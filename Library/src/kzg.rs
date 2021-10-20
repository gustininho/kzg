#![allow(non_camel_case_types)]
use ark_poly_commit::*;
use ark_poly_commit::kzg10::{KZG10, UniversalParams, Powers, VerifierKey, Randomness, Proof, Commitment};

use super::{Fp, Fr as BlstFr, P1, P2};
use super::utils::{
    blst_fr_into_pc_fr, blst_p1_into_pc_g1projective, blst_poly_into_pc_poly, pc_fr_into_blst_fr,
    pc_g1projective_into_blst_p1, polydata,
};
use ark_bls12_381::{g1, Bls12_381, Fr};
use ark_ec::{models::short_weierstrass_jacobian::GroupAffine, PairingEngine, AffineCurve, ProjectiveCurve};
use ark_ec::msm::{VariableBaseMSM};
use ark_poly::univariate::DensePolynomial as DensePoly;
use ark_std::{marker::PhantomData, ops::Div, vec, test_rng, Zero,};
use ark_ff::{Field, PrimeField};
use rand::rngs::StdRng;
use std::ops::{Neg,};

use super::fft::SCALE2_ROOT_OF_UNITY;

type UniPoly_381 = DensePoly<<Bls12_381 as PairingEngine>::Fr>;
type KZG_Bls12_381 = KZG10<Bls12_381, UniPoly_381>;

pub fn trim(
    pp: &UniversalParams<Bls12_381>,
    mut supported_degree: usize,
) -> Result<(Powers<Bls12_381>, VerifierKey<Bls12_381>), Error> {
    if supported_degree == 1 {
        supported_degree += 1;
    }
    let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
    let powers_of_gamma_g = (0..=supported_degree)
        .map(|i| pp.powers_of_gamma_g[&i])
        .collect();

    let powers = Powers {
        powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
        powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
    };
    let vk = VerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[&0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    Ok((powers, vk))
}

pub struct KZG<E: PairingEngine, P: UVPolynomial<E::Fr>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E, P> KZG<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{

    pub fn open<'a>(
        powers: &Powers<Bls12_381>,
        p: &DensePoly<Fr>,
        point: Fr,
        rand: &Randomness<Fr, DensePoly<Fr>>,
    ) -> Result<Proof<Bls12_381>, Error> {
        Self::check_degree_is_too_large(p.degree(), powers.size())?;

        let (witness_poly, hiding_witness_poly) = KZG_Bls12_381::compute_witness_polynomial(p, point, rand)?;

        let proof = Self::open_with_witness_polynomial(
            powers,
            point,
            rand,
            &witness_poly,
            hiding_witness_poly.as_ref(),
        );

        proof
    }

    pub(crate) fn open_with_witness_polynomial<'a>(
        powers: &Powers<Bls12_381>,
        point: Fr,
        randomness: &Randomness<Fr, DensePoly<Fr>>,
        witness_polynomial: &DensePoly<Fr>,
        hiding_witness_polynomial: Option<&DensePoly<Fr>>,
    ) -> Result<Proof<Bls12_381>, Error> {
        Self::check_degree_is_too_large(witness_polynomial.degree(), powers.size())?;
        let (num_leading_zeros, witness_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(witness_polynomial);

        let mut w = VariableBaseMSM::multi_scalar_mul(
            &powers.powers_of_g[num_leading_zeros..],
            &witness_coeffs,
        );

        let random_v = if let Some(hiding_witness_polynomial) = hiding_witness_polynomial {
            let blinding_p = &randomness.blinding_polynomial;
            let blinding_evaluation = blinding_p.evaluate(&point);

            let random_witness_coeffs = convert_to_bigints(&hiding_witness_polynomial.coeffs());
            let witness_comm_time =
            w += &VariableBaseMSM::multi_scalar_mul(
                &powers.powers_of_gamma_g,
                &random_witness_coeffs,
            );
            Some(blinding_evaluation)
        } else {
            None
        };

        Ok(Proof {
            w: w.into_affine(),
            random_v,
        })
    }

    pub(crate) fn check_degree_is_too_large(degree: usize, num_powers: usize) -> Result<(), Error> {
        let num_coefficients = degree + 1;
        if num_coefficients > num_powers {
            Err(Error::TooManyCoefficients {
                num_coefficients,
                num_powers,
            })
        } else {
            Ok(())
        }
    }
}
fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField, P: UVPolynomial<F>>(
    p: &P,
) -> (usize, Vec<F::BigInt>) {
    let mut num_leading_zeros = 0;
    while num_leading_zeros < p.coeffs().len() && p.coeffs()[num_leading_zeros].is_zero() {
        num_leading_zeros += 1;
    }
    let coeffs = convert_to_bigints(&p.coeffs()[num_leading_zeros..]);
    (num_leading_zeros, coeffs)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let coeffs = p.iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    coeffs
}

pub(crate) const fr_one: BlstFr = BlstFr {
    l: [
        0xc999e990f3f29c6d,
        0x2b6cedcb87925c23,
        0x5d314967254398f,
        0x748d9d99f59ff11,
    ],
};

pub(crate) const fr_zero: BlstFr = BlstFr {
    l: [0x0, 0x0, 0x0, 0x0],
};

pub(crate) const g1_identity: P1 = P1 {
    x: Fp {
        l: [0, 0, 0, 0, 0, 0],
    },
    y: Fp {
        l: [0, 0, 0, 0, 0, 0],
    },
    z: Fp {
        l: [0, 0, 0, 0, 0, 0],
    },
};

pub struct FFTSettings {
    pub max_width: u64,
    pub root_of_unity: BlstFr,
    pub expanded_roots_of_unity: Vec<BlstFr>,
    pub reverse_roots_of_unity: Vec<BlstFr>,
}

impl Default for FFTSettings {
    fn default() -> FFTSettings {
        FFTSettings {
            max_width: 0,
            root_of_unity: BlstFr { l: [0, 0, 0, 0] },
            expanded_roots_of_unity: Vec::new(),
            reverse_roots_of_unity: Vec::new(),
        }
    }
}

impl FFTSettings {
    pub fn from_scale(max_scale: usize) -> Result<FFTSettings, String> {
        if max_scale >= SCALE2_ROOT_OF_UNITY.len() {
            return Err(String::from("Scale is expected to be within root of unity matrix row size"));
        }
        let max_width: u64 = 1 << max_scale;

        Ok(FFTSettings {
            max_width: max_width,
            ..Default::default()
        })
    }
}

pub(crate) struct KZGSettings {
    pub(crate) fs: FFTSettings,
    pub(crate) secret_g1: Vec<P1>,
    pub(crate) secret_g2: Vec<P2>,
    pub(crate) length: u64,
    pub(crate) params: UniversalParams<Bls12_381>,
    pub(crate) rand: StdRng,
    pub(crate) rand2: Randomness<Fr, UniPoly_381>,
}

impl Default for KZGSettings {
    fn default() -> KZGSettings {
        KZGSettings {
            fs: FFTSettings::default(),
            secret_g1: Vec::new(),
            secret_g2: Vec::new(),
            length: 0,
            params: KZG_Bls12_381::setup(1, false, &mut test_rng()).unwrap(),
            rand: test_rng(),
            rand2: Randomness::empty(),
        }
    }
}

pub(crate) fn generate_trusted_setup(len: usize) -> Result<(Vec<P1>, Vec<P2>), Error> {
    Ok((Vec::new(), Vec::new()))
}

pub(crate) fn new_fft_settings(max_scale: u64) -> FFTSettings {
    FFTSettings::default()
}

pub(crate) fn new_kzg_settings(
    secret_g1: &Vec<P1>,
    secret_g2: &Vec<P2>,
    length: u64,
    fs: FFTSettings,
) -> KZGSettings {
    KZGSettings {
        length: length,
        params: KZG_Bls12_381::setup(length as usize, true, &mut test_rng()).unwrap(),
        ..Default::default()
    }
}

pub(crate) fn fr_from_uint64(num: u64) -> BlstFr {
    let fr = Fr::from(num);
    pc_fr_into_blst_fr(fr)
}

pub(crate) fn new_poly(len: usize) -> polydata {
    polydata {
        coeffs: vec![BlstFr { l: [0, 0, 0, 0] }; len],
    }
}

#[derive(Debug, PartialEq)]
pub(crate) struct TooLongPoly;

pub(crate) fn commit_to_poly(p: &polydata, ks: &mut KZGSettings) -> Result<P1, TooLongPoly> {
    if p.coeffs.len() > ks.length as usize {
        Err(TooLongPoly)
    } else if blst_poly_into_pc_poly(&p).unwrap().is_zero() {
        Ok(g1_identity)
    } else {
        let (powers, _) = trim(&ks.params, &ks.params.max_degree() - 1).unwrap();
        let (com, rand) = KZG_Bls12_381::commit(
            &powers,
            &blst_poly_into_pc_poly(&p).unwrap(),
            None,
            Some(&mut ks.rand),
        )
        .unwrap();
        ks.rand2 = rand;
        Ok(pc_g1projective_into_blst_p1(com.0.into_projective()).unwrap())
    }
}

pub(crate) fn compute_proof_single(p: &polydata, x: &BlstFr, ks: &KZGSettings) -> Proof<Bls12_381> {
    let (powers, _) = trim(&ks.params, &ks.params.max_degree() - 1).unwrap();
    
    let proof = KZG::<Bls12_381, UniPoly_381>::open(
        &powers,
        &blst_poly_into_pc_poly(p).unwrap(),
        blst_fr_into_pc_fr(x),
        &ks.rand2,
    )
    .unwrap();
    proof
}

pub(crate) fn eval_poly(p: &polydata, x: &BlstFr) -> BlstFr {
    let poly = blst_poly_into_pc_poly(p).unwrap();
    pc_fr_into_blst_fr(poly.evaluate(&blst_fr_into_pc_fr(x)))
}

pub(crate) fn check_proof_single(
    com: &P1,
    proof: &Proof<Bls12_381>,
    x: &BlstFr,
    value: &BlstFr,
    ks: &KZGSettings,
) -> bool {
    let (powers, vk) = trim(&ks.params, &ks.params.max_degree() - 1).unwrap();
    let projective = blst_p1_into_pc_g1projective(com).unwrap();
    let affine = GroupAffine::<g1::Parameters>::from(projective.clone());
    let mut com = Commitment::empty();
    com.0 = affine;
    KZG_Bls12_381::check(
        &vk,
        &com,
        blst_fr_into_pc_fr(x),
        blst_fr_into_pc_fr(value),
        &proof,
    )
    .unwrap()
}

pub(crate) fn fr_add(x: &BlstFr, y: &BlstFr) -> BlstFr {
    let pcx = blst_fr_into_pc_fr(x);
    let pcy = blst_fr_into_pc_fr(y);
    let sum = pcx + pcy;
    pc_fr_into_blst_fr(sum)
}

pub(crate) fn fr_mul(x: &BlstFr, roots: &BlstFr) -> BlstFr {
    let PCx = blst_fr_into_pc_fr(x);
    let PCroots = blst_fr_into_pc_fr(roots);
    let mul = PCx * PCroots;
    pc_fr_into_blst_fr(mul)
}

pub(crate) fn compute_proof_multi(
    p: &polydata,
    x: &BlstFr,
    n: u64,
    ks: &mut KZGSettings,
) -> Proof<Bls12_381> {
    let rng = &mut test_rng();
    let mut PCdivisor = DensePoly::rand(n as usize, rng);

    let PCx = blst_fr_into_pc_fr(x);

    let x_pow_n = PCx.pow(&[n]);
    let fr_neg = x_pow_n.neg();
    let temp = std::mem::replace(&mut PCdivisor.coeffs[0], blst_fr_into_pc_fr(&fr_zero));
    for x in 1..n - 1 {
        let temp = std::mem::replace(&mut PCdivisor.coeffs[x as usize], fr_neg);
    }
    let temp = std::mem::replace(
        &mut PCdivisor.coeffs[n as usize],
        blst_fr_into_pc_fr(&fr_one),
    );

    let q = &blst_poly_into_pc_poly(p).unwrap() / &PCdivisor;

    let (powers, _) = trim(&ks.params, &ks.params.max_degree() - 1).unwrap();
    let (com, rand) = KZG_Bls12_381::commit(&powers, &q, None, Some(&mut ks.rand)).unwrap();
    ks.rand2 = rand;
    Proof {
        w: com.0,
        random_v: None,
    }
}

//TODO:
pub(crate) fn check_proof_multi(
    com: &P1,
    proof: &Proof<Bls12_381>,
    x: &BlstFr,
    ys: &Vec<BlstFr>,
    n: u64,
    ks: &KZGSettings,
) -> bool {
    // let rng = &mut test_rng();
    // let mut interp: DensePoly<Fr> = DensePoly::rand(n as usize, rng);
    true
}
