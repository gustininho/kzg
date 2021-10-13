use ark_bls12_381::{g1, g2};
use ark_bls12_381::Fq;
use ark_bls12_381::Fr;
use ark_ff::{biginteger::BigInteger384, biginteger::BigInteger256, One, PrimeField, UniformRand, Zero, Fp2, Fp384};
use ark_poly::univariate::DensePolynomial as DensePoly;
use ark_poly::UVPolynomial;
use num_bigint::BigUint;
use super::*;

pub struct polydata {
    coeffs: Vec<u64>,
}

pub fn blst_poly_into_pc_poly(pd: polydata) -> Result<DensePoly<Fr>, Error> {
    let mut poly = Vec::new();
    for x in pd.coeffs {
        poly.push(Fr::from(x))
    }

    let p = DensePoly::from_coefficients_slice(&poly);
    Ok(p)
}

pub fn pc_poly_into_blst_poly(poly: DensePoly<Fr>) -> Result<polydata, Error> {
    let mut bls_pol = polydata { coeffs: Vec::new() };

    for x in poly.coeffs {
        let bigU: BigUint = x.into();
        bls_pol.coeffs.push(bigU.to_u64_digits()[0])
    }

    Ok(bls_pol)
}

pub fn pc_affine_into_blst_affine(
    affine: GroupAffine<g1::Parameters>,
) -> Result<super::P1Affine, Error> {
    let bl_aff = super::P1Affine {
        x: blst_fp { l: affine.x.0 .0 },
        y: blst_fp { l: affine.y.0 .0 },
    };

    Ok(bl_aff)
}

pub fn blst_fp_into_pc_fq(fp: super::Fp) -> Fq {
    let big_int = BigInteger384::new(fp.l);
    Fq::new(big_int)
}

pub fn blst_fr_into_pc_fr(fr: super::Fr) -> Fr {
    let big_int = BigInteger256::new(fr.l);
    Fr::new(big_int)
}

pub fn pc_fr_into_blst_fr(fr: Fr) -> super::Fr {
    super::Fr { l: fr.0 .0 }
}

pub fn blst_affine_into_pc_affine(
    affine: super::P1Affine
) -> Result<GroupAffine<g1::Parameters>, Error> {
    let pc_affine = GroupAffine::new(
        blst_fp_into_pc_fq(affine.x),
        blst_fp_into_pc_fq(affine.y),
        false,
    );
    Ok(pc_affine)
}

pub fn blst_p1_into_pc_g1projective(
    p1: super::P1
) -> Result<GroupProjective<g1::Parameters>, Error> {
    let pc_projective = GroupProjective::new(
        blst_fp_into_pc_fq(p1.x),
        blst_fp_into_pc_fq(p1.y),
        blst_fp_into_pc_fq(p1.z),
    );
    Ok(pc_projective)
}

pub fn blst_p2_into_pc_g2projective(
    p2: super::P2
) -> Result<GroupProjective<g2::Parameters>, Error> {
    let pc_projective = GroupProjective::new(
        Fp2::new(Fp384::new(BigInteger384::new(p2.x.fp[0].l)), Fp384::new(BigInteger384::new(p2.x.fp[1].l))),
        Fp2::new(Fp384::new(BigInteger384::new(p2.y.fp[0].l)), Fp384::new(BigInteger384::new(p2.y.fp[1].l))),
        Fp2::new(Fp384::new(BigInteger384::new(p2.z.fp[0].l)), Fp384::new(BigInteger384::new(p2.z.fp[1].l))),
    );
    Ok(pc_projective)
}