use ark_bls12_381::{g1, g2};
use ark_bls12_381::Fq;
use ark_bls12_381::Fr;
use ark_ec::{
     models::short_weierstrass_jacobian::GroupAffine, models::short_weierstrass_jacobian::GroupProjective
};
use ark_ff::{biginteger::BigInteger384, biginteger::BigInteger256, Fp2, Fp384};
use ark_poly::univariate::DensePolynomial as DensePoly;
use ark_poly::UVPolynomial;
use super::{Fp, Fr as BlstFr, P1, P2, P1Affine};
use crate::*;

#[derive(Debug, PartialEq)]
pub struct Error;

pub(crate) struct polydata {
    pub(crate) coeffs: Vec<BlstFr>,
}

pub(crate) fn pc_poly_into_blst_poly(poly: DensePoly<Fr>) -> Result<polydata, Error> {
        let mut bls_pol = polydata { coeffs: Vec::new() };

        for x in poly.coeffs {
            bls_pol.coeffs.push(BlstFr{l:x.0.0});
        }

        Ok(bls_pol)
}

pub(crate) fn blst_poly_into_pc_poly(pd: &polydata) -> Result<DensePoly<Fr>, Error> {
        let mut poly = Vec::new();
        let x = pd.coeffs.clone();
        for x in x {
            poly.push(Fr::from(BigInteger256::new(x.l)))
        }

        let p = DensePoly::from_coefficients_slice(&poly);
        Ok(p)
}

pub(crate) fn pc_fr_into_blst_fr(fr: Fr) -> BlstFr {
        BlstFr { l: fr.0 .0 }
}

pub(crate) fn pc_fq_into_blst_fp(fq: Fq) -> Fp {
        Fp { l: fq.0.0}
}

pub(crate) fn blst_fr_into_pc_fr(fr: &BlstFr) -> Fr {
        let big_int = BigInteger256::new(fr.l);
        Fr::new(big_int)
}

pub(crate) fn blst_fp_into_pc_fq(fp: &Fp) -> Fq {
        let big_int = BigInteger384::new(fp.l);
        Fq::new(big_int)
}

pub(crate) fn pc_g1projective_into_blst_p1(
        gp: GroupProjective<g1::Parameters>
    ) -> Result<P1, Error> {
        let p1 = P1{
            x: pc_fq_into_blst_fp(gp.x),
            y: pc_fq_into_blst_fp(gp.y),
            z: pc_fq_into_blst_fp(gp.z),
        };
        Ok(p1)
}

pub(crate) fn blst_p1_into_pc_g1projective(
        p1: &P1
    ) -> Result<GroupProjective<g1::Parameters>, Error> {
        let pc_projective = GroupProjective::new(
            blst_fp_into_pc_fq(&p1.x),
            blst_fp_into_pc_fq(&p1.y),
            blst_fp_into_pc_fq(&p1.z),
        );
        Ok(pc_projective)
}

pub(crate) fn blst_p2_into_pc_g2projective(
    p2: &P2
) -> Result<GroupProjective<g2::Parameters>, Error> {
    let pc_projective = GroupProjective::new(
        Fp2::new(Fp384::new(BigInteger384::new(p2.x.fp[0].l)), Fp384::new(BigInteger384::new(p2.x.fp[1].l))),
        Fp2::new(Fp384::new(BigInteger384::new(p2.y.fp[0].l)), Fp384::new(BigInteger384::new(p2.y.fp[1].l))),
        Fp2::new(Fp384::new(BigInteger384::new(p2.z.fp[0].l)), Fp384::new(BigInteger384::new(p2.z.fp[1].l))),
    );
    Ok(pc_projective)
}

pub(crate) fn pc_affine_into_blst_affine(
    affine: GroupAffine<g1::Parameters>,
) -> Result<P1Affine, Error> {
    let bl_aff = P1Affine {
        x: Fp { l: affine.x.0 .0 },
        y: Fp { l: affine.y.0 .0 },
    };

    Ok(bl_aff)
}

pub(crate) fn blst_affine_into_pc_affine(
    affine: &P1Affine
) -> Result<GroupAffine<g1::Parameters>, Error> {
    let pc_affine = GroupAffine::new(
        blst_fp_into_pc_fq(&affine.x),
        blst_fp_into_pc_fq(&affine.y),
        false,
    );
    Ok(pc_affine)
}