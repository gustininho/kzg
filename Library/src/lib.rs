pub type Pairing = blst::Pairing;
pub type Fp = blst::blst_fp;
pub type Fp12 = blst::blst_fp12;
pub type Fp6 = blst::blst_fp6;
pub type Fr = blst::blst_fr;
pub type P1 = blst::blst_p1;
pub type P1Affine = blst::blst_p1_affine;
pub type P2 = blst::blst_p2;
pub type P2Affine = blst::blst_p2_affine;
pub type Scalar = blst::blst_scalar;
pub type Uniq = blst::blst_uniq;

//pub mod finite;
mod utils;
mod tests;
mod fft;
mod main;
mod kzg;

pub(crate) trait Eq:  {
    fn equals(&self, other: &P1) -> bool;
}

impl Eq for P1 {
    fn equals(&self, other: &P1) -> bool {
        if self.x.l.eq(&other.x.l) && self.y.l.eq(&other.x.l) && self.z.l.eq(&other.x.l) {
            true
        } else{
            false
        }
    }
}