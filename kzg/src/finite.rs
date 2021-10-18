use super::Fr;

pub fn add_fr(first: Fr, second: Fr) -> Fr {
    let mut sum = Fr::default();
    unsafe {
        blst::blst_fr_add(&mut sum.0, &first.0, &second.0);
    }
    sum
}
