use std::{mem, vec};
use mcl_rust::old::*;
use mcl_rust::utilities::*;
use mcl_rust::CurveType;
use mcl_rust::data_types::fr::Fr;
use mcl_rust::mlc_methods::init;

#[test]
fn alog_2_byte_works() {
    assert!(0 == log_2(0x01));
    assert!(7 == log_2(0x80));
    assert!(7 == log_2(0xff));
    assert!(4 == log_2(0x10));
}

#[test]
fn fr_is_zero_works() {
    assert!(init(CurveType::BLS12_381));
    let x = Fr::from_int(0);
    assert!(x.is_zero());
}

#[test]
fn fr_is_one_works() {
    assert!(init(CurveType::BLS12_381));
    let x = Fr::from_int(1);
    assert!(x.is_one());
}