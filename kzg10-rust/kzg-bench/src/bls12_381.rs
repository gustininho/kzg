use mcl_rust::old::*;
use mcl_rust::CurveType;
use mcl_rust::data_types::fr::Fr;
use mcl_rust::mlc_methods::*;
use mcl_rust::utilities::*;

#[test]
fn log_2_byte_works() {
    assert!(0 == log_2(0x01));
    assert!(7 == log_2(0x80));
    assert!(7 == log_2(0xff));
    assert!(4 == log_2(0x10));
}

