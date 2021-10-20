
use crate::kzg_types::{FFTSettings, Poly, create_fr_zero, negate_fr, fr_is_zero, create_fr_u64};
use blst::{
    blst_fr_eucl_inverse,
    blst_fr_mul,
    blst_fr_sub,
    blst_fr_add,
};
use kzg::{Fr};

use crate::utils::{next_power_of_two, log2_pow2, log2_u64, min_u64};
use crate::fft_fr::fft_fr;

pub fn fr_div(a: &Fr, b: &Fr) -> Result<Fr, String> {
    let mut tmp = Fr::default();
    let mut out = Fr::default();

    unsafe {
        blst_fr_eucl_inverse(&mut tmp, b);
        blst_fr_mul(&mut out, a, &tmp);
    }
    Ok(out)
}

pub fn poly_norm(p: &Poly) -> Result<Poly, String> {
    let mut ret = p.clone();
    let mut temp_len: usize = ret.coeffs.len();
    while temp_len > 0 && fr_is_zero(&ret.coeffs[temp_len - 1]) {
        temp_len -= 1;
    }
    if temp_len == 0 {
        ret.coeffs = Vec::default();
    }

    Ok(ret)
}

pub fn poly_quotient_length(dividend: &Poly, divisor: &Poly) -> Result<usize, String>{
    if dividend.coeffs.len() >= divisor.coeffs.len() {
        return Ok(dividend.coeffs.len() - divisor.coeffs.len() + 1);
    }

    Ok(0)
}

pub fn poly_long_div(dividend: &Poly, divisor: &Poly) -> Result<Poly, String> {
    // Dividing by zero is undefined
    assert!(divisor.coeffs.len() > 0);
    // The divisor's highest coefficient must be non-zero
    //CHECK(!fr_is_zero(&divisor->coeffs[divisor->length - 1]));
    assert!(!fr_is_zero(&divisor.coeffs[divisor.coeffs.len() - 1]));

    let mut out: Poly = Poly { coeffs: Vec::default()};
    //uint64_t a_pos = dividend->length - 1;
    //uint64_t b_pos = divisor->length - 1;
    //uint64_t diff = a_pos - b_pos;
    let mut a_pos = dividend.coeffs.len();
    let b_pos = divisor.coeffs.len();
    let mut diff = a_pos - b_pos;
    
    // Deal with the size of the output polynomial
    // uint64_t out_length = poly_quotient_length(dividend, divisor);
    let result = poly_quotient_length(&dividend, &divisor);
    assert!(result.is_ok());
    let out_length = result.unwrap();
    
    // CHECK(out->length >= out_length);
    assert!(out.coeffs.len() >= out_length);
    
    // If the divisor is larger than the dividend, the result is zero-length
    if out_length == 0 {
        return Ok(out);
    }
    
    //fr_t *a;
    // TRY(new_fr_array(&a, dividend->length));
    let mut a = vec![Fr::default(); dividend.coeffs.len()];
    for i in 0..dividend.coeffs.len() {
        //a[i] = dividend->coeffs[i];
        a.push(dividend.coeffs[i]);
    }

    while diff > 0 {
        // fr_div(&out->coeffs[diff], &a[a_pos], &divisor->coeffs[b_pos]);
        let result = fr_div(&a[a_pos], &divisor.coeffs[b_pos]);
        assert!(result.is_ok());
        out.coeffs[diff] = result.unwrap();
    
        unsafe {
            for i in 0..(b_pos+1) {
                // fr_t tmp;
                let mut tmp = Fr::default();
                // a[diff + i] -= b[i] * quot
                blst_fr_mul(&mut tmp, &out.coeffs[diff], &divisor.coeffs[i]);
                blst_fr_sub(&mut a[diff + i], &a[diff + i], &tmp);
            }
        }
        diff -= 1;
        a_pos -= 1;
    }
    // fr_div(&out->coeffs[0], &a[a_pos], &divisor->coeffs[b_pos]);
    let result = fr_div(&a[a_pos], &divisor.coeffs[b_pos]);
    assert!(result.is_ok());
    out.coeffs[0] = result.unwrap();

    Ok(out)
}

pub fn poly_flip(input: &Poly) -> Result<Poly, String> {
    let mut output = Poly { coeffs: Vec::default() };
    for i in 0..input.coeffs.len() {
        output.coeffs.push(input.coeffs[input.coeffs.len() - i - 1]);
    }
    Ok(output)
}

pub fn poly_inverse(b: &Poly, output_len: usize) -> Result<Poly, String> {
    assert!(b.coeffs.len() > 0);
    assert!(!fr_is_zero(&b.coeffs[0]));

    let mut output = Poly { coeffs: Vec::default() };
    // If the input polynomial is constant, the remainder of the series is zero
    if b.coeffs.len() == 1 {
        unsafe {
            blst_fr_eucl_inverse(&mut output.coeffs[0], &b.coeffs[0]);
        }
        for i in 1..output_len {
            output.coeffs[i] = create_fr_zero();
        }
        return Ok(output);
    }
    
    //poly tmp0, tmp1;
    //fr_t fr_two;

    // uint64_t length = out->length;
    // uint64_t maxd = length - 1;
    let maxd = output_len - 1;

    //uint64_t d = 0;

    // Max space for multiplications is (2 * length - 1)
    // int scale = log2_pow2(next_power_of_two(2 * length - 1));
    let scale: usize = log2_pow2(next_power_of_two(2 * output_len - 1));
    // FFTSettings fs;
    // TRY(new_fft_settings(&fs, scale));

    let fs: FFTSettings = FFTSettings::from_scale(scale).unwrap();

    // To store intermediate results
    // TRY(new_poly(&tmp0, length));
    // TRY(new_poly(&tmp1, length));
    let mut tmp0 = Poly { coeffs: Vec::default() };
    let mut tmp1 = Poly { coeffs: Vec::default() };

    // Base case for d == 0
    unsafe {
        blst_fr_eucl_inverse(&mut output.coeffs[0], &b.coeffs[0]);
    }
    //output->length = 1;

    let mut d: usize = 0;

    // uint64_t mask = (uint64_t)1 << log2_u64(maxd);
    let mut mask: usize = 1 << log2_u64(maxd);
    while mask != 0 {

        d = 2 * d + (if (maxd & mask) != 0 {1} else {0});
        mask = mask >> 1;

        // b.c -> tmp0 (we're using out for c)
        // tmp0.length = min_u64(d + 1, b->length + output->length - 1);
        let len_temp;
        if d + 1 < b.coeffs.len() + output_len - 1 {
            len_temp = d + 1;
        } else {
            len_temp = b.coeffs.len() + output_len - 1
        }

        // TRY(poly_mul_(&tmp0, b, output, &fs));
        tmp0 = poly_mul_(&b, &output, &fs, len_temp).unwrap();

        // 2 - b.c -> tmp0
        for i in 0..tmp0.coeffs.len() {
            // fr_negate(&tmp0.coeffs[i], &tmp0.coeffs[i]);
            let cloned_fr = tmp0.coeffs[i].clone();
            negate_fr(&mut tmp0.coeffs[i], &cloned_fr);
        }
        // fr_from_uint64(&fr_two, 2);
        let fr_two: Fr = create_fr_u64(2);
        unsafe {
            blst_fr_add(&mut tmp0.coeffs[0], &tmp0.coeffs[0], &fr_two);
        }

        // c.(2 - b.c) -> tmp1;
        // tmp1.length = d + 1;
        // TRY(poly_mul_(&tmp1, output, &tmp0, &fs));
        tmp1 = poly_mul_(&output, &tmp0, &fs, d + 1).unwrap();

        // tmp1 -> c
        // output->length = tmp1.length;
        for i in 0..tmp1.coeffs.len() {
            output.coeffs.push(tmp1.coeffs[i]);
        }
    }
    assert!(d + 1 == output_len);
    Ok(output)
}

pub fn poly_fast_div(dividend: &Poly, divisor: &Poly) -> Result<Poly, String> {

    // Dividing by zero is undefined
    assert!(divisor.coeffs.len() > 0);

    // The divisor's highest coefficient must be non-zero
    assert!(!fr_is_zero(&divisor.coeffs[divisor.coeffs.len() - 1]));

    let m: usize = dividend.coeffs.len() - 1;
    let n: usize = divisor.coeffs.len() - 1;

    // If the divisor is larger than the dividend, the result is zero-length
    if n > m {
        return Ok(Poly { coeffs: Vec::default() });
    }

    // Ensure the output poly has enough space allocated
    //CHECK(out->length >= m - n + 1);

    // Ensure that the divisor is well-formed for the inverse operation
    assert!(!fr_is_zero(&divisor.coeffs[divisor.coeffs.len() - 1]));

    let mut out = Poly { coeffs: Vec::default() };
    // Special case for divisor.length == 1 (it's a constant)
    if divisor.coeffs.len() == 1 {
        //out->length = dividend->length;
        for i in 0..dividend.coeffs.len() {
            out.coeffs.push(fr_div(&dividend.coeffs[i], &divisor.coeffs[0]).unwrap());
        }
        return Ok(out);
    }

    // poly a_flip, b_flip;
    let mut a_flip = Poly { coeffs: Vec::default() };
    let mut b_flip = Poly { coeffs: Vec::default() };

    // TRY(new_poly(&a_flip, dividend->length));
    // TRY(new_poly(&b_flip, divisor->length));
    // TRY(poly_flip(&a_flip, dividend));
    // TRY(poly_flip(&b_flip, divisor));
    a_flip = poly_flip(&dividend).unwrap();
    b_flip = poly_flip(&divisor).unwrap();

    // poly inv_b_flip;
    let mut inv_b_flip = Poly { coeffs: Vec::default() };
    // TRY(new_poly(&inv_b_flip, m - n + 1));
    // TRY(poly_inverse(&inv_b_flip, &b_flip));
    inv_b_flip = poly_inverse(&b_flip, m - n + 1).unwrap();

    // poly q_flip;
    let mut q_flip = Poly { coeffs: Vec::default() };
    // We need only m - n + 1 coefficients of q_flip
    // TRY(new_poly(&q_flip, m - n + 1));
    // TRY(poly_mul(&q_flip, &a_flip, &inv_b_flip));

    q_flip = poly_mul(&a_flip, &inv_b_flip, m - n + 1).unwrap();

    // out->length = m - n + 1;
    // TRY(poly_flip(out, &q_flip));
    out = poly_flip(&q_flip).unwrap();

    Ok(out)
}

pub fn poly_mul_direct(a: &Poly, b: &Poly, output_len: usize) -> Result<Poly, String>{

    //uint64_t a_degree = a->length - 1;
    //uint64_t b_degree = b->length - 1;
    let a_degree: usize = a.coeffs.len() - 1;
    let b_degree: usize = b.coeffs.len() - 1;
    let mut output = Poly { coeffs: Vec::default() };

    for _ in 0..output_len {
        output.coeffs.push(create_fr_zero());
    }

    // Truncate the output to the length of the output polynomial
    for i in 0..(a_degree + 1) {
        let mut j: usize = 0;
        unsafe {
            while j <= b_degree && i + j < output.coeffs.len() {
                let mut tmp: Fr = Fr::default();
                blst_fr_mul(&mut tmp, &a.coeffs[i], &b.coeffs[j]);
                blst_fr_add(&mut output.coeffs[i + j], &output.coeffs[i + j], &tmp);
                
                j += 1;
            }
        }
    }

    Ok(output)
}

pub fn pad(input: &Vec<Fr>, n_in: usize, n_out: usize) -> Result<Vec<Fr>, String> {
    // uint64_t num = min_u64(n_in, n_out);
    let num: usize = min_u64(n_in, n_out).unwrap();
    let mut output: Vec<Fr> = Vec::default();
    for i in 0..num {
        output[i] = input[i].clone();
    }
    for i in num..n_out {
        output[i] = create_fr_zero();
    }
    Ok(output)
}

pub fn poly_mul_fft(a: &Poly, b: &Poly, fs_: Option<&FFTSettings>, output_len: usize) -> Result<Poly, String>{

    // Truncate a and b so as not to do excess work for the number of coefficients required.
    // uint64_t a_len = min_u64(a->length, out->length);
    // uint64_t b_len = min_u64(b->length, out->length);
    // uint64_t length = next_power_of_two(a_len + b_len - 1);

    let a_len = min_u64(a.coeffs.len(), output_len).unwrap();
    let b_len = min_u64(a.coeffs.len(), output_len).unwrap();
    let length =  next_power_of_two(a_len + b_len - 1);

    // If the FFT settings are NULL then make a local set, otherwise use the ones passed in.
    // FFTSettings fs, *fs_p;

    let mut fs: FFTSettings = FFTSettings::from_scale(0).unwrap();
    match fs_ {
        Some(x) => fs = x.clone(),
        None    => {
            let scale: usize = log2_pow2(length);
            fs = FFTSettings::from_scale(scale).unwrap();
        }
    }

    // if (fs_ != NULL) {
    //     fs_p = fs_;
    // } else {
    //     fs_p = &fs;
    //     int scale = log2_pow2(length); // TODO only good up to length < 32 bits
    //     TRY(new_fft_settings(fs_p, scale));
    // }

    // CHECK(length <= fs_p->max_width);
    assert!(length <= fs.max_width);

    // fr_t *a_pad, *b_pad, *a_fft, *b_fft;
    let mut a_pad: Vec<Fr> = Vec::default();
    let mut b_pad: Vec<Fr> = Vec::default();
    let mut a_fft: Vec<Fr> = Vec::default();
    let mut b_fft: Vec<Fr> = Vec::default();

    // TRY(new_fr_array(&a_pad, length));
    // TRY(new_fr_array(&b_pad, length));
    a_pad = pad(&a.coeffs, a_len, length).unwrap();
    b_pad = pad(&b.coeffs, b_len, length).unwrap();

    // TRY(new_fr_array(&a_fft, length));
    // TRY(new_fr_array(&b_fft, length));
    // TRY(fft_fr(a_fft, a_pad, false, length, fs_p));
    // TRY(fft_fr(b_fft, b_pad, false, length, fs_p));
    a_fft = fft_fr(&a_pad, false, &fs).unwrap();
    b_fft = fft_fr(&b_pad, false, &fs).unwrap();

    // fr_t *ab_fft = a_pad; // reuse the a_pad array
    // fr_t *ab = b_pad;     // reuse the b_pad array
    let mut ab_fft: Vec<Fr> = a_pad.clone();
    // let mut ab: &Vec<Fr> = &b_pad;
    
    unsafe {
        for i in 0..length {
            blst_fr_mul(&mut ab_fft[i], &a_fft[i], &b_fft[i]);
        }
    }
    // TRY(fft_fr(ab, ab_fft, true, length, fs_p));
    let ab = &fft_fr(&ab_fft, true, &fs).unwrap();

    // Copy result to output
    // uint64_t data_len = min_u64(out->length, length);
    let mut output = Poly { coeffs: Vec::default() };
    let data_len = min_u64(output_len, length).unwrap();
    for i in 0..data_len {
        output.coeffs.push(ab[i]);
    }
    for _ in data_len..output_len {
        output.coeffs.push(create_fr_zero());
    }

    return Ok(output);
}

pub fn poly_mul_(a: &Poly, b: &Poly, fs: &FFTSettings, output_len: usize) -> Result<Poly, String>{
    if a.coeffs.len() < 64 || b.coeffs.len() < 64 || output_len < 128 { // Tunable parameter
        return poly_mul_direct(&a, &b, output_len);
    } else {
        return poly_mul_fft(a, b, Some(&fs), output_len);
    }
}

pub fn poly_mul(a: &Poly, b: &Poly, output_len: usize) -> Result<Poly, String> {
    let fft_settings = FFTSettings::from_scale(0).unwrap();
    return poly_mul_(&a, &b, &fft_settings, output_len);
}

pub fn new_poly_div(dividend_: &Poly, divisor_: &Poly) -> Result<Poly, String>{

    //poly dividend = poly_norm(dividend_);
    let result = poly_norm(dividend_);
    assert!(result.is_ok());
    let dividend: Poly = result.unwrap();
    
    //poly divisor = poly_norm(divisor_);
    let result = poly_norm(divisor_);
    assert!(result.is_ok());
    let divisor: Poly = result.unwrap();

    //TRY(new_poly(out, poly_quotient_length(&dividend, &divisor)));
    // let newLength = poly_quotient_length(&dividend, &divisor).unwrap();

    if divisor.coeffs.len() >= dividend.coeffs.len() || divisor.coeffs.len() < 128 { // Tunable paramter
        let result = poly_long_div(&dividend, &divisor);
        assert!(result.is_ok());
        result
    } else {
        let result = poly_fast_div(&dividend, &divisor);
        assert!(result.is_ok());
        result
    }
}
