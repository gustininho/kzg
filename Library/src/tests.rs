use crate::fft::fr_are_equal;
use crate::fft::fft_fr;
use blst::blst_fr_from_uint64;
use crate::kzg::{generate_trusted_setup, new_kzg_settings, KZGSettings, FFTSettings, new_fft_settings, fr_from_uint64,
    eval_poly, fr_add, fr_one, commit_to_poly, check_proof_single, new_poly, compute_proof_single, compute_proof_multi,
    fr_mul, check_proof_multi, g1_identity, TooLongPoly};
use crate::{Fr, P1, Eq};

#[test]
pub(crate) fn proof_single() {
        // Our polynomial: degree 15, 16 coefficients
        let coeffs = vec![1, 2, 3, 4, 7, 7, 7, 7, 13, 13, 13, 13, 13, 13, 13, 13];
        let poly_len = coeffs.len();
        let secrets_len = poly_len + 1;


        // Create the polynomial
        let mut p = new_poly(poly_len);
        for x in 0..poly_len{
                let temp = std::mem::replace(&mut p.coeffs[x], fr_from_uint64(coeffs[x]));
        }

        // Initialise the secrets and data structures
        let (s1,s2) = generate_trusted_setup(secrets_len).unwrap();
        let fs: FFTSettings = new_fft_settings(4);
        let mut ks: KZGSettings = new_kzg_settings(&s1, &s2, secrets_len as u64, fs);

        // Compute the proof for x = 25
        let x = fr_from_uint64(25);
        let commitment = commit_to_poly(&p, &mut ks).unwrap();
        let proof = compute_proof_single(&p, &x, &ks);
        let mut value = eval_poly(&p, &x);

        // Verify the proof that the (unknown) polynomial has y = value at x = 25
        assert!(check_proof_single(&commitment, &proof, &x, &value, &ks));

        // Change the value and check that the proof fails
        
        value = fr_add(&value, &fr_one);
        assert!(!check_proof_single(&commitment, &proof, &x, &value, &ks));
}

#[test]
pub(crate) fn proof_multi() {
        // Our polynomial: degree 15, 16 coefficients
        let coeffs = vec![1, 2, 3, 4, 7, 7, 7, 7, 13, 13, 13, 13, 13, 13, 13, 13];
        let poly_len = coeffs.len();
        let secrets_len = poly_len + 1;

        // Compute proof at 2^coset_scale points
        let coset_scale = 3;
        let coset_len = 1 << coset_scale;
        //     fr_t y[coset_len];
        let mut y: Vec<Fr> = Vec::new();


        // Create the polynomial
        let mut p = new_poly(poly_len);
        for x in 0..poly_len{
                let temp = std::mem::replace(&mut p.coeffs[x], fr_from_uint64(coeffs[x]));
        }

        // Initialise the secrets and data structures
        let (s1,s2) = generate_trusted_setup(secrets_len).unwrap();
        let fs1: FFTSettings = new_fft_settings(4);
        let mut ks1: KZGSettings = new_kzg_settings(&s1, &s2, secrets_len as u64, fs1);

        // Commit to the polynomial
        let commitment = commit_to_poly(&p, &mut ks1).unwrap();

        let fs2: FFTSettings = new_fft_settings(coset_len);
        let mut ks2: KZGSettings = new_kzg_settings(&s1, &s2, secrets_len as u64, fs2);

        // Compute proof at the points [x * root_i] 0 <= i < coset_len
        let x = fr_from_uint64(5431);
        let proof = compute_proof_multi(&p, &x, coset_len, &mut ks2);

        // y_i is the value of the polynomial at each x_i
        for i in 0..coset_len {
                let tmp = fr_mul(&x, &ks2.fs.expanded_roots_of_unity[i as usize]);
                y.push(eval_poly(&p, &tmp));
        }

        //     Verify the proof that the (unknown) polynomial has value y_i at x_i
        let result = check_proof_multi(&commitment, &proof, &x, &y, coset_len, &ks2);
        assert_eq!(result, true);

        // Change a value and check that the proof fails
        let temp = fr_add(&y[4], &fr_one);
        let temp = std::mem::replace(&mut y[4], temp);
        let result = check_proof_multi(&commitment, &proof, &x, &y, coset_len, &ks2);
        assert_eq!(result, false);
}

#[test]
pub(crate) fn commit_to_nil_poly() {
        let secrets_len = 16;

        // Initialise the (arbitrary) secrets and data structures
        let (s1,s2) = generate_trusted_setup(secrets_len).unwrap();
        let fs: FFTSettings = new_fft_settings(4);
        let mut ks: KZGSettings = new_kzg_settings(&s1, &s2, secrets_len as u64, fs);

        let a = new_poly(0);
        let result: P1 = commit_to_poly(&a, &mut ks).unwrap();
        assert!(result.equals(&g1_identity));
}


#[test]
pub(crate) fn commit_to_too_long_poly() {
        let secrets_len = 16;
        let poly_len = 32; // poly is longer than secrets!

        // Initialise the (arbitrary) secrets and data structures
        let (s1,s2) = generate_trusted_setup(secrets_len).unwrap();
        let fs: FFTSettings = new_fft_settings(4);
        let mut ks: KZGSettings = new_kzg_settings(&s1, &s2, secrets_len as u64, fs);

        let a = new_poly(poly_len);
        assert_eq!(commit_to_poly(&a, &mut ks), Err(TooLongPoly));
}


// Check that computing FFT and inverse FFT results in the starting data
#[test]
pub(crate) fn roundtrip_fft() {
    let size: usize = 12;

    let fft_settings = FFTSettings::from_scale(size).unwrap();

    let mut starting_data = vec![Fr::default(); fft_settings.max_width as usize];
    for i in 0..fft_settings.max_width {
        unsafe {
            blst_fr_from_uint64(&mut starting_data[i as usize], [i, 0, 0, 0].as_ptr());
        }
    }

    // Forward and inverse FFT
    let forward_result = fft_fr(&starting_data, false, &fft_settings).unwrap();
    let inverse_result = fft_fr(&forward_result, true, &fft_settings).unwrap();

    for i in 0..fft_settings.max_width {
        assert!(fr_are_equal(&starting_data[i as usize], &inverse_result[i as usize]));
    }
}

#[test]
pub(crate) fn stride_fft() {
    let size1: usize = 9;
    let size2: usize = 12;

    let width: usize = 1 << size1;

    let fft_settings1 = FFTSettings::from_scale(size1).unwrap();
    let fft_settings2 = FFTSettings::from_scale(size2).unwrap();

    let mut data = vec![Fr::default(); width];
    for i in 0..width {
        unsafe {
            blst_fr_from_uint64(&mut data[i], [i as u64, 0, 0, 0].as_ptr());
        }
    }

    let result1 = fft_fr(&data, false, &fft_settings1).unwrap();
    let result2 = fft_fr(&data, false, &fft_settings2).unwrap();

    for i in 0..width {
        println!("expected {:?}",result1[i]);
        println!("forward_result {:?}",result2[i]);
        assert!(fr_are_equal(&result1[i], &result2[i]));
    }
}


#[test]
pub(crate)fn inverse_fft() {
    let inv_fft_expected: [[u64; 4]; 16] =
        [
            [0x7fffffff80000008, 0xa9ded2017fff2dff, 0x199cec0404d0ec02, 0x39f6d3a994cebea4],
            [0xef296e7ffb8ca216, 0xd5b902cbcef9c1b6, 0xf06dfe5c7fca260d, 0x13993b7d05187205],
            [0xe930fdda2306c7d4, 0x40e02aff48e2b16b, 0x83a712d1dd818c8f, 0x5dbc603bc53c7a3a],
            [0xf9925986d0d25e90, 0xcdf85d0a339d7782, 0xee7a9a5f0410e423, 0x2e0d216170831056],
            [0x80007fff80000000, 0x1fe05202bb00adff, 0x6045d26b3fd26e6b, 0x39f6d3a994cebea4],
            [0x27325dd08ac4cee9, 0xcbb94f168ddacca9, 0x6843be68485784b1, 0x5a6faf9039451673],
            [0xe92ffdda2306c7d4, 0x54dd2afcd2dfb16b, 0xf6554603677e87be, 0x5dbc603bc53c7a39],
            [0x1cc772c9b57f126f, 0xfb73f4d33d3116dd, 0x4f9388c8d80abcf9, 0x3ffbc9abcdda7821],
            [0x7fffffff80000000, 0xa9ded2017fff2dff, 0x199cec0404d0ec02, 0x39f6d3a994cebea4],
            [0xe3388d354a80ed91, 0x5849af2fc2cd4521, 0xe3a64f3f31971b0b, 0x33f1dda75bc30526],
            [0x16d00224dcf9382c, 0xfee079062d1eaa93, 0x3ce49204a2235046, 0x163147176461030e],
            [0xd8cda22e753b3117, 0x880454ec72238f55, 0xcaf6199fc14a5353, 0x197df7c2f05866d4],
            [0x7fff7fff80000000, 0x33dd520044fdadff, 0xd2f4059cc9cf699a, 0x39f6d3a994cebea3],
            [0x066da6782f2da170, 0x85c546f8cc60e47c, 0x44bf3da90590f3e1, 0x45e085f1b91a6cf1],
            [0x16cf0224dcf9382c, 0x12dd7903b71baa93, 0xaf92c5362c204b76, 0x163147176461030d],
            [0x10d6917f04735dea, 0x7e04a13731049a48, 0x42cbd9ab89d7b1f7, 0x60546bd624850b42]
        ];

    let fft_settings = FFTSettings::from_scale(4).unwrap();

    let mut data = vec![Fr::default(); fft_settings.max_width as usize];
    for i in 0..fft_settings.max_width {
        unsafe {
            blst_fr_from_uint64(&mut data[i as usize], [i, 0, 0, 0].as_ptr());
        }
    }

    let forward_result = fft_fr(&data, true, &fft_settings).unwrap();

    assert_eq!(inv_fft_expected.len(), fft_settings.max_width as usize);
    
    for i in 0..inv_fft_expected.len() {
        let mut expected: Fr = Fr::default();
        unsafe {
            blst_fr_from_uint64(&mut expected, inv_fft_expected[i].as_ptr());
        }
        assert!(fr_are_equal(&expected, &forward_result[i]));
    }
}