// use criterion::{criterion_group, criterion_main, Criterion};
// use kzg::Poly;

// fn poly_division_in_finite_field(c: &mut Criterion) {
//     c.bench_function(
//         "poly_division_in_finite_field",
//         |b| b.iter(|| {
//             Poly::bench_divide_in_finite_field(8)
//         })
//     );
// }

// criterion_group!(benches, poly_division_in_finite_field);
// criterion_main!(benches);
