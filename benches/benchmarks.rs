#[macro_use]
extern crate criterion;
use criterion::Criterion;

use rand::distr::{Distribution, Uniform};
use rand::rng;

#[allow(unused_must_use)]
fn bench_encode(c: &mut Criterion) {
    let mut rng = rng();
    // These coordinates cover London, approximately
    let between_lon = Uniform::new(-6.379880, 1.768960).unwrap();
    let between_lat = Uniform::new(49.871159, 55.811741).unwrap();
    let mut lon_vec: Vec<f64> = vec![];
    let mut lat_vec: Vec<f64> = vec![];
    (0..100000).for_each(|_| {
        lon_vec.push(between_lon.sample(&mut rng));
        lat_vec.push(between_lat.sample(&mut rng));
    });
    c.bench_function("bench encode: 100000 coordinates", move |b| {
        b.iter(|| {
            lonlat_bng::convert_to_bng_threaded_vec(&mut lon_vec, &mut lat_vec);
        })
    });
}

criterion_group!(benches, bench_encode,);
criterion_main!(benches);
