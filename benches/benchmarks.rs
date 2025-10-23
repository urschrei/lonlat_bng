#[macro_use]
extern crate criterion;
use criterion::Criterion;

use rand::distributions::Distribution;
use rand::distributions::Uniform;
use rand::thread_rng;

#[allow(unused_must_use)]
fn bench_encode(c: &mut Criterion) {
    let mut rng = thread_rng();
    // These coordinates cover London, approximately
    let between_lon = Uniform::from(-6.379880..1.768960);
    let between_lat = Uniform::from(49.871159..55.811741);
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
