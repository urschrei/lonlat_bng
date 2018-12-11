#![feature(test)]

extern crate test;
use test::Bencher;

use lonlat_bng;
use rand::distributions::{IndependentSample, Range};

#[bench]
fn bench_threads(b: &mut Bencher) {
    let num_coords = 100000;
    let between_lon = Range::new(-6.379880, 1.768960);
    let between_lat = Range::new(49.871159, 55.811741);
    let mut rng = rand::thread_rng();
    let mut lon_vec = vec![between_lon.ind_sample(&mut rng); num_coords];
    let mut lat_vec = vec![between_lat.ind_sample(&mut rng); num_coords];
    b.iter(|| {
        lonlat_bng::convert_to_bng_threaded_vec(&mut lon_vec, &mut lat_vec);
    });
}
