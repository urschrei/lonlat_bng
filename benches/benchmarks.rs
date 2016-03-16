#![feature(test)]

extern crate test;
use test::Bencher;
extern crate lonlat_bng;
use lonlat_bng::convert_to_bng_threaded_vec;
use lonlat_bng::utils::get_ostn_ref;

extern crate rand;
use rand::distributions::{IndependentSample, Range};

#[bench]
fn bench_threads(b: &mut Bencher) {
    let num_coords = 100000;
    let between_lon = Range::new(-6.379880, 1.768960);
    let between_lat = Range::new(49.871159, 55.811741);
    let mut rng = rand::thread_rng();
    let mut lon_vec = vec![between_lon.ind_sample(&mut rng); num_coords];
    let mut lat_vec = vec![between_lat.ind_sample(&mut rng); num_coords];
    b.iter(||{
        convert_to_bng_threaded_vec(
            &mut lon_vec,
            &mut lat_vec
        );
    });
}

#[bench]
fn bench_hex_lookup(b: &mut Bencher) {
    let eastings = 651;
    let northings = 313;
    b.iter(||{
        get_ostn_ref(&eastings, &northings);
    });
}
