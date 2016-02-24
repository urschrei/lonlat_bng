// this is just a test binary that can be hooked up to Instruments
extern crate lonlat_bng;
use lonlat_bng::{convert_to_bng_threaded, Array};

extern crate rand;
use rand::distributions::{IndependentSample, Range};

extern crate libc;

fn main() {
    let num_coords = 1000000;
    let between_lon = Range::new(7.5600, 1.7800);
    let between_lat = Range::new(49.9600, 60.8400);
    let mut rng = rand::thread_rng();
    let lon_vec = vec![between_lon.ind_sample(&mut rng); num_coords];
    let lat_vec = vec![between_lat.ind_sample(&mut rng); num_coords];
    let lon_arr = Array {
        data: lon_vec.as_ptr() as *const libc::c_void,
        len: lon_vec.len() as libc::size_t,
    };
    let lat_arr = Array {
        data: lat_vec.as_ptr() as *const libc::c_void,
        len: lat_vec.len() as libc::size_t,
    };

    convert_to_bng_threaded(lon_arr, lat_arr);
    // convert_to_osgb36_threaded(lon_arr, lat_arr);
}
