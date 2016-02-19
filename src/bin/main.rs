// this is just a test binary that can be hooked up to Instruments
extern crate lonlat_bng;
use lonlat_bng::{convert_to_bng_threaded, convert_to_bng_threaded_nocopy, Array};

extern crate rand;
use rand::distributions::{IndependentSample, Range};

extern crate libc;

fn main() {
    let num_coords = 10000000;
    let between_lon = Range::new(-6.379880, 1.768960);
    let between_lat = Range::new(49.871159, 55.811741);
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

    convert_to_bng_threaded_nocopy(lon_arr, lat_arr);
    // convert_to_osgb36_threaded(lon_arr, lat_arr);
}
