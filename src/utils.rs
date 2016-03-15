#![doc(html_root_url = "https://urschrei.github.io/lonlat_bng/")]
//! This module provides utilities to the conversions module
use conversions::MIN_X_SHIFT;
use conversions::MIN_Y_SHIFT;
use conversions::MIN_Z_SHIFT;

use std;
use std::fmt;
use ostn02_phf::ostn02_lookup;

extern crate phf;
// include!("hexkeys.rs");

// fn helmert(lon_vec: [&f64], lat_vec: [&f64]) -> (Vec<f64>, Vec<f64>) {
//     let t_array = Vec3::new(TX, TY, TZ);
//     let params = Mat3::new(1. + s, RZS, -RYS, -RZS, 1. + s, RXS, RYS, -RXS, 1. + s);
//     // let zeros: Vec<&f64> = vec![&1.; lon_vec.len()];
//     let test = vec![1, 2, 3, 4, 5, 6];
//     // let mut combined: Vec<&f64> = lon_vec.extend(lat_vec.iter()).collect();
//     // combined.extend(zeros.iter().cloned());
//     // let tmat = DMat::from_row_vec(3, zeros.len(), &combined);
//     let tmat = DMat::from_row_vec(3, 2, &test);
//     // build input vector of x, y, and z columns
//     // let inp = DMat::from_row_vec(3, 1, vec![lon_vec, lat_vec, h_vec]);
//     (vec![1.], vec![2.])
// }

/// Bounds checking for input values
pub fn check<T>(to_check: T, bounds: (T, T)) -> Result<T, ()>
    where T: std::cmp::PartialOrd + fmt::Display + Copy
{
    match to_check {
        to_check if bounds.0 <= to_check && to_check <= bounds.1 => Ok(to_check),
        _ => Err(()),
    }
}

/// Round Easting, Northing, and Orthometric height to the nearest millimetre
pub fn round_to_nearest_mm(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let new_x = (x * 1000.).round() as f64 / 1000.;
    let new_y = (y * 1000.).round() as f64 / 1000.;
    let new_z = (z * 1000.).round() as f64 / 1000.;
    (new_x, new_y, new_z)
}

/// Round a float to nine decimal places
pub fn round_to_eight(x: f64, y: f64) -> (f64, f64) {
    let new_x = (x * 100000000.).round() as f64 / 100000000.;
    let new_y = (y * 100000000.).round() as f64 / 100000000.;
    (new_x, new_y)
}

/// Try to get OSTN02 shift parameters, and calculate offsets
pub fn get_ostn_ref(x: &i16, y: &i16) -> Result<(f64, f64, f64), ()> {
    let key = format!("{:03x}{:03x}", y, x);
    // let key = [hexkeys.get(y).unwrap().cloned(), hexkeys.get(x).unwrap().cloned()].concat();
    // Some or None, so convert to Result, which we can try!
    let result = try!(ostn02_lookup(&*key).ok_or(()));
    Ok((result.0 as f64 / 1000. + MIN_X_SHIFT,
        result.1 as f64 / 1000. + MIN_Y_SHIFT,
        result.2 as f64 / 1000. + MIN_Z_SHIFT))

}

// Input values must be valid ETRS89 grid references
// See p20 of the transformation user guide at
// https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html
pub fn ostn02_shifts(x: &f64, y: &f64) -> Result<(f64, f64, f64), ()> {
    let e_index = (*x / 1000.) as i16;
    let n_index = (*y / 1000.) as i16;

    // eastings and northings of the south-west corner of the cell
    let x0 = e_index as i32 * 1000;
    let y0 = n_index as i32 * 1000;

    // The easting, northing and geoid shifts for the four corners of the cell
    // any of these could be Err, so use try!

    // bottom-left grid intersection
    let s0: (f64, f64, f64) = try!(get_ostn_ref(&(e_index), &(n_index)));
    // bottom-right
    let s1: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 1), &(n_index + 0)));
    // top-left
    let s2: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 0), &(n_index + 1)));
    // top-right
    let s3: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 1), &(n_index + 1)));

    // offset within square
    let dx = x - (x0 as f64);
    let dy = y - (y0 as f64);

    let t = dx / 1000.;
    let u = dy / 1000.;

    let f0 = (1. - t as f64) * (1. - u as f64);
    let f1 = t as f64 * (1. - u as f64);
    let f2 = (1. - t as f64) * u as f64;
    let f3 = t as f64 * u as f64;

    // bilinear interpolation, to obtain the actual shifts
    let se = f0 * s0.0 + f1 * s1.0 + f2 * s2.0 + f3 * s3.0;
    let sn = f0 * s0.1 + f1 * s1.1 + f2 * s2.1 + f3 * s3.1;
    let sg = f0 * s0.2 + f1 * s1.2 + f2 * s2.2 + f3 * s3.2;

    Ok(round_to_nearest_mm(se, sn, sg))

}

#[cfg(test)]
mod tests {
    use super::get_ostn_ref;
    use super::ostn02_shifts;
    use super::check;
    // use super::HEXKEYS;
    use super::phf;
    include!("hexkeys.rs");


    #[test]
    // original coordinates are 651307.003, 313255.686
    fn test_ostn_hashmap_retrieval() {
        let eastings = 651;
        let northings = 313;
        let expected = (102.775, -78.244, 44.252);
        assert_eq!(expected, get_ostn_ref(&eastings, &northings).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_failed_ostn_hashmap_retrieval() {
        // we're try!ing this in the shift calculation function, so an Err is fine
        let eastings = 999;
        let northings = 999;
        let expected = (1., 1., 1.);
        assert_eq!(expected, get_ostn_ref(&eastings, &northings).unwrap());
    }

    #[test]
    fn test_ostn02_shift_incorporation() {
        // these are the input values and corrections on p20-21
        let eastings = 651307.003;
        let northings = 313255.686;
        let expected = (102.789, -78.238, 44.244);
        assert_eq!(expected, ostn02_shifts(&eastings, &northings).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_min_lon_extents() {
        let max_lon = 1.768960;
        let min_lon = -6.379880;
        // below min_lon
        check(&-6.379881, (&min_lon, &max_lon)).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_min_lat_extents() {
        let max_lat = 55.811741;
        let min_lat = 49.871159;
        // below min lat
        check(&49.871158, (&min_lat, &max_lat)).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_max_lon_extents() {
        let max_lon = 1.768960;
        let min_lon = -6.379880;
        // above max lon
        check(&1.768961, (&min_lon, &max_lon)).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_max_lat_extents() {
        let max_lat = 55.811741;
        let min_lat = 49.871159;
        // above max lat
        check(&55.811742, (&min_lat, &max_lat)).unwrap();
    }

    #[test]
    fn test_hexkeys() {
        let val = 711;
        assert_eq!(HEXKEYS[&val], "2c7");
    }

}
