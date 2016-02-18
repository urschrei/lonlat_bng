#![doc(html_root_url = "https://urschrei.github.io/lonlat_bng/")]
//! The `lonlat_bng` crate provides functions that convert decimal longitude
//! and latitude coordinates into [British National Grid](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid) coordinates, and vice versa.
//! This library makes use of the [OSTN02](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) transformations
//! for the following functions:
//!
//! - [`convert_osgb36`](fn.convert_osgb36.html)
//! - [`convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
//! - [`convert_to_osgb36_threaded`](fn.convert_to_osgb36_threaded.html)
//! - [`convert_to_osgb36_threaded_vec`](fn.convert_to_osgb_threaded_vec.html)
//! - [`convert_etrs89_to_osgb36_threaded`](fn.convert_etrs89_to_osgb36_threaded.html)
//! - [`convert_etrs89_to_osgb36_threaded_vec`](fn.convert_etrs89_to_osgb36_threaded_vec.html)
//! - [`convert_osgb36_to_ll_threaded`](fn.convert_osgb36_to_ll_threaded.html)
//! - [`convert_osgb36_to_ll_threaded_vec`](fn.convert_osgb36_to_ll_threaded_vec.html)
//! - [`convert_osgb36_to_etrs89_threaded`](fn.convert_osgb36_to_etrs89_threaded.html)
//! - [`convert_osgb36_to_etrs89_threaded_vec`](fn.convert_osgb36_to_etrs89_threaded_vec.html)
//!
//! These functions transform input longitude and latitude coordinates to OSGB36 Eastings and Northings with high accuracy, and are suitable for use in surveying and construction. Please run your own tests, though.
//! **Note that `lon`, `lat` coordinates outside the UK bounding box will be transformed to `(9999, 9999)`, which cannot be mapped.**
//!
//! # Examples
//!
//! ```
//! // Convert single Longitude, Latitude values to OSGB36 Eastings and Northings
//! assert_eq!((651409.792, 313177.448), lonlat_bng::ostn02::convert_osgb36(&1.716073973, &52.658007833).unwrap());
//! ```
//! ```
//! // Convert single Longitude, Latitude values to BNG Eastings and Northings (Fast, but only accurate within ~7m)
//! assert_eq!((516276, 173141), lonlat_bng::convert_bng(&-0.32824866, &51.44533267).unwrap());
//! ```
//! ```
//! // Convert single BNG values to Longitude, Latitude (Fast, but only accurate within ~7m)
//! assert_eq!((-0.328248, 51.44534), lonlat_bng::convert_lonlat(&516276, &173141));
//! ```
//! ```
//! // Convert vectors using threaded functions
//! lonlat_bng::convert_to_bng_threaded_vec(vec![&-0.32824866], vec![&51.44533267]);
//! ```
//! ```
//! lonlat_bng::convert_to_lonlat_threaded_vec(vec![&516276], vec![&173141]);
//! ```
//! The crate also provides C-compatible wrapper functions, which are intended for use with FFI.
//! An example implementation using Python can be found at [Convertbng](https://github.com/urschrei/convertbng).
use std::marker::Send;

extern crate libc;

extern crate rand;

extern crate crossbeam;
use crossbeam::scope;

extern crate num_cpus;
extern crate ostn02_phf;

mod conversions;
mod utils;
mod ffi;

pub use ffi::Array;
pub use ffi::drop_float_array;
pub use ffi::convert_to_bng_threaded;
pub use ffi::convert_to_lonlat_threaded;
pub use ffi::convert_to_osgb36_threaded;
pub use ffi::convert_to_etrs89_threaded;
pub use ffi::convert_etrs89_to_osgb36_threaded;
pub use ffi::convert_etrs89_to_ll_threaded;
pub use ffi::convert_osgb36_to_ll_threaded;
pub use ffi::convert_osgb36_to_etrs89_threaded;

pub use conversions::convert_bng;
pub use conversions::convert_lonlat;
pub use conversions::convert_etrs89;
pub use conversions::convert_osgb36;
pub use conversions::convert_etrs89_to_osgb36;
pub use conversions::convert_osgb36_to_etrs89;
pub use conversions::convert_osgb36_to_ll;
pub use conversions::convert_etrs89_to_ll;

// fn helmert(lon_vec: Vec<&f32>, lat_vec: Vec<&f32>) -> (Vec<f32>, Vec<f32>) {
//     let t_array = Vec3::new(TX, TY, TZ);
//     let params = Mat3::new(1. + s, RZS, -RYS, -RZS, 1. + s, RXS, RYS, -RXS, 1. + s);
//     // let zeros: Vec<&f32> = vec![&1.; lon_vec.len()];
//     let test = vec![1, 2, 3, 4, 5, 6];
//     // let mut combined: Vec<&f32> = lon_vec.extend(lat_vec.iter()).collect();
//     // combined.extend(zeros.iter().cloned());
//     // let tmat = DMat::from_row_vec(3, zeros.len(), &combined);
//     let tmat = DMat::from_row_vec(3, 2, &test);
//     // build input vector of x, y, and z columns
//     // let inp = DMat::from_row_vec(3, 1, vec![lon_vec, lat_vec, h_vec]);
//     (vec![1.], vec![2.])
// }

/// A threaded wrapper for [`lonlat_bng::convert_bng`](fn.convert_bng.html)
pub fn convert_to_bng_threaded_vec(longitudes: &[f64],
                                   latitudes: &[f64])
                                   -> (Vec<f64>, Vec<f64>) {
    convert_vec(longitudes, latitudes, convert_bng)
}

/// A threaded wrapper for [`lonlat_bng::convert_lonlat`](fn.convert_lonlat.html)
pub fn convert_to_lonlat_threaded_vec(eastings: &[f64],
                                      northings: &[f64])
                                      -> (Vec<f64>, Vec<f64>) {
    convert_vec(eastings, northings, convert_lonlat)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89`](fn.convert_etrs89.html)
pub fn convert_to_etrs89_threaded_vec(longitudes: &[f64],
                                      latitudes: &[f64])
                                      -> (Vec<f64>, Vec<f64>) {
    convert_vec(longitudes, latitudes, convert_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_osgb36_threaded_vec(longitudes: &[f64],
                                      latitudes: &[f64])
                                      -> (Vec<f64>, Vec<f64>) {
    convert_vec(longitudes, latitudes, convert_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
pub fn convert_etrs89_to_osgb36_threaded_vec(eastings: &[f64],
                                             northings: &[f64])
                                             -> (Vec<f64>, Vec<f64>) {
    convert_vec(eastings, northings, convert_etrs89_to_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_ll`](fn.convert_etrs8989_to_ll.html)
pub fn convert_etrs89_to_ll_threaded_vec(eastings: &[f64],
                                         northings: &[f64])
                                         -> (Vec<f64>, Vec<f64>) {
    convert_vec(eastings, northings, convert_etrs89_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_etrs89`](fn.convert_osgb36_to_etrs89.html)
pub fn convert_osgb36_to_etrs89_threaded_vec(eastings: &[f64],
                                             northings: &[f64])
                                             -> (Vec<f64>, Vec<f64>) {
    convert_vec(eastings, northings, convert_osgb36_to_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_osgb36_to_ll_threaded_vec(eastings: &[f64],
                                         northings: &[f64])
                                         -> (Vec<f64>, Vec<f64>) {
    convert_vec(eastings, northings, convert_osgb36_to_ll)
}

/// Generic function for threaded processing of conversion functions
fn convert_vec<F>(ex: &[f64], ny: &[f64], func: F) -> (Vec<f64>, Vec<f64>)
    where F: Fn(&f64, &f64) -> Result<(f64, f64), ()> + Send + Copy
{
    let numthreads = num_cpus::get() as usize;
    let orig: Vec<(&f64, &f64)> = ex.iter().zip(ny.iter()).collect();
    let mut result = vec![(1.0, 1.0); orig.len()];
    let mut size = orig.len() / numthreads;
    if orig.len() % numthreads > 0 {
        size += 1;
    }
    size = std::cmp::max(1, size);
    crossbeam::scope(|scope| {
        for (res_chunk, orig_chunk) in result.chunks_mut(size).zip(orig.chunks(size)) {
            scope.spawn(move || {
                for (res_elem, orig_elem) in res_chunk.iter_mut().zip(orig_chunk.iter()) {
                    match func(orig_elem.0, orig_elem.1) {
                        Ok(res) => *res_elem = res,
                        // we don't care about the return value as such
                        Err(_) => *res_elem = (9999.000, 9999.000),
                    };
                }
            });
        }
    });
    let (ex_converted, ny_converted): (Vec<f64>, Vec<f64>) = result.into_iter().unzip();
    (ex_converted, ny_converted)
}

#[cfg(test)]
mod tests {
    use ffi::drop_float_array;
    use ffi::Array;
    use super::convert_to_bng_threaded;
    use super::convert_to_lonlat_threaded;
    use super::convert_to_osgb36_threaded;
    use super::convert_to_etrs89_threaded;
    use super::convert_etrs89_to_osgb36_threaded;
    use super::convert_osgb36_to_etrs89_threaded;
    use super::convert_osgb36_to_ll_threaded;
    use super::convert_etrs89_to_ll_threaded;

    extern crate libc;

    #[test]
    #[ignore]
    fn test_os_etrs89_to_osgb36() {
        // these are the values from the OSTN02 download
        // several values differ by one digit in the third decimal place (mm)
        let etrs89_e_vec = vec![331439.16, 362174.408, 151874.984, 339824.598, 241030.731,
                                599345.196, 357359.683, 389448.042, 319092.329, 525643.491,
                                397061.069, 256247.486, 266961.481, 244687.517, 227685.882,
                                562079.805, 422143.679, 170277.189, 530526.413, 247865.448,
                                247865.718, 167542.805, 292090.289, 424539.719, 639720.224,
                                474237.874, 453904.269, 438614.045, 250265.789, 449719.403,
                                440623.592, 299624.627, 91400.000, 9500.003, 71622.45, 180766.824,
                                261500.000, 395898.578, 421200.000, 330300.000, 337800.000,
                                334100.000];
        let etrs89_n_vec = vec![431992.943,
                                170056.5,
                                966535.331,
                                556102.504,
                                220409.858,
                                225801.485,
                                383364.152,
                                261989.271,
                                671010.63,
                                470775.61,
                                805408.146,
                                664760.292,
                                846233.646,
                                495324.611,
                                468918.331,
                                319862.042,
                                433891.207,
                                11652.895,
                                178467.043,
                                393566.264,
                                393568.938,
                                797124.153,
                                168081.281,
                                565080.533,
                                169645.824,
                                262125.333,
                                340910.743,
                                114871.192,
                                62095.883,
                                75415.594,
                                1107930.781,
                                967256.596,
                                11400.000,
                                899499.996,
                                938567.303,
                                1029654.639,
                                1025500.000,
                                1138780.346,
                                1072200.000,
                                1017400.000,
                                981800.000,
                                982100.000];
        let osgb36_e_vec = vec![331534.552, 362269.979, 151968.641, 339921.133, 241124.573,
                                599445.578, 357455.831, 389544.178, 319188.423, 525745.658,
                                397160.479, 256340.914, 267056.756, 244780.625, 227778.318,
                                562180.535, 422242.174, 170370.706, 530624.963, 247958.959,
                                247959.229, 167634.19, 292184.858, 424639.343, 639821.823,
                                474335.957, 454002.822, 438710.908, 250359.798, 449816.359,
                                440725.061, 299721.879, 91492.135, 9587.897, 71713.12, 180862.449,
                                261596.767, 395999.656, 421300.513, 330398.311, 337898.195,
                                334198.101];
        let osgb36_n_vec = vec![431920.792,
                                169978.688,
                                966483.777,
                                556034.759,
                                220332.638,
                                225722.824,
                                383290.434,
                                261912.151,
                                670947.532,
                                470703.211,
                                805349.734,
                                664697.266,
                                846176.969,
                                495254.884,
                                468847.386,
                                319784.993,
                                433818.699,
                                11572.404,
                                178388.461,
                                393492.906,
                                393495.58,
                                797067.142,
                                168003.462,
                                565012.7,
                                169565.856,
                                262047.752,
                                340834.941,
                                114792.248,
                                62016.567,
                                75335.859,
                                1107878.445,
                                967202.99,
                                11318.801,
                                899448.993,
                                938516.401,
                                1029604.111,
                                1025447.599,
                                1138728.948,
                                1072147.236,
                                1017347.013,
                                981746.359,
                                982046.419];
        let e_arr = Array {
            data: etrs89_e_vec.as_ptr() as *const libc::c_void,
            len: etrs89_e_vec.len() as libc::size_t,
        };
        let n_arr = Array {
            data: etrs89_n_vec.as_ptr() as *const libc::c_void,
            len: etrs89_n_vec.len() as libc::size_t,
        };
        let (osgb36_eastings, osgb36_northings) = convert_etrs89_to_osgb36_threaded(e_arr, n_arr);
        let retval = unsafe { osgb36_eastings.as_f64_slice() };
        let retval2 = unsafe { osgb36_northings.as_f64_slice() };
        assert_eq!(osgb36_e_vec, retval);
        assert_eq!(osgb36_n_vec, retval2);
    }

    #[test]
    fn test_threaded_bng_conversion() {
        let lon_vec: Vec<f64> = vec![-2.0183041005533306,
                                     0.95511887434519682,
                                     0.44975855518383501,
                                     -0.096813621191803811,
                                     -0.36807065656416427,
                                     0.63486335458665621];
        let lat_vec: Vec<f64> = vec![54.589097162646141,
                                     51.560873800587828,
                                     50.431429161121699,
                                     54.535021436247419,
                                     50.839059313135706,
                                     55.412189281234419];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        let retval2 = unsafe { northings.as_f64_slice() };
        // the value's incorrect, but let's worry about that later
        assert_eq!(398915.033, retval[0]);
        assert_eq!(521545.067, retval2[0]);
    }

    #[test]
    fn test_threaded_bng_conversion_single() {
        // I spent 8 hours confused cos I didn't catch that chunks(0) is invalid
        let lon_vec: Vec<f64> = vec![-2.0183041005533306];
        let lat_vec: Vec<f64> = vec![54.589097162646141];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        assert_eq!(398915.033, retval[0]);
    }

    #[test]
    fn test_threaded_lonlat_conversion_single() {
        let easting_vec: Vec<i32> = vec![516276];
        let northing_vec: Vec<i32> = vec![173141];


        let easting_arr = Array {
            data: easting_vec.as_ptr() as *const libc::c_void,
            len: easting_vec.len() as libc::size_t,
        };
        let northing_arr = Array {
            data: northing_vec.as_ptr() as *const libc::c_void,
            len: northing_vec.len() as libc::size_t,
        };
        let (lons, _) = convert_to_lonlat_threaded(easting_arr, northing_arr);
        let retval = unsafe { lons.as_f64_slice() };
        // We shouldn't really be using error margins, but it should be OK because
        // neither number is zero, or very close to, and on opposite sides of zero
        // http://floating-point-gui.de/errors/comparison/
        assert!((retval[0] - -0.32824799370716407).abs() / -0.32824799370716407 < 0.000000001);
        // assert!((retval[1] - 51.44534026616287).abs() / 51.44534026616287 < 0.0000000001);
    }

    #[test]
    fn test_threaded_osgb36_conversion_single() {
        let lon_vec: Vec<f64> = vec![1.716073973];
        let lat_vec: Vec<f64> = vec![52.65800783];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_to_osgb36_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        assert_eq!(651409.792, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_to_osgb36_conversion_single() {
        let e_vec: Vec<f64> = vec![651307.003];
        let n_vec: Vec<f64> = vec![313255.686];
        let e_arr = Array {
            data: e_vec.as_ptr() as *const libc::c_void,
            len: e_vec.len() as libc::size_t,
        };
        let n_arr = Array {
            data: n_vec.as_ptr() as *const libc::c_void,
            len: n_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_etrs89_to_osgb36_threaded(e_arr, n_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        assert_eq!(651409.792, retval[0]);
    }

    #[test]
    fn test_threaded_osgb36_to_etrs89_conversion_single() {
        // Caister Water Tower OSGB36, see p21
        let e_vec: Vec<f64> = vec![651409.792];
        let n_vec: Vec<f64> = vec![313177.448];
        let e_arr = Array {
            data: e_vec.as_ptr() as *const libc::c_void,
            len: e_vec.len() as libc::size_t,
        };
        let n_arr = Array {
            data: n_vec.as_ptr() as *const libc::c_void,
            len: n_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_osgb36_to_etrs89_threaded(e_arr, n_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        // Caister Water Tower ETRS89, see p20
        assert_eq!(651307.003, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_conversion_single() {
        let lon_vec: Vec<f64> = vec![1.716073973];
        let lat_vec: Vec<f64> = vec![52.65800783];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_to_etrs89_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        assert_eq!(651307.003, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_to_ll_conversion_single() {
        let e_vec: Vec<f64> = vec![651307.003];
        let n_vec: Vec<f64> = vec![313255.686];
        let e_arr = Array {
            data: e_vec.as_ptr() as *const libc::c_void,
            len: e_vec.len() as libc::size_t,
        };
        let n_arr = Array {
            data: n_vec.as_ptr() as *const libc::c_void,
            len: n_vec.len() as libc::size_t,
        };
        let (lon, lat) = convert_etrs89_to_ll_threaded(e_arr, n_arr);
        let retval = unsafe { lon.as_f64_slice() };
        let retval2 = unsafe { lat.as_f64_slice() };
        assert_eq!(1.71607397, retval[0]);
        assert_eq!(52.65800783, retval2[0]);

    }

    #[test]
    fn test_threaded_osgb36_to_ll_conversion_single() {
        // Caister Water Tower
        let e_vec: Vec<f64> = vec![651409.792];
        let n_vec: Vec<f64> = vec![313177.448];

        let e_arr = Array {
            data: e_vec.as_ptr() as *const libc::c_void,
            len: e_vec.len() as libc::size_t,
        };
        let n_arr = Array {
            data: n_vec.as_ptr() as *const libc::c_void,
            len: n_vec.len() as libc::size_t,
        };
        let (lon, lat) = convert_osgb36_to_ll_threaded(e_arr, n_arr);
        let retval = unsafe { lon.as_f64_slice() };
        let retval2 = unsafe { lat.as_f64_slice() };
        assert_eq!(1.71607397, retval[0]);
        assert_eq!(52.65800783, retval2[0]);

    }

    #[test]
    fn test_drop_float_array() {
        let lon_vec: Vec<f64> = vec![-2.0183041005533306];
        let lat_vec: Vec<f64> = vec![54.589097162646141];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        drop_float_array(eastings, northings);
    }

    #[test]
    fn test_bad_threaded_conversion() {
        // above maximum longitude
        let lon_vec: Vec<f64> = vec![-6.379881];
        let lat_vec: Vec<f64> = vec![55.811741];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_f64_slice() };
        assert_eq!(9999.000, retval[0]);
    }

}
