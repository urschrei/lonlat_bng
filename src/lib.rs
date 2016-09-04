#![doc(html_logo_url = "https://raw.githubusercontent.com/urschrei/lonlat_bng/master/points.png",
       html_root_url = "https://urschrei.github.io/lonlat_bng/")]
//! The `lonlat_bng` crate provides functions that convert decimal (WGS84 / ETRS89) longitude
//! and latitude coordinates into [British National Grid](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid) coordinates, and vice versa.
//! This library makes use of the [OSTN02](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) transformations
//!
//! These functions transform input longitude and latitude coordinates to OSGB36 Eastings and Northings with high accuracy, and are suitable for use in surveying and construction. Please run your own tests, though.
//! **Note that `lon`, `lat` coordinates outside the UK bounding box will be transformed to `(NAN, NAN)`, which cannot be mapped.**
//!
//! # Examples
//!
//! ```
//! // Convert single Longitude, Latitude values to OSGB36 Eastings and Northings
//! assert_eq!((651409.792, 313177.448), lonlat_bng::convert_osgb36(&1.716073973, &52.658007833).unwrap());
//! ```
//! ```
//! // Convert vectors or slices using multi-threaded functions
//! lonlat_bng::convert_to_osgb36_threaded_vec(vec![&-0.32824866], vec![&51.44533267]);
//! lonlat_bng::convert_osgb36_to_lonlat_threaded_vec(vec![&516276], vec![&173141]);
//! ```
//! The crate also provides C-compatible wrapper functions which are intended for use with FFI.
//!
//! **An example FFI implementation using Python can be found at [Convertbng](https://github.com/urschrei/convertbng)**.
//!
//! # Note
//!
//! FFI implementations **must** implement [`lonlat_bng::drop_float_array`](fn.drop_float_array.html), in order to free memory which has been allocated across the FFI boundary. Failure to do so may result in memory leaks.

use std::marker::Send;

extern crate rand;
extern crate phf;
extern crate ostn15_phf;
extern crate crossbeam;
use crossbeam::scope;
extern crate num_cpus;

mod conversions;
pub mod utils;
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
pub use ffi::convert_epsg3857_to_wgs84_threaded;

pub use conversions::convert_etrs89;
pub use conversions::convert_osgb36;
pub use conversions::convert_etrs89_to_osgb36;
pub use conversions::convert_osgb36_to_etrs89;
pub use conversions::convert_osgb36_to_ll;
pub use conversions::convert_etrs89_to_ll;
pub use conversions::convert_epsg3857_to_wgs84;

use std::f64;
pub const NAN: f64 = f64::NAN;

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_bng_threaded_vec<'a>(longitudes: &'a mut [f64],
                                       latitudes: &'a mut [f64])
                                       -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_to_lonlat_threaded_vec<'a>(eastings: &'a mut [f64],
                                          northings: &'a mut [f64])
                                          -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89`](fn.convert_etrs89.html)
pub fn convert_to_etrs89_threaded_vec<'a>(longitudes: &'a mut [f64],
                                          latitudes: &'a mut [f64])
                                          -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_osgb36_threaded_vec<'a>(longitudes: &'a mut [f64],
                                          latitudes: &'a mut [f64])
                                          -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
pub fn convert_etrs89_to_osgb36_threaded_vec<'a>(eastings: &'a mut [f64],
                                                 northings: &'a mut [f64])
                                                 -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_etrs89_to_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_ll`](fn.convert_etrs8989_to_ll.html)
pub fn convert_etrs89_to_ll_threaded_vec<'a>(eastings: &'a mut [f64],
                                             northings: &'a mut [f64])
                                             -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_etrs89_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_etrs89`](fn.convert_osgb36_to_etrs89.html)
pub fn convert_osgb36_to_etrs89_threaded_vec<'a>(eastings: &'a mut [f64],
                                                 northings: &'a mut [f64])
                                                 -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_osgb36_to_ll_threaded_vec<'a>(eastings: &'a mut [f64],
                                             northings: &'a mut [f64])
                                             -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_epsg3857_to_wgs84`](fn.convert_epsg3857_to_wgs84.html)
pub fn convert_epsg3857_to_wgs84_threaded_vec<'a>(x: &'a mut [f64],
                                                  y: &'a mut [f64])
                                                  -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(x, y, convert_epsg3857_to_wgs84)
}

// Generic function which applies conversion functions to vector or slice chunks within threads
// As opposed to the earlier convert_vec, we're directly modifying and returning the
// input vectors here, at the cost of having to use lifetime annotations
fn convert_vec_direct<'a, F>(ex: &'a mut [f64],
                             ny: &'a mut [f64],
                             func: F)
                             -> (&'a mut [f64], &'a mut [f64])
    where F: Fn(&f64, &f64) -> Result<(f64, f64), ()> + Send + Copy + Sync
{
    let numthreads = num_cpus::get() as usize;
    let mut size = ex.len() / numthreads;
    if ex.len() % numthreads > 0 {
        size += 1;
    }
    size = std::cmp::max(1, size);
    crossbeam::scope(|scope| {
        // chunks_mut returns chunks of "size"
        // e.g. 20 items / 4 (numthreads) = 5, resulting in 4 chunks, 1 per thread
        for (ex_chunk, ny_chunk) in ex.chunks_mut(size).zip(ny.chunks_mut(size)) {
            scope.spawn(move || {
                for (ex_elem, ny_elem) in ex_chunk.iter_mut().zip(ny_chunk.iter_mut()) {
                    match func(ex_elem, ny_elem) {
                        // mutate values, or assign default error values
                        Ok(res) => {
                            *ex_elem = res.0;
                            *ny_elem = res.1;
                        }
                        Err(_) => {
                            *ex_elem = NAN;
                            *ny_elem = NAN;
                        }
                    };
                }
            });
        }
    });
    (ex, ny)
}

#[cfg(test)]
mod tests {
    use super::*;
    use conversions::convert_bng;
    use super::convert_vec_direct;

    extern crate libc;
    use std::ptr;

    #[test]
    // Test Google/Bing Maps to WGS84 conversion
    fn test_epsg3857() {
        let x: &mut[f64] = &mut[-626172.1357121646];
        let y: &mut[f64] = &mut[6887893.4928337997];
        let x_arr = Array::from(x);
        let y_arr = Array::from(y);
        let (lon, lat) = convert_epsg3857_to_wgs84_threaded(x_arr, y_arr);
        let retval: &mut [f64] = lon.into();
        let retval2: &mut [f64] = lat.into();
        let expected = (-5.625000000783013, 52.48278022732355);
        assert_eq!(expected, (retval[0], retval2[0]));
    }

    #[test]
    // this test verifies that we aren't mangling memory inside our threads
    fn test_threading() {
        let mut lon = [-2.0183041005533306,
                       0.95511887434519682,
                       0.44975855518383501,
                       -0.096813621191803811,
                       -0.36807065656416427,
                       0.63486335458665621];
        let mut lat = [54.589097162646141,
                       51.560873800587828,
                       50.431429161121699,
                       54.535021436247419,
                       50.839059313135706,
                       55.412189281234419];
        convert_vec_direct(&mut lon, &mut lat, convert_bng);
        assert_eq!(398915.033, lon[0]);
        assert_eq!(521545.067, lat[0]);
    }

    #[test]
    fn test_os_etrs89_to_osgb36() {
        // these are the values from the OSTN02 download
        // they exclude the two values which fall outside the bounding box
        // several values differ by one digit in the third decimal place (mm)
        let etrs89_e_vec: &mut[f64] =
            &mut[331439.16, 362174.408, 151874.984, 339824.598, 241030.731, 599345.196,
                 357359.683, 389448.042, 319092.329, 525643.491, 397061.069, 256247.486,
                 266961.481, 244687.517, 227685.882, 562079.805, 422143.679, 170277.189,
                 530526.413, 247865.448, 247865.718, 167542.805, 292090.289, 424539.719,
                 639720.224, 474237.874, 453904.269, 438614.045, 250265.789, 449719.403,
                 440623.592, 299624.627, 91400.000, 9500.003, 71622.45, 180766.824, 261500.000,
                 395898.578, 421200.000, 330300.000, 337800.000, 334100.000];
        let etrs89_n_vec: &mut[f64] = &mut[431992.943,
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
        let osgb36_e_vec =
            [331534.552, 362269.979, 151968.641, 339921.133, 241124.573, 599445.578,
                 357455.831, 389544.178, 319188.423, 525745.658, 397160.479, 256340.914,
                 267056.756, 244780.625, 227778.318, 562180.535, 422242.174, 170370.706,
                 530624.963, 247958.959, 247959.229, 167634.19, 292184.858, 424639.343,
                 639821.823, 474335.957, 454002.822, 438710.908, 250359.798, 449816.359,
                 440725.061, 299721.879, 91492.135, 9587.897, 71713.12, 180862.449, 261596.767,
                 395999.656, 421300.513, 330398.311, 337898.195, 334198.101];
        let osgb36_n_vec = [431920.792,
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
        let e_arr = Array::from(etrs89_e_vec);
        let n_arr = Array::from(etrs89_n_vec);
        let (osgb36_eastings, osgb36_northings) = convert_etrs89_to_osgb36_threaded(e_arr, n_arr);
        let retval: &mut [f64] = osgb36_eastings.into();
        let retval2: &mut [f64] = osgb36_northings.into();
        // We're using absolute error margins here, but it should be OK
        // test eastings
        for (expect, result) in osgb36_e_vec.iter().zip(retval.iter()) {
            assert!(((expect - result) / result).abs() < 0.0000001)
        }
        // test northings
        for (expect, result) in osgb36_n_vec.iter().zip(retval2.iter()) {
            assert!(((expect - result) / result).abs() < 0.0000001)
        }
        // This will fail because 11 results differ by .001 (so 1 mm)
        // assert_eq!(osgb36_e_vec, retval);
        // assert_eq!(osgb36_n_vec, retval2);
    }

    #[test]
    #[should_panic]
    // This should panic because the OSTN02 test suite specifies it as outside the polygon: Outside#1
    // Lon: 4, 51, 3.503091
    // Lat: 53,20, 49.312599
    fn test_ostn_invalid_outside_1() {
        let _ = convert_osgb36(&4.850973, &53.347031).unwrap();
    }

    #[test]
    #[should_panic]
    // This should panic because the OSTN02 test suite specifies it as outside the polygon: Outside#2
    // Lon: 2, 22, 31.048596
    // Lat: 56, 10, 31.115299
    fn test_ostn_invalid_outside_2() {
        let _ = convert_osgb36(&2.375291, &56.17531).unwrap();

    }

    #[test]
    fn test_threaded_bng_conversion_single() {
        let lon_vec: &mut[f64] = &mut[1.716073973];
        let lat_vec: &mut[f64] = &mut[52.65800783];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (eastings, _) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval: &mut [f64] = eastings.into();
        assert_eq!(651409.792, retval[0]);
    }

    #[test]
    fn test_threaded_lonlat_conversion_single() {
        // Caister Water Tower OSGB36 coords
        let easting_vec: &mut[f64] = &mut[651409.792];
        let northing_vec: &mut[f64] = &mut[313177.448];


        let easting_arr = Array::from(easting_vec);
        let northing_arr = Array::from(northing_vec);
        let (lons, _) = convert_to_lonlat_threaded(easting_arr, northing_arr);
        let retval: &mut [f64] = lons.into();
        // We shouldn't really be using error margins, but it should be OK because
        // neither number is zero, or very close to, and on opposite sides of zero
        // http://floating-point-gui.de/errors/comparison/
        assert_eq!(1.71607397, retval[0]);
    }

    #[test]
    fn test_threaded_osgb36_conversion_single() {
        let lon_vec: &mut[f64] = &mut[1.716073973];
        let lat_vec: &mut[f64] = &mut[52.65800783];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (eastings, _) = convert_to_osgb36_threaded(lon_arr, lat_arr);
        let retval: &mut [f64] = eastings.into();
        assert_eq!(651409.792, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_to_osgb36_conversion_single() {
        let e_vec: &mut[f64] = &mut[651307.003];
        let n_vec: &mut[f64]= &mut[313255.686];
        let e_arr = Array::from(e_vec);
        let n_arr = Array::from(n_vec);
        let (eastings, _) = convert_etrs89_to_osgb36_threaded(e_arr, n_arr);
        let retval: &mut [f64] = eastings.into();
        assert_eq!(651409.792, retval[0]);
    }

    #[test]
    fn test_threaded_osgb36_to_etrs89_conversion_single() {
        // Caister Water Tower OSGB36, see p21
        let e_vec: &mut[f64] = &mut[651409.792];
        let n_vec: &mut[f64] = &mut[313177.448];
        let e_arr = Array::from(e_vec);
        let n_arr = Array::from(n_vec);
        let (eastings, _) = convert_osgb36_to_etrs89_threaded(e_arr, n_arr);
        let retval: &mut [f64] = eastings.into();
        // Caister Water Tower ETRS89, see p20
        assert_eq!(651307.003, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_conversion_single() {
        let lon_vec: &mut[f64] = &mut[1.716073973];
        let lat_vec: &mut[f64] = &mut[52.65800783];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (eastings, _) = convert_to_etrs89_threaded(lon_arr, lat_arr);
        let retval: &mut [f64] = eastings.into();
        assert_eq!(651307.003, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_to_ll_conversion_single() {
        let e_vec: &mut[f64] = &mut[651307.003];
        let n_vec: &mut[f64] = &mut[313255.686];
        let e_arr = Array::from(e_vec);
        let n_arr = Array::from(n_vec);
        let (lon, lat) = convert_etrs89_to_ll_threaded(e_arr, n_arr);
        let retval: &mut [f64] = lon.into();
        let retval2: &mut [f64] = lat.into();
        assert_eq!(1.71607397, retval[0]);
        assert_eq!(52.65800783, retval2[0]);

    }

    #[test]
    fn test_threaded_osgb36_to_ll_conversion_single() {
        // Caister Water Tower
        let e_vec: &mut[f64] = &mut[651409.792];
        let n_vec: &mut[f64] = &mut[313177.448];

        let e_arr = Array::from(e_vec);
        let n_arr = Array::from(n_vec);
        let (lon, lat) = convert_osgb36_to_ll_threaded(e_arr, n_arr);
        let retval: &mut [f64] = lon.into();
        let retval2: &mut [f64] = lat.into();
        assert_eq!(1.71607397, retval[0]);
        assert_eq!(52.65800783, retval2[0]);

    }

    #[test]
    fn test_drop_float_array() {
        let lon_vec: &mut[f64] = &mut[-2.0183041005533306];
        let lat_vec: &mut[f64] = &mut[54.589097162646141];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        drop_float_array(eastings, northings);
    }

    #[test]
    fn test_empty_lon_array() {
        let lon_vec: &mut[f64] = &mut[];
        let lat_vec: &mut[f64] = &mut[];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (mut eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        eastings.data = ptr::null();
        drop_float_array(eastings, northings);
    }

    #[test]
    fn test_empty_lat_array() {
        let lon_vec: &mut[f64] = &mut[];
        let lat_vec: &mut[f64] = &mut[];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (eastings, mut northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        northings.data = ptr::null();
        drop_float_array(eastings, northings);
    }

    #[test]
    fn test_bad_threaded_conversion() {
        // above maximum longitude
        let lon_vec: &mut[f64] = &mut[1.85];
        let lat_vec: &mut[f64] = &mut[55.811741];
        let lon_arr = Array::from(lon_vec);
        let lat_arr = Array::from(lat_vec);
        let (eastings, _) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval: &mut [f64] = eastings.into();
        assert!(retval[0].is_nan());
    }

}
