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
        // these are the values from the OSTN15 NTv2 download
        // they exclude the two values which fall outside the bounding box
        // several values differ by one digit in the third decimal place (mm)
        let etrs89_lons_vec: &mut[f64] = &mut[
            -6.29977752,
            -5.2030461,
            -4.108645636,
            -1.297822772,
            -1.450514337,
            -3.551283492,
            1.444547304,
            -2.544076183,
            -0.119925572,
            -4.30852477,
            0.89724327,
            -2.154586144,
            -0.91248957,
            0.401535471,
            -1.197476559,
            -2.640493208,
            -4.289180698,
            -4.289177929,
            -3.040454907,
            -1.663791682,
            -4.634521682,
            -0.077731332,
            -4.388491181,
            -2.938277411,
            -1.616576852,
            -4.296490163,
            -3.294792193,
            -5.828366919,
            -2.048560307,
            -4.219263986,
            -8.578544561,
            -7.592555606,
            -6.260914555,
            -3.726310221,
            -3.214540011,
            -4.417576746,
            -5.827993398,
            -1.625169661,
            -1.274869104,
            -2.073828228
        ];
        let etrs89_lats_vec: &mut[f64] = &mut[
            49.92226394,
            49.96006138,
            50.43885826,
            50.57563665,
            50.93127938,
            51.4007822,
            51.37447026,
            51.42754743,
            51.48936565,
            51.85890896,
            51.89436637,
            52.25529382,
            52.25160951,
            52.75136687,
            52.96219109,
            53.3448028,
            53.41628516,
            53.41630925,
            53.77911026,
            53.8002152,
            54.08666318,
            54.11685144,
            54.32919541,
            54.8954234,
            54.97912274,
            55.85399953,
            55.92478266,
            57.00606696,
            57.13902519,
            57.48625001,
            57.81351838,
            58.21262247,
            58.51560361,
            58.58120461,
            59.03743871,
            59.09335035,
            59.09671617,
            59.53470794,
            59.85409914,
            60.13308092
        ];

        let osgb36_e_vec = [
            91492.146,
            170370.718,
            250359.811,
            449816.371,
            438710.92,
            292184.87,
            639821.835,
            362269.991,
            530624.974,
            241124.584,
            599445.59,
            389544.19,
            474335.969,
            562180.547,
            454002.834,
            357455.843,
            247958.971,
            247959.241,
            331534.564,
            422242.186,
            227778.33,
            525745.67,
            244780.636,
            339921.145,
            424639.355,
            256340.925,
            319188.434,
            167634.202,
            397160.491,
            267056.768,
            9587.909,
            71713.132,
            151968.652,
            299721.891,
            330398.323,
            261596.778,
            180862.461,
            421300.525,
            440725.073,
            395999.668
        ];
        let osgb36_n_vec = [
            11318.804,
            11572.405,
            62016.569,
            75335.861,
            114792.25,
            168003.465,
            169565.858,
            169978.69,
            178388.464,
            220332.641,
            225722.826,
            261912.153,
            262047.755,
            319784.995,
            340834.943,
            383290.436,
            393492.909,
            393495.583,
            431920.794,
            433818.701,
            468847.388,
            470703.214,
            495254.887,
            556034.761,
            565012.703,
            664697.269,
            670947.534,
            797067.144,
            805349.736,
            846176.972,
            899448.996,
            938516.404,
            966483.78,
            967202.992,
            1017347.016,
            1025447.602,
            1029604.114,
            1072147.239,
            1107878.448,
            1138728.951
        ];
        let e_arr = Array::from(etrs89_lons_vec);
        let n_arr = Array::from(etrs89_lats_vec);
        let (osgb36_eastings, osgb36_northings) = convert_to_bng_threaded(e_arr, n_arr);
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
        assert_eq!(651409.804, retval[0]);
    }

    #[test]
    fn test_threaded_lonlat_conversion_single() {
        // Caister Water Tower OSGB36 coords
        let easting_vec: &mut[f64] = &mut[651409.804];
        let northing_vec: &mut[f64] = &mut[313177.450];


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
        assert_eq!(651409.804, retval[0]);
    }

    #[test]
    fn test_threaded_etrs89_to_osgb36_conversion_single() {
        let e_vec: &mut[f64] = &mut[651307.003];
        let n_vec: &mut[f64]= &mut[313255.686];
        let e_arr = Array::from(e_vec);
        let n_arr = Array::from(n_vec);
        let (eastings, _) = convert_etrs89_to_osgb36_threaded(e_arr, n_arr);
        let retval: &mut [f64] = eastings.into();
        assert_eq!(651409.804, retval[0]);
    }

    #[test]
    fn test_threaded_osgb36_to_etrs89_conversion_single() {
        // Caister Water Tower OSGB36, see p21
        let e_vec: &mut[f64] = &mut[651409.804];
        let n_vec: &mut[f64] = &mut[313177.450];
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
        let e_vec: &mut[f64] = &mut[651409.804];
        let n_vec: &mut[f64] = &mut[313177.450];

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
