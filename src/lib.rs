#![deny(
    clippy::cast_slice_from_raw_parts,
    clippy::cast_slice_different_sizes,
    invalid_null_arguments,
    clippy::ptr_as_ptr,
    clippy::transmute_ptr_to_ref
)]
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/urschrei/lonlat_bng/master/points.png",
    html_root_url = "https://urschrei.github.io/lonlat_bng/"
)]
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
use std::f64;

use rayon::prelude::*;

mod conversions;
mod ffi;
mod tests;
pub mod utils;

pub use crate::ffi::*;

pub use crate::conversions::convert_epsg3857_to_wgs84;
pub use crate::conversions::convert_etrs89;
pub use crate::conversions::convert_etrs89_to_ll;
pub use crate::conversions::convert_etrs89_to_osgb36;
pub use crate::conversions::convert_osgb36;
pub use crate::conversions::convert_osgb36_to_etrs89;
pub use crate::conversions::convert_osgb36_to_ll;

pub const NAN: f64 = f64::NAN;

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_bng_threaded_vec<'a>(
    longitudes: &'a mut [f64],
    latitudes: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_to_lonlat_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89`](fn.convert_etrs89.html)
pub fn convert_to_etrs89_threaded_vec<'a>(
    longitudes: &'a mut [f64],
    latitudes: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_osgb36_threaded_vec<'a>(
    longitudes: &'a mut [f64],
    latitudes: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
pub fn convert_etrs89_to_osgb36_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_etrs89_to_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_ll`](fn.convert_etrs8989_to_ll.html)
pub fn convert_etrs89_to_ll_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_etrs89_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_etrs89`](fn.convert_osgb36_to_etrs89.html)
pub fn convert_osgb36_to_etrs89_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_osgb36_to_ll_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_epsg3857_to_wgs84`](fn.convert_epsg3857_to_wgs84.html)
pub fn convert_epsg3857_to_wgs84_threaded_vec<'a>(
    x: &'a mut [f64],
    y: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(x, y, convert_epsg3857_to_wgs84)
}

// Generic function which applies conversion functions to vector or slice chunks within threads
// As opposed to the earlier convert_vec, we're directly modifying and returning the
// inputs here, at the cost of having to use lifetime annotations
fn convert_vec_direct<'a>(
    ex: &'a mut [f64],
    ny: &'a mut [f64],
    func: impl Fn(f64, f64) -> Result<(f64, f64), ()> + Sync,
) -> (&'a mut [f64], &'a mut [f64]) {
    ex.par_iter_mut().zip(ny.par_iter_mut()).for_each(|p| {
        match func(*p.0, *p.1) {
            // mutate values, or assign default error values
            Ok(res) => {
                *p.0 = res.0;
                *p.1 = res.1;
            }
            Err(_) => {
                *p.0 = NAN;
                *p.1 = NAN;
            }
        }
    });
    (ex, ny)
}
