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

mod conversions;
mod ffi;
mod tests;
mod threaded;
pub(crate) mod utils;

pub use crate::ffi::*;

pub use crate::conversions::convert_epsg3857_to_wgs84;
pub use crate::conversions::convert_etrs89;
pub use crate::conversions::convert_etrs89_to_ll;
pub use crate::conversions::convert_etrs89_to_osgb36;
pub use crate::conversions::convert_osgb36;
pub use crate::conversions::convert_osgb36_to_etrs89;
pub use crate::conversions::convert_osgb36_to_ll;

pub use crate::threaded::convert_epsg3857_to_wgs84_threaded_vec;
pub use crate::threaded::convert_etrs89_to_ll_threaded_vec;
pub use crate::threaded::convert_etrs89_to_osgb36_threaded_vec;
pub use crate::threaded::convert_osgb36_to_etrs89_threaded_vec;
pub use crate::threaded::convert_osgb36_to_ll_threaded_vec;
pub use crate::threaded::convert_to_bng_threaded_vec;
pub use crate::threaded::convert_to_etrs89_threaded_vec;
pub use crate::threaded::convert_to_lonlat_threaded_vec;
pub use crate::threaded::convert_to_osgb36_threaded_vec;

pub const NAN: f64 = f64::NAN;
