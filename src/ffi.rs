use std::slice;

use libc::c_void;

#[repr(C)]
pub struct Array {
    pub data: *const c_void,
    pub len: libc::size_t,
}

/// A simple struct that can be leaked across the FFI boundary in lieu of an actual tuple
#[repr(C)]
pub struct ResultTuple {
    pub e: Array,
    pub n: Array,
}

use super::convert_epsg3857_to_wgs84_threaded_vec;
use super::convert_etrs89_to_ll_threaded_vec;
use super::convert_etrs89_to_osgb36_threaded_vec;
use super::convert_osgb36_to_etrs89_threaded_vec;
use super::convert_osgb36_to_ll_threaded_vec;
use super::convert_to_bng_threaded_vec;
use super::convert_to_etrs89_threaded_vec;
use super::convert_to_lonlat_threaded_vec;
use super::convert_to_osgb36_threaded_vec;

/// Free memory which Rust has allocated across the FFI boundary (f64 values)
///
/// # Examples
///
/// ```
/// # use libc;
/// let lon_vec: Vec<f64> = vec![-2.0183041005533306];
/// let lat_vec: Vec<f64> = vec![54.589097162646141];
/// let lon_arr = Array {
///     data: lon_vec.as_ptr() as *const libc::c_void,
///     len: lon_vec.len() as libc::size_t,
/// };
/// let lat_arr = Array {
///     data: lat_vec.as_ptr() as *const libc::c_void,
///     len: lat_vec.len() as libc::size_t,
/// };
/// let rtp = convert_to_bng_threaded(lon_arr, lat_arr);
/// drop_float_array(rtp.e, rtp.n);
/// ```
/// An example FFI implementation is available at [Convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py), specifically in the `_void_array_to_list` function.
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn drop_float_array(lons: Array, lats: Array) {
    if lons.data.is_null() {
        return;
    }
    if lats.data.is_null() {
        return;
    }
    let _: &mut [f64] = lons.into();
    let _: &mut [f64] = lats.into();
}

// Build an Array so it can be leaked across the FFI boundary
impl<'a> From<&'a mut [f64]> for Array {
    fn from(sl: &mut [f64]) -> Self {
        Array {
            data: sl.as_mut_ptr().cast::<libc::c_void>(),
            len: sl.len() as libc::size_t,
        }
    }
}

// Build &mut[f64] from an Array, so it can be dropped
impl<'a> From<Array> for &'a mut [f64] {
    fn from(arr: Array) -> Self {
        unsafe { slice::from_raw_parts_mut(arr.data as *mut f64, arr.len) }
    }
}

// Slice tuple to ResultTuple
impl<T> From<(T, T)> for ResultTuple
where
    T: Into<Array>,
{
    fn from(tup: (T, T)) -> Self {
        ResultTuple {
            e: tup.0.into(),
            n: tup.1.into(),
        }
    }
}

/// A threaded, FFI-compatible wrapper for `lonlat_bng::convert_osgb36`
///
/// # Examples
///
/// ```
/// use libc;
/// let lon_vec: Vec<f64> = vec![-2.0183041005533306,
///                              0.95511887434519682,
///                              0.44975855518383501,
///                              -0.096813621191803811,
///                              -0.36807065656416427,
///                              0.63486335458665621];
/// let lat_vec: Vec<f64> = vec![54.589097162646141,
///                              51.560873800587828,
///                              50.431429161121699,
///                              54.535021436247419,
///                              50.839059313135706,
///                              55.412189281234419];
/// let lon_arr = Array {
///     data: lon_vec.as_ptr() as *const libc::c_void,
///     len: lon_vec.len() as libc::size_t,
/// };
/// let lat_arr = Array {
///     data: lat_vec.as_ptr() as *const libc::c_void,
///     len: lat_vec.len() as libc::size_t,
/// };
/// let rtp = convert_to_bng_threaded(lon_arr, lat_arr);
/// ```
/// For an FFI implementation, see the code at [Convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_to_bng_threaded(longitudes: Array, latitudes: Array) -> ResultTuple {
    convert_to_bng_threaded_vec(longitudes.into(), latitudes.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
///
/// # Examples
///
/// See [`lonlat_bng::convert_to_bng_threaded`](fn.convert_to_bng_threaded.html)
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_to_lonlat_threaded(eastings: Array, northings: Array) -> ResultTuple {
    convert_to_lonlat_threaded_vec(eastings.into(), northings.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
///
/// # Examples
///
/// See [`lonlat_bng::convert_to_bng_threaded`](fn.convert_to_bng_threaded.html) for examples, substituting f64 vectors
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_to_osgb36_threaded(longitudes: Array, latitudes: Array) -> ResultTuple {
    convert_to_osgb36_threaded_vec(longitudes.into(), latitudes.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_etrs89`](fn.convert_etrs89.html)
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples, substituting f64 vectors
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_to_etrs89_threaded(longitudes: Array, latitudes: Array) -> ResultTuple {
    convert_to_etrs89_threaded_vec(longitudes.into(), latitudes.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_etrs89_to_osgb36_threaded(
    eastings: Array,
    northings: Array,
) -> ResultTuple {
    convert_etrs89_to_osgb36_threaded_vec(eastings.into(), northings.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_etrs89_to_ll`](fn.convert_etrs89_to_ll.html)
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_etrs89_to_ll_threaded(eastings: Array, northings: Array) -> ResultTuple {
    convert_etrs89_to_ll_threaded_vec(eastings.into(), northings.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_osgb36_to_ll_threaded(eastings: Array, northings: Array) -> ResultTuple {
    convert_osgb36_to_ll_threaded_vec(eastings.into(), northings.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36_to_etrs89`](fn.convert_osgb36_to_etrs89.html)
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_osgb36_to_etrs89_threaded(
    eastings: Array,
    northings: Array,
) -> ResultTuple {
    convert_osgb36_to_etrs89_threaded_vec(eastings.into(), northings.into()).into()
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_epsg3857_to_wgs84`](fn.convert_epsg3857_to_wgs84.html)
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_epsg3857_to_wgs84_threaded(x: Array, y: Array) -> ResultTuple {
    convert_epsg3857_to_wgs84_threaded_vec(x.into(), y.into()).into()
}
