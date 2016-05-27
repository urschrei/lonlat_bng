use std::mem;
use std::slice;

extern crate libc;
use self::libc::c_void;

#[repr(C)]
pub struct Array {
    pub data: *const c_void,
    pub len: libc::size_t,
}

use super::convert_to_bng_threaded_vec;
use super::convert_to_lonlat_threaded_vec;
use super::convert_to_osgb36_threaded_vec;
use super::convert_to_etrs89_threaded_vec;
use super::convert_etrs89_to_osgb36_threaded_vec;
use super::convert_etrs89_to_ll_threaded_vec;
use super::convert_osgb36_to_ll_threaded_vec;
use super::convert_osgb36_to_etrs89_threaded_vec;
use super::convert_epsg3857_to_wgs84_threaded_vec;

/// Free memory which Rust has allocated across the FFI boundary (f64 values)
///
/// # Examples
///
/// ```
/// # extern crate libc;
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
/// let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
/// drop_float_array(eastings, northings);
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

// Convert Array to mutable [f64] slice
impl Array {
    pub unsafe fn as_f64_slice(&self) -> &mut [f64] {
        assert!(!self.data.is_null());
        slice::from_raw_parts_mut(self.data as *mut f64, self.len as usize)
    }
}

// Build an Array from &mut[T], so it can be leaked across the FFI boundary
impl<'a, T> From<&'a mut [T]> for Array {
    fn from(sl: &mut [T]) -> Self {
        let array = Array {
            data: sl.as_ptr() as *const libc::c_void,
            len: sl.len() as libc::size_t,
        };
        mem::forget(sl);
        array
    }
}

// Build &mut[f64] from an Array, so it can be dropped
impl<'a> From<Array> for &'a mut [f64] {
    fn from(arr: Array) -> Self {
        unsafe { slice::from_raw_parts_mut(arr.data as *mut f64, arr.len) }
    }
}

// Build an Array from a Vec, so it can be leaked across the FFI boundary
impl<T> From<Vec<T>> for Array {
    fn from(vec: Vec<T>) -> Self {
        let array = Array {
            data: vec.as_ptr() as *const libc::c_void,
            len: vec.len() as libc::size_t,
        };
        mem::forget(vec);
        array
    }
}

// Build a Vec from an Array, so it can be dropped
impl From<Array> for Vec<f64> {
    fn from(arr: Array) -> Self {
        unsafe { Vec::from_raw_parts(arr.data as *mut f64, arr.len, arr.len) }
    }
}

/// A threaded, FFI-compatible wrapper for `lonlat_bng::convert_osgb36`
///
/// # Examples
///
/// ```
/// extern crate libc;
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
/// let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
/// ```
/// For an FFI implementation, see the code at [Convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data
#[no_mangle]
pub extern "C" fn convert_to_bng_threaded(longitudes: Array, latitudes: Array) -> (Array, Array) {
    let mut longitudes_v: &mut [f64] = longitudes.into();
    let mut latitudes_v: &mut [f64] = latitudes.into();
    convert_to_bng_threaded_vec(&mut longitudes_v, &mut latitudes_v);
    (Array::from(longitudes_v), Array::from(latitudes_v))
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
pub extern "C" fn convert_to_lonlat_threaded(eastings: Array, northings: Array) -> (Array, Array) {
    let mut eastings_vec: &mut [f64] = eastings.into();
    let mut northings_vec: &mut [f64] = northings.into();
    convert_to_lonlat_threaded_vec(&mut eastings_vec, &mut northings_vec);
    (Array::from(eastings_vec), Array::from(northings_vec))
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
pub extern "C" fn convert_to_osgb36_threaded(longitudes: Array,
                                             latitudes: Array)
                                             -> (Array, Array) {
    let mut longitudes_v: &mut [f64] = longitudes.into();
    let mut latitudes_v: &mut [f64] = latitudes.into();
    convert_to_osgb36_threaded_vec(&mut longitudes_v, &mut latitudes_v);
    (Array::from(longitudes_v), Array::from(latitudes_v))
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
pub extern "C" fn convert_to_etrs89_threaded(longitudes: Array,
                                             latitudes: Array)
                                             -> (Array, Array) {
    let mut longitudes_v: &mut [f64] = longitudes.into();
    let mut latitudes_v: &mut [f64] = latitudes.into();
    convert_to_etrs89_threaded_vec(&mut longitudes_v, &mut latitudes_v);
    (Array::from(longitudes_v), Array::from(latitudes_v))
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
pub extern "C" fn convert_etrs89_to_osgb36_threaded(eastings: Array,
                                                    northings: Array)
                                                    -> (Array, Array) {
    let mut eastings_vec: &mut [f64] = eastings.into();
    let mut northings_vec: &mut [f64] = northings.into();
    convert_etrs89_to_osgb36_threaded_vec(&mut eastings_vec, &mut northings_vec);
    (Array::from(eastings_vec), Array::from(northings_vec))
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
pub extern "C" fn convert_etrs89_to_ll_threaded(eastings: Array,
                                                northings: Array)
                                                -> (Array, Array) {
    let mut eastings_vec: &mut [f64] = eastings.into();
    let mut northings_vec: &mut [f64] = northings.into();
    convert_etrs89_to_ll_threaded_vec(&mut eastings_vec, &mut northings_vec);
    (Array::from(eastings_vec), Array::from(northings_vec))
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
pub extern "C" fn convert_osgb36_to_ll_threaded(eastings: Array,
                                                northings: Array)
                                                -> (Array, Array) {
    let mut eastings_vec: &mut [f64] = eastings.into();
    let mut northings_vec: &mut [f64] = northings.into();
    convert_osgb36_to_ll_threaded_vec(&mut eastings_vec, &mut northings_vec);
    (Array::from(eastings_vec), Array::from(northings_vec))
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
pub extern "C" fn convert_osgb36_to_etrs89_threaded(eastings: Array,
                                                    northings: Array)
                                                    -> (Array, Array) {
    let mut eastings_vec: &mut [f64] = eastings.into();
    let mut northings_vec: &mut [f64] = northings.into();
    convert_osgb36_to_etrs89_threaded_vec(&mut eastings_vec, &mut northings_vec);
    (Array::from(eastings_vec), Array::from(northings_vec))
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
pub extern "C" fn convert_epsg3857_to_wgs84_threaded(x: Array, y: Array) -> (Array, Array) {
    let mut x_vec: &mut [f64] = x.into();
    let mut y_vec: &mut [f64] = y.into();
    convert_epsg3857_to_wgs84_threaded_vec(&mut x_vec, &mut y_vec);
    (Array::from(x_vec), Array::from(y_vec))
}
