use std::mem;
use std::slice;

extern crate libc;
use libc::{c_void, c_double};

#[repr(C)]
pub struct Array {
    pub data: *const c_void,
    pub len: libc::size_t,
}

use super::convert_to_bng_threaded_vec;
use super::convert_to_bng_threaded_arr;
use super::convert_to_lonlat_threaded_vec;
use super::convert_to_osgb36_threaded_vec;
use super::convert_to_etrs89_threaded_vec;
use super::convert_etrs89_to_osgb36_threaded_vec;
use super::convert_etrs89_to_ll_threaded_vec;
use super::convert_osgb36_to_ll_threaded_vec;
use super::convert_osgb36_to_etrs89_threaded_vec;

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
    unsafe { Vec::from_raw_parts(lons.data as *mut c_double, lons.len, lons.len) };
    unsafe { Vec::from_raw_parts(lats.data as *mut c_double, lats.len, lats.len) };
}

impl Array {
    pub unsafe fn as_f64_slice(&self) -> &[f64] {
        assert!(!self.data.is_null());
        slice::from_raw_parts(self.data as *const f64, self.len as usize)
    }

    pub unsafe fn as_f64_slice_mut(&self) -> &mut[f64] {
        assert!(!self.data.is_null());
        slice::from_raw_parts_mut(self.data as *mut f64, self.len as usize)
    }

    pub fn from_vec<T>(mut vec: Vec<T>) -> Array {
        // Important to make length and capacity match
        // A better solution is to track both length and capacity
        vec.shrink_to_fit();

        let array = Array {
            data: vec.as_ptr() as *const libc::c_void,
            len: vec.len() as libc::size_t,
        };

        // Leak the memory, and now the raw pointer (and
        // eventually the FFI parent process) is the owner
        mem::forget(vec);

        array
    }

    pub fn from_slice(s: &mut[f64]) -> Array {
        let array = Array {
            data: s.as_ptr() as *const libc::c_void,
            len: s.len() as libc::size_t,
        };
        mem::forget(s);

        array
    }
}

/// A threaded, FFI-compatible wrapper for `lonlat_bng::convert_bng`
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
    let longitudes_vec = unsafe { longitudes.as_f64_slice().to_vec() };
    let latitudes_vec = unsafe { latitudes.as_f64_slice().to_vec() };
    let (eastings, northings) = convert_to_bng_threaded_vec(&longitudes_vec, &latitudes_vec);
    (Array::from_vec(eastings), Array::from_vec(northings))
    
}
#[no_mangle]
pub extern "C" fn convert_to_bng_threaded_slice(longitudes: Array, latitudes: Array) -> (Array, Array) {
    let longitudes = unsafe { longitudes.as_f64_slice_mut() };
    let latitudes = unsafe { latitudes.as_f64_slice_mut() };
    convert_to_bng_threaded_arr(longitudes, latitudes);
    (Array::from_slice(longitudes), Array::from_slice(latitudes))
}

/// A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_lonlat`](fn.convert_lonlat.html)
///
/// # Examples
///
/// See [`lonlat_bng::convert_to_bng_threaded`](fn.convert_to_bng_threaded.html) , substituting i32 vectors
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn convert_to_lonlat_threaded(eastings: Array, northings: Array) -> (Array, Array) {
    let eastings_vec = unsafe { eastings.as_f64_slice().to_vec() };
    let northings_vec = unsafe { northings.as_f64_slice().to_vec() };
    let (eastings_shifted, northings_shifted) = convert_to_lonlat_threaded_vec(&eastings_vec,
                                                                               &northings_vec);
    (Array::from_vec(eastings_shifted),
     Array::from_vec(northings_shifted))
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
    let longitudes_vec = unsafe { longitudes.as_f64_slice().to_vec() };
    let latitudes_vec = unsafe { latitudes.as_f64_slice().to_vec() };
    let (eastings, northings) = convert_to_osgb36_threaded_vec(&longitudes_vec, &latitudes_vec);
    (Array::from_vec(eastings), Array::from_vec(northings))
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
    let longitudes_vec = unsafe { longitudes.as_f64_slice().to_vec() };
    let latitudes_vec = unsafe { latitudes.as_f64_slice().to_vec() };
    let (eastings, northings) = convert_to_etrs89_threaded_vec(&longitudes_vec, &latitudes_vec);
    (Array::from_vec(eastings), Array::from_vec(northings))
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
    let eastings_vec = unsafe { eastings.as_f64_slice().to_vec() };
    let northings_vec = unsafe { northings.as_f64_slice().to_vec() };
    let (eastings_shifted, northings_shifted) =
        convert_etrs89_to_osgb36_threaded_vec(&eastings_vec, &northings_vec);
    (Array::from_vec(eastings_shifted),
     Array::from_vec(northings_shifted))
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
    let eastings_vec = unsafe { eastings.as_f64_slice().to_vec() };
    let northings_vec = unsafe { northings.as_f64_slice().to_vec() };
    let (eastings_shifted, northings_shifted) = convert_etrs89_to_ll_threaded_vec(&eastings_vec,
                                                                                  &northings_vec);
    (Array::from_vec(eastings_shifted),
     Array::from_vec(northings_shifted))
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
    let eastings_vec = unsafe { eastings.as_f64_slice().to_vec() };
    let northings_vec = unsafe { northings.as_f64_slice().to_vec() };
    let (eastings_shifted, northings_shifted) = convert_osgb36_to_ll_threaded_vec(&eastings_vec,
                                                                                  &northings_vec);
    (Array::from_vec(eastings_shifted),
     Array::from_vec(northings_shifted))
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
    let eastings_vec = unsafe { eastings.as_f64_slice().to_vec() };
    let northings_vec = unsafe { northings.as_f64_slice().to_vec() };
    let (eastings_shifted, northings_shifted) =
        convert_osgb36_to_etrs89_threaded_vec(&eastings_vec, &northings_vec);
    (Array::from_vec(eastings_shifted),
     Array::from_vec(northings_shifted))
}
