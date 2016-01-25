//! The `lonlat_bng` crate provides functions that convert decimal longitude
//! and latitude coordinates into [British National Grid](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid) coordinates, and vice versa.
//!
//! # Examples
//!
//! ```
//! // Convert single values
//! assert_eq!((516276, 173141), lonlat_bng::convert_bng(&-0.32824866, &51.44533267));
//! ```
//! ```
//! assert_eq!((-0.328248, 51.44534), lonlat_bng::convert_lonlat(&516276, &173141));
//! ```
//! ```
//! // Convert vectors using threaded functions
//! lonlat_bng::convert_to_bng_threaded_vec(vec![&-0.32824866], vec![&51.44533267]);
//! ```
//! ```
//! lonlat_bng::convert_to_lonlat_threaded_vec(vec![&516276], vec![&173141]);
//! ```
//! The crate also provides threaded versions of each function, which accept vectors of values, and wrapper functions for these, which are intended for use with FFI.
//! An example implementation using Python can be found at [Convertbng](https://github.com/urschrei/convertbng).
use std::f64;
use std::mem;
use std::slice;

extern crate libc;
use libc::{c_void, c_float, uint32_t};

extern crate rand;

extern crate crossbeam;
use crossbeam::scope;

extern crate num_cpus;

// Constants used for coordinate conversions
//
// Ellipsoids
const AIRY_1830_SEMI_MAJOR: f64 = 6377563.396;
const AIRY_1830_SEMI_MINOR: f64 = 6356256.909;
const GRS80_SEMI_MAJOR: f64 = 6378137.000;
const GRS80_SEMI_MINOR: f64 = 6356752.3141;
// Northing & easting of true origin (m)
const TRUE_ORIGIN_NORTHING: f64 = -100000.;
const TRUE_ORIGIN_EASTING: f64 = 400000.;
// For Helmert Transform to OSGB36, translations along the x, y, z axes
// When transforming to WGS84, reverse the signs
const TX: f64 = -446.448;
const TY: f64 = 125.157;
const TZ: f64 = -542.060;
// Rotations along the x, y, z axes, in seconds
const RXS: f64 = -0.1502;
const RYS: f64 = -0.2470;
const RZS: f64 = -0.8421;
// etc
const PI: f64 = f64::consts::PI;
//

#[repr(C)]
pub struct IntTuple {
    a: uint32_t,
    b: uint32_t,
}

#[repr(C)]
pub struct FloatTuple {
    a: c_float,
    b: c_float,
}

#[repr(C)]
pub struct Array {
    pub data: *const c_void,
    pub len: libc::size_t,
}

/// Free memory "leaked" by rust by sending it back across the FFI boundary
/// and reconstituting it (i32 values).
///
/// # Examples
///
/// ```
/// # extern crate libc;
/// let lon_vec: Vec<f32> = vec![-2.0183041005533306];
/// let lat_vec: Vec<f32> = vec![54.589097162646141];
/// let lon_arr = Array {
///     data: lon_vec.as_ptr() as *const libc::c_void,
///     len: lon_vec.len() as libc::size_t,
/// };
/// let lat_arr = Array {
///     data: lat_vec.as_ptr() as *const libc::c_void,
///     len: lat_vec.len() as libc::size_t,
/// };
/// let converted = convert_to_bng_threaded(lon_arr, lat_arr);
/// drop_int_array(converted);
/// ```
/// # Safety
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn drop_int_array(arr: Array) {
    if arr.data.is_null() {
        return;
    }
    unsafe { Vec::from_raw_parts(arr.data as *mut i32, arr.len, arr.len) };
}

/// Free memory "leaked" by rust by sending it back across the FFI boundary
/// and reconstituting it (f32 values).
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn drop_float_array(arr: Array) {
    if arr.data.is_null() {
        return;
    }
    unsafe { Vec::from_raw_parts(arr.data as *mut f32, arr.len, arr.len) };
}

impl Array {
    unsafe fn as_f32_slice(&self) -> &[f32] {
        assert!(!self.data.is_null());
        slice::from_raw_parts(self.data as *const f32, self.len as usize)
    }
    unsafe fn as_i32_slice(&self) -> &[i32] {
        assert!(!self.data.is_null());
        slice::from_raw_parts(self.data as *const i32, self.len as usize)
    }

    fn from_vec<T>(mut vec: Vec<T>) -> Array {
        // Important to make length and capacity match
        // A better solution is to track both length and capacity
        vec.shrink_to_fit();

        let array = Array {
            data: vec.as_ptr() as *const libc::c_void,
            len: vec.len() as libc::size_t,
        };

        // Whee! Leak the memory, and now the raw pointer (and
        // eventually Python) is the owner.
        mem::forget(vec);

        array
    }
}

/// Calculate the meridional radius of curvature
#[allow(non_snake_case)]
fn curvature(a: f64, F0: f64, e2: f64, lat: f64) -> f64 {
    a * F0 * (1. - e2) * (1. - e2 * lat.sin().powi(2)).powf(-1.5)
}

/// Perform Longitude, Latitude to British National Grid conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_bng;
/// assert_eq!((516276, 173141), convert_bng(&-0.32824866, &51.44533267));
#[allow(non_snake_case)]
pub fn convert_bng(longitude: &f32, latitude: &f32) -> (i32, i32) {
    // input is restricted to the UK bounding box
    assert!(-6.379880 <= *longitude && *longitude <= 1.768960,
            "Out of bounds! Longitude must be between -6.379880 and 1.768960: {}",
            longitude);
    assert!(49.871159 <= *latitude && *latitude <= 55.811741,
            "Out of bounds! Latitude must be between 49.871159 and 55.811741: {}",
            latitude);
    // Convert input to degrees
    let lat_1: f64 = *latitude as f64 * PI / 180.;
    let lon_1: f64 = *longitude as f64 * PI / 180.;
    // The GRS80 semi-major and semi-minor axes used for WGS84 (m)
    let a_1 = GRS80_SEMI_MAJOR;
    let b_1 = GRS80_SEMI_MINOR;
    // The eccentricity (squared) of the GRS80 ellipsoid
    let e2_1 = 1. - (b_1.powi(2)) / (a_1.powi(2));
    // Transverse radius of curvature
    let nu_1 = a_1 / (1. - e2_1 * lat_1.sin().powi(2)).sqrt();
    // Third spherical coordinate is 0, in this case
    let H: f64 = 0.;
    let x_1 = (nu_1 + H) * lat_1.cos() * lon_1.cos();
    let y_1 = (nu_1 + H) * lat_1.cos() * lon_1.sin();
    let z_1 = ((1. - e2_1) * nu_1 + H) * lat_1.sin();

    // Perform Helmert transform (to go between Airy 1830 (_1) and GRS80 (_2))
    let cst: f64 = 20.4894;
    let s = cst * (10 as f64).powi(-6);
    // The translations along x, y, z axes respectively
    let tx = TX;
    let ty = TY;
    let tz = TZ;
    // The rotations along x, y, z respectively, in seconds
    let rxs = RXS;
    let rys = RYS;
    let rzs = RZS;
    // In radians
    let rx = rxs * PI / (180. * 3600.);
    let ry = rys * PI / (180. * 3600.);
    let rz = rzs * PI / (180. * 3600.);
    let x_2 = tx + (1. + s) * x_1 + -rz * y_1 + ry * z_1;
    let y_2 = ty + rz * x_1 + (1. + s) * y_1 + -rx * z_1;
    let z_2 = tz + -ry * x_1 + rx * y_1 + (1. + s) * z_1;

    // The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
    let a = AIRY_1830_SEMI_MAJOR;
    let b = AIRY_1830_SEMI_MINOR;
    // The eccentricity of the Airy 1830 ellipsoid
    let e2 = 1. - b.powi(2) / a.powi(2);
    let p = (x_2.powi(2) + y_2.powi(2)).sqrt();
    // Initial value
    let mut lat = z_2.atan2((p * (1. - e2)));
    let mut latold = 2. * PI;
    // this is cheating, but not sure how else to initialise nu
    let mut nu: f64 = 1.;
    // Latitude is obtained by iterative procedure
    while (lat - latold).abs() > (10 as f64).powi(-16) {
        mem::swap(&mut lat, &mut latold);
        nu = a / (1. - e2 * latold.sin().powi(2)).sqrt();
        lat = (z_2 + e2 * nu * latold.sin()).atan2(p);
    }
    let lon = y_2.atan2(x_2);
    // Scale factor on the central meridian
    let F0: f64 = 0.9996012717;
    // Latitude of true origin (radians)
    let lat0 = 49. * PI / 180.;
    // Longitude of true origin and central meridian (radians)
    let lon0 = -2. * PI / 180.;
    // Northing & easting of true origin (m)
    let N0 = TRUE_ORIGIN_NORTHING;
    let E0 = TRUE_ORIGIN_EASTING;
    let n = (a - b) / (a + b);
    // Meridional radius of curvature
    let rho = curvature(a, F0, e2, lat);
    let eta2 = nu * F0 / rho - 1.;

    let M1 = (1. + n + (5. / 4.) * n.powi(2) + (5. / 4.) * n.powi(3)) * (lat - lat0);
    let M2 = (3. * n + 3. * n.powi(2) + (21. / 8.) * n.powi(3)) *
             ((lat.sin() * lat0.cos()) - (lat.cos() * lat0.sin())).ln_1p().exp_m1() *
             (lat + lat0).cos();
    let M3 = ((15. / 8.) * n.powi(2) + (15. / 8.) * n.powi(3)) * (2. * (lat - lat0)).sin() *
             (2. * (lat + lat0)).cos();
    let M4 = (35. / 24.) * n.powi(3) * (3. * (lat - lat0)).sin() * (3. * (lat + lat0)).cos();
    let M = b * F0 * (M1 - M2 + M3 - M4);

    let I = M + N0;
    let II = nu * F0 * lat.sin() * lat.cos() / 2.;
    let III = nu * F0 * lat.sin() * lat.cos().powi(3) * (5. - lat.tan().powi(2) + 9. * eta2) / 24.;
    let IIIA = nu * F0 * lat.sin() * lat.cos().powi(5) *
               (61. - 58. * lat.tan().powi(2) + lat.tan().powi(4)) / 720.;
    let IV = nu * F0 * lat.cos();
    let V = nu * F0 * lat.cos().powi(3) * (nu / rho - lat.tan().powi(2)) / 6.;
    let VI = nu * F0 * lat.cos().powi(5) *
             (5. - 18. * lat.tan().powi(2) + lat.tan().powi(4) + 14. * eta2 -
              58. * eta2 * lat.tan().powi(2)) / 120.;
    let N = I + II * (lon - lon0).powi(2) + III * (lon - lon0).powi(4) +
            IIIA * (lon - lon0).powi(6);
    let E = E0 + IV * (lon - lon0) + V * (lon - lon0).powi(3) + VI * (lon - lon0).powi(5);
    (E.round() as i32, N.round() as i32)
}


/// Perform British National Grid Eastings, Northings to Longitude, Latitude conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_lonlat;
/// assert_eq!((-0.328248, 51.44534), convert_lonlat(&516276, &173141));
#[allow(non_snake_case)]
pub fn convert_lonlat(easting: &i32, northing: &i32) -> (f32, f32) {
    // The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
    let a = AIRY_1830_SEMI_MAJOR;
    let b = AIRY_1830_SEMI_MINOR;
    // Scale factor on the central meridian
    let F0: f64 = 0.9996012717;
    // Latitude of true origin (radians)
    let lat0 = 49. * PI / 180.;
    // Longitude of true origin and central meridian (radians)
    let lon0 = -2. * PI / 180.;
    // Northing & easting of true origin (m)
    let N0 = TRUE_ORIGIN_NORTHING;
    let E0 = TRUE_ORIGIN_EASTING;
    // Eccentricity squared
    let e2 = 1. - b.powi(2) / a.powi(2);
    let n = (a - b) / (a + b);

    let mut lat = lat0;
    let mut M: f64 = 0.0;
    while (*northing as f64 - N0 - M) >= 0.00001 {
        lat = (*northing as f64 - N0 - M) / (a * F0) + lat;
        let M1 = (1. + n + (5. / 4.) * n.powi(3) + (5. / 4.) * n.powi(3)) * (lat - lat0);
        let M2 = (3. * n + 3. * n.powi(2) + (21. / 8.) * n.powi(3)) *
                 ((lat.sin() * lat0.cos()) - (lat.cos() * lat0.sin())).ln_1p().exp_m1() *
                 (lat + lat0).cos();
        let M3 = ((15. / 8.) * n.powi(2) + (15. / 8.) * n.powi(3)) * (2. * (lat - lat0)).sin() *
                 (2. * (lat + lat0)).cos();
        let M4 = (35. / 24.) * n.powi(3) * (3. * (lat - lat0)).sin() * (3. * (lat + lat0)).cos();
        // Meridional arc!
        M = b * F0 * (M1 - M2 + M3 - M4);
    }
    // Transverse radius of curvature
    let nu = a * F0 / (1. - e2 * lat.sin().powi(2)).sqrt();
    // Meridional radius of curvature
    let rho = curvature(a, F0, e2, lat);
    let eta2 = nu / rho - 1.;

    let secLat = 1. / lat.cos();
    let VII = lat.tan() / (2. * rho * nu);
    let VIII = lat.tan() / (24. * rho * nu.powi(3)) *
               (5. + 3. * lat.tan().powi(2) + eta2 - 9. * lat.tan().powi(2) * eta2);
    let IX = lat.tan() / (720. * rho * nu.powi(5)) *
             (61. + 90. * lat.tan().powi(2) + 45. * lat.tan().powi(4));
    let X = secLat / nu;
    let XI = secLat / (6. * nu.powi(3)) * (nu / rho + 2. * lat.tan().powi(2));
    let XII = secLat / (120. * nu.powi(5)) *
              (5. + 28. * lat.tan().powi(2) + 24. * lat.tan().powi(4));
    let XIIA = secLat / (5040. * nu.powi(7)) *
               (61. + 662. * lat.tan().powi(2) + 1320. * lat.tan().powi(4) +
                720. * lat.tan().powi(6));
    let dE = *easting as f64 - E0;
    // These are on the wrong ellipsoid currently: Airy1830 (Denoted by _1)
    let lat_1 = lat - VII * dE.powi(2) + VIII * dE.powi(4) - IX * dE.powi(6);
    let lon_1 = lon0 + X * dE - XI * dE.powi(3) + XII * dE.powi(5) - XIIA * dE.powi(7);

    // We Want to convert to the GRS80 ellipsoid
    // First, convert to cartesian from spherical polar coordinates
    let H = 0.;
    let x_1 = (nu / F0 + H) * lat_1.cos() * lon_1.cos();
    let y_1 = (nu / F0 + H) * lat_1.cos() * lon_1.sin();
    let z_1 = ((1. - e2) * nu / F0 + H) * lat_1.sin();

    // Perform Helmert transform (to go between Airy 1830 (_1) and GRS80 (_2))
    let s = -20.4894 * (10. as f64).powi(-6); // The scale factor -1
    // The translations along x, y, z axes respectively
    let tx = TX.abs();
    let ty = TY * -1.;
    let tz = TZ.abs();
    // The rotations along x, y, z respectively, in seconds
    let rxs = RXS * -1.;
    let rys = RYS * -1.;
    let rzs = RZS * -1.;

    let rx = rxs * PI / (180. * 3600.);
    let ry = rys * PI / (180. * 3600.);
    let rz = rzs * PI / (180. * 3600.); // In radians
    let x_2 = tx + (1. + s) * x_1 + (-rz) * y_1 + (ry) * z_1;
    let y_2 = ty + (rz) * x_1 + (1. + s) * y_1 + (-rx) * z_1;
    let z_2 = tz + (-ry) * x_1 + (rx) * y_1 + (1. + s) * z_1;

    // Back to spherical polar coordinates from cartesian
    // Need some of the characteristics of the new ellipsoid
    // The GRS80 semi-major and semi-minor axes used for WGS84(m)
    let a_2 = GRS80_SEMI_MAJOR;
    let b_2 = GRS80_SEMI_MINOR;
    // The eccentricity of the GRS80 ellipsoid
    let e2_2 = 1. - b_2.powi(2) / a_2.powi(2);
    let p = (x_2.powi(2) + y_2.powi(2)).sqrt();

    // Lat is obtained by iterative procedure
    // Initial value
    let mut lat = z_2.atan2((p * (1. - e2_2)));
    let mut latold = 2. * PI;
    let mut nu_2: f64;
    while (lat - latold).abs() > (10. as f64).powi(-16) {
        mem::swap(&mut lat, &mut latold);
        nu_2 = a_2 / (1. - e2_2 * latold.sin().powi(2)).sqrt();
        lat = (z_2 + e2_2 * nu_2 * latold.sin()).atan2(p);
    }

    let mut lon = y_2.atan2(x_2);
    lat = lat * 180. / PI;
    lon = lon * 180. / PI;
    (lon as f32, lat as f32)
}

/// A C-compatible wrapper for `lonlat_bng::convert_bng`
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn convert_vec_c(longitudes: Array, latitudes: Array) -> Array {
    // we're receiving floats
    let lon = unsafe { longitudes.as_f32_slice() };
    let lat = unsafe { latitudes.as_f32_slice() };
    // copy values and combine
    let orig = lon.iter()
                  .zip(lat.iter());
    // carry out the conversion
    let result = orig.map(|elem| convert_bng(elem.0, elem.1));
    // convert back to vector of unsigned integer Tuples
    let nvec = result.map(|ints| {
                         IntTuple {
                             a: ints.0 as u32,
                             b: ints.1 as u32,
                         }
                     })
                     .collect();
    Array::from_vec(nvec)
}

/// A threaded, FFI-compatible wrapper for `lonlat_bng::convert_bng`
///
/// # Examples
///
/// ```
/// extern crate libc;
/// let lon_vec: Vec<f32> = vec![-2.0183041005533306,
///                              0.95511887434519682,
///                              0.44975855518383501,
///                              -0.096813621191803811,
///                              -0.36807065656416427,
///                              0.63486335458665621];
/// let lat_vec: Vec<f32> = vec![54.589097162646141,
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
/// let converted = convert_to_bng_threaded(lon_arr, lat_arr);
/// ```
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn convert_to_bng_threaded(longitudes: Array, latitudes: Array) -> Array {
    let lons = unsafe { longitudes.as_f32_slice().to_vec() };
    let lats = unsafe { latitudes.as_f32_slice().to_vec() };
    let result = convert_to_bng_threaded_vec(&lons, &lats);
    Array::from_vec(result)
}

/// A threaded wrapper for `lonlat_bng::convert_bng`
pub fn convert_to_bng_threaded_vec(longitudes: &Vec<f32>, latitudes: &Vec<f32>) -> Vec<(i32, i32)> {
    let numthreads = num_cpus::get() as usize;
    let orig: Vec<(&f32, &f32)> = longitudes.iter().zip(latitudes.iter()).collect();
    let mut result = vec![(1, 1); orig.len()];
    let mut size = orig.len() / numthreads;
    if orig.len() % numthreads > 0 {
        size += 1;
    }
    size = std::cmp::max(1, size);
    crossbeam::scope(|scope| {
        for (res_chunk, orig_chunk) in result.chunks_mut(size).zip(orig.chunks(size)) {
            scope.spawn(move || {
                for (res_elem, orig_elem) in res_chunk.iter_mut().zip(orig_chunk.iter()) {
                    *res_elem = convert_bng(orig_elem.0, orig_elem.1);
                }
            });
        }
    });
    result
}

/// A threaded, FFI-compatible wrapper for `lonlat_bng::convert_lonlat`
///
/// # Examples
///
/// See `lonlat_bng::convert_to_bng_threaded` for examples, substituting i32 vectors
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn convert_to_lonlat_threaded(eastings: Array, northings: Array) -> Array {
    let eastings_vec = unsafe { eastings.as_i32_slice().to_vec() };
    let northings_vec = unsafe { northings.as_i32_slice().to_vec() };
    let result = convert_to_lonlat_threaded_vec(&eastings_vec, &northings_vec);
    Array::from_vec(result)
}

/// A threaded wrapper for `lonlat_bng::convert_lonlat`
pub fn convert_to_lonlat_threaded_vec(eastings: &Vec<i32>, northings: &Vec<i32>) -> Vec<(f32, f32)> {
    let numthreads = num_cpus::get() as usize;
    let orig: Vec<(&i32, &i32)> = eastings.iter().zip(northings.iter()).collect();
    let mut result: Vec<(f32, f32)> = vec![(1.0, 1.0); orig.len()];
    let mut size = orig.len() / numthreads;
    if orig.len() % numthreads > 0 {
        size += 1;
    }
    size = std::cmp::max(1, size);
    crossbeam::scope(|scope| {
        for (res_chunk, orig_chunk) in result.chunks_mut(size).zip(orig.chunks(size)) {
            scope.spawn(move || {
                for (res_elem, orig_elem) in res_chunk.iter_mut().zip(orig_chunk.iter()) {
                    *res_elem = convert_lonlat(orig_elem.0, orig_elem.1);
                }
            });
        }
    });
    result
}

#[cfg(test)]
mod tests {
    use super::drop_int_array;
    use super::convert_bng;
    use super::convert_lonlat;
    use super::convert_vec_c;
    use super::convert_to_bng_threaded;
    use super::convert_to_lonlat_threaded;
    use super::Array;

    extern crate libc;

    #[test]
    fn test_threaded_bng_conversion() {
        let lon_vec: Vec<f32> = vec![-2.0183041005533306,
                                     0.95511887434519682,
                                     0.44975855518383501,
                                     -0.096813621191803811,
                                     -0.36807065656416427,
                                     0.63486335458665621];
        let lat_vec: Vec<f32> = vec![54.589097162646141,
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
        let converted = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { converted.as_i32_slice() };
        // the value's incorrect, but let's worry about that later
        assert_eq!(398915, retval[0]);
        assert_eq!(521545, retval[1]);
    }

    #[test]
    fn test_threaded_bng_conversion_single() {
        // I spent 8 hours confused cos I didn't catch that chunks(0) is invalid
        let lon_vec: Vec<f32> = vec![-2.0183041005533306];
        let lat_vec: Vec<f32> = vec![54.589097162646141];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let converted = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { converted.as_i32_slice() };
        assert_eq!(398915, retval[0]);
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
        let converted = convert_to_lonlat_threaded(easting_arr, northing_arr);
        let retval = unsafe { converted.as_f32_slice() };
        // We shouldn't really be using error margins, but it should be OK because
        // neither number is zero, or very close to, and on opposite sides of zero
        // http://floating-point-gui.de/errors/comparison/
        assert!((retval[0] - -0.32824799370716407).abs() / -0.32824799370716407 < 0.0000000001);
        // assert!((retval[1] - 51.44534026616287).abs() / 51.44534026616287 < 0.0000000001);
    }

    #[test]
    fn test_nonthreaded_vector_bng_conversion() {
        let lon_vec: Vec<f32> = vec![-2.0183041005533306,
                                     0.95511887434519682,
                                     0.44975855518383501,
                                     -0.096813621191803811,
                                     -0.36807065656416427,
                                     0.63486335458665621];
        let lat_vec: Vec<f32> = vec![54.589097162646141,
                                     51.560873800587828,
                                     50.431429161121699,
                                     54.535021436247419,
                                     50.839059313135706,
                                     55.412189281234419];

        // from http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        let correct_values = vec![398915, 521545, 604932, 188804, 574082, 61931, 523242, 517193,
                                  515004, 105661, 566898, 616298];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let converted = convert_vec_c(lon_arr, lat_arr);
        let retval = unsafe { converted.as_i32_slice() };
        let combined: Vec<(&i32, &i32)> = retval.iter()
                                                .zip(correct_values.iter())
                                                .collect();
        for val in combined.iter() {
            assert_eq!(val.0, val.1);
        }
    }

    #[test]
    fn test_drop_int_array() {
        let lon_vec: Vec<f32> = vec![-2.0183041005533306];
        let lat_vec: Vec<f32> = vec![54.589097162646141];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let converted = convert_to_bng_threaded(lon_arr, lat_arr);
        drop_int_array(converted);
    }

    #[test]
    fn test_bng_conversion() {
        // verified to be correct at http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        assert_eq!((516276, 173141), convert_bng(&-0.32824866, &51.44533267));
    }

    #[test]
    fn test_lonlat_conversion() {
        let res = convert_lonlat(&516276, &173141);
        // We shouldn't really be using error margins, but it should be OK because
        // neither number is zero, or very close to, and on opposite sides of zero
        // epsilon is .000001 here, because BNG coords are 6 digits, so
        // we should be fine if the error is in the 7th digit (i.e. < epsilon)
        // http://floating-point-gui.de/errors/comparison/
        assert!(((res.0 - -0.328248269313) / -0.328248269313).abs() < 0.000001);
        assert!(((res.1 - 51.4453318435) / 51.4453318435).abs() < 0.000001);
    }

    #[test]
    #[should_panic]
    fn test_bad_lon() {
        assert_eq!((516276, 173141), convert_bng(&181., &51.44533267));
    }

    #[test]
    #[should_panic]
    fn test_bad_lat() {
        assert_eq!((516276, 173141), convert_bng(&-0.32824866, &-90.01));
    }
}
