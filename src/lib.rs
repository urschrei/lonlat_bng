#![doc(html_root_url = "https://urschrei.github.io/lonlat_bng/")]
//! The `lonlat_bng` crate provides functions that convert decimal longitude
//! and latitude coordinates into [British National Grid](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid) coordinates, and vice versa.
//! Please note that this library does not make use of the [OSTN02](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) transformations, and should not be used for surveying work
//! or other applications requiring accuracy greater than ~5m.
//!
//! **Note that `lon`, `lat` coordinates outside the UK bounding box will be transformed to `(9999, 9999)`, which cannot be mapped.**
//!
//! # Examples
//!
//! ```
//! // Convert single values
//! assert_eq!((516276, 173141), lonlat_bng::convert_bng(&-0.32824866, &51.44533267).unwrap());
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
//! The crate also provides C-compatible wrapper functions, which are intended for use with FFI.
//! An example implementation using Python can be found at [Convertbng](https://github.com/urschrei/convertbng).
use std::f64;
use std::mem;
use std::slice;
use std::fmt;
use std::collections::HashMap;

extern crate libc;
use libc::{c_void, c_float, c_int};

extern crate rand;

extern crate crossbeam;
use crossbeam::scope;

extern crate num_cpus;

extern crate nalgebra;
use nalgebra::{Vec3, Mat3, DMat};

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

const s: f64 = 20.4894 * 0.000001;
// etc
const PI: f64 = f64::consts::PI;

const MAX_EASTING: i32  =  700000;
const MAX_NORTHING: i32 = 1250000;

const MIN_X_SHIFT: f64 =  86.275;
const MIN_Y_SHIFT: f64 = -81.603;
const MIN_Z_SHIFT: f64 =  43.982;

#[repr(C)]
pub struct Array {
    pub data: *const c_void,
    pub len: libc::size_t,
}


fn helmert(lon_vec: Vec<&f32>, lat_vec: Vec<&f32>) -> (Vec<f32>, Vec<f32>) {
    let t_array = Vec3::new(TX, TY, TZ);
    let params = Mat3::new(1. + s, RZS, -RYS, -RZS, 1. + s, RXS, RYS, -RXS, 1. + s);
    // let zeros: Vec<&f32> = vec![&1.; lon_vec.len()];
    let test = vec![1, 2, 3, 4, 5, 6];
    // let mut combined: Vec<&f32> = lon_vec.extend(lat_vec.iter()).collect();
    // combined.extend(zeros.iter().cloned());
    // let tmat = DMat::from_row_vec(3, zeros.len(), &combined);
    let tmat = DMat::from_row_vec(3, 2, &test);
    // build input vector of x, y, and z columns
    // let inp = DMat::from_row_vec(3, 1, vec![lon_vec, lat_vec, h_vec]);
    (vec![1.], vec![2.])
}

/// Free memory which Rust has allocated across the FFI boundary (i32 values)
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
/// let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
/// drop_int_array(eastings, northings);
/// ```
/// # Safety
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn drop_int_array(eastings: Array, northings: Array) {
    if eastings.data.is_null() {
        return;
    }
    if northings.data.is_null() {
        return;
    }
    unsafe { Vec::from_raw_parts(eastings.data as *mut c_int, eastings.len, eastings.len) };
    unsafe { Vec::from_raw_parts(northings.data as *mut c_int, northings.len, northings.len) };
}

/// Free memory which Rust has allocated across the FFI boundary (f32 values)
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
    unsafe { Vec::from_raw_parts(lons.data as *mut c_float, lons.len, lons.len) };
    unsafe { Vec::from_raw_parts(lats.data as *mut c_float, lats.len, lats.len) };
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

/// Bounds checking for input values
fn check<T>(to_check: T, bounds: (T, T)) -> Result<T, T>
    where T: std::cmp::PartialOrd + fmt::Display + Copy
{
    match to_check {
        to_check if bounds.0 <= to_check && to_check <= bounds.1 => Ok(to_check),
        _ => Err(to_check),
    }
}

/// Perform Longitude, Latitude to British National Grid conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_bng;
/// assert_eq!((516276, 173141), convert_bng(&-0.32824866, &51.44533267).unwrap());
#[allow(non_snake_case)]
pub fn convert_bng(longitude: &f32, latitude: &f32) -> Result<(c_int, c_int), f32> {
    // input is restricted to the UK bounding box
    let max_lon = 1.768960;
    let min_lon = -6.379880;
    let max_lat = 55.811741;
    let min_lat = 49.871159;
    // Convert bounds-checked input to degrees, or return an Err
    let lon_1: f64 = try!(check(*longitude, (min_lon, max_lon)).map_err(|e| e)) as f64 * PI / 180.;
    let lat_1: f64 = try!(check(*latitude, (min_lat, max_lat)).map_err(|e| e)) as f64 * PI / 180.;
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

    // TODO solve this for all lat and lon using matrices in an intermediate step?
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
    Ok((E.round() as c_int, N.round() as c_int))
}


/// Perform British National Grid Eastings, Northings to Longitude, Latitude conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_lonlat;
/// assert_eq!((-0.328248, 51.44534), convert_lonlat(&516276, &173141));
#[allow(non_snake_case)]
pub fn convert_lonlat(easting: &i32, northing: &i32) -> (c_float, c_float) {
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
    let minus_s = -s; // The scale factor -1
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
    let x_2 = tx + (1. + minus_s) * x_1 + (-rz) * y_1 + (ry) * z_1;
    let y_2 = ty + (rz) * x_1 + (1. + minus_s) * y_1 + (-rx) * z_1;
    let z_2 = tz + (-ry) * x_1 + (rx) * y_1 + (1. + minus_s) * z_1;

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
    (lon as c_float, lat as c_float)
}


// Input values should be valid ETRS89 grid references
fn ostn02_shifts(x: &i32, y :&i32) -> (f64, f64, f64) {
    let e_index: i32 = *x / 1000;
    let n_index: i32 = *y / 1000;

    // any of these could be Err, so use try!
    let s0_ref: (f64, f64, f64) = get_ostn_ref(&(e_index + 0), &(n_index + 0));
    let s1_ref: (f64, f64, f64) = get_ostn_ref(&(e_index + 1), &(n_index + 0));
    let s2_ref: (f64, f64, f64) = get_ostn_ref(&(e_index + 0), &(n_index + 1));
    let s3_ref: (f64, f64, f64) = get_ostn_ref(&(e_index + 1), &(n_index + 1));
    // only continue if we get Results for the above

    let x0 = e_index * 1000;
    let y0 = n_index * 1000;
    // offset within square
    let dx = *x - x0;
    let dy = *y - y0;

    let t = dx / 1000;
    let u = dy / 1000;

    let f0 = (1. - t as f64) * (1. - u as f64);
    let f1 = t as f64 * (1. - u as f64);
    let f2 = (1. - t as f64) * u as f64;
    let f3 = t as f64 * u as f64;

    let se = f0 * s0_ref.0 + f1 * s1_ref.0 + f2 * s2_ref.0 + f3 * s3_ref.0;
    let sn = f0 * s0_ref.1 + f1 * s1_ref.1 + f2 * s2_ref.1 + f3 * s3_ref.1;
    let sg = f0 * s0_ref.2 + f1 * s1_ref.2 + f2 * s2_ref.2 + f3 * s3_ref.2;

    (se, sn, sg)

}

fn get_ostn_ref(x: &i32, y: &i32) -> (f64, f64, f64) {

    // TODO populate ostn02 with the full OSTN02 data
    let mut keys = vec!["0a0266", "09f266", "09f267", "0a0267"];
    let mut values:Vec<(_, _, _)> = vec![(14192, 1769, 643), (14197, 1737, 641), (14247, 1732, 627), (14249, 1759, 628)];
    let ostn02 = keys.drain(..).zip(values.drain(..)).collect::<HashMap<_, (_, _, _)>>();
    let key = format!("{:03x}{:03x}", y, x);
    // some or None, so try! this
    let result = ostn02.get(&*key).unwrap();
    // if we get a hit
    let data2 = (result.0 as f64  / 1000. + MIN_X_SHIFT, result.1 as f64 / 1000. + MIN_Y_SHIFT, result.2 as f64 / 1000. + MIN_Z_SHIFT);
    data2
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
/// let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
/// ```
/// For an FFI implementation, see the code at [Convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).
///
/// # Safety
///
/// This function is unsafe because it accesses a raw pointer which could contain arbitrary data 
#[no_mangle]
pub extern "C" fn convert_to_bng_threaded(longitudes: Array, latitudes: Array) -> (Array, Array) {
    let lons = unsafe { longitudes.as_f32_slice().to_vec() };
    let lats = unsafe { latitudes.as_f32_slice().to_vec() };
    let (eastings, northings): (Vec<i32>, Vec<i32>) = convert_to_bng_threaded_vec(&lons, &lats);
    (Array::from_vec(eastings), Array::from_vec(northings))
}

/// A threaded wrapper for `lonlat_bng::convert_bng`
pub fn convert_to_bng_threaded_vec(longitudes: &Vec<f32>,
                                   latitudes: &Vec<f32>)
                                   -> (Vec<i32>, Vec<i32>) {
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
                    match convert_bng(orig_elem.0, orig_elem.1) {
                        Ok(res) => *res_elem = res,
                        // we don't care about the return value as such
                        Err(_) => *res_elem = (9999, 9999),
                    };
                }
            });
        }
    });
    let (eastings, northings): (Vec<i32>, Vec<i32>) = result.into_iter().unzip();
    (eastings, northings)
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
pub extern "C" fn convert_to_lonlat_threaded(eastings: Array, northings: Array) -> (Array, Array) {
    let eastings_vec = unsafe { eastings.as_i32_slice().to_vec() };
    let northings_vec = unsafe { northings.as_i32_slice().to_vec() };
    let (lons, lats) = convert_to_lonlat_threaded_vec(&eastings_vec, &northings_vec);
    (Array::from_vec(lons), Array::from_vec(lats))
}

/// A threaded wrapper for `lonlat_bng::convert_lonlat`
pub fn convert_to_lonlat_threaded_vec(eastings: &Vec<i32>,
                                      northings: &Vec<i32>)
                                      -> (Vec<f32>, Vec<f32>) {
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
    let (lons, lats): (Vec<f32>, Vec<f32>) = result.into_iter().unzip();
    (lons, lats)
}

#[cfg(test)]
mod tests {
    use super::get_ostn_ref;
    use super::ostn02_shifts;
    use super::drop_int_array;
    use super::convert_bng;
    use super::convert_lonlat;
    use super::convert_to_bng_threaded;
    use super::convert_to_lonlat_threaded;
    use super::check;
    use super::Array;

    extern crate libc;

    #[test]
    // original coordinates are 614300, 159900
    fn test_ostn_hashmap_retrieval() {
        let eastings = 614;
        let northings = 160;
        let expected = (100.46700000000001, -79.83399999999999, 44.625);
        assert_eq!(expected, get_ostn_ref(&eastings, &northings));

    }

    #[test]
    fn test_ostn02_shift_incorporation() {
        let eastings = 614300;
        let northings = 159900;
        // The offsets result in the following hex values for the E and N above:
        // "0a0266", "09f266", "09f267", "0a0267" (00, 10, 01, 11)
        let expected = (100.47200000000001, -79.866, 44.623);
        assert_eq!(expected, ostn02_shifts(&eastings, &northings));
    }

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
        let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_i32_slice() };
        let retval2 = unsafe { northings.as_i32_slice() };
        // the value's incorrect, but let's worry about that later
        assert_eq!(398915, retval[0]);
        assert_eq!(521545, retval2[0]);
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
        let (eastings, _) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_i32_slice() };
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
        let (lons, _) = convert_to_lonlat_threaded(easting_arr, northing_arr);
        let retval = unsafe { lons.as_f32_slice() };
        // We shouldn't really be using error margins, but it should be OK because
        // neither number is zero, or very close to, and on opposite sides of zero
        // http://floating-point-gui.de/errors/comparison/
        assert!((retval[0] - -0.32824799370716407).abs() / -0.32824799370716407 < 0.0000000001);
        // assert!((retval[1] - 51.44534026616287).abs() / 51.44534026616287 < 0.0000000001);
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
        let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
        drop_int_array(eastings, northings);
    }

    #[test]
    fn test_bng_conversion() {
        // verified to be correct at http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        assert_eq!((516276, 173141),
                   convert_bng(&-0.32824866, &51.44533267).unwrap());
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
        assert_eq!((516276, 173141), convert_bng(&181., &51.44533267).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_bad_lat() {
        assert_eq!((516276, 173141),
                   convert_bng(&-0.32824866, &-90.01).unwrap());
    }

    #[test]
    fn test_bad_threaded_conversion() {
        // above maximum longitude
        let lon_vec: Vec<f32> = vec![-6.379881];
        let lat_vec: Vec<f32> = vec![55.811741];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t,
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t,
        };
        let (eastings, _) = convert_to_bng_threaded(lon_arr, lat_arr);
        let retval = unsafe { eastings.as_i32_slice() };
        assert_eq!(9999, retval[0]);
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
}
