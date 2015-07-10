//! The `lonlat_bng` crate provides a function that converts decimal longitude
//! and latitude coordinates into British National Grid Coordinates
//!
//! Examples
//!
//!```
//! assert_eq!((516276, 173141), lonlat_bng::convert(-0.32824866, 51.44533267));
//!```

use std::{f32, f64};
// use std::f32::consts;
use std::mem;
use std::slice;
use std::thread::{self, JoinHandle};

extern crate libc;
use libc::{size_t, c_void, c_float, uint32_t};

extern crate rand;

const NUMTHREADS: usize = 7;

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
    data: *const c_void,
    len: libc::size_t,
}

#[no_mangle]
pub extern fn drop_array(p: *const Array) {
    if p.is_null() { return }
    unsafe { drop(p) };
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
            len: vec.len() as libc::size_t
        };

        // Whee! Leak the memory, and now the raw pointer (and
        // eventually C) is the owner.
        mem::forget(vec);

        array
    }
}

// http://stackoverflow.com/a/28124775/155423
fn round(x: f64) -> f64 {
    let y = x.floor();
    if x == y {
        x
    } else {
        let z = (2.0*x-y).floor();
        z * x.signum() // Should use copysign, but not stably-available
    }
}

/// This function performs lon, lat to BNG conversion
///
/// Examples
///
/// ```
/// use lonlat_bng::convert_bng;
/// assert_eq!((516276, 173141), convert_bng(-0.32824866, 51.44533267));
#[allow(non_snake_case)]
#[no_mangle]
pub extern fn convert_bng(input_lon: f64, input_lat: f64) -> (i32, i32) {
    // match is restricted to the UK bounding box
    match input_lon {
        -6.379880...1.768960 => input_lon,
        _ => panic!("Out of bounds! Longitude must be between -180.00 and 180.00: {:?}", input_lon)
    };
    match input_lat {
        49.871159...55.811741 => input_lat,
        _ => panic!("Out of bounds! Latitude must be between -90.00 and 90.00: {:?}", input_lat)
    };
    let pi: f64 = f64::consts::PI;
    //Convert input to degrees
    let lat_1: f64 = input_lat * pi / 180.;
    let lon_1: f64 = input_lon * pi / 180.;
    // The GSR80 semi-major and semi-minor axes used for WGS84 (m)
    let a_1: f64 = 6378137.;
    let b_1: f64 = 6356752.3141;
    // The eccentricity (squared) of the GRS80 ellipsoid
    let e2_1: f64 = 1. - (b_1.powi(2)) / (a_1.powi(2));
    // Transverse radius of curvature
    let nu_1: f64 = a_1 / (1. - e2_1 * lat_1.sin().powi(2)).sqrt();
    // Third spherical coordinate is 0, in this case
    let H: f64 = 0.;
    let x_1: f64 = (nu_1 + H) * lat_1.cos() * lon_1.cos();
    let y_1: f64 = (nu_1 + H) * lat_1.cos() * lon_1.sin();
    let z_1: f64 = ((1. - e2_1) * nu_1 + H) * lat_1.sin();

    // Perform Helmert transform (to go between Airy 1830 (_1) and GRS80 (_2))
    let small: f64 = 10.;
    let cst: f64 = 20.4894;
    let s: f64 =  cst * small.powi(-6);
    // The translations along x, y, z axes respectively
    let tx: f64 = -446.448;
    let ty: f64 = 125.157;
    let tz: f64 = -542.060;
    // The rotations along x, y, z respectively, in seconds
    let rxs: f64 = -0.1502;
    let rys: f64 = -0.2470;
    let rzs: f64 = -0.8421;
    // In radians
    let rx: f64 = rxs * pi / (180. * 3600.);
    let ry: f64 = rys * pi / (180. * 3600.);
    let rz: f64 = rzs * pi / (180. * 3600.);
    let x_2: f64 = tx + (1. + s) * x_1 + -rz * y_1 + ry * z_1;
    let y_2: f64 = ty + rz * x_1 + (1. + s) * y_1 + -rx * z_1;
    let z_2: f64 = tz + -ry * x_1 + rx * y_1 + (1. + s) * z_1;

    // The GSR80 semi-major and semi-minor axes used for WGS84 (m)
    let a: f64 = 6377563.396;
    let b: f64 = 6356256.909;
    // The eccentricity of the Airy 1830 ellipsoid
    let e2: f64 = 1. - b.powi(2) / a.powi(2);
    let p: f64 = (x_2.powi(2) + y_2.powi(2)).sqrt();
    // Initial value
    let mut lat: f64 = z_2.atan2((p * (1. - e2)));
    let mut latold: f64 = 2. * pi;
    // this is cheating, but not sure how else to initialise nu
    let mut nu: f64 = 1.;
    // Latitude is obtained by iterative procedure
    while (lat - latold).abs() > small.powi(-16) {
        mem::swap(&mut lat, &mut latold);
        nu = a / (1. - e2 * latold.sin().powi(2)).sqrt();
        lat = (z_2 + e2 * nu * latold.sin()).atan2(p);
    };
    let lon: f64 = y_2.atan2(x_2);
    // Scale factor on the central meridian
    let F0: f64 = 0.9996012717;
    // Latitude of true origin (radians)
    let lat0: f64 = 49. * pi / 180.;
    // Longtitude of true origin and central meridian (radians)
    let lon0: f64 = -2. * pi / 180.;
    // Northing & easting of true origin (m)
    let N0: f64 = -100000.;
    let E0: f64 = 400000.;
    let n: f64 = (a - b) / (a + b);
    // Meridional radius of curvature
    let rho: f64 = a * F0 * (1. - e2) * (1. - e2 * lat.sin().powi(2)).powf(-1.5);
    let eta2: f64 = nu * F0 / rho - 1.;

    let M1: f64 = (1. + n + (5. / 4.) * n.powi(2)
        + (5. / 4.) * n.powi(3)) * (lat - lat0);
    let M2: f64 = (3. * n + 3. * n.powi(2) + (21. / 8.)
        * n.powi(3)) * (lat - lat0).sin() * (lat + lat0).cos();
    let M3: f64 = ((15. / 8.) * n.powi(2) + (15. / 8.)
        * n.powi(3)) * (2. * (lat-lat0)).sin() * (2. * (lat + lat0)).cos();
    let M4: f64 = (35. / 24.) * n.powi(3) * (3. * (lat - lat0)).sin()
        * (3. * (lat + lat0)).cos();
    let M: f64 = b * F0 * (M1 - M2 + M3 - M4);

    let I: f64 = M + N0;
    let II: f64 = nu * F0 * lat.sin() * lat.cos() / 2.;
    let III: f64 = nu * F0 * lat.sin() * lat.cos().powi(3)
        * (5. - lat.tan().powi(2) + 9. * eta2) / 24.;
    let IIIA: f64 = nu * F0 * lat.sin() * lat.cos().powi(5)
        * (61. - 58. * lat.tan().powi(2)
        + lat.tan().powi(4)) / 720.;
    let IV: f64 = nu * F0 * lat.cos();
    let V: f64 = nu * F0 * lat.cos().powi(3)
        * (nu / rho - lat.tan().powi(2)) / 6.;
    let VI: f64 = nu * F0 * lat.cos().powi(5)
        * (5. - 18. * lat.tan().powi(2) + lat.tan().powi(4)
        + 14. * eta2 - 58. * eta2 * lat.tan().powi(2)) / 120.;
    let N: f64 = I + II * (lon - lon0).powi(2)
        + III * (lon - lon0).powi(4) + IIIA * (lon - lon0).powi(6);
    let E: f64 = E0 + IV * (lon - lon0) + V * (lon - lon0).powi(3)
        + VI * (lon - lon0).powi(5);
    (round(E) as i32, round(N) as i32)
}


/// This function performs BNG Eastings, Northings to lon, lat conversion
///
/// Examples
///
/// ```
/// use lonlat_bng::convert_lonlat;
/// assert_eq!((-0.328248, 51.44534), convert_lonlat(516276, 173141)));
#[allow(non_snake_case)]
#[no_mangle]
pub extern fn convert_lonlat(input_e: i32, input_n: i32) -> (f64, f64) {
    let pi: f64 = f64::consts::PI;
    // The Airy 180 semi-major and semi-minor axes used for OSGB36 (m)
    let a: f64 = 6377563.396;
    let b: f64 = 6356256.909;
    // Scale factor on the central meridian
    let F0: f64 = 0.9996012717;
    // Latitude of true origin (radians)
    let lat0: f64 = 49. * pi / 180.;
    // Longtitude of true origin and central meridian (radians)
    let lon0: f64 = -2. * pi / 180.;
    // Northing & easting of true origin (m)
    let N0 = -100000.;
    let E0 = 400000.;
    // Eccentricity squared
    let e2 = 1. - (b * b) / (a * a);
    let n = (a - b) / (a + b);

    let mut lat = lat0;
    let mut M = 0.0;
    while (input_n as f64 - N0 - M) >= 0.00001 {
        lat = (input_n as f64 - N0 - M) / (a * F0) + lat;
        let M1 = (1. + n + (5. / 4.) * n.powi(3)
            + (5. / 4.) * n.powi(3)) * (lat - lat0);
        let M2 = (3. * n + 3. * n.powi(2) + (21. / 8.)
            * n.powi(3)) * (lat - lat0).sin() * (lat + lat0).cos();
        let M3 = ((15. / 8.) * n.powi(2) + (15. / 8.) * n.powi(3))
            * (2. * (lat - lat0)).sin() * (2. * (lat + lat0)).cos();
        let M4 = (35. / 24.) * n.powi(3) * (3. * (lat - lat0)).sin()
            * (3. * (lat + lat0)).cos();
        // Meridional arc!
        M = b * F0 * (M1 - M2 + M3 - M4);
    }
    // Transverse radius of curvature
    let nu = a * F0 / (1. - e2 * lat.sin().powi(2)).sqrt();
    // Meridional radius of curvature
    let rho = a * F0 * (1. - e2) * (1. - e2 * lat.sin().powi(2)).powf(-1.5);
    let eta2 = nu / rho - 1.;

    let secLat = 1. / lat.cos();
    let VII = lat.tan() / (2. * rho * nu);
    let VIII = lat.tan() / (24. * rho * nu.powi(3))
        * (5. + 3. * lat.tan().powi(2)
        + eta2 - 9. * lat.tan().powi(2) * eta2);
    let IX = lat.tan() / (720. * rho * nu.powi(5))
        * (61. + 90. * lat.tan().powi(2)
        + 45. * lat.tan().powi(4));
    let X = secLat / nu;
    let XI = secLat / (6. * nu.powi(3)) * (nu / rho + 2. * lat.tan().powi(2));
    let XII = secLat / (120. * nu.powi(5))
        * (5. + 28. * lat.tan().powi(2) + 24.
        * lat.tan().powi(4));
    let XIIA = secLat / (5040. * nu.powi(7))
        * (61. + 662. * lat.tan().powi(2) + 1320.
        * lat.tan().powi(4) + 720. * lat.tan().powi(6));
    let dE = input_e as f64 - E0;
    //  These are on the wrong ellipsoid currently: Airy1830 (Denoted by _1)
    let lat_1 = lat - VII * dE.powi(2)
        + VIII * dE.powi(4) - IX * dE.powi(6);
    let lon_1 = lon0 + X * dE - XI * dE.powi(3)
        + XII * dE.powi(5) - XIIA * dE.powi(7);

    // We Want to convert to the GRS80 ellipsoid
    // First, convert to cartesian from spherical polar coordinates
    let H = 0.;
    let x_1 = (nu / F0 + H) * lat_1.cos() * lon_1.cos();
    let y_1 = (nu / F0 + H) * lat_1.cos() * lon_1.sin();
    let z_1 = ((1. - e2) * nu / F0 + H) * lat_1.sin();

    // Perform Helmert transform (to go between Airy 1830 (_1) and GRS80 (_2))
    let ten: f64 = 10.;
    let s = -20.4894 * ten.powi(-6); // The scale factor -1
    let tx = 446.448;
    let ty = -125.157;
    let tz = 542.060; // The translations along x,y,z axes respectively
    let rxs =0.1502;
    let rys= 0.2470;
    let rzs = 0.8421; // The rotations along x,y,z respectively, in seconds
    let rx = rxs * pi / (180.*3600.);
    let ry = rys * pi / (180.*3600.);
    let rz = rzs * pi / (180.*3600.); // In radians
    let x_2 = tx + (1. + s) * x_1 + (-rz) * y_1 + (ry) * z_1;
    let y_2 = ty + (rz) * x_1 + (1. + s) * y_1 + (-rx) * z_1;
    let z_2 = tz + (-ry) * x_1 + (rx) * y_1 + (1. + s) * z_1;

    // Back to spherical polar coordinates from cartesian
    // Need some of the characteristics of the new ellipsoid
    let a_2 = 6378137.00;
    let b_2 = 6356752.3141; // The GSR80 semi-major and semi-minor axes used for WGS84(m)
    let e2_2 = 1. - (b_2 * b_2) / (a_2 * a_2); // The eccentricity of the GRS80 ellipsoid
    let p = (x_2.powi(2) + y_2.powi(2)).sqrt();

    // Lat is obtained by iterative procedure
    let mut lat = z_2.atan2((p * (1. - e2_2))); // Initial value
    let mut latold = 2. *pi;
    let mut nu_2: f64;
    while (lat - latold).abs() > ten.powi(-16) {
        mem::swap(&mut lat, &mut latold);
        nu_2 = a_2 / (1. - e2_2 * latold.sin().powi(2)).sqrt();
        lat = (z_2 + e2_2 * nu_2 * latold.sin()).atan2(p);
    }

    let mut lon = y_2.atan2(x_2);
    lat = lat * 180. / pi;
    lon = lon * 180. / pi;
    return (lon, lat)
}

/// A safer C-compatible wrapper for convert_bng()
#[no_mangle]
pub extern fn convert_vec_c(lon: Array, lat: Array) -> Array {
    // we're receiving floats
    let lon = unsafe { lon.as_f32_slice() };
    let lat = unsafe { lat.as_f32_slice() };
    // copy values and combine
    let orig = lon
        .iter()
        .cloned()
    .zip(lat
        .iter()
        .cloned());
    // carry out the conversion 
    let result = orig
        .map(|elem| convert_bng(elem.0 as f64, elem.1 as f64));
    // convert back to vector of unsigned integer Tuples
    let nvec = result
        .map(|ints| IntTuple { a: ints.0 as u32, b: ints.1 as u32 })
        .collect();
    Array::from_vec(nvec)
}

/// A threaded version of the C-compatible wrapper for convert_bng()
#[no_mangle]
pub extern fn convert_to_bng(lon: Array, lat: Array) -> Array {
    // we're receiving floats
    let lon = unsafe { lon.as_f32_slice() };
    let lat = unsafe { lat.as_f32_slice() };
    // copy values and combine
    let orig: Vec<(f32, f32)> = lon
        .iter()
        .cloned()
    .zip(lat
        .iter()
        .cloned())
    .collect();

    let mut guards: Vec<JoinHandle<Vec<(i32, i32)>>> = vec!();
    // split into slices
    let mut size = orig.len() / NUMTHREADS;
    if orig.len() % NUMTHREADS > 0 { size += 1; }
    // if orig.len() == 0, we need another adjustment
    size = std::cmp::max(1, size);
    for chunk in orig.chunks(size) {
        let chunk = chunk.to_owned();
        let g = thread::spawn(move || chunk
            .into_iter()
            .map(|elem| convert_bng(elem.0 as f64, elem.1 as f64))
            .collect());
        guards.push(g);
    }
    let mut result: Vec<IntTuple> = Vec::with_capacity(orig.len());
    for g in guards {
        result.extend(g.join().unwrap().into_iter()
                       .map(|ints| IntTuple { a: ints.0 as u32, b: ints.1 as u32 }));
    }
    Array::from_vec(result)
}

/// A threaded version of the C-compatible wrapper for convert_lonlat()
#[no_mangle]
pub extern fn convert_to_lonlat(eastings: Array, northings: Array) -> Array {
    // we're receiving floats
    let lon = unsafe { eastings.as_i32_slice() };
    let lat = unsafe { northings.as_i32_slice() };
    // copy values and combine
    let orig: Vec<(i32, i32)> = lon
        .iter()
        .cloned()
    .zip(lat
        .iter()
        .cloned())
    .collect();

    let mut guards: Vec<JoinHandle<Vec<(f64, f64)>>> = vec!();
    // split into slices
    let mut size = orig.len() / NUMTHREADS;
    if orig.len() % NUMTHREADS > 0 { size += 1; }
    // if orig.len() == 0, we need another adjustment
    size = std::cmp::max(1, size);
    for chunk in orig.chunks(size) {
        let chunk = chunk.to_owned();
        let g = thread::spawn(move || chunk
            .into_iter()
            .map(|elem| convert_lonlat(elem.0, elem.1))
            .collect());
        guards.push(g);
    }
    let mut result: Vec<FloatTuple> = Vec::with_capacity(orig.len());
    for g in guards {
        result.extend(g.join().unwrap().into_iter()
                       .map(|floats| FloatTuple { a: floats.0 as f32, b: floats.1 as f32 }));
    }
    Array::from_vec(result)
}

#[cfg(test)]
mod tests {
    use super::convert_bng;
    use super::convert_lonlat;
    use super::convert_vec_c;
    use super::convert_to_bng;
    use super::convert_to_lonlat;
    use super::Array;

    extern crate libc;
    use libc::{size_t};

    #[test]
    fn test_threaded_vector_conversion() {
        let lon_vec: Vec<f32> = vec!(
            -2.0183041005533306,
            0.95511887434519682,
            0.44975855518383501,
            -0.096813621191803811,
            -0.36807065656416427,
            0.63486335458665621);
        let lat_vec: Vec<f32> = vec!(
            54.589097162646141,
            51.560873800587828,
            50.431429161121699,
            54.535021436247419,
            50.839059313135706,
            55.412189281234419);
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t
        };
        let converted = convert_to_bng(lon_arr, lat_arr);
        let retval = unsafe{ converted.as_i32_slice() };
        // the value's incorrect, but let's worry about that later
        assert_eq!(398915, retval[0]);
        assert_eq!(521545, retval[1]);
    }

    #[test]
    fn test_threaded_bng_conversion_single() {
        // I spent 8 hours confused cos I didn't catch that chunks(0) is invalid
        let lon_vec: Vec<f32> = vec!(
            -2.0183041005533306);
        let lat_vec: Vec<f32> = vec!(
            54.589097162646141);
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t
        };
        let converted = convert_to_bng(lon_arr, lat_arr);
        let retval = unsafe{ converted.as_i32_slice() };
        assert_eq!(398915, retval[0]);
    }

    #[test]
    fn test_threaded_lonlat_conversion_single() {
        let easting_vec: Vec<i32> = vec!(
            516276);
        let northing_vec: Vec<i32> = vec!(
            173141);
        let easting_arr = Array {
            data: easting_vec.as_ptr() as *const libc::c_void,
            len: easting_vec.len() as libc::size_t
        };
        let northing_arr = Array {
            data: northing_vec.as_ptr() as *const libc::c_void,
            len: northing_vec.len() as libc::size_t
        };
        let converted = convert_to_lonlat(easting_arr, northing_arr);
        let retval = unsafe{ converted.as_f32_slice() };
        assert_eq!(-0.328248, retval[0]);
    }

    #[test]
    fn test_nonthreaded_vector_conversion() {
         let lon_vec: Vec<f32> = vec!(
            -2.0183041005533306,
            0.95511887434519682,
            0.44975855518383501,
            -0.096813621191803811,
            -0.36807065656416427,
            0.63486335458665621);
        let lat_vec: Vec<f32> = vec!(
            54.589097162646141,
            51.560873800587828,
            50.431429161121699,
            54.535021436247419,
            50.839059313135706,
            55.412189281234419);

        // from http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        let correct_values = vec![
            398915,
            521545,
            604932,
            188804,
            574082,
            61931,
            523242,
            517193,
            515004,
            105661,
            566898,
            616298
            ];
        let lon_arr = Array {
            data: lon_vec.as_ptr() as *const libc::c_void,
            len: lon_vec.len() as libc::size_t
        };
        let lat_arr = Array {
            data: lat_vec.as_ptr() as *const libc::c_void,
            len: lat_vec.len() as libc::size_t
        };
        let converted = convert_vec_c(lon_arr, lat_arr);
        let retval = unsafe{ converted.as_i32_slice() };
        let combined: Vec<(&i32, &i32)> = retval.iter()
            .zip(correct_values.iter())
            .collect();
        for val in combined.iter() {
            assert_eq!(val.0, val.1);
        };
    }

    #[test]
    fn test_bng_conversion() {
        // verified to be correct at http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        assert_eq!((516276, 173141), convert_bng(-0.32824866, 51.44533267));
    }

    #[test]
    fn test_lonlat_conversion() {
        assert_eq!((-0.32824799370716407, 51.44534026616287), convert_lonlat(516276, 173141));
    }

    #[test]
    #[should_panic]
    fn test_bad_lon() {
        assert_eq!((516276, 173141), convert_bng(181., 51.44533267));
    }

    #[test]
    #[should_panic]
    fn test_bad_lat() {
        assert_eq!((516276, 173141), convert_bng(-0.32824866, -90.01));
    }
}
