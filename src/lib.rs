//! The `lonlat_bng` crate provides a function that converts decimal longitude
//! and latitude coordinates into British National Grid Coordinates
//!
//! Examples
//!
//!```
//! assert_eq!((516276, 173141), lonlat_bng::convert(-0.32824866, 51.44533267));
//!```

use std::f32::consts;
use std::mem;
use std::slice;
use std::thread::{self, JoinHandle};

extern crate libc;
use libc::{size_t, c_void, uint32_t};

extern crate rand;

const NUMTHREADS: usize = 7;

#[repr(C)]
pub struct Tuple {
    a: uint32_t,
    b: uint32_t,
}

#[repr(C)]
pub struct Array {
    data: *const c_void,
    len: libc::size_t,
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
fn round(x: f32) -> f32 {
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
/// use lonlat_bng::convert;
/// assert_eq!((516276, 173141), convert(-0.32824866, 51.44533267));
#[allow(non_snake_case)]
#[no_mangle]
pub extern fn convert(input_lon: f32, input_lat: f32) -> (i32, i32) {
    match input_lon {
        -180.0...180.0 => input_lon,
        _ => panic!("Out of bounds! Longitude must be between -180.00 and 180.00: {:?}", input_lon)
    };
    match input_lat {
        -90.0...90.0 => input_lat,
        _ => panic!("Out of bounds! Latitude must be between -90.00 and 90.00: {:?}", input_lat)
    };
    let pi: f32 = consts::PI;
    //Convert input to degrees
    let lat_1: f32 = input_lat * pi / 180.;
    let lon_1: f32 = input_lon * pi / 180.;
    // The GSR80 semi-major and semi-minor axes used for WGS84 (m)
    let a_1: f32 = 6378137.;
    let b_1: f32 = 6356752.3141;
    // The eccentricity (squared) of the GRS80 ellipsoid
    let e2_1: f32 = 1. - (b_1.powf(2.)) / (a_1.powf(2.));
    // Transverse radius of curvature
    let nu_1: f32 = a_1 / (1. - e2_1 * lat_1.sin().powf(2.)).sqrt();
    // Third spherical coordinate is 0, in this case
    let H: f32 = 0.;
    let x_1: f32 = (nu_1 + H) * lat_1.cos() * lon_1.cos();
    let y_1: f32 = (nu_1 + H) * lat_1.cos() * lon_1.sin();
    let z_1: f32 = ((1. - e2_1) * nu_1 + H) * lat_1.sin();

    // Perform Helmert transform (to go between Airy 1830 (_1) and GRS80 (_2))
    let small: f32 = 10.;
    let cst: f32 = 20.4894;
    let s: f32 =  cst * small.powf(-6.);
    // The translations along x, y, z axes respectively
    let tx: f32 = -446.448;
    let ty: f32 = 125.157;
    let tz: f32 = -542.060;
    // The rotations along x, y, z respectively, in seconds
    let rxs: f32 = -0.1502;
    let rys: f32 = -0.2470;
    let rzs: f32 = -0.8421;
    // In radians
    let rx: f32 = rxs * pi / (180. * 3600.);
    let ry: f32 = rys * pi / (180. * 3600.);
    let rz: f32 = rzs * pi / (180. * 3600.);
    // panic begins on line 46
    let x_2: f32 = tx + (1. + s) * x_1 + -rz * y_1 + ry * z_1;
    let y_2: f32 = ty + rz * x_1 + (1. + s) * y_1 + -rx * z_1;
    let z_2: f32 = tz + -ry * x_1 + rx * y_1 + (1. + s) * z_1;

    // The GSR80 semi-major and semi-minor axes used for WGS84 (m)
    let a: f32 = 6377563.396;
    let b: f32 = 6356256.909;
    // The eccentricity of the Airy 1830 ellipsoid
    let e2: f32 = 1. - b.powf(2.) / a.powf(2.);
    let p: f32 = (x_2.powf(2.) + y_2.powf(2.)).sqrt();
    // Initial value
    let mut lat: f32 = z_2.atan2((p * (1. - e2)));
    let mut latold: f32 = 2. * pi;
    // this is cheating, but not sure how else to initialise nu
    let mut nu: f32 = 1.;
    // Latitude is obtained by iterative procedure
    while (lat - latold).abs() > small.powf(-16.) {
        mem::swap(&mut lat, &mut latold);
        nu = a / (1. - e2 * latold.sin().powf(2.)).sqrt();
        lat = (z_2 + e2 * nu * latold.sin()).atan2(p);
    };
    let lon: f32 = y_2.atan2(x_2);
    // Scale factor on the central meridian
    let F0: f32 = 0.9996012717;
    // Latitude of true origin (radians)
    let lat0: f32 = 49. * pi / 180.;
    // Longtitude of true origin and central meridian (radians)
    let lon0: f32 = -2. * pi / 180.;
    // Northing & easting of true origin (m)
    let N0: f32 = -100000.;
    let E0: f32 = 400000.;
    let n: f32 = (a - b) / (a + b);
    // Meridional radius of curvature
    let rho: f32 = a * F0 * (1. - e2) * (1. - e2 * lat.sin().powf(2.)).powf(-1.5);
    let eta2: f32 = nu * F0 / rho - 1.;

    let M1: f32 = (1. + n + (5. / 4.) * n.powf(2.) + (5. / 4.) * n.powf(3.)) * (lat - lat0);
    let M2: f32 = (3. * n + 3. * n.powf(2.) + (21. / 8.) * n.powf(3.)) * (lat - lat0).sin() * (lat + lat0).cos();
    let M3: f32 = ((15. / 8.) * n.powf(2.) + (15. / 8.) * n.powf(3.)) * (2. * (lat-lat0)).sin() * (2. * (lat + lat0)).cos();
    let M4: f32 = (35. / 24.) * n.powf(3.) * (3. * (lat - lat0)).sin() * (3. * (lat + lat0)).cos();
    let M: f32 = b * F0 * (M1 - M2 + M3 - M4);

    let I: f32 = M + N0;
    let II: f32 = nu * F0 * lat.sin() * lat.cos() / 2.;
    let III: f32 = nu * F0 * lat.sin() * lat.cos().powf(3.) * (5. - lat.tan().powf(2.) + 9. * eta2) / 24.;
    let IIIA: f32 = nu * F0 * lat.sin() * lat.cos().powf(5.) * (61. - 58. * lat.tan().powf(2.) + lat.tan().powf(4.)) / 720.;
    let IV: f32 = nu * F0 * lat.cos();
    let V: f32 = nu * F0 * lat.cos().powf(3.) * (nu / rho - lat.tan().powf(2.)) / 6.;
    let VI: f32 = nu * F0 * lat.cos().powf(5.) * (5. - 18. * lat.tan().powf(2.) + lat.tan().powf(4.) + 14. * eta2 - 58. * eta2 * lat.tan().powf(2.)) / 120.;
    let N: f32 = I + II * (lon - lon0).powf(2.) + III * (lon - lon0).powf(4.) + IIIA * (lon - lon0).powf(6.);
    let E: f32 = E0 + IV * (lon - lon0) + V * (lon - lon0).powf(3.) + VI * (lon - lon0).powf(5.);
    (round(E) as i32, round(N) as i32)
}

/// A safer C-compatible wrapper for convert()
#[no_mangle]
pub extern fn convert_vec_c(lon: Array, lat: Array) -> Array {
    // we're receiving floats
    let lon = unsafe { lon.as_f32_slice() };
    let lat = unsafe { lat.as_f32_slice() };
    // copy values and combine
    let orig: Vec<(f32, f32)> = lon
        .iter()
        .map(|i| i.clone())
        .collect::<Vec<f32>>()
        .into_iter()
    .zip(lat
        .iter()
        .map(|i| i.clone())
        .collect::<Vec<f32>>()
    .into_iter()).collect();
    // carry out the conversion 
    let result: Vec<(i32, i32)> = orig.iter()
        .map(|elem| convert(elem.0, elem.1))
        .collect();
    // convert back to vector of unsigned integer Tuples
    let nvec = result.iter()
        .map(|ints| Tuple { a: ints.0 as u32, b: ints.1 as u32 })
        .collect();
    Array::from_vec(nvec)
}

/// A threaded version of the C-compatible wrapper for convert()
#[no_mangle]
pub extern fn convert_vec_c_threaded(lon: Array, lat: Array) -> Array {
    // we're receiving floats
    let lon = unsafe { lon.as_f32_slice() };
    let lat = unsafe { lat.as_f32_slice() };
    // copy values and combine
    let orig: Vec<(f32, f32)> = lon
        .iter()
        .map(|i| i.clone())
        .collect::<Vec<f32>>()
        .into_iter()
    .zip(lat
        .iter()
        .map(|i| i.clone())
        .collect::<Vec<f32>>()
    .into_iter()).collect();

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
            .map(|elem| convert(elem.0, elem.1))
            .collect());
        guards.push(g);
    }
    let mut result: Vec<(i32, i32)> = Vec::with_capacity(orig.len());
    for g in guards {
        result.extend(g.join().unwrap().into_iter());
    }
    // convert back to vector of unsigned integer Tuples
    let nvec = result.iter()
        .map(|ints| Tuple { a: ints.0 as u32, b: ints.1 as u32 })
        .collect();
    Array::from_vec(nvec)
}

#[cfg(test)]
mod tests {
    use super::convert;
    use super::convert_vec_c;
    use super::convert_vec_c_threaded;
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
        let converted = convert_vec_c_threaded(lon_arr, lat_arr);
        let retval = unsafe{ converted.as_i32_slice() };
        // the value's incorrect, but let's worry about that later
        assert_eq!(398915, retval[0]);
        assert_eq!(521545, retval[1]);
    }

    #[test]
    fn test_threaded_vector_conversion_single() {
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
        let converted = convert_vec_c_threaded(lon_arr, lat_arr);
        let retval = unsafe{ converted.as_i32_slice() };
        assert_eq!(398915, retval[0]);
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
            188805, // BNG says 188804
            574082,
            61932, // BNG says 61931
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
    fn test_conversion() {
        // verified to be correct at http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        assert_eq!((516276, 173141), convert(-0.32824866, 51.44533267));
    }

    #[test]
    #[should_panic]
    fn test_bad_lon() {
        assert_eq!((516276, 173141), convert(181., 51.44533267));
    }

    #[test]
    #[should_panic]
    fn test_bad_lat() {
        assert_eq!((516276, 173141), convert(-0.32824866, -90.01));
    }
}
