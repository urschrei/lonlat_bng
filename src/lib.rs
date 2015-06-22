//! The `lonlat_bng` crate provides a function that converts decimal longitude
//! and latitude coordinates into British National Grid Coordinates
//!
//! Examples
//!
//!```
//! assert_eq!((516276, 173141), lonlat_bng::convert(-0.32824866, 51.44533267));
//!```

#[allow(non_snake_case)]
use std::f32::consts;
use std::mem;

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
#[no_mangle]
pub extern fn convert(input_lon: f32, input_lat: f32) -> (i32, i32) {

    // if not all([0 <= input_lat <= 90, -180 <= input_lon <= 180]):
    match input_lon {
        -180.0...180.0 => input_lon,
        _ => panic!("Out of bounds! Longitude must be between -180 and 180.")
    };
      match input_lat {
        -90.0...90.0 => input_lat,
        _ => panic!("Out of bounds! Latitude must be between -90 and 90.")
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
    return (round(E) as i32, round(N) as i32);
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
