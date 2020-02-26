#![doc(html_root_url = "https://urschrei.github.io/lonlat_bng/")]
//! This module provides all the conversion algorithms

// Constants used for coordinate conversions
//
// Ellipsoids
pub const AIRY_1830_SEMI_MAJOR: f64 = 6377563.396;
pub const AIRY_1830_SEMI_MINOR: f64 = 6356256.909;
pub const GRS80_SEMI_MAJOR: f64 = 6378137.000;
pub const GRS80_SEMI_MINOR: f64 = 6356752.3141;
// Northing & easting of true origin (m)
pub const TRUE_ORIGIN_NORTHING: f64 = -100000.;
pub const TRUE_ORIGIN_EASTING: f64 = 400000.;
// For Helmert Transform to OSGB36, translations along the x, y, z axes
// When transforming to WGS84, reverse the signs
pub const TX: f64 = -446.448;
pub const TY: f64 = 125.157;
pub const TZ: f64 = -542.060;
// Rotations along the x, y, z axes, in seconds
pub const RXS: f64 = -0.1502;
pub const RYS: f64 = -0.2470;
pub const RZS: f64 = -0.8421;

// scale factor, referred to as "S" in OSTN15, but some C compilers do not like that
pub const CS: f64 = 20.4894 * 0.000001;
// etc
pub const PI: f64 = f64::consts::PI;
pub const RAD: f64 = PI / 180.;
pub const MAX_EASTING: f64 = 700000.000;
pub const MAX_NORTHING: f64 = 1250000.000;

pub const MIN_LONGITUDE: f64 = -8.5790;
pub const MAX_LONGITUDE: f64 = 1.7800;
pub const MIN_LATITUDE: f64 = 49.922;
pub const MAX_LATITUDE: f64 = 60.8400;

// lon and lat of true origin
const LAM0: f64 = RAD * -2.0;
const PHI0: f64 = RAD * 49.0;

// Easting and Northing of origin
const E0: f64 = TRUE_ORIGIN_EASTING;
const N0: f64 = TRUE_ORIGIN_NORTHING;
// convergence factor
const F0: f64 = 0.9996012717;

use libc::c_double;
use std::f64;
use std::mem;

use crate::utils::check;
use crate::utils::ostn15_shifts;
use crate::utils::round_to_eight;
use crate::utils::ToMm;

/// Calculate the meridional radius of curvature
#[allow(non_snake_case)]
fn curvature(a: f64, f0: f64, e2: f64, lat: f64) -> f64 {
    a * f0 * (1. - e2) * (1. - e2 * lat.sin().powi(2)).powf(-1.5)
}

/// Perform Longitude, Latitude to ETRS89 conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_etrs89
/// assert_eq!((651307.003, 313255.686), convert_etrs89(&1.716073973, &52.658007833).unwrap());
#[allow(non_snake_case)]
// See Annexe B (p23) of the transformation user guide for instructions
pub fn convert_etrs89(longitude: f64, latitude: f64) -> Result<(f64, f64), ()> {
    // Input is restricted to the UK bounding box
    // Convert bounds-checked input to degrees, or return an Err
    let lon_1: f64 = check(longitude, (MIN_LONGITUDE, MAX_LONGITUDE))?.to_radians();
    let lat_1: f64 = check(latitude, (MIN_LATITUDE, MAX_LATITUDE))?.to_radians();
    // ellipsoid squared eccentricity constant
    let e2 = (GRS80_SEMI_MAJOR.powi(2) - GRS80_SEMI_MINOR.powi(2)) / GRS80_SEMI_MAJOR.powi(2);
    let n = (GRS80_SEMI_MAJOR - GRS80_SEMI_MINOR) / (GRS80_SEMI_MAJOR + GRS80_SEMI_MINOR);
    let phi = lat_1;
    let lambda = lon_1;

    let sp2 = phi.sin().powi(2);
    let nu = GRS80_SEMI_MAJOR * F0 * (1. - e2 * sp2).powf(-0.5); // v
    let rho = GRS80_SEMI_MAJOR * F0 * (1. - e2) * (1. - e2 * sp2).powf(-1.5);
    let eta2 = nu / rho - 1.;

    let m = compute_m(phi, GRS80_SEMI_MINOR, n);

    let cp = phi.cos();
    let sp = phi.sin();
    let tp = phi.tan();
    let tp2 = tp.powi(2);
    let tp4 = tp.powi(4);

    let I = m + N0;
    let II = nu / 2. * sp * cp;
    let III = nu / 24. * sp * cp.powi(3) * (5. - tp2 + 9. * eta2);
    let IIIA = nu / 720. * sp * cp.powi(5) * (61. - 58. * tp2 + tp4);

    let IV = nu * cp;
    let V = nu / 6. * cp.powi(3) * (nu / rho - tp2);
    let VI = nu / 120. * cp.powi(5) * (5. - 18. * tp2 + tp4 + 14. * eta2 - 58. * tp2 * eta2);

    let l = lambda - LAM0;
    let north = I + II * l.powi(2) + III * l.powi(4) + IIIA * l.powi(6);
    let east = E0 + IV * l + V * l.powi(3) + VI * l.powi(5);
    Ok((east.round_to_mm(), north.round_to_mm()))
}

/// Perform ETRS89 to OSGB36 conversion, using [OSTN15](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html) data
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_ETRS89_to_OSGB36
/// assert_eq!((651409.792, 313177.448), convert_ETRS89_to_OSGB36(&651307.003, &313255.686).unwrap());
#[allow(non_snake_case)]
pub fn convert_etrs89_to_osgb36(eastings: f64, northings: f64) -> Result<(f64, f64), ()> {
    // ensure that we're within the boundaries
    check(eastings, (0.000, MAX_EASTING))?;
    check(northings, (0.000, MAX_NORTHING))?;
    // obtain OSTN15 corrections, and incorporate
    let (e_shift, n_shift, _) = ostn15_shifts(eastings, northings)?;
    println!("{:?}", (e_shift, n_shift));
    Ok((
        (eastings + e_shift).round_to_mm(),
        (northings + n_shift).round_to_mm(),
    ))
}

/// Perform Longitude, Latitude to OSGB36 conversion, using [OSTN15](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html) data
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_osgb36
/// assert_eq!((651409.792, 313177.448), convert_etrs89(&1.716073973, &52.658007833).unwrap());
#[allow(non_snake_case)]
pub fn convert_osgb36(longitude: f64, latitude: f64) -> Result<(f64, f64), ()> {
    // convert input to ETRS89
    let (eastings, northings) = convert_etrs89(longitude, latitude)?;
    // obtain OSTN15 corrections, and incorporate
    let (e_shift, n_shift, _) = ostn15_shifts(eastings, northings)?;
    Ok((
        (eastings + e_shift).round_to_mm(),
        (northings + n_shift).round_to_mm(),
    ))
}

// Intermediate calculation used for lon, lat to ETRS89 and reverse conversion
fn compute_m(phi: f64, b: f64, n: f64) -> f64 {
    let p_plus = phi + PHI0;
    let p_minus = phi - PHI0;

    b * F0
        * ((1. + n * (1. + 5. / 4. * n * (1. + n))) * p_minus
            - 3. * n * (1. + n * (1. + 7. / 8. * n)) * p_minus.sin() * p_plus.cos()
            + (15. / 8. * n * (n * (1. + n))) * (2. * p_minus).sin() * (2. * p_plus).cos()
            - 35. / 24. * n.powi(3) * (3. * p_minus).sin() * (3. * p_plus).cos())
}

// Easting and Northing to Lon, Lat conversion using a Helmert transform
// Note that either GRS80 or Airy 1830 ellipsoids can be passed
#[allow(non_snake_case)]
fn convert_to_ll(eastings: f64, northings: f64, ell_a: f64, ell_b: f64) -> Result<(f64, f64), ()> {
    // ensure that we're within the boundaries
    check(eastings, (0.000, MAX_EASTING))?;
    check(northings, (0.000, MAX_NORTHING))?;
    // ellipsoid squared eccentricity constant
    let a = ell_a;
    let b = ell_b;
    let e2 = (a.powi(2) - b.powi(2)) / a.powi(2);
    let n = (a - b) / (a + b);

    let dN = northings - N0;
    let mut phi = PHI0 + dN / (a * F0);
    let mut m = compute_m(phi, b, n);
    while (dN - m) >= 0.00001 {
        m = compute_m(phi, b, n);
        phi += (dN - m) / (a * F0);
    }
    let sp2 = phi.sin().powi(2);
    let nu = a * F0 * (1. - e2 * sp2).powf(-0.5);
    let rho = a * F0 * (1. - e2) * (1. - e2 * sp2).powf(-1.5);
    let eta2 = nu / rho - 1.;

    let tp = phi.tan();
    let tp2 = tp.powi(2);
    let tp4 = tp.powi(4);

    let VII = tp / (2. * rho * nu);
    let VIII = tp / (24. * rho * nu.powi(3)) * (5. + 3. * tp2 + eta2 - 9. * tp2 * eta2);
    let IX = tp / (720. * rho * nu.powi(5)) * (61. + 90. * tp2 + 45. * tp4);

    let sp = 1.0 / phi.cos();
    let tp6 = tp4 * tp2;

    let X = sp / nu;
    let XI = sp / (6. * nu.powi(3)) * (nu / rho + 2. * tp2);
    let XII = sp / (120. * nu.powi(5)) * (5. + 28. * tp2 + 24. * tp4);
    let XIIA = sp / (5040. * nu.powi(7)) * (61. + 662. * tp2 + 1320. * tp4 + 720. * tp6);

    let e = eastings - E0;

    phi = phi - VII * e.powi(2) + VIII * e.powi(4) - IX * e.powi(6);
    let mut lambda = LAM0 + X * e - XI * e.powi(3) + XII * e.powi(5) - XIIA * e.powi(7);

    phi = phi.to_degrees();
    lambda = lambda.to_degrees();
    Ok(round_to_eight(lambda, phi))
}

/// Convert ETRS89 coordinates to Lon, Lat
#[allow(non_snake_case)]
pub fn convert_etrs89_to_ll(E: f64, N: f64) -> Result<(f64, f64), ()> {
    // ETRS89 uses the WGS84 / GRS80 ellipsoid constants
    convert_to_ll(E, N, GRS80_SEMI_MAJOR, GRS80_SEMI_MINOR)
}

/// Convert OSGB36 coordinates to Lon, Lat using OSTN15 data
#[allow(non_snake_case)]
pub fn convert_osgb36_to_ll(E: f64, N: f64) -> Result<(f64, f64), ()> {
    // Apply reverse OSTN15 adustments
    let epsilon = 0.009;
    let (mut dx, mut dy, _) = ostn15_shifts(E, N)?;
    let (mut x, mut y) = (E - dx, N - dy);
    let (mut last_dx, mut last_dy) = (dx, dy);
    let mut res;
    loop {
        res = ostn15_shifts(x, y)?;
        dx = res.0;
        dy = res.1;
        x = E - dx;
        y = N - dy;
        // If the difference […] is more than 0.00010m (User Guide, p15)
        // TODO: invert this logic
        if (dx - last_dx).abs() < epsilon && (dy - last_dy).abs() < epsilon {
            break;
        }
        last_dx = dx;
        last_dy = dy;
    }
    let x = (E - dx).round_to_mm();
    let y = (N - dy).round_to_mm();
    // We've converted to ETRS89, so we need to use the WGS84/ GRS80 ellipsoid constants
    convert_to_ll(x, y, GRS80_SEMI_MAJOR, GRS80_SEMI_MINOR)
}

/// Convert OSGB36 coordinates to ETRS89 using OSTN15 data
#[allow(non_snake_case)]
pub fn convert_osgb36_to_etrs89(E: f64, N: f64) -> Result<(f64, f64), ()> {
    // Apply reverse OSTN15 adustments
    let epsilon = 0.00001;
    let (mut dx, mut dy, _) = ostn15_shifts(E, N)?;
    let (mut x, mut y) = (E - dx, N - dy);
    let (mut last_dx, mut last_dy) = (dx, dy);
    let mut res;
    loop {
        res = ostn15_shifts(x, y)?;
        dx = res.0;
        dy = res.1;
        x = E - dx;
        y = N - dy;
        if (dx - last_dx).abs() < epsilon && (dy - last_dy).abs() < epsilon {
            break;
        }
        last_dx = dx;
        last_dy = dy;
    }
    let x = (E - dx).round_to_mm();
    let y = (N - dy).round_to_mm();
    Ok((x, y))
}

/// **THIS FUNCTION IS DEPRECATED**
///
/// Perform Longitude, Latitude to British National Grid conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_bng;
/// assert_eq!((516276.000, 173141.000), convert_bng(&-0.32824866, &51.44533267).unwrap());
#[allow(non_snake_case)]
#[allow(dead_code)]
pub fn convert_bng(longitude: f64, latitude: f64) -> Result<(c_double, c_double), ()> {
    // input is restricted to the UK bounding box
    // Convert bounds-checked input to degrees, or return an Err
    let lon_1: f64 = check(longitude, (MIN_LONGITUDE, MAX_LONGITUDE))?.to_radians();
    let lat_1: f64 = check(latitude, (MIN_LATITUDE, MAX_LATITUDE))?.to_radians();
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
    let x_2 = tx + (1. + CS) * x_1 + -rz * y_1 + ry * z_1;
    let y_2 = ty + rz * x_1 + (1. + CS) * y_1 + -rx * z_1;
    let z_2 = tz + -ry * x_1 + rx * y_1 + (1. + CS) * z_1;

    // The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
    let a = AIRY_1830_SEMI_MAJOR;
    let b = AIRY_1830_SEMI_MINOR;
    // The eccentricity of the Airy 1830 ellipsoid
    let e2 = 1. - b.powi(2) / a.powi(2);
    let p = (x_2.powi(2) + y_2.powi(2)).sqrt();
    // Initial value
    let mut lat = z_2.atan2(p * (1. - e2));
    let mut latold = 2. * PI;
    // this is cheating, but not sure how else to initialise nu
    let mut nu: f64 = 1.;
    // Latitude is obtained by iterative procedure
    while (lat - latold).abs() > (10. as f64).powi(-16) {
        mem::swap(&mut lat, &mut latold);
        nu = a / (1. - e2 * latold.sin().powi(2)).sqrt();
        lat = (z_2 + e2 * nu * latold.sin()).atan2(p);
    }
    let lon = y_2.atan2(x_2);
    // Latitude of true origin (radians)
    let lat0 = 49. * PI / 180.;
    // Longitude of true origin and central meridian (radians)
    let lon0 = -2. * PI / 180.;
    // Northing & easting of true origin (m)
    let n = (a - b) / (a + b);
    // Meridional radius of curvature
    let rho = curvature(a, F0, e2, lat);
    let eta2 = nu * F0 / rho - 1.;

    let M1 = (1. + n + (5. / 4.) * n.powi(2) + (5. / 4.) * n.powi(3)) * (lat - lat0);
    let M2 = (3. * n + 3. * n.powi(2) + (21. / 8.) * n.powi(3))
        * ((lat.sin() * lat0.cos()) - (lat.cos() * lat0.sin()))
            .ln_1p()
            .exp_m1() * (lat + lat0).cos();
    let M3 = ((15. / 8.) * n.powi(2) + (15. / 8.) * n.powi(3))
        * (2. * (lat - lat0)).sin()
        * (2. * (lat + lat0)).cos();
    let M4 = (35. / 24.) * n.powi(3) * (3. * (lat - lat0)).sin() * (3. * (lat + lat0)).cos();
    let M = b * F0 * (M1 - M2 + M3 - M4);

    let I = M + N0;
    let II = nu * F0 * lat.sin() * lat.cos() / 2.;
    let III = nu * F0 * lat.sin() * lat.cos().powi(3) * (5. - lat.tan().powi(2) + 9. * eta2) / 24.;
    let IIIA = nu
        * F0
        * lat.sin()
        * lat.cos().powi(5)
        * (61. - 58. * lat.tan().powi(2) + lat.tan().powi(4)) / 720.;
    let IV = nu * F0 * lat.cos();
    let V = nu * F0 * lat.cos().powi(3) * (nu / rho - lat.tan().powi(2)) / 6.;
    let VI = nu * F0 * lat.cos().powi(5)
        * (5. - 18. * lat.tan().powi(2) + lat.tan().powi(4) + 14. * eta2
            - 58. * eta2 * lat.tan().powi(2)) / 120.;
    let N =
        I + II * (lon - lon0).powi(2) + III * (lon - lon0).powi(4) + IIIA * (lon - lon0).powi(6);
    let E = E0 + IV * (lon - lon0) + V * (lon - lon0).powi(3) + VI * (lon - lon0).powi(5);

    Ok((E.round_to_mm(), N.round_to_mm()))
}

/// **THIS FUNCTION IS DEPRECATED**
///
/// Perform British National Grid Eastings, Northings to Longitude, Latitude conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_lonlat;
/// assert_eq!((-0.328248, 51.44534), convert_lonlat(&516276, &173141));
#[allow(non_snake_case)]
#[allow(dead_code)]
pub fn convert_lonlat(easting: f64, northing: f64) -> Result<(f64, f64), ()> {
    // The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
    let a = AIRY_1830_SEMI_MAJOR;
    let b = AIRY_1830_SEMI_MINOR;
    // Scale factor on the central meridian
    // let F0: f64 = 0.9996012717;
    // Latitude of true origin (radians)
    let lat0 = 49. * PI / 180.;
    // Longitude of true origin and central meridian (radians)
    let lon0 = -2. * PI / 180.;
    // Northing & easting of true origin (m)
    // Eccentricity squared
    let e2 = 1. - b.powi(2) / a.powi(2);
    let n = (a - b) / (a + b);

    let mut lat = lat0;
    let mut M: f64 = 0.0;
    while (northing - N0 - M) >= 0.00001 {
        lat += (northing - N0 - M) / (a * F0);
        let M1 = (1. + n + (5. / 4.) * n.powi(3) + (5. / 4.) * n.powi(3)) * (lat - lat0);
        let M2 = (3. * n + 3. * n.powi(2) + (21. / 8.) * n.powi(3))
            * ((lat.sin() * lat0.cos()) - (lat.cos() * lat0.sin()))
                .ln_1p()
                .exp_m1() * (lat + lat0).cos();
        let M3 = ((15. / 8.) * n.powi(2) + (15. / 8.) * n.powi(3))
            * (2. * (lat - lat0)).sin()
            * (2. * (lat + lat0)).cos();
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
    let VIII = lat.tan() / (24. * rho * nu.powi(3))
        * (5. + 3. * lat.tan().powi(2) + eta2 - 9. * lat.tan().powi(2) * eta2);
    let IX = lat.tan() / (720. * rho * nu.powi(5))
        * (61. + 90. * lat.tan().powi(2) + 45. * lat.tan().powi(4));
    let X = secLat / nu;
    let XI = secLat / (6. * nu.powi(3)) * (nu / rho + 2. * lat.tan().powi(2));
    let XII =
        secLat / (120. * nu.powi(5)) * (5. + 28. * lat.tan().powi(2) + 24. * lat.tan().powi(4));
    let XIIA = secLat / (5040. * nu.powi(7))
        * (61. + 662. * lat.tan().powi(2) + 1320. * lat.tan().powi(4) + 720. * lat.tan().powi(6));
    let dE = easting - E0;
    // These are on the wrong ellipsoid currently: Airy1830 (Denoted by _1)
    let lat_1 = lat - VII * dE.powi(2) + VIII * dE.powi(4) - IX * dE.powi(6);
    let lon_1 = lon0 + X * dE - XI * dE.powi(3) + XII * dE.powi(5) - XIIA * dE.powi(7);

    // We want to convert to the GRS80 ellipsoid
    // First, convert to cartesian from spherical polar coordinates
    let H = 0.;
    let x_1 = (nu / F0 + H) * lat_1.cos() * lon_1.cos();
    let y_1 = (nu / F0 + H) * lat_1.cos() * lon_1.sin();
    let z_1 = ((1. - e2) * nu / F0 + H) * lat_1.sin();

    // Perform Helmert transform (to go between Airy 1830 (_1) and GRS80 (_2))
    let minus_s = -CS; // The scale factor -1
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
    let mut lat = z_2.atan2(p * (1. - e2_2));
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
    Ok(round_to_eight(lon, lat))
}

/// Convert Web Mercator (from Google Maps or Bing Maps) to WGS84
// from https://alastaira.wordpress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection/
pub fn convert_epsg3857_to_wgs84(x: f64, y: f64) -> Result<(f64, f64), ()> {
    let lon = (x / 20037508.34) * 180.;
    let mut lat = (y / 20037508.34) * 180.;
    lat = 180. / PI * (2. * (lat * PI / 180.).exp().atan() - PI / 2.);
    Ok((lon, lat))
}

#[cfg(test)]
mod tests {

    use super::convert_bng;
    use super::convert_epsg3857_to_wgs84;
    use super::convert_etrs89;
    use super::convert_etrs89_to_ll;
    use super::convert_etrs89_to_osgb36;
    use super::convert_lonlat;
    use super::convert_osgb36;
    use super::convert_osgb36_to_ll;

    #[test]
    fn test_gmaps_to_wgs() {
        let x = -626172.1357121646;
        let y = 6887893.4928337997;
        let expected = (-5.625000000783013, 52.48278022732355);
        assert_eq!(expected, convert_epsg3857_to_wgs84(x, y).unwrap());
    }

    #[test]
    fn test_convert_osgb36_to_ll() {
        // Caister Water Tower, with OSTN15 corrections applied. See p23
        // Final Lon, Lat rounded to eight decimal places
        // p20 gives the correct lon, lat as (1.716073973, 52.658007833)
        let easting = 651409.804;
        let northing = 313177.450;
        let expected = (1.71607397, 52.65800783);
        assert_eq!(expected, convert_osgb36_to_ll(easting, northing).unwrap());
    }

    #[test]
    fn test_convert_etrs89_to_ll() {
        // Caister Water Tower, ETRS89. See p20
        let easting = 651307.003;
        let northing = 313255.686;
        assert_eq!(
            (1.71607397, 52.65800783),
            convert_etrs89_to_ll(easting, northing).unwrap()
        );
    }

    #[test]
    fn test_etrs89_conversion() {
        // these are the input values and intermediate result in the example on p20–23
        let longitude = 1.716073973;
        let latitude = 52.658007833;
        let expected = (651307.003, 313255.686);
        assert_eq!(expected, convert_etrs89(longitude, latitude).unwrap());
    }

    #[test]
    fn test_osgb36_conversion() {
        // these are the input values and final result in the example on p20–23
        let longitude = 1.716073973;
        let latitude = 52.658007833;
        let expected = (651409.804, 313177.450);
        assert_eq!(expected, convert_osgb36(longitude, latitude).unwrap());
    }

    #[test]
    fn test_etrs89_to_osgb36_conversion() {
        // these are the input values and final result in the example on p20–23
        let eastings = 651307.003;
        let northings = 313255.686;
        let expected = (651409.804, 313177.450);
        assert_eq!(
            expected,
            convert_etrs89_to_osgb36(eastings, northings).unwrap()
        );
    }

    #[test]
    #[should_panic]
    fn test_bad_max_easting() {
        let max_easting = 700001.000;
        let max_northing = 1250000.000;
        // above max lat
        convert_etrs89_to_osgb36(max_easting, max_northing).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_bad_max_northing() {
        let max_easting = 700000.000;
        let max_northing = 1250000.001;
        // above max lat
        convert_etrs89_to_osgb36(max_easting, max_northing).unwrap();
    }

    #[test]
    fn test_bng_conversion() {
        // verified to be correct at http://www.bgs.ac.uk/data/webservices/convertForm.cfm
        assert_eq!(
            (516275.973, 173141.092),
            convert_bng(-0.32824866, 51.44533267).unwrap()
        );
    }

    #[test]
    fn test_lonlat_conversion() {
        let res = convert_lonlat(516276.000, 173141.000).unwrap();
        // We shouldn't really be using error margins, but it should be OK because
        // neither number is zero, or very close to, and on opposite sides of zero
        // epsilon is .000001 here, because BNG coords are 6 digits, so
        // we should be fine if the error is in the 7th digit (i.e. < epsilon)
        // http://floating-point-gui.de/errors/comparison/
        assert!(((res.0 - -0.328248269313) / -0.328248269313).abs() < 0.000001);
        assert!(((res.1 - 51.4453318435) / 51.4453318435).abs() < 0.000001);
    }

    #[test]
    // TrainTrick reported that this coordinate doesn't converge at an epsilon of 0.00001
    fn test_traintrick() {
        let res = convert_osgb36_to_ll(515415.0, 202612.0).unwrap();
        assert_eq!(res.0, -0.33093489);
        assert_eq!(res.1, 51.71038497);
    }

    #[test]
    #[should_panic]
    fn test_bad_lon() {
        assert_eq!(
            (516276.000, 173141.000),
            convert_bng(181., 51.44533267).unwrap()
        );
    }

    #[test]
    #[should_panic]
    fn test_bad_lat() {
        assert_eq!(
            (516276.000, 173141.000),
            convert_bng(-0.32824866, -90.01).unwrap()
        );
    }

}
