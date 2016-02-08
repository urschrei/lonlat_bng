#![doc(html_root_url = "https://urschrei.github.io/lonlat_bng/")]
//! This module provides high-quality transforms to OSGB36 using OSTN02 data
//! As such, it should be suitable for use in surveying and construction applications
//!
//!
//!
//!
//!
use super::GRS80_SEMI_MAJOR;
use super::GRS80_SEMI_MINOR;

use super::RAD;
use super::DAR;
use super::MIN_LONGITUDE;
use super::MAX_LONGITUDE;
use super::MIN_LATITUDE;
use super::MAX_LATITUDE;

use super::TRUE_ORIGIN_EASTING;
use super::TRUE_ORIGIN_NORTHING;

const WGS84_A: f64 = GRS80_SEMI_MAJOR;
const WGS84_B: f64 = GRS80_SEMI_MINOR;

// lon and lat of true origin
const LAM0: f64 = RAD * -2.0;
const PHI0: f64 = RAD * 49.0;

// Easting and Northing of origin
const E0: f64 = TRUE_ORIGIN_EASTING;
const N0: f64 = TRUE_ORIGIN_NORTHING;
// convergence factor
const F0: f64 = 0.9996012717;

const MIN_X_SHIFT: f64 = 86.275;
const MIN_Y_SHIFT: f64 = -81.603;
const MIN_Z_SHIFT: f64 = 43.982;

use ostn02_phf::ostn02_lookup;
use super::check;

// TODO Herbie's going to have a field day with this
/// Round Easting, Northing, and Orthometric height to the nearest millimetre
pub fn round_to_nearest_mm(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let new_x = (x * 1000.).round() as f64 / 1000.;
    let new_y = (y * 1000.).round() as f64 / 1000.;
    let new_z = (z * 1000.).round() as f64 / 1000.;
    (new_x, new_y, new_z)
}

// TODO Herbie's going to have a field day with this
/// Round a float to six decimal places
pub fn round_to_six(x: f64, y:f64) -> (f64, f64) {
    let new_x = (x * 1000000.).round() as f64 / 1000000.;
    let new_y = (y * 1000000.).round() as f64 / 1000000.;
    (new_x, new_y)
}

// Try to get OSTN02 shift parameters, and calculate offsets
pub fn get_ostn_ref(x: &i32, y: &i32) -> Result<(f64, f64, f64), ()> {
    let key = format!("{:03x}{:03x}", y, x);
    // some or None, so try! this
    let lookup = ostn02_lookup(&*key);
    let intermediate = lookup.ok_or(());
    let result = intermediate.unwrap();
    Ok((result.0 as f64 / 1000. + MIN_X_SHIFT,
        result.1 as f64 / 1000. + MIN_Y_SHIFT,
        result.2 as f64 / 1000. + MIN_Z_SHIFT))

}

// Input values must be valid ETRS89 grid references
// See p20 of the transformation user guide at
// https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html
pub fn ostn02_shifts(x: &f64, y: &f64) -> Result<(f64, f64, f64), ()> {
    let e_index = (*x / 1000.) as i32;
    let n_index = (*y / 1000.) as i32;

    // eastings and northings of the south-west corner of the cell
    let x0 = e_index * 1000;
    let y0 = n_index * 1000;

    // The easting, northing and geoid shifts for the four corners of the cell
    // any of these could be Err, so use try!
    let s0: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 0), &(n_index + 0)));
    let s1: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 1), &(n_index + 0)));
    let s2: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 0), &(n_index + 1)));
    let s3: (f64, f64, f64) = try!(get_ostn_ref(&(e_index + 1), &(n_index + 1)));

    // offset within square
    let dx = *x - (x0 as f64);
    let dy = *y - (y0 as f64);

    let t = dx / 1000.;
    let u = dy / 1000.;

    let f0 = (1. - t as f64) * (1. - u as f64);
    let f1 = t as f64 * (1. - u as f64);
    let f2 = (1. - t as f64) * u as f64;
    let f3 = t as f64 * u as f64;

    // bilinear interpolation, to obtain the actual shifts
    let se = f0 * s0.0 + f1 * s1.0 + f2 * s2.0 + f3 * s3.0;
    let sn = f0 * s0.1 + f1 * s1.1 + f2 * s2.1 + f3 * s3.1;
    let sg = f0 * s0.2 + f1 * s1.2 + f2 * s2.2 + f3 * s3.2;

    Ok(round_to_nearest_mm(se, sn, sg))

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
pub fn convert_etrs89(longitude: &f64, latitude: &f64) -> Result<(f64, f64), ()> {
    // input is restricted to the UK bounding box
    // Convert bounds-checked input to degrees, or return an Err
    let lon_1: f64 = try!(check(*longitude, (MIN_LONGITUDE, MAX_LONGITUDE))) as f64 * RAD;
    let lat_1: f64 = try!(check(*latitude, (MIN_LATITUDE, MAX_LATITUDE))) as f64 * RAD;
    // ellipsoid squared eccentricity constant
    let e2 = (WGS84_A.powf(2.) - WGS84_B.powf(2.)) / WGS84_A.powf(2.);
    let n = (WGS84_A - WGS84_B) / (WGS84_A + WGS84_B);
    let phi = lat_1;
    let lambda = lon_1;

    let sp2 = phi.sin().powf(2.);
    let nu = WGS84_A * F0 * (1. - e2 * sp2).powf(-0.5); // v
    let rho = WGS84_A * F0 * (1. - e2) * (1. - e2 * sp2).powf(-1.5);
    let eta2 = nu / rho - 1.;

    let m = compute_m(&phi, &WGS84_B, &n);

    let cp = phi.cos();
    let sp = phi.sin();
    let tp = phi.tan();
    let tp2 = tp.powf(2.);
    let tp4 = tp.powf(4.);

    let I = m + N0;
    let II = nu / 2. * sp * cp;
    let III = nu / 24. * sp * cp.powf(3.) * (5. - tp2 + 9. * eta2);
    let IIIA = nu / 720. * sp * cp.powf(5.) * (61. - 58. * tp2 + tp4);

    let IV = nu * cp;
    let V = nu / 6. * cp.powf(3.) * (nu / rho - tp2);
    let VI = nu / 120. * cp.powf(5.) * (5. - 18. * tp2 + tp4 + 14. * eta2 - 58. * tp2 * eta2);

    let l = lambda - LAM0;
    let north = I + II * l.powf(2.) + III * l.powf(4.) + IIIA * l.powf(6.);
    let east = E0 + IV * l + V * l.powf(3.) + VI * l.powf(5.);
    let (rounded_eastings, rounded_northings, _) = round_to_nearest_mm(east, north, 1.00);
    Ok((rounded_eastings, rounded_northings))
}

/// Perform ETRS89 to OSGB36 conversion, using [OSTN02](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html) data
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_ETRS89_to_OSGB36
/// assert_eq!((651409.792, 313177.448), convert_ETRS89_to_OSGB36(&651307.003, &313255.686).unwrap());
#[allow(non_snake_case)]
pub fn convert_ETRS89_to_OSGB36(eastings: &f64, northings: &f64) -> Result<(f64, f64), ()> {
    // obtain OSTN02 corrections, and incorporate
    let (e_shift, n_shift, _) = try!(ostn02_shifts(&eastings, &northings));
    let (shifted_e, shifted_n) = (eastings + e_shift, northings + n_shift);
    Ok((shifted_e, shifted_n))

}

/// Perform Longitude, Latitude to OSGB36 conversion, using [OSTN02](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html) data
///
/// # Examples
///
/// ```
/// use lonlat_bng::convert_osgb36
/// assert_eq!((651409.792, 313177.448), convert_etrs89(&1.716073973, &52.658007833).unwrap());
#[allow(non_snake_case)]
pub fn convert_osgb36(longitude: &f64, latitude: &f64) -> Result<(f64, f64), ()> {
    // convert input to ETRS89
    let (eastings, northings) = try!(convert_etrs89(longitude, latitude));
    // obtain OSTN02 corrections, and incorporate
    let (e_shift, n_shift, _) = try!(ostn02_shifts(&eastings, &northings));
    let (shifted_e, shifted_n) = (eastings + e_shift, northings + n_shift);
    Ok((shifted_e, shifted_n))
}

fn compute_m(phi: &f64, b: &f64, n: &f64) -> f64 {
    let p_plus = *phi + PHI0;
    let p_minus = *phi - PHI0;

    let result = *b * F0 *
                 ((1. + *n * (1. + 5. / 4. * *n * (1. + *n))) * p_minus -
                  3. * *n * (1. + *n * (1. + 7. / 8. * *n)) * p_minus.sin() * p_plus.cos() +
                  (15. / 8. * *n * (*n * (1. + *n))) * (2. * p_minus).sin() * (2. * p_plus).cos() -
                  35. / 24. * n.powf(3.) * (3. * p_minus).sin() * (3. * p_plus).cos());
    result
}

#[allow(non_snake_case)]
fn convert_ETRS89_to_ll(eastings: &f64, northings: &f64) -> Result<(f64, f64), ()> {
    // ellipsoid squared eccentricity constant
    let e2 = (WGS84_A.powf(2.) - WGS84_B.powf(2.)) / WGS84_A.powf(2.);
    let n = (WGS84_A - WGS84_B) / (WGS84_A + WGS84_B);

    let dN = *northings - N0;
    let mut phi = PHI0 + dN / (WGS84_A * F0);
    let mut m = compute_m(&phi, &WGS84_B, &n);
    while (dN - m) >= 0.001 {
        phi = phi + (dN - m) / (WGS84_A * F0);
        m = compute_m(&phi, &WGS84_B, &n);    
    }
    let sp2 = phi.sin().powf(2.);
    let nu = WGS84_A * F0 * (1. - e2 * sp2).powf(-0.5);
    let rho = WGS84_A * F0 * (1. - e2) * (1. - e2 * sp2).powf(-1.5);
    let eta2 = nu / rho - 1.;

    let tp = phi.tan();
    let tp2 = tp.powf(2.);
    let tp4 = tp.powf(4.);

    let VII = tp / (2. * rho * nu);
    let VIII = tp / (24. * rho * nu.powf(3.)) * (5. + 3. * tp2 + eta2 - 9. * tp2 * eta2);
    let IX = tp / (720. * rho * nu.powf(5.)) * (61. + 90. * tp2 + 45. * tp4);

    let sp = 1.0 / phi.cos();
    let tp6 = tp4 * tp2;

    let X = sp / nu;
    let XI = sp / (6. * nu.powf(3.)) * (nu / rho + 2. * tp2);
    let XII = sp / ( 120. * nu.powf(5.)) * (5.+ 28. * tp2 + 24.*tp4);
    let XIIA = sp / (5040. * nu.powf(7.)) * (61. + 662. * tp2 + 1320. * tp4 + 720. * tp6);

    let e = *eastings - E0;

    phi = phi - VII * e.powf(2.) + VIII * e.powf(4.) - IX * e.powf(6.);
    let mut lambda = LAM0 + X  *e - XI * e.powf(3.) + XII * e.powf(5.) - XIIA * e.powf(7.);

    phi = phi * DAR;
    lambda = lambda * DAR;
    Ok(round_to_six(lambda, phi))
}

/// Convert OSGB36 coordinates to Lon, Lat using OSTN02 shifts
// pub fn shift_osgb36_to_etrs89(E: &f64, N: &f64) -> Result<(f64, f64), ()> {
//     let z0 = 0.000;
//     let epsilon = 0.00001;
//     let (dx, dy, dz) = try!(ostn02_shifts(&E, &N));
//     let (x, y, z) = (E - dx, N - dy, dz);
//     let (last_dx, last_dy) = (dx.clone(), dy.clone());
//     while (dx - last_dx).abs() < epsilon && (dy - last_dy).abs() < epsilon {
//         let (dx, dy, dz) = try!(ostn02_shifts(&x, &y));
//         let (x, y) = (E - dx, N - dy);
//         let (last_dx, last_dy) = (dx, dy);
//     }
//     let (x, y, _) = round_to_nearest_mm(E - dx, N - dy, z0 + dz);
//     // this function returns a Result
//     convert_ETRS89_to_ll(&x, &y)
// }

#[cfg(test)]
mod tests {
    use super::get_ostn_ref;
    use super::ostn02_shifts;
    use super::convert_etrs89;
    use super::convert_osgb36;
    use super::convert_ETRS89_to_OSGB36;
    use super::convert_ETRS89_to_ll;
    // use super::shift_osgb36_to_etrs89;

    // #[test]
    // fn test_shift_osgb36_to_etrs89() {
    //     let easting = 651409.792;
    //     let northing = 313177.448;
    //     assert_eq!((1.716074, 52.658008), shift_osgb36_to_etrs89(&easting, &northing).unwrap());
    // }

    #[test]
    fn test_convert_ETRS89_to_ll() {
        let easting = 651307.003;
        let northing = 313255.686;
        assert_eq!((1.716074, 52.658008), convert_ETRS89_to_ll(&easting, &northing).unwrap());
    }

    #[test]
    fn test_etrs89_conversion() {
        // these are the input values and intermediate result in the example on p20–21
        let longitude = 1.716073973;
        let latitude = 52.658007833;
        let expected = (651307.003, 313255.686);
        assert_eq!(expected, convert_etrs89(&longitude, &latitude).unwrap());
    }

    #[test]
    fn test_osgb36_conversion() {
        // these are the input values and final result in the example on p20–21
        let longitude = 1.716073973;
        let latitude = 52.658007833;
        let expected = (651409.792, 313177.448);
        assert_eq!(expected, convert_osgb36(&longitude, &latitude).unwrap());
    }

    #[test]
    fn test_etrs89_to_osgb36_conversion() {
        // these are the input values and final result in the example on p20–21
        let eastings = 651307.003;
        let northings = 313255.686;
        let expected = (651409.792, 313177.448);
        assert_eq!(expected,
                   convert_ETRS89_to_OSGB36(&eastings, &northings).unwrap());
    }

    #[test]
    // original coordinates are 651307.003, 313255.686
    fn test_ostn_hashmap_retrieval() {
        let eastings = 651;
        let northings = 313;
        let expected = (102.775, -78.244, 44.252);
        assert_eq!(expected, get_ostn_ref(&eastings, &northings).unwrap());
    }

    #[test]
    fn test_ostn02_shift_incorporation() {
        // these are the input values and corrections on p20-21
        let eastings = 651307.003;
        let northings = 313255.686;
        let expected = (102.789, -78.238, 44.244);
        assert_eq!(expected, ostn02_shifts(&eastings, &northings).unwrap());
    }

}
