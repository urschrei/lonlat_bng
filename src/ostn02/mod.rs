// mod shifts;
// use shifts::get_shifts_hashmap;
use super::GRS80_SEMI_MAJOR;
use super::GRS80_SEMI_MINOR;

use super::RAD;
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

// Herbie's going to have a field day with this
pub fn round_to_nearest_mm(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let new_x = (x * 1000.).round() as f64 / 1000.;
    let new_y = (y * 1000.).round() as f64 / 1000.;
    let new_z = (z * 1000.).round() as f64 / 1000.;
    (new_x, new_y, new_z)
}

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
    // only continue if we get Results for the above

    // offset within square
    let dx = *x - (x0 as f64);
    let dy = *y - (y0 as f64);

    // the python script divides by an int here, which = 0 (e.g. 300 / 1000)
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

    let (r_se, r_sn, r_sg) = round_to_nearest_mm(se, sn, sg);
    Ok((r_se, r_sn, r_sg))

}

/// Perform Longitude, Latitude to ETRS89 conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::ostn02::convert_etrs89
/// assert_eq!((651307.003, 313255.686), convert_etrs89(&1.716073973, &52.658007833).unwrap());
#[allow(non_snake_case)]
// See Annexe B (p23) of the transformation user guide for instructions
pub fn convert_etrs89(longitude: &f64, latitude: &f64) -> Result<(f64, f64), ()> {
    // input is restricted to the UK bounding box
    // Convert bounds-checked input to degrees, or return an Err
    let lon_1: f64 = try!(check(*longitude as f32, (MIN_LONGITUDE, MAX_LONGITUDE))) as f64 * RAD;
    let lat_1: f64 = try!(check(*latitude as f32, (MIN_LATITUDE, MAX_LATITUDE))) as f64 * RAD;
    let alt = 0.0;
    // ellipsoid squared eccentricity constant
    let e2 = (WGS84_A.powf(2.) - WGS84_B.powf(2.)) / WGS84_A.powf(2.);
    let n = (WGS84_A - WGS84_B) / (WGS84_A + WGS84_B);
    let phi = RAD * *latitude;
    let lambda = RAD * *longitude;

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

/// Perform Longitude, Latitude to OSGB36 conversion
///
/// # Examples
///
/// ```
/// use lonlat_bng::ostn02::convert_osgb36
/// assert_eq!((651409.792, 313177.448), convert_etrs89(&1.716073973, &52.658007833).unwrap());
#[allow(non_snake_case)]
pub fn convert_osgb36(longitude: &f64, latitude: &f64) -> Result<(f64, f64), ()> {
    // convert input to ETRS89
    let (eastings, northings) = try!(convert_etrs89(longitude, latitude));
    let alt = 0.0;
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
// Convert OSGB36 coordinates to ETRS89 using OSTN02 shifts
// pub fn shift_osgb36_to_etrs89(E: &f64, N: &f64) -> (f64, f64) {
//     let z0 = 0.000;
//     let epsilon = 0.00001;
//     let (dx, dy, dz) = ostn02_shifts(&E, &N);
//     let (x, y, z) = (E - dx, N - dy, dz);
//     let (last_dx, last_dy) = (dx.clone(), dy.clone());
//     while (dx - last_dx).abs() < epsilon && (dy - last_dy).abs() < epsilon {
//         let (dx, dy, dz) = ostn02_shifts(&x, &y);
//         let (x, y) = (E - dx, N - dy);
//         let (last_dx, last_dy) = (dx, dy);
//     }
//     let (x, y, z) = round_to_nearest_mm(E - dx, N - dy, z0 + dz);
//     (x, y)
// }
