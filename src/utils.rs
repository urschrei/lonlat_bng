//! This module provides utilities to the conversions module
use ostn15_phf::ostn15_lookup;

use crate::error::{Axis, TransformError};

/// Bounds checking for an input coordinate on a given axis.
///
/// Returns the value unchanged if it lies within the inclusive `bounds`,
/// otherwise a [`TransformError::OutOfBounds`] naming the axis and range.
pub(crate) fn check(value: f64, bounds: (f64, f64), axis: Axis) -> Result<f64, TransformError> {
    if bounds.0 <= value && value <= bounds.1 {
        Ok(value)
    } else {
        Err(TransformError::OutOfBounds {
            axis,
            value,
            min: bounds.0,
            max: bounds.1,
        })
    }
}

/// Round an Easting or Northing coordinate to the nearest millimetre
pub(crate) trait ToMm {
    fn round_to_mm(self) -> f64;
}

impl ToMm for f64 {
    fn round_to_mm(self) -> f64 {
        (self * 1000.).round() / 1000.
    }
}

/// Kahan compensated summation algorithm
/// Reduces accumulated floating-point rounding errors when summing a series of terms
/// See: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
// Only used by the Redfearn inverse series; the Karney path sums differently.
#[cfg_attr(feature = "karney_tm", allow(dead_code))]
pub(crate) fn kahan_sum(terms: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut c = 0.0; // Compensation for lost low-order bits

    for &term in terms {
        let y = term - c; // Subtract the compensation
        let t = sum + y; // Add to running sum
        c = (t - sum) - y; // Recover the low-order part that was lost
        sum = t;
    }

    sum
}

/// Round a float to eight decimal places
pub(crate) fn round_to_eight(x: f64, y: f64) -> (f64, f64) {
    let new_x = (x * 100000000.).round() / 100000000.;
    let new_y = (y * 100000000.).round() / 100000000.;
    (new_x, new_y)
}

/// Look up the OSTN15 shift parameters for a grid intersection, if present.
pub(crate) fn get_ostn_ref(x: i32, y: i32) -> Option<(f64, f64, f64)> {
    let key = x + (y * 701) + 1;
    let result = ostn15_lookup(&key)?;
    Some((result.0, result.1, result.2))
}

// Input values must be valid ETRS89 grid references
// See p20 of the transformation user guide at
// https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html
/// Calculate OSTN15 shifts for a given coordinate
pub(crate) fn ostn15_shifts(x: f64, y: f64) -> Result<(f64, f64, f64), TransformError> {
    let e_index = (x / 1000.) as i32;
    let n_index = (y / 1000.) as i32;

    // eastings and northings of the south-west corner of the cell
    let x0 = e_index * 1000;
    let y0 = n_index * 1000;

    // The shifts for the four corners of the cell. A missing intersection means
    // the point lies outside the OSTN15 grid coverage.
    let coverage_err = || TransformError::OutsideOstn15Coverage {
        easting: x,
        northing: y,
    };

    // intersections
    // this is a 3 x 4 matrix (using column-major order)
    // bottom-left grid intersection
    let s0: (f64, f64, f64) = get_ostn_ref(e_index, n_index).ok_or_else(coverage_err)?;
    // bottom-right
    let s1: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index).ok_or_else(coverage_err)?;
    // top-left
    let s2: (f64, f64, f64) = get_ostn_ref(e_index, n_index + 1).ok_or_else(coverage_err)?;
    // top-right
    let s3: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index + 1).ok_or_else(coverage_err)?;

    // offset within square
    let dx = x - f64::from(x0);
    let dy = y - f64::from(y0);

    let t = dx / 1000.;
    let u = dy / 1000.;

    // Calculation of the weights for each intersection (W)
    // this is a 4 x 1 matrix
    let f0 = (1. - t) * (1. - u);
    let f1 = t * (1. - u);
    let f2 = (1. - t) * u;
    let f3 = t * u;

    // bilinear interpolation, to obtain the actual shifts
    // We could also do a dot product:
    // weights.dot(offsets)
    let se = f0 * s0.0 + f1 * s1.0 + f2 * s2.0 + f3 * s3.0;
    let sn = f0 * s0.1 + f1 * s1.1 + f2 * s2.1 + f3 * s3.1;
    // this isn't needed for this library, since it's a height offset
    let sg = f0 * s0.2 + f1 * s1.2 + f2 * s2.2 + f3 * s3.2;
    Ok((se, sn, sg))
}

/// Detailed shift information including corner shifts and interpolated values
/// Used for testing and validation against OSTN15 reference data
#[cfg(test)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ShiftDetails {
    // Corner shifts (SW, SE, NW, NE)
    pub se0: f64,
    pub sn0: f64,
    pub sg0: f64,
    pub se1: f64,
    pub sn1: f64,
    pub sg1: f64,
    pub se2: f64,
    pub sn2: f64,
    pub sg2: f64,
    pub se3: f64,
    pub sn3: f64,
    pub sg3: f64,
    // Interpolated shifts
    pub se: f64,
    pub sn: f64,
    pub sg: f64,
}

/// Calculate OSTN15 shifts with detailed corner and interpolation data
/// Returns all intermediate values for validation against reference data
#[cfg(test)]
pub(crate) fn ostn15_shifts_detailed(x: f64, y: f64) -> Result<ShiftDetails, TransformError> {
    let e_index = (x / 1000.) as i32;
    let n_index = (y / 1000.) as i32;

    // eastings and northings of the south-west corner of the cell
    let x0 = e_index * 1000;
    let y0 = n_index * 1000;

    // The shifts for the four corners of the cell. A missing intersection means
    // the point lies outside the OSTN15 grid coverage.
    let coverage_err = || TransformError::OutsideOstn15Coverage {
        easting: x,
        northing: y,
    };

    // intersections
    // this is a 3 x 4 matrix (using column-major order)
    // bottom-left grid intersection (SW corner)
    let s0: (f64, f64, f64) = get_ostn_ref(e_index, n_index).ok_or_else(coverage_err)?;
    // bottom-right (SE corner)
    let s1: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index).ok_or_else(coverage_err)?;
    // top-left (NW corner)
    let s2: (f64, f64, f64) = get_ostn_ref(e_index, n_index + 1).ok_or_else(coverage_err)?;
    // top-right (NE corner)
    let s3: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index + 1).ok_or_else(coverage_err)?;

    // offset within square
    let dx = x - f64::from(x0);
    let dy = y - f64::from(y0);

    let t = dx / 1000.;
    let u = dy / 1000.;

    // Calculation of the weights for each intersection (W)
    // this is a 4 x 1 matrix
    let f0 = (1. - t) * (1. - u);
    let f1 = t * (1. - u);
    let f2 = (1. - t) * u;
    let f3 = t * u;

    // bilinear interpolation, to obtain the actual shifts
    let se = f0 * s0.0 + f1 * s1.0 + f2 * s2.0 + f3 * s3.0;
    let sn = f0 * s0.1 + f1 * s1.1 + f2 * s2.1 + f3 * s3.1;
    let sg = f0 * s0.2 + f1 * s1.2 + f2 * s2.2 + f3 * s3.2;

    Ok(ShiftDetails {
        // Corner shifts
        se0: s0.0,
        sn0: s0.1,
        sg0: s0.2,
        se1: s1.0,
        sn1: s1.1,
        sg1: s1.2,
        // Note: CSV uses different corner ordering - swap corners 2 and 3
        se2: s3.0,
        sn2: s3.1,
        sg2: s3.2,
        se3: s2.0,
        sn3: s2.1,
        sg3: s2.2,
        // Interpolated shifts
        se,
        sn,
        sg,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    // original coordinates are 651307.003, 313255.686
    fn test_ostn_hashmap_retrieval() {
        let eastings = 651;
        let northings = 313;
        let expected = (102.787, -78.242, 44.236);
        assert_eq!(expected, get_ostn_ref(eastings, northings).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_failed_ostn_hashmap_retrieval() {
        // we're try!ing this in the shift calculation function, so an Err is fine
        let eastings = 999;
        let northings = 999;
        let expected = (1., 1., 1.);
        assert_eq!(expected, get_ostn_ref(eastings, northings).unwrap());
    }

    #[test]
    fn test_ostn15_shift_incorporation() {
        // these are the input values and corrections on p20-21
        let eastings = 651307.003;
        let northings = 313255.686;
        // Expected values from guide (rounded to mm)
        let expected = (102.801, -78.236, 44.228);
        let result = ostn15_shifts(eastings, northings).unwrap();

        // Round calculated shifts to millimetre precision and compare with published guide values
        assert_eq!(result.0.round_to_mm(), expected.0, "e_shift mismatch");
        assert_eq!(result.1.round_to_mm(), expected.1, "n_shift mismatch");
        assert_eq!(result.2.round_to_mm(), expected.2, "geoid_shift mismatch");
    }

    #[test]
    #[should_panic]
    fn test_min_lon_extents() {
        let max_lon = 1.768960;
        let min_lon = -6.379880;
        // below min_lon
        check(-6.379881, (min_lon, max_lon), Axis::Longitude).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_min_lat_extents() {
        let max_lat = 55.811741;
        let min_lat = 49.871159;
        // below min lat
        check(49.871158, (min_lat, max_lat), Axis::Latitude).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_max_lon_extents() {
        let max_lon = 1.768960;
        let min_lon = -6.379880;
        // above max lon
        check(1.768961, (min_lon, max_lon), Axis::Longitude).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_max_lat_extents() {
        let max_lat = 55.811741;
        let min_lat = 49.871159;
        // above max lat
        check(55.811742, (min_lat, max_lat), Axis::Latitude).unwrap();
    }

    #[test]
    fn test_out_of_bounds_error_detail() {
        // The error names the axis and the violated range.
        let err = check(2.0, (-8.5790, 1.7800), Axis::Longitude).unwrap_err();
        assert_eq!(
            err,
            TransformError::OutOfBounds {
                axis: Axis::Longitude,
                value: 2.0,
                min: -8.5790,
                max: 1.7800,
            }
        );
    }
}
