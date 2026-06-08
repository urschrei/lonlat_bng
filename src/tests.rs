//! Tests module for lonlat_bng
//!
//! This module contains unit tests for the lonlat_bng library.

#![cfg(test)]

use std::collections::HashMap;

#[cfg(not(feature = "karney_tm"))]
use crate::conversions::GRS80_SEMI_MAJOR;
#[cfg(not(feature = "karney_tm"))]
use crate::conversions::GRS80_SEMI_MINOR;
#[cfg(not(feature = "karney_tm"))]
use crate::conversions::convert_to_ll;
use crate::conversions::osgb36_to_etrs89_iterative_detailed;

use super::*;

/// Represents whether this row is an iteration step or the final result
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IterationValue {
    /// An iteration step (1, 2, 3, ...)
    Iteration(i32),
    /// The final result row with Lat/Long coordinates
    Result,
}

/// Test data row from OSTN15_OSGM15_TestOutput_OSGBtoETRS.csv
/// Represents a single iteration or final result of the OSGB36→ETRS89 conversion
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Osgb36ToEtrs89TestRow {
    /// Test point identifier (e.g., "TP01")
    pub point_id: String,

    /// Iteration number or final result indicator
    pub iteration: IterationValue,

    /// ETRS89 Easting (metres) for iteration rows, or Latitude (decimal degrees) for RESULT rows
    pub etrs_east_or_lat: f64,

    /// ETRS89 Northing (metres) for iteration rows, or Longitude (decimal degrees) for RESULT rows
    pub etrs_north_or_long: f64,

    /// ETRS89 Height (metres)
    pub etrs_height: f64,

    // OSTN15 grid shift values at the four corners of the grid cell:
    // Corner 0: South-West, Corner 1: South-East, Corner 2: North-West, Corner 3: North-East
    /// Easting shift at SW corner (metres)
    pub se0: f64,
    /// Northing shift at SW corner (metres)
    pub sn0: f64,
    /// Geoid-ellipsoid separation shift at SW corner (metres)
    pub sg0: f64,

    /// Easting shift at SE corner (metres)
    pub se1: f64,
    /// Northing shift at SE corner (metres)
    pub sn1: f64,
    /// Geoid-ellipsoid separation shift at SE corner (metres)
    pub sg1: f64,

    /// Easting shift at NW corner (metres)
    pub se2: f64,
    /// Northing shift at NW corner (metres)
    pub sn2: f64,
    /// Geoid-ellipsoid separation shift at NW corner (metres)
    pub sg2: f64,

    /// Easting shift at NE corner (metres)
    pub se3: f64,
    /// Northing shift at NE corner (metres)
    pub sn3: f64,
    /// Geoid-ellipsoid separation shift at NE corner (metres)
    pub sg3: f64,

    /// Bilinearly interpolated easting shift (metres)
    pub se: f64,
    /// Bilinearly interpolated northing shift (metres)
    pub sn: f64,
    /// Geoid-ellipsoid separation shift (metres)
    pub sg: f64,
    osgb36_northing: f64,
    osgb36_easting: f64,
}

/// Load test data from the revised.csv file
/// Returns a vector of all data rows (iterations and RESULT rows), skipping empty rows
pub(crate) fn load_ostn15_test_data() -> Vec<Osgb36ToEtrs89TestRow> {
    let csv_data = include_str!("../test_inputs/OSGB36_to_ETRS.csv");
    let mut rows = Vec::new();

    for (line_num, line) in csv_data.lines().enumerate() {
        // Skip header row (line 0)
        if line_num == 0 {
            continue;
        }

        // Skip empty or whitespace-only lines
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with(',') {
            continue;
        }

        // Parse the CSV line
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 20 {
            panic!(
                "Expected 20 fields on line {}, got {}",
                line_num + 1,
                fields.len()
            );
        }

        // Parse point_id (may have leading space or BOM)
        let point_id = fields[0].trim().trim_start_matches('\u{feff}').to_string();

        // Parse iteration
        let iteration_str = fields[1].trim();
        let iteration = if iteration_str == "RESULT" {
            IterationValue::Result
        } else {
            let num = iteration_str.parse::<i32>().unwrap_or_else(|_| {
                panic!(
                    "Failed to parse iteration '{}' on line {}",
                    iteration_str,
                    line_num + 1
                )
            });
            IterationValue::Iteration(num)
        };

        // Parse all f64 fields
        let parse_f64 = |idx: usize, name: &str| -> f64 {
            fields[idx].trim().parse::<f64>().unwrap_or_else(|_| {
                panic!(
                    "Failed to parse {} '{}' on line {}",
                    name,
                    fields[idx],
                    line_num + 1
                )
            })
        };

        rows.push(Osgb36ToEtrs89TestRow {
            point_id,
            osgb36_easting: f64::NAN,
            osgb36_northing: f64::NAN,
            iteration,
            etrs_east_or_lat: parse_f64(2, "ETRSEast or Lat"),
            etrs_north_or_long: parse_f64(3, "ETRSNorth or Long"),
            etrs_height: parse_f64(4, "ETRSHeight"),
            se0: parse_f64(5, "Se0"),
            sn0: parse_f64(6, "Sn0"),
            sg0: parse_f64(7, "Sg0"),
            se1: parse_f64(8, "Se1"),
            sn1: parse_f64(9, "Sn1"),
            sg1: parse_f64(10, "Sg1"),
            se2: parse_f64(11, "Se2"),
            sn2: parse_f64(12, "Sn2"),
            sg2: parse_f64(13, "Sg2"),
            se3: parse_f64(14, "Se3"),
            sn3: parse_f64(15, "Sn3"),
            sg3: parse_f64(16, "Sg3"),
            se: parse_f64(17, "Se"),
            sn: parse_f64(18, "Sn"),
            sg: parse_f64(19, "Sg"),
        });
    }

    rows
}

/// Strip a leading UTF-8 BOM and surrounding whitespace from a CSV field.
fn clean(field: &str) -> &str {
    field.trim().trim_start_matches('\u{feff}')
}

/// OSTN15 developer-pack ETRS89 -> OSGB36 inputs: (PointID, latitude, longitude).
fn dev_pack_etrs_inputs() -> Vec<(String, f64, f64)> {
    let csv = include_str!("../test_inputs/OSTN15_OSGM15_TestInput_ETRStoOSGB.csv");
    csv.lines()
        .skip(1)
        .filter_map(|line| {
            let f: Vec<&str> = line.split(',').collect();
            if f.len() < 3 {
                return None;
            }
            Some((
                clean(f[0]).to_string(),
                clean(f[1]).parse().unwrap(),
                clean(f[2]).parse().unwrap(),
            ))
        })
        .collect()
}

/// OSTN15 developer-pack OSGB36 eastings/northings keyed by PointID.
///
/// Both the ETRS->OSGB *output* file and the OSGB->ETRS *input* file share the
/// first three columns (PointID, easting, northing), so this reads either.
fn dev_pack_osgb_en(csv: &str) -> HashMap<String, (f64, f64)> {
    csv.lines()
        .skip(1)
        .filter_map(|line| {
            let f: Vec<&str> = line.split(',').collect();
            if f.len() < 3 {
                return None;
            }
            Some((
                clean(f[0]).to_string(),
                (clean(f[1]).parse().unwrap(), clean(f[2]).parse().unwrap()),
            ))
        })
        .collect()
}

#[test]
/// Forward ETRS89 lon/lat -> OSGB36 against the developer-pack reference output.
///
/// The published values are the OS exact-method results rounded to the millimetre.
/// Both the Redfearn and Karney projections agree to last-digit precision; the
/// residual is bounded by the millimetre rounding of the published coordinates.
fn dev_pack_forward_etrs_to_osgb() {
    let expected = dev_pack_osgb_en(include_str!(
        "../test_inputs/OSTN15_OSGM15_TestOutput_ETRStoOSGB.csv"
    ));
    // 1 mm rounding tolerance per axis, with a little slack for fp subtraction.
    const TOL: f64 = 0.0015;
    let mut failures = Vec::new();
    for (id, lat, lon) in dev_pack_etrs_inputs() {
        let (e, n) = convert_osgb36(lon, lat).unwrap();
        let (ee, en) = expected[&id];
        if (e - ee).abs() > TOL || (n - en).abs() > TOL {
            failures.push(format!(
                "{id}: got ({e:.3}, {n:.3}) expected ({ee:.3}, {en:.3})"
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} forward mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn test_load_ostn15_test_data() {
    let data = load_ostn15_test_data();

    // Should have 40 test points with multiple iterations each
    // Expected: TP01-TP40, each with iterations + 1 RESULT row
    // TP01: 3 iterations + 1 RESULT = 4 rows
    // TP02: 2 iterations + 1 RESULT = 3 rows
    // etc.
    assert!(!data.is_empty(), "Should load at least some data");

    // Check first row (TP01, iteration 1)
    let first = &data[0];
    assert_eq!(first.point_id, "TP01");
    assert_eq!(first.iteration, IterationValue::Iteration(1));
    assert_eq!(first.etrs_east_or_lat, 91399.9984);
    assert_eq!(first.etrs_north_or_long, 11399.9999);

    // Check a RESULT row exists
    let result_rows: Vec<_> = data
        .iter()
        .filter(|r| r.iteration == IterationValue::Result)
        .collect();
    assert_eq!(
        result_rows.len(),
        40,
        "Should have 40 RESULT rows (one per test point)"
    );

    // Check first RESULT row (TP01)
    let tp01_result = result_rows
        .iter()
        .find(|r| r.point_id == "TP01")
        .expect("Should have TP01 RESULT row");
    assert_eq!(tp01_result.etrs_east_or_lat, 49.92226394); // Latitude
    assert_eq!(tp01_result.etrs_north_or_long, -6.29977752); // Longitude
}

#[test]
/// Reverse OSGB36 -> ETRS89 pipeline validated against the developer-pack vectors.
///
/// The OSGB36 eastings/northings are taken from the OSGB->ETRS *input* file, which
/// is the basis OS used to generate `OSGB36_to_ETRS.csv`. For every point the
/// per-iteration ETRS89 coordinates and the twelve corner shifts plus the
/// interpolated shifts are checked against the published values; these depend only
/// on the OSTN15 grid interpolation, so they hold for either projection.
///
/// The final lon/lat is compared against the published RESULT row, but only on the
/// default (Redfearn) build: the published RESULT is itself a product of the OS
/// truncated inverse, and the Karney projection deliberately diverges from it at
/// the far west of the grid (validated instead by the `dev_pack_karney_*` tests).
fn test_osgb36_to_etrs89_iterations_detailed() {
    // Helper functions to round values to match CSV precision
    let round_to_3dp = |x: f64| -> f64 { (x * 1000.).round() / 1000. };
    let round_to_4dp = |x: f64| -> f64 { (x * 10000.).round() / 10000. };

    // OSGB36 input coordinates keyed by PointID (OSGB->ETRS reverse-test input).
    let osgb_in = dev_pack_osgb_en(include_str!(
        "../test_inputs/OSTN15_OSGM15_TestInput_OSGBtoETRS.csv"
    ));

    let data = load_ostn15_test_data();

    // Group test data by point_id
    let mut test_points: HashMap<String, Vec<Osgb36ToEtrs89TestRow>> = HashMap::new();
    for row in data {
        test_points
            .entry(row.point_id.clone())
            .or_default()
            .push(row);
    }

    let mut failures = Vec::new();

    // Process all 40 test points
    for idx in 0..40 {
        let point_id = format!("TP{:02}", idx + 1);
        let (osgb36_e, osgb36_n) = match osgb_in.get(&point_id) {
            Some(&en) => en,
            None => {
                failures.push(format!("{}: No OSGB->ETRS input found", point_id));
                continue;
            }
        };
        let rows = match test_points.get(&point_id) {
            Some(rows) => rows,
            None => {
                failures.push(format!("{}: No test data found", point_id));
                continue;
            }
        };

        // Separate iteration rows from RESULT row
        let iteration_rows: Vec<_> = rows
            .iter()
            .filter(|r| matches!(r.iteration, IterationValue::Iteration(_)))
            .collect();
        let result_row = rows.iter().find(|r| r.iteration == IterationValue::Result);

        let result_row = match result_row {
            Some(r) => r,
            None => {
                failures.push(format!("{}: No RESULT row found", point_id));
                continue;
            }
        };

        // Run the detailed conversion
        let conversion_result = osgb36_to_etrs89_iterative_detailed(osgb36_e, osgb36_n);
        let (final_etrs89_e, final_etrs89_n, iterations) = match conversion_result {
            Ok(result) => result,
            Err(_) => {
                failures.push(format!("{}: Conversion failed", point_id));
                continue;
            }
        };

        // The number of iterations to convergence is an implementation detail of
        // the convergence test and may legitimately differ from the published
        // table by one step, so the count itself is not asserted; only the
        // overlapping iterations are compared.

        // Validate each iteration
        let min_iterations = iterations.len().min(iteration_rows.len());
        for iter_idx in 0..min_iterations {
            let expected = &iteration_rows[iter_idx];
            let actual = &iterations[iter_idx];
            let iter_num = iter_idx + 1;

            // Compare ETRS89 coordinates. These intermediates carry cross-file
            // rounding from the dev pack (the reverse-test inputs are rounded
            // differently from the coordinates OS used to generate the iteration
            // table), so allow a 1 mm tolerance rather than exact equality.
            const ITER_TOL: f64 = 0.001;
            if (actual.etrs89_e - expected.etrs_east_or_lat).abs() > ITER_TOL {
                failures.push(format!(
                    "{} Iteration {}: ETRS89 Easting mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_4dp(actual.etrs89_e),
                    expected.etrs_east_or_lat
                ));
            }
            if (actual.etrs89_n - expected.etrs_north_or_long).abs() > ITER_TOL {
                failures.push(format!(
                    "{} Iteration {}: ETRS89 Northing mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_4dp(actual.etrs89_n),
                    expected.etrs_north_or_long
                ));
            }

            // Compare all 12 corner shifts (round to 3 decimal places to match CSV precision)
            if round_to_3dp(actual.shifts.se0) != expected.se0 {
                failures.push(format!(
                    "{} Iteration {}: se0 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.se0),
                    expected.se0
                ));
            }
            if round_to_3dp(actual.shifts.sn0) != expected.sn0 {
                failures.push(format!(
                    "{} Iteration {}: sn0 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sn0),
                    expected.sn0
                ));
            }
            if round_to_3dp(actual.shifts.sg0) != expected.sg0 {
                failures.push(format!(
                    "{} Iteration {}: sg0 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sg0),
                    expected.sg0
                ));
            }
            if round_to_3dp(actual.shifts.se1) != expected.se1 {
                failures.push(format!(
                    "{} Iteration {}: se1 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.se1),
                    expected.se1
                ));
            }
            if round_to_3dp(actual.shifts.sn1) != expected.sn1 {
                failures.push(format!(
                    "{} Iteration {}: sn1 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sn1),
                    expected.sn1
                ));
            }
            if round_to_3dp(actual.shifts.sg1) != expected.sg1 {
                failures.push(format!(
                    "{} Iteration {}: sg1 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sg1),
                    expected.sg1
                ));
            }
            if round_to_3dp(actual.shifts.se2) != expected.se2 {
                failures.push(format!(
                    "{} Iteration {}: se2 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.se2),
                    expected.se2
                ));
            }
            if round_to_3dp(actual.shifts.sn2) != expected.sn2 {
                failures.push(format!(
                    "{} Iteration {}: sn2 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sn2),
                    expected.sn2
                ));
            }
            if round_to_3dp(actual.shifts.sg2) != expected.sg2 {
                failures.push(format!(
                    "{} Iteration {}: sg2 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sg2),
                    expected.sg2
                ));
            }
            if round_to_3dp(actual.shifts.se3) != expected.se3 {
                failures.push(format!(
                    "{} Iteration {}: se3 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.se3),
                    expected.se3
                ));
            }
            if round_to_3dp(actual.shifts.sn3) != expected.sn3 {
                failures.push(format!(
                    "{} Iteration {}: sn3 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sn3),
                    expected.sn3
                ));
            }
            if round_to_3dp(actual.shifts.sg3) != expected.sg3 {
                failures.push(format!(
                    "{} Iteration {}: sg3 mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_3dp(actual.shifts.sg3),
                    expected.sg3
                ));
            }

            // Compare interpolated shifts (round to 4 decimal places to match CSV precision)
            if round_to_4dp(actual.shifts.se) != expected.se {
                failures.push(format!(
                    "{} Iteration {}: se mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_4dp(actual.shifts.se),
                    expected.se
                ));
            }
            if round_to_4dp(actual.shifts.sn) != expected.sn {
                failures.push(format!(
                    "{} Iteration {}: sn mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_4dp(actual.shifts.sn),
                    expected.sn
                ));
            }
            if round_to_4dp(actual.shifts.sg) != expected.sg {
                failures.push(format!(
                    "{} Iteration {}: sg mismatch: {} != {}",
                    point_id,
                    iter_num,
                    round_to_4dp(actual.shifts.sg),
                    expected.sg
                ));
            }
        }

        // Final lon/lat vs the published RESULT row. The published RESULT is a
        // product of the OS truncated inverse, so only the default (Redfearn)
        // projection is expected to reproduce it; the Karney build is validated by
        // the round-trip tests below.
        #[cfg(not(feature = "karney_tm"))]
        {
            let (lon, lat) = match convert_to_ll(
                final_etrs89_e,
                final_etrs89_n,
                GRS80_SEMI_MAJOR,
                GRS80_SEMI_MINOR,
            ) {
                Ok(result) => result,
                Err(_) => {
                    failures.push(format!("{}: Lon/Lat conversion failed", point_id));
                    continue;
                }
            };
            let expected_lat = result_row.etrs_east_or_lat;
            let expected_lon = result_row.etrs_north_or_long;
            // The published RESULT is given to 8 dp (~0.6 mm); compare on the ground.
            const R: f64 = 6_378_137.0;
            let dlat_mm = (lat - expected_lat).to_radians() * R * 1000.0;
            let dlon_mm =
                (lon - expected_lon).to_radians() * R * expected_lat.to_radians().cos() * 1000.0;
            let ground_mm = (dlat_mm * dlat_mm + dlon_mm * dlon_mm).sqrt();
            if ground_mm > 1.5 {
                failures.push(format!(
                    "{}: final lon/lat off by {:.4} mm",
                    point_id, ground_mm
                ));
            }
        }
        #[cfg(feature = "karney_tm")]
        {
            let _ = (final_etrs89_e, final_etrs89_n, result_row);
        }
    }

    // Assert at the end if there were any failures
    assert!(
        failures.is_empty(),
        "{} validation failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
