//! Error types for coordinate transformations.

use std::error::Error;
use std::fmt;

/// The coordinate axis that failed a bounds check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// Geographic longitude (degrees).
    Longitude,
    /// Geographic latitude (degrees).
    Latitude,
    /// Projected easting (metres).
    Easting,
    /// Projected northing (metres).
    Northing,
}

impl fmt::Display for Axis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Axis::Longitude => "longitude",
            Axis::Latitude => "latitude",
            Axis::Easting => "easting",
            Axis::Northing => "northing",
        };
        f.write_str(name)
    }
}

/// An error arising during a coordinate transformation.
///
/// The crate's conversion functions are fallible because their inputs may lie
/// outside the area the OSTN15 transformation covers, or because the reverse
/// transformation may fail to converge. Each variant carries the offending
/// coordinate(s) so the cause can be reported or logged.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum TransformError {
    /// An input coordinate fell outside the valid range for the conversion.
    OutOfBounds {
        /// Which coordinate axis was out of range.
        axis: Axis,
        /// The offending value.
        value: f64,
        /// The inclusive lower bound of the valid range.
        min: f64,
        /// The inclusive upper bound of the valid range.
        max: f64,
    },
    /// The coordinate lies outside the area covered by the OSTN15 grid (for
    /// example, offshore), so no shift could be interpolated.
    OutsideOstn15Coverage {
        /// Easting of the query point (metres).
        easting: f64,
        /// Northing of the query point (metres).
        northing: f64,
    },
    /// The iterative OSGB36 to ETRS89 conversion did not converge within the
    /// iteration limit.
    NonConvergent {
        /// Easting of the OSGB36 point (metres).
        easting: f64,
        /// Northing of the OSGB36 point (metres).
        northing: f64,
    },
}

impl fmt::Display for TransformError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TransformError::OutOfBounds {
                axis,
                value,
                min,
                max,
            } => write!(
                f,
                "{axis} {value} is outside the valid range [{min}, {max}]"
            ),
            TransformError::OutsideOstn15Coverage { easting, northing } => write!(
                f,
                "coordinate ({easting}, {northing}) lies outside OSTN15 grid coverage"
            ),
            TransformError::NonConvergent { easting, northing } => write!(
                f,
                "OSGB36 to ETRS89 conversion did not converge for ({easting}, {northing})"
            ),
        }
    }
}

impl Error for TransformError {}
