use rayon::prelude::*;

use super::NAN;
use super::convert_epsg3857_to_wgs84;
use super::convert_etrs89;
use super::convert_etrs89_to_ll;
use super::convert_etrs89_to_osgb36;
use super::convert_osgb36;
use super::convert_osgb36_to_etrs89;
use super::convert_osgb36_to_ll;

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_bng_threaded_vec<'a>(
    longitudes: &'a mut [f64],
    latitudes: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_to_osgb36_threaded_vec(longitudes, latitudes)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_to_lonlat_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_osgb36_to_ll_threaded_vec(eastings, northings)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89`](fn.convert_etrs89.html)
pub fn convert_to_etrs89_threaded_vec<'a>(
    longitudes: &'a mut [f64],
    latitudes: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
pub fn convert_to_osgb36_threaded_vec<'a>(
    longitudes: &'a mut [f64],
    latitudes: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(longitudes, latitudes, convert_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
pub fn convert_etrs89_to_osgb36_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_etrs89_to_osgb36)
}

/// A threaded wrapper for [`lonlat_bng::convert_etrs89_to_ll`](fn.convert_etrs89_to_ll.html)
pub fn convert_etrs89_to_ll_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_etrs89_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_etrs89`](fn.convert_osgb36_to_etrs89.html)
pub fn convert_osgb36_to_etrs89_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_etrs89)
}

/// A threaded wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
pub fn convert_osgb36_to_ll_threaded_vec<'a>(
    eastings: &'a mut [f64],
    northings: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(eastings, northings, convert_osgb36_to_ll)
}

/// A threaded wrapper for [`lonlat_bng::convert_epsg3857_to_wgs84`](fn.convert_epsg3857_to_wgs84.html)
pub fn convert_epsg3857_to_wgs84_threaded_vec<'a>(
    x: &'a mut [f64],
    y: &'a mut [f64],
) -> (&'a mut [f64], &'a mut [f64]) {
    convert_vec_direct(x, y, convert_epsg3857_to_wgs84)
}

// Generic function which applies conversion functions to paired slices within threads
fn convert_vec_direct<'a>(
    ex: &'a mut [f64],
    ny: &'a mut [f64],
    func: impl Fn(f64, f64) -> Result<(f64, f64), ()> + Sync,
) -> (&'a mut [f64], &'a mut [f64]) {
    ex.par_iter_mut().zip(ny.par_iter_mut()).for_each(|p| {
        match func(*p.0, *p.1) {
            // mutate values, or assign default error values
            Ok(res) => {
                *p.0 = res.0;
                *p.1 = res.1;
            }
            Err(_) => {
                *p.0 = NAN;
                *p.1 = NAN;
            }
        }
    });
    (ex, ny)
}
