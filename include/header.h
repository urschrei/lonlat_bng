/* Generated with cbindgen:0.6.8 */

/* Warning, this file is autogenerated by cbindgen. Don't modify this manually. */

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define AIRY_1830_SEMI_MAJOR 6377563.396

#define AIRY_1830_SEMI_MINOR 6356256.909

#define F0 0.9996012717

#define GRS80_SEMI_MAJOR 6378137

#define GRS80_SEMI_MINOR 6356752.3141

#define MAX_EASTING 700000

#define MAX_LATITUDE 60.84

#define MAX_LONGITUDE 1.78

#define MAX_NORTHING 1250000

#define MIN_LATITUDE 49.96

#define TRUE_ORIGIN_EASTING 400000

#define TY 125.157

typedef struct Array {
    const void *data;
    size_t len;
} Array;

/**
 * A simple struct that can be leaked across the FFI boundary in lieu of an actual tuple
 */
typedef struct ResultTuple {
    Array e;
    Array n;
} ResultTuple;

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_epsg3857_to_wgs84`](fn.convert_epsg3857_to_wgs84.html)
 * # Examples
 * See `lonlat_bng::convert_to_bng_threaded` for examples
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_epsg3857_to_wgs84_threaded(Array x,
                                               Array y);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_etrs89_to_ll`](fn.convert_etrs89_to_ll.html)
 * # Examples
 * See `lonlat_bng::convert_to_bng_threaded` for examples
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_etrs89_to_ll_threaded(Array eastings,
                                          Array northings);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_etrs89_to_osgb36`](fn.convert_etrs89_to_osgb36.html)
 * # Examples
 * See `lonlat_bng::convert_to_bng_threaded` for examples
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_etrs89_to_osgb36_threaded(Array eastings,
                                              Array northings);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36_to_etrs89`](fn.convert_osgb36_to_etrs89.html)
 * # Examples
 * See `lonlat_bng::convert_to_bng_threaded` for examples
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_osgb36_to_etrs89_threaded(Array eastings,
                                              Array northings);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
 * # Examples
 * See `lonlat_bng::convert_to_bng_threaded` for examples
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_osgb36_to_ll_threaded(Array eastings,
                                          Array northings);

/**
 * A threaded, FFI-compatible wrapper for `lonlat_bng::convert_osgb36`
 * # Examples
 * ```
 * use libc;
 * let lon_vec: Vec<f64> = vec![-2.0183041005533306,
 * 0.95511887434519682,
 * 0.44975855518383501,
 * -0.096813621191803811,
 * -0.36807065656416427,
 * 0.63486335458665621];
 * let lat_vec: Vec<f64> = vec![54.589097162646141,
 * 51.560873800587828,
 * 50.431429161121699,
 * 54.535021436247419,
 * 50.839059313135706,
 * 55.412189281234419];
 * let lon_arr = Array {
 * data: lon_vec.as_ptr() as *const libc::c_void,
 * len: lon_vec.len() as libc::size_t,
 * };
 * let lat_arr = Array {
 * data: lat_vec.as_ptr() as *const libc::c_void,
 * len: lat_vec.len() as libc::size_t,
 * };
 * let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
 * ```
 * For an FFI implementation, see the code at [Convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_to_bng_threaded(Array longitudes,
                                    Array latitudes);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_etrs89`](fn.convert_etrs89.html)
 * # Examples
 * See `lonlat_bng::convert_to_bng_threaded` for examples, substituting f64 vectors
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_to_etrs89_threaded(Array longitudes,
                                       Array latitudes);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36_to_ll`](fn.convert_osgb36_to_ll.html)
 * # Examples
 * See [`lonlat_bng::convert_to_bng_threaded`](fn.convert_to_bng_threaded.html)
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_to_lonlat_threaded(Array eastings,
                                       Array northings);

/**
 * A threaded, FFI-compatible wrapper for [`lonlat_bng::convert_osgb36`](fn.convert_osgb36.html)
 * # Examples
 * See [`lonlat_bng::convert_to_bng_threaded`](fn.convert_to_bng_threaded.html) for examples, substituting f64 vectors
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
ResultTuple convert_to_osgb36_threaded(Array longitudes,
                                       Array latitudes);

/**
 * Free memory which Rust has allocated across the FFI boundary (f64 values)
 * # Examples
 * ```
 * # use libc;
 * let lon_vec: Vec<f64> = vec![-2.0183041005533306];
 * let lat_vec: Vec<f64> = vec![54.589097162646141];
 * let lon_arr = Array {
 * data: lon_vec.as_ptr() as *const libc::c_void,
 * len: lon_vec.len() as libc::size_t,
 * };
 * let lat_arr = Array {
 * data: lat_vec.as_ptr() as *const libc::c_void,
 * len: lat_vec.len() as libc::size_t,
 * };
 * let (eastings, northings) = convert_to_bng_threaded(lon_arr, lat_arr);
 * drop_float_array(eastings, northings);
 * ```
 * An example FFI implementation is available at [Convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py), specifically in the `_void_array_to_list` function.
 * # Safety
 * This function is unsafe because it accesses a raw pointer which could contain arbitrary data
 */
void drop_float_array(Array lons,
                      Array lats);
