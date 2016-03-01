[![Build Status](https://travis-ci.org/urschrei/lonlat_bng.png?branch=master)](https://travis-ci.org/urschrei/lonlat_bng) [![](https://img.shields.io/crates/v/lonlat_bng.svg)](https://crates.io/crates/lonlat_bng) [![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](license.txt)  

# Introduction
An attempt at speeding up the conversion between decimal longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates, using a Rust binary with FFI. Conversions use a standrad Helmert transform, with the addition of OSTN02 corrections for [accuracy](#accuracy), where appropriate.

# Motivation
Python (etc.) is relatively slow; this type of conversion is usually carried out in bulk, so an order-of-magnitude improvement using FFI could save precious minutes.

# Accuracy
The Helmert transforms and their threaded and vectorised versions are accurate to within around 5 metres, and are **not suitable** for calculations or conversions used in e.g. surveying.    
Thus, we use the OSTN02 transform, which adjusts for local variation within the Terrestrial Reference Frame by incorporating OSTN02 data. [See here](http://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) for more information.  

The OSTN02-enabled functions are:

- convert_osgb36
- convert_etrs89_to_osgb36
- convert_to_osgb36_threaded ← FFI
- convert_to_osgb36_threaded_vec
- convert_etrs89_to_osgb36_threaded ← FFI
- convert_etrs89_to_osgb36_threaded_vec
- convert_osgb36_to_ll_threaded ← FFI
- convert_osgb36_to_ll_threaded_vec
- convert_osgb36_to_etrs89_threaded ← FFI
- convert_osgb36_to_etrs89_threaded_vec

- convert_bng_threaded (an alias for convert_osgb36_threaded)
- convert_bng_threaded_vec ← FFI version of the above

- convert_lonlat_threaded (an alias for convert_osgb36_to_ll)
- convert_lonlat_threaded_vec ← FFI version of the above

[![OSTN02](ostn002_s.gif)]( "OSTN02")

# Library Use
## As a Rust Library
Add the following to your `Cargo.toml` (the latest version is displayed on the second badge at the top of this screen)  

    lonlat_bng = "x.x.x"

Full library documentation is available [here](http://urschrei.github.io/lonlat_bng/)  

**Note that `lon`, `lat` coordinates outside the [UK bounding box](http://spatialreference.org/ref/epsg/27700/) will be transformed to `(NAN, NAN)`, which cannot be mapped.**  

The functions exposed by the library can be found [here](http://urschrei.github.io/lonlat_bng/lonlat_bng/index.html#functions)

## As an FFI Library
The FFI C-compatible functions exposed by the library are:  
`convert_to_bng_threaded(Array, Array) -> Array`  
`convert_to_lonlat_threaded(Array, Array) -> Array`  

`convert_to_osgb36_threaded(Array, Array) -> Array`  
`convert_to_etrs89_threaded(Array, Array) -> Array)`  
`convert_osgb36_to_ll_threaded(Array, Array) -> Array`  
`convert_etrs89_to_ll_threaded(Array, Array) -> Array`  

`convert_etrs89_to_osgb36_threaded(Array, Array) -> Array`  
`convert_osgb36_to_etrs89_threaded(Array, Array) -> Array`  

`convert_epsg3857_to_wgs84_threaded(Array, Array) -> Array`  

### FFI and Memory Management
If your library, module, or script uses the FFI functions, it **must** implement `drop_float_array`. **Failing to do so may result in memory leaks**.  
Its signature: `drop_float_array(ar1: Array, ar2: Array)`  

The Array structs you pass to `drop_float_array` must be those you receive from the FFI function. For examples, see the `Array` struct and tests in [ffi.rs](src/ffi.rs), and the `_FFIArray` class in [convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).

### Building the Shared Library
Running `cargo build --release` will build an artefact called `liblonlat_bng.dylib` on OSX, and `liblonlat_bng.a` on `*nix` systems. Note that you'll have to generate `liblonlat_bng.so` for `*nix` hosts using the following steps:

- `ar -x target/release/liblonlat_bng.a`
- `gcc -shared *.o -o target/release/liblonlat_bng.so -lrt` 

## As a Python Package
`convert_bng` is [available](https://pypi.python.org/pypi/convertbng/) from PyPI for OSX and *nix:  
`pip install convertbng`  
More information is available in its [repository](https://github.com/urschrei/rust_bng)

# Benchmark
An IPython (sorry, *Jupyter*) notebook with some benchmarks is [here](rust_BNG.ipynb)

## Installation Instructions
- Ensure you have Rust 1.x installed
- If you want to run the benchmark yourself, install Pandas, IPython, Numpy, Pyproj, and their dependencies
- Clone this repository
- run `cargo build --release` from the repo root
- run `ipython notebook`, and open `rust_BNG`.

## Results
Test machine:  
- Late-2012 27" iMac
- 3.4 GHz Intel Core i7
- 4 cores (8 logical CPUs)
- 16GB 1600 MHz DDR3 RAM  

| Coordinates    | Method | Time (ms) | Speedup |
|:---------------|:------:|:---------:|--------:|
| 100k (50 runs) | Python | 547       | N/A     |
|                |**Rust**| 68.2      |**8x**   |
|                | Pyproj | 57.2      | 9.5x    |
| 1mm (50 runs)  | Python | 5560      | N/A     |
|                |**Rust**| 624       |**9.1x** |
|                | Pyproj | 537       | 10.9x   |
| 10mm (50 runs) | Python | 54500     | N/A     |
|                |**Rust**| 6150      |**8.9x** |
|                | Pyproj | 4710      | 11.5    | 


## Conclusion
Using multithreading gives excellent performance (Pyproj – which is a compiled [Cython](http://cython.org) binary – is now only ~20% faster than Rust, on average). Not bad, considering the relative youth of Rust *as a language* (let alone this library), and the maturity of the [PROJ.4](https://en.wikipedia.org/wiki/PROJ.4) project.

# Comparing Crossbeam and Rayon
Comparing how varying threads and weights affects overall speed, using [`cargo bench`](benches/benchmarks.rs)  
On both 2-core i5 and 8-core i7 machines, running `convert_bng_threaded_vec` using one thread per core gives optimum performance, whereas Rayon does a good job at choosing its own optimum weight.

<img src="crossbeam_v_rayon.png" alt="Comparison" style="width: 789px;"/>

# License
[MIT](license.txt)  

This software makes use of OSTN02 data, which is © Crown copyright, Ordnance Survey and the Ministry of Defence (MOD) 2002. All rights reserved. Provided under the BSD 2-clause [license](OSTN02_license.txt).

† Really, pyproj?  
[![mediocre](mediocre.png)]( "MEDIOCRE")
