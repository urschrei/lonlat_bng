[![Build Status](https://travis-ci.org/urschrei/lonlat_bng.png?branch=master)](https://travis-ci.org/urschrei/lonlat_bng) [![](https://img.shields.io/crates/v/lonlat_bng.svg)](https://crates.io/crates/lonlat_bng) [![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](license.txt)  

## Introduction
An attempt at speeding up the conversion between decimal longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates, using an external Rust binary and Python FFI.

## Motivation
Python is relatively slow; this type of conversion is usually carried out in bulk, so an order-of-magnitude improvement could save precious minutes

## Library Use
### As a Rust Library
Add the following to your `Cargo.toml` (the latest version is displayed on the second badge at the top of this screen)  

    lonlat_bng = "x.x.x"

Full library documentation is available [here](http://urschrei.github.io/lonlat_bng/)  

**Note that `lon`, `lat` coordinates outside the UK bounding box will be transformed to `(9999, 9999)`, which cannot be mapped.**  

The native functions exposed by the library are:

`lonlat_bng::convert_bng(&f32, &f32) -> (i32, i32)`  
`lonlat_bng::ostn02::convert_osgb36(&f64, &f64) -> (f64, f64)`  
`lonlat_bng::convert_lonlat(&i21 &i32) -> (f32, f32)`  

`lonlat_bng::convert_to_bng_threaded_vec(Vec<(&f32, &f32)>) -> Vec<(i32, i32)>`  
`lonlat_bng::convert_to_osgb36_threaded_vec(Vec<(&f64, &f64)>) -> Vec<(f64, f64)>`  
`lonlat_bng::convert_to_lonlat_threaded_vec(Vec<(&i32, &i32)>) -> Vec<(f32, f32)>`  

### As an FFI Library
The FFI C-compatible functions exposed by the library are:  
`convert_to_bng_threaded(Array, Array) -> Array`  
`convert_to_osgb36_threaded(Array, Array) -> Array`  
`convert_to_lonlat_threaded(Array, Array) -> Array`  

And for freeing the memory allocated by the above  
`drop_int_array(Array) -> Null`  
`drop_float_array(Array) -> Null` (for `convert_to_lonlat_threaded` and `convert_to_osgb36_threaded`)  

The `Array`s must contain 32-bit `Float`s and 32-bit `Int`s, respectively. `Arrays` used in `convert_osgb36_threaded` must contain 64-bit `Float`s. For examples, see the `Array` struct and tests in [lib.rs](src/lib.rs), and the `_BNG_FFIArray` class in [convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).  

#### FFI and Memory Management
If your FFI library implements `convert_to_bng_threaded`, it **must** also implement `drop_int_array`, and if it implements `convert_to_lonlat_threaded` or `convert_to_osgb36_threaded`, it **must** implement `drop_float_array`. **Failing to do so will result in memory leaks**. 

#### Building the Shared Library
Running `cargo build --release` will build an artefact called `liblonlat_bng.dylib` on OSX, and `liblonlat_bng.a` on `*nix` systems. Note that you'll have to generate `liblonlat_bng.so` for `*nix` hosts using the following steps:

- `ar -x target/release/liblonlat_bng.a`
- `gcc -shared *.o -o target/release/liblonlat_bng.so -lrt` 

### As a Python Package
`convert_bng` is [available](https://pypi.python.org/pypi/convertbng/) from PyPI for OSX and *nix:  
`pip install convertbng`  
More information is available in its [repository](https://github.com/urschrei/rust_bng)

## Benchmark
An IPython (sorry, *Jupyter*) notebook with some benchmarks is [here](rust_BNG.ipynb)

### Installation Instructions
- Ensure you have Rust 1.x installed
- If you want to run the benchmark yourself, install Pandas, IPython, Numpy, Pyproj, and their dependencies
- Clone this repository
- run `cargo build --release` from the repo root
- run `ipython notebook`, and open `rust_BNG`.

### Results
Test machine:  
- Late-2012 27" iMac
- 3.4 GHz Intel Core i7
- 4 cores (8 logical CPUs)
- 16GB 1600 MHz DDR3 RAM  

| Coordinates    | Method | Time (ms) | Speedup |
|:---------------|:------:|:---------:|--------:|
| 100k (50 runs) | Python | 547       | N/A     |
|                |**Rust**| 73.9      |**7.4x** |
|                | Pyproj | 57.2      | 9.5x    |
| 1mm (50 runs)  | Python | 5560      | N/A     |
|                |**Rust**| 624       |**8.9x** |
|                | Pyproj | 510       | 10.9x   |
| 10mm (50 runs) | Python | 54500     | N/A     |
|                |**Rust**| 6360      |**8.6x** |
|                | Pyproj | 4710      | 11.5    | 


### Conclusion
Using multithreading gives excellent performance (Pyproj – which is a compiled [Cython](http://cython.org) binary – is now only ~20% faster than Rust, on average). Not bad, considering the relative youth of Rust *as a language* (let alone this library), and the maturity of the [PROJ.4](https://en.wikipedia.org/wiki/PROJ.4) project.

## Accuracy
The Helmert transform used in `convert_bng` and its threaded and vectorised versions is accurate to within 7 metres on average, and is **not suitable** for calculations used in e.g. surveying.  

**If higher accuracy is required, please use `convert_osgb36` and `convert_osgb36_threaded`, which adjust for local variation within the Terrestrial Reference Frame by incorporating OSTN02 data**. [See here](http://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) for more information.

## License
[MIT](license.txt)  

This software makes use of OSTN02 data, which is © Crown copyright, Ordnance Survey and the Ministry of Defence (MOD) 2002. All rights reserved. Provided under the BSD 2-clause [license](OSTN02_license.txt).

† Really, pyproj?  
[![mediocre](mediocre.png)]( "MEDIOCRE")
