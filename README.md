[![Build Status](https://travis-ci.org/urschrei/lonlat_bng.png?branch=master)](https://travis-ci.org/urschrei/lonlat_bng)  
## Introduction

An attempt at speeding up the conversion between decimal longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates, using an external Rust binary and Python FFI.

## Motivation

Python is relatively slow; this type of conversion is usually carried out in bulk, so an order-of-magnitude improvement could save precious minutes

## Instructions

- Ensure you have Rust 1.x installed
- If you want to run the benchmark yourself, install Pandas, IPython, Numpy, Pyproj, and their dependencies
- Clone this repository
- run `cargo build --release` from the repo root
- run `ipython notebook`, and open `rust_BNG`.

## Benchmark
An IPython (sorry, *Jupyter*) notebook with some benchmarks is [here](rust_BNG.ipynb)

## Results
### Simple Test
Python: 10000 loops, best of 10: **31 µs** per loop  
Rust: 100000 loops, best of 10: **2.04 µs** per loop* 💅  
Pyproj: 100000 loops, best of 10: **11.8 µs** per loop<sup>†</sup>  
<sup>*</sup>Test warns that intermediate results may have been cached  

### Real-world Test
Convert 10,000000 sets of random coordinates  

Python: 1 loops, best of 10: **804 ms** per loop  
Rust: 1 loops, best of 10: **204 ms** per loop  
Pyproj: 10 loops, best of 10: **99.5 ms** per loop 💅  
Rust (threaded): 10 loops, best of 10: **162.5 ms** per loop  


## Conclusion
Using multithreading, we can get much closer (pyproj is now only 65% faster). Not bad, considering the relative youth of Rust *as a language* (let alone this library), and the maturity of the [PROJ.4](https://en.wikipedia.org/wiki/PROJ.4) project.

## Accuracy
The Helmert transform used is accurate to within 4 – 5 metres, so this library is **not suitable** for calculations used in e.g. surveying. If higher accuracy is required, please use a product which incorporates the OSTN02 calculations, which adjust for local variation within the Terrestrial Reference Frame. [See here](http://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) for more information.

## Library Use
### As a Rust Library
Add the following to your `Cargo.toml` (You'll have to look up the latest version on [crates.io](https://crates.io/crates/lonlat_bng/))  

    lonlat_bng = "0.1.8"

Full library documentation is available [here](http://urschrei.github.io/lonlat_bng/)  

The native functions exposed by the library are:

`lonlat_bng::convert_bng(&f32, &f32) -> (i32, i32)`  
`lonlat_bng::convert_lonlat(&i21 &i32) -> (f32, f32)`  

`lonlat_bng::convert_to_bng_threaded_vec(Vec<(&f32, &f32)>) -> Vec<(i32, i32)>`  
`lonlat_bng::convert_to_lonlat_threaded_vec(Vec<(&i32, &i32)>) -> Vec<(f32, f32)>`  

### As an FFI Library
The FFI C-compatible functions exposed by the library are:

`convert_to_bng_threaded(Array, Array) -> Array`  
`convert_to_bng_threaded(Array, Array) -> Array`  

The `Array`s must contain 32-bit `Float`s and 32-bit `Int`s, respectively. For examples, see the `Array` struct and tests in [lib.rs](src/lib.rs), and the `_BNG_FFIArray` class in [convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py)

### As a Python Package
`convert_bng` is [available](https://pypi.python.org/pypi/convertbng/) from PyPI for OSX and *nix:  
`pip install convertbng`  
More information is available in its [repository](https://github.com/urschrei/rust_bng)

## Benchmark machine spec:

- Mid-2011 Macbook Air
- 1.8 GHz Intel Core i7
- OSX 10.10
- Rust 1.0 (installed using Homebrew)
- Python 2.7.9
- Numpy 1.9.2
- Pandas 0.16.2
- Pyproj 1.9.4

## License
MIT

† Really, pyproj?  
[![mediocre](mediocre.png)]( "MEDIOCRE")
