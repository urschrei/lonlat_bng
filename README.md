[![Build Status](https://travis-ci.org/urschrei/lonlat_bng.png?branch=master)](https://travis-ci.org/urschrei/lonlat_bng)  
## Introduction
An attempt at speeding up the conversion between decimal longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates, using an external Rust binary and Python FFI.

## Motivation
Python is relatively slow; this type of conversion is usually carried out in bulk, so an order-of-magnitude improvement could save precious minutes

## Library Use
### As a Rust Library
Add the following to your `Cargo.toml` (You'll have to look up the latest version on [crates.io](https://crates.io/crates/lonlat_bng/))  

    lonlat_bng = "0.1.9"

Full library documentation is available [here](http://urschrei.github.io/lonlat_bng/)  

The native functions exposed by the library are:

`lonlat_bng::convert_bng(&f32, &f32) -> (i32, i32)`  
`lonlat_bng::convert_lonlat(&i21 &i32) -> (f32, f32)`  

`lonlat_bng::convert_to_bng_threaded_vec(Vec<(&f32, &f32)>) -> Vec<(i32, i32)>`  
`lonlat_bng::convert_to_lonlat_threaded_vec(Vec<(&i32, &i32)>) -> Vec<(f32, f32)>`  

### As an FFI Library
The FFI C-compatible functions exposed by the library are:  
`convert_to_bng_threaded(Array, Array) -> Array`  
`convert_to_lonlat_threaded(Array, Array) -> Array`  

And for freeing the memory allocated by the above  
`drop_int_array(Array) -> Null`  
`drop_float_array(Array) -> Null`  

The `Array`s must contain 32-bit `Float`s and 32-bit `Int`s, respectively. For examples, see the `Array` struct and tests in [lib.rs](src/lib.rs), and the `_BNG_FFIArray` class in [convertbng](https://github.com/urschrei/convertbng/blob/master/convertbng/util.py).  

#### FFI and Memory Management
If your FFI library implements `convert_to_bng_threaded`, it **must** also implement `drop_int_array`, and if it implements `convert_to_lonlat_threaded`, it **must** implement `drop_float_array`. **Failing to do so will result in memory leaks**. 

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

| Coordinates | Method | Time (ms) | Speedup |
|:------------|:------:|:---------:|--------:|
| 100k        | Python | 539       | N/A     |
|             |**Rust**| 104       |**5.2x** |
|             | Pyproj | 61.1      | 8.8x    |
| 1mm         | Python | 5430      | N/A     |
|             |**Rust**| 1070      |**5.1x** |
|             | Pyproj | 479       | 11.3    |
| 10mm        | Python | 54500     | N/A     |
|             |**Rust**| 10800     |**5.1x** |
|             | Pyproj | 4710      | 11.5    | 


### Conclusion
Using multithreading gives excellent performance (Pyproj – which is a compiled [Cython](http://cython.org) binary – is now only around 130% faster, on average). Not bad, considering the relative youth of Rust *as a language* (let alone this library), and the maturity of the [PROJ.4](https://en.wikipedia.org/wiki/PROJ.4) project.

## Accuracy
The Helmert transform used is accurate to within 7 metres on average, so this library is **not suitable** for calculations used in e.g. surveying. If higher accuracy is required, please use a product which incorporates the OSTN02 calculations, which adjust for local variation within the Terrestrial Reference Frame. [See here](http://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) for more information.

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
