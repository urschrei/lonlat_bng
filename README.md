## Introduction

An attempt at speeding up the conversion between decimal longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates, using an external Rust binary and Python FFI.

## Motivation

Python is relatively slow; this type of conversion is usually carried out in bulk, so an order-of-magnitude improvement could save precious minutes

## Instructions

- Ensure you have Rust 1.x installed
- If you want to run the benchmark youself, install Pandas, IPython, Numpy, Pyproj, and their dependencies
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
Convert 100,000 sets of random coordinates  

Python: 1 loops, best of 10: **804 ms** per loop  
Rust: 1 loops, best of 10: **204 ms** per loop  
Pyproj: 10 loops, best of 10: **99.5 ms** per loop 💅  
Rust (threaded): 10 loops, best of 10: **162.5 ms** per loop  


## Conclusion
Using multithreading, we can get much closer (pyproj is now only 65% faster). Not bad, considering the relative youth of Rust *as a language* (let alone this library), and the maturity of the [PROJ.4](https://en.wikipedia.org/wiki/PROJ.4) project.

## Package
`convert_bng` is [available](https://pypi.python.org/pypi/convertbng/) from PyPI:  
`pip install convertbng`  
Usage:

    from convertbng.util import convertbng
    convertbng(lon, lat)

## Benchmark machine spec:

- Mid-2011 Macbook Air
- 1.8 GHz Intel Core i7
- OSX 10.10
- Rust 1.0 (installed using Homebrew)
- Python 2.7.9
- Numpy 1.9.2
- Pandas 0.16.2
- Pyproj 1.9.4

## TODO

- Write a better real-world test
- Allow the converter to accept Numpy arrays

## License
MIT

† Really, pyproj?  
[![mediocre](mediocre.png)]( "MEDIOCRE")
