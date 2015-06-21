## Introduction

An attempt at speeding up the conversion between decimal longitude and latitude and British national grid coordinates, using an external Rust binary and Python FFI.

## Motivation

Python is relatively slow; this type of conversion is usually carried out in bulk, so a 10x improvement could save precious minutes

## Instructions

- Ensure you have Rust 1.0 installed
- If you want to run the benchmark youself, install Pandas, IPython, Numpy, Pyproj, and their dependencies
- Clone this repository
- run `cargo build --release` from the repo root
- run `ipython notebook`, and open `rust_BNG`.

## Benchmark
An IPython (sorry, *Jupyter*) notebook with some benchmarks is [here](rust_BNG.ipynb)

## Results:

# Simple Test
Python: 10000 loops, best of 10: **31 µs** per loop  
Rust: 100000 loops, best of 10: **2.04 µs** per loop*  
Pyproj: 100000 loops, best of 10: **11.8** µs per loop*  
<sup>*</sup>Test warns that intermediate results may have been cached  

An approximately 15x improvement on the simple test, a 5x improvement on Pyproj, and around a 2x improvement on the "real-world" test.

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

