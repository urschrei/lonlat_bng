#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Standalone benchmark runner
"""

import cProfile
import pstats
import profile
import numpy as np

print("Running Rust (Ctypes + Cython) and Pyproj benchmarks, 1mm conversions x 50 runs\n")

# calibrate
pr = profile.Profile()
calibration = np.mean([pr.calibrate(100000) for x in xrange(5)])
# add the bias
profile.Profile.bias = calibration

cProfile.run(open('benches/cprofile_rust.py', 'rb'), 'benches/output_stats_rust')
rust = pstats.Stats('benches/output_stats_rust')

cProfile.run(open('benches/cprofile_rust_cython.py', 'rb'), 'benches/output_stats_rust_cython')
rust_cython = pstats.Stats('benches/output_stats_rust_cython')

cProfile.run(open('benches/cprofile_pyproj.py', 'rb'), 'benches/output_stats_pyproj')
pyproj_ = pstats.Stats('benches/output_stats_pyproj')

print("Rust Benchmark\n")
rust.sort_stats('cumulative').print_stats(5)
print("Rust Cython Benchmark\n")
rust_cython.sort_stats('cumulative').print_stats(5)
print("Pyproj Benchmark\n")
pyproj_.sort_stats('cumulative').print_stats(5)
