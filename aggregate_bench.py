#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A wrapper for cargo bench
Its numeric output is parsed and dumped to a csv
Pass an an optional dependent variable from the command line
(C) Stephan Hügel 2016
License: MIT
"""
import os
import sys
import csv
from subprocess import check_output
import re
pattern = re.compile(
    r"bench:\s+([0-9,]*)\D+([0-9,]*)"
)

def dump_benchmark(pattern, filepath=None, headers=None, dep_var=None, **kwargs):
    """ If I have to append benchmark output to a CSV once more I'm going
    to drown the world in a bath of fire. This should just work.
    Customise with your own output path and header row.
    dep_var is an optional dependent variable.
    """
    if not filepath:
        filepath = "measurements.csv"
    if not headers:
        headers = ['time', 'error']
    # run cargo bench in cwd, capture output
    result = re.search(pattern, check_output(["cargo", "bench"]))
    output = [int(group.translate(None, ',')) for group in result.groups()]
    # this one's special because wtf are we measuring without a dependent variable
    if dep_var:
        headers.append("dependent_variable")
        output.append(dep_var)
    # anything else will get written as a CSV header row and value
    # nothing prevents you from writing rows that don't have a header
    for k, v in kwargs.items():
        headers.append(k),
        output.append(v)
    # check that path and file exist, or create them
    path_wrangle(filepath, headers)
    # get rid of nasty commas, convert to int, and write
    with open(filepath, 'a') as handle:
        wr = csv.writer(handle)
        wr.writerow(output)

def path_wrangle(filepath, headers):
    """ Check for or create path and output file
    There's no error handling, because noisy failure's probably a good thing
    """
    # check for or create directory path
    directory = os.path.split(filepath)[0]
    if not os.path.exists(directory):
            os.makedirs(directory)
    # if that worked, see if the CSV exists, and create a new one if necessary
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as newhandle:
            wr = csv.writer(newhandle)
            wr.writerow(headers)

if __name__ == "__main__":
    dep_var = None
    # So brittle. Shhh.
    if sys.argv[1] is not None:
        dep_var = sys.argv[1]
    dump_benchmark(pattern, filepath="benches/measurements.csv", dep_var=dep_var)
