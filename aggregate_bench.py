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

def dump_benchmark(pattern, filepath=None, headers=None, dep_var=None):
    """ If I have to append benchmark output to a CSV once more I'm going
    to drown the world in a bath of fire. This should just work.
    Customise with your own output path and header row.
    dep_var is an optional dependent variable.
    """
    if not filepath:
        filepath = "benches/measurements.csv"
    if not headers:
        headers = ['time', 'error']
    # run cargo bench in cwd, capture output
    result = re.search(pattern, check_output(["cargo", "bench"]))
    output = [int(group.translate(None, ',')) for group in result.groups()]
    if dep_var:
        headers.append("dependent_variable")
        output.append(dep_var)
    # TODO this should be way more robust
    # does the CSV exist?
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as newhandle:
            wr = csv.writer(newhandle)
            wr.writerow(headers)
    # get rid of nasty commas, convert to int, and write
    with open(filepath, 'a') as handle:
        wr = csv.writer(handle)
        wr.writerow(output)

if __name__ == "__main__":
    dep_var = None
    # So brittle. Shhh.
    if sys.argv[1] is not None:
        dep_var = sys.argv[1]
    dump_benchmark(pattern, dep_var=dep_var)
