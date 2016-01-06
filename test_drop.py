import numpy as np
import pandas as pd
import math
from ctypes import cdll, c_float, Structure, ARRAY, POINTER, c_int32, c_uint32, c_size_t, c_void_p, cast
from sys import platform
from bng import bng
import pyproj


if platform == "darwin":
    ext = "dylib"
else:
    ext = "so"
    
lib = cdll.LoadLibrary('target/release/liblonlat_bng.' + ext)


class BNG_FFITuple(Structure):
    _fields_ = [("a", c_uint32),
                ("b", c_uint32)]

class BNG_FFIArray(Structure):
    _fields_ = [("data", c_void_p),
                ("len", c_size_t)]

    # Allow implicit conversions from a sequence of 32-bit unsigned
    # integers.
    @classmethod
    def from_param(cls, seq):
        return seq if isinstance(seq, cls) else cls(seq)

    # Wrap sequence of values. You can specify another type besides a
    # 32-bit unsigned integer.
    def __init__(self, seq, data_type = c_float):
        array_type = data_type * len(seq)
        raw_seq = array_type(*seq)
        self.data = cast(raw_seq, c_void_p)
        self.len = len(seq)

# A conversion function that cleans up the result value to make it
# nicer to consume.
def bng_void_array_to_tuple_list(array, _func, _args):
    res = cast(array.data, POINTER(BNG_FFITuple * array.len))[0]
    drop_array(array)
    return res


class LONLAT_FFITuple(Structure):
    _fields_ = [("a", c_float),
                ("b", c_float)]

class LONLAT_FFIArray(Structure):
    _fields_ = [("data", c_void_p),
                ("len", c_size_t)]

    # Allow implicit conversions from a sequence of 32-bit unsigned
    # integers.
    @classmethod
    def from_param(cls, seq):
        return seq if isinstance(seq, cls) else cls(seq)

    # Wrap sequence of values. You can specify another type besides a
    # 32-bit unsigned integer.
    def __init__(self, seq, data_type = c_uint32):
        array_type = data_type * len(seq)
        raw_seq = array_type(*seq)
        self.data = cast(raw_seq, c_void_p)
        self.len = len(seq)

# A conversion function that cleans up the result value to make it
# nicer to consume.
def lonlat_void_array_to_tuple_list(array, _func, _args):
    res = cast(array.data, POINTER(LONLAT_FFITuple * array.len))[0]
    return res

# Multi-threaded
convert_bng = lib.convert_to_bng
convert_bng.argtypes = (BNG_FFIArray, BNG_FFIArray)
convert_bng.restype = BNG_FFIArray
convert_bng.errcheck = bng_void_array_to_tuple_list

convert_lonlat = lib.convert_to_lonlat
convert_lonlat.argtypes = (LONLAT_FFIArray, LONLAT_FFIArray)
convert_lonlat.restype = LONLAT_FFIArray
convert_lonlat.errcheck = lonlat_void_array_to_tuple_list

# cleanup
drop_array = lib.drop_array 
drop_array.argtypes = (BNG_FFIArray,)


def convertbng_threaded(lons, lats):
    """ Multi-threaded lon lat to BNG wrapper """
    return [(i.a, i.b) for i in iter(convert_bng(lons, lats))]

def convertlonlat_threaded(eastings, northings):
    """ Multi-threaded BNG to lon, lat wrapper """
    return [(i.a, i.b) for i in iter(convert_lonlat(eastings, northings))]

# actually test the thing
print convertbng_threaded([-0.32824866], [51.44533267])