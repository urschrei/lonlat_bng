# -*- coding: utf-8 -*-

from ctypes import cdll, c_float, Structure, POINTER, c_uint32, c_size_t, c_void_p, cast
from sys import platform
# from bng import bng
import numpy as np
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
    res_list = [(i.a, i.b) for i in iter(res)]
    drop_bng_array(array)
    return res_list


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
    res_list = [(i.a, i.b) for i in iter(res)]
    drop_ll_array(array)
    return res_list

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
drop_bng_array = lib.drop_int_array
drop_bng_array.argtypes = (BNG_FFIArray,)
drop_bng_array.restype = None
drop_ll_array = lib.drop_float_array
drop_ll_array.argtypes = (LONLAT_FFIArray,)
drop_ll_array.restype = None 


def convertbng_threaded(lons, lats):
    """ Multi-threaded lon lat to BNG wrapper """
    return convert_bng(lons, lats)

def convertlonlat_threaded(eastings, northings):
    """ Multi-threaded BNG to lon, lat wrapper """
    return convert_lonlat(eastings, northings)

# UK bounding box
N = 55.811741
E = 1.768960
S = 49.871159
W = -6.379880

pp_bng = pyproj.Proj(init='epsg:27700')
wgs84 = pyproj.Proj(init='epsg:4326')

num_coords = 1000000
lon_ls = list(np.random.uniform(W, E, [num_coords]))
lat_ls = list(np.random.uniform(S, N, [num_coords]))

# actually test the thing
# print("Pyproj")
# proj_res = zip(*pyproj.transform(wgs84, pp_bng, lon_ls, lat_ls))
print("Threaded lon, lat --> BNG")
convertbng_threaded([-0.32824866], [51.44533267])
print("Threaded BNG --> lon, lat")
convertlonlat_threaded([516276], [173141])
print("Threaded conversion of %s lon, lat coords --> BNG…" % num_coords)
convertbng_threaded(lon_ls, lat_ls)
print("done.")
