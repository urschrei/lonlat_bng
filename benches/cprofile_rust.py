import numpy as np
from ctypes import cdll, c_double, Structure, POINTER, c_size_t, c_void_p, cast
from sys import platform
import pyproj
from array import array

if platform == "darwin":
    ext = "dylib"
else:
    ext = "so"
    
lib = cdll.LoadLibrary('target/release/liblonlat_bng.' + ext)

class _FFIArray(Structure):
    """ Convert sequence of floats to a C-compatible void array """
    _fields_ = [("data", c_void_p),
                ("len", c_size_t)]

    @classmethod
    def from_param(cls, seq):
        """  Allow implicit conversions from a sequence of 64-bit floats."""
        return seq if isinstance(seq, cls) else cls(seq)

    def __init__(self, seq, data_type = c_double):
        """
        Convert sequence of values into array, then ctypes Structure

        Rather than checking types (bad), we just try to blam seq
        into a ctypes object using from_buffer. If that doesn't work,
        we try successively more conservative approaches:
        numpy array -> array.array -> read-only buffer -> CPython iterable
        """
        if isinstance(seq, float):
            seq = array('d', [seq])
        try:
            len(seq)
        except TypeError:
             # we've got an iterator or a generator, so consume it
            seq = array('d', seq)
        array_type = data_type * len(seq)
        try:
            raw_seq = array_type.from_buffer(seq.astype(np.float64))
        except (TypeError, AttributeError):
            try:
                raw_seq = array_type.from_buffer_copy(seq.astype(np.float64))
            except (TypeError, AttributeError):
                # it's a list or a tuple
                raw_seq = array_type.from_buffer(array('d', seq))
        self.data = cast(raw_seq, c_void_p)
        self.len = len(seq)
        

class _Result_Tuple(Structure):
    """ Container for returned FFI data """
    _fields_ = [("e", _FFIArray),
                ("n", _FFIArray)]


def _void_array_to_list(restuple, _func, _args):
    """ Convert the FFI result to Python data structures """
    # eastings = POINTER(c_double * restuple.e.len).from_buffer_copy(restuple.e)[0]
    # northings = POINTER(c_double * restuple.n.len).from_buffer_copy(restuple.n)[0]
    res_list = [
        list(POINTER(c_double * restuple.e.len).from_buffer_copy(restuple.e)[0]),
        list(POINTER(c_double * restuple.n.len).from_buffer_copy(restuple.n)[0])
    ]
    drop_array(restuple.e, restuple.n)
    return res_list

# Multi-threaded FFI functions
convert_bng = lib.convert_to_bng_threaded
convert_bng.argtypes = (_FFIArray, _FFIArray)
convert_bng.restype = _Result_Tuple
convert_bng.errcheck = _void_array_to_list
convert_bng.__doc__ = """
    Multi-threaded lon, lat --> BNG conversion
    Returns a list of two lists containing Easting and Northing floats,
    respectively
    Uses the Helmert transform
    """

# Free FFI-allocated memory
drop_array = lib.drop_float_array
drop_array.argtypes = (_FFIArray, _FFIArray)
drop_array.restype = None

# UK bounding box
N = 55.811741
E = 1.768960
S = 49.871159
W = -6.379880

bng = pyproj.Proj(init='epsg:27700')
wgs84 = pyproj.Proj(init='epsg:4326')

num_coords = 1000000
lon_ls = list(np.random.uniform(W, E, [num_coords]))
lat_ls = list(np.random.uniform(S, N, [num_coords]))

if __name__ == "__main__":
    for x in xrange(50):
        convert_bng(lon_ls, lat_ls)