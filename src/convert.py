from ctypes import cdll, c_float, Structure, ARRAY, c_int32
from sys import platform

if platform == "darwin":
    ext = "dylib"
else:
    ext = "so"

# hacky
# http://stackoverflow.com/a/30789980/416626
class Int32_2(Structure):
    _fields_ = [("array", ARRAY(c_int32, 2))]

lib = cdll.LoadLibrary('target/debug/liblonlat_bng.' + ext)
convert = lib.convert
convert.restype = Int32_2

lon = -0.32824866
lat = 51.44533267
print [i for i in convert(c_float(lon), c_float(lat)).array]
