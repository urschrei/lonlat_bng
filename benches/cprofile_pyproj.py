import numpy as np
from pyproj import Proj, transform

# London bounding box
N = 51.691874116909894
E = 0.3340155643740321
S = 51.28676016315085
W = -0.5103750689005356

# osgb36 = Proj(init='epsg:27700')
osgb36 = Proj('+init=EPSG:27700 +nadgrids=OSTN15_NTv2.gsb')
wgs84 = Proj(init='epsg:4326')
num_coords = 1000000
lon_ls = np.random.uniform(W, E, [num_coords])
lat_ls = np.random.uniform(S, N, [num_coords])

if __name__ == "__main__":
    for x in xrange(50):
        transform(wgs84, osgb36, lon_ls, lat_ls)
