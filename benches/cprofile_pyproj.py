import numpy as np
import pyproj

# UK bounding box
N = 55.811741
E = 1.768960
S = 49.871159
W = -6.379880

bng = pyproj.Proj(init='epsg:27700')
wgs84 = pyproj.Proj(init='epsg:4326')

num_coords = 1000000
lon_ls = np.random.uniform(W, E, [num_coords])
lat_ls = np.random.uniform(S, N, [num_coords])

if __name__ == "__main__":
    for x in xrange(50):
        pyproj.transform(wgs84, bng, lon_ls, lat_ls)
