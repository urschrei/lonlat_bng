import numpy as np
from convertbng.cutil import convert_bng

# UK bounding box
N = 55.811741
E = 1.768960
S = 49.871159
W = -6.379880

num_coords = 1000000
lon_ls = np.random.uniform(W, E, [num_coords])
lat_ls = np.random.uniform(S, N, [num_coords])

if __name__ == "__main__":
    for x in xrange(50):
        convert_bng(lon_ls, lat_ls)
