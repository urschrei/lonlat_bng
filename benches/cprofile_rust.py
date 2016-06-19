import numpy as np
from convertbng.util import convert_bng

# London bounding box
N = 51.691874116909894
E = 0.3340155643740321
S = 51.28676016315085
W = -0.5103750689005356

num_coords = 1000000
lon_ls = np.random.uniform(W, E, [num_coords])
lat_ls = np.random.uniform(S, N, [num_coords])

if __name__ == "__main__":
    for x in xrange(50):
        convert_bng(lon_ls, lat_ls)
