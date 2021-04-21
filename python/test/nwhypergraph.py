import nwhy
import numpy as np
col = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4,
         4, 4, 4, 4, 5, 5, 6])
row = np.array([ 0,  1,  2,  5,  9, 10,  0,  1,  2,  3,  4,  5,  0,  1,  3,  5,  6,
          7,  1,  2,  5,  0,  1,  2,  3,  5,  4,  8,  5])
data = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1])
h = nwhy.NWHypergraph(row, col, data, collapse=False)


nh = h.collapse_edges(return_equavilance_class=False)
print(nh)