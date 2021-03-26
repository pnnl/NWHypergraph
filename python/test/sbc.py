import nwhy
import numpy as np
srow = np.array([4592, 4601, 4601, 4605, 4605, 4609, 4611, 4612, 4612, 4613, 4613,
       4615, 4615, 4616, 4616, 4617])
scol = np.array([4611, 4617, 4618, 4617, 4618, 4618, 4617, 4617, 4618, 4616, 4617,
       4617, 4618, 4617, 4618, 4618])
sdata = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

linegraph = nwhy.Slinegraph(srow, scol, sdata, 15, True)
sbc = linegraph.s_betweenness_centrality(normalized=False)
assert sbc[0] == 0