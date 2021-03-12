import numpy as np
import nwhy

col = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4,
         4, 4, 4, 4, 5, 5, 6])
row = np.array([ 0,  1,  2,  5,  9, 10,  0,  1,  2,  3,  4,  5,  0,  1,  3,  5,  6,
          7,  1,  2,  5,  0,  1,  2,  3,  5,  4,  8,  5])
data = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1])
'''
# test weighted line graph
newx, newy, newz, oldx, oldy, oldz = nwhy.convert_to_s_overlap(row, col, data, 1)
print("newx:", newx)
print("newy:", newy)
print("newz:", newz)
print("oldx:", oldx)
print("oldy:", oldy)
print("oldz:", oldz)

# test unweighted line graph
emptydata = np.array([])
newu, newv, news, oldu, oldv, olds = nwhy.convert_to_s_overlap(row, col, emptydata, 1)
print("newu:", newu)
print("newv:", newv)
print("news:", news)
print("oldu:", oldu)
print("oldv:", oldv)
print("olds:", olds)
'''
# declare NWHypergraph
h = nwhy.NWHypergraph(row, col, data, collapse=True)
#print(h)
# 
print("=====s_connected_component=====")
s1linegraph = h.s_linegraph(s=1, edges=True)
ccs0 = s1linegraph.s_connected_components(return_singleton=True)
#s1linegraph.s_distance(s=1, source=0, edges=True)
print("ccs0:", ccs0)

# s_distance
print("=====s_distance=====")
source = 1
destination = 2
dist2 = s1linegraph.s_distance(src=source, dest=destination)
print("dist from", source, "to", destination, "is", dist2)
#dis1 = infy or a value

# s_neighbor
print("=====s_neighbor=====")
vertex=0
neighborsofv2 = s1linegraph.s_neighbors(v=vertex)
print("neighbors of vertex", vertex, "are:", neighborsofv2)

# s_path


# s_betweenness_centrality
print("=====s_betweenness_centrality=====")
#s1edgelinegraph = h.s_linegraph(s=1, edges=False)
print("normalized:")
print(s1linegraph.s_betweenness_centrality(normalized=True))
print("unnormalized:")
s3edgelinegraph = h.s_linegraph(s=3, edges=True)
bc = s1linegraph.s_betweenness_centrality(normalized=False)
print(bc)
#for k, v in bc.items():
#    print(k, ' : ', v)

# s_closeness_centrality
print("=====s_closeness_centrality=====")
outofboundvertex = 11
print(s1linegraph.s_closeness_centrality(vertex))
print(s1linegraph.s_closeness_centrality(outofboundvertex))
print(s1linegraph.s_closeness_centrality())

# s_harmonic_closeness_centrality
print("=====s_harmonic_closeness_centrality=====")
print(s1linegraph.s_harmonic_closeness_centrality(vertex))
print(s1linegraph.s_harmonic_closeness_centrality(outofboundvertex))
print(s1linegraph.s_harmonic_closeness_centrality())

# s_eccentricity
print("=====s_eccentricity=====")
print(s1linegraph.s_eccentricity(2))
print(s1linegraph.s_eccentricity(outofboundvertex))
print(s1linegraph.s_eccentricity())

import gc
gc.collect()