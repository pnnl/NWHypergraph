import numpy as np
import nwhy

col = np.array([0, 3, 1, 0, 3, 3, 3])
row = np.array([0, 3, 1, 2, 1, 1, 1])
data = np.array([4, 5, 7, 9, 2, 2, 2])
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
g = nwhy.NWHypergraph(row, col, data, collapse=True)
#print(g)
# 
print("=====s_connected_component=====")
s1linegraph = g.s_linegraph(s=1, edges=True)
ccs = g.s_connected_components(s1linegraph, return_singleton=True)
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

# TODO s_centrality