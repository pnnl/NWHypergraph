import numpy as np
import nwhy

col = np.array([0, 3, 1, 0, 3])
row = np.array([0, 3, 1, 2, 1])
data = np.array([4, 5, 7, 9, 2])
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
g = nwhy.NWHypergraph(row, col, data)
#print(g)
# 
print("=====s_connected_component=====")
s1linegraph = g.s_linegraph(s=1, edges=True)
ccs = g.s_connected_component(s1linegraph, return_singleton=True)
ccs0 = s1linegraph.s_connected_component(return_singleton=True)
#s1linegraph.s_distance(s=1, source=0, edges=True)
print("ccs0:", ccs0)

# return python::list of python::set
ccs1 = g.s_connected_component(s=1, edges=True, return_singleton=True)
print("ccs1:", ccs1)
# [{0}, {1,2}, {3}]
ccs2 = g.s_connected_component(s=2, edges=True, return_singleton=True)
print("ccs2:", ccs2)
# [{1,2}]

# s_distance
print("=====s_distance=====")
source = 1
destination = 2
dist0 = g.s_distance(s1linegraph, src=source, dest=destination)

dist1 = g.s_distance(src=0, dest=2, s=1, edges=True)
print("dist from", 0, "to", 2, "is", dist1)
dist2 = s1linegraph.s_distance(src=source, dest=destination)
print("dist from", source, "to", destination, "is", dist2)
#dis1 = infy or a value

# s_neighbo
print("=====s_neighbor=====")
vertex=0
neighborsofv0 = g.s_neighbor(s1linegraph, v=vertex)
neighborsofv1 = g.s_neighbor(v=1, s=1, edges=True)
print("neighbors of vertex", 1, "are:", neighborsofv1)
neighborsofv2 = s1linegraph.s_neighbor(v=vertex)
print("neighbors of vertex", vertex, "are:", neighborsofv2)

# TODO s_centrality