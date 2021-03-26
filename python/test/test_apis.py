import numpy as np
import nwhy

col = np.array([0, 0, 0, 1, 1, 1])
row = np.array([0, 1, 2, 0, 1, 2])
weight = np.array([1, 1, 1, 1, 1, 1])
#create a hypergraph hg 
hg = nwhy.NWHypergraph(row, col, weight, collapse=False)
#compute the s-line graph of hg with s=2
s2lg = hg.s_linegraph(s=2, edges=True)
#query whether the 2-line graph is 2-connected
tmp = s2lg.is_s_connected()
#query the neighboring hyperedges of hyperedge 0
sn = s2lg.s_neighbors(v=0)
#query the s-degree of hyperedge 0
sd = s2lg.s_degree(v=0)
#compute the s-connected components of line graph
scc = s2lg.s_connected_components()
#compute the s-distance between hyperedge 0 and 1
sdist = s2lg.s_distance(src=0, dest=1)
#compute the s-path between hyperedge 0 and 1
sp = s2lg.s_path(src=0, dest=1)
#compute the s-betweenness centrality of line graph
sbc = s2lg.s_betweenness_centrality(normalized=True)
#compute the s-closeness centrality of every hyperedge 
#in line graph
sc = s2lg.s_closeness_centrality(v=None)
#compute the s-harmonic closeness centrality of every 
#hyperedge in line graph
shc = s2lg.s_harmonic_closeness_centrality(v=None)
#compute s-eccentricity of every hyperedges in line graph
se = s2lg.s_eccentricity(v=None)