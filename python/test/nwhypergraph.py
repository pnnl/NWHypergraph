import nwhy
import numpy as np

row = np.array([0, 0, 0, 1, 1, 1, 2])
col = np.array([ 0,  1,  2,  0,  1,  2,  0])
data = np.array([1, 1, 1, 1, 1, 1, 1])

print('-- without collapsing')
h = nwhy.NWHypergraph(row, col, data)
print(h)

print('-- collapsing edges without returning equal class')
equal_class = h.collapse_edges()
print(equal_class)

print('-- collapsing nodes without returning equal class')
equal_class = h.collapse_nodes()
print(equal_class)

print('-- collapsing nodes and edges without returning equal class')
equal_class = h.collapse_nodes_and_edges()
print(equal_class)

print('-- collapsing edges with returning equal class')
equal_class = h.collapse_edges(return_equivalence_class=True)
print(equal_class)

print('-- collapsing nodes with returning equal class')
equal_class = h.collapse_nodes(return_equivalence_class=True)
print(equal_class)

print('-- collapsing nodes and edges with returning equal class')
equal_class = h.collapse_nodes_and_edges(return_equivalence_class=True)
print(equal_class)