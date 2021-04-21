import nwhy
import numpy as np

row = np.array([0, 0, 0, 1, 1, 1, 2])
col = np.array([ 0,  1,  2,  0,  1,  2,  0])
data = np.array([1, 1, 1, 1, 1, 1, 1])

print('-- without collapsing')
h = nwhy.NWHypergraph(row, col, data)
print(h)

print('-- collapsing edges')
newe, equal_class = h.collapse_edges()
print(newe)
print(equal_class)

print('-- collapsing nodes')
newn, equal_class = h.collapse_nodes()
print(newn)
print(equal_class)

print('-- collapsing nodes and edges')
new, equal_class = h.collapse_nodes_and_edges()
print(new)
print(equal_class)