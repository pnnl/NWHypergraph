import numpy as np
import nwhy

col = np.array([0, 3, 1, 0])
row = np.array([0, 3, 1, 2])
data = np.array([4, 5, 7, 9])

# test weighted line graph
x, y, z = nwhy.convert_to_s_overlap(row, col, data, 1)
print("row:", row)
print("x:", x)
print("y:", y)
print("z:", z)

# test unweighted line graph
emptydata = np.array([])
u, v, s = nwhy.convert_to_s_overlap(row, col, emptydata, 1)
print("u:", x)
print("v:", y)
print("s:", s)
