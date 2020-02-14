# need to be in the etc conda env
from planar import Polygon
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib

fig = plt.figure()
ax = Axes3D(fig)

h1 = np.asarray(Polygon.regular(6, radius=5.5, angle=30))
h1 = np.hstack([h1, np.zeros(6).reshape([6,1])])
h2 = h1 
h2[:,2] =- 3

v1 = [list(zip(h1.tolist()))]
v2 = [list(zip(h2.tolist()))]
ax.add_collection3d(Poly3DCollection(v1))
ax.add_collection3d(Poly3DCollection(v2))

plt.show()
