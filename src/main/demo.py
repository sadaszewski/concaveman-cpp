from concaveman import concaveman2d
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import json
import numpy as np
from scipy.spatial import ConvexHull


with open('../../data/points-1k.json', 'r') as f:
    pts = json.load(f)
    pts = np.array(pts)

h = ConvexHull(pts)
cc = concaveman2d(pts, h.vertices, 2, 0.005)

plt.plot(pts[:,0], pts[:,1], 'b*',
    cc[:,0].tolist() + [ cc[0,0] ],
    cc[:,1].tolist() + [ cc[0,1] ], 'r-')
plt.show()
