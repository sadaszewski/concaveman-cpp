from concaveman import concaveman2d
import numpy as np

def twospirals(n_points, noise=.5):
    """
     Returns the two spirals dataset.
    """
    n = np.sqrt(np.random.rand(n_points,1)) * 780 * (2*np.pi)/360
    d1x = -np.cos(n)*n + np.random.rand(n_points,1) * noise
    d1y = np.sin(n)*n + np.random.rand(n_points,1) * noise
    return (np.vstack((np.hstack((d1x,d1y)),np.hstack((-d1x,-d1y)))),
            np.hstack((np.zeros(n_points),np.ones(n_points))))

X, y = twospirals(1000)

import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import json
import numpy as np
from scipy.spatial import ConvexHull

with open('../../../data/points-1k.json', 'r') as f:
    pts = json.load(f)
    pts = np.array(pts)

h = ConvexHull(pts)
cc = concaveman2d(pts, h.vertices, 2, 0.005)

plt.plot(pts[:,0], pts[:,1], 'b*',
    cc[:,0].tolist() + [ cc[0,0] ],
    cc[:,1].tolist() + [ cc[0,1] ], 'r-')
plt.show()
