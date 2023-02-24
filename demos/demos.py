#!/usr/bin/env python

"""
demos / checks for concaveman
"""

import numpy as np

import matplotlib.pyplot as plt

from concaveman import ConcaveHull

def circular_points(N):
    """
    create a bunch of random points in a circle
    """
    # create a bunch of random points
    points_in_circle = []

    while len(points_in_circle) < N:
        n = max(100, N // 2)
#        print(f"picking {n} random points")
        points = (np.random.random((n, 2)) - 0.5) * 100

        # check which ones are in the circle
        points_in_circle.extend(points[np.hypot(points[:,0], points[:,1]) <= 50.0])
        # print(f"there are now {len(points_in_circle)} in the circle")
    return np.array(points_in_circle[:N])

# get some points
points = circular_points(100)

# compute the hulls
hull = ConcaveHull(points, concavity=1.0)
hull_1 = hull.cc_hull

hull.concavity = 3.0
hull.compute()
hull_2 = hull.cc_hull


hull.concavity = np.inf
hull.compute()
hull_inf = hull.cc_hull


fig, ax = plt.subplots()
ax.plot(points[:, 0], points[:, 1], '.')
ax.plot(hull_1[:, 0], hull_1[:, 1], label="alpha=1.0")
ax.plot(hull_2[:, 0], hull_2[:, 1], label="alpha=3.0")
ax.plot(hull_inf[:, 0], hull_inf[:, 1], label="alpha=Inf")
ax.legend()

ax.set_xlim(-80, 80)
ax.set_ylim(-80, 80)
ax.set_aspect('equal', 'box')

plt.show()

