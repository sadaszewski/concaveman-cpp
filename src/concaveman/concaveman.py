import numpy as np
import scipy.spatial

from .build_hull import concave_hull


class ConcaveHull():
    """
    class modeled after the scipy.spatial.ConvexHull implementation

    One of these days, we could add area and other computations.
    """

    def __init__(self, points, hull=None, concavity=2.0, length_threshold=0.0):

        self.points = points

        if hull is not None:
            if len(hull.shape) != 1:
                raise ValueError('hull must be a 1-D array')

            if np.any(hull >= len(points)) or np.any(hull < 0):
                raise ValueError('hull indices out of bounds')
            self.cvx_hull = hull
            self.cvx_hull = scipy.spatial.ConvexHull(points).vertices

        self.concavity = concavity
        self.length_threshold = length_threshold

        self.compute()

    @property
    def concavity(self):
        return self._concavity
    @concavity.setter
    def concavity(self, value):
        self._concavity = value
        self._cc_hull = None

    @property
    def length_threshold(self):
        return self._length_threshold
    @length_threshold.setter
    def length_threshold(self, value):
        self._length_threshold = value
        self._cc_hull = None

    @property
    def points(self):
        return self._points
    @points.setter
    def points(self, points):
        points = np.asarray(points, dtype=np.float64)
        if len(points.shape) != 2 or points.shape[1] != 2:
            raise ValueError('points must be an Nx2 array')
        self._points = points
        self.cvx_hull = None
        self._cc_hull = None

    @property
    def cc_hull(self):
        if self._cc_hull is None:
            self.compute()
        return self._cc_hull


    def compute(self):
        """
        compute the concave hull -- should be called after changing any parameters
        """
        if self.cvx_hull is None:
            self.cvx_hull = scipy.spatial.ConvexHull(self.points).vertices

        self._cc_hull = concave_hull(self.points,
                                     self.cvx_hull,
                                     self.concavity,
                                     self.length_threshold)





# def concaveman2d(points, hull=None, concavity=2.0, length_threshold=0.0):
#     points = np.asrray(points, dytpe=np.float64)

#     if len(points.shape) != 2 or points.shape[0] != 2:
#         raise ValueError('points must be an Nx2 array')

#     if hull is not None:
#         if len(hull.shape) != 1:
#             raise ValueError('hull must be a 1-D array')

#         if np.any(hull >= len(points)) or np.any(hull < 0):
#             raise ValueError('hull indices out of bounds')
#     else:
#         # compute the convex hull
#         # compute the ConvexHull
#         hull = scipy.spatial.ConvexHull(points).vertices

#     return concave_hull(points, hull, concavity, length_threshold)
