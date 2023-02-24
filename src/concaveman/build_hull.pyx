# cython: language_level=3

"""
wrapper of concaveman

NOTE: after working on this for a while, I realized that concaveman.cpp was
mostly a wrapper that was set up to call the templates from C, with
regular C arrays.

But we should be able to call the Templates directly from Cython, so save that
translation step.

but other than some wasted data copying, this works.

NOTE: I removed the function pointer for the free() call, as this code is now using that
      malloc'ed pointer to build the returning numpy array -- so should be freed by
      python when the numpy array is deleted.
"""

import numpy as np
cimport numpy as cnp

from libc.stdlib cimport free
from libc.string cimport memcpy


cdef extern from "pyconcaveman.h":

    void pyconcaveman2d(double *points_c,
                        size_t num_points,
                        int *hull_points_c,
                        size_t num_hull_points,
                        double concavity,
                        double lengthThreshold,
                        double **concave_points_c,
                        size_t *num_concave_points,
                        )

cpdef concave_hull(cnp.float64_t[:,:] points,
                   int[:] hull,
                   double concavity=2.0,
                   double length_threshold=0.0):
    """
    compute the concave hull of a bunch of points

    :param points: The points to make the hull out of.
    :type points: Nx2 array of float64 type

    :param hull: The convex hull of the points
    :type points: 1D array of int32

    :param concavity=2.0: concavity parameter: large means smoother
    :type concavity: python float

    :param length_threshold:
    :type length_threshold: python float
    """

    if points.shape[1] != 2:
        raise ValueError('points must be an Nx2 array')

    # now a memoryview, so any isn't helpful here.
    # should this check be at a higher level?
    # if np.any(hull >= len(points)) or np.any(hull < 0):
    #     raise ValueError('hull indices out of bounds')

    cdef cnp.float64_t* p_concave_points = NULL
    cdef size_t num_concave_points = 0

    # print("in cython: about to call pyconcaveman2d")

    pyconcaveman2d(&points[0, 0],
                   len(points),
                   &hull[0],
                   len(hull),
                   concavity,
                   length_threshold,
                   &p_concave_points,
                   &num_concave_points,
                   )

    # print("cpp concave hull returned")
    # print("num concave points:", num_concave_points)

    # create an array to return
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] arr_concave_points
    arr_concave_points = np.zeros((num_concave_points, 2), dtype=np.float64)

    # loping through element by element
    # cdef cnp.float64_t* temp = &arr_concave_points[0, 0]
    # for i in range(num_concave_points * 2):
    #     temp[i] = p_concave_points[i]

    # copy the data to the ndarray's memory block directly
    memcpy(&arr_concave_points[0, 0],
           p_concave_points,
           sizeof(cnp.float64_t) * num_concave_points * 2)

    # print('in cython again: concave_points:\n', arr_concave_points)

    # As we are compiling the C++ code all together, the
    # plain free() should work.
    # on OS-X, it doesn't crash nor leak memory.
    free(p_concave_points)

    return arr_concave_points
