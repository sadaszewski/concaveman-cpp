# cython: language_level=3

"""
wrapper of concaveman

NOTE: after working on this for a while, I realized that concaveman.cpp was
mostly a wrapper that was set up to call the templates from C, with
regular C arrays.

But we should be able to call the Templates directly from Cython, so save that tranlsation step.

Though maybe getting it to stop segfaulting first would be a good step
"""

import numpy as np
cimport numpy as cnp

from libc.stdlib cimport free

cdef extern from "pyconcaveman.h":

    void pyconcaveman2d(double *points_c,
                        size_t num_points,
                        int *hull_points_c,
                        size_t num_hull_points,
                        double concavity,
                        double lengthThreshold,
                        double **concave_points_c,
                        size_t *num_concave_points,
                        # void(**p_free)(void*),
                        )

# defining the function pointer type
# ctypedef void (**f_type)(void*)

# def concave_hull(cnp.ndarray[double, ndim=2, mode="c"] points,
#                  cnp.ndarray[cnp.int32_t, ndim=2, mode="c"] hull,
#                  double concavity=2.0,
#                  double length_threshold=0.0):

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

    # points = np.array(points).astype(np.double)
    # hull = np.array(hull).astype(np.int32)

    if points.shape[1] != 2:
        raise ValueError('points must be an Nx2 array')

    # if len(hull.shape) != 1:
    #     raise ValueError('hull must be a 1-D array')

    # if np.any(hull >= len(points)) or np.any(hull < 0):
    #     raise ValueError('hull indices out of bounds')

    cdef double* p_concave_points = NULL
    cdef size_t[1] num_concave_points
    num_concave_points[0] = 2
    # cdef f_type p_free = NULL

    print("num concave points:", num_concave_points[0])
    print("in cython: about to call pyconcaveman2d")

    pyconcaveman2d(&points[0, 0],
                   len(points),
                   &hull[0],
                   len(hull),
                   concavity,
                   length_threshold,
                   &p_concave_points,
                   num_concave_points,
                   # p_free). # should't need this, as we're c omiling with same lib.
                   )

    print("cpp concave hull returned")
    print("num concave points:", num_concave_points[0])

    #cdef cnp.float64_t[:, :] concave_points_mview = p_concave_points
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] arr_concave_points
    arr_concave_points = np.zeros((num_concave_points[0], 2), dtype=np.float64)
    cdef unsigned int i
    print(p_concave_points[0])
    for i in range(num_concave_points[0]):
    #     arr_concave_points[i, 0] = p_concave_points[i]
    #     arr_concave_points[i, 1] = p_concave_points[i]
        arr_concave_points[i, 0] = p_concave_points[i * 2]
        arr_concave_points[i, 1] = p_concave_points[i * 2 + 1]

    print('in cython again: concave_points:', arr_concave_points)

    # fixme: need to make sure this isn't a memory leak!
    # as we are compiling the C++ code all together, the
    # plan free() should work. I think.
    # free(p_concave_points)

    return arr_concave_points

