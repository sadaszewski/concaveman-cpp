
import numpy as np
cimport numpy as cnp

# from libc.stdlib cimport malloc, free

cdef extern from "concaveman.h":

    void pyconcaveman2d(double *points_c,
                        size_t num_points,
                        int *hull_points_c,
                        size_t num_hull_points,
                        double concavity,
                        double lengthThreshold,
                        double **concave_points_c,
                        size_t *num_concave_points,
                        void(**p_free)(void*),
                        )

# defining the function pointer type
ctypedef void (**f_type)(void*)

# def concave_hull(cnp.ndarray[double, ndim=2, mode="c"] points,
#                  cnp.ndarray[cnp.int32_t, ndim=2, mode="c"] hull,
#                  double concavity=2.0,
#                  double length_threshold=0.0):
def concave_hull(cnp.float64_t [:,:] points,
                 int [:] hull,
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

    cdef double** p_concave_points
    cdef size_t p_num_concave_points
    cdef f_type p_free

    # p_concave_points_c = _ffi.new('double**')
    # p_num_concave_points = _ffi.new('size_t*')
    # p_free = _ffi.new('void(**)(void*)')


    # points_c = _ffi.cast('double*', points.ctypes.data)
    # hull_c = _ffi.cast('int*', hull.ctypes.data)
    pyconcaveman2d(&points[0,0],
                   len(points),
                   &hull[0],
                   len(hull),
                   concavity,
                   length_threshold,
                   p_concave_points,
                   &p_num_concave_points,
                   p_free)

    # _lib.pyconcaveman2d(points_c, len(points),
    #     hull_c, len(hull),
    #     concavity, lengthThreshold,
    #     p_concave_points_c, p_num_concave_points,
    #     p_free)

    num_concave_points = p_num_concave_points
    concave_points_c = p_concave_points

    # buffer = _ffi.buffer(concave_points_c, 8 * 2 * num_concave_points)

    # concave_points = np.frombuffer(buffer, dtype=np.double)
    # concave_points = concave_points.reshape((num_concave_points, 2))
    # concave_points = concave_points.copy()

    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] concave_points
    concave_points = np.empty((p_num_concave_points, 2), dtype=np.float64)
    cdef unsigned int i
    for i in range(len(concave_points)):
        concave_points[i, 0] = p_concave_points[i][0]
        concave_points[i, 1] = p_concave_points[i][1]

    print('concave_points:', concave_points)



    # fixme: need to make sure this isn't a memory leak!
    # p_free[0](concave_points_c)

    return concave_points
