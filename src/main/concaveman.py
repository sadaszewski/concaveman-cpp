import numpy as np

import cppimport
concaveman = cppimport.imp("cpp.concaveman")


def concaveman2d(points, hull, concavity=2.0, lengthThreshold=0.0):
    points = np.array(points).astype(np.double)
    hull = np.array(hull).astype(np.int32)

    if len(points.shape) != 2:
        raise ValueError('points must be a 2-D array')

    if len(hull.shape) != 1:
        raise ValueError('hull must be a 1-D array')

    if np.any(hull >= len(points)) or np.any(hull < 0):
        raise ValueError('hull indices out of bounds')

    concave_points = concaveman.concaveman(
        points, hull, 
        concavity, lengthThreshold)

    return concave_points


if __name__ == '__main__':
    print(concaveman2d([[0, 0], [.25, .15], [1, 0], [1, 1]], [0, 2, 3]))
