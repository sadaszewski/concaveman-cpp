"""
test_build_hull.py
"""
import sys
import numpy as np
from concaveman.build_hull import concave_hull


def test_simple():
    """
    a simple call of the function
    """
    points = np.array([[0, 0], [.25, .15], [1, 0], [1, 1]], dtype=np.float64)
    hull = np.array([0, 2, 3], dtype=np.int32)

    print("about to call concave_hull")
    cvx_hull = concave_hull(points, hull)
    print("concave_hull returned")
    print(cvx_hull)

    assert np.alltrue(cvx_hull == [[0.25, 0.15],
                                   [0., 0.],
                                   [1., 0.],
                                   [1., 1.],
                                   ])


if __name__ == "__main__":
    test_simple()

