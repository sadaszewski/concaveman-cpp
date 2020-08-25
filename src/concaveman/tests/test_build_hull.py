"""
test_build_hull.py
"""
import numpy as np
from concaveman.build_hull import build_hull


def test_simple():
    """
    a simple call of the function
    """
    points = np.array([[0, 0], [.25, .15], [1, 0], [1, 1]], dtype=np.float64)
    hull = np.array([0, 2, 3], dtype=np.int32)

    cvx_hull = build_hull(points, hull)

    print(cvx_hull)

