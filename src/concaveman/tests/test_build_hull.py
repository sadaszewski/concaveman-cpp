"""
test_build_hull.py
"""
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

    sys.stdout.flush()
    print(cvx_hull)

if __name__ == "__main__":
    test_simple()

