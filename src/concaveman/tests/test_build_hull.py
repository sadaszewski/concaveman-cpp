"""
test_build_hull.py
"""
import sys
import numpy as np
import scipy.spatial
from concaveman.build_hull import concave_hull
from concaveman import ConcaveHull, EXAMPLE_SQUARE


def test_build_hull():
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


def test_ConcaveHull():
    """
    test with defaults
    """
    points = EXAMPLE_SQUARE
    hull = ConcaveHull(EXAMPLE_SQUARE, concavity=2.0, length_threshold=0.0)

    print(hull.cvx_hull)

    # convex hull should be correct
    assert sorted(hull.cvx_hull) == [0, 1, 2, 3]

    print(hull.cc_hull)


    # same results as first got -- correct?
    assert hull.cc_hull.tolist() == [[5., 20.],
                                     [6., 13.],
                                     [5., 10.],
                                     [10., 10.],
                                     [10., 20.],
                                     ]

def test_ConcaveHull_inf_concavity():
    """
    test with defaults
    """
    points = EXAMPLE_SQUARE
    hull = ConcaveHull(EXAMPLE_SQUARE,
                       concavity=float("inf"),
                       length_threshold=0.0)

    print(hull.cvx_hull)

    # convex hull should be correct
    assert sorted(hull.cvx_hull) == [0, 1, 2, 3]

    print(hull.cc_hull)


    # should be the outer square
    assert hull.cc_hull.tolist() == [[5., 20.],
                                     [5., 10.],
                                     [10., 10.],
                                     [10., 20.],
                                     ]

def test_ConcaveHull_small_concavity():
    """
    test with defaults
    """
    points = EXAMPLE_SQUARE
    hull = ConcaveHull(EXAMPLE_SQUARE,
                       concavity= -5.0,
                       length_threshold=0.0)

    print(points)
    print(hull.cvx_hull)

    # convex hull should be correct
    assert sorted(hull.cvx_hull) == [0, 1, 2, 3]

    print(hull.cc_hull)


    # should be the outer square
    assert hull.cc_hull.tolist() == [[5., 20.],
                                     [5., 10.],
                                     [10., 10.],
                                     [10., 20.],
                                     ]


def check_memory():
    """
    not a real test, but something to run to make sure that there isn't a memory leak

    runs an infinite loop, so you can see if the memory use is growing.

    I ran this for a good while, and memory use was rock-steady :-)

    """
    # create a whole bunch of points
    N = 100000

    while True:
        print(f"computing a hull of {N} points")
        points = np.random.random((N, 2)) * 1000

        # compute the ConvexHull
        cvx_hull = scipy.spatial.ConvexHull(points)

        print(f"convex hull has {len(cvx_hull.vertices)} points")
        # print(cvx_hull.vertices)

        # call concave_hull
        cc_hull = concave_hull(points, cvx_hull.vertices)
        print(f"concave hull has {len(cc_hull)} points")

        # print(cc_hull)


if __name__ == "__main__":
    # test_simple()

    check_memory()




