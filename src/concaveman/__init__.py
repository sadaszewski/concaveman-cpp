
__version__ = "0.1.0"

import numpy as np

from .concaveman import ConcaveHull


# an example square with a couple interior points
EXAMPLE_SQUARE = np.array(((5.0, 10.0),
                           (5.0, 20.0),
                           (10.0, 20.0),
                           (10.0, 10.0),
                           (6.0, 13.0),
                           (8.0, 18.0),
                           ))
