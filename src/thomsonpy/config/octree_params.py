"""
.. module:: octree_params
        :platform: Unix
        :synopsis: parameters for managing the octree data structure
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>

.. note::
    The space is often divided in quadrants so we define our quadrants as:
        * Quadrant 1: x >= 0 & y >= 0
        * Quadrant 2: x <= 0 & y >= 0
        * Quadrant 3: x <= 0 & y <= 0
        * Quadrant 4: x >= 0 & y <= 0
"""

import numpy as np
import thomsonpy.constants.units as units

MAX_LEVEL = 6
""" Maximum deep level of the octree. """

MAX_DATA = 1000
""" Maximum amount of data in each leaf node. """

MAX_R = 5 * units.RSOL_TO_METERS
""" Maximum and minimum distance of octree data from the center of the Sun in :math:`m`. """

MAX_COORD = MAX_R
""" Maximum bounds for octree according a ``MAX_R`` radius in :math:`m`. """

MIN_COORD = - MAX_R
""" Minimum bounds for octree according a ``MAX_R`` radius in :math:`m`. """

# Quadrant 1: x >= 0 & y >= 0
MAX_1 = np.array([MAX_COORD, MAX_COORD, MAX_COORD / 2])
""" Maximum coordinates for the first quadrant in :math:`m`. """
MIN_1 = np.array([0, 0, MIN_COORD / 2])
""" Minimum coordinates for the first quadrant in :math:`m`. """

# Quadrant 2: x <= 0 & y >= 0
MAX_2 = np.array([0, MAX_COORD, MAX_COORD / 2])
""" Maximum coordinates for the second quadrant in :math:`m`. """
MIN_2 = np.array([MIN_COORD, 0, MIN_COORD / 2])
""" Minimum coordinates for the second quadrant in :math:`m`. """

# Quadrant 3: x <= 0 & y <= 0
MAX_3 = np.array([0, 0, MAX_COORD / 2])
""" Maximum coordinates for the third quadrant in :math:`m`. """
MIN_3 = np.array([MIN_COORD, MIN_COORD, MIN_COORD / 2])
""" Minimum coordinates for the third quadrant in :math:`m`. """

# Quadrant 4: x >= 0 & y <= 0
MAX_4 = np.array([MAX_COORD, 0, MAX_COORD / 2])
""" Maximum coordinates for the fourth quadrant in :math:`m`. """
MIN_4 = np.array([0, MIN_COORD, MIN_COORD / 2])
""" Minimum coordinates for the fourth quadrant in :math:`m`. """
