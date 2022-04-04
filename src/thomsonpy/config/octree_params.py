import numpy as np
import thomsonpy.constants.units as units

# Maximum deep level of the octree.
MAX_LEVEL = 8
# Maximum amount of data in each leaf node.
MAX_DATA = 1000

# Maximum and minimum distance of octree data from the center of the Sun in RSol.
MAX_R = 3.1 * units.RSOL_TO_METERS# m

# Maximum bounds for octree according a MAX_R radius.
MAX_COORD = MAX_R
MIN_COORD = -MAX_R

# Quadrant 1: x >= 0 & y >= 0
MAX_1 = np.array([MAX_COORD, MAX_COORD, MAX_COORD])
MIN_1 = np.array([0, 0, MIN_COORD])

# Quadrant 2: x <= 0 & y >= 0
MAX_2 = np.array([0, MAX_COORD, MAX_COORD])
MIN_2 = np.array([MIN_COORD, 0, MIN_COORD])

# Quadrant 3: x <= 0 & y <= 0
MAX_3 = np.array([0, 0, MAX_COORD])
MIN_3 = np.array([MIN_COORD, MIN_COORD, MIN_COORD])

# Quadrant 4: x >= 0 & y <= 0
MAX_4 = np.array([MAX_COORD, 0, MAX_COORD])
MIN_4 = np.array([0, MIN_COORD, MIN_COORD])
