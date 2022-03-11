import numpy as np
import thomsonpy.constants.units as units
import thomsonpy.config.octree_params as op

"""
Coordinate system centered on the center of the Sun defined by the right-hand rule.
    - Y: vertical axis (|).
    - X: horizontal axis (---).
    - Z: axis aligned with the center of the Sun and the Observer point (·).

All parameters are in the International System of Units.
"""
SUN_CENTER = np.array([0, 0, 0]) # m
OBSERVER = np.array([0, 0, -units.AU_TO_METERS]) # m
SOLAR_RADIUS = units.RSOL_TO_METERS # m

"""
Image parameters.

All parameters are in the International System of Units.
"""
MAX_VIS_R = 1.9223684 * units.RSOL_TO_METERS # m
MAX_COORD = np.sqrt(MAX_VIS_R**2 / 2) 
MIN_COORD = -np.sqrt(MAX_VIS_R**2 / 2)
IMAGE_MAX = np.array([MAX_COORD, MAX_COORD])
IMAGE_MIN = np.array([MIN_COORD, MIN_COORD])
IMAGE_RESOLUTION = 2000000 #  (2000 km)
IMAGE_NUM_POINTS = (MAX_COORD - MIN_COORD) / IMAGE_RESOLUTION