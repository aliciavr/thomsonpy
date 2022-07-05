"""
.. module:: solar_imager_params
        :platform: Unix
        :synopsis: parameters to configure the solar imager
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>

.. note::
    The coordinates system is centered on the center of the Sun and defined by the right-hand rule.
        * **Y**: vertical axis (|).
        * **X**: horizontal axis (---).
        * **Z**: axis aligned with the center of the Sun and the Observer point (·).
"""

import numpy as np
import thomsonpy.constants.units as units
import thomsonpy.config.octree_params as op
import thomsonpy.thomson_scattering.thomson_scattering_tools as tsp

SUN_CENTER = np.array([0, 0, 0])
""" Position of the centre of the Sun in :math:`m`. """
OBSERVER = np.array([0, 0, units.AU_TO_METERS])
""" Position of the observer (the Earth) in :math:`m`."""

MAX_VIS_R = 3 * units.RSOL_TO_METERS
""" Maximum radius of image visualization in :math:`m`. """
MAX_COORD = MAX_VIS_R# * np.sqrt(1/2)
""" Maximum common value for the coordinates of the image in :math:`m`. """
MIN_COORD = - MAX_VIS_R
""" Minimum common value for the coordinates of the image in :math:`m`. """

IMAGE_MAX = np.array([MAX_COORD, MAX_COORD])
""" Maximum coordinates of the image in :math:`m`. """

IMAGE_MIN = np.array([MIN_COORD, MIN_COORD])
""" Minimum coordinates of the image in :math:`m`. """

IMAGE_RESOLUTION = 7250000 #  (7250 km)
""" Image resolution in :math:`m`. """

IMAGE_NUM_POINTS = int(np.rint((MAX_COORD - MIN_COORD) / IMAGE_RESOLUTION))
""" Number of points of the image according to the coordinates (``MAX_COORD`` and ``MIN_COORD``) and the ``IMAGE_RESOLUTION``. """