"""
.. module:: thomson_scattering_params
        :platform: Unix
        :synopsis: parameters to set up and configure the Thomson scattering tools.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

from math import radians
import thomsonpy.constants.units as units
import thomsonpy.data_management.octree.octree as octr
import thomsonpy.config.paths as paths

SIGMA_E = 7.95E-30
""" Differential cross section for perpendicular Thomson scattering in :math:`m^2 sr^{-1}`. """

WAVELENGTH = 5000E-10
""" Observed wavelength in :math:`m`. """

T_SOL = 5778
""" Solar temperature in :math:`K`. """

X = 1 * units.AU_TO_METERS
""" Distance between the Observer (the Earth) and the Sun (distance :math:`OS`) in :math:`m`. """

EPSILON = (radians(1))
""" Elongation angle in :math:`rad`. """

SOLAR_RADIUS = units.RSOL_TO_METERS 
""" Solar radius in :math:`m`. """

INI_Z = 213 * units.RSOL_TO_METERS
""" Initial/Lower limit for numerical integration of the Thomson scattering in :math:`m`. """

FIN_Z = 217 * units.RSOL_TO_METERS
""" Final/Upper limit for numerical integration of the Thomson scattering in :math:`m`. """

INCR_Z = 0.25 * units.RSOL_TO_METERS
""" Step for numerical integration of the Thomson scattering in :math:`m`. """

NUM_Z = (FIN_Z - INI_Z) / INCR_Z
""" Number of points to be computed in numerical integration of the Thomson scattering in :math:`m`. """