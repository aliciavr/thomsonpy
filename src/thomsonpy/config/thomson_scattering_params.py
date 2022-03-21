# -*- coding: utf-8 -*-
from math import radians
import thomsonpy.constants.units as units
import thomsonpy.data_management.octree.octree as octr
import thomsonpy.config.paths as paths

SIGMA_E = 7.95E-30 # Differential cross section for perpendicular scattering in m²·sr⁻¹.
WAVELENGTH = 5000E-10 # Target wavelength (5000A) in m.
T_SOL = 5778 # Solar temperature in K.
X = 1 * units.AU_TO_METERS # Distance between the Observer (Earth) and the Sun: O-S.
EPSILON = (radians(10)) # Elongation in rad.
SOLAR_RADIUS = units.RSOL_TO_METERS # m

INI_Z = 210 * units.RSOL_TO_METERS # m
FIN_Z = 220 * units.RSOL_TO_METERS # m
INCR_Z = 0.1 * units.RSOL_TO_METERS # m
NUM_Z = (FIN_Z - INI_Z) / 0.1