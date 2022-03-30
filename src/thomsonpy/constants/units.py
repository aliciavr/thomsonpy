# -*- coding: utf-8 -*-
"""
.. module:: units
        :platform: Unix
        :synopsis: fundamental unit conversion factors for computations 
.. moduleauthor:: 
"""

RSOL_TO_METERS = 696340000
METERS_TO_RSOL = 1/696340000
AU_TO_RSOL = 215
RSOL_TO_AU = 1/215
AU_TO_METERS = AU_TO_RSOL * RSOL_TO_METERS
METERS_TO_AU = METERS_TO_RSOL * RSOL_TO_AU
ANGSTROM_TO_METERS = 1E-10
METERS_TO_ANGSTROM = 1E10

NE_MAS_FACTOR = 1E14 # From MAS to m⁻³