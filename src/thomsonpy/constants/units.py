"""
.. module:: units
        :platform: Unix.
        :synopsis: fundamental unit conversion factors needed for computations. 
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>.

.. note::
    The constants:
        * NE_MAS_FACTOR.
        * ``NE_MAS_FACTOR``.
        * :py:const:`NE_MAS_FACTOR`.
        * :py:var:`NE_MAS_FACTOR`.
    Are needed for data formatting processes when the data have been obtained from Predictive Science Inc. simulations which are in raw format (MAS units).

"""

RSOL_TO_METERS = 696340000
""" 
:var RSOL_TO_METERS: unit conversion factor: from Solar radius (:math:`RSol`) to meters (:math:`m`)."""
METERS_TO_RSOL = 1/696340000
""" Unit conversion factor: from Solar radius to meters. """
AU_TO_RSOL = 215
""" Unit conversion factor: from Solar radius to meters. """
RSOL_TO_AU = 1/215
""" Unit conversion factor: from Solar radius to meters. """
AU_TO_METERS = AU_TO_RSOL * RSOL_TO_METERS
""" Unit conversion factor: from Solar radius to meters. """
METERS_TO_AU = METERS_TO_RSOL * RSOL_TO_AU
""" Unit conversion factor: from Solar radius to meters. """
ANGSTROM_TO_METERS = 1E-10
""" Unit conversion factor: from Solar radius to meters. """
METERS_TO_ANGSTROM = 1E10
""" Unit conversion factor: from Solar radius to meters. """


NE_MAS_FACTOR = 1E14 # From MAS to m⁻³
""" Unit conversion factor: from ne raw data (:math:`rho`) of Predictive Science Inc. to meters (:math:`m^{-3}`). """