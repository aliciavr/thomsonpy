"""
.. module:: units
        :platform: Unix.
        :synopsis: fundamental unit conversion factors needed for computations. 
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>.

.. note:: The constant ``NE_MAS_FACTOR`` is needed for data formatting processes when the data have been obtained
from Predictive Science Inc. simulations which are in raw format (MAS units).

"""

RSOL_TO_METERS = 696340000
""" Unit conversion factor: from Solar radius (:math:`RSol`) to meters (:math:`m`)."""
METERS_TO_RSOL = 1 / 696340000
""" Unit conversion factor: from meters (:math:`m`) to Solar radius (:math:`RSol`)."""
AU_TO_RSOL = 215
""" Unit conversion factor: from :math:`AU` to  Solar radius (:math:`RSol`). """
RSOL_TO_AU = 1 / 215
""" Unit conversion factor: from Solar radius (:math:`RSol`) to :math:`AU`. """
AU_TO_METERS = AU_TO_RSOL * RSOL_TO_METERS
""" Unit conversion factor: from :math:`AU` to meters (:math:`m`). """
METERS_TO_AU = METERS_TO_RSOL * RSOL_TO_AU
""" Unit conversion factor: from meters (:math:`m`) to :math:`AU`. """
ANGSTROM_TO_METERS = 1E-10
""" Unit conversion factor: from ángstrom (:math:`\mathring{A}`) to meters (:math:`m`). """
METERS_TO_ANGSTROM = 1E+10
""" Unit conversion factor: from meters (:math:`m`) to Ángstroms (:math:`\mathring{A}`)."""

NE_MAS_FACTOR = 1E14  # From MAS to m⁻³
"""Unit conversion factor: from ne raw coordinates (:math:`rho`) of Predictive Science Inc. to meters (:math:`m^{
-3}`). """
T_MAS_FACTOR = 28.07067E+6  # From MAS to K.
