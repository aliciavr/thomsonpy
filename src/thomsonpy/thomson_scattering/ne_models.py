"""
.. module:: ne_models
        :platform: Unix
        :synopsis: this module contains functions that manages electron density models readable by :py:mod:`thomsonpy`.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

import thomsonpy.config.thomson_scattering_params as tsp
import thomsonpy.config.solar_imager_params as sip
import thomsonpy.constants.units as units

def crammer_model(d):
    '''
    This model retrieves the electron number density profile obtained by Crammer et al. (1999) using the UVCS/WLC aboard SOHO satellite. Initially, the units are :math:`cm^{-3}`. Therefore a units conversion factor have been applied to convert from :math:`cm^{-3}` to :math:`m^{-3}`.

    :param d: distance QS.
    :type d: float
    
    .. warning:: 
        Be careful with the units of parameter ``d``. Its must be in :math:`RSol`. Otherwise, the model will not work properly.
    
    :return: electron number density profile at the given distance :math:`d`.
    :rtype: float
    '''
    # en cm⁻³ --> m⁻³
    return 1E8 * (3.89 * d**-10.5 + 0.00869 * d**-2.57) * 1E6

def predictive_science_model(z, TG, NE_MODEL):
    """
    This model retrieves the electron number density from a simulation made and published by **Predictive Science Inc.** for the **solar eclipse of 4th December, 2021** for a given distance :math:`z` from the point of view of the observer (:math:`O`) to the scattering point (:math:`Q`) with a predefined ``ThomsonGeometry`` object that manages the coherence of the rest of the geometric components for the ray-tracing process.
    
    .. warning:: 
        This model will not work without a model (``NE_MODEL``) structured in the ``Octree`` data structure managed in the :py:mod:`thomsonpy` package.

    :param z: distance OQ.
    :type z: float
    :param TG: an object managing the internal coherence of the geometric components of the Thomson scattering ray-tracing process.
    :type TG: :py:class:`thomsonpy.thomson_scattering.thomson_scattering_tools.ThomsonGeometry`
    :param NE_MODEL: an ``Octree`` containing a cloud of ``Data``` objects representing the cloud of points with physical measures provided by Predictive Science Inc.
    :type NE_MODEL: :py:class:`thomsonpy.data_management.octree.octree.Octree`
    
    :return: electron number density for the solar eclipse of 4th December, 2021.
    :rtype: float
    
    
    """
    target = TG.get_target(z)
    data = NE_MODEL.search_nearest(target)
    return data.get_ne() 