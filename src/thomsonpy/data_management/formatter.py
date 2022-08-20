# -*- coding: utf-8 -*-
"""
.. module:: formatter
        :platform: Unix
        :synopsis: tools for manage data format, usage and storing 
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

import pickle
import numpy as np
import thomsonpy.constants.units as units
import thomsonpy.data_management.octree.octree as octr
import thomsonpy.data_management.spherical_mesh.spherical_mesh as sphmesh


def dump(filepath, obj):
    """
    Auxiliary function for storing data.

    :param filepath: name for the file storing the data
    :type filepath: string
    :param obj: object containing the data that will be stored
    :type obj: any 
    """
    f = open(filepath, "wb")
    pickle.dump(obj, f)
    f.close()


def load(filepath):
    """
    Auxiliary function for loading stored data.
    
    :param filepath: name for the file storing the data
    :type filepath: string

    :return: the object stored in filepath
    :rtype: any
    """
    f = open(filepath, "rb")
    obj = pickle.load(f)
    f.close()
    return obj


def cartesian_to_spherical(coords):
    x = coords[0]
    y = coords[1]
    z = coords[2]
    radial = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = np.arctan(x / z)
    theta = np.arccos(y / radial)
    return np.array([phi, theta, radial])


def spherical_to_cartesian(r, theta, phi):
    """
    tangential_intensity gets cartesian coordinates from spherical coordinates.

    .. math::
        x &= \sin{\\theta} \sin{\phi} 
        
        y &= \cos{\\theta} 
        
        z &= \sin{\\theta} \cos{\phi}

    :param r: radius of the spherical coordinates
    :type r: float
    :param theta: latitude, starting from the North [0 - PI]
    :type theta: float
    :param phi: longitude: Carrington Longitude [0 - 2PI]
    :type phi: float

    :return: coordinates in the Cartesian Coordinate System
    :rtype: numpy.ndarray[float, float, float]
    """
    x = r * np.sin(theta) * np.sin(phi)
    y = r * np.cos(theta)
    z = r * np.sin(theta) * np.cos(phi)
    return np.array([x, y, z])


def apply_octree_data_format(r, theta, phi, ne_mas):
    """
    tangential_intensity gives the required data format by the octree structure.
    
    The original data is stored in hdf format in the Spherical Coordinates System
    in mas units. The octree structure needs the data stored as a list of Data 
    objects where each one stores its coordinates in the Cartesian Coordinates 
    System and its physical magnitudes associated in the International System of 
    Units. 
    
    :param r: radius coordinate of the Spherical Coordinates System in RSol
    :type r: float
    :param theta: theta coordinate of the Spherical Coordinates System in rad
    :type theta: float
    :param phi: phi coordinate of the Spherical Coordinates System in rad
    :type phi: float
    :param ne_mas: ne value for the given coordinates in MAS units.
    :type ne_mas: float
    
    :return: a Data object with coordinates in the Cartesian Coordinates System and
            with units in the International System of Units.
    :rtype: :class:`thomsonpy.data_management.octree.octree.Data`
    
    """
    coords = spherical_to_cartesian(r, theta, phi) * units.RSOL_TO_METERS  # From RSol to m
    ne = ne_mas * units.NE_MAS_FACTOR  # From MAS to m⁻³.
    data = octr.Data(coords, ne)
    return data


def apply_spherical_mesh_data_format(r, theta, phi, ne_mas):
    coords = spherical_to_cartesian(r, theta, phi) * units.RSOL_TO_METERS  # From RSol to m
    ne = ne_mas * units.NE_MAS_FACTOR  # From MAS to m⁻³.
    data = sphmesh.Data(coords, ne)
    return data
