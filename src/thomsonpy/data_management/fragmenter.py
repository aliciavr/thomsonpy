# -*- coding: utf-8 -*-
"""
.. module:: fragmenter
        :platform: Unix
        :synopsis: it has tools for selecting and fragmenting the original data.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

import numpy as np
from pyhdf.SD import SD
import thomsonpy.data_management.formatter as formatter


def get_ne_raw(filepath):
    """
    tangential_intensity loads and the raw model of electron density provided by Predictive Science Inc.
    
    :param filepath: filepath to the file storing the raw model.
    :type filepath: string
    
    :return: the electron density raw model.
    :rtype: np.array([int, int, int])
    """
    hdf = SD(filepath)
    print("\nLoading datasets from", filepath, ":\n", hdf.datasets())
    data = hdf.select(3).get()  # DATA CUBE OF ELECTRON DENSITY (RHO)

    num_points = data.size
    print("# POINTS = ", num_points)

    return data


def get_ne_raw_coords(filepath, opt):
    """
    tangential_intensity gets the values for the coordinates indicated at ``opt`` parameter of the raw data cube provided by Predictive
    Science and that can be obtained with the function :py:func:`thomsonpy.data_management.fragmenter`.
    
    :param filepath: filepath to the file storing the raw model.
    :type filepath: string :param opt: it accepts the values ``phi``, ``theta`` or ``radial``, indicating the coordinate
    in the Spherical Coordinates Systems to be loaded.
    :type opt: string
    
    :return: the values for the chosen spherical coordinate. ``None`` if the ``opt`` parameter has an invalid value.
    :rtype: np.array([float])
    """

    hdf = SD(filepath)
    print("\nLoading datasets from", filepath, ":\n", hdf.datasets())

    coords = None
    if opt == "phi":
        coords = hdf.select(0).get()  # PHI (rad, 0-2PI) 699
        num_phi = coords.size
        print("\nPHI (rad, 0-2PI). # phi =", num_phi)
    elif opt == "theta":
        coords = hdf.select(1).get()  # THETA (rad, 0-PI) 327
        num_theta = coords.size
        print("THETA (rad, 0-PI). # theta =", num_theta)
    elif opt == "radial":
        coords = hdf.select(2).get()  # RADIAL (RSOL, 1-30) 288
        num_radial = coords.size
        print("RADIAL (RSOL, 1-30). # radial =", num_radial)

    return coords


def selection(r, theta, phi, ne):
    """
    
    :param r: radial coordinate in the Spherical System of Coordinates.
    :type r: float
    :param theta: theta coordinate in the Spherical System of Coordinates.
    :type theta: float
    :param phi: phi coordinate in the Spherical System of Coordinates.
    :type phi: float
    :param ne: electron density model.
    :type ne: float
    
    :return:
    :rtype: boolean
    
    """
    cartesian_coords = formatter.spherical_to_cartesian(r, theta, phi)
    x = cartesian_coords[0]
    y = cartesian_coords[1]
    z = cartesian_coords[2]
    limit = 3.5
    return x <= limit and x >= -limit and y <= limit and y >= -limit and z <= limit and z >= -limit


def fragment(selection_func, format_func, ne_raw, data_radial, data_theta, data_phi, progress_step=1e6):
    """
    
    :param selection_func: function defining the subspace of points.
    :type selection_func: function
    :param format_func: function defining the format of the output of this fragmentation.
    :type format_func: function
    :param ne_raw: original electron density model from Predictive Science Inc.
    :type ne_raw: 3d array
    :param data_radial: values for the radial coordinate.
    :type data_radial: array of float
    :param data_theta: values for the theta coordinate.
    :type data_theta: array of float
    :param data_phi: values for the phi coordinate.
    :type data_phi: array of float
    :param progress_step: int
    :type progress_step: int
    """

    it = np.nditer(ne_raw, flags=['multi_index'])
    num_points = ne_raw.size
    print(num_points)
    data = []
    # Progress and auxiliary variables
    progress = 0
    # Fragmentation process.
    for ne_mas in it:
        i, j, k = it.multi_index
        r = data_radial[k]
        theta = data_theta[j]
        phi = data_phi[i]

        if selection_func(r, theta, phi, ne_mas):
            d = format_func(r=r, theta=theta, phi=phi, ne_mas=ne_mas)
            data.append(d)
        if progress % progress_step == 0:
            print(progress / num_points * 100, "%")
        progress += 1

    return data
