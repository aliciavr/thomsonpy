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
    It loads and the raw model of electron density provided by Predictive Science Inc.
    
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
    It gets the values for the coordinates indicated at ``opt`` parameter of the raw data cube provided by Predictive Science and that can be obtained with the function :py:func:`thomsonpy.data_management.fragmenter`.
    
    :param filepath: filepath to the file storing the raw model.
    :type filepath: string
    :param opt: it accepts the values ``phi``, ``theta`` or ``radial``, indicating the coordinate in the Spherical Coordinates Systems to be loaded.
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
    
    :param r:
    :type r: float
    :param theta:
    :type theta: float
    :param phi:
    :type phi: float
    :param ne:
    :type ne: float
    
    :return:
    :rtype: boolean
    
    """
    cartesian_coords = formatter.spherical_to_cartesian(r, theta, phi)
    return cartesian_coords[0] >= 0 and cartesian_coords[1] >= 0 and r <= 3

def fragment(selection_func, ne_raw, data_radial, data_theta, data_phi, progress_step = 1e6): 
    """
    
    :param selection_func:
    :type:
    
    """
    it = np.nditer(ne_raw, flags=['multi_index'])
    num_points = ne_raw.size
    print(num_points)
    octree_data = []
    # Progress and auxiliar variables
    max_r = -1 
    progress = 0
    # Fragmentation process.
    for ne_mas in it:
        i, j, k = it.multi_index
        r = data_radial[k]
        theta = data_theta[j]
        phi = data_phi[i]

        if selection(r, theta, phi, ne_mas) == True:
            d = formatter.apply_octree_data_format(r, theta, phi, ne_mas)
            octree_data.append(d)

            if r > max_r:
                max_r = r
        if progress % 1e6 == 0:
            print(progress / num_points * 100, "%")
        progress += 1
        
    return octree_data