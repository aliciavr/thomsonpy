#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: fragmenter
        :platform: Unix
        :synopsis: selects and applies octree format to original data 
.. moduleauthor:: 
"""
import os
import numpy as np
import pickle
from pyhdf.SD import SD
import time
import thomsonpy.data_management.formatter as formatter
import thomsonpy.config.paths as paths

# DATA LOADING FROM PREDICTIVE SCIENCE FILES

filepath = paths.PREDSCI_DATA_PATH + paths.PREDSCI_FILENAME

hdf = SD(filepath)
print("\nLoading datasets from", filepath, ":\n", hdf.datasets())
data_phi = hdf.select(0).get()  # PHI (rad, 0-2PI) 699
data_theta = hdf.select(1).get()  # THETA (rad, 0-PI) 327
data_radial = hdf.select(2).get()  # RADIAL (RSOL, 1-30) 288
data = hdf.select(3).get()  # DATA CUBE OF ELECTRON DENSITY (RHO)

num_phi = data_phi.size
num_theta = data_theta.size
num_radial = data_radial.size                                                                               
num_points = data.size
print("\nPHI (rad, 0-2PI). # phi =", num_phi)
print("THETA (rad, 0-PI). # theta =", num_theta)
print("RADIAL (RSOL, 1-30). # radial =", num_radial)
print("# POINTS = ", num_points)

# DATA FRAGMENTING
print("\nStarting data formatting and fragmenting...")

it = np.nditer(data, flags=['multi_index'])

octree_data = []
# Progress and auxiliar variables
max_r = -1 
progress = 0
ini_time = time.perf_counter()

def selection(r, theta, phi, ne):
    cartesian_coords = formatter.spherical_to_cartesian(r, theta, phi)
    return cartesian_coords[0] >= 0 and cartesian_coords[1] >= 0 and r <= 3

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

fin_time = time.perf_counter()
print("Data formatting and fragmentation in", fin_time - ini_time, "seconds.")

formatter.dump(f"{paths.OCTREE_DATA_PATH}{os.path.splitext(paths.PREDSCI_FILENAME)[0]}.data", octree_data)
print(f"Stored data for octree at {paths.OCTREE_DATA_PATH}")

print("Octree will have", len(octree_data), "points")
