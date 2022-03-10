#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 11:52:57 2021

@author: aliciavr
"""
import numpy as np
import pickle
from pyhdf.SD import SD
import time
import thomsonpy.data_management.formatter as formatter
import thomsonpy.config.paths as paths

# DATA LOADING FROM PREDICTIVE SCIENCE FILES

rho_file = paths.NE_PREDSCI_DATA_PATH

hdf = SD(rho_file)
print("\nLoading datasets from", paths.NE_PREDSCI_DATA_PATH, ":\n", hdf.datasets())
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
print("\nStarting data fragmenting...")

it = np.nditer(data, flags=['multi_index'])

points_1 = list()
ne_1 = list()
points_2 = list()
ne_2 = list()
points_3 = list()
ne_3 = list()
points_4 = list()
ne_4 = list()

# Progress and auxiliar variables
max_r = -1 
progress = 0
ini_time = time.perf_counter()

# Fragmentation process.
for ne in it:
    i, j, k = it.multi_index
    if ne > 1e-1:
        r = data_radial[k]
        theta = data_theta[j]
        phi = data_phi[i]
        coords = formatter.spherical_to_cartesian(r, theta, phi)
        if coords[0] >= 0 and coords[1] >= 0:
            # Octree 1: x >=0 & y >= 0
            points_1.append(coords)
            ne_1.append(ne)
        elif coords[0] <= 0 and coords[1] >= 0:
            # Octree 2: x <= 0 & y >= 0
            points_2.append(coords)
            ne_2.append(ne)
        elif coords[0] <= 0 and coords[1] <= 0:
            # Octree 3: x <= 0 & y <= 0
            points_3.append(coords)
            ne_3.append(ne)
        elif coords[0] >= 0 and coords[1] <= 0:
            # Octree 4: x >= 0 & y <= 0
            points_4.append(coords)
            ne_4.append(ne)
        
        if r > max_r:
            max_r = r
    if progress % 3000000 == 0:
        print(progress / num_points * 100, "%")
    progress += 1

fin_time = time.perf_counter()
print("Data fragmentation in", fin_time - ini_time, "seconds.")

print("\n# POINTS =", len(ne_1) + len(ne_2) + len(ne_3) + len(ne_4), " with ne > 1e-1 and max_r =", max_r, "RSol")
print("Octree 1 has", len(ne_1), "points")
print("Octree 2 has", len(ne_2), "points")
print("Octree 3 has", len(ne_3), "points")
print("Octree 4 has", len(ne_4), "points")

# STORAGE
# Octree 1
formatter.dump(paths.OCTREE_DATA_PATH + "points_1.obj", np.array(points_1))
formatter.dump(paths.OCTREE_DATA_PATH + "ne_1.obj", np.array(ne_1))

# Octree 2
formatter.dump(paths.OCTREE_DATA_PATH + "points_2.obj", np.array(points_2))
formatter.dump(paths.OCTREE_DATA_PATH + "ne_2.obj", np.array(ne_2))

# Octree 3
formatter.dump(paths.OCTREE_DATA_PATH + "points_3.obj", np.array(points_3))
formatter.dump(paths.OCTREE_DATA_PATH + "ne_3.obj", np.array(ne_3))

# Octree 4
formatter.dump(paths.OCTREE_DATA_PATH + "points_4.obj", np.array(points_4))
formatter.dump(paths.OCTREE_DATA_PATH + "ne_4.obj", np.array(ne_4))

print("Stored points and ne values for all octrees at", paths.OCTREE_DATA_PATH)

