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
import matplotlib as mpl
import matplotlib.pyplot as plt
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

model = np.log(data[:, :, 100])

plt.figure(figsize=(10, 10))
plt.title("Thomson Scattering with Cramer Ne model")
plt.ylabel("Y (RSol)")
plt.xlabel("X (RSol)")
plt.imshow(model)

plt.colorbar()
plt.show()


fin_time = time.perf_counter()
print("Data fragmentation in", fin_time - ini_time, "seconds.")
