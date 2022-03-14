#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:06:09 2022

@author: aliciavr
"""

import pickle
import numpy as np
import open3d as o3d
from octree import Octree, Data
import thomsonpy.config.octree_params as op
import thomsonpy.config.paths as paths
import thomsonpy.data_management.formatter as formatter
import time

## OCTREE CREATION AND STORAGE
min_v = [op.MIN_1, op.MIN_2, op.MIN_3, op.MIN_4]
max_v = [op.MAX_1, op.MAX_2, op.MAX_3, op.MAX_4]

for i in range(4):  
    if i != 0:
        # Loading data...
        octree_data = formatter.load(paths.OCTREE_DATA_PATH + "octree_data_" + str(i + 1) + ".obj")
        print("Loaded octree data:", len(octree_data), "points.")
        # Creating octree...
        print("Creating octree with params:")
        print("MAX_LEVEL =", op.MAX_LEVEL)
        print("MAX_DATA =", op.MAX_DATA)
        print("MIN_V =", min_v[i])
        print("MAX_V =", max_v[i])
        ini_time = time.perf_counter()
        octree = Octree(op.MAX_LEVEL, op.MAX_DATA, octree_data, min_v[i], max_v[i])
        fin_time = time.perf_counter()
        print("Octree " + str(i + 1) + " built in", str((fin_time - ini_time) / 60), "minutes.")
        Octree.save(octree, paths.OCTREE_OBJECTS_PATH + "octree_" + str(i + 1) + ".oct")
        print("Octree " + str(i + 1) + " created at", paths.OCTREE_OBJECTS_PATH)
"""
# Loading data...
f = open('../octree_data_2.obj', 'rb')
octree_data_2 = pickle.load(f)
f.close()
# Creating octree...
octree_2 = Octree(op.MAX_LEVEL, op.MAX_DATA, octree_data_2, op.MIN_2, op.MAX_2)
Octree.save(octree_2, "octree_2.obj")
# Deleting resources...
del(octree_data_2)
del(octree_2)
print("Octree 2 created.")

# Loading data...
f = open('../octree_data_3.obj', 'rb')
octree_data_3 = pickle.load(f)
f.close()
# Creating octree...
octree_3 = Octree(op.MAX_LEVEL, op.MAX_DATA, octree_data_3, op.MIN_3, op.MAX_3)
Octree.save(octree_3, "octree_3.obj")
# Deleting resources...
del(octree_data_3)
del(octree_3)
print("Octree 3 created.")

# Loading data...
f = open('../octree_data_4.obj', 'rb')
octree_data_4 = pickle.load(f)
f.close()
# Creating octree...
octree_4 = Octree(op.MAX_LEVEL, op.MAX_DATA, octree_data_4, op.MIN_4, op.MAX_4)
Octree.save(octree_4, "octree_4.obj")
# Deleting resources...
del(octree_data_4)
del(octree_4)
print("Octree 4 created.")
"""
"""
for d in octree_data_1:
    node = octree_1.search(d)
    if node == None:
        print(d, "NOT FOUND")
    else:
        print(node)
"""