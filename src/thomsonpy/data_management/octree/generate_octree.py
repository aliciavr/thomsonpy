#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: generate_octree
        :platform: Unix
        :synopsis: generates and stores an octree
.. moduleauthor:: 
"""
import os
import pickle
import numpy as np
import open3d as o3d
import time
import thomsonpy.data_management.octree.octree as octr
import thomsonpy.config.octree_params as op
import thomsonpy.config.paths as paths
import thomsonpy.data_management.formatter as formatter


## OCTREE CREATION AND STORAGE
min_v = [op.MIN_1, op.MIN_2, op.MIN_3, op.MIN_4]
max_v = [op.MAX_1, op.MAX_2, op.MAX_3, op.MAX_4]

i = 0
# Loading data...
octree_data = formatter.load(f"{paths.OCTREE_DATA_PATH}{os.path.splitext(paths.PREDSCI_FILENAME)[0]}.data")
print("Loaded octree data:", len(octree_data), "points.")
octree = 0

# Creating octree...
print("Creating octree with params:")
print("MAX_LEVEL =", op.MAX_LEVEL)
print("MAX_DATA =", op.MAX_DATA)
print("MIN_V =", min_v[i])
print("MAX_V =", max_v[i])
ini_time = time.perf_counter()
octree = octr.Octree(op.MAX_LEVEL, op.MAX_DATA, octree_data, min_v[i], max_v[i])
fin_time = time.perf_counter()
print("Octree " + str(i + 1) + " built in", str((fin_time - ini_time) / 60), "minutes.")
octr.Octree.save(octree, paths.OCTREES_PATH + "octree_" + str(i + 1) + ".oct")
print("Octree " + str(i + 1) + " created at", paths.OCTREES_PATH)
