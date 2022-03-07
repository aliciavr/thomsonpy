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

## OCTREE CREATION AND STORAGE

# Loading data...
f = open('../octree_data_1.obj', 'rb')
octree_data_1 = pickle.load(f)
f.close()
# Creating octree...
octree_1 = Octree(op.MAX_LEVEL, op.MAX_DATA, octree_data_1, op.MIN_1, op.MAX_1)
Octree.save(octree_1, "octree_1.obj")
# Deleting resources...
del(octree_data_1)
del(octree_1)
print("Octree 1 created.")

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
