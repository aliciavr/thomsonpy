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


"""
## CREATION

points = np.random.randint(0,1000,(1000,3)) 

cloud = list()

for p in points:
    cloud.append(Data(p, np.random.randint(50, 100)))

my_octree = Octree(20, 100, cloud, np.array([0, 0, 0]), np.array([1000, 1000, 1000]))
#my_octree = Octree.load('my_octree.obj')

## OUTPUT

print("Cloud data:")    
for d in cloud:
    print(d)
    
print("Searching each point:")
for d in cloud:
    print(my_octree.search(d))
"""

"""
## DATA LOADING
f = open('../octree_data', 'rb')
data_cloud = pickle.load(f)
f.close()         

print("Data cloud loaded.")

f = open('../xyz', 'rb')
point_cloud = pickle.load(f)
f.close()         

print("Point cloud loaded.")
"""
## OCTREE CREATION
point_cloud = np.random.uniform(op.MIN_COORD, op.MAX_COORD,(1000,3)) 

data_cloud = list()

for p in point_cloud:
    print(p)
    data_cloud.append(Data(p, np.random.uniform(50, 100)))

octree_1 = Octree(op.MAX_LEVEL, op.MAX_DATA, data_cloud, op.MIN_1, op.MAX_1)
octree_2 = Octree(op.MAX_LEVEL, op.MAX_DATA, data_cloud, op.MIN_2, op.MAX_2)
octree_3 = Octree(op.MAX_LEVEL, op.MAX_DATA, data_cloud, op.MIN_3, op.MAX_3)
octree_4 = Octree(op.MAX_LEVEL, op.MAX_DATA, data_cloud, op.MIN_4, op.MAX_4)

print("My octrees created.")
    
## STORAGE
Octree.save(octree_1, "octree_1.obj")
Octree.save(octree_2, "octree_2.obj")
Octree.save(octree_3, "octree_3.obj")
Octree.save(octree_4, "octree_4.obj")

print("My octrees stored.")

## VISUALIZATION
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(point_cloud)
pcd.paint_uniform_color([1, 0.5, 0])
geoms = octree_1.get_visual_octree() + octree_2.get_visual_octree() + octree_3.get_visual_octree() + octree_4.get_visual_octree()
geoms.append(pcd)
o3d.visualization.draw_geometries(geoms)
