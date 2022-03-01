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

## CREATION

points = np.random.randint(0,1000,(1000,3)) 

cloud = list()

for p in points:
    cloud.append(Data(p, np.random.randint(50, 100)))

#my_octree = Octree(20, 100, cloud, np.array([0, 0, 0]), np.array([1000, 1000, 1000]))
my_octree = Octree.load('my_octree.obj')

print("hello")

## OUTPUT

print("Cloud data:")    
for d in cloud:
    print(d)
    
print("Searching each point:")
for d in cloud:
    print(my_octree.search(d))
    
    
## STORAGE

#Octree.save(my_octree, 'my_octree.obj')

## VISUALIZATION

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(points)
pcd.paint_uniform_color([1, 0.5, 0])
geoms = my_octree.get_visual_octree()
geoms.append(pcd)
o3d.visualization.draw_geometries(geoms)
