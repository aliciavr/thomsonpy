from thomsonpy.data_management.octree.octree import Data
import numpy as np

def from_numpy_to_octree_data(xyz, ne):
    """
    Adapts the numpy data to the octree data format.
    
    Parameters
    ----------
    xyz : numpy.array[float, float, float]
        Numpy array of points with data.
    ne : float
        Electronic density for each point stored at xyz in cm⁻³.

    Returns
    -------
    data : list
        List of data in octree data format.
        """
    progress = 0
    total = xyz.size
    data = list()
    for i, ne_val in enumerate(ne):
        ne_val = ne_val * 1e6 # conversion to m⁻³.
        d = Data(xyz[i], ne_val)
        data.append(d)
        
        # progress...
        if progress % 500000 == 0:
            print(progress / total * 100, "%")
        progress += 1
        
    return data

def spherical_to_cartesian(r, theta, phi):
    """
    Gets cartesian coordinates form spherical coordinates.

    Parameters
    ----------
    r : float
        Radius of the spherical coordinates.
    theta : float
        Latitude: starting from the North [0 - PI].
    phi : float
        Longitude: Carrington Longitude [0 - 2PI].

    Returns
    -------
    x : float
        x coordinate in cartesian coordinates.
    y : float
        y coordinate in cartesian coordinates.
    z : float
        z coordinate in cartesian coordinates.

    """
    x = r * np.sin(theta) * np.cos(phi) 
    y = r * np.cos(theta)
    z = r * np.sin(theta) * np.sin(phi)
    return np.array([x, y, z])

def coords_system_change(coords, distance = 215):
    """
    Description
    -------------
    Changes the coordinate system between Sun Centre as System S and
    Observer Point as System O:
    From System S to System O  
    Xo = Xs
    Yo = Ys
    Zo = 215RSol - Zs

    From System O to System S
    Xs = Xo
    Ys = Yo
    Zs = 215RSol - Zo

    Parameters
    -------------
    coords: (float, float, float)
        Coordinates to be changed in either both sytems.
    distance: float
        Distance to be translated (the system of units dependes on the 
        user). It is set by default to 215 in RSun units, equivalent to a 1 AU.

    Returns
    -------------
    new_coords: (float, float, float)
        New coordinates according to the other system involved.
    """

    x, y, z = coords[0], coords[1], coords[2]
    z = distance - z
    new_coords = (x, y, z)
    return new_coords


import pickle
# FORMATTING AND STORAGE
# Octree 1
f = open('points_1', 'rb')
points_1 = pickle.load(f)
f.close()

f = open('ne_1', 'rb')
ne_1 = pickle.load(f)
f.close()     

octree_data_1 = from_numpy_to_octree_data(points_1, ne_1)
f = open('octree_data_1', 'wb')
pickle.dump(octree_data_1, f)
f.close()



import open3d as o3d
## VISUALIZATION
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(points_1)
geoms = pcd
o3d.visualization.draw_geometries(geoms)

"""
octree_data_2 = formatter.from_numpy_to_octree_data(points_2, ne_2)
f = open('octree_data_2', 'wb')
pickle.dump(octree_data_2, f)
f.close()

octree_data_3 = formatter.from_numpy_to_octree_data(points_3, ne_3)
f = open('octree_data_3', 'wb')
pickle.dump(octree_data_3, f)
f.close()

octree_data_4 = formatter.from_numpy_to_octree_data(points_4, ne_4)
f = open('octree_data_4', 'wb')
pickle.dump(octree_data_4, f)
f.close()


import open3d as o3d
## VISUALIZATION
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(points_1)
geoms = pcd
o3d.visualization.draw_geometries(geoms)

## VISUALIZATION
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(points_2)
geoms = pcd
o3d.visualization.draw_geometries(geoms)

## VISUALIZATION
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(points_3)
geoms = pcd
o3d.visualization.draw_geometries(geoms)

## VISUALIZATION
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(points_4)
geoms = pcd
o3d.visualization.draw_geometries(geoms)
"""