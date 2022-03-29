import thomsonpy.config.paths as paths
from thomsonpy.data_management.octree.octree import Data
import pickle
import numpy as np

def dump(filename, obj):
    f = open(filename, "wb")
    pickle.dump(obj, f)
    f.close()

def load(filename):
    f = open(filename, "rb")
    obj = pickle.load(f)
    f.close()
    return obj

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
    total = ne.size
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
    x = r * np.sin(theta) * np.sin(phi) 
    y = r * np.cos(theta)
    z = r * np.sin(theta) * np.cos(phi)
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

def format_data():
    # FORMATTING AND STORAGE
    # Octree 1
    points_1 = load(paths.OCTREE_DATA_PATH + "points_1.obj")
    print("Loaded points for octree 1 from", paths.OCTREE_DATA_PATH + "points_1.obj")

    ne_1 = load(paths.OCTREE_DATA_PATH + "ne_1.obj")
    print("Loaded ne values for octree 1 from", paths.OCTREE_DATA_PATH + "ne_1.obj")

    octree_data_1 = from_numpy_to_octree_data(points_1, ne_1)
    dump(paths.OCTREE_DATA_PATH + "octree_data_1.obj", octree_data_1)
    print("Stored data for octree 1 at", paths.OCTREE_DATA_PATH + "octree_data_1.obj")

    del points_1
    del ne_1
    del octree_data_1

    # Octree 2
    points_2 = load(paths.OCTREE_DATA_PATH + "points_2.obj")
    print("Loaded points for octree 2 from", paths.OCTREE_DATA_PATH + "points_2.obj")

    ne_2 = load(paths.OCTREE_DATA_PATH + "ne_2.obj")
    print("Loaded ne values for octree 2 from", paths.OCTREE_DATA_PATH + "ne_2.obj")

    octree_data_2 = from_numpy_to_octree_data(points_2, ne_2)
    dump(paths.OCTREE_DATA_PATH + "octree_data_2.obj", octree_data_2)
    print("Stored data for octree 2 at", paths.OCTREE_DATA_PATH + "octree_data_2.obj")

    del points_2
    del ne_2
    del octree_data_2
    
    # Octree 3
    points_3 = load(paths.OCTREE_DATA_PATH + "points_3.obj")
    print("Loaded points for octree 3 from", paths.OCTREE_DATA_PATH + "points_3.obj")

    ne_3 = load(paths.OCTREE_DATA_PATH + "ne_3.obj")
    print("Loaded ne values for octree 3 from", paths.OCTREE_DATA_PATH + "ne_3.obj")

    octree_data_3 = from_numpy_to_octree_data(points_3, ne_3)
    dump(paths.OCTREE_DATA_PATH + "octree_data_3.obj", octree_data_3)
    print("Stored data for octree 3 at", paths.OCTREE_DATA_PATH + "octree_data_3.obj")

    del points_3
    del ne_3
    del octree_data_3
    
    # Octree 4
    points_4 = load(paths.OCTREE_DATA_PATH + "points_4.obj")
    print("Loaded points for octree 4 from", paths.OCTREE_DATA_PATH + "points_4.obj")

    ne_4 = load(paths.OCTREE_DATA_PATH + "ne_4.obj")
    print("Loaded ne values for octree 4 from", paths.OCTREE_DATA_PATH + "ne_4.obj")

    octree_data_4 = from_numpy_to_octree_data(points_4, ne_4)
    dump(paths.OCTREE_DATA_PATH + "octree_data_4.obj", octree_data_4)
    print("Stored data for octree 4 at", paths.OCTREE_DATA_PATH + "octree_data_4.obj")

    del points_4
    del ne_4
    del octree_data_4
    
