import pickle
import numpy as np
from pyhdf.SD import SD
import thomsonpy.constants.units as units
from thomsonpy.data_management.data_model import Data


def store_object(obj, filename, path):
    """
    Auxiliary function for storing coordinates.

    :param obj: object to be stored.
    :type obj: any
    :param filename: name for the file storing the object.
    :type filename: string
    :param path: path of the folder where the object will be stored.
    :type path: string
    """
    f = open(f"{path}{filename}.dstruct", 'wb')
    pickle.dump(obj, f)
    f.close()


def load_object(filepath):
    """
    Auxiliary function for loading stored coordinates.

    :param filepath: name for the file storing the coordinates.
    :type filepath: string

    :return: the object stored in filepath passed as parameter.
    :rtype: any
    """
    f = open(filepath, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj


def cartesian_to_spherical(coordinates):
    x = coordinates[0]
    y = coordinates[1]
    z = coordinates[2]
    radial = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = np.arctan(x / z)
    theta = np.arccos(y / radial)
    return np.array([phi, theta, radial])


def spherical_to_cartesian(r, theta, phi):
    """
    It gets cartesian coordinates from spherical coordinates.

    .. math::
        x &= \sin{\\theta} \sin{\phi}

        y &= \cos{\\theta}

        z &= \sin{\\theta} \cos{\phi}

    :param r: radius of the spherical coordinates
    :type r: float
    :param theta: latitude, starting from the North [0 - PI]
    :type theta: float
    :param phi: longitude: Carrington Longitude [0 - 2PI]
    :type phi: float

    :return: coordinates in the Cartesian Coordinate System
    :rtype: numpy.ndarray[float, float, float]
    """
    x = r * np.sin(theta) * np.sin(phi)
    y = r * np.cos(theta)
    z = r * np.sin(theta) * np.cos(phi)
    return np.array([x, y, z])


def apply_data_format(r, theta, phi, ne_mas, t_mas):
    """
    It gives the required coordinates format by the octree structure.

    The original coordinates is stored in hdf format in the Spherical Coordinates System
    in mas units. The octree structure needs the coordinates stored as a list of Data
    objects where each one stores its coordinates in the Cartesian Coordinates
    System and its physical magnitudes associated in the International System of
    Units.

    :param r: radius coordinate of the Spherical Coordinates System in RSol
    :type r: float
    :param theta: theta coordinate of the Spherical Coordinates System in rad
    :type theta: float
    :param phi: phi coordinate of the Spherical Coordinates System in rad
    :type phi: float
    :param ne_mas: ne value for the given coordinates in MAS units.
    :type ne_mas: float

    :return: a Data object with coordinates in the Cartesian Coordinates System and
            with units in the International System of Units.
    :rtype: :class:`thomsonpy.data_management.octree.octree.Data`

    """
    coordinates = spherical_to_cartesian(r, theta, phi) * units.RSOL_TO_METERS  # From RSol to m
    ne = ne_mas * units.NE_MAS_FACTOR  # From MAS to m⁻³.
    temperature = t_mas * units.T_MAS_FACTOR # From MAS to K.
    data = Data(coordinates, (ne, temperature))
    return data


def get_predsci_raw(filepath):
    """
    It loads a raw model of the solar corona provided by Predictive Science Inc.

    :param filepath:
    :type filepath:

    :return: a raw model of one of the physical magnitudes of the solar corona.
    :rtype: np.ndarray(3D)
    """
    hdf = SD(filepath)
    print("\nLoading datasets from", filepath, ":\n", hdf.datasets())
    data = hdf.select(3).get()  # DATA CUBE

    num_points = data.size
    print("# POINTS = ", num_points)

    return data


def get_raw_coords(filepath, opt):
    """
    It gets the values for the coordinates indicated at ``opt`` parameter of the raw coordinates cube provided by Predictive
    Science and that can be obtained with the function :py:func:`thomsonpy.data_management.fragmenter`.

    :param filepath: filepath to the file storing the raw model.
    :type filepath: string :param opt: it accepts the values ``phi``, ``theta`` or ``radial``, indicating the coordinate
    in the Spherical Coordinates Systems to be loaded.
    :type opt: string

    :return: the values for the chosen spherical coordinate. ``None`` if the ``opt`` parameter has an invalid value.
    :rtype: np.array([float])
    """

    hdf = SD(filepath)
    print("\nLoading datasets from", filepath, ":\n", hdf.datasets())

    coords = None
    if opt == "phi":
        coords = hdf.select(0).get()  # PHI (rad, 0-2PI) 699
        num_phi = coords.size
        print("\nPHI (rad, 0-2PI). # phi =", num_phi)
    elif opt == "theta":
        coords = hdf.select(1).get()  # THETA (rad, 0-PI) 327
        num_theta = coords.size
        print("THETA (rad, 0-PI). # theta =", num_theta)
    elif opt == "radial":
        coords = hdf.select(2).get()  # RADIAL (RSOL, 1-30) 288
        num_radial = coords.size
        print("RADIAL (RSOL, 1-30). # radial =", num_radial)

    return coords


def selection(r, theta, phi, ne, temperature):
    """

    :param r: radial coordinate in the Spherical System of Coordinates.
    :type r: float
    :param theta: theta coordinate in the Spherical System of Coordinates.
    :type theta: float
    :param phi: phi coordinate in the Spherical System of Coordinates.
    :type phi: float
    :param ne: electron density model.
    :type ne: float
    :param temperature: temperature model.
    :type temperature: float

    :return:
    :rtype: boolean

    """
    cartesian_coords = spherical_to_cartesian(r, theta, phi)
    x = cartesian_coords[0]
    y = cartesian_coords[1]
    z = cartesian_coords[2]
    limit = 3.5
    return x <= limit and x >= -limit and y <= limit and y >= -limit and z <= limit and z >= -limit


def fragment(selection_func, format_func, raw_model, data_radial, data_theta, data_phi, progress_step=1e6):
    """

    :param selection_func: function defining the subspace of points.
    :type selection_func: function
    :param format_func: function defining the format of the output of this fragmentation.
    :type format_func: function
    :param raw_model: original electron density model from Predictive Science Inc.
    :type raw_model: 3d array
    :param data_radial: values for the radial coordinate.
    :type data_radial: array of float
    :param data_theta: values for the theta coordinate.
    :type data_theta: array of float
    :param data_phi: values for the phi coordinate.
    :type data_phi: array of float
    :param progress_step: int
    :type progress_step: int
    """

    it = np.nditer(raw_model, flags=['multi_index'])
    num_points = raw_model.size
    print(num_points)
    data = []
    # Progress and auxiliary variables
    progress = 0
    # Fragmentation process.
    for ne_mas in it:
        i, j, k = it.multi_index
        r = data_radial[k]
        theta = data_theta[j]
        phi = data_phi[i]

        if selection_func(r, theta, phi, ne_mas):
            d = format_func(r=r, theta=theta, phi=phi, ne_mas=ne_mas)
            data.append(d)
        if progress % progress_step == 0:
            print(progress / num_points * 100, "%")
        progress += 1

    return data
