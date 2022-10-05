import pickle
import numpy as np
from pyhdf.SD import SD
import thomsonpy.constants.units as units
from thomsonpy.data_management.data_model import Data


def store_object(obj, filename, extension, path):
    """
    Auxiliary function for storing coordinates.

    :param obj: object to be stored.
    :type obj: any
    :param filename: name for the file storing the object.
    :type filename: string
    :param extension: extension for the file storing the object.
    :type extension: string
    :param path: path of the folder where the object will be stored.
    :type path: string
    """
    f = open(f"{path}{filename}.{extension}", 'wb')
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


def apply_data_format_to_predsci(r, theta, phi, mas_physical_magnitudes):
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
    :param mas_physical_magnitudes: tuple with the physical magnitudes in mas.
    :type mas_physical_magnitudes: tuple

    :return: a Data object with coordinates in the Cartesian Coordinates System and
            with units in the International System of Units.
    :rtype: :class:`thomsonpy.data_management.octree.octree.Data`

    """

    # It transforms the coordinates into the Cartesian System and in I.S. of units.
    coordinates = spherical_to_cartesian(r, theta, phi) * units.RSOL_TO_METERS  # From RSol to m

    # It transforms the physical magnitudes into the I.S. of units.
    r = range(0, len(mas_physical_magnitudes))
    conv_factors = [units.NE_MAS_FACTOR, units.T_MAS_FACTOR]
    physical_magnitudes = list()
    for m, c in zip(mas_physical_magnitudes, conv_factors):
        physical_magnitudes.append(m * c)

    data = Data(coordinates, physical_magnitudes)
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


def predsci_selection(r, theta, phi):
    """

    :param r: radial coordinate in the Spherical System of Coordinates.
    :type r: float
    :param theta: theta coordinate in the Spherical System of Coordinates.
    :type theta: float
    :param phi: phi coordinate in the Spherical System of Coordinates.
    :type phi: float

    :return:
    :rtype: boolean

    """
    cartesian_coords = spherical_to_cartesian(r, theta, phi)
    x = cartesian_coords[0]
    y = cartesian_coords[1]
    z = cartesian_coords[2]
    limit = 3.5
    return x <= limit and x >= -limit and y <= limit and y >= -limit and z <= limit and z >= -limit


def predsci_fragmentation(ori_filepath, selection_func, format_func, progress_step=1e6):
    """
    It creates a preliminary model view of an original Predictive Science raw model with a common
    and useful structure for further computations.

    :param ori_filepath: path of the original Predictive Science raw model.
    :param selection_func: function defining the subspace of points.
    :type selection_func: function
    :param format_func: function defining the format of the output of this fragmentation.
    :type format_func: function
    :param progress_step: int
    :type progress_step: int
    """

    raw_models = list()
    for f in ori_filepath:
        raw_model = get_predsci_raw(f)
        raw_models.append(raw_model)
        print("Size of the raw model", raw_model.size)

    if raw_models:
        data_radial = get_raw_coords(ori_filepath[0], "radial")
        data_theta = get_raw_coords(ori_filepath[0], "theta")
        data_phi = get_raw_coords(ori_filepath[0], "phi")

        iterators = list()
        for raw in raw_models:
            it = np.nditer(raw, flags=["multi_index"])
            iterators.append(it)
            num_points = raw_models[0].size

        data = list()
        # Progress and auxiliary variables
        progress = 0
        # Fragmentation process.
        for p in zip(*iterators):
            i, j, k = iterators[0].multi_index
            r = data_radial[k]
            theta = data_theta[j]
            phi = data_phi[i]
            if selection_func(r=r, theta=theta, phi=phi):
                d = format_func(r=r, theta=theta, phi=phi, mas_physical_magnitudes=p)
                data.append(d)
            if progress % progress_step == 0:
                print(progress / num_points * 100, "%")
            progress += 1

    return data
