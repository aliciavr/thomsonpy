import argparse
import time
import numpy as np

from thomsonpy.data_management.data_model import DataModel
import thomsonpy.data_management.utils as utils
import thomsonpy.config.paths as paths




def arguments_set_up():
    # Parser initialization
    parser = argparse.ArgumentParser(prog="SPGCoronalModeler",
                                     description="It creates images modelling the Thomson scattering of the solar "
                                                 "corona.",
                                     epilog="Further information at https://aliciavr.github.io/thomsonpy")

    # Arguments creation
    parser.add_argument("-n", "--name", help="Name of the new model")
    parser.add_argument("-m", "--model-path", help="Path to the coronal model to be used.")
    parser.add_argument("echo")

    # Subparsers
    subparsers = parser.add_subparsers()
    # Model creator arguments
    model_creator_parser = subparsers.add_parser(name="Model creator",
                                                 description="Creates a pre-model for the solar image.")
    # Imager arguments
    imager_parser = subparsers.add_parser(name="Imager",
                                          description="Generates solar image.")

    # It reads the arguments from the command line
    arguments = parser.parse_args()

    return arguments


def create_model(ori_filepath, des_path, model_name, fragment_func, selection_func, format_func):
    """
    It creates a view of the raw data model ready to be used to create de image model and data
    structures that can accelerate the imaging process.
    This transformation implies filtering the original raw data to obtain a subspace to study, a
    reorganization of the points with the needed physical magnitudes in another structure and a
    normalization of the units to the I.S.

    :param ori_filepath: path of the original raw model.
    :type ori_filepath: string
    :param des_path: path where the new model will be stored.
    :type des_path: string
    :param model_name:
    :type model_name: string
    :param fragment_func: function to extract a subspace of the original raw model.
    :type fragment_func: `function`
    :param selection_func: filter for selecting a subspace of the original raw model.
    :type selection_func: `function`
    :param format_func: function to apply a common format.
    :type format_func: `function`
    """

    # Creating a preliminary view of the raw model.
    pre_data = fragment_func(   ori_filepath=ori_filepath,
                                selection_func=selection_func,
                                format_func=format_func)

    data_model = DataModel(pre_data)

    # Storing the data model.
    utils.store_object(obj=data_model, filename=model_name,
                       extension=paths.DATA_MODELS_EXTENSION, path=des_path)

    return data_model


"""

def create_data_structure(name_dstruct, data_model):
    coordinates = data_model.get_all_coordinates()
    dstruct = None
    if name_dstruct == "octree":
       dstruct = Octree(level_limit=,
                        data_limit=,
                        coordinates=coordinates,
                        min_v=,
                        max_v=)
    elif name_dstruct == "kdtree":
       dstruct = KDTree(coordinates=coordinates)
    else:
        raise Exception("Non valid data structure")
    return dstruct

def compute_thomson_scattering_image(num_points, min_coord, max_coord, data_model, data_structure):
    resolution = (max_coord - min_coord) / num_points
    print("Imaging between", (min_coord, min_coord) * units.METERS_TO_RSOL, "RSol and", (max_coord, max_coord) * units.METERS_TO_RSOL,
          "RSol.")
    print("Resolution = ", resolution / 1000, " km (", num_points, "x", num_points, ").")
    x_values = np.linspace(min_coord, max_coord, num_points)  # de - a +
    # print(x_values * units.METERS_TO_RSOL)
    y_values = np.linspace(min_coord, max_coord, num_points)[::-1]  # de + a -
    # print(y_values * units.METERS_TO_RSOL)
    image_model = np.zeros((num_points, num_points))
    print("# Numeric integral steps =", tsp.NUM_Z)

    num_points = num_points ** 2
    print(f"Num points = {num_points}")
    count = 0
    ini_time = time.perf_counter()
    for y in range(min_coord, max_coord):
        for x in range(min_coord, max_coord):
            if (x_values[x] ** 2 + y_values[y] ** 2) > tsp.SOLAR_RADIUS ** 2:
                # Coordinates with the center of the Sun as Origin of the Reference System.
                target = (x_values[x], y_values[y], 0)
                # Creates a ThomsonGeometry object to manage the ray - tracing across the Corona.
                TG = thtools.ThomsonGeometry(sip.SUN_CENTER, sip.OBSERVER, target, tsp.SOLAR_RADIUS)

                # Line of sight integration generating a value for the scattered light model.
                scattered_light = thtools.get_scattered_light(tsp.WAVELENGTH, tsp.T_SOL, tsp.X,
                                                              TG.get_elongation(), tsp.INI_Z, tsp.FIN_Z,
                                                              tsp.INCR_Z, TG, NE_MODEL)
                model[y][x] = scattered_light
            else:
                model[y][x] = 0

            # Progress
            count += 1
            if count % 10000 == 0:
                print("Progress", count / num_points * 100, "%")

    fin_time = time.perf_counter()
    print("Model built in", str((fin_time - ini_time) / 60 / 60), "hours.")
    return image_model

"""

if __name__ == "__main__":

    # It manages the console input parameters.
    args = arguments_set_up()
    print(args.name, args.model_path, args.echo)

    model_name = "default"
    model = "default"
    echo = "default"

    if args.name:
        model_name = args.name
    if args.model_path:
        model = args.model_path
    if args.echo:
        echo = args.echo
    print(f"my name is {model_name}")
    print(f"my model is in {model}")

    # It creates the data model.
    data_model = None
    if args.data_model:
        data_model = utils.load_object(f"{paths.MODELS_PATH}{args.name}.{paths.DATA_MODELS_EXTENSION}")
    else:
        data_model = create_model(ori_path=paths.PREDSCI_FILENAME,
                                  dest_path=paths.MODELS_PATH,
                                  model_name=args.name,
                                  fragment_func=utils.predsci_fragmentation,
                                  selection_func=utils.predsci_selection,
                                  format_func=utils.apply_data_format_to_predsci)

    """
    # It creates the data structure.
    data_structure = None
    if args.data_structure:
        data_structure = utils.load_object(f"{dest_path}{data_struct_name}.dstruct")
    else:
        data_structure = create_data_structure(name_dstruct, data_model)

    # It creates the solar corona model to be imaged.
    if data_model & data_structure:
        compute_thomson_scattering_image(data_model, data_structure)
        utils.store_object(f"{num_points}p_{image_model}", model_name, "model", models_path)
    """