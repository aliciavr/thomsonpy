"""
.. module:: paths
        :platform: Unix
        :synopsis: package relative paths for loading and storing models.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""
PREFIX = "../../"
RAW_MODELS_EXTENSION = "hdf"
DATA_MODELS_EXTENSION = "model"
DATA_STRUCTURES_EXTENSION = "dstruct"
IMAGE_MODELS_EXTENSION = "imodel"

PROJECT_NAME = "cmag"
""" Path to the base folder of the project in execution. """

RAW_MODELS_PATH = f"{PREFIX}data/input/eclipse2021_mhd_final/corona/"
""" Path to the base folder of the raw models provided by Predictive Science Inc. """

RAW_FILENAME = [f"rho002.{RAW_MODELS_EXTENSION}", f"t002.{RAW_MODELS_EXTENSION}"]
""" Filename of the raw model. """

DATA_MODELS_PATH = f"{PREFIX}data/output/{PROJECT_NAME}/data_models/"
""" Path to the formatted data. """

STRUCTURES_PATH = f"{PREFIX}data/output/{PROJECT_NAME}/structures/"
""" Path to the folder containing the binary objects with the data structures. """

IMODELS_PATH = f"{PREFIX}data/output/{PROJECT_NAME}/imodels/"
""" Path to the folder containing the final image models of the Thomson scattering computation. """

