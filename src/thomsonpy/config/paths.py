"""
.. module:: paths
        :platform: Unix
        :synopsis: package relative paths for loading and storing models.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

PROJECT_NAME = "first-quadrant"
""" Path to the base folder of the project in execution. """

PREDSCI_FILENAME = "rho002.hdf"
""" Filename of the raw model. """

PREDSCI_DATA_PATH = "data/predictive_science/eclipse2021_mhd_final/corona/"
""" Path to the base folder of the raw models provided by Predictive Science Inc. """

OCTREE_DATA_PATH = f"data/projects/{PROJECT_NAME}/octree_data_format/"
""" Path to the octree formated data. """

OCTREES_PATH = f"data/projects/{PROJECT_NAME}/octrees/"
""" Path to the folder containing the binary objects with the octrees. """

MODELS_PATH = f"data/projects/{PROJECT_NAME}/models/"
""" Path to the folder containing the final image models of the Thomson scattering computation. """