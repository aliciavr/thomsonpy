# -*- coding: utf-8 -*-
"""
.. module:: octree
        :platform: Unix
        :synopsis: it manages an octree and its functionalities.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

import numpy as np
import open3d as o3d
import sys

from thomsonpy.data_management.data_structures.data_structure import DataStructure
from thomsonpy.data_management.data_model import Data
from thomsonpy.main import PHYSICAL_MODEL


class Node:
    """
    Node coordinates structure for constructing the ``Octree`` (:py:class:`thomsonpy.data_management.octree.octree.Octree`).
    It contains points as ``Data`` (:py:class:`thomsonpy.data_management.octree.octree.Data`).
    """

    def __init__(self, level, ant, min_v, max_v, octree):
        """
        Constructor of ``Node`` class.
        
        :param level: depth level where the node is created in the ``Octree`` coordinates structure.
        :type level: int
        :param ant: reference to the ``Node`` in the previous depth level which spatially contains this one.
        :type ant: int
        :param min: minimum coordinates of the spatial region enclosed by this ``Node``.
        :type min: numpy.array[(float), (float), (float)]
        :param max: maximum coordinates of the spatial region enclosed by this ``Node``.
        :type max: numpy.array[(float), (float), (float)]
        :param octree: reference to the ``Octree`` containing this node in its internal structure.
        :type octree: :py:class:`thomsonpy.data_management.octree.octree.Octree`
        
        :return: an empty octree node.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Node`
        """
        self.__predecessor = ant
        self.__level = level
        self.__children = list()
        self.__octree = octree
        self.__min = min_v
        self.__max = max_v
        self.__indices = list()

    def __str__(self):
        return 'Node at level {} with {} coordinates at coordinates {} and {}. Total coordinates:{}.'\
            .format(self.__level, len(self.__indices), self.__min, self.__max, len(self.__indices))

    def __create_children(self):
        """
        It creates the eight children of this ``Node``.
        """
        self.__children = [None] * 8

        min_x = self.__min[0]
        min_y = self.__min[1]
        min_z = self.__min[2]
        max_x = self.__max[0]
        max_y = self.__max[1]
        max_z = self.__max[2]
        mid_x = (max_x + min_x) / 2
        mid_y = (max_y + min_y) / 2
        mid_z = (max_z + min_z) / 2

        self.__children[0] = Node(self.__level + 1, self, np.array([min_x, min_y, min_z]),
                                  np.array([mid_x, mid_y, mid_z]), self.__octree)

        self.__children[1] = Node(self.__level + 1, self, np.array([mid_x, min_y, min_z]),
                                  np.array([max_x, mid_y, mid_z]), self.__octree)

        self.__children[2] = Node(self.__level + 1, self, np.array([mid_x, mid_y, min_z]),
                                  np.array([max_x, max_y, mid_z]), self.__octree)

        self.__children[3] = Node(self.__level + 1, self, np.array([min_x, mid_y, min_z]),
                                  np.array([mid_x, max_y, mid_z]), self.__octree)

        self.__children[4] = Node(self.__level + 1, self, np.array([min_x, min_y, mid_z]),
                                  np.array([mid_x, mid_y, max_z]), self.__octree)

        self.__children[5] = Node(self.__level + 1, self, np.array([mid_x, min_y, mid_z]),
                                  np.array([max_x, mid_y, max_z]), self.__octree)

        self.__children[6] = Node(self.__level + 1, self, np.array([mid_x, mid_y, mid_z]),
                                  np.array([max_x, max_y, max_z]), self.__octree)

        self.__children[7] = Node(self.__level + 1, self, np.array([min_x, mid_y, mid_z]),
                                  np.array([mid_x, max_y, max_z]), self.__octree)

    def add(self, coordinates):
        """
        It adds a new ``Data`` to the ``Octree``.
        
        :param coordinates: new spatial ``Data`` to be added to the ``Octree``.
        :type coordinates: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        index = PHYSICAL_MODEL.get_index(coordinates)  # Index identifying the coordinates passed as a parameter.
        if not self.__children:  # If it does not have children.
            if self.__level == self.__octree.get_level_limit() or len(self.__indices) < self.__octree.get_data_limit():
                # If it has reached the maximum level or there is free space in this node it inserts the index
                # associated with the coordinates passed as a parameter.
                self.__indices.append(index)
            else:
                # If it has not reached the maximum level and there is no free space in this node it creates children
                # nodes.
                self.__indices.append(index)
                self.__create_children()
                for i in self.__indices:
                    for c in self.__children:
                        if c.contains(PHYSICAL_MODEL.get_coordinates(i)):
                            c.add(i)
                self.__indices.clear()
        else:  # If it has children.
            for c in self.__children:
                # It searches in each node containing the coordinates passed as a parameter.
                if c.contains(coordinates):
                    c.add(coordinates)

    def search(self, coordinates):
        """
        It searches for a specific coordinates previously stored in the ``Octree`` coordinates structure.
        
        :param coordinates: coordinates of the coordinates to be found with information about physical measures of the solar
        corona.
        :type coordinates: numpy.ndarray([float, float, float])
        
        :return: a ``Data`` obj representing the result of the searching process. If the search is successful,
        it returns the ``Data`` obj related with the 3D coordinates ``coordinates`` passed as parameter, otherwise it returns a
        ``None`` obj. :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        if self.contains(coordinates):
            # If this ``Node`` contains the coordinates passed as a parameter.
            if self.__children:
                # If it has children.
                for c in self.__children:
                    # It searches which child contains the coordinates passed as a parameter.
                    if c.contains(coordinates):
                        return c.search(coordinates)
            else:
                # It is a leaf node.
                return self
        return None

    def search_nearest_data_index(self, coordinates):
        """
        It searches for the nearest coordinates previously stored in the ``Octree`` coordinates structure.
        
        :param coordinates: coordinates to be located for retrieving the nearest coordinates information physical coordinates
        of the solar corona.
        :type coordinates: numpy.ndarray([float, float, float])

        :return: -1 if the coordinates are not located inside this ``Octree``, the index associated by the
        ``PHYSICAL_MODEL`` to the nearest coordinates to the ones passed as a parameter.
        :rtype: int
        """
        if self.contains(coordinates):
            # If this ``Node`` contains the coordinates passed as a parameter.
            if self.__children:
                # If it has children.
                for c in self.__children:
                    # It searches which child contains the coordinates passed as a parameter.
                    if c.contains(coordinates):
                        return c.search_nearest(coordinates)
            else:
                # If it is a leaf node it searches for the nearest coordinates.
                nearest = -1
                min_distance = sys.float_info.max
                for i in self.__indices:
                    i_coordinates = PHYSICAL_MODEL.get_coordinates(i)
                    distance = np.linalg.norm(coordinates - i_coordinates)
                    if distance < min_distance:
                        nearest = i
                        min_distance = distance
                return nearest
        return -1

    def contains(self, coordinates):
        """
        It checks if this ``Node`` is containing the coordinates ``coordinates`` passed as parameter.

        :param coordinates: coordinates in :math:`\mathbb{R}^3`.
        :type coordinates: numpy.ndarray([float, float, float])

        :return: ``True`` if ``coordinates`` is inside this ``Node``, ``False`` otherwise.
        :rtype: boolean
        """
        x = coordinates[0]
        y = coordinates[1]
        z = coordinates[2]
        return x <= self.__max[0] and y <= self.__max[1] and z <= self.__max[2] and x >= self.__min[0] and y >= \
               self.__min[1] and z >= self.__min[2]

    def get_min(self):
        """
        It gets the minimum coordinates of the spatial region enclosed by this ``Node``.
        
        :return: the minimum coordinates of the spatial region enclosed by this ``Node``.  
        :rtype: numpy.ndarray([float, float, float])
        """
        return self.__min

    def get_max(self):
        """
        It gets the maximum coordinates of the spatial region enclosed by this ``Node``.
        
        :return: the maximum coordinates of the spatial region enclosed by this ``Node``.
        :rtype: numpy.ndarray([float, float, float])
        """
        return self.__max

    def get_data(self):
        """
        It gets the ``list`` of ``Data`` objects inside this ``Node``.
        
        :return: the ``Data`` objects inside this 
        :rtype: Data
        """
        return self.__indices

    def get_num_points(self):
        """
        It gets the number of points inside this ``Node``
        
        :return: number of points inside this ``Node``.
        :rtype: float
        """
        return len(self.__indices)

    def has_children(self):
        """
        It checks if this ``Node`` has children nodes.
        
        :return: ``True`` if it has children nodes created, ``False`` otherwise. 
        :rtype: boolean
        """
        return len(self.__children) != 0

    def has_data(self):
        """
        It checks if this ``Node`` has any coordinates stored.
        
        :return: ``True`` if it has coordinates stored, ``False`` otherwise.
        :rtype: boolean
        """
        return len(self.__indices) != 0

    def get_level(self):
        """
        It gets the depth level where this ``Node`` has been created in the internal ``Octree`` structure.
        
        :return: the depth level of this ``Node``.
        :rtype: int
        """
        return self.__level

    def get_representation(self):
        """
        It generates geometric structures for the ``open3D`` library for rendering the octree. tangential_intensity adds the geometric
        structures to the attribute ``visual_nodes`` of the ``Octree`` coordinates structure.
        """
        min_x = self.__min[0]
        min_y = self.__min[1]
        min_z = self.__min[2]
        max_x = self.__max[0]
        max_y = self.__max[1]
        max_z = self.__max[2]
        points = [
            [min_x, min_y, min_z],
            [max_x, min_y, min_z],
            [min_x, max_y, min_z],
            [max_x, max_y, min_z],
            [min_x, min_y, max_z],
            [max_x, min_y, max_z],
            [min_x, max_y, max_z],
            [max_x, max_y, max_z],
        ]
        lines = [
            [0, 1],
            [0, 2],
            [1, 3],
            [2, 3],
            [4, 5],
            [4, 6],
            [5, 7],
            [6, 7],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
        ]
        line_set = o3d.geometry.LineSet(
            points=o3d.utility.Vector3dVector(points),
            lines=o3d.utility.Vector2iVector(lines),
        )

        if self.__level >= 0:
            # if (self.has_data()):
            self.__octree.visual_nodes.append(line_set)
        for c in self.__children:
            c.see_node()


class Octree(DataStructure):
    """
    It manages the internal and root structure of an octree coordinates structure.
    """

    def __init__(self, level_limit, data_limit, coordinates, min_v, max_v):
        """
        Constructor of the ``Octree`` data structure implemented manually in this library adapting its structure to the
        requirements of the computations.

        :param level_limit: the maximum level of tree depth to be allowed.
        :type level_limit: int
        :param data_limit: the maximum number of indices to be stored by internal nodes.
        :type data_limit: int
        :param coordinates: list of coordinates to structured in this ``Octree``.
        :type coordinates: list(numpy.ndarray([float, float, float]))
        :param min_v: vector representing the lower bounds of the 3D space represented by this ``Octree``.
        :type min_v: numpy.ndarray([float, float, float])
        :param max_v: vector representing the upper bounds of the 3D space represented by this ``Octree``.
        :type max_v: numpy.ndarray([float, float, float)]

        :return: an ``Octree`` with the specified parameters.
        :rtype: :py:class:`thomsonpy.data_management.octree.Octree`
        """
        self.__min = min_v
        self.__max = max_v
        self.__level_limit = level_limit
        self.__data_limit = data_limit
        self.__num_data = len(coordinates)
        self.__root = Node(0, None, min_v, max_v, self)
        progress = 0
        total = len(coordinates)
        for c in coordinates:
            self.__root.add(c)
            # progress...
            if progress % 100000 == 0:
                print(progress / total * 100, "%")
            progress += 1
        self.visual_nodes = list()

    def search_nearest_data_index(self, coordinates):
        """
        It calls the ``search_nearest_data_index`` method from the root ``Node`` of this ``Octree``.

        .. seealso:: The method :py:meth:`thomsonpy.data_management.octree.octree.Node.search_nearest` using the
        recursion to search for the ``coordinates`` passed as parameter.

        :return: -1 if the coordinates are not located inside this ``Octree``, the index associated by the
        ``PHYSICAL_MODEL`` to the nearest coordinates to the ones passed as a parameter.
        :rtype: int
        """
        return self.__root.search_nearest_data_index(coordinates)

    def search_nearest(self, coordinates):
        """
        It calls the ``search_nearest`` method from the root ``Node`` of this ``Octree``.
        It searches for the nearest coordinates previously stored in the ``Octree`` coordinates structure.

        .. seealso:: The method :py:meth:`thomsonpy.data_management.octree.octree.Node.search_nearest` using the
        recursion to search for the ``coordinates`` passed as parameter.

        :param coordinates: coordinates to be located for retrieving the nearest coordinates information physical coordinates
        of the solar corona.
        :type coordinates: numpy.ndarray([float, float, float])

        :return: a ``Data`` obj representing the result of the nearest searching process. If the search is
        successful, it returns the nearest ``Data`` obj to the 3D coordinates ``coordinates`` passed as parameter,
        otherwise it returns a ``Data`` with the default values.
        :rtype: :py:class:`thomsonpy.data_management.data_model.Data`
        """

        nearest_index = self.__root.search_nearest_data_index(coordinates)
        if nearest_index == -1:
            data = Data(None, 0)
        else:
            data = PHYSICAL_MODEL.get_data_from_index(nearest_index)
        return data

    def search(self, coordinates):
        """
        It calls the ``search`` method from the root ``Node`` of this ``Octree``.

        .. seealso:: The method :py:meth:`thomsonpy.data_management.octree.octree.Node.search` using the recursion
        to search for the ``coordinates`` passed as parameter.

        :return: a ``Data`` obj representing the result of the searching process. If the search is successful,
        it returns the ``Data`` obj related with the 3D coordinates ``coordinates`` passed as parameter, otherwise it returns a
        ``None`` obj. :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        return self.__root.search(coordinates)

    def get_representation(self):
        """
        It gets a ``list`` with all the geometric objects representing the internal structure of this ``Octree``.

        :return: a list with geometric objects readable by the ``Open3D`` library.
        :rtype: ``list``
        """
        self.__root.get_representation()
        return self.visual_nodes

    def get_level_limit(self):
        """
        It gets the level limit used as a restriction for the construction process. If this level is reached the
        ``Node`` objects at this level do not create more children nodes anymore.
        
        :return: the level limit.
        :rtype: int
        """
        return self.__level_limit

    def get_data_limit(self):
        """
        It gets the coordinates limit used as a restriction for the construction process. If this quantity of ``Data``
        objects in a particular ``Node`` is reached and a new ``Data`` obj is inserted in that ``Node``,
        then that ``Node`` creates eight children if ``level_limit`` has not been reached yet.
        
        :return: the limit of coordinates allowed by node when constructing the octree.
        :rtype: int
        """
        return self.__data_limit

    def get_num_data(self):
        """
        It gets the number of ``Data`` objects stored in this ``Octree``.
        
        :return: the number of coordinates stored in the octree.
        :rtype: int
        """
        return self.__num_data
