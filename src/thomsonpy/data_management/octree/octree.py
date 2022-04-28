# -*- coding: utf-8 -*-
"""
.. module:: octree
        :platform: Unix
        :synopsis: it manages an octree and its functionalities.
.. moduleauthor:: Alicia Vázquez Ramos (SPG - IAA) <aliciavr@iaa.es>
"""

import numpy as np
import open3d as o3d
import pickle
import sys

class Data:
    """
    Class representing structured object for the Octree class :py:class:`thomsonpy.data_management.octree.octree.Octree`.
    It relates the electron density with a point in :math:`\mathbb{R}^3`.
    """
    def __init__(self, coordinates, ne):
        """
        Constructor of ``Data`` class.
        
        :param coordinates: coordinates of the data in :math:`\mathbb{R}^3`.
        :type coordinates: numpy.array([float, float, float])
        :param ne: electron density at the point defined by the ``coordinates`` parameter.
        :type ne: float
        
        :return: a data with spatial and electron density information.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        self.__coordinates = coordinates
        self.__ne = ne

    def __str__(self):
        return 'Data = {} at ({}, {}, {})'.format(self.__ne, self.__coordinates[0], self.__coordinates[1], self.__coordinates[2])

    def get_coordinates(self):
        """
        It gets the coordinates of the point represented by this ``Data`` object.
        
        :return: the coordinates of this ``Data`` object.
        :rtype: numpy.array([float, float, float])
        """
        return self.__coordinates

    def get_x(self):
        """
        It gets the :math:`x` coordinate of the point represented by this ``Data`` object.
        
        :return: the coordinate :math:`x` of this ``Data`` object.
        :rtype: float
        """
        return self.__coordinates[0]

    def get_y(self):
        """
        It gets the :math:`y` coordinate of the point represented by this ``Data`` object.
        
        :return: the coordinate :math:`y` of this ``Data`` object.
        :rtype: float
        """
        return self.__coordinates[1]

    def get_z(self):
        """
        It gets the :math:`z` coordinate of the point represented by this ``Data`` object.
        
        :return: the coordinate :math:`z` of this ``Data`` object.
        :rtype: float
        """
        return self.__coordinates[2]

    def get_ne(self):
        """
        It gets the electron density :math:`\\rho` related with the point represented by this ``Data`` object.
        
        :return: the electron density :math:`\\rho` of this ``Data`` object.
        :rtype: float
        """
        return self.__ne

class Node:
    """
    Node data structure for constructing the ``Octree`` (:py:class:`thomsonpy.data_management.octree.octree.Octree`). It contains points as ``Data`` (:py:class:`thomsonpy.data_management.octree.octree.Data`).
    """

    def __init__(self, level, ant, min, max, octree):
        """
        Constructor of ``Node`` class.
        
        :param level: depth level where the node is created in the ``Octree`` data structure.
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
        self.__antecessor = ant
        self.__level = level
        self.__children = list()
        self.__octree = octree
        self.__min = min
        self.__max = max
        self.__data = list()

    def __str__(self):
        return 'Node at level {} with {} data at coordinates {} and {}. Total data:{}.'.format(self.__level, len(self.__data), self.__min, self.__max, len(self.__data))

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
    
    def add(self, data):
        """
        It adds a new ``Data`` to the ``Octree``. 
        
        :param data: new spatial ``Data`` to be added to the ``Octree``.
        :type data: :py:class:`thomsonpy.data_mangement.octree.octree.Data`
        """
        if len(self.__children) == 0: # If it does not has children
            if self.__level == self.__octree.get_level_limit() or len(self.__data) < self.__octree.get_data_limit():
                # If it has reached the maximum level or there is free space in this node it inserts data
                self.__data.append(data)
            else:
                # If it has not reached the maximum level and there is not free space in this node
                self.__data.append(data)
                self.__create_children()
                for p in self.__data:
                    for c in self.__children:
                        if c.contains(p.get_coordinates()):
                            c.add(p)
                self.__data.clear()
        else: # If it has children
            for c in self.__children:
                # Searchs in each node containing this data
                if c.contains(data.get_coordinates()):
                    c.add(data)

    def search(self, p):
        """
        It searchs for a specific point previously stored in the ``Octree`` data structure.
        
        :param p: point to be found with information about electron density :math:`\\rho`.
        :type p: numpy.array([float], [float], [float])
        
        :return: a ``Data`` object representing the result of the searching process. If the search is successful, it returns the ``Data`` object related with the 3D point ``p`` passed as parameter, otherwise it returns a ``None`` object.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        if self.contains(p):
            if len(self.__children) != 0:
                # Has children and contains the point
                for c in self.__children:
                    # Searchs the point in each children which contains the point
                    if c.contains(p):
                        return c.search(p)
            else:
                # Is a leaf node
                return self
        return None
    
    def search_nearest(self, p):
        """
        It searchs for a specific point previously stored in the ``Octree`` data structure.
        
        :param p: point to be located for retrieving the nearest data information about electron density :math:`\\rho`.
        :type p: numpy.array([float], [float], [float])
        
        :return: a ``Data`` object representing the result of the nearest searching process. If the search is successful, it returns the nearest ``Data`` object to the 3D point ``p`` passed as parameter, otherwise it returns a ``None`` object.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        if self.contains(p):
            # If contains the point
            if len(self.__children) != 0:
                # If has children
                for c in self.__children:
                    if c.contains(p):
                        return c.search_nearest(p)
            else:
                # Is a leave node
                nearest = Data(None, 0)
                min_distance = sys.float_info.max
                for data in self.__data:
                    distance = np.linalg.norm(p - data.get_coordinates())
                    if distance < min_distance:
                        nearest = data
                        min_distance = distance
                return nearest
        return Data(None, 0)

    def contains(self, p):
        """
        It checks if this ``Node`` is containing the point ``p`` passed as parameter.
        
        :param p: a point in :math:`\mathbb{R}^3`.
        :type p: numpy.array([float], [float], [float])
        
        :return: ``True`` if ``p`` is inside this ``Node``, ``False`` otherwise.
        :rtype: boolean
        """
        x = p[0]
        y = p[1]
        z = p[2]        
        return x <= self.__max[0] and y <= self.__max[1] and z <= self.__max[2] and x >= self.__min[0] and y >= self.__min[1] and z >= self.__min[2]

    def get_min(self):
        """
        It gets the minimum coordinates of the spatial region enclosed by this ``Node``.
        
        :return: the minimum coordinates of the spatial region enclosed by this ``Node``.  
        :rtype: numpy.array([float][float][float])
        """
        return self.__min

    def get_max(self):
        """
        It gets the maximum coordinates of the spatial region enclosed by this ``Node``.
        
        :return: the maximum coordinates of the spatial region enclosed by this ``Node``.
        :rtype: numpy.float([float][float][float])
        """
        return self.__max

    def get_data(self):
        """
        It gets the ``list`` of ``Data`` objects inside this ``Node``.
        
        :return: the ``Data`` objects inside this 
        :rtype: 
        """
        return self.__data

    def get_num_points(self):
        """
        It gets the number of points inside this ``Node``
        
        :return: number of points inside this ``Node``.
        :rtype: float
        """
        return len(self.__data)

    def has_children(self):
        """
        It checks if this ``Node`` has children nodes.
        
        :return: ``True`` if it has children nodes created, ``False`` otherwise. 
        :rtype: boolean
        """
        return len(self.__children) != 0

    def has_data(self):
        """
        It checks if this ``Node`` has any data stored.
        
        :return: ``True`` if it has data stored, ``False`` otherwise.
        :rtype: boolean
        """
        return len(self.__data) != 0

    def get_level(self):
        """
        It gets the depth level where this ``Node`` has been created in the internal ``Octree`` structure.
        
        :return: the depth level of this ``Node``.
        :rtype: int
        """
        return self.__level

    def see_node(self):
        """
        It generates geometric structures for the ``open3D`` library for rendering the octree. It adds the geometric structures to the attrubute ``visual_nodes`` of the ``Octree`` data structure.
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

        if (self.__level >= 0):
        #if (self.has_data()):
            self.__octree.visual_nodes.append(line_set)
        for c in self.__children:
            c.see_node()

class Octree:
    """
    It manages the internal and root structure of an octree data structure.
    """
    
    def __init__(self, level_limit, data_limit, data, min_v, max_v):
        """
        """
        self.__min = min_v
        self.__max = max_v
        self.__level_limit = level_limit
        self.__data_limit = data_limit
        self.__num_data = len(data)
        self.__root = Node(0, None, self.__min, self.__max, self)
        progress = 0
        total = len(data)
        for p in data:
            self.__root.add(p)
            # progress...
            if progress % 100000 == 0:
                print(progress / total * 100, "%")
            progress += 1
        self.visual_nodes = list()

    def search(self, p):
        """
        It calls the ``search`` method from the root ``Node`` of this ``Octree``.

        .. seealso::
            The method :py:meth:`thomsonpy.data_management.octree.octree.Node.search` using the recursivity to search for the ``p`` passed as parameter.
        
        :return: a ``Data`` object representing the result of the searching process. If the search is successful, it returns the ``Data`` object related with the 3D point ``p`` passed as parameter, otherwise it returns a ``None`` object.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        return self.__root.search(p)
    
    def search_nearest(self, p):
        """
        It calls the ``search_nearest`` method from the root ``Node`` of this ``Octree``.

        .. seealso::
            The method :py:meth:`thomsonpy.data_management.octree.octree.Node.search_nearest` using the recursivity to search for the ``p`` passed as parameter.
        
        :return: a ``Data`` object representing the result of the nearest searching process. If the search is successful, it returns the nearest ``Data`` object to the 3D point ``p`` passed as parameter, otherwise it returns a ``None`` object.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Data`
        """
        return self.__root.search_nearest(p)

    def get_level_limit(self):
        """
        It gets the level limit used as a restriction for the construction process. If this level is reached the ``Node`` objects at this level do not create more children nodes anymore.
        
        :return: the level limit.
        :rtype: int
        """
        return self.__level_limit

    def get_data_limit(self):
        """
        It gets the data limit used as a restriction for the construction process. If this quantity of ``Data`` objects in a particular ``Node`` is reached and a new ``Data`` object is inserted in that ``Node``, then that ``Node`` creates eight children if ``level_limit`` has not been reached yet.
        
        :return: the data_limit.
        :rtype: int
        """
        return self.__data_limit

    def get_num_data(self):
        """
        It gets the number of ``Data`` objects stored in this ``Octree``.
        
        :return: 
        :rtype: int
        """
        return self.__num_data

    def get_root(self):
        """
        It gets the root of this ``Octree`` from where the rest of the nodes are processed recursively.
        
        :return: the root ``Node`` of this ``Octree``.
        :rtype: :py:class:`thomsonpy.data_management.octree.octree.Node`
        """
        return self.__root

    def get_visual_octree(self):
        """
        It gets a ``list`` with all the geometric objects representing the internal structure of this ``Octree``.
        
        :return: a list with geometric objects readable by the ``Open3D`` library.
        :rtype: ``list``
        """
        self.__root.see_node()
        return self.visual_nodes

    @staticmethod
    def save(octree, filename):
        """
        It saves the octree structure in a binary file.
        """
        f = open(filename, 'wb')
        pickle.dump(octree, f)
        f.close()

    @staticmethod
    def load(filename):
        """
        It loads an octree structure from a binary file.
        """
        f = open(filename, 'rb')
        obj = pickle.load(f)
        f.close()
        return obj