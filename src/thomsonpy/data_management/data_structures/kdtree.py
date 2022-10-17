from thomsonpy.data_management.data_structures.data_structure import DataStructure
import thomsonpy.data_management.utils as utils
from sklearn.neighbors import KDTree


class KDTree(DataStructure):

    def __init__(self, coordinates, return_distance=True,
                 dual_tree=False,
                 breadth_first=False,
                 sort_results=True,
                 metric='minkowski',
                 leaf_size=2,
                 filename='kdtree',
                 path='.obj'):
        self.__return_distance = return_distance
        self.__dual_tree = dual_tree
        self.__breadth_first = breadth_first
        self.__sort_results = sort_results
        self.__metric = metric
        self.__leaf_size = leaf_size
        self.__kdtree = KDTree(X=coordinates, leaf_size=self.__leaf_size, metric=self.__metric)
        utils.store_object(self, filename, path)

    def search_nearest_data_index(self, coordinates):
        """
        It searches for the index with the nearest coordinates associated to the ones passed as a parameter.
        It assumes receiving the coordinates of a single point.

        :param coordinates: coordinates of a single point in 3D space.
        :type coordinates: numpy.ndarray([float, float, float]).

        :return: the index which identifies the nearest ``Data`` object to the coordinates passed as a parameter.
        :rtype: int
        """
        if self.__return_distance:
            distance, index = self.__kdtree.query(X=coordinates, k=1,
                                                  return_distance=self.__return_distance,
                                                  dualtree=self.__dual_tree,
                                                  breadth_first=self.__breadth_first)
            return distance, index
        else:
            index = self.__kdtree.query(X=coordinates, k=1)
            return index

    def search_nearest(self, coordinates, data_model):
        """
        It searches for the data object associated with the index retrieved by
        :python:func:`thomsonpy.data_management.data_structures.kdtree.search_nearest_index`.
        It assumes receiving the coordinates of a single point.

        :param coordinates: coordinates of a single point in 3D space.
        :type coordinates: numpy.ndarray([float, float, float]).

        :return: the nearest ``Data`` object to the coordinates passed as a parameter.
        :rtype: ``thomsonpy.data_management.data_model.Data``

        """
        if self.__return_distance:
            distance, nearest_index = self.search_nearest_data_index(coordinates)
        else:
            nearest_index = self.search_nearest_data_index(coordinates)

        data = data_model.get_data_from_index(nearest_index)
        return data

    def search_k_nearest_indices(self, coordinates, k):
        """
        It searches for the :math:`k` ``Data`` objects.
        """
        if self.__return_distance:

            distances, nearest_indices = self.__kdtree.query(X=coordinates, k=k,
                                                             return_distance=self.__return_distance,
                                                             dualtree=self.__dual_tree,
                                                             breadth_first=self.__breadth_first)

        else:
            nearest_indices = self.__kdtree.query(X=coordinates, k=k,
                                                  return_distance=self.__return_distance,
                                                  dualtree=self.__dual_tree,
                                                  breadth_first=self.__breadth_first)

        return nearest_indices

    def search(self, coordinates):
        data = self.search_nearest(coordinates)
        if data.get_x() != coordinates[0] | data.get_y() != coordinates[1] | data.get_z() != coordinates[2]:
            data = None
        return data

    def get_representation(self):
        return self.__kdtree.get_tree_stats()
