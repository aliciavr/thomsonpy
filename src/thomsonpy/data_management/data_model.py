class DataModel:

    def __init__(self, coordinates, physical_data):
        self.__model = dict()
        self.__model_coordinates = dict()
        num_points = len(coordinates)
        for i, c in zip(num_points, coordinates):
            self.__model[i] = Data(coordinates, physical_data)
            self.__model_coordinates[c] = i

    def __init__(self, data):
        self.__model = dict()
        self.__model_coordinates = dict()
        num_points = len(data)
        for i, d in zip(num_points, data):
            self.__model[i] = d
            self.__model_coordinates[d.get_coordinates()] = i

    def get_all_indices(self):
        return self.__model.keys()

    def get_all_data(self):
        return self.__model.values()

    def get_all_coordinates(self):
        return self.__model_coordinates.keys()

    def get_data_from_index(self, index):
        return self.__model.get(index)

    def get_data_from_coordinates(self, coordinates):
        """
        It gets the index that identifies the coordinates of a ``Data`` element in the physical model.

        :param coordinates: the Cartesian coordinates which identify the ``Data`` elements.
        :type coordinates: numpy.ndarray(float, float, float)

        :return:
        :rtype:
        """
        index = self.get_data_index(coordinates)
        return self.__model.get(index)

    def get_index(self, coordinates):
        """
        It gets the index that identifies the coordinates of a ``Data`` element in the physical model.

        :param coordinates:
        :type coordinates:

        :return:
        :rtype:
        """
        return self.__model_coordinates[coordinates]

    def get_coordinates(self, index):
        data = self.get_data_from_index(index)
        return data.get_coordinates()

    def get_x(self, index):
        data = self.get_data_from_index(index)
        return data.get_x()

    def get_y(self, index):
        data = self.get_data_from_index(index)
        return data.get_y()

    def get_z(self, index):
        data = self.get_data_from_index(index)
        return data.get_z()

    def get_ne(self, index):
        data = self.get_data_from_index(index)
        return data.get_ne()

    def get_temperature(self, index):
        data = self.get_data_from_index(index)
        return data.get_temperature()

    def get_physical_data(self, index):
        data = self.get_data_from_index(index)
        return data.get_physical_data()


class Data:
    __NE = 0
    __T = 1

    def __init__(self, coordinates, physical_data):
        self.__coordinates = coordinates
        self.__physical_data = physical_data

    def get_x(self):
        return self.__coordinates[0]

    def get_y(self):
        return self.__coordinates[1]

    def get_z(self):
        return self.__coordinates[2]

    def get_coordinates(self):
        return self.__coordinates

    def get_ne(self):
        return self.__physical_data[self.__NE]

    def get_temperature(self):
        return self.__physical_data[self.__T]

    def get_physical_data(self):
        return self.__physical_data
