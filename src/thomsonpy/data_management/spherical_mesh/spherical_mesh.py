import numpy as np
import sys

class SphericalMesh:
    
    def __init__(self, nd_theta, nd_phi, nd_radial, radial_max):
        self.__nd_phi = nd_phi
        self.__nd_theta = nd_theta
        self.__nd_radial = nd_radial
        self.__radial_max = radial_max
        
        self.__length_phi = np.pi / nd_phi
        self.__length_theta = 2 * np.pi / nd_theta
        self.__length_radial = radial_max / nd_radial

        self.__mesh =  [[[None]*nd_radial]*nd_theta]*nd_phi
        for i in range(nd_phi):
            for j in range (nd_theta):
                for k in range (nd_radial):
                    self.__mesh[i][j][k] = Sector3D()

    def add_data(self, data):
        sector = self.get_sector(data.get_phi(), data.get_theta(), data.get_radial())
        sector.__add_data(data)
        
    def get_sector(self, phi, theta, radial):
        i = int((phi//self.__length_phi)) % self.__nd_phi
        j = int((theta//self.__length_theta)) % self.__nd_theta
        k = int((radial//self.__length_radial)) % self.__nd_radial
        print(i, j, k)
        sector = self.__mesh[i][j][k] 
        return sector
    
    def search_nearest(self, coords):
        
        phi = 
        theta = 
        radial = 
        sector = get_sector(phi, theta, radial)
        nearest = sector.__search_nearest()
        return nearest
    
    @staticmethod
    def save(sphmesh, filename):
        """
        It saves the spherical mesh structure in a binary file.
        """
        f = open(filename, 'wb')
        pickle.dump(octree, f)
        f.close()

    @staticmethod
    def load(filename):
        """
        It loads an spherical mesh structure from a binary file.
        """
        f = open(filename, 'rb')
        obj = pickle.load(f)
        f.close()
        return obj
    
class Sector3D:
    def __init__(self):
        self.__data = []
        
    def __search_nearest(self):
        nearest = Data(None, 0)
        min_distance = sys.float_info.max
        for data in self.__data:
            distance = np.linalg.norm(p - data.get_coordinates())
            if distance < min_distance:
                nearest = data
                distance = min_distance
        return nearest
    
    def __add_data(self, data):
        self.__data.append(data)
        
class Data:

    def __init__(self, coordinates, ne):

        self.__coordinates = coordinates
        self.__ne = ne

    def __str__(self):
        return 'Data = {} at ({}, {}, {})'.format(self.__ne, self.__coordinates[0], self.__coordinates[1], self.__coordinates[2])

    def get_coordinates(self):

        return self.__coordinates

    def get_phi(self):
        return self.__coordinates[0]

    def get_theta(self):
        return self.__coordinates[1]

    def get_radial(self):
        return self.__coordinates[2]

    def get_ne(self):
        """
        It gets the electron density :math:`\\rho` related with the point represented by this ``Data`` object.
        
        :return: the electron density :math:`\\rho` of this ``Data`` object.
        :rtype: float
        """
        return self.__ne
        
"""
nd_phi = 4
nd_theta = 3
nd_radial = 2
radial_max = 10
mesh = SphericalMesh(nd_phi, nd_theta, nd_radial, radial_max)


phi = np.pi
theta = np.pi / 2
radial = 10
print(mesh.add_data())
"""