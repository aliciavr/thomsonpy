import numpy as np
import sys
import pickle
import thomsonpy.data_management.formatter as fmt

class SphericalMesh:
    
    def __init__(self, nd_phi, nd_theta, nd_radial, radial_max):
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
        sector = self.get_sector(fmt.cartesian_to_spherical(data.get_coordinates()))
        sector.add_data(data)
        
    def get_sector(self, coords):
        i = int((coords[0]//self.__length_phi)) % self.__nd_phi
        j = int((coords[1]//self.__length_theta)) % self.__nd_theta
        k = int((coords[2]//self.__length_radial)) % self.__nd_radial
        sector = self.__mesh[i][j][k] 
        return sector
    
    def get_sector_and_neighbors(self, coords):
        i = int((coords[0]//self.__length_phi)) % self.__nd_phi
        j = int((coords[1]//self.__length_theta)) % self.__nd_theta
        k = int((coords[2]//self.__length_radial)) % self.__nd_radial
        sectors = []
        num_data = 0
        r = [-1,0,1]
        for ii in r:
            for jj in r:
                for kk in r:
                    sector = self.__mesh[(i + ii) % self.__nd_phi][(j + jj) % self.__nd_theta][(k + kk) % self.__nd_radial]
                    sectors.append(sector)
                    num_data += sector.get_num_data()
        
        print(num_data)
        return sectors
    
    def search_nearest(self, coords):
        sectors = self.get_sector_and_neighbors(fmt.cartesian_to_spherical(coords))
        nearest = Sector3D.search_nearest(coords, sectors)
        return nearest
    
    @staticmethod
    def save(sphmesh, filename):
        """
        It saves the spherical mesh structure in a binary file.
        """
        f = open(filename, 'wb')
        pickle.dump(sphmesh, f)
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
        
    @staticmethod
    def search_nearest(coords, sectors):
        nearest = Data(None, 0)
        min_distance = sys.float_info.max
        for sector in sectors:      
            for data in sector.get_all_data():
                distance = np.linalg.norm(coords - data.get_coordinates())
                if distance < min_distance:
                    nearest = data
                    distance = min_distance
        return nearest
    
    def add_data(self, data):
        self.__data.append(data)
        
    def get_all_data(self):
        return self.__data
        
    def get_num_data(self):
        return len(self.__data)
        
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