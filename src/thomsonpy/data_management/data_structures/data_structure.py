from abc import ABC, abstractmethod


class DataStructure(ABC):

    @abstractmethod
    def search(self, coordinates):
        raise NotImplementedError

    @abstractmethod
    def search_nearest(self, coordinates):
        raise NotImplementedError

    @abstractmethod
    def get_representation(self):
        raise NotImplementedError
