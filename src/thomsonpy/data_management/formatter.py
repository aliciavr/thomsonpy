from thomsonpy.data_management.octree.octree import Data

def from_numpy_to_octree_data(xyz, ne):
    progress = 0
    data = list()
    for i, ne_val in enumerate(ne):
        d = Data(xyz[i], ne_val)
        data.append(d)
        
        # progress...
        if progress % 2000000 == 0:
            print(progress)
        
        progress += 1
    return data
