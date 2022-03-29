from thomsonpy.data_management.octree.octree import *
import thomsonpy.config.paths as paths
import thomsonpy.data_management.formatter as formatter
import time

my_octree = Octree.load(paths.OCTREE_OBJECTS_PATH + "octree_1.oct")
my_data = formatter.load(paths.OCTREE_DATA_PATH + "octree_data_1.obj")
count = 0
print("Start searching")

"""
for p in my_data:
    init_time = time.perf_counter()
    my_octree.search_nearest(p.get_coordinates())
    fin_time = time.perf_counter()
    print(str(fin_time - init_time), "seconds")
"""

"""
for p in my_data:
    init_time = time.perf_counter()
    nearest = None
    for q in my_data:
        min_distance = sys.float_info.max
        distance = np.linalg.norm(p.get_coordinates() - q.get_coordinates())
        if distance < min_distance:
            nearest = q
            distance = min_distance   
    fin_time = time.perf_counter()
    print(str(fin_time - init_time), "seconds") 
    
"""
listaFOUND = list()
listaNOTFOUND = list()

count = 0
found = 0
not_found = 0
for p in my_data:
    if my_octree.search(p.get_coordinates()) is None:
        not_found += 1
        listaNOTFOUND.append(p.get_coordinates())
    else:
        found += 1
        listaFOUND.append(p.get_coordinates())
    if count % 100000 == 0:
        print(count / len(my_data) * 100)
        print("FOUND =", found)
        print("NOT FOUND =", not_found)
    count += 1
       
sphere = o3d.geometry.TriangleMesh.create_sphere(radius=1.0)
sphere.paint_uniform_color([0.8, 0.5, 0.0])

pcdFOUND = o3d.geometry.PointCloud()
pcdFOUND.points = o3d.utility.Vector3dVector(listaFOUND)
pcdFOUND.paint_uniform_color([1, 0, 0])

pcdNOTFOUND = o3d.geometry.PointCloud()
pcdNOTFOUND.points = o3d.utility.Vector3dVector(listaNOTFOUND)
pcdNOTFOUND.paint_uniform_color([0, 1, 0])

viewer = o3d.visualization.Visualizer()
viewer.create_window()
viewer.add_geometry(pcdFOUND)
viewer.add_geometry(pcdNOTFOUND)
viewer.add_geometry(sphere)
opt = viewer.get_render_option()
opt.show_coordinate_frame = True
opt.background_color = np.asarray([0.0, 0.2, 0.2])
viewer.run()
viewer.destroy_window()
