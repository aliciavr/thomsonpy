"""
.. module:: visualizer
        :platform: Unix
        :synopsis: tools for manage data format, usage and storing 
.. moduleauthor:: 
"""
import os
import open3d as o3d
import numpy as np
import thomsonpy.data_management.octree.octree as octr
import thomsonpy.data_management.formatter as formatter
import thomsonpy.config.paths as paths
import thomsonpy.constants.units as units

def vis_points_and_ne(i):
    points_filenames = ['points_1.obj', 'points_2.obj', 'points_3.obj', 'points_4.obj']
    pclouds = list()
    ne_filenames = ['ne_1.obj', 'ne_2.obj', 'ne_3.obj', 'ne_4.obj']
    ne_clouds = list()

    pcloud = formatter.load(paths.OCTREE_DATA_PATH + points_filenames[i])
    pclouds.append(pcloud)
    ne_cloud = formatter.load(paths.OCTREE_DATA_PATH + ne_filenames[i])
    ne_clouds.append(ne_cloud)
    sphere = o3d.geometry.TriangleMesh.create_sphere(radius=1.0)
    sphere.paint_uniform_color([0.8, 0.5, 0.0])

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pclouds[0])
    colors = list()
    for ne_val in ne_clouds[0]:
        color = ne_val
        colors.append(np.array([color, color, color]))
    colors = np.array(colors)
    pcd.colors = o3d.utility.Vector3dVector(colors)

    viewer = o3d.visualization.Visualizer()
    viewer.create_window()
    viewer.add_geometry(pcd)
    viewer.add_geometry(sphere)
    opt = viewer.get_render_option()
    opt.show_coordinate_frame = True
    opt.background_color = np.asarray([0.0, 0.2, 0.2])
    viewer.run()
    viewer.destroy_window()
    del pcloud
    del ne_cloud
    del pcd
    
def vis_octree(i):

    octrees_filenames = ["octree_1.oct", "octree_2.oct", "octree_3.oct", "octree_4.oct"]
    octree = Octree.load(paths.OCTREE_OBJECTS_PATH + octrees_filenames[i])
    octree_geom = octree.get_visual_octree()
    sphere = o3d.geometry.TriangleMesh.create_sphere(radius=1.0)
    sphere.paint_uniform_color([0.8, 0.5, 0.0])
    viewer = o3d.visualization.Visualizer()
    viewer.create_window()
    for i in octree_geom:
            viewer.add_geometry(i)
    viewer.add_geometry(sphere)
    opt = viewer.get_render_option()
    opt.show_coordinate_frame = True
    opt.background_color = np.asarray([0.0, 0.2, 0.2])
    viewer.run()
    viewer.destroy_window()

def vis_octree_data():
    my_data = formatter.load(f"{paths.OCTREE_DATA_PATH}{ os.path.splitext(paths.PREDSCI_FILENAME)[0]}.data")

    lista = list()
    for p in my_data:
        lista.append(p.get_coordinates() * units.METERS_TO_RSOL)
        
    sphere = o3d.geometry.TriangleMesh.create_sphere(radius=1.0)
    sphere.paint_uniform_color([0.8, 0.5, 0.0])

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(lista)
    pcd.paint_uniform_color([1, 0, 0])

    viewer = o3d.visualization.Visualizer()
    viewer.create_window()
    viewer.add_geometry(pcd)

    viewer.add_geometry(sphere)
    opt = viewer.get_render_option()
    opt.show_coordinate_frame = True
    opt.background_color = np.asarray([0.0, 0.2, 0.2])
    viewer.run()
    viewer.destroy_window()

def vis_found_and_not_found(i):
    my_octree = Octree.load(paths.OCTREE_OBJECTS_PATH + "octree_" + str(i) + ".oct")
    my_data = formatter.load(paths.OCTREE_DATA_PATH + "octree_data_" + str(i) + ".obj")

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
vis_octree_data()