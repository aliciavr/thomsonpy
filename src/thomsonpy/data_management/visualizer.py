import open3d as o3d
import numpy as np
import pickle

filenames = ['points_1.obj', 'points_2.obj', 'points_3.obj', 'points_4.obj']

f = open(filenames[3], 'rb')
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(pickle.load(f))
f.close()
viewer = o3d.visualization.Visualizer()
viewer.create_window()
viewer.add_geometry(pcd)
opt = viewer.get_render_option()
opt.show_coordinate_frame = True
opt.background_color = np.asarray([0.5, 0.5, 0.5])
viewer.run()
viewer.destroy_window()
del(pcd)