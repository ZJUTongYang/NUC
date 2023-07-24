#!/usr/bin/python3
import nuc
import numpy as np
import trimesh

stl_data = trimesh.load_mesh('remeshed_saddle.stl')
mesh_tri = stl_data.faces
mesh_ver = stl_data.vertices

mesh_tri_row = np.reshape(mesh_tri, (mesh_tri.shape[0]*mesh_tri.shape[1], -1))
mesh_ver_row = np.reshape(mesh_ver, (mesh_ver.shape[0]*mesh_ver.shape[1], -1))

result = nuc.add(mesh_tri_row, mesh_ver_row)

result = np.reshape(result, (int(len(result)/3), 3))

## We output the data (for the visualisation in MATLAB)
mesh_tri_f = open("mesh_tri.txt", "w+")
mesh_ver_f = open("mesh_ver.txt", "w+")
path_f = open("result_path.txt", "w+")

np.savetxt(mesh_tri_f, mesh_tri)
np.savetxt(mesh_ver_f, mesh_ver)
np.savetxt(path_f, result)

mesh_tri_f.close()
mesh_ver_f.close()
path_f.close()

