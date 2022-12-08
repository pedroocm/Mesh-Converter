import pymesh
import sys
import numpy as np

if(len(sys.argv) != 2):
    exit(1)

path = sys.argv[1]

mesh = pymesh.load_mesh(path)

debug_file_path = f"{path[:-4]}_collapsed_edges.txt"
file = open(debug_file_path, "w")
collapsed_voxels = []

for cell_id, voxel_array in enumerate(mesh.voxels):
    collapsed_edges = {}
    has_collapsed_edge = False
    voxel = list(voxel_array)

    for vertex_id in voxel:
        if(vertex_id not in collapsed_edges and voxel.count(vertex_id) > 1):
            has_collapsed_edge = True
            collapsed_edges[vertex_id] = [i for i in range(len(voxel)) if voxel[i] == vertex_id]
    
    if(has_collapsed_edge):
        file.write(f"Cell Id: {cell_id}\n")
        file.write(f"Vertices Id: {voxel}\n")
        file.write(f"Collapsed Vertices:\n")
        for v_id, v_array in collapsed_edges.items():
            file.write(f"{v_id}: {[i + 1 for i in v_array]}\n")
        file.write("\n")

print(f"Created file: {debug_file_path}")
file.close() 