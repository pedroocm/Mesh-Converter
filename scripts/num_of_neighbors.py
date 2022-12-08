import pymesh
import sys


THRESHOLD = 50
PATH = sys.argv[1]

mesh = pymesh.load_mesh(PATH)
max_num_of_neighbors = 0

mesh.enable_connectivity()
for v_id, coords in enumerate(mesh.vertices):
    num_neighbors = len(mesh.get_vertex_adjacent_vertices(v_id))
    if(num_neighbors > THRESHOLD):
        print(f"{v_id}: {num_neighbors}")
    if(num_neighbors > max_num_of_neighbors):
        max_num_of_neighbors = num_neighbors

print(f"Max number of neighbors: {max_num_of_neighbors}")
        
