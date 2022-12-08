# coding: utf-8
import pymesh
import random
mesh = pymesh.load_mesh('uniao_tetrav.msh')
mesh.enable_connectivity()

counter = 0
elementsWithTwoNeigh = []

for i in range(mesh.num_elements):
    if len(mesh.get_voxel_adjacent_voxels(i)) == 2:
    	elementsWithTwoNeigh.append(i + 1)


print(f"Number of elements: {len(elementsWithTwoNeigh)}")

chosenElement = elementsWithTwoNeigh[ random.randint(0, len(elementsWithTwoNeigh) - 1) ]
print(f"Check element #{chosenElement} on GMSH. Its neighbours are: {[elem + 1 for elem in mesh.get_voxel_adjacent_voxels(chosenElement - 1)]}")
        
