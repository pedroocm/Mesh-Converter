#!/usr/bin/env python

import pymesh
import numpy as np
import math
import sys
import os
from lib.transformData import map_indexes

DEBUG = False

def generatePocosInput(filename, args):
    nx, ny, nz = args

    with open(filename, "w") as file:
        fileString = f'''{nx} {ny} {nz}
1 1 1:{nz}
{nx} 1 1:{nz}
1 {ny} 1:{nz}
{nx} {ny} 1:{nz}'''

        file.write(fileString)

def read_input_file(file_name):
    HEADER = True
    indexes = []
    dimensions = []
    zRange_array = []
    with open(file_name) as file:
        for line in file:
            if(HEADER):
                dimensions = line.split()
                HEADER = False
            else:
                (i, j, zRange_string) = line.split()
                (zBeg, zEnd) = zRange_string.split(":")
                indexes.append((i, j))
                zRange_array.append((zBeg, zEnd))

    return (indexes, dimensions, zRange_array)

def lua_table_to_python_list(string):
    # String format -> var = {1, 2, 3, 4, 5, 6, 7, ...}
    numbers_list = string.split('=')[-1].strip()
    py_list = numbers_list[1:-1].split(',')
    return py_list

def read_input_lua_file(file_name):
    header = []
    LXW = []
    LYW = []
    LZWF = []
    LZWL = []
    
    with open(file_name, 'r') as file:
        for linha in file:
            if linha.lstrip().startswith('NX'):
                header.append(linha.strip()[:-1])
                header.append(file.readline().strip()[:-1])
                header.append(file.readline().strip()[:-1])
            elif(linha.lstrip().startswith('LXW')):
                LXW = lua_table_to_python_list(linha.strip()[:-1])
                LYW = lua_table_to_python_list(file.readline().strip()[:-1])
                LZWF = lua_table_to_python_list(file.readline().strip()[:-1])
                LZWL = lua_table_to_python_list(file.readline().strip()[:-1])
    
    dimensions = [c.split('=')[-1].strip() for c in header]
    indexes = [t for t in zip(LXW, LYW)]
    zRange_array = [t for t in zip(LZWF, LZWL)]

    return (indexes, dimensions, zRange_array)

def write_output(filename, well_vertices_index, indexes, zRange, nwbcm):
    with open(filename, 'w') as file:
        file.write("NWBCM\n")
        file.write(f"\t{nwbcm}\n")
        for index, well in enumerate(well_vertices_index):
            file.write(f"({indexes[index][0]}, {indexes[index][1]}, {zRange[index][0]})/")
            file.write(f"({indexes[index][0]}, {indexes[index][1]}, {zRange[index][1]}): ")
            for vertex in well:
                file.write(str(vertex))
                file.write(", ")
            file.write('\n')


def generatePocosOutput(filename, original_mesh_name, converted_mesh_name, inactiveCells = []):
    
    if(filename.endswith('lua')):
        indexes, dimensions, zRange_array = read_input_lua_file(filename)
    else:
        indexes, dimensions, zRange_array = read_input_file(filename)
    
    print(f"Dimensions: {dimensions}\n")
    print(f"Indexes: {indexes}\n")
    print(f"Z Range: {zRange_array}\n")

    original_mesh = pymesh.load_mesh(original_mesh_name)
    converted_mesh = pymesh.load_mesh(converted_mesh_name)
    locked_vertices = {}


    (nx, ny, nz) = (int(dimensions[0]), int(dimensions[1]), int(dimensions[2]) )
    well_indexes = []
    well_enumerate = 0
    for index, zRange in zip(indexes, zRange_array):
        well_enumerate += 1
        print(f'Well {well_enumerate}')
        i, j = (int(index[0]) - 1, int(index[1]) - 1)
        zBeg, zEnd = (int(zRange[0]) - 1, int(zRange[1]) - 1)
        
        column_coord = j * nz + i * nz * ny
        first_depth = -1
        first_element_coord, second_element_coord = -1, -1
        for k in range(zBeg, zEnd + 1):
            if(inactiveCells[column_coord + k]):
                first_element_coord = column_coord + k
                first_depth = k
                break
        
        if(first_element_coord == -1):
            print("Column of inactive cells")
            continue

        for k in range(zEnd, first_depth-1, -1):
            if(inactiveCells[column_coord + k]):
                second_element_coord = column_coord + k
                break
        
        if(DEBUG):
            print('First Element:', first_element_coord, 'Mapped:', map_indexes(inactiveCells, first_element_coord))
            print('Second Element:', second_element_coord, 'Mapped:', map_indexes(inactiveCells, second_element_coord))
        
        first_selectedElement = original_mesh.elements[map_indexes(inactiveCells, first_element_coord)]
        second_selectedElement = original_mesh.elements[map_indexes(inactiveCells, second_element_coord)]

        # ray
        vector_endpoints = (original_mesh.vertices[ second_selectedElement[0] ], original_mesh.vertices[ first_selectedElement[3] ])
        ray_origin = vector_endpoints[1]
        ray_vector_direction =  vector_endpoints[0] - ray_origin
        
        if(DEBUG):
            print("Ray Origin:", vector_endpoints[1])
            print("Ray End:", vector_endpoints[0])
            print("Ray Direction:", ray_vector_direction    )
        
        vertices_index = []
        for (index, vertex) in enumerate(converted_mesh.vertices):
            t_z = (vertex[2] - ray_origin[2]) / ray_vector_direction[2]

            if(0 <= t_z <= 1):
            
                estimated_x = t_z * ray_vector_direction[0] + ray_origin[0]
                estimated_y = t_z * ray_vector_direction[1] + ray_origin[1]

                if( np.isclose(estimated_x, vertex[0]) and np.isclose(estimated_y, vertex[1]) ): 
                    vertices_index.append(index + 1)
                    locked_vertices[index] = True

        vertices_index.sort(key=lambda i: converted_mesh.vertices[i - 1][2])
        well_indexes.append(vertices_index)

    nwbcm = 0
    for well in well_indexes:
        if(len(well) > nwbcm):
            nwbcm = len(well)
    
    path, _ = os.path.split(filename)
    output_path = os.path.join(path, "pocos_output.txt")
    
    write_output(output_path, well_indexes, indexes, zRange_array, nwbcm)
    return locked_vertices
         

if __name__ == "__main__":
    
    if len(sys.argv) != 4:
        print("Número incorreto de argumentos")
        sys.exit(1)
    
    generatePocosOutput(sys.argv[1], sys.argv[2], sys.argv[3])
