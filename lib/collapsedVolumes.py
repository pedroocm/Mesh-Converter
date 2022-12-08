import pymesh
import numpy as np
import sys
import os

threshold = 10**-1
DEBUG = False
#DEBUG = True

face_indices = {
    'TOP': [0, 1, 4, 5],
    'BOTTOM': [2, 3, 6, 7],
}

def avgVolumeCritiria(idV, hexa_volumes, avg_volume):
    return np.divide(hexa_volumes[idV],avg_volume) <= threshold

def thicknessCritiria(idV, voxels, vertices):
    voxel = voxels[idV]

    avg_z_top = np.mean([vertices[voxel[i]][2] for i in face_indices["TOP"]])
    avg_z_bottom = np.mean([vertices[voxel[i]][2] for i in face_indices["BOTTOM"]])

    val = abs(avg_z_bottom - avg_z_top) <= threshold
    if(val):
        print(f"Bloco: {idV}   Value: {val}")
    return val

def getCollapsedVolumes(hexa_volumes, num_hexa, voxels, vertices, pinchArray, volume_critiria_is_enabled):
    avg_volume = np.mean(hexa_volumes)
    print(f"Average volume: {avg_volume}")
    print(f"Number of hexa: {num_hexa}")
    collapsedVolumes = []

    for idV in range(num_hexa):
        #if(not pinchArray[idV]):
        #    collapsedVolumes.append(idV)
        #elif
        if(thicknessCritiria(idV, voxels, vertices)):
            collapsedVolumes.append(idV)
        elif(volume_critiria_is_enabled and avgVolumeCritiria(idV, hexa_volumes, avg_volume)):
            collapsedVolumes.append(idV)
    return collapsedVolumes

def getDirection(id_voxel, id_neighbor, voxels):
    
    shared_face_indices = []
    voxel = voxels[id_voxel]
    neighbor = voxels[id_neighbor]

    for i in range(8):
        if(voxel[i] in neighbor):
            shared_face_indices.append(i)
    
    for direction, index_list in face_indices.items():
        if(index_list == shared_face_indices):
            return direction
    
    return None

def collapseVolume(id_collapsed, id_neighbor, voxels, direction, collapsed_blocks, mesh):
    
    collapsed = False

    for block in collapsed_blocks:

        if(block[0] == id_collapsed and direction == "BOTTOM"):
            collapsed = True
            block[0], block[1] = block[1], block[0]
            neighbors = mesh.get_voxel_adjacent_voxels(block[0])
            for neighbor in neighbors:
                if(getDirection(block[0], neighbor, voxels) == "BOTTOM"):
                    block[0] = neighbor

        if(block[0] == id_collapsed):
            collapsed = True
            block[0] = id_neighbor
        

    if(not collapsed):
        block = [id_neighbor, id_collapsed]
        collapsed_blocks.append(block)


def expand_active_block(active_cell_id, collapsed_cell_id, voxels, vertices, threshold=0.05):
    collapsed_cell = voxels[collapsed_cell_id]
    active_cell = voxels[active_cell_id]

    expand_direction = ''

    if(vertices[active_cell[0]][2] <= vertices[collapsed_cell[0]][2]):
        expand_direction = 'TOP'
    else:
        expand_direction = 'BOTTOM'

    for i in face_indices[expand_direction]:
        distance = abs(vertices[active_cell[i]][2] - vertices[collapsed_cell[i]][2])
        if(distance > threshold):
            return

    for i in face_indices[expand_direction]:
        active_cell_vertex_id = active_cell[i]
        collapsed_cell_vertex_id = collapsed_cell[i]
        
        vertices[active_cell_vertex_id] = vertices[collapsed_cell_vertex_id]

def removeCollapsedVolumes(collapsedVolumes, voxels):
    return np.delete(voxels, collapsedVolumes, 0)

def writeMSH(voxels, points, numHexa, nome):
    with open(nome + ".msh", "w") as newFile:
        newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        newFile.write("$Nodes\n" + str(len(points)) + "\n")
        for (index, point) in enumerate(points):
            newFile.write(str(
                index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
        newFile.write("$EndNodes\n")
        newFile.write("$Elements\n" + str(numHexa) + "\n")

        tagNumber = 0

        for (idV, voxel) in enumerate(voxels):
            newFile.write(
                str(idV + 1) + f" 5 2 0 {tagNumber} " +
                str(voxel[0] + 1) + " " + str(voxel[1] + 1) + " "
                + str(voxel[2] + 1) + " " + str(voxel[3] + 1) + " "
                + str(voxel[4] + 1) + " " + str(voxel[5] + 1) + " "
                + str(voxel[6] + 1) + " " + str(voxel[7] + 1) + "\n"
            )
        newFile.write("$EndElements")

def highlightCollpsedVolumes(collapsedVolumes, voxels, points, numHexa, nome):

    with open(nome + ".msh", "w") as newFile:
        newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        newFile.write("$Nodes\n" + str(len(points)) + "\n")
        for (index, point) in enumerate(points):
            newFile.write(str(
                index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
        newFile.write("$EndNodes\n")
        newFile.write("$Elements\n" + str(numHexa) + "\n")

        tagNumber = 0

        for (idV, voxel) in enumerate(voxels):
            if(idV in collapsedVolumes):
                tagNumber = 1
            else:
                tagNumber = 0
            newFile.write(
                str(idV + 1) + f" 5 2 0 {tagNumber} " +
                str(voxel[0] + 1) + " " + str(voxel[1] + 1) + " "
                + str(voxel[2] + 1) + " " + str(voxel[3] + 1) + " "
                + str(voxel[4] + 1) + " " + str(voxel[5] + 1) + " "
                + str(voxel[6] + 1) + " " + str(voxel[7] + 1) + "\n"
            )
        newFile.write("$EndElements")

def main(path_file, inactiveCells = [], pinchArray = [], volume_critiria_is_enabled=False):
    mesh = pymesh.load_mesh(path_file)
    mesh.add_attribute("voxel_volume")
    
    voxels = mesh.voxels.copy()
    vertices = mesh.vertices.copy()
    hexa_volumes = mesh.get_attribute("voxel_volume")
    
    collapsed_volumes = getCollapsedVolumes(
        hexa_volumes, 
        mesh.num_voxels, 
        voxels, 
        vertices,
        pinchArray,
        volume_critiria_is_enabled
    )
    collapsed_blocks = []
    print(f"collapsed_volumes: {collapsed_volumes}")

    mesh.enable_connectivity()
    for idV in collapsed_volumes:
        
        neighbors = [cell for cell in mesh.get_voxel_adjacent_voxels(idV) 
                        if getDirection(cell, idV, voxels)]
        
        neighbors.reverse()
        if(neighbors == []):
            continue
        else:
            direction = getDirection(idV, neighbors[0], voxels)
            collapseVolume(idV, neighbors[0], voxels, direction, collapsed_blocks, mesh)
            
    # for block in collapsed_blocks:
    #     active_cell_id, collapsed_cell_id = block
    #     expand_active_block(
    #         active_cell_id,
    #         collapsed_cell_id,
    #         voxels,
    #         vertices,
    #         threshold=0.05
    #     )

    voxels = removeCollapsedVolumes(collapsed_volumes, voxels)

    print("  Mapping Indexes")
    index = 0
    mapped_index = 0
    mapped_indexes = []

    for idx in collapsed_volumes:
        while(index != idx + 1):
            if(inactiveCells[mapped_index]):
                index += 1
            mapped_index += 1
        mapped_indexes.append(mapped_index - 1)

    for idx in mapped_indexes:
        inactiveCells[idx] = 0

    if(DEBUG):
        print(collapsed_volumes)
        print(collapsed_blocks)
    
    collapsed_mesh = pymesh.form_mesh(vertices, mesh.faces, voxels)
    
    #print("  Removing Isolated Vertices")
    #collapsed_mesh, _ = pymesh.remove_isolated_vertices(collapsed_mesh)
    
    dir_name = os.path.dirname(path_file)
    pymesh.save_mesh(path_file, collapsed_mesh, ascii=True)
    if(DEBUG): 
        highlightCollpsedVolumes(collapsed_volumes, mesh.voxels, mesh.vertices, mesh.num_voxels, os.path.join(dir_name, 'highlightedVolumes'))
    
    return inactiveCells

if __name__ == "__main__":
    if(len(sys.argv) < 2):
        print("Wrong numbers of arguments")
    path = sys.argv[1]
    main(path)
