import pymesh
import sys
import os

DEBUG = False
THRESHOLD = 11

def areVolumesUnstable(vol1: float, vol2: float):
    if (vol1 > vol2):
        return (vol1 / vol2) >= THRESHOLD
    else:
        return (vol2 / vol1) >= THRESHOLD

def checkVolumes(checkedVolumes, tetraIdx:int, unstableNeighbours: list):
    checkedVolumes[tetraIdx] = True
    for neigh in unstableNeighbours:
        checkedVolumes[neigh] = True

def isVolumeChecked(checkedVolumes, tetraIdx: int):
    return checkedVolumes[tetraIdx] 

def getUnstableNeighbours(mesh, checkedVolumes, tetra_volumes,tetraIdx: int):
    unstableNeighbours = list()
    neighbours = mesh.get_voxel_adjacent_voxels(tetraIdx)
    for neigh in neighbours:
        if (not isVolumeChecked(checkedVolumes, neigh) ) and (areVolumesUnstable(tetra_volumes[tetraIdx], tetra_volumes[neigh]) ):
            unstableNeighbours.append(neigh)
    return unstableNeighbours

def printTetraInfo(mesh, tetra_volumes, tetraIdx: int):
    print(f"{tetraIdx} {mesh.elements[tetraIdx]} {tetra_volumes[tetraIdx] }\n")

def printTetraInfoWithNeighbours(mesh, tetra_volumes, tetraIdx: int, unstableNeighbours: list):
    print('ID     VERTICES        VOLUME\n')
    printTetraInfo(mesh, tetra_volumes, tetraIdx)
    for neigh in unstableNeighbours:
        printTetraInfo(mesh, tetra_volumes, neigh)
    print("==================================")

def printTable(mesh, checkedVolumes, tetra_volumes):
    unstableRegions = 0
    for tetraIdx in range(mesh.num_elements):
        unstableNeighbours = getUnstableNeighbours(mesh, checkedVolumes, tetra_volumes, tetraIdx)
        if (len(unstableNeighbours) > 0):
            printTetraInfoWithNeighbours(mesh, tetra_volumes, tetraIdx, unstableNeighbours)
            checkVolumes(checkedVolumes, tetraIdx, unstableNeighbours)

def unstableVolumes(mesh, checkedVolumes, tetra_volumes, THRESHOLD, DEBUG=False):
    if(DEBUG): print(f"Threshold: {THRESHOLD}")
    unstableRegions = 0
    for tetraIdx in range(mesh.num_elements):
        unstableNeighbours = getUnstableNeighbours(mesh, checkedVolumes, tetra_volumes, tetraIdx)
        if (len(unstableNeighbours) > 0):
            checkVolumes(checkedVolumes, tetraIdx, unstableNeighbours)
            if(DEBUG): 
                print(f"Found Unstable volumes for Threshold {THRESHOLD}")
                printTetraInfoWithNeighbours(tetraIdx, unstableNeighbours)
            return True
    return False

def writeMSH(checkedVolumes, voxels, points, numTetra, nome):
    with open(nome + ".msh", "w") as newFile:
        newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        newFile.write("$Nodes\n" + str(len(points)) + "\n")
        for (index, point) in enumerate(points):
            newFile.write(str(
                index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
        newFile.write("$EndNodes\n")
        newFile.write("$Elements\n" + str(numTetra) + "\n")

        tagNumber = 0

        for (idV, voxel) in enumerate(voxels):
            if (isVolumeChecked(checkedVolumes, idV)):
                tagNumber = idV + 1
            else:
                tagNumber = 0
            newFile.write(
                str(idV + 1) + f" 4 2 0 {tagNumber} " +
                str(voxel[0] + 1) + " " + str(voxel[1] + 1) + " "
                + str(voxel[2] + 1) + " " + str(voxel[3] + 1) + "\n"
            )
        newFile.write("$EndElements")

def loadData(filename):
    mesh = pymesh.load_mesh(filename)
    mesh.add_attribute("voxel_volume")
    tetra_volumes = mesh.get_attribute("voxel_volume")
    mesh.enable_connectivity()
    checkedVolumes = {tIdx: False for tIdx in range(mesh.num_elements)}
    
    return mesh, checkedVolumes, tetra_volumes

def findThreshold(filename):
    global THRESHOLD

    mesh, checkedVolumes, tetra_volumes = loadData(filename)

    max_threshold = 100
    min_threshold = 2
    while(min_threshold <= max_threshold):
        THRESHOLD = min_threshold + (max_threshold - min_threshold)//2

        if(unstableVolumes(mesh, checkedVolumes, tetra_volumes, THRESHOLD)):
            min_threshold = THRESHOLD + 1
        else:
            max_threshold = THRESHOLD - 1
    
    if(THRESHOLD >= 100):
        print(f"Mesh has unstable volumes of ratio {THRESHOLD}")
        unstableVolumes(unstableVolumes(mesh, checkedVolumes, tetra_volumes, THRESHOLD, DEBUG=True))
        writeMSH(chekedVolumes, mesh.voxels, mesh.vertices, mesh.num_voxels, "unstableRegions")
    else:
        print(f"No unstable volumes for threshold: {THRESHOLD}")

def main(path):
    mesh, checkedVolumes, tetra_volumes = loadData(path)
    printTable(mesh, checkedVolumes, tetra_volumes)
    output_path, _ = os.path.split(path)
    output_path = os.path.join(output_path,"unstableRegions")
    writeMSH(checkedVolumes, mesh.voxels, mesh.vertices, mesh.num_voxels, output_path)

if __name__ == "__main__":    
    path = input("Type path to mesh file: ")
    THRESHOLD = int(input("Type the threshold for the volumes: "))
    main(path)