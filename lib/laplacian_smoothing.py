import numpy as np
import pymesh 
import sys
from os import path

DEBUG = False

def computeWeightBetweenVertices(vert1: int, vert2: int):
    return 1



## Expects mesh to have connectivity enabled
def computeAverageNeighbourVertices(vert: int, mesh: list, vertices: list):
    phi = 0.5
    sumVertices = 0
    sumWeights = 0

    for neighbourIndex in mesh.get_vertex_adjacent_vertices(vert):
        sumWeights = sumWeights + computeWeightBetweenVertices(vert, neighbourIndex)
        sumVertices = sumVertices + computeWeightBetweenVertices(vert, neighbourIndex) * (vertices[neighbourIndex] - vertices[vert]) 
    return phi * (sumVertices / sumWeights)

def isRatioTooBig(vol1: float, vol2: float, limit: int):
    return min(vol1, vol2) == 0 or max(vol1, vol2)/min(vol1, vol2) > limit

# Expects mesh to have connectivity enabled. 
def isQualityBetter(vert: int, mesh: pymesh.Mesh, new_mesh: pymesh.Mesh):
    vertAdjacentTetra = mesh.get_vertex_adjacent_voxels(vert)
    new_mesh.add_attribute("voxel_volume")
    volumes = new_mesh.get_attribute("voxel_volume") 
    for adjacent in vertAdjacentTetra:
        if volumes[adjacent] < 0: 
            return False
        for adjacentVolume in mesh.get_voxel_adjacent_voxels(adjacent):
            if isRatioTooBig(volumes[adjacent], volumes[adjacentVolume], 22):
                return False
    return True

def getSharedFaces(tetra: int, tetraNeighbours: list, mesh: pymesh.Mesh) -> list:
    tetraVertices = set(mesh.elements[tetra])
    sharedFaces = list()
    for neighbour in tetraNeighbours:
        neighbourVertices = set(mesh.elements[neighbour])
        sharedFace = tetraVertices.intersection(neighbourVertices)
        sharedFaces.append(sharedFace)
    return sharedFaces

def isLockedVertex(vertIndex: int, lockedVertices: dict) -> bool:
    try:
        return lockedVertices[vertIndex]
    except KeyError:
        return False

def isVertexExternal(vertIndex: int, mesh: pymesh.Mesh, vertices: list)-> bool:
    
    yMin, zMin, yMax, zMax = mesh.bbox[0][1], mesh.bbox[0][2], mesh.bbox[1][1], mesh.bbox[1][2] 
    vertex = vertices[vertIndex]

    if True in np.isclose([vertex[1], vertex[2], vertex[1], vertex[2]], [yMin, zMin, yMax, zMax]):
        return True

    vertexAdjacentTetrahedra = mesh.get_vertex_adjacent_voxels(vertIndex)

    for tetra in vertexAdjacentTetrahedra:
        numOfNeighbours = len(mesh.get_voxel_adjacent_voxels(tetra))
        if numOfNeighbours <= 2:
            return True
        elif numOfNeighbours == 3:
            neighboursTetra = mesh.get_voxel_adjacent_voxels(tetra)
            sharedFaces = getSharedFaces(tetra, neighboursTetra, mesh)
            for face in sharedFaces:
                if vertIndex not in face:
                    return True
        else:
            continue

    return False

def checkInputArguments():
    if len(sys.argv) != 2:
        print("Incorrect number of arguments. Usage: laplacian_smoothing name_of_mesh")
        sys.exit(1)

    if not path.exists(sys.argv[1]):
        print(f" File {sys.argv[1]} does not exist! ")
        sys.exit(1)

    return sys.argv[1]

def laplacian_smoothing(mesh: pymesh.Mesh, locked_vertices={}):
    mesh.enable_connectivity()
    
    changedVertices = np.copy(mesh.vertices)
    for i in range(5):
        for vertexIndex in range(mesh.num_vertices):
            if isVertexExternal(vertexIndex, mesh, changedVertices):
                continue
            elif isLockedVertex(vertexIndex, locked_vertices):
                continue
            averageTerm = computeAverageNeighbourVertices(vertexIndex, mesh, changedVertices)
            oldValue = changedVertices[vertexIndex]
            changedVertices[vertexIndex] = oldValue + averageTerm
            new_mesh = pymesh.form_mesh(changedVertices, mesh.faces, mesh.voxels)
            if not isQualityBetter(vertexIndex, mesh, new_mesh):
                if(DEBUG): print("Quality is not better")
                changedVertices[vertexIndex] = oldValue

    mesh = pymesh.form_mesh(changedVertices, mesh.faces, mesh.voxels)

def main(path_to_mesh: str, locked_vertices={}):
    mesh = pymesh.load_mesh(path_to_mesh)
    laplacian_smoothing(mesh, locked_vertices)
    output_path = f"{path_to_mesh[:-4]}_laplacian.msh"
    pymesh.save_mesh(output_path, mesh, ascii=True)
    print(f"Created File: {output_path}")

if __name__ == "__main__":
    mesh_name = checkInputArguments()
    main(mesh_name)
