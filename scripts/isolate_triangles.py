# coding: utf-8

import pymesh

mesh = pymesh.load_mesh("triang_apos_remove.msh")
mesh.enable_connectivity()

face_level = dict()

def isolateTriangles(maxDepth, curFace):
    depth = 0
    current = [curFace]
    neighbours = []
    faces = []
    seen = []
    while depth < maxDepth:
        for f in current:
            if f not in seen:
                faces.append(mesh.elements[f])
                seen.append(f)
                face_level[f] = depth + 1
                for n in mesh.get_face_adjacent_faces(f):
                    neighbours.append(n)
            
        current = neighbours
        neighbours = []
        depth += 1
    return faces

def writeMSH(voxels, points, numTetra, nome):
    with open(nome + ".msh", "w") as newFile:
        newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        newFile.write("$Nodes\n" + str(len(points)) + "\n")
        for (index, point) in enumerate(points):
            newFile.write(str(index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
        newFile.write("$EndNodes\n")
        newFile.write("$Elements\n" + str(numTetra) + "\n")

        for (idV, voxel) in enumerate(voxels):
            tag = 0
            if idV in face_level:
                tag = face_level[idV]
            newFile.write(
                str(idV + 1) + f" 2 2 0 {tag} " + str(voxel[0] + 1) + " " + str(voxel[1] + 1) + " " 
                + str(voxel[2] + 1) + " " + "\n"
                )


        newFile.write("$EndElements")
    
new_faces = isolateTriangles(10, 9501)
import numpy
new_faces = numpy.asarray(new_faces)
new_mesh = pymesh.form_mesh(mesh.vertices, new_faces)
pymesh.save_mesh("isolated_triangles.msh", new_mesh, ascii=True)
writeMSH(mesh.faces, mesh.vertices, mesh.num_elements, "isolated_triangles_colorized")


