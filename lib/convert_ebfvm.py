import pymesh
import sys
import os

def writeCFX5(voxels, points, num_voxels, type_of_mesh, nome):
        with open(nome + ".txt", "w") as newMesh:
            newMesh.write("1128683573\n")
            newMesh.write("Version number: 5.6D\n")
            
            if(type_of_mesh.upper() == 'HEXA'):
                newMesh.write(str(len(points)) + " 0 0 " + str(num_voxels) + " 0 1 1\n") #numberVertices, numberTetra, numberPrism, numberHexa, numberPiramid
            elif(type_of_mesh.upper() == 'TETRA'):
                newMesh.write(str(len(points)) + " " + str(num_voxels) + " 0 0" + " 0 1 1\n") #numberVertices, numberTetra, numberPrism, numberHexa, numberPiramid
            
            for point in points:
                newMesh.write(str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
            if(type_of_mesh.upper() == 'HEXA'):
                for voxel in voxels:
                    newMesh.write(f"{voxel[0] + 1} {voxel[1] + 1} {voxel[3] + 1} {voxel[2] + 1} {voxel[4] + 1} {voxel[5] + 1} {voxel[7] + 1} {voxel[6] + 1}")
                    newMesh.write("\n")
            elif(type_of_mesh.upper() == 'TETRA'):
                for voxel in voxels:
                    for idx in voxel:
                        newMesh.write(str(idx + 1) + " ")
                    newMesh.write("\n")
def main(file_path, type_of_mesh="TETRA"):
    mesh = pymesh.load_mesh(file_path)

    path, filename = os.path.split(file_path)
    filename = f"CFX5_{filename[:-4]}"
    file_path = os.path.join(path,filename)
    
    writeCFX5(mesh.voxels, mesh.vertices, mesh.num_voxels, type_of_mesh, file_path)
    print(f"Created file: {file_path}.txt")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Número incorreto de argumentos")
        os.exit(1)

    main(sys.argv[1])
