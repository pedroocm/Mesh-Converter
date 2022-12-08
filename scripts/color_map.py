import pymesh


def get_tag_number(volume: float):
    for i in range(6):
        if (10 ** i) <= volume < (10 ** (i+1)):
            return i


def get_tag_array(volumes: list):
    tag_array = []

    for vol in volumes:
        tag_array.append(get_tag_number(vol))

    return tag_array


mesh = pymesh.load_mesh("uniao_tetrav.msh")
mesh.add_attribute("voxel_volume")

tag_array = get_tag_array(mesh.get_attribute("voxel_volume"))


def writeMSH(voxels, points, numTetra, nome):
    with open(nome + ".msh", "w") as newFile:
        newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        newFile.write("$Nodes\n" + str(len(points)) + "\n")
        for (index, point) in enumerate(points):
            newFile.write(str(
                index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
        newFile.write("$EndNodes\n")
        newFile.write("$Elements\n" + str(numTetra) + "\n")

        for (idV, voxel) in enumerate(voxels):
            newFile.write(
                str(idV + 1) + f" 4 2 0 {tag_array[idV]} " +
                str(voxel[0] + 1) + " " + str(voxel[1] + 1) + " "
                + str(voxel[2] + 1) + " " + str(voxel[3] + 1) + "\n"
            )
        newFile.write("$EndElements")

writeMSH(mesh.voxels, mesh.vertices, mesh.num_voxels, "uniao_colormap")
