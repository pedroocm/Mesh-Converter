import pymesh

def writeCFX5(voxels, points, numTetra, nome):
	with open(nome + ".txt", "w") as newMesh:
		newMesh.write("1128683573\n")
		newMesh.write("Version number: 5.6D\n")
		newMesh.write(str(len(points)) + " " + str(numTetra) + " 0 0" + " 0 1 1\n") #numberVertices, numberTetra, numberPrism, numberHexa, numberPiramid
		for point in points:
			newMesh.write(str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
		for voxel in voxels:
			for idx in voxel:
				newMesh.write(str(idx + 1) + " ")
			newMesh.write("\n")


mesh = pymesh.load_mesh("uniao_tetra.msh")

mesh = pymesh.resolve_self_intersection(mesh)

tet = pymesh.tetgen()
tet.points = mesh.vertices
tet.triangles = mesh.faces
tet.verbosity = 4
tet.split_boundary = False
tet.run()
new_mesh = tet.mesh



pymesh.save_mesh("second_tetrahedralization.msh", new_mesh, ascii=True)
writeCFX5(new_mesh.voxels, new_mesh.vertices, new_mesh.num_voxels, "second_tetrahedralization")
