import pymesh
from enum import Enum
import numpy
import sys
from lib.transformData import map_indexes
from lib.extractData2 import getInactiveCells
import faulthandler

faulthandler.enable()
DEBUG = False

class Direction(Enum):
	NORTHWEST = 0
	NORTHEAST = 1
	SOUTHWEST = 2
	SOUTHEAST = 3

directionEdge = {
	Direction.NORTHWEST: (3, 0),
	Direction.NORTHEAST: (2, 1),
	Direction.SOUTHWEST: (7, 4),
	Direction.SOUTHEAST: (6, 5)
}

class fDirection(Enum):
	NORTH = 0
	SOUTH = 1
	EAST = 2
	WEST = 3

facesEdges = {
	fDirection.NORTH: (Direction.NORTHWEST, Direction.NORTHEAST),
	fDirection.SOUTH: (Direction.SOUTHWEST, Direction.SOUTHEAST),
	fDirection.WEST: (Direction.NORTHWEST, Direction.SOUTHWEST),
	fDirection.EAST: (Direction.NORTHEAST, Direction.SOUTHEAST)
}


neighborEdge = {
	Direction.NORTHWEST: {fDirection.WEST: Direction.SOUTHWEST, fDirection.NORTH: Direction.NORTHEAST},
	Direction.NORTHEAST: {fDirection.EAST: Direction.SOUTHEAST, fDirection.NORTH: Direction.NORTHWEST},
	Direction.SOUTHWEST: {fDirection.WEST: Direction.NORTHWEST, fDirection.SOUTH: Direction.SOUTHEAST},
	Direction.SOUTHEAST: {fDirection.EAST: Direction.NORTHEAST, fDirection.SOUTH: Direction.SOUTHWEST}
}

def is_collapsed(cellIdx, tol):
	return avg_depth(cellIdx) < tol

def avg_depth(cellIdx):
	ne_depth = getCellDepth(cellIdx, directionEdge[Direction.NORTHEAST])
	nw_depth = getCellDepth(cellIdx, directionEdge[Direction.NORTHWEST])
	se_depth = getCellDepth(cellIdx, directionEdge[Direction.SOUTHEAST])
	sw_depth = getCellDepth(cellIdx, directionEdge[Direction.SOUTHWEST])
	
	return (ne_depth + nw_depth + se_depth + sw_depth) / 4

def getCellDepth(cellIdx, direction):
	return abs( getVertex(cellIdx, direction[1])[2] - getVertex(cellIdx, direction[0])[2] )

def is_face_collapsed(face):
	return (face[1][2] - face[0][2] == 0) and (face[3][2] - face[2][2] == 0)

def generate_mesh(tetraedros):
	points = dict()
	index = 0
	faces = list()

	for triangles in tetraedros:
		for triangle in triangles.faces:
			face = list()
			for vertex_index in triangle:
				vertex = triangles.vertices[vertex_index]
				if tuple(vertex) not in points:
					points[tuple(vertex)] = index
					index += 1
				face.append(points[tuple(vertex)])
			faces.append(face)

	vertices = [numpy.asarray(point) for point in points]
	vertices = numpy.asarray(vertices)

	faces = numpy.asarray(faces)

	return pymesh.form_mesh(vertices, faces)


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

def writeMSH(voxels, points, numTetra, nome):
	with open(nome + ".msh", "w") as newFile:
		newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
		newFile.write("$Nodes\n" + str(len(points)) + "\n")
		for (index, point) in enumerate(points):
			newFile.write(str(index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
		newFile.write("$EndNodes\n")
		newFile.write("$Elements\n" + str(numTetra) + "\n")

		for (idV, voxel) in enumerate(voxels):
			newFile.write(
				str(idV + 1) + f" 4 2 0 {idV + 1} " + str(voxel[0] + 1) + " " + str(voxel[1] + 1) + " " 
				+ str(voxel[2] + 1) + " " + str(voxel[3] + 1)  + "\n"
				)

		newFile.write("$EndElements")

def tetraConversor(meshFilename: str, numCellsPerAxis: list, inactiveCells: list, optimization=False, refinement=0.125):
	print(f"meshFileName: {meshFilename}")

	mesh = pymesh.load_mesh(meshFilename)

	nz = numCellsPerAxis[2]
	ny = numCellsPerAxis[1]
	nx = numCellsPerAxis[0]
	corrector = 0
	columns = {}

	if(optimization):
		for vertice in mesh.vertices:
			x, y, z = vertice
			if((x, y) not in columns):
				columns[(x, y)] = [vertice]
			else:
				columns[(x, y)].append(vertice)
	
	def getNeighborColumn(col_index, direction, faceDir):
		col = []
				
		edge_points_idx = directionEdge[neighborEdge[direction][faceDir]]

		first_col_index = -1
		last_col_index = -1
		
		for i in range(col_index, col_index + nz):
			if(inactiveCells[i]):
				first_col_index = i
				break
		
		if(first_col_index == -1):
			return col
		
		for i in range(col_index + nz - 1, first_col_index + 1, -1):
			if(inactiveCells[i]):
				last_col_index = i
				break

		if(last_col_index == -1):
			new_end_index = first_col_index

		first_element_current_column = mesh.elements[map_indexes(inactiveCells, first_col_index)]
		last_element_current_column = mesh.elements[map_indexes(inactiveCells, last_col_index)]

		first_vertex_id = first_element_current_column[edge_points_idx[0]]
		last_vertex_id = last_element_current_column[edge_points_idx[1]]

		first_point = mesh.vertices[ first_vertex_id ] 
		last_point = mesh.vertices[ last_vertex_id ]

		if(optimization):
			return columns[(first_point[0], first_point[1])]
		
		for vertex_id in columns:
			vertex_coords = mesh.vertices[ vertex_id ]
			if(pymesh.is_colinear(vertex_coords, first_point, last_point)):
				return columns[vertex_id]
		
		for vertex_coords in mesh.vertices:
			if(pymesh.is_colinear(vertex_coords, first_point, last_point)):
				col.append(vertex_coords)
		
		columns[first_vertex_id] = col
		return col

	def getNeighborColIdx(xy_index):
		return [xy_index - nz, xy_index + nz, xy_index - ny * nz, xy_index + ny * nz] #N, S, W, E

	def hasEdge(xy_index, direction, curCoord):
		
		x, y = curCoord
		idxArray = getNeighborColIdx(xy_index)
		hasNorth = (y - 1) >= 0
		hasWest = (x - 1) >= 0
		hasEast = (x + 1) < nx
		hasSouth = (y + 1) < ny

		truthValuesDirection = {
			Direction.NORTHWEST: (hasWest, hasNorth, 2, 0),
			Direction.NORTHEAST: (hasEast, hasNorth, 3, 0),
			Direction.SOUTHWEST: (hasWest, hasSouth, 2, 1),
			Direction.SOUTHEAST: (hasEast, hasSouth, 3, 1)
		}   

		indexes = [-1, -1]

		truthValue = truthValuesDirection[direction]
		if truthValue[0]: indexes[0] = idxArray[truthValue[2]]
		if truthValue[1]: indexes[1] = idxArray[truthValue[3]]

		return indexes

	def getVertex(cellIdx, vertexIdx, corrector=-1):
		return mesh.vertices[ mesh.elements[ map_indexes(inactiveCells, cellIdx, corrector) ][ vertexIdx ] ]


	celulas = {}

	print("Extrating points coords")
	for x in range(nx):
		for y in range(ny):
			xy_index = y*nz + x*nz*ny
			dictCols = {
				Direction.NORTHWEST: {fDirection.NORTH:[], fDirection.WEST: []},
				Direction.NORTHEAST: {fDirection.NORTH:[], fDirection.EAST: []},
				Direction.SOUTHWEST: {fDirection.SOUTH:[], fDirection.WEST: []}, 
				Direction.SOUTHEAST: {fDirection.SOUTH:[], fDirection.EAST: []}
			}

			for d in dictCols:
				idx = hasEdge(xy_index, d, [x, y])
				for (i, fd) in enumerate(dictCols[d]):
					if (idx[i] != -1): dictCols[d][fd] = getNeighborColumn(idx[i], d, fd)
			for z in range(nz):
				cell_id = xy_index + z
				
				if(not inactiveCells[cell_id]):
					corrector += 1
					continue
				
				cellFaces = []
				for face in facesEdges:
					facePoints = []
					for edge in facesEdges[face]:
						(topVertex, botVertex) = directionEdge[edge]
						
						top_vertex_coords = getVertex(cell_id, topVertex, corrector)
						bot_vertex_coords = getVertex(cell_id, botVertex, corrector)
						
						facePoints.append(top_vertex_coords)
						facePoints.append(bot_vertex_coords)

						for fDir in dictCols[edge]:
							for point in dictCols[edge][fDir]:
								if top_vertex_coords[2] < point[2] < bot_vertex_coords[2]:
									if ((point[2] - top_vertex_coords[2]) >= 0.001) and ((bot_vertex_coords[2] - point[2]) >= 0.001):
										facePoints.append(point)
					
					cellFaces.append(facePoints)
				
				cellFaces.append([getVertex(cell_id, directionEdge[Direction.NORTHEAST][0], corrector), getVertex(cell_id, directionEdge[Direction.NORTHWEST][0], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHEAST][0], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHWEST][0], corrector) ])
				cellFaces.append([getVertex(cell_id, directionEdge[Direction.NORTHEAST][1], corrector), getVertex(cell_id, directionEdge[Direction.NORTHWEST][1], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHEAST][1], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHWEST][1], corrector) ])

				celulas[cell_id] = cellFaces

	tetraedros = []
	def tri_tetra(cell, index):
		tri_meshs = []

		for i, face in enumerate(cell):
			tri = pymesh.triangle()
			tri.points = face
			tri.max_num_steiner_points = 0
			tri.verbosity = 0

			if(i > len(cell) - 3):
				faces = numpy.array([[3, 1, 0], [0, 2, 3]])
				triangulation = pymesh.form_mesh(face, faces)
			elif(len(face) == 4 and numpy.array_equal(face[0], face[1]) and numpy.array_equal(face[2], face[3])):
				continue
			else:
				tri.run()
				triangulation = tri.mesh
			tri_meshs.append(triangulation)
		return tri_meshs

	del columns
	
	print("Generating triangles")
	tetraedros = (tri_face for cel in celulas 
					for tri_face in tri_tetra(celulas[cel], cel))

	print("Generating triangular mesh")
	uniao = generate_mesh(tetraedros)
	
	print(  "Removing duplicated faces and vertices")
	uniao, info = pymesh.remove_duplicated_vertices(uniao, 0.1)
	uniao, info = pymesh.remove_duplicated_faces(uniao)
	print("  Removing degenerated triangles")
	uniao, info = pymesh.remove_degenerated_triangles(uniao, num_iterations=1500)
	print("  Collapsing short edges")
	uniao, info = pymesh.collapse_short_edges(uniao, rel_threshold=refinement)
	if(DEBUG): pymesh.save_mesh("malha_triangular.msh", uniao)


	print("  Removing self intersections")
	while( len(pymesh.detect_self_intersection(uniao)) != 0):
		if(DEBUG): print(f' Number of self intersections before = {len(pymesh.detect_self_intersection(uniao))}')
		uniao = pymesh.resolve_self_intersection(uniao)
		uniao, info = pymesh.remove_duplicated_vertices(uniao, 0.1) #no caso anterior, esta linha estava abaixo
		uniao, info = pymesh.remove_duplicated_faces(uniao)
		if(DEBUG): print(f' Number of self intersections after = {len(pymesh.detect_self_intersection(uniao))}')

	
	if(DEBUG): print(f' Number of self intersections after = {len(pymesh.detect_self_intersection(uniao))}')
	uniao = pymesh.resolve_self_intersection(uniao)
	if(DEBUG): pymesh.save_mesh("triang_apos_remove.msh", uniao)

	print("Generating tetra mesh")
	tetgen = pymesh.tetgen()
	tetgen.points = uniao.vertices # Input points.
	tetgen.triangles = uniao.faces # Input triangles
	#tetgen.max_tet_volume = 0.01;
	tetgen.verbosity = 4 if DEBUG else 0
	tetgen.split_boundary = False
	tetgen.coarsening = True
	#tetgen.optimization_level = 10
	#tetgen.max_num_steiner_points = 0
	tetgen.run() # Execute tetgen
	uniao = tetgen.mesh

	pymesh.save_mesh(f"{meshFilename[:-4]}_tetra.msh", uniao, ascii=True)
	if(DEBUG): writeCFX5(uniao.voxels, uniao.vertices, uniao.num_voxels, "uniao_tetra")
	if(DEBUG): writeMSH(uniao.voxels, uniao.vertices, uniao.num_voxels, "uniao_tetrav")

def noHNTetraConversor(meshFilename, numCellsPerAxis, inactiveCells):
	print(f"meshFileName: {meshFilename}")

	mesh = pymesh.load_mesh(meshFilename)


	'''
	nz = 3
	ny = 1
	nx = 2
	'''
	nz = numCellsPerAxis[2] #4
	ny = numCellsPerAxis[1] #17
	nx = numCellsPerAxis[0] #27
	corrector = 0

	def getVertex(cellIdx, vertexIdx, corrector=-1):
		return mesh.vertices[ mesh.elements[ map_indexes(inactiveCells, cellIdx, corrector) ][ vertexIdx ] ]

	celulas = {}

	print("Extrating points coords")
	for x in range(nx):
		for y in range(ny):
			xy_index = y*nz + x*nz*ny

			for z in range(nz):
				cell_id = xy_index + z
				
				if(not inactiveCells[cell_id]):
					corrector += 1
					continue

				cellFaces = []
				for face in facesEdges:
					facePoints = []
					for edge in facesEdges[face]:
						(topVertex, botVertex) = directionEdge[edge]
						facePoints.append(getVertex(cell_id, topVertex, corrector))
						facePoints.append(getVertex(cell_id, botVertex, corrector))

					cellFaces.append(facePoints)
				
				cellFaces.append([getVertex(cell_id, directionEdge[Direction.NORTHEAST][0], corrector), getVertex(cell_id, directionEdge[Direction.NORTHWEST][0], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHEAST][0], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHWEST][0], corrector) ])
				cellFaces.append([getVertex(cell_id, directionEdge[Direction.NORTHEAST][1], corrector), getVertex(cell_id, directionEdge[Direction.NORTHWEST][1], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHEAST][1], corrector), getVertex(cell_id, directionEdge[Direction.SOUTHWEST][1], corrector) ])

				celulas[cell_id] = cellFaces


	tetraedros = []

	def tri_tetra(cell, index):
		tri_meshs = []

		for i, face in enumerate(cell):
			vertices = face
			tri = pymesh.triangle()
			tri.points = vertices
			tri.max_num_steiner_points = 0
			tri.verbosity = 0
			if(i > len(cell) - 3):
				faces = numpy.array([[3, 1, 0], [0, 2, 3]])
				triangulation = pymesh.form_mesh(face, faces)
				tri_meshs.append(triangulation)
			else: 
				tri.run()
				triangulation = tri.mesh
				tri_meshs.append(triangulation)
			
		#uniao = pymesh.merge_meshes(tri_meshs)
	
		return tri_meshs #uniao

	print("Generating triangles")
	for cel in celulas:
		for tri_face in tri_tetra(celulas[cel], cel):
			tetraedros.append(tri_face)

	print("Generating triangular mesh")	
	uniao = generate_mesh(tetraedros)
	
	uniao, info = pymesh.remove_duplicated_vertices(uniao, 0.1) #no caso anterior, esta linha estava abaixo
	uniao, info = pymesh.remove_duplicated_faces(uniao)
	uniao, info = pymesh.remove_degenerated_triangles(uniao, num_iterations=1500)
	uniao, info = pymesh.collapse_short_edges(uniao, rel_threshold=0.125)
	if(DEBUG): pymesh.save_mesh(f"{meshFilename[:-4]}malha_triangular.msh", uniao)

	print("Removing self intersections")
	while( len(pymesh.detect_self_intersection(uniao)) != 0):
		if(DEBUG): print(f' Number of self intersections before = {len(pymesh.detect_self_intersection(uniao))}')
		uniao = pymesh.resolve_self_intersection(uniao)
		uniao, info = pymesh.remove_duplicated_vertices(uniao, 0.1) #no caso anterior, esta linha estava abaixo
		uniao, info = pymesh.remove_duplicated_faces(uniao)
		if(DEBUG): print(f' Number of self intersections after = {len(pymesh.detect_self_intersection(uniao))}')

	
	if(DEBUG): print(f' Number of self intersections after = {len(pymesh.detect_self_intersection(uniao))}')
	uniao = pymesh.resolve_self_intersection(uniao)
	if(DEBUG): pymesh.save_mesh(f"{meshFilename[:-4]}_triang_apos_remove.msh", uniao)

	#pymesh.tetrahedralize(uniao, 100, engine='cgal')
	print("Generating tetrahedralized mesh")
	tetgen = pymesh.tetgen()
	tetgen.points = uniao.vertices # Input points.
	tetgen.triangles = uniao.faces # Input triangles
	#tetgen.max_tet_volume = 0.01;
	tetgen.verbosity = 4 if DEBUG else 0
	tetgen.split_boundary = False
	tetgen.coarsening = True
	# tetgen.coplanar_tolerance = 10**-6
	#tetgen.optimization_level = 10
	#tetgen.max_num_steiner_points = 0
	tetgen.run() # Execute tetgen
	uniao = tetgen.mesh

	pymesh.save_mesh(f"{meshFilename[:-4]}_tetra.msh", uniao, ascii=True)
	if(DEBUG): writeCFX5(uniao.voxels, uniao.vertices, uniao.num_voxels, f"{meshFilename[:-4]}_tetra.txt")
	if(DEBUG): writeMSH(uniao.voxels, uniao.vertices, uniao.num_voxels, "uniao_tetrav")

def debugTetralization(path_to_mesh):

	uniao = pymesh.load_mesh(path_to_mesh)

	print("Generating tetrahedralized mesh")
	tetgen = pymesh.tetgen();
	tetgen.points = uniao.vertices; # Input points.
	tetgen.triangles = uniao.faces; # Input triangles
	#tetgen.max_tet_volume = 0.01;
	tetgen.verbosity = 4 if DEBUG else 0
	tetgen.split_boundary = False
	tetgen.coarsening = True
	tetgen.coplanar_tolerance = 10**-4
	#tetgen.optimization_level = 10
	tetgen.max_num_steiner_points = 0
	error = True
	tetgen.run()
	uniao = tetgen.mesh
	
	
	'''
	tetgen = pymesh.tetgen();
	tetgen.points = uniao.vertices; # Input points.
	tetgen.triangles = uniao.faces; # Input triangles
	tetgen.tetrahedra = uniao.voxels;
	#tetgen.max_tet_volume = 0.01;
	tetgen.verbosity = 4
	tetgen.split_boundary = False
	tetgen.coarsening = True
	#tetgen.optimization_level = 10
	#tetgen.max_num_steiner_points = 0
	tetgen.run(); # Execute tetgen
	uniao = tetgen.mesh;
	'''

	pymesh.save_mesh(f"{meshFilename[:-4]}_tetra.msh", uniao, ascii=True)
