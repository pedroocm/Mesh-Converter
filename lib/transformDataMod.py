# This function transforms the COORD and ZCORN lists obtained by the extracData module, and transforms them into a tridimensional list, in which each element is a Cell. 
#import rayIntersection as ri
#from rayIntersection import *
from os import nice
from lib.generic_binary import binarySearch
import pymesh
from time import time
import re
from sys import exit


class Vector:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	@classmethod
	def fromList(cls, lista):
		return cls(lista[0], lista[1], lista[2])
	def __add__(self, other):
		x = self.x + other.x
		y = self.y + other.y
		z = self.z + other.z
		return Vector(x, y, z)
	def __sub__(self, other):
		x = self.x - other.x
		y = self.y - other.y
		z = self.z - other.z
		return Vector(x, y, z)
	def __mul__(self, scalar):
		x = self.x * scalar
		y = self.y * scalar
		z = self.z * scalar
		return Vector(x, y, z)
	def __str__(self):
		return "({0}, {1}, {2})".format(self.x, self.y, self.z)

def map_indexes(inactiveCells, index, offset=-1):
	
	if(offset != -1):
		return index - offset
	
	return sum(inactiveCells[:index])


def getPoint(coord, zCorn):
	dirVector = Vector(coord[3], coord[4], coord[5]) - Vector(coord[0], coord[1], coord[2])

	scaleFactor = 0
	if dirVector.z == 0:
		#print("[Error - getPoint] z-coordinates cannot be equal: ({0}, {1}, {2}), ({3}, {4}, {5})".format(coord[0], coord[1], coord[2], coord[3], coord[4], coord[5]))
		#print(f"ZCorn: {zCorn}")
		scaleFactor = 0
		#exit(1)
	else:
		scaleFactor = (zCorn - coord[2])/dirVector.z

	xCorn = scaleFactor * dirVector.x + coord[0]
	yCorn = scaleFactor * dirVector.y + coord[1]

	return (xCorn, yCorn, zCorn)



def hexLaplacianOld(matrixCells, arrayOfPoints, numberOfCellValues):
	(nx, ny, nz) = numberOfCellValues
	phi = 0.5
	for i in range(nx):
		for j in range(ny):
			for k in range(1, nz):
				cell = matrixCells[i][j][k]
				topCell = matrixCells[i][j][k-1]
				for x in range(4):
					topVertex = arrayOfPoints[cell[x]]
					botVertex = arrayOfPoints[cell[x+4]]
					if topVertex == botVertex and topVertex == arrayOfPoints[topCell[x+4]]:
						topNeighbor = arrayOfPoints[topCell[x]]
						newPoint = Vector.fromList(botVertex)*(1-phi) + Vector.fromList(topNeighbor)*phi
						arrayOfPoints[cell[x]] = (newPoint.x, newPoint.y, newPoint.z)
	return arrayOfPoints

def hasSmallThickness(cell, points, tol):
	topMean = 0
	botMean = 0
	for i in range(4):
		topMean += points[cell[i]][2]
		botMean += points[cell[i+4]][2]
	#print("top:", topMean, "   bot =", botMean)
	#print("abs: ", abs((topMean - botMean)/4) <= tol)
	return abs((topMean - botMean)/4) <= tol



def checkEdgeHeight(cell, column, arrayOfPoints, pinchTol):
	height = abs(arrayOfPoints[cell[column+4]][2] - arrayOfPoints[cell[column]][2])
	return height < pinchTol

def needsLaplacian(cell, cellIndex, column, arrayOfPoints, inactiveCells, pinchArray, pinchTol):
	if not checkEdgeHeight(cell, column, arrayOfPoints, pinchTol):
		return False
	if pinchArray[cellIndex] == 0:
		return True
	if inactiveCells[cellIndex] == 0:
		return False
	return True
	



def vertexHasUpperConnectivity(i, j, k, columnIndex, matrixCells, arrayOfPoints):
	if k == 0:
		return True
	cell = matrixCells[i][j][k]
	upperCell = matrixCells[i][j][k-1]
	point = arrayOfPoints[cell[columnIndex]]
	pointAbove = arrayOfPoints[upperCell[columnIndex+4]]
	if point != pointAbove:
		return False
	return True




def getUpperPointIndex(i, j, k, index, pointIndex, cellIndex, matrixCells, arrayOfPoints, inactiveCells):
	if index >= 4:
		return matrixCells[i][j][k][index-4]
	if k == 0 or inactiveCells[cellIndex - 1] == 0:
		return pointIndex
	upperCell = matrixCells[i][j][k-1]
	if pointIndex != upperCell[index+4]:
		return pointIndex

	return upperCell[index]

def getPointAndUpperPoint(i, j, k, matrixCells, arrayOfPoints, columnIndex):
	point = arrayOfPoints[matrixCells[i][j][k][columnIndex]]
	if columnIndex < 4:
		upperPoint = arrayOfPoints[matrixCells[i][j][k-1][columnIndex+4]]
	else:
		upperPoint = arrayOfPoints[matrixCells[i][j][k][columnIndex-4]]
	return (point, upperPoint)
	
def vertexNeedsLaplacian(i, j, k, matrixCells, arrayOfPoints, columnIndex, tol):
	(point, upperPoint) = getPointAndUpperPoint(i, j, k, matrixCells, arrayOfPoints, columnIndex+4)
	return abs(point[2] - upperPoint[2]) <= tol

def cellHasUpperConnectivity(i, j, k, matrixCells, arrayOfPoints):
	if k == 0:
		return True
	for columnIndex in range(4):
		(point, pointAbove) = getPointAndUpperPoint(i, j, k, matrixCells, arrayOfPoints, columnIndex)
		if point != pointAbove:
			return False
	return True

def hexLaplacian(matrixCells, arrayOfPoints, numberOfCellValues, inactiveCells, phi, tol = 0.1):
	(nx, ny, nz) = numberOfCellValues
	for i in range(nx):
		for j in range(ny):
			for columnIndex in range(4):
				k = 0
				while k < nz:
					cellIndex = i*nz*ny + j*nz + k
					if inactiveCells[cellIndex] == 0:
						k += 1
						continue

					firstCell = matrixCells[i][j][k]

					topFirstIndex = firstCell[columnIndex]
					firstPoint = arrayOfPoints[topFirstIndex]
					upperPoint = firstPoint
					if k > 0:
						upperCell = matrixCells[i][j][k-1]
						if inactiveCells[cellIndex - 1] != 0 and firstPoint == arrayOfPoints[upperCell[columnIndex+4]]:
							upperPoint = arrayOfPoints[upperCell[columnIndex]]
							
					pinchVertices = [topFirstIndex]
					lowerCell = firstCell
					while(k < nz):
						currCell = matrixCells[i][j][k]
						if inactiveCells[cellIndex] == 0 or not cellHasUpperConnectivity(i, j, k, matrixCells, arrayOfPoints):
							lowerCell = matrixCells[i][j][k-1]
							break
						if not vertexNeedsLaplacian(i, j, k, matrixCells, arrayOfPoints, columnIndex, tol):
							lowerCell = currCell
							break

						pinchVertices.append(currCell[columnIndex+4])
						k += 1
						cellIndex += 1

					nPinch = len(pinchVertices)
					if nPinch <= 1:
						k += 1
						continue
						
					lastCell = matrixCells[i][j][k-1]
					lastPoint = arrayOfPoints[lastCell[columnIndex+4]]
					lowerPoint = arrayOfPoints[lowerCell[columnIndex+4]]

					topLimitPoint = Vector.fromList(firstPoint)*(1 - phi) + Vector.fromList(upperPoint)*(phi)
					botLimitPoint = Vector.fromList(lastPoint) *(1 - phi) + Vector.fromList(lowerPoint)*(phi)
					for pinchIndex in range(nPinch):
						delta = pinchIndex / (nPinch - 1)
						newPoint = topLimitPoint*(1-delta) + botLimitPoint*delta
						arrayOfPoints[pinchVertices[pinchIndex]] = (newPoint.x, newPoint.y, newPoint.z)

					k += 1

	
	return arrayOfPoints

def applyAdaptedSmoothing(pinchVertices, topVertexIndex, botVertexIndex, needsSmoothing, arrayOfPoints, phi):
	if pinchVertices == []:
		return
	nPinch = len(pinchVertices)
	if nPinch <= 1:
		return

	topLimit = Vector.fromList(arrayOfPoints[pinchVertices[0]])*(1 - phi) + Vector.fromList(arrayOfPoints[topVertexIndex])*(phi)
	botLimit = Vector.fromList(arrayOfPoints[pinchVertices[nPinch-1]])*(1 - phi) + Vector.fromList(arrayOfPoints[botVertexIndex])*(phi)
	for index, pinchVertexIndex in enumerate(pinchVertices):
		delta = index / (nPinch - 1)
		newPoint = topLimit*(1-delta) + botLimit*delta
		arrayOfPoints[pinchVertexIndex] = (newPoint.x, newPoint.y, newPoint.z)
		needsSmoothing[pinchVertexIndex] = 0

	return (needsSmoothing, arrayOfPoints)

def enablesFlux1(cell, index, points, inactive, pinchout, porosity, tol):
	if pinchout != [] and pinchout[index] == 0:
		return True
	if inactive != [] and inactive[index] == 0:
		return False
	if hasSmallThickness(cell, points, tol):
		return True
	if porosity != [] and porosity[index] == 0:
		return False
	return True

def enablesFlux(cell, index, points, inactive, pinchout, porosity, tol):
	if pinchout != [] and pinchout[index] == 0:
		return True
	if inactive != [] and inactive[index] == 0:
		return False
	if hasSmallThickness(cell, points, tol):
		return True
	if porosity != [] and porosity[index] == 0:
		return False
	return True
	
def determineFluxCells(matrixCells, numberOfCellValues, arrayOfPoints, qtdCells, inactiveCells, pinchArray, porosity, tol = 0.1):
	(nx, ny, nz) = numberOfCellValues
	newQtdCells = nx*ny*nz
	fluxCells = [1 for i in range(newQtdCells)]
	for i in range(nx):
		for j in range(ny):
			k = 0
			while k < nz:
				cell = matrixCells[i][j][k]
				cellIndex = i*nz*ny + j*nz + k
				if not enablesFlux(cell, cellIndex) or hasSmallThickness(cell, arrayOfPoints, tol):
					fluxCells[cellIndex] = 0
					newQtdCells -= 1
					k += 1
				else:
					intermediateCells = []
					k += 1
					while k < nz:
						currCell = matrixCells[i][j][k]
						currCellIndex = i*nz*ny + j*nz + k
						enablesFlux = enablesFlux(currCell, currCellIndex, arrayOfPoints, inactiveCells, pinchArray, porosity, tol)
						upperConnected = cellHasUpperConnectivity(i, j, k, matrixCells, arrayOfPoints)
						if not enablesFlux or not upperConnected:
							for index in intermediateCells:
								newQtdCells -= 1
								fluxCells[index] = 0
							k += 1
							break
						elif pinchArray[currCellIndex] != 0 and not hasSmallThickness(currCell, arrayOfPoints, tol):
							break
						intermediateCells.append(currCellIndex)
						k += 1

					
	return (fluxCells, newQtdCells)
	

def hexLaplacianBig(matrixCells, arrayOfPoints, numberOfCellValues, qtdCells, inactiveCells, pinchArray, phi, pinchTol = 0.1):
	(nx, ny, nz) = numberOfCellValues
	needsSmoothing = [0  for i in range(len(arrayOfPoints))]
	pinchSequences = [ [ [[], [], [], []] for y in range(ny)] for x in range(nx) ]
	for i in range(nx):
		for j in range(ny):
			for column in range(4):
				k = 0
				while k < nz:
					#Checagem da necessidade de aplicar o laplaciano:
					firstCell = matrixCells[i][j][k]
					firstCellIndex = i*nz*ny + j*nz + k
					if k == 0 or not needsLaplacian(firstCell, firstCellIndex, column, arrayOfPoints, inactiveCells, pinchArray, pinchTol):
						#print("notneeds -> i:", i, "j:", j, "k:", k)
						k += 1
						continue

					upperCell = matrixCells[i][j][k-1]
					upperCellIndex = firstCellIndex - 1
					if inactiveCells[upperCellIndex] == 0:
						k += 1
						continue

					#Definicao dos pontos superiores a regiao suavizada
					firstPointIndex = firstCell[column]
					#print("First point:", firstPointIndex)
					#print("i:", i, "j:", j, "k:", k)
					firstPoint = arrayOfPoints[firstPointIndex]
					upperPointIndex = firstPointIndex
					if firstPoint == arrayOfPoints[upperCell[column+4]]:
							upperPointIndex = upperCell[column]

					#Descobrir sequencia de vertices que serao suavizadas
					firstTopIndex = firstCell[column]
					firstBotIndex = firstCell[column+4]
					pinchVertices = [upperPointIndex, firstTopIndex, firstBotIndex]
					needsSmoothing[firstTopIndex] = 1
					needsSmoothing[firstBotIndex] = 1
					flux = True
					k += 1
					cellIndex = firstCellIndex + 1
					while(k < nz):
						currCell = matrixCells[i][j][k]
						#Se a celula nao obedece o criterio de fluxo, cancelamos a suavizacao
						if inactiveCells[cellIndex] == 0 or \
						   not checkUpperConnectivity(i, j, k, column, matrixCells, arrayOfPoints):
							flux = False
							break
						#Se a celula obedece mas a aresta nao eh fina o suficiente, terminamos a seq. de vertices colapsados
						if not checkEdgeHeight(currCell, column, arrayOfPoints, pinchTol):
							break
						#print("i:", i, "j:", j, "k:", k)
						#Adicionamos o vertice ao vetor
						pinchIndex = currCell[column+4]
						#print("Next point:", pinchIndex)
						pinchVertices.append(pinchIndex)
						needsSmoothing[pinchIndex] = 1
						k += 1
						#print("k -->", k)
						cellIndex += 1
					
					if k == nz:
						break
					if not flux:
						k += 1
						continue

					#Definicao dos pontos da camada abaixo da regiao suavizada
					lastCell = matrixCells[i][j][k-1]
					lastPointIndex = lastCell[column+4]
					#print("Last point:", lastPointIndex)
					lastPoint = arrayOfPoints[lastPointIndex]
					lowerPointIndex = lastPointIndex
					lowerCell = matrixCells[i][j][k]
					if lastPoint == arrayOfPoints[lowerCell[column]]:
						lowerPointIndex = lowerCell[column+4]

					pinchVertices.append(lowerPointIndex)
					#print("Pinch Vertices =", pinchVertices)
					pinchSequences[i][j][column].append(pinchVertices)

	if phi == 0.0:
		return (arrayOfPoints, inactiveCells, qtdCells)
	
	for i in range(nx):
		for j in range(ny):
			for column in range(4):
				for pinchSequence in pinchSequences[i][j][column]:
					debug = False
					for vertex in pinchSequence:
						if vertex == 8603 or vertex == 8604:
							debug = True

					if debug:
						print("Pinches =", pinchSequence)
					#print("pinchSeq", pinchSequence)
					topVertexIndex = pinchSequence[0]
					pinchVertices = []
					for pinchVertexIndex in pinchSequence[1:len(pinchSequence)-1]:
						#print("pinchSeq", pinchSequence)
						#print("pinchV", pinchVertices)
						#print("topIndex", topVertexIndex)
						#print("currIndex", pinchVertexIndex)
						#print("needs", needsSmoothing)
						#print("points", arrayOfPoints)
						if not needsSmoothing[pinchVertexIndex]:
							if pinchVertices == [] or len(pinchVertices) <= 1:
								topVertexIndex = pinchVertexIndex
								continue
							#print("smoothing: top vertex =", topVertexIndex, "vertices =", pinchVertices, "; botVertex =", pinchVertexIndex)
							(needsSmoothing, arrayOfPoints) = applyAdaptedSmoothing(pinchVertices, topVertexIndex, pinchVertexIndex, needsSmoothing, arrayOfPoints, phi)
							topVertexIndex = pinchVertexIndex
							pinchVertices = []
							continue
						pinchVertices.append(pinchVertexIndex)
					if pinchVertices == [] or len(pinchVertices) <= 1:
						continue
					botVertexIndex = pinchSequence[len(pinchSequence)-1]
					#print("smoothing: top vertex =", topVertexIndex, "; botVertex =", botVertexIndex)
					(needsSmoothing, arrayOfPoints) = applyAdaptedSmoothing(pinchVertices, topVertexIndex, botVertexIndex, needsSmoothing, arrayOfPoints, phi)

			
	
	return (arrayOfPoints, inactiveCells, qtdCells)

def removeInactivePoints(matrixCells, arrayOfPoints, numberOfCellValues, out_inactive):
	(nx, ny, nz) = (int(numberOfCellValues[0]), int(numberOfCellValues[1]), int(numberOfCellValues[2]))
	newMatrix = [ [ [ [] for z in range(nz) ] for y in range(ny)] for x in range(nx) ]
	newArrayOfPoints = []
	newIndices = [ -1 for i in range(len(arrayOfPoints)) ]
	newQtdPontos = 0
	
	for k in range(nz):
			for j in range(ny):
				for i in range(nx):
					if out_inactive[i*nz*ny + j*nz + k] == 0:
						continue
					cell = matrixCells[i][j][k]
					newCell = []
					for pointCount in range(8):
						pointIndex = cell[pointCount]
						newIndex = newIndices[pointIndex]
						if newIndex == -1:
							newArrayOfPoints.append(arrayOfPoints[pointIndex])
							newIndex = newQtdPontos
							newIndices[pointIndex] = newIndex
							newQtdPontos += 1
						newCell.append(newIndex)
					
					newMatrix[i][j][k] = newCell
	print("Qtd nova de pontos: ", newQtdPontos)
	return (newMatrix, newArrayOfPoints)

def removeIsolatedPoints(meshPath):
	mesh = pymesh.load_mesh(meshPath)
	mesh, _ = pymesh.remove_isolated_vertices(mesh)
	pymesh.meshio.save_mesh(meshPath, mesh)
	return


def transformTo3dList(coordValues, zcornValues, numberOfCellValues, inactiveCells = [], pinchArray = [], porosity = []):
	#Create bidimensional list with line coordinates
		lineCoords = [] # The number of inner brackets is equal to the number of line coordinates in the Y-direction
		lineCoordsTuples = zip(*[iter(coordValues)] * 6) # Each tuple contains 6 numbers (x1, y1, z1, x2, y2, z2). (x1, y1, z1) represents the first point of a coord line; the other 3 number, the second point

		for i in range(int(numberOfCellValues[1]) + 1):
			xDirLineCoords = [next(lineCoordsTuples) for x in range(int(numberOfCellValues[0]) + 1)]
			lineCoords.append(xDirLineCoords)

		(nx, ny, nz) = (int(numberOfCellValues[0]), int(numberOfCellValues[1]), int(numberOfCellValues[2]))
		matrixCells = [ [ [] for y in range(ny)] for x in range(nx) ]
		arrayOfPoints = []
		#Dicionario sera da forma:
		#  - Chave: Ponto (coordenadas)
		#  - Valor: Lista de indices de pontos nessa coordenada (em ordem crescente de profundidade na malha)
		dicArrayOfPoints = {}
		qtdCelulas = nx * ny * nz
		qtdPontos = 0
		out_inactive = [1 for i in range(qtdCelulas)]
		out_pinchArray = [1 for i in range(qtdCelulas)]
		out_porosity = [-1 for i in range(qtdCelulas)]
		count = 0

		print(inactiveCells)
		
		for k in range(nz):
			for j in range(ny):
				for i in range(nx):
					
					if(pinchArray != [] and not int(pinchArray[i + j*nx + k*nx*ny])):
						out_pinchArray[i*nz*ny + j*nz + k] = 0
					elif (inactiveCells != [] and not int(inactiveCells[i + j*nx + k*nx*ny])):
						qtdCelulas -= 1
						out_inactive[i*nz*ny + j*nz + k] = 0
					
					if porosity != []:
						out_porosity[i*nz*ny + j*nz + k] = porosity[i + j*nx + k*nx*ny]

					beginningPoint = 8*nx*ny*k + 4*nx*j + 2*i

					topNorthWestPoint = getPoint(lineCoords[j][i], zcornValues[beginningPoint] )
					topNorthEastPoint = getPoint(lineCoords[j][i+1], zcornValues[beginningPoint + 1])
					topSouthWestPoint = getPoint(lineCoords[j+1][i], zcornValues[beginningPoint + 2*nx])
					topSouthEastPoint = getPoint(lineCoords[j+1][i+1], zcornValues[beginningPoint + 2*nx + 1])

					botNorthWestPoint = getPoint(lineCoords[j][i], zcornValues[beginningPoint + 4*nx*ny])
					botNorthEastPoint = getPoint(lineCoords[j][i+1], zcornValues[beginningPoint + 4*nx*ny + 1])
					botSouthWestPoint = getPoint(lineCoords[j+1][i], zcornValues[beginningPoint + 4*nx*ny + 2*nx])
					botSouthEastPoint = getPoint(lineCoords[j+1][i+1], zcornValues[beginningPoint + 4*nx*ny + 2*nx +1])

					newCell = (topNorthWestPoint, topNorthEastPoint, topSouthWestPoint, topSouthEastPoint, botNorthWestPoint, botNorthEastPoint, botSouthWestPoint, botSouthEastPoint)
					newCellPointIndex = []

					for counter, point in enumerate(newCell):
						if point in dicArrayOfPoints:
							#Se ponto esta no topo da celula e ja esta no dicionario:
							if counter < 4:
								#Se ele esta na camada do topo (k = 0), entao ele tera o 1o indice do ponto no dicionario
								if k == 0:
									newCellPointIndex.append(dicArrayOfPoints[point][0])
									continue

								supCell = matrixCells[i][j][k-1]
								supPointIndex = supCell[counter + 4]
								#Se ele for igual ao ponto superior a ele, ele tera seu mesmo indice
								if point == arrayOfPoints[supPointIndex]:
									newCellPointIndex.append(supPointIndex)
									continue
								#Senao, ele tera o indice do ponto no dicionario
								else:
									newCellPointIndex.append(dicArrayOfPoints[point][0])
									continue
							
							#Se o ponto esta na base da celula e ja esta no dicionario:
							#  - Senao, ele tera o indice do primeiro na lista do dicionario
							else:
								supPointIndex = newCellPointIndex[counter - 4]
								#Se ele for igual ao ponto superior a ele na celula, ele tera o indice sucessor
								#na lista do dicionario. Se o sucessor nao existir, o geramos
								if point == arrayOfPoints[supPointIndex]:
									idx = dicArrayOfPoints[point].index(supPointIndex) + 1
									if idx < len(dicArrayOfPoints[point]):
										newCellPointIndex.append(dicArrayOfPoints[point][idx])
										continue
									else:
										dicArrayOfPoints[point].append(qtdPontos)
								else:
									newCellPointIndex.append(dicArrayOfPoints[point][0])
									continue
						else:
							#Se o ponto nao esta no dicionario, adicionamos nele o elemento (ponto, [indice])
							dicArrayOfPoints[point] = [qtdPontos]

						#Ponto novo gerado e adicionado a array de pontos
						arrayOfPoints.append(point)
						newCellPointIndex.append(qtdPontos)
						qtdPontos += 1
					#tmp = 1
					#for index in newCellPointIndex:
					#	print("Ponto ", tmp, " = ", index)
					#	tmp += 1
					
					matrixCells[i][j].append(newCellPointIndex)
					count += 1

		return matrixCells, arrayOfPoints, qtdCelulas, out_inactive, out_pinchArray, out_porosity


#pressão, temperatura, altura, gps
def transformTo3dList2(coordValues, zcornValues, numberOfCellValues, inactiveCells):
	#Create bidimensional list with line coordinates
		lineCoords = [] # The number of inner brackets is equal to the number of line coordinates in the Y-direction
		lineCoordsTuples = zip(*[iter(coordValues)] * 6) # Each tuple contains 6 numbers (x1, y1, z1, x2, y2, z2). (x1, y1, z1) represents the first point of a coord line; the other 3 number, the second point

		for i in range(int(numberOfCellValues[1]) + 1):
			xDirLineCoords = [next(lineCoordsTuples) for x in range(int(numberOfCellValues[0]) + 1)]
			lineCoords.append(xDirLineCoords)
		
		#return lineCoords #This return statement is here just for testing.

		(nx, ny, nz) = (int(numberOfCellValues[0]), int(numberOfCellValues[1]), int(numberOfCellValues[2]))
		matrixCells = [ [ [] for y in range(ny)] for x in range(nx) ]
		arrayOfPoints = []
		dicArrayOfPoints = {}
		qtdCelulas = nx * ny * nz
		qtdPontos = 0

		
		for k in range(nz):
			for j in range(ny):
				for i in range(nx):
					
					
					if not int(inactiveCells[i + j*nx + k*nx*ny]):
						qtdCelulas -= 1
						continue

					beginningPoint = 8*nx*ny*k + 4*nx*j + 2*i
					
					topNorthWestPoint = (lineCoords[j][i][0], lineCoords[j][i][1], zcornValues[beginningPoint] )
					topNorthEastPoint = (lineCoords[j][i+1][0], lineCoords[j][i+1][1], zcornValues[beginningPoint + 1])
					topSouthWestPoint = (lineCoords[j+1][i][0], lineCoords[j+1][i][1], zcornValues[beginningPoint + 2*nx])
					topSouthEastPoint = (lineCoords[j+1][i+1][0], lineCoords[j+1][i+1][1], zcornValues[beginningPoint + 2*nx + 1])

					botNorthWestPoint = (lineCoords[j][i][0], lineCoords[j][i][1], zcornValues[beginningPoint + 4*nx*ny])
					botNorthEastPoint = (lineCoords[j][i+1][0], lineCoords[j][i+1][1], zcornValues[beginningPoint + 4*nx*ny + 1])
					botSouthWestPoint = (lineCoords[j+1][i][0], lineCoords[j+1][i][1], zcornValues[beginningPoint + 4*nx*ny + 2*nx])
					botSouthEastPoint = (lineCoords[j+1][i+1][0], lineCoords[j+1][i+1][1], zcornValues[beginningPoint + 4*nx*ny + 2*nx +1])
					


					newCell = (topNorthWestPoint, topNorthEastPoint, topSouthWestPoint, topSouthEastPoint, botNorthWestPoint, botNorthEastPoint, botSouthWestPoint, botSouthEastPoint)
					
					#tempArrayPoints = [point for point in newCell if point not in arrayOfPoints]
					#arrayOfPoints.extend(tempArrayPoints)
					#del(tempArrayPoints)
					#print("Addin newCell...")
					for point in newCell:
						if point not in dicArrayOfPoints:
							arrayOfPoints.append(point)
							dicArrayOfPoints[point] = qtdPontos
							qtdPontos += 1
					#print("New Cell Added!")


					#newCellPointIndex = [arrayOfPoints.index(i) for i in newCell]
					newCellPointIndex = [dicArrayOfPoints[i] for i in newCell]

					matrixCells[i][j].append(newCellPointIndex)
					#matrixCells[i][j].append(newCell)

		return matrixCells, arrayOfPoints, qtdCelulas


def dtopTransformTo3dList(width, numberOfCellValues, dtop):
	(nx, ny, nz) = (numberOfCellValues[0], numberOfCellValues[1], numberOfCellValues[2])
	startingPoint = (0.0, 0.0, 0.0)
	matrixCells = [ [ [] for y in range(ny)] for x in range(nx) ]
	arrayOfPoints = []
	dicArrayOfPoints = {}
	numberOfPoints = 0

	for j in range(ny):
		for i in range(nx):
			currentDTOP = j*nx + i
			for k in range(nz):
				topNorthWest = (startingPoint[0] + width[0]*i, startingPoint[1] + width[1] * j, startingPoint[2] + width[2] * k + dtop[currentDTOP])
				topNorthEast = (startingPoint[0] + width[0]*(i+1), startingPoint[1] + width[1] * j, startingPoint[2] + width[2] * k + dtop[currentDTOP])
				topSouthWest = (startingPoint[0] + width[0]*i, startingPoint[1] + width[1] * (j+1), startingPoint[2] + width[2] * k + dtop[currentDTOP])
				topSouthEast = (startingPoint[0] + width[0]*(i+1), startingPoint[1] + width[1] * (j+1), startingPoint[2] + width[2] * k + dtop[currentDTOP])

				botNorthWest = (startingPoint[0] + width[0]*i, startingPoint[1] + width[1] * j, startingPoint[2] + width[2] * (k+1) + dtop[currentDTOP])
				botNorthEast = (startingPoint[0] + width[0]*(i+1), startingPoint[1] + width[1] * j, startingPoint[2] + width[2] * (k+1) + dtop[currentDTOP])
				botSouthWest = (startingPoint[0] + width[0]*i, startingPoint[1] + width[1] * (j+1), startingPoint[2] + width[2] * (k+1) + dtop[currentDTOP])
				botSouthEast = (startingPoint[0] + width[0]*(i+1), startingPoint[1] + width[1] * (j+1), startingPoint[2] + width[2] * (k+1) + dtop[currentDTOP])

				newCell = (topNorthWest, topNorthEast, topSouthWest, topSouthEast, botNorthWest, botNorthEast, botSouthWest, botSouthEast)

				newCellArray = []
				for point in newCell:
					if point not in dicArrayOfPoints:
						dicArrayOfPoints[point] = numberOfPoints
						numberOfPoints += 1
						arrayOfPoints.append(point)
					newCellArray.append(dicArrayOfPoints[point])
				matrixCells[i][j].append(newCellArray)
	return matrixCells, arrayOfPoints

#@profile
def createNewCells(cell, extPoint, arrayOfPoints, pointsDict, orientation="LEFT"):
	if orientation == "LEFT":
		northWesternPoint = (Vector(*arrayOfPoints[cell[0]]) - Vector(*arrayOfPoints[cell[1]])) + Vector(*extPoint)
		northEasternPoint = Vector(*extPoint) 


	if orientation == "RIGHT":
		northEasternPoint = (Vector(*arrayOfPoints[cell[1]]) - Vector(*arrayOfPoints[cell[0]])) + Vector(*extPoint)
		northWesternPoint = Vector(*extPoint)

	
	southWesternPoint = (Vector(*arrayOfPoints[cell[2]]) - Vector(*arrayOfPoints[cell[0]])) + Vector(northWesternPoint.x, northWesternPoint.y, northWesternPoint.z) #Adding the point found to the vector that goes from topNorthWest to topSouthWest
	southEasternPoint = (Vector(*arrayOfPoints[cell[3]]) - Vector(*arrayOfPoints[cell[1]])) + Vector(northEasternPoint.x, northEasternPoint.y, northEasternPoint.z) #Adding the point found to the vector that goes from topNorthEast to topSouthEast

	topCell = (arrayOfPoints[cell[0]], arrayOfPoints[cell[1]], arrayOfPoints[cell[2]], arrayOfPoints[cell[3]], (northWesternPoint.x, northWesternPoint.y, northWesternPoint.z), (northEasternPoint.x, northEasternPoint.y, northEasternPoint.z), 

				(southWesternPoint.x, southWesternPoint.y, southWesternPoint.z), (southEasternPoint.x, southEasternPoint.y, southEasternPoint.z))

	botCell = ((northWesternPoint.x, northWesternPoint.y, northWesternPoint.z), (northEasternPoint.x, northEasternPoint.y, northEasternPoint.z), 

				(southWesternPoint.x, southWesternPoint.y, southWesternPoint.z), (southEasternPoint.x, southEasternPoint.y, southEasternPoint.z), arrayOfPoints[cell[4]], arrayOfPoints[cell[5]], arrayOfPoints[cell[6]], arrayOfPoints[cell[7]])

	'''
	for point in topCell:
		if point not in arrayOfPoints:
			arrayOfPoints.append(point)

	topCellIndexes = [arrayOfPoints.index(point) for point in topCell]

	for point in botCell:
		if point not in arrayOfPoints:
			arrayOfPoints.append(point)

	
	botCellIndexes = [arrayOfPoints.index(point) for point in botCell]

	'''
	for point in topCell:
		if point not in pointsDict:
			arrayOfPoints.append(point)
			pointsDict[point] = len(arrayOfPoints) - 1

	topCellIndexes = [pointsDict[point] for point in topCell]

	for point in botCell:
		if point not in pointsDict:
			arrayOfPoints.append(point)
			pointsDict[point] = len(arrayOfPoints) - 1

	botCellIndexes = [pointsDict[point] for point in botCell]


	return (topCellIndexes, botCellIndexes)

def triggerRight(nova, rightColumn):
	def eq(a, b):
		return b[0][2] <= a[2] <= b[1][2]
	def lt(a, b):
		return a[2] < b[0][2]
	def gt(a, b):
		return a[2] > b[0][2]
	
	compFuncs = (eq, lt, gt)

	arrayOfIntervalsRightColumn = [(rightCell[0], rightCell[4]) for rightCell in rightColumn[i-1][j]]
	novaTopWestPoint, novaTopEastPoint = nova[0], nova[1]
	index = binarySearch(arrayOfIntervalsRightColumn, currentTopEast, compFuncs)


	if index != -1 and  (arrayOfIntervalsLeftColumn[index][0] < currentTopEastPoint < arrayOfIntervalsLeftColumn[index][1]):
		rightCell = rightColumn[index]
		newCells = createNewCells(rightCell, (novaTopWestPoint, novaTopEastPoint), orientation="RIGHT")
		rightColumn.remove(rightCell)
		rightColumn.insert(index, newCells[0])
		rightColumn.insert(index + 1, newCells[1])

		return newCells
 
	return ()
#@profile
def extend(nova, column, arrayOfPoints, pointsDict, edge="TOP", orientation="LEFT"): #Orientation say in which way we're going to extend our edge; edge, if we are going to extend the top or bot edge. Column is the column in which we are doing the splits.
	#nova is the cell from which we're extending the edges.
	def eq(a, b):
		return b[0][2] <= a[2] <= b[1][2]
	def lt(a, b):
		return a[2] < b[0][2]
	def gt(a, b):
		return a[2] > b[1][2]

	compFuncs = (eq, lt, gt)


	if(orientation == "LEFT"):
		arrayOfIntervalsColumn = [(arrayOfPoints[cell[1]], arrayOfPoints[cell[5]]) for cell in column]
		if(edge=="TOP"):
			novaSearchPoint, novaAdjacentPoint = arrayOfPoints[nova[0]], arrayOfPoints[nova[1]]
		elif(edge=="BOT"):
			novaSearchPoint, novaAdjacentPoint = arrayOfPoints[nova[4]], arrayOfPoints[nova[5]]
	elif(orientation == "RIGHT"):
		arrayOfIntervalsColumn = [(arrayOfPoints[cell[0]], arrayOfPoints[cell[4]]) for cell in column]
		if(edge=="TOP"):
			novaSearchPoint, novaAdjacentPoint = arrayOfPoints[nova[1]], arrayOfPoints[nova[0]]
		elif(edge=="BOT"):
			novaSearchPoint, novaAdjacentPoint = arrayOfPoints[nova[5]], arrayOfPoints[nova[4]]

	
	index = binarySearch(arrayOfIntervalsColumn, novaSearchPoint, compFuncs)

	if index != -1 and  (arrayOfIntervalsColumn[index][0] < novaSearchPoint < arrayOfIntervalsColumn[index][1]):
		adjacentCell = column[index]
		newCells = createNewCells(adjacentCell, novaSearchPoint, arrayOfPoints, pointsDict, orientation)
		column.remove(adjacentCell)
		column.insert(index, newCells[0])
		#column[index] = newCells[0]
		column.insert(index + 1, newCells[1])

		return newCells

	return ()




def destroyHN(matrixCells, arrayOfPoints, numberOfCellValues):
	
	'''
	def eq(a, b):
	 	return b[0][2] <= a[2] <= b[1][2]
	def lt(a, b):
	 	return a[2] < b[0][2]
	def gt(a, b):
	 	return a[2] > b[0][2]
	 compFuncs = (eq, lt, gt)
	'''

	(nx, ny) = (int(numberOfCellValues[0]), int(numberOfCellValues[1])) 

	qtdCelulas = int(numberOfCellValues[0] * numberOfCellValues[1] * numberOfCellValues[2])

	pointsDict = {}
	for numberPoints in range(len(arrayOfPoints)):
		point = arrayOfPoints[numberPoints]
		pointsDict[point] = numberPoints

	def trigger(listaCelulasAtuais, currentX, j):
		listaCelulasFilhas = []
		addNumberCells = 0
		while len(listaCelulasAtuais) > 0 and currentX < nx:
			for celulasAtuais in listaCelulasAtuais:
				for celulaAtual in celulasAtuais:
					celulasFilhasTop = extend(celulaAtual, matrixCells[currentX][j], arrayOfPoints, pointsDict, edge="TOP", orientation="RIGHT") #currentX coluna na qual a divisão está sendo feita
					celulasFilhasBot = extend(celulaAtual, matrixCells[currentX][j], arrayOfPoints, pointsDict, edge="BOT", orientation="RIGHT")
					if len(celulasFilhasTop) > 0:
						addNumberCells += 1
						listaCelulasFilhas.append(celulasFilhasTop)
					if len(celulasFilhasBot) > 0:
						addNumberCells += 1
						listaCelulasFilhas.append(celulasFilhasBot)
			listaCelulasAtuais = listaCelulasFilhas
			listaCelulasFilhas = []
			currentX += 1
		return addNumberCells

	 	#Possibilidade: extend retornar índice em que células criadas estão na coluna. Depois de obter todos os índices, pegar células na coluna com essa informação e colocá-las num SET. 
	 	


	inicio = time()	
	for j in range(ny):
		for i in range(nx-1, -1, -1):
			currentColumn = matrixCells[i][j].copy()
			for cell in currentColumn:
				if i > 0:

					'''
					arrayOfIntervalsLeftColumn = [(leftCell[1], leftCell[5]) for leftCell in matrixCells[i-1][j]]
					currentTopWestPoint, currentTopEastPoint = cell[0], cell[1]

					index = binarySearch(arrayOfIntervalsLeftColumn, currentTopWestPoint, compFuncs)

					if index != -1 and  (arrayOfIntervalsLeftColumn[index][0] < currentTopWestPoint < arrayOfIntervalsLeftColumn[index][1]):
						lefCell = matrixCells[i-1][j][index]
						newCells = createNewCells(leftCell, (currentTopWestPoint, currentTopEastPoint))
						matrixCells[i-1][j].remove(leftCell)
						matrixCells[i-1][j].insert(index, newCells[0]) #Inserting top cell in the position his parent was in
						matrixCells[i-1][j].insert(index + 1, newCells[1]) #Inserting bot cell one position below top cell 
				
						listaCelulasAtuais = [newCells]
						listaCelulasFilhas = []
						currentX = i
					'''

					newCellsTop = extend(cell, matrixCells[i-1][j], arrayOfPoints, pointsDict, edge="TOP", orientation="LEFT")
					listaCelulas = []
					if len(newCellsTop) > 0:
						qtdCelulas += 1
						#print("A quantidade atual é: " + str(qtdCelulas))
						listaCelulas.append(newCellsTop)
						qtdCelulas += trigger(listaCelulas, i, j)

					newCellsBot = extend(cell, matrixCells[i-1][j], arrayOfPoints, pointsDict, edge="BOT", orientation="LEFT")
					listaCelulas = []
					if len(newCellsBot) > 0:
						qtdCelulas += 1
						#print("A quantidade atual é: " + str(qtdCelulas))
						listaCelulas.append(newCellsBot)
						qtdCelulas += trigger(listaCelulas, i, j)

						
						'''
						while len(listaCelulasAtuais) > 0 and currentX < nx:
							for celulasAtuais in listaCelulasAtuais:
								for celulaAtual in celulasAtuais:
									celulasFilhas = triggerRight(celulaAtual, matrixCells[currentX][j]) #currentX coluna na qual a divisão está sendo feita
									if len(celulasFilhas) > 0:
										listaCelulasFilhas.append(celulasFilhas)
							listaCelulasAtuais = listaCelulasFilhas
							listaCelulasFilhas = []
							currentX += 1
						'''
				
				if i < nx - 1:
					

					listaCelulas = []
					newCellsTop = extend(cell, matrixCells[i+1][j], arrayOfPoints, pointsDict, edge="TOP", orientation="RIGHT")
					#print("Extended top right")
					if len(newCellsTop) > 0:
						qtdCelulas += 1
						#print("A quantidade atual é: " + str(qtdCelulas))
						listaCelulas.append(newCellsTop)
						qtdCelulas += trigger(listaCelulas, i+1, j)


					listaCelulas = []
					newCellsBot = extend(cell, matrixCells[i+1][j], arrayOfPoints, pointsDict, edge="BOT", orientation="RIGHT")
					if len(newCellsBot) > 0:
						qtdCelulas += 1
						#print("A quantidade atual é: " + str(qtdCelulas))
						listaCelulas.append(newCellsBot)
						qtdCelulas += trigger(listaCelulas, i+1, j)

					'''
					listaCelulasFilhas = []
					currentX = i+1
								
					while len(listaCelulasAtuais) > 0 and currentX < nx:
						for celulasAtuais in listaCelulasAtuais:
							for celulaAtual in celulasAtuais:
								celulasFilhas = triggerRight(celulaAtual, matrixCells[currentX][j]) #currentX coluna na qual a divisão está sendo feita
								if len(celulasFilhas) > 0:
									listaCelulasFilhas.append(celulasFilhas)
						listaCelulasAtuais = listaCelulasFilhas
						listaCelulasFilhas = []
						currentX += 1
					'''

				#print(matrixCells)
		print("Progress: {0}".format((j/ny) * 100))
	
	fim = time()	

	print("Time spent: {0} seconds".format(fim-inicio))

	#print("A quantidade de células é: " + str(qtdCelulas))
	return qtdCelulas








	'''
	TO-DO:
	1) Percorrer matriz de células, crescendo em y, diminuindo em x, crescendo em z
	2) Para cada célula, descobrir se há hanging node com célula vizinha, à direita e à esquerda.
	Se houver, criar duas novas células, reinseri-las na lista, e iniciar trigger na direção contrária.
	3) Implementar trigger (talvez seja melhor fazer como função separada) em que cada nova célula criada estende arestas para a direção contrária. 
	'''



def transformToEbFVM(matrixCells, arrayOfPoints, numberOfCellValues):
	with open("newEbFVM.txt", "w") as newMesh:
		newMesh.write("1128683573\n")
		newMesh.write(str(len(arrayOfPoints)) + " 0 0 " + str(numberOfCellValues[0] * numberOfCellValues[1] * numberOfCellValues[2]) + " 0 1 1\n") #numberVertices, numberTetra, numberPrism, numberHexa, numberPiramid
		for point in arrayOfPoints:
			newMesh.write(str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
		for x in matrixCells:
			for y in x:
				for z in y:
					for cell in z:
						newMesh.write(str(cell + 1) + " ")
					newMesh.write("\n")


def transformToEbFVM2(matrixCells, arrayOfPoints, qtdCelulas, nome):
	with open(nome + ".txt", "w") as newMesh:
		newMesh.write("1128683573\n")
		newMesh.write("Version number: 5.6D\n")
		newMesh.write(str(len(arrayOfPoints)) + " 0 0 " + str(qtdCelulas) + " 0 1 1\n") #numberVertices, numberTetra, numberPrism, numberHexa, numberPiramid
		for point in arrayOfPoints:
			newMesh.write(str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
		for x in matrixCells:
			for y in x:
				for z in y:
					for cell in z:
						newMesh.write(str(cell + 1) + " ")
					newMesh.write("\n")


def transformToMSH(matrixCells, arrayOfPoints, numberOfCellValues):
	with open("newMSH.msh", "w") as newFile:
		newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
		newFile.write("$Nodes\n" + str(len(arrayOfPoints)) + "\n")
		for (index, point) in enumerate(arrayOfPoints):
			newFile.write(str(index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
		newFile.write("$EndNodes\n")
		newFile.write("$Elements\n" + str(numberOfCellValues[0] * numberOfCellValues[1] * numberOfCellValues[2]) + "\n")

		idCell = 0

		for xDirCells in matrixCells:
			for yDirCells in xDirCells:
				for zDirCells in yDirCells:
					idCell += 1
					newFile.write(
						str(idCell) + " 5 0 " + str(zDirCells[4] + 1) + " " + str(zDirCells[5] + 1) + " " + str(zDirCells[1] + 1) + " " + str(zDirCells[0] + 1) 
						+ " " + str(zDirCells[6] + 1) + " " + str(zDirCells[7] + 1) + " " + str(zDirCells[3] + 1) + " " + str(zDirCells[2] + 1) + "\n"
						)


		newFile.write("$EndElements")

def writeProperty(property, path):
	if property == []:
		return
	with open(path, 'w') as newFile:
		for value in property:
			newFile.write(str(value) + " ")

def outputProperties(permeability, permeabiliy_path, porosity, porosity_path):
	writeProperty(porosity, porosity_path)
	writeProperty(permeability[0], permeabiliy_path + "_x.txt")
	writeProperty(permeability[1], permeabiliy_path + "_y.txt")
	writeProperty(permeability[2], permeabiliy_path + "_z.txt")

def transformToMSH2(matrixCells, arrayOfPoints, qtdCelulas, numberOfCellValues, inactiveCells, pinchArray, permeability, porosity, nome):
	out_permeability = [[], [], []]
	out_porosity = []
	out_pinchArray = []
	with open(nome + ".msh", "w") as newFile:
		newFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
		newFile.write("$Nodes\n" + str(len(arrayOfPoints)) + "\n")
		for (index, point) in enumerate(arrayOfPoints):
			newFile.write(str(index + 1) + " " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
		newFile.write("$EndNodes\n")
		newFile.write("$Elements\n" + str(int(qtdCelulas)) + "\n")

		(nx, ny, nz) = numberOfCellValues
		idCell = 0
		print("qtd =", qtdCelulas)
		print("inac = ", inactiveCells)
		for i in range(nx):
			for j in range(ny):
				for k in range(nz):
					if inactiveCells[i*nz*ny + j*nz + k] == 0:
						continue
					zDirCells = matrixCells[i][j][k]
					idCell += 1


					if pinchArray != []:
						out_pinchArray.append(pinchArray[i + j*nx + k*nx*ny])
					for l in range(3):
						if(permeability[l] != []):
							out_permeability[l].append(permeability[l][i + j*nx + k*nx*ny])
					if porosity != [] and porosity != "":
						out_porosity.append(porosity[i + j*nx + k*nx*ny])

					newFile.write(
						str(idCell) + " 5 0 " + str(zDirCells[4] + 1) + " " + str(zDirCells[5] + 1) + " " + str(zDirCells[1] + 1) + " " + str(zDirCells[0] + 1) 
						+ " " + str(zDirCells[6] + 1) + " " + str(zDirCells[7] + 1) + " " + str(zDirCells[3] + 1) + " " + str(zDirCells[2] + 1) + "\n"
						)

		newFile.write("$EndElements")
	return (out_permeability, out_porosity)


