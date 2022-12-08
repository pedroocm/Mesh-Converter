#!/usr/bin/env python

import lib.extractData2 as extractData2 
import lib.transformData as transformData
import lib.laplacian_smoothing as laplacian
from sys import exit

coordValues = []
zcornValues = []
numberCellsPerAxis = []

cornerPointFileName = input("Enter the file name of the corner point mesh you wish to convert. It must be a .txt file. Type only the name, without the ext: ")
print("Type 1 to use DTOP values.\nType 2 to use ZCORN values: ")
dtop = int(input("> "))

if(dtop != 1 and dtop != 2):
	print("Valor inválido! Saindo...")
	exit()

print("Hanging nodes removal method:\n Type 1 for TETRAHEDRALIZATION.\n Type 2 for HEXAHEDRA SUBDIVISION")
method = int(input("> "))

print("Type 1 to do not use smoothing.\nType 2 to use Laplacian Smoothing:")
smoothing = int(input("> "))

try:
	with open( cornerPointFileName + ".txt") as corner:
		print("Extracting data from file...")
		if(dtop == 2):
			(coordValues, zcornValues, numberCellsPerAxis) = extractData2.extractData3(corner)
		else:
			(numberCellsPerAxis, dtop, width) = extractData2.extractDataDTOP(corner)
		print("Data extracted!\n")
except IOError as e:
	print("Error reading file: {0}".format(e.strerror))
	exit(1)

#Cria matriz tridimensional com todas as células da malha corner point

print("Converting...")

if (dtop == 2):
	(matrixOfCells, arrayOfPoints) = transformData.transformTo3dList(coordValues, zcornValues, numberCellsPerAxis)
else:
	(matrixOfCells, arrayOfPoints) = transformData.dtopTransformTo3dList(width, numberCellsPerAxis, dtop)

#Destroi os hangind nodes presentes na malha

print("Chegou aqui")
if(method == 2):
	print("Removing hanging nodes...")
	qtdCelulasDepois = transformData.destroyHN(matrixOfCells, arrayOfPoints, numberCellsPerAxis)
	#qtdCelulasDepois = numberCellsPerAxis[0] * numberCellsPerAxis[1] * numberCellsPerAxis[2]
	print("Hanging nodes removed!")

	#Cria mesh para EbFVM:

	transformData.transformToEbFVM2(matrixOfCells, arrayOfPoints, qtdCelulasDepois, cornerPointFileName+"EbFVM")
	print(cornerPointFileName+"EbFVM.txt file created!")

	#Cria arquivo em formato MSH para visualização no GMSH:
	transformData.transformToMSH2(matrixOfCells, arrayOfPoints, qtdCelulasDepois, cornerPointFileName+"MSH")
	print(cornerPointFileName+"MSH.msh file created!")

elif(method == 1):
	print("Entrou")
	import lib.tetraconversor as tetraconversor
	qtdCelulasDepois = numberCellsPerAxis[0] * numberCellsPerAxis[1] * numberCellsPerAxis[2]
	transformData.transformToMSH2(matrixOfCells, arrayOfPoints, qtdCelulasDepois, cornerPointFileName+"MSH")
	tetraconversor.tetraConversor(cornerPointFileName+"MSH.msh", [int(numberCellsPerAxis[0]), int(numberCellsPerAxis[1]), int(numberCellsPerAxis[2])])

mesh_file = f"{cornerPointFileName}MSH.msh"

if(smoothing == 2):
	print(f"Applying Laplacian Smoothing to {mesh_file}")
	laplacian.main(mesh_file)
	print("Laplacian Smoothing - DONE")