#!/usr/bin/env python

from re import X
import lib.extractData2 as extractData2
import lib.transformData as transformData
import lib.laplacian_smoothing as laplacian
import lib.tetraconversor as tetraconversor
import lib.unstableVolumes as unstableVolumes
import lib.pocos as pocos
import lib.convert_ebfvm as convertEbfvm
import lib.collapsedVolumes as collapsedVolumes
import lib.collapseShortEdges as collapseShortEdges
import os
import argparse
from sys import exit

DEBUG = {
	"hexa": False,
	"volumes": False,
	"tetra": False,
	"wells": False,
	"laplacian": False
}

def is_float(value):
	try:
		float(value)
		return True
	except:
		return False

def check_path(path):
	if((path != "") and not os.path.exists(path)):
		raise FileNotFoundError(f'"{path}" not found')
	return path

def check_permeability_path(path):
	if is_float(path):
		return float(path)

	return check_path(path)

def get_filename_arg(prompt, required=False):
	filename = input(prompt)

	if((required or filename != "") and not os.path.exists(filename)):
		raise FileNotFoundError(f"\"{filename}\" not found")

	return filename

def get_option_arg(prompt, options=[]):
	option = int(input(prompt))

	if(option not in options):
		raise ValueError("Invalid input option")

	return option

coordValues = []
zcornValues = []
numberCellsPerAxis = []
inactiveCells = []
pinchArray = []
permeability = [ [], [], [] ]
porosity = []
have_hanging_nodes = True
tetra_optimization = False
apply_laplacian = False
output_mesh_path = ""
qtdCells = -1

arg_parser = argparse.ArgumentParser(
	prog="converter",
	description='Convert hexa mesh to tetra mesh',
	usage="%(prog)s path type_of_file [options]"
)
arg_parser.add_argument('path', type=str, help='Path to mesh file')
arg_parser.add_argument(
	'file_type',
	type=int,
	choices=[1, 2],
	help='1 for DTOP and 2 for ZCORN'
)
arg_parser.add_argument(
	'--hn-removal-method',
	action='store',
	type=int,
	default=1,
	choices=[0, 1],
	help='0 to not remove HN and 1 to TETRAHEDRALIZATION. Default = 1'
)
arg_parser.add_argument(
	'--wells-path',
	action='store',
	metavar='wells_path',
	type=str,
	default="",
	help='Path to wells configuration file'
)
arg_parser.add_argument(
	'--inactive-cells-path',
	action='store',
	metavar='inactive_cells_path',
	type=str,
	default="",
	help='Path to inactive cells file'
)
arg_parser.add_argument(
	'--pinchout-array-path',
	action='store',
	metavar='pinchout_array_path',
	type=str,
	default="",
	help='Path to pinchout array file'
)
arg_parser.add_argument(
	'--x-permeability',
	action='store',
	metavar='x_permeability',
	type=str,
	default="",
	help='X-permeability value or path to x-permeability array file'
)
arg_parser.add_argument(
	'--y-permeability',
	action='store',
	metavar='y_permeability',
	type=str,
	default="",
	help='Y-permeability value or path to y-permeability array file'
)
arg_parser.add_argument(
	'--z-permeability',
	action='store',
	metavar='z_permeability',
	type=str,
	default="",
	help='Z-permeability value or path to z-permeability array file'
)
arg_parser.add_argument(
	'--porosity',
	action='store',
	metavar='porosity',
	type=str,
	default="",
	help='Porosity value or path to porosity array file'
)
arg_parser.add_argument(
	'--tetra-optimization',
	action='store_true',
	default=False,
	help='Enable optimization for mesh with vertical columns'
)
arg_parser.add_argument(
	'--apply-laplacian-smoothing',
	action='store_true',
	default=False,
	help='Enable laplacian smoothing'
)
arg_parser.add_argument(
	'--enable-volume-critiria',
	action='store_true',
	default=False,
	help='Enable volume critiria to identify collapsed volumes'
)
arg_parser.add_argument(
	'--refinement-threshold',
	type=float,
	default=0.125,
	help='Generates more refined meshs the LOWER the threshold. DEFAULT = 0.125'
)
arg_parser.add_argument(
	'--deactivate-hexa-collapsing',
	action='store_true',
	default=False,
	help='Deactivate the removal of collapsed volumes in the hexa mesh'
)
arg_parser.add_argument(
	'--smoothing-factor',
	type=float,
	default=0,
	help='Generates smoother meshs the HIGHER the factor. DEFAULT = 0.0'
)

args = arg_parser.parse_args()
file_path = check_path(args.path)
corner_point_filename = file_path[:-4]

file_type = args.file_type

hn_removal_method = args.hn_removal_method

if(hn_removal_method == 1):

	# Tetra optimization can only be applied to mesh with vertical columns
	tetra_optimization = args.tetra_optimization

	apply_laplacian = args.apply_laplacian_smoothing

refinement = args.refinement_threshold
smooth_factor = args.smoothing_factor
hexa_volumes_collapsing = not args.deactivate_hexa_collapsing
hexa_volumes_critiria = args.enable_volume_critiria

wells_input_path = check_path(args.wells_path)
active_cells_path = check_path(args.inactive_cells_path)
pinch_array_path = check_path(args.pinchout_array_path)
input_perm = []
input_perm.append(check_permeability_path(args.x_permeability))
input_perm.append(check_permeability_path(args.y_permeability))
input_perm.append(check_permeability_path(args.z_permeability))
input_porosity = check_permeability_path(args.porosity)
print()

print("Extracting data from file...")

with open(file_path) as corner:
	if(file_type == 1):
		(numberCellsPerAxis, dtop, width) = extractData2.extractDataDTOP(corner)
		print(f"Nxyz: $numberCellsPerAxis")
		print(f"DTOP: $dtop")
		print(f"Width: $width")
	elif(file_type == 2):
		(coordValues, zcornValues, numberCellsPerAxis) = extractData2.extractData3(corner)

print("Data extracted!\n")


print("Generating Hexa mesh...")

# Cria matriz tridimensional com todas as celulas da malha corner point
if (file_type == 1):
	(matrixOfCells, arrayOfPoints) = transformData.dtopTransformTo3dList(width, numberCellsPerAxis, dtop)

elif(file_type == 2):

	if(active_cells_path != ""):
		inactiveCells = extractData2.getInactiveCells(active_cells_path)

	if(pinch_array_path != ""):
		pinchArray = extractData2.getInactiveCells(pinch_array_path)
	
	for i in range(3):
		if(input_perm[i] != "" and not is_float(input_perm[i])):
			permeability[i] = extractData2.getPermeability(input_perm[i])

	if(input_porosity != "" and not is_float(input_porosity)):
		porosity = extractData2.getPermeability(input_porosity)

	(matrixOfCells, arrayOfPoints, qtdCells, inactiveCells, pinchArray, porosity, permeability) = transformData.transformTo3dList(coordValues, zcornValues, numberCellsPerAxis, inactiveCells, pinchArray, porosity, permeability)

for i in range(3):
	if(permeability[i] == [] and is_float(input_perm[i])):
		perm = float(input_perm[i])
		permeability[i] = [perm for j in range(qtdCells)]

if is_float(input_porosity):
	poro = float(input_porosity)
	porosity = [poro for i in range(qtdCells)]

(inactiveCells, qtdCells) = transformData.determineFluxCells(matrixOfCells, numberCellsPerAxis, arrayOfPoints, qtdCells, inactiveCells, pinchArray, porosity)

if smooth_factor != 0:
	arrayOfPoints = transformData.hexLaplacian(matrixOfCells, arrayOfPoints, numberCellsPerAxis, inactiveCells, smooth_factor)

(matrixOfCells, arrayOfPoints) = transformData.removeInactivePoints(matrixOfCells, arrayOfPoints, numberCellsPerAxis, inactiveCells)

if(qtdCells < 0):
	qtdCells = numberCellsPerAxis[0] * numberCellsPerAxis[1] * numberCellsPerAxis[2]


if(inactiveCells == []):
	inactiveCells = [1 for i in range(qtdCells)]


#print("Permeabilidades: ", porosity)

# Gera malha de hexa
qtdCelulasDepois = qtdCells
print("Num de cells = ", qtdCells)
(out_permeability, out_porosity) =     transformData.transformToMSH2(matrixOfCells, arrayOfPoints, qtdCelulasDepois, numberCellsPerAxis, inactiveCells, pinchArray, permeability, porosity, corner_point_filename+"MSH")
(out_permAllVol, out_porosityAllVol) = transformData.transformToMSH2(matrixOfCells, arrayOfPoints, qtdCelulasDepois, numberCellsPerAxis, inactiveCells, pinchArray, permeability, porosity, corner_point_filename+"MSH_all_volumes")

#print("Poro:", porosity)


if(pinchArray == []):
	pinchArray = [1 for i in range(qtdCells)]

#print("Permeabilidades Apos: ", porosity)

hexa_mesh_filename = corner_point_filename + "MSH.msh"
output_mesh_path = hexa_mesh_filename

print("Hexa mesh generated\n")

if(DEBUG['hexa']): exit()

print("Removing Collapsed Volumes...")
pocos_inactiveCells = inactiveCells[:]
if(hexa_volumes_collapsing):
	print("NAME: ", hexa_mesh_filename)
	#Esta deletando pinch cells
	inactiveCells = collapsedVolumes.main(hexa_mesh_filename, inactiveCells, pinchArray, hexa_volumes_critiria)

print("Collapsed volumes removed\n")
if(DEBUG['volumes']): exit()

#collapseShortEdges.collapseShortEdges(hexa_mesh_filename)

#transformData.removeIsolatedPoints(hexa_mesh_filename)

# Destroi os hanging nodes presentes na malha
if(hn_removal_method == 1):
	print("\nConverting from Hexa to Tetra...")

	if(have_hanging_nodes):
		tetraconversor.tetraConversor(hexa_mesh_filename, numberCellsPerAxis, inactiveCells, tetra_optimization, refinement)
	else:
		tetraconversor.noHNTetraConversor(hexa_mesh_filename, numberCellsPerAxis, inactiveCells)

	print("Tetra Mesh Generated\n")
	tetra_mesh_filename = f"{corner_point_filename}MSH_tetra.msh"
	output_mesh_path = tetra_mesh_filename

if(DEBUG["tetra"]): exit()

# Generate pocos_input file
pocos_input_args = [
	f"{corner_point_filename}MSH_all_volumes.msh",
	f"{output_mesh_path}",
	int(numberCellsPerAxis[0]),
	int(numberCellsPerAxis[1]),
	int(numberCellsPerAxis[2]),
]

path, _ = os.path.split(corner_point_filename)

permeability_path = os.path.join(path, "out_permeability")
porosity_path = os.path.join(path, "out_porosity.txt")

transformData.outputProperties(out_permeability, permeability_path, out_porosity, porosity_path)

if(wells_input_path == ""):
	wells_input_path = os.path.join(path, "pocos_input.txt")

	print(f"Generating pocos_input.txt")
	pocos.generatePocosInput(wells_input_path, pocos_input_args[2:])
	print(f"Created file: {wells_input_path}\n")

# Generate the pocos_output file

print(f"Generating pocos_output.txt")
locked_vertices = pocos.generatePocosOutput(wells_input_path, pocos_input_args[0], pocos_input_args[1], pocos_inactiveCells)
print(f"Created file: {path}/pocos_output.txt\n")

if(DEBUG["wells"]): exit()

if(apply_laplacian and hn_removal_method == 1):
	print(f"Applying Laplacian Smoothing...")
	laplacian.main(tetra_mesh_filename, locked_vertices)
	print("Laplacian Smoothing - DONE\n")
	output_mesh_path = f"{tetra_mesh_filename[:-4]}_laplacian.msh"

if(DEBUG['laplacian']): exit()

# Convert from MSH to CFX 5

print(f"Converting to CFX5 format")
if(hn_removal_method == 1):
	convertEbfvm.main(output_mesh_path, type_of_mesh="TETRA")
else:
	convertEbfvm.main(output_mesh_path, type_of_mesh="HEXA")

