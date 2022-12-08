import pymesh
import threading

mesh = pymesh.load_mesh("./outputs/output_mesh2/converted_meshs/first_tetrahedralization.msh")
mesh.enable_connectivity()

def isNumberNeighborsCorrect(tetraIndex, tetra_neighborsList):
	tetraVerticesIndex = set(mesh.elements[tetraIndex])
	tetra_neighbors = set(tetra_neighborsList)

	for i in range(mesh.num_elements):
		if (i != tetraIndex) and (i not in tetra_neighbors):
			intersection = tetraVerticesIndex.intersection(set(mesh.elements[i]))
			if len(intersection) == 3:
				return False
	return True

def checkMeshCorrectness(start, end):
	print("Started")
	for i in range(start, end):
		#print(f"Checking tetrahedra {i + 1}...")
		tetra_neighborsList = mesh.get_voxel_adjacent_voxels(i)
		number_of_neighbors = len(tetra_neighborsList)
		if (number_of_neighbors < 4 and (not isNumberNeighborsCorrect(i, tetra_neighborsList) ) ):
			print(f"[PROBLEM] Problem at tetrahedrum {i + 1} with neighbours  {[k + 1 for k in tetra_neighborsList]} ")
		#else:
			#print("No problem")
	print("Finished")


threads = list()
numOfThreads = 33
numElemsPerThread = mesh.num_elements // numOfThreads

print(numElemsPerThread)


for t in range(numOfThreads):
	x = threading.Thread(target=checkMeshCorrectness, args=(t * numElemsPerThread, t * numElemsPerThread + numElemsPerThread))
	threads.append(x)
	x.start()

for thread in threads:
	thread.join()




