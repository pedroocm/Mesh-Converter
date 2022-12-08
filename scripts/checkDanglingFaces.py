import pymesh 

mesh = pymesh.load_mesh('uniao_tetra.msh');

def isFaceDangling(face: list):
	for tetra in mesh.elements:
		intersection = set(face).intersection(set(tetra))
		if ( len(intersection) != 0 ):
			return False
	return True

print("reading faces...")
for i, face in enumerate(mesh.faces):
	if isFaceDangling(face):
		print(f"[PROBLEM] at face {face} index {i}")
	else:
		print("OK")
print("process ended")


