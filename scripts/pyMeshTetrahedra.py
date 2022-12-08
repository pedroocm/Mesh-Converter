import pymesh
import os
import sys

if len(sys.argv) < 2:
    print("Deu pau")
    sys.exit(1)

mesh = pymesh.load_mesh(sys.argv[1])

new_mesh = pymesh.hex_to_tet(mesh)

pymesh.save_mesh("pymesh_tetrahedra.msh", new_mesh)
