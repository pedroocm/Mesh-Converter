#!/usr/bin/env python

import sys
import polyscope as ps
import pymesh as pm

if len(sys.argv) == 1:
    print("Please enter the mesh file")
    sys.exit(1)

MESH_NAMES = sys.argv[1:]

MESHES = []
for name in MESH_NAMES:
    MESHES.append(pm.load_mesh(name))

ps.init()

for i, mesh in enumerate(MESHES):
    vertices = mesh.vertices
    faces = mesh.faces
    ps.register_surface_mesh(f"my mesh {i}", vertices, faces)

ps.show()
