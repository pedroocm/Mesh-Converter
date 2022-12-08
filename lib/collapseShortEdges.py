import pymesh


def collapseShortEdges(meshPath):
    mesh = pymesh.load_mesh(meshPath)
    mesh = pymesh.collapse_short_edges(mesh, rel_threshold=0.05)
    pymesh.save_mesh(meshPath, mesh)