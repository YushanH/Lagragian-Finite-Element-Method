import grid_mesh
import elastic

"""Example 3*3 grid, mesh"""
#   *-*-*
#   6 - 7 - 8
#     /   /
#   3 - 4 - 5
#     /   /
#   0 - 1 - 2

N = 3
d = 2
dx = 1/(N-1)
npt = N*N
config = grid_mesh.Config(N,d,dx,npt)
grid = grid_mesh.create_grid(config)
mesh = grid_mesh.create_mesh(config)
#print(grid
# print(grid[0,:])
# print(mesh.size)
# print(grid_mesh.incident_element(config, grid, mesh))

efem = elastic.efem(config, grid, mesh)
# print(efem.Dm)
# print(efem.vol)