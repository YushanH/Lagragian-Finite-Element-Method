import grid_mesh
import elastic

"""Example 3*3 grid, mesh"""
#   *-*-*
#   6 - 7 - 8
#   | / | / |
#   3 - 4 - 5
#   | / | / |
#   0 - 1 - 2

N = 3
d = 2
dx = 1/(N-1)
npt = N*N
T = 10
num_tpt = 100
config = grid_mesh.Config(N,d,dx,npt,T,num_tpt)
grid = grid_mesh.create_grid(config)
mesh = grid_mesh.create_mesh(config)
dirichlet_bc = lambda x,y: x in [0,1] or y in [0,1]    # fix all boundary
dirichlet_mapping = lambda x,y: (1.2*x, y)
# print(grid)
# print(grid[0,:])
# print(mesh.size)
# print(grid_mesh.incident_element(config, grid, mesh))

efem = elastic.efem(config, grid, mesh, dirichlet_bc, dirichlet_mapping)


# print(efem.Dm)
# print(efem.Dminv)
# print(efem.vol)
# print(efem.nodalmass)
# print(efem.dirchlet_pts)
# print(efem.non_dirichlet_pts)
efem.run()
#print(efem.deformed_grid)
print(efem.M)