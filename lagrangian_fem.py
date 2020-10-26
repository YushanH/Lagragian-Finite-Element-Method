import numpy as np
import grid_mesh
import elastic
d = 2
N = 5
npt = N*N
dx = 1/(N-1)
T = 10
num_tpt = 100
config = grid_mesh.Config(N,d,dx,npt,T,num_tpt)
grid = grid_mesh.create_grid(config)
mesh = grid_mesh.create_mesh(config)
incident_element = grid_mesh.incident_element(config,grid,mesh)

dirichlet_bc = lambda x,y: x in [0,1] or y in [0,1]    # fix all boundary
dirichlet_bc = lambda x,y: x in [0,1]    # fix left/right boundary
dirichlet_bc = lambda x,y: False           # No dirichlet, Pure Neumann boundary

dirichlet_mapping = lambda x,y: (1.2*x, y)
g = 9.8

efem = elastic.efem(config, grid, mesh, dirichlet_bc, dirichlet_mapping)
efem.initialize() #volume, Dm, Dminv


