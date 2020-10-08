import numpy as np
import grid_mesh
import elastic
d = 2
N = 5
npt = N*N
dx = 1/(N-1)

config = grid_mesh.Config(N,d,dx,npt)
grid = grid_mesh.create_grid(config)
mesh = grid_mesh.create_mesh(config)
incident_element = grid_mesh.incident_element(config,grid,mesh)

dirichlet_bc = lambda x,y: x == 0 or x == 1     # fix left and right side
g = 9.8

efem = elastic.efem(config, grid, mesh)
efem.initialize() #volume, Dm, Dminv


#nodal mass