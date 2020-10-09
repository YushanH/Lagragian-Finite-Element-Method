import numpy as np
import math
class efem:
    def __init__(self, config, grid, mesh):
        self.config = config
        self.grid = grid
        self.mesh = mesh
        self.initialize()
        self.initialize_nodalmass()
    def initialize(self):
        N, d, npt = self.config.N, self.config.d, self.config.npt
        num_triangles = self.mesh.size//(d+1)
        vol = np.zeros(num_triangles)
        self.Dm = np.zeros((num_triangles,d,d))
        self.Dminv = np.zeros((num_triangles, d, d))
        self.vol = np.zeros(num_triangles)
        for e in range(num_triangles):
            index_list = []
            point_list = []
            vector_list = []
            for j in range(d+1):
                index = self.mesh[(d+1)*e+j]
                index_list.append(index)
                point_list.append(self.grid[index,:])
            for j in range(1,d+1):
                vector_list.append(point_list[j]-point_list[0])
            D_m = np.stack(vector_list).transpose()
            self.Dm[e,:,:] = D_m
            self.Dminv[e,:,:] = np.linalg.inv(D_m)
            self.vol[e] = (1/math.factorial(d))*abs(np.linalg.det(D_m))
    def initialize_nodalmass(self):
        """grid point i, nodalmass[i] = mass matrix M_ii, after mass lumping"""
        npt = self.config.npt
        d = self.config.d
        rho = self.config.rho
        self.nodalmass = np.zeros(npt)
        for e in range(self.mesh.size//(d+1)):
            for i in range(d+1):
                mass = self.vol[e]*rho/3
                self.nodalmass[self.mesh[(d+1)*e+i]] += mass



