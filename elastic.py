import numpy as np
import math
class efem:
    def __init__(self, config, grid, mesh):
        self.config = config
        self.grid = grid
        self.mesh = mesh
        self.calculate_volume()
    def initialize(self):
        N, d, npt = self.config.N, self.config.d, self.config.npt
        num_triangles = self.mesh.size//(d+1)
        vol = np.zeros(num_triangles)
        self.Dm = np.zeros((num_triangles,d,d))
        self.Dminv = np.zeros((num_triangles, d, d))
        self.vol = np.zeros(num_triangles)
        for i in range(num_triangles):
            i0,i1,i2 = self.mesh[(d+1)*i], self.mesh[(d+1)*i+1],self.mesh[(d+1)*i+2]
            x0 = self.grid[i0,:]
            x1 = self.grid[i1,:]
            x2 = self.grid[i2,:]
            D_m = np.stack((x1-x0,x2-x0)).transpose()
            self.Dm[i,:,:] = D_m
            self.Dminv[i,:,:] = np.linalg.inv(D_m)
            self.vol[i] = (1/math.factorial(d))*abs(np.linalg.det(D_m))


