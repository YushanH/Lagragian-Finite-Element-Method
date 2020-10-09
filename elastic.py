import numpy as np
import math
import sys
sys.path.append('../')
# import PyMesh.python.pymesh
class efem:
    def __init__(self, config, grid, mesh, dirichlet_bc, dirichlet_mapping):
        self.config = config
        self.grid = grid
        self.deformed_grid = grid
        self.mesh = mesh
        self.dirchlet_bc = dirichlet_bc
        self.dirichlet_mapping = dirichlet_mapping
        self.initialize()
        self.initialize_nodalmass()
        # self.dirchlet_pts
        # self.non_dirichlet_pts
        # self.nodalmass
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
        self.build_dirichlet_pts()

    def build_dirichlet_pts(self):
        dirchlet_pts = []
        non_dirichlet_pts = []
        for i in range(self.config.npt):
            x,y = self.grid[i,:]
            if self.dirchlet_bc(x,y):
                dirchlet_pts.append(i)
            else:
                non_dirichlet_pts.append(i)
        self.dirchlet_pts = dirchlet_pts
        self.non_dirichlet_pts = non_dirichlet_pts

    def map_dirichlet(self):
        for i in self.dirchlet_pts:
            x,y = self.grid[i,:]
            self.deformed_grid[i,:] = self.dirichlet_mapping(x,y)
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

    def run(self):
        d = self.config.d
        self.map_dirichlet()
        num_inside_pts = len(self.non_dirichlet_pts)
        phi_0 = np.zeros(d*num_inside_pts)
        phi = phi_0
        for i in range(num_inside_pts):
            phi_0[i*d: (i+1)*d] = self.grid[self.dirchlet_pts[i], :]
        for timestep in range(self.config.num_tpt):
            phi += self.dphi(phi)
            pass
    def dphi(self, phi):
        return 0
        pass

    def write_obj(self, filename):
        pass


