import numpy as np
import math
import sys
sys.path.append('../')
import pymesh

class efem:
    def __init__(self, config, grid, mesh, dirichlet_bc, dirichlet_mapping):
        self.config = config
        self.grid = grid
        self.deformed_grid = grid
        self.mesh = mesh
        self.faces = np.reshape(mesh, (-1,3))
        self.dirchlet_bc = dirichlet_bc
        self.dirichlet_mapping = dirichlet_mapping
        self.initialize()
        self.initialize_nodalmass()
        self.grad_N = self.interpolant_gradient()
        # self.dirchlet_pts
        # self.non_dirichlet_pts
        # self.nodalmass
        # self.M
        # self.num_inside_pts
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
        self.num_inside_pts = len(self.non_dirichlet_pts)

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

        self.M = np.zeros((d*len(self.non_dirichlet_pts),d*len(self.non_dirichlet_pts)))
        for i in range(self.num_inside_pts):
            index = self.non_dirichlet_pts[i]
            for j in range(d):
                self.M[i+j,i+j] = self.nodalmass[index]

    def run(self):
        d = self.config.d
        self.map_dirichlet()
        num_inside_pts = self.num_inside_pts
        phi_0 = np.zeros(d*num_inside_pts)
        for i in range(num_inside_pts):
            phi_0[i*d: (i+1)*d] = self.deformed_grid[i, :]
        phi_1 = phi_0
        phi_2 = phi_0
        for timestep in range(self.config.num_tpt):
            # given phi_(n-1), phi_n, find phi_(n+1) using newton's method
            self.advance_one_step(phi_0, phi_1)
            phi_1, phi_0 = phi_1 + self.dphi(phi_0, phi_1), phi_1

    def advance_one_step(self,phi_0, phi_1):
        dt = self.config.dt
        g = lambda phi: M*phi - 2*M*phi_1 + M*phi_0 - dt*dt*self.external_force(phi)

    def dphi(self, phi_0, phi_1):
        return 0
        pass

    def write_obj(self, filename):
        pymesh.save_mesh(filename, pymesh.form_mesh(self.deformed_grid, self.faces))
    def external_force(self, phi):
        d, num_inside_pts = self.config.d, self.num_inside_pts
        F
        f = np.zeros(d * num_inside_pts)
        for i in range(num_inside_pts):
            for beta in range(d):
                pass
    def interpolant_gradient(self):
        def grad_N(i, alpha, e):
            pass
        return grad_N
