import numpy as np
import math
import pymesh
import grid_mesh
from scipy.sparse.linalg import cg, LinearOperator
from cons_model import Corotated
from GEO import writeGEO
class efem:
    def __init__(self, config, grid, mesh, dirichlet_bc, dirichlet_mapping, g):
        self.config = config
        self.grid = np.copy(grid)
        self.deformed_grid = np.copy(grid)
        self.mesh = np.copy(mesh)
        self.incident_element = grid_mesh.incident_element(config, grid, mesh)
        self.faces = np.reshape(mesh, (-1,config.d+1))
        self.Ne = mesh.size//(config.d+1)
        self.F = np.zeros((self.Ne, config.d, config.d))
        self.dirchlet_bc = dirichlet_bc
        self.dirichlet_mapping = dirichlet_mapping
        self.g = g
        self.initialize()
        self.initialize_nodalmass()
        self.initialize_gravity()
        self.initalize_interpolant_gradient()

        # self.dirchlet_pts
        # self.non_dirichlet_pts
        # self.nodalmass
        # self.M
        # self.num_inside_pts
    def initialize(self):
        """Initialize Dm, Dminv, Ds, vol, dirichlet_pts, non_dirichlet_pts"""
        N, d, npt = self.config.N, self.config.d, self.config.npt
        Ne = self.Ne
        vol = np.zeros(Ne)
        self.Dm = np.zeros((Ne,d,d))
        self.Dminv = np.zeros((Ne, d, d))
        self.Ds = np.zeros((Ne,d,d))
        self.vol = np.zeros(Ne)
        for e in range(Ne):
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
        self.map_dirichlet()
        self.updateDs_F()

    def build_dirichlet_pts(self):
        dirchlet_pts = []
        non_dirichlet_pts = []
        self.i_to_vectori = []
        vectori = 0
        for i in range(self.config.npt):
            x,y = self.grid[i,:]
            if self.dirchlet_bc(x,y):
                dirchlet_pts.append(i)
                self.i_to_vectori.append(None)
            else:
                non_dirichlet_pts.append(i)
                self.i_to_vectori.append(vectori)
                vectori += 1
        self.dirchlet_pts = dirchlet_pts
        self.non_dirichlet_pts = non_dirichlet_pts
        self.num_inside_pts = len(self.non_dirichlet_pts)

    def map_dirichlet(self):
        for i in self.dirchlet_pts:
            x,y = self.grid[i,:]
            self.deformed_grid[i,:] = self.dirichlet_mapping(x,y)
        # x, y = self.grid[4, :]
        # self.deformed_grid[4, :] = [1,0.5]
        self.updateDs_F()

    def initialize_nodalmass(self):
        """grid point i, nodalmass[i] = mass matrix M_ii, after mass lumping"""
        npt = self.config.npt
        d = self.config.d
        rho = self.config.rho
        self.nodalmass = np.zeros(npt)
        for e in range(self.mesh.size//(d+1)):
            for i in range(d+1):
                mass = self.vol[e]*rho/(d+1)
                self.nodalmass[self.mesh[(d+1)*e+i]] += mass
        self.M = np.zeros((d*self.num_inside_pts,d*self.num_inside_pts))
        for i in range(self.num_inside_pts):
            index = self.non_dirichlet_pts[i]
            for j in range(d):
                self.M[i*d+j,i*d+j] = self.nodalmass[index]

    def initialize_gravity(self):
        d, num_inside_pts = self.config.d, self.num_inside_pts
        f = np.zeros(d * num_inside_pts)
        for vectori,i in enumerate(self.non_dirichlet_pts):
            for beta in range(d):
                f[d*vectori+beta] += self.g[beta]*self.nodalmass[i]
        self.gravity_force = f

    def initalize_interpolant_gradient(self):
        d = self.config.d
        self.grad_N = np.zeros((self.Ne * (d + 1), d))  # (len(mesh), 2)
        self.canonical_grad_N = np.array([[-1, -1], [1, 0], [0, 1]])
        for e in range(self.Ne):
            for alpha in range(d+1):
                self.grad_N[(d+1)*e+alpha] = self.Dminv[e].transpose().dot(self.canonical_grad_N[alpha])

    def run(self, verbose = True):
        d = self.config.d
        num_inside_pts = self.num_inside_pts
        phi = np.zeros(d*num_inside_pts)
        v = np.zeros(d*num_inside_pts)
        for i in range(num_inside_pts):
            phi[i*d: (i+1)*d] = self.deformed_grid[self.non_dirichlet_pts[i], :]
        self.deformed_to_obj("output/frame_0.obj")
        # solve for phi_1
        # print("Initial BE energy =", self.BE_energy(phi_1, phi_0, v_0))
        # print("Initial energy =", self.config.dt ** 2 * self.energy())
        for timestep in range(self.config.num_tpt):
            # given phi_(n-1), phi_n, find phi_(n+1) s.t. g(phi_(n+1)) = 0 using newton's method
            # g = lambda phi: self.M.dot(phi) - 2*self.M.dot(phi_1) + self.M.dot(phi_0) - dt*dt*self.external_force(phi)
            phi, v = self.advance_one_step(phi, v, verbose, timestep)
            # print(f"phi = {phi}, v = {v}")
            self.deformed_to_obj(f"output/frame_{timestep+1}.obj")
            # self.deformed_to_geo(f"frame_{timestep+2}.vtk")


    def advance_one_step(self, phi_prev, v_prev, verbose, timestep, tol_newton = 10e-6, tol_cg = 10e-4):
        dt = self.config.dt
        dim = self.num_inside_pts*self.config.d
        tol_newton *= dt**2
        const = - self.M.dot(phi_prev) - dt*self.M.dot(v_prev)
        # print(self.deformed_grid)
        # print(self.internal_force())
        g = lambda phi: self.M.dot(phi) - dt*dt*self.internal_force() + const
        Dg = lambda phi: self.M - dt*dt*self.Df()   # returns matrix size dim*dim
        phi = np.copy(phi_prev)
        max_iter = 100
        iter = 0
        self.update_phi(phi)
        self.updateDs_F()
        while np.linalg.norm(g(phi)) > tol_newton and iter < max_iter:
            # Df_fast = self.del_f()
            # Dg_fast1 = lambda phi: self.M.dot(phi) - dt * dt * Df_fast(phi)
            # Dg_fast = LinearOperator((dim, dim), matvec=Dg_fast1)
            dphi, res = cg(Dg(phi), -g(phi))
            # dphi, res = cg(Dg_fast, -g(phi))

            phi += dphi
            iter += 1
            # print(phi)
            # print("BE energy =", self.BE_energy(phi, phi_1, phi_0))
            # print("cg residual =", res)
            # print("residual g =", np.linalg.norm(g(phi)), "internal force =", dt * dt * self.internal_force())
            self.update_phi(phi)
            self.updateDs_F()
        if verbose:
            print(f"Newton iter {timestep}, BE energy =", self.BE_energy(phi, phi_prev, v_prev))
            print(f"Newton iter {timestep}, energy =", self.config.dt**2*self.energy())
            print(f"Newton iter {timestep}, residual g =", np.linalg.norm(g(phi)))
        return phi, (phi-phi_prev)/dt

    def internal_force(self):
        """assembly external force vector based on current deformed_grid, Ds, F"""
        """f_int = -grad(e)"""
        d, num_inside_pts = self.config.d, self.num_inside_pts
        f = np.zeros(d * num_inside_pts)
        mu, lambd = self.config.mu, self.config.lambd
        for e in range(self.Ne):
            model = Corotated(mu, lambd, self.F[e, :, :])
            P = model.P()
            for p in range(d + 1):
                i = self.mesh[(d+1)*e+p]
                vectori = self.i_to_vectori[i]
                if vectori != None:
                    for beta in range(d):
                        for gamma in range(d):
                            f[d*vectori+beta] -= self.vol[e]*P[beta,gamma]*self.grad_N[(d+1)*e+p, gamma]
        # Add gravity
        f += self.gravity_force
        return f

    def Df(self):
        """Df = grad(f_int) = -Hessian(e)"""
        d, num_inside_pts = self.config.d, self.num_inside_pts
        df_dphi = np.zeros((d * num_inside_pts, d * num_inside_pts))
        mu, lambd = self.config.mu, self.config.lambd
        for e in range(self.Ne):
            model = Corotated(mu, lambd, self.F[e, :, :])
            dPdF = model.dPdF()
            for p in range(d + 1):
                mesh_i = (d+1)*e+p
                i = self.mesh[mesh_i]
                for q in range(d+1):
                    mesh_j = (d+1)*e+q
                # for mesh_j in self.incident_element[i]:
                    j = self.mesh[mesh_j]
                    vectori, vectorj = self.i_to_vectori[i], self.i_to_vectori[j]
                    if vectori != None and vectorj != None:
                        for beta in range(d):
                            for gamma in range(d):
                                for alpha in range(d):
                                    for epsilon in range(d):
                                        df_dphi[vectori*d+beta, vectorj*d+alpha] -= \
                                            self.vol[e]*dPdF[beta*d+gamma, alpha*d+epsilon]*self.grad_N[mesh_j,epsilon]*self.grad_N[mesh_i,gamma]
        return df_dphi

    def del_f(self):
        """Blackbox function for CG"""
        d = self.config.d
        dP_function_list = []
        for e in range(self.Ne):
            model = Corotated(self.config.mu, self.config.lambd, self.F[e, :, :])
            dP_function_list.append(model.dP)
        def calculate_del_f(del_phi):
            ans = np.zeros(self.num_inside_pts*d)
            del_F_vector = self.del_F(del_phi)
            for e in range(self.Ne):
                dP = dP_function_list[e](del_F_vector[e])
                for p in range(d + 1):
                    i = self.mesh[(d+1)*e+p]
                    vectori = self.i_to_vectori[i]
                    if vectori != None:
                        for beta in range(d):
                            for gamma in range(d):
                                ans[d*vectori+beta] -= self.vol[e] * dP[beta, gamma] * self.grad_N[(d+1)*e+p,gamma]
            return ans
        df = LinearOperator((self.num_inside_pts*d,self.num_inside_pts*d), matvec=calculate_del_f)
        return df

    def del_F(self, del_phi):
        """calculates del_F (difference in F) with del_phi"""
        """For 3*3 left/right dirichlet BC, del_phi = [d1,d4,d7]"""
        d = self.config.d
        result = np.zeros((self.Ne,d,d))
        for e in range(self.Ne):
            dDs = np.zeros((d,d))
            i0 = self.i_to_vectori[self.mesh[(d+1)*e]]
            dx0 = del_phi[d*i0:d*i0+d] if i0 != None else np.zeros(d)
            for j in range(1,d+1):
                i = self.i_to_vectori[self.mesh[(d+1)*e+j]]
                dDs[:,j-1] -= dx0
                if i != None:
                    dDs[:,j-1] += del_phi[d*i:d*i+d]
            result[e,:,:] = dDs.dot(self.Dminv[e])
        return result

    def update_phi(self, phi):
        """updates self.deformed_grid based on phi"""
        for i in range(self.num_inside_pts):
            index = self.non_dirichlet_pts[i]
            for beta in range(self.config.d):
                self.deformed_grid[index, beta] = phi[self.config.d*i+beta]

    def BE_energy(self, phi, phi0, v0):
        dt = self.config.dt
        return 1/2*phi.transpose().dot(self.M).dot(phi)- phi.transpose().dot(self.M).dot(phi0+dt*v0) + self.config.dt**2*self.energy()

    def energy(self):
        ans = 0
        for e in range(self.Ne):
            model = Corotated(self.config.mu, self.config.lambd, self.F[e,:,:])
            ans += self.vol[e]*model.psi()
        return ans

    def updateDs_F(self):
        d = self.config.d
        for e in range(self.Ne):
            if d == 2:
                i0,i1,i2 = self.mesh[(d+1)*e:(d+1)*e+(d+1)]
                self.Ds[e,:,:] = np.stack([self.deformed_grid[i1]-self.deformed_grid[i0], self.deformed_grid[i2]-self.deformed_grid[i0]]).transpose()

        # F = Ds*Dminv
        for e in range(self.Ne):
            self.F[e] = np.dot(self.Ds[e],self.Dminv[e])

    def deformed_to_obj(self, filename):
        deformed = np.append(self.deformed_grid, np.zeros((self.config.npt,1)), 1)

        pymesh.save_mesh(filename, pymesh.form_mesh(deformed, self.faces))

    def deformed_to_geo(self, filename):
        deformed = np.append(self.deformed_grid, np.zeros((self.config.npt,1)), 1)
        writeGEO(deformed, self.mesh, filename)