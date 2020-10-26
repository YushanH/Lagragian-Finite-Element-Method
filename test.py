import grid_mesh
import elastic
import time
from cons_model import Corotated
import numpy as np
"""Example 3*3 grid, mesh"""
#   *-*-*
#   6 - 7 - 8       6  -  7  -  8
#   | / | / |       |  / /  /  |
#   3 - 4 - 5       3 - 4 ----- 5
#   | / | / |       | /  \   /  |
#   0 - 1 - 2       0  -  1  -  2

N = 5
d = 2
dx = 1/(N-1)
npt = N*N
T = 2
num_tpt = 10
E, nu = 10, 0.3 # Young's modulus/Poisson's ratio
mu, lambd = E/(2*(1+nu)), E*nu/((1+nu)*(1-2*nu))
# mu, lambd = 1,1

rho = 10
# g = np.array([0, 10*-0.98])
g = np.array([0, 0])
config = grid_mesh.Config(N,d,dx,npt,T,num_tpt,rho=rho,mu=mu,lambd=lambd)
grid = grid_mesh.create_grid(config)
mesh = grid_mesh.create_mesh(config)
dirichlet_bc = lambda x,y: x in [0,1] or y in [0,1]    # fix all boundary
# dirichlet_bc = lambda x,y: x in [0,1]    # fix left/right boundary
# dirichlet_bc = lambda x,y: False           # No dirichlet, Pure Neumann boundary
dirichlet_mapping = lambda x,y: (2*x, y)
# dirichlet_mapping = lambda x,y: (x, y)
# print(grid)
# print(grid[0,:])
# print(mesh.size)
# print(grid_mesh.incident_element(config, grid, mesh))

efem = elastic.efem(config, grid, mesh, dirichlet_bc, dirichlet_mapping,g)
# for i in [2,5,8]:
#     efem.deformed_grid[i,0] *= 2
# efem.updateDs_F()
# for i,m in enumer               ate(mesh):
#     print(efem.grad_N[i,:], m)
# # print(efem.F)
# print("START")
# for f in efem.F:
#     print(Corotated(mu,lambd,f).P())
# print("END")
# print(efem.internal_force())

# print(efem.Dminv)
# print(efem.vol)
# print(efem.nodalmass)
# print(efem.dirchlet_pts)
# print(efem.non_dirichlet_pts)

start = time.time()
efem.run(verbose=True)
end = time.time()
print(f"Elapsed time = {end-start} seconds")
# print(efem.deformed_grid)
# print(efem.M)
# print(efem.faces)
# efem.deformed_to_obj("example.obj")
# print(efem.F)
# model = cons_model.Corotated(mu,lambd, efem.F[0])
# model = cons_model.Corotated(mu,lambd, np.diag([1,2]))
# efem.updateDs()
# efem.updateF()
# print(efem.Ds[:3,:,:], efem.Dm[:3,:,:])
# print(efem.F[:3,:,:])
# print(model.dPdF())


# Debugging e ~ f ~ Df
def test_e_f_Df():
    mu, lambd = 1, 1
    g = np.array([0, -0.98])
    g = np.array([0, 0])
    config = grid_mesh.Config(N, d, dx, npt, T, num_tpt)
    grid = grid_mesh.create_grid(config)
    mesh = grid_mesh.create_mesh(config)
    dirichlet_bc = lambda x, y: x in [0, 1] or y in [0, 1]  # fix all boundary
    dirichlet_mapping = lambda x, y: (2 * x, y)
    direction = np.random.rand(2)
    direction_F = np.random.rand(2,2)
    eps = 0.1
    epsilon = [eps]
    iter = 5
    for i in range(iter):
        eps /= 2
        epsilon.append(eps)
    print("epsilon =", epsilon)
    error_ef = []
    error_fDf = []

    for i in range(5):
        efem = elastic.efem(config, grid, mesh, dirichlet_bc, dirichlet_mapping, g)
        f = efem.internal_force()
        del_f = efem.del_f()
        Df = efem.Df()
        efem.deformed_grid[4] += epsilon[i]*direction
        efem.updateDs_F()
        # print(efem.deformed_grid[4])
        # print(efem.internal_force())
        # print(efem.F[0])
        energy_plus = efem.energy()
        force_plus = efem.internal_force()

        efem.deformed_grid[4] -= 2*epsilon[i]*direction
        efem.updateDs_F()
        # print(efem.deformed_grid[4])
        # print(efem.F[0])
        energy_minus = efem.energy()
        force_minus = efem.internal_force()
        error_ef.append(-(energy_plus-energy_minus)/2/epsilon[i]-f.dot(direction))
        error_fDf.append(np.linalg.norm((force_plus-force_minus)/2/epsilon[i]- Df.dot(direction)))
        # print((force_plus-force_minus)/2, del_f(epsilon[i]*direction))
    print("error e vs.f\n", error_ef)
    print("error f vs.Df\n", error_fDf)

# test_e_f_Df()