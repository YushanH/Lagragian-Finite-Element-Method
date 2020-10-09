import numpy as np
"""Example 3*3 grid, mesh"""
#   *-*-*
#   6 - 7 - 8
#   | / | / |
#   3 - 4 - 5
#   | / | / |
#   0 - 1 - 2
#   043 : triangle 0
#   014 : triangle 1


class Config:
    def __init__(self, N, d, dx, npt, rho = 1):
        self.N = N
        self.d = d
        self.dx = dx
        self.npt = npt
        self.rho = rho

def create_grid(config):
    N, d, dx, npt = config.N, config.d, config.dx, config.npt
    grid = np.zeros((npt, d))
    for i in range(npt):
        if d == 2:
            grid[i, 0] = i % N * dx
            grid[i, 1] = i // N * dx
    return grid

def _flat_index(i,j,N):
    return i+j*N

def create_mesh(config):
    # mesh = [0,4,3, 0,1,4,1,5,4,1,2,5,...]
    # Only need left bottom (N-1)*(N-1) points
    N, d, npt = config.N, config.d, config.npt
    mesh = np.zeros((2*(d+1)*(N-1)*(N-1)), dtype=int)
    for i in range(N-1):
        for j in range(N-1):
            num = j*(N-1)+i
            flat_i = _flat_index(i,j,N)
            mesh[2*(d+1)*num] = flat_i
            mesh[2*(d+1)*num+1] = flat_i + N + 1
            mesh[2*(d+1)*num+2] = flat_i + N
            mesh[2*(d+1)*num+3] = flat_i
            mesh[2*(d+1)*num+4] = flat_i + 1
            mesh[2*(d+1)*num+5] = flat_i + N + 1
    return mesh

def incident_element(config, grid, mesh) -> [set]:
    """inci_ele[i] returns list [j] such that mesh[j] == i
    inci_ele[0] = [0,3], inci_ele[1] = [4,6,9]"""

    N, d, npt = config.N, config.d, config.npt
    inci_ele = []
    for i in range(npt):
        inci_ele.append(set())
    for j in range(mesh.size):
        inci_ele[mesh[j]].add(j)
    return inci_ele
if __name__ == "__main__":
    N = 3
    d = 2
    dx = 1/(N-1)
    npt = N*N
    example_config = Config(N,d,dx,npt)
    grid = create_grid(example_config)
    mesh = create_mesh(example_config)
    #print(grid)
    print(grid[0,:])
    print(mesh)
    print(incident_element(example_config, grid, mesh))
