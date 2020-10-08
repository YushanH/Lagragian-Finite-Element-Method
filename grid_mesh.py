import numpy as np
"""Example 3*3 grid, mesh"""
#   *-*-*
#   6 - 7 - 8
#     /   /
#   3 - 4 - 5
#     /   /
#   0 - 1 - 2

class Config:
    def __init__(self, N, d, dx, npt):
        self.N = N
        self.d = d
        self.dx = dx
        self.npt = npt


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
    #   3 - 4 - 5
    #     /   /
    #   0 - 1 - 2
    # mesh = [0,1,4,0,4,3,1,2,5,1,5,4,...]
    # Only need left bottom (N-1)*(N-1) points
    N, d, npt = config.N, config.d, config.npt
    mesh = np.zeros((2*(d+1)*(N-1)*(N-1)), dtype=int)
    for i in range(N-1):
        for j in range(N-1):
            num = j*(N-1)+i
            flat_i = _flat_index(i,j,N)
            mesh[2*(d+1)*num] = flat_i
            mesh[2*(d+1)*num+1] = flat_i + 1
            mesh[2*(d+1)*num+2] = flat_i + N + 1
            mesh[2*(d+1)*num+3] = flat_i
            mesh[2*(d+1)*num+4] = flat_i + N + 1
            mesh[2*(d+1)*num+5] = flat_i + N
    return mesh

def incident_element(config, grid, mesh) -> [set]:
    """returns neighbors of point i in mesh"""
    N, d, dx, npt = config.N, config.d, config.dx, config.npt
    inci_ele = []
    for i in range(npt):
        inci_ele.append(set())
        for e in range(mesh.size//((d+1))):
            if i in mesh[(d+1)*e:(d+1)*e+(d+1)]:
                for neighbor in mesh[(d+1)*e:(d+1)*e+(d+1)]:
                    inci_ele[i].add(neighbor)
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
    print(mesh.size)
    print(incident_element(example_config, grid, mesh))
