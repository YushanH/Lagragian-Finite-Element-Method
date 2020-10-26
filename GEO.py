import pygmsh
from py2gmsh import Mesh, Entity
def writeGEO(grid, mesh, filename):
    """Example: grid = np.array([[0,0],[0.5,0.5],[0,0.5]]), mesh = np.array([1,4,3])"""
    # with pygmsh.geo.Geometry() as geom:
    #     for e in range(mesh.shape[0] // 3):
    #         p1 = grid[mesh[3*e]]
    #         p2 = grid[mesh[3*e+1]]
    #         p3 = grid[mesh[3*e+2]]
    #
    #         geom.add_polygon(
    #             [
    #                 list(p1),list(p2),list(p3)
    #             ],
    #         )
    #     mesh = geom.generate_mesh()
    #
    # mesh.write(filename)

    mesh_geo = Mesh()
    pt_list = []
    for i in range(grid.shape[0]):
        p = Entity.Point(grid[i,:])
        pt_list.append(p)
        mesh_geo.addEntity(p)
    for e in range(mesh.shape[0]//3):
        p1 = pt_list[mesh[3*e]]
        p2 = pt_list[mesh[3*e+1]]
        p3 = pt_list[mesh[3*e+2]]
        l1 = Entity.Curve([p1, p2])
        l2 = Entity.Curve([p2, p3])
        l3 = Entity.Curve([p3, p1])
        # entities can also be added in a batch
        mesh_geo.addEntities([l1, l2, l3])
    mesh_geo.writeGeo(filename)

