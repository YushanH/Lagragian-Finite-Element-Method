import numpy as np

def writeOBJ(filename, points, faces):
    # Tri mesh to obj
    with open(filename, "w") as file:
        for p in points:
            file.write(f"v {p[0]} {p[1]} 0\n")
        for f in faces:
            file.write(f"f {f[0]+1} {f[1]+1} {f[2]+1}\n")

if __name__ == "__main__":
    points = np.array([[0.1,0.1],[0.3,0.5]])
    writeOBJ("test_geometry.obj", points)