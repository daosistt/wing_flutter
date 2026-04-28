import numpy as np
import meshio

msh = meshio.read("ellipse.msh")

points_2d = msh.points[:, :2]   # убираем Z=0

cells_tri = msh.cells_dict.get("triangle")
if cells_tri is None:
    raise RuntimeError("No triangle cells found in circle.msh")

try:
    tri_tags = msh.cell_data_dict["gmsh:physical"]["triangle"]
except KeyError:
    tri_tags = np.ones(len(cells_tri), dtype=np.int32)

mesh_tri = meshio.Mesh(
    points=points_2d,
    cells=[("triangle", cells_tri)],
    cell_data={"markers": [tri_tags]},
)
meshio.write("ellipse.xdmf", mesh_tri)
print("Written: ellipse.xdmf + ellipse.h5")
