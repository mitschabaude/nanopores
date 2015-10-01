from nanopores import *
from checksolve import check_solve

name3D = "H_cyl_geo"
x0 = [0.0, 0.0, 0.0]
#PNPProblem.method["iterative"] = False

# 3D with molecule
generate_mesh(15.0, name3D, x0=x0)
geo = geo_from_name(name3D, x0=x0)
p = PNPS(geo)
check_solve(p)

