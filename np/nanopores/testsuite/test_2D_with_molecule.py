from nanopores import *
from checksolve import check_solve

geo_name = "H_geo"
x0 = [0.0, 0.0, 0.0e-9]

# 2D with molecule
generate_mesh(2.0, geo_name, x0=x0)
geo = geo_from_name(geo_name, x0=x0)
p = PNPSAxisym(geo)
check_solve(p)
