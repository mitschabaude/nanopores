from nanopores import *
from .checksolve import check_solve

name = "H_geo"

# 2D without molecule
generate_mesh(2.0, name, dim=2)
geo = geo_from_name(name)

# hack to prevent "Warning: found no facets matching domain for bc"
# TODO: find better way to resolve this (dynamic synonymes in .subdomains)
exec("from nanopores.%s.subdomains import synonymes" %name)
geo.import_synonymes({"moleculeb":set()})
geo.import_synonymes(synonymes)

p = PNPSAxisym(geo)
check_solve(p)

