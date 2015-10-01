from nanopores import *
from dolfin import *

name = "H_cyl_geo"

geo = geo_from_name(name)

# ions on dna?
geo.import_synonymes({"ions":{"fluid","dna"}})
IllposedNonlinearSolver.newtondamp = 0.9

PNPS.maxcells = 100000

pnps = PNPS(geo)

pnps.solve(refinement=False)
#pnps.save_mesh()
#pnps.print_results()
pnps.visualize("center")

