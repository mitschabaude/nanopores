#!/usr/bin/env python
''' create mesh and midpoints file for force calculations '''

from nanopores import *
from dolfin import *
from numpy import *

from protocol import Data
import os, sys

DIR = DATADIR + "/sim"
simname = "sim"

if not os.path.exists(DIR):
    os.makedirs(DIR)
    
if len(sys.argv) > 1:
    simname = sys.argv[1]

FILENAME = DIR + "/" + simname + ".dat"
ndigits = 5

def create_data_file(mesh, lscale):
    midpoints = [CellFunction("double", mesh, 0.0) for i in [0,1]]
    for c in cells(mesh):
        # be careful with precision of x0, therefore round
        # otherwise there are problems with mesh generation and subdomain evaluation
        midpoints[0][c] = round(c.midpoint()[0]*lscale, ndigits)
        midpoints[1][c] = round(c.midpoint()[1]*lscale, ndigits)

    N = mesh.num_cells()

    data = Data(FILENAME, N=N, overwrite=True)
    data["r"] = midpoints[0].array()
    data["z"] = midpoints[1].array()

    return data.filename

name = "H_geo"
geoparams = dict(
x0=None,
lcCenter=0.32e-9,
Ry = 9*nm,
Rx = 4*nm,
)
lscale = 1e9

geo_dict = generate_mesh(.4, name, **geoparams)
geo1 = geo_from_name(name, **geoparams)
#print geo1._physical_boundary
#geo1.import_synonymes({"chargeddnab":{"chargeddnainb","unchargeddnab"}})
#from nanopores.tools.geometry import _invert_dict_nonunique
#geo1._bou2phys = _invert_dict_nonunique(geo1._physical_boundary)
#print geo1._bou2phys
#phys = Physics("pore_molecule", geo1, Membraneqs = 0.)
#print phys.charge
#pde = LinearPBAxisym(geo1, phys)
#pde.maxcells = 5000
#pde.marking_fraction = 0.2
#pde.solve(refinement=True)
mesh = geo1.mesh

geo = geo_from_subdomains(mesh, "inside", **geoparams)
mesh_o = geo.submesh("moleculeinfluid")
geo_i = geo_from_subdomains(mesh_o, "inside-outside", **geoparams)
mesh_i = geo_i.submesh("moleculeinfluid")

plot(mesh_o)
plot(mesh_i)
interactive()

File(DIR + "/" + simname + "_mesh_i.xml") << mesh_i
File(DIR + "/" + simname + "_mesh_o.xml") << mesh_o

print "Created data file in", create_data_file(mesh_i, lscale)
print "Number of cells:", mesh_i.num_cells()
print "Rounded midpoints to",ndigits ,"digits"
