"""save/load/modify PNPS force field."""

import numpy, dolfin
import nanopores
from nanopores.models import Howorka

nanopores.add_params(Howorka,
    z0 = None,
    bV = 0.,
    dnaqsdamp = 0.5,
    save = False,
)
# save and load implicit force field
FNAME = "howorka2D_implicit"

def save_force_field(**params):
    F, Fel, Fdrag = Howorka.F_field_implicit(**params)
    mesh = Howorka.geo.mesh
    nanopores.save_functions(FNAME, mesh, meta=params, F=F, Fel=Fel, Fdrag=Fdrag)
    
if save:
    save_force_field(**PARAMS)
    exit()
    
def load_force_field():
    forces, mesh, params = nanopores.load_vector_functions(FNAME)
    #mesh = forces["F"].function_space().mesh()
    return forces["F"], forces["Fel"], forces["Fdrag"], mesh, params

# get forces and accompanying function spaces
F, Fel, Fdrag, mesh, params = load_force_field()
geo, phys = Howorka.setup2D(mesh=mesh, **params)
F0, Fel0, Fdrag0 = F, Fel, Fdrag
mesh0 = geo.mesh
V0 = dolfin.FunctionSpace(mesh, "CG", 1)
VV0 = dolfin.VectorFunctionSpace(mesh, "CG", 1)

# interpolate on fluid submesh
submesh = geo.submesh("fluid")
mesh = submesh
geo, phys = Howorka.setup2D(mesh=mesh, **params)
V = dolfin.FunctionSpace(mesh, "CG", 1)
VV = dolfin.VectorFunctionSpace(mesh, "CG", 1)
F = dolfin.interpolate(F, VV)
Fel = dolfin.interpolate(Fel, VV)
Fdrag = dolfin.interpolate(Fdrag, VV)

v2d = dolfin.vertex_to_dof_map(V)
coord = mesh.coordinates() # numpy array of 2-1 arrays

Rx = params["Rx"]
Ry = params["Ry"]

def function_from_values(values):
    u = dolfin.Function(V)
    u.vector()[v2d] = numpy.array(values)
    return u    
def function_from_lambda(f):
    values = [f(x) for x in coord]
    return function_from_values(values)

if __name__ == "__main__":
    print "params:", params
    print geo
    dolfin.plot(F)
    dolfin.plot(mesh)
    dolfin.interactive()
