import dolfin
import nanopores
import Howorka
from nanopores.physics.convectiondiffusion import ConvectionDiffusion

nanopores.add_params(
    save = False,
    **Howorka.PARAMS)
PARAMS.pop("z0")

# save and load implicit force field
FNAME = "howorka2D_implicit"

def save_force_field(**params):
    F, Fel, Fdrag = Howorka.F_field_implicit(**params)
    mesh = Howorka.geo.mesh
    nanopores.save_functions(FNAME, mesh, F=F, Fel=Fel, Fdrag=Fdrag)
    
if save:
    save_force_field(**PARAMS)
    exit()
    
def load_force_field():
    forces, mesh = nanopores.load_vector_functions(FNAME)
    #mesh = forces["F"].function_space().mesh()
    return forces["F"], forces["Fel"], forces["Fdrag"], mesh
        
def convect(geo, phys, Fel, Fdrag):
    #F, Fel, Fdrag, mesh = load_force_field()
    #geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **params)
    #dolfin.plot(F)
    #frac = .01
    #t = 1e-9    
    #dt = t*frac
    t = 5e-7
    dt = 1e-9 
    
    bc = dict(upperb=dolfin.Constant(1.), lowerb=dolfin.Constant(0.))
    u0 = dolfin.Constant(0.)
    #F = dolfin.Constant((0.,0.))
    F = Fel + dolfin.Constant(3.)*Fdrag
    pde = ConvectionDiffusion(geo, phys, dt=dt, F=F, u0=u0, bc=bc, cyl=True)
    #pde.timerange = nanopores.logtimerange(t, levels=10, frac=frac, change_dt=pde.change_dt)
    pde.solve(t=t, visualize=True)
    
    #dolfin.interactive()

F, Fel, Fdrag, mesh = load_force_field()
geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **PARAMS)
convect(geo, phys, Fel, Fdrag)
  
"""  
class ForceCorrection(object):
    
    
def construct_correction():
    load_force_field()
    mesh = Howorka.loadedmesh
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    VV = dolfin.VectorFunctionSpace(mesh, "CG", 1)
    
    v2d = dolfin.vertex_to_dof_map(VV)
    coord = mesh.coordinates() # numpy array of 2-1 arrays
    linspace = numpy.linspace(*space)
"""