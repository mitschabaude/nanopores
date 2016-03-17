''' this script is designed to determine the electric field induced near the nanopore by a cylindrical upper liquid reservoir filled with charged particles (fluorophores), given the concentration, single molecule charge, and dimensions of the reservoir.
the charge is smeared out across the liquid.
'''

from nanopores import *
from dolfin import *

# ----- physical parameters -----
mm = 1e-3
lx = 10*mm # radius [m]
ly = 2*10*mm # height [m]
uM = 1e-3
bulkcon = 10*uM # bulk concentration [mM] [mol/m**3]
Qmol = -3*qq # charge per molecule

# ----- discretization parameters -----
N = 60 # number of points in every space direction

def surfcharge(Qmol, bulkcon=10e-3, h=10e-3):
    return bulkcon*mol*Qmol*h/2
    
print "surf charge",surfcharge(Qmol, bulkcon, ly/2)
if not __name__=="__main__":
    exit()

# volume charge [C/m**3]
volcharge0 = bulkcon*mol*Qmol
print "volume charge [C/m**3]:",volcharge0

mesh = RectangleMesh(Point(0.,0.), Point(lx, ly), N, N)

class Top(SubDomain):
    def inside(self, x, on):
        return near(x[1], ly)
        
class Bottom(SubDomain):
    def inside(self, x, on):
        return near(x[1], 0.)
        
class CenterLine(SubDomain):
    def inside(self, x, on):
        return near(x[1], ly/2.)
    
class Left(SubDomain):
    def inside(self, x, on):
        return near(x[0], 0.)
        
class Center(SubDomain):
    def inside(self, x, on):
        return near(x[1], ly/2) and near(x[0], 0.)
        
class CenterBlock(SubDomain):
    def inside(self, x, on):
        return between(x[0], (0.,lx/2.)) and between(x[1], (ly/3., 2.*ly/3.))
        
class TopBlock(SubDomain):
    def inside(self, x, on):
        return between(x[1], (ly/2., ly))

subdomains = CellFunction("size_t", mesh, 0)       
boundaries = FacetFunction("size_t", mesh, 0)
TopBlock().mark(subdomains, 1)
CenterLine().mark(boundaries, 1)
Top().mark(boundaries, 2)
Bottom().mark(boundaries, 3)
physical_domain = {"center":(1,), "rest":(0,)}
physical_boundary = {"centerline":(1,), "top":(2,), "bottom":(3,)}
synonymes = {
    "fluid":{"center","rest"},
    "ground":set(),
    "bV":set()
    }
volcharge = {
    "center": volcharge0,
    "rest": 0.
    }

geo = Geometry(
    mesh = mesh,
    boundaries = boundaries,
    subdomains = subdomains,
    physical_boundary = physical_boundary,
    physical_domain = physical_domain,
    synonymes = synonymes,
    )
    
phys = Physics(
    geo = geo,
    surfcharge = {},
    volcharge = volcharge,
    bV = None)


area = assemble(Expression("x[0]*2*pi")("+")*geo.dS("centerline"))    
print "Total volume charge: %g [C]" %assemble(geo.pwconst("volcharge")*Expression("x[0]*2*pi")*geo.dx())
print "Induced surface charge: %g [C/m**2]" %(assemble(volcharge0*Expression("x[0]*2*pi")*geo.dx("center"))/area/2,)
   
PoissonProblemPureNeumannAxisym.method["iterative"] = False
#V = PoissonProblemAxisym.space(mesh)
#bc = DirichletBC(V, Constant(0.), Center(), method="pointwise")
#pde = LinearPDE(geo, PoissonProblemAxisym, phys=phys, bcs=bc)
#pde = LinearPDE(geo, PoissonProblemAxisym, phys=phys)
pde = LinearPDE(geo, PoissonProblemPureNeumannAxisym, phys=phys) #, bcs=bc)
pde.solve()
(u, p) = pde.solution.split()
#u = pde.solution

eps = geo.pwconst("permittivity")

Eline = assemble((-grad(u)[1]*Expression("x[0]*2*pi"))("+")*geo.dS("centerline"))/area
h = 1e-1*mm
Epoint = -(u([0.,ly/2+h]) - u([0.,ly/2-h]))/(2.*h)

Cline = assemble((eps*grad(u)[1]*Expression("x[0]*2*pi"))("+")*geo.dS("centerline"))/area
h = 1e-1*mm
Cpoint = phys.permittivity["fluid"]*(u([0.,ly/2+h]) - u([0.,ly/2-h]))/(2.*h)

print "Electric field average over center line: %g [V/m]" %Eline
print "Electric field at center point: %g [V/m]" %Epoint

print "Equivalent surface charge, average over center line: %g [C/m**2]" %Cline
print "Equivalent surface charge at center point: %g [C/m**2]" %Cpoint

'''
plot(-grad(u), title="E")
plot(u, title="u")
interactive()
'''
