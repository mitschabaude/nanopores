#FIXME: the current still does not converge, but the potential seems to do.
# plot the other solutions as well! where is the solution not accurate?

" analytical test problem to validate 2D and 3D solvers "
import math
from collections import OrderedDict
from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *
from nanopores.geometries.curved import Cylinder

# --- define parameters ---
add_params(
bV = -0.1, # [V]
rho = -0.05, # [C/m**2]
h2D = 1.,
h3D = 1.,
nmax2D = 2e3,
nmax3D = 2e4,
damp = 1.,
bulkcon = 300.,
iterative = True,
inewton = 1,
ipicard = 30,
)
print PARAMS
set_log_level(100)

# --- create 2D geometry ---
Rz = 2. # [nm] length in z direction of channel part
R = 1. # [nm] pore radius
hcross = 1. # [nm] height of crossection

domain2D = Box([0., -Rz], [R, Rz])
cross = Box([0., 0.], [R, hcross])
domain2D.addsubdomains(
    main = domain2D - cross,
    cross = cross,
)
domain2D.addboundaries(
    lowerb = domain2D.boundary("bottom"),
    upperb = domain2D.boundary("top"),
    wall = domain2D.boundary("right"),
    cross = cross.boundary("bottom"),
    #center = domain2D.boundary("left")
)
domain2D.params["lscale"] = 1e9
domain2D.synonymes = dict(
    fluid = {"main", "cross"},
    pore = "fluid",
    #chargedmembraneb = "wall",
    noslip = "wall",
    nopressure = "center",
    bulk = {"upperb", "lowerb"},
    #nocbc = {"lowerb"},
)

geo2D = domain2D.create_geometry(lc=h2D)
print "Number of cells (2D):", geo2D.mesh.num_cells()
#domain2D.plot()

# --- create 3D geometry by rotating ---
domain3D = rotate_z(domain2D)
geo3D = domain3D.create_geometry(lc=h3D)
print "Number of cells (3D):", geo3D.mesh.num_cells()
#domain3D.plot()
plot(geo3D.subdomains, title="initial")

# define cylinder
cyl = Cylinder(R=R, L=2.*Rz)
# this causes geo to automatically snap boundaries when adapting
geo3D.curved = dict(wall = cyl.snap)

# --- create geometry for 1D crossection ---
# TODO: it would be cool if the 1D domain could be extracted from the 2D one
# (should be pretty easy)
domain1D = Interval(0., R)
domain1D.addsubdomain(domain1D, "fluid")
domain1D.addboundaries(
    wall = domain1D.boundary("right"),
    center = domain1D.boundary("left")
)
domain1D.params["lscale"] = 1e9
domain1D.synonymes = dict(
    pore = "fluid",
    chargedmembraneb = "wall",
)
geo1D = domain1D.create_geometry(lc=.001) 

# --- define physical parameters for 1D problem ---
phys_params = dict(
    Membraneqs = rho,
    bulkcon = bulkcon,
    v0 = {}
)
phys = Physics("pore", geo1D, **phys_params)

# --- solve 1D problem for "exact" solution ---
pb = solve_pde(SimplePBProblem, geo1D, phys, cyl=True, iterative=False, tolnewton=1e-10)

# define expressions for interpolation into 2D/3D
phi = pb.solution
UT = phys.UT
c0 = phys.bulkcon
D = phys.DPore
lscale = phys.lscale
E0 = -lscale*bV/(2.*Rz)
eps = phys.permittivity["water"]
eta = phys.eta

print "Diffusion constant in pore:", D*1e9, "[nm**2/ns]"
print "Constant electric field:", E0*1e-9, "[V/nm]"
print "Debye length:", phys.debye*1e9, "[nm]"

def cpPB(x):
    return c0*math.exp(-phi(x)/UT)
def cmPB(x):
    return c0*math.exp(phi(x)/UT)
def pPB(x):
    return -2.*c0*cFarad*UT*(math.cosh(phi(x)/UT) - math.cosh(phi(0.)/UT))
def uPB(x):
    return eps*E0/eta*(phi(x) - phi(R))
    
def r(x): # radius for 2D AND 3D:
    return sqrt(sum(t**2 for t in x[:-1]))
    
class vPB(Expression):
    def eval(self, value, x):
        value[0] = bV*x[-1]/(2.*Rz) + phi(r(x))
class JpPB(Expression):
    def eval(self, value, x):
        value[0] = (+D/UT*E0 + uPB(r(x)))*cpPB(r(x))
class JmPB(Expression):
    def eval(self, value, x):
        value[0] = (-D/UT*E0 + uPB(r(x)))*cmPB(r(x))
class cpPBEx(Expression):
    def eval(self, value, x):
        value[0] = cpPB(r(x))
class cmPBEx(Expression):
    def eval(self, value, x):
        value[0] = cmPB(r(x))
class pPBEx(Expression):
    def eval(self, value, x):
        value[0] = pPB(r(x))

# compute "exact" current
r2pi = Expression("2*pi*x[0]")
u_PB = Constant(eps/eta)*(phi - Constant(phi(R)))
J_el = Constant(D/UT)*(exp(-phi/UT) + exp(phi/UT))
J_u = u_PB*(exp(-phi/UT) - exp(phi/UT))
J_PB_el = assemble(Constant(cFarad*c0*E0/lscale**2)*J_el*r2pi*dx)
J_PB_u = assemble(Constant(cFarad*c0*E0/lscale**2)*J_u*r2pi*dx)
J_PB = J_PB_el + J_PB_u
print "J (PB): %s [A]" % J_PB
print "   J_el: %s [A]" % J_PB_el
print "   J_u : %s [A]" % J_PB_u

# --- define physical parameters and customized BCs of 2D problem ---

# constant Dirichlet BCs for v, cp, cm on wall,
# non-zero flux BCs on top/bottom
# non-standard pressure BC
lscale = Constant(phys.lscale)
phys_params.update(
    cp0 = dict(
        wall = cpPBEx(), #c0*exp(-phi(R)/UT),
    #    bulk = cpPBEx(),
        ),
    cm0 = dict(
        wall = cmPBEx(), #c0*exp(+phi(R)/UT),
    #    bulk = cmPBEx(),
        ),
    v0 = dict(wall = vPB()),
    cpflux = dict(upperb = +JpPB(),
                  lowerb = -JpPB(),),
    cmflux = dict(upperb = +JmPB(),
                  lowerb = -JmPB(),),
    pressure = dict(bulk = pPBEx()),
    surfcharge = dict(
        wall = rho,
        upperb = lscale*eps*bV/(2.*Rz),
        lowerb = -lscale*eps*bV/(2.*Rz),)
)
phys2D = Physics("pore", geo2D, **phys_params)
phys3D = Physics("pore", geo3D, **phys_params)

# --- define goal functional: current through crosssection ---
grad = phys2D.grad

def J(U, geo):
    v, cp, cm, u, p = U
    #u, p, v, cp, cm = U
    dim = geo.physics.dim
    coeff = Constant(1.) if dim==3 else r2pi
    Jp = Constant(D)*(-grad(cp) - Constant(1/UT)*cp*grad(v)) + cp*u
    Jm = Constant(D)*(-grad(cm) + Constant(1/UT)*cm*grad(v)) + cm*u
    Jsurf = avg(Constant(cFarad/lscale**2)*(Jp - Jm)[dim-1] * coeff) * geo.dS("cross")
    Jvol = Constant(cFarad/lscale**2/hcross)*(Jp - Jm)[dim-1] * coeff * geo.dx("cross")
    return dict(Jsurf=Jsurf, Jvol=Jvol)

def saveJ(self):
    i = sum(u.function_space().dim() for u in self.solutions())
    (v, cp, cm, u, p) = self.solutions()
    if self.geo.physics.dim == 3:
        #plot1D({"phi (3D), N=%s" %i: v}, (0., R, 101), "x", dim=3, axlabels=("r [nm]", "potential [V]"), newfig=False)
        #plot1D({"p (3D), N=%s" %i:p}, (0., R, 101), "x", dim=3, axlabels=("r [nm]", "pressure [Pa s]"), newfig=False)
        #plot1D({"uz (3D), N=%s" %i:u[2]}, (0., R, 101), "x", dim=3, axlabels=("r [nm]", "velocity [m/s]"), newfig=False)
        plot1D({"c+ (3D), N=%s" %i: cp}, (-Rz, Rz, 101), "z", dim=3, 
            axlabels=("z [nm]", "concentration [mol/m**3]"), newfig=False)

    self.save_estimate("(J*h - J)/J" + Dstr, abs((self.functionals["Jsurf"].evaluate()-J_PB)/J_PB), N=i)
    self.save_estimate("(Jh - J)/J" + Dstr, abs((self.functionals["Jvol"].evaluate()-J_PB)/J_PB), N=i)
    print "     rel. error Jv:", abs((self.functionals["Jvol"].value()-J_PB)/J_PB)
    print "     rel. error Js:", abs((self.functionals["Jsurf"].value()-J_PB)/J_PB)

# --- set up PNP+Stokes problem ---

problems = OrderedDict([
    ("pnp", SimplePNPProblem),
    ("stokes", SimpleStokesProblem),])

def couple_pnp(ustokes):
    return dict(ustokes = ustokes.sub(0))

def couple_stokes(upnp, phys):
    v, cp, cm = upnp.split()
    f = -phys.cFarad*(cp - cm)*grad(v)
    return dict(f = f)

couplers = dict(
    pnp = couple_pnp,
    stokes = couple_stokes
)

# --- solve 2D problem ---
Dstr = " (2D)"
problem = CoupledProblem(problems, couplers, geo2D, phys2D, cyl=True, conservative=False, ku=1, beta=10.)
pnps2D = CoupledSolver(problem, goals=[J], damp=damp, inewton=1, ipicard=20, tolnewton=1e-2, verbose=False)
#pnps2D.uniform_refinement = True
pnps2D.marking_fraction = 1.
pnps2D.maxcells = nmax2D
pnps2D.solve(refinement=True, inside_loop=saveJ)

# --- solve 3D problem ---
Dstr = " (3D)"
problem = CoupledProblem(problems, couplers, geo3D, phys3D, cyl=False, conservative=False, ku=1, beta=1.)
problem.problems["pnp"].method["iterative"] = iterative
problem.problems["stokes"].method["iterative"] = iterative
pnps3D = CoupledSolver(problem, goals=[J], damp=damp, inewton=inewton, ipicard=ipicard, tolnewton=1e-2, verbose=True, nverbose=True)
#pnps3D.uniform_refinement = True
pnps3D.marking_fraction = 1.
pnps3D.maxcells = nmax3D
#plot1D({"phi PB":phi}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "potential [V]"))
#plot1D({"p PB":pPB}, (0., R, 101), "x", dim=1)
#plot1D({"uz PB":uPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
plot1D({"c+ PB":lambda x: cpPB(0.)}, (-Rz, Rz, 101), "y", dim=2, 
            axlabels=("z [nm]", "concentration [mol/m**3]"))
pnps3D.solve(refinement=True, inside_loop=saveJ)

# --- visualization ---
pnps2D.estimators["(Jh - J)/J (2D)"].plot(rate=-1.)
pnps3D.estimators["(Jh - J)/J (3D)"].plot(rate=-2./3, fig=False)

pnps3D.visualize("cross")
#saveplots("anaPNPSrefine", meta=PARAMS)
showplots()

