" analytical test problem to validate 2D and 3D solvers "
import math
from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

# --- define parameters ---
bV = -0.05 # [V]
rho = -0.025 # [C/m**2]
initialh = .1
Nmax = 1e1

# --- create 2D geometry ---
Rz = 2. # [nm] length in z direction of channel part
R = 2. # [nm] pore radius
hcross = .2

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
    center = domain2D.boundary("left")
)
domain2D.params["lscale"] = 1e9
domain2D.synonymes = dict(
    fluid = {"main", "cross"},
    pore = "fluid",
    chargedmembraneb = "wall",
    noslip = "wall",
    nopressure = "center",
    bulk = {"upperb", "lowerb"},
    #nocbc = {"lowerb"},
)
geo2D = domain2D.create_geometry(lc=initialh)
#mesh = geo2D.mesh
#boundary = MeshFunction("size_t", mesh, 1)
#boundary.set_all(0)
#DomainBoundary().mark(boundary, 2)
#plot(boundary)
#domain2D.plot()

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
    bulkcon = 300,
    v0 = {}
)
phys = Physics("pore", geo1D, **phys_params)

# --- solve 1D problem for "exact" solution ---
pb = solve_pde(SimplePBProblem, geo1D, phys, cyl=True, iterative=False, tolnewton=1e-10)

# define expression for interpolation into 2D
phi = pb.solution
UT = phys.UT
c0 = phys.bulkcon
D = phys.DPore
lscale = phys.lscale
E0 = -lscale*bV/(2.*Rz)
eps = phys.permittivity["water"]
eta = phys.eta

print "Diffusion constant in pore:", D*1e9, "[nm**2/ns]"
print "Constant electric field:", E0, "[V/m]"
def cpPB(x):
    return c0*exp(-phi(x)/UT)
def cmPB(x):
    return c0*exp(phi(x)/UT)
def pPB(x):
    return -2.*c0*cFarad*UT*(math.cosh(phi(x)/UT) - math.cosh(phi(0.)/UT))
def uPB(x):
    return eps*E0/eta*(phi(x[0]) - phi(R))
    
class vPB(Expression):
    def eval(self, value, x):
        value[0] = bV*x[1]/(2.*Rz) + phi(x[0])
class JpPB(Expression):
    def eval(self, value, x):
        value[0] = D/UT*E0*cpPB(x[0])
class JmPB(Expression):
    def eval(self, value, x):
        value[0] = -D/UT*E0*cmPB(x[0])

# compute current
r2pi = Expression("2*pi*x[0]")
J_PB = assemble(Constant(cFarad*D/UT*E0*c0/lscale**2)*(exp(-phi/UT) + exp(phi/UT))*r2pi*dx)
print "J (PB): %s [A]" % J_PB


# --- define physical parameters and customized BCs of 2D problem ---

# constant Dirichlet BCs for v, cp, cm on wall,
# non-zero flux BCs on top/bottom
n = FacetNormal(geo2D.mesh)
lscale = Constant(phys.lscale)
phys_params.update(
    cp0 = dict(wall = c0*exp(-phi(R)/UT)),
    cm0 = dict(wall = c0*exp(+phi(R)/UT)),
    v0 = dict(wall = vPB()),
    cpflux = dict(bulk = JpPB()*n[1]),
    cmflux = dict(bulk = JmPB()*n[1]),
)
phys = Physics("pore", geo2D, **phys_params)
phys.surfcharge.update(
    upperb = lscale*eps*bV/(2.*Rz),
    lowerb = -lscale*eps*bV/(2.*Rz),
)

# --- solve 2D PNP+Stokes problem ---

# the goal functional: current through crosssection
grad = phys.grad
def J(U, geo):
    v, cp, cm = U
    Jp = Constant(D)*(-grad(cp) - Constant(1/UT)*cp*grad(v))
    Jm = Constant(D)*(-grad(cm) + Constant(1/UT)*cm*grad(v))
    Jsurf = avg(Constant(cFarad/lscale**2)*(Jp - Jm)[1] * r2pi) * geo.dS("cross")
    Jvol = Constant(cFarad/lscale**2/hcross)*(Jp - Jm)[1] * r2pi * geo.dx("cross")
    return dict(Jsurf=Jsurf, Jvol=Jvol)

def saveJ(self):
    self.save_estimate("(Jsing_h - J)/J", abs((self.functionals["Jsurf"].value()-J_PB)/J_PB), 
    N=self.solution.function_space().dim())
    self.save_estimate("(J_h - J)/J", abs((self.functionals["Jvol"].value()-J_PB)/J_PB), 
    N=self.solution.function_space().dim())

# solve    
pnps = solve_pde(SimplePNPProblem, geo2D, phys, cyl=True, newtondamp=1., goals=[J], inside_loop=saveJ, 
    refinement=True, marking_fraction=.5, maxcells=Nmax, iterative=False)
v, cp, cm = pnps.solutions()
stokes = solve_pde(SimpleStokesProblem, geo2D, phys, cyl=True, conservative=False, f=-cFarad*(cp-cm)*grad(v), ku=1, beta=10.)

# --- visualization ---
#plot(-cFarad*(cp-cm)*grad(v)[1]/(lscale**2*eta), title="electroosmotic forcing [m/s]")
#pnps.visualize()
#stokes.visualize()
(u, p) = stokes.solutions()

#fig = plot1D({"phi PB":phi}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "potential [V]"))
#plot1D({"phi PNP (2D)": v}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "potential [V]"), fig=fig)

#fig = plot1D({"c+ PB":cpPB, "c- PB":cmPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "concentration [mol/m**3]"))
#plot1D({"c+ PNP (2D)": cp, "c- PNP (2D)": cm}, (0., R, 101), "x", origin=(0.,-Rz), dim=2, axlabels=("r [nm]", "concentration [mol/m**3]"), fig=fig)

#plot1D({"c+ PNP (2D)": cp, "c- PNP (2D)": cm, "c+ PB":lambda x: cpPB(0.), "c- PB":lambda x: cmPB(0.)}, 
#    (-Rz, Rz, 101), "y", origin=(.0*R, 0.), dim=2, axlabels=("z [nm]", "concentration [mol/m**3]"))

fig = plot1D({"uz PB":uPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
fig = plot1D({"uz PNP (2D)":u[1]}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "velocity [m/s]"), fig=fig)
#fig = plot1D({"ur PB":lambda x:0.}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
#fig = plot1D({"ur PNP (2D)":u[0]}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "velocity [m/s]"), fig=fig)
fig = plot1D({"p PB":pPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
fig = plot1D({"p PNP (2D)":p}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "velocity [m/s]"), fig=fig)

#pnps.estimators["(Jsing_h - J)/J"].plot(rate=-1.)
#pnps.estimators["(J_h - J)/J"].plot(rate=-1.)
showplots()

