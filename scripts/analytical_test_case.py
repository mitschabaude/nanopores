" analytical test problem to validate 2D solver "
import math
from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

# --- define parameters ---
bV = -0.5 # [V]
rho = -0.025 # [C/m**2]

# --- create 2D geometry ---
Rz = 2. # [nm] length in z direction of channel part
R = 2. # [nm] pore radius

domain2D = Box([0., -Rz], [R, Rz])
cross = Box([0., 0.], [R, 0.])
domain2D.addsubdomains(fluid = domain2D)
domain2D.addboundaries(
    lowerb = domain2D.boundary("bottom"),
    upperb = domain2D.boundary("top"),
    wall = domain2D.boundary("right"),
    cross = cross,
)
domain2D.params["lscale"] = 1e9
domain2D.synonymes = dict(
    pore = "fluid",
    chargedmembraneb = "wall",
)
geo2D = domain2D.create_geometry(lc=.1)
domain2D.plot()

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
geo1D = domain1D.create_geometry(lc=.01) 

# --- define physical parameters for 1D problem ---
phys_params = dict(
    Membraneqs = rho,
    bulkcon = 300,
    v0 = {}
)
phys = Physics("pore", geo1D, **phys_params)

# --- solve 1D problem for "exact" solution and interpolate onto 2D mesh ---
pb = solve_pde(SimplePBProblem, geo1D, phys, cyl=True, iterative=False, tolnewton=1e-10)
phi = pb.solution

UT = phys.UT
c0 = phys.bulkcon
D = phys.DPore
lscale = phys.lscale
E0 = -lscale*bV/(2.*Rz)
print("Diffusion constant in pore:",D)
print("Constant electric field:",E0)

def cpPB(x):
    return c0*exp(-phi(x)/UT)
def cmPB(x):
    return c0*exp(phi(x)/UT)

class Potential(Expression):
    def eval(self, value, x):
        value[0] = phi(x[0]) + bV*x[1]/(2.*Rz)
class Jp(Expression):
    def eval(self, value, x):
        value[0] = D/UT*E0*cpPB(x[0])
class Jm(Expression):
    def eval(self, value, x):
        value[0] = -D/UT*E0*cmPB(x[0])
"""
phi1 = Function(FunctionSpace(geo2D.mesh, 'CG', 2))
phi1.interpolate(Potential())
plot(phi1, interactive=True)
phi1.interpolate(Jp())
plot(phi1, interactive=True)
phi1.interpolate(Jm())
plot(phi1, interactive=True)
"""

# --- define physical parameters and non-standard BCs for 2D problem ---

#v0 = dict(wall = phi(R))
v0ex = Potential()
v0 = dict(
    upperb = v0ex,
    lowerb = v0ex,
    wall = v0ex,
)
cp0 = dict(wall = c0*exp(-phi(R)/UT))
cm0 = dict(wall = c0*exp(+phi(R)/UT))

phys_params.update(
    cp0 = cp0,
    cm0 = cm0,
    v0 = v0,
)
phys = Physics("pore", geo2D, **phys_params)

V = SimplePNPProblem.space(geo2D.mesh)
bcs = geo2D.pwBC(V.sub(0), "v0")
bcs = bcs + geo2D.pwconstBC(V.sub(1), "cp0")
bcs = bcs + geo2D.pwconstBC(V.sub(2), "cm0")

# --- create customized problem with non-zero flux BCs for nernst-planck ---
problem = SimplePNPProblem(geo2D, phys, cyl=True, newtondamp=0.8, bcs=bcs)

w, dp, dm = split(problem.a.arguments()[0])
r2pi = Expression("2*pi*x[0]")
lscale = Constant(phys.lscale)
Lp = -lscale*Jp()*dp*r2pi*geo2D.ds("upperb") + lscale*Jp()*dp*r2pi*geo2D.ds("lowerb")
Lm = -lscale*Jm()*dm*r2pi*geo2D.ds("upperb") + lscale*Jm()*dm*r2pi*geo2D.ds("lowerb")
problem.addRHS(- Lp - Lm)

# --- solve 2D problem
pnps = solve_problem(problem, geo2D)
pnps.visualize()
(v, cp, cm) = pnps.solutions()

#print geo2D
#domain2D.plot()

#print geo1D

fig = plot1D({"PB":phi}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "potential [V]"))
plot1D({"PNP (2D)": v}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "potential [V]"), fig=fig)

fig = plot1D({"c+ PB":cpPB, "c- PB":cmPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "concentration [mol/m**3]"))
plot1D({"c+ PNP (2D)": cp, "c- PNP (2D)": cm}, (0., R, 101), "x", origin=(0.,-Rz), dim=2, axlabels=("r [nm]", "concentration [mol/m**3]"), fig=fig)
showplots()



