" analytical test problem to validate 2D and 3D solvers "
import math
from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

# --- define parameters ---
bV = -0.05 # [V]
rho = -0.025 # [C/m**2]
Nmax = 1e5

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
    #nocbc = {"upperb", "lowerb"},
)
geo2D = domain2D.create_geometry(lc=.2)
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
class Jp(Expression):
    def eval(self, value, x):
        value[0] = D/UT*E0*cpPB(x[0])
class Jm(Expression):
    def eval(self, value, x):
        value[0] = -D/UT*E0*cmPB(x[0])

        
#jm = Function(FunctionSpace(geo2D.mesh, 'CG', 4))
#jm.interpolate(Jm())
#jp = Function(FunctionSpace(geo2D.mesh, 'CG', 4))
#jp.interpolate(Jp())        
jp = Jp()
jm = Jm()

# compute current
r2pi = Expression("2*pi*x[0]")
J_PB = assemble(Constant(cFarad*D/UT*E0*c0/lscale**2)*(exp(-phi/UT) + exp(phi/UT))*r2pi*dx)
print "J (PB): %s [A]" % J_PB


# --- define physical parameters and non-standard BCs for 2D problem ---

#v0 = dict(wall = phi(R))
v0 = dict(wall = vPB())
cp0 = dict(wall = c0*exp(-phi(R)/UT))
cm0 = dict(wall = c0*exp(+phi(R)/UT))

phys_params.update(
    cp0 = cp0,
    cm0 = cm0,
    v0 = v0,
)
phys = Physics("pore", geo2D, **phys_params)
phys.surfcharge.update(
    upperb = lscale*eps*bV/(2.*Rz),
    lowerb = -lscale*eps*bV/(2.*Rz),
)

V = SimplePNPProblem.space(geo2D.mesh)
bcs = geo2D.pwBC(V.sub(0), "v0")
bcs = bcs + geo2D.pwconstBC(V.sub(1), "cp0")
bcs = bcs + geo2D.pwconstBC(V.sub(2), "cm0")

# --- solve customized 2D PNP problem with non-zero flux BCs for nernst-planck ---
problem = SimplePNPProblem(geo2D, phys, cyl=True, newtondamp=0.5, bcs=bcs)
"""
# neumann bc
w, dp, dm = split(problem.a.arguments()[0])
lscale = Constant(phys.lscale)
Lp = lscale*jp*dp*r2pi*geo2D.ds("upperb") - lscale*jp*dp*r2pi*geo2D.ds("lowerb")
Lm = lscale*jm*dm*r2pi*geo2D.ds("upperb") - lscale*jm*dm*r2pi*geo2D.ds("lowerb")
problem.addRHS(- Lp - Lm)
"""

# the goal functional: current through crosssection
v, cp, cm = problem.solution().split()
grad = phys.grad
Jp = Constant(D)*(-grad(cp) - Constant(1/UT)*cp*grad(v))
Jm = Constant(D)*(-grad(cm) + Constant(1/UT)*cm*grad(v))
J = avg(Constant(cFarad/lscale**2)*(Jp - Jm)[1] * r2pi) * geo2D.dS("cross")
Jvol = Constant(cFarad/lscale**2/hcross)*(Jp - Jm)[1] * r2pi * geo2D.dx("cross")
def saveJ(self):
    self.save_estimate("(Jsing_h - J)/J", (self.functionals["J"].value()-J_PB)/J_PB, N=self.solution.function_space().dim())
    self.save_estimate("(J_h - J)/J", (self.functionals["Jvol"].value()-J_PB)/J_PB, N=self.solution.function_space().dim())

# solve    
pnps = solve_problem(problem, geo2D, goals={"J": J, "Jvol": Jvol}, inside_loop=saveJ, 
    refinement=True, marking_fraction=1., maxcells=Nmax, iterative=False)
    
# --- solve 2D Stokes Problem --
plot(-cFarad*(cp-cm)*grad(v)/(lscale**2*eta), title="electroosmotoic forcing [m/s]")
stokes = solve_pde(SimpleStokesProblem, geo2D, phys, cyl=True, conservative=False, f=-cFarad*(cp-cm)*grad(v), ku=2, beta=0.)

# --- visualization ---
pnps.visualize()
stokes.visualize()
(v, cp, cm) = pnps.solutions()
(u, p) = stokes.solutions()

#fig = plot1D({"phi PB":phi}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "potential [V]"))
#plot1D({"phi PNP (2D)": v}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "potential [V]"), fig=fig)

fig = plot1D({"c+ PB":cpPB, "c- PB":cmPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "concentration [mol/m**3]"))
plot1D({"c+ PNP (2D)": cp, "c- PNP (2D)": cm}, (0., R, 101), "x", origin=(0.,-Rz), dim=2, axlabels=("r [nm]", "concentration [mol/m**3]"), fig=fig)

plot1D({"c+ PNP (2D)": cp, "c- PNP (2D)": cm, "c+ PB":lambda x: cpPB(0.), "c- PB":lambda x: cmPB(0.)}, 
    (-Rz, Rz, 101), "y", origin=(.0*R, 0.), dim=2, axlabels=("z [nm]", "concentration [mol/m**3]"))

#fig = plot1D({"uz PB":uPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
#fig = plot1D({"uz PNP (2D)":u[1]}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "velocity [m/s]"), fig=fig)
#fig = plot1D({"ur PB":lambda x:0.}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
#fig = plot1D({"ur PNP (2D)":u[0]}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "velocity [m/s]"), fig=fig)
#fig = plot1D({"p PB":pPB}, (0., R, 101), "x", dim=1, axlabels=("r [nm]", "velocity [m/s]"))
#fig = plot1D({"p PNP (2D)":p}, (0., R, 101), "x", dim=2, axlabels=("r [nm]", "velocity [m/s]"), fig=fig)

#pnps.estimators["(Jsing_h - J)/J"].plot(rate=-1.)
#pnps.estimators["(J_h - J)/J"].plot(rate=-1.)
showplots()

