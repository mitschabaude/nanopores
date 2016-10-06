import dolfin
import nanopores as nano
import nanopores.tools.box as box
import nanopores.geometries.pughpore as pughpore
from nanopores.geometries.curved import Sphere
#from nanopores.models.mysolve import pbpnps
import nanopores.physics.simplepnps as simplepnps
import solvers

up = nano.user_params(
#up = nano.Params( # for iPython
    h = 2.,
    Qmol = -1.,
    Nmax = 1e5,
    R = 35.,
)

geop = nano.Params(
    R = up.R,
    H = 70.,
    x0 = [0.,0.,15.]
)
physp = nano.Params(
    Qmol = up.Qmol,
    bulkcon = 300.,
    dnaqsdamp = .5,
    bV = -0.025,
)
solverp = nano.Params(
    h = up.h,
    frac = 0.2,
    Nmax = up.Nmax,  
    imax = 30,
    tol = 1e-2,
    cheapest = False,
)

# setup geometry and physics
box.set_tol(None)
geo = pughpore.get_geo(solverp.h, **geop)
#print geo
#nano.plot_sliced(geo)
#dolfin.interactive()
molec = Sphere(R=geo.params["rMolecule"], center=geo.params["x0"])
geo.curved = dict(moleculeb = molec.snap)

phys = nano.Physics("pore_mol", geo, **physp)
solvers.set_sideBCs(phys, geop, physp)

R, H = float(geo.params["R"]), float(geo.params["H"])
print R, H
mesh2D = nano.RectangleMesh([-R,-H/2.], [R, H/2.], int(8*R), int(4*H))

#simplepnps.SimpleStokesProblem.method["kparams"].update(
#    monitor_convergence = False,
#    relative_tolerance = 1e-5,
#    absolute_tolerance = 1e-5,
#    maximum_iterations = 1000,
#    nonzero_initial_guess = True,
#    ) 
#simplepnps.SimpleLinearPBProblem.method["iterative"] = False
#simplepnps.SimpleLinearPBProblem.method["kparams"].update(
#    relative_tolerance = 1e-10,
#    absolute_tolerance = 1e-6,
#    monitor_convergence = False,
#    #nonzero_initial_guess = True,
#)

# prerefine with PB
dolfin.tic()
goal = phys.CurrentPB
pb = simplepnps.SimpleLinearPBGO(geo, phys, goal=goal,
                                 cheapest=solverp.cheapest)
print pb.solvers["primal"].problem.bcs
print [bc.g([-R,0.,0.]) for bc in pb.solvers["primal"].problem.bcs]

for i in pb.adaptive_loop(solverp.Nmax, solverp.frac):
    nano.plot_cross(pb.solution, mesh2D,title="pb potential", key="pb")

print "CPU time (PB): %.3g s" %(dolfin.toc(),)
nano.plot1D(dict(pbx=pb.solution), (-R, R, 1001), axis="x",
            dim=3, axlabels=("x [nm]", "pb potential [V]"))
nano.plot1D(dict(pby=pb.solution), (-R, R, 1001), axis="y",
            dim=3, axlabels=("x [nm]", "pb potential [V]"), newfig=False)
            
nano.plot1D(dict(pbleft=pb.solution), (-H, H, 1001), axis="z",
            origin=(-R-0.1,0.,0.),
            dim=3, axlabels=("x [nm]", "pb potential [V]"))
nano.plot1D(dict(pbfront=pb.solution), (-H, H, 1001), axis="z",
            origin=(0,-R-0.1,0.),
            dim=3, axlabels=("x [nm]", "pb potential [V]"), newfig=False)

pb.visualize("fluid")
nano.showplots()
#dolfin.interactive()

pnps = simplepnps.PNPSFixedPoint(geo, phys, ipicard=solverp.imax,
           #taylorhood=True,
           stokesiter=True, v0=pb.solution, 
           tolnewton=solverp.tol, verbose=True, iterative=True)          

print "Number of cells:", geo.mesh.num_cells()
print "DOFs:", pnps.dofs()

# solve PNPS
dolfin.tic()
for i in pnps.fixedpoint(ipnp=5):
    v, cp, cm, u, p = pnps.solutions()
#    if isinstance(pnps, simplepnps.PNPSHybrid):
#        v, cp, cm = pnps.functions["pnp"].split()
#    else:
#        v = pnps.functions["poisson"]
    nano.plot_cross(v, mesh2D, title="potential", key="u")
    nano.plot_cross(cm, mesh2D, title="negative ions", key="cm")
    #u = pnps.functions["stokes"].sub(0)
    #nano.plot_cross_vector(u, mesh2D, title="u", key="uu")
    print "v0 =", v([0., 0., -25.])
    
print "CPU time (solve): %.3g s" %(dolfin.toc(),)

print pnps.evaluate(phys.CurrentPNPS)
print pnps.evaluate(phys.ForcesPNPS)

# visualize
v, cp, cm, u, p = pnps.solutions()
#nano.plot_cross(v, mesh2D, title="potential")
nano.plot_cross(cm, mesh2D, title="cm")
nano.plot_cross(cp, mesh2D, title="cp")
nano.plot_cross(p, mesh2D, title="p")
nano.plot_cross_vector(u, mesh2D, title="u")

#if not solverp.cheapest:
#    pb.estimators["rep"].plot()
#    pb.estimators["err"].plot(rate=-2./3., fig=False)

#pnps.visualize("pore")
nano.showplots()
#dolfin.interactive()
