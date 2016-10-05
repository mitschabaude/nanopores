import dolfin
import nanopores as nano
import nanopores.tools.box as box
import nanopores.geometries.pughpore as pughpore
from nanopores.geometries.curved import Sphere
#from nanopores.models.mysolve import pbpnps
import nanopores.physics.simplepnps as simplepnps
from nanopores.models.mysolve import mesh_quality

up = nano.user_params(
    h = 2.,
    Qmol = -1.,
    Nmax = 1e6,
)

geop = nano.Params(
    R = 35.,
    H = 70.,
    x0 = [0.,0.,8.]
)
physp = nano.Params(
    Qmol = up.Qmol,
    bulkcon = 300.,
    dnaqsdamp = 0.2,
    bV = -0.05,
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
print geo
form = dolfin.Constant(1.)*geo.dx("molecule")
print form
print dolfin.assemble(form)
#nano.plot_sliced(geo)
#dolfin.interactive()
molec = Sphere(R=geo.params["rMolecule"], center=geo.params["x0"])
geo.curved = dict(moleculeb = molec.snap)

phys = nano.Physics("pore_mol", geo, **physp)
print phys.surfcharge
print phys.volcharge
R, H = float(geo.params["R"]), float(geo.params["H"])
mesh2D = nano.RectangleMesh([0,-H/2.], [R, H/2.], int(4*R), int(4*H))

simplepnps.SimpleStokesProblem.method["kparams"].update(
    monitor_convergence = False,
    relative_tolerance = 1e-5,
    absolute_tolerance = 1e-5,
    maximum_iterations = 1000,
    nonzero_initial_guess = True,
    ) 
#simplepnps.SimpleLinearPBProblem.method["iterative"] = False
simplepnps.SimpleLinearPBProblem.method["kparams"].update(
    relative_tolerance = 1e-10,
    absolute_tolerance = 1e-6,
    monitor_convergence = False,
    #nonzero_initial_guess = True,
)

dolfin.tic()
# TODO
goal = phys.CurrentPB
pb = simplepnps.SimpleLinearPBGO(geo, phys, goal=goal,
                                 cheapest=solverp.cheapest)
pb.maxcells = solverp.Nmax
pb.marking_fraction = solverp.frac
refined = True
i = 0

print "Number of cells:", pb.geo.mesh.num_cells()
while refined:
    i += 1
    print "\nAssessing mesh quality."
    mesh_quality(pb.geo.mesh, ratio=0.01, geo=pb.geo, plothist=False)
    print "\nSolving PB."
    pb.single_solve()
    nano.plot_cross(pb.functions["primal"], mesh2D, title="pb potential", key="pb")
    print "\nError estimation."
    (ind, err) = pb.estimate()
    print "\nMesh refinement."
    refined = pb.refine(ind)
    if not refined:
        print "Maximal number of cells reached."
    else:
        print "New total number of cells:", pb.geo.mesh.num_cells()

#pb = nano.solve_pde(simplepnps.SimplePBProblem, geo, phys,
#                    iterative=True, tolnewton=1e-2)
print "CPU time (PB): %.3g s" %(dolfin.toc(),)
#pb.visualize("pore")
#dolfin.interactive()
#pnps = simplepnps.PNPSHybrid(geo, phys, v0=pb.solution, pscale=1e7,
#         iterative=True, inewton=1, ipicard=30, tolnewton=solverp.tol,
#         verbose=True)

pnps = simplepnps.PNPSFixedPointbV(geo, phys, ipicard=solverp.imax,
           #ku=2, kp=1, beta=0.,
           pscale=1e7, stokesiter=True, v0=pb.solution, 
           tolnewton=solverp.tol, verbose=True, iterative=True)          

print "Number of cells:", geo.mesh.num_cells()
print "DOFs:", pnps.dofs()

# solve PNPS
dolfin.tic()
for i in pnps.fixedpoint(ipnp=5):
    if isinstance(pnps, simplepnps.PNPSHybrid):
        v = pnps.functions["pnp"].sub(0)
    else:
        v = pnps.functions["poisson"]
    u = pnps.functions["stokes"].sub(0)
    nano.plot_cross(v, mesh2D, title="potential", key="u")
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

pnps.visualize("pore")
nano.showplots()
#dolfin.interactive()
