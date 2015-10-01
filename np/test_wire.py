from nanopores import generate_mesh, geo_from_name, Poisson, PoissonProblem, qq, PB, PBProblem, import_vars, params_physical, get_mesh, u_to_matlab, IllposedNonlinearSolver
from dolfin import *
from random import random

t = Timer("time")
geo_name = "Nanowire"
globals().update(import_vars("nanopores.%s.params_geo" % geo_name))

def random_dopant_positions(N, **params):
    # create N random points in unit cube
    x = [[random() for i in range(3)] for i in range(N)]
    # affine transformation on transducer
    rr = params["r_eff"]
    ll = params["l"]
    T = lambda x : [x[0]*(w_core - 2*rr) + x_core + rr, x[1]*(ll - 2*rr) + rr,
                    x[2]*(h_core - 2*rr) + z_core + rr]
    return map(T, x)

N = 5 # number of impurity atoms in pore
l = 10 # length of wire [nm]
clscale = 1.0 # mesh width scale
Nsamples = 1 # number of samples
r_eff = 3 # effective radius [nm]
voltage = -3 # voltage drop across 1000 nm wire [V]

# change physical parameters
#params_physical.permanent_charge["impurity"] = 1e6
#params_physical.permittivity["impurity"] = ???
#params_physical.ion_concentration["impurity"] = ???

geo_params = dict(l = l*1e-9, r_eff = r_eff*1e-9, lc = 10e-9)
generate_mesh(clscale, geo_name, **geo_params)  # constructs the mesh
mesh = get_mesh(geo_name)

print "CPU Time (mesh generation):",t.stop()
print "Number of cells:",mesh.num_cells()

PB.adapt = PB.rebuild
PB.maxcells = 5000000
PB.marking_fraction = 0.9
PB.tolnewton = 1e-5
#IllposedNonlinearSolver.newtondamp = 0.9

u0 = Function(PBProblem.space(mesh))
    
for k in range(Nsamples):
    print "\n --- Sample %d of %d" %(k+1, Nsamples)
    
    # xi = list of dopant positions
    xi = random_dopant_positions(N, **geo_params)
    geo = geo_from_name(geo_name, mesh=mesh, check_midpoint=True, xi=xi, **geo_params)
    t.start()
    pde = PB(geo, bV=-1.0)
    #rho = geo.pwconst("permanent_charge")
    #pde = Poisson(geo, bV=-0.0, f=rho)
    
    pde.solve(refinement = False)
    print "CPU Time (solving):",t.stop()
    
    (u,) = pde.solutions()
    u0.vector()[:] = u0.vector()[:] + u.vector()[:]/Nsamples
    
#mesh_core = geo.submesh("core")
#geo_core = geo_from_name("Nanowire", mesh=mesh_core, check_midpoint=True, xi=xi, **geo_params)
#plot(geo_core.boundaries)

# compute current
from nanopores.physics.params_physical import UT, qq
ni = 1e16
mup = 100*1e4
mun = 1000*1e4
phiF = params_physical.fermi_level["sil"]

E = -voltage/1000e-9
dx = geo.dx("core")
p = ni*exp(-(u0 - phiF)/UT)
n = ni*exp((u0 - phiF)/UT)

Jp = qq*mup*p*E
Jn = qq*mun*n*E
j0 = Jp + Jn

I = assemble(j0*dx)/l
print "current: ", I

#j1 = Constant(qq*E*ni/l)*(mup*exp(-(u0 - phiF)/UT) + mun*exp((u0 - phiF)/UT))
#I1 = assemble(j1*dx)


#u_to_matlab(mesh, u0, "potential")

#plot(u, title="example potential")
plot(u0, title="average potential")
#pde.functions["PBProblem"] = u0
#pde.visualize("core")
interactive()


