from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *

# @Benjamin, Gregor TODO:
# -) check permittivity and surface charge of ahem
# -) what biased voltage to use?

geo_params = dict(
    l4 = 15.,
    R = 20.,
    x0 = [2., 4., 6.],
)
phys_params = dict(
    bV = 1.0,
    ahemqs = 0.01,
    rpermProtein = 2.,
)    

t = Timer("meshing")
meshdict = generate_mesh(5., "aHem", **geo_params)

print "Mesh generation time:",t.stop()
print "Mesh file:",meshdict["fid_xml"]
print "Mesh metadata:"
for item in meshdict["meta"].items():
    print "%s = %s" %item
print 

t = Timer("reading geometry")
geo = geo_from_xml("aHem")

print "Geo generation time:",t.stop()
print "Geo params:", geo.params
print "Geo physical domains:", geo._physical_domain
print "Geo physical boundaries:", geo._physical_boundary

plot(geo.boundaries)
plot(geo.submesh("solid"))
#plot(geo.submesh("exittime"))

phys = Physics("pore_molecule", geo, **phys_params)

#print "Physics:"
#for item in phys.__dict__.items():
#    print "%s = %s" %item
    
pde = PNPS(geo, phys)
pde.solve()
pde.print_results()

(v, cp, cm, u, p) = pde.solutions(deepcopy=True)
E = -phys.grad(v)
Q = -qq
pi = 3.141592
r = 0.5*nm
Fel = Q*E
Fdrag = 6*pi*eta*r*u
F = Fel + Fdrag

'''
VV = VectorFunctionSpace(geo.mesh, "CG", 1)
Fdrag = project(Fdrag, VV)
F = project(F, VV)
plot(F, title="F")
plot(Fdrag, title="Fdrag")
interactive()
'''
# solve exit time problem
phys.Dtarget = phys.kT/(6*pi*eta*r)

etp = LinearPDE(geo, ExitTimeProblem, phys, F=F)
etp.solve()
etp.visualize("exittime")

etp_noF = LinearPDE(geo, ExitTimeProblem, phys, F=Constant((0.0,0.0,0.0)))
etp_noF.solve()
etp_noF.visualize("exittime")

