''' record IV curve for open Howorka pore and compare with experiment '''

from nanopores import *
from dolfin import *
import json
import numpy
import matplotlib.pyplot as plt

def save(data, fname):
    with open('%s.txt' % fname, 'w') as f:
        f.write('\n'.join([str(s) for s in data]))

def load(fname):
    with open('%s.txt' % fname, 'r') as f:
        data = f.read().split('\n')
    return data

# a crucial parameter is R, i.e. the distance
# where the potential bias is applied
R = 30*nm
maxcells = 80e3

# get experimental data from csv file
csvfile = '../howorkaIV'
IV = numpy.genfromtxt(csvfile+'.csv', delimiter=',')
V0 = IV[:,0]
I0 = IV[:,1]

geo_name = "H_geo"
nm = 1e-9
geo_params = dict(
Ry = R,
Rx = R,
x0 = None,
)
phys_params = {"Membraneqs": -0.03, "dnaqsdamp":0.1}

# IMPORTANT: diffusion and stokes damping in pore
phys_params.update({"stokesdampPore":1.0, "rDPore":0.4,})

meshgen_dict = generate_mesh(0.8, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore_molecule", geo, **phys_params)

if geo.parameter("x0") is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)

#mesh = Mesh("/".join([DATADIR, geo_name, "mesh", "last_adapted_mesh.xml"]))
#geo = geo_from_name(geo_name, mesh=mesh, **geo_params)
PNPSAxisym.tolnewton = 1e-1
PNPSAxisym.maxcells = maxcells
PNPProblemAxisym.method["iterative"] = False
PNPProblemAxisym.method["kparams"]["monitor_convergence"] = False
IllposedNonlinearSolver.newtondamp = 1.0

phys.bV = 0.05
print(phys.charge)
goal = (lambda v : phys.Fbare(v, 1) + phys.CurrentPB(v)) if geo_params["x0"] else (lambda v : phys.CurrentPB(v))
pb = LinearPBAxisymGoalOriented(geo, phys, goal=goal)
pb.maxcells = maxcells
pb.marking_fraction = 0.5
pb.solve(refinement=True)
v0 = pb.solution

geo.mesh = pb.geo.mesh
#plot(geo.subdomains)
N = geo.mesh.num_cells()
hmin = geo.mesh.hmin()

I = []
F = []
G = []
for bVmV in V0:
    bV = bVmV*1e-3
    print("\nbV = %.0f [mV]\n" % bVmV)

    phys.bV = bV
    pnps = PNPSAxisym(geo, phys)
    pnps.solve(refinement=False, save_mesh=False, print_functionals=False, visualize=False)

    (v,cp,cm,u,p) = pnps.solutions()
    f_dict = pnps.get_functionals()

    Jcomponents = ["Jzdiff", "Jzdrift", "Jzstokes"]
    Jzdiff = phys.cFarad*phys.D*phys.rDPore*phys.grad(-cp+cm)[1] * phys.r2pi/geo.params["l0"]*geo.dx("pore")
    Jzdrift = phys.cFarad*phys.mu*phys.rDPore*(-cp-cm)*phys.grad(v)[1] * phys.r2pi/geo.params["l0"]*geo.dx("pore")
    Jzstokes = phys.cFarad*phys.stokesdampPore*(cp-cm)*u[1] * phys.r2pi/geo.params["l0"]*geo.dx("pore")
    for j in Jcomponents:
        f_dict.update({j: 1e12*assemble(locals()[j])})

    print(json.dumps(f_dict, indent=4, sort_keys=True))

    I_ = f_dict["Javgbtm"]
    I.append(I_)
    F.append(f_dict)

    V = (v([0.0, -10*nm]) - v([0.0, 10*nm]))
    print("I (current through pore center):",I_,"[pA]")
    print("V (transmembrane potential):",V,"[V]")
    print("conductance I/V:",I_/bV,"[pS]")

# plot(geo.mesh)
# plot(v)
# plot(cm-cp, title="cdiff")
# plot(cm+cp)
# interactive()
# exit()


#for s in ["Z","I","V"]:
#    save(vars()[s], s)

plt.figure(1)
plt.plot(V0,I0, 's-', label="Howorka et. al")
plt.plot(V0,I, 'o-', label="PNPS simulation")
plt.xlabel("voltage bias [mV]")
plt.ylabel("current [pA]")
plt.legend(loc='lower right')
fname = ("data/%s_N%sdnad%.2frD%.2f" %(csvfile, N, phys.dnaqsdamp, phys.rDPore) )
plt.savefig(fname+".eps", bbox_inches='tight')

plt.figure(2)
Jcomp = numpy.asarray([[f[j] for f in F] for j in Jcomponents])
for j in range(len(Jcomponents)):
     plt.plot(V0, Jcomp[j], 'x-', label=Jcomponents[j])
plt.xlabel("voltage bias [mV]")
plt.ylabel("current [pA]")
plt.legend(loc='lower right')
plt.grid(True)

plt.savefig(fname+"add.eps", bbox_inches='tight')

plt.show()
