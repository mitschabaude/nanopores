''' record IV curve for Wei pore and compare with experimental data
Experimental IV data of Wei is fitted linearly to pore resistances of 11.5 MOhm without SAM and 16.0 MOhm with SAM
'''

from nanopores import *
from dolfin import *
import numpy
import json
import matplotlib.pyplot as plt

def drange(start, stop, step):
    r = start
    while min(start,stop) <= r <= max(start,stop):
        yield r
        r += step

def save(data, fname):
    with open('%s.txt' % fname, 'w') as f:
        f.write('\n'.join([str(s) for s in data]))

def load(fname):
    with open('%s.txt' % fname, 'r') as f:
        data = f.read().split('\n')
    return data


geo_name = "W_2D_geo"
geo_dict = import_vars("nanopores.geometries.%s.params_geo" %geo_name)
config = "2b"

# no sam and no molecule causes a "found no facets..."
lsam = geo_dict["lsam"]
maxcells = 10e3
r0factor = 1.25

if config=="1f":
    csvfile = "data/weiIV/weiIV1f"
    # config Fig 1f: sam, bV=0.2V, c0 = 1M, pH=7.6
    # dp=26nm was est via resistance -> J~25nA
    geo_params = {"x0": None, "lsam": lsam, "r0": (13*nm -lsam)*r0factor,}
    phys_params = {"bV": 0.2, "bulkcon":1e3, "SAMqs": -0.032, "SiNqs": -0.022}

elif config=="2b":
    csvfile = "data/weiIV/weiIVs2bsam" if lsam else "data/weiIV/weiIVs2b"
    # config Fig SI2b: c0 = 400mM, pH=9
    # dp=23nm was est via resistance (0.2V -> J~
    geo_params = {"x0": None, "lsam": lsam, "r0": (11.5*nm -lsam)*r0factor,}
    phys_params = {"bulkcon":4e2, "SAMqs": -0.0, "SiNqs": -0.022, "bulkconduct": 4.9}
else:
    print("config %s is not a valid configuration" %config)
    exit()

l0 = geo_dict["lsin"]+geo_dict["lau"]+lsam
geo_params.update({"Rx":50*nm, "Rz":200*nm, "l0":l0, "boxfields":True,
                   #"membraneblayer":True, "moleculeblayer":True ,
               })

# IMPORTANT: diffusion and stokes damping in pore
phys_params.update({"stokesdampPore":1.0, "rDPore": 1.0,})

# get experimental data from csv file
IV = numpy.genfromtxt(csvfile+".csv", delimiter=',')
V0 = IV[:,0]
I0 = IV[:,1]*1e3
IV_dict = dict(list(zip(V0,I0)))

PNPSAxisym.tolnewton = 1e-1
PNPSAxisym.maxcells = maxcells
PNPProblemAxisym.method["iterative"] = False
PNPProblemAxisym.method["kparams"]["monitor_convergence"] = False
PNPSAxisym.newtondamp = 1.0

meshgen_dict = generate_mesh(0.8, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore_molecule", geo, **phys_params)

if not geo.parameter("x0"):
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)
if geo.parameter("lsam") < 1e-12:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"samb":set(),"chargedsamb":set()})
    geo.import_synonymes(synonymes)

phys.bV = 0.1
goal = (lambda v : phys.Fbare(v, 1) + phys.CurrentPB(v)) if geo_params["x0"] else (lambda v : phys.CurrentPB(v))
pb = LinearPBAxisymGoalOriented(geo, phys, goal=goal)
pb.maxcells = maxcells
pb.marking_fraction = 0.5
pb.solve(refinement=True)

geo.mesh = pb.geo.mesh
plot(geo.subdomains)
N = geo.mesh.num_cells()
hmin = geo.mesh.hmin()

print(phys.diffusion_factor)
I = []
F = []
G = []
Gnaive = []
Gexp = []
for bVmV in V0:
    bV = bVmV*1e-3
    print("\nbV = %.0f [mV]\n" % bVmV)
    phys.bV = bV if abs(bV) > 1e-7 else None

    phys.bV = bV
    pnps = PNPSAxisym(geo, phys)
    pnps.solve(refinement=False, save_mesh=True, visualize=False)
    #pnps.print_results()
    (v,cp,cm,u,p) = pnps.solutions()

    f_dict = pnps.get_functionals()
    if geo_params["x0"]:
        J = "Javgtop" if geo_params["x0"][2] <=0 else "Javgbtm"
        Javg = f_dict[J]
    else:
        Javg = 1/2.0*(f_dict["Javgtop"] + f_dict["Javgbtm"])
    f_dict.update(Javg = Javg)

    data = [N, hmin, bV, Javg]

    Jcomponents = ["Jzdiff", "Jzdrift", "Jzstokes"]
    Jzdiff = phys.cFarad*phys.D*phys.rDPore*phys.grad(-cp+cm)[1] * phys.r2pi/geo.params["l0"]*geo.dx("pore")
    Jzdrift = phys.cFarad*phys.mu*phys.rDPore*(-cp-cm)*phys.grad(v)[1] * phys.r2pi/geo.params["l0"]*geo.dx("pore")
    Jzstokes = phys.cFarad*phys.stokesdampPore*(cp-cm)*u[1] * phys.r2pi/geo.params["l0"]*geo.dx("pore")
    for j in Jcomponents:
        f_dict.update({j: 1e12*assemble(locals()[j])})

    # naive pore conductance in pS via geometric mean approximation
    AgmeanPore = pi*geo_params["r0"]*(geo_params["r0"] + geo_params["l0"]/numpy.tan((90-geo_dict["angle"]/2)*pi/180))
    Gnaive_ = 1e12*phys.bulkconduct/geo_params["l0"]*AgmeanPore
    Jnaive = bV*Gnaive_

    Gexp_ = IV_dict[bVmV]/bV
    Gexp.append(Gexp_)

    if geo_params["x0"]:
        Feff = f_dict["Fp_z"] + f_dict["Fshear_z"] + f_dict["Fbare_z"]
        f_dict.update(Feff = Feff)
    else: pass

    F.append(f_dict)
    I_ = f_dict["Javg"]
    I.append(I_)

    V = (v([0.0, -l0/2*1.5]) - v([0.0, l0/2*1.5]))
    print("l0: ", l0*1e9, "[nm] ", "lsam: ", lsam*1e9, "[nm]", "r0factor: ", r0factor)
    print("naive current :", Jnaive, "[pA]")
    print("current in experiment: ", IV_dict[bVmV], "[pA]")
    print("I (current through pore center):", I_, "[pA]")
    print("V (transmembrane potential):", V, "[V]")
    print("bV: ", phys.bV, "[V]")
    print("simulated resistance bV/I:", bV/I_*1e6, "[MOhm]")
    print("experimental resistance :", 1./Gexp_ *1e6, "[MOhm]")
    print("naive resistance:", 1./Gnaive_*1e6, "[MOhm]")
    print("simulated conductance I/bV:", I_/bV*1e-3, "[nS]")
    print("experimental conductance:", Gexp_ *1e-3, "[nS]")
    print("naive conductance:", Gnaive_*1e-3, "[nS]")
    Gnaive.append(Gnaive_)
    G.append(I_/V)

print(json.dumps(F, indent=4, sort_keys=True))

plot(v)
plot(cm-cp, title="cdiff")
plot(cm+cp)
plot(u)
interactive()

plt.figure(1)
plt.plot(V0,I0, 's-', label="Wei et. al")
plt.plot(V0,I, 'o-', label="PNPS simulation")
plt.xlabel("voltage bias [mV]")
plt.ylabel("current [pA]")
plt.legend(loc='lower right')

fname = ("%s_N%s_lsam%.1f_r0factor%.2f" %(csvfile, N, lsam*1e9, r0factor) )
plt.savefig(fname+".eps") # bbox_inches='tight')

Jcomp = numpy.asarray([[f[j] for f in F] for j in Jcomponents])
plt.figure(2)
plt.subplot(211)
plt.plot(V0,Gnaive, '--', label="naive conductance [pS]")
plt.plot(V0,Gexp, '--', label="experimental conductance [pS]")
plt.plot(V0,G, 'x-', label="simulated conductance [pS]")
plt.xlabel("voltage bias [mV]")
plt.ylabel("conductance [pS]")
plt.legend(loc='lower right')
plt.subplot(212)
for j in range(len(Jcomponents)):
     plt.plot(V0, Jcomp[j], 'x-', label=Jcomponents[j])
plt.xlabel("voltage bias [mV]")
plt.ylabel("current [pA]")
plt.legend(loc='lower right')
plt.grid(True)

plt.savefig(fname+"add.eps") #, bbox_inches='tight')

plt.show()
