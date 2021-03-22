from nanopores import *
from dolfin import *
import numpy
import json, os
from datetime import datetime
from nanopores.tools.geometry import make_domain,make_boundary
from importlib import import_module

'''
Script to compare simulations with electrodes to simulations
with an applied potential near the pore and to find the 
correct biased voltage for "spatially truncated" simulations
'''

geo_name = "W_2D_geo"
geo_dict = import_vars("nanopores.%s.params_geo" %geo_name)
physical_dict = import_vars("nanopores.physics.params_physical")
default_dict = dict(geo_dict = geo_dict, physical_dict = physical_dict)
print("Default parameters: \n", json.dumps(default_dict, indent=4, sort_keys=True))
nm = geo_dict["nm"]

lsam = geo_dict["lsam"]
#config Fig 1f: sam, dp=26nm, bV=0.2V, c0 = 1M, pH=7.6  -> J~25nA
config1f = [{"x0": None, "lsam": lsam, "r0": 13*nm -lsam,},
            {"bV": 0.2, "bulkcon":1e3, "SAMqs": -0.031, "SiNqs": -0.022}]
#config Fig SI2b: dp=23nm, c0 = 400mM, pH=9
configs2b= [{"x0": None, "lsam": lsam, "r0": 11.5*nm -lsam,},
            {"bulkcon":4e2, "SAMqs": -0.078, "SiNqs": -0.022},]

configl = [geo_dict,]
configl.extend(configs2b)
defparams = dict((k, v) for d in configl for k, v in list(d.items()))
defparams.update({"name":geo_name, "R":250*nm, "bV":0.1,
                  "lowerblayer":True, "membraneblayer": False, })

logdict = dict(params = defparams)
qois = "Javg"  # quantity of interest
qoidict = {}

# varying params
vparams = dict(#Membraneqs = numpy.linspace(0.0, -0.2, 11).tolist(),
               #x0 = [None,] # 0,0,0*nm],],
               #maxcells = numpy.linspace(1e4, 11e4, 3).tolist(),
               #bV = numpy.linspace(0.0, 0.4, 5).tolist(),
               Rz = numpy.linspace(250*nm, 2500*nm, 11).tolist(),
)

def compvmean(rz, v, R, n):
    vsample = []
    for i in range(n):
        vsample.append(v(i*R/n,rz))
    vmean = numpy.mean(vsample)
    vdiff = numpy.max(vsample) - numpy.min(vsample)
    return (vmean,vdiff)
Rz_newbV = []

module = "nanopores.%s.subdomains" %geo_name
subd = import_module(module)
t = Timer('Start')
# for now only one one parameter is varying
for (ps,pran) in list(vparams.items()):
    ps_qoi = []
    aparams = defparams.copy()  # actual params for calculation
    for i in range(len(pran)):
        print(("---\nLoop %s of %s \nparameter: %s = %r"
               %(i+1,len(pran), ps,pran[i])))
        aparams.update({ps: pran[i]})

        pid = "e"
        mesh_dict = generate_mesh(0.5, geo_name, pid=pid, **aparams)
        mesh = Mesh("%s/%s/mesh/mesh%s.xml" %(DATADIR,geo_name,pid))
        # plot(geo.boundaries, title="geometry")

        (subdomains, physical_domain) = make_domain(mesh, subd.subdomain_list(**aparams), True)
        (boundaries, physical_boundary) = make_boundary(mesh, subd.boundaries_list(**aparams), True)
        synonymes = subd.__dict__.get("synonymes")
        synonymesnew = {
            # for testing biased Voltage with electrodes
            "bulk":"rightfluidb",
            "noslip":{"membraneb","moleculeb","upperb","lowerb"},
        }
        synonymes.update(synonymesnew)

        geo = Geometry(None, mesh, subdomains, boundaries, physical_domain, physical_boundary, synonymes, aparams)

        PNPSAxisym.imax = 50
        if "maxcells" in aparams:
            PNPSAxisym.maxcells = aparams["maxcells"]
        else:
            PNPSAxisym.maxcells = 8e4
        PNPSAxisym.marking_fraction = 0.5
        PNPSAxisym.tolnewton = 1e-2
        PNPProblem.method["iterative"] = True

        pnps = PNPSAxisym(geo, **aparams)
        pnps.solvers.pop("Stokes")
        pnps.solve(refinement=True, save_mesh=False, visualize=False)
        pnps.print_results()

        N = geo.mesh.num_cells()
        hmin = geo.mesh.hmin()
        f_dict = pnps.get_functionals()
        if aparams["x0"] is not None:
            if aparams["x0"][2] <=0:
                Javg = f_dict["Javgtop"]
            else:
                Javg = f_dict["Javgbtm"]
        else:
            Javg = 1/2.0*(f_dict["Javgtop"] + f_dict["Javgbtm"])
        f_dict.update(Javg = Javg)
        if aparams["x0"] is not None:
            ps_qoi.append([N, hmin, pran[i], f_dict[qois], f_dict["Feff"],])
        else:
            ps_qoi.append([N, hmin, pran[i], f_dict[qois]])

        debye_length = 1.0/sqrt(2*physical_dict["cFarad"]*aparams["bulkcon"]/\
                            (physical_dict["rpermw"]*physical_dict["eperm"]*physical_dict["UT"]))
        if ps == "bulkcon":
            debye_l = [pran[i], debye_length]
        else:
            debye_l = [aparams["bulkcon"], debye_length]
            
        (v,cp,cm,u,p) = pnps.solutions()
        Rz = geo_dict["Rz"]
        rzvp = Rz
        rzvm = -Rz
        
        (vmeanp,vdiffp) = compvmean(rzvp, v, aparams["R"], 100)
        print(("vmean @ %s" %rzvp), vmeanp, "    vdiff: ", vdiffp)
        (vmeanm,vdiffm) = compvmean(rzvm, v, aparams["R"], 100)
        print(("vmean @ %s" %rzvm), vmeanm, "    vdiff: ", vdiffm)
        newbV = vmeanm - vmeanp
        print("new biased Voltage @ %s :" %rzvm, newbV)
        Rz_newbV.append([N, hmin, aparams["Rz"], newbV, rzvp])
            
    qoidict.update({"N_hmin_%s_%s" %(ps,qois) : ps_qoi})
    qoidict.update({"N_hmin_Rz_bV_rzv": Rz_newbV})
    logdict.update(qoidict = qoidict, bulkcon_debye = debye_l)

t.stop()

logdict.update(dict(
    vparams = vparams,
    solvers = list(pnps.solvers.keys()),
    geo_name = geo_name,
    comptime = t.value(),
    datetime = str(datetime.now()),
))
print("total time: ", t.value())


'''
solve corresponding poisson problem with biased voltage nearer to pore
'''
params1 = defparams.copy()
aparams1 = {"name":geo_name, "membraneblayer":False,
            "Rz": Rz, "bV": newbV,
}
params1.update(aparams1)
#print json.dumps(params1, indent=2, sort_keys=True)
print("Rz: ", params1["Rz"])

pid1 = "e1"
mesh_dict1 = generate_mesh(0.5, geo_name, pid=pid1, **params1)
mesh1 = Mesh("%s/%s/mesh/mesh%s.xml" %(DATADIR,geo_name,pid1))
(subdomains1, physical_domain1) = make_domain(mesh1, subd.subdomain_list(**params1), True)
(boundaries1, physical_boundary1) = make_boundary(mesh1, subd.boundaries_list(**params1), True)
synonymesnew = {
    # for testing biased Voltage without electrodes
    "bulk":{"rightfluidb","upperb","lowerb"},
    "noslip":{"membraneb","moleculeb",},
}
synonymes.update(synonymesnew)

geo1 = Geometry(None, mesh1, subdomains1, boundaries1, physical_domain1, physical_boundary1, synonymes, params1)

pnps1 = PNPSAxisym(geo1, **params1)
pnps1.solvers.pop("Stokes")
pnps1.solve(refinement=True, save_mesh=False, visualize=False)

#compare
v1 = pnps1.solutions()[0]
von1 = interpolate(v, FunctionSpace(geo1.mesh,"CG",1))
plot(von1, title="potential w electrodes")
plot(v1, title="potential1 w bV %s" %newbV)
f_vmeanp = Constant(vmeanp)
diff = project(v1+f_vmeanp-von1, FunctionSpace(geo1.mesh, "CG",1))
print("L2-norm of potential difference: ", norm(diff,"l2"))
print("H1-norm of potential difference: ", norm(diff,"h1"))
print("linf-norm of potential difference: ", norm(diff.vector(),"linf"))
print("new biased Voltage @ %s :" %rzvm, newbV)
plot(diff, title="potential difference", interactive=True)

logdict.update(dict(rzvm=rzvm, newbV=newbV))
print(json.dumps(logdict, indent=4, sort_keys=True))

logdict.update(_default = default_dict, _method=PNPProblem.method)
logdir = os.path.join(DATADIR, geo_name, "log","")
timeidstr = datetime.now().strftime("%Y-%m-%d_%Hh%M")
filestr = timeidstr+".json"
if not os.path.exists(logdir):
    os.makedirs(logdir)
with open(logdir + filestr, "w") as fobj:
    json.dump(logdict, fobj, indent=4, sort_keys = True, ensure_ascii=False)
print("\n+++++\n\n", qoidict)
print("filestr = '%s'" %filestr)


pnps.visualize("fluid")
