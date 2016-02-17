from nanopores import *
from dolfin import *

geo_name = "H_geo"
nm = 1e-9
params = dict(
Ry = 30*nm,
Rx = 15*nm,
rMolecule = 0.77*nm,
r0 = 1.2*nm,
)
phys_params = {"Membraneqs": -0.03, "bV": -0.1, "Qmol": -0*qq,}

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
    
Ry_ = params["Ry"]-2*params["rMolecule"]
Z = []
V = []
I = []
ndigits = 9+4

for z in drange(Ry_, -Ry_, -1*nm):
    z = round(z, ndigits)
    print "\nz = ",z,"\n"
    Z.append(z*1e9)
    x0 = [0, 0, z]
    meshgen_dict = generate_mesh(0.5, geo_name, x0 = x0, **params)
    geo = geo_from_name(geo_name, x0 = x0, **params)
    PNPSAxisym.tolnewton = 1e-1
    pnps = PNPSAxisym(geo, **phys_params)
    pnps.solve(refinement=False, save_mesh=False)

    (v,cp,cm,u,p) = pnps.solutions()
    I0 = -pnps.get_functionals()["Javgbtm"] if z>0 else -pnps.get_functionals()["Javgtop"]
    I.append(I0)
    V.append(v([0.0, 10*nm]) - v([0.0, -10*nm]))
    #print "I (current through pore center):",I,"[pA]"
    #print "V (transmembrane potential):",V,"[V]"
    #print "conductance I/V:",I/V,"[pS]"

from numpy import isnan

for i,x in enumerate(V):
    if isnan(x):
        V[i] = V[i-1]
        I[i] = I[i-1]
        
#for s in ["Z","I","V"]:
#    save(vars()[s], s)
    
import matplotlib.pyplot as plt

label = "r = %.2fnm" %(params["rMolecule"]*1e9,)
fname = "current_%dmV_%de_%.2fnm.eps" %(int(phys_params["bV"]*1000), int(phys_params["Qmol"]/qq), params["rMolecule"]*1e9)

plt.plot(Z,I, label=label)
plt.xlabel("z-coordinate of molecule center [nm]")
plt.ylabel("current [pA]")
plt.legend(loc='lower right')

plt.savefig(fname, bbox_inches='tight')

