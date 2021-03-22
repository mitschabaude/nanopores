from nanopores import *
from dolfin import *

geo_name = "H_geo"
nm = 1e-9
z0 = -10.*nm

params = dict(
#x0 = None,
#x0 = [0.,0.,z0],
Ry = 15*nm,
Rx = 15*nm,
rMolecule = 0.4*nm,
r0 = 1*nm,
)

phys_params = dict(
Membraneqs = -0.03,
bV = -0.0,
Qmol = -3.*qq,
bulkcon = 3e2,
#uppermbias = 1.,
lowermbias = -.03,
)

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
    
Ry_ = 12*nm
Z = []
F = []
I = []
ndigits = 9+4
PNPSAxisym.tolnewton = 1e0
IllposedNonlinearSolver.newtondamp = 1.
step = .5*nm    

for z in drange(Ry_, -Ry_, -step):
    z = round(z, ndigits)
    print("\nz = ",z*1e9,"[nm]\n")
    Z.append(z*1e9)
    x0 = [0, 0, z]
    meshgen_dict = generate_mesh(1.0, geo_name, x0 = x0, **params)
    geo = geo_from_name(geo_name, x0 = x0, **params)
    phys = Physics("pore_molecule", geo, **phys_params)
    
    pb = LinearPBAxisym(geo, **phys_params)
    pb.maxcells = 12e4
    pb.marking_fraction = 0.8
    pb.solve(refinement=True)
    geo = pb.geo
    v0 = pb.solutions()[0]
    pnps = PNPSAxisym(geo, v0=v0, **phys_params)
    pnps.solve(refinement=False)

    fs = pnps.get_functionals()
    Fdrag = fs["Fp1"] + fs["Fshear1"]
    Fel = fs["Fbare1"]
    F.append(Fdrag + Fel)

    I0 = fs["Javgbtm"] if z>0 else fs["Javgtop"]
    I.append(I0)
    
    print(F)
    
from numpy import isnan

for i,x in enumerate(I):
    if isnan(x):
        F[i] = F[i-1]
        I[i] = I[i-1]

U = []
for f in F:
    if 0 < len(U):
        U.append(U[-1] + 0.5e-12*(f+f0)*step/(kB*T))
    else:
        U.append(0.5e-12*f*step/(kB*T))
    f0 = f
        
for s in ["Z","F","I","U"]:
    save(vars()[s], s + "%.2fnm" %(params["rMolecule"]*1e9,))
    
import matplotlib.pyplot as plt

label = "r = %.2fnm, q = %.1fe" %(params["rMolecule"]*1e9, phys_params["Qmol"]/qq)
fname = "PMF2D%.2fnm.eps" %(params["rMolecule"]*1e9,)

plt.plot(Z,U, label=label)
plt.xlabel("z-coordinate of molecule center [nm]")
plt.ylabel("PMF [kT]")
plt.legend(loc = "lower center")

plt.savefig(fname, bbox_inches='tight')
#plt.show()

