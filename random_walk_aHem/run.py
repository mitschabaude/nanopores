from random import gauss
from math import sqrt, pi
import math
import numpy as np
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
# CALCULATE FORCE:


geo_params = dict(
    l3 = 60.,
    l4 = 10.,
    R = 60.,
    x0 = [5., 0., 10.], # |x0| > 2.2
    exit_i = 1,
)
phys_params = dict(
    bV = .5,
    ahemqs = 0.01,
    rTarget = 0.5*nm,
    bulkcon = 1000.,
)

badexit = {"upperbulkb"}
goodexit = {"exit"}
skip_stokes = True


t = Timer("meshing")
meshdict = generate_mesh(5., "aHem", **geo_params)

print "Mesh generation time:",t.stop()

t = Timer("reading geometry")
geo = geo_from_xml("aHem")

print "Geo generation time:",t.stop()

phys = Physics("pore_molecule", geo, **phys_params)

x0 = geo.params["x0"]
r0 = math.sqrt(sum(x**2 for x in x0))
rnear = r0 - geo.params["rMolecule"]
rfar = r0 + geo.params["rMolecule"]
xnear = map(lambda x: rnear/r0*x, x0)
xfar = map(lambda x: rfar/r0*x, x0)

def avg(u, meas):
    return assemble(u*meas)/assemble(Constant(1.0)*meas)

def exit_times(tau):
    Tmin = tau(xnear)
    Tmax = tau(xfar)
    Tavg = avg(tau, geo.dS("moleculeb"))
    return (Tmin, Tavg, Tmax)
    
pde = PNPS(geo, phys)
if skip_stokes:
    pde.solvers.pop("Stokes")
pde.solve()

(v, cp, cm, u, p) = pde.solutions(deepcopy=True)
F = phys.Feff(v, u)

for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average F in %s:"%domain, assemble(F[2]*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain))

VV = VectorFunctionSpace(geo.mesh, "CG", 1)
F = project(F, VV)



def radius(x,y):
    return sqrt(x**2+y**2)

# def F(x,y,z):
#     return [0,0,-5e-3]

def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

geo = geo_from_xml("aHem")
indicator_ahem = geo.indicator("ahem",callable=True)
indicator_molecule = geo.indicator("molecule",callable=True)
indicator_poretop = geo.indicator("poretop",callable=True)
indicator_membrane_geo = geo.indicator("membrane",callable=True)
def indicator_membrane(vec): #"infinite" large membrane
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>=60.:
        return indicator_membrane_geo(argument(10,0,z))
    else:
        return indicator_membrane_geo(vec)


kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
damp = 1 # diffusion should be much lower in pore than in bulk
D=(kb*T)/(6*pi*visc)*damp*1e9  #diffusion[m^2/s]
gamma=(6*pi*visc) #friction [kg/s]
tau=0.1 #-1 # [ns]
steps=1e3 # 
C=1/(gamma)*tau
coeff=sqrt(2*D*1e9*tau)


counter = np.array([0,0,0,0,0]) # ahem, molecule, poretop, membrane, bulk
EXIT_X, EXIT_Y, EXIT_Z = np.array([]), np.array([]), np.array([])

for index in range(1):
    print index
    X=np.zeros((steps))
    Y=np.zeros((steps))
    Z=np.zeros((steps))

    Z[0] = 10
    i=0
    while i<=steps-2:
        xi_x=gauss(0,1)
        xi_y=gauss(0,1)
        xi_z=gauss(0,1)
        
        X[i+1]=X[i] + coeff*xi_x + C*F(argument(X[i],Y[i],Z[i]))[0]
        Y[i+1]=Y[i] + coeff*xi_y + C*F(argument(X[i],Y[i],Z[i]))[1]
        Z[i+1]=Z[i] + coeff*xi_z + C*F(argument(X[i],Y[i],Z[i]))[2]
        if indicator_membrane(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[3] += 1
            i+=1
            break
        elif indicator_ahem(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[0] += 1
            i+=1
            break
        elif indicator_poretop(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[2] += 1
            i+=1
            break
        elif indicator_molecule(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[1] += 1
            i+=1
            break
        i+=1
    if i>steps-2:
        counter[4] += 1
    else:
        EXIT_X = np.append(EXIT_X, np.array([exit_x]))
        EXIT_Y = np.append(EXIT_Y, np.array([exit_y]))
        EXIT_Z = np.append(EXIT_Z, np.array([exit_z]))

np.save('exit_x',EXIT_X)
np.save('exit_y',EXIT_Y)
np.save('exit_z',EXIT_Z)

# X=X[:i]
# Y=Y[:i]
# Z=Z[:i]

# np.save('X',X)
# np.save('Y',Y)
# np.save('Z',Z)

import plot
