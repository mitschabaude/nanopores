import random,math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab

import sys

# bad hack to allow arbitrary nm in params_geo
from nanopores import import_vars
params = import_vars("nanopores.W_3D_geo.params_geo")
for x in params:
    exec("%s = %s*1e0/%s" %(x, params[x], params['nm']))


#from collision import *
from nanopores.force import F
from additionalforce import *
nm_force = nm

def radius(x,y):
    return sqrt(x**2+y**2)

################
kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
#rnm = 0.5*1.1 #radius particle [nm]
r = rMolecule*nm #rnm*1e-9
visc=1e-3 #[Pa*s]
damp = 0.1 # diffusion should be much lower in pore than in bulk

D=(kb*T)/(6*pi*r*visc)*damp  #diffusion[m^2/s]
if nm==1:
    D=D*1e9
#F=-1e-12 #force[N]
m= 1e-21 #mass[kg]

gamma=(6*pi*r*visc) #friction [kg/s]

tau=0.5 #-1 # [ns]
steps=1e4 # 


C=1e-12/(gamma)*tau
coeff=math.sqrt(2*D*1e9*tau)
pi=math.pi

X=np.zeros((steps))
Y=np.zeros((steps))
Z=np.zeros((steps))

Z[0] = 3
X[0] = 0

time=0
i=0
j=0
S=np.zeros(6)
pitch=1.
delta=0.

def cyl2cart(F, x):
    # transform vector in cylindrical coordinate system to cartesian
    r = radius(x[0],x[1])
    if r==0.0:
        return [0.0, 0.0, F[2]]
        
    Fx = 1/r * (x[0]*F[0] - x[1]*F[1])
    Fy = 1/r * (x[1]*F[0] + x[0]*F[1])
    return [Fx, Fy, F[2]]
    
print "\n"*2

while i<=steps-2 and Z[i] <=R and radius(X[i],Y[i])<=R:
    time+=1
    xi_x=random.gauss(0,1)
    xi_y=random.gauss(0,1)
    xi_z=random.gauss(0,1)
    Fi = F([sqrt(X[i]**2 + Y[i]**2)*nm, Z[i]*nm])
    
    Fi = cyl2cart(Fi, (X[i], Y[i]))
    X[i+1]=X[i] + coeff*xi_x + Fi[0]*C + additionalforce(X[i],Y[i],Z[i],delta=delta,pitch=pitch)[0]
    Y[i+1]=Y[i] + coeff*xi_y + Fi[1]*C + additionalforce(X[i],Y[i],Z[i],delta=delta,pitch=pitch)[1]
    Z[i+1]=Z[i] + coeff*xi_z + Fi[2]*C + additionalforce(X[i],Y[i],Z[i],delta=delta,pitch=pitch)[2]
    
    print "\x1b[A"*2,"\r",
    print "Z position:",Z[i]
    print "Z forcing:",-Fi[2]*C/tau

    i+=1


X=X[:i]
Y=Y[:i]
Z=Z[:i]

np.save('X',X)
np.save('Y',Y)
np.save('Z',Z)

print 'collisions: ',j
print 'number of steps:',i
print 'time [microsec]=',time*1e-3
import plot
