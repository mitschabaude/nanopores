import random
from math import sqrt, pi
import numpy as np
from nanopores import *
from dolfin import *

def radius(x,y):
    return sqrt(x**2+y**2)

def F(x,y,z):
    return [0,0,-5e-3]

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

counter = np.array([0,0,0,0,0]) # ahem, molecule, poretop, membrane, bulk

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

X=np.zeros((steps))
Y=np.zeros((steps))
Z=np.zeros((steps))

Z[0] = 10
X[0] = 10

time=0
i=0
j=0
S=np.zeros(6)
exit_x = None
exit_y = None
exit_z = None
while i<=steps-2:
    time+=1
    xi_x=random.gauss(0,1)
    xi_y=random.gauss(0,1)
    xi_z=random.gauss(0,1)
    
    X[i+1]=X[i] + coeff*xi_x + C*F(X[i],Y[i],Z[i])[0]
    Y[i+1]=Y[i] + coeff*xi_y + C*F(X[i],Y[i],Z[i])[1]
    Z[i+1]=Z[i] + coeff*xi_z + C*F(X[i],Y[i],Z[i])[2]
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
print counter


X=X[:i]
Y=Y[:i]
Z=Z[:i]

np.save('X',X)
np.save('Y',Y)
np.save('Z',Z)

print 'number of steps:',i
# print 'time [microsec]=',time*1e-3
import plot
