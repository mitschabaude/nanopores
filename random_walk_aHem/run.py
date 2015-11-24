import random,math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
import sys

def radius(x,y):
    return sqrt(x**2+y**2)

def F(x,y,z):
    return [0,0,-5e-3]

geo = geo_from_xml("aHem")
chi_fluid = geo.indicator("fluid_bulk_top",callable=True)

# "flui0_bulk_top"
# "fluid_bulk_bottom"
# "porebottom"
# "porecenter"
# "poretop"
# "pore"



kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
damp = 1 # diffusion should be much lower in pore than in bulk

D=(kb*T)/(6*pi*visc)*damp*1e9  #diffusion[m^2/s]


gamma=(6*pi*visc) #friction [kg/s]

tau=0.1 #-1 # [ns]
steps=1e3 # 


C=1/(gamma)*tau
coeff=math.sqrt(2*D*1e9*tau)

X=np.zeros((steps))
Y=np.zeros((steps))
Z=np.zeros((steps))

Z[0] = 10
X[0] = 0

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
    if chi_fluid(np.array([float(X[i+1]),float(Y[i+1]),float(Z[i+1])]))==0:
    	exit_x = X[i+1]
    	exit_y = Y[i+1]
    	exit_z = Z[i+1]
    	i+=1
    	break
    i+=1
# print 'end'


# X=X[:i]
# Y=Y[:i]
# Z=Z[:i]

# np.save('X',X)
# np.save('Y',Y)
# np.save('Z',Z)
if exit_x != None:
	file=open("exit_x.txt","a")
	file.write("%f\n" %exit_x)
	file.close()
if exit_y != None:
	file=open("exit_y.txt","a")
	file.write("%f\n" %exit_y)
	file.close()
if exit_z != None:
	file=open("exit_z.txt","a")
	file.write("%f\n" %exit_z)
	file.close()

# print 'number of steps:',i
# print 'time [microsec]=',time*1e-3
# import plot
