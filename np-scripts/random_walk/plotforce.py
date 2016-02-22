from math import cos, sin, pi
import numpy as np
from numpy import linalg as LA
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
import sys
from calculateforce import *
from aHem_array import *
import matplotlib as mpl
import matplotlib.pyplot as plt
F = calculateforce(clscale=6., tol=5e-3) # 6. 5e-3
#def F(x):
#    return [0.,0.,-1e-12]
def radius(x,y):
    return sqrt(x**2+y**2)
def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

z=-4.9
rad=np.linspace(5,500,200)
size=rad.shape[0]
f=np.zeros(size)
fz=np.zeros(size)
temp=np.zeros(8)
tempz=np.zeros(8)
for index in range(size):
    for i in range(8):
        temp[i]=LA.norm(np.array(F(argument(rad[index]*cos(float(i)*(pi*2)/8.),rad[index]*sin(float(i)*(pi*2)/8.),z))))
        tempz[i]=abs(F(argument(rad[index]*cos(float(i)*(pi*2)/8.),rad[index]*sin(float(i)*(pi*2)/8.),z))[2])
    f[index]=np.mean(temp)
    fz[index]=np.mean(tempz)
plt.plot(rad,f)
plt.plot(rad,fz)
plt.show()
