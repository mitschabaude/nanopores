#http://nbviewer.ipython.org/github/jrjohansson/scientific-python-lectures/blob/master/Lecture-4-Matplotlib.ipynb

import random,math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab

import sys

from nanopores import *
from nanopores.H_cyl_geo.params_geo import *
nm=1e-9
from collision import *

from nanopores.force import F

################
kb=1.3806488e-23 #boltzmann [J/K]
T=293 #temp [K]
Md=2.0 #diameter particle [nm]
rnm=Md/2.0
r=rnm*1e-9
visc=1e-3 #[Pa*s]
D=(kb*T)/(6*pi*r*visc) #diffusion[m^2/s]
#F=-1e-12 #force[N]
m=1e-21 #mass[kg]

gamma=(6*pi*r*visc)/m #friction [1/s]



tau=1e-1 # [ns]
steps=1e7 # 


C=F/(m*gamma)*tau
coeff=math.sqrt(2*D*1e9*tau)
pi=math.pi

X=np.zeros((steps))
Y=np.zeros((steps))
Z=np.zeros((steps))
#X_C=np.zeros((steps))
#Y_C=np.zeros((steps))
#Z_C=np.zeros((steps))

#Z[0]=l0*5e8


time=0
i=0
j=0
S=np.zeros(6)

while i<=steps-2 and Z[i]>=-l0*5e8 and Z[i] <=15 and radius(X[i],Y[i])<=7.5:
    time+=1
    xi_x=random.gauss(0,1)
    xi_y=random.gauss(0,1)
    xi_z=random.gauss(0,1)
    X[i+1]=X[i]+coeff*xi_x
    Y[i+1]=Y[i]+coeff*xi_y
    Z[i+1]=Z[i]+coeff*xi_z+C

    while True:
        ax, ay, az, bx, by, bz=X[i], Y[i], Z[i], X[i+1], Y[i+1], Z[i+1]
        S[0]=s_innercylinder(ax,ay,az,bx,by,bz)
        S[1]=s_dnatop(ax,ay,az,bx,by,bz)
        S[2]=s_innertorus(ax,ay,az,bx,by,bz)
        S[3]=s_outertorus(ax,ay,az,bx,by,bz)
        S[4]=s_outercylinder(ax,ay,az,bx,by,bz)
        S[5]=s_membran(ax,ay,az,bx,by,bz)
        s=np.min(S)
        if s>1.0 or s==0.0:
            break
        pos=np.argmin(S)
        if pos==0:
            newpoint=newpoint_innercylinder
            #s_function=s_innercylinder
        elif pos==1:
            newpoint=newpoint_dnatop
            #s_function=s_dnatop
        elif pos==2:
            newpoint=newpoint_innertorus
            #s_function=s_innertorus
        elif pos==3:
            newpoint=newpoint_outertorus
            #s_function=s_outertorus
        elif pos==4:
            newpoint=newpoint_outercylinder
            #s_function=s_outercylinder_membran  
        elif pos==5:
            newpoint=newpoint_membran
            #s_function      
    
        j+=1
        if True:
            print 's=', s
            print 'pos=',pos
        print j
        X[i+1]=ax+s*(bx-ax)
        Y[i+1]=ay+s*(by-ay)
        Z[i+1]=az+s*(bz-az)
        i+=1
        array=newpoint(ax,ay,az,bx,by,bz,s)
        X[i+1], Y[i+1], Z[i+1]=array[0], array[1], array[2]

    i+=1
if Z[i]<-l0*5e8:
    print 'unten'





###################
#X_C=X_C[:j]
#Y_C=Y_C[:j]
#Z_C=Z_C[:j]

X=X[:i]
Y=Y[:i]
Z=Z[:i]

np.save('X',X)
np.save('Y',Y)
np.save('Z',Z)
#np.save('X_C',X_C)
#np.save('Y_C',Y_C)
#np.save('Z_C',Z_C)
x1=np.zeros(100)+r1*1e9+rMolecule*1e9
y1=np.linspace(l1*5e8+rMolecule*1e9,l0*5e8,100)
x2=np.linspace(r1*1e9+rMolecule*1e9,7.5,100)
y2=np.zeros(100)+l1*5e8+rMolecule*1e9
x3=np.linspace(r0*1e9,r1*1e9,100)
y3=np.zeros(100)+l0*5e8+rMolecule*1e9
x4=np.zeros(100)+r0*1e9-rMolecule*1e9
y4=np.linspace(-l0*5e8,l0*5e8,100)
x=np.linspace(0,pi/2,100)
y=np.linspace(0,pi/2,100)
x_1=np.cos(x)*rMolecule*1e9+r1*1e9
y_1=np.sin(y)*rMolecule*1e9+l0*5e8
x_=np.linspace(pi/2,pi,100)
y_=np.linspace(pi/2,pi,100)
x_2=np.cos(x_)*rMolecule*1e9+r0*1e9
y_2=np.sin(y_)*rMolecule*1e9+l0*5e8
plt.plot(x_2,y_2)
plt.plot(-x_2,y_2)
plt.plot(x_1,y_1)
plt.plot(-x_1,y_1)
plt.plot(-x1,y1)
plt.plot(-x2,y2)
plt.plot(-x3,y3)
plt.plot(-x4,y4)
plt.plot(x4,y4)
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x3,y3)
X_norm=np.zeros(X.shape[0])
for index in range(X.shape[0]):
    X_norm[index]=radius(X[index],Y[index])

plt.plot(X,Z)
plt.scatter(X,Z)
#plt.plot(np.linspace(0,time*1e-3,Z.shape[0]),Z)
#plt.plot(np.linspace(0,time*1e-3,Z.shape[0]),X_norm)
#plt.plot(np.linspace(0,time*1e-3,Z.shape[0]),np.full(Z.shape[0],r0*1e9-rMolecule*1e9))
#plt.savefig('plot_z.png')
print 'collisions: ',j
print 'time [microsec]=',time*1e-3
import plot
