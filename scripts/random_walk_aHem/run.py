import matplotlib.path as mplPath
from random import gauss
from math import sqrt, pi, acos, cos, sin, exp
import math
import numpy as np
from numpy import linalg as LA
from nanopores import *
from dolfin import *
import sys
from calculateforce import loadforces
F, Fel, Fdrag = loadforces()
hbond=np.load('hbond.npy')
bbPath=mplPath.Path(hbond) #bigger aHem pore for H-Bond


def radius(x,y):
    return sqrt(x**2+y**2)
def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])
def FF(array):
    if radius(array[0],array[1])>5. or array[2]>5.:
        return [0.,0.,0.]
    return F(array)


geo = geo_from_xml("aHem")
indicator_poretop_geo = geo.indicator("poretop",callable=True)
def indicator_poretop(vec):
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>5.:
    	return 0
    elif z>0.:
    	return 0
    else:
    	return indicator_poretop_geo(vec)

beta=1000
kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
D=(kb*T)/(6*pi*0.5e-9*visc) #diffusion[m^2/s]
gamma=(6*pi*0.5*visc) #friction [microgramm/s]
tau=0.05 # [ns]
steps=1e8# 5 milliseconds = 1e8*tau
C=1/(gamma)*tau # [s^2/kg]==>multiply force with 1e9 to convert from N to kg*nm/s^2
coeff=sqrt(2*D*1e9*tau) # [nm]

TIME=np.load('timer.npy')
counter=np.load('counter.npy')

#######
timeend=5e6
#######
sims=1000
done=np.sum(counter)
left=sims-done

r_h=np.array([])
z_h=np.array([])

#clock=np.zeros(5)
Range = range(done,sims)
start=time()
for index in Range:
    print str(index)+" out of "+str(sims)
    X=np.zeros(steps)
    Y=np.zeros(steps)
    Z=np.zeros(steps)
    Z[0] = 1.0
    timer = 0.

    i=0
    boolexit=False
    newhbond=0
    while timer<timeend and Z[i]<5e2 and X[i]**2+Y[i]**2<2.5e5:
        faraway=False
        timefac=1.
        timefacsq = 1.
        timeadd = tau
        rad=radius(X[i],Y[i])
        xi_x=gauss(0,1)
        xi_y=gauss(0,1)
        xi_z=gauss(0,1)
        [Fx, Fy, Fz] = FF(argument(X[i],Y[i],Z[i]))
        if Z[i]>-5.41+12. and ( rad>10. or Z[i]>7.):
            timefac = 20.
            timefacsq = 4.47213
            timeadd = 1.
            faraway=True
        if Z[i]>100.:
            timefac = 200.
            timefacsq = 14.142136
            timeadd = 10.
        if Z[i]>1000.:
            timefac = 2000.
            timefacsq = 44.72135
            timeadd = 100.
        if Z[i]>3000.:
            timefac = 20000.
            timefacsq = 141.42135
            timeadd = 1000.
        X[i+1] = X[i] + coeff*xi_x*timefacsq + timefac*1e9*C*(Fx)# + fsurfx + fmemx)
        Y[i+1] = Y[i] + coeff*xi_y*timefacsq + timefac*1e9*C*(Fy)# + fsurfy + fmemy)
        Z[i+1] = Z[i] + coeff*xi_z*timefacsq + timefac*1e9*C*(Fz)# + fsurfz + fmemz)
        if bbPath.contains_point((radius(X[i+1],Y[i+1]),Z[i+1])) or Z[i+1]<-5.41+.5:
            i-=1
            newhbond+=1
            timer+=np.random.exponential(beta,size=1)[0]
            r_h=np.append(r_h,np.array([rad]))
            z_h=np.append(z_h,np.array([Z[i+1]]))
        timer += timeadd
        i+=1
        if indicator_poretop(argument(X[i],Y[i],Z[i]))==1:# and Z[i]<-1.0:
            counter[0] += 1
            TIME = np.append(TIME,np.array([timer]))
            boolexit=True
            break
    if not boolexit:
        counter[1] += 1
    np.save('timer',TIME)
    np.save('counter',counter)
end=time()
work=end-start
workh=work/3600.
print 'work = %.1f hours'%workh
print newhbond
from plot_probability import *
