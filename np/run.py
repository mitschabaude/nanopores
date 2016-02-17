from random import gauss
from math import sqrt, pi, acos, cos, sin, exp
import math
import numpy as np
from numpy import linalg as LA
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
import sys
from calculateforce import *
from aHem_array import *
#F = calculateforce(clscale=10., tol=5e-1) # 6. 5e-3
def F(vec):
    return [0.,0.,-1e-12]
def radius(x,y):
    return sqrt(x**2+y**2)
def normal(ax,ay,bx,by,px,py):
    AP2=(ax-px)**2+(ay-py)**2
    BP2=(bx-px)**2+(by-py)**2
    AB2=(ax-bx)**2+(ay-by)**2
    AB=sqrt(AB2)
    c = (AP2-BP2+AB2)/(2*AB)
    qx = ax+c*(bx-ax)/AB
    qy = ay+c*(by-ay)/AB
    if c>0. and c<AB:
        if AP2<=c**2:
            return [0.,qx,qy]
        else:
            return [sqrt(AP2-c**2),qx,qy]
    else:
        return [100.,qx,qy]

size=X_aHem.shape[0]
def dist(rad,z):
    if z>3.:
        return [100.,0.,0.]
    elif rad>8.:
        return [100.,0.,0.]
    else:
        D = np.zeros(size)
        E = np.zeros(size*3).reshape(3,size)
        for index in range(size):
            D[index] = radius(rad-X_aHem[index][0],z-X_aHem[index][2])
            vec = normal(X_aHem[(index-1)][0],X_aHem[(index-1)][2],X_aHem[(index)][0],X_aHem[(index)][2],rad,z)
            E[0][index] = vec[0]
            E[1][index] = vec[1]
            E[2][index] = vec[2]
        di = np.argmin(D)
        ei = np.argmin(E[0])
        if D[di]<E[0][ei]:
            return [D[di],X_aHem[di][0],X_aHem[di][2]]
        else:
            return [E[0][ei],E[1][ei],E[2][ei]]
def surfdiff(r):
    if r>=3.:
        return 0.
    else:
        return 12*1e-14/(r**13)

def F_surf(x,y,z):
    rad=radius(x,y)
    [d,xs,ys] = dist(rad,z)
    if d>3.:
        return np.array([0.,0.,0.,d])
    if rad==0.:
        alpha=0.
    else:
        alpha = acos(x/rad)
    if y<0:
        alpha = 2*pi-alpha
    vec_proj = np.array([rad-xs,0,z-ys])
    vec_proj = vec_proj*(1./LA.norm(vec_proj))
    co=cos(alpha)
    si=sin(alpha)
    A=np.array([[co,-si,0],[si,co,0],[0,0,1]])
    vec = A.dot(vec_proj)
    force = surfdiff(d)
    return np.append(force*vec,np.array([d]))

def F_membrane(x,y,z):
    surf = X_aHem[18][2]
    rad=radius(x,y)
    if rad<4. or z>surf+4.5:
        return np.array([0.,0.,0.,10.])
    r=z-surf
    force = surfdiff(r)
    return np.append(force*np.array([0,0,1]),np.array([r]))

def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

geo = geo_from_xml("aHem")
indicator_ahem_geo = geo.indicator("ahem",callable=True)
indicator_molecule = geo.indicator("molecule",callable=True)
indicator_porecenter_geo = geo.indicator("porecenter",callable=True)
indicator_membrane_geo = geo.indicator("membrane",callable=True)
indicator_poretop_geo = geo.indicator("poretop",callable=True)
def indicator_membrane(vec): #"infinite" large membrane
    x, y, z = vec[0], vec[1], vec[2]
    if z>60.:
    	return 0
    elif radius(x,y)>=30.:
        return indicator_membrane_geo(argument(10,0,z))
    else:
        return indicator_membrane_geo(vec)
def indicator_ahem(vec):
	x, y, z = vec[0], vec[1], vec[2]
	if radius(x,y)>=30. or z>60.:
		return 0
	else:
		return indicator_ahem_geo(vec)
def indicator_porecenter(vec):
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>50.:
    	return 0
    elif z>0.:
    	return 0
    else:
    	return indicator_porecenter_geo(vec)
def indicator_poretop(vec):
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>50.:
    	return 0
    elif z>0.:
    	return 0
    else:
    	return indicator_poretop_geo(vec)

def oor(x,y,z):
    return radius(x,y)>500 or z>500

kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
D=(kb*T)/(6*pi*0.5e-9*visc)  #diffusion[m^2/s]
gamma=(6*pi*0.5*visc) #friction [microgramm/s]
tau=0.05 # [ns]
steps=1e8# 5 milliseconds = 1e8*tau
C=1/(gamma)*tau # [s^2/kg]==>multiply force with 1e9 to convert from N to kg*nm/s^2
coeff=sqrt(2*D*1e9*tau) # [nm]

counter = np.array([0,0,0,0,0]) # ahem, molecule, poretop, membrane, bulk
EXIT_X, EXIT_Y, EXIT_Z = np.array([]), np.array([]), np.array([])
Range = range(1)
for index in Range:
    print str(index)+" out of "+str(len(Range))
    X=np.zeros((steps))
    Y=np.zeros((steps))
    Z=np.zeros((steps))
    Z[0]=10.
    time = 0.
    redos = 0
    hbonds = 0

    i=0
    while time<5e6:
        timefac=1.
        timefacsq = 1.
        timeadd = tau
        b=time/5e4
        rad=radius(X[i],Y[i])
        print
        print ">>>>>>>>>>>>>>>> progress %0.2f percent" %b
        print 'RAD = %.0f , Z = %.1f '%(rad, Z[i])
        sys.stdout.write("\033[F") # Cursor up one line
        sys.stdout.write("\033[F") # Cursor up one line
        sys.stdout.write("\033[F") # Cursor up one line
        xi_x=gauss(0,1)
        xi_y=gauss(0,1)
        xi_z=gauss(0,1)
        if indicator_poretop(argument(X[i],Y[i],Z[i]))==1: #Targetmolecule in Pore=>diffusion damp factor 1/10
            xi_x *= 0.316
            xi_y *= 0.316
            xi_z *= 0.316
        [fsurfx,fsurfy,fsurfz,dsurf] = F_surf(X[i],Y[i],Z[i])
        [fmemx, fmemy, fmemz, dmem] = F_membrane(X[i],Y[i],Z[i])
        [Fx, Fy, Fz] = F(argument(X[i],Y[i],Z[i]))
        if dsurf<1.5 or dmem<1.5:
            hbonds+=1
            time+=1e3*np.random.poisson(10.,1)[0]
            Z[i+1]=Z[i]+5.
            if rad<7. and rad>3.8:
                if rad==0.:
                    alpha=0.
                else:
                    alpha = acos(X[i]/rad)
                if Y[i]<0:
                    alpha = 2*pi-alpha
                co=cos(alpha)
                si=sin(alpha)
                A=np.array([[co,-si],[si,co]])
                vec = A.dot(np.array([5.,0]))
                X[i+1]=X[i]+vec[0]
                Y[i+1]=Y[i]+vec[1]
            else:
                X[i+1]=X[i]
                Y[i+1]=Y[i]
            i+=1
            [fsurfx,fsurfy,fsurfz,dsurf] = F_surf(X[i],Y[i],Z[i])
            [fmemx, fmemy, fmemz, dmem] = F_membrane(X[i],Y[i],Z[i])
            [Fx, Fy, Fz] = F(argument(X[i],Y[i],Z[i]))
        if Z[i]>X_aHem[18][2]+10. and ( rad>15. or Z[i]>10.):
            timefac = 20.
            timefacsq = 4.472135
            timeadd = 1.
        X[i+1]=X[i] + coeff*xi_x*timefacsq + timefac*1e9*C*(Fx + fsurfx + fmemx)
        Y[i+1]=Y[i] + coeff*xi_y*timefacsq + timefac*1e9*C*(Fy + fsurfy + fmemy)
        Z[i+1]=Z[i] + coeff*xi_z*timefacsq + timefac*1e9*C*(Fz + fsurfz + fmemz)
        if LA.norm(F_surf(X[i+1],Y[i+1],Z[i+1])[0:3])+LA.norm(F_membrane(X[i+1],Y[i+1],Z[i+1])[0:3])>5e-10:
            time-=timeadd
            i-=1
            redos+=1
        time += timeadd
#        if indicator_porecenter(argument(X[i+1],Y[i+1],Z[i+1]))==1:
#            exit_x = X[i+1]
#            exit_y = Y[i+1]
#            exit_z = Z[i+1]
#            counter[2] += 1
#            i+=1
#            break
        i+=1

X=X[:i]
Y=Y[:i]
Z=Z[:i]

print 'redos: ',redos
print 'hbonds: ',hbonds
np.save('exit_x',EXIT_X)
np.save('exit_y',EXIT_Y)
np.save('exit_z',EXIT_Z)
np.save('counter',counter)
np.save('x',X)
np.save('y',Y)
np.save('z',Z)
import plot3
#import plot2
