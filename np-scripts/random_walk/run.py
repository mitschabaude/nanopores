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
F = calculateforce(clscale=6., tol=5e-3) # 6. 5e-3
#def F(vec):
#    return [0.,0.,-1e-12]
def radius(x,y):
    return sqrt(x**2+y**2)
def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])
def FF(array):
    if radius(array[0],array[1])>50. or array[2]>50.:
        return [0.,0.,0.]
    return F(array)
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
        return 12*1e-14/((r+.3)**13)

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


geo = geo_from_xml("aHem")
indicator_porecenter_geo = geo.indicator("porecenter",callable=True)
indicator_poretop_geo = geo.indicator("poretop",callable=True)
indicator_aHem_geo = geo.indicator("ahem",callable=True)
def indicator_aHem(vec):
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>5. or z>0. or z<-10.5:
        return 0
    else:
        return indicator_aHem_geo(vec)
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

kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
D=(kb*T)/(6*pi*0.5e-9*visc) #diffusion[m^2/s]
gamma=(6*pi*0.5*visc) #friction [microgramm/s]
tau=0.05 # [ns]
steps=1e8# 5 milliseconds = 1e8*tau
C=1/(gamma)*tau # [s^2/kg]==>multiply force with 1e9 to convert from N to kg*nm/s^2
coeff=sqrt(2*D*1e9*tau) # [nm]

counter = np.array([0,0])
EXIT_X, EXIT_Y, EXIT_Z, TIME = np.array([]), np.array([]), np.array([]), np.array([])
Range = range(1)
for index in Range:
    print str(index)+" out of "+str(len(Range))
    X=np.zeros(steps)
    Y=np.zeros(steps)
    Z=np.zeros(steps)
    xia_x=np.zeros(steps)
    xia_y=np.zeros(steps)
    xia_z=np.zeros(steps)
    Z[0] = 2.
    time = 0.
    redos = 0
    hbonds = 0

    i=0
    timeend=5e6
    mean_hbond = 1e3 #1 microsec
    lambda_poisson = 10.
    boolexit=False
    while time<timeend and Z[i]<1e6 and X[i]**2+Y[i]**2<1e12:
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
            xi_x *= 0.316228
            xi_y *= 0.316228
            xi_z *= 0.316228
        [fsurfx,fsurfy,fsurfz,dsurf] = F_surf(X[i],Y[i],Z[i])
        [fmemx, fmemy, fmemz, dmem] = F_membrane(X[i],Y[i],Z[i])
        [Fx, Fy, Fz] = FF(argument(X[i],Y[i],Z[i]))
        if dsurf<0.5:
            hbonds+=1
            time+=mean_hbond*np.random.poisson(lambda_poisson,1)[0]
            vec = np.array([fsurfx,fsurfy,fsurfz])
            vec *= 1./(LA.norm(vec))
            X[i+1]=X[i]+vec[0]
            Y[i+1]=Y[i]+vec[1]
            Z[i+1]=Z[i]+vec[2]
            i+=1
            rad=radius(X[i],Y[i])
            [fsurfx,fsurfy,fsurfz,dsurf] = F_surf(X[i],Y[i],Z[i])
            [fmemx, fmemy, fmemz, dmem] = F_membrane(X[i],Y[i],Z[i])
            [Fx, Fy, Fz] = FF(argument(X[i],Y[i],Z[i]))
        if dmem<0.5:
            hbonds+=1
            time+=mean_hbond*np.random.poisson(lambda_poisson,1)[0]
            X[i+1]=X[i]
            Y[i+1]=Y[i]
            Z[i+1]=Z[i]+1.
            i+=1
            rad=radius(X[i],Y[i])
            [fsurfx,fsurfy,fsurfz,dsurf] = F_surf(X[i],Y[i],Z[i])
            [fmemx, fmemy, fmemz, dmem] = F_membrane(X[i],Y[i],Z[i])
            [Fx, Fy, Fz] = FF(argument(X[i],Y[i],Z[i]))
        if Z[i]>X_aHem[18][2]+12. and ( rad>10. or Z[i]>7.):
            timefac = 20.
            timefacsq = 4.47213
            timeadd = 1.
        if Z[i]>100. and rad>100.:
            timefac = 200.
            timefacsq = 14.142136
            timeadd = 10.
        X[i+1]=X[i] + coeff*xi_x*timefacsq + timefac*1e9*C*(Fx + fsurfx + fmemx)
        Y[i+1]=Y[i] + coeff*xi_y*timefacsq + timefac*1e9*C*(Fy + fsurfy + fmemy)
        Z[i+1]=Z[i] + coeff*xi_z*timefacsq + timefac*1e9*C*(Fz + fsurfz + fmemz)
        xia_x[i]=xi_x
        xia_y[i]=xi_y
        xia_z[i]=xi_z
        if LA.norm(F_surf(X[i+1],Y[i+1],Z[i+1])[0:3])+LA.norm(F_membrane(X[i+1],Y[i+1],Z[i+1])[0:3])>5e-10 or indicator_aHem(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            time-=timeadd
            i-=1
            redos+=1
        time += timeadd
        i+=1
        if indicator_poretop(argument(X[i],Y[i],Z[i]))==1 and Z[i]<-1.0:
            exit_x = X[i]
            exit_y = Y[i]
            exit_z = Z[i]
            counter[0] += 1
            TIME = np.append(TIME,np.array([time]))
            boolexit=True
            break
        if Z[i]<-5.41 or indicator_aHem(argument(X[i],Y[i],Z[i]))==1:
            np.save('x',X)
            np.save('y',Y)
            np.save('z',Z)
            np.save('xia_x',xia_x)
            np.save('xia_y',xia_y)
            np.save('xia_z',xia_z)
            sys.exit()
    if not boolexit:
        counter[1] += 1
        exit_x = X[i]
        exit_y = Y[i]
        exit_z = Z[i]
    EXIT_X = np.append(EXIT_X,np.array([exit_x]))
    EXIT_Y = np.append(EXIT_Y,np.array([exit_y]))
    EXIT_Z = np.append(EXIT_Z,np.array([exit_z]))

X=X[:i]
Y=Y[:i]
Z=Z[:i]

#print 'redos: ',redos
#print 'hbonds: ',hbonds
np.save('x',X)
np.save('y',Y)
np.save('z',Z)
np.save('exit_x',EXIT_X)
np.save('exit_y',EXIT_Y)
np.save('exit_z',EXIT_Z)
np.save('time',TIME)
np.save('counter',counter)
import plot3
