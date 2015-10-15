import random, math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab

import sys

# bad hack to allow arbitrary nm in params_geo
from nanopores import import_vars
params = import_vars("nanopores.geometries.H_cyl_geo.params_geo")
for x in params:
    if params[x] is not None:
        exec("%s = %s*1e0/%s" %(x, params[x], params['nm']))


from collision import *
from nanopores.force import F, Jz, fmaps
#F = fmaps["Fp"]
#import dolfin
#dolfin.plot(F)

nm_force = 1e-9

################
rMolecule = 0.55

################
kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
#rnm = 0.5*1.1 #radius particle [nm]
r = rMolecule*1e-9 #rnm*1e-9
visc=1e-3 #[Pa*s]
damp = 1e-0 # diffusion should be much lower in pore than in bulk

D=(kb*T)/(6*pi*r*visc)*damp  #diffusion[m^2/s]
#F=-1e-12 #force[N]
m= 1e-21 #mass[kg]

gamma=(6*pi*r*visc)/damp #divided by mass - already reduced #friction [1/s]

tau=0.1 # [ns]
steps=1e6


C = 1e9*1e-12/(gamma)*(tau*1e-9) #gamma times mass - already reduced
coeff = 1e9*math.sqrt(2*D*(tau*1e-9))
pi=math.pi

X=np.zeros((steps))
Y=np.zeros((steps))
Z=np.zeros((steps))
J = np.zeros((steps))
TT = np.zeros((steps))
#X_C=np.zeros((steps))
#Y_C=np.zeros((steps))
#Z_C=np.zeros((steps))

Z[0] = l0/2+1.
X[0] = 0

time=0
i=0
j=0
S=np.zeros(11)

def cyl2cart(F, x):
    # transform vector in cylindrical coordinate system to cartesian
    r = math.sqrt(x[0]**2 + x[1]**2)
    if r==0.0:
        return [0.0, 0.0, F[2]]
        
    Fx = 1/r * (x[0]*F[0] - x[1]*F[1])
    Fy = 1/r * (x[1]*F[0] + x[0]*F[1])
    return [Fx, Fy, F[2]]
    
print "\nTypical random displacement |dx|_L2 [nm] ",coeff
print "Displacement caused by force of 1 pN [nm] ",C
print "\n"*4

def lowerradius(x,y,z):
    if z>= -l0/2-rMolecule:
        return 0.
    else:
        return math.sqrt((z+l0/2+rMolecule)**2+x**2+y**2)

while i<=steps-2 and Z[i]>=-l0*0.5-5 and Z[i] <=15 and radius(X[i],Y[i])<=7.5 and lowerradius(X[i],Y[i],Z[i])<=3.:
    print "\x1b[A"*4,"\r",
    timer_a=0
    if association(radius(X[i],Y[i]),Z[i]):
        timer_a+=1
        print 'Association!'
        '''while dissociation()==False:
            i+=1
            X[i]=X[i-1]
            Y[i]=Y[i-1]
            Z[i]=Z[i-1]
            time+=1
            timer_b+=1
            TT[i-1]=time-1
            J[i-1]=J[i-2]
        print 'Dissociation after ',timer_b*1e-3'''
    timer_b=0
    timer_a=0
    
    TT[i] = time
    time+=tau
    xi_x=random.gauss(0,1)
    xi_y=random.gauss(0,1)
    xi_z=random.gauss(0,1)
    Fi = F([sqrt(X[i]**2 + Y[i]**2)*nm_force, Z[i]*nm_force])
    J[i] = Jz([sqrt(X[i]**2 + Y[i]**2)*nm_force, Z[i]*nm_force])
    
    #print "R force:",Fi[0]
    Fi = cyl2cart(Fi, (X[i], Y[i]))
    
    #print [sqrt(X[i]**2 + Y[i]**2)*nm_force, Z[i]*nm_force]
    #print Fi
    
    X[i+1]=X[i] + coeff*xi_x + Fi[0]*C
    Y[i+1]=Y[i] + coeff*xi_y + Fi[1]*C
    Z[i+1]=Z[i] + coeff*xi_z + Fi[2]*C #- 1*C
    
    
    print "Z force:", -Fi[2]*C, "  Z step:", coeff*xi_z - Fi[2]*C 
    print "Z position [nm]:",Z[i]
    print "Z forcing [nm]:",-Fi[2]*C
    print "Time:",time,"[ns]"

    while True:
        ax, ay, az, bx, by, bz=X[i], Y[i], Z[i], X[i+1], Y[i+1], Z[i+1]
        S[0]=s_innercylinder(ax,ay,az,bx,by,bz)
        S[1]=s_dnatop(ax,ay,az,bx,by,bz)
        S[2]=s_innertorus(ax,ay,az,bx,by,bz)
        S[3]=s_outertorus(ax,ay,az,bx,by,bz)
        S[4]=s_outercylinder(ax,ay,az,bx,by,bz)
        S[5]=s_membran(ax,ay,az,bx,by,bz)
        S[6]=s_dnabottom(ax,ay,az,bx,by,bz)
        S[7]=s_innertorus_bottom(ax,ay,az,bx,by,bz)
        S[8]=s_outertorus_bottom(ax,ay,az,bx,by,bz)
        S[9]=s_outercylinder_bottom(ax,ay,az,bx,by,bz)
        S[10]=s_membran_bottom(ax,ay,az,bx,by,bz)
        s=np.min(S)
        if s>1.0 or s==0.0:
            break
        pos=np.argmin(S)
        if pos==0:
            newpoint=newpoint_innercylinder
        elif pos==1:
            newpoint=newpoint_dnatop
        elif pos==2:
            newpoint=newpoint_innertorus
        elif pos==3:
            newpoint=newpoint_outertorus
        elif pos==4 or pos==9:
            newpoint=newpoint_outercylinder
        elif pos==5 or pos==10:
            newpoint=newpoint_membran
        elif pos==6:
            newpoint=newpoint_dnabottom
        elif pos==7:
            newpoint=newpoint_innertorus_bottom
        elif pos==8:
            newpoint=newpoint_outertorus_bottom
        elif pos==10:
            newpoint=newpoint_membran_bottom
    
        j+=1
        if False:
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

file=open("results.txt","a")
if Z[i]<-l0*0.5-rMolecule:
    print 'bottom'
    file.write("{:.3f}".format(time*1e-3)+"\n") 
elif Z[i]>0:
    print 'top'
    file.write("-1.00\n")
elif i==steps-1:
    'ran out of steps'
    file.write("-2.00\n")
file.close()

'''
X=X[:i]
Y=Y[:i]
Z=Z[:i]
J=J[:i]
TT=TT[:i]
#print TT[:], J[:]
#print J[J.nonzero()[0]]
plt.plot(J[J.nonzero()[0]])
plt.show()

np.save('X',X)
np.save('Y',Y)
np.save('Z',Z)
np.save('TT',Y)
np.save('J',Z)

x1=np.zeros(100)+r1+rMolecule
y1=np.linspace(l1*0.5+rMolecule,l0*0.5,100)
x2=np.linspace(r1+rMolecule,7.5,100)
y2=np.zeros(100)+l1*0.5+rMolecule
x3=np.linspace(r0,r1,100)
y3=np.zeros(100)+l0*0.5+rMolecule
x4=np.zeros(100)+r0-rMolecule
y4=np.linspace(-l0*0.5,l0*0.5,100)
x=np.linspace(0,pi/2,100)
y=np.linspace(0,pi/2,100)
x_1=np.cos(x)*rMolecule+r1
y_1=np.sin(y)*rMolecule+l0*0.5
x_=np.linspace(pi/2,pi,100)
y_=np.linspace(pi/2,pi,100)
x_2=np.cos(x_)*rMolecule+r0
y_2=np.sin(y_)*rMolecule+l0*0.5
x__=np.linspace(pi,4*pi/2,100)
y__=np.linspace(pi,4*pi/2,100)
x_3=np.cos(x__)*rMolecule+r0
y_3=np.sin(y__)*rMolecule-l0*0.5
x5=np.linspace(r0,r1,100)
y5=np.zeros(100)-l0*0.5-rMolecule
x6_1=np.linspace(3*pi/2,2*pi,100)
y6_1=x6_1[:]
x6_2=np.cos(x6_1)*rMolecule+r1
y6_2=np.sin(y6_1)*rMolecule-l0/2
x7=np.zeros(100)+r1+rMolecule
y7=np.linspace(-l0/2,-l1/2-rMolecule,100)
x8=np.linspace(r1+rMolecule,7.5,100)
y8=np.zeros(100)-l1/2-rMolecule
plt.plot(x_2,y_2,color='black')
plt.plot(-x_2,y_2,color='black')
plt.plot(x_1,y_1,color='black')
plt.plot(-x_1,y_1,color='black')
plt.plot(-x1,y1,color='black')
plt.plot(-x2,y2,color='black')
plt.plot(-x3,y3,color='black')
plt.plot(-x4,y4,color='black')
plt.plot(x4,y4,color='black')
plt.plot(x1,y1,color='black')
plt.plot(x2,y2,color='black')
plt.plot(x3,y3,color='black')
plt.plot(x_3,y_3,color='black')
plt.plot(-x_3,y_3,color='black')
plt.plot(x5,y5,color='black')
plt.plot(-x5,y5,color='black')
plt.plot(x6_2,y6_2,color='black')
plt.plot(x7,y7,color='black')
plt.plot(x8,y8,color='black')
X_norm=np.zeros(X.shape[0])
for index in range(X.shape[0]):
    X_norm[index]=radius(X[index],Y[index])

plt.plot(X_norm,Z)
plt.scatter(X_norm,Z)

print 'collisions: ',j
print 'number of steps:',i
print 'time [microsec]=',time*1e-3
import plot
'''
