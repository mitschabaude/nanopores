import random,math
from math import sqrt, pi, atan2
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
import sympy as sp

from nanopores import *

# bad hack to allow arbitrary nm in params_geo
params = import_vars("nanopores.geometries.H_cyl_geo.params_geo")
for x in params:
    if params[x] is not None:
        exec("%s = %s*1e0/%s" %(x, params[x], params['nm']))

################
rMolecule = 0.55
################
tau=1e-1 #[ns]
ra=20.9e6 #[1/(Mol*sec)]
Conc=400.e-9 #Mol
rd=0.27 #1/sec
pa=0.#ra*Conc*tau*1e-9
pd=0.#rd*tau*1e-9
def association(rad,z):
    delta=1.0*rMolecule
    if not (z<=l0*0.5+delta and z>=-l0*0.5-delta and rad<=r1 and rad>=r0-delta):
        return False
    if np.random.binomial(1,pa)==0:
        return False
    else: return True
def dissociation():
    if np.random.binomial(1,pd)==0:
        return False
    else: return True

def sign(x):
    if x>0.:
        return 1
    elif x<0.:
        return -1
    else: return 0

def radius(x,y):
    return sqrt(x**2+y**2)

def collision_2(x,y,z):
    return sqrt(pow(x,2)+pow(y,2))>=r0-rMolecule

def para(xa,ya,xb,yb):
    return (-(xa*(xb-xa)+ya*(yb-ya))+sqrt(pow(xa*(xb-xa)+ya*(yb-ya),2)-(pow(xa,2)+pow(ya,2)-pow(r0-rMolecule,2))*(pow(xb-xa,2)+pow(yb-ya,2))))/(pow(xb-xa,2)+pow(yb-ya,2))

def para2(xa,ya,xb,yb):
    return (-(xa*(xb-xa)+ya*(yb-ya))+sqrt(pow(xa*(xb-xa)+ya*(yb-ya),2)-(pow(xa,2)+pow(ya,2)-pow(r1+rMolecule,2))*(pow(xb-xa,2)+pow(yb-ya,2))))/(pow(xb-xa,2)+pow(yb-ya,2))

def newpoint(xa,xb,s):
    return xa+s*(xb-xa)

def mirror_x(xc,yc,px,py):
    return px-2*((px*xc+py*yc)/(pow(xc,2)+pow(yc,2)))*xc

def mirror_y(xc,yc,px,py):
    return py-2*((px*xc+py*yc)/(pow(xc,2)+pow(yc,2)))*yc

def pos(x,y,z):
    rad=radius(x,y)
    if rad<=r0-rMolecule:
        if z>l0*0.5+rMolecule:
            return 1
        else:
            return 2
    elif rad>r0-rMolecule and rad<r1+rMolecule and z>l0*0.5+rMolecule:
        return 3
    elif rad>=r1+rMolecule:
        if z>l0*0.5+rMolecule:
            return 4
        else:
            return 5

def collision(x,y,z):
    rad=radius(x,y)
    return rad>r0-rMolecule and rad<r1+rMolecule and z<l0*0.5+rMolecule

def derivative(f):
    def compute(x, dx):
        return (f(x+dx) - f(x))/dx
    return compute
def newton(f, x, dx=0.000001, tolerance=0.000001):
    '''f is the function f(x)'''
    df = derivative(f)
    while True:
        x1 = x - f(x)/df(x, dx)
        t = abs(x1 - x)
        if t < tolerance:
            break
        x = x1
    return x

def mirror(a,b,c,d,dx,dy,dz):
    s=-2*(a*dx+b*dy+c*dz-d)/(a**2+b**2+c**2)
    return np.array([dx+s*a,dy+s*b,dz+s*c])

def s_innercylinder(ax,ay,az,bx,by,bz):
    if radius(ax,ay)<=round(r0-rMolecule,2) and radius(bx,by)>r0-rMolecule:
        s=para(ax,ay,bx,by)
        z_new=az+s*(bz-az)
        if z_new<=l0*0.5 and z_new>=-l0*0.5:
            return s
        else:
            return 2.0
    else: return 2.0

def newpoint_innercylinder(ax,ay,az,bx,by,bz,s):
    a, b, c=ax+s*(bx-ax), ay+s*(by-ay), 0
    z_c=az+s*(bz-az)
    dx=bx-a
    dy=by-b
    dz=bz-z_c
    D=mirror(a,b,c,0,dx,dy,dz)
    ex=D[0]+a
    ey=D[1]+b
    ez=D[2]+z_c
    return np.array([ex,ey,ez])

def s_dnatop(ax,ay,az,bx,by,bz):
    if az>l0*0.5+rMolecule and bz<=l0*0.5+rMolecule:
        s=(az-(l0*0.5+rMolecule))/(az-bz)
        x_new, y_new=ax+s*(bx-ax), ay+s*(by-ay)
        if radius(x_new,y_new)<=r1 and radius(x_new,y_new)>=r0:
            return s
        else: return 2.0
    else: return 2.0

def newpoint_dnatop(ax,ay,az,bx,by,bz,s):
    a, b, c=0, 0, 1
    dx=(1-s)*(bx-ax)
    dy=(1-s)*(by-ay)
    dz=(1-s)*(bz-az)
    D=mirror(a,b,c,0,dx,dy,dz)
    ex=D[0]+ax+s*(bx-ax)
    ey=D[1]+ay+s*(by-ay)
    ez=D[2]+az+s*(bz-az)
    return np.array([ex,ey,ez])


def s_innertorus(ax,ay,az,bx,by,bz):
    def function(s):
        return (r0-sqrt((ax+s*(bx-ax))**2+(ay+s*(by-ay))**2))**2+(az+s*(bz-az)-l0*0.5)**2-(rMolecule)**2
    S=np.linspace(0,1,200)
    vz=sign(function(S[0]))
    s_1=2.0
    s_2=2.0
    index=0
    while True:
        if vz != sign(function(S[index+1])):
            s_1=S[index+1]
        index+=1
        if index==S.shape[0]-1 or s_1<2.0:
            break
    if s_1<2.0:
        s_2=s_innertorus_new(ax,ay,az,bx,by,bz)
    return s_2
        

def s_outertorus(ax,ay,az,bx,by,bz):
    def function(s):
        return (r1-sqrt((ax+s*(bx-ax))**2+(ay+s*(by-ay))**2))**2+(az+s*(bz-az)-l0*0.5)**2-(rMolecule)**2
    S=np.linspace(0,1,200)
    vz=sign(function(S[0]))
    s_1=2.0
    s_2=2.0
    index=0
    while True:
        if vz != sign(function(S[index+1])):
            s_1=S[index+1]
        index+=1
        if index==S.shape[0]-1 or s_1<2.0:
            break
    if s_1<2.0:
        s_2=s_outertorus_new(ax,ay,az,bx,by,bz)
    return s_2

def s_outercylinder(ax,ay,az,bx,by,bz):
    if radius(ax,ay)>r1+rMolecule and radius(bx,by)<=r1+rMolecule:
        s=1-para2(bx,by,ax,ay)
        cz=az+s*(bz-az)
        if s>=0 and s<=1 and cz<=l0*0.5 and cz>l1*0.5:
            return s
        else: return 2.0
    else: return 2.0

def newpoint_innertorus(ax,ay,az,bx,by,bz,s):
    cx, cy, cz=ax+s*(bx-ax), ay+s*(by-ay), az+s*(bz-az)
    rad=radius(cx,cy)
    dx=cx/rad*r0
    dy=cy/rad*r0
    dz=l0*0.5
    a, b, c=cx-dx, cy-dy, cz-dz
    d=cx*a+cy*b+cz*c
    D=mirror(a,b,c,d,bx,by,bz)
    return np.array([D[0],D[1],D[2]])

def newpoint_outertorus(ax,ay,az,bx,by,bz,s):
    cx, cy, cz=ax+s*(bx-ax), ay+s*(by-ay), az+s*(bz-az)
    rad=radius(cx,cy)
    dx=cx/rad*r1
    dy=cy/rad*r1
    dz=l0*0.5
    a, b, c=cx-dx, cy-dy, cz-dz
    d=cx*a+cy*b+cz*c
    D=mirror(a,b,c,d,bx,by,bz)
    return np.array([D[0],D[1],D[2]])

def newpoint_outercylinder(ax,ay,az,bx,by,bz,s):
    a, b, c=ax+s*(bx-ax), ay+s*(by-ay), 0
    z_c=az+s*(bz-az)
    dx=bx-a
    dy=by-b
    dz=bz-z_c
    D=mirror(a,b,c,0,dx,dy,dz)
    ex=D[0]+a
    ey=D[1]+b
    ez=D[2]+z_c
    return np.array([ex,ey,ez])

def s_membran(ax,ay,az,bx,by,bz):
    if az>l1*0.5+rMolecule and bz<=l1*0.5+rMolecule:
        s=(az-(l1*0.5+rMolecule))/(az-bz)
        x_new, y_new=ax+s*(bx-ax), ay+s*(by-ay)
        if radius(x_new,y_new)>=r1+rMolecule:
            return s
        else: return 2.0
    else: return 2.0
        
def newpoint_membran(ax,ay,az,bx,by,bz,s):
    a, b, c=0, 0, 1
    dx=(1-s)*(bx-ax)
    dy=(1-s)*(by-ay)
    dz=(1-s)*(bz-az)
    D=mirror(a,b,c,0,dx,dy,dz)
    ex=D[0]+ax+s*(bx-ax)
    ey=D[1]+ay+s*(by-ay)
    ez=D[2]+az+s*(bz-az)
    return np.array([ex,ey,ez])

def s_innertorus_new(ax,ay,az,bx,by,bz):
    az=az-l0*0.5
    bz=bz-l0*0.5
    p=np.array([ax,ay,az])
    d=np.array([bx-ax,by-ay,bz-az])
    R=r0
    r=rMolecule
    a_4=np.inner(d,d)**2
    a_3=4*np.inner(d,d)*np.inner(p,d)
    a_2=4*np.inner(p,d)**2+2*np.inner(d,d)*(np.inner(p,p)-r**2-R**2)+4*R**2*(bz-az)**2
    a_1=4*np.inner(p,d)*(np.inner(p,p)-r**2-R**2)+8*R**2*az*(bz-az)
    a_0=(np.inner(p,p)-r**2-R**2)**2+4*R**2*az**2-4*R**2*r**2
    s=sp.symbols('s')
    roots=np.array(sp.solve(a_4*s**4+a_3*s**3+a_2*s**2+a_1*s+a_0,s))
    t=2.0
    for index in range(roots.shape[0]):
        if abs(sp.im(roots[index]))<=1e-20:
            f=sp.re(roots[index])
            if f>=0 and f<=1:
                t=min(t,f)
    return t

def s_outertorus_new(ax,ay,az,bx,by,bz):
    az=az-l0*0.5
    bz=bz-l0*0.5
    p=np.array([ax,ay,az])
    d=np.array([bx-ax,by-ay,bz-az])
    R=r1
    r=rMolecule
    a_4=np.inner(d,d)**2
    a_3=4*np.inner(d,d)*np.inner(p,d)
    a_2=4*np.inner(p,d)**2+2*np.inner(d,d)*(np.inner(p,p)-r**2-R**2)+4*R**2*(bz-az)**2
    a_1=4*np.inner(p,d)*(np.inner(p,p)-r**2-R**2)+8*R**2*az*(bz-az)
    a_0=(np.inner(p,p)-r**2-R**2)**2+4*R**2*az**2-4*R**2*r**2
    s=sp.symbols('s')
    roots=np.array(sp.solve(a_4*s**4+a_3*s**3+a_2*s**2+a_1*s+a_0,s))
    t=2.0
    for index in range(roots.shape[0]):
        if sp.im(roots[index])<=1e-20:
            f=sp.re(roots[index])
            if f>=0 and f<=1:
                t=min(t,f)
    return t

def s_innertorus_bottom(ax,ay,az,bx,by,bz):
    def function(s):
        return (r0-sqrt((ax+s*(bx-ax))**2+(ay+s*(by-ay))**2))**2+(az+s*(bz-az)+l0*0.5)**2-(rMolecule)**2
    S=np.linspace(0,1,200)
    vz=sign(function(S[0]))
    s_1=2.0
    s_2=2.0
    index=0
    while True:
        if vz != sign(function(S[index+1])):
            s_1=S[index+1]
        index+=1
        if index==S.shape[0]-1 or s_1<2.0:
            break
    if s_1<2.0:
        s_2=s_innertorus_new_bottom(ax,ay,az,bx,by,bz)
    return s_2

def s_innertorus_new_bottom(ax,ay,az,bx,by,bz):
    az=az+l0*0.5
    bz=bz+l0*0.5
    p=np.array([ax,ay,az])
    d=np.array([bx-ax,by-ay,bz-az])
    R=r0
    r=rMolecule
    a_4=np.inner(d,d)**2
    a_3=4*np.inner(d,d)*np.inner(p,d)
    a_2=4*np.inner(p,d)**2+2*np.inner(d,d)*(np.inner(p,p)-r**2-R**2)+4*R**2*(bz-az)**2
    a_1=4*np.inner(p,d)*(np.inner(p,p)-r**2-R**2)+8*R**2*az*(bz-az)
    a_0=(np.inner(p,p)-r**2-R**2)**2+4*R**2*az**2-4*R**2*r**2
    s=sp.symbols('s')
    roots=np.array(sp.solve(a_4*s**4+a_3*s**3+a_2*s**2+a_1*s+a_0,s))
    t=2.0
    for index in range(roots.shape[0]):
        if sp.im(roots[index])<=1e-20:
            f=sp.re(roots[index])
            if f>=0 and f<=1:
                t=min(t,f)
    return t

def s_outertorus_bottom(ax,ay,az,bx,by,bz):
    def function(s):
        return (r1-sqrt((ax+s*(bx-ax))**2+(ay+s*(by-ay))**2))**2+(az+s*(bz-az)+l0*0.5)**2-(rMolecule)**2
    S=np.linspace(0,1,200)
    vz=sign(function(S[0]))
    s_1=2.0
    s_2=2.0
    index=0
    while True:
        if vz != sign(function(S[index+1])):
            s_1=S[index+1]
        index+=1
        if index==S.shape[0]-1 or s_1<2.0:
            break
    if s_1<2.0:
        s_2=s_outertorus_new_bottom(ax,ay,az,bx,by,bz)
    return s_2

def s_outertorus_new_bottom(ax,ay,az,bx,by,bz):
    az=az+l0*0.5
    bz=bz+l0*0.5
    p=np.array([ax,ay,az])
    d=np.array([bx-ax,by-ay,bz-az])
    R=r1
    r=rMolecule
    a_4=np.inner(d,d)**2
    a_3=4*np.inner(d,d)*np.inner(p,d)
    a_2=4*np.inner(p,d)**2+2*np.inner(d,d)*(np.inner(p,p)-r**2-R**2)+4*R**2*(bz-az)**2
    a_1=4*np.inner(p,d)*(np.inner(p,p)-r**2-R**2)+8*R**2*az*(bz-az)
    a_0=(np.inner(p,p)-r**2-R**2)**2+4*R**2*az**2-4*R**2*r**2
    s=sp.symbols('s')
    roots=np.array(sp.solve(a_4*s**4+a_3*s**3+a_2*s**2+a_1*s+a_0,s))
    t=2.0
    for index in range(roots.shape[0]):
        if sp.im(roots[index])<=1e-20:
            f=sp.re(roots[index])
            if f>=0 and f<=1:
                t=min(t,f)
    return t

def newpoint_innertorus_bottom(ax,ay,az,bx,by,bz,s):
    cx, cy, cz=ax+s*(bx-ax), ay+s*(by-ay), az+s*(bz-az)
    rad=radius(cx,cy)
    dx=cx/rad*r0
    dy=cy/rad*r0
    dz=-l0*0.5
    a, b, c=cx-dx, cy-dy, cz-dz
    d=cx*a+cy*b+cz*c
    D=mirror(a,b,c,d,bx,by,bz)
    return np.array([D[0],D[1],D[2]])

def newpoint_outertorus_bottom(ax,ay,az,bx,by,bz,s):
    cx, cy, cz=ax+s*(bx-ax), ay+s*(by-ay), az+s*(bz-az)
    rad=radius(cx,cy)
    dx=cx/rad*r1
    dy=cy/rad*r1
    dz=-l0*0.5
    a, b, c=cx-dx, cy-dy, cz-dz
    d=cx*a+cy*b+cz*c
    D=mirror(a,b,c,d,bx,by,bz)
    return np.array([D[0],D[1],D[2]])

def s_dnabottom(ax,ay,az,bx,by,bz):
    if az<-l0*0.5-rMolecule and bz>=-l0*0.5-rMolecule:
        s=(-l0/2.-rMolecule-az)/(bz-az)
        x_new, y_new=ax+s*(bx-ax), ay+s*(by-ay)
        if radius(x_new,y_new)<=r1 and radius(x_new,y_new)>=r0:
            return s
        else: return 2.0
    else: return 2.0

def newpoint_dnabottom(ax,ay,az,bx,by,bz,s):
    a, b, c=0, 0, 1
    dx=(1-s)*(bx-ax)
    dy=(1-s)*(by-ay)
    dz=(1-s)*(bz-az)
    D=mirror(a,b,c,0,dx,dy,dz)
    ex=D[0]+ax+s*(bx-ax)
    ey=D[1]+ay+s*(by-ay)
    ez=D[2]+az+s*(bz-az)
    return np.array([ex,ey,ez])

def s_outercylinder_bottom(ax,ay,az,bx,by,bz):
    if radius(ax,ay)>r1+rMolecule and radius(bx,by)<=r1+rMolecule:
        s=1-para2(bx,by,ax,ay)
        cz=az+s*(bz-az)
        if s>=0 and s<=1 and cz>=-l0*0.5 and cz<-l1*0.5:
            return s
        else: return 2.0
    else: return 2.0

def s_membran_bottom(ax,ay,az,bx,by,bz):
    if az<-l1*0.5-rMolecule and bz>=-l1*0.5-rMolecule:
        s=(-(l1*0.5+rMolecule)-az)/(bz-az)
        x_new, y_new=ax+s*(bx-ax), ay+s*(by-ay)
        if radius(x_new,y_new)>=r1+rMolecule:
            return s
        else: return 2.0
    else: return 2.0
