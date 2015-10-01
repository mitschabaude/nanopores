import random,math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab

from nanopores import *

# bad hack to allow arbitrary nm in params_geo
params = import_vars("nanopores.W_3D_geo.params_geo")
for x in params:
    exec("%s = %s*1e0/%s" %(x, params[x], params['nm']))
    
height=(lsin+lau+lsam)*nm
r1=height/math.tan(angle*pi/180.)

mpl.rcParams['legend.fontsize']=10

fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax = fig.add_subplot(111, projection='3d')
ax=Axes3D(fig)
#plt.hold(True)


X=np.load('X.npy')
Y=np.load('Y.npy')
Z=np.load('Z.npy')


'''
X_r=np.load('X.npy')
Y_r=np.load('Y.npy')
Z_r=np.load('Z.npy')
#X_C=np.load('X_C.npy')
#Y_C=np.load('Y_C.npy')
#Z_C=np.load('Z_C.npy')
Z=np.zeros(Z_r.shape[0])
Y=np.zeros(Z_r.shape[0])
X=np.zeros(Z_r.shape[0])
for index in range(Z_r.shape[0]):
    Z[index]=Z_r[Z_r.shape[0]-1-index]
    X[index]=X_r[X_r.shape[0]-1-index]
    Y[index]=Y_r[Y_r.shape[0]-1-index]
'''

#################
show_all=True
show_schichten=False
alpha=1.0
color='#d3d3d3'
i=0
while i<=0:#X.shape[0]-1:
    t = float(i)/float(X.shape[0]-1)
    ax.view_init((1-t)*32 + t*6 , (1-t)*(-104) + t*(-99)) # from 32,-104 to 6,-99
    #ax.view_init(6,-99)
    ax.axis('off')
    lim=90
    ax.set_xlim3d(-lim,lim)
    ax.set_ylim3d(-lim,lim)
    ax.set_zlim3d(-80,80)

    # Sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = rMolecule*nm*np.outer(np.cos(u), np.sin(v))
    y = rMolecule*nm*np.outer(np.sin(u), np.sin(v))
    z = rMolecule*nm*np.outer(np.ones(np.size(u)), np.cos(v))
    #ax.plot_surface(x+X[i], y+Y[i], z+Z[i],  rstride=10, cstride=10, color='r',alpha=1.0)


    # lowerTop = lowerMembranNEW

    u = np.linspace(0, np.pi, 100)
    v = np.linspace(r0*nm, R*nm, 10)

    x = np.outer(v, np.cos(u))
    y = np.outer(v, np.sin(u))
    z = np.zeros(shape=x.shape)-height/2.0*nm
    #ax.plot_surface(x, y, z,  rstride=15, cstride=15, color=color, alpha=alpha,linewidth=0)


    #background

    x=np.linspace(-R*nm, R*nm, 100)
    z=np.linspace(-height/2.0, height/2.0, 100)
    z_sin=np.linspace(-height/2.0,-height/2.0+lsin*nm,100)
    z_au=np.linspace(-height/2.0+lsin*nm,-height/2.0+lsin*nm+lau*nm,100)
    z_sam=np.linspace(-height/2.0+lsin*nm+lau*nm,height/2.0,100)
    Xic, Zic=np.meshgrid(x, z)
    Xic, Zic_sin=np.meshgrid(x, z_sin)
    Xic, Zic_au=np.meshgrid(x, z_au)
    Xic, Zic_sam=np.meshgrid(x, z_sam)
    Yic = np.sqrt((R*nm)**2-Xic**2)

    # Draw parameters
    #ax.plot_surface(Xic, Yic, Zic, alpha=alpha, rstride=10, cstride=3, color=color, linewidth=0)
    if show_all:
        ax.plot_surface(Xic, -Yic, Zic_sin, alpha=alpha, rstride=10, cstride=3, color='#00cbff', linewidth=0)
        ax.plot_surface(Xic, -Yic, Zic_au, alpha=alpha, rstride=10, cstride=3, color='#ff7f00', linewidth=0)
        ax.plot_surface(Xic, -Yic, Zic_sam, alpha=alpha, rstride=10, cstride=3, color='#d3d3d3', linewidth=0)





    # innercone

    # Set up the grid in polar
    theta = np.linspace(0,np.pi,90)
    r = np.linspace(r0*nm,r1*nm,60)
    T, R_cone = np.meshgrid(theta, r)

    # Then calculate X, Y, and Z
    X_cone = -R_cone * np.cos(T)
    Y_cone = R_cone * np.sin(T)
    Z_cone = np.sqrt(X_cone**2 + Y_cone**2)-height/2.0*nm-r0*nm

    # Set the Z values outside your range to NaNs so they aren't plotted
    Z_cone[Z_cone < -height/2.0*nm] = np.nan
    Z_cone[Z_cone > +height/2.0*nm] = np.nan
    ax.plot_surface(X_cone, Y_cone, Z_cone, rstride=3, cstride=1, color=color,alpha=alpha,linewidth=0.2)
    if show_all:
        ax.plot_surface(X_cone, -Y_cone, Z_cone, rstride=3, cstride=1, color=color,alpha=alpha,linewidth=0.2)



    # upperTop = upperMembranNEW

    u = np.linspace(0, np.pi, 100)
    v = np.linspace(r1*nm-3, R*nm, 50)

    x = np.outer(v, np.cos(u))
    y = np.outer(v, np.sin(u))
    z = np.zeros(shape=x.shape)+height/2.0*nm
    ax.plot_surface(x, y, z,  rstride=15, cstride=15, color=color, alpha=alpha,linewidth=0.2)
    if show_all:
        ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color=color, alpha=alpha,linewidth=0.2)



    #schichten
    #SIN
    x=[r0*nm+25,r1*nm-3,R*nm,R*nm,r1*nm-3,r0*nm+25]
    x_=[-r0*nm-25,-r1*nm+3,-R*nm,-R*nm,-r1*nm+3,-r0*nm-25]
    y=[0,0,0,0,0,0]
    z=[-height/2.,-height/2.,-height/2.,-height/2.+lsin*nm,-height/2.+lsin*nm,-height/2.]
    verts1=[zip(x,y,z)]
    verts2=[zip(x_,y,z)]
    if show_schichten:
        ax.add_collection3d(Poly3DCollection(verts1,color='#00cbff'))
        ax.add_collection3d(Poly3DCollection(verts2,color='#00cbff'))
    #AU
    x=[r0*nm+25,r1*nm-3,R*nm,R*nm,r1*nm-3,r0*nm+25-22,r0*nm+25]
    x_=[-r0*nm-25,-r1*nm+3,-R*nm,-R*nm,-r1*nm+3,-r0*nm-25+22,-r0*nm-25]
    y=[0,0,0,0,0,0,0]
    z=[-height/2.,-height/2.+lsin*nm,-height/2.+lsin*nm,height/2.-lsam*nm,height/2.-lsam,-height/2.,-height/2.]
    verts1=[zip(x,y,z)]
    verts2=[zip(x_,y,z)]
    if show_schichten:
        ax.add_collection3d(Poly3DCollection(verts1,color='#ff7f00'))
        ax.add_collection3d(Poly3DCollection(verts2,color='#ff7f00'))

    #SAM
    x=[r0*nm+25-22,r1*nm-3,R*nm,R*nm,r1*nm-3,r0*nm,r0*nm+25-22]
    x_=[-r0*nm-25+22,-r1*nm+3,-R*nm,-R*nm,-r1*nm+3,-r0*nm,-r0*nm-25+22]
    y=[0,0,0,0,0,0,0]
    z=[-height/2.,height/2.-lsam*nm,height/2.-lsam*nm,height/2.,height/2.,-height/2.]
    verts1=[zip(x,y,z)]
    verts2=[zip(x_,y,z)]
    if show_schichten:
        ax.add_collection3d(Poly3DCollection(verts1,color='#d3d3d3'))
        ax.add_collection3d(Poly3DCollection(verts2,color='#d3d3d3'))

    # Sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = rMolecule*nm*np.outer(np.cos(u), np.sin(v))
    y = rMolecule*nm*np.outer(np.sin(u), np.sin(v))
    z = rMolecule*nm*np.outer(np.ones(np.size(u)), np.cos(v))
    #ax.plot_surface(x+X[i], y+Y[i], z+Z[i],  rstride=10, cstride=10, color='r',alpha=1.0)
    #############################################
    plt.show()
    #fig.savefig('movie/%i.jpg'%i)
    fig.savefig('start2.jpg')
    print "\x1b[A"*1,"\r",
    print "%d out of %d"% (i, X.shape[0]-1)
    i+=1
    plt.cla()

