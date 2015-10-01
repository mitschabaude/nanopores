"""
    A simple example of an animated plot... In 3D!
    """
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

import sys

# bad hack to allow arbitrary nm in params_geo
from nanopores import import_vars
params = import_vars("nanopores.H_cyl_geo.params_geo")
for x in params:
    exec("%s = %s*1e-9/%s" %(x, params[x], params['nm']))

import scipy as sp
from array import array

X=np.load('X.npy')
Y=np.load('Y.npy')
Z=np.load('Z.npy')
l=X.shape[0]
X=X.tolist()
Y=Y.tolist()
Z=Z.tolist()

show_all=False
#innercylinder
alphb=0.5
#rest
alpha=0.1

def update_lines(num, dataLines, lines) :
    for line, data in zip(lines, dataLines) :
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2,:num])
    return lines

def update(num, dataLines, lines):
    return update_lines(num, dataLines, lines)


fig = plt.figure()
ax = p3.Axes3D(fig)



data = [np.array([X,Y,Z])]


# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1],c='green')[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([-7.5, 7.5])
ax.set_xlabel('X')

ax.set_ylim3d([-7.5, 7.5])
ax.set_ylabel('Y')

ax.set_zlim3d([-7.5, 7.5])
ax.set_zlabel('Z')

ax.set_title('Random Walk')





# innerCylinder

x=np.linspace(-r0*1e9, r0*1e9, 100)
z=np.linspace(-l0/2.0*1e9, l0/2.0*1e9, 100)
Xic, Zic=np.meshgrid(x, z)
Yic = np.sqrt((r0*1e9)**2-Xic**2)

# Draw parameters
rstride = 40
cstride = 40
ax.plot_surface(Xic, Yic, Zic, alpha=alphb, rstride=rstride, cstride=cstride)
if show_all:
    ax.plot_surface(Xic, -Yic, Zic, alpha=alphb, rstride=rstride, cstride=cstride)

# outerupperCylinder

x=np.linspace(-r1*1e9, r1*1e9, 100)
z=np.linspace(l1/2.0*1e9, l0/2.0*1e9, 100)
Xouc, Zouc=np.meshgrid(x, z)
Youc = np.sqrt((r1*1e9)**2-Xouc**2)

# Draw parameters
rstride = 50
cstride = 50
ax.plot_surface(Xouc, Youc, Zouc, alpha=alpha, rstride=rstride, cstride=cstride)
if show_all:
    ax.plot_surface(Xouc, -Youc, Zouc, alpha=alpha, rstride=rstride, cstride=cstride)



# outerlowerCylinder

x=np.linspace(-r1*1e9, r1*1e9, 100)
z=np.linspace(-l0/2.0*1e9, -l1/2.0*1e9, 100)
Xodc, Zodc=np.meshgrid(x, z)
Yodc = np.sqrt((r1*1e9)**2-Xodc**2)

# Draw parameters
#rstride = 40
#cstride = 40
ax.plot_surface(Xodc, Yodc, Zodc, alpha=alpha, rstride=rstride, cstride=cstride)
if show_all:
    ax.plot_surface(Xodc, -Yodc, Zodc, alpha=alpha, rstride=rstride, cstride=cstride)


# upperTop

u = np.linspace(0, np.pi, 100)
v = np.linspace(r0*1e9, r1*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)+l0/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)

# lowerTop

u = np.linspace(0, np.pi, 100)
v = np.linspace(r0*1e9, r1*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)-l0/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)

# upperMembran

u = np.linspace(0, np.pi, 100)
v = np.linspace(r1*1e9, Rz*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)+l1/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)

# lowerMembran

u = np.linspace(0, np.pi, 100)
v = np.linspace(r1*1e9, Rz*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)-l1/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)



# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update, l, fargs=(data, lines),
                                   interval=1, blit=False)
#line_ani.save('animation_path.mp4', fps=50, writer="avconv", codec="libx264")

#plt.close(fig)
plt.show()
