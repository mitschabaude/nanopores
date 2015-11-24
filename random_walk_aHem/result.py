import random,math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
# from nanopores import *

exit_x = np.loadtxt('exit_x.txt',delimiter='\n')
exit_y = np.loadtxt('exit_y.txt',delimiter='\n')
exit_z = np.loadtxt('exit_z.txt',delimiter='\n')

fig = plt.figure()

ax=Axes3D(fig)
ax.scatter(exit_x,exit_y,exit_z,color='red')

show_half=False
fac=2
if show_half:
    fac=1

X_aHem = np.array([[2.16, 0.0, 0.0],
                      [2.77, 0.0, -0.19],
                      [3.24, 0.0, -0.1 ],
                      [3.59, 0.0, -0.1 ],
                      [3.83, 0.0, -0.35],
                      [3.84, 0.0, -0.8 ],
                      [3.67, 0.0, -1.34],
                      [3.73, 0.0, -1.96],
                      [3.93, 0.0, -2.31],
                      [4.23, 0.0, -2.67],
                      [4.44, 0.0, -2.81],
                      [4.33, 0.0, -3.25],
                      [4.01, 0.0, -3.5 ],
                      [3.99, 0.0, -3.67],
                      [4.11, 0.0, -3.94],
                      [4.39, 0.0, -4.12],
                      [4.44, 0.0, -4.52],
                      [4.73, 0.0, -4.86],
                      [4.96, 0.0, -5.41],
                      [4.89, 0.0, -5.87],
                      [4.63, 0.0, -6.44],
                      [4.43, 0.0, -6.96],
                      [4.07, 0.0, -7.32],
                      [3.71, 0.0, -7.51],
                      [3.46, 0.0, -7.36],
                      [3.41, 0.0, -7.1 ],
                      [3.31, 0.0, -6.9 ],
                      [3.04, 0.0, -6.87],
                      [2.73, 0.0, -6.73],
                      [2.41, 0.0, -6.6 ],
                      [2.17, 0.0, -6.41],
                      [1.97, 0.0, -6.23],
                      [1.84, 0.0, -6.03],
                      [1.76, 0.0, -5.87],
                      [1.54, 0.0, -5.87],
                      [1.4 , 0.0, -5.96],
                      [1.31, 0.0, -6.16],
                      [1.39, 0.0, -6.57],
                      [1.6 , 0.0, -6.81],
                      [1.71, 0.0, -7.09],
                      [1.76, 0.0, -7.32],
                      [1.67, 0.0, -7.65],
                      [1.44, 0.0, -7.81],
                      [1.49, 0.0, -8.06],
                      [1.56, 0.0, -8.36],
                      [1.44, 0.0, -8.61],
                      [1.43, 0.0, -8.79],
                      [1.44, 0.0, -9.1 ],
                      [1.6 , 0.0, -9.48],
                      [1.74, 0.0, -9.84],
                      [1.63, 0.0, -10.0],
                      [1.47, 0.0, -10.19],
                      [1.26, 0.0, -10.21],
                      [1.07, 0.0, -10.05],
                      [1.03, 0.0, -9.76],
                      [1.09, 0.0, -9.44],
                      [1.07, 0.0, -9.02],
                      [0.86, 0.0, -8.79],
                      [0.64, 0.0, -8.68],
                      [0.63, 0.0, -8.36],
                      [0.8 , 0.0, -8.22],
                      [0.81, 0.0, -7.93],
                      [0.89, 0.0, -7.71],
                      [1.04, 0.0, -7.51],
                      [1.1 , 0.0, -7.25],
                      [0.91, 0.0, -7.02],
                      [0.91, 0.0, -6.76],
                      [0.91, 0.0, -6.48],
                      [0.69, 0.0, -6.25],
                      [0.69, 0.0, -6.  ],
                      [0.66, 0.0, -5.68],
                      [0.59, 0.0, -5.36],
                      [0.53, 0.0, -5.12],
                      [0.54, 0.0, -4.92],
                      [0.79, 0.0, -4.84],
                      [1.03, 0.0, -4.89],
                      [1.21, 0.0, -4.7 ],
                      [1.36, 0.0, -4.42],
                      [1.49, 0.0, -4.16],
                      [1.66, 0.0, -3.92],
                      [1.66, 0.0, -3.7 ],
                      [1.8 , 0.0, -3.41],
                      [2.  , 0.0, -3.22],
                      [1.91, 0.0, -2.93],
                      [1.8 , 0.0, -2.71],
                      [1.56, 0.0, -2.55],
                      [1.46, 0.0, -2.38],
                      [1.3 , 0.0, -2.19],
                      [1.21, 0.0, -1.93],
                      [1.09, 0.0, -1.64],
                      [0.9 , 0.0, -1.45],
                      [0.8 , 0.0, -1.28],
                      [0.84, 0.0, -1.  ],
                      [1.  , 0.0, -0.8 ],
                      [1.26, 0.0, -0.64],
                      [1.7 , 0.0, -0.31]])

x_aHem, z_aHem = np.zeros_like(X_aHem), np.zeros_like(X_aHem)
for index in range(x_aHem.shape[0]):
  x_aHem[index], z_aHem[index] = X_aHem[index][0], X_aHem[index][2]
x_aHem = np.append(x_aHem,x_aHem[0])
z_aHem = np.append(z_aHem,z_aHem[0])
angle = np.linspace(0,fac*pi,50)
theta, mesh_x = np.meshgrid(angle,x_aHem)
theta, mesh_z = np.meshgrid(angle,z_aHem)

X_plot = mesh_x*np.cos(theta)
Y_plot = mesh_x*np.sin(theta)
Z_plot = mesh_z
ax.plot_wireframe(X_plot, Y_plot, Z_plot,color='blue', rstride = 1, cstride = 1, alpha=1.0)

if show_half:
  x=x_aHem.tolist()
  left=-x_aHem
  x_=left.tolist()
  y=np.zeros_like(x_aHem).tolist()
  z=z_aHem.tolist()
  verts1=[zip(x,y,z)]
  verts2=[zip(x_,y,z)]
  ax.add_collection3d(Poly3DCollection(verts1,color='blue',alpha=1.0))
  ax.add_collection3d(Poly3DCollection(verts2,color='blue',alpha=1.0))
ax.set_xlim3d(-5,5)
ax.set_ylim3d(-5,5)
ax.set_zlim3d(-5,5)
plt.show()
