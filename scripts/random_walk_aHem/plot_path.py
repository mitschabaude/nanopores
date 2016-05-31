from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from aHem_array_2d import *
def radius(x,y):
    return sqrt(x**2+y**2)

a=np.load('hbond.npy')
X=np.load('x.npy')
Y=np.load('y.npy')
Z=np.load('z.npy')
r_h=np.load('r_h.npy')
z_h=np.load('z_h.npy')
def radius(x,y):
    return sqrt(x**2+y**2)
Rad=np.array([radius(X[i],Y[i]) for i in range(X.shape[0])])
size=X_aHem_2d.shape[0]

leftend=10.
x_mem=np.linspace(X_aHem_2d[18][0],leftend,100)
y_mem=np.zeros(x_mem.shape[0])+X_aHem_2d[18][1]
X=np.zeros(size+1)
Y=np.zeros(size+1)
X_=np.zeros(a.shape[0]+1)
Y_=np.zeros(a.shape[0]+1)
for index in range(size):
    X[index]=X_aHem_2d[index][0]
    Y[index]=X_aHem_2d[index][1]
for index in range(a.shape[0]):
    X_[index]=a[index][0]
    Y_[index]=a[index][1]
X[size]=X[0]
Y[size]=Y[0]
X_[a.shape[0]]=X_[0]
Y_[a.shape[0]]=Y_[0]
#fig=plt.figure(figsize=(5.7,11), dpi=400)
plt.plot(X,Y,linewidth='2',color='blue')
plt.scatter(X,Y,50,color='blue')
plt.plot(X_,Y_,linewidth=2,color='green')
plt.scatter(X_,Y_,50,color='green')
plt.plot(x_mem,y_mem,color='black',linewidth=1)
plt.plot(Rad,Z,color='red')
plt.scatter(r_h,z_h,50,color='red')
#plt.scatter([0.,0.,2.,2.],[0.,-2.,-2.,0.],50,color='red')
#plt.savefig('array.png')
plt.show()
