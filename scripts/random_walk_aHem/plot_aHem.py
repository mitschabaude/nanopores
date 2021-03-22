import matplotlib.pyplot as plt
import numpy as np
from aHem_array_2d import *
a=np.load('hbond.npy')


def radius(x,y):
    return sqrt(x**2+y**2)
def det(a,b,c,d):
	return a*d-b*c

size=X_aHem_2d.shape[0]
hbond=X_aHem_2d

A=0
for i in range(size):
    A+=(hbond[i-1][0]*hbond[i][1]-hbond[i][0]*hbond[i-1][1])
A*=-0.5
print(A)
Cx=0
Cy=0
for i in range(size):
    Cx+=(hbond[i-1][0]+hbond[i][0])*(hbond[i-1][0]*hbond[i][1]-hbond[i][0]*hbond[i-1][1])
    Cy+=(hbond[i-1][1]+hbond[i][1])*(hbond[i-1][0]*hbond[i][1]-hbond[i][0]*hbond[i-1][1])
Cx*=1./(6*A)
Cy*=1./(6*A)
shift=np.array([[Cx,Cy] for i in range(size)])

hbond=hbond+shift
hbond2=hbond*1.2
#hbond=hbond-shift


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
#plt.scatter([0.,0.,2.,2.],[0.,-2.,-2.,0.],50,color='red')
#plt.savefig('array.png')
plt.show()
