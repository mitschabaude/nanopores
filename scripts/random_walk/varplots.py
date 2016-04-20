from scipy import special
from scipy import integrate
import numpy as np
from math import pi, sqrt,exp
from matplotlib import pyplot as plt
sims=100
time=6e5


kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
D=(kb*T)/(6*pi*0.5e-9*visc)*1e9 #diffusion[nm^2/ns]
d=2.+5.41
def f(x,t):
    return exp(-x**2/(4*D*t))/sqrt(4*D*pi*t)
def g(x,t):
    return x*(f(x+d,t)+f(x-d,t))
def g2(x,t):
    return x*x*(f(x+d,t)+f(x-d,t))
def expected(t):
    h = lambda x:g(x,t)
    return integrate.quad(h,0,np.inf)[0]
def variance(t):
    h = lambda x:g2(x,t)
    return integrate.quad(h,0,np.inf)[0]


xaxis = np.linspace(0,time,sims) # time [ns]
y1=np.zeros(xaxis.shape[0])
y2=np.array([2*D*t for t in xaxis])
ex=np.array([expected(xaxis[i]) for i in range(1,xaxis.shape[0])])
ex=np.insert(ex,0,d)
var=np.array([variance(xaxis[i])-ex[i]**2 for i in range(1,xaxis.shape[0])])
var=np.insert(var,0,0)

Vecx = np.load('Vecx.npy')/float(sims)
Vecy = np.load('Vecy.npy')/float(sims)
Vecz = np.load('Vecz.npy')/float(sims)
Vecx2 = np.load('Vecx2.npy')/float(sims)
Vecy2 = np.load('Vecy2.npy')/float(sims)
Vecz2 = np.load('Vecz2.npy')/float(sims)
Vecx2=Vecx2-Vecx**2
Vecy2=Vecy2-Vecy**2
Vecz2=Vecz2-Vecz**2
plt.plot(xaxis,Vecx)
plt.plot(xaxis,y1)
plt.title('E(X)')
plt.savefig('compareplots/exp_x.png')
plt.close()
plt.plot(xaxis,Vecy)
plt.plot(xaxis,y1)
plt.title('E(Y)')
plt.savefig('compareplots/exp_y.png')
plt.close()
plt.plot(xaxis,Vecz)
plt.plot(xaxis,ex)
plt.title('E(Z)')
plt.savefig('compareplots/exp_z.png')
plt.close()
plt.plot(xaxis,np.sqrt(Vecx2))
plt.plot(xaxis,np.sqrt(y2))
plt.title('Var(X)')
plt.savefig('compareplots/var_x.png')
plt.close()
plt.plot(xaxis,np.sqrt(Vecy2))
plt.plot(xaxis,np.sqrt(y2))
plt.title('Var(Y)')
plt.savefig('compareplots/var_y.png')
plt.close()
plt.plot(xaxis,np.sqrt(Vecz2))
plt.plot(xaxis,np.sqrt(var))
plt.title('Var(Z)')
plt.savefig('compareplots/var_z.png')
