from sys import path
from scipy import special
from scipy import integrate
import numpy as np
from math import pi, sqrt,exp
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
sims=np.sum(np.load('counter.npy'))
time=5e6
path.append('/home/benjamin/projekt/texfiles/')
from colors import *


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


xaxis = np.linspace(0,time,100) # time [ns]
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
############################################################
fig=plt.figure(num=None, figsize=(12,7), dpi=100)
label1='Sample mean'
label2='Calculated mean'
label3='Sample standard deviation'
label4='Calculated standard deviation'
line1, = plt.plot(xaxis*1e-3,Vecx, label=label1, color=color1, linewidth=2)
line2, = plt.plot(xaxis*1e-3,y1, label=label2, color=color1, linestyle='--')
line3, = plt.plot(xaxis*1e-3,np.sqrt(Vecx2), label=label3, color=color2, linewidth=2)
line4, = plt.plot(xaxis*1e-3,np.sqrt(y2), label=label4, color=color2, linestyle='--')
ax=plt.gca()
ax.set_xlabel('Time [microsec]')
ax.set_ylabel('Distance [nm]')
first_legend = plt.legend(handles=[line1,line3, line2, line4], loc=2)
plt.title('x')
plt.show()
#plt.savefig('compareplots/x.png')

fig=plt.figure(num=None, figsize=(12,7), dpi=100)
label1='Sample mean'
label2='Calculated mean'
label3='Sample standard deviation'
label4='Calculated standard deviation'
line1, = plt.plot(xaxis*1e-3,Vecy, label=label1, color=color1, linewidth=2)
line2, = plt.plot(xaxis*1e-3,y1, label=label2, color=color1, linestyle='--')
line3, = plt.plot(xaxis*1e-3,np.sqrt(Vecy2), label=label3, color=color2, linewidth=2)
line4, = plt.plot(xaxis*1e-3,np.sqrt(y2), label=label4, color=color2, linestyle='--')
ax=plt.gca()
ax.set_xlabel('Time [microsec]')
ax.set_ylabel('Distance [nm]')
first_legend = plt.legend(handles=[line1,line3, line2, line4], loc=2)
plt.title('y')
plt.show()
#plt.savefig('compareplots/y.png')

fig=plt.figure(num=None, figsize=(12,7), dpi=100)
label1='Sample mean'
label2='Calculated mean'
label3='Sample standard deviation'
label4='Calculated standard deviation'
line1, = plt.plot(xaxis*1e-3,Vecz, label=label1, color=color1, linewidth=2)
line2, = plt.plot(xaxis*1e-3,ex, label=label2, color=color1, linestyle='--')
line3, = plt.plot(xaxis*1e-3,np.sqrt(Vecz2), label=label3, color=color2, linewidth=2)
line4, = plt.plot(xaxis*1e-3,np.sqrt(var), label=label4, color=color2, linestyle='--')
ax=plt.gca()
ax.set_xlabel('Time [microsec]')
ax.set_ylabel('Distance [nm]')
first_legend = plt.legend(handles=[line1,line3, line2, line4], loc=2)
plt.title('z')
plt.show()
#plt.savefig('compareplots/z.png')
