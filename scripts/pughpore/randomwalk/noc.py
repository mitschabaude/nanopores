from math import exp, factorial, pow
import matplotlib.pyplot as plt
import numpy as np
import os
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

fieldsname='number_of_collisions'
params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=0.,P_bind2=0.,z0=23.)
color1 = 'blue'
color2 = 'green'
color3 = 'red'
alpha=.2

data=f.get_fields(fieldsname,**params)
data2=f.get_fields(fieldsname+'2',**params)
data3=f.get_fields(fieldsname+'3',**params)
Nc=np.array(data["Nc"])
Nc2=np.array(data2["Nc"])
Nc3=np.array(data3["Nc"])
lam=np.mean(Nc)
lam2=np.mean(Nc2)
lam3=np.mean(Nc3)
print 'len = 14 - lam = '+str(lam)
print 'len =  8 - lam = '+str(lam2)
print 'len =  3 - lam = '+str(lam3)
plt.hist(Nc,20,normed=1,alpha=alpha,color=color1)
plt.hist(Nc2,15,normed=1,alpha=alpha,color=color2)
plt.hist(Nc3,10,normed=1,alpha=alpha,color=color3)
k=np.arange(0,21)
def P(k):
    return pow(lam,k)/(factorial(k))*exp(-lam)
P=np.array([P(x) for x in k])
plt.scatter(k,P,color=color1)
def P2(k):
    return pow(lam2,k)/(factorial(k))*exp(-lam2)
P2=np.array([P2(x) for x in k])
plt.scatter(k,P2,color=color2)
def P3(k):
    return pow(lam3,k)/(factorial(k))*exp(-lam3)
P3=np.array([P3(x) for x in k])
plt.scatter(k,P3,color=color3)
plt.show()
