import sys
#import matplotlib.pyplot as plt
from random import gauss
from math import sqrt, ceil, floor, pi
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore

#eps = 5e-2
#L = 1.

geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")
kT = physp.kT
eta = physp.eta
l3 =        geop.l3
rMolecule = geop.rMolecule
Dmol = kT/(6.*pi*eta*rMolecule*1e-9) # [m^2/s]
D = Dmol*1e9 # from [m^2/s] to [nm^2/ns]
fac = sqrt(2*D) # correction factor

L = l3-rMolecule
eps = L*1e-1

taufacinv = int(sys.argv[1])
tau = eps**2/(2*D)*5e-2*(1./taufacinv)
iter = 50000000
len=tau*iter
sqrttau = sqrt(tau)



attempt = 0

try:
    np.load('does_not_exist.npy')
    X = np.load('X.npy')
    Y = np.load('Y.npy')
except:
    Y = np.zeros(iter)
    for i in range(iter-1):
        if i%int(iter/100.)==0:
            print('%.0f %%'%(100*float(i)/float(iter)))
        xi = gauss(0.,1.)
        Y[i+1]=Y[i]+sqrttau*fac*xi
    X = np.linspace(0.,len,iter)
    np.save('X',X)
    np.save('Y',Y)

maxL = ceil(np.max(Y))
minL = floor(np.min(Y))
L_ = np.concatenate((-np.arange(L,abs(minL)+L,L)[::-1],np.arange(0.,maxL+L,L)))
#iold = 0
sig = 0.
ffa = False
try:
    np.load('does_not_exist.npy')
    Exp=np.load('Exp.npy')
    XA = np.load('XA.npy')
    YA = np.load('YA.npy')
except:
    XA = np.array([])
    YA = np.array([])
    Exp = np.zeros(100)
    exptau=len/100.
    j=1
    for i in range(1,X.shape[0]):
        if i%int(X.shape[0]/100.)==0 and i!=0:
            print('%.0f %%'%(100*float(i)/float(X.shape[0])))
        dist = np.min(np.absolute(L_-Y[i]))
        Li = np.where(abs(L_-Y[i])==dist)[0][0]
        if dist<=eps and np.sign(Y[i]-L_[Li])!=sig and ffa:
            attempt+=1
            XA = np.append(XA,X[i])
            YA = np.append(YA,Y[i])
            ffa = False
        if dist>eps:
            sig = np.sign(Y[i]-L_[Li])
#            if np.sign(Y[iold]-L_[Li])!=sig:
#                attempt+=1
#                plt.plot([X[iold],X[i]],[Y[iold],Y[i]],color='#00ff00',linewidth=2.)
            ffa = True
#            iold = i
        if X[i]>=exptau*j:
            Exp[j-1]=attempt/(j*exptau)
            j+=1
    Exp[-1]=attempt/(len)
    np.save('Exp',Exp)
    np.save('XA',XA)
    np.save('YA',YA)
attemptrate = Exp[-1]
print('L = %.2e; eps = %.2e, D = %.2e'%(L,eps,D))
print('tau = %.3e'%tau)
theo = 2*D/(L*eps)
print('analytic attempt rate = %.2e'%theo)
print('numeric attemptrate   = %.3e'%attemptrate)

# plot stuff
#for i in L_:
#    plt.plot([0.,iter*tau],[i,i],'r--')
#    plt.fill_between([0.,iter*tau],[float(i)+eps,float(i)+eps],[float(i)-eps,float(i)-eps],color='#ff0000',alpha=.2)
#plt.plot(X,Y,color='#000000')
#plt.tight_layout()
#plt.plot(XA,YA,'ro')
#plt.show()
#plt.plot(np.linspace(1,len,100),Exp,color='#0000ff')
#plt.plot([0.,len],[attemptrate,attemptrate],color='#0000ff')
#plt.show()
