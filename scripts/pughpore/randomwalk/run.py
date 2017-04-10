from dolfin import *
import sys
from random import gauss, expovariate
import math
from math import atan, pi, atan2, sqrt
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
from get_F import Force, Current
from get_D import Dx, Dy, Dz, dxDx, dyDy, dzDz, dis
import os
from time import time as timer


HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)

def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")

kT = physp.kT
eta = physp.eta

l0 =        geop.l0
l1 =        geop.l1
l2 =        geop.l2
l3 =        geop.l3
l4 =        geop.l4
hpore =     geop.hpore
hmem =      geop.hmem
h2 =        geop.h2
h1 =        geop.h1
h4 =        geop.h4
rMolecule = geop.rMolecule
beps = (l3 - rMolecule)*1e-1


Dmol = kT/(6.*math.pi*eta*rMolecule*1e-9) # [m^2/s]
gamma = (6.*math.pi*eta*rMolecule) #friction [microgramm/s]
maxiter = 1e6 # [ns]
tau = 1. # [ns]
C = tau/gamma*1e9 # [s^2/kg * 1e9 nm/m]
coeff = math.sqrt(2*Dmol*1e9*tau) # [nm]

def run(params,fieldsname,outcome,outside,b1,b2):
    def area1(x,y,z):
        for seg in b1:
            h=np.array([p[1] for p in seg])
            if np.min(h)<=z and z<=np.max(h):
                return True
        return False
    def area2(x,y,z):
        for seg in b2:
            h=np.array([p[1] for p in seg])
            if np.min(h)<=z and z<=np.max(h):
                return True
        return False
    z0 = params["z0"]
    avgbind1=params["avgbind1"]
    P_bind1=params["P_bind1"]
    avgbind2=params["avgbind2"]
    P_bind2=params["P_bind2"]
    should_restart = True
    while should_restart:
        should_restart = False
        X = np.array([0.])
        Y = np.array([0.])
        Z = np.array([z0])
        J1 = np.array([])
        T = np.array([])
        Nc = 0
        ffa = True
        i=0
        ood = False
        while i<maxiter and Z[-1]>=-hpore/2.-4.:
            add=tau
            xi_x=gauss(0.,1.)
            xi_y=gauss(0.,1.)
            xi_z=gauss(0.,1.)
            arg = argument(X[-1],Y[-1],Z[-1])
            F = Force(X[-1],Y[-1],Z[-1])
            D = [Dx(arg)*1e9,Dy(arg)*1e9,Dz(arg)*1e9]
            dD = [dxDx(arg)*1e9,dyDy(arg)*1e9,dzDz(arg)*1e9]
    #        x_new = X[-1] + coeff*xi_x*math.sqrt(abs(Dxfac)) + C*Force[0]*Dxfac + DDx*tau*Dmol
    #        y_new = Y[-1] + coeff*xi_y*math.sqrt(abs(Dyfac)) + C*Force[1]*Dyfac + DDy*tau*Dmol
    #        z_new = Z[-1] + coeff*xi_z*math.sqrt(abs(Dzfac)) + C*Force[2]*Dzfac + DDz*tau*Dmol
    #        x_new = X[-1] + coeff*xi_x + C*Force[0]
    #        y_new = Y[-1] + coeff*xi_y + C*Force[1]
    #        z_new = Z[-1] + coeff*xi_z + C*Force[2]
            x_new = X[-1] + sqrt(2*D[0]*tau)*xi_x + F[0]*D[0]*1e-9*tau/kT+dD[0]*tau
            y_new = Y[-1] + sqrt(2*D[1]*tau)*xi_y + F[1]*D[1]*1e-9*tau/kT+dD[1]*tau
            z_new = Z[-1] + sqrt(2*D[2]*tau)*xi_z + F[2]*D[2]*1e-9*tau/kT+dD[2]*tau
            if dis(argument(x_new,y_new,z_new)) < rMolecule:
                x_new = X[-1]
                y_new = Y[-1]
                z_new = Z[-1]
                if ffa and area2(0.,0.,Z[-1]): Nc+=1
                if ffa and np.random.binomial(1,P_bind1)==1 and area2(0.,0.,Z[-1]):
                    add+=expovariate(lambd=1./avgbind1)
                elif ffa and np.random.binomial(1,P_bind2)==1 and area1(0.,0.,Z[-1]):
                    add+=expovariate(lambd=1./avgbind2)
                else:
                    add+=0.
                ffa = False
            elif dis(argument(x_new,y_new,z_new)) < rMolecule + beps:
                pass
            else:
                ffa = True
            X = np.append(X,x_new)
            Y = np.append(Y,y_new)
            Z = np.append(Z,z_new)
            if abs(Z[-1])>35. or abs(X[-1])>10. or abs(Y[-1])>10.:
                print 'Out of domain!'
                ood = True
                if not outside:
                    should_restart = True
                    print 'restart!'
                break
            Jx=Current(X[-1],Y[-1],Z[-1])
            if math.isnan(Jx):
                if add<=tau:
                    Jx = J1[-1]
                else:
                    print 'current at binding position is NaN!!!'
                    print 'current = %.1e A'%Jx
                    print 'X = %.8f'%X[-1]
                    print 'Y = %.8f'%Y[-1]
                    print 'Z = %.8f'%Z[-1]
                    print 'add = %.2f nanoseconds'%add
                    exit()
            J1=np.append(J1,Jx)
            T =np.append(T,add)
            i+=1
    if i>=maxiter:
        print 'randomwalk: more than 1e6 steps!'
    fields.save_fields(fieldsname,params,Nc=[Nc])
    if outcome=='type' or outcome=='both':
        tau_off = np.sum(T)*1e-6
        curr = 7.523849e-10
        amp = (curr-np.inner(T*1e-6,J1)/tau_off)/curr*100.
        if math.isnan(amp):
            np.save('T',T)
            np.save('J1',J1)
            file=open('nanerror.txt','w')
            file.write('tau_off = %.10f\n'% tau_off)
            file.write('amp = %.10f\n'% amp)
            file.close()
            exit()
        t=[tau_off]
        a=[amp]
        if ood:
            ood=[1]
        else:
            ood=[0]
            
        fields.save_fields(fieldsname,params,t=t,a=a,ood=ood)
    if outcome=='traj' or outcome=='both':
        X=[X]
        Y=[Y]
        Z=[Z]
        T=[T]
        J1=[J1]
        fields.save_fields(fieldsname,params,X=X, Y=Y, Z=Z, T=T, J=J1)
