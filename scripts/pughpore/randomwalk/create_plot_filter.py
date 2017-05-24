# -*- coding: utf-8 -*-
from matplotlib import gridspec
import math
import matplotlib
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import os
import sys
import nanopores.tools.fields as f

curr = 752.3849
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

number=False

geop = nano.Params(pughpore.params)
hpore=geop.hpore
#fieldsname='eventsnew_both_1_'
#params=dict(avgbind1=23e6,avgbind2=3e4,P_bind1=0.035,P_bind2=3e-1,z0=hpore/2.+0.)
#i=350


def save_fig_filter(params,fieldsname,i):
    data=f.get_fields(fieldsname,**params)
    T_=data["T"][i].load()

    J_=data["J"][i].load()
    T=np.array([])
    J=np.array([])

    for k in range(T_.shape[0]):
        time=T_[k]
        if time==1.:
            T=np.append(T,1.)
            J=np.append(J,J_[k])
        else:
            vec=np.ones(int(time))
            T=np.append(T,vec)
            J=np.append(J,J_[k]*vec)

    T=np.append(T,np.ones(1000000))
    J=np.append(J,np.ones(1000000)*J[0])


    for k in range(1,T.shape[0]):
        T[k]+=T[k-1]

#    T*=1e-9
    J*=1e12

    zero=J[0]
    J-=zero
    from scipy import signal
    b, a = signal.bessel(2, 1./15.*1e-3, analog=False)
    sig_ff = signal.lfilter(b, a, J)
    amp=data["a"][i]
    ftime=np.min(np.where(sig_ff>1e-10))-1
    ampf=abs(np.sum(sig_ff[:ftime]))/float(ftime)/curr*100
    print 'without filter: %.1f' %amp
    print 'with    filter: %.1f' %ampf
    tau_off=float(ftime)
    figname = fieldsname+'_filter_'+'%.8f'%(tau_off*1e-6)+'_%04d'%i+'_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"])+str(params["z0"])
    af=[ampf]
    tf=[tau_off]
    f.save_fields(fieldsname,params,tf3=tf,af3=af)


    plt.figure(figsize=(6,4),dpi=80)
    plt.title('Current signal')
    ax=plt.gca()
    if tau_off<1e3:
        t = np.linspace(0.,tau_off,3)
        fac=1.
        ax.set_xlabel('time [$ns$]')
    elif tau_off<1e6 and tau_off>=1e3:
        t = np.linspace(0.,tau_off*1e-3,3)
        fac = 1e-3
        ax.set_xlabel(r'time [$\mu s$]')
    else:
        t = np.linspace(0.,tau_off*1e-6,3)
        fac = 1e-6
        ax.set_xlabel('time [$ms$]')
    T=T*fac
    plt.plot(T,J+zero,color='#000000',label='Original')
    plt.plot(T,sig_ff+zero,linewidth=2.,color='#ff6600',label='Filtered')
    plt.legend(loc='best')
    yt = np.linspace(580.,760,4)
    ax.set_ylabel(r'A [$pA$]')
    ax.set_yticks(yt)
    ax.set_xticks(t)
    xfmt=FormatStrFormatter('%.1f')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim([-4e-2*tau_off*fac,(1.+4e-2)*tau_off*fac])

    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
    print 'savefig: %s'%figname
    plt.close("all")
#    plt.show()

#save_fig_filter(params,fieldsname,i)
