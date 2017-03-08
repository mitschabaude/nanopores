from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
from nanopores.models.pughpore import polygon
from nanopores.models.pughpoints import plot_polygon
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

up = nano.Params(pughpore.params, k=3)
hpore=up.hpore

def save_fig(params, fieldsname, i):
    data=f.get_fields(fieldsname,**params)
    X = np.array(data["X"][i])
    Y = np.array(data["Y"][i])
    Z = np.array(data["Z"][i])
    T = np.array(data["T"][i])
    J = np.array(data["J"][i])
    curr = 7.523849e-10
    bind1 = np.where(T>1e6)
    bind2 = np.intersect1d(np.where(T<=1e6),np.where(T>100.))
    amplitude = curr-np.inner(J,T)/np.sum(T)
    for i in range(1,T.shape[0]):
        T[i]=T[i]+T[i-1]
    tau_off=T[-1]
    figname = fieldsname+'_traj_'+'%.3f'%(tau_off*1e-6)+'_%.1e_%.1e_%.1e_%.1e_'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"],)+str(params["z0"])
    #print 'tau_off = %.3f ms'% (tau_off*1e-6)
    #print 'amplitude = %.0f pA'% amplitude
    #J=np.append(np.array([curr,curr]),J)
    #J=np.append(J,np.array([curr,curr]))
    J=J*1e12
    #T=np.append(np.array([-1e6,0.]),T)
    #T=np.append(T,np.array([tau_off,1e6+tau_off]))
    T=T*1e-6

    fig=plt.figure(figsize=(10,5),dpi=80)
    gs=gridspec.GridSpec(1, 2, width_ratios=[6e-1,1])

    ax0=plt.subplot(gs[0])
    plt.title('Projected trajectory')
    ax0.plot(X,Z,linewidth=1.,c='#0000ff')
    ax0.scatter(X[bind1],Z[bind1],s=200,marker='h',c='#ff0000',linewidth=0.)
    ax0.scatter(X[bind2],Z[bind2],s=100,marker='h',c='#00ff00',linewidth=0.)
    ax=plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-25.,25.])
    ax.set_ylim([-30.,40.])
    ax.set_xticks([-20,0,20])
    ax.set_yticks([-30,0,30])
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    plot_polygon(ax,polygon())


    ax1=plt.subplot(gs[1])
    plt.title('Current Signal')
    ax1.plot(T,J,linewidth=2.,color='#0000ff')
    ax=plt.gca()
    t = np.linspace(0.,tau_off*1e-6,5)
    ax.set_xlabel('Time [ms]')
    ax.set_ylabel('Current [pA]')
    ax.set_xticks(t)
    xfmt=FormatStrFormatter('%.1e')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim([-.05*tau_off*1e-6,1.05*tau_off*1e-6])


    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
    print 'savefig: %s'%figname
    plt.close("all")
