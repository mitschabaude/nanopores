from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
#from matplotlib import rc
#rc('text', usetex=True)
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

def save_fig(params,fieldsname,i):
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
    figname = fieldsname+'_traj_'+'%.3f'%(tau_off*1e-6)+'_%i'%i+'_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["avgbind3"],params["P_bind1"],params["P_bind2"],params["P_bind3"])+str(params["z0"])
    #print 'tau_off = %.3f ms'% (tau_off*1e-6)
    #print 'amplitude = %.0f pA'% amplitude
    #J=np.append(np.array([curr,curr]),J)
    #J=np.append(J,np.array([curr,curr]))
    J=J*1e12
    #T=np.append(np.array([-1e6,0.]),T)
    #T=np.append(T,np.array([tau_off,1e6+tau_off]))
    #T=T*1e-6

    fig=plt.figure(figsize=(10,5),dpi=80)
    gs=gridspec.GridSpec(1, 2, width_ratios=[6e-1,1])

    ax0=plt.subplot(gs[0])
    plt.title('Trajectory')
    ax0.plot(X,Z,linewidth=1.,c='#0000ff')
    longer = ax0.scatter(X[bind1],Z[bind1],s=200,marker='h',c='#ff0000',linewidth=0.)
    shorter = ax0.scatter(X[bind2],Z[bind2],s=100,marker='h',c='#ffa500',linewidth=0.)
    start = ax0.scatter([X[0]],[Z[0]],s=200,marker='x',c='#00cc00',linewidth=2.)
    ax0.legend([start,longer,shorter],['Start','Long binding','Short binding'],scatterpoints=1,loc=(-.17,.3))
    ax=plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-35.,15.])
    ax.set_ylim([-30.,40.])
    ax.set_xticks([])
    ax.set_yticks([])
    plot_polygon(ax,polygon(rmem=40.))


    ax1=plt.subplot(gs[1])
    plt.title('Current signal')
    ax=plt.gca()
    if tau_off<1e3:
        t = np.linspace(0.,tau_off,5)
        fac=1.
        ax.set_xlabel('Time [$ns$]')
    elif tau_off<1e6 and tau_off>=1e3:
        t = np.linspace(0.,tau_off*1e-3,5)
        fac = 1e-3
        ax.set_xlabel(r'Time [$\mu s$]')
    else:
        t = np.linspace(0.,tau_off*1e-6,5)
        fac = 1e-6
        ax.set_xlabel('Time [$ms$]')
    T=T*fac
    ax1.plot(T,J,linewidth=2.,color='#0000ff')
    yt = np.linspace(580.,760,7)
    ax.set_ylabel(r'Current [$pA$]')
    ax.set_yticks(yt)
    ax.set_xticks(t)
    xfmt=FormatStrFormatter('%.1f')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim([-1e-2*tau_off*fac,(1.+1e-2)*tau_off*fac])


    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
    print 'savefig: %s'%figname
    plt.close("all")
