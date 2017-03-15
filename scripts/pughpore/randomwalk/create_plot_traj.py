from matplotlib.ticker import FormatStrFormatter
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
l0 =        up.l0
l1 =        up.l1
l2 =        up.l2
l3 =        up.l3
l4 =        up.l4
hpore =     up.hpore
hmem =      up.hmem
h2 =        up.h2
h1 =        up.h1
h4 =        up.h4

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
    figname = fieldsname+'_traj_'+'%f'%(tau_off*1e-6)+'_%i'%i+'_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["avgbind3"],params["P_bind1"],params["P_bind2"],params["P_bind3"])+str(params["z0"])
    J=J*1e12

    fig=plt.figure(figsize=(8,5),dpi=80)
    cmap = matplotlib.cm.get_cmap('viridis')
    color2=cmap(0.6)
    color1=cmap(0.8)
    color3=cmap(0.3)

    plt.plot(X,Z,linewidth=1.,c=cmap(0.0),alpha=.5)
    plt.plot([l3/2.-.5,l3/2.-.5],[-3.,11.],color=color2,linewidth=2.)
    plt.plot([-l3/2.+.5,-l3/2.+.5],[-3.,11.],color=color2,linewidth=2.)
    plt.plot([l1/2.,l1/2.],[17.,19.],color=color1,linewidth=2.)
    plt.plot([-l1/2.,-l1/2.],[17.,19.],color=color1,linewidth=2.)
    plt.plot([l3/2.,l3/2.,l2/2.,l2/2.],[-hpore/2.,hpore/2.-h2,hpore/2.-h2,14.],color=color1,linewidth=2.)
    plt.plot([-l3/2.,-l3/2.,-l2/2.,-l2/2.],[-hpore/2.,hpore/2.-h2,hpore/2.-h2,14.],color=color1,linewidth=2.)
    longer = plt.scatter(X[bind1],Z[bind1],s=200,marker='h',c=color2,linewidth=0.)
    shorter = plt.scatter(X[bind2],Z[bind2],s=100,marker='h',c=color1,linewidth=0.)
    start = plt.scatter([X[0]],[Z[0]],s=200,marker='x',c=color3,linewidth=2.)
    plt.legend([start,longer,shorter],['Start','Long binding','Short binding'],scatterpoints=1,loc=(.42,.15))
    ax=plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([20.,-55.])
    ax.set_ylim([-25.,40.])
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')
    plot_polygon(ax,polygon(rmem=60.))


    plt.axes([.55,.5,.2,.3])
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
    plt.plot(T,J,color='#000000')
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
