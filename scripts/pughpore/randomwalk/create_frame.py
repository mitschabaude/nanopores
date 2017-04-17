from matplotlib.ticker import FormatStrFormatter
import matplotlib
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
fieldsname='eventsnew_onlyone_4_'
params=dict(avgbind1=23e6,avgbind2=3e4,P_bind1=0.035,P_bind2=0*3e-1,z0=hpore/2.+0.)
i=311
fs=30
label_size = 30
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
data=f.get_fields(fieldsname,**params)
b1 =data["b1"]
b2 =data["b2"]
X = data["X"][i]
Y = data["Y"][i]
Z = data["Z"][i]
T = data["T"][i]
J = data["J"][i]
J=J.load()
T=T.load()
curr = 7.523849e-10
bind1 = np.where(T>1e6)
bind2 = np.intersect1d(np.where(T<=1e6),np.where(T>100.))
amplitude = curr-np.inner(J,T)/np.sum(T)
for k in range(1,T.shape[0]):
    T[k]=T[k]+T[k-1]
tau_off=T[-1]
J=J*1e12
len=X.shape[0]

def save_frame(r,X=X,Z=Z,T=T,J=J):
    figname = 'frames/'+'%08d'%r
    suffix='.png'
    figname+=suffix
    fig=plt.figure(figsize=(32,18),dpi=80)
    color2='#ff0000'
    color1='#ff9900'
    color3='#00ff00'

    #b1 = [[[l1/2.,17.],[l1/2.,19.]],[[l3/2.,-hpore/2.],[l3/2.,hpore/2.-h2],[l2/2.,hpore/2.-h2],[l2/2.,14.]]]
    for seq in b1:
        x= [p[0] for p in seq]
        xm=[-p[0] for p in seq]
        y= [p[1] for p in seq]
        plt.plot(x,y,color=color1,linewidth=2.)
        plt.plot(xm,y,color=color1,linewidth=2.)
    #b2 = [[[l3/2.-.5,-3.],[l3/2.-.5,11.]]]
    for seq in b2:
        x= [p[0] for p in seq]
        xm=[-p[0] for p in seq]
        y= [p[1] for p in seq]
        plt.plot(x,y,color=color2,linewidth=2.)
        plt.plot(xm,y,color=color2,linewidth=2.)
    plt.plot(X[:r],Z[:r],linewidth=1.,c='#0000ff')
    molecule = plt.Circle((X[r-1], Z[r-1]),2.0779,color='#33ff99')
    longer = plt.scatter(X[bind1],Z[bind1],s=200,marker='h',c=color2,linewidth=0.)
    shorter = plt.scatter(X[bind2],Z[bind2],s=100,marker='h',c=color1,linewidth=0.)
    start = plt.scatter([X[0]],[Z[0]],s=400,marker='x',c=color3,linewidth=2.)
    patches=[start]
    labels=['Start']
    plt.legend(patches,labels,scatterpoints=1,loc=(.42,.15),fontsize=fs)
    ax=plt.gca()
    ax.add_artist(molecule)
    ax.set_aspect('equal')
    ax.set_xlim([20.,-55.])
    ax.set_ylim([-25.,40.])
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')
    plot_polygon(ax,polygon(rmem=60.))


    plt.axes([.55,.5,.2,.3])
    plt.title('Current signal',fontsize=fs)
    ax=plt.gca()
    if tau_off<1e3:
        t = np.linspace(0.,tau_off,3)
        fac=1.
        ax.set_xlabel('time [$ns$]',fontsize=fs)
    elif tau_off<1e6 and tau_off>=1e3:
        t = np.linspace(0.,tau_off*1e-3,3)
        fac = 1e-3
        ax.set_xlabel(r'time [$\mu s$]',fontsize=fs)
    else:
        t = np.linspace(0.,tau_off*1e-6,3)
        fac = 1e-6
        ax.set_xlabel('time [$ms$]',fontsize=fs)
    T=T*fac
    plt.plot(T[:r],J[:r],color='#000000')
    yt = np.linspace(580.,760,4)
    ax.set_ylabel(r'A [$pA$]',fontsize=fs)
    ax.set_yticks(yt)
    ax.set_xticks(t)
    xfmt=FormatStrFormatter('%.1f')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_xlim([-4e-2*tau_off*fac,(1.+4e-2)*tau_off*fac])


    plt.tight_layout()
    #plt.show()
    plt.savefig(figname)
    print 'save: '+figname+' from %08d'%len
    plt.close('all')

for r in range(1,len):
    save_frame(r)
