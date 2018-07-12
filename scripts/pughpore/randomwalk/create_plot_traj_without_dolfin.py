from matplotlib.ticker import FormatStrFormatter
import matplotlib
#import nanopores as nano
#import nanopores.geometries.pughpore as pughpore
#from nanopores.models.pughpore import polygon
#from nanopores.models.pughpoints import plot_polygon
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#import nanopores.tools.fields as f
#HOME = os.path.expanduser("~")
#PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
#FIGDIR = os.path.join(PAPERDIR, "figures", "")
#DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
#f.set_dir(DATADIR)
###############################
# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

default = dict(
    R = 20.,
    H = 70.,
    l0 = 18., #22.5,
    l1 = 14., #17.5,
    l2 = 10., #12.5,
    l3 = 6., #7.5,
    l4 = 14., #17.5,
    hpore = 46,
    hmem = 2.2,

    h2 = 46.-35., # 11.
    h1 = 46.-35.-2.5, # 8.5
    h4 = 10.,

    hnear = 12,

    rMolecule = 2.0779, # molecular radius of protein trypsin
    x0 = [0., 0., 0.],
    lcMolecule = 0.4, # relative to global mesh size
    center_at_x0 = False,
    center_z_at_x0 = False,
)

class Params(dict):
    "for writing params.Qmol instead of params['Qmol']"
    def __getattr__(self, key):
        return self[key]
    def __or__(self, other):
        new = Params(self)
        new.update(other)
        return new

def plot_polygon(ax, polygon):
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=1.,
                    edgecolor="black")
    polygon = np.array(polygon)
    polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    patch = patches.Polygon(polygon, **settings)
    patchm = patches.Polygon(polygon_m, **settings)
    #patch.set_zorder(10)
    #patchm.set_zorder(10)
    ax.add_patch(patch)
    ax.add_patch(patchm)

def polygon(rmem=20., **params):
    "polygon of pore + membrane for plotting"
    #setup = SetupNoGeo(**params)
    #params = nano.Params(pughpore.params) | setup.geop
    params = Params(default, **params)

    r = [0.5*params.l3, 0.5*params.l2, 0.5*params.l1, 0.5*params.l0,
         0.5*params.l4, rmem]
    ztop = params.hpore/2.
    zbot = -ztop
    z = [zbot, ztop - params.h2, ztop - params.h1, ztop, zbot + params.h4,
         zbot + params.hmem]
    # indices: [(0,0), (0,1), (1,1), (1,2), ..., (5,5), (5,0)]
    return [(r[i / 2 % 6], z[(i+1) / 2 % 6]) for i in range(12)]
###############################

up = Params(default, k=3)
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
#fieldsname='events_onlyone_2'
#params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=8.e-2,P_bind2=0*3e-1,z0=hpore/2.+0.)
#i=15
#showtraj = True

#def save_fig_traj(params,fieldsname,i,showtraj):
def assertdir(DIR):
    if not os.path.exists(DIR):
        os.makedirs(DIR)

def savefigs(name="fig", DIR="/tmp/", size=None, pdf=False):
    if not DIR.endswith("/"): DIR = DIR + "/"
    assertdir(DIR)
    if len(plt.get_fignums()) == 1:
        fig = plt.figure(plt.get_fignums()[0])
        if size is not None:
            fig.set_size_inches(size)
        if not pdf: suffix='.eps'
        else: suffix='.pdf'
        fig.savefig(DIR + name + suffix, bbox_inches="tight")
        return
    for num in plt.get_fignums():
        fig = plt.figure(num)
        label = fig.get_label()
        label = str(num) if label=="" else label
        if size is not None:
            fig.set_size_inches(size)
        if not pdf: suffix='.eps'
        else: suffix='.pdf'
        fig.savefig(DIR + name + "_" + label + suffix, bbox_inches="tight")
def save_fig_traj():
    showtraj = False
#    data=f.get_fields(fieldsname,**params)
    b1 = [[[l3/2.,-hpore/2.],[l3/2.,hpore/2.-h2],[l2/2.,hpore/2.-h2],[l2/2.,hpore/2.-h1],[l1/2.,hpore/2.-h1],[l1/2.,hpore/2.]]]
#    b2 = [[[l3/2.-.5,11.],[l3/2.-.5,-3.]]]
#    b2 = [[[l3/2.,-11.],[l3/2.,3.]]]
   # b1 =data["b1"]
   # b2 =data["b2"]
    if showtraj:
        X = data["X"][0]
        Y = data["Y"][0]
        Z = data["Z"][0]
        T = data["T"][0]
        J = data["J"][0]
        J=J.load()
        T=T.load()
        curr = 7.523849e-10
        bind1 = np.where(T>1e6)[0]
        bind2 = np.intersect1d(np.where(T<=1e6)[0],np.where(T>100.)[0])
        amplitude = curr-np.inner(J,T)/np.sum(T)
        for k in range(1,T.shape[0]):
            T[k]=T[k]+T[k-1]
        tau_off=T[-1]
        J=J*1e12
        figname = fieldsname+'_traj_'+'%.8f'%(tau_off*1e-6)+'_%04d'%i+'_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"])+str(params["z0"])
    else:
#        figname = fieldsname+'_bindzones'+'_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"])+str(params["z0"])
        figname = 'bindzones_both.eps'

    if showtraj:
        fig=plt.figure(figsize=(8,5),dpi=80)
    else:
        fig=plt.figure(figsize=(3,5),dpi=80)
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
    b2 = [[[l3/2.-.5,-4.],[l3/2.-.5,4.]]]
    for seq in b2:
        x= [p[0] for p in seq]
        xm=[-p[0] for p in seq]
        y= [p[1] for p in seq]
        plt.plot(x,y,color=color2,linewidth=2.)
        plt.plot(xm,y,color=color2,linewidth=2.)
    if showtraj:
        plt.plot(X,Z,linewidth=1.,c='#0000ff')
        longer = plt.scatter(X[bind1],Z[bind1],s=200,marker='h',c=color2,linewidth=0.)
        shorter = plt.scatter(X[bind2],Z[bind2],s=100,marker='h',c=color1,linewidth=0.)
        start = plt.scatter([X[0]],[Z[0]],s=200,marker='x',c=color3,linewidth=2.)
        patches=[start]
        labels=['Start']
    if showtraj and len(bind1)>0:
        patches=patches+[longer]
        labels+=['Longer bindings']
    if showtraj and len(bind2)>0:
        patches=patches+[shorter]
        labels+=['Shorter bindings']
    if showtraj:
        plt.legend(patches,labels,scatterpoints=1,loc=(.42,.15))
    ax=plt.gca()
    ax.set_aspect('equal')
    if showtraj:
        ax.set_xlim([20.,-55.])
        ax.set_ylim([-25.,40.])
    else:
        ax.set_xlim([20.,-20.])
        ax.set_ylim([-25.,40.])
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')
    plot_polygon(ax,polygon(rmem=60.))


    if showtraj:
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
    savefigs(name=figname,DIR='/home/bstadlbau/plots/')
#    plt.show()
#    print 'savefig: %s'%figname
    plt.close("all")

save_fig_traj()
