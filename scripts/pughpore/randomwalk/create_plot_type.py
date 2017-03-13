import math
import matplotlib
matplotlib.use("Agg")
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import os
import sys
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

number=False

geop = nano.Params(pughpore.params)
hpore=geop.hpore
#params=dict(avgbind=8.7e6,P_bind=5.e-3,z0=hpore/2.+5.)
def save_fig(params,fieldsname):
    figname = fieldsname+'_%.1e_%.1e_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["avgbind3"],params["P_bind1"],params["P_bind2"],params["P_bind3"])+str(params["z0"])
    data=f.get_fields(fieldsname,**params)
    t1 = data["t1"]
    a1 = data["a1"]
    t2 = data["t2"]
    a2 = data["a2"]
    t0 = data["t0"]
    a0 = data["a0"]
    print len(t1)
    print len(t2)
    lendata=len(t1)+len(t2)
    for i in range(len(t1)):
            if math.isnan(a1[i]): print 'a1: NaN %i'%i
    plt.figure(figsize=(6,4),dpi=80)
    plt.plot([3e-6,.7],[0.,0.],linewidth=2,color='lightgreen')
    plt.plot([1.,1e2],[0.,0.],linewidth=2,color='green')
    ax=plt.gca()
    plt.scatter(t1,a1,color='lightgreen',s=8)
    plt.scatter(t2,a2,color='green',s=8)
    plt.scatter(t0,a0,marker='o',s=50,facecolors='none',edgecolors='#ff0000')
    #ax.text(t1[k],a1[k],'%i'%k,fontsize=9)
    xfmt=FormatStrFormatter('%g')
    ax.set_xlim([1e-6,500.])
    ax.set_ylim([-2.,25.0])
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(xfmt)
    ax.invert_yaxis()
    ax.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15)
    ax.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
    ax.text(.001,-0.03,'I',fontsize=15)
    ax.text(5.,-0.03,'II',fontsize=15)
    #ax.text(1.,5.,lendata,fontsize=25)
    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
    print 'savefig:'
    print figname
