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
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

number=False

geop = nano.Params(pughpore.params)
hpore=geop.hpore
#fieldsname='events_onlyone_1'
#params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=8.e-2,P_bind2=0*3e-1,z0=hpore/2.+0.)

def save_fig_type(fieldsname,params):
    data=f.get_fields(fieldsname,**params)
    figname = fieldsname+'_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"])+str(params["z0"])
    t = data["t"]
    a = data["a"]
    ood = data["ood"]
    lendata=len(t)
    fac=1.
    if max(t)<1e-2:
        fac=1e6
        t = [x*1e6 for x in t]


    color2='green'
    color1='lightgreen'
    color3='red'

    plt.figure(figsize=(10,5),dpi=80)
    gs = gridspec.GridSpec(2,3,width_ratios=[4,2,1],height_ratios=[1,2.5])
    gs.update(wspace=0.,hspace=0.)



    plt1=plt.subplot(gs[1,0])
    for k in range(lendata):
        if t[k]<1.*fac:
            type1 = plt1.scatter([t[k]],[a[k]],color=color1,s=8)
        elif t[k]>=1.*fac:
            type2 = plt1.scatter([t[k]],[a[k]],color=color2,s=8)
        if ood[k]==1:
            type0 = plt1.scatter([t[k]],[a[k]],marker='o',s=50,facecolors='none',edgecolors=color3)
    try: plt.legend([type1,type2,type0],[r'$\tau_{off}$ shorter than 1ms',r'$\tau_{off}$ longer than 1ms','did not translocate'],scatterpoints=4,loc=(.4,1.02))
    except: plt.legend([type1,type0],[r'$\tau_{off}$ shorter than 1ms','did not translocate'],scatterpoints=4,loc=(.4,1.02))
    xfmt=FormatStrFormatter('%g')
    plt1.set_xlim([.2*min(t),max(t)*5.])
    plt1.set_ylim([-2.,25.0])
    plt1.set_xscale('log')
    plt1.xaxis.set_major_formatter(xfmt)
    plt1.invert_yaxis()
    if fac is not None:
        plt1.set_xlabel(r'$\tau_{off}$ [ns]',fontsize=15,x=.76)
    else:
        plt1.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15,x=.76)
    plt1.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
    if fac==1.:
        plt1.plot([3e-6,.7],[0.,0.],linewidth=3,color=color1)
        plt1.plot([1.,1e2],[0.,0.],linewidth=3,color=color2)
        plt1.text(.001,-0.03,'I',fontsize=15)
        plt1.text(5.,-0.03,'II',fontsize=15)
    plt2=plt.subplot(gs[1,1])
    for k in range(lendata):
        if t[k]<1.*fac:
            plt2.scatter([t[k]],[a[k]],color=color1,s=8)
        elif t[k]>=1.*fac:
            plt2.scatter([t[k]],[a[k]],color=color2,s=8)
    plt2.invert_yaxis()
    plt2.set_ylim([25.0,-2.])
    plt2.set_xlim([-2e-2*max(t),max(t)*(1.+2e-2)])
    plt2.axes.get_yaxis().set_visible(False)
    plt2.axes.get_xaxis().major.locator.set_params(nbins=6)

    plt3=plt.subplot(gs[1,2])
    n, bins, patches = plt3.hist(np.array(a),20,normed=1,orientation='horizontal',color=color1,alpha=.5)
    plt3.invert_yaxis()
    plt3.set_xlim([0.,max(n)*1.2])
    plt3.set_ylim([25.0,-2.])
    plt3.axes.get_xaxis().set_visible(False)
    plt3.axes.get_yaxis().set_visible(False)



    plt4=plt.subplot(gs[0,1])
    n, bins, patches = plt4.hist(np.array(t),20,normed=1,color=color1,alpha=.5)
    plt4.set_xlim([-2e-2*max(t),max(t)*(1.+2e-2)])
    plt4.axes.get_xaxis().set_visible(False)
    plt4.axes.get_yaxis().set_visible(False)


    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
    print 'savefig:'
    print figname
