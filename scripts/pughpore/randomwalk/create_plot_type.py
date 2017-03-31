# -*- coding: utf-8 -*-
from matplotlib import gridspec
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
#fieldsname='events3_onlyone_5'
#params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=8.e-2,P_bind2=0*3e-1,z0=hpore/2.+0.)

def save_fig_type(params,fieldsname):
    cmap=matplotlib.cm.get_cmap('viridis')
    data=f.get_fields(fieldsname,**params)
    figname = fieldsname+'_%.1e_%.1e_%.1e_%.1e'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"])+str(params["z0"])
    t = data["t"]
    a = data["a"]
    ood = data["ood"]
    lendata=len(t)
    fac=1.
    if max(t)<1e-2:
        fac=1e3
        t = [x*1e3 for x in t]

    P_bind1=params["P_bind1"]
    P_bind2=params["P_bind2"]
    avgbind1=params["avgbind1"]*1e-6
    avgbind2=params["avgbind2"]*1e-6


    color2='green'
    color1='lightgreen'
    color3='red'

    plt.figure(figsize=(10,5),dpi=80)
    gs = gridspec.GridSpec(2,3,width_ratios=[4,2,1],height_ratios=[1,2.5])
    gs.update(wspace=0.,hspace=0.)



    minperc=0.
    maxperc=30.
    plt1=plt.subplot(gs[1,0])
    for k in range(lendata):
        if ood[k]==0:
            type1 = plt1.scatter([t[k]],[a[k]],color=color2,s=8)
        else:
            type0 = plt1.scatter([t[k]],[a[k]],color=color3,s=8)
    plt.legend([type1,type0],['successfull translocation','did not translocate'],scatterpoints=4,loc=(.4,1.02))
    xfmt=FormatStrFormatter('%g')
    plt1.set_xlim([.2*min(t),max(t)*5.])
    plt1.set_ylim([minperc,maxperc])
    plt1.set_xscale('log')
    plt1.xaxis.set_major_formatter(xfmt)
    plt1.invert_yaxis()
    plt1.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
    if fac==1.:
        if P_bind1!=0.:
            plt1.text(avgbind1*.5,27.,'Long binding',fontsize=9,horizontalalignment='center')
            k=1.0
            plt1.add_patch(matplotlib.patches.Rectangle((avgbind1*10**(-k*2),0.),avgbind1*(10**(k)-10**(-k)),maxperc,facecolor=cmap(.7),alpha=.15))
        if P_bind2!=0.:
            plt1.text(avgbind2*.5,27.,'Short binding',fontsize=9,horizontalalignment='center')
            k=1.0
            plt1.add_patch(matplotlib.patches.Rectangle((avgbind2*10**(-k),0.),avgbind2*(10**(k)-10**(-k)),maxperc,facecolor=cmap(.4),alpha=.15))
        plt1.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15,x=.76)
    else:
        plt1.set_xlabel(ur'$\tau_{off}$ [µs]',fontsize=15,x=.76)
    plt2=plt.subplot(gs[1,1])
    for k in range(lendata):
        if ood[k]==0:
            type1 = plt2.scatter([t[k]],[a[k]],color=color2,s=8)
        else:
            type0 = plt2.scatter([t[k]],[a[k]],color=color3,s=8)
    plt2.invert_yaxis()
    plt2.set_ylim([maxperc,minperc])
    plt2.set_xlim([-2e-2*max(t),max(t)*(1.+2e-2)])
    plt2.axes.get_yaxis().set_visible(False)
    plt2.axes.get_xaxis().major.locator.set_params(nbins=6)

    plt3=plt.subplot(gs[1,2])
    n, bins, patches = plt3.hist(np.array(a),15,normed=1,orientation='horizontal',color=color1,alpha=.5)
    plt3.invert_yaxis()
    plt3.set_xlim([0.,max(n)*1.2])
    plt3.set_ylim([maxperc,minperc])
    plt3.axes.get_xaxis().set_visible(False)
    plt3.axes.get_yaxis().set_visible(False)



    plt4=plt.subplot(gs[0,1])
    n, bins, patches = plt4.hist(np.array(t),20,normed=1,color=color1,alpha=.5)
    plt4.set_xlim([-2e-2*max(t),max(t)*(1.+2e-2)])
    plt4.axes.get_xaxis().set_visible(False)
    plt4.axes.get_yaxis().set_visible(False)


    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/',pdf=True)
    print 'savefig:'
    print figname
