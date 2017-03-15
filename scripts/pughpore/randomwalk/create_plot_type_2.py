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
#fieldsname='test_new_bind_5'
#params=dict(avgbind1=2e7,avgbind2=3e4,avgbind3=2e5,P_bind1=8.e-2,P_bind2=3e-1,P_bind3=5e-2,z0=hpore/2.+0.)
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



    cmap = matplotlib.cm.get_cmap('viridis')
    color2=cmap(0.6)
    color1=cmap(0.8)
    color3=cmap(0.3)

    plt.figure(figsize=(10,5),dpi=80)
    gs = gridspec.GridSpec(2,3,width_ratios=[4,2,1],height_ratios=[1,2.5])
    gs.update(wspace=0.,hspace=0.)



    plt1=plt.subplot(gs[1,0])
    plt1.plot([3e-6,.7],[0.,0.],linewidth=3,color=color1)
    plt1.plot([1.,1e2],[0.,0.],linewidth=3,color=color2)
    type1 = plt1.scatter(t1,a1,color=color1,s=8)
    type2 = plt1.scatter(t2,a2,color=color2,s=8)
    type0 = plt1.scatter(t0,a0,marker='o',s=50,facecolors='none',edgecolors=color3)
    plt.legend([type1,type2,type0],[r'$\tau_{off}$ shorter than 1ms',r'$\tau_{off}$ longer than 1ms','did not translocate'],scatterpoints=4,loc=(.4,1.02))
    #ax1.text(t1[k],a1[k],'%i'%k,fontsize=9)
    xfmt=FormatStrFormatter('%g')
    plt1.set_xlim([1e-6,500.])
    plt1.set_ylim([-2.,25.0])
    plt1.set_xscale('log')
    plt1.xaxis.set_major_formatter(xfmt)
    plt1.invert_yaxis()
    ticks=[10**(r) for r in [-5,-4,-3,-2,-1,0.,1.,2.]]
    #plt1.get_xaxis().set_tick_params(direction='out',pad=-20)
    plt1.set_xticks(ticks)
    plt1.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15,x=.76)
    plt1.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
    plt1.text(.001,-0.03,'I',fontsize=15)
    plt1.text(5.,-0.03,'II',fontsize=15)
    #plt1.text(1.,5.,lendata,fontsize=25)

    plt2=plt.subplot(gs[1,1])
    plt2.scatter(t1,a1,color=color1,s=8)
    plt2.scatter(t2,a2,color=color2,s=8)
    plt2.invert_yaxis()
    plt2.set_ylim([25.0,-2.])
    plt2.set_xlim([-2,max(max(t1),max(t2))])
    plt2.axes.get_yaxis().set_visible(False)
    ticks_2=[20.,40.,60.,80.,100.]
    plt2.set_xticks(ticks_2)
    #plt2.get_xaxis().set_tick_params(direction='out',pad=-20)
    #plt2.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15)

    plt3=plt.subplot(gs[1,2])
    n, bins, patches = plt3.hist(np.array(a1+a2),35,normed=1,orientation='horizontal',color=color1,alpha=.5)
    plt3.invert_yaxis()
    plt3.set_xlim([0.,max(n)*1.2])
    plt3.set_ylim([25.0,-2.])
    plt3.axes.get_xaxis().set_visible(False)
    plt3.axes.get_yaxis().set_visible(False)



    plt4=plt.subplot(gs[0,1])
    n, bins, patches = plt4.hist(np.array(t1+t2),20,normed=1,color=color1,alpha=.5)
    #plt4.invert_yaxis()
    #plt4.set_xlim([0.,max(n)*1.2])
    #plt4.set_ylim([25.0,-2.])
    plt4.axes.get_xaxis().set_visible(False)
    plt4.axes.get_yaxis().set_visible(False)


    plt.tight_layout()
    nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
    print 'savefig:'
    print figname
