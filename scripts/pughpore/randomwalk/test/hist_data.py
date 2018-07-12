# -*- coding: utf-8 -*-
from scipy.stats import gamma
from matplotlib import gridspec
import math
from math import exp
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


drop, th = f.get("events_pugh_experiment", "drop", "t")
len=th.load().shape[0]

la = 1./23.
a = 1.074

color2='green'
orange='lightgreen'
color3='red'
gray='#888888'
orange='#ff6600'
#orange=gray

log=True
#log=False

if log:
    plt.figure(figsize=(10,5),dpi=80)
    gs = gridspec.GridSpec(2,3,width_ratios=[4,2,1],height_ratios=[1,2.5])
else:
    plt.figure(figsize=(4,6),dpi=80)
    gs = gridspec.GridSpec(2,1,height_ratios=[1,2.5])
gs.update(wspace=0.,hspace=0.)

thleft=[x for x in th if x<2.]
thright=[x for x in th if x>=2.]
dropleft=[drop[i] for i in range(len) if th[i]<2.]
dropright=[drop[i] for i in range(len) if th[i]>=2.]



minperc=0.
maxperc=40.
if log:
    plt1=plt.subplot(gs[1,0])
    plt1.scatter(thright,dropright,s=8,color=orange)
    plt1.scatter(thleft,dropleft,s=8,color=gray)
    xfmt=FormatStrFormatter('%g')
    plt1.set_xlim([.2*min(th),max(th)*5.])
    plt1.set_ylim([minperc,maxperc])
    plt1.set_xscale('log')
    plt1.xaxis.set_major_formatter(xfmt)
    plt1.invert_yaxis()
    plt1.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
    plt1.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15,x=.76)

if log:
    plt2=plt.subplot(gs[1,1])
else:
    plt2=plt.subplot(gs[1,0])
plt2.scatter(thright,dropright,s=8,color=orange)
plt2.scatter(thleft,dropleft,s=8,color=gray)
plt2.invert_yaxis()
plt2.set_ylim([maxperc,minperc])
plt2.set_xlim([-2e-2*max(th),max(th)*(1.+2e-2)])
if log:
    plt2.axes.get_yaxis().set_visible(False)
    plt2.axes.get_xaxis().major.locator.set_params(nbins=6)
plt2.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
if not log:
    plt2.axes.get_xaxis().major.locator.set_params(nbins=7)
    plt2.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15)

alpha=.3

if log:
    plt3=plt.subplot(gs[1,2])
    n, bins, patches = plt3.hist(np.array(dropright),5,normed=1,orientation='horizontal',color=orange,alpha=alpha)
    plt3.invert_yaxis()
    plt3.set_xlim([0.,max(n)*1.2])
    plt3.set_ylim([maxperc,minperc])
    plt3.axes.get_xaxis().set_visible(False)
    plt3.axes.get_yaxis().set_visible(False)



if log:
    plt4=plt.subplot(gs[0,1])
else:
    plt4=plt.subplot(gs[0,0])
plt4.plot(np.linspace(1,100,100),np.array([gamma.pdf(x,a)*la**a*exp(x*(1-la)) for x in np.linspace(1,100,100)]),color=orange,linewidth=1.5)
n, bins, patches = plt4.hist(np.array(thright),15,normed=1,color=orange,alpha=alpha)
plt4.set_xlim([-2e-2*max(th),max(th)*(1.+2e-2)])
plt4.axes.get_xaxis().set_visible(False)
plt4.axes.get_yaxis().set_visible(False)


plt.tight_layout()
if log:
    pass
#    plt.savefig('hist_data1.pdf')
else:
    pass
#    plt.savefig('hist_data2.pdf')
plt.show()
