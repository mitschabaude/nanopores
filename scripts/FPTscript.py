result = {'tmolavg': [1.3859278792745953e-05,
  5.620329856283472e-05,
  0.00012694031170803607,
  0.0002257161225403611,
  0.00035293731252042414],
 'tmolmax': [1.5096023956506486e-05,
  6.11560324796202e-05,
  0.0001380757058487864,
  0.00024564587998346154,
  0.00038420164505047227],
 'tmolmin': [8.922562751966426e-06,
  3.635761294185034e-05,
  8.316546349049934e-05,
  0.00014571394626902007,
  0.00023023670497108691],
 'tporeavg': [1.5574475812731874e-06,
  6.311337654507861e-06,
  1.4305761006503032e-05,
  2.5421526895897307e-05,
  3.973465937932169e-05],
 'tporemax': [7.115944951898106e-06,
  2.891679075712449e-05,
  6.539653615875759e-05,
  0.00011615400759909166,
  0.00018149779368301853],
 'tporemin': [2.414590480759685e-07,
  1.004871012098735e-06,
  2.2949392849512985e-06,
  4.114123090460408e-06,
  6.455691024101006e-06]}

from numpy import *
from matplotlib.pyplot import *

def fit(x, y, k): return polyval(polyfit(x, y, k), x)
def err(y, y_): return linalg.norm(y-y_)/linalg.norm(y)

def plotfit(x, y, k=1):
    x = array(x)
    y = array(y)
    p = polyfit(x,y,k)
    C = p[0] # leading constant
    y_ = polyval(p, x)
    #plot(x, y_, label="degree %s fit" %k)
    plot(x, C*x**k, label="cx^%s fit (c=%s)" %(k,C))
    
def plot_inversefit(x, y, k=1, plot=plot):
    x = array(x)
    y = array(y)
    p = polyfit(1./x, y, k)
    #print p
    y_ = polyval(p, 1./x)
    c = p[-1] # constant term
    plot(x, y_, "+-", label="degree %s fit in 1/x" %k)
    #plot(x, C*x**k, label="cx^%s fit (c=%s)" %(k,C))
    plot(x, [c for i in x], "--", label="c = %.5f (from the fit)" %c)
    return y_

def plotMFPT(x, y, title=""):
    plot(x,y, label=title)
    plotfit(x, y, 2)
    legend(loc="upper left")
    xlabel("domain scaling")
    ylabel("time [s]")
    savefig("np/data/ahem/%s.eps" %title, bbox_inches='tight')
    show()
    
def plotSER(x, y, k=1, title=""):
    plot = semilogx
    L = 15. # length scale [nm]
    plot(x*L, y, "s-", label=title)
    y_ = plot_inversefit(x*L, y, k=k, plot=plot)
    #print err(y, y_)
    legend(loc="lower right")
    xlabel("upper reservoir size [nm]")
    ylabel("probability [%]")
    #savefig("np/data/ahem/%s.eps" %title, bbox_inches='tight')
    #show()
    return err(y, y_)
    
x = array([  1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 20., 40.])
y = array([ 1.43748472,  2.85920221,  3.34647272,  3.67373476,  3.84163926,    3.93641666,  4.02680935,  4.09577207,  4.13149867,  4.16848892, 4.267684, 4.318430])
#x = array([ 2., 5., 10., 20., 40.])
#y = array([ 2.85920221,  3.84163926, 4.16848892, 4.267684, 4.318430])
l = []
ran = list(range(11))
for k in ran: l.append(plotSER(x, y, k, title="Successful exit rate from molecule"))
figure()
semilogy(ran, l, "s-")
show()
        

    
