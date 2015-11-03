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
from nanopores import 

def plotfit(x, y, k=1):
    x = numpy.array(x)
    y = numpy.array(y)
    p = numpy.polyfit(x,y,k)
    C = p[0] # leading constant
    y_ = numpy.polyval(p, x)
    #plot(x, y_, label="degree %s fit" %k)
    plot(x, C*x**k, label="cx^%s fit (c=%s)" %(k,C))

def plotMFPT(x, y, title=""):
    plot(x,y, label=title)
    plotfit(x, y, 2)
    legend(loc="upper left")
    xlabel("domain scaling")
    ylabel("time [s]")
    savefig("np/data/ahem/%s.eps" %title, bbox_inches='tight')
    show()
    
x = numpy.arange(1,6)
