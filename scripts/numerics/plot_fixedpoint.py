import os
import matplotlib.pyplot as plt
from mysolve import load_estimators as load

DIR = os.path.expanduser("~") + "/papers/pnps-numerics/figures/"

fixedpoint = load("fixedpoint")
hybrid = load("hybrid")
newton = load("newton")

fixedpoint["err hybrid i"].newtonplot()
hybrid["err hybrid i"].newtonplot(fig=False, style="rs-")
newton["err newton i"].newtonplot(fig=False, style="gs-")
plt.ylim(ymax=2.)
plt.legend(bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
fig = plt.gcf()
fig.set_size_inches(4, 3.5)
fig.savefig(DIR + "fixedhybridnewton_i.eps", bbox_inches='tight')

fixedpoint["err hybrid time"].newtonplot()
hybrid["err hybrid time"].newtonplot(fig=False, style="rs-")
newton["err newton time"].newtonplot(fig=False, style="gs-")
plt.ylim(ymax=2.)
plt.legend(bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
plt.xlabel("time [s]")
plt.xscale("log")
fig = plt.gcf()
fig.set_size_inches(4, 3.5)
fig.savefig(DIR + "fixedhybridnewton_t.eps", bbox_inches='tight')


#plt.xlim(xmin=estimators_unif[F].pairs[0][0]*0.75, xmax=estimators_unif[F].pairs[-1][0]*1.33)
#plt.ylim(ymin=estimators[F].pairs[-1][1]*0.75, ymax=max([x[1] for x in estimators[F].pairs])*1.33)

