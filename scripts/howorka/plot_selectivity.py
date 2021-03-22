import numpy
import nanopores
import matplotlib.pyplot as plt
from .selectivity import selectivity, default

Qs = [-1.,-3.]
label = r"$Q = %.0f$"

default = nanopores.user_params(**default)

def avg(J):
    n = len(J)
    J0 = list(numpy.array(J)[n*0.2:n*0.5])
    return sum(J0)/len(J0)

Javg = []

for Q in Qs:
    params = dict(default)
    params["Qmol"] = Q
    results = selectivity(**params)
    t = results["time"]
    J = results["current"]
    rel = results["release"]
    params0 = results["params"]

    plt.figure(0)
    plt.semilogx(t, rel, "x-", label=label % Q)
    plt.xlabel("time [s]")
    plt.ylabel("% release")
    plt.title("reservoir size: %.0f nm" % (params0["Ry"],))
    plt.ylim(ymin=0.)

    plt.figure(1)
    plt.semilogx(t, J, "x-", label=label % Q)
    plt.xlabel("time [s]")
    plt.ylabel("current through pore [1/ms]")
    plt.ylim(ymin=0.)

    # determine average current
    Javg.append(avg(J))
    plt.plot(t, [Javg[-1]]*len(t), "k--", label="%.1f" % Javg[-1])

plt.title("selectivity: %.1f / %.1f = %.1f" % (
    max(Javg), min(Javg), max(Javg)/min(Javg)))


plt.figure(0)
plt.legend(loc="best")

plt.figure(1)
#plt.plot(t, [0.]*len(t), "k-")
plt.legend(loc="best")

plt.show()