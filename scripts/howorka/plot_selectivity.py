import nanopores, numpy
import matplotlib.pyplot as plt

Qs = [-1.,-3.]
NAME = "howorka2D_selectivity_Q%.0f"
label = r"$Q = %.0f$"

         
def avg(J):
    n = len(J)
    J0 = list(numpy.array(J)[n*0.2:n*0.5])
    return sum(J0)/len(J0)
    
Javg = []

for Q in Qs:
    results, params = nanopores.load_stuff(NAME % Q)
    t = results["time"]
    J = results["current"]
    rel = results["release"]
    
    plt.figure(0)
    plt.semilogx(t, rel, "x-", label=label % Q)
    plt.xlabel("time [s]")
    plt.ylabel("% release")
    plt.title("reservoir size: %.0f nm" % (params["Ry"],))
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