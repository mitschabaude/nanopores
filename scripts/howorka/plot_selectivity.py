import nanopores
import matplotlib.pyplot as plt

Qs = [-1.,-3.]
NAME = "howorka2D_selectivity_Q%.0f"
label = r"$Q = %.0f$"

for Q in Qs:
    results, params = nanopores.load_stuff(NAME % Q)
    t = results["time"]
    J = results["current"]
    rel = results["release"]
    
    plt.figure(0)
    plt.plot(t, rel, "x-", label=label % Q)
    plt.xlabel("time [s]")
    plt.ylabel("% release")
    plt.title("reservoir size: %.0f nm" % (params["Ry"],))
    
    plt.figure(1)
    plt.plot(t, J, "x-", label=label % Q)    
    plt.xlabel("time [s]")
    plt.ylabel("current through pore [1/ms]")
    plt.title("reservoir size: %.0f nm" % (params["Ry"],))
    
plt.figure(0)
plt.legend(loc="best")

plt.figure(1)
plt.legend(loc="best")

plt.show()