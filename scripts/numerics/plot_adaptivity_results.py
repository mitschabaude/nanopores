import os
import nanopores
import matplotlib.pyplot as plt
from mysolve import load_estimators as load

nanopores.add_params(est = "F")
F = "err ref" if est == "ref" else est
DIR = os.path.expanduser("~") + "/papers/pnps-numerics/figures/"

for dim in "3D", "2D":
    estimators = load("adap%s" %dim)
    estimators_cheap = load("adap%scheap" %dim)
    estimators_unif = load("adap%suniform" %dim)
    
    for est in estimators, estimators_cheap, estimators_unif:
        est["Fabs"] = nanopores.Estimator("Fabs")
        est["Fabs"].pairs = [(N, abs(drag) + abs(el))
            for (N, drag), (N1, el) in zip(est["Fdrag"].pairs, est["Fel"].pairs)]
    
    estimators[F].name = "adaptive"
    estimators_cheap[F].name = "adaptive (cheap)"
    estimators_unif[F].name = "uniform"
    
    rate = -2./float(dim[0])
    estimators[F].plot()
    estimators_cheap[F].plot(fig=False)
    estimators_unif[F].plot(fig=False, rate=rate)
    
    ax = plt.gca()
    line = ax.lines[3]
    label = line.get_label()
    label = label.replace(r"0.67", r"2/3")
    line.set_label(label)
    
    #plt.title(dim)
    plt.legend(bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
    
    plt.xlim(xmin=estimators_unif[F].pairs[0][0]*0.75, xmax=estimators_unif[F].pairs[-1][0]*1.33)
    plt.ylim(ymin=estimators[F].pairs[-1][1]*0.75, ymax=max([x[1] for x in estimators[F].pairs])*1.33)
    # save
    fig = plt.gcf()
    fig.set_size_inches(4, 3.5)
    #fig.savefig(DIR + "adap%scomp.eps" %dim, bbox_inches='tight')
    
    estimators["err ref"].plot(rate=rate)
    estimators["rep"].plot(fig=False)
    estimators["err"].plot(fig=False)

    
nanopores.showplots()
