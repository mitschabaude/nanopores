import os
import nanopores
import matplotlib.pyplot as plt
from mysolve import load_estimators as load

nanopores.add_params(est = "F")
F = "err ref" if est == "ref" else est
DIR = os.path.expanduser("~") + "/papers/pnps-numerics/figures/"

for dim in "2D", "3D":
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
    
    estimators[F].plot()
    estimators_cheap[F].plot(fig=False)
    estimators_unif[F].plot(fig=False, rate=-2./float(dim[0]))
    #plt.title(dim)
    plt.legend(bbox_to_anchor=(0.7, 1.1), loc="upper left", borderaxespad=0.,)
    
    plt.xlim(xmin=estimators[F].pairs[0][0]*0.75)
    # save
    fig = plt.gcf()
    fig.set_size_inches(4, 3.5)
    fig.savefig(DIR + "adap%scomp.eps" %dim, bbox_inches='tight')
    
#nanopores.showplots()
