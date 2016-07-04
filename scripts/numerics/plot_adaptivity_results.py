import nanopores
from mysolve import load_estimators as load

nanopores.add_params(est = "F")
F = "err ref" if est == "ref" else est

estimators = load("adap2D")
estimators_cheap = load("adap2Dcheap")
estimators_unif = load("adap2Duniform")

estimators[F].name = "adaptive"
estimators_cheap[F].name = "adaptive (cheap estimator)"
estimators_unif[F].name = "uniform"

estimators[F].plot()
estimators_cheap[F].plot(fig=False)
estimators_unif[F].plot(fig=False, rate=-0.75)

nanopores.showplots()
