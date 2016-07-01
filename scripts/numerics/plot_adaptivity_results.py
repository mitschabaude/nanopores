import nanopores
import mysolve

estimators = mysolve.load_estimators("adap2D")
estimators_cheap = mysolve.load_estimators("adap2Dcheap")
estimators_unif = mysolve.load_estimators("adap2Duniform")

F = "F"

estimators[F].name = "adaptive"
estimators_cheap[F].name = "adaptive (cheap estimator)"
estimators_unif[F].name = "uniform"

estimators[F].plot()
estimators_cheap[F].plot(fig=False)
estimators_unif[F].plot(fig=False, rate=-0.75)

nanopores.showplots()
