# (c) 2017 Gregor Mitscha-Baude
import numpy as np

def solve(C, f, f1, x0=1., n=20):
    "solve f(x) == C"
    x = x0 # initial value
    print("Newton iteration:")
    for i in range(n):
        dx = -(f(x) - C)/f1(x)
        x = x + dx
        print(i, "Residual", f(x) - C, "Value", x)
    print()
    return x

def poisson_from_positiveK(mean):
    # solve x/(1 - exp(-x)) == mean
    def f(x):
        return x/(1. - np.exp(-x))
    def f1(x):
        return (np.expm1(x) - x)/(2.*np.cosh(x) - 2.)

    x = solve(mean, f, f1, mean, n=10)
    return x
