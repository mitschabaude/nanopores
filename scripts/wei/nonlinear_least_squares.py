# (c) 2017 Gregor Mitscha-Baude
import numpy as np

def functional(xi, yi):
    def g(t):
        return np.sum((np.exp(-np.exp(xi)/t) - (1. - yi))**2)
    return g

def newton(f, x0=0., tol=1e-14):
    "solve f(x) = 0 with initial guess x = x0"
    from tangent import grad
    df = grad(f)
    x = x0
    while np.abs(f(x)) > tol:
        x -= f(x)/df(x)
        print "|f(x)|", np.abs(f(x))
    return x

def minimize(F, x0=0., tol=1e-14):
    "find local minimum of F near initial guess x=x0"
    # solve dF(x) = 0 with newton
    from tangent import grad
    return newton(grad(F), x0=x0, tol=tol)

def NLS(ti, yi, t0=0., tol=1e-14):
    "nonlinear least squares to fit exponential distribution with mean exp(x)"
    xi = np.log(ti)
    
    def f(x, xi, yi):
        return np.sum((1. - np.exp(-np.exp(xi - x)) - yi)**2)
    
    # minimize f by solving df(x) = 0 with newton method
    from tangent import grad
    df = grad(f)
    ddf = grad(df)
    x = np.log(t0)
    while np.abs(df(x, xi, yi)) > tol:
        x -= df(x, xi, yi)/ddf(x, xi, yi)
        #print "|f(x)|", np.abs(df(x, xi, yi))
    return np.exp(x)

def NLS_(F, xi, yi, x0=0., tol=1e-14):
    "nonlinear least squares to find parameter x so that F(xi, x) \approx yi"
    
    def f(x, xi, yi):
        return np.sum((F(xi, x) - yi)**2)
    
    # minimize f by solving df(x) = 0 with newton method
    from tangent import grad
    df = grad(f)
    ddf = grad(df)
    x = x0
    while np.abs(df(x, xi, yi)) > tol:
        x -= df(x, xi, yi)/ddf(x, xi, yi)
        print "|df(x)|", np.abs(df(x, xi, yi))
    return x

if __name__ == "__main__":
    x = 0.1234
    xi = np.linspace(-3, 3, 100)
    
    def F(xi, x):
        return 1. - np.exp(-np.exp(xi - x))
        
    yi = F(xi, x) + 1e-6*np.random.randn(100)
    
    print NLS_(F, xi, yi), x
