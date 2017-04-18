# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt

def Dt_plane(h, r):
    x = r/h
    return 1. - 9./16.*x + 1./8.*x**3 - 45./256.*x**4 - 1/16.*x**5

def Dn_plane(l, r, N=20):
    alpha = np.arccosh(l/r)
    sinh = np.sinh
    s = 0.
    for n in range(1, N):
        n = float(n)
        K = n*(n+1)/(2*n-1)/(2*n+3)
        s += K*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2) - 1)
    return 1./((4./3.)*sinh(alpha)*s)
    
def Dn_sphere(l, r, R): # happel brenner
    l = l + R
    return 1. / (1. + 9./4.*r*R/l**2 + 3./4./l**4*(-2.*r**3*R + 27./4.*r**2*R**2 + 3.*r*R**3))
    
def Dt_sphere(l, r, R): # happel brenner
    l = l + R
    return 1. / (1. + 9./16.*r*R/l**2 + 3./8./l**4*(r**3*R + 27./32.*r**2*R**2 + 3.*r*R**3))

r = 1.
# TODO: for large R, approximation breaks down completely..
# better use plane in that case!
R = 1. #2.0779 / 0.11
x = np.logspace(1e-10, 2, 100)

plt.plot(x, Dn_sphere(x, r, R), label="D sphere normal")
plt.plot(x, Dt_sphere(x, r, R), label="D sphere parallel")
plt.plot(x, Dn_plane(x, r), label="D plane parallel")
plt.plot(x, Dt_plane(x, r), label="D plane parallel")
plt.xscale("log")
ticks = [1, 2, 3, 5, 10, 20, 30, 50, 100]
plt.xticks(ticks, ticks)
plt.grid()
plt.legend(loc="lower right")