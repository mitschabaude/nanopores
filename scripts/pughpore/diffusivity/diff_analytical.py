# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#from math import sinh, acosh

def Dt_plane(h, r):
    x = r/h
    return 1. - 9./16.*x + 1./8.*x**3 - 45./256.*x**4 - 1/16.*x**5

sinh = np.sinh
acosh = np.arccosh
def Dn_plane(l, r, N=10):
    alpha = acosh(l/r)
    s = 0.
    for n in range(1, N):
        n = float(n)
        K = n*(n+1)/(2*n-1)/(2*n+3)
        s += K*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2) - 1)
    return 1./((4./3.)*sinh(alpha)*s)


def Dt_equal_planes(R, r):
    x = r/R
    return 1. - 1.004*x + 0.418*x**3 + 0.21*x**4 - 0.169*x**5

def Dt_3R_planes(R, r):
    x = r/R
    return 1. - 0.6526*x + 0.1475*x**3 - 0.131*x**4 - 0.0644*x**5

# TODO: general two planes and sphere
# infinite cylinder

# tabulated values from Happel, Brenner
_x = [.00, .01, .03, .05, .1, .2, .3, .35, .37, .39, .4, .41, .43, .45, .5, .6,
     .7, .8, .9]
_y = [2.10444, 2.10436, 2.10381, 2.10270, 2.09763, 2.07942, 2.05691, 2.04805,
     2.04567, 2.04426, 2.04401, 2.04407, 2.04530, 2.04825, 2.06566, 2.16965,
     2.45963, 3.2316, 5.905]
x = _x + [1.]
y1 = [(1 - t)*ft for t, ft in zip(_x, _y)] + [9./16.]
f1 = interp1d(x, y1)
def f(t):
    return f1(t) / (1. - t)

#eta = 1e-3
#abs_friction = 6.*np.pi*eta*r
#rel_friction = 1. / (1. - f(t)*r/R)
def Dt_cyl(R, r, ecc=0.): # Happel, Brenner 1965 (book), valid for r/R small
    return 1. - f(ecc)*r/R

def Dt_cyl_alt(R, r, ecc=0.):
    a = f(ecc)*r/R
    return 1./ (1. + a + a**2)

def Dt_cyl_ax_large_dist(R, r):
    t = r/R
    return 1. - 2.10444*t + 2.08877*t**3 -  0.94813*t**5 - 1.372*t**6 + 3.87*t**8 - 4.19*t**10

def Dt_cyl_ax(R, r): # Bungay, Brenner 1973, Motion of a sphere in a fluid-filled tube
    t = r/R
    #print t
    tm = 1. - t
    brac = 1. - 73./60.*tm + 77.293/50.400*tm**2
    pol = -22.5083 - 5.6117*t - 0.3363*t**2 - 1.216*t**3 + 1.647*t**4
    fric = 9./4.*np.pi**2*np.sqrt(2.)*tm**(-2.5)*brac + pol
    return 6.*np.pi/fric

def D_sphere(R, r): # Happel, Brenner 1965 (book), valid for r/R small
    return 1. - 9./4.*r/R

r = 1.
eps = 1e-10
Rmax = 100.
x = np.logspace(eps, np.log10(Rmax), 100)

def plot_dashed(x, y, thresh=2., **kwargs):
    label = kwargs.pop("label", None)
    x0 = x[x < thresh]
    y0 = y[x < thresh]
    plt.plot(x0, y0, "--", **kwargs)
    x0 = x[x >= thresh]
    y0 = y[x >= thresh]
    plt.plot(x0, y0, "-", label=label, **kwargs)

blue = ["#000099", "#0000ff", "#9999ff"][::-1]
plt.plot(x, Dt_plane(x, r), color=blue[0], label="D one plane parallel")
N = 10
#plt.plot(x, Dn_plane(x, r, N=N), label="D plane normal")
plot_dashed(x, Dt_equal_planes(x, r), color=blue[1], label="D two planes")
plot_dashed(x, Dt_3R_planes(x, r), color=blue[2], label="D two planes ecc=0.50")

plt.plot(x, Dt_cyl_ax(x, r), color="g", label="D cylinder parallel")
#plt.plot(x, Dt_cyl_ax_large_dist(x, r), label="D cylinder parallel far")

colors = ["#99ff99", "#990000", "#ff0000", "#ff9999"]
for i, eccentricity in enumerate([0., 0.5, 0.75, 0.99]):
    plot_dashed(x, Dt_cyl(x/(1. - eccentricity), r, eccentricity),
        color=colors[i], label="D cylinder parallel ecc=%.2f" % eccentricity)

plot_dashed(x, D_sphere(x, r), color="black", thresh=5., label="D sphere")

plt.xlim(1., Rmax)
plt.ylim(0., 1.)
plt.xscale("log")
ticks = [1, 2, 3, 5, 10, 20, 30, 50, 100]
plt.xticks(ticks, ticks)
plt.grid()
plt.xlabel("distance to wall in multiples of particle radius")
plt.ylabel("relative diffusivity")
plt.legend(loc="lower left", bbox_to_anchor=(0.6, 0.))
#plt.legend(loc="lower right")
#plt.show()

from nanopores import DROPBOX, assertdir
#from matplotlib.backends.backend_pdf import PdfPages
DIR = DROPBOX + "/figures/diffusivity/"
assertdir(DIR)
#pp = PdfPages(DIR + "analytical.pdf")
#pp.savefig(plt.gcf())
#plt.gcf().set_size_inches(11.69, 8.27) # A4
plt.savefig(DIR + "analytical.pdf", papertype = 'a4', orientation = 'portrait', bbox_inches="tight")