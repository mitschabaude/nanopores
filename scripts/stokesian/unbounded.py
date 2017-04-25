# (c) 2017 Gregor Mitscha-Baude
"""stokesian dynamics simulation of two kinds of particles

framework: given list of particle centers xi and radii ri, for each timestep, calculate:
1.) external force on each particle (and torque?) F(xi, ri) or even F(x1,...,xn; r1,...,rn)
2.) brownian (random) force (and torque?) Frand_i
3.) velocities via ui = sum_j ( M_ij * F_j ) where M_ij(x_i, x_j) is the
    3x3 translational - translational mobility matrix.
4.) next position x_i = x_i + ui * delta t

rotational velocities can be ignored in this framework because they do not
change the particle position."""
import numpy as np
from nanopores import kT, eta

I = np.eye(3)

class Particle(object):

    def __init__(self, x, a):
        self.x = np.array(x)
        self.a = a

# external force: gravity + repulsive potential at the bottom
def f_gravity(p):
    a = 1.
    b = 1e-2
    f = np.zeros_like(p.x)
    f[:, 3] = -a*p.a**3 + b*p.x[:, 3]**(-7)
    return f

def f_brownian(p):
    return np.random.randn(p.x.shape)

def RR(R):
    R = np.reshape(R, [1, -1])
    r = np.linalg.norm(R)
    return R.T * R / r**2

def mobility(p1, p2): # D/kT
    a1 = p1.a
    a2 = p2.a
    R = p1.x - p2.x
    r = np.linalg.norm(R)

    if r <= abs(a1 - a2):
        return 1./(6.*np.pi*eta*max(a1, a2)) * I
    elif r <= a1 + a2:
        A = 0.5*(a1 + a2) - ((a1-a2)**2 + 3.*r**2)**2/(32.*r**3)
        B = 3.*((a1-a2)**2 - r**2)**2/(32.*r**3)
        return 1./(6.*np.pi*eta*a1*a2) * (A*I + B*RR(R))
    else:
        a = a1**2 + a2**2
        A = 1. + a/(3.*r**2)
        B = 1. - a/r**2
        return 1./(8.*np.pi*eta*r) * (A*I + B*RR(R))


