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
from nanopores import kT, eta, eperm, qq, rpermw
from matplotlib import pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation
from matplotlib import collections

class Particle(object):

    def __init__(self, x, a, charge=0., color="blue"):
        self.x = np.reshape(np.array(x), (-1, 1))
        self.a = a
        self.circle = None
        self.charge = charge
        self.color = color
#
#    def add_patch(self, ax):
#        self.circle = patches.Circle(xy=self.x[::2], radius=self.a)
#        ax.add_patch(self.circle)
#
#    def patch(self):
#        self.circle = patches.Circle(xy=self.x[::2], radius=self.a)
#        return self.circle

    def move(self, dx):
        self.x = self.x + np.array(dx)
#        if self.circle is not None:
#            self.circle.center = self.x[::2]

def msqrt(M): # matrix square root
    #U, S, V = np.linalg.svd(M)
    #return np.dot(U, np.dot(np.diag(np.sqrt(S)), V))
    return np.linalg.cholesky(M)

# external force: constant electrical field
def f_extern(p):
    #return np.zeros_like(p.x)
    a = 0e-14
    b = 1e-26
    f = np.zeros_like(p.x)
    f[2, :] = -a*p.a**3 #+ b*(p.x[2, :] - 1.)**(-7)
    return f

def f_brownian(p):
    rand = np.random.randn(*p.x.shape)
    rand[1] = 0.
    return rand

def RR(R):
    #R = np.reshape(R, (-1, 1))
    r = np.linalg.norm(R)
    return np.matrix(R * R.T / r**2)

def mobility_pp(p1, p2): # D/kT
    a1 = 1e-9*p1.a
    a2 = 1e-9*p2.a
    R = 1e-9*p1.x - 1e-9*p2.x
    r = np.linalg.norm(R)
    I = np.matrix(np.eye(3))

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

def mobility_vectorized(P):
    # return matrix mobility(pi, pj)_kl where kl are the space coordinates
    # is, at the moment, list of particles
    # TODO: need pure numerical representation !!!!!!!!!!!!1111
    n = len(P)
    M = np.matrix(np.zeros((3*n, 3*n)))
    for i, p in enumerate(P):
        for j, q in enumerate(P):
            M[i::n, j::n] = mobility_pp(p, q)
    return M

def mobility(P):
    # nice: M is of the form Mijkl = Aij * Ikl + Bij * Rijkl
    # where Ikl, Rijkl are easily formed
    # obtain vector representations of positions, radii
    n = len(P)
    a = 1e-9*np.array([p.a for p in P])
    x = 1e-9*np.array([p.x.flatten() for p in P])
    rr = np.sum((x[:, None, :] - x[None, :, :])**2, 2) + 1e-100
    r = np.sqrt(rr)

    A = np.zeros_like(r)
    B = np.zeros_like(r)

    ama = a[:, None] - a[None, :]
    apa = a[:, None] + a[None, :]
    asq = a[:, None]**2 + a[None, :]**2

    case1 = r <= np.abs(ama)
    A[case1] = (1./(6.*np.pi*eta*np.maximum(a[:,None], a[None,:])))[case1]
    B[case1] = 0.

    case2 = ~case1 & (r <= apa)
    C = 1./(6.*np.pi*eta*a[:, None]*a[None, :])
    A[case2] = (C*(0.5*apa - (ama**2 + 3.*r**2)**2/(32.*r**3)))[case2]
    B[case2] = (C*(3.*(ama**2 - r**2)**2/(32.*r**3)))[case2]

    case3 = ~(case1 | case2) # else
    C = 1./(8.*np.pi*eta*r)
    A[case3] = (C*(1. + asq/(3.*r**2)))[case3]
    B[case3] = (C*(1. - asq/r**2))[case3]

    I = np.eye(3)
    R = (x[:, None, :] - x[None, :, :])
    RR = (R[:, :, :, None] * R[:, :, None, :]).transpose(2, 0, 3, 1)
    RR = RR / (rr[None, :, None, :])

    M = (A[None, :, None, :] * I[:, None, :, None]).reshape(3*n, 3*n) + \
        (B[None, :, None, :] * RR).reshape(3*n, 3*n)
    return np.matrix(M)


def f_vectorized(P, f):
    n = len(P)
    F = np.zeros((3*n, 1))
    for i, p in enumerate(P):
        F[i::n] = f(p)
    return F

def f_electric(P):
    n = len(P)
    a = 1e-9*np.array([p.a for p in P])
    apa = a[:, None] + a[None, :]


    x = 1e-9*np.array([p.x.flatten() for p in P])

    R = x[:, None, :] - x[None, :, :]
    r = np.sqrt(np.sum(R**2, 2) + 1e-100)
    R0 = R / (r**3)[:, :, None]

    q = np.array([float(p.charge) for p in P])
    const = qq**2 / (4.*np.pi*eperm*rpermw)

    QQ = q[:, None] * q[None, :]
    F = const * QQ[:, :, None] * R0
    #F[np.diag_indices_from(r)] = 0.
    tooclose = r <= apa
    R0i = R / (np.maximum(a[:, None], a[None, :])**3)[:, :, None]
    F[tooclose] = (const * QQ[:, :, None] * R0i)[tooclose]

    f = np.sum(F, 1).T.reshape(3*n, 1)

    return f


def simulation(P, T=100., dt=1.):
    t = 0.
    while t < T:
        # calculate forces
        force = f_electric(P)
        brownian = f_vectorized(P, f_brownian)
        M = mobility(P)
        sqM = np.matrix(msqrt(M))

        # determine resulting velocities
        U = M*force + np.sqrt(2.*kT/dt*1e9)*sqM*brownian

        n = len(P)
        # move particles
        for i, p in enumerate(P):
            u = U[i::n]
            p.move(dt*u)

        yield t
        t += dt

# static particles: no forces act on them, are not affected by
from time import time

def move(P, dt=1.):
    # calculate forces
    force = f_electric(P)
    brownian = f_vectorized(P, f_brownian)
    #t = time()
    M = mobility(P)
    #print "forming mobility", time() - t
    #t = time()
    sqM = np.matrix(msqrt(M))
    #print "matrix square root", time() - t

    # determine resulting velocities
    Udet = M*force
    Ubro = np.sqrt(2.*kT/dt*1e9)*sqM*brownian
    #U = M*force + np.sqrt(2.*kT/dt*1e9)*sqM*brownian
    U = Udet + Ubro
    #print "P0:", Udet[0], Ubro[0]
    #print "P1:", Udet[1], Ubro[1]
    #print

    n = len(P)
    # move particles
    for i, p in enumerate(P):
        u = U[i::n]
        p.move(dt*u)

# animation
dt = .001

fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

L = 5
rlarge = 2.1
rsmall = .11
Qlarge = 50
N = 100

Np = (N - Qlarge)/2
Nm = (N + Qlarge)/2

ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-L*1.1, L*1.1), ylim=(-L*1.1, L*1.1))

# rect is the box edge
rect = plt.Rectangle([-L, -L], 2*L, 2*L, ec='k', lw=2, fc='none')

def random_particle(L, *args, **kwargs):
    x = [(2*np.random.rand() - 1.)*L, 0., (2*np.random.rand() - 1.)*L]
    a = rsmall
    return Particle(x, a, *args, **kwargs)

large = Particle([0.,0.,0.], rlarge, 50., "#9999ff")

particles = [large]
ccolors = {1: "blue", -1: "red"}
for q, N in (-1, Nm), (1, Np):
    newparticles = []
    while len(newparticles) < N:
        p0 = random_particle(L, charge=q, color=ccolors[q])
        if all(np.linalg.norm(q.x - p0.x) > p0.a + q.a for q in particles + newparticles):
            newparticles.append(p0)
    particles = particles + newparticles

xy = np.array([p.x[::2].flatten() for p in particles])
sizes = np.array([p.a for p in particles])
coll = collections.EllipseCollection(sizes, sizes, np.zeros_like(sizes),
                                     offsets=xy, units='x',
                                     facecolors=[p.color for p in particles],
                                     transOffset=ax.transData, alpha=0.7)

def init():
    return ()

def animate(i):
    global rect, coll
    if i == 0:
        ax.add_patch(rect)
        ax.add_collection(coll)
    else:
        move(particles, dt)
        xy = np.array([p.x[::2].flatten() for p in particles])
        coll.set_offsets(xy)
    return rect, coll

ani = animation.FuncAnimation(fig, animate, frames=600,
                              interval=10, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
