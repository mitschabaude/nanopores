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
from nanopores import kT, eta, eperm, qq, rpermw, HOME
from matplotlib import pyplot as plt
#from matplotlib import patches
import matplotlib.animation as animation
from matplotlib import collections
from functools import partial
#import matplotlib
#matplotlib.use("Agg")

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
        self.x = self.x + np.array(dx).reshape(-1, 1)
#        if self.circle is not None:
#            self.circle.center = self.x[::2]

class Plane(object):

    def __init__(self, p, n): # outer normal vextor and one point on plane
        self.n = np.array(n).reshape(1, 3)
        self.p = np.array(p).reshape(1, 3)

    def reflect(self, P):
        # reflect all points in P downwards that lie above plane
        # ("above" is direction of normal vector)
        x = np.array([p.x.flatten() for p in P])
        # n * (x - p) > 0
        excess = np.dot(x - self.p, self.n.T)
        above = (excess > 0).flatten()
        dx = np.zeros_like(x)
        dx[above, :] = (-2.*excess*self.n)[above, :]

        for i, p in enumerate(P):
            if above[i]:
                p.move(dx[i, :])

class Box(object):
    # reflect to be inside box
    def __init__(self, x0, Lx, Ly, Lz):
        ex = np.eye(1, 3, 0)
        ey = np.eye(1, 3, 1)
        ez = np.eye(1, 3, 2)
        x0 = np.array(x0).reshape(1, 3)
        self.planes = [Plane(x0, -ex), Plane(x0 + Lx*ex, ex),
                       Plane(x0, -ey), Plane(x0 + Ly*ey, ey),
                       Plane(x0, -ez), Plane(x0 + Lz*ez, ez)]

    def reflect(self, P):
        for plane in self.planes:
            plane.reflect(P)


def msqrt(M): # matrix square root
    #U, S, V = np.linalg.svd(M)
    #return np.dot(U, np.dot(np.diag(np.sqrt(S)), V))
    return np.linalg.cholesky(M)

# external force: constant electrical field
def f_extern(p):
    #return np.zeros_like(p.x)
    a = 0e-14
    #b = 1e-26
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

def f_shortrange(P):
    n = len(P)
    a = 1e-9*np.array([p.a for p in P])
    apa = a[:, None] + a[None, :]

    x = 1e-9*np.array([p.x.flatten() for p in P])
    R = x[:, None, :] - x[None, :, :]
    r = np.sqrt(np.sum(R**2, 2)) + 1e-100

    #E0 = apa*1e9*10.*kT # total energy required for r = apa*1.1 --> r = 0
    #E = E0*((r/apa/1.1 - 1.)**2)
    E0 = 1e9*1e1*kT
    cutoff = 1.05
    f = 2./cutoff*E0*np.maximum(1. - r/apa/cutoff, 0.)

    R0 = R / r[:, :, None]
    F = f[:, :, None] * R0
    ff = np.sum(F, 1).T.reshape(3*n, 1)
    return ff


def simulation(P, T=100., dt=1.):
    t = 0.
    while t < T:
        # calculate forces
        force = f_electric(P) + f_shortrange(P)
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
#from time import time

def move(P, dt=1., boundaries=()):
    # calculate forces
    force = f_electric(P) + f_shortrange(P)
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

    for b in boundaries:
        b.reflect(P)

class ParticleCollection(object):
    def __init__(self, particles=None):
        if particles is None:
            particles = []
        self.P = particles

    def add(self, generate, N):
        i = 0
        while i < N:
            p0 = generate()
            if all(np.linalg.norm(q.x - p0.x) > p0.a + q.a for q in self.P):
                self.P.append(p0)
                i += 1

def ellipse_collection(ax, P):
    xy = np.array([p.x[::2].flatten() for p in P])
    sizes = np.array([p.a for p in P])
    coll = collections.EllipseCollection(sizes, sizes, np.zeros_like(sizes),
               offsets=xy, units='x', facecolors=[p.color for p in P],
               transOffset=ax.transData, alpha=0.7)
    return coll

def panimate(ax, dt, pcoll, patches, boundaries, **kwargs):
    particles = pcoll.P
    coll = ellipse_collection(ax, particles)
    def init():
        return ()
        #ax.clear()
        #ax.add_patch(rect)
        #return (rect, )

    def animate(i):
        if i == 0:
            for patch in patches:
                ax.add_patch(patch)
            ax.add_collection(coll)
        else:
            move(particles, dt, boundaries)
            xy = np.array([p.x[::2].flatten() for p in particles])
            coll.set_offsets(xy)
        return tuple(patches + [coll])

    kwargs = dict(dict(frames=1800, interval=10, blit=True), **kwargs)
    ani = animation.FuncAnimation(ax.figure, animate, init_func=init, **kwargs)
    return ani

def maybe_save(do_save, ani, name):
    if do_save:
        ani.save(HOME + "/presentations/anaday17/" + name, fps=30, dpi=200,
                 writer="ffmpeg_file",
                 savefig_kwargs={"bbox_inches":0, "pad_inches":0},
                 extra_args=['-vcodec', 'libx264'])
    else:
        plt.show()


def random_particle(L, *args, **kwargs):
    x = [(2*np.random.rand() - 1.)*L, 0., (2*np.random.rand() - 1.)*L]
    return Particle(x, *args, **kwargs)

def setup_box(L):
    # ax, patches, boundaries = setup_box(L)
    fig = plt.figure()
    fig.set_size_inches(6, 6)

    ax = plt.axes([0,0,1,1], autoscale_on=False,
                     xlim=(-L, L), ylim=(-L, L))
    ax.set_axis_off()
    rect = plt.Rectangle([-L, -L], 2*L, 2*L, ec='k', lw=3, fc='none')
    box = Box([-L, -L, -L], 2*L, 2*L, 2*L)
    return ax, [rect], [box]
