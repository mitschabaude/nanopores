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
from nanopores import kT, eta, collect
from matplotlib import pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

I = np.matrix(np.eye(3))

class Particle(object):

    def __init__(self, x, a):
        self.x = np.reshape(np.array(x), (-1, 1))
        self.a = a
        self.circle = None

    def add_patch(self, ax):
        self.circle = patches.Circle(xy=self.x[::2], radius=self.a)
        ax.add_patch(self.circle)

    def move(self, dx):
        self.x = self.x + dx
        if self.circle is not None:
            self.circle.center = self.x[::2]

def msqrt(M): # matrix square root
    U, S, V = np.linalg.svd(M)
    return np.dot(U, np.dot(np.diag(np.sqrt(S)), V))

# external force: gravity + repulsive potential at the bottom
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

def mobility(p1, p2): # D/kT
    a1 = 1e-9*p1.a
    a2 = 1e-9*p2.a
    R = 1e-9*p1.x - 1e-9*p2.x
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

def simulation(particles, T=100., dt=1.):
    t = 0.
    while t < T:
        # calculate forces
        force = [f_extern(p) for p in particles]
        brownian = [f_brownian(p) for p in particles]

        # determine resulting velocities
        for p, up in collect(particles):
            #u = np.zeros(1, 3)
            for (iq, q), uq in collect(enumerate(particles)):
                M = mobility(p, q)
                sqM = np.matrix(msqrt(M))
                uq.new = M*force[iq] + np.sqrt(2.*kT/dt*1e9)*sqM*brownian[iq]

            up.new = reduce(lambda a, b: a + b, uq)

        # move particles
        for p, u in zip(particles, up):
            #p.x = p.x + dt*u
            p.move(dt*u)

        yield t
        t += dt

def move(particles, dt=1.):
    # calculate forces
    force = [f_extern(p) for p in particles]
    brownian = [f_brownian(p) for p in particles]

    # determine resulting velocities
    for p, up in collect(particles):
        #u = np.zeros(1, 3)
        for (iq, q), uq in collect(enumerate(particles)):
            M = mobility(p, q)
            sqM = np.matrix(msqrt(M))
            uq.new = M*force[iq] + np.sqrt(2.*kT/dt*1e9)*sqM*brownian[iq]

        up.new = reduce(lambda a, b: a + b, uq)

    # move particles
    for p, u in zip(particles, up):
        #p.x = p.x + dt*u
        p.move(dt*u)

# example simulation
#p = Particle([-1, 0, -1], 1)
#q = Particle([1, 0, 1], 0.5)
#particles = [p, q]

#for t in simulation(particles, T=100.):
#    print t
#    print p.x
#    print q.x

# animation
dt = 10.

fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

L = 100
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-L*1.1, L*1.1), ylim=(-L*1.1, L*1.1))

# rect is the box edge
rect = plt.Rectangle([-L, -L], 2*L, 2*L, ec='none', lw=2, fc='none')
ax.add_patch(rect)

def random_particle(L):
    x = [np.random.randn()*L/5., 0., np.random.randn()**2*L/5.]
    a = np.random.rand()
    return Particle(x, a)

particles = [random_particle(L) for i in range(50)]

def init():
    """initialize animation"""
#    global rect
#    rect.set_edgecolor('none')
#    for p in particles:
#        p.circle.set_color("none")
    return () #tuple([rect] + [p.circle for p in particles])

def animate(i):
    """perform animation step"""
    global particles, rect, dt, ax, fig
    if i == 0:
        ax.add_patch(rect)
        for p in particles:
            p.add_patch(ax)
    else:
        move(particles, dt)

    # update pieces of the animation
#    rect.set_edgecolor('k')
#    for p in particles:
#        p.circle.set_color("blue")
#    particles.set_data(box.state[:, 0], box.state[:, 1])
#    particles.set_markersize(ms)
    return tuple([rect] + [p.circle for p in particles])

ani = animation.FuncAnimation(fig, animate, frames=600,
                              interval=10, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
