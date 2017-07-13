# (c) 2017 Gregor Mitscha-Baude
"random walk of many particles in cylindrical pore"
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections

import dolfin
import nanopores
from nanopores.tools import fields
from nanopores.geometries.allpores import get_pore
fields.set_dir_dropbox()
dolfin.parameters["allow_extrapolation"] = False
params = nanopores.user_params(
    geoname = "wei",
    rMolecule = 6.,
    Qmol = -1.,
    bV = -0.5,

    N = 10,
    dt = 1.,
)

# external forces
def load_externals(**params):
    F, = fields.get_functions("wei_force_ps", "F", **params)
    D, = fields.get_functions("wei_D_2D", "D", **params)
    return F, D

# initial positions: uniformly distributed over disc
def initial(R, z, N=10):
    # create uniform polar coordinates r**2, theta
    r2 = R**2*np.random.rand(N)
    theta = 2.*np.pi*np.random.rand(N)
    r = np.sqrt(r2)
    x = np.zeros((N, 3))
    x[:, 0] = r*np.cos(theta)
    x[:, 1] = r*np.sin(theta)
    x[:, 2] = z
    return x, r, x[:, 2]

class RandomWalk(object):

    def __init__(self, pore, N=10, dt=1.):
        # dt is timestep in nanoseconds
        self.params = pore.params
        R = self.params.R
        ztop = pore.protein.zmax()[1]
        x, r, z = initial(R, ztop, N)
        self.N = N
        self.x = x
        self.r = r
        self.z = z
        self.dt = dt

    def walk(self):
        "one step of random walk"

    def ellipse_collection(self, ax):
        "for matplotlib plotting"
        xz = self.x[:,[0,2]]
        sizes = self.params.rMolecule*np.ones(self.N)
        colors = ["b"]*self.N
        coll = collections.EllipseCollection(sizes, sizes, np.zeros_like(sizes),
                   offsets=xz, units='x', facecolors=colors,
                   transOffset=ax.transData, alpha=0.7)
        return coll

    def move_ellipses(self, coll):
        xz = self.x[:,[0,2]]
        coll.set_offsets(xz)

def animate(rw, **aniparams):
    R = rw.params.R
    H = rw.params.H

    fig = plt.figure()
    fig.set_size_inches(6, 6)
    ax = plt.axes([0,0,1,1], autoscale_on=False, xlim=(-R, R), ylim=(-H, H))
    coll = rw.ellipse_collection(ax)

    def init():
        return ()

    def animate(i):
        if i == 0:
            ax.add_collection(coll)
        else:
            rw.walk()
            rw.move_ellipses(coll)
        return tuple(coll)

    aniparams = dict(dict(frames=1800, interval=10, blit=True), **aniparams)
    animation.FuncAnimation(ax.figure, animate, init_func=init, **aniparams)
    plt.show()

pore = get_pore(**params)
rw = RandomWalk()

dt = 1e-9
ax, patches, boundaries = setup_box(L)
ani = panimate(ax, dt, pcoll, patches, boundaries, frames=1800)
maybe_save(True, ani, "two_ions_only.mp4")