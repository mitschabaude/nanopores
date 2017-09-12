# (c) 2017 Gregor Mitscha-Baude
"random walk of many particles in cylindrical pore"
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections
import matplotlib.patches as mpatches
from scipy.stats import poisson, gamma

import dolfin
import nanopores
from nanopores import get_pore
from nanopores.tools.polygons import Polygon, Ball, isempty
from nanopores.models import nanopore
from nanopores.tools.poreplots import streamlines
from nanopores.tools import fields

dolfin.parameters["allow_extrapolation"] = True

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 6.,
    h = 5.,
    Nmax = 4e4,
    Qmol = -1.,
    bV = -0.2,
    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = 10., # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 20.,
    margbot = 10.,
)

# domains are places where molecule can bind
# and/or be reflected after collision
domain_params = dict(
    cyl = False, # determines whether rz or xyz coordinates are passed to .inside
    walldist = 1.5, # multiple of radius that determines what counts as collision

    exclusion = True,
    minsize = 0.01, # accuracy when performing reflection

    binding = False,
    eps = 1., # margin in addition to walldist, determines re-attempting
    p = 0.1, # binding probability for one attempt
    t = 1e6, # mean of exponentially distributed binding duration [ns]
    dx = 0.4, # width of bond energy barrier [nm]
    use_force = True, # if True, t_mean = t*exp(-|F|*dx/kT)
)

class Domain(object):
    """based on existing domain object which need only support the
    methods .inside and .inside_single, which get either xyz or rz
    coordinates depending on the cyl attribute set in init."""

    def __init__(self, domain, **params):
        self.domain = domain
        self.__dict__.update(domain_params)
        self.__dict__.update(params)
        if isinstance(domain, Polygon):
            self.cyl = True

    def collide(self, rw):
        "compute collisions and consequences with RandomWalk instance"
        radius = rw.params.rMolecule * self.walldist

        # determine collisions and then operate only on collided particles
        X = rw.rz if self.cyl else rw.x[rw.alive]
        collided = self.domain.inside(X, radius=radius)

        # "reflect" particles by shortening last step
        if self.exclusion:
            X0, X1 = rw.xold[rw.alive], rw.x[rw.alive]
            for i in np.nonzero(collided)[0]:
                x = self.binary_search_inside(X0[i], X1[i], radius)
                rw.update_one(i, x)

        # attempt binding for particles that can bind
        if self.binding:
            can_bind = rw.can_bind[rw.alive]
            attempt = collided & can_bind
            # bind with probability p
            bind = np.random.rand(np.sum(attempt)) <= self.p
            # draw exponentially distributed binding time
            duration = self.draw_binding_durations(attempt, bind, rw)
            # update can_bind and bind_times of random walk
            iattempt = rw.i[rw.alive][attempt]
            ibind = iattempt[bind]
            rw.can_bind[iattempt] = False
            rw.bind_times[ibind] += duration
            # some statistics
            rw.attempts[rw.i[rw.alive][attempt]] += 1
            rw.bindings[rw.i[rw.alive][attempt][bind]] += 1

            # unbind particles that can not bind and are out of nobind zone
            X_can_not_bind = X[~can_bind]
            rnobind = radius + self.eps
            unbind = ~self.domain.inside(X_can_not_bind, radius=rnobind)
            iunbind = rw.i[rw.alive][~can_bind][unbind]
            rw.can_bind[iunbind] = True

    def binary_search_inside(self, x0, x1, radius):
        if self.domain.inside_single(x0, radius=radius):
            print "ERROR: x0 is in domain despite having been excluded before."
            print "x0", x0, "x1", x1
            #raise Exception
        if np.sum((x0 - x1)**2) < self.minsize**2:
            return x0
        x05 = .5*(x0 + x1)
        if self.domain.inside_single(x05, radius=radius):
            x1 = x05
        else:
            x0 = x05
        return self.binary_search_inside(x0, x1, radius)

    def draw_binding_durations(self, attempt, bind, rw):
        if self.use_force and np.sum(bind) > 0:
            # evaluate force magnitude at binding particles
            ibind = np.nonzero(attempt)[0][bind]
            F = np.array([rw.F(x) for x in rw.rz[ibind]])
            F = np.sqrt(np.sum(F**2, 1))
            # create array of mean times
            kT = rw.phys.kT
            dx = 1e-9*self.dx
            t = self.t * np.exp(-F*dx/kT)
        else:
            t = self.t
        return np.random.exponential(t, np.sum(bind))

# external forces
def load_externals(**params):
    return nanopore.force_diff(**params)

class RandomWalk(object):

    def __init__(self, pore, N=10, dt=1., walldist=2.,
                 margtop=20., margbot=10., xstart=None, zstart=None,
                 rstart=None, **params):
        # dt is timestep in nanoseconds
        self.pore = pore
        self.params = pore.params
        self.params.update(params, margtop=margtop, margbot=margbot,
                           walldist=walldist, dt=dt,
                           rstart=rstart, xstart=xstart, zstart=zstart)

        # initialize some parameters and create random walkers at entrance
        self.rtop = pore.protein.radiustop() - self.params.rMolecule
        self.ztop = pore.protein.zmax()[1]
        self.rbot = pore.protein.radiusbottom() - self.params.rMolecule
        self.zbot = pore.protein.zmin()[1]
        self.zmid = .5*(self.ztop + self.zbot) if isempty(pore.membrane) else pore.params["zmem"]
        self.N = N

        x, r, z = self.initial()

        self.x = x
        self.xold = x
        self.rz = np.column_stack([r, z])
        self.dt = dt
        self.t = 0.
        self.i = np.arange(N)

        # load force and diffusivity fields
        F, D, divD = load_externals(**params)
        self.F = F #self.ood_evaluation(F)
        self.D = D #self.ood_evaluation(D)
        self.divD = divD #self.ood_evaluation(divD)
        self.phys = nanopores.Physics("pore_mol", **params)

        self.alive = np.full((N,), True, dtype=bool)
        self.success = np.full((N,), False, dtype=bool)
        self.fail = np.full((N,), False, dtype=bool)
        self.can_bind = np.full((N,), True, dtype=bool)
        self.times = np.zeros(N)
        self.bind_times = np.zeros(N)
        self.attempts = np.zeros(N, dtype=int)
        self.bindings = np.zeros(N, dtype=int)

        self.domains = []
        self.add_domain(pore.protein, binding=False, exclusion=True,
                        walldist=walldist)
        if not isempty(pore.membrane):
            self.add_domain(pore.membrane, binding=False, exclusion=True,
                        walldist=walldist)

    # initial positions: uniformly distributed over disc
    def initial(self):
        rstart = self.params.rstart
        xstart = self.params.xstart
        zstart = self.params.zstart
        if rstart is None:
            rstart = self.rtop - self.params.rMolecule*(self.params.walldist - 1.)
        if xstart is None:
            xstart = 0.
        if zstart is None:
            zstart = self.ztop

        # create uniform polar coordinates r, theta
        r = rstart * np.sqrt(np.random.rand(self.N))
        theta = 2.*np.pi * np.random.rand(self.N)

        x = np.zeros((self.N, 3))
        x[:, 0] = xstart + r*np.cos(theta)
        x[:, 1] = r*np.sin(theta)
        x[:, 2] = zstart
        return x, np.sqrt(x[:, 0]**2 + x[:, 1]**2), x[:, 2]

    def add_domain(self, domain, **params):
        """add domain where particles can bind and/or are excluded from.
        domain only has to implement the .inside(x, radius) method.
        params can be domain_params"""
        self.domains.append(Domain(domain, **params))

    def add_wall_binding(self, **params):
        self.domains[0].__dict__.update(params, binding=True)

    def ood_evaluation(self, f):
        dim = self.params.dim
        def newf(x):
            try:
                return f(x)
            except RuntimeError:
                print "ood:", x
                return np.zeros(dim)
        return newf

    def evaluate(self, function):
        return np.array([function(rz) for rz in self.rz])

    def evaluate_vector_cyl(self, function):
        r = self.rz[:, 0] + 1e-30
        R = self.x[self.alive] / r[:, None]
        F = self.evaluate(function)
        return np.column_stack([F[:, 0]*R[:, 0], F[:, 0]*R[:, 1], F[:, 1]])

    def evaluate_D_cyl_matrix(self):
        # approximation based on Dn \sim Dt
        D = self.evaluate(self.D)
        Dn = D[:, 0]
        Dt = D[:, 1]
        r = self.rz[:, 0] + 1e-30
        xbar = self.x[:, 0]/r
        ybar = self.x[:, 1]/r

        Dmatrix = np.zeros((self.N, 3, 3))
        Dmatrix[:, 0, 0] = Dn*xbar**2 + Dt*(1.-xbar**2)
        Dmatrix[:, 1, 1] = Dn*ybar**2 + Dt*(1.-ybar**2)
        Dmatrix[:, 2, 2] = Dt
        return Dmatrix

    def evaluate_D_cyl(self):
        # approximation based on Dn \sim Dt
        D = self.evaluate(self.D)
        Dn = D[:, 0]
        Dt = D[:, 1]
        r = self.rz[:, 0]
        xbar = self.x[self.alive, 0]/r
        ybar = self.x[self.alive, 1]/r
        Dx = Dn*xbar**2 + Dt*(1.-xbar**2)
        Dy = Dn*ybar**2 + Dt*(1.-ybar**2)
        return np.column_stack([Dx, Dy, Dt])

    def evaluate_D_simple(self):
        # just take D = Dzz
        D = self.evaluate(self.D)
        return D[:, 1, None]

    def brownian(self, D):
        n = np.count_nonzero(self.alive)
        zeta = np.random.randn(n, 3)
        return np.sqrt(2.*self.dt*1e9*D) * zeta

    def update(self, dx):
        self.xold = self.x.copy()
        #self.rzold = self.rz.copy()
        self.x[self.alive] = self.x[self.alive] + dx
        self.update_alive()
        r = np.sqrt(np.sum(self.x[self.alive, :2]**2, 1))
        self.rz = np.column_stack([r, self.x[self.alive, 2]])

    def is_success(self, r, z):
        return (z < self.zbot - self.params.margbot) | (
               (r > self.rbot + self.params.margbot) & (z < self.zbot))

    def is_fail(self, r, z):
        return (z > self.ztop + self.params.margtop) | (
               (r > self.rtop + self.params.margtop) & (z > self.ztop))

    def update_alive(self):
        alive = self.alive
        z = self.x[alive, 2]
        r = np.sqrt(np.sum(self.x[alive, :2]**2, 1))
        self.success[alive] = self.is_success(r, z)
        self.fail[alive] = self.is_fail(r, z)
        died = self.fail[alive] | self.success[alive]
        self.alive[alive] = ~died
        self.times[alive] = self.t

    def update_one(self, i, xnew):
        self.x[np.nonzero(self.alive)[0][i]] = xnew
        self.rz[i, 0] = np.sqrt(xnew[0]**2 + xnew[1]**2)
        self.rz[i, 1] = xnew[2]

    def step(self):
        "one step of random walk"
        # evaluate F and D
        D = self.evaluate_D_cyl()
        F = self.evaluate_vector_cyl(self.F)
        divD = 1e9*self.evaluate_vector_cyl(self.divD)
        kT = self.phys.kT
        dt = self.dt
        self.t += self.dt

        # get step
        dW = self.brownian(D)
        dx = dW + dt*divD + dt*D/kT*F
        #print "%.2f (dx) = %.2f (dW) + %.2f (divD) + %.2f (F)" % (
        #        abs(dx[0, 2]), abs(dW[0, 2]), abs(dt*divD[0, 2]), abs((dt*D/kT*F)[0, 2]))
        #print ("t = %.2f microsec" % (self.t*1e-3))

        # update position and time and determine which particles are alive
        self.update(dx)

        # correct particles that collided with pore wall
        #self.simple_reflect()
        for domain in self.domains:
            domain.collide(self)

    def walk(self):
        with nanopores.Log("Running..."):
            yield self.t
            while np.any(self.alive):
                #with nanopores.Log("%.0f ns, cpu time:" % self.t):
                self.step()
                yield self.t

        self.finalize()

    def finalize(self):
        print "finished!"
        print "mean # of attempts:", self.attempts.mean()
        print "mean # of bindings:", self.bindings.mean()
        print "mean dwell time with binding: %.3f mus" % (
            1e-3*(self.bind_times + self.times).mean())
        print "mean dwell time without binding: %.3f mus" % (
            1e-3*self.times.mean())
        self.times += self.bind_times

    def save(self, name="rw"):
        if "N" in self.params:
            self.params.pop("N")
        fields.save_fields(name, self.params,
            times = self.times,
            success = self.success,
            fail = self.fail,
            bind_times = self.bind_times,
            attempts = self.attempts,
            bindings = self.bindings,
            )
        fields.update()

    def ellipse_collection(self, ax):
        "for matplotlib plotting"
        xz = self.x[:, [0,2]]
        #xz = self.rz
        sizes = self.params.rMolecule*np.ones(self.N)
        colors = ["b"]*self.N
        coll = collections.EllipseCollection(sizes, sizes, np.zeros_like(sizes),
                   offsets=xz, units='x', facecolors=colors,
                   transOffset=ax.transData, alpha=0.7)
        return coll

    def move_ellipses(self, coll, cyl=False):
        xz = self.x[:, ::2] if not cyl else np.column_stack(
           [np.sqrt(np.sum(self.x[:, :2]**2, 1)), self.x[:, 2]])
        coll.set_offsets(xz)
        #inside = self.inside_wall()
        #margin = np.nonzero(self.alive)[0][self.inside_wall(2.)]
        colors = np.full((self.N,), "b", dtype=str)
        #colors[margin] = "r"
        colors[self.success] = "k"
        colors[self.fail] = "k"
        colors[self.alive & ~self.can_bind] = "r"
        #colors = [("r" if inside[i] else "g") if margin[i] else "b" for i in range(self.N)]
        coll.set_facecolors(colors)
        #y = self.x[:, 1]
        #d = 50.
        #sizes = self.params.rMolecule*(1. + y/d)
        #coll.set(widths=sizes, heights=sizes)

    def polygon_patches(self, cyl=False):
        poly_settings = dict(closed=True, facecolor="#eeeeee", linewidth=1.,
                        edgecolor="k")
        ball_settings = dict(facecolor="#aaaaaa", linewidth=1., edgecolor="k",
                             alpha=0.5)
        patches = []
        for dom in self.domains:
            dom = dom.domain
            if isinstance(dom, Polygon):
                polygon = dom.nodes
                polygon = np.array(polygon)
                patches.append(mpatches.Polygon(polygon, **poly_settings))
                if not cyl:
                    polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])
                    patches.append(mpatches.Polygon(polygon_m, **poly_settings))
            elif isinstance(dom, Ball):
                xy = dom.x0[0], dom.x0[2]
                p = mpatches.Circle(xy, dom.r, **ball_settings)
                p.set_zorder(200)
                patches.append(p)

        return patches
    def plot_streamlines(self, cyl=False):
        # TODO:
        pass

def load_results(name, **params):
    data = fields.get_fields(name, **params)
    load = lambda a: a if isinstance(a, np.ndarray) else a.load()
    data = nanopores.Params({k: load(data[k]) for k in data})
    print "Found %d simulated events." % len(data.times)
    return data

def video(rw, cyl=False, **aniparams):
    R = rw.params.R
    Htop = rw.params.Htop
    Hbot = rw.params.Hbot

    #fig = plt.figure()
    #fig.set_size_inches(6, 6)
    #ax = plt.axes([0,0,1,1], autoscale_on=False, xlim=(-R, R), ylim=(-H, H))
    xlim = (-R, R) if not cyl else (0., R)
    ax = plt.axes(xlim=xlim, ylim=(-Hbot, Htop))
    #streamlines(rx=R, ry=Htop, Nx=100, Ny=100, maxvalue=None, F=rw.F)
    coll = rw.ellipse_collection(ax)
    patches = rw.polygon_patches(cyl)

    def init():
        return ()

    def animate(t):
        if t == 0:
            ax.add_collection(coll)
            for p in patches:
                ax.add_patch(p)

        rw.move_ellipses(coll, cyl=cyl)
        return tuple([coll] + patches)

    aniparams = dict(dict(interval=10, blit=True, save_count=5000), **aniparams)
    ani = animation.FuncAnimation(ax.figure, animate, frames=rw.walk(),
                                  init_func=init, **aniparams)
    return ani

def integrate_hist(hist, cutoff):
    n, bins, _ = hist
    I, = np.nonzero(bins > cutoff)
    return np.dot(n[I[:-1]], np.diff(bins[I]))

def integrate_values(T, fT, cutoff):
    values = 0.5*(fT[:-1] + fT[1:])
    I, = np.nonzero(T > cutoff)
    return np.dot(values[I[:-1]], np.diff(T[I]))

def exponential_hist(times, a, b, **params):
    cutoff = 0.03 # cutoff frequency in ms
    if len(times) == 0:
        return
    bins = np.logspace(a, b, 100)
    hist = plt.hist(times, bins=bins, alpha=0.5, **params)
    plt.xscale("log")
    params.pop("label")
    color = params.pop("color")
    total = integrate_hist(hist, cutoff)
    if sum(times > cutoff) == 0:
        return
    tmean = times[times > cutoff].mean()
    T = np.logspace(a-3, b, 1000)
    fT = np.exp(-T/tmean)*T/tmean
    fT *= total/integrate_values(T, fT, cutoff)
    plt.plot(T, fT, label="exp. fit, mean = %.2f ms" % (tmean,),
             color="dark" + color, **params)
    plt.xlim(10**a, 10**b)

def histogram(rw, a=-3, b=3, scale=1e-3):
    t = rw.times * 1e-9 / scale # assuming times are in nanosaconds

    exponential_hist(t[rw.success], a, b, color="green", label="translocated")
    exponential_hist(t[rw.fail], a, b, color="red", label="did not translocate")

    plt.xlabel(r"$\tau$ off [ms]")
    plt.ylabel("count")
    plt.legend(loc="best")

def hist_poisson(rw, name="attempts", ran=None, n=10):
    attempts = getattr(rw, name)
    if ran is None:
        n0 = 0
        n1 = n
    else:
        n0, n1 = ran

    k = np.arange(n0, n1 + 1)

    plt.figure()
    plt.hist(attempts, bins=(k - 0.5), label=name, color="#aaaaff")
    # poisson fit
    if n0 >= 1:
        mean = attempts[attempts >= n0].mean()
        a = poisson_from_positiveK(mean)
        print "Poisson fit, mean of K>0:", mean
        print "Inferred total mean:", a
        print "Mean of all K:", attempts.mean()
        K = len(attempts[attempts > 0])/(1.-np.exp(-a))
    else:
        a = attempts.mean()
        K = len(attempts)
    k0 = np.linspace(n0, n1, 500)
    plt.plot(k0, K*gamma.pdf(a, k0 + 1), "-",
             label="Poisson fit, mean = %.4f" %a, color="b")
    plt.plot(k, K*poisson.pmf(k, a), "o", color="b")
    plt.xlim(n0 - 0.5, n1 + 0.5)
    plt.xticks(k, k)
    plt.xlabel("# %s" % name)
    plt.ylabel("count")
    plt.legend()

def solve_newton(C, f, f1, x0=1., n=20):
    "solve f(x) == C"
    x = x0 # initial value
    print "Newton iteration:"
    for i in range(n):
        dx = -(f(x) - C)/f1(x)
        x = x + dx
        res = abs(f(x) - C)
        print i, "Residual", res, "Value", x
        if res < 1e-12:
            break
    print
    return x

def poisson_from_positiveK(mean):
    # solve x/(1 - exp(-x)) == mean
    def f(x):
        return x/(1. - np.exp(-x))
    def f1(x):
        return (np.expm1(x) - x)/(2.*np.cosh(x) - 2.)

    x = solve_newton(mean, f, f1, mean, n=10)
    return x

def save(ani, name="rw"):
    ani.save(nanopores.HOME + "/presentations/nanopores/%s.mp4" % name,
                 fps=30, dpi=200, writer="ffmpeg_file",
                 extra_args=['-vcodec', 'libx264'])

# convenience function for interactive experiments
def run(rw, name="rw", plot=True, a=-1, b=4, **aniparams):
    params = nanopores.user_params(video=False, save=False, cyl=False)
    if params.video:
        ani = video(rw, cyl=params.cyl, **aniparams)
        if params.save:
            save(ani, name=name)
        else:
            plt.show()
    else:
        for t in rw.walk(): pass
    if plot and ((not params.video) or (not params.save)):
        histogram(rw, a, b)
        hist_poisson(rw, "attempts")
        hist_poisson(rw, "bindings")
        plt.show()

if __name__ == "__main__":
    pore = get_pore(**params)
    rw = RandomWalk(pore, **params)
    receptor = Ball([9., 0., -30.], 8.)
    rw.add_domain(receptor, exclusion=True, walldist=1.,
                  binding=True, eps=1., t=1.5e6, p=0.14)
    run(rw, name="wei")

