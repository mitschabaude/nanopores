# (c) 2017 Gregor Mitscha-Baude
"random walk of many particles in cylindrical pore"
# TODO: bug: walldist = 2 works like what would be expected for = 1
# have hardcoded *2 in .contains_point
# => properly investigate matplotlib.path.contains_point
# and make code transparent
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import animation
from matplotlib import collections
from scipy.stats import poisson, gamma

import dolfin
import nanopores
from nanopores import get_pore
from nanopores.tools.polygons import Polygon, Ball, isempty
from nanopores.models import nanopore
from nanopores.tools.poreplots import streamlines
from nanopores.tools import fields

dolfin.parameters["allow_extrapolation"] = True #False #TODO:

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 6.,
    h = 5.,
    Nmax = 4e4,
    Qmol = -1.,
    bV = -0.2,
    posDTarget = True,
    x0 = None,
    
    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = 10., # time step [ns]
    walldist = 1., # in multiples of radius, should be >= 1
    margtop = 40.,
    margbot = 20.,
    initial = "disc",  # oder "sphere"
)

# domains are places where molecule can bind
# and/or be reflected after collision
domain_params = dict(
    cyl = False, # determines whether rz or xyz coordinates are passed to .inside
    walldist = 1., # multiple of radius that determines what counts as collision

    exclusion = True,
    minsize = 0.01, # accuracy when performing reflection

    binding = False,
    bind_type = "collision", # or "zone"
    # "collision" parameters:
    eps = 1., # margin in addition to walldist, determines re-attempting
    p = 0.1, # binding probability for one attempt
    # "zone" parameters:
    ka = 1e5, # (bulk) association rate constant [1/Ms]
    ra = 1, # radius of the association zone (w/o rMolecule) [nm]
    # in collect_stats_mode, which applies only to "zone" binding, we separate
    # the random walk and sampling phases; during rw, only stats about attempt
    # time and (if use_force = True) position are collected
    collect_stats_mode = False,
    
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
            
    def initialize_binding_zone(self, rw):
        if not self.binding or not self.bind_type == "zone":
            return
        # calculate binding rate in binding zone
        self.rbind = rw.params.rMolecule + self.ra
        Vbind = 4./3.*np.pi*(self.rbind**3 - rw.params.rMolecule**3) # [nm**3]
        Vbind *= (1e-8)**3 * nanopores.mol # [dm**3/mol = 1/M]
        kbind = 1e-9 * self.ka / Vbind # [1/ns]
        # mean no. bindings during this step
        self.nbind = kbind * rw.dt
        self.Vbind = Vbind
        self.kbind = kbind
        
    def collide(self, rw):
        "compute collisions and consequences with RandomWalk instance"
        radius = rw.params.rMolecule * self.walldist

        # determine collisions and then operate only on collided particles
        X = rw.rz if self.cyl else rw.x[rw.alive]
        if self.exclusion or (self.binding and self.bind_type == "collision"):
            collided = self.domain.inside(X, radius=radius)

        # "reflect" particles by shortening last step
        if self.exclusion:
            X0, X1 = rw.xold[rw.alive], rw.x[rw.alive]
            for i in np.nonzero(collided)[0]:
                x = self.binary_search_inside(X0[i], X1[i], radius)
                rw.update_one(i, x)

        # attempt binding for particles that can bind
        if self.binding:
            if self.bind_type == "collision":
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
            elif self.bind_type == "zone":
                rzone = self.rbind - self.domain.r if isinstance(
                        self.domain, Ball) else self.rbind
                # determine particles in binding zone
                attempt = self.domain.inside(X, radius=rzone)
                iattempt = rw.i[rw.alive][attempt]
                rw.attempt_times[iattempt] += rw.dt
                
                if not self.collect_stats_mode:
                    # draw poisson distributed number of bindings
                    bindings = np.random.poisson(self.nbind, size=np.sum(attempt))
                    # draw gamma distributed binding durations and add to time
                    duration = self.draw_zone_binding_durations(bindings, rw)
                    #duration = np.random.gamma(bindings, scale=self.t)
                    rw.bind_times[iattempt] += duration
                    # statistics
                    rw.bindings[iattempt] += bindings
                elif self.use_force:
                    self.collect_forces(attempt, iattempt, rw)
                
                # update can_bind for video
                # actually can_bind should be called can_not_bind here
                rw.can_bind[iattempt] = False
                can_bind = rw.can_bind[rw.alive]
                X_can_not_bind = X[~can_bind]
                unbind = ~self.domain.inside(X_can_not_bind, radius=rzone)
                iunbind = rw.i[rw.alive][~can_bind][unbind]
                rw.can_bind[iunbind] = True
                

    def binary_search_inside(self, x0, x1, radius):
        if self.domain.inside_single(x0, radius=radius):
            print("ERROR: x0 is in domain despite having been excluded before.")
            print("x0", x0, "x1", x1)
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
    
    def draw_zone_binding_durations(self, bindings, rw):
        if self.use_force and np.sum(bindings) > 0:
            # evaluate force magnitude at binding particles
            ibind = np.nonzero(bindings)[0]
            F = np.array([rw.F(x) for x in rw.rz[ibind]])
            F = np.sqrt(np.sum(F**2, 1))
            # create array of mean times
            kT = rw.phys.kT
            dx = 1e-9*self.dx
            t = np.zeros(bindings.shape)
            t[ibind] = self.t * np.exp(-F*dx/kT)
        else:
            t = self.t
        return np.random.gamma(bindings, scale=t)
    
    def collect_forces(self, attempt, iattempt, rw):
        if not hasattr(rw, "binding_zone_forces"):
            rw.binding_zone_forces = [[] for i in range(rw.N)]
        if np.any(attempt):
            F = np.array([rw.F(x) for x in rw.rz[attempt]])
            F = np.sqrt(np.sum(F**2, 1))
            for i, f in enumerate(F):
                rw.binding_zone_forces[iattempt[i]].append(f)

# external forces
def load_externals(**params):
    return nanopore.force_diff(**params)

class RandomWalk(object):

    def __init__(self, pore, N=10, dt=1., walldist=2.,
                 margtop=20., margbot=10., xstart=None, zstart=None,
                 rstart=None, record_positions=False, initial="sphere",
                 **params):
        # dt is timestep in nanoseconds
        self.pore = pore
        self.params = pore.params
        self.params.update(params, margtop=margtop, margbot=margbot,
                           walldist=walldist, dt=dt,
                           rstart=rstart, xstart=xstart, zstart=zstart,
                           initial=initial)
        self.sim_params = params

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
        self.attempt_times = np.zeros(N)
        self.bindings = np.zeros(N, dtype=int)

        self.domains = []
        self.add_domain(pore.protein, binding=False, exclusion=True,
                        walldist=walldist)
        if not isempty(pore.membrane):
            self.add_domain(pore.membrane, binding=False, exclusion=True,
                        walldist=walldist)
        
        self.record_positions = record_positions
        if self.record_positions:
            self.timetraces = [[] for i in range(N)]
            self.positions = [[] for i in range(N)]
            self.update_positions_record()

    # initial positions: uniformly distributed over disc
    def initial(self):
        if self.params.initial == "sphere":
            return self.initial_half_sphere()
        else:
            return self.initial_disc()
        
    def initial_disc(self):
        rstart = self.params.rstart
        xstart = self.params.xstart
        zstart = self.params.zstart
        if rstart is None:
            rstart = self.rtop - self.params.rMolecule*(self.params.walldist - 1.)
        if xstart is None:
            xstart = 0.
        if zstart is None:
            zstart = self.ztop
        self.rstart = rstart

        # create uniform polar coordinates r, theta
        r = rstart * np.sqrt(np.random.rand(self.N))
        theta = 2.*np.pi * np.random.rand(self.N)

        x = np.zeros((self.N, 3))
        x[:, 0] = xstart + r*np.cos(theta)
        x[:, 1] = r*np.sin(theta)
        x[:, 2] = zstart
        return x, np.sqrt(x[:, 0]**2 + x[:, 1]**2), x[:, 2]
    
    def initial_half_sphere(self):
        rstart = self.params.rstart
        xstart = self.params.xstart
        zstart = self.params.zstart
        if rstart is None:
            rstart = 2.*self.pore.protein.radiustop()
        if xstart is None:
            xstart = 0.
        if zstart is None:
            zstart = self.ztop + self.params.rMolecule*self.params.walldist
        self.rstart = rstart
            
        # draw 3D gaussian points, project to half-sphere and
        # only accept if above channel
        x = np.random.randn(self.N, 3)
        x[:, 2] = np.abs(x[:, 2])
        R = np.sqrt(np.sum(x**2, 1))
        m = np.array([[xstart, 0., zstart]])
        x = m + rstart*x/R[:, None]
        
        return x, np.sqrt(x[:, 0]**2 + x[:, 1]**2), x[:, 2]
        
    def add_domain(self, domain, **params):
        """add domain where particles can bind and/or are excluded from.
        domain only has to implement the .inside(x, radius) method.
        params can be domain_params"""
        dom = Domain(domain, **params)
        self.domains.append(dom)
        dom.initialize_binding_zone(self)

    def add_wall_binding(self, **params):
        dom = self.domains[0]
        dom.__dict__.update(params, binding=True)
        dom.initialize_binding_zone(self)
        
    def set_stopping_criteria(self, success=None, fail=None):
        if success is not None:
            self.is_success = success.__get__(self)
        if fail is not None:
            self.is_fail = fail.__get__(self)

    was_ood = False
    def ood_evaluation(self, f):
        dim = self.params.dim
        def newf(x):
            try:
                return f(x)
            except RuntimeError:
                if not self.was_ood:
                    print("\nFirst particle out of domain:", x)
                    self.was_ood = True
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

#    def is_success(self, r, z):
#        return (z < self.zbot - self.params.margbot) | (
#               (r > self.rbot + self.params.margbot) & (z < self.zbot))
#
#    def is_fail(self, r, z):
#        return (z > self.ztop + self.params.margtop) | (
#               (r > self.rtop + self.params.margtop) & (z > self.ztop))
        
    def is_success(self, r, z):
        return (r**2 + (z - self.zbot)**2 > self.params.margbot**2) & (
               self.below_channel(r, z))

    def is_fail(self, r, z):
        return (r**2 + (z - self.ztop)**2 > self.params.margtop**2) & (
               self.above_channel(r, z))
        
    def in_channel(self, r, z):
        if not hasattr(self, "_channel"):
            self._channel = self.pore.get_subdomain("pore")
        return self._channel.inside_winding(r, z)
    
    def above_channel(self, r, z):
        if not hasattr(self, "_above_channel"):
            self._above_channel = self.pore.get_subdomain(
                    {"bulkfluid_top", "poreregion_top"})
        return self._above_channel.inside_winding(r, z)
        
    def below_channel(self, r, z):
        if not hasattr(self, "_below_channel"):
            self._below_channel = self.pore.get_subdomain(
                    {"bulkfluid_bottom", "poreregion_bottom"})
        return self._below_channel.inside_winding(r, z)

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
        
    def update_positions_record(self):
        for i in range(self.N):
            if self.alive[i]:
                t = self.times[i] + self.bind_times[i]
                x = np.copy(self.x[i])
                self.timetraces[i].append(t)
                self.positions[i].append(x)
    

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
            
        if self.record_positions:
            self.update_positions_record()

    def walk(self):
        with nanopores.Log("Running..."):
            yield self.t
            while np.any(self.alive):
                #with nanopores.Log("%.0f ns, cpu time:" % self.t):
                self.step()
                yield self.t

        self.finalize()
        
    def draw_bindings(self, idomain=None, **bind_params):
        "(re-)draw bindings after running or loading rw"
        # get either specified domain or unique binding domain
        if idomain is None:
            domains = [dom for dom in self.domains if dom.binding]
            assert len(domains) == 1
            domain = domains[0]
        else:
            domain = self.domains[idomain]
            
        domain.__dict__.update(bind_params)
        
        # draw number of bindings
        if domain.bind_type == "zone":
            domain.initialize_binding_zone(self)
            ta = self.attempt_times
            ka = domain.kbind
            self.bindings = np.random.poisson(ka*ta)
        else:
            raise NotImplementedError("currently only for zone binding")
            
        # draw binding durations
        t = self.t
        if domain.use_force:
            raise NotImplementedError("currently no force dependency")
        self.bind_times = 1e-9*np.random.gamma(self.bindings, scale=t)
        self.times = self.walk_times + self.bind_times

    def finalize(self):
        print("finished!")
        print("mean # of attempts:", self.attempts.mean())
        tdwell = self.times.mean()
        tbind = self.bind_times.mean()
        ta = self.attempt_times.mean()
        print("mean attempt time: %s ns (fraction of total time: %s)" % (
                ta, ta/tdwell))
        print("mean # of bindings:", self.bindings.mean())
        print("mean dwell time with binding: %.3f mus"%(1e-3*(tbind + tdwell)))
        print("mean dwell time without binding: %.3f mus" % (1e-3*tdwell))
        self.walk_times = np.copy(self.times)
        self.times += self.bind_times
        
        for domain in self.domains:
            if not domain.binding or not domain.bind_type == "zone":
                continue
            # calculate effective association rate in pore
            phys = self.phys
            Dbulk = phys.DTargetBulk
            r = 1e-9*self.rstart # radius of arrival zone
            karr = 2.*self.phys.pi*r*Dbulk*1e3*phys.mol # events/Ms
            ka = domain.ka # bulk association rate [1/Ms]
            Vbind = domain.Vbind # binding volume per mole [1/M]
            cchar = 1./Vbind # [M], receptor concentration at which all targets
            # are in binding zone, so that ka * cchar = kbind
            kbind = ka * cchar # binding zone assoc. rate [1/s]
            ta_ = 1e-9*ta # mean attempt time := time in binding zone [s]
            nbind = ta_ * kbind # mean no. bindings per event
            
            keff = karr * ta_ * cchar * ka # = karr * nbind = bindings / Ms = 
            # effective association rate
            frac = karr * ta_ / Vbind # fraction of time spent in binding zone
            # in simulation, relative to bulk = keff / ka
            print("Dbulk", Dbulk)
            print("karr", karr)
            #print("ta", ta)
            #print("Vbind", Vbind)
            #print "kbind", kbind
            print("nbind: %.3f (bindings per event)" % nbind)
            print("ka [1/Ms]. Effective: %.3g, bulk: %.3g, fraction: %.3g" % (
                    keff, ka, frac))
        print('\n==============\n\n')
        
        if self.record_positions:
            for i in range(self.N):
                self.positions[i] = np.array(self.positions[i])
                self.timetraces[i] = np.array(self.timetraces[i])
                
        if hasattr(self, "binding_zone_forces"):
            for i in range(self.N):
                self.binding_zone_forces[i] = np.array(self.binding_zone_forces[i])

    def save(self, name="rw"):
        if "N" in self.params:
            self.params.pop("N")
            
        optional = dict()
        if self.record_positions:
            optional.update(positions = self.positions,
                            timetraces = self.timetraces)
        if hasattr(self, "binding_zone_forces"):
            optional.update(binding_zone_forces = self.binding_zone_forces)
        
        fields.save_fields(name, self.params,
            times = self.times,
            success = self.success,
            fail = self.fail,
            bind_times = self.bind_times,
            attempts = self.attempts,
            bindings = self.bindings,
            attempt_times = self.attempt_times,
            **optional)
        fields.update()
        
    def ellipse_collection(self, ax):
        "for matplotlib plotting"
        xz = self.x[:, [0,2]]
        #xz = self.rz
        sizes = 2.*self.params.rMolecule*np.ones(self.N)
        colors = ["b"]*self.N
        coll = collections.EllipseCollection(sizes, sizes, np.zeros_like(sizes),
                   offsets=xz, units="xy", facecolors=colors,
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
        ball_bind_zone_settings = dict(facecolor="#ffaaaa", linewidth=0.,
                                       alpha=0.5)
        patches = []
        for dom in self.domains:
            domp = dom.domain
            if isinstance(domp, Polygon):
                polygon = domp.nodes
                polygon = np.array(polygon)
                patches.append(mpatches.Polygon(polygon, **poly_settings))
                if not cyl:
                    polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])
                    patches.append(mpatches.Polygon(polygon_m, **poly_settings))
            elif isinstance(domp, Ball):
                xy = domp.x0[0], domp.x0[2]
                if dom.binding and dom.bind_type == "zone":
                    p1 = mpatches.Circle(xy, domp.r + dom.ra,
                                         **ball_bind_zone_settings)
                    p1.set_zorder(-100)
                    patches.append(p1)
                
                p = mpatches.Circle(xy, domp.r, **ball_settings)
                p.set_zorder(200)
                patches.append(p)
        return patches
    
    def plot_streamlines(self, both=False, Hbot=None, Htop=None, R=None, **params):
        R = self.params.R if R is None else R
        Htop = self.params.Htop if Htop is None else Htop
        Hbot = self.params.Hbot if Hbot is None else Hbot
        #ax = plt.axes(xlim=(-R, R), ylim=(-Hbot, Htop))
        dolfin.parameters["allow_extrapolation"] = True
        if both:
            Fel, Fdrag = fields.get_functions("force_pointsize",
                                              "Fel", "Fdrag", **self.sim_params)
            streamlines(patches=[self.polygon_patches(), self.polygon_patches()],
                        R=R, Htop=Htop, Hbot=Hbot,
                        Nx=100, Ny=100, Fel=Fel, Fdrag=Fdrag, **params)
        else:
            streamlines(patches=[self.polygon_patches()],
                        R=R, Htop=Htop, Hbot=Hbot,
                        Nx=100, Ny=100, F=self.F, **params)
        dolfin.parameters["allow_extrapolation"] = False
        
#        for p in patches:
#            p.set_zorder(100)
#            plt.gca().add_patch(p)
        plt.xlim(-R, R)
        plt.ylim(-Hbot, Htop)
        
    def plot_pore(self, cyl=False):
        R = self.params.R
        Htop = self.params.Htop
        Hbot = self.params.Hbot
        
        xlim = (-R, R) if not cyl else (0., R)
        ax = plt.axes(xlim=xlim, ylim=(-Hbot, Htop))
        patches = self.polygon_patches(cyl)
        for p in patches:
            ax.add_patch(p)
        return ax        
        
    def plot_path(self, i=0, cyl=False, **plot_params):
        self.plot_pore(cyl)
        path = self.positions[i]
        x = np.sqrt(path[:, 0]**2 + path[:, 1]**2) if cyl else path[:, 0]
        z = path[:, 2]
        plt.plot(x, z, **plot_params)
        
def setup_default(params):
    pore = get_pore(**params)
    return RandomWalk(pore, **params)
        
def _load(a):
    return a.load() if isinstance(a, fields.NpyFile) else a
    
def load_results(name, **params):
    data = fields.get_fields(name, **params) 
    data = nanopores.Params({k: _load(data[k]) for k in data})
    print("Found %d simulated events." % len(data.times))
    return data

def get_results(name, params, setup=setup_default, calc=True):
    # setup is function rw = setup(params) that sets up rw
    # check existing saved rws
    data = None
    if fields.exists(name, **params):
        data = load_results(name, **params)
        N = len(data.times)
    else:
        N = 0
    # determine number of missing rws and run
    N_missing = params["N"] - N
    if N_missing > 0 and calc:
        new_params = nanopores.Params(params, N=N_missing)
        rw = setup(new_params)
        run(rw, name)
        rw.save(name)
        data = load_results(name, **params)
    # return results
    elif data is None:
        data = load_results(name, **params)
    return data

def reconstruct_rw(data, params, setup=setup_default, finalize=True,
                   **setup_args):
    """if positions were recorded, the complete random walk instance
    can IN PRINCIPLE be reconstructed and restarted.
    this is a simple first attempt."""
    # TODO: recreate rw.rz, rw.t, rw.can_bind    
    rw = setup(params, **setup_args)
    rw.__dict__.update(
        times = data.times - data.bind_times,
        success = data.success,
        fail = data.fail,
        bind_times = data.bind_times,
        attempts = data.attempts,
        bindings = data.bindings,
        attempt_times = data.attempt_times,
    )
    if rw.record_positions:
        rw.positions = [[x for x in _load(X)] for X in data.positions]
        rw.timetraces = [[t for t in _load(T)] for T in data.timetraces]
        rw.x = np.array([X[-1] for X in data.positions])
        rw.alive = ~(data.success | data.fail)
    if hasattr(data, "binding_zone_forces"):
        rw.binding_zone_forces = [list(_load(X)) for X in data.binding_zone_forces]
        
    if finalize:
        rw.finalize()
    return rw
    
def get_rw(name, params, setup=setup_default, calc=True, finalize=True,
           **setup_args):
    # setup_args are only relevant to reconstructed rw
    # everything that influences data have to be in params
    data = get_results(name, params, setup, calc)
    return reconstruct_rw(data, params, setup, finalize, **setup_args)

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

def histogram(rw, a=-3, b=3, scale=1e-0):
    t = rw.times * 1e-9 / scale # assuming times are in nanosaconds

    exponential_hist(t[rw.success], a, b, color="green", label="translocated")
    exponential_hist(t[rw.fail], a, b, color="red", label="did not translocate")

    plt.xlabel(r"$\tau$ off [s]")
    plt.ylabel("count")
    plt.legend(loc="best")

def hist_poisson(rw, name="attempts", ran=None, n=10, pfit=True, mpfit=True, lines=True):
    attempts = getattr(rw, name)
    if ran is None:
        n0 = 0
        n1 = n
    else:
        n0, n1 = ran

    k = np.arange(n0, n1 + 1)
    bins = np.arange(n0 - 0.5, n1 + 1.5)
    
    astr = "ap = %.3f" if name == "bindings" else "a = %.1f"

    plt.hist(attempts, bins=bins, label="Simulated "+name, color="#aaaaff")
    # poisson fit
    a = attempts.mean()
    K = len(attempts)
    a0 = attempts[attempts >= 1].mean()
    a1 = poisson_from_positiveK(a0)
    print("Mod. Poisson fit, mean of K>0:", a0)
    print("Inferred total mean:", a1)
    print("Standard Poisson fit, mean:", a)
    p1 = a/a1
    K1 = len(attempts[attempts > 0])/(1.-np.exp(-a1))
    
    pdf = K*poisson.pmf(k, a)
    pdf1 = K*p1*poisson.pmf(k, a1)
    if n0 == 0:
        pdf1[0] += K*(1. - p1)
        
    k0 = np.linspace(n0, n1, 500)
    if pfit:
        if lines:
            plt.plot(k0, K*gamma.pdf(a, k0 + 1), "-", color="C1")        
        plt.plot(k, pdf, "s",
                 label=("Poisson fit, "+astr)%(a), color="C1")
    if mpfit:
        if lines:
            plt.plot(k0, K1*gamma.pdf(a1, k0 + 1), "-", color="C2")
        plt.plot(k, pdf1, "v",
                 label=("Mod. Poisson fit, "+astr)%(a1), color="C2")
    plt.xlim(n0 - 0.5, n1 + 0.5)
    plt.xticks(k, k)
    plt.yscale("log")
    plt.ylim(ymin=1.)
    plt.xlabel("# %s" % name)
    plt.ylabel("Count")
    plt.legend()

def solve_newton(C, f, f1, x0=1., n=20):
    "solve f(x) == C"
    x = x0 # initial value
    print("Newton iteration:")
    for i in range(n):
        dx = -(f(x) - C)/f1(x)
        x = x + dx
        res = abs(f(x) - C)
        print(i, "Residual", res, "Value", x)
        if res < 1e-12:
            break
    print()
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
                 extra_args=["-vcodec", "libx264"])

# convenience function for interactive experiments
def run(rw, name="rw", plot=False, a=-3, b=3, **aniparams):
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
        #plt.figure()
        #hist_poisson(rw, "attempts")
        #plt.figure()
        #hist_poisson(rw, "bindings")
        plt.show()

if __name__ == "__main__":
    pore = get_pore(**params)
    rw = RandomWalk(pore, **params)
    receptor = Ball([15., 0., 30.], 8.)
    rw.add_domain(receptor, exclusion=True, walldist=1.,
                  binding=True, eps=1., t=1.5e6, p=0.14)
    run(rw, name="wei")
