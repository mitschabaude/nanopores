# (c) 2017 Gregor Mitscha-Baude
"random walk of many particles in cylindrical pore"
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections
import matplotlib.patches as mpatches

import dolfin
import nanopores
from nanopores.tools import fields
from nanopores.geometries.allpores import get_pore
fields.set_dir_dropbox()
dolfin.parameters["allow_extrapolation"] = False
params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 6.,
    Qmol = -1.,
    bV = -0.5,
    # random walk params
    N = 10, # number of (simultaneous) random walks
    dt = 1.,
)

# external forces
def load_externals(**params):
    F, = fields.get_functions("wei_force_ps", "F", **params)
    D, = fields.get_functions("wei_D_2D", "D", **params)
    V = D.function_space()
    divD = dolfin.project(dolfin.as_vector([
               dolfin.grad(D[0])[0], dolfin.grad(D[1])[1]]), V)
    return F, D, divD

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

    def __init__(self, pore, N=10, dt=1., **params):
        # dt is timestep in nanoseconds
        self.pore = pore
        self.params = pore.params
        r1 = pore.protein.radiustop() - self.params.rMolecule - 10.
        ztop = pore.protein.zmax()[1]
        x, r, z = initial(r1, ztop, N)
        self.N = N
        self.x = x
        self.xold = x
        self.rz = np.column_stack([r, z])
        self.dt = dt
        
        self.t = 0.
        F, D, divD = load_externals(**params)
        self.F = F #self.ood_evaluation(F)
        self.D = D #self.ood_evaluation(D)
        self.divD = divD #self.ood_evaluation(divD)
        self.phys = nanopores.Physics("pore_mol", **params)
        
        self.alive = np.full((N,), True, dtype=bool)
        self.success = np.full((N,), False, dtype=bool)
        self.fail = np.full((N,), False, dtype=bool)

        #self._Dx = np.zeros((N, params["dim"]))
        #self._Fx = np.zeros((N, params["dim"]))
        
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
        self.finish()
        r = np.sqrt(np.sum(self.x[self.alive, :2]**2, 1))
        self.rz = np.column_stack([r, self.x[self.alive, 2]])
        self.t += self.dt
        
    def update_one(self, i, xnew):
        self.x[np.nonzero(self.alive)[0][i]] = xnew
        self.rz[i, 0] = np.sqrt(xnew[0]**2 + xnew[1]**2)
        self.rz[i, 1] = xnew[2]
        
    def inside_wall(self, factor=1., x=None):
        if x is None:
            return self.pore.protein.inside(self.rz,
                   radius=self.params.rMolecule*factor)
        else:
            return self.pore.protein.inside_single(x,
                   radius=self.params.rMolecule*factor)

    def simple_reflect(self, factor=2., minsize=0.5):
        radius = self.params.rMolecule*factor
        inside = self.inside_wall(factor)
        alive = self.alive
        X0, X1 = self.xold[alive], self.x[alive]
        for i in np.nonzero(inside)[0]:
            x0, x1 = X0[i], X1[i]
            x = self.binary_search_inside(x0, x1, radius, minsize)
            #t = (x - x0)/(x1 - x0)
            #print "update", t
            #print "update", self.xold[i], ",", x, ",", self.x[i]
            self.update_one(i, x)
            
    def binary_search_inside(self, x0, x1, radius, minsize=0.5):
        if self.pore.protein.inside_single(x0, radius=radius):
            print "ERROR: something wrong"
            print x0, x1
            raise Exception
        if np.sum((x0 - x1)**2) < minsize**2:
            return x0
        x05 = .5*(x0 + x1)
        if self.pore.protein.inside_single(x05, radius=radius):
            x1 = x05
        else:
            x0 = x05
        return self.binary_search_inside(x0, x1, radius, minsize)
    
    def finish(self):
        alive = self.alive
        z = self.x[alive, 2]
        r = np.sqrt(np.sum(self.x[alive, :2]**2, 1))
        Htop = self.params.Htop - 1.
        Hbot = self.params.Hbot - 1.
        R = self.params.R - 1.
        zmem = self.params.zmem
        self.success[alive] = (z < -Hbot) | ((r > R) & (z < zmem))
        self.fail[alive] = (z > Htop) | ((r > R) & (z > zmem))
        self.alive[alive] = ~(self.fail[alive] | self.success[alive])
        #print self.alive

    def walk(self):
        "one step of random walk"
        # evaluate F and D
        D = self.evaluate_D_cyl()
        F = self.evaluate_vector_cyl(self.F)
        divD = 1e9*self.evaluate_vector_cyl(self.divD)
        kT = self.phys.kT
        dt = self.dt
        
        # get step
        dW = self.brownian(D)
        dx = dW + dt*divD + dt*D/kT*F*0.5
        #print "%.2f (dx) = %.2f (dW) + %.2f (divD) + %.2f (F)" % (
        #        abs(dx[0, 2]), abs(dW[0, 2]), abs(dt*divD[0, 2]), abs((dt*D/kT*F)[0, 2]))
        #print ("t = %.2f microsec" % (self.t*1e-3))
        
        # update position and time and discard finished particles
        self.update(dx)
        
        # correct particles that collided with pore wall
        self.simple_reflect(factor=1.8, minsize=0.01)
        
        
#        inside = self.inside_wall(factor=2.)
#        if any(inside):
#            print self.x[inside]
#            print [self.inside_wall(2., x) for x in self.x[inside]]
#            raise Exception("ERROR")

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

    def move_ellipses(self, coll):
        xz = self.x[:, [0,2]]
        #xz = self.rz
        coll.set_offsets(xz)
        #inside = self.inside_wall()
        margin = np.nonzero(self.alive)[0][self.inside_wall(2.)]
        colors = np.full((self.N,), "b", dtype=str)
        colors[margin] = "g"
        colors[self.success] = "w"
        colors[self.fail] = "k"
        #colors = [("r" if inside[i] else "g") if margin[i] else "b" for i in range(self.N)]
        coll.set_facecolors(colors)
        #y = self.x[:, 1]
        #d = 50.
        #sizes = self.params.rMolecule*(1. + y/d)
        #coll.set(widths=sizes, heights=sizes)
        
        
def polygon_patches(rw, ax):
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=1.,
                    edgecolor="black")
    polygon = rw.pore.protein.nodes
    polygon = np.array(polygon)
    polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    patch1 = mpatches.Polygon(polygon, **settings)
    patch2 = mpatches.Polygon(polygon_m, **settings)
    #patch.set_zorder(10)
    #patchm.set_zorder(10)
    return patch1, patch2
    #ax.add_patch(patch)
    #ax.add_patch(patchm)

def panimate(rw, **aniparams):
    R = rw.params.R
    Htop = rw.params.Htop
    Hbot = rw.params.Hbot

    fig = plt.figure()
    fig.set_size_inches(6, 6)
    #ax = plt.axes([0,0,1,1], autoscale_on=False, xlim=(-R, R), ylim=(-H, H))
    ax = plt.axes(xlim=(-R, R), ylim=(-Hbot, Htop))
    coll = rw.ellipse_collection(ax)
    patch1, patch2 = polygon_patches(rw, ax)

    def init():
        return ()

    def animate(i):
        if i == 0:
            ax.add_collection(coll)
            ax.add_patch(patch1)
            ax.add_patch(patch2)
        else:
            rw.walk()
            rw.move_ellipses(coll)
        return coll, patch1, patch2

    aniparams = dict(dict(frames=1800, interval=10, blit=True), **aniparams)
    ani = animation.FuncAnimation(ax.figure, animate, init_func=init, **aniparams)
    return ani

if __name__ == "__main__":
    pore = get_pore(**params)
    rw = RandomWalk(pore, **params)
    
    #while rw.t < 1e4:
    #    rw.walk()
    
    ani = panimate(rw)
    plt.show()