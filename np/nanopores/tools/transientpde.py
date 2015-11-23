""" PDE class for time-dependent systems """

from dolfin import plot
from .pdesystem import PDESystem
from .illposed import IllposedLinearSolver
import matplotlib.pyplot as plt
import numpy as np

__all__ = ["TransientLinearPDE", "timerange", "logtimerange", "TimeDependentPlotter"]

def timerange(T, dt):
    t = dt
    while t <= T:
        yield t
        t += dt
        
def logtimerange(T0, levels, frac=0.01, mult=10, change_dt=None):
    t = 0
    for l in range(levels):
        dt = T0*frac
        if change_dt is not None:
            change_dt(dt)
        while t < T0:
            t += dt
            yield t
        T0 *= mult

class TransientLinearPDE(PDESystem):
    dt = 1 # default time step [s]
    time = []
    
    def __init__(self, Problem, geo=None, phys=None,
                 dt=None, **problem_params):
        if dt is not None:
            self.dt = dt
        
        # convention: Problem sets solution to initial condition
        # and uses it at previous timestep in the forms
        # i don't think there's any need for storing two functions
        problem = Problem(geo=geo, phys=phys, dt=dt, **problem_params)
        solver = IllposedLinearSolver(problem)

        self.geo = geo
        self.functions = {Problem.__name__: problem.solution()}
        self.solution = problem.solution()
        self.solvers = {Problem.__name__: solver}
        self.solver = solver
        self.functionals = {}
        
    def change_dt(self, dt):
        self.dt = dt
        # needs problem.update_forms
        for solver in self.solvers.values():
            solver.problem.update_forms(dt=dt)
            solver.assemble_A()
        
    def timestep(self, **params):
        PDESystem.solve(self, verbose=False, **params)
        
    def solve(self, t=0, verbose=True, visualize=False, record_functionals=True, **params):
        if not hasattr(self, "timerange"):
            self.timerange = timerange(t, self.dt)
    
        if verbose:
            print "\n"
        # evolve system up to time t
        for t_ in self.timerange:
            self.time.append(t_)
            self.timestep(**params)
            if verbose:
                print "\x1b[A","\r",
                print "t = %s [s]" %t_
            if record_functionals:
                self.record_functionals()
            if visualize:
                self.visualize()
                
        if visualize:
            self.finish_plots()
            
    def visualize(self):
        u = self.solution
        plot(u, title="solution at time %s [s]" %(self.time[-1],))
        
    def finish_plots(self):
        interactive()
        
    def record_functionals(self):
        for J in self.functionals.values():
            J.evaluate()
            
    def plot_functionals(self, plot="plot", title=""):
        plt.figure()
        for Jstr, J in self.functionals.items():
            getattr(plt, plot)(self.time, J.values, "-x", label=Jstr)
            plt.xlabel("time [s]")
            plt.title(title)
            plt.legend(loc="lower right")
            #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show(block=True)
        
class TimeDependentPlotter(object):
    "plot a function that changes with time over 1D range."
    
    title = "t = %s [s]"
    xlabel = "z [nm]"
    ylabel = "p(x, t)"
    #ylim = [0.,1.]
    style = "-x"
    
    def __init__(self, ft, xran, dt):
        # ft is a function that is supposed to change with t
        # xran = (a, b, step)
        plt.ion()
        self.x = np.linspace(xran[0], xran[1], xran[2])
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ft = ft
        self.dt = dt
        self.t = 0.
        self.isclear = True
        
    def plot(self, t=None):
        self.t += self.dt
        if t is not None:
            self.t = t
        if not self.isclear:
            self.ax.clear()
        y = np.array([self.ft(z) for z in self.x])
        self.ax.plot(self.x, y, self.style)
        self.ax.set_title(self.title %(self.t,))
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        #self.ax.set_ylim(self.ylim)
        self.fig.canvas.draw()
        self.isclear = False
        
    def finish(self):
        plt.show(block=True)
          
        
