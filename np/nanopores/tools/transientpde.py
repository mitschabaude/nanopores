""" PDE class for time-dependent systems """

from dolfin import plot
from .pdesystem import PDESystem
from .illposed import IllposedLinearSolver
import matplotlib.pyplot as plt
import numpy as np

__all__ = ["TransientLinearPDE"]

def timerange(T, dt):
    t = dt
    while t <= T:
        yield t
        t += dt

class TransientLinearPDE(PDESystem):
    dt = 1 # default time step [s]
    
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
        
    def timestep(self, **params):
        PDESystem.solve(self, verbose=False, **params)
        
    def solve(self, t=0, verbose=True, visualize=True, **params):
        if verbose:
            print "\n"
        # evolve system up to time t
        for t_ in timerange(t, self.dt):
            self.timestep(**params)
            if verbose:
                print "\x1b[A","\r",
                print "t = %s [s]" %t_
            if visualize:
                self.visualize(t_)
        if visualize:
            plt.show(block=True)
            
    def visualize(self, t):
        u = self.solution
        #plot(u, title="solution at time %s [s]" %t)
        
        if t < 1.5*self.dt: # first step
            # FIXME: softcode this
            #self.x = np.array([[0., 0., float(z)] for z in range(-15, 10, 1)])
            self.x = np.linspace(-25., 60., 200)
            plt.ion()
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
        else:
            self.ax.clear()
            
        # TODO: efficient?
        y = np.array([u([0., 0., z]) for z in self.x])
        self.ax.plot(self.x, y, 'r-')
        self.ax.set_title("t = %s [s]" %t)
        self.ax.set_xlabel("z [nm]")
        self.ax.set_ylabel("p(x, t)")
        self.ax.set_ylim([0.,1.])
        #self.ax.draw()
        #self.line.set_ydata(y)
        self.fig.canvas.draw()
        #plt.plot(self.x, y)
          
        
