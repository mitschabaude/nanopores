""" PDE class for time-dependent systems """

from dolfin import plot
from .pdesystem import PDESystem
from .illposed import IllposedLinearSolver
import matplotlib.pyplot as pyplot

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
        
    def timestep(self, t, verbose=True, visualize=False, **params):
        PDESystem.solve(self, verbose=False, **params)
        if verbose:
            print "Time:",t
        if visualize:
            self.visualize(t)
        
    def solve(self, t=0, **params):
        # evolve system up to time t
        for t_ in timerange(t, self.dt):
            self.timestep(t_, **params)
            
    def visualize(self, t):
        plot(self.solution, title="solution at time %s" %t)
        '''
        x = range(-15, 5, 1)
        
        pyplot.ion()

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        line1, = ax.plot(x, y, 'r-') # Returns a tuple of line objects, thus the comma

        for phase in np.linspace(0, 10*np.pi, 500):
            line1.set_ydata(np.sin(x + phase))
            fig.canvas.draw()
        
        
            
        sol = self.solutions(deepcopy=True)
        mesh = self.geo.mesh
        on = ""

        if subdomain:
            on = " on " + subdomain
            mesh = self.geo.submesh(subdomain)
            for f in sol:
                adaptfunction(f,mesh,assign=True)

        plot(mesh, title="final mesh"+on)
        for f in sol:
            plot(f, title = str(f)+on)
        interactive()
        '''    
        
