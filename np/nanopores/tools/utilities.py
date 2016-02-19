''' some utility functions for global use after importing * from nanopores '''

from importlib import import_module
import inspect
import matplotlib.pyplot as plt
import numpy as np

__all__ = ["import_vars", "get_mesh", "u_to_matlab", "plot_on_sub", "save_dict",
           "crange", "plot1D", "showplots"]

def crange(a, b, N): # continuous range with step 1/N
    return [x/float(N) for x in range(a*N, b*N+1)]

def import_vars(mod):
    d = vars(import_module(mod))
    return {k:d[k] for k in d if not k.startswith("_")}
    
def get_mesh(geo_name, mesh_name = "mesh/mesh.xml"):
    from nanopores import DATADIR
    from dolfin import Mesh
    return Mesh("/".join([DATADIR, geo_name, mesh_name]))
    
def u_to_matlab(mesh, u, name="u"):
    # save real-valued P1 Function u at the mesh nodes to matlab arrays X, U
    # X = mesh.coordinates() = [x1_1, .., x1_d; ...; xN_1, .., xN_d]
    # U = u(X)
    from dolfin import vertex_to_dof_map
    from numpy import array
    X = mesh.coordinates()
    v2d = vertex_to_dof_map(u.function_space())
    U = u.vector()[v2d]
    from scipy.io import savemat
    dic = {"X": X, "U": U[:,None]}
    savemat("%s.mat" %name, dic)
    
def plot_on_sub(u, geo, sub, expr=None, title=""):
    from nanopores.tools.illposed import adaptfunction
    from dolfin import plot
    submesh = geo.submesh(sub)
    adaptfunction(u, submesh, assign=True)
    u0 = u if expr is None else expr
    plot(u0, title=title)
    
def save_dict(data, dir=".", name="file"):
    with open('%s/%s.txt' % (dir,name), 'w') as f:
        f.write(repr(data))
        
def _call(f, params):
    # call f without knowing its arguments --
    # they just have to be subset of dict params.
    argnames = inspect.getargspec(f).args
    args = {k: params[k] for k in argnames if k in params}
    return f(**args)
    
def plot1D(functions, rng=(0.,1.,101), axis=None, dim=3, axlabels=("",""),
        line=[1.,0.,0.], origin=[0.,0.,0.], plot="plot", title="", legend="upper right",
        show=False):
    # functions is dict(string = callable)
    # output of functions should be real-valued
    
    # create array of points
    x = np.linspace(*rng)
    if axis is not None:
        index = dict(x=0, y=1, z=2)
        line = [0. for i in range(dim)]
        line[index[axis]] = 1.
    origin = np.array(origin)[:dim]
    line = np.array(line)[:dim]
    z = np.array([t*line + origin for t in x]) # could also be done in a numpy way
    
    # plot functions
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for fstr, f in functions.items():
        y = np.array([f(list(t)) for t in z])
        getattr(ax, plot)(x, y, "-x", label=fstr)
        ax.set_xlabel(axlabels[0])
        ax.set_ylabel(axlabels[1])
        ax.set_title(title)
        ax.legend(loc=legend)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if show: plt.show()
    return fig
    
def showplots():
    plt.show()
    
