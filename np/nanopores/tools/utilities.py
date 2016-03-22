''' some utility functions for global use after importing * from nanopores '''

from importlib import import_module
import inspect, os, sys
import matplotlib.pyplot as plt
import numpy as np
import dolfin
from nanopores.dirnames import DATADIR
from nanopores.tools.protocol import unique_id

__all__ = ["import_vars", "get_mesh", "u_to_matlab", "plot_on_sub", "save_dict", "plot_sliced",
           "crange", "plot1D", "showplots", "saveplots", "loadplots", "add_params",
           "plot_cross", "plot_cross_vector", "load_dict"]

def crange(a, b, N): # continuous range with step 1/N
    return [x/float(N) for x in range(a*N, b*N+1)]

def import_vars(mod):
    d = vars(import_module(mod))
    return {k:d[k] for k in d if not k.startswith("_")}
    
def get_mesh(geo_name, mesh_name = "mesh/mesh.xml"):

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
    u = adaptfunction(u, submesh, assign=False, interpolate=True)
    u0 = u if expr is None else expr
    plot(u0, title=title)
    
def plot_sliced(geo):
    tol = 1e-5
    class Back(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return x[1] >= -tol
    back = dolfin.CellFunction("size_t", geo.mesh, 0)
    Back().mark(back, 1)
    submesh = dolfin.SubMesh(geo.mesh, back, 1)
    #plot(submesh)
    bb = geo.mesh.bounding_box_tree()
    subsub = dolfin.CellFunction("size_t", submesh, 0)
    sub = geo.subdomains
    for cell in dolfin.cells(submesh):
        iparent = bb.compute_first_entity_collision(cell.midpoint())
        subsub[cell] = sub[int(iparent)]
    dolfin.plot(subsub)
    
class uCross(dolfin.Expression):
    def __init__(self, u, axis=1):
        self.u = u
        self.i = axis
        dolfin.Expression.__init__(self)
    def eval(self, value, x):
        y = list(x)
        y.insert(self.i, 0.)
        value[0] = self.u(y)
        
class uCrossVector(dolfin.Expression):
    def __init__(self, u, axis=1):
        self.u = u
        self.i = axis
        dolfin.Expression.__init__(self)
    def eval(self, value, x):
        y = list(x)
        i = self.i
        y.insert(i, 0.)
        ux = self.u(y)
        value[0] = ux[0] # for axis=1
        value[1] = ux[2]
    def value_shape(self):
        return (2,)

def plot_cross(u, mesh2D, title="", axis=1):
    # create Expression to evaluate u on a hyperplane
    ucross = uCross(u=u, axis=axis)
    # interpolate u onto the 2D mesh
    V = dolfin.FunctionSpace(mesh2D, "CG", 1)
    u2D = dolfin.Function(V)
    u2D.interpolate(ucross)
    dolfin.plot(u2D, title=title)
    
def plot_cross_vector(u, mesh2D, title="", axis=1):
    # create Expression to evaluate u on a hyperplane
    ucross = uCrossVector(u=u, axis=axis)
    # interpolate u onto the 2D mesh
    V = dolfin.VectorFunctionSpace(mesh2D, "CG", 1)
    u2D = dolfin.Function(V)
    u2D.interpolate(ucross)
    dolfin.plot(u2D, title=title)
    
def save_dict(data, dir=".", name="file"):
    # works for all data that can be recovered from their repr()
    if not os.path.exists(dir):
        os.makedirs(dir)
    with open(os.path.join(dir, name + ".txt"), 'w') as f:
        f.write(repr(data))
        
def load_dict(dir, name):
    from numpy import array
    fname = os.path.join(dir, name + ".txt")
    with open(fname, 'r') as f:
        data = f.read()
    return eval("".join(data.split("\n")))
        
def _call(f, params):
    # call f without knowing its arguments --
    # they just have to be subset of dict params.
    argnames = inspect.getargspec(f).args
    args = {k: params[k] for k in argnames if k in params}
    return f(**args)
    
def plot1D(functions, rng=(0.,1.,101), axis=None, dim=3, axlabels=("",""),
        line=[1.,0.,0.], origin=[0.,0.,0.], plot="plot", title="", legend="upper right",
        show=False, newfig=True):
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
    if newfig:
        plt.figure()
    #ax = fig.add_subplot(111)
    lines = []
    for fstr, f in functions.items():
        y = np.array([f(list(t)) for t in z])
        lines.append(getattr(plt, plot)(x, y, "-x", label=fstr))
        plt.xlabel(axlabels[0])
        plt.ylabel(axlabels[1])
        plt.title(title)
        plt.legend(loc=legend)
        #ax.set_xlabel(axlabels[0])
        #ax.set_ylabel(axlabels[1])
        #ax.set_title(title)
        #ax.legend(loc=legend)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if show: plt.show()
    #return lines
    
def showplots():
    plt.show()
    
def saveplots(name="plot", meta=None, uid=False):
    # collect data from every open figure
    plots = []
    figs = map(plt.figure, plt.get_fignums())
    for fig in figs:
        for ax in fig.axes:
            plot = dict(
                xlabel = ax.get_xlabel(),
                ylabel = ax.get_ylabel(),
                xscale = ax.get_xscale(),
                yscale = ax.get_yscale(),
                data = [])
            for line in ax.lines:
                x, y = line.get_data()
                marker = line.get_marker()
                if marker == "None":
                    marker = ""
                data = dict(
                    x = x,
                    y = y,
                    style = marker + line.get_linestyle(),
                    label = line.get_label())
                plot["data"].append(data)
            plots.append(plot)
    # add metadata if provided
    meta = {} if meta is None else meta
    data = dict(plots=plots, meta=meta)
    # save to txt file in DATADIR
    DIR = os.path.join(DATADIR, "plots")
    name = name + "_" + str(unique_id()) if uid else name
    save_dict(data, DIR, name)
    return plots
    
def loadplots(name, show=True):
    defaultstyle = "-x"
    DIR = os.path.join(DATADIR, "plots")
    dic = load_dict(DIR, name)
    meta = dic["meta"]
    plots = dic["plots"]
    for plot in plots:
        plt.figure()
        for line in plot["data"]:
            style = line["style"] if "style" in line else defaultstyle
            plt.plot(line["x"], line["y"], style, label=line["label"])
        plt.xlabel(plot["xlabel"])
        plt.ylabel(plot["ylabel"])
        if "xscale" in plot:
            plt.xscale(plot["xscale"])
        if "yscale" in plot:
            plt.yscale(plot["yscale"])
        plt.legend()
    if show: plt.show()
    return meta
        
    
def add_params(**params):
    # this is some magic to attach a set of parameters to a module
    # but be able to use them in the same module without any ugly complication
    # TODO: connect with parsed command line args
    args = _argparse()
    params.update({key: args[key] for key in args if key in params})
    frm = inspect.stack()[1]
    mod = inspect.getmodule(frm[0])
    mod.__dict__.update(params)
    if not hasattr(mod, "PARAMS"):
        mod.PARAMS = dict()
    mod.PARAMS.update(params)

def _argparse():
    # parse arguments only of the form 'arg1 val1 arg2 val2 ...'
    # the values are tried to be converted with eval
    dic = {}
    for name, val in zip(sys.argv[1::2], sys.argv[2::2]):
        try:
            val = eval(val)
        except NameError:
            pass
        dic[name] = val
    return dic
            
    
