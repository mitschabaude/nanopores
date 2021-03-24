''' some utility functions for global use after importing * from nanopores '''

from importlib import import_module
import inspect, os, sys, glob, time
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
import json
import dolfin
from nanopores.dirnames import DATADIR
from nanopores.tools.protocol import unique_id
from functools import reduce

__all__ = ["import_vars", "get_mesh", "u_to_matlab", "plot_on_sub", "save_dict", "plot_sliced",
           "crange", "plot1D", "showplots", "saveplots", "loadplots", "add_params",
           "plot_cross", "plot_cross_vector", "load_dict", "save_stuff", "load_stuff",
           "save_functions", "load_functions", "load_vector_functions", "load_mesh",
           "convert3D", "convert2D", "RectangleMesh", "savefigs", "Params",
           "user_params", "user_param", "dict_union", "union", "plot_sliced_mesh",
           "smooth", "collect", "collect_dict", "assertdir", "any_params", "tic", "toc"]

def crange(a, b, N): # continuous range with step 1/N
    return [x/float(N) for x in range(a*N, b*N+1)]

def smooth(a, k=3):
    a = np.array(a)
    # range of kernel
    start = -(k // 2)
    end = (k // 2) + (k % 2)
    N = a.shape[0]
    b = [a[max(0, i + start) : min(N, i + end)].mean() for i in range(N)]
    return np.array(b)

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
    X = mesh.coordinates()
    v2d = dolfin.vertex_to_dof_map(u.function_space())
    U = u.vector()[v2d]
    from scipy.io import savemat
    dic = {"X": X, "U": U[:,None]}
    savemat("%s.mat" %name, dic)

def plot_on_sub(u, geo, sub, expr=None, title=""):
    from nanopores.tools.illposed import adaptfunction
    submesh = geo.submesh(sub)
    u = adaptfunction(u, submesh, assign=False, interpolate=True)
    u0 = u if expr is None else expr
    dolfin.plot(u0, title=title)

def plot_sliced(geo, **params):
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
    plot_params = dict(title="sliced geometry with subdomains",
                       elevate=-90., **params)
    dolfin.plot(subsub, **plot_params)

def plot_sliced_mesh(geo, **kwargs):
    tol = 1e-5
    class Back(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return x[1] >= -tol
    back = dolfin.CellFunction("size_t", geo.mesh, 0)
    Back().mark(back, 1)
    submesh = dolfin.SubMesh(geo.mesh, back, 1)
    dolfin.plot(submesh, **kwargs)

class uCross(dolfin.Expression):
    def __init__(self, u, axis=1, **kw):
        self.u = u
        self.i = axis
    def eval(self, value, x):
        y = list(x)
        y.insert(self.i, 0.)
        value[0] = self.u(y)

class uCrossVector(dolfin.Expression):
    def __init__(self, u, axis=1, **kw):
        self.u = u
        self.i = axis
    def eval(self, value, x):
        y = list(x)
        i = self.i
        y.insert(i, 0.)
        ux = self.u(y)
        value[0] = ux[0] # for axis=1
        value[1] = ux[2]
    def value_shape(self):
        return (2,)

def plot_cross(u, mesh2D, axis=1, **kwargs):
    # create Expression to evaluate u on a hyperplane
    ucross = uCross(u=u, axis=axis, degree=1)
    # interpolate u onto the 2D mesh
    V = dolfin.FunctionSpace(mesh2D, "CG", 1)
    u2D = dolfin.Function(V)
    u2D.interpolate(ucross)
    dolfin.plot(u2D, **kwargs)

def plot_cross_vector(u, mesh2D, axis=1, **kwargs):
    # create Expression to evaluate u on a hyperplane
    ucross = uCrossVector(u=u, axis=axis, degree=1)
    # interpolate u onto the 2D mesh
    V = dolfin.VectorFunctionSpace(mesh2D, "CG", 1)
    u2D = dolfin.Function(V)
    u2D.interpolate(ucross)
    dolfin.plot(u2D, **kwargs)

def save_dict(data, dir=".", name="file"):
    # works for all data that can be recovered from their repr()
    if not os.path.exists(dir):
        os.makedirs(dir)
    with open(os.path.join(dir, name + ".txt"), 'w') as f:
        f.write(repr(data))

def load_dict(dir, name):
    fname = os.path.join(dir, name + ".txt")
    with open(fname, 'r') as f:
        data = f.read()
    return eval("".join(data.split("\n")))

def _open(name, folder, mode):
    DIR = os.path.join(DATADIR, folder)
    if not os.path.exists(DIR):
        os.makedirs(DIR)
    FILE = os.path.join(DIR, name + ".txt")
    return open(FILE, mode)

def save_stuff(name, *stuff):
    with _open(name, "stuff", "w") as f:
        json.dump(stuff, f)

def load_stuff(name):
    # Stuff API:
    # save_stuff("name", stuff) --> stuff = load_stuff("name")
    # save_stuff("name", stuff1, stuff2) --> stuff1, stuff2 = load_stuff("name")
    try:
        file_ = _open(name, "stuff", "r")
    except IOError:
        raise IOError("nanopores.load_stuff: Nothing named '%s' yet exists." % name)
    with file_ as f:
        stuff = json.load(f)
    if len(stuff) == 1:
        return stuff[0]
    else:
        return tuple(stuff)

def save_functions(name, mesh, meta=None, DIR=None, **functions):
    # save dolfin functions and mesh
    if DIR is None:
        DIR = os.path.join(DATADIR, "functions", "")
    if not os.path.exists(DIR):
        os.makedirs(DIR)
    dolfin.File(DIR + name + "_mesh.xml") << mesh
    for fname, f in list(functions.items()):
        assert f.function_space().mesh().id() == mesh.id()
        dolfin.File(DIR + name + "_" + fname + ".xml") << f
    with _open(name + "_meta", "functions", "w") as f:
        json.dump(meta, f)

def _find_names_in_files(pre, post):
    return [string[len(pre):-len(post)] for string in glob.glob(pre + "*" + post)]

def load_functions(name, space=None, DIR=None):
    # load dolfin functions with pattern matching approach
    # space is a lambda with input mesh
    if DIR is None:
        DIR = os.path.join(DATADIR, "functions", "")
    mesh = dolfin.Mesh(DIR + name + "_mesh.xml")
    if space is None:
        space = lambda mesh: dolfin.FunctionSpace(mesh, "CG", 1)
    V = space(mesh)
    functions = {}
    fnames = _find_names_in_files(DIR + name + "_", ".xml")
    for fname in fnames:
        if fname == "mesh": continue
        #print "name:", fname
        #print "file:", DIR + name + "_" + fname + ".xml"
        f = dolfin.Function(V, DIR + name + "_" + fname + ".xml")
        functions[fname] = f
    with _open(name + "_meta", "functions", "r") as f:
        meta = json.load(f)
    return functions, mesh, meta

def load_mesh(name, DIR=None):
    if DIR is None:
        DIR = os.path.join(DATADIR, "functions", "")
    mesh = dolfin.Mesh(DIR + name + "_mesh.xml")
    return mesh

def load_vector_functions(name, DIR=None):
    def space(mesh):
        return dolfin.VectorFunctionSpace(mesh, "CG", 1)
    return load_functions(name, space, DIR)

def _call(f, params):
    # call f without knowing its arguments --
    # they just have to be subset of dict params.
    argnames = inspect.getargspec(f).args
    args = {k: params[k] for k in argnames if k in params}
    return f(**args)

def plot1D(functions, rng=(0.,1.,101), axis=None, dim=3, axlabels=("",""),
        line=[1.,0.,0.], origin=[0.,0.,0.], plot="plot", title="", legend="upper right",
        show=False, newfig=True, style="-x"):
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
    for fstr, f in list(functions.items()):
        y = np.array([f(list(t)) for t in z])
        lines.append(getattr(plt, plot)(x, y, style, label=fstr))
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

def assertdir(DIR):
    if not os.path.exists(DIR):
        os.makedirs(DIR)

def savefigs(name="fig", DIR="/tmp/", size=None, bbox_inches="tight",
             ending=".eps", **params):
    if not DIR.endswith("/"): DIR = DIR + "/"
    assertdir(DIR)
    if len(plt.get_fignums()) == 1:
        num = plt.get_fignums()[0]
        fig = plt.figure(num)
        label = fig.get_label()
        label = "" if label=="" else "_"+label
        if size is not None:
            fig.set_size_inches(size)
        fig.savefig(DIR + name + label + ending, bbox_inches=bbox_inches,
                    **params)
        return
    for num in plt.get_fignums():
        fig = plt.figure(num)
        label = fig.get_label()
        label = str(num) if label=="" else label
        if size is not None:
            fig.set_size_inches(size)
        fig.savefig(DIR + name + "_" + label + ending, bbox_inches=bbox_inches,
                    **params)

def saveplots(name="plot", meta=None, uid=False):
    # collect data from every open figure
    plots = []
    figs = list(map(plt.figure, plt.get_fignums()))
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

# TODO: to make this truly magical, we could recursively modify all parents
def add_params(PARENT=None, **params):
    # this is some magic to attach a set of parameters to a module
    # but be able to use them in the same module without any ugly complication
    if PARENT is not None:
        pparams = PARENT.PARAMS
        params.update({key: pparams[key] for key in pparams if not key in params})
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

def user_params(default=None, **params):
    "cleaner version of add_params without secret module mangling"
    if default is not None:
        params.update({key: default[key] for key in default if not key in params})
    args = _argparse()
    params.update({key: args[key] for key in args if key in params})
    #frm = inspect.stack()[1]
    #mod = inspect.getmodule(frm[0])
    #if not hasattr(mod, "PARAMS"):
    #    mod.PARAMS = dict()
    #mod.PARAMS.update(params)
    return Params(params)

def user_param(**params):
    "expect len(params) = 1"
    return user_params(**params)[list(params.keys())[0]]

def any_params(**params):
    args = _argparse()
    params.update(args)
    return Params(params)

class Params(dict):
    "for writing params.Qmol instead of params['Qmol']"
    def __getattr__(self, key):
        return self[key]
    def __or__(self, other):
        new = Params(self)
        new.update(other)
        return new

def union(*seq):
    return reduce(lambda x, y: x | y, seq)

def dict_union(*seq):
    def union(a, b):
        new = dict(a)
        new.update(b)
        return new
    return reduce(union, seq)

def RectangleMesh(a, b, nx, ny):
    return dolfin.RectangleMesh(dolfin.Point(array(a)), dolfin.Point(array(b)), nx, ny)

def convert3D(mesh3D, *forces):
    "convert force from axisymmetric 2D simulation to 3D vector function"
    def rad(x, y):
        return dolfin.sqrt(x**2 + y**2)

    class Convert3DExpression(dolfin.Expression):
        def __init__(self, F, **kwargs):
            self.F = F
        def value_shape(self):
            return (3,)
        def eval(self, value, x):
            r = rad(x[0], x[1])
            F = self.F([r, x[2]])
            if r==0.:
                value[0] = 0.
                value[1] = 0.
                value[2] = F[1]
            else:
                value[0] = x[0]/r*F[0]
                value[1] = x[1]/r*F[0]
                value[2] = F[1]

    #U = dolfin.FunctionSpace(mesh3D, "CG", 1)
    #V = dolfin.MixedFunctionSpace([U, U, U])
    V = dolfin.VectorFunctionSpace(mesh3D, "CG", 1)
    def to3D(F):
        F2 = dolfin.project(Convert3DExpression(F, degree=1), V)
        #F2 = dolfin.Function(V)
        #F2.interpolate(Convert3DExpression(F))
        return F2
    return tuple(map(to3D, forces))

def convert2D(mesh2D, *forces):
    "convert force from axisymmetric 2D simulation to 2D vector function"
    dolfin.parameters['allow_extrapolation'] = False
    class Convert2DExpression(dolfin.Expression):
        def __init__(self, F, **kwargs):
            self.F = F
        def value_shape(self):
            return (2,)
        def eval(self, value, x):
            r = abs(x[0])
            F = self.F([r, x[1]])
            if r==0.:
                value[0] = 0.
                value[1] = F[1]
            else:
                value[0] = x[0]/r*F[0]
                value[1] = F[1]

    #U = dolfin.FunctionSpace(mesh2D, "CG", 1)
    #V = dolfin.MixedFunctionSpace([U, U])
    V = dolfin.VectorFunctionSpace(mesh2D, "CG", 1)
    def to2D(F):
        F2 = dolfin.project(Convert2DExpression(F, degree=1), V)
        #F2 = dolfin.Function(V)
        #F2.interpolate(Convert3DExpression(F))
        return F2
    return tuple(map(to2D, forces))

class Collector(list):
    new = None

def collect(iterator):
    lst = Collector([])
    for i in iterator:
        yield i, lst
        lst.append(lst.new)

class CollectorDict(dict):
    new = None

def collect_dict(iterator):
    result = CollectorDict({})
    for i, obj in enumerate(iterator):
        yield obj, result
        if i==0:
            for key in result.new:
                result[key] = [result.new[key]]
        else:
            for key in result:
                result[key].append(result.new[key])

def print_dict_difference(first, second):
    print("First w/o second:", {k: first[k] for k in first if not k in second or first[k]!=second[k]})
    print("Second w/o first:", {k: second[k] for k in second if not k in first or first[k]!=second[k]})

def printnow(s):
    print(s, end=' ')
    sys.stdout.flush()

tics = []
def tic():
    tics.append(time.time())
def toc():
    return time.time() - tics.pop()

class Log(object):
    def __init__(self, msg):
        self.msg = msg
    def __enter__(self):
        printnow(self.msg)
        tic()
    def __exit__(self, *args):
        #print("")
        print("%.2g s" %(toc(),))
