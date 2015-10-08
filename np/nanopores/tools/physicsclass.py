''' Physics class '''

from .utilities import import_vars
from collections import defaultdict
from importlib import import_module
import inspect

__all__ = ["Physics"]

d = defaultdict(lambda :"other", {
type(lambda: None) : "function",
dict : "dict",
str : "str",
float : "float",
int : "int",
type(inspect): "module",
})
def typestr(t): return d[type(t)]

class Physics(object):

    def __init__(self, name="default", geo=None, **params):
        # initialize from module nanopores.physics.name
        name = "nanopores.physics."+name
        mod = import_module(name)
        mod = reload(mod)
        var = import_vars(name)

        #print var

        # override with user-specified parameters
        var.update(params)

        # sort according to method of calling
        self.base = defaultdict(lambda : None)
        self.functions = {}
        self.maps = {}
        for k,v in var.items():
            if typestr(v) is "function":
                self.functions[k] = v
            elif typestr(v) is "dict":
                self.maps[k] = v
            elif typestr(v) is "module":
                pass
            else:
                self.base[k] = v
                setattr(mod, k, v)

        # calculate all dependent parameters
        if geo:
            self.base["geo"] = geo
        self.precalculate(mod) # could be optional at some point, for performance

        # crazy way to finish
        self.__dict__ = self.base

        # change physical parameters in geo to self
        # TODO: not ideal.. should we move all "physics functionality" from Geometry to Physics?
        if geo:
            geo.physics = self

    def precalculate(self, mod):
        for fstr, f in self.functions.items():
            argnames = inspect.getargspec(f).args
            args = [self.base[k] for k in argnames]
            self.base[fstr] = f(*args)
            setattr(mod, fstr, self.base[fstr])

        for mstr, m in self.maps.items():
            for k in m:
                if isinstance(m[k], str):
                    m[k] = self.base[m[k]]
            self.base[mstr] = m
            setattr(mod, mstr, m)
