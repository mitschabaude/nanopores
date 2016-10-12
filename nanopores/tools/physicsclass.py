''' Physics class '''

#from .utilities import import_vars
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

    def __init__(self, name="default", geo=None, module=None, **params):
        # initialize from module nanopores.physics.name
        if module is not None:
            mod = module
        else:
            name = "nanopores.physics."+name
            mod = import_module(name)
        mod = reload(mod)
        #var = import_vars(name)
        vardic = vars(mod)
        var = {k: vardic[k] for k in vardic if not k.startswith("_")}
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
                setattr(mod, k, v) # TODO ????

        #self.precalculate(mod) # could be optional at some point, for performance
        self.mod = mod
        self.__name__ = mod.__name__

        # crazy way to finish
        #self.__dict__ = self.base

        # change physical parameters in geo to self
        # TODO: not ideal.. should we move all "physics functionality" from Geometry to Physics?
        if geo:
            self.base["geo"] = geo
            geo.physics = self
            # conservative import of default synonymes
            # TODO: i'm not totally sure about this...
            # should some synonymes be a *physical* property rather than *geometrical*?
            if "synonymes" in self.maps:
                syns = self.maps.pop("synonymes")
                self.base["synonymes"] = syns
                geo.import_synonymes(syns, conservative=True)
            
            self.toadapt = self.base.pop("toadapt") if "toadapt" in self.base else []

    def precalculate(self, mod):
        for fstr, f in self.functions.items():
            argnames = inspect.getargspec(f).args
            args = [self.base[k] for k in argnames]
            self.base[fstr] = f(*args)
            #setattr(mod, fstr, self.base[fstr])

        for mstr, m in self.maps.items():
            for k in m:
                if isinstance(m[k], str):
                    m[k] = self.base[m[k]]
            self.base[mstr] = m
            #setattr(mod, mstr, m)
            
    def __getattr__(self, name):
        if name in self.base:
            return self.base[name]
            
        elif name in self.functions:
            f = self.functions[name] #.pop(name)
            argnames = inspect.getargspec(f).args
            args = [getattr(self, k) for k in argnames]
            result = f(*args)
            self.base[name] = result
            #setattr(self.mod, name, result)
            return result
            
        elif name in self.maps:
            m0 = self.maps[name] #.pop(name)
            self.base[name] = m = dict(m0)
            for k in m:
                if isinstance(m[k], str):
                    m[k] = getattr(self, m[k])
            self.base[name] = m
            #setattr(self.mod, name, m)
            return m
            
    # FIXME this next two functions were never tried out because we did not actually need them
    # could be useful at some time though
    """
    def _getfunction(self, name):
        f = self.functions[name]
        argnames = inspect.getargspec(f).args
        args = [getattr(self, k) for k in argnames]
        result = f(*args)
        self.base[name] = result
        setattr(self.mod, name, result)
        return result
        
    def _getmap(self, name):
        m0 = self.maps[name]
        self.base[name] = m = dict(m0)
        for k in m:
            if isinstance(m[k], str):
                m[k] = getattr(self, m[k])
        self.base[name] = m
        setattr(self.mod, name, m)
        return m
          
    def adapt(self):
        # simply recalculate stuff in self.toadapt
        print "Adapting Physics."
        for name in self.toadapt:
            if name in self.functions:
                res = self._getfunction(name)
            elif name in self.maps:
                res = self._getmap(name)
            print "   %s: %s" %(name, res)
    """    
            
    def __str__(self):
        output = ["","Physics:"]
        s = "   %s = %s"
        for dstr in ["base", "functions", "maps"]:
            d = getattr(self, dstr)
            output.append("%s (%d):" % (dstr, len(d)))
            for item in d.items():
                output.append(s % item)
            output.append("")
        return "\n".join(output)
        
