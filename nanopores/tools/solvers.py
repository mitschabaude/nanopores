# (c) 2016 Gregor Mitscha-Baude
"handle complex PDE solvers and parallel force evaluation"
import traceback
from nanopores.scripts.simulation2D import iterate_in_parallel
from nanopores.tools.utilities import Params
from nanopores.tools import fields
__all__ = ["Setup", "calculate_forcefield"]

class Setup(object):
    "handle input parameters and setup geometry"
    default = {}
    
    def __init__(self, geop=None, physp=None, solverp=None, **params):
        self.init_params(params, geop=geop, physp=physp, solverp=solverp)
        self.init_geo()
        self.init_phys()
        
    def init_params(self, params, **paramsets):
        for p in paramsets:
            setattr(self, p, Params(self.default[p]))
            dic = getattr(self, p)
            if paramsets[p] is not None:
                dic.update(paramsets[p])
            for k in dic:
                if k in params:
                    dic[k] = params[k]
                    
    # subclasses have to overwrite
    def init_geo(self): 
        self.geo = None
    def init_phys(self):
        self.phys = None

def calculate_forcefield(name, X, calculate, params={}, default={}, nproc=1):
    "assuming function calculate(x0, **params)"
    #fields.update()
    if "x0" in default: default.pop("x0")
    save_params = dict(default, **params)
    # TODO calculate in parallel
    N = len(X)
    if fields.exists(name, **save_params):
        Xdone = fields.get_field(name, "x", **save_params)
        X = [x0 for x0 in X if x0 not in Xdone]
        print "Existing force file found, %d/%d points remaining." % (
            len(X), N)
    Xfailed = []
    iter_params = dict(x0=X)
    
    def run(x0=None):
        try:
            result = calculate([x0], **params)
            #result = {k: [v] for k, v in result.items()}
            fields.save_fields(name, save_params, x=[x0], **result)
        except: # Exception, RuntimeError:
            print "Error occured, continuing without saving."
            print traceback.print_exc()
            Xfailed.append(x0)
            result = None
        return result
    
    results, _ = iterate_in_parallel(run, nproc, **iter_params)
    
    if nproc == 1:
        print "%d of %d force calculations failed." % (len(Xfailed), len(X))
    fields.update()
    return results
    
class cache_forcefield(fields.CacheBase):
    "caching decorator for function calculate(X, **params) --> dict()"
    def __init__(self, name, default={}, nproc=1):
        self.name = name
        self.default = default
        self.nproc = nproc
        
    def __call__(self, f):
        def wrapper(X, nproc=self.nproc, name=self.name, **params):
            # calculate remaining points (in parallel)
            calculate_forcefield(name, X, f, params,
                                 self.default, nproc)
            # load requested data points
            try:
                result = fields.get_fields(name, **params)
                I = [i for i, x in enumerate(result["x"]) if x in X]
            except KeyError:
                result = {}
                I = []
            result = {key: [val[i] for i in I] for key, val in result.items()}
            return result
        return wrapper

    

    