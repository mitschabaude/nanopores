# (c) 2016 Gregor Mitscha-Baude
"handle complex PDE solvers and parallel force evaluation"
from nanopores import iterate_in_parallel
from .utilities import Params
from . import fields
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
        
#class cache_force_field(fields.CacheBase):

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
            result = calculate(x0, **params)
            result = {k: [v] for k, v in result.items()}
            fields.save_fields(name, save_params, x=[x0], **result)
        except: # Exception, RuntimeError:
            print "Error occured, continuing without saving."
            Xfailed.append(x0)
            result = None
        return result
    
    results, _ = iterate_in_parallel(run, nproc, **iter_params)
    
    print "failed:"       
    print Xfailed
    print "%d of %d force calculations failed." % (len(Xfailed), len(X))
    fields.update()
    return results