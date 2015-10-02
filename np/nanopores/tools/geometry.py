from dolfin import *
import nanopores
from nanopores.tools.illposed import AdaptableBC,adaptmeshfunction,adaptfunction
from nanopores.physics import params_physical

from warnings import warn
from importlib import import_module
import types

__all__ = ["Geometry", "PhysicalBC", "geo_from_name", "geo_from_subdomains", "geo_from_xml"]

class Geometry(object):
    """ Interface between numerical routines and files describing the geometry.

        Holds information such as Mesh, MeshFunctions, etc. and
        supplies BCs and Measures corresponding to physical meanings.
        Example use:

        After initialization somewhere in the program,

        >>> geo = Geometry(mesh=mesh,subdomains=subdomains,...)

        possibly depending on an external geometry module,

        >>> import nanopores.geometry
        >>> geo = Geometry(nanopores.geometry)

        a subroutine is passed the geo object and can do stuff like this:

        >>> # Set up Poisson problem with piecewise constant permittivity
        >>> # and mixed boundary conditions
        >>> V = FunctionSpace(geo.mesh,'CG',1)
        >>> dx = geo.dx()
        >>> ds = geo.ds('charged_boundary')
        >>> eps = geo.pwconst('permittivity')
        >>> #phi = geo.pwconst('boundary_charge') # TODO: doesn't work yet
        >>> # Define form
        >>> u = TrialFunction(V)
        >>> v = TestFunction(V)
        >>> F = inner(eps*grad(u),grad(v))*dx - phi*v*ds
        >>> # Define boundary condition
        >>> bc = geo.BC(V, Constant(0.0), 'ground')

        Note that the problem definition here only knows about the *physical*
        meanings of coefficients, boundary conditions, etc., and the translation
        of these into quantities with *geometrical* content is provided by the
        external geometry module. The Geometry class should enable easy
        communication between the two.
        """

    def __init__(self, module=None, mesh=None,
                 #TODO: names are stupid
                 subdomains=None, boundaries=None,
                 physical_domain=None, physical_boundary=None,
                 synonymes=None, params=None):
        if module:
            exec('from %s import *' %module.__name__)
        self.mesh = mesh
        self.subdomains = subdomains
        self.boundaries = boundaries
        self.params = params
        self.physics = params_physical

        self._physical_domain = physical_domain
        self._physical_boundary = physical_boundary

        #print synonymes
        if synonymes:
            for i in range(3):
                self.import_synonymes(synonymes)
        #print self._physical_domain
        #print self._physical_boundary

        self._dom2phys = _invert_dict_nonunique(physical_domain)
        self._bou2phys = _invert_dict_nonunique(physical_boundary)
        self.dg = {}

    def physicalboundary(self, string):
        try:
            return self._physical_boundary[string]
        except (TypeError,KeyError):
            dolfin_error(__name__+".py",
                         "interprete physical boundary description",
                         "This geometry has not implemented '%s'" %string)

    def physicaldomain(self, string):
        try:
            return self._physical_domain[string]
        except (TypeError,KeyError):
            dolfin_error(__name__+".py",
                "interprete physical domain description",
                "This geometry has not implemented '%s'" %string)

    def interprete(self, string):
        if self.subdomains and self._physical_domain and self._physical_domain.has_key(string):
            return (self.subdomains,self._physical_domain[string])
        elif self.boundaries and self._physical_boundary and self._physical_boundary.has_key(string):
            return (self.boundaries,self._physical_boundary[string])
        else:
            dolfin_error(__name__+".py",
                "interprete string description",
                "This geometry has not implemented '%s'" %string)

    def dx(self, string=None):
        if string:
            return Measure("dx", domain=self.mesh, subdomain_data=self.subdomains,
                    subdomain_id=self.physicaldomain(string))
        else:
            return Measure("dx", domain=self.mesh, subdomain_data=self.subdomains)

    def ds(self, string=None):
        if string:
            return Measure("ds", domain=self.mesh, subdomain_data=self.boundaries,
                    subdomain_id=self.physicalboundary(string))
        else:
            return Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    def dS(self, string=None):
        if string:
            return Measure("dS", domain=self.mesh, subdomain_data=self.boundaries,
                    subdomain_id=self.physicalboundary(string))
        else:
            return Measure("dS", domain=self.mesh, subdomain_data=self.boundaries)

    def BC(self, V, g, string=None):
        if string:
            return PhysicalBC(V, g, string, self)
        else: # TODO: implement AdaptableBC(V,g,SubDomain)
            return DirichletBC(V, g, DomainBoundary())

    def NeumannRHS(self, v, string, value=None, axisym=False):
        # L = geo.NeumannRHS(v, "surfcharge") == charge("dna")*v*dS("dna") +
        #                                    charge("mol")*v*dS("mol") + ...
        # thus v should be TestFunction (or Function)
        if value is None:
            try:
                value = vars(self.physics)[string]
            except KeyError:
                dolfin_error(__name__+".py",
                             "interprete string description",
                             "The module %s has not implemented '%s'" % (self.physics.__name__, string))

        bou2value = self._neumann_lookup(self._bou2phys, value)
        dS = self.dS()
        return sum([avg(inner(bou2value[i], v)) * dS(i) for i in bou2value])

    def linearRHS(self, v, string, value=None, axisym = False):
        # L = geo.linearRHS(v, "volcharge") == charge("mol")*v*dx("mol") + ...
        # thus v should be TestFunction (or Function)
        if value is None:
            try:
                value = vars(self.physics)[string]
            except KeyError:
                dolfin_error(__name__+".py",
                             "interprete string description",
                             "The module %s has not implemented '%s'" % (self.physics.__name__, string))

        dom2value = self._pwconst_lookup(self._dom2phys,value)
        dx = self.dx()
        return sum([inner(dom2value[i], v) * dx(i) for i in dom2value])

    def pwconst(self, string, value=None, DG=True): #TODO: should also work as in docstring
        if DG and self.dg.has_key(string):
            return self.dg[string]
        if value is None:
            try:
                value = vars(self.physics)[string]
            except KeyError:
                dolfin_error(__name__+".py",
                             "interprete string description",
                             "The module %s has not implemented '%s'" % (self.physics.__name__, string))

        dom2value = self._pwconst_lookup(self._dom2phys,value)
        #print value
        #print dom2value
        if DG:
            dgfun = self._dict_to_DG(dom2value)
            self.dg[string] = dgfun
            return dgfun
        else:
            return value

    def parameter(self, string):
        try:
            return self.params[string]
        except (TypeError,KeyError):
            dolfin_error(__name__+".py",
                "find parameter",
                "This geometry has not implemented '%s'" %string)

    def adapt(self,mesh):
        self.mesh = mesh
        #print "subdomain id (geo):",self.subdomains.id()
        self.subdomains = adaptmeshfunction(self.subdomains,mesh)
        self.boundaries = adaptmeshfunction(self.boundaries,mesh)
        #plot(self.boundaries, interactive=True)
        #print "subdomain id (geo):",self.subdomains.id()
        for f in self.dg.values():
            adaptfunction(f, mesh, interpolate=True, assign=True)

    # alternative to adapt, should be overwritten dynamically
    rebuild = adapt

    def import_synonymes(self,syns):
        # philosophy behind synonymes: the more, the better!
        for syn in syns:
            if not isinstance(syns[syn],set):
                syns[syn] = {syns[syn]}

        for syn in syns:
            for dic in [self._physical_domain,
                        self._physical_boundary]:
                t = set()
                for phys in syns[syn]:
                    if dic.has_key(phys):
                        t = t | set(dic[phys])
                    else:
                        t = None
                        break
                if t is not None:
                    dic[syn] = tuple(t)

    def submesh(self, string):
        # assumes volume submesh
        t = self.physicaldomain(string)
        if len(t)==0:
            return SubMesh(self.mesh, self.subdomains, 1234)
        elif len(t)==1:
            return SubMesh(self.mesh, self.subdomains, t[0])
        else:
            return SubMesh(self.mesh, self.indicator(string), 1)

    def indicator(self, string):
        # return "indicator" CellFunction for subdomain
        t = self.physicaldomain(string)
        sub = self.subdomains
        chi = CellFunction("size_t", self.mesh, 0)
        for cell in cells(self.mesh):
            if sub[cell] in t:
                chi[cell] = 1
        return chi

    def _neumann_lookup(self, bou2phys, value):
        bou2value = {}
        for i in bou2phys:
            for s in bou2phys[i]:
                if s in value:
                    if i in bou2value and not bou2value[i] == value[s]:
                        dolfin_error(__name__+".py",
                            "create Neumann RHS",
                            "The value on '%s' is ambigous, check %s"%(s,self.physics.__name__))
                    else:
                        bou2value[i] = value[s]
        # partially defined is ok => don't need "if not i in bou2value:..."
        return bou2value

    def _pwconst_lookup(self, dom2phys, value):
        dom2value = {}
        for i in dom2phys:
            for s in dom2phys[i]:
                if value.has_key(s):
                    if dom2value.has_key(i) and (not value[s] == dom2value[i]):
                        dolfin_error(__name__+".py",
                            "create piecewise constant",
                            "The value on '%s' is ambigous, check %s" %(s,self.physics.__name__))
                    else:
                        dom2value[i] = value[s]
            if not dom2value.has_key(i):
                warn("While creating piecewise constant, there was no value assigned to the %s domain" %str(dom2phys[i]))
        return dom2value

    def _dict_to_DG(self, dom2value): #TODO: not assume domain
        expr = Dict2Expression(dom2value, self.subdomains)
        dgfun = Function(FunctionSpace(self.mesh,'DG',0))
        dgfun.interpolate(expr)
        return dgfun

class Dict2Expression(Expression): #TODO: too slow... --> compiled expr??
    def __init__(self, dom2value, subdomains):
        self.subdomains = subdomains
        self.dom2value = dom2value
    def eval_cell(self, values, x, cell):
        values[0] = self.dom2value[self.subdomains[cell.index]]


class PhysicalBC(object):
    """ boundary condition defined by its physical meaning
        together with a geometry object to interprete that.
        can also be in an 'unrealized' state without geometry object """
    def __init__(self, V, g, description, geo=None):
        self.real = False
        self.V = V
        self.g = g
        self.description = description
        if geo:
            self.realize(geo)

    def realize(self, geo):
        boundaries = geo.boundaries
        numbers = geo.physicalboundary(self.description)
        self.bcs = [AdaptableBC(self.V, self.g, boundaries, i) for i in numbers]
        self.boundaries = boundaries
        self.real = True

    def apply(self, *args, **kwargs):
        for bc in self.bcs:
            bc.apply(*args,**kwargs)

    def adapt(self, *args, **kwargs):
        self.bcs = [bc.adapt(*args, **kwargs) for bc in self.bcs]
        return self

    def value(self):
        return self.bcs[0].value() if self.bcs else self.g

    def function_space(self):
        return self.bcs[0].function_space() if self.bcs else self.V

    def homogenize(self):
        for bc in self.bcs:
            bc.homogenize()


def _invert_dict(d):
    if isinstance(d,dict):
        return {i:s for s,t in d.items() for i in t}

# TODO:
def _invert_dict_nonunique(dom):
    d = {}
    for s in dom:
        for i in dom[s]:
            if d.has_key(i):
                d[i].append(s)
            else:
                d[i] = [s]
    return d


# get Geometry
def make_domain(mesh, list, check_midpoint=False):
    subdomains = CellFunction("size_t", mesh, 0)
    physical_domain = {}
    for i,sub in enumerate(list):
        sub.mark(subdomains, i, check_midpoint)
        name = type(sub).__name__.lower()
        physical_domain[name] = (i,)
    return (subdomains, physical_domain)

def make_boundary(mesh, list, check_midpoint=False):
    i0 = 0
    subdomains = FacetFunction("size_t", mesh, i0)
    physical_domain = {}
    for i,sub in enumerate(list):
        if hasattr(sub, "check_midpoint"):
            sub.mark(subdomains, i+i0+1, sub.check_midpoint)
        else:
            sub.mark(subdomains, i+i0+1, check_midpoint)
        name = type(sub).__name__.lower()
        physical_domain[name] = (i+i0+1,)
    return (subdomains, physical_domain)

# DEPRECATED: this is unnessecary, since the same could be accomplished by
#             providing a SubDomain with its mark() function overwritten.
def mark_domains_with_function(subdomains, physical_domain, flist, params):
    j = len(physical_domain.values())
    for i,f in enumerate(flist):
        f(subdomains, i+j, **params)
        name = f.__name__.lower()
        physical_domain[name] = (i+j,)

def geo_from_subdomains(mesh, module, check_midpoint=False, **params):
    subd = import_module(module)

    (subdomains, physical_domain) = make_domain(mesh, subd.subdomain_list(**params), check_midpoint)
    (boundaries, physical_boundary) = make_boundary(mesh, subd.boundaries_list(**params), check_midpoint)
    synonymes = subd.__dict__.get("synonymes")
    if hasattr(subd, "mark_directly"):
        mark_domains_with_function(subdomains, physical_domain, subd.mark_directly(), params)

    return Geometry(None, mesh, subdomains, boundaries, physical_domain, physical_boundary, synonymes, params)

def geo_from_name(name, mesh=None, check_midpoint=False, **params):

    if not mesh:
        mesh = Mesh("%s/%s/mesh/mesh.xml" %(nanopores.DATADIR,name))

    module = "nanopores.geometries.%s.subdomains" %name
    tmp = vars(import_module('nanopores.geometries.%s.params_geo' %name))
    params.update({key:tmp[key] for key in tmp if (not key in params and not key.startswith("__"))})
    params["name"] = name

    def rebuild(self, othermesh):
        other = geo_from_subdomains(othermesh, module, check_midpoint=check_midpoint, **params)
        self.__dict__.update(other.__dict__)
        return other

    geo = geo_from_subdomains(mesh, module, check_midpoint=check_midpoint, **params)
    #setattr(Geometry, 'rebuild', rebuild)        #<-- modifies all instances
    #geo.rebuild = rebuild.__get__(geo, Geometry) #<-- equivalent
    geo.rebuild = types.MethodType(rebuild, geo)
    return geo
    
    
def geo_from_xml(name, mesh=mesh):
    DIR = "%s/%s/mesh/" %(nanopores.DATADIR,name)
    mesh = Mesh(DIR+"mesh.xml")
    subdomains = MeshFunction("size_t", mesh, DIR+"mesh_physical_region.xml")
    boundaries = MeshFunction("size_t", mesh, DIR+"mesh_facet_region.xml")
    
    # quick hack for testing aHem geo
    #physical_domain = {"fluid":(4,), "membrane":(5,), "ahem":(6,)}
    #physical_boundary = {"ext":(1,), "ahem":(2,), "membrane":(3,)}
    
    with open(DIR+"meta.txt", "r") as f:
        meta = eval(f.read())
        
    print meta
    
    if not isinstance(meta, dict):
        meta = {}
    
    return Geometry(None, mesh, subdomains, boundaries, **meta)
    
