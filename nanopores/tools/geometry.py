from dolfin import *
import ufl
import nanopores
from nanopores.tools.illposed import AdaptableBC, adaptmeshfunction, adaptfunction
#from nanopores.tools.physicsclass import Physics
from nanopores.physics import params_physical
import dolfin
import numpy as np

from importlib import import_module
import types

__all__ = ["Geometry", "PhysicalBC", "geo_from_name", "geo_from_subdomains",
           "geo_from_xml", "geo_from_xml_threadsafe", "PointBC"]

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
            exec 'from %s import *' %module.__name__ in globals(), locals()

        self.mesh = mesh

        # trivial default subdomains/boundaries
        if subdomains is None:
            subdomains = CellFunction("size_t", mesh, 0)
            physical_domain = {"domain":(0,)}
        if boundaries is None:
            boundaries = FacetFunction("size_t", mesh, 0)
            AutoSubDomain(lambda x, on_boundary : on_boundary).mark(boundaries, 1)
            physical_boundary = {"boundary":(1,)}

        self.subdomains = subdomains
        self.boundaries = boundaries
        self.params = params if params else {}
        self.synonymes = {}
        self.physics = params_physical

        self._physical_domain = physical_domain
        self._physical_boundary = physical_boundary

        self.import_synonymes((synonymes if synonymes else {}))

        self.dg = {}
        self.constants = {}
        self.volumes = {}

    def physicalboundary(self, string):
        try:
            return self._physical_boundary[string]
        except (TypeError,KeyError):
            dolfin_error(__name__+".py",
                         "interprete physical boundary description",
                         "This geometry has not implemented '%s'" %string)
        return tup

    def physicaldomain(self, string):
        try:
            return self._physical_domain[string]
        except (TypeError,KeyError):
            dolfin_error(__name__+".py",
                "interprete physical domain description",
                "This geometry has not implemented '%s'" %string)

    def _notempty(self, tup):
        "assure that assemble evaluates to 0. if subdomain is empty"
        return tup if len(tup)>0 else (9999,)

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
                    subdomain_id=self._notempty(self.physicaldomain(string)))
        else:
            return Measure("dx", domain=self.mesh, subdomain_data=self.subdomains)

    def ds(self, string=None):
        if string:
            return Measure("ds", domain=self.mesh, subdomain_data=self.boundaries,
                    subdomain_id=self._notempty(self.physicalboundary(string)))
        else:
            return Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    def dS(self, string=None):
        if string:
            return Measure("dS", domain=self.mesh, subdomain_data=self.boundaries,
                    subdomain_id=self._notempty(self.physicalboundary(string)))
        else:
            return Measure("dS", domain=self.mesh, subdomain_data=self.boundaries)

    def BC(self, V, g, string=None):
        if string:
            return PhysicalBC(V, g, string, self)
        else: # TODO: implement AdaptableBC(V,g,SubDomain)
            return DirichletBC(V, g, DomainBoundary())

    # TODO: pwconstBC is DEPRECATED since pwBC has strictly more general functionality
    def pwconstBC(self, V, string, homogenize=False, value=None):
        " piecewise constant boundary condition from dict(boundary = value) where value is a number (or None) "
        value = self._getvalue(string, value)
        #print "DEBUG:", value
        if isinstance(value, dict):
            if homogenize:
                value = {key: 0. for key in value}
            return [self.BC(V, Constant(value[key]), key) for key in value if (value[key] is not None)]
        else: # assume value is number and apply on whole boundary
            if homogenize:
                value = 0.
            return [self.BC(V, Constant(value))] if value is not None else []

    # TODO: OMG this is ugly code. readability = -1000000
    def pwBC(self, V, string, homogenize=False, value=None):
        """ piecewise boundary condition from dict(boundary = value) where value is a dolfin.GenericFunction.
            will wrap with dolfin.Constant is input is a number/tuple, and do nothin if value is None. """
        value = self._getvalue(string, value)
        #print "DEBUG:", value
        if isinstance(value, dict):
            # TODO: get rid of homogenize (it is only used for PB and is ugly, because assumes scalar pde)
            if homogenize:
                value = {key: Constant(0.) for key in value}
            bcs = []
            for key, val in value.items():
                if val is not None:
                    if isinstance(val, GenericFunction):
                        #print type(val)
                        #print val.__class__
                        bcs.append(self.BC(V, val, key))
                    elif isinstance(val, float) or isinstance(val, tuple):
                        bcs.append(self.BC(V, Constant(val), key))
                    else:
                        dolfin_error(__name__+".py",
                            "assign boundary condition",
                            "Value on '%s' for the BC '%s' is of unexpected type %s" % (key, string, type(val)))
            return bcs
        else: # try to assign value on whole boundary
            if homogenize:
                value = Constant(0.)
            if value is None:
                return []
            elif isinstance(value, GenericFunction):
                bc = self.BC(V, value)
            elif isinstance(val, float) or isinstance(val, tuple):
                bc = self.BC(V, Constant(value))
            else:
                dolfin_error(__name__+".py",
                    "assign boundary condition",
                    "Value for the BC '%s' is of unexpected type '%s'." % (string, type(val)))
            return [bc]

    def VolumeBC(self, V, name, f):
        return _VolumeBC(V, self, name, f)

    def PointBC(self, V, points, values):
        return PointBC(V, points, values)

    def _getvalue(self, string, value):
        if value is None:
            try:
                value = getattr(self.physics, string)
            except KeyError:
                dolfin_error(__name__+".py",
                    "interprete string description",
                    "The module %s has not implemented '%s'" % (self.physics.__name__, string))
        return value

    def NeumannRHS(self, v, string=None, value=None):
        # L = geo.NeumannRHS(v, "surfcharge") == charge("dna")*v*dS("dna") +
        #                                    charge("mol")*v*dS("mol") + ...
        # thus v should be TestFunction (or Function)
        value = self._getvalue(string, value)
        bou2value = self._neumann_lookup(self._bou2phys, value)
        #print bou2value
        dS = self.dS()
        ds = self.ds()
        return sum([avg(inner(_wrapf(bou2value[i]), v)) * dS(i) for i in bou2value]) \
             + sum([inner(_wrapf(bou2value[i]), v) * ds(i) for i in bou2value])

    def linearRHS(self, v, string=None, value=None):
        # L = geo.linearRHS(v, "volcharge") == charge("mol")*v*dx("mol") + ...
        # thus v should be TestFunction (or Function)
        # value can be dict or float
        value = self._getvalue(string, value)
        dx = self.dx()

        if isinstance(value, dict):
            dom2value = self._neumann_lookup(self._dom2phys, value)
            return sum([inner(_wrapf(dom2value[i]), v) * dx(i) for i in dom2value])
        else:
            return inner(_wrapf(value), v) * dx

    def pwconst(self, string, value=None, DG=True): #TODO: should also work as in docstring
        if DG and self.dg.has_key(string):
            return self.dg[string]
        value = self._getvalue(string, value)
        dom2value = self._pwconst_lookup(self._dom2phys, value)
        #print value
        #print self._dom2phys
        #print dom2value
        if DG:
            dgfun = self._dict_to_DG(dom2value)
            self.dg[string] = dgfun
            return dgfun
        else:
            return value

    def constant(self, name, compute=None):
        # compute has to be provided the first time the constant is requested
        if not name in self.constants:
            self.constants[name] = GeometricConstant(name, compute, self)
            print "Computed %s." %(self.constants[name],)
        return self.constants[name].function


    def volume(self, string, dx="dx", cyl=True):
        if dx in self.volumes and string in self.volumes[dx]:
            return self.volumes[dx][string]
        dmu = getattr(self, dx)(string)
        r = Constant(1.0) #if not cyl else Expression("2*pi*x[0]")
        vol = assemble(r*dmu)
        cvol = Constant(vol)
        if dx in self.volumes:
            self.volumes[dx][string] = cvol
        else:
            self.volumes[dx] = {string: cvol}
        return cvol

    def avg(self, u, string, dx="dx"):
        meas = getattr(self, dx)(string)
        return assemble(u*meas)/assemble(Constant(1.0)*meas)

    def parameter(self, string):
        try:
            return self.params[string]
        except (TypeError,KeyError):
            dolfin_error(__name__+".py",
                "find parameter",
                "This geometry has not implemented '%s'" %string)

    def adapt(self,mesh):
        # save history of past meshes and meshfunctions
        # --> to prevent segfaults because the old stuff gets garbage-collected!
        if not hasattr(self, "old"):
            self.old = []
        self.old.append((self.mesh, self.subdomains, self.boundaries))

        self.mesh = mesh
        #print "subdomain id (geo):",self.subdomains.id()
        self.subdomains = adaptmeshfunction(self.subdomains, mesh)
        self.boundaries = adaptmeshfunction(self.boundaries, mesh)
        # adapt functions and constants
        for f in self.dg.values():
            adaptfunction(f, mesh, interpolate=True, assign=True)

        # if curved boundaries are defined, snap back those
        if hasattr(self, "curved"):
            for boundary, snap in self.curved.items():
                #print "Adapting curved boundary '%s'." % boundary
                self.snap_to_boundary(boundary, snap)

        for const in self.constants.values():
            const.recompute()
            #print "Recomputed %s." %const

        for meas in self.volumes:
            for name in self.volumes[meas]:
                dmu = getattr(self, meas)(name)
                vol = assemble(Constant(1.0)*dmu)
                print "DEBUG New volume:", vol
                self.volumes[meas][name].assign(vol)
        #print self.volumes

        # TODO maybe needed some time
        # adapt self.Physics if we have one
        #if isinstance(self.physics, Physics):
        #    self.physics.adapt()


    # alternative to adapt, should be overwritten dynamically
    rebuild = adapt

    def import_synonymes(self, synonymes, conservative=False):
        if conservative:
            synonymes.update(self.synonymes)
        self.synonymes.update(synonymes)
        for i in range(3):
            self._import_synonymes(self.synonymes)
        #print self._physical_domain
        #print self._physical_boundary
        self._dom2phys = _invert_dict_nonunique(self._physical_domain)
        self._bou2phys = _invert_dict_nonunique(self._physical_boundary)

    def _import_synonymes(self,syns):
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

    def add_subdomain(self, string, marker):
        i = max(self._dom2phys.keys()) + 1
        marker.mark(self.subdomains, i)
        self._physical_domain[string] = (i,)
        self._dom2phys[i] = [string]

    def submesh(self, string):
        # assumes volume submesh
        t = self.physicaldomain(string)
        if len(t)==0:
            return SubMesh(self.mesh, self.subdomains, 1234)
        elif len(t)==1:
            return SubMesh(self.mesh, self.subdomains, t[0])
        else:
            return SubMesh(self.mesh, self.indicator(string), 1)

    def indicator(self, string, DG=False, callable=False):
        # return "indicator" CellFunction for subdomain
        # set either DG or callable to True to use as function (callable with points)
        t = self.physicaldomain(string)
        sub = self.subdomains

        if DG: # <-- slow. callable is better
            D = FunctionSpace(self.mesh, "DG", 0)
            dofmap = D.dofmap()
            chi = Function(D)
            chi.vector()[:] = 0
            for cell in cells(self.mesh):
                if sub[cell] in t:
                    i = dofmap.cell_dofs(cell.index())
                    chi.vector()[i] = 1
            return chi

        chi = CellFunction("size_t", self.mesh, 0)
        for cell in cells(self.mesh):
            if sub[cell] in t:
                chi[cell] = 1

        if callable:
            return CallableMeshFunction(chi)
        return chi

    def snap_to_boundary(self, name, snap, smooth=False):
        mesh = self.mesh
        # get vertices that lie on the boundary (via BC for CG1 function)
        V = FunctionSpace(mesh, 'CG', 1)
        bc = self.BC(V, 1., name)
        u = Function(V)
        bc.apply(u.vector())
        d2v = dof_to_vertex_map(V)
        vertices_on_boundary = d2v[u.vector() == 1.0]

        # DEBUG plot for check
        '''
        testf = VertexFunction("bool", mesh, False)
        testf.array()[vertices_on_boundary] = True
        R = self.params["rMolecule"]
        C = self.params["x0"][::2]
        '''

        # snap those vertices
        for v in vertices_on_boundary:
            x = mesh.coordinates()[v]
            #r0 = sqrt((x[0]-C[0])**2 + (x[1]-C[1])**2)
            snap(x)
            #r1 = sqrt((x[0]-C[0])**2 + (x[1]-C[1])**2)
            #print "DEBUG Radius: %s --> %s" %(r0, r1)
            mesh.geometry().set(v, x)
        if smooth:
            mesh.smooth(1)
        #plot(testf, title="after snap")
        #interactive()

    def _neumann_lookup(self, bou2phys, value):
        bou2value = {}
        for i in bou2phys:
            for s in bou2phys[i]:
                if s in value:
                    if i in bou2value and not bou2value[i] == value[s]:
                        dolfin_error(__name__+".py",
                            "create Neumann or volume RHS",
                            "The value on '%s' is ambigous, check %s" %(s, self.physics.__name__))
                    else:
                        bou2value[i] = value[s]
        # partially defined is ok => don't need "if not i in bou2value:..."
        return bou2value

    def _pwconst_lookup(self, dom2phys, value):
        default = value["default"] if "default" in value else None
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
            if not i in dom2value:
                if default is not None:
                    dom2value[i] = default
                else:
                    #warning("pwconst: no value specified for '%s' (used 0 as default)." %str(dom2phys[i][0]))
                    raise Exception("pwconst: no value specified on (part of) '%s'." %str(dom2phys[i][0]))
                    dom2value[i] = 0.
        return dom2value

    def _dict_to_DG(self, dom2value): #TODO: not assume domain
        expr = Dict2Expression(dom2value, self.subdomains, degree=1)
        dgfun = Function(FunctionSpace(self.mesh,'DG',0))
        dgfun.interpolate(expr)
        return dgfun

    def __str__(self):
        return "Boundaries:\n%s\nSubdomains:\n%s\n" % (self._physical_boundary, self._physical_domain)


class Dict2Expression(Expression): #TODO: too slow... --> compiled expr??
    def __init__(self, dom2value, subdomains, **kwargs):
        self.subdomains = subdomains
        self.dom2value = dom2value
    def eval_cell(self, values, x, cell):
        values[0] = self.dom2value[self.subdomains[cell.index]]


from numpy import array
class CallableMeshFunction(object):
    def __init__(self, f):
        self.f = f
        self.btree = f.mesh().bounding_box_tree()
    def __call__(self, x):
        i = self.btree.compute_first_entity_collision(Point(array(x)))
        return self.f[int(i)]


class GeometricConstant(object):
    def __init__(self, name, compute, geo):
        c = compute(geo)
        self.name = name
        self.function = Constant(c)
        self.value = c
        self.compute = compute
        self.geo = geo
    def recompute(self):
        c = self.compute(self.geo)
        self.function.assign(c)
        self.value = c
    def __str__(self):
        return "%s = %s" %(self.name, self.value)


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
            self.geo = geo
            self.realize(geo)

    def realize(self, geo):
        boundaries = geo.boundaries
        #if True:
        try:
            numbers = geo.physicalboundary(self.description)
            self.bcs = [AdaptableBC(self.V, self.g, boundaries, i) for i in numbers]
        #else:
        except RuntimeError:
            warning("PhysicalBC for '%s' could not be assigned; maybe the boundary was not found." % (self.description,))
            self.bcs = []
        self.boundaries = boundaries
        self.real = True

    def apply(self, *args, **kwargs):
        for bc in self.bcs:
            bc.apply(*args,**kwargs)

    def adapt(self, *args, **kwargs):
        self.bcs = [bc.adapt(*args, **kwargs) for bc in self.bcs]
        return self

    def value(self):
        return self.g

    def function_space(self):
        return self.bcs[0].function_space() if self.bcs else self.V

    def homogenized(self): # TODO: arbitrary tensors
        if dolfin.__version__ == "1.6.0":
            shape = self.g.shape()
        else:
            shape = self.g.ufl_shape
        c = Constant(tuple(0. for i in range(shape[0])) if shape else 0.)
        return PhysicalBC(self.V, c, self.description, self.geo)

    def damp(self, scalar):
        if scalar == 0.:
            raise Exception("PhysicalBC: Cannot damp with 0.")
        new = scalar/self.damping if hasattr(self, "damping") else scalar
        self.damping = scalar
        self.bcs = [bc.damped(new) for bc in self.bcs]

class PointBC(object):

    def __init__(self, V, points, values, tol=1e-5):
        if V.component().shape[0] > 0: # V is subspace
            V = V.collapse()
        self.V = V
        if callable(values):
            self.values = [values(p) for p in points]
        else:
            self.values = values
        self.points = points
        self.bc_f = dolfin.Function(V)
        mesh = V.mesh()

        co = mesh.coordinates()
        dim = co.shape[1]
        dof_map = dolfin.vertex_to_dof_map(V)
        node_set = set()
        bc_values = self.bc_f.vector().array()
        for p, v in zip(points, self.values):
            wh = np.where(sum((co[:,j] - p[j])**2 for j in range(dim)) < tol)[0]
            if wh.shape[0]:
                i = wh[0]
                bc_values[dof_map[i]] = v
                node_set.add(i)
                #print "found:", p
            else:
                pass
                #print "not found:", p
        print "Found %d of %d points." %(len(node_set), len(points))

        self.bc_f.vector().set_local(bc_values)
        self.bc_f.vector().apply("insert") # TODO: what does this do?
        self.dof_set = np.array(dof_map[list(node_set)], dtype="intc")

    def apply(self, a):
        # Manual application of bcs
        if isinstance(a, dolfin.Matrix):
            A = a
            # Modif A: zero bc row & set diagonal to 1
            A.ident_local(self.dof_set)
            A.apply("insert")

        elif isinstance(a, dolfin.GenericVector):
            b = a
            # Modif b: entry in the bc row is taken from bc_f
            bc_values = self.bc_f.vector().array()
            b_values = b.array()
            b_values[self.dof_set] = bc_values[self.dof_set]
            b.set_local(b_values)
            b.apply("insert")

        else:
            dolfin.warning("Could not apply Point BC.")

class _VolumeBC(PointBC):

    def __init__(self, V, geo, name, f):
        self.V = V
        # get dofs lying in subdomain
        dofmap = V.dofmap()
        tup = geo.physicaldomain(name)
        sub = geo.subdomains
        mesh = geo.mesh

        subdofs = set()
        for i, cell in enumerate(dolfin.cells(mesh)):
            if sub[cell] in tup:
                celldofs = dofmap.cell_dofs(i)
                subdofs.update(celldofs)

        subdofs = np.array(list(subdofs), dtype="intc")
        d2v = dolfin.dof_to_vertex_map(V)
        co = mesh.coordinates()

        # create function with desired values
        # could also be implemented with Expression.eval_cell like pwconst
        bc_f = dolfin.Function(V)
        for dof in subdofs:
            x = co[d2v[dof]]
            bc_f.vector()[dof] = f(x)
        self.bc_f = bc_f
        self.dof_set = subdofs


def _wrapf(f):
# takes either Function/uflExpression or float/tuple and wraps with Constant in the latter case
# for easy specification of Dirichlet or Neumann data
    if isinstance(f, GenericFunction) or isinstance(f, ufl.core.expr.Expr):
        return f
    elif isinstance(f, float) or isinstance(f, int) or isinstance(f, tuple):
        return Constant(f)
    else:
        dolfin_error(__name__+".py", "use given function",
            "Function is of unexpected type '%s'." % (type(f),))

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


def geo_from_xml(name):
    DIR = "%s/%s/mesh/" %(nanopores.DATADIR, name)
    mesh = Mesh(DIR+"mesh.xml")
    subdomains = MeshFunction("size_t", mesh, DIR+"mesh_physical_region.xml")
    boundaries = MeshFunction("size_t", mesh, DIR+"mesh_facet_region.xml")

    with open(DIR+"meta.txt", "r") as f:
        meta = eval(f.read())

    physdom = meta.pop("physical_domain")
    physbou = meta.pop("physical_boundary")

    module = "nanopores.geometries.%s.params_geo" %name
    params = nanopores.import_vars(module)
    params.update(meta)
    params["name"] = name
    syn = params.pop("synonymes")

    return Geometry(None, mesh, subdomains, boundaries, physdom, physbou, syn, params)

def geo_from_xml_threadsafe(name, reuse_mesh=False, clscale=None, **params):
    import os, mpi4py
    comm = mpi_comm_self()

    if mpi4py.MPI.COMM_WORLD.Get_size() > 1:
        pid = str(mpi4py.MPI.COMM_WORLD.Get_rank())
    else:
        pid = str(os.getpid())

    if not reuse_mesh:
        nanopores.generate_mesh(clscale, name, pid=pid, **params)

    DIR = "%s/%s/mesh/" %(nanopores.DATADIR, name)
    mesh = Mesh(comm, DIR+"mesh%s.xml" %pid)
    subdomains = MeshFunction("size_t", mesh, DIR+"mesh%s_physical_region.xml" %pid)
    boundaries = MeshFunction("size_t", mesh, DIR+"mesh%s_facet_region.xml" %pid)

    with open(DIR+"meta%s.txt" %pid, "r") as f:
        meta = eval(f.read())

    physdom = meta.pop("physical_domain")
    physbou = meta.pop("physical_boundary")

    module = "nanopores.geometries.%s.params_geo" %name
    params = nanopores.import_vars(module)
    params.update(meta)
    params["name"] = name
    syn = params.pop("synonymes")

    return Geometry(None, mesh, subdomains, boundaries, physdom, physbou, syn, params)



