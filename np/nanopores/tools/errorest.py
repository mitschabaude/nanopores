from dolfin import *
import math

__all__ = ["edge_residual_indicator", "poisson_indicator", "zz_indicator",
           "pb_indicator", "Estimator", "pb_indicator_GO"]

class Estimator(object):
    ''' object consisting of pairs (N, f(N)) describing convergence of an error or similar '''
    def __init__(self, name):
        self.name = name
        self.pairs = []
    def __iadd__(self, pair):
        self.pairs.append(pair)
        return self

    # implement list-like behavior:
    def __len__(self):
        return len(self.pairs)
    def __getitem__(self, i):
        return self.pairs[i]
    def __setitem__(self, i, val):
        self.pairs[i] = val
    def __delitem__(self, i):
        del self.pairs[i]
    def __iter__(self):
        return iter(self.pairs)

    def split(self):
        return ([x for x,y in self.pairs], [y for x,y in self.pairs])

    # convergence rates:
    def rates(self):
        # assuming a rate f(N) = c*N^a, return list of a
        tmp = zip(self.pairs[1:],self.pairs[:-1])
        return [math.log(fN/fK)/math.log(float(N)/K) for (N,fN),(K,fK) in tmp]

    # for plotting with matlab:
    def save_to_matlab(self):
        from scipy.io import savemat
        from numpy import array
        dic = {"N": array(self.pairs)[:,0], "err": array(self.pairs)[:,1]}
        savemat("est_%s.mat" %self.name, dic)

    # errorplot using matplotlib
    def plot(self, rate=None, dim=2):
        from matplotlib.pyplot import loglog, xlabel, ylabel, legend, show
        N, err = self.split()
        loglog(N, err, 's-', label="Estimator")
        if rate and N[0] != 0:
            alg = [err[0]/(N[0]**rate)*n**rate for n in N]
            loglog(N, alg, '--', label="N^{%s}" %rate)
        xlabel("# Elements")
        ylabel("Error estimator")
        legend(loc='upper right')
        show()

    # newtonplot using matplotlib
    def newtonplot(self):
        from matplotlib.pyplot import semilogy, xlabel, ylabel, legend, show
        i, err = self.split()
        semilogy(i, err, 's-', label="Estimator")
        xlabel(str(self.name))
        ylabel("Error estimator")
        legend(loc='upper right')
        show()



def zz_indicator(v,flux=None,dx=None):
    """ v is assumed to be scalar and of piece-wise polynomial degree >= 1 """
    V = v.function_space()
    mesh = V.mesh()
    DV = VectorFunctionSpace(mesh, 'DG', V.ufl_element().degree()-1)
    #DV_e = VectorFunctionSpace(mesh, 'DG', V.ufl_element().degree())
    DV_e = VectorFunctionSpace(mesh, 'CG', V.ufl_element().degree())
    DG = FunctionSpace(mesh, 'DG', 0)

    if not flux:
        flux = grad(v)

    # flux recovery
    # TODO: is there a better way to do this??
    # (project is slow and in theory unnessecary)

    g = project(flux,DV)
    g_e = project(g,DV_e)
    #g_e = Function(DV_e)
    #g_e.extrapolate(g)

    if not dx:
        dx = Measure("dx")

    w = TestFunction(DG)
    r = w*inner(g-g_e,g-g_e)*dx

    ind = Function(DG)
    assemble(r,tensor=ind.vector())
    err = errornorm(g_e,g,'L2')/norm(g,'L2')
    return ind,err

def edge_residual_indicator(mesh,flux,force=Constant(0.0)):
    residual = div(flux) + force
    W = FunctionSpace(mesh,"DG",0)
    w = TestFunction(W)
    n = FacetNormal(mesh)
    h = CellSize(mesh)
    r = w*(h*residual)**2*dx + avg(w)*avg(h)*jump(flux,n)**2*dS
    indicators = Function(W)
    b = indicators.vector()
    assemble(r,tensor=b)
    return indicators

def poisson_indicator(geo, u, f=None, cyl=False):
    mesh = geo.mesh
    W = FunctionSpace(mesh,"DG",0)
    w = TestFunction(W)
    n = FacetNormal(mesh)
    h = CellSize(mesh)/geo.physics.lscale

    dS = geo.dS()
    dx = geo.dx()

    Aperm = geo.pwconst("permittivity")
    volcharge = geo.pwconst("volcharge")
    flux = Aperm*geo.physics.grad(u)
    r = Expression("2*pi*x[0]") if cyl else Constant(1.0)

    residual = geo.physics.div(flux) + volcharge
    if f:
        residual = residual + f

    lscale = geo.physics.lscale
    def Clscale(i):
        return Constant( lscale**i )

    res = Clscale(1)*avg(w)*avg(h)*jump(flux,n)**2*r('+')*dS + w*h**2*residual**2*r*dx
    restot = Clscale(1)*avg(h)*jump(flux,n)**2*r('+')*dS + h**2*residual**2*r*dx
    energynorm = inner(flux, geo.physics.grad(u))*r*dx

    indicators = Function(W)
    error = sqrt(assemble(restot)/assemble(energynorm))
    assemble(res,tensor=indicators.vector())
    return indicators, error

def pb_indicator(geo, phys, u, cyl=False):
    c0 = phys.bulkcon
    chi = geo.pwconst("ions", value={"ions":1.,"solid":0.})
    f = Constant(-phys.cFarad*2*c0/phys.UT)*u*chi
    return poisson_indicator(geo, u, f=f, cyl=cyl)

def pb_indicator_GO(geo, phys, u, z, cyl=False):
    # u .. primal solution
    # z .. dual solution
    import numpy

    mesh = geo.mesh
    V = FunctionSpace(mesh, "DG", 0)
    W = FunctionSpace(mesh, "CG", 1)
    EW = FunctionSpace(mesh, "CG", 2)
    Ez = Function(EW)
    Ez.extrapolate(z)
    w = Ez - interpolate(Ez, W)

    v = TestFunction(V)
    n = FacetNormal(mesh)
    h = CellSize(mesh)

    r = Expression("2*pi*x[0]") if cyl else Constant(1.)
    dS = geo.dS() # interior facets
    ds = geo.ds() # exterior facets
    dx0 = geo.dx("ions")
    dx = geo.dx()

    c0 = phys.bulkcon
    cFarad = phys.cFarad
    UT = phys.UT
    eps = geo.pwconst('permittivity')
    def Clscale(i):
        return Constant( (phys.lscale)**i )

    flux = eps*phys.grad(u)

    # local residuals
    rform = -v*phys.div(flux)*w*r*dx \
        +v*Constant(cFarad*2*c0/UT)*u*w*r*dx0 \
        -geo.linearRHS(v*w*r, "volcharge") \
        +Clscale(1)*(
            +avg(v)*jump(n, flux*w)*r('+')*dS \
            +v*inner(n, flux)*w*r*ds \
            -geo.NeumannRHS(v*w*r, "surfcharge"))

    # global residual
    def R(w):
        return assemble(
            inner(flux, phys.grad(w))*r*dx \
            +Constant(cFarad*2*c0/UT)*u*w*r*dx0 \
            -geo.linearRHS(w*r, "volcharge")
            -Clscale(1)*geo.NeumannRHS(w*r, "surfcharge"))

    indicators = Function(V)
    vec = indicators.vector()
    assemble(rform, tensor=vec)
    vec[:] = numpy.abs(vec[:]) #vec[:]*vec[:] #
    #plotind = plot(indicators, title="indicator pb GO", elevate=0.0)

    errormult = 1e12  # Fmult
    error_res = R(z)*(1./phys.lscale**3)
    error_rep = R(Ez)*errormult*(1./phys.lscale**3)
    error_sum = sum(vec)*errormult*(1./phys.lscale**3)

    # FIXME ?
    print "This should be zero (dual global residual) :", error_res
    print "Extrapolated dual residual :", error_rep
    print "indicator sum :", error_sum

    # return indicators, error_rep, error_sum
    return indicators, error_sum
