from dolfin import *
import math, numpy

__all__ = ["edge_residual_indicator", "poisson_indicator", "zz_indicator",
           "pb_indicator", "Estimator", "pb_indicator_GO", "pb_indicator_GO_cheap"]

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
    def plot(self, rate=None, fig=True, style="s-"):
        from matplotlib.pyplot import figure, loglog, xlabel, ylabel, legend, show
        if fig is True:
            figure()
        N, err = self.split()
        loglog(N, err, style, label=self.name)
        if rate and N[0] != 0:
            alg = [err[0]/(N[0]**rate)*n**rate for n in N]
            loglog(N, alg, 'k--', label=r"$O(N^{%.2g})$" %rate)
        #xlabel("# Elements")
        xlabel("degrees of freedom")
        ylabel("rel. error")
        legend(loc='upper right')

    # newtonplot using matplotlib
    def newtonplot(self, fig=True, style="s-"):
        from matplotlib.pyplot import semilogy, xlabel, ylabel, legend, show, figure
        if fig is True:
            figure()
        i, err = self.split()
        semilogy(i, err, style, label=str(self.name))
        xlabel("# iterations")
        ylabel("rel. error")
        legend(loc='upper right')
        #show()



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

    mesh = geo.mesh
    V = FunctionSpace(mesh, "DG", 0)
    W = FunctionSpace(mesh, "CG", 1)
    EW = FunctionSpace(mesh, "CG", 2)
    Ez = Function(EW)
    Ez.extrapolate(z)
    w = Ez - z #interpolate(Ez, W)

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
    def rform(w):
        return -v*phys.div(flux)*w*r*dx \
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
            
    # global functional value
    def J(w):
        return assemble(
            inner(flux, phys.grad(w))*r*dx \
            +Constant(cFarad*2*c0/UT)*u*w*r*dx0)

    indicators = Function(V)
    vec = indicators.vector()
    assemble(rform(w), tensor=vec)
    vec[:] = numpy.abs(vec[:])

    goal = J(z)
    goal_ex = J(Ez)
    # precise relevant scale for error (abs value of functional)
    scale = abs(1./goal) if not goal == 0. else 1e12*(1./phys.lscale**3)
    #scale = 1e12*(1./phys.lscale**3) # rough scale (cheaper)
    
    error_res = abs(R(z))*scale
    error_rep = abs(R(Ez))*scale
    error_sum = sum(vec)*scale
    
    # cheaper estimator without extrapolation
    indicators2 = Function(V)
    vec2 = indicators2.vector()
    assemble(rform(z), tensor=vec2)
    vec2[:] = numpy.abs(vec2[:])
    cheap_sum = sum(vec2)*scale
    #plotind = plot(indicators2, title="indicator pb GO", elevate=0.0, interactive=True)

    # FIXME ?
    print "Goal (dual):", goal
    print "Goal (extrapolated dual):", goal_ex
    print "This should be zero (dual global residual):", error_res
    print "Extrapolated dual residual:", error_rep
    print "indicator sum:", error_sum
    print "ind sum w/o extrap:", cheap_sum

    # return indicators, error_rep, error_sum
    return indicators, error_sum, error_rep, cheap_sum, goal, goal_ex
    
def pb_indicator_GO_cheap(geo, phys, u, z, cyl=False):
    # u .. primal solution
    # z .. dual solution

    mesh = geo.mesh
    V = FunctionSpace(mesh, "DG", 0)
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
    def rform(w):
        return -v*phys.div(flux)*w*r*dx \
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
            
    # global functional value
    def J(w):
        return assemble(
            inner(flux, phys.grad(w))*r*dx \
            +Constant(cFarad*2*c0/UT)*u*w*r*dx0)


    goal = J(z)
    scale = abs(1./goal) # precise relevant scale for error (abs value of functional)
    #scale = 1e12*(1./phys.lscale**3) # rough scale (cheaper)
    
    error_res = abs(R(z))*scale
    
    # cheap estimator without extrapolation
    indicators = Function(V)
    vec = indicators.vector()
    assemble(rform(z), tensor=vec)
    vec[:] = numpy.abs(vec[:])
    error_sum = sum(vec)*scale
    #plotind = plot(indicators2, title="indicator pb GO", elevate=0.0, interactive=True)

    print "Goal (dual):", goal
    print "This should be zero (dual global residual):", error_res
    print "indicator sum (does not make sense as error estimate):", error_sum

    # return indicators, error_rep, error_sum
    return indicators, error_sum, goal
    
def simple_pb_indicator_GO(geo, phys, u, z):
    # u .. primal solution
    # z .. dual solution
    mesh = geo.mesh
    V = FunctionSpace(mesh, "DG", 0)
    EW = FunctionSpace(mesh, "CG", 2)
    Ez = Function(EW)
    Ez.extrapolate(z)
    w = Ez - z #interpolate(Ez, W)

    v = TestFunction(V)
    n = FacetNormal(mesh)

    r = phys.r2pi
    dS = geo.dS() # interior facets
    ds = geo.ds() # exterior facets
    dx0 = geo.dx("ions")
    dx = geo.dx()

    c0 = phys.bulkcon
    cFarad = phys.cFarad
    UT = phys.UT
    eps = geo.pwconst('permittivity')
    k = Constant(cFarad*2.*c0/UT)
    def Clscale(i):
        return Constant( (phys.lscale)**i )

    flux = eps*phys.grad(u)

    # local residuals
    def rform(w):
        return -v*phys.div(flux)*w*r*dx + v*k*u*w*r*dx0 \
        -geo.linearRHS(v*w*r, "volcharge") \
        +Clscale(1)*(
            +avg(v)*jump(n, flux*w)*r('+')*dS \
            +v*inner(n, flux)*w*r*ds \
            -geo.NeumannRHS(v*w*r, "surfcharge"))

    # global residual
    def R(w):
        return assemble(
            inner(flux, phys.grad(w))*r*dx + k*u*w*r*dx0 \
            -geo.linearRHS(w*r, "volcharge")
            -Clscale(1)*geo.NeumannRHS(w*r, "surfcharge"))
            
    # global functional value
    def J(w):
        return assemble(inner(flux, phys.grad(w))*r*dx + k*u*w*r*dx0)

    indicators = Function(V)
    vec = indicators.vector()
    assemble(rform(w), tensor=vec)
    vec[:] = numpy.abs(vec[:])

    # precise relevant scale for error (abs value of functional)
    goal = J(z)
    scale = abs(1./goal) if not goal == 0. else 1e12*(1./phys.lscale**3)
    error_rep = abs(R(Ez))*scale

    return indicators, error_rep
    

