"simple application-ready framework for 1D meshes, functions, plotting etc."
import dolfin, numpy, nanopores
import matplotlib.pyplot as plt

# FIXME: make Geometry1D subclass Geometry, so we don't have to write geo1D.geo etc.
# that shouldn't be too hard, as there is no overlap, init is the only issue:
# either try assigning to self.__dict__, or write proper copy constructor, or adapt domain.create_geometry()
class Geometry1D(object):
    
    def __init__(self, a=None, b=None, N=100, domain=None):
        if domain is None:
            domain = nanopores.Interval(a, b)
            domain.addsubdomains(fluid = domain)
            domain.addboundaries(upperb = domain.boundary("right"),
                                 lowerb = domain.boundary("left"))
        else:
            a = domain.a[0]
            b = domain.b[0]
        self.a = a
        self.b = b
        self.N = N
        self.domain = domain
        geo = domain.create_geometry(lc=(b-a)/N)
        self.geo = geo
        mesh = geo.mesh
        self.V = dolfin.FunctionSpace(mesh, "CG", 1)
        self.v2d = dolfin.vertex_to_dof_map(self.V)
        self.Z = geo.mesh.coordinates()[:, 0] # numpy array
        self.Zlin = numpy.linspace(a, b, N+1)
        
    def function_from_values(self, values):
        u = dolfin.Function(self.V)
        u.vector()[self.v2d] = numpy.array(values)
        return u
    
    def function_from_lambda(self, f):
        values = [f(z) for z in self.Z]
        return self.function_from_values(values)
    
    def plot(self, u, *args, **kwargs):
        Zlin = self.Zlin
        plt.plot(Zlin, [u(z) for z in Zlin], *args, **kwargs)
    
    def extend_from(self, f, a0, b0, left=0., right=0.):
        "extend function from smaller interval (a0, b0) to (a,b) by constants."
        def ff(z):
            if z > b0:
                return right
            elif z < a0:
                return left
            else:
                return f(z)
        return self.function_from_lambda(ff)
        
if __name__ == "__main__":
    geo = Geometry1D(0, 2*numpy.pi, 50)
    f = geo.function_from_lambda(numpy.sin)
    geo.plot(f, "s-")
    plt.show()
