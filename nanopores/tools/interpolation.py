# (c) 2016 Gregor Mitscha-Baude
"harmonic interpolation."
import numpy as np
import dolfin

# TODO: add points in 3D
# TODO: also interpolate from given values on subdomains and boundaries

class PointBC(object):
    
    def __init__(self, V, points, values, tol=1e-5):
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
        
        self.bc_f.vector().set_local(bc_values)
        self.bc_f.vector().apply("insert") # TODO: what does this do?
        self.dof_set = np.array(dof_map[list(node_set)], dtype="intc")
        #self.bc_f = dolfin.interpolate(function, V)
    
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
            
def harmonic_interpolation(mesh, points, values):
    # Laplace equation
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    a = dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dolfin.dx
    L = dolfin.Constant(0.)*v*dolfin.dx(domain=mesh)
    
    # Point-wise boundary condition
    bc = PointBC(V, points, values)
    
    # Assemble, apply bc and solve
    A, b = dolfin.assemble_system(a, L)
    bc.apply(A)
    bc.apply(b)
    
    u = dolfin.Function(V)
    dolfin.solve(A, u.vector(), b)
    return u

if __name__ == "__main__":    
    from nanopores.geometries import pughpore
    from nanopores.tools.utilities import user_params
    
    h = 1.
    domain = pughpore.get_domain_cyl(h)
    
    # create points
    N = user_params(N=1000).N
    R, H = pughpore.square2circle(domain.params["R"]), domain.params["H"]
    px = R*np.random.random(N)
    py = H*np.random.random(N) - H/2.
    points = zip(px, py)
    
    # prepare mesh containing points
    domain.write_gmsh_code(h)
    domain.insert_points_2D(points, h, forbidden=["dna", "membrane"])
    geo = domain.code_to_mesh()
    mesh = geo.mesh
    
    # interpolate function from point values
    values = np.sin(np.array(points)[:,1]/5.)
    u = harmonic_interpolation(mesh, points, values)
    
    dolfin.plot(mesh)
    dolfin.plot(u)
    dolfin.interactive()
