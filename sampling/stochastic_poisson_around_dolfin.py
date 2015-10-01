from numpy import *
from dolfin import *
from numpy import dot, array, exp, sqrt
from irlb import irlb
import matplotlib.pyplot as plt

#mesh = Mesh("dolfin_fine.xml")
mesh = Mesh("dolfin_coarse.xml")
visualize = True

# mesh vertices
N = mesh.num_vertices()
print "Number of vertices:",N
nodes = mesh.coordinates()[:,:,None]

# (translation-invariant) correlation function
def c(Z):
    return exp(-0.5*sqrt(sum(Z**2, 1)))
    #return exp(-1000*sum(Z**2, 1))

# assembly of covariance matrix C
print "Assembling and factorizing covariance matrix..."
t = Timer("assemble C")
C = c(nodes - nodes.T)
t.stop()

# truncated svd
t = Timer("truncated svd")
N = min(200, N)
X = irlb(C, N)
U, E, _ = X[0], X[1], X[2]
K = U*sqrt(E)
t.stop()
print "Smallest / largest EV of reduced model:",E[-1]/E[0]

# cholesky
#t = Timer("cholesky")
#K = linalg.cholesky(C)
#t.stop()

# set up pde problem
V = FunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Function(V)
a = inner(grad(u), grad(v))*dx
L = f*v*dx
u = Function(V)
def boundary(x):
    return near(x[0], 0) or near(x[0], 1) or near(x[1], 0) or near(x[1], 1) 
bc = DirichletBC(V, Constant(0.0), boundary)
#bc = DirichletBC(V, Constant(0.0), DomainBoundary())

# assemble system matrix and define solver
print "Assembling FEM matrix..."
t = Timer("assemble FEM matrix")
A = assemble(a)
bc.apply(A)
pde = LUSolver(A) #, method="default")
pde.parameters["reuse_factorization"] = True
t.stop()

# solve with different samples
n = 100
v2d = vertex_to_dof_map(V)
u_avg = Function(V)
f_avg = Function(V)
u0 = Function(V) # reference avg solution
error = []

for i in range(n):
    if i==0:
        print "Computing LU factorization...\n"
    print "\x1b[A","\r",
    print "Solving... sample %5d" %i

    # create sample
    t = Timer("create sample")
    z = random.randn(N, 1)
    z = K.dot(z)
    t.stop()
    
    # interpolate sample on dolfin space and assemble rhs
    t = Timer("assemble rhs")
    f.vector()[v2d] = z[:,:]
    b = assemble(L)
    bc.apply(b)
    t.stop()
    
    # solve
    t = Timer("%ssolve" % ("LU decomposition + " if i==0 else ""))
    pde.solve(u.vector(), b)
    t.stop()
        
    u_avg.vector()[:] = u_avg.vector()[:] + u.vector()[:]/n
    f_avg.vector()[:] = f_avg.vector()[:] + f.vector()[:]/n
    
    error.append(norm(u_avg, "H1")*n/(i+1))

    
list_timings()
    
# plot average, eigenfunctions
if visualize:
    plot(f, title="example noise")
    plot(u, title="example solution")
    plot(f_avg, title="average noise")
    plot(u_avg, title="average solution")
    if vars().has_key("U"):
        v = []
        for i in range(3):
            v.append(Function(V))
            v[i].vector()[v2d] = ascontiguousarray(U[:,10*i])
            plot(v[i], title="eigenfunction #%d" %(10*i,))
    interactive()
    
plt.loglog(error)
plt.loglog(error[0]*2/sqrt(arange(1,n+1)))
plt.show()

