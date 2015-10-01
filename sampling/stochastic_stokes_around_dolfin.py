from numpy import *
from dolfin import *
from numpy import dot, array, exp, sqrt
from irlb import irlb
import matplotlib.pyplot as plt

mesh = Mesh("dolfin_fine.xml")
#mesh = Mesh("dolfin_coarse.xml")
visualize = True
plot_efs = False

# mesh vertices
N = mesh.num_vertices()
print "Number of vertices:",N
nodes = mesh.coordinates()[:,:,None]

# (translation-invariant) correlation function
def c(Z):
    #return 10*exp(-0.5*sqrt(sum(Z**2, 1)))
    return 10*exp(-10*sum(Z**2, 1))

# assembly of covariance matrix C
print "Assembling and factorizing covariance matrix..."
t = Timer("assemble C")
C = c(nodes - nodes.T)
t.stop()

# truncated svd
t = Timer("truncated svd")
N = min(50, N)
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
W = VectorFunctionSpace(mesh, "CG", 2)
VV = VectorFunctionSpace(mesh, "CG", 1)
V1 = FunctionSpace(mesh, "CG", 1)
Q = FunctionSpace(mesh, "CG", 1)
V = W * Q
(u, p) = TrialFunctions(V)
(v, q) = TestFunctions(V)
f = Function(VV)
f1 = Function(V1)
f2 = Function(V1)

a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
L = (f1*v[0] + f2*v[1])*dx
w = Function(V)

class noslip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (
               (between(x[0], (0.1,0.9)) and between(x[1], (0.1,0.9))))
               #or near(x[1], 0) or near(x[1], 1) or \
def inflow(x):
    return near(x[0], 1)
def nopressure(x):
    return near(x[0], 0)
    
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0)), noslip()),
       DirichletBC(V.sub(0), Constant((-0.1, 0.0)), inflow),
       DirichletBC(V.sub(1), Constant(0.0), nopressure)]

# assemble system matrix and define solver
print "Assembling FEM matrix..."
t = Timer("assemble FEM matrix")
A = assemble(a)
for bc in bcs: bc.apply(A)
pde = LUSolver(A) #, method="default")
pde.parameters["reuse_factorization"] = True
t.stop()

# compute deterministic solution (f=0)
print "Computing LU factorization...","\n"
w0 = Function(V)
b = assemble(L)
for bc in bcs: bc.apply(b)
Timer("LU decomposition + solve")
pde.solve(w0.vector(), b)
(u0, p0) = w0.split()
t.stop()

# solve with different samples
n = 1000
v2d = vertex_to_dof_map(V1)
w_avg = Function(V)
f_avg = Function(VV)
f1_avg = Function(V1)
f2_avg = Function(V1)
w_diff = Function(V)
(u_diff, p_diff) = w_diff.split()
error = []

for i in range(n):
    print "\x1b[A","\r",
    print "Solving... sample %5d" %i

    # create sample
    t = Timer("create sample")
    z = random.randn(N, 1)
    z = K.dot(z)
    t.stop()
    
    # interpolate sample on dolfin space and assemble rhs
    t = Timer("assemble rhs")
    f2.vector()[v2d] = z
    b = assemble(L)
    for bc in bcs: bc.apply(b)
    t.stop()
    
    # solve
    t = Timer("solve")
    pde.solve(w.vector(), b)
    t.stop()
        
    w_avg.vector()[:] = w_avg.vector()[:] + w.vector()[:]/n
    f1_avg.vector()[:] = f1_avg.vector()[:] + f1.vector()[:]/n
    f2_avg.vector()[:] = f2_avg.vector()[:] + f2.vector()[:]/n
    
    w_diff.vector()[:] = w_avg.vector()[:]*n/(i+1) - w0.vector()[:]
    error.append(norm(u_diff, "L2"))
    
assign(f, [f1,f2])
assign(f_avg, [f1_avg,f2_avg])
(u, p) = w.split()
(u_avg, p_avg) = w_avg.split()
    
list_timings()

# plot average, eigenfunctions
if visualize:
    plot(f, title="example noise")
    plot(u, title="example velocity")
    plot(f_avg, title="average noise")
    plot(u_avg, title="average velocity")
    plot(u0, title="deterministic velocity")
    if vars().has_key("U") and plot_efs:
        v = []
        for i in range(5):
            v.append(Function(V1))
            v[i].vector()[v2d] = ascontiguousarray(U[:,10*i])
            plot(v[i], title="eigenfunction #%d" %(10*i,))
    interactive()
    
plt.loglog(error)
plt.loglog(error[0]*2/sqrt(arange(1,n+1)))
plt.show()


