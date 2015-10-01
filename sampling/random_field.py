from numpy import *
from irlb import irlb
from dolfin import Timer, list_timings

# create simple point cloud
t = Timer("initial stuff")
a = arange(0, 1 + 1e-10, 0.025)
A = array([[[x], [y]] for x in a for y in a])
n = a.size
N = A.shape[0]
print N

# translation-invariant correlation function
def c(Z):
    return exp(-5*sqrt(sum(Z**2, 1)))
    #return exp(-100*sum(Z**2, 1))

t.stop()    
# correlation matrix
t = Timer("assemble C")
C = c(A - A.T)
t.stop()
t = Timer("numpy cholesky")
L = linalg.cholesky(C)
t.stop()
t = Timer("numpy full svd")
U0, E0, V0 = linalg.svd(C)
t.stop()
N0 = min(200,N)
t = Timer("truncated svd (%.3f)" %(float(N0)/N))
X = irlb(C, N0)
U, E, V = U0[:,:N0], E0[:N0], (V0.T)[:,:N0]
#U, E, V = X[0], X[1], X[2]
t.stop()
print "Smallest / largest EV of reduced model:",E[-1]/E[0]
print "Sqrt of leftover variance fraction:",sqrt(1 - sum(E)/sum(E0))
#print C
#print (U0[:,N0:]*E0[N0:]).dot((V0)[:,N0:].T)
#print "Sqrt of relative error in matrix norm:",sqrt(linalg.norm(C-(U*E).dot(V.T))/linalg.norm(C))
del C
#print linalg.norm(U0.dot(diag(E0)).dot(V0)-C)/linalg.norm(C)
#print linalg.norm(U.dot(diag(E)).dot(V.T)-C)/linalg.norm(C)

# create sample
t = Timer("create correlated sample")
#z = random.standard_normal((N0,1))
z0 = random.standard_normal((N,1))
z = z0[:N0,:]
u = U.dot(diag(sqrt(E)).dot(z))
u0 = U0.dot(diag(sqrt(E0)).dot(z0))
print "Relative error of sample:",linalg.norm(u-u0)/linalg.norm(u0)
t.stop()

# create white noise to simulate difference
W = sqrt(E[-1])/2*random.randn(N,1)

list_timings()

# plot the field
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

# assumption: E[n] = C*q^n
#q = (E[50]/E[25])**(1./25)
#plt.semilogy(sqrt(E))
#plt.semilogy([q**(k*0.5) for k in range(E.size)])
# assumption: E[n] = C*(n+1)^j, j < 0
j = log(E0[49]/E0[24])/log(2)
j = -1.5 # -(d+1)/d
plt.loglog(sqrt(E0))
plt.loglog([(k+1)**(j*0.5) for k in range(E0.size)])
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = meshgrid(a, a)
#Z0 = U0[:,2].reshape(n,n)
#Z = U[:,10].reshape(n,n)
Z0 = u0.reshape(n,n)
Z = u.reshape(n,n)
Zw = W.reshape(n,n)

# plot reduced and full model sample
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
fig = plt.figure(); ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z0, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
fig = plt.figure(); ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z + Zw, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
# plot difference, noise   
fig = plt.figure(); ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, u.reshape(n,n)-u0.reshape(n,n), rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
fig = plt.figure(); ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Zw, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)

#surf = ax.plot_surface(A[:,0], A[:,1], u, rstride=1, cstride=1, #cmap=cm.coolwarm,
#        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()


