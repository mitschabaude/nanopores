''' test file to analyze influence of axisymmetrical confinement on the validity of stoke's law.

stoke's law relates the drag force on a particle moving in a fluid to its velocity v, radius r and the viscosity eta as follows:

    f_drag = 6*pi*eta*r*v

we solve the stokes equation in a idealized cylindrical pore with the given v as boundary condition on the particle and noslip conditions at the walls. from the resulting velocity field around the particle, we compute f_drag and compare it to the theoretical value.
'''

from nanopores import *
from dolfin import *

# *** PARAMETERS

eta = 1e-3  # fluid viscosity [Pa*s]
v = (0.,-0.01) # particle velocity [m/s] in cylindrical coordinates (r,z)

l = 1e-6  # length scale [m]
r = 1*l   # particle radius [m]
lx = 2*l   # horizontal confinement [m]
ly = 2*l   # vertical confinement [m]

# ***
# set up and plot geometry
name = "Cyl2D"
generate_mesh(0.9, name, r=r, l=ly*2, R=lx)
geo = geo_from_name(name, r=r, l=ly*2, R=lx)
#plot(geo.boundaries)
print("Number of elements: %d" %geo.mesh.num_cells())
print("hmin/r: %f" %(geo.mesh.hmin()/r,))

# set up and solve pde model
from nanopores import params_physical as phys
phys.eta = eta
phys.inflow = v

pde = LinearPDE(geo, StokesProblemAxisymEqualOrder, phys=phys)
pde.solve(verbose=False)

# obtain drag force from the solutions
(u, p) = pde.solutions()
dS = geo.dS("moleculeb")
n = FacetNormal(geo.mesh)
r2pi = Expression("2*pi*x[0]")
strain = eta*2.0*dot(sym(grad(u)),n)
Fp = l**2*assemble(Constant(1./l**2)*(-r2pi*p*n[1])('-')*dS)
Fs = l**2*assemble(Constant(1./l**2)*(-r2pi*strain[1])('-')*dS)
Fd = Fp + Fs
print()
print("Forces obtained numerically [N]:")
print("   Fp = %.4e" %(Fp,))
print("   Fs = %.4e" %(Fs,))
print("   Fd = %.4e" %(Fd,))
print("   Fs/Fp = %.4f" %(Fs/Fp,))

# calculate theoretical drag force from stoke's law
v0 = abs(v[1])
Fd0 = 6*pi*eta*r*v0
print()
print("Forces according to Stokes' Law [N]:")
print("   Fp = %.4e" %(Fd0/3.,))
print("   Fs = %.4e" %(Fd0*2./3.,))
print("   Fd = %.4e" %(Fd0,))
print("   Fs/Fp = 2")

print()
print("Diffusion coefficient is reduced due to confinement by the factor")
print("   D_conf/D = %.4f" %(Fd0/Fd,))
print()

# plot solutions
pde.visualize()
#interactive()
