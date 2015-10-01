''' test file to analyze influence of axisymmetrical confinement on the validity of stoke's law.

stoke's law relates the drag force on a particle moving in a fluid to its velocity v, radius r and the viscosity eta as follows:

    f_drag = 6*pi*eta*r*v

we solve the stokes equation in a idealized cylindrical pore with the given v as boundary condition on the particle and noslip conditions at the walls. from the resulting velocity field around the particle, we compute f_drag and compare it to the theoretical value.

this script creates plots of the reduced diffusion constant for varying confinements.
'''

from nanopores import *
from dolfin import *

# *** PARAMETERS

eta = 1e-3  # fluid viscosity [Pa*s]
v = (0.,-0.1) # particle velocity [m/s] in cylindrical coordinates (r,z)

l = 1e-6   # length scale [m]
r = 1*l   # particle radius [m]

# range of horizontal/vertical confinement to be analyzed [1/l*m]
#lrange = [10**((n+1.5)/8.5) for n in range(20)] 
lrange = [1.15, 1.23, 1.31, 1.4, 1.5, 1.7, 2, 2.25, 2.6, 3.4, 4.4, 5.8, 7.6, 10, 13, 17, 23, 30, 40, 50, 67, 87, 115, 150, 200]

# fixed horizontal/vertical distance modeling no confinement [1/l*m]
l0 = 200

# ***
# function that performs calculation for given confinement
def stokeslaw(lx,ly):
    print "lx = %.1f, ly = %.1f" %(lx,ly)
    lx *= l
    ly *= l
    
    # set up and plot geometry
    name = "Cyl2D"
    generate_mesh(0.9, name, r=r, l=ly*2, R=lx)
    geo = geo_from_name(name, r=r, l=ly*2, R=lx)
    
    #plot(geo.boundaries)
    print "Number of elements: %d" %geo.mesh.num_cells()
    print "hmin/r: %f" %(geo.mesh.hmin()/r,)

    # set up and solve pde model
    from nanopores import params_physical
    params_physical.eta = eta
    pde = LinearPDE(geo, StokesProblemAxisym, inflow = tuple(-t for t in v), eta=eta)
    pde.solve(verbose=False)

    # obtain drag force from the solutions
    (u, p) = pde.solutions()
    dS = geo.dS("moleculeb")
    n = FacetNormal(geo.mesh)
    r2pi = Expression("2*pi*x[0]")
    strain = eta*2.0*dot(sym(grad(u)),n)
    Fp = assemble((-r2pi*p*n[1])('-')*dS)
    Fs = assemble((-r2pi*strain[1])('-')*dS)
    Fd = Fp + Fs

    # calculate theoretical drag force from stoke's law
    v0 = abs(v[1])
    Fd0 = 6*pi*eta*r*v0
    
    # return factor by which diffusion constant is reduced
    print "D_red = %.4f" %(Fd0/Fd,)
    print
    #pde.visualize()
    return Fd0/Fd
    
    
Dx = []
Dy = []
Dxy = []

for l1 in lrange:
    Dx.append(stokeslaw(l1, l0))
    Dy.append(stokeslaw(l0, l1))
    Dxy.append(stokeslaw(l1, l1))
    
import matplotlib.pyplot as plt

D = {"r":Dx, "z":Dy, "rz":Dxy}

for s in D:
    fname = "confinement_%s.eps" %s
    plt.plot(lrange, D[s], 's-')
    plt.xlabel("distance of confining wall compared to particle radius")
    plt.ylabel("dimensionless diffusion constant D/D_bulk")
    plt.savefig(fname, bbox_inches='tight')
    plt.show()
    
    fname = "confinement_%s_semilogx.eps" %s
    plt.semilogx(lrange, D[s], 's-')
    plt.xlabel("distance of confining wall compared to particle radius")
    plt.ylabel("dimensionless diffusion constant D/D_bulk")
    plt.savefig(fname, bbox_inches='tight')
    plt.show()
    
    fname = "confinement_%s_loglog.eps" %s
    plt.loglog(lrange, D[s], 's-')
    plt.xlabel("distance of confining wall compared to particle radius")
    plt.ylabel("dimensionless diffusion constant D/D_bulk")
    plt.savefig(fname, bbox_inches='tight')
    plt.show()
    #plt.legend(loc='lower right')

