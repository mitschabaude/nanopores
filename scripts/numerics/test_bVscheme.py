from nanopores import *
import HoworkaSimple as H

add_params(
bV = -1.,
bVstep = 0.025,
tol = 1e-2,
h = .5,
dnaqsdamp = 1.,
)

geo, phys = H.setup2D(z0=None, bV=bV, h=h)
converged0 = H.solve2D_fixedpoint_bVscheme(geo, phys, imax=20, bVstep=bVstep, tol=tol)
converged1 = H.solve2D_hybrid_PB(geo, phys, imax=15, tol=tol)

print "fixed point converged:", converged0
print "hybrid converged:", converged1