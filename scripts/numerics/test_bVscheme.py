from nanopores import *
import HoworkaSimple as H

add_params(
bV = -.2,
bVstep = 0.025,
tol = 1e-2,
h = .5,
dnaqsdamp = 1.,
)

geo, phys = H.setup2D(z0=None, bV=bV, h=h, dnaqsdamp=dnaqsdamp)
converged0 = H.solve2D_fixedpoint(geo, phys, imax=20, tol=tol)
converged1 = H.solve2D_hybrid(geo, phys, imax=10, tol=tol)
converged2 = H.solve2D_fixedpoint_bVscheme(geo, phys, imax=20, bVstep=bVstep, tol=tol)
converged3 = H.solve2D_hybrid_PB(geo, phys, imax=10, tol=tol)

print "fixed point converged:", converged0
print "hybrid converged:", converged1
print "fixed point with bV damping converged:", converged2
print "hybrid with PB initial guess converged:", converged3