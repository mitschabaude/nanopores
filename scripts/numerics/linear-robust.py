import itertools, numpy
import nanopores
import HoworkaSimple as H

nanopores.add_params(
bVstep = 0.025,
tol = 1e-2,
h = .5,
)

bVs = [ -.01, -.02, -.05, -.1, -.2, -.5, -1., -2.]
qs = [0.1, 0.25, 0.5, 1., 2.]
params = itertools.product(bVs, qs)
M, N = len(bVs), len(qs)

geo, _ = H.setup2D(z0=None, h=h)
solve = ( lambda phys: H.solve2D_fixedpoint(geo, phys, imax=20, tol=tol),
    lambda phys: H.solve2D_fixedpoint_bVscheme(geo, phys, imax=20, bVstep=bVstep, tol=tol),
    lambda phys: H.solve2D_fixedpoint_bVscheme(geo, phys, imax=100, bVstep=bVstep, tol=tol),
    lambda phys: H.solve2D_hybrid(geo, phys, imax=10, tol=tol),
    lambda phys: H.solve2D_hybrid_PB(geo, phys, imax=10, tol=tol)
    )
data = tuple(numpy.zeros([M, N], dtype=bool) for k in range(5))

#for bV, dnaqs in params:
for i in range(M):
    for j in range(N):
        phys = H.phys2D(geo, bV=bVs[i], dnaqsdamp=qs[j])
        for k in range(5):
            data[k][i, j] = solve[k](phys)
            
jsondata = tuple(d.tolist() for d in data)          
nanopores.save_stuff("robustness_fixedpoint", jsondata)

for d in data:
    print
    print d