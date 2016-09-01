import nanopores, dolfin
import nanopores.models.Howorka as H

#geo, phys = H.setup3D(h3D=8., lcpore=0.1, x0=[0.2,0,7], Rx = 8.)
#pb, pnps = H.solve3D(geo, phys, Nmax3D=2e3)
#nanopores.plot_sliced(geo)
#dolfin.interactive()

params = dict(
dnaqsdamp = 0.25,
bV = -0.05,
h3D=8., lcpore=0.1, Rx = 8., Nmax3D=2e3
)

print H.F_explicit3D(x=[[0.2,0,7], [0,0,0]], **params)
