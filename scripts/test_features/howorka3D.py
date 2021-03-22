#import nanopores, dolfin
import nanopores.models.Howorka as H

#geo, phys = H.setup3D(h3D=8., lcpore=0.1, x0=[0.2,0,7], Rx = 8.)
#pb, pnps = H.solve3D(geo, phys, Nmax3D=2e3)
#nanopores.plot_sliced(geo)
#dolfin.interactive()

params = dict(
dnaqsdamp = 0.5,
bV = -0.0,
h3D = 8.,
h = 1.,
lcpore = 0.1,
Rx = 8.,
Nmax3D = 3e4, # 8 GB => max. 3e4, 32 GB => 
Nmax = 1e4,
stokesLU = True,
)

z = [7, 4.2, 0]
x = [[0, 0, z0] for z0 in z]

F2, Fel2, Fdrag2 = H.F_explicit(z, **params)
F, Fel, Fdrag = H.F_explicit3D(x, **params)

for i in range(len(z)):
    print("F_z(z=%s)" % z[i])
    print("2D: %.3f, %.3f, %.3f" % (F2[i], Fel2[i], Fdrag2[i]))
    print("3D: %.3f, %.3f, %.3f" % (F[i][2], Fel[i][2], Fdrag[i][2]))
