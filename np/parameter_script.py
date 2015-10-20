import nanopores

params = dict(
z0 = 2.,
rMolecule = 0.55,
dnaqsdamp = 0.25,
bV0 = 0.01,
Qmol = nanopores.crange(-3, 3, 1),

refinement = False,
clscale = 0.5,
maxcells = 3e4,
newtondamp = 1.
)

nanopores.simulation2D(**params)


