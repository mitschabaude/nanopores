import nanopores

nanopores.simulation2D(
z = [0.,7.5],
r = 0.4,
dnaqsdamp = 0.1,
bV = 0.01,
q = -3,

refinement = False,
clscale = 2.,
maxcells = 5e4,
newtondamp = .8,

write_files = False,
new_mesh = False,
)
