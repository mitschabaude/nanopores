import nanopores

nanopores.simulation2D(
z = nanopores.crange(-10, 10, 4),
r = [0.4, 0.55],
dnaqsdamp = 0.2,
bV0 = 0.01,
q = nanopores.crange(-3, 1, 1),

refinement = True,
clscale = .8,
maxcells = 5e4,
newtondamp = .8,

plot = "z",
outputs = ["F", "J", "Fdrag", "Fbare"],
nproc = 3,
)

