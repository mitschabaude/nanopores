import nanopores
'''
nanopores.simulation2D(
z = nanopores.crange(0, 11, 5),
r = [0.4, 0.55],
dnaqsdamp = 0.2,
bV = [0.002, 0.003],
q = [-3, 1],

refinement = False,
clscale = .5,
maxcells = 5e4,
newtondamp = .8,

plot = "z",
outputs = ["F", "J", "Fdrag", "Fbare"],
nproc = 4,
)
'''

nanopores.simulate("ahemFPT",
nproc = 1,
plot = "domscale",

domscale = 40, #nanopores.crange(1, 10, 1),
r0 = 5.,
z0 = 10.,

ahemqs = 0.01,
#rTarget = [3e-10, 10e-10, 1e-9]
bV = .5,
bulkcon = 1000.,

clscale = 12.,
skip_stokes = True, #[True, False],
iterative = True,
)


