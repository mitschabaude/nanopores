import nanopores

nanopores.simulation2D(
z = 7.5, #nanopores.crange(0, 11, 5),
r = 0.4,
dnaqsdamp = 0.2,
bV = [0.003, 0.03, 0.3],
q = -3,

refinement = True,
clscale = .8,
maxcells = 5e4,
newtondamp = .8,

plot = "bV",
outputs = ["F", "J", "Fdrag", "Fbare"],
nproc = 1,
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
'''

