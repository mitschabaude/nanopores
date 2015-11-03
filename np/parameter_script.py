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

result = nanopores.simulate("ahemFPT",
nproc = 1,
domscale = [1, 2, 3, 4, 5],
z0 = 3.,
plot = "domscale",
)


