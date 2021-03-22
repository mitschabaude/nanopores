from nanopores.tools.balls import Box, Ball

# form building blocks

#........................R.............................
#                                                     .
#                                                     . 
#              .........l0..........                  .
#              .                   .                  .
#              ._ _______________ _...............    .
#              |D|               |D|     .   .   .    .
#              |D|......l1.......|D|    h1   .   .    .
#              |D|_ ____l2_____ _|D|......   h2  .    .
#              |DDD|_ _______ _|DDD|..........   .    .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|       |DDDDD|             .    .
#       DNA--->|DDDDD|       |DDDDD|           hpore  .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|..l3...|DDDDD|             .    .
#   MEMBRANE   |DDDDD|       |DDDDD|             .    H
#      |       |DDDDD|       |DDDDD|             .    .
#      |       |DDDDD|       |DDDDD|....h4       .    .
#______V_________|DDD|       |DDD|_____.________ .___ .......
#MMMMMMMMMMMMMMMM|DDD|       |DDD|MMMMM.MMMMMMMMM.MMMM.    hmem
#MMMMMMMMMMMMMMMM|DDD|_______|DDD|MMMMM.MMMMMMMMM.MMMM.......
#                .               .                    .
#                .......l4........                    .
#                                                     .
#                                                     .
#                                                     .
#......................................................

zero = [0, 0, 0]

R = 40
H = 70
l0 = 22.5
l1 = 17.5
l2 = 12.5
l3 = 7.5
l4 = 17.5
hpore = 46
hmem = 2.2
h2 = hpore-35. 
h1 = h2-2.5
h4 = 10.

cmem = [0.,0.,-.5*(hpore-hmem)]
c1 = [0.,0.,.5*(hpore-h1)]
c2 = [0.,0.,hpore*.5-h1-(h2-h1)*.5]
cpore = [0.,0.,-.5*h2]
c4 = [0.,0.,-.5*(hpore-h4)]

rMolecule = 2.
x0 = [0.,1.,.5*hpore - rMolecule - 0.2]
lcMolecule = 0.2

reservoir = Box(center=zero, l=R, w=R, h=H)
upperhalf = Box([-R, -R, cmem[2]], [R, R, 0.5*H])

closed_membrane = Box(center=cmem, l=R, w=R, h=hmem)
closed_dna = Box(center=zero, l=l0, w=l0, h=hpore)
enter_1 = Box(center=c1, l=l1, w=l1, h=h1)
enter_2 = Box(center=c2, l=l2, w=l2, h=h2-h1)
enter_3 = Box(center=cpore, l=l3, w=l3, h=hpore-h2)
substract_mem = Box(center=c4, l=l4, w=l4, h=h4)
substract_mem_spanning = Box(center=c4, l=l0, w=l0, h=h4)
substract_dna = substract_mem_spanning - substract_mem
add_bulkfluid = substract_dna - closed_membrane

domain = reservoir
membrane = closed_membrane - substract_mem
membrane_boundary = closed_membrane - closed_dna
dna = closed_dna - enter_1 - enter_2 - enter_3 - substract_dna
pore = enter_1 | enter_2 | enter_3

bulkfluid = reservoir - membrane - closed_dna | add_bulkfluid
bulkfluid_top = bulkfluid & upperhalf
bulkfluid_bottom = bulkfluid - upperhalf

domain.addsubdomains(
    membrane = membrane,
    dna = dna,
    pore = pore,
    bulkfluid_top = bulkfluid_top,
    bulkfluid_bottom = bulkfluid_bottom,
)

molecule = Ball(x0, r=rMolecule, lc=lcMolecule)
domain.addball(molecule, "molecule", "moleculeb")

dnainnerb = enter_1.boundary("front", "back", "left", "right") | enter_2.boundary("front", "back", "left", "right") | enter_3.boundary("front", "back", "left", "right")
dnaupperb = closed_dna.boundary("top") - enter_1
dnaupperb = dnaupperb | enter_1.boundary("bottom") - enter_2
dnaupperb = dnaupperb | enter_2.boundary("bottom") - enter_3
dnaouterb = closed_dna.boundary("front", "back", "left", "right") - \
            substract_mem_spanning.boundary("front", "back", "left", "right")
dnaouterb = dnaouterb | (substract_mem_spanning.boundary("top") - substract_mem.boundary("top"))
dnaouterb = dnaouterb | (substract_mem.boundary("front", "back", "right", "left") - membrane.boundary())
dnalowerb = substract_mem.boundary("bottom") - enter_3
memb = closed_membrane.boundary("top", "bottom") - substract_mem
outermemb = closed_membrane.boundary("front", "back", "left", "right")
sideb = reservoir.boundary("front", "back", "left", "right") - outermemb
upperb = reservoir.boundary("top")
lowerb = reservoir.boundary("bottom")

domain.addboundaries(
    dnainnerb = dnainnerb,
    dnaupperb = dnaupperb,
    dnaouterb = dnaouterb,
    dnalowerb = dnalowerb,
    memb = memb,
    sideb = sideb,
    upperb = upperb,
    lowerb = lowerb,
)

# add synonymes for overlapping subdomains and boundaries
domain.synonymes = dict(
    #subdomains
    bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
    fluid = {"bulkfluid", "pore"},
    solid = {"membrane", "dna", "molecule"},

    #boundaries
    chargeddnab = {"dnaouterb", "dnainnerb"},
    dnab = {"chargeddnab", "dnaupperb", "dnalowerb"},
    noslip = {"dnab", "memb", "moleculeb"},
    bV = "lowerb",
    ground = "upperb",
    nopressure = "upperb",
)

# add parameters (this should include params needed by physics module)
domain.params = dict(
    dim = 3,
    nm = 1.,
    lscale = 1e9,
)

if __name__ == "__main__":
    import dolfin
    from nanopores import plot_sliced
    
    solid = membrane | dna | molecule
    solid.addsubdomains(dna=dna, membrane=membrane)
    solid.addball(molecule, "molecule", "moleculeb")
    solid.addboundaries(
        dnainnerb = dnainnerb,
        dnaupperb = dnaupperb,
        dnaouterb = dnaouterb,
        dnalowerb = dnalowerb,
        memb = memb,
    )  
    print("COMPUTING SOLID")
    solidgeo = solid.create_geometry(lc=2.)
    print(solidgeo)
    
    print("COMPUTING DOMAIN")
    geo = domain.create_geometry(lc=2.)
    print(geo)
    print(geo.params)
    
    plot_sliced(geo)
    dolfin.plot(solidgeo.boundaries, title="boundaries")
    dolfin.interactive()

