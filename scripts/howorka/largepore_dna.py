from nanopores import Box, plot_sliced

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
#      |       |DDDDD|       |DDDDD|             .    .
#______V_______|DDDDD|       |DDDDD|____________ .___ .......
#MMMMMMMMMMMMMMMM|DDD|       |DDD|MMMMMMMMMMMMMMM.MMMM.    hmem
#MMMMMMMMMMMMMMMM|DDD|_______|DDD|MMMMMMMMMMMMMMM.MMMM.......
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

cmem = [0.,0.,-.5*(hpore-hmem)]
c1 = [0.,0.,.5*(hpore-h1)]
c2 = [0.,0.,hpore*.5-h1-(h2-h1)*.5]
cpore = [0.,0.,-.5*h2]

reservoir = Box(center=zero, l=R, w=R, h=H)

closed_membrane = Box(center=cmem, l=R, w=R, h=hmem)
closed_dna = Box(center=zero, l=l0, w=l0, h=hpore)
enter_1 = Box(center=c1, l=l1, w=l1, h=h1)
enter_2 = Box(center=c2, l=l2, w=l2, h=h2-h1)
enter_3 = Box(center=cpore, l=l3, w=l3, h=hpore-h2)
substract_mem = Box(center=cmem, l=l4, w=l4, h=hmem)

domain = reservoir
membrane = closed_membrane - substract_mem
membrane_boundary = closed_membrane - closed_dna
dna = closed_dna - enter_1 - enter_2 - enter_3 - membrane

domain.addsubdomains(
    membrane = membrane,
    dna = dna,
    pore = enter_1 | enter_2 | enter_3,
    bulkfluid = reservoir - membrane - closed_dna,
)

dnainnerb = enter_1.boundary("front", "back", "left", "right") | enter_2.boundary("front", "back", "left", "right") | enter_3.boundary("front", "back", "left", "right")
dnaupperb = closed_dna.boundary("top") - enter_1
dnaupperb = dnaupperb | enter_1.boundary("bottom") - enter_2
dnaupperb = dnaupperb | enter_2.boundary("bottom") - enter_3
dnaouterb = closed_dna.boundary("front", "back", "left", "right") - membrane
dnalowerb = substract_mem.boundary("bottom") - enter_3
uppermemb = closed_membrane.boundary("top") - closed_dna
lowermemb = closed_membrane.boundary("bottom") - substract_mem
memb = uppermemb | lowermemb
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
    fluid = {"bulkfluid", "pore"},
    solid = {"membrane", "dna"},

    #boundaries
    chargeddnab = {"dnaouterb", "dnainnerb"},
    dnab = {"chargeddnab", "dnaupperb", "dnalowerb"},
    noslip = {"dnab", "memb"},
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
    geo = domain.create_geometry(lc=1.)
    print geo
    
    plot_sliced(geo)
    
    import dolfin
    dolfin.plot(geo.submesh("solid"))
    dolfin.interactive()
    
    # TODO: doesnt work
    # pure solid domain with boundaries for visual verification
#    solid = membrane | dna
#    solid.addsubdomains(dna=dna, membrane=membrane)
#    solid.addboundaries(
#        dnaouterb = dnaouterb,
#        dnainnerb = dnainnerb,
#        dnaedgeb = dnaedgeb,
#        memb = memb,
#    )
#    geo = solid.create_geometry(lc=1.)
#    solid.plot()
