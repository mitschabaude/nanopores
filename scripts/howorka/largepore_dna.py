from nanopores import Box, plot_sliced

# form building blocks

#...........................R...................................
#                                                              .
#                                                              . 
#                 .........l0..........                        .
#                 .                   .                        .
#                 ._ _______________ _.....................    .
#                 |S|               |S|         .   .     .    .
#                 |S|......l1.......|S|        h1   .     .    .
#                 |S|_ ____l2_____ _|S|..........   h2    .    .
#                 |SSS|_ _______ _|SSS|..............     .    .
#                 |SSSSS|       |SSSSS|                   .    .
#                 |SSSSS|       |SSSSS|                   .    .
#                 |SSSSS|       |SSSSS|                 hpore  .
#                 |SSSSS|       |SSSSS|                   .    .
#                 |SSSSS|..l3...|SSSSS|                   .    .
#                 |SSSSS|       |SSSSS|                   .    H
#                 |SSSSS|       |SSSSS|                   .    .
#                 |SSSSS|       |SSSSS|                   .    .
#         ________|SSSSS|       |SSSSS|________......     .    .
#         SSSSSSSSSSSSSS|       |SSSSSSSSSSSSSS    hmem   .    .
#         SSSSSSSSSSSSSS|_______|SSSSSSSSSSSSSS............    .
#         .                                   .                .
#         .................lmem................                .
#                                                              .
#                                                              .
#                                                              .
#...............................................................

zero = [0, 0, 0]

R = 40
H = 70
l0 = 22.5
l1 = 17.5
l2 = 12.5
l3 = 7.5
lmem = 35
hpore = 46
hmem = 2.2
h2 = hpore-35. 
h1 = h2-2.5

cmem = [0.,0.,-.5*(hpore-hmem)]
c1 = [0.,0.,.5*(hpore-h1)]
c2 = [0.,0.,hpore*.5-h1-(h2-h1)*.5]
cpore = [0.,0.,-.5*h2]

reservoir = Box(center=zero, l=R, w=R, h=H)

closed_membrane = Box(center=cmem, l=lmem, w=lmem, h=hmem)
closed_dna = Box(center=zero, l=l0, w=l0, h=hpore)
enter_1 = Box(center=c1, l=l1, w=l1, h=h1)
enter_2 = Box(center=c2, l=l2, w=l2, h=h2-h1)
pore = Box(center=cpore, l=l3, w=l3, h=hpore-h2)

domain = reservoir
membrane = closed_membrane - closed_dna
dna = closed_dna - pore - enter_1 - enter_2

domain.addsubdomains(
    membrane = membrane,
    dna=dna,
    pore=pore,
    bulkfluid = reservoir - membrane - closed_dna,
)


domain.params = dict(
    dim = 3,
    nm = 1.,
    lscale = 1e9,
)

domain.synonymes = dict(
    fluid = {"bulkfluid", "pore"},
    solid = {"membrane", "dna"},
)

##########################################


# build disjoint boundaries
#dnaouterb = closed_dna.boundary("front", "back", "left", "right")
#dnamemb = dnaouterb & closed_membrane
#dnaouterb = dnaouterb - dnamemb
#dnainnerb = pore.boundary()
#dnaedgeb = closed_dna.boundary("front", "back", "left", "right") - pore
#
#outermemb = closed_membrane.boundary("front", "back", "left", "right")
#sideb = reservoir.boundary("front", "back", "left", "right") - outermemb
#memb = closed_membrane.boundary("top", "bottom") - closed_dna
#
#upperb = reservoir.boundary("top")
#lowerb = reservoir.boundary("bottom")
#
#domain.addboundaries(
#    dnaouterb = dnaouterb,
#    dnainnerb = dnainnerb,
#    dnaedgeb = dnaedgeb,
#    memb = memb,
#    sideb = sideb,
#    upperb = upperb,
#    lowerb = lowerb,
#)
#
## add synonymes for overlapping subdomains and boundaries
#domain.synonymes = dict(
#    # subdomains
#    fluid = {"bulkfluid", "pore"},
#    solid = {"membrane", "dna"},
#
#    # boundaries
#    chargeddnab = {"dnaouterb", "dnainnerb"},
#    dnab = {"chargeddnab", "dnaedgeb"},
#    noslip = {"dnab", "memb"},
#    bV = "lowerb",
#    ground = "upperb",
#    nopressure = "upperb",
#)
#
## add parameters (this should include params needed by physics module)
#domain.params = dict(
#    dim = 3,
#    nm = 1.,
#    lscale = 1e9,
#)

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
