from nanopores import Box, plot_sliced

# form building blocks
zero = [0, 0, 0]

R = 30
H = 40

hmem = 2
closed_membrane = Box(center=zero, l=R, w=R, h=hmem)

rdna = 20
hdna = 20
cdna = [0, 0, 0.5*(hdna-hmem)]

reservoir = Box(center=cdna, l=R, w=R, h=H)
closed_dna = Box(center=cdna, l=rdna, w=rdna, h=hdna)

rpore = 10
pore = Box(center=cdna, l=rpore, w=rpore, h=hdna)

# build domain with disjoint subdomains
domain = reservoir
membrane = closed_membrane - closed_dna
dna = closed_dna - pore

domain.addsubdomains(
    membrane = membrane,
    dna = dna,
    pore = pore,
    bulkfluid = reservoir - membrane - dna - pore,
)

# build disjoint boundaries
dnaouterb = closed_dna.boundary("front", "back", "left", "right")
dnamemb = dnaouterb & closed_membrane
dnaouterb = dnaouterb - dnamemb
dnainnerb = pore.boundary()
dnaedgeb = closed_dna.boundary("front", "back", "left", "right") - pore

outermemb = closed_membrane.boundary("front", "back", "left", "right")
sideb = reservoir.boundary("front", "back", "left", "right") - outermemb
memb = closed_membrane.boundary("top", "bottom") - closed_dna

upperb = reservoir.boundary("top")
lowerb = reservoir.boundary("bottom")

domain.addboundaries(
    dnaouterb = dnaouterb,
    dnainnerb = dnainnerb,
    dnaedgeb = dnaedgeb,
    memb = memb,
    sideb = sideb,
    upperb = upperb,
    lowerb = lowerb,
)

# add synonymes for overlapping subdomains and boundaries
domain.synonymes = dict(
    # subdomains
    fluid = {"bulkfluid", "pore"},
    solid = {"membrane", "dna"},

    # boundaries
    chargeddnab = {"dnaouterb", "dnainnerb"},
    dnab = {"chargeddnab", "dnaedgeb"},
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