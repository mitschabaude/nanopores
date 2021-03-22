from nanopores import Box, plot_sliced

# form building blocks
zero = [0, 0, 0]

R = 30
H = 40

hmem = 2
closed_membrane = Box(center=zero, l=R, w=R, h=hmem)

rdna = 25
hdna = 20
cdna = [0, 0, 0.5*(hdna-hmem)]

reservoir = Box(center=cdna, l=R, w=R, h=H)
upperhalf = Box([-R, -R, 0], [R, R, H])

closed_dna = Box(center=cdna, l=rdna, w=rdna, h=hdna)

rpore = 22
pore = Box(center=cdna, l=rpore, w=rpore, h=hdna)

# build domain with disjoint subdomains
domain = reservoir
membrane = closed_membrane - closed_dna
dna = closed_dna - pore

bulkfluid = domain - (dna | membrane | pore)
bulkfluid_top = bulkfluid & upperhalf
bulkfluid_bottom = bulkfluid - upperhalf

domain.addsubdomains(
    membrane = membrane,
    dna = dna,
    pore = pore,
    bulkfluid_top = bulkfluid_top,
    bulkfluid_bottom = bulkfluid_bottom,
    #bulkfluid = reservoir - membrane - dna - pore,
)

# build disjoint boundaries
dnaouterb = closed_dna.boundary("front", "back", "left", "right")
dnamemb = dnaouterb & closed_membrane
dnaouterb = dnaouterb - dnamemb
dnainnerb = pore.boundary("front", "back", "left", "right")
dnaedgeb = closed_dna.boundary("top", "bottom") - pore

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
    bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
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
    merge = True
    
    solid = membrane | dna
    solid.addsubdomains(dna=dna, membrane=membrane)
    solid.addboundaries(
        dnaouterb = dnaouterb,
        dnainnerb = dnainnerb,
        dnaedgeb = dnaedgeb,
        memb = memb,
    )
    
    #print solid
    from dolfin import tic, toc
    tic()
    geo = solid.create_geometry(lc=2., merge=merge)
    print("time:", toc())
    print(geo)
    
    #print dnaedgeb.indexset()
    solid.plot()    

    tic()
    geo = domain.create_geometry(lc=2., merge=merge)
    print("time:", toc())
    #domain.plot()
    #print geo
    
    plot_sliced(geo)
    
    import dolfin
    dolfin.plot(geo.submesh("solid"))
    dolfin.interactive()
    
#    print domain.indexsets
#    
#    for i,en in enumerate(domain.entities[3]):
#        print i,":",en
#    for sub in domain.subdomains:
#        if isinstance(sub, Box):
#            print sub.name, ":", sub.indexsets
#        else:
#            print sub.name, ":", sub.indexset
#    for i,en in enumerate(domain.entities[2]):
#        print i,":",en
#    for sub in domain.boundaries:
#        if isinstance(sub, Box):
#            print sub.name, ":", sub.indexsets
#        else:
#            print sub.name, ":", sub.indexset
            
    
    # TODO: doesnt work
    # pure solid domain with boundaries for visual verification

#    solid.addsubdomains(dna=dna, membrane=membrane)
#    solid.addboundaries(
#        dnaouterb = dnaouterb,
#        dnainnerb = dnainnerb,
#        dnaedgeb = dnaedgeb,
#        memb = memb,
#    )
