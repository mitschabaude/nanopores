# (c) 2016 Gregor Mitscha-Baude and Benjamin Stadlbauer
"Large DNA pore from Pugh et al. 2016"

from nanopores.tools.balls import Box, Ball, EmptySet, set_tol

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

params = dict(
    R = 20.,
    H = 70.,
    l0 = 22.5,
    l1 = 17.5,
    l2 = 12.5,
    l3 = 7.5,
    l4 = 17.5,
    hpore = 46,
    hmem = 2.2,
    
    h2 = 46.-35.,
    h1 = 46.-35.-2.5,
    h4 = 10.,
    
    rMolecule = 2.0779, # molecular radius of protein trypsin
    x0 = [0., 0., -20.],
    lcMolecule = 0.4, # relative to global mesh size
)
# change global settings for mesh generation
#set_tol(None) # faster, may lead to degenerate elements or even gmsh error
set_tol(0.5) # more robust

def set_params(**newparams):
    params.update(newparams)
    
def get_geo(lc=1., **newparams):
    domain = get_domain(lc, **newparams)
    geo = domain.create_geometry(lc=lc)
    return geo
    
def get_domain(lc=1., **newparams):
    _params = dict(params, **newparams)
    zero = [0., 0., 0.]
    
    R = _params["R"]
    H = _params["H"]
    l0 = _params["l0"]
    l1 = _params["l1"]
    l2 = _params["l2"]
    l3 = _params["l3"]
    l4 = _params["l4"]
    hpore = _params["hpore"]
    hmem = _params["hmem"]
    h2 = _params["h2"]
    h1 = _params["h1"]
    h4 = _params["h4"]
    
    cmem = [0.,0.,-.5*(hpore-hmem)]
    c1 = [0.,0.,.5*(hpore-h1)]
    c2 = [0.,0.,hpore*.5-h1-(h2-h1)*.5]
    cpore = [0.,0.,-.5*h2]
    c4 = [0.,0.,-.5*(hpore-h4)]
    
    rMolecule = _params["rMolecule"]
    x0 = _params["x0"]
    lcMolecule = lc*_params["lcMolecule"]
    
    hcenter = hpore - h2    
    lporecurrent = hcenter/3. # for current calculation
    
    # form building blocks
    reservoir = Box(center=zero, l=2.*R, w=2.*R, h=H)
    upperhalf = Box([-2.*R, -2.*R, cmem[2]], [2.*R, 2.*R, 0.5*H])
    
    closed_membrane = Box(center=cmem, l=2.*R, w=2.*R, h=hmem)
    closed_dna = Box(center=zero, l=l0, w=l0, h=hpore)
    enter_1 = Box(center=c1, l=l1, w=l1, h=h1)
    enter_2 = Box(center=c2, l=l2, w=l2, h=h2-h1)

    hporetop = closed_dna.b[2]
    hporebot = closed_dna.a[2]
    hcross0 = enter_2.a[2]    
    hcross1 = cpore[2] + hcenter/6.
    hcross2 = cpore[2] - hcenter/6.
    poretop = Box([-l3/2, -l3/2, hcross1], [l3/2, l3/2, hcross0])
    porectr = Box([-l3/2, -l3/2, hcross2], [l3/2, l3/2, hcross1])
    porebot = Box([-l3/2, -l3/2, hporebot], [l3/2, l3/2, hcross2])
    poreenter = enter_1 | enter_2
    pore = poreenter | poretop | porectr | porebot
    
    enter_3 = Box(center=cpore, l=l3, w=l3, h=hcenter)
    substract_mem = Box(center=c4, l=l4, w=l4, h=h4)
    substract_mem_spanning = Box(center=c4, l=l0, w=l0, h=h4)
    substract_dna = substract_mem_spanning - substract_mem
    add_bulkfluid = substract_dna - closed_membrane
    
    domain = reservoir
    membrane = closed_membrane - substract_mem
    dna = closed_dna - pore - substract_dna
    
    bulkfluid = (reservoir - (membrane | closed_dna)) | add_bulkfluid
    bulkfluid_top = bulkfluid & upperhalf
    bulkfluid_bottom = bulkfluid - upperhalf
    
    if x0 is not None:
        # add molecule
        molecule = Ball(x0, r=rMolecule, lc=lcMolecule)
        domain.addball(molecule, "molecule", "moleculeb")
        # if molecule intersects pore boundary, change pore domain
        epsi = min(lcMolecule, .5) # buffer width of mesh around molecule
        
        if abs(x0[2] - hporetop) <= rMolecule + epsi:
            bulkentry = poreenter & Box(a=[-l1/2,-l1/2, x0[2]-rMolecule-epsi],
                                   b=[ l1/2, l1/2, hporetop])
            bulkfluid_top |= bulkentry
            poreenter -= bulkentry
            
        elif abs(x0[2] - hporebot) <= rMolecule + epsi:
            bulkentry = porebot & Box(a=[-l3/2,-l3/2, hporebot],
                                   b=[ l3/2, l3/2, x0[2]+rMolecule+epsi])
            bulkfluid_bottom |= bulkentry
            porebot -= bulkentry
            
        porecurrent = "porebot" if x0[2] >= cpore[2] else "poretop"
    else:
        domain.addsubdomain(EmptySet(), "molecule")
        domain.addboundary(EmptySet(), "moleculeb")
        porecurrent = "poretop"
    
    domain.addsubdomains(
        membrane = membrane,
        dna = dna,
        poreenter = poreenter,
        poretop = poretop,
        porectr = porectr,
        porebot = porebot,
        bulkfluid_top = bulkfluid_top,
        bulkfluid_bottom = bulkfluid_bottom,
    )

    dnainnerb = enter_1.boundary("front", "back", "left", "right") |\
                enter_2.boundary("front", "back", "left", "right") |\
                enter_3.boundary("front", "back", "left", "right")
    dnaupperb = closed_dna.boundary("top") - enter_1
    dnaupperb = dnaupperb | enter_1.boundary("bottom") - enter_2
    dnaupperb = dnaupperb | enter_2.boundary("bottom") - enter_3
    dnaouterb = closed_dna.boundary("front", "back", "left", "right") - \
                substract_mem_spanning.boundary("front", "back", "left", "right")
    dnaouterb = dnaouterb | (substract_mem_spanning.boundary("top") - \
                substract_mem.boundary("top"))
    dnaouterb |= (substract_mem.boundary("front", "back", "right", "left")\
                  - membrane.boundary())
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
        pore = {"poreenter", "poretop", "porectr", "porebot"},
        bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
        fluid = {"bulkfluid", "pore"},
        solid = {"membrane", "dna", "molecule"},
        ions = "fluid",
        porecurrent = porecurrent,
    
        #boundaries
        chargeddnab = {"dnaouterb", "dnainnerb", "dnaupperb", "dnalowerb"},
        dnab = {"chargeddnab"},
        noslip = {"dnab", "memb", "moleculeb", "sideb"},
        bV = "lowerb",
        ground = "upperb",
        bulk = {"lowerb", "upperb"},
        nopressure = "upperb",
    )
    
    # add parameters (this should include params needed by physics module)
    domain.params = dict(_params,
        name = "pughpore",
        dim = 3,
        nm = 1.,
        lscale = 1e9,
        lporecurrent = lporecurrent,
        lpore = hpore,
    )
    return domain
    
def get_geo1D(lc=0.01, **newparams):
    _params = dict(params, **newparams)
    H = _params["H"]
    hmem = _params["hmem"]
    hpore = _params["hpore"]
    cmem = -.5*(hpore-hmem)
    
    domain = Box([-.5*H], [.5*H])
    upperhalf = Box([cmem], [.5*H])
    membrane = Box([cmem - .5*hmem], [cmem + .5*hmem])
    bulkfluid = domain - membrane
    
    lowerb = domain.boundary("left")
    upperb = domain.boundary("right")
    
    domain.addsubdomains(
        membrane = membrane,
        bulkfluid_top = bulkfluid & upperhalf,
        bulkfluid_bottom = bulkfluid - upperhalf,
    )
    domain.addboundaries(
        lowerb = lowerb,
        upperb = upperb,
        memb = membrane.boundary(),
    )
    domain.params["lscale"] = 1e9
    domain.synonymes = dict(
        #subdomains
        bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
        fluid = {"bulkfluid"},
        solid = {"membrane"},
        ions = "fluid",
    
        #boundaries
        noslip = "memb", # "upperb", "sideb", "lowerb"},
        bV = "lowerb",
        ground = "upperb",
        bulk = {"lowerb", "upperb"},
        nopressure = "upperb",
    )
    domain.params = dict(_params,
        name = "pughpore",
        dim = 1,
        nm = 1.,
        lscale = 1e9,
    )
    geo = domain.create_geometry(lc=lc)
    return geo
    
    
if __name__ == "__main__":
    import dolfin
    from nanopores import plot_sliced
    
    domain = get_domain()
    membrane = domain.getsubdomain("membrane")
    dna = domain.getsubdomain("dna")
    molecule = domain.getsubdomain("molecule")
    
    solid = membrane | dna | molecule
    solid.addsubdomains(dna=dna, membrane=membrane)
    solid.addball(molecule, "molecule", "moleculeb")
    solid.addboundaries(
        dnainnerb = domain.getboundary("dnainnerb"),
        dnaupperb = domain.getboundary("dnaupperb"),
        dnaouterb = domain.getboundary("dnaouterb"),
        dnalowerb = domain.getboundary("dnalowerb"),
        memb = domain.getboundary("memb"),
    )  
    print "COMPUTING SOLID"
    solidgeo = solid.create_geometry(lc=4.)
    print solidgeo
    
    print "COMPUTING DOMAIN"
    geo = get_geo(lc=4.)
    print geo
    print geo.params
    
    plot_sliced(geo)
    dolfin.plot(solidgeo.boundaries, title="boundaries")
    dolfin.interactive()

