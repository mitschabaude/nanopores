# (c) 2016 Gregor Mitscha-Baude and Benjamin Stadlbauer
"Large DNA pore from Pugh et al. 2016"

from nanopores.tools.balls import Box, Ball, EmptySet, set_tol
import nanopores.tools.box as box
import nanopores.py4gmsh as gmsh

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
    l0 = 18., #22.5,
    l1 = 14., #17.5,
    l2 = 10., #12.5,
    l3 = 6., #7.5,
    l4 = 14., #17.5,
    hpore = 46,
    hmem = 2.2,

    h2 = 46.-35.,
    h1 = 46.-35.-2.5,
    h4 = 10.,

    hnear = 10.,

    rMolecule = 2.0779, # molecular radius of protein trypsin
    x0 = [0., 0., 0.],
    lcMolecule = 0.4, # relative to global mesh size
    center_at_x0 = False,
    center_z_at_x0 = False,
)
# change global settings for mesh generation
#set_tol(None) # faster, may lead to degenerate elements or even gmsh error
set_tol(0.1) # more robust

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
    hnear = _params["hnear"]

    cmem = [0.,0.,-.5*(hpore-hmem)]
    c1 = [0.,0.,.5*(hpore-h1)]
    c2 = [0.,0.,hpore*.5-h1-(h2-h1)*.5]
    cpore = [0.,0.,-.5*h2]
    c4 = [0.,0.,-.5*(hpore-h4)]

    rMolecule = _params["rMolecule"]
    x0 = _params["x0"]
    lcMolecule = lc*_params["lcMolecule"]
    center = zero if (not _params["center_at_x0"] or x0 is None) else x0
    if _params["center_z_at_x0"]:
        center[2] = x0[2]

    hcenter = hpore - h2
    lporecurrent = hcenter/3. # for current calculation

    # form building blocks

    reservoir = Box(center=center, l=2.*R, w=2.*R, h=H)
    upperhalf = Box([-2.*R, -2.*R, cmem[2]], [2.*R, 2.*R, center[2]+0.5*H])

    closed_membrane = Box(center=cmem, l=2.*R, w=2.*R, h=hmem)
    closed_dna = Box(center=zero, l=l0, w=l0, h=hpore)

    enter_1 = Box(center=c1, l=l1, w=l1, h=h1)
    enter_2 = Box(center=c2, l=l2, w=l2, h=h2-h1)

    hporetop = closed_dna.b[2]
    hporebot = closed_dna.a[2]
    nearpore_top = Box([-l0/2, -l0/2, hporetop], [l0/2, l0/2, hporetop+hnear])
    nearpore_bot = Box([-l0/2, -l0/2, hporebot-hnear], [l0/2, l0/2, hporebot])

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
    bulkfluid_top = (bulkfluid & upperhalf) - nearpore_top
    bulkfluid_bottom = (bulkfluid - upperhalf) - nearpore_bot

    if x0 is not None:
        # add molecule
        molecule = Ball(x0, r=rMolecule, lc=lcMolecule)
        domain.addball(molecule, "molecule", "moleculeb")
        # if molecule intersects pore boundary, change pore domain
        epsi = min(lcMolecule, .5) # buffer width of mesh around molecule

        if abs(x0[2] - hporetop) <= rMolecule + epsi:
            bulkentry = poreenter & Box(a=[-l1/2,-l1/2, x0[2]-rMolecule-epsi],
                                   b=[ l1/2, l1/2, hporetop])
            nearpore_top |= bulkentry
            poreenter -= bulkentry

        elif abs(x0[2] - hporebot) <= rMolecule + epsi:
            bulkentry = porebot & Box(a=[-l3/2,-l3/2, hporebot],
                                   b=[ l3/2, l3/2, x0[2]+rMolecule+epsi])
            nearpore_bot |= bulkentry
            porebot -= bulkentry

        if x0[2] >= cpore[2]:
            porecurrent = porebot
            porerest = poreenter | poretop | porectr
            poreenter = EmptySet()
        else:
            porecurrent = poretop
            porerest = porebot | porectr
    else:
        domain.addsubdomain(EmptySet(), "molecule")
        domain.addboundary(EmptySet(), "moleculeb")
        porecurrent = porebot
        porerest = poreenter | poretop | porectr
        poreenter = EmptySet()

    domain.addsubdomains(
        membrane = membrane,
        dna = dna,
        poreenter = poreenter,
        porerest = porerest,
        porecurrent = porecurrent,
        bulkfluid_top = bulkfluid_top,
        bulkfluid_bottom = bulkfluid_bottom,
        nearpore_top = nearpore_top,
        nearpore_bottom = nearpore_bot,
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
        pore = {"porerest", "porecurrent", "poreenter"},
        bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
        nearpore = {"nearpore_bottom", "nearpore_top"},
        poreregion = {"nearpore", "pore"},
        fluid = {"bulkfluid", "pore", "nearpore"},
        solid = {"membrane", "dna", "molecule"},
        ions = {"bulkfluid", "pore", "nearpore"},

        #boundaries
        chargeddnab = {"dnaouterb", "dnainnerb", "dnaupperb", "dnalowerb"},
        dnab = {"chargeddnab"},
        noslip = {"dnab", "memb", "moleculeb"},# "sideb"},
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
        lc = lc,
        lscale = 1e9,
        lporecurrent = lporecurrent,
        lpore = hpore,
    )
    return domain

def square2circle(r):
    "obtain radius of circle from half-sidelength of square with same area"
    from math import pi, sqrt
    return r*2./sqrt(pi)

def get_domain_cyl(lc=1., **newparams):
    # TODO
    set_tol(1e-2)
    _params = dict(params, **newparams)
    zero = [0., 0.]

    R = square2circle(_params["R"])
    H = _params["H"]
    l0 = square2circle(_params["l0"])
    l1 = square2circle(_params["l1"])
    l2 = square2circle(_params["l2"])
    l3 = square2circle(_params["l3"])
    l4 = square2circle(_params["l4"])
    hpore = _params["hpore"]
    hmem = _params["hmem"]
    h2 = _params["h2"]
    h1 = _params["h1"]
    h4 = _params["h4"]
    hnear = _params["hnear"]

    cmem = [0.,-.5*(hpore-hmem)]
    c1 = [0.,.5*(hpore-h1)]
    c2 = [0.,hpore*.5-h1-(h2-h1)*.5]
    cpore = [0.,-.5*h2]
    c4 = [0.,-.5*(hpore-h4)]

    rMolecule = _params["rMolecule"]
    x0 = _params["x0"]
    lcMolecule = lc*_params["lcMolecule"]
    zcenter = 0. if (not _params["center_at_x0"] and\
                     not _params["center_z_at_x0"] or x0 is None) else x0[2]

    hcenter = hpore - h2
    lporecurrent = hcenter/3. # for current calculation

    # form building blocks
    reservoir = Box([0., zcenter-.5*H], [R, zcenter+.5*H])
    upperhalf = Box([-2.*R, cmem[1]], [2.*R, zcenter+0.5*H])

    closed_membrane = Box(center=cmem, l=2.*R, w=hmem)
    closed_dna = Box(center=zero, l=l0, w=hpore)
    enter_1 = Box(center=c1, l=l1, w=h1)
    enter_2 = Box(center=c2, l=l2, w=h2-h1)

    hporetop = closed_dna.b[1]
    hporebot = closed_dna.a[1]
    nearpore_top = Box([-l0/2, hporetop], [l0/2, hporetop+hnear])
    nearpore_bot = Box([-l0/2, hporebot-hnear], [l0/2, hporebot])

    hcross0 = enter_2.a[1]
    hcross1 = cpore[1] + hcenter/6.
    hcross2 = cpore[1] - hcenter/6.
    poretop = Box([-l3/2, hcross1], [l3/2, hcross0])
    porectr = Box([-l3/2, hcross2], [l3/2, hcross1])
    porebot = Box([-l3/2, hporebot], [l3/2, hcross2])
    poreenter = enter_1 | enter_2
    pore = poreenter | poretop | porectr | porebot

    substract_mem = Box(center=c4, l=l4, w=h4)
    substract_mem_spanning = Box(center=c4, l=l0, w=h4)
    substract_dna = substract_mem_spanning - substract_mem
    add_bulkfluid = substract_dna - closed_membrane

    domain = reservoir
    membrane = closed_membrane - substract_mem
    dna = closed_dna - pore - substract_dna

    bulkfluid = (reservoir - (membrane | closed_dna)) | add_bulkfluid
    bulkfluid_top = (bulkfluid & upperhalf) - nearpore_top
    bulkfluid_bottom = (bulkfluid - upperhalf) - nearpore_bot

    if x0 is not None:
        # molecule will be added later with hack
        # if molecule intersects pore boundary, change pore domain
        epsi = min(lcMolecule, .5) # buffer width of mesh around molecule

        if abs(x0[2] - hporetop) <= rMolecule + epsi:
            bulkentry = poreenter & Box(a=[0., x0[2]-rMolecule-epsi],
                                   b=[l1/2, hporetop])
            nearpore_top |= bulkentry
            poreenter -= bulkentry

        elif abs(x0[2] - hporebot) <= rMolecule + epsi:
            bulkentry = porebot & Box(a=[0., hporebot],
                                   b=[ l3/2, x0[2]+rMolecule+epsi])
            nearpore_bot |= bulkentry
            porebot -= bulkentry

        if x0[2] >= cpore[1]:
            porecurrent = porebot
            porerest = poreenter | poretop | porectr
            poreenter = EmptySet()
        else:
            porecurrent = poretop
            porerest = porebot | porectr
    else:
        domain.addsubdomain(EmptySet(), "molecule")
        domain.addboundary(EmptySet(), "moleculeb")
        porecurrent = porebot
        porerest = poreenter | poretop | porectr
        poreenter = EmptySet()

    domain.addsubdomains(
        membrane = membrane,
        dna = dna,
        poreenter = poreenter,
        porerest = porerest,
        porecurrent = porecurrent,
        bulkfluid_top = bulkfluid_top,
        bulkfluid_bottom = bulkfluid_bottom,
        nearpore_top = nearpore_top,
        nearpore_bottom = nearpore_bot,
    )

    dnab = dna.boundary() - membrane.boundary()
    outermemb = closed_membrane.boundary("right")
    memb = membrane.boundary() - dna.boundary() - outermemb
    sideb = reservoir.boundary("right") - outermemb
    upperb = reservoir.boundary("top")
    lowerb = reservoir.boundary("bottom")

    domain.addboundaries(
        dnab = dnab,
        memb = memb,
        sideb = sideb,
        upperb = upperb,
        lowerb = lowerb,
    )

    # add synonymes for overlapping subdomains and boundaries
    domain.synonymes = dict(
        #subdomains
        pore = {"porerest", "porecurrent", "poreenter"},
        bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
        nearpore = {"nearpore_bottom", "nearpore_top"},
        poreregion = {"nearpore", "pore"},
        fluid = {"bulkfluid", "pore", "nearpore"},
        solid = {"membrane", "dna", "molecule"},
        ions = "fluid",

        #boundaries
        chargeddnab = {"dnab"},
        noslip = {"dnab", "memb", "moleculeb", "sideb"},
        bV = "lowerb",
        ground = "upperb",
        bulk = {"lowerb", "upperb"},
        nopressure = "upperb",
    )

    # add parameters (this should include params needed by physics module)
    domain.params = dict(_params,
        name = "pughpore",
        lc = lc,
        dim = 3,
        nm = 1.,
        lscale = 1e9,
        lporecurrent = lporecurrent,
        lpore = hpore,
    )
    return domain

def entity2box(ent):
    intervals = [(f if isinstance(f, tuple) else (f,f)) for f in ent]
    return Box(intervals=intervals)

def get_geo_cyl(lc=1., **newparams):
    domain = get_domain_cyl(lc, **newparams)
    if domain.params["x0"] is not None:
        geo = add_molecule(domain, lc)
    else:
        geo = domain.create_geometry(lc=lc)
    return geo

def add_molecule(dom, lc):
    "hack to insert half molecule at left boundary"
    x0 = [0., dom.params["x0"][2]]
    r = dom.params["rMolecule"]
    lcMolecule = dom.params["lcMolecule"]

    #dom = get_domain_cyl()
    dom.addboundaries(leftb=dom.boundary("left"))
    left = dom.getboundary("leftb")

    mol = EmptySet()
    dom.addsubdomain(mol, "molecule")
    dom.addboundary(mol.boundary() - left, "moleculeb")

    dom.compute_entities()
    dom.compute_boundaries(True)

    edgeinds = list(left.indexsets[1])
    #edgeents = [dom.entities[1][i] for i in edgeinds]
    #print edgeents
    edge = [entity2box(dom.entities[1][i]) for i in edgeinds]
    points = [(x0[0], x0[1]-r), tuple(x0), (x0[0], x0[1]+r)]
    circle = [Box(points[i], points[i+1]) for i in range(len(points)-1)]
    N = len(edge)

    dic = box.multi_box_union(edge + circle)

    # add additional point entities
    for p in dic["entities"][0]:
        if not p in dom.entities[0]:
            dom.entities[0].append(p)

    # add new edge entities and compute replacement
    replace = {i:[] for i in edgeinds}
    circleb = []
    for s, ent in zip(dic["esets"][1], dic["entities"][1]):
        for j in s:
            if j < len(edgeinds): # is old edge
                i = edgeinds[j]
                replace[i].append(ent)
            if j >= len(edgeinds): # belongs to circle
                #print j
                circleb.append(ent)
    #print replace
    for k in replace.keys():
        for i, ent in enumerate(replace[k]):
            if ent in dom.entities[1]:
                j = dom.entities[1].index(ent)
            else:
                dom.entities[1].append(ent)
                j = len(dom.entities[1]) - 1
                #print j, ent
            replace[k][i] = j
    for k, v in replace.items():
        if len(v)==1 and k==v[0]:
            replace.pop(k)
    #print replace
    old = set(replace.keys())
    new = box.union(set(v) for v in replace.values())
    # replace edge indices in boundary
    left.indexsets[1] = left.indexsets[1] - old | new

    # compute left circle boundary
    for i, ent in enumerate(circleb):
        if ent in dom.entities[1]:
            j = dom.entities[1].index(ent)
        else:
            dom.entities[1].append(ent)
            j = len(dom.entities[1]) - 1
        circleb[i] = j
    #print "circle:", circleb

    # gmsh circle
    lcCirc = lcMolecule*lc
    m0, m1 = x0[0], x0[1]
    pcirc = [(m0, m1), (m0, m1-r), (m0+r, m1), (m0, m1+r)]
    dom.entities[0].append(pcirc[2])

    dom.gmsh_entities = [[None for e in k] for k in dom.entities]
    pcirc = [dom.entity_to_gmsh(p, 0, lcCirc) for p in pcirc]

    surfs = [gmsh.Circle([pcirc[1], pcirc[0], pcirc[2]]),
             gmsh.Circle([pcirc[2], pcirc[0], pcirc[3]])]
    dom.gmsh_entities[1] += surfs
    N = len(dom.gmsh_entities[1])
    circlearc = [N-2, N-1]

    for k, v in replace.items():
        removed = False
        for j in list(v):
            #print "adding", j,"to replace"
            if j in circleb:
                replace[k].remove(j)
                removed = True
        if removed:
            replace[k].extend(circlearc)
    for j in circleb:
        if not j in new and not j in replace:
            #print "adding", j,"to replace"
            replace[j] = circlearc

    #print replace
        # replace edge indices sub.boundaries
    for sub in dom.subdomains + dom.boundarysubs:
        iset = sub.bdry().indexset
        orients = sub.bdry().orients
        for i in iset & set(replace.keys()):
            iset.remove(i)
            for j in replace[i]:
                iset.add(j)
                if j in circlearc:
                    orients[j] = -1
                else:
                    orients[j] = orients[i]
                #print sub.name, i, j, orients[j]

    # add edge indices to molecule boundary
    mol.bdry().indexset = set(circleb + circlearc)
    mol.bdry().indexsets[1] = set(circleb + circlearc)
    for i in circleb:
        mol.bdry().orients[i] = -1
    for i in circlearc:
        mol.bdry().orients[i] = 1

    dom.entities_to_gmsh_merge(lc)
    # rebuild boundaries involving balls
    for bou in dom.boundaries:
        bou.indexset = bou.csg.evalsets()[1]

    dom.physical_to_gmsh(True)
    #print gmsh.basic._PHYSSURF
    #print gmsh.basic._PHYSVOL
    dom.geo = box.to_mesh()
    dom.geo.params = dom.params
    if hasattr(dom, "synonymes"):
        dom.geo.import_synonymes(dom.synonymes)

    return dom.geo


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
    from nanopores import plot_sliced, user_params
    up = user_params(params, h=4.)

    geo2D = get_geo_cyl(lc=1., **up)
    dolfin.plot(geo2D.subdomains, title="subdomains")
    dolfin.plot(geo2D.boundaries, title="boundaries")
    print geo2D

    domain = get_domain(**up)
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
    solidgeo = solid.create_geometry(lc=up.h)
    print solidgeo

    print "COMPUTING DOMAIN"
    geo = get_geo(lc=up.h, **up)
    print geo
    print geo.params

    plot_sliced(geo, scalarbar=False)
    dolfin.plot(solidgeo.boundaries, title="boundaries", scalarbar=False)
    dolfin.interactive()

