"simple base class for pore geometries compatible with PNPS"

from nanopores.physics.electrolyte import *

bV = None # voltage bias across pore [V] (None for no enforced bias)

Membraneqs = -0.0 # membrane surface charge density [C/m**2]
DNAqsPure = -qq/nm**2 # = -0.16 # DNA surface charge density [C/m**2]
dnaqsdamp = 1. # DNA charge damping
SiNqs = -0.022
SAMqs = -0.078
ahemqs = 0.

rpermPore = rpermw
rpermProtein = 2. # TODO ?????
# rel. diffusivity in pore depends on pore radius
# => good practice to override default value
rDPore = 0.5

DNAqs = lambda DNAqsPure, dnaqsdamp: DNAqsPure*dnaqsdamp
permPore = lambda eperm, rPermPore: eperm*rpermPore
permProtein = lambda eperm, rpermProtein: eperm*rpermProtein
DPore = lambda D, rDPore: D*rDPore

# default bc for diffusivity
bulkbc = True

# piece-wise boundary conditions
v0 = dict(
    upperb = 0.,
    lowerb = "bV"
)
c0 = dict(
    upperb = "bulkcon",
    lowerb = "bulkcon"
)
cp0 = cm0 = c0

permittivity.update(
    #default = eperm*rpermw,
    bulkfluid = eperm*rpermw,
    nearpore = eperm*rpermw,
    pore = "permPore",
    protein = "permProtein", # for protein pores
    alphahem = "permProtein",
    membrane = eperm*rpermLipid,
)

surfcharge = dict( # surface charge densities for Neumann RHS
    memb = "Membraneqs",
    chargeddnab = "DNAqs",
    sinb = "SiNqs",
    samb = "SAMqs",
    ahemb = "ahemqs",
    alphahemb = "ahemqs",
)

Dpdict = dict(
    nearpore = "D",
    bulkfluid = "D",
    pore = "DPore",
    solid = 0.,
)
Dmdict = Dpdict

# functionals
def CurrentPB(geo, r2pi, bulkcon, mu, rDPore, UT, lscale, cFarad, invscale):
    "approximation of current at 100 mV as linear functional of PB solution"
    bV0 = 0.1
    def J0(v):
        L = geo.params["lpore"]/lscale
        E = bV0/L
        dx = geo.dx("pore")
        Jz = 2*cFarad*bulkcon*mu*rDPore*v/UT*E* r2pi/L*dx
        return Jz
    return J0

def CurrentPNPS(geo, cFarad, UT, grad, r2pi, dim, invscale, Dp, Dm):
    def _current(U):
        v, cp, cm, u, p = U
        L = dolfin.Constant(geo.params["lporecurrent"])
        cUT = dolfin.Constant(UT)
        F = dolfin.Constant(cFarad)

        jm = -Dm*grad(cm) + Dm/cUT*cm*grad(v) + cm*u
        jp = -Dp*grad(cp) - Dp/cUT*cp*grad(v) + cp*u
        jz = F*(jp - jm)[dim-1]

        J = -jz/L * r2pi*invscale(2)*geo.dx("porecurrent")
        J = dolfin.assemble(J)
        return dict(J=J)
    return _current

def CurrentPNPSDetail(geo, cFarad, UT, grad, r2pi, dim, invscale, Dp, Dm):
    def _current(U):
        v, cp, cm, u, p = U
        L = dolfin.Constant(geo.params["lporecurrent"])
        cUT = dolfin.Constant(UT)
        F = dolfin.Constant(cFarad)

        j = dict(
        Jmdif = -Dm*grad(cm),
        Jmmig = Dm/cUT*cm*grad(v),
        Jmconv = cm*u,
        Jpdif = Dp*grad(cp),
        Jpmig = Dp/cUT*cp*grad(v),
        Jpconv = -cp*u
        )

        jform = lambda jz: -F*jz[dim-1]/L * r2pi*invscale(2)*geo.dx("porecurrent")
        Jdict = {key: dolfin.assemble(jform(j[key])) for key in j}
        Jdict["Jm"] = sum(Jdict[s] for s in ["Jmdif", "Jmmig", "Jmconv"])
        Jdict["Jp"] = sum(Jdict[s] for s in ["Jpdif", "Jpmig", "Jpconv"])
        Jdict["Jdif"] = Jdict["Jmdif"] + Jdict["Jpdif"]
        Jdict["Jmig"] = Jdict["Jmmig"] + Jdict["Jpmig"]
        Jdict["Jconv"] = Jdict["Jmconv"] + Jdict["Jpconv"]
        Jdict["J"] = Jdict["Jm"] + Jdict["Jp"]
        return Jdict
    return _current

# surface charges of alphahemolysin
ahemqstotal = [-3.815, 3.1, 16.6, -8.885][::-1]
ahemuniformqs = False
ahemqs0 = lambda ahemqsmulti: ahemqsmulti[0]
ahemqs1 = lambda ahemqsmulti: ahemqsmulti[1]
ahemqs2 = lambda ahemqsmulti: ahemqsmulti[2]
ahemqs3 = lambda ahemqsmulti: ahemqsmulti[3]
surfcharge.update({"alphahemb%d" % i: "ahemqs%d" % i for i in range(4)})

def ahemqsmulti(geo, ahemqstotal, r2pi, invscale, dim, qq, ahemuniformqs):
    if not "alphahemb0" in geo._physical_boundary:
        return [0.]*4
    def area(i):
        return dolfin.assemble(r2pi*invscale(2)*geo.dS("alphahemb%d" % i))
    if not ahemuniformqs:
        qs = [ahemqstotal[i]*qq/area(i) for i in range(4)]
    else:
        qs0 = qq*sum(ahemqstotal)/sum(area(i) for i in range(4))
        qs = [qs0]*4
    print "Calculated alphahem surface charges:", qs
    return qs
