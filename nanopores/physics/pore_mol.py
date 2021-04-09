# (c) 2016 Gregor Mitscha-Baude
"""parameters for pore with spherical molecule(s) inside
(3D or axisymmetric 2D geometries are assumed)"""

#import dolfin
from nanopores.physics.pore import *

Qmol = 0. # molecule charge [q]
Qmolq = lambda Qmol, qq: Qmol*qq
qTarget = Qmolq
rpermMol = rpermDNA
permMol = lambda eperm, rpermMol: eperm*rpermMol
cyl = lambda dim: True if dim==2 else False
posDTarget = True # if True, position-dep. D used for target mols

# params for unspecific binding
bind_prob = 0.1
bind_time = .04e6 # binding duration [ns]

volcharge.update(
    molecule = "Moleculeqv",
)
permittivity.update(
    molecule = "permMol",
    molecules = "permMol",
)

# geometry-dependent parameters
def Moleculeqv(geo, Qmolq, lscale, r2pi): #
    "Molecule volume charge density [C/nm**3], adapted to discrete volume"
    scale = dolfin.Constant(1.0/lscale**3)
    r = scale*r2pi
    def compute(geo):
        vol = dolfin.assemble(r*geo.dx("molecule"))
        return Qmolq/vol if vol > 0. else 0.
    const = geo.constant("Moleculeqv", compute)
    return const

def rMolecule(geo):
    return geo.params["rMolecule"]

def rTarget(rMolecule, lscale):
    return rMolecule/lscale

def DTargetBulk(rTarget, kT, eta, pi):
    "Stokes-Einstein relation"
    return kT/(6.*pi*eta*rTarget)

def ForceField(geo, grad, eta, qTarget, rTarget, pi):
    def Forces0(v, u, subdomain=None):
        E = -grad(v)
        Fel = dolfin.Constant(qTarget)*E
        Fdrag = dolfin.Constant(6.*pi*eta*rTarget)*u
        F = Fel + Fdrag

        # projecting
        if subdomain is not None:
            mesh = geo.submesh(subdomain)
        else:
            mesh = geo.mesh
        V = dolfin.VectorFunctionSpace(mesh, "CG", 1)
        Fel = dolfin.project(Fel, V)
        Fdrag = dolfin.project(Fdrag, V)
        F = dolfin.project(F, V)
        return F, Fel, Fdrag
    return Forces0

# goal functionals of various continuous quantities
def Fbare(geo, r2pi, Moleculeqv, grad, invscale):
    def _Fel(v, i):
        dx = geo.dx("molecule")
        #scale = dolfin.Constant(lscale**(-3))
        return Moleculeqv*(-r2pi*grad(v)[i])*dx
    return _Fel

# results functionals
def ForcesPNPS(geo, Moleculeqv, div, grad, eta, r2pi,
               invscale, dim, cFarad, pscale):
    "electric and drag forces in 3D based on volume integrals"
    def _forces(U):
        v, cp, cm, u, p = U
        p *= dolfin.Constant(pscale)
        dx = geo.dx("molecule")
        dxf = geo.dx("fluid")
        rho0 = Moleculeqv
        eta2 = dolfin.Constant(2.*eta)
        sym = dolfin.sym
        inner = dolfin.inner
        V = dolfin.VectorFunctionSpace(geo.mesh, "CG", 1)

        Farad = dolfin.Constant(cFarad)
        fstokes = -Farad*(cp - cm)*grad(v)

        F_dict = dict(Fel=[], Fdrag=[])
        for i in range(dim):
            Fbarevol = rho0*(-grad(v)[i]) * r2pi*invscale(3)*dx

            ei = tuple((1. if j==i else 0.) for j in range(dim))
            ei = dolfin.Constant(ei)
            uaux = dolfin.Function(V)
            geo.BC(V, ei, "moleculeb").apply(uaux.vector())

            Fdragvol = -(-inner(fstokes, uaux) + \
                eta2*inner(sym(grad(u)), sym(grad(uaux))) + \
                div(uaux)*p)* r2pi*invscale(3)*dxf

            F_dict["Fel"].append(dolfin.assemble(Fbarevol))
            F_dict["Fdrag"].append(dolfin.assemble(Fdragvol))

        F_dict["F"] = [Fe+Fd for Fe, Fd in zip(F_dict["Fel"], F_dict["Fdrag"])]
        return F_dict
    return _forces

def Fdrag(geo, div, grad, eta, r2pi, invscale, dim, pscale):
    "Drag force on Stokes solution with zero RHS"
    def _force(U):
        u, p = U
        p *= dolfin.Constant(pscale)
        dxf = geo.dx("fluid")
        eta2 = dolfin.Constant(2.*eta)
        sym = dolfin.sym
        inner = dolfin.inner
        V = dolfin.VectorFunctionSpace(geo.mesh, "CG", 1)

        F_dict = dict(Fdrag=[])
        for i in range(dim):
            ei = tuple((1. if j==i else 0.) for j in range(dim))
            ei = dolfin.Constant(ei)
            uaux = dolfin.Function(V)
            geo.BC(V, ei, "moleculeb").apply(uaux.vector())

            Fdragvol = -(eta2*inner(sym(grad(u)), sym(grad(uaux))) + \
                         div(uaux)*p)* r2pi*invscale(3)*dxf

            F_dict["Fdrag"].append(dolfin.assemble(Fdragvol))
        return F_dict
    return _force

def AverageDiffusivity(geo, r2pi, dim, invscale, Dp, Dm, DTargetBulk):
    def _avgdiff(U):
        v, cp, cm, u, p = U

        totalCp = dolfin.assemble(cp * r2pi*invscale(3)*geo.dx("pore"))
        avgDp = dolfin.assemble(Dp[dim-1, dim-1] * cp * r2pi*invscale(3)*geo.dx("pore")) / totalCp

        totalCm = dolfin.assemble(cm * r2pi*invscale(3)*geo.dx("pore"))
        avgDm = dolfin.assemble(Dm[dim-1, dim-1] * cm * r2pi*invscale(3)*geo.dx("pore")) / totalCm

        return dict(avgDp=avgDp/DTargetBulk, avgDm=avgDm/DTargetBulk, D0=DTargetBulk)
    return _avgdiff