# (c) 2016 Gregor Mitscha-Baude
"""parameters for pore with spherical molecule(s) inside
(3D or axisymmetric 2D geometries are assumed)"""

#import dolfin
from nanopores.physics.pore import *

Qmol = 0. # molecule charge [q]
Qmolq = lambda Qmol, qq: Qmol*qq
qTarget = Qmolq
rpermMol = rpermDNA
permMol = lambda rpermMol: eperm*rpermMol

def r2pi(dim):
    return dolfin.Expression("2*pi*x[0]") if dim==2 else dolfin.Constant(1.)

volcharge.update(
    molecule = "Moleculeqv",
)
permittivity.update(
    molecule = "permMol",
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
    
def rTarget(geo, lscale):
    return geo.params["rMolecule"]/lscale
    
def invscale(lscale):
    return lambda i: dolfin.Constant(1./lscale**i)

# goal functionals of various continuous quantities
def Fbare(geo, r2pi, Moleculeqv, grad, invscale):
    def _Fel(v, i):
        dx = geo.dx("molecule")
        #scale = dolfin.Constant(lscale**(-3))
        return Moleculeqv*(-r2pi*grad(v)[i])*dx
    return _Fel
    
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
    
# results functionals
def CurrentPNPS(geo, cFarad, UT, grad, r2pi, dim, invscale):
    def _current(U):
        v, cp, cm, u, p = U
        L = dolfin.Constant(geo.params["lpore"])
        Dp = geo.pwconst("Dp")
        Dm = geo.pwconst("Dm")
        cUT = dolfin.Constant(UT)
        F = dolfin.Constant(cFarad)
        
        jm = -Dm*grad(cm) + Dm/cUT*cm*grad(v) + cm*u
        jp = -Dp*grad(cp) - Dp/cUT*cp*grad(v) + cp*u
        jz = F*(jp - jm)[dim-1]
        
        J = -jz/L * r2pi*invscale(2)*geo.dx("pore")
        J = dolfin.assemble(J)
        return dict(J=J)
    return _current

def ForcesPNPS(geo, Moleculeqv, div, grad, r2pi, eta,
               invscale, dim, cFarad, pscale):
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
            Fbarevol = rho0*(-grad(v)[i]) * invscale(3)*dx
            
            ei = tuple((1. if j==i else 0.) for j in range(dim))
            ei = dolfin.Constant(ei)
            uaux = dolfin.Function(V)
            geo.BC(V, ei, "moleculeb").apply(uaux.vector())
            
            Fdragvol = -(-inner(fstokes, uaux) + \
                eta2*inner(sym(grad(u)), sym(grad(uaux))) + \
                div(uaux)*p) *invscale(3)*dxf
                
            F_dict["Fel"].append(dolfin.assemble(Fbarevol))
            F_dict["Fdrag"].append(dolfin.assemble(Fdragvol))
            
        return F_dict
    return _forces

                
    