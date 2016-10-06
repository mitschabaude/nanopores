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
cyl = lambda dim: True if dim==2 else False

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

                
    