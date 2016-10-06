# (c) 2016 Gregor Mitscha-Baude
"PNPS solvers and visualization for pugh pore"

import dolfin
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import nanopores.physics.simplepnps as simplepnps

def solve1D(geop, physp):
    geo = pughpore.get_geo1D(lc=.001, **geop)
    phys = nano.Physics("pore", geo, **physp)
    pnp = nano.solve_pde(simplepnps.SimplePNPProblem, geo, phys)
    return geo, pnp
    
def visualize1D(geo, pnp):
    v, cp, cm = pnp.solutions()
    h = geo.params["H"]
    nano.plot1D({"potential": v}, (-h/2, h/2, 1001),
                "x", dim=1, axlabels=("z [nm]", "potential [V]"))
    nano.plot1D({"c+": cp, "c-":cm},  (-h/2, h/2, 1001),
                "x", dim=1, axlabels=("z [nm]", "concentrations [mol/m^3]"))
                
class u1D(dolfin.Expression):
    def __init__(self, u):
        self.u = u
        dolfin.Expression.__init__(self)
    def eval(self, value, x):
        value[0] = self.u(x[2])
        
def set_sideBCs(phys, geop, physp):
    geo, pnp = solve1D(geop, physp)
    v, cp, cm = pnp.solutions()
    phys.v0["sideb"] = u1D(v)
    phys.cp0["sideb"] = u1D(cp)
    phys.cm0["sideb"] = u1D(cm)
    