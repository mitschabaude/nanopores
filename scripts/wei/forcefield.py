# (c) 2017 Gregor Mitscha-Baude
import dolfin
import nanopores as nano
from nanopores.tools import fields
from nanopores.models import nanopore
from nanopores.tools.poreplots import streamlines
from nanopores.geometries.allpores import WeiPore
from nanopores.tools.polygons import MultiPolygon, Polygon
from nanopores.models.diffusion_interpolation import diffusivity_field
params = nano.user_params(
    geoname = "wei",
    h = 5.,
    Nmax = 4e4,
    dim = 2,
    x0 = None,
    rMolecule = 6.,
    Qmol = -1.,
    bV = -0.5,
)
fields.set_dir_dropbox()

name = "wei_force_ps"
if not fields.exists(name, **params):
    setup = nanopore.Setup(**params)
    _, pnps = nanopore.solve(setup, True)
    v, cp, cm, u, p = pnps.solutions()

    F, Fel, Fdrag = setup.phys.ForceField(v, u, "fluid")
    fields.save_functions(name, params, F=F, Fel=Fel, Fdrag=Fdrag)
    fields.update()

name_D = "wei_D_2D"
if not fields.exists(name_D, **params):
    setup = nanopore.Setup(**params)
    dic = diffusivity_field(setup, r=params["rMolecule"], boundary="poresolidb")
    fields.save_functions(name_D, params, **dic)
    fields.update()

F, Fel, Fdrag = fields.get_functions(name, "F", "Fel", "Fdrag", **params)
D, dist = fields.get_functions(name_D, "D", "dist", **params)

if __name__ == "__main__":
    Fplot = nano.user_param(Fplot=False)
    Dplot = nano.user_param(Dplot=False)

    if Fplot:
        plt = nanopore.Plotter(dim=2)
        plt.plot_vector(F, title="F", interactive=True)
        #plt.plot_vector(Fel, title="Fel")
        #plt.plot_vector(Fdrag, title="Fdrag", interactive=True)
        pore = WeiPore()
        params = nano.Params(pore.default, **params)
        sam, au, sin, _ = pore.polygons(params)
        poly = MultiPolygon(Polygon(sam), Polygon(au), Polygon(sin)).nodes
        streamlines(polygon=poly, rx=100., ry=100., Nx=100, Ny=100,
                            maxvalue=None, F=F)
        nano.showplots()

    if Dplot:
        plt = nanopore.Plotter(dim=2)
        #plt.plot(dist, "dist")
        plt.plot(D[1], "D", interactive=True)

for x in [0.,0.], [0.,1.]:
    print "x", x
    print "D", D(x)
    print "F", F(x)