import nanopores.models.nanopore as nanopore
from nanopores.tools import fields
fields.set_dir_mega()

ddcoupled = dict(name="Dalphahem-coupled", dim=2, Nmax=1e5, h=.5, ahemqsuniform=True, rMolecule=0.11)

params = dict(geoname="alphahem", dim=2, h=1., Nmax=2e4, rDPore=0.3, ahemqs=-0.1, rMolecule=0.11)
params["ahemuniformqs"] = False
params["diffusivity_data"] = ddcoupled
params["x0"] = None
params["bV"] = 0.
setup = nanopore.Setup(**params)
pb, pnps = nanopore.solve(setup, visualize=False)

result = pnps.evaluate(setup.phys.AverageDiffusivity)
print(result)