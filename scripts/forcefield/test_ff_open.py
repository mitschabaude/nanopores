# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
import folders

# calculate 2D currents
params = dict(name="pughopen", h=1., Nmax=.2e5, dim=2)

data = dict(name="Dpugh", Nmax=1e5, dim=2, r=0.11, h=1.0)
func, mesh = folders.fields.get_functions(**data)
D = func["D"]
rDPore = D([0.,0.])[1]/pugh.nano.Physics().D
print "rel. D", rDPore

J0 = pugh.F_explicit([None], diffusivity_data=data, **params)
J1 = pugh.F_explicit([None], diffusivity_data=None, rDPore=rDPore, **params)
J2 = pugh.F_explicit([None], diffusivity_data=None, rDPore=1., **params)

print "\nOPEN PORE CURRENT (2D):"
print "fully position-dependent diff.:", J0["J"][0]
print "pw. const. diff.:", J1["J"][0]
print "constant diff.:", J2["J"][0]

# calculate 3D currents
params = dict(name="pughopen", h=2., dim=3, stokesiter=False)

data = dict(name="Dpugh", Nmax=2e6, dim=3, r=0.11, h=2.0)
func, mesh = folders.fields.get_functions(**data)
D = func["D"]
rDPore = D([0.,0.,0.])[2]/pugh.nano.Physics().D
print "rel. D", rDPore

J0 = pugh.F_explicit([None], diffusivity_data=data, Nmax=3e5, **params)
J1 = pugh.F_explicit([None], diffusivity_data=None, Nmax=2e5, rDPore=rDPore, **params)
J2 = pugh.F_explicit([None], diffusivity_data=None, Nmax=2e5, rDPore=1., **params)

print "\nOPEN PORE CURRENT (3D):"
print "fully position-dependent diff.:", J0["J"][0]
print "pw. const. diff.:", J1["J"][0]
print "constant diff.:", J2["J"][0]
