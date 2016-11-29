# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
import folders

data = dict(name="Dpugh", Nmax=100000.0, dim=2, r=0.11, h=1.0)

J = pugh.F_explicit([None], nproc=1, name="pughopen", h=6.,
                    Nmax=.5e5, cache=False, dim=2, diffusivity_data=data)

print J