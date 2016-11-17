# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
J = pugh.F_explicit([[0.,0.,0.]], nproc=1, name="pughopen", h=6.,
                    Nmax=.5e5, cache=False, dim=2)

print J