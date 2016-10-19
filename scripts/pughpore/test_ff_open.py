# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
J = pugh.F_explicit([None], nproc=1, name="pughopen", h=1., Nmax=1e5)

print J