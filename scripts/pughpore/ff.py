# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
import nanopores.tools.fields as fields

X = fields.get_entry("pughx", "x", k=3)
result = pugh.F_explicit(X, nproc=6, name="pugh")
print result