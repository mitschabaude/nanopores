# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
from folders import fields

X = fields.get_entry("pughx", "x", k=3)
result = pugh.F_explicit(X, nproc=4, name="pugh")
print result
import nanopores
nanopores.save_stuff("pugh_backup", result)