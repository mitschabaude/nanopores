# (c) 2016 Gregor Mitscha-Baude
from dolfin import *
from nanopores.tools import fields

# set save/load directory
fields.set_dir("/tmp/nanopores/")

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
u.interpolate(Expression("sin(x[0]*x[1]*4*pi)"))

if not fields.exists("test_save"):
    fields.save_functions("test_save", {}, utest=u)
    fields.update()

functions, mesh = fields.get_functions("test_save")
u1 = functions["utest"]
plot(u1, interactive=True)
