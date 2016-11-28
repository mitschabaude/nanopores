# (c) 2016 Gregor Mitscha-Baude
"provide simple dolfin classes for code inspection"
import dolfin

mesh = dolfin.UnitSquareMesh(10, 10)
V = dolfin.FunctionSpace(mesh, "CG", 1)
VV = V * V

u = dolfin.Function(V)