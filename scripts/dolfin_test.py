# (c) 2016 Gregor Mitscha-Baude
"provide simple dolfin classes for code inspection"
import dolfin

mesh = dolfin.UnitSquareMesh(10, 10)
V = dolfin.FunctionSpace(mesh, "CG", 1)

u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

a = dolfin.inner(dolfin.grad(u), dolfin.grad(v)) * dolfin.dx
L = dolfin.Constant(1.) * v * dolfin.dx

u = dolfin.Function(V)