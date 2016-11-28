"testing new dolfin version 2016.1.0"
from dolfin import *
mesh = UnitSquareMesh(2, 2)
V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
l = Constant(1.)*v*dx

u = Function(V)
bc = DirichletBC(V, Constant(0.), DomainBoundary())

solve(a == l, u, bc)
plot(u)
for i in range(5):
    oldmesh = mesh
    mesh = adapt(mesh)

u = adapt(u, mesh, False)
bc = adapt(bc, mesh, V)

import dolfin.cpp as cpp
from dolfin.fem.assembling import _create_tensor

def assemble_cppform(form, tensor=None):
    # Create tensor
    comm = form.mesh().mpi_comm()
    tensor = _create_tensor(comm, form, form.rank(), None, tensor)
    # Call C++ assemble function
    assembler = cpp.Assembler()
    assembler.assemble(tensor, form)
    return tensor

def adapt_form(form, mesh):
    form = Form(form)
    adapted_form = adapt(form, mesh)
    adapted_form._compiled_form = form._compiled_form
    return adapted_form
#
#a = adapt(Form(a), mesh)
#l = adapt(Form(l), mesh)
#A = assemble_cppform(a)
#b = assemble_cppform(l)

a = adapt_form(a, mesh)
l = adapt_form(l, mesh)
A = assemble(a)
b = assemble(l)

bc.apply(A)
bc.apply(b)

solve(A, u.vector(), b)
plot(u)
interactive()
