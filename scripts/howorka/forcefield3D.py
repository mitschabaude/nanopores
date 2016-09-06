"calculate force on given points; clever save/continue depending on params."
from nanopores.models import Howorka
from nanopores.tools import fields
from nanopores import add_params

add_params(
    rMolecule = 0.5,
    Qmol = -1.,
)

params = dict(
    bV = 0.,
    dnaqsdamp = 0.5,
    rMolecule = rMolecule,
    Qmol = Qmol,
    bulkcon = 3e2,
)

solver_params = dict(params,
    Rx = 8.,
    h3D = 8.,
    h = 1.,
    lcpore = 0.1,
    Nmax3D = 2.5e4, # UMFPACK: max. ca. 3e4
    Nmax = 1e4,
    stokesLU = True,
)
fields.update()

X = fields.get_entry("xforce", "X")
N = len(X)
if fields.exists("force3D", **params):
    Xdone = fields.get_field("force3D", "x", **params)
    X = [x0 for x0 in X if x0 not in Xdone]
    print "Existing force file found, %d/%d points remaining." % (len(X), N)
Xfailed = []

for x in X:
    x0 = [x[0], 0., x[1]]
    try:
        F, Fel, Fdrag = Howorka.F_explicit3D([x0], **solver_params)
        fields.save_fields("force3D", params, x=[x], F=F, Fel=Fel, Fdrag=Fdrag)
    except RuntimeError:
        print "RuntimeError occured, continuing without saving."
        Xfailed.append(x)

print "failed:"       
print Xfailed
print "%d of %d force calculations failed." % (len(Xfailed), len(X))

fields.update()