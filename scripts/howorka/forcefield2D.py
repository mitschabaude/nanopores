"calculate 2D implicit forcefield; clever save/continue depending on params."
from nanopores.models import Howorka
from nanopores.tools import fields
import nanopores
   
nanopores.add_params(
    rMolecule = 0.5,
    Qmol = -1.,
    h = 1.,
    Nmax = 2e4,
)

params = dict(
    bV = 0.,
    dnaqsdamp = 0.5,
    rMolecule = rMolecule,
    Qmol = Qmol,
    bulkcon = 3e2,
    Nmax = Nmax,
    h = h,
    Rx = 12.,
    Ry = 12.,
)

# save and load implicit force field
NAME = "force2Dimp"

def save_forcefield_implicit(**params):
    F, Fel, Fdrag = Howorka.F_field_implicit(**params)
    mesh = Howorka.geo.mesh
    uid = fields._unique_id()
    FNAME = NAME + uid
    nanopores.save_functions(FNAME, mesh, meta=params, F=F, Fel=Fel, Fdrag=Fdrag)
    fields.save_entries(NAME, params, FNAME=FNAME)
    fields.update()
    
def load_forcefield_implicit(**params):
    FNAME = fields.get_entry(NAME, "FNAME", **params)
    forces, mesh, params = nanopores.load_vector_functions(str(FNAME))
    return forces["F"], forces["Fel"], forces["Fdrag"], mesh, params
    
def maybe_calculate(**newparams):
    params.update(newparams)

    if fields.exists(NAME, **params):
        print "Existing force field found."
    else:
        print "Calculating force field."
        save_forcefield_implicit(**params)
        
    return load_forcefield_implicit(**params)

if __name__ == "__main__":
    import dolfin
    from plot_forcefield import porestreamlines
    F, Fel, Fdrag, mesh, params = maybe_calculate(**params)
    dolfin.plot(mesh)
    porestreamlines(Howorka.polygon(), 6., 8., F=F)
    nanopores.showplots()
    