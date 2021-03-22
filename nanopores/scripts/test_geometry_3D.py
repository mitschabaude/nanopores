from nanopores import *
from dolfin import *

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]
z0 = 5.0*nm

params = dict(
#x0 = None,
x0 = [0.,0.,z0],
R = 15*nm,
Rz = 15*nm,
rMolecule = 0.55*nm,
r0 = 1*nm,
#lcMolecule = nm*0.1,
#moleculeblayer = True,
)

t = Timer("Mesh Generation")
meshgen_dict = generate_mesh(8.0, geo_name, **params)
geo = geo_from_name(geo_name, **params)
'''
for subd in geo._physical_domain:
    submesh = geo.submesh(subd)
    geo_sub = geo_from_subdomains(submesh,
                "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)
    plot(geo_sub.subdomains, title=("subdomains on %s" %subd), elevate=-3e1)
    #plot(submesh, title=("initial mesh on %s" %subd), wireframe=True, elevate=-3e1)
interactive()
'''
print("Boundaries:")
for i in geo._bou2phys:
    print("%d: %s" %(i, str(geo._bou2phys[i])))
    
for subd in geo._physical_domain:
    submesh = geo.submesh(subd)
    geo_sub = geo_from_subdomains(submesh,
                "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)
    plot(geo_sub.boundaries, title=("boundaries on %s" %subd), elevate=-3e1)
    #plot(submesh, title=("initial mesh on %s" %subd), wireframe=True, elevate=-3e1)
interactive()
exit()
