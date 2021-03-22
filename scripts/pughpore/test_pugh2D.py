# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:42:54 2016

@author: gregor
"""
import os
import dolfin
import nanopores.geometries.pughpore as pugh
import nanopores

nanopores.add_params(
    z = 0.,
    r = 1.,
)
print(os.getpid())
geo = pugh.get_geo_cyl(x0=[1.,1.,z], rMolecule=r) #rMolecule=0.15)
dolfin.plot(geo.subdomains)
dolfin.plot(geo.boundaries)
dolfin.interactive()
