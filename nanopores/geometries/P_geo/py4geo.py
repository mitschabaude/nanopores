"""
python script that generates mesh for Howorka geometry
"""

import numpy
from importlib import import_module
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
from warnings import warn
import os
from nanopores import INSTALLDIR

def get_geo(**params):
    """
    writes a 2d geo file for an axissymmetric geometry for PRE paper
    'Effective driving force on DNA inside a solid-state nanopore'
    _________
    |        |
    |  _     |
    | | |____|
    | | |____|
    | |_|    |
    |        |
    |________|

    """


    # for now
    for key, value in params.iteritems():
        if value is not None:
            warn("%s is not yet considered for mesh generation" %key)
    geofid = os.path.join(INSTALLDIR, "geometries","P_geo","mesh.geo")
    geofile = open(geofid, "r")
    geo_code = geofile.read()
    geo_dict = {"gmsh mesh generating sript": __name__,
                "used the following static geo": geofid,
                "geo_code": geo_code,
            }
    return geo_dict

# -----
if __name__ == '__main__':
    print(get_geo())
    print('\n - This is the sample code for the geo file')
