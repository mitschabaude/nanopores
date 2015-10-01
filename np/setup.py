#!/usr/bin/env python

from distutils.core import setup

scripts = ["np-test", "np-test-local", "np-sim"]
packages = ["nanopores", "nanopores.tools", "nanopores.physics", "nanopores.testsuite",
    "nanopores.py4gmsh", "nanopores.geometries", "nanopores.geometries.H_cyl_geo", "nanopores.geometries.H_geo", "nanopores.geometries.P_geo",
    "nanopores.geometries.Nanowire", "nanopores.geometries.W_3D_geo", "nanopores.geometries.W_2D_geo"]

setup(name = "fenics-np",
      version = "1.0",
      description = "Nanopore Simulations with FEniCS",
      author = "Gregor Mitscha-Eibl and Andreas Buttinger-Kreuzhuber",
      author_email = "gregor.mitscha-eibl@tuwien.ac.at",
      url = "https://gitlab.asc.tuwien.ac.at/",
      packages = packages,
      scripts = scripts,
      data_files = [])
