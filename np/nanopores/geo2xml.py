import subprocess
from importlib import import_module
import os
import nanopores

def generate_mesh(clscale, gid, xml=True, pid="", dim=3, **params):
    """
    python function that writes geo for given geometry and xml for fenics

    Input: clscale... scaling of characteristic length in gmsh [float]
           gid ... geometry identifier [string]
           pid ... optional process id to prevent file access clashes
           ...
    Out: geo_dict... file identifier dictionary + geo_dict
    """

    inputfile = "input%s.geo" %pid
    outfile = "out%s.msh" %pid
    meshfile = "mesh%s.xml" %pid

    py4geo = "geometries.%s.py4geo" %gid
    exec('from %s import get_geo' %py4geo)

    # create path/to/nanoporesdata/gid/mesh if not already there
    meshdir = os.path.join(nanopores.DATADIR, gid, "mesh")
    if not os.path.exists(meshdir):
        os.makedirs(meshdir)

    fid_dict = {"fid_geo": os.path.join(meshdir, inputfile),
                "fid_msh": os.path.join(meshdir, outfile),
    }

    # save code to .geo file
    geo_dict = get_geo(**params)
    fobj = open(fid_dict["fid_geo"], "w")
    fobj.write(geo_dict["geo_code"])
    fobj.close()
    del geo_dict["geo_code"]

    # after writing the geo file, call gmsh
    gmsh_out = subprocess.call(["gmsh", "-%s" %dim, "-v", "1","-clscale", "%f" %clscale,
                     fid_dict["fid_geo"], "-o", fid_dict["fid_msh"]])

    if gmsh_out < 0:
        raise RuntimeError('Gmsh failed in generating this geometry')
    if xml:
        fid_dict["fid_xml"] = os.path.join(meshdir, meshfile)
        subprocess.check_output(["dolfin-convert", fid_dict["fid_msh"], fid_dict["fid_xml"]])
        
    # optionally, write metadata to file ("meta" should be dict-like)
    if "meta" in geo_dict:
        save(geo_dict["meta"], meshdir, "meta")

    geo_dict.update(fid_dict)
    return geo_dict
    
    
def save(data, dir=".", name="file"):
    with open('%s/%s.txt' % (dir,name), 'w') as f:
        f.write(repr(data))


# -----
# to test script run '>> python -m nanopores.geo2xml'
if __name__ == '__main__':
    params = {"x0": None}
    print(generate_mesh(
        clscale=7.0, gid="W_2D_geo", xml=False, **params)
    )
