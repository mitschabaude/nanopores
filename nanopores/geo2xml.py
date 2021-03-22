import subprocess
from importlib import import_module
import os
import dolfin
import nanopores
from nanopores.tools.utilities import Log
#FIXME: deprecated because of license conflict -> import from dolfin
#from nanopores.meshconvert import convert2xml
MESHDIR = "/tmp/nanopores"

def geofile2geo(code, meta, name=None, clscale=1.):
    pid = str(os.getpid())
    meshdir = (MESHDIR + "/" + name) if name is not None else MESHDIR

    if not os.path.exists(meshdir):
        os.makedirs(meshdir)

    inputfile = "%s/input%s.geo" % (meshdir, pid)
    outfile = "%s/out%s.msh" % (meshdir, pid)
    meshfile = "%s/mesh%s.xml" % (meshdir, pid)

    xml_sub = "%s/mesh%s_physical_region.xml" % (meshdir, pid)
    xml_bou = "%s/mesh%s_facet_region.xml" % (meshdir, pid)
    if os.path.exists(xml_sub): os.remove(xml_sub)
    if os.path.exists(xml_bou): os.remove(xml_bou)

    with Log("executing gmsh..."):
        # save code to .geo file
        with open(inputfile, "w") as f:
            f.write(code)
        # after writing the geo file, call gmsh
        gmsh_out = subprocess.call(["gmsh", "-3", "-v", "1",
            "-clscale", "%f" %clscale, inputfile, "-o", outfile, "-optimize"])

        if gmsh_out != 0:
            raise RuntimeError("Gmsh failed in generating this geometry")

    with Log("converting to dolfin..."):
        subprocess.check_output(["dolfin-convert", outfile, meshfile])
        # for debugging:
        # convert2xml(outfile, meshfile)
        mesh = dolfin.Mesh(meshfile)

    with open('%s/meta%s.txt' % (meshdir, pid), 'w') as f:
        f.write(repr(meta))

    pdom = meta.pop("physical_domain")
    pbou = meta.pop("physical_boundary")
    subdomains = dolfin.MeshFunction("size_t", mesh, xml_sub) if pdom else None
    boundaries = dolfin.MeshFunction("size_t", mesh, xml_bou) if pbou else None
    geo = nanopores.Geometry(None, mesh, subdomains, boundaries, pdom, pbou)
    return geo

def reconstructgeo(name=None, pid=None, params=None):
    # if pid is None, simply take latest mesh
    # if params is not None, check if they agree with meta["params"]
    # throw error if no matching mesh is available
    meshdir = (MESHDIR + "/" + name) if name is not None else MESHDIR
    if not os.path.exists(meshdir):
        raise EnvironmentError("Geometry folder does not exist yet.")
    if pid is None:
        # get pid of latest mesh
        files = os.listdir(meshdir)
        mfiles = [f for f in files if f.startswith("input")]
        if not mfiles:
            raise EnvironmentError("No existing mesh files found.")
        latest = max(mfiles, key=lambda f: os.path.getmtime(meshdir + "/" + f))
        pid = latest.lstrip("input").rstrip(".geo")

    meshfile = "%s/mesh%s.xml" % (meshdir, pid)
    if not os.path.exists(meshfile):
        raise EnvironmentError(
                  "No existing mesh files found with pid %s." % pid)
    print("Found existing mesh file with pid %s." % pid)

    with open('%s/meta%s.txt' % (meshdir, pid), "r") as f:
        meta = eval(f.read())

    if params is not None:
        if not params == meta["params"]:
            #mparams = meta["params"]
            #print {k: v for k, v in params.items() if k not in mparams or mparams[k] != v}
            #print {k: v for k, v in mparams.items() if k not in params or params[k] != v}
            raise EnvironmentError(
                "Mesh file does not have compatible parameters.")
        print("Mesh file has compatible parameters.")

    print("Reconstructing geometry from %s." % meshfile)

    xml_sub = "%s/mesh%s_physical_region.xml" % (meshdir, pid)
    xml_bou = "%s/mesh%s_facet_region.xml" % (meshdir, pid)

    mesh = dolfin.Mesh(meshfile)
    pdom = meta.pop("physical_domain")
    pbou = meta.pop("physical_boundary")
    subdomains = dolfin.MeshFunction("size_t", mesh, xml_sub) if pdom else None
    boundaries = dolfin.MeshFunction("size_t", mesh, xml_bou) if pbou else None
    geo = nanopores.Geometry(None, mesh, subdomains, boundaries, pdom, pbou)
    return geo


def generate_mesh(clscale, gid, xml=True, pid="", dim=3, optimize=True, **params):
    """
    python function that writes geo for given geometry and xml for fenics

    Input: clscale... scaling of characteristic length in gmsh [float]
           gid ... geometry identifier [string]
           pid ... optional process id to prevent file access clashes
           ...
    Out: geo_dict... file identifier dictionary + geo_dict
    """
    pid = str(os.getpid())
    inputfile = "input%s.geo" %pid
    outfile = "out%s.msh" %pid
    meshfile = "mesh%s.xml" %pid

    py4geo = "nanopores.geometries.%s.py4geo" %gid
    #exec('from %s import get_geo' %py4geo)
    mod = import_module(py4geo)
    get_geo = mod.get_geo

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
    callstr = ["gmsh", "-%s" %dim, "-v", "1","-clscale", "%f" %clscale,
                     fid_dict["fid_geo"], "-o", fid_dict["fid_msh"]]
    if optimize:
        callstr.append("-optimize")
    gmsh_out = subprocess.call(callstr)

    if gmsh_out != 0:
        raise RuntimeError('Gmsh failed in generating this geometry')
    if xml:
        fid_dict["fid_xml"] = os.path.join(meshdir, meshfile)
        subprocess.check_output(["dolfin-convert", fid_dict["fid_msh"], fid_dict["fid_xml"]])
        # for debugging:
        #convert2xml(fid_dict["fid_msh"], fid_dict["fid_xml"])


    # optionally, write metadata to file ("meta" should be dict)
    if "meta" in geo_dict:
        save(geo_dict["meta"], meshdir, "meta%s" %pid)

    geo_dict.update(fid_dict)
    return geo_dict


def save(data, dir=".", name="file"):
    with open('%s/%s.txt' % (dir,name), 'w') as f:
        f.write(repr(data))


# -----
# to test script run '>> python -m nanopores.geo2xml'
if __name__ == '__main__':
    params = {"x0": None}
    print((generate_mesh(
        clscale=7.0, gid="W_2D_geo", xml=False, **params)
    ))
