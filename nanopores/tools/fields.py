"""
storing/retrieving of (position-dependent) fields obtained from simulations.

a field is a collection of functions x -> f(x; params), and will be stored
simply as a number of lists with matching lengths, including the list of
positions x, alongside with a name and a uniquely identifying
parameter dict params.

essentially, this is a self-hacked implementation of a database
for a very specific purpose.

achievements of this module:
-) save fields in hidden database
-) add additional data points x, f(x) at any time
-) automatically identify field based on the given name and params
-) if parameters are superset of existing parameter set,
   assume they match and update existing field
-) retrieve field values and handle non-existence
-) display list of available fields
-) reading AND saving fields can be done in parallel, only update() step not

TODO: if demanded, overwrite already calculated values
TODO: in remove/save_functions, remove obsolete dolfin and .txt files
TODO: fields.list_all (with indices)
TODO: fields.list[i]
"""
import os, json
from nanopores.dirnames import DATADIR, INSTALLDIR, HOME
DIR = os.path.join(DATADIR, "fields")
#DIR = os.path.join(os.path.dirname(INSTALLDIR), "data", "fields")
HEADER = "header.txt"
SUFFIX = ".field.txt"

def set_dir(NEWDIR):
    global DIR
    DIR = NEWDIR
    # assert directory and header exists
    if not os.path.exists(DIR):
        os.makedirs(DIR)
    if not HEADER in os.listdir(DIR):
        _save(dict(_flist=[]), HEADER)
    
# user interface that wraps Header object
def update():
    Header().update()
    
def load_file(name, **params):
    return Header().load_file(name, params)
    
def get_fields(name, **params):
    return Header().get_fields(name, **params)
    

    
def get_field(name, field, **params):
    return Header().get_field(name, field, **params)
    
def save_fields(name, params, **fields):
    # fields has to be a dictionary of lists, x a list
    data = dict(name = name, params = params, fields = fields)
    FILE = name + _unique_id() + SUFFIX
    _save(data, FILE)
    
def save_entries(name, params, **entries):
    data = dict(name = name, params = params, **entries)
    FILE = name + _unique_id() + SUFFIX
    _save(data, FILE)        

def remove(name, **params):
    Header().remove(name, params)
    
def purge(name, **params):
    while exists(name, **params):
        remove(name, **params)
    
def get_entry(name, entry, **params):
    return Header().get_entry(name, entry, **params)
    
def set_entries(name, params, **entries):
    # TODO could make sense also to change name/params
    assert all(k not in entries for k in ("name", "params"))
    Header().set_entries(name, params, **entries)
    
def exists(name, **params):
    try:
        Header().get_file_params(name, params)
    except KeyError:
        return False
    return True     
    
# core class that communicates with txt files
class Header(object):
    """access/organize data files in directory
    usage:
    Header.get_XXX(name, **params) -- retrieve simulation data
    Header.update() -- incorporate new data in folder
    
    header file:
    dict name : list paramsets
         "_flist" : list of files in folder
    paramset also contains filename"""
    
    def __init__(self):
        self.header = _load(HEADER)
    
    # TODO: could add pretty print of header content to view existing data
    def list_files(self):
        # TODO: could include global file names and sizes in list
        return self.header["_flist"]
            
    def get_file_params(self, name, params):
        "take first file compatible with params"
        if name not in self.header:
            raise KeyError("Header: Name '%s' not found." %name)
        for params0 in self.header[name]:
            if _compatible(params0, params):
                FILE = params0["FILE"]
                break
        else:
            raise KeyError("Header: No matching parameter set.")
        return FILE, params0
        
    def get_file(self, name, params):
        "take first file compatible with params"
        FILE, _ = self.get_file_params(name, params)
        return FILE
        
    def load_file(self, name, params):
        FILE = self.get_file(name, params)
        return _load(FILE)
        
    def get_fields(self, name, **params):
        fdict = self.load_file(name, params)
        return fdict["fields"]
        
    def get_field(self, name, field, **params):
        fdict = self.load_file(name, params)
        return fdict["fields"][field]
        
    # TODO: this would be slightly more elegant and much more
    # efficient if it was possible to get handle on one file
    def get_entry(self, name, entry, **params):
        fdict = self.load_file(name, params)
        return fdict[entry]
        
    def set_entries(self, name, params, **entries):
        FILE = self.get_file(name, params)
        f = _load(FILE)
        for entry, value in entries.items():
            f[entry] = value
        _save(f, FILE)
        
    def update(self):
        "update file list, merge all mutually compatible files"
        new = self._update_list()
        N = len(new)
        n = 0
        # inspect new files and create entries in header[name]
        for FILE in new:
            f = _load(FILE)
            name, params = f["name"], f["params"].copy()
            try:
                MATCH, params0 = self.get_file_params(name, params)
            except KeyError:
                MATCH = None
            if MATCH is None:
                # create new entry
                params["FILE"] = FILE
                self._add_entry(name, params)
            else:
                # merge file into existing file
                params0.update(f["params"])
                mergefile(f, MATCH)
                self._delete_file(FILE)
                n += 1
        self._write()
        if N>0:
            print ("Found %d new files, merged %d of them into "
               "existing files.") % (N, n)
        else: print "Nothing to be updated."
        
    def reread(self):
        "completely clear existing information and read again"
        self.header = dict(_flist=[])
        self.update()
        
    def remove(self, name, params):
        FILE, params = self.get_file_params(name, params)
        self._delete_file(FILE)
        self.header[name].remove(params)
        if len(self.header[name]) == 0:
            self.header.pop(name)
        self._write()
    
    def _update_list(self):
        "update filelist field in header and return list of new files"
        flist = [f for f in os.listdir(DIR) if f.endswith(SUFFIX)]
        old = self.header["_flist"]
        self.header["_flist"] = flist
        new = [f for f in flist if f not in old]
        return new
        
    def _delete_file(self, FILE):
        if FILE in self.header["_flist"]:
            self.header["_flist"].remove(FILE)
        os.remove(os.path.join(DIR, FILE))
        
    def _add_entry(self, name, params):
        if not name in self.header:
            self.header[name] = []
        self.header[name].append(params)
        
    def _write(self):
        _save(self.header, HEADER)
    
def mergefile(f, FILE):
    "merge file content f into FILE, knowing they are compatible"
    # read FILE
    f0 = _load(FILE)
    # update f0 by f. keys "params", fields" get special treatment
    f0["params"].update(f["params"])
    dontupdate = ["params"] 
        
    if "fields" in f0 and "fields" in f:
        dontupdate.append("fields")
        for F in f["fields"]:
            if F in f0["fields"]:
                f0["fields"][F].extend(f["fields"][F])
            else:
                f0["fields"][F] = f["fields"][F]        
    # other keys are simply overwritten
    f0.update({k:f[k] for k in f if not k in dontupdate})
    # write back to file
    _save(f0, FILE)
    
# helper functions    
def _compatible(a, b):
    "check whether to dicts have no conflicting fields"
    return all([a[k] == b[k] for k in a if k in b])

def _merge(a, b):
    "merge two dicts"
    c = a.copy()
    c.update(b)
    return c

def _unique_id():
    "unique ID based on cputime"
    import numpy as np
    from time import time
    return str(np.int64(time()*1e6))
    
# basic file IO with json
# saving twice in the same file means overwriting
def _save_global(data, FILE):
    with open(FILE, "w") as f:
        json.dump(data, f)    
    
def _save(data, FILE):
    _save_global(data, os.path.join(DIR, FILE))

def _load_global(FILE):
    with open(FILE, "r") as f:
        data = json.load(f)
    return data
    
def _load(FILE):
    return _load_global(os.path.join(DIR, FILE))
    
# assert directory and header exists
if not os.path.exists(DIR):
    os.makedirs(DIR)
if not HEADER in os.listdir(DIR):
    _save(dict(_flist=[]), HEADER)

"""
automatic caching decorator

use like:
@cache("name", default=default_params)
def calculate(params):
    return results
    
=> calculate(**params) gives cached results if they exist

exists should be generic
save and load could be medium generic.
several types of return values could be covered:
-) simple json-able objects
-) fields, where calculation appends to previous one, e.g. 3D force field
-) file links e.g. for continuous 2D force fields
"""
    
class Params(dict):
    "for writing params.Qmol instead of params['Qmol']"
    def __getattr__(self, key):
        return self[key]
    
class CacheBase(object):
    def __init__(self, name, default={}, overwrite=False):
        self.name = name
        self.default = default
        self.overwrite = overwrite
        
    def __call__(self, f):
        def wrapper(*args, **params):
            params = Params(self.default, **params)
    
            if self.overwrite or not self._exists(params):
                out = f(params, *args)
                self.save(out, params)
            
            return self.load(params)
        return wrapper
        
    def _exists(self, params):
        return exists(self.name, **params)
            
class cache(CacheBase):
    "default -- for cashing simple json-able output"
    
    def save(self, result, params):
        save_entries(self.name, params, result=result)
        update()
            
    def load(self, params):
        return get_entry(self.name, "result", **params)
        
"caching discrete dolfin functions"
import dolfin
    
def _save_dolfin(data, FILE):
    dolfin.File(os.path.join(DIR, FILE)) << data

def save_functions(name, params, **functions):
    # create common file prefix
    PREFIX = name + _unique_id()
    FILE = PREFIX + SUFFIX
    data = dict(name=name, params=params, prefix=PREFIX)
    keys = functions.keys()
    data["functions"] = keys
    
    # if no functions are given, only save metadata 
    if len(keys)==0:
        data["empty"] = True
        _save(data, FILE)
        return
    data["empty"] = False
    
    # save mesh (we assume functions are defined on the same mesh!)
    mesh = functions.values()[0].function_space().mesh()
    for f in functions.values():
        assert f.function_space().mesh().id() == mesh.id()
    MESHFILE = PREFIX + "_mesh.xml"
    _save_dolfin(mesh, MESHFILE)

    # save functions
    # currently only scalar/vector CG1 is supported
    ranks = []
    for fname in keys:
        f = functions[fname]
        ranks.append(f.rank())
        FFILE = PREFIX + "_" + fname + ".xml"
        _save_dolfin(f, FFILE)
    data["ranks"] = ranks
        
    # save metadata
    _save(data, FILE)
    
def _space(mesh, rank):
    if rank==0:
        return dolfin.FunctionSpace(mesh, "CG", 1)
    elif rank==1:
        return dolfin.VectorFunctionSpace(mesh, "CG", 1)
    else:
        raise Exception("Loading Functions of rank > 1 is not supported.")
            
def _load_mesh(FILE):
    return dolfin.Mesh(str(os.path.join(DIR, FILE)))
    
def _load_function(FILE, mesh, rank):
    V = _space(mesh, rank)
    return dolfin.Function(V, str(os.path.join(DIR, FILE)))

def get_functions(name, **params):
    data = load_file(name, **params)
    if data["empty"]:
        return dict(), None
    
    PREFIX = data["prefix"]
    mesh = _load_mesh(PREFIX + "_mesh.xml")
    functions = dict()
    for fname, rank in zip(data["functions"], data["ranks"]):
        FFILE = PREFIX + "_" + fname + ".xml"
        functions[fname] = _load_function(FFILE, mesh, rank)
        
    return functions, mesh
    
# print information
def show():
    h = Header().header
    h.pop("_flist")
    for key in h:
        print "\n%s" %key
        lst = h[key]
        for i, dic in enumerate(lst):
            dic.pop("FILE")
            print "%d)" %(i+1,),
            print ", ".join(["%s=%s" %x for x in dic.items()])
            
def showfields():
    h = Header().header
    h.pop("_flist")
    for key in h:
        print "\n%s" %key
        lst = h[key]
        for i, dic in enumerate(lst):
            FILE = dic.pop("FILE")
            content = _load(FILE)
            if not "fields" in content or not content["fields"]:
                continue
            n = len(content["fields"].values()[0])
            print "-) %d field values, params:" %(n,),
            print ", ".join(["%s=%s" %x for x in dic.items()])
        