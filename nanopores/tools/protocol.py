from dictio import read_dict, write_dict
import numpy as np
import copy
import os

__all__ = ["Data"]

class Data(object):
    ''' This class acts as a wrapper between a data file and the scripts that fill it with results.
        Compatible data dictionaries must support two keys:
            -) "i" - indices
            -) "status" - current status:
                0 ... remains to be calculated
                1 ... already calculated
                uid ... currently occupied by calculating process
                (the uid identifies the process uniquely)
                
        The crucial methods are __init__, get_points and write_values.
        Attributes are data (dict containing data), filename and N (number of data points).
        Example usage:
        
        >> data = Data("sim.dat")  # to read existing data file, or
        >> # data = Data("sim.dat", N=100, overwrite=True) # to create new one
        >> subdata, uid = data.get_points(20) # get subdata for 20 free data points and unique id        
        >> subdata["force"] = ... # calculate force on the 20 points
        >> data.write_values(subdata, uid) # write the changes back to data file
    '''
    
    def __init__(self, filename, N=1, overwrite=False):
        ''' Either initialize with existing file, or create new one '''
        self.filename = filename

        if os.path.exists(filename) and not overwrite:
            self.data = read_dict(filename)
            self.N = self.data["i"].shape[0]
        else:
            self.data = {"i":np.arange(N), "status":np.zeros(N, dtype="int64")}
            self.N = N
            self.write()
            
    # be careful with using "magic" access functions, these are inefficient.  
    def __getitem__(self, I): # implements Data[I]
        #self.read() <-- questionable
        if isinstance(I,str):
            return self.data[I]
        else:
            return get_subindices(self.data, I)
        
    def __setitem__(self, I, val): # Data[key] = foo
        self.read()
        self.data[I] = val
        self.write()
            
    def read(self):
        self.data = read_dict(self.filename)
            
    def write(self):
        return write_dict(self.data, self.filename)
        
    def reset(self):
        # reset status to 0
        # this will cause every calculated value to be overwritten again
        self.data["status"][:] = 0
        self.write()
        
    def newkeys(self, keys):
        # much more efficient than self[key] = np.zeros(self.N)
        self.read()
        for key in keys:
            self.data[key] = np.zeros(self.N)
        self.write()
            
    def remaining(self):
        return np.nonzero(self.data["status"] == 0)[0]
        #return np.nonzero(np.logical_not(self.data["status"]))[0]
        
    def get_points(self, N, read=True):
        ''' get max. N data points (vertices) for calculation  '''
        if read: self.read()
        
        # get indices of N remaining vertices
        vertices = self.remaining()[:N]
    
        # set status of these to occupied by unique ID and write back
        uid = unique_id()
        self.data["status"][vertices] = uid
        self.write()
        
        # return subdata
        return (self[vertices], uid)
        
    def write_values(self, other, uid):
        ''' write back values, update status '''
        # only change status to 1 if index gets returned, else reset to 0
        # nothing can happen to data if no correct uid is known
        # RETURNS remaining points
        
        # first of course we have to re-read
        # because somebody else could have changed our data
        self.read()
        data = self.data
        
        # check for valid id
        if uid in [0,1]:
            print "Error: %s is not a valid ID, returning." %uid
            return get_remaining(data)
        
        # get indices of vertices to update by id
        vertices = np.nonzero(data["status"] == uid)[0]
        
        if vertices.size == 0:
            print "Warning: ID not found, returning."
            return get_remaining(data)
        
        # check which of the relevant vertices are present in other 
        nonpresent = np.setdiff1d(vertices, other["i"])
        present = np.intersect1d(vertices, other["i"])
            
        # if contribution is non-empty, initialize keys not in data
        if present.size > 0:
            for key in other:
                if not key in data:
                    d = type(other[key][0])
                    data[key] = np.zeros(self.N, dtype=d)
        
        # update vertices present in other
        for i in present:
            ii = np.where(other["i"] == i)[0][0]
            for key in other:
                data[key][i] = other[key][ii]
            
        # reset non-present indices to 0, present to 1
        data["status"][nonpresent] = 0
        data["status"][present] = 1
        
        self.write()
        return get_remaining(data)
        
    def set_free(self, uid=None):
        ''' reset status of occupied points to zero '''
        self.read()
        data = self.data
        
        # check for valid id
        if uid in [0,1]:
            print "Error: %s is not a valid ID, returning." %uid
            return None
            
        if uid is not None:
            # get indices of vertices by id
            vertices = np.nonzero(data["status"] == uid)[0]
        else:
            vertices = np.nonzero(np.logical_and(data["status"] != 0, data["status"] != 1))[0]
            
        # reset vertices to 0
        data["status"][vertices] = 0
        self.write()
        
def get_subindices(data, I):
    return {key: a[I] for key, a in data.items()}
        
def get_remaining(data):
    ''' number of remaining vertices for calculation '''
    return np.count_nonzero(np.logical_not(data["status"]))
        
def unique_id():
    ''' unique ID based on cputime '''
    from time import time
    return np.int64(time()*1e6)

