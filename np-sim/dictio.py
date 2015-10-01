from collections import defaultdict
from numpy import *
import os

__all__ = ["write_dict","read_dict"]

def write_dict(data, filename):
    ''' write dict {key:array} in file '''
    if not data.has_key("i"):
        raise SystemExit("Data dict has no index column!")
    
    f = open(filename,"w")
    
    # write keys
    for key in data:
        print >>f, key,
    print >>f
    
    # write array types
    for key in data:
        print >>f, type(data[key][0]),
    print >>f
    
    # write array lengths
    for key in data:
        print >>f, data[key].shape[0],
    print >>f
    
    # write actual data
    for i in data["i"]:
        for key in data:
            if data[key].shape[0] > i:
                print >>f, data[key][i],
            else:
                print >>f, "X",
        print >>f
    
    f.close()
    return filename

def read_dict(filename):
    ''' read defaultdict {key:array} from file written by write_dict '''      
    f = open(filename,"r")
    key = f.next().strip().split(" ")
    dtype = [t[1:-2].strip("numpy.") for t in f.next().strip().split(" ")[1::2]]
    length = [int(s) for s in f.next().strip().split(" ")]
    
    data = {k : zeros(l, dtype=d) for k, d, l in zip(key,dtype,length)}
    cols = range(len(key))
    index = key.index('i')
    N = length[index]
    
    for line in f:
        l = line.strip().split(" ")
        i = eval(l[index])
        for j in cols:
            if length[j]>i:
                data[key[j]][i] = eval(l[j])
    
    return data
    #return defaultdict(lambda:zeros(N), **data) #TODO: useful?
    
#if not os.path.exists(filename)

