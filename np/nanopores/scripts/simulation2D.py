""" attempt at a more general-purpose parallel simulation script using the 2D solver.

should do the following: simulate forces in the pore for a given list of ranges of parameter values.
distribute this simulation to a given number of processors.
create a data and metadata file for every range.
if data file already exists and metadata match, attempt to finish the data file.
if simulation is finished at the end, save plot of the range in same DIR.
"""

from .calculate_forces import calculate2D
from pathos.helpers import mp # mp = fork of multiprocessing package

__all__ = ["iterate_in_parallel"]

# params is expected to be a collection of parameter sets.

def iterate_in_parallel(method, nproc=1, **params):
    ''' evaluate a given method for a given parameter set and save input as well as output.
    
        params is a dict and EXACTLY ONE of its values is expected to be iterable.
        the method is expected to return a dict with the SAME KEYS for every parameter in one iteration.
        an exception ocurring in the method is handled without stopping the iteration.
    '''
    # find the one parameter to be iterated through
    iterkey = None
    for key, value in params.items():
        if hasattr(value, "__iter__"):
            if iterkey is not None:
                raise Exception("Only one of the parameters was expected to be iterable.")
            iterkey = key
    if iterkey is None:
        raise Exception("At least one of the parameters was expected to be iterable.")
    
    # create the iterator and function to be mapped with
    iterator = params.pop(iterkey)
    def f(x):
        params[iterkey] = x
        return method(**params)
        
    # map iterator using multiprocessing.Pool
    pool = mp.Pool(nproc)
    result = pool.map(f, iterator)
    pool.close()
    pool.join()
    print result
    
if __name__=="__main__":
    iterate_in_parallel(calculate2D, nproc=4, Qmol=[1.,2.,3.,4.], refinement=False)
    
