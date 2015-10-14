""" attempt at a more general-purpose parallel simulation script using the 2D solver.

should do the following: simulate forces in the pore for a given list of ranges of parameter values.
distribute this simulation to a given number of processors.
create a data and metadata file for every range.
#if data file already exists and metadata match, attempt to finish the data file.
if simulation is finished at the end, save plot of the range in same DIR.
"""

from ..tools.protocol import Data, unique_id
from ..tools.utilities import save_dict
from pathos.helpers import mp # mp = fork of multiprocessing package
from .. import DATADIR
import numpy, os
from .calculate_forces import calculate2D

__all__ = ["iterate_in_parallel", "post_iteration", "simulation2D"]

# directory where data are saved
savedir = DATADIR + "/sim/stamps/"
if not os.path.exists(savedir):
        os.makedirs(savedir)

# two quite general function for running simulations and saving output

def iterate_in_parallel(method, nproc=1, **params):
    ''' evaluate a given method for a given parameter set.
    
        params is a dict and EXACTLY ONE of its values is expected to be iterable.
        the method is expected to return a dict with the SAME KEYS for every parameter in one iteration.
        an exception ocurring in the method is NOT handled without stopping the iteration.
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
        
    # create stamp of the input
    stamp = dict(params)
    stamp["iterated"] = iterkey
    stamp["method"] = method.__name__
    print stamp
    
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
    print {key:[dic[key] for dic in result] for key in result[0]}
    return result, stamp

def post_iteration(result, stamp, showplot=False):
    ''' in case method output is a dict, put result of iterate_in_parallel
        into nicer form, save in .dat file and create plots '''
    # create unique id for filenames
    uid = str(unique_id()) 
    
    # save stamp to file
    save_dict(stamp, dir=savedir, name=("stamp"+uid))
    
    # result = [{"F":1.0}, {"F":2.0}, ...] --> {"F":[1.0, 2.0, ...]}
    result = {key:[dic[key] for dic in result] for key in result[0]}
    
    # save iterated parameter and result to data file
    iterkey = stamp["iterated"]
    parameter = stamp[iterkey]
    N = len(parameter)
    data = Data(savedir+"result"+uid+".dat", N=N, overwrite=True)
    data.data["status"][:] = 1
    data.data[iterkey] = numpy.array(list(parameter))
    for key in result:
        data.data[key] = numpy.array(list(result[key]))
    data.write()
    
    # create plot for every result column
    from matplotlib.pyplot import plot, xlabel, ylabel, figure, savefig, show
    for key in result:
        figure()
        plot(parameter, result[key], '-x')
        xlabel(iterkey)
        ylabel(key)
        savefig(savedir+"plot"+uid+key+".eps", bbox_inches='tight')
    if showplot: show()

# simulation in 2D    
def simulation2D(nproc=1, **params):
    result, stamp = iterate_in_parallel(calculate2D, nproc=nproc, **params)
    post_iteration(result, stamp, showplot=False)
    
    
    
