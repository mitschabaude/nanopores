""" attempt at a more general-purpose parallel simulation script using the 2D solver.

should do the following: simulate forces in the pore for a given list of ranges of parameter values.
distribute this simulation to a given number of processors.
create a data and metadata file for every range.
#if data file already exists and metadata match, attempt to finish the data file.
if simulation is finished at the end, save plot of the range in same DIR.
"""

from ..tools.protocol import Data, unique_id
from ..tools.utilities import save_dict
from ..tools.mpipool import mpimap
from mpi4py import MPI
from pathos.helpers import mp # mp = fork of multiprocessing package
from ..dirnames import DATADIR
import numpy, os
#from .calculate_forces import calculate2D

__all__ = ["iterate_in_parallel", "post_iteration", "simulate",
           "parallel_output"]

# directory where data are saved
savedir = DATADIR + "/sim/stamps/"
if not os.path.exists(savedir):
        os.makedirs(savedir)

# general functions for running simulations and saving output
# TODO: this is a bit messy.. need class for storing and iterating
#       through parameter sets (while holding some fixed)
def iterate_in_parallel(method, nproc=1, iterkeys=None, **params):
    ''' evaluate a given method for a given parameter set.

        params is a dict and some of its values are allowed to be iterable.
        the method is expected to return a dict with the SAME KEYS for every parameter in one iteration.
        an exception ocurring in the method is NOT handled without stopping the iteration.
    '''
    # find the parameters to be iterated through
    iterkeys2 = [key for key in params if hasattr(params[key], "__iter__")]
    
    if iterkeys is None:
        iterkeys = iterkeys2
    elif set(iterkeys) <= set(iterkeys2):
        for key in iterkeys:
            iterkeys2.remove(key)
        iterkeys = iterkeys + iterkeys2
    else:
        print("I'm ignoring your iterkeys.")
        iterkeys = iterkeys2

    # create stamp of the input
    stamp = dict(params)
    stamp["iterkeys"] = iterkeys
    stamp["method"] = method.__name__

    # create list of params instances to be mapped
    iterator = combinations(params, iterkeys)

    # create the function to be mapped with
    def f(params): return method(**params)

    # map iterator using mpi4py
    # FIXME: doesn't work if some dolfin function are used, e.g. Function.extrapolate
    if MPI.COMM_WORLD.Get_size() > 1:
        result = mpimap(f, iterator)
    # map iterator using multiprocessing.Pool
    # FIXME: this approach of distributing across multiple processors is inconvenient
    #        since a single error kills the whole simulation.
    #        (not necessarily, error can be catched and displayed by method)
    #        also it's not supposed to be appropriate for HPC architectures
    elif nproc>1:
        pool = mp.Pool(nproc)
        result = pool.map(f, iterator)
        pool.close()
        pool.join()
    # map in serial
    else:
        result = list(map(f, iterator))

    return join_dicts(result), stamp


def combinations(dic, iterkeys):
    # Input: dict of iterables and/or single values, list of iterable keys to provide order
    # Output: list of dicts with all possible combinations of single values
    P = [{k:dic[k] for k in dic if k not in iterkeys}]
    for key in iterkeys:
        P2 = []
        for val in dic[key]:
            for p in P:
                p2 = dict(p)
                p2[key] = val
                P2.append(p2)
        P = P2
        #print P
    return P

def join_dicts(lst):
    # [{"F":1.0}, {"F":2.0}, ...] --> {"F":[1.0, 2.0, ...]}
    keys = []
    for dic in lst:
        if dic is not None:
            keys = list(dic.keys())
            
    return {key:[dic[key] for dic in lst if dic is not None] for key in keys}

def post_iteration(result, stamp, showplot=False):
    ''' in case method output is a dict, put result of iterate_in_parallel
        into nicer form, save in .dat file and create plots '''
    # create unique id for filenames
    uid = str(unique_id())

    # save stamp to file
    save_dict(stamp, dir=savedir, name=("stamp"+uid))

    # put result and input into form
    #result = join_dicts(result)
    iterkeys = stamp.pop("iterkeys")
    # create combinations only of relevant (iterable) parameters
    input_params = {k:stamp[k] for k in iterkeys} # can be empty
    input = join_dicts(combinations(input_params, iterkeys))

    # save iterated parameters and result to data file
    N = len(list(result.values())[0])
    data = Data(savedir+"result"+uid+".dat", N=N, overwrite=True)
    data.data["status"][:] = 1
    for key in input:
        data.data[key] = numpy.array(input[key])
    for key in result:
        data.data[key] = numpy.array(result[key])
    data.write()

    # no plot if not at least two different parameter sets
    if len(iterkeys) == 0:
        return

    # create plot for every result column
    # TODO: for the moment i assume that iterkeys[0] is the one to be plotted
    # thus i can use numpy indexing to get the right chunks of the results
    #if plotkey is None:
    plotkey = iterkeys[0]
    iterkeys.remove(plotkey)
    x = stamp.pop(plotkey)
    nx = len(x)
    # create combinations only of relevant (iterable) parameters
    input_params = {k:stamp[k] for k in iterkeys}
    params = combinations(input_params, iterkeys)

    from matplotlib.pyplot import plot, xlabel, ylabel, legend, savefig, show
    import matplotlib.pyplot as plt
    plots = {}

    # for every result column
    for key, rescol in list(result.items()):
        i = 0
        # create new figure
        fig, ax = plt.subplots()
        plots[key] = ax
        # for every non-axis input parameter set held fixed
        for pset in params:
            # get the corresponding chunk of length nx of result column
            chunk = slice(i*nx, (i+1)*nx)
            i += 1
            y = rescol[chunk]
            # create fitting label using the fixed params
            label = ", ".join("%s=%s" % t for t in list(pset.items()))
            # add x,y to plot and label axis with keys
            #print x,y
            plot(x, y, '-x', label=label)
            xlabel(plotkey)
            ylabel(key)
            legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        savefig(savedir+"plot"+uid+key+".eps", bbox_inches='tight')
    if showplot: show()
    return plots
    #else: close()

# general simulation (for modules with calculate() function)
# optionally, module can also provide post_calculate() function that receives the result
def simulate(name, nproc=1, outputs=None, plot=None,
             write_files=True, **params):
    script = __import__("nanopores.scripts."+name, fromlist=["calculate"])
    calculate = script.calculate

    if outputs is not None:
        def f(**x):
            res = calculate(**x)
            return {key:res[key] for key in outputs if key in res}
    else:
        f = calculate
    if plot is not None:
        result, stamp = iterate_in_parallel(f, nproc=nproc, iterkeys=[plot], **params)
    else:
        result, stamp = iterate_in_parallel(f, nproc=nproc, **params)
    if MPI.COMM_WORLD.Get_rank() > 0 or not write_files:
        return

    stamp["script"] = name
    print(result, stamp)
    if hasattr(script, "post_calculate"):
        script.post_calculate(result, stamp)
    else:
        post_iteration(result, stamp, showplot=False)
    return result
    
def parallel_output(calculate, nproc=1, plot=None, showplot=False, **params):
    if plot is not None:
        result, stamp = iterate_in_parallel(calculate, nproc=nproc, iterkeys=[plot], **params)
    else:
        result, stamp = iterate_in_parallel(calculate, nproc=nproc, **params)
    if MPI.COMM_WORLD.Get_rank() > 0:
        return
    plots = post_iteration(result, stamp, showplot=showplot)
    return plots
                 

# simulation in 2D (script for howorka pore)
#def simulation2D(nproc=1, outputs=None, plot=None, write_files=True, **params):
#    if outputs is not None:
#        def f(**x):
#            res = calculate2D(**x)
#            return {key:res[key] for key in outputs}
#    else:
#        f = calculate2D
#    if plot is not None:
#        result, stamp = iterate_in_parallel(f, nproc=nproc, iterkeys=[plot], **params)
#    else:
#        result, stamp = iterate_in_parallel(f, nproc=nproc, **params)
#    if MPI.COMM_WORLD.Get_rank() > 0 or not write_files:
#        return
#    post_iteration(result, stamp, showplot=False)
#    return result
