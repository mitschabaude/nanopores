"""
this module specifies a number of often-occuring routines and enables to call them easily for different parameter ranges

a script is defined as a .py module containing a function calculate(). this function is expected to accept arbitrary keyword arguments and no positional arguments. this allows the script to be called by a generic routine with a range of different parameter sets, possibly in parallel, where parameter sets are implemented as dicts.
the script only has to take care that every relevant parameter to its execution gets recognized correctly by calculate(), which can be assured most conveniently by defining default values for all those parameters globally at the module level and do something like
    def calculate(**params):
        globals().update(params)
        ...
the output of calculate is expected to be a dict again -- preferably with string keys and float or integer values (since the values should be recognized by plot engines).

the script can then be recognized by its name NAME.py and be called conveniently in parallel for a range of parameter sets with the following simple code:

from nanopores import simulate
simulate("NAME",
param1 = value1,
...
paramN = valueN,
)

where all the values are allowed to be either single values or lists of values.
"""

from simulation2D import *

