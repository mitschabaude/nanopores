import numpy as np

raw_input('THIS SCRIPT MAY DELETE EXISTING NUMPY ARRAYS OF DATA!! Warning 1/3')
raw_input('THIS SCRIPT MAY DELETE EXISTING NUMPY ARRAYS OF DATA!! Warning 2/3')
raw_input('THIS SCRIPT MAY DELETE EXISTING NUMPY ARRAYS OF DATA!! Warning 3/3')

counter = np.array([0,0])
EXIT_X, EXIT_Y, EXIT_Z, TIME = np.array([]), np.array([]), np.array([]), np.array([])
np.save('exit_x',EXIT_X)
np.save('exit_y',EXIT_Y)
np.save('exit_z',EXIT_Z)
np.save('timer',TIME)
np.save('counter',counter)

