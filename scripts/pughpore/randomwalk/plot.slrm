#!/bin/bash

#SBATCH -J events3_two
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<bstadlbauer@posteo.de>

# when srun is used, you need to set:
export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/current/lib/libpmi.so

module load intel/16 intel-mpi/5 eigen/3.2.4 intel-mkl/11 cmake/3.6.2 parmetis/4.0.3 python/2.7 swig/3.0.10 boost/1.58.0 mpi4py/2.0.0 numpy/1.9.1 petsc/3.7.5 petsc4py/3.7.0 fenics/2016.2.0
module load matplotlib/1.5.3
module load scipy/0.18.0

srun -q -l python plot.py both 400 events10_onlyone
