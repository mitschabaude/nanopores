#!/bin/bash

sudo pip install git+https://github.com/uqfoundation/dill.git@master
sudo pip install git+https://github.com/uqfoundation/pathos.git@master

sudo apt-get install python-mpi4py
sudo apt-get install cython
sudo apt-get install python-h5py
sudo apt-get install hdf5-tools

wget https://pypi.python.org/packages/source/h/h5py/h5py-2.5.0.tar.gz
tar xzf h5py-2.5.0.tar.gz
cd h5py-2.5.0/
export CC=mpicc
python setup.py configure --mpi
python setup.py build
sudo python setup.py install
