#!/bin/bash

# auf false setzen wenn du noch keinen ssh-key am repository hast:
ssh_key=true

if $ssh_key; then
    repo="git@gitlab.asc.tuwien.ac.at:Clemens/nanopores.git"
else
    repo="https://gitlab.asc.tuwien.ac.at/Clemens/nanopores.git"
fi

dir=.git_nanopores_data

cd ~
if [ -d $dir ]; then
    sudo rm -R $dir
fi

mkdir $dir   
cd $dir

git clone "$repo"
cd nanopores
git checkout data
cd np
sudo python setup.py install
cd ..
np-sim init $PWD
cd np-sim
python -c "from nanopores.force import F"

if [ $? == 0 ]; then
    echo "Installation successful."
else
    echo "Something went wrong"
fi

