# Simulation of Nanopores

### Run with fenics & python2 inside docker

```sh
# once, in the project root (this folder)
sudo docker pull quay.io/fenicsproject/stable:2016.2.0
sudo docker run -ti -v $(pwd):/home/fenics/shared --name fenics quay.io/fenicsproject/stable:2016.2.0

# to restart
sudo docker start -i fenics

# set up container
sudo apt update
sudo apt install gmsh # and paste more recent gmsh bin in this folder

# in the container, to execute scripts
cd shared
export PYTHONPATH=.

python <YOUR SCRIPT>

# since the project root is shared, you can iterate on your scripts
# from outside the container while running them from inside


# sometimes you'll need a second shell in the running container:
sudo docker exec -it -u fenics fenics /bin/bash
```
