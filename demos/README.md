# The 3rdHand visualization slice of the MLR code base

This code is a slice of the full MLR code base containing code relative to the
3rdHand project, e.g. to visualize the data collected by MLR and Inria.  This
is the reason for the unusual directory layout.

## Installation

Generically, the following should work:

cd share
make
cd projects/mocap/<project_name>
make
./x.exe

On problems, maybe you have to install some more Ubuntu
packages or adapt to other OS. Please refer to the files

* README.mlr.md
* install/INSTALL_ALL_UBUNTU_PACKAGES.sh

These two refer to the whole MLR code base. You probably need to install only a
fraction of these packages.  The current packages support the installation for
Ubuntu 14.04; Ubuntu 12.04 (and possibly others) are supported by installing
older versions of the indicated packages.

## Documentation

Read the README file relative to each project to learn the details of each
individual program.

## Contact

Andrea Baisero
andrea.baisero@ipvs.uni-stuttgart.de

