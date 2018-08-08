#!/bin/bash

# mpi/mpich-3.3_gcc-8

FC=mpif90 cmake ../ \
  -DCMAKE_BUILD_TYPE=DEBUG  \
  -DCMAKE_Fortran_FLAGS="-O3 -g -cpp" \
  -DCMAKE_INSTALL_PREFIX=$HOME/doestnmatter
