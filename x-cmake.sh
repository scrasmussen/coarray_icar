#!/bin/bash

# Step 1. Load needed modules
# Step 2. Make sure below script matches user's need
# Step 3. Run `./x-cmake.sh` from command line

# ----- GNU -----
FC=caf cmake ../ \
  -DCMAKE_Fortran_FLAGS="-g -O3 -Wall -cpp -DUSE_ASSERTIONS=.false." \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON


# ----- CRAY SYSTEM -----
# FC=caf cmake ../ \
# 	-DCMAKE_Fortran_COMPILER_ID="Cray" \
# 	-DCMAKE_Fortran_COMPILER=ftn \
# 	-DCMAKE_Fortran_FLAGS="-O3 -DUSE_ASSERTIONS=.false. -e F -h caf" \
# 	-DCMAKE_BUILD_TYPE=Release \
# 	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
