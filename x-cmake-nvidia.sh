#!/bin/bash

# Step 1. Load needed modules
# Step 2. Make sure below script matches user's need
# Step 3. Create build directory and `cd` to it
# Step 4. Run `cp ../x-cmake.sh ./` from command line
# Step 4. Run `./x-cmake.sh` from command line

# ----- GNU -----
# FC=caf cmake ../ \
#   -DCMAKE_Fortran_FLAGS="-g -O3 -Wall -cpp -DUSE_ASSERTIONS=.false." \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON

loop="DO_LOOP"
if [ "$#" -ge 0 ]; then
    loop="$1_LOOP"
fi
loop="DO_LOOP"

set -e -x
# ----- NVidia SYSTEM -----
FC=nvfortran cmake ../ \
        -DCMAKE_Fortran_COMPILER_ID="NVidia" \
        -DCMAKE_Fortran_COMPILER=nvfortran \
        -DCMAKE_Fortran_FLAGS="\
-g -O3 \
-D${loop}=true \
-cpp" \
        -DCMAKE_BUILD_TYPE=Release \
				-DNO_ASSERTIONS=1 \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON

# redefinition, is this a bug?
# -DUSE_ASSERTIONS=.false. \


# nvidia has no coarrays



# -DDO_LOOP=true \
# -DOMP_LOOP=true -homp \
# -DDC_LOOP=true -hthread_do_concurrent \
