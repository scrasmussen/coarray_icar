#!/bin/bash

set -e -x

# --- gnu ---
# COMPILER=gnu \
# 				make \
# 				USE_ASSERTIONS=.false.

# loop options are DO, DC, OMP
loop=OMP

# --- nvidia ---
COMPILER=nvidia \
LOOP=${loop}_LOOP \
USE_ASSERTIONS=.false. \
DEBUGVAL=.false. \
				make test-ideal
cp test-ideal run/test-ideal-${loop}
