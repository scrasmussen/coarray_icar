#!/bin/bash


# --- gnu ---
# COMPILER=gnu \
# 				make \
# 				USE_ASSERTIONS=.false.

# --- nvidia ---
COMPILER=nvidia \
LOOP=DO_LOOP \
USE_ASSERTIONS=.false. \
				make test-ideal
