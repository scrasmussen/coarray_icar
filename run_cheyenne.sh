#!/bin/bash
#PBS -N job_name
#PBS -A project_code
#PBS -q queue_name
#PBS -j oe
#PBS -m n
#PBS -M gutmann@ucar.edu
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32


# --- fill in parameters ---
nodes=1
strong_scaling=false
weak_scaling=true
exe=./test-ideal
# --- parameters automatically calculated ---
n=$((${nodes} * 18))


# --- define user functions ---
file=input-parameters.txt
function create_input {
    rm -f $file
    echo "&grid nx=$1,ny=$2,nz=30 /" > $file
    echo "" >> $file
    echo "&options preferred_ratio=1/" >> $file
    cat $file
}


set -e

# ---- strong scaling ----
if [[ "$strong_scaling" = true ]]; then
    if [[ ${nodes} == 1 ]]; then
	mpiexec -np 1 ${exe}
	mpiexec -np 2 ${exe}
	mpiexec -np 4 ${exe}
	mpiexec -np 8 ${exe}
	mpiexec -np 16 ${exe}
	mpiexec -np 36 ${exe}
    else
	mpiexec -np ${n} ${exe}
    fi
fi
#
# | Nodes | Num Images |
# |     1 |          1 |
# |     1 |          2 |
# |     1 |          4 |
# |     1 |          8 |
# |     1 |         16 |
# |     1 |         36 |
# |     2 |         72 |
# |     5 |        180 |
# |    10 |        360 |
# |    20 |        720 |
# ---- end strong scaling ----



# ---- weak scaling ----
if [[ "$weak_scaling" = true ]]; then
    if [[ ${nodes} == 1 ]]; then
	mpiexec -np 1 ${exe}
	mpiexec -np 36 ${exe}
    else
	mpiexec -np ${n} ${exe}
    fi
fi
#
# 12000 per image
# | Nodes | Num Images |   x |   y |  z |
# |     1 |          1 |  20 |  20 | 30 |
# |     1 |         36 | 120 | 120 | 30 |
# |     2 |         72 | 160 | 180 | 30 |
# |     5 |        180 | 250 | 288 | 30 |
# |    10 |        360 | 375 | 384 | 30 |
# |    20 |        720 | 500 | 576 | 30 |


# 768000 per image
# | Nodes | Num Images |    x |    y |  z |
# |     1 |          1 |  160 |  160 | 30 |
# |     1 |         36 |  960 |  960 | 30 |
# |     2 |         72 | 1280 | 1440 | 30 |
# |     5 |        180 | 2250 | 2048 | 30 |
# |    10 |        360 | 3000 | 3072 | 30 |
# |    20 |        720 | 4096 | 4500 | 30 |
