#!/bin/bash
#PBS -N job_name
#PBS -A project_code
#PBS -q queue_name
#PBS -j oe
#PBS -m n
#PBS -M gutmann@ucar.edu
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32

set -e  # stop running on error

# --- fill in parameters ---
nodes=2
strong_scaling=true
weak_scaling=true
exe=./test-ideal
# --- parameters automatically calculated ---
n=$((${nodes} * 18))
echo "--- VARIABLES CHOSEN ---"
echo " nodes=${nodes}"
echo " processes=${n}"
echo " strong scaling=${strong_scaling}"
echo " weak scaling=${weak_scaling}"
echo " executable=${exe}"
echo "------------------------"



# --- define user functions ---
file=input-parameters.txt
function create_input {
    rm -f ${file}
    echo "&grid nx=$1,ny=$2,nz=30 /" > ${file}
    echo "" >> ${file}
    echo "&options preferred_ratio=1/" >> ${file}
    echo "created ${1}x${2}x30 parameter file"
}

function run_n_images {
    echo "mpiexec -np ${1} ${exe}"
    mpiexec -np ${1} ${exe}
}



# ---- run strong scaling ----
if [[ "$strong_scaling" = true ]]; then
    echo "- running strong scaling tests on ${nodes} nodes"
    if [[ ${nodes} == 1 ]]; then
	create_input 500 500
	run_n_images 1
	run_n_images 2
	run_n_images 4
	run_n_images 8
	run_n_images 16
	run_n_images 36

	create_input 2000 2000
	run_n_images 1
	run_n_images 2
	run_n_images 4
	run_n_images 8
	run_n_images 16
	run_n_images 36

    else
	create_input 500 500
	run_n_images ${n}
	create_input 2000 2000
	run_n_images ${n}
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


# ---- run weak scaling ----
if [[ "$weak_scaling" = true ]]; then
    echo "- running weak scaling tests on ${nodes} nodes"
    if [[ ${nodes} == 1 ]]; then
	# 12k per image
	create_input 20 20
	run_n_images 1
	create_input 120 120
	run_n_images 36
	# 768k per image
	create_input 160 160
	run_n_images 1
	create_input 960 960
	run_n_images 36
    elif [[ ${nodes} == 2 ]] ; then
	# 12000 per image
	create_input 160 180
	run_n_images ${n}
	# 768k per image
	create_input 1280 1440
	run_n_images ${n}
    elif [[ ${nodes} == 5 ]] ; then
	# 12000 per image
	create_input 250 288
	run_n_images ${n}
	# 768k per image
	create_input 2250 2048
	run_n_images ${n}
    elif [[ ${nodes} == 10 ]] ; then
	# 12000 per image
	create_input 375 384
	run_n_images ${n}
	# 768k per image
	create_input 3000 3072
	run_n_images ${n}
    elif [[ ${nodes} == 20 ]] ; then
	# 12000 per image
	create_input 500 576
	run_n_images ${n}
	# 768k per image
	create_input 4096 4500
	run_n_images ${n}
    fi
fi

# Problem size table for weak scaling
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
