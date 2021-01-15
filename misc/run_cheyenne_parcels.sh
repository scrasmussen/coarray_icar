#!/bin/bash

# --------------------------
# --- FILL IN PARAMETERS ---
# --------------------------
# values should be 1, 2, 5, 10 or 20
nodes=2
#
strong_scaling=true
weak_scaling=true
parcel_is_dry=.true.
total_parcels=-1 # parcels are being set in cheyenne_parcels.sub
np=36


function run_qsub {
    echo "qsub -l select=${nodes}:ncpus=${np}:mpiprocs=${np},walltime=${1} -v STRONG_SCALING=${strong_scaling},WEAK_SCALING=${weak_scaling},NODES=${nodes},TOTAL_PARCELS=${total_parcels},PARCEL_IS_DRY=${parcel_is_dry} cheyenne_parcels.sub"
    qsub -l select=${nodes}:ncpus=${np}:mpiprocs=${np},walltime=${1} -v STRONG_SCALING=${strong_scaling},WEAK_SCALING=${weak_scaling},NODES=${nodes},TOTAL_PARCELS=${total_parcels},PARCEL_IS_DRY=${parcel_is_dry} cheyenne_parcels.sub

		# -- for testing cheyenne_parcels.sub -
    # STRONG_SCALING=${strong_scaling} \
    #     	  WEAK_SCALING=${weak_scaling} \
    #     	  PARCEL_IS_DRY=${parcel_is_dry} \
    #     	  TOTAL_PARCELS=${total_parcels} \
    #     	  NODES=${nodes} \
    #     	  bash cheyenne_parcels.sub
}


if [ $nodes -eq 1 ]; then
		time=00:10:00
		run_qsub $time
fi

if [ $nodes -eq 2 ]; then
		time=00:30:00
		run_qsub $time
fi

if [ $nodes -eq 5 ]; then
		time=00:10:00
		run_qsub $time
fi

if [ $nodes -eq 10 ]; then
		time=00:10:00
		run_qsub $time
fi

if [ $nodes -eq 20 ]; then
		time=00:10:00
		run_qsub $time
fi
