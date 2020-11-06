#!/bin/bash

# --------------------------
# --- FILL IN PARAMETERS ---
# --------------------------
# values should be 1, 2, 5, 10 or 20
nodes=1
#
strong_scaling=true
weak_scaling=true
parcel_is_dry=.true.
parcels_per_image=0


function run_qsub {
		echo "qsub -l select=${nodes}:ncpus=36:mpiprocs=36,walltime=${1} \
		-v STRONG_SCALING=${strong_scaling},WEAK_SCALING=${weak_scaling},NODES=${nodes}, \
		PARCELS_PER_IMAGE=${parcels_per_image},PARCEL_IS_DRY=${parcel_is_dry}\
		cheyenne.sub"

		# -- for testing cheyenne.sub -
		STRONG_SCALING=${strong_scaling} \
									WEAK_SCALING=${weak_scaling} \
									PARCEL_IS_DRY=${parcel_is_dry} \
									PARCELS_PER_IMAGE=${parcels_per_image} \
									NODES=${nodes} \
									bash cheyenne.sub
}


if [ $nodes -eq 1 ]; then
		time=00:10:00
		run_qsub $time
fi

if [ $nodes -eq 2 ]; then
		time=00:10:00
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
