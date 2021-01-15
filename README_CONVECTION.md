# Convection Readme

## Turning on air parcels

Within convection_exchangeable_implementation.f90 are various parameters to
tune the behaviour of the air parcels.
The file parcel-parameter.txt is used to control the total number of air
parcels, the wind speed, and the saturation level.
If the variable parcel_is_dry is .true. than the parcels are created without
water vapor and if the variable is .false. the parcels are created saturated.


The file test-ideal.f90 was used to test the air parcels and would offer a
good starting point.


## Data
The timing data produced by this project is kept within the data directory.
The data directory also contains the python files used to produce the graphs.

Additionally the git repository https://github.com/scrasmussen/icar_data/ is
used to save the larger data files, such as the ones that track particle
movement over time.