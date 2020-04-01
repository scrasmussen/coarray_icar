# Confection Readme

convected parcel, computationally particle         if multiple poitns its an object
 particle fine

adapt MP driver or implement new subroutine that just calls MP for particle,
here is one grid cell, calc the MP physics for that single grid cell

my convected particles in their own system,
  routines calculates how it moves
  pototenial temp of parcel, should have its own stuff, needs to do its own calculations

will there

MP expects to performs on a column of air, rain falls and if it falls through the bottom it treats it as air
  if water particle falls through partcle need to add to MP

need to spend more time looking at MP to think how it handle the MP party



to the first order, dont even need the complex MP, only one subroutine to get the saturation of air, if saturated create cloud particles and warm air
    pressure caluation for cahne in elevtion, calc of bouncy tendancy, the force, velocity calc,

atmospheric glossary, tricky theres a little bit of iteration



need to  add movement of particles within local system

## program flow : setup

Use domain module to setup and call `halo_{send,retrieve}`
                  has and initializes water_vapor, rain and snow mass, etc, etc.
 BUT
Add optional arguments to indicate particle array or linked list, and array size



Use microphysics to with domain argument to call `mp_gt_driver`, which does all the fancy stuff
 BUT
Add to/modify the microphysics driver (`mp_driver.f90` file) to handle the calculations


test-program sets up domain
 - `domain%initialize_from_file`

domain_interface: has all the variable, initializes them, calls halo send, recv, exchange
 - `domain_t`
   - `type(exchangeable_t) :: water_vapor`

exchangeable_interface: has all the halo regions, defines put and retrieve subroutines


## program flow : new stuff

=going to add=
-------------
-handles creation, destruction or particles, constantly knows its surroundings
`convection_particle_interface.f90`
`convection_particle_implementation.f90`

-handles the exchange of the halo regions for the array and linked list types, all the puts and recieves
`convection_exchangeable_interface.f90`
`convection_exchangeable_implementation.f90`

-has the different types, convection particle, array, node and linked list, etc
-defines the linked list functions
`convection_type_interface.f90`
`convection_type_implementation.f90`

## given goals
seperate interface that handles creation, destruction of convection particles, constantly knows its surroundings
switch between array and linked list at runtime
