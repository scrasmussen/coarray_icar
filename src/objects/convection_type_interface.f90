module convection_type_interface
  implicit none
  private
  public :: convection_particle

  type convection_particle
     integer :: particle_id
     logical :: exists = .false., moved
     real :: x, y, z
     real :: u, v, w
     real :: z_meters, z_interface
     real :: pressure, temperature, potential_temp
     real :: velocity, water_vapor, cloud_water
   contains
     procedure :: move_to
  end type convection_particle

  interface convection_particle
     module procedure :: constructor ! not being used right now
  end interface convection_particle

contains
  subroutine move_to(from,to)
    class(convection_particle), intent(inout) :: from
    type(convection_particle), intent(inout)  :: to
    ! handle the from
    from%exists = .false.
    from%moved  = .false.

    ! handle the to
    to%exists = .true.
    to%moved  = .true.
    to%particle_id = from%particle_id
    to%x = from%x; to%y = from%y; to%z = from%z
    to%u = from%u; to%v = from%v; to%w = from%w
    to%z_meters = from%z_meters
    to%z_interface = from%z_interface
    to%pressure = from%pressure
    to%temperature = from%temperature
    to%potential_temp = from%potential_temp
    to%velocity = from%velocity
    to%water_vapor = from%water_vapor
    to%cloud_water = from%cloud_water
  end subroutine move_to

  ! not being used right now
  function constructor() result(this)
      type(convection_particle) :: this
  end function constructor
end module convection_type_interface
