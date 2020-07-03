module convection_type_interface
  ! use domain_interface, only : pressure_at_elevation, exner_function, sat_mr
  implicit none
  private
  public :: convection_particle, convection_array
  public :: convection_node, convection_list
  public :: convection_particle_e, convection_array_e, &
      convection_linked_list_e

  enum, bind(C)
     enumerator :: convection_particle_e
     enumerator :: convection_array_e
     enumerator :: convection_linked_list_e
  end enum

  type convection_particle
     ! when adding element, remember about constructor
     integer :: particle_id
     logical :: exists = .false., moved
     real :: x, y, z
     ! real :: xd, yd, zd  ! ARTLESS these are really the domain xyz
     real :: u, v, w
     real :: z_meters
     real :: pressure, temperature, potential_temp
     real :: velocity, water_vapor, relative_humidity
   contains
     procedure :: move_to
     procedure :: send_particle
  end type convection_particle
  interface convection_particle
     module procedure :: constructor
     ! module procedure :: constructor2
  end interface convection_particle


  type convection_array
     type(convection_particle), allocatable :: array(:)
  end type convection_array

  type convection_node
     type(convection_particle), allocatable :: val
     type(convection_node), pointer :: next => null()
     type(convection_node), pointer :: prev => null()
  end type convection_node

  type convection_list
     private
     integer :: num_nodes
     type(convection_node), pointer :: head => null()
     type(convection_node), pointer :: tail => null()
  end type convection_list

  ! interface
  !     module subroutine test()
  !     end subroutine
  ! end interface
contains
  subroutine send_particle(this,new_ijk)
    class(convection_particle), intent(inout) :: this
    integer, dimension(3), intent(in)         :: new_ijk
    print *, "TRYING TO SEND ---------"
  end subroutine send_particle

  subroutine move_to(from,to)
    class(convection_particle), intent(inout) :: from
    type(convection_particle), intent(inout)  :: to
    ! Logic for array of particles
    ! if (to%moved .eqv. .true.) then
    !   error stop &
    !       "TO Particle has been Moved, need to combine with FROM Particle"
    ! else if (to%exists .eqv. .true.) then
    !   error stop "TO Particle Exists, need to move it"
    ! end if

    ! handle the from
    from%exists = .false.
    from%moved  = .false.

    ! handle the to
    to%exists = .true.
    to%moved  = .true.
    to%x = from%x; to%y = from%y; to%z = from%z
    to%u = from%u; to%v = from%v; to%w = from%w
    to%pressure = from%pressure
    to%temperature = from%temperature
  end subroutine move_to

  ! function constructor2(x, z, y, dxyz, u, v, w, z_meters, potential_temp, temp, &
  !     water_vapor, particle_id) &
  !     result(this)
  !   type(convection_particle) :: this
  !   integer :: particle_id
  !   real :: x, z, y, dxyz, u, v, w, z_meters
  !   real :: pressure, potential_temp, temp, water_vapor, exner_val
  !   this%particle_id = particle_id
  !   this%exists = .true.
  !   this%moved = .false.

  !   this%pressure = pressure
  !   this%potential_temp = potential_temp
  !   this%temperature = temp
  !   this%water_vapor = water_vapor

  ! end function constructor2

  function constructor(particle_id,x,y,z,u,v,w,pressure,temperature, &
      water_vapor, relative_humidity) &
      result(this)
    type(convection_particle) :: this
    integer :: particle_id
    real :: x, y, z
    real :: u, v, w, pressure, temperature, water_vapor, relative_humidity
    real :: exner
    this%particle_id = particle_id
    this%exists = .true.
    this%moved = .false.
    this%x = x
    this%y = y
    this%z = z

    this%u = u
    this%v = v
    this%w = w
    this%pressure = pressure
    this%temperature = temperature * 1.01
    ! print *, "------IN CONSTRUCTOR: currently adding 0 K to temp"

    this%water_vapor = water_vapor
    this%relative_humidity = relative_humidity

    ! Amontons's law
    ! this%k = pressure / temperature

    associate(po=>100000, Rd=>287.058, cp=>1003.5)
      exner = (pressure / po) ** (Rd/cp)
    end associate
    this%potential_temp = temperature / exner
    this%velocity = 0

  end function constructor
end module convection_type_interface
