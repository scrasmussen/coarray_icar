module convection_type_interface
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
     logical :: exists, moved
     real :: x, y, z
     real :: u, v, w
     real :: pressure, temperature
   contains
     procedure :: move_particle
  end type convection_particle
  interface convection_particle
     module procedure constructor
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
  subroutine move_particle(from,to)
    class(convection_particle), intent(inout) :: from
    type(convection_particle), intent(inout)  :: to
    if (to%moved .eqv. .true.) then
      error stop &
          "TO Particle has been Moved, need to combine with FROM Particle"
    else if (to%exists .eqv. .true.) then
      error stop "TO Particle Exists, need to move it"
    end if
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
  end subroutine move_particle

  function constructor(x,y,z,u,v,w,pressure,temperature) result(this)
    type(convection_particle) :: this
    integer :: x, y, z
    real :: u, v, w, pressure, temperature
    this%exists = .true.
    this%moved = .false.
    this%x = x
    this%y = y
    this%z = z

    this%u = u
    this%v = v
    this%w = w
    this%pressure = pressure
    this%temperature = temperature + 30
    print *, "------IN CONSTRUCTOR: currently adding 30 K to temp"
  end function constructor
end module convection_type_interface
