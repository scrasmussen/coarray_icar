module convection_type_interface
  implicit none
  private
  public :: convection_particle, convection_array
  public :: convection_node, convection_list

  type convection_particle
     ! this will contain 5-10 scalars, remember to edit constructor
     real :: pressure
     real :: temp
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
  function constructor(build) result(this)
    type(convection_particle) :: this
    logical :: build
    this%pressure = 0.0
    this%temp = 0.0
  end function constructor
end module convection_type_interface
