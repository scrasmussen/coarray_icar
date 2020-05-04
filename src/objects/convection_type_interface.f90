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
     ! this will contain 5-10 scalars, remember to edit constructor
     logical :: exist
     real :: i, j, k
     real :: u, v, w
     real :: pressure, temperature
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
  function constructor(u,v,w,pressure,temp) result(this)
    type(convection_particle) :: this
    real :: u, v, w, pressure, temp
    this%u = u
    this%v = v
    this%w = w
    print *, "------IN CONSTRUCTOR"
    ! print *, "k = ", this%k
  end function constructor
  ! function constructor(build,i,j,k,u,v,w) result(this)
  !   type(convection_particle) :: this
  !   logical :: build
  !   real :: i, j, k
  !   real :: u, v, w
  !   this%i = i
  !   this%j = j
  !   this%k = k
  !   this%u = u
  !   this%v = v
  !   this%w = w
  !   ! print *, "------IN CONSTRUCTOR"
  !   ! print *, "k = ", this%k
  ! end function constructor
end module convection_type_interface
