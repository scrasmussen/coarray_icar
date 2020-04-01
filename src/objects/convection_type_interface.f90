module convection_type_interface
  public :: convection_particle

  type convection_particle
     ! this will contain 5-10 scalars
     real :: pressure
     real :: temp
  end type convection_particle

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

  interface

     module subroutine test()
       implicit none
     end subroutine

  end interface
end module convection_type_interface
