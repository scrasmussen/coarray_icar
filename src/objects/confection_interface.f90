module convection_interface
  implicit none

  private
  public :: convection_t

  type convection_obj
     real :: pressure
     real :: temp
     ! add more as needed
  end type convection_obj

  type convection_node
     type(convection_obj), allocatable :: val
     type(convection_node), pointer :: next => null()
     type(convection_node), pointer :: prev => null()
  end type convection_node

  type convection_list
     private
     integer :: num_nodes
     type(convection_node), pointer :: head => null()
     type(convection_node), pointer :: tail => null()
   contains
     ! procedure :: add, remove
  end type convection_list

  type convection_t
     private
     type(convection_list), allocatable, public :: local
     type(convection_obj), allocatable :: halo_south_in[:]
     type(convection_obj), allocatable :: halo_north_in[:]
     type(convection_obj), allocatable :: halo_west_in[:]
     type(convection_obj), allocatable :: halo_east_in[:]
     type(convection_obj), allocatable :: halo_southwest_in[:]
     type(convection_obj), allocatable :: halo_southeast_in[:]
     type(convection_obj), allocatable :: halo_northwest_in[:]
     type(convection_obj), allocatable :: halo_northeast_in[:]

     logical, allocatable :: halo_south_done[:]
     logical, allocatable :: halo_north_done[:]
     logical, allocatable :: halo_west_done[:]
     logical, allocatable :: halo_east_done[:]
     logical, allocatable :: halo_southwest_done[:]
     logical, allocatable :: halo_southeast_done[:]
     logical, allocatable :: halo_northwest_done[:]
     logical, allocatable :: halo_northeast_done[:]


     logical :: north_boundary=.false.
     logical :: south_boundary=.false.
     logical :: east_boundary=.false.
     logical :: west_boundary=.false.
     logical :: northeast_boundary=.false.
     logical :: northwest_boundary=.false.
     logical :: southeast_boundary=.false.
     logical :: southwest_boundary=.false.

   contains
     private
     ! SIMILAR TO EXCHANGE_T
     ! procedure :: const, send, retrieve, exchange, initialize procedures
     ! put_{north,south,west,east,northeast,northwest,southeast,southwest}
     ! retrieve_{north,south,west,east,northeast,northwest,southeast,southwest}
     ! etc. etc.
  end type convection_t

end module convection_interface
