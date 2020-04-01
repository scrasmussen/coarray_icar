module convection_exchange_interface
  use, intrinsic :: iso_fortran_env, only : event_type
  use convection_type_interface
  implicit none

  private
  public :: convection_exchange_array_t, convection_exchange_list_t

  type convection_exchange_array_t
     private
     type(convection_array), allocatable, public :: local(:,:,:)
     type(convection_particle), allocatable :: halo_south_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_north_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_west_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_east_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_southwest_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_southeast_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_northwest_in(:,:,:)[:]
     type(convection_particle), allocatable :: halo_northeast_in(:,:,:)[:]

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
  end type convection_array_t


  type convection_exchange_list_t
  !    private
  !    type(convection_list), allocatable, public :: local
  !    type(convection_particle), allocatable :: halo_south_in[:]
  !    type(convection_particle), allocatable :: halo_north_in[:]
  !    type(convection_particle), allocatable :: halo_west_in[:]
  !    type(convection_particle), allocatable :: halo_east_in[:]
  !    type(convection_particle), allocatable :: halo_southwest_in[:]
  !    type(convection_particle), allocatable :: halo_southeast_in[:]
  !    type(convection_particle), allocatable :: halo_northwest_in[:]
  !    type(convection_particle), allocatable :: halo_northeast_in[:]



     ! still figuring out how these will be used
     ! event to declare if there are more particles in the linked list to move
     ! type(event_type) :: halo_south_type[*]
     ! type(event_type) :: halo_north_type[*]
     ! type(event_type) :: halo_west_type[*]
     ! type(event_type) :: halo_east_type[*]
     ! type(event_type) :: halo_southwest_type[*]
     ! type(event_type) :: halo_southeast_type[*]
     ! type(event_type) :: halo_northwest_type[*]
     ! type(event_type) :: halo_northeast_type[*]

  !    logical :: north_boundary=.false.
  !    logical :: south_boundary=.false.
  !    logical :: east_boundary=.false.
  !    logical :: west_boundary=.false.
  !    logical :: northeast_boundary=.false.
  !    logical :: northwest_boundary=.false.
  !    logical :: southeast_boundary=.false.
  !    logical :: southwest_boundary=.false.

  !  contains
  !    private
  !    ! SIMILAR TO EXCHANGE_T
  !    ! procedure :: const, send, retrieve, exchange, initialize procedures
  !    ! put_{north,south,west,east,northeast,northwest,southeast,southwest}
  !    ! retrieve_{north,south,west,east,northeast,northwest,southeast,southwest}
  !    ! etc. etc.
  end type convection_exchange_list_t


end module convection_exchange_interface
