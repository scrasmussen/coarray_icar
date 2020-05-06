module convection_exchangeable_interface
  ! use, intrinsic :: iso_fortran_env, only : event_type
  use grid_interface, only : grid_t
  use iso_c_binding, only: c_int
  use convection_type_interface
  implicit none

  private
  public :: convection_exchangeable_array_t, convection_exchangeable_list_t
  public :: convection_exchangeable_t

  type convection_t
  end type convection_t

  type, extends(convection_t) :: convection_exchangeable_array_t
     type(convection_array), allocatable, public :: local(:,:,:)
  end type convection_exchangeable_array_t
  ! interface convection_exchangeable_array_t
  !    module procedure initialize_convection_array
  ! end interface convection_exchangeable_array_t

  type, extends(convection_t) :: convection_exchangeable_list_t
     type(convection_list), allocatable, public :: local(:,:,:)
  end type convection_exchangeable_list_t

  type convection_exchangeable_t
     private
     type(convection_particle), allocatable, public :: local(:,:,:)
     type(convection_particle), allocatable :: buf_south_in(:)[:]
     type(convection_particle), allocatable :: buf_north_in(:)[:]
     type(convection_particle), allocatable :: buf_west_in(:)[:]
     type(convection_particle), allocatable :: buf_east_in(:)[:]
     type(convection_particle), allocatable :: buf_southwest_in(:)[:]
     type(convection_particle), allocatable :: buf_southeast_in(:)[:]
     type(convection_particle), allocatable :: buf_northwest_in(:)[:]
     type(convection_particle), allocatable :: buf_northeast_in(:)[:]

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
     procedure, public :: const
     procedure, public :: send
     procedure, public :: load_buf
     procedure, public :: retrieve
     procedure, public :: exchange
     procedure, public :: process
     generic,   public :: initialize=>const

     procedure :: put_north
     procedure :: put_south
     procedure :: put_west
     procedure :: put_east
     procedure :: retrieve_north_buf
     procedure :: retrieve_south_buf
     procedure :: retrieve_west_buf
     procedure :: retrieve_east_buf
     ! TODO
     ! put_{northeast,northwest,southeast,southwest}
     ! retrieve_{northeast,northwest,southeast,southwest}
  end type convection_exchangeable_t

  ! type convection_object_t
  !    private
  !    type(convection_exchangeable_t), public :: convection_array
  ! end type convection_object_t

  interface
     module subroutine process(this, dt, its,ite, jts,jte, kts,kte, &
         temperature)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       real,           intent(in)    :: dt
       integer,        intent(in)    :: its,ite, jts,jte, kts,kte
       real, dimension(:,:,:), intent(in) :: temperature
     end subroutine

     module subroutine const(this, convection_type_enum, grid, input_buf_size, &
         halo_width, u_in, v_in, w_in, temperature, pressure)
       class(convection_exchangeable_t), intent(inout) :: this
       type(grid_t) :: grid
       integer, intent(in), optional :: input_buf_size
       integer, intent(in), optional :: halo_width
       integer(c_int), intent(in) :: convection_type_enum
       real, optional, intent(in) :: u_in,v_in,w_in
       real, dimension(:,:,:), intent(in) :: temperature, pressure
     end subroutine


     ! module subroutine const(this, grid, initial_value, buf_width, u, v, w)
     !   implicit none
     !   class(convection_exchangeable_t), intent(inout)  :: this
     !   type(grid_t),              intent(in)     :: grid
     !   type(convection_particle), intent(in), optional :: initial_value
     !   integer,                   intent(in), optional :: buf_width
     !   real,                      intent(in), optional :: u,v,w
     ! end subroutine

     module subroutine send(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine load_buf(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine retrieve(this, no_sync)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       logical,               intent(in),   optional :: no_sync
     end subroutine

     module subroutine exchange(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine put_north(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine put_south(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine retrieve_north_buf(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine retrieve_south_buf(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine


     module subroutine put_east(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine put_west(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine retrieve_east_buf(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine retrieve_west_buf(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

  end interface


  ! type, extends(convection_exchange_t) :: convection_exchange_list_t
  !    private
  !    type(convection_list), allocatable, public :: local(:,:,:)

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
  ! end type convection_exchange_list_t


end module convection_exchangeable_interface
