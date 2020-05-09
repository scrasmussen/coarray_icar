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
     type(convection_particle), allocatable, public :: local(:)
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

     integer :: north_i=1
     integer :: south_i=1
     integer :: east_i=1
     integer :: west_i=1


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
     procedure :: retrieve_buf
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

     module subroutine send(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine load_buf(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine retrieve_buf(this, buf)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: buf(:)[*]
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

     module subroutine put_south(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine put_east(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine put_west(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

  end interface
end module convection_exchangeable_interface
