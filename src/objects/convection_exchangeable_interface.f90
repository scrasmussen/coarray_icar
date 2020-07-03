module convection_exchangeable_interface
  ! use, intrinsic :: iso_fortran_env, only : event_type
  use grid_interface, only : grid_t
  use iso_c_binding, only: c_int
  use convection_type_interface
  use exchangeable_interface, only : exchangeable_t
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
     type(convection_particle), allocatable :: buf_north_in(:)[:]
     type(convection_particle), allocatable :: buf_south_in(:)[:]
     type(convection_particle), allocatable :: buf_east_in(:)[:]
     type(convection_particle), allocatable :: buf_west_in(:)[:]
     type(convection_particle), allocatable :: buf_northeast_in(:)[:]
     type(convection_particle), allocatable :: buf_northwest_in(:)[:]
     type(convection_particle), allocatable :: buf_southeast_in(:)[:]
     type(convection_particle), allocatable :: buf_southwest_in(:)[:]

     logical :: north_boundary=.false.
     logical :: south_boundary=.false.
     logical :: east_boundary=.false.
     logical :: west_boundary=.false.
     logical :: northeast_boundary=.false.
     logical :: northwest_boundary=.false.
     logical :: southeast_boundary=.false.
     logical :: southwest_boundary=.false.
     logical :: wrapped_north=.false.
     logical :: wrapped_south=.false.
     logical :: wrapped_east=.false.
     logical :: wrapped_west=.false.


     integer :: north_i=1
     integer :: south_i=1
     integer :: east_i=1
     integer :: west_i=1
     integer :: northeast_i=1
     integer :: northwest_i=1
     integer :: southeast_i=1
     integer :: southwest_i=1

     integer :: particle_id_count=-1

   contains
     private
     procedure, public :: const, const2
     procedure, public :: send
     procedure, public :: load_buf
     procedure, public :: retrieve
     procedure, public :: exchange
     procedure, public :: process
     generic,   public :: initialize=>const2

     procedure :: put_north
     procedure :: put_south
     procedure :: put_east
     procedure :: put_west
     procedure :: put_northeast
     procedure :: put_northwest
     procedure :: put_southeast
     procedure :: put_southwest
     procedure :: retrieve_buf
     procedure :: create_particle_id
     procedure :: setup_neighbors
  end type convection_exchangeable_t

  ! type convection_object_t
  !    private
  !    type(convection_exchangeable_t), public :: convection_array
  ! end type convection_object_t

  interface
     module subroutine process(this, nx_global, ny_global, grid, &
         dt, dz, temperature)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       integer, intent(in) :: nx_global, ny_global
       real, intent(in)    :: dt, dz
       type(grid_t), intent(in) :: grid
       real, dimension(:,:,:), intent(in) :: temperature
     end subroutine

  module subroutine const2(this, potential_temp, u_in, v_in, w_in, grid, z_m, &
      ims, ime, kms, kme, jms, jme, dz_value, &
      input_buf_size, halo_width)
    class(convection_exchangeable_t), intent(inout) :: this
    class(exchangeable_t), intent(in)    :: potential_temp
    class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
    type(grid_t), intent(in)      :: grid
    real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
    integer, intent(in)           :: ims, ime, kms, kme, jms, jme
    real, intent(in)              :: dz_value
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width
  end subroutine


     module subroutine const(this, convection_type_enum, grid, tims,time,tkms,tkme,&
         tjms,tjme, input_buf_size, &
         halo_width, u_in, v_in, w_in, temperature, pressure, water_vapor)
       class(convection_exchangeable_t), intent(inout) :: this
       type(grid_t) :: grid
       integer, intent(in) :: tims,time,tkms,tkme,tjms,tjme
       integer, intent(in), optional :: input_buf_size
       integer, intent(in), optional :: halo_width
       integer(c_int), intent(in) :: convection_type_enum
       real, optional, intent(in) :: u_in,v_in,w_in
       real, dimension(:,:,:), intent(in) :: temperature, pressure, water_vapor
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

     module subroutine put_northeast(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine put_northwest(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine put_southeast(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine put_southwest(this, particle)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(convection_particle), intent(inout) :: particle
     end subroutine

     module subroutine create_particle_id(this)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
     end subroutine

     module subroutine setup_neighbors(this,grid)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       type(grid_t), intent(in) :: grid
     end subroutine
  end interface
end module convection_exchangeable_interface
