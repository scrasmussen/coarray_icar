module convection_exchangeable_interface
  use grid_interface, only : grid_t
  use iso_c_binding, only: c_int
  use convection_type_interface
  use exchangeable_interface, only : exchangeable_t
  implicit none

  private
  public :: convection_exchangeable_t, num_particles

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
     procedure, public :: const
     procedure, public :: send
     procedure, public :: retrieve
     procedure, public :: exchange
     procedure, public :: process
     generic,   public :: initialize=>const

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


  interface
     module subroutine process(this, nx_global, ny_global, &
         ims, ime, kms, kme, jms, jme, dt, dz, temperature, z_interface, &
         its, ite, kts, kte, jts, jte)
       implicit none
       class(convection_exchangeable_t), intent(inout) :: this
       integer, intent(in) :: nx_global, ny_global
       real, intent(in)    :: dt, dz
       integer, intent(in) :: ims, ime, kms, kme, jms, jme
       integer, intent(in) :: its, ite, kts, kte, jts, jte
       real, intent(in) :: temperature(ims:ime,kms:kme,jms:jme)
       real, intent(in) :: z_interface(ims:ime,jms:jme)
     end subroutine

     module subroutine const(this, potential_temp, u_in, v_in, w_in, grid, z_m, &
         z_interface, ims, ime, kms, kme, jms, jme, dz_value, &
         its, ite, kts, kte, jts, jte, input_buf_size, halo_width)
       class(convection_exchangeable_t), intent(inout) :: this
       class(exchangeable_t), intent(in)    :: potential_temp
       class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
       type(grid_t), intent(in)      :: grid
       real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
       real, intent(in)              :: z_interface(ims:ime,jms:jme)
       integer, intent(in)           :: ims, ime, kms, kme, jms, jme
       integer, intent(in)           :: its, ite, kts, kte, jts, jte
       real, intent(in)              :: dz_value
       integer, intent(in), optional :: input_buf_size
       integer, intent(in), optional :: halo_width
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

     module function num_particles()
       integer :: num_particles
     end function num_particles
  end interface
end module convection_exchangeable_interface
