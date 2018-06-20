module exchangeable_interface
  use grid_interface, only : grid_t
  use mpi, only : MPI_REQUEST_NULL
  implicit none

  type direction_t
    integer :: north = 1, south = 2, east = 3 , west = 4
  end type
  type sendrecv_t
    integer :: send = 0, recv = 1
  end type

  private
  public :: exchangeable_t

  type exchangeable_t
    private
    real, allocatable, public :: local(:,:,:)
    real, allocatable :: halo_south_in(:,:,:)
    real, allocatable :: halo_north_in(:,:,:)
    real, allocatable :: halo_west_in(:,:,:)
    real, allocatable :: halo_east_in(:,:,:)

    integer :: north_request
    integer :: south_request
    integer :: west_request
    integer :: east_request
    ! integer :: north_request = MPI_REQUEST_NULL
    ! integer :: south_request = MPI_REQUEST_NULL
    ! integer :: west_request = MPI_REQUEST_NULL
    ! integer :: east_request = MPI_REQUEST_NULL

    integer :: rank
    integer :: north_tag
    integer :: south_tag
    integer :: west_tag
    integer :: east_tag

    logical :: north_boundary=.false.
    logical :: south_boundary=.false.
    logical :: east_boundary=.false.
    logical :: west_boundary=.false.

  contains
    private
    procedure, public :: const
    procedure, public :: send
    procedure, public :: retrieve
    procedure, public :: exchange
    generic,   public :: initialize=>const
    procedure, public :: save_request
    procedure, public :: get_tag

    procedure :: put_north
    procedure :: put_south
    procedure :: put_west
    procedure :: put_east
    procedure :: retrieve_north_halo
    procedure :: retrieve_south_halo
    procedure :: retrieve_west_halo
    procedure :: retrieve_east_halo
  end type

  integer, parameter :: space_dim=3

  interface

    module subroutine const(this, grid, initial_value, halo_width)
      implicit none
      class(exchangeable_t), intent(inout)  :: this
      type(grid_t),          intent(in)     :: grid
      real,                  intent(in)     :: initial_value
      integer,               intent(in), optional :: halo_width
    end subroutine

    module subroutine send(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve(this, no_sync)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      logical,               intent(in),   optional :: no_sync
    end subroutine

    module subroutine exchange(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine save_request(this, request, direction)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer, intent(in) :: request
      integer, intent(in) :: direction
    end subroutine

    module subroutine get_tag(this, sendrecv, from, to, tag)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer, intent(in)  :: sendrecv, from, to
      integer, intent(out) :: tag
    end subroutine

    module subroutine put_north(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine put_south(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_north_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_south_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine


    module subroutine put_east(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine put_west(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_east_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_west_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine


  end interface

end module
