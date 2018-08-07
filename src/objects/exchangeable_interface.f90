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

    integer :: send_north_request
    integer :: send_south_request
    integer :: send_west_request
    integer :: send_east_request
    integer :: recv_north_request
    integer :: recv_south_request
    integer :: recv_west_request
    integer :: recv_east_request
    integer :: send_north_type = 0
    integer :: send_south_type = 0
    integer :: send_west_type  = 0
    integer :: send_east_type  = 0
    integer :: recv_north_type = 0
    integer :: recv_south_type = 0
    integer :: recv_west_type  = 0
    integer :: recv_east_type  = 0



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
    procedure, public :: retrieve
    procedure, public :: exchange
    generic,   public :: initialize=>const
    procedure, public :: get_tag

    procedure :: exchange_north
    procedure :: exchange_south
    procedure :: exchange_west
    procedure :: exchange_east
    procedure :: send_north
    procedure :: send_south
    procedure :: send_west
    procedure :: send_east
    procedure :: recv_north
    procedure :: recv_south
    procedure :: recv_west
    procedure :: recv_east
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

    module subroutine retrieve(this, no_sync)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      logical,               intent(in),   optional :: no_sync
    end subroutine

    module subroutine exchange(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine get_tag(this, sendrecv, from, to, tag)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer, intent(in)  :: sendrecv, from, to
      integer, intent(out) :: tag
    end subroutine

    module subroutine retrieve_north_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_south_halo(this)
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

    module subroutine exchange_north(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine exchange_south(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine exchange_east(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine exchange_west(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine send_north(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine send_south(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine send_east(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine send_west(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine recv_north(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine recv_south(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine recv_east(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine recv_west(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine


  end interface

end module
