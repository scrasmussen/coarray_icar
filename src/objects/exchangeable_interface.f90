module exchangeable_interface
  use mpi_f08, only : MPI_Request
  use grid_interface, only : grid_t
  implicit none

  private
  public :: exchangeable_t

  type exchangeable_t
    private
    real, allocatable, public :: local(:,:,:)
    real, allocatable :: halo_south_in(:,:,:)
    real, allocatable :: halo_north_in(:,:,:)
    real, allocatable :: halo_west_in(:,:,:)
    real, allocatable :: halo_east_in(:,:,:)

    type(MPI_Request) :: request(4)
    integer :: num_request=0

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
