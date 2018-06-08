#define ARTLESS 1
submodule(exchangeable_interface) exchangeable_implementation
  use mpi_f08, only : MPI_COMM_WORLD, MPI_Comm_rank, MPI_Barrier, &
                      MPI_Real, MPI_STATUSES_IGNORE, MPI_Request
  use assertions_interface, only : assert, assertions
  use grid_interface, only : grid_t
  implicit none

  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_neighbor, south_neighbor, halo_size
  integer, save :: east_neighbor, west_neighbor

contains

  module subroutine const(this,grid,initial_value,halo_width)
    class(exchangeable_t), intent(inout) :: this
    type(grid_t),          intent(in)    :: grid
    real,                  intent(in)    :: initial_value
    integer,               intent(in), optional :: halo_width

    integer :: n_neighbors, current, rank, ierr
    integer :: ims,ime,kms,kme,jms,jme

    if (present(halo_width)) then
        halo_size = halo_width
    else
        halo_size = default_halo_size
    end if

    if (allocated(this%local)) deallocate(this%local)
    this%north_boundary = (grid%yimg == grid%yimages)
    this%south_boundary = (grid%yimg == 1)
    this%east_boundary  = (grid%ximg == grid%ximages)
    this%west_boundary  = (grid%ximg == 1)


    associate( halo_south => merge(0,halo_size,this%south_boundary), &
               halo_north => merge(0,halo_size,this%north_boundary), &
               halo_east  => merge(0,halo_size,this%east_boundary), &
               halo_west  => merge(0,halo_size,this%west_boundary))
      ims = grid%ims - halo_east
      ime = grid%ime + halo_west
      jms = grid%jms - halo_south
      jme = grid%jme + halo_north
      kms = grid%kms
      kme = grid%kme

      allocate(this%local(ims:ime,kms:kme,jms:jme), source=initial_value)
    end associate

    allocate( this%halo_south_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    ), source=initial_value)
    allocate( this%halo_north_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    ), source=initial_value)
    allocate( this%halo_east_in(    halo_size,     grid%halo_nz, grid%ew_halo_ny+halo_size*2), source=initial_value)
    allocate( this%halo_west_in(    halo_size,     grid%halo_nz, grid%ew_halo_ny+halo_size*2), source=initial_value)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! set up the neighbors array so we can sync with our neighbors when needed
    if (.not.allocated(neighbors)) then
      associate(me=>rank)
        south_neighbor = me - grid%ximages
        north_neighbor = me + grid%ximages
        east_neighbor  = me + 1
        west_neighbor  = me - 1

        n_neighbors = merge(0,1,this%south_boundary)  &
                     +merge(0,1,this%north_boundary)  &
                     +merge(0,1,this%east_boundary)   &
                     +merge(0,1,this%west_boundary)
        n_neighbors = max(1, n_neighbors)

        allocate(neighbors(n_neighbors))

        current = 1
        if (.not. this%south_boundary) then
            neighbors(current) = south_neighbor
            current = current+1
        endif
        if (.not. this%north_boundary) then
            neighbors(current) = north_neighbor
            current = current+1
        endif
        if (.not. this%east_boundary) then
            neighbors(current) = east_neighbor
            current = current+1
        endif
        if (.not. this%west_boundary) then
            neighbors(current) = west_neighbor
            current = current+1
        endif
        ! if current = 1 then all of the boundaries were set, just store ourself as our "neighbor"
        if (current == 1) then
            neighbors(current) = me
        endif

      end associate
    endif

  end subroutine

  module subroutine send(this)
    class(exchangeable_t), intent(inout) :: this
    ! each put has to be called by both? or just the send
    ! have get call the receive?? Doesnt make sense for isends
    ! do both immedielty, let implementation complete when able
    if (.not. this%north_boundary) call this%put_north  !! ARTLESS
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west
  end subroutine

  module subroutine retrieve(this, no_sync)
    class(exchangeable_t), intent(inout) :: this
    logical,               intent(in),   optional :: no_sync
    integer :: ierr

    if (.not. present(no_sync)) then
        call this%waitall
    else
        if (.not. no_sync) then
            call this%waitall
        endif
    endif

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary) call this%retrieve_east_halo
    if (.not. this%west_boundary) call this%retrieve_west_halo
  end subroutine

  module subroutine exchange(this)
    class(exchangeable_t), intent(inout) :: this
    if (.not. this%north_boundary) call this%put_north
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west

    call this%waitall

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary)  call this%retrieve_east_halo
    if (.not. this%west_boundary)  call this%retrieve_west_halo
  end subroutine

  module subroutine waitall(this)
    class(exchangeable_t), intent(inout) :: this
    integer :: ierr
    call MPI_Waitall(this%num_request, this%request, &
                             MPI_STATUSES_IGNORE, ierr)
    this%num_request = 0
  end subroutine

  module subroutine save_request(this, request)
    implicit none
    class(exchangeable_t), intent(inout) :: this
    type(MPI_Request), intent(in) :: request
    this%num_request = this%num_request + 1
    this%request(this%num_request) = request
  end subroutine

  module subroutine put_north(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: n, nx, len, ierr
      n = ubound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        write (*,*) "implementation does not support assert because it accesses [north_neighbor] shape"
        stop
        ! call assert( shape(this%halo_south_in(:nx,:,1:halo_size)[north_neighbor]) &
        !              == shape(this%local(:,:,n-halo_size+1:n)),         &
        !              "put_north: conformable halo_south_in and local " )
      end if


#ifdef ARTLESS
      if (size(this%local(:,:,n-halo_size*2+1:n-halo_size)) .NE. &
          size(this%halo_south_in(:nx,:,1:halo_size))) then
        print*, "ERROR: the number of elements in send and recv arrays differs"
      end if
#endif

      !dir$ pgas defer_sync
      ! this%halo_south_in(1:nx,:,1:halo_size)[north_neighbor] = this%local(:,:,n-halo_size*2+1:n-halo_size)
      len = size(this%local(:,:,n-halo_size*2+1:n-halo_size))
      if (.not. this%north_boundary) then ! send to north_neighbor
        call MPI_Isend(this%local(:,:,n-halo_size*2+1:n-halo_size), &
                       len, MPI_Real, north_neighbor, -1, &
                       MPI_COMM_WORLD, request, ierr)
        call this%save_request(request)
      end if
      if (.not. this%south_boundary) then ! get from south_neighbor
        call MPI_Irecv(this%halo_south_in(:nx,:,1:halo_size), len, &
                       MPI_Real, south_neighbor, -1, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request)
      end if

  end subroutine

  module subroutine put_south(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: start, nx, len, ierr
      start = lbound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        write (*,*) "implementation does not support assert because it accesses [south_neighbor] shape"
        stop
        ! call assert( shape(this%halo_north_in(:nx,:,1:halo_size)[south_neighbor]) &
        !              == shape(this%local(:,:,start:start+halo_size-1)), &
        !              "put_south: conformable halo_north_in and local " )
      end if
      !dir$ pgas defer_sync
      ! this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
      len = size(this%halo_north_in(1:nx,:,1:halo_size))
      if (.not. this%south_boundary) then ! send to south_neighbor
        call MPI_Isend(this%local(:,:,start+halo_size:start+halo_size*2-1), &
                       len, MPI_Real, south_neighbor, -1, &
                       MPI_COMM_WORLD, request, ierr)
        call this%save_request(request)
      end if
      if (.not. this%north_boundary) then ! get from north_neighbor
        call MPI_Irecv(this%halo_north_in(1:nx,:,1:halo_size), len, &
                       MPI_Real, north_neighbor, -1, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request)
      end if

  end subroutine

  module subroutine retrieve_north_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: n, nx
      n = ubound(this%local,3)
      nx = size(this%local,1)

      this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine retrieve_south_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: start, nx
      start = lbound(this%local,3)
      nx = size(this%local,1)

      this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine put_east(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: n, ny, len, ierr
      n = ubound(this%local,1)
      ny = size(this%local,3)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        write (*,*) "implementation does not support assert because it accesses [east_neighbor] shape"
        stop
        ! call assert( shape(this%halo_west_in(1:halo_size,:,:ny)[east_neighbor])       &
        !              == shape(this%local(n-halo_size*2+1:n-halo_size,:,:)), &
        !              "put_east: conformable halo_west_in and local " )
      end if

      !dir$ pgas defer_sync
      ! this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)

      len = size(this%halo_west_in(1:halo_size,:,1:ny))
      if (.not. this%east_boundary) then ! send to east_neighbor
        call MPI_Isend(this%local(n-halo_size*2+1:n-halo_size,:,:), &
                       len, MPI_Real, east_neighbor, -1, &
                       MPI_COMM_WORLD, request, ierr)
        call this%save_request(request)
      end if
      if (.not. this%west_boundary) then ! get from west_neighbor
        call MPI_Irecv(this%halo_west_in(1:halo_size,:,1:ny), len, &
                       MPI_Real, west_neighbor, -1, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request)
      end if
  end subroutine

  module subroutine put_west(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: start, ny, len, ierr
      start = lbound(this%local,1)
      ny = size(this%local,3)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        write (*,*) "implementation does not support assert because it accesses [west_neighbor] shape"
        stop

        ! call assert( shape(this%halo_east_in(1:halo_size,:,:ny)[west_neighbor])               &
        !              == shape(this%local(start+halo_size:start+halo_size*2-1,:,:)), &
        !              "put_west: conformable halo_east_in and local " )
      end if

      !dir$ pgas defer_sync
      ! this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)

      len = size(this%halo_east_in(1:halo_size,:,1:ny))
      if (.not. this%west_boundary) then ! send to west_neighbor
        call MPI_Isend(this%local(start+halo_size:start+halo_size*2-1,:,:), &
                       MPI_Real, west_neighbor, -1, &
                       MPI_COMM_WORLD, request, ierr)
        call this%save_request(request)
      end if
      if (.not. this%east_boundary) then ! get from east_neighbor
        call MPI_Irecv(this%halo_east_in(1:halo_size,:,1:ny), len, &
                       MPI_Real, east_neighbor, -1, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request)
      end if
  end subroutine

  module subroutine retrieve_east_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: n, ny
      n = ubound(this%local,1)
      ny = size(this%local,3)

      this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_halo(this)
      class(exchangeable_t), intent(inout) :: this

      integer :: start, ny
      start = lbound(this%local,1)
      ny = size(this%local,3)

      this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
  end subroutine

end submodule
