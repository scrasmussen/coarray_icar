!! #define ARTLESS 1
!! #define ARTLESSSEND 1
submodule(exchangeable_interface) exchangeable_implementation
  use mpi_f08, only : MPI_Comm_rank, MPI_Request, MPI_Barrier, &
                      MPI_STATUSES_IGNORE, MPI_COMM_WORLD, MPI_Status, &
                      MPI_Type_size, MPI_Double_precision
  use assertions_interface, only : assert, assertions
  use grid_interface, only : grid_t
  implicit none

  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_neighbor, south_neighbor, halo_size
  integer, save :: east_neighbor, west_neighbor
  integer, save :: north_tag, south_tag, east_tag, west_tag
  integer :: north = 1, south = 2, east = 3, west = 4

contains

  module subroutine const(this,grid,initial_value,halo_width)
    class(exchangeable_t), intent(inout) :: this
    type(grid_t),          intent(in)    :: grid
    real,                  intent(in)    :: initial_value
    integer,               intent(in), optional :: halo_width

    integer :: n_neighbors, current, ierr
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

    call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, ierr)
    this%rank = this%rank + 1
    ! set up the neighbors array so we can sync with our neighbors when needed
    if (.not.allocated(neighbors)) then
      associate(me=>this%rank)
        south_neighbor = me - grid%ximages
        north_neighbor = me + grid%ximages
        east_neighbor  = me + 1
        west_neighbor  = me - 1
        south_tag = me * 10 + south_neighbor
        east_tag  = me * 10 + east_neighbor
        north_tag = north_neighbor * 10 + me
        west_tag  = west_neighbor * 10 + me
        ! print *, this%rank, ":north=",north_neighbor, "south=",south_neighbor,"east=",east_neighbor,"west=",west_neighbor
        ! print *, this%rank, ":north_tag=",north_tag, "south_tag=",south_tag,"east_tag=",east_tag,"west_tag=",west_tag


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
    integer :: ierr
    ! have get call the receive?? Doesnt make sense for isends
    ! do both immedielty, let implementation complete when able

#ifdef ARTLESSSEND
    print *, "@@@put_north"
#endif
    call this%put_north  !! ARTLESS
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr) !! ARTLESS
! #ifdef ARTLESSSEND
!     print *, "@@@put_south"
! #endif
    call this%put_south

#ifdef ARTLESSEASTWEST
    call this%put_east
    call this%put_west
#endif
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr) !! ARTLESS
! #ifdef ARTLESSSEND
!     print *, "      === FIN SEND ===  "
! #endif
    ! if (.not. this%north_boundary) call this%put_north
    ! if (.not. this%south_boundary) call this%put_south
    ! if (.not. this%east_boundary)  call this%put_east
    ! if (.not. this%west_boundary)  call this%put_west
  end subroutine

  module subroutine retrieve(this, no_sync)
    class(exchangeable_t), intent(inout) :: this
    logical,               intent(in),   optional :: no_sync
    integer :: ierr

! #ifdef ARTLESS
    print *, "======================================="
! #endif
    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    print *, "==========COMM COMPLETE=========="
#ifdef ARTLESSEASTWEST
    if (.not. this%east_boundary) call this%retrieve_east_halo
    if (.not. this%west_boundary) call this%retrieve_west_halo
#endif
    if (.not. present(no_sync)) then
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    else
        if (.not. no_sync) then
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        endif
    endif

  end subroutine

  module subroutine exchange(this) ! ARTLESS: TODO CHANGE THIS
    class(exchangeable_t), intent(inout) :: this
    integer :: ierr
    if (.not. this%north_boundary) call this%put_north
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west

    call MPI_Barrier(MPI_Comm_world,ierr)

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary)  call this%retrieve_east_halo
    if (.not. this%west_boundary)  call this%retrieve_west_halo
  end subroutine

  ! module subroutine waitall(this)
  !   class(exchangeable_t), intent(inout) :: this
  !   integer :: ierr
  !   call MPI_Waitall(this%num_request, this%request, &
  !                            MPI_STATUSES_IGNORE, ierr)
  !   this%num_request = 0
  ! end subroutine

  module subroutine save_request(this, request, direction)
    implicit none
    class(exchangeable_t), intent(inout) :: this
    type(MPI_Request), intent(in) :: request
    integer, intent(in) :: direction
    type(direction_t) :: d
    if (direction == d%north) then
    this%north_request = request
    else if (direction == d%south) then
    this%south_request = request
    else if (direction == d%west) then
    this%west_request = request
    else if (direction == d%east) then
    this%east_request = request
    end if
  end subroutine

  module subroutine put_north(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: n, nx, len, ierr
      integer :: mpir
      real :: r
      type(MPI_Status) :: status

      ! ARTLESS
      print*, this%rank, "BOUNDARIES NESW:", this%north_boundary,this%east_boundary,this%south_boundary,this%west_boundary


      n = ubound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
#ifdef VERBOSE
        write (*,*) "implementation does not support assert because it accesses [north_neighbor] shape"
#endif
        ! stop
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
! #ifdef ARTLESSthis%rank,": and this%north_boundary =", this%north_boundary, ": north_neighbor =", north_neighbor
! #endif
      if (.not. this%north_boundary) then ! send to north_neighbor
#ifdef ARTLESSISEND
        print *, "MPI_Isend:",this%rank,"to",north_neighbor - 1
#endif
        len = size(this%local(:,:,n-halo_size*2+1:n-halo_size))
        print *, "MPI_Isend:",this%rank,"to",north_neighbor, " OF SIZE", len
        call MPI_Isend(this%local(:,:,n-halo_size*2+1:n-halo_size), &
                       len, MPI_Double_precision, north_neighbor-1, this%north_tag, &
                       MPI_COMM_WORLD, request, ierr)
      end if
      if (.not. this%south_boundary) then ! get from south_neighbor
#ifdef ARTLESSISEND
        print *, "MPI_Irecv:",this%rank,"gets from",south_neighbor - 1
#endif
        len = size(this%halo_south_in(:nx,:,1:halo_size))
        call MPI_Irecv(this%halo_south_in(:nx,:,1:halo_size), len, &
                       MPI_Double_precision, south_neighbor-1, this%south_tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, south)
        ! call MPI_Wait(this%north_request, status, ierr)
      end if
#ifdef ARTLESS
  print *, "end of put_north in rank", this%rank
#endif
    ! call MPI_Finalize(ierr)
    ! call exit
  end subroutine

  module subroutine put_south(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: start, nx, len, ierr
      start = lbound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
#ifdef VERBOSE
        write (*,*) "implementation does not support assert because it accesses [south_neighbor] shape"
#endif
        ! stop
        ! call assert( shape(this%halo_north_in(:nx,:,1:halo_size)[south_neighbor]) &
        !              == shape(this%local(:,:,start:start+halo_size-1)), &
        !              "put_south: conformable halo_north_in and local " )
      end if
      !dir$ pgas defer_sync
      ! this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
      len = size(this%halo_north_in(1:nx,:,1:halo_size))
#ifdef ARTLESS
      print *, this%rank, ": PUT_SOUTH, len =", len
#endif
      if (.not. this%south_boundary) then ! send to south_neighbor
        call MPI_Isend(this%local(:,:,start+halo_size:start+halo_size*2-1), &
                       len, MPI_Double_precision, south_neighbor-1, this%south_tag, &
                       MPI_COMM_WORLD, request, ierr)
      end if
      if (.not. this%north_boundary) then ! get from north_neighbor
        call MPI_Irecv(this%halo_north_in(1:nx,:,1:halo_size), len, &
                       MPI_Double_precision, north_neighbor-1, this%north_tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, north)
      end if

#ifdef ARTLESS
      print *, this%rank, ": end PUT_SOUTH"
#endif
  end subroutine

  module subroutine retrieve_north_halo(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Status) :: status
      integer :: n, nx, ierr
      print *, "north start"
      n = ubound(this%local,3)
      nx = size(this%local,1)
      call MPI_Wait(this%north_request, status, ierr)
      this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
      print *, "retrieve_north_halo complete for north_neighbor", north_neighbor
  end subroutine

  module subroutine retrieve_south_halo(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Status) :: status
      integer :: start, nx, ierr
      print *, "south start"
      start = lbound(this%local,3)
      nx = size(this%local,1)
      call MPI_Wait(this%south_request, status, ierr)
      this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
      print *, "retrieve_south_halo complete for south_neighbor", south_neighbor
  end subroutine

  module subroutine put_east(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Request) :: request
      integer :: n, ny, len, ierr
      n = ubound(this%local,1)
      ny = size(this%local,3)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
#ifdef VERBOSE
        write (*,*) "implementation does not support assert because it accesses [east_neighbor] shape"
#endif
        ! stop
        ! call assert( shape(this%halo_west_in(1:halo_size,:,:ny)[east_neighbor])       &
        !              == shape(this%local(n-halo_size*2+1:n-halo_size,:,:)), &
        !              "put_east: conformable halo_west_in and local " )
      end if

      !dir$ pgas defer_sync
      ! this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)

      len = size(this%halo_west_in(1:halo_size,:,1:ny))
      if (.not. this%east_boundary) then ! send to east_neighbor
        call MPI_Isend(this%local(n-halo_size*2+1:n-halo_size,:,:), &
                       len, MPI_Double_precision, east_neighbor-1, this%east_tag, &
                       MPI_COMM_WORLD, request, ierr)
      end if
      if (.not. this%west_boundary) then ! get from west_neighbor
        call MPI_Irecv(this%halo_west_in(1:halo_size,:,1:ny), len, &
                       MPI_Double_precision, west_neighbor-1, this%west_tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, west)
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
#ifdef VERBOSE
        write (*,*) "implementation does not support assert because it accesses [west_neighbor] shape"
#endif
        ! stop

        ! call assert( shape(this%halo_east_in(1:halo_size,:,:ny)[west_neighbor])               &
        !              == shape(this%local(start+halo_size:start+halo_size*2-1,:,:)), &
        !              "put_west: conformable halo_east_in and local " )
      end if

      !dir$ pgas defer_sync
      ! this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)
      len = size(this%halo_east_in(1:halo_size,:,1:ny))
      if (.not. this%west_boundary) then ! send to west_neighbor
        call MPI_Isend(this%local(start+halo_size:start+halo_size*2-1,:,:), &
                       MPI_Double_precision, west_neighbor-1, this%west_tag, &
                       MPI_COMM_WORLD, request, ierr)
      end if
      if (.not. this%east_boundary) then ! get from east_neighbor
        call MPI_Irecv(this%halo_east_in(1:halo_size,:,1:ny), len, &
                       MPI_Double_precision, east_neighbor-1, this%east_tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, east)
      end if
  end subroutine

  module subroutine retrieve_east_halo(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Status) :: status
      integer :: n, ny, ierr

      n = ubound(this%local,1)
      ny = size(this%local,3)
      call MPI_Wait(this%east_request, status, ierr)
      this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_halo(this)
      class(exchangeable_t), intent(inout) :: this
      type(MPI_Status) :: status
      integer :: start, ny, ierr

      start = lbound(this%local,1)
      ny = size(this%local,3)
      call MPI_Wait(this%west_request, status, ierr)
      this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
  end subroutine

end submodule
