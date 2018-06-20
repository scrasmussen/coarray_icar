#define ARTLESSISEND 1
#define ARTLESSNORTHSOUTH 1
#define ARTLESSEASTWEST 1
submodule(exchangeable_interface) exchangeable_implementation
  use mpi, only : MPI_Comm_rank, MPI_Barrier, &
                      MPI_STATUSES_IGNORE, MPI_COMM_WORLD, MPI_STATUS_SIZE, &
                      MPI_Type_size, MPI_Real
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

        print *, this%rank, ":north=",north_neighbor, "south=",south_neighbor,"east=",east_neighbor,"west=",west_neighbor

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
    call this%put_north
    call this%put_south
#ifdef ARTLESSEASTWEST
    call this%put_east
    call this%put_west
#endif
  end subroutine

  module subroutine retrieve(this, no_sync)
    class(exchangeable_t), intent(inout) :: this
    logical,               intent(in),   optional :: no_sync
    integer :: ierr

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    print *, "======AFTER RETRIEVALS======="
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

    call this%exchange_north

    ! call this%put_north
    ! call this%put_south
    ! call this%put_east
    ! call this%put_west
    ! print *, "======++THIS SHOULD NOT BE TOUCHED++++======="
    ! if (.not. this%north_boundary) call this%retrieve_north_halo
    ! if (.not. this%south_boundary) call this%retrieve_south_halo
    ! if (.not. this%east_boundary)  call this%retrieve_east_halo
    ! if (.not. this%west_boundary)  call this%retrieve_west_halo

    ! call MPI_Barrier(MPI_Comm_world,ierr)
  end subroutine

  module subroutine save_request(this, request, direction)
    implicit none
    class(exchangeable_t), intent(inout) :: this
    integer, intent(in) :: request, direction
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

  module subroutine get_tag(this, sendrecv, from, to, tag)
    implicit none
    class(exchangeable_t), intent(inout) :: this
    integer, intent(in)  :: sendrecv, from, to
    integer, intent(out) :: tag

    if (sendrecv == 0) then !send
      tag = from * 10 + to
    else if (sendrecv == 1) then !recv
      tag = to * 10 + from
    else
      print *, "ERROR in get_tag: sendrecv out of bounds"
      tag = -1
    end if
  end subroutine


  module subroutine put_north(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: request
      integer :: n, nx, len, ierr
      integer :: mpir, tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      real :: r

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
      ! print *, "@@@@@@@@@@@@@BEFORE NORTH ISEND?RECV@@@@@@@@", this%north_boundary, this%south_boundary
      len = size(this%local(:,:,n-halo_size*2+1:n-halo_size))
      if (.not. this%north_boundary) then ! send to north_neighbor
        len = size(this%local(:,:,n-halo_size*2+1:n-halo_size))
        call this%get_tag(sr%send, this%rank, north_neighbor, tag)

        print *, "~~~~~MPI_Isend:",this%rank-1,"to north wit tag ",tag
        call MPI_Isend(this%local(:,:,n-halo_size*2+1:n-halo_size), &
                       len, MPI_Real, north_neighbor-1, tag, &
                       MPI_COMM_WORLD, request, ierr)
        if (ierr .ne. 0) print *, this%rank-1,":*****ERROR MPI_Isend***** 285", ierr
      end if
      if (.not. this%south_boundary) then ! get from south_neighbor
        len = size(this%halo_south_in(:nx,:,1:halo_size))
        call this%get_tag(sr%recv, this%rank, south_neighbor, tag)
        print *, "~~~~~MPI_Irecv:",this%rank-1,"from south wit tag ",tag
        call MPI_Irecv(this%halo_south_in(:nx,:,1:halo_size), len, &
                       MPI_Real, south_neighbor-1, tag, MPI_COMM_WORLD, &
                       request, ierr)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 296", ierr
        call this%save_request(request, south)
        ! call MPI_Wait(this%north_request, status, ierr)
      end if

  ! print *, "end of put_north in rank", this%rank
  end subroutine

  module subroutine put_south(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: request
      integer :: start, nx, len, ierr, tag
      type(sendrecv_t) :: sr
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
#ifdef ARTLESSISEND
        print *, "MPI_Isend:",this%rank-1,"to south",south_neighbor - 1
#endif
        call this%get_tag(sr%send, this%rank, south_neighbor, tag)
        call MPI_Isend(this%local(:,:,start+halo_size:start+halo_size*2-1), &
                       len, MPI_Real, south_neighbor-1, tag, &
                       MPI_COMM_WORLD, request, ierr)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Isend***** 338", ierr
      end if
      if (.not. this%north_boundary) then ! get from north_neighbor
#ifdef ARTLESSISEND
        print *, "MPI_Irecv:",this%rank-1,"from north",north_neighbor - 1
#endif
        call this%get_tag(sr%recv, this%rank, north_neighbor, tag)
        call MPI_Irecv(this%halo_north_in(1:nx,:,1:halo_size), len, &
                       MPI_Real, north_neighbor-1, tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, north)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 350", ierr
      end if

#ifdef ARTLESS
      print *, this%rank, ": end PUT_SOUTH"
#endif
  end subroutine

  module subroutine retrieve_north_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: n, nx, ierr, status(MPI_STATUS_SIZE)
      ! print *, "north start"
      n = ubound(this%local,3)
      nx = size(this%local,1)
      print *, this%rank, ": has north_request", this%north_request
      call MPI_Wait(this%north_request, status, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 366", ierr
      this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
      print *, "retrieve_north_halo complete for north_neighbor", north_neighbor
  end subroutine

  module subroutine retrieve_south_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: start, nx, ierr, status(MPI_STATUS_SIZE)
      ! print *, "south start"
      start = lbound(this%local,3)
      nx = size(this%local,1)
      print *, this%rank, ": has south_request", this%south_request
      call MPI_Wait(this%south_request, status, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 379", ierr
      this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
      print *, "retrieve_south_halo complete for south_neighbor", south_neighbor
  end subroutine

  module subroutine put_east(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: request
      integer :: n, ny, len, ierr, tag
      type(sendrecv_t) :: sr
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
        call this%get_tag(sr%send, this%rank, east_neighbor, tag)
#ifdef ARTLESSISEND
        print *, "MPI_Isend:",this%rank-1,"to east",east_neighbor - 1, "tag", tag
#endif
        call MPI_Isend(this%local(n-halo_size*2+1:n-halo_size,:,:), &
                       len, MPI_Real, east_neighbor-1, tag, &
                       MPI_COMM_WORLD, request, ierr)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Isend***** 411", ierr
      end if
      if (.not. this%west_boundary) then ! get from west_neighbor
        call this%get_tag(sr%recv, this%rank, west_neighbor, tag)
#ifdef ARTLESSISEND
        print *, "MPI_Irecv:",this%rank-1,"from west",west_neighbor - 1, "tag", tag
#endif
        call MPI_Irecv(this%halo_west_in(1:halo_size,:,1:ny), len, &
                       MPI_Real, west_neighbor-1, tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, west)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 424", ierr
      end if
  end subroutine

  module subroutine put_west(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: request
      integer :: start, ny, len, ierr, tag
      type(sendrecv_t) :: sr
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
        call this%get_tag(sr%send, this%rank, west_neighbor, tag)
#ifdef ARTLESSISEND
        print *, "MPI_Isend:",this%rank-1,"to west",west_neighbor - 1, "tag", tag
#endif
        call MPI_Isend(this%local(start+halo_size:start+halo_size*2-1,:,:), &
                       len, MPI_Real, west_neighbor-1, tag, &
                       MPI_COMM_WORLD, request, ierr)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Isend***** 456", ierr
      end if
      if (.not. this%east_boundary) then ! get from east_neighbor
        call this%get_tag(sr%recv, this%rank, east_neighbor, tag)
#ifdef ARTLESSISEND
        print *, "MPI_Irecv:",this%rank-1,"from east",east_neighbor - 1, "tag", tag
#endif
        call MPI_Irecv(this%halo_east_in(1:halo_size,:,1:ny), len, &
                       MPI_Real, east_neighbor-1, tag, MPI_COMM_WORLD, &
                       request, ierr)
        call this%save_request(request, east)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 470", ierr
      end if
  end subroutine

  module subroutine retrieve_east_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: n, ny, ierr, status(MPI_STATUS_SIZE)

      n = ubound(this%local,1)
      ny = size(this%local,3)
      call MPI_Wait(this%east_request, status, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 484", ierr
      this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: start, ny, ierr, status(MPI_STATUS_SIZE)

      start = lbound(this%local,1)
      ny = size(this%local,3)
      call MPI_Wait(this%west_request, status, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 496", ierr
      this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
  end subroutine

    module subroutine exchange_north(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, nx, len, ierr
      integer :: mpir, tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr

      n = ubound(this%local,3)
      nx = size(this%local,1)

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
        call this%get_tag(sr%send, this%rank, north_neighbor, tag)
        ! print *, "~~~~~MPI_send:",this%rank-1,"to north wit tag ",tag
        call MPI_Isend(this%local(:,:,n-halo_size*2+1:n-halo_size), &
                       len, MPI_Real, north_neighbor-1, tag, &
                       MPI_COMM_WORLD, send_request, ierr)
        if (ierr .ne. 0) print *, this%rank-1,":*****ERROR MPI_Isend***** 285", ierr
      end if
      if (.not. this%south_boundary) then ! get from south_neighbor
        call this%get_tag(sr%recv, this%rank, south_neighbor, tag)
        ! print *, "~~~~~MPI_recv:",this%rank-1,"from south wit tag ",tag
        call MPI_Irecv(this%halo_south_in(:nx,:,1:halo_size), len, &
                       MPI_Real, south_neighbor-1, tag, MPI_COMM_WORLD, &
                       recv_request, ierr)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 296", ierr
      end if

      ! Wait for complete of exchange
      if (.not. this%north_boundary) then ! send to north_neighbor
        call MPI_Wait(send_request, status, ierr)
      end if
      if (.not. this%south_boundary) then ! get from south_neighbor
        call MPI_Wait(recv_request, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (.not. this%south_boundary) then ! get from south_neighbor
        this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
        start = lbound(this%local,3)
        this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
      end if
    end subroutine

    module subroutine exchange_south(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, nx, len, ierr
      integer :: mpir, tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      start = lbound(this%local,3)
      nx = size(this%local,1)

      !dir$ pgas defer_sync
      ! this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
      len = size(this%halo_north_in(1:nx,:,1:halo_size))
      if (.not. this%south_boundary) then ! send to north_neighbor
        call this%get_tag(sr%send, this%rank, south_neighbor, tag)
        call MPI_Isend(this%local(:,:,start+halo_size:start+halo_size*2-1), &
                       len, MPI_Real, south_neighbor-1, tag, &
                       MPI_COMM_WORLD, send_request, ierr)
        if (ierr .ne. 0) print *, this%rank-1,":*****ERROR MPI_Isend***** 285", ierr
      end if
      if (.not. this%north_boundary) then ! get from north_neighbor
        call this%get_tag(sr%recv, this%rank, north_neighbor, tag)
        call MPI_Irecv(this%halo_north_in(1:nx,:,1:halo_size), len, &
                       MPI_Real, north_neighbor-1, tag, MPI_COMM_WORLD, &
                       recv_request, ierr)
        if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 296", ierr
      end if

      ! Wait for complete of exchange
      if (.not. this%south_boundary) then ! send to north_neighbor
        call MPI_Wait(send_request, status, ierr)
      end if
      if (.not. this%north_boundary) then ! get from south_neighbor
        call MPI_Wait(recv_request, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (.not. this%north_boundary) then ! get from south_neighbor
        this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
      end if
    end subroutine

    module subroutine exchange_east(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine exchange_west(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

end submodule
