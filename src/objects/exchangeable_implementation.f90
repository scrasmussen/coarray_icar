! #define DEBUG
#define ALTSTART
submodule(exchangeable_interface) exchangeable_implementation
  ! use mpi, only : MPI_Comm_rank, MPI_Barrier, &
  !                     MPI_STATUSES_IGNORE, MPI_COMM_WORLD, MPI_STATUS_SIZE, &
  !                     MPI_Type_size, MPI_Real, MPI_FORTRAN_ORDER
  use mpi
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

    call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, ierr)
    this%rank = this%rank + 1
#ifdef DEBUG
print *,this%rank,":",ims,":",ime,",",kms,":",kme,",",jms,":",jme
#endif
    end associate

    allocate( this%halo_south_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    ), source=initial_value)
    allocate( this%halo_north_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    ), source=initial_value)
    allocate( this%halo_east_in(    halo_size,     grid%halo_nz, grid%ew_halo_ny+halo_size*2), source=initial_value)
    allocate( this%halo_west_in(    halo_size,     grid%halo_nz, grid%ew_halo_ny+halo_size*2), source=initial_value)

    ! call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, ierr)
    ! this%rank = this%rank + 1
! #ifdef DEBUG
! print *,this%rank,":",grid%ns_halo_nx+halo_size*2,grid%halo_nz,halo_size,"EAST",halo_size,grid%halo_nz, grid%ew_halo_ny+halo_size*2
! #endif

    ! set up the neighbors array so we can sync with our neighbors when needed
    if (.not.allocated(neighbors)) then
      associate(me=>this%rank)
        south_neighbor = me - grid%ximages
        north_neighbor = me + grid%ximages
        east_neighbor  = me + 1
        west_neighbor  = me - 1

! #ifdef DEBUG
!         print *, this%rank, ":north=",north_neighbor, "south=",south_neighbor,"east=",east_neighbor,"west=",west_neighbor
! #endif

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


  module subroutine retrieve(this, no_sync)
    class(exchangeable_t), intent(inout) :: this
    logical,               intent(in),   optional :: no_sync
    integer :: ierr
    ! print *, "WARNING: retrieve is not yet implemented"

!@@@@#    ARTLESS
    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary) call this%retrieve_east_halo
    if (.not. this%west_boundary) call this%retrieve_west_halo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine

  module subroutine exchange(this)
    class(exchangeable_t), intent(inout) :: this
    integer :: ierr

    call this%exchange_north
    call this%exchange_south
    call this%exchange_east
    call this%exchange_west

    call MPI_Barrier(MPI_Comm_world,ierr)
  end subroutine

  module subroutine get_tag(this, sendrecv, from, to, tag)
    implicit none
    class(exchangeable_t), intent(inout) :: this
    integer, intent(in)  :: sendrecv, from, to
    integer, intent(out) :: tag

    if (sendrecv == 0) then !send
      tag = from * 10000 + to
    else if (sendrecv == 1) then !recv
      tag = to * 10000 + from
    else
      print *, "ERROR in get_tag: sendrecv out of bounds"
      tag = -1
    end if
  end subroutine

  module subroutine retrieve_north_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: n, nx, ierr, status(MPI_STATUS_SIZE)
      n = ubound(this%local,3)
      nx = size(this%local,1)

        this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine retrieve_south_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: start, nx, ierr, status(MPI_STATUS_SIZE)
      start = lbound(this%local,3)
      nx = size(this%local,1)

        this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
  end subroutine


  module subroutine retrieve_east_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: n, ny, ierr, status(MPI_STATUS_SIZE)
      n = ubound(this%local,1)
      ny = size(this%local,3)
      this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_halo(this)
      class(exchangeable_t), intent(inout) :: this
      integer :: start, ny, ierr, status(MPI_STATUS_SIZE)
      start = lbound(this%local,1)
      ny = size(this%local,3)
      this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
  end subroutine

    module subroutine exchange_north(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: ierr, status(MPI_STATUS_SIZE)
      if (.not. this%south_boundary) then ! get from south_neighbor
         call this%recv_south
      end if
      if (.not. this%north_boundary) then ! send to north_neighbor
         call this%send_north
      end if
      if (.not. this%south_boundary) then ! get from south_neighbor
         call MPI_Wait(this%recv_south_request, status, ierr)
      end if
      if (.not. this%north_boundary) then ! send to north_neighbor
         call MPI_Wait(this%send_north_request, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine

    module subroutine exchange_south(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: ierr, status(MPI_STATUS_SIZE)
      if (.not. this%north_boundary) then ! get from north_neighbor
         call this%recv_north
      end if
      if (.not. this%south_boundary) then ! send to south_neighbor
         call this%send_south
      end if
      if (.not. this%north_boundary) then ! get from north_neighbor
         call MPI_Wait(this%recv_north_request, status, ierr)
      end if
      if (.not. this%south_boundary) then ! send to south_neighbor
         call MPI_Wait(this%send_south_request, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine

    module subroutine exchange_east(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: ierr, status(MPI_STATUS_SIZE)
      if (.not. this%west_boundary) then ! get from west_neighbor
         call this%recv_west
      end if
      if (.not. this%east_boundary) then ! send to east_neighbor
         call this%send_east
      end if
      if (.not. this%west_boundary) then ! get from west_neighbor
         call MPI_Wait(this%recv_west_request, status, ierr)
      end if
      if (.not. this%east_boundary) then ! send to east_neighbor
         call MPI_Wait(this%send_east_request, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine

    module subroutine exchange_west(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: ierr, status(MPI_STATUS_SIZE)
      if (.not. this%east_boundary) then ! get from east_neighbor
         call this%recv_east
      end if
      if (.not. this%west_boundary) then ! send to west_neighbor
         call this%send_west
      end if
      if (.not. this%east_boundary) then ! get from east_neighbor
         call MPI_Wait(this%recv_east_request, status, ierr)
      end if
      if (.not. this%west_boundary) then ! send to west_neighbor
         call MPI_Wait(this%send_west_request, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine


    module subroutine send_north(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, nx, len, ierr
      integer :: tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      integer :: n2, n3, n1, starts2(3) !ARTLESS


      n = ubound(this%local,3)
      n1 = ubound(this%local,1)
      n2 = ubound(this%local,2)
      n3 = ubound(this%local,3)
      nx = size(this%local,1)

#ifdef DEBUG
      if (size(this%local(:,:,n-halo_size*2+1:n-halo_size)) .NE. &
          size(this%halo_south_in(:nx,:,1:halo_size))) then
        print*, "ERROR: the number of elements in send and recv arrays differs"
      end if
#endif

      !dir$ pgas defer_sync
      ! this%halo_south_in(1:nx,:,1:halo_size)[north_neighbor] = this%local(:,:,n-halo_size*2+1:n-halo_size)


      call this%get_tag(sr%send, this%rank, north_neighbor, tag)
      
      if (this%send_north_type == 0) then 
         sizes = shape(this%local(:,:,:))
         subsizes = shape(this%local(:,:,n-halo_size*2+1:n-halo_size))
         starts = (/1-1,1-1,n-halo_size*2+1-1/)
#ifdef ALTSTART
         starts = (/1-1,1-1,sizes(3)-halo_size*2+1-1/)
#endif
#ifdef DEBUG
         print *, "sendnorth307 rank:", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts, " n =", n
         ! print *, "rank:", this%rank, "n1 = ", n1, "n2 = ", n2, "n3 = ", n3
         ! print *, this%rank, "::SHAPE=",sizes, "SIZE=", size(this%local)
#endif
         call MPI_Type_create_subarray(3, sizes, subsizes, &
              starts, MPI_ORDER_FORTRAN, MPI_Real, this%send_north_type, ierr)
         call MPI_Type_commit(this%send_north_type, ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray succes 307"
#endif
      end if
      
      call MPI_Isend(this%local, &
                     1, this%send_north_type, north_neighbor-1, tag, &
                     MPI_COMM_WORLD, send_request, ierr)
      ! call MPI_Isend(this%local(:,:,n-halo_size*2+1:n-halo_size), &
      !                len, MPI_Real, north_neighbor-1, tag, &
      !                MPI_COMM_WORLD, send_request, ierr)
      if (ierr .ne. 0) print *, this%rank-1,":*****ERROR MPI_Isend***** 285", ierr

      this%send_north_request = send_request
    end subroutine

    module subroutine recv_south(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: recv_request
      integer :: start, n, nx, len, len2, ierr
      integer :: tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      !dir$ pgas defer_sync
      ! this%halo_south_in(1:nx,:,1:halo_size)[north_neighbor] = this%local(:,:,n-halo_size*2+1:n-halo_size)

      n = ubound(this%local,3)
      nx = size(this%local,1)
      len = size(this%halo_south_in(:nx,:,1:halo_size))

      if (this%recv_south_type == 0) then 
         sizes = shape(this%halo_south_in(:,:,:))
         subsizes = shape(this%halo_south_in(:nx,:,1:halo_size))
         starts = (/1-1,1-1,1-1/)
#ifdef DEBUG
         print *, "revsouth345", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
              MPI_Real,this%recv_south_type,ierr)
         call MPI_Type_commit(this%recv_south_type,ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray success 345"
#endif
      end if

      ! print *, this%rank-1,":", len, "and", len2   ARTLESS
      call this%get_tag(sr%recv, this%rank, south_neighbor, tag)
      call MPI_Irecv(this%halo_south_in, 1, this%recv_south_type, & 
                     south_neighbor-1, tag, MPI_COMM_WORLD, &
                     recv_request, ierr)
      ! call MPI_Irecv(this%halo_south_in(:nx,:,1:halo_size), len, &
      !                MPI_Real, south_neighbor-1, tag, MPI_COMM_WORLD, &
      !                recv_request, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 296", ierr
      this%recv_south_request = recv_request
    end subroutine

    module subroutine send_south(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, nx, len,len2, ierr
      integer :: tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      start = lbound(this%local,3)
      nx = size(this%local,1)

      !dir$ pgas defer_sync
      ! this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
      len = size(this%halo_north_in(1:nx,:,1:halo_size))
      ! len2 = size(this%local(:,:,start+halo_size:start+halo_size*2-1))
      ! print *, this%rank, ":",len," to ", len2

      if (this%send_south_type == 0) then 
         sizes = shape(this%local(:,:,:))
         subsizes = shape(this%local(:,:,start+halo_size:start+halo_size*2-1))
         starts = (/1-1,1-1,start+halo_size-1/)
#ifdef ALTSTART
         starts = (/1-1,1-1,halo_size-1/)
#endif
#ifdef DEBUG
         print *, "sendsouth402", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3, sizes, subsizes, &
              starts, MPI_ORDER_FORTRAN, MPI_Real, this%send_south_type, ierr)
         call MPI_Type_commit(this%send_south_type, ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray success 386"
#endif
      end if

      call this%get_tag(sr%send, this%rank, south_neighbor, tag)
      call MPI_Isend(this%local, &
                     1, this%send_south_type, south_neighbor-1, tag, &
                     MPI_COMM_WORLD, send_request, ierr)
      ! call MPI_Isend(this%local(:,:,start+halo_size:start+halo_size*2-1), &
      !                len, MPI_Real, south_neighbor-1, tag, &
      !                MPI_COMM_WORLD, send_request, ierr)
      if (ierr .ne. 0) print *, this%rank-1,":*****ERROR MPI_Isend***** 285", ierr
      this%send_south_request = send_request
    end subroutine

    module subroutine recv_north(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, nx, len, ierr
      integer :: tag, status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      start = lbound(this%local,3)
      nx = size(this%local,1)

      !dir$ pgas defer_sync
      ! this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
      len = size(this%halo_north_in(1:nx,:,1:halo_size))

      if (this%recv_north_type == 0) then 
         sizes = shape(this%halo_north_in(:,:,:))
         subsizes = shape(this%halo_north_in(1:nx,:,1:halo_size))
         starts = (/1-1,1-1,1-1/)
#ifdef DEBUG
         print *, "recvnorth424", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
              MPI_Real,this%recv_north_type,ierr)
         call MPI_Type_commit(this%recv_north_type,ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray success 424"
#endif
      end if


      call this%get_tag(sr%recv, this%rank, north_neighbor, tag)
      ! print *, this%rank, ": tag =", tag, "nx=",nx,"halo_size=",halo_size,"len",len
      call MPI_Irecv(this%halo_north_in, 1, &
                     this%recv_north_type, north_neighbor-1, tag, MPI_COMM_WORLD, &
                     recv_request, ierr)
      ! call MPI_Irecv(this%halo_north_in(1:nx,:,1:halo_size), len, &
      !                MPI_Real, north_neighbor-1, tag, MPI_COMM_WORLD, &
      !                recv_request, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 296", ierr
      this%recv_north_request = recv_request
    end subroutine

    module subroutine send_east(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, ny, len, ierr, tag
      integer :: status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      start = lbound(this%local,1)
      n = ubound(this%local,1)
      ny = size(this%local,3)

      !dir$ pgas defer_sync
      ! this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)
      len = size(this%halo_west_in(1:halo_size,:,1:ny))

      if (this%send_east_type == 0) then 
         sizes = shape(this%local(:,:,:))
         subsizes = shape(this%local(n-halo_size*2+1:n-halo_size,:,:))
         starts = (/n-halo_size*2+1-1,1-1,1-1/)
#ifdef ALTSTART
         starts = (/sizes(1)-halo_size*2+1-1,1-1,1-1/)
#endif
#ifdef DEBUG
         print *, "sendeast490", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3, sizes, subsizes, &
              starts, MPI_ORDER_FORTRAN, MPI_Real, this%send_east_type, ierr)
         call MPI_Type_commit(this%send_east_type, ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray success 465"
#endif
      end if

      call this%get_tag(sr%send, this%rank, east_neighbor, tag)
      call MPI_Isend(this%local, &
                     1, this%send_east_type, east_neighbor-1, tag, &
                     MPI_COMM_WORLD, send_request, ierr)
      ! call MPI_Isend(this%local(n-halo_size*2+1:n-halo_size,:,:), &
      !                len, MPI_Real, east_neighbor-1, tag, &
      !                MPI_COMM_WORLD, send_request, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Isend***** 411", ierr

      this%send_east_request = send_request
      ! this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
    end subroutine

    module subroutine recv_west(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, ny, len, ierr, tag
      integer :: status(MPI_STATUS_SIZE)
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      start = lbound(this%local,1)
      n = ubound(this%local,1)
      ny = size(this%local,3)

      !dir$ pgas defer_sync
      ! this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)
      len = size(this%halo_west_in(1:halo_size,:,1:ny))

      if (this%recv_west_type == 0) then 
         sizes = shape(this%halo_west_in(:,:,:))
         subsizes = shape(this%halo_west_in(1:halo_size,:,1:ny))
         starts = (/1-1,1-1,1-1/)
#ifdef DEBUG
         print *, "recvwest506", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3, sizes, subsizes, &
              starts, MPI_ORDER_FORTRAN, MPI_Real, this%recv_west_type, ierr)
         call MPI_Type_commit(this%recv_west_type, ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray success 506"
#endif
      end if


      call this%get_tag(sr%recv, this%rank, west_neighbor, tag)
      call MPI_Irecv(this%halo_west_in, 1, &
                     this%recv_west_type, west_neighbor-1, tag, MPI_COMM_WORLD, &
                     recv_request, ierr)
      ! call MPI_Irecv(this%halo_west_in(1:halo_size,:,1:ny), len, &
      !                MPI_Real, west_neighbor-1, tag, MPI_COMM_WORLD, &
      !                recv_request, ierr)
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 424", ierr

      this%recv_west_request = recv_request
      ! this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
    end subroutine

    module subroutine send_west(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, ny, len, ierr, tag
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      start = lbound(this%local,1)
      ny = size(this%local,3)

      !dir$ pgas defer_sync
      ! this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)
      len = size(this%halo_east_in(1:halo_size,:,1:ny))

      if (this%send_west_type == 0) then 
         sizes = shape(this%local(:,:,:))
         subsizes = shape(this%local(start+halo_size:start+halo_size*2-1,:,:))
         starts = (/start-halo_size-1,1-1,1-1/)
#ifdef ALTSTART
         starts = (/1-1,1-1,1-1/)
#endif

#ifdef DEBUG
         print *, "sendwest490", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3, sizes, subsizes, &
              starts, MPI_ORDER_FORTRAN, MPI_Real, this%send_west_type, ierr)
         call MPI_Type_commit(this%send_west_type, ierr)
! #ifdef DEBUG
!          print *, "MPI_Type_create_subarray success 546"
! #endif
      end if

      call this%get_tag(sr%send, this%rank, west_neighbor, tag)
      call MPI_Isend(this%local, &
                     1, this%send_west_type, west_neighbor-1, tag, &
                     MPI_COMM_WORLD, send_request, ierr)

      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 470", ierr

      this%send_west_request = send_request
      ! n = ubound(this%local,1)
      ! this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
    end subroutine

    module subroutine recv_east(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      integer :: send_request, recv_request
      integer :: start, n, ny, len, ierr, tag
      type(sendrecv_t) :: sr
      integer :: starts(3), sizes(3), subsizes(3)
      start = lbound(this%local,1)
      ny = size(this%local,3)

      !dir$ pgas defer_sync
      ! this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)
      len = size(this%halo_east_in(1:halo_size,:,1:ny))

      if (this%recv_east_type == 0) then 
         sizes = shape(this%halo_east_in(:,:,:))
         subsizes = shape(this%halo_east_in(1:halo_size,:,1:ny))
         starts = (/1-1,1-1,1-1/)
#ifdef DEBUG
         print *, "recveast584", this%rank, "sizes =", sizes, "subsizes =", subsizes, "starts =", starts
#endif
         call MPI_Type_create_subarray(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
              MPI_Real,this%recv_east_type,ierr)
         call MPI_Type_commit(this%recv_east_type,ierr)
#ifdef DEBUG
         print *, "MPI_Type_create_subarray success 584"
#endif
      end if

      call this%get_tag(sr%recv, this%rank, east_neighbor, tag)
      call MPI_Irecv(this%halo_east_in, 1, &
                     this%recv_east_type, east_neighbor-1, tag, MPI_COMM_WORLD, &
                     recv_request, ierr)
      ! call MPI_Irecv(this%halo_east_in(1:halo_size,:,1:ny), 1, &
      if (ierr .ne. 0) print *, this%rank-1, ":*****ERROR MPI_Irecv***** 470", ierr

      this%recv_east_request = recv_request
      ! n = ubound(this%local,1)
      ! this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
    end subroutine

end submodule
