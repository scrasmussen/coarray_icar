submodule(convection_exchangeable_interface) &
    convection_exchangeable_implementation
  use assertions_interface, only : assert, assertions
  use grid_interface, only : grid_t
  implicit none

  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_neighbor, south_neighbor, halo_size
  integer, save :: east_neighbor, west_neighbor
  integer, save :: northeast_neighbor, northwest_neighbor
  integer, save :: southeast_neighbor, southwest_neighbor

contains
  ! function initialize_convection_array_t(this)
  !   class(convection_exchangeable_array_t), intent(inout) :: this
  ! end function initialize_convection_array_t

  ! constructor
  module subroutine const(this, convection_type_enum, grid, halo_width, &
      u_in, v_in, w_in, temperature, pressure)
    use iso_c_binding, only: c_int
    class(convection_exchangeable_t), intent(inout) :: this
    type(grid_t) :: grid
    integer, intent(in), optional :: halo_width
    integer(c_int), intent(in) :: convection_type_enum
    real, optional, intent(in) :: u_in,v_in,w_in
    real, dimension(:,:,:), intent(in) :: temperature, pressure

    integer :: n_neighbors, current
    integer :: ims,ime,kms,kme,jms,jme,i,j,k
    real :: u,v,w
    if (present(u_in)) then
      u = u_in
    else
      u = 0.0
    end if
    if (present(v_in)) then
      v = v_in
    else
      v = 0.0
    end if
    if (present(w_in)) then
      w = w_in
    else
      w = 0.0
    end if
    print *, "=============++HIHI HI H IH HIHI ============="
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

      allocate(this%local(ims:ime,kms:kme,jms:jme))
      do i=ims,ime
        do j=jms,jme
          do k=kms,kme
              if (i.eq.5 .and. k.eq.5 .and. j.eq.5) then
                print*,"ARTLESS change this to percentage once exchange working"
                ! print*, "~~~~~~~~~~~~~", i, k, j, temperature(i,k,j)
                this%local(i,k,j) = convection_particle(i,j,k, &
                    0.5,0.5,0.0,pressure(i,k,j), temperature(i,k,j))
              else
                this%local(i,k,j)%exists = .false.
              end if

            end do
          end do
        end do
        print *, "-----", ims, ime, kms, kme, jms, jme
        print *, "----- shape = ", shape(this%local)
    end associate

    allocate( this%halo_south_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    )[*])
    allocate( this%halo_north_in( grid%ns_halo_nx+halo_size*2, grid%halo_nz,   halo_size    )[*])
    allocate( this%halo_east_in( halo_size, grid%halo_nz, grid%ew_halo_ny+halo_size*2)[*])
    allocate( this%halo_west_in( halo_size, grid%halo_nz, grid%ew_halo_ny+halo_size*2)[*])

    if (this_image() .eq. 1) then
      print *, "===--- ARTLESS NEED TO DOUBLE CHECK NEIGHBORS ---==="
    end if
    ! set up the neighbors array so we can sync with our neighbors when needed
    if (.not.allocated(neighbors)) then
      associate(me=>this_image())
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
        print *, this_image(), ": north", north_neighbor, ": east", &
             east_neighbor, ": south", south_neighbor, &
             ": west", west_neighbor
      end associate
    endif



    print *, "=============++YBE YBE YBYE BYE  ============="
  end subroutine


  module subroutine send(this)
    class(convection_exchangeable_t), intent(inout) :: this
    if (.not. this%north_boundary) call this%put_north
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west
  end subroutine

  module subroutine retrieve(this, no_sync)
    class(convection_exchangeable_t), intent(inout) :: this
    logical,               intent(in),   optional :: no_sync

    if (.not. present(no_sync)) then
        sync images( neighbors )
    else
        if (.not. no_sync) then
            sync images( neighbors )
        endif
    endif

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary) call this%retrieve_east_halo
    if (.not. this%west_boundary) call this%retrieve_west_halo
  end subroutine

  module subroutine exchange(this)
    class(convection_exchangeable_t), intent(inout) :: this
    if (.not. this%north_boundary) call this%put_north
    if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west

    sync images( neighbors )

    if (.not. this%north_boundary) call this%retrieve_north_halo
    if (.not. this%south_boundary) call this%retrieve_south_halo
    if (.not. this%east_boundary)  call this%retrieve_east_halo
    if (.not. this%west_boundary)  call this%retrieve_west_halo
  end subroutine

  module subroutine put_north(this)
      class(convection_exchangeable_t), intent(inout) :: this
      integer :: n, nx
      n = ubound(this%local,3)
      nx = size(this%local,1)
      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_south_in(:nx,:,1:halo_size)[north_neighbor]) &
                     == shape(this%local(:,:,n-halo_size+1:n)),         &
                     "put_north convection: conformable halo_south_in and local " )
      end if

      !dir$ pgas defer_sync
      this%halo_south_in(1:nx,:,1:halo_size)[north_neighbor] = this%local(:,:,n-halo_size*2+1:n-halo_size)
  end subroutine

  module subroutine put_south(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: start, nx
      start = lbound(this%local,3)
      nx = size(this%local,1)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_north_in(:nx,:,1:halo_size)[south_neighbor]) &
                     == shape(this%local(:,:,start:start+halo_size-1)), &
                     "put_south convection: conformable halo_north_in and local " )
      end if
      !dir$ pgas defer_sync
      this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
  end subroutine

  module subroutine retrieve_north_halo(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: n, nx
      n = ubound(this%local,3)
      nx = size(this%local,1)

      this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine retrieve_south_halo(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: start, nx
      start = lbound(this%local,3)
      nx = size(this%local,1)

      this%local(:,:,start:start+halo_size-1) = this%halo_south_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine put_east(this)
      class(convection_exchangeable_t), intent(inout) :: this
      integer :: n, ny
      n = ubound(this%local,1)
      ny = size(this%local,3)

      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_west_in(1:halo_size,:,:ny)[east_neighbor])       &
                     == shape(this%local(n-halo_size*2+1:n-halo_size,:,:)), &
                     "put_east: conformable halo_west_in and local " )
      end if

      !dir$ pgas defer_sync
      this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)
  end subroutine

  module subroutine put_west(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: start, ny
      start = lbound(this%local,1)
      ny = size(this%local,3)
      ! print *, "START = ", start, "when shape is ", shape(this%local(:,:,:)),
      if (assertions) then
        !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
        call assert( shape(this%halo_east_in(1:halo_size,:,:ny)[west_neighbor])               &
                     == shape(this%local(start+halo_size:start+halo_size*2-1,:,:)), &
                     "put_west: conformable halo_east_in and local " )
      end if

      !dir$ pgas defer_sync
      this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)
  end subroutine

  module subroutine retrieve_east_halo(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: n, ny
      n = ubound(this%local,1)
      ny = size(this%local,3)

      this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_halo(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: start, ny
      start = lbound(this%local,1)
      ny = size(this%local,3)

      this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
    end subroutine

    module subroutine process(this, dt, its,ite, jts,jte, kts,kte, &
        temperature)
      implicit none
      class(convection_exchangeable_t), intent(inout) :: this
      real,           intent(in)    :: dt
      integer,        intent(in)    :: its,ite, jts,jte, kts,kte
      real, dimension(:,:,:), intent(in) :: temperature
      real, parameter :: gravity = 9.80665
      real, parameter :: artless_density = 1.003
      real :: a_prime, displacement, t, t_prime
      integer :: i,j,k, l_bound(3), dif(3), new_ijk(3)

      l_bound = lbound(this%local)
      dif  = l_bound - lbound(temperature)

      do i=its,ite
        do j=jts,jte
          do k=kts,kte
            associate (particle=>this%local(i,k,j))
            if (particle%exists .eqv. .true. .and. &
                particle%moved .eqv. .false.) then
              T = temperature(i-dif(0), k-dif(2), j-dif(1))
              T_prime = particle%temperature
              a_prime = (T_prime - T) / T * gravity
              displacement =  0    + 0.5 * a_prime * 1 * 1
              particle%z = particle%z + displacement

              print *, "z = ", particle%z  , "with displacement", displacement

              ! checking to see if need to move particle
              new_ijk = floor((/particle%x,particle%y,particle%z/))

              if ( .not. all( (/i,j,k/) .eq. &
                  new_ijk)) then
                ! print *, "----- MOVING PARTICLE: from",(/i,j,k/), "to",new_ijk
                call particle%move_particle( &
                    this%local(new_ijk(1),new_ijk(3),new_ijk(2)))
              end if
              ! ARTLESS DIF BETWEEN TEMP IMAGES WORRISOME? maybe explained with init
              ! print *, me, ":", temperature(i-dif(0), k-dif(2), j-dif(1))
              ! print *, me, ":", T, "||||", this%local(i,k,j)%temperature
            end if
            end associate
          end do
        end do
      end do

      ! after moving everything need to reset movement flags to false
      do concurrent (i=its:ite, j=jts:jte, k=kts:kte)
        this%local(i,k,j)%moved = .false.
      end do
    end subroutine

end submodule
