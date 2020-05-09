submodule(convection_exchangeable_interface) &
    convection_exchangeable_implementation
  use assertions_interface, only : assert, assertions
  use grid_interface, only : grid_t
  implicit none

  integer, parameter :: default_buf_size=1
  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_neighbor, south_neighbor, buf_size, halo_size
  integer, save :: east_neighbor, west_neighbor
  integer, save :: northeast_neighbor, northwest_neighbor
  integer, save :: southeast_neighbor, southwest_neighbor

contains
  ! function initialize_convection_array_t(this)
  !   class(convection_exchangeable_array_t), intent(inout) :: this
  ! end function initialize_convection_array_t

  ! constructor
  module subroutine const(this, convection_type_enum, grid, input_buf_size, &
      halo_width, u_in, v_in, w_in, temperature, pressure)
    use iso_c_binding, only: c_int
    class(convection_exchangeable_t), intent(inout) :: this
    type(grid_t) :: grid
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width
    integer(c_int), intent(in) :: convection_type_enum
    real, optional, intent(in) :: u_in,v_in,w_in
    real, dimension(:,:,:), intent(in) :: temperature, pressure
    real :: random_start(3)

    integer :: n_neighbors, current
    integer :: ims,ime,kms,kme,jms,jme,i,j,k
    real :: x,y,z
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
    if (present(input_buf_size)) then
        buf_size = input_buf_size
    else
        buf_size = default_buf_size
    end if
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
      ims = grid%ims - halo_east  ! ARTLESS: no halos right now
      ime = grid%ime + halo_west
      jms = grid%jms - halo_south
      jme = grid%jme + halo_north
      kms = grid%kms
      kme = grid%kme

! -      allocate(this%local(ims:ime,kms:kme,jms:jme))
! -      do i=ims,ime
! -        do j=jms,jme
! -          do k=kms,kme
! -              if (i.eq.5 .and. k.eq.5 .and. j.eq.5) then
! -                print*,"ARTLESS change this to percentage once exchange working"
! -                ! print*, "~~~~~~~~~~~~~", i, k, j, temperature(i,k,j)
! -                this%local(i,k,j) = convection_particle(i,j,k, &
! -                    0.5,0.5,0.0,pressure(i,k,j), temperature(i,k,j))
! -              else
! -                this%local(i,k,j)%exists = .false.
! -              end if
! -
! -            end do
! -          end do
! -        end do


      allocate(this%local(input_buf_size * 4))
      call random_number(random_start)
      if (this_image() .eq. 1) then
        print*,"ARTLESS change this to percentage once exchange working"
        ! pick random place in image to start
        call random_number(random_start)
        x = (ims+halo_east) + (random_start(1) * (ime-ims-halo_west))
        z = kms + (random_start(3) * (kme-kms))
        y = (jms+halo_south) + (random_start(2) * (jme-jms-halo_north))
        print *, "BOUNDS::", ims, ime, kms, kme, jms, jme
        print *, "HALO::", halo_north ,halo_south ,halo_east ,halo_west
        ! print *, x,z,y

        this%local(1) = convection_particle(x,y,z, 0.5,0.5,0.0,&
            pressure(floor(x),floor(z),floor(y)), &
            temperature(floor(x),floor(z),floor(y)))
        ! need to set exists to false
      end if

      if (this_image() .eq. 3) then
        print*,"ARTLESS change this to percentage once exchange working"
        call random_number(random_start)
        x = (ims+halo_east) + (random_start(1) * (ime-ims-halo_west))
        z = kms + (random_start(3) * (kme-kms))
        y = (jms+halo_south) + (random_start(2) * (jme-jms-halo_north))

        this%local(1) = convection_particle(x,y,z, 0.5,0.5,0.0,&
            pressure(floor(x),floor(y),floor(z)), &
            temperature(floor(x),floor(y),floor(z)))
      ! need to set exists to false
        ! print *, "-----", ims, ime, kms, kme, jms, jme
        ! print *, "----- shape = ", shape(this%local)
      end if

    end associate
    print *, "ALLOCATING BUFFERS OF SIZE", buf_size
    allocate( this%buf_south_in(buf_size)[*])
    allocate( this%buf_north_in(buf_size)[*])
    allocate( this%buf_east_in(buf_size)[*])
    allocate( this%buf_west_in(buf_size)[*])

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
        print *, this_image(), ": boundry", this%north_boundary, &
            this%south_boundary, this%east_boundary, this%west_boundary
      end associate
    endif



    print *, "=============++YBE YBE YBYE BYE  ============="
  end subroutine


  module subroutine send(this)
    class(convection_exchangeable_t), intent(inout) :: this
    ! if (.not. this%north_boundary) call this%put_north
    ! if (.not. this%south_boundary) call this%put_south
    ! if (.not. this%east_boundary)  call this%put_east
    ! if (.not. this%west_boundary)  call this%put_west
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

    if (.not. this%north_boundary) call this%retrieve_north_buf
    if (.not. this%south_boundary) call this%retrieve_south_buf
    if (.not. this%east_boundary) call this%retrieve_east_buf
    if (.not. this%west_boundary) call this%retrieve_west_buf
  end subroutine

  module subroutine exchange(this)
    class(convection_exchangeable_t), intent(inout) :: this
    ! if (.not. this%north_boundary) call this%put_north
    ! if (.not. this%south_boundary) call this%put_south
    if (.not. this%east_boundary)  call this%put_east
    if (.not. this%west_boundary)  call this%put_west

    sync images( neighbors )

    if (.not. this%north_boundary) call this%retrieve_north_buf
    if (.not. this%south_boundary) call this%retrieve_south_buf
    if (.not. this%east_boundary)  call this%retrieve_east_buf
    if (.not. this%west_boundary)  call this%retrieve_west_buf
  end subroutine

  module subroutine load_buf(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    integer :: i
    particle%exists = .false.
    do i=1,ubound(this%buf_north_in, dim=1)
    end do

    print *, "LOADING BUF SUBROUTINE"
  end subroutine


  module subroutine put_north(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle
      integer :: n, nx

      if (this%north_boundary) return
      print *, "PUTTING NORTH!!!!!!!!!!!!!!!!!!!!!!!!!!", "and ", this%north_i
      ! this%buf_north_in(this%north_i)[south_neighbor] = particle
      this%buf_south_in(this%south_i)[north_neighbor] = particle
      particle%exists=.false.


      ! n = ubound(this%local,3)
      ! nx = size(this%local,1)
      ! if (assertions) then
      !   !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
      !   call assert( shape(this%halo_south_in(:nx,:,1:halo_size)[north_neighbor]) &
      !                == shape(this%local(:,:,n-halo_size+1:n)),         &
      !                "put_north convection: conformable halo_south_in and local " )
      ! end if

      ! !dir$ pgas defer_sync
      ! this%halo_south_in(1:nx,:,1:halo_size)[north_neighbor] = this%local(:,:,n-halo_size*2+1:n-halo_size)
  end subroutine

  module subroutine put_south(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      integer :: start, nx
      print *, "PUTTING SOUTH to ", south_neighbor, "and ", this%north_i
      particle%exists=.false.
      this%buf_north_in(this%north_i)[south_neighbor] = particle
      ! this%north_i = this%north_i + 1


      ! start = lbound(this%local,3)
      ! nx = size(this%local,1)

      ! if (assertions) then
      !   !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
      !   call assert( shape(this%halo_north_in(:nx,:,1:halo_size)[south_neighbor]) &
      !                == shape(this%local(:,:,start:start+halo_size-1)), &
      !                "put_south convection: conformable halo_north_in and local " )
      ! end if
      ! !dir$ pgas defer_sync
      ! this%halo_north_in(1:nx,:,1:halo_size)[south_neighbor] = this%local(:,:,start+halo_size:start+halo_size*2-1)
  end subroutine

  module subroutine retrieve_north_buf(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: n, nx
      ! n = ubound(this%local,3)
      ! nx = size(this%local,1)

      ! this%local(:,:,n-halo_size+1:n) = this%halo_north_in(:nx,:,1:halo_size)
  end subroutine

  module subroutine retrieve_south_buf(this)
      class(convection_exchangeable_t), intent(inout) :: this
      integer :: i,n
      n = ubound(this%buf_south_in, dim=1)
      do i=1,n
        if (this%buf_south_in(i)%exists .eqv. .true.) then
          print* , "RETRIEVE SOUTH BUF--"
          this%buf_south_in(i)%exists=.false.
        end if
      end do


      ! this%buf_south_in(this%south_i) = particle
  end subroutine

  module subroutine put_east(this)
      class(convection_exchangeable_t), intent(inout) :: this
      integer :: n, ny
      ! n = ubound(this%local,1)
      ! ny = size(this%local,3)

      ! if (assertions) then
      !   !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
      !   call assert( shape(this%halo_west_in(1:halo_size,:,:ny)[east_neighbor])       &
      !                == shape(this%local(n-halo_size*2+1:n-halo_size,:,:)), &
      !                "put_east: conformable halo_west_in and local " )
      ! end if

      ! !dir$ pgas defer_sync
      ! this%halo_west_in(1:halo_size,:,1:ny)[east_neighbor] = this%local(n-halo_size*2+1:n-halo_size,:,:)
  end subroutine

  module subroutine put_west(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: start, ny
      ! start = lbound(this%local,1)
      ! ny = size(this%local,3)
      ! ! print *, "START = ", start, "when shape is ", shape(this%local(:,:,:)),
      ! if (assertions) then
      !   !! gfortran 6.3.0 doesn't check coarray shape conformity with -fcheck=all so we verify with an assertion
      !   call assert( shape(this%halo_east_in(1:halo_size,:,:ny)[west_neighbor])               &
      !                == shape(this%local(start+halo_size:start+halo_size*2-1,:,:)), &
      !                "put_west: conformable halo_east_in and local " )
      ! end if

      ! !dir$ pgas defer_sync
      ! this%halo_east_in(1:halo_size,:,1:ny)[west_neighbor] = this%local(start+halo_size:start+halo_size*2-1,:,:)
  end subroutine

  module subroutine retrieve_east_buf(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: n, ny
      ! n = ubound(this%local,1)
      ! ny = size(this%local,3)

      ! this%local(n-halo_size+1:n,:,:) = this%halo_east_in(1:halo_size,:,1:ny)
  end subroutine

  module subroutine retrieve_west_buf(this)
      class(convection_exchangeable_t), intent(inout) :: this

      integer :: start, ny
      ! start = lbound(this%local,1)
      ! ny = size(this%local,3)

      ! this%local(start:start+halo_size-1,:,:) = this%halo_west_in(1:halo_size,:,1:ny)
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
      real :: ws
      integer :: i,j,k, l_bound(1), dif(1), new_ijk(3), me
      !integer :: l_bound(3), dif(3), new_ijk(3)  ! old three dimensional

      ! l_bound = lbound(this%local)
      ! dif  = l_bound - lbound(temperature)
! -      do i=its,ite
! -        do j=jts,jte
! -          do k=kts,kte
! -            associate (particle=>this%local(i,k,j))
! -            if (particle%exists .eqv. .true. .and. &
      me = this_image()
      do i=1,ubound(this%local,1)
        associate (particle=>this%local(i))
          if (particle%exists .eqv. .true.) then
            ! print *, "HERE"
            if (particle%x .lt. its .or. particle%x .gt. ite .or. &
                particle%z .lt. kts .or. particle%z .gt. kte .or. &
                particle%y .lt. jts .or. particle%y .gt. jte) then
              print *, "x:", its, "<", particle%x, "<", ite
              print *, "z:", kts, "<", particle%z, "<", kte
              print *, "y:", jts, "<", particle%y, "<", jte
              stop "x,y,z is out of bounds"
            end if

            !-----------------------------------------------------------------
            ! Handle Buoyancy
            !-----------------------------------------------------------------
            ! print *, me,":xzy=", particle%x, particle%z,particle%y
            ! ! print *, me,":",  its,ite, ":",kts,kte, ":",jts,jte
            ! print *, me,":", shape(temperature), ":", lbound(temperature), &
            !     ":", ubound(temperature)
            T = temperature(floor(particle%x), floor(particle%z), &
                floor(particle%y))

            T_prime = particle%temperature
            a_prime = (T_prime - T) / T * gravity
            displacement =  0    + 0.5 * a_prime * 1 * 1
            particle%z = particle%z + displacement

            ! print *,me, "z = ", particle%z  , "with displacement", displacement
            ! print *, "limits????????", ite, jte, kte

            ! checking to see if need to move particle
            new_ijk = floor((/particle%x,particle%y,particle%z/))
            if (particle%z .gt. kte) then
              particle%exists=.false.
            else if (particle%z .lt. kts) then ! kts will be 1
              particle%exists=.false.
            end if

            !-----------------------------------------------------------------
            ! Handle Windfield
            !-----------------------------------------------------------------
            ws = sqrt(particle%u * particle%u + particle%v * particle%v)

            particle%y = particle%y + ws
            print *, "windspeed is", ws
            print *, "jts  <   y    <  jte"
            print *,  jts, particle%y, jte
            if (particle%y .lt. jts) then
              call this%put_south(particle)
            else if (particle%y .gt. jte) then ! jts will be 1
              call this%put_north(particle)
            end if

            !     call this%load_buf(particle)
            !     ! call this%put_north(particle)


          end if
        end associate
      end do

      ! OLD after moving everything need to reset movement flags to false
      ! do concurrent (i=its:ite, j=jts:jte, k=kts:kte)
      !   this%local(i,k,j)%moved = .false.
      ! end do
    end subroutine

end submodule
