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

    if (this%north_boundary) then
      this%northeast_boundary = .true.
      this%northwest_boundary = .true.
    end if
    if (this%south_boundary) then
      this%southeast_boundary = .true.
      this%southwest_boundary = .true.
    end if
    if (this%east_boundary) then
      this%northeast_boundary = .true.
      this%southeast_boundary = .true.
    end if
    if (this%west_boundary) then
      this%northwest_boundary = .true.
      this%southwest_boundary = .true.
    end if

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

      ! --------- this is old code from when particle arrays are used
! -      allocate(this%local(ims:ime,kms:kme,jms:jme))
! -      do i=ims,ime
! -        do j=jms,jme
! -          do k=kms,kme
! -              if (i.eq.5 .and. k.eq.5 .and. j.eq.5) then
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
    allocate( this%buf_north_in(buf_size)[*])
    allocate( this%buf_south_in(buf_size)[*])
    allocate( this%buf_east_in(buf_size)[*])
    allocate( this%buf_west_in(buf_size)[*])
    allocate( this%buf_northeast_in(buf_size)[*])
    allocate( this%buf_northwest_in(buf_size)[*])
    allocate( this%buf_southeast_in(buf_size)[*])
    allocate( this%buf_southwest_in(buf_size)[*])

    if (this_image() .eq. 1) then
      print *, "===--- ARTLESS NEED TO DOUBLE CHECK NEIGHBORS ---==="
    end if
    ! set up the neighbors array so we can sync with our neighbors when needed
    if (.not.allocated(neighbors)) then
      associate(me=>this_image())
        north_neighbor = me + grid%ximages
        south_neighbor = me - grid%ximages
        east_neighbor  = me + 1
        west_neighbor  = me - 1
        northeast_neighbor = me + grid%ximages + 1
        northwest_neighbor = me + grid%ximages - 1
        southeast_neighbor = me - grid%ximages + 1
        southwest_neighbor = me - grid%ximages - 1

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
        print *, this_image(), ": boundry         ", this%north_boundary, &
            this%south_boundary, this%east_boundary, this%west_boundary
        print *, this_image(), ": diagonal boundry", this%northeast_boundary, &
            this%northwest_boundary, this%southeast_boundary, &
            this%southwest_boundary
      end associate
    endif
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

    if (.not. this%north_boundary) call this%retrieve_buf(this%buf_north_in)
    if (.not. this%south_boundary) call this%retrieve_buf(this%buf_south_in)
    if (.not. this%east_boundary) call this%retrieve_buf(this%buf_east_in)
    if (.not. this%west_boundary) call this%retrieve_buf(this%buf_west_in)
    if (.not. this%northeast_boundary) &
        call this%retrieve_buf(this%buf_northeast_in)
    if (.not. this%northwest_boundary) &
        call this%retrieve_buf(this%buf_northwest_in)
    if (.not. this%southeast_boundary) &
        call this%retrieve_buf(this%buf_southeast_in)
    if (.not. this%southwest_boundary) &
        call this%retrieve_buf(this%buf_southwest_in)

    this%north_i = 1
    this%south_i = 1
    this%east_i  = 1
    this%west_i  = 1
    this%northeast_i = 1
    this%northwest_i = 1
    this%southeast_i = 1
    this%southwest_i = 1
  end subroutine

  module subroutine exchange(this)
    class(convection_exchangeable_t), intent(inout) :: this
    ! if (.not. this%north_boundary) call this%put_north
    ! if (.not. this%south_boundary) call this%put_south
    ! if (.not. this%east_boundary)  call this%put_east
    ! if (.not. this%west_boundary)  call this%put_west

    sync images( neighbors )

    if (.not. this%north_boundary) call this%retrieve_buf(this%buf_north_in)
    if (.not. this%south_boundary) call this%retrieve_buf(this%buf_south_in)
    if (.not. this%east_boundary) call this%retrieve_buf(this%buf_east_in)
    if (.not. this%west_boundary) call this%retrieve_buf(this%buf_west_in)
    if (.not. this%northeast_boundary) &
        call this%retrieve_buf(this%buf_northeast_in)
    if (.not. this%northwest_boundary) &
        call this%retrieve_buf(this%buf_northwest_in)
    if (.not. this%southeast_boundary) &
        call this%retrieve_buf(this%buf_southeast_in)
    if (.not. this%southwest_boundary) &
        call this%retrieve_buf(this%buf_southwest_in)

    this%north_i = 1
    this%south_i = 1
    this%east_i  = 1
    this%west_i  = 1
    this%northeast_i = 1
    this%northwest_i = 1
    this%southeast_i = 1
    this%southwest_i = 1
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

  module subroutine retrieve_buf(this, buf)
    implicit none
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: buf(:)[*]
    integer :: i, n, local_i, local_n
    n = ubound(buf, dim=1)
    local_n = ubound(this%local, dim=1)
    do i=1,n
      associate (particle=>buf(i))
        if (particle%exists .eqv. .true.) then
          do local_i=1,local_n
            if (this%local(local_i)%exists .eqv. .false.) then
              call particle%move_to(this%local(local_i))
              exit
            end if
          end do
        end if
      end associate
    end do
  end subroutine

  module subroutine put_north(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%north_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_south_in(this%south_i)[north_neighbor] = particle
      particle%exists=.false.
      this%south_i = this%south_i + 1
  end subroutine

  module subroutine put_south(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%south_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_north_in(this%north_i)[south_neighbor] = particle
      particle%exists=.false.
      this%north_i = this%north_i + 1
  end subroutine

  module subroutine put_east(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%east_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_west_in(this%west_i)[east_neighbor] = particle
      particle%exists=.false.
      this%west_i = this%west_i + 1
  end subroutine

  module subroutine put_west(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%west_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_east_in(this%east_i)[west_neighbor] = particle
      particle%exists=.false.
      this%east_i = this%east_i + 1
  end subroutine

  module subroutine put_northeast(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%northeast_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_southwest_in(this%southwest_i)[northeast_neighbor] = particle
      particle%exists=.false.
      this%southwest_i = this%southwest_i + 1
  end subroutine

  module subroutine put_northwest(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%northwest_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_southeast_in(this%southeast_i)[northwest_neighbor] = particle
      particle%exists=.false.
      this%southeast_i = this%southeast_i + 1
  end subroutine

  module subroutine put_southeast(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%southeast_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_northwest_in(this%northwest_i)[southeast_neighbor] = particle
      particle%exists=.false.
      this%northwest_i = this%northwest_i + 1
  end subroutine

  module subroutine put_southwest(this, particle)
      class(convection_exchangeable_t), intent(inout) :: this
      type(convection_particle), intent(inout) :: particle

      if (this%southwest_boundary) then
        particle%exists = .false.
        return
      end if

      this%buf_northeast_in(this%northeast_i)[southwest_neighbor] = particle
      particle%exists=.false.
      this%northeast_i = this%northeast_i + 1
  end subroutine


  module subroutine process(this, dt, its,ite, jts,jte, kts,kte, &
        temperature)
      implicit none
      class(convection_exchangeable_t), intent(inout) :: this
      real,           intent(in)    :: dt
      integer,        intent(in)    :: its,ite, jts,jte, kts,kte
      real, dimension(:,:,:), intent(in) :: temperature
      real, parameter :: gravity = 9.80665
      real :: a_prime, displacement, t, t_prime
      real :: ws
      integer :: i,j,k, l_bound(1), dif(1), new_ijk(3), me

      me = this_image()
      do i=1,ubound(this%local,1)
        associate (particle=>this%local(i))
          if (particle%exists .eqv. .true.) then
            ! Check if particle is out of bounds
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
            ! Buoyancy (B), Temperature (T), Surrounding Tmperature (T')
            ! B = (T - T') / T'
            ! acceleration'_z = B * g                                    (11)
            !-----------------------------------------------------------------
            ! acceleration in terms of potential temperature             (13)
            ! theta = T(p_0 / p') ^ (R_d / c_p)
            ! acceleration'_z = B * g = (theta' - theta) / theta * gravity
            !-----------------------------------------------------------------
            ! displacement (s), velocity (u)
            ! s = u*t + 1/2 a*t^2
            !-----------------------------------------------------------------
            T = temperature(floor(particle%x), floor(particle%z), &
                floor(particle%y))

            T_prime = particle%temperature
            a_prime = (T_prime - T) / T * gravity

            ! time step is 1 so t and t^2 go away
            displacement = particle%velocity + 0.5 * a_prime
            particle%z = particle%z + displacement
            particle%velocity = displacement

            ! print *,me, "z = ", particle%z  , "with displacement",displacement

            if (particle%z .gt. kte) then
              particle%exists=.false.
            else if (particle%z .lt. kts) then ! kts will be 1
              particle%exists=.false.
            end if

            !-----------------------------------------------------------------
            ! Handle Windfield
            !-----------------------------------------------------------------
            ! ARTLESS, this is fake math right now
            ws = sqrt(particle%u * particle%u + particle%v * particle%v)
            ! u: zonal velocity, wind towards the east
            particle%x = particle%x + particle%u
            ! v: meridional velocity, wind towards north
            particle%y = particle%y + particle%v
            ! w: tangential velocity: ARTLESS

            if (particle%y .lt. jts) then      ! "jts  <   y    <  jte"
              if (particle%x .lt. its) then
                call this%put_southwest(particle)
              else if (particle%x .gt. ite) then
                call this%put_southeast(particle)
              else
                call this%put_south(particle)
              end if
            else if (particle%y .gt. jte) then ! jts will be 1
              if (particle%x .lt. its) then
                call this%put_northwest(particle)
              else if (particle%x .gt. ite) then
                call this%put_northeast(particle)
              else
                call this%put_north(particle)
              endif
            else if (particle%x .lt. its) then ! "its  <   x    <  ite"
              call this%put_west(particle)      ! need to double check this!
            else if (particle%x .gt. ite) then
              call this%put_east(particle)
            end if
          end if
        end associate
      end do

      ! Only needed if array of particles
      ! after moving everything need to reset movement flags to false
      ! do concurrent (i=its:ite, j=jts:jte, k=kts:kte)
      !   this%local(i,k,j)%moved = .false.
      ! end do
    end subroutine

end submodule
