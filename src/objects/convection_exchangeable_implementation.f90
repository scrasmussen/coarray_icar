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
  module subroutine const(this, convection_type_enum, grid,ims,ime,kms,kme,jms,&
      jme, input_buf_size, &
      halo_width, u_in, v_in, w_in, temperature, pressure, water_vapor)
    use iso_c_binding, only: c_int
    class(convection_exchangeable_t), intent(inout) :: this
    type(grid_t) :: grid
    integer, intent(in) :: ims,ime,kms,kme,jms,jme
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width
    integer(c_int), intent(in) :: convection_type_enum
    real, optional, intent(in) :: u_in,v_in,w_in
    real, dimension(:,:,:), intent(in) :: temperature, pressure, water_vapor
    real :: random_start(3)

    logical :: wrap_neighbors
    integer :: n_images
    integer :: n_neighbors, current, id_range, particle_id, i, j, k, me
    integer :: create, num_create
    real :: x,y,z,fx,fy,fz,t0,wv0
    real :: u,v,w
    me = this_image()
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
      ! ARTLESS:: removing ims, ime, etc, trying to get proper range
      ! ims = grid%ims - halo_east
      ! ime = grid%ime + halo_west
      ! jms = grid%jms - halo_south
      ! jme = grid%jme + halo_north
      ! kms = grid%kms
      ! kme = grid%kme


      ! ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
      ! ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
      ! its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

      ! print *, "----------------------------------"
      ! print *, "me = ", this_image()
      ! ! model domain
      ! print *, "ids,ide, jds,jde, kds,kde"
      ! print *, grid%ids,grid%ide, grid%jds,grid%jde, grid%kds,grid%kde
      ! ! memory in arrays
      ! print *, "ims,ime, jms,jme, kms,kme"
      ! print *, grid%ims,grid%ime, grid%jms,grid%jme, grid%kms,grid%kme
      ! ! data tile to process
      ! print *, "its,ite, jts,jte, kts,kte"
      ! print *, grid%its,grid%ite, grid%jts,grid%jte, grid%kts,grid%kte
      ! print *, "----------------------------------"


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

      num_create = 4
      if (me == 1) print*, "Creating", num_create, "parcels per image"

      do create=1,num_create
        call random_number(random_start)
        x = (ims) + (random_start(1) * (ime-ims))
        z = kms   + (random_start(3) * (kme-kms))
        y = (jms) + (random_start(2) * (jme-jms))
        ! x=4; y=2; z=3
        ! print *, "BOUNDS::", ims, ime, kms, kme, jms, jme, "CREATED::", x,z,y, " ON IMAGE", this_image()
        call this%create_particle_id()
        ! print *, me, ": particle_id = ", this%particle_id_count
        fx = floor(x)
        fy = floor(y)
        fz = floor(z)
        t0 = temperature(fx,fz,fy)
        wv0 = water_vapor(fx,fz,fy)

        this%local(create) = convection_particle( &
            this%particle_id_count, x,y,z, 0.5,0.5,0.0,&
            pressure(fx,fz,fy), &
            t0, wv0, e_sat_mr(t0,wv0))
      end do

      ! call random_number(random_start)
      ! if (this_image() .eq. 1) then
      !   print*,"ARTLESS change this to percentage once exchange working"
      !   ! pick random place in image to start
      !   call random_number(random_start)
      !   x = (ims) + (random_start(1) * (ime-ims))
      !   z = kms   + (random_start(3) * (kme-kms))
      !   y = (jms) + (random_start(2) * (jme-jms))
      !   ! x = (ims+halo_east) + (random_start(1) * (ime-ims-halo_west))
      !   ! z = kms + (random_start(3) * (kme-kms))
      !   ! y = (jms+halo_south) + (random_start(2) * (jme-jms-halo_north))
      !   print *, "BOUNDS::", ims, ime, kms, kme, jms, jme
      !   ! print *, "HALO::", halo_north ,halo_south ,halo_east ,halo_west
      !   print *, "CREATED::", x,z,y, " ON IMAGE", this_image()
      !   call this%create_particle_id()
      !   this%local(1) = convection_particle( &
      !       this%particle_id_count, x,y,z, 0.5,0.5,0.0,&
      !       pressure(floor(x),floor(z),floor(y)), &
      !       temperature(floor(x),floor(z),floor(y)))
      !   ! need to set exists to false
      ! end if

      ! if (this_image() .eq. 3) then
      !   print*,"ARTLESS change this to percentage once exchange working"
      !   call random_number(random_start)
      !   x = ims + (random_start(1) * (ime-ims))
      !   z = kms + (random_start(3) * (kme-kms))
      !   y = jms + (random_start(2) * (jme-jms))
      !   print *, "BOUNDS::", ims, ime, kms, kme, jms, jme
      !   ! print *, "HALO::", halo_north ,halo_south ,halo_east ,halo_west
      !   print *, "CREATED::", x,z,y, " ON IMAGE", this_image()

      !   call this%create_particle_id()
      !   this%local(1) = convection_particle( &
      !       this%particle_id_count, x,y,z, 0.5,0.5,0.0,&
      !       pressure(floor(x),floor(y),floor(z)), &
      !       temperature(floor(x),floor(y),floor(z)))
      ! ! need to set exists to false
      !   ! print *, "-----", ims, ime, kms, kme, jms, jme
      !   ! print *, "----- shape = ", shape(this%local)
      ! end if

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

        ! Wrap boundaries: so particles in the xy direction don't disappear
        wrap_neighbors = .true.
        if (wrap_neighbors .eqv. .true.) then
          n_images = num_images()
          if (this%north_boundary .eqv. .true.) then
            north_neighbor = north_neighbor - n_images
            this%north_boundary = .false.
            this%wrapped_north = .true.
          end if
          if (this%south_boundary .eqv. .true.) then
            south_neighbor = south_neighbor + n_images
            this%south_boundary = .false.
            this%wrapped_south = .true.
          end if
          if (this%east_boundary .eqv. .true.) then
            east_neighbor = east_neighbor - grid%ximages
            this%east_boundary = .false.
            this%wrapped_east = .true.
          end if
          if (this%west_boundary .eqv. .true.) then
            west_neighbor = west_neighbor + grid%ximages
            this%west_boundary = .false.
            this%wrapped_west = .true.
          end if
        end if

        ! if (this_image() .eq. 1) then
        !   print *, "grid x y ", grid%ximages, grid%yimages
        ! end if
        ! call flush()
        ! sync all
        ! print *, this_image(), ": north", north_neighbor, ": east", &
        !      east_neighbor, ": south", south_neighbor, &
        !      ": west", west_neighbor
        ! call flush()
        ! sync all
        ! if (this_image() .eq. 1) then
        !   print *, "                               n s e w"
        ! end if
        ! call flush()
        ! sync all
        ! print *, this_image(), ": boundry         ", this%north_boundary, &
        !     this%south_boundary, this%east_boundary, this%west_boundary
        ! print *, this_image(), ": diagonal boundry", this%northeast_boundary, &
        !     this%northwest_boundary, this%southeast_boundary, &
        !     this%southwest_boundary
      end associate
    endif

    sync all
    call exit
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
        temperature, dz)
      implicit none
      class(convection_exchangeable_t), intent(inout) :: this
      real,           intent(in)    :: dt, dz
      integer,        intent(in)    :: its,ite, jts,jte, kts,kte
      real, dimension(:,:,:), intent(in) :: temperature
      real, parameter :: gravity = 9.80665, Gamma = 0.01
      real :: a_prime, z_displacement, t, t_prime, buoyancy
      real :: ws, fake_wind_correction, delta_z, rcp
      integer :: i,j,k, l_bound(1), dif(1), new_ijk(3), me
      real :: new_pressure, R_s, p0, exponent, alt_pressure, alt_pressure2
      real :: alt_pressure3, alt_pressure4, mixing_ratio, sat_mr
      real :: vapor_p, sat_p, T_C, mr

      me = this_image()

      do i=1,ubound(this%local,1)
        associate (particle=>this%local(i))
          if (particle%exists .eqv. .true.) then
            ! Check if particle is out of bounds
            if (floor(particle%x) .lt. its .or. &
                floor(particle%x) .gt. ite .or. &
                floor(particle%z) .lt. kts .or. &
                floor(particle%z) .gt. kte .or. &
                floor(particle%y) .lt. jts .or. &
                floor(particle%y) .gt. jte) then
              print *, "particle", particle%particle_id, "on image", me
              print *, "x:", its, "<", particle%x, "<", ite
              print *, "z:", kts, "<", particle%z, "<", kte
              print *, "y:", jts, "<", particle%y, "<", jte
              stop "x,y,z is out of bounds"
            end if

            !-----------------------------------------------------------------
            ! Handle Buoyancy
            !-----------------------------------------------------------------
            ! Buoyancy (B), Temperature (T), Surrounding Tmperature (T')
            ! B = (T - T') / T'                                         (4.12)
            ! acceleration'_z = B * g                                    (11)
            !-----------------------------------------------------------------
            ! acceleration in terms of potential temperature             (13)
            ! theta = T(p_0 / p') ^ (R_d / c_p)
            ! acceleration'_z = B * g = (theta' - theta) / theta * gravity
            !-----------------------------------------------------------------
            ! displacement (s), initial velocity (u), final velocity (v)
            ! v = u + a * t
            ! s = u*t + 1/2 a*t^2
            !-----------------------------------------------------------------
            ! orig: properites of air parcel are prime
            ! new: properites of environment are prime
            T = particle%temperature
            T_prime = temperature(floor(particle%x), floor(particle%z), &
                floor(particle%y))

            ! 4.12 in Rogers and Yao, simple for now
            buoyancy = (T - T_prime) / T_prime
            a_prime = buoyancy * gravity

            ! time step is 1 so t and t^2 go away
            z_displacement = particle%velocity + 0.5 * a_prime
            ! number from dz_interface, currently always 500
            delta_z = z_displacement / dz
            particle%z = particle%z + delta_z
            particle%velocity = z_displacement



            ! print *, ":::::ARTLESS:::::"
            ! print *, "z_displacement = ", z_displacement
            ! print *, "temp     =" ,particle%temperature
            ! print *, "potential temp     =" ,particle%potential_temp
            ! print *, "pressure =" ,particle%pressure

            rcp = 0.286

            ! OLD
            !-----------------------------------------------------------------
            ! ! barometric formula for an adiabatic atmosphere
            ! ! p(h) = p_0 * (1 + Gamma / T_0 * h) ^ -g/R_s*Gamma
            ! new_pressure = particle%pressure * (1 + Gamma / )
            !-----------------------------------------------------------------
            ! dp = - g / (R_s * T) * p * dh
            ! R_s is specific gas constant 287


            !-----------------------------------------------------------------
            ! p = p_0 * e^( ((9.81/287.058)*dz) / t_mean )
            !-----------------------------------------------------------------
            p0 = particle%pressure
            particle%pressure = p0 - z_displacement * &
                gravity / (287.05 * particle%temperature) * p0


            ! ! barometric formula for an adiabatic atmosphere
            ! exponent = gravity / (287.05 * Gamma)
            ! alt_pressure = p0 * (1 - Gamma * particle%z / T)

            ! ! Exner
            ! ! T/potential =
            ! associate(po=>100000, Rd=>287.058, cp=>1003.5)
            !   alt_pressure2= p0 / ((particle%potential_temp / T)**(cp/Rd))
            ! end associate



            !-----------------------------------------------------------------
            ! Update the parcel temperature for dry air
            ! Gamma is the dry adiabatic lapse rate
            ! Temperature is reduced T - Gamma * delta_z  (3.8)
            ! print *, this_image(), "z_displacement", z_displacement
            ! Gamma in Kelvin per meter

            particle%temperature = particle%temperature - Gamma * z_displacement



            ! this is close to the potential temp given by the Exner func
            ! particle%potential_temp = particle%temperature * &
            !     (p0 / particle%pressure) ** (0.286)
            ! print *, "new potential temp", particle%potential_temp

            ! potential_temp from exner
            associate(po=>100000, Rd=>287.058, cp=>1003.5)
              particle%potential_temp = particle%temperature / &
                  ((particle%pressure / p0) ** (Rd/cp))
            end associate

            ! Equaion Ethan gave me, same as pressure
            ! alt_pressure = p0 * exp( -((9.81/287.058)*z_displacement) / &
            !     particle%temperature)
            ! alt_pressure4 = (p0 * (-Gamma * z_displacement)) / (-0.286*T)




            !-----------------------------------------------------------------
            ! Better update
            ! T_0 + gamma U dt - Gamma_s U d t
            ! T_0 is initial temperature
            ! Gamma_s is the pseudoadiabaitc lapse rate
            ! gamme is the ambient lapse rate


            !-----------------------------------------------------------------
            ! WANTED: LCL, lifted condensation level:
            ! approximation: h_LCL = 125(T-T_d)
            ! T_d = dew-point Temperature

            !-----------------------------------------------------------------
            ! Time to figure outclear





            ! print *,me,"z = ",particle%z,"with z_displacement",z_displacement

            if (particle%z .gt. kte) then
              particle%exists=.false.
            else if (particle%z .lt. kts) then ! kts will be 1
              particle%exists=.false.
            end if

            !-----------------------------------------------------------------
            ! Handle Windfield
            !-----------------------------------------------------------------
            ! ! ARTLESS, this fake math right now
            ! ws = sqrt(particle%u * particle%u + particle%v * particle%v)
            ! ARTLESS below might be correctish math
            ! u: zonal velocity, wind towards the east
            ! dz: 500, supposes all dz_interface values are the same
            fake_wind_correction = (1.0 / dz) ! real correction
            fake_wind_correction = (1.0 / 4) ! ARTLESS
            particle%x = particle%x + (particle%u * fake_wind_correction)
            ! v: meridional velocity, wind towards north
            particle%y = particle%y + (particle%v * fake_wind_correction)
            ! w: tangential velocity: ARTLESS

            !-----------------------------------------------------------------
            ! Handle Saturated Mixing Ratio
            !
            ! saturate mixing ratio: max amount of water vapor parcel can hold
            !                        without condensation
            !-----------------------------------------------------------------
            ! water_vapor_p =
            sat_mr = e_sat_mr(particle%temperature, particle%pressure)


            ! ! Mixing ratio: r, or mr
            ! ! r = R_d/ R_v  *  e/(p-e)
            ! ! from
            ! ! http://snowball.millersville.edu/~adecaria/ESCI340/esci340-Lesson02-Thermodynamics.pdf
            ! !
            ! ! Specific gas constants for dry air and water vapor R_d, R_v
            ! ! R_d = 287.058, R_v = 461.5, R_d / R_v = 0.622
            ! ! p, total air pressure
            ! ! e, vapor pressure, vapor_[]

            ! ! Antoine equation to find vapor pressure
            ! if (T_C .lt. 100) then
            !   vapor_p = 10 ** (8.07131 - 1730.63 / (233.426 + T_C))
            ! else
            !   vapor_p = 10 ** (8.14019 - 1810.94 / (244.485 + T_C))
            ! end if

            ! mr =  0.622 * (vapor_p / particle%pressure - vapor_p)

            ! print *, "-------"
            ! print *, "vapor_p =", vapor_p
            ! print *, "pressure =", particle%pressure
            ! print *, "mr =", mr
            ! ! print *, "saturation vapor pressure =", sat_p
            ! print *, "sat_mr = ",sat_mr
            ! print *, "relative humidity =", mr / sat_mr * 100
            ! print *, "-------"

            ! ! T in celcius
            ! T_C = particle%temperature - 273.15


            ! sat_p = 611**((17.27*T_C) / (237.3+T_C))
            ! print *, "-------"
            ! print *, "vapor pressure =", vapor_p
            ! print *, "saturation vapor pressure =", sat_p
            ! print *, "relative humidity =", vapor_p / sat_p * 100
            ! print *, "sat_mr = ",sat_mr
            ! print *, "-------"
            ! print *, "sat_mr" = sat_mr

            ! particle%mixing_ratio = ! change this to relative humidity



            associate (x => floor(particle%x), y => floor(particle%y), &
                z => floor(particle%z))
              if (x .lt. its .or. x .gt. ite .or. &
                  z .lt. kts .or. z .gt. kte .or. &
                  y .lt. jts .or. y .gt. jte) then
                print *, "PUTTING", particle%x, particle%z,  particle%y, &
                    "FROM", this_image(), "REMINDER!!! FAKE WIND"
                print *, its,ite, kts, kte, jts, jte
                print *, "------"
              end if

              if (y .lt. jts) then      ! "jts  <   y    <  jte"
                if (x .lt. its) then
                  call this%put_southwest(particle)
                else if (x .gt. ite) then
                  call this%put_southeast(particle)
                else
                  call this%put_south(particle)
                end if
              else if (y .gt. jte) then ! jts will be 1
                if (x .lt. its) then
                  call this%put_northwest(particle)
                else if (x .gt. ite) then
                  call this%put_northeast(particle)
                else
                  call this%put_north(particle)
                endif
              else if (x .lt. its) then ! "its  <   x    <  ite"
                call this%put_west(particle)      ! need to double check this!
              else if (x .gt. ite) then
                call this%put_east(particle)
              end if
            end associate

        end if
        end associate
      end do

      ! Only needed if array of particles
      ! after moving everything need to reset movement flags to false
      ! do concurrent (i=its:ite, j=jts:jte, k=kts:kte)
      !   this%local(i,k,j)%moved = .false.
      ! end do
    end subroutine

    module subroutine create_particle_id(this)
      use iso_fortran_env, only : int32
      implicit none
      class(convection_exchangeable_t), intent(inout) :: this
      integer :: id_range, h
      if (this%particle_id_count .eq. -1) then
        id_range = huge(int32) / num_images()
        this%particle_id_count = (this_image()-1) * id_range
      else
        this%particle_id_count = this%particle_id_count + 1
      end if
    end subroutine


    ! THIS IS COPIED FROM DOMAIN_IMPLEMENTATION.F90
    !>----------------------------------------------------------
    !!  Calculate the saturated mixing ratio for a given temperature and pressure
    !!
    !!  If temperature > 0C: returns the saturated mixing ratio with respect to liquid
    !!  If temperature < 0C: returns the saturated mixing ratio with respect to ice
    !!
    !!  @param temperature  Air Temperature [K]
    !!  @param pressure     Air Pressure [Pa]
    !!  @retval sat_mr      Saturated water vapor mixing ratio [kg/kg]
    !!
    !!  @see http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
    !!   Lowe, P.R. and J.M. Ficke., 1974: The Computation of Saturation Vapor Pressure
    !!   Environmental Prediction Research Facility, Technical Paper No. 4-74
    !!
    !!----------------------------------------------------------
    elemental function e_sat_mr(temperature,pressure)
    ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
        implicit none
        real,intent(in) :: temperature,pressure
        real :: e_s,a,b
        real :: e_sat_mr
        ! from http://www.dtic.mil/dtic/tr/fulltext/u2/778316.pdf
        !   Lowe, P.R. and J.M. Ficke., 1974: THE COMPUTATION OF SATURATION VAPOR PRESSURE
        !       Environmental Prediction Research Facility, Technical Paper No. 4-74
       ! which references:
        !   Murray, F. W., 1967: On the computation of saturation vapor pressure.
        !       Journal of Applied Meteorology, Vol. 6, pp. 203-204.
        ! Also notes a 6th order polynomial and look up table as viable options.
        if (temperature < 273.15) then
            a = 21.8745584
            b = 7.66
        else
            a = 17.2693882
            b = 35.86
        endif

        e_s = 610.78 * exp(a * (temperature - 273.16) / (temperature - b)) !(Pa)

        ! alternate formulations
        ! Polynomial:
        ! e_s = ao + t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t))))) a0-6 defined separately for water and ice
        ! e_s = 611.2*exp(17.67*(t-273.15)/(t-29.65)) ! (Pa)
        ! from : http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
        ! e_s = 611.0*10.0**(7.5*(t-273.15)/(t-35.45))


        if ((pressure - e_s) <= 0) then
            e_s = pressure * 0.99999
        endif
        ! from : http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf
        e_sat_mr = 0.6219907 * e_s / (pressure - e_s) !(kg/kg)
    end function e_sat_mr
end submodule
