submodule(convection_exchangeable_interface) &
    convection_exchangeable_implementation
  use assertions_interface, only : assert, assertions
  use exchangeable_interface, only : exchangeable_t
  use domain_interface, only : pressure_at_elevation, exner_function, sat_mr

  use grid_interface, only : grid_t
  implicit none

  integer, parameter :: default_buf_size=1
  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_neighbor, south_neighbor, buf_size, halo_size
  integer, save :: east_neighbor, west_neighbor
  integer, save :: northeast_neighbor, northwest_neighbor
  integer, save :: southeast_neighbor, southwest_neighbor
  logical, parameter :: wrap_neighbors = .true.
  logical, parameter :: advection = .false.
  logical, parameter :: wind = .true.
  logical, parameter :: caf_comm_message = .false.
  logical, parameter :: particle_create_message = .false.
  integer, parameter :: particles_per_image = 1
  integer, parameter :: local_buf_size = particles_per_image * 4


contains
  ! function initialize_convection_array_t(this)
  !   class(convection_exchangeable_array_t), intent(inout) :: this
  ! end function initialize_convection_array_t

  ! ----- STEPS -----
  ! CREATE PARTICLE
  ! input, z elevation, potential_temp,
  ! get z elevation
  !  -> pressure from pressure_at_elevation(sealevel_pressure, z_meters)
  !     VARY THE PRESSURE BY AN AMOUNT
  !  -> exner from exner_function(pressure)
  !  -> temp  from exner * potential_temp
  !  -> water_vapor = sat_mr(temp, pressure)

  ! Advection
  !  potential_temp, potential_temp` -> boyancy -> z displacement
  !    -> aka change in z
  !    -> change in temp,
  !       -> change in pressure
  !       -> change in potential temp ! NOT WHEN dry adiabatic
  !


  module subroutine const2(this, potential_temp, u_in, v_in, w_in, grid, z_m, &
      ims, ime, kms, kme, jms, jme, dz_value, &
      input_buf_size, halo_width)
    class(convection_exchangeable_t), intent(inout) :: this
    class(exchangeable_t), intent(in)    :: potential_temp
    class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
    type(grid_t), intent(in)      :: grid
    real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
    integer, intent(in)           :: ims, ime, kms, kme, jms, jme
    real, intent(in)              :: dz_value
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width

    integer :: me, create
    real :: random_start(3), x, z, y, north_adjust, east_adjust
    real :: z_meters, z_floor, z_ceiling
    real :: theta_val, theta_floor, theta_ceiling
    real :: pressure_val, exner_val, temp_val, water_vapor_val
    me = this_image()
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

    ! --- setup boundaries ---
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

    allocate(this%local(particles_per_image * 4))
    if (particle_create_message .eqv. .true.) then
      if (me == 1) then
        print*, "Creating", particles_per_image, "parcels per image"
      end if
    end if

    do create=1,particles_per_image
      call random_number(random_start)
      north_adjust = 0.99
      east_adjust = 0.99
      x = ims + (random_start(1) * (ime+east_adjust-ims))
      z = kms + (random_start(3) * (kme-kms))
      y = jms + (random_start(2) * (jme+north_adjust-jms))

      print *, x, z, y

      z_floor = z_m(floor(x),floor(z),floor(y))
      z_ceiling = z_m(ceiling(x),ceiling(z),ceiling(y))
      z_meters = z_floor + (z_ceiling - z_floor) * (modulo(z,1.0))

      print *, "z_meter =", z_meters
      ! sealevel_pressure => 100000.0
      pressure_val = pressure_at_elevation(100000.0, z_meters)
      exner_val = exner_function(pressure_val)

      theta_floor = potential_temp%local(floor(x),floor(z),floor(y))
      theta_ceiling = potential_temp%local(ceiling(x),ceiling(z),ceiling(y))
      theta_val = theta_floor + (theta_ceiling - theta_floor) * (modulo(z,1.0))

      temp_val = exner_val * theta_val
      water_vapor_val = sat_mr(temp_val, pressure_val)

      call this%create_particle_id()

      this%local(create) = convection_particle(x, z, y, dz_value, u_in, v_in, &
          w_in, z_meters, theta_val, water_vapor_val, this%particle_id_count)
      ! -- old --
      ! this%local(create) = convection_particle( &
      !     this%particle_id_count, x,y,z, 0.5,0.5,0.0,&
      !     pressure_val, t0, wv0, wv0 / e_sat_mr(t0,wv0))
    end do


    if (particle_create_message .eqv. .true.) then
      print *, "ALLOCATING BUFFERS OF SIZE", buf_size
    end if
    allocate( this%buf_north_in(buf_size)[*])
    allocate( this%buf_south_in(buf_size)[*])
    allocate( this%buf_east_in(buf_size)[*])
    allocate( this%buf_west_in(buf_size)[*])
    allocate( this%buf_northeast_in(buf_size)[*])
    allocate( this%buf_northwest_in(buf_size)[*])
    allocate( this%buf_southeast_in(buf_size)[*])
    allocate( this%buf_southwest_in(buf_size)[*])


    ! input, z elevation, potential_temp :: potential is exchangable, z is real
    print *, "--- fin - ish ---"
  end subroutine

  ! constructor
  module subroutine const(this, convection_type_enum, grid,tims,time,tkms,tkme,tjms,&
      tjme, input_buf_size, &
      halo_width, u_in, v_in, w_in, temperature, pressure, water_vapor)
    use iso_c_binding, only: c_int
    class(convection_exchangeable_t), intent(inout) :: this
    type(grid_t) :: grid
    integer, intent(in) :: tims,time,tkms,tkme,tjms,tjme
    integer :: ims,ime,kms,kme,jms,jme
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width
    integer(c_int), intent(in) :: convection_type_enum
    real, optional, intent(in) :: u_in,v_in,w_in
    real, dimension(:,:,:), intent(in) :: temperature, pressure, water_vapor
    real :: random_start(3)

    integer :: n_images
    integer :: n_neighbors, current, id_range, particle_id, i, j, k, me
    integer :: create
    real :: north_adjust, east_adjust
    real :: x,y,z,t0,wv0
    real :: u,v,w,pressure_val
    integer :: fx,fy,fz
    logical, parameter :: broken=.true.
    integer :: u_bound(3), broken_fx, broken_fy, broken_fz

    me = this_image()

    if (particle_create_message .eqv. .true.) then
      sync all
      print *, me,":grid = ", grid
      sync all
    end if

    ims = grid%ims
    ime = grid%ime
    kms = grid%kms
    kme = grid%kme
    jms = grid%jms
    jme = grid%jme

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

    ! ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
    ! ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
    ! its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)
    ! ARTLESS:: removing ims, ime, etc, trying to get proper range

    allocate(this%local(local_buf_size))

    if (particle_create_message .eqv. .true.) then
      if (me == 1) then
        print*, "Creating", particles_per_image, "parcels per image"
        print*, "buffer based on particles_per_image value not input_buf_size"
      end if
    end if

    do i=1, num_images()
      call flush()
      sync all
      print *, me, "---l=", lbound(temperature), "u=", ubound(temperature)
      call flush()
      sync all
    end do



    do create=1,particles_per_image
      call random_number(random_start)
      north_adjust = 0.99
      east_adjust = 0.99
      x = ims + (random_start(1) * (ime+east_adjust-ims))
      z = kms + (random_start(3) * (kme-kms))
      y = jms + (random_start(2) * (jme+north_adjust-jms))



      ! if (particle_create_message .eqv. .true. &
      !     .and. me .eq. 16 &
      !     ) then
      !   print *, "BOUNDS::", ims, ime, kms, kme, jms, jme, &
      !       "CREATED::", x,z,y, " ON IMAGE", this_image()
      !   print *, "---l=", lbound(temperature), "u=", ubound(temperature)
      ! end if

      call this%create_particle_id()
    ! print *, me, ": particle_id = ", this%particle_id_count

      ! fx = floor(x)
      ! fy = floor(y)
      ! fz = floor(z)
      fx = floor(x)
      fy = floor(y)
      fz = floor(z)
      t0 = temperature(fx,fz,fy)
      wv0 = water_vapor(fx,fz,fy)
      pressure_val = pressure(fx,fz,fy)

      ! doing this because bounds aren't being passed
      if (broken .eqv. .true.) then
        u_bound = ubound(temperature)
        broken_fx = mod((fx-1),u_bound(1)+1)
        broken_fz = mod((fz-1),u_bound(2)+1)
        broken_fy = mod((fy-1),u_bound(3)+1)
        t0 = temperature(broken_fx,broken_fz,broken_fy)
        wv0 = water_vapor(broken_fx,broken_fz,broken_fy)
        pressure_val = pressure(broken_fx,broken_fz,broken_fy)
        ! print *, me, "BROKEN:::::::::", broken_fx, broken_fz, broken_fy, t0, wv0, pressure_val
      end if

      this%local(create) = convection_particle( &
          this%particle_id_count, x,y,z, 0.5,0.5,0.0,&
          pressure_val, t0, wv0, wv0 / e_sat_mr(t0,wv0))


            ! sometimes we create a particle that will grab an out of bounds
      ! temperature, this fixes that
      ! do while (t0 .eq. 0)
      !   call random_number(random_start)
      !   x = ims + (random_start(1) * (ime+east_adjust-ims))
      !   z = kms + (random_start(3) * (kme-kms))
      !   y = jms + (random_start(2) * (jme+north_adjust-jms))
      !   fx = floor(x)
      !   fy = floor(y)
      !   fz = floor(z)
      !   t0 = temperature(fx,fz,fy)
      ! end do

    end do

    if (particle_create_message .eqv. .true.) then
      print *, "ALLOCATING BUFFERS OF SIZE", buf_size
    end if
    allocate( this%buf_north_in(buf_size)[*])
    allocate( this%buf_south_in(buf_size)[*])
    allocate( this%buf_east_in(buf_size)[*])
    allocate( this%buf_west_in(buf_size)[*])
    allocate( this%buf_northeast_in(buf_size)[*])
    allocate( this%buf_northwest_in(buf_size)[*])
    allocate( this%buf_southeast_in(buf_size)[*])
    allocate( this%buf_southwest_in(buf_size)[*])

    ! if (this_image() .eq. 1) then
    !   print *, "===--- ARTLESS NEED TO DOUBLE CHECK NEIGHBORS ---==="
    ! end if
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
        associate(nx => grid%ximages, ny => grid%yimages, &
            nimages => grid%ximages * grid%yimages)
        if (wrap_neighbors .eqv. .true.) then
          n_images = num_images()
          ! --- handle diagonals
          if (this%north_boundary .eqv. .true.) then
            northeast_neighbor = modulo( (me+nx+1)-nimages, (nx+1))
            if (northeast_neighbor .eq. 0) northeast_neighbor = 1
            northwest_neighbor = (me+nx-1)-nimages
            if (northwest_neighbor .eq. 0) northwest_neighbor = nx
          else  if (this%east_boundary .eqv. .true.) then
            northeast_neighbor = me + 1
          else  if (this%west_boundary .eqv. .true.) then
            northwest_neighbor = me + nx * 2 - 1
          end if

          if (this%south_boundary .eqv. .true.) then
            southeast_neighbor = me + nimages - nx + 1
            if (southeast_neighbor > nimages) then
              southeast_neighbor = southeast_neighbor - nx
            end if
            southwest_neighbor = me + nimages - nx - 1
            if (modulo(southwest_neighbor,nx) == 0) then
              southwest_neighbor = nimages
            end if
            ! southwest_neighbor = - 2
          else  if (this%east_boundary .eqv. .true.) then
            southeast_neighbor = me - nx * 2 + 1
          else  if (this%west_boundary .eqv. .true.) then
            southwest_neighbor = me - 1
            ! southwest_neighbor = - 1
          end if
          ! this%northeast_boundary = .false.
          ! this%northwest_boundary = .false.
          ! this%southeast_boundary = .false.
          ! this%southwest_boundary = .false.

          ! --- handle up/down/left/right
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
        end associate


        do i=1,num_images()
          call flush()
          sync all
          if (this_image() .eq. i) then
          print *, "                               n s e w"
          print *, this_image(), ": boundry         ", this%north_boundary, &
              this%south_boundary, this%east_boundary, this%west_boundary
          print *, "                               nw ne se sw"
          print *, this_image(), ": diagonal boundry", this%northwest_boundary, &
              this%northeast_boundary, this%southeast_boundary, &
              this%southwest_boundary
          end if
          sync all
        end do

        ! if (this_image() .eq. 1) then
        !   print *, "grid x y ", grid%ximages, grid%yimages
        ! end if
        ! call flush()
        ! sync all
        ! ! print *, this_image(), ": w", west_neighbor, ": n", &
        ! !      north_neighbor, ": e", east_neighbor, &
        ! !      ": s", south_neighbor
        ! call flush()
        ! sync all
        ! print *, this_image(), ": nw", northwest_neighbor, ": ne", &
        !      northeast_neighbor, ": se", southeast_neighbor, &
        !      ": sw", southwest_neighbor
        ! call flush()
        ! sync all
        ! ! if (this_image() .eq. 1) then
        ! !   ! print *, "                               n s e w"
        ! !   print *, "                               nw ne se sw"
        ! ! end if
        ! ! call flush()
        ! ! sync all
        ! ! ! print *, this_image(), ": boundry         ", this%north_boundary, &
        ! ! !     this%south_boundary, this%east_boundary, this%west_boundary
        ! ! print *, this_image(), ": diagonal boundry", this%northwest_boundary, &
        ! !     this%northeast_boundary, this%southeast_boundary, &
        ! !     this%southwest_boundary
        ! sync all
        ! call exit

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


  module subroutine process(this, nx_global, ny_global, grid, &
         dt, dz, temperature)
      implicit none
      class(convection_exchangeable_t), intent(inout) :: this
      type(grid_t), intent(in) :: grid
      integer, intent(in) :: nx_global, ny_global
      real, intent(in)    :: dt, dz
      real, dimension(:,:,:), intent(in) :: temperature
      real, parameter :: gravity = 9.80665, Gamma = 0.01
      integer :: ims,ime, jms,jme, kms,kme
      real :: a_prime, z_displacement, t, t_prime, buoyancy
      real :: ws, fake_wind_correction, delta_z
      integer :: i,j,k, l_bound(1), dif(1), new_ijk(3), me
      real :: new_pressure, R_s, p0, exponent, alt_pressure, alt_pressure2
      real :: alt_pressure3, alt_pressure4, mixing_ratio, sat_mr
      real :: vapor_p, sat_p, T_C, mr, T_squared, tmp
      logical, parameter :: broken=.true.
      integer :: u_bound(3), broken_fx, broken_fy, broken_fz
      real :: T_broken


      me = this_image()
      ims = grid%ims
      ime = grid%ime
      jms = grid%jms
      jme = grid%jme
      kms = grid%kms
      kme = grid%kme

      do i=1,ubound(this%local,1)
        associate (particle=>this%local(i))
          if (particle%exists .eqv. .true.) then
            ! Check if particle is out of bounds
            if (floor(particle%x) .lt. ims .or. &
                floor(particle%x) .gt. ime .or. &
                floor(particle%z) .lt. kms .or. &
                floor(particle%z) .gt. kme .or. &
                floor(particle%y) .lt. jms .or. &
                floor(particle%y) .gt. jme) then
              print *, "particle", particle%particle_id, "on image", me
              print *, "x:", ims, "<", particle%x, "<", ime
              print *, "z:", kms, "<", particle%z, "<", kme
              print *, "y:", jms, "<", particle%y, "<", jme
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

            if (broken .eqv. .true.) then
              u_bound = ubound(temperature)
              broken_fx = mod((floor(particle%x)-1),u_bound(1)+1)
              broken_fz = mod((floor(particle%z)-1),u_bound(2)+1)
              broken_fy = mod((floor(particle%y)-1),u_bound(3)+1)
              T_broken = temperature(broken_fx,broken_fz,broken_fy)
              ! print *, "BROKEN:::::::::", broken_fx, broken_fz, broken_fy, &
              !     T, T_broken, T_prime
              T_prime = T_broken
            end if



            ! print *, "------------ ubound---", ubound(temperature), T_prime

            ! if (me .eq. 2) then
            !   print *,  "ME ================", 2, "|", ubound(this%local,1)
            ! end if

            if (advection .eqv. .true.) then
            ! if (advection .eqv. .true. .and. me .eq. 2) then
              ! 4.12 in Rogers and Yao, simple for now
              ! print *, "id = ", particle%particle_id

              buoyancy = (T - T_prime) / T_prime
              a_prime = buoyancy * gravity
              ! time step is 1 so t and t^2 go away
              z_displacement = particle%velocity + 0.5 * a_prime
              ! number from dz_interface, currently always 500
              delta_z = z_displacement / dz
              ! print *, "bouyancy::", buoyancy, T,  T_prime
              ! print *, "dispalce::", z_displacement, particle%velocity, a_prime
              ! print *, "delta_z ::", delta_z, z_displacement, dz
              if (z_displacement /= z_displacement) then
                ! when T=0 it causes division by 0. This problem should be fixed
                print *, me, ":: ------------NAN--------------"
                particle%z = -1
              else
                particle%z = particle%z + delta_z
                particle%velocity = z_displacement
              end if
            else
              z_displacement = 0.0
            end if

            ! print *, ":::::ARTLESS:::::"
            ! print *, "z_displacement = ", z_displacement
            ! print *, "temp     =" ,particle%temperature
            ! print *, "potential temp     =" ,particle%potential_temp
            ! print *, "pressure =" ,particle%pressure

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
            if (particle%relative_humidity .ge. 1.01) then
              ! equation from
              ! wikipedia.org/wiki/Lapse_rate#Moist_adiabatic_lapse_rate
              T_C = particle%temperature - 273.15
              T_squared = T_C * T_C
              mixing_ratio = particle%pressure

              tmp = gravity * (286 * T_squared + &
                  2501000 * mixing_ratio * T_C ) / &
                  (1003.5 *  287 * T_squared + &
                  0.622 * mixing_ratio * 2501000 * 2501000)
              particle%temperature = tmp + 273.15
              ! particle%temperature = gravity * (286 * T_squared + &
              !     2501000 * mixing_ratio * T_C ) / &
              !     (1003.5 *  287 * T_squared + &
              !      0.622 * mixing_ratio * 2501000 * 2501000) + 273.15
              ! ---ARTLESS---
              print *, "MALR: from", T_C+273.15, "to", particle%temperature, &
                  "diff", T_C - tmp

              ! other possible equation
              ! http://www.theweatherprediction.com/habyhints/161/
              ! MALR = dT/dz = DALR / (1 + L/Cp*dWs/dT)
              ! dWs/dT is the change in saturation mixing ratio with change in T
            else
              particle%temperature = particle%temperature - Gamma * z_displacement
            end if



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
            ! remove particle if beyond the z axis
            !-----------------------------------------------------------------
            if (particle%z .gt. kme) then
              particle%exists=.false.
              ! print *, "particle", particle%particle_id, "gone out the top"
              cycle
            else if (particle%z .lt. kms) then ! kts will be 1
              particle%exists=.false.
              ! print *, "particle", particle%particle_id, "gone out the bottom"
              cycle
            end if

            !-----------------------------------------------------------------
            ! Handle Windfield
            !-----------------------------------------------------------------
            ! ! ARTLESS, this fake math right now
            ! ws = sqrt(particle%u * particle%u + particle%v * particle%v)
            ! ARTLESS below might be correctish math
            ! dz: 500, supposes all dz_interface values are the same
            fake_wind_correction = 0
            ! fake_wind_correction = (1.0 / dz) ! real correction
            ! fake_wind_correction = (1.0 / 10) ! ARTLESS
            fake_wind_correction = (1.0 / 4) ! ARTLESS
            fake_wind_correction = (1.0 / 4) ! ARTLESS
            fake_wind_correction = (1.0) ! ARTLESS

            if (wind .eqv. .true.) then
              ! u: zonal velocity, wind towards the east
              particle%x = particle%x + (particle%u * fake_wind_correction)
              ! v: meridional velocity, wind towards north
              particle%y = particle%y + (particle%v * fake_wind_correction)
            end if

            !-----------------------------------------------------------------
            ! Handle Saturated Mixing Ratio
            !
            ! saturate mixing ratio: max amount of water vapor parcel can hold
            !                        without condensation
            !-----------------------------------------------------------------
            ! water_vapor = water vapor mixing ratio (w)
            ! relative humidity = w / w_s

            sat_mr = e_sat_mr(particle%temperature, particle%pressure)
            particle%relative_humidity = particle%water_vapor / sat_mr

            ! print *, me, "RH = ", particle%relative_humidity
            ! if (particle%relative_humidity .ge. 1) then
            !   print *, "particle", particle%particle_id, &
            !       particle%relative_humidity, "=", &
            !       particle%water_vapor, "/", sat_mr
            ! end if

            ! ! Antoine equation to find vapor pressure
            ! if (T_C .lt. 100) then
            !   vapor_p = 10 ** (8.07131 - 1730.63 / (233.426 + T_C))
            ! else
            !   vapor_p = 10 ** (8.14019 - 1810.94 / (244.485 + T_C))
            ! end if
            ! mr =  0.622 * (vapor_p / particle%pressure - vapor_p)

            associate (x => floor(particle%x), y => floor(particle%y), &
                z => floor(particle%z))
              if (x .lt. ims .or. x .gt. ime .or. &
                  z .lt. kms .or. z .gt. kme .or. &
                  y .lt. jms .or. y .gt. jme) then
                if (caf_comm_message .eqv. .true.) then
                  print *, "PUTTING", particle%x, particle%z,  particle%y, &
                      "FROM", this_image(), "id:", particle%particle_id
                  print *, ims,ime, kms, kme, jms, jme
                  print *, "------"
                end if
              end if


              ! If particle is getting wrapped the x and y values need to be
              ! properly updated
              if (wrap_neighbors .eqv. .true.) then
                if (particle%x > nx_global + 1) then
                  if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                  particle%x = particle%x - nx_global
                else if (particle%x < 1) then
                  if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                  particle%x = particle%x + nx_global - 1
                end if

                if (particle%y > ny_global + 1) then
                  if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                  particle%y = particle%y - ny_global
                else if (particle%y < 1) then
                  if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                  particle%y = particle%y + ny_global - 1
                end if
              end if

              ! Check values to know where to send particle
              if (y .lt. jms) then      ! "jts  <   y    <  jte"
                if (x .lt. ims) then
                  call this%put_southwest(particle)
                else if (x .gt. ime) then
                  call this%put_southeast(particle)
                else
                  call this%put_south(particle)
                end if
              else if (y .gt. jme) then ! jts will be 1
                if (x .lt. ims) then
                  call this%put_northwest(particle)
                else if (x .gt. ime) then
                  call this%put_northeast(particle)
                else
                  call this%put_north(particle)
                endif
              else if (x .lt. ims) then ! "its  <   x    <  ite"
                call this%put_west(particle)      ! need to double check this!
              else if (x .gt. ime) then
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
