#ifdef __NVCOMPILER
#ifdef DC_LOOP
#define EXPAND_FUNC
#endif
#endif

#ifdef __NVCOMPILER
#ifdef DO_LOOP
#define EXPAND_FUNC
#endif
#endif

#ifdef __NVCOMPILER
#ifdef OMP_LOOP
#define EXPAND_FUNC
#endif
#endif

! #define EXPAND_FUNC
#define NO_RANDOM_NUMBER



submodule(convection_exchangeable_interface) &
    convection_exchangeable_implementation
  use assertions_interface, only : assert, assertions
  ! use exchangeable_interface, only : exchangeable_t
  use domain_interface, only : pressure_at_elevation, exner_function, sat_mr
  ! use cudafor, only : random_number
  ! use grid_interface, only : grid_t
  implicit none

  ! ----- PARAMETERS TO TUNE CONVECTION MODEL -----
  logical, parameter :: debug = .false.
  logical, parameter :: wrap_neighbors = .true.
  logical, parameter :: convection = .true.
  logical, parameter :: wind = .true.
  logical, parameter :: fake_wind_correction = .false.
  logical, parameter :: use_input_wind = .true.
  logical, parameter :: caf_comm_message = .false.
  logical, parameter :: particle_create_message = .false.
  logical, parameter :: brunt_vaisala_data = .false.
  logical, parameter :: replacement = .true.
  logical, parameter :: replacement_message = .false.
  logical, parameter :: init_theta = .false.
  logical, parameter :: init_velocity = .true.
  logical, parameter :: count_p_comm = .false.

#ifndef NO_COARRAYS
  integer, save      :: particles_communicated[*]
#else
  integer, save      :: particles_communicated
#endif

! #ifnotdef NO_COARRAYS
! #else
! #endif
  integer, save      :: particles_per_image
  integer, save      :: local_buf_size
  real,    save      :: input_wind
  logical, save      :: dry_air_particles
  ! integer, parameter :: particles_per_image=1
  ! integer, parameter :: local_buf_size=4*particles_per_image
  ! logical, parameter :: dry_air_particles=.true.
  ! -----------------------------------------------

  integer, parameter :: default_buf_size=1
  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_con_neighbor, south_con_neighbor, buf_size, halo_size
  integer, save :: east_con_neighbor, west_con_neighbor
  integer, save :: northeast_con_neighbor, northwest_con_neighbor
  integer, save :: southeast_con_neighbor, southwest_con_neighbor
  integer, save :: current_max_local_particles

contains
  module subroutine convect_const(this, potential_temp, u_in, v_in, w_in, grid,&
      z_m, z_interface, ims, ime, kms, kme, jms, jme, dz_val, &
      its, ite, kts, kte, jts, jte, pressure, input_buf_size, halo_width)
    class(convection_exchangeable_t), intent(inout) :: this
    class(exchangeable_t), intent(in)    :: potential_temp
    class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
    type(grid_t), intent(in)      :: grid
    real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
    real, intent(in)              :: pressure(ims:ime,kms:kme,jms:jme)
    real, intent(in)              :: z_interface(ims:ime,jms:jme)
    integer, intent(in)           :: ims, ime, kms, kme, jms, jme
    integer, intent(in)           :: its, ite, kts, kte, jts, jte
    real, intent(in)              :: dz_val
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width

    integer :: me, create, seed(34)
    real :: random_start(3), x, z, y
    real :: z_meters, z_interface_val, theta_val
    real :: pressure_val, exner_val, temp_val, water_vapor_val
    real :: u_val, v_val, w_val
    integer :: x0, x1, z0, z1, y0, y1
    logical :: calc

#if NO_COARRAYS
    me = 1
#else
    me = this_image()
#endif
    call initialize_from_file()
    particles_communicated = 0
    if (particles_per_image .eq. 0) then
       if (me .eq. 1) print *, "No air parcels used"
       return
    end if

    ! if (present(input_buf_size)) then
    !   buf_size = input_buf_size
    !   print *, 'using input_buf_size'
    ! else
    !   buf_size = default_buf_size
    !   print *, 'using default_buf_size'
    ! end if
    buf_size = ceiling(particles_per_image * 1.25) + 10  ! artless
    ! buf_size = particles_per_image * 2

    if (present(halo_width)) then
      halo_size = halo_width
    else
      halo_size = default_halo_size
    end if

    ! if (me .eq. 1) print *, "Buf size is", buf_size


    if (allocated(this%local)) deallocate(this%local)
    this%north_boundary = (grid%yimg == grid%yimages)
    this%south_boundary = (grid%yimg == 1)
    this%east_boundary  = (grid%ximg == grid%ximages)
    this%west_boundary  = (grid%ximg == 1)

    allocate(this%local(buf_size))
    if (particle_create_message .eqv. .true.) then
      if (me == 1) then
        print*, "Creating", particles_per_image, "parcels per image"
      end if
    end if

    seed = -1
    call random_seed(PUT=seed)
#ifndef __NVCOMPILER
    call random_init(.true.,.true.)
#endif

    current_max_local_particles = particles_per_image
    do create=1,particles_per_image
      call this%create_particle_id()
      this%local(create) = create_particle(this%particle_id_count, &
          its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, &
          z_m, potential_temp, z_interface, pressure, u_in, v_in, w_in)

      ! call exit ! artless
      ! this%local(create) = convection_particle(this%particle_id_count, .true., &
      !     .false., x, y, z, u_val, v_val, w_val, z_meters, z_interface_val, &
      !     pressure_val, temp_val, theta_val, 0, water_vapor_val, 0)
    end do


    if (particle_create_message .eqv. .true.) then
      print *, "ALLOCATING BUFFERS OF SIZE", buf_size
    end if
#if NO_COARRAYS
    allocate( this%buf_north_in(buf_size))
    allocate( this%buf_south_in(buf_size))
    allocate( this%buf_east_in(buf_size))
    allocate( this%buf_west_in(buf_size))
    allocate( this%buf_northeast_in(buf_size))
    allocate( this%buf_northwest_in(buf_size))
    allocate( this%buf_southeast_in(buf_size))
    allocate( this%buf_southwest_in(buf_size))
#else
    allocate( this%buf_north_in(buf_size)[*])
    allocate( this%buf_south_in(buf_size)[*])
    allocate( this%buf_east_in(buf_size)[*])
    allocate( this%buf_west_in(buf_size)[*])
    allocate( this%buf_northeast_in(buf_size)[*])
    allocate( this%buf_northwest_in(buf_size)[*])
    allocate( this%buf_southeast_in(buf_size)[*])
    allocate( this%buf_southwest_in(buf_size)[*])
#endif

    call this%setup_neighbors(grid)
  end subroutine convect_const

  module function create_particle(particle_id, its, ite, kts, kte, jts, jte, &
      ims, ime, kms, kme, jms, jme, z_m, potential_temp, z_interface, &
      pressure, u_in, v_in, w_in, times_moved) result(particle)
    integer :: particle_id
    type(convection_particle) :: particle
    integer, intent(in)           :: ims, ime, kms, kme, jms, jme
    integer, intent(in)           :: its, ite, kts, kte, jts, jte
    real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
    real, intent(in)              :: pressure(ims:ime,kms:kme,jms:jme)
    class(exchangeable_t), intent(in)    :: potential_temp
    real, intent(in)              :: z_interface(ims:ime,jms:jme)
    class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
    integer, intent(in), optional :: times_moved
    real :: relative_humidity_in
    real :: z_meters, z_interface_val, theta_val
    real :: random_start(3), x, z, y
    integer :: x0, x1, z0, z1, y0, y1, times_moved_val
    real :: pressure_val, exner_val, temp_val, water_vapor_val
    real :: u_val, v_val, w_val, velocity, cloud_water, r1

#ifdef NO_RANDOM_NUMBER
    x = its + ((1.0/(particle_id+2.0)) * (ite-its))
    z = kts + ((2.0/(particle_id+3.0)) * (kte-kts))
    y = jts + ((3.0/(particle_id+4.0)) * (jte-jts))
    ! x = curandCreateGenerator() 100
    ! call random_number(random_start)
    ! print *, random_start
    ! call random_number(r1)
    ! print*, r1
#else
    call random_number(random_start)
    x = its + (random_start(1) * (ite-its))
    z = kts + (random_start(3) * (kte-kts))
    y = jts + (random_start(2) * (jte-jts))
#endif
    ! print *, x,y,z

    if (x .lt. its .or. &
        x .gt. ite .or. &
        z .lt. kts .or. &
        z .gt. kte .or. &
        y .lt. jts .or. &
        y .gt. jte) then
      print *, "x:", its, "<", x, "<", ite
      print *, "z:", kts, "<", z, "<", kte
      print *, "y:", jts, "<", y, "<", jte
      stop "x,y,z is out of bounds"
    end if

    x0 = floor(x); z0 = floor(z); y0 = floor(y);
    x1 = ceiling(x); z1 = ceiling(z); y1 = ceiling(y);

    associate (A => z_m)
      z_meters = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
          A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
          A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
    end associate

    associate (A => potential_temp%local)
      theta_val = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
          A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
          A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
    end associate

    associate (A => pressure)
      pressure_val = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
          A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
          A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
    end associate

    ! print *, "val and func", pressure_val, pressure_at_elevation(100000.0, z_meters)
    ! call exit

    associate (A => z_interface)
      z_interface_val = bilinear_interpolation(x, x0, x1, y, y0, y1, &
          A(x0,y0), A(x0,y1), A(x1,y0), A(x1,y1))
    end associate

    ! print *, "z meters =" , z_meters, "z =", z, "zm =", z_m(x0,z0,y0)
    ! print *, "z_interface_val =", z_interface_val, "try =", z*dz_val + z_interface_val -250
    ! print *, "----"

    ! call flush()
    ! sync all
    ! call exit

    ! sealevel_pressure => 100000.0
    ! pressure_val = pressure_at_elevation(100000.0, z_meters)
    exner_val = exner_function(pressure_val)

    ! call random_number(rand)
    if (init_theta .eqv. .true.) then
      theta_val = theta_val * (1 + 1.0 / 100) ! random 0-1% change
    else
      theta_val = theta_val ! * (1 + 0.01) ! increase by 1%
    end if

    if (init_velocity .eqv. .true.) then
      velocity = 5
    else
      velocity = 0
    end if


    temp_val = exner_val * theta_val

    if (dry_air_particles .eqv. .true.) then
       water_vapor_val = 0
    else
       water_vapor_val = sat_mr(temp_val, pressure_val) *  1.0 !0.99
    end if
    relative_humidity_in = water_vapor_val

    ! ARTLESS: Testing, set to 0
    ! relative_humidity_in = 0

    ! wind is constant in this system, ignoring w wind (aka z-direction)
    u_val = u_in%local(x0, z0, y0)
    v_val = v_in%local(x0, z0, y0)
    w_val = 0


    cloud_water = 0
    if (present(times_moved) .eqv. .true.) then
       times_moved_val = times_moved
    else
       times_moved_val = 0
    end if


    particle = convection_particle(particle_id, .true., times_moved_val, &
        x, y, z, u_val, v_val, w_val, z_meters, z_interface_val, &
        pressure_val, temp_val, theta_val, velocity, water_vapor_val, &
        cloud_water, relative_humidity_in)

  end function

  module subroutine replace_particle(particle_id, its, ite, kts, kte, jts, jte, &
      ims, ime, kms, kme, jms, jme, z_m, potential_temp, z_interface, &
      pressure, u_in, v_in, w_in, times_moved, particle)
    type(convection_particle), intent(inout) :: particle
    integer, intent(in) :: particle_id
    integer, intent(in)           :: ims, ime, kms, kme, jms, jme
    integer, intent(in)           :: its, ite, kts, kte, jts, jte
    real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
    real, intent(in)              :: pressure(ims:ime,kms:kme,jms:jme)
    class(exchangeable_t), intent(in)    :: potential_temp
    real, intent(in)              :: z_interface(ims:ime,jms:jme)
    class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
    integer, intent(in), optional :: times_moved
    real :: relative_humidity_in
    real :: z_meters, z_interface_val, theta_val
    real :: random_start(3), x, z, y
    integer :: x0, x1, z0, z1, y0, y1, times_moved_val
    real :: pressure_val, exner_val, temp_val, water_vapor_val
    real :: u_val, v_val, w_val, velocity, cloud_water

    call random_number(random_start)

    x = its + (random_start(1) * (ite-its))
    z = kts + (random_start(3) * (kte-kts))
    y = jts + (random_start(2) * (jte-jts))

    ! if (x .lt. its .or. &
    !     x .gt. ite .or. &
    !     z .lt. kts .or. &
    !     z .gt. kte .or. &
    !     y .lt. jts .or. &
    !     y .gt. jte) then
    !   ! print *, "x:", its, "<", x, "<", ite
    !   ! print *, "z:", kts, "<", z, "<", kte
    !   ! print *, "y:", jts, "<", y, "<", jte
    !    ! print*, "x,y,z is out of bounds"
    !   return
    ! end if

    x0 = floor(x); z0 = floor(z); y0 = floor(y);
    x1 = ceiling(x); z1 = ceiling(z); y1 = ceiling(y);

    associate (A => z_m)
      z_meters = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
          A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
          A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
    end associate

    associate (A => potential_temp%local)
      theta_val = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
          A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
          A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
    end associate

    associate (A => pressure)
      pressure_val = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
          A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
          A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
    end associate

    ! print *, "val and func", pressure_val, pressure_at_elevation(100000.0, z_meters)
    ! call exit

    associate (A => z_interface)
      z_interface_val = bilinear_interpolation(x, x0, x1, y, y0, y1, &
          A(x0,y0), A(x0,y1), A(x1,y0), A(x1,y1))
    end associate

    ! sealevel_pressure => 100000.0
    ! pressure_val = pressure_at_elevation(100000.0, z_meters)
    ! exner_val =

    ! call random_number(rand)
    ! if (init_theta .eqv. .true.) then
    !   theta_val = theta_val * (1 + 1.0 / 100) ! random 0-1% change
    ! else
    !   theta_val = theta_val ! * (1 + 0.01) ! increase by 1%
    ! end if

    ! if (init_velocity .eqv. .true.) then
    velocity = 5
    ! else
    !   velocity = 0
    ! end if


    temp_val = exner_function(pressure_val) * theta_val

    ! if (dry_air_particles .eqv. .true.) then
    !    water_vapor_val = 0
    ! else
       water_vapor_val = sat_mr(temp_val, pressure_val) *  1.0 !0.99
    ! end if
    relative_humidity_in = water_vapor_val

    ! ARTLESS: Testing, set to 0
    ! relative_humidity_in = 0

    ! wind is constant in this system, ignoring w wind (aka z-direction)
    u_val = u_in%local(x0, z0, y0)
    v_val = v_in%local(x0, z0, y0)
    w_val = 0


    cloud_water = 0
    if (present(times_moved) .eqv. .true.) then
       times_moved_val = times_moved
    else
       times_moved_val = 0
    end if

    ! ARTLESS
    ! particle = convection_particle(particle_id, .true., times_moved_val, &
    !     x, y, z, u_val, v_val, w_val, z_meters, z_interface_val, &
    !     pressure_val, temp_val, theta_val, velocity, water_vapor_val, &
    !     cloud_water, relative_humidity_in)
    particle%particle_id = particle_id
    particle%exists = .true.
    particle%moved = times_moved_val
    particle%x = x
    particle%y = y
    particle%z = z
    particle%u = u_val
    particle%v = v_val
    particle%w = w_val
    particle%z_meters = z_meters
    particle%z_interface = z_interface_val
    particle%pressure = pressure_val
    particle%temperature = temp_val
    particle%potential_temp = theta_val
    particle%velocity = velocity
    particle%water_vapor = water_vapor_val
    particle%cloud_water = cloud_water
    particle%relative_humidity = relative_humidity_in
  end subroutine


  module subroutine process(this, nx_global, ny_global, &
      ims, ime, kms, kme, jms, jme, dt, dz, temperature, z_interface, &
      its, ite, kts, kte, jts, jte, z_m, potential_temp, pressure, u_in, v_in, &
      w_in, timestep)
    implicit none
    class(convection_exchangeable_t), intent(inout) :: this
    integer, intent(in) :: nx_global, ny_global
    real, intent(in)    :: dt, dz
    real, intent(in)    :: temperature(ims:ime,kms:kme,jms:jme)
    real, intent(in)    :: pressure(ims:ime,kms:kme,jms:jme)
    real, intent(in)    :: z_interface(ims:ime,jms:jme)
    integer, intent(in) :: ims, ime, kms, kme, jms, jme
    integer, intent(in) :: its, ite, kts, kte, jts, jte
    real, intent(in) :: z_m(ims:ime,kms:kme,jms:jme)
    class(exchangeable_t), intent(in) :: potential_temp
    class(exchangeable_t), intent(in), optional :: u_in, v_in, w_in
    integer, intent(in), optional :: timestep

#ifdef __NVCOMPILER
    real :: gravity = 9.80665
#else
    real, parameter :: gravity = 9.80665
#endif
    type(convection_particle) :: new_particle
    ! real :: Gamma
    real :: a_prime, z_displacement, T, t_prime, buoyancy
    real :: wind_correction, delta_z, z_interface_val, z_wind_change

    integer :: i,j,k, me
    real :: new_pressure, R_s, exponent, alt_pressure, alt_pressure2
    real :: alt_pressure3, alt_pressure4, mixing_ratio, sat_mr_val
    real :: vapor_p, sat_p, T_C, T_K, mr, T_squared, T_original, tmp
    real :: xx, yy, z_0, z_1
    real :: rate_of_temp_change, bv(local_buf_size), input_wind_val
    integer :: x0, x1, z0, z1, y0, y1, bv_i, image, particle_id(local_buf_size)
    integer :: u_bound(3)!, n_local_particles
    logical :: calc, calc_x, calc_y, exist
    ! logical, allocatable :: truth_a(:)
    character(len=32) :: filename
    logical :: dry_air_particles_val


    !-----------------------------------------
    !  block memory
    !-----------------------------------------
    real :: saturate, condensate, vapor, vapor_needed, RH
    ! specific latent heat values, calculating using formula
#ifdef __NVCOMPILER
    real :: condensation_lh = 2600!000 ! 2.5 x 10^6 J/kg
    real :: vaporization_lh = -2600
#else
    real, parameter :: condensation_lh = 2600!000 ! 2.5 x 10^6 J/kg
    real, parameter :: vaporization_lh = -condensation_lh
#endif
    ! specific heat of water vapor at constant volume
    ! real, parameter :: C_vv = 1.0 / 1390.0
    real :: C_vv, c_p
    real :: specific_latent_heat, Q_heat, delta_t, T0, T1, &
         p0, p1, q_dry, q_wet, potential_temp0, q_new_dry, q_dif
    real :: water_vapor0
    real :: cloud_water0,  cloud_water1
    real :: potential_T0, potential_T1
    integer :: iter
    logical :: bv_data_val, replacement_message_val, replacement_val
    logical :: debug_val, wrap_neighbors_val, convection_val

    !-----------------------------------------
    !  associate
    !-----------------------------------------
    real :: x, y, z
    type(convection_particle) :: particle

    ! extra variables needed because procedure calls aren't allowed
    type(convection_particle),allocatable :: local(:)

#ifdef EXPAND_FUNC
    real :: e_s,a,b
    real :: c000, c001, c010, c100, c011, c101, c110, c111
    real :: xd, zd, yd, c00, c01, c10, c11, c0, c1, c
    real :: gravity1, exner_function_val
    ! integer, allocatable :: test_all(:)
#endif




    input_wind_val = input_wind
    dry_air_particles_val = dry_air_particles

    replacement_message_val = replacement_message
    replacement_val = replacement
! #ifdef __NVCOMPILER
    debug_val = .false.
    bv_data_val = .false.
! #else
!     ! bv_data_val = brunt_vaisala_data
!     ! debug_val = debug
! #endif
    wrap_neighbors_val = wrap_neighbors
    convection_val = convection

    ! print *, "ARTLESS: PARTICLE PROCESSING"
    if (debug .eqv. .true.) print*, "start particle processing"
    if (bv_data_val .eqv. .true.) then
       bv_i = 1
       bv = 0
       particle_id = 0
    end if


    ! allocate(truth_a(ubound(this%local,1)))

    ! yy, xx, ny_global, nx_global, q_dif, q_new_dry, condensate,
    ! potential_temp0, delta_t, c_p, q_wet, q_heat, vapor, vapor_needed
    ! specific_latent_heat, cloud_water0, water_vapor0, potential_t1
    ! rh, saturate, iter, q_dry, p1, t1, only_dry, potential_t0, p0
    ! t0

#if NO_COARRAYS
    me = 1
#else
    me = this_image()
#endif

    allocate(local(size(this%local)))
    local = this%local

    ! allocate(test_all(size(this%local)))
    ! test_all = 0

#ifdef DC_LOOP
    ! ARTLESS
    do concurrent (i=1:current_max_local_particles) & ! this is 14% faster
         local(particle, bv_i,x0,z0,y0,x1,z1,y1,T,T_prime,x,y,z, &
#ifdef EXPAND_FUNC
         e_s,a,b,c000, c001, c010, c100, c011, c101, c110, c111, &
         xd, zd, yd, c00, c01, c10, c11, c0, c1, c, exner_function_val, &
#endif
         buoyancy,a_prime,z_displacement, delta_z, wind_correction, &
         z_interface_val,z_wind_change,z_0,z_1,T0,p0,potential_T0, &
         t1, p1, q_dry, iter,saturate, RH, potential_T1, water_vapor0, &
         cloud_water0, specific_latent_heat, vapor_needed, &
         vapor, Q_heat, q_wet, c_p, delta_t, condensate) &
         shared(local, temperature, dt, dz, z_interface, kte, kts) &
         default(none)
#endif
#ifdef OMP_LOOP
       ! !$omp parallel loop
       !!      !$omp parallel do simd schedule(auto)
!!!$omp parallel
!!$omp parallel do simd schedule(auto) &
!$omp  do simd schedule(auto) &
!$omp  private(T,T_prime,x0,z0,y0,x1,z1,y1,buoyancy,a_prime,z_displacement) &
!$omp  private(delta_z,wind_correction,z_interface_val,z_wind_change,z_0,z_1) &
!$omp  private(bv_i) &
!$omp  private(saturate, condensate, vapor, vapor_needed, RH, c_p) &
!$omp  private(specific_latent_heat, Q_heat, delta_t, T0, T1) &
!$omp  private(p0, p1, q_dry, q_wet) &
!$omp  private(water_vapor0, cloud_water0) &
!$omp  private(potential_T0, potential_T1, iter, x,y,z) &
#ifdef EXPAND_FUNC
!$omp  private(e_s,a,b,c000, c001, c010, c100, c011, c101, c110, c111) &
!$omp  private(xd, zd, yd, c00, c01, c10, c11, c0, c1, c, exner_function_val) &
#endif
#ifdef ONLY_IF_CHECKING_BV
!$omp  private(particle_id, bv) &
#endif
!$omp  private(particle)
!!$omp  shared(me,dt,dz,local,temperature,z_interface) &
!!$omp  shared(kte, kts) &
!!$omp  shared(z_m, pressure, potential_temp,input_wind_val,nx_global,ny_global) &
!!$omp  shared(w_in, v_in, u_in, dry_air_particles_val) &
!!$omp  shared(bv_data_val, replacement_message_val, current_max_local_particles) &
!!$omp  shared(replacement_val, debug_val, wrap_neighbors_val, convection_val) &
!!$omp  default(none)  ! remove gravity vaporization_lh condensation_lh
    do i=1,current_max_local_particles
       ! do i=1,ubound(this%local,1)
       ! print *, ""
#endif
#ifdef DO_LOOP
    do i=1,current_max_local_particles
    ! do i=1,ubound(this%local,1)
#endif

       particle = local(i) ! particle = this%local(i)

       ! print *, particle%particle_id, particle%exists
       ! if (i .eq. 1) print*, omp_in_parallel(), "omp_num_threads", omp_get_num_threads()

       if (particle%exists .neqv. .true.) cycle

! #ifndef __NVCOMPILER
!         if (particle%x .lt. its-1 .or. &
!              particle%z .lt. kts-1 .or. &
!              particle%y .lt. jts-1 .or. &
!              particle%x .gt. ite+1 .or. &
!              particle%z .gt. kte+1 .or. &
!              particle%y .gt. jte+1) then
!            print *, "x:", its, "<", particle%x, "<", ite, "with halo 2"
!            print *, "z:", kts, "<", particle%z, "<", kte, "with halo 2"
!            print *, "y:", jts, "<", particle%y, "<", jte, "with halo 2"
!            ! stop "x,y,z is out of bounds" ! can't be in DC
!            print *, "ERROR: x,y,z is out of bounds", " particle_id:", &
!                 particle%particle_id, "on image", me
!            ! print *, "particle_id:", &
!            !      particle%particle_id,"E=", particle%exists
!            particle%exists = .false.
! #ifdef DO_LOOP
!            stop "x,y,z is out of bounds" ! can't be in DC
! #endif
!            cycle
!         end if
! #endif


! #ifdef __NVCOMPILER
        if (.true. .eqv. .true.) bv_i = 1
! #else
!         if (bv_data_val .eqv. .true.) bv_i = 1
! #endif


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
        ! new: properites of environment are prime
        T = particle%temperature


#ifdef EXPAND_FUNC
        ! associate (A => temperature,
        !      y => particle%y)
        x = particle%x
        z = particle%z
        y = particle%y
        x0 = floor(x); z0 = floor(z); y0 = floor(y);
        x1 = ceiling(x); z1 = ceiling(z); y1 = ceiling(y);
    ! start trilinear function
          c000 = temperature(x0,z0,y0); c001 = temperature(x0,z0,y1); c010 = temperature(x0,z1,y0)
          c100 = temperature(x1,z0,y0)
          c011 = temperature(x0,z1,y1); c101 = temperature(x1,z0,y1); c110 = temperature(x1,z1,y0)
          c111 = temperature(x1,z1,y1)
    ! start bilinear check
    if ((x0 .eq. x1) .and. (z0 .eq. z1) .and. (y0 .eq. y1)) then
       c = c000
    else if ((x0 .eq. x1) .and. (z0 .eq. z1)) then
       c = c000 + (c011-c000) * (y-y0) / (y1-y0)
    else if ((x0 .eq. x1) .and. (y0 .eq. y1)) then
       c = c000 + (c011-c000) * (z-z0) / (z1-z0)
    else if ((y0 .eq. y1) .and. (z0 .eq. z1)) then
       c = c000 + (c110-c000) * (x-x0) / (x1-x0)

    else if ((x0 .eq. x1)) then
       c00 = c000; c01 = c001; c10 = c010; c11 = c011
       c0 = ((z1 - z)/(z1 - z0)) * c00 + ((z - z0)/(z1 - z0)) * c10
       c1 = ((z1 - z)/(z1 - z0)) * c01 + ((z - z0)/(z1 - z0)) * c11
       c  = ((y1 - y)/(y1 - y0)) * c0 + ((y - y0)/(y1 - y0)) * c1

    else if (y0 .eq. y1) then
       c00 = c000; c01 = c010; c10 = c100; c11 = c110
       c0 = ((x1 - x)/(x1 - x0)) * c00 + ((x - x0)/(x1 - x0)) * c10
       c1 = ((x1 - x)/(x1 - x0)) * c01 + ((x - x0)/(x1 - x0)) * c11
       c = ((z1 - z)/(z1 - z0)) * c0 + ((z - z0)/(z1 - z0)) * c1
    else if (z0 .eq. z1 ) then
       c00 = c000; c01 = c001; c10 = c100; c11 = c101
       c0 = ((x1 - x)/(x1 - x0)) * c00 + ((x - x0)/(x1 - x0)) * c10
       c1 = ((x1 - x)/(x1 - x0)) * c01 + ((x - x0)/(x1 - x0)) * c11
       c = ((y1 - y)/(y1 - y0)) * c0 + ((y - y0)/(y1 - y0)) * c1
    ! end bilinear
    else
       xd = (x - x0) / (x1 - x0)
       zd = (z - z0) / (z1 - z0)
       yd = (y - y0) / (y1 - y0)

       c00 = c000 * (1 - xd) + c100 * xd
       c01 = c001 * (1 - xd) + c101 * xd
       c10 = c010 * (1 - xd) + c110 * xd
       c11 = c011 * (1 - xd) + c111 * xd

       c0 = c00 * (1 - zd) + c10 * zd
       c1 = c01 * (1 - zd) + c11 * zd

       c = c0 * (1 - yd) + c1 * yd
       T_prime = c
    end if
    ! end trilinear
#else
    x0 = floor(particle%x); z0 = floor(particle%z)
    y0 = floor(particle%y)
    x1 = ceiling(particle%x); z1 = ceiling(particle%z)
    y1 = ceiling(particle%y)
    T_prime = trilinear_interpolation(particle%x, x0, x1, &
         particle%z, z0, z1, &
         particle%y, y0,y1,&
         temperature(x0,z0,y0), temperature(x0,z0,y1), &
         temperature(x0,z1,y0), temperature(x1,z0,y0), &
         temperature(x0,z1,y1), temperature(x1,z0,y1), &
         temperature(x1,z1,y0), temperature(x1,z1,y1))
#endif
     ! checked, the same


! #ifdef __NVCOMPILER
        if (.true. .eqv. .true.) then
! #else
!         if (convection_val .eqv. .true.) then
! #endif
           buoyancy = (T - T_prime) / T_prime
           a_prime = buoyancy * 9.80665  ! gravity
           ! d = v_0 * t + 1/2 * a * t^2
           z_displacement = particle%velocity * dt + 0.5 * a_prime * dt * dt
           ! z_displacement = particle%velocity + 0.5 * a_prime


           ! number from dz_interface, currently always 500
           delta_z = z_displacement / dz

           if (z_displacement /= z_displacement) then
#ifdef DO_LOOP
              print *, me, ":: ---------NAN ERROR---------", &
                   particle%particle_id, T, T_prime
              particle%exists = .false.
              cycle
              particle%z = -1
              print *, p0, z_displacement, gravity, particle%temperature
              ! call exit
#endif
           else
              particle%z = particle%z + delta_z
              particle%velocity = z_displacement
           end if
        else
           z_displacement = 0.0
        end if
     ! end do


        !-----------------------------------------------------------------
        ! Orographic lift
        ! Find dz change from change in environment
        !-----------------------------------------------------------------
        if (wind .eqv. .true.) then
           wind_correction = (1.0 / dz) ! real correction
           ! print *, "wind", wind_correction
           ! wind_correction = (1.0 / 10) ! real correction
           if (fake_wind_correction .eqv. .true.) then
              wind_correction = wind_correction * 0.0
           end if

           ! u: zonal velocity, wind towards the east
           ! v: meridional velocity, wind towards north
! #ifdef __NVCOMPILER
           if (.true. .eqv. .true.) then
              particle%x = particle%x + (8.0 * wind_correction)
              particle%y = particle%y + (8.0 * wind_correction)
! #else
!            if (use_input_wind .eqv. .true.) then
!               particle%x = particle%x + (input_wind_val * wind_correction)
!               particle%y = particle%y + (input_wind_val * wind_correction)
! #endif
           else
              particle%x = particle%x + (particle%u * wind_correction)
              particle%y = particle%y + (particle%v * wind_correction)
           end if

#ifdef EXPAND_FUNC
    ! start bilinear check
    x0 = floor(particle%x); x1 = ceiling(particle%x)
    y0 = floor(particle%y); y1 = ceiling(particle%y)

    c00 = z_interface(x0,y0); c01 = z_interface(x0,y1)
    c10 = z_interface(x1,y0); c11 = z_interface(x1,y1)
    if ((x0 .eq. x1) .and. (y0 .eq. y1)) then
       c = c00
    else if (x0 .eq. x1) then
       c = c00 + (c11-c00) * (particle%y-y0) / (y1-y0)
    else if (y0 .eq. y1) then
       c = c00 + (c11-c00) * (particle%x-x0) / (x1-x0)
    else
       c0 = ((x1 - particle%x)/(x1 - x0)) * c00 + ((particle%x - x0)/(x1 - x0)) * c10
       c1 = ((x1 - particle%x)/(x1 - x0)) * c01 + ((particle%x - x0)/(x1 - x0)) * c11
       c = ((y1 - particle%y)/(y1 - y0)) * c0 + ((particle%y - y0)/(y1 - y0)) * c1
    end if
    ! end bilinear
             z_interface_val = c

#else
             x0 = floor(particle%x); x1 = ceiling(particle%x)
             y0 = floor(particle%y); y1 = ceiling(particle%y)
             z_interface_val = bilinear_interpolation(particle%x, x0, x1, &
                  particle%y, y0, y1, &
                  z_interface(x0,y0), z_interface(x0,y1), &
                  z_interface(x1,y0), z_interface(x1,y1))
#endif
           ! checked

           z_wind_change = z_interface_val - particle%z_interface
           particle%z_interface = z_interface_val
           z_displacement = z_displacement + z_wind_change
        end if


        !-----------------------------------------------------------------
        ! Move particle, remove particle if beyond the z axis
        !-----------------------------------------------------------------
        z_0 = particle%z_meters
        z_1 = particle%z_meters + z_displacement
        particle%z_meters = z_1

        ! out of bounds is handled in data movement do loop
        if (particle%z .lt. kts) then
           cycle
        else if (particle%z .gt. kte) then
           cycle
        end if


#ifdef ONLY_IF_CHECKING_BV
        if (bv_data_val .eqv. .true.) then
           x0 = floor(particle%x)
           z0 = floor(particle%z)
           y0 = floor(particle%y)
           x1 = ceiling(particle%x)
           z1 = ceiling(particle%z)
           y1 = ceiling(particle%y);

           bv(bv_i) = trilinear_interpolation( &
                particle%x, x0, x1, &
                particle%z, z0, z1, &
                particle%y, y0, y1, &
                potential_temp%local(x0,z0,y0), &
                potential_temp%local(x0,z0,y1), &
                potential_temp%local(x0,z1,y0), &
                potential_temp%local(x1,z0,y0), &
                potential_temp%local(x0,z1,y1), &
                potential_temp%local(x1,z0,y1), &
                potential_temp%local(x1,z1,y0), &
                potential_temp%local(x1,z1,y1))
           particle_id(bv_i) = particle%particle_id
           bv_i = bv_i + 1
        end if
#endif
    ! end do
    ! print *, "DONE"

        !-----------------------------------------------------------------
        ! Dry Lapse Rate
        !-----------------------------------------------------------------
        ! p = p_0 * e^( ((9.81/287.058)*dz) / t_mean )
        !-----------------------------------------------------------------
        ! Method one: physics
        !        a) change pressure         b) update temperature
        !-----------------------------------------------------------------
! #ifdef __NVCOMPILER
        if (.false. .eqv. .true.) then !artless, parcels wet
           if (DEBUG_VAL .eqv. .true.) print *, "-- only dry air parcels --"
! #else
!         if (dry_air_particles_val .eqv. .true.) then
!            if (debug_val .eqv. .true.) print *, "-- only dry air parcels --"
! #endif
#ifdef EXPAND_FUNC
           ! --- dry_lapse_rate function ---
           ! particle%pressure = particle%pressure - z_displacement * &
           !      gravity / (287.05 * particle%temperature) * particle%pressure
           ! particle%temperature = exner_function(particle%pressure) * &
           !      particle%potential_temp
           ! -------------------------------
             particle%pressure = particle%pressure - z_displacement * &
                  9.80665 / (287.05 * particle%temperature) * particle&
                  &%pressure
             exner_function_val = (particle%pressure/100000)**(287.058/1003.5)
             particle%temperature = exner_function_val * &
                  particle%potential_temp
           ! -------------

#else
             call dry_lapse_rate(particle%pressure, particle%temperature, &
                  particle%potential_temp, z_displacement)
#endif
        else
           !-----------------------------------------------------------------
           ! Relative Humidity and physics of Moist Adiabatic Lapse Rate
           !-----------------------------------------------------------------
           ! | Mixing Ratio |
           ! saturate = sat_mr(t,p)
           ! if (water_vapor > saturated)
           !   condensate = water_vapor - satured
           !   water_vapor -= condensate
           !   clouds += condensate
           !   Q_heat  = specific_latent_heat * condensate
           !   delta_T = Q_heat / c_p   ! c_p is specific heat capacity
           !   temperature += delta_T
           !-----------------------------------------------------------------

           ! ---- this was start of block ----
             T0 = particle%temperature
             p0 = particle%pressure
             potential_T0 = particle%potential_temp

#ifdef EXPAND_FUNC
             particle%pressure = particle%pressure - z_displacement * &
                  9.80665 / (287.05 * particle%temperature) * particle&
                  &%pressure
             exner_function_val = (particle%pressure/100000)**(287.058/1003.5)
             particle%temperature = exner_function_val * &
                  particle%potential_temp
#else
             call dry_lapse_rate(particle%pressure, particle%temperature, &
                  particle%potential_temp, z_displacement)
#endif

             T1 = particle%temperature
             p1 = particle%pressure
             q_dry = 1004 * 1 * (abs(t1-t0))  ! q = c_p x m x delta_T


             do iter = 1,4 ! ARTLESS
#ifdef EXPAND_FUNC
                if (particle%temperature < 273.15) then
                   a = 21.8745584
                   b = 7.66
                else
                   a = 17.2693882
                   b = 35.86
                endif
                e_s = 610.78 * exp(a * (particle%temperature - 273.16) / &
                     (particle%temperature - b)) !(Pa)
                if ((particle%pressure - e_s) <= 0) then
                   e_s = particle%pressure * 0.99999
                endif
                saturate = 0.6219907 * e_s / (particle%pressure - e_s) !(kg/kg)
#else
                saturate = sat_mr(particle%temperature, particle%pressure)
#endif

                RH = particle%water_vapor / saturate
                potential_T1 = particle%potential_temp
                water_vapor0 = particle%water_vapor
                cloud_water0 = particle%cloud_water
                particle%relative_humidity = RH
                ! https://en.wikipedia.org/wiki/Latent_heat#Specific_latent_heat
                ! specific latent heat for condensation and evaporation
                specific_latent_heat = 2.5 * 10**6 ! J kg^-1
                ! Particle is falling, using evaporation of cloud water to keep
                ! the relative humidity at 1 if possible
                vapor_needed = saturate - particle%water_vapor
                ! print *, ":",saturate, "-", particle%water_vapor
                ! print *, "----",particle%cloud_water, RH, vapor_needed

                if (particle%cloud_water .gt. 0.0 .and. RH .lt. 1.0 &
                     .and. vapor_needed .gt. 0.0000001) then
#ifdef DO_LOOP
                   if (debug_val .eqv. .true.) &
                        print*, "==== cloud_water .gt. 0, rh .lt. 1, wet ===="
#endif
                   ! test_all(i) = 1

                   if (vapor_needed > particle%cloud_water) then
                      vapor = particle%cloud_water
                      vapor = vapor / (6-iter)                ! Saturated
                      particle%cloud_water = 0
                   else
                      vapor = vapor_needed
                      vapor = vapor / (6-iter)                ! Saturated
                      particle%cloud_water = particle%cloud_water - vapor
                   end if


                   particle%water_vapor = particle%water_vapor + vapor

                   ! heat required by phase change
                   Q_heat = specific_latent_heat * vapor ! kJ
                   q_wet = Q_heat
                   c_p = (1004 * (1 + 1.84 * vapor)) ! 3.3
                   ! c_p = 1004 ! specific heat of dry air at 0C
                   delta_t = Q_heat / c_p   ! 3.2c
                   particle%temperature = T1 - delta_t

#ifdef DO_LOOP
                   if (debug_val .eqv. .true.) &
                        potential_temp0 = particle%potential_temp
#endif

                   ! update potential temperature, assumming pressure is
                   ! constant
#ifdef EXPAND_FUNC
             exner_function_val = (particle%pressure/100000)**(287.058/1003.5)
             particle%potential_temp = particle%temperature /&
                     exner_function_val
#else
                   particle%potential_temp = particle%temperature / exner_function(particle%pressure)
#endif
                   ! Particle is raising and condensation is occurring


                else if (RH .gt. 1.0) then
                   ! test_all(i) = 2
                   ! ==== rh .gt. 1 ====
                   ! water_vapor 5.271018017E-4 saturate 5.271018017E-4 cloud water 0.
                   ! condensate 0. new water_vapor 5.271018017E-4 new cloud water 0.
                   ! Q_heat 0.
                   ! ============= process done ===============

#ifdef DO_LOOP
                   if (debug_val .eqv. .true.) print*, "==== rh .gt. 1, wet ===="
#endif

                   condensate = particle%water_vapor - saturate
                   condensate = condensate / (6-iter) ! Saturated
                   particle%water_vapor = particle%water_vapor - condensate
                   particle%cloud_water = particle%cloud_water + condensate

                   !--------------------------------------------------------------
                   ! a different way to calculate specific latent heat
                   !--------------------------------------------------------------
                   ! if (T1 .lt. 248.15) then
                   !   specific_latent_heat = 2600
                   ! else if (T1 .gt. 314.15) then
                   !   specific_latent_heat = 2400
                   ! else
                   !   T_C = T1 - 273.15
                   !   specific_latent_heat = 2500.8 - 2.36*T_C + 0.0016 * T_C**2 - &
                   !       0.00006 * T_C**3
                   ! end if

                   ! heat from phase change
                   Q_heat = specific_latent_heat * condensate ! kJ
                   q_wet = Q_heat
                   ! print *, "Q_heat", Q_heat
                   ! exit
                   !--------------------------------------------------------------
                   ! calculate change in heat using specific heat (c_p)
                   !--------------------------------------------------------------
                   ! Stull: Practical Meteorology
                   c_p = (1004 * (1 + 1.84 * condensate)) ! 3.3
                   ! c_p = 1004 ! specific heat of dry air at 0C
                   ! print *, "c_p", c_p
                   delta_T = Q_heat / c_p   ! 3.2c
                   ! print *, "delta_T", delta_t
                   particle%temperature = T1 + delta_T


                   ! update potential temperature, assumming pressure is constant
                   ! particle%potential_temp = exner_function(particle%pressure) / &
                   !      particle%temperature
#ifdef EXPAND_FUNC
             exner_function_val = (particle%pressure/100000)**(287.058/1003.5)
             particle%potential_temp = particle%temperature /&
                  exner_function_val

#else
                   particle%potential_temp = particle%temperature /&
                        & exner_function(particle%pressure)
#endif
                   ! call dry_lapse_rate(particle%pressure, particle%temperature, &
                   !      particle%potential_temp, z_displacement)

                   ! stop
                   ! -- is pressure constant during this process?
                   ! using Poisson's equation
                   ! particle%pressure = p0/ ((T0/particle%temperature)**(1/0.286))

                else if (iter .eq. 1) then
                   ! test_all(i) = 3
#ifdef DO_LOOP
                   if (debug_val .eqv. .true.) &
                      print*, "==== was dry process ===="
#endif
                   ! cycle
                   exit
                end if
                ! test_all(i) = 4


! #ifdef __NVCOMPILER
!                 if (DEBUG_VAL .eqv. .true.) then
! #else
!                 if (debug_val .eqv. .true.) then
! #endif
!                    print *, "     pressure  |  temp      |   ~heat    | potential"
!                    print *, "pre ", p0, t0, ", -none-       ,", potential_t0

!                    q_new_dry = &
!                         1004 * (1) * (abs(particle%temperature-t0))
!                    q_dif = abs(q_dry-(q_new_dry+q_wet))
!                    print *, q_dry, q_new_dry, q_wet, q_dif, &
!                         1004 * (1) * q_dif
!                 end if
!                 ! if ((debug_val .eqv. .true.) .and. (only_dry .eqv. .false.)) then
!                 !    print *, "post", p1, t1, q_dry, potential_t1
!                 !    print *, "new ", &
!                 !         particle%pressure, particle%temperature, &
!                 !         q_wet,  particle%potential_temp


!                 ! print *, "heat: q_dry = q_new_dry + q_wet"
!                 ! print *, " ", q_dry, "=", q_new_dry, "+", q_wet
!                 ! print *, "dif         =", q_dif, "which is ~ temp diff", &
!                 !      1004 * (1) * q_dif

!                 ! print *, "--test--"
!                 ! print *, " t? ", t0 - (q_new_dry / (1004 * (1)) )

!                 ! print *, "vapor_needed", vapor_needed, &
!                 !      "new water_vapor", particle%water_vapor, &
!                 !      "new cloud water", particle%cloud_water, &
!                 !      "water_vapor0", water_vapor0, &
!                 !      "cloud water0", cloud_water0
!                 ! end if

             end do
!              ! saturate = sat_mr(particle%temperature, particle%pressure)
!              ! RH = particle%water_vapor / saturate
!              ! particle%relative_humidity = RH


        end if ! ---- end of relative humidity seciton ----
#ifdef OMP_LOOP
     end do
#ifdef __NVCOMPILER
     !$omp end do
#else
     !$omp end do simd
     !!$omp end parallel do simd
#endif

   !!
   !!!$omp end parallel loop
#else
     end do
#endif


    ! print *, test_all(1)
    ! print *, "DONE"

!     !     !-----------------------------------------------------------------
!     !     ! Move particle if needed
!         !-----------------------------------------------------------------
    do i=1,current_max_local_particles
    ! ! do i=1,ubound(this%local,1)     !
      associate (particle=>this%local(i))        !
        if (particle%exists .neqv. .true.) cycle !
        x = particle%x
        y = particle%y
        z = particle%z


        if (particle%z .lt. kts) then
           particle%exists = .false.
           if (replacement_val .eqv. .true.) then
              call replace_particle(particle%particle_id, &
                   its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, &
                   z_m, potential_temp, z_interface, pressure, u_in, v_in, w_in,&
                   particle%moved, particle)
              if (replacement_message_val .eqv. .true.) then
                 print *, me,":",particle%particle_id, "hit the ground"
              end if
           end if
           cycle
        else if (particle%z .gt. kte) then
           particle%exists = .false.
           if (replacement_val .eqv. .true.) then
              call replace_particle(particle%particle_id, &
                   its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, &
                   z_m, potential_temp, z_interface, pressure, u_in, v_in, w_in,&
                   particle%moved, particle)
              if (replacement_message_val .eqv. .true.) then
                 print *, me,":",particle%particle_id, "went off the top"
              end if
           end if
           cycle
        end if


#ifdef DO_LOOP
        if (caf_comm_message .eqv. .true.) then
           if (  x .lt. its-1 .or. x .gt. ite+1 .or. &
                z .lt. kms-1 .or. z .gt. kme   .or. &
                y .lt. jts-1 .or. y .gt. jte+1 .or. &
                x .lt. 1 .or. x .gt. nx_global .or. &
                y .lt. 1 .or. y .gt. ny_global &
                ) then
              print *, "PUTTING", particle%x, particle%y, particle%z_meters, &
                   "FROM", me, "id:", particle%particle_id, &
                   "M", ims,ime, kms, kme, jms, jme, "T", its,ite,jts,jte
           end if
        end if
#endif
        xx = x
        yy = y
        ! If particle is getting wrapped the x and y values need to be
        ! properly updated
        if (wrap_neighbors_val .eqv. .true.) then
           if (x > nx_global) then
              ! if (caf_comm_message .eqv. .true.) print *, "WRAPPED" &
              !     , particle%particle_id
              x = x - nx_global + 1
              xx = xx + 2
           else if (x < 1) then
              ! if (caf_comm_message .eqv. .true.) print *, "WRAPPED" &
              !     , particle%particle_id
              x = x + nx_global - 1
              xx = xx - 2
           end if

           if (y > ny_global) then
              ! if (caf_comm_message .eqv. .true.) print *, "WRAPPED" &
              !     , particle%particle_id
              y = y - ny_global + 1
              yy = yy + 2
           else if (y < 1) then
              ! if (caf_comm_message .eqv. .true.) print *, "WRAPPED" &
              !     , particle%particle_id
              y = y + ny_global - 1
              yy = yy - 2
           end if
        end if

        ! Check values to know where to send particle
#if NO_COARRAYS
        if (1 .gt. 1) then
#else
        if (num_images() .gt. 1) then
#endif
        if (yy .lt. jts-1) then      ! "jts  <   y    <  jte"
           if (xx .lt. its-1) then
              ! particle%exists = .false.
              call this%put_southwest(particle)
              ! if (count_p_comm .eqv. .true.) &
              ! particles_communicated = particles_communicated + 1
           else if (xx .gt. ite+1) then
              ! particle%exists = .false.
              call this%put_southeast(particle)
              ! if (count_p_comm .eqv. .true.) &
              ! particles_communicated = particles_communicated + 1
           else
              ! particle%exists = .false.
              call this%put_south(particle)
              ! if (count_p_comm .eqv. .true.) &
              ! particles_communicated = particles_communicated + 1
           end if
        else if (yy .gt. jte+1) then ! jts will be 1
           if (xx .lt. its-1) then
              ! particle%exists = .false.
              call this%put_northwest(particle)
              ! if (count_p_comm .eqv. .true.) &
              ! particles_communicated = particles_communicated + 1
           else if (xx .gt. ite+1) then
              ! particle%exists = .false.
              call this%put_northeast(particle)
              ! if (count_p_comm .eqv. .true.) &
              ! particles_communicated = particles_communicated + 1
           else
              ! particle%exists = .false.
              call this%put_north(particle)
              ! if (count_p_comm .eqv. .true.) &
              ! particles_communicated = particles_communicated + 1
           endif
        else if (xx .lt. its-1) then ! "its  <   x    <  ite"
           ! particle%exists = .false.
           call this%put_west(particle) ! need to double check this!
           ! if (count_p_comm .eqv. .true.) &
           ! particles_communicated = particles_communicated + 1
        else if (xx .gt. ite+1) then
           ! particle%exists = .false.
           call this%put_east(particle)
           ! if (count_p_comm .eqv. .true.) &
           ! particles_communicated = particles_communicated + 1
        end if
        end if

        ! end if
      end associate
        ! print *, "end loop"
   end do
   return

!     if (brunt_vaisala_data .eqv. .true.) then
! #if NO_COARRAYS
!       do image=1,1
! #else
!       do image=1,num_images()
! #endif
!       if (me .eq. image) then
!          write (filename,"(A17)") "brunt_vaisala.txt"
!          inquire(file=filename, exist=exist)
!          if (exist) then
!             open(unit=me, file=filename, status='old', position='append')
!          else
!             open(unit=me, file=filename, status='new')
!          end if
!          do i=1,bv_i-1
!             write(me,*) me, timestep, particle_id(i), bv(i)
!          end do
!          close(me)
!       end if
! #ifndef NO_COARRAYS
!       sync all
! #endif
!       end do
 ! ARTLESS
!    end if
!    ! if (debug_val .eqv. .true.) print *, "============= process done ==============="
  end subroutine process



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

#ifndef NO_COARRAYS
    if (.not. present(no_sync)) then
      ! sync images (neighbors) ! sync neighbors currently broken
      sync all
    else
      if (.not. no_sync) then
        ! sync images (neighbors) ! sync neighbors currently broken
        sync all
      endif
    endif
#endif
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
#ifndef NO_COARRAYS
    sync all
#endif
  end subroutine

  module subroutine exchange(this)
    class(convection_exchangeable_t), intent(inout) :: this
    ! if (.not. this%north_boundary) call this%put_north
    ! if (.not. this%south_boundary) call this%put_south
    ! if (.not. this%east_boundary)  call this%put_east
    ! if (.not. this%west_boundary)  call this%put_west

#ifndef NO_COARRAYS
    sync images( neighbors )
#endif

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

  module subroutine retrieve_buf(this, buf)
    implicit none
    class(convection_exchangeable_t), intent(inout) :: this
#if NO_COARRAYS
    type(convection_particle), intent(inout) :: buf(:)
#else
    type(convection_particle), intent(inout) :: buf(:)[*]
#endif
    integer :: i, buf_n, local_i, local_n
    buf_n = ubound(buf, dim=1)
    local_n = ubound(this%local, dim=1)
    local_i = 1
    do i=1,buf_n
       associate (particle=>buf(i))
         do while (particle%exists .eqv. .true.)
            if (local_i .gt. local_n) then
               ! print *, this_image(), ": id",particle%particle_id," i" , &
               !      local_i, "but buf size", local_n
               ! print *, "ERROR: retrieve_buf is out of bounds"
               particle%exists = .false.
            end if

            if (this%local(local_i)%exists .eqv. .false.) then
               call particle%move_to(this%local(local_i))
               ! print *, "------", particle%exists
               local_i = local_i + 1
               exit
            end if
            local_i = local_i + 1
         end do
       end associate
    end do
    if (local_i > current_max_local_particles) then
       current_max_local_particles = local_i
    end if
  end subroutine

  module subroutine put_north(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", north_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", north_con_neighbor, particle%particle_id
#endif
    end if

    if (this%north_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%south_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_south_in(this%south_i)[north_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%south_i = this%south_i + 1
  end subroutine

  module subroutine put_south(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", south_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", south_con_neighbor, particle%particle_id
#endif
    end if

    if (this%south_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%north_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_north_in(this%north_i)[south_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%north_i = this%north_i + 1
  end subroutine

  module subroutine put_east(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", east_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", east_con_neighbor, particle%particle_id
#endif
    end if

    if (this%east_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%west_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_west_in(this%west_i)[east_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%west_i = this%west_i + 1
  end subroutine

  module subroutine put_west(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", west_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", west_con_neighbor, particle%particle_id
#endif
    end if

    if (this%west_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%east_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_east_in(this%east_i)[west_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%east_i = this%east_i + 1
  end subroutine

  module subroutine put_northeast(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", northeast_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", northeast_con_neighbor, particle%particle_id
#endif
    end if

    if (this%northeast_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%southwest_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_southwest_in(this%southwest_i)[northeast_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%southwest_i = this%southwest_i + 1
  end subroutine

  module subroutine put_northwest(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", northwest_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", northwest_con_neighbor, particle%particle_id
#endif
    end if

    if (this%northwest_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%southeast_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_southeast_in(this%southeast_i)[northwest_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%southeast_i = this%southeast_i + 1
  end subroutine

  module subroutine put_southeast(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", southeast_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", southeast_con_neighbor, particle%particle_id
#endif
    end if

    if (this%southeast_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%northwest_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_northwest_in(this%northwest_i)[southeast_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%northwest_i = this%northwest_i + 1
  end subroutine

  module subroutine put_southwest(this, particle)
    class(convection_exchangeable_t), intent(inout) :: this
    type(convection_particle), intent(inout) :: particle
    if (caf_comm_message .eqv. .true.) then
#if NO_COARRAYS
       print*, "from", 1, "to", southwest_con_neighbor, particle%particle_id
#else
       print*, "from", this_image(), "to", southwest_con_neighbor, particle%particle_id
#endif
    end if

    if (this%southwest_boundary) then
      particle%exists = .false.
      return
    end if

    call check_buf_size(this%northeast_i)
#ifndef NO_COARRAYS
    !dir$ pgas defer_sync
    this%buf_northeast_in(this%northeast_i)[southwest_con_neighbor] = particle
#endif
    particle%exists = .false.
    this%northeast_i = this%northeast_i + 1
  end subroutine

  ! pure module subroutine put_north_p(this, particle,i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%north_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_south_in(i)[north_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine


  ! pure module subroutine put_south_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%south_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_north_in(i)[south_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine

  ! pure module subroutine put_east_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%east_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_west_in(i)[east_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine

  ! pure module subroutine put_west_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%west_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_east_in(i)[west_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine

  ! pure module subroutine put_northeast_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%northeast_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_southwest_in(i)[northeast_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine

  ! pure module subroutine put_northwest_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%northwest_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_southeast_in(i)[northwest_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine

  ! pure module subroutine put_southeast_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%southeast_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_northwest_in(i)[southeast_con_neighbor] = particle
  !   particle%exists = .false.
  ! end subroutine

  ! pure module subroutine put_southwest_p(this, particle, i)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   type(convection_particle), intent(inout) :: particle
  !   integer, intent(in) :: i
  !   if (this%southwest_boundary) then
  !     particle%exists = .false.
  !     return
  !   end if

  !   !dir$ pgas defer_sync
  !   this%buf_northeast_in(this%northeast_i)[southwest_con_neighbor] = particle
  !   particle%exists = .false.
  !   this%northeast_i = this%northeast_i + 1
  ! end subroutine


  module subroutine create_particle_id(this)
    use iso_fortran_env, only : int32
    implicit none
    class(convection_exchangeable_t), intent(inout) :: this
    integer :: id_range, h, me, n_images
#if NO_COARRAYS
    n_images = 1
    me = 1
#else
    n_images = num_images()
    me = this_image()
#endif
    if (this%particle_id_count .eq. -1) then
      id_range = huge(int32) / n_images
      this%particle_id_count = (me-1) * id_range
    else
      this%particle_id_count = this%particle_id_count + 1
    end if
  end subroutine

  module subroutine setup_neighbors(this, grid)
    implicit none
    class(convection_exchangeable_t), intent(inout) :: this
    type(grid_t), intent(in)      :: grid
    integer :: current, n_neighbors, n_images

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


    if (.not.allocated(neighbors)) then
#if NO_COARRAYS
      associate(me=>1)
#else
      associate(me=>this_image())
#endif
        north_con_neighbor = me + grid%ximages
        south_con_neighbor = me - grid%ximages
        east_con_neighbor  = me + 1
        west_con_neighbor  = me - 1
        northeast_con_neighbor = me + grid%ximages + 1
        northwest_con_neighbor = me + grid%ximages - 1
        southeast_con_neighbor = me - grid%ximages + 1
        southwest_con_neighbor = me - grid%ximages - 1

        n_neighbors = &
             merge(0,1,this%south_boundary) &
            +merge(0,1,this%north_boundary) &
            +merge(0,1,this%east_boundary)  &
            +merge(0,1,this%west_boundary)  &
            +merge(0,1,this%northeast_boundary) &
            +merge(0,1,this%northwest_boundary) &
            +merge(0,1,this%southeast_boundary) &
            +merge(0,1,this%southwest_boundary)
        n_neighbors = max(1, n_neighbors)

        allocate(neighbors(n_neighbors))

        current = 1
        if (.not. this%south_boundary) then
          neighbors(current) = south_con_neighbor
          current = current+1
        endif
        if (.not. this%north_boundary) then
          neighbors(current) = north_con_neighbor
          current = current+1
        endif
        if (.not. this%east_boundary) then
          neighbors(current) = east_con_neighbor
          current = current+1
        endif
        if (.not. this%west_boundary) then
          neighbors(current) = west_con_neighbor
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
#if NO_COARRAYS
            n_images = 1
#else
            n_images = num_images()
#endif
            ! --- handle diagonals
            if (this%north_boundary .eqv. .true.) then
              northeast_con_neighbor = modulo( (me+nx+1)-nimages, (nx+1))
              if (northeast_con_neighbor .eq. 0) northeast_con_neighbor = 1
              northwest_con_neighbor = (me+nx-1)-nimages
              if (northwest_con_neighbor .eq. 0) northwest_con_neighbor = nx
            else  if (this%east_boundary .eqv. .true.) then
              northeast_con_neighbor = me + 1
            else  if (this%west_boundary .eqv. .true.) then
              northwest_con_neighbor = me + nx * 2 - 1
            end if

            if (this%south_boundary .eqv. .true.) then
              southeast_con_neighbor = me + nimages - nx + 1
              if (southeast_con_neighbor > nimages) then
                southeast_con_neighbor = southeast_con_neighbor - nx
              end if
              southwest_con_neighbor = me + nimages - nx - 1
              if (modulo(southwest_con_neighbor,nx) == 0) then
                southwest_con_neighbor = nimages
              end if
              ! southwest_con_neighbor = - 2
            else  if (this%east_boundary .eqv. .true.) then
              southeast_con_neighbor = me - nx * 2 + 1
            else  if (this%west_boundary .eqv. .true.) then
              southwest_con_neighbor = me - 1
              ! southwest_con_neighbor = - 1
            end if
            this%northeast_boundary = .false.
            this%northwest_boundary = .false.
            this%southeast_boundary = .false.
            this%southwest_boundary = .false.

            ! --- handle up/down/left/right
            if (this%north_boundary .eqv. .true.) then
              north_con_neighbor = north_con_neighbor - n_images
              this%north_boundary = .false.
              this%wrapped_north = .true.
            end if
            if (this%south_boundary .eqv. .true.) then
              south_con_neighbor = south_con_neighbor + n_images
              this%south_boundary = .false.
              this%wrapped_south = .true.
            end if
            if (this%east_boundary .eqv. .true.) then
              east_con_neighbor = east_con_neighbor - grid%ximages
              this%east_boundary = .false.
              this%wrapped_east = .true.
            end if
            if (this%west_boundary .eqv. .true.) then
              west_con_neighbor = west_con_neighbor + grid%ximages
              this%west_boundary = .false.
              this%wrapped_west = .true.
            end if

            ! fix neighbors for when parcels are wrapped
            deallocate(neighbors)
            allocate(neighbors(8))
            neighbors = (/ &
                 north_con_neighbor, &
                 south_con_neighbor, &
                 east_con_neighbor, &
                 west_con_neighbor, &
                 northeast_con_neighbor, &
                 southeast_con_neighbor, &
                 northwest_con_neighbor, &
                 southwest_con_neighbor /)
          end if
        end associate
      end associate
    end if

  end subroutine


  pure function bilinear_interpolation(x, x0, x1, y, y0, y1, &
      c00, c01, c10, c11) result(c)
    real, intent(in) :: x, y
    real, intent(in) :: c00, c01, c10, c11
    integer, intent(in) :: x0, x1, y0, y1
    real :: xd, yd, c0, c1, c
    if (x0 .eq. x1) then
       if (y0 .eq. y1) then
          c = c00
          return
       end if
      c = c00 + (c11-c00) * (y-y0) / (y1-y0)
    else if (y0 .eq. y1) then
      c = c00 + (c11-c00) * (x-x0) / (x1-x0)
    else
      c0 = ((x1 - x)/(x1 - x0)) * c00 + ((x - x0)/(x1 - x0)) * c10
      c1 = ((x1 - x)/(x1 - x0)) * c01 + ((x - x0)/(x1 - x0)) * c11
      c = ((y1 - y)/(y1 - y0)) * c0 + ((y - y0)/(y1 - y0)) * c1
    end if
  end function bilinear_interpolation

  pure function trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
      c000, c001, c010, c100, c011, c101, c110, c111) result(c)
    real, intent(in) :: x, z, y
    real, intent(in) :: c000, c001, c010, c100, c011, c101, c110, c111
    integer, intent(in) :: x0, x1, z0, z1, y0, y1
    real :: xd, zd, yd, c00, c01, c10, c11, c0, c1, c
    if (x0 .eq. x1) then
      c = bilinear_interpolation(z, z0, z1, y, y0, y1, &
          c000, c001, c010, c011)
    else if (y0 .eq. y1) then
      c = bilinear_interpolation(x, x0, x1, z, z0, z1, &
          c000, c010, c100, c110)
    else if (z0 .eq. z1 ) then
      c = bilinear_interpolation(x, x0, x1, y, y0, y1, &
          c000, c001, c100, c101)
    else
      xd = (x - x0) / (x1 - x0)
      zd = (z - z0) / (z1 - z0)
      yd = (y - y0) / (y1 - y0)

      c00 = c000 * (1 - xd) + c100 * xd
      c01 = c001 * (1 - xd) + c101 * xd
      c10 = c010 * (1 - xd) + c110 * xd
      c11 = c011 * (1 - xd) + c111 * xd

      c0 = c00 * (1 - zd) + c10 * zd
      c1 = c01 * (1 - zd) + c11 * zd

      c = c0 * (1 - yd) + c1 * yd
   end if

    ! ! expanded version of trilinear interpolation
    ! if ((x0 .eq. x1) .and. (z0 .eq. z1) .and. (y0 .eq. y1)) then
    !    c = c000
    ! else if ((x0 .eq. x1) .and. (z0 .eq. z1)) then
    !    c = c000 + (c011-c000) * (y-y0) / (y1-y0)
    ! else if ((x0 .eq. x1) .and. (y0 .eq. y1)) then
    !    c = c000 + (c011-c000) * (z-z0) / (z1-z0)
    ! else if ((y0 .eq. y1) .and. (z0 .eq. z1)) then
    !    c = c000 + (c110-c000) * (x-x0) / (x1-x0)
    ! ! end linear

    ! else if ((x0 .eq. x1)) then
    !    c00 = c000; c01 = c001; c10 = c010; c11 = c011
    !    c0 = ((z1 - z)/(z1 - z0)) * c00 + ((z - z0)/(z1 - z0)) * c10
    !    c1 = ((z1 - z)/(z1 - z0)) * c01 + ((z - z0)/(z1 - z0)) * c11
    !    c  = ((y1 - y)/(y1 - y0)) * c0 + ((y - y0)/(y1 - y0)) * c1

    ! else if (y0 .eq. y1) then
    !    c00 = c000; c01 = c010; c10 = c100; c11 = c110
    !    c0 = ((x1 - x)/(x1 - x0)) * c00 + ((x - x0)/(x1 - x0)) * c10
    !    c1 = ((x1 - x)/(x1 - x0)) * c01 + ((x - x0)/(x1 - x0)) * c11
    !    c = ((z1 - z)/(z1 - z0)) * c0 + ((z - z0)/(z1 - z0)) * c1
    ! else if (z0 .eq. z1 ) then
    !    c00 = c000; c01 = c001; c10 = c100; c11 = c101
    !    c0 = ((x1 - x)/(x1 - x0)) * c00 + ((x - x0)/(x1 - x0)) * c10
    !    c1 = ((x1 - x)/(x1 - x0)) * c01 + ((x - x0)/(x1 - x0)) * c11
    !    c = ((y1 - y)/(y1 - y0)) * c0 + ((y - y0)/(y1 - y0)) * c1
    ! ! end bilinear
    ! else
    !    xd = (x - x0) / (x1 - x0)
    !    zd = (z - z0) / (z1 - z0)
    !    yd = (y - y0) / (y1 - y0)

    !    c00 = c000 * (1 - xd) + c100 * xd
    !    c01 = c001 * (1 - xd) + c101 * xd
    !    c10 = c010 * (1 - xd) + c110 * xd
    !    c11 = c011 * (1 - xd) + c111 * xd

    !    c0 = c00 * (1 - zd) + c10 * zd
    !    c1 = c01 * (1 - zd) + c11 * zd

    !    c = c0 * (1 - yd) + c1 * yd
    ! ! end trilinear
    ! end if
  end function trilinear_interpolation

  module function total_num_particles()
    integer :: total_num_particles
#if NO_COARRAYS
    total_num_particles = particles_per_image * 1
#else
    total_num_particles = particles_per_image * num_images()
#endif
  end function total_num_particles

  module function num_particles_per_image()
    integer :: num_particles_per_image
    num_particles_per_image = particles_per_image
  end function num_particles_per_image

  module function num_particles_communicated()
    integer :: num_particles_communicated
    if (count_p_comm .eqv. .true.) then
       call co_sum(particles_communicated)
       num_particles_communicated = particles_communicated
    else
       num_particles_communicated = -1
       print *, "WARNING: number of particles communicated wasn't calculated"
    end if
  end function num_particles_communicated

  module function are_particles_dry()
    logical :: are_particles_dry
    are_particles_dry = dry_air_particles
  end function are_particles_dry

  module function get_wind_speed()
    real :: get_wind_speed
    if (use_input_wind .eqv. .true.) then
      get_wind_speed = input_wind
    end if
  end function get_wind_speed

  module function current_num_particles(convection_obj)
    integer :: current_num_particles
    type(convection_exchangeable_t),intent(in) :: convection_obj
    integer :: buf_size, i
    buf_size = size(convection_obj%local)
    current_num_particles = 0
    do i=1,buf_size
       if (convection_obj%local(i)%exists .eqv. .true.) &
            current_num_particles = current_num_particles + 1
    end do
  end function current_num_particles


  ! module subroutine create_particle(this, index)
  !   class(convection_exchangeable_t), intent(inout) :: this
  !   integer :: index
  !   integer :: a
  !   a = -1
  !   ! print *, "hi!"
  ! end subroutine create_particle
  module subroutine check_buf_size(i)
    integer, intent(in) :: i
    if (i .gt. particles_per_image) then
#if NO_COARRAYS
       print *, 1, ": ERROR put buffer overflow"
#else
       print *, this_image(), ": ERROR put buffer overflow"
#endif
       call exit
    end if
  end subroutine check_buf_size


  module subroutine initialize_from_file()
    integer :: total_parcels
    logical :: parcel_is_dry
    real    :: wind_speed
    namelist/parcel_parameters/ total_parcels, parcel_is_dry, wind_speed

    character(len=*), parameter :: file = 'parcel-parameters.txt'
    integer :: unit, rc, me, n_images
    logical :: exists

#if NO_COARRAYS
    me = 1
    n_images = 1
#else
    me = this_image()
    n_images = num_images()
#endif

    if (me .eq. 1) &
         print *, "Initializing convection parcels from file"
    inquire(file=file, exist=exists)

    if (exists .neqv. .true.) then
       if (me .eq. 1) &
            print*, trim(file), " does not exist, please create file"
       call exit
    end if

    unit = 10
    open(unit, file=file, action='read', status='old', iostat=rc)
    read(unit=unit, nml=parcel_parameters, iostat=rc)
    close(unit)

    particles_per_image = nint(total_parcels / real(n_images))
    local_buf_size = particles_per_image * 4
    dry_air_particles = parcel_is_dry
    input_wind = wind_speed
  end subroutine

  !-----------------------------------------------------------------
  ! p = p_0 * e^( ((9.81/287.058)*dz) / t_mean )
  !-----------------------------------------------------------------
  ! Method one: physics
  !        a) change pressure         b) update temperature
  !-----------------------------------------------------------------

  pure module subroutine dry_lapse_rate(pressure, temperature, potential_temp, &
          z_displacement)
    real, intent(inout) :: pressure, temperature, potential_temp
    real, intent(in) :: z_displacement

    real, parameter :: gravity = 9.80665
    real ::  p0, T_original
    p0 = pressure
    pressure = p0 - z_displacement * &
         gravity / (287.05 * temperature) * p0

    T_original = temperature
    temperature = exner_function(pressure) * potential_temp
  end subroutine
  ! if (temperature /= temperature) then
  !    print *, this_image(), ":: ~~~~~~NAN ERROR~~~~~~", &
  !         p0, pressure, z_displacement, T_original, temperature, &
  !         exner_function(pressure), potential_temp
  !    ! call exit
  !    return
  ! end if
  ! end procedure !pure dry_lapse_rate
  !-----------------------------------------------------------------
  ! Method two: adiabatic lapse rate, not using
  !-----------------------------------------------------------------
  ! if (particle%cloud_water .gt. 0.0) then
  !    gamma = 0.006 ! 6 celsius, not sure if this conversion is correct
  ! else
  !    gamma = 0.01 ! 10 celcius
  ! end if

end submodule
