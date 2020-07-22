submodule(convection_exchangeable_interface) &
    convection_exchangeable_implementation
  use assertions_interface, only : assert, assertions
  use exchangeable_interface, only : exchangeable_t
  use domain_interface, only : pressure_at_elevation, exner_function, sat_mr

  use grid_interface, only : grid_t
  implicit none

  ! ----- PARAMETERS TO TUNE CONVECTION MODEL -----
  logical, parameter :: wrap_neighbors = .true.
  logical, parameter :: advection = .true.
  logical, parameter :: wind = .true.
  logical, parameter :: relative_humidity = .true.
  logical, parameter :: caf_comm_message = .false.
  logical, parameter :: particle_create_message = .false.
  integer, parameter :: particles_per_image = 30
  integer, parameter :: local_buf_size = particles_per_image * 4

  integer, parameter :: default_buf_size=1
  integer, parameter :: default_halo_size=1
  integer, save, allocatable :: neighbors(:)
  integer, save :: north_neighbor, south_neighbor, buf_size, halo_size
  integer, save :: east_neighbor, west_neighbor
  integer, save :: northeast_neighbor, northwest_neighbor
  integer, save :: southeast_neighbor, southwest_neighbor

contains
  module subroutine const(this, potential_temp, u_in, v_in, w_in, grid, z_m, &
      z_interface, ims, ime, kms, kme, jms, jme, dz_value, &
      its, ite, kts, kte, jts, jte, input_buf_size, halo_width)
    class(convection_exchangeable_t), intent(inout) :: this
    class(exchangeable_t), intent(in)    :: potential_temp
    class(exchangeable_t), intent(in)    :: u_in, v_in, w_in
    type(grid_t), intent(in)      :: grid
    real, intent(in)              :: z_m(ims:ime,kms:kme,jms:jme)
    real, intent(in)              :: z_interface(ims:ime,jms:jme)
    integer, intent(in)           :: ims, ime, kms, kme, jms, jme
    integer, intent(in)           :: its, ite, kts, kte, jts, jte
    real, intent(in)              :: dz_value
    integer, intent(in), optional :: input_buf_size
    integer, intent(in), optional :: halo_width

    integer :: me, create
    real :: random_start(3), x, z, y
    real :: z_meters, z_interface_val, theta_val
    real :: pressure_val, exner_val, temp_val, water_vapor_val
    real :: u_val, v_val, w_val, rand
    integer :: x0, x1, z0, z1, y0, y1
    logical :: calc
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

    allocate(this%local(particles_per_image * 4))
    if (particle_create_message .eqv. .true.) then
      if (me == 1) then
        print*, "Creating", particles_per_image, "parcels per image"
      end if
    end if

    do create=1,particles_per_image
      call random_number(random_start)
      x = its + (random_start(1) * (ite-its))
      z = kts + (random_start(3) * (kte-kts))
      y = jts + (random_start(2) * (jte-jts))

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

      associate (A => z_interface)
        z_interface_val = bilinear_interpolation(x, x0, x1, y, y0, y1, &
            A(x0,y0), A(x0,y1), A(x1,y0), A(x1,y1))
      end associate

      ! sealevel_pressure => 100000.0
      pressure_val = pressure_at_elevation(100000.0, z_meters)
      exner_val = exner_function(pressure_val)

      call random_number(rand)
      theta_val = theta_val * (1 + rand / 10) ! random 1-10% change
      temp_val = exner_val * theta_val

      water_vapor_val = sat_mr(temp_val, pressure_val)

      ! wind is constant in this system, ignoring w wind (aka z-direction)
      u_val = u_in%local(x0, z0, y0)
      v_val = v_in%local(x0, z0, y0)
      w_val = 0

      call this%create_particle_id()
      this%local(create) = convection_particle(this%particle_id_count, .true., &
          .false., x, y, z, u_val, v_val, w_val, z_meters, z_interface_val, &
          pressure_val, temp_val, theta_val, 0, water_vapor_val, 0)
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

    call this%setup_neighbors(grid)
  end subroutine

  module subroutine process(this, nx_global, ny_global, &
      ims, ime, kms, kme, jms, jme, dt, dz, temperature, z_interface, &
      its, ite, kts, kte, jts, jte)
    implicit none
    class(convection_exchangeable_t), intent(inout) :: this
    integer, intent(in) :: nx_global, ny_global
    real, intent(in)    :: dt, dz
    real, intent(in) :: temperature(ims:ime,kms:kme,jms:jme)
    real, intent(in) :: z_interface(ims:ime,jms:jme)
    integer, intent(in) :: ims, ime, kms, kme, jms, jme
    integer, intent(in) :: its, ite, kts, kte, jts, jte
    real, parameter :: gravity = 9.80665
    real :: Gamma
    real :: a_prime, z_displacement, t, t_prime, buoyancy
    real :: ws, fake_wind_correction, delta_z, z_interface_val, z_wind_change
    integer :: i,j,k, l_bound(1), dif(1), new_ijk(3), me
    real :: new_pressure, R_s, p0, exponent, alt_pressure, alt_pressure2
    real :: alt_pressure3, alt_pressure4, mixing_ratio, sat_mr_val
    real :: vapor_p, sat_p, T_C, T_K, mr, T_squared, T_original, tmp
    real :: xx, yy
    real :: rate_of_temp_change
    integer :: x0, x1, z0, z1, y0, y1
    integer :: u_bound(3)
    logical :: calc, calc_x, calc_y

    me = this_image()

    do i=1,ubound(this%local,1)
      associate (particle=>this%local(i))
        if (particle%exists .eqv. .true.) then
          if (particle%x .lt. its-1 .or. &
              particle%z .lt. kts-1 .or. &
              particle%y .lt. jts-1 .or. &
              particle%x .gt. ite+1 .or. &
              particle%z .gt. kte+1 .or. &
              particle%y .gt. jte+1) then
            print *, "particle", particle%particle_id, "on image", me
            print *, "x:", its, "<", particle%x, "<", ite, "with halo 2"
            print *, "z:", kts, "<", particle%z, "<", kte, "with halo 2"
            print *, "y:", jts, "<", particle%y, "<", jte, "with halo 2"
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
          ! new: properites of environment are prime
          T = particle%temperature

          associate (A => temperature, x => particle%x, z => particle%z , &
              y => particle%y)
            x0 = floor(x); z0 = floor(z); y0 = floor(y);
            x1 = ceiling(x); z1 = ceiling(z); y1 = ceiling(y);

            T_prime = trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0,y1,&
                A(x0,z0,y0), A(x0,z0,y1), A(x0,z1,y0), A(x1,z0,y0), &
                A(x0,z1,y1), A(x1,z0,y1), A(x1,z1,y0), A(x1,z1,y1))
          end associate

          if (advection .eqv. .true.) then
            buoyancy = (T - T_prime) / T_prime
            a_prime = buoyancy * gravity
            ! time step is 1 so t and t^2 go away
            z_displacement = particle%velocity + 0.5 * a_prime
            ! number from dz_interface, currently always 500
            delta_z = z_displacement / dz
            if (z_displacement /= z_displacement) then
              ! when T=0 it causes division by 0. This problem should be fixed
              print *, me, ":: ------------NAN--------------", &
                  particle%particle_id, T, T_prime
              particle%z = -1
              print *, p0, z_displacement, gravity, particle%temperature
              cycle
              ! call flush()
              ! call exit
            else
              particle%z = particle%z + delta_z
              particle%velocity = z_displacement
            end if
          else
            z_displacement = 0.0
          end if

          ! print *, "delta_z ::", delta_z, "z_displacement",z_displacement
          ! print *, me, ":", particle%particle_id, z_displacement
          !-----------------------------------------------------------------
          ! Orographic lift
          ! Find dz change from change in environment
          !-----------------------------------------------------------------
          if (wind .eqv. .true.) then
            ! print *, "old wind", particle%x, particle%y, "z", particle%z_meters
            fake_wind_correction = (1.0 / dz) ! real correction
            ! fake_wind_correction = (1.0) ! correction
            ! u: zonal velocity, wind towards the east
            particle%x = particle%x + (particle%u * fake_wind_correction)
            ! v: meridional velocity, wind towards north
            particle%y = particle%y + (particle%v * fake_wind_correction)
            associate (A => z_interface, x=> particle%x, y=>particle%y)
              x0 = floor(x); x1 = ceiling(x)
              y0 = floor(y); y1 = ceiling(y)
              z_interface_val = bilinear_interpolation(x, x0, x1, y, y0, y1, &
                  A(x0,y0), A(x0,y1), A(x1,y0), A(x1,y1))
            end associate
            z_wind_change = z_interface_val - particle%z_interface
            particle%z_interface = z_interface_val
            z_displacement = z_displacement + z_wind_change
          end if

          !-----------------------------------------------------------------
          ! Move particle, remove particle if beyond the z axis
          !-----------------------------------------------------------------
          particle%z_meters = particle%z_meters + z_displacement
          if (particle%z .lt. kts .or. particle%z .gt. kte) then
            particle%exists=.false.
            cycle
          end if


          !-----------------------------------------------------------------
          ! p = p_0 * e^( ((9.81/287.058)*dz) / t_mean )
          !-----------------------------------------------------------------
          ! Method one, change pressure, update temperature
          !-----------------------------------------------------------------
          if (1 .eq. 1) then
            p0 = particle%pressure
            particle%pressure = p0 - z_displacement * &
                gravity / (287.05 * particle%temperature) * p0

            T_original = particle%temperature
            particle%temperature = exner_function(particle%pressure) * &
                particle%potential_temp

            if (particle%temperature /= particle%temperature) then
              print *, me, ":: ~~~~~~NAN~~~~~~", &
                  particle%particle_id, particle%temperature, &
                  exner_function(particle%pressure), particle%pressure, &
                  particle%potential_temp
              print*, p0, z_displacement, gravity ,287.05,particle%temperature
              particle%exists = .false.
              ! call flush()
              ! call exit
            end if
          end if
          !-----------------------------------------------------------------
          ! Method two: heuristics, not using
          !-----------------------------------------------------------------
          if (0 .eq. 1) then
            if (particle%cloud_water .gt. 0.0) then
              gamma = 0.006 ! 6 celsius, not sure if this conversion is correct
            else
              gamma = 0.01 ! 10 celcius
            end if
          end if


          !-----------------------------------------------------------------
          ! Mixing Ratio
          ! saturate = sat_mr(t,p)
          ! if (water_vapor > saturated)
          !   condensate = water_vapor - satured
          !   water_vapor -= condensate
          !   clouds += condensate
          !   temperature += latent_heat * condensate
          !-----------------------------------------------------------------
          if (relative_humidity .eqv. .true.) then
            block
              real :: saturate, condensate, vapor, vapor_needed, RH
              ! latent heat values, calculating using formula
              real, parameter :: condensation_lh = 2600!000 ! 2.5 x 10^6 J/kg
              real, parameter :: vaporization_lh = -condensation_lh
              ! specific heat of water vapor at constant volume
              ! real, parameter :: C_vv = 1.0 / 1390.0
              real :: C_vv, c_p
              real :: latent_heat, Q_heat, temp_c, delta_t, T0
              integer :: repeat

              saturate = sat_mr(particle%temperature, particle%pressure)
              RH = particle%water_vapor / saturate

              if (particle%cloud_water .gt. 0.0 .and. RH .lt. 1.0) then
                vapor_needed = 1.0 - RH
                if (vapor_needed > particle%cloud_water) then
                  vapor = particle%cloud_water
                  particle%cloud_water = 0
                else
                  vapor = vapor_needed
                  particle%cloud_water = particle%cloud_water - vapor_needed
                end if

                particle%water_vapor = particle%water_vapor + vapor
                latent_heat = 2.5 * 10**6 ! J kg^-1  for condensation
                Q_heat = latent_heat * vapor ! kJ
                c_p = (1004 * (1 + 1.84 * vapor)) ! 3.3
                delta_t = Q_heat / c_p   ! 3.2c
                particle%temperature = particle%temperature - delta_t
                ! print *, "VAPOR WAS NEEDED"

              else if (RH .ge. 1.0) then    ! if (0.0001 > saturate) then
                condensate = particle%water_vapor - saturate
                particle%water_vapor = particle%water_vapor - condensate
                particle%cloud_water = particle%cloud_water + condensate

                T0 = particle%temperature
                p0 = particle%pressure

                !--------------------------------------------------------------
                ! calculate specific latent heat
                !--------------------------------------------------------------
                ! if (T0 .lt. 248.15) then
                !   latent_heat = 2600
                ! else if (T0 .gt. 314.15) then
                !   latent_heat = 2400
                ! else
                !   T_C = T0 - 273.15
                !   latent_heat = 2500.8 - 2.36*T_C + 0.0016 * T_C**2 - &
                !       0.00006 * T_C**3
                ! end if

                ! specific latent heat for condensation
                latent_heat = 2.5 * 10**6 ! J kg^-1
                ! https://en.wikipedia.org/wiki/Latent_heat#Specific_latent_heat
                Q_heat = latent_heat * condensate ! kJ


                !--------------------------------------------------------------
                ! calculate specific heat
                !--------------------------------------------------------------
                ! Stull: Practical Meteorology
                c_p = (1004 * (1 + 1.84 * condensate)) ! 3.3
                delta_t = Q_heat / c_p   ! 3.2c
                particle%temperature = T0 + delta_t


                ! -- is pressure constant during this process?
                ! using Poisson's equation
                ! particle%pressure = p0/ ((T0/particle%temperature)**(1/0.286))

                ! call flush()
                ! call exit

              else
                exit
              end if
            end block
          end if ! ---- end of relative humidity seciton ----



          !-----------------------------------------------------------------
          ! Move particle if needed
          !-----------------------------------------------------------------
          associate (x => particle%x, y => particle%y, z => particle%z)
            if (x .lt. its .or. x .gt. ite .or. &
                z .lt. kms .or. z .gt. kme .or. &
                y .lt. jts .or. y .gt. jte) then
              if (caf_comm_message .eqv. .true.) then
                print *, "PUTTING", particle%x, particle%y, particle%z_meters, &
                    "FROM", this_image(), "id:", particle%particle_id, &
                    "M", ims,ime, kms, kme, jms, jme, "T", its,ite,jts,jte
                print *, "------"
              end if
            end if

            xx = x
            yy = y
            ! If particle is getting wrapped the x and y values need to be
            ! properly updated
            if (wrap_neighbors .eqv. .true.) then
              if (x > nx_global) then
                if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                x = x - nx_global + 1
              else if (x < 1) then
                if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                x = x + nx_global - 1
              end if

              if (y > ny_global) then
                if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                y = y - ny_global + 1
              else if (y < 1) then
                if (caf_comm_message .eqv. .true.) print *, "WRAPPED"
                y = y + ny_global - 1
              end if
            end if

            ! Check values to know where to send particle
            if (yy .lt. jts) then      ! "jts  <   y    <  jte"
              if (xx .lt. its) then
                call this%put_southwest(particle)
              else if (xx .gt. ite) then
                call this%put_southeast(particle)
              else
                call this%put_south(particle)
              end if
            else if (yy .gt. jte) then ! jts will be 1
              if (xx .lt. its) then
                call this%put_northwest(particle)
              else if (xx .gt. ite) then
                call this%put_northeast(particle)
              else
                call this%put_north(particle)
              endif
            else if (xx .lt. its) then ! "its  <   x    <  ite"
              call this%put_west(particle)      ! need to double check this!
            else if (xx .gt. ite) then
              call this%put_east(particle)
            end if
          end associate
        end if
      end associate
    end do
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
      end associate
    end if

  end subroutine

  function bilinear_interpolation(x, x0, x1, y, y0, y1, &
      c00, c01, c10, c11) result(c)
    real, intent(in) :: x, y
    real, intent(in) :: c00, c01, c10, c11
    integer, intent(in) :: x0, x1, y0, y1
    real :: xd, yd, c00, c01, c10, c11, c0, c1, c
    if (x0 .eq. x1) then
      c = c00 + (c11-c00) * (y-y0) / (y1-y0)
    else if (y0 .eq. y1) then
      c = c00 + (c11-c00) * (x-x0) / (x1-x0)
    else
      c0 = ((x1 - x)/(x1 - x0)) * c00 + ((x - x0)/(x1 - x0)) * c10
      c1 = ((x1 - x)/(x1 - x0)) * c01 + ((x - x0)/(x1 - x0)) * c11
      c = ((y1 - y)/(y1 - y0)) * c0 + ((y - y0)/(y1 - y0)) * c1
    end if
  end function bilinear_interpolation

  function trilinear_interpolation(x, x0, x1, z, z0, z1, y, y0, y1, &
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
  end function trilinear_interpolation

  function num_particles()
    integer :: num_particles
    num_particles = particles_per_image
  end function num_particles
end submodule
