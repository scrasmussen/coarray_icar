program main
  use iso_fortran_env, only : input_unit
  use domain_interface, only : domain_t
  use assertions_interface, only : assert
  use module_mp_driver, only: microphysics
  use timer_interface, only: timer_t
  use convection_exchangeable_interface, only : num_particles
  implicit none

  integer :: me
  me = this_image()
  if (me==1) print *,"Number of images = ",num_images()

  block
    type(domain_t), save :: domain

    ! parameters to setup test
    integer, parameter :: timesteps = 500
    logical            :: report = .true.
    logical, parameter :: convection = .true.
    logical, parameter :: print_timestep = .false.

    integer :: i,nz, ypos,xpos, n_particles
    type(timer_t) :: timer
    character*32 hostname

    if (convection .eqv. .false.) report = .false.

    ! call get_environment_variable('HOSTNAME',hostname)
    ! print *, me, ': has hostname ', trim(hostname)
    ! sync all

    if (me==1) print *,me,"domain%initialize_from_file('input-parameters.txt')"
    call domain%initialize_from_file('input-parameters.txt', convection)

    if (me==1) then
        nz = size(domain%pressure,2)
        print *, " Layer height       Pressure        Temperature      Water Vapor"
        print *, "     [m]              [hPa]             [K]            [kg/kg]"
        do i=nz,1,-4
            print *,domain%z(1,i,1), domain%pressure(1,i,1)/100, domain%temperature(1,i,1), domain%water_vapor%local(1,i,1)
        end do
    endif

    ypos = (ubound(domain%accumulated_precipitation,2)-lbound(domain%accumulated_precipitation,2))/2
    ypos = ypos + lbound(domain%accumulated_precipitation,2)

    ! initialize microphysics before starting the timer
    call microphysics(domain, dt = 20.0)

    if (me==1) print*, ""
    if (me==1) print*, "Beginning simulation..."

    if (report .eqv. .true.) call domain%report_convection(0)


    sync all
    call timer%start()
    do i=1, timesteps
        if ((print_timestep .eqv. .true.) .and. me==1) then
          print *, &
              "________________________timestep ", i,"________________________"
        end if
        ! note should this be wrapped into the domain object(?)
        call microphysics(domain, dt = 20.0, halo=1, &
            convected_particles = convection)
        call domain%halo_send()
        call microphysics(domain, dt = 20.0, subset=1, &
            convected_particles = convection)
        ! call domain%halo_send()
        call domain%halo_retrieve(convection)

        if (report .eqv. .true.) call domain%report_convection(i)

        call domain%advect(dt = 1.0)
    end do
    sync all
    call timer%stop()

    if (me==1) then
        n_particles = num_particles()
        print *, "For", timesteps, "timesteps"
        if (convection .eqv. .true.) then
            print *, "With", n_particles, "particles per image for a total of", &
                n_particles * num_images()
        else
            print *, "With no particles"
        end if
        print *,"Model run time:",timer%as_string('(f8.3," seconds")')
    endif

    ypos = (domain%jde-domain%jds)/2 + domain%jds
    do i=1,num_images()
        sync all
        if (me==i) then
            if ((ypos>=domain%jts).and.(ypos<=domain%jte)) then
                xpos = (domain%ite-domain%its)/2 + domain%its
                ! print*, me, " : ", domain%accumulated_precipitation(domain%its:domain%ite:2,ypos)
            endif
        endif
    enddo

  end block


  if (me==1) print *,"Test passed."

end program
