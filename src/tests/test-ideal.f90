program main
  use iso_fortran_env, only : input_unit
  use domain_interface, only : domain_t
  use assertions_interface, only : assert
  use module_mp_driver, only: microphysics
  use timer_interface, only: timer_t
  use convection_exchangeable_interface, only : num_particles, &
       are_particles_dry, num_particles_communicated, get_wind_speed, &
       current_num_particles
  implicit none

  integer :: me, ierrr
  me = this_image()
  ! call TAU_PROFILE_SET_NODE(me)
#if _CRAYFTN
  call assign('assign -S on -y on p:%.txt', ierrr)
#endif

  if (me==1) print *,"Number of images = ",num_images()

  block
    type(domain_t), save :: domain

    ! parameters to setup test
    integer, parameter :: timesteps = 200
    logical            :: report = .false.
    logical, parameter :: use_sounding   = .false.
    logical, parameter :: print_timestep = .true.

    integer :: i,nz, ypos,xpos, n_particles, particles_communicated, p
    integer, allocatable :: current_n_particles[:]
    type(timer_t) :: timer
    logical :: exist
    character(len=32) :: filename
    integer :: len, ierr
    character*32 processorname

    ! call get_environment_variable('HOSTNAME',hostname)
    ! call MPI_Get_processor_name( processorname, len, ierr );
    ! print *, me, ': has processor name ', trim(processorname)
    ! sync all
    ! call exit

    if (me==1) print *,me,"domain%initialize_from_file('input-parameters.txt')"
    call domain%initialize_from_file('input-parameters.txt', use_sounding)

    if (me==1) then
        nz = size(domain%pressure,2)
        print *, " Layer height       Pressure        Temperature      Water Vapor"
        print *, "     [m]              [hPa]             [K]            [kg/kg]"
        do i=nz,1,-4
            print *,domain%z(1,i,1), domain%pressure(1,i,1)/100, &
                  domain%temperature(1,i,1), domain%water_vapor%local(1,i,1)
        end do
    endif

    ypos = (ubound(domain%accumulated_precipitation,2) - &
          lbound(domain%accumulated_precipitation,2))/2
    ypos = ypos + lbound(domain%accumulated_precipitation,2)

    ! initialize microphysics before starting the timer
    if (me==1) print*, "Initializing microphysics..."
    call microphysics(domain, dt = 20.0)

    if (me==1) print*, ""
    if (me==1) print*, "Beginning simulation..."

    n_particles = num_particles()
    if (n_particles .eq. 0) report = .false.
    if (report .eqv. .true.) call domain%report_convection(0)

    sync all
    call timer%start()
    do i=1, timesteps
        if ((print_timestep .eqv. .true.) .and. me==1) then
            print *, &
                  "_____________________timestep ", i,"________________________"
        end if
        ! note should this be wrapped into the domain object(?)
        call microphysics(domain, dt = 20.0, halo=1)

        call domain%halo_send()
        call microphysics(domain, dt = 20.0, subset=1, t=i,&
              report_convection = report)
        ! call domain%halo_send()
        call domain%halo_retrieve()

        call domain%advect(dt = 1.0)
    end do
    sync all
    call timer%stop()

    particles_communicated = num_particles_communicated()
    allocate(current_n_particles[*])
    current_n_particles = current_num_particles(domain%convection_obj)
    call co_sum(current_n_particles)

    if (me==1) then
        print *, "For", timesteps, "timesteps"
        if (n_particles .gt. 0) then
            print *, "With", n_particles, &
                  "particles per image for a total of", &
                  n_particles * num_images()
            if (current_n_particles .ne. (n_particles * num_images())) &
               print *, "ERROR: final n particles .ne. orignal n particles",&
               current_n_particles, "vs.", n_particles
        else
            print *, "With no particles"
        end if
        print *,"Model run time:",timer%as_string('(f8.3," seconds")')
        print *,"Model get_time():", timer%get_time()
        me = this_image()

        ! handle opening of file and report
        write (filename,"(A18)") "timing_results.txt"
        inquire(file=filename, exist=exist)
        if (exist) then
           open(unit=me, file=filename, status='old', position='append')
        else
           open(unit=me, file=filename, status='new')
        end if
        write(me,*) domain%nx_global, domain%nz, domain%ny_global, &
             num_images(), domain%ximages, domain%yimages, &
             n_particles, timesteps, timer%get_time(), &
             are_particles_dry(), get_wind_speed(), &
             particles_communicated
        close(me)
    endif

    do i=1,num_images()
        sync all
        if (me==i) then
           write (filename,"(A13)") "p_results.txt"
           inquire(file=filename, exist=exist)
           if (exist) then
              open(unit=me, file=filename, status='old', position='append')
           else
              open(unit=me, file=filename, status='new')
           end if

           do p=1,size(domain%convection_obj%local)
              if (domain%convection_obj%local(p)%exists .eqv. .true.) then
                 write(me,*) domain%convection_obj%local(p)%particle_id, &
                      domain%convection_obj%local(p)%moved
              end if
           end do
        end if
     end do


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

end program main
