program main
  use omp_lib
  use iso_fortran_env, only : input_unit
  use domain_interface, only : domain_t
  use assertions_interface, only : assert
  use module_mp_driver, only: microphysics
  use timer_interface, only: timer_t
  use convection_exchangeable_interface, only : total_num_particles, &
       are_particles_dry, num_particles_communicated, get_wind_speed, &
       current_num_particles, num_particles_per_image
  implicit none

  integer :: me, n_images, ierrr
#if NO_COARRAYS
  me = 1
  n_images = 1
#else
  me = this_image()
  n_images = num_images()
#endif
  ! call TAU_PROFILE_SET_NODE(me)
#if _CRAYFTN
  call assign('assign -S on -y on p:%.txt', ierrr)
#endif

  if (me==1) print *,"Number of images = ",n_images

  block
    type(domain_t), save :: domain

    ! parameters to setup test
    logical, parameter :: just_particles = .true.
    integer, parameter :: timesteps = 50! 200
    logical            :: report = .false.
    logical, parameter :: use_sounding   = .false.
    logical, parameter :: print_timestep = .false.
    logical, parameter :: count_p_comm = .false.
    ! if count_p_comm is true, need to edit other convect_exhange_implementation
    logical, parameter :: save_particles_moved = .false.

    integer :: i,nz, ypos,xpos, particles_communicated, p
#if NO_COARRAYS
    integer, allocatable :: current_n_particles
#else
    integer, allocatable :: current_n_particles[:]
#endif
    type(timer_t) :: timer
    logical :: exist
    character(len=32) :: filename
    integer :: len, ierr
    integer :: dz_lb(3), ii
    real    :: dt


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
    dt = 20.0

    if (me==1) print*, ""
    if (me==1) print*, "Beginning simulation..."
#ifdef DC_LOOP
    if (me==1) print *, "--- do concurrent ---"
    print *, "max threads =", omp_get_max_threads()
#endif
#ifdef OMP_LOOP
    if (me==1) then
       print *, "--- OpenMP ---"
       print *, "max threads =", omp_get_max_threads()
    end if
#endif
#ifdef DO_LOOP
    if (me==1) print *, "--- do ---"
#endif


    if (total_num_particles() .eq. 0) report = .false.
    if (report .eqv. .true.) call domain%report_convection(0)

    ! print *, me,":", domain%nx, domain%ny, domain%nx_global, domain%ny_global

#ifndef NO_COARRAYS
    sync all
#endif
    call timer%start()

    do i=1, timesteps
        if (just_particles .eqv. .true.) then
           do ii=1,int(dt)
              if ((print_timestep .eqv. .true.) .and. me==1) then
                 print *, &
                      "_____________________timestep ", &
                      int((i-1)*dt + ii), "________________________"
                 ! print *, me,":", current_num_particles(domain%convection_obj)
              end if

              dz_lb = lbound(domain%dz_interface)
              call domain%convection_obj%process( &
                   domain%nx_global, domain%ny_global, &
                   domain%ims, domain%ime, domain%kms, &
                   domain%kme, domain%jms, domain%jme, &
                   1.0, domain%dz_interface(dz_lb(1),dz_lb(2),dz_lb(3)), &
                   domain%temperature, domain%z_interface(:,1,:), &
                   domain%its, domain%ite, domain%kts, domain%kte, domain%jts, &
                   domain%jte, domain%z, domain%potential_temperature, &
                   domain%pressure, domain%u, domain%v, domain%w, &
                   int((i-1)*dt + ii) )
              if (report .eqv. .true.) &
                   call domain%report_convection(int((i-1)*dt + ii))
              call domain%convection_obj%retrieve(no_sync=.false.)
           end do
           cycle
        end if

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

#ifndef NO_COARRAYS
    sync all
#endif
    call timer%stop()

    if (me==1) print *, "Finished"


    if (count_p_comm .eqv. .true.) then
       particles_communicated = num_particles_communicated()
#if NO_COARRAYS
       allocate(current_n_particles)
#else
       allocate(current_n_particles[*])
#endif
       current_n_particles = current_num_particles(domain%convection_obj)
#ifndef NO_COARRAYS
       call co_sum(current_n_particles)
#endif
    else
       particles_communicated = -1
    end if

    if (me==1) then
        print *, "For", timesteps, "timesteps"
        if (total_num_particles() .gt. 0) then
            print *, "With", num_particles_per_image(), &
                  "particles per image for a total of", &
                  total_num_particles()
        else
            print *, "With no particles"
        end if
        print *,"Model run time:",timer%as_string('(f8.3," seconds")')
        print *,"Model get_time():", timer%get_time()

        ! handle opening of file and report
        write (filename,"(A18)") "timing_results.txt"
        inquire(file=filename, exist=exist)
        if (exist) then
           open(unit=me, file=filename, status='old', position='append')
        else
           open(unit=me, file=filename, status='new')
        end if
        write(me,*) ceiling(n_images/44.0), &
             domain%nx_global, domain%nz, domain%ny_global, &
             n_images, domain%ximages, domain%yimages, &
             num_particles_per_image(), timesteps, timer%get_time(), &
             are_particles_dry(), get_wind_speed(), &
#ifdef OMP_LOOP
             particles_communicated, "omp", omp_get_max_threads()
#endif
#ifdef DC_LOOP
             particles_communicated, 'dc', omp_get_max_threads()
#endif
#ifdef DO_LOOP
             particles_communicated, "do", 0

#endif
        close(me)
    endif


    if (save_particles_moved .eqv. .true.) then
       do i=1,n_images
#ifndef NO_COARRAYS
          sync all
#endif
          if (me==i) then
             write (filename,"(A19)") "particles-moved.txt"
             inquire(file=filename, exist=exist)
             if (exist) then
                open(unit=me, file=filename, status='old', position='append')
             else
                open(unit=me, file=filename, status='new')
             end if
             do p=1,size(domain%convection_obj%local)
                if (domain%convection_obj%local(p)%exists .eqv. .true.) then
                   write(me,*) domain%convection_obj%local(p)%moved
                end if
             end do
             close(unit=me)
          end if
       end do
    end if

    ypos = (domain%jde-domain%jds)/2 + domain%jds
    do i=1,n_images
#ifndef NO_COARRAYS
        sync all
#endif
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
