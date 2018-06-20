program main
  use iso_fortran_env, only : input_unit
  use domain_interface, only : domain_t
  use assertions_interface, only : assert
  use module_mp_driver, only: microphysics
  use timer_interface, only: timer_t
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  implicit none

  integer :: num_ranks, rank, ierr
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
  print *, "  ________I AM RANK", rank, "__________"
  print *, "-------------------------------"
  rank = rank + 1
  if (rank==1) print *,"Number of images = ", num_ranks
  if (num_ranks == 1) then
    call exit
  endif
  call exit
  block
    type(domain_t), save :: domain
    integer :: i,nz, ypos,xpos
    type(timer_t) :: timer

    if (rank==1) print *,rank, &
                 "domain%initialize_from_file('input-parameters.txt')"
    call domain%initialize_from_file('input-parameters.txt')

    if (rank==1) then
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
    if (rank==1) print*, ""
    if (rank==1) print*, "Beginning simulation..."
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call timer%start()
    do i=1,20 !200
        ! note should this be wrapped into the domain object(?)
        call microphysics(domain, dt = 20.0, halo=1)
        call domain%halo_send()
        call microphysics(domain, dt = 20.0, subset=1)

        call domain%halo_retrieve()

        call domain%advect(dt = 1.0)

        ! if (this_image()==(num_images()/2)) then
        !     print*, domain%accumulated_precipitation(::3,ypos)
        ! endif

    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call timer%stop()


    if (rank==1) then
        print *,"Model run time:",timer%as_string('(f8.3," seconds")')
    endif

    ypos = (domain%jde-domain%jds)/2 + domain%jds
    do i=1,num_ranks
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        if (rank==i) then
            if ((ypos>=domain%jts).and.(ypos<=domain%jte)) then
                xpos = (domain%ite-domain%its)/2 + domain%its
                print*, rank, " : ", domain%accumulated_precipitation(domain%its:domain%ite:2,ypos)
            endif
        endif
    enddo
  end block

  if (rank==1) print *,"Test passed."
  call MPI_Finalize(ierr)
end program
