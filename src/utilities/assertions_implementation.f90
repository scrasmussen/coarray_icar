submodule(assertions_interface) assertions_implementation
  implicit none
contains
  elemental impure module subroutine assert(assertion,description,success)
    use iso_fortran_env, only : error_unit
    use mpi_f08, only : MPI_Comm_rank, MPI_COMM_WORLD
    logical, intent(in) :: assertion
    character(len=*), intent(in) :: description
    logical, intent(out), optional :: success
    integer :: rank, ierr
    if (present(success)) success=assertion
    if (.not.assertion) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      rank = rank + 1

      write(error_unit,*) 'Assertion "',description,'" failed on image ', rank
      if (.not. present(success)) error stop
    end if
  end subroutine
end submodule
