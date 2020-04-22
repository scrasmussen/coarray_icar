module convection_object
  use convection_exchangeable_interface
  use convection_type_interface
  use iso_c_binding, only: c_int
  implicit none
  private

  type, public :: convection_object_t
     integer(c_int) :: convection_type
     type(convection_particle) :: particle
     type(convection_array) :: array
     type(convection_list) :: list
     ! type(convection_exchangeable_t) :: exchangeable
     ! type(convection_exchangeable_particle_t) :: particle
   contains
     procedure :: initialize => init
  end type convection_object_t

contains

  subroutine init(this, convection_type_enum,u_in,v_in,w_in)
    use convection_type_interface
    use iso_c_binding, only : c_int
    class(convection_object_t), intent(inout) :: this
    integer(c_int), intent(in) :: convection_type_enum
    real, optional, intent(in) :: u_in,v_in,w_in
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

    this%convection_type = convection_type_enum
    select case (convection_type_enum)
    case(convection_particle_e)
      this%particle = convection_particle(u,v,w)
    case(convection_array_e)
      ! this%particle = convection_particle(u,v,w)
    case(convection_linked_list_e)
    end select
  end subroutine init

end module convection_object
