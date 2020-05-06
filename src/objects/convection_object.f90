module convection_object
  use convection_type_interface
  use grid_interface
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
     procedure :: process => simple_math
     procedure :: send
     procedure :: retrieve
     ! procedure :: exchange => exchange_data  ! Artless
  end type convection_object_t

contains
  subroutine simple_math(this, temperature)
    class(convection_object_t), intent(inout) :: this
    real, dimension(:,:,:), intent(in) :: temperature
    real, parameter :: gravity = 9.80665
    real, parameter :: artless_density = 1.003
    real :: a_prime, displacement, t, t_prime
    ! print*, " -- set -- "
    ! what type
    ! print *, "convection object type is ", this%convection_type
    ! print*, this%particle%u, this%particle%v, this%particle%w
    ! these are initialized in domain_implementation

    ! print *, "ARTLESS-----------", this%particle%temperature, &
    !     this%particle%pressure

    ! u,v,w wind field to control advection
    ! first just concerned with w, which is initially 0

    ! Parcel theory
    ! temp, pressure, density  :: T, p, rho

    !!              calculate acceleration equation
    !!! a_prime_z =  -1/rho_prime * partial_p/partial_z - g
    !!
    !!! partial_p_partial_z =
    ! a_prime = -(1 / artless_density) * partial_p_partial_z - gravity
    ! print *, "~~~~~~~ a_prime = ", a_prime

    !!!
    !!! a'_z = (T' - T) / T * g
    !!!
    T = temperature(this%particle%x, this%particle%y, this%particle%z)
    T_prime =this%particle%temperature
    a_prime = (T_prime - T) / T * gravity
    ! print *, "----------- variables -------------"
    ! print *, T, T_prime, a_prime

    !!! move particle
    !!! displacement
    !!! d = v_0*t + 1/2*a*t^2
    displacement =  0        + 0.5 * a_prime * 1 * 1
    ! print *, "~~~~~~~ displacement = ", displacement
    ! print *, " to be added to ", this%particle%k
    this%particle%z = this%particle%z + displacement
    print *, "z = ", this%particle%z  , "with displacement", displacement
    ! print*, " -- fin -- "
  end subroutine simple_math


  ! Choose type of initialization based on type of particle
  subroutine init(this, convection_type_enum, grid, u_in, v_in, w_in, &
      temperature, pressure)
    use convection_type_interface
    class(convection_object_t), intent(inout) :: this
    type(grid_t) :: grid
    integer(c_int), intent(in) :: convection_type_enum
    real, optional, intent(in) :: u_in,v_in,w_in
    real, dimension(:,:,:) :: temperature, pressure
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
      ! ! ARTLESS :: TODO
      ! ! NOW STARTING IN CENTER, FUTURE RANDOM
      ! this%particle = convection_particle(.true., &
      !     (grid%ims + grid%ime) / 2.0, &
      !     (grid%jms + grid%jme) / 2.0, &
      !     (grid%kms + grid%kme) / 2.0, &
      !     u,v,w, &
      !     -0.0,-0.0)

      ! ! ---- calculate pressure
      ! ! ARTLESS :: better way to get cordinates
      ! this%particle%temperature = temperature(floor(this%particle%x), &
      !     floor(this%particle%z), &
      !     floor(this%particle%y))
      ! this%particle%pressure = pressure(floor(this%particle%x), &
      !     floor(this%particle%z), &
      !     floor(this%particle%y))

      ! ! print *, "------ RANGE ----- ", &
      ! !     temperature(floor(this%particle%i), &
      ! !     floor(this%particle%k), &
      ! !     floor(this%particle%j)), &
      ! !     temperature(ceiling(this%particle%i), &
      ! !     ceiling(this%particle%k), &
      ! !     ceiling(this%particle%j))

      ! print *, "!!! I'm add 10k to temperature !!! For Testing"
      ! this%particle%temperature = this%particle%temperature + 10

      ! ! print *, "====", temperature(5,:,5)
      ! ! print *, "SHAPE OF ", shape(temperature)
      ! ! print *, "ARTLESS--- TEMPERATURE ----------", this%particle%temperature
      ! !     this%particle%pressure

    case(convection_array_e)
      ! this%array = convection_array(u,v,w)
    case(convection_linked_list_e)
    end select
  end subroutine init

  subroutine send(this)
    use convection_type_interface
    class(convection_object_t), intent(in) :: this
    print *, "---sending---"

  end subroutine send

  subroutine retrieve(this)
    use convection_type_interface
    class(convection_object_t), intent(in) :: this
    print *, "---retrieving---"
  end subroutine retrieve

end module convection_object
