module domain_interface
  use configuration_interface, only : configuration_t
  use exchangeable_interface, only : exchangeable_t
  use grid_interface, only: grid_t
  ! use convection_object
  use convection_exchangeable_interface, only : &
      convection_exchangeable_t
  implicit none

  private
  public :: domain_t
  public :: pressure_at_elevation, exner_function, sat_mr

  type domain_t
    ! private
    ! core model species to be advected
    type(exchangeable_t) :: water_vapor
    type(exchangeable_t) :: potential_temperature
    type(exchangeable_t) :: cloud_water_mass
    type(exchangeable_t) :: cloud_ice_mass
    type(exchangeable_t) :: cloud_ice_number
    type(exchangeable_t) :: rain_mass
    type(exchangeable_t) :: rain_number
    type(exchangeable_t) :: snow_mass
    type(exchangeable_t) :: graupel_mass
    ! type(convection_object_t) :: convection_obj
    type(convection_exchangeable_t) :: convection_obj

    ! core model variables (not advected)
    real, public, allocatable :: exner(:,:,:)
    real, public, allocatable :: pressure(:,:,:)
    real, public, allocatable :: temperature(:,:,:)
    real, public, allocatable :: z(:,:,:)
    real, public, allocatable :: dz_interface(:,:,:)
    real, allocatable :: z_interface(:,:,:)
    real, allocatable :: dz_mass(:,:,:)
    real, allocatable :: accumulated_precipitation(:,:)
    real, allocatable :: accumulated_snowfall(:,:)

    ! these coarrays are used to send all data back to a master image for IO... feels very inefficient
    real, public, allocatable :: transfer_array_3d(:,:,:)[:]
    real, public, allocatable :: transfer_array_2d(:,:)[:]

    ! wind field to control advection
    type(exchangeable_t) :: u
    type(exchangeable_t) :: v
    type(exchangeable_t) :: w

    ! contains the size of the domain (or the local tile?)
    integer :: nx, ny, nz, nx_global, ny_global
    integer :: ximg, ximages, yimg, yimages

    ! store the start (s) and end (e) for the i,j,k dimensions
    integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

    type(configuration_t) :: domain_configuration
  contains
    procedure :: default_initialize
    procedure :: get_grid_dimensions
    procedure :: initialize_from_file
    procedure :: advect
    procedure :: halo_send
    procedure :: halo_retrieve
    procedure :: halo_exchange
    procedure :: enforce_limits
    procedure :: domain_decomposition
    procedure :: report_convection
   !generic :: read(formatted)=>initialize_from_file
   !procedure, private :: initialize_with_configuration
   !procedure :: update_boundary
  end type

  integer, parameter :: space_dimension=3

  interface


    ! Set default component values
    module subroutine default_initialize(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    ! MPDATA algorithm
    module subroutine advect(this, dt)
      implicit none
      class(domain_t), intent(inout) :: this
      real,            intent(in)    :: dt
    end subroutine


    module subroutine halo_send(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine halo_retrieve(this, convected_particles)
      implicit none
      class(domain_t), intent(inout) :: this
      logical, intent(in) :: convected_particles
    end subroutine


    ! Exchange subdomain boundary information
    module subroutine halo_exchange(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    ! Make sure no hydrometeors are getting below 0
    module subroutine enforce_limits(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine domain_decomposition(this, nx, ny, nimages, ratio)
      implicit none
      class(domain_t), intent(inout) :: this
      integer,         intent(in)    :: nx, ny, nimages
      real,            intent(in), optional :: ratio
    end subroutine

    ! Return x, y, z dimensions of grid
    module function get_grid_dimensions(this, nx_extra, ny_extra) result(grid)
      implicit none
      class(domain_t), intent(in) :: this
      integer,         intent(in), optional :: nx_extra, ny_extra
      type(grid_t) :: grid
    end function

    ! Input domain_t object from file
    module subroutine initialize_from_file(this,file_name,convected_particles)
      implicit none
      class(domain_t), intent(inout) :: this
      character(len=*), intent(in) :: file_name
      logical, intent(in) :: convected_particles
    end subroutine

    module subroutine report_convection(this, step)
      implicit none
      class(domain_t), intent(inout) :: this
      integer, intent(in) :: step
    end subroutine

    elemental module function pressure_at_elevation(sealevel_pressure, elevation) result(pressure)
      implicit none
      real, intent(in) :: sealevel_pressure, elevation
      real :: pressure
    end function pressure_at_elevation

    elemental module function exner_function(pressure) result(exner)
      implicit none
      real, intent(in) :: pressure
      real :: exner
    end function exner_function

    elemental module function sat_mr(temperature,pressure)
      ! Calculate the saturated mixing ratio at a temperature (K), pressure (Pa)
      implicit none
      real,intent(in) :: temperature,pressure
      real :: sat_mr
    end function sat_mr

    subroutine load_sounding(sounding, sounding_n, sounding_start, sounding_max)
      real, allocatable, intent(out) :: sounding(:)
      integer, intent(out) :: sounding_n, sounding_start, sounding_max
    end subroutine load_sounding
  end interface


end module
