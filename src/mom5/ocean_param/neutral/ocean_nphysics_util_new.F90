module ocean_nphysics_util_new_mod
  !
  !<CONTACT EMAIL="tim.leslie@gmail.com"> Tim Leslie
  !</CONTACT>
  !
  !<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
  !</REVIEWER>
  !
  !<OVERVIEW>
  ! Utilities for neutral physics modules. 
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! Utilities for neutral physics modules. 
  !</DESCRIPTION>
  !
  ! <INFO>
  !
  ! <REFERENCE>
  ! D.B. Chelton,  R.A. deSzoeke, M.G. Schlax, K.E. Naggar, N. Siwertz
  ! Geographical Variability of the First Baroclinic Rossby Radius of Deformation
  ! Journal of Physical Oceanography (1998) vol 28 pages 433-460 
  ! </REFERENCE>
  !
  ! <REFERENCE>
  ! K. Eden and R. Greatbatch, 2008: Towards a mesoscale eddy closure,
  ! Ocean Modelling, vol. 20, pages 223-239
  ! </REFERENCE>
  !
  ! <REFERENCE>
  ! S.M. Griffies "Elements of MOM (2012) (EoM)"
  ! </REFERENCE>
  !
  ! <REFERENCE>
  ! S.M. Griffies 
  ! Fundamentals of Ocean Climate Models (2004)
  ! Princeton University Press 
  ! </REFERENCE>
  !  
  ! <NOTE>
  ! </NOTE>
  !
  ! </INFO>
  !
  !<NAMELIST NAME="ocean_nphysics_util_new_nml">
  !
  !  <DATA NAME="num_121_passes" TYPE="integer">
  !  For number of 1-2-1 passes through to smooth drhodz or 
  !  eady_rate in vertical. Default num_121_passes=1. 
  !  </DATA>
  !
  !</NAMELIST>

#define COMP isc:iec,jsc:jec
#define COMPXLYL isc-1:iec,jsc-1:jec
#define COMPXRYL isc:iec+1,jsc-1:jec
#define COMPXLYR isc-1:iec,jsc:jec+1

  use constants_mod,    only: epsln
  use diag_manager_mod, only: register_diag_field
  use fms_mod,          only: FATAL, open_namelist_file, check_nml_error, close_file
  use mpp_mod,          only: mpp_error, stdout, stdlog
  use mpp_mod,          only: input_nml_file
  use mpp_domains_mod,  only: mpp_update_domains
  
  use ocean_domains_mod,    only: get_local_indices
  use ocean_parameters_mod, only: omega_earth
  use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_time_type

  implicit none

  public ocean_nphysics_util_new_init
  public stencil_centre_to_vert
  public stencil_centre_to_horiz
  public vert_smooth
  public horz_smooth 
  public compute_rossby_radius

  private check_init
  private

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics_util_new.F90,v 20.0 2013/12/14 00:14:54 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


  logical :: module_is_initialized = .false.

  ! absolute value of the Coriolis parameter (sec^-1)
  real, dimension(:,:), allocatable :: coriolis_param

  ! beta = d(Coriolis)/dy (m^-1 sec^-1)
  real, dimension(:,:), allocatable :: beta_param

  !**************nml settings**************
  ! vert_smooth()
  integer :: num_121_passes = 1
  !**************end of nml settings**************

  namelist /ocean_nphysics_util_new_nml/ num_121_passes

contains

  !#######################################################################
  ! <SUBROUTINE NAME="stencil_centre_to_vert">
  !
  ! <DESCRIPTION>
  ! Initialise the grid indices and constants.
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_util_new_init(Domain, Grid)

    type(ocean_domain_type), intent(in) :: Domain
    type(ocean_grid_type),   intent(in) :: Grid

    integer :: ioun, ierr, io_status
    integer :: stdoutunit, stdlogunit
    stdoutunit = stdout()
    stdlogunit = stdlog()

    if (module_is_initialized) then
       call mpp_error(FATAL, &
            '==> Error from ocean_nphysics_util_new (ocean_nphysics_util_new_init): already initialized.')
    endif
    module_is_initialized = .true.

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_nphysics_util_new_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_util_new_nml')
#else
    ioun = open_namelist_file()
    read (ioun, ocean_nphysics_util_new_nml, IOSTAT=io_status)
    ierr = check_nml_error(io_status, 'ocean_nphysics_util_new_nml')
    call close_file(ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_nphysics_util_new_nml)
    write (stdlogunit, ocean_nphysics_util_new_nml)

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    ! Coriolis parameter and beta parameter
    allocate (coriolis_param(isd:ied,jsd:jed))
    allocate (beta_param(isd:ied,jsd:jed))
    coriolis_param(COMP) = 2.0*omega_earth*abs(sin(Grid%phit(COMP)))
    beta_param(COMP) = max(2.28e-11*abs(cos(Grid%phit(COMP))), epsln)

  end subroutine ocean_nphysics_util_new_init
  ! </SUBROUTINE> NAME="ocean_nphysics_util_new_init"


  !#######################################################################
  ! <SUBROUTINE NAME="stencil_centre_to_vert">
  !
  ! <DESCRIPTION>
  ! Take an array defined over the 4 quadrants centred on a tracer cell
  ! and shift them vertically to the four quadrants spanning the bottom
  ! cell face. The top and bottom half cells are set to zero.
  ! _____________       _____________
  ! |     :     |       |     :     |
  ! | 0,0 : 1,0 |       |     :     |
  ! |     :     !       |     :     |
  ! |-----T-----|  ==>  |-----T-----|
  ! |     :     !       |     :     |
  ! | 0,1 : 1,1 !       | 0,0 : 1,0 |
  ! |_____:_____|       |_____:_____|
  !                     |     :     |
  !                     | 0,1 : 1,1 |
  !                     |     :     |
  !                     |-----o-----|
  !                     |     :     |
  !                     |     :     |
  !                     |_____:_____|
  ! </DESCRIPTION>
  !
  subroutine stencil_centre_to_vert(array_xz, array_yz, array_v_xz, array_v_yz)

    real, dimension(isd:,jsd:,:,0:,0:), intent(in)    :: array_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in)    :: array_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: array_v_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: array_v_yz

    integer :: k

    call check_init("stencil_centre_to_vert")

    do k = 1, nk-1
       array_v_xz(COMP,k,:,0) = array_xz(COMP,k,  :,1)
       array_v_xz(COMP,k,:,1) = array_xz(COMP,k+1,:,0)
       array_v_yz(COMP,k,:,0) = array_yz(COMP,k,  :,1)
       array_v_yz(COMP,k,:,1) = array_yz(COMP,k+1,:,0)
    enddo
    array_v_xz(COMP,nk,:,0) = array_xz(COMP,nk,:,1)
    array_v_yz(COMP,nk,:,0) = array_yz(COMP,nk,:,1)
    array_v_xz(COMP,nk,:,1) = 0.0
    array_v_yz(COMP,nk,:,1) = 0.0
    
  end subroutine stencil_centre_to_vert
  ! </SUBROUTINE> NAME="stencil_centre_to_vert"


  !#######################################################################
  ! <SUBROUTINE NAME="stencil_centre_to_horiz">
  !
  ! <DESCRIPTION>
  ! Take an array defined over the 4 quadrants centred on a tracer cell
  ! and shift them horizontally to the four quadrants spanning the right hand
  ! cell face.
  ! _____________       _________________________
  ! |     :     |       |     :     |     :     |
  ! | 0,0 : 1,0 |       |     : 0,0 | 1,0 :     |
  ! |     :     !       |     :     |     :     |
  ! |-----T-----|  ==>  |-----T-----|-----o-----|
  ! |     :     !       |     :     |     :     |
  ! | 0,1 : 1,1 !       |     : 0,1 | 1,1 :     |
  ! |_____:_____|       |_____:_____|_____:_____|
  ! </DESCRIPTION>
  !
  subroutine stencil_centre_to_horiz(array_xz, array_yz, array_h_xz, array_h_yz)

    real, dimension(isd:,jsd:,:,0:,0:), intent(in)    :: array_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in)    :: array_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: array_h_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: array_h_yz

    call check_init("stencil_centre_to_horz")

    array_h_xz(COMPXLYL,:,0,:) = array_xz(COMPXLYL,:,1,:)
    array_h_xz(COMPXLYL,:,1,:) = array_xz(COMPXRYL,:,0,:)

    array_h_yz(COMPXLYL,:,0,:) = array_yz(COMPXLYL,:,1,:)
    array_h_yz(COMPXLYL,:,1,:) = array_yz(COMPXLYR,:,0,:)

  end subroutine stencil_centre_to_horiz
  ! </SUBROUTINE> NAME="stencil_centre_to_horiz"


  !#######################################################################
  ! <SUBROUTINE NAME="vert_smooth">
  !
  ! <DESCRIPTION>
  ! Apply a vertical smoothing using a weighted 3 point stencil.
  !
  !                    1
  ! Stencil weights :  2
  !                    1
  ! 
  ! Smoothing is applied to all points in the column which admit
  ! the full stencil (so vertical boundaries are not modified).
  ! 
  ! The smoothing is applied multiple times, controlled by the namelist
  ! parameter num_121_passes.
  ! </DESCRIPTION>
  !
  subroutine vert_smooth(kmt, array)
    integer, dimension(isd:,jsd:), intent(in)    :: kmt
    real, dimension(isd:,jsd:,:),  intent(inout) :: array

    integer :: m, i, j, kbot

    call check_init("vert_smooth")

    do m=1,num_121_passes
       do j=jsd,jed
          do i=isd,ied
             kbot = kmt(i,j)
             if (kbot > 1) then
                array(i,j,2:kbot-1) = (array(i,j,1:kbot-2) + 2*array(i,j,2:kbot-1) + array(i,j,3:kbot))/4.0
             endif
          enddo
       enddo
    enddo

  end subroutine vert_smooth
  ! </SUBROUTINE> NAME="vert_smooth"


  !#######################################################################
  ! <SUBROUTINE NAME="horz_smooth">
  !
  ! <DESCRIPTION>
  ! Apply a horizontal smoothing using a weighted five point stencil.
  !
  !                     1
  ! Stencil weights : 1 4 1
  !                     1
  ! 
  ! Only active cells, as defined by the tmask, are used in the calculation.
  ! Operates over the COMP domain.
  ! </DESCRIPTION>
  !
  subroutine horz_smooth(Domain, tmask, array)

    type(ocean_domain_type),      intent(inout) :: Domain
    real, dimension(isd:,jsd:,:), intent(in)    :: tmask
    real, dimension(isd:,jsd:,:), intent(inout) :: array

    real, dimension(isd:ied,jsd:jed) :: tmp
    integer :: i, j, k
    real :: active_cells

    call check_init("horz_smooth")

    array(COMP,:) = tmask(COMP,:)*array(COMP,:)
    call mpp_update_domains(array(:,:,:), Domain%domain2d)
    do k = 1, nk
       tmp(:,:) = 0.0
       do j=jsc,jec
          do i=isc,iec
             if (tmask(i,j,k) == 1.0) then
                active_cells = 4.0   +&
                     tmask(i-1,j,k)  +&
                     tmask(i+1,j,k)  +&
                     tmask(i,j-1,k)  +&
                     tmask(i,j+1,k)
                if (active_cells > 4.0) then
                   tmp(i,j) = &
                        (4.0*array(i,j,k) +&
                        array(i-1,j,k)    +&
                        array(i+1,j,k)    +&
                        array(i,j-1,k)    +&
                        array(i,j+1,k)) / active_cells
                else
                   tmp(i,j) = array(i,j,k)
                endif
             endif
          enddo
       enddo
       array(COMP,k) = tmp(COMP)
    enddo

  end subroutine horz_smooth
  ! </SUBROUTINE> NAME="horz_smooth"



  !#######################################################################
  ! <SUBROUTINE NAME="compute_rossby_radius">
  !
  ! <DESCRIPTION>
  ! Subroutine computes the first baroclinic Rossby radius of deformation.
  ! Employ WKB approach described by Chelton et al.  In particular,
  ! use formulae (2.2), (2.3a) and (2.3b) from their paper.  
  ! </DESCRIPTION>
  !
  subroutine compute_rossby_radius(gravity_wave_speed, rossby_radius)

    real, dimension(isd:,jsd:), intent(in)    :: gravity_wave_speed
    real, dimension(isd:,jsd:), intent(inout) :: rossby_radius

    real, dimension(isd:ied,jsd:jed) :: rossby_non_equator
    real, dimension(isd:ied,jsd:jed) :: rossby_equator

    ! EoM: 16.62, G(2004) 14.81
    rossby_non_equator(COMP) = gravity_wave_speed(COMP)/(coriolis_param(COMP) + epsln)
    
    ! EoM: 16.63, G(2004) 14.82    
    rossby_equator(COMP)     = sqrt(gravity_wave_speed(COMP)/(2.0*beta_param(COMP)))

    ! Eden and Greatbatch section 4.2    
    rossby_radius(COMP)      = min(rossby_non_equator(COMP), rossby_equator(COMP))          

  end subroutine compute_rossby_radius
  ! </SUBROUTINE> NAME="compute_rossby_radius"

  !#######################################################################
  ! <SUBROUTINE NAME="check_init">
  !
  ! <DESCRIPTION>
  ! Helper function to ensure the module is initialised when calling the
  ! module's functions.
  ! </DESCRIPTION> 
  !
  subroutine check_init(function_name)

    character(len=*), intent(in) :: function_name

    if (.not. module_is_initialized) then 
       call mpp_error(FATAL, &
            '==> Error '//FILENAME//' (' // function_name // '): needs initialization.')
    endif

  end subroutine check_init
  ! </SUBROUTINE> NAME="check_init"


end module ocean_nphysics_util_new_mod
