module ocean_lapgen_friction_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes the thickness weighted time tendency for  
! horizontal velocity arising from horizontal Laplacian friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted time tendency for
! horizontal velocity arising from horizontal Laplacian friction. 
! The viscosity used to determine the strength of the tendency 
! can be a general function of space and time as specified by 
! the Smagorinsky approach as well as a grid-scale dependent
! background viscosity.  The form of the friction operator 
! can be isotropic or anisotropic in the horizontal plane. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies and R.W. Hallberg, 2000: 
! Biharmonic friction with a Smagorinsky viscosity for use in large-scale
! eddy-permitting ocean models
! Monthly Weather Review, vol. 128, pages 2935-2946
! </REFERENCE>
!
! <REFERENCE>
! R. D. Smith and J. C. McWilliams, 2003:
! Anisotropic horizontal viscosity for ocean models,
! Ocean Modelling, vol. 5, pages 129-156.
! </REFERENCE>
!
! <REFERENCE>
! Maltrud and Holloway, 2008: Implementing biharmonic neptune in a
! global eddying ocean model, Ocean Modelling, vol. 21, pages 22-34.
! </REFERENCE>
!
! <REFERENCE>
! Deremble, Hogg, Berloff, and Dewar, 2011:
! On the application of no-slip lateral boundary conditions to coarsely
! resolved ocean models, Ocean Modelling. 
! </REFERENCE>
!
! <REFERENCE>
! Griffies: Elements of MOM (2012)
! </REFERENCE>
!
! <NOTE>
! The ocean model can generally run with both Laplacian and biharmonic friction
! enabled at the same time.  Such has been found useful for some eddying 
! ocean simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_lapgen_friction_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA>
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging by printing checksums.  
!  </DATA> 
!
!  <DATA NAME="viscosity_scale_by_rossby" TYPE="logical">
!  To scale down the laplacian viscosity according to the relative scale of the 
!  horizontal grid and the first baroclinic Rossby radius. This is a useful 
!  scheme for models that resolve the Rossby radius in the lower latitudes, and so
!  presumably do not wish to have much laplacian friction, whereas the higher latitudes
!  need more friction.  Default viscosity_scale_by_rossby=.false.
!  </DATA> 
!
!  <DATA NAME="viscosity_scale_by_rossby_power" TYPE="real">
!  The power used to determine the viscosity scaling function. 
!  Default viscosity_scale_by_rossby_power=2.0.
!  </DATA> 
!
!  <DATA NAME="divergence_damp" TYPE="logical">
!  To damp the divergence field.  
!  </DATA> 
!  <DATA NAME="divergence_damp_vel_micom" TYPE="real" UNITS="m/s">
!  Velocity scale to set the viscosity used with divergence damping. 
!  </DATA> 
!
!  <DATA NAME="k_smag_iso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="k_smag_aniso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky anisotropic viscosity. 
!  </DATA> 
!  <DATA NAME="viscosity_ncar" TYPE="logical">
!  Anisotropic background viscosities used by NCAR. 
!  </DATA> 
!  <DATA NAME="viscosity_ncar_2000" TYPE="logical">
!  Anisotropic background viscosities used by NCAR, using the 
!  formulation as of 2000.  Default viscosity_ncar_2000=.true.
!  </DATA> 
!  <DATA NAME="viscosity_ncar_2007" TYPE="logical">
!  Anisotropic background viscosities used by NCAR, using the 
!  formulation as of 2007.  Default viscosity_ncar_2007=.false.
!  </DATA> 
!  <DATA NAME="vel_micom_iso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_aniso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity. 
!  </DATA> 
!
!  <DATA NAME="equatorial_zonal" TYPE="logical">
!  Orient the anisotropic friction within a latitudinal band according to zonal direction. 
!  </DATA> 
!  <DATA NAME="equatorial_zonal_lat" TYPE="real">
!  Latitudinal band to use the zonal friction orientation. 
!  </DATA> 
!  <DATA NAME="ncar_isotropic_off_equator" TYPE="logical">
!  Polewards of equatorial_zonal_lat, revert NCAR scheme to isotropic 
!  </DATA> 
!  <DATA NAME="equatorial_no_smag" TYPE="logical">
!  Turn smag off within equatorial_zonal_lat region. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_iso" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity within
!  a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_aniso" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity within
!  a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_lat_micom" TYPE="real">
!  Equatorial latitude band (degrees) within which the MICOM viscosity is set according 
!  to eq_vel_micom_iso and eq_vel_micom_aniso.
!  </DATA> 
!
!  <DATA NAME="restrict_polar_visc" TYPE="logical">
!  For restricting the background viscosity poleward of a 
!  latitude.  This method may be useful for coupling to an ice model
!  in which case the horizontal viscosity may need to be a bit 
!  smaller to maintain time step constraints.  This is because the 
!  effective friction is larger than that just within the ocean.  
!  </DATA> 
!  <DATA NAME="restrict_polar_visc_lat" TYPE="real">
!  Latitude poleward of which we restrict the viscosity.
!  </DATA> 
!  <DATA NAME="restrict_polar_visc_ratio" TYPE="real">
!  Ratio of the normal critical value that we limit the 
!  viscosity to be no greater than.  If restrict_polar_visc_ratio=1.0
!  then there is no special limitation of the viscosity beyond that 
!  of the one-dimensional stability constraint.  
!  </DATA> 
!
!  <DATA NAME="bottom_5point" TYPE="logical">
!  To alleviate problems with small partial cells, it is often necessary to reduce the 
!  operator to the traditional 5-point Laplacian at the ocean bottom.  This logical 
!  implements this mixing. Default bottom_5point=.false.
!  </DATA> 
!
!  <DATA NAME="neptune" TYPE="logical">
!  Set to true for computing friction relative to Neptune barotropic velocity. 
!  Default neptune=.false. 
!  </DATA> 
!  <DATA NAME="neptune_length_eq" UNITS="m" TYPE="real">
!  Length scale used to compute Neptune velocity at equator.  
!  </DATA> 
!  <DATA NAME="neptune_length_pole" UNITS="m" TYPE="real">
!  Length scale used to compute Neptune velocity at pole. 
!  </DATA> 
!  <DATA NAME="neptune_depth_min" UNITS="m" TYPE="real">
!  Minimum depth scale used for computing Neptune velocity.
!  Default neptune_depth_min=100.0
!  </DATA> 
!  <DATA NAME="neptune_smooth" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the diagnosed  
!  neptune velocity scale. 
!  Default neptune_smooth=.true.
!  </DATA> 
!  <DATA NAME="neptune_smooth_num" TYPE="integer">
!  Number of smoothing passes for neptune velocity.
!  Default neptune_smooth_num=1.
!  </DATA> 
!
!  <DATA NAME="vconst_1" UNITS="cm^2/sec" TYPE="real">
!  Background viscosity for NCAR algorithm.
!  </DATA> 
!  <DATA NAME="vconst_2" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_3" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_4" UNITS="1/cm" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_5" TYPE="integer">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_6" UNITS="cm^2/sec" TYPE="real">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_7" UNITS="cm/sec">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="vconst_8" UNITS="degrees">
!  For NCAR viscosity algorithm.
!  </DATA> 
!  <DATA NAME="visc_vel_scale_length" UNITS="cm">
!  For NCAR viscosity algorithm: efolding depth for 
!  depth dependent background viscosity.  
!  Default visc_vel_scale_length=1500.e2 cm
!  </DATA> 
!  <DATA NAME="ncar_isotropic_at_depth" TYPE="logical">
!  Sets the NCAR scheme to be isotropic beneath a chosen depth.
!  </DATA> 
!  <DATA NAME="ncar_isotropic_depth" TYPE="real" UNITS="m">
!  Sets the NCAR scheme to be isotropic beneath this chosen depth.
!  </DATA> 
!  <DATA NAME="ncar_isotropic_at_depth_visc" TYPE="real" UNITS="m2/sec">
!  Sets the NCAR scheme to be isotropic beneath this chosen depth, with 
!  minimum viscosity set according to this value. 
!  </DATA> 
!  <DATA NAME="debug_ncar_A" TYPE="logical">
!  Sets f_perp=f_para for debugging purposes with the NCAR scheme.
!  </DATA> 
!  <DATA NAME="debug_ncar_B" TYPE="logical">
!  Sets f_para=f_perp for debugging purposes with the NCAR scheme.
!  </DATA> 
!
!  <DATA NAME="use_side_drag_friction" TYPE="logical">
!  For converting friction at U-cells next to walls into 
!  a drag law, as per Deremble et al. Use cdbot_array
!  from ocean_core/ocean_bbc.F90 to compute drag force. 
!  Default use_side_drag_friction=.false.
!  </DATA> 
!  <DATA NAME="side_drag_friction_scaling" TYPE="real">
!  Dimensionless scaling used for cdbot_array when setting
!  side drag friction. So the effective side dragy coefficient
!  is side_drag_friction_scaling*cdbot_array.  
!  Default side_drag_friction_scaling=1.0.
!  </DATA> 
!  <DATA NAME="side_drag_friction_uvmag_max" UNITS="m/s" TYPE="real">
!  Maximum magnitude of horizontal velocity used to compute the 
!  side drag friction. This parameter can be useful especially
!  for pressure models where the bottom cells can be quite thin 
!  and subject to sporadic large magnitudes.  We do the same thing with 
!  bottom drag calculations. 
!  Default side_drag_friction_uvmag_max=10.0.
!  </DATA> 
!  <DATA NAME="side_drag_friction_max" UNITS="N/m^2" TYPE="real">
!  Maximum magnitude of the side drag induced friction. 
!  This parameter can be useful especially for pressure models 
!  where the bottom cells can be quite thin and subject to sporadic
!  large magnitudes.  We do the same thing with bottom drag calculations. 
!  Default side_drag_friction_max=1.0.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: pi, radius, epsln, radian
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: open_namelist_file, check_nml_error, write_version_number, close_file
use mpp_domains_mod,  only: mpp_update_domains, BGRID_NE, mpp_global_field, XUPDATE 
use mpp_domains_mod,  only: mpp_start_update_domains, mpp_complete_update_domains
use mpp_mod,          only: input_nml_file, mpp_sum, mpp_pe, mpp_error, mpp_max
use mpp_mod,          only: FATAL, NOTE, stdout, stdlog

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: BAY, BAX, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_operators_mod,  only: FAY, FAX, FDX_NT, FDY_ET, FDX_T, FDY_T 
use ocean_parameters_mod, only: missing_value, omega_earth, rho0
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type
use ocean_types_mod,      only: ocean_domain_type, ocean_adv_vel_type 
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,      only: ocean_options_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d_u, diagnose_3d_u, diagnose_2d, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_2d, wrk1_v, wrk2_v, wrk3_v, wrk1_v2d  

implicit none

private

public ocean_lapgen_friction_init
public lapgen_friction
public lapgen_viscosity_check
public lapgen_reynolds_check

private BDX_EU_smag
private BDY_NU_smag
private compute_neptune_velocity

! length of forward time step  
real :: dtime

! for diagnostics 
integer :: id_aiso                 =-1
integer :: id_aaniso               =-1
integer :: id_along_back           =-1
integer :: id_across_back          =-1
integer :: id_aiso_diverge         =-1
integer :: id_neptune_lap_u        =-1
integer :: id_neptune_lap_v        =-1
integer :: id_neptune_psi          =-1
integer :: id_neptune_ft           =-1
integer :: id_along                =-1
integer :: id_across               =-1
integer :: id_visc_crit_lap        =-1
integer :: id_lap_fric_u           =-1
integer :: id_lap_fric_v           =-1
integer :: id_horz_lap_diss        =-1
integer :: id_viscosity_scaling    =-1
integer :: id_umask_next_to_land   =-1
integer :: id_lap_plus_side_fric_u =-1
integer :: id_lap_plus_side_fric_v =-1
integer :: id_side_drag_friction_u =-1
integer :: id_side_drag_friction_v =-1

logical :: used

real :: k_smag_iso            = 2.0     ! smag scaling coeff for isotropic viscosity (dimensionless) 
real :: k_smag_aniso          = 0.0     ! smag scaling coeff for anisotripic viscosity (dimensionless) 
real :: vel_micom_iso         = 0.0     ! background scaling velocity for isotropic viscosity (m/sec) 
real :: vel_micom_aniso       = 0.0     ! background scaling velocity for anisotropic viscosity (m/sec) 
real :: eq_vel_micom_iso      = 0.0     ! background scaling velocity (m/sec) within equatorial band for isotropic visc
real :: eq_vel_micom_aniso    = 0.0     ! background scaling velocity (m/sec) within equatorial band for anisotropic visc
real :: eq_lat_micom          = 0.0     ! equatorial latitude band for micom (degrees)
real :: equatorial_zonal_lat  = 0.0     ! latitudinal band for orienting the friction along zonal direction
logical :: equatorial_zonal   = .false. ! zonally orient anisotropic friction w/i equatorial_zonal_lat
logical :: equatorial_no_smag = .false. ! remove smag within equatorial_zonal_lat
logical :: bottom_5point      = .false. ! for bottom Laplacian 5point mixing to avoid problems with thin partial cells. 

! for scaling the viscosity down according to the Rossby radius 
logical :: viscosity_scale_by_rossby=.false. 
real    :: viscosity_scale_by_rossby_power=2.0

! for restricting critical viscosity in high latitudes.
! this has been found of use for coupling to ice models, where 
! the effective friction is larger than just that within the ocean.   
logical :: restrict_polar_visc=.false.  
real    :: restrict_polar_visc_lat=60.0
real    :: restrict_polar_visc_ratio=0.35

! for divergence damping 
logical  :: divergence_damp           = .false.
real     :: divergence_damp_vel_micom = 0.0
                                 
! From NCAR viscosity algorithm.  Defaults are from CCSM2.0 ocean component gx1v3 
integer :: nx, ny                         ! number of global points in generalized x and y directions
real :: vconst_1              = 1.e7      ! (cm^2/s)
real :: vconst_2              = 0.0       ! (dimensionless)
real :: vconst_3              = 0.16      ! (dimensionless)
real :: vconst_4              = 2.e-8     ! (1/cm) 
integer :: vconst_5           = 3         ! (number of grid points)
real :: vconst_6              = 1.e7      ! (cm^2/sec)
real :: vconst_7              = 100.0     ! (cm/sec)
real :: vconst_8              = 45.0      ! (degrees on a circle)
real :: visc_vel_scale_length = 1500.e2   ! cm efolding depth scale 
logical :: viscosity_ncar     = .false.   ! to get the (x,y,z) dependent isotropic 
                                          ! and anisotropic viscosities used by NCAR.  
logical :: viscosity_ncar_2000= .true.    ! to get the (x,y,z) dependent isotropic 
                                          ! and anisotropic viscosities used by NCAR as of 2000.  
logical :: viscosity_ncar_2007= .false.   ! to get the (x,y,z) dependent isotropic 
                                          ! and anisotropic viscosities used by NCAR as of 2007.  
logical :: debug_ncar_A     = .false.
logical :: debug_ncar_B     = .false.
logical :: ncar_isotropic_off_equator = .false.
logical :: ncar_only_equatorial = .false.  ! for use of the ncar scheme only within the tropical band 

logical :: ncar_isotropic_at_depth      = .false. ! to set the NCAR scheme to be isotropic beneath a depth
real    :: ncar_isotropic_depth         = 4000.0  ! depth (m) beneath which set NCAR scheme to isotropic 
real    :: ncar_isotropic_at_depth_visc = 1e4     ! m2/sec minimum viscosity for ncar_isotropic_at_depth

! for neptune scheme  
logical :: neptune              = .false. ! for computing friction relative to barotropic Neptune velocity
logical :: neptune_smooth       = .true.  ! for smoothing diagnosed neptune velocity
integer :: neptune_smooth_num   = 1       ! number of smoothing passes for neptune smoother 
real    :: neptune_length_eq    = 1.2e3   ! (metres)
real    :: neptune_length_pole  = 3.0e3   ! (metres)
real    :: neptune_depth_min    = 100.0   ! (metres)

! for side-boundary drag law
logical :: use_side_drag_friction       = .false.
real    :: side_drag_friction_scaling   = 1.0
real    :: side_drag_friction_max       = 1.0
real    :: side_drag_friction_uvmag_max = 10.0
real, dimension(:,:,:), allocatable :: umask_next_to_land


#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed)    :: aiso_diverge ! (m^2/s) viscosity for divergence damping
 
real, dimension(isd:ied,jsd:jed)    :: fsmag_iso   ! (m^2) combination of terms for computing isotropic smag visc
real, dimension(isd:ied,jsd:jed)    :: fsmag_aniso ! (m^2) combination of terms for computing anisotropic smag visc
real, dimension(isd:ied,jsd:jed,nk) :: aiso_back   ! minimum allowable isotropic viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk) :: aaniso_back ! minimum allowable anisotropic viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk) :: aiso        ! isotropic viscosity (m^2/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed)    :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed)    :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed,nk) :: massqc      ! (1/4)*dxu*dyu*rho_dzu (kg) for time-independent quarter-cell mass
real, dimension(isd:ied,jsd:jed)    :: visc_crit   ! critical value of the viscosity (m^2/sec) for linear stability 

real, dimension(isd:ied,jsd:jed,0:1,0:1,2) :: stress   !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(isd:ied,jsd:jed,2)         :: tmpfdelx
real, dimension(isd:ied,jsd:jed,2)         :: tmpfdely
real, dimension(isd:ied,jsd:jed)           :: tmp
real, dimension(isd:ied,jsd:jed)           :: cos2theta
real, dimension(isd:ied,jsd:jed)           :: sin2theta

#else

real, dimension(:,:), allocatable   :: aiso_diverge ! (m^2/s) viscosity for divergence damping
real, dimension(:,:), allocatable   :: fsmag_iso    ! (m^2) combination of terms for computing isotropic smag visc
real, dimension(:,:), allocatable   :: fsmag_aniso  ! (m^2) combination of terms for computing anisotropic smag visc
real, dimension(:,:,:), allocatable :: aiso_back    ! minimum allowable isotropic viscosity (m^2/s)
real, dimension(:,:,:), allocatable :: aaniso_back  ! minimum allowable anisotropic viscosity (m^2/s)
real, dimension(:,:,:), allocatable :: aiso         ! isotropic viscosity (m^2/s) averaged to U cell for diagnostics
real, dimension(:,:), allocatable   :: daur_dxur    ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(:,:), allocatable   :: daur_dyur    ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(:,:,:), allocatable :: massqc       ! (1/4)*dxu*dyu*rho_dzu (kg) for time-independent quarter-cell mass 
real, dimension(:,:), allocatable   :: visc_crit    ! critical value of the viscosity (m^2/sec) for linear stability 

real, dimension(:,:,:,:,:), allocatable        :: stress   !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(:,:,:), allocatable            :: tmpfdelx
real, dimension(:,:,:), allocatable            :: tmpfdely
real, dimension(:,:), allocatable              :: tmp
real, dimension(:,:), allocatable              :: cos2theta
real, dimension(:,:), allocatable              :: sin2theta


#endif


real, dimension(:,:,:), allocatable :: neptune_velocity  ! barotropic velocity (m/s) from Neptune 
real, dimension(:,:),   allocatable :: grid_length       ! horizontal grid length (m)
real, dimension(:,:),   allocatable :: viscosity_scaling ! dimensionless function of grid scale and Rossby radius

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

character(len=256) :: version=&
     '=>Using: ocean_lapgen_friction.F90 ($Id: ocean_lapgen_friction.F90,v 20.0 2013/12/14 00:14:26 fms Exp $)'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc              = .FALSE.
logical :: async_domain_update   = .false.
integer :: blocksize             = 10

namelist /ocean_lapgen_friction_nml/ use_this_module, debug_this_module,                &
    bottom_5point, k_smag_iso, k_smag_aniso,                                            &
    vel_micom_iso, vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom, &
    equatorial_zonal, equatorial_zonal_lat, equatorial_no_smag,                         &
    viscosity_ncar, viscosity_ncar_2000, viscosity_ncar_2007,                           &
    ncar_isotropic_off_equator, ncar_only_equatorial,                                   &
    vconst_1, vconst_2, vconst_3, vconst_4, vconst_5, vconst_6, vconst_7, vconst_8,     &
    debug_ncar_A, debug_ncar_B, visc_vel_scale_length,                                  &
    neptune, neptune_length_eq, neptune_length_pole, neptune_depth_min,                 &
    neptune_smooth, neptune_smooth_num,                                                 &
    restrict_polar_visc, restrict_polar_visc_lat, restrict_polar_visc_ratio,            &
    ncar_isotropic_at_depth, ncar_isotropic_depth, ncar_isotropic_at_depth_visc,        &
    divergence_damp, divergence_damp_vel_micom,                                         &
    viscosity_scale_by_rossby, viscosity_scale_by_rossby_power, async_domain_update,    &
    blocksize, use_side_drag_friction, side_drag_friction_scaling,                      &
    side_drag_friction_uvmag_max, side_drag_friction_max 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_lapgen_friction_init">

! <DESCRIPTION>
! Initialize the horizontal Laplacian friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_lapgen_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                      obc, use_lapgen_friction, debug)

  type(ocean_grid_type),    target, intent(in)    :: Grid
  type(ocean_domain_type),  target, intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: d_time
  logical,                          intent(in)    :: obc
  logical,                          intent(inout) :: use_lapgen_friction
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real :: coeff_iso, coeff_aniso, dxdy

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapgen_friction_mod (ocean_lapgen_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif 

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_lapgen_friction_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_lapgen_friction_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_lapgen_friction_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_lapgen_friction_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_lapgen_friction_nml)  
  write (stdlogunit,ocean_lapgen_friction_nml)

  if(use_this_module) then 
      call mpp_error(NOTE, '==> NOTE: USING ocean_lapgen_friction_mod.')
      Ocean_options%horz_lap_friction = 'Used general horizontal Laplacian friction.'
      use_lapgen_friction = .true. 
  else 
      call mpp_error(NOTE, '==> NOTE: NOT using ocean_lapgen_friction_mod.')
      Ocean_options%horz_lap_friction = 'Did NOT use horizontal Laplacian friction.'
      use_lapgen_friction = .false. 
      return 
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  if(async_domain_update) then 
    write(stdoutunit,'(a)') &
    '==>Note: using asynchronous domain update in the vertical loop.'  
  else
    write(stdoutunit,'(a)') &
    '==>Note: not using asynchronous domain update in the vertical loop. This may be slow.'  
  endif 

  dtime = d_time
  write(stdoutunit,'(/a,f10.2)') &
  '==> Note from ocean_lapgen_friction_mod: using forward time step of (secs)', dtime 

  have_obc       = obc
  if (have_obc) write(stdoutunit,'(/a,f10.2)') &
   '==> Note from ocean_lapgen_friction_mod: considering obc' 

  if(viscosity_scale_by_rossby) then 
    write(stdoutunit,'(1x,a)') '==> Note: Scaling the laplacian viscosity according to grid scale and Rossby radius.'
  endif 

  if(bottom_5point) then 
    write(stdoutunit,'(1x,a)') '==> Note: Will reduce horizontal friction &
         &to a 5point Laplacian on the bottom'
    write(stdoutunit,'(5x,a)') 'This helps to alleviate numerical problems&
         & with thin bottom partial cells.'
  endif 

  if(neptune) then 
    write(stdoutunit,'(1x,a)') '==> Note: Computing Laplacian friction relative to Neptune'
    write(stdoutunit,'(5x,a,e10.5)') &
    'Equatorial length (m) for computing Neptune velocity is   ',neptune_length_eq 
    write(stdoutunit,'(5x,a,e10.5)') &
    'Polar length (m) for computing Neptune velocity is        ',neptune_length_pole 
    write(stdoutunit,'(5x,a,e10.5)') &
    'Minimum depth scale (m) for computing Neptune velocity is ',neptune_depth_min
  endif 

  if(viscosity_ncar)  then 

    write( stdoutunit,'(1x,a)')'==> NOTE: USING background horz viscosities &
         &according to NCAR CCSM2.0 algorithm.' 
    if(viscosity_ncar_2000) then 
       write( stdoutunit,'(1x,a)')'==> NOTE: USING NCAR viscosity as formulated in 2000.' 
    endif 
    if(viscosity_ncar_2007) then 
       write( stdoutunit,'(1x,a)')'==> NOTE: USING NCAR viscosity as formulated in 2007.' 
    endif 
    if(viscosity_ncar_2000 .and. viscosity_ncar_2007) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_lapgen_friction_mod: cannot choose both viscosity_ncar_2000 and viscosity_ncar_2007')
    endif 
    if(.not. viscosity_ncar_2000 .and. .not. viscosity_ncar_2007) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_lapgen_friction_mod: viscosity_ncar requires one of viscosity_ncar_2000 or viscosity_ncar_2007')
    endif 

    write( stdoutunit,'(a,e15.4)')'  NCAR vconst_1 (cm^2/sec)  = ',vconst_1
    write( stdoutunit,*)' NCAR vconst_2             = ',vconst_2
    write( stdoutunit,*)' NCAR vconst_3             = ',vconst_3
    write( stdoutunit,*)' NCAR vconst_4 (1/cm)      = ',vconst_4
    write( stdoutunit,*)' NCAR vconst_5             = ',vconst_5
    write( stdoutunit,'(a,e15.4)')'  NCAR vconst_6 (cm^2/sec)  = ',vconst_6
    write( stdoutunit,*)' NCAR vconst_7 (cm/sec)    = ',vconst_7
    write( stdoutunit,*)' NCAR vconst_8 (degrees)   = ',vconst_8
    if(debug_ncar_A) then 
      write( stdoutunit,'(1x,a)')'==> NOTE: USING debug_ncar_A=.true., &
           &so will set f_perp=f_para'
    endif 
    if(debug_ncar_B) then 
      write( stdoutunit,'(1x,a)')'==> NOTE: USING debug_ncar_B=.true.,&
           & so will set f_para=f_perp'
    endif 
    if(ncar_isotropic_off_equator) then 
      write( stdoutunit,'(1x,a,f10.4)')&
      '==>ncar_isotropic_off_equator=.true. =>f_para=f_perp poleward of lat',equatorial_zonal_lat
    endif 
    if(ncar_only_equatorial) then 
      write( stdoutunit,'(1x,a,f10.4)')&
      '==>ncar_only_equatorial=.true. =>ncar scheme only in band +/- lat',equatorial_zonal_lat
    endif 

  endif  

  if(k_smag_iso > 0.0) then 
     write( stdoutunit,'(1x,a)')'==> NOTE: USING horz isotropic viscosity &
          &via Smagorinsky.'
  endif 
  if(k_smag_aniso > 0.0) then 
     write( stdoutunit,'(1x,a)')'==> NOTE: USING horz anisotropic viscosity &
          &via Smagorinsky.'
  endif      
  if(k_smag_iso==0.0)    write( stdoutunit,'(1x,a)')'==> NOTE: Setting horz &
       &isotropic Smagorinsky viscosity to zero.'
  if(k_smag_aniso==0.0)  write( stdoutunit,'(1x,a)')'==> NOTE: Setting horz &
       &anisotropic Smagorinsky viscosity to zero.'

  if(.not. viscosity_ncar) then 
    write( stdoutunit,'(1x,a)') &
    '==> NOTE: NOT using NCAR scheme for computing viscosities. '
    if(vel_micom_iso > 0.0) then 
       write( stdoutunit,'(1x,a)') &
       '==> NOTE: USING background horz isotropic viscosity via MICOM.'
    endif
    if(vel_micom_aniso > 0.0) then 
       write( stdoutunit,'(1x,a)') &
       '==> NOTE: USING background horz anisotropic viscosity via MICOM.'
    endif 
    if(vel_micom_iso==0.0)  write( stdoutunit,'(1x,a)')  '==> NOTE: USING zero &
         &background horz isotropic viscosity.'
    if(vel_micom_aniso==0.0)  write( stdoutunit,'(1x,a)')'==> NOTE: USING zero &
         &background horz anisotropic viscosity.' 
  endif  
  if(eq_lat_micom > 0) write( stdoutunit,'(1x,a)')'==> NOTE: USING different &
       &background horz viscosity within equatorial zone.'

  if(k_smag_aniso > 2.0*k_smag_iso) then
    write( stdoutunit,'(1x,a)')'==> ERROR: Smag horz anisotropic visc too large. &
         &Must be less than twice Smagorinsky isotropic visc.'
  endif  
  if(vel_micom_aniso > 2.0*vel_micom_iso) then
    write( stdoutunit,'(1x,a)')'==> ERROR: Background anisotropic visc too large.&
         & Must be less than twice back isotropic visc.'
  endif  

  if(equatorial_zonal) then
    write( stdoutunit,'(/a)') &
    'If using anisotropic friction, zonally orient the friction within a latitudinal band'
    write( stdoutunit,'(a,f12.5,a)') &
    ' of width ',equatorial_zonal_lat,' degrees.'
  endif 

  if(restrict_polar_visc) then
    write( stdoutunit,'(/a,f12.5)')&
    ' Using restrict_polar_visc to lower visc_crit poleward of (deg)',restrict_polar_visc_lat
    write( stdoutunit,'(a,f12.5)')&
    ' by an amount given by the fraction ',restrict_polar_visc_ratio
    write( stdoutunit,'(a/)')&
    ' This approach is useful when coupling to ice, where effective (ocn+ice) visc > ocn visc.'
  endif 

  Grd => Grid
  Dom => Domain
  nx  =  Grd%ni
  ny  =  Grd%nj

#ifndef MOM_STATIC_ARRAYS
  allocate (aiso_diverge(isd:ied,jsd:jed))
  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (aiso_back(isd:ied,jsd:jed,nk))
  allocate (aaniso_back(isd:ied,jsd:jed,nk))
  allocate (aiso(isd:ied,jsd:jed,nk)) 
  allocate (daur_dxur(isd:ied,jsd:jed))
  allocate (daur_dyur(isd:ied,jsd:jed))
  allocate (massqc(isd:ied,jsd:jed,nk))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (stress(isd:ied,jsd:jed,0:1,0:1,2))  !stress tensor: last index 1=stress_xx, 2=stress_xy
  allocate (tmpfdelx(isd:ied,jsd:jed,2))    
  allocate (tmpfdely(isd:ied,jsd:jed,2))   
  allocate (tmp(isd:ied,jsd:jed))
  allocate (cos2theta(isd:ied,jsd:jed))
  allocate (sin2theta(isd:ied,jsd:jed))
#endif

  allocate (grid_length(isd:ied,jsd:jed))
  allocate (viscosity_scaling(isd:ied,jsd:jed))
  grid_length(:,:)       = Grd%umask(:,:,1)*2.0*Grd%dxu(:,:)*Grd%dyu(:,:)/(epsln+Grd%dxu(:,:)+Grd%dyu(:,:))
  viscosity_scaling(:,:) = Grd%umask(:,:,1)

  tmp(:,:)          = 0.0
  aiso_diverge(:,:) = 0.0
  
  ! Some commonly used grid arrays 
  daur_dxur(:,:) = 0.0
  daur_dyur(:,:) = 0.0
  daur_dxur(:,:) = Grd%daur(:,:)*Grd%dxur(:,:)
  daur_dyur(:,:) = Grd%daur(:,:)*Grd%dyur(:,:)

  ! critical value of viscosity, above which 2-dim linear stability is not satisfied
  ! visc_crit taken from equation (18.17) of Griffies book.
  visc_crit(:,:) = 0.0
  do j=jsd,jed
    do i=isd,ied
      dxdy = 1.0/( 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j)) )
      visc_crit(i,j) = 0.5*dxdy/(dtime + epsln)
    enddo
  enddo
  
  ! may need a more conservative restriction in high latitudes when coupling to ice models 
  if(restrict_polar_visc) then 
     do j=jsd,jed
       do i=isd,ied
          if (abs(Grd%yu(i,j)) > restrict_polar_visc_lat) then    
            visc_crit(i,j) = restrict_polar_visc_ratio*visc_crit(i,j)
          endif
       enddo
     enddo
  endif

  id_visc_crit_lap = register_static_field ('ocean_model', 'visc_crit_lap', &
                     Grd%vel_axes_uv(1:2), 'critical viscosity',            &
                     'm^2/sec',missing_value=missing_value, range=(/0.0,1.e20/))
  call diagnose_2d_u(Time, Grd, id_visc_crit_lap, visc_crit(:,:))


  ! compute some smag fields and parameters 

  ! f_smag = (delta s * k_smag/pi)^2 where (delta s) is for velocity cell
  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = (k_smag_iso/pi)**2
  coeff_aniso      = (k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = coeff_iso  *((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**2
  fsmag_aniso(:,:) = coeff_aniso*((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**2

  ! remove smag within equatorial_zonal_lat band 
  if(equatorial_no_smag) then
    do j=jsd,jed
      do i=isd,ied
        if(abs(Grd%yu(i,j)) < equatorial_zonal_lat) then
          fsmag_iso(i,j)   = 0.0
          fsmag_aniso(i,j) = 0.0
        endif 
      enddo
    enddo 
  endif 

  ! compute time independent background viscosities 
  aiso(:,:,:)        = 0.0
  aiso_back(:,:,:)   = 0.0
  aaniso_back(:,:,:) = 0.0

  ! grid scale dependent "MICOM" background viscosities
  if(.not. viscosity_ncar)  then
      do j=jsd,jed
         do i=isd,ied
            aiso_back(i,j,:)     = vel_micom_iso  *(2.0*Grd%dxu(i,j)*Grd%dyu(i,j)) &
                                                   /(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
            aaniso_back(i,j,:)   = vel_micom_aniso*(2.0*Grd%dxu(i,j)*Grd%dyu(i,j)) &
                                                  /(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
            if(abs(Grd%yu(i,j)) < eq_lat_micom) then
                aiso_back(i,j,:)   = eq_vel_micom_iso  *(2.0*Grd%dxu(i,j)*Grd%dyu(i,j)) &
                                                       /(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
                aaniso_back(i,j,:) = eq_vel_micom_aniso*(2.0*Grd%dxu(i,j)*Grd%dyu(i,j)) &
                                                       /(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
            endif
         enddo
      enddo
  endif

  ! for divergence damping 
  do j=jsd,jed
     do i=isd,ied
        aiso_diverge(i,j) = divergence_damp_vel_micom*(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)) &
                            /(epsln + Grd%dxt(i,j) + Grd%dyt(i,j))
     enddo
  enddo


  ! (x,y,z) viscosities according to NCAR CCSM2.0 algorithm
  ! if use ncar only for equatorial region, revert to micom for off-equatorial 
  if(viscosity_ncar)  then
     call anisotropic_ncar
     if(ncar_only_equatorial) then    
        do j=jsd,jed
           do i=isd,ied  
              if(abs(Grd%yu(i,j)) > equatorial_zonal_lat) then
                 aiso_back(i,j,:) = vel_micom_iso  *(2.0*Grd%dxu(i,j)*Grd%dyu(i,j)) &
                                                   /(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
                 aaniso_back(i,j,:) = vel_micom_aniso*(2.0*Grd%dxu(i,j)*Grd%dyu(i,j)) &
                                                     /(epsln + Grd%dxu(i,j) + Grd%dyu(i,j))
              endif
           enddo
        enddo
     endif  
  endif 

  ! enhance background viscosity near open boundaries if needed
  if (have_obc) call ocean_obc_enhance_visc_back(aiso_back, aaniso_back)

  ! ensure background viscosities are not too large 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(aiso_back(i,j,k)   > visc_crit(i,j))  aiso_back(i,j,k)   = visc_crit(i,j)
           if(aaniso_back(i,j,k) > visc_crit(i,j))  aaniso_back(i,j,k) = visc_crit(i,j)
        enddo
     enddo
  enddo

  ! ensure viscosities maintain proper relative values so friction dissipates 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
            if(aaniso_back(i,j,k) >= 2.0*aiso_back(i,j,k) .and. aaniso_back(i,j,k) > 0.0)  then 
              write(stdoutunit,'(a,i4,a,i4,a,i3,a)') &
              'Violating lap iso/aniso constraint at (',i,',',j,',',k,')'
              aaniso_back(i,j,k) = 1.9*aiso_back(i,j,k)
            endif 
         enddo
      enddo
  enddo

  call compute_neptune_velocity(Time)


  if(use_side_drag_friction) then 
     write(stdoutunit,'(a)') '==> Note from ocean_lapgen_friction: using side drag friction for cells next to land.'
     allocate (umask_next_to_land(isd:ied,jsd:jed,nk))
     umask_next_to_land(:,:,:) = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              umask_next_to_land(i,j,k) = (1.0-Grd%umask(i+1,j,k)) + (1.0-Grd%umask(i-1,j,k)) &
                                         +(1.0-Grd%umask(i,j+1,k)) + (1.0-Grd%umask(i,j-1,k)) 
              umask_next_to_land(i,j,k) = Grd%umask(i,j,k)*min(1.0,umask_next_to_land(i,j,k))
           enddo
        enddo
     enddo 
     id_umask_next_to_land = register_static_field('ocean_model', 'umask_next_to_land_lap',&
            Grd%vel_axes_uv(1:3),'U-cell mask for cells next to land for lapgen module',   &
            'dimensionless', missing_value=missing_value, range=(/-10.0,10.0/))
     call diagnose_3d_u(Time, Grd, id_umask_next_to_land, umask_next_to_land(:,:,:))
  endif 


  ! send static fields to diag_manager 
  id_along_back   = register_static_field('ocean_model', 'along_lap_back',      &
                    Grd%vel_axes_uv(1:3),'U-cell background along-stream visc', &
                    'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))
  id_across_back = register_static_field('ocean_model', 'across_lap_back',      &
                   Grd%vel_axes_uv(1:3), 'U-cell background cross-stream visc', &
                   'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))

  if (id_along_back > 0) then 
     call diagnose_3d_u(Time, Grd, id_along_back, aiso_back(:,:,:)+0.5*aaniso_back(:,:,:))
  endif 
  if (id_across_back > 0) then 
     call diagnose_3d_u(Time, Grd, id_across_back, aiso_back(:,:,:)-0.5*aaniso_back(:,:,:))
  endif 

  id_aiso_diverge = register_static_field('ocean_model', 'aiso_lap_diverge', &
              Grd%tracer_axes(1:2), 'U-cell visc for divergence damping',    &
              'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_2d(Time, Grd, id_aiso_diverge, aiso_diverge(:,:))

  ! other viscosities for diagnostic output 

  id_aiso   = register_diag_field ('ocean_model', 'aiso_lap', Grd%vel_axes_uv(1:3), &
              Time%model_time, 'U-cell isotropic visc', 'm^2/sec',                  &
              missing_value=missing_value, range=(/-10.0,1.e10/),                   &
              standard_name='ocean_momentum_xy_laplacian_diffusivity')
  id_aaniso = register_diag_field ('ocean_model', 'aaniso_lap', Grd%vel_axes_uv(1:3), &
              Time%model_time, 'U-cell anisotropic visc', 'm^2/sec',                  &
              missing_value=missing_value, range=(/-10.0,1.e10/))
  id_along  = register_diag_field ('ocean_model', 'along', Grd%vel_axes_uv(1:3), &
              Time%model_time, 'U-cell along-stream visc', 'm^2/sec',            &
              missing_value=missing_value, range=(/-10.0,1.e10/))
  id_across = register_diag_field ('ocean_model', 'across', Grd%vel_axes_uv(1:3), &
              Time%model_time, 'U-cell cross-stream visc', 'm^2/sec',             &
              missing_value=missing_value, range=(/-10.0,1.e10/))
  id_lap_fric_u = register_diag_field ('ocean_model', 'lap_fric_u', Grd%vel_axes_uv(1:3),              &
                  Time%model_time, 'Thick & rho wghtd horz lap frict on u-zonal','(kg/m^3)*(m^2/s^2)', &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_lap_fric_v = register_diag_field ('ocean_model', 'lap_fric_v', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thick & rho wghtd horz lap frict on v-merid', '(kg/m^3)*(m^2/s^2)', &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_horz_lap_diss = register_diag_field ('ocean_model', 'horz_lap_diss', Grd%vel_axes_uv(1:3),     &
                  Time%model_time, 'Energy dissipation from horizontal Laplacian friction', 'W/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_viscosity_scaling = register_diag_field ('ocean_model', 'viscosity_scaling', Grd%vel_axes_uv(1:2),  &
              Time%model_time, 'Grid/Rossby radius scaling for the laplacian viscosity', 'dimensionless',&
              missing_value=missing_value, range=(/-1.0,10.0/))

  id_side_drag_friction_u = register_diag_field ('ocean_model', 'side_drag_friction_lap_u', Grd%vel_axes_uv(1:3),&
                  Time%model_time, 'u-friction due to side drag from lapgen module', 'N/m^2',                    &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_side_drag_friction_v = register_diag_field ('ocean_model', 'side_drag_friction_lap_v', Grd%vel_axes_uv(1:3),&
                  Time%model_time, 'v-friction due to side drag from lapgen module', 'N/m^2',                    &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_lap_plus_side_fric_u = register_diag_field ('ocean_model', 'lap_plus_side_fric_u', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz lap frict + side friction on u-zonal',       &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_lap_plus_side_fric_v = register_diag_field ('ocean_model', 'lap_plus_side_fric_v', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz lap frict + side friction on v-merid',       &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))


end subroutine ocean_lapgen_friction_init
! </SUBROUTINE>  NAME="ocean_lapgen_friction_init"


!#######################################################################
! <SUBROUTINE NAME="lapgen_friction">
!
! <DESCRIPTION>
! This routine computes thickness weighted and density weighted 
! time tendency for horizontal velocity arising from horizontal 
! Laplacian friction.  
!
! The algorithm is derived from a functional approach that ensures
! kinetic energy is consistenty dissipated for all flow configurations. 
! The triad do-loops are expanded in order to enhance the 
! ability of cache-based machines to keep most of the variables 
! on-cache.  
! 
! Fundamental to the scheme are the rates of horizontal deformation  <BR/> 
! horizontal tension = DT = (dy)(u/dy)_x - (dx)(v/dx)_y              <BR/>    
! horizontal strain  = DS = (dx)(u/dx)_y + (dy)(v/dy)_x              <BR/> 
! Units of the tension and strain are sec^-1.
!
! Four tensions and four strains are computed for each velocity point, <BR/>
! corresponding to the four triads surrounding the point.              <BR/>  
! The following notation is used to distinguish the triads:            <BR/>  
! (0,1)=northwest triad  (1,1)=northeast triad,                        <BR/> 
! (0,0)=southwest triad, (1,0)=southeast triad
!
! A triad contributes when at least one of its velocities is            
! not a land point.  In order to obtain the correct tension          
! and strain next to boundaries, tension and strain should not be   
! masked with umask. 
!
! </DESCRIPTION>
!
subroutine lapgen_friction(Time, Thickness, Adv_vel, Velocity, &
                           lap_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step
  real, dimension(isd:,jsd:), intent(inout) :: lap_viscosity

  real, dimension(0:1,0:1) :: aiso_smag
  real, dimension(0:1,0:1) :: aaniso_smag
  real, dimension(0:1,0:1) :: dissipate

  integer :: i, j, k, n
  integer :: ip, jq, jjq, iip
  integer :: taum1, tau

  real :: u1_m10, u2_m10
  real :: u1_10,  u2_10
  real :: u1_00,  u2_00
  real :: u1_01,  u2_01
  real :: u1_0m1, u2_0m1 
  real :: dxuer_ip0, dxuer_ip1
  real :: dyunr_jq0, dyunr_jq1
  real :: usqrd, vsqrd, umagr
  real :: uvmag, umag, vmag, velocity_gamma
  real :: tension, strain, deform, delta 
  real :: tension_metric, strain_metric
  real,dimension (isd:ied,jsd:jed,8,nk) :: stress_B

  integer :: stdoutunit, ibl
  integer, dimension(nk) :: id_update, kstart, kend

  stdoutunit=stdout() 

  if(.not. use_this_module) then

      if(energy_analysis_step) then 
          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      Velocity%wrkv(i,j,k,n) = 0.0
                   enddo
                enddo
             enddo
          enddo
      endif

    return 
  endif 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapgen_friction_mod (lapgen_friction): module not initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  ! to scale down the viscosity if the Rossby radius is resolved by the grid 
  if(viscosity_scale_by_rossby) then 
     do j=jsc,jec
        do i=isc,iec
           wrk1_2d(i,j) = Velocity%rossby_radius(i,j)
        enddo
     enddo
     call mpp_update_domains (wrk1_2d(:,:), Dom%domain2d)    
     do j=jsd,jed
        do i=isd,ied
           viscosity_scaling(i,j) = grid_length(i,j)**viscosity_scale_by_rossby_power &
           /(grid_length(i,j)**viscosity_scale_by_rossby_power + wrk1_2d(i,j)**viscosity_scale_by_rossby_power + epsln)
        enddo
     enddo
  else 
     viscosity_scaling(:,:) = Grd%umask(:,:,1)
  endif 


  stress(:,:,:,:,:) = 0.0
  wrk1_v(:,:,:,:)   = 0.0
  wrk2_v(:,:,:,:)   = 0.0
  wrk3_v(:,:,:,:)   = 0.0
  wrk1(:,:,:)       = 0.0 
  wrk2(:,:,:)       = 0.0 

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           massqc(i,j,k) = 0.25*Grd%dau(i,j)*Thickness%rho_dzu(i,j,k,tau)
        enddo
     enddo
  enddo 

  kstart = 1
  kend   = nk
  ibl    = 1
  do k=1,nk

     do j=jsd,jec
        do i=isd,iec

           tension_metric = -(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))*Grd%dh2dx(i,j) &
                            +(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*Grd%dh1dy(i,j)  

           strain_metric  = -(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))*Grd%dh1dy(i,j) &
                            -(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*Grd%dh2dx(i,j)  

           if(equatorial_zonal .and. abs(Grd%yu(i,j)) <= equatorial_zonal_lat) then  
             sin2theta(i,j) = 0.0
             cos2theta(i,j) = 1.0
           else 
             usqrd = (Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))**2
             vsqrd = (Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))**2
             umagr = 1.0/(epsln + usqrd + vsqrd)
             sin2theta(i,j) = 2.0*(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1)) &
                                 *(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*umagr
             cos2theta(i,j) = (usqrd-vsqrd)*umagr
           endif  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are quicker. 

           ip=0 ; iip = max(isd,i+ip-1)
           dxuer_ip0 = Grd%dxuer(iip,j)

           ip=1 
           dxuer_ip1 = Grd%dxuer(i+ip-1,j)

           jq=0 ; jjq = max(jsd,j+jq-1)
           dyunr_jq0 = Grd%dyunr(i,jjq)

           jq=1 
           dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; iip = max(i+ip,isd)
           u1_m10 = Velocity%u(iip,j+jq,k,1,taum1) - neptune_velocity(iip,j+jq,1) 
           u2_m10 = Velocity%u(iip,j+jq,k,2,taum1) - neptune_velocity(iip,j+jq,2) 

           ip=1  ; jq=0  
           u1_10  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1) 
           u2_10  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)

           ip=0  ; jq=0  
           u1_00  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1) 
           u2_00  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)

           ip=0  ; jq=1  
           u1_01  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)  
           u2_01  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2) 

           ip=0  ; jq=-1;  jjq = max(j+jq,jsd)
           u1_0m1 = Velocity%u(i+ip,jjq,k,1,taum1) - neptune_velocity(i+ip,jjq,1) 
           u2_0m1 = Velocity%u(i+ip,jjq,k,2,taum1) - neptune_velocity(i+ip,jjq,2)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = viscosity_scaling(i,j)*min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = viscosity_scaling(i,j)*min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)
           dissipate(ip,jq)    = -aiso_smag(ip,jq)*deform**2 + 2.0*aaniso_smag(ip,jq)*delta**2

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = viscosity_scaling(i,j)*min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = viscosity_scaling(i,j)*min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)
           dissipate(ip,jq)    = -aiso_smag(ip,jq)*deform**2 + 2.0*aaniso_smag(ip,jq)*delta**2

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = viscosity_scaling(i,j)*min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = viscosity_scaling(i,j)*min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)
           dissipate(ip,jq)    = -aiso_smag(ip,jq)*deform**2 + 2.0*aaniso_smag(ip,jq)*delta**2

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq)    = viscosity_scaling(i,j)*min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(ip,jq)  = viscosity_scaling(i,j)*min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           stress(i,j,ip,jq,1) = aiso_smag(ip,jq)*tension + aaniso_smag(ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(ip,jq)*strain  - aaniso_smag(ip,jq)*delta*cos2theta(i,j)
           dissipate(ip,jq)    = -aiso_smag(ip,jq)*deform**2 + 2.0*aaniso_smag(ip,jq)*delta**2

           ! viscosities and dissipation for diagnostics 
           aiso(i,j,k) = 0.25*(aiso_smag(0,0)   + aiso_smag(1,0)   + aiso_smag(0,1)   + aiso_smag(1,1))
           wrk1(i,j,k) = 0.25*(aaniso_smag(0,0) + aaniso_smag(1,0) + aaniso_smag(0,1) + aaniso_smag(1,1))
           wrk2(i,j,k) = massqc(i,j,k)*(dissipate(0,0) + dissipate(1,0) + dissipate(0,1) + dissipate(1,1))

        enddo
     enddo
     stress_B(:,:,1,k)  = stress(:,:,0,0,1)
     stress_B(:,:,2,k)  = stress(:,:,1,0,1)
     stress_B(:,:,3,k)  = stress(:,:,0,1,1)
     stress_B(:,:,4,k)  = stress(:,:,1,1,1)
     stress_B(:,:,5,k)  = stress(:,:,0,0,2)
     stress_B(:,:,6,k)  = stress(:,:,1,0,2)
     stress_B(:,:,7,k)  = stress(:,:,0,1,2)
     stress_B(:,:,8,k)  = stress(:,:,1,1,2)
     if( async_domain_update ) then
        if(mod(k,blocksize) == 0 .or. k == nk)   then
           kend(ibl) = k
           id_update(ibl) = mpp_start_update_domains (stress_B(:,:,:,kstart(ibl):kend(ibl)), Dom%domain2d)  
           ibl = ibl+1
           if(ibl<nk) kstart(ibl) = k+1
        endif
     endif
  enddo

  if( .not. async_domain_update ) then
     call mpp_update_domains (stress_B(:,:,:,:), Dom%domain2d)    
  endif  

  ibl = 1
  do k=1,nk
     if( async_domain_update ) then
        if(k == kstart(ibl))   then
           call mpp_complete_update_domains (id_update(ibl), stress_B(:,:,:,kstart(ibl):kend(ibl)), Dom%domain2d) 
           ibl = ibl+1
        endif
     endif
     stress(:,:,0,0,1) = stress_B(:,:,1,k)
     stress(:,:,1,0,1) = stress_B(:,:,2,k)
     stress(:,:,0,1,1) = stress_B(:,:,3,k)
     stress(:,:,1,1,1) = stress_B(:,:,4,k)
     stress(:,:,0,0,2) = stress_B(:,:,5,k)
     stress(:,:,1,0,2) = stress_B(:,:,6,k)
     stress(:,:,0,1,2) = stress_B(:,:,7,k)
     stress(:,:,1,1,2) = stress_B(:,:,8,k)
     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains

     tmpfdelx(:,:,:) = 0.0
     do j=jsc,jec
        do i=isd,iec
           tmpfdelx(i,j,1) = (stress(i,j,1,0,1)   + stress(i,j,1,1,1)  )*massqc(i,j,k)    &
                            +(stress(i+1,j,0,0,1) + stress(i+1,j,0,1,1))*massqc(i+1,j,k)
           tmpfdelx(i,j,2) = (stress(i,j,1,0,2)   + stress(i,j,1,1,2)  )*massqc(i,j,k)    &  
                            +(stress(i+1,j,0,0,2) + stress(i+1,j,0,1,2))*massqc(i+1,j,k)              
        enddo
     enddo

     tmpfdely(:,:,:) = 0.0
     do j=jsd,jec
        do i=isc,iec  
           tmpfdely(i,j,1) = (stress(i,j,0,1,2)   + stress(i,j,1,1,2)  )*massqc(i,j,k)    &
                            +(stress(i,j+1,0,0,2) + stress(i,j+1,1,0,2))*massqc(i,j+1,k) 
           tmpfdely(i,j,2) = (stress(i,j,0,1,1)   + stress(i,j,1,1,1)  )*massqc(i,j,k)    &
                            +(stress(i,j+1,0,0,1) + stress(i,j+1,1,0,1))*massqc(i,j+1,k)
        enddo
     enddo


     ! thickness weighted and density weighted frictional acceleration (kg/m^3)*(m/sec^2) 
     do n=1,2

        tmp(:,:) =  BDX_EU_smag(tmpfdelx(:,:,n))+(3-2*n)*BDY_NU_smag(tmpfdely(:,:,n))    
        do j=jsc,jec
           do i=isc,iec
              wrk1_v(i,j,k,n) = tmp(i,j)*Grd%umask(i,j,k)
           enddo
        enddo

        ! reduce to 5-point laplacian at bottom to avoid problems with thin bottom partial cells
        if (bottom_5point) then
            tmp(:,:) = BDX_EU(aiso_back(:,:,k)*FMX(Thickness%rho_dzu(:,:,k,tau)) &
                      *FDX_U(Velocity%u(:,:,k,n,taum1)-neptune_velocity(:,:,n))) &
                      +BDY_NU(aiso_back(:,:,k)*FMY(Thickness%rho_dzu(:,:,k,tau)) &
                      *FDY_U(Velocity%u(:,:,k,n,taum1)-neptune_velocity(:,:,n)))
            do j=jsc,jec
               do i=isc,iec
                  if(k==Grd%kmu(i,j)) wrk1_v(i,j,k,n) = tmp(i,j)*Grd%umask(i,j,k)*viscosity_scaling(i,j)
               enddo
            enddo
        endif

     enddo ! end of n-loop

     ! use side drag for cells next to side boundaries 
     if(use_side_drag_friction) then
        do j=jsc,jec
           do i=isc,iec
              uvmag           = sqrt(Velocity%u(i,j,k,1,taum1)**2 + Velocity%u(i,j,k,2,taum1)**2)
              uvmag           = min(uvmag,side_drag_friction_uvmag_max)
              velocity_gamma  = rho0*side_drag_friction_scaling*Velocity%cdbot_array(i,j) &
                                *uvmag*umask_next_to_land(i,j,k)
              umag            = abs(Velocity%u(i,j,k,1,taum1))
              vmag            = abs(Velocity%u(i,j,k,2,taum1))
              wrk2_v(i,j,k,1) = -min(side_drag_friction_max, velocity_gamma*umag)*sign(1.0,Velocity%u(i,j,k,1,taum1))
              wrk2_v(i,j,k,2) = -min(side_drag_friction_max, velocity_gamma*vmag)*sign(1.0,Velocity%u(i,j,k,2,taum1))
              wrk3_v(i,j,k,1) =  wrk1_v(i,j,k,1)*(1.0-umask_next_to_land(i,j,k))
              wrk3_v(i,j,k,2) =  wrk1_v(i,j,k,2)*(1.0-umask_next_to_land(i,j,k))
              wrk1_v(i,j,k,1) =  wrk2_v(i,j,k,1) + wrk3_v(i,j,k,1)
              wrk1_v(i,j,k,2) =  wrk2_v(i,j,k,2) + wrk3_v(i,j,k,2)              
            enddo
         enddo
     endif 

  enddo  !end of k-loop


  ! damping on the divergence field   
  if(divergence_damp) then 
      do n=1,2 
         do k=1,nk
            wrk1_v2d(:,:,:) = 0.0
            wrk1_v2d(:,:,1) = FAY(FDX_T(aiso_diverge(:,:)*Adv_vel%diverge_t(:,:,k)))
            wrk1_v2d(:,:,2) = FAX(FDY_T(aiso_diverge(:,:)*Adv_vel%diverge_t(:,:,k)))
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v(i,j,k,n) = wrk1_v(i,j,k,n) + wrk1_v2d(i,j,n)*Grd%umask(i,j,k)
               enddo
            enddo
         enddo
      enddo
  endif


  if(energy_analysis_step) then 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%wrkv(i,j,k,n) = wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

  else 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      ! vertically averaged isotropic portion of the viscosity 
      lap_viscosity(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               lap_viscosity(i,j) = lap_viscosity(i,j) &
                 + Grd%umask(i,j,k)*Thickness%rho_dzu(i,j,k,tau)*aiso(i,j,k)
            enddo
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            lap_viscosity(i,j) = &
            Grd%umask(i,j,1)*lap_viscosity(i,j)/(epsln+Thickness%mass_u(i,j,tau))
         enddo
      enddo
      call mpp_update_domains (lap_viscosity(:,:), Dom%domain2d)    


      call diagnose_3d_u(Time, Grd, id_aiso, aiso(:,:,:))
      call diagnose_3d_u(Time, Grd, id_aaniso, wrk1(:,:,:))
      if (id_along > 0)  call diagnose_3d_u(Time, Grd, id_along,  aiso(:,:,:)+0.5*wrk1(:,:,:))
      if (id_across > 0) call diagnose_3d_u(Time, Grd, id_across, aiso(:,:,:)-0.5*wrk1(:,:,:))

      call diagnose_3d_u(Time, Grd, id_lap_fric_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_lap_fric_v, wrk1_v(:,:,:,2))

      if (id_horz_lap_diss > 0) then 
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk2(i,j,k) = wrk2(i,j,k)*Grd%daur(i,j)
                enddo
             enddo
          enddo
          call diagnose_3d_u(Time, Grd, id_horz_lap_diss, wrk2(:,:,:))
      endif 

      call diagnose_2d_u(Time, Grd, id_viscosity_scaling, viscosity_scaling(:,:))

      if(use_side_drag_friction) then 
         call diagnose_3d_u(Time, Grd, id_lap_fric_u, wrk3_v(:,:,:,1))
      else 
         call diagnose_3d_u(Time, Grd, id_lap_fric_u, wrk1_v(:,:,:,1))
      endif

      if(use_side_drag_friction) then 
         call diagnose_3d_u(Time, Grd, id_lap_fric_v, wrk3_v(:,:,:,2))
      else 
         call diagnose_3d_u(Time, Grd, id_lap_fric_v, wrk1_v(:,:,:,2))
      endif

      call diagnose_3d_u(Time, Grd, id_lap_plus_side_fric_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_lap_plus_side_fric_v, wrk1_v(:,:,:,2))
      call diagnose_3d_u(Time, Grd, id_side_drag_friction_u, wrk2_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_side_drag_friction_v, wrk2_v(:,:,:,2))

  endif

  if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_lapgen_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('lapgen friction(1)', wrk1_v(COMP,:,1)*Grd%umask(COMP,:))
      call write_chksum_3d('lapgen friction(2)', wrk1_v(COMP,:,2)*Grd%umask(COMP,:))
  endif


end subroutine lapgen_friction
! </SUBROUTINE> NAME="lapgen_friction"


!#######################################################################
! <FUNCTION NAME="BDX_EU_smag">
!
! <DESCRIPTION>
! Compute backwards derivative in X of a quantity defined on the east 
! face of a U-cell. Slightly modified version of BDX_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i-1/2,j).
!
! BDX_EU_smag(a) has dimensions of a*m^-3 
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! field defined on the east face of a U-cell
! </IN>
!
function BDX_EU_smag(a)

  real, dimension(isd:,jsd:), intent(in) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDX_EU_smag
  integer :: i, j

  do j=jsd,jed
    do i=isd+1,ied
      BDX_EU_smag(i,j) = &
      (Grd%dyue_dxuer(i,j)*a(i,j) - Grd%dyue_dxuer(i-1,j)*a(i-1,j))*daur_dyur(i,j)
    enddo
    BDX_EU_smag(isd,j) = 0.0
  enddo

end function BDX_EU_smag
! </FUNCTION> NAME="BDX_EU_smag"


!#######################################################################
! <FUNCTION NAME="BDY_NU_smag">
!
! <DESCRIPTION>
! Compute backwards derivative in Y of a quantity defined on the north
! face of a U-cell. Slightly modified version of BDY_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i,j-1/2)
!
! BDY_NU_smag(a) has dimensions of a*m^-3 
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! field defined on the north face of a U-cell
! </IN>
!
function BDY_NU_smag(a)

  real, dimension(isd:,jsd:), intent(in) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDY_NU_smag
  integer :: i, j

  do j=jsd+1,jed
    do i=isd,ied            
      BDY_NU_smag(i,j) = &
      (Grd%dxun_dyunr(i,j)*a(i,j) - Grd%dxun_dyunr(i,j-1)*a(i,j-1))*daur_dxur(i,j)
    enddo
  enddo
  BDY_NU_smag(:,jsd) = 0.0

end function BDY_NU_smag
! </FUNCTION> NAME="BDY_EU_smag"


!#######################################################################
! <SUBROUTINE NAME="anisotropic_ncar">
!
! <DESCRIPTION>
!
!     Spatially-varying anisotropic viscosity initialization
!
!     This routine defines NCOM-like spatial distributions of
!     viscosity coefficients F_PARA and F_PERP.
!     Uses NCAR CCSM2.0 algorithm with cm^2/sec --> m^2/sec.  
!
!     written by:     Stephen Yeager 3/2000                           <BR/> 
!
!     modified by:    Gokhan Danabasoglu (08/2001)                    <BR/>
!
!     port to mom4:   Stephen.Griffies (9/2002)  
!
!     update to mom4p1 based on new tunes from NCAR
!                    Stephen.Griffies (7/2007)  
!
!
!   "A_viscosity" = F_PARA = Along = viscosity parallel to flow 
!                  = max{0.5*visc_vel_scale(z)*A*max[dx,dy],vconst_6}
!
!   where                                                                <BR/>    
!          A = 0.425 * cos(pi*y*radian/30) + 0.575   for |y*radian| < 30 <BR/>
!          A = 0.15                                  otherwise 
!
!   Here, A provides a horizontal variation for visc_vel_scale.
!
!   "B_viscosity" = F_PERP = Across = viscosity perpendicular to flow = max( bu, bv)
!
!   and                                                                  <BR/>    
!        F_PARA = min(F_PARA, AMAX_CFL),                                 <BR/> 
!        F_PERP = min(F_PERP, AMAX_CFL),                                 <BR/>   
!        F_PARA = max(F_PARA, F_PERP)                                    <BR/> 
!   are enforced 
!
!   In the above equations, 
!
!        bu  = vconst_1 * ( 1 + vconst_2  * ( 1 + cos( 2*y + pi ) ) )        <BR/>
!        bv  = vconst_3 * beta_f * dx^3   * exp( - (vconst_4 * distance)^2 ) <BR/>
!
!   with                                                                     <BR/>        
!        beta_f         (x,y)   = 2 * omega_earth* cos(ULAT(i,j)) / radius   <BR/>  
!        distance       (x,y,z) = actual distance to "vconst_5" points       <BR/>  
!                                 west of the nearest western boundary       <BR/>  
!        dx             (x,y)   = DXU(i,j)                                   <BR/>  
!        dy             (x,y)   = DYU(i,j)                                   <BR/> 
!        visc_vel_scale (z)     = vconst_7 * exp(-zt(k)/visc_vel_scale_length)  <BR/> 
!        visc_vel_scale_length  = e-folding scale ( default = 1500.0e2 cm)      <BR/> 
!        y              (x,y)   = ULAT(i,j), latitude of "u/v" grid pts in radians   <BR/> 
!        In MOM, ULAT(radians) = xu*pi/180 with xu(i,j) the longitude of U grid points in degrees
!
!   "vconst_#" are input parameters defined in namelist ocean_lapgen_friction_general_nml. 
!   "vconst_1", "vconst_6", and "vconst_4" have dimensions of cm^2/s,
!   cm^2/s, and 1/cm, respectively. "vconst_5" is an INTEGER.
!
!   NOTE: The nearest western boundary computations are done along the
!         model i-grid lines. Therefore, viscosity based on these are 
!         only approximate in the high Northern Hemisphere when using 
!         generalized coordinates with coordinate pole(s) shifted onto land. 
!
! </DESCRIPTION>
!
subroutine anisotropic_ncar

  integer :: i, j, k
  integer :: n, ncount
  integer :: is, ie, ip1, ii
  integer :: index, indexo
  integer, dimension(nx) :: iwp
  integer, dimension(nx) :: nwbp_global

  real, dimension(nx,jsc:jec)            :: dxtn_global
  real, dimension(nx,jsc:jec)            :: kmu_global
  real, dimension(isd:ied,jsd:jed)  :: kmu_tmp
  real, dimension(nx)               :: dist_g 

  real :: dist_max, vvsl  
  real :: a_para, bu, bv 
  real :: f_para, f_perp
  real :: visc_vel_scale

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapgen_friction_mod (anisotropic_ncar): module must be initialized')
  endif 

  dist_max = 1.e10                  ! distance for ACC region (cm)
  vvsl     = visc_vel_scale_length  ! exponential decay scale (cm) for decay away from ocean surface

  dxtn_global(:,:) = 0.0
  kmu_global(:,:)  = 0.0
  kmu_tmp(:,:)     = Grd%kmu(:,:)

  call mpp_global_field(Dom%domain2d,kmu_tmp,kmu_global, flags = XUPDATE)  
  call mpp_global_field(Dom%domain2d,Grd%dxtn,dxtn_global, flags = XUPDATE)

  do k=1,nk

     ! 100.0 factor in exponential is to convert MOM zt from (m) to (cm)
     ! visc_vel_scale is in cm/sec.  It is used for along-stream viscosity. 
     visc_vel_scale = vconst_7 * exp(-100.0*Grd%zt(k)/vvsl)   

     do j=jsc,jec


        ! determine nearest western boundary
        ncount         = 0
        nwbp_global(:) = 0
        iwp(:)         = 0
        dist_g(:)      = 0.0

        do i=1,Grd%ni
           ip1 = i+1
           if(i==Grd%ni) then 
               if(Grd%cyclic_x) then 
                   ip1=1
               else
                   ip1=i  
               endif
           endif
           if(kmu_global(i,j)<k .and. kmu_global(ip1,j) >= k) then 
               ncount      = ncount+1
               iwp(ncount) = i
           endif
        enddo

        if (ncount > 0) then
            do n=1,ncount-1
               is = iwp(n)
               ie = iwp(n+1)-1
               do i=is,ie
                  nwbp_global(i)=is
               enddo
            enddo
            do i=1,Grd%ni
               if(nwbp_global(i)==0) nwbp_global(i) = iwp(ncount)
            enddo
        endif


        ! determine distance (cm) to nearest western boundary
        do i=1,Grd%ni

           index  = nwbp_global(i)
           indexo = index + vconst_5
           if ( index .eq. 0) then
               dist_g(i) = dist_max
           elseif ( i .ge. index  .and.  i .le. indexo ) then
               dist_g(i) = 0.0
           elseif ( (i .gt. indexo) ) then
               dist_g(i) = dxtn_global(i,j)+dist_g(i-1)
           elseif ( i .lt. index ) then
               if (indexo .le. Grd%ni) then
                   if (i .eq. 1) then
                       dist_g(i) = 0.0
                       do ii=indexo+1,Grd%ni
                          dist_g(i)=dxtn_global(ii,j) + dist_g(i)
                       enddo
                       dist_g(i) = dxtn_global(i,j)+dist_g(i)
                   else
                       dist_g(i) = dxtn_global(i,j)+dist_g(i-1)
                   endif
               else
                   if (i .le. (indexo - Grd%ni)) then
                       dist_g(i) = 0.0
                   else
                       dist_g(i) = dxtn_global(i,j)+dist_g(i-1)
                   endif
               endif
           endif

        enddo

        ! convert distance in MOM (m) to POP1.4 (cm)
        dist_g(:) = 100.0*dist_g(:)

 
        !calculate viscosity over comp domain 
        do i=isc,iec

           
           if(viscosity_ncar_2000) then 

               ! f_para is viscosity parallel to flow 
               ! 100.0 factor in f_para converts (m) to (cm)
               a_para = 0.15
               if ( abs(Grd%yu(i,j)) < 30.0) then  ! yu is in degrees, not radians 
                   a_para = 0.425*cos(pi*Grd%yu(i,j)/30.0) + 0.575
               endif
               f_para = &
                    max(0.5 * visc_vel_scale * a_para * 100.0*max(Grd%dxu(i,j),Grd%dyu(i,j)), vconst_6 )

               ! f_perp is viscosity perpendicular to flow 
               ! 1e4 factor in bv converts (m) to (cm) for dxu and beta_f 
               bu = vconst_1*(1.0 + vconst_2*(1.0 + cos((2.0*Grd%yu(i,j)/radian)+pi)))
               bv = 1e4*vconst_3*Grd%beta(i,j)*Grd%dxu(i,j)**3
               bv = bv*exp(-(vconst_4*dist_g(i+Dom%ioff))**2)
               f_perp = min(max(bu,bv),0.99*f_para)

           elseif(viscosity_ncar_2007) then 

               ! bv here is a temporary 
               bv = (min(abs(Grd%yu(i,j)),vconst_8)*90.0 / vconst_8) / radian 
               bu = vconst_1 * ( 1.0 + vconst_2*(1.0-cos(2.0*bv)))

               bv = 1e4*vconst_3*Grd%beta(i,j)*Grd%dxu(i,j)**3
               bv = bv*exp(-(vconst_4*dist_g(i+Dom%ioff))**2)
               f_para = max(bv,vconst_6)
               f_perp = min(max(bu,bv),0.99*f_para)

           endif

           if(debug_ncar_A) then 
             f_perp=f_para
           endif 
           if(debug_ncar_B) then 
             f_para=f_perp
           endif 

           ! revert to isotropic background viscosity outside equatorial_zonal_lat
           if(ncar_isotropic_off_equator .and. abs(Grd%yu(i,j)) >= equatorial_zonal_lat) then 
             f_perp=f_para
           endif   

           ! 1.e-4 converts from cm^2/sec to m^2/sec
           aiso_back(i,j,k)   = 1.e-4*0.5*(f_para + f_perp) 
           aaniso_back(i,j,k) = 1.e-4*(f_para - f_perp)

           ! revert to isotropic background viscosity beneath some level
           if(ncar_isotropic_at_depth .and. Grd%zt(k) >= ncar_isotropic_depth) then 
             aaniso_back(i,j,k) = 0.0
             aiso_back(i,j,k)   = max(aiso_back(i,j,k), ncar_isotropic_at_depth_visc)
           endif   


        enddo !i

     enddo ! j
  enddo    ! k

  ! need these viscosities in the halo regions
  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


end subroutine anisotropic_ncar
! </SUBROUTINE>  NAME="anisotropic_ncar"


!#######################################################################
! <SUBROUTINE NAME="lapgen_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the Laplacian 
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine lapgen_viscosity_check

  integer       :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: fudge

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL,&
     '==>Error in ocean_lapgen_friction_mod (lapgen_viscosity_check): needs initialization')
  endif 

  fudge=1.001  ! to eliminate roundoffs causing aiso=visc_crit+epsln

  write (stdoutunit,'(/60x,a/)') &
  ' Excessive horizontal Laplacian friction summary:'
  write (stdoutunit,'(1x,a/)') &
  'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then 
          if (aiso(i,j,k) > fudge*visc_crit(i,j) .and. num < max_num) then
            num = num + 1
            write (unit,9600) &
            'aiso(',aiso(i,j,k), visc_crit(i,j), i, j, k, &
            Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum(num)
  if (num > max_num) write (stdoutunit,*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es16.8,' m^2/s) exceeds max value (',es16.8,') at (i,j,k) = ', &
            '(',i4,',',i4,',',i4,'),',' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')


end subroutine lapgen_viscosity_check
! </SUBROUTINE> NAME="lapgen_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="lapgen_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the LLaplacian grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
! <IN NAME="u" TYPE="real" DIM="(isd:ied,jsd:jed,nk,2)">
! Horizontal velocity field at time tau
! </IN>
!
subroutine lapgen_reynolds_check(Time, Velocity)

  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  real :: rame, ramn, reyx, reyy
  real :: reynx0, reyny0, reynx,reyny, reynu,reynv, reynmu,reynmv
  integer :: ireynx,jreynx,kreynx, ireyny,jreyny,kreyny
  integer :: i, j, k, tau
  integer, save :: unit=6

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapgen_friction_mod (lapgen_reynolds_check): needs initialization')
  endif 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(aiso(i,j,k) + epsln)
        ramn = 1.0/(aiso(i,j,k) + epsln)
        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j))*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j))*ramn
        if (reyy > reyny) then
          ireyny = i
          jreyny = j
          kreyny = k
          reyny  = reyy
          reynv  = Velocity%u(i,j,k,2,tau)
          reynmv = 1.0/ramn
        endif
      enddo
    enddo
  enddo
  write (stdoutunit,'(/60x,a/)') ' Horizontal Laplacian Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0 .and. reynx > 1.e-6) then
    write (unit,10300) &
    reynx, ireynx, jreynx, kreynx, Grd%xu(ireynx,jreynx), &
    Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0 .and. reyny > 1.e-6) then
    write (unit,10400) &
    reyny, ireyny, jreyny, kreyny, Grd%xu(ireyny,jreyny), &
    Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine lapgen_reynolds_check
! </SUBROUTINE> NAME="lapgen_reynolds_check"


!#######################################################################
! <SUBROUTINE NAME="compute_neptune_velocity">
!
! <DESCRIPTION>
! Compute Neptune velocity.  
!
! Method follows that used in MOM2 and MOM3 as implemented by 
! Greg Holloway (zounds@ios.bc.ca) and Michael Eby (eby@uvic.ca) 
! Coded in mom4 by Stephen.Griffies 
!
! Neptune is calculated as an equilibrium streamfunction given by 
! pnep = -f*snep*snep*ht and is applied through friction whereby 
! the solution is damped towards the equilibrium streamfunction 
! rather than being damped towards zero kinetic energy. 
!
! ht    = depth of tracer cells 
! snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
!
! Neptune length scale snep has a value of senep at the
! equator and smoothly changes to spnep at the poles
!
! Reference:
! Holloway, G., 1992: Representing topographic stress for large
! scale ocean models, J. Phys. Oceanogr., 22, 1033-1046
!
! Eby and Holloway, 1994: Sensitivity of a large scale ocean model
! to a parameterization of topographic stress.  JPO, vol. 24,
! pages 2577-2588
!
! March 2012
! Stephen.Griffies 
! Algorithm updated to Eby and Holloway (1994)
!
! </DESCRIPTION>
!
subroutine compute_neptune_velocity(Time)

  type(ocean_time_type), intent(in) :: Time
  real, dimension(isd:ied,jsd:jed)  :: neptune_ft
  real, dimension(isd:ied,jsd:jed)  :: neptune_psi
  real                              :: neptune_length
  real                              :: active_cells
  integer                           :: i,j,k
  integer                           :: num_smooth

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapgen_friction_mod (compute_neptune_velocity): needs initialization')
  endif 

  allocate (neptune_velocity(isd:ied,jsd:jed,2))
  neptune_velocity = 0.0
  if(.not. neptune) return 

  ! neptune_psi is at T-cell point 
  do j=jsd,jed
    do i=isd,ied
      neptune_ft(i,j)  = 2.0*omega_earth*sin(Grd%phit(i,j))
      neptune_length   = neptune_length_pole  &
      + 0.5*(neptune_length_eq-neptune_length_pole)*(1.0 + cos(2.0*Grd%phit(i,j)))
      neptune_psi(i,j) = -neptune_ft(i,j)*neptune_length*neptune_length*Grd%ht(i,j)
    enddo
  enddo

  ! compute velocity 
  neptune_velocity(:,:,1) = -FDY_ET(FAX(neptune_psi(:,:)))
  neptune_velocity(:,:,2) =  FDX_NT(FAY(neptune_psi(:,:)))
  do j=jsc,jec
     do i=isc,iec
        neptune_velocity(i,j,1) = neptune_velocity(i,j,1)*Grd%umask(i,j,1)/(neptune_depth_min+Grd%hu(i,j))
        neptune_velocity(i,j,2) = neptune_velocity(i,j,2)*Grd%umask(i,j,1)/(neptune_depth_min+Grd%hu(i,j)) 
     enddo
  enddo
  call mpp_update_domains(neptune_velocity(:,:,1), neptune_velocity(:,:,2),&
                          Dom%domain2d,gridtype=BGRID_NE)


  ! optional smoothing
  if(neptune_smooth) then 
      do num_smooth=1,neptune_smooth_num

         wrk1_v2d(:,:,:) = 0.0
         k=1
         do j=jsc,jec
            do i=isc,iec

               if(Grd%tmask(i,j,k)==1.0) then 

                   active_cells = 4.0       +&
                        Grd%umask(i-1,j,k)  +&
                        Grd%umask(i+1,j,k)  +&
                        Grd%umask(i,j-1,k)  +&
                        Grd%umask(i,j+1,k)

                   if (active_cells > 4.0) then
                       wrk1_v2d(i,j,1) =                  &  
                         (4.0*neptune_velocity(i,j,1)    +&
                            neptune_velocity(i-1,j,1)    +&
                            neptune_velocity(i+1,j,1)    +&
                            neptune_velocity(i,j-1,1)    +&
                            neptune_velocity(i,j+1,1)) / active_cells
                       wrk1_v2d(i,j,2) =                  &  
                         (4.0*neptune_velocity(i,j,2)    +&
                            neptune_velocity(i-1,j,2)    +&
                            neptune_velocity(i+1,j,2)    +&
                            neptune_velocity(i,j-1,2)    +&
                            neptune_velocity(i,j+1,2)) / active_cells
                   else
                       wrk1_v2d(i,j,1) = neptune_velocity(i,j,1)
                       wrk1_v2d(i,j,2) = neptune_velocity(i,j,2)
                   endif

               endif

            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               neptune_velocity(i,j,1) = wrk1_v2d(i,j,1)*Grd%umask(i,j,k)
               neptune_velocity(i,j,2) = wrk1_v2d(i,j,2)*Grd%umask(i,j,k)
            enddo
         enddo

         call mpp_update_domains(neptune_velocity(:,:,1), neptune_velocity(:,:,2),&
                                 Dom%domain2d,gridtype=BGRID_NE)

      enddo  ! enddo for number of smoothing iterations 
  endif      ! endif for smooth



  ! diagnostics 
  id_neptune_lap_u = register_static_field('ocean_model', 'neptune_lap_u',               &
                   Grd%vel_axes_uv(1:2), 'Zonal velocity from neptune Laplacian scheme', &
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d_u(Time, Grd, id_neptune_lap_u, neptune_velocity(:,:,1))

  id_neptune_lap_v = register_static_field('ocean_model', 'neptune_lap_v',                    &
                   Grd%vel_axes_uv(1:2), 'Meridional velocity from neptune Laplacian scheme', &
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d_u(Time, Grd, id_neptune_lap_v, neptune_velocity(:,:,2))

  id_neptune_psi   = register_static_field('ocean_model', 'neptune_psi',             &
                     Grd%tracer_axes(1:2), 'Transport for neptune parameterization', &
                     'm^3/sec', missing_value=missing_value, range=(/-1.e12,1.e12/))
  call diagnose_2d_u(Time, Grd, id_neptune_psi, neptune_psi(:,:))

  id_neptune_ft   = register_static_field('ocean_model', 'neptune_ft',                                &
                    Grd%tracer_axes(1:2),'Coriolis parameter used for  neptune parameterization Psi', &
                    '1/sec', missing_value=missing_value, range=(/-1.e3,1.e3/))
  call diagnose_2d_u(Time, Grd, id_neptune_ft, neptune_ft(:,:))


end subroutine compute_neptune_velocity
! </SUBROUTINE> NAME="compute_neptune_velocity"

end module ocean_lapgen_friction_mod
      
      




