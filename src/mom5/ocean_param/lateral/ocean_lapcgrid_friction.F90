module ocean_lapcgrid_friction_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes the thickness weighted time tendency for  
! horizontal velocity arising from horizontal Laplacian friction. 
! Friction is formulated for the C-grid here. 
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
!
!
! Friction is formulated for the C-grid here. 
!
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
!<NAMELIST NAME="ocean_lapcgrid_friction_nml">
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
!
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
!  large magnitudes.  We do the same thing with bottom drag calculations
!  in ocean_bbc. Default side_drag_friction_max=1.0.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: pi, radius, epsln, radian
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: open_namelist_file, check_nml_error, write_version_number, close_file
use mpp_domains_mod,  only: mpp_update_domains, CGRID_NE, mpp_global_field
use mpp_domains_mod,  only: mpp_start_update_domains, mpp_complete_update_domains
use mpp_mod,          only: input_nml_file, mpp_sum, mpp_pe, mpp_error, mpp_max
use mpp_mod,          only: FATAL, NOTE, stdout, stdlog

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: BAY, BAX, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_operators_mod,  only: FAY, FAX, FDX_NT, FDY_ET, FDX_T, FDY_T 
use ocean_parameters_mod, only: missing_value, onehalf, omega_earth, rho0
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type
use ocean_types_mod,      only: ocean_domain_type, ocean_adv_vel_type 
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,      only: ocean_options_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d, diagnose_2d_u, diagnose_3d_u, diagnose_2d, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4
use ocean_workspace_mod,  only: wrk1_2d, wrk1_v2d  
use ocean_workspace_mod,  only: wrk1_v, wrk2_v, wrk3_v

implicit none

private

public ocean_lapcgrid_friction_init
public lapcgrid_friction
public lapcgrid_viscosity_check
public lapcgrid_reynolds_check

private compute_neptune_velocity

! length of forward time step  
real :: dtime

! for diagnostics 
integer :: id_aiso                 =-1
integer :: id_aaniso               =-1
integer :: id_aiso_smag            =-1
integer :: id_aaniso_smag          =-1
integer :: id_along_back           =-1
integer :: id_across_back          =-1
integer :: id_aiso_diverge         =-1
integer :: id_sin2theta            =-1
integer :: id_cos2theta            =-1
integer :: id_neptune_lap_u        =-1
integer :: id_neptune_lap_v        =-1
integer :: id_neptune_psi          =-1
integer :: id_along                =-1
integer :: id_across               =-1
integer :: id_visc_crit_lap        =-1
integer :: id_lap_fric_u           =-1
integer :: id_lap_fric_v           =-1
integer :: id_horz_lap_diss        =-1
integer :: id_viscosity_scaling    =-1
integer :: id_tmask_next_to_land   =-1
integer :: id_stress_xx_lap        =-1
integer :: id_stress_xy_lap        =-1
integer :: id_stress_yx_lap        =-1

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
real, dimension(:,:,:), allocatable :: tmask_next_to_land


#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS
real, dimension(isd:ied,jsd:jed,nk)  :: aiso_back      ! minimum allowable isotropic viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk)  :: aaniso_back    ! minimum allowable anisotropic viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk)  :: aiso_xx        ! isotropic viscosity (m^2/s) for stress_xx 
real, dimension(isd:ied,jsd:jed,nk)  :: aiso_xy        ! isotropic viscosity (m^2/s) for stress_xy 
real, dimension(isd:ied,jsd:jed,nk)  :: aaniso_xx      ! anisotropic viscosity (m^2/s) for stress_xx 
real, dimension(isd:ied,jsd:jed,nk)  :: aaniso_xy      ! anisotropic viscosity (m^2/s) for stress_xy 
real, dimension(isd:ied,jsd:jed,nk)  :: horz_ten       ! horizontal rate of deformation (1/s)
real, dimension(isd:ied,jsd:jed,nk)  :: horz_str       ! horizontal rate of strain (1/s)
real, dimension(isd:ied,jsd:jed,nk)  :: stress_xx      ! stress tensor component
real, dimension(isd:ied,jsd:jed,nk)  :: stress_xy      ! stress tensor component
real, dimension(isd:ied,jsd:jed,nk)  :: stress_yx      ! stress tensor component
real, dimension(isd:ied,jsd:jed,nk)  :: dissipate      ! diagnosed lateral dissipation from friction 
#else
real, dimension(:,:,:),  allocatable  :: aiso_back      ! minimum allowable isotropic viscosity (m^2/s)
real, dimension(:,:,:),  allocatable  :: aaniso_back    ! minimum allowable anisotropic viscosity (m^2/s)
real, dimension(:,:,:),  allocatable  :: aiso_xx        ! isotropic viscosity (m^2/s) for stress_xx 
real, dimension(:,:,:),  allocatable  :: aiso_xy        ! isotropic viscosity (m^2/s) for stress_xy 
real, dimension(:,:,:),  allocatable  :: aaniso_xx      ! anisotropic viscosity (m^2/s) for stress_xx 
real, dimension(:,:,:),  allocatable  :: aaniso_xy      ! anisotropic viscosity (m^2/s) for stress_xy 
real, dimension(:,:,:),  allocatable  :: horz_ten       ! horizontal rate of tension (1/s)
real, dimension(:,:,:),  allocatable  :: horz_str       ! horizontal rate of strain (1/s)
real, dimension(:,:,:),  allocatable  :: stress_xx      ! stress tensor component
real, dimension(:,:,:),  allocatable  :: stress_xy      ! stress tensor component
real, dimension(:,:,:),  allocatable  :: stress_yx      ! stress tensor component
real, dimension(:,:,:),  allocatable  :: dissipate      ! diagnosed lateral dissipation from friction 
#endif

real, dimension(:,:),   allocatable  :: fsmag_iso    ! (m^2) combination of terms for computing isotropic smag visc
real, dimension(:,:),   allocatable  :: fsmag_aniso  ! (m^2) combination of terms for computing anisotropic smag visc
real, dimension(:,:),   allocatable  :: aiso_diverge ! (m^2/s) viscosity for divergence damping
real, dimension(:,:),   allocatable  :: visc_crit    ! critical value of the viscosity (m^2/sec) for linear stability 
real, dimension(:,:),   allocatable  :: dxu_dyur          ! (dxu/dyu)
real, dimension(:,:),   allocatable  :: dyu_dxur          ! (dyu/dxu)
real, dimension(:,:),   allocatable  :: dxt_dytr          ! (dxt/dyt)
real, dimension(:,:),   allocatable  :: dyt_dxtr          ! (dyt/dxt)
real, dimension(:,:),   allocatable  :: grid_length       ! horizontal grid length (m)
real, dimension(:,:),   allocatable  :: viscosity_scaling ! dimensionless function of grid scale and Rossby radius
real, dimension(:,:,:), allocatable  :: neptune_velocity  ! barotropic velocity (m/s) from Neptune 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

character(len=256) :: version=&
     '=>Using: ocean_lapcgrid_friction.F90 ($Id: ocean_lapcgrid_friction.F90,v 20.0 2013/12/14 00:14:22 fms Exp $)'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc              = .FALSE.

namelist /ocean_lapcgrid_friction_nml/ use_this_module, debug_this_module,              &
    k_smag_iso, k_smag_aniso,                                                           &
    vel_micom_iso, vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom, &
    equatorial_zonal, equatorial_zonal_lat, equatorial_no_smag,                         &
    neptune, neptune_length_eq, neptune_length_pole, neptune_depth_min,                 &
    neptune_smooth, neptune_smooth_num,                                                 &
    restrict_polar_visc, restrict_polar_visc_lat, restrict_polar_visc_ratio,            &
    divergence_damp, divergence_damp_vel_micom,                                         &
    viscosity_scale_by_rossby, viscosity_scale_by_rossby_power,                         &
    use_side_drag_friction, side_drag_friction_scaling,                                 &
    side_drag_friction_uvmag_max, side_drag_friction_max 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_lapcgrid_friction_init">

! <DESCRIPTION>
! Initialize the horizontal Laplacian friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_lapcgrid_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                        obc, use_lapcgrid_friction, debug)
 
  type(ocean_grid_type),    target, intent(in)    :: Grid
  type(ocean_domain_type),  target, intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: d_time
  logical,                          intent(in)    :: obc
  logical,                          intent(inout) :: use_lapcgrid_friction
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real :: coeff_iso, coeff_aniso, dxdy

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapcgrid_friction_mod (ocean_lapcgrid_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif 

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_lapcgrid_friction_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_lapcgrid_friction_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_lapcgrid_friction_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_lapcgrid_friction_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_lapcgrid_friction_nml)  
  write (stdlogunit,ocean_lapcgrid_friction_nml)

  if(use_this_module) then 
      call mpp_error(NOTE, '==> NOTE: USING ocean_lapcgrid_friction_mod.')
      Ocean_options%horz_lap_friction = 'Used general horizontal Laplacian friction for Cgrid.'
      use_lapcgrid_friction = .true. 
  else 
      call mpp_error(NOTE, '==> NOTE: NOT using ocean_lapcgrid_friction_mod.')
      Ocean_options%horz_lap_friction = 'Did NOT use horizontal Laplacian friction for Cgrid.'
      use_lapcgrid_friction = .false. 
      return 
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_lapcgrid_friction with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  dtime = d_time
  write(stdoutunit,'(/a,f10.2)') &
  '==> Note from ocean_lapcgrid_friction_mod: using forward time step of (secs)', dtime 

  have_obc = obc
  if (have_obc) write(stdoutunit,'(/a,f10.2)') &
   '==> Note from ocean_lapcgrid_friction_mod: considering obc' 

  if(viscosity_scale_by_rossby) then 
    write(stdoutunit,'(1x,a)') '==> Note: Scaling the laplacian viscosity according to grid scale and Rossby radius.'
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

#ifndef MOM_STATIC_ARRAYS
  allocate (aiso_back(isd:ied,jsd:jed,nk))
  allocate (aaniso_back(isd:ied,jsd:jed,nk))
  allocate (aiso_xx(isd:ied,jsd:jed,nk))
  allocate (aiso_xy(isd:ied,jsd:jed,nk))
  allocate (aaniso_xx(isd:ied,jsd:jed,nk))
  allocate (aaniso_xy(isd:ied,jsd:jed,nk))
  allocate (horz_ten(isd:ied,jsd:jed,nk))
  allocate (horz_str(isd:ied,jsd:jed,nk))
  allocate (stress_xx(isd:ied,jsd:jed,nk))  
  allocate (stress_xy(isd:ied,jsd:jed,nk))  
  allocate (stress_yx(isd:ied,jsd:jed,nk))  
  allocate (dissipate(isd:ied,jsd:jed,nk))  
#endif
  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (aiso_diverge(isd:ied,jsd:jed))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (dxu_dyur(isd:ied,jsd:jed))
  allocate (dyu_dxur(isd:ied,jsd:jed))
  allocate (dxt_dytr(isd:ied,jsd:jed))
  allocate (dyt_dxtr(isd:ied,jsd:jed))
  allocate (grid_length(isd:ied,jsd:jed))
  allocate (viscosity_scaling(isd:ied,jsd:jed))

  grid_length(:,:)       = Grd%tmask(:,:,1)*2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(epsln+Grd%dxt(:,:)+Grd%dyt(:,:))
  viscosity_scaling(:,:) = Grd%tmask(:,:,1)

  aiso_back(:,:,:)   = 0.0
  aaniso_back(:,:,:) = 0.0
  aiso_xx(:,:,:)     = 0.0
  aiso_xy(:,:,:)     = 0.0
  aaniso_xx(:,:,:)   = 0.0
  aaniso_xy(:,:,:)   = 0.0
  aiso_diverge(:,:)  = 0.0
  horz_ten(:,:,:)    = 0.0
  horz_str(:,:,:)    = 0.0
  stress_xx(:,:,:)   = 0.0
  stress_xy(:,:,:)   = 0.0
  stress_yx(:,:,:)   = 0.0
  dissipate(:,:,:)   = 0.0
  
  ! dimensionless ratios used for friction operator 
  dxu_dyur(:,:) = Grd%dxu(:,:)*Grd%dyur(:,:)
  dyu_dxur(:,:) = Grd%dyu(:,:)*Grd%dxur(:,:)
  dxt_dytr(:,:) = Grd%dxt(:,:)*Grd%dytr(:,:)
  dyt_dxtr(:,:) = Grd%dyt(:,:)*Grd%dxtr(:,:)

  ! critical value of viscosity, above which 2-dim linear stability is not satisfied
  ! visc_crit taken from equation (18.17) of Griffies book.
  visc_crit(:,:) = 0.0
  do j=jsd,jed
    do i=isd,ied
      dxdy = 1.0/( 1.0/(Grd%dxt(i,j)*Grd%dxt(i,j)) + 1.0/(Grd%dyt(i,j)*Grd%dyt(i,j)) )
      visc_crit(i,j) = 0.5*dxdy/(dtime + epsln)
    enddo
  enddo
  
  ! may need a more conservative restriction in high latitudes when coupling to ice models 
  if(restrict_polar_visc) then 
     do j=jsd,jed
       do i=isd,ied
          if (abs(Grd%yt(i,j)) > restrict_polar_visc_lat) then    
            visc_crit(i,j) = restrict_polar_visc_ratio*visc_crit(i,j)
          endif
       enddo
     enddo
  endif

  id_visc_crit_lap = register_static_field ('ocean_model', 'visc_crit_lap', &
                     Grd%tracer_axes(1:2), 'critical viscosity',            &
                     'm^2/sec',missing_value=missing_value, range=(/0.0,1.e20/))
  call diagnose_2d(Time, Grd, id_visc_crit_lap, visc_crit(:,:))

  ! compute some smag fields and parameters 

  ! f_smag = (delta s * k_smag/pi)^2 where (delta s) is for velocity cell
  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = (k_smag_iso/pi)**2
  coeff_aniso      = (k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = coeff_iso  *((2.0*Grd%dxt(:,:)*Grd%dyt(:,:))/(epsln + Grd%dxt(:,:) + Grd%dyt(:,:)))**2
  fsmag_aniso(:,:) = coeff_aniso*((2.0*Grd%dxt(:,:)*Grd%dyt(:,:))/(epsln + Grd%dxt(:,:) + Grd%dyt(:,:)))**2

  ! remove smag within equatorial_zonal_lat band 
  if(equatorial_no_smag) then
    do j=jsd,jed
      do i=isd,ied
        if(abs(Grd%yt(i,j)) < equatorial_zonal_lat) then
          fsmag_iso(i,j)   = 0.0
          fsmag_aniso(i,j) = 0.0
        endif 
      enddo
    enddo 
  endif 

  ! compute time independent background viscosities 
  aiso_back(:,:,:)   = 0.0
  aaniso_back(:,:,:) = 0.0

  ! grid scale dependent "MICOM" background viscosities
  do j=jsd,jed
     do i=isd,ied
        aiso_back(i,j,:)   = vel_micom_iso  *(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)) &
                                            /(epsln + Grd%dxt(i,j) + Grd%dyt(i,j))
        aaniso_back(i,j,:) = vel_micom_aniso*(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)) &
                                            /(epsln + Grd%dxt(i,j) + Grd%dyt(i,j))
        if(abs(Grd%yt(i,j)) < eq_lat_micom) then
            aiso_back(i,j,:)   = eq_vel_micom_iso  *(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)) &
                                                   /(epsln + Grd%dxt(i,j) + Grd%dyt(i,j))
            aaniso_back(i,j,:) = eq_vel_micom_aniso*(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)) &
                                                   /(epsln + Grd%dxt(i,j) + Grd%dyt(i,j))
        endif
     enddo
  enddo

  ! for divergence damping 
  do j=jsd,jed
     do i=isd,ied
        aiso_diverge(i,j) = divergence_damp_vel_micom*(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)) &
                            /(epsln + Grd%dxt(i,j) + Grd%dyt(i,j))
     enddo
  enddo


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

  ! compute static neptune viscosity 
  call compute_neptune_velocity(Time)


  if(use_side_drag_friction) then 
     write(stdoutunit,'(a)') '==> Note from ocean_lapcgrid_friction: using side drag friction for cells next to land.'
     allocate (tmask_next_to_land(isd:ied,jsd:jed,nk))
     tmask_next_to_land(:,:,:) = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              tmask_next_to_land(i,j,k) = (1.0-Grd%tmask(i+1,j,k)) + (1.0-Grd%tmask(i-1,j,k)) &
                                         +(1.0-Grd%tmask(i,j+1,k)) + (1.0-Grd%tmask(i,j-1,k)) 
              tmask_next_to_land(i,j,k) = Grd%tmask(i,j,k)*min(1.0,tmask_next_to_land(i,j,k))
           enddo
        enddo
     enddo 
     id_tmask_next_to_land = register_static_field('ocean_model', 'tmask_next_to_land_lap',&
            Grd%tracer_axes(1:3),'T-cell mask for cells next to land for lapcgrid module', &
            'dimensionless', missing_value=missing_value, range=(/-10.0,10.0/))
     call diagnose_3d(Time, Grd, id_tmask_next_to_land, tmask_next_to_land(:,:,:))
  endif 


  ! send static fields to diag_manager 
  id_along_back   = register_static_field('ocean_model', 'along_lap_back',      &
                    Grd%tracer_axes(1:3),'T-cell background along-stream visc', &
                    'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))
  id_across_back = register_static_field('ocean_model', 'across_lap_back',      &
                   Grd%tracer_axes(1:3), 'T-cell background cross-stream visc', &
                   'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))

  if (id_along_back > 0) then 
     call diagnose_3d(Time, Grd, id_along_back, aiso_back(:,:,:)+0.5*aaniso_back(:,:,:))
  endif 
  if (id_across_back > 0) then 
     call diagnose_3d(Time, Grd, id_across_back, aiso_back(:,:,:)-0.5*aaniso_back(:,:,:))
  endif 

  id_aiso_diverge = register_static_field('ocean_model', 'aiso_lap_diverge', &
              Grd%tracer_axes(1:2), 'T-cell visc for divergence damping',    &
              'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_2d(Time, Grd, id_aiso_diverge, aiso_diverge(:,:))

  ! other viscosities for diagnostic output 

  id_aiso   = register_diag_field ('ocean_model', 'aiso_lap', Grd%tracer_axes(1:3), &
              Time%model_time, 'T-cell isotropic visc', 'm^2/sec',                  &
              missing_value=missing_value, range=(/-10.0,1.e10/),                   &
              standard_name='ocean_momentum_xy_laplacian_diffusivity')
  id_aaniso = register_diag_field ('ocean_model', 'aaniso_lap', Grd%tracer_axes(1:3), &
              Time%model_time, 'T-cell anisotropic visc', 'm^2/sec',                  &
              missing_value=missing_value, range=(/-10.0,1.e10/))
  id_aiso_smag = register_diag_field ('ocean_model', 'aiso_lap_smag', Grd%tracer_axes(1:3), &
                 Time%model_time, 'T-cell isotropic lap visc from smag', 'm^2/sec',         &
                 missing_value=missing_value, range=(/-10.0,1.e10/))
  id_aaniso_smag = register_diag_field ('ocean_model', 'aaniso_lap_smag', Grd%tracer_axes(1:3), &
                   Time%model_time, 'T-cell anisotropic lap visc from smag', 'm^2/sec',         &
                   missing_value=missing_value, range=(/-10.0,1.e10/))
  id_along  = register_diag_field ('ocean_model', 'along', Grd%tracer_axes(1:3), &
              Time%model_time, 'T-cell along-stream visc', 'm^2/sec',            &
              missing_value=missing_value, range=(/-10.0,1.e10/))
  id_across = register_diag_field ('ocean_model', 'across', Grd%tracer_axes(1:3), &
              Time%model_time, 'T-cell cross-stream visc', 'm^2/sec',             &
              missing_value=missing_value, range=(/-10.0,1.e10/))
  id_sin2theta = register_diag_field ('ocean_model', 'sin2theta_lap', Grd%tracer_axes(1:3),               &
              Time%model_time, 'sin2theta for orientation angle with laplacian friction', 'dimensionless',&
              missing_value=missing_value, range=(/-10.0,10.0/))
  id_cos2theta = register_diag_field ('ocean_model', 'cos2theta_lap', Grd%tracer_axes(1:3),               &
              Time%model_time, 'cos2theta for orientation angle with laplacian friction', 'dimensionless',&
              missing_value=missing_value, range=(/-10.0,10.0/))
 
  id_lap_fric_u = register_diag_field ('ocean_model', 'lap_fric_u', Grd%tracer_axes(1:3),              &
                  Time%model_time, 'Thick & rho wghtd horz lap frict on u-zonal','(kg/m^3)*(m^2/s^2)', &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_lap_fric_v = register_diag_field ('ocean_model', 'lap_fric_v', Grd%tracer_axes(1:3), &
                  Time%model_time, 'Thick & rho wghtd horz lap frict on v-merid', '(kg/m^3)*(m^2/s^2)', &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_horz_lap_diss = register_diag_field ('ocean_model', 'horz_lap_diss', Grd%tracer_axes(1:3),     &
                  Time%model_time, 'Energy dissipation from horizontal Laplacian friction', 'W/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_viscosity_scaling = register_diag_field ('ocean_model', 'viscosity_scaling', Grd%tracer_axes(1:2),  &
              Time%model_time, 'Grid/Rossby radius scaling for the laplacian viscosity', 'dimensionless',&
              missing_value=missing_value, range=(/-1.0,10.0/))
  id_stress_xx_lap = register_diag_field ('ocean_model', 'stress_xx_lap', Grd%tracer_axes(1:3),  &
                  Time%model_time, 'Stress tensor xx component from Laplacian friction', 'N/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_stress_xy_lap = register_diag_field ('ocean_model', 'stress_xy_lap', Grd%vel_axes_uv(1:3),  &
                  Time%model_time, 'Stress tensor xy component from Laplacian friction', 'N/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_stress_yx_lap = register_diag_field ('ocean_model', 'stress_yx_lap', Grd%vel_axes_uv(1:3),  &
                  Time%model_time, 'Stress tensor xy component from Laplacian friction', 'N/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))


end subroutine ocean_lapcgrid_friction_init
! </SUBROUTINE>  NAME="ocean_lapcgrid_friction_init"


!#######################################################################
! <SUBROUTINE NAME="lapcgrid_friction">
!
! <DESCRIPTION>
! This routine computes thickness weighted and density weighted 
! time tendency for horizontal velocity arising from horizontal 
! Laplacian friction.  
!
! The algorithm is derived from a functional approach that ensures
! kinetic energy is consistenty dissipated for all flow configurations. 
! The stencil is far simpler than the B-grid approach. In particular, 
! there are no triads here for the C-grid.  
! 
! Fundamental to the scheme are the rates of horizontal deformation 
! horizontal tension = DT = (dy)(u/dy)_x - (dx)(v/dx)_y             
! horizontal strain  = DS = (dx)(u/dx)_y + (dy)(v/dy)_x             
! Units of the tension and strain are sec^-1.
!
! </DESCRIPTION>
!
subroutine lapcgrid_friction(Time, Thickness, Velocity, lap_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step
  real, dimension(isd:,jsd:), intent(inout) :: lap_viscosity

  real :: deform, delta 
  real :: usqrd, vsqrd, umagr
  real :: uvmag, umag, vmag, velocity_gamma  
  real :: term1, term2 

  integer :: i, j, k, n
  integer :: taum1, tau
  integer :: stdoutunit
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
    '==>Error in ocean_lapcgrid_friction_mod (lapcgrid_friction): module not initialized')
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
     do j=jsd,jed
        do i=isd,ied
           viscosity_scaling(i,j) = Grd%tmask(i,j,1)
        enddo
     enddo
  endif 

  ! store velocity relative to neptune 
  wrk1_v(:,:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk1_v(i,j,k,1) = Velocity%u(i,j,k,1,taum1) - neptune_velocity(i,j,1)
           wrk1_v(i,j,k,2) = Velocity%u(i,j,k,2,taum1) - neptune_velocity(i,j,2)
        enddo
     enddo
  enddo

  ! compute the horizontal deformation rates
  ! horz_ten = (dy)(u/dy)_x - (dx)(v/dx)_y 
  !          = located at tracer points 
  ! horz_str = (dx)(u/dx)_y + (dy)(v/dy)_x
  !          = located at corner vorticity points 
  !          = zero at corner points via umask 
  horz_ten(:,:,:) = 0.0
  horz_str(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           horz_ten(i,j,k) = dyt_dxtr(i,j)*(wrk1_v(i,j,k,1)  *Grd%dyter(i,j)  -wrk1_v(i-1,j,k,1)*Grd%dyter(i-1,j)) &
                            -dxt_dytr(i,j)*(wrk1_v(i,j,k,2)  *Grd%dxtnr(i,j)  -wrk1_v(i,j-1,k,2)*Grd%dxtnr(i,j-1))
           horz_str(i,j,k) = dxu_dyur(i,j)*(wrk1_v(i,j+1,k,1)*Grd%dxter(i,j+1)-wrk1_v(i,j,k,1)  *Grd%dxter(i,j))   &
                            +dyu_dxur(i,j)*(wrk1_v(i+1,j,k,2)*Grd%dytnr(i+1,j)-wrk1_v(i,j,k,2)  *Grd%dytnr(i,j))
           horz_ten(i,j,k) = horz_ten(i,j,k)*Grd%tmask(i,j,k) 
           horz_str(i,j,k) = horz_str(i,j,k)*Grd%umask(i,j,k) 
        enddo
     enddo
  enddo
  call mpp_update_domains (horz_ten(:,:,:), Dom%domain2d)    
  call mpp_update_domains (horz_str(:,:,:), Dom%domain2d)    


  ! compute viscosities according to Smagorinsky scheme 
  wrk1 = 0.0   ! aiso_xx_smag 
  wrk2 = 0.0   ! aiso_xy_smag   
  wrk3 = 0.0   ! aaniso_xx_smag
  wrk4 = 0.0   ! aaniso_xy_smag
  if(k_smag_iso > 0.0) then 
     do k=1,nk 
        do j=jsc,jec
           do i=isc,iec
              deform      = sqrt(horz_ten(i,j,k)**2 + horz_str(i,j,k)**2)
              wrk1(i,j,k) = fsmag_iso(i,j)*deform
              wrk3(i,j,k) = fsmag_aniso(i,j)*deform

              deform      = sqrt(horz_ten(i,j,k)**2 + horz_str(i+1,j+1,k)**2)
              wrk2(i,j,k) = fsmag_iso(i,j)*deform
              wrk4(i,j,k) = fsmag_aniso(i,j)*deform
           enddo
        enddo
     enddo
  endif 
  ! write diagnostics here since use wrk arrays, some of which are reused later 
  if(.not. energy_analysis_step) then 
     if (id_aiso_smag > 0) then
        call diagnose_3d(Time, Grd, id_aiso_smag, onehalf*(wrk1(:,:,:)+wrk2(:,:,:)))
     endif 
     if (id_aaniso_smag > 0) then
        call diagnose_3d(Time, Grd, id_aaniso_smag, onehalf*(wrk3(:,:,:)+wrk4(:,:,:)))
     endif  
  endif 

  ! compute the full viscosities 
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec
           aiso_xx(i,j,k)   = viscosity_scaling(i,j)*min(visc_crit(i,j), aiso_back(i,j,k)   + wrk1(i,j,k))
           aaniso_xx(i,j,k) = viscosity_scaling(i,j)*min(visc_crit(i,j), aaniso_back(i,j,k) + wrk3(i,j,k))
           aiso_xy(i,j,k)   = viscosity_scaling(i,j)*min(visc_crit(i,j), aiso_back(i,j,k)   + wrk2(i,j,k))
           aaniso_xy(i,j,k) = viscosity_scaling(i,j)*min(visc_crit(i,j), aaniso_back(i,j,k) + wrk4(i,j,k))
        enddo
     enddo
  enddo

  ! compute orientation angle; ignore offset of u,v 
  wrk1 = 0.0    ! sin2theta
  wrk2 = 0.0    ! cos2theta
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec
           if(equatorial_zonal .and. abs(Grd%yt(i,j)) <= equatorial_zonal_lat) then
             wrk1(i,j,k) = 0.0
             wrk2(i,j,k) = 1.0
           else
             usqrd = wrk1_v(i,j,k,1)**2
             vsqrd = wrk1_v(i,j,k,2)**2
             umagr = 1.0/(epsln + usqrd + vsqrd)
             wrk1(i,j,k) = 2.0*wrk1_v(i,j,k,1)*wrk1_v(i,j,k,2)*umagr
             wrk2(i,j,k) = (usqrd-vsqrd)*umagr
           endif
        enddo
     enddo
  enddo
  if(.not. energy_analysis_step) then 
     call diagnose_3d(Time, Grd, id_sin2theta, wrk1(:,:,:))
     call diagnose_3d(Time, Grd, id_cos2theta, wrk2(:,:,:))
  endif 


  ! compute the stress tensor components (N/m2)
  ! stress_xx(i,j) = rho*
  !   [aiso_xx(i,j)*horz_ten(i,j) + aaniso_xx(i,j)*sin2theta*0.5*(horz_str(i,j)*cos2theta - horz_ten(i,j)*sin2theta)] 
  ! stress_xy(i,j) = 0.5*rho*
  !   [aiso_xy(i,j)*horz_str(i,j) - aaniso_xy(i,j)*cos2theta*0.5*(horz_str(i,j)*cos2theta - horz_ten(i+1,j+1)*sin2theta)]
  ! set rho = rho0 even for non-Boussinesq 
  ! stress_yy = -stress_xx
  ! stress_xy =  stress_yx except for possible differences next to boundaries with partial slip 
  stress_xx(:,:,:) = 0.0
  stress_xy(:,:,:) = 0.0
  stress_yx(:,:,:) = 0.0
  
  ! umask applied to stress_xy corresponds to free-slip side boundaries 
  wrk1 = 0.0
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec

           delta = onehalf*(horz_str(i,j,k)*wrk2(i,j,k)-horz_ten(i+1,j+1,k)*wrk1(i,j,k))
           stress_xy(i,j,k) = rho0*Grd%umask(i,j,k) &
                              *(aiso_xy(i,j,k)*horz_str(i,j,k) - aaniso_xy(i,j,k)*wrk2(i,j,k)*delta)
           stress_yx(i,j,k) = stress_xy(i,j,k) 

           delta = onehalf*(horz_str(i,j,k)*wrk2(i,j,k)-horz_ten(i,j,k)*wrk1(i,j,k))
           stress_xx(i,j,k) = rho0*grd%tmask(i,j,k) &
                              *(aiso_xx(i,j,k)*horz_ten(i,j,k) + aaniso_xx(i,j,k)*wrk1(i,j,k)*delta)

           ! estimate the dissipation (W/m2) based on diagonal terms 
           wrk1(i,j,k) = -Thickness%rho_dzt(i,j,k,tau) &
                        *(aiso_xx(i,j,k)*(horz_ten(i,j,k)**2 + horz_str(i,j,k)**2) - 2.0*aaniso_xx(i,j,k)*delta**2)

        enddo
     enddo
   enddo

   ! side drag law for partial slip on stress_xy and stress_yx (N/m^2)
   wrk2_v(:,:,:,:) = 0.0
   if(use_side_drag_friction) then
      do k=1,nk    
         do j=jsc,jec
            do i=isc,iec
               uvmag            = sqrt(wrk1_v(i,j,k,1)**2 + wrk1_v(i,j,k,2)**2)
               uvmag            = min(uvmag,side_drag_friction_uvmag_max)
               velocity_gamma   = rho0*side_drag_friction_scaling*Velocity%cdbot_array(i,j)*uvmag
               umag             = abs(wrk1_v(i,j,k,1))
               vmag             = abs(wrk1_v(i,j,k,2))
               wrk2_v(i,j,k,1)  = -min(side_drag_friction_max, velocity_gamma*umag)*sign(1.0,wrk1_v(i,j,k,1))
               wrk2_v(i,j,k,2)  = -min(side_drag_friction_max, velocity_gamma*vmag)*sign(1.0,wrk1_v(i,j,k,2))
               stress_xy(i,j,k) = stress_xy(i,j,k) + wrk2_v(i,j,k,1)*tmask_next_to_land(i,j,k)   
               stress_yx(i,j,k) = stress_yx(i,j,k) + wrk2_v(i,j,k,2)*tmask_next_to_land(i,j,k)   
            enddo
         enddo
      enddo
   endif 

   ! need stress components in halos 
   call mpp_update_domains (stress_xx(:,:,:), Dom%domain2d)    
   call mpp_update_domains (stress_xy(:,:,:), Dom%domain2d)    
   call mpp_update_domains (stress_yx(:,:,:), Dom%domain2d)    


   ! acceleration from lateral friction kg/(m*s^2) = N/m^2
   ! the factor of dzt ensures proper units, and scales down friction for thin cells.  
   do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            term1 = Grd%dyter(i,j)*(dyt_dxtr(i+1,j)*stress_xx(i+1,j,k)-dyt_dxtr(i,j)  *stress_xx(i,j,k))  
            term2 = Grd%dxter(i,j)*(dxu_dyur(i,j)  *stress_xy(i,j,k)  -dxu_dyur(i,j-1)*stress_xy(i,j-1,k))  
            wrk3_v(i,j,k,1) = Thickness%dzt(i,j,k)*(term1 + term2)

            term1 = -Grd%dxtnr(i,j)*(dxt_dytr(i,j+1)*stress_xx(i,j+1,k) - dxt_dytr(i,j)  *stress_xx(i,j,k))  
            term2 =  Grd%dytnr(i,j)*(dyu_dxur(i,j)  *stress_yx(i,j,k)   - dyu_dxur(i-1,j)*stress_yx(i-1,j,k))  
            wrk3_v(i,j,k,2) = Thickness%dzt(i,j,k)*(term1 + term2)
         enddo
      enddo
   enddo

   if(energy_analysis_step) then 
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%wrkv(i,j,k,n) = wrk3_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
      return 
   endif 

   do n=1,2
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk3_v(i,j,k,n)
            enddo
         enddo
      enddo
   enddo

   ! vertically averaged isotropic portion of the viscosity 
   lap_viscosity(:,:) = 0.0
   wrk1_2d(:,:)       = 0.0
   do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            lap_viscosity(i,j) = lap_viscosity(i,j) + Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau) &
                                                      *onehalf*(aiso_xx(i,j,k) + aiso_xy(i,j,k))
            wrk1_2d(i,j)       = wrk1_2d(i,j) + Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
         enddo
      enddo
   enddo
   do j=jsc,jec
      do i=isc,iec
         lap_viscosity(i,j) = Grd%tmask(i,j,1)*lap_viscosity(i,j)/(epsln+wrk1_2d(i,j))
      enddo
   enddo
   call mpp_update_domains (lap_viscosity(:,:), Dom%domain2d)    


   if (id_aiso > 0) then 
      wrk1(:,:,:) = onehalf*(aiso_xx(:,:,:)+aiso_xy(:,:,:))    
      call diagnose_3d(Time, Grd, id_aiso, wrk1(:,:,:))
   endif 
   if (id_aaniso > 0) then 
      wrk1(:,:,:) = onehalf*(aaniso_xx(:,:,:)+aaniso_xy(:,:,:))      
      call diagnose_3d(Time, Grd, id_aaniso, wrk1(:,:,:))
   endif 
   if (id_along > 0) then 
      wrk1(:,:,:) = onehalf*(aiso_xx(:,:,:)+aiso_xy(:,:,:) + onehalf*(aaniso_xx(:,:,:)+aaniso_xy(:,:,:)))
      call diagnose_3d(Time, Grd, id_along,  wrk1(:,:,:))
   endif 
   if (id_across > 0) then 
      wrk1(:,:,:) = onehalf*(aiso_xx(:,:,:)+aiso_xy(:,:,:) - onehalf*(aaniso_xx(:,:,:)+aaniso_xy(:,:,:)))
      call diagnose_3d(Time, Grd, id_across,  wrk1(:,:,:))
   endif 
   call diagnose_3d_u(Time, Grd, id_lap_fric_u, wrk3_v(:,:,:,1))
   call diagnose_3d_u(Time, Grd, id_lap_fric_v, wrk3_v(:,:,:,2))
   call diagnose_2d_u(Time, Grd, id_viscosity_scaling, viscosity_scaling(:,:))
   call diagnose_3d(Time, Grd, id_horz_lap_diss, dissipate(:,:,:))
   call diagnose_3d(Time, Grd, id_stress_xx_lap, stress_xx(:,:,:))
   call diagnose_3d_u(Time, Grd, id_stress_xy_lap, stress_xy(:,:,:))
   call diagnose_3d_u(TIme, Grd, id_stress_yx_lap, stress_yx(:,:,:))

   if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_lapcgrid_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('lapcgrid friction(1)', wrk1_v(COMP,:,1)*Grd%umask(COMP,:))
      call write_chksum_3d('lapcgrid friction(2)', wrk1_v(COMP,:,2)*Grd%umask(COMP,:))
   endif


end subroutine lapcgrid_friction
! </SUBROUTINE> NAME="lapcgrid_friction"



!#######################################################################
! <SUBROUTINE NAME="lapcgrid_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the Laplacian 
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine lapcgrid_viscosity_check

  integer       :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: fudge

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL,&
     '==>Error in ocean_lapcgrid_friction_mod (lapcgrid_viscosity_check): needs initialization')
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
        if (Grd%tmask(i,j,k) > 0.0) then 
          if (aiso_xx(i,j,k) > fudge*visc_crit(i,j) .and. num < max_num) then
            num = num + 1
            write (unit,9600) &
            'aiso(',aiso_xx(i,j,k), visc_crit(i,j), i, j, k, &
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
9601  format(/' Warning: ',a,es16.8,' m^2/s) exceeds max value (',es16.8,') at (i,j,k) = ', &
            '(',i4,',',i4,',',i4,'),',' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')


end subroutine lapcgrid_viscosity_check
! </SUBROUTINE> NAME="lapcgrid_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="lapcgrid_reynolds_check">
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
subroutine lapcgrid_reynolds_check(Time, Velocity)

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
    '==>Error in ocean_lapcgrid_friction_mod (lapcgrid_reynolds_check): needs initialization')
  endif 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(aiso_xx(i,j,k) + epsln)
        ramn = 1.0/(aiso_xx(i,j,k) + epsln)
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


end subroutine lapcgrid_reynolds_check
! </SUBROUTINE> NAME="lapcgrid_reynolds_check"


!#######################################################################
! <SUBROUTINE NAME="compute_neptune_velocity">
!
! <DESCRIPTION>
! Compute Neptune velocity for c-grid.  
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
! hu    = depth of B-grid velocity corner 
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
! May 2012
! Stephen.Griffies 
! upgraded to Cgrid 
!
! </DESCRIPTION>
!
subroutine compute_neptune_velocity(Time)

  type(ocean_time_type), intent(in) :: Time
  real, dimension(isd:ied,jsd:jed)  :: neptune_psi
  real                              :: neptune_length
  real                              :: active_cells_u
  real                              :: active_cells_v
  integer                           :: i,j,k
  integer                           :: num_smooth

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapcgrid_friction_mod (compute_neptune_velocity): needs initialization')
  endif 

  allocate (neptune_velocity(isd:ied,jsd:jed,2))
  neptune_velocity = 0.0
  if(.not. neptune) return 

  ! neptune_psi defined at U-cell point 
  do j=jsd,jed
     do i=isd,ied
        neptune_length   = neptune_length_pole  &
        + onehalf*(neptune_length_eq-neptune_length_pole)*(1.0 + cos(2.0*Grd%phiu(i,j)))
        neptune_psi(i,j) = -Grd%f(i,j)*neptune_length*neptune_length*Grd%hu(i,j)
     enddo
  enddo

  ! compute neptune velocity on C-grid 
  do j=jsc,jec
    do i=isc,iec
       neptune_velocity(i,j,1) = -(neptune_psi(i,j)-neptune_psi(i,j-1))*Grd%dyunr(i,j-1)
       neptune_velocity(i,j,2) =  (neptune_psi(i,j)-neptune_psi(i-1,j))*Grd%dxuer(i-1,j)
    enddo
  enddo
  call mpp_update_domains(neptune_velocity(:,:,1), neptune_velocity(:,:,2),&
                          Dom%domain2d,gridtype=CGRID_NE)


  ! optional smoothing
  if(neptune_smooth) then 
      do num_smooth=1,neptune_smooth_num

         wrk1_v2d(:,:,:) = 0.0
         k=1
         do j=jsc,jec
            do i=isc,iec

               if(Grd%tmask(i,j,k)==1.0) then 

                   active_cells_u = 4.0         +&
                        Grd%tmasken(i-1,j,k,1)  +&
                        Grd%tmasken(i+1,j,k,1)  +&
                        Grd%tmasken(i,j-1,k,1)  +&
                        Grd%tmasken(i,j+1,k,1)

                   active_cells_v = 4.0         +&
                        Grd%tmasken(i-1,j,k,2)  +&
                        Grd%tmasken(i+1,j,k,2)  +&
                        Grd%tmasken(i,j-1,k,2)  +&
                        Grd%tmasken(i,j+1,k,2)

                   if (active_cells_u > 4.0) then
                       wrk1_v2d(i,j,1) =                  &  
                         (4.0*neptune_velocity(i,j,1)    +&
                              neptune_velocity(i-1,j,1)  +&
                              neptune_velocity(i+1,j,1)  +&
                              neptune_velocity(i,j-1,1)  +&
                              neptune_velocity(i,j+1,1)) / active_cells_u
                   else
                       wrk1_v2d(i,j,1) = neptune_velocity(i,j,1)
                   endif 

                   if (active_cells_v > 4.0) then
                       wrk1_v2d(i,j,2) =                    &  
                         (4.0*neptune_velocity(i,j,2)      +&
                              neptune_velocity(i-1,j,2)    +&
                              neptune_velocity(i+1,j,2)    +&
                              neptune_velocity(i,j-1,2)    +&
                              neptune_velocity(i,j+1,2)) / active_cells_v
                   else
                       wrk1_v2d(i,j,2) = neptune_velocity(i,j,2)
                   endif

               endif

            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               neptune_velocity(i,j,1) = wrk1_v2d(i,j,1)*Grd%tmasken(i,j,k,1)
               neptune_velocity(i,j,2) = wrk1_v2d(i,j,2)*Grd%tmasken(i,j,k,2)
            enddo
         enddo

         call mpp_update_domains(neptune_velocity(:,:,1), neptune_velocity(:,:,2),&
                                 Dom%domain2d,gridtype=CGRID_NE)

      enddo  ! enddo for number of smoothing iterations 
  endif      ! endif for smooth



  ! diagnostics; output onto T-cell is fine  
  id_neptune_lap_u = register_static_field('ocean_model', 'neptune_lap_u',                &
                     Grd%tracer_axes(1:2), 'Zonal velocity from neptune Laplacian scheme',&
                     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d_u(Time, Grd, id_neptune_lap_u, neptune_velocity(:,:,1))

  id_neptune_lap_v = register_static_field('ocean_model', 'neptune_lap_v',                   &
                   Grd%tracer_axes(1:2), 'Meridional velocity from neptune Laplacian scheme',&
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d_u(Time, Grd, id_neptune_lap_v, neptune_velocity(:,:,2))

  id_neptune_psi   = register_static_field('ocean_model', 'neptune_psi',            &
                     Grd%vel_axes_uv(1:2), 'Transport for neptune parameterization',&
                     'm^3/sec', missing_value=missing_value, range=(/-1.e12,1.e12/))
  call diagnose_2d_u(Time, Grd, id_neptune_psi, neptune_psi(:,:))


end subroutine compute_neptune_velocity
! </SUBROUTINE> NAME="compute_neptune_velocity"

end module ocean_lapcgrid_friction_mod
      
      




