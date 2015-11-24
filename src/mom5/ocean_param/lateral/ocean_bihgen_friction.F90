module ocean_bihgen_friction_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes thickness weighted and density weighted
! time tendency for horizontal velocity arising from 
! biharmonic friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes thickness weighted and density weighted
! time tendency for horizontal velocity arising from biharmonic
! friction. 
!
! The viscosity used to determine the strength of the tendency 
! can be a general function of space and time as specified by 
! the Smagorinsky approach as well as a grid-scale dependent
! background viscosity.  The form of the friction operator 
! can be isotropic or anisotropic in the lateral plane. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies and R.W. Hallberg, 2000:
! Biharmonic friction with a Smagorinsky viscosity for use
! in large-scale eddy-permitting ocean models
! Monthly Weather Review, vol. 128, pages 2935-2946.
! </REFERENCE>
!
! <REFERENCE>
! R.D. Smith and J.C. McWilliams, 2003:
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
! S.M. Griffies, 2004: 
! Fundamentals of Ocean Climate Models 
! Princeton University Press
! </REFERENCE>
!
! <REFERENCE>
!  Griffies, Elements of MOM (2012)
! </REFERENCE>
!
! <NOTE>
! The ocean model can generally run with both Laplacian and biharmonic
! friction enabled at the same time.  Such has been found useful 
! for some eddying ocean simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bihgen_friction_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging by printing checksums.  
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
!  <DATA NAME="vel_micom_iso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_aniso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity. 
!  </DATA> 
!
!  <DATA NAME="visc_crit_scale" UNITS="dimensionless" TYPE="real">
!  Scaling factor used to determine the critical viscosity, above which 
!  the viscosity is not allowed to reach. 
!  Use visc_crit_scale < 1.0 for cases where the visc_crit from linear stability 
!  allows for still too large of a viscosity.  Use visc_crit_scale>1.0 when wish 
!  to allow for larger viscosity. Default is visc_crit_scale=1.0.
!  </DATA> 
!
!  <DATA NAME="equatorial_zonal" TYPE="real">
!  Orient the anisotropic friction within a latitudinal band according 
!  to zonal direction. 
!  </DATA> 
!  <DATA NAME="equatorial_zonal_lat" TYPE="real">
!  Latitudinal band to use the zonal friction orientation. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_iso" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity 
!  within a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_aniso" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic 
!  viscosity within a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_lat_micom" TYPE="real">
!  Equatorial latitude band (degrees) within which the MICOM viscosity
!  is set according to eq_vel_micom_iso and eq_vel_micom_aniso.
!  </DATA> 
!  <DATA NAME="bottom_5point" TYPE="logical">
!  To alleviate problems with small partial cells, it is often necessary
!  to reduce the operator to the traditional 5-point Laplacian at the 
!  ocean bottom.  This logical implements this mixing. 
!  Default bottom_5point=.false. 
!  </DATA> 
!  <DATA NAME="vel_micom_bottom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity for 
!  5point Laplacian at the bottom. 
!  </DATA> 
!
!  <DATA NAME="ncar_boundary_scaling" TYPE="logical">
!  To enhance the velocity scale used in western boundaries 
!  for the isotropic and anisotropic  background viscosities, 
!  we compute a scaling using the algorithm from the laplacian
!  NCAR anisotropic scheme.
!  Default ncar_boundary_scaling=.false.
!  </DATA> 
!  <DATA NAME="ncar_boundary_scaling_read" TYPE="logical">
!  To read in the ncar boundary scaling field rather than 
!  generating it during initialization.  Generating during 
!  initialization can be a bottle-neck on fine resolution models
!  since there are some global 2d fields needed.  So if the 
!  rescaling is produced once and then saved, it can be read
!  in during subsequent runs without incurring the slowdown 
!  of re-generating the scalings.  
!  Default ncar_boundary_scaling_read=.false.
!  </DATA> 
!  <DATA NAME="ncar_rescale_power" UNITS="dimensionless" TYPE="integer">
!  For determining rescaling of the viscosity so to enhance the 
!  friction near the western boundaries. Default ncar_rescale_power=1.
!  </DATA> 
!  <DATA NAME="ncar_vconst_4" UNITS="1/cm" TYPE="real">
!  Inverse damping length for exponential falloff of the velocity scale 
!  as move eastward away from western boundary. Default ncar_vconst_4=2.e-8.
!  </DATA> 
!  <DATA NAME="ncar_vconst_5" UNITS="dimensionless" TYPE="integer">
!  For determining number of grid points in boundary calculation.
!  Default ncar_vconst_5=3.
!  </DATA> 
!
!  <DATA NAME="neptune" TYPE="logical">
!  Set to true for computing friction relative to Neptune barotropic velocity. 
!  Default neptune=.false. 
!  </DATA> 
!  <DATA NAME="neptune_length_eq" UNITS="m" TYPE="real">
!  Length scale used to compute Neptune velocity at equator.  
!  Default neptune_length_eq= 4.2e3 from Maltrud and Holloway.
!  </DATA> 
!  <DATA NAME="neptune_length_pole" UNITS="m" TYPE="real">
!  Length scale used to compute Neptune velocity at pole. 
!  Default neptune_length_pole= 17.0e3 from Maltrud and Holloway.
!  </DATA> 
!  <DATA NAME="neptune_depth_min" UNITS="m" TYPE="real">
!  Minimum depth scale used for computing Neptune velocity.
!  Default neptune_depth_min=100.0 from Maltrud and Holloway.
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
!  <DATA NAME="neptune_scaling" TYPE="real">
!  Overall scaling parameter to help tune neptune.  
!  Default neptune_scaling=1.0
!  </DATA> 
!
!  <DATA NAME="visc_diverge_scaling" TYPE="real">
!  Dimensionless scaling used for divergence based viscosity.  
!  Default visc_diverge_scaling=0.0 turns off the scheme. 
!  visc_diverge_scaling=10.0 produces sensible viscosities for 
!  1-degree model.  May need tuning for different resolutions.  
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
!
use constants_mod,       only: pi, radius, epsln
use diag_manager_mod,    only: register_diag_field, register_static_field
use fms_mod,             only: open_namelist_file, check_nml_error
use fms_mod,             only: write_version_number, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog, read_data
use mpp_domains_mod,     only: mpp_update_domains, mpp_global_field, XUPDATE 
use mpp_domains_mod,     only: mpp_global_min, mpp_global_max, BGRID_NE
use mpp_mod,             only: input_nml_file, mpp_error, mpp_max, mpp_sum, mpp_pe

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: BAY, BAX, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_parameters_mod, only: missing_value, oneeigth, omega_earth, rho0, rho0r
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type, ocean_adv_vel_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type, ocean_options_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d_u, diagnose_2d_u, diagnose_2d, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_v, wrk2_v, wrk3_v, wrk1_v2d

implicit none

private

public ocean_bihgen_friction_init
public bihgen_friction
public bihgen_viscosity_check
public bihgen_reynolds_check

private ncar_boundary_scale_create 
private ncar_boundary_scale_read
private BDX_EU_smag
private BDY_NU_smag
private compute_neptune_velocity

real :: k_smag_iso          = 0.0    ! smag scaling coeff for isotropic viscosity (dimensionless) 
real :: k_smag_aniso        = 0.0    ! smag scaling coeff for anisotripic viscosity (dimensionless) 
real :: vel_micom_iso       = 0.0    ! background scaling velocity for isotropic viscosity (m/sec) 
real :: vel_micom_aniso     = 0.0    ! background scaling velocity for anisotropic viscosity (m/sec) 
real :: eq_vel_micom_iso    = 0.2    ! background scaling velocity (m/sec) within equatorial band for isotropic visc
real :: eq_vel_micom_aniso  = 0.0    ! background scaling velocity (m/sec) within equatorial band for anisotropic visc
real :: eq_lat_micom        = 0.0    ! equatorial latitude band for micom (degrees)
real :: vel_micom_bottom    = 0.01   ! velocity scale for determining viscosity at the bottom 
real :: equatorial_zonal_lat= 0.0    ! latitudinal band for orienting the friction along zonal direction
real :: visc_crit_scale     = 1.0    ! constant for setting the critical viscosity on the instability constraint 
logical :: equatorial_zonal = .false.! zonally orient anisotropic friction w/i equatorial_zonal_lat
logical :: bottom_5point    =.false.  ! for bottom Laplacian 5point mixing to avoid problems with thin partial cells. 

! for the western boundary scaling of the background viscosity
! nx and ny are number of global points in generalized x and y 
! directions
logical :: ncar_boundary_scaling      = .false.
logical :: ncar_boundary_scaling_read = .false.
real    :: ncar_vconst_4              = 2.e-8
integer :: ncar_vconst_5              = 3  
integer :: ncar_rescale_power         = 1
integer :: nx, ny 

! for neptune scheme  
logical :: neptune              = .false. ! for computing friction relative to barotropic Neptune velocity
logical :: neptune_smooth       = .true.  ! for smoothing diagnosed neptune velocity
integer :: neptune_smooth_num   = 1       ! number of smoothing passes for neptune smoother 
real    :: neptune_length_eq    = 4.2e3   ! (metres; from Maltrud and Holloway)
real    :: neptune_length_pole  = 17.0e3  ! (metres; from Maltrud and Holloway)
real    :: neptune_depth_min    = 100.0   ! (metres; from Maltrud and Holloway)
real    :: neptune_scaling      = 1.0     ! dimensionless scaling used for tuning 

! for divergence-based viscosity
real :: visc_diverge_scaling = 0.0  ! dimensionless
real, dimension(:,:),   allocatable :: visc_diverge_factor ! (m^4) for computing visc_diverge
real, dimension(:,:,:), allocatable :: divergence_t        ! array to hold diverge_t 

! for side-boundary drag law
logical :: use_side_drag_friction       = .false.
real    :: side_drag_friction_scaling   = 1.0
real    :: side_drag_friction_max       = 1.0
real    :: side_drag_friction_uvmag_max = 10.0
real, dimension(:,:,:), allocatable :: umask_next_to_land

! for diagnostics 
logical :: used
integer :: id_aiso                 =-1
integer :: id_aaniso               =-1
integer :: id_along                =-1
integer :: id_across               =-1
integer :: id_aiso_back            =-1
integer :: id_aaniso_back          =-1
integer :: id_along_back           =-1
integer :: id_across_back          =-1
integer :: id_visc_diverge         =-1
integer :: id_visc_crit_bih        =-1
integer :: id_bih_fric_u           =-1
integer :: id_bih_fric_v           =-1
integer :: id_horz_bih_diss        =-1
integer :: id_ncar_rescale         =-1
integer :: id_neptune_bih_u        =-1
integer :: id_neptune_bih_v        =-1
integer :: id_neptune_fu           =-1
integer :: id_umask_next_to_land   =-1
integer :: id_bih_plus_side_fric_u =-1
integer :: id_bih_plus_side_fric_v =-1
integer :: id_side_drag_friction_u =-1
integer :: id_side_drag_friction_v =-1

! time step 
real ::  dtime = 0.0  ! time step used for friction tendency (2*dtuv if threelevel, dtuv if twolevel) 

#include <ocean_memory.h>


#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed)    :: fsmag_iso       ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(isd:ied,jsd:jed)    :: fsmag_aniso     ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(isd:ied,jsd:jed,nk) :: aiso_back       ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk) :: aaniso_back     ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk) :: aiso            ! isotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: visc_diverge    ! viscosity (m^4/s) arising from |grad(diverge_t)|
real, dimension(isd:ied,jsd:jed,nk) :: ncar_rescale    ! viscosity scale (cm^2/s) based on NCAR western boundary scaling 

real, dimension(isd:ied,jsd:jed)    :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed)    :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed,nk) :: massqc      ! (1/4)*dxu*dyu*rho_dzu (kg) time-independent quarter-cell mass
real, dimension(isd:ied,jsd:jed)    :: visc_crit   ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(isd:ied,jsd:jed)    :: visc_bottom ! Laplacian viscosity at the bottom (m^2/sec)

real, dimension(isd:ied,jsd:jed,nk) :: stress001   !stress tensor: last index stress_xx
real, dimension(isd:ied,jsd:jed,nk) :: stress011   !stress tensor: last index stress_xx
real, dimension(isd:ied,jsd:jed,nk) :: stress101   !stress tensor: last index stress_xx
real, dimension(isd:ied,jsd:jed,nk) :: stress111   !stress tensor: last index stress_xx
real, dimension(isd:ied,jsd:jed,nk) :: stress002   !stress tensor: last index stress_xy
real, dimension(isd:ied,jsd:jed,nk) :: stress012   !stress tensor: last index stress_xy
real, dimension(isd:ied,jsd:jed,nk) :: stress102   !stress tensor: last index stress_xy
real, dimension(isd:ied,jsd:jed,nk) :: stress112   !stress tensor: last index stress_xy
real, dimension(isd:ied,jsd:jed,nk) :: tmplap1
real, dimension(isd:ied,jsd:jed,nk) :: tmplap2
real, dimension(isd:ied,jsd:jed)    :: tmpfdelx1
real, dimension(isd:ied,jsd:jed)    :: tmpfdelx2
real, dimension(isd:ied,jsd:jed)    :: tmpfdely1
real, dimension(isd:ied,jsd:jed)    :: tmpfdely2
real, dimension(isd:ied,jsd:jed)    :: tmp1
real, dimension(isd:ied,jsd:jed)    :: tmp2
real, dimension(isc:iec,jsc:jec,nk) :: cos2theta
real, dimension(isc:iec,jsc:jec,nk) :: sin2theta

#else

real, dimension(:,:), allocatable   :: fsmag_iso    ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(:,:), allocatable   :: fsmag_aniso  ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(:,:,:), allocatable :: aiso_back    ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(:,:,:), allocatable :: aaniso_back  ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(:,:,:), allocatable :: aiso         ! isotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(:,:,:), allocatable :: visc_diverge ! viscosity (m^4/s) arising from |grad(diverge_t)|
real, dimension(:,:,:), allocatable :: ncar_rescale ! viscosity scale (cm^2/s) based on NCAR western boundary scaling 

real, dimension(:,:), allocatable   :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(:,:), allocatable   :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(:,:,:), allocatable :: massqc      ! (1/4)*dxu*dyu*rho_dzu (kg) time-independent quarter-cell mass
real, dimension(:,:), allocatable   :: visc_crit   ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(:,:), allocatable   :: visc_bottom ! Laplacian viscosity at the bottom (m^2/sec)

real, dimension(:,:,:), allocatable :: stress001  !stress tensor: last index stress_xx
real, dimension(:,:,:), allocatable :: stress011  !stress tensor: last index stress_xx
real, dimension(:,:,:), allocatable :: stress101  !stress tensor: last index stress_xx
real, dimension(:,:,:), allocatable :: stress111  !stress tensor: last index stress_xx
real, dimension(:,:,:), allocatable :: stress002  !stress tensor: last index stress_xy
real, dimension(:,:,:), allocatable :: stress012  !stress tensor: last index stress_xy
real, dimension(:,:,:), allocatable :: stress102  !stress tensor: last index stress_xy
real, dimension(:,:,:), allocatable :: stress112  !stress tensor: last index stress_xy
real, dimension(:,:,:), allocatable :: tmplap1
real, dimension(:,:,:), allocatable :: tmplap2
real, dimension(:,:), allocatable   :: tmpfdelx1
real, dimension(:,:), allocatable   :: tmpfdelx2
real, dimension(:,:), allocatable   :: tmpfdely1
real, dimension(:,:), allocatable   :: tmpfdely2
real, dimension(:,:), allocatable   :: tmp1
real, dimension(:,:), allocatable   :: tmp2
real, dimension(:,:,:), allocatable   :: cos2theta
real, dimension(:,:,:), allocatable   :: sin2theta

#endif

real, dimension(:,:), allocatable :: aiso_back_obc   ! for enhancing visc next to OBC boundaries 
real, dimension(:,:), allocatable :: aaniso_back_obc ! for enhancing visc next to OBC boundaries 

! barotropic velocity (m/s) from Neptune 
real, dimension(:,:,:), allocatable :: neptune_velocity 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

character(len=256) :: version=&
     '=>Using: ocean_bihgen_friction.F90 ($Id: ocean_bihgen_friction.F90,v 20.0 2013/12/14 00:14:16 fms Exp $)'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc              = .FALSE.
logical :: read_aiso_bih_back    = .false.

namelist /ocean_bihgen_friction_nml/ use_this_module, debug_this_module,           &
          k_smag_iso, k_smag_aniso, vel_micom_iso,                                 &
          vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom,     &
          vel_micom_bottom, bottom_5point, equatorial_zonal, equatorial_zonal_lat, &
          visc_crit_scale, read_aiso_bih_back,                                     &
          ncar_boundary_scaling, ncar_rescale_power, ncar_vconst_4, ncar_vconst_5, &
          ncar_boundary_scaling_read,                                              &
          neptune, neptune_length_eq, neptune_length_pole, neptune_depth_min,      &
          neptune_scaling, neptune_smooth, neptune_smooth_num,                     &
          visc_diverge_scaling,                                                    &
          use_side_drag_friction, side_drag_friction_scaling,                      &
          side_drag_friction_uvmag_max, side_drag_friction_max 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bihgen_friction_init">

! <DESCRIPTION>
! Initialize the horizontal biharmonic friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_bihgen_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                      obc, use_bihgen_friction, debug)

  type(ocean_grid_type),    target, intent(in)    :: Grid
  type(ocean_domain_type),  target, intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: d_time
  logical,                          intent(in)    :: obc
  logical,                          intent(inout) :: use_bihgen_friction
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real    :: coeff_iso, coeff_aniso, dxdy

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_bihgen_friction_mod(ocean_bihgen_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_bihgen_friction_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bihgen_friction_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_bihgen_friction_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_bihgen_friction_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_bihgen_friction_nml)  
  write (stdlogunit,ocean_bihgen_friction_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  if(use_this_module) then 
      call mpp_error(NOTE, '==> NOTE: USING ocean_bihgen_friction_mod.')
      Ocean_options%horz_bih_friction = 'Used general horizontal biharmonic friction.'
      use_bihgen_friction=.true.
  else 
      call mpp_error(NOTE, '==> NOTE: NOT using ocean_bihgen_friction_mod.')
      Ocean_options%horz_bih_friction = 'Did NOT use horizontal biharmonic friction.'
      use_bihgen_friction=.false.
      return 
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  Grd => Grid
  Dom => Domain
  nx  =  Grd%ni
  ny  =  Grd%nj

  dtime = d_time

  write(stdoutunit,'(a)') ' ' 

  write(stdoutunit,'(/a,f10.2)') &
   '==> Note from ocean_bihgen_friction_mod: using forward time step of (secs)', dtime 

  have_obc = obc
  if (have_obc) then 
     write(stdoutunit,'(/a,f10.2)') '==> Note from ocean_bihgen_friction_mod: considering obc' 
  endif 

  if (Dom%xhalo > 1 .or. Dom%yhalo > 1) then
    write (stdoutunit,'(/1x,a)') ' ==> NOTE: General biharmonic friction allows xhalo=yhalo=1.'
  endif

  if(bottom_5point) then
    write(stdoutunit,'(/)')  
    write(stdoutunit,'(/1x,a)') ' ==> NOTE: Will make horizontal friction to a 5point Laplacian on the bottom'
    write(stdoutunit,'(a)')     '     This helps alleviate numerical problems with thin bottom partial cells.'
  endif 

  if(neptune) then 
    write(stdoutunit,'(1x,a)') '==> Note: Computing biharmonic friction relative to Neptune'
    write(stdoutunit,'(5x,a,e10.5)') &
    'Equatorial length (m) for computing Neptune velocity is   ',neptune_length_eq 
    write(stdoutunit,'(5x,a,e10.5)') &
    'Polar length (m) for computing Neptune velocity is        ',neptune_length_pole 
    write(stdoutunit,'(5x,a,e10.5)') &
    'Minimum depth scale (m) for computing Neptune velocity is ',neptune_depth_min
  endif 

  if(k_smag_iso > 0.0) then 
     write( stdoutunit,'(a)')' Computing horziontal isotropic biharmonic viscosity via Smagorinsky.'
  endif
  if(k_smag_iso==0.0)  then 
    write( stdoutunit,*)' Setting horzizontal isotropic biharmonic Smagorinsky viscosity to zero.'
  endif 
  if(k_smag_aniso > 0.0) then 
     write( stdoutunit,'(a)')' Computing horzizontal anisotropic biharmonic viscosity via Smagorinsky.'
     write( stdoutunit,'(1x,a)')' ==>WARNING: anisotropic biharmonic Smagorinsky not well tested.'
  endif
  if(k_smag_aniso==0.0) then 
    write( stdoutunit,*)' Setting horzizontal anisotropic biharmonic Smagorinsky viscosity to zero.'
  endif 
  if(vel_micom_iso > 0.0) then 
     write( stdoutunit,'(a)')' Computing background horzizontal biharmonic isotropic viscosity via MICOM.'
  endif
  if(vel_micom_iso==0.0) then 
    write( stdoutunit,*)' Setting background horzizontal biharmonic isotropic viscosity to zero.'
  endif 
  if(vel_micom_aniso > 0.0) then
     write( stdoutunit,'(a)')' Computing background horzizontal anisotropic biharmonic viscosity via MICOM.'
     write( stdoutunit,'(1x,a)')' ==> WARNING: anisotropic biharmonic background viscosity not well tested.'
  endif 
  if(vel_micom_aniso==0.0) then
     write(stdoutunit,'(a)')' Setting background horziontal biharmonic anisotropic viscosity to zero.' 
  endif 
  if(k_smag_aniso > 2.0*k_smag_iso) then
     call mpp_error(FATAL, &
     '==>Error: Smag horizontal bih anisotropic visc too large. Must be < 2*isotropic-visc.')
  endif  
  if(vel_micom_aniso > 2.0*vel_micom_iso) then
    call mpp_error(FATAL, &
     '==>Error: Background bih horz aniso visc too large. Must be < 2.0*background iso visc.')
  endif  
  if(eq_lat_micom > 0) then 
    write( stdoutunit,'(a)')'Setting different background horz bih viscosity w/i equatorial zone.'
  endif 
  if(equatorial_zonal) then
    write( stdoutunit,'(a)') &
    'If using anisotropic biharmonic friction, zonally orient friction w/i latitudinal'
    write( stdoutunit,'(a)') &
    'band of width ',equatorial_zonal_lat,' degrees.'
  endif 


#ifndef MOM_STATIC_ARRAYS
  
  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (ncar_rescale(isd:ied,jsd:jed,nk))
  allocate (aiso_back(isd:ied,jsd:jed,nk))
  allocate (aaniso_back(isd:ied,jsd:jed,nk))
  allocate (aiso(isd:ied,jsd:jed,nk)) 
  allocate (visc_diverge(isd:ied,jsd:jed,nk)) 
  allocate (daur_dxur(isd:ied,jsd:jed))
  allocate (daur_dyur(isd:ied,jsd:jed))
  allocate (massqc(isd:ied,jsd:jed,nk))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (visc_bottom(isd:ied,jsd:jed))
  allocate (stress001(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xx
  allocate (stress011(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xx
  allocate (stress101(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xx
  allocate (stress111(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xx
  allocate (stress002(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xy
  allocate (stress012(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xy
  allocate (stress102(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xy
  allocate (stress112(isd:ied,jsd:jed,nk))  !stress tensor: last index stress_xy
  allocate (tmplap1(isd:ied,jsd:jed,nk))    
  allocate (tmplap2(isd:ied,jsd:jed,nk))    
  allocate (tmpfdelx1(isd:ied,jsd:jed))    
  allocate (tmpfdelx2(isd:ied,jsd:jed))    
  allocate (tmpfdely1(isd:ied,jsd:jed))   
  allocate (tmpfdely2(isd:ied,jsd:jed))   
  allocate (tmp1(isd:ied,jsd:jed))
  allocate (tmp2(isd:ied,jsd:jed))
  allocate (cos2theta(isc:iec,jsc:jec,nk))
  allocate (sin2theta(isc:iec,jsc:jec,nk))
  
#endif

  allocate (aiso_back_obc(isd:ied,jsd:jed))
  allocate (aaniso_back_obc(isd:ied,jsd:jed))
  allocate (visc_diverge_factor(isd:ied,jsd:jed))
  allocate (divergence_t(isd:ied,jsd:jed,nk))
 
  tmp1(:,:) = 0.0
  tmp2(:,:) = 0.0

  ! Some commonly used grid arrays 
  daur_dxur(:,:) = 0.0
  daur_dyur(:,:) = 0.0
  daur_dxur(:,:) = Grd%daur(:,:)*Grd%dxur(:,:)
  daur_dyur(:,:) = Grd%daur(:,:)*Grd%dyur(:,:)

  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = 0.125*(k_smag_iso/pi)**2
  coeff_aniso      = 0.125*(k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = Grd%umask(:,:,1)*coeff_iso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4
  fsmag_aniso(:,:) = Grd%umask(:,:,1)*coeff_aniso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4

  aiso_back_obc(:,:)      = 0.0
  aaniso_back_obc(:,:)    = 0.0
  visc_bottom(:,:)        = 0.0
  do j=jsd,jed
    do i=isd,ied 
      visc_bottom(i,j)   = Grd%umask(i,j,1)*vel_micom_bottom* &
                           ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))
      aiso_back_obc(i,j)     = Grd%umask(i,j,1)*vel_micom_iso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      aaniso_back_obc(i,j)   = Grd%umask(i,j,1)*vel_micom_aniso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      if(abs(Grd%yu(i,j)) < eq_lat_micom) then
        aiso_back_obc(i,j)   = Grd%umask(i,j,1)*eq_vel_micom_iso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
        aaniso_back_obc(i,j) = Grd%umask(i,j,1)*eq_vel_micom_aniso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      endif
    enddo  
  enddo

  aiso(:,:,:)         = 0.0
  aiso_back(:,:,:)    = 0.0
  aaniso_back(:,:,:)  = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied 

           aiso_back(i,j,k)     = Grd%umask(i,j,k)*vel_micom_iso* &
                ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           aaniso_back(i,j,k)   = Grd%umask(i,j,k)*vel_micom_aniso* &
                ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3

           if(abs(Grd%yu(i,j)) < eq_lat_micom) then
               aiso_back(i,j,k)   = Grd%umask(i,j,k)*eq_vel_micom_iso* &
                    ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
               aaniso_back(i,j,k) = Grd%umask(i,j,k)*eq_vel_micom_aniso* &
                    ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           endif

        enddo
     enddo
  enddo

  ! factor used for computing viscosity based on divergence 
  visc_diverge(:,:,:) = 0.0
  if(visc_diverge_scaling > 0.0) then 
     write(stdoutunit,'(/a,f10.2)') &
     '==> Note from ocean_bihgen_friction_mod: using visc_diverge_scaling for viscosity, with a value = ',visc_diverge_scaling 
  endif 
  do j=jsc,jec
     do i=isc,iec  
        visc_diverge_factor(i,j) = Grd%umask(i,j,1)*visc_diverge_scaling* &
                                   ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**4
      enddo
  enddo
  call mpp_update_domains (visc_diverge_factor(:,:),Dom%domain2d)    


  ! rescale the background viscosities according to NCAR boundary scaling  
  if(ncar_boundary_scaling) then 
      write(stdoutunit,'(a)') ' Note: rescaling background bih viscosities so they are larger in western boundaries.'
      if(ncar_boundary_scaling_read) then 
          call ncar_boundary_scale_read(Time) 
      else 
          call ncar_boundary_scale_create(Time)
      endif
  endif

  if (read_aiso_bih_back) then
      call read_data('INPUT/aiso_bih_back.nc','aiso_bih_back',aiso_back, Dom%domain2d)
  endif  

  ! enhance background viscosity near open boundaries if needed
  if (have_obc) then 
      call ocean_obc_enhance_visc_back(aiso_back_obc, 'bih')
      call ocean_obc_enhance_visc_back(aaniso_back_obc, 'bih')
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if(aiso_back_obc(i,j)   > aiso_back(i,j,k))   aiso_back(i,j,k)   = aiso_back_obc(i,j)
               if(aaniso_back_obc(i,j) > aaniso_back(i,j,k)) aaniso_back(i,j,k) = aaniso_back_obc(i,j)
            enddo
         enddo
      enddo
  endif

  call compute_neptune_velocity(Time)

  ! critical value of viscosity, above which 2-dim linear stability is not satisfied.
  ! visc_crit taken from equation (18.26) of Griffies book. 
  visc_crit(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
      dxdy = 1.0/( 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j)) )
      visc_crit(i,j) =  visc_crit_scale*oneeigth*dxdy**2/(dtime+epsln)
    enddo
  enddo

  if(use_side_drag_friction) then 
     write(stdoutunit,'(a)') '==> Note from ocean_bihgen_friction: using side drag friction for cells next to land.'
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
     id_umask_next_to_land = register_static_field('ocean_model', 'umask_next_to_land_bih',&
            Grd%vel_axes_uv(1:3),'U-cell mask for cells next to land for bihgen module',   &
            'dimensionless', missing_value=missing_value, range=(/-10.0,10.0/))
     call diagnose_3d_u(Time, Grd, id_umask_next_to_land, umask_next_to_land(:,:,:))
  endif 

  id_visc_crit_bih = register_static_field ('ocean_model', 'visc_crit_bih', &
                     Grd%vel_axes_uv(1:2), 'critical viscosity', 'm^4/sec', &
                     missing_value=missing_value, range=(/0.0,1.e20/))
  call diagnose_2d_u(Time, Grd, id_visc_crit_bih, visc_crit(:,:))


  ! ensure that background viscosities are not too large 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(aiso_back(i,j,k)   > visc_crit(i,j))  aiso_back(i,j,k)   = visc_crit(i,j)
           if(aaniso_back(i,j,k) > visc_crit(i,j))  aaniso_back(i,j,k) = visc_crit(i,j)
        enddo
     enddo
  enddo

  ! ensure viscosities maintain proper relative values so that friction dissipates 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(aaniso_back(i,j,k) >= 2.0*aiso_back(i,j,k) .and. aaniso_back(i,j,k) > 0.0)  then 
               write(stdoutunit,'(a,i4,a,i4,a,i3,a)')'Violating bih iso/aniso constraint at (',i,',',j,')'
               aaniso_back(i,j,k) = 1.9*aiso_back(i,j,k)
           endif
        enddo
     enddo
  enddo

  ! apply U-cell masking to background viscosities; needed in particular 
  ! for the case when read_aiso_back is enabled, to remove any spurious 
  ! missing values.  
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           aiso_back(i,j,k)   = aiso_back(i,j,k)*Grd%umask(i,j,k)
           aaniso_back(i,j,k) = aaniso_back(i,j,k)*Grd%umask(i,j,k)
        enddo
     enddo
  enddo


  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


  ! diagnostic output 
  id_aiso    = register_diag_field ('ocean_model', 'aiso_bih', Grd%vel_axes_uv(1:3), &
               Time%model_time,'U-cell isotropic bih visc', 'm^4/sec',               &
               missing_value=-10.0, range=(/-10.0,1.e20/),                           &
               standard_name='ocean_momentum_xy_biharmonic_diffusivity')
  id_aaniso  = register_diag_field ('ocean_model', 'aaniso_bih', Grd%vel_axes_uv(1:3), &
               Time%model_time,'U-cell  anisotropic bih visc', 'm^4/sec',              &
               missing_value=missing_value, range=(/-10.0,1.e20/))
  id_along  = register_diag_field ('ocean_model', 'along_bih', Grd%vel_axes_uv(1:3),   &
              Time%model_time,'U-cell along-stream bih visc', 'm^4/sec',               &
              missing_value=missing_value, range=(/-10.0,1.e20/))
  id_across  = register_diag_field ('ocean_model', 'across_bih', Grd%vel_axes_uv(1:3), &
               Time%model_time, 'U-cell cross-stream bih visc', 'm^4/sec',             &
               missing_value=missing_value, range=(/-10.0,1.e20/))
  id_visc_diverge = register_diag_field ('ocean_model', 'visc_diverge_bih', Grd%vel_axes_uv(1:3),    &
               Time%model_time,'U-cell biharmonic visc based on horz_diverge_t gradients', 'm^4/sec',&
               missing_value=-10.0, range=(/-10.0,1.e20/))

  id_bih_fric_u = register_diag_field ('ocean_model', 'bih_fric_u', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict on u-zonal',   &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_bih_fric_v = register_diag_field ('ocean_model', 'bih_fric_v', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict on v-merid',   &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_aiso_back   = register_static_field('ocean_model', 'aiso_bih_back', Grd%vel_axes_uv(1:3),   &
                    'U-cell background aiso bih visc', 'm^4/sec',                                &
                    missing_value=missing_value, range=(/-10.0,1.e20/))
  id_aaniso_back = register_static_field('ocean_model', 'aaniso_bih_back', Grd%vel_axes_uv(1:3), &
                    'U-cell background aaniso bih visc', 'm^4/sec',                              &
                    missing_value=missing_value, range=(/-10.0,1.e20/))
  id_along_back   = register_static_field('ocean_model', 'along_bih_back', Grd%vel_axes_uv(1:3), &
                    'U-cell background along-stream bih visc', 'm^4/sec',                        &
                    missing_value=missing_value, range=(/-10.0,1.e20/))
  id_across_back = register_static_field('ocean_model', 'across_bih_back', Grd%vel_axes_uv(1:3), &
                   'U-cell background cross-stream bih visc', 'm^4/sec',                         &
                   missing_value=missing_value, range=(/-10.0,1.e20/))

  call diagnose_3d_u(Time, Grd, id_aiso_back, aiso_back(:,:,:))
  call diagnose_3d_u(Time, Grd, id_aaniso_back, aaniso_back(:,:,:))
  if (id_along_back > 0) then 
     call diagnose_3d_u(Time, Grd, id_along_back, aiso_back(:,:,:)+0.5*aaniso_back(:,:,:))
  endif 
  if (id_across_back > 0) then 
     call diagnose_3d_u(Time, Grd, id_across_back, aiso_back(:,:,:)-0.5*aaniso_back(:,:,:))
  endif 

  id_horz_bih_diss = register_diag_field ('ocean_model', 'horz_bih_diss', Grd%vel_axes_uv(1:3),      &
                  Time%model_time, 'Energy dissipation from horizontal biharmonic friction', 'W/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_side_drag_friction_u = register_diag_field ('ocean_model', 'side_drag_friction_bih_u', Grd%vel_axes_uv(1:3),&
                  Time%model_time, 'u-friction due to side drag from bihgen module', 'N/m^2',                    &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_side_drag_friction_v = register_diag_field ('ocean_model', 'side_drag_friction_bih_v', Grd%vel_axes_uv(1:3),&
                  Time%model_time, 'v-friction due to side drag from bihgen module', 'N/m^2',                    &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_bih_plus_side_fric_u = register_diag_field ('ocean_model', 'bih_plus_side_fric_u', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict + side friction on u-zonal',       &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_bih_plus_side_fric_v = register_diag_field ('ocean_model', 'bih_plus_side_fric_v', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict + side friction on v-merid',       &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))

end subroutine ocean_bihgen_friction_init
! </SUBROUTINE>  NAME="ocean_bihgen_friction_init"


!#######################################################################
! <SUBROUTINE NAME="bihgen_friction">
!
! <DESCRIPTION>
! This subroutine computes the time tendency for horizontal 
! velocity (i.e., the acceleration) from horizontal biharmonic friction.  
! The algorithm is derived from a functional approach that ensures kinetic 
! energy is consistenty dissipated for all flow configurations. 
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
! As shown in Griffies and Hallberg (2000), 
! a biharmonic operator with a nonconstant viscosity is guaranteed to 
! dissipate kinetic energy *only* when using the sqrt of the biharmonic
! viscosity at each of the two stages of the algorithm. 
! The sqrt approach is employed here.  
!
! </DESCRIPTION>
!
subroutine bihgen_friction(Time, Thickness, Adv_vel, Velocity, bih_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  real, dimension(isd:,jsd:), intent(inout) :: bih_viscosity
  logical,                    intent(in)    :: energy_analysis_step

  real, dimension(0:1,0:1,isc:iec,jsc:jec,nk) :: aiso_smag
  real, dimension(0:1,0:1,isc:iec,jsc:jec,nk) :: aaniso_smag
  real    :: uvmag, umag, vmag, velocity_gamma
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
  real :: tension, strain, deform, delta 
  real :: tension_metric, strain_metric
  real :: grad_diverge_sq 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihgen_friction_mod (horz_bih_friction): needs to be initialized')
  endif 

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

  stress001(:,:,:) = 0.0
  stress011(:,:,:) = 0.0
  stress101(:,:,:) = 0.0
  stress111(:,:,:) = 0.0
  stress002(:,:,:) = 0.0
  stress012(:,:,:) = 0.0
  stress102(:,:,:) = 0.0
  stress112(:,:,:) = 0.0
  tmplap1(:,:,:)   = 0.0
  tmplap2(:,:,:)   = 0.0
  wrk1_v(:,:,:,:)  = 0.0
  wrk2_v(:,:,:,:)  = 0.0
  wrk3_v(:,:,:,:)  = 0.0
  wrk1(:,:,:)      = 0.0
  wrk2(:,:,:)      = 0.0

  taum1 = Time%taum1
  tau   = Time%tau
  
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           massqc(i,j,k) = 0.25*Grd%dau(i,j)*Thickness%rho_dzu(i,j,k,tau)
        enddo
     enddo
  enddo

  ! divergence_t has units m/s
  divergence_t(:,:,:) = 0.0
  if(visc_diverge_scaling > 0.0) then 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              divergence_t(i,j,k) = rho0r*Adv_vel%diverge_t(i,j,k)
           enddo
        enddo
     enddo
     call mpp_update_domains (divergence_t(:,:,:), Dom%domain2d)    
  endif 

  ! big k-loop 
  do k=1,nk

     ! Laplacian part of algorithm
     do j=jsc,jec
        do i=isc,iec

           ! compute visc_diverge using gradient that straddles U-cell northeast sides  
           grad_diverge_sq     = ((divergence_t(i+1,j,k)-divergence_t(i,j,k))*Grd%dxter(i,j))**2 &
                                +((divergence_t(i,j+1,k)-divergence_t(i,j,k))*Grd%dytnr(i,j))**2 
           visc_diverge(i,j,k) = Grd%umask(i,j,k)*visc_diverge_factor(i,j)*sqrt(grad_diverge_sq)

           tension_metric = -(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))*Grd%dh2dx(i,j) &
                            +(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*Grd%dh1dy(i,j)  

           strain_metric  = -(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))*Grd%dh1dy(i,j) &
                            -(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*Grd%dh2dx(i,j)  

           if(equatorial_zonal .and. abs(Grd%yu(i,j)) <= equatorial_zonal_lat) then  
             sin2theta(i,j,k) = 0.0
             cos2theta(i,j,k) = 1.0
           else 
             usqrd = (Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1))**2
             vsqrd = (Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))**2
             umagr = 1.0/(epsln + usqrd + vsqrd)
             sin2theta(i,j,k) = 2.0*(Velocity%u(i,j,k,1,taum1)-neptune_velocity(i,j,1)) &
                                 *(Velocity%u(i,j,k,2,taum1)-neptune_velocity(i,j,2))*umagr
             cos2theta(i,j,k) = (usqrd-vsqrd)*umagr
           endif  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are faster.
          ! dimensions[aiso_smag]=m^2/s^0.5
          ! dimensions[tension and strain]=1/s
          ! dimensions[stress]=m^2/s^1.5

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  
           u1_m10 = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)  
           u2_m10 = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)  

           ip=1  ; jq=0  
           u1_10 = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)  
           u2_10 = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2)  

           ip=0  ; jq=0  
           u1_00 = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)  
           u2_00 = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2) 

           ip=0  ; jq=1  
           u1_01  = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)   
           u2_01  = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2) 

           ip=0  ; jq=-1 
           u1_0m1 = Velocity%u(i+ip,j+jq,k,1,taum1) - neptune_velocity(i+ip,j+jq,1)    
           u2_0m1 = Velocity%u(i+ip,j+jq,k,2,taum1) - neptune_velocity(i+ip,j+jq,2) 

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq,i,j,k)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform   + visc_diverge(i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform + visc_diverge(i,j,k))
           aiso_smag(ip,jq,i,j,k)    = sqrt(aiso_smag(ip,jq,i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = sqrt(aaniso_smag(ip,jq,i,j,k))
           stress001(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress002(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq,i,j,k)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform   + visc_diverge(i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform + visc_diverge(i,j,k))
           aiso_smag(ip,jq,i,j,k)    = sqrt(aiso_smag(ip,jq,i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = sqrt(aaniso_smag(ip,jq,i,j,k))
           stress101(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress102(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq,i,j,k)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform   + visc_diverge(i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform + visc_diverge(i,j,k))
           aiso_smag(ip,jq,i,j,k)    = sqrt(aiso_smag(ip,jq,i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = sqrt(aaniso_smag(ip,jq,i,j,k))
           stress011(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress012(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(ip,jq,i,j,k)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform   + visc_diverge(i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform + visc_diverge(i,j,k))
           aiso_smag(ip,jq,i,j,k)    = sqrt(aiso_smag(ip,jq,i,j,k))
           aaniso_smag(ip,jq,i,j,k)  = sqrt(aaniso_smag(ip,jq,i,j,k))
           stress111(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress112(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

        enddo
     enddo
  enddo  !end of k-loop

  call mpp_update_domains (stress001(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress011(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress101(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress111(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress002(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress012(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress102(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress112(:,:,:), Dom%domain2d, complete=.true. )    

  do k=1,nk
     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains
     ! dimensions[tmpfdelx and tmpfdely]=kg*m^2/s^1.5

     tmpfdelx1(:,:) = 0.0
     tmpfdelx2(:,:) = 0.0
     do j=jsc,jec
        do i=isd,iec
           tmpfdelx1(i,j) = (stress101(i,j,k)   + stress111(i,j,k)  )*massqc(i,j,k)    &
                +(stress001(i+1,j,k) + stress011(i+1,j,k))*massqc(i+1,j,k)
           tmpfdelx2(i,j) = (stress102(i,j,k)   + stress112(i,j,k)  )*massqc(i,j,k)    &  
                +(stress002(i+1,j,k) + stress012(i+1,j,k))*massqc(i+1,j,k)              
        enddo
     enddo

     tmpfdely1(:,:) = 0.0
     tmpfdely2(:,:) = 0.0
     do j=jsd,jec
        do i=isc,iec  
           tmpfdely1(i,j) = (stress012(i,j,k)   + stress112(i,j,k)  )*massqc(i,j,k)    &                      
                +(stress002(i,j+1,k) + stress102(i,j+1,k))*massqc(i,j+1,k) 
           tmpfdely2(i,j) = (stress011(i,j,k)   + stress111(i,j,k)  )*massqc(i,j,k)    &                      
                +(stress001(i,j+1,k) + stress101(i,j+1,k))*massqc(i,j+1,k)
        enddo
     enddo

     ! compute laplacian operator
     ! dimensions [tmplap=m/s^1.5]
     tmp1(:,:) = BDX_EU_smag(tmpfdelx1(:,:))+BDY_NU_smag(tmpfdely1(:,:))    
     tmp2(:,:) = BDX_EU_smag(tmpfdelx2(:,:))-BDY_NU_smag(tmpfdely2(:,:))    

     ! do not use rho_dzur since that is at taup1 and need to divide by rho_dau(tau)
     do j=jsc,jec
        do i=isc,iec
           tmplap1(i,j,k) = tmp1(i,j)*Grd%umask(i,j,k)/Thickness%rho_dzu(i,j,k,tau)
           tmplap2(i,j,k) = tmp2(i,j)*Grd%umask(i,j,k)/Thickness%rho_dzu(i,j,k,tau)
        enddo
     enddo
  enddo  !end of k-loop

  call mpp_update_domains (tmplap1(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (tmplap2(:,:,:), Dom%domain2d, complete=.true. )    


  do k=1,nk
     ! Second part of the iteration
     ! tmplap[m/s^1.5] replaces velocity[m/s]
     ! dimensions[aiso_smag]=m^2/s^0.5
     ! dimensions[tension and strain]=1/s^1.5
     ! dimensions[stress]=m^2/s^2
     do j=jsc,jec
        do i=isc,iec

           tension_metric = -tmplap1(i,j,k)*Grd%dh2dx(i,j) + tmplap2(i,j,k)*Grd%dh1dy(i,j)  
           strain_metric  = -tmplap1(i,j,k)*Grd%dh1dy(i,j) - tmplap2(i,j,k)*Grd%dh2dx(i,j)  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are quicker. 

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; u1_m10 = tmplap1(i+ip,j+jq,k) ; u2_m10 = tmplap2(i+ip,j+jq,k)
           ip=1  ; jq=0  ; u1_10  = tmplap1(i+ip,j+jq,k) ; u2_10  = tmplap2(i+ip,j+jq,k) 
           ip=0  ; jq=0  ; u1_00  = tmplap1(i+ip,j+jq,k) ; u2_00  = tmplap2(i+ip,j+jq,k) 
           ip=0  ; jq=1  ; u1_01  = tmplap1(i+ip,j+jq,k) ; u2_01  = tmplap2(i+ip,j+jq,k)
           ip=0  ; jq=-1 ; u1_0m1 = tmplap1(i+ip,j+jq,k) ; u2_0m1 = tmplap2(i+ip,j+jq,k)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           stress001(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress002(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           stress101(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress102(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           stress011(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress012(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j,k) - tension*sin2theta(i,j,k))
           stress111(i,j,k) = aiso_smag(ip,jq,i,j,k)*tension + aaniso_smag(ip,jq,i,j,k)*delta*sin2theta(i,j,k)
           stress112(i,j,k) = aiso_smag(ip,jq,i,j,k)*strain  - aaniso_smag(ip,jq,i,j,k)*delta*cos2theta(i,j,k)

        enddo
     enddo
  enddo  !end of k-loop

  call mpp_update_domains (stress001(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress011(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress101(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress111(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress002(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress012(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress102(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (stress112(:,:,:), Dom%domain2d, complete=.true. )    

  do k=1,nk
     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains
     ! dimensions[tmpfdelx and tmpfdely]=kg*m^2/s^2

     tmpfdelx1(:,:) = 0.0
     tmpfdelx2(:,:) = 0.0
     do j=jsc,jec
        do i=isd,iec
           tmpfdelx1(i,j) = (stress101(i,j,k)   + stress111(i,j,k)  )*massqc(i,j,k)    &
                +(stress001(i+1,j,k) + stress011(i+1,j,k))*massqc(i+1,j,k)
           tmpfdelx2(i,j) = (stress102(i,j,k)   + stress112(i,j,k)  )*massqc(i,j,k)    &  
                +(stress002(i+1,j,k) + stress012(i+1,j,k))*massqc(i+1,j,k)              
        enddo
     enddo

     tmpfdely1(:,:) = 0.0
     tmpfdely2(:,:) = 0.0
     do j=jsd,jec
        do i=isc,iec  
           tmpfdely1(i,j) = (stress012(i,j,k)   + stress112(i,j,k)  )*massqc(i,j,k)    &                      
                +(stress002(i,j+1,k) + stress102(i,j+1,k))*massqc(i,j+1,k) 
           tmpfdely2(i,j) = (stress011(i,j,k)   + stress111(i,j,k)  )*massqc(i,j,k)    &                      
                +(stress001(i,j+1,k) + stress101(i,j+1,k))*massqc(i,j+1,k)
        enddo
     enddo

    ! compute acceleration from horizontal biharmonic friction
    ! dimensions[tmp]=kg*/(m*s^2)
    ! dimensions[friction]=(kg/m^3)*(m^2/s^2)=N/m^2
     tmp1(:,:) = BDX_EU_smag(tmpfdelx1(:,:))+BDY_NU_smag(tmpfdely1(:,:)) 
     tmp2(:,:) = BDX_EU_smag(tmpfdelx2(:,:))-BDY_NU_smag(tmpfdely2(:,:)) 
     do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,1) = -tmp1(i,j)*Grd%umask(i,j,k)
           wrk1_v(i,j,k,2) = -tmp2(i,j)*Grd%umask(i,j,k)
        enddo
     enddo

     ! reduce to 5-point laplacian at bottom to avoid problems with thin bottom partial cells
     if (bottom_5point) then
         tmp1(:,:) = BDX_EU(visc_bottom(:,:)*FMX(Thickness%rho_dzu(:,:,k,tau)) &
                    *FDX_U(Velocity%u(:,:,k,1,taum1)-neptune_velocity(:,:,1))) &
                    +BDY_NU(visc_bottom(:,:)*FMY(Thickness%rho_dzu(:,:,k,tau)) &
                    *FDY_U(Velocity%u(:,:,k,1,taum1)-neptune_velocity(:,:,1)))

         tmp2(:,:) = BDX_EU(visc_bottom(:,:)*FMX(Thickness%rho_dzu(:,:,k,tau)) &
                    *FDX_U(Velocity%u(:,:,k,2,taum1)-neptune_velocity(:,:,2))) &
                    +BDY_NU(visc_bottom(:,:)*FMY(Thickness%rho_dzu(:,:,k,tau)) &
                    *FDY_U(Velocity%u(:,:,k,2,taum1)-neptune_velocity(:,:,2)))

         do j=jsc,jec
            do i=isc,iec
               if(k==Grd%kmu(i,j)) then
                  wrk1_v(i,j,k,1) = tmp1(i,j)*Grd%umask(i,j,k)
                  wrk1_v(i,j,k,2) = tmp2(i,j)*Grd%umask(i,j,k)
               endif
            enddo
         enddo
     endif

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

  else ! not energy_analysis_step

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      ! vertically averaged viscosity at U-cell centre
      bih_viscosity(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               aiso(i,j,k) = 0.25*(  aiso_smag(0,0,i,j,k)**2   + aiso_smag(1,0,i,j,k)**2    &
                                   + aiso_smag(0,1,i,j,k)**2   + aiso_smag(1,1,i,j,k)**2)
               wrk1(i,j,k) = 0.25*(  aaniso_smag(0,0,i,j,k)**2 + aaniso_smag(1,0,i,j,k)**2  &
                                   + aaniso_smag(0,1,i,j,k)**2 + aaniso_smag(1,1,i,j,k)**2)
               bih_viscosity(i,j) = bih_viscosity(i,j)  &
               + Grd%umask(i,j,k)*Thickness%rho_dzu(i,j,k,tau)*aiso(i,j,k)
            enddo
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            bih_viscosity(i,j) = Grd%umask(i,j,1)*bih_viscosity(i,j)/(epsln+Thickness%mass_u(i,j,tau))
         enddo
      enddo
      call mpp_update_domains (bih_viscosity(:,:), Dom%domain2d)    

      call diagnose_3d_u(Time, Grd, id_aiso, aiso(:,:,:))
      call diagnose_3d_u(Time, Grd, id_visc_diverge, visc_diverge(:,:,:))
      call diagnose_3d_u(Time, Grd, id_aaniso, wrk1(:,:,:))

      if (id_along > 0)  call diagnose_3d_u(Time, Grd, id_along,  aiso(:,:,:)+0.5*wrk1(:,:,:))
      if (id_across > 0) call diagnose_3d_u(Time, Grd, id_across, aiso(:,:,:)-0.5*wrk1(:,:,:))

      if (id_horz_bih_diss > 0)  then
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk2(i,j,k) = wrk2(i,j,k) - 4.0*massqc(i,j,k)*tmplap1(i,j,k)**2  &
                                             - 4.0*massqc(i,j,k)*tmplap2(i,j,k)**2 
                   wrk2(i,j,k) = wrk2(i,j,k)*Grd%daur(i,j)
                enddo
             enddo
          enddo
          call diagnose_3d_u(Time, Grd, id_horz_bih_diss, wrk2(:,:,:))
      endif 

      if(use_side_drag_friction) then
         call diagnose_3d_u(Time, Grd, id_bih_fric_u, wrk3_v(:,:,:,1))
         call diagnose_3d_u(Time, Grd, id_bih_fric_v, wrk3_v(:,:,:,2))
      else 
         call diagnose_3d_u(Time, Grd, id_bih_fric_u, wrk1_v(:,:,:,1))
         call diagnose_3d_u(Time, Grd, id_bih_fric_v, wrk1_v(:,:,:,2))
      endif

      call diagnose_3d_u(Time, Grd, id_bih_plus_side_fric_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_bih_plus_side_fric_v, wrk1_v(:,:,:,2))
      call diagnose_3d_u(Time, Grd, id_side_drag_friction_u, wrk2_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_side_drag_friction_v, wrk2_v(:,:,:,2))

  endif    ! endif for energy_analysis_step


  if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_bihgen_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('bihgen friction(1)', wrk1_v(COMP,:,1)*Grd%umask(COMP,:))
      call write_chksum_3d('bihgen friction(2)', wrk1_v(COMP,:,2)*Grd%umask(COMP,:))
  endif 


end subroutine bihgen_friction
! </SUBROUTINE> NAME="bihgen_friction"


!#######################################################################
! <SUBROUTINE NAME="ncar_boundary_scale_read">
!
! <DESCRIPTION>
!
! Read in the 3d ncar boundary scaling field and use this to 
! rescale the background viscosities. 
! 
! To use this routine, we need to already have generated the field
! ncar_rescale using the routine ncar_boundary_scale_create.
!
! The advantage of reading ncar_rescale is that we do not need to 
! introduce any global 2d arrays required for ncar_boundary_scale_create.     
! So the idea is to pay the price once by running ncar_boundary_scale_create,
! save ncar_rescale, then read that field in during subsequent runs through 
! ncar_boundary_scale_read.
!
! Here are the steps:
! 1/ run one time with ncar_boundary_scaling_read=.false.
! and ncar_boundary_scaling=.true. 
! Be sure that the field ncar_rescale is saved in diagnostic table.
! To ensure answers agree whether reading ncar_rescale or creating it
! during initialization, it is necessary to save ncar_rescale using the
! double precision option in the diagnostic table (packing=1). 
!
! 2/ extract field ncar_rescale from the diagnostics output
! and place into its own file INPUT/ncar_rescale.nc
! example extraction using ncks:
! ncks -v ncar_rescale 19900101.ocean_month.nc ncar_rescale.nc
!
! 3/ set ncar_boundary_scaling_read=.true. 
! and ncar_boundary_scaling=.true., and now run the model 
! reading in ncar_rescale rather than regenerating
! it during each initialization (which can be a bottleneck 
! for large models on huge processor counts).   
!
! 4/ As a check that all is fine, save ncar_rescale as a diagnostic
! for both the create and the read stage and make sure they agree.  
! Also, all checksums should agree whether reading in ncar_rescale
! or creating it each initialization, so long as the ncar_rescale.nc
! was saved with double precision  (see step 1/ above).
!
! </DESCRIPTION>
!
subroutine ncar_boundary_scale_read(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: i, j, k
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihgen_friction_mod (ncar_boundary_scale_read): module must be initialized')
  endif 

  call read_data('INPUT/ncar_rescale.nc','ncar_rescale',ncar_rescale, Dom%domain2d)

  ! rescale the background viscosities 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           aiso_back(i,j,k)   = aiso_back(i,j,k)  *ncar_rescale(i,j,k)
           aaniso_back(i,j,k) = aaniso_back(i,j,k)*ncar_rescale(i,j,k)
        enddo
     enddo
  enddo

  id_ncar_rescale = register_static_field('ocean_model', 'ncar_rescale', Grd%vel_axes_uv(1:3),&
                    'rescaling used for the background viscosity', 'dimensionless',           &
                     missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_3d_u(Time, Grd, id_ncar_rescale, ncar_rescale(:,:,:))

  ! need these viscosities in the halo regions
  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


end subroutine ncar_boundary_scale_read
! </SUBROUTINE>  NAME="ncar_boundary_scale_read"


!#######################################################################
! <SUBROUTINE NAME="ncar_boundary_scale_create">
!
! <DESCRIPTION>
!
!     Recale the background viscosities to be larger in the western 
!     boundary regions.  The algorithm is taken directly from the 
!     anisotropic_ncar routine in ocean_lapgen_friction.F90.
!
!   NOTE: The nearest western boundary computations are done along the
!         model i-grid lines. Therefore, viscosity based on these are 
!         only approximate in the high Northern Hemisphere when using 
!         generalized coordinates with coordinate pole(s) shifted onto 
!         land. 
!
! </DESCRIPTION>
!
subroutine ncar_boundary_scale_create(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: i, j, k
  integer :: n, ncount
  integer :: is, ie, ip1, ii
  integer :: index, indexo
  integer, dimension(nx) :: iwp
  integer, dimension(nx) :: nwbp_global

  real, dimension(nx,jsc:jec)       :: dxtn_global
  real, dimension(nx,jsc:jec)       :: kmu_global
  real, dimension(isd:ied,jsd:jed)  :: kmu_tmp
  real, dimension(isd:ied,jsd:jed)  :: ncar_rescale2d
  real, dimension(nx)               :: dist_g 

  real :: dist_max
  real :: ncar_rescale_min
  real :: huge=1e20

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihgen_friction_mod (ncar_velocity_scale): module must be initialized')
  endif 

  if(ncar_boundary_scaling_read) return 

  dist_max         = 1.e10  ! distance for ACC region (cm)
  dxtn_global(:,:) = 0.0
  kmu_global(:,:)  = 0.0
  kmu_tmp(:,:)     = Grd%kmu(:,:)
  ncar_rescale(:,:,:) = 0.0
  ncar_rescale2d(:,:) = 0.0

  call mpp_global_field(Dom%domain2d,kmu_tmp,kmu_global, flags=XUPDATE)  
  call mpp_global_field(Dom%domain2d,Grd%dxtn,dxtn_global, flags=XUPDATE)

  do k=1,nk
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
           indexo = index + ncar_vconst_5
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
        do i=isc,iec
           ncar_rescale(i,j,k) = Grd%umask(i,j,k)*(1.0+exp(-(ncar_vconst_4*dist_g(i+Dom%ioff))**2))
        enddo

     enddo ! j
  enddo    ! k

  ! the only depth dependence to ncar_rescale arises from absence or presence of 
  ! topography.  so to compute the min ncar_rescale, we only need to look at k=1.  
  ! also, to avoid getting ncar_rescale=0.0 due to land, we set ncar_rescale2d
  ! to a huge number if over land. 
  do j=jsc,jec
     do i=isc,iec
        if(Grd%umask(i,j,1)== 0.0) then 
           ncar_rescale2d(i,j) = huge
        else 
           ncar_rescale2d(i,j) = ncar_rescale(i,j,1) 
        endif 
     enddo
  enddo
  ncar_rescale_min = mpp_global_min(Dom%domain2d,ncar_rescale2d(:,:))

  ! check for error
  if(ncar_rescale_min==0.0) then
    call mpp_error(FATAL, &
    '==>Error from ocean_bihgen_friction_mod(ncar_boundary_scale): ncar_rescale_min=0.  Something wrong')
  else
    write(stdoutunit,'(a,e14.6)') 'From ncar_boundary_scale, minimum ncar_rescale =  ',ncar_rescale_min
  endif 

  ! rescale the background viscosities 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           ncar_rescale(i,j,k) = (ncar_rescale(i,j,k)/ncar_rescale_min)**ncar_rescale_power
           aiso_back(i,j,k)    = aiso_back(i,j,k)  *ncar_rescale(i,j,k)
           aaniso_back(i,j,k)  = aaniso_back(i,j,k)*ncar_rescale(i,j,k)
        enddo
     enddo
  enddo

  id_ncar_rescale = register_static_field('ocean_model', 'ncar_rescale', Grd%vel_axes_uv(1:3),&
                    'rescaling used for the background viscosity', 'dimensionless',           &
                     missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_3d_u(Time, Grd, id_ncar_rescale, ncar_rescale(:,:,:))

  ! need these viscosities in the halo regions
  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


end subroutine ncar_boundary_scale_create
! </SUBROUTINE>  NAME="ncar_boundary_scale_create"



!#######################################################################
! <FUNCTION NAME="BDX_EU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in X of a quantity defined on the east 
! face of a U-cell. Slightly modified version of BDX_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i-1/2,j).
!
! BDX_EU_smag changes dimensions by m^-3 
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
      BDX_EU_smag(i,j) = (Grd%dyue_dxuer(i,j)*a(i,j) - Grd%dyue_dxuer(i-1,j)*a(i-1,j))*daur_dyur(i,j)
    enddo
    BDX_EU_smag(isd,j) = 0.0
  enddo

end function BDX_EU_smag
! </FUNCTION> NAME="BDX_EU_smag"


!#######################################################################
! <FUNCTION NAME="BDY_NU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in Y of a quantity defined on the north
! face of a U-cell. Slightly modified version of BDY_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i,j-1/2).
!
! BDY_EU_smag changes dimensions by m^-3 
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
      BDY_NU_smag(i,j) = (Grd%dxun_dyunr(i,j)*a(i,j) - Grd%dxun_dyunr(i,j-1)*a(i,j-1))*daur_dxur(i,j)
    enddo
  enddo
  BDY_NU_smag(:,jsd) = 0.0

end function BDY_NU_smag
! </FUNCTION> NAME="BDY_EU_smag"


!#######################################################################
! <SUBROUTINE NAME="bihgen_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the biharmonic
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine bihgen_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: fudge

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihgen_friction_mod (bihgen_viscosity_check): needs initialization')
  endif 

  fudge=1.001  ! to eliminate roundoffs causing aiso=visc_crit+epsln

  write (stdoutunit,'(//60x,a/)') ' Excessive horizontal biharmonic friction summary:'
  write (stdoutunit,'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then            
          if (aiso(i,j,k) > fudge*visc_crit(i,j) .and. num < max_num ) then
            num = num + 1
            write (unit,9600) 'aiso(',aiso(i,j,k), visc_crit(i,j),i,j,k,Grd%xu(i,j),Grd%yu(i,j),Grd%zt(k),mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum (num)
  if (num > max_num) write (stdoutunit,*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es16.8,' m^4/s) exceeds max value (',es16.8,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')


end subroutine bihgen_viscosity_check
! </SUBROUTINE> NAME="bihgen_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="bihgen_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the biharmonic grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
subroutine bihgen_reynolds_check(Time, Velocity)

  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  real :: rame, ramn, reyx, reyy
  real :: reynx0, reyny0, reynx,reyny, reynu,reynv, reynmu,reynmv
  integer :: ireynx,jreynx,kreynx, ireyny,jreyny,kreyny
  integer :: i, j, k, tau
  integer, save :: unit=6

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihgen_friction_mod (gen_bih_reynolds_check): nees to be initialized')
  endif 

  if(.not. use_this_module) return 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(aiso(i,j,k) + epsln)
        ramn = 1.0/(aiso(i,j,k) + epsln)

        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j)**3)*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j)**3)*ramn
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
  write (stdoutunit,'(/60x,a/)') ' Horizontal biharmonic Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0) then
    write (unit,10300) reynx, ireynx, jreynx, kreynx, &
                       Grd%xu(ireynx,jreynx), Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, &
                       Grd%xu(ireyny,jreyny), Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine bihgen_reynolds_check
! </SUBROUTINE> NAME="bihgen_reynolds_check"


!#######################################################################
! <SUBROUTINE NAME="compute_neptune_velocity">
!
! <DESCRIPTION>
! Compute Neptune velocity.  
!
! Method follows that of 
! Maltrud and Holloway, 2008: Implementing biharmonic neptune in a 
! global eddying ocean model, Ocean Modelling, vol. 21, pages 22-34. 
!
! March 2012
! Stephen.Griffies 
!
! </DESCRIPTION>
!
subroutine compute_neptune_velocity(Time)

  type(ocean_time_type), intent(in) :: Time
  real, dimension(isd:ied,jsd:jed)  :: neptune_fu
  real, dimension(isd:ied,jsd:jed)  :: neptune_length
  real                              :: active_cells
  integer                           :: i,j,k
  integer                           :: num_smooth

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihgen_friction_mod (compute_neptune_velocity): needs initialization')
  endif 

  allocate (neptune_velocity(isd:ied,jsd:jed,2))
  neptune_velocity = 0.0
  if(.not. neptune) return 

  ! length scale on U-cell point 
  neptune_length(:,:) = 0.0
  neptune_fu(:,:)     = 0.0
  do j=jsd,jed
    do i=isd,ied
      neptune_fu(i,j)     = 2.0*omega_earth*sin(Grd%phiu(i,j))
      neptune_length(i,j) = neptune_scaling                                        &
        *( neptune_length_pole                                                     &
      + 0.5*(neptune_length_eq-neptune_length_pole)*(1.0 + cos(2.0*Grd%phiu(i,j))) &
         )
    enddo
  enddo

  ! compute velocity 
  do j=jsc,jec
     do i=isc,iec
        neptune_velocity(i,j,1) = -neptune_fu(i,j)*neptune_length(i,j)*neptune_length(i,j) &
                                   *Grd%dht_dy(i,j)*Grd%umask(i,j,1)/(neptune_depth_min+Grd%hu(i,j))
        neptune_velocity(i,j,2) =  neptune_fu(i,j)*neptune_length(i,j)*neptune_length(i,j) &
                                   *Grd%dht_dx(i,j)*Grd%umask(i,j,1)/(neptune_depth_min+Grd%hu(i,j))
     enddo
  enddo
  call mpp_update_domains(neptune_velocity(:,:,1), neptune_velocity(:,:,2),&
                          Dom%domain2d,gridtype=BGRID_NE)

  ! optional smoothing
  if(neptune_smooth) then 
      do num_smooth=1,neptune_smooth_num

         k=1
         wrk1_v2d(:,:,:) = 0.0
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
  id_neptune_bih_u = register_static_field('ocean_model', 'neptune_bih_u',                &
                   Grd%vel_axes_uv(1:2), 'Zonal velocity from neptune biharmonic scheme', &
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d_u(Time, Grd, id_neptune_bih_u, neptune_velocity(:,:,1))

  id_neptune_bih_v = register_static_field('ocean_model', 'neptune_bih_v',                     &
                   Grd%vel_axes_uv(1:2), 'Meridional velocity from neptune biharmonic scheme', &
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d_u(Time, Grd, id_neptune_bih_v, neptune_velocity(:,:,2))

  id_neptune_fu   = register_static_field('ocean_model', 'neptune_fu',                            &
                    Grd%tracer_axes(1:2),'Coriolis parameter used for  neptune parameterization', &
                    '1/sec', missing_value=missing_value, range=(/-1.e3,1.e3/))
  call diagnose_2d(Time, Grd, id_neptune_fu, neptune_fu(:,:))

end subroutine compute_neptune_velocity
! </SUBROUTINE> NAME="compute_neptune_velocity"


end module ocean_bihgen_friction_mod
      
      




