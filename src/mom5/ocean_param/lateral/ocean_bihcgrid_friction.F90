module ocean_bihcgrid_friction_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes the thickness weighted time tendency for  
! horizontal velocity arising from horizontal biharmonic friction. 
! Friction is formulated for the C-grid here. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted time tendency for
! horizontal velocity arising from horizontal biharmonic friction. 
! The viscosity used to determine the strength of the tendency 
! can be a general function of space and time as specified by 
! the Smagorinsky approach; a grid-scale dependent
! background viscosity; or other options.  
! The form of the friction operator can be isotropic or 
! anisotropic in the horizontal plane. 
!
! Friction is formulated for the C-grid in this module. 
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
! enabled at the same time.  Such has been found useful for some simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bihcgrid_friction_nml">
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
use fms_mod,          only: read_data
use mpp_domains_mod,  only: mpp_update_domains, CGRID_NE
use mpp_domains_mod,  only: mpp_start_update_domains, mpp_complete_update_domains
use mpp_domains_mod,  only: mpp_global_min, mpp_global_max, mpp_global_field, XUPDATE 
use mpp_mod,          only: input_nml_file, mpp_sum, mpp_pe, mpp_error, mpp_max
use mpp_mod,          only: FATAL, NOTE, stdout, stdlog

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: BAY, BAX, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_operators_mod,  only: FAY, FAX, FDX_NT, FDY_ET, FDX_T, FDY_T 
use ocean_parameters_mod, only: missing_value, onehalf, oneeigth
use ocean_parameters_mod, only: omega_earth, rho0, rho0r
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type
use ocean_types_mod,      only: ocean_domain_type, ocean_adv_vel_type 
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,      only: ocean_options_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d, diagnose_3d_u, diagnose_2d, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4
use ocean_workspace_mod,  only: wrk1_2d, wrk1_v2d  
use ocean_workspace_mod,  only: wrk1_v, wrk2_v, wrk3_v

implicit none

private

public ocean_bihcgrid_friction_init
public bihcgrid_friction
public bihcgrid_viscosity_check
public bihcgrid_reynolds_check

private compute_neptune_velocity

! length of forward time step  
real :: dtime

! for diagnostics 
integer :: id_aiso                 =-1
integer :: id_aaniso               =-1
integer :: id_aiso_smag            =-1
integer :: id_aaniso_smag          =-1
integer :: id_aiso_back            =-1
integer :: id_aaniso_back          =-1
integer :: id_along_back           =-1
integer :: id_across_back          =-1
integer :: id_aiso_diverge         =-1
integer :: id_sin2theta            =-1
integer :: id_cos2theta            =-1
integer :: id_neptune_bih_u        =-1
integer :: id_neptune_bih_v        =-1
integer :: id_ncar_rescale         =-1
integer :: id_along                =-1
integer :: id_across               =-1
integer :: id_visc_diverge         =-1
integer :: id_visc_crit_bih        =-1
integer :: id_bih_fric_u           =-1
integer :: id_bih_fric_v           =-1
integer :: id_horz_bih_diss        =-1
integer :: id_tmask_next_to_land   =-1
integer :: id_stress_xx_bih        =-1
integer :: id_stress_xy_bih        =-1
integer :: id_stress_yx_bih        =-1

logical :: used

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


! for divergence-based viscosity
real :: visc_diverge_scaling = 0.0  ! dimensionless
real, dimension(:,:),   allocatable :: visc_diverge_factor ! (m^4) for computing visc_diverge
real, dimension(:,:,:), allocatable :: divergence_t        ! array to hold diverge_t 
real, dimension(:,:,:), allocatable :: visc_diverge        ! viscosity based on horz divergence

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

! for side-boundary drag law
logical :: use_side_drag_friction       = .false.
real    :: side_drag_friction_scaling   = 1.0
real    :: side_drag_friction_max       = 1.0
real    :: side_drag_friction_uvmag_max = 10.0
real, dimension(:,:,:), allocatable :: tmask_next_to_land


#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,nk)  :: aiso_back      ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk)  :: aaniso_back    ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk)  :: aiso           ! isotropic viscosity (m^4/s) for stability check
real, dimension(isd:ied,jsd:jed,nk)  :: aiso_xx        ! isotropic viscosity (m^4/s) for stress_xx 
real, dimension(isd:ied,jsd:jed,nk)  :: aiso_xy        ! isotropic viscosity (m^4/s) for stress_xy 
real, dimension(isd:ied,jsd:jed,nk)  :: aaniso_xx      ! anisotropic viscosity (m^4/s) for stress_xx 
real, dimension(isd:ied,jsd:jed,nk)  :: aaniso_xy      ! anisotropic viscosity (m^4/s) for stress_xy 
real, dimension(isd:ied,jsd:jed,nk)  :: horz_ten       ! horizontal rate of deformation (1/s)
real, dimension(isd:ied,jsd:jed,nk)  :: horz_str       ! horizontal rate of strain (1/s)
real, dimension(isd:ied,jsd:jed,nk)  :: sin2theta      ! orientation angle 
real, dimension(isd:ied,jsd:jed,nk)  :: cos2theta      ! orientation angle 
real, dimension(isd:ied,jsd:jed,nk)  :: stress_xx      ! stress tensor component
real, dimension(isd:ied,jsd:jed,nk)  :: stress_xy      ! stress tensor component
real, dimension(isd:ied,jsd:jed,nk)  :: stress_yx      ! stress tensor component
real, dimension(isd:ied,jsd:jed,nk)  :: tmplap1        ! for holding the Laplacian friction 
real, dimension(isd:ied,jsd:jed,nk)  :: tmplap2        ! for holding the Laplacian friction 
real, dimension(isd:ied,jsd:jed,nk)  :: ncar_rescale   ! viscosity scale (cm^2/s) based on NCAR western boundary scaling 


#else

real, dimension(:,:,:),  allocatable  :: aiso_back      ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(:,:,:),  allocatable  :: aaniso_back    ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(:,:,:),  allocatable  :: aiso           ! isotropic viscosity (m^4/s) for diagnostic 
real, dimension(:,:,:),  allocatable  :: aiso_xx        ! isotropic viscosity (m^4/s) for stress_xx 
real, dimension(:,:,:),  allocatable  :: aiso_xy        ! isotropic viscosity (m^4/s) for stress_xy 
real, dimension(:,:,:),  allocatable  :: aaniso_xx      ! anisotropic viscosity (m^4/s) for stress_xx 
real, dimension(:,:,:),  allocatable  :: aaniso_xy      ! anisotropic viscosity (m^4/s) for stress_xy 
real, dimension(:,:,:),  allocatable  :: horz_ten       ! horizontal rate of tension (1/s)
real, dimension(:,:,:),  allocatable  :: horz_str       ! horizontal rate of strain (1/s)
real, dimension(:,:,:),  allocatable  :: sin2theta      ! orientation angle 
real, dimension(:,:,:),  allocatable  :: cos2theta      ! orientation angle 
real, dimension(:,:,:),  allocatable  :: stress_xx      ! stress tensor component
real, dimension(:,:,:),  allocatable  :: stress_xy      ! stress tensor component
real, dimension(:,:,:),  allocatable  :: stress_yx      ! stress tensor component
real, dimension(:,:,:),  allocatable  :: tmplap1        ! for holding the Laplacian friction 
real, dimension(:,:,:),  allocatable  :: tmplap2        ! for holding the Laplacian friction 
real, dimension(:,:,:),  allocatable  :: ncar_rescale   ! viscosity scale (cm^2/s) based on NCAR western boundary scaling 

#endif

real, dimension(:,:),   allocatable  :: aiso_back_obc     ! (m^4/s) for OBC 
real, dimension(:,:),   allocatable  :: aaniso_back_obc   ! (m^4/s) for OBC 
real, dimension(:,:),   allocatable  :: fsmag_iso         ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(:,:),   allocatable  :: fsmag_aniso       ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(:,:),   allocatable  :: visc_crit         ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(:,:),   allocatable  :: dxu_dyur          ! (dxu/dyu)
real, dimension(:,:),   allocatable  :: dyu_dxur          ! (dyu/dxu)
real, dimension(:,:),   allocatable  :: dxt_dytr          ! (dxt/dyt)
real, dimension(:,:),   allocatable  :: dyt_dxtr          ! (dyt/dxt)
real, dimension(:,:),   allocatable  :: grid_length       ! horizontal grid length (m)
real, dimension(:,:,:), allocatable  :: neptune_velocity  ! barotropic velocity (m/s) from Neptune 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

character(len=256) :: version=&
     '=>Using: ocean_bihcgrid_friction.F90 ($Id: ocean_bihcgrid_friction.F90,v 20.0 2013/12/14 00:14:12 fms Exp $)'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc              = .FALSE.
logical :: read_aiso_bih_back    = .false.


namelist /ocean_bihcgrid_friction_nml/ use_this_module, debug_this_module,              &
    k_smag_iso, k_smag_aniso,                                                           &
    vel_micom_iso, vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom, &
    equatorial_zonal, equatorial_zonal_lat,                                             &
    neptune, neptune_length_eq, neptune_length_pole, neptune_depth_min,                 &
    neptune_smooth, neptune_smooth_num,                                                 &
    ncar_boundary_scaling, ncar_rescale_power, ncar_vconst_4, ncar_vconst_5,            &
    ncar_boundary_scaling_read,                                                         &
    use_side_drag_friction, side_drag_friction_scaling,                                 &
    side_drag_friction_uvmag_max, side_drag_friction_max 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bihcgrid_friction_init">

! <DESCRIPTION>
! Initialize the lateral biharmonic friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_bihcgrid_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                        obc, use_bihcgrid_friction, debug)
 
  type(ocean_grid_type),    target, intent(in)    :: Grid
  type(ocean_domain_type),  target, intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: d_time
  logical,                          intent(in)    :: obc
  logical,                          intent(inout) :: use_bihcgrid_friction
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real :: coeff_iso, coeff_aniso, dxdy

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihcgrid_friction_mod (ocean_bihcgrid_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)


  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_bihcgrid_friction_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bihcgrid_friction_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_bihcgrid_friction_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_bihcgrid_friction_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_bihcgrid_friction_nml)  
  write (stdlogunit,ocean_bihcgrid_friction_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  if(use_this_module) then 
      call mpp_error(NOTE, '==> NOTE: USING ocean_bihcgrid_friction_mod.')
      Ocean_options%horz_bih_friction = 'Used general horizontal biharmonic friction for Cgrid.'
      use_bihcgrid_friction=.true.
  else 
      call mpp_error(NOTE, '==> NOTE: NOT using ocean_bihcgrid_friction_mod.')
      Ocean_options%horz_bih_friction = 'Did NOT use horizontal biharmonic friction for Cgrid.'
      use_bihcgrid_friction=.false.
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
   '==> Note from ocean_bihcgrid_friction_mod: using forward time step of (secs)', dtime 

  have_obc = obc
  if (have_obc) then 
     write(stdoutunit,'(/a,f10.2)') '==> Note from ocean_bihcgrid_friction_mod: considering obc' 
  endif 

  if (Dom%xhalo > 1 .or. Dom%yhalo > 1) then
    write (stdoutunit,'(/1x,a)') ' ==> NOTE: General biharmonic friction allows xhalo=yhalo=1.'
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
  allocate (aiso_back(isd:ied,jsd:jed,nk))
  allocate (aaniso_back(isd:ied,jsd:jed,nk))
  allocate (aiso(isd:ied,jsd:jed,nk))
  allocate (aiso_xx(isd:ied,jsd:jed,nk))
  allocate (aiso_xy(isd:ied,jsd:jed,nk))
  allocate (aaniso_xx(isd:ied,jsd:jed,nk))
  allocate (aaniso_xy(isd:ied,jsd:jed,nk))
  allocate (horz_ten(isd:ied,jsd:jed,nk))
  allocate (horz_str(isd:ied,jsd:jed,nk))
  allocate (sin2theta(isd:ied,jsd:jed,nk))
  allocate (cos2theta(isd:ied,jsd:jed,nk))
  allocate (stress_xx(isd:ied,jsd:jed,nk))  
  allocate (stress_xy(isd:ied,jsd:jed,nk))  
  allocate (stress_yx(isd:ied,jsd:jed,nk))  
  allocate (tmplap1(isd:ied,jsd:jed,nk))    
  allocate (tmplap2(isd:ied,jsd:jed,nk))    
  allocate (ncar_rescale(isd:ied,jsd:jed,nk))
#endif

  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (dxu_dyur(isd:ied,jsd:jed))
  allocate (dyu_dxur(isd:ied,jsd:jed))
  allocate (dxt_dytr(isd:ied,jsd:jed))
  allocate (dyt_dxtr(isd:ied,jsd:jed))
  allocate (grid_length(isd:ied,jsd:jed))

  allocate (aiso_back_obc(isd:ied,jsd:jed))
  allocate (aaniso_back_obc(isd:ied,jsd:jed))


  grid_length(:,:)    = Grd%tmask(:,:,1)*2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(epsln+Grd%dxt(:,:)+Grd%dyt(:,:))
  aiso_back(:,:,:)    = 0.0
  aaniso_back(:,:,:)  = 0.0
  aiso(:,:,:)         = 0.0
  aiso_xx(:,:,:)      = 0.0
  aiso_xy(:,:,:)      = 0.0
  aaniso_xx(:,:,:)    = 0.0
  aaniso_xy(:,:,:)    = 0.0
  horz_ten(:,:,:)     = 0.0
  horz_str(:,:,:)     = 0.0
  sin2theta(:,:,:)    = 0.0
  cos2theta(:,:,:)    = 0.0
  stress_xx(:,:,:)    = 0.0
  stress_xy(:,:,:)    = 0.0
  stress_yx(:,:,:)    = 0.0
  tmplap1(:,:,:)      = 0.0
  tmplap2(:,:,:)      = 0.0
  aiso_back_obc(:,:)  = 0.0
  aaniso_back_obc(:,:)= 0.0


  ! dimensionless ratios used for friction operator 
  dxu_dyur(:,:) = Grd%dxu(:,:)*Grd%dyur(:,:)
  dyu_dxur(:,:) = Grd%dyu(:,:)*Grd%dxur(:,:)
  dxt_dytr(:,:) = Grd%dxt(:,:)*Grd%dytr(:,:)
  dyt_dxtr(:,:) = Grd%dyt(:,:)*Grd%dxtr(:,:)


  ! smag related arrays 
  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = 0.125*(k_smag_iso/pi)**2
  coeff_aniso      = 0.125*(k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = Grd%umask(:,:,1)*coeff_iso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4
  fsmag_aniso(:,:) = Grd%umask(:,:,1)*coeff_aniso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4


  aiso_back_obc(:,:)   = 0.0
  aaniso_back_obc(:,:) = 0.0
  do j=jsd,jed
    do i=isd,ied 
      aiso_back_obc(i,j)     = Grd%tmask(i,j,1)*vel_micom_iso* &
                               ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
      aaniso_back_obc(i,j)   = Grd%tmask(i,j,1)*vel_micom_aniso* &
                               ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
      if(abs(Grd%yu(i,j)) < eq_lat_micom) then
        aiso_back_obc(i,j)   = Grd%tmask(i,j,1)*eq_vel_micom_iso* &
                               ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
        aaniso_back_obc(i,j) = Grd%tmask(i,j,1)*eq_vel_micom_aniso* &
                               ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
      endif
    enddo  
  enddo

  aiso_back(:,:,:)    = 0.0
  aaniso_back(:,:,:)  = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied 

           aiso_back(i,j,k)     = Grd%tmask(i,j,k)*vel_micom_iso* &
                ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
           aaniso_back(i,j,k)   = Grd%tmask(i,j,k)*vel_micom_aniso* &
                ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3

           if(abs(Grd%yt(i,j)) < eq_lat_micom) then
               aiso_back(i,j,k)   = Grd%tmask(i,j,k)*eq_vel_micom_iso* &
                    ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
               aaniso_back(i,j,k) = Grd%tmask(i,j,k)*eq_vel_micom_aniso* &
                    ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**3
           endif

        enddo
     enddo
  enddo

  ! for computing viscosity based on horizontal divergence 
  allocate (visc_diverge(isd:ied,jsd:jed,nk))
  visc_diverge(:,:,:) = 0.0
  if(visc_diverge_scaling > 0.0) then 
      allocate (visc_diverge_factor(isd:ied,jsd:jed))
      allocate (divergence_t(isd:ied,jsd:jed,nk))
      visc_diverge_factor(:,:)= 0.0
      divergence_t(:,:,:)     = 0.0 
      write(stdoutunit,'(/a,f10.2)') &
       '==> Note from ocean_bihcgrid_friction_mod: using visc_diverge_scaling for viscosity, with a value = ',visc_diverge_scaling 
      do j=jsc,jec
         do i=isc,iec  
            visc_diverge_factor(i,j) = Grd%tmask(i,j,1)*visc_diverge_scaling* &
                                       ((2.0*Grd%dxt(i,j)*Grd%dyt(i,j))/(epsln + Grd%dxt(i,j) + Grd%dyt(i,j)))**4
         enddo
      enddo
      call mpp_update_domains (visc_diverge_factor(:,:),Dom%domain2d)    
  endif


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
      dxdy = 1.0/( 1.0/(Grd%dxt(i,j)*Grd%dxt(i,j)) + 1.0/(Grd%dyt(i,j)*Grd%dyt(i,j)) )
      visc_crit(i,j) =  visc_crit_scale*oneeigth*dxdy**2/(dtime+epsln)
    enddo
  enddo

  if(use_side_drag_friction) then 
     write(stdoutunit,'(a)') '==> Note from ocean_bihcgrid_friction: using side drag friction for cells next to land.'
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
     id_tmask_next_to_land = register_static_field('ocean_model', 'tmask_next_to_land_bih',&
            Grd%tracer_axes(1:3),'mask for cells next to land for bihcgrid module',        &
            'dimensionless', missing_value=missing_value, range=(/-10.0,10.0/))
     call diagnose_3d(Time, Grd, id_tmask_next_to_land, tmask_next_to_land(:,:,:))
  endif 

  id_visc_crit_bih = register_static_field ('ocean_model', 'visc_crit_bih',&
                     Grd%tracer_axes(1:2), 'critical viscosity', 'm^4/sec',&
                     missing_value=missing_value, range=(/0.0,1.e20/))
  call diagnose_2d(Time, Grd, id_visc_crit_bih, visc_crit(:,:))


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

  ! apply T-cell masking to background viscosities; needed in particular 
  ! for the case when read_aiso_back is enabled, to remove any spurious 
  ! missing values.  
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           aiso_back(i,j,k)   = aiso_back(i,j,k)*Grd%tmask(i,j,k)
           aaniso_back(i,j,k) = aaniso_back(i,j,k)*Grd%tmask(i,j,k)
        enddo
     enddo
  enddo


  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


  ! diagnostic output 
  id_aiso    = register_diag_field ('ocean_model', 'aiso_bih', Grd%tracer_axes(1:3), &
               Time%model_time,'T-cell isotropic bih visc', 'm^4/sec',               &
               missing_value=-10.0, range=(/-10.0,1.e20/),                           &
               standard_name='ocean_momentum_xy_biharmonic_diffusivity')
  id_aaniso  = register_diag_field ('ocean_model', 'aaniso_bih', Grd%tracer_axes(1:3), &
               Time%model_time,'T-cell anisotropic bih visc', 'm^4/sec',               &
               missing_value=missing_value, range=(/-10.0,1.e20/))
  id_aiso_smag = register_diag_field ('ocean_model', 'aiso_bih_smag', Grd%tracer_axes(1:3), &
                 Time%model_time, 'T-cell isotropic bih visc from smag', 'm^4/sec',         &
                 missing_value=missing_value, range=(/-10.0,1.e10/))
  id_aaniso_smag = register_diag_field ('ocean_model', 'aaniso_bih_smag', Grd%tracer_axes(1:3), &
                   Time%model_time, 'T-cell anisotropic bih visc from smag', 'm^4/sec',         &
                   missing_value=missing_value, range=(/-10.0,1.e10/))
  id_along  = register_diag_field ('ocean_model', 'along_bih', Grd%tracer_axes(1:3),   &
              Time%model_time,'T-cell along-stream bih visc', 'm^4/sec',               &
              missing_value=missing_value, range=(/-10.0,1.e20/))
  id_across  = register_diag_field ('ocean_model', 'across_bih', Grd%tracer_axes(1:3), &
               Time%model_time, 'T-cell cross-stream bih visc', 'm^4/sec',             &
               missing_value=missing_value, range=(/-10.0,1.e20/))
  id_visc_diverge = register_diag_field ('ocean_model', 'visc_diverge_bih', Grd%tracer_axes(1:3),    &
               Time%model_time,'T-cell biharmonic visc based on horz_diverge_t gradients', 'm^4/sec',&
               missing_value=-10.0, range=(/-10.0,1.e20/))

  id_sin2theta = register_diag_field ('ocean_model', 'sin2theta_bih', Grd%tracer_axes(1:3),                &
              Time%model_time, 'sin2theta for orientation angle with biharmonic friction', 'dimensionless',&
              missing_value=missing_value, range=(/-10.0,10.0/))
  id_cos2theta = register_diag_field ('ocean_model', 'cos2theta_bih', Grd%tracer_axes(1:3),                &
              Time%model_time, 'cos2theta for orientation angle with biharmonic friction', 'dimensionless',&
              missing_value=missing_value, range=(/-10.0,10.0/))

  id_bih_fric_u = register_diag_field ('ocean_model', 'bih_fric_u', Grd%tracer_axes(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict on u-zonal',   &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_bih_fric_v = register_diag_field ('ocean_model', 'bih_fric_v', Grd%tracer_axes(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict on v-merid',   &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_horz_bih_diss = register_diag_field ('ocean_model', 'horz_bih_diss', Grd%tracer_axes(1:3),      &
                  Time%model_time, 'Energy dissipation from horizontal biharmonic friction', 'W/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_stress_xx_bih = register_diag_field ('ocean_model', 'stress_xx_bih', Grd%tracer_axes(1:3),   &
                  Time%model_time, 'Stress tensor xx component from biharmonic friction', 'N/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_stress_xy_bih = register_diag_field ('ocean_model', 'stress_xy_bih', Grd%vel_axes_uv(1:3),   &
                  Time%model_time, 'Stress tensor xy component from biharmonic friction', 'N/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_stress_yx_bih = register_diag_field ('ocean_model', 'stress_yx_bih', Grd%vel_axes_uv(1:3),   &
                  Time%model_time, 'Stress tensor xy component from biharmonic friction', 'N/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_aiso_back   = register_static_field('ocean_model', 'aiso_bih_back', Grd%tracer_axes(1:3),   &
                    'U-cell background aiso bih visc', 'm^4/sec',                                &
                    missing_value=missing_value, range=(/-10.0,1.e20/))
  id_aaniso_back = register_static_field('ocean_model', 'aaniso_bih_back', Grd%tracer_axes(1:3), &
                    'U-cell background aaniso bih visc', 'm^4/sec',                              &
                    missing_value=missing_value, range=(/-10.0,1.e20/))
  id_along_back   = register_static_field('ocean_model', 'along_bih_back', Grd%tracer_axes(1:3), &
                    'U-cell background along-stream bih visc', 'm^4/sec',                        &
                    missing_value=missing_value, range=(/-10.0,1.e20/))
  id_across_back = register_static_field('ocean_model', 'across_bih_back', Grd%tracer_axes(1:3), &
                   'U-cell background cross-stream bih visc', 'm^4/sec',                         &
                   missing_value=missing_value, range=(/-10.0,1.e20/))

  call diagnose_3d(Time, Grd, id_aiso_back, aiso_back(:,:,:))
  call diagnose_3d(Time, Grd, id_aaniso_back, aaniso_back(:,:,:))
  if (id_along_back > 0) then 
     call diagnose_3d(Time, Grd, id_along_back, aiso_back(:,:,:)+0.5*aaniso_back(:,:,:))
  endif 
  if (id_across_back > 0) then 
     call diagnose_3d(Time, Grd, id_across_back, aiso_back(:,:,:)-0.5*aaniso_back(:,:,:))
  endif 

end subroutine ocean_bihcgrid_friction_init
! </SUBROUTINE>  NAME="ocean_bihcgrid_friction_init"


!#######################################################################
! <SUBROUTINE NAME="bihcgrid_friction">
!
! <DESCRIPTION>
! This routine computes thickness weighted and density weighted 
! time tendency for horizontal velocity arising from horizontal 
! biharmonic friction.  
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
! As shown in Griffies and Hallberg (2000), 
! a biharmonic operator with a nonconstant viscosity is guaranteed to 
! dissipate kinetic energy *only* when using the sqrt of the biharmonic
! viscosity at each of the two stages of the algorithm. 
! The sqrt approach is employed here.  
!
! </DESCRIPTION>
!
subroutine bihcgrid_friction(Time, Thickness, Adv_vel, Velocity, bih_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step
  real, dimension(isd:,jsd:), intent(inout) :: bih_viscosity

  real :: deform, delta 
  real :: usqrd, vsqrd, umagr
  real :: uvel, vvel
  real :: uvmag, umag, vmag, velocity_gamma  
  real :: grad_diverge_sq
  real :: term1, term2 
  real :: volume 

  integer :: i, j, k, n
  integer :: taum1, tau
  integer :: stdoutunit, ibl
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
    '==>Error in ocean_bihcgrid_friction_mod (bihcgrid_friction): module not initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau



  ! divergence_t has units m/s.
  visc_diverge(:,:,:) = 0.0
  if(visc_diverge_scaling > 0.0) then 

     divergence_t(:,:,:) = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              divergence_t(i,j,k) = rho0r*Adv_vel%diverge_t(i,j,k)
           enddo
        enddo
     enddo
     call mpp_update_domains (divergence_t(:,:,:), Dom%domain2d)    

     ! ignore offset of u,v 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              grad_diverge_sq     = ((divergence_t(i+1,j,k)-divergence_t(i,j,k))*Grd%dxter(i,j))**2 &
                                   +((divergence_t(i,j+1,k)-divergence_t(i,j,k))*Grd%dytnr(i,j))**2 
              visc_diverge(i,j,k) = Grd%tmask(i,j,k)*visc_diverge_factor(i,j)*sqrt(grad_diverge_sq)
           enddo
        enddo
     enddo 

  endif 

  !-----------------------------------
  ! Laplacian portion of the algorithm 


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

  ! compute the horizontal deformation rates (1/sec)
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
  ! write diagnostics here since wrk arrays are reused later 
  if(.not. energy_analysis_step) then 
     if (id_aiso_smag > 0) then
        call diagnose_3d(Time, Grd, id_aiso_smag, onehalf*(wrk1(:,:,:)+wrk2(:,:,:)))
     endif 
     if (id_aaniso_smag > 0) then
        call diagnose_3d(Time, Grd, id_aaniso_smag, onehalf*(wrk3(:,:,:)+wrk4(:,:,:)))
     endif  
  endif 

  ! compute the full viscosities with dimension sqrt(m^4/s).
  ! take sqrt to ensure that non-constant biharmonic viscosity dissipates 
  ! (see Section 17.9.2 of Griffies book). 
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec
           aiso_xx(i,j,k)   = min(visc_crit(i,j), aiso_back(i,j,k)   + wrk1(i,j,k) + visc_diverge(i,j,k))
           aaniso_xx(i,j,k) = min(visc_crit(i,j), aaniso_back(i,j,k) + wrk3(i,j,k) + visc_diverge(i,j,k))
           aiso_xy(i,j,k)   = min(visc_crit(i,j), aiso_back(i,j,k)   + wrk2(i,j,k) + visc_diverge(i,j,k))
           aaniso_xy(i,j,k) = min(visc_crit(i,j), aaniso_back(i,j,k) + wrk4(i,j,k) + visc_diverge(i,j,k))
           aiso(i,j,k)      = aiso_xx(i,j,k)
           aiso_xx(i,j,k)   = sqrt(aiso_xx(i,j,k))
           aaniso_xx(i,j,k) = sqrt(aaniso_xx(i,j,k))
           aiso_xy(i,j,k)   = sqrt(aiso_xy(i,j,k))
           aaniso_xy(i,j,k) = sqrt(aaniso_xy(i,j,k))
        enddo
     enddo
  enddo

  ! compute orientation angle; ignore offset of u,v 
  sin2theta = 0.0  
  cos2theta = 0.0   
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec
           if(equatorial_zonal .and. abs(Grd%yt(i,j)) <= equatorial_zonal_lat) then
             sin2theta(i,j,k) = 0.0
             cos2theta(i,j,k) = 1.0
           else
             usqrd = wrk1_v(i,j,k,1)**2
             vsqrd = wrk1_v(i,j,k,2)**2
             umagr = 1.0/(epsln + usqrd + vsqrd)
             sin2theta(i,j,k) = 2.0*wrk1_v(i,j,k,1)*wrk1_v(i,j,k,2)*umagr
             cos2theta(i,j,k) = (usqrd-vsqrd)*umagr
           endif
        enddo
     enddo
  enddo


  ! compute the stress tensor components (kg/m^3)*(m^2/s^1.5)
  ! stress_xy(i,j) = 0.5*rho*
  !   [aiso_xy(i,j)*horz_str(i,j) - aaniso_xy(i,j)*cos2theta*0.5*(horz_str(i,j)*cos2theta - horz_ten(i+1,j+1)*sin2theta)]
  ! stress_xx(i,j) = rho*
  !   [aiso_xx(i,j)*horz_ten(i,j) + aaniso_xx(i,j)*sin2theta*0.5*(horz_str(i,j)*cos2theta - horz_ten(i,j)*sin2theta)] 
  ! set rho = rho0 even for non-Boussinesq 
  ! stress_yy = -stress_xx
  ! stress_xy =  stress_yx except for possible differences next to boundaries with partial slip 
  stress_xx(:,:,:) = 0.0
  stress_xy(:,:,:) = 0.0
  stress_yx(:,:,:) = 0.0
  
  ! umask applied to stress_xy corresponds to free-slip side boundaries 
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec

           ! note use of different indices for horz_ten
           delta = onehalf*(horz_str(i,j,k)*cos2theta(i,j,k)-horz_ten(i+1,j+1,k)*sin2theta(i,j,k))
           stress_xy(i,j,k) = rho0*Grd%umask(i,j,k) &
                              *(aiso_xy(i,j,k)*horz_str(i,j,k) - aaniso_xy(i,j,k)*cos2theta(i,j,k)*delta)
           stress_yx(i,j,k) = stress_xy(i,j,k) 

           delta = onehalf*(horz_str(i,j,k)*cos2theta(i,j,k)-horz_ten(i,j,k)*sin2theta(i,j,k))
           stress_xx(i,j,k) = rho0*grd%tmask(i,j,k) &
                              *(aiso_xx(i,j,k)*horz_ten(i,j,k) + aaniso_xx(i,j,k)*sin2theta(i,j,k)*delta)

        enddo
     enddo
  enddo

  ! need stress components in halos 
  call mpp_update_domains (stress_xx(:,:,:), Dom%domain2d)    
  call mpp_update_domains (stress_xy(:,:,:), Dom%domain2d)    
  call mpp_update_domains (stress_yx(:,:,:), Dom%domain2d)    


  ! compute Laplacian operator with units (kg/m^3)*(m/s^1.5)
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           term1 = Grd%dyter(i,j)*(dyt_dxtr(i+1,j)*stress_xx(i+1,j,k)-dyt_dxtr(i,j)  *stress_xx(i,j,k))  
           term2 = Grd%dxter(i,j)*(dxu_dyur(i,j)  *stress_xy(i,j,k)  -dxu_dyur(i,j-1)*stress_xy(i,j-1,k))  
           tmplap1(i,j,k) = term1 + term2

           term1 = -Grd%dxtnr(i,j)*(dxt_dytr(i,j+1)*stress_xx(i,j+1,k) - dxt_dytr(i,j)  *stress_xx(i,j,k))  
           term2 =  Grd%dytnr(i,j)*(dyu_dxur(i,j)  *stress_yx(i,j,k)   - dyu_dxur(i-1,j)*stress_yx(i-1,j,k))  
           tmplap2(i,j,k) = term1 + term2
        enddo
     enddo
  enddo
  call mpp_update_domains (tmplap1(:,:,:), Dom%domain2d, complete=.false.)    
  call mpp_update_domains (tmplap2(:,:,:), Dom%domain2d, complete=.true. )    


  !-------------------------------------- 
  ! second portion of the algorithm, with velocity replaced by tmplap

  ! store laplacian operator (kg/m^3)*(m/s^1.5)
  wrk1_v(:,:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk1_v(i,j,k,1) = tmplap1(i,j,k)
           wrk1_v(i,j,k,2) = tmplap2(i,j,k)
        enddo
     enddo
  enddo

  ! compute the horizontal deformation rates with units (kg/m^3)*(1/s^1.5)
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


  ! compute the stress tensor components with units N/m^2 = (kg/m^3)*(m^2/s^2)
  ! umask applied to stress_xy corresponds to free-slip side boundaries 
  stress_xx(:,:,:) = 0.0
  stress_xy(:,:,:) = 0.0
  stress_yx(:,:,:) = 0.0
  
  do k=1,nk 
     do j=jsc,jec
        do i=isc,iec

           delta = onehalf*(horz_str(i,j,k)*cos2theta(i,j,k)-horz_ten(i+1,j+1,k)*sin2theta(i,j,k))
           stress_xy(i,j,k) = Grd%umask(i,j,k) &
                              *(aiso_xy(i,j,k)*horz_str(i,j,k) - aaniso_xy(i,j,k)*cos2theta(i,j,k)*delta)
           stress_yx(i,j,k) = stress_xy(i,j,k) 

           delta = onehalf*(horz_str(i,j,k)*cos2theta(i,j,k)-horz_ten(i,j,k)*sin2theta(i,j,k))
           stress_xx(i,j,k) = grd%tmask(i,j,k) &
                              *(aiso_xx(i,j,k)*horz_ten(i,j,k) + aaniso_xx(i,j,k)*sin2theta(i,j,k)*delta)

        enddo
     enddo
  enddo

   ! side drag law for partial slip on stress_xy and stress_yx (N/m^2)
  wrk2_v(:,:,:,:) = 0.0
  if(use_side_drag_friction) then
     do k=1,nk    
        do j=jsc,jec
           do i=isc,iec
              uvel             = Velocity%u(i,j,k,1,taum1)
              vvel             = Velocity%u(i,j,k,2,taum1)
              uvmag            = sqrt(uvel**2 + vvel**2)
              uvmag            = min(uvmag,side_drag_friction_uvmag_max)
              velocity_gamma   = rho0*side_drag_friction_scaling*Velocity%cdbot_array(i,j)*uvmag
              umag             = abs(uvel)
              vmag             = abs(vvel)
              wrk2_v(i,j,k,1)  = -min(side_drag_friction_max, velocity_gamma*umag)*sign(1.0,uvel)
              wrk2_v(i,j,k,2)  = -min(side_drag_friction_max, velocity_gamma*vmag)*sign(1.0,vvel)
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


  ! compute biharmonic operator with units (kg/m^3)*(m/sec^2) = N/m^2.
  ! the factor of dzt ensures proper units, and scales down friction for thin cells.  
  wrk3_v(:,:,:,:) = 0.0
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



  ! compute some diagnostics 

  bih_viscosity(:,:) = 0.0
  wrk1_2d(:,:)       = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%rho_dzt(i,j,k,tau)
           bih_viscosity(i,j) = bih_viscosity(i,j)  &
                              + Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*aiso_xx(i,j,k)
           aiso_xx(i,j,k)   = aiso_xx(i,j,k)**2
           aiso_xy(i,j,k)   = aiso_xy(i,j,k)**2
           aaniso_xx(i,j,k) = aaniso_xx(i,j,k)**2
           aaniso_xy(i,j,k) = aaniso_xy(i,j,k)**2
        enddo
     enddo
  enddo


  ! for sending back to driver 
  do j=jsc,jec
     do i=isc,iec
        bih_viscosity(i,j) = Grd%tmask(i,j,1)*bih_viscosity(i,j)/(epsln+wrk1_2d(i,j))
     enddo
  enddo
  call mpp_update_domains (bih_viscosity(:,:), Dom%domain2d)    
  
  call diagnose_3d(Time, Grd, id_sin2theta, sin2theta(:,:,:))
  call diagnose_3d(Time, Grd, id_cos2theta, cos2theta(:,:,:))

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

   call diagnose_3d_u(Time, Grd, id_bih_fric_u, wrk3_v(:,:,:,1))
   call diagnose_3d_u(Time, Grd, id_bih_fric_v, wrk3_v(:,:,:,2))

   if (id_horz_bih_diss > 0) then 
      wrk1(:,:,:) = 0.0  
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               volume = Grd%tmask(i,j,k)*Grd%dat(i,j)*Thickness%dzt(i,j,k)
               wrk1(i,j,k) = wrk1(i,j,k) - rho0r*volume*(tmplap1(i,j,k)**2 + tmplap2(i,j,k))**2 
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_horz_bih_diss, wrk1(:,:,:))
   endif 

   call diagnose_3d(Time, Grd, id_stress_xx_bih, stress_xx(:,:,:))
   call diagnose_3d_u(Time, Grd, id_stress_xy_bih, stress_xy(:,:,:))
   call diagnose_3d_u(Time, Grd, id_stress_yx_bih, stress_yx(:,:,:))

   if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_bihcgrid_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('bihcgrid friction(1)', wrk3_v(COMP,:,1)*Grd%tmask(COMP,:))
      call write_chksum_3d('bihcgrid friction(2)', wrk3_v(COMP,:,2)*Grd%tmask(COMP,:))
   endif


end subroutine bihcgrid_friction
! </SUBROUTINE> NAME="bihcgrid_friction"


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

  ncar_rescale(:,:,:) = 0.0

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihcgrid_friction_mod (ncar_boundary_scale_read): module must be initialized')
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

  id_ncar_rescale = register_static_field('ocean_model', 'ncar_rescale', Grd%tracer_axes(1:3),&
                    'rescaling used for the background viscosity', 'dimensionless',           &
                     missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_3d(Time, Grd, id_ncar_rescale, ncar_rescale(:,:,:))

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

  real, dimension(nx,jsc:jec)       :: dxt_global
  real, dimension(nx,jsc:jec)       :: kmt_global
  real, dimension(isd:ied,jsd:jed)  :: kmt_tmp
  real, dimension(isd:ied,jsd:jed)  :: ncar_rescale2d
  real, dimension(nx)               :: dist_g 

  real :: dist_max
  real :: ncar_rescale_min
  real :: huge=1e20

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihcgrid_friction_mod (ncar_velocity_scale): module must be initialized')
  endif 

  if(ncar_boundary_scaling_read) return 

  dist_max         = 1.e10  ! distance for ACC region (cm)
  dxt_global(:,:)  = 0.0
  kmt_global(:,:)  = 0.0
  kmt_tmp(:,:)     = Grd%kmt(:,:)
  ncar_rescale(:,:,:) = 0.0
  ncar_rescale2d(:,:) = 0.0

  call mpp_global_field(Dom%domain2d,kmt_tmp,kmt_global, flags = XUPDATE)  
  call mpp_global_field(Dom%domain2d,Grd%dxt,dxt_global, flags = XUPDATE)

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
           if(kmt_global(i,j)<k .and. kmt_global(ip1,j) >= k) then 
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
               dist_g(i) = dxt_global(i,j)+dist_g(i-1)
           elseif ( i .lt. index ) then
               if (indexo .le. Grd%ni) then
                   if (i .eq. 1) then
                       dist_g(i) = 0.0
                       do ii=indexo+1,Grd%ni
                          dist_g(i)=dxt_global(ii,j) + dist_g(i)
                       enddo
                       dist_g(i) = dxt_global(i,j)+dist_g(i)
                   else
                       dist_g(i) = dxt_global(i,j)+dist_g(i-1)
                   endif
               else
                   if (i .le. (indexo - Grd%ni)) then
                       dist_g(i) = 0.0
                   else
                       dist_g(i) = dxt_global(i,j)+dist_g(i-1)
                   endif
               endif
           endif

        enddo

        ! convert distance in MOM (m) to POP1.4 (cm)
        dist_g(:) = 100.0*dist_g(:)
        do i=isc,iec
           ncar_rescale(i,j,k) = Grd%tmask(i,j,k)*(1.0+exp(-(ncar_vconst_4*dist_g(i+Dom%ioff))**2))
        enddo

     enddo ! j
  enddo    ! k

  ! the only depth dependence to ncar_rescale arises from absence or presence of 
  ! topography.  so to compute the min ncar_rescale, we only need to look at k=1.  
  ! also, to avoid getting ncar_rescale=0.0 due to land, we set ncar_rescale2d
  ! to a huge number if over land. 
  do j=jsc,jec
     do i=isc,iec
        if(Grd%tmask(i,j,1)== 0.0) then 
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
    '==>Error from ocean_bihcgrid_friction_mod(ncar_boundary_scale): ncar_rescale_min=0.  Something wrong')
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

  id_ncar_rescale = register_static_field('ocean_model', 'ncar_rescale', Grd%tracer_axes(1:3),&
                    'rescaling used for the background viscosity', 'dimensionless',           &
                     missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_3d(Time, Grd, id_ncar_rescale, ncar_rescale(:,:,:))

  ! need these viscosities in the halo regions
  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


end subroutine ncar_boundary_scale_create
! </SUBROUTINE>  NAME="ncar_boundary_scale_create"



!#######################################################################
! <SUBROUTINE NAME="bihcgrid_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the biharmonic
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine bihcgrid_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: fudge

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihcgrid_friction_mod (bihcgrid_viscosity_check): needs initialization')
  endif 

  fudge=1.001  ! to eliminate roundoffs causing aiso=visc_crit+epsln

  write (stdoutunit,'(//60x,a/)') ' Excessive horizontal biharmonic friction summary:'
  write (stdoutunit,'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%tmask(i,j,k) > 0.0) then            
          if (aiso(i,j,k) > fudge*visc_crit(i,j) .and. num < max_num ) then
            num = num + 1
            write (unit,9600) 'aiso(',aiso(i,j,k), visc_crit(i,j),i,j,k,Grd%xt(i,j),Grd%yt(i,j),Grd%zt(k),mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum (num)
  if (num > max_num) write (stdoutunit,*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es16.8,' m^4/s) exceeds max value (',es16.8,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')


end subroutine bihcgrid_viscosity_check
! </SUBROUTINE> NAME="bihcgrid_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="bihcgrid_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the biharmonic grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
subroutine bihcgrid_reynolds_check(Time, Velocity)

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
    call mpp_error(FATAL, '==>Error from ocean_bihcgrid_friction_mod (gen_bih_reynolds_check): nees to be initialized')
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

        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxt(i,j)**3)*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyt(i,j)**3)*ramn
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
                       Grd%xt(ireynx,jreynx), Grd%yt(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, &
                       Grd%xt(ireyny,jreyny), Grd%yt(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine bihcgrid_reynolds_check
! </SUBROUTINE> NAME="bihcgrid_reynolds_check"


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
! This approach is slightly different than the Eby and Holloway
! method implemented in the laplacian module. There is no fundamental
! reason to favor one versus the other.  We use the Maltrud and 
! Holloway method here sinc they implemented it for biharmonic. 
!
! May 2012
! Stephen.Griffies 
!
! </DESCRIPTION>
!
subroutine compute_neptune_velocity(Time)

  type(ocean_time_type), intent(in) :: Time

  real, dimension(isd:ied,jsd:jed)  :: dhu_dx
  real, dimension(isd:ied,jsd:jed)  :: dhu_dy
  real, dimension(isd:ied,jsd:jed)  :: neptune_length
  real                              :: active_cells_u
  real                              :: active_cells_v
  integer                           :: i,j,k
  integer                           :: num_smooth

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihcgrid_friction_mod (compute_neptune_velocity): needs initialization')
  endif 

  allocate (neptune_velocity(isd:ied,jsd:jed,2))
  neptune_velocity = 0.0
  if(.not. neptune) return 

  ! length scale on U-cell point 
  neptune_length(:,:) = 0.0
  do j=jsd,jed
    do i=isd,ied
      neptune_length(i,j) = neptune_scaling                                        &
        *( neptune_length_pole                                                     &
      + 0.5*(neptune_length_eq-neptune_length_pole)*(1.0 + cos(2.0*Grd%phiu(i,j))) &
         )
    enddo
  enddo

  ! gradient of U-cell bottom depth 
  do j=jsc,jec
     do i=isc,iec
        dhu_dx(i,j) = (Grd%hu(i,j)-Grd%hu(i-1,j))*Grd%dxuer(i-1,j)    ! positioned at v(i,j)
        dhu_dy(i,j) = (Grd%hu(i,j)-Grd%hu(i,j-1))*Grd%dyunr(i-1,j)    ! positioned at u(i,j)
     enddo
  enddo

  ! computed neptune velocity 
  do j=jsc,jec
     do i=isc,iec
        neptune_velocity(i,j,1) = -Grd%f(i,j)*neptune_length(i,j)*neptune_length(i,j) &
                                   *dhu_dy(i,j)*Grd%tmasken(i,j,1,1)/(neptune_depth_min+Grd%hu(i,j))
        neptune_velocity(i,j,2) =  Grd%f(i,j)*neptune_length(i,j)*neptune_length(i,j) &
                                   *dhu_dx(i,j)*Grd%tmasken(i,j,1,2)/(neptune_depth_min+Grd%hu(i,j))
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



  ! diagnostics 
  id_neptune_bih_u = register_static_field('ocean_model', 'neptune_bih_u',                &
                   Grd%tracer_axes(1:2), 'Zonal velocity from neptune biharmonic scheme', &
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d(Time, Grd, id_neptune_bih_u, neptune_velocity(:,:,1))

  id_neptune_bih_v = register_static_field('ocean_model', 'neptune_bih_v',                     &
                   Grd%tracer_axes(1:2), 'Meridional velocity from neptune biharmonic scheme', &
                   'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  call diagnose_2d(Time, Grd, id_neptune_bih_v, neptune_velocity(:,:,2))


end subroutine compute_neptune_velocity
! </SUBROUTINE> NAME="compute_neptune_velocity"



end module ocean_bihcgrid_friction_mod
      
      




