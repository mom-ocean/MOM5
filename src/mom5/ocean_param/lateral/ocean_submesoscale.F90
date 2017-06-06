module ocean_submesoscale_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes a streamfunction within
! the upper surface boundary layer, and applies this
! streamfunction to all tracers. It also optionally 
! applies horizontal diffusion in the surface layer 
! as determined by the strength of the streamfunction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes a streamfunction within
! the upper surface boundary layer, and applies this
! streamfunction to all tracers.  It also optionally 
! applies horizontal diffusion in the surface layer 
! as determined by the strength of the streamfunction. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Fox-Kemper, Ferrari, and Hallberg 2008: Parameterization of 
! mixed layer eddies. Part I: theory and diagnosis
! Journal of Physical Oceanography, vol. 38, pages 1145-1165. 
! </REFERENCE>
!
! <REFERENCE>
! Fox-Kemper, Danabasoglu, Ferrari, and Hallberg 2008: 
! Parameterizing submesoscale physics in global models.
! Clivar Exchanges, vol 13, no.1,  Jan2008. pages 3-5.
! </REFERENCE>
!
! <REFERENCE>
! Fox-Kemper, Danabasoglu, Ferrari, Griffies, Hallberg,
! Holland, Peacock, Samuels, 2011: Parameterization of 
! Mixed Layer Eddies. III: Global Implementation and 
! Impact on Ocean Climate Simulations, Ocean Modelling, 
! vol. 39, pages 61-78.
! </REFERENCE>
!
! <REFERENCE>
! Griffies, 2012: Elements of MOM 
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_submesoscale_nml">
!  <DATA NAME="use_this_module=" TYPE="logical">
!  Must be .true. to use this module.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!  <DATA NAME="diag_step" TYPE="integer">
!  Number of time steps between computing max bottom value for
!  wrho_bt_submeso. Default diag_step=1200.
!  </DATA> 
!
!  <DATA NAME="submeso_skew_flux" TYPE="logical">
!  For computing the tendency as convergence of skew flux.
!  This is the recommended method.
!  Default submeso_skew_flux=.true.
!  </DATA> 
!
!  <DATA NAME="submeso_advect_flux" TYPE="logical">
!  For computing the tendency as convergence of advective flux.
!  This approach uses either a flux limited sweby advection or 
!  first order upwind, both of which ensure that the resulting 
!  tendency will not create extrema in the tracer field.  
!  Default submeso_advect_flux=.false.
!  </DATA> 
!  <DATA NAME="submeso_advect_upwind" TYPE="logical">
!  For computing the tendency as convergence of a first order
!  advective flux. 
!  Default submeso_advect_upwind=.true.
!  </DATA> 
!  <DATA NAME="submeso_advect_sweby" TYPE="logical">
!  For computing the tendency as convergence of a sweby 
!  advective flux. This routine is incomplete and has a bug.  
!  Default submeso_advect_sweby=.false.
!  </DATA> 
!  <DATA NAME="submeso_advect_limit" TYPE="logical">
!  For limiting the value of the horizontal transports 
!  to be less than a velocity scale set by limit_psi_velocity_scale.
!  This option is not needed if limit_psi=.true.
!  Default submeso_advect_limit=.false.
!  </DATA> 
!  <DATA NAME="submeso_advect_zero_bdy" TYPE="logical">
!  For removing the advective transport next to boundaries. 
!  This is useful since computation of the advective transport
!  velocity components can be problematic next to boundaries. 
!  Default submeso_advect_zero_bdy=.false.
!  </DATA> 
!  <DATA NAME="smooth_advect_transport" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the diagnosed  
!  uhrho_et_submeso and vhrho_nt_submeso fields.  
!  Default smooth_advect_transport=.true.
!  </DATA> 
!  <DATA NAME="smooth_advect_transport_num" TYPE="integer">
!  Number of iterations for the smooothing of horizontal transport. 
!  Default smooth_advect_transport_num=2.
!  </DATA> 
!
!  <DATA NAME="submeso_diffusion" TYPE="logical">
!  For computing a horizontal diffusive flux in the boundary layer
!  as determined by the strength of the vector streamfunction.
!  Default submeso_diffusion=.false.
!  </DATA> 
!  <DATA NAME="submeso_diffusion_biharmonic" TYPE="logical">
!  The default submeso diffusion is Laplacian. However, one may wish to
!  use a biharmonic mixing operator instead.  
!  Default submeso_diffusion_biharmonic=.false. 
!  </DATA> 
!  <DATA NAME="submeso_diffusion_scale" UNITS="dimensionless" TYPE="real">
!  A dimensionless scaling to be used for scaling up or down the effects from 
!  horizontal diffusion in the boundary layer. Default submeso_diffusion_scale=1.0.
!  </DATA> 
!
!  <DATA NAME="use_hblt_constant" TYPE="logical">
!  For running with a constant boundary layer depth. This for the case when 
!  not using a realistic mixed layer scheme.  Default use_hblt_constant=.false.
!  </DATA> 
!  <DATA NAME="constant_hblt" UNITS="metre" TYPE="real">
!  The boundary layer depth for the case when use_hblt_constant=.true.
!  Default constant_hblt=100.0.
!  </DATA> 
!  <DATA NAME="use_hblt_equal_mld" TYPE="logical">
!  For using the diagnosed mld as the hblt for submeso.  
!  This is useful for those test models that do not have a mixed layer
!  scheme enabled, such as KPP, where the mixed layer scheme provides a
!  boundary layer depth.  In this case, it is sensible to employ the diagnosed
!  mixed layer depth for the submeso scheme. Additionally, in general it is 
!  more physical to use the mld than the KPP hblt as the depth over which 
!  the submesoscale eddies act.  Hence, default use_hblt_equal_mld=.true.
!  </DATA> 
!  <DATA NAME="min_kblt" UNITS="dimensionless" TYPE="integer">
!  The minimum number of vertical cells in the surface boundary layer 
!  that are required in order to compute the submesoscale streamfunction.
!  Default min_kblt=4.  Need at least three to fit a parabola with zero 
!  streamfunction at the top and bottom of the boundary layer.  
!  </DATA> 
!
!  <DATA NAME="minimum_hblt" TYPE="real" UNITS="metre">
!  For setting a floor to the hblt used for submesoscale scheme. 
!  Default minimum_hblt=0.0.
!  </DATA> 
!
!  <DATA NAME="smooth_hblt" TYPE="logical">
!  For smoothing on the submeso bldepth field. This is useful 
!  since the bldepth obtained from KPP or diagnosed mld can 
!  have some grid noise. 
!  Default smooth_hblt=.false. since this agrees with legacy.
!  Note that this scheme fails to reproduce across
!  processor layout, so it remains broken.  
!  </DATA> 
!  <DATA NAME="smooth_hblt_num" TYPE="integer">
!  Number of iterations for the smooothing of bldepth. 
!  Default smooth_hblt_num=1.
!  </DATA> 
!
!  <DATA NAME="use_psi_legacy" TYPE="logical">
!  For computing psi using older legacy methods. 
!  These methods are not ideal, and can be problematic
!  depending on nml settings for the limiters and smoothers.
!  This option is retained only for legacy purposes. 
!  Default use_psi_legacy=.false.
!  </DATA> 
!  <DATA NAME="smooth_psi" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the 
!  psix_horz and psiy_horz fields. 
!  Default smooth_psi=.true.
!  </DATA> 
!  <DATA NAME="smooth_psi_num" TYPE="integer">
!  Number of iterations for the smooothing of psi. 
!  Default smooth_psi_num=2.
!  </DATA> 
!  <DATA NAME="limit_psi" TYPE="logical">
!  For limiting the magnitude of psi in order to reduce possibility of 
!  model crashes. Rescales the full psi to maintain vertical structure
!  but to keep overall magnitude within bounds.  
!  Default limit_psi=.false.
!  </DATA> 
!  <DATA NAME="limit_psi_velocity_scale" UNITS="metre/sec" TYPE="real">
!  Velocity scale used to limit the value of psi when limit_psi=.true.
!  Default limit_psi_velocity_scale=5.0
!  </DATA> 
!
!  <DATA NAME="submeso_limit_flux" TYPE="logical">
!  For limiting the fluxes arising from submeso scheme, according to 
!  tmask_limit. When reach a point where tmask_limit=1.0, then set
!  the submeso flux for this cell to zero. 
!  Default submeso_limit_flux=.true.
!  </DATA> 
!
! <DATA NAME="coefficient_ce" UNITS="dimensionless" TYPE="real">
!  The dimensionless coefficient from the Fox-Kemper etal scheme. 
!  They recommend setting coefficient_ce between 0.06 and 0.08.  
!  Default coefficient_ce=0.07.
!  </DATA> 
! <DATA NAME="time_constant" UNITS="seconds" TYPE="real">
!  Timescale to mix momentum across the mixed layer.  
!  Default time_constant=86400.0 = 1day. 
!  </DATA> 
! <DATA NAME="front_length_const" UNITS="metre" TYPE="real">
!  Take constant horizontal length scale of submesoscale front. 
!  Default front_length_const=5e3.
!  </DATA> 
! <DATA NAME="front_length_deform_radius" TYPE="logical">
!  To compute the front length using the mixed layer deformation 
!  radius. Default front_length_deform_radius=.true.  Note, 
!  will have a floor on the variable front length set by the
!  nml setting for front_length_const.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field, register_static_field, need_data, send_data
use fms_mod,           only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,           only: stdout, stdlog, read_data, NOTE, FATAL, WARNING
use mpp_domains_mod,   only: mpp_update_domains, XUPDATE, YUPDATE, CGRID_NE
use mpp_domains_mod,   only: mpp_global_sum, NON_BITWISE_EXACT_SUM
use mpp_mod,           only: input_nml_file, mpp_error, mpp_max, mpp_pe
use time_manager_mod,  only: time_type, increment_time

use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_operators_mod,  only: FMX, FMY, FDX_T, FDY_T, BDX_ET, BDY_NT
use ocean_parameters_mod, only: missing_value, DEPTH_BASED, omega_earth, grav
use ocean_parameters_mod, only: rho0, rho0r, onehalf, onefourth, onesixth, oneeigth
use ocean_tracer_diag_mod,only: calc_mixed_layer_depth, diagnose_eta_tend_3dflux 
use ocean_types_mod,      only: tracer_2d_type, tracer_3d_0_nk_type, tracer_3d_1_nk_type
use ocean_types_mod,      only: ocean_time_type, ocean_domain_type, ocean_grid_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_thickness_type, ocean_density_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_sum, write_chksum_3d
use ocean_tracer_util_mod,only: diagnose_3d_rho
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk1_2d, wrk2_2d, wrk1_v

implicit none

private

! for diagnostics 
real, dimension(:,:,:),   allocatable :: advect_tendency 

! internally set for computing watermass diagnstics
logical :: compute_watermass_diag      = .false. 
logical :: compute_watermass_diff_diag = .false. 

! for diagnostics 
integer :: id_kblt_submeso         =-1
integer :: id_hblt_submeso         =-1
integer :: id_mu_submeso           =-1
integer :: id_dmu_submeso          =-1
integer :: id_psix_submeso         =-1
integer :: id_psiy_submeso         =-1
integer :: id_tx_trans_submeso     =-1
integer :: id_ty_trans_submeso     =-1
integer :: id_tz_trans_submeso     =-1
integer :: id_tx_trans_submeso_adv =-1
integer :: id_ty_trans_submeso_adv =-1
integer :: id_tz_trans_submeso_adv =-1
integer :: id_tx_trans_nrho_submeso=-1
integer :: id_ty_trans_nrho_submeso=-1
integer :: id_tz_trans_nrho_submeso=-1
integer :: id_front_length_submeso =-1
integer :: id_buoy_freq_ave_submeso=-1 
integer :: id_uhrho_et_submeso     =-1
integer :: id_vhrho_nt_submeso     =-1
integer :: id_wrho_bt_submeso      =-1
integer :: id_u_et_submeso         =-1
integer :: id_v_nt_submeso         =-1
integer :: id_w_bt_submeso         =-1
integer :: id_subdiff_diffusivity  =-1

integer :: id_eta_tend_submeso_flx      =-1
integer :: id_eta_tend_submeso_flx_glob =-1
integer :: id_eta_tend_submeso_tend     =-1
integer :: id_eta_tend_submeso_tend_glob=-1

integer :: id_eta_tend_subdiff_flx      =-1
integer :: id_eta_tend_subdiff_flx_glob =-1
integer :: id_eta_tend_subdiff_tend     =-1
integer :: id_eta_tend_subdiff_tend_glob=-1

integer :: id_neut_rho_submeso          =-1
integer :: id_neut_rho_submeso_on_nrho  =-1
integer :: id_wdian_rho_submeso         =-1
integer :: id_wdian_rho_submeso_on_nrho =-1
integer :: id_tform_rho_submeso         =-1
integer :: id_tform_rho_submeso_on_nrho =-1

integer :: id_neut_temp_submeso          =-1
integer :: id_neut_temp_submeso_on_nrho  =-1
integer :: id_wdian_temp_submeso         =-1
integer :: id_wdian_temp_submeso_on_nrho =-1
integer :: id_tform_temp_submeso         =-1
integer :: id_tform_temp_submeso_on_nrho =-1

integer :: id_neut_salt_submeso          =-1
integer :: id_neut_salt_submeso_on_nrho  =-1
integer :: id_wdian_salt_submeso         =-1
integer :: id_wdian_salt_submeso_on_nrho =-1
integer :: id_tform_salt_submeso         =-1
integer :: id_tform_salt_submeso_on_nrho =-1

integer :: id_neut_rho_subdiff          =-1
integer :: id_neut_rho_subdiff_on_nrho  =-1
integer :: id_wdian_rho_subdiff         =-1
integer :: id_wdian_rho_subdiff_on_nrho =-1
integer :: id_tform_rho_subdiff         =-1
integer :: id_tform_rho_subdiff_on_nrho =-1

integer :: id_neut_temp_subdiff          =-1
integer :: id_neut_temp_subdiff_on_nrho  =-1
integer :: id_wdian_temp_subdiff         =-1
integer :: id_wdian_temp_subdiff_on_nrho =-1
integer :: id_tform_temp_subdiff         =-1
integer :: id_tform_temp_subdiff_on_nrho =-1

integer :: id_neut_salt_subdiff          =-1
integer :: id_neut_salt_subdiff_on_nrho  =-1
integer :: id_wdian_salt_subdiff         =-1
integer :: id_wdian_salt_subdiff_on_nrho =-1
integer :: id_tform_salt_subdiff         =-1
integer :: id_tform_salt_subdiff_on_nrho =-1


integer, dimension(:), allocatable :: id_xflux_submeso       ! i-directed flux 
integer, dimension(:), allocatable :: id_yflux_submeso       ! j-directed flux 
integer, dimension(:), allocatable :: id_zflux_submeso       ! k-directed flux 
integer, dimension(:), allocatable :: id_xflux_submeso_int_z ! vertically integrated i-flux
integer, dimension(:), allocatable :: id_yflux_submeso_int_z ! vertically integrated j-flux
integer, dimension(:), allocatable :: id_submeso             ! tendency from streamfunction portion of submesoscale param 
integer, dimension(:), allocatable :: id_submeso_on_nrho     ! tendency from streamfunction portion of submesoscale param binned to neutral density

integer, dimension(:), allocatable :: id_xflux_subdiff       ! i-directed horz diffusive flux 
integer, dimension(:), allocatable :: id_yflux_subdiff       ! j-directed horz diffusive flux 
integer, dimension(:), allocatable :: id_xflux_subdiff_int_z ! vertically integrated horz diffusive i-flux
integer, dimension(:), allocatable :: id_yflux_subdiff_int_z ! vertically integrated horz diffusive j-flux
integer, dimension(:), allocatable :: id_subdiff             ! tendency from horz diffusion associated with submesoscale 

logical :: used


#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()
type(ocean_domain_type), save    :: Dom_flux_sub

character(len=128)  :: version='$$'
character (len=128) :: tagname = '$Name: tikal $'


#ifdef MOM_STATIC_ARRAYS
real, dimension(isd:ied,jsd:jed,nk)   :: uhrho_et_submeso ! i-component of transport for submeso 
real, dimension(isd:ied,jsd:jed,nk)   :: vhrho_nt_submeso ! j-component of transport for submeso 
real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bt_submeso  ! k-component of transport for submeso 

real, dimension(isd:ied,jsd:jed,nk,0:1) :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(isd:ied,jsd:jed,0:nk)   :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(isd:ied,jsd:jed,0:1)    :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(isd:ied,jsd:jed,nk,0:1) :: psix             ! streamfunction x-component (m^2/sec) 
real, dimension(isd:ied,jsd:jed,nk,0:1) :: psiy             ! streamfunction y-component (m^2/sec) 
real, dimension(isd:ied,jsd:jed)        :: hblt             ! boundary layer depth (m) 
real, dimension(isd:ied,jsd:jed)        :: grid_length      ! grid length scale (m)
real, dimension(isd:ied,jsd:jed)        :: front_length_inv ! inverse front length (1/m)
real, dimension(isd:ied,jsd:jed)        :: buoy_freq_ave    ! buoyancy frequency averaged over mixed layer depth (1/sec)
real, dimension(isd:ied,jsd:jed)        :: time_factor      ! time factor (sec) for computing streamfunction
real, dimension(isd:ied,jsd:jed)        :: coriolis_param   ! absolute value of the Coriolis parameter (sec^-1) on T-grid 
integer, dimension(isd:ied,jsd:jed)     :: kblt             ! k-level encompassing hblt 

real, dimension(isd:ied,jsd:jed,nk)     :: flux_x      ! i-component to tracer flux
real, dimension(isd:ied,jsd:jed,nk)     :: flux_y      ! j-component to tracer flux
real, dimension(isd:ied,jsd:jed,0:nk)   :: flux_z      ! k-component to tracer flux

#else 

real, dimension(:,:,:), allocatable :: uhrho_et_submeso  ! i-component of transport for submeso 
real, dimension(:,:,:), allocatable :: vhrho_nt_submeso  ! j-component of transport for submeso 
real, dimension(:,:,:), allocatable :: wrho_bt_submeso   ! k-component of transport for submeso 

real, dimension(:,:,:,:), allocatable :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(:,:,:),   allocatable :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(:,:,:),   allocatable :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),   allocatable :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),   allocatable :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(:,:,:),   allocatable :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(:,:,:,:), allocatable :: psix             ! streamfunction x-component (m^2/sec) 
real, dimension(:,:,:,:), allocatable :: psiy             ! streamfunction y-component (m^2/sec) 
real, dimension(:,:),     allocatable :: hblt             ! boundary layer depth (m) 
real, dimension(:,:),     allocatable :: grid_length      ! grid length scale (m)
real, dimension(:,:),     allocatable :: front_length_inv ! inverse front length (1/m)
real, dimension(:,:),     allocatable :: buoy_freq_ave    ! buoyancy frequency averaged over mixed layer depth (1/sec)
real, dimension(:,:),     allocatable :: time_factor      ! time factor (sec) for computing streamfunction
real, dimension(:,:),     allocatable :: coriolis_param   ! absolute value of the Coriolis parameter (sec^-1) on T-grid 
integer, dimension(:,:),  allocatable :: kblt             ! k-level encompassing hblt 

real, dimension(:,:,:),   allocatable :: flux_x      ! i-component to tracer flux
real, dimension(:,:,:),   allocatable :: flux_y      ! j-component to tracer flux
real, dimension(:,:,:),   allocatable :: flux_z      ! k-component to tracer flux


#endif 

! for saving the horizontally dependent portion of psix and psiy
real, dimension(:,:,:), allocatable :: psix_horz ! streamfunction x-component sans vertical structure (m^2/sec) 
real, dimension(:,:,:), allocatable :: psiy_horz ! streamfunction y-component sans vertical structure (m^2/sec) 

! for eta_tend diagnostics 
real, dimension(:,:,:), allocatable :: flux_x_temp
real, dimension(:,:,:), allocatable :: flux_y_temp
real, dimension(:,:,:), allocatable :: flux_z_temp
real, dimension(:,:,:), allocatable :: flux_x_salt
real, dimension(:,:,:), allocatable :: flux_y_salt
real, dimension(:,:,:), allocatable :: flux_z_salt

! for advecting tracers with sweby submesoscale advection 
real, dimension(:,:,:), allocatable :: tmask_mdfl
real, dimension(:,:,:), allocatable :: tracer_mdfl
type  :: tracer_mdfl_type
  real, dimension(:,:,:), pointer :: field => NULL()
end type tracer_mdfl_type
type(tracer_mdfl_type),  dimension(:), allocatable  :: tracer_mdfl_all ! tracer array for sweby advection   
type(ocean_domain_type), save :: Dom_mdfl


type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdx ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdy ! tracer partial derivative (tracer/m)
type(tracer_3d_0_nk_type), dimension(:), allocatable  :: dTdz ! tracer partial derivative (tracer/m)

type(tracer_3d_1_nk_type), save :: dSdx ! density-salinity partial derivative (tracer/m)
type(tracer_3d_1_nk_type), save :: dSdy ! density-salinity partial derivative (tracer/m)
type(tracer_3d_0_nk_type), save :: dSdz ! density-salinity partial derivative (tracer/m)

public ocean_submesoscale_init 
public submeso_restrat
private tracer_derivs 
private salinity_derivs 
private compute_flux_x
private compute_flux_y
private compute_flux_z
private compute_psi
private compute_psi_legacy
private compute_bldepth 
private compute_advect_transport 
private compute_submeso_skewsion
private compute_submeso_upwind 
private compute_submeso_sweby
private maximum_bottom_w_submeso
private transport_on_nrho_submeso
private transport_on_nrho_submeso_adv
private watermass_diag_init
private watermass_diag
private watermass_diag_diffusion

integer :: index_temp=-1
integer :: index_salt=-1 
integer :: num_prog_tracers=0

! for diagnosing fluxes
real    :: flux_sign

! vertical coordinate 
integer :: vert_coordinate_class=1

! for output
integer :: unit=6

real    :: fiveover21
real    :: eightover21
real    :: dtime 
real    :: ce_grav_rho0r 
real    :: time_constant2_r
real    :: grav_rho0r
real    :: front_length_const_inv
real    :: front_length_max=1e10

! for global area normalization
real    :: cellarea_r

logical :: module_is_initialized=.false.

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 

! for eta_tend diagnostics (internally set)
logical :: diagnose_eta_tend_submeso_flx=.false. 
logical :: diagnose_eta_tend_subdiff_flx=.false.

! internally determined 
logical :: diag_advect_transport = .false.

! nml parameters 
logical :: use_this_module               = .false.
logical :: debug_this_module             = .false. 
logical :: submeso_skew_flux             = .true.
logical :: submeso_advect_flux           = .false.
logical :: submeso_advect_upwind         = .true. 
logical :: submeso_advect_sweby          = .false. 
logical :: submeso_advect_limit          = .false. 
logical :: submeso_advect_zero_bdy       = .false.
logical :: submeso_limit_flux            = .true.
logical :: use_hblt_constant             = .false.    
logical :: use_hblt_equal_mld            = .true. 
logical :: smooth_hblt                   = .false.    
logical :: smooth_psi                    = .true.    
logical :: smooth_advect_transport       = .true. 
logical :: front_length_deform_radius    = .true.
logical :: limit_psi                     = .true.
logical :: use_psi_legacy                = .false.
logical :: submeso_diffusion             = .false.
logical :: submeso_diffusion_biharmonic  = .false.
real    :: limit_psi_velocity_scale      = 0.5
real    :: coefficient_ce                = 0.07
real    :: time_constant                 = 86400.0 
real    :: front_length_const            = 5e3 
real    :: constant_hblt                 = 100.0
real    :: minimum_hblt                  = 0.0
real    :: submeso_diffusion_scale       = 1.0 
integer :: min_kblt                      = 4
integer :: diag_step                     = 1200
integer :: smooth_psi_num                = 2
integer :: smooth_advect_transport_num   = 2
integer :: smooth_hblt_num               = 2

namelist /ocean_submesoscale_nml/ use_this_module, debug_this_module, diag_step,        &
                                  use_hblt_constant, use_hblt_equal_mld,                &
                                  smooth_hblt, smooth_hblt_num, constant_hblt,          &
                                  coefficient_ce, time_constant, front_length_const,    &
                                  min_kblt, minimum_hblt, smooth_psi, smooth_psi_num,   &
                                  front_length_deform_radius,                           &
                                  limit_psi, use_psi_legacy, limit_psi_velocity_scale,  &
                                  submeso_limit_flux, smooth_advect_transport,          &
                                  smooth_advect_transport_num,                          &
                                  submeso_skew_flux, submeso_advect_flux,               &
                                  submeso_advect_upwind, submeso_advect_sweby,          &
                                  submeso_advect_limit, submeso_advect_zero_bdy,        &
                                  submeso_diffusion, submeso_diffusion_scale,           &
                                  submeso_diffusion_biharmonic  

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_submesoscale_init">
!
! <DESCRIPTION>
! Initialization for the ocean_submesoscale module.
! </DESCRIPTION>
  subroutine ocean_submesoscale_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, dtime_t, &
                                     ver_coordinate_class, cmip_units, debug)
  
    type(ocean_grid_type),        intent(in), target   :: Grid
    type(ocean_domain_type),      intent(in), target   :: Domain
    type(ocean_time_type),        intent(in)           :: Time
    type(ocean_density_type),     intent(in)           :: Dens
    type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
    type(ocean_options_type),     intent(inout)        :: Ocean_options 
    real,                         intent(in)           :: dtime_t
    integer,                      intent(in)           :: ver_coordinate_class 
    logical,                      intent(in)           :: cmip_units
    logical,                      intent(in), optional :: debug

    integer :: unit, io_status, ierr
    integer :: i,j,k,n
    integer :: num_methods=0

    integer :: stdoutunit,stdlogunit 
    stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
      call mpp_error(FATAL,&
      '==>Error from ocean_submesoscale_init_mod: module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_submesoscale_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_submesoscale_nml')
#else
    unit = open_namelist_file()
    read(unit, ocean_submesoscale_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_submesoscale_nml')
    call close_file(unit)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_submesoscale_nml)    
    write (stdlogunit,ocean_submesoscale_nml)

    Dom => Domain
    Grd => Grid

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif

    if(cmip_units) then
        transport_convert=1.0
        transport_dims   = 'kg/s'
    else
        transport_convert=1.0e-9 
        transport_dims   = 'Sv (10^9 kg/s)'
    endif

    if (PRESENT(debug) .and. .not. debug_this_module) then
        debug_this_module = debug
    endif
    if(debug_this_module) then 
        write(stdoutunit,'(a)') '==>Note: running ocean_submesoscale_mod with debug_this_module=.true.'  
    endif

    if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING ocean_submesoscale_mod')
      Ocean_options%submesoscale = 'Used submesoscale closure for surface restratification.'
    else 
      call mpp_error(NOTE, '==>Note: NOT using ocean_submesoscale_mod')
      Ocean_options%submesoscale = 'Did NOT use submesoscale closure for surface restratification.'
      return 
    endif 

    if(use_hblt_equal_mld) then 
      write(stdoutunit,'(a)') &
      '==>Note: For ocean_submesoscale, setting bldepth equal to diagnosed mld.'
    endif 
    if(use_hblt_constant) then 
      write(stdoutunit,'(a)') &
      '==>Note: For ocean_submesoscale, setting bldepth equal to prescribed constant.'
    endif 
    if(use_hblt_equal_mld .and. use_hblt_constant) then 
      call mpp_error(FATAL, &
      '==>Error: in ocean_submesoscale, use_hblt_equal_mld & use_hblt_constant cannot both be true.')
    endif 

    if(submeso_advect_flux) then 
        num_methods=num_methods+1
        flux_sign = 1.0
        if(submeso_advect_upwind) then 
            write(stdoutunit,'(a)') &
                 '==>Note: For ocean_submesoscale, tendency computed as convergence of upwind advective flux.'
        elseif(submeso_advect_sweby) then 
            write(stdoutunit,'(a)') &
                 '==>Note: For ocean_submesoscale, tendency computed as convergence of mdfl_sweby advective flux; has known bugs!'
        else
            call mpp_error(FATAL, &
            '==>Error: in ocean_submesoscale: submeso_advect_flux=.true. yet no advection scheme chosen')
        endif
        if(use_psi_legacy) then
            call mpp_error(FATAL, &
            '==>Error: in ocean_submesoscale: submeso_advect_flux=.true. is not compatible with use_psi_legacy=.true.')
        endif 
    endif
    if(submeso_skew_flux) then 
       write(stdoutunit,'(a)') &
       '==>Note: For ocean_submesoscale, tendency computed as skew flux convergence.'
       num_methods=num_methods+1
       flux_sign = -1.0
    endif 
    if(num_methods > 1) then 
       call mpp_error(FATAL, &
       '==>Error: in ocean_submesoscale, can choose only one method for computing tendency.')
    endif 
    if(num_methods == 0) then 
       call mpp_error(FATAL, &
       '==>Error: in ocean_submesoscale, must choose one method for computing tendency.')
    endif 

    if(use_psi_legacy) then
       write(stdoutunit,'(a)') &
       '==>Note: in ocean_submesoscale: use_psi_legacy=.true. This method has problems, and is retained solely for legacy.'  
    endif 

    if(submeso_diffusion) then
        if(submeso_diffusion_biharmonic) then 
            write(stdoutunit,'(a)') &
            '==>Note: in ocean_submesoscale: submeso_diffusion=.true. Adding horizontal biharmonic mixing in blayer.'  
        else 
            write(stdoutunit,'(a)') &
            '==>Note: in ocean_submesoscale: submeso_diffusion=.true. Adding horizontal Laplacian diffusion in blayer.'  
        endif
    endif


    fiveover21             = 5.0/21.0
    eightover21            = 8.0/21.0
    grav_rho0r             = grav*rho0r
    ce_grav_rho0r          = coefficient_ce*grav*rho0r
    time_constant2_r       = 1.0/(time_constant**2) 
    dtime                  = dtime_t 
    vert_coordinate_class  = ver_coordinate_class
    front_length_const_inv = 1.0/front_length_const
    cellarea_r             = 1.0/(epsln + Grd%tcellsurf)

    call set_ocean_domain(Dom_mdfl,Grd,xhalo=2,yhalo=2,name='mdfl',maskmap=Dom%maskmap)

    call set_ocean_domain(Dom_flux_sub, Grd,xhalo=Dom%xhalo,yhalo=Dom%yhalo,name='flux dom submeso',&
                          maskmap=Dom%maskmap)

    ! for diagnostics 
    allocate( advect_tendency(isd:ied,jsd:jed,nk) )
    advect_tendency(:,:,:) = 0.0

    num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo

    allocate( dTdx(num_prog_tracers) )
    allocate( dTdy(num_prog_tracers) )
    allocate( dTdz(num_prog_tracers) )

#ifndef MOM_STATIC_ARRAYS
    allocate (uhrho_et_submeso(isd:ied,jsd:jed,nk))
    allocate (vhrho_nt_submeso(isd:ied,jsd:jed,nk))
    allocate (wrho_bt_submeso(isd:ied,jsd:jed,0:nk))

    allocate (delqc(isd:ied,jsd:jed,nk,0:1))
    allocate (dzwtr(isd:ied,jsd:jed,0:nk))
    allocate (dtew(isd:ied,jsd:jed,0:1))
    allocate (dtns(isd:ied,jsd:jed,0:1))
    allocate (dtwedyt(isd:ied,jsd:jed,0:1))
    allocate (dxtdtsn(isd:ied,jsd:jed,0:1))

    allocate (psix(isd:ied,jsd:jed,nk,0:1))
    allocate (psiy(isd:ied,jsd:jed,nk,0:1))
    allocate (hblt(isd:ied,jsd:jed))
    allocate (grid_length(isd:ied,jsd:jed))
    allocate (front_length_inv(isd:ied,jsd:jed))
    allocate (buoy_freq_ave(isd:ied,jsd:jed))
    allocate (time_factor(isd:ied,jsd:jed))
    allocate (coriolis_param(isd:ied,jsd:jed))
    allocate (kblt(isd:ied,jsd:jed))

    allocate (flux_x(isd:ied,jsd:jed,nk) )
    allocate (flux_y(isd:ied,jsd:jed,nk) )
    allocate (flux_z(isd:ied,jsd:jed,0:nk) )

    do n=1,num_prog_tracers
       allocate ( dTdx(n)%field(isd:ied,jsd:jed,nk) )
       allocate ( dTdy(n)%field(isd:ied,jsd:jed,nk) )
       allocate ( dTdz(n)%field(isd:ied,jsd:jed,0:nk) )
    enddo
    allocate ( dSdx%field(isd:ied,jsd:jed,nk) )
    allocate ( dSdy%field(isd:ied,jsd:jed,nk) )
    allocate ( dSdz%field(isd:ied,jsd:jed,0:nk) )

#endif 

    allocate (psix_horz(isd:ied,jsd:jed,0:1))
    allocate (psiy_horz(isd:ied,jsd:jed,0:1))

    uhrho_et_submeso = 0.0
    vhrho_nt_submeso = 0.0
    wrho_bt_submeso  = 0.0

    do n=1,num_prog_tracers
       dTdx(n)%field(:,:,:) = 0.0
       dTdy(n)%field(:,:,:) = 0.0
       dTdz(n)%field(:,:,:) = 0.0
    enddo
    dSdx%field(:,:,:) = 0.0
    dSdy%field(:,:,:) = 0.0
    dSdz%field(:,:,:) = 0.0

    psix      = 0.0
    psiy      = 0.0
    psix_horz = 0.0
    psiy_horz = 0.0

    kblt = 0
    hblt = 0.0

    flux_x = 0.0
    flux_y = 0.0
    flux_z = 0.0

    dzwtr = 0.0
    delqc = 0.0

    grid_length(:,:)      = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))
    front_length_inv(:,:) = front_length_const_inv
    buoy_freq_ave(:,:)    = 0.0
    coriolis_param(:,:)   = 2.0*omega_earth*abs(sin(Grd%phit(:,:)))

    dtew(:,:,0) = Grd%dtw(:,:)
    dtew(:,:,1) = Grd%dte(:,:)
    dtns(:,:,0) = Grd%dts(:,:)
    dtns(:,:,1) = Grd%dtn(:,:)

    dtwedyt(:,:,:) = 0.0
    dtwedyt(:,:,0) = Grd%dte(:,:)*Grd%dyt(:,:)
    do i=isc-1,iec
       dtwedyt(i,:,1) = Grd%dtw(i+1,:)*Grd%dyt(i+1,:)
    enddo

    dxtdtsn(:,:,:) = 0.0
    dxtdtsn(:,:,0) = Grd%dxt(:,:)*Grd%dtn(:,:)
    do j=jsc-1,jec
       dxtdtsn(:,j,1) = Grd%dxt(:,j+1)*Grd%dts(:,j+1)
    enddo

    do j=jsd,jed
       do i=isd,ied
          time_factor(i,j) = 1.0/(sqrt( time_constant2_r + coriolis_param(i,j)**2 ))
       enddo
    enddo


    ! for the sweby advection scheme 
    if(submeso_advect_flux .and. submeso_advect_sweby) then 
        allocate(tracer_mdfl_all(num_prog_tracers))
        allocate(tmask_mdfl (isc-2:iec+2,jsc-2:jec+2,nk))    
        allocate(tracer_mdfl(isc-2:iec+2,jsc-2:jec+2,nk))    
        do n=1,num_prog_tracers
           allocate (tracer_mdfl_all(n)%field(isc-2:iec+2,jsc-2:jec+2,nk))
        enddo
        tmask_mdfl  = 0.0
        tracer_mdfl = 0.0 
        do n=1,num_prog_tracers
           tracer_mdfl_all(n)%field(:,:,:) = 0.0
        enddo
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 tmask_mdfl(i,j,k) = Grd%tmask(i,j,k)
              enddo
           enddo
        enddo
        call mpp_update_domains(tmask_mdfl,Dom_mdfl%domain2d)  
    endif


    ! diagnostics 
    id_kblt_submeso = register_diag_field ('ocean_model', 'kblt_submeso', &
         Grid%tracer_axes(1:2), Time%model_time,                          &
         'Number of k-levels in boundary layer for submesoscale closure', &
         'dimensionless', missing_value=missing_value, range=(/-1.0,1e6/))
    id_hblt_submeso = register_diag_field ('ocean_model', 'hblt_submeso', &
         Grid%tracer_axes(1:2), Time%model_time,                          &
         'Boundary layer depth used for submesoscale closure',            &
         'metre', missing_value=missing_value, range=(/-1.0,1e6/))
    id_front_length_submeso = register_diag_field ('ocean_model', 'front_length_submeso', &
         Grid%tracer_axes(1:2), Time%model_time,                                          &
         'Front length used for submesoscale closure',                                    &
         'metre', missing_value=missing_value, range=(/-1.0,1e6/))
    id_buoy_freq_ave_submeso = register_diag_field ('ocean_model', 'buoy_freq_ave_submeso',&
         Grid%tracer_axes(1:2), Time%model_time,                                           &
         'Buoyancy frequency averaged over depth of mixed layer for submesoscale closure', &
         '1/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_mu_submeso = register_diag_field ('ocean_model', 'mu_submeso',   &
         Grd%tracer_axes(1:3), Time%model_time,                         &
         'vertical structure function for submesoscale streamfunction', &
         'dimensionless', missing_value=missing_value, range=(/-1.e2,1e2/))
    id_dmu_submeso = register_diag_field ('ocean_model', 'dmu_submeso',   &
         Grd%tracer_axes(1:3), Time%model_time,                                                &
         'vertical derivative of vertical structure function for submesoscale streamfunction', &
         '1/metre', missing_value=missing_value, range=(/-1.e6,1e6/))
    id_psix_submeso = register_diag_field ('ocean_model', 'psix_submeso', &
         Grd%tracer_axes_flux_y(1:3), Time%model_time,                    &
         'i-comp of submesoscale streamfunction',                         &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))
    id_psiy_submeso = register_diag_field ('ocean_model', 'psiy_submeso', &
         Grd%tracer_axes_flux_x(1:3), Time%model_time,                    &
         'j-comp of submesoscale streamfunction',                         &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))
    id_tx_trans_submeso = register_diag_field ('ocean_model', 'tx_trans_submeso', &
         Grd%tracer_axes_flux_x(1:3), Time%model_time,                            &
         'T-cell mass i-transport from submesoscale param',                       &
         trim(transport_dims), missing_value=missing_value, range=(/-1.e20,1e20/))
    id_ty_trans_submeso = register_diag_field ('ocean_model', 'ty_trans_submeso', &
         Grd%tracer_axes_flux_y(1:3), Time%model_time,                            &
         'T-cell mass j-transport from submesoscale param',                       &
         trim(transport_dims), missing_value=missing_value, range=(/-1.e20,1e20/))
    id_tz_trans_submeso = register_diag_field ('ocean_model', 'tz_trans_submeso', &
         Grd%tracer_axes_wt(1:3), Time%model_time,                                &
         'T-cell mass k-transport from submesoscale param',                       &
         trim(transport_dims), missing_value=missing_value, range=(/-1.e20,1e20/))

    id_tx_trans_nrho_submeso = register_diag_field ('ocean_model','tx_trans_nrho_submeso',&
         Dens%neutralrho_axes_flux_x(1:3),Time%model_time,                                &
         'T-cell i-mass transport from submesoscale param on neutral rho'                 &
         ,trim(transport_dims),missing_value=missing_value, range=(/-1e20,1e20/))
    id_ty_trans_nrho_submeso = register_diag_field ('ocean_model','ty_trans_nrho_submeso',&
         Dens%neutralrho_axes_flux_y(1:3),Time%model_time,                                &
         'T-cell j-mass transport from submesoscale param on neutral rho'                 &
         ,trim(transport_dims),missing_value=missing_value, range=(/-1e20,1e20/))
    id_tz_trans_nrho_submeso = register_diag_field ('ocean_model','tz_trans_nrho_submeso',&
         Dens%neutralrho_axes(1:3),Time%model_time,                                       &
         'T-cell k-mass transport from submesoscale param on neutral rho'                 &
         ,trim(transport_dims),missing_value=missing_value, range=(/-1e20,1e20/))

    id_u_et_submeso = register_diag_field ('ocean_model', 'u_et_submeso', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                     &
       'i-component of submesoscale transport velocity',                  &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    id_v_nt_submeso = register_diag_field ('ocean_model', 'v_nt_submeso', &
        Grd%tracer_axes_flux_y(1:3), Time%model_time,                     &
       'j-component of submesoscale transport velocity',                  &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    id_w_bt_submeso = register_diag_field ('ocean_model', 'w_bt_submeso', &
        Grd%tracer_axes_wt(1:3), Time%model_time,                         &
       'vertical component of submesoscale transport velocity',           &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    id_uhrho_et_submeso = register_diag_field ('ocean_model', 'uhrho_et_submeso', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                             &
       'i-component of submesoscale advective mass transport',                    &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    if(id_uhrho_et_submeso > 0) diag_advect_transport = .true. 

    id_vhrho_nt_submeso = register_diag_field ('ocean_model', 'vhrho_nt_submeso', &
        Grd%tracer_axes_flux_y(1:3), Time%model_time,                             &
       'j-component of submesoscale advective mass transport',                    &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    if(id_vhrho_nt_submeso > 0) diag_advect_transport = .true. 

    id_wrho_bt_submeso = register_diag_field ('ocean_model', 'wrho_bt_submeso', &
        Grd%tracer_axes_wt(1:3), Time%model_time,                               &
       'k-component of submesoscale advective mass transport',                  &
       '(kg/m^3)*m/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    if(id_wrho_bt_submeso > 0) diag_advect_transport = .true. 

    id_tx_trans_submeso_adv = register_diag_field ('ocean_model', 'tx_trans_submeso_adv', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                                     &
       'i-component of submesoscale advective mass transport',                            &
       'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    if(id_tx_trans_submeso_adv > 0) diag_advect_transport = .true. 

    id_ty_trans_submeso_adv = register_diag_field ('ocean_model', 'ty_trans_submeso_adv', &
        Grd%tracer_axes_flux_y(1:3), Time%model_time,                                     &
       'j-component of submesoscale advective mass transport',                            &
       'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    if(id_ty_trans_submeso_adv > 0) diag_advect_transport = .true. 
    
    id_tz_trans_submeso_adv = register_diag_field ('ocean_model', 'tz_trans_submeso_adv', &
        Grd%tracer_axes_wt(1:3), Time%model_time,                                         &
       'k-component of submesoscale advective mass transport',                            &
       'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
    if(id_tz_trans_submeso_adv > 0) diag_advect_transport = .true. 

    if(submeso_diffusion_biharmonic) then 
        id_subdiff_diffusivity = register_diag_field ('ocean_model', 'subdiff_diffusivity', &
             Grd%tracer_axes(1:3), Time%model_time,                                         &
             'horizontal diffusivity used for submesoscale biharmonic mixing scheme',       &
             'm^4/sec', missing_value=missing_value, range=(/-1.e10,1e10/))
    else 
        id_subdiff_diffusivity = register_diag_field ('ocean_model', 'subdiff_diffusivity', &
             Grd%tracer_axes(1:3), Time%model_time,                                         &
             'horizontal diffusivity used for submesoscale laplacian diffusion scheme',     &
             'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))
    endif


    call watermass_diag_init(Time, Dens)


    id_eta_tend_submeso_flx= -1          
    id_eta_tend_submeso_flx= register_diag_field ('ocean_model','eta_tend_submeso_flx',&
         Grd%tracer_axes(1:2), Time%model_time,                                        &
         'non-Bouss steric sea level tendency from submesoscale fluxes', 'm/s',        &
         missing_value=missing_value, range=(/-1e10,1.e10/))
    if(id_eta_tend_submeso_flx > 0) diagnose_eta_tend_submeso_flx=.true.

    id_eta_tend_submeso_flx_glob= -1          
    id_eta_tend_submeso_flx_glob= register_diag_field ('ocean_model', 'eta_tend_submeso_flx_glob',&
         Time%model_time,                                                                         &
         'global mean non-bouss steric sea level tendency from submesoscale fluxes',              &
         'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
    if(id_eta_tend_submeso_flx_glob > 0) diagnose_eta_tend_submeso_flx=.true.

    id_eta_tend_subdiff_flx= -1          
    id_eta_tend_subdiff_flx= register_diag_field ('ocean_model','eta_tend_subdiff_flx',&
         Grd%tracer_axes(1:2), Time%model_time,                                        &
         'non-Bouss steric sea level tendency from subdiff fluxes', 'm/s',             &
         missing_value=missing_value, range=(/-1e10,1.e10/))
    if(id_eta_tend_subdiff_flx > 0) diagnose_eta_tend_subdiff_flx=.true.

    id_eta_tend_subdiff_flx_glob= -1          
    id_eta_tend_subdiff_flx_glob= register_diag_field ('ocean_model', 'eta_tend_subdiff_flx_glob',&
         Time%model_time,                                                                         &
         'global mean non-bouss steric sea level tendency from subdiff fluxes',                   &
         'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
    if(id_eta_tend_subdiff_flx_glob > 0) diagnose_eta_tend_subdiff_flx=.true.

    if(diagnose_eta_tend_submeso_flx .or. diagnose_eta_tend_subdiff_flx) then 
       allocate( flux_x_temp(isd:ied,jsd:jed,nk) )
       allocate( flux_y_temp(isd:ied,jsd:jed,nk) )
       allocate( flux_z_temp(isd:ied,jsd:jed,nk) )
       allocate( flux_x_salt(isd:ied,jsd:jed,nk) )
       allocate( flux_y_salt(isd:ied,jsd:jed,nk) )
       allocate( flux_z_salt(isd:ied,jsd:jed,nk) )
       flux_x_temp(:,:,:) = 0.0
       flux_y_temp(:,:,:) = 0.0
       flux_z_temp(:,:,:) = 0.0
       flux_x_salt(:,:,:) = 0.0
       flux_y_salt(:,:,:) = 0.0
       flux_z_salt(:,:,:) = 0.0
    endif 


    allocate (id_xflux_submeso(num_prog_tracers))
    allocate (id_yflux_submeso(num_prog_tracers))
    allocate (id_zflux_submeso(num_prog_tracers))
    allocate (id_xflux_submeso_int_z(num_prog_tracers))
    allocate (id_yflux_submeso_int_z(num_prog_tracers))
    allocate (id_submeso(num_prog_tracers))
    allocate (id_submeso_on_nrho(num_prog_tracers))

    allocate (id_xflux_subdiff(num_prog_tracers))
    allocate (id_yflux_subdiff(num_prog_tracers))
    allocate (id_xflux_subdiff_int_z(num_prog_tracers))
    allocate (id_yflux_subdiff_int_z(num_prog_tracers))
    allocate (id_subdiff(num_prog_tracers))

    do n=1,num_prog_tracers

       if(n == index_temp) then 

           id_xflux_submeso(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_xflux_submeso',              &
                Grd%tracer_axes_flux_x(1:3), Time%model_time,        &
                'cp*submeso_xflux*dyt*rho_dzt*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_yflux_submeso(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_yflux_submeso',              &
                Grd%tracer_axes_flux_y(1:3), Time%model_time,        &
                'cp*submeso_yflux*dxt*rho_dzt*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_zflux_submeso(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_zflux_submeso',              &
                Grd%tracer_axes_wt(1:3), Time%model_time,            &
                'cp*submeso_zflux*dxt*dyt*rho*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_xflux_submeso_int_z(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_xflux_submeso_int_z',              &
                Grd%tracer_axes_flux_x(1:2), Time%model_time,              &
                'z-integral cp*submeso_xflux*dyt*rho_dzt*temp',            &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_yflux_submeso_int_z(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_yflux_submeso_int_z',              &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,              &
                'z-integral cp*submeso_yflux*dxt*rho_dzt*temp',            &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_submeso(n) = register_diag_field ('ocean_model',           &
                trim(T_prog(n)%name)//'_submeso',                        &              
                Grd%tracer_axes(1:3), Time%model_time,                   &
                'rho*dzt*cp*submesoscale tendency (heating)',            &
                trim(T_prog(n)%flux_units), missing_value=missing_value, &
                range=(/-1.e10,1.e10/))
           id_submeso_on_nrho(n) = register_diag_field ('ocean_model',   &
                trim(T_prog(n)%name)//'_submeso_on_nrho',                &
                Dens%neutralrho_axes(1:3), Time%model_time,              &
                'rho*dzt*cp*submesoscale tendency (heating) binned to neutral density',&
                trim(T_prog(n)%flux_units), missing_value=missing_value, &
                range=(/-1.e20,1.e20/))

           id_xflux_subdiff(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_xflux_subdiff',              &
                Grd%tracer_axes_flux_x(1:3), Time%model_time,        &
                'cp*xflux_subdiff*dyt*rho_dzt*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_yflux_subdiff(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_yflux_subdiff',              &
                Grd%tracer_axes_flux_y(1:3), Time%model_time,        &
                'cp*yflux_subdiff*dxt*rho_dzt*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_xflux_subdiff_int_z(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_xflux_subdiff_int_z',              &
                Grd%tracer_axes_flux_x(1:2), Time%model_time,              &
                'z-integral cp*xflux_subdiff*dyt*rho_dzt*temp',            &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_yflux_subdiff_int_z(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_yflux_subdiff_int_z',              &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,              &
                'z-integral cp*yflux_subdiff*dxt*rho_dzt*temp',            &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_subdiff(n) = register_diag_field ('ocean_model',           &
                trim(T_prog(n)%name)//'_subdiff',                        &              
                Grd%tracer_axes(1:3), Time%model_time,                   &
                'rho*dzt*cp*submesoscale_diff tendency (heating)',       &
                trim(T_prog(n)%flux_units), missing_value=missing_value, &
                range=(/-1.e10,1.e10/))

       else

           id_xflux_submeso(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_xflux_submeso',                       &
                Grd%tracer_axes_flux_x(1:3), Time%model_time,                 &
                'submeso_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
                'kg/sec', missing_value=missing_value,                        &
                range=(/-1.e18,1.e18/))
           id_yflux_submeso(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_yflux_submeso',                       &
                Grd%tracer_axes_flux_y(1:3), Time%model_time,                 &
                'submeso_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
                'kg/sec', missing_value=missing_value,                        &
                range=(/-1.e18,1.e18/))
           id_zflux_submeso(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_zflux_submeso',                       &
                Grd%tracer_axes_wt(1:3), Time%model_time,                     &
                'submeso_yflux*dxt*dyt*rho*tracer for'//trim(T_prog(n)%name), &
                'kg/sec', missing_value=missing_value,                        &
                range=(/-1.e18,1.e18/))
           id_xflux_submeso_int_z(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_xflux_submeso_int_z',                           &
                Grd%tracer_axes_flux_x(1:2), Time%model_time,                           &
                'z-integral submeso_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
                'kg/sec', missing_value=missing_value,                                  &
                range=(/-1.e18,1.e18/))
           id_yflux_submeso_int_z(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_yflux_submeso_int_z',                           &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,                           &
                'z-integral submeso_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
                'kg/sec', missing_value=missing_value,                                  &
                range=(/-1.e18,1.e18/))
           id_submeso(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_submeso',                           &
                Grd%tracer_axes(1:3), Time%model_time,                      &
                'rho*dzt*submesoscale tendency for '//trim(T_prog(n)%name), &
                trim(T_prog(n)%flux_units), missing_value=missing_value,    &
                range=(/-1.e10,1.e10/))
           id_submeso_on_nrho(n) = register_diag_field ('ocean_model',      &
                trim(T_prog(n)%name)//'_submeso_on_nrho',                   &
                Dens%neutralrho_axes(1:3), Time%model_time,                 &
                'rho*dzt*submesoscale tendency for '//trim(T_prog(n)%name)//' binned to neutral density', &
                trim(T_prog(n)%flux_units), missing_value=missing_value,    &
                range=(/-1.e20,1.e20/))


           id_xflux_subdiff(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_xflux_subdiff',                       &
                Grd%tracer_axes_flux_x(1:3), Time%model_time,                 &
                'xflux_subdiff*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
                'kg/sec', missing_value=missing_value,                        &
                range=(/-1.e18,1.e18/))
           id_yflux_subdiff(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_yflux_subdiff',                       &
                Grd%tracer_axes_flux_y(1:3), Time%model_time,                 &
                'yflux_subdiff*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
                'kg/sec', missing_value=missing_value,                        &
                range=(/-1.e18,1.e18/))
           id_xflux_subdiff_int_z(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_xflux_subdiff_int_z',                           &
                Grd%tracer_axes_flux_x(1:2), Time%model_time,                           &
                'z-integral xflux_subdiff*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
                'kg/sec', missing_value=missing_value,                                  &
                range=(/-1.e18,1.e18/))
           id_yflux_subdiff_int_z(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_yflux_subdiff_int_z',                           &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,                           &
                'z-integral yflux_subdiff*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
                'kg/sec', missing_value=missing_value,                                  &
                range=(/-1.e18,1.e18/))
           id_subdiff(n) = register_diag_field ('ocean_model',                   &
                trim(T_prog(n)%name)//'_subdiff',                                &
                Grd%tracer_axes(1:3), Time%model_time,                           &
                'rho*dzt*submesoscale_diff tendency for '//trim(T_prog(n)%name), &
                trim(T_prog(n)%flux_units), missing_value=missing_value,         &
                range=(/-1.e10,1.e10/))

       endif

    enddo


end subroutine ocean_submesoscale_init
! </SUBROUTINE> NAME="ocean_submesoscale_init"


!#######################################################################
! <SUBROUTINE NAME="submeso_restrat">
!
! <DESCRIPTION>
! This routine computes a thickness and density weighted time tendency
! for each tracer, arising from the effects of parameterized 
! submesoscale eddies acting in the surface mixed layer.  
! </DESCRIPTION>
!
  subroutine submeso_restrat(Time, Thickness, Dens, T_prog, surf_blthick)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(in)    :: surf_blthick

  integer :: i,j,k
  integer :: tau, taum1

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_submesoscale_mode (ocean_submeso): needs initialization')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1 

  ! time dependent delqc geometric factor 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           delqc(i,j,k,0) = Grd%fracdz(k,0)*Thickness%rho_dzt(i,j,k,tau)
           delqc(i,j,k,1) = Grd%fracdz(k,1)*Thickness%rho_dzt(i,j,k,tau)
        enddo
     enddo
  enddo

  ! time dependent inverse dzwt 
  do k=0,nk
     do j=jsd,jed
        do i=isd,ied
           dzwtr(i,j,k) = 1.0/Thickness%dzwt(i,j,k) 
        enddo
     enddo
  enddo

  call compute_bldepth(Time, Thickness, Dens, T_prog, surf_blthick)
  call tracer_derivs(taum1, T_prog) 
  call salinity_derivs(taum1, Dens) 

  if(use_psi_legacy) then 
     call compute_psi_legacy(Time, Dens, Thickness)
  else 
     call compute_psi(Time, Dens, Thickness)
  endif 

  if(diag_advect_transport .or. submeso_advect_flux) then 
     call compute_advect_transport(Time, Dens, Thickness)
  endif 

  ! compute tracer flux components and their convergence
  if(submeso_skew_flux) then 
      call compute_submeso_skewsion(Thickness, Dens, Time, T_prog)
  elseif(submeso_advect_flux) then
      if(submeso_advect_upwind) then  
         call compute_submeso_upwind(Time, Dens, T_prog)
      elseif(submeso_advect_sweby) then 
         call compute_submeso_sweby(Thickness, Time, Dens, T_prog)
      endif 
  endif
  call watermass_diag(Time, T_prog, Dens)


  if(submeso_diffusion) then 
      call compute_submeso_diffusion(Thickness, Dens, Time, T_prog)
  endif 
  call watermass_diag_diffusion(Time, T_prog, Dens)


end subroutine submeso_restrat
! </SUBROUTINE> NAME="submeso_restrat"



!#######################################################################
! <SUBROUTINE NAME="compute_bldepth">
!
! <DESCRIPTION>
! Compute the boundary layer depth and kblt.
! </DESCRIPTION>
!
subroutine compute_bldepth(Time, Thickness, Dens, T_prog, surf_blthick) 

  type(ocean_time_type),       intent(in) :: Time
  type(ocean_thickness_type),  intent(in) :: Thickness
  type(ocean_density_type),    intent(in) :: Dens
  type(ocean_prog_tracer_type),intent(in) :: T_prog(:)
  real, dimension(isd:,jsd:),  intent(in) :: surf_blthick

  integer :: i, j, k, tau
  integer :: num_smooth
  real    :: mld_thickness
  real    :: active_cells 

  tau  = Time%tau 
  hblt = 0.0

  if(use_hblt_constant) then 

      do j=jsc,jec
         do i=isc,iec
            hblt(i,j) = Grd%tmask(i,j,1)*min(constant_hblt, Grd%ht(i,j))
         enddo
      enddo

  elseif(use_hblt_equal_mld) then 

      call calc_mixed_layer_depth(Thickness,                &
           Dens%rho_salinity(isd:ied,jsd:jed,:,tau),        &
           T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau), &
           Dens%rho(isd:ied,jsd:jed,:,tau),                 &
           Dens%pressure_at_depth(isd:ied,jsd:jed,:),       &
           hblt, smooth_mld_input=.false.)

      do j=jsc,jec
         do i=isc,iec
            hblt(i,j) = Grd%tmask(i,j,1)*min(hblt(i,j), Grd%ht(i,j))
         enddo
      enddo

  else 

      do j=jsc,jec
         do i=isc,iec
            hblt(i,j) = Grd%tmask(i,j,1)*min(surf_blthick(i,j), Grd%ht(i,j))
         enddo
      enddo

  endif

  if(smooth_hblt) then 
      do num_smooth=1,smooth_hblt_num

         wrk1_2d(:,:) = 0.0
         call mpp_update_domains(hblt(:,:), Dom%domain2d) 
         k=1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%tmask(i,j,k)==1.0) then 
                   active_cells = 4.0       +&
                        Grd%tmask(i-1,j,k)  +&
                        Grd%tmask(i+1,j,k)  +&
                        Grd%tmask(i,j-1,k)  +&
                        Grd%tmask(i,j+1,k)
                   if (active_cells > 4.0) then
                       wrk1_2d(i,j) = & 
                            (4.0*hblt(i,j) +&
                            hblt(i-1,j)    +&
                            hblt(i+1,j)    +&
                            hblt(i,j-1)    +&
                            hblt(i,j+1)) / active_cells
                   else
                       wrk1_2d(i,j) = hblt(i,j)
                   endif
               endif
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               hblt(i,j) = wrk1_2d(i,j)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
            enddo
         enddo

      enddo  ! enddo for number of smoothing iterations 
  endif  ! endif for smooth_hblt

  ! set floor to hblt 
  do j=jsc,jec
     do i=isc,iec
        hblt(i,j) = Grd%tmask(i,j,1)*max(minimum_hblt, hblt(i,j))
     enddo
  enddo

  ! halo values needed 
  call mpp_update_domains(hblt(:,:), Dom%domain2d) 


  ! k-index at bottom of hblt
  ! also move the hblt to equal depth_zwt(kblt); this is 
  ! necessary to get a full extent of the streamfunction 
  ! contained in the boundary layer. 
  kblt=0
  do j=jsd,jed
     do i=isd,ied
        if(Grd%kmt(i,j) > 1) then  
            kloop: do k=1,nk
               if(Thickness%depth_zwt(i,j,k) >= hblt(i,j)) then
                   kblt(i,j) = k
                   hblt(i,j) = Thickness%depth_zwt(i,j,k)
                   exit kloop 
               endif
            enddo kloop
        endif
     enddo
  enddo

  ! inverse front length = f/(<N> H), with H=hblt and <N> ave buoyancy freq over hblt
  if(front_length_deform_radius) then 
      do j=jsd,jed
         do i=isd,ied
            buoy_freq_ave(i,j)    = 0.0
            mld_thickness         = epsln
            front_length_inv(i,j) = front_length_const_inv 
            if(kblt(i,j) >= min_kblt) then  
                buoy_freq_ave(i,j)= epsln
                mld_thickness     = epsln
                do k=1,kblt(i,j)
                   mld_thickness      = mld_thickness      + Thickness%dzt(i,j,k)  
                   buoy_freq_ave(i,j) = buoy_freq_ave(i,j) - grav_rho0r*Thickness%dzt(i,j,k)*Dens%drhodz_zt(i,j,k)  
                enddo
                buoy_freq_ave(i,j)    = sqrt(abs(buoy_freq_ave(i,j))/mld_thickness) 
                front_length_inv(i,j) = min(front_length_const_inv, coriolis_param(i,j)/(epsln+mld_thickness*buoy_freq_ave(i,j))) 
            endif 
         enddo
      enddo
  endif


  ! diagnostics  
  if (id_front_length_submeso > 0) then 
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = min(front_length_max, 1.0/(front_length_inv(i,j)+epsln))
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_front_length_submeso, wrk1_2d(:,:))
  endif 
  call diagnose_2d(Time, Grd, id_buoy_freq_ave_submeso, buoy_freq_ave(:,:))
  call diagnose_2d(Time, Grd, id_hblt_submeso, hblt(:,:))
  if (id_kblt_submeso > 0) then 
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = kblt(i,j)
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_kblt_submeso, wrk1_2d(:,:))
  endif 


end subroutine compute_bldepth
! </SUBROUTINE> NAME="compute_bldepth"


!#######################################################################
! <SUBROUTINE NAME="tracer_derivs">
!
! <DESCRIPTION>
! Compute the tracer derivatives, with the  
! lateral derivatives computed along constant k-level.
! </DESCRIPTION>
!
subroutine tracer_derivs(taum1, T_prog)

  integer,                      intent(in) :: taum1
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  integer :: i, j, k, n
  integer :: kp1
  real    :: tmaski, tmaskj

  do n=1,num_prog_tracers

     dTdx(n)%field(:,:,:) = 0.0
     dTdy(n)%field(:,:,:) = 0.0
     dTdz(n)%field(:,:,:) = 0.0

     do j=jsc-1,jec
        do i=isc-1,iec

           if(kblt(i,j) >= min_kblt) then  
               do k=1,kblt(i,j)
                  tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
                  tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  dTdx(n)%field(i,j,k) = (T_prog(n)%field(i+1,j,k,taum1)-T_prog(n)%field(i,j,k,taum1)) &
                                          *Grd%dxter(i,j)*tmaski
                  dTdy(n)%field(i,j,k) = (T_prog(n)%field(i,j+1,k,taum1)-T_prog(n)%field(i,j,k,taum1)) &
                                          *Grd%dytnr(i,j)*tmaskj
               enddo
           endif

        enddo
     enddo

     do j=jsd,jed
        do i=isd,ied
           if(kblt(i,j) >= min_kblt) then  
               do k=1,kblt(i,j)
                  kp1 = min(k+1,nk)
                  dTdz(n)%field(i,j,k) = (T_prog(n)%field(i,j,k,taum1)-T_prog(n)%field(i,j,kp1,taum1)) &
                                          *Grd%tmask(i,j,kp1)*dzwtr(i,j,k)
               enddo
           endif
        enddo
     enddo

  enddo

end subroutine tracer_derivs
! </SUBROUTINE> NAME="tracer_derivs"


!#######################################################################
! <SUBROUTINE NAME="salinity_derivs">
!
! <DESCRIPTION>
! Compute the density-salinity derivatives, with lateral 
! derivative computed on constant k-level. 
! </DESCRIPTION>
!
subroutine salinity_derivs(taum1, Dens)

  integer,                      intent(in) :: taum1
  type(ocean_density_type),     intent(in) :: Dens

  integer :: i, j, k
  integer :: kp1
  real    :: tmaski, tmaskj

  dSdx%field(:,:,:) = 0.0
  dSdy%field(:,:,:) = 0.0
  dSdz%field(:,:,:) = 0.0

  do j=jsc-1,jec
     do i=isc-1,iec

        if(kblt(i,j) >= min_kblt) then  
            do k=1,kblt(i,j)
               tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
               dSdx%field(i,j,k) = (Dens%rho_salinity(i+1,j,k,taum1)-Dens%rho_salinity(i,j,k,taum1)) &
                                       *Grd%dxter(i,j)*tmaski
               dSdy%field(i,j,k) = (Dens%rho_salinity(i,j+1,k,taum1)-Dens%rho_salinity(i,j,k,taum1)) &
                                       *Grd%dytnr(i,j)*tmaskj
            enddo
        endif

     enddo
  enddo

  do j=jsd,jed
     do i=isd,ied
        if(kblt(i,j) >= min_kblt) then  
            do k=1,kblt(i,j)
               kp1 = min(k+1,nk)
               dSdz%field(i,j,k) = (Dens%rho_salinity(i,j,k,taum1)-Dens%rho_salinity(i,j,kp1,taum1)) &
                                       *Grd%tmask(i,j,kp1)*dzwtr(i,j,k)
            enddo
        endif
     enddo
  enddo


end subroutine salinity_derivs
! </SUBROUTINE> NAME="salinity_derivs"


!#######################################################################
! <SUBROUTINE NAME="compute_psi">
!
! <DESCRIPTION>
! Compute the vector streamfunction from parameterized 
! submesoscale restratification. 
!
! Units of psi are m^2/sec
!
! psix is defined on north face of tracer cell for jq=0,1.
! psiy is defined on east  face of tracer cell for ip=0,1.
!
! NOTE: the mpp updates for psix and psiy are treated as a 
! scalar, whereas they are actually components to a pseudo-vector.
! Some further thought is required for the tripolar grid. We ignore
! this detail in the present implementation.  
!
! </DESCRIPTION>
!
subroutine compute_psi(Time, Dens, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness

  real    :: gradx, grady
  real    :: tmaski, tmaskj
  real    :: gradxrho(0:1), gradyrho(0:1)
  real    :: factor, coefficient
  real    :: active_cells 
  real    :: max_psi
  real    :: abs_psi
  real    :: rescale_psi
  real    :: hblt_r
  integer :: i, j, k
  integer :: ip, jq
  integer :: tau
  integer :: num_smooth

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau = Time%tau

  ! compute depth-independent portion of vector streamfunction (m^2/sec)
  psix_horz = 0.0
  psiy_horz = 0.0
  do j=jsc,jec
     do i=isc,iec

        if(kblt(i,j) >= min_kblt) then  

            gradxrho(:) = 0.0
            gradyrho(:) = 0.0

            do k=1,kblt(i,j)
               do ip=0,1 
                  jq=ip
                  gradx        =   Dens%drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k) &
                                  +Dens%drhodS(i+ip,j,k)*dSdx%field(i,j,k) 
                  grady        =   Dens%drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k) &
                                  +Dens%drhodS(i,j+jq,k)*dSdy%field(i,j,k)
                  gradxrho(ip) = gradxrho(ip) + gradx*Thickness%dzt(i,j,k)  
                  gradyrho(jq) = gradyrho(jq) + grady*Thickness%dzt(i,j,k)  
               enddo
            enddo

            ! normalization by depth of boundary layer 
            hblt_r = 1.0/hblt(i,j)
            do ip=0,1 
               jq=ip
               gradxrho(ip) = gradxrho(ip)*hblt_r
               gradyrho(jq) = gradyrho(jq)*hblt_r
            enddo

            ! coefficient has units m^6/(sec*kg)
            coefficient = time_factor(i,j)*ce_grav_rho0r*grid_length(i,j)*front_length_inv(i,j)*hblt(i,j)**2

            k=1
            tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
            tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
            do ip=0,1 
               jq=ip
               psix_horz(i,j,jq) = -coefficient*gradyrho(jq)*tmaskj
               psiy_horz(i,j,ip) =  coefficient*gradxrho(ip)*tmaski
            enddo

        endif ! endif for kblt(i,j) >= min_kblt

     enddo  ! i-loop
  enddo     ! j-loop


  ! rescale psi if its maximum is larger than limited value. 
  ! to be conservative, base the maximum on the upper ocean thickness.
  ! this option is useful to avoid instabilities in presence of 
  ! strong lateral density gradients.
  if(limit_psi) then 
      do j=jsc,jec
         do i=isc,iec
            if(Grd%tmask(i,j,1) > 0.0) then 

               max_psi = limit_psi_velocity_scale*Thickness%dzt(i,j,1)   

               do ip=0,1 
                  jq=ip

                  abs_psi = abs(psix_horz(i,j,jq))
                  if(abs_psi > max_psi) then 
                     rescale_psi = max_psi/(epsln+abs_psi)
                     psix_horz(i,j,jq) = rescale_psi*psix_horz(i,j,jq)
                  endif 

                  abs_psi = abs(psiy_horz(i,j,ip))
                  if(abs_psi > max_psi) then 
                     rescale_psi = max_psi/(epsln+abs_psi)
                     psiy_horz(i,j,ip) = rescale_psi*psiy_horz(i,j,ip)
                  endif 

               enddo

            endif ! endif for tmask > 0
         enddo    ! enddo for i
      enddo       ! enddo for j
  endif

  ! spatially smooth psi_horz
  if(smooth_psi) then 
      do num_smooth=1,smooth_psi_num

         call mpp_update_domains(psix_horz(:,:,:), Dom%domain2d) 
         call mpp_update_domains(psiy_horz(:,:,:), Dom%domain2d) 

         do ip=0,1
            jq=ip
            do k=1,1     ! to keep loop structure as in other 3d smoothers 

               wrk1_2d(:,:) = 0.0
               wrk2_2d(:,:) = 0.0
               do j=jsc,jec
                  do i=isc,iec

                     if(Grd%tmask(i,j,k)==1.0) then 

                         active_cells = 4.0       +&
                              Grd%tmask(i-1,j,k)  +&
                              Grd%tmask(i+1,j,k)  +&
                              Grd%tmask(i,j-1,k)  +&
                              Grd%tmask(i,j+1,k)

                         if (active_cells > 4.0) then
                             wrk1_2d(i,j) =               &  
                                  (4.0*psix_horz(i,j,jq) +&
                                  psix_horz(i-1,j,jq)    +&
                                  psix_horz(i+1,j,jq)    +&
                                  psix_horz(i,j-1,jq)    +&
                                  psix_horz(i,j+1,jq)) / active_cells
                             wrk2_2d(i,j) =               &
                                  (4.0*psiy_horz(i,j,ip) +&
                                  psiy_horz(i-1,j,ip)    +&
                                  psiy_horz(i+1,j,ip)    +&
                                  psiy_horz(i,j-1,ip)    +&
                                  psiy_horz(i,j+1,ip)) / active_cells

                         else
                             wrk1_2d(i,j) = psix_horz(i,j,jq)
                             wrk2_2d(i,j) = psiy_horz(i,j,ip)
                         endif

                     endif

                  enddo
               enddo

               do j=jsc,jec
                  do i=isc,iec
                     psix_horz(i,j,jq) = wrk1_2d(i,j)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                     psiy_horz(i,j,ip) = wrk2_2d(i,j)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
                  enddo
               enddo

            enddo  ! k=1,1 loop
         enddo     ! ip=0,1 loop

         call mpp_update_domains(psix_horz(:,:,:), Dom%domain2d) 
         call mpp_update_domains(psiy_horz(:,:,:), Dom%domain2d) 

      enddo  ! enddo for number of smoothing iterations 
  endif      ! endif for smooth_psi


  ! compute 3d vector streamfunction components by applying 
  ! dimensionless vertical structure function 0 <= wrk3 <= 1.
  psix = 0.0
  psiy = 0.0
  wrk3 = 0.0
  do j=jsc,jec
      do i=isc,iec
          if(hblt(i,j) > 0.0 .and. kblt(i,j) > min_kblt) then 
              hblt_r = 1.0/hblt(i,j)
              do k=1,kblt(i,j) 
                  factor      = (1.0 - 2.0*Thickness%depth_zt(i,j,k)*hblt_r)**2
                  wrk3(i,j,k) = (1.0 - factor)*(1.0 + fiveover21*factor) 
                  tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
                  tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  do ip=0,1 
                     jq=ip
                     psix(i,j,k,jq) = wrk3(i,j,k)*psix_horz(i,j,jq)*tmaskj
                     psiy(i,j,k,ip) = wrk3(i,j,k)*psiy_horz(i,j,ip)*tmaski
                  enddo
              enddo
          endif
      enddo
  enddo
  call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
  call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 


  ! send diagnostics 
  if (id_psix_submeso > 0) then 
     call diagnose_3d(Time, Grd, id_psix_submeso, onehalf*(psix(:,:,:,0)+psix(:,:,:,1)))
  endif 
  if (id_psiy_submeso > 0) then
     call diagnose_3d(Time, Grd, id_psiy_submeso, onehalf*(psiy(:,:,:,0)+psiy(:,:,:,1)))
  endif 
  call diagnose_3d(Time, Grd, id_mu_submeso, wrk3(:,:,:))


  if(debug_this_module) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_submeso_mod: chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('psix(:,:,:,0)', psix(COMP,:,0)*Grd%tmask(COMP,:))
      call write_chksum_3d('psix(:,:,:,1)', psix(COMP,:,1)*Grd%tmask(COMP,:))
      call write_chksum_3d('psiy(:,:,:,0)', psiy(COMP,:,0)*Grd%tmask(COMP,:))
      call write_chksum_3d('psiy(:,:,:,1)', psiy(COMP,:,1)*Grd%tmask(COMP,:))
  endif


end subroutine compute_psi
! </SUBROUTINE> NAME="compute_psi"



!#######################################################################
! <SUBROUTINE NAME="compute_psi_legacy">
!
! <DESCRIPTION>
! Compute the vector streamfunction 
!
! Units of psi are m^2/sec
!
! If computing skewsion tendency, then need psi at depth_zt. 
! If computing advection tendency, then need psi at depth_zwt. 
!
! Jan2012: This scheme has problems with the limiters and smoothers.
! These problems become particularly egregious when trying to compute
! an advective flux rather than a skew flux.  This routine is 
! retained only for legacy purposes. 
! Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine compute_psi_legacy(Time, Dens, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness
  real    :: gradx, grady
  real    :: tmaski, tmaskj
  real    :: gradxrho(0:1), gradyrho(0:1)
  real    :: factor, coefficient
  real    :: active_cells 
  real    :: max_psi
  real    :: hblt_r
  integer :: i, j, k
  integer :: ip, jq
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau = Time%tau

  wrk1 = 0.0   ! depth 
  wrk2 = 0.0
  wrk3 = 0.0   ! mu_submeso

  ! vector streamfunction components
  psix = 0.0
  psiy = 0.0


  ! for holding the depths 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk1(i,j,k) = Thickness%depth_zt(i,j,k) 
        enddo
     enddo
  enddo

  do j=jsd,jec
     do i=isd,iec

        if(kblt(i,j) >= min_kblt) then  

            gradxrho(:) = 0.0
            gradyrho(:) = 0.0

            do k=1,kblt(i,j)
               do ip=0,1 
                  jq=ip
                  gradx        =   Dens%drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k) &
                                  +Dens%drhodS(i+ip,j,k)*dSdx%field(i,j,k) 
                  grady        =   Dens%drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k) &
                                  +Dens%drhodS(i,j+jq,k)*dSdy%field(i,j,k)
                  gradxrho(ip) = gradxrho(ip) + gradx*Thickness%dzt(i,j,k)  
                  gradyrho(jq) = gradyrho(jq) + grady*Thickness%dzt(i,j,k)  
               enddo

            enddo   ! enddo for k=1,kblt(i,j)

            ! normalization by depth of boundary layer 
            do ip=0,1 
               jq=ip
               gradxrho(ip) = gradxrho(ip)/hblt(i,j)
               gradyrho(jq) = gradyrho(jq)/hblt(i,j)
            enddo

            ! coefficient has units m^6/(sec*kg)
            coefficient = time_factor(i,j)*ce_grav_rho0r*grid_length(i,j)*front_length_inv(i,j)*hblt(i,j)**2
            hblt_r      = 1.0/hblt(i,j)

            ! compute the vector streamfunction (m^2/sec)
            do k=1,kblt(i,j) 

               tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)

               factor      = (1.0 - 2.0*wrk1(i,j,k)*hblt_r)**2
               wrk3(i,j,k) = (1.0 - factor)*(1.0 + fiveover21*factor) 

               do ip=0,1 
                  jq=ip
                  psix(i,j,k,jq) = -coefficient*wrk3(i,j,k)*gradyrho(jq)*tmaskj
                  psiy(i,j,k,ip) =  coefficient*wrk3(i,j,k)*gradxrho(ip)*tmaski
               enddo

            enddo

        endif ! endif for kblt(i,j) >= min_kblt

     enddo
  enddo


  ! limit magnitude of psi to reduce potential for crashes 
  ! Caution: this limiter method is really wrong, as it corrupts the vertical structure
  ! function and leads to spurious jumps in the horizontal eddy induced transport.  
  if(limit_psi) then 
      do j=jsd,jec
         do i=isd,iec
            do k=1,kblt(i,j) 
               max_psi = limit_psi_velocity_scale*Thickness%dzt(i,j,k)
               do ip=0,1 
                  jq=ip
                    psix(i,j,k,jq) = sign(1.0,psix(i,j,k,jq))*min(max_psi,abs(psix(i,j,k,jq))) 
                    psiy(i,j,k,ip) = sign(1.0,psiy(i,j,k,ip))*min(max_psi,abs(psiy(i,j,k,ip))) 
               enddo
            enddo
         enddo
      enddo
  endif


  ! spatially smooth psi
  ! Caution: this smoother method is really wrong, as it corrupts the vertical structure
  ! function and leads to spurious jumps in the horizontal eddy induced transport.  
  if(smooth_psi) then 

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

      do ip=0,1
         jq=ip
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec

                  if(Grd%tmask(i,j,k)==1.0) then 

                      active_cells = 4.0       +&
                           Grd%tmask(i-1,j,k)  +&
                           Grd%tmask(i+1,j,k)  +&
                           Grd%tmask(i,j-1,k)  +&
                           Grd%tmask(i,j+1,k)

                      if (active_cells > 4.0) then
                          wrk1(i,j,k) = & 
                               (4.0*psix(i,j,k,jq) +&
                               psix(i-1,j,k,jq)    +&
                               psix(i+1,j,k,jq)    +&
                               psix(i,j-1,k,jq)    +&
                               psix(i,j+1,k,jq)) / active_cells
                          wrk2(i,j,k) =  &
                               (4.0*psiy(i,j,k,ip) +&
                               psiy(i-1,j,k,ip)    +&
                               psiy(i+1,j,k,ip)    +&
                               psiy(i,j-1,k,ip)    +&
                               psiy(i,j+1,k,ip)) / active_cells

                      else
                          wrk1(i,j,k) = psix(i,j,k,jq)
                          wrk2(i,j,k) = psiy(i,j,k,ip)
                      endif

                  endif

               enddo
            enddo

            do j=jsc,jec
               do i=isc,iec
                  psix(i,j,k,jq) = wrk1(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  psiy(i,j,k,ip) = wrk2(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               enddo
            enddo

         enddo
      enddo

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

  endif


  ! send diagnostics 
  if (id_psix_submeso > 0) then 
     call diagnose_3d(Time, Grd, id_psix_submeso, onehalf*(psix(:,:,:,0)+psix(:,:,:,1)))
  endif 

  if (id_psiy_submeso > 0) then 
     call diagnose_3d(Time, Grd, id_psiy_submeso, onehalf*(psiy(:,:,:,0)+psiy(:,:,:,1)))
  endif 
  call diagnose_3d(Time, Grd, id_mu_submeso, wrk3(:,:,:))


  if(debug_this_module) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_submeso_mod: chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('psix(:,:,:,0)', psix(COMP,:,0)*Grd%tmask(COMP,:))
      call write_chksum_3d('psix(:,:,:,1)', psix(COMP,:,1)*Grd%tmask(COMP,:))
      call write_chksum_3d('psiy(:,:,:,0)', psiy(COMP,:,0)*Grd%tmask(COMP,:))
      call write_chksum_3d('psiy(:,:,:,1)', psiy(COMP,:,1)*Grd%tmask(COMP,:))
  endif


end subroutine compute_psi_legacy
! </SUBROUTINE> NAME="compute_psi_legacy"


!#######################################################################
! <SUBROUTINE NAME="compute_advect_transport">
!
! <DESCRIPTION>
! Compute the mass transport from submeso.  
!
! This routine is a diagnostic routine if skewsion, and 
! part of the calculation of the eddy-induced velocity if 
! advective approach used. 
!
! Comments on the scheme:
! 
! 1/ compute vertical component from convergence of horizontal, just
! as for the vertical velocity component for the Eulerian transport.
! 
! 2/ wrho_bt_submeso(:,:,k=0) = 0.0 by definition
!
! 3/ expand the BDX_ET and BDY_NT operators for efficiency.  
!
! 4/ mask to zero those regions where the horizontal divergence 
! vanishes, as these are regions beneath the submeso boundary 
! layer.  Base the mask on horz divergence rather than kblt(i,j), 
! since any smoothing performed to uhrho_et_submeso and vhrho_nt_submeso
! will modify the region of nonzero submesoscale advection so that 
! it reaches potentially to below kblt.  
!
! </DESCRIPTION>
!
subroutine compute_advect_transport(Time, Dens, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness

  real    :: active_cells
  real    :: deriv_x, deriv_y
  real    :: divergence_mask, divergence
  real    :: tmask_bdy
  real    :: term1, term2, hblt_r
  integer :: i, j, k, kp1
  integer :: tau
  integer :: num_smooth

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau = Time%tau

  ! compute rho*Psi according to vertical coordinate class.   
  ! compute analytic vertical derivative of structure function
  ! and multipy by psix_horz and psiy_horz.
  ! ignore vertical derivative of density for pressure vertical coordinates. 

  uhrho_et_submeso = 0.0
  vhrho_nt_submeso = 0.0
  wrk1(:,:,:) = 0.0   ! dmu/dz = vertical derivative of structure function 
  wrk2(:,:,:) = 0.0   ! to hold the rho_dzt field
  do j=jsc,jec
     do i=isc,iec
        if(kblt(i,j) > min_kblt) then 
            hblt_r = Grd%tmask(i,j,1)/(epsln+hblt(i,j))
            do k=1,kblt(i,j)
               wrk2(i,j,k) = Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
               term1       = (8.0/21.0)*hblt_r*(1.0 - 2.0*Thickness%depth_zt(i,j,k)*hblt_r)
               term2       = 8.0 + 5.0*(1.0-2.0*Thickness%depth_zt(i,j,k)*hblt_r)**2
               wrk1(i,j,k) = -term1*term2  
               uhrho_et_submeso(i,j,k) = -wrk2(i,j,k)*wrk1(i,j,k)*onehalf*(psiy_horz(i,j,1)+psiy_horz(i,j,0))
               vhrho_nt_submeso(i,j,k) =  wrk2(i,j,k)*wrk1(i,j,k)*onehalf*(psix_horz(i,j,1)+psix_horz(i,j,0))
            enddo
        endif  
     enddo
  enddo

  ! this option is typically not needed if limit_psi is enabled.
  if(submeso_advect_limit) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               uhrho_et_submeso(i,j,k) = sign(1.0,uhrho_et_submeso(i,j,k)) &
               *min(limit_psi_velocity_scale*Thickness%rho_dzt(i,j,k,tau),abs(uhrho_et_submeso(i,j,k))) 
               vhrho_nt_submeso(i,j,k) = sign(1.0,vhrho_nt_submeso(i,j,k)) &
               *min(limit_psi_velocity_scale*Thickness%rho_dzt(i,j,k,tau),abs(vhrho_nt_submeso(i,j,k))) 
            enddo
         enddo
      enddo
  endif


  ! update domains as per C-grid fluxes 
  call mpp_update_domains(uhrho_et_submeso(:,:,:), vhrho_nt_submeso(:,:,:),&
                          Dom_flux_sub%domain2d,gridtype=CGRID_NE)


  ! apply horizontal smoother to the transport.
  ! this option is typically not be needed if smooth_psi option enabled. 
  if(smooth_advect_transport) then 
      do num_smooth = 1,smooth_advect_transport_num

         wrk1_v(:,:,:,:) = 0.0

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec

                  if(Grd%tmask(i,j,k)==1.0) then 

                      active_cells = 4.0       +&
                           Grd%tmask(i-1,j,k)  +&
                           Grd%tmask(i+1,j,k)  +&
                           Grd%tmask(i,j-1,k)  +&
                           Grd%tmask(i,j+1,k)

                      if (active_cells > 4.0) then
                          wrk1_v(i,j,k,1) = & 
                          (4.0*uhrho_et_submeso(i,j,k)      +&
                               uhrho_et_submeso(i-1,j,k)    +&
                               uhrho_et_submeso(i+1,j,k)    +&
                               uhrho_et_submeso(i,j-1,k)    +&
                               uhrho_et_submeso(i,j+1,k)) / active_cells
                          wrk1_v(i,j,k,2) =  &
                          (4.0*vhrho_nt_submeso(i,j,k)      +&
                               vhrho_nt_submeso(i-1,j,k)    +&
                               vhrho_nt_submeso(i+1,j,k)    +&
                               vhrho_nt_submeso(i,j-1,k)    +&
                               vhrho_nt_submeso(i,j+1,k)) / active_cells

                      else
                          wrk1_v(i,j,k,1) = uhrho_et_submeso(i,j,k)
                          wrk1_v(i,j,k,2) = vhrho_nt_submeso(i,j,k)
                      endif

                  endif

               enddo
            enddo

            do j=jsc,jec
               do i=isc,iec
                  uhrho_et_submeso(i,j,k) = wrk1_v(i,j,k,1)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
                  vhrho_nt_submeso(i,j,k) = wrk1_v(i,j,k,2)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
               enddo
            enddo

         enddo

         call mpp_update_domains(uhrho_et_submeso(:,:,:), vhrho_nt_submeso(:,:,:),&
                                 Dom_flux_sub%domain2d,gridtype=CGRID_NE)

      enddo  ! enddo for number of smoothing iterations 
  endif      ! endif for smooth_advect_transport

  ! set transports to zero next to any land points, including corners.
  ! this option aims to remove what can be some "noisy" features from 
  ! the advective tendencies the occurs near land.  
  if(submeso_advect_zero_bdy) then 
      do k=1,nk
         kp1 = min(k+1,nk)
         do j=jsc,jec
            do i=isc,iec
               tmask_bdy = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)*Grd%tmask(i-1,j,k)         &
                                           *Grd%tmask(i,j+1,k)*Grd%tmask(i,j-1,k)         &
                                           *Grd%tmask(i+1,j+1,k)*Grd%tmask(i-1,j-1,k)     &
                          *Grd%tmask(i,j,kp1)*Grd%tmask(i+1,j,kp1)*Grd%tmask(i-1,j,kp1)   &
                                           *Grd%tmask(i,j+1,kp1)*Grd%tmask(i,j-1,kp1)     &
                                           *Grd%tmask(i+1,j+1,kp1)*Grd%tmask(i-1,j-1,kp1)
               uhrho_et_submeso(i,j,k) = tmask_bdy*uhrho_et_submeso(i,j,k) 
               vhrho_nt_submeso(i,j,k) = tmask_bdy*vhrho_nt_submeso(i,j,k) 
            enddo
         enddo
      enddo
      call mpp_update_domains(uhrho_et_submeso(:,:,:), vhrho_nt_submeso(:,:,:),&
                              Dom_flux_sub%domain2d,gridtype=CGRID_NE)
  endif


  ! diagnose vertical transport. read comments given at
  ! top of this subroutine for details.
  wrho_bt_submeso(:,:,:) = 0.0
  do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            deriv_x = (Grd%dyte(i,j)*uhrho_et_submeso(i,j,k) - Grd%dyte(i-1,j)*uhrho_et_submeso(i-1,j,k))*Grd%datr(i,j)
            deriv_y = (Grd%dxtn(i,j)*vhrho_nt_submeso(i,j,k) - Grd%dxtn(i,j-1)*vhrho_nt_submeso(i,j-1,k))*Grd%datr(i,j)
            divergence = deriv_x + deriv_y 
            divergence_mask = 0.0
            if(divergence /= 0.0) divergence_mask = 1.0
            wrho_bt_submeso(i,j,k) = divergence_mask*Grd%tmask(i,j,k)*(wrho_bt_submeso(i,j,k-1) + divergence)
         enddo
     enddo
  enddo


  ! send diagnostics 

  call diagnose_3d(Time, Grd, id_dmu_submeso, wrk1(:,:,:))
  call diagnose_3d(Time, Grd, id_uhrho_et_submeso, uhrho_et_submeso(:,:,:))
  call diagnose_3d(Time, Grd, id_vhrho_nt_submeso, vhrho_nt_submeso(:,:,:))
  call diagnose_3d(Time, Grd, id_wrho_bt_submeso, wrho_bt_submeso(:,:,1:nk))

  if (id_u_et_submeso > 0) then
      wrk1 = 0.0 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = Grd%tmask(i,j,k)*uhrho_et_submeso(i,j,k) &
                             /(epsln+Thickness%rho_dzt(i,j,k,tau))
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_u_et_submeso, wrk1(:,:,:))
  endif
  if (id_v_nt_submeso > 0) then
      wrk1 = 0.0 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = Grd%tmask(i,j,k)*vhrho_nt_submeso(i,j,k) &
                             /(epsln+Thickness%rho_dzt(i,j,k,tau))
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_v_nt_submeso, wrk1(:,:,:))
  endif
  if (id_w_bt_submeso > 0) then
      wrk1 = 0.0 
      do k=1,nk
         kp1 = min(k+1,nk)
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)*Thickness%dzt(i,j,k)*wrho_bt_submeso(i,j,k) &
                             /(epsln+Thickness%rho_dzt(i,j,k,tau))
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_w_bt_submeso, wrk1(:,:,:))
  endif

  
  !  for submeso_advect_flux, tx_trans_submeso_adv = tx_trans_submeso, as for other components.
  !  for submeso_shew_flux, they are distinct.  it is useful when running with skew approach to
  !  diagnose the advective mass transports from submeso. GM has an analogous diagnostic.  
  if (id_tx_trans_submeso_adv > 0 .or. id_ty_trans_submeso_adv > 0 .or. id_tz_trans_submeso_adv > 0) then 

      wrk1_v(:,:,:,:) = 0.0
      wrk1(:,:,:)     = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = transport_convert*uhrho_et_submeso(i,j,k)*Grd%dyte(i,j)
               wrk1_v(i,j,k,2) = transport_convert*vhrho_nt_submeso(i,j,k)*Grd%dxtn(i,j)
               wrk1(i,j,k)     = transport_convert*wrho_bt_submeso(i,j,k)*Grd%dat(i,j)
            enddo
         enddo
      enddo
      if (id_tx_trans_submeso_adv > 0) then 
          used = send_data (id_tx_trans_submeso_adv, wrk1_v(:,:,:,1), &
               Time%model_time, rmask=Grd%tmask(:,:,:),               &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_ty_trans_submeso_adv > 0) then 
          used = send_data (id_ty_trans_submeso_adv, wrk1_v(:,:,:,2), &
               Time%model_time, rmask=Grd%tmask(:,:,:),               &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_tz_trans_submeso_adv > 0) then 
          used = send_data (id_tz_trans_submeso_adv, wrk1(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:),           &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

  endif


  ! mass transports (skew approach diagnoses transports in compute_submeso_skewsion).
  ! to compute overturning streamfunction in Ferret, do so just like the Eulerian overturning. 
  if(submeso_advect_flux) then 
      if (id_tx_trans_submeso > 0 .or. id_ty_trans_submeso > 0 .or. id_tz_trans_submeso > 0) then 

          wrk1_v = 0.0
          wrk1   = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1_v(i,j,k,1) = transport_convert*uhrho_et_submeso(i,j,k)*Grd%dyte(i,j)
                   wrk1_v(i,j,k,2) = transport_convert*vhrho_nt_submeso(i,j,k)*Grd%dxtn(i,j)
                   wrk1(i,j,k)     = transport_convert*wrho_bt_submeso(i,j,k)*Grd%dat(i,j)
                enddo
             enddo
          enddo
          call diagnose_3d(Time, Grd, id_tx_trans_submeso, wrk1_v(:,:,:,1))
          call diagnose_3d(Time, Grd, id_ty_trans_submeso, wrk1_v(:,:,:,2))
          call diagnose_3d(Time, Grd, id_tz_trans_submeso, wrk1(:,:,:))

      endif
      call transport_on_nrho_submeso_adv(Time, Dens, wrk1_v(:,:,:,1), wrk1_v(:,:,:,2))  
  endif

  if (diag_step > 0) then
      if (mod(Time%itt,diag_step) == 0) then
          call maximum_bottom_w_submeso(wrho_bt_submeso(:,:,1:nk),Grd%xt,Grd%yt,Thickness%depth_zwt,rho0r,Grd%kmt)
      endif
  endif

  if(debug_this_module) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_submeso_mod: chksums for transport'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('uhrho_et_submeso', uhrho_et_submeso(COMP,:)*Grd%tmask(COMP,:))
      call write_chksum_3d('vhrho_nt_submeso', vhrho_nt_submeso(COMP,:)*Grd%tmask(COMP,:))
      call write_chksum_3d('wrho_bt_submeso', wrho_bt_submeso(COMP,1:nk)*Grd%tmask(COMP,:))
  endif

end subroutine compute_advect_transport
! </SUBROUTINE> NAME="compute_advect_transport"


!#######################################################################
! <SUBROUTINE NAME="compute_submeso_skewsion">
!
! <DESCRIPTION>
! Compute tendency from submeso skewsion. 
! </DESCRIPTION>
!
subroutine compute_submeso_skewsion(Thickness, Dens, Time, T_prog)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real, dimension(isd:ied,jsd:jed) :: eta_tend
  real    :: eta_tend_glob
  integer :: i,j,k,n,tau

  tau = Time%tau 
  
  do n=1,num_prog_tracers

     call compute_flux_x(Time,n,T_prog(n))
     call compute_flux_y(Time,n,T_prog(n))
     call compute_flux_z(Time,n,T_prog(n))

     if (Grd%tripolar) then 
         call mpp_update_domains(flux_x(:,:,:), flux_y(:,:,:), Dom_flux_sub%domain2d, &
              gridtype=CGRID_NE) 
     endif

     ! tracer tendency (units rho*dzt * tracer concentration/sec)
     T_prog(n)%wrk1 = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k)*(flux_z(i,j,k-1)-flux_z(i,j,k)           &
                   +(flux_x(i,j,k)-flux_x(i-1,j,k)+flux_y(i,j,k)-flux_y(i,j-1,k))*Grd%datr(i,j) &
                   )
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
           enddo
        enddo
     enddo

     if(id_submeso(n) > 0) then
        call diagnose_3d(Time, Grd, id_submeso(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion)
     endif
     if(id_submeso_on_nrho(n) > 0) then
        call diagnose_3d_rho(Time, Dens, id_submeso_on_nrho(n), T_prog(n)%wrk1*T_prog(n)%conversion)
     endif

  enddo ! enddo for n=1,num_prog_tracers


  ! diagnose mass transports for sending to diagnostics.
  wrk1_v = 0.0
  if(vert_coordinate_class==DEPTH_BASED) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = &
                    -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*rho0
               wrk1_v(i,j,k,2) = &
                    transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*rho0
            enddo
         enddo
      enddo
  else
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = &
                    -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*Dens%rho(i,j,k,tau)
               wrk1_v(i,j,k,2) = &
                    transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*Dens%rho(i,j,k,tau)
            enddo
         enddo
      enddo
  endif

  ! change signs to agree with convention used for ty_trans_gm 
  if (id_tx_trans_submeso > 0) then 
     call diagnose_3d(Time, Grd, id_tx_trans_submeso, -1.0*wrk1_v(:,:,:,1))
  endif
  if (id_ty_trans_submeso > 0) then 
     call diagnose_3d(Time, Grd, id_ty_trans_submeso, -1.0*wrk1_v(:,:,:,2))
  endif

  call transport_on_nrho_submeso(Time, Dens,-1.0*wrk1_v(:,:,:,1),-1.0*wrk1_v(:,:,:,2))  

  if(diagnose_eta_tend_submeso_flx) then 
      call diagnose_eta_tend_3dflux (Time, Thickness, Dens,&
           flux_x_temp(:,:,:),                             &
           flux_y_temp(:,:,:),                             &
           flux_z_temp(:,:,:),                             &
           flux_x_salt(:,:,:),                             &
           flux_y_salt(:,:,:),                             &
           flux_z_salt(:,:,:),                             &    
           eta_tend, eta_tend_glob)   

      call diagnose_2d(Time, Grd, id_eta_tend_submeso_flx, eta_tend(:,:))
      if(id_eta_tend_submeso_flx_glob > 0) then  
          used  = send_data (id_eta_tend_submeso_flx_glob, eta_tend_glob, Time%model_time)
      endif
  endif


end subroutine compute_submeso_skewsion
! </SUBROUTINE> NAME="compute_submeso_skewsion"



!#######################################################################
! <SUBROUTINE NAME="compute_flux_x">
!
! <DESCRIPTION>
!
! Subroutine computes the zonal submesoscale tracer skew flux component.
!
! fx has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine compute_flux_x(Time,n,Tracer)

  type(ocean_time_type),        intent(in) :: Time
  integer,                      intent(in) :: n
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  integer :: i, j, k
  integer :: ip, kr
  real    :: tensor_13(isd:ied,jsd:jed,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)

  flux_x = 0.0

  do k=1,nk

     ! tracer-independent part of the calculation 
     tensor_13(:,:,:) = 0.0
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec
              tensor_13(i,j,ip) = -psiy(i,j,k,ip)
           enddo
        enddo
     enddo

     ! tracer-dependent part of the calculation
     sumz(:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                      *tensor_13(i,j,ip)*dTdz(n)%field(i+ip,j,k-1+kr)                              &
                      *min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc,jec
        do i=isc-1,iec
           flux_x(i,j,k) = Grd%dxter(i,j)*(sumz(i,j,0)+sumz(i,j,1))
        enddo
     enddo

     if(submeso_limit_flux) then 
         do j=jsc,jec
            do i=isc-1,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then 
                   flux_x(i,j,k) = 0.0
               endif
            enddo
         enddo
     endif

  enddo   ! k-loop


  ! send fluxes to diag_manager 
  if(id_xflux_submeso(n) > 0) then 
     call diagnose_3d(Time, Grd, id_xflux_submeso(n), flux_sign*Tracer%conversion*flux_x(:,:,:))
  endif

  if(id_xflux_submeso_int_z(n) > 0) then 
      wrk1_2d = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + flux_x(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_xflux_submeso_int_z(n), flux_sign*Tracer%conversion*wrk1_2d(:,:))
  endif

  ! save for eta_tend diagnostics 
  if(diagnose_eta_tend_submeso_flx) then 
      if(n==index_temp) then
          flux_x_temp(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   flux_x_temp(i,j,k) = flux_x(i,j,k)
                enddo
             enddo
          enddo
      elseif(n==index_salt) then
          flux_x_salt(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   flux_x_salt(i,j,k) = flux_x(i,j,k)
                enddo
             enddo
          enddo
      endif
  endif


end subroutine compute_flux_x
! </SUBROUTINE> NAME="compute_flux_x"



!#######################################################################
! <SUBROUTINE NAME="compute_flux_y">
!
! <DESCRIPTION>
!
! Subroutine computes the meridional submesoscale tracer skew flux component.
!
! fy has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine compute_flux_y(Time,n,Tracer)

  type(ocean_time_type),        intent(in) :: Time
  integer,                      intent(in) :: n
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  integer :: i, j, k
  integer :: jq, kr
  real    :: tensor_23(isd:ied,jsd:jed,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)

  flux_y = 0.0

  do k=1,nk

     ! tracer-independent part of the calculation 
     tensor_23(:,:,:) = 0.0
     do jq=0,1  
        do j=jsc-1,jec
           do i=isc,iec
              tensor_23(i,j,jq) = psix(i,j,k,jq)
           enddo
        enddo
     enddo

     ! tracer-dependent part of the calculation
     sumz(:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                      *tensor_23(i,j,jq)*dTdz(n)%field(i,j+jq,k-1+kr)                              &
                      *min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc-1,jec
        do i=isc,iec
           flux_y(i,j,k) = Grd%dytnr(i,j)*(sumz(i,j,0)+sumz(i,j,1))
        enddo
     enddo

     if(submeso_limit_flux) then 
         do j=jsc-1,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then 
                   flux_y(i,j,k) = 0.0
               endif
            enddo
         enddo
     endif

  enddo   ! k-loop


  ! diagnostics 
  if(id_yflux_submeso(n) > 0) then 
     call diagnose_3d(Time, Grd, id_yflux_submeso(n), flux_sign*Tracer%conversion*flux_y(:,:,:))
  endif

  if(id_yflux_submeso_int_z(n) > 0) then 
      wrk1_2d = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + flux_y(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_yflux_submeso_int_z(n), flux_sign*Tracer%conversion*wrk1_2d(:,:))
  endif

  ! save for eta_tend diagnostics 
  if(diagnose_eta_tend_submeso_flx) then 
      if(n==index_temp) then
          flux_y_temp(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   flux_y_temp(i,j,k) = flux_y(i,j,k)
                enddo
             enddo
          enddo
      elseif(n==index_salt) then
          flux_y_salt(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   flux_y_salt(i,j,k) = flux_y(i,j,k)
                enddo
             enddo
          enddo
      endif
  endif


end subroutine compute_flux_y
! </SUBROUTINE> NAME="compute_flux_y"


!#######################################################################
! <SUBROUTINE NAME="compute_flux_z">
!
! <DESCRIPTION>
!
! Subroutine computes the vertical submeso tracer skew flux component.
!
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
!
! fz has physical dimensions (density*diffusivity*tracer gradient)
!
! </DESCRIPTION>
!
subroutine compute_flux_z(Time,n,Tracer)

  type(ocean_time_type),        intent(in) :: Time
  integer,                      intent(in) :: n
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  integer :: i, j, k, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1
  real :: temparray31(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32(isc:iec,jsc:jec,0:1,0:1)
  real :: tensor_31(isd:ied,jsd:jed,0:1)
  real :: tensor_32(isd:ied,jsd:jed,0:1)

  flux_z      = 0.0
  tensor_31   = 0.0
  tensor_32   = 0.0
  temparray31 = 0.0
  temparray32 = 0.0

  do k=1,nk-1

     do ip=0,1
        jq=ip  
        do j=jsc-1,jec
           do i=isc-1,iec
              tensor_31(i,j,ip) =  psiy(i,j,k,ip)
              tensor_32(i,j,jq) = -psix(i,j,k,jq)
           enddo
        enddo
     enddo

     do kr=0,1

        do ip=0,1
           do j=jsc,jec
              do i=isc,iec
                 temparray31(i,j,ip,kr) = tensor_31(i,j,ip)*dtew(i,j,ip) &
                      *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
              enddo
           enddo
        enddo
        do jq=0,1
           do j=jsc,jec
              do i=isc,iec
                 temparray32(i,j,jq,kr) = tensor_32(i,j,jq)*dtns(i,j,jq) &
                      *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
              enddo
           enddo
        enddo

     enddo

     do j=jsc,jec
        do i=isc,iec
           sumx_0 =  temparray31(i,j,0,0)*dTdx(n)%field(i-1,j,k) &
                  +  temparray31(i,j,0,1)*dTdx(n)%field(i-1,j,k+1)
           sumx_1 =  temparray31(i,j,1,0)*dTdx(n)%field(i,j,k)   &
                  +  temparray31(i,j,1,1)*dTdx(n)%field(i,j,k+1)
           sumy_0 =  temparray32(i,j,0,0)*dTdy(n)%field(i,j-1,k) &
                  +  temparray32(i,j,0,1)*dTdy(n)%field(i,j-1,k+1)
           sumy_1 =  temparray32(i,j,1,0)*dTdy(n)%field(i,j,k)   &
                  +  temparray32(i,j,1,1)*dTdy(n)%field(i,j,k+1)

           flux_z(i,j,k) = Grd%tmask(i,j,k+1)                                      &
                *( Grd%dxtr(i,j)*(sumx_0+sumx_1) + Grd%dytr(i,j)*(sumy_0+sumy_1) ) &
                *dzwtr(i,j,k)

        enddo
     enddo

     if(submeso_limit_flux) then 
         do j=jsc,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then 
                   flux_z(i,j,k) = 0.0
               endif
            enddo
         enddo
     endif


  enddo   ! end of k=1,nk-1 loop


  ! diagnostics 
  if(id_zflux_submeso(n) > 0) then 
      wrk2 = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk2(i,j,k) = Grd%dat(i,j)*flux_z(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_zflux_submeso(n), flux_sign*Tracer%conversion*wrk2(:,:,1:nk))
  endif

  ! save for eta_tend diagnostics 
  if(diagnose_eta_tend_submeso_flx) then 
      if(n==index_temp) then
          flux_z_temp(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   flux_z_temp(i,j,k) = flux_z(i,j,k)
                enddo
             enddo
          enddo
      elseif(n==index_salt) then
          flux_z_salt(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   flux_z_salt(i,j,k) = flux_z(i,j,k)
                enddo
             enddo
          enddo
      endif
  endif

end subroutine compute_flux_z
! </SUBROUTINE> NAME="compute_flux_z"


!#######################################################################
! <SUBROUTINE NAME="compute_submeso_upwind">
!
! <DESCRIPTION>
! First order upwind to compute the tendency from submeso advection.
!
! Although this method adds diffusion, some of the mixing 
! is physically relevant. Absent this mixing, the submesoscale 
! parameterization is incomplete. The submesoscale parameterization is, 
! afterall, active only in the mixed layer, where there is lot of 
! physical mixing.  
!
! Use of first order upwind ensures that the tendency computed
! from submesoscale parameterization will not, in principle, 
! introduce extrema. However, there remain some issues with large
! tendencies appearing near boundaries that may compromise this 
! monotonicity property.   
!
! Apply masks so that there is no flux leaving cell next to bottom.  
!
! </DESCRIPTION>
!
subroutine compute_submeso_upwind(Time, Dens, T_prog)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real,dimension(isc:iec,jsc:jec) :: ft1
  real,dimension(isc:iec,jsc:jec) :: ft2
  real,dimension(isc:iec,jsc:jec) :: fe
  real,dimension(isc:iec,jsc:jec) :: fn

  integer  :: i, j, k, n, kp1, kp2
  integer  :: tau, taum1
  real     :: velocity, upos, uneg, wpos, wneg

  tau     = Time%tau
  taum1   = Time%taum1

  ! do loop over all tracers 
  do n=1,num_prog_tracers

     ! initialize some arrays 
     advect_tendency(:,:,:) = 0.0
     flux_x(:,:,:)          = 0.0
     flux_y(:,:,:)          = 0.0
     flux_z(:,:,:)          = 0.0

     ! vertical flux divergence 
     ft1 = 0.0
     do k=1,nk
        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        do j=jsc,jec
           do i=isc,iec
              velocity = 0.5*wrho_bt_submeso(i,j,k)
              wpos     = velocity + abs(velocity) 
              wneg     = velocity - abs(velocity) 
              ft2(i,j) = (wneg*T_prog(n)%field(i,j,k,taum1) + wpos*T_prog(n)%field(i,j,kp1,taum1)) &
                         *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)*Grd%tmask(i,j,kp2)
              advect_tendency(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
              ft1(i,j) = ft2(i,j)
              flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
           enddo
        enddo
     enddo

     ! horizontal fluxes
     do k=1,nk 
        kp1 = min(k+1,nk)
 
        ! i-flux
        do j=jsc,jec
           do i=isc-1,iec
              velocity = 0.5*uhrho_et_submeso(i,j,k)
              upos     = velocity + abs(velocity)
              uneg     = velocity - abs(velocity)
              fe(i,j)  = Grd%dyte(i,j)*(upos*T_prog(n)%field(i,j,k,taum1) + uneg*T_prog(n)%field(i+1,j,k,taum1)) &
                         *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)*Grd%tmask(i,j,kp1)*Grd%tmask(i+1,j,kp1)
              flux_x(i,j,k) = fe(i,j)
           enddo
        enddo

        ! j-flux
        do j=jsc-1,jec
           do i=isc,iec
              velocity = 0.5*vhrho_nt_submeso(i,j,k)
              upos     = velocity + abs(velocity)
              uneg     = velocity - abs(velocity)
              fn(i,j)  = Grd%dxtn(i,j)*(upos*T_prog(n)%field(i,j,k,taum1) + uneg*T_prog(n)%field(i,j+1,k,taum1)) &
                         *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)*Grd%tmask(i,j,kp1)*Grd%tmask(i,j+1,kp1)
              flux_y(i,j,k) = fn(i,j)
           enddo
        enddo

     enddo  ! k-loop completed

     ! needed for tripolar 
     call mpp_update_domains(flux_x(:,:,:), flux_y(:,:,:), Dom_flux_sub%domain2d, gridtype=CGRID_NE) 

 
     ! advect_tendency is a divergence; minus sign produces a convergence needed for tendency 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              advect_tendency(i,j,k) = advect_tendency(i,j,k) + &
              Grd%tmask(i,j,k)*(flux_x(i,j,k)-flux_x(i-1,j,k)+flux_y(i,j,k)-flux_y(i,j-1,k))*Grd%datr(i,j)
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) - advect_tendency(i,j,k)
           enddo
        enddo
     enddo


     ! send fluxes to diag_manager 
     if(id_xflux_submeso(n) > 0) then 
        call diagnose_3d(Time, Grd, id_xflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_x(:,:,:))
     endif
     if(id_xflux_submeso_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_x(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_xflux_submeso_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:))
     endif

     if(id_yflux_submeso(n) > 0) then 
        call diagnose_3d(Time, Grd, id_yflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_y(:,:,:))
     endif
     if(id_yflux_submeso_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_y(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_yflux_submeso_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:))
     endif

     ! sans the (1:nk) argument in flux_z, the diagnostics output spuriously has flux_z(k=1)=0.0.
     if(id_zflux_submeso(n) > 0) then 
        call diagnose_3d(Time, Grd, id_zflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_z(:,:,1:nk))
     endif
     
     ! minus sign produces a convergence rather than divergence, 
     ! which is consistent with upwind computed in tracer_advect module.       
     if(id_submeso(n) > 0) then 
        call diagnose_3d(Time, Grd, id_submeso(n), -advect_tendency(:,:,:)*T_prog(n)%conversion)
     endif
     if(id_submeso_on_nrho(n) > 0) then
        call diagnose_3d_rho(Time, Dens, id_submeso_on_nrho(n), -advect_tendency*T_prog(n)%conversion)
     endif


  enddo ! end of tracer n-loop
  
end subroutine compute_submeso_upwind
! </SUBROUTINE> NAME="compute_submeso_upwind"



!#######################################################################
! <SUBROUTINE NAME="compute_submeso_sweby">
!
! <DESCRIPTION>
! Sweby scheme to compute the tendency from submeso advection.
! Algorithm taken after advect_tracer_sweby_all in the module
! ocean_tracers/ocean_tracer_advect.F90.
!
! Jan 2012: Stephen.Griffies
! This scheme has known bugs; it is not meant for general use. 
!
! </DESCRIPTION>
!
subroutine compute_submeso_sweby(Thickness, Time, Dens, T_prog)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real,dimension(isc:iec,jsc:jec) :: ftp
  real,dimension(isc:iec,jsc:jec) :: fbt
  real,dimension(isc:iec,jsc:jec) :: wkm1

  integer  :: i, j, k, n
  integer  :: kp1, kp2, km1
  integer  :: tau, taum1
  real     :: Rjm, Rj, Rjp, cfl, massflux
  real     :: d0, d1, thetaP, psiP 
  real     :: thetaM, psiM

  tau     = Time%tau
  taum1   = Time%taum1
  flux_x  = 0.0
  flux_y  = 0.0
  flux_z  = 0.0
  wrk1    = 0.0

  ! calculate flux at bottom face of the T-cells
  do n=1,num_prog_tracers

     ftp  = 0.0
     fbt  = 0.0
     wkm1 = 0.0

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(n)%wrk1(i,j,k) = 0.0
           enddo
        enddo
     enddo

     do k=1,nk

        do j=jsc,jec
           do i=isc,iec
              tracer_mdfl_all(n)%field(i,j,k) = T_prog(n)%field(i,j,k,taum1)
           enddo
        enddo

        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        km1 = max(k-1,1)

        do j=jsc,jec
           do i=isc,iec

              Rjp = (T_prog(n)%field(i,j,km1,taum1) - T_prog(n)%field(i,j,k,taum1))    &
                   *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
              Rj  = (T_prog(n)%field(i,j,k,taum1) - T_prog(n)%field(i,j,kp1,taum1))    &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
              Rjm = (T_prog(n)%field(i,j,kp1,taum1) - T_prog(n)%field(i,j,kp2,taum1))  &
                   *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

              massflux = Grd%dat(i,j) * wrho_bt_submeso(i,j,k)
              cfl = abs(wrho_bt_submeso(i,j,k) * dtime  / Thickness%rho_dzt(i,j,k,tau))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              fbt(i,j) =  0.5 * ( ( massflux + abs(massflux) )         &
                   * ( T_prog(n)%field(i,j,kp1,taum1) + psiP * Rj )    &
                   + ( massflux - abs(massflux) )                      &
                   * ( T_prog(n)%field(i,j, k ,taum1) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

              wrk1(i,j,k) = Grd%datr(i,j)*( fbt(i,j) - ftp(i,j))     &
                   + T_prog(n)%field(i,j,k,taum1)*(wkm1(i,j) - wrho_bt_submeso(i,j,k)) 

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)  &
                   + wrk1(i,j,k) * dtime/Thickness%rho_dzt(i,j,k,tau) 

              flux_z(i,j,k) = fbt(i,j)
              ftp(i,j)      = fbt(i,j)

           enddo
        enddo

        ! update vertical velocity for next k-level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = wrho_bt_submeso(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                              flags=XUPDATE, complete=T_prog(n)%complete)

  enddo ! end of n-loop for tracers  


  ! calculate flux at the eastern face of the T-cells
  do n=1,num_prog_tracers

     do k=1,nk

        do j=jsc,jec
           do i=isc-1,iec

              Rjp = (tracer_mdfl_all(n)%field(i+2,j,k) - tracer_mdfl_all(n)%field(i+1,j,k))    &
                   *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
              Rj  = (tracer_mdfl_all(n)%field(i+1,j,k) - tracer_mdfl_all(n)%field( i ,j,k) )   &
                   *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k) - tracer_mdfl_all(n)%field(i-1,j,k))      &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

              massflux = Grd%dyte(i,j) * uhrho_et_submeso(i,j,k)
              cfl = abs(uhrho_et_submeso(i,j,k) * dtime * 2.0    &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau)) * Grd%dxte(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_x(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )      &
                   * ( tracer_mdfl_all(n)%field( i ,j,k) + psiP * Rj )   &
                   + ( massflux - abs(massflux) )                        &
                   * ( tracer_mdfl_all(n)%field(i+1,j,k) - psiM * Rj ) ) &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
           enddo
        enddo

        ! update the tracer
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = tmask_mdfl(i,j,k) * Grd%datr(i,j)       &
                   * (flux_x(i-1,j,k) - flux_x(i,j,k)               &
                   + T_prog(n)%field(i,j,k,taum1)* (                &
                     Grd%dyte(i,j)   * uhrho_et_submeso( i ,j,k)    &
                   - Grd%dyte(i-1,j) * uhrho_et_submeso(i-1,j,k) ) )

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k) &
                   + wrk1(i,j,k)*dtime/Thickness%rho_dzt(i,j,k,tau)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                              flags=YUPDATE, complete=T_prog(n)%complete)

  enddo ! end of n-loop for tracers  


  ! calculate flux at the northern face of the T-cells 
  do n=1,num_prog_tracers 

     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = 0.0
        enddo
     enddo

     do k=1,nk

        do j=jsc-1,jec
           do i=isc,iec

              Rjp = (tracer_mdfl_all(n)%field(i,j+2,k) - tracer_mdfl_all(n)%field(i,j+1,k))   &
                   *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
              Rj  = (tracer_mdfl_all(n)%field(i,j+1,k) - tracer_mdfl_all(n)%field(i,j,k))     &
                   *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k)   - tracer_mdfl_all(n)%field(i,j-1,k))   &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

              massflux = Grd%dxtn(i,j) * vhrho_nt_submeso(i,j,k)
              cfl = abs(vhrho_nt_submeso(i,j,k) * dtime * 2.0              &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau)) * Grd%dytn(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_y(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )       &
                   * ( tracer_mdfl_all(n)%field(i,j,k) + psiP * Rj )      &
                   + ( massflux - abs(massflux) )                         &
                   * ( tracer_mdfl_all(n)%field(i,j+1,k) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

           enddo
        enddo

        ! calculate the overall tendency (advection convergence) and update tracer 
        do j=jsc,jec
           do i=isc,iec

              wrk1(i,j,k) = tmask_mdfl(i,j,k)*Grd%datr(i,j)*(flux_y(i,j-1,k)-flux_y(i,j,k))  &
                            + T_prog(n)%field(i,j,k,taum1) *                                 &
                   (wrho_bt_submeso(i,j,k) - wkm1(i,j)                                       &
                   + Grd%datr(i,j)*                                                          &
                   ( Grd%dyte(i-1,j) * uhrho_et_submeso(i-1,j,k)                             &
                   - Grd%dyte( i ,j) * uhrho_et_submeso( i ,j,k)))

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k) &
                                              + wrk1(i,j,k)*dtime/Thickness%rho_dzt(i,j,k,tau)

              T_prog(n)%wrk1(i,j,k) = &
                Thickness%rho_dzt(i,j,k,tau)*(tracer_mdfl_all(n)%field(i,j,k)-T_prog(n)%field(i,j,k,taum1))/dtime  &
                   *tmask_mdfl(i,j,k) 

           enddo
        enddo

        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
          enddo
        enddo

        ! update vertical velocity for next level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = wrho_bt_submeso(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop



     ! send fluxes to diag_manager 
     if(id_xflux_submeso(n) > 0) then 
        call diagnose_3d(Time, Grd, id_xflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_x(:,:,:))
     endif
     if(id_xflux_submeso_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_x(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_xflux_submeso_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:))
     endif

     if(id_yflux_submeso(n) > 0) then 
        call diagnose_3d(Time, Grd, id_yflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_y(:,:,:))
     endif
     if(id_yflux_submeso_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_y(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_yflux_submeso_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:))
     endif

     if(id_zflux_submeso(n) > 0) then 
         wrk2 = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk2(i,j,k) = Grd%dat(i,j)*flux_z(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_3d(Time, Grd, id_zflux_submeso(n), flux_sign*T_prog(n)%conversion*wrk2(:,:,1:nk))
     endif
     
     if(id_submeso(n) > 0) then
        call diagnose_3d(Time, Grd, id_submeso(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion)
     endif

     if(id_submeso_on_nrho(n) > 0) then
        call diagnose_3d_rho(Time, Dens, id_submeso_on_nrho(n), T_prog(n)%wrk1*T_prog(n)%conversion)
     endif


  enddo ! end of n-loop
  
end subroutine compute_submeso_sweby
! </SUBROUTINE> NAME="compute_submeso_sweby"


!#######################################################################
! <SUBROUTINE NAME="compute_submeso_diffusion">
!
! <DESCRIPTION>
! Compute tendency from submeso horizontal diffusion.
! </DESCRIPTION>
!
subroutine compute_submeso_diffusion(Thickness, Dens, Time, T_prog)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real, dimension(isd:ied,jsd:jed)    :: eta_tend
  real, dimension(isd:ied,jsd:jed)    :: fe
  real, dimension(isd:ied,jsd:jed)    :: fn
  real, dimension(isd:ied,jsd:jed,nk) :: del2_tracer 
  real    :: eta_tend_glob
  integer :: i,j,k,n,tau

  tau = Time%tau 

  ! compute effective horizontal diffusivity (m^2/sec)
  wrk1(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = onefourth*Grd%tmask(i,j,k) &
                         *(psix(i,j,k,0)**2 + psix(i,j,k,1)**2 + psiy(i,j,k,0)**2 + psiy(i,j,k,1)**2) 
           wrk1(i,j,k) = submeso_diffusion_scale*(wrk1(i,j,k)**0.5)
        enddo
     enddo
  enddo

  ! biharmonic mixing coefficient (m^4/sec) = laplacian(m^2/sec) * horz_area(m^2) 
  if(submeso_diffusion_biharmonic) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = wrk1(i,j,k)*grid_length(i,j)**2
            enddo
         enddo
      enddo
  endif


  ! compute tendency from horizontal mixing in submeso boundary layer 
  do n=1,num_prog_tracers


     if(submeso_diffusion_biharmonic) then 

         ! del2_tracer = laplacian operator with unit diffusivity.
         ! del2_tracer has units (tracer concentration / m^2) 
         del2_tracer(:,:,:) = 0.0 
         do k=1,nk
            fe(:,:) = dTdx(n)%field(:,:,k)*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
            fn(:,:) = dTdy(n)%field(:,:,k)*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
            del2_tracer(:,:,k) = &
            Grd%tmask(:,:,k)*(BDX_ET(fe(:,:)) + BDY_NT(fn(:,:)))/(epsln+Thickness%rho_dzt(:,:,k,tau))
         enddo
         call mpp_update_domains (del2_tracer, Dom%domain2d)

         ! compute biharmonic tracer flux components. 
         ! fluxes have dimensions (area*diffusivity*density*tracer gradient)
         flux_x = 0.0
         flux_y = 0.0
         do k=1,nk
            flux_x(:,:,k) = &
              -wrk1(:,:,k)*Grd%dyte(:,:)*FDX_T(del2_tracer(:,:,k))*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
            flux_y(:,:,k) = &
              -wrk1(:,:,k)*Grd%dxtn(:,:)*FDY_T(del2_tracer(:,:,k))*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         enddo

     else  

         ! laplacian tracer flux components. 
         ! fluxes have dimensions (area*diffusivity*density*tracer gradient)
         flux_x = 0.0
         flux_y = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  flux_x(i,j,k) = wrk1(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*Grd%dyte(i,j)*dTdx(n)%field(i,j,k)
                  flux_y(i,j,k) = wrk1(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*Grd%dxtn(i,j)*dTdy(n)%field(i,j,k)
               enddo
            enddo
         enddo

     endif

     call mpp_update_domains(flux_x(:,:,:), flux_y(:,:,:), Dom_flux_sub%domain2d, gridtype=CGRID_NE) 


     ! tracer tendency (units rho*dzt * tracer concentration/sec)
     T_prog(n)%wrk1(:,:,:) = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k) &
               *(flux_x(i,j,k)-flux_x(i-1,j,k)+flux_y(i,j,k)-flux_y(i,j-1,k))*Grd%datr(i,j)
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
           enddo
        enddo
     enddo


     ! diagnostics

     if(id_subdiff(n) > 0) then 
        call diagnose_3d(Time, Grd, id_subdiff(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion)
     endif
     if(id_xflux_subdiff(n) > 0) then 
        call diagnose_3d(Time, Grd, id_xflux_subdiff(n), flux_sign*T_prog(n)%conversion*flux_x(:,:,:))
     endif

     if(id_xflux_subdiff_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_x(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_xflux_subdiff_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:))
     endif

     if(id_yflux_subdiff(n) > 0) then 
        call diagnose_3d(Time, Grd, id_yflux_subdiff(n), flux_sign*T_prog(n)%conversion*flux_y(:,:,:))
     endif
     if(id_yflux_subdiff_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_y(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_yflux_subdiff_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:))
     endif


     ! save fluxes for eta_tend diagnostics 
     if(diagnose_eta_tend_subdiff_flx) then 
         if(n==index_temp) then
             flux_x_temp(:,:,:) = 0.0
             flux_y_temp(:,:,:) = 0.0
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      flux_x_temp(i,j,k) = flux_x(i,j,k)
                      flux_y_temp(i,j,k) = flux_y(i,j,k)
                   enddo
                enddo
             enddo
         elseif(n==index_salt) then
             flux_x_salt(:,:,:) = 0.0
             flux_y_salt(:,:,:) = 0.0
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      flux_x_salt(i,j,k) = flux_x(i,j,k)
                      flux_y_salt(i,j,k) = flux_y(i,j,k)
                   enddo
                enddo
             enddo
         endif
     endif


  enddo ! enddo for n=1,num_prog_tracers


  ! save diffusivity 
  call diagnose_3d(Time, Grd, id_subdiff_diffusivity, wrk1(:,:,:))


  ! diagnose effects on sea level 
  if(diagnose_eta_tend_subdiff_flx) then 
      flux_z_temp(:,:,:) = 0.0
      flux_z_salt(:,:,:) = 0.0
      call diagnose_eta_tend_3dflux (Time, Thickness, Dens,&
           flux_x_temp(:,:,:),                             &
           flux_y_temp(:,:,:),                             &
           flux_z_temp(:,:,:),                             &
           flux_x_salt(:,:,:),                             &
           flux_y_salt(:,:,:),                             &
           flux_z_salt(:,:,:),                             &    
           eta_tend, eta_tend_glob)   

      call diagnose_2d(Time, Grd, id_eta_tend_subdiff_flx, eta_tend(:,:))
      if(id_eta_tend_subdiff_flx_glob > 0) then  
          used  = send_data (id_eta_tend_subdiff_flx_glob, eta_tend_glob, Time%model_time)
      endif
  endif


end subroutine compute_submeso_diffusion
! </SUBROUTINE> NAME="compute_submeso_diffusion"



!#######################################################################
! <SUBROUTINE NAME="maximum_bottom_w_general">
!
! <DESCRIPTION>
! Compute maximum vertical velocity from submeso.
! </DESCRIPTION>
!
subroutine maximum_bottom_w_submeso(wb,x_array,y_array,depth_array,scale,km_array)

  real,    dimension(isd:,jsd:,:) , intent(in) :: wb
  real,    dimension(isd:,jsd:)   , intent(in) :: x_array
  real,    dimension(isd:,jsd:)   , intent(in) :: y_array
  real,    dimension(isd:,jsd:,:) , intent(in) :: depth_array
  real,                             intent(in) :: scale
  integer, dimension(isd:,jsd:)   , intent(in) :: km_array

  real    :: wbot, wbot0, fudge
  integer :: i, j, k, iwbot, jwbot, kwbot

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_submesoscale_mod (maximum_bottom_w_submeso): module needs initialization ')
  endif 

  wbot=0.0; iwbot=isc; jwbot=jsc; kwbot=1
  do j=jsc,jec
    do i=isc,iec
      k = km_array(i,j)
      if (k /= 0 .and. (abs(wb(i,j,k)) > abs(wbot))) then
        wbot  = wb(i,j,k)
        iwbot = i
        jwbot = j
        kwbot = k
      endif
    enddo
  enddo
  wbot = scale*abs(wbot)


  write (stdoutunit,'(//60x,a/)') ' Summary of bottom vertical velocity from submeso:'
  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when tracer is independent of processor
  wbot  = wbot*fudge
  wbot0 = wbot
  wbot  = abs(wbot)
  call mpp_max(wbot)

  if (abs(wbot0) == wbot) then
    wbot = wbot0/fudge
    write (unit,9912) wbot, iwbot+Dom%ioff, jwbot+Dom%joff, kwbot, &
                      x_array(iwbot,jwbot), y_array(iwbot,jwbot), depth_array(iwbot,jwbot,kwbot)
  endif
  
9912  format(/' Maximum submeso bottom velocity (',es10.3,' m/s){error}  at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m)')

end subroutine maximum_bottom_w_submeso
! </SUBROUTINE> NAME="maximum_bottom_wrho_submeso"


!#######################################################################
! <SUBROUTINE NAME="transport_on_nrho_submeso">
!
! <DESCRIPTION>
! Classify horizontal submeso mass transport according to neutral
! density classes. 
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_rho_submeso, a remapping is 
! performed, rather than the binning done for transport_on_nrho_submeso_adv.
!
! This is the same algorithm as used for GM skew fluxes on rho surfaces. 
!
! Caveat: Since the submeso scheme operates only in the mixed layer,
! there are difficulties mapping this transport to neutral density 
! layers.  The user should be mindful of the problems with this 
! remapping.  An alternative that may be more suitable is to use 
! Ferret to remap the time mean submeso transport to the time mean
! neutral density surfaces.  There are missing correlations, but for
! many purposes, the Ferret remapping may be preferable.  
!
! Briefly, the Ferret command is the following:
! 
! let ty_trans_nrho_submeso_new = ZAXREPLACE(TY_TRANS_SUBMESO,NEUTRAL_RHO,TY_TRANS_NRHO)
! where TY_TRANS_SUBMESO is the level-space transport
!       NEUTRAL_RHO is the level space version of the neutral density 
!       TY_TRANS_NRHO is any density space field whose vertical coordinates 
!                     are accessed for the remapping. 
! </DESCRIPTION>
!
subroutine transport_on_nrho_submeso (Time, Dens, tx_trans_lev, ty_trans_lev)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_lev
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_lev
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_rho, neutralrho_nk
  real    :: work(isd:ied,jsd:jed,size(Dens%neutralrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error transport_on_nrho_submeso (transport_on_nrho_submeso): module needs initialization ')
  endif 

  next_time = increment_time(Time%model_time, int(dtime), 0)

  if (need_data(id_tx_trans_nrho_submeso,next_time) .or. need_data(id_ty_trans_nrho_submeso,next_time)) then

      neutralrho_nk = size(Dens%neutralrho_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with neutralrho_ref < neutralrho(k=1),   work=0
      ! for (i,j) points with neutralrho_ref > neutralrho(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_lev from k-levels to neutralrho_nk-levels
      do k_rho=1,neutralrho_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%neutralrho_ref(k_rho) >  Dens%neutralrho(i,j,k)  ) then
                      if( Dens%neutralrho_ref(k_rho) <= Dens%neutralrho(i,j,k+1)) then 
                          W1= Dens%neutralrho_ref(k_rho)- Dens%neutralrho(i,j,k)
                          W2= Dens%neutralrho(i,j,k+1)  - Dens%neutralrho_ref(k_rho)
                          work(i,j,k_rho,1) = (tx_trans_lev(i,j,k+1)*W1 +tx_trans_lev(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                          work(i,j,k_rho,2) = (ty_trans_lev(i,j,k+1)*W1 +ty_trans_lev(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_rho=1,neutralrho_nk
         do j=jsc,jec
            do i=isc,iec
               work(i,j,k_rho,1) = work(i,j,k_rho,1)*Grd%tmask(i,j,1)
               work(i,j,k_rho,2) = work(i,j,k_rho,2)*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_nrho_submeso > 0) then 
          used = send_data (id_tx_trans_nrho_submeso, work(:,:,:,1), Time%model_time, &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif
      if (id_ty_trans_nrho_submeso > 0) then 
          used = send_data (id_ty_trans_nrho_submeso, work(:,:,:,2), Time%model_time, &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif

  endif

end subroutine transport_on_nrho_submeso
! </SUBROUTINE> NAME="transport_on_nrho_submeso"



!#######################################################################
! <SUBROUTINE NAME="transport_on_nrho_submeso_adv">
!
! <DESCRIPTION>
! Classify horizontal transport according to neutral density classes. 
!
! Based on transport_on_nrho in ocean_diag/ocean_adv_diag.F90.
!
! </DESCRIPTION>
!
subroutine transport_on_nrho_submeso_adv (Time, Dens, utrans, vtrans)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: utrans
  real, dimension(isd:,jsd:,:), intent(in) :: vtrans
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_rho, neutralrho_nk
  real    :: work1(isc:iec,jsc:jec,size(Dens%neutralrho_ref))
  real    :: work2(isc:iec,jsc:jec,size(Dens%neutralrho_ref))
  real    :: tmp(2,isc:iec,jsc:jec)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error transport_on_nrho_submeso_adv: module needs initialization ')
  endif 

  neutralrho_nk = size(Dens%neutralrho_ref(:))
  next_time = increment_time(Time%model_time, int(dtime), 0)

  if (need_data(id_tx_trans_nrho_submeso,next_time) .or. need_data(id_ty_trans_nrho_submeso,next_time)) then

      work1(:,:,:) = 0.0
      work2(:,:,:) = 0.0

      do k_rho=1,neutralrho_nk

         tmp(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if (k_rho == 1) then
                     if(Dens%neutralrho(i,j,k) < Dens%neutralrho_bounds(k_rho)) then 
                         tmp(1,i,j) = tmp(1,i,j) + utrans(i,j,k)
                         tmp(2,i,j) = tmp(2,i,j) + vtrans(i,j,k)
                     endif
                  elseif(k_rho < neutralrho_nk) then
                     if( (Dens%neutralrho_bounds(k_rho) <= Dens%neutralrho(i,j,k)) .and.  &
                         (Dens%neutralrho(i,j,k)        <  Dens%neutralrho_bounds(k_rho+1)) ) then 
                           tmp(1,i,j) = tmp(1,i,j) + utrans(i,j,k)
                           tmp(2,i,j) = tmp(2,i,j) + vtrans(i,j,k)
                     endif
                  else    ! if (k_rho == neutralrho_nk) then
                     if(Dens%neutralrho_bounds(k_rho) <= Dens%neutralrho(i,j,k)) then 
                         tmp(1,i,j) = tmp(1,i,j) + utrans(i,j,k)             
                         tmp(2,i,j) = tmp(2,i,j) + vtrans(i,j,k)
                     endif
                  endif
               enddo
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               work1(i,j,k_rho) = (tmp(1,i,j)+work1(i,j,k_rho))*Grd%dyte(i,j)*transport_convert*Grd%tmask(i,j,1)
               work2(i,j,k_rho) = (tmp(2,i,j)+work2(i,j,k_rho))*Grd%dxtn(i,j)*transport_convert*Grd%tmask(i,j,1)
            enddo
         enddo

      enddo

      if (id_tx_trans_nrho_submeso > 0) then 
          used = send_data (id_tx_trans_nrho_submeso, work1, Time%model_time, &
                            ks_in=1, ke_in=neutralrho_nk)
      endif
      if (id_ty_trans_nrho_submeso > 0) then 
          used = send_data (id_ty_trans_nrho_submeso, work2, Time%model_time, &
                            ks_in=1, ke_in=neutralrho_nk)
      endif

  endif

end subroutine transport_on_nrho_submeso_adv
! </SUBROUTINE> NAME="transport_on_nrho_submeso_adv"



!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag = .false. 

  
  ! for the submeso streamfunction parameterization 

  id_neut_rho_submeso = register_diag_field ('ocean_model', 'neut_rho_submeso',&
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'update of locally referenced potential density from submeso param',       &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_submeso > 0) compute_watermass_diag = .true.

  id_wdian_rho_submeso = register_diag_field ('ocean_model', 'wdian_rho_submeso',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'dianeutral mass transport due to submeso closure',                          &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_submeso > 0) compute_watermass_diag = .true.

  id_tform_rho_submeso = register_diag_field ('ocean_model', 'tform_rho_submeso',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'watermass transform due to submeso on levels (pre-layer binning)',          &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_submeso > 0) compute_watermass_diag = .true.

  id_neut_rho_submeso_on_nrho = register_diag_field ('ocean_model',                            &
    'neut_rho_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                    &
    'update of locally ref potential density from submeso as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_rho_submeso_on_nrho = register_diag_field ('ocean_model',                &
     'wdian_rho_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
     'dianeutral mass transport due to submeso as binned to neutral density layers',&
     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_rho_submeso_on_nrho = register_diag_field ('ocean_model',          &
     'tform_rho_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
     'watermass transform due to submeso as binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_eta_tend_submeso_tend= -1          
  id_eta_tend_submeso_tend= register_diag_field ('ocean_model','eta_tend_submeso_tend',&
       Grd%tracer_axes(1:2), Time%model_time,                                          &
       'non-Bouss steric sea level tendency from submesoscale tendency', 'm/s',        &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_submeso_tend > 0) compute_watermass_diag=.true.

  id_eta_tend_submeso_tend_glob= -1          
  id_eta_tend_submeso_tend_glob= register_diag_field ('ocean_model', 'eta_tend_submeso_tend_glob',&
       Time%model_time,                                                                           &
       'global mean non-bouss steric sea level tendency from submesoscale tendency',              &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_submeso_tend_glob > 0) compute_watermass_diag=.true.

  id_eta_tend_subdiff_tend= -1          
  id_eta_tend_subdiff_tend= register_diag_field ('ocean_model','eta_tend_subdiff_tend',&
       Grd%tracer_axes(1:2), Time%model_time,                                          &
       'non-Bouss steric sea level tendency from subdiff tendency', 'm/s',             &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_subdiff_tend > 0) compute_watermass_diag=.true.

  id_eta_tend_subdiff_tend_glob= -1          
  id_eta_tend_subdiff_tend_glob= register_diag_field ('ocean_model', 'eta_tend_subdiff_tend_glob',&
       Time%model_time,                                                                           &
       'global mean non-bouss steric sea level tendency from subdiff tendency',                   &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_subdiff_tend_glob > 0) compute_watermass_diag=.true.


  id_neut_temp_submeso = register_diag_field ('ocean_model', 'neut_temp_submeso',    &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'temp related update of locally referenced potential density from submeso param',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_submeso > 0) compute_watermass_diag = .true.

  id_wdian_temp_submeso = register_diag_field ('ocean_model', 'wdian_temp_submeso',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'temp related dianeutral mass transport due to submeso closure',               &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_submeso > 0) compute_watermass_diag = .true.

  id_tform_temp_submeso = register_diag_field ('ocean_model', 'tform_temp_submeso', &
    Grd%tracer_axes(1:3), Time%model_time,                                          &
    'temp related watermass transform due to submeso on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_submeso > 0) compute_watermass_diag = .true.

  id_neut_temp_submeso_on_nrho = register_diag_field ('ocean_model',                                        &
    'neut_temp_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                &
    'temp related update of locally ref potential density from submeso as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_temp_submeso_on_nrho = register_diag_field ('ocean_model',                            &
     'wdian_temp_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
     'temp related dianeutral mass transport due to submeso as binned to neutral density layers',&
     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_temp_submeso_on_nrho = register_diag_field ('ocean_model',                       &
     'tform_temp_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
     'temp related watermass transform due to submeso as binned to neutral density layers', &
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_submeso_on_nrho > 0) compute_watermass_diag = .true.



  id_neut_salt_submeso = register_diag_field ('ocean_model', 'neut_salt_submeso',    &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'salt related update of locally referenced potential density from submeso param',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_submeso > 0) compute_watermass_diag = .true.

  id_wdian_salt_submeso = register_diag_field ('ocean_model', 'wdian_salt_submeso',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'salt related dianeutral mass transport due to submeso closure',               &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_submeso > 0) compute_watermass_diag = .true.

  id_tform_salt_submeso = register_diag_field ('ocean_model', 'tform_salt_submeso', &
    Grd%tracer_axes(1:3), Time%model_time,                                          &
    'salt related watermass transform due to submeso on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_submeso > 0) compute_watermass_diag = .true.

  id_neut_salt_submeso_on_nrho = register_diag_field ('ocean_model',                                        &
    'neut_salt_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                &
    'salt related update of locally ref potential density from submeso as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_salt_submeso_on_nrho = register_diag_field ('ocean_model',                            &
     'wdian_salt_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
     'salt related dianeutral mass transport due to submeso as binned to neutral density layers',&
     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_submeso_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_salt_submeso_on_nrho = register_diag_field ('ocean_model',                      &
     'tform_salt_submeso_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,             &
     'salt related watermass transform due to submeso as binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_submeso_on_nrho > 0) compute_watermass_diag = .true.


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_submesoscale_mod w/ compute_watermass_diag=.true.'  
  endif 



  ! for the submeso horizontal diffusion parameterization 

  id_neut_rho_subdiff = register_diag_field ('ocean_model', 'neut_rho_subdiff',&
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'update of locally referenced potential density from subdiff param',       &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_subdiff > 0) compute_watermass_diff_diag = .true.

  id_wdian_rho_subdiff = register_diag_field ('ocean_model', 'wdian_rho_subdiff',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'dianeutral mass transport due to subdiff closure',                          &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_subdiff > 0) compute_watermass_diff_diag = .true.

  id_tform_rho_subdiff = register_diag_field ('ocean_model', 'tform_rho_subdiff',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'watermass transform due to subdiff on levels (pre-layer binning)',          &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_subdiff > 0) compute_watermass_diff_diag = .true.

  id_neut_rho_subdiff_on_nrho = register_diag_field ('ocean_model',                            &
    'neut_rho_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                    &
    'update of locally ref potential density from subdiff as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_wdian_rho_subdiff_on_nrho = register_diag_field ('ocean_model',                &
     'wdian_rho_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
     'dianeutral mass transport due to subdiff as binned to neutral density layers',&
     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_tform_rho_subdiff_on_nrho = register_diag_field ('ocean_model',          &
     'tform_rho_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
     'watermass transform due to subdiff as binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_eta_tend_subdiff_tend= -1          
  id_eta_tend_subdiff_tend= register_diag_field ('ocean_model','eta_tend_subdiff_tend',&
       Grd%tracer_axes(1:2), Time%model_time,                                          &
       'non-Bouss steric sea level tendency from subdiffscale tendency', 'm/s',        &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_subdiff_tend > 0) compute_watermass_diff_diag=.true.

  id_eta_tend_subdiff_tend_glob= -1          
  id_eta_tend_subdiff_tend_glob= register_diag_field ('ocean_model', 'eta_tend_subdiff_tend_glob',&
       Time%model_time,                                                                           &
       'global mean non-bouss steric sea level tendency from subdiffscale tendency',              &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_subdiff_tend_glob > 0) compute_watermass_diff_diag=.true.


  id_neut_temp_subdiff = register_diag_field ('ocean_model', 'neut_temp_subdiff',    &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'temp related update of locally referenced potential density from subdiff param',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_subdiff > 0) compute_watermass_diff_diag = .true.

  id_wdian_temp_subdiff = register_diag_field ('ocean_model', 'wdian_temp_subdiff',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'temp related dianeutral mass transport due to subdiff closure',               &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_subdiff > 0) compute_watermass_diff_diag = .true.

  id_tform_temp_subdiff = register_diag_field ('ocean_model', 'tform_temp_subdiff', &
    Grd%tracer_axes(1:3), Time%model_time,                                          &
    'temp related watermass transform due to subdiff on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_subdiff > 0) compute_watermass_diff_diag = .true.

  id_neut_temp_subdiff_on_nrho = register_diag_field ('ocean_model',                                        &
    'neut_temp_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                &
    'temp related update of locally ref potential density from subdiff as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_wdian_temp_subdiff_on_nrho = register_diag_field ('ocean_model',                            &
     'wdian_temp_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
     'temp related dianeutral mass transport due to subdiff as binned to neutral density layers',&
     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_tform_temp_subdiff_on_nrho = register_diag_field ('ocean_model',                       &
     'tform_temp_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
     'temp related watermass transform due to subdiff as binned to neutral density layers', &
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.



  id_neut_salt_subdiff = register_diag_field ('ocean_model', 'neut_salt_subdiff',    &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'salt related update of locally referenced potential density from subdiff param',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_subdiff > 0) compute_watermass_diff_diag = .true.

  id_wdian_salt_subdiff = register_diag_field ('ocean_model', 'wdian_salt_subdiff',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'salt related dianeutral mass transport due to subdiff closure',               &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_subdiff > 0) compute_watermass_diff_diag = .true.

  id_tform_salt_subdiff = register_diag_field ('ocean_model', 'tform_salt_subdiff', &
    Grd%tracer_axes(1:3), Time%model_time,                                          &
    'salt related watermass transform due to subdiff on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_subdiff > 0) compute_watermass_diff_diag = .true.

  id_neut_salt_subdiff_on_nrho = register_diag_field ('ocean_model',                                        &
    'neut_salt_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                &
    'salt related update of locally ref potential density from subdiff as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_wdian_salt_subdiff_on_nrho = register_diag_field ('ocean_model',                            &
     'wdian_salt_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
     'salt related dianeutral mass transport due to subdiff as binned to neutral density layers',&
     'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.

  id_tform_salt_subdiff_on_nrho = register_diag_field ('ocean_model',                      &
     'tform_salt_subdiff_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,             &
     'salt related watermass transform due to subdiff as binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_subdiff_on_nrho > 0) compute_watermass_diff_diag = .true.


  if(compute_watermass_diff_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_submesoscale_mod w/ compute_watermass_diff_diag=.true.'  
  endif 



end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from submesoscale on watermass transformation.  
! </DESCRIPTION>
!
subroutine watermass_diag(Time, T_prog, Dens)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ocean_density_type),       intent(in) :: Dens

  integer :: i,j,k,tau
  real,  dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_submesoscale (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return 

  tau=Time%tau 

  ! rho related diagnostics 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) &
                        +Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_submeso, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_submeso, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_submeso, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_submeso_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_submeso_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_submeso_on_nrho, wrk4)

  ! sea level tendency 
  if(id_eta_tend_submeso_tend > 0 .or. id_eta_tend_submeso_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_submeso_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_submeso_tend_glob, eta_tend, cellarea_r)
  endif


  ! temp related diagnostics 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) 
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_submeso, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_submeso, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_submeso, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_submeso_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_submeso_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_submeso_on_nrho, wrk4)

  ! salinity related diagnostics 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_submeso, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_submeso, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_submeso, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_submeso_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_submeso_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_submeso_on_nrho, wrk4)
 
end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"



!#######################################################################
! <SUBROUTINE NAME="watermass_diag_diffusion">
!
! <DESCRIPTION>
! Diagnose effects from submesoscale horizontal diffusion 
! on watermass transformation.  
! </DESCRIPTION>
!
subroutine watermass_diag_diffusion(Time, T_prog, Dens)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ocean_density_type),       intent(in) :: Dens

  integer :: i,j,k,tau
  real,  dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_submesoscale (watermass_diag_diffusion): module needs initialization ')
  endif 

  if(.not. compute_watermass_diff_diag) return 

  tau=Time%tau 

  ! rho related diagnostics 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) &
                        +Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_subdiff, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_subdiff, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_subdiff, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_subdiff_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_subdiff_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_subdiff_on_nrho, wrk4)

  ! sea level tendency 
  if(id_eta_tend_subdiff_tend > 0 .or. id_eta_tend_subdiff_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_subdiff_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_subdiff_tend_glob, eta_tend, cellarea_r)
  endif

  ! temp related diagnostics 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) 
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_subdiff, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_subdiff, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_subdiff, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_subdiff_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_subdiff_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_subdiff_on_nrho, wrk4)

  ! salinity related diagnostics 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_subdiff, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_subdiff, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_subdiff, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_subdiff_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_subdiff_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_subdiff_on_nrho, wrk4)

 
end subroutine watermass_diag_diffusion
! </SUBROUTINE>  NAME="watermass_diag_diffusion"


end module ocean_submesoscale_mod
