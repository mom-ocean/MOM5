module ocean_nphysics_util_mod
#define COMP isc:iec,jsc:jec
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<REVIEWER EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</REVIEWER>
!
!<OVERVIEW>
! Utilities for neutral physics, including the code to compute 
! space-time dependent diffusivities.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Utilities for neutral physics, including the code to compute 
! space-time dependent diffusivities and many diagnostics.  
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
! G. Danabasoglu and J. C. McWilliams
! Sensitivity of the global ocean circulation to 
! parameterizations of mesoscale tracer transports
! Journal of Climate (1995) vol 8 pages 2967--2987 
! </REFERENCE>
!
! <REFERENCE>
! Held and Larichev
! A scaling theory for horizontally homogeneous baroclinically 
! unstable flow on a beta plane
! Journal of Atmospheric Sciences (1996) vol 53 pages 946-952
! </REFERENCE>
!
! <REFERENCE>
! M. Visbeck, J.C. Marshall, T. Haine and M. Spall
! Specification of eddy transfer coefficients in coarse resolution ocean
! circulation models
! Journal of Physical Oceanography (1997) vol 27 pages 381--402
! </REFERENCE>
!
! <REFERENCE>
! D. Ferreira, J. Marshall, and P. Heimbach, 
! Estimating eddy stresses by fitting dynamics to observations 
! using a residual-mean ocean circulation omdel and its adjoint. 
! Journal of Physical Oceanography (2005) vol 35 pages 1891-1910.
! </REFERENCE>
!
! <REFERENCE>
! K. Eden, Eddy length scales in the North Atlantic, 2007,
! Preprint. 
! </REFERENCE>
!
! <REFERENCE>
! K. Eden and R. Greatbatch, 2008: Towards a mesoscale eddy closure,
! Ocean Modelling, vol. 20, pages 223-239
! </REFERENCE>
!
! <NOTE>
! Diffusivities can be determined in a number of manners
!
! TIME INDEPENDENT
!
! Various methods are available for specifying a time 
! independent diffusivity, either globally uniform or 
! with selections of spatial dependence.  
!
! TIME DEPENDENT (as a function of the flow)
!
! Various methods are available for determining the 
! diffusivity that changes in time according to the 
! properties of the fluid.  There are various means 
! for specifying the length and time scales needed
! to set the diffusivity. 
!
! LENGTH SCALES 
!
! 1. First baroclinic Rossby radius (estimated as per Chelton etal).  
! Equatorial Rossby radius is used within 5deg of the equator.
! 
! 2. Width of the baroclinic zone, as done in the Hadley Centre 
! model and documented in the MOM3 Manual.
!
! 3. Specified length scale set independent of the flow. 
!
! TIME SCALE 
!
! When using either of the above for the length scale, 
! the time scale is determined by the Eady growth rate.    
!
! COMBINED LENGTH/TIME SCALE 
! 
! Another option, used in the GFDL CM2.X coupled climate models,
! is to set the diffusivity proportional to the depth averaged
! absolute value of the horizontal density gradient.  
! </NOTE> 
!
! </INFO>
!
!
!<NAMELIST NAME="ocean_nphysics_util_nml">
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!
!  <DATA NAME="nphysics_util_zero_init" TYPE="logical">
!  For Time%init=.true. and wishing to ensure starting with a clean 
!  suite of nphysics_util fields, even if ocean_neutral.res.nc exists.  
!  </DATA> 
!
!  <DATA NAME="epsln_drhodz" UNITS="kg/m^3"  TYPE="real">
!  For computing drhodz used in slope calculation.
!  We must keep drhodz < 0 in order to maintain integrity of the 
!  quasi-Stokes streamfunction as well as computation of buoyancy frequency.  
!  Default epsln_drhodz=1e-30.
!  </DATA> 
!  <DATA NAME="drhodz_mom4p1" TYPE="logical">
!  For computing the vertical deriviative of locally referenced
!  potrho as in the preferred MOM algorithm rather than the 
!  earlier mom4p0 approach. Default drhodz_mom4p1=.true.
!  </DATA> 
!  <DATA NAME="drhodz_smooth_horz" TYPE="logical">
!  For horizontal laplacian smoothing the vertical derivative 
!  of density prior to its use in computing the neutral slopes. 
!  This smoothing helps to produce regularized slopes.  
!  Note that this option breaks the integrity of the triads
!  and is thus NOT generally recommended.  
!  Default drhodz_smooth_horz=.false.
!  </DATA> 
!  <DATA NAME="drhodz_smooth_vert" TYPE="logical">
!  For vertical 1-2-1 smoothing the vertical derivative of 
!  density prior to its use in computing the neutral slopes.
!  This smoothing helps to produce regularized slopes.  
!  Note that this option breaks the integrity of the triads
!  and is thus NOT generally recommended.  
!  Default drhodz_smooth_vert=.false.
!  </DATA> 
!
!  <DATA NAME="vel_micom_smooth" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in smoothing of drhodzb. 
!  Default vel_micom_smooth=0.2. 
!  </DATA>
!  <DATA NAME="num_121_passes" TYPE="integer">
!  For number of 1-2-1 passes through to smooth drhodz or 
!  eady_rate in vertical. Default num_121_passes=1. 
!  </DATA>
!
!  <DATA NAME="aredi" TYPE="real">
!  Neutral diffusivity used for experiments using a constant diffusivity. 
!  </DATA> 
!  <DATA NAME="agm" TYPE="real">
!  GM-skew diffusivity used for experiments using a constant diffusivity. 
!  </DATA> 
!  <DATA NAME="aredi_equal_agm" TYPE="logical">
!  Will set aredi_array=agm_array, over-riding any other specification 
!  of aredi_array. Default aredi_equal_agm=.true. 
!  </DATA> 
!
!  <DATA NAME="smax" TYPE="real">
!  Value of the maximum neutral direction slope above which the neutral fluxes are 
!  either tapered to zero or saturated.  Typical value is smax=0.01 or smaller. 
!  </DATA> 
!  <DATA NAME="swidth" TYPE="real">
!  Width in slope over which use tanh with dm_taper scheme to taper fluxes in 
!  steep sloped regions. Typical value swidth=0.1*smax
!  </DATA> 
!
!  <DATA NAME="smax_grad_gamma_scalar" TYPE="real">
!  For calculation of gradients of scalars along a neutral direction, then 
!  when abs(slope) > smax_grad_gamma_scalar, will compute the gradient using 
!  only the vertical scalar gradient, since the slopes are so large they are 
!  effectively infinite. 
!  Default smax_grad_gamma_scalar=.01 
!  </DATA> 

!  <DATA NAME="neutral_horz_mix_bdy" TYPE="logical">
!  If .true., then use a horizontal diffusivity in the neutral boundary layer. 
!  </DATA> 
!  <DATA NAME="vel_micom_bdy" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM horizontal diffusivity 
!  within the neutral boundary layer. 
!  </DATA> 
!  <DATA NAME="ah_bdy" UNITS="m^2/sec" TYPE="real">
!  Constant horizontal diffusivity for the boundary layer.  
!  Default ah_bdy=0.0.
!  </DATA> 
!
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the GM-skew diffusivity is set according to a velocity scale 
!  times the grid spacing. 
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!
!  <DATA NAME="agm_lat_zones" TYPE="logical">
!  If true, will set agm_array as constant within two latitudinal zones.  
!  The idea is that one may wish to use a larger agm in the ACC than 
!  elsewhere. 
!  </DATA> 
!  <DATA NAME="agm_lat_zones_boundary" TYPE="real">
!  Boundary between agm in the south and north zones. 
!  </DATA> 
!  <DATA NAME="agm_lat_zones_ratio" TYPE="real">
!  Ratio between the large agm used in the southern latitudinal zone
!  to that used in the north.  
!  agm_array(north) = agm
!  agm_array(south) = agm*agm_lat_zones_ratio
!  </DATA> 
!
!  <DATA NAME="bryan_lewis_aredi" TYPE="logical">
!  Set bryan_lewis_aredi=.true. when want to have aredi a function of depth
!  according to the Bryan and Lewis (1979) profile. Maintained for legacy 
!  purposes, and not recommended for new models.                                 
!  </DATA>
!  <DATA NAME="ahs" TYPE="real">
!  ahs = adjustable parameter at the surface for bryan_lewis_aredi   
!  </DATA>
!  <DATA NAME="ahb" TYPE="real">
!  ahb = adjustable parameter at the bottom for bryan_lewis_aredi   
!  </DATA>
!
!  <DATA NAME="agm_read_restart" TYPE="logical">
!  For those cases with agm_closure=.false. where we wish to read in 
!  the agm_array from restart files and keep the value from the restart.
!  This approach allows us to read in a spatially dependent agm_array 
!  that may have been computed from another integration, but to leave
!  the coefficient static in time.  
!  Default agm_read_restart=.false.
!  </DATA> 
!
!  <DATA NAME="agm_closure" TYPE="logical">
!  If .true. then will compute the GM-skew diffusivity as a function of the flow.
!  The length scale is determined by the Rossby radius and the time scale is 
!  determined by the Eady growth rate.  Diffusivities are depth independent.  
!  </DATA> 
!  <DATA NAME="agm_closure_max" UNITS="m^2/sec" TYPE="real">
!  Maximum GM diffusivity allowed when using agm_closure=.true. 
!  </DATA> 
!  <DATA NAME="agm_closure_min" UNITS="m^2/sec" TYPE="real">
!  Minimum GM diffusivity allowed when using agm_closure=.true. 
!  </DATA> 
!  <DATA NAME="agm_closure_scaling" UNITS="dimensionless" TYPE="logical">
!  Dimensionless tuning parameter for computing flow dependent diffusivities. 
!  </DATA> 
!  <DATA NAME="agm_closure_upper_depth" UNITS="m" TYPE="real">
!  Upper depth where start the depth integration to compute the Eady 
!  growth rate and/or baroclinicity.
!  </DATA> 
!  <DATA NAME="agm_closure_lower_depth" UNITS="m" TYPE="real">
!  Deeper depth where finish the depth integration to compute the Eady 
!  growth rate and/or baroclinicity.
!  </DATA> 
!
!  <DATA NAME="agm_closure_n2_scale" TYPE="logical">
!  For computing the agm coefficient using a 3-dimensional 
!  scaling by (N/Nref)^2, with N=buoyancy frequency and 
!  Nref the buoyancy frequency at the base of the neutral 
!  blayer. Default agm_closure_n2_scale=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_n2_scale_coeff" UNITS="m^2/s" TYPE="real">
!  Coefficient setting the scale for the diffusivity computed from 
!  agm_closure_n2_scale. 
!  Default agm_closure_n2_scale_coeff=1e3. 
!  </DATA> 
!  <DATA NAME="agm_closure_n2_scale_nref_cst" TYPE="logical">
!  For taking the reference buoyancy frequency as agm_closure_buoy_freq
!  for the (N/Nref)^2 scaling.  
!  Default agm_closure_n2_scale_nref_cst=.false.
!  </DATA> 
!
!  <DATA NAME="agm_closure_baroclinic" TYPE="logical">
!  For computing the agm coefficient using only the vertically
!  averaged magnitude of the horizontal density gradient 
!  (i.e., the "baroclinicity").
!  </DATA> 
!  <DATA NAME="agm_closure_buoy_freq" UNITS="sec^-1" TYPE="real">
!  For computing the agm coefficient using only the vertically
!  averaged horizontal density gradient, we need to specify a 
!  buoyancy frequency,  which is taken to be fixed over all space-time.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_cap" TYPE="logical">
!  For setting a maximum length scale for the agm_closure calculation.
!  Default agm_closure_length_cap=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_length_max" UNITS="metre" TYPE="real">
!  Maximum length scale used for computing agm_closure.  
!  Default agm_closure_length_max=50e3.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_fixed" TYPE="logical">
!  Use fixed length scale for computing agm_closure diffusivity 
!  </DATA> 
!  <DATA NAME="agm_closure_length" UNITS="meter" TYPE="real">
!  Fixed length scale for use with agm_closure_fixed_length
!  </DATA>
!  
!  <DATA NAME="agm_closure_length_rossby" TYPE="logical">
!  For computing the agm_closure length scale according to Rossby radius. 
!  </DATA> 
!  <DATA NAME="rossby_radius_max" UNITS="meter" TYPE="real">
!  Maximum Rossby radius used for agm_closure_length_rossby and 
!  the neutral_sine_taper. Default = 100e3 m.
!  </DATA> 
!  <DATA NAME="rossby_radius_min" UNITS="meter" TYPE="real">
!  Minimum Rossby Radius used for agm_closure_length_rossby and 
!  the neutral_sine_taper. Default = 15e3 m.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_bczone" TYPE="logical">
!  For computing the agm_closure length scale according to radius of baroclinic zone. 
!  </DATA> 
!  <DATA NAME="bczone_max_pts" TYPE="integer">
!  Max number of horizontal grid points for use in computing the baroclinic zone radius.  
!  </DATA> 
!  <DATA NAME="agm_closure_bczone_crit_rate" UNITS="sec^-1" TYPE="real">
!  Critical growth rate for determining width of the baroclinic zone. 
!  </DATA> 
!  <DATA NAME="agm_closure_growth_scale" UNITS="dimensionless" TYPE="real">
!  Dimensionless scaling used to set a maximum for agm_growth. 
!  </DATA> 
!
!  <DATA NAME="agm_closure_eden_greatbatch" TYPE="logical">
!  For computing the agm_closure length scale according to minimum 
!  of the Rhines scale and the Rossby radius, and using 3d Eady 
!  growth rate. 
!  </DATA> 
!  <DATA NAME="agm_closure_eden_gamma" TYPE="real" UNITS="dimensionless">
!  For use in regularizing the growth rate used in the eden/greatbatch approach.
!  Default agm_closure_eden_gamma=200. Setting to zero removes the regularization.  
!  </DATA> 
!  <DATA NAME="agm_closure_eden_length_const" TYPE="logical">
!  To set the length scale for agm_closure_eden_greatbatch to constant. 
!  Default agm_closure_eden_length_const=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eden_length" TYPE="real" UNITS="metre">
!  Length scale for use with agm_closure_eden_length_const=.true.
!  Default agm_closure_eden_length=10e3. 
!  </DATA> 
!
!  <DATA NAME="agm_closure_eady_smooth_vert" TYPE="logical">
!  For vertical 1-2-1 smoothing the eady_rate 
!  Default agm_closure_eady_smooth=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eady_smooth_horz" TYPE="logical">
!  For horizontal Laplacian smoothing of eady growth rate. 
!  Default agm_closure_eady_smooth_horz=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eady_cap" TYPE="logical">
!  For capping the eady growth rate to avoid huge values.  
!  Default agm_closure_eady_cap=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eady_ave_mixed" TYPE="logical">
!  To set the Eady growth rate to its average within mixed layer region.  
!  This is used to avoid spuriously large values which often appear just 
!  in the upper regions of the ocean mixed layer.  
!  Default agm_closure_eady_ave_mixed=.false.
!  </DATA> 
!
!  <DATA NAME="agm_smooth_space" TYPE="logical">
!  For smoothing the agm diffusivity in space when nonconstant diffusivity used. 
!  Default is agm_smooth_space=.false.
!  </DATA> 
!  <DATA NAME="agm_smooth_time" TYPE="logical">
!  For smoothing the agm diffusivity in time when nonconstant diffusivity used. 
!  Default is agm_smooth_time=.false.
!  </DATA> 
!  <DATA NAME="agm_damping_time" UNITS="days" TYPE="real">
!  The damping time used for time smoothing agm_array.
!  Default agm_damping_time=10days.  
!  </DATA> 
!
!  <DATA NAME="agm_closure_grid_scaling" TYPE="logical">
!  For an overall scaling of the agm coefficient, according to
!  the relative resolution of the grid and deformation radius. 
!  Default is agm_closure_grid_scaling=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_grid_scaling_power" TYPE="real">
!  Power used to scale the agm_closure diffusivity. 
!  Default is agm_closure_grid_scaling_power=2.0
!  </DATA> 
!  <DATA NAME="aredi_diffusivity_grid_scaling" TYPE="logical">
!  For an overall scaling of the aredi coefficient, according to
!  the relative resolution of the grid and deformation radius.
!  This option is used only when aredi_equal_agm=.false. 
!  Default is aredi_diffusivity_grid_scaling=.false.
!  </DATA> 
!
!  <DATA NAME="epsln_drhodz_diagnostics" UNITS="kg/m4" TYPE="real">
!  For drhodz used in calculation of dianeutral velocity component
!  from cabbeling and thermobaricity.
!  Default epsln_drhodz_diagnostics=1e-7.
!  </DATA>
!  <DATA NAME="wdianeutral_smooth" TYPE="logical">
!  For smoothing the diagnosed dianeutral velocity component using a 
!  horizontal 1-2-1 smoother. Default is wdianeutral_smooth=.true.
!  </DATA> 
!
!  <DATA NAME="smooth_eta_tend_gm90" TYPE="logical">
!  For smoothing the diagnosed contribution to steric sea level 
!  time tendency associated with the GM90 scheme. 
!  Default is smooth_eta_tend_gm90=.false.
!  </DATA> 
!</NAMELIST>

use constants_mod,    only: epsln, pi
use diag_manager_mod, only: register_diag_field, register_static_field, need_data
use fms_mod,          only: FATAL, file_exist
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: restart_file_type
use mpp_mod,          only: input_nml_file, mpp_pe, mpp_min, mpp_error, stdout, stdlog
use mpp_mod,          only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use mpp_domains_mod,  only: mpp_update_domains, EUPDATE, NUPDATE
use time_manager_mod, only: time_type, increment_time

use ocean_domains_mod,           only: get_local_indices, set_ocean_domain
use ocean_density_mod,           only: calc_cabbeling_thermobaricity
use ocean_operators_mod,         only: FDX_T, FDY_T, FMX, FMY, LAP_T, S2D
use ocean_parameters_mod,        only: missing_value, onehalf, onefourth
use ocean_parameters_mod,        only: rho0, rho0r, grav
use ocean_tracer_diag_mod,       only: calc_mixed_layer_depth
use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,             only: ocean_prog_tracer_type, ocean_density_type
use ocean_types_mod,             only: ocean_time_type, ocean_thickness_type, ocean_time_steps_type
use ocean_types_mod,             only: tracer_3d_0_nk_type, tracer_3d_1_nk_type 
use ocean_util_mod,              only: write_timestamp, write_chksum_3d, write_chksum_2d, write_chksum_header
use ocean_util_mod,              only: diagnose_2d, diagnose_3d, diagnose_sum
use ocean_util_mod,              only: write_line, write_note, write_warning, write_chksum_3d
use ocean_tracer_util_mod,       only: diagnose_3d_rho
use ocean_workspace_mod,         only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk6
use ocean_workspace_mod,         only: wrk1_v, wrk2_v, wrk3_v, wrk1_2d, wrk2_2d, wrk3_2d


implicit none

public ocean_nphysics_util_init
public ocean_nphysics_coeff_init
public ocean_nphysics_coeff_end
public tracer_derivs 
public neutral_slopes
public compute_eady_rate
public compute_baroclinicity
public compute_rossby_radius
public compute_bczone_radius
public compute_diffusivity
public ocean_nphysics_util_restart
public transport_on_rho_gm
public transport_on_nrho_gm
public transport_on_theta_gm
public cabbeling_thermob_tendency
public compute_eta_tend_gm90
public watermass_diag_init
public watermass_diag
public watermass_diag_ndiffuse
public watermass_diag_sdiffuse

private calc_gradx_gamma_scalar
private calc_grady_gamma_scalar
private pressure_derivs 
private cabbeling_thermob_init

private 

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: BCzone_domain   

! to write some output to all processors 
integer :: unit=6

! clock ids
integer :: id_clock_tracer_derivs
integer :: id_clock_neutral_slopes
integer :: id_clock_compute_eady_rate
integer :: id_clock_compute_baroclinicity
integer :: id_clock_compute_rossby_radius
integer :: id_clock_compute_bczone_radius
integer :: id_clock_compute_diffusivity 
integer :: id_clock_cabbeling_thermob
integer :: id_clock_eta_tend_gm90
integer :: id_transport_on_nrho_gm
integer :: id_transport_on_rho_gm
integer :: id_transport_on_theta_gm

! for diagnostics manager 
integer :: id_slope31                  =-1
integer :: id_slope32                  =-1
integer :: id_aredi                    =-1
integer :: id_agm                      =-1
integer :: id_aredi_3d                 =-1
integer :: id_agm_3d                   =-1
integer :: id_eady_mld                 =-1
integer :: id_eady_rate                =-1
integer :: id_eady_rate_zave           =-1
integer :: id_rossby                   =-1
integer :: id_rossby_radius            =-1
integer :: id_rossby_equator           =-1
integer :: id_rossby_nonequator        =-1
integer :: id_bczone                   =-1
integer :: id_gw_speed                 =-1
integer :: id_sqrt2betaCr              =-1
integer :: id_growth_rate_baroclinic   =-1
integer :: id_baroclinicity            =-1
integer :: id_agm_growth_rate          =-1
integer :: id_agm_length               =-1
integer :: id_rhines_length            =-1
integer :: id_agm_qg                   =-1
integer :: id_N2slope                  =-1
integer :: id_N2_nblayer_base          =-1
integer :: id_N2_for_agm               =-1
integer :: id_coriolis_param           =-1
integer :: id_beta_param               =-1
integer :: id_grid_length              =-1
integer :: id_agm_grid_scaling         =-1
integer :: id_aredi_grid_scaling       =-1
integer :: id_ksurf_blayer             =-1
integer :: id_cabbeling_tend           =-1
integer :: id_thermobaric_tend         =-1
integer :: id_cabbeling_speed          =-1
integer :: id_thermobaric_speed        =-1
integer :: id_cabbeling_tend_intz      =-1
integer :: id_thermobaric_tend_intz    =-1
integer :: id_eta_tend_gm90            =-1
integer :: id_eta_tend_gm90_glob       =-1

integer :: id_cabbeling_tend_intz_glob =-1
integer :: id_thermobaric_tend_intz_glob =-1

integer :: id_neut_rho_cabbeling         =-1
integer :: id_neut_rho_cabbeling_on_nrho =-1
integer :: id_wdian_cabbeling            =-1
integer :: id_wdian_cabbeling_on_nrho    =-1
integer :: id_tform_rho_cabbel_on_nrho   =-1

integer :: id_neut_rho_thermob         =-1
integer :: id_neut_rho_thermob_on_nrho =-1
integer :: id_wdian_thermob            =-1
integer :: id_wdian_thermob_on_nrho    =-1
integer :: id_tform_rho_thermb_on_nrho =-1

! for watermass transformation diagnostics
integer :: id_neut_rho_nphysics
integer :: id_wdian_rho_nphysics

integer :: id_neut_rho_ndiff
integer :: id_wdian_rho_ndiff
integer :: id_tform_rho_ndiff
integer :: id_neut_rho_ndiff_on_nrho
integer :: id_wdian_rho_ndiff_on_nrho
integer :: id_tform_rho_ndiff_on_nrho

integer :: id_eta_tend_ndiff_tend     
integer :: id_eta_tend_ndiff_tend_glob

integer :: id_neut_temp_ndiff
integer :: id_wdian_temp_ndiff
integer :: id_tform_temp_ndiff
integer :: id_neut_temp_ndiff_on_nrho
integer :: id_wdian_temp_ndiff_on_nrho
integer :: id_tform_temp_ndiff_on_nrho

integer :: id_neut_salt_ndiff
integer :: id_wdian_salt_ndiff
integer :: id_tform_salt_ndiff
integer :: id_neut_salt_ndiff_on_nrho
integer :: id_wdian_salt_ndiff_on_nrho
integer :: id_tform_salt_ndiff_on_nrho

integer :: id_neut_rho_gm
integer :: id_wdian_rho_gm
integer :: id_tform_rho_gm
integer :: id_neut_rho_gm_on_nrho
integer :: id_wdian_rho_gm_on_nrho
integer :: id_tform_rho_gm_on_nrho

integer :: id_eta_tend_gm_tend     
integer :: id_eta_tend_gm_tend_glob

integer :: id_neut_temp_gm
integer :: id_wdian_temp_gm
integer :: id_tform_temp_gm
integer :: id_neut_temp_gm_on_nrho
integer :: id_wdian_temp_gm_on_nrho
integer :: id_tform_temp_gm_on_nrho

integer :: id_neut_salt_gm
integer :: id_wdian_salt_gm
integer :: id_tform_salt_gm
integer :: id_neut_salt_gm_on_nrho
integer :: id_wdian_salt_gm_on_nrho
integer :: id_tform_salt_gm_on_nrho

! for transport_on_nrho_gm
integer :: id_tx_trans_nrho_gm=-1
integer :: id_ty_trans_nrho_gm=-1

! for transport_on_rho_gm
integer :: id_tx_trans_rho_gm=-1
integer :: id_ty_trans_rho_gm=-1

! for transport_on_theta_gm
integer :: id_tx_trans_theta_gm=-1
integer :: id_ty_trans_theta_gm=-1

! for restart
type(restart_file_type), save :: Nphysics_util_restart

integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1
integer :: neutralrho_nk

logical :: used

! internally set for computing watermass diagnstics
logical :: compute_watermass_diag = .false. 

! 2D array of micom gm diffusivities (m^2/sec)
real, dimension(:,:), allocatable :: agm_micom 

! 2D array of deformation radius scaling
real, dimension(:,:), allocatable :: agm_grid_scaling
real, dimension(:,:), allocatable :: aredi_grid_scaling

! 3D array of gm length scale (metre)
real, dimension(:,:,:), allocatable :: agm_length 

! 3D array of agm_array for local purposes (m^2/sec)
real, dimension(:,:,:), allocatable :: agm_array_local

! 3D array of aredi_array for local purposes (m^2/sec)
real, dimension(:,:,:), allocatable :: aredi_array_local

! 3D array of gm growth rate scale (sec^-1)
real, dimension(:,:,:), allocatable :: agm_growth_rate

! 2D array of gm growth rate max (sec^-1)
real, dimension(:,:), allocatable :: agm_growth_rate_max

! grid length array (metre)
real, dimension(:,:),  allocatable :: grid_length

! absolute value of the Coriolis parameter (sec^-1)
real, dimension(:,:),  allocatable :: coriolis_param 

! beta = d(Coriolis)/dy (m^-1 sec^-1)
real, dimension(:,:),  allocatable :: beta_param  

! for determining baroclinic zone radius (metre) 
real, dimension(:,:), allocatable :: bczone_dxt  
real, dimension(:,:), allocatable :: bczone_dyt  
real, dimension(:,:), allocatable :: bczone_rate 

! for Eady growth rate and baroclinicity calculation 
real, dimension(:,:),  allocatable :: count_x 
real, dimension(:,:),  allocatable :: count_y 

! Eady growth rate (sec^-1)
real, dimension(:,:,:),  allocatable :: eady_rate
real, dimension(:,:),    allocatable :: eady_rate_zave

! 1st baroclinic gravity wave speed (m/s)
real, dimension(:,:),  allocatable :: gravity_wave_speed

! for Eden and Greatbatch near equator (m)
real, dimension(:,:),  allocatable :: sqrt2betaCr

!vertically ave abs(horiz density gradient) (kg/m^4)
real, dimension(:,:),  allocatable :: baroclinicity

! (m^2/sec) for smoothing
real, dimension(:,:),  allocatable :: smooth_lap

! cabbeling and thermobaricity tendencies 
real, dimension(:,:,:), allocatable :: cabbeling_param
real, dimension(:,:,:), allocatable :: thermobaric_param
real, dimension(:,:,:), allocatable :: gradx_gamma_temp
real, dimension(:,:,:), allocatable :: grady_gamma_temp
real, dimension(:,:,:), allocatable :: gradx_gamma_press
real, dimension(:,:,:), allocatable :: grady_gamma_press
real, dimension(:,:,:), allocatable :: deriv_x_press
real, dimension(:,:,:), allocatable :: deriv_y_press
real, dimension(:,:,:), allocatable :: deriv_z_press

! gm contribution to steric sea level tendency 
real, dimension(:,:,:), allocatable :: taper_fcn
real, dimension(:,:,:), allocatable :: slopex
real, dimension(:,:,:), allocatable :: agm_rho_slopex
real, dimension(:,:,:), allocatable :: dslopex_dx
real, dimension(:,:,:), allocatable :: dslopex_dz
real, dimension(:,:,:), allocatable :: gradx_gamma_slopex
real, dimension(:,:,:), allocatable :: slopey
real, dimension(:,:,:), allocatable :: agm_rho_slopey
real, dimension(:,:,:), allocatable :: dslopey_dy
real, dimension(:,:,:), allocatable :: dslopey_dz
real, dimension(:,:,:), allocatable :: grady_gamma_slopey

! some useful constants 
real :: grav_r
real :: gravrho0r
real :: gravrho0r_buoyr
real :: gamma_damp  
real :: dtime
real :: cellarea_r

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 


!*************************************
! nml settings 
!
! for how to compute drhodz prior to computing slopes 
logical :: drhodz_mom4p1      =.true.
logical :: drhodz_smooth_horz =.false.
logical :: drhodz_smooth_vert =.false.
integer :: num_121_passes     = 1
real    :: vel_micom_smooth   = 0.2   ! m/sec
real    :: epsln_drhodz       = 1e-30    

! globally constant diffusivities 
real    :: aredi = 1.0e3 ! constant neutral diffusion tracer diffusivity (m^2/sec)
real    :: agm   = 1.0e3 ! constant gent-mcwilliams skew-diffusion diffusivity (m^2/sec)

! will set aredi_array=agm_array 
logical :: aredi_equal_agm =.true.

! maximum neutral direction slope allowed before tapering fluxes
real    :: smax=0.01

! slope width over which fluxes tapered using tanh function using dm_taper scheme 
real    :: swidth = 0.05*.01

! if true, apply horizontal diffusivity in neutral boundary layer for nphysicsA module
logical :: neutral_horz_mix_bdy =.false. 
real    :: vel_micom_bdy        = 0.0  ! velocity scale (m/s) for horizontal diffusivity w/i boundary 
real    :: ah_bdy               = 0.0  ! constant horiz diffusivity in neutral bdy

! time independent, vertically dependent neutral diffusivity according to Bryan-Lewis. 
! maintained in MOM only for legacy purposes.  not recommended for new models. 
logical :: bryan_lewis_aredi = .false.  
real    :: ahs               = 0.0      
real    :: ahb               = 0.0      

! for setting agm according to latitude zones 
logical :: agm_lat_bands          = .false. 
real    :: agm_lat_bands_boundary = -999.   ! boundary between agm in the south and north zones
real    :: agm_lat_bands_ratio    = 1.0     ! ratio agm(south)/agm(north)

! set diffusivity according to local horizontal area 
! of grid and a specified velocity scale (m/s)
logical :: tracer_mix_micom =.false. 
real    :: vel_micom        = 0.0    

! determining how to read in restart diffusivity
logical :: agm_read_restart=.false.

! for setting diffusivity according to flow properties 
logical :: agm_closure         =.false. 
real    :: agm_closure_scaling = 2.0    ! dimensionless tuning parameter 
real    :: agm_closure_max     = 2.e3   ! maximum diffusivity allowed when agm_closure=.true.
real    :: agm_closure_min     = 2.e2   ! minimum diffusivity allowed when agm_closure=.true.

! for fixed length scale (metres) set by agm_closure_length
logical :: agm_closure_length_fixed =.false. 
real    :: agm_closure_length       = 50.e3  ! metres

! for length scale set according to estimate of first baroclinic Rossby radius
logical :: agm_closure_length_rossby =.false. 
real    :: rossby_radius_max=100e3  ! metres 
real    :: rossby_radius_min=15e3   ! metres 

! for length scale set according to radius of baroclinic zone
logical :: agm_closure_length_bczone    =.false. 
integer :: bczone_max_pts               =10      ! max # points searched for determining baroclinic zone width 
real    :: agm_closure_bczone_crit_rate =1.4e-6  ! critical growth rate for determining baroclinic zone (sec^-1)

! for capping the agm_length scale 
logical :: agm_closure_length_cap = .false.
real    :: agm_closure_length_max = 50e3

! for scaling diffusivity by (N/Nref)^2
logical :: agm_closure_n2_scale=.false. 
logical :: agm_closure_n2_scale_nref_cst=.false.
real    :: agm_closure_n2_scale_coeff=1e3

! for diffusity according to Eden and Greatbatch (2008) and Eden(2007)
logical :: agm_closure_eden_greatbatch  =.false. 
logical :: agm_closure_eden_length_const=.false.
real    :: agm_closure_eden_gamma  = 200.0
real    :: agm_closure_eden_length = 10e3

! for diffusivity set proportional to vertically averaged horiz density gradient  
logical :: agm_closure_baroclinic = .true. 
! sec^-1 buoyancy frequency for use with agm_closure_baroclinic
real    :: agm_closure_buoy_freq  = 0.004  

! depths over which compute the Eady growth and/or horizontal density gradient 
real    :: agm_closure_upper_depth  =100.0  ! metre
real    :: agm_closure_lower_depth  =2000.0 ! metre 

! for smoothing and capping eady growth rate
logical :: agm_closure_eady_smooth_vert=.false.
logical :: agm_closure_eady_smooth_horz=.false.
logical :: agm_closure_eady_cap        =.false.
logical :: agm_closure_eady_ave_mixed  =.false.

! dimensionless number to set maximum value of agm_growth. 
! agm_closure_growth_scale =0.5 yields agm_growth max=coriolis parameter.
real    :: agm_closure_growth_scale =0.5    

! for smoothing agm_array in space and time
logical :: agm_smooth_space=.false.
logical :: agm_smooth_time =.false.
real    :: agm_damping_time=10.0    

! for overall scaling of the agm coefficient 
logical :: agm_closure_grid_scaling=.false.
real    :: agm_closure_grid_scaling_power=2.0

! for scaling aredi_array according to size of grid and Rossby radius  
logical :: aredi_diffusivity_grid_scaling=.false.

! for Tim%init=.true. and ocean_neutral.res.nc exists, but still wish
! to start from initial fields.
logical :: nphysics_util_zero_init=.false.

! for computing tendencies from thermobaricity and cabbeling 
logical :: diagnose_cabbeling_thermob=.false.  ! internally set  
logical :: wdianeutral_smooth=.true.
real    :: epsln_drhodz_diagnostics=1e-7
real    :: smax_grad_gamma_scalar=.01
real    :: swidthr                     ! inverse swidth  
real    :: smax_swidthr                ! useful combination of terms 

! for computing contribution to steric sea level evolution from GM
logical :: smooth_eta_tend_gm90   =.false.
logical :: diagnose_eta_tend_gm90 =.false.  ! internally set  

! for debugging 
logical :: debug_this_module = .false.

!*************************************

character(len=128) :: version=&
     '$Id: ocean_nphysics_util.F90,v 20.0 2013/12/14 00:14:52 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


logical :: module_is_initialized = .FALSE.

namelist /ocean_nphysics_util_nml/ debug_this_module, nphysics_util_zero_init,      &
          smax, swidth, epsln_drhodz, drhodz_mom4p1,                                &
          drhodz_smooth_horz, drhodz_smooth_vert, num_121_passes,                   &
          aredi, agm, aredi_equal_agm, tracer_mix_micom, vel_micom,                 &
          bryan_lewis_aredi, ahs, ahb,                                              &
          neutral_horz_mix_bdy, vel_micom_bdy, ah_bdy,                              &
          agm_lat_bands, agm_lat_bands_boundary, agm_lat_bands_ratio,               &  
          rossby_radius_max, rossby_radius_min,  agm_read_restart,                  &
          agm_closure, agm_closure_scaling, agm_closure_max, agm_closure_min,       &
          agm_closure_growth_scale, agm_closure_length_fixed, agm_closure_length,   &
          agm_closure_length_rossby, agm_closure_length_bczone,                     &
          bczone_max_pts, agm_closure_bczone_crit_rate,                             &
          agm_closure_eden_greatbatch, agm_closure_eden_gamma,                      &
          agm_closure_eden_length_const, agm_closure_eden_length,                   &
          agm_closure_eady_smooth_vert, agm_closure_eady_smooth_horz,               &
          agm_closure_eady_ave_mixed, agm_closure_eady_cap,                         &
          agm_closure_baroclinic, agm_closure_buoy_freq,                            &   
          agm_closure_upper_depth, agm_closure_lower_depth,                         &
          agm_closure_length_cap, agm_closure_length_max,                           &
          agm_smooth_space, vel_micom_smooth, agm_smooth_time, agm_damping_time,    &
          agm_closure_grid_scaling, agm_closure_grid_scaling_power,                 &
          aredi_diffusivity_grid_scaling,                                           &
          agm_closure_n2_scale, agm_closure_n2_scale_coeff,                         &
          agm_closure_n2_scale_nref_cst,                                            &
          smax_grad_gamma_scalar,                                                   &
          epsln_drhodz_diagnostics, wdianeutral_smooth,                             &
          smooth_eta_tend_gm90 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_util_init">
!
! <DESCRIPTION>
! Initialize the utility module for neutral physics.
! </DESCRIPTION>
!
subroutine ocean_nphysics_util_init(Grid, Domain, Time, Time_steps, Dens, T_prog, &
           agm_closure_lower_dept, agm_closure_upper_dept, agm_closure_buoy_frq,  &
           smx, swidt, cmip_units, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  real,                         intent(inout)        :: agm_closure_lower_dept
  real,                         intent(inout)        :: agm_closure_upper_dept
  real,                         intent(inout)        :: agm_closure_buoy_frq
  real,                         intent(inout)        :: smx 
  real,                         intent(inout)        :: swidt
  logical,                      intent(in)           :: cmip_units
  logical,                      intent(in), optional :: debug

  integer :: ioun, io_status, ierr  
  integer :: i,j,n
  integer :: num_schemes 
  real    :: param 
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysics_util_mod: already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  Dom => Domain
  Grd => Grid

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_nphysics_util_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_nphysics_util_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysics_util_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_nphysics_util_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysics_util_nml)  
  write (stdlogunit,ocean_nphysics_util_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
     call write_note(FILENAME,&
     'running ocean_nphysics_util_mod with debug_this_module=.true.')
  endif

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  ! for setting up temp and salt tracers
  num_prog_tracers = size(T_prog(:))
  index_temp=-1;index_salt=-1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then 
     call mpp_error(FATAL, &
     '==>Error: temp and/or salt not identified in call to ocean_nphysics_util_init')
  endif 

  ! some useful constants 
  dtime           = Time_steps%dtime_t 
  grav_r          = 1.0/grav 
  gravrho0r       = grav*rho0r                       !for buoyancy frequency calculation  
  gravrho0r_buoyr = gravrho0r/agm_closure_buoy_freq  !for agm_closure_baroclinic
  gamma_damp      = dtime/(86400.0*agm_damping_time) !for damping time dependent agm_array 
  swidthr         = 1.0/(swidth + epsln)             !for slope taper function for cabbeling/thermob diagnostic
  smax_swidthr    = smax_grad_gamma_scalar*swidthr   !for slope taper function for cabbeling/thermob diagnostic
  neutralrho_nk   = size(Dens%neutralrho_ref(:))
  cellarea_r      = 1.0/(epsln + Grd%tcellsurf)


  if(cmip_units) then
      transport_convert=1.0
      transport_dims   = 'kg/s'
  else
      transport_convert=1.0e-9 
      transport_dims   = 'Sv (10^9 kg/s)'
  endif

  ! for sending back to ocean_nphysics_mod 
  agm_closure_lower_dept = agm_closure_lower_depth
  agm_closure_upper_dept = agm_closure_upper_depth
  agm_closure_buoy_frq   = agm_closure_buoy_freq
  smx                    = smax
  swidt                  = swidth 

  ! Coriolis parameter and beta parameter 
  allocate (coriolis_param(isd:ied,jsd:jed))
  allocate (beta_param(isd:ied,jsd:jed))
  coriolis_param(:,:) = 2.0*7.292e-5*abs(sin(Grd%phit(:,:)))
  beta_param(:,:)     = 2.28e-11*abs(cos(Grd%phit(:,:)))

  allocate (gravity_wave_speed(isd:ied,jsd:jed))
  allocate (sqrt2betaCr(isd:ied,jsd:jed))
  allocate (eady_rate(isd:ied,jsd:jed,nk))
  allocate (eady_rate_zave(isd:ied,jsd:jed))
  allocate (baroclinicity(isd:ied,jsd:jed))
  sqrt2betaCr(:,:)        = 0.0
  gravity_wave_speed(:,:) = 0.0
  eady_rate(:,:,:)        = 0.0
  eady_rate_zave(:,:)     = 0.0
  baroclinicity(:,:)      = 0.0

  ! for computing the baroclinic zone radius using Hadley Centre search algorithm 
  call set_ocean_domain(BCzone_domain,Grd,xhalo=bczone_max_pts,yhalo=bczone_max_pts,name='bczone',maskmap=Dom%maskmap)
  allocate (bczone_rate(isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  allocate (bczone_dxt (isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  allocate (bczone_dyt (isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  bczone_dxt=0.0 ; bczone_dyt=0.0 ; bczone_rate=0.0
  bczone_dxt(isc:iec,jsc:jec) = Grd%dxt(isc:iec,jsc:jec)
  bczone_dyt(isc:iec,jsc:jec) = Grd%dyt(isc:iec,jsc:jec)
  call mpp_update_domains (bczone_dxt(:,:), BCzone_domain%domain2d)
  call mpp_update_domains (bczone_dyt(:,:), BCzone_domain%domain2d)

  allocate (grid_length(isd:ied,jsd:jed))
  grid_length(:,:) = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))

  allocate (agm_micom(isd:ied,jsd:jed))
  if(tracer_mix_micom) then 
       agm_micom(:,:) = vel_micom*grid_length(:,:)
  else 
       agm_micom(:,:) = 0.0
  endif 

  allocate (agm_length(isd:ied,jsd:jed,nk))
  agm_length(:,:,:) = agm_closure_length*Grd%tmask(:,:,:)

  allocate (agm_grid_scaling(isd:ied,jsd:jed))
  agm_grid_scaling(:,:) = Grd%tmask(:,:,1)

  allocate (aredi_grid_scaling(isd:ied,jsd:jed))
  aredi_grid_scaling(:,:) = Grd%tmask(:,:,1)

  allocate (agm_growth_rate(isd:ied,jsd:jed,nk))
  agm_growth_rate(:,:,:) = 0.0

  allocate (agm_array_local(isd:ied,jsd:jed,nk))
  agm_array_local(:,:,:) = agm*Grd%tmask(:,:,:)

  allocate (aredi_array_local(isd:ied,jsd:jed,nk))
  aredi_array_local(:,:,:) = aredi*Grd%tmask(:,:,:)


  !max agm_growth_rate is 2.0*agm_closure_growth_scale*param
  !if agm_closure_growth_scale=0.5, then agm_growth_rate is <= param
  allocate (agm_growth_rate_max(isd:ied,jsd:jed))
  do j=jsd,jed
     do i=isd,ied
        param = max(coriolis_param(i,j), sqrt2betaCr(i,j))
        agm_growth_rate_max(i,j) = agm_closure_growth_scale*param
     enddo
  enddo

  if(Time_steps%aidif /= 1.0 .and. aredi /= 0) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_nphysics_util_mod: stability requires aidif=1.0 to handle K33 implicitly.')
  endif

  if(agm_closure) then

      num_schemes = 0 
      if(agm_closure_length_fixed) then 
          num_schemes = num_schemes + 1
          call write_note(FILENAME,&
          'Computing 2d flow-dependent tracer diffusivity with agm_closure_length_fixed.')
      endif
      if(agm_closure_length_rossby) then 
          num_schemes = num_schemes + 1
          call write_note(FILENAME,&
          'Computing 2d flow-dependent tracer diffusivity with agm_closure_length_rossby.')
      endif
      if(agm_closure_length_bczone) then 
          num_schemes = num_schemes + 1
          call write_note(FILENAME,&
          'Computing 2d flow-dependent tracer diffusivity with agm_closure_length_bczone.')
      endif
      if(agm_closure_eden_greatbatch) then 
          num_schemes = num_schemes + 1
          call write_note(FILENAME,&
          'Computing 3d flow-dependent tracer diffusivity with agm_closure_eden_greatbatch.')
      endif
      if(agm_closure_baroclinic) then 
          num_schemes = num_schemes + 1
          call write_note(FILENAME,&
          'Computing 2d flow-dependent tracer diffusivity with agm_closure_baroclinic.')
      endif
      if(agm_closure_n2_scale) then 
          num_schemes = num_schemes + 1
          call write_note(FILENAME,&
          'Computing 3d flow-dependent tracer diffusivity with agm_closure_n2_scale.')
      endif

      if(num_schemes == 0) then 
          call mpp_error(FATAL, &
          '==>Error: with agm_closure=.true., must choose one of the "agm_closure" methods')
      endif
      if(num_schemes > 1) then 
          call mpp_error(FATAL, &
          '==>Error: with agm_closure=.true., can choose only one of the available agm_closure methods')
      endif

      write(stdoutunit,'(a,e10.5)') &
      '        The maximum allowable diffusivity (m^2/s) is given by ',agm_closure_max
      write(stdoutunit,'(a,e10.5)') &
      '        The minimum allowable diffusivity (m^2/s) is given by ',agm_closure_min

      write(stdoutunit,'(a,e10.5,a,e10.5)') &
      '  Depths (m) between which compute eady growth and baroclinicity = ', &
         agm_closure_upper_depth, ' ',agm_closure_lower_depth


      if(agm_closure_n2_scale) then 
         call write_note(FILENAME,&
         'Diffusivity will be computed as agm = coeff*(N/Nref)^2.')
      endif 

      if(agm_closure_baroclinic) then 
         call write_note(FILENAME,&
         'Length and time scales set by vertically averaged baroclinicity |grad(rho)|,')
          write(stdoutunit,'(a,e10.5)') &
               '        as well as the constant buoyancy freq(sec^-1) = ',agm_closure_buoy_freq
          write(stdoutunit,'(a,e10.5)') &
               '        and the constant length scale (m) = ',agm_closure_length
      else 
         call write_note(FILENAME,&
         'Eady growth rate gives inverse time scale.')
          if(agm_closure_length_fixed) then 
              write(stdoutunit,'(a,e12.5,a)') &
                   '        Length scale set by nml parameter agm_closure_length =',agm_closure_length,' metre'
          elseif(agm_closure_length_rossby) then 
             call write_line('First baroclinic Rossby radius gives the length scale.')
          elseif(agm_closure_length_bczone) then 
             call write_line('Radius of baroclinic zone gives the length scale.')
          elseif(agm_closure_eden_greatbatch) then
              if(agm_closure_eden_length_const) then 
                  write(stdoutunit,'(a,e12.5,a)') &
                       '        Take constant length scale of ',agm_closure_eden_length, 'metre'
              else 
                 call write_line('Min(rossby,rhines) gives the length scale.')
              endif
          endif
      endif

  endif   ! endif for agm_closure 

  ! for smoothing
  allocate (smooth_lap(isd:ied,jsd:jed))
  smooth_lap(:,:) = vel_micom_smooth*grid_length(:,:)

  ! for Eady growth rate and baroclinicity calculation 
  allocate (count_x(isd:ied,jsd:jed))
  allocate (count_y(isd:ied,jsd:jed))
  count_x(:,:) = 0.0
  count_y(:,:) = 0.0
  do i=isc,iec
    do j=jsc,jec
      count_x(i,j) =  &
      min(Grd%tmask(i+1,j,1) + Grd%tmask(i-1,j,1), 1.0/(2.0*(Grd%tmask(i+1,j,1) + Grd%tmask(i-1,j,1) + epsln)))
      count_y(i,j) =  &
      min(Grd%tmask(i,j+1,1) + Grd%tmask(i,j-1,1), 1.0/(2.0*(Grd%tmask(i,j+1,1) + Grd%tmask(i,j-1,1) + epsln)))
    enddo
  enddo

  ! initialize clock ids 
  id_clock_tracer_derivs         = mpp_clock_id('(Ocean neutral: tracer derivs)'  ,grain=CLOCK_ROUTINE)
  id_clock_neutral_slopes        = mpp_clock_id('(Ocean neutral: slopes)'         ,grain=CLOCK_ROUTINE)
  id_clock_compute_eady_rate     = mpp_clock_id('(Ocean neutral: eady rate)'      ,grain=CLOCK_ROUTINE)
  id_clock_compute_baroclinicity = mpp_clock_id('(Ocean neutral: baroclinic)'     ,grain=CLOCK_ROUTINE)
  id_clock_compute_rossby_radius = mpp_clock_id('(Ocean neutral: rossby radius)'  ,grain=CLOCK_ROUTINE)
  id_clock_compute_bczone_radius = mpp_clock_id('(Ocean neutral: bc zone radius)' ,grain=CLOCK_ROUTINE)
  id_clock_compute_diffusivity   = mpp_clock_id('(Ocean neutral: diffusivity)'    ,grain=CLOCK_ROUTINE)
  id_clock_cabbeling_thermob     = mpp_clock_id('(Ocean neutral: cabbel/thermob)' ,grain=CLOCK_ROUTINE)
  id_clock_eta_tend_gm90         = mpp_clock_id('(Ocean neutral: eta_tend_gm90)'  ,grain=CLOCK_ROUTINE)
  id_transport_on_nrho_gm        = mpp_clock_id('(Ocean neutral: nrho-trans)'     ,grain=CLOCK_ROUTINE)
  id_transport_on_rho_gm         = mpp_clock_id('(Ocean neutral: rho-trans)'      ,grain=CLOCK_ROUTINE)
  id_transport_on_theta_gm       = mpp_clock_id('(Ocean neutral: theta-trans)'    ,grain=CLOCK_ROUTINE)


  ! for diagnostics manager 

  id_ksurf_blayer= -1            
  id_ksurf_blayer= register_diag_field ('ocean_model', 'ksurf_blayer',     &
                Grd%tracer_axes(1:2), Time%model_time,                     &
                'k-value at base of surf nblayer region', 'dimensionless', &
                missing_value=missing_value, range=(/-1.0,1.e1/))

  id_N2slope = -1
  id_N2slope = register_diag_field ('ocean_model', 'N2slope',               & 
            Grd%tracer_axes_wt(1:3), Time%model_time,                       &
            'Squared buoyancy frequency used in neutral slope calculation', &
            '1/s^2', missing_value=missing_value, range=(/-1e3,1.e10/))

  id_N2_for_agm = -1
  id_N2_for_agm = register_diag_field ('ocean_model', 'N2_for_agm',           & 
            Grd%tracer_axes(1:3), Time%model_time,                            &
            'Squared buoyancy frequency used for (N/Nref)^2 agm calculation', &
            '1/s^2', missing_value=missing_value, range=(/-1e3,1.e10/))

  id_N2_nblayer_base= -1
  id_N2_nblayer_base= register_diag_field ('ocean_model', 'N2_nblayer_base',& 
            Grd%tracer_axes(1:2), Time%model_time,                          &
            'Squared buoyancy frequency at base of nblayer',                &
            '1/s^2', missing_value=missing_value, range=(/-1e3,1.e10/))

  id_slope31 = -1
  id_slope31 = register_diag_field ('ocean_model', 'slope31',   &
               Grd%tracer_axes(1:3), Time%model_time,           &
               'neutral slope -(rho_x/rho_z)', 'dimensionless', &
               missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_slope32 = -1
  id_slope32 = register_diag_field ('ocean_model', 'slope32',   &
               Grd%tracer_axes(1:3), Time%model_time,           &
               'neutral slope -(rho_y/rho_z)', 'dimensionless', &
               missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_eady_mld = -1          
  id_eady_mld = register_diag_field ('ocean_model', 'eady_mld',        &
            Grd%tracer_axes(1:2), Time%model_time,                     &
            'Mixed layer depth inside of which compute ave Eady rate', &
            's^-1', missing_value=missing_value, range=(/-10.0,1.e10/))

  id_eady_rate_zave = -1          
  id_eady_rate_zave = register_diag_field ('ocean_model', 'eady_rate_zave',&
            Grd%tracer_axes(1:2), Time%model_time,                         &
            'Eady growth rate averaged over depth range for nphysics ',    &
            's^-1', missing_value=missing_value, range=(/-10.0,1.e10/))

  id_eady_rate = -1          
  id_eady_rate = register_diag_field ('ocean_model', 'eady_rate',&
            Grd%tracer_axes(1:3), Time%model_time,               &
            'Eady growth rate used in neutral physics ', 's^-1', &
            missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby = -1      
  id_rossby = register_diag_field ('ocean_model', 'rossby', &
              Grd%tracer_axes(1:2), Time%model_time,        &
              'Rossby radius used in neutral physics', 'm', &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby_radius = -1      
  id_rossby_radius = register_diag_field ('ocean_model', 'rossby_radius', &
              Grd%tracer_axes(1:2), Time%model_time,                      &
              'Rossby radius computed without min/max bounds', 'm',       &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby_equator = -1      
  id_rossby_equator = register_diag_field ('ocean_model', 'rossby_equator', &
              Grd%tracer_axes(1:2), Time%model_time,                        &
              'Equatorial Rossby radius', 'm',                              &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby_nonequator = -1      
  id_rossby_nonequator = register_diag_field ('ocean_model', 'rossby_nonequator', &
              Grd%tracer_axes(1:2), Time%model_time,                              &
              'Rossby radius outside equatorial region', 'm',                     &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_bczone = -1      
  id_bczone = register_diag_field ('ocean_model', 'bczone', &
              Grd%tracer_axes(1:2), Time%model_time,        &
              'radius of baroclinic zone', 'm',             &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_sqrt2betaCr = -1          
  id_sqrt2betaCr = register_diag_field ('ocean_model', 'sqrt2betaCr',         &
            Grd%tracer_axes(1:2), Time%model_time,                            &
            'sqrt(2 * 1st b/c wave speed *beta) in neutral physics ', '1/sec',&
            missing_value=missing_value, range=(/-1e0,1e10/))

  id_gw_speed = -1          
  id_gw_speed = register_diag_field ('ocean_model', 'gw_speed',  &
            Grd%tracer_axes(1:2), Time%model_time,               &
            'Gravity wave speed used in neutral physics ', 'm/s',&
            missing_value=missing_value, range=(/-1e2,1e2/))

  id_baroclinicity = -1          
  id_baroclinicity = register_diag_field ('ocean_model', 'baroclinicity',   &
                     Grd%tracer_axes(1:2), Time%model_time,                 &
                     'vertically averaged horz density gradient', 'kg/m^4', &
                     missing_value=missing_value, range=(/-1e10,1.e10/))

  id_growth_rate_baroclinic = -1          
  id_growth_rate_baroclinic = register_diag_field ('ocean_model', 'growth_rate_baroclinic', &
                              Grd%tracer_axes(1:3), Time%model_time,                        &
                              'growth rate using baroclinicity', 's^-1',                    &
                              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm_growth_rate = -1          
  id_agm_growth_rate = register_diag_field ('ocean_model', 'agm_growth_rate',&
                  Grd%tracer_axes(1:3), Time%model_time,                     &
                  'effective growth rate for agm', 's^-1',                   &
                  missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm_length = -1          
  id_agm_length = register_diag_field ('ocean_model', 'agm_length',  &
                  Grd%tracer_axes(1:3), Time%model_time,             &
                  'effective length scale for agm', 'm',             &
                  missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rhines_length = -1          
  id_rhines_length = register_diag_field ('ocean_model', 'rhines_length', &
                  Grd%tracer_axes(1:3), Time%model_time,                  &
                  'Rhines length approximated as eady_rate/beta', 'm',    &
                  missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm_qg = -1          
  id_agm_qg = register_diag_field ('ocean_model', 'agm_qg',  &
              Grd%tracer_axes(1:3), Time%model_time,         &
              'agm from QG theory', 'm^2/s',                 &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm = -1
  id_agm = register_diag_field ('ocean_model', 'agm',         &
           Grd%tracer_axes(1:2), Time%model_time,             &
           'GM diffusivity at surface', 'm^2/sec',            &
           missing_value=missing_value, range=(/-10.0,1.e10/),&
           standard_name='ocean_tracer_bolus_laplacian_diffusivity')

  id_agm_grid_scaling = register_diag_field ('ocean_model','agm_grid_scaling',          &
                        Grd%tracer_axes(1:2), Time%model_time,                          &
                       'Scaling of AGM according to Delta(s)^2/(Delta(s)^2 + Rossby^2)',&
                        'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))

  id_aredi_grid_scaling = register_diag_field ('ocean_model','aredi_grid_scaling',        &
                        Grd%tracer_axes(1:2), Time%model_time,                            &
                       'Scaling of Aredi according to Delta(s)^2/(Delta(s)^2 + Rossby^2)',&
                        'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))

  id_agm_3d = -1
  id_agm_3d = register_diag_field ('ocean_model', 'agm_3d',      &
              Grd%tracer_axes(1:3), Time%model_time,             &
              '3d GM diffusivity', 'm^2/sec',                    &
              missing_value=missing_value, range=(/-10.0,1.e10/),&
              standard_name='ocean_tracer_bolus_laplacian_diffusivity')

  id_aredi = -1            
  id_aredi = register_diag_field ('ocean_model', 'aredi',       &
             Grd%tracer_axes(1:2), Time%model_time,             &
             'neutral diffusivity at k=1', 'm^2/sec',           &
             missing_value=missing_value, range=(/-10.0,1.e20/),&
             standard_name='ocean_tracer_epineutral_laplacian_diffusivity')

  id_aredi_3d = -1            
  id_aredi_3d = register_diag_field ('ocean_model', 'aredi_3d', &
                Grd%tracer_axes(1:3), Time%model_time,          &
                '3d neutral diffusivity', 'm^2/sec',            &
                missing_value=missing_value, range=(/-10.0,1.e20/))

  id_tx_trans_nrho_gm = register_diag_field ('ocean_model','tx_trans_nrho_gm', Dens%neutralrho_axes_flux_x(1:3),&
                        Time%model_time, 'T-cell i-mass transport from GM on neutral rho',trim(transport_dims), &
                        missing_value=missing_value, range=(/-1e20,1e20/))
  id_ty_trans_nrho_gm = register_diag_field ('ocean_model','ty_trans_nrho_gm', Dens%neutralrho_axes_flux_y(1:3),&
                        Time%model_time, 'T-cell j-mass transport from GM on neutral rho',trim(transport_dims), &
                        missing_value=missing_value, range=(/-1e20,1e20/))

  id_tx_trans_rho_gm = register_diag_field ('ocean_model','tx_trans_rho_gm', Dens%potrho_axes_flux_x(1:3),&
                       Time%model_time, 'T-cell i-mass transport from GM on pot_rho',trim(transport_dims),&
                       missing_value=missing_value, range=(/-1e20,1e20/))
  id_ty_trans_rho_gm = register_diag_field ('ocean_model','ty_trans_rho_gm', Dens%potrho_axes_flux_y(1:3),&
                       Time%model_time, 'T-cell j-mass transport from GM on pot_rho',trim(transport_dims),&
                       missing_value=missing_value, range=(/-1e20,1e20/))

  id_tx_trans_theta_gm = register_diag_field ('ocean_model','tx_trans_theta_gm', Dens%theta_axes_flux_x(1:3),&
                       Time%model_time, 'T-cell i-mass transport from GM on theta',trim(transport_dims),     &
                       missing_value=missing_value, range=(/-1e20,1e20/))
  id_ty_trans_theta_gm = register_diag_field ('ocean_model','ty_trans_theta_gm', Dens%theta_axes_flux_y(1:3),&
                       Time%model_time, 'T-cell j-mass transport from GM on theta',trim(transport_dims),     &
                       missing_value=missing_value, range=(/-1e20,1e20/))

  ! static fields 
  id_coriolis_param = register_static_field ('ocean_model', 'coriolis_param', Grd%tracer_axes(1:2), &
                                        'Coriolis frequency on T-cell for nphysics', '1/s',         &
                                         missing_value=missing_value, range=(/-10.0,10.0/))
  call diagnose_2d(Time, Grd, id_coriolis_param, coriolis_param(:,:))

  id_beta_param = register_static_field ('ocean_model', 'beta_param', Grd%tracer_axes(1:2),&
                                        'Beta=df/dy on T-cell for nphysics', '1/s',        &
                                         missing_value=missing_value, range=(/-10.0,10.0/))
  call diagnose_2d(Time, Grd, id_beta_param, beta_param(:,:))

  id_grid_length = register_static_field ('ocean_model', 'grid_length', Grd%tracer_axes(1:2),   &
                                        'Grid length scale used for nphysics calculations', 'm',&
                                         missing_value=missing_value, range=(/-10.0,1e10/))
  call diagnose_2d(Time, Grd, id_grid_length, grid_length(:,:))


  ! initialize cabbeling and thermobaricity related diagnostics 
  call cabbeling_thermob_init(Time, Dens)


  ! eta tendency calculations  
  id_eta_tend_gm90= -1          
  id_eta_tend_gm90= register_diag_field ('ocean_model','eta_tend_gm90',                         &
    Grd%tracer_axes(1:2), Time%model_time,                                                      &
    'analytic (and less accurate) form of non-Bouss steric sea level tendency from GM90', 'm/s',&
    missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_gm90 > 0) diagnose_eta_tend_gm90=.true.

  id_eta_tend_gm90_glob= -1          
  id_eta_tend_gm90_glob= register_diag_field ('ocean_model', 'eta_tend_gm90_glob',                  &
    Time%model_time,                                                                                &
   'global mean analytic (and less accurate) form of non-bouss steric sea level tendency from GM90',&
   'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_gm90_glob > 0) diagnose_eta_tend_gm90=.true.

  if(diagnose_eta_tend_gm90) then 
     call write_note(FILENAME, &
      'diagnose_eta_tend_gm90==.true., diagnosing global sea level tendency from GM90.')
     allocate (taper_fcn(isd:ied,jsd:jed,nk))
     allocate (slopex(isd:ied,jsd:jed,nk))
     allocate (agm_rho_slopex(isd:ied,jsd:jed,nk))
     allocate (dslopex_dx(isd:ied,jsd:jed,nk))
     allocate (dslopex_dz(isd:ied,jsd:jed,nk))
     allocate (gradx_gamma_slopex(isd:ied,jsd:jed,nk))
     allocate (slopey(isd:ied,jsd:jed,nk))
     allocate (agm_rho_slopey(isd:ied,jsd:jed,nk))
     allocate (dslopey_dy(isd:ied,jsd:jed,nk))
     allocate (dslopey_dz(isd:ied,jsd:jed,nk))
     allocate (grady_gamma_slopey(isd:ied,jsd:jed,nk))
     taper_fcn(:,:,:)           = 0.0
     slopex(:,:,:)              = 0.0
     agm_rho_slopex(:,:,:)      = 0.0
     dslopex_dx(:,:,:)          = 0.0
     dslopex_dz(:,:,:)          = 0.0
     gradx_gamma_slopex(:,:,:)  = 0.0
     slopey(:,:,:)              = 0.0
     agm_rho_slopey(:,:,:)      = 0.0
     dslopey_dy(:,:,:)          = 0.0
     dslopey_dz(:,:,:)          = 0.0
     grady_gamma_slopey(:,:,:)  = 0.0
  endif 


end subroutine ocean_nphysics_util_init
! </SUBROUTINE>  NAME="ocean_nphysics_util_init"



!#######################################################################
! <SUBROUTINE NAME="cabbeling_thermob_init">
!
! <DESCRIPTION>
! Initialize the cabbeling and thermobaricity related diagnostic 
! fields and register the diagnostics to the diag manager. 
! </DESCRIPTION>
!

subroutine cabbeling_thermob_init(Time, Dens)
  type(ocean_time_type),   intent(in) :: Time
  type(ocean_density_type),intent(in) :: Dens


  id_cabbeling_tend= -1          
  id_cabbeling_tend = register_diag_field ('ocean_model', 'cabbeling_tend',&
                  Grd%tracer_axes(1:3), Time%model_time,                   &
                  '(1/rho)*d(rho)/dt from cabbeling', 's^-1',              &
                  missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_cabbeling_tend > 0) diagnose_cabbeling_thermob =.true. 

  id_cabbeling_speed= -1          
  id_cabbeling_speed= register_diag_field ('ocean_model', 'cabbeling_speed',&
                  Grd%tracer_axes(1:3), Time%model_time,                    &
                  '(dz/rho)*d(rho)/dt from cabbeling', 'm/s',               &
                  missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_cabbeling_tend > 0) diagnose_cabbeling_thermob =.true. 

  id_cabbeling_tend_intz= -1          
  id_cabbeling_tend_intz= register_diag_field ('ocean_model', 'cabbeling_tend_intz',&
                  Grd%tracer_axes(1:2), Time%model_time,                            &
                  'vertical sum of (dz/rho)*d(rho)/dt from cabbeling', 'm/s',       &
                  missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_cabbeling_tend_intz > 0) diagnose_cabbeling_thermob =.true. 

  id_cabbeling_tend_intz_glob= -1          
  id_cabbeling_tend_intz_glob= register_diag_field ('ocean_model',                &
                  'cabbeling_tend_intz_glob', Time%model_time,                    &
                  'global mean vertical sum of (dz/rho)*d(rho)/dt from cabbeling',&
                  'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_cabbeling_tend_intz_glob > 0) diagnose_cabbeling_thermob =.true. 

  id_thermobaric_tend= -1          
  id_thermobaric_tend= register_diag_field ('ocean_model', 'thermobaric_tend',&
                  Grd%tracer_axes(1:3), Time%model_time,                      &
                  '(1/rho)*d(rho)/dt from thermobaricity', 's^-1',            &
                  missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_thermobaric_tend > 0) diagnose_cabbeling_thermob =.true. 

  id_thermobaric_speed= -1          
  id_thermobaric_speed= register_diag_field ('ocean_model', 'thermobaric_speed',&
                  Grd%tracer_axes(1:3), Time%model_time,                        &
                  '(dz/rho)*d(rho)/dt from thermobaricity', 'm/s',              &
                  missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_thermobaric_speed > 0) diagnose_cabbeling_thermob =.true. 

  id_thermobaric_tend_intz= -1          
  id_thermobaric_tend_intz= register_diag_field ('ocean_model', 'thermobaric_tend_intz',&
                  Grd%tracer_axes(1:2), Time%model_time,                                &
                  'vertical sum of (dz/rho)*d(rho)/dt from thermobaricity', 'm/s',      &  
                  missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_thermobaric_tend_intz > 0) diagnose_cabbeling_thermob =.true. 

  id_thermobaric_tend_intz_glob= -1          
  id_thermobaric_tend_intz_glob= register_diag_field ('ocean_model',                   &
                  'thermobaric_tend_intz_glob', Time%model_time,                       &
                  'global mean vertical sum of (dz/rho)*d(rho)/dt from thermobaricity',&
                  'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_thermobaric_tend_intz_glob > 0) diagnose_cabbeling_thermob =.true. 

  id_neut_rho_cabbeling= -1          
  id_neut_rho_cabbeling= register_diag_field ('ocean_model', 'neut_rho_cabbeling', &
     Grd%tracer_axes(1:3), Time%model_time,                                        &
    'update of locally referenced potential density from cabbeling', '(kg/m^3)/s', &  
     missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_neut_rho_cabbeling > 0) diagnose_cabbeling_thermob =.true. 

  id_wdian_cabbeling= -1          
  id_wdian_cabbeling= register_diag_field ('ocean_model', 'wdian_cabbeling',&
    Grd%tracer_axes(1:3), Time%model_time,                                  &
   'dianeutral mass transport due to cabbeling', 'kg/s',                    &  
    missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_wdian_cabbeling > 0) diagnose_cabbeling_thermob =.true. 

  id_neut_rho_thermob= -1          
  id_neut_rho_thermob= register_diag_field ('ocean_model', 'neut_rho_thermob',   &
     Grd%tracer_axes(1:3), Time%model_time,                                      &
     'update of locally referenced potential density from thermob', '(kg/m^3)/s',&  
     missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_neut_rho_thermob > 0) diagnose_cabbeling_thermob =.true. 

  id_wdian_thermob= -1          
  id_wdian_thermob= register_diag_field ('ocean_model', 'wdian_thermob',&
     Grd%tracer_axes(1:3), Time%model_time,                             &
     'dianeutral mass transport due to thermobaricity', 'kg/s',         &  
      missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_wdian_thermob > 0) diagnose_cabbeling_thermob =.true. 

  id_neut_rho_cabbeling_on_nrho= -1          
  id_neut_rho_cabbeling_on_nrho= register_diag_field ('ocean_model', 'neut_rho_cabbeling_on_nrho',  &
     Dens%neutralrho_axes(1:3), Time%model_time,                                                    &
     'update of locally referenced potential density due to cabbeling as binned to neutral density',&
     '(kg/m^3)/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_neut_rho_cabbeling_on_nrho > 0) diagnose_cabbeling_thermob =.true. 

  id_wdian_cabbeling_on_nrho= -1          
  id_wdian_cabbeling_on_nrho= register_diag_field ('ocean_model', 'wdian_cabbeling_on_nrho', &
     Dens%neutralrho_axes(1:3), Time%model_time,                                             &
     'dianeutral mass transport due to cabbeling as binned to neutral density', 'kg/s',      &  
     missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_wdian_cabbeling_on_nrho > 0) diagnose_cabbeling_thermob =.true. 

  id_neut_rho_thermob_on_nrho= -1          
  id_neut_rho_thermob_on_nrho= register_diag_field ('ocean_model', 'neut_rho_thermob_on_nrho',           &
     Dens%neutralrho_axes(1:3), Time%model_time,                                                         &
     'update of locally referenced potential density due to thermobaricity as binned to neutral density',&
     '(kg/m^3)/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_neut_rho_thermob_on_nrho > 0) diagnose_cabbeling_thermob =.true. 

  id_wdian_thermob_on_nrho= -1          
  id_wdian_thermob_on_nrho= register_diag_field ('ocean_model', 'wdian_thermob_on_nrho',    &
     Dens%neutralrho_axes(1:3), Time%model_time,                                            &
     'dianeutral mass transport due to thermobaricity as binned to neutral density', 'kg/s',&  
     missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_wdian_thermob_on_nrho > 0) diagnose_cabbeling_thermob =.true. 

  id_tform_rho_cabbel_on_nrho =-1
  id_tform_rho_cabbel_on_nrho = register_diag_field ('ocean_model', 'tform_rho_cabbel_on_nrho',&
   Dens%neutralrho_axes(1:3), Time%model_time,                                                 &
   'watermass transform due to cabbeling as binned to neutral rho layers',                     &
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_cabbel_on_nrho > 0) diagnose_cabbeling_thermob =.true. 

  id_tform_rho_thermb_on_nrho =-1
  id_tform_rho_thermb_on_nrho = register_diag_field ('ocean_model', 'tform_rho_thermb_on_nrho',&
   Dens%neutralrho_axes(1:3), Time%model_time,                                                 &
   'watermass transform due to thermobaricity as binned to neutral rho layers',                &
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_thermb_on_nrho > 0) diagnose_cabbeling_thermob =.true. 


  if(diagnose_cabbeling_thermob) then 
     call write_note(FILENAME, &
     'diagnose_cabbeling_thermob==.true., diagnosing tendencies from cabbeling and thermobaricity.')
     allocate (cabbeling_param(isd:ied,jsd:jed,nk))
     allocate (thermobaric_param(isd:ied,jsd:jed,nk))
     allocate (gradx_gamma_temp(isd:ied,jsd:jed,nk))
     allocate (grady_gamma_temp(isd:ied,jsd:jed,nk))
     allocate (gradx_gamma_press(isd:ied,jsd:jed,nk))
     allocate (grady_gamma_press(isd:ied,jsd:jed,nk))
     allocate (deriv_x_press(isd:ied,jsd:jed,nk))
     allocate (deriv_y_press(isd:ied,jsd:jed,nk))
     allocate (deriv_z_press(isd:ied,jsd:jed,nk))
     cabbeling_param(:,:,:)   = 0.0
     thermobaric_param(:,:,:) = 0.0
     gradx_gamma_temp(:,:,:)  = 0.0
     grady_gamma_temp(:,:,:)  = 0.0
     gradx_gamma_press(:,:,:) = 0.0
     grady_gamma_press(:,:,:) = 0.0
     deriv_x_press(:,:,:)     = 0.0
     deriv_y_press(:,:,:)     = 0.0
     deriv_z_press(:,:,:)     = 0.0
  endif 


end subroutine cabbeling_thermob_init
! </SUBROUTINE>  NAME="cabbeling_thermob_init"



!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_coeff_init">
!
! <DESCRIPTION>
! Initialize the diffusivities used in neutral physics.
! Need to initialize them after the ocean_nphysics_util_init routine,
! since need to have the domain parameters known for passing the 
! array size information into the ocean_nphysics_coeff_init routine. 
! </DESCRIPTION>
!

subroutine ocean_nphysics_coeff_init(Time, Thickness, rossby_radius, rossby_radius_raw, &
                      bczone_radius, agm_array, aredi_array, ah_array)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  real, dimension(isd:,jsd:),   intent(inout) :: rossby_radius
  real, dimension(isd:,jsd:),   intent(inout) :: rossby_radius_raw
  real, dimension(isd:,jsd:),   intent(inout) :: bczone_radius
  real, dimension(isd:,jsd:,:), intent(inout) :: agm_array
  real, dimension(isd:,jsd:,:), intent(inout) :: aredi_array
  real, dimension(isd:,jsd:),   intent(inout) :: ah_array

  integer :: i,j,k
  integer :: i_delta, j_delta, k_delta
  integer :: id_restart
  real    :: ft, deltaX, deltaY, delta, delta_iso, delta_iso0, delta_min, A_max, A_max0
  character(len=64) :: file_name

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! initialization based on constant coefficients 
  agm_array(:,:,:)         = agm*Grd%tmask(:,:,:)
  agm_array_local(:,:,:)   = agm*Grd%tmask(:,:,:)
  aredi_array(:,:,:)       = aredi*Grd%tmask(:,:,:)
  aredi_array_local(:,:,:) = aredi*Grd%tmask(:,:,:)

  ! horizontal diffusivity in neutral boundary layer region 
  ah_array(:,:) = 0.0 
  if(neutral_horz_mix_bdy) then
     call write_note(FILENAME,&
      'Adding horizontal diffusivity in neutral boundary.')
     call write_line('This method is implemented only for the case with neutral_physics_simple=.false.')
      if(vel_micom_bdy > 0.0) then
         ah_array(:,:) = vel_micom_bdy*grid_length(:,:)
      else 
         ah_array(:,:) = ah_bdy
      endif 
      do j=jsc,jec
         write (stdoutunit,'(a,i4,a,e14.7,a)') &
         ' ah_array in neutral bdy layer at (isc,',j,',1)= ',ah_array(isc,j),' m^2/s'
      enddo
  endif


  ! Bryan-Lewis profile for aredi_array
  ! this profile is implemented for legacy purposes  
  if(bryan_lewis_aredi) then
     call write_note(FILENAME,&
      'Using Bryan-Lewis depth profile for the Redi diffusivity.')
     call write_line('The Bryan-Lewis profile has been tested ONLY with agm==0.')
     call write_line('Vertical dependence to the neutral diffusivity has not been implemented')
     call write_line('according to a physical theory giving a vertical structure depending')
     call write_line('on the flow field. More theoretical work is required.')
     do k=1,nk
        aredi_array(:,:,k) = (ahb + (ahs - ahb)*exp(-Grd%zt(k)/500.0))
     enddo
     do k=1,nk
        write (stdoutunit,'(a,i3,e16.8)') '  k, diffusivity = ', k, aredi_array(isc,jsc,k)
     enddo
  endif

  ! Set diffusivity according to agm_lat_bands
  if(agm_lat_bands) then
     call write_note(FILENAME,&
      'Setting agm_array according to agm_lat_bands.')
      write(stdoutunit,'(a,e10.5)') &
      '      The ratio agm(south)/agm(north) is given by ',agm_lat_bands_boundary
      write(stdoutunit,'(a,e10.5)') &
      '      with the latitude separating the bands given by ',agm_lat_bands_ratio
      if(agm_lat_bands_boundary <= -90.0) then 
          write(stdoutunit,'(1x,a)') &
           '      Since agm_lat_bands_boundary <= -90, will default to globally constant agm value' 
      endif
      agm_array_local(:,:,:) = agm*Grd%tmask(:,:,:)
      if(agm_lat_bands_boundary > -90.) then 
          do j=jsc-Dom%yhalo,jec+Dom%yhalo
             do i=isc-Dom%xhalo,iec+Dom%xhalo
                if(Grd%yt(i,j) <= agm_lat_bands_boundary) then
                  agm_array_local(i,j,:) = agm*agm_lat_bands_ratio*Grd%tmask(i,j,:)
                endif 
             enddo
          enddo
      endif
  endif

  ! grid-scale dependent diffusivity suggested by that commonly used in MICOM
  ! vel_micom (m/s) sets velocity scale.
  ! space scale is set by grid size. 
  if(tracer_mix_micom) then
      do k=1,nk
        agm_array_local(:,:,k) = agm_micom(:,:)
      enddo
      do j=jsc,jec
         write (stdoutunit,'(a,i4,a,e14.7,a)') &
         ' Micom agm_array at (isc,',j,',1) = ',agm_array_local(isc,j,1),' m^2/s'
      enddo
      if(aredi_equal_agm) then
          aredi_array(:,:,:) = agm_array_local(:,:,:) 
          call write_note(FILENAME,&
          'aredi_array = agm_array as given by Micom grid scale dependent diffusivity')
      else
         call write_note(FILENAME,&
          'Since (agm/=aredi), aredi_array=aredi')
      endif
  endif

  ! to have the halos filled 
  agm_array(:,:,:) = agm_array_local(:,:,:)

  ! make mandatory=.false. to facilitate backward compatibility with older restart files. 
  file_name = 'ocean_neutral.res.nc'
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'agm_array', agm_array, &
       domain=Dom%domain2d, mandatory=.false.) 
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'aredi_array', aredi_array, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'rossby_radius', rossby_radius, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'rossby_radius_raw', rossby_radius_raw, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'bczone_radius', bczone_radius, &
       domain=Dom%domain2d, mandatory=.false.)

  if(Time%init .and. nphysics_util_zero_init) then 
     call write_note(FILENAME,&
     'Starting ocean_nphysics_util fields from raw initialization.')
  elseif(.NOT. file_exist('INPUT/ocean_neutral.res.nc')) then
     if (.NOT. Time%init) call mpp_error(FATAL,'Expecting file INPUT/ocean_neutral.res.nc to exist.&
             &This file was not found and Time%init=.false.')
  else
      call restore_state(Nphysics_util_restart)

      call mpp_update_domains(agm_array,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'agm_array being read from restart.'
      call write_chksum_3d('checksum start agm_array', agm_array(COMP,:)*Grd%tmask(COMP,:))

      call mpp_update_domains(aredi_array,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'aredi_array being read from restart.'
      call write_chksum_3d('checksum start aredi_array', aredi_array(COMP,:)*Grd%tmask(COMP,:))

      call mpp_update_domains(rossby_radius,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'rossby_radius being read from restart.'
      call write_chksum_2d('checksum start rossby_radius', rossby_radius(COMP)*Grd%tmask(COMP,1))

      call mpp_update_domains(rossby_radius_raw,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'rossby_radius_raw being read from restart.'
      call write_chksum_2d('checksum start rossby_radius_raw', rossby_radius_raw(COMP)*Grd%tmask(COMP,1))

      call mpp_update_domains(bczone_radius,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'bczone_radius being read from restart.'
      call write_chksum_2d('checksum start bczone_radius', bczone_radius(COMP)*Grd%tmask(COMP,1))

  endif 

  ! for those cases with agm_array determined by 
  ! time invariant initialization methods.  
  if(.not. agm_closure) then
      if(.not. agm_read_restart) then
         call write_note(FILENAME,&
         'agm_closure=.false. and agm_read_restart=.false.')
         call write_line('=> agm_array set to static profiles.')
         agm_array(:,:,:) = agm_array_local(:,:,:)
      else
         call write_note(FILENAME,&
          'agm_closure=.false. and agm_read_restart=.true.')
         call write_line('=> agm_array set to restart values.')
      endif
  endif

  ! allow aredi value to be read from namelist rather than from a restart
  if(.not. agm_read_restart .and. .not. aredi_equal_agm) then
      write(stdoutunit,'(1x,a)') &
     'aredi_equal_agm=.false. and agm_read_restart=.false. => aredi_array set to static profiles.'
      aredi_array(:,:,:) = aredi_array_local(:,:,:)
  endif

  if(aredi_equal_agm) then 
     call write_note(FILENAME,&
     'aredi_equal_agm=.true. will force aredi_array to equal agm_array')
     aredi_array(:,:,:) = agm_array(:,:,:)
  else 
     call write_note(FILENAME,&
     'aredi_equal_agm=.false. allows aredi_array to differ from agm_array')
  endif 

  ! compute maximum stable neutral slope available for neutral diffusion
  i_delta = isc; j_delta = jsc; k_delta = 1
  delta_iso = 1e30
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,k) /= 0.0) then 
          ft     = 0.5/(aredi_array(i,j,k)*dtime + epsln)
          deltaX = Grd%dxt(i,j)*Thickness%dzt(i,j,k)*ft
          deltaY = Grd%dyt(i,j)*Thickness%dzt(i,j,k)*ft
          if (delta_iso >= deltaX .or. delta_iso >= deltaY) then
            i_delta = i; j_delta = j; k_delta = k
            delta_iso = min(deltaX,deltaY)
          endif
        endif  
      enddo
    enddo
  enddo
  delta_iso  = delta_iso + 1.e-6*mpp_pe() ! to separate redundancies
  delta_iso0 = delta_iso
  call mpp_min (delta_iso)

  ! show most unstable location
  if (delta_iso == delta_iso0) then
     write(unit,'(a)')
     write(unit,'(a)')'---Neutral direction slope check I for linear stability of neutral diffusion---'
     write(unit,'(a,e14.7)')'With a neutral physics time step (secs) of ', dtime
     write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
     write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                   = ',Grd%xt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                   = ',Grd%yt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,') = ', &
                                           Thickness%dzt(i_delta,j_delta,k_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'aredi(',i_delta,',',j_delta,',',k_delta,') = ', &
                                           aredi_array(i_delta,j_delta,k_delta)
     write(unit,'(a,e14.7,a)')'delta_iso           = ',delta_iso,' is the maximum neutral direction slope' 
     write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
     write(unit,'(a)')'The namelist parameter smax should conservatively be <= delta_iso.'
     if(smax >= delta_iso) then 
        write(unit,'(a,f10.5,a)')'==> Warning: The namelist parameter smax= ',smax, ' is >= to delta_iso.'
        write(unit,'(a)')'Linear stability of the neutral diffusion scheme may be compromised.'
     endif
     write(unit,'(a)')
  endif


  ! Compute maximum diffusivity available given a maximum slope of smax
  i_delta = isc; j_delta = jsc; k_delta = 1
  ft = 0.5/(smax*dtime + epsln)
  delta_min = Thickness%dzt(i_delta,j_delta,k_delta)*Grd%dxt(i_delta,j_delta)
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,k) /= 0.0) then        
          delta = min(Grd%dxt(i,j),Grd%dyt(i,j))*Thickness%dzt(i,j,k)
          if (delta_min > delta) then
            i_delta = i; j_delta = j; k_delta = k
            delta_min = delta
          endif
        endif   
      enddo
    enddo
  enddo
  A_max  = ft*delta_min + 1.e-6*mpp_pe() ! to separate redundancies
  A_max0 = A_max
  call mpp_min (A_max)

  ! show most unstable location
  if (A_max == A_max0) then
     write(unit,'(a)')
     write(unit,'(a)')'---Neutral direction slope check II for linear stability of neutral diffusion---'
     write(unit,'(a,e14.7)')'Assuming maximum Redi neutral diffusion slope of ', smax
     write(unit,'(a,e14.7)')'and neutral physics time step (secs) of ', dtime
     write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
     write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                  = ',&
     Grd%xt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                  = ',&
     Grd%yt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,')= ',&
     Thickness%dzt(i_delta,j_delta,k_delta)
     write(unit,'(a,e14.7,a)')'A_max      = ',A_max,' (m^2/sec) is the maximum neutral diffusivity'
     write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
     write(unit,'(a)')'Conservatively, neutral diffusivities used in the model should be less than A_max.'
     write(unit,'(a)')'--------------------------------------------------------------------------------'
     write(unit,'(a)')
  endif


end subroutine ocean_nphysics_coeff_init
! </SUBROUTINE>  NAME="ocean_nphysics_coeff_init"



!#######################################################################
! <SUBROUTINE NAME="tracer_derivs">
!
! <DESCRIPTION>
! Compute the tracer derivatives.
!
! Horizontal derivatives are taken along surfaces of 
! constant vertical coordinate (constant k-level)
!
! This approach ensures that when neutral physics defaults to "horizontal" physics
! next to boundaries, it will do so as horizontal, defined along surfaces of constant 
! s-surfaces, and so will not generate spurious extrema.  
!
! Additionally, when using generalized vertical coordinates, the neutral diffusion
! slope should be computed relative to the s-surfaces.  The skew diffusion slope 
! should ideally be computed with respect to z-surfaces, as z-surfaces define
! available potential energy. However, when s and z surfaces are reasonably close, 
! as they are in the interior for zstar and pstar vertical coordinates, then we 
! choose to to dissipate thickness as defined relative to the zstar or pstar surfaces. 
! This should not be such a big deal, and it is certainly easier computationally than
! worrying about computing two separate sets of slopes.  More on this detail is 
! discussed in "Elements of MOM".
! 
! NOTE: This approach is not appropriate for sigma-models. Indeed, many assumptions
! in the neutral physics modules need to be rethought for terrain following vertical
! coordinates. 
!
! Vertical neutral density derivative for use in fz_terms
! and fz_flux, and for use in fx_flux and fy_flux. 
! Note that the derivative at k=nk vanishes by definition
! since these derivatives are at the bottom of tracer cell. 
! also note the use of -epsln_drhodz ensures the vertical 
! derivative is always < 0.  We also support the same 
! approach used in the mom4p0d code for legacy purposes. 
!
! Comments about smoothing drhodz:
!
! 1/ Tests in coupled 1-degree model showed extreme sensitivity 
! of MOC to smoothing.  GFDL users generally do NOT smooth, hence
! the default drhodz_smooth_vert=drhodz_smooth_horz=.false. 
!
! 2/ Smoothing the vertical derivative of drhodzb and drhodzh helps  
! is greatly needed for producing a regularized (i.e., well behaved)
! neutral slope vector.  
!
! 3/ An attempt was made to smooth dTdz and dSdz rather 
! than drhodz.  The resulting slope was smooth, but not as 
! smooth as when acting on drhodz itself.
!
! </DESCRIPTION>
!
subroutine tracer_derivs(Time, taum1, dtime, drhodT, drhodS, T_prog, Dens, dzwtr, &
                         dTdx, dTdy, dTdz, dSdx, dSdy, dSdz, drhodzb, drhodzh)

  type(ocean_time_type),           intent(in)    :: Time
  integer,                         intent(in)    :: taum1
  real,                            intent(in)    :: dtime
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodS
  real, dimension(isd:,jsd:,0:),   intent(in)    :: dzwtr
  type(ocean_prog_tracer_type),    intent(in)    :: T_prog(:)
  type(ocean_density_type),        intent(in)    :: Dens

  type(tracer_3d_1_nk_type),       intent(inout) :: dTdx(:)
  type(tracer_3d_1_nk_type),       intent(inout) :: dTdy(:)
  type(tracer_3d_0_nk_type),       intent(inout) :: dTdz(:)
  type(tracer_3d_1_nk_type),       intent(inout) :: dSdx
  type(tracer_3d_1_nk_type),       intent(inout) :: dSdy
  type(tracer_3d_0_nk_type),       intent(inout) :: dSdz
  real, dimension(isd:,jsd:,:,0:), intent(inout) :: drhodzb
  real, dimension(isd:,jsd:,:,0:), intent(inout) :: drhodzh


  integer :: i, j, k, m, n
  integer :: kp1, kbot
  integer :: kr, kpkr, km1pkr
  real    :: dTdz_ijk1, dTdz_ijk2
  real    :: drhodzb0_prev, drhodzb1_prev
  real    :: drhodzh0_prev, drhodzh1_prev
  real    :: tmpb0, tmpb1, tmph0, tmph1

  call mpp_clock_begin(id_clock_tracer_derivs)


  ! Horizontal derivatives taken along surfaces of 
  ! constant vertical coordinate (constant k-level)
  do n=1,num_prog_tracers
     do k=1,nk
        kp1 = min(k+1,nk)
        wrk1(:,:,1)          = T_prog(n)%field(:,:,k,taum1)
        dTdx(n)%field(:,:,k) = FDX_T(wrk1(:,:,1))*FMX(Grd%tmask(:,:,k))
        dTdy(n)%field(:,:,k) = FDY_T(wrk1(:,:,1))*FMY(Grd%tmask(:,:,k))
     enddo
  enddo
  do k=1,nk
     kp1 = min(k+1,nk)
     wrk1(:,:,1)       = Dens%rho_salinity(:,:,k,taum1)
     dSdx%field(:,:,k) = FDX_T(wrk1(:,:,1))*FMX(Grd%tmask(:,:,k))
     dSdy%field(:,:,k) = FDY_T(wrk1(:,:,1))*FMY(Grd%tmask(:,:,k))
  enddo

  ! vertical derivative 
  do n=1,num_prog_tracers
     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsd,jed
           do i=isd,ied
              dTdz(n)%field(i,j,k) = Grd%tmask(i,j,kp1)*dzwtr(i,j,k) &
                   *(T_prog(n)%field(i,j,k,taum1)-T_prog(n)%field(i,j,kp1,taum1))
           enddo
        enddo
     enddo
  enddo
  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsd,jed
        do i=isd,ied
           dSdz%field(i,j,k) = Grd%tmask(i,j,kp1)*dzwtr(i,j,k) &
                *(Dens%rho_salinity(i,j,k,taum1)-Dens%rho_salinity(i,j,kp1,taum1))
        enddo
     enddo
  enddo


  ! vertical neutral density derivative
  if(drhodz_mom4p1) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied

               dTdz_ijk1 = dTdz(index_temp)%field(i,j,k) ; dTdz_ijk2 = dSdz%field(i,j,k)

               kr=0
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 + drhodS(i,j,kpkr)*dTdz_ijk2 
               drhodzb(i,j,k,kr) = min(drhodzb(i,j,k,kr), -epsln_drhodz)

               kr=1
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 + drhodS(i,j,kpkr)*dTdz_ijk2
               drhodzb(i,j,k,kr) = min(drhodzb(i,j,k,kr), -epsln_drhodz)


               kr=0
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dSdz%field(i,j,km1pkr)
               drhodzh(i,j,k,kr) = min(drhodzh(i,j,k,kr), -epsln_drhodz)

               kr=1
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dSdz%field(i,j,km1pkr)
               drhodzh(i,j,k,kr) = min(drhodzh(i,j,k,kr), -epsln_drhodz)

            enddo
         enddo
      enddo

  else 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied

               dTdz_ijk1 = dTdz(index_temp)%field(i,j,k) ; dTdz_ijk2 = dSdz%field(i,j,k)

               kr=0
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 &
                                 + drhodS(i,j,kpkr)*dTdz_ijk2 -epsln

               kr=1
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 &
                                 + drhodS(i,j,kpkr)*dTdz_ijk2 -epsln

               kr=0
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dSdz%field(i,j,km1pkr) -epsln

               kr=1
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dSdz%field(i,j,km1pkr) -epsln

            enddo
         enddo
      enddo
  endif

  ! vertically smooth the vertical derivative of density to 
  ! produce a smooth neutral slope vector for all flux components.
  if(drhodz_smooth_vert) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied

               drhodzb0_prev = onefourth*drhodzb(i,j,1,0)
               drhodzb1_prev = onefourth*drhodzb(i,j,1,1)
               drhodzh0_prev = onefourth*drhodzh(i,j,1,0)
               drhodzh1_prev = onefourth*drhodzh(i,j,1,1)

               kbot=Grd%kmt(i,j)
               if(kbot>1) then
                   do k=2,kbot-2

                      tmpb0            = drhodzb(i,j,k,0)
                      drhodzb(i,j,k,0) = drhodzb0_prev + onehalf*drhodzb(i,j,k,0) &
                                                       + onefourth*drhodzb(i,j,k+1,0)
                      drhodzb0_prev    = onefourth*tmpb0

                      tmpb1            = drhodzb(i,j,k,1)
                      drhodzb(i,j,k,1) = drhodzb1_prev + onehalf*drhodzb(i,j,k,1) &
                                                       + onefourth*drhodzb(i,j,k+1,1)
                      drhodzb1_prev    = onefourth*tmpb1

                      tmph0            = drhodzh(i,j,k,0)
                      drhodzh(i,j,k,0) = drhodzh0_prev + onehalf*drhodzh(i,j,k,0) &
                                                       + onefourth*drhodzh(i,j,k+1,0)
                      drhodzh0_prev    = onefourth*tmph0

                      tmph1            = drhodzh(i,j,k,1)
                      drhodzh(i,j,k,1) = drhodzh1_prev + onehalf*drhodzh(i,j,k,1) &
                                                       + onefourth*drhodzh(i,j,k+1,1)
                      drhodzh1_prev     = onefourth*tmph1

                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! horizontally smooth the vertical derivative of density to 
  ! produce a smooth neutral slope vector for all flux components.
  ! since drhodzb is needed only on the computational domain, there 
  ! is no need to update its values in the halos.  
  if(drhodz_smooth_horz) then 
      do k=1,nk-1
         drhodzb(:,:,k,0) = drhodzb(:,:,k,0) + dtime*LAP_T(drhodzb(:,:,k,0),smooth_lap(:,:))
         drhodzb(:,:,k,1) = drhodzb(:,:,k,1) + dtime*LAP_T(drhodzb(:,:,k,1),smooth_lap(:,:))
         drhodzh(:,:,k,0) = drhodzh(:,:,k,0) + dtime*LAP_T(drhodzh(:,:,k,0),smooth_lap(:,:))
         drhodzh(:,:,k,1) = drhodzh(:,:,k,1) + dtime*LAP_T(drhodzh(:,:,k,1),smooth_lap(:,:))
      enddo
      call mpp_update_domains(drhodzh(:,:,:,0), Dom%domain2d, complete=.false.) 
      call mpp_update_domains(drhodzh(:,:,:,1), Dom%domain2d, complete=.true.) 
  endif

  ! compute squared buoyancy frequency 
  wrk2(:,:,:)=0.0
  do k=1,nk
     wrk2(:,:,k) = -onehalf*grav*rho0r*(drhodzb(:,:,k,0)+drhodzb(:,:,k,1))
  enddo
  call diagnose_3d(Time, Grd, id_N2slope, wrk2(:,:,:))

  call mpp_clock_end(id_clock_tracer_derivs)

end subroutine tracer_derivs
! </SUBROUTINE> NAME="tracer_derivs"



!#######################################################################
! <SUBROUTINE NAME="neutral_slopes">
!
! <DESCRIPTION>
! Subroutine computes the neutral slopes for the triads associated 
! with the vertical flux component.  
!
! Array tensor_31 initially holds the x-slope used for flux component fz.
! Array tensor_32 initially holds the y-slope used for flux component fz.
!
! In subsequent calculations, these arrays will be multipied by the
! diffusivities.  
!
! No slope tapering is applied in this routine. 
!
! slopes are computed over k=1,nk-1, since the slope at k=nk 
! should be zero. 
!
! </DESCRIPTION>
!
subroutine neutral_slopes(Time, dTdx, dTdy, dSdx, dSdy, drhodT, drhodS, drhodzb, tensor_31, tensor_32)

  type(ocean_time_type),              intent(in)    :: Time
  type(tracer_3d_1_nk_type),          intent(in)    :: dTdx(:)
  type(tracer_3d_1_nk_type),          intent(in)    :: dTdy(:)
  type(tracer_3d_1_nk_type),          intent(in)    :: dSdx
  type(tracer_3d_1_nk_type),          intent(in)    :: dSdy
  real, dimension(isd:,jsd:,:),       intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),       intent(in)    :: drhodS
  real, dimension(isd:,jsd:,:,0:),    intent(in)    :: drhodzb
  real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: tensor_31
  real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: tensor_32

  integer :: i, j, k, kp1
  real :: drhodT_ijk, drhodS_ijk, drhodT_ijkp1, drhodS_ijkp1
  real :: tmask_ijkp1

  real :: dTdx_ijk1, dTdx_ijk2, dTdx_im1jk1, dTdx_im1jk2
  real :: dTdx_ijkp11, dTdx_ijkp12, dTdx_im1jkp11, dTdx_im1jkp12
  real :: dTdy_ijk1, dTdy_ijk2, dTdy_ijm1k1, dTdy_ijm1k2
  real :: dTdy_ijkp11, dTdy_ijkp12, dTdy_ijm1kp11, dTdy_ijm1kp12
  real :: drhodzbr_ijk0, drhodzbr_ijk1

  call mpp_clock_begin(id_clock_neutral_slopes)

  tensor_31(:,:,:,:,:) = 0.0
  tensor_32(:,:,:,:,:) = 0.0

  do k=1,nk-1
     kp1 = k+1
     do j=jsc,jec
        do i=isc,iec

           tmask_ijkp1    = Grd%tmask(i,j,kp1)

           drhodT_ijk     = drhodT(i,j,k)
           drhodS_ijk     = drhodS(i,j,k)
           drhodT_ijkp1   = drhodT(i,j,kp1)
           drhodS_ijkp1   = drhodS(i,j,kp1)

           dTdx_ijk1     = dTdx(index_temp)%field(i,j,k)
           dTdx_ijk2     = dSdx%field(i,j,k)
           dTdx_im1jk1   = dTdx(index_temp)%field(i-1,j,k)
           dTdx_im1jk2   = dSdx%field(i-1,j,k)
           dTdx_ijkp11   = dTdx(index_temp)%field(i,j,kp1)
           dTdx_ijkp12   = dSdx%field(i,j,kp1)
           dTdx_im1jkp11 = dTdx(index_temp)%field(i-1,j,kp1)
           dTdx_im1jkp12 = dSdx%field(i-1,j,kp1)

           dTdy_ijk1     = dTdy(index_temp)%field(i,j,k)
           dTdy_ijk2     = dSdy%field(i,j,k)
           dTdy_ijm1k1   = dTdy(index_temp)%field(i,j-1,k)
           dTdy_ijm1k2   = dSdy%field(i,j-1,k)
           dTdy_ijkp11   = dTdy(index_temp)%field(i,j,kp1)
           dTdy_ijkp12   = dSdy%field(i,j,kp1)
           dTdy_ijm1kp11 = dTdy(index_temp)%field(i,j-1,kp1)
           dTdy_ijm1kp12 = dSdy%field(i,j-1,kp1)

           drhodzbr_ijk0 = 1.0/drhodzb(i,j,k,0)
           drhodzbr_ijk1 = 1.0/drhodzb(i,j,k,1)


           ! ip=jq=0

           !   kr=0
           tensor_31(i,j,k,0,0) = -tmask_ijkp1                     &
                *(drhodT_ijk*dTdx_im1jk1 + drhodS_ijk*dTdx_im1jk2) &
                *drhodzbr_ijk0
           tensor_32(i,j,k,0,0) = -tmask_ijkp1                     &
                *(drhodT_ijk*dTdy_ijm1k1 + drhodS_ijk*dTdy_ijm1k2) &
                *drhodzbr_ijk0

           !   kr=1 
           tensor_31(i,j,k,0,1) = -tmask_ijkp1                             &
                *(drhodT_ijkp1*dTdx_im1jkp11 + drhodS_ijkp1*dTdx_im1jkp12) &
                *drhodzbr_ijk1

           tensor_32(i,j,k,0,1) = -tmask_ijkp1                             &
                *(drhodT_ijkp1*dTdy_ijm1kp11 + drhodS_ijkp1*dTdy_ijm1kp12) &
                *drhodzbr_ijk1


           ! ip=jq=1

           !   kr=0 
           tensor_31(i,j,k,1,0) = -tmask_ijkp1                 &
                *(drhodT_ijk*dTdx_ijk1 + drhodS_ijk*dTdx_ijk2) &
                *drhodzbr_ijk0

           tensor_32(i,j,k,1,0) = -tmask_ijkp1                 &
                *(drhodT_ijk*dTdy_ijk1 + drhodS_ijk*dTdy_ijk2) &
                *drhodzbr_ijk0
           !   kr=1 
           tensor_31(i,j,k,1,1) = -tmask_ijkp1                         &
                *(drhodT_ijkp1*dTdx_ijkp11 + drhodS_ijkp1*dTdx_ijkp12) &
                *drhodzbr_ijk1

           tensor_32(i,j,k,1,1) = -tmask_ijkp1                         &
                *(drhodT_ijkp1*dTdy_ijkp11 + drhodS_ijkp1*dTdy_ijkp12) &
                *drhodzbr_ijk1

        enddo
     enddo
  enddo


  ! send to diagnostic manager 
  if (id_slope31 > 0) then 
       wrk1(isc:iec,jsc:jec,:) = onefourth*                                     &
         (tensor_31(isc:iec,jsc:jec,:,0,0) + tensor_31(isc:iec,jsc:jec,:,0,1) + &
          tensor_31(isc:iec,jsc:jec,:,1,0) + tensor_31(isc:iec,jsc:jec,:,1,1))
       call diagnose_3d(Time, Grd, id_slope31, wrk1(:,:,:))
  endif 
  if (id_slope32 > 0) then 
       wrk1(isc:iec,jsc:jec,:) = onefourth*                                     &
         (tensor_32(isc:iec,jsc:jec,:,0,0) + tensor_32(isc:iec,jsc:jec,:,0,1) + &
          tensor_32(isc:iec,jsc:jec,:,1,0) + tensor_32(isc:iec,jsc:jec,:,1,1))
       call diagnose_3d(Time, Grd, id_slope32, wrk1(:,:,:))
  endif  

  call mpp_clock_end(id_clock_neutral_slopes)

end subroutine neutral_slopes
! </SUBROUTINE> NAME="neutral_slopes"



!#######################################################################
! <SUBROUTINE NAME="compute_eady_rate">
!
! <DESCRIPTION>
! Finish computing eady growth rate.
! </DESCRIPTION>
!
subroutine compute_eady_rate(Time, Thickness, T_prog, Dens, eady_termx, eady_termy)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  real,  dimension(isd:,jsd:,:),intent(in) :: eady_termx
  real,  dimension(isd:,jsd:,:),intent(in) :: eady_termy

  integer :: i, j, k, m, kbot, tau
  real    :: eady_rate_prev, tmp
  real    :: eady_mld(isd:ied,jsd:jed)

  if(.not. agm_closure) return 

  call mpp_clock_begin(id_clock_compute_eady_rate)

  tau           = Time%tau 
  wrk1(:,:,:)   = Grd%tmask(:,:,:)   
  eady_mld(:,:) = 0.0


  ! raw Eady growth rate, with regularization not yet applied 
  do k=1,nk
     eady_rate(COMP,k) = Grd%tmask(COMP,k) &
          *sqrt((eady_termx(COMP,k)*count_x(COMP)+eady_termy(COMP,k)*count_y(COMP))*gravrho0r)
  enddo

  ! apply cap to the Eady growth rate
  if(agm_closure_eady_cap) then 
      do k=1,nk
         eady_rate(COMP,k) = 2.0*eady_rate(COMP,k)*agm_growth_rate_max(COMP) &
              /(epsln+eady_rate(COMP,k)+agm_growth_rate_max(COMP))
      enddo
  endif


  ! vertically average Eady within surface mixed layer 
  if(agm_closure_eady_ave_mixed) then

      call calc_mixed_layer_depth(Thickness,                &
           T_prog(index_salt)%field(isd:ied,jsd:jed,:,tau), &
           Dens%rho_salinity(isd:ied,jsd:jed,:,tau),        &
           Dens%rho(isd:ied,jsd:jed,:,tau),                 &
           Dens%pressure_at_depth(isd:ied,jsd:jed,:),       &
           eady_mld)

      do j=jsc,jec
         do i=isc,iec
            eady_mld(i,j) = Grd%tmask(i,j,1)*min(eady_mld(i,j), Grd%ht(i,j))
         enddo
      enddo

      ! do not believe eady_rate at k=1, so always skip its contribution 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      do k=2,nk
         do j=jsc,jec
            do i=isc,iec
               if(Thickness%depth_zt(i,j,k) < eady_mld(i,j)) then 
                   wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%dzt(i,j,k)
                   wrk2_2d(i,j) = wrk2_2d(i,j) + eady_rate(i,j,k)*Thickness%dzt(i,j,k)
               endif
            enddo
         enddo
      enddo
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if(Thickness%depth_zt(i,j,k) <= eady_mld(i,j)) then 
                   eady_rate(i,j,k) = Grd%tmask(i,j,k)*wrk2_2d(i,j)/(wrk1_2d(i,j)+epsln)
               endif
            enddo
         enddo
      enddo

  endif ! endif for agm_closure_eady_ave_mixed


  ! apply vertical 1-2-1 smoothing 
  if(agm_closure_eady_smooth_vert) then 
      do m=1,num_121_passes
         do j=jsc,jec
            do i=isc,iec
               eady_rate_prev = onefourth*eady_rate(i,j,1)
               kbot=Grd%kmt(i,j)
               if(kbot>1) then
                   do k=2,kbot-2
                      tmp              = eady_rate(i,j,k)
                      eady_rate(i,j,k) = eady_rate_prev + onehalf*eady_rate(i,j,k) &
                                                        + onefourth*eady_rate(i,j,k+1)
                      eady_rate_prev   = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif 

  ! apply horizontal 1-2-1 smoothing 
  if(agm_closure_eady_smooth_horz) then 
      call mpp_update_domains (eady_rate(:,:,:), Dom%domain2d) 
      do k=1,nk
         eady_rate(:,:,k) = S2D(eady_rate(:,:,k))
      enddo
  endif 

  ! compute vertical average over specified depth 
  eady_rate_zave(:,:) = 0.0
  wrk1_2d(:,:)        = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(Thickness%depth_zwt(i,j,k) >= agm_closure_upper_depth .and. &
              Thickness%depth_zwt(i,j,k) <= agm_closure_lower_depth) then 
              eady_rate_zave(i,j) = eady_rate_zave(i,j) + eady_rate(i,j,k)*Thickness%dzwt(i,j,k)
              wrk1_2d(i,j)        = wrk1_2d(i,j) + Thickness%dzwt(i,j,k)
           endif
        enddo
     enddo
  enddo
  eady_rate_zave(COMP) = Grd%tmask(COMP,1)*eady_rate_zave(COMP)/(wrk1_2d(COMP)+epsln)

  call diagnose_3d(Time, Grd, id_eady_rate, eady_rate(:,:,:))
  call diagnose_2d(Time, Grd, id_eady_rate_zave, eady_rate_zave(:,:))
  call diagnose_2d(Time, Grd, id_eady_mld, eady_mld(:,:))

  call mpp_clock_end(id_clock_compute_eady_rate)

end subroutine compute_eady_rate
! </SUBROUTINE> NAME="compute_eady_rate"



!#######################################################################
! <SUBROUTINE NAME="compute_baroclinicity">
!
! <DESCRIPTION>
! Finish computing baroclinicity, which is defined to be the vertically
! averaged magnitude of the horizontal density gradient.
! </DESCRIPTION>
!
subroutine compute_baroclinicity(Time, baroclinic_termx, baroclinic_termy)
  type(ocean_time_type),        intent(in) :: Time
  real,  dimension(isd:,jsd:),  intent(in) :: baroclinic_termx
  real,  dimension(isd:,jsd:),  intent(in) :: baroclinic_termy

  real :: vertical_range

  if(.not. agm_closure) return 

  call mpp_clock_begin(id_clock_compute_baroclinicity)
  
  vertical_range = epsln + agm_closure_lower_depth - agm_closure_upper_depth

  baroclinicity(COMP) = Grd%tmask(COMP,1)*                                       &
       (baroclinic_termx(COMP)*count_x(COMP)+baroclinic_termy(COMP)*count_y(COMP)) &
       /vertical_range

  call diagnose_2d(Time, Grd, id_baroclinicity, baroclinicity(:,:))

  call mpp_clock_end(id_clock_compute_baroclinicity)

end subroutine compute_baroclinicity
! </SUBROUTINE> NAME="compute_baroclinicity"


!#######################################################################
! <SUBROUTINE NAME="compute_rossby_radius">
!
! <DESCRIPTION>
! Subroutine computes the first baroclinic Rossby radius of deformation. 
! Employ WKB approach described by Chelton et al.  In particular, 
! use formulae (2.2), (2.3a) and (2.3b) from their paper. 
!
! Place a max and min value on the Rossby radius.
!
! Compute buoyancy frequency in terms of vertical gradient of 
! locally referenced potential density.  Place the reference point
! at the interface between the tracer cells, which is also where 
! the vertical derivative of neutral density is located.  This amounts 
! to a centered difference computation similar to that used by 
! Chelton et al. equation (B.4). 
! </DESCRIPTION>
!
subroutine compute_rossby_radius(Thickness, dTdz, dSdz, Time, drhodT, drhodS, &
                                 rossby_radius, rossby_radius_raw)

  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(tracer_3d_0_nk_type),     intent(in)    :: dTdz(:)
  type(tracer_3d_0_nk_type),     intent(in)    :: dSdz
  type(ocean_time_type),         intent(in)    :: Time
  real, dimension(isd:,jsd:,:),  intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),  intent(in)    :: drhodS
  real, dimension(isd:,jsd:),    intent(inout) :: rossby_radius
  real, dimension(isd:,jsd:),    intent(inout) :: rossby_radius_raw

  real     :: drhodzb_speed
  integer  :: i,j,k

  call mpp_clock_begin(id_clock_compute_rossby_radius)
 
  gravity_wave_speed(:,:) = 0.0
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           drhodzb_speed  = 0.5*( (drhodT(i,j,k) + drhodT(i,j,k+1))*dTdz(index_temp)%field(i,j,k) &
                                 +(drhodS(i,j,k) + drhodS(i,j,k+1))*dSdz%field(i,j,k) )
           gravity_wave_speed(i,j) = gravity_wave_speed(i,j) &
                                   + Thickness%dzwt(i,j,k)*sqrt(abs(gravrho0r*drhodzb_speed))
        enddo
     enddo
  enddo
  gravity_wave_speed(:,:) = gravity_wave_speed(:,:)/pi

  ! for Eden near the equator 
  sqrt2betaCr(:,:) = sqrt(2.0*gravity_wave_speed(:,:)*beta_param(:,:))

  wrk1_2d(:,:)           = 0.0
  wrk2_2d(:,:)           = 0.0
  rossby_radius_raw(:,:) = 0.0
  wrk1_2d(COMP)           = gravity_wave_speed(COMP)/(coriolis_param(COMP) + epsln)
  wrk2_2d(COMP)           = sqrt(gravity_wave_speed(COMP)/(2.0*beta_param(COMP)+epsln))
  rossby_radius_raw(COMP) = min(wrk1_2d(COMP), wrk2_2d(COMP))

  agm_grid_scaling(:,:) = 1.0
  if(agm_closure_grid_scaling) then 

      ! for backward bitwise compatibility 
      if(agm_closure_grid_scaling_power==2.0) then 
         agm_grid_scaling(COMP) =  Grd%tmask(COMP,1) &
              *grid_length(COMP)**2/(grid_length(COMP)**2 + rossby_radius_raw(COMP)**2)
      else  
         agm_grid_scaling(COMP) =  Grd%tmask(COMP,1)                 &
              *grid_length(COMP)**agm_closure_grid_scaling_power    &
              /( grid_length(COMP)**agm_closure_grid_scaling_power  &
              + rossby_radius_raw(COMP)**agm_closure_grid_scaling_power)
      endif

  endif

  aredi_grid_scaling(:,:) = 1.0
  if(aredi_diffusivity_grid_scaling) then 
     aredi_grid_scaling(COMP) = Grd%tmask(COMP,1)              &
          *grid_length(COMP)**agm_closure_grid_scaling_power &
          /( grid_length(COMP)**agm_closure_grid_scaling_power &
          + rossby_radius_raw(COMP)**agm_closure_grid_scaling_power)
      call mpp_update_domains (aredi_grid_scaling(:,:), Dom%domain2d)
  endif

  rossby_radius(COMP) = min(rossby_radius_max,max(rossby_radius_min,rossby_radius_raw(COMP)))
   
  call diagnose_2d(Time, Grd, id_rossby_radius,rossby_radius_raw(:,:))
  call diagnose_2d(Time, Grd, id_rossby, rossby_radius(:,:))
  call diagnose_2d(Time, Grd, id_rossby_nonequator, wrk1_2d(:,:))
  call diagnose_2d(Time, Grd, id_rossby_equator, wrk2_2d(:,:))
  call diagnose_2d(Time, Grd, id_gw_speed, gravity_wave_speed(:,:))
  call diagnose_2d(Time, Grd, id_sqrt2betaCr, sqrt2betaCr(:,:))
  call diagnose_2d(Time, Grd, id_agm_grid_scaling, agm_grid_scaling(:,:))
  call diagnose_2d(Time, Grd, id_aredi_grid_scaling, aredi_grid_scaling(:,:))

  if(debug_this_module) then 
      call write_chksum_header(FILENAME, ' ', Time%model_time)
      call write_chksum_2d('checksum rossby_radius', rossby_radius(COMP)*Grd%tmask(COMP,1))
      call write_chksum_2d('checksum rossby_radius_raw', rossby_radius_raw(COMP)*Grd%tmask(COMP,1))
  endif

  call mpp_clock_end(id_clock_compute_rossby_radius)

end subroutine compute_rossby_radius
! </SUBROUTINE> NAME="compute_rossby_radius"



!#######################################################################
! <SUBROUTINE NAME="compute_bczone_radius">
!
! <DESCRIPTION>
! Subroutine computes the radius of the baroclinic zone in a manner 
! suggested by the Hadley Centre approach (Malcolm Roberts, personal 
! communication).  
!
! Algorithm is used in MOM3 and documented in the MOM3 Manual.
!
! </DESCRIPTION>
!
subroutine compute_bczone_radius(Time, bczone_radius)

  type(ocean_time_type),      intent(in)    :: Time
  real, dimension(isd:,jsd:), intent(inout) :: bczone_radius

  integer :: i,j
  real    :: n_zone, e_zone, s_zone, w_zone, fract, nstot, ewtot
  integer :: ip, jq

  if(.not. agm_closure) return
  if(.not. agm_closure_length_bczone) return 

  call mpp_clock_begin(id_clock_compute_bczone_radius)

  bczone_rate(COMP) = eady_rate_zave(COMP)
  call mpp_update_domains (bczone_rate(:,:), BCzone_domain%domain2d)

  do j=jsc,jec
    do i=isc,iec
      if (bczone_rate(i,j) > agm_closure_bczone_crit_rate) then

        ! search northward 
        n_zone = bczone_dyt(i,j)
        do jq=j+1,j+bczone_max_pts
          if (bczone_rate(i,jq) > agm_closure_bczone_crit_rate) then
            n_zone = n_zone + bczone_dyt(i,jq)
          else
            exit
          endif
        enddo

        ! search southward 
        s_zone = bczone_dyt(i,j)
        do jq=j-1,j-bczone_max_pts,-1
          if (bczone_rate(i,jq) > agm_closure_bczone_crit_rate) then
            s_zone = s_zone + bczone_dyt(i,jq)
          else
            exit
          endif
        enddo

        ! search eastward 
        e_zone = bczone_dxt(i,j)
        do ip=i+1,i+bczone_max_pts
          if (bczone_rate(ip,j) > agm_closure_bczone_crit_rate) then
            e_zone = e_zone + bczone_dxt(ip,j)
          else
            exit
          endif
        enddo

        ! search westward 
        w_zone = bczone_dxt(i,j)
        do ip=i-1,i-bczone_max_pts,-1
          if (bczone_rate(ip,j) > agm_closure_bczone_crit_rate) then
            w_zone = w_zone + bczone_dxt(ip,j)
          else
            exit
          endif
        enddo

        ! total radius (subtraction accounts for double-counting central point)
        nstot=n_zone+s_zone-bczone_dyt(i,j)
        ewtot=e_zone+w_zone-bczone_dxt(i,j)

        if (nstot < ewtot) then 
          fract=min(n_zone,s_zone)/max(n_zone,s_zone) 
          bczone_radius(i,j)=fract*nstot 
        else   
          fract=min(e_zone,w_zone)/max(e_zone,w_zone)
          bczone_radius(i,j)=fract*ewtot
        endif
      endif
          
    enddo
  enddo

  call diagnose_2d(Time, Grd, id_bczone, bczone_radius(:,:))

  if(debug_this_module) then 
      call write_chksum_header(FILENAME, ' ', Time%model_time)
      call write_chksum_2d('checksum bczone_radius', bczone_radius(COMP)*Grd%tmask(COMP,1))
  endif


  call mpp_clock_end(id_clock_compute_bczone_radius)

end subroutine compute_bczone_radius
! </SUBROUTINE> NAME="compute_bczone_radius"



!#######################################################################
! <SUBROUTINE NAME="compute_diffusivity">
!
! <DESCRIPTION>
! Subroutine computes flow dependent diffusivity.
! Allow for an added dimensionless tuning factor as well as a 
! minimum and maximum diffusivity. 
! </DESCRIPTION>
!
subroutine compute_diffusivity(Time, ksurf_blayer, drhodz_zt, rossby_radius, &
                               bczone_radius, agm_array, aredi_array)

  type(ocean_time_type),        intent(in)    :: Time
  integer, dimension(isd:,jsd:),intent(in)    :: ksurf_blayer
  real, dimension(isd:,jsd:,:), intent(in)    :: drhodz_zt
  real, dimension(isd:,jsd:),   intent(in)    :: rossby_radius
  real, dimension(isd:,jsd:),   intent(in)    :: bczone_radius
  real, dimension(isd:,jsd:,:), intent(inout) :: agm_array
  real, dimension(isd:,jsd:,:), intent(inout) :: aredi_array

  integer :: i, j, k
  real    :: denom, param, active_cells, tmp

  call mpp_clock_begin(id_clock_compute_diffusivity)

  if(agm_closure) then 


     agm_growth_rate = 0.0
     agm_length      = 0.0
     wrk1            = 0.0 !growth_rate
     wrk2            = 0.0 !rhines length 
     wrk3            = 0.0 !agm_qg
     wrk4            = 0.0 !agm_fast
     wrk2_2d(:,:)    = 0.0 !drhodz_ref

     ! set the length scale 
     if(agm_closure_length_fixed .or. agm_closure_baroclinic) then 
         agm_length(:,:,:) = Grd%tmask(:,:,:)*agm_closure_length

     elseif(agm_closure_length_rossby) then 
         do k=1,nk
            agm_length(COMP,k) = Grd%tmask(COMP,k)* &
                 2.0*grid_length(COMP)*rossby_radius(COMP)/(grid_length(COMP)+rossby_radius(COMP)+epsln)
         enddo

     elseif(agm_closure_length_bczone) then 
         do k=1,nk
            agm_length(COMP,k) = Grd%tmask(COMP,k)* &
                 2.0*grid_length(COMP)*bczone_radius(COMP)/(grid_length(COMP)+bczone_radius(COMP)+epsln)
         enddo

     elseif(agm_closure_eden_greatbatch) then
         if(agm_closure_eden_length_const) then 
             do k=1,nk
                agm_length(COMP,k) = agm_closure_eden_length*Grd%tmask(COMP,k)
             enddo
         else 
             do k=1,nk
                wrk2(COMP,k)       = eady_rate(COMP,k)/(beta_param(COMP)+epsln)
                agm_length(COMP,k) = min(wrk2(COMP,k),rossby_radius(COMP))
                agm_length(COMP,k) = Grd%tmask(COMP,k)* &
                     2.0*grid_length(COMP)*agm_length(COMP,k)/(grid_length(COMP)+agm_length(COMP,k)+epsln)
             enddo
         endif
     endif

     if(agm_closure_length_cap) then 
         do k=1,nk
            agm_length(COMP,k) = Grd%tmask(COMP,k)*min(agm_closure_length_max, agm_length(COMP,k))
         enddo
     endif


     ! set the growth rates 
     if(agm_closure_baroclinic) then 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp = gravrho0r_buoyr*baroclinicity(i,j)
                  wrk1(i,j,k) = 2.0*tmp*agm_growth_rate_max(i,j) &
                                /(tmp+agm_growth_rate_max(i,j)+epsln)  
               enddo
            enddo
         enddo

     elseif(agm_closure_eden_greatbatch) then 
         if(agm_closure_eden_gamma /= 0.0) then 
            do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      param = max(coriolis_param(i,j), sqrt2betaCr(i,j))
                      denom = sqrt(param**2 + agm_closure_eden_gamma*eady_rate(i,j,k)**2)+epsln
                      wrk1(i,j,k) = param*eady_rate(i,j,k)/denom
                   enddo
               enddo
            enddo
         else 
            do k=1,nk
               wrk1(COMP,k) = eady_rate(COMP,k)
            enddo
         endif 
     else
         do k=1,nk
            wrk1(COMP,k) = eady_rate_zave(COMP)
         enddo
     endif


     ! diffusivity computed as coeff*(N/Nref)**2
     if(agm_closure_n2_scale) then 

         if(agm_closure_n2_scale_nref_cst) then 
             wrk2_2d(:,:) = rho0*grav_r*agm_closure_buoy_freq**2
         else 
             ! reference drhodz taken one level beneath blayer base,
             ! and no deeper than one cell from bottom.
             do j=jsc,jec
                do i=isc,iec
                   if(ksurf_blayer(i,j) > 0 .and. Grd%kmt(i,j) > 1) then 
                       k = min(Grd%kmt(i,j)-1, ksurf_blayer(i,j)+1)
                       wrk2_2d(i,j) = abs(drhodz_zt(i,j,k))
                   endif
                enddo
             enddo
         endif

         do k=1,nk
            agm_growth_rate(COMP,k) = Grd%tmask(COMP,k)*wrk1(COMP,k)   ! computed just for diagnostics 
            wrk3(COMP,k) = agm_closure_n2_scale_coeff*Grd%tmask(COMP,k) &
                 *abs(drhodz_zt(COMP,k))/(epsln+wrk2_2d(COMP))
            wrk4(COMP,k) = agm_grid_scaling(COMP)*(wrk3(COMP,k) + agm_micom(COMP)) 
            wrk4(COMP,k) = Grd%tmask(COMP,k)*max(agm_closure_min, min(agm_closure_max, wrk4(COMP,k)))      
         enddo


     ! diffusivity computed as growth_rate*length_scale**2
     else 

         do k=1,nk
            agm_growth_rate(COMP,k) = Grd%tmask(COMP,k)*wrk1(COMP,k)
            wrk3(COMP,k) = agm_closure_scaling*agm_growth_rate(COMP,k)*agm_length(COMP,k)**2 
            wrk4(COMP,k) = agm_grid_scaling(COMP)*(wrk3(COMP,k) + agm_micom(COMP)) 
            wrk4(COMP,k) = Grd%tmask(COMP,k)*max(agm_closure_min, min(agm_closure_max, wrk4(COMP,k)))      
         enddo

     endif


     ! time damping to get slowly evolving diffusivity 
     if(agm_smooth_time) then 
         do k=1,nk
            agm_array(COMP,k) = &
                 Grd%tmask(COMP,k)*(agm_array(COMP,k) - gamma_damp*(agm_array(COMP,k)-wrk4(COMP,k)))
         enddo
     else
         do k=1,nk
            agm_array(COMP,k) = wrk4(COMP,k)
         enddo
     endif
  

     ! spatial smoothing 
     if(agm_smooth_space) then 
         wrk1_2d(:,:) = 0.0
         call mpp_update_domains (agm_array(:,:,:), Dom%domain2d) 
         do k=1,nk

            do j=jsc,jec
               do i=isc,iec
                  if(Grd%tmask(i,j,k)==1.0) then 
                      active_cells = 4.0    +&
                        Grd%tmask(i-1,j,k)  +&
                        Grd%tmask(i+1,j,k)  +&
                        Grd%tmask(i,j-1,k)  +&
                        Grd%tmask(i,j+1,k)
                      if (active_cells > 4.0) then
                          wrk1_2d(i,j) = & 
                            (4.0*agm_array(i,j,k) +&
                            agm_array(i-1,j,k)    +&
                            agm_array(i+1,j,k)    +&
                            agm_array(i,j-1,k)    +&
                            agm_array(i,j+1,k)) / active_cells
                      else
                          wrk1_2d(i,j) = agm_array(i,j,k)
                      endif
                  endif
               enddo
            enddo
            agm_array(COMP,k) = wrk1_2d(COMP)*Grd%tmask(COMP,k)
         enddo
     endif

     ! need agm_array on full data domain 
     call mpp_update_domains (agm_array(:,:,:), Dom%domain2d) 

     call diagnose_3d(Time, Grd, id_agm_growth_rate, agm_growth_rate(:,:,:))
     call diagnose_3d(Time, Grd, id_agm_length, agm_length(:,:,:))
     call diagnose_3d(Time, Grd, id_rhines_length, wrk2(:,:,:))
     call diagnose_3d(Time, Grd, id_agm_qg, wrk3(:,:,:))
     call diagnose_3d(Time, Grd, id_growth_rate_baroclinic, wrk1(:,:,:))

     if(id_N2_nblayer_base > 0) then 
        wrk3_2d(:,:) = grav*rho0r*wrk2_2d(:,:)
        call diagnose_2d(Time, Grd, id_N2_nblayer_base, wrk2_2d(:,:))
     endif

     if(id_N2_for_agm > 0) then 
         wrk4(:,:,:) = 0.0
         do k=1,nk
            wrk4(:,:,k) = grav*rho0r*abs(drhodz_zt(:,:,k))
         enddo
         call diagnose_3d(Time, Grd, id_N2_for_agm, wrk4(:,:,:))
     endif
  endif ! endif of agm_closure 


  ! set redi diffusivity 
  if(aredi_diffusivity_grid_scaling) then 
      do k=1,nk
         aredi_array(:,:,k) = aredi*aredi_grid_scaling(:,:)*Grd%tmask(:,:,k) 
      enddo
  endif
  if(aredi_equal_agm) then 
      aredi_array(:,:,:) = agm_array(:,:,:) 
  endif


  ! diagnostics 
  call diagnose_2d(Time, Grd, id_agm, agm_array(:,:,1))
  call diagnose_3d(Time, Grd, id_agm_3d, agm_array(:,:,:))
  call diagnose_2d(Time, Grd, id_aredi, agm_array(:,:,1))
  call diagnose_3d(Time, Grd, id_aredi_3d, agm_array(:,:,:))

  if (id_ksurf_blayer >  0) then 
     wrk1_2d(:,:) = ksurf_blayer(:,:)
     call diagnose_2d(Time, Grd, id_ksurf_blayer, wrk1_2d(:,:))
  endif 

  if(debug_this_module) then 
      call write_chksum_header(FILENAME, ' ', Time%model_time)
      call write_chksum_3d('checksum agm_array', agm_array(COMP,:)*Grd%tmask(COMP,:))
  endif


  call mpp_clock_end(id_clock_compute_diffusivity)

end subroutine compute_diffusivity
! </SUBROUTINE> NAME="compute_diffusivity"


!#######################################################################
! <SUBROUTINE NAME="transport_on_nrho_gm">
!
! <DESCRIPTION>
! Classify horizontal GM mass transport according to neutral density classes. 
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_nrho_gm, a remapping is 
! performed, rather than the binning done for trans_rho.  
!
! Code history 
! 2008: algorithm based (incorrectly) on transport_on_rho 
! 2009: algorithm corrected to be consistent with remapping 
!       used in tracer_on_rho algorithm
! </DESCRIPTION>
!
subroutine transport_on_nrho_gm (Time, Dens, tx_trans_gm, ty_trans_gm)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_gm
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_gm
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_rho, neutralrho_nk
  real    :: work(isd:ied,jsd:jed,size(Dens%neutralrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (transport_on_nrho_gm): module needs initialization ')
  endif 

  next_time = increment_time(Time%model_time, int(dtime), 0)

  if (need_data(id_tx_trans_nrho_gm,next_time) .or. need_data(id_ty_trans_nrho_gm,next_time)) then

      neutralrho_nk = size(Dens%neutralrho_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with neutralrho_ref < neutralrho(k=1),   work=0
      ! for (i,j) points with neutralrho_ref > neutralrho(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_gm from k-levels to neutralrho_nk-levels
      do k_rho=1,neutralrho_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%neutralrho_ref(k_rho) >  Dens%neutralrho(i,j,k)  ) then
                      if( Dens%neutralrho_ref(k_rho) <= Dens%neutralrho(i,j,k+1)) then 
                          W1= Dens%neutralrho_ref(k_rho)- Dens%neutralrho(i,j,k)
                          W2= Dens%neutralrho(i,j,k+1)  - Dens%neutralrho_ref(k_rho)
                          work(i,j,k_rho,1) = (tx_trans_gm(i,j,k+1)*W1 +tx_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                          work(i,j,k_rho,2) = (ty_trans_gm(i,j,k+1)*W1 +ty_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_rho=1,neutralrho_nk
         work(COMP,k_rho,1) = work(COMP,k_rho,1)*Grd%tmask(COMP,1)
         work(COMP,k_rho,2) = work(COMP,k_rho,2)*Grd%tmask(COMP,1)
      enddo

      if (id_tx_trans_nrho_gm > 0) then 
         call diagnose_3d(Time, Grd, id_tx_trans_nrho_gm, work(:,:,:,1), &
                          nk_lim=neutralrho_nk, use_mask=.false.)
      endif
      if (id_ty_trans_nrho_gm > 0) then 
         call diagnose_3d(Time, Grd, id_ty_trans_nrho_gm, work(:,:,:,2), &
                          nk_lim=neutralrho_nk, use_mask=.false.)
      endif

  endif

end subroutine transport_on_nrho_gm
! </SUBROUTINE> NAME="transport_on_nrho_gm"


!#######################################################################
! <SUBROUTINE NAME="transport_on_rho_gm">
!
! <DESCRIPTION>
! Classify horizontal GM mass transport according to potential density classes. 
!
! Algorithm based on linear interpolation of function on s-surfaces to 
! function on rho-surfaces.  
!
! Diagnostic makes sense when potrho is monotonically increasing with 
! depth, although the algorithm does not explicitly make this assumption.  
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_rho_gm, a remapping is 
! performed, rather than the binning done for trans_rho.  
!
! Code history 
!
! 2008: algorithm based (incorrectly) on transport_on_rho 
! 2009: algorithm corrected to be consistent with remapping 
!       used in tracer_on_rho algorithm
!
! </DESCRIPTION>
!
subroutine transport_on_rho_gm (Time, Dens, tx_trans_gm, ty_trans_gm)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_gm
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_gm
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_rho, potrho_nk
  real    :: work(isd:ied,jsd:jed,size(Dens%potrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (transport_on_rho_gm): module needs initialization ')
  endif 

  next_time = increment_time(Time%model_time, int(dtime), 0)
  if (need_data(id_tx_trans_rho_gm,next_time) .or. need_data(id_ty_trans_rho_gm,next_time)) then

      potrho_nk = size(Dens%potrho_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with potrho_ref < potrho(k=1),   work=0
      ! for (i,j) points with potrho_ref > potrho(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_gm from k-levels to krho-levels
      do k_rho=1,potrho_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%potrho_ref(k_rho) >  Dens%potrho(i,j,k)  ) then
                      if( Dens%potrho_ref(k_rho) <= Dens%potrho(i,j,k+1)) then 
                          W1= Dens%potrho_ref(k_rho)- Dens%potrho(i,j,k)
                          W2= Dens%potrho(i,j,k+1)  - Dens%potrho_ref(k_rho)
                          work(i,j,k_rho,1) = (tx_trans_gm(i,j,k+1)*W1 +tx_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                          work(i,j,k_rho,2) = (ty_trans_gm(i,j,k+1)*W1 +ty_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_rho=1,potrho_nk
         work(COMP,k_rho,1) = work(COMP,k_rho,1)*Grd%tmask(COMP,1)
         work(COMP,k_rho,2) = work(COMP,k_rho,2)*Grd%tmask(COMP,1)
      enddo

      call diagnose_3d(Time, Grd, id_tx_trans_rho_gm, work(:,:,:,1), &
                       nk_lim=potrho_nk, use_mask=.false.)
      call diagnose_3d(Time, Grd, id_ty_trans_rho_gm, work(:,:,:,2), &
                       nk_lim=potrho_nk, use_mask=.false.)
  endif


end subroutine transport_on_rho_gm
! </SUBROUTINE> NAME="transport_on_rho_gm"


!#######################################################################
! <SUBROUTINE NAME="transport_on_theta_gm">
!
! <DESCRIPTION>
! Classify horizontal GM mass transport according to potential temp classes. 
!
! Algorithm based on linear interpolation of function on s-surfaces to 
! function on rho-surfaces.  
!
! Diagnostic makes sense when potential temp is monotonically increasing 
! with depth, although the algorithm does not explicitly make this assumption.  
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_theta_gm, a remapping is 
! performed, rather than the binning done for trans_rho.  
!
! Code history 
!
! 2009: algorithm based on transport_on_rho_gm
!
! </DESCRIPTION>
!
subroutine transport_on_theta_gm (Time, Dens, Theta, tx_trans_gm, ty_trans_gm)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_prog_tracer_type), intent(in) :: Theta
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_gm
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_gm
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_theta, theta_nk, tau
  real    :: work(isd:ied,jsd:jed,size(Dens%potrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (transport_on_theta_gm): module needs initialization ')
  endif 

  tau = Time%tau
  next_time = increment_time(Time%model_time, int(dtime), 0)
  if (need_data(id_tx_trans_theta_gm,next_time) .or. need_data(id_ty_trans_theta_gm,next_time)) then

      theta_nk = size(Dens%theta_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with theta_ref > theta(k=1),   work=0
      ! for (i,j) points with theta_ref < theta(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_gm from k-levels to theta-levels
      ! note sign change in the if-tests relative to transport on rho and nrho
      do k_theta=1,theta_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%theta_ref(k_theta) <  Theta%field(i,j,k,tau) ) then
                      if( Dens%theta_ref(k_theta) >= Theta%field(i,j,k+1,tau)) then 
                          W1= -Dens%theta_ref(k_theta)  + Theta%field(i,j,k,tau)
                          W2= -Theta%field(i,j,k+1,tau) + Dens%theta_ref(k_theta)
                          work(i,j,k_theta,1) = (tx_trans_gm(i,j,k+1)*W1 +tx_trans_gm(i,j,k)*W2) &
                                                /(W1 + W2 + epsln)
                          work(i,j,k_theta,2) = (ty_trans_gm(i,j,k+1)*W1 +ty_trans_gm(i,j,k)*W2) &
                                                /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_theta=1,theta_nk
         work(COMP,k_theta,1) = work(COMP,k_theta,1)*Grd%tmask(COMP,1)
         work(COMP,k_theta,2) = work(COMP,k_theta,2)*Grd%tmask(COMP,1)
      enddo

      call diagnose_3d(Time, Grd, id_tx_trans_theta_gm, work(:,:,:,1), &
                       nk_lim=theta_nk, use_mask=.false.)
      call diagnose_3d(Time, Grd, id_ty_trans_theta_gm, work(:,:,:,2), &
                       nk_lim=theta_nk, use_mask=.false.)
  endif


end subroutine transport_on_theta_gm
! </SUBROUTINE> NAME="transport_on_theta_gm"


!#######################################################################
! <SUBROUTINE NAME="compute_eta_tend_gm90">
!
! <DESCRIPTION>
! Diagnose contribution to global mean sea level evolution arising  
! from analytic form of Gent-McWilliams scheme.
!
! This routine computes the diagnostic based on an analytic form of the 
! GM90 contribution.  The raw numerical form is computed inside the respective
! nphysics modules.  The raw numerical form is more accurate and thus
! recomended for purposes of sea level budgets.  
!
! Compute an averaged slope using tensor_31 and tensor_32.
! Then compute the neutral divergence of this slope vector, just
! as if each component of the slope was a scalar.  
!
! To avoid stencil issues with bottom cells, mask to zero 
! contributions from cells next to the bottom in either of the
! three directions.  
!
! Send output to diagnostic manager. 
! 
! Subroutine history: 
! Feb2010 version 1.0: Stephen.Griffies 
! 
! </DESCRIPTION>
!
subroutine compute_eta_tend_gm90 (Time, Thickness, Dens,  &
                                  dtheta_dx, dtheta_dy,   & 
                                  dsalt_dx , dsalt_dy ,   & 
                                  drhodzh,                &
                                  dtwedyt, dzwtr, delqc,  &
                                  tensor_31, tensor_32, agm_array)

  type(ocean_time_type),              intent(in) :: Time
  type(ocean_thickness_type),         intent(in) :: Thickness
  type(ocean_density_type),           intent(in) :: Dens
  real, dimension(isd:,jsd:,:),       intent(in) :: dtheta_dx
  real, dimension(isd:,jsd:,:),       intent(in) :: dtheta_dy
  real, dimension(isd:,jsd:,:),       intent(in) :: dsalt_dx
  real, dimension(isd:,jsd:,:),       intent(in) :: dsalt_dy
  real, dimension(isd:,jsd:,:,0:),    intent(in) :: drhodzh
  real, dimension(isd:,jsd:,0:),      intent(in) :: dtwedyt
  real, dimension(isd:,jsd:,0:),      intent(in) :: dzwtr
  real, dimension(isd:,jsd:,:,0:),    intent(in) :: delqc
  real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: tensor_31
  real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: tensor_32
  real, dimension(isd:,jsd:,:),       intent(in) :: agm_array

  integer :: i, j, k, kp1, tau
  real    :: extended_mask
  real    :: absslope, absslopex, absslopey
  real    :: drho_dz, drho_dz_rhor2, divergence

  if(.not. diagnose_eta_tend_gm90) return 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (compute_eta_tend_gm90): module needs initialization ')
  endif 

  call mpp_clock_begin(id_clock_eta_tend_gm90)

  tau = Time%tau
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  ! approximate the neutral slope at the T-point as the average 
  ! of the previously computed triad slopes contained in tensor_31 and tensor_32. 
  ! taper to zero in regions of steep slope.  
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           slopex(i,j,k) =  onefourth*  &
           (tensor_31(i,j,k,0,0) + tensor_31(i,j,k,0,1) + tensor_31(i,j,k,1,0) + tensor_31(i,j,k,1,1))
           slopey(i,j,k) =  onefourth*  &
           (tensor_32(i,j,k,0,0) + tensor_32(i,j,k,0,1) + tensor_32(i,j,k,1,0) + tensor_32(i,j,k,1,1))

           taper_fcn(i,j,k) = Grd%tmask(i,j,k)
           absslopex = abs(slopex(i,j,k))
           absslopey = abs(slopey(i,j,k))
           absslope  = max(absslopex,absslopey)
           if(absslope > smax_grad_gamma_scalar) then 
               taper_fcn(i,j,k) = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               slopex(i,j,k) = slopex(i,j,k)*taper_fcn(i,j,k)  
               slopey(i,j,k) = slopey(i,j,k)*taper_fcn(i,j,k)  
           endif

           agm_rho_slopex(i,j,k) = slopex(i,j,k)*Dens%rho(i,j,k,tau)*agm_array(i,j,k)            
           agm_rho_slopey(i,j,k) = slopey(i,j,k)*Dens%rho(i,j,k,tau)*agm_array(i,j,k)            

        enddo
     enddo
  enddo
  call mpp_update_domains(agm_rho_slopex(:,:,:), Dom%domain2d) 
  call mpp_update_domains(agm_rho_slopey(:,:,:), Dom%domain2d) 
 

  ! compute derivatives of agm_rho_slope vector  
  do k=1,nk
     kp1 = min(k+1,nk)
     wrk1(:,:,1)       = agm_rho_slopex(:,:,k)
     dslopex_dx(:,:,k) = FDX_T(wrk1(:,:,1))*FMX(Grd%tmask(:,:,k))
     wrk2(:,:,1)       = agm_rho_slopey(:,:,k)
     dslopey_dy(:,:,k) = FDY_T(wrk2(:,:,1))*FMY(Grd%tmask(:,:,k))
  enddo
  do k=1,nk
     kp1 = min(k+1,nk)
     dslopex_dz(:,:,k) = Grd%tmask(:,:,kp1)*dzwtr(:,:,k) &
          *(agm_rho_slopex(:,:,k)-agm_rho_slopex(:,:,kp1))
     dslopey_dz(:,:,k) = Grd%tmask(:,:,kp1)*dzwtr(:,:,k) &
          *(agm_rho_slopey(:,:,k)-agm_rho_slopey(:,:,kp1))
  enddo

  ! compute the i-gradient of agm_rho_slopex along a neutral direction  
  ! units are rho*diffusivity/length  
  call calc_gradx_gamma_scalar(Thickness,Dens%drhodT(:,:,:),Dens%drhodS(:,:,:), &
                               dtheta_dx,dsalt_dx,                              &
                               drhodzh,dtwedyt,delqc,                           &
                               dslopex_dx,dslopex_dz,gradx_gamma_slopex)

  ! compute the j-gradient of agm_rho_slopey along a neutral direction   
  ! units are rho*diffusivity/length  
  call calc_grady_gamma_scalar(Thickness,Dens%drhodT(:,:,:),Dens%drhodS(:,:,:), &
                               dtheta_dy,dsalt_dy,                              &
                               drhodzh,dtwedyt,delqc,                           &
                               dslopey_dy,dslopey_dz,grady_gamma_slopey)
  
  ! accumulate the vertical integral for the tendency 
  ! unconcerned with diagnosing contributions where slope is steep, 
  ! so taper the contributions to zero in those regions.  
  wrk1_2d(:,:) = 0.0  
  do k=1,nk
     kp1=min(nk,k+1)
     do j=jsc,jec
        do i=isc,iec

           extended_mask   = Grd%tmask(i,j,k)*              &
                             (Grd%tmask(i,j,kp1)+epsln)*    &
                             (Grd%tmask(i+1,j,k)+epsln)*    &
                             (Grd%tmask(i-1,j,k)+epsln)*    &
                             (Grd%tmask(i,j+1,k)+epsln)*    &
                             (Grd%tmask(i,j-1,k)+epsln)*    &
                             (Grd%tmask(i+1,j,kp1)+epsln)*  &
                             (Grd%tmask(i-1,j,kp1)+epsln)*  &
                             (Grd%tmask(i,j+1,kp1)+epsln)*  &
                             (Grd%tmask(i,j-1,kp1)+epsln)

           drho_dz         = abs(Dens%drhodz_zt(i,j,k))
           drho_dz_rhor2   = extended_mask*drho_dz/(Dens%rho(i,j,k,tau)**2 + epsln)
           divergence      = gradx_gamma_slopex(i,j,k) + grady_gamma_slopey(i,j,k)
           wrk1_2d(i,j)    = wrk1_2d(i,j) + taper_fcn(i,j,k)*Thickness%dzt(i,j,k)*drho_dz_rhor2*divergence 
 
        enddo
     enddo
  enddo

  if (id_eta_tend_gm90 > 0 .or. id_eta_tend_gm90_glob > 0) then 
      if(smooth_eta_tend_gm90) then 
          call mpp_update_domains (wrk1_2d(:,:), Dom%domain2d) 
          wrk1_2d(:,:) = S2D(wrk1_2d(:,:))
      endif
      if(id_eta_tend_gm90 > 0) then   
         call diagnose_2d(Time, Grd, id_eta_tend_gm90, -1.0*wrk1_2d(:,:))
      endif 
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_gm90_glob, wrk1_2d, -1.0*cellarea_r)
  endif

  call mpp_clock_end(id_clock_eta_tend_gm90)


end subroutine compute_eta_tend_gm90
! </SUBROUTINE> NAME="compute_eta_tend_gm90"



!#######################################################################
! <SUBROUTINE NAME="cabbeling_thermob_tendency">
!
! <DESCRIPTION>
! Compute tendencies from cabbeling and thermobaricity.
!
! To avoid stencil issues with bottom cells, mask to zero 
! contributions from cells next to the bottom in either of the
! three directions.  
!
! Set cabbeling and thermobaricity to zero in regions where 
! vertical stratification is gravitationally unstable.  The idea is 
! that in such regions, convective mixing will wash away the impact of 
! along-neutral diffusion, in which case cabbeling and thermobaricity 
! are not relevant anyhow.  
!
! Set the vertical stratification drho_dz to a negative number 
! with lower bound on its magnitude as set by epsln_drhodz_diagnostics.
!
! Send output to diagnostic manager. 
! 
! Subroutine history: 
! Jan2010 version 1.0: initial coding by Stephen.Griffies  
! Mar2010 version 2.0: tweaks on the division by drho_dz and strat_mask
! 
! </DESCRIPTION>
!
subroutine cabbeling_thermob_tendency (Time, Thickness, T_prog, Dens,   &
                                       dtheta_dx, dtheta_dy, dtheta_dz, & 
                                       dsalt_dx, dsalt_dy,              & 
                                       drhodzh,                         &
                                       dxtdtsn, dtwedyt, dzwtr, delqc,  &
                                       aredi_array)

  type(ocean_time_type),              intent(in) :: Time
  type(ocean_thickness_type),         intent(in) :: Thickness
  type(ocean_prog_tracer_type),       intent(in) :: T_prog(:)
  type(ocean_density_type),           intent(in) :: Dens
  real, dimension(isd:,jsd:,:),       intent(in) :: dtheta_dx
  real, dimension(isd:,jsd:,:),       intent(in) :: dtheta_dy
  real, dimension(isd:,jsd:,0:),      intent(in) :: dtheta_dz
  real, dimension(isd:,jsd:,:),       intent(in) :: dsalt_dx
  real, dimension(isd:,jsd:,:),       intent(in) :: dsalt_dy
  real, dimension(isd:,jsd:,:,0:),    intent(in) :: drhodzh
  real, dimension(isd:,jsd:,0:),      intent(in) :: dxtdtsn
  real, dimension(isd:,jsd:,0:),      intent(in) :: dtwedyt
  real, dimension(isd:,jsd:,0:),      intent(in) :: dzwtr
  real, dimension(isd:,jsd:,:,0:),    intent(in) :: delqc
  real, dimension(isd:,jsd:,:),       intent(in) :: aredi_array

  integer :: i, j, k, kp1, tau
  real    :: extended_mask
  real    :: drho_dz, strat_mask
  real    :: temporary 

  if(.not. diagnose_cabbeling_thermob) return 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (cabbeling_thermob_tendency): module needs initialization ')
  endif 

  call mpp_clock_begin(id_clock_cabbeling_thermob)

  tau = Time%tau

  ! get the cabbeling and thermobaricity parameters from ocean_density.F90 
  call calc_cabbeling_thermobaricity(Dens%rho(:,:,:,tau), T_prog(index_salt)%field(:,:,:,tau),          &
                                    T_prog(index_temp)%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:), &
                                    Dens%drhodT(:,:,:), Dens%drhodS(:,:,:), Dens%drhodP(:,:,:),   &
                                    cabbeling_param, thermobaric_param)

  ! compute gradient of the pressure field 
  call pressure_derivs(Dens%pressure_at_depth(:,:,:), dzwtr)

  ! compute the i-gradient of temperature and pressure along a neutral direction   
  call calc_gradx_gamma_scalar(Thickness,Dens%drhodT(:,:,:),Dens%drhodS(:,:,:), &
                               dtheta_dx,dsalt_dx,                              &
                               drhodzh,dtwedyt,delqc,                           &
                               dtheta_dx,dtheta_dz,gradx_gamma_temp)
  call calc_gradx_gamma_scalar(Thickness,Dens%drhodT(:,:,:),Dens%drhodS(:,:,:), &
                               dtheta_dx,dsalt_dx,                              &
                               drhodzh,dtwedyt,delqc,                           &
                               deriv_x_press,deriv_z_press,gradx_gamma_press)

  ! compute the j-gradient of temperature and pressure along a neutral direction   
  call calc_grady_gamma_scalar(Thickness,Dens%drhodT(:,:,:),Dens%drhodS(:,:,:), &
                               dtheta_dy,dsalt_dy,                              &
                               drhodzh,dxtdtsn,delqc,                           &
                               dtheta_dy,dtheta_dz,grady_gamma_temp)
  call calc_grady_gamma_scalar(Thickness,Dens%drhodT(:,:,:),Dens%drhodS(:,:,:), &
                               dtheta_dy,dsalt_dy,                              &
                               drhodzh,dxtdtsn,delqc,                           &
                               deriv_y_press,deriv_z_press,grady_gamma_press)


  wrk1(:,:,:)     = 0.0   ! tendency from cabbeling (1/sec)
  wrk2(:,:,:)     = 0.0   ! tendency from thermobaricity (1/sec)
  wrk3(:,:,:)     = 0.0   ! speed from cabbeling (m/sec)
  wrk4(:,:,:)     = 0.0   ! speed from thermobaricity (m/sec)
  wrk5(:,:,:)     = 0.0   ! (drho/dt) from cabbeling (kg/m^3)/s
  wrk6(:,:,:)     = 0.0   ! (drho/dt) from thermobaricity (kg/m^3)/s
  wrk1_2d(:,:)    = 0.0   ! vertical sum of cabbeling speed (m/s)
  wrk2_2d(:,:)    = 0.0   ! vertical sum of thermobaricity speed (m/s)
  wrk1_v(:,:,:,:) = 0.0   ! dianeutral mass transport (kg/s) from (cabbeling,thermobaricity)
  do k=1,nk
     kp1=min(nk,k+1)
     do j=jsc,jec
        do i=isc,iec

           extended_mask   = Grd%tmask(i,j,k)*              &
                             (Grd%tmask(i,j,kp1)+epsln)*    &
                             (Grd%tmask(i+1,j,k)+epsln)*    &
                             (Grd%tmask(i-1,j,k)+epsln)*    &
                             (Grd%tmask(i,j+1,k)+epsln)*    &
                             (Grd%tmask(i,j-1,k)+epsln)*    &
                             (Grd%tmask(i+1,j,kp1)+epsln)*  &
                             (Grd%tmask(i-1,j,kp1)+epsln)*  &
                             (Grd%tmask(i,j+1,kp1)+epsln)*  &
                             (Grd%tmask(i,j-1,kp1)+epsln)

           ! zero those points where stratification is gravitationally unstable
           strat_mask = 1.0
           if(Dens%drhodz_zt(i,j,k) > 0.0) strat_mask=0.0

           wrk1(i,j,k)     = extended_mask*strat_mask*aredi_array(i,j,k)*cabbeling_param(i,j,k) &
            *(gradx_gamma_temp(i,j,k)*gradx_gamma_temp(i,j,k)+grady_gamma_temp(i,j,k)*grady_gamma_temp(i,j,k))
           wrk2(i,j,k)     = extended_mask*strat_mask*aredi_array(i,j,k)*thermobaric_param(i,j,k) &
            *(gradx_gamma_temp(i,j,k)*gradx_gamma_press(i,j,k)+grady_gamma_temp(i,j,k)*grady_gamma_press(i,j,k))

           ! provide lower bound on abs(drho_dz), and ensure drho_dz < 0
           drho_dz = min(-epsln_drhodz_diagnostics,Dens%drhodz_zt(i,j,k))

           temporary      =  Grd%tmask(i,j,k)*Grd%dat(i,j)*Dens%rho(i,j,k,tau) &
                             *Thickness%rho_dzt(i,j,k,tau)/(epsln+drho_dz*Thickness%dzt(i,j,k))

           wrk1_v(i,j,k,1) = wrk1(i,j,k)*temporary
           wrk1_v(i,j,k,2) = wrk2(i,j,k)*temporary

           wrk3(i,j,k)     = wrk1(i,j,k)*Thickness%dzt(i,j,k)
           wrk4(i,j,k)     = wrk2(i,j,k)*Thickness%dzt(i,j,k)

           wrk5(i,j,k)     = wrk1(i,j,k)*Dens%rho(i,j,k,tau)
           wrk6(i,j,k)     = wrk2(i,j,k)*Dens%rho(i,j,k,tau)

           wrk1_2d(i,j)    = wrk1_2d(i,j) + wrk3(i,j,k) 
           wrk2_2d(i,j)    = wrk2_2d(i,j) + wrk4(i,j,k) 

        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_cabbeling_tend,        wrk1(:,:,:))
  call diagnose_3d(Time, Grd, id_thermobaric_tend,      wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_cabbeling_speed,       wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_thermobaric_speed,     wrk4(:,:,:))
  call diagnose_2d(Time, Grd, id_cabbeling_tend_intz,   wrk1_2d(:,:))
  call diagnose_2d(Time, Grd, id_thermobaric_tend_intz, wrk2_2d(:,:))

  call diagnose_sum(Time, Grd, Dom, id_cabbeling_tend_intz_glob, wrk1_2d, cellarea_r)
  call diagnose_sum(Time, Grd, Dom, id_thermobaric_tend_intz_glob, wrk2_2d, cellarea_r)

  if(wdianeutral_smooth) then 
     if(id_wdian_cabbeling > 0 .or. id_wdian_cabbeling_on_nrho > 0 .or. &
        id_wdian_thermob   > 0 .or. id_wdian_thermob_on_nrho   > 0) then   
          call mpp_update_domains (wrk1_v(:,:,:,1), Dom%domain2d) 
          call mpp_update_domains (wrk1_v(:,:,:,2), Dom%domain2d) 
          do k=1,nk
             wrk1_v(:,:,k,1) = S2D(wrk1_v(:,:,k,1))
             wrk1_v(:,:,k,2) = S2D(wrk1_v(:,:,k,2))
          enddo
     endif
  endif 

  call diagnose_3d(Time, Grd, id_neut_rho_cabbeling, wrk5(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_cabbeling_on_nrho, wrk5)

  call diagnose_3d(Time, Grd, id_wdian_cabbeling, wrk1_v(:,:,:,1))
  call diagnose_3d_rho(Time, Dens, id_wdian_cabbeling_on_nrho, wrk1_v(:,:,:,1))

  call diagnose_3d(Time, Grd, id_neut_rho_thermob, wrk6(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_thermob_on_nrho, wrk6)

  call diagnose_3d(Time, Grd, id_wdian_thermob, wrk1_v(:,:,:,2))
  call diagnose_3d_rho(Time, Dens, id_wdian_thermob_on_nrho, wrk1_v(:,:,:,2))

  if (id_tform_rho_cabbel_on_nrho > 0) then 
      wrk3(:,:,:)      = 0.0
      do k=1,nk
         wrk3(COMP,k) = wrk1(COMP,k)*Dens%watermass_factor(COMP,k) &
              *Dens%rho(COMP,k,tau)*Thickness%rho_dzt(COMP,k,tau)*Grd%dat(COMP)
      enddo 
      call diagnose_3d_rho(Time, Dens, id_tform_rho_cabbel_on_nrho, wrk3)
  endif
  if (id_tform_rho_thermb_on_nrho > 0) then 
      wrk3(:,:,:)      = 0.0
      do k=1,nk
         wrk3(COMP,k) = wrk2(COMP,k)*Dens%watermass_factor(COMP,k) &
              *Dens%rho(COMP,k,tau)*Thickness%rho_dzt(COMP,k,tau)*Grd%dat(COMP)
      enddo 
      call diagnose_3d_rho(Time, Dens, id_tform_rho_thermb_on_nrho, wrk3)
  endif
    

  call mpp_clock_end(id_clock_cabbeling_thermob)


end subroutine cabbeling_thermob_tendency
! </SUBROUTINE> NAME="cabbeling_thermob_tendency"


!#######################################################################
! <SUBROUTINE NAME="pressure_derivs">
!
! <DESCRIPTION>
! Compute the pressure derivatives for use in thermobaricity diagnostic.
! </DESCRIPTION>
!
subroutine pressure_derivs(pressure, dzwtr)

  real, dimension(isd:,jsd:,:),  intent(in) :: pressure 
  real, dimension(isd:,jsd:,0:), intent(in) :: dzwtr

  integer :: k
  integer :: kp1

  do k=1,nk
     kp1 = min(k+1,nk)
     wrk1(:,:,1)          = pressure(:,:,k)
     deriv_x_press(:,:,k) = FDX_T(wrk1(:,:,1))*FMX(Grd%tmask(:,:,k))
     deriv_y_press(:,:,k) = FDY_T(wrk1(:,:,1))*FMY(Grd%tmask(:,:,k))
  enddo

  do k=1,nk
     kp1 = min(k+1,nk)
     deriv_z_press(:,:,k) = Grd%tmask(:,:,kp1)*dzwtr(:,:,k)*(pressure(:,:,k) - pressure(:,:,kp1))
  enddo

end subroutine pressure_derivs
! </SUBROUTINE> NAME="pressure_derivs"



!#######################################################################
! <SUBROUTINE NAME="calc_gradx_gamma_scalar">
!
! <DESCRIPTION>
! Subroutine computes the i-gradient along neutral surfaces
! for a scalar field. For use in the cabbeling and thermobaricity 
! diagnostic.  Thus, when slope is steep, the full derivative is 
! tapered, rather than just the off-diagonal term. 
!
! We only compute neutral i-gradient for interior regions,
! since will be setting thermobaricity and cabbeling to zero at 
! top and bottom levels; do not trust the calculation at top 
! and bottom boundaries due to stencil truncations. 
!
! Algorithm slightly modified from ocean_nphysC calculation of 
! neutral diffusion x-flux.  
! 
! </DESCRIPTION>
!
subroutine calc_gradx_gamma_scalar(Thickness,drhodT,drhodS,         &
                                   dtheta_dx,dsalt_dx,    &
                                   drhodzh,dtwedyt,delqc,           &
                                   dscalar_dx,dscalar_dz,gradx_gamma_scalar)

  type(ocean_thickness_type),      intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodS
  real, dimension(isd:,jsd:,:),    intent(in)    :: dtheta_dx
  real, dimension(isd:,jsd:,:),    intent(in)    :: dsalt_dx
  real, dimension(isd:,jsd:,:,0:), intent(in)    :: drhodzh
  real, dimension(isd:,jsd:,0:),   intent(in)    :: dtwedyt
  real, dimension(isd:,jsd:,:,0:), intent(in)    :: delqc
  real, dimension(isd:,jsd:,:),    intent(in)    :: dscalar_dx
  real, dimension(isd:,jsd:,0:),   intent(in)    :: dscalar_dz
  real, dimension(isd:,jsd:,:),    intent(inout) :: gradx_gamma_scalar

  integer :: i, j, k, ip
  integer :: kr, kpkr
  real    :: tensor_11(isd:ied,jsd:jed,0:1,0:1)
  real    :: tensor_13(isd:ied,jsd:jed,0:1,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)
  real    :: taperA
  real    :: slope, absslope
  real    :: drhodx, drhodz 
  real    :: cell_mass_inv

  gradx_gamma_scalar(:,:,:) = 0.0

  do k=2,nk-1

     ! initialize arrays to zero 
     tensor_11(:,:,:,:) = 0.0
     tensor_13(:,:,:,:) = 0.0

     ! scalar field independent part of the calculation 
     do kr=0,1
        kpkr = min(k+kr,nk)
        do ip=0,1 
           do j=jsc,jec
              do i=isc-1,iec

                 drhodz   = drhodzh(i+ip,j,k,kr)
                 drhodx   = Grd%tmask(i+ip,j,kpkr)               &
                            *( drhodT(i+ip,j,k)*dtheta_dx(i,j,k) &
                              +drhodS(i+ip,j,k)*dsalt_dx(i,j,k))
                 slope    = -drhodx/drhodz
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax_grad_gamma_scalar) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                 else 
                     taperA = 1.0
                 endif

                 ! fill tensor components; note tapering applied to both terms  
                 tensor_11(i,j,ip,kr) = taperA
                 tensor_13(i,j,ip,kr) = taperA*slope

              enddo
           enddo
        enddo
     enddo

     ! scalar field dependent part of the calculation
     sumz(:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)          &
                 *(tensor_11(i,j,ip,kr)*dscalar_dx(i,j,k) + tensor_13(i,j,ip,kr)*dscalar_dz(i+ip,j,k-1+kr)) &
                 *min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo

     ! need to divide by mass of cell to get to a scalar gradient 
     do j=jsc,jec
        do i=isc-1,iec
           cell_mass_inv             = Grd%tmask(i,j,k)*Thickness%rho_dztr(i,j,k)*Grd%dater(i,j)
           gradx_gamma_scalar(i,j,k) = (sumz(i,j,0)+sumz(i,j,1))*cell_mass_inv 
        enddo
     enddo


  enddo ! enddo for k-loop

end subroutine calc_gradx_gamma_scalar
! </SUBROUTINE> NAME="calc_gradx_gamma_scalar"


!#######################################################################
! <SUBROUTINE NAME="calc_grady_gamma_scalar">
!
! <DESCRIPTION>
! Subroutine computes the y-gradient along neutral surfaces
! for a scalar field. For use in the cabbeling and thermobaricity 
! diagnostic.  Thus, when slope is steep, the full derivative is 
! tapered, rather than just the off-diagonal term. 
!
! We only compute neutral y-gradient for interior regions,
! since will be setting thermobaricity and cabbeling to zero at 
! top and bottom levels; do not trust the calculation at top 
! and bottom boundaries due to stencil truncations.  
!
! Algorithm slightly modified from ocean_nphysC calculation of 
! neutral diffusion y-flux.  
! 
! </DESCRIPTION>
!
subroutine calc_grady_gamma_scalar(Thickness,drhodT,drhodS,     &
                                   dtheta_dy,dsalt_dy,          &
                                   drhodzh,dxtdtsn,delqc, &
                                   dscalar_dy,dscalar_dz,grady_gamma_scalar)

  type(ocean_thickness_type),      intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodS
  real, dimension(isd:,jsd:,:),    intent(in)    :: dtheta_dy
  real, dimension(isd:,jsd:,:),    intent(in)    :: dsalt_dy
  real, dimension(isd:,jsd:,:,0:), intent(in)    :: drhodzh
  real, dimension(isd:,jsd:,0:),   intent(in)    :: dxtdtsn
  real, dimension(isd:,jsd:,:,0:), intent(in)    :: delqc
  real, dimension(isd:,jsd:,:),    intent(in)    :: dscalar_dy
  real, dimension(isd:,jsd:,0:),   intent(in)    :: dscalar_dz
  real, dimension(isd:,jsd:,:),    intent(inout) :: grady_gamma_scalar

  integer :: i, j, k, jq
  integer :: kr, kpkr
  real :: tensor_22(isd:ied,jsd:jed,0:1,0:1)
  real :: tensor_23(isd:ied,jsd:jed,0:1,0:1)
  real :: sumz(isd:ied,jsd:jed,0:1)
  real :: taperA
  real :: slope, absslope
  real :: drhody, drhodz 
  real :: cell_mass_inv

  grady_gamma_scalar(:,:,:) = 0.0

  do k=2,nk-1

     ! initialize arrays to zero 
     tensor_22(:,:,:,:) = 0.0
     tensor_23(:,:,:,:) = 0.0


     do kr=0,1
        kpkr = min(k+kr,nk)
        do jq=0,1  
           do j=jsc-1,jec
              do i=isc,iec

                 drhodz   = drhodzh(i,j+jq,k,kr) 
                 drhody   = Grd%tmask(i,j+jq,kpkr)               &
                            *( drhodT(i,j+jq,k)*dtheta_dy(i,j,k) &
                              +drhodS(i,j+jq,k)*dsalt_dy(i,j,k)) 
                 slope    = -drhody/drhodz
                 absslope = abs(slope)  

                 ! taper for steep slope regions 
                 if(absslope > smax_grad_gamma_scalar) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                 else 
                     taperA = 1.0
                 endif

                 ! fill tensor components; note tapering applied to both terms  
                 tensor_22(i,j,jq,kr) = taperA
                 tensor_23(i,j,jq,kr) = taperA*slope

              enddo
           enddo
        enddo
     enddo


     ! scalar field dependent part of the calculation
     sumz(:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                      *(tensor_22(i,j,jq,kr)*dscalar_dy(i,j,k) + tensor_23(i,j,jq,kr)*dscalar_dz(i,j+jq,k-1+kr)) &
                      * min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo

     ! need to divide by mass of cell to get to a scalar gradient 
     do j=jsc,jec
        do i=isc-1,iec
           cell_mass_inv             = Grd%tmask(i,j,k)*Thickness%rho_dztr(i,j,k)*Grd%datnr(i,j)
           grady_gamma_scalar(i,j,k) = (sumz(i,j,0)+sumz(i,j,1))*cell_mass_inv 
        enddo
     enddo


  enddo ! enddo for k-loop


end subroutine calc_grady_gamma_scalar
! </SUBROUTINE> NAME="calc_grady_gamma_scalar"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens, diagnose_gm_redi)

  type(ocean_time_type),    intent(in)    :: Time
  type(ocean_density_type), intent(in)    :: Dens
  logical                 , intent(inout) :: diagnose_gm_redi

  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag = .false. 


  id_neut_rho_nphysics = -1
  id_neut_rho_nphysics = register_diag_field ('ocean_model', 'neut_rho_nphysics', &
    Grd%tracer_axes(1:3), Time%model_time,                                        &
    'update of locally referenced potential density from time-explicit  nphysics',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_nphysics > 0) compute_watermass_diag=.true.

  id_wdian_rho_nphysics = -1
  id_wdian_rho_nphysics = register_diag_field ('ocean_model', 'wdian_rho_nphysics',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'dianeutral mass transport due to time-explicit nphysics',                     &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_nphysics > 0) compute_watermass_diag=.true.


  ! full contributions to neutral diffusion  
  id_neut_rho_ndiff = -1
  id_neut_rho_ndiff = register_diag_field ('ocean_model', 'neut_rho_ndiff',          &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'update of locally referenced potential density time-explicit neutral diffusion',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_ndiff > 0) compute_watermass_diag=.true.

  id_neut_rho_ndiff_on_nrho = -1
  id_neut_rho_ndiff_on_nrho = register_diag_field ('ocean_model', 'neut_rho_ndiff_on_nrho',                &
       Dens%neutralrho_axes(1:3), Time%model_time,                                                         &
       'update of locally ref potrho from time-explicit neutral diffusion as binned to neutral rho layers',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_ndiff = -1
  id_wdian_rho_ndiff = register_diag_field ('ocean_model', 'wdian_rho_ndiff',&
    Grd%tracer_axes(1:3), Time%model_time,                                   &
    'dianeutral mass transport due to time-explicit neutral diffusion',      &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_ndiff > 0) compute_watermass_diag=.true.

  id_wdian_rho_ndiff_on_nrho = -1
  id_wdian_rho_ndiff_on_nrho = register_diag_field ('ocean_model', 'wdian_rho_ndiff_on_nrho',         &
   Dens%neutralrho_axes(1:3), Time%model_time,                                                        &
   'dianeutral mass transport due to time-explicit neutral diffusion as binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_ndiff = -1
  id_tform_rho_ndiff = register_diag_field ('ocean_model', 'tform_rho_ndiff',                  & 
    Grd%tracer_axes(1:3), Time%model_time,                                                     &
    'watermass transform due to time-explicit neutral diffusion on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_ndiff > 0) compute_watermass_diag=.true.

  id_tform_rho_ndiff_on_nrho = -1
  id_tform_rho_ndiff_on_nrho = register_diag_field ('ocean_model', 'tform_rho_ndiff_on_nrho',   &
   Dens%neutralrho_axes(1:3), Time%model_time,                                                  &
   'watermass transform due to time-explicit neutral diffusion as binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_eta_tend_ndiff_tend= -1          
  id_eta_tend_ndiff_tend= register_diag_field ('ocean_model','eta_tend_ndiff_tend',          &
    Grd%tracer_axes(1:2), Time%model_time,                                                   &
    'numerical form of non-Bouss steric sea level tend from time explicit ndiffusion', 'm/s',&
    missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_ndiff_tend > 0) compute_watermass_diag=.true.

  id_eta_tend_ndiff_tend_glob= -1          
  id_eta_tend_ndiff_tend_glob= register_diag_field ('ocean_model', 'eta_tend_ndiff_tend_glob',   &
   Time%model_time,                                                                              &
   'global mean numerical form of non-bouss steric sea level tend from time explicit ndiffusion',&
   'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_ndiff_tend_glob > 0) compute_watermass_diag=.true.


  ! temp contributions to neutral diffusion  
  id_neut_temp_ndiff = -1
  id_neut_temp_ndiff = register_diag_field ('ocean_model', 'neut_temp_ndiff',   &
    Grd%tracer_axes(1:3), Time%model_time,                                      &
    'temp related update of locally ref potrho time-explicit neutral diffusion',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_ndiff > 0) compute_watermass_diag=.true.

  id_neut_temp_ndiff_on_nrho = -1
  id_neut_temp_ndiff_on_nrho = register_diag_field ('ocean_model', 'neut_temp_ndiff_on_nrho',                &
    Dens%neutralrho_axes(1:3), Time%model_time,                                                              &
    'temp related update of locally ref potrho from time-explicit neutral diff binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_temp_ndiff = -1
  id_wdian_temp_ndiff = register_diag_field ('ocean_model', 'wdian_temp_ndiff',   &
    Grd%tracer_axes(1:3), Time%model_time,                                        &
    'temp related dianeutral mass transport from time-explicit neutral diffusion',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_ndiff > 0) compute_watermass_diag=.true.

  id_wdian_temp_ndiff_on_nrho = -1
  id_wdian_temp_ndiff_on_nrho = register_diag_field ('ocean_model', 'wdian_temp_ndiff_on_nrho',               &
   Dens%neutralrho_axes(1:3), Time%model_time,                                                                &
   'temp related dianeutral mass transport from time-explicit neutral diffusion binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_temp_ndiff = -1
  id_tform_temp_ndiff = register_diag_field ('ocean_model', 'tform_temp_ndiff',         &
    Grd%tracer_axes(1:3), Time%model_time,                                              &
    'temp related watermass transform due to time-explicit neutral diffusion on levels',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_ndiff > 0) compute_watermass_diag=.true.

  id_tform_temp_ndiff_on_nrho = -1
  id_tform_temp_ndiff_on_nrho = register_diag_field ('ocean_model', 'tform_temp_ndiff_on_nrho',           &
   Dens%neutralrho_axes(1:3), Time%model_time,                                                            &
   'temp related watermass transform due to time-explicit neutral diffusion binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_ndiff_on_nrho > 0) compute_watermass_diag=.true.


  ! salinity contributions to neutral diffusion  
  id_neut_salt_ndiff = -1
  id_neut_salt_ndiff = register_diag_field ('ocean_model', 'neut_salt_ndiff',   &
    Grd%tracer_axes(1:3), Time%model_time,                                      &
    'salt related update of locally ref potrho time-explicit neutral diffusion',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_ndiff > 0) compute_watermass_diag=.true.

  id_neut_salt_ndiff_on_nrho = -1
  id_neut_salt_ndiff_on_nrho = register_diag_field ('ocean_model', 'neut_salt_ndiff_on_nrho',                &
    Dens%neutralrho_axes(1:3), Time%model_time,                                                              &
    'salt related update of locally ref potrho from time-explicit neutral diff binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_salt_ndiff = -1
  id_wdian_salt_ndiff = register_diag_field ('ocean_model', 'wdian_salt_ndiff',   &
    Grd%tracer_axes(1:3), Time%model_time,                                        &
    'salt related dianeutral mass transport from time-explicit neutral diffusion',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_ndiff > 0) compute_watermass_diag=.true.

  id_wdian_salt_ndiff_on_nrho = -1
  id_wdian_salt_ndiff_on_nrho = register_diag_field ('ocean_model', 'wdian_salt_ndiff_on_nrho',               &
   Dens%neutralrho_axes(1:3), Time%model_time,                                                                &
   'salt related dianeutral mass transport from time-explicit neutral diffusion binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_ndiff_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_salt_ndiff = -1
  id_tform_salt_ndiff = register_diag_field ('ocean_model', 'tform_salt_ndiff',         &
    Grd%tracer_axes(1:3), Time%model_time,                                              &
    'salt related watermass transform due to time-explicit neutral diffusion on levels',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_ndiff > 0) compute_watermass_diag=.true.

  id_tform_salt_ndiff_on_nrho = -1
  id_tform_salt_ndiff_on_nrho = register_diag_field ('ocean_model', 'tform_salt_ndiff_on_nrho',           &
   Dens%neutralrho_axes(1:3), Time%model_time,                                                            &
   'salt related watermass transform due to time-explicit neutral diffusion binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_ndiff_on_nrho > 0) compute_watermass_diag=.true.


  ! full contributions to GM transport 
  id_neut_rho_gm = -1
  id_neut_rho_gm = register_diag_field ('ocean_model', 'neut_rho_gm',  &
    Grd%tracer_axes(1:3), Time%model_time,                             &
    'update of locally referenced potential density from GM transport',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_gm > 0) compute_watermass_diag=.true.

  id_neut_rho_gm_on_nrho = -1
  id_neut_rho_gm_on_nrho = register_diag_field ('ocean_model', 'neut_rho_gm_on_nrho',&
    Dens%neutralrho_axes(1:3), Time%model_time,                                      &
    'update of locally ref potrho from GM transport as binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_gm = -1
  id_wdian_rho_gm = register_diag_field ('ocean_model', 'wdian_rho_gm',&
    Grd%tracer_axes(1:3), Time%model_time,                             &
    'dianeutral mass transport due to GM transport',                   &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_gm > 0) compute_watermass_diag=.true.

  id_wdian_rho_gm_on_nrho = -1
  id_wdian_rho_gm_on_nrho = register_diag_field ('ocean_model', 'wdian_rho_gm_on_nrho',&
    Dens%neutralrho_axes(1:3), Time%model_time,                                        &
    'dianeutral mass transport due to GM transport as binned to neutral rho layers',   &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_gm = -1
  id_tform_rho_gm = register_diag_field ('ocean_model', 'tform_rho_gm',     &
    Grd%tracer_axes(1:3), Time%model_time,                                  &
    'watermass transform due to GM transport on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_gm > 0) compute_watermass_diag=.true.

  id_tform_rho_gm_on_nrho = -1
  id_tform_rho_gm_on_nrho = register_diag_field ('ocean_model', 'tform_rho_gm_on_nrho',&
    Dens%neutralrho_axes(1:3), Time%model_time,                                        &
   'watermass transform due to GM transport as binned to neutral rho layers',          &
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_eta_tend_gm_tend= -1          
  id_eta_tend_gm_tend= register_diag_field ('ocean_model','eta_tend_gm_tend',         &
    Grd%tracer_axes(1:2), Time%model_time,                                            &
    'numerical form of non-Bouss steric sea level tend from GM-based tendency', 'm/s',&
    missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_gm_tend > 0) compute_watermass_diag=.true.

  id_eta_tend_gm_tend_glob= -1          
  id_eta_tend_gm_tend_glob= register_diag_field ('ocean_model', 'eta_tend_gm_tend_glob',  &
   Time%model_time,                                                                       &
   'global mean numerical form of non-bouss steric sea level tend from GM-based tendency',&
   'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_gm_tend_glob > 0) compute_watermass_diag=.true.


  ! temp contributions to GM transport
  id_neut_temp_gm = -1
  id_neut_temp_gm = register_diag_field ('ocean_model', 'neut_temp_gm',             &
    Grd%tracer_axes(1:3), Time%model_time,                                          &
    'temp related update of locally referenced potential density from GM transport',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_gm > 0) compute_watermass_diag=.true.

  id_neut_temp_gm_on_nrho = -1
  id_neut_temp_gm_on_nrho = register_diag_field ('ocean_model', 'neut_temp_gm_on_nrho',        &
    Dens%neutralrho_axes(1:3), Time%model_time,                                                &
    'temp related update of locally ref potrho from GM transport binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_temp_gm = -1
  id_wdian_temp_gm = register_diag_field ('ocean_model', 'wdian_temp_gm',&
    Grd%tracer_axes(1:3), Time%model_time,                               &
    'temp related dianeutral mass transport due to GM transport',        &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_gm > 0) compute_watermass_diag=.true.

  id_wdian_temp_gm_on_nrho = -1
  id_wdian_temp_gm_on_nrho = register_diag_field ('ocean_model', 'wdian_temp_gm_on_nrho',     &
    Dens%neutralrho_axes(1:3), Time%model_time,                                               & 
    'temp related dianeutral mass transport due to GM transport binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_temp_gm = -1
  id_tform_temp_gm = register_diag_field ('ocean_model', 'tform_temp_gm',&
    Grd%tracer_axes(1:3), Time%model_time,                               &
    'temp related watermass transform due to GM transport on levels',    &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_gm > 0) compute_watermass_diag=.true.

  id_tform_temp_gm_on_nrho = -1
  id_tform_temp_gm_on_nrho = register_diag_field ('ocean_model', 'tform_temp_gm_on_nrho',&
    Dens%neutralrho_axes(1:3), Time%model_time,                                          &
   'temp related watermass transform due to GM transport binned to neutral rho layers',  &
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_gm_on_nrho > 0) compute_watermass_diag=.true.


  ! salinity contributions to GM transport
  id_neut_salt_gm = -1
  id_neut_salt_gm = register_diag_field ('ocean_model', 'neut_salt_gm',             &
    Grd%tracer_axes(1:3), Time%model_time,                                          &
    'salt related update of locally referenced potential density from GM transport',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_gm > 0) compute_watermass_diag=.true.

  id_neut_salt_gm_on_nrho = -1
  id_neut_salt_gm_on_nrho = register_diag_field ('ocean_model', 'neut_salt_gm_on_nrho',        &
    Dens%neutralrho_axes(1:3), Time%model_time,                                                &
    'salt related update of locally ref potrho from GM transport binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_salt_gm = -1
  id_wdian_salt_gm = register_diag_field ('ocean_model', 'wdian_salt_gm',&
    Grd%tracer_axes(1:3), Time%model_time,                               &
    'salt related dianeutral mass transport due to GM transport',        &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_gm > 0) compute_watermass_diag=.true.

  id_wdian_salt_gm_on_nrho = -1
  id_wdian_salt_gm_on_nrho = register_diag_field ('ocean_model', 'wdian_salt_gm_on_nrho',     &
    Dens%neutralrho_axes(1:3), Time%model_time,                                               & 
    'salt related dianeutral mass transport due to GM transport binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_gm_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_salt_gm = -1
  id_tform_salt_gm = register_diag_field ('ocean_model', 'tform_salt_gm',&
    Grd%tracer_axes(1:3), Time%model_time,                               &
    'salt related watermass transform due to GM transport on levels',    &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_gm > 0) compute_watermass_diag=.true.

  id_tform_salt_gm_on_nrho = -1
  id_tform_salt_gm_on_nrho = register_diag_field ('ocean_model', 'tform_salt_gm_on_nrho',&
    Dens%neutralrho_axes(1:3), Time%model_time,                                          &
   'salt related watermass transform due to GM transport binned to neutral rho layers',  &
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_gm_on_nrho > 0) compute_watermass_diag=.true.


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_nphysics_util_mod w/ compute_watermass_diag=.true.'  
    diagnose_gm_redi=.true.
  endif 


end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from nphysics on the watermass transformation.
! For use with nphysicsA and nphysicsB. 
! </DESCRIPTION>
!
subroutine watermass_diag(Time, T_prog, Dens, tendency_redi_temp, tendency_redi_salt)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_density_type),      intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),  intent(in)    :: tendency_redi_temp
  real, dimension(isd:,jsd:,:),  intent(in)    :: tendency_redi_salt

  integer :: i,j,k,tau
  real    :: tmp1, tmp2
  real,  dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysics_util (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return 
  tau = Time%tau

  ! impacts from full neutral physics operator 
  if(id_neut_rho_nphysics > 0 .or. id_wdian_rho_nphysics > 0) then 
       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) &
                             +Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
             enddo
          enddo
       enddo
       call diagnose_3d(Time, Grd, id_neut_rho_nphysics, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_rho_nphysics, wrk3(:,:,:))
   endif

   ! impacts from neutral diffusion 
   if(id_neut_rho_ndiff      > 0 .or. id_neut_rho_ndiff_on_nrho   > 0 .or. &
      id_wdian_rho_ndiff     > 0 .or. id_wdian_rho_ndiff_on_nrho  > 0 .or. &
      id_tform_rho_ndiff     > 0 .or. id_tform_rho_ndiff_on_nrho  > 0 .or. &
      id_eta_tend_ndiff_tend > 0 .or. id_eta_tend_ndiff_tend_glob > 0) then 

       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       wrk4(:,:,:) = 0.0
       wrk5(:,:,:) = 0.0

       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk1(i,j,k) = Dens%drhodT(i,j,k)*tendency_redi_temp(i,j,k)  &
                             +Dens%drhodS(i,j,k)*tendency_redi_salt(i,j,k) 
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
                wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) ! for eta_tend 
             enddo
          enddo
       enddo

       call diagnose_3d(Time, Grd, id_neut_rho_ndiff, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_rho_ndiff, wrk3(:,:,:))
       call diagnose_3d(Time, Grd, id_tform_rho_ndiff, wrk4(:,:,:))
       call diagnose_3d_rho(Time, Dens, id_neut_rho_ndiff_on_nrho, wrk2)
       call diagnose_3d_rho(Time, Dens, id_wdian_rho_ndiff_on_nrho, wrk3)
       call diagnose_3d_rho(Time, Dens, id_tform_rho_ndiff_on_nrho, wrk4)
       if(id_eta_tend_ndiff_tend > 0 .or. id_eta_tend_ndiff_tend_glob > 0) then
           eta_tend(:,:) = 0.0
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
                 enddo
              enddo
           enddo
           call diagnose_2d(Time, Grd, id_eta_tend_ndiff_tend, eta_tend(:,:))
           call diagnose_sum(Time, Grd, Dom, id_eta_tend_ndiff_tend_glob, eta_tend, cellarea_r)
       endif

   endif   ! endif for ndiffuse


   ! temperature contributions to neutral diffusion 
   if(id_neut_temp_ndiff  > 0 .or. id_neut_temp_ndiff_on_nrho  > 0 .or. &
      id_wdian_temp_ndiff > 0 .or. id_wdian_temp_ndiff_on_nrho > 0 .or. &
      id_tform_temp_ndiff > 0 .or. id_tform_temp_ndiff_on_nrho > 0) then 

       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       wrk4(:,:,:) = 0.0

       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk1(i,j,k) = Dens%drhodT(i,j,k)*tendency_redi_temp(i,j,k)
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
             enddo
          enddo
       enddo

       call diagnose_3d(Time, Grd, id_neut_temp_ndiff, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_temp_ndiff, wrk3(:,:,:))
       call diagnose_3d(Time, Grd, id_tform_temp_ndiff, wrk4(:,:,:))
       call diagnose_3d_rho(Time, Dens, id_neut_temp_ndiff_on_nrho, wrk2)
       call diagnose_3d_rho(Time, Dens, id_wdian_temp_ndiff_on_nrho, wrk3)
       call diagnose_3d_rho(Time, Dens, id_tform_temp_ndiff_on_nrho, wrk4)
   endif   ! endif for temp_ndiffuse


   ! salinity contributions to neutral diffusion 
   if(id_neut_salt_ndiff  > 0 .or. id_neut_salt_ndiff_on_nrho  > 0 .or. &
      id_wdian_salt_ndiff > 0 .or. id_wdian_salt_ndiff_on_nrho > 0 .or. &
      id_tform_salt_ndiff > 0 .or. id_tform_salt_ndiff_on_nrho > 0) then 

       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       wrk4(:,:,:) = 0.0

       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk1(i,j,k) = Dens%drhodS(i,j,k)*tendency_redi_salt(i,j,k) 
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
             enddo
          enddo
       enddo

       call diagnose_3d(Time, Grd, id_neut_salt_ndiff, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_salt_ndiff, wrk3(:,:,:))
       call diagnose_3d(Time, Grd, id_tform_salt_ndiff, wrk4(:,:,:))
       call diagnose_3d_rho(Time, Dens, id_neut_salt_ndiff_on_nrho, wrk2)
       call diagnose_3d_rho(Time, Dens, id_wdian_salt_ndiff_on_nrho, wrk3)
       call diagnose_3d_rho(Time, Dens, id_tform_salt_ndiff_on_nrho, wrk4)
   endif   ! endif for salt_ndiffuse



   ! impacts from GM skewsion 
   if(id_neut_rho_gm      > 0 .or. id_neut_rho_gm_on_nrho   > 0 .or. &
      id_wdian_rho_gm     > 0 .or. id_wdian_rho_gm_on_nrho  > 0 .or. &
      id_tform_rho_gm     > 0 .or. id_tform_rho_gm_on_nrho  > 0 .or. &
      id_eta_tend_gm_tend > 0 .or. id_eta_tend_gm_tend_glob > 0) then 

       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       wrk4(:,:,:) = 0.0
       wrk5(:,:,:) = 0.0

       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                tmp1 = T_prog(index_temp)%wrk1(i,j,k) - tendency_redi_temp(i,j,k)
                tmp2 = T_prog(index_salt)%wrk1(i,j,k) - tendency_redi_salt(i,j,k)
                wrk1(i,j,k) = Dens%drhodT(i,j,k)*tmp1 + Dens%drhodS(i,j,k)*tmp2
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
                wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) ! for eta_tend 
             enddo
          enddo
       enddo

       call diagnose_3d(Time, Grd, id_neut_rho_gm, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_rho_gm, wrk3(:,:,:))
       call diagnose_3d(Time, Grd, id_tform_rho_gm, wrk4(:,:,:))
       call diagnose_3d_rho(Time, Dens, id_neut_rho_gm_on_nrho, wrk2)
       call diagnose_3d_rho(Time, Dens, id_wdian_rho_gm_on_nrho, wrk3)
       call diagnose_3d_rho(Time, Dens, id_tform_rho_gm_on_nrho, wrk4)
       if(id_eta_tend_gm_tend > 0 .or. id_eta_tend_gm_tend_glob > 0) then
           eta_tend(:,:) = 0.0
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
                 enddo
              enddo
           enddo
           call diagnose_2d(Time, Grd, id_eta_tend_gm_tend, eta_tend(:,:))
           call diagnose_sum(Time, Grd, Dom, id_eta_tend_gm_tend_glob, eta_tend, cellarea_r)
       endif


   endif   ! endif for GM 



   ! impacts from GM skewsion due to temperature contributions 
   if(id_neut_temp_gm  > 0 .or. id_neut_temp_gm_on_nrho  > 0 .or. &
      id_wdian_temp_gm > 0 .or. id_wdian_temp_gm_on_nrho > 0 .or. &
      id_tform_temp_gm > 0 .or. id_tform_temp_gm_on_nrho > 0) then 

       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       wrk4(:,:,:) = 0.0

       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                tmp1 = T_prog(index_temp)%wrk1(i,j,k) - tendency_redi_temp(i,j,k)
                wrk1(i,j,k) = Dens%drhodT(i,j,k)*tmp1 
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
             enddo
          enddo
       enddo

       call diagnose_3d(Time, Grd, id_neut_temp_gm, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_temp_gm, wrk3(:,:,:))
       call diagnose_3d(Time, Grd, id_tform_temp_gm, wrk4(:,:,:))
       call diagnose_3d_rho(Time, Dens, id_neut_temp_gm_on_nrho, wrk2)
       call diagnose_3d_rho(Time, Dens, id_wdian_temp_gm_on_nrho, wrk3)
       call diagnose_3d_rho(Time, Dens, id_tform_temp_gm_on_nrho, wrk4)
   endif   ! endif for GM from temperature contributions 


   ! impacts from GM skewsion due to salinity contributions 
   if(id_neut_salt_gm  > 0 .or. id_neut_salt_gm_on_nrho  > 0 .or. &
      id_wdian_salt_gm > 0 .or. id_wdian_salt_gm_on_nrho > 0 .or. &
      id_tform_salt_gm > 0 .or. id_tform_salt_gm_on_nrho > 0) then 

       wrk1(:,:,:) = 0.0
       wrk2(:,:,:) = 0.0
       wrk3(:,:,:) = 0.0
       wrk4(:,:,:) = 0.0

       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                tmp2        = T_prog(index_salt)%wrk1(i,j,k) - tendency_redi_salt(i,j,k)
                wrk1(i,j,k) = Dens%drhodS(i,j,k)*tmp2
                wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
             enddo
          enddo
       enddo

       call diagnose_3d(Time, Grd, id_neut_salt_gm, wrk2(:,:,:))
       call diagnose_3d(Time, Grd, id_wdian_salt_gm, wrk3(:,:,:))
       call diagnose_3d(Time, Grd, id_tform_salt_gm, wrk4(:,:,:))
       call diagnose_3d_rho(Time, Dens, id_neut_salt_gm_on_nrho, wrk2)
       call diagnose_3d_rho(Time, Dens, id_wdian_salt_gm_on_nrho, wrk3)
       call diagnose_3d_rho(Time, Dens, id_tform_salt_gm_on_nrho, wrk4)
   endif   ! endif for GM from salinity contributions 



end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_ndiffuse">
! <DESCRIPTION>
!
!  Diagnose watermass transformation from neutral diffusion. 
!  For use with the nphysicsC scheme.  
!
! </DESCRIPTION>
subroutine watermass_diag_ndiffuse(Time, Dens, tendency_redi_temp, tendency_redi_salt)

  type(ocean_time_type),         intent(in)  :: Time
  type(ocean_density_type),      intent(in)  :: Dens
  real, dimension(isd:,jsd:,:),  intent(in)  :: tendency_redi_temp
  real, dimension(isd:,jsd:,:),  intent(in)  :: tendency_redi_salt

  integer :: i,j,k,tau
  real,  dimension(isd:ied,jsd:jed) :: eta_tend

  if(.not. compute_watermass_diag) return 

  tau = Time%tau

  if (id_neut_rho_ndiff      > 0 .or. id_neut_rho_ndiff_on_nrho   > 0  .or.  &
      id_wdian_rho_ndiff     > 0 .or. id_wdian_rho_ndiff_on_nrho  > 0  .or.  & 
      id_tform_rho_ndiff     > 0 .or. id_tform_rho_ndiff_on_nrho  > 0  .or.  &
      id_eta_tend_ndiff_tend > 0 .or. id_eta_tend_ndiff_tend_glob > 0) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*(Dens%drhodT(i,j,k)*tendency_redi_temp(i,j,k) &
                                              +Dens%drhodS(i,j,k)*tendency_redi_salt(i,j,k))
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
               wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_ndiff,wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_ndiff, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_ndiff, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_ndiff_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_ndiff_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_ndiff_on_nrho, wrk4)

     if(id_eta_tend_ndiff_tend > 0 .or. id_eta_tend_ndiff_tend_glob > 0) then
         eta_tend(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_eta_tend_ndiff_tend, eta_tend(:,:))
         call diagnose_sum(Time, Grd, Dom, id_eta_tend_ndiff_tend_glob, eta_tend, cellarea_r)
     endif

  endif !endif for the neut_rho, wdian, and tform diagnostics iftest  


  ! ndiffuse scheme contributions from temperature 
  if (id_neut_temp_ndiff  > 0 .or. id_neut_temp_ndiff_on_nrho  > 0  .or.  &
      id_wdian_temp_ndiff > 0 .or. id_wdian_temp_ndiff_on_nrho > 0  .or.  & 
      id_tform_temp_ndiff > 0 .or. id_tform_temp_ndiff_on_nrho > 0) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*tendency_redi_temp(i,j,k)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_temp_ndiff,wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_temp_ndiff, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_temp_ndiff, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_temp_ndiff_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_temp_ndiff_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_temp_ndiff_on_nrho, wrk4)
  endif !endif for the neut_rho, wdian, and tform diagnostics from temperature effects 


  ! ndiffuse scheme contributions from salinity 
  if (id_neut_salt_ndiff  > 0 .or. id_neut_salt_ndiff_on_nrho  > 0  .or.  &
      id_wdian_salt_ndiff > 0 .or. id_wdian_salt_ndiff_on_nrho > 0  .or.  & 
      id_tform_salt_ndiff > 0 .or. id_tform_salt_ndiff_on_nrho > 0) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodS(i,j,k)*tendency_redi_salt(i,j,k)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_salt_ndiff,wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_salt_ndiff, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_salt_ndiff, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_salt_ndiff_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_salt_ndiff_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_salt_ndiff_on_nrho, wrk4)
  endif !endif for the neut_rho, wdian, and tform diagnostics from salinity effects 



end subroutine watermass_diag_ndiffuse
! </SUBROUTINE> NAME="watermass_diag_ndiffuse"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_sdiffuse">
! <DESCRIPTION>
!
!  Diagnose watermass transformation from skew diffusion. 
!  For use with the neutral physics C scheme.  
!
! </DESCRIPTION>
subroutine watermass_diag_sdiffuse(Time, Dens, tendency_gm_temp, tendency_gm_salt) 

  type(ocean_time_type),         intent(in)  :: Time
  type(ocean_density_type),      intent(in)  :: Dens
  real, dimension(isd:,jsd:,:),  intent(in)  :: tendency_gm_temp
  real, dimension(isd:,jsd:,:),  intent(in)  :: tendency_gm_salt

  integer :: i,j,k,tau
  real,  dimension(isd:ied,jsd:jed) :: eta_tend

  if(.not. compute_watermass_diag) return 

  tau = Time%tau  

  if (id_neut_rho_gm      > 0 .or. id_neut_rho_gm_on_nrho   > 0  .or. &
      id_wdian_rho_gm     > 0 .or. id_wdian_rho_gm_on_nrho  > 0  .or. & 
      id_tform_rho_gm     > 0 .or. id_tform_rho_gm_on_nrho  > 0  .or. &
      id_eta_tend_gm_tend > 0 .or. id_eta_tend_gm_tend_glob > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*(Dens%drhodT(i,j,k)*tendency_gm_temp(i,j,k) &
                                              +Dens%drhodS(i,j,k)*tendency_gm_salt(i,j,k))
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
               wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend 
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_gm,wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_gm, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_gm, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_gm_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_gm_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_gm_on_nrho, wrk4)
     if(id_eta_tend_gm_tend > 0 .or. id_eta_tend_gm_tend_glob > 0) then
         eta_tend(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_eta_tend_gm_tend, eta_tend(:,:))
         call diagnose_sum(Time, Grd, Dom, id_eta_tend_gm_tend_glob, eta_tend, cellarea_r)
     endif

  endif !endif for the neut_rho, wdian, and tform diagnostics iftest  



  ! temperature contributions
  if (id_neut_temp_gm  > 0 .or. id_neut_temp_gm_on_nrho  > 0  .or.  &
      id_wdian_temp_gm > 0 .or. id_wdian_temp_gm_on_nrho > 0  .or.  & 
      id_tform_temp_gm > 0 .or. id_tform_temp_gm_on_nrho > 0) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*tendency_gm_temp(i,j,k)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_temp_gm, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_temp_gm, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_temp_gm, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_temp_gm_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_temp_gm_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_temp_gm_on_nrho, wrk4)
  endif !endif for the neut_rho, wdian, and tform diagnostics from temperature effects  



  ! salinity contributions
  if (id_neut_salt_gm  > 0 .or. id_neut_salt_gm_on_nrho  > 0  .or.  &
      id_wdian_salt_gm > 0 .or. id_wdian_salt_gm_on_nrho > 0  .or.  & 
      id_tform_salt_gm > 0 .or. id_tform_salt_gm_on_nrho > 0) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodS(i,j,k)*tendency_gm_salt(i,j,k)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_salt_gm, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_salt_gm, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_salt_gm, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_salt_gm_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_salt_gm_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_salt_gm_on_nrho, wrk4)
  endif !endif for the neut_rho, wdian, and tform diagnostics from salinity effects  


end subroutine watermass_diag_sdiffuse
! </SUBROUTINE> NAME="watermass_diag_sdiffuse"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_util_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_nphysics_util_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  call save_restart(Nphysics_util_restart, time_stamp)

end subroutine ocean_nphysics_util_restart
! </SUBROUTINE> NAME="ocean_nphysics_util_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_coeff_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysics_coeff_end(Time, agm_array, aredi_array, &
                    rossby_radius, rossby_radius_raw, bczone_radius)

  type(ocean_time_type),        intent(in) :: Time
  real, dimension(isd:,jsd:,:), intent(in) :: agm_array
  real, dimension(isd:,jsd:,:), intent(in) :: aredi_array
  real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
  real, dimension(isd:,jsd:),   intent(in) :: rossby_radius_raw
  real, dimension(isd:,jsd:),   intent(in) :: bczone_radius 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  write(stdoutunit,*)' ' 
  write(stdoutunit,*) 'From ocean_nphysics_coeff_end: ending chksum'
  call write_timestamp(Time%model_time)

  call write_chksum_3d('ending agm_array', agm_array(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('ending aredi_array', aredi_array(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_2d('ending rossby_radius', rossby_radius(COMP)*Grd%tmask(COMP,1))
  call write_chksum_2d('ending rossby_radius_raw', rossby_radius_raw(COMP)*Grd%tmask(COMP,1))
  call write_chksum_2d('ending bczone_radius', bczone_radius(COMP)*Grd%tmask(COMP,1))

end subroutine ocean_nphysics_coeff_end
! </SUBROUTINE> NAME="ocean_nphysics_coeff_end"


end module ocean_nphysics_util_mod
