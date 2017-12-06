module ocean_density_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> Russell Fiedler 
!</CONTACT>
!
!<OVERVIEW>
! Compute density and related quantities.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module computes the in-situ density and its partial derivatives with 
! respect to conservative temperature or potential temperature, and with 
! respect to salinity.  
!
! There are three basic means for performing this calculation.
!
! A/ Linear equation for use in idealized studies
! 
! This equation renders density a linear function of potential 
! temperature and salinity.  All nonlinearities are ignored, as are  
! pressure effects. 
!
! The valid range for theta and salinity arbitrary for the 
! linear equation of state. 
!
! B/ pre-TEOS10 method: this method uses density as a rational 
! polynomial function of potential temperature, practical salinity,
! and gauge pressure.  There is also an implementation that computes
! density as a function of conservative temperature rather than 
! potential temperature.  The equation of state is based on that 
! documented in Jackett, McDougall, Feistel, Wright, and Griffies(2006).  
!
! This equation of state is valid over a "cone-shaped" range 
! corresponding to 
!
! 0psu <= salinity <= 40 psu
!
! -3C <= theta <= 40C    "theta" = either conservative or potential temp	 
!
! 0dbar <= pressure <= 8000dbar 
!
! with the cone getting smaller in the deeper ocean where 
! theta and salinity vary over a smaller range.  
!
!  Input variables are the following:
!
!  --salinity in psu or g/kg 
!  --conservative temperature or potential temperature (theta) in deg C
!  --pressure in dbars  (1bar = 10dbar = 10^5 Newton/m^2 = 10^5 Pascals). 
!
!  Note that in the ocean, pressure increases roughly by 1dbar for each meter depth.
!  Also note that pressure is the "sea pressure", which is the absolute pressure
!  minus the pressure of a standard atmosphere, which is 10.1325 dbars.
!
! check values                                                        <BR/>
!
!  for "theta" = conservative temperature 
!  rho(s=20psu,theta=20C,p=1000dbar)   = 1017.842890411975 (kg/m^3)   <BR/>
!  alpha(s=20psu,theta=20C,p=1000dbar) = 2.436057013634649e-4 (1/C)   <BR/>
!  beta(s=20psu,theta=20C,p=1000dbar)  = 7.314818108935248e-4 (1/psu) <BR/>
!
!  for "theta" = potential temperature 
!  rho(s=20psu,theta=20C,p=1000dbar)   = 1017.728868019642 (kg/m^3)   <BR/>
!  alpha(s=20psu,theta=20C,p=1000dbar) = 2.525481286927133e-4 (1/C)   <BR/>
!  beta(s=20psu,theta=20C,p=1000dbar)  = 7.379638527217575e-4 (1/psu) <BR/> 
!
! This equation of state should be suitable for purposes of realistic 
! global ocean climate modeling. 
!
! C/ TEOS10 method: this method makes use of the recommendations from 
! the SCOR working group on seawater thermodynamics, 2010.  Here, density
! is a function of conservative temperature and absolute salinity. 
! The equation is valid from 0 g/kg salinty to a very large value.
!
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
!  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2011:  A
!   computationally efficient 48-term expression for the density of
!   seawater in terms of Conservative Temperature, and related properties
!   of seawater.  To be submitted to Ocean Science.
! </REFERENCE>
!
! <REFERENCE>
! IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
! seawater 2010: Calculation and use of thermodynamic properties. 
! Intergovernmental Oceanographic Commission, Manuals and Guides No.
! 56, UNESCO (English), 196 pp.
! </REFERENCE>
!
! <REFERENCE>
! D. Iudicone, G. Madec, and T.J. McDougall (2008)
! Water-mass transformations in a neutral density framework and the
! key role of light penetration.  JPO vol 38, pages 1357-1376.
! </REFERENCE>
!
! <REFERENCE>
! Jackett, McDougall, Feistel, Wright, and Griffies (2006)
! Algorithms for density, potential temperature, conservative
! temperature, and freezing temperature of seawater.  
! Journal of Atmospheric and Oceanic Technology, 2006, 
! volume 23, pages 1709-1728.
! </REFERENCE>
!
! <REFERENCE>
! McDougall and Jackett (2005)
! The material derivative of neutral density
! Journal of Marine Research, vol 63, pages 159-185.
! </REFERENCE>
!
! <REFERENCE>
! Feistel (2003), A new extended Gibbs thermodynamic potential 
! of seawater. Progress in Oceanography. vol 58, pages 43-114.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison,  R.C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! <REFERENCE>
! T. McDougall (1987)
! Cabbeling, Thermobaricity, and water mass conversion
! JGR vol 92, pages 5448-5464
! </REFERENCE>
!
! <NOTE>
!
! Density is computed as a function of conservative temperature (degC) 
! or potential temperature (degC), salinity (psu or g/kg), and pressure (dbar).
! The pressure contribution includes that from the free surface height 
! and the applied atmospheric and/or sea ice pressure.  However, it is referenced
! to standard atmosphere, so that we use the "gauge" pressure rather than the 
! full in-situ pressure.  
!
! For vert_coordinate==GEOPOTENTIAL, ZSTAR, or ZSIGMA, baroclinic component of
! hydrostatic pressure is not known until density is known.  In this case,
! the baroclinic pressure contribution to density is lagged by a time step.  
! rho(tau) = rho[theta(tau),s(tau), p_atm(tau) + p_fs(tau) + p_baroclinic(tau-1)].  
! This issue does not arise when using vert_coordinate=PRESSURE, PSTAR, or PSIGMA.
!
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_density_nml">
!
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="press_standard" UNITS="dbar" TYPE="real">
!  Standard atmospheric pressure (dbar).  The realistic 
!  EOS used in MOM requires "sea pressure" as an argument
!  rather than absolute pressure.  Sea pressure is 
!  absolute pressure minus a standard atmospheric pressure 
!  of 10.1325dbar.  
!
!  For models that do have a realistic atmospheric loading, then it
!  is appropriate to remove 10.1325dbar prior to computing the EOS.
!  For those cases with zero atmospheric pressure, then it is not
!  necessary to remove the standard atmosphere.  The default for the 
!  press_standard is 0.0dbar.   
!  </DATA> 
!
!  <DATA NAME="t_test" UNITS="C" TYPE="real"> 			
!  Conservative temperature or potential temperature for 
!  testing the EOS.
!  </DATA> 
!  <DATA NAME="s_test" UNITS="psu or g/kg" TYPE="real">
!  Salinity for testing the EOS.
!  </DATA> 
!  <DATA NAME="p_test" UNITS="dbar" TYPE="real">
!  Sea pressure for testing the EOS.
!  </DATA> 
!
!  <DATA NAME="tn_test" UNITS="C" TYPE="real"> 			
!  Conservative temperature or potential temperature for 
!  testing the equation for neutral density.
!  </DATA> 
!  <DATA NAME="sn_test" UNITS="psu or g/kg" TYPE="real">
!  Salinity the equation for neutral density.
!  </DATA> 
!
!  <DATA NAME="eos_teos10" TYPE="logical">
!  Set to true to use TEOS-10 equation of state, which 
!  is a function of conservative temperature and absolute
!  salinity.  
!  Default eos_teos10=.false.
!  </DATA>
!  <DATA NAME="eos_preteos10" TYPE="logical">
!  Set to true to use pre-TEOS-10 equation of state, which 
!  is a function of potential temperature and practical salinity, 
!  or conservative temperature and practical salinity.
!  Default eos_preteos10=.false.
!  </DATA>
!  <DATA NAME="eos_linear" TYPE="logical">
!  Set to true to use an idealized linear equation of state, which 
!  has no pressure dependence, and is a linear function of salinity
!  and temperature. 
!  Default eos_linear=.false.
!  </DATA>
!  <DATA NAME="alpha_linear_eos" TYPE="real">
!  Constant "thermal expansion coefficient" for linear EOS 
!  rho = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
!  </DATA>
!  <DATA NAME="beta_linear_eos" TYPE="real">
!  Constant "saline contraction coefficient" for linear EOS 
!  rho = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
!  </DATA>
!
!  <DATA NAME="potrho_press" UNITS="dbar" TYPE="real">
!  Reference sea pressure for computing diagnostic potential density
!  of use for computing diagnostics with potential density.  
!  Default potrho_press=2000.0
!  </DATA> 
!  <DATA NAME="potrho_min" UNITS="kg/m^3" TYPE="real">
!  Minimum potential density used to partition vertical according 
!  to potential density.  
!  </DATA> 
!  <DATA NAME="potrho_max" UNITS="kg/m^3" TYPE="real">
!  Maximum potential density used to partition vertical according 
!  to potential density.
!  </DATA> 
!
!  <DATA NAME="neutral_density_omega" TYPE="logical">
!  Set to true to compute the neutral density according to 
!  the omega method based on Klocker and McDougall. 
!  This approach has not yet been coded. Presently as a 
!  placeholder we use potential density referenced to 2000dbar.
!  Default neutral_density_omega=.false.
!  </DATA>
!  <DATA NAME="neutral_density_potrho" TYPE="logical">
!  Set to true to compute the neutral density as just 
!  a selected potential density, set according to potrho_press. 
!  Since the neutral_density_omega approach has yet to be coded,
!  we only have the neutral_density_potrho option to choose from
!  at this time.  
!  Default neutral_density_potrho=.true.
!  </DATA>
!  <DATA NAME="neutral_density_theta" TYPE="logical">
!  Set to true to use temperature instead of neutral density as the
!  binning variable for water-mass diagnostics.
!  Default neutral_density_theta=.false.
!  </DATA>
!  <DATA NAME="neutralrho_min" UNITS="kg/m^3" TYPE="real">
!  Minimum neutral density used to partition vertical according 
!  to rational polynomial approximation to neutral density.  
!  </DATA> 
!  <DATA NAME="neutralrho_max" UNITS="kg/m^3" TYPE="real">
!  Maximum neutral density used to partition vertical according 
!  to rational polynomial approximation to neutral density.  
!  </DATA> 
!
!  <DATA NAME="theta_min" UNITS="C" TYPE="real">
!  Minimum conservative temperature or potential temperature used to 
!  partition vertical according to temperature.  
!  </DATA> 
!  <DATA NAME="theta_max" UNITS="C" TYPE="real">
!  Maximum conservative temperature or potential temperature used to
!  partition vertical according to temperature. 
!  </DATA> 
!
!  <DATA NAME="layer_nk" TYPE="integer">
!  Number of classes used to partition vertical according to potential density,
!  conservative temperature, or potential temperature. Used for diagnostics. 
!  </DATA> 
!
!  <DATA NAME="buoyfreq_smooth_vert" TYPE="logical">
!  To smooth the vertical temp and salt derivative for diagnosing 
!  the buoyancy frequency. Default buoyfreq_smooth_vert=.true.
!  </DATA>
!
!  <DATA NAME="epsln_drhodz" UNITS="kg/m4" TYPE="real">
!  To normalize the inverse vertical derivative of
!  locally referenced potential density for computing the
!  buoyancy frequency. Default epsln_drhodz=1e-10.
!  </DATA>
!
!  <DATA NAME="epsln_drhodz_diag" UNITS="kg/m4" TYPE="real">
!  To normalize the inverse vertical derivative of locally referenced
!  potential density for computing neutral_rho and wdian diagnostics. 
!  Default epsln_drhodz_diag=1e-10.
!  </DATA>
!
!  <DATA NAME="smax_diag" UNITS="dimensionless" TYPE="real">
!  A diagnostic maximum neutral slope for use in computing which direction 
!  is deemed the most stratified.  For use in computing the stratification_factor
!  which is then used to diagnose the dianeutral mass transport.  
!  smax_diag should corresond to the choice used in neutral diffusion scheme.  
!  Should have 0 <= smax_diag <= 1.0. 
!  Default smax_diag=-1.0, in which case we compute the smax according to 
!  the vertical to horizontal grid aspect ratio.  This method ensures that 
!  the slope is adequately "resolved" by the grid.
!  </DATA>
!  <DATA NAME="smax_min_in_column" TYPE="logical">
!  To compute the diagnostic maximum neutral slope within a column as the minimum 
!  vertical to horizontal grid aspect ratio.  This method ensures that 
!  the slope is adequately "resolved" by the grid, and that all depths use the 
!  same definition of "resolved", even if presumably thicker grid cells can 
!  "resolve" larger neutral slopes.  This approach is not very useful generally,
!  so it is retained only for testing purposes.  
!  Default smax_min_in_column=.false.  
!  </DATA>
!
!  <DATA NAME="mask_domain_restart" TYPE="logical">
!  For cases where use the domain masking, it is necessary to initialize the field 
!  denominator_r to nonzero in order to avoid NaNs in the case when change processor
!  layout in between restarts.  Note that when use solid wall boundary conditions, 
!  this logical should remain false in order to bitwise reproduce across restarts.
!  Default mask_domain_restart=.false. 
!  </DATA>
!
!  <DATA NAME="drhodz_diag_stable" TYPE="logical">
!  When computing drhodz_diag, we can enforce that it is negative,
!  thus reflecting a stable stratification.  The field drhodz_diag 
!  is used for many water mass transformation diagnostics, such as 
!  wdian_rho.  Allowing for unstable profiles can bias the wdian_rho 
!  calculation in an improper way, since the magnitude of drhodz_diag
!  is very small when it is positive, whereas it is larger magnitude when
!  negative.  Default drhodz_diag_stable=.true. 
!  </DATA>
!
!  <DATA NAME="grad_nrho_lrpotrho_compute" TYPE="logical">
!  To perform the diagnostic calculation of grad_nrho_lrpotrho
!  for analysis diagnostics.  This factor is not well constrained,
!  and can be problematic in certain regions.  So presently we do 
!  not recommend computing it, so that the default 
!  is grad_nrho_lrpotrho_compute=.false. 
!  </DATA>
!  <DATA NAME="grad_nrho_lrpotrho_max" UNITS="dimensionless" TYPE="real">
!  Maximum value used for grad_nrho_lrpotrho.  
!  Default grad_nrho_lrpotrho_max=10.
!  </DATA>
!  <DATA NAME="grad_nrho_lrpotrho_min" UNITS="dimensionless" TYPE="real">
!  Minimum value used for grad_nrho_lrpotrho.  
!  Default grad_nrho_lrpotrho_min=1.
!  </DATA>
!
!  <DATA NAME="smooth_stratification_factor" TYPE="logical">
!  For doing an S2D smoothing of the stratification factor used 
!  for diagnostic purposes.  Requires an extra call to mpp update.
!  Default smooth_stratification_factor=.false. since the smoothing 
!  incurs a cost that should be borne only when desired.  
!  </DATA>
!
!  <DATA NAME="update_diagnostic_factors" TYPE="logical">
!  To update the watermass_factor and stratification_factor
!  for use in the water mass transformation diagnostics.  
!  Default update_diagnostic_factors=.false. 
!  </DATA>
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging nonlinear equation of state 
!  </DATA>
!
!  <DATA NAME="rho0_density" TYPE="logical">
!  For debugging, it is often useful to have rho=rho0 uniform.
!  </DATA>
!
!  <DATA NAME="density_equal_potrho" TYPE="logical">
!  For idealized tests, set the in situ density equal to the 
!  potential density referenced to potrho_press.  All density 
!  derivatives will also be computed with respect to constant
!  potrho_press pressure.  
!  Default density_equal_potrho=.false.
!  </DATA>
!
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase 
!  efficiency.
!  default: do_bitwise_exact_sum=.false.
!  </DATA>
!
!</NAMELIST>

#include <fms_platform.h>

use constants_mod,       only: epsln, c2dbars
use diag_manager_mod,    only: register_diag_field, register_static_field, diag_axis_init
use diag_manager_mod,    only: need_data, send_data
use fms_mod,             only: write_version_number, mpp_error
use fms_mod,             only: field_exist, FATAL, WARNING
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_io_mod,          only: register_restart_field, save_restart, restore_state
use fms_io_mod,          only: reset_field_pointer, restart_file_type
use mpp_domains_mod,     only: mpp_update_domains, mpp_global_sum
use mpp_domains_mod,     only: BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,             only: input_nml_file, stdout, stdlog
use platform_mod,        only: i8_kind
use time_manager_mod,    only: time_type, increment_time
use field_manager_mod,   only: fm_get_index

use ocean_domains_mod,    only: get_local_indices
use ocean_operators_mod,  only: FDX_T, FDY_T, FMX, FMY, S2D 
use ocean_parameters_mod, only: GEOPOTENTIAL, ZSTAR, ZSIGMA, PRESSURE, PSTAR, PSIGMA 
use ocean_parameters_mod, only: CONSERVATIVE_TEMP, POTENTIAL_TEMP, PRACTICAL_SALT, PREFORMED_SALT
use ocean_parameters_mod, only: missing_value, onefourth, onehalf, rho0r, rho0, grav
use ocean_pressure_mod,   only: pressure_in_dbars, pressure_in_dbars_blob
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_density_type, ocean_external_mode_type
use ocean_types_mod,      only: ocean_diag_tracer_type, ocean_lagrangian_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d, diagnose_3d, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk6
use ocean_workspace_mod,  only: wrk1_v, wrk2_v, wrk3_v, wrk4_v
use ocean_workspace_mod,  only: wrk1_2d, wrk2_2d, wrk3_2d, wrk4_2d

implicit none

private

real, dimension(:), allocatable     :: a, b                   ! polynomial coefficients in the realistic EOS
real, dimension(:,:), allocatable   :: rhodz_tau              ! vertical integral of rho*dzt at time tau
real, dimension(:,:), allocatable   :: rhodz_taup1            ! vertical integral of rho*dzt at time taup1
real, dimension(:,:,:), allocatable :: denominator_r          ! reciprocal of denominator for relastic EOS

integer  :: vert_coordinate        ! for setting which vertical coordinate to use
integer  :: global_sum_flag        ! for determing type of global sum BITWISE/NON

! polynomial coefficients in the preTEOS10 EOS
real :: a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11      
real :: b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12

! multiplied constants for density partial derivs in the preTEOS10 EOS
real :: two_a2, three_a3, two_a6, six_a3, two_a8, two_a10, two_a11, four_a11             
real :: two_b2, three_b3, six_b3, four_b4, twelve_b4, three_b7, six_b7             
real :: onep5_b8, onep5_b9, two_b9, three_b9, two_b11, three_b11, six_b11, three_b12   

! polynomial coefficients in the TEOS10 eos 
real :: v01,v02,v03,v04,v05,v06,v07,v08,v09,v10
real :: v11,v12,v13,v14,v15,v16,v17,v18,v19,v20
real :: v21,v22,v23,v24,v25,v26,v27,v28,v29,v30
real :: v31,v32,v33,v34,v35,v36,v37,v38,v39,v40
real :: v41,v42,v43,v44,v45,v46,v47,v48

! multiplied constants for density partial derivs in the TEOS10 EOS
real :: two_v03, two_v07, two_v10, two_v14, two_v17, two_v19, two_v18, two_v20, two_v23, two_v28, two_v33
real :: two_v36, two_v39, two_v43, two_v44, two_v45, two_v46
real :: three_v04, three_v11, three_v24, three_v29, three_v34, three_v40, three_v47, three_v48
real :: four_v25, four_v30, four_v35
real :: onep5_v08, onep5_v09, onep5_v10, onep5_v11, onep5_v31, onep5_v32, onep5_v33, onep5_v34, onep5_v35
real :: six_v04, six_v11, six_v24, six_v29, six_v34, six_v40, six_v47, six_v48
real :: twelve_v25, twelve_v30, twelve_v35
real :: p75_v08, p75_v09, p75_v10, p75_v11, p75_v31, p75_v32, p75_v33, p75_v34, p75_v35
real :: three_v10, three_v33
real :: fourp5_v11, fourp5_v34
real :: four_v19, four_v45
real :: six_v35

! polynomial coefficients in neutral density 
real :: a0n, a1n, a2n, a3n, a4n, a5n, a6n
real :: b0n, b1n, b2n, b3n, b4n, b5n, b6n, b7n, b8n, b9n

! time steps 
real :: dtts=0.0

! inverse area (1/m) of k=1 tracer cells 
real :: area_total_r  

! inverse neutralrho_interval
real :: neutralrho_interval_r

! inverse gravity constant 
real :: grav_r 

! for diagnostics manager 
logical :: used
integer :: id_eos_salinity   =-1
integer :: id_drhodtheta     =-1
integer :: id_drhodsalt      =-1
integer :: id_dpotrhodtheta  =-1
integer :: id_dpotrhodsalt   =-1
integer :: id_dpotrhodpress  =-1
integer :: id_drhodpress     =-1
integer :: id_thermal_expand =-1
integer :: id_haline_contract=-1
integer :: id_sound_speed2   =-1
integer :: id_press          =-1
integer :: id_rho            =-1
integer :: id_rhoE           =-1
integer :: id_rho0           =-1
integer :: id_neutral_rho    =-1
integer :: id_pot_rho        =-1
integer :: id_pot_rho_0      =-1
integer :: id_pot_rho_1      =-1
integer :: id_pot_rho_2      =-1
integer :: id_pot_rho_3      =-1
integer :: id_pot_rho_4      =-1
integer :: id_pot_rho_et     =-1
integer :: id_pot_rho_nt     =-1
integer :: id_pot_rho_wt     =-1
integer :: id_rho_average    =-1
integer :: id_mass_level     =-1
integer :: id_volume_level   =-1
integer :: id_drhodx_zt      =-1  
integer :: id_drhody_zt      =-1  
integer :: id_drhodz_zt      =-1  
integer :: id_drhodz_wt      =-1  
integer :: id_drhodz_diag    =-1  
integer :: id_buoyfreq2_zt   =-1  
integer :: id_buoyfreq2_wt   =-1  
integer :: id_cabbeling      =-1
integer :: id_thermobaricity =-1
integer :: id_int_rhodz      =-1

integer :: id_rhoave              =-1
integer :: id_pbot_adjust         =-1
integer :: id_eta_adjust          =-1
integer :: id_eta_adjust_approx   =-1
integer :: id_eta_adjust_cstvolume=-1

integer :: id_grad_nrho             =-1
integer :: id_grad_lrpotrho         =-1
integer :: id_grad_nrho_lrpotrho    =-1
integer :: id_watermass_factor      =-1
integer :: id_stratification_factor =-1
integer :: id_stratification_axis   =-1
integer :: id_smax_dianeutral       =-1

!for restart
integer                       :: id_restart_rho = 0
integer                       :: id_restart_rho_s = 0
type(restart_file_type), save :: Den_restart

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type),   pointer :: Grd =>NULL()

character(len=128) :: version=&
       '$Id: ocean_density.F90,v 20.0 2013/12/14 00:10:40 fms Exp $'
character (len=128) :: tagname = &
       '$Name: tikal $'

public ocean_density_init
public ocean_density_end
public ocean_density_diag
public update_ocean_density
public update_ocean_density_salinity
public density
public potential_density
public neutral_density
public density_sfc
public density_derivs
public density_delta_z
public density_delta_sfc
public ocean_density_restart
public calc_cabbeling_thermobaricity
public buoyfreq2
public compute_buoyfreq

private compute_drhodxy
private ocean_density_chksum
private cabbeling_thermobaricity
private compute_diagnostic_factors
private compute_density_diagnostics 
private density_coeffs_init
private density_diagnostics_init

! interfaces for density 
interface density
   module procedure density_field
   module procedure density_level
   module procedure density_line
   module procedure density_point
end interface

! interfaces for neutral density 
interface neutral_density
   module procedure neutral_density_field
   module procedure neutral_density_point
end interface


! interfaces for density derivatives 
interface density_derivs
   module procedure density_derivs_field
   module procedure density_derivs_level
   module procedure density_derivs_point
end interface

integer :: index_temp   =-1
integer :: index_salt   =-1
integer :: index_Fdelta =-1
integer :: num_prog_tracers
integer :: salt_variable
integer :: temp_variable

! indices for diagnostic tracer density elements
integer :: ind_rho        = -1
integer :: ind_neutralrho = -1
integer :: ind_potrho_0   = -1
integer :: ind_potrho_1   = -1
integer :: ind_potrho_2   = -1
integer :: ind_potrho_3   = -1
integer :: ind_potrho_4   = -1

! nml settings 

! for linear EOS 
logical :: eos_linear=.false.    
real    :: alpha_linear_eos=0.255
real    :: beta_linear_eos =0.0

! for teos10 or preteos10 eos
logical :: eos_preteos10 = .false.
logical :: eos_teos10    = .false.

! some test settings for density 
real :: s_test=20.0
real :: t_test=20.0
real :: p_test=1000.0
real :: jmfwg_rho 
real :: jmfwg_alpha 
real :: jmfwg_beta 
real :: mbfj_rho 
real :: mbfj_alpha 
real :: mbfj_beta 

! some test settings for neutral density 
real :: sn_test=35.0
real :: tn_test=20.0
real :: mj_neutralrho 
real :: mb_neutralrho 

! for diagnostic partitioning of vertical according to 
! potential density, neutral density, or theta classes 
integer :: layer_nk = 80 ! # of classes used to compute diagnostics associated with 
                         ! potential density, potential temperature, or conservative
                         ! temperature classes. 

! for diagnostic partitioning of vertical according to potential density classes 
real :: potrho_press = 2000.0  ! sea pressure (dbars) for potential density computed for diagnostics
real :: potrho_min   = 1028.0  ! (kg/m^3)         
real :: potrho_max   = 1038.0  ! (kg/m^3)           

! for diagnostic partitioning of vertical according to 
! rational polynomial approximation to neutral density
logical :: neutral_density_omega  = .false.
logical :: neutral_density_potrho = .true.
logical :: neutral_density_theta  = .false.
real    :: neutralrho_min = 1020.0  ! (kg/m^3)         
real    :: neutralrho_max = 1030.0  ! (kg/m^3)           

! for diagnostic partitioning of vertical according 
! to potential temperature or conservative temperature classes 
real :: theta_min = -2.0  ! (degrees C)         
real :: theta_max = 30.0  ! (degrees C)           

! standard atmospheric pressure (dbar).
! Should be set to 0.0 if assume zero pressure for the 
! overlying atmosphere.  But if have a realistic atmospheric
! pressure loading, then set press_standard=10.1325.  
real :: press_standard = 0.0 

! for smoothing the vertical density gradient 
! when diagnosing the buoyancy frequency. 
logical :: buoyfreq_smooth_vert=.true.
integer :: num_121_passes      =1
real    :: epsln_drhodz        =1.e-10

! for enforcing negative drhodz_diag
logical :: drhodz_diag_stable = .true.
real    :: epsln_drhodz_diag = 1.e-10

! for setting grad_nrho_lrpotrho=1.0 
logical :: grad_nrho_lrpotrho_compute=.false.
real    :: grad_nrho_lrpotrho_max=10.0
real    :: grad_nrho_lrpotrho_min=1.0

! for stratification factor used for diagnostics 
logical :: smooth_stratification_factor=.false.
logical :: smax_min_in_column=.false.
real    :: smax_diag=-1.0
real, dimension(:,:,:), allocatable :: smax
real, dimension(:,:),   allocatable :: dhorz_r

! to compute diagnostic factors for 
! water mass transformation diagnostics 
logical :: update_diagnostic_factors=.false. 

! for case when use domain masking, 
logical :: mask_domain_restart = .false. 

logical :: debug_this_module     = .false. 
logical :: rho0_density          = .false.
logical :: density_equal_potrho  = .false. 
logical :: module_is_initialized = .false.
logical :: write_a_restart       = .true. 
logical :: do_bitwise_exact_sum  = .false.

namelist /ocean_density_nml/ s_test, t_test, p_test, press_standard,                  &
                             sn_test, tn_test,                                        &
                             eos_linear, alpha_linear_eos, beta_linear_eos,           & 
                             eos_preteos10, eos_teos10,                               &
                             potrho_press, potrho_min, potrho_max,                    &
                             neutralrho_min, neutralrho_max,                          &
                             layer_nk, theta_min, theta_max,                          &
                             debug_this_module, write_a_restart,                      &
                             rho0_density, density_equal_potrho,                      &
                             buoyfreq_smooth_vert, num_121_passes, epsln_drhodz,      &
                             mask_domain_restart, do_bitwise_exact_sum,               &
                             drhodz_diag_stable, epsln_drhodz_diag,                   &
                             grad_nrho_lrpotrho_compute, grad_nrho_lrpotrho_max,      &
                             grad_nrho_lrpotrho_min,                                  &
                             neutral_density_omega, neutral_density_potrho,           &
                             neutral_density_theta,                                   &
                             smooth_stratification_factor, update_diagnostic_factors, &
                             smax_diag, smax_min_in_column 
contains

!#######################################################################
! <SUBROUTINE NAME="ocean_density_init">
!
! <DESCRIPTION>
! Initialize the density module
! </DESCRIPTION>
!
  subroutine ocean_density_init (Grid, Domain, Time, Time_steps, Thickness, T_prog, T_diag, &
                                 Ocean_options, Dens, ver_coordinate, use_blobs, debug)

    type(ocean_grid_type),        target,  intent(in)    :: Grid
    type(ocean_domain_type),      target,  intent(in)    :: Domain
    type(ocean_time_type),                 intent(in)    :: Time
    type(ocean_time_steps_type),           intent(in)    :: Time_steps
    type(ocean_thickness_type),            intent(in)    :: Thickness
    type(ocean_prog_tracer_type), target,  intent(in)    :: T_prog(:)
    type(ocean_diag_tracer_type),          intent(inout) :: T_diag(:)
    type(ocean_options_type),              intent(inout) :: Ocean_options
    type(ocean_density_type),     target,  intent(inout) :: Dens

    integer,   intent(in)            :: ver_coordinate
    logical,   intent(in)            :: use_blobs
    logical,   intent(in), optional  :: debug
    
    integer :: ioun, io_status, ierr
    integer :: i, j, k, n
    integer :: tau, taup1, taum1
    integer :: num_eos
    real    :: rho0_nk(nk)

    real, allocatable, dimension(:) :: rho_average
    real, allocatable, dimension(:) :: mass_level
    real, allocatable, dimension(:) :: volume_level 

    integer :: id_grid_xt, id_grid_yt
    integer :: id_grid_xu, id_grid_yu

    integer :: id_potrho_bounds,     id_potrho_axis
    integer :: id_neutralrho_bounds, id_neutralrho_axis
    integer :: id_theta_bounds,      id_theta_axis

    integer :: id_restart(5)
    integer :: neutral_rho_method=0

    real, allocatable, dimension(:) :: potrho_bounds
    real, allocatable, dimension(:) :: neutralrho_bounds
    real, allocatable, dimension(:) :: theta_bounds

    ! salinity and temperature pointers to aid readability 
    real, dimension(:,:,:), pointer :: salinity    => NULL()
    real, dimension(:,:,:), pointer :: temperature => NULL()

    real    :: potrho_interval 
    real    :: neutralrho_interval 
    real    :: theta_interval 
    character*128 filename

    integer :: stdoutunit,stdlogunit 
    stdoutunit=stdout();stdlogunit=stdlog() 


    if ( module_is_initialized ) then 
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod (ocean_density_init): module already initialized.')
    endif 

    module_is_initialized = .TRUE.

    tau             = Time%tau
    taup1           = Time%taup1
    taum1           = Time%taum1
    vert_coordinate = ver_coordinate
 
    num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'temp')   index_temp   = n
       if (T_prog(n)%name == 'salt')   index_salt   = n
       if (T_prog(n)%name == 'fdelta') index_Fdelta = n
       if (T_prog(n)%longname == 'Conservative temperature') temp_variable = CONSERVATIVE_TEMP
       if (T_prog(n)%longname == 'Potential temperature')    temp_variable = POTENTIAL_TEMP
       if (T_prog(n)%longname == 'Practical Salinity')       salt_variable = PRACTICAL_SALT
       if (T_prog(n)%longname == 'Preformed Salinity')       salt_variable = PREFORMED_SALT
    enddo

    call write_version_number(version, tagname)

    ! provide for namelist override of defaults
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_density_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_density_nml')
#else
    ioun = open_namelist_file()
    read (ioun,ocean_density_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status, 'ocean_density_nml')
    call close_file(ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_density_nml)
    write (stdlogunit,ocean_density_nml)

    num_eos = 0
    if(eos_linear) then 
      num_eos = num_eos+1
    elseif(eos_preteos10) then 
      num_eos = num_eos+1
    elseif(eos_teos10) then 
      num_eos = num_eos+1
    endif        
    if(num_eos == 0) then
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: No EOS chosen. Set one to true: eos_linear, eos_preteos10, eos_teso10.')
    endif 
    if(num_eos > 1) then
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: More than one EOS chosen. Set just one true: eos_linear, eos_preteos10, eos_teso10.')
    endif 
    Dens%use_teos10 = eos_teos10


    if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
    endif 
    if(debug_this_module) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_density with debug_this_module=.true.'  
    endif 

    if(.not. write_a_restart) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_density with write_a_restart=.false.'
      write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
    endif 

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

    if(eos_linear) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: USING linear EOS designed for idealized simulations.'
        write (stdoutunit,'(7x,a)') &
        'It is a linear function of potential temperature and salinity.'
        write (stdoutunit,'(7x,a)') &
        'There is no pressure dependence.'
        Ocean_options%equation_of_state = 'Linear equation of state used for density.'
    elseif(eos_preteos10) then 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: USING preTEOS10 EOS, as relevant for realistic ocean climate simulations.' 
        write (stdoutunit,'(1x,a,f12.6,a)') &
        ' Subtracting standard atmosphere of ',press_standard,' dbar for EOS calculation.'  
        Ocean_options%equation_of_state = 'preTEOS10 equation of state used for density.'
    elseif(eos_teos10) then 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: USING TEOS10 EOS, as relevant for realistic ocean climate simulations.' 
        write (stdoutunit,'(1x,a,f12.6,a)') &
        ' Subtracting standard atmosphere of ',press_standard,' dbar for EOS calculation.'  
        Ocean_options%equation_of_state = 'TEOS10 equation of state used for density.'
    endif

    if(temp_variable==CONSERVATIVE_TEMP) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing EOS assuming prognostic temp = conservative temperature.'
    elseif(temp_variable==POTENTIAL_TEMP) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing EOS assuming prognostic temp = potential temperature.'
    else
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: model temperature variable remains unspecified')
    endif 

    if(salt_variable==PRACTICAL_SALT) then
        if(eos_teos10) then 
            call mpp_error(FATAL, &
            '==>Error in ocean_density_mod: eos_teos10 requires preformed salinity as prognostic salinity.')
        endif
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing EOS assuming prognostic salinity = practical salinity.'
    elseif(salt_variable==PREFORMED_SALT) then
        if(eos_linear .or. eos_preteos10) then 
            call mpp_error(FATAL, &
            '==>Error in ocean_density_mod: preformed salinity as prognostic salinity requires eos_teos10.')
        endif
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing TEOS10 EOS assuming prognostic salinity = preformed salinity.'
    else
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: model salinity variable remains unspecified')
    endif 

    write(stdoutunit,*) &
    '==>Note: The Boussinesq rho0 density has a value of (kg/m^3) ',rho0
    if(rho0_density) then 
        write(stdoutunit,*) &
        '==>Warning: rho0_density=.true, so rho=rho0 everywhere. Is this really what you wish to do?'
    endif

    if(neutral_density_omega) then 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing diagnostic neutral_rho according to Klocker/McDougall (not yet coded).'
        neutral_rho_method = neutral_rho_method + 1 
        call mpp_error(WARNING, &
       '==>Warning in ocean_density_mod: neutral_density_omega not coded. Will use potential density 2000dbar.')
    endif 
    if(neutral_density_potrho) then 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing diagnostic neutral_rho as potential density referenced to pressure potrho_press.'
        neutral_rho_method = neutral_rho_method + 1 
    endif 
    if(neutral_density_theta) then 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing diagnostic neutral_rho as conservative temperature.'
        neutral_rho_method = neutral_rho_method + 1 
    endif
    if(neutral_rho_method==0) then 
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: set one of the options "neutral_density_omega" or "neutral_density_potrho" or "neutral_density_theta" true.')
    endif   
    if(neutral_rho_method > 1) then 
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: set only **one** of the options "neutral_density_omega" or "neutral_density_potrho" or "neutral_density_theta" true.')
    endif   

    if(drhodz_diag_stable) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Enforcing drhodz_diag < 0, so to use stable stratification for certain diagnostic purposes.'
    endif 

    if(update_diagnostic_factors) then 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing diagnostic factors in ocean_density.F90 for watermass diagnostics.'
    else 
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Diagnostic factors are NOT computed. So if enable watermass diagnostics, they will be corrupted.'
    endif 

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
    allocate(Dens%rho(isd:ied,jsd:jed,nk,3))
    allocate(Dens%rho_fresh(isd:ied,jsd:jed,nk))
    allocate(Dens%rho_salinity(isd:ied,jsd:jed,nk,3))
    allocate(Dens%rho_dztr_tau(isd:ied,jsd:jed,nk))
    allocate(Dens%potrho(isd:ied,jsd:jed,nk))
    allocate(Dens%neutralrho(isd:ied,jsd:jed,nk))
    allocate(Dens%pressure_at_depth(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodT(isd:ied,jsd:jed,nk))
    allocate(Dens%dpotrhodT(isd:ied,jsd:jed,nk))
    allocate(Dens%dpotrhodS(isd:ied,jsd:jed,nk))
    allocate(Dens%dpotrhodP(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodS(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodP(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodz_wt(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodz_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodx_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%drhody_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodz_diag(isd:ied,jsd:jed,nk))
    allocate(Dens%watermass_factor(isd:ied,jsd:jed,nk))
    allocate(Dens%stratification_factor(isd:ied,jsd:jed,nk))
    allocate(Dens%dTdz_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%dSdz_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%mld_subduction(isd:ied,jsd:jed))
    if (use_blobs) allocate(Dens%rhoT(isd:ied,jsd:jed,nk))
#endif

    do k=1,nk 
       Dens%rho_dztr_tau(:,:,k)      = rho0r*Grid%tmask(:,:,k)
       Dens%rho_fresh(:,:,k)         = 1000.0*Grid%tmask(:,:,k)
       Dens%potrho(:,:,k)            = rho0*Grid%tmask(:,:,k)
       Dens%neutralrho(:,:,k)        = rho0*Grid%tmask(:,:,k)
       Dens%pressure_at_depth(:,:,k) = rho0*grav*Thickness%depth_zt(:,:,k)*c2dbars
    enddo
    Dens%rho(:,:,:,1)            = Grid%tmask(:,:,:)*rho0 
    Dens%rho(:,:,:,2)            = Grid%tmask(:,:,:)*rho0 
    Dens%rho(:,:,:,3)            = Grid%tmask(:,:,:)*rho0 
    Dens%rho_salinity(:,:,:,1)   = Grid%tmask(:,:,:)*35 
    Dens%rho_salinity(:,:,:,2)   = Grid%tmask(:,:,:)*35 
    Dens%rho_salinity(:,:,:,3)   = Grid%tmask(:,:,:)*35 
    Dens%drhodT(:,:,:)           = 0.0
    Dens%dpotrhodT(:,:,:)        = 0.0
    Dens%dpotrhodS(:,:,:)        = 0.0
    Dens%dpotrhodP(:,:,:)        = 0.0
    Dens%drhodS(:,:,:)           = 0.0
    Dens%drhodP(:,:,:)           = 0.0
    Dens%drhodz_wt(:,:,:)        = Grid%tmask(:,:,:)*epsln_drhodz
    Dens%drhodz_zt(:,:,:)        = Grid%tmask(:,:,:)*epsln_drhodz
    Dens%drhodx_zt(:,:,:)        = 0.0
    Dens%drhody_zt(:,:,:)        = 0.0
    Dens%drhodz_diag(:,:,:)      = Grid%tmask(:,:,:)*epsln_drhodz_diag
    Dens%dTdz_zt(:,:,:)          = 0.0
    Dens%dSdz_zt(:,:,:)          = 0.0
    Dens%watermass_factor(:,:,:) = Grid%tmask(:,:,:)
    Dens%mld_subduction(:,:)     = 0.0
    do k=1,nk   
       Dens%stratification_factor(:,:,k) = Grid%tmask(:,:,k)*Grid%dzt(k)*Grid%dat(:,:)
    enddo
    if (use_blobs) Dens%rhoT(:,:,:) = 0.0

    Dom  => Domain
    Grd  => Grid
    dtts = Time_steps%dtts

    ! for diagnostics    
    area_total_r = mpp_global_sum(Dom%domain2d,Grd%dat(:,:)*Grd%tmask(:,:,1), global_sum_flag)
    if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) area_total_r = real(area_total_r,kind=FLOAT_KIND)
    area_total_r = 1.0/(epsln+area_total_r) 

    ! for diagnostic calculation of dianeutral transport 
    allocate(smax(isd:ied,jsd:jed,nk))
    allocate(dhorz_r(isd:ied,jsd:jed))
    dhorz_r(:,:) = 0.0
    smax(:,:,:)  = 0.0
    do j=jsd,jed
       do i=isd,ied
          dhorz_r(i,j) = Grd%tmask(i,j,1)  &
          *( Grd%dxt(i,j) + Grd%dyt(i,j) ) / (epsln + 2.0*Grd%dxt(i,j)*Grd%dyt(i,j) ) 
       enddo
    enddo

    grav_r = 1.0/grav
    allocate(rhodz_tau(isd:ied,jsd:jed))
    allocate(rhodz_taup1(isd:ied,jsd:jed))
    rhodz_tau   = 0.0
    rhodz_taup1 = 0.0


   ! send as 1-dim field, due to problems with diag manager 
   id_rho0 = register_static_field('ocean_model','rho0', Grid%tracer_axes(3:3), &
             'reference density for Boussinesq approximation','kg/m^3',         &
             missing_value=missing_value, range=(/-1e0,1e10/),                  &
             standard_name='reference_sea_water_density_for_boussinesq_approximation')
   rho0_nk(:) = rho0 
   if (id_rho0 > 0) then 
      used = send_data(id_rho0, rho0_nk(:), Time%model_time)
   endif 

   allocate(denominator_r(isd:ied,jsd:jed,nk))
   denominator_r(:,:,:) = 0.0

   ! specify coefficients for the polynomical equations of state 
   ! and perform some pointwise tests.  
   call density_coeffs_init()

   filename = 'ocean_density.res.nc'    
   id_restart_rho = register_restart_field(Den_restart, filename, 'rho', Dens%rho(:,:,:,taup1), &
                   domain=Dom%domain2d )
   id_restart_rho_s = register_restart_field(Den_restart, filename, 'rho_salinity', Dens%rho_salinity(:,:,:,taup1), &
                   domain=Dom%domain2d )
   id_restart(1)  = register_restart_field(Den_restart, filename, 'pressure_at_depth', Dens%pressure_at_depth(:,:,:), &
                   domain=Dom%domain2d )
   id_restart(2)  = register_restart_field(Den_restart, filename, 'denominator_r', denominator_r(:,:,:), &
                   domain=Dom%domain2d )
   id_restart(3)  = register_restart_field(Den_restart, filename, 'drhodT', Dens%drhodT(:,:,:), &
                   domain=Dom%domain2d )
   id_restart(4)  = register_restart_field(Den_restart, filename, 'drhodS', Dens%drhodS(:,:,:), &
                   domain=Dom%domain2d )
   id_restart(5)  = register_restart_field(Den_restart, filename, 'drhodz_zt', Dens%drhodz_zt(:,:,:), &
                   domain=Dom%domain2d )

   ! pointers to relevant salinity and temperature to enhance readability    
   salinity    => Dens%rho_salinity(:,:,:,tau)
   temperature => T_prog(index_temp)%field(:,:,:,tau)

   filename = 'INPUT/ocean_density.res.nc'
   if (.NOT.file_exist(trim(filename)) )then

      write (stdoutunit,'(/a)') '  Initialising salinity for use in density calculation'
      call update_ocean_density_salinity(T_prog,taum1,Dens)
      call update_ocean_density_salinity(T_prog,tau,Dens)
      call update_ocean_density_salinity(T_prog,taup1,Dens)

      Dens%rho(:,:,:,taum1) = density(salinity, temperature, Dens%pressure_at_depth(:,:,:))
      Dens%rho(:,:,:,tau)   = Dens%rho(:,:,:,taum1)
      Dens%rho(:,:,:,taup1) = Dens%rho(:,:,:,taum1)
      if (use_blobs) Dens%rhoT(:,:,:) = Dens%rho(:,:,:,taup1)

   else

       write (stdoutunit,'(/a)') '  Reading restart for density from INPUT/ocean_density.res.nc'
       call restore_state( Den_restart, id_restart_rho )
       call mpp_update_domains(Dens%rho(:,:,:,taup1),        Dom%domain2d)
       call restore_state( Den_restart, id_restart(1) )
       call mpp_update_domains(Dens%pressure_at_depth(:,:,:),Dom%domain2d)
       call restore_state( Den_restart, id_restart(2) )

       ! initialize to epsln to eliminate NaNs when mask 
       ! out processors and then change layout upon a restart 
       if(mask_domain_restart) then 
          where (denominator_r==0.)  denominator_r=rho0r
       endif 
       call mpp_update_domains(denominator_r(:,:,:),Dom%domain2d)

       ! determine whether to read density derivs at the initial time step of integration.
       ! early versions of MOM4p1 did not contain drhodT and drhodS in the restart file, so 
       ! to allow older restart files to be used with newer releases of MOM4p1, we 
       ! provide for the possibility that the restarts do not contain drhodT and drhodS.
       if(field_exist(filename, 'drhodT')) then
           call restore_state( Den_restart, id_restart(3) )
           call restore_state( Den_restart, id_restart(4) )
           call mpp_update_domains(Dens%drhodT(:,:,:),Dom%domain2d)
           call mpp_update_domains(Dens%drhodS(:,:,:),Dom%domain2d)
       else 
           write(stdoutunit,'(a)') '==>Note: ocean_density_mod: did not read density derivatives from restart.' 
       endif 

       ! determine whether to read rho_salinity at the initial time step of
       ! integration.  Early versions of MOM4p1 did not contain the rho_salinity
       ! in the restart file.  This allows older restart files to be used with
       ! newer releases of MOM.
       if(field_exist(filename, 'rho_salinity')) then
          write (stdoutunit,'(/a)') '  Initializing salinity for use in density calculation from restart'
          call restore_state( Den_restart, id_restart_rho_s )
          Dens%rho_salinity(:,:,:,taum1) = Dens%rho_salinity(:,:,:,taup1)
          Dens%rho_salinity(:,:,:,tau) = Dens%rho_salinity(:,:,:,taup1)
          call mpp_update_domains(Dens%rho_salinity(:,:,:,taum1),Dom%domain2d)
          call mpp_update_domains(Dens%rho_salinity(:,:,:,tau),  Dom%domain2d)
          call mpp_update_domains(Dens%rho_salinity(:,:,:,taup1),Dom%domain2d)
       else
          write (stdoutunit,'(/a)') '  Initialising salinity for use in density calculation'
          call update_ocean_density_salinity(T_prog,taum1,Dens)
          call update_ocean_density_salinity(T_prog,tau,Dens)
          call update_ocean_density_salinity(T_prog,taup1,Dens)
       endif

   endif

   ! compute buoyancy frequency for diagnostic purposes 
   ! note: Dens%rhoT and the buoyancy frequency is initialised in ocean_blob_init
   ! when use_blobs=.true.
   if (.not. use_blobs) then
      call compute_buoyfreq(Time, Thickness, salinity, temperature, Dens, use_blobs) 
   endif

   ! over-write drhodz_zt from compute_buoyfreq with drhodz_zt from restart. 
   ! Note that we determine whether to read drhodz_zt at the initial time 
   ! step of integration. early versions of MOM4p1 did not contain drhodz_zt 
   ! in the restart file, so to allow older restart files to be used with 
   ! newer releases of MOM4p1, we provide for the possibility that the 
   ! restarts do not contain drhodz_zt.
   if (file_exist(trim(filename)) )then 
       if(field_exist(filename, 'drhodz_zt')) then
           call restore_state( Den_restart, id_restart(5) )
          call mpp_update_domains(Dens%drhodz_zt(:,:,:),Dom%domain2d)
       else 
          write(stdoutunit,'(a)') '==>Note: ocean_density_mod: did not read drhodz_zt from restart.' 
       endif
   endif

   ! compute neutral density for diagnostic purposes 
   Dens%neutralrho(:,:,:) = Grd%tmask(:,:,:)*neutral_density(salinity, temperature)
   ind_neutralrho = fm_get_index('/ocean_mod/diag_tracers/rho_neutral')
   if (ind_neutralrho .gt. 0) then
     T_diag(ind_neutralrho)%field(:,:,:) = Dens%neutralrho(:,:,:)
   endif

   ! compute some density related diagnostic factors
   ! initialize neutralrho_interval_r first
   neutralrho_interval  = (neutralrho_max-neutralrho_min)/(epsln+layer_nk)
   neutralrho_interval_r = 1.0/neutralrho_interval
   call compute_diagnostic_factors(Time, Thickness, Dens, salinity, temperature)

   ! compute potential density for diagnostic purposes 
   Dens%potrho(:,:,:) = Grd%tmask(:,:,:)*potential_density(salinity, temperature, potrho_press)

   ! compute various potential densities for diagnostic purposes 
   ind_potrho_0 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_0')
   if (ind_potrho_0 .gt. 0) then
     T_diag(ind_potrho_0)%field(:,:,:) = Grd%tmask(:,:,:)*potential_density(salinity, temperature, 0.0)
   endif
   ind_potrho_1 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_1')
   if (ind_potrho_1 .gt. 0) then
     T_diag(ind_potrho_1)%field(:,:,:) = Grd%tmask(:,:,:)*potential_density(salinity, temperature, 1000.0)
   endif
   ind_potrho_2 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_2')
   if (ind_potrho_2 .gt. 0) then
     T_diag(ind_potrho_2)%field(:,:,:) = Grd%tmask(:,:,:)*potential_density(salinity, temperature, 2000.0)
   endif
   ind_potrho_3 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_3')
   if (ind_potrho_3 .gt. 0) then
     T_diag(ind_potrho_3)%field(:,:,:) = Grd%tmask(:,:,:)*potential_density(salinity, temperature, 3000.0)
   endif
   ind_potrho_4 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_4')
   if (ind_potrho_4 .gt. 0) then
     T_diag(ind_potrho_4)%field(:,:,:) = Grd%tmask(:,:,:)*potential_density(salinity, temperature, 4000.0)
   endif

   if(rho0_density) then 
       Dens%rho(:,:,:,taum1) = Grd%tmask(:,:,:)*rho0
       Dens%rho(:,:,:,tau)   = Grd%tmask(:,:,:)*rho0
       Dens%rho(:,:,:,taup1) = Grd%tmask(:,:,:)*rho0
   endif
   ind_rho = fm_get_index('/ocean_mod/diag_tracers/rho_in_situ')
   if (ind_rho .gt. 0) then
     T_diag(ind_rho)%field(:,:,:) = Dens%rho(:,:,:,tau)
   endif

   ! compute initial level mass, volume, and density 
   allocate(rho_average(0:nk))
   allocate(mass_level(0:nk))
   allocate(volume_level(0:nk))
   rho_average(:)  = 0.0
   mass_level(:)   = 0.0
   volume_level(:) = 0.0
   wrk1_2d(:,:)    = 0.0
   wrk2_2d(:,:)    = 0.0
   wrk3_2d(:,:)    = 0.0
   do k=1,nk
      wrk1_2d(:,:)    = Grd%dat(:,:)*Grd%tmask(:,:,k)  
      wrk2_2d(:,:)    = wrk1_2d(:,:)*Thickness%dzt(:,:,k)
      wrk3_2d(:,:)    = wrk1_2d(:,:)*Thickness%dzt(:,:,k)*Dens%rho(:,:,k,tau)
      volume_level(k) = mpp_global_sum(Dom%domain2d,wrk2_2d(:,:),global_sum_flag)
      mass_level(k)   = mpp_global_sum(Dom%domain2d,wrk3_2d(:,:),global_sum_flag)
      if (global_sum_flag == NON_BITWISE_EXACT_SUM) then
        volume_level(k) = real(volume_level(k),kind=FLOAT_KIND)
        mass_level(k)   = real(mass_level(k),kind=FLOAT_KIND)
      endif
      if(volume_level(k) > 0.0) rho_average(k) = mass_level(k)/volume_level(k)
   enddo
   do k=1,nk
      volume_level(0) = volume_level(0) + volume_level(k)
      mass_level(0)   = mass_level(0)   + mass_level(k)
   enddo
   if(volume_level(0) > 0.0) rho_average(0) = mass_level(0)/volume_level(0)

   write(stdoutunit,'(/a,e24.12)') &
   ' ==>Note: From ocean_density_mod: Boussinesq reference density rho0(kg/m3) = ',rho0
   write(stdoutunit,'(a,e24.12)')  &
   '          Initial rho_average(kg/m3) = ',rho_average(0)
   if(rho_average(0) /= rho0 .and. Time%init) then 
     write(stdoutunit,'(a)')  &
     '          Since rho0 .ne. rho_average, consider changing rho0 in' 
     write(stdoutunit,'(a/)') &
     '          ocean_parameters.F90 to be equal to rho_average for better accuracy.'
   endif 
 
   id_rho_average = register_static_field('ocean_model','rho_average', Grid%tracer_axes(3:3), &
                   'level average density at initial step of segment','kg/m^3',               &
                   missing_value=missing_value, range=(/-1e0,1e10/))
   id_mass_level  = register_static_field('ocean_model','mass_level', Grid%tracer_axes(3:3), &
                   'mass of levels at initial step of segment','kg',                         &
                   missing_value=missing_value, range=(/-1e0,1e10/))
   id_volume_level= register_static_field('ocean_model','volume_level', Grid%tracer_axes(3:3), &
                   'volume of levels at initial step of segment','m^3',                        &
                   missing_value=missing_value, range=(/-1e0,1e10/))

   if (id_rho_average > 0) then 
      used = send_data(id_rho_average, rho_average(1:nk), Time%model_time)
   endif 
   if (id_mass_level > 0) then 
      used = send_data(id_mass_level, mass_level(1:nk), Time%model_time)
   endif 
   if (id_volume_level > 0) then 
      used = send_data(id_volume_level, volume_level(1:nk), Time%model_time)
   endif 

  ! define vertical axes according to potential density classes  
  allocate ( Dens%potrho_ref(layer_nk))
  allocate ( Dens%potrho_bounds(layer_nk+1))
  allocate ( potrho_bounds(layer_nk+1))
  potrho_bounds(1) = potrho_min
  potrho_interval  = (potrho_max-potrho_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    potrho_bounds(k)=potrho_bounds(k-1)+potrho_interval
  enddo 
  do k=1,layer_nk
    Dens%potrho_ref(k)=potrho_bounds(k)+0.5*potrho_interval
  enddo 
  do k=1,layer_nk+1
    Dens%potrho_bounds(k) = potrho_bounds(k)
  enddo 

  ! define vertical axes according to neutral density classes  
  allocate ( Dens%neutralrho_ref(layer_nk))
  allocate ( Dens%neutralrho_bounds(layer_nk+1))
  allocate ( neutralrho_bounds(layer_nk+1))
  neutralrho_bounds(1) = neutralrho_min
  do k=2,layer_nk+1
    neutralrho_bounds(k)=neutralrho_bounds(k-1)+neutralrho_interval
  enddo 
  do k=1,layer_nk
    Dens%neutralrho_ref(k)=neutralrho_bounds(k)+0.5*neutralrho_interval
  enddo 
  do k=1,layer_nk+1
    Dens%neutralrho_bounds(k) = neutralrho_bounds(k)
  enddo 

  ! define vertical axes according to 
  ! potential temperature or conservative temperature classes
  allocate ( Dens%theta_ref(layer_nk))
  allocate ( Dens%theta_bounds(layer_nk+1))
  allocate ( theta_bounds(layer_nk+1))
  theta_bounds(1)  = theta_min
  theta_interval   = (theta_max-theta_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    theta_bounds(k) =theta_bounds(k-1) +theta_interval
  enddo 
  do k=1,layer_nk
    Dens%theta_ref(k) =theta_bounds(k) +0.5*theta_interval
  enddo 
  do k=1,layer_nk+1
    Dens%theta_bounds(k) = theta_bounds(k)
  enddo 

  ! assuming that if the Grd%name = 'ocean' then we will
  ! be using this diagnostic axis. Otherwise it is not used.
  ! This is done to prevent duplicate axis names for multiple
  ! grid objects.
  if (trim(Grd%name) == 'ocean') then

    id_grid_xt = diag_axis_init ('grid_xt_'//trim(Grd%name),Grd%grid_x_t,'degrees_E','x', &
                                   'tcell longitude',set_name='ocean', Domain2=Dom%domain2d)

    id_grid_yt = diag_axis_init ('grid_yt_'//trim(Grd%name),Grd%grid_y_t,'degrees_N','y', &
                                   'tcell latitude',set_name='ocean', Domain2=Dom%domain2d)

    id_grid_xu = diag_axis_init ('grid_xu_'//trim(Grd%name),Grd%grid_x_u,'degrees_E','x', &
                                   'ucell longitude',set_name='ocean', Domain2=Dom%domain2d)

    id_grid_yu = diag_axis_init ('grid_yu_'//trim(Grd%name),Grd%grid_y_u,'degrees_N','y', &
                                   'ucell latitude',set_name='ocean', Domain2=Dom%domain2d)


    id_potrho_bounds = diag_axis_init                                        &
     ( 'potrho_edges',potrho_bounds, 'kg/m^3', 'z','potential density edges',&
       direction=-1, set_name='ocean')

    id_potrho_axis   = diag_axis_init                           &
     ( 'potrho',Dens%potrho_ref,                                &
       'kg/m^3', 'z','potential density',edges=id_potrho_bounds,&
        direction=-1,set_name='ocean')

    Dens%potrho_axes        = (/ id_grid_xt, id_grid_yt, id_potrho_axis /)
    Dens%potrho_axes_flux_x = (/ id_grid_xu, id_grid_yt, id_potrho_axis /)
    Dens%potrho_axes_flux_y = (/ id_grid_xt, id_grid_yu, id_potrho_axis /)


    id_neutralrho_bounds = diag_axis_init                                          &
     ( 'neutralrho_edges',neutralrho_bounds, 'kg/m^3', 'z','neutral density edges',&
       direction=-1, set_name='ocean')

    id_neutralrho_axis   = diag_axis_init                         &
     ( 'neutral',Dens%neutralrho_ref,                             &
       'kg/m^3', 'z','neutral density',edges=id_neutralrho_bounds,&
        direction=-1,set_name='ocean')

    Dens%neutralrho_axes        = (/ id_grid_xt, id_grid_yt, id_neutralrho_axis /)
    Dens%neutralrho_axes_flux_x = (/ id_grid_xu, id_grid_yt, id_neutralrho_axis /)
    Dens%neutralrho_axes_flux_y = (/ id_grid_xt, id_grid_yu, id_neutralrho_axis /)


    id_theta_bounds = diag_axis_init                                                       &
     ( 'theta_edges',potrho_bounds, 'C', 'z','potential or conservative temperature edges',&
       direction=1, set_name='ocean')

    id_theta_axis   = diag_axis_init                                            &
     ( 'theta',Dens%theta_ref, 'C', 'z','potential or conservative temperature',&
       edges=id_theta_bounds, direction=1,set_name='ocean')

    Dens%theta_axes        = (/ id_grid_xt,  id_grid_yt,  id_theta_axis  /)
    Dens%theta_axes_flux_x = (/ id_grid_xu,  id_grid_yt,  id_theta_axis  /)
    Dens%theta_axes_flux_y = (/ id_grid_xt,  id_grid_yu,  id_theta_axis  /)

  endif

  ! register the diagnostic fields 
  call density_diagnostics_init(Time)

  write(stdoutunit,'(/a)') ' From ocean_density_mod: density chksums from ocean_density_init'
  call write_timestamp(Time%model_time)
  call ocean_density_chksum(Time, Dens, use_blobs)

  ! clean up pointers used in this subroutine 
  nullify(salinity)
  nullify(temperature)

  end subroutine ocean_density_init
! </SUBROUTINE> NAME="ocean_density_init"


!#######################################################################
! <SUBROUTINE NAME="density_diagnostics_init">
!
! <DESCRIPTION>
! Register the diagnostic fields. 
! </DESCRIPTION>
!
  subroutine density_diagnostics_init(Time)

  type(ocean_time_type),  intent(in) :: Time

  id_press = register_diag_field ('ocean_model', 'press', Grd%tracer_axes(1:3),&
             Time%model_time, 'absolute pressure', 'dbar',                     &
             missing_value=missing_value, range=(/-10.,1e6/))

  id_rho = register_diag_field ('ocean_model', 'rho', Grd%tracer_axes(1:3), &
           Time%model_time, 'in situ density', 'kg/m^3',                    &
           missing_value=missing_value, range=(/-10.0,1e5/))

  id_rhoE = register_diag_field ('ocean_model', 'rhoE', Grd%tracer_axes(1:3), &
            Time%model_time, 'E system in situ density', 'kg/m^3',            &
            missing_value=missing_value, range=(/-10.0,1e5/))

  id_pot_rho = register_diag_field ('ocean_model', 'pot_rho', Grd%tracer_axes(1:3), &
               Time%model_time, 'potential density', 'kg/m^3',                      &
               missing_value=missing_value, range=(/-10.0,1e5/))    

  id_int_rhodz = register_diag_field ('ocean_model', 'int_rhodz', Grd%tracer_axes(1:2), &
               Time%model_time, 'vertical integral of density', 'm*(kg/m^3)',           &
               missing_value=missing_value, range=(/-10.0,1e5/))    

  if(neutral_density_omega) then   
      id_neutral_rho = register_diag_field ('ocean_model', 'neutral_rho', Grd%tracer_axes(1:3), &
                       Time%model_time, 'polynomial estimate of neutral density', 'kg/m^3',     &
                       missing_value=missing_value, range=(/-10.0,1e5/))    
  elseif(neutral_density_potrho) then 
      id_neutral_rho = register_diag_field ('ocean_model', 'neutral_rho', Grd%tracer_axes(1:3),   &
                       Time%model_time, 'potential density estimate of neutral density', 'kg/m^3',&
                       missing_value=missing_value, range=(/-10.0,1e5/))    
  elseif(neutral_density_theta) then 
      id_neutral_rho = register_diag_field ('ocean_model', 'neutral_rho', Grd%tracer_axes(1:3),   &
                       Time%model_time, 'conservative temperature is neutral density', 'degC',&
                       missing_value=missing_value, range=(/-10.0,1e5/))    
  endif 

  id_eos_salinity = register_diag_field ('ocean_model', 'eos_salinity', Grd%tracer_axes(1:3),&
           Time%model_time, 'absolute salinity as used for TEOS10 equation of state', 'g/kg',&
           missing_value=missing_value, range=(/0.0,100.0/))

  id_pot_rho_0 = register_diag_field ('ocean_model', 'pot_rho_0', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 0 dbar', 'kg/m^3',   &
                 missing_value=missing_value, range=(/-10.0,1e5/),                      &
                 standard_name='sea_water_potential_density')    
  id_pot_rho_1 = register_diag_field ('ocean_model', 'pot_rho_1', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 1000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_2 = register_diag_field ('ocean_model', 'pot_rho_2', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 2000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_3 = register_diag_field ('ocean_model', 'pot_rho_3', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 3000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_4 = register_diag_field ('ocean_model', 'pot_rho_4', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 4000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_et = register_diag_field ('ocean_model', 'pot_rho_et', Grd%tracer_axes_flux_x(1:3), &
                  Time%model_time, 'potential density at east face of tracer cell', 'kg/m^3',    &
                  missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_nt = register_diag_field ('ocean_model', 'pot_rho_nt', Grd%tracer_axes_flux_y(1:3), &
                  Time%model_time, 'potential density at north face of tracer cell', 'kg/m^3',   &
                  missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_wt = register_diag_field ('ocean_model', 'pot_rho_wt', Grd%tracer_axes_wt(1:3), &
                  Time%model_time, 'potential density at wt points', 'kg/m^3',               &
                  missing_value=missing_value, range=(/-10.0,1e5/))    
  id_thermal_expand = register_diag_field ('ocean_model', 'thermal_expand', Grd%tracer_axes(1:3), &
                  Time%model_time, 'thermal expansion coefficient -(1/rho)d(rho)/d(theta)', '1/C',&
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_haline_contract= register_diag_field ('ocean_model', 'haline_contract', Grd%tracer_axes(1:3),        & 
                  Time%model_time, 'haline contraction coefficient (1/rho)d(rho)/d(salinity)', '1/(g/kg)',&
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_drhodtheta = register_diag_field ('ocean_model', 'drhodtheta', Grd%tracer_axes(1:3), &
                  Time%model_time, 'd(rho)/d(theta)', 'kg/m^3/C',                         &
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_drhodsalt   = register_diag_field ('ocean_model', 'drhodsalinity', Grd%tracer_axes(1:3), &
                   Time%model_time, 'd(rho)/d(salinity)', 'kg/m^3/psu',                       &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_dpotrhodtheta = register_diag_field ('ocean_model', 'dpotrhodtheta', Grd%tracer_axes(1:3), &
                  Time%model_time, 'd(potrho)/d(theta)', 'kg/m^3/C',                            &
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_dpotrhodsalt   = register_diag_field ('ocean_model', 'dpotrhodsalt', Grd%tracer_axes(1:3), & 
                   Time%model_time, 'd(potrho)/d(salinity)', 'kg/m^3/psu',                      &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_dpotrhodpress  = register_diag_field ('ocean_model', 'dpotrhodpress', Grd%tracer_axes(1:3), &
                   Time%model_time, 'd(potrho)/d(press)', 'kg/m^3/psu',                          &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_drhodpress  = register_diag_field ('ocean_model', 'drhodpress', Grd%tracer_axes(1:3), &
                   Time%model_time, 'd(rho)/d(press)', '(kg/m^3)/dbar',                    &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_sound_speed2 = register_diag_field ('ocean_model', 'sound_speed2', Grd%tracer_axes(1:3),&
                   Time%model_time, 'squared sound speed = 1/[d(rho)/d(press)]', '(m/s)^2',  &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_cabbeling   = register_diag_field ('ocean_model', 'cabbeling', Grd%tracer_axes(1:3), &
                  Time%model_time, 'cabbeling parameter', '(1/degC)^2',                   &
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_thermobaricity = register_diag_field ('ocean_model', 'thermobaricity', Grd%tracer_axes(1:3), &
                  Time%model_time, 'thermobaricity parameter', '1/(dbar*degC)',                   &
                  missing_value=-1e10, range=(/-1e10,1e10/))

  id_buoyfreq2_zt = register_diag_field ('ocean_model','buoyfreq2_zt',  &
                    Grd%tracer_axes(1:3), Time%model_time,              &
                    'Squared buoyancy frequency at T-point', '1/s^2',   &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_buoyfreq2_wt = register_diag_field ('ocean_model','buoyfreq2_wt',      &
                    Grd%tracer_axes_wt(1:3), Time%model_time,               &
                    'Squared buoyancy frequency at T-cell bottom', '1/s^2', &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhodz_zt = register_diag_field ('ocean_model','drhodz_zt',&
                    Grd%tracer_axes(1:3), Time%model_time,      &
                    'd(neutral rho)/dz at T-point', 'kg/m^4',   &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhodx_zt = register_diag_field ('ocean_model','drhodx_zt',     &
                    Grd%tracer_axes(1:3), Time%model_time,           &
                    'd(neutral rho)/dx at nominal T-point', 'kg/m^4',&
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhody_zt = register_diag_field ('ocean_model','drhody_zt',     &
                    Grd%tracer_axes(1:3), Time%model_time,           &
                    'd(neutral rho)/dy at nominal T-point', 'kg/m^4',&
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhodz_diag = register_diag_field ('ocean_model','drhodz_diag',     &
                    Grd%tracer_axes(1:3), Time%model_time,               &
                    'regularized d(neutral rho)/dz at T-point', 'kg/m^4',&
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhodz_wt = register_diag_field ('ocean_model','drhodz_wt',  &
                    Grd%tracer_axes_wt(1:3), Time%model_time,     &
                    'd(neutral rho)/dz at W-point', 'kg/m^4',     &
                    missing_value=missing_value, range=(/-1e6,1e6/))

  id_rhoave    = register_diag_field('ocean_model','rhoave', Time%model_time, &
                   'global mean ocean in-situ density from ocean_density_mod',&
                   'kg/m^3', missing_value=missing_value, range=(/-1e10,1e10/))
  id_pbot_adjust = register_diag_field('ocean_model','pbot_adjust',Time%model_time,         &
                   'pbot adjustment to counteract spurious mass source in Boussinesq fluid',&
                   'dbar', missing_value=missing_value, range=(/-1e10,1e10/))
  id_eta_adjust = register_diag_field('ocean_model','eta_adjust',Time%model_time,       &
                   'global eta adjustment to include steric effect in Boussinesq fluid',&
                   'm', missing_value=missing_value, range=(/-100.0,100.0/))
  id_eta_adjust_approx = register_diag_field('ocean_model','eta_adjust_approx',Time%model_time,     &
                   'approximate global eta adjustment to include steric effect in Boussinesq fluid',&
                   'm', missing_value=missing_value, range=(/-100.0,100.0/))
  id_eta_adjust_cstvolume = register_diag_field('ocean_model','eta_adjust_cstvolume',Time%model_time,              &
                   'global eta adjustment (assuming constant volume) to include steric effect in Boussinesq fluid',&
                   'm', missing_value=missing_value, range=(/-100.0,100.0/))

  id_grad_nrho     = register_diag_field ('ocean_model','grad_nrho', &
                    Grd%tracer_axes(1:3), Time%model_time,           &
                    '|grad neutral rho|', 'dimensionless',           &
                    missing_value=missing_value, range=(/-1e1,1e10/))
  id_grad_lrpotrho = register_diag_field ('ocean_model','grad_lrpotrho', &
                    Grd%tracer_axes(1:3), Time%model_time,               &
                    '|grad local ref potential rho|', 'dimensionless',   &
                    missing_value=missing_value, range=(/-1e1,1e10/))
  id_grad_nrho_lrpotrho = register_diag_field ('ocean_model','grad_nrho_lrpotrho',         &
                    Grd%tracer_axes(1:3), Time%model_time,                                 &
                    '|grad neutral rho| / |grad local ref potential rho|', 'dimensionless',&
                    missing_value=missing_value, range=(/-1e1,1e10/))
  id_watermass_factor = register_diag_field ('ocean_model','watermass_factor',                       &
                    Grd%tracer_axes(1:3), Time%model_time,                                           &
                    '(|grad neutral rho|/|grad local ref potential rho|)/(delta(gamma))', '(m^3/kg)',&
                    missing_value=missing_value, range=(/-1e1,1e15/))
  id_stratification_factor = register_diag_field ('ocean_model','stratification_factor',         &
                    Grd%tracer_axes(1:3), Time%model_time, 'abs ( rho*Area(h)/gamma_h )', 'm^3', &
                    missing_value=missing_value, range=(/-1e1,1e22/))
  id_stratification_axis = register_diag_field ('ocean_model','stratification_axis',Grd%tracer_axes(1:3), &
                    Time%model_time, 'direction of strongest stratification', 'dimensionless',            &
                    missing_value=missing_value, range=(/-1e2,1e2/))
  id_smax_dianeutral = register_diag_field ('ocean_model','smax_dianeutral',Grd%tracer_axes(1:3),      &
                    Time%model_time, 'slope factor used for computing dianeutral transport diagnostic',&
                    'dimensionless', missing_value=missing_value, range=(/-1e2,1e2/))


  return 
  end subroutine density_diagnostics_init
! </SUBROUTINE> NAME="density_diagnostics_init"


!#######################################################################
! <SUBROUTINE NAME="density_coeffs_init">
!
! <DESCRIPTION>
! Initialize the EOS coefficients, and write some test values.  
! </DESCRIPTION>
!
  subroutine density_coeffs_init()

    integer :: stdoutunit
    real    :: rho_neutralrho
    real    :: neutralrho_test, rho_test, alpha_test, beta_test, speed2_test
    real    :: drho_dtheta_test, drho_dsal_test, drho_dpress_test
    real    :: diff_test

    stdoutunit=stdout()

    ! for the TESO10 EOS 
    mbfj_rho   = 1.017775176234136d+3
    mbfj_alpha = 2.435473441547041d-4
    mbfj_beta  = 7.284367916939847d-4
    mb_neutralrho=1033.093610463980

    v01 =  9.998420897506056d+2
    v02 =  2.839940833161907
    v03 = -3.147759265588511d-2
    v04 =  1.181805545074306d-3
    v05 = -6.698001071123802
    v06 = -2.986498947203215d-2
    v07 =  2.327859407479162d-4
    v08 = -3.988822378968490d-2
    v09 =  5.095422573880500d-4
    v10 = -1.426984671633621d-5
    v11 =  1.645039373682922d-7
    v12 = -2.233269627352527d-2
    v13 = -3.436090079851880d-4
    v14 =  3.726050720345733d-6
    v15 = -1.806789763745328d-4
    v16 =  6.876837219536232d-7
    v17 = -3.087032500374211d-7
    v18 = -1.988366587925593d-8
    v19 = -1.061519070296458d-11
    v20 =  1.550932729220080d-10
    v21 =  1.0
    v22 =  2.775927747785646d-3
    v23 = -2.349607444135925d-5
    v24 =  1.119513357486743d-6
    v25 =  6.743689325042773d-10
    v26 = -7.521448093615448d-3
    v27 = -2.764306979894411d-5
    v28 =  1.262937315098546d-7
    v29 =  9.527875081696435d-10
    v30 = -1.811147201949891d-11
    v31 = -3.303308871386421d-5
    v32 =  3.801564588876298d-7
    v33 = -7.672876869259043d-9
    v34 = -4.634182341116144d-11
    v35 =  2.681097235569143d-12
    v36 =  5.419326551148740d-6
    v37 = -2.742185394906099d-5
    v38 = -3.212746477974189d-7
    v39 =  3.191413910561627d-9
    v40 = -1.931012931541776d-12
    v41 = -1.105097577149576d-7
    v42 =  6.211426728363857d-10
    v43 = -1.119011592875110d-10
    v44 = -1.941660213148725d-11
    v45 = -1.864826425365600d-14
    v46 =  1.119522344879478d-14
    v47 = -1.200507748551599d-15
    v48 =  6.057902487546866d-17

    ! Save some multiples
    two_v03 = 2.0*v03
    two_v07 = 2.0*v07
    two_v10 = 2.0*v10
    two_v14 = 2.0*v14
    two_v17 = 2.0*v17
    two_v18 = 2.0*v18
    two_v19 = 2.0*v19
    two_v20 = 2.0*v20
    two_v23 = 2.0*v23
    two_v28 = 2.0*v28
    two_v33 = 2.0*v33
    two_v36 = 2.0*v36
    two_v39 = 2.0*v39
    two_v43 = 2.0*v43
    two_v44 = 2.0*v44
    two_v45 = 2.0*v45
    two_v46 = 2.0*v46

    three_v04 = 3.0*v04
    three_v11 = 3.0*v11
    three_v24 = 3.0*v24
    three_v29 = 3.0*v29
    three_v34 = 3.0*v34
    three_v40 = 3.0*v40
    three_v47 = 3.0*v47
    three_v48 = 3.0*v48

    four_v25 = 4.0*v25
    four_v30 = 4.0*v30
    four_v35 = 4.0*v35

    onep5_v08 = 1.5*v08
    onep5_v09 = 1.5*v09
    onep5_v10 = 1.5*v10
    onep5_v11 = 1.5*v11
    onep5_v10 = 1.5*v10
    onep5_v31 = 1.5*v31
    onep5_v32 = 1.5*v32
    onep5_v33 = 1.5*v33
    onep5_v34 = 1.5*v34
    onep5_v35 = 1.5*v35

    ! Double derivs
    six_v04 = 6.0*v04
    six_v11 = 6.0*v11
    six_v24 = 6.0*v24
    six_v29 = 6.0*v29
    six_v34 = 6.0*v34
    six_v40 = 6.0*v40
    six_v47 = 6.0*v47
    six_v48 = 6.0*v48

    twelve_v25 = 12.0*v25
    twelve_v30 = 12.0*v30
    twelve_v35 = 12.0*v35

    p75_v08 = 0.75*v08
    p75_v09 = 0.75*v09
    p75_v10 = 0.75*v10
    p75_v11 = 0.75*v11
    p75_v10 = 0.75*v10
    p75_v31 = 0.75*v31
    p75_v32 = 0.75*v32
    p75_v33 = 0.75*v33
    p75_v34 = 0.75*v34
    p75_v35 = 0.75*v35

    ! Cross derivs
    three_v10 = 3.0*v10
    three_v33 = 3.0*v33

    fourp5_v11 = 4.5*v11
    fourp5_v34 = 4.5*v34

    four_v19 = 4.0*v19
    four_v45 = 4.0*v45

    six_v35 = 6.0*v35


    ! 25 coefficients in the preTEOS10 equation of state 
    if(temp_variable==CONSERVATIVE_TEMP) then

        jmfwg_rho   = 1017.842890411975d0
        jmfwg_alpha = 2.436057013634649d-4
        jmfwg_beta  = 7.314818108935248d-4

        a0  =  9.9983912878771446d+02
        a1  =  7.0687133522652896d+00
        a2  = -2.2746841916232965d-02
        a3  =  5.6569114861400121d-04
        a4  =  2.3849975952593345d+00
        a5  =  3.1761924314867009d-04
        a6  =  1.7459053010547962d-03
        a7  =  1.2192536310173776d-02
        a8  =  2.4643435731663949d-07
        a9  =  4.0525405332794888d-06
        a10 = -2.3890831309113187d-08
        a11 = -5.9016182471196891d-12

        b0  =  1.0000000000000000d+00 
        b1  =  7.0051665739672298d-03
        b2  = -1.5040804107377016d-05 
        b3  =  5.3943915288426715d-07
        b4  =  3.3811600427083414d-10
        b5  =  1.5599507046153769d-03
        b6  = -1.8137352466500517d-06
        b7  = -3.3580158763335367d-10
        b8  =  5.7149997597561099d-06
        b9  =  7.8025873978107375d-10
        b10 =  7.1038052872522844d-06
        b11 = -2.1692301739460094d-17
        b12 = -8.2564080016458560d-18

        ! Coefficients for neutral density based on McDougall/Jackett (2005).
        ! To be replaced by Klocker/McDougall approach in near future. 
        rho_neutralrho=1024.43863927763d0
        
        a0n =  1.0022048243661291d+003
        a1n =  2.0634684367767725d-001
        a2n =  8.0483030880783291d-002
        a3n = -3.6670094757260206d-004
        a4n = -1.4602011474139313d-003
        a5n = -2.5860953752447594d-003
        a6n = -3.0498135030851449d-007

        b0n =  1.0000000000000000d+000 
        b1n =  4.4946117492521496d-005
        b2n =  7.9275128750339643d-005
        b3n = -1.2358702241599250d-007
        b4n = -4.1775515358142458d-009
        b5n = -4.3024523119324234d-004
        b6n =  6.3377762448794933d-006
        b7n = -7.2640466666916413d-010
        b8n = -5.1075068249838284d-005
        b9n = -5.8104725917890170d-009


    elseif(temp_variable==POTENTIAL_TEMP) then 
    
        jmfwg_rho   =  1017.728868019642d0
        jmfwg_alpha = 2.525481286927133d-4
        jmfwg_beta  = 7.379638527217575d-4

        a0  =  9.9984085444849347d+02
        a1  =  7.3471625860981584d+00
        a2  = -5.3211231792841769d-02
        a3  =  3.6492439109814549d-04
        a4  =  2.5880571023991390d+00
        a5  = -6.7168282786692355d-03
        a6  =  1.9203202055760151d-03
        a7  =  1.1798263740430364d-02
        a8  =  9.8920219266399117d-08
        a9  =  4.6996642771754730d-06
        a10 = -2.5862187075154352d-08
        a11 = -3.2921414007960662d-12

        b0  =  1.0000000000000000d+00 
        b1  =  7.2815210113327091d-03
        b2  = -4.4787265461983921d-05 
        b3  =  3.3851002965802430d-07
        b4  =  1.3651202389758572d-10
        b5  =  1.7632126669040377d-03
        b6  = -8.8066583251206474d-06
        b7  = -1.8832689434804897d-10
        b8  =  5.7463776745432097d-06
        b9  =  1.4716275472242334d-09
        b10 =  6.7103246285651894d-06
        b11 = -2.4461698007024582d-17
        b12 = -9.1534417604289062d-18

        ! Coefficients for neutral density based on McDougall/Jackett (2005).
        ! To be replaced by Klocker/McDougall approach in near future. 
        rho_neutralrho=1024.59416751197d0

        a0n =  1.0023063688892480d+003
        a1n =  2.2280832068441331d-001
        a2n =  8.1157118782170051d-002
        a3n = -4.3159255086706703d-004
        a4n = -1.0304537539692924d-004
        a5n = -3.1710675488863952d-003
        a6n = -1.7052298331414675d-007

        b0n =  1.0000000000000000d+000 
        b1n =  4.3907692647825900d-005
        b2n =  7.8717799560577725d-005
        b3n = -1.6212552470310961d-007
        b4n = -2.3850178558212048d-009
        b5n = -5.1268124398160734d-004
        b6n =  6.0399864718597388d-006
        b7n = -2.2744455733317707d-009
        b8n = -3.6138532339703262d-005
        b9n = -1.3409379420216683d-009
       
    endif 

    ! save some multiples of the coefficients 
    two_a2   = 2.0*a2
    three_a3 = 3.0*a3
    six_a3   = 6.0*a3
    two_a6   = 2.0*a6
    two_a8   = 2.0*a8
    two_a10  = 2.0*a10
    two_a11  = 2.0*a11
    four_a11 = 4.0*a11

    two_b2    = 2.0*b2
    three_b3  = 3.0*b3
    six_b3    = 6.0*b3
    four_b4   = 4.0*b4
    twelve_b4 = 12.0*b4
    three_b7  = 3.0*b7
    six_b7    = 6.0*b7
    onep5_b8  = 1.5*b8
    onep5_b9  = 1.5*b9
    two_b9    = 2.0*b9
    three_b9  = 3.0*b9
    two_b11   = 2.0*b11
    three_b11 = 3.0*b11
    six_b11   = 6.0*b11
    three_b12 = 3.0*b12


   ! Test values for preTEOS10 EOS 
   if(eos_preteos10) then 

      write (stdoutunit,'(/,a)') 'preTEOS10 EQUATION OF STATE TEST VALUES'

      write (stdoutunit,'(a,f6.2,a,f6.2,a,f8.2)') &
      's_test(psu) = ',s_test,', t_test(C) = ',t_test,', p_test(dbar) = ',p_test  

      rho_test = density(s_test,t_test,p_test) 
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'rho  (',s_test,',',t_test,',',p_test,') = ',rho_test,' kg/m^3'

      diff_test = rho_test-jmfwg_rho
      write (stdoutunit,' (a,e22.16,a)') 'diff from JMFWG = ',diff_test,' kg/m^3'

      call density_derivs(rho_test, s_test, t_test, p_test, &
                          drho_dtheta_test, drho_dsal_test, drho_dpress_test)

      alpha_test = -drho_dtheta_test/(epsln+rho_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'alpha(',s_test,',',t_test,',',p_test,') = ',alpha_test,' 1/C'

      diff_test = alpha_test-jmfwg_alpha
      write (stdoutunit,' (a,e22.16,a)') 'diff from JMFWG = ',diff_test,' 1/C'

      beta_test  =  drho_dsal_test  /(epsln+rho_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'beta (',s_test,',',t_test,',',p_test,') = ',beta_test,' 1/psu'

      speed2_test  =  1.0/(epsln + drho_dpress_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'squared sound speed (',s_test,',',t_test,',',p_test,') = ',speed2_test,' (m/s)^2'

      diff_test = beta_test-jmfwg_beta
      write (stdoutunit,' (a,e22.16,a)') 'diff from JMFWG = ',diff_test,' 1/psu'

      write (stdoutunit,'(/,a)') 'NEUTRAL DENSITY EQUATION TEST VALUES'

      write (stdoutunit,'(a,f6.2,a,f6.2)') &
      'sn_test(psu) = ',sn_test,', tn_test(C) = ',tn_test

      neutralrho_test = neutral_density(sn_test,tn_test) 
      write (stdoutunit,' (a,f6.2,a,f6.2,a,e22.16,a)') &
      'rho  (',sn_test,',',tn_test,') = ',neutralrho_test,' kg/m^3'

      diff_test = neutralrho_test-rho_neutralrho
      write (stdoutunit,' (a,e22.16,a)') 'diff from Klocker and McDougall test = ',diff_test,' kg/m^3'

  endif 

  ! Test values for TEOS10 EOS 
  if (eos_teos10) then

      write (stdoutunit,'(/,a)') 'TEOS10 EQUATION OF STATE TEST VALUES'

      write (stdoutunit,'(a,f6.2,a,f6.2,a,f8.2)') &
      's_test(ppt) = ',s_test,', t_test(C) = ',t_test,', p_test(dbar) = ',p_test

      rho_test = density(s_test,t_test,p_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'rho  (',s_test,',',t_test,',',p_test,') = ',rho_test,' kg/m^3'

      diff_test = rho_test-mbfj_rho
      write (stdoutunit,' (a,e22.16,a)') 'diff from MBFJ = ',diff_test,' kg/m^3'

      call density_derivs(rho_test, s_test, t_test, p_test, &
                          drho_dtheta_test, drho_dsal_test, drho_dpress_test)

      alpha_test = -drho_dtheta_test/(epsln+rho_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'alpha(',s_test,',',t_test,',',p_test,') = ',alpha_test,' 1/C'

      diff_test = alpha_test-mbfj_alpha
      write (stdoutunit,' (a,e22.16,a)') 'diff from MBFJ = ',diff_test,' 1/C'

      beta_test  =  drho_dsal_test  /(epsln+rho_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'beta (',s_test,',',t_test,',',p_test,') = ',beta_test,' 1/ppt'

      diff_test = beta_test-mbfj_beta
      write (stdoutunit,' (a,e22.16,a)') 'diff from MBFJ = ',diff_test,' 1/ppt'

      speed2_test  =  1.0/(epsln + drho_dpress_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'squared sound speed (',s_test,',',t_test,',',p_test,') = ',speed2_test,' (m/s)^2'

      write (stdoutunit,'(/,a)') 'NEUTRAL DENSITY EQUATION TEST VALUES'

      write (stdoutunit,'(a,f6.2,a,f6.2)') &
      'sn_test(ppt) = ',sn_test,', tn_test(C) = ',tn_test

      neutralrho_test = neutral_density(sn_test,tn_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,e22.16,a)') &
      'rho  (',sn_test,',',tn_test,') = ',neutralrho_test,' kg/m^3'

      diff_test = neutralrho_test-mb_neutralrho
      write (stdoutunit,' (a,e22.16,a)') 'diff from McDougall and Barker = ',diff_test,' kg/m^3'

  endif 


  return 
  end subroutine density_coeffs_init
! </SUBROUTINE> NAME="density_coeffs_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_diag">
!
! <DESCRIPTION>
! Diagnostic ocean density fields: neutral density and potential density.  
! Also send some diagnostics to diagnostic manager.  
! </DESCRIPTION>
!
  subroutine ocean_density_diag(Time, Temp, Thickness, Dens, T_diag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_prog_tracer_type),   intent(in)    :: Temp
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(inout) :: Dens
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  
  integer :: tau

  tau   = Time%tau

  if (ind_rho .gt. 0) then
    T_diag(ind_rho)%field(:,:,:) = Dens%rho(:,:,:,tau)
  endif

  ! compute neutral density for diagnostic purposes 
  Dens%neutralrho(:,:,:) = Grd%tmask(:,:,:) &
   *neutral_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau))
  if (ind_neutralrho .gt. 0) then
    T_diag(ind_neutralrho)%field(:,:,:) = Dens%neutralrho(:,:,:)
  endif

  ! compute potential density for diagnostic purposes 
  Dens%potrho(:,:,:) = Grd%tmask(:,:,:) &
   *potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), potrho_press)

  ! compute some diagnostic factors 
  call compute_diagnostic_factors(Time, Thickness, Dens,       &
                                  Dens%rho_salinity(:,:,:,tau),&
                                  Temp%field(:,:,:,tau))

  ! compute potential densities for diagnostic purposes 
  if (ind_potrho_0 .gt. 0) then
    T_diag(ind_potrho_0)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 0.0)
  endif
  if (ind_potrho_1 .gt. 0) then
    T_diag(ind_potrho_1)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 1000.0)
  endif
  if (ind_potrho_2 .gt. 0) then
    T_diag(ind_potrho_2)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 2000.0)
  endif
  if (ind_potrho_3 .gt. 0) then
    T_diag(ind_potrho_3)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 3000.0)
  endif
  if (ind_potrho_4 .gt. 0) then
    T_diag(ind_potrho_4)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 4000.0)
  endif

  call diagnose_3d(Time, Grd, id_neutral_rho, Dens%neutralrho(:,:,:))
  call diagnose_3d(Time, Grd, id_pot_rho, Dens%potrho(:,:,:))


  if (id_pot_rho_0 > 0) then 
    if (ind_potrho_0 .gt. 0) then
       call diagnose_3d(Time, Grd, id_pot_rho_0, T_diag(ind_potrho_0)%field(:,:,:))
    else
      wrk1(:,:,:) = potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 0.0)
      call diagnose_3d(Time, Grd, id_pot_rho_0, wrk1(:,:,:))
    endif
  endif

  if (id_pot_rho_1 > 0) then 
    if (ind_potrho_1 .gt. 0) then
       call diagnose_3d(Time, Grd, id_pot_rho_1, T_diag(ind_potrho_1)%field(:,:,:))
    else
      wrk1(:,:,:) = potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 1000.0)
      call diagnose_3d(Time, Grd, id_pot_rho_1, wrk1(:,:,:))
    endif
  endif

  if (id_pot_rho_2 > 0) then 
    if (ind_potrho_2 .gt. 0) then
       call diagnose_3d(Time, Grd, id_pot_rho_2, T_diag(ind_potrho_2)%field(:,:,:))
    else
      wrk1(:,:,:) = potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 2000.0)
      call diagnose_3d(Time, Grd, id_pot_rho_2, wrk1(:,:,:))
    endif
  endif

  if (id_pot_rho_3 > 0) then 
    if (ind_potrho_3 .gt. 0) then
       call diagnose_3d(Time, Grd, id_pot_rho_3, T_diag(ind_potrho_3)%field(:,:,:))
    else
      wrk1(:,:,:) = potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 3000.0)
      call diagnose_3d(Time, Grd, id_pot_rho_3, wrk1(:,:,:))
    endif
  endif

  if (id_pot_rho_4 > 0) then 
    if (ind_potrho_4 .gt. 0) then
       call diagnose_3d(Time, Grd, id_pot_rho_4, T_diag(ind_potrho_4)%field(:,:,:))
    else
      wrk1(:,:,:) = potential_density(Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), 4000.0)
      call diagnose_3d(Time, Grd, id_pot_rho_4, wrk1(:,:,:))
    endif
  endif


end subroutine ocean_density_diag
! </SUBROUTINE> NAME="ocean_density_diag"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_density_salinity">
!
! <DESCRIPTION>
!
! If TEOS-10 is being used then we need to multiply Preformed Salinity and the
! Salinity factor to obtain absolute salinity for use in the TEOS10 EOS. 
!
! If not using TEOS10 EOS, then copy the practical salinity into the
! density_salinity field for use in the preTEOS10 EOS or the linear EOS.  
!
! Note that halo values are not generally valid for taup1 until the halos
! for index_salt and index_Fdelta have been updated inside of ocean_model.F90.
!
! </DESCRIPTION>
  subroutine update_ocean_density_salinity(T_prog, time_level, Dens)

  type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
  integer,                        intent(in)    :: time_level
  type(ocean_density_type),       intent(inout) :: Dens

  integer i,j,k

  if (Dens%use_teos10) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%rho_salinity(i,j,k,time_level) = T_prog(index_salt)%field(i,j,k,time_level)* &
                    (1.0 + T_prog(index_Fdelta)%field(i,j,k,time_level))

            enddo
         enddo
      enddo
  else
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%rho_salinity(i,j,k,time_level) = T_prog(index_salt)%field(i,j,k,time_level)
            enddo
         enddo
      enddo
  endif

  end subroutine update_ocean_density_salinity
! </SUBROUTINE> NAME="update_ocean_density_salinity"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_density">
!
! <DESCRIPTION>
! Diagnose pressure_at_depth and ocean density.  
! Also send some diagnostics to diagnostic manager.  
! </DESCRIPTION>
!
  subroutine update_ocean_density(Time, Thickness, Temp, Salt, Ext_mode, Dens, L_system, use_blobs)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(in)    :: Temp
  type(ocean_prog_tracer_type),   intent(in)    :: Salt
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(inout) :: Dens
  type(ocean_lagrangian_type),    intent(in)    :: L_system
  logical,                        intent(in)    :: use_blobs
  type(time_type)                               :: next_time 
  
  integer :: i, j, k
  integer :: tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1  

  ! compute pressure at depth (dbars).
  ! do so consistently with all fields at time tau.  
  if(vert_coordinate==GEOPOTENTIAL) then 
     if (use_blobs) then
        Dens%pressure_at_depth(:,:,:) = pressure_in_dbars_blob(Thickness, Dens%rho(:,:,:,tau), &
                                        L_system%rho_dztup(:,:,:), L_system%rho_dztlo(:,:,:),  &
                                        Ext_mode%patm_t(:,:,tau)+grav*Dens%rho(:,:,1,tau)*Ext_mode%eta_t(:,:,tau))
     else
        Dens%pressure_at_depth(:,:,:) = pressure_in_dbars(Thickness, Dens%rho(:,:,:,tau), &
                                        Ext_mode%patm_t(:,:,tau)+grav*Dens%rho(:,:,1,tau)*Ext_mode%eta_t(:,:,tau))
     endif
  elseif(vert_coordinate==ZSTAR .or. vert_coordinate==ZSIGMA) then 
     if (use_blobs) then
        Dens%pressure_at_depth(:,:,:) = pressure_in_dbars_blob(Thickness, Dens%rho(:,:,:,tau), &
                                        L_system%rho_dztup(:,:,:), L_system%rho_dztlo(:,:,:),  &
                                        Ext_mode%patm_t(:,:,tau))
     else
        Dens%pressure_at_depth(:,:,:) = pressure_in_dbars(Thickness, Dens%rho(:,:,:,tau), &
                                        Ext_mode%patm_t(:,:,tau))
     endif
  elseif(vert_coordinate==PRESSURE) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%pressure_at_depth(i,j,k) = c2dbars*Thickness%depth_st(i,j,k)
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSTAR) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%pressure_at_depth(i,j,k) = c2dbars*(                   &
                      Ext_mode%patm_t(i,j,tau)                             &
                    + Thickness%depth_st(i,j,k)*Thickness%pbot0r(i,j)      &
                      *(Ext_mode%pbot_t(i,j,tau)-Ext_mode%patm_t(i,j,tau)) &
                    )
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSIGMA) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%pressure_at_depth(i,j,k) = c2dbars*(                   &
                      Ext_mode%patm_t(i,j,tau)                             &
                    + Thickness%depth_st(i,j,k)                            &
                      *(Ext_mode%pbot_t(i,j,tau)-Ext_mode%patm_t(i,j,tau)) &
                    )
            enddo
         enddo
      enddo
  endif
  
  ! diagnose in situ density 
  Dens%rho(:,:,:,taup1) = Grd%tmask(:,:,:) &
   *density(Dens%rho_salinity(:,:,:,taup1), Temp%field(:,:,:,taup1), Dens%pressure_at_depth(:,:,:))

  if(rho0_density) then 
     Dens%rho(:,:,:,taup1) = Grd%tmask(:,:,:)*rho0
  endif 

  ! compute drhodT and drhodS and drhodP
  ! do so with tau values for salinity, temperature, and pressure
  call density_derivs(Dens%rho(:,:,:,tau), Dens%rho_salinity(:,:,:,tau),    &
                      Temp%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:), &
                      Time, Dens%drhodT(:,:,:), Dens%drhodS(:,:,:), Dens%drhodP(:,:,:)) 

  ! compute dpotrhodT and dpotrhodS
  ! do so with tau values for salinity, temperature, and pressure
  call density_derivs_potrho(Time, Dens%potrho(:,:,:), Dens%rho_salinity(:,:,:,tau),&
                             Temp%field(:,:,:,tau), Dens%dpotrhodT(:,:,:),          &
                             Dens%dpotrhodS(:,:,:), Dens%dpotrhodP(:,:,:))

  ! compute buoyancy frequency for diagnostic purposes 
  if (use_blobs) then
     call compute_buoyfreq(Time, Thickness, Salt%fieldT(:,:,:), Temp%fieldT(:,:,:), Dens, use_blobs) 
  else
     ! compute buoyancy frequency for diagnostic purposes 
     call compute_buoyfreq(Time, Thickness, Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), Dens, use_blobs) 
  endif

  call compute_density_diagnostics(Time, Thickness, Temp, Ext_mode, Dens, use_blobs)
  call compute_drhodxy(Time, Dens%rho_salinity(:,:,:,tau), Temp%field(:,:,:,tau), Dens)   

end subroutine update_ocean_density
! </SUBROUTINE> NAME="update_ocean_density"


!#######################################################################
! <FUNCTION NAME="density_field">
!
! <DESCRIPTION>
! Compute density for all grid points.  
!
! Note that pressure here is 
!
! sea pressure = absolute pressure - press_standard (dbars) 
!
! and salinity is in model units (psu or g/kg).  
!
! </DESCRIPTION>
!
  function density_field (salinity, theta, press)

    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta
    real, dimension(isd:,jsd:,:), intent(in) :: press

    real, dimension(isd:ied,jsd:jed,nk) :: density_field
    real, dimension(isd:ied,jsd:jed,nk) :: pressure 

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_field): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press
    else
       pressure(:,:,:) = press(:,:,:)
    endif 


    if(eos_linear) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_field(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
             enddo
          enddo
       enddo

    elseif(eos_preteos10) then 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,k)
                t2  = t1*t1

                s1  = salinity(i,j,k)
                sp5 = sqrt(s1) 

                p1   = pressure(i,j,k) - press_standard
                p1t1 = p1*t1

                num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                     + s1*(a4 + a5*t1  + a6*s1)      &
                     + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                     + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                     + p1*(b10 + p1t1*(b11*t2 + b12*p1))
 
                denominator_r(i,j,k) = 1.0/(epsln+den) 

                density_field(i,j,k) = num*denominator_r(i,j,k)

             enddo
          enddo
       enddo

    else   ! eos_teos10

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,k)
                t2  = t1*t1

                s1  = salinity(i,j,k)
                sp5 = sqrt(s1)

                p1   = pressure(i,j,k) - press_standard
                p1t1 = p1*t1

                num = v01 + t1*(v02 + t1*(v03 + v04*t1))                &
                      + s1*(v05 + t1*(v06 + v07*t1)                     &
                      + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                      + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                      + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

                den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))                &
                       + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                       + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                       + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                       + s1*(v41 + v42*t1)                                          &
                       + p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                       &
                       + p1*(v47 + v48*t1)))

                denominator_r(i,j,k) = 1.0/(epsln+den)

                density_field(i,j,k) = num*denominator_r(i,j,k)

             enddo
          enddo
       enddo

    endif

  end function density_field
! </FUNCTION> NAME="density_field"


!#######################################################################
! <FUNCTION NAME="density_level">
!
! <DESCRIPTION>
! Compute density at a particular k-level. 
!
! Note that pressure here is 
!
! sea pressure = absolute pressure - press_standard (dbars)
!
! </DESCRIPTION>
!
  function density_level(salinity, theta, press)

    real, dimension(isd:,jsd:), intent(in) :: salinity
    real, dimension(isd:,jsd:), intent(in) :: theta
    real, dimension(isd:,jsd:), intent(in) :: press

    real, dimension(isd:ied,jsd:jed) :: density_level 
    real, dimension(isd:ied,jsd:jed) :: pressure 

    integer :: i, j
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_level): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:,:) = Grd%tmask(:,:,1)*potrho_press
    else 
       pressure(:,:) = press(:,:)
    endif 


    if(eos_linear) then

       do j=jsd,jed
          do i=isd,ied
             density_level(i,j) = rho0 - alpha_linear_eos*theta(i,j) + beta_linear_eos*salinity(i,j)
          enddo
       enddo

    elseif(eos_preteos10) then 

       do j=jsd,jed
          do i=isd,ied

             t1  = theta(i,j)
             t2  = t1*t1

             s1  = salinity(i,j)
             sp5 = sqrt(s1) 

             p1   = pressure(i,j) - press_standard
             p1t1 = p1*t1

             num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                  + s1*(a4 + a5*t1  + a6*s1)        &
                  + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

             den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                  + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                  + p1*(b10 + p1t1*(b11*t2 + b12*p1))

             density_level(i,j) = num/(epsln+den)

          enddo
       enddo

    else   ! eos_teos10

       do j=jsd,jed
          do i=isd,ied

             t1  = theta(i,j)
             t2  = t1*t1

             s1  = salinity(i,j)
             sp5 = sqrt(s1)

             p1   = pressure(i,j) - press_standard
             p1t1 = p1*t1

             num = v01 + t1*(v02 + t1*(v03 + v04*t1))                &
                   + s1*(v05 + t1*(v06 + v07*t1)                     &
                   + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                   + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                   + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

             den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))                &
                    + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                    + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                    + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                    + s1*(v41 + v42*t1)                                          &
                    + p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                       &
                    + p1*(v47 + v48*t1)))


             density_level(i,j) = num/(epsln+den)

          enddo
       enddo

    endif

  end function density_level
! </FUNCTION> NAME="density_level"


!#######################################################################
! <FUNCTION NAME="density_line">
!
! <DESCRIPTION>
! Compute density at a particular k-level and j index.  This scheme
! is used in the vectorized version of the full convection scheme. 
!
! Note that pressure here is
!
! sea pressure = absolute pressure - press_standard
!
! </DESCRIPTION>
!
  function density_line(salinity, theta, press)

    real, dimension(isd:), intent(in) :: salinity
    real, dimension(isd:), intent(in) :: theta
    real, dimension(isd:), intent(in) :: press

    real, dimension(isd:ied) :: density_line 
    real, dimension(isd:ied) :: pressure

    integer :: i
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_line): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:) = potrho_press
    else 
       pressure(:) = press(:)
    endif 


    if(eos_linear) then

        do i=isd,ied
           density_line(i) = rho0 - alpha_linear_eos*theta(i) + beta_linear_eos*salinity(i)
        enddo

    elseif(eos_preteos10) then 

        do i=isd,ied

           t1  = theta(i)
           t2  = t1*t1

           s1  = salinity(i)
           sp5 = sqrt(s1) 

           p1   = pressure(i) - press_standard
           p1t1 = p1*t1

           num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                + s1*(a4 + a5*t1  + a6*s1)        &
                + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

           den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                  + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                  + p1*(b10 + p1t1*(b11*t2 + b12*p1))

           density_line(i) = num/(epsln+den)
       enddo

    else  ! eos_teos10

        do i=isd,ied

           t1  = theta(i)
           t2  = t1*t1

           s1  = salinity(i)
           sp5 = sqrt(s1)

           p1   = pressure(i) - press_standard
           p1t1 = p1*t1

           num = v01 + t1*(v02 + t1*(v03 + v04*t1))                &
                 + s1*(v05 + t1*(v06 + v07*t1)                     &
                 + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                 + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                 + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

           den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))                &
                  + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                  + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                  + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                  + s1*(v41 + v42*t1)                                          &
                  + p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                       &
                  + p1*(v47 + v48*t1)))


           density_line(i) = num/(epsln+den)

        enddo

    endif

  end function density_line
! </FUNCTION> NAME="density_line"


!#######################################################################
! <FUNCTION NAME="neutral_density_field">
!
! <DESCRIPTION>
! Compute neutral density for use in various layer diagnostics.
!
! Two options are presently available:
!
! A/ use rational polynomial (to be done)
!
! B/ use potential density referenced to pressure potrho_press
!    McDougall recommends potential density referenced to 2000dbar,
!    since the rational polynomial is not too good.   
!
! Note that presently, the rational polynomial method defaults to 
! potential density referenced to 2000dbar.
! The polynomial approximation from McDougall and Jackett (2005)
! is not recommended (as per Trevor McDougall, 2011). A new polynomial 
! is being constructed and should be ready end of 2011.  
!
! </DESCRIPTION>
!
  function neutral_density_field (salinity, theta)

    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta

    real, dimension(isd:ied,jsd:jed,nk) :: neutral_density_field

    integer :: i, j, k

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (neutral_density_field): module must be initialized')
    endif 

    if(neutral_density_theta) then
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 neutral_density_field(i,j,k) = theta(i,j,k)
              enddo
           enddo
        enddo
    elseif(eos_linear) then
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 neutral_density_field(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
              enddo
           enddo
        enddo
    else
         neutral_density_field = potential_density(salinity,theta,potrho_press)
    endif 

  end function neutral_density_field
! </FUNCTION> NAME="neutral_density_field"


!#######################################################################
! <FUNCTION NAME="neutral_density_point">
!
! <DESCRIPTION>
! Compute neutral density for use in various layer diagnostics.
!
! Only test here the rational polynomial 
! approximation given by McDougall and Jackett (2005).
! This test needs to be updated.  
! 
! </DESCRIPTION>
!
  function neutral_density_point (salinity, theta)

    real, intent(in) :: salinity
    real, intent(in) :: theta

    real :: neutral_density_point
    real :: t1, t2, s1, sp5
    real :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (neutral_density_point): module must be initialized')
    endif 

    if(neutral_density_theta) then

        neutral_density_point = theta

    elseif(eos_linear) then

        neutral_density_point = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity

    else

        t1  = theta
        t2  = t1*t1

        s1  = salinity
        sp5 = sqrt(s1) 

        num = a0n + t1*(a1n + t1*(a2n+a3n*t1) )    &
             + s1*(a4n + a5n*t1  + a6n*s1)       

        den = b0n + t1*(b1n + t1*(b2n + t1*(b3n + t1*b4n)))      &
             + s1*(b5n + t1*(b6n + b7n*t2) + sp5*(b8n + b9n*t2)) 

        neutral_density_point = num/(epsln+den)

    endif

  end function neutral_density_point
! </FUNCTION> NAME="neutral_density_point"


!#######################################################################
! <FUNCTION NAME="potential_density">
!
! <DESCRIPTION>
! Compute potential density referenced to some given sea pressure. 
!
! Note that potential density referenced to the surface (i.e., sigma_0)
! has a zero sea pressure, so pressure=0.0 should be the argument
! to the function. 
!
! Note that pressure here is 
! sea pressure = absolute pressure - press_standard  (dbars)
!
! input pressure < 0 is an error, and model is brought down.  
! 
! </DESCRIPTION>
!
  function potential_density (salinity, theta, pressure)

    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta
    real,                         intent(in) :: pressure

    real, dimension(isd:ied,jsd:jed,nk) :: potential_density

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (potential_density): module must be initialized')
    endif 

    if(pressure < 0.0) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: potential density at pressure < 0 is not defined')
    endif 


    if(eos_linear) then

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 potential_density(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
              enddo
           enddo
        enddo

    elseif(eos_preteos10) then 

        p1 = pressure
        if(p1 > 0.0) then 

            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     t1  = theta(i,j,k)
                     t2  = t1*t1

                     s1  = salinity(i,j,k)
                     sp5 = sqrt(s1) 

                     p1t1 = p1*t1

                     num = a0 + t1*(a1 + t1*(a2+a3*t1) )    &
                          + s1*(a4 + a5*t1  + a6*s1)        &
                          + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                     den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                          + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                          + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                     potential_density(i,j,k) = num/(epsln+den)
                  enddo
               enddo
            enddo

        elseif(p1==0.0) then ! simplification for sigma_0

            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     t1  = theta(i,j,k)
                     t2  = t1*t1

                     s1  = salinity(i,j,k)
                     sp5 = sqrt(s1) 

                     num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                          + s1*(a4 + a5*t1  + a6*s1)        

                     den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                          + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) 

                     potential_density(i,j,k) = num/(epsln+den)
                  enddo
               enddo
            enddo

        endif  ! endif for value of pressure 

    else   ! eos_teos10

        p1 = pressure
        if(p1 > 0.0) then 

            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     t1  = theta(i,j,k)
                     t2  = t1*t1

                     s1  = salinity(i,j,k)
                     sp5 = sqrt(s1) 

                     p1t1 = p1*t1

                     num = v01 + t1*(v02 + t1*(v03 + v04*t1))               &
                          + s1*(v05 + t1*(v06 + v07*t1)                     &
                          + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                          + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                          + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

                     den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))              &
                          + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                          + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                          + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                          +     s1*(v41 + v42*t1)                                      &
                          +     p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                   &
                          +     p1*(v47 + v48*t1)))

                     potential_density(i,j,k) = num/(epsln+den)
                  enddo
               enddo
            enddo

        elseif(p1==0.0) then ! for sigma_0

            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     t1  = theta(i,j,k)
                     t2  = t1*t1

                     s1  = salinity(i,j,k)
                     sp5 = sqrt(s1) 

                     num = v01 + t1*(v02 + t1*(v03 + v04*t1))         &
                          + s1*(v05 + t1*(v06 + v07*t1)               &
                          + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1)))) 

                     den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))              &
                          + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                          + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))


                     potential_density(i,j,k) = num/(epsln+den)
                  enddo
               enddo
            enddo

        endif  ! endif for pressure 

    endif  ! endif for eos choice 


  end function potential_density
! </FUNCTION> NAME="potential_density"


!#######################################################################
! <SUBROUTINE NAME="compute_density_diagnostics">
!
! <DESCRIPTION>
! Diagnostics related to density. 
! </DESCRIPTION>
!
  subroutine compute_density_diagnostics(Time, Thickness, Temp, Ext_mode, Dens, use_blobs)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(in)    :: Temp
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(in)    :: Dens
  logical,                        intent(in)    :: use_blobs
  type(time_type)                               :: next_time 
  
  logical :: used
  integer :: i, j, k
  integer :: tau, taup1
  real    :: mass, massip1, massjp1
  real    :: density_tau         =0.0
  real    :: density_taup1       =0.0
  real    :: volume_tau          =0.0
  real    :: volume_taup1        =0.0
  real    :: volume_taup12       =0.0
  real    :: mass_tau            =0.0
  real    :: mass_taup1          =0.0
  real    :: pbot_adjust         =0.0
  real    :: eta_adjust          =0.0
  real    :: eta_adjust_approx   =0.0
  real    :: eta_adjust_cstvolume=0.0

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1  

  ! compute cabbeling and thermobaricity parameters for diagnostics 
  if(id_cabbeling > 0 .or. id_thermobaricity > 0) then 
      call cabbeling_thermobaricity(Dens%rho(:,:,:,tau), Dens%rho_salinity(:,:,:,tau),    &
                                    Temp%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:), &
                                    Time, Dens%drhodT(:,:,:), Dens%drhodS(:,:,:), Dens%drhodP(:,:,:))
  endif 

  call diagnose_3d(Time, Grd, id_press, Dens%pressure_at_depth(:,:,:))

  if (use_blobs) then
     call diagnose_3d(Time, Grd, id_rho, Dens%rhoT(:,:,:))
     call diagnose_3d(Time, Grd, id_rhoE, Dens%rho(:,:,:,tau))
  else
     call diagnose_3d(Time, Grd, id_rho, Dens%rho(:,:,:,tau))
  endif

  call diagnose_3d(Time, Grd, id_eos_salinity, Dens%rho_salinity(:,:,:,tau))

  if(id_int_rhodz > 0) then 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%dzt(i,j,k)*Dens%rho(i,j,k,tau) 
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_int_rhodz, wrk1_2d(:,:))
  endif

  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_pot_rho_wt,next_time)) then 
      wrk1(:,:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,nk) = Dens%potrho(i,j,nk)
         enddo
      enddo
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%tmask(i,j,k) > 0.0) then   
                   wrk1(i,j,k) = ( Grd%tmask(i,j,k+1)*Dens%potrho(i,j,k+1)*Thickness%rho_dzt(i,j,k+1,tau) &
                        +Dens%potrho(i,j,k)  *Thickness%rho_dzt(i,j,k,tau))                               &
                        /( epsln + Thickness%rho_dzt(i,j,k,tau) + Grd%tmask(i,j,k+1)*Thickness%rho_dzt(i,j,k+1,tau) ) 
               endif
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_pot_rho_wt, wrk1(:,:,:))
  endif

  if (need_data(id_pot_rho_et,next_time)) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               mass      = Thickness%rho_dzt(i,j,k,tau)  *Grd%dat(i,j)  *Grd%tmask(i,j,k) 
               massip1   = Thickness%rho_dzt(i+1,j,k,tau)*Grd%dat(i+1,j)*Grd%tmask(i+1,j,k) 
               wrk1(i,j,k) = ( Dens%potrho(i,j,k)*mass + Dens%potrho(i+1,j,k)*massip1 ) &
                    /( epsln + mass + massip1 ) 
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_pot_rho_et, wrk1(:,:,:))
  endif

  if (need_data(id_pot_rho_nt,next_time)) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               mass      = Thickness%rho_dzt(i,j,k,tau)  *Grd%dat(i,j)  *Grd%tmask(i,j,k) 
               massjp1   = Thickness%rho_dzt(i,j+1,k,tau)*Grd%dat(i,j+1)*Grd%tmask(i,j+1,k) 
               wrk1(i,j,k) = ( Dens%potrho(i,j,k)*mass + Dens%potrho(i,j+1,k)*massjp1 ) &
                    /( epsln + mass + massjp1 ) 
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_pot_rho_nt, wrk1(:,:,:))
  endif


  if (id_pbot_adjust>0       .or. id_eta_adjust>0           .or. &
      id_eta_adjust_approx>0 .or. id_eta_adjust_cstvolume>0 .or. id_rhoave>0) then 
    ! These diagnostics are used to adjust the bottom pressure and surface height computed
    ! from a Boussinesq model.  The adjustments correct, globally, for the spurious 
    ! mass sources appearing in the Boussinesq fluid and the missing steric effect.  
    ! The adjustments are NOT needed with non-Boussinesq models (pressure based vert coord),
    ! since the non-Boussinesq fluid correctly computes the bottom pressure and 
    ! surface height according to mass conserving kinematics.
    !
    ! In a Bouss model (depth based vert coordinates), when global mean density increases, 
    ! so does the Boussinesq mass. So the pbot_adjust is negative in this case, in 
    ! order to globally counteract effects from the spurious mass source.
    !
    ! In a non-Bouss model, when global mean density increases, global mean sea level decreases. 
    ! This global steric effect is absent from the Boussinesq fluid.  To incorporate this 
    ! this missing contribution to sea level, eta_adjust will be negative in this case. 
    !
    ! note that dzt has already been updated to its taup1 value. that is why we use  
    ! rho0r*rho_dzt = dzt (when depth-based vertical coordinates).

      mass_tau         = 0.0
      mass_taup1       = 0.0
      volume_tau       = 0.0
      volume_taup1     = 0.0
      volume_taup12    = 0.0
      rhodz_tau(:,:)   = 0.0
      rhodz_taup1(:,:) = 0.0
      wrk1_2d(:,:)     = 0.0
      wrk2_2d(:,:)     = 0.0
      wrk3_2d(:,:)     = 0.0
      wrk4_2d(:,:)     = 0.0

      ! do the vertical sum of rho*dzt 
      if(vert_coordinate==GEOPOTENTIAL .or. vert_coordinate==ZSTAR .or. vert_coordinate==ZSIGMA) then 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   rhodz_tau(i,j)  = rhodz_tau(i,j)  + Grd%tmask(i,j,k)*rho0r*Thickness%rho_dzt(i,j,k,tau)   &
                                                       *Dens%rho(i,j,k,tau)
                   rhodz_taup1(i,j)= rhodz_taup1(i,j)+ Grd%tmask(i,j,k)*rho0r*Thickness%rho_dzt(i,j,k,taup1) &
                                                       *Dens%rho(i,j,k,taup1)
                enddo
             enddo
          enddo
      else 
          k=1
          do j=jsc,jec
             do i=isc,iec
                rhodz_tau(i,j)   = Grd%tmask(i,j,k)*grav_r*(Ext_mode%pbot_t(i,j,tau)  -Ext_mode%patm_t(i,j,tau))
                rhodz_taup1(i,j) = Grd%tmask(i,j,k)*grav_r*(Ext_mode%pbot_t(i,j,taup1)-Ext_mode%patm_t(i,j,taup1))
             enddo
          enddo
      endif

      ! compute the mass and volume of a seawater column
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*rhodz_tau(i,j)
            wrk2_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*rhodz_taup1(i,j)
            wrk3_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau))
            wrk4_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,taup1))
         enddo
      enddo

      ! perform global sums to get global mass and volume for ocean 
      mass_tau    = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:),NON_BITWISE_EXACT_SUM)
      volume_tau  = mpp_global_sum(Dom%domain2d,wrk3_2d(:,:),NON_BITWISE_EXACT_SUM)
      density_tau = mass_tau/(epsln+volume_tau)

      mass_taup1    = mpp_global_sum(Dom%domain2d,wrk2_2d(:,:),NON_BITWISE_EXACT_SUM)
      volume_taup1  = mpp_global_sum(Dom%domain2d,wrk4_2d(:,:),NON_BITWISE_EXACT_SUM)
      density_taup1 = mass_taup1/(epsln+volume_taup1)

      volume_taup12 = 0.5*(volume_taup1+volume_tau)

      if(id_pbot_adjust > 0) then 
         pbot_adjust = -grav*c2dbars*area_total_r*volume_taup1*(density_taup1-density_tau) 
         used = send_data(id_pbot_adjust, pbot_adjust, Time%model_time)
      endif 
      if(id_eta_adjust > 0) then 
         eta_adjust = area_total_r*volume_taup12*log((density_tau+epsln)/(density_taup1+epsln))
         used = send_data(id_eta_adjust, eta_adjust, Time%model_time)
      endif 
      if(id_eta_adjust_approx > 0) then 
         eta_adjust_approx = area_total_r*volume_taup1*(density_tau/(epsln+density_taup1)-1.0)
         used = send_data(id_eta_adjust_approx, eta_adjust_approx, Time%model_time)
      endif 
      if(id_eta_adjust_cstvolume > 0) then 
         eta_adjust_cstvolume = area_total_r*Grd%tcellv*(density_tau/(epsln+density_taup1)-1.0)
         used = send_data(id_eta_adjust_cstvolume, eta_adjust_cstvolume, Time%model_time)
      endif 
      if(id_rhoave > 0) then 
         used = send_data(id_rhoave, density_tau, Time%model_time)
      endif    

  endif 

  if(debug_this_module) then 
      write(stdoutunit,'(a)') ' ' 
      write(stdoutunit,*) 'From ocean_density_mod: density chksums from compute_density_diagnostics'
      call write_timestamp(Time%model_time)
      call ocean_density_chksum(Time, Dens, use_blobs)
  endif

end subroutine compute_density_diagnostics
! </SUBROUTINE> NAME="compute_density_diagnostics"



!#######################################################################
! <SUBROUTINE NAME="compute_diagnostic_factors">
!
! <DESCRIPTION>
!
! 1/ Compute ratio |grad neutral rho| / |grad local ref pot rho|
! for use in tform water mass analysis as per Iudicone et al. (2008).
!
! 2/ Compute rho*Area(h)/gamma_h, where "h" is the direction where
! gamma has the largest stratification, and where gamma is the 
! locally referenced potential density.
!
! </DESCRIPTION>
!
  subroutine compute_diagnostic_factors (Time, Thickness, Dens, salinity, theta)

    type(ocean_time_type),         intent(in)    :: Time
    type(ocean_thickness_type),    intent(in)    :: Thickness 
    type(ocean_density_type),      intent(inout) :: Dens
    real, dimension(isd:,jsd:,:),  intent(in)    :: salinity
    real, dimension(isd:,jsd:,:),  intent(in)    :: theta

    integer :: i, j, k, kp1, km1, tau
    real    :: tmp, active_cells 

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (compute_diagnostic_factors): module must be initialized')
    endif 

    tau = Time%tau 

    if (.not. update_diagnostic_factors) then 
        do k=1,nk
            do j=jsd,jed
               do i=isd,ied                
                  Dens%watermass_factor(i,j,k)      = neutralrho_interval_r*Grd%tmask(i,j,k)
                  Dens%stratification_factor(i,j,k) = Grd%dxt(i,j)*Grd%dyt(i,j)*Thickness%dzt(i,j,k)
                  Dens%rho_dztr_tau(i,j,k)          = Thickness%rho_dztr(i,j,k)
               enddo
            enddo
        enddo
        return 
    endif 

    ! based on thinnest cell in a column, which is typically at k=1 or partial bottom 
    if(smax_min_in_column) then 
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 smax(i,j,k) = Grd%tmask(i,j,k)*min(Thickness%dzt(i,j,k),Thickness%dzt(i,j,1))*dhorz_r(i,j)
              enddo
           enddo
        enddo
    ! straightforward aspect ratio calculation (generally the recommended method)
    else  
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 smax(i,j,k) = Grd%tmask(i,j,k)*Thickness%dzt(i,j,k)*dhorz_r(i,j)
              enddo
           enddo
        enddo
    endif

    ! over-ride the above using a prescribed smax 
    if(smax_diag > 0.0) then  
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 smax(i,j,k) = min(1.0,smax_diag) 
              enddo
           enddo
        enddo
    endif 

    call diagnose_3d(Time, Grd, id_smax_dianeutral, smax(:,:,:))


    ! some initialization 
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied                
             wrk1_v(i,j,k,1) = 0.0
             wrk1_v(i,j,k,2) = 0.0
             wrk2_v(i,j,k,1) = 0.0
             wrk2_v(i,j,k,2) = 0.0
             wrk3_v(i,j,k,1) = 0.0
             wrk3_v(i,j,k,2) = 0.0
             wrk4_v(i,j,k,1) = 0.0
             wrk4_v(i,j,k,2) = 0.0
             wrk1(i,j,k)     = 0.0
             wrk2(i,j,k)     = 0.0
             wrk3(i,j,k)     = 0.0
             Dens%watermass_factor(i,j,k) = neutralrho_interval_r*Grd%tmask(i,j,k)
             wrk4(i,j,k)                  = Grd%tmask(i,j,k)/Thickness%dzwt(i,j,k) 
             Dens%rho_dztr_tau(i,j,k)     = Thickness%rho_dztr(i,j,k)
          enddo
       enddo
    enddo

    ! spatial derivatives of locally referenced potential density  
    do k=1,nk
       kp1 = min(k+1,nk)
       wrk1_v(:,:,k,1) = FDX_T(theta(:,:,k))   *FMX(Grd%tmask(:,:,k))
       wrk1_v(:,:,k,2) = FDX_T(salinity(:,:,k))*FMX(Grd%tmask(:,:,k))
       wrk2_v(:,:,k,1) = FDY_T(theta(:,:,k))   *FMY(Grd%tmask(:,:,k))
       wrk2_v(:,:,k,2) = FDY_T(salinity(:,:,k))*FMY(Grd%tmask(:,:,k))
       do j=jsc,jec
          do i=isc,iec
             wrk3_v(i,j,k,1) = Grd%tmask(i,j,kp1)*wrk4(i,j,k)*(theta(i,j,k)-theta(i,j,kp1))
             wrk3_v(i,j,k,2) = Grd%tmask(i,j,kp1)*wrk4(i,j,k)*(salinity(i,j,k)-salinity(i,j,kp1))
             wrk1(i,j,k)     = abs( Dens%drhodT(i,j,k)*wrk1_v(i,j,k,1) + Dens%drhodS(i,j,k)*wrk1_v(i,j,k,2) )
             wrk2(i,j,k)     = abs( Dens%drhodT(i,j,k)*wrk2_v(i,j,k,1) + Dens%drhodS(i,j,k)*wrk2_v(i,j,k,2) )
             wrk3(i,j,k)     = abs( Dens%drhodT(i,j,k)*wrk3_v(i,j,k,1) + Dens%drhodS(i,j,k)*wrk3_v(i,j,k,2) )
             wrk4_v(i,j,k,2) = sqrt(wrk1(i,j,k)**2 + wrk2(i,j,k)**2 + wrk3(i,j,k)**2) + epsln_drhodz
          enddo
       enddo
    enddo

    ! compute |grad neutral density|
    wrk5(:,:,:) = 0.0 
    wrk6(:,:,:) = 0.0 
    do k=1,nk
       kp1 = min(k+1,nk)
       wrk5(:,:,k) = FDX_T(Dens%neutralrho(:,:,k))*FMX(Grd%tmask(:,:,k))
       wrk6(:,:,k) = FDY_T(Dens%neutralrho(:,:,k))*FMY(Grd%tmask(:,:,k))
       do j=jsc,jec
          do i=isc,iec
             tmp             = Grd%tmask(i,j,kp1)*wrk4(i,j,k)*(Dens%neutralrho(i,j,k)-Dens%neutralrho(i,j,kp1))
             wrk4_v(i,j,k,1) = sqrt(wrk5(i,j,k)**2 + wrk6(i,j,k)**2 + tmp**2)
          enddo
       enddo
    enddo
  
    ! compute watermass factor 
    if(.not. eos_linear .and. grad_nrho_lrpotrho_compute) then
        wrk5(:,:,:) = 0.0 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk5(i,j,k) = max(grad_nrho_lrpotrho_min, &
                      min(grad_nrho_lrpotrho_max,Grd%tmask(i,j,k)*wrk4_v(i,j,k,1)/wrk4_v(i,j,k,2)))
                 Dens%watermass_factor(i,j,k) = wrk5(i,j,k)*Dens%watermass_factor(i,j,k)
              enddo
           enddo
        enddo
    endif

    call diagnose_3d(Time, Grd, id_grad_nrho, wrk4_v(:,:,:,1))
    call diagnose_3d(Time, Grd, id_grad_lrpotrho, wrk4_v(:,:,:,2))
    call diagnose_3d(Time, Grd, id_grad_nrho_lrpotrho, wrk5(:,:,:))
    call diagnose_3d(Time, Grd, id_watermass_factor, Dens%watermass_factor(:,:,:))


    ! For stratification_factor, move derivatives to centre
    ! of a grid cell through 2-point averaging. 
    call mpp_update_domains (wrk1(:,:,:), Dom%domain2d, complete=.false.)
    call mpp_update_domains (wrk2(:,:,:), Dom%domain2d, complete=.true.)
    wrk4(:,:,:) = 0.0
    wrk5(:,:,:) = 0.0
    wrk6(:,:,:) = 0.0
    do k=1,nk
       kp1 = min(k+1,nk)
       km1 = max(k-1,1)
       do j=jsc,jec
          do i=isc,iec
             active_cells = epsln + Grd%tmask(i-1,j,k) + Grd%tmask(i+1,j,k)
             wrk4(i,j,k)  = Grd%tmask(i,j,k)*(wrk1(i-1,j,k)+wrk1(i,j,k))/active_cells 

             active_cells = epsln + Grd%tmask(i,j-1,k) + Grd%tmask(i,j+1,k)
             wrk5(i,j,k)  = Grd%tmask(i,j,k)*(wrk2(i,j-1,k)+wrk2(i,j,k))/active_cells 

             active_cells = epsln + Grd%tmask(i,j,km1) + Grd%tmask(i,j,kp1)
             wrk6(i,j,k)  = Grd%tmask(i,j,k)*(wrk3(i,j,km1)+wrk3(i,j,k))/active_cells 
          enddo
       enddo
    enddo

    ! smoothing may be useful in some cases
    if(smooth_stratification_factor) then 
       call mpp_update_domains (wrk4(:,:,:), Dom%domain2d, complete=.false.)
       call mpp_update_domains (wrk5(:,:,:), Dom%domain2d, complete=.false.)
       call mpp_update_domains (wrk6(:,:,:), Dom%domain2d, complete=.true.)
       do k=1,nk
          wrk4(:,:,k) = S2D(wrk4(:,:,k))  
          wrk5(:,:,k) = S2D(wrk5(:,:,k))  
          wrk6(:,:,k) = S2D(wrk6(:,:,k))  
       enddo 
    endif 

    ! compute stratification factor. 
    ! Vertical stratification is scaled down by smax in the if-tests 
    ! to account for grid aspect ratio, and/or specific setting for smax_diag.
    wrk1(:,:,:) = 0.0
    do k=1,nk
       do j=jsc,jec
          do i=isc,iec
             if(wrk4(i,j,k) >= wrk5(i,j,k) .and. wrk4(i,j,k) >= smax(i,j,k)*wrk6(i,j,k) .and. wrk4(i,j,k) > 0.0) then 
                Dens%stratification_factor(i,j,k) = Dens%rho(i,j,k,tau)*Grd%dyt(i,j)*Thickness%dzt(i,j,k)/wrk4(i,j,k)
                wrk1(i,j,k) = 1.0
             elseif(wrk5(i,j,k) >= smax(i,j,k)*wrk6(i,j,k) .and. wrk5(i,j,k) >= wrk4(i,j,k) .and. wrk5(i,j,k) > 0.0) then 
                Dens%stratification_factor(i,j,k) = Dens%rho(i,j,k,tau)*Grd%dxt(i,j)*Thickness%dzt(i,j,k)/wrk5(i,j,k)
                wrk1(i,j,k) = 2.0
             elseif(smax(i,j,k)*wrk6(i,j,k) >= wrk4(i,j,k) .and. smax(i,j,k)*wrk6(i,j,k) >= wrk5(i,j,k) .and. wrk6(i,j,k) > 0.0) then 
                Dens%stratification_factor(i,j,k) = Dens%rho(i,j,k,tau)*Grd%dxt(i,j)*Grd%dyt(i,j)/wrk6(i,j,k)
                wrk1(i,j,k) = 3.0
             else
                Dens%stratification_factor(i,j,k) = 0.0
                wrk1(i,j,k) = 0.0
             endif  
          enddo
       enddo
    enddo

    call diagnose_3d(Time, Grd, id_stratification_factor, Dens%stratification_factor(:,:,:))
    call diagnose_3d(TIme, Grd, id_stratification_axis, wrk1(:,:,:))


  end subroutine compute_diagnostic_factors
! </SUBROUTINE> NAME="compute_diagnostic_factors"


!#######################################################################
! <FUNCTION NAME="density_sfc">
!
! <DESCRIPTION>
! Compute density as a function of surface salinity, surface theta, 
! and in situ sea pressure. 
!
! Note that pressure here is 
! sea pressure = absolute pressure - press_standard  (dbars)
! 
! For use in KPP mixed layer scheme 
! </DESCRIPTION>
!
  function density_sfc (salinity, theta, press)

    real, intent(in), dimension(isd:,jsd:,:) :: salinity
    real, intent(in), dimension(isd:,jsd:,:) :: theta
    real, intent(in), dimension(isd:,jsd:,:) :: press

    real, dimension(isd:ied,jsd:jed,nk) :: density_sfc
    real, dimension(isd:ied,jsd:jed,nk) :: pressure 

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_density_mod (density_sfc): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press
    else 
       pressure(:,:,:) = press(:,:,:)
    endif 


    if(eos_linear) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_sfc(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,1) + beta_linear_eos*salinity(i,j,1)
             enddo
          enddo
       enddo

    elseif(eos_preteos10) then 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,1)
                t2  = t1*t1

                s1  = salinity(i,j,1)
                sp5 = sqrt(s1) 

                p1   = pressure(i,j,k) - press_standard 
                p1t1 = p1*t1

                num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                     + s1*(a4 + a5*t1  + a6*s1)        &
                     + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                     + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                     + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                density_sfc(i,j,k) = num/(epsln + den)

             enddo
          enddo
       enddo

    else ! eos_teos10 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,1)
                t2  = t1*t1

                s1  = salinity(i,j,1)
                sp5 = sqrt(s1)

                p1   = pressure(i,j,k) - press_standard
                p1t1 = p1*t1

                num = v01 + t1*(v02 + t1*(v03 + v04*t1))                 &
                       + s1*(v05 + t1*(v06 + v07*t1)                     &
                       + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                       + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                       + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

                 den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))                &
                        + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                        + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                        + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                        +     s1*(v41 + v42*t1)                                      &
                        +     p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                   &
                        +     p1*(v47 + v48*t1)))
                density_sfc(i,j,k) = num/(epsln + den)

             enddo
          enddo
       enddo

    endif

  end function density_sfc
! </FUNCTION> NAME="density_sfc"


!#######################################################################
! <FUNCTION NAME="density_point">
!
! <DESCRIPTION>
! Compute density at a single model grid point. 
!
! Note that pressure here is 
!
! sea pressure = absolute pressure - press_standard  (dbars)
!
! </DESCRIPTION>
!
  function density_point (s1, t1, p1_dbars)

    real, intent(in) :: s1
    real, intent(in) :: t1
    real, intent(in) :: p1_dbars

    real :: t2, sp5, p1, p1t1
    real :: num, den
    real :: density_point
    real :: pressure 

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_point): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure = potrho_press
    else 
       pressure = p1_dbars 
    endif 


    if(eos_linear) then

       density_point = rho0 - alpha_linear_eos*t1 + beta_linear_eos*s1

    elseif(eos_preteos10) then 

       t2  = t1*t1
       sp5 = sqrt(s1) 

       p1   = pressure - press_standard 
       p1t1 = p1*t1

       num = a0 + t1*(a1 + t1*(a2+a3*t1))                   &
            + s1*(a4 + a5*t1  + a6*s1)                      &
            + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

       den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4 )))      &
            + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2))  &  
            + p1*(b10 + p1t1*(b11*t2 + b12*p1))

       density_point = num/(epsln+den)


    else ! eos_teos10 

       t2  = t1*t1
       sp5 = sqrt(s1) 

       p1   = pressure - press_standard
       p1t1 = p1*t1

       num = v01 + t1*(v02 + t1*(v03 + v04*t1))                &
             + s1*(v05 + t1*(v06 + v07*t1)                     &
             + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
             + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
             + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

       den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))                &
              + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
              + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
              + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
              +     s1*(v41 + v42*t1)                                      &
              +     p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                   &
              +     p1*(v47 + v48*t1)))


       density_point = num/(epsln+den)


    endif

  end function density_point
! </FUNCTION> NAME="density_point"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_field">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to 
! temperature and with respect to salinity.  Hold pressure constant.  
!
! Pressure here is 
!
! sea pressure = absolute press - press_standard (dbars)
!
! </DESCRIPTION>
!
  subroutine density_derivs_field (rho, salinity, theta, press, Time, &
             density_theta, density_salinity, density_press)

    type(ocean_time_type),        intent(in)    :: Time
    real, dimension(isd:,jsd:,:), intent(in)    :: rho
    real, dimension(isd:,jsd:,:), intent(in)    :: salinity
    real, dimension(isd:,jsd:,:), intent(in)    :: theta
    real, dimension(isd:,jsd:,:), intent(in)    :: press
    real, dimension(isd:,jsd:,:), intent(inout) :: density_theta
    real, dimension(isd:,jsd:,:), intent(inout) :: density_salinity
    real, dimension(isd:,jsd:,:), intent(inout) :: density_press
    real, dimension(isd:ied,jsd:jed,nk)         :: pressure 

    integer :: i, j, k
    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta 
    real :: dnum_dsalinity, dden_dsalinity
    real :: dnum_dpress, dden_dpress

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_derivs_field): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press
    else 
       pressure(:,:,:) = press(:,:,:)
    endif 

    if(eos_linear) then

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 density_theta(i,j,k)    = -alpha_linear_eos
                 density_salinity(i,j,k) =  beta_linear_eos
                 density_press(i,j,k)    =  0.0
              enddo
           enddo
        enddo

    elseif(eos_preteos10) then  

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied

                 t1  = theta(i,j,k)
                 t2  = t1*t1
                 s1  = salinity(i,j,k)
                 sp5 = sqrt(s1) 

                 p1   = pressure(i,j,k) - press_standard 
                 p1t1 = p1*t1

                 dnum_dtheta = a1 + t1*(two_a2 + three_a3*t1) &
                      + a5*s1                                 &
                      + p1t1*(two_a8 + two_a11*p1)    
                 dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                      + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                      + p1*p1*(three_b11*t2 + b12*p1)

                 dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
                 dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

                 dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
                 dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

                 density_theta(i,j,k)    = denominator_r(i,j,k)*(dnum_dtheta    - rho(i,j,k)*dden_dtheta)
                 density_salinity(i,j,k) = denominator_r(i,j,k)*(dnum_dsalinity - rho(i,j,k)*dden_dsalinity)
                 density_press(i,j,k)    = denominator_r(i,j,k)*(dnum_dpress    - rho(i,j,k)*dden_dpress)

              enddo
           enddo
        enddo

    else  ! eos_teos10 

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied

                 t1  = theta(i,j,k)
                 t2  = t1*t1
                 s1  = salinity(i,j,k)
                 sp5 = sqrt(s1)

                 p1   = pressure(i,j,k) - press_standard
                 p1t1 = p1*t1

                 dnum_dtheta = v02 + t1*(two_v03 + three_v04*t1)     &
                      + s1*((v06+two_v07*t1)                         &
                      +     sp5*(v09 + t1*(two_v10 + three_v11*t1))) &
                      + p1*((v13 + two_v14*t1) + s1*v16              &
                      +     p1*(v18+two_v19*t1))

                 dden_dtheta = v22 + t1*(two_v23 +t1*(three_v24+four_v25*t1))      &
                      + s1*(v27 + t1*(two_v28 + t1*(three_v29  + four_v30*t1))     &
                      +     sp5*(v32+t1*(two_v33 + t1*(three_v34 + four_v35*t1)))) &
                      + p1*((v38 + t1*(two_v39 + three_v40*t1))                    &
                      +     s1*v42                                                 &
                      +     p1*((v44 + two_v45*t1 + v46*s1)                        &
                      +     p1*v48))

                 dnum_dsalinity =                                    &
                      + (v05 + t1*(v06 + v07*t1)                     &
                      + 1.5*sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))&
                      + p1*((v15 + v16*t1) &
                      +      p1*v20)

                 dden_dsalinity =                                                 &
                      (v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + two_v36*s1 &
                      + 1.5*sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))  &
                      + p1*( v41 + v42*t1 + p1*t1*v46 )

                 dnum_dpress = v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1)    &
                      + p1*(two_v17 + t1*(two_v18 + two_v19*t1) + two_v20*s1) 
                 dden_dpress = (v37 + t1*(v38 + t1*(v39 + v40*t1))            &
                      + s1*(v41 + v42*t1)                                     &
                      + p1*(two_v43 + t1*(two_v44 + two_v45*t1 + two_v46*s1)  &
                      + p1*(three_v47 + three_v48*t1)))

                 density_theta(i,j,k)    = denominator_r(i,j,k)*(dnum_dtheta    - rho(i,j,k)*dden_dtheta)
                 density_salinity(i,j,k) = denominator_r(i,j,k)*(dnum_dsalinity - rho(i,j,k)*dden_dsalinity)
                 density_press(i,j,k)    = denominator_r(i,j,k)*(dnum_dpress    - rho(i,j,k)*dden_dpress)*c2dbars

              enddo
           enddo
        enddo

    endif

    call diagnose_3d(TIme, Grd, id_drhodtheta, density_theta(:,:,:))
    call diagnose_3d(Time, Grd, id_drhodsalt, density_salinity(:,:,:))
    call diagnose_3d(Time, Grd, id_drhodpress, density_press(:,:,:))
    if (id_sound_speed2 > 0) then 
        wrk1(:,:,:) = 0.0  
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if(density_press(i,j,k) /= 0.0) then 
                     wrk1(i,j,k) = 1.0/(c2dbars*density_press(i,j,k))
                 endif
              enddo
           enddo
        enddo
        call diagnose_3d(Time, Grd, id_sound_speed2, wrk1(:,:,:))
    endif
    if(id_thermal_expand > 0) then 
        wrk1(:,:,:) = 0.0  
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = -Grd%tmask(i,j,k)*density_theta(i,j,k)/(epsln+rho(i,j,k))
              enddo
           enddo
        enddo
        call diagnose_3d(Time, Grd, id_thermal_expand, wrk1(:,:,:))
    endif 
    if(id_haline_contract > 0) then 
        wrk1(:,:,:) = 0.0  
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = Grd%tmask(i,j,k)*density_salinity(i,j,k)/(epsln+rho(i,j,k))
              enddo
           enddo
        enddo
        call diagnose_3d(Time, Grd, id_haline_contract, wrk1(:,:,:))
    endif 


  end subroutine density_derivs_field
! </SUBROUTINE> NAME="density_derivs_field"

!#######################################################################
! <SUBROUTINE NAME="density_derivs_pothro">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to
! temperature and with respect to salinity.  Use potrho_press
! reference depth for potential density rather than neutral.
!
! These coefficients are useful for certain diagnostics, such as APE
! calculations.
!
! The algorithm is identical to density_derivs_field, only here we set
! the pressure to a constant reference pressure value, potrho_press.
!
! Pressure here is
! sea pressure = absolute press - press_standard (dbars)
!
! </DESCRIPTION>
!
  subroutine density_derivs_potrho (Time, rho, salinity, theta, &
             density_theta, density_salinity, density_press)

    type(ocean_time_type),        intent(in)    :: Time
    real, dimension(isd:,jsd:,:), intent(in)    :: rho
    real, dimension(isd:,jsd:,:), intent(in)    :: salinity
    real, dimension(isd:,jsd:,:), intent(in)    :: theta
    real, dimension(isd:,jsd:,:), intent(inout) :: density_theta
    real, dimension(isd:,jsd:,:), intent(inout) :: density_salinity
    real, dimension(isd:,jsd:,:), intent(inout) :: density_press
    real, dimension(isd:ied,jsd:jed,nk)         :: pressure

    integer :: i, j, k
    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta
    real :: dnum_dsalinity, dden_dsalinity
    real :: dnum_dpress, dden_dpress

    if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_derivs_field): module must be initialized')
    endif

    pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press

    if(eos_linear) then

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 density_theta(i,j,k)    = -alpha_linear_eos
                 density_salinity(i,j,k) =  beta_linear_eos
                 density_press(i,j,k)    =  0.0
              enddo
           enddo
        enddo

    elseif(eos_preteos10) then

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied

                 t1  = theta(i,j,k)
                 t2  = t1*t1
                 s1  = salinity(i,j,k)
                 sp5 = sqrt(s1)

                 p1   = pressure(i,j,k) - press_standard
                 p1t1 = p1*t1

                 dnum_dtheta = a1 + t1*(two_a2 + three_a3*t1) &
                      + a5*s1                                 &
                      + p1t1*(two_a8 + two_a11*p1)
                 dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                      + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                      + p1*p1*(three_b11*t2 + b12*p1)

                 dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
                 dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

                 dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
                 dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1)

                 density_theta(i,j,k)    = denominator_r(i,j,k)*(dnum_dtheta    - rho(i,j,k)*dden_dtheta)
                 density_salinity(i,j,k) = denominator_r(i,j,k)*(dnum_dsalinity - rho(i,j,k)*dden_dsalinity)
                 density_press(i,j,k)    = denominator_r(i,j,k)*(dnum_dpress    - rho(i,j,k)*dden_dpress)

              enddo
           enddo
        enddo

    else  ! eos_teos10

        do k=1,nk
           do j=jsd,jed
              do i=isd,ied

                 t1  = theta(i,j,k)
                 t2  = t1*t1
                 s1  = salinity(i,j,k)
                 sp5 = sqrt(s1)

                 p1   = pressure(i,j,k) - press_standard
                 p1t1 = p1*t1

                 dnum_dtheta = v02 + t1*(two_v03 + three_v04*t1)     &
                      + s1*((v06+two_v07*t1)                         &
                      +     sp5*(v09 + t1*(two_v10 + three_v11*t1))) &
                      + p1*((v13 + two_v14*t1) + s1*v16              &
                      +     p1*(v18+two_v19*t1))

                 dden_dtheta = v22 + t1*(two_v23 +t1*(three_v24+four_v25*t1))       &
                      + s1*(v27 + t1*(two_v28 + t1*(three_v29  + four_v30*t1))      &
                      +     sp5*(v32+t1*(two_v33 + t1*(three_v34 +four_v35*t1))))   &
                      + p1*((v38 + t1*(two_v39 + three_v40*t1))                     &
                      +     s1*v42                                                  &
                      +     p1*((v44 + two_v45*t1 + v46*s1)                         &
                      +     p1*v48))

                 dnum_dsalinity =                                    &
                      + (v05 + t1*(v06 + v07*t1)                     &
                      + 1.5*sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))&
                      + p1*((v15 + v16*t1) &
                      +      p1*v20)

                 dden_dsalinity =                                                 &
                      (v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + two_v36*s1 &
                      + 1.5*sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))  &
                      + p1*( v41 + v42*t1 + p1*t1*v46 )

                 dnum_dpress = v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1)    &
                      + p1*(two_v17 + t1*(two_v18 + two_v19*t1) + two_v20*s1)
                 dden_dpress = (v37 + t1*(v38 + t1*(v39 + v40*t1))            &
                      + s1*(v41 + v42*t1)                                     &
                      + p1*(two_v43 + t1*(two_v44 + two_v45*t1 + two_v46*s1)  &
                      + p1*(three_v47 + three_v48*t1)))

                 density_theta(i,j,k)    = denominator_r(i,j,k)*(dnum_dtheta    - rho(i,j,k)*dden_dtheta)
                 density_salinity(i,j,k) = denominator_r(i,j,k)*(dnum_dsalinity - rho(i,j,k)*dden_dsalinity)
                 density_press(i,j,k)    = denominator_r(i,j,k)*(dnum_dpress    - rho(i,j,k)*dden_dpress)*c2dbars

              enddo
           enddo
        enddo

    endif

    call diagnose_3d(Time, Grd, id_dpotrhodtheta, density_theta(:,:,:))
    call diagnose_3d(Time, Grd, id_dpotrhodsalt, density_salinity(:,:,:))
    call diagnose_3d(Time, Grd, id_dpotrhodpress, density_press(:,:,:))

  end subroutine density_derivs_potrho
! </SUBROUTINE> NAME="density_derivs_potrho"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_level">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to 
! temperature and with respect to salinity.  Hold pressure constant.  
!
! Pressure here is sea pressure = absolute press - press_standard
!
! </DESCRIPTION>
!
  subroutine density_derivs_level (rho, salinity, theta, press, Time, klevel, &
                                   density_theta, density_salinity, density_press)

    type(ocean_time_type),      intent(in)    :: Time
    real, dimension(isd:,jsd:), intent(in)    :: rho
    real, dimension(isd:,jsd:), intent(in)    :: salinity
    real, dimension(isd:,jsd:), intent(in)    :: theta
    real, dimension(isd:,jsd:), intent(in)    :: press
    integer,                    intent(in)    :: klevel 
    real, dimension(isd:,jsd:), intent(inout) :: density_theta
    real, dimension(isd:,jsd:), intent(inout) :: density_salinity
    real, dimension(isd:,jsd:), intent(inout) :: density_press
    real, dimension(isd:ied,jsd:jed)          :: pressure 

    integer :: i, j
    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta 
    real :: dnum_dsalinity, dden_dsalinity
    real :: dnum_dpress, dden_dpress

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_derivs_level): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:,:) = Grd%tmask(:,:,1)*potrho_press
    else 
       pressure(:,:) = press(:,:)
    endif 

    if(eos_linear) then

        do j=jsd,jed
           do i=isd,ied
              density_theta(i,j)    = -alpha_linear_eos 
              density_salinity(i,j) =  beta_linear_eos
              density_press(i,j)    =  0.0
           enddo
        enddo

    elseif(eos_preteos10) then 

        do j=jsd,jed
           do i=isd,ied


              t1  = theta(i,j)
              t2  = t1*t1
              s1  = salinity(i,j)
              sp5 = sqrt(s1) 

              p1   = pressure(i,j) - press_standard 
              p1t1 = p1*t1

              dnum_dtheta = a1 + t1*(two_a2 + three_a3*t1) &
                   + a5*s1                                 &
                   + p1t1*(two_a8 + two_a11*p1)    
              dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                   + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                   + p1*p1*(three_b11*t2 + b12*p1)

              dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
              dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

              dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
              dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

              density_theta(i,j)    = denominator_r(i,j,klevel)*(dnum_dtheta    - rho(i,j)*dden_dtheta)
              density_salinity(i,j) = denominator_r(i,j,klevel)*(dnum_dsalinity - rho(i,j)*dden_dsalinity)
              density_press(i,j)    = denominator_r(i,j,klevel)*(dnum_dpress    - rho(i,j)*dden_dpress)

           enddo
        enddo

    else ! eos_teos10 

        t1  = theta(i,j)
        t2  = t1*t1
        s1  = salinity(i,j)
        sp5 = sqrt(s1) 

        p1   = pressure(i,j) - press_standard 
        p1t1 = p1*t1

        dnum_dtheta = v02 + t1*(two_v03 + three_v04*t1)     &
             + s1*((v06+two_v07*t1)                         &
             +     sp5*(v09 + t1*(two_v10 + three_v11*t1))) &
             + p1*((v13 + two_v14*t1) + s1*v16              &
             +     p1*(v18+two_v19*t1))

        dden_dtheta = v22 + t1*(two_v23 +t1*(three_v24+four_v25*t1))      &
             + s1*(v27 + t1*(two_v28 + t1*(three_v29  + four_v30*t1))     &
             +     sp5*(v32+t1*(two_v33 + t1*(three_v34 + four_v35*t1)))) &
             + p1*((v38 + t1*(two_v39 + three_v40*t1))                    &
             +     s1*v42                                                 &
             +     p1*((v44 + two_v45*t1 + v46*s1)                        &
             +     p1*(v48)))

        dnum_dsalinity =                                     &
             + (v05 + t1*(v06 + v07*t1)                      &
             + 1.5*sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1)))) &
             + p1*((v15 + v16*t1)                            &
             +      p1*(v20))

        dden_dsalinity =                                                 &
             (v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + two_v36*s1 &
             + 1.5*sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))  &
             + p1*( v41 + v42*t1 + p1*t1*v46 )

        dnum_dpress = v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
             + 2.0*p1*(v17 + t1*(v18 + v19*t1) + v20*s1)
        dden_dpress = (v37 + t1*(v38 + t1*(v39 + v40*t1)) &
             + s1*(v41 + v42*t1)                          &
             + 2.*p1*(v43 + t1*(v44 + v45*t1 + v46*s1)    &
             + 1.5*p1*(v47 + v48*t1)))

        density_theta(i,j)    = denominator_r(i,j,klevel)*(dnum_dtheta    - rho(i,j)*dden_dtheta)
        density_salinity(i,j) = denominator_r(i,j,klevel)*(dnum_dsalinity - rho(i,j)*dden_dsalinity)
        density_press(i,j)    = denominator_r(i,j,klevel)*(dnum_dpress    - rho(i,j)*dden_dpress)


    endif   ! endif for eos choice


  end subroutine density_derivs_level
! </SUBROUTINE> NAME="density_derivs_level"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_point">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to
! temperature and with respect to salinity.  Do so here for a point. 
!
! Pressure here is 
!
! sea pressure = absolute pressure - press_standard  (dbars)
!
! </DESCRIPTION>
!
  subroutine density_derivs_point (rho, salinity, theta, press,&
                                   density_theta, density_salinity, density_press)

    real, intent(in)    :: rho
    real, intent(in)    :: salinity
    real, intent(in)    :: theta
    real, intent(in)    :: press
    real, intent(out)   :: density_theta
    real, intent(out)   :: density_salinity
    real, intent(out)   :: density_press

    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta 
    real :: dnum_dsalinity, dden_dsalinity
    real :: dnum_dpress, dden_dpress
    real :: denominator_point, denrecip
    real :: pressure 

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_density_mod (density_derivs_point): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure = potrho_press
    else
       pressure = press
    endif 


    if(eos_linear) then

        density_theta    = -alpha_linear_eos 
        density_salinity =  beta_linear_eos
        density_press    =  0.0

    elseif(eos_preteos10) then 

        t1  = theta
        t2  = t1*t1
        s1  = salinity
        sp5 = sqrt(s1) 

        p1   = pressure - press_standard 
        p1t1 = p1*t1

        denominator_point = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4 ))) &
             + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2))           &  
             + p1*(b10 + p1t1*(b11*t2 + b12*p1))

        denrecip = 1.0/(denominator_point + epsln)

        dnum_dtheta = a1 + t1*(2.0*a2 + 3.0*a3*t1) &
             + a5*s1                               &
             + 2.0*p1t1*(a8 + a11*p1)    
        dden_dtheta = b1 + t1*(2.0*b2 + t1*(3.0*b3 + 4.0*b4*t1)) &
             + s1*(b6 + t1*(3.0*b7*t1 + 2.0*b9*sp5))             &
             + p1*p1*(3.0*b11*t2 + b12*p1)

        dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
        dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

        dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
        dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

        density_theta    = denrecip*(dnum_dtheta    - rho*dden_dtheta)
        density_salinity = denrecip*(dnum_dsalinity - rho*dden_dsalinity)
        density_press    = denrecip*(dnum_dpress    - rho*dden_dpress)


    else  ! eos_teos10 

        t1  = theta
        t2  = t1*t1
        s1  = salinity
        sp5 = sqrt(s1) 

        p1   = pressure - press_standard 
        p1t1 = p1*t1

        denominator_point =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))&
             + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
             + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
             + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
             +     s1*(v41 + v42*t1)                                      &
             +     p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                   &
             +     p1*(v47 + v48*t1)))
        denrecip = 1./(denominator_point + epsln )

        dnum_dtheta = v02 + t1*(two_v03 + three_v04*t1)     &
             + s1*(v06+two_v07*t1                           &
             +     sp5*(v09 + t1*(two_v10 + three_v11*t1))) &
             + p1*((v13 + two_v14*t1) + s1*v16              &
             +     p1*(v18+two_v19*t1))

        dden_dtheta = v22 + t1*(two_v23 +t1*(three_v24+four_v25*t1))      &
             + s1*     (v27+t1*(two_v28 + t1*(three_v29 + four_v30*t1))   &
             +     sp5*(v32+t1*(two_v33 + t1*(three_v34 + four_v35*t1)))) &
             + p1*(v38 + t1*(two_v39 + three_v40*t1)+  s1*v42             &
             +     p1*(v44 + two_v45*t1 + v46*s1                          &
             +     p1*v48))

        dnum_dsalinity =                                    &
             + (v05 + t1*(v06 + v07*t1)                     &
             + 1.5*sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))&
             + p1*((v15 + v16*t1) &
             +      p1*(v20))

        dden_dsalinity =                                                 &
             (v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + two_v36*s1 &
             + 1.5*sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))  &
             + p1*( v41 + v42*t1 + p1*t1*v46 )

        dnum_dpress = v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
             + 2.0*p1*(v17 + t1*(v18 + v19*t1) + v20*s1)
        dden_dpress = (v37 + t1*(v38 + t1*(v39 + v40*t1)) &
             + s1*(v41 + v42*t1)                          &
             + 2.*p1*(v43 + t1*(v44 + v45*t1 + v46*s1)    &
             + 1.5*p1*(v47 + v48*t1)))
        density_theta    = denrecip*(dnum_dtheta    - rho*dden_dtheta)
        density_salinity = denrecip*(dnum_dsalinity - rho*dden_dsalinity)
        density_press    = denrecip*(dnum_dpress    - rho*dden_dpress)

    endif   ! endif for eos choice 


  end subroutine density_derivs_point
! </SUBROUTINE> NAME="density_derivs_point"



!#######################################################################
! <SUBROUTINE NAME="cabbeling_thermobaricity">
!
! <DESCRIPTION>
! Diagnostic sends for cabbeling and thermobaricity parameters.
!
! Pressure here is 
! sea pressure = absolute press - press_standard (dbars)
!
! </DESCRIPTION>
!
  subroutine cabbeling_thermobaricity(rho, salinity, theta, press, Time, &
             density_theta, density_salinity, density_press)

    real, dimension(isd:,jsd:,:), intent(in) :: rho
    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta
    real, dimension(isd:,jsd:,:), intent(in) :: press
    type(ocean_time_type),        intent(in) :: Time
    real, dimension(isd:,jsd:,:), intent(in) :: density_theta
    real, dimension(isd:,jsd:,:), intent(in) :: density_salinity
    real, dimension(isd:,jsd:,:), intent(in) :: density_press

    real, dimension(isd:ied,jsd:jed,nk) :: cabbeling_param
    real, dimension(isd:ied,jsd:jed,nk) :: thermobaric_param
    real, dimension(isd:ied,jsd:jed,nk) :: pressure 


    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (cabbeling_thermobaricity): module must be initialized')
    endif 

    if(density_equal_potrho) then 
       pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press
    else 
       pressure(:,:,:) = press(:,:,:)
    endif 

    call calc_cabbeling_thermobaricity(rho, salinity, theta, pressure, &
    density_theta, density_salinity, density_press, cabbeling_param, thermobaric_param)

    call diagnose_3d(Time, Grd, id_cabbeling, cabbeling_param(:,:,:))
    call diagnose_3d(Time, Grd, id_thermobaricity, thermobaric_param(:,:,:))


  end subroutine cabbeling_thermobaricity
! </SUBROUTINE> NAME="cabbeling_thermobaricity"


!#######################################################################
! <SUBROUTINE NAME="calc_cabbeling_thermobaricity">
!
! <DESCRIPTION>
! Compute cabbeling and thermobaricity parameters, as defined in 
! McDougall (1987).
!
! Pressure here is 
! sea pressure = absolute press - press_standard (dbars)
!
! </DESCRIPTION>
!
  subroutine calc_cabbeling_thermobaricity(rho, salinity, theta, press, &
             density_theta, density_salinity, density_press, cabbeling_param, thermobaric_param)

    real, dimension(isd:,jsd:,:), intent(in)    :: rho
    real, dimension(isd:,jsd:,:), intent(in)    :: salinity
    real, dimension(isd:,jsd:,:), intent(in)    :: theta
    real, dimension(isd:,jsd:,:), intent(in)    :: press
    real, dimension(isd:,jsd:,:), intent(in)    :: density_theta
    real, dimension(isd:,jsd:,:), intent(in)    :: density_salinity
    real, dimension(isd:,jsd:,:), intent(in)    :: density_press
    real, dimension(isd:,jsd:,:), intent(inout) :: cabbeling_param
    real, dimension(isd:,jsd:,:), intent(inout) :: thermobaric_param

    integer :: i, j, k
    real :: t1, t2, s1, sp5, spm5, p1, p2, p1t1

    real :: rhotheta_rhosalinity
    real :: rho_inv
    real :: d2rho_dtheta2, d2rho_dsalin2, d2rho_dsalin_dtheta
    real :: d2rho_dsalin_dpress, d2rho_dtheta_dpress

    real :: dden_dtheta, dden_dsalinity, dden_dpress

    real :: d2num_dtheta2, d2num_dsalin2, d2num_dtheta_dsalin
    real :: d2num_dsalin_dpress, d2num_dtheta_dpress
    real :: d2den_dtheta2, d2den_dsalin2, d2den_dsalin_dtheta 
    real :: d2den_dsalin_dpress, d2den_dtheta_dpress

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (calc_cabbeling_thermobaricity): module must be initialized')
    endif 

    ! zero initialization, and values for when eos_linear is used
    cabbeling_param(:,:,:)   = 0.0   ! cabbeling parameter 
    thermobaric_param(:,:,:) = 0.0   ! thermobaricity parameter 

    if(eos_preteos10) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1   = theta(i,j,k)
                t2   = t1*t1
                s1   = salinity(i,j,k)
                sp5  = sqrt(s1) 
                p1   = press(i,j,k) - press_standard 
                p2   = p1*p1
                p1t1 = p1*t1

                if(s1 > 0.0) then 
                   spm5 = 1.0/sp5
                else 
                   spm5 = 0.0
                endif
                rhotheta_rhosalinity = Grd%tmask(i,j,k)*density_theta(i,j,k)/(epsln+density_salinity(i,j,k))

                dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                     + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                     + p1*p1*(three_b11*t2 + b12*p1)
                dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)
                dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

                d2num_dtheta2 = two_a2 + six_a3*t1 + two_a8*p1 + two_a11*p2
                d2num_dsalin2 = two_a6
                d2num_dtheta_dsalin = a5
                d2num_dsalin_dpress = a9
                d2num_dtheta_dpress = two_a8*t1 + four_a11*p1*t1

                d2den_dtheta2 = two_b2 + six_b3*t1 + twelve_b4*t2  &
                              + six_b7*s1*t1 + two_b9*s1*sp5 + six_b11*p2*t1
                d2den_dsalin2 = .75*spm5*(b8 + b9*t2)
                d2den_dsalin_dtheta = b6 + three_b7*t2 + three_b9*sp5*t1
                d2den_dsalin_dpress = 0.0
                d2den_dtheta_dpress = six_b11*p1*t2 + three_b12*p2

                d2rho_dtheta2 = denominator_r(i,j,k) &
                 *(d2num_dtheta2-2.0*density_theta(i,j,k)*dden_dtheta-rho(i,j,k)*d2den_dtheta2)

                d2rho_dsalin2 = denominator_r(i,j,k) &
                 *(d2num_dsalin2-2.0*density_salinity(i,j,k)*dden_dsalinity-rho(i,j,k)*d2den_dsalin2)

                d2rho_dsalin_dtheta = denominator_r(i,j,k)  &
                 *( d2num_dtheta_dsalin                     &
                   -density_salinity(i,j,k)*dden_dtheta     &
                   -density_theta(i,j,k)*dden_dsalinity     &
                   -rho(i,j,k)*d2den_dsalin_dtheta)

                d2rho_dsalin_dpress = denominator_r(i,j,k)  &
                 *( d2den_dsalin_dpress                     &
                   -density_press(i,j,k)*dden_dsalinity     &
                   -density_salinity(i,j,k)*dden_dpress     &
                   -rho(i,j,k)*d2den_dsalin_dpress)

                d2rho_dtheta_dpress = denominator_r(i,j,k)  &
                 *( d2den_dtheta_dpress                     &
                   -density_press(i,j,k)*dden_dtheta        &
                   -density_theta(i,j,k)*dden_dpress        &
                   -rho(i,j,k)*d2den_dtheta_dpress)
 
                rho_inv = Grd%tmask(i,j,k)/(rho(i,j,k) + epsln)

                cabbeling_param(i,j,k) = -rho_inv*               &
                  ( d2rho_dtheta2                                &
                   -2.0*d2rho_dsalin_dtheta*rhotheta_rhosalinity &
                   +d2rho_dsalin2*rhotheta_rhosalinity**2 )
                   
                thermobaric_param(i,j,k) = -rho_inv*(d2rho_dtheta_dpress-d2rho_dsalin_dpress*rhotheta_rhosalinity)

             enddo
          enddo
       enddo

    elseif(eos_teos10) then 

      do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1   = theta(i,j,k)
                t2   = t1*t1
                s1   = salinity(i,j,k)
                sp5  = sqrt(s1)
                p1   = press(i,j,k) - press_standard
                p2   = p1*p1
                p1t1 = p1*t1

                if(s1 > 0.0) then
                   spm5 = 1.0/sp5
                else
                   spm5 = 0.0
                endif
                rhotheta_rhosalinity = Grd%tmask(i,j,k)*density_theta(i,j,k)/(epsln+density_salinity(i,j,k))

                dden_dtheta = v22 + t1*(two_v23 +t1*(three_v24+four_v25*t1))      &
                     + s1*(v27 + t1*(two_v28 + t1*(three_v29  + four_v30*t1))     &
                     +     sp5*(v32+t1*(two_v33 + t1*(three_v34 + four_v35*t1)))) &
                     + p1*((v38 + t1*(two_v39 + three_v40*t1))                    &
                     +     s1*v42                                                 &
                     +     p1*((v44 + two_v45*t1 + v46*s1)                        &
                     +     p1*(v48)))


                dden_dsalinity =                                                    &
                        (v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + two_v36*s1 &
                       + 1.5*sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))   &
                       + p1*( v41 + v42*t1 + p1*t1*v46 )

               dden_dpress = (v37 + t1*(v38 + t1*(v39 + v40*t1))                    &
                       + s1*(v41 + v42*t1)                                          &
                       + 2.*p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                    &
                       + 1.5*p1*(v47 + v48*t1)))

                d2num_dtheta2 = (two_v03 + six_v04*t1)    &
                     + s1*(two_v07                        &
                     +     sp5*( (two_v10 + six_v11*t1))) &
                     + p1*(two_v14                        &
                     +     p1*two_v19)
                d2num_dsalin2 = 0.75*(v08 + t1*(v09 + t1*(v10 + v11*t1)))*spm5
                d2num_dtheta_dsalin = v06 + two_v07*t1              &
                      + 1.5*sp5*(v09 + t1*(two_v10 + three_v11*t1)) &
                      + p1*v16 
                d2num_dsalin_dpress = (v15 + v16*t1) &
                      + 2.0*p1*v20
                d2num_dtheta_dpress = (v13 + two_v14*t1) + s1*v16 &
                      + 2.0*p1*(v18 + two_v19*t1 )



                d2den_dtheta2 = (two_v23 +t1*(six_v24 + twelve_v25*t1))    &
                     + s1*(two_v28 + t1*(six_v29  + twelve_v30*t1)         &
                     +     sp5*((two_v33 + t1*(six_v34 + twelve_v35*t1)))) &
                     + p1*(((two_v39 + six_v40*t1))                        &
                     +     p1*two_v45       )
                d2den_dsalin2 = two_v36  &
                       + 0.75*spm5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1))))
                d2den_dsalin_dtheta = ((v27 + t1*(two_v28 + t1*(three_v29 + four_v30*t1))) &
                       + 1.5*sp5*((v32 + t1*(two_v33 + t1*(three_v34 + four_v35*t1)))))    &
                       + p1*( v42 + p1*v46 )
                d2den_dsalin_dpress = v41 + v42*t1 + p1*t1*two_v46 
                d2den_dtheta_dpress = (v38 + t1*(two_v39 + three_v40*t1)) &
                     +     s1*v42                                         &
                     +     2.*p1*((v44 + two_v45*t1 + v46*s1)             &
                     +             1.5*p1*(v48))


                d2rho_dtheta2 = denominator_r(i,j,k) &
                 *(d2num_dtheta2-2.0*density_theta(i,j,k)*dden_dtheta-rho(i,j,k)*d2den_dtheta2)

                d2rho_dsalin2 = denominator_r(i,j,k) &
                 *(d2num_dsalin2-2.0*density_salinity(i,j,k)*dden_dsalinity-rho(i,j,k)*d2den_dsalin2)

                d2rho_dsalin_dtheta = denominator_r(i,j,k)  &
                 *( d2num_dtheta_dsalin                     &
                   -density_salinity(i,j,k)*dden_dtheta     &
                   -density_theta(i,j,k)*dden_dsalinity     &
                   -rho(i,j,k)*d2den_dsalin_dtheta)

                d2rho_dsalin_dpress = denominator_r(i,j,k)  &
                 *( d2den_dsalin_dpress                     &
                   -density_press(i,j,k)*dden_dsalinity     &
                   -density_salinity(i,j,k)*dden_dpress     &
                   -rho(i,j,k)*d2den_dsalin_dpress)

                d2rho_dtheta_dpress = denominator_r(i,j,k)  &
                 *( d2den_dtheta_dpress                     &
                   -density_press(i,j,k)*dden_dtheta        &
                   -density_theta(i,j,k)*dden_dpress        &
                   -rho(i,j,k)*d2den_dtheta_dpress)

                rho_inv = Grd%tmask(i,j,k)/(rho(i,j,k) + epsln)

                cabbeling_param(i,j,k) = -rho_inv*               &
                  ( d2rho_dtheta2                                &
                   -2.0*d2rho_dsalin_dtheta*rhotheta_rhosalinity &
                   +d2rho_dsalin2*rhotheta_rhosalinity**2 )

                thermobaric_param(i,j,k) = -rho_inv*(d2rho_dtheta_dpress-d2rho_dsalin_dpress*rhotheta_rhosalinity)

             enddo
          enddo
       enddo

    endif ! endif for eos choice 



  end subroutine calc_cabbeling_thermobaricity
! </SUBROUTINE> NAME="calc_cabbeling_thermobaricity"



!#######################################################################
! <FUNCTION NAME="density_delta_z">
!
! <DESCRIPTION>
! rho(k)-rho(k+1) for all i,j with both temperatures referenced to the 
! deeper pressure depth.  
!
! Of use for KPP scheme. 
! </DESCRIPTION>
!
function density_delta_z (rho_initial, salinity, theta, press)

  real, intent(in), dimension(isd:,jsd:,:) :: rho_initial
  real, intent(in), dimension(isd:,jsd:,:) :: salinity
  real, intent(in), dimension(isd:,jsd:,:) :: theta
  real, intent(in), dimension(isd:,jsd:,:) :: press

  real, dimension(isd:ied,jsd:jed,nk) :: density_delta_z
  real, dimension(isd:ied,jsd:jed,nk) :: pressure 
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_density_mod (density_delta_z): module must be initialized')
  endif 

  if(density_equal_potrho) then 
     pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press
  else 
     pressure(:,:,:) = press(:,:,:)
  endif 

  do k=1,nk-1
    wrk1(:,:,k) = pressure(:,:,k+1)
  enddo
  wrk1(:,:,nk) = pressure(:,:,nk)

  wrk2(:,:,:) = density(salinity(:,:,:), theta(:,:,:), wrk1(:,:,:))

  do k=1,nk-1
    density_delta_z(:,:,k) = (wrk2(:,:,k) - rho_initial(:,:,k+1))*Grd%tmask(:,:,k+1)
  enddo
  density_delta_z(:,:,nk) = (wrk2(:,:,nk) - rho_initial(:,:,nk))*Grd%tmask(:,:,nk)

end function density_delta_z
! </FUNCTION> NAME="density_delta_z"



!#######################################################################
! <FUNCTION NAME="density_delta_sfc">
!
! <DESCRIPTION>
! rho(1)-rho(k+1) for all i,j. 
!
! Of use for KPP scheme. 
! </DESCRIPTION>
!
function density_delta_sfc (rho_initial, salinity, theta, press)

  real, intent(in), dimension(isd:,jsd:,:) :: rho_initial
  real, intent(in), dimension(isd:,jsd:,:) :: salinity
  real, intent(in), dimension(isd:,jsd:,:) :: theta
  real, intent(in), dimension(isd:,jsd:,:) :: press

  real, dimension(isd:ied,jsd:jed,nk) :: density_delta_sfc
  real, dimension(isd:ied,jsd:jed,nk) :: pressure
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_density_mod (density_delta_sfc): module must be initialized')
  endif 

  if(density_equal_potrho) then 
     pressure(:,:,:) = Grd%tmask(:,:,:)*potrho_press
  else 
     pressure(:,:,:) = press(:,:,:)
  endif 

  do k=1,nk-1
    wrk1(:,:,k) = pressure(:,:,k+1)
  enddo
  wrk1(:,:,nk) = pressure(:,:,nk)

  wrk2(:,:,:) = density_sfc(salinity, theta, wrk1)

  do k=1,nk-1
    density_delta_sfc(:,:,k) = (wrk2(:,:,k) - rho_initial(:,:,k+1))*Grd%tmask(:,:,k+1)
  enddo
  density_delta_sfc(:,:,nk) = (wrk2(:,:,nk) - rho_initial(:,:,nk))*Grd%tmask(:,:,nk)

end function density_delta_sfc
! </FUNCTION> NAME="density_delta_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_end">
!
! <DESCRIPTION>
!
! Write density and pressure_at_depth to a restart.
!
! </DESCRIPTION>
subroutine ocean_density_end(Time, Dens, use_blobs)

  type(ocean_time_type),    intent(in)  :: Time
  type(ocean_density_type), intent(in)  :: Dens
  logical,                  intent(in)  :: use_blobs
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error in ocean_density_mod (ocean_density_end): module must be initialized')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') '==>Warning from ocean_density_mod (ocean_density_end): NO restart written.'
    call mpp_error(WARNING,'==>Warning from ocean_density_mod (ocean_density_end): NO restart written.')
    return
  endif 

  call ocean_density_restart(Time, Dens)

  write(stdoutunit,'(/a)') ' From ocean_density_mod: ending density chksums'
  call write_timestamp(Time%model_time)
  call ocean_density_chksum(Time, Dens, use_blobs)

  module_is_initialized = .FALSE.

end subroutine ocean_density_end
! </SUBROUTINE> NAME="ocean_density_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_density_restart(Time, Dens, time_stamp)
  type(ocean_time_type),    intent(in)           :: Time
  type(ocean_density_type), intent(in)           :: Dens
  character(len=*),         intent(in), optional :: time_stamp
   integer :: taup1

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error in ocean_density_mod (ocean_density_end): module must be initialized')
  endif 

  taup1  = Time%taup1

  call reset_field_pointer(Den_restart, id_restart_rho, Dens%rho(:,:,:,taup1) )
  call reset_field_pointer(Den_restart, id_restart_rho_s, Dens%rho_salinity(:,:,:,taup1) )

  call save_restart(Den_restart, time_stamp)

end subroutine ocean_density_restart
! </SUBROUTINE> NAME="ocean_density_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_chksum">
!
! <DESCRIPTION>
! Compute checksums for density. 
! </DESCRIPTION>
subroutine ocean_density_chksum(Time, Dens, use_blobs)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens
  logical,                  intent(in) :: use_blobs

  integer :: taup1

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_density_mod (ocean_density_chksum): module not yet initialized ')
  endif 

  taup1 = Time%taup1 

  call write_chksum_3d('rho(taup1)', Dens%rho(COMP,:,taup1)*Grd%tmask(COMP,:))
  if (use_blobs) then
     call write_chksum_3d('rhoT(taup1)', Dens%rhoT(COMP,:)*Grd%tmask(COMP,:))
  endif
  call write_chksum_3d('pressure_at_depth', Dens%pressure_at_depth(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('denominator_r', denominator_r(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('drhodT', Dens%drhodT(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('drhodS', Dens%drhodS(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('drhodz_zt', Dens%drhodz_zt(COMP,:)*Grd%tmask(COMP,:))

end subroutine ocean_density_chksum
! </SUBROUTINE>  NAME="ocean_density_chksum"


!#######################################################################
! <SUBROUTINE NAME="compute_buoyfreq">
! <DESCRIPTION>
!
! Diagnose the buoyancy frequency, both at T-points and at 
! vertical interfaces of T-cells. 
!
! Author: Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine compute_buoyfreq(Time, Thickness, salinity, theta, Dens, use_blobs)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity
  real, dimension(isd:,jsd:,:), intent(in)    :: theta
  type(ocean_density_type),     intent(inout) :: Dens
  logical,                      intent(in)    :: use_blobs

  integer   :: i,j,k,kp1,kbot,m
  integer   :: tau, taup1
  real      :: drhodz_prev, tmpdrhodz

  if (.not. module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_density (compute_buoyfreq): module needs initialization')
  endif

  tau   = Time%tau
  taup1 = Time%taup1
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  if (use_blobs) then
     ! vertical derivative of temp and salt
     do k=1,nk-1
        kp1 = k+1
        do j=jsd,jed
           do i=isd,ied
              wrk1(i,j,k) = (theta(i,j,k)-theta(i,j,kp1))/Thickness%dzwtT(i,j,k)
              wrk2(i,j,k) = (salinity(i,j,k)-salinity(i,j,kp1))/Thickness%dzwtT(i,j,k)
           enddo
        enddo
     enddo
  else
     ! vertical derivative of temp and salt
     do k=1,nk-1
        kp1 = k+1
        do j=jsd,jed
           do i=isd,ied
              wrk1(i,j,k) = (theta(i,j,k)-theta(i,j,kp1))/Thickness%dzwt(i,j,k)
              wrk2(i,j,k) = (salinity(i,j,k)-salinity(i,j,kp1))/Thickness%dzwt(i,j,k)
           enddo
        enddo
     enddo
  endif

  ! vertical derivative of locally referenced potential density.
  ! vanishes at the bottom of a column.
  do k=1,nk-1
     kp1=k+1
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = Dens%drhodT(i,j,k)*wrk1(i,j,k) + Dens%drhodS(i,j,k)*wrk2(i,j,k)  
           if(abs(wrk3(i,j,k)) < epsln_drhodz) then 
              wrk3(i,j,k) = sign(1.0,wrk3(i,j,k))*epsln_drhodz
           endif 
           wrk3(i,j,k) = Grd%tmask(i,j,kp1)*wrk3(i,j,k)
        enddo
     enddo
  enddo

  ! vertically smooth drhodz; otherwise can get noisy results 
  if(buoyfreq_smooth_vert) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               drhodz_prev = onefourth*wrk3(i,j,1) 
               kbot=Grd%kmt(i,j)
               if(kbot>3) then
                   do k=2,kbot-2
                      tmpdrhodz   = wrk3(i,j,k)
                      wrk3(i,j,k) = drhodz_prev + onehalf*wrk3(i,j,k) + onefourth*wrk3(i,j,k+1)
                      drhodz_prev = onefourth*tmpdrhodz
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! drhodz at bottom of T-cell
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Dens%drhodz_wt(i,j,k) = wrk3(i,j,k)
        enddo
     enddo
  enddo
  if(id_buoyfreq2_wt > 0) then 
     call diagnose_3d(Time, Grd, id_buoyfreq2_wt, -grav*rho0r*wrk3(:,:,:))
  endif
  call diagnose_3d(Time, Grd, id_drhodz_wt, wrk3(:,:,:))

  ! drhodz and squared buoyancy frequency at T-cell point 
  wrk4(:,:,:) = 0.0
  k=1
  do j=jsd,jed
     do i=isd,ied
        Dens%drhodz_zt(i,j,k)   = Dens%drhodz_wt(i,j,k)
        wrk4(i,j,k)             = Dens%drhodz_wt(i,j,k)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           wrk4(i,j,k) = wrk3(i,j,k-1)*(Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k))   &
                        +wrk3(i,j,k)  *(Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)) 
           wrk4(i,j,k) = wrk4(i,j,k)/(epsln+Thickness%dzt(i,j,k))
           Dens%drhodz_zt(i,j,k)   = wrk4(i,j,k)
        enddo
     enddo
  enddo

  if(drhodz_diag_stable) then

      ! enforce stable stratification 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               if(abs(Dens%drhodz_zt(i,j,k)) < epsln_drhodz_diag .or. Dens%drhodz_zt(i,j,k) > 0.0) then 
                   Dens%drhodz_diag(i,j,k) = -epsln_drhodz_diag
               else 
                   Dens%drhodz_diag(i,j,k) = Dens%drhodz_zt(i,j,k)
               endif
            enddo
         enddo
      enddo
  else 

      ! allow for unstable stratification 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               if(abs(Dens%drhodz_zt(i,j,k)) < epsln_drhodz_diag) then 
                   Dens%drhodz_diag(i,j,k) = sign(1.0,Dens%drhodz_zt(i,j,k))*epsln_drhodz_diag
               else 
                   Dens%drhodz_diag(i,j,k) = Dens%drhodz_zt(i,j,k)
               endif
            enddo
         enddo
      enddo
  endif


  if(id_buoyfreq2_zt > 0) then 
     call diagnose_3d(Time, Grd, id_buoyfreq2_zt, -grav*rho0r*wrk4(:,:,:))
  endif 
  call diagnose_3d(Time, Grd, id_drhodz_zt, wrk4(:,:,:))
  call diagnose_3d(Time, Grd, id_drhodz_diag, Dens%drhodz_diag(:,:,:))

  ! dTdz and dSdz  at T-cell point 
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  k=1
  do j=jsd,jed
     do i=isd,ied
        Dens%dTdz_zt(i,j,k) = wrk1(i,j,k)
        Dens%dSdz_zt(i,j,k) = wrk2(i,j,k)
     enddo
  enddo

  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = wrk1(i,j,k-1)*(Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k))   &
                        +wrk1(i,j,k)  *(Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)) 
           wrk4(i,j,k) = wrk2(i,j,k-1)*(Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k))   &
                        +wrk2(i,j,k)  *(Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)) 
           wrk3(i,j,k) = wrk3(i,j,k)/(epsln+Thickness%dzt(i,j,k))
           wrk4(i,j,k) = wrk4(i,j,k)/(epsln+Thickness%dzt(i,j,k))
           Dens%dTdz_zt(i,j,k) = wrk3(i,j,k)
           Dens%dSdz_zt(i,j,k) = wrk4(i,j,k)
        enddo
     enddo
  enddo


end subroutine compute_buoyfreq
! </SUBROUTINE>  NAME="compute_buoyfreq"


!#######################################################################
! <FUNCTION NAME="buoyfreq2">
! <DESCRIPTION>
!
! Diagnose the square of the buoyancy frequency at the bottom of 
! T-cells, NOT at T-points.
! The algorithm follows that used by the private function
! compute_buoyfreq in the density module.
!
! We take the square of the buoyancy frequency as is, we do not
! smooth or force it to be positive.  This allows us to search for
! instabilities.
!
! Authors: m.bates
! </DESCRIPTION>
!
function buoyfreq2(Time, Thickness, Dens, salinity, theta, rho, use_blobs)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: salinity
  real, dimension(isd:,jsd:,:), intent(in) :: theta
  real, dimension(isd:,jsd:,:), intent(in) :: rho
  logical,                      intent(in) :: use_blobs

  real, dimension(isd:ied,jsd:jed,nk) :: buoyfreq2
  real, dimension(isd:ied,jsd:jed,nk) :: drhodT
  real, dimension(isd:ied,jsd:jed,nk) :: drhodS

  integer   :: i,j,k,kp1
  integer   :: taup1
  real      :: tmp, grav_rho0r

  if (.not. module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_density (buoyfreq2): module needs initialization')
  endif

  taup1 = Time%taup1

  grav_rho0r = grav*rho0r

  wrk4(:,:,:)      = 0.0
  buoyfreq2(:,:,:) = 0.0

  !get the partial derivatives drho/dT and drho/dS
  call density_derivs( rho(:,:,:), salinity(:,:,:), theta(:,:,:), &
                      Dens%pressure_at_depth(:,:,:),              &
                      Time, drhodT(:,:,:), drhodS(:,:,:),         &
                      wrk4(:,:,:))

  wrk3(:,:,:) = 0.0
  ! We don't need drho/dp and so can overwrite it.
  wrk4(:,:,:) = 0.0

  ! vertical derivative of temp and salt
  if (use_blobs) then
     do k=1,nk-1
        kp1 = k+1
        do j=jsd,jed
           do i=isd,ied
              tmp = Grd%tmask(i,j,kp1)/Thickness%dzwtT(i,j,k)
              wrk3(i,j,k) = tmp*(theta(i,j,k)    - theta(i,j,kp1)   )
              wrk4(i,j,k) = tmp*(salinity(i,j,k) - salinity(i,j,kp1))
           enddo
        enddo
     enddo
  else
     do k=1,nk-1
        kp1 = k+1
        do j=jsd,jed
           do i=isd,ied
              tmp = Grd%tmask(i,j,kp1)/Thickness%dzwt(i,j,k)
              wrk3(i,j,k) = tmp*(theta(i,j,k)    - theta(i,j,kp1)   )
              wrk4(i,j,k) = tmp*(salinity(i,j,k) - salinity(i,j,kp1))
           enddo
        enddo
     enddo
  endif

  ! vertical derivative of locally referenced potential density
  ! at the bottom of a cell.
  ! N**2==0 at the bottom of a column.
  do k=1,nk-1
     do j=jsd,jed
        do i=isd,ied
           buoyfreq2(i,j,k) = -grav_rho0r*( drhodT(i,j,k)*wrk3(i,j,k) + drhodS(i,j,k)*wrk4(i,j,k) )
        enddo
     enddo
  enddo

end function buoyfreq2
! </FUNCTION>  NAME="buoyfreq2"



!#######################################################################
! <SUBROUTINE NAME="compute_drhodxy">
! <DESCRIPTION>
!
! Diagnose horiz derivatives of locally referenced potential density.
! For use in diagnosing balanced Ertel PV.  
!
! For simplicity, do no spatial averaging. Hence, the placement
! of the derivative is at the cell face, not cell center. But to 
! simplify matters, we consider the derivs to be at the cell center,
! which simplifies how we take |grad b| for the PV calculation, meaning
! there will be no averaging operations. This assumption may need to
! be revisited.  
!
! Author: Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine compute_drhodxy(Time, salinity, theta, Dens)

  type(ocean_time_type),        intent(in)    :: Time
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity
  real, dimension(isd:,jsd:,:), intent(in)    :: theta
  type(ocean_density_type),     intent(inout) :: Dens

  integer   :: i,j,k
  real      :: tmp1,tmp2

  if (.not. module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_density (compute_drhodxy): module needs initialization')
  endif

  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0

  ! x-derivative 
  do k=1,nk
     do j=jsd,jed
        do i=isd,iec
           tmp1 = (theta(i+1,j,k)   -theta(i,j,k)   )*Dens%drhodT(i,j,k)
           tmp2 = (salinity(i+1,j,k)-salinity(i,j,k))*Dens%drhodS(i,j,k)
           wrk1(i,j,k)           = Grd%tmask(i+1,j,k)*(tmp1+tmp2)*Grd%dxter(i,j)
           Dens%drhodx_zt(i,j,k) = wrk1(i,j,k)
        enddo
     enddo
  enddo
  
  ! y-derivative 
  do k=1,nk
     do j=jsd,jec
        do i=isd,ied
           tmp1 = (theta(i,j+1,k)   -theta(i,j,k)   )*Dens%drhodT(i,j,k)
           tmp2 = (salinity(i,j+1,k)-salinity(i,j,k))*Dens%drhodS(i,j,k)
           wrk2(i,j,k)           = Grd%tmask(i,j+1,k)*(tmp1+tmp2)*Grd%dytnr(i,j)
           Dens%drhody_zt(i,j,k) = wrk2(i,j,k)
        enddo
     enddo
  enddo
  
  if(id_drhodx_zt > 0) then 
     call diagnose_3d(Time, Grd, id_drhodx_zt, wrk1(:,:,:))
  endif
  if(id_drhody_zt > 0) then 
     call diagnose_3d(Time, Grd, id_drhody_zt, wrk2(:,:,:))
  endif



end subroutine compute_drhodxy
! </SUBROUTINE>  NAME="compute_drhodxy"



end module ocean_density_mod





