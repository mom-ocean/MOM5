#include <fms_platform.h>

! ----------------------------------------------------------------
!                   GNU General Public License                        
! This file is a part of MOM.                                                                 
!                                                                      
! MOM is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! 
! <CONTACT EMAIL="Eric.Galbraith@mcgill.ca"> Eric Galbraith
! </CONTACT>
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Anand Gnanandesikan
! </CONTACT>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Rick Slater
! </REVIEWER>
!
!<OVERVIEW>
!                    IDEALIZED OCEAN BIOGEOCHEMISTRY                        *
!                                                                           *
! An idealized biogeochemical cycling module, in which the ecosystem and    * 
! particles are treated implicitly. This is intended as a first-order       *
! complexity complement to more intricate ecosystem/biogeochemical models   *
! such as BLING, or the GFDL TOPAZ model.                                   *
!
! iBGC is not intended to simulate marine biogeochemistry in the most       *
! accurate sense possible - that is left to the more complex models.        *
! Instead, it is intended as a tool for basic research, being easily        *
! modified, simply understood, and readily used in idealized configurations.*
!
! Note that the model parameters are not particularly well tuned for any    *
! ocean model, and the user should feel free to alter anything as they see  * 
! fit.                                                                      *
! 
!</OVERVIEW>
!
!<DESCRIPTION>
!
!* Inorganic nutrients are incorporated into organic matter in the presence *
!* of light. Once in the organic pool, the nutrient elements are recycled   *
!* between organisms and labile reservoirs as long as light is present in   *
!* sufficient abundance to fuel continued growth. In the absence of light,  *
!* inorganic nutrients gradually accumulate once again. The ideal nutrient  *
!* tracers are intended to capture the fundamental dynamics of these        *
!* light-dependent inorganic-organic transitions, but in a very simple,     *
!* transparent way that depends only on the physical circulation and the    *
!* availability of PAR. Thus, they ignore the complexities of ecosystems,   *
!* in terms of uptake dynamics and remineralization pathways. Although this *
!* will certainly decrease their ability to reproduce the observed          *
!* distributions of nutrients in the ocean, they are computationally        *
!* inexpensive and easy to interpret, features which will hopefully give    *
!* them some utility.                                                       *

!* The simplest ideal nutrient tracer, ideal_n, is non-conservative,        *
!* designed to quickly approach  equilibrium (order 10 years), avoiding the *
!* need for long integrations. The slightly less straightforward,           *
!* conservative tracers are remineralized through the water column and      *
!* therefore approach steady state on much longer, multi-centennial         *
!* timescales.                                                              *
!*
!* The core behaviour of the ideal nutrients are given by the exponential,  *
!* exp(-IRR/IRRk). IRR is the available light (average SW radiation, in     *
!* W/m2, within the grid cell) and IRRk (in W/m2) determines the light      *
!* level at which the phytoplankton growth rate approaches saturation.      *
!* Thus, this term is 1 when IRR is 0, and approaches 0 as IRR increases. A *
!* larger value for IRRk causes the exponential to approach 0 more slowly,  *
!* ie. to saturate at a higher light level.                                 *
!*                                                                          
!* Uptake is determined by the product of the 1 minus the exponential term, *
!* a maximum uptake velocity (Vmax, in s-1), and limitation caused by the   *
!* concentration of the ideal nutrient itself, N (in mol kg-1):             *
!*                                                                          *
!* -Vmax * (1 - exp( -IRR / IRRk )) * (N / N + k)                           *
!*                                                                          *
!* The last term causes nutrient uptake to slow as the nutrient             *
!* concentration approaches 0, according to Michaelis-Menten dynamics.      *
!*                                                                          *
!* In ideal_n, regeneration is simply a constant rate, R.                   *
!*
!* Starting from this simple foundation, biogeochemical functionality is    *
!* added by creating a macronutrient (iPO4) with identical uptake to        *
!* ideal_n, but for which mass is conserved within the ocean. This is       *
!* achieved by separating the uptake into dissolved and particulate         *
!* organic matter, the latter of which is instantaneously remineralized     *
!* throughout the water column below according to a variant of the OCMPIP2  *
!* protocol. This calculates the remineralization rate at each level as a   *
!* function of the temperature, base remineralization rate, and sinking     *
!* rate (itself a function of depth). Any flux reaching the bottom box is   *
!* instantly remineralized there to conserve mass within the ocean.         *
!* Meanwhile, Dissolved Organic Phosphorus (iDOP) is transported within the *
!* ocean, decaying to PO4 at a temperature-dependent rate.                  *
!*
!* Given the apparent importance of iron (Fe) as a limiting micronutrient   *
!* in the oceans, it seemed interesting to try making another PO4 tracer,   * 
!* limited by Fe, in addition to the basic limitations given above: this is *
!* called iPO4f. The iron input is highly idealized, and the scavenging     *
!* function very simple, but iron-limitation develops and modulates         *
!* growth rate through iron-light colimitation in a way that may not be     *
!* totally unlike the real ocean. The Fe-limited PO4 cycle is completely    *
!* independent of the non-Fe-limited PO4, and they can be run on their own  *
!* or in parallel for comparison.
!*
!* A full biogeochemical simulation can then be made by selecting one of    *
!* the iPO4s as the central bgc variable (default is iPO4). Biological      * 
!* consumption and production of other quantities, including a number of    *
!* gases and idealized gas tracers, are then determined from the generally  *
!* reliable Redfield stoichiometries.                                       *
!*
!* Isotopic fractionations are simulated by multiplying the uptake rate of  *
!* the chemical species by the isotopic rate ratio, alpha. The tracer       *
!* concentrations of the heavy isotope species are actually                 *
!* equal to the true concentration divided by the isotopic ratio of the     *
!* standard, e.g. 15NO3 = [15NO3] * [14Nair] / [15Nair]. The true value in  *
!* permil is given by (15NO3 / NO3 - 1.) * 1000.                            *
!*
!* Dissolved gases are handled according to the OCMIP2 protocol. The air-   *
!* sea exchange code was taken from the ocmip2_abiotic and ocmip2_biotic    *
!* modules, with negligible modification. Note that the 14C implementation  *
!* differs from that of abiotic, in that biological uptake and remineraliz- *
!* ation of carbon impact the distribution of 14C, analagous to 15N.        *
!*
!* Chlorophyll is estimated diagnostically from the net bgc iPO4 (or iPO4f) *
!* uptake rate, with an adjustment for photoadaptation.
!
!
!    A general note on nomenclature
! Variables starting with 'j' are source/sink terms; jprod_x is the 
! biological production term for the quantity 'x'; jremin_x is the 
! remineralization rate of the quantity 'x'. Variables starting with 'fp' 
! are sinking particulate flux terms. Most of the 'i's stand for idealized,
! as in iBGC, and help to avoid duplication of arrays and tracer names used
! by TOPAZ and BLING, so that ibgc can be run in parallel with these other 
! models. Most of the 'p's stand for particulate, though a lot of them 
! stand for phosphorus. So, for example, jprod_pip is the production rate 
! of particulate idealized phosphorus through the uptake of iPO4.
!
! Tracer units are generally in mol kg-1, with some exceptions (suntan, 
! ideal_n, Fe).
!
! A general note on parameter values:
! Pretty much all of the parameter values were selected ad hoc, from very
! loose principles. They were not particularly well 'tuned' to match any
! given circulation model, though they gave reasonable results with the 
! 3-degree ocean model, om1p7. As such, the user should feel free to make 
! changes to the default values given here.
! Exceptions to this are: gas exchange parameters, radiocarbon decay, 
! nutrient stoichiometry, fractionation factors.

!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! The basic nutrient uptake equations will be described in an upcoming paper
! by Galbraith et al. (in prep.). The remainder of the biogeochemistry 
! module may be documented more officially elsewhere. All are welcome to use
! ibgc model output for the purposes of publication; however, please contact
! Eric Galbraith (Eric.Galbraith@mcgill.ca) for the appropriate reference.
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_ibgc_nml">
!
!  <DATA NAME="do_ideal" TYPE="logical">
!  If true, then do ideal_n and suntan. This does not require any other
! part of the model.
!  </DATA> 
!
!  <DATA NAME="do_po4" TYPE="logical">
!  If true, then do the non-Fe-limited P cycle, iPO4 and iDOP. If either 
! this or do_po4f is true, PO4_pre and chl will be calculated as well. If 
! both do_po4 and do_po4f are true, po4 will be the master variable,
! unless bgc_felim is true.
!  </DATA> 
!
!  <DATA NAME="do_po4f" TYPE="logical">
!  If true, then do the Fe-limited P cycle, with iFe, iPO4f and DOP. If  
! either this or do_po4 is true, PO4_pre and chl will be calculated as well.
! If both do_po4 and do_po4f are true, po4 will be the master variable,
! unless do_bgc_felim is true.
!  </DATA> 
!
!  <DATA NAME="do_bgc_felim" TYPE="logical">
!  If true, then use PO4f as the master variable for biogeochemical 
! calculations (gases, PO4_pre, chl, isotopes). Requires that do_PO4f be true.
!  </DATA> 
!
!  <DATA NAME="do_gasses" TYPE="logical">
!  If true, then do the gases Dissolved Inorganic Carbon (iDIC) and oxygen 
! (iO2). Requires that do_po4 and/or do_po4f be true.
!  </DATA> 
!
!  <DATA NAME="do_carbon_comp" TYPE="logical">
!  If true, then do the dissolved inorganic carbon component tracers,  
! Saturation DIC (iDIC_sat) and preformed DIC (iDIC_pre). Requires
! that do_po4 and/or do_po4f be true, and that do_gasses be true.
!  </DATA> 
!
!  <DATA NAME="do_radiocarbon" TYPE="logical">
!  If true, then do the radiocarbon tracers (iDI14C and iDO14C). Requires
! that do_po4 and/or do_po4f be true, and that do_gasses be true. 
!  </DATA> 
!
!  <DATA NAME="do_isio4" TYPE="logical">
!  If true, then do silica cycle and silicon isotopes, iSiO4 and i30SiO4.
!  </DATA> 
!
!  <DATA NAME="do_no3_iso" TYPE="logical">
!  If true, then do NO3 isotopes, i15NO3, iN18O3 and iDO15N. Requires that
! do_po4 and/or do_po4f be true.
!  </DATA> 
!
!</NAMELIST>
!
!------------------------------------------------------------------
!
!       Module ocean_ibgc_mod
!
!------------------------------------------------------------------

module  ocean_ibgc_mod

use diag_manager_mod,         only: register_diag_field, diag_axis_init
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value, fm_get_index
use mpp_mod,                  only: input_nml_file, stdout, stdlog, mpp_error, mpp_sum, mpp_chksum, FATAL, NOTE
use fms_mod,                  only: field_exist, file_exist
use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
use fms_io_mod,               only: register_restart_field, save_restart, restore_state
use fms_io_mod,               only: restart_file_type
use time_manager_mod,         only: get_date, time_type, days_in_year, days_in_month
use time_interp_external_mod, only: time_interp_external, init_external_field
use mpp_domains_mod,          only: domain2d
use constants_mod,            only: WTMCO2, WTMO2

use ocean_tpm_util_mod,       only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer
use fm_util_mod,              only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,              only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,              only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,              only: fm_util_start_namelist, fm_util_end_namelist
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,          only: ocean_thickness_type, ocean_density_type, ocean_time_type, ocean_grid_type
use ocmip2_co2calc_mod,       only: ocmip2_co2calc
use coupler_types_mod,        only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use atmos_ocean_fluxes_mod,   only: aof_set_coupler_flux
use ocean_util_mod,           only: diagnose_2d, diagnose_2d_comp, diagnose_3d_comp

implicit none

private

public  :: ocean_ibgc_bbc
public  :: ocean_ibgc_end
public  :: ocean_ibgc_init
public  :: ocean_ibgc_flux_init
public  :: ocean_ibgc_sbc
public  :: ocean_ibgc_source
public  :: ocean_ibgc_start
public  :: ocean_ibgc_init_sfc
public  :: ocean_ibgc_avg_sfc
public  :: ocean_ibgc_sum_sfc
public  :: ocean_ibgc_zero_sfc
public  :: ocean_ibgc_sfc_end
public  :: ocean_ibgc_tracer
public  :: ocean_ibgc_restart

private :: allocate_arrays

character(len=32), parameter              :: package_name = 'ocean_ibgc'
character(len=48), parameter              :: mod_name = 'ocean_ibgc_mod'
character(len=48), parameter              :: diag_name = 'ocean_ibgc'
character(len=fm_string_len), parameter   :: default_restart_file =  'ocean_ibgc.res.nc'
character(len=fm_string_len), parameter   :: default_local_restart_file =  'ocean_ibgc_local.res.nc'
character(len=fm_string_len), parameter   :: default_ice_restart_file   =    'ice_ibgc.res.nc'
character(len=fm_string_len), parameter   :: default_ocean_restart_file =  'ocean_ibgc_airsea_flux.res.nc'

! coefficients for O2 saturation
real, parameter :: a_0 = 2.00907
real, parameter :: a_1 = 3.22014
real, parameter :: a_2 = 4.05010
real, parameter :: a_3 = 4.94457
real, parameter :: a_4 = -2.56847e-01
real, parameter :: a_5 = 3.88767
real, parameter :: b_0 = -6.24523e-03
real, parameter :: b_1 = -7.37614e-03
real, parameter :: b_2 = -1.03410e-02
real, parameter :: b_3 = -8.17083e-03
real, parameter :: c_0 = -4.88682e-07

! Options

  logical :: do_ideal       = .true.  ! ideal_n and suntan
  logical :: do_po4         = .true.  ! incl PO4, DOP, PO4_pre
  logical :: do_gasses      = .true.  ! DIC, DI14C, O2 etc, requires po4 and/or po4f
  logical :: do_carbon_comp = .true.  ! use po4f for gases, PO4_pre, isotopes
  logical :: do_radiocarbon = .true.  ! use po4f for gases, PO4_pre, isotopes
  logical :: do_bgc_felim   = .true.  ! use po4f for gases, PO4_pre, isotopes
  logical :: do_po4f        = .true.  ! iFe and iron-limited PO4, DOP 
  logical :: do_isio4       = .true.  ! silica and d30Si
  logical :: do_no3_iso     = .true.  ! 15N and 18O of NO3

namelist /ocean_ibgc_nml/ do_ideal, do_po4, do_gasses, do_carbon_comp, do_radiocarbon, &
  do_bgc_felim, do_po4f, do_isio4, do_no3_iso

type ibgc_type

  real                                  :: kappa_eppley
  real                                  :: kappa_eppley_remin
  real                                  :: irrk

  real                                  :: ideal_n_vmax
  real                                  :: ideal_n_k
  real                                  :: ideal_n_r

  real                                  :: suntan_fade

  real                                  :: ipo4_vmax
  real                                  :: ipo4_k
  real                                  :: idop_frac
  real                                  :: gamma_idop
  real                                  :: gamma_ip 
  real                                  :: gamma_isi 
  real                                  :: sinking_init
  real                                  :: sinking_acc

  real                                  :: uptake_2_chl
  real                                  :: ichl_lag
  real                                  :: ichl_highlight

  real                                  :: ife_supply
  real                                  :: ligand_tot 
  real                                  :: ligand_k 
  real                                  :: scav_k 
  real                                  :: ipo4f_vmax
  real                                  :: ife_irrsuf 

!  real                                  :: alpha_12c_13c
  real                                  :: alpha_14n_15n
  real                                  :: alpha_28si_30si

  real                                  :: o2_2_ip
  real                                  :: io2_min

  real                                  :: half_life_14c
  real                                  :: lambda_14c
  real                                  :: alkbar
  real                                  :: si_2_ip
  real                                  :: c_2_ip
  real                                  :: sal_global
  real                                  :: frac_14catm_const = 0.0
  character(len=128)                    :: frac_14catm_file = ' '
  character(len=128)                    :: frac_14catm_name = ' '
  integer                               :: frac_14catm_id = 0

  real, _ALLOCATABLE, dimension(:,:,:)  :: jideal_n  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jsuntan  _NULL

  real, _ALLOCATABLE, dimension(:,:,:)  :: remin_ip  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: remin_isi  _NULL

  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pisi  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpisi  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pisi  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pi30si  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpi30si  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pi30si  _NULL

  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pip  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_idop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_idop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpip  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pip  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jife_new  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: felim_irrk  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pipf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_idopf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_idopf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpipf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pipf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pife  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpife  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pife  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: ife_free  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jife_scav  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jipo4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jipo4f  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jipo4_bgc  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: ipo4_bgc  _NULL

  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pi15n  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_ido15n  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_ido15n  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpi15n  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pi15n  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: ji15no3  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_i18o  _NULL

  real, _ALLOCATABLE, dimension(:,:,:)  :: ichl_new  _NULL

  real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: csatsurf  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: htotal  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: htotal_sat  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alk  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: c14surf  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: frac_14catm  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: pco2surf  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: pco2satsurf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdecay_idi14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdecay_ido14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pi14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_ido14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_ido14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jremin_pi14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpi14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jidi14c  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jio2  _NULL

  integer                               :: id_sc_co2 = -1
  integer                               :: id_sc_o2 = -1
  real, _ALLOCATABLE, dimension(:,:)    :: sc_co2  _NULL
  real                                  :: sc_co2_0
  real                                  :: sc_co2_1
  real                                  :: sc_co2_2
  real                                  :: sc_co2_3
  real, _ALLOCATABLE, dimension(:,:)    :: sc_o2  _NULL
  real                                  :: sc_o2_0
  real                                  :: sc_o2_1
  real                                  :: sc_o2_2
  real                                  :: sc_o2_3

  integer                               :: ind_irr = -1
  integer                               :: id_jideal_n = -1
  integer                               :: id_jsuntan = -1
  
  integer                               :: id_remin_ip = -1

  integer                               :: id_remin_isi = -1
  integer                               :: id_jprod_pisi = -1
  integer                               :: id_fpisi = -1
  integer                               :: id_jremin_pisi = -1
  integer                               :: id_jprod_pi30si = -1
  integer                               :: id_fpi30si = -1
  integer                               :: id_jremin_pi30si = -1

  integer                               :: id_jprod_pip = -1
  integer                               :: id_jprod_idop = -1
  integer                               :: id_jremin_idop = -1
  integer                               :: id_fpip = -1
  integer                               :: id_jremin_pip = -1
  integer                               :: id_jife_new = -1
  integer                               :: id_felim_irrk = -1
  integer                               :: id_jprod_pipf = -1
  integer                               :: id_jprod_idopf = -1
  integer                               :: id_jremin_idopf = -1
  integer                               :: id_fpipf = -1
  integer                               :: id_jremin_pipf = -1
  integer                               :: id_jprod_pife = -1
  integer                               :: id_fpife = -1
  integer                               :: id_jremin_pife = -1
  integer                               :: id_ife_free = -1
  integer                               :: id_jife_scav = -1
  integer                               :: id_jipo4 = -1
  integer                               :: id_jipo4f = -1
  integer                               :: id_jipo4_bgc = -1
  integer                               :: id_ipo4_bgc = -1

  integer                               :: id_jprod_pi15n = -1
  integer                               :: id_jprod_ido15n = -1
  integer                               :: id_jremin_ido15n = -1
  integer                               :: id_fpi15n = -1
  integer                               :: id_jremin_pi15n = -1
  integer                               :: id_ji15no3 = -1
  integer                               :: id_jprod_i18o = -1

  integer                               :: ind_ichl = -1
  integer                               :: ind_irrmem = -1
  integer                               :: id_ichl_new = -1

  integer                               :: id_alpha = -1
  integer                               :: id_csurf = -1
  integer                               :: id_csatsurf = -1
  integer                               :: id_htotal = -1
  integer                               :: id_htotal_sat = -1
  integer                               :: id_alk = -1
  integer                               :: id_c14surf = -1
  integer                               :: id_frac_14catm = -1
  integer                               :: id_pco2surf = -1
  integer                               :: id_pco2satsurf = -1
  integer                               :: id_sfc_flux_co2 = -1
  integer                               :: id_sfc_flux_co2sat = -1
  integer                               :: id_sfc_flux_14co2 = -1
  integer                               :: id_jdecay_idi14c = -1
  integer                               :: id_jdecay_ido14c = -1
  integer                               :: id_jprod_pi14c = -1
  integer                               :: id_jprod_ido14c = -1
  integer                               :: id_jremin_pi14c = -1
  integer                               :: id_jremin_ido14c = -1
  integer                               :: id_fpi14c = -1
  integer                               :: id_jidi14c = -1
  integer                               :: id_jio2 = -1
  integer                               :: id_sfc_flux_o2 = -1

  integer                               :: ind_ideal_n = -1
  integer                               :: ind_suntan = -1
  integer                               :: ind_ipo4_pre = -1
  integer                               :: ind_isio4 = -1
  integer                               :: ind_i30sio4 = -1
  integer                               :: ind_ipo4 = -1
  integer                               :: ind_idop = -1
  integer                               :: ind_i15no3 = -1
  integer                               :: ind_in18o3 = -1
  integer                               :: ind_ido15n = -1
  integer                               :: ind_ipo4f = -1
  integer                               :: ind_idopf = -1
  integer                               :: ind_ife = -1

  integer                               :: ind_idi14c = -1
  integer                               :: ind_ido14c = -1
  integer                               :: ind_idic = -1
  integer                               :: ind_idic_pre = -1
  integer                               :: ind_idic_sat = -1
  integer                               :: ind_co2_flux = -1
  integer                               :: ind_co2sat_flux = -1
  integer                               :: ind_14co2_flux = -1
  integer                               :: ind_io2 = -1
  integer                               :: ind_o2_flux = -1

  character(len=fm_string_len)          :: local_restart_file
  character(len=fm_field_name_len)      :: name

end type ibgc_type

logical, public :: do_ocean_ibgc
integer, public :: indsal
integer, public :: indtemp

integer                                 :: package_index
logical                                 :: module_initialized = .false.

character(len=128) :: version = '$Id: ocean_ibgc.F90,v 20.0 2013/12/14 00:09:32 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!       Input parameters:
!
!  htotal_in            = default value for htotal for an initial run
!  htotal_scale_lo      = scaling parameter to chose htotallo
!  htotal_scale_hi      = scaling parameter to chose htotalhi
real                                    :: htotal_in
real, allocatable, dimension(:,:)       :: htotal_scale_hi
real                                    :: htotal_scale_hi_in
real, allocatable, dimension(:,:)       :: htotal_scale_lo
real                                    :: htotal_scale_lo_in

! Calculated parameters (with possible initial input values):
integer                                         :: id_o2_sat = -1
real, allocatable, dimension(:,:)               :: sc_no_term
type(ibgc_type), allocatable, dimension(:) :: ibgc
integer                                         :: instances
real, allocatable, dimension(:,:,:)             :: o2_sat
real, allocatable, dimension(:)                 :: tk
real, allocatable, dimension(:)                 :: ts
real, allocatable, dimension(:)                 :: ts2
real, allocatable, dimension(:)                 :: ts3
real, allocatable, dimension(:)                 :: ts4
real, allocatable, dimension(:)                 :: ts5
real, allocatable, dimension(:)                 :: tt

! for restart
integer                              :: num_restart = 0
type(restart_file_type), allocatable :: restart(:)


contains

!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays for quantities with unknown dimensions.
! These are arrays that only exist temporarily.
! </DESCRIPTION>
!
subroutine allocate_arrays(isc, iec, jsc, jec, nk)

integer, intent(in)     :: isc
integer, intent(in)     :: iec
integer, intent(in)     :: jsc
integer, intent(in)     :: jec
integer, intent(in)     :: nk


integer :: i, j, k, n

allocate( sc_no_term(isc:iec,jsc:jec) )
allocate( htotal_scale_lo(isc:iec,jsc:jec) )
allocate( htotal_scale_hi(isc:iec,jsc:jec) )
allocate( o2_sat(isc:iec,jsc:jec,nk) )
allocate( tt(isc:iec) )
allocate( tk(isc:iec) )
allocate( ts(isc:iec) )
allocate( ts2(isc:iec) )
allocate( ts3(isc:iec) )
allocate( ts4(isc:iec) )
allocate( ts5(isc:iec) )

! initialize some arrays
sc_no_term(:,:)      = 0.0
htotal_scale_lo(:,:) = 0.0
htotal_scale_hi(:,:) = 0.0
o2_sat(:,:,:) = 0.0
tt(:) = 0.0
tk(:) = 0.0
ts(:) = 0.0
ts2(:) = 0.0
ts3(:) = 0.0
ts4(:) = 0.0
ts5(:) = 0.0

! allocate ibgc array elements
do n = 1, instances

  allocate( ibgc(n)%sc_co2(isc:iec,jsc:jec) )
  allocate( ibgc(n)%sc_o2(isc:iec,jsc:jec) )
  allocate( ibgc(n)%jideal_n(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jsuntan(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%remin_ip(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%remin_isi(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pisi(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpisi(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pisi(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pi30si(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpi30si(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pi30si(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pip(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_idop(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_idop(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpip(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pip(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jife_new(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%felim_irrk(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pipf(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_idopf(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_idopf(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpipf(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pipf(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pife(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpife(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pife(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%ife_free(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jife_scav(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jipo4(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jipo4f(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jipo4_bgc(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%ipo4_bgc(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pi15n(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_ido15n(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_ido15n(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpi15n(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pi15n(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%ji15no3(isc:iec,jsc:jec,nk) )  
  allocate( ibgc(n)%jprod_i18o(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%ichl_new(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%alpha(isc:iec,jsc:jec) )
  allocate( ibgc(n)%csurf(isc:iec,jsc:jec) )
  allocate( ibgc(n)%csatsurf(isc:iec,jsc:jec) )
  allocate( ibgc(n)%htotal(isc:iec,jsc:jec) )
  allocate( ibgc(n)%htotal_sat(isc:iec,jsc:jec) )
  allocate( ibgc(n)%alk(isc:iec,jsc:jec) )
  allocate( ibgc(n)%c14surf(isc:iec,jsc:jec) )
  allocate( ibgc(n)%frac_14catm(isc:iec,jsc:jec) )
  allocate( ibgc(n)%pco2surf(isc:iec,jsc:jec) )
  allocate( ibgc(n)%pco2satsurf(isc:iec,jsc:jec) )
  allocate( ibgc(n)%jdecay_idi14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jdecay_ido14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_pi14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jremin_ido14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_pi14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jprod_ido14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%fpi14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jidi14c(isc:iec,jsc:jec,nk) )
  allocate( ibgc(n)%jio2(isc:iec,jsc:jec,nk) )

enddo

do n = 1, instances

  do j = jsc, jec
    do i = isc, iec
     ibgc(n)%sc_co2(i,j) = 0.0
     ibgc(n)%sc_o2(i,j) = 0.0
     ibgc(n)%alpha(i,j) = 0.0
     ibgc(n)%csurf(i,j) = 0.0
     ibgc(n)%csatsurf(i,j) = 0.0
     ibgc(n)%htotal(i,j) = 0.0
     ibgc(n)%htotal_sat(i,j) = 0.0
     ibgc(n)%alk(i,j) = 0.0
     ibgc(n)%c14surf(i,j) = 0.0
     ibgc(n)%frac_14catm(i,j) = ibgc(n)%frac_14catm_const
     ibgc(n)%pco2surf(i,j) = 0.0
     ibgc(n)%pco2satsurf(i,j) = 0.0
      do k = 1, nk
       ibgc(n)%jideal_n(i,j,k) = 0.0
       ibgc(n)%jsuntan(i,j,k) = 0.0
       ibgc(n)%remin_ip(i,j,k) = 0.0
       ibgc(n)%remin_isi(i,j,k) = 0.0
       ibgc(n)%jprod_pisi(i,j,k) = 0.0
       ibgc(n)%fpisi(i,j,k) = 0.0
       ibgc(n)%jremin_pisi(i,j,k) = 0.0
       ibgc(n)%jprod_pi30si(i,j,k) = 0.0
       ibgc(n)%fpi30si(i,j,k) = 0.0
       ibgc(n)%jremin_pi30si(i,j,k) = 0.0
       ibgc(n)%jprod_pip(i,j,k) = 0.0
       ibgc(n)%jprod_idop(i,j,k) = 0.0
       ibgc(n)%jremin_idop(i,j,k) = 0.0
       ibgc(n)%fpip(i,j,k) = 0.0
       ibgc(n)%jremin_pip(i,j,k) = 0.0
       ibgc(n)%jife_new(i,j,k) = 0.0
       ibgc(n)%felim_irrk(i,j,k) = 0.0
       ibgc(n)%jprod_pipf(i,j,k) = 0.0
       ibgc(n)%jprod_idopf(i,j,k) = 0.0
       ibgc(n)%jremin_idopf(i,j,k) = 0.0
       ibgc(n)%fpipf(i,j,k) = 0.0
       ibgc(n)%jremin_pipf(i,j,k) = 0.0
       ibgc(n)%jprod_pife(i,j,k) = 0.0
       ibgc(n)%fpife(i,j,k) = 0.0
       ibgc(n)%jremin_pife(i,j,k) = 0.0
       ibgc(n)%ife_free(i,j,k) = 0.0
       ibgc(n)%jife_scav(i,j,k) = 0.0
       ibgc(n)%jipo4(i,j,k) = 0.0
       ibgc(n)%jipo4f(i,j,k) = 0.0
       ibgc(n)%jipo4_bgc(i,j,k) = 0.0
       ibgc(n)%ipo4_bgc(i,j,k) = 0.0
       ibgc(n)%jprod_pi15n(i,j,k) = 0.0
       ibgc(n)%jprod_ido15n(i,j,k) = 0.0
       ibgc(n)%jremin_ido15n(i,j,k) = 0.0
       ibgc(n)%fpi15n(i,j,k) = 0.0
       ibgc(n)%jremin_pi15n(i,j,k) = 0.0
       ibgc(n)%ji15no3(i,j,k) = 0.0
       ibgc(n)%jprod_i18o(i,j,k) = 0.0
       ibgc(n)%ichl_new(i,j,k) = 0.0
       ibgc(n)%jdecay_idi14c(i,j,k) = 0.0
       ibgc(n)%jdecay_ido14c(i,j,k) = 0.0
       ibgc(n)%jremin_pi14c(i,j,k) = 0.0
       ibgc(n)%jremin_ido14c(i,j,k) = 0.0
       ibgc(n)%jprod_pi14c(i,j,k) = 0.0
       ibgc(n)%jprod_ido14c(i,j,k) = 0.0
       ibgc(n)%fpi14c(i,j,k) = 0.0
       ibgc(n)%jidi14c(i,j,k) = 0.0
       ibgc(n)%jio2(i,j,k) = 0.0
      enddo
    enddo
  enddo

enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_bbc">
!
! <DESCRIPTION>
!     This sets up the boundary conditions at the bottom of the water 
! column. Here, this does nothing and is just a placeholder.
! </DESCRIPTION>
!

subroutine ocean_ibgc_bbc

end subroutine  ocean_ibgc_bbc
! </SUBROUTINE> NAME="ocean_ibgc_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_end">
!
! <DESCRIPTION>
!     Clean up various quantities for this run. This includes writing out 
! additional information to ensure reproduction across restarts.
! </DESCRIPTION>

subroutine ocean_ibgc_end()

character(len=64), parameter    :: sub_name = 'ocean_ibgc_end'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

integer :: i, j, k, n

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! save out additional information for a restart

write(stdoutunit,*)
call ocean_ibgc_restart

do n = 1, instances

  write(stdoutunit,*) trim(note_header),                         &
       'Writing additional restart information for instance ', &
       trim(ibgc(n)%name)

  write (stdoutunit,*) trim(note_header),                        &
       'Done writing additional restart information for instance ',&
       trim(ibgc(n)%name)

enddo

return
end subroutine  ocean_ibgc_end
! </SUBROUTINE> NAME="ocean_ibgc_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_ibgc_restart(time_stamp)
  character(len=*),             intent(in), optional :: time_stamp
  integer :: n

character(len=64), parameter    :: sub_name = 'ocean_ibgc_restart'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
  integer :: stdoutunit 
  stdoutunit=stdout() 


  do n=1, instances
     call save_restart(restart(n), time_stamp)
  end do

  ! perform checksums on restart fields
write (stdoutunit,*)
write (stdoutunit,*) trim(note_header), ' Saved check sums for extra variables'
do n = 1, instances
  write(stdoutunit,*) 'htotal chksum = ', mpp_chksum(ibgc(n)%htotal)
  write(stdoutunit,*) 'htotal_sat chksum = ', mpp_chksum(ibgc(n)%htotal_sat)
enddo
write (stdoutunit,*)

end subroutine ocean_ibgc_restart
! </SUBROUTINE> NAME="ocean_ibgc_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions. This includes things
! like gas exchange, atmospheric deposition, and riverine inputs.
! </DESCRIPTION>

subroutine ocean_ibgc_sbc(isc, iec, jsc, jec, &
     isc_bnd, jsc_bnd, T_prog,     &
     Grid, Time, ice_ocean_boundary_fluxes)

integer, intent(in)                                       :: isc
integer, intent(in)                                       :: iec
integer, intent(in)                                       :: jsc
integer, intent(in)                                       :: jec
integer, intent(in)                                       :: isc_bnd
integer, intent(in)                                       :: jsc_bnd
type(ocean_prog_tracer_type), intent(inout), dimension(:) :: T_prog
type(ocean_grid_type), intent(in)                         :: Grid
type(ocean_time_type), intent(in)                         :: Time
type(coupler_2d_bc_type), intent(in)                      :: ice_ocean_boundary_fluxes

integer :: i, j, k, n, m
integer :: i_bnd_off
integer :: j_bnd_off
integer :: kz
integer :: hour
integer :: day
real    :: days_in_this_year
integer :: minute
integer :: month
integer :: num_days
integer :: second
integer :: year
logical :: used

!     use the surface fluxes from the coupler
!       stf is in mol/m^2/s, flux from coupler is positive upwards
i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

if (do_gasses) then

do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
  t_prog(ibgc(n)%ind_idic)%stf(i,j) =                                     &
   -ice_ocean_boundary_fluxes%bc(ibgc(n)%ind_co2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
  t_prog(ibgc(n)%ind_io2)%stf(i,j) =                      &
   -ice_ocean_boundary_fluxes%bc(ibgc(n)%ind_o2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d(Time, Grid, ibgc(n)%id_sfc_flux_co2, t_prog(ibgc(n)%ind_idic)%stf(:,:))
   call diagnose_2d(Time, Grid, ibgc(n)%id_sfc_flux_o2, t_prog(ibgc(n)%ind_io2)%stf(:,:))
enddo


if (do_carbon_comp) then
do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
  t_prog(ibgc(n)%ind_idic_sat)%stf(i,j) =               &
   -ice_ocean_boundary_fluxes%bc(ibgc(n)%ind_co2sat_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo
do n = 1, instances
   call diagnose_2d(Time, Grid, ibgc(n)%id_sfc_flux_co2sat, t_prog(ibgc(n)%ind_idic_sat)%stf(:,:))
enddo
endif

if (do_radiocarbon) then
do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
  t_prog(ibgc(n)%ind_idi14c)%stf(i,j) =                  &
   -ice_ocean_boundary_fluxes%bc(ibgc(n)%ind_14co2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo
do n = 1, instances
   call diagnose_2d(Time, Grid, ibgc(n)%id_sfc_flux_14co2, t_prog(ibgc(n)%ind_idi14c)%stf(:,:))
enddo
endif

endif

return

end subroutine  ocean_ibgc_sbc
! </SUBROUTINE> NAME="ocean_ibgc_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocean_ibgc_flux_init

character(len=64), parameter    :: sub_name = 'ocean_ibgc_flux_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=256)                                      :: caller_str

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       First, perform some initialization if this module has not been
!       initialized because the normal initialization routine will
!       not have been called as part of the normal ocean model
!       initialization if this is an Atmosphere pe of a coupled
!       model running in concurrent mode
if (.not. module_initialized) then

   ! Initialize the package
  package_index = otpm_set_tracer_package(package_name,            &
       restart_file = default_restart_file,                        &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')

  ! Check whether to use this package
  path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
  instances = fm_get_length(path_to_names)
  if (instances .lt. 0) then
    call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
  endif

  write (stdoutunit,*)
  if (instances .eq. 0) then
    write (stdoutunit,*) trim(note_header), ' No instances'
    do_ocean_ibgc = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocean_ibgc = .true.
  endif

  module_initialized = .true.

endif

! return if we do not want to use this package
if (.not. do_ocean_ibgc) then
  return
endif

if (.not. allocated(ibgc)) then

   ! allocate storage for ibgc array
  allocate ( ibgc(instances) )

  ! loop over the names, saving them into the ibgc array
  do n = 1, instances

    if (fm_get_value(path_to_names, name, index = n)) then
      ibgc(n)%name = name
    else
      write (name,*) n
      call mpp_error(FATAL, trim(error_header) //        &
           'Bad field name for index ' // trim(name))
    endif

  enddo

endif

if (do_gasses) then
   ! Set up the ocean-atmosphere gas flux fields
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = ibgc(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // name
  endif

  ! Coupler fluxes
  ibgc(n)%ind_co2_flux = aof_set_coupler_flux('ico2_flux' // suffix,            &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',               &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),                     &
       ice_restart_file = default_ice_restart_file,                             &
       ocean_restart_file = default_ocean_restart_file,                         &
       caller = caller_str)

  ibgc(n)%ind_o2_flux = aof_set_coupler_flux('io2_flux' // suffix,              &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',               &
       mol_wt = WTMO2, param = (/ 9.36e-07, 9.7561e-06 /),                      &
       ice_restart_file = default_ice_restart_file,                             &
       ocean_restart_file = default_ocean_restart_file,                         &
       caller = caller_str)

if (do_carbon_comp) then
  ibgc(n)%ind_co2sat_flux = aof_set_coupler_flux('ico2sat_flux' // suffix,      &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',               &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),                     &
       ice_restart_file = default_ice_restart_file,                             &
       ocean_restart_file = default_ocean_restart_file,                         &
       caller = caller_str)
endif

if (do_radiocarbon) then
  ibgc(n)%ind_14co2_flux = aof_set_coupler_flux('i14co2_flux' // suffix,        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',               &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),                     &
       ice_restart_file = default_ice_restart_file,                             &
       ocean_restart_file = default_ocean_restart_file,                         &
       caller = caller_str)
endif

enddo

endif

return

end subroutine  ocean_ibgc_flux_init
! </SUBROUTINE> NAME="ocean_ibgc_flux_init"

!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocean_ibgc_init

character(len=64), parameter    :: sub_name = 'ocean_ibgc_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

real, parameter :: sperd = 24.0 * 3600.0

integer                                                 :: ioun
integer                                                 :: n
integer                                                 :: ierr
integer                                                 :: io_status
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_string_len)                            :: string
character(len=fm_field_name_len+3)                      :: long_suffix
logical, dimension(12)                                  :: t_mask
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  ! Initialize the package
package_index = otpm_set_tracer_package(package_name,                 &
     restart_file = default_restart_file,                             &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

! Check whether to use this package
path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then
  call mpp_error(FATAL, trim(error_header) // 'Could not get number of instances')
endif

write (stdoutunit,*)
if (instances .eq. 0) then
  write (stdoutunit,*) trim(note_header), ' No instances'
  do_ocean_ibgc = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocean_ibgc = .true.
endif

module_initialized = .true.

! Return if we do not want to use this package,
! after changing the list back
if (.not. do_ocean_ibgc) then
  return
endif

! provide for namelist over-ride
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_ibgc_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_ibgc_nml')
#else
ioun = open_namelist_file()
read  (ioun, ocean_ibgc_nml,iostat=io_status)
ierr = check_nml_error(io_status,'ocean_ibgc_nml')
call close_file (ioun)
#endif
write (stdoutunit,'(/)')
write (stdoutunit, ocean_ibgc_nml)  
write (stdlogunit, ocean_ibgc_nml)

do n = 1, instances

  if ((do_gasses) .and. ((do_po4) .or. (do_po4f))) then
    write (stdoutunit,*) trim(note_header), 'Doing iBGC gasses'
  else if ((do_gasses) .and. .not. ((do_po4) .or. (do_po4f))) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_gasses requires do_po4 and/or do_po4f' // trim(name))
  endif

  if ((do_carbon_comp) .and. ((do_po4) .or. (do_po4f)) .and. (do_gasses)) then
    write (stdoutunit,*) trim(note_header), 'Doing DIC components'
  else if ((do_gasses) .and. .not. ((do_po4) .or. (do_po4f)) .and. .not. (do_gasses)) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_carbon_comp requires do_po4 and/or do_po4f as well as do_carbon' // trim(name))
  endif
  
  if ((do_radiocarbon) .and. ((do_po4) .or. (do_po4f)) .and. (do_gasses)) then
    write (stdoutunit,*) trim(note_header), 'Doing radiocarbon'
  else if ((do_gasses) .and. .not. ((do_po4) .or. (do_po4f)) .and. .not. (do_gasses)) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_radiocarbon requires do_po4 and/or do_po4f as well as do_carbon' // trim(name))
  endif
  
  if ((do_no3_iso) .and. ((do_po4) .or. (do_po4f))) then
    write (stdoutunit,*) trim(note_header), 'Doing NO3 isotopes'
  else if ((do_gasses) .and. .not. ((do_po4) .or. (do_po4f))) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_no3_iso requires do_po4 and/or do_po4f' // trim(name))
  endif
  
  if ((do_bgc_felim) .and. (do_po4f)) then
    write (stdoutunit,*) trim(note_header), 'Using Fe-limited P for ibgc'
  else if ((do_bgc_felim) .and. .not. (do_po4f)) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_bgc_felim requires do_po4f' // trim(name))
  endif
  
enddo

!       Otherwise, go ahead to allocate storage for ibgc array.
allocate ( ibgc(instances) )

! Loop over the names, saving them into the ibgc array.
do n = 1, instances

  if (fm_get_value(path_to_names, name, index = n)) then
    ibgc(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = ibgc(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

if (do_ideal) then
   ! Ideal_N
  ibgc(n)%ind_ideal_n = otpm_set_prog_tracer('ideal_n' // suffix, package_name,   &
       longname = 'Ideal Nutrient' // trim(long_suffix),                          &
       units = 'mol kg-1', flux_units = 'mol m-2 s-1',                            &
       caller = caller_str)

  ! Suntan
  ibgc(n)%ind_suntan = otpm_set_prog_tracer('suntan' // suffix, package_name,     &
       longname = 'Suntan' // trim(long_suffix),                                  &
       units = 'J kg-1', flux_units = 'W/kg-1',                                   &
       caller = caller_str)
endif

if (do_po4) then
   ! PO4
  ibgc(n)%ind_ipo4 = otpm_set_prog_tracer('ipo4' // suffix, package_name,         &
       longname = 'Phosphate' // trim(long_suffix),                               &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! DOP
  ibgc(n)%ind_idop = otpm_set_prog_tracer('idop' // suffix, package_name,         &
       longname = 'DOP' // trim(long_suffix),                                     &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

endif

if (do_po4f) then
   ! PO4 Fe/Irr lim 
  ibgc(n)%ind_ipo4f = otpm_set_prog_tracer('ipo4f' // suffix, package_name,       &
       longname = 'Phosphate Fe/Irr lim' // trim(long_suffix),                    &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! DOPf
  ibgc(n)%ind_idopf = otpm_set_prog_tracer('idopf' // suffix, package_name,       &
       longname = 'DOPf' // trim(long_suffix),                                    &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! Fe Fe/Irr lim 
  ibgc(n)%ind_ife = otpm_set_prog_tracer('ife' // suffix, package_name,           &
       longname = 'Iron Fe/Irr lim' // trim(long_suffix),                         &
       units = 'mmol/kg', flux_units = 'mmol m-2 s-1',                            &
       caller = caller_str)

endif

if ((do_po4f) .or. (do_po4)) then
   ! Preformed PO4
  ibgc(n)%ind_ipo4_pre = otpm_set_prog_tracer('ipo4_pre' // suffix, package_name, &
       longname = 'Preformed iPO4' // trim(long_suffix),                          &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! Irradiance memory
  ibgc(n)%ind_irrmem = otpm_set_diag_tracer('irrmem' // suffix, package_name,   &
       longname = 'Irradiance memory' // trim(long_suffix),               &
       restart_file = default_restart_file,                         &
       units = 'Watts/m^2',const_init_tracer=.true.,    &
       const_init_value=0.0)

  ! Chlorophyll
  ibgc(n)%ind_ichl = otpm_set_diag_tracer('ichl' // suffix, package_name,         &
       longname = 'Chlorophyll' // trim(long_suffix),                             &
       restart_file = default_restart_file,                                       &
       units = 'ug/kg', const_init_tracer=.true.,                                 &
       const_init_value=0.08)

endif

if (do_gasses) then
   ! DIC
  ibgc(n)%ind_idic = otpm_set_prog_tracer('idic' // suffix, package_name,         &
       longname = 'DIC' // trim(long_suffix),                                     &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! O2
  ibgc(n)%ind_io2 = otpm_set_prog_tracer('io2' // suffix, package_name,           &
       longname = 'O2' // trim(long_suffix),                                      &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

if (do_carbon_comp) then
   ! Saturation DIC
  ibgc(n)%ind_idic_sat = otpm_set_prog_tracer('idic_sat' // suffix, package_name, &
       longname = 'Saturation DIC' // trim(long_suffix),                          &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! Preformed DIC
  ibgc(n)%ind_idic_pre = otpm_set_prog_tracer('idic_pre' // suffix, package_name, &
       longname = 'Preformed iDIC' // trim(long_suffix),                          &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)
endif

if (do_radiocarbon) then
  ! DI14C
  ibgc(n)%ind_idi14c = otpm_set_prog_tracer('idi14c' // suffix, package_name,     &
       longname = 'DI14C' // trim(long_suffix),                                   &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! DO14C
  ibgc(n)%ind_ido14c = otpm_set_prog_tracer('ido14c' // suffix, package_name,     &
       longname = 'DO14C' // trim(long_suffix),                                   &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)
endif

endif

if (do_no3_iso) then
   ! 15no3
  ibgc(n)%ind_i15no3 = otpm_set_prog_tracer('i15no3' // suffix, package_name,     &
       longname = '15NO3' // trim(long_suffix),                                   &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! DO15N
  ibgc(n)%ind_ido15n = otpm_set_prog_tracer('ido15n' // suffix, package_name,     &
       longname = 'DO15N' // trim(long_suffix),                                   &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! n18o3
  ibgc(n)%ind_in18o3 = otpm_set_prog_tracer('in18o3' // suffix, package_name,     &
       longname = 'N18O3' // trim(long_suffix),                                   &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

endif

if (do_isio4) then
   ! SiO4
  ibgc(n)%ind_isio4 = otpm_set_prog_tracer('isio4' // suffix, package_name,       &
       longname = 'Silicate' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

  ! 30SiO4
  ibgc(n)%ind_i30sio4 = otpm_set_prog_tracer('i30sio4' // suffix, package_name,   &
       longname = '30SiO4' // trim(long_suffix),                                  &
       units = 'mol/kg', flux_units = 'mol m-2 s-1',                              &
       caller = caller_str)

endif

enddo

! Process the namelists

! Add the package name to the list of good namelists, to be used
! later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

! Set up the *global* namelist
call fm_util_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call fm_util_set_value('htotal_scale_lo_in', 0.001 )
call fm_util_set_value('htotal_scale_hi_in', 1000.0)
call fm_util_set_value('htotal_in', 1.0e-08)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

! Set up the instance namelists
t_mask(:) = .true.

do n = 1, instances

   !       create the instance namelist
  call fm_util_start_namelist(package_name, ibgc(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('local_restart_file', default_local_restart_file)

  ! Global constants
  call fm_util_set_value('kappa_eppley', 0.0378)         ! deg C-1. 
  call fm_util_set_value('kappa_eppley_remin', 0.0249)   ! deg C-1. 
  call fm_util_set_value('irrk', 20.)                    !Light half-saturation constant W m-2

  ! Ideal Nutrient constants
  call fm_util_set_value('ideal_n_vmax', 0.007 / sperd)  ! Maximum uptake rate mol kg-1 s-1
  call fm_util_set_value('ideal_n_k', 0.2)               ! Nutrient half-saturation constant mol kg-1
  call fm_util_set_value('ideal_n_r', 0.00025 / sperd)   ! regeneration rate s-1

  ! Suntan constant
  call fm_util_set_value('suntan_fade', 0.01 )           ! Suntan fading rate W m-3 s-1

  ! Ideal Phosphate constants
  call fm_util_set_value('ipo4_vmax', 0.010 * 2.16e-6 / sperd) ! Maximum uptake rate mol kg-1 s-1 at 0C
  call fm_util_set_value('ipo4_k', 0.2 * 2.16e-6)        ! Nutrient half-saturation constant mol kg-1
  call fm_util_set_value('idop_frac', 0.65)              ! fraction of total production converted to dissolved pool
  call fm_util_set_value('gamma_idop', 0.007 / sperd)    ! remineralization rate for iDOP at 0C
  call fm_util_set_value('gamma_ip', 0.04 / sperd)       ! remineralization rate for soft tissue s-1
  call fm_util_set_value('sinking_init', 12. / sperd)    ! m s-1
  call fm_util_set_value('sinking_acc', .009 / sperd)    ! m s-2

  call fm_util_set_value('gamma_isi', 0.02 / sperd)      ! remineralization rate for silica

  ! Ideal Chlorophyll constants
  call fm_util_set_value('uptake_2_chl', 3.e13)          ! convert uptake rate to chl concentration mg-chl s molP-1
  call fm_util_set_value('ichl_lag', 0.5 / sperd)        ! Smoothing timescale for ichl adjustment s-1
  call fm_util_set_value('ichl_highlight', 0.1)          ! Fractional ichl:iP reduction under intense light

  ! Ideal Iron constants
  call fm_util_set_value('ipo4f_vmax', 0.012 * 2.16e-6 / sperd) ! Maximum uptake rate for PO4f mol kg-1 s-1 at 0C
  call fm_util_set_value('ligand_tot', 1. * 1.e-6)       !Ligand concentration mmol kg-1, Parekh 2005
  call fm_util_set_value('ligand_k', 1.e8)               ! Ligand stability constant, kg mmol-1 Dutkiewicz 2005
  call fm_util_set_value('scav_k', 0.001 / sperd)        ! Scavenging rate of free Fe, s-1, Dutkiewicz 2005 
  call fm_util_set_value('ife_irrsuf', 0.25 * 1.e-6)     !Iron-sufficient-ish Fe concentration for IRRk calc mmol kg-1
  call fm_util_set_value('ife_supply', 0.02  * 1.e-3 / sperd) ! Average global Fe input to surface mmol m-2 s-1

  ! Ideal isotope fractionations
!  call fm_util_set_value('alpha_12c_13c', 1.04)         ! 14C isotope fractionation factor (=2x13c)
  call fm_util_set_value('alpha_14n_15n', 1.005)         ! Nitrogen isotope fractionation factor
  call fm_util_set_value('alpha_28si_30si', 1.001)       ! Silicon isotope fractionation factor
  
  ! Gases 
  call fm_util_set_value('o2_2_ip', 170.)                ! Molar O2 to P ratio
  call fm_util_set_value('io2_min', 2.e-6)               ! Minimum O2 concentration

  call fm_util_set_value('sal_global', 35.0)             ! Average global salinity PSU
  call fm_util_set_value('half_life_14c', 5730.0)        ! a
  call fm_util_set_value('alkbar', 2.310e-03)            ! Average alkalinity eq/kg
  call fm_util_set_value('si_2_ip', 92. / 2.16)          ! Global average Si/PO4 concentration (to convert iSiO4 to SiO4)
  call fm_util_set_value('c_2_ip', 117.)                 ! Molar Redfield C to P ratio
  call fm_util_set_value('frac_14catm_file', ' ')      
  call fm_util_set_value('frac_14catm_name', ' ')
  call fm_util_set_value('frac_14catm_const', ibgc(n)%frac_14catm_const)

  !New Wanninkhof numbers
  call fm_util_set_value('sc_co2_0', 2068.9)
  call fm_util_set_value('sc_co2_1', -118.63)
  call fm_util_set_value('sc_co2_2', 2.9311)
  call fm_util_set_value('sc_co2_3', -0.027)
  call fm_util_set_value('sc_o2_0', 1929.7)
  call fm_util_set_value('sc_o2_1', -117.46)
  call fm_util_set_value('sc_o2_2', 3.116)
  call fm_util_set_value('sc_o2_3', -0.0306)

  call fm_util_end_namelist(package_name, ibgc(n)%name, check = .true., caller = caller_str)

enddo


! Check for any errors in the number of fields in the namelists for this package
good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then
  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif

return

end subroutine ocean_ibgc_init
! </SUBROUTINE> NAME="ocean_ibgc_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_init_sfc">
!
! <DESCRIPTION>
!
! SURFACE GAS FLUXES
! 
! This subroutine coordinates the calculation of gas solubilities in the surface layer, 
! and sends the appropriate values to the coupler.
!
! First, for CO2 and 14CO2, the carbon solubility and speciation are calculated by the
! subroutine co2calc, following the OCMIP2 protocol. These calculations are both made
! using total CO2, following which the surface CO2 concentration (CO2*, also known as
! H2CO3*) is scaled by the DI14C/DIC ratio to give the surface 14CO2 concentration.
! The speciation calculation uses in situ temperature, salinity, idealized PO4 and 
! idealized SiO4 (which must be scaled to real SiO4 units, since I'm using PO4 units
! for iSiO4). 
! Oxygen solubility is calculated here, using in situ temperature and salinity.
!
! The actual gas fluxes will be calculated in the coupler using a piston velocity (Kw),
!    Flux = Kw * (alpha - csurf)
! and returned as elements of ice_ocean_boundary_fluxes%bc in the ocean_ibgc_sbc 
! subroutine.

! </DESCRIPTION>

subroutine ocean_ibgc_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,   &
     isc_bnd, jsc_bnd,                                         &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer         :: i, j, k, n
integer         :: i_bnd_off
integer         :: j_bnd_off
integer         :: ind

real            :: epsln = 1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

if (do_gasses) then

do n = 1, instances

   ! CO2 flux
  ind = ibgc(n)%ind_co2_flux

  ! Have ocmip2_co2calc calculate alpha (which is the saturation co2star per atmosphere CO2, 
  ! in mol kg-1 atm-1), co2star (the actual co2star, in mol kg-1) and pco2surf (uatm)
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),  &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,            &
         grid_tmask(isd:ied,jsd:jed,1),                                    &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                   &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                    &
         t_prog(ibgc(n)%ind_idic)%field(isd:ied,jsd:jed,1,taum1),          &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1) *                   &
              ibgc(n)%alkbar / ibgc(n)%sal_global,                    &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1),     &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1) *    &
              ibgc(n)%si_2_ip,                                        &
         htotal_scale_lo, htotal_scale_hi, ibgc(n)%htotal,            &
         co2star = ibgc(n)%csurf, alpha = ibgc(n)%alpha,              &
         pco2surf = ibgc(n)%pco2surf)


    !  Compute the Schmidt number of CO2 in seawater using the 
    !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    !  7373-7382). This parameterizes the temperature-dependence of the
    !  gas exchange rate.
    do j = jsc, jec
      do i = isc, iec
        ibgc(n)%sc_co2(i,j) =                                              &
             ibgc(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *       &
             (ibgc(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *      &
              (ibgc(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *     &
               ibgc(n)%sc_co2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (ibgc(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)

        ! Send surface CO2* solubility (alpha) to coupler in units of mol m-3 atm-1
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) = &
             ibgc(n)%alpha(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
        ! Send surface CO2* concentration (csurf) to coupler in units of mol m-3 
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) = &
             ibgc(n)%csurf(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo

  endif

  ! O2 flux
  ind = ibgc(n)%ind_o2_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file), &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

!  Compute the oxygen saturation concentration at 1 atm total
!  pressure in mol/m^3 given the temperature (t, in deg C) and
!  the salinity (s, in permil)
!
!  From Garcia and Gordon (1992), Limnology and Oceonography.
!  The formula used is from page 1310, eq (8).
!
!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!  *** It shouldn't be there.                                ***
!
!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol/m^3
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        tt(i) = 298.15 - t_prog(indtemp)%field(i,j,k,taum1)
        tk(i) = 273.15 + t_prog(indtemp)%field(i,j,k,taum1)
        ts(i) = log(tt(i) / tk(i))
        ts2(i) = ts(i) * ts(i)
        ts3(i) = ts2(i) * ts(i)
        ts4(i) = ts3(i) * ts(i)
        ts5(i) = ts4(i) * ts(i)
        o2_sat(i,j,k) =                                                 &
             exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                         &
                 a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +                 &
                 t_prog(indsal)%field(i,j,1,taum1) *                    &
                 (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +           &
                  c_0*t_prog(indsal)%field(i,j,1,taum1)))
      enddo
    enddo
  enddo

  ! convert from ml/l to mol/m^3
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        o2_sat(i,j,k) = o2_sat(i,j,k) * (1000.0/22391.6)
      enddo
    enddo
  enddo
  
  !  Compute the Schmidt number of O2 in seawater using the 
  !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
  !  Cycles, 12, 141-163).
    do j = jsc, jec
      do i = isc, iec
        ibgc(n)%sc_o2(i,j) =                                               &
             ibgc(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *        &
             (ibgc(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *       &
              (ibgc(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *      &
               ibgc(n)%sc_o2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (ibgc(n)%sc_o2(i,j) + epsln)) * grid_tmask(i,j,1)
        
        ! Send the O2 solubility term (alpha) to the coupler, in units of mol m-3.
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) = &
             o2_sat(i,j,1) * sc_no_term(i,j)

        ! Send the in situ O2 concentration to the coupler, in units of mol m-3.
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             t_prog(ibgc(n)%ind_io2)%field(i,j,1,taum1) * rho(i,j,1,taum1) *       &
             sc_no_term(i,j)
      enddo
    enddo

  endif
  
if (do_carbon_comp) then
   !     DIC Sat CO2 flux
   !  This is similar to the CO2 flux above, but the Schmidt number term is increased  
   !  by 30x to cause rapid approach to saturation DIC concentrations.
  ind = ibgc(n)%ind_co2sat_flux

  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file), &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,                &
         grid_tmask(isd:ied,jsd:jed,1),                                        &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                       &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                        &
         t_prog(ibgc(n)%ind_idic_sat)%field(isd:ied,jsd:jed,1,taum1),          &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1) *                       &
              ibgc(n)%alkbar / ibgc(n)%sal_global,                             &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1),              &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1) *            &
              ibgc(n)%si_2_ip,                                                 &
         htotal_scale_lo, htotal_scale_hi, ibgc(n)%htotal_sat,                 &
         co2star = ibgc(n)%csatsurf,                                           &
         pco2surf = ibgc(n)%pco2satsurf)


    do j = jsc, jec
      do i = isc, iec

         ! Send surface CO2* solubility (alpha) to coupler in units of mol m-3 atm-1
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) = &
             ibgc(n)%alpha(i,j) * sc_no_term(i,j) * 30. * rho(i,j,1,taum1)
        ! Send surface CO2* concentration (csatsurf) to coupler in units of mol m-3.
        ! Schmidt number term is removed, replaced by constant multiplier of 1e4.
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) = &
             ibgc(n)%csatsurf(i,j) * sc_no_term(i,j) * 30. * rho(i,j,1,taum1)

      enddo
    enddo

  endif

endif

if (do_radiocarbon) then
    ! 14CO2 flux
    ind = ibgc(n)%ind_14co2_flux

    ! Calculate interpolated frac_14catm (fractionation of atmospheric 14CO2)
    if (ibgc(n)%frac_14catm_file .ne. ' ') then
      call time_interp_external(ibgc(n)%frac_14catm_id, model_time, ibgc(n)%frac_14catm)
    endif

    do j = jsc, jec
      do i = isc, iec

         ! Calculate surface 14CO2* concentration in mol kg-1, by scaling the total CO2* concentration 
         ! by 14C/12C
        ibgc(n)%c14surf(i,j) = ibgc(n)%csurf(i,j) *                              &
             t_prog(ibgc(n)%ind_idi14c)%field(i,j,1,taum1) /                     &
             (t_prog(ibgc(n)%ind_idic)%field(i,j,1,taum1) + epsln)

        ! Send surface 14CO2* solubility (alpha) to coupler in units of mol m-3 
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =  &
             ibgc(n)%alpha(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1) *           &
             (1.0 + ibgc(n)%frac_14catm(i,j) * 1.0e-03)

        ! Send surface 14CO2* concentration (csurf) to coupler in units of mol m-3 
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =  &
             ibgc(n)%c14surf(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo
endif

enddo

endif

return

end subroutine ocean_ibgc_init_sfc
! </SUBROUTINE> NAME="ocean_ibgc_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations. 
! </DESCRIPTION>
subroutine ocean_ibgc_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, jsc_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer         :: i, j, k, n
integer         :: i_bnd_off
integer         :: j_bnd_off
integer         :: ind

real            :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

if (do_gasses) then

do n = 1, instances

  ind = ibgc(n)%ind_co2_flux

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,           &
         grid_tmask(isd:ied,jsd:jed,1),                                    &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                   &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                    &
         t_prog(ibgc(n)%ind_idic)%field(isd:ied,jsd:jed,1,taum1),     &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1) *                   &
              ibgc(n)%alkbar / ibgc(n)%sal_global,               &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1),     &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1) *   &
              ibgc(n)%si_2_ip,                                        &
         htotal_scale_lo, htotal_scale_hi, ibgc(n)%htotal,            &
         co2star = ibgc(n)%csurf, alpha = ibgc(n)%alpha,         &
         pco2surf = ibgc(n)%pco2surf)

    ! Compute the Schmidt number of CO2 in seawater using the 
    ! formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    ! 7373-7382).
    do j = jsc, jec
      do i = isc, iec
        ibgc(n)%sc_co2(i,j) =                                              &
             ibgc(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *       &
             (ibgc(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *      &
              (ibgc(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *     &
               ibgc(n)%sc_co2_3)) * grid_tmask(i,j,1)
      sc_no_term(i,j) = sqrt(660.0 / (ibgc(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           ibgc(n)%alpha(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           ibgc(n)%csurf(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo

    ! O2 flux
    ind = ibgc(n)%ind_o2_flux

!  Compute the oxygen saturation concentration at 1 atm total
!  pressure in mol/m^3 given the temperature (t, in deg C) and
!  the salinity (s, in permil)
!
!  From Garcia and Gordon (1992), Limnology and Oceonography.
!  The formula used is from page 1310, eq (8).
!
!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!  *** It shouldn't be there.                                ***
!
!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol/m^3
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        tt(i) = 298.15 - t_prog(indtemp)%field(i,j,k,taum1)
        tk(i) = 273.15 + t_prog(indtemp)%field(i,j,k,taum1)
        ts(i) = log(tt(i) / tk(i))
        ts2(i) = ts(i) * ts(i)
        ts3(i) = ts2(i) * ts(i)
        ts4(i) = ts3(i) * ts(i)
        ts5(i) = ts4(i) * ts(i)
        o2_sat(i,j,k) =                                         &
             exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                 &
                 a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +         &
                 t_prog(indsal)%field(i,j,1,taum1) *            &
                 (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +   &
                  c_0*t_prog(indsal)%field(i,j,1,taum1)))
      enddo
    enddo
  enddo

  ! convert from ml/l to mol/m^3
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        o2_sat(i,j,k) = o2_sat(i,j,k) * (1000.0/22391.6)
      enddo
    enddo
  enddo

  !  Compute the Schmidt number of O2 in seawater using the 
  !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
  !  Cycles, 12, 141-163).
    do j = jsc, jec
      do i = isc, iec
        ibgc(n)%sc_o2(i,j) =                                               &
             ibgc(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *        &
             (ibgc(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *       &
              (ibgc(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *      &
               ibgc(n)%sc_o2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (ibgc(n)%sc_o2(i,j) + epsln)) * grid_tmask(i,j,1)

        ! Send the O2 solubility term (alpha) to the coupler, in units of mol m-3.
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +    &
             o2_sat(i,j,1) * sc_no_term(i,j)

        ! Send the in situ O2 concentration to the coupler, in units of mol m-3.
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +    &
             t_prog(ibgc(n)%ind_io2)%field(i,j,1,taum1) * rho(i,j,1,taum1) *       &
             sc_no_term(i,j)
      enddo
    enddo

if (do_carbon_comp) then

   !     DIC Sat CO2 flux
   !  This is similar to the CO2 flux above, but the Schmidt number term is increased  
   !  by 30x to cause rapid approach to saturation DIC concentrations.
  ind = ibgc(n)%ind_co2sat_flux

      call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,          &
         grid_tmask(isd:ied,jsd:jed,1),                                    &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                   &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                    &
         t_prog(ibgc(n)%ind_idic_sat)%field(isd:ied,jsd:jed,1,taum1),     &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1) *                   &
              ibgc(n)%alkbar / ibgc(n)%sal_global,               &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1),     &
         t_prog(ibgc(n)%ind_ipo4)%field(isd:ied,jsd:jed,1,taum1) *   &
              ibgc(n)%si_2_ip,                                        &
         htotal_scale_lo, htotal_scale_hi, ibgc(n)%htotal_sat,            &
         co2star = ibgc(n)%csatsurf,          &
         pco2surf = ibgc(n)%pco2satsurf)

    do j = jsc, jec
      do i = isc, iec
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           ibgc(n)%alpha(i,j) * sc_no_term(i,j) * 30. * rho(i,j,1,taum1)
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           ibgc(n)%csatsurf(i,j) * sc_no_term(i,j) * 30. * rho(i,j,1,taum1)
      enddo
    enddo

endif

if (do_radiocarbon) then
   ! 14CO2 flux
  ind = ibgc(n)%ind_14co2_flux

  ! Calculate interpolated frac_14catm (fractionation of atmospheric 14CO2)
  if (ibgc(n)%frac_14catm_file .ne. ' ') then
    call time_interp_external(ibgc(n)%frac_14catm_id, model_time, ibgc(n)%frac_14catm)
  endif

  do j = jsc, jec
    do i = isc, iec

       ! Calculate surface 14CO2* concentration in mol kg-1, by scaling the total CO2* concentration 
       ! by 14C/12C
      ibgc(n)%c14surf(i,j) = ibgc(n)%csurf(i,j) *                             &
           t_prog(ibgc(n)%ind_idi14c)%field(i,j,1,taum1) /                         &
           (t_prog(ibgc(n)%ind_idic)%field(i,j,1,taum1) + epsln)

      ! Send surface 14CO2* solubility (alpha) to coupler in units of mol m-3 
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           ibgc(n)%alpha(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1) *               &
           (1.0 + ibgc(n)%frac_14catm(i,j) * 1.0e-03)

      ! Send surface 14CO2* concentration (csurf) to coupler in units of mol m-3 
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           ibgc(n)%c14surf(i,j) * sc_no_term(i,j)  * rho(i,j,1,taum1)
    enddo
  enddo

endif

enddo

endif

return

end subroutine ocean_ibgc_sum_sfc
! </SUBROUTINE> NAME="ocean_ibgc_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations. 
! </DESCRIPTION>

subroutine ocean_ibgc_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

if (do_gasses) then

do n = 1, instances

  ind = ibgc(n)%ind_co2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = ibgc(n)%ind_o2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

if (do_carbon_comp) then
  ind = ibgc(n)%ind_co2sat_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0
endif

if (do_radiocarbon) then
  ind = ibgc(n)%ind_14co2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0
endif


enddo

endif

return

end subroutine ocean_ibgc_zero_sfc
! </SUBROUTINE> NAME="ocean_ibgc_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations. 
! </DESCRIPTION>

subroutine ocean_ibgc_avg_sfc(isc, iec, jsc, jec, isd, jsd,    &
     isc_bnd, jsc_bnd, Ocean_fields, Ocean_avg_kount, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
integer                                                 :: Ocean_avg_kount
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer         :: i_bnd_off
integer         :: j_bnd_off
integer         :: i, j, n
integer         :: ind
real            :: divid

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

divid = 1./float(Ocean_avg_kount)

if (do_gasses) then

do n = 1, instances

  ind = ibgc(n)%ind_co2_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo

  ind = ibgc(n)%ind_o2_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo

if (do_carbon_comp) then
  ind = ibgc(n)%ind_co2sat_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo
endif

if (do_radiocarbon) then
  ind = ibgc(n)%ind_14co2_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo
endif

enddo

endif

return

end subroutine ocean_ibgc_avg_sfc
! </SUBROUTINE> NAME="ocean_ibgc_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_sfc_end">
!
! <DESCRIPTION>
!       Finish up stuff for surface fields for flux calculations. 
! </DESCRIPTION>

subroutine ocean_ibgc_sfc_end

end subroutine ocean_ibgc_sfc_end
! </SUBROUTINE> NAME="ocean_ibgc_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_source">
!
! <DESCRIPTION>
!     Compute the source terms for the ideal nutrients, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>

subroutine ocean_ibgc_source(isc, iec, jsc, jec, nk, isd, jsd,       &
     T_prog, T_diag, taum1, grid_tmask, Grid, Time, grid_kmt, depth_zt,   &
     rho, rho_dzt, dzt, hblt_depth, dtts)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: jsd
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
type(ocean_diag_tracer_type), intent(inout), dimension(:)       :: T_diag
integer, intent(in)                                             :: taum1
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt
real, dimension(isd:,jsd:,:), intent(in)                        :: depth_zt
real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho
real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho_dzt
real, dimension(isd:,jsd:,:), intent(in)                        :: dzt
real, dimension(isd:,jsd:), intent(in)                          :: hblt_depth
real, intent(in)                                                :: dtts
integer :: i, j, k, n

integer :: ind_ideal_n
integer :: ind_suntan
integer :: ind_ipo4_pre
integer :: ind_idic_pre
integer :: ind_isio4
integer :: ind_i30sio4
integer :: ind_ipo4
integer :: ind_idop
integer :: ind_i15no3
integer :: ind_in18o3
integer :: ind_ido15n
integer :: ind_ichl
integer :: ind_ipo4f
integer :: ind_idopf
integer :: ind_ife
integer :: ind_irr
integer :: ind_irrmem
integer :: ind_io2
integer :: ind_idic
integer :: ind_idi14c
integer :: ind_ido14c

integer :: kblt

real :: epsln = 1.0e-30
real :: tmp_irr
real :: tmp_hblt

logical :: used

real,dimension(isc:iec,jsc:jec,1:nk)             :: expkT
real,dimension(isc:iec,jsc:jec,1:nk)             :: expkT_remin
real,dimension(isc:iec,jsc:jec,1:nk)             :: irr_mix

real,dimension(isc:iec,jsc:jec,1:nk)             :: wsink

real,dimension(isc:iec,jsc:jec,1:nk)             :: jprod_ip
real,dimension(isc:iec,jsc:jec,1:nk)             :: jprod_ipf
real,dimension(isc:iec,jsc:jec,1:nk)             :: jprod_ip_bgc
real,dimension(isc:iec,jsc:jec,1:nk)             :: jprod_i15n
real,dimension(isc:iec,jsc:jec,1:nk)             :: jprod_i14c

real,dimension(isc:iec,jsc:jec,1:nk)             :: fe_ligand
real,dimension(isc:iec,jsc:jec,1:nk)             :: ligand
real,dimension(isc:iec,jsc:jec,1:nk)             :: scav_k

real,dimension(isc:iec,jsc:jec,1:nk)             :: cideal_n
real,dimension(isc:iec,jsc:jec,1:nk)             :: cisio4
real,dimension(isc:iec,jsc:jec,1:nk)             :: ci30sio4
real,dimension(isc:iec,jsc:jec,1:nk)             :: cipo4
real,dimension(isc:iec,jsc:jec,1:nk)             :: cidop
real,dimension(isc:iec,jsc:jec,1:nk)             :: cido15n
real,dimension(isc:iec,jsc:jec,1:nk)             :: cichl
real,dimension(isc:iec,jsc:jec,1:nk)             :: cipo4f
real,dimension(isc:iec,jsc:jec,1:nk)             :: cidopf
real,dimension(isc:iec,jsc:jec,1:nk)             :: cife
real,dimension(isc:iec,jsc:jec,1:nk)             :: cidic
real,dimension(isc:iec,jsc:jec,1:nk)             :: cido14c

! SOURCES

! Calculate the source component fluxes for each of the ideal 
! nutrient tracers.

!  Loop over multiple instances
do n = 1, instances

! Use shortened naming convention for indices
  ind_irr = fm_get_index('/ocean_mod/diag_tracers/irr')
if (do_ideal) then
  ind_ideal_n = ibgc(n)%ind_ideal_n
  ind_suntan = ibgc(n)%ind_suntan
endif
if (do_po4) then
  ind_ipo4 = ibgc(n)%ind_ipo4
  ind_idop = ibgc(n)%ind_idop
endif
if (do_po4f) then
  ind_ipo4f = ibgc(n)%ind_ipo4f
  ind_idopf = ibgc(n)%ind_idopf
  ind_ife = ibgc(n)%ind_ife
endif
if ((do_po4f) .or. (do_po4)) then
  ind_irrmem = ibgc(n)%ind_irrmem
  ind_ichl = ibgc(n)%ind_ichl
endif
if (do_gasses) then
  ind_io2 = ibgc(n)%ind_io2
  ind_idic = ibgc(n)%ind_idic
endif
if (do_radiocarbon) then
  ind_idi14c = ibgc(n)%ind_idi14c
  ind_ido14c = ibgc(n)%ind_ido14c
endif
if (do_no3_iso) then
  ind_i15no3 = ibgc(n)%ind_i15no3
  ind_in18o3 = ibgc(n)%ind_in18o3
  ind_ido15n = ibgc(n)%ind_ido15n
endif
if (do_isio4) then
  ind_isio4 = ibgc(n)%ind_isio4
  ind_i30sio4 = ibgc(n)%ind_i30sio4
endif

! Calculate positive tracer concentrations for relevant variables
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        irr_mix(i,j,k) = t_diag(ind_irr)%field(i,j,k)
  enddo; enddo ; enddo
if (do_ideal) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cideal_n(i,j,k)  = max(0.0,t_prog(ind_ideal_n)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif
if (do_po4) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cipo4(i,j,k)  = max(0.0,t_prog(ind_ipo4)%field(i,j,k,taum1))
        cidop(i,j,k)  = max(0.0,t_prog(ind_idop)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif
if (do_po4f) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cipo4f(i,j,k)  = max(0.0,t_prog(ind_ipo4f)%field(i,j,k,taum1))
        cidopf(i,j,k)  = max(0.0,t_prog(ind_idopf)%field(i,j,k,taum1))
        cife(i,j,k)  = max(0.0,t_prog(ind_ife)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif
if ((do_po4f) .or. (do_po4)) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cichl(i,j,k)  = max(0.0,t_diag(ind_ichl)%field(i,j,k))
  enddo; enddo ; enddo
endif
if (do_gasses) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cidic(i,j,k)  = max(0.0,t_prog(ind_idic)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif
if (do_radiocarbon) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cido14c(i,j,k)  = max(0.0,t_prog(ind_ido14c)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif
if (do_no3_iso) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cido15n(i,j,k)  = max(0.0,t_prog(ind_ido15n)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif
if (do_isio4) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
        cisio4(i,j,k)  = max(0.0,t_prog(ind_isio4)%field(i,j,k,taum1))
        ci30sio4(i,j,k)  = max(0.0,t_prog(ind_i30sio4)%field(i,j,k,taum1))
  enddo; enddo ; enddo
endif


! The light used for uptake-limitation is the average irradiance within 
! the actively mixed layer, as defined in the KPP routine, plus an
! additional layer to account for mixing directly below the boundary.
  do j = jsc, jec
    do i = isc, iec
      kblt=1
      tmp_irr = irr_mix(i,j,1) * dzt(i,j,1)
      tmp_hblt = dzt(i,j,1)
      do k = 2, nk
        if (tmp_hblt .lt. hblt_depth(i,j)) then
          kblt=kblt+1
          tmp_irr = tmp_irr + irr_mix(i,j,k) * dzt(i,j,k)
          tmp_hblt = tmp_hblt + dzt(i,j,k)
        endif
      enddo
      irr_mix(i,j,1:kblt) = tmp_irr / max(1.0e-6,tmp_hblt)
    enddo
  enddo

  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
     ! Temperature functionality of export and remineralization

     ! After Eppley (1972). The eppley constant in use is 0.6 of the normal
     ! value (0.0378 rather than 0.063). In addition, the remineralization 
     ! temperature dependence is only 2/3 of that used for growth. This was 
     ! an arbitrary change, to reduce hypercycling in tropics.

       expkT(i,j,k) = exp(ibgc(n)%kappa_eppley *                           &
         t_prog(indtemp)%field(i,j,k,taum1))
       
       expkT_remin(i,j,k) = exp(ibgc(n)%kappa_eppley_remin *               &
          t_prog(indtemp)%field(i,j,k,taum1))

  enddo; enddo ; enddo

if (do_ideal) then
   ! NON-CONSERVATIVE TRACERS

   ! Ideal Nutrient

   ! Source function: 
   ! dN/dT = -Vmax * expkT * (1 - exp(-IRR/IRRk)) * (N / (N + k)) + expkT * R
   ! Do not allow the concentration of Ideal Nutrient to go above 1.
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

      if (cideal_n(i,j,k) .lt. 1.) then
        ibgc(n)%jideal_n(i,j,k) = -ibgc(n)%ideal_n_vmax *                  &
        expkT(i,j,k) * (1.0 - exp(- irr_mix(i,j,k) /                       &
        ibgc(n)%irrk)) * cideal_n(i,j,k) / (cideal_n(i,j,k) +              &
        ibgc(n)%ideal_n_k) + expkT(i,j,k) * ibgc(n)%ideal_n_r *            &
        grid_tmask(i,j,k)                
      else
        ibgc(n)%jideal_n(i,j,k) = -ibgc(n)%ideal_n_vmax * (1.0 -           &
        exp(- irr_mix(i,j,k) / ibgc(n)%irrk)) *                            &
        cideal_n(i,j,k) / (cideal_n(i,j,k) + ibgc(n)%ideal_n_k) *          &
        grid_tmask(i,j,k)
      endif

      ! Suntan

      ! Source function: dS/dT = IRR / (layer thickness) - f
      ! This needs to be fixed - not to be used now.
       if (t_prog(ind_suntan)%field(i,j,k,taum1) .gt. 0.0) then
         ibgc(n)%jsuntan(i,j,k) = (t_diag(ind_irr)%field(i,j,k)            &
          - ibgc(n)%suntan_fade) / dzt(i,j,k)  * grid_tmask(i,j,k)
       else
         ibgc(n)%jsuntan(i,j,k) = t_diag(ind_irr)%field(i,j,k)             &
         / dzt(i,j,k) * grid_tmask(i,j,k)
       endif

  enddo; enddo ; enddo

endif

! CONSERVATIVE TRACERS

! UPTAKE

! Ideal Phosphate
! Option provides 2 types of P: normal, and iron-limited.

if (do_po4) then
   ! Ideal Phosphate with DOP

   ! Ideal phosphate is done as before, but with a dissolved organic 
   ! phosphorus pool, after Yamanaka and Tajika (1997). A fraction (idop_frac) 
   ! of the total production (jprod_ip) is instantly converted to the dissolved 
   ! pool, while the remainder goes into the sinking particulate pool. This
   ! 'dissolved' organic phosphorus should be thought of as the sum of all 
   ! non-sinking organic matter, from small molecules to plankton to whales, 
   ! and is therefore not directly comparable to real DOP. Thus, the DOP
   ! concentrations and distributions should not necessarily match those of 
   ! the ocean - iBGC should probably have higher-than-observed concentrations
   ! in general, especially in nutrient-rich regions.

 do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

        jprod_ip(i,j,k) = ibgc(n)%ipo4_vmax * expkT(i,j,k) *               &
          (1.0 - exp(- irr_mix(i,j,k) / ibgc(n)%irrk)) *                   &         
          cipo4(i,j,k) / (cipo4(i,j,k) + ibgc(n)%ipo4_k) * grid_tmask(i,j,k)                              

        ibgc(n)%jprod_idop(i,j,k) = jprod_ip(i,j,k) * ibgc(n)%idop_frac                             

        ibgc(n)%jprod_pip(i,j,k) = jprod_ip(i,j,k) -                       &
          ibgc(n)%jprod_idop(i,j,k)

  enddo; enddo ; enddo

endif

if (do_po4f) then

   ! Iron cycling and growth limitation

   ! Iron is input in the upper 100m of the water column at a globally uniform 
   ! rate, intended to represent the sum of aeolian deposition and sediment-
   ! derived iron, without imposing geographical source distributions.

   ! Growth is limited by iron and phosphate concentrations, with iron 
   ! limitation expressed as an accentuator of light limitation, in 
   ! accordance with iron's role as an accessory in photosynthesis (e.g. 
   ! Strzepek et al., 2005). Iron uptake is a product of p uptake and  
   ! the iron concentration, to allow flexible stoichiometry and luxury 
   ! uptake. 

   ! The associated PO4 and DOP are called ipo4f and ipo4d, respectively.

   ! Iron taken up by organisms is assumed to be bound within organic matter
   ! until remineralized from a sinking particle. Sinking and remineralization 
   ! are identical to ipo4. Iron is also scavenged, further below.

   ! Iron source: spread over upper 100m and convert to mmol kg-1 s-1 
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

        if (depth_zt(i,j,k) .lt. 100.) then
          ibgc(n)%jife_new(i,j,k) = ibgc(n)%ife_supply * 0.01 /            &
            (rho(i,j,k,taum1) + epsln) * grid_tmask(i,j,k)
        else
          ibgc(n)%jife_new(i,j,k) = 0.0
        endif
        
        ! Growth and uptake
        ! dP/dT = -Vmax * expkT * (1 - exp(-IRR/IRRFek) * (P / P + kP)
        ! where IRRFek = IRRk * Fesuf / (Fe + x) where x is a small number to 
        ! prevent 'infinite' iron limitation.
        ! The idea here is that Fe helps phytoplankton harvest light, so that iron 
        ! limitation hampers the phytoplankton's ability to grow at a given light 
        ! level. With brighter light less iron is required, while with abundant iron 
        ! less light is required.
        ibgc(n)%felim_irrk(i,j,k) = ibgc(n)%irrk *                         &
           ibgc(n)%ife_irrsuf / (cife(i,j,k) + 1.e-8) 

        jprod_ipf(i,j,k) = ibgc(n)%ipo4f_vmax * expkT(i,j,k) *             &
          (1.0 - exp(- irr_mix(i,j,k) / ibgc(n)%felim_irrk(i,j,k)))*       &         
          min(cipo4f(i,j,k) / (cipo4f(i,j,k) + ibgc(n)%ipo4_k),            &
          cife(i,j,k) / ibgc(n)%ife_irrsuf)*                               &
          grid_tmask(i,j,k)                              

        ibgc(n)%jprod_idopf(i,j,k) = jprod_ipf(i,j,k) * ibgc(n)%idop_frac                             

        ibgc(n)%jprod_pipf(i,j,k) = jprod_ipf(i,j,k) -                     &
          ibgc(n)%jprod_idopf(i,j,k)

        ! Biological uptake        
        ibgc(n)%jprod_pife(i,j,k) =   ibgc(n)%jprod_pipf(i,j,k) *          &
          cife(i,j,k) / ibgc(n)%ife_irrsuf

  enddo; enddo ; enddo

endif

if ((do_po4f) .or. (do_po4)) then
! Does BGC use Fe-limited P, or non-Fe-limited P?
! Default is non-Fe-limited P (iPO4).

! This part sets the jprod_ip and po4 to be used for the fluxes of 
! other elements (chl, isotopes, O2, DIC etc), called *_bgc.
if (do_bgc_felim) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      jprod_ip_bgc(i,j,k) = jprod_ipf(i,j,k)
      ibgc(n)%ipo4_bgc(i,j,k) = t_prog(ind_ipo4f)%field(i,j,k,taum1)
  enddo; enddo ; enddo
else
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      jprod_ip_bgc(i,j,k) = jprod_ip(i,j,k)
      ibgc(n)%ipo4_bgc(i,j,k) = t_prog(ind_ipo4)%field(i,j,k,taum1)
  enddo; enddo ; enddo
endif

! Ideal Chlorophyll

! Chlorophyll is an important determinant of the absorption of shortwave 
! radiation in seawater, thereby influencing the circulation. In order to 
! capture this effect in a prognostic, interactive sense with the minimum 
! possible complexity, we link the formation of chlorophyll to the uptake 
! of inorganic nutrients, modulated by the available light (with 
! phytoplankton producing relatively less chlorophyll under very intense 
! light) following Behrenfeld et al., 2005. Because chlorophyll cannot 
! react instantaneously to changes in light and uptake, but would be 
! expected to respond only on the timescale of phytoplankton lifetimes 
! (days), we introduce two lag terms. First, an 'irradiance memory' term 
! that smoothes the available light, representing photoadaptation; second, 
! a chlorophyll lag that smoothes the diagnostic ichl concentration, 
! representing the delay of community growth and mortality, reducing
! short-term excursions (and preventing zero values appearing overnight!).
!
! A better way to do this would be to diagnose a biomass from some 
! nonlinear function of the growth rate, calculate the theta from the
! irrmem and iron (if used), and multiply them together to get a
! chlorophyll.
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

        t_diag(ind_irrmem)%field(i,j,k) = t_diag(ind_irrmem)%field(i,j,k)+ &
         (irr_mix(i,j,k) - t_diag(ind_irrmem)%field(i,j,k)) * min(1.0,     &
         ibgc(n)%ichl_lag * dtts) * grid_tmask(i,j,k)
         
        ibgc(n)%ichl_new(i,j,k) = jprod_ip_bgc(i,j,k) *                    &
          ibgc(n)%uptake_2_chl * (ibgc(n)%ichl_highlight + (1. -           &
          ibgc(n)%ichl_highlight) * exp (-0.05 *                           &
          t_diag(ind_irrmem)%field(i,j,k))) * grid_tmask(i,j,k)

        t_diag(ind_ichl)%field(i,j,k) = cichl(i,j,k) +                     & 
          (ibgc(n)%ichl_new(i,j,k) - cichl(i,j,k)) *                       &
          min(1.0, ibgc(n)%ichl_lag * dtts)         

  enddo; enddo ; enddo

endif

if (do_no3_iso) then
   ! Ideal Nitrate isotopes - with DOM

   ! Identical to implementation without DOM, except that a fraction of the 
   ! 15N taken up goes into the DOM pool. 
   ! The total uptake fractionation is the same, and there is no fractionation
   ! between particulate and dissolved organic forms.
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

       jprod_i15n(i,j,k) = jprod_ip_bgc(i,j,k) *                           &
         t_prog(ind_i15no3)%field(i,j,k,taum1) /                           &
         (epsln + ibgc(n)%ipo4_bgc(i,j,k)) /                               &
         ibgc(n)%alpha_14n_15n * grid_tmask(i,j,k)

       ibgc(n)%jprod_ido15n(i,j,k) = jprod_i15n(i,j,k) * ibgc(n)%idop_frac                             

       ibgc(n)%jprod_pi15n(i,j,k) = jprod_i15n(i,j,k) -                   &
         ibgc(n)%jprod_ido15n(i,j,k)

       ibgc(n)%jprod_i18o(i,j,k) = jprod_ip_bgc(i,j,k) *                  &
         t_prog(ind_in18o3)%field(i,j,k,taum1) /                          &
         (epsln + ibgc(n)%ipo4_bgc(i,j,k)) /                              &
         ibgc(n)%alpha_14n_15n * grid_tmask(i,j,k)

  enddo; enddo ; enddo
  
endif

if (do_radiocarbon) then
   ! Ideal Carbon isotopes - with DOM

   ! The uptake is done in the same way as N and  Si isotopes, to account for 
   ! the biological impact on 14C due to uptake and remineralization of DIC. 
   ! Because the convention for D14C notation includes an implicit correction 
   ! for mass-dependent fractionation (calculated from the observed d13C of 
   ! the same sample), uptake fractionation is not included.
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

       jprod_i14c(i,j,k) = jprod_ip_bgc(i,j,k) * ibgc(n)%c_2_ip *          &
         grid_tmask(i,j,k) * t_prog(ind_idi14c)%field(i,j,k,taum1) /       &
         (epsln + t_prog(ind_idic)%field(i,j,k,taum1))
         
       ibgc(n)%jprod_ido14c(i,j,k) = jprod_i14c(i,j,k) *                   &
         ibgc(n)%idop_frac                             

       ibgc(n)%jprod_pi14c(i,j,k) = jprod_i14c(i,j,k) -                    & 
         ibgc(n)%jprod_ido14c(i,j,k)

  enddo; enddo ; enddo

endif

if (do_isio4) then
   ! Ideal Silicate

   ! Uptake function: dN/dT = -Vmax * (1 - exp(-IRR/IRRk) * (N / N + k)
   ! Identical to Ideal Phosphate (without DOP), with different constants.
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

       ibgc(n)%jprod_pisi(i,j,k) = ibgc(n)%ipo4_vmax * expkT(i,j,k) *      &
          (1.0 - exp(- irr_mix(i,j,k) / ibgc(n)%irrk)) *                   &
          cisio4(i,j,k)  / (cisio4(i,j,k)                                  &
          + ibgc(n)%ipo4_k) * grid_tmask(i,j,k)                              

       ! Ideal Silicon isotopes
       ! The uptake is given by the uptake of SiO4 (total silicate), 
       ! multiplied by the rate ratio alpha.
       ibgc(n)%jprod_pi30si(i,j,k) = ibgc(n)%jprod_pisi(i,j,k) *           &
         t_prog(ind_i30sio4)%field(i,j,k,taum1) /                          &
         (epsln + t_prog(ind_isio4)%field(i,j,k,taum1)) /                  &
         ibgc(n)%alpha_28si_30si * grid_tmask(i,j,k)                     

  enddo; enddo ; enddo

endif

  ! SINKING AND REMINERALIZATION

  ! First, convert all ideal particulate production in the first layer to 
  ! particulate flux, no remineralization.
  do j=jsc,jec
    do i=isc,iec

      ibgc(n)%fpisi(i,j,1) = ibgc(n)%jprod_pisi(i,j,1) *                   &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pisi(i,j,1) = 0.0

      ibgc(n)%fpi30si(i,j,1) = ibgc(n)%jprod_pi30si(i,j,1) *               &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pi30si(i,j,1) = 0.0

      ibgc(n)%fpi14c(i,j,1) = ibgc(n)%jprod_pi14c(i,j,1) *                 &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pi14c(i,j,1) = 0.0

      ibgc(n)%fpip(i,j,1) = ibgc(n)%jprod_pip(i,j,1) *                     &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pip(i,j,1) = 0.0

      ibgc(n)%fpi15n(i,j,1) = ibgc(n)%jprod_pi15n(i,j,1) *               &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pi15n(i,j,1) = 0.0

      ibgc(n)%fpipf(i,j,1) = ibgc(n)%jprod_pipf(i,j,1) *                   &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pipf(i,j,1) = 0.0

      ibgc(n)%fpife(i,j,1) = ibgc(n)%jprod_pife(i,j,1) *                   &
        rho_dzt(i,j,1,taum1) * grid_tmask(i,j,1)
      ibgc(n)%jremin_pife(i,j,1) = 0.0

    enddo
  enddo

 do k = 2, nk
    do j = jsc, jec
      do i = isc, iec
         ! Rest of water column

         ! Calculate the remineralization lengthscale matrix, remin, a function of 
         ! T and z. Sinking rate (wsink) increases linearly with depth from an 
         ! initial value. Everything uses the remin_ip (no DOM) or remin_ip (with 
         ! DOM) values except for silica.
        wsink(i,j,k) = (ibgc(n)%sinking_acc * depth_zt(i,j,k) +            &
          ibgc(n)%sinking_init) * grid_tmask(i,j,k)

        ibgc(n)%remin_ip(i,j,k) = ibgc(n)%gamma_ip * expkT_remin(i,j,k) /  &
          (wsink(i,j,k) + epsln)
          
        ibgc(n)%remin_isi(i,j,k) = ibgc(n)%gamma_isi * expkT_remin(i,j,k)/ &
          (wsink(i,j,k) + epsln)

        ! Calculate the flux at level k from the flux in the overlying layer and
        ! the local remineralization lengthscale.
        ibgc(n)%fpisi(i,j,k) = ibgc(n)%fpisi(i,j,k-1) /                    &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_isi(i,j,k)) * grid_tmask(i,j,k)

        ibgc(n)%fpi30si(i,j,k) = ibgc(n)%fpi30si(i,j,k-1) /                &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_isi(i,j,k)) * grid_tmask(i,j,k)

        ibgc(n)%fpi14c(i,j,k) = ibgc(n)%fpi14c(i,j,k-1) /                  &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_ip(i,j,k)) * grid_tmask(i,j,k)

         ibgc(n)%fpip(i,j,k) = ibgc(n)%fpip(i,j,k-1) /                     &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_ip(i,j,k)) * grid_tmask(i,j,k)

        ibgc(n)%fpi15n(i,j,k) = ibgc(n)%fpi15n(i,j,k-1) /                  &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_ip(i,j,k)) * grid_tmask(i,j,k)

        ibgc(n)%fpipf(i,j,k) = ibgc(n)%fpipf(i,j,k-1) /                    &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_ip(i,j,k)) * grid_tmask(i,j,k)

        ibgc(n)%fpife(i,j,k) = ibgc(n)%fpife(i,j,k-1) /                    &
          (1.0 + dzt(i,j,k) * ibgc(n)%remin_ip(i,j,k)) * grid_tmask(i,j,k)

        ! Calculate regeneration term assuming flux through bottom of grid cell
        ibgc(n)%jremin_pisi(i,j,k) = (ibgc(n)%fpisi(i,j,k-1) -             &
          ibgc(n)%fpisi(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ibgc(n)%jremin_pi30si(i,j,k) = (ibgc(n)%fpi30si(i,j,k-1) -         &
          ibgc(n)%fpi30si(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ibgc(n)%jremin_pi14c(i,j,k) = (ibgc(n)%fpi14c(i,j,k-1) -           &
          ibgc(n)%fpi14c(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ibgc(n)%jremin_pip(i,j,k) = (ibgc(n)%fpip(i,j,k-1) -               &
          ibgc(n)%fpip(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ibgc(n)%jremin_pi15n(i,j,k) = (ibgc(n)%fpi15n(i,j,k-1) -           &
          ibgc(n)%fpi15n(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ibgc(n)%jremin_pipf(i,j,k) = (ibgc(n)%fpipf(i,j,k-1) -             &
          ibgc(n)%fpipf(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ibgc(n)%jremin_pife(i,j,k) = (ibgc(n)%fpife(i,j,k-1) -             &
          ibgc(n)%fpife(i,j,k)) * grid_tmask(i,j,k) / rho_dzt(i,j,k,taum1)

        ! Add production within box to flux assuming flux through bottom of
        ! grid cell
        ibgc(n)%fpisi(i,j,k) = (ibgc(n)%fpisi(i,j,k) +                     & 
          ibgc(n)%jprod_pisi(i,j,k) * rho_dzt(i,j,k,taum1)) *              & 
          grid_tmask(i,j,k)

        ibgc(n)%fpi30si(i,j,k) = (ibgc(n)%fpi30si(i,j,k) +                 & 
          ibgc(n)%jprod_pi30si(i,j,k) * rho_dzt(i,j,k,taum1)) *            &
          grid_tmask(i,j,k)

        ibgc(n)%fpi14c(i,j,k) = (ibgc(n)%fpi14c(i,j,k) +                   & 
          ibgc(n)%jprod_pi14c(i,j,k) * rho_dzt(i,j,k,taum1)) *             &  
          grid_tmask(i,j,k)

        ibgc(n)%fpip(i,j,k) = (ibgc(n)%fpip(i,j,k) +                       & 
          ibgc(n)%jprod_pip(i,j,k) * rho_dzt(i,j,k,taum1)) *               & 
          grid_tmask(i,j,k)

        ibgc(n)%fpi15n(i,j,k) = (ibgc(n)%fpi15n(i,j,k) +                   &  
          ibgc(n)%jprod_pi15n(i,j,k) * rho_dzt(i,j,k,taum1)) *             & 
          grid_tmask(i,j,k)

        ibgc(n)%fpipf(i,j,k) = (ibgc(n)%fpipf(i,j,k) +                     & 
          ibgc(n)%jprod_pipf(i,j,k) * rho_dzt(i,j,k,taum1)) *              &  
          grid_tmask(i,j,k)

        ibgc(n)%fpife(i,j,k) = (ibgc(n)%fpife(i,j,k) +                     & 
          ibgc(n)%jprod_pife(i,j,k) * rho_dzt(i,j,k,taum1)) *              & 
          grid_tmask(i,j,k)
 
      enddo
    enddo
  enddo

  ! Remineralize everything remaining in the bottom cell
  do j = jsc, jec
    do i = isc, iec
      k = grid_kmt(i,j)

      if (k .gt. 0) then

        ibgc(n)%jremin_pisi(i,j,k) = (ibgc(n)%jremin_pisi(i,j,k) +         &
          ibgc(n)%fpisi(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

        ibgc(n)%jremin_pi30si(i,j,k) = (ibgc(n)%jremin_pi30si(i,j,k) +     &
          ibgc(n)%fpi30si(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

        ibgc(n)%jremin_pi14c(i,j,k) = (ibgc(n)%jremin_pi14c(i,j,k) +       &
          ibgc(n)%fpi14c(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

        ibgc(n)%jremin_pip(i,j,k) = (ibgc(n)%jremin_pip(i,j,k) +           &
          ibgc(n)%fpip(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

        ibgc(n)%jremin_pi15n(i,j,k) = (ibgc(n)%jremin_pi15n(i,j,k) +       &
          ibgc(n)%fpi15n(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

        ibgc(n)%jremin_pipf(i,j,k) = (ibgc(n)%jremin_pipf(i,j,k) +         &
          ibgc(n)%fpipf(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

        ibgc(n)%jremin_pife(i,j,k) = (ibgc(n)%jremin_pife(i,j,k) +         &
          ibgc(n)%fpife(i,j,k) / rho_dzt(i,j,k,taum1)) * grid_tmask(i,j,k)

      endif

    enddo
  enddo

   ! Remineralize the dissolved organic matter, as a function of the first 
   ! order decay constant, and sum the remineralization and uptake sources. 
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

      ibgc(n)%jremin_idop(i,j,k) = (ibgc(n)%gamma_idop *                   &
        expkT_remin(i,j,k) * cidop(i,j,k)) * grid_tmask(i,j,k)
      
      ibgc(n)%jipo4(i,j,k) = (ibgc(n)%jremin_idop(i,j,k) -                 &
          ibgc(n)%jprod_pip(i,j,k) - ibgc(n)%jprod_idop(i,j,k)) *          &
          grid_tmask(i,j,k)

      ibgc(n)%jremin_idopf(i,j,k) = (ibgc(n)%gamma_idop *                  &
        expkT_remin(i,j,k) * cidopf(i,j,k)) * grid_tmask(i,j,k)
      
      ibgc(n)%jipo4f(i,j,k) = (ibgc(n)%jremin_idopf(i,j,k) -               &
          ibgc(n)%jprod_pipf(i,j,k) - ibgc(n)%jprod_idopf(i,j,k)) *        &
          grid_tmask(i,j,k)

      ibgc(n)%jremin_ido15n(i,j,k) = (ibgc(n)%gamma_idop *                 &
        expkT_remin(i,j,k) * cido15n(i,j,k)) * grid_tmask(i,j,k)
      
      ibgc(n)%ji15no3(i,j,k) = (ibgc(n)%jremin_ido15n(i,j,k) -             &
          ibgc(n)%jprod_pi15n(i,j,k) - ibgc(n)%jprod_ido15n(i,j,k)) *      &
          grid_tmask(i,j,k)
          
      ibgc(n)%jremin_ido14c(i,j,k) = (ibgc(n)%gamma_idop *                 &
        expkT_remin(i,j,k) * cido14c(i,j,k)) * grid_tmask(i,j,k)
      
      ibgc(n)%jidi14c(i,j,k) = (ibgc(n)%jremin_ido14c(i,j,k) -             &
          ibgc(n)%jprod_pi14c(i,j,k) - ibgc(n)%jprod_ido14c(i,j,k)) *      &
          grid_tmask(i,j,k)
          
  enddo; enddo ; enddo

if ((do_po4f) .or. (do_po4)) then
! Does BGC use Fe-limited P, or non-Fe-limited P?
! Default is non-Fe-limited P (iPO4).
! This is done now for jipo4_bgc.
if (do_bgc_felim) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      ibgc(n)%jipo4_bgc(i,j,k) = ibgc(n)%jipo4f(i,j,k)
  enddo; enddo ; enddo
else
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      ibgc(n)%jipo4_bgc(i,j,k) = ibgc(n)%jipo4(i,j,k)
  enddo; enddo ; enddo
endif

endif

if (do_po4f) then

   ! Iron scavenging
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

! Calculate the Ligand complexation throughout the full domain.
! This is used for scavenging only; it's assumed that phytoplankton are 
! able to access all iron, regardless of complexation.
! The calculation used here came from Parekh.

        ligand(i,j,k) = (-ibgc(n)%ligand_k * cife(i,j,k) +                 &
          ibgc(n)%ligand_k * ibgc(n)%ligand_tot - 1. +                     &
          ((ibgc(n)%ligand_k * cife(i,j,k) - ibgc(n)%ligand_k *            &
          ibgc(n)%ligand_tot + 1.) **2 + 4. * ibgc(n)%ligand_k *           &
          ibgc(n)%ligand_tot) ** 0.5) / (2. * ibgc(n)%ligand_k) *          &
          grid_tmask(i,j,k)

        fe_ligand(i,j,k) = (ibgc(n)%ligand_tot - ligand(i,j,k)) *          &
          grid_tmask(i,j,k)

        ibgc(n)%ife_free(i,j,k) = cife(i,j,k) - fe_ligand(i,j,k)

        ! Calculate the loss of Fe due to scavenging. This is done in the
        ! absolute simplest way, as a first order rate constant.
        ibgc(n)%jife_scav(i,j,k) = ibgc(n)%scav_k * ibgc(n)%ife_free(i,j,k)

  enddo; enddo ; enddo

endif

if (do_gasses) then
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
! O2

! Oxygen is produced by photosynthesis and consumed by respiration, the two 
! biological processes which simultaneously cause, respectively, the 
! consumption and remineralization of inorganic nutrients. Assuming a 
! constant stoichiometric ratio between oxygen and iPO4 allows oxygen to be 
! calculated simply as a function of the iPO4 source terms. 
! Oxygen consumption is not allowed to occur below a very low threshold 
! concentration (o2_min), in order to avoid negative tracer values. Note 
! that reminerization continues at the same rate irrespective of oxygen 
! concentrations. This is not realistic, and will lead to very large O2
! minimum zones.

        if (t_prog(ind_io2)%field(i,j,k,taum1) .gt. ibgc(n)%io2_min) then
          ibgc(n)%jio2(i,j,k) = -ibgc(n)%o2_2_ip * grid_tmask(i,j,k) *      &
          ibgc(n)%jipo4_bgc(i,j,k)
        else
          ibgc(n)%jio2(i,j,k) = ibgc(n)%o2_2_ip * grid_tmask(i,j,k) *      &
          jprod_ip_bgc(i,j,k)
        endif
  enddo; enddo ; enddo

if (do_radiocarbon) then
    ! Compute DI14C decay. This depends on the decay rate constant, lambda_14c.
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec

        ibgc(n)%jdecay_idi14c(i,j,k) =                                     &
          t_prog(ind_idi14c)%field(i,j,k,taum1) * ibgc(n)%lambda_14c *     &
          grid_tmask(i,j,k)

        ibgc(n)%jdecay_ido14c(i,j,k) =                                     &
          t_prog(ind_ido14c)%field(i,j,k,taum1) * ibgc(n)%lambda_14c *     &
          grid_tmask(i,j,k)
             
  enddo; enddo ; enddo
endif

endif

! ACCUMULATE SOURCES IN TENDENCY TERMS
!
! Calculate source/sink tendency terms for each prognostic tracer as the
! sum of fluxes generated above. These tendency terms are passed to the
! ocean model in units of mol m-3.

! All prognostic terms go inside this loop
if (do_ideal) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_ideal_n)%th_tendency(i,j,k) =                             &
        t_prog(ind_ideal_n)%th_tendency(i,j,k) + ibgc(n)%jideal_n(i,j,k)   &
        * rho_dzt(i,j,k,taum1)

      t_prog(ind_suntan)%th_tendency(i,j,k) =                              &
        t_prog(ind_suntan)%th_tendency(i,j,k) + ibgc(n)%jsuntan(i,j,k) *   &
        dzt(i,j,k)  
  enddo; enddo ; enddo
endif

if (do_po4) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_ipo4)%th_tendency(i,j,k) =                                &
        t_prog(ind_ipo4)%th_tendency(i,j,k) + ibgc(n)%jipo4(i,j,k) *       &
        rho_dzt(i,j,k,taum1)

      t_prog(ind_idop)%th_tendency(i,j,k) =                                &
        t_prog(ind_idop)%th_tendency(i,j,k) + (ibgc(n)%jprod_idop(i,j,k) - &
        ibgc(n)%jremin_idop(i,j,k) + ibgc(n)%jremin_pip(i,j,k)) *          &
         rho_dzt(i,j,k,taum1)
  enddo; enddo ; enddo
endif
 
if (do_po4f) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_ipo4f)%th_tendency(i,j,k) =                               & 
        t_prog(ind_ipo4f)%th_tendency(i,j,k) + ibgc(n)%jipo4f(i,j,k) *     &
        rho_dzt(i,j,k,taum1)

      t_prog(ind_ife)%th_tendency(i,j,k) =                                 &
        t_prog(ind_ife)%th_tendency(i,j,k) + (ibgc(n)%jremin_pife(i,j,k) + &
        ibgc(n)%jife_new(i,j,k) - ibgc(n)%jprod_pife(i,j,k) -              &
         ibgc(n)%jife_scav(i,j,k)) * rho_dzt(i,j,k,taum1)

      t_prog(ind_idopf)%th_tendency(i,j,k) =                               & 
        t_prog(ind_idopf)%th_tendency(i,j,k) +                             &
        (ibgc(n)%jprod_idopf(i,j,k) - ibgc(n)%jremin_idopf(i,j,k) +        &
        ibgc(n)%jremin_pipf(i,j,k)) * rho_dzt(i,j,k,taum1)
  enddo; enddo ; enddo
endif
 
if (do_gasses) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_io2)%th_tendency(i,j,k) =                                 & 
        t_prog(ind_io2)%th_tendency(i,j,k) + ibgc(n)%jio2(i,j,k) *         &
        rho_dzt(i,j,k,taum1)

      t_prog(ind_idic)%th_tendency(i,j,k) =                                & 
        t_prog(ind_idic)%th_tendency(i,j,k) + rho_dzt(i,j,k,taum1) *       &
        ibgc(n)%c_2_ip *  ibgc(n)%jipo4_bgc(i,j,k)
  enddo; enddo ; enddo
endif

if (do_radiocarbon) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_idi14c)%th_tendency(i,j,k) =                              & 
        t_prog(ind_idi14c)%th_tendency(i,j,k) + (ibgc(n)%jidi14c(i,j,k) -  &
        ibgc(n)%jdecay_idi14c(i,j,k)) * rho_dzt(i,j,k,taum1)

      t_prog(ind_ido14c)%th_tendency(i,j,k) =                              &
        t_prog(ind_ido14c)%th_tendency(i,j,k) +                            &
        (ibgc(n)%jprod_ido14c(i,j,k) - ibgc(n)%jremin_ido14c(i,j,k) -      &
        ibgc(n)%jdecay_ido14c(i,j,k)+ibgc(n)%jremin_pi14c(i,j,k)) *        &
        rho_dzt(i,j,k,taum1)
  enddo; enddo ; enddo
endif

if (do_no3_iso) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_i15no3)%th_tendency(i,j,k) =                             &
        t_prog(ind_i15no3)%th_tendency(i,j,k) +                           &
        ibgc(n)%ji15no3(i,j,k) * rho_dzt(i,j,k,taum1)

      t_prog(ind_ido15n)%th_tendency(i,j,k) =                              &
        t_prog(ind_ido15n)%th_tendency(i,j,k) +                            &
        (ibgc(n)%jprod_ido15n(i,j,k) - ibgc(n)%jremin_ido15n(i,j,k) +      &
        ibgc(n)%jremin_pi15n(i,j,k)) * rho_dzt(i,j,k,taum1)
          
      t_prog(ind_in18o3)%th_tendency(i,j,k) =                             &
        t_prog(ind_in18o3)%th_tendency(i,j,k) +                           &
        (ibgc(n)%jremin_idop(i,j,k) -                                      &
        ibgc(n)%jprod_i18o(i,j,k)) * rho_dzt(i,j,k,taum1)
  enddo; enddo ; enddo
endif

if (do_isio4) then
   do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
      t_prog(ind_isio4)%th_tendency(i,j,k) =                               &
        t_prog(ind_isio4)%th_tendency(i,j,k) +                             &
        (ibgc(n)%jremin_pisi(i,j,k) - ibgc(n)%jprod_pisi(i,j,k)) *         &
        rho_dzt(i,j,k,taum1)

      t_prog(ind_i30sio4)%th_tendency(i,j,k) =                             &
        t_prog(ind_i30sio4)%th_tendency(i,j,k) +                           &
        (ibgc(n)%jremin_pi30si(i,j,k) - ibgc(n)%jprod_pi30si(i,j,k)) *     &
        rho_dzt(i,j,k,taum1)
  enddo; enddo ; enddo
endif

enddo

! Save variables for diagnostics
call diagnose_3d_comp(Time, Grid, id_o2_sat, o2_sat(:,:,:))

do n = 1, instances
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_sc_co2, ibgc(n)%sc_co2(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_sc_o2, ibgc(n)%sc_o2(:,:))
   if (ibgc(n)%id_jideal_n .gt. 0) then
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jideal_n, ibgc(n)%jideal_n(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jsuntan .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jsuntan, ibgc(n)%jsuntan(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_remin_ip, ibgc(n)%remin_ip(:,:,:))
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_remin_isi, ibgc(n)%remin_isi(:,:,:))
   if (ibgc(n)%id_jprod_pisi .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pisi, ibgc(n)%jprod_pisi(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpisi, ibgc(n)%fpisi(:,:,:))
   if (ibgc(n)%id_jremin_pisi .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pisi, ibgc(n)%jremin_pisi(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_pi30si .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pi30si, ibgc(n)%jprod_pi30si(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpi30si, ibgc(n)%fpi30si(:,:,:))
   if (ibgc(n)%id_jremin_pi30si .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pi30si, ibgc(n)%jremin_pi30si(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_pip .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pip, ibgc(n)%jprod_pip(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_idop .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_idop, ibgc(n)%jprod_idop(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jremin_idop .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_idop, ibgc(n)%jremin_idop(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpip, ibgc(n)%fpip(:,:,:))
   if (ibgc(n)%id_jremin_pip .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pip, ibgc(n)%jremin_pip(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jife_new .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jife_new, ibgc(n)%jife_new(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_felim_irrk, ibgc(n)%felim_irrk(:,:,:))
   if (ibgc(n)%id_jprod_pipf .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pipf, ibgc(n)%jprod_pipf(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_idopf .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_idopf, ibgc(n)%jprod_idopf(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jremin_idopf .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_idopf, ibgc(n)%jremin_idopf(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpipf, ibgc(n)%fpipf(:,:,:))
   if (ibgc(n)%id_jremin_pipf .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pipf, ibgc(n)%jremin_pipf(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_pife .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pife, ibgc(n)%jprod_pife(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpife, ibgc(n)%fpife(:,:,:))
   if (ibgc(n)%id_jremin_pife .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pife, ibgc(n)%jremin_pife(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_ife_free, ibgc(n)%ife_free(:,:,:))
   if (ibgc(n)%id_jife_scav .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jife_scav, ibgc(n)%jife_scav(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jipo4 .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jipo4, ibgc(n)%jipo4(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jipo4f .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jipo4f, ibgc(n)%jipo4f(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jipo4_bgc .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jipo4_bgc, ibgc(n)%jipo4_bgc(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_ipo4_bgc, ibgc(n)%ipo4_bgc(:,:,:))
   if (ibgc(n)%id_jprod_pi15n .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pi15n, ibgc(n)%jprod_pi15n(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_ido15n .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_ido15n, ibgc(n)%jprod_ido15n(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jremin_ido15n .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_ido15n, ibgc(n)%jremin_ido15n(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpi15n, ibgc(n)%fpi15n(:,:,:))
   if (ibgc(n)%id_jremin_pi15n .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pi15n, ibgc(n)%jremin_pi15n(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_ji15no3 .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_ji15no3, ibgc(n)%ji15no3(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_i18o .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_i18o, ibgc(n)%jprod_i18o(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_ichl_new, ibgc(n)%ichl_new(:,:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_alpha, ibgc(n)%alpha(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_csurf, ibgc(n)%csurf(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_csatsurf, ibgc(n)%csatsurf(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_c14surf, ibgc(n)%c14surf(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_pco2surf, ibgc(n)%pco2surf(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_pco2satsurf, ibgc(n)%pco2satsurf(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_htotal, ibgc(n)%htotal(:,:))
   call diagnose_2d_comp(Time, Grid, ibgc(n)%id_htotal_sat, ibgc(n)%htotal_sat(:,:))
   if (ibgc(n)%id_alk .gt. 0) then
      call diagnose_2d(Time, Grid, ibgc(n)%id_alk, t_prog(indsal)%field(:,:,1,taum1) * ibgc(n)%alkbar / ibgc(n)%sal_global)
   endif
   if (ibgc(n)%id_jdecay_idi14c .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jdecay_idi14c, ibgc(n)%jdecay_idi14c(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jdecay_ido14c .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jdecay_ido14c, ibgc(n)%jdecay_ido14c(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jremin_pi14c .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_pi14c, ibgc(n)%jremin_pi14c(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jremin_ido14c .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jremin_ido14c, ibgc(n)%jremin_ido14c(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_pi14c .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_pi14c, ibgc(n)%jprod_pi14c(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   if (ibgc(n)%id_jprod_ido14c .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jprod_ido14c, ibgc(n)%jprod_ido14c(:,:,:) * rho_dzt(:,:,:,taum1))
   endif
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_fpi14c, ibgc(n)%fpi14c(:,:,:))
   call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jidi14c, ibgc(n)%jidi14c(:,:,:))
   if (ibgc(n)%id_jio2 .gt. 0) then                     
      call diagnose_3d_comp(Time, Grid, ibgc(n)%id_jio2, ibgc(n)%jio2(:,:,:) * rho_dzt(:,:,:,taum1))
   endif

enddo

return

end subroutine  ocean_ibgc_source
! </SUBROUTINE> NAME="ocean_ibgc_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants for a given run
! and allocate diagnostic arrays
! </DESCRIPTION>

subroutine ocean_ibgc_start(isc, iec, jsc, jec, nk, isd, jsd, &
     model_time, grid_tmask,            &
     grid_tracer_axes, &
     mpp_domain2d)

real, parameter :: spery = 365.25 * 24.0 * 3600.0

character(len=64), parameter    :: sub_name = 'ocean_ibgc_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
integer, dimension(3), intent(in)                       :: grid_tracer_axes
type(domain2d), intent(in)                              :: mpp_domain2d

  real                                  :: ideal_n_vmax
  real                                  :: ideal_n_k
  real                                  :: ideal_n_r

  real                                  :: suntan_fade

  real                                  :: ipo4_vmax
  real                                  :: ipo4_k
  real                                  :: idop_frac
  real                                  :: gamma_idop
  real                                  :: gamma_ip 
  real                                  :: gamma_isi 
  real                                  :: sinking_init
  real                                  :: sinking_acc
 
  real                                  :: uptake_2_chl
  real                                  :: ichl_lag
  real                                  :: ichl_highlight

  real                                  :: ligand_tot 
  real                                  :: ligand_k 
  real                                  :: scav_k 

  real                                  :: ipo4f_vmax
  real                                  :: ife_irrsuf 
  real                                  :: ife_supply

!  real                                  :: alpha_12c_13c
  real                                  :: alpha_14n_15n
  real                                  :: alpha_28si_30si

  real                                  :: o2_2_ip
  real                                  :: io2_min

  real                                  :: half_life_14c
  real                                  :: lambda_14c
  real                                  :: alkbar
  real                                  :: si_2_ip
  real                                  :: c_2_ip
  real                                  :: sal_global

integer                                         :: i, j, l, n
character(len=fm_field_name_len+1)              :: suffix
character(len=fm_field_name_len+3)              :: long_suffix
character(len=256)                              :: caller_str
character(len=fm_string_len), allocatable       :: local_restart_file(:)
integer                                         :: ind
logical                                         :: fld_exist
integer                                 :: id_restart

  integer :: stdoutunit 
  stdoutunit=stdout() 

! Determine indices for temperature and salinity
indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif

indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif

! dynamically allocate the global Ideal Nut arrays
call allocate_arrays(isc, iec, jsc, jec, nk)

! save the *global* namelist values
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str)

htotal_scale_lo_in =  fm_util_get_real   ('htotal_scale_lo_in', scalar = .true.)
htotal_scale_hi_in =  fm_util_get_real   ('htotal_scale_hi_in', scalar = .true.)
htotal_in          =  fm_util_get_real   ('htotal_in', scalar = .true.)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str)
      
! set default values for htotal_scale bounds
htotal_scale_lo(:,:) = htotal_scale_lo_in
htotal_scale_hi(:,:) = htotal_scale_hi_in

! read in the namelists for each instance
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  call fm_util_start_namelist(package_name, ibgc(n)%name, caller = caller_str)

  ibgc(n)%local_restart_file      = fm_util_get_string ('local_restart_file', scalar = .true.)

  ibgc(n)%kappa_eppley            = fm_util_get_real ('kappa_eppley', scalar = .true.)
  ibgc(n)%kappa_eppley_remin      = fm_util_get_real ('kappa_eppley_remin', scalar = .true.)
  ibgc(n)%irrk                    = fm_util_get_real ('irrk', scalar = .true.)

  ibgc(n)%ideal_n_vmax            = fm_util_get_real ('ideal_n_vmax', scalar = .true.)
  ibgc(n)%ideal_n_k               = fm_util_get_real ('ideal_n_k', scalar = .true.)
  ibgc(n)%ideal_n_r               = fm_util_get_real ('ideal_n_r', scalar = .true.)

  ibgc(n)%suntan_fade             = fm_util_get_real ('suntan_fade', scalar = .true.)
 
  ibgc(n)%ipo4_vmax               = fm_util_get_real ('ipo4_vmax', scalar = .true.)
  ibgc(n)%ipo4_k                  = fm_util_get_real ('ipo4_k', scalar = .true.)
  ibgc(n)%idop_frac               = fm_util_get_real ('idop_frac', scalar = .true.)
  ibgc(n)%gamma_idop              = fm_util_get_real ('gamma_idop', scalar = .true.)
  ibgc(n)%gamma_ip                = fm_util_get_real ('gamma_ip', scalar = .true.)
  ibgc(n)%gamma_isi               = fm_util_get_real ('gamma_isi', scalar = .true.)
  ibgc(n)%sinking_init            = fm_util_get_real ('sinking_init', scalar = .true.)
  ibgc(n)%sinking_acc             = fm_util_get_real ('sinking_acc', scalar = .true.)

  ibgc(n)%uptake_2_chl            = fm_util_get_real ('uptake_2_chl', scalar = .true.)
  ibgc(n)%ichl_lag                = fm_util_get_real ('ichl_lag', scalar = .true.)
  ibgc(n)%ichl_highlight          = fm_util_get_real ('ichl_highlight', scalar = .true.)

  ibgc(n)%ligand_tot              = fm_util_get_real ('ligand_tot', scalar = .true.)
  ibgc(n)%ligand_k                = fm_util_get_real ('ligand_k', scalar = .true.)
  ibgc(n)%scav_k                  = fm_util_get_real ('scav_k', scalar = .true.)

  ibgc(n)%ipo4f_vmax              = fm_util_get_real ('ipo4f_vmax', scalar = .true.)
  ibgc(n)%ife_irrsuf              = fm_util_get_real ('ife_irrsuf', scalar = .true.)
  ibgc(n)%ife_supply              = fm_util_get_real ('ife_supply', scalar = .true.)

!  ibgc(n)%alpha_12c_13c           = fm_util_get_real ('alpha_12c_13c', scalar = .true.)
  ibgc(n)%alpha_14n_15n           = fm_util_get_real ('alpha_14n_15n', scalar = .true.)
  ibgc(n)%alpha_28si_30si         = fm_util_get_real ('alpha_28si_30si', scalar = .true.)

  ibgc(n)%o2_2_ip                 = fm_util_get_real ('o2_2_ip', scalar = .true.)
  ibgc(n)%io2_min                 = fm_util_get_real ('io2_min', scalar = .true.)

  ibgc(n)%frac_14catm_file        = fm_util_get_string ('frac_14catm_file', scalar = .true.)
  ibgc(n)%frac_14catm_name        = fm_util_get_string ('frac_14catm_name', scalar = .true.)
  ibgc(n)%frac_14catm_const       = fm_util_get_real   ('frac_14catm_const', scalar = .true.)
  ibgc(n)%half_life_14c           = fm_util_get_real ('half_life_14c', scalar = .true.)
  ibgc(n)%alkbar                  = fm_util_get_real ('alkbar', scalar = .true.)
  ibgc(n)%si_2_ip                 = fm_util_get_real ('si_2_ip', scalar = .true.)
  ibgc(n)%c_2_ip                  = fm_util_get_real ('c_2_ip', scalar = .true.)
  ibgc(n)%sal_global              = fm_util_get_real ('sal_global', scalar = .true.)

  ibgc(n)%sc_co2_0                = fm_util_get_real ('sc_co2_0', scalar = .true.)
  ibgc(n)%sc_co2_1                = fm_util_get_real ('sc_co2_1', scalar = .true.)
  ibgc(n)%sc_co2_2                = fm_util_get_real ('sc_co2_2', scalar = .true.)
  ibgc(n)%sc_co2_3                = fm_util_get_real ('sc_co2_3', scalar = .true.)
  ibgc(n)%sc_o2_0                 = fm_util_get_real ('sc_o2_0', scalar = .true.)
  ibgc(n)%sc_o2_1                 = fm_util_get_real ('sc_o2_1', scalar = .true.)
  ibgc(n)%sc_o2_2                 = fm_util_get_real ('sc_o2_2', scalar = .true.)
  ibgc(n)%sc_o2_3                 = fm_util_get_real ('sc_o2_3', scalar = .true.)

  call fm_util_end_namelist(package_name, ibgc(n)%name, caller = caller_str)


enddo

if (do_gasses) then

do n = 1, instances

if (do_radiocarbon) then
   ! Open the frac_14catm (fractionation of atmospheric 14CO2) file
   ! If the file name is blank, then the 14C fractionation is assumed to
   ! be added to the atmospheric concentration
  if (ibgc(n)%frac_14catm_file .ne. ' ') then
    ibgc(n)%frac_14catm_id = init_external_field(ibgc(n)%frac_14catm_file, &
      ibgc(n)%frac_14catm_name, domain = mpp_domain2d,                     &
      use_comp_domain = .true.)
    if (ibgc(n)%frac_14catm_id .eq. 0) then
      call mpp_error(FATAL, trim(error_header) //                          &
        ' Could not open frac_14catm_file file: ' //                       &
        trim(ibgc(n)%frac_14catm_file) // ' for ' //                       &
        trim(ibgc(n)%frac_14catm_name))
    endif
  else
    call mpp_error(NOTE, trim(error_header) //                             &
      ' Using constant field for atmospheric 14C for instance '            &
      // trim(ibgc(n)%name))
    do j = jsc, jec
      do i = isc, iec
        ibgc(n)%frac_14catm(i,j) = ibgc(n)%frac_14catm_const *             &
        grid_tmask(i,j,1)
      enddo
    enddo
  endif
endif

enddo

!       Read in additional information for a restart.
!
!       We must process all of the instances before restoring any files
!       as all fields must be registered before the fields are
!       restored, and fields from different instances may be in the
!       same file.
!
!       Note that the restart file names here must be different from
!       those for the tracer values.

allocate(restart(instances))
allocate(local_restart_file(instances))

write(stdoutunit,*)

do n = 1, instances

   ! Set the suffix for this instance (if instance name is "_",
   ! then use a blank suffix).
  if (ibgc(n)%name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // ibgc(n)%name
  endif

  ! Check whether we are already using this restart file, if so,
  ! we do not want to duplicate it in the list of restart files
  ! since we only read each restart file once.
  ind = 0
  do l = 1, num_restart
    if (ibgc(n)%local_restart_file == local_restart_file(l)) then
      ind = l
      exit
    endif
  end do

  if (ind .eq. 0) then
    num_restart = num_restart + 1
    ind = num_restart
    local_restart_file(ind) = trim(ibgc(n)%local_restart_file)
  end if

  ! Check whether the field already exists in the restart file.
  ! If not, then set a default value.
  fld_exist = field_exist('INPUT/' // trim(ibgc(n)%local_restart_file), 'htotal' // trim(suffix) )

  if ( fld_exist ) then
    write (stdoutunit,*) trim(note_header),                       &
         ' Reading additional information for instance ',       &
         trim(ibgc(n)%name)
  else
    write (stdoutunit,*) trim(note_header),                       &
         ' Initializing instance ', trim(ibgc(n)%name)
    ibgc(n)%htotal(:,:) = htotal_in
    ibgc(n)%htotal_sat(:,:) = htotal_in
  endif

  ! Register the field for restart
  id_restart = register_restart_field(restart(ind), ibgc(n)%local_restart_file,         &
                    'htotal' // trim(suffix), ibgc(n)%htotal,                           &
                    domain=mpp_domain2d, mandatory=fld_exist)

if (do_carbon_comp) then
  id_restart = register_restart_field(restart(ind), ibgc(n)%local_restart_file,         &
                    'htotal_sat' // trim(suffix), ibgc(n)%htotal_sat,                   &
                    domain=mpp_domain2d, mandatory=fld_exist)
endif

enddo

! Restore the restart fields if the file exists
do l = 1, num_restart
  if (file_exist('INPUT/' // trim(local_restart_file(l)))) then
    call restore_state(restart(l))
  end if
end do

deallocate(local_restart_file)

! perform checksums on restart fields
write (stdoutunit,*)
write (stdoutunit,*) trim(note_header), ' Starting check sums for extra variables'
do n = 1, instances
  write(stdoutunit,*) 'htotal chksum = ', mpp_chksum(ibgc(n)%htotal)
if (do_carbon_comp) then
  write(stdoutunit,*) 'htotal_sat chksum = ', mpp_chksum(ibgc(n)%htotal_sat)
endif
enddo
write (stdoutunit,*)

! initialize some arrays which are held constant for this
! simulation
if (do_radiocarbon) then
do n = 1, instances
  if (ibgc(n)%half_life_14c .gt. 0.0) then
    ibgc(n)%lambda_14c = log(2.0) / (ibgc(n)%half_life_14c * spery)
  else
    call mpp_error(FATAL,trim(error_header) // ' Half-life <= 0')
  endif
enddo
endif

endif

!       Set up diagnostic fields

! Register the global fields

suffix = '_' // package_name
long_suffix = ' (' // trim(package_name) // ')'

! Register the instance fields
do n = 1, instances

  if (instances .eq. 1) then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // ibgc(n)%name
    long_suffix = ' (' // trim(ibgc(n)%name) // ')'
  endif

if (do_ideal) then
  ibgc(n)%id_jideal_n = register_diag_field(trim(diag_name),                                  &
       'jideal_n' // trim(suffix), grid_tracer_axes(1:3),                                     &
       model_time, 'Ideal Nutrient source, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jsuntan = register_diag_field(trim(diag_name),                                   &
       'jsuntan' // trim(suffix), grid_tracer_axes(1:3),                                      &
       model_time, 'Suntan source, layer integral' // trim(long_suffix), 'mol m-2 s-1',       &
       missing_value = -1.0e+10)
endif

if (do_po4) then
  ibgc(n)%id_remin_ip = register_diag_field(trim(diag_name),                                  &
       'remin_ip' // trim(suffix), grid_tracer_axes(1:3),                                     &
       model_time, 'Remineralization inverse lengthscale for fpip' // trim(long_suffix), 'm-1',       &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_pip = register_diag_field(trim(diag_name),                                 &
       'jprod_pip' // trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Particulate ideal P production, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_idop = register_diag_field(trim(diag_name),                                &
       'jprod_idop' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'DOP source, layer integral' // trim(long_suffix), 'mol m-2 s-1',          &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_idop = register_diag_field(trim(diag_name),                               &
       'jremin_idop' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'DOP remineralization, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpip = register_diag_field(trim(diag_name),                                      &
       'fpip' // trim(suffix), grid_tracer_axes(1:3),                                         &
       model_time, 'Particulate ideal P sinking flux' // trim(long_suffix), 'mol m-2 s-1',    &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_pip = register_diag_field(trim(diag_name),                                &
       'jremin_pip' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Ideal PO4 remin fr particles, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jipo4 = register_diag_field(trim(diag_name),                                     &
       'jipo4' // trim(suffix), grid_tracer_axes(1:3),                                        &
       model_time, 'Ideal PO4 total source, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

endif

if (do_po4f) then

 ibgc(n)%id_jife_new = register_diag_field(trim(diag_name),                                    &
       'jife_new' // trim(suffix), grid_tracer_axes(1:3),                                      &
       model_time, 'iFe surface source, layer integral' // trim(long_suffix), 'mol m-2 s-1',   &
       missing_value = -1.0e+10)

  ibgc(n)%id_felim_irrk = register_diag_field(trim(diag_name),                                 &
       'felim_irrk' // trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Fe-limited light half-saturation' // trim(long_suffix), 'W m-2',           &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_pipf = register_diag_field(trim(diag_name),                                 &
       'jprod_pipf' // trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Particulate iPf production, layer integral' // trim(long_suffix), 'mmol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_idopf = register_diag_field(trim(diag_name),                                &
       'jprod_idopf' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'DOPf source, layer integral' // trim(long_suffix), 'mol m-2 s-1',          &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_idopf = register_diag_field(trim(diag_name),                               &
       'jremin_idopf' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'DOPf remineralization, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpipf = register_diag_field(trim(diag_name),                                      &
       'fpipf' // trim(suffix), grid_tracer_axes(1:3),                                         &
       model_time, 'Particulate iPf sinking flux' // trim(long_suffix), 'mol m-2 s-1',         &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_pipf = register_diag_field(trim(diag_name),                                &
       'jremin_pipf' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'iPO4f remineralization, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_pife = register_diag_field(trim(diag_name),                                 &
       'jprod_pife' // trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Particulate iFe production, layer integral' // trim(long_suffix), 'mmol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpife = register_diag_field(trim(diag_name),                                      &
       'fpife' // trim(suffix), grid_tracer_axes(1:3),                                         &
       model_time, 'Particulate iFe sinking flux' // trim(long_suffix), 'mmol m-2 s-1',        &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_pife = register_diag_field(trim(diag_name),                                &
       'jremin_pife' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'iFe remineralization, layer integral' // trim(long_suffix), 'mmol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_ife_free = register_diag_field(trim(diag_name),                                   &
       'ife_free' // trim(suffix), grid_tracer_axes(1:3),                                      &
       model_time, 'Free iFe concentration' // trim(long_suffix), 'mmol kg-1',                 &
       missing_value = -1.0e+10)

  ibgc(n)%id_jife_scav = register_diag_field(trim(diag_name),                                  &
       'jife_scav' // trim(suffix), grid_tracer_axes(1:3),                                     &
       model_time, 'iFe scavenging, layer integral' // trim(long_suffix), 'mmol m-2 s-1',      &
       missing_value = -1.0e+10)

  ibgc(n)%id_jipo4f = register_diag_field(trim(diag_name),                                     &
       'jipo4f' // trim(suffix), grid_tracer_axes(1:3),                                        &
       model_time, 'Ideal PO4f net source, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

endif

if ((do_po4f) .or. (do_po4)) then

  ibgc(n)%id_ichl_new = register_diag_field(trim(diag_name),                                   &
       'ichl_new' // trim(suffix), grid_tracer_axes(1:3),                                      &
       model_time, 'New ichl' // trim(long_suffix), 'mg kg-1',                                 &
       missing_value = -1.0e+10)

  ibgc(n)%id_ipo4_bgc = register_diag_field(trim(diag_name),                                   &
       'ipo4_bgc' // trim(suffix), grid_tracer_axes(1:3),                                      &
       model_time, 'PO4 used for bgc calculations' // trim(long_suffix), 'mol kg-1',           &
       missing_value = -1.0e+10)

  ibgc(n)%id_jipo4_bgc = register_diag_field(trim(diag_name),                                  &
       'jipo4_bgc' // trim(suffix), grid_tracer_axes(1:3),                                     &
       model_time, 'Net PO4 source for bgc calculations, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

endif

if (do_gasses) then

  id_o2_sat = register_diag_field(trim(diag_name),                                             &
       'o2_sat' // trim(suffix), grid_tracer_axes(1:3),                                        &
       model_time, 'O2 saturation' // trim(long_suffix), 'mol m-3',                            &
       missing_value = -1.0e+10)

  ibgc(n)%id_sc_co2 = register_diag_field(trim(diag_name),                                     &
       'sc_ico2' // trim(suffix), grid_tracer_axes(1:2),                                       &
       model_time, 'Schmidt number - CO2' // trim(long_suffix), ' ',                           &
       missing_value = -1.0e+10)

  ibgc(n)%id_sc_o2 = register_diag_field(trim(diag_name),                                      &
       'sc_io2' // trim(suffix), grid_tracer_axes(1:2),                                        &
       model_time, 'Schmidt number - O2' // trim(long_suffix), ' ',                            &
       missing_value = -1.0e+10)

  ibgc(n)%id_alpha = register_diag_field(trim(diag_name),                                      &
       'ialpha' // trim(suffix), grid_tracer_axes(1:2),                                        &
       model_time, 'CO2* saturation' // trim(long_suffix), 'mol kg-1',                         &
       missing_value = -1.0e+10)

  ibgc(n)%id_csurf = register_diag_field(trim(diag_name),                                      &
       'icsurf' // trim(suffix), grid_tracer_axes(1:2),                                        &
       model_time, 'CO2* water' // trim(long_suffix), 'mol kg-1',                              &
       missing_value = -1.0e+10)

  ibgc(n)%id_pco2surf = register_diag_field(trim(diag_name),                                   &
       'ipco2surf' // trim(suffix), grid_tracer_axes(1:2),                                     &
       model_time, 'Oceanic pCO2' // trim(long_suffix), 'ppm',                                 &
       missing_value = -1.0e+10)

  ibgc(n)%id_sfc_flux_co2 = register_diag_field(trim(diag_name),                               &
       'sfc_flux_ico2' // trim(suffix), grid_tracer_axes(1:2),                                 &
       model_time, 'CO2 surface flux' // trim(long_suffix), 'mol m-2 s-1',                     &
       missing_value = -1.0e+10)

  ibgc(n)%id_htotal = register_diag_field(trim(diag_name),                                     &
       'ihtotal' // trim(suffix), grid_tracer_axes(1:2),                                       &
       model_time, 'H+ ion concentration' // trim(long_suffix), 'mol kg-1',                    &
       missing_value = -1.0e+10)

  ibgc(n)%id_alk = register_diag_field(trim(diag_name),                                        &
       'ialk' // trim(suffix), grid_tracer_axes(1:2),                                          &
       model_time, 'ALK' // trim(long_suffix), 'mol eq kg-1',                                  &
       missing_value = -1.0e+10)

  ibgc(n)%id_sfc_flux_o2 = register_diag_field(trim(diag_name),                                &
       'sfc_flux_io2' // trim(suffix), grid_tracer_axes(1:2),                                  &
       model_time, 'iO2 surface flux' // trim(long_suffix), 'mol m-2 s-1',                     &
       missing_value = -1.0e+10)

  ibgc(n)%id_jio2 = register_diag_field(trim(diag_name),                                       &
       'jio2' // trim(suffix), grid_tracer_axes(1:3),                                          &
       model_time, 'Total iO2 biological source, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

if (do_carbon_comp) then

  ibgc(n)%id_csatsurf = register_diag_field(trim(diag_name),                                   &
       'icsatsurf' // trim(suffix), grid_tracer_axes(1:2),                                     &
       model_time, 'CO2* saturation DIC' // trim(long_suffix), 'mol kg-1',                     &
       missing_value = -1.0e+10)

  ibgc(n)%id_htotal_sat = register_diag_field(trim(diag_name),                                 &
       'ihtotal_sat' // trim(suffix), grid_tracer_axes(1:2),                                   &
       model_time, 'H+ ion concentration for DIC_sat' // trim(long_suffix), 'mol kg-1',        &
       missing_value = -1.0e+10)

  ibgc(n)%id_pco2satsurf = register_diag_field(trim(diag_name),                                &
       'ipco2satsurf' // trim(suffix), grid_tracer_axes(1:2),                                  &
       model_time, 'Oceanic pCO2 saturation DIC' // trim(long_suffix), 'ppm',                  &
       missing_value = -1.0e+10)

  ibgc(n)%id_sfc_flux_co2sat = register_diag_field(trim(diag_name),                            &
       'sfc_flux_ico2sat' // trim(suffix), grid_tracer_axes(1:2),                              &
       model_time, 'Saturation CO2 surface flux' // trim(long_suffix), 'mol m-2 s-1',          &
       missing_value = -1.0e+10)

endif

if (do_radiocarbon) then

  ibgc(n)%id_c14surf = register_diag_field(trim(diag_name),                                    &
       'ic14surf' // trim(suffix), grid_tracer_axes(1:2),                                      &
       model_time, '14CO2* water' // trim(long_suffix), 'mol kg-1',                            &
       missing_value = -1.0e+10)

  ibgc(n)%id_sfc_flux_14co2 = register_diag_field(trim(diag_name),                             &
       'sfc_flux_i14co2' // trim(suffix), grid_tracer_axes(1:2),                               &
       model_time, '14CO2 surface flux' // trim(long_suffix), 'mol m-2 s-1',                   &
       missing_value = -1.0e+10)

  ibgc(n)%frac_14catm = register_diag_field(trim(diag_name),                                   &
       'frac_14catm' // trim(suffix), grid_tracer_axes(1:2),                                   &
       model_time, '14C fraction' // trim(long_suffix), '',                                    &
       missing_value = -1.0e+10)

  ibgc(n)%id_jdecay_idi14c = register_diag_field(trim(diag_name),                              &
       'jdecay_idi14c' // trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'iDI14C decay, layer integral' // trim(long_suffix), 'mol m-2 s-1',         &
       missing_value = -1.0e+10)

  ibgc(n)%id_jdecay_ido14c = register_diag_field(trim(diag_name),                              &
       'jdecay_ido14c' // trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'iDO14C decay, layer integral' // trim(long_suffix), 'mol m-2 s-1',         &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_pi14c = register_diag_field(trim(diag_name),                               &
       'jremin_pi14c' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'Remineralization of Pi14C, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_ido14c = register_diag_field(trim(diag_name),                              &
       'jremin_ido14c' // trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'Remineralization of iDO14C, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_pi14c = register_diag_field(trim(diag_name),                                &
       'jprod_pi14c' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Pi14C production, layer integral' // trim(long_suffix), 'mol m-2 s-1',     &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_ido14c = register_diag_field(trim(diag_name),                               &
       'jprod_ido14c' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'iDO14C production, layer integral' // trim(long_suffix), 'mol m-2 s-1',    &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpi14c = register_diag_field(trim(diag_name),                                     &
       'fpi14c' // trim(suffix), grid_tracer_axes(1:3),                                        &
       model_time, 'Particulate i14C sinking flux' // trim(long_suffix), 'mol m-2 s-1',        &
       missing_value = -1.0e+10)

  ibgc(n)%id_jidi14c = register_diag_field(trim(diag_name),                                    &
       'jidi14c' // trim(suffix), grid_tracer_axes(1:3),                                       &
       model_time, 'Total iDI14C source, layer integral' // trim(long_suffix), 'mol m-2 s-1',  &
       missing_value = -1.0e+10)

endif

endif

if (do_no3_iso) then

   ibgc(n)%id_jprod_pi15n = register_diag_field(trim(diag_name),                              &
       'jprod_pi15n' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'Particulate ideal 15N production, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_ido15n = register_diag_field(trim(diag_name),                               &
       'jprod_ido15n' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'Dissolved ideal 15N source, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_ido15n = register_diag_field(trim(diag_name),                              &
       'jremin_ido15n' // trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'Dissolved ideal 15N remineralization' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpi15n = register_diag_field(trim(diag_name),                                     &
       'fpi15n' // trim(suffix), grid_tracer_axes(1:3),                                        &
       model_time, 'Particulate ideal 15N sinking flux' // trim(long_suffix), 'mol m-2 s-1',   &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_pi15n = register_diag_field(trim(diag_name),                               &
       'jremin_pi15n' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'Ideal 15no3 remin fr particles, layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_ji15no3 = register_diag_field(trim(diag_name),                                    &
       'ji15no3' // trim(suffix), grid_tracer_axes(1:3),                                       &
       model_time, '15no3 net source, layer integral' // trim(long_suffix), 'mol m-2 s-1',     &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_i18o = register_diag_field(trim(diag_name),                                 &
       'jprod_i18o' // trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'N18O3 uptake source, layer integral' // trim(long_suffix), 'mol m-2 s-1',  &
       missing_value = -1.0e+10)

endif

if (do_isio4) then

  ibgc(n)%id_remin_isi = register_diag_field(trim(diag_name),                                  &
       'remin_isi' // trim(suffix), grid_tracer_axes(1:3),                                     &
       model_time, 'Remineralization inverse lengthscale for fpisi' // trim(long_suffix), 'm-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_pisi = register_diag_field(trim(diag_name),                                 &
       'jprod_pisi' // trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Particulate ideal Si production layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpisi = register_diag_field(trim(diag_name),                                      &
       'fpisi' // trim(suffix), grid_tracer_axes(1:3),                                         &
       model_time, 'Particulate ideal Si sinking flux' // trim(long_suffix), 'mol m-2 s-1',    &
       missing_value = -1.0e+10)

  ibgc(n)%id_jremin_pisi = register_diag_field(trim(diag_name),                                &
       'jremin_pisi' // trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Ideal SiO4 remineralization layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_jprod_pi30si = register_diag_field(trim(diag_name),                               &
       'jprod_pi30si' // trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'Particulate ideal 30Si production layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

  ibgc(n)%id_fpi30si = register_diag_field(trim(diag_name),                                    &
       'fpi30si' // trim(suffix), grid_tracer_axes(1:3),                                       &
       model_time, 'Particulate ideal 30Si sinking flux' // trim(long_suffix), 'mol m-2 s-1',  &
       missing_value = -1.0e+10)
       
  ibgc(n)%id_jremin_pi30si = register_diag_field(trim(diag_name),                              &
       'jremin_pi30si' // trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'Ideal 30SiO4 remineralization layer integral' // trim(long_suffix), 'mol m-2 s-1', &
       missing_value = -1.0e+10)

endif

enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), 'Tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocean_ibgc_start
! </SUBROUTINE> NAME="ocean_ibgc_start"


!#######################################################################
! <SUBROUTINE NAME="ocean_ibgc_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode.
!
! Here, it is used only to set the values of preformed iPO4 and iDIC in 
! the mixed layer.
!
! </DESCRIPTION>
!
subroutine ocean_ibgc_tracer(isc, iec, jsc, jec, isd, jsd, nk,   &
                                Time, t_prog,              &
                                depth_zt, hblt_depth)
integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: nk
type(ocean_time_type), intent(in)                          :: Time
type(ocean_prog_tracer_type), intent(inout), dimension(:)  :: t_prog
real, dimension(isd:,jsd:,:), intent(in)                   :: depth_zt
real, dimension(isd:,jsd:), intent(in)                     :: hblt_depth

integer :: i, j, k, n
integer :: taup1

!     set Preformed phosphate, DIC
!  Considered doing them differently, with DIC_pre set only in the very 
! surface layer, where gas exchange takes place...but for ease of 
! interpretation, made them equivalent.
taup1 = Time%taup1

if ((do_po4f) .or. (do_po4)) then

do n = 1, instances

if (do_bgc_felim) then
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
         if (depth_zt(i,j,k) <= hblt_depth(i,j)) then
           t_prog(ibgc(n)%ind_ipo4_pre)%field(i,j,k,taup1) =               &
             t_prog(ibgc(n)%ind_ipo4f)%field(i,j,k,taup1) 
         endif
  enddo; enddo ; enddo
else
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec
         if (depth_zt(i,j,k) <= hblt_depth(i,j)) then
           t_prog(ibgc(n)%ind_ipo4_pre)%field(i,j,k,taup1) =               &
             t_prog(ibgc(n)%ind_ipo4)%field(i,j,k,taup1) 
         endif
  enddo; enddo ; enddo
endif

enddo

endif

if (do_carbon_comp) then

do n = 1, instances

  do j = jsc, jec
    do i = isc, iec
      do k = 1,nk
         if (depth_zt(i,j,k) <= hblt_depth(i,j)) then
           t_prog(ibgc(n)%ind_idic_pre)%field(i,j,k,taup1) =               &
             t_prog(ibgc(n)%ind_idic)%field(i,j,k,taup1) 
         endif
      enddo
    enddo
  enddo

enddo

endif

return

end subroutine  ocean_ibgc_tracer
! </SUBROUTINE> NAME="ocean_ibgc_tracer"

!#######################################################################
end module  ocean_ibgc_mod
