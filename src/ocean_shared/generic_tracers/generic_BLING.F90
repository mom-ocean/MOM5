!----------------------------------------------------------------
! <CONTACT EMAIL="Eric.Galbraith@mcgill.ca"> Eric D. Galbraith
! </CONTACT>
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
! </CONTACT>
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Anand Gnanandesikan
! </CONTACT>
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Rick Slater
! </REVIEWER>
!
! <OVERVIEW>
!   This module contains the generic version of BLING.
!   It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
!
!   WARNING: although the core components of the model (PO4, Fed, DOP, O2) 
!   have been quite well tested, the other components should be viewed as
!   developmental at this point. There may still be some bugs, particularly
!   in the CaCO3 burial scheme. EDG June 4, 2009
!
! </OVERVIEW>
!<DESCRIPTION>
!   Biogeochemistry with Light, Iron, Nutrient and Gas (BLING) includes an
!   implicit ecological model of growth limitation by light,
!   temperature, phosphate and iron, along with dissolved organic
!   phosphorus and O2 pools.
!   Food web processing in the euphotic zone and remineralization/
!   dissolution through the ocean interior are handled as in Dunne et al. 
!   (2005).  O2 equilibria and gas exchange follow OCMIP2 protocols.
!   Additional functionality comes from an optional carbon cycle that is 
!   non-interactive, i.e. does not change the core BLING behaviour, as
!   well as tracers for radiocarbon (14c), a decomposition of carbon 
!   components by gas exchange and remineralization (carbon_pre), and a 
!   decomposition of phosphate as preformed and remineralized (po4_pre).
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! This model is available for public use. 
! The current version is BLING.0. The version number refers to the core
! model behaviour; additional tracers exist in different iterations of the
! module. In publications it should be referenced as:
! Galbraith, E.D., Gnanadesikan, A., Dunne, J. and Hiscock, M. 2009.
! Regional impacts of iron-light colimitation in a global 
! biogeochemical model. Biogeosciences Discussions, 6, 1-47.
!
! All parameter values are as described in this paper.
! Note that this reference is only for the core model components, and
! does not include any of the additional functionalities, which remain
! undocumented. Please contact Eric Galbraith (eric.galbraith@mcgill.ca)
! for more information.
! </REFERENCE>
!
! <DEVELOPER_NOTES>
! This code was developed based on the template of Perth generic TOPAZ code.
! </DEVELOPER_NOTES>
! </INFO>
!
!<NAMELIST NAME="generic_bling_nml">
!
!  <DATA NAME="do_14c" TYPE="logical">
!  If true, then simulate radiocarbon. Includes 2 prognostic tracers, DI14C
! and DO14C. Requires that do_carbon = .true. Note that 14C is not taken up
! by CaCO3 at the current time, but cycles only through the soft tissue.
! This is a mistake that will be fixed later.
!  </DATA> 
!
!  <DATA NAME="do_carbon" TYPE="logical">
!  If true, then simulate the carbon cycle based on strict stoichiometry
! of C:P. Includes 1 prognostic tracer, DIC.
!  </DATA> 
!
!  <DATA NAME="do_carbon_pre" TYPE="logical">
!  If true, then simulate the carbon cycle based on strict stoichiometry
!  of C:P. Includes 3 prognostic tracers, DIC, ALK and ALK_pre. Requires 
! that do_carbon = .true.
!  </DATA> 
!
!  <DATA NAME="do_po4_pre" TYPE="logical">
!  If true, then simulate preformed PO4, a useful theoretical construct
! equal to PO4 in the surface layer and subject only to passive transport 
! everywhere else. 1 prognostic tracer, PO4_pre.
!  </DATA> 
!
!  <DATA NAME="bury_caco3" TYPE="logical">
!  If true, then allow CaCO3 to be buried in sediments as a function
! of sinking CaCO3 flux and bottom water saturation state, and allow
! a river input of alkalinity to compensate. Caution: this will cause
! the alkalinity to have a long term drift, which will produce a long
! term drift in the carbon cycle. Should be considered highly
! experimental at this point. Requires that do_carbon=.true.
!  </DATA> 
!
!  <DATA NAME="use_bling_chl" TYPE="logical">
!  If true, names the BLING diagnostic chlorophyll field 'chl'. Thus,
! if read_chl=.false. in the shortwave namelist, BLING chlorophyll will
! be used to calculate sw absorption. Set this to false in order to 
! have the chlorophyll field named chl_b, preventing conflict with the 
! TOPAZ chlorophyll field, so that the two modules can run
! simultaneously.
!  </DATA> 
!
!</NAMELIST>
!
!----------------------------------------------------------------

module generic_BLING

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len, fm_path_name_len
  use mpp_mod,           only: input_nml_file, mpp_error, stdlog, NOTE, WARNING, FATAL, stdout, mpp_chksum
  use fms_mod,           only: write_version_number, open_namelist_file, check_nml_error, close_file
  use fms_mod,           only: field_exist, file_exist
  use time_manager_mod,  only: time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  
  use diag_manager_mod,  only: register_diag_field, send_data 
  use constants_mod,     only: WTMCO2, WTMO2

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer
  use g_tracer_utils, only : g_tracer_get_common,g_tracer_set_common 
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector

  implicit none ; private

  character(len=128) :: version = '$Id: generic_BLING.F90,v 20.0 2013/12/14 00:18:02 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  character(len=fm_string_len), parameter :: mod_name       = 'generic_BLING'
  character(len=fm_string_len), parameter :: package_name   = 'generic_bling'

  public do_generic_BLING
  public generic_BLING_register
  public generic_BLING_init
  public generic_BLING_register_diag
  public generic_BLING_update_from_coupler
  public generic_BLING_update_from_source
  public generic_BLING_update_from_bottom
  public generic_BLING_set_boundary_values
  public generic_BLING_end

  !The following logical for using this module is overwritten 
  ! generic_tracer_nml namelist
  logical, save :: do_generic_BLING = .false.

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30

! Namelist Options

  logical :: do_14c             = .true.  ! Requires do_carbon = .true.
  logical :: do_carbon          = .true.  
  logical :: do_carbon_pre      = .true.  ! Requires do_carbon = .true.
  logical :: do_po4_pre         = .true.
  logical :: bury_caco3         = .false. ! Requires do_carbon = .true.
  logical :: use_bling_chl      = .false. ! Must be false to run with TOPAZ

namelist /generic_bling_nml/ do_14c, do_carbon, do_carbon_pre, &
  do_po4_pre, bury_caco3, use_bling_chl

  !
  !The following two types contain all the parameters and arrays used in this module.

  type generic_BLING_type

     logical  ::       &
          init,             &                  ! If tracers should be initializated
          prevent_neg_o2,   &
          tracer_debug

     real  ::          &
          alpha_max,        &                  ! Quantum yield under low light, Fe-replete
          alpha_min,        &                  ! Quantum yield under low light, Fe-limited
          c_2_p,            &                  ! Carbon to Phosphorus ratio
          ca_2_p,           &                  ! CaCO3 to Phosphorus ratio (of small phytoplankton)
          ca_remin_depth,   &                  ! CaCO3 dissolution length scale (subject to omega)
          chl_min,          &                  ! Minimum chl concentration allowed (for numerical stability)
          def_fe_min,       &                  ! Minimum value for iron deficiency term
          fast_gasex,       &                  ! Gas exchange rate multiplier for DIC_sat
          fe_2_p_max,       &                  ! Iron to Phosphate uptake ratio scaling
          fe_2_p_sed,       &                  ! Iron to Phosphorus ratio in sediments
          felig_bkg,        &                  ! Iron ligand concentration
          gamma_biomass,    &                  ! Biomass adjustment timescale
          gamma_dop,        &                  ! Dissolved organic phosphorus decay
          gamma_irr_mem,    &                  ! Photoadaptation timescale
          gamma_pop,        &                  ! Patriculate Organic Phosphorus decay
          gamma_tag,        &                  ! Restoring time constant for nutrient source tracers
          half_life_14c,    &                  ! Radiocarbon half-life
          k_fe_2_p,         &                  ! Fe:P half-saturation constant
          k_fe_uptake,      &                  ! Iron half-saturation concentration
          k_o2,             &                  ! Oxygen half-saturation concentration
          k_po4,            &                  ! Phosphate half-saturation concentration
          kappa_eppley,     &                  ! Temperature dependence
          kappa_remin,      &                  ! Temperature dependence for particle fractionation
          kfe_inorg,        &                  ! Iron scavenging, 2nd order
          kfe_eq_lig_max,   &                  ! Maximum light-dependent iron ligand stability constant
          kfe_eq_lig_min,   &                  ! Minimum light-dependent iron ligand stability constant
          kfe_eq_lig_irr,   &                  ! Irradiance scaling for iron ligand stability constant
          kfe_eq_lig_femin, &                  ! Low-iron threshold for ligand stability constant
          kfe_org,          &                  ! Iron scavenging, 1st order
          lambda0,          &                  ! Total mortality rate constant
          lambda_14c,       &                  ! Radiocarbon decay rate
          resp_frac,        &                  ! Biomass maintenace requirement as fraction of pc_0
          mass_2_p,         &                  ! Organic matter mass to Phosphorus ratio
          n_2_p,            &                  ! Nitrogen to Phosphorus ratio
          o2_2_p,           &                  ! Oxygen to Phosphorus ratio
          o2_min,           &                  ! Anaerobic respiration threshold
          P_star,           &                  ! Pivotal phytoplankton concentration
          pc_0,             &                  ! Maximum carbon-specific growth rate
          phi_dop,          &                  ! Dissolved organic phosphorus fraction of uptake
          phi_lg,           &                  ! Fraction of small phytoplankton converted to detritus
          phi_sm,           &                  ! Fraction of large phytoplankton converted to detritus
          remin_min,        &                  ! Minimum remineralization under low O2
          rho_dense,        &                  ! Deep boundary density for exposure tracers
          rho_light,        &                  ! Shallow boundary density for exposure tracers
          sed_flux,         &                  ! Seafloor sedimentary flux (globally constant)
          thetamax_hi,      &                  ! Maximum Chl:C ratio when iron-replete
          thetamax_lo,      &                  ! Maximum Chl:C ratio when iron-limited
          wsink_acc,        &                  ! Sinking rate acceleration with depth
          wsink0,           &                  ! Sinking rate at surface
          wsink0_z,         &                  ! Depth to which sinking rate remains constant
          z_sed                                ! Thickness of active sediment layer

     real    :: htotal_scale_lo, htotal_scale_hi, htotal_in
     real    :: Rho_0, a_0, a_1, a_2, a_3, a_4, a_5, b_0, b_1, b_2, b_3, c_0
     real    :: a1_co2, a2_co2, a3_co2, a4_co2, a1_o2, a2_o2, a3_o2, a4_o2

!
! The prefixes "f_" refers to a "field" and "j" to a volumetric rate, and 
! "b_" to a bottom flux. The prefix "p_" refers to a "pointer".
!
     real, dimension(:,:,:), ALLOCATABLE ::  &
          alpha,&
          biomass_p_ts,&
          def_fe,&
          expkT,&
          f_biomass_p,&
          f_cased,&
          f_chl,&
          f_dop,&
          f_fed,&
          f_htotal,&
          f_htotal_sat,&
          f_irr_mem,&
          f_o2,&
          f_po4,&
          f_po4_pre,&
          fe_2_p_uptake,&
          feprime,&
          fpofe,&
          fpop,&
          frac_lg,&
          frac_pop,&
          irr_inst,&
          irr_mix,&
          irrk,&
          jdop,&
          jfe_ads_inorg,&
          jfe_ads_org,&
          jfe_recycle,&
          jfe_reminp,&
          jfe_uptake,&
          jo2,&
          jp_recycle,&
          jp_reminp,&
          jp_uptake,&
          jpo4,&
          jpofe,&
          jpop,&
          kfe_eq_lig,&
          pc_m,&
          mu,&
          pc_tot,&
          theta,&
          thetamax_fe,&
          wsink,&
          zremin,&
          zt

     real, dimension(:,:,:), ALLOCATABLE ::  &
          f_dop_n,&
          f_dop_nx,&
          f_dop_s,&
          f_dop_sx,&
          f_po4_n,&
          f_po4_nx,&
          f_po4_s,&
          f_po4_sx,&
          f_po4_pre_n,&
          f_po4_pre_s,&
          fpop_n,&
          fpop_nx,&
          fpop_s,&
          fpop_sx,&
          jpn_reminp,&
          jpnx_reminp,&
          jps_reminp,&
          jpsx_reminp

     real, dimension(:,:,:), ALLOCATABLE ::  &
          c14_2_p,&
          co3_solubility,&
          f_alk,&
          f_alk_pre,&
          f_co3_ion,&
          f_di14c,&
          f_dic,&
          f_dic_pre,&
          f_dic_sat,&
          f_do14c,&
          fcaco3,&
          fpo14c,&
          j14c_decay_dic,&
          j14c_decay_doc,&
          j14c_reminp,&
          jca_reminp,&
          jca_uptake,&
          jdi14c,&
          jdo14c,&
          zremin_caco3

     real, dimension(:,:), ALLOCATABLE :: &
          fe_burial,&
          ffe_sed,&
          b_fed,&
          b_o2,&
          b_po4

     real, dimension(:,:), ALLOCATABLE :: &
          b_po4_n,&
          b_po4_nx,&
          b_po4_s,&
          b_po4_sx

     real, dimension(:,:), ALLOCATABLE ::  &
          co2_csurf,pco2_surf,co2_alpha,&
          c14o2_csurf,c14o2_alpha,co2_sat_csurf,pco2_sat_surf,&
          htotallo, htotalhi,&
          htotal_satlo, htotal_sathi,&
          b_alk,&
          b_di14c,&
          b_dic,&
          fcaco3_to_sed,&
          fcased_burial,&
          fcased_redis
     
     real, dimension(:,:,:,:), pointer :: &
          p_dop,&
          p_fed,&
          p_o2,&
          p_po4,&
          p_po4_pre

     real, dimension(:,:,:,:), pointer :: &
          p_dop_n,&
          p_dop_nx,&
          p_dop_s,&
          p_dop_sx,&
          p_po4_n,&
          p_po4_nx,&
          p_po4_s,&
          p_po4_sx,&
          p_po4_pre_n,&
          p_po4_pre_s

     real, dimension(:,:,:,:), pointer :: &
          p_alk,&
          p_alk_pre,&
          p_di14c,&
          p_dic,&
          p_dic_pre,&
          p_do14c
     
     integer :: nkml
     character(len=fm_string_len) :: file
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file

     integer                   :: &
      id_alpha         = -1,  & ! Iron-limited initial slope of P-I curve 
      id_b_alk         = -1,  & ! Bottom flux of alkalinity
      id_b_dic         = -1,  & ! Bottom flux of DIC
      id_b_di14c       = -1,  & ! Bottom flux of DI14C
      id_b_fed         = -1,  & ! Bottom flux of Fe
      id_b_o2          = -1,  & ! Bottom flux of O2
      id_b_po4         = -1,  & ! Bottom flux of PO4
      id_b_po4_n       = -1,  & ! Bottom flux of PO4_n
      id_b_po4_nx      = -1,  & ! Bottom flux of PO4_nx
      id_b_po4_s       = -1,  & ! Bottom flux of PO4_s
      id_b_po4_sx      = -1,  & ! Bottom flux of PO4_sx
      id_biomass_p_ts  = -1,  & ! Instantaneous P concentration in biomass
      id_c14_2_p       = -1,  & ! DI14C to PO4 uptake ratio
      id_c14o2_csurf   = -1,  & ! Surface water 14CO2*
      id_c14o2_alpha   = -1,  & ! Surface water 14CO2* solubility 
      id_co2_csurf     = -1,  & ! Surface water CO2*
      id_co2_alpha     = -1,  & ! Surface water CO2* solubility
      id_co2_sat_csurf = -1,  & ! Surface water CO2* for DIC_sat
      id_co3_solubility= -1,  & ! Calcite solubility
      id_def_fe        = -1,  & ! Iron deficiency term 
      id_expkT         = -1,  & ! Temperature dependence 
      id_fcaco3        = -1,  & ! CaCO3 sinking flux
      id_fcaco3_to_sed = -1,  & ! CaCO3 sinking flux in bottom layer
      id_fcased_burial = -1,  & ! CaCO3 permanent burial flux
      id_fcased_redis  = -1,  & ! CaCO3 dissolution flux from active sediment layer
      id_fe_2_p_uptake = -1,  & ! Fed:PO4 of instantaneous uptake
      id_feprime       = -1,  & ! Free (unbound) iron concentration
      id_fe_burial     = -1,  & ! Flux of iron to sediment as particulate
      id_ffe_sed       = -1,  & ! Sediment iron efflux
      id_fpofe         = -1,  & ! POFe sinking flux
      id_fpo14c        = -1,  & ! PO14C sinking flux
      id_fpop          = -1,  & ! POP sinking flux
      id_fpop_n        = -1,  & ! POP_N sinking flux
      id_fpop_nx       = -1,  & ! POP_Nx sinking flux
      id_fpop_s        = -1,  & ! POP_S sinking flux
      id_fpop_sx       = -1,  & ! POP_Sx sinking flux
      id_frac_lg       = -1,  & ! Fraction of production by large phytoplankton
      id_frac_pop      = -1,  & ! Fraction of uptake converted to particulate
      id_irr_inst      = -1,  & ! Instantaneous irradiance 
      id_irr_mix       = -1,  & ! Mixed layer irradiance 
      id_irrk          = -1,  & ! Effective susceptibility to light limitation 
      id_j14c_decay_dic= -1,  & ! Radioactive decay of DI14C
      id_j14c_decay_doc= -1,  & ! Radioactive decay of DO14C
      id_j14c_reminp   = -1,  & ! 14C particle remineralization layer integral
      id_jca_reminp    = -1,  & ! CaCO3 dissolution layer integral
      id_jca_uptake    = -1,  & ! CaCO3 formation layer integral
      id_jdi14c        = -1,  & ! DI14C source layer integral
      id_jdo14c        = -1,  & ! Semilabile DO14C source layer integral
      id_jdop          = -1,  & ! Semilabile DOP source layer integral
      id_jfe_ads_inorg = -1,  & ! Iron adsorption (2nd order) layer integral
      id_jfe_ads_org   = -1,  & ! Iron adsorption to fpop layer integral
      id_jfe_recycle   = -1,  & ! Iron fast recycling layer integral
      id_jfe_reminp    = -1,  & ! Iron particle remineralization layer integral
      id_jfe_uptake    = -1,  & ! Iron uptake layer integral
      id_jo2           = -1,  & ! O2 source layer integral
      id_jp_recycle    = -1,  & ! Phosphorus fast recycling layer integral
      id_jp_reminp     = -1,  & ! Phosphorus particle remineralization layer integral
      id_jpn_reminp    = -1,  & ! Northern phosphorus particle remineralization layer integral
      id_jpnx_reminp   = -1,  & ! Northern exposure phosphorus particle remineralization layer integral
      id_jps_reminp    = -1,  & ! Southern phosphorus particle remineralization layer integral
      id_jpsx_reminp   = -1,  & ! Southern exposure phosphorus particle remineralization layer integral
      id_jp_uptake     = -1,  & ! Phosphorus uptake layer integral
      id_jpo4          = -1,  & ! PO4 source layer integral
      id_jpofe         = -1,  & ! Particulate organic iron source layer integral
      id_jpop          = -1,  & ! Particulate organic phosphorus source layer integral
      id_kfe_eq_lig    = -1,  & ! Iron-ligand stability constant
      id_pc_m          = -1,  & ! Light-saturated maximum photosynthesis rate (carbon specific)
      id_mu            = -1,  & ! Growth rate after respiratory loss(carbon specific)
      id_pc_tot        = -1,  & ! Fully-limited photosynthesis rate (carbon specific)
      id_pco2_surf     = -1,  & ! Surface water pCO2
      id_pco2_sat_surf = -1,  & ! Surface water pCO2 for DIC_sat
      id_theta         = -1,  & ! Chl:C ratio
      id_thetamax_fe   = -1,  & ! Iron-limited maximum Chl:C ratio
      id_wsink         = -1,  & ! Sinking rate
      id_zremin        = -1,  & ! Remineralization length scale
      id_zremin_caco3  = -1,  & ! CaCO3 remineralization length scale
      id_alk           = -1,  & ! Alkalinity Prognostic tracer
      id_alk_pre       = -1,  & ! Preformed Alkalinity Prognostic tracer
      id_di14c         = -1,  & ! Dissolved inorganic radiocarbon Prognostic tracer
      id_dic           = -1,  & ! Dissolved inorganic carbon Prognostic tracer
      id_dic_pre       = -1,  & ! Preformed DIC Prognostic tracer
      id_dic_sat       = -1,  & ! Saturation DIC Prognostic tracer
      id_do14c         = -1,  & ! Semi-labile DO radiocarbon Prognostic tracer
      id_dop           = -1,  & ! Semi-labile dissolved organic P Prognostic tracer
      id_dop_n         = -1,  & ! Semi-labile dissolved organic P Prognostic tracer
      id_dop_nx         = -1,  & ! Semi-labile dissolved organic P Prognostic tracer
      id_dop_s         = -1,  & ! Semi-labile dissolved organic P Prognostic tracer
      id_dop_sx         = -1,  & ! Semi-labile dissolved organic P Prognostic tracer
      id_fed           = -1,  & ! Dissolved Iron Prognostic tracer
      id_o2            = -1,  & ! Oxygen Prognostic tracer
      id_po4           = -1,  & ! Phosphate Prognostic tracer
      id_po4_n         = -1,  & ! Phosphate Prognostic tracer
      id_po4_nx        = -1,  & ! Phosphate Prognostic tracer
      id_po4_s         = -1,  & ! Phosphate Prognostic tracer
      id_po4_sx        = -1,  & ! Phosphate Prognostic tracer
      id_po4_pre       = -1,  & ! Preformed Phosphate Prognostic tracer
      id_po4_pre_n     = -1,  & ! Preformed Phosphate Prognostic tracer
      id_po4_pre_s     = -1,  & ! Preformed Phosphate Prognostic tracer
      id_htotal        = -1,  & ! Hydrogen ion Diagnostic tracer
      id_htotal_sat    = -1,  & ! Hydrogen ion for DIC_sat Diagnostic tracer
      id_co3_ion       = -1,  & ! CO3= ion Diagnostic tracer
      id_cased         = -1,  & ! Active sediment CaCO3 concentration Diagnostic tracer
      id_chl           = -1,  & ! Chlorophyll Diagnostic tracer
      id_biomass_p     = -1,  & ! Biomass Diagnostic tracer
      id_irr_mem       = -1     ! Irradiance Memory Diagnostic tracer

  end type generic_BLING_type

  !An auxiliary type for storing varible names
  type, public :: vardesc
     character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
     character(len=fm_string_len) :: longname ! The long name of that variable.
     character(len=1)  :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
     character(len=1)  :: z_grid   ! The vert. grid:  L, i, or 1.
     character(len=1)  :: t_grid   ! The time description: s, a, m, or 1.
     character(len=fm_string_len) :: units    ! The dimensions of the variable.
     character(len=1)  :: mem_size ! The size in memory: d or f.
  end type vardesc

  type(generic_BLING_type), save :: bling

  type(CO2_dope_vector) :: CO2_dope_vec

  
contains

!#######################################################################
  subroutine generic_BLING_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

integer                                                 :: ioun
integer                                                 :: ierr
integer                                                 :: io_status
character(len=fm_string_len)                            :: name
integer                                                 :: stdoutunit,stdlogunit
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!
character(len=fm_string_len), parameter :: sub_name = 'generic_bling_register'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '


! provide for namelist over-ride
! This needs to go before the add_tracers in order to allow the namelist 
! settings to switch tracers on and off.
!
stdoutunit=stdout();stdlogunit=stdlog()

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=generic_bling_nml, iostat=io_status)
ierr = check_nml_error(io_status,'generic_bling_nml')
#else
ioun = open_namelist_file()
read  (ioun, generic_bling_nml,iostat=io_status)
ierr = check_nml_error(io_status,'generic_bling_nml')
call close_file (ioun)
#endif

write (stdoutunit,'(/)')
write (stdoutunit, generic_bling_nml)
write (stdlogunit, generic_bling_nml)

  if ((do_14c) .and. (do_carbon)) then
    write (stdoutunit,*) trim(note_header), 'Simulating radiocarbon'
  else if ((do_14c) .and. .not. (do_carbon)) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_14c requires do_carbon' // trim(name))
  endif

  if ((do_carbon_pre) .and. (do_carbon)) then
    write (stdoutunit,*) trim(note_header), 'Calculating DIC_pre and DIC_sat'
  else if ((do_carbon_pre) .and. .not. (do_carbon)) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Do_carbon_pre requires do_carbon' // trim(name))
  endif

  if ((bury_caco3) .and. (do_carbon)) then
    write (stdoutunit,*) trim(note_header), &
         'CAUTION: burying CaCO3, are you sure you want to do this?'
  else if ((bury_caco3) .and. .not. (do_carbon)) then
    call mpp_error(FATAL, trim(error_header) //        &
         'Bury_caco3 requires do_carbon' // trim(name))
  endif

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

    
    
  end subroutine generic_BLING_register

!#######################################################################
  ! <SUBROUTINE NAME="generic_BLING_init">
  !  <OVERVIEW>
  !   Initialize the generic BLING module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the BLING Tracers to the list of generic Tracers passed 
  !  to it via utility subroutine g_tracer_add(). Adds all the parameters 
  !  used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_BLING_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_BLING_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    call write_version_number( version, tagname )

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate all the private work arrays used by this module.
    call user_allocate_arrays

  end subroutine generic_BLING_init

!#######################################################################
  !   Register diagnostic fields to be used in this module. 
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !
  subroutine generic_BLING_register_diag
    real,parameter :: missing_value1=-1.0e+10
    type(vardesc)  :: vardesc_temp
    integer        :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, axes(3)
    type(time_type):: init_time 

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes=axes,init_time=init_time) 


    !   The following vardesc types contain a package of metadata about each tracer,
    ! including, in order, the following elements: name; longname; horizontal
    ! staggering ('h') for collocation with thickness points ; vertical staggering
    ! ('L') for a layer variable ; temporal staggering ('s' for snapshot) ; units ;
    ! and precision in non-restart output files ('f' for 32-bit float or 'd' for
    ! 64-bit doubles). For most tracers, only the name, longname and units should
    ! be changed.  

    !
    ! Register Diagnostics
    !===========================================================
    !
    ! Core diagnostics

    vardesc_temp = vardesc&
    ("alpha","Fe-limitated initial slope of P-I curve",'h','L','s','g C g Chl-1 m2 W-1 s-1','f')
    bling%id_alpha = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("b_fed","Bottom flux of Fe into sediment",'h','1','s','mol m-2 s-1','f')
    bling%id_b_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("b_o2","Bottom flux of O2 into sediment",'h','1','s','mol m-2 s-1','f')
    bling%id_b_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("b_po4","Bottom flux of PO4 into sediment",'h','1','s','mol m-2 s-1','f')
    bling%id_b_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("biomass_p_ts","Instantaneous P concentration in biomass",'h','L','s','mol kg-1','f')
    bling%id_biomass_p_ts = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("def_Fe","Iron deficiency term",'h','L','s','unitless','f')
    bling%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("expkT","Temperature dependence",'h','L','s','unitless','f')
    bling%id_expkT = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fe_2_p_uptake","Uptake ratio of Fed:PO4",'h','L','s','mol Fe mol P-1','f')
    bling%id_fe_2_p_uptake = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fe_burial","Sedimenting iron flux",'h','1','s','mol m-2 s-1','f')
    bling%id_fe_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("feprime","Concentration of free, unbound iron",'h','L','s','mol kg-1','f')
    bling%id_feprime = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("ffe_sed","Sediment iron efflux",'h','1','s','mol m-2 s-1','f')
    bling%id_ffe_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fpofe","POFe sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
    bling%id_fpofe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fpop","POP sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
    bling%id_fpop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("frac_lg","Fraction of production by large phytoplankton",'h','L','s','unitless','f')
    bling%id_frac_lg = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("frac_pop","Particulate fraction of total uptake",'h','L','s','unitless','f')
    bling%id_frac_pop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("irr_inst","Instantaneous light",'h','L','s','W m-2','f')
    bling%id_irr_inst = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("irr_mix","Mixed layer light",'h','L','s','W m-2','f')
    bling%id_irr_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("irrk","Tendency to light limitation",'h','L','s','W m-2','f')
    bling%id_irrk = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jdop","DOP source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jfe_ads_inorg","Iron adsorption (2nd order) layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jfe_ads_inorg = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jfe_ads_org","Iron adsorption to FPOP layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jfe_ads_org = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jfe_recycle","Fast recycling of iron layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jfe_recycle = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jfe_reminp","Sinking particulate Fe decay layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jfe_reminp = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jfe_uptake","Iron production layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jfe_uptake = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jp_recycle","Fast recycling of PO4 layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jp_recycle = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jp_reminp","Sinking particulate P decay layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jp_reminp = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jp_uptake","PO4 uptake layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jp_uptake = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jo2","O2 source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jpo4","PO4 source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jpo4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jpop","Particulate P source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jpop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jpofe","Particulate Fe source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jpofe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("kfe_eq_lig","Iron ligand stability constant",'h','L','s','mol-1 kg','f')
    bling%id_kfe_eq_lig = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("pc_m","Light-saturated photosynthesis rate (carbon specific)",'h','L','s','s-1','f')
    bling%id_pc_m = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("mu","Net growth rate after respiratory loss",'h','L','s','s-1','f')
    bling%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("pc_tot","Photosynthesis rate (carbon specific)",'h','L','s','s-1','f')
    bling%id_pc_tot = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("theta","Chl:C ratio",'h','L','s','g Chl g C-1','f')
    bling%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("thetamax_fe","Fe-limited max Chl:C",'h','L','s','g Chl g C-1','f')
    bling%id_thetamax_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("wsink","Sinking rate",'h','L','s','m s-1','f')
    bling%id_wsink = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("zremin","Remineralization lengthscale",'h','L','s','m','f')
    bling%id_zremin = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    if (do_carbon) then                                      !<<CARBON CYCLE
    vardesc_temp = vardesc&
    ("b_alk","Bottom flux of Alk into sediment",'h','1','s','mol m-2 s-1','f')
    bling%id_b_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("b_dic","Bottom flux of DIC into sediment",'h','1','s','mol m-2 s-1','f')
    bling%id_b_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("co2_alpha","Saturation surface CO2* per uatm",'h','1','s','mol kg-1 atm-1','f')
    bling%id_co2_alpha = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("co2_csurf","CO2* concentration at surface",'h','1','s','mol kg-1','f')
    bling%id_co2_csurf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("co3_solubility","Calcite solubility",'h','L','s','mol kg-1','f')
    bling%id_co3_solubility = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fcaco3","CaCO3 sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
    bling%id_fcaco3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fcaco3_to_sed","CaCO3 sinking flux at ocean bottom",'h','1','s','mol m-2 s-1','f')
    bling%id_fcaco3_to_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jca_reminp","Sinking CaCO3 dissolution layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jca_reminp = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jca_uptake","CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jca_uptake = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("pco2_surf","Seawater pCO2 in surface layer",'h','1','s','uatm','f')
    bling%id_pco2_surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("zremin_caco3","CaCO3 Remineralization lengthscale",'h','L','s','m','f')
    bling%id_zremin_caco3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

      if (bury_caco3) then                                     !<<BURY CACO3  
    vardesc_temp = vardesc&
    ("fcased_burial","CaCO3 permanent burial flux",'h','1','s','mol m-2 s-1','f')
    bling%id_fcased_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fcased_redis","CaCO3 redissolution from sediments",'h','1','s','mol m-2 s-1','f')
    bling%id_fcased_redis = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      endif                                                    !BURY CACO3>>
         
      if (do_14c) then                                        !<<RADIOCARBON
    vardesc_temp = vardesc&
    ("b_di14c","Bottom flux of DI14C into sediment",'h','1','s','mol m-2 s-1','f')
    bling%id_b_di14c = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("c14_2_p","Ratio of DI14C to PO4",'h','L','s','mol kg-1','f')
    bling%id_c14_2_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("c14o2_alpha","Saturation surface 14CO2* per uatm",'h','1','s','mol kg-1 atm-1','f')
    bling%id_c14o2_alpha = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("c14o2_csurf","14CO2* concentration at surface",'h','1','s','mol kg-1','f')
    bling%id_c14o2_csurf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fpo14c","PO14C sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
    bling%id_fpo14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("j14c_decay_dic","DI14C radioactive decay layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_j14c_decay_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("j14c_decay_doc","DO14C radioactive decay layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_j14c_decay_doc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("j14c_reminp","Sinking PO14C remineralization layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_j14c_reminp = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jdi14c","DI14C source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jdi14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jdo14c","DO14C source layer integral",'h','L','s','mol m-2 s-1','f')
    bling%id_jdo14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      endif                                                   !RADIOCARBON>>

      if (do_carbon_pre) then                                     !<<DIC_PRE  
    vardesc_temp = vardesc&
    ("co2_sat_csurf","CO2* concentration for DIC_sat at surface",'h','1','s','mol kg-1','f')
    bling%id_co2_sat_csurf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("pco2_sat_surf","Seawater pCO2 in surface layer",'h','1','s','uatm','f')
    bling%id_pco2_sat_surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      endif                                                       !DIC_PRE>>

    endif                                                    !CARBON CYCLE>>

  end subroutine generic_BLING_register_diag

!#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_BLING_params type
    !==============================================================    

    !=============
    !Block Starts: g_tracer_add_param
    !=============
    !Add the known experimental parameters used for calculations
    !in this module.
    !All the g_tracer_add_param calls must happen between 
    !g_tracer_start_param_list and g_tracer_end_param_list  calls.
    !This implementation enables runtime overwrite via field_table.
    !
    call g_tracer_start_param_list(package_name)

    !     g_tracer_add_param(name   , variable   ,  default_value)
    !call g_tracer_add_param('', bling%,  )
    !
    call g_tracer_add_param('init', bling%init, .false. )
    !
    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    call g_tracer_add_param('RHO_0', bling%Rho_0, 1035.0)
    call g_tracer_add_param('NKML' , bling%nkml, 1)
    !
    !-----------------------------------------------------------------------
    ! Gas exchange
    !-----------------------------------------------------------------------
    !       coefficients for O2 saturation
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_0', bling%a_0, 2.00907)
    call g_tracer_add_param('a_1', bling%a_1, 3.22014)
    call g_tracer_add_param('a_2', bling%a_2, 4.05010)
    call g_tracer_add_param('a_3', bling%a_3, 4.94457)
    call g_tracer_add_param('a_4', bling%a_4, -2.56847e-01)
    call g_tracer_add_param('a_5', bling%a_5, 3.88767)
    call g_tracer_add_param('b_0', bling%b_0, -6.24523e-03)
    call g_tracer_add_param('b_1', bling%b_1, -7.37614e-03)
    call g_tracer_add_param('b_2', bling%b_2, -1.03410e-02 )
    call g_tracer_add_param('b_3', bling%b_3, -8.17083e-03)
    call g_tracer_add_param('c_0', bling%c_0, -4.88682e-07)
    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !-----------------------------------------------------------------------
    !  Compute the Schmidt number of CO2 in seawater using the 
    !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    !  7373-7382).
    !-----------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_co2', bling%a1_co2,  2068.9)
    call g_tracer_add_param('a2_co2', bling%a2_co2, -118.63)
    call g_tracer_add_param('a3_co2', bling%a3_co2,  2.9311)
    call g_tracer_add_param('a4_co2', bling%a4_co2, -0.027)
    !---------------------------------------------------------------------
    !  Compute the Schmidt number of O2 in seawater using the 
    !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
    !  Cycles, 12, 141-163).
    !---------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_o2', bling%a1_o2, 1929.7)
    call g_tracer_add_param('a2_o2', bling%a2_o2, -117.46)
    call g_tracer_add_param('a3_o2', bling%a3_o2, 3.116)
    call g_tracer_add_param('a4_o2', bling%a4_o2, -0.0306)

    call g_tracer_add_param('htotal_scale_lo', bling%htotal_scale_lo, 0.01)
    call g_tracer_add_param('htotal_scale_hi', bling%htotal_scale_hi, 100.0)

    ! Fast gas exchange multiplier for DIC_sat
    ! Should be >30, but can be reduced during spinup to prevent instability
    call g_tracer_add_param('fast_gasex', bling%fast_gasex, 30.)   ! mol Fed mol PO4-1
    !-----------------------------------------------------------------------
    ! Uptake
    !-----------------------------------------------------------------------
    !
    ! Phytoplankton growth altered from Geider et al (1997)
    ! and Moore et al (2002). 
    ! The factor of 6.022e17 is to convert
    ! from umol to quanta and 2.77e18 to convert from quanta/sec
    ! to Watts given the average energy spectrum for underwater
    ! PAR from the Seabird sensor.  
    ! 
    call g_tracer_add_param('alpha_max', bling%alpha_max, 1.6e-5 *2.77e18/6.022e17) ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('alpha_min', bling%alpha_min, 0.4e-5 *2.77e18/6.022e17) ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('kappa_eppley', bling%kappa_eppley, 0.063)              ! deg C-1
    call g_tracer_add_param('resp_frac', bling%resp_frac, 0.0)                      ! dimensionless
    call g_tracer_add_param('pc_0',     bling%pc_0, 1.0e-5)                         ! s-1
    call g_tracer_add_param('thetamax_hi', bling%thetamax_hi, 0.040)                ! g Chl g C-1
    call g_tracer_add_param('thetamax_lo', bling%thetamax_lo, 0.010)                ! g Chl g C-1
    !
    ! Chl:C response rate constant for phytoplankton calibrated to 1 d-1
    ! after Owens et al (1980, Diel Periodicity in cellular Chlorophyll
    ! content of marine diatoms, Mar. Biol, 59, 71-77).
    !
    call g_tracer_add_param('gamma_irr_mem', bling%gamma_irr_mem, 1.0 / sperd)   ! s-1

    ! Introduce a minimum chlorophyll concentration for numerical stability.
    ! Value is an order of magnitude less than the minimum produced in topaz.
    !
    call g_tracer_add_param('chl_min', bling%chl_min, 1.e-5)                     ! ug kg-1
    !
    ! The biomass reponds to changes in growth rate with an arbitrary 2 day lag.
    !
    call g_tracer_add_param('gamma_biomass', bling%gamma_biomass, 0.5 / sperd)   ! s-1

    !-----------------------------------------------------------------------
    ! Monod half saturation coefficient for phosphate. Value of Aumont (JGR, 2002) 
    ! used for large phytoplankton.

    call g_tracer_add_param('k_po4', bling%k_po4,  1.0e-7)                       ! mol PO4 kg-1

    !-----------------------------------------------------------------------
    ! Fe uptake and limitation.
    ! The uptake ratio of Fe:P is determined from a Monod constant and a
    ! scaling factor. 
    ! The k_Fe_uptake is high, to provide luxury uptake of iron as a 
    ! relatively linear function of iron concentrations under open-ocean
    ! conditions, consistent with the results of Sunda and Huntsman (Fig 1, 
    ! Nature, 1997).

    call g_tracer_add_param('k_fe_uptake', bling%k_fe_uptake,  0.8e-9)           ! mol Fe kg-1

    ! This Monod term, which is nearly linear with [Fe], is multiplied by a
    ! scaling term to provide the actual Fe:P uptake ratio such that, at
    ! [Fe] = k_fe_uptake, Fe:P = fe_2_p_max / 2.
    ! This maximum value was set in accordance with the range of 
    ! open-ocean Fe:C ratios summarized by Boyd et al. (Science, 2007) and
    ! converted to a Fe:P ratio.
    ! As a tuning parameter, it affects the amount of Fe that cycles via the 
    ! organic matter pathway, and its ratio to k_fe_2_p determines the
    ! degree of iron limitation (the larger this ratio, the less iron
    ! limitation there will be).

    call g_tracer_add_param('fe_2_p_max', bling%fe_2_p_max, 28.e-6 * 106.)   ! mol Fed mol PO4-1

    ! The k_fe_2_p is the Fe:P at which the iron-limitation term has a 
    ! value of 0.5, chosen according to Sunda and Huntsman (Fig. 2, 
    ! Nature, 1997). Converted from Fe:C ratio.

    call g_tracer_add_param('k_fe_2_p', bling%k_fe_2_p,  7.e-6 * 106.)           ! mol Fe mol P-1

    ! In order to represent enzymatic plasticity and the ability of plankton
    ! to make do with very low Fe supply - including by liberating recalcitrant
    ! iron - the def_fe should not approach zero, but instead some small 
    ! positive number << 1, def_fe_min. This prevents unrealistically-strong
    ! limitation of phytoplankton under very low iron concentrations.

    call g_tracer_add_param('def_fe_min', bling%def_fe_min, 0.)           ! mol Fe mol P-1

    !-----------------------------------------------------------------------
    ! Mortality & Remineralization
    !-----------------------------------------------------------------------
    !
    ! T=0 phytoplankton specific total-mortality rate from the global
    ! synthesis of Dunne et al. (2005)
    !
    call g_tracer_add_param('lambda0', bling%lambda0, 0.19 / sperd)              ! s-1
    !
    ! Pivot phytoplankton concentration for grazing-based
    ! variation in ecosystem structure from the global
    ! synthesis of Dunne et al. (2005). Converted from mol C m-3.
    !
    call g_tracer_add_param('P_star', bling%P_star, 1.9e-3 / 1028. / 106.0)      ! mol P kg-1
    !
    ! Temperature-dependence of fractional detritus production
    ! from the global synthesis of Dunne et al. (2005)
    !
    call g_tracer_add_param('kappa_remin', bling%kappa_remin, -0.032)            ! deg C-1

    ! Phytoplankton fractional detritus production by size class,
    ! from the global synthesis of Dunne et al. (2005)
    call g_tracer_add_param('phi_lg', bling%phi_lg, 1.0)             ! unitless
    call g_tracer_add_param('phi_sm', bling%phi_sm, 0.18)            ! unitless
    !
    !-----------------------------------------------------------------------
    ! Dissolved Organic Material remineralization rate constants
    ! and fractional production ratios, all to be consistent
    ! with the work of Abell et al. (2000, Distributions of TOP, TON and TOC
    ! in the North pacific subtropical gyre: Implications for nutrient supply
    ! in the surface ocean and remineralization in the upper thermocline,
    ! J. Mar. Res., 58, 203-222).  
    ! Phi_dop is the fraction of non-sinking OM converted to DOP, and 
    ! gamma_dop is the first-order decay rate constant for DOP.
    ! 
    call g_tracer_add_param('phi_dop'  ,  bling%phi_dop, 0.1)                    ! dimensionless
    call g_tracer_add_param('gamma_dop',  bling%gamma_dop, 1.0 / (4.0 * spery))  ! s-1
    !
    !
    !-----------------------------------------------------------------------
    ! Remineralization
    !-----------------------------------------------------------------------
    !
    ! Stoichiometric ratios taken from Anderson (1995) as discussed in
    ! Sarmiento and Gruber (2008), and Sarmiento et al. (2002) for Ca:P.
    !
    call g_tracer_add_param('c_2_p', bling%c_2_p, 106.0 )                        ! mol C mol P-1
    call g_tracer_add_param('ca_2_p', bling%ca_2_p, 106.0 * .015 )                        ! mol C mol P-1
    call g_tracer_add_param('n_2_p', bling%n_2_p, 16. )                        ! mol C mol P-1
    call g_tracer_add_param('o2_2_p', bling%o2_2_p, 150.0 )                      ! mol O2 mol P-1
    ! Convert from mol P m-3 to mg C l-1
    call g_tracer_add_param('mass_2_p', bling%mass_2_p, 106. * 12.001 )          ! g C mol P-1

    ! Radiocarbon
    call g_tracer_add_param('half_life_14c', bling%half_life_14c, 5730.0 )       ! a

    !
    !-----------------------------------------------------------------------
    ! Remineralization length scales
    !
    ! Values of parameters to approximate upper e-folding of the globally-tuned
    ! "Martin curve" used in the OCMIP-II Biotic configuration of (z/75)^-0.9
    ! that gives a value of exp(-1) at 228 m from 75 m for an e-folding scale
    ! of 188 m.
    ! Here these are given as a linear function of depth, 
    !   wsink = wsink0 + wsink_acc * (z - wsink0_z)

    call g_tracer_add_param('wsink_acc', bling%wsink_acc, 0.05 / sperd)          ! s-1 
    call g_tracer_add_param('wsink0', bling%wsink0, 16.0 / sperd)                ! m s-1
    call g_tracer_add_param('wsink0_z', bling%wsink0_z, 80. )                    ! m
    call g_tracer_add_param('gamma_pop', bling%gamma_pop, 0.12 / sperd )         ! s-1

    ! Half saturation oxygen concentration for oxic remineralization rate.
    !
    call g_tracer_add_param('k_o2', bling%k_o2, 20.0e-6)                         ! mol O2 kg-1
    !
    ! Remineralization rate under suboxic/anoxic conditions, as a fraction of the rate under
    ! fully oxidized conditions.
    !
    call g_tracer_add_param('remin_min', bling%remin_min, 0.3)                   ! dimensionless
    !
    ! Minimum oxygen concentration for oxic remineralization.
    ! At O2 less than this, anaerobic remineralization occurs at remin_min rate.
    !
    call g_tracer_add_param('o2_min', bling%o2_min, 1.0e-06)                     ! mol O2 kg-1
    !
    ! Prevent oxygen from becoming negative. Setting to false allows negative
    ! oxygen in anoxic zones, which can be thought of as equivalent to 
    ! denitrification plus H2S production.
    !
    call g_tracer_add_param('prevent_neg_o2', bling%prevent_neg_o2, .true. ) 

    ! CaCO3 remineralization length scale. From Station P calibration (Dunne).
    !
    call g_tracer_add_param('ca_remin_depth', bling%ca_remin_depth, 1343. )      ! m

    ! Thickness of active sediment layer for  CaCO3 dissolution calculation
    !
    call g_tracer_add_param('z_sed',  bling%z_sed, 0.1 )                         ! m

    ! Global constant background (non-CaCO3) sedimentation flux, for CaCO3
    ! burial calculation. Given as a flux in g cm-2 ky-1, but then converted 
    ! to a 'mol lith m-2 a-1' flux, by dividing by 10, for consistency with 
    ! the calculations (from TOPAZ).
    !
    call g_tracer_add_param('sed_flux',  bling%sed_flux, 0.14 * .1 )             ! mol m-2 a-1 

    !-----------------------------------------------------------------------
    !       Iron Cycling
    !
    ! Global uniform iron ligand concentration.
    ! Taken from Parekh, P., M. J. Follows and E. A. Boyle (2005) Decoupling of iron
    ! and phosphate in the global ocean. Glob. Biogeochem. Cycles, 19, 
    ! doi: 10.1029/2004GB002280.
    !
    call g_tracer_add_param('felig_bkg', bling%felig_bkg, 1.0e-9)                ! mol ligand kg-1
    !
    ! Ratio of iron efflux from bottom sediment boundaries to the sedimenting phosphorus flux.
    ! From Elrod et al. (2004), 0.68 mmol Fe mol C-1, after Moore et al (2008):
    !
    call g_tracer_add_param('fe_2_p_sed', bling%fe_2_p_sed, 1.e-4 * 106.0 )      ! mol Fe mol P-1
    !
    ! 1.5-order iron scavenging in order to prevent high iron
    ! accumulations in high deposition regions (like the tropical
    ! Atlantic). This also helps prevent Fe accumulating in oligotrophic gyres and in 
    ! the abyssal ocean, where organic fluxes are low.  
    !
    call g_tracer_add_param('kfe_inorg', bling%kfe_inorg, 1.e3/sperd)            ! mol.5 Fe-.5 kg s-1
    !
    ! Equilibrium constant for (free and inorganically bound) iron binding with organic
    ! ligands taken from range similar to Parekh, P., M. J. Follows and E. A. Boyle 
    ! (2005) Decoupling of iron and phosphate in the global ocean. Glob. Biogeochem. 
    ! Cycles, 19, doi: 10.1029/2004GB002280.
    !
    call g_tracer_add_param('kfe_eq_lig_max', bling%kfe_eq_lig_max, 8.e10)       ! mol lig-1 kg
    !
    ! Minimum ligand strength under high light, to represent photodissociation of 
    ! ligand-Fe complexes.
    !
    call g_tracer_add_param('kfe_eq_lig_min', bling%kfe_eq_lig_min, 0.8e10)      ! mol lig-1 kg
    !
    ! Photodecay irradiance scaling.
    !
    call g_tracer_add_param('kfe_eq_lig_irr', bling%kfe_eq_lig_irr, 0.1)         ! W m-2
    !
    ! Iron concentration near which photodecay is compensated by enhanced siderophore
    ! production.
    !
    call g_tracer_add_param('kfe_eq_lig_femin', bling%kfe_eq_lig_femin, 0.05e-9) ! W m-2
    !
    ! Adsorption rate coefficient for detrital organic material.
    !
    call g_tracer_add_param('kfe_org', bling%kfe_org, 0.5/sperd)                 ! g org-1 m3 s-1
    !
    !-----------------------------------------------------------------------
    ! Miscellaneous
    !-----------------------------------------------------------------------
    !
    ! Debug flag to calculate global integrals for tracers
    !
    call g_tracer_add_param('tracer_debug',  bling%tracer_debug, .false.)
    !
    call g_tracer_end_param_list(package_name)
    !===========
    !Block Ends: g_tracer_add_param
    !===========

  end subroutine user_add_params

!#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Add all the tracers to be used in this module. 
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'

    !Add here only the parameters that are required at the time of registeration 
    !(to make flux exchanging Ocean tracers known for all PE's) 
    !
    call g_tracer_start_param_list(package_name)

    call g_tracer_add_param('ice_restart_file'   , bling%ice_restart_file  , 'ice_bling.res.nc')
    call g_tracer_add_param('ocean_restart_file' , bling%ocean_restart_file, 'ocean_bling.res.nc')
    call g_tracer_add_param('IC_file'       , bling%IC_file       , '')

    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file    = bling%ice_restart_file,&
         ocean_restart_file  = bling%ocean_restart_file )

    !All tracer fields shall be registered for diag output.

    !=====================================================
    !Specify all prognostic tracers of this modules.
    !=====================================================
    !User adds one call for each prognostic tracer below!
    !User should specify if fluxes must be extracted from boundary 
    !by passing one or more of the following methods as .true.  
    !and provide the corresponding parameters array
    !methods: flux_gas,flux_runoff,flux_wetdep,flux_drydep  
    !
    !Pass an init_value arg if the tracers should be initialized to a nonzero value everywhere
    !otherwise they will be initialized to zero.
    !
    !===========================================================
    !Prognostic Tracers
    !===========================================================
    !

    !       Dissolved Fe 
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'fed_b',                  &
         longname   = 'Dissolved Iron',         & 
         units      = 'mol/kg',                 &
         prog       = .true.,                   &
         flux_runoff    = .false.,              &
         flux_wetdep    = .true.,               &
         flux_drydep    = .true.,               &
         flux_param     = (/ 55.847e-03 /),     &
         flux_bottom    = .true. )


    !       DOP (Dissolved organic phosphorus)
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'dop_b',                  &
         longname   = 'DOP',                    &
         units      = 'mol/kg',                 &
         prog       = .true.)

    !
    !       O2
    !
    !NOTE: flux_gas_type = 'air_sea_gas_flux' is needed since the calculated alpha and csurf 
    !      in this module include the Schmidt number.
    !      If you want to pass Schmidt number separately refer to the calculation in generic_TOPAZ.F90
    !      which uses flux_gas_type = 'air_sea_gas_flux_generic' (which is default).
    !
    call g_tracer_add(tracer_list,package_name,                     &
         name       = 'o2_b',                                       &
         longname   = 'Oxygen',                                     &
         units      = 'mol/kg',                                     &
         prog       = .true.,                                       &
         flux_gas       = .true.,                                   &
         flux_bottom    = .true.,                                   &
         flux_gas_name  = 'o2_b_flux',                              &
         flux_gas_type  = 'air_sea_gas_flux',                       &
         flux_gas_molwt = WTMO2,                                    &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),               &
         flux_gas_restart_file  = 'ocean_bling_airsea_flux.res.nc' )

    !
    !       PO4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'po4_b',                   &
         longname   = 'Phosphate',             &
         units      = 'mol/kg',                &
         prog       = .true.,                  &
         flux_bottom    = .true.     )

    if (do_po4_pre) then                                          !<<PO4_PRE
    !
    !       PO4_pre
    !
    call g_tracer_add(tracer_list,package_name,  &
         name       = 'po4_pre',                 &
         longname   = 'Preformed phosphate',     &
         units      = 'mol/kg',                  &
         prog       = .true.)
    endif                                                         !PO4_PRE>>
    
    !===========================================================
    !Diagnostic Tracers
    !===========================================================
    !
    if (use_bling_chl) then    
    !       Chl (Chlorophyll)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'chl',                   &
         longname   = 'Chlorophyll',           &
         units      = 'ug kg-1',               &
         prog       = .false.,                 &
         init_value = 0.08          )
    else
    call g_tracer_add(tracer_list,package_name,&
         name       = 'chl_b',                 &
         longname   = 'Chlorophyll',           &
         units      = 'ug kg-1',               &
         prog       = .false.,                 &
         init_value = 0.08          )
    endif    
    !
    !       Biomass
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'biomass_p',             &
         longname   = 'Biomass in P units',    &
         units      = 'mol P kg-1',            &
         prog       = .false.)

    !       Irr_mem (Irradiance Memory)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'irr_mem_b',             &
         longname   = 'Irradiance memory',     &
         units      = 'Watts/m^2',             &
         prog       = .false.)


    if (do_carbon) then                                      !<<CARBON CYCLE
      if (bury_caco3) then                                     !<<BURY CACO3  
      !
      !       ALK (Total carbonate alkalinity)
      ! If CaCO3 is being buried, alkalinity is restored
      ! by riverine input (flux_runoff = true).
      !
      call g_tracer_add(tracer_list,package_name,&
           name       = 'alk_b',                   &
           longname   = 'Alkalinity',            &
           units      = 'mol/kg',                &
           prog       = .true.,                  &
           flux_runoff    = .true.,              &
           flux_param     = (/ 1.0e-03 /),       &
           flux_bottom    = .true.           )
      !
      !    Cased (CaCO3 concentration in active sediment layer)   
      !
      call g_tracer_add(tracer_list,package_name,&
           name       = 'cased_b',          &
           longname   = 'Sediment CaCO3 concentration', &
           units      = 'mol m-3',       &
           prog       = .false.          )
      else
      !
      !       ALK (Total carbonate alkalinity)
      !
      call g_tracer_add(tracer_list,package_name,&
           name       = 'alk_b',                   &
           longname   = 'Alkalinity',            &
           units      = 'mol/kg',                &
           prog       = .true.,                  &
           flux_runoff    = .false.,             &
           flux_param     = (/ 1.0e-03 /),       &
           flux_bottom    = .true.           )
      endif                                                    !BURY CACO3>>
    !
    !       DIC (Dissolved inorganic carbon)
    !
    call g_tracer_add(tracer_list,package_name,      &
         name       = 'dic_b',                       &
         longname   = 'Dissolved Inorganic Carbon',  &
         units      = 'mol/kg',                      &
         prog       = .true.,                        &
         flux_gas       = .true.,                    &
         flux_gas_name  = 'co2_b_flux',              &
         flux_gas_type  = 'air_sea_gas_flux',        &
         flux_gas_molwt = WTMCO2,                    &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),&
         flux_gas_restart_file  = 'ocean_bling_airsea_flux.res.nc',  &
         flux_runoff    = .false.,                   &
         flux_param     = (/12.011e-03  /),          &
         flux_bottom    = .true.,                    &
         init_value     = 0.001)
    !     
    !Diagnostic Tracers:
    !
    !       CO3_ion (Carbonate ion)
    !
    call g_tracer_add(tracer_list,package_name,    &
         name       = 'co3_ion_b',                   &
         longname   = 'Carbonate ion',             &
         units      = 'mol/kg',                    &
         prog       = .false. )
    !
    !       htotal (H+ ion concentration)
    !
    call g_tracer_add(tracer_list,package_name,  &
         name       = 'htotal_b',                  &
         longname   = 'H+ ion concentration',    &
         units      = 'mol/kg',                  &
         prog       = .false.,                   &
         init_value = bling%htotal_in)

      if (do_carbon_pre) then                                     !<<DIC_PRE  
      !       ALK_pre
      !
      call g_tracer_add(tracer_list,package_name,  &
           name       = 'alk_pre',                 &
           longname   = 'Preformed ALK',           &
           units      = 'mol/kg',                  &
           prog       = .true.)
      !          
      !       DIC_pre
      !
      call g_tracer_add(tracer_list,package_name,  &
           name       = 'dic_pre',                 &
           longname   = 'Preformed DIC',           &
           units      = 'mol/kg',                  &
           prog       = .true.)
      !   
      !       DIC_sat (Saturation Dissolved inorganic carbon)
      !
      call g_tracer_add(tracer_list,package_name,      &
           name       = 'dic_sat',                     &
           longname   = 'Saturation Dissolved Inorganic Carbon', &
           units      = 'mol/kg',                      &
           prog       = .true.,                        &
           flux_gas       = .true.,                    &
           flux_gas_name  = 'co2_sat_flux',            &
           flux_gas_type  = 'air_sea_gas_flux',        &
           flux_gas_molwt = WTMCO2,                    &
           flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),&
           flux_gas_restart_file  = 'ocean_bling_airsea_flux.res.nc',  &
           flux_param     = (/12.011e-03  /),          &
           init_value     = 0.001)
      !
      !Diagnostic tracers
      !
      !       htotal_sat (H+ ion concentration for DIC_sat)
      !
      call g_tracer_add(tracer_list,package_name,       &
           name       = 'htotal_sat',                   &
           longname   = 'H+ ion concentration for DIC_sat',    &
           units      = 'mol/kg',                       &
           prog       = .false.,                        &
           init_value = bling%htotal_in)
      endif                                                       !DIC_PRE>>
 
 
      if (do_14c) then                                        !<<RADIOCARBON
      !       D14IC (Dissolved inorganic radiocarbon)
      !
      call g_tracer_add(tracer_list,package_name,       &
         name       = 'di14c',                          &
         longname   = 'Dissolved Inorganic Radiocarbon',&
         units      = 'mol/kg',                         &
         prog       = .true.,                           &
         flux_gas       = .true.,                       &
         flux_gas_name  = 'c14o2_flux',                 &
         flux_gas_type  = 'air_sea_gas_flux',           &
         flux_gas_molwt = WTMCO2,                       &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),   &
         flux_gas_restart_file  = 'ocean_bling_airsea_flux.res.nc',  &
         flux_param     = (/14.e-03  /),                &
         flux_bottom    = .true.,                       &
         init_value     = 0.001)
      !
      !       DO14C (Dissolved organic radiocarbon)
      !
      call g_tracer_add(tracer_list,package_name, &
         name       = 'do14c',                    &
         longname   = 'DO14C',                    &
         units      = 'mol/kg',                   &
         prog       = .true.)
      endif                                                   !RADIOCARBON>>

    endif                                                    !CARBON CYCLE>>

  end subroutine user_add_tracers

!#######################################################################
! <SUBROUTINE NAME="generic_BLING_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracer fields could be modified after values are obtained from the 
  !  coupler. This subroutine is the place for specific tracer manipulations.
  !  BLING currently does not use this.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_BLING_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_BLING_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_BLING_update_from_coupler'

  end subroutine generic_BLING_update_from_coupler
!#######################################################################
  ! <SUBROUTINE NAME="generic_BLING_update_from_bottom">
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracers could have bottom fluxes and reservoirs. 
  !   This subroutine is the place for specific tracer manipulations.
  !   BLING currently does not use this.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_BLING_update_from_bottom(tracer_list,dt, tau) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment 
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_BLING_update_from_bottom(tracer_list, dt, tau)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    real, dimension(:,:,:),pointer :: grid_tmask

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

  end subroutine generic_BLING_update_from_bottom

!#######################################################################
  ! <SUBROUTINE NAME="generic_BLING_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This is the subroutine to contain most of the biogeochemistry for calculating the 
  !   interaction of tracers with each other and with outside forcings.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_BLING_update_from_source(tracer_list,Temp,Salt,dzt,hblt_depth,&
  !                                         ilb,jlb,tau,dt, grid_dat,sw_pen,opacity) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="Temp" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean temperature   
  !  </IN>
  !  <IN NAME="Salt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean salinity
  !  </IN>
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean layer thickness (meters)
  !  </IN>
  !  <IN NAME="opacity" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean opacity
  !  </IN>
  !  <IN NAME="sw_pen" TYPE="real, dimension(ilb:,jlb:)">
  !   Shortwave peneteration
  !  </IN>
  !  <IN NAME="hblt_depth" TYPE="real, dimension(ilb:,jlb:)">
  !   
  !  </IN>
  !  <IN NAME="grid_dat" TYPE="real, dimension(ilb:,jlb:)">
  !   Grid area
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_BLING_update_from_source(tracer_list,Temp,Salt,&
       rho_dzt,dzt,hblt_depth,ilb,jlb,tau,dt,grid_dat,model_time,nbands,   &
       max_wavelength_band,sw_pen_band,opacity_band)

    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp,Salt,rho_dzt,dzt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    integer,                        intent(in) :: ilb,jlb,tau
    real,                           intent(in) :: dt
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat
    type(time_type),                intent(in) :: model_time

    integer,                        intent(in) :: nbands
    real, dimension(:),             intent(in) :: max_wavelength_band
    real, dimension(:,ilb:,jlb:),   intent(in) :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band


    character(len=fm_string_len), parameter :: sub_name = 'generic_BLING_update_from_source'
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau, i, j, k , kblt,n
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast,grid_kmt

    integer :: nb
    logical :: used
    real :: tmp_hblt, tmp_Irrad, tmp_irrad_ML, tmp_opacity
    real, dimension(:), Allocatable :: tmp_irr_band 
    real :: s_over_p, fe_2_p
    real :: TK, PRESS, PKSPC

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

  ! SURFACE GAS FLUXES
  ! 
  ! This subroutine coordinates the calculation of gas concentrations and solubilities 
  ! in the surface layer. The concentration of a gas is written as csurf, while the
  ! solubility (in mol kg-1 atm-1 or mol m-3 atm-1) is written as alpha. These two
  ! quantities are passed to the coupler, which multiplies their difference by the
  ! gas exchange piston velocity over the mixed layer depth to provide the gas
  ! exchange flux,
  !    Flux = Kw/dz * (alpha - csurf)
  ! In order to simplify code flow, the Schmidt number parameters, which are part of 
  ! the piston velocity, are calculated here and applied to each of csurf and alpha 
  ! before being sent to the coupler.
  !
  ! For CO2 and 14CO2, the carbon solubility and speciation are calculated by the
  ! subroutine co2calc, following the OCMIP2 protocol. These calculations are both made
  ! using total CO2, following which the surface CO2 concentration (CO2*, also known as
  ! H2CO3*) is scaled by the DI14C/DIC ratio to give the surface 14CO2 concentration.
  ! The speciation calculation uses in situ temperature, salinity, ALK, PO4 and SiO4.
  !
  ! Oxygen solubility is calculated here, using in situ temperature and salinity.  

    !---------------------------------------------------------------------
    ! Get positive tracer concentrations for carbon calculation
    !---------------------------------------------------------------------

    call g_tracer_get_values(tracer_list,'po4_b' ,'field', bling%f_po4,isd,jsd,ntau=tau,positive=.true.)

    if (do_carbon) then                                      !<<CARBON CYCLE

    call g_tracer_get_values(tracer_list,'htotal_b','field', bling%f_htotal,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'alk_b'   ,'field', bling%f_alk   ,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dic_b'   ,'field', bling%f_dic   ,isd,jsd,ntau=tau)

    !---------------------------------------------------------------------
    !Calculate co3_ion.
    !Also calculate co2 fluxes csurf and alpha for the next round of exchange
    ! Note a scaled value of the PO4, rather than SiOH3, is used for all 
    ! calculations since there is no prognostic silica cycle 
    !---------------------------------------------------------------------
   
    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       bling%htotallo(i,j) = bling%htotal_scale_lo * bling%f_htotal(i,j,k)
       bling%htotalhi(i,j) = bling%htotal_scale_hi * bling%f_htotal(i,j,k)
    enddo; enddo ; !} i, j
 
    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         bling%f_dic(:,:,k),                          &
         bling%f_po4(:,:,k),                          &  
         bling%f_po4(:,:,k),                          &
         bling%f_alk(:,:,k),                          &
         bling%htotallo, bling%htotalhi,&
                                !InOut
         bling%f_htotal(:,:,k),                       & 
                                !OUT
         co2star=bling%co2_csurf(:,:), alpha=bling%co2_alpha(:,:), &
         pCO2surf=bling%pco2_surf(:,:), &
         co3_ion=bling%f_co3_ion(:,:,k))

    do k = 2, nk
       do j = jsc, jec ; do i = isc, iec  !{
          bling%htotallo(i,j) = bling%htotal_scale_lo * bling%f_htotal(i,j,k)
          bling%htotalhi(i,j) = bling%htotal_scale_hi * bling%f_htotal(i,j,k)
       enddo; enddo ; !} i, j
  
       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
            Temp(:,:,k), Salt(:,:,k),                    &
            bling%f_dic(:,:,k),                          &
            bling%f_po4(:,:,k),                          &  
            bling%f_po4(:,:,k),                          &
            bling%f_alk(:,:,k),                          &
            bling%htotallo, bling%htotalhi,&
                                !InOut
            bling%f_htotal(:,:,k),                       & 
                                !OUT
            co3_ion=bling%f_co3_ion(:,:,k))
    enddo

    call g_tracer_set_values(tracer_list,'htotal_b','field',bling%f_htotal  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion_b','field',bling%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'dic_b','alpha',bling%co2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic_b','csurf',bling%co2_csurf    ,isd,jsd)

      if (do_carbon_pre) then                                     !<<DIC_PRE  

    call g_tracer_get_values(tracer_list,'htotal_sat' ,'field',bling%f_htotal_sat,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'dic_sat'    ,'field',bling%f_dic_sat ,isd,jsd,ntau=tau)
    ! not needed for gas calc, but get them now for later
    call g_tracer_get_values(tracer_list,'alk_pre'    ,'field',bling%f_alk_pre ,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dic_pre'    ,'field',bling%f_dic_pre ,isd,jsd,ntau=tau)

    !---------------------------------------------------------------------
    !Calculate co2 fluxes csurf and alpha for the next round of exchange
    !---------------------------------------------------------------------
   
    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       bling%htotal_satlo(i,j) = bling%htotal_scale_lo * bling%f_htotal_sat(i,j,k)
       bling%htotal_sathi(i,j) = bling%htotal_scale_hi * bling%f_htotal_sat(i,j,k)
    enddo; enddo ; !} i, j
 
    ! The saturation carbon uses the same solubility as the regular DIC, so the
    ! alpha for dic_sat is simply set equal to the alpha for DIC. However, it
    ! has its own concentrations of things (and will use a faster gas exchange).
    ! The Htotal_sat is calculated for the whole ocean, even though it's only
    ! required in the surface layer (since the CO3= concentration related to 
    ! DIC_sat is not used).
 
    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         bling%f_dic_sat(:,:,k),                      &
         bling%f_po4(:,:,k),                          &  
         bling%f_po4(:,:,k),                          &  
         bling%f_alk(:,:,k),                          &
         bling%htotal_satlo, bling%htotal_sathi,&
                                !InOut
         bling%f_htotal_sat(:,:,k),                   & 
                                !OUT
         co2star=bling%co2_sat_csurf(:,:),            &
         pCO2surf=bling%pco2_sat_surf(:,:))

    do k = 2, nk
       do j = jsc, jec ; do i = isc, iec  !{
          bling%htotal_satlo(i,j) = bling%htotal_scale_lo * bling%f_htotal_sat(i,j,k)
          bling%htotal_sathi(i,j) = bling%htotal_scale_hi * bling%f_htotal_sat(i,j,k)
       enddo; enddo ; !} i, j
  
       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
            Temp(:,:,k), Salt(:,:,k),                    &
            bling%f_dic_sat(:,:,k),                          &
            bling%f_po4(:,:,k),                          &  
            bling%f_po4(:,:,k),                          &  
            bling%f_alk(:,:,k),                          &
            bling%htotal_satlo, bling%htotal_sathi,&
                                !InOut
            bling%f_htotal_sat(:,:,k))
    enddo

    call g_tracer_set_values(tracer_list,'htotal_sat','field',bling%f_htotal_sat  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'dic_sat','alpha',bling%co2_alpha        ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic_sat','csurf',bling%co2_sat_csurf    ,isd,jsd)

      endif                                                       !DIC_PRE>>
    
      if (do_14c) then                                        !<<RADIOCARBON
      
      ! Normally, the alpha would be multiplied by the atmospheric 14C/12C ratio. However,
      ! here that is set to 1, so that alpha_14C = alpha_12C. This needs to be changed!

    call g_tracer_get_values(tracer_list,'di14c' ,'field', bling%f_di14c,isd,jsd,ntau=tau)
    ! This is not used until later, but get it now
    call g_tracer_get_values(tracer_list,'do14c' ,'field', bling%f_do14c,isd,jsd,ntau=tau,positive=.true.)
    
       do j = jsc, jec ; do i = isc, iec  !{
       bling%c14o2_csurf(i,j) =  bling%co2_csurf(i,j) *                &
         bling%f_di14c(i,j,1) / (bling%f_dic(i,j,1) + epsln)
       bling%c14o2_alpha(i,j) =  bling%co2_alpha(i,j)
       enddo; enddo ; !} i, j

    call g_tracer_set_values(tracer_list,'di14c','alpha',bling%c14o2_alpha      ,isd,jsd)
    call g_tracer_set_values(tracer_list,'di14c','csurf',bling%c14o2_csurf      ,isd,jsd)

      endif                                                   !RADIOCARBON>>
      
    endif                                                    !CARBON CYCLE>>
    
    !---------------------------------------------------------------------
    ! Get positive concentrations for core tracers
    !---------------------------------------------------------------------
    call g_tracer_get_values(tracer_list,'fed_b'    ,'field',bling%f_fed       ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'dop_b'    ,'field',bling%f_dop       ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'o2_b'     ,'field',bling%f_o2        ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'biomass_p','field',bling%f_biomass_p ,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'irr_mem_b','field',bling%f_irr_mem   ,isd,jsd,ntau=1)

    if (do_po4_pre) &
    call g_tracer_get_values(tracer_list,'po4_pre'  ,'field',bling%f_po4_pre   ,isd,jsd,ntau=tau,positive=.true.)
    
    if (do_carbon) then
    call g_tracer_get_values(tracer_list,'co3_ion_b','field',bling%f_co3_ion   ,isd,jsd,ntau=1,positive=.true.)
      if (bury_caco3) &
      call g_tracer_get_values(tracer_list,'cased_b'    ,'field',bling%f_cased     ,isd,jsd,ntau=1)
    endif

    bling%zt = 0.0
    s_over_p = 0.0

 !--------------------------------------------------------------------------
 ! NUTRIENT UPTAKE
 !--------------------------------------------------------------------------
   
    ! Available light calculation
    !-----------------------------------------------------------------------
    ! There are multiple types of light.
    !   irr_inst is the instantaneous irradiance field.
    !   irr_mix is the same, but with the irr_inst averaged throughout the  
    ! mixed layer as defined in the KPP routine plus one more vertical box 
    ! to account for mixing directly below the boundary layer. This quantity  
    ! is intended to represent the light to which phytoplankton subject to
    ! turbulent transport in the mixed-layer would be exposed.
    !   irr_mem is a temporally smoothed field carried between timesteps, to 
    ! represent photoadaptation.
    !-----------------------------------------------------------------------

    allocate(tmp_irr_band(nbands))
    do j = jsc, jec ; do i = isc, iec   !{

       do nb=1,nbands !{
          if (max_wavelength_band(nb) .lt. 710) then !{
             tmp_irr_band(nb) = sw_pen_band(nb,i,j)
          else
             tmp_irr_band(nb) = 0.0
          endif !}
       enddo !}

       kblt = 0 ; tmp_irrad_ML = 0.0 ; tmp_hblt = 0.0
       do k = 1, nk !{
          tmp_Irrad = 0.0
          do nb=1,nbands !{
             tmp_opacity = opacity_band(nb,i,j,k)
             tmp_Irrad = tmp_Irrad + tmp_irr_band(nb) * exp(-tmp_opacity * &
               dzt(i,j,k) * 0.5)
             ! Change tmp_irr_band from being the value atop layer k to the 
             ! value at the bottom of layer k.
             tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity *      &
               dzt(i,j,k))
          enddo !}
          bling%irr_inst(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
          bling%irr_mix(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
          if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then !{
             kblt = kblt+1
             tmp_irrad_ML = tmp_irrad_ML + bling%irr_mix(i,j,k) * dzt(i,j,k)
             tmp_hblt = tmp_hblt + dzt(i,j,k)
          endif !}
       enddo !} k-loop
       bling%irr_mix(i,j,1:kblt) = tmp_irrad_ML / max(1.0e-6,tmp_hblt)

    enddo;  enddo !} i,j

    deallocate(tmp_irr_band)

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        

       !--------------------------------------------------------------------
       ! Phytoplankton photoadaptation timescale
       
       bling%f_irr_mem(i,j,k) = (bling%f_irr_mem(i,j,k) +                   &
            (bling%irr_mix(i,j,k) - bling%f_irr_mem(i,j,k)) * min( 1.0 ,   &
            bling%gamma_irr_mem * dt)) * grid_tmask(i,j,k)

       !--------------------------------------------------------------------
       ! Temperature functionality of growth and grazing
       ! NB The temperature effect of Eppley (1972) is used instead
       !    of that in Geider et al (1997) for both simplicity and
       !    to incorporate combined effects on uptake, incorporation
       !    into organic matter and photorespiration.  Values of PCmax
       !    are normalized to 0C rather than 20C in Geider et al. (1997)
        
       bling%expkT(i,j,k) = exp(bling%kappa_eppley * Temp(i,j,k))
     
    enddo; enddo ; enddo !} i,j,k

    !-----------------------------------------------------------------------
    ! Phytoplankton are assumed to grow according to the general properties 
    ! described in Geider (1997). This formulation gives a biomass-specific 
    ! growthrate as a function of light, nutrient limitation, and 
    ! temperature. We modify this relationship slightly here, as described 
    ! below, and also use the assumption of steady state growth vs. loss to 
    ! derive a simple relationship between growth rate, biomass and uptake.
    !
    !-----------------------------------------------------------------------
    ! First, we calculate the limitation terms for PO4 and Fe, and the 
    ! Fe-limited Chl:C maximum.
    ! The light-saturated maximal photosynthesis rate term (pc_m) is simply 
    ! the product of a prescribed maximal photosynthesis rate (pc_0), the 
    ! Eppley temperature dependence, and a Liebig limitation (the minimum
    ! of Michaelis-Menton PO4-limitation, or iron-limitation).
    ! The iron limitation term has a lower limit of def_fe_min 
    ! and is scaled by (k_fe_2_p + fe_2_p_max) / fe_2_p_max
    ! so that it approaches 1 as fed approaches infinity. Thus, 
    ! it's of comparable magnitude to the PO4 limitation term.
    !
    ! Fe limitation acts in two additional mechanisms:
    ! 1. By reducing the maximum achievable Chl:C ratio 
    ! (theta) below a prescribed, Fe-replete maximum value (thetamax), to 
    ! approach a prescribed minimum Chl:C (thetamin) under extreme
    ! Fe-limitation.
    ! 2. By reducing alpha (the initial slope of the P-I curve) under Fe-
    ! limitation.
    !-----------------------------------------------------------------------
    
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{      
         bling%fe_2_p_uptake(i,j,k) = bling%fe_2_p_max *                   &
           bling%f_fed(i,j,k) / (bling%k_fe_uptake + bling%f_fed(i,j,k))
         bling%def_fe(i,j,k) = bling%def_fe_min +                          &
           (1. - bling%def_fe_min) * bling%fe_2_p_uptake(i,j,k) /          & 
           (bling%k_fe_2_p + bling%fe_2_p_uptake(i,j,k)) *                 &
           (bling%k_fe_2_p + bling%fe_2_p_max) / bling%fe_2_p_max 
         bling%pc_m(i,j,k) = bling%pc_0 * bling%expkT(i,j,k) * min(      &
           bling%f_po4(i,j,k) / (bling%k_po4 + bling%f_po4(i,j,k)) ,       &
           bling%def_fe(i,j,k))
         bling%thetamax_fe(i,j,k) = bling%thetamax_lo +                    &
           (bling%thetamax_hi - bling%thetamax_lo) * bling%def_fe(i,j,k)
         bling%alpha(i,j,k) = bling%alpha_min +                            &
           (bling%alpha_max - bling%alpha_min) * bling%def_fe(i,j,k)
            
    !-----------------------------------------------------------------------
    ! Next, the nutrient-limited efficiency of algal photosystems, Irrk, is
    ! calculated. This requires a prescribed quantum yield, alpha.
    ! The iron deficiency term is included here as a multiplier of the 
    ! thetamax_fe to represent the importance of Fe in forming chlorophyll
    ! accessory antennae, which do not affect the Chl:C but still affect the
    ! phytoplankton ability to use light (eg Stzrepek & Harrison Nature 
    ! 2004).

          bling%irrk(i,j,k) = (bling%pc_m(i,j,k) / ( epsln +               &
            bling%alpha(i,j,k) * bling%thetamax_fe(i,j,k) )) +             &
            bling%f_irr_mem(i,j,k) * 0.5

    !-----------------------------------------------------------------------
    ! We also calculate the Chl:C ratio here, although it does not enter  
    ! into the uptake calculation and is only used for the diagnostic
    ! chlorophyll concentration, below.

          bling%theta(i,j,k) = bling%thetamax_fe(i,j,k) / (1. +            &
            bling%thetamax_fe(i,j,k) * bling%alpha(i,j,k) *                &
            bling%f_irr_mem(i,j,k) / (epsln + 2. * bling%pc_m(i,j,k)))
            
    !-----------------------------------------------------------------------
    ! Now we can calculate the carbon-specific photosynthesis rate, pc_tot.
    
          bling%pc_tot(i,j,k) = bling%pc_m(i,j,k) *                        &
            (1. - exp(-bling%irr_mix(i,j,k) / (epsln + bling%irrk(i,j,k))))

    !-----------------------------------------------------------------------
    ! Next, we account for the maintenance effort that phytoplankton must 
    ! exert in order to combat decay. This is prescribed as a fraction of the
    ! light-saturated photosynthesis rate, resp_frac. The result of this is 
    ! to set a level of energy availability below which net growth (and 
    ! therefore nutrient uptake) is zero, given by resp_frac * pc_m.

          bling%mu(i,j,k) = max (0.,                                       &              
            bling%pc_tot(i,j,k) - bling%resp_frac * bling%pc_m(i,j,k))

    !-----------------------------------------------------------------------
    ! We now must convert this net carbon-specific growth rate to nutrient 
    ! uptake rates, the quantities we are interested in. Since we have no 
    ! explicit biomass tracer, we use the result of Dunne et al. (GBC, 2005) 
    ! to calculate an implicit biomass from the uptake rate through the  
    ! application of a simple idealized grazing law. This has the effect of 
    ! reducing uptake in low growth-rate regimes and increasing uptake in 
    ! high growth-rate regimes - essentially a non-linear amplification of 
    ! the growth rate variability. The result is:
    
          bling%biomass_p_ts(i,j,k) =                                      &
            ((bling%mu(i,j,k)/(bling%lambda0 * bling%expkT(i,j,k)))**3     &
            + (bling%mu(i,j,k)/(bling%lambda0 * bling%expkT(i,j,k))))      &
            * bling%p_star

          bling%f_biomass_p(i,j,k) = bling%f_biomass_p(i,j,k) +            &
            (bling%biomass_p_ts(i,j,k) - bling%f_biomass_p(i,j,k)) *       &
            min(1.0, bling%gamma_biomass * dt) * grid_tmask(i,j,k)
          
          bling%jp_uptake(i,j,k) = bling%f_biomass_p(i,j,k) *              &
            bling%mu(i,j,k)

    ! We can now use the diagnostic biomass to calculate the chlorophyll
    ! concentration:
    
          bling%f_chl(i,j,k) = max(bling%chl_min, bling%f_biomass_p(i,j,k) &
            * bling%c_2_p * 12.011e6 * bling%theta(i,j,k)) *               &
            grid_tmask(i,j,k)

    !-----------------------------------------------------------------------
    ! Iron is then taken up as a function of PO4 uptake and iron limitation,
    ! with a maximum Fe:P uptake ratio of fe2p_max:
    
          bling%jfe_uptake(i,j,k) = bling%jp_uptake(i,j,k) *               &
            bling%fe_2_p_uptake(i,j,k)

    enddo; enddo ; enddo !} i,j,k
    
  !-------------------------------------------------------------------------
  ! PARTITIONING BETWEEN ORGANIC POOLS
  !-------------------------------------------------------------------------
   
  ! The uptake of nutrients is assumed to contribute to the growth of
  ! phytoplankton, which subsequently die and are consumed by heterotrophs.
  ! This can involve the transfer of nutrient elements between many
  ! organic pools, both particulate and dissolved, with complex histories.
  ! We take a simple approach here, partitioning the total uptake into two
  ! fractions - sinking and non-sinking - as a function of temperature, 
  ! following Dunne et al. (2005). 
  ! Then, the non-sinking fraction is further subdivided, such that the 
  ! majority is recycled instantaneously to the inorganic nutrient pool,
  ! representing the fast turnover of labile dissolved organic matter via
  ! the microbial loop, and the remainder is converted to semi-labile
  ! dissolved organic matter. Iron and phosphorus are treated identically 
  ! for the first step, but all iron is recycled instantaneously in the
  ! second step (i.e. there is no dissolved organic iron pool).
  !-------------------------------------------------------------------------
    
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        

          bling%frac_pop(i,j,k) = (bling%phi_sm + bling%phi_lg *           &
            (bling%mu(i,j,k)/(bling%lambda0*bling%expkT(i,j,k)))**2.)/     &
            (1. +                                                          &
            (bling%mu(i,j,k)/(bling%lambda0*bling%expkT(i,j,k)))**2.)*     &
            exp(bling%kappa_remin * Temp(i,j,k))

          bling%jpop(i,j,k) = bling%frac_pop(i,j,k) * bling%jp_uptake(i,j,k)

          bling%jpofe(i,j,k) = bling%frac_pop(i,j,k)*bling%jfe_uptake(i,j,k)
 
     ! As a helpful diagnostic, the implied fraction of production by large 
     ! phytoplankton is calculated, also following Dunne et al. 2005. This
     ! could be done more simply, but is done here in a complicated way as
     ! a sanity check. Looks fine.
     ! Note the calculation is made in P units, rather than C.
     ! This is also used for the CaCO3 production.

          s_over_p = ( -1. + ( 1. + 4. * bling%jp_uptake(i,j,k) / &
            (bling%expkT(i,j,k) * bling%lambda0 * bling%p_star))**0.5) * .5

          bling%frac_lg(i,j,k) = s_over_p / (1 + s_over_p)
          
    !-----------------------------------------------------------------------
    ! Then the remainder is divided between instantaneously recycled and
    ! long-lived dissolved organic matter,
            
          bling%jdop(i,j,k) = bling%phi_dop * (bling%jp_uptake(i,j,k) - &
            bling%jpop(i,j,k))
            
          bling%jp_recycle(i,j,k) = bling%jp_uptake(i,j,k) - &
            bling%jpop(i,j,k) - bling%jdop(i,j,k)

          bling%jfe_recycle(i,j,k) = bling%jfe_uptake(i,j,k) -  &
            bling%jpofe(i,j,k)

    enddo; enddo ; enddo !} i,j,k
    
    if (do_carbon) then                                     !<<CARBON CYCLE

    !-----------------------------------------------------------------------
    ! Calcium carbonate production
    ! Alkalinity is consumed through the production of CaCO3. Here, this is
    ! simply a linear function of the implied growth rate of small
    ! phytoplankton, which gave a reasonably good fit to the global 
    ! observational synthesis of Dunne (in prep., 2009). This is consistent
    ! with the findings of Jin et al. (GBC,2006).
    !-----------------------------------------------------------------------

         do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
            bling%jca_uptake(i,j,k) = (1. - bling%frac_lg(i,j,k)) *        &
              bling%jp_uptake(i,j,k) * bling%ca_2_p
         enddo; enddo ; enddo !} i,j,k
       
      if (do_14c) then                                        !<<RADIOCARBON
         do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

          bling%c14_2_p(i,j,k) = bling%c_2_p *                             &
            bling%f_di14c(i,j,k) / (epsln + bling%f_dic(i,j,k))

          bling%jdo14c(i,j,k) = bling%phi_dop * (bling%jp_uptake(i,j,k) -  &
            bling%jpop(i,j,k)) * bling%c14_2_p(i,j,k)
            
         enddo; enddo ; enddo !} i,j,k
      endif                                                   !RADIOCARBON>>
       
    endif                                                    !CARBON CYCLE>>

  !-------------------------------------------------------------------------
  ! SINKING AND REMINERALIZATION
  !-------------------------------------------------------------------------
    ! Calculate the depth of each grid cell (needs to be 3d for use with
    ! isopycnal co-ordinate model).

    do j = jsc, jec ;      do i = isc, iec   !{
       bling%zt(i,j,1) = dzt(i,j,1)
    enddo; enddo !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       bling%zt(i,j,k) = bling%zt(i,j,k-1) + dzt(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    !-----------------------------------------------------------------------
    ! Calculate the remineralization lengthscale matrix, zremin, a function 
    ! of z. Sinking rate (wsink) is constant over the upper wsink0_z metres,
    ! then  increases linearly with depth.
    ! The remineralization rate is a function of oxygen concentrations,
    ! following a Holling type 2 dependence, decreasing to a minimum value
    ! of remin_min. This is ad hoc, following work by Bianchi, Sarmiento,
    ! Galbraith and Kwon (unpublished).

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
      
        if (bling%zt(i,j,k) .lt. bling%wsink0_z) then !{
          bling%wsink(i,j,k) = bling%wsink0
        else
          bling%wsink(i,j,k) = (bling%wsink_acc * (bling%zt(i,j,k) -       &
            bling%wsink0_z) + bling%wsink0)
        endif  !}

        bling%zremin(i,j,k) = bling%gamma_pop * (bling%f_o2(i,j,k)**2 /    &
          (bling%k_o2**2 + bling%f_o2(i,j,k)**2) * (1. - bling%remin_min)+ &
          bling%remin_min) / (bling%wsink(i,j,k) + epsln)

    enddo; enddo ; enddo !} i,j,k

    if (do_carbon) then                                       !<<CARBON CYCLE
     do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
     
     ! Using Sayles for solubility (will change to Mucci later)
       bling%co3_solubility(i,j,k) = max(4.95e-7 * exp ( 0.05021 / &
            (Temp(i,j,k) + 273.15) * bling%zt(i,j,k)) * 3.42031e-3 *       &
            bling%Rho_0 * bling%Rho_0 / max(epsln, Salt(i,j,k)),epsln)
      
     ! CaCO3 dissolution lengthscale is a function of the saturation state,
     ! CO3 / CO3solubility, such that the lengthscale decreases from
     ! infinity at CO3 > CO3solubility to zero at CO3 = 0.
     
     bling%zremin_caco3(i,j,k) = 1. / bling%ca_remin_depth * (1.0 -        &
       min(1., bling%f_co3_ion(i,j,k) / (bling%co3_solubility(i,j,k) +     &
       epsln)))
         
     enddo; enddo ; enddo !} i,j,k
     
    ! Generate CaCO3 sinking flux, and dissolve it through the water column.
    ! Same as for other elements - see below for more detailed explanation.

    do j = jsc, jec ;      do i = isc, iec   !{

      bling%fcaco3(i,j,1) = bling%jca_uptake(i,j,1) * rho_dzt(i,j,1) /     &
        (1.0 + dzt(i,j,1) * bling%zremin_caco3(i,j,1)) 

      bling%jca_reminp(i,j,1) =                                            &
       (bling%jca_uptake(i,j,1) * rho_dzt(i,j,1) - bling%fcaco3(i,j,1)) /  &
       (epsln + rho_dzt(i,j,1))
       
    enddo; enddo !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{

      bling%fcaco3(i,j,k) = (bling%fcaco3(i,j,k-1) +                       &
        bling%jca_uptake(i,j,k) * rho_dzt(i,j,k)) /                        &
        (1.0 + dzt(i,j,k) * bling%zremin_caco3(i,j,k)) 

      bling%jca_reminp(i,j,k) = (bling%fcaco3(i,j,k-1) +                   &
       bling%jca_uptake(i,j,k) * rho_dzt(i,j,k) - bling%fcaco3(i,j,k)) /   &
       (epsln + rho_dzt(i,j,k))
       
    enddo; enddo ; enddo !} i,j,k

      if (do_14c) then                                        !<<RADIOCARBON

      ! Sinking particulate 14C is generated in the local ratio of 14C/12C
      ! to sinking 12C, which itself is strictly tied to P through a fixed
      ! C:P. Therefore, jpop can be used to calculate fpo14c.

      do j = jsc, jec ;      do i = isc, iec   !{
        bling%fpo14c(i,j,1) = bling%jpop(i,j,1) * bling%c14_2_p(i,j,1) *   &
          rho_dzt(i,j,1) /                                                 &
          (1.0 + dzt(i,j,1) * bling%zremin(i,j,1)) 

        bling%j14c_reminp(i,j,1) = (bling%jpop(i,j,1) *                    &
          bling%c14_2_p(i,j,1) * rho_dzt(i,j,1) - bling%fpo14c(i,j,1)) /   &
         (epsln + rho_dzt(i,j,1))
      enddo; enddo !} i,j

      do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
        bling%fpo14c(i,j,k) = (bling%fpo14c(i,j,k-1) +                     &
          bling%jpop(i,j,k) * bling%c14_2_p(i,j,k) * rho_dzt(i,j,k)) /     &
          (1.0 + dzt(i,j,k) * bling%zremin(i,j,k)) 

        bling%j14c_reminp(i,j,k) = (bling%fpo14c(i,j,k-1) +                &
         bling%jpop(i,j,k) * bling%c14_2_p(i,j,k) * rho_dzt(i,j,k) -       &
         bling%fpo14c(i,j,k)) / (epsln + rho_dzt(i,j,k))
      enddo; enddo ; enddo !} i,j,k

     ! Decay the radiocarbon in both DIC and DOC
     
     bling%lambda_14c = log(2.0) / (bling%half_life_14c * spery)

      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        

        bling%j14c_decay_dic(i,j,k) = bling%f_di14c(i,j,k) *               &
          bling%lambda_14c 

        bling%j14c_decay_doc(i,j,k) = bling%f_do14c(i,j,k) *               &
          bling%lambda_14c 

      enddo; enddo ; enddo !} i,j,k
      endif                                                   !RADIOCARBON>>

    endif                                                    !CARBON CYCLE>>

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        

     !---------------------------------------------------------------------
     ! Calculate free and inorganically associated iron concentration for
     ! scavenging.
     ! We assume that there is a 
     ! spectrum of iron ligands present in seawater, with varying binding
     ! strengths and whose composition varies with light and iron 
     ! concentrations. For example, photodissocation of ligand complexes 
     ! occurs under bright light, weakening the binding strength 
     ! (e.g. Barbeau et al., Nature 2001), while at very low iron 
     ! concentrations (order kfe_eq_lig_femin), siderophores are thought
     ! to be produced as a response to extreme
     ! iron stress.
     ! In anoxic waters, iron should be reduced, and therefore mostly 
     ! immune to scavenging. Easiest way to do this is to skip the feprime
     ! calculation if oxygen is less than 0.
     
     if (bling%f_o2(i,j,k) .gt. bling%o2_min) then  !{
     bling%kfe_eq_lig(i,j,k) = bling%kfe_eq_lig_max -                      &
       (bling%kfe_eq_lig_max - bling%kfe_eq_lig_min) *                     &
       (bling%irr_inst(i,j,k)**2. / (bling%irr_inst(i,j,k)**2. +       &
       bling%kfe_eq_lig_irr **2.)) * max(0., min(1., (bling%f_fed(i,j,k) - &
       bling%kfe_eq_lig_femin) / (epsln + bling%f_fed(i,j,k)) * 1.2))

     bling%feprime(i,j,k) = 1.0 + bling%kfe_eq_lig(i,j,k) *                &
       (bling%felig_bkg - bling%f_fed(i,j,k))
     
     bling%feprime(i,j,k) = (-bling%feprime(i,j,k) +(bling%feprime(i,j,k)* &
       bling%feprime(i,j,k) + 4.0 * bling%kfe_eq_lig(i,j,k) *              &
       bling%f_fed(i,j,k))**(0.5)) /(2.0 * bling%kfe_eq_lig(i,j,k))
     else
     bling%feprime(i,j,k) = 0.
     endif !}

     bling%jfe_ads_inorg(i,j,k) = min(0.5/dt, bling%kfe_inorg * &
       bling%feprime(i,j,k) ** 0.5) * bling%feprime(i,j,k) 
     
    enddo; enddo ; enddo  !} i,j,k

    !---------------------------------------------------------------------
    ! In general, the flux at the bottom of a grid cell should equal
    ! Fb = (Ft + Prod*dz) / (1 + zremin*dz)
    ! where Ft is the flux at the top, and prod*dz is the integrated 
    ! production of new sinking particles within the layer.
    ! Since Ft=0 in the first layer,

    do j = jsc, jec ;      do i = isc, iec   !{

      bling%fpop(i,j,1) = bling%jpop(i,j,1) * rho_dzt(i,j,1) /             &
        (1.0 + dzt(i,j,1) * bling%zremin(i,j,1)) 

    !-----------------------------------------------------------------------
    ! Now, calculate the Fe adsorption using this fpop:
    ! The absolute first order rate constant is calculated from the 
    ! concentration of organic particles, after Parekh et al. (2005). Never
    !  allowed to be greater than 1/2dt for numerical stability.

     bling%jfe_ads_org(i,j,1) = min (0.5/dt,                               &
       bling%kfe_org * (bling%fpop(i,j,1) / (epsln + bling%wsink(i,j,1)) * &
       bling%mass_2_p) ** 0.58) * bling%feprime(i,j,1)
       
      bling%fpofe(i,j,1) = (bling%jpofe(i,j,1) +bling%jfe_ads_inorg(i,j,1) &
        + bling%jfe_ads_org(i,j,1)) * rho_dzt(i,j,1) /                     &
        (1.0 + dzt(i,j,1) * bling%zremin(i,j,1)) 
      
    !-----------------------------------------------------------------------
    ! Calculate remineralization terms

      bling%jp_reminp(i,j,1) =                                             &
       (bling%jpop(i,j,1) * rho_dzt(i,j,1) - bling%fpop(i,j,1)) /          &
       (epsln + rho_dzt(i,j,1))
       
      bling%jfe_reminp(i,j,1) =                                            &
       ((bling%jpofe(i,j,1) + bling%jfe_ads_org(i,j,1) +                   &
       bling%jfe_ads_inorg(i,j,1)) * rho_dzt(i,j,1) -                      &
       bling%fpofe(i,j,1)) / (epsln + rho_dzt(i,j,1))

     enddo; enddo !} i,j

    !-----------------------------------------------------------------------
    ! Then, for the rest of water column, include Ft:

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{

      bling%fpop(i,j,k) = (bling%fpop(i,j,k-1) +                           &
        bling%jpop(i,j,k) * rho_dzt(i,j,k)) /                              &
        (1.0 + dzt(i,j,k) * bling%zremin(i,j,k)) 

    !-----------------------------------------------------------------------
    ! Again, calculate the Fe adsorption using this fpop:

     bling%jfe_ads_org(i,j,k) = min (0.5/dt,                               &
       bling%kfe_org * (bling%fpop(i,j,k) / (epsln + bling%wsink(i,j,k)) * &
       bling%mass_2_p) ** 0.58) * bling%feprime(i,j,k)
       
      bling%fpofe(i,j,k) = (bling%fpofe(i,j,k-1) +                         &
        (bling%jfe_ads_org(i,j,k) + bling%jfe_ads_inorg(i,j,k) +           &
        bling%jpofe(i,j,k)) *rho_dzt(i,j,k)) /                             &
        (1.0 + dzt(i,j,k) * bling%zremin(i,j,k)) 

    !---------------------------------------------------------------------
    ! Calculate remineralization terms

      bling%jp_reminp(i,j,k) = (bling%fpop(i,j,k-1) +                      &
       bling%jpop(i,j,k) * rho_dzt(i,j,k) - bling%fpop(i,j,k)) /           &
       (epsln + rho_dzt(i,j,k))
       
      bling%jfe_reminp(i,j,k) = (bling%fpofe(i,j,k-1) +                    &
       (bling%jfe_ads_org(i,j,k) + bling%jfe_ads_inorg(i,j,k) +            & 
         bling%jpofe(i,j,k)) * rho_dzt(i,j,k) -                            &
       bling%fpofe(i,j,k)) / (epsln + rho_dzt(i,j,k))

    enddo; enddo ; enddo !} i,j,k

    !---------------------------------------------------------------------
    ! BOTTOM LAYER 
    ! Account for remineralization in bottom box, and bottom fluxes

    do j = jsc, jec ; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{

      !---------------------------------------------------------------------
      ! Calculate iron addition from sediments as a function of organic
      ! matter supply.

        bling%ffe_sed(i,j) = bling%fe_2_p_sed * bling%fpop(i,j,k)

      ! Added the burial flux of sinking particulate iron here as a 
      ! diagnostic, needed to calculate mass balance of iron.

        bling%fe_burial(i,j) = bling%fpofe(i,j,k)

      !---------------------------------------------------------------------
      ! Calculate external bottom fluxes for tracer_vertdiff. Positive fluxes
      ! are from the water column into the seafloor. For P, the bottom flux  
      ! puts the sinking flux reaching the bottom cell into the water column 
      ! through diffusion. For iron, the sinking flux disappears into the 
      ! sediments if bottom waters are oxic (assumed adsorbed as oxides),
      ! while an efflux of dissolved iron occurs dependent on the supply of
      ! reducing organic matter (scaled by the org-P sedimentation rate).
      ! If bottom waters are anoxic, the sinking flux of Fe is returned to
      ! the water column. Note this is not appropriate for very long runs
      ! with an anoxic ocean (iron will keep accumulating forever).
      ! For oxygen, the consumption of oxidant required to respire  
      ! the settling flux of organic matter (in support of the
      ! PO4 bottom flux) diffuses from the bottom water into the sediment.

        bling%b_po4(i,j) = - bling%fpop(i,j,k)

        if (bling%f_o2(i,j,k) .gt. bling%o2_min) then  !{
          bling%b_fed(i,j) = - bling%ffe_sed(i,j) 
        else
          bling%b_fed(i,j) = - bling%ffe_sed(i,j) - bling%fpofe(i,j,k)
        endif !}      

        bling%b_o2(i,j) = bling%o2_2_p * bling%fpop(i,j,k)

      endif !}
    enddo; enddo  !} i, j

    call g_tracer_set_values(tracer_list,'fed_b', 'btf', bling%b_fed ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4_b', 'btf', bling%b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2_b',  'btf', bling%b_o2 ,isd,jsd)

    if (do_carbon) then                                !<<CARBON CYCLE
    ! Do bottom box calcs for carbon cycle
    
      do j = jsc, jec ; do i = isc, iec  !{
         k = grid_kmt(i,j)
         if (k .gt. 0) then !{
           ! CaCO3 dissolved in the sediment is returned as an alkalinity 
           ! bottom flux. This can allow for some burial of CaCO3, or return
           ! all CaCO3 to the water column to preserve the alkalinity 
           ! inventory at the user's preference.
           ! In addition, the flux of NO3 out of the sediment, inferred from 
           ! the PO4 flux, causes a negative flux of alkalinity.
           ! As a diagnostic,
           bling%fcaco3_to_sed(i,j) = bling%fcaco3(i,j,k)
           
           if (bury_caco3) then 
           !-----------------------------------------------------------------
           ! Determine the flux of CaCO3 retained in sediment using the Dunne
           ! et al (in prep) metamodel calibrated to the Hales (2003) steady 
           ! state 1-d sediment model. In short, this uses the bottom water 
           ! calcite saturation state (omega), as well as additional
           ! dissolution caused by the remineralization of organic matter, to
           ! calculate dissolution within a 10cm-thick active sediment layer. 
           ! The sediment is assumed to accumulate at a rate determined by
           ! the sum of 
           ! so that, assuming a density of 2.7 g cm-3,
           ! a porosity of 0.7, and an average molecular weight of 100 for 
           ! sedimentary components to give: 2.7e6*(1-0.7)=8.1e3 mol m-3.  
           !
           ! The permanent burial flux of CaCO3 is then given by the product
           ! of the concentration of CaCO3 in the active sediment layer and 
           ! the mass accumulation rate, and a new concentration of CaCO3 in 
           ! the active sediment layer is calcuated from the sum of the 
           ! fluxes.
           !
           ! Note that, at steady state, the burial efficiency, 
           !     f = fcased_burial/f_cadet_btf
           ! reduces to (in TOPAZ nomenclature):
           !
           !  f = min(f_cadet_btf,0.0653 * f_ndet_btf * topaz%c_2_n) / 
           !   (7.34e-02 / spery * max(0.0, 1.0 - f_co3_ion / co3_solubility 
           !   + 1.64 * f_ndet_btf * c_2_n * spery)**2.83 * (f_lithdet
           !   * spery + f_cadet_btf * 100.0 * spery)**-1.84 + f_cadet_btf)
           !
           ! which can be used to generate a useful initial condition. 
           !----------------------------------------------------------------
           !
           bling%fcased_redis(i,j) = min(bling%fcaco3(i,j,k),              &
             0.0653 * bling%fpop(i,j,k) * bling%c_2_p) + 7.34e-02 / spery  &
             * max(0.0, 1.0 - bling%f_co3_ion(i,j,k) /                     &
             bling%co3_solubility(i,j,k) + 1.64 * bling%fpop(i,j,k) *      &
             bling%c_2_p * spery)**(2.83) * (bling%sed_flux +              &
             bling%fcaco3(i,j,k) * 100.0 * spery)**(-1.84) *               &
             bling%f_cased(i,j,1)
           bling%fcased_burial(i,j) = max(0.0, bling%fcaco3(i,j,k) *       &
             bling%f_cased(i,j,1) / 8.1e3)
           bling%f_cased(i,j,1) = bling%f_cased(i,j,1) +                   &
             (bling%fcaco3(i,j,k) - bling%fcased_redis(i,j) -              &
             bling%fcased_burial(i,j)) / bling%z_sed * dt *                &
             grid_tmask(i,j,k)
           bling%b_alk(i,j) = - 2. * bling%fcased_redis(i,j) +             &
             bling%fpop(i,j,k) * bling%n_2_p
           ! Do not bury any C-org - all goes back to water column  
           bling%b_dic(i,j) = - bling%fpop(i,j,k) * bling%c_2_p -          &
             bling%fcased_redis(i,j)
           else
           !----------------------------------------------------------------
           !   Conserve alkalinity in the ocean
           !
           bling%b_alk(i,j) = - 2. * bling%fcaco3(i,j,k) +                 &
             bling%fpop(i,j,k) * bling%n_2_p
           ! Do not bury any C - all goes back to water column  
           bling%b_dic(i,j) = - bling%fpop(i,j,k) * bling%c_2_p -          &
             bling%fcaco3(i,j,k)
           endif
         endif  
      enddo; enddo  !} i, j

      call g_tracer_set_values(tracer_list,'alk_b', 'btf', bling%b_alk ,isd,jsd)
      call g_tracer_set_values(tracer_list,'dic_b', 'btf', bling%b_dic ,isd,jsd)

      if (do_14c) then                                        !<<RADIOCARBON
      do j = jsc, jec ; do i = isc, iec  !{
         k = grid_kmt(i,j)
         if (k .gt. 0) then !{
           bling%b_di14c(i,j) = - bling%fpo14c(i,j,k)
         endif  
      enddo; enddo  !} i, j

     call g_tracer_set_values(tracer_list,'di14c','btf',bling%b_di14c,isd,jsd)
      endif                                                   !RADIOCARBON>>
      
    endif                                                    !CARBON CYCLE>>


  !-------------------------------------------------------------------------
  !     CALCULATE SOURCE/SINK TERMS FOR EACH TRACER
  !-------------------------------------------------------------------------

    !Update the prognostics tracer fields via their pointers.

    call g_tracer_get_pointer(tracer_list,'fed_b'    ,'field',bling%p_fed    )
    call g_tracer_get_pointer(tracer_list,'dop_b'    ,'field',bling%p_dop    )
    call g_tracer_get_pointer(tracer_list,'o2_b'     ,'field',bling%p_o2     )
    call g_tracer_get_pointer(tracer_list,'po4_b'    ,'field',bling%p_po4    )

    if (do_po4_pre) &
    call g_tracer_get_pointer(tracer_list,'po4_pre','field',bling%p_po4_pre)
    
    if (do_carbon) then
    call g_tracer_get_pointer(tracer_list,'alk_b','field',bling%p_alk)
    call g_tracer_get_pointer(tracer_list,'dic_b','field',bling%p_dic)
    endif 

    if (do_14c) then
    call g_tracer_get_pointer(tracer_list,'di14c','field',bling%p_di14c)
    call g_tracer_get_pointer(tracer_list,'do14c','field',bling%p_do14c)
    endif 

    if (do_carbon_pre) then 
    call g_tracer_get_pointer(tracer_list,'alk_pre','field',bling%p_alk_pre)
    call g_tracer_get_pointer(tracer_list,'dic_pre','field',bling%p_dic_pre)
    endif

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       
       !
       ! PO4
       ! Sum of fast recycling, decay of sinking POP, and decay of DOP,
       ! less uptake.
       !    
       bling%jpo4(i,j,k) = bling%jp_recycle(i,j,k) +                       &
         (1. - bling%phi_dop) * bling%jp_reminp(i,j,k) +                   &
         (bling%gamma_dop * bling%f_dop(i,j,k))  - bling%jp_uptake(i,j,k)
  
       bling%p_po4(i,j,k,tau) = bling%p_po4(i,j,k,tau) +                   &
          bling%jpo4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! Fed
       !
       bling%p_fed(i,j,k,tau) = bling%p_fed(i,j,k,tau) +                   &
         (bling%jfe_recycle(i,j,k) + bling%jfe_reminp(i,j,k) -             &
         bling%jfe_uptake(i,j,k) - bling%jfe_ads_org(i,j,k) -              &
         bling%jfe_ads_inorg(i,j,k) ) * dt * grid_tmask(i,j,k)

       !
       ! Dissolved Organic Phosphorus
       ! Add first-order decay and generation from sinking particles
       ! to jprod, initially produced in Uptake above.
       !
       bling%jdop(i,j,k) = bling%jdop(i,j,k) + bling%phi_dop *             &
         bling%jp_reminp(i,j,k) - bling%gamma_dop * bling%f_dop(i,j,k)
         
       bling%p_dop(i,j,k,tau) = bling%p_dop(i,j,k,tau)+bling%jdop(i,j,k) * &
        dt * grid_tmask(i,j,k)

    !-----------------------------------------------------------------------
    !     O2
    ! Assuming constant P:O ratio.
    ! Optional prevention of negative oxygen (does not conserve ocean 
    ! redox potential) or alternatively it can be allowed to go negative, 
    ! keeping track of an implicit nitrate deficit 
    ! plus sulfate reduction.
    !-----------------------------------------------------------------------

       if ( (bling%prevent_neg_o2) .and.                                   &
            (bling%f_o2(i,j,k) .lt. bling%o2_min) ) then !{
         bling%jo2(i,j,k) = 0. * grid_tmask(i,j,k)
       else
         bling%jo2(i,j,k) = - bling%o2_2_p * bling%jpo4(i,j,k)              &
           * grid_tmask(i,j,k)
       endif !}

       bling%p_o2(i,j,k,tau) = bling%p_o2(i,j,k,tau) + bling%jo2(i,j,k) *  &
        dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

    if (do_carbon) then                                      !<<CARBON CYCLE
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
    
       bling%p_alk(i,j,k,tau) = bling%p_alk(i,j,k,tau) + (2. *             &
         (bling%jca_reminp(i,j,k) - bling%jca_uptake(i,j,k)) -             &
         bling%jpo4(i,j,k) * bling%n_2_p) * dt * grid_tmask(i,j,k)
         
       bling%p_dic(i,j,k,tau) = bling%p_dic(i,j,k,tau) +                   &
         (bling%jpo4(i,j,k) * bling%c_2_p + bling%jca_reminp(i,j,k) -      &
         bling%jca_uptake(i,j,k)) * dt * grid_tmask(i,j,k)
       
      if (do_14c) then                                        !<<RADIOCARBON
       bling%jdi14c(i,j,k) = (bling%jp_recycle(i,j,k) -                    &
         bling%jp_uptake(i,j,k)) * bling%c14_2_p(i,j,k) +                  &
         (1. - bling%phi_dop) * bling%j14c_reminp(i,j,k) +                 &
         (bling%gamma_dop * bling%f_do14c(i,j,k)) 
  
       bling%p_di14c(i,j,k,tau) = bling%p_di14c(i,j,k,tau) +               &
         (bling%jdi14c(i,j,k) - bling%j14c_decay_dic(i,j,k)) * dt          &
         * grid_tmask(i,j,k)

       bling%jdo14c(i,j,k) = bling%jdo14c(i,j,k) + bling%phi_dop *         &
         bling%j14c_reminp(i,j,k) - bling%gamma_dop * bling%f_do14c(i,j,k)
         
       bling%p_do14c(i,j,k,tau) = bling%p_do14c(i,j,k,tau) +               &
         (bling%jdo14c(i,j,k) - bling%j14c_decay_doc(i,j,k)) * dt          &
         * grid_tmask(i,j,k)
      endif                                                   !RADIOCARBON>>
 
    enddo; enddo ; enddo  !} i,j,k
    endif                                                    !CARBON CYCLE>>
       

    !
    ! PO4_pre
    ! Set equal to PO4 in surface layer only.
    if (do_po4_pre) then
      do j = jsc, jec ; do i = isc, iec  !{
        bling%p_po4_pre(i,j,1,tau) = bling%f_po4(i,j,1) 
      enddo; enddo ;  !} i,j
    endif
    
    !
    ! DIC_pre and ALK_pre
    ! Set equal to DIC and ALK, respectively, in surface layer.
    if (do_carbon_pre) then 
      do j = jsc, jec ; do i = isc, iec  !{
        bling%p_alk_pre(i,j,1,tau) = bling%f_alk(i,j,1) 
        bling%p_dic_pre(i,j,1,tau) = bling%f_dic(i,j,1) 
      enddo; enddo ;  !} i,j
    endif

    !
    !Set the diagnostics tracer fields.
    !
    if (use_bling_chl) then    
    call g_tracer_set_values(tracer_list,'chl',       'field',bling%f_chl       ,isd,jsd,ntau=1)
    else
      call g_tracer_set_values(tracer_list,'chl_b',       'field',bling%f_chl       ,isd,jsd,ntau=1)
    endif
    call g_tracer_set_values(tracer_list,'biomass_p', 'field',bling%f_biomass_p ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'irr_mem_b' ,  'field',bling%f_irr_mem   ,isd,jsd,ntau=1)
    if (bury_caco3) &
      call g_tracer_set_values(tracer_list,'cased_b',   'field',bling%f_cased     ,isd,jsd,ntau=1)

    !-----------------------------------------------------------------------
    !       Save variables for diagnostics
    !-----------------------------------------------------------------------
    !

    if (bling%id_alpha .gt. 0)                                                   &
         used = send_data(bling%id_alpha,          bling%alpha,                  &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_b_fed .gt. 0)                                                   &
         used = send_data(bling%id_b_fed,          bling%b_fed,                  &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_o2 .gt. 0)                                                    &
         used = send_data(bling%id_b_o2,           bling%b_o2,                   &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_po4 .gt. 0)                                                   &
         used = send_data(bling%id_b_po4,          bling%b_po4,                  &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_po4_n .gt. 0)                                                 &
         used = send_data(bling%id_b_po4_n,          bling%b_po4_n,              &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_po4_nx .gt. 0)                                                &
         used = send_data(bling%id_b_po4_nx,          bling%b_po4_nx,            &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_po4_s .gt. 0)                                                 &
         used = send_data(bling%id_b_po4_s,          bling%b_po4_s,              &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_po4_sx .gt. 0)                                                &
         used = send_data(bling%id_b_po4_sx,          bling%b_po4_sx,            &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_biomass_p_ts .gt. 0)                                            &
         used = send_data(bling%id_biomass_p_ts,   bling%biomass_p_ts,           &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_def_fe .gt. 0)                                                  &
         used = send_data(bling%id_def_fe,         bling%def_fe,                 &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_expkT .gt. 0)                                                   &
         used = send_data(bling%id_expkT,          bling%expkT,                  &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fe_2_p_uptake .gt. 0)                                           &
         used = send_data(bling%id_fe_2_p_uptake,  bling%fe_2_p_uptake,          &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_feprime .gt. 0)                                                 &
         used = send_data(bling%id_feprime,        bling%feprime,                &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fe_burial .gt. 0)                                               &
         used = send_data(bling%id_fe_burial,      bling%fe_burial,              &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_ffe_sed .gt. 0)                                                 &
         used = send_data(bling%id_ffe_sed,        bling%ffe_sed,                &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_fpofe .gt. 0)                                                   &
         used = send_data(bling%id_fpofe,          bling%fpofe,                  &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fpop .gt. 0)                                                    &
         used = send_data(bling%id_fpop,           bling%fpop,                   &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fpop_n .gt. 0)                                                  &
         used = send_data(bling%id_fpop_n,           bling%fpop_n,               &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fpop_nx .gt. 0)                                                 &
         used = send_data(bling%id_fpop_nx,           bling%fpop_nx,             &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fpop_s .gt. 0)                                                  &
         used = send_data(bling%id_fpop_s,           bling%fpop_s,               &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fpop_sx .gt. 0)                                                 &
         used = send_data(bling%id_fpop_sx,           bling%fpop_sx,             &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_frac_lg .gt. 0)                                                 &
         used = send_data(bling%id_frac_lg,        bling%frac_lg,                &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_frac_pop .gt. 0)                                                &
         used = send_data(bling%id_frac_pop,       bling%frac_pop,               &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_irr_inst .gt. 0)                                                &
         used = send_data(bling%id_irr_inst,       bling%irr_inst,               &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_irr_mix .gt. 0)                                                 &
         used = send_data(bling%id_irr_mix,        bling%irr_mix,                &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_irrk .gt. 0)                                                    &
         used = send_data(bling%id_irrk,           bling%irrk,                   &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jdop .gt. 0)                                                    &
         used = send_data(bling%id_jdop,           bling%jdop*rho_dzt,           &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jfe_ads_inorg .gt. 0)                                           &
         used = send_data(bling%id_jfe_ads_inorg,  bling%jfe_ads_inorg*rho_dzt,  &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jfe_ads_org .gt. 0)                                             &
         used = send_data(bling%id_jfe_ads_org,    bling%jfe_ads_org*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jfe_recycle .gt. 0)                                             &
         used = send_data(bling%id_jfe_recycle,    bling%jfe_recycle*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jfe_reminp .gt. 0)                                              &
         used = send_data(bling%id_jfe_reminp,     bling%jfe_reminp*rho_dzt,     &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jfe_uptake .gt. 0)                                              &
         used = send_data(bling%id_jfe_uptake,     bling%jfe_uptake*rho_dzt,     &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jo2 .gt. 0)                                                     &
         used = send_data(bling%id_jo2,            bling%jo2*rho_dzt,            &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jp_recycle .gt. 0)                                              &
         used = send_data(bling%id_jp_recycle,     bling%jp_recycle*rho_dzt,     &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jp_reminp .gt. 0)                                               &
         used = send_data(bling%id_jp_reminp,      bling%jp_reminp*rho_dzt,      &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jpn_reminp .gt. 0)                                              &
         used = send_data(bling%id_jpn_reminp,      bling%jpn_reminp*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jpnx_reminp .gt. 0)                                             &
         used = send_data(bling%id_jpnx_reminp,    bling%jpnx_reminp*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jps_reminp .gt. 0)                                              &
         used = send_data(bling%id_jps_reminp,      bling%jps_reminp*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jpsx_reminp .gt. 0)                                             &
         used = send_data(bling%id_jpsx_reminp,    bling%jpsx_reminp*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jp_uptake .gt. 0)                                               &
         used = send_data(bling%id_jp_uptake,      bling%jp_uptake*rho_dzt,      &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jpo4 .gt. 0)                                                    &
         used = send_data(bling%id_jpo4,           bling%jpo4*rho_dzt,           &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jpofe .gt. 0)                                                   &
         used = send_data(bling%id_jpofe,          bling%jpofe*rho_dzt,          &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jpop .gt. 0)                                                    &
         used = send_data(bling%id_jpop,           bling%jpop*rho_dzt,           &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_kfe_eq_lig .gt. 0)                                              &
         used = send_data(bling%id_kfe_eq_lig,     bling%kfe_eq_lig,             &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_pc_m .gt. 0)                                                    &
         used = send_data(bling%id_pc_m,           bling%pc_m,                   &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_mu .gt. 0)                                                      &
         used = send_data(bling%id_mu,         bling%mu,                         &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_pc_tot .gt. 0)                                                  &
         used = send_data(bling%id_pc_tot,         bling%pc_tot,                 &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_theta .gt. 0)                                                   &
         used = send_data(bling%id_theta,          bling%theta,                  &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_thetamax_fe .gt. 0)                                             &
         used = send_data(bling%id_thetamax_fe,    bling%thetamax_fe,            &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_wsink .gt. 0)                                                   &
         used = send_data(bling%id_wsink,          bling%wsink,                  &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_zremin .gt. 0)                                                  &
         used = send_data(bling%id_zremin,         bling%zremin,                 &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_b_alk .gt. 0)                                                   &
         used = send_data(bling%id_b_alk,          bling%b_alk,                  &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_dic .gt. 0)                                                   &
         used = send_data(bling%id_b_dic,          bling%b_dic,                  &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_co2_csurf .gt. 0)                                               &
         used = send_data(bling%id_co2_csurf,      bling%co2_csurf,              &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_co2_alpha .gt. 0)                                               &
         used = send_data(bling%id_co2_alpha,      bling%co2_alpha,              &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_co3_solubility .gt. 0)                                          &
         used = send_data(bling%id_co3_solubility, bling%co3_solubility,         &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fcaco3 .gt. 0)                                                  &
         used = send_data(bling%id_fcaco3,         bling%fcaco3,                 &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_fcaco3_to_sed .gt. 0)                                           &
         used = send_data(bling%id_fcaco3_to_sed,  bling%fcaco3_to_sed,          &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_fcased_burial .gt. 0)                                           &
         used = send_data(bling%id_fcased_burial,  bling%fcased_burial,          &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_fcased_redis .gt. 0)                                            &
         used = send_data(bling%id_fcased_redis,   bling%fcased_redis,           &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_jca_reminp .gt. 0)                                              &
         used = send_data(bling%id_jca_reminp,     bling%jca_reminp*rho_dzt,     &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jca_uptake .gt. 0)                                              &
         used = send_data(bling%id_jca_uptake,     bling%jca_uptake*rho_dzt,     &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_pco2_surf .gt. 0)                                               &
         used = send_data(bling%id_pco2_surf,      bling%pco2_surf,              &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_zremin_caco3 .gt. 0)                                            &
         used = send_data(bling%id_zremin_caco3,   bling%zremin_caco3,           &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_co2_sat_csurf .gt. 0)                                           &
         used = send_data(bling%id_co2_sat_csurf,  bling%co2_sat_csurf,          &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_pco2_sat_surf .gt. 0)                                           &
         used = send_data(bling%id_pco2_sat_surf,  bling%pco2_sat_surf,          &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_b_di14c .gt. 0)                                                 &
         used = send_data(bling%id_b_di14c,        bling%b_di14c,                &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_c14_2_p .gt. 0)                                                 &
         used = send_data(bling%id_c14_2_p,        bling%c14_2_p,                &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_c14o2_csurf .gt. 0)                                             &
         used = send_data(bling%id_c14o2_csurf,    bling%c14o2_csurf,            &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_c14o2_alpha .gt. 0)                                             &
         used = send_data(bling%id_c14o2_alpha,    bling%c14o2_alpha,            &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (bling%id_fpo14c .gt. 0)                                                  &
         used = send_data(bling%id_fpo14c,         bling%fpo14c,                 &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_j14c_decay_dic .gt. 0)                                          &
         used = send_data(bling%id_j14c_decay_dic, bling%j14c_decay_dic*rho_dzt, &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_j14c_decay_doc .gt. 0)                                          &
         used = send_data(bling%id_j14c_decay_doc, bling%j14c_decay_doc*rho_dzt, &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_j14c_reminp .gt. 0)                                             &
         used = send_data(bling%id_j14c_reminp,    bling%j14c_reminp*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jdi14c .gt. 0)                                                  &
         used = send_data(bling%id_jdi14c,         bling%jdi14c*rho_dzt,         &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bling%id_jdo14c .gt. 0)                                                  &
         used = send_data(bling%id_jdo14c,         bling%jdo14c*rho_dzt,         &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

  end subroutine generic_BLING_update_from_source

!#######################################################################
  ! <SUBROUTINE NAME="generic_BLING_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_BLING_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="SST" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Temperature   
  !  </IN>
  !  <IN NAME="SSS" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Salinity
  !  </IN>
  !  <IN NAME="rho" TYPE="real, dimension(ilb:,jlb:,:,:)">
  !   Ocean density
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  ! </SUBROUTINE>

  !User must provide the calculations for these boundary values.
  subroutine generic_BLING_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),   intent(in)   :: SST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,tau

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: sal,ST,sc_co2,sc_o2,sc_no_term,o2_saturation
    real    :: tt,tk,ts,ts2,ts3,ts4,ts5
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: o2_field
    real, dimension(:,:), ALLOCATABLE :: co2_alpha,co2_csurf,o2_alpha,o2_csurf
    real, dimension(:,:), ALLOCATABLE :: co2_sat_alpha,co2_sat_csurf,c14o2_alpha,c14o2_csurf
    character(len=fm_string_len), parameter :: sub_name = 'generic_BLING_set_boundary_values'

    !Get the necessary properties
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    call g_tracer_get_pointer(tracer_list,'o2_b' ,'field',  o2_field)

    allocate(    co2_alpha(isd:ied, jsd:jed)); co2_alpha=0.0
    allocate(    co2_csurf(isd:ied, jsd:jed)); co2_csurf=0.0
    allocate(co2_sat_alpha(isd:ied, jsd:jed)); co2_sat_alpha=0.0
    allocate(co2_sat_csurf(isd:ied, jsd:jed)); co2_sat_csurf=0.0
    allocate(  c14o2_alpha(isd:ied, jsd:jed)); c14o2_alpha=0.0
    allocate(  c14o2_csurf(isd:ied, jsd:jed)); c14o2_csurf=0.0
    allocate(     o2_alpha(isd:ied, jsd:jed)); o2_alpha=0.0
    allocate(     o2_csurf(isd:ied, jsd:jed)); o2_csurf=0.0

    do j=jsc,jec ; do i=isc,iec
       !This calculation needs an input of SST and SSS
       sal = SSS(i,j) ; ST = SST(i,j)

       !nnz: 
       !Note: In the following calculations in order to get results for o2 
       !      identical with bling code in MOM bling%Rho_0 must be replaced with rho(i,j,1,tau)
       !      This is achieved by uncommenting the following if desired.
       !! bling%Rho_0 = rho(i,j,1,tau)
       !      But since bling%Rho_0 plays the role of a unit conversion factor in this module
       !      it may be safer to keep it as a constant (1035.0) rather than the actual variable
       !      surface density rho(i,j,1,tau)
       !---------------------------------------------------------------------
       !     O2
       !---------------------------------------------------------------------
       !  Compute the oxygen saturation concentration at 1 atm total
       !  pressure in mol/kg given the temperature (t, in deg C) and
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
       !              o2_saturation = 0.282015 mol m-3
       !---------------------------------------------------------------------
       !
       tt = 298.15 - ST
       tk = 273.15 + ST
       ts = log(tt / tk)
       ts2 = ts  * ts
       ts3 = ts2 * ts
       ts4 = ts3 * ts
       ts5 = ts4 * ts

       o2_saturation = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
            exp( bling%a_0 + bling%a_1*ts + bling%a_2*ts2 + bling%a_3*ts3 + bling%a_4*ts4 + bling%a_5*ts5 + &
            (bling%b_0 + bling%b_1*ts + bling%b_2*ts2 + bling%b_3*ts3 + bling%c_0*sal)*sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the 
       !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
       !  Cycles, 12, 141-163).
       !---------------------------------------------------------------------

       sc_o2  = bling%a1_o2  + ST * (bling%a2_o2  + ST * (bling%a3_o2  + ST * bling%a4_o2 )) * &
            grid_tmask(i,j,1)
       sc_no_term = sqrt(660.0 / (sc_o2 + epsln)) 
     
       o2_alpha(i,j) = o2_saturation       * sc_no_term * bling%Rho_0
       o2_csurf(i,j) = o2_field(i,j,1,tau) * sc_no_term * bling%Rho_0 !nnz: MOM has rho(i,j,1,tau)

    enddo; enddo
    !
    !Set %csurf and %alpha for these tracers. This will mark them for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'o2_b', 'alpha',o2_alpha, isd,jsd)
    call g_tracer_set_values(tracer_list,'o2_b', 'csurf',o2_csurf, isd,jsd)

    if (do_carbon) then                                !<<CARBON CYCLE
    
    call g_tracer_get_values(tracer_list,'dic_b','alpha', co2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'dic_b','csurf', co2_csurf ,isd,jsd)

    do j=jsc,jec ; do i=isc,iec
       !This calculation needs an input of SST and SSS
       sal = SSS(i,j) ; ST = SST(i,j)
       !nnz: 
       !Note: In the following calculations in order to get results for co2 and o2 
       !      identical with bling code in MOM bling%Rho_0 must be replaced with rho(i,j,1,tau)
       !      This is achieved by uncommenting the following if desired.
       !! bling%Rho_0 = rho(i,j,1,tau)
       !      But since bling%Rho_0 plays the role of a unit conversion factor in this module
       !      it may be safer to keep it as a constant (1035.0) rather than the actual variable
       !      surface density rho(i,j,1,tau)
       !---------------------------------------------------------------------
       !     CO2
       !---------------------------------------------------------------------
       !---------------------------------------------------------------------
       !  Compute the Schmidt number of CO2 in seawater using the 
       !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
       !  7373-7382).
       !---------------------------------------------------------------------
       sc_co2 = bling%a1_co2 + ST * (bling%a2_co2 + ST * (bling%a3_co2 + ST * bling%a4_co2)) * &
            grid_tmask(i,j,1)
       sc_no_term = sqrt(660.0 / (sc_co2 + epsln)) 
       
       co2_alpha(i,j) = co2_alpha(i,j)* sc_no_term * bling%Rho_0 !nnz: MOM has rho(i,j,1,tau)  
       co2_csurf(i,j) = co2_csurf(i,j)* sc_no_term * bling%Rho_0 !nnz: MOM has rho(i,j,1,tau) 

    enddo; enddo

    call g_tracer_set_values(tracer_list,'dic_b','alpha',co2_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic_b','csurf',co2_csurf,isd,jsd)


      if (do_carbon_pre) then
 
        call g_tracer_get_values(tracer_list,'dic_sat','alpha', co2_sat_alpha ,isd,jsd)
        call g_tracer_get_values(tracer_list,'dic_sat','csurf', co2_sat_csurf ,isd,jsd)

        do j=jsc,jec ; do i=isc,iec
         !---------------------------------------------------------------------
         !     CO2_sat - schmidt number is calculated same as before, and then
         ! multiplied by 50 to get the fast gas exchange.
         ! NOTE: if unstable during initial spinup, can decrease multiplier to
         ! 1x for a year, then go back to 50
         !---------------------------------------------------------------------

         sal = SSS(i,j) ; ST = SST(i,j)
         sc_co2 = bling%a1_co2 + ST * (bling%a2_co2 + ST * (bling%a3_co2 +  &
           ST * bling%a4_co2)) * grid_tmask(i,j,1)
         sc_no_term = sqrt(660.0 / (sc_co2 + epsln)) 
       
         co2_sat_alpha(i,j) = co2_sat_alpha(i,j) * sc_no_term * bling%fast_gasex * &
           bling%Rho_0 
         co2_sat_csurf(i,j) = co2_sat_csurf(i,j) * sc_no_term * bling%fast_gasex * &
           bling%Rho_0 

        enddo; enddo

        call g_tracer_set_values(tracer_list,'dic_sat','alpha',co2_sat_alpha,isd,jsd)
        call g_tracer_set_values(tracer_list,'dic_sat','csurf',co2_sat_csurf,isd,jsd)

      endif

      if (do_14c) then
      
        call g_tracer_get_values(tracer_list,'di14c','alpha', c14o2_alpha ,isd,jsd)
        call g_tracer_get_values(tracer_list,'di14c','csurf', c14o2_csurf ,isd,jsd)

        do j=jsc,jec ; do i=isc,iec
         !---------------------------------------------------------------------
         !     14CO2 - schmidt number is calculated same as before, as is alpha.
         !   Need to multiply alpha by frac_14catm to get the real alpha.
         !   NOTE: FOR NOW, D14C fixed at 0 permil!! Need to fix this later.
         !---------------------------------------------------------------------

         sal = SSS(i,j) ; ST = SST(i,j)
         sc_co2 = bling%a1_co2 + ST * (bling%a2_co2 + ST * (bling%a3_co2 +  &
           ST * bling%a4_co2)) * grid_tmask(i,j,1)
         sc_no_term = sqrt(660.0 / (sc_co2 + epsln)) 
       
         c14o2_alpha(i,j) = c14o2_alpha(i,j) * 1. * sc_no_term * bling%Rho_0 
         c14o2_csurf(i,j) = c14o2_csurf(i,j) * sc_no_term * bling%Rho_0 

        enddo; enddo

        call g_tracer_set_values(tracer_list,'di14c','alpha',c14o2_alpha,isd,jsd)
        call g_tracer_set_values(tracer_list,'di14c','csurf',c14o2_csurf,isd,jsd)

      endif
      
    endif                                                    !CARBON CYCLE>>


    deallocate(                    &
      co2_alpha,co2_csurf,         &
      co2_sat_alpha,co2_sat_csurf, &
      c14o2_alpha,c14o2_csurf,     &
      o2_alpha,o2_csurf)

  end subroutine generic_BLING_set_boundary_values

!#######################################################################
  ! <SUBROUTINE NAME="generic_BLING_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_BLING_end
  !  </TEMPLATE>
  ! </SUBROUTINE>


  subroutine generic_BLING_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_BLING_end'
    call user_deallocate_arrays
  end subroutine generic_BLING_end

!#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Allocate all the work arrays to be used in this module.
  !
  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    allocate(bling%alpha(isd:ied, jsd:jed, 1:nk));          bling%alpha=0.0
    allocate(bling%biomass_p_ts(isd:ied, jsd:jed, 1:nk));   bling%biomass_p_ts=0.0
    allocate(bling%def_fe(isd:ied, jsd:jed, 1:nk));         bling%def_fe=0.0
    allocate(bling%expkT(isd:ied, jsd:jed, 1:nk));          bling%expkT=0.0
    allocate(bling%f_biomass_p(isd:ied, jsd:jed, 1:nk));    bling%f_biomass_p=0.0
    allocate(bling%f_chl(isd:ied, jsd:jed, 1:nk));          bling%f_chl=0.0
    allocate(bling%f_dop(isd:ied, jsd:jed, 1:nk));          bling%f_dop=0.0
    allocate(bling%f_fed(isd:ied, jsd:jed, 1:nk));          bling%f_fed=0.0
    allocate(bling%f_irr_mem(isd:ied, jsd:jed, 1:nk));      bling%f_irr_mem=0.0
    allocate(bling%f_o2(isd:ied, jsd:jed, 1:nk));           bling%f_o2=0.0
    allocate(bling%f_po4(isd:ied, jsd:jed, 1:nk));          bling%f_po4=0.0
    allocate(bling%fe_2_p_uptake(isd:ied, jsd:jed, 1:nk));  bling%fe_2_p_uptake=0.0
    allocate(bling%feprime(isd:ied, jsd:jed, 1:nk));        bling%feprime=0.0
    allocate(bling%fpofe(isd:ied, jsd:jed, 1:nk));          bling%fpofe=0.0
    allocate(bling%fpop(isd:ied, jsd:jed, 1:nk));           bling%fpop=0.0
    allocate(bling%frac_lg(isd:ied, jsd:jed, 1:nk));        bling%frac_lg=0.0
    allocate(bling%frac_pop(isd:ied, jsd:jed, 1:nk));       bling%frac_pop=0.0
    allocate(bling%irr_inst(isd:ied, jsd:jed, 1:nk));       bling%irr_inst=0.0
    allocate(bling%irr_mix(isd:ied, jsd:jed, 1:nk));        bling%irr_mix=0.0
    allocate(bling%irrk(isd:ied, jsd:jed, 1:nk));           bling%irrk=0.0
    allocate(bling%jdop(isd:ied, jsd:jed, 1:nk));           bling%jdop=0.0
    allocate(bling%jfe_ads_org(isd:ied, jsd:jed, 1:nk));    bling%jfe_ads_org=0.0
    allocate(bling%jfe_ads_inorg(isd:ied, jsd:jed, 1:nk));  bling%jfe_ads_inorg=0.0
    allocate(bling%jfe_recycle(isd:ied, jsd:jed, 1:nk));    bling%jfe_recycle=0.0
    allocate(bling%jfe_reminp(isd:ied, jsd:jed, 1:nk));     bling%jfe_reminp=0.0
    allocate(bling%jfe_uptake(isd:ied, jsd:jed, 1:nk));     bling%jfe_uptake=0.0
    allocate(bling%jo2(isd:ied, jsd:jed, 1:nk));            bling%jo2=0.0
    allocate(bling%jp_recycle(isd:ied, jsd:jed, 1:nk));     bling%jp_recycle=0.0
    allocate(bling%jp_reminp(isd:ied, jsd:jed, 1:nk));      bling%jp_reminp=0.0
    allocate(bling%jp_uptake(isd:ied, jsd:jed, 1:nk));      bling%jp_uptake=0.0
    allocate(bling%jpo4(isd:ied, jsd:jed, 1:nk));           bling%jpo4=0.0
    allocate(bling%jpofe(isd:ied, jsd:jed, 1:nk));          bling%jpofe=0.0
    allocate(bling%jpop(isd:ied, jsd:jed, 1:nk));           bling%jpop=0.0
    allocate(bling%kfe_eq_lig(isd:ied, jsd:jed, 1:nk));     bling%kfe_eq_lig=0.0
    allocate(bling%pc_m(isd:ied, jsd:jed, 1:nk));           bling%pc_m=0.0
    allocate(bling%mu(isd:ied, jsd:jed, 1:nk));         bling%mu=0.0
    allocate(bling%pc_tot(isd:ied, jsd:jed, 1:nk));         bling%pc_tot=0.0
    allocate(bling%theta(isd:ied, jsd:jed, 1:nk));          bling%theta=0.0
    allocate(bling%thetamax_fe(isd:ied, jsd:jed, 1:nk));    bling%thetamax_fe=0.0
    allocate(bling%wsink(isd:ied, jsd:jed, 1:nk));          bling%wsink=0.0
    allocate(bling%zremin(isd:ied, jsd:jed, 1:nk));         bling%zremin=0.0
    allocate(bling%zt(isd:ied, jsd:jed, 1:nk));             bling%zt=0.0
    allocate(bling%fe_burial    (isd:ied, jsd:jed));        bling%fe_burial=0.0
    allocate(bling%ffe_sed      (isd:ied, jsd:jed));        bling%ffe_sed=0.0
    allocate(bling%b_fed        (isd:ied, jsd:jed));        bling%b_fed=0.0
    allocate(bling%b_o2         (isd:ied, jsd:jed));        bling%b_o2=0.0
    allocate(bling%b_po4        (isd:ied, jsd:jed));        bling%b_po4=0.0

    if (do_po4_pre) then                                          !<<PO4_PRE
    allocate(bling%f_po4_pre(isd:ied, jsd:jed, 1:nk));      bling%f_po4_pre=0.0
    endif                                                         !PO4_PRE>>

    if (do_carbon) then                                      !<<CARBON CYCLE
    !Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc ; CO2_dope_vec%iec = iec 
    CO2_dope_vec%jsc = jsc ; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd ; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd ; CO2_dope_vec%jed = jed    
    allocate(bling%htotallo(isd:ied, jsd:jed))
    allocate(bling%htotalhi(isd:ied, jsd:jed))
    allocate(bling%co3_solubility(isd:ied, jsd:jed, 1:nk)); bling%co3_solubility=0.0
    allocate(bling%f_alk(isd:ied, jsd:jed, 1:nk));          bling%f_alk=0.0
    allocate(bling%f_co3_ion(isd:ied, jsd:jed, 1:nk));      bling%f_co3_ion=0.0
    allocate(bling%f_dic(isd:ied, jsd:jed, 1:nk));          bling%f_dic=0.0
    allocate(bling%f_htotal(isd:ied, jsd:jed, 1:nk));       bling%f_htotal=0.0
    allocate(bling%fcaco3(isd:ied, jsd:jed, 1:nk));         bling%fcaco3=0.0
    allocate(bling%jca_reminp(isd:ied, jsd:jed, 1:nk));     bling%jca_reminp=0.0
    allocate(bling%jca_uptake(isd:ied, jsd:jed, 1:nk));     bling%jca_uptake=0.0
    allocate(bling%zremin_caco3(isd:ied, jsd:jed, 1:nk));   bling%zremin_caco3=0.0
    allocate(bling%co2_csurf    (isd:ied, jsd:jed));        bling%co2_csurf=0.0
    allocate(bling%pco2_surf    (isd:ied, jsd:jed));        bling%pco2_surf=0.0
    allocate(bling%co2_alpha    (isd:ied, jsd:jed));        bling%co2_alpha=0.0
    allocate(bling%fcaco3_to_sed(isd:ied, jsd:jed));        bling%fcaco3_to_sed=0.0
    allocate(bling%b_alk        (isd:ied, jsd:jed));        bling%b_alk=0.0
    allocate(bling%b_dic        (isd:ied, jsd:jed));        bling%b_dic=0.0
      if (bury_caco3) then                                     !<<BURY CACO3  
    allocate(bling%f_cased(isd:ied, jsd:jed, 1:nk));        bling%f_cased=0.0
    allocate(bling%fcased_burial(isd:ied, jsd:jed));        bling%fcased_burial=0.0
    allocate(bling%fcased_redis (isd:ied, jsd:jed));        bling%fcased_redis=0.0
      endif                                                    !BURY CACO3>>
      if (do_carbon_pre) then                                     !<<DIC_PRE  
    allocate(bling%htotal_satlo(isd:ied, jsd:jed))
    allocate(bling%htotal_sathi(isd:ied, jsd:jed))
    allocate(bling%f_alk_pre(isd:ied, jsd:jed, 1:nk));      bling%f_alk_pre=0.0
    allocate(bling%f_dic_pre(isd:ied, jsd:jed, 1:nk));      bling%f_dic_pre=0.0
    allocate(bling%f_dic_sat(isd:ied, jsd:jed, 1:nk));      bling%f_dic_sat=0.0
    allocate(bling%f_htotal_sat(isd:ied, jsd:jed, 1:nk));   bling%f_htotal_sat=0.0
    allocate(bling%co2_sat_csurf(isd:ied, jsd:jed));        bling%co2_sat_csurf=0.0
    allocate(bling%pco2_sat_surf(isd:ied, jsd:jed));        bling%pco2_sat_surf=0.0
      endif                                                       !DIC_PRE>>
      if (do_14c) then                                        !<<RADIOCARBON
    allocate(bling%c14_2_p(isd:ied, jsd:jed, 1:nk));        bling%c14_2_p=0.0
    allocate(bling%f_di14c(isd:ied, jsd:jed, 1:nk));        bling%f_di14c=0.0
    allocate(bling%f_do14c(isd:ied, jsd:jed, 1:nk));        bling%f_do14c=0.0
    allocate(bling%fpo14c(isd:ied, jsd:jed, 1:nk));         bling%fpo14c=0.0
    allocate(bling%j14c_decay_dic(isd:ied, jsd:jed, 1:nk)); bling%j14c_decay_dic=0.0
    allocate(bling%j14c_decay_doc(isd:ied, jsd:jed, 1:nk)); bling%j14c_decay_doc=0.0
    allocate(bling%j14c_reminp(isd:ied, jsd:jed, 1:nk));    bling%j14c_reminp=0.0
    allocate(bling%jdi14c(isd:ied, jsd:jed, 1:nk));         bling%jdi14c=0.0
    allocate(bling%jdo14c(isd:ied, jsd:jed, 1:nk));         bling%jdo14c=0.0
    allocate(bling%c14o2_csurf  (isd:ied, jsd:jed));        bling%c14o2_csurf=0.0
    allocate(bling%c14o2_alpha  (isd:ied, jsd:jed));        bling%c14o2_alpha=0.0
    allocate(bling%b_di14c      (isd:ied, jsd:jed));        bling%b_di14c=0.0
      endif                                                   !RADIOCARBON>>
    endif                                                    !CARBON CYCLE>>

  end subroutine user_allocate_arrays

!#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays

    deallocate(&
         bling%alpha,&
         bling%biomass_p_ts,&
         bling%def_fe,&
         bling%expkT,&
         bling%f_biomass_p,&
         bling%f_chl,&
         bling%f_dop,&
         bling%f_fed,&
         bling%f_irr_mem,&
         bling%f_o2,&
         bling%f_po4,&
         bling%fe_2_p_uptake,&
         bling%feprime,&
         bling%fpofe,&
         bling%fpop,&
         bling%frac_lg,&
         bling%frac_pop,&
         bling%irr_inst,&
         bling%irr_mix,&
         bling%irrk,&
         bling%jdop,&
         bling%jfe_ads_inorg,&
         bling%jfe_ads_org,&
         bling%jfe_recycle,&
         bling%jfe_reminp,&
         bling%jfe_uptake,&
         bling%jo2,&
         bling%jp_recycle,&
         bling%jp_reminp,&
         bling%jp_uptake,&
         bling%jpo4,&
         bling%jpofe,&
         bling%jpop,&
         bling%kfe_eq_lig,&
         bling%pc_m,&
         bling%mu,&
         bling%pc_tot,&
         bling%theta,&
         bling%thetamax_fe,&
         bling%wsink,&
         bling%zremin,&
         bling%zt,&
         bling%fe_burial,&
         bling%ffe_sed,&
         bling%b_fed,&
         bling%b_o2,&
         bling%b_po4 )

    if (do_po4_pre) then                                          !<<PO4_PRE
    deallocate(&
         bling%f_po4_pre )
    endif                                                         !PO4_PRE>>

    if (do_carbon) then                                      !<<CARBON CYCLE
    deallocate(&
         bling%f_alk,&
         bling%co3_solubility,&
         bling%f_co3_ion,&
         bling%f_dic,&
         bling%f_htotal,&
         bling%fcaco3,&
         bling%jca_reminp,&
         bling%jca_uptake,&
         bling%zremin_caco3,&
         bling%co2_csurf,&
         bling%pco2_surf,&
         bling%co2_alpha,&
         bling%htotallo,&
         bling%htotalhi,&
         bling%fcaco3_to_sed,&
         bling%b_alk,&
         bling%b_dic )
      if (bury_caco3) then                                     !<<BURY CACO3  
      deallocate(&
           bling%f_cased,&
           bling%fcased_burial,&
           bling%fcased_redis )
      endif                                                    !BURY CACO3>>
      if (do_carbon_pre) then                                     !<<DIC_PRE  
      deallocate(&
           bling%f_alk_pre,&
           bling%f_dic_pre,&
           bling%f_dic_sat,&
           bling%f_htotal_sat,&
           bling%co2_sat_csurf,&
           bling%pco2_sat_surf,&
           bling%htotal_satlo,&
           bling%htotal_sathi )
      endif                                                       !DIC_PRE>>
      if (do_14c) then                                        !<<RADIOCARBON
      deallocate(&
           bling%c14_2_p,&
           bling%f_di14c,&
           bling%f_do14c,&
           bling%fpo14c,&
           bling%j14c_decay_dic,&
           bling%j14c_decay_doc,&
           bling%j14c_reminp,&
           bling%jdi14c,&
           bling%jdo14c,&
           bling%c14o2_alpha,&
           bling%c14o2_csurf,&
           bling%b_di14c )
      endif                                                   !RADIOCARBON>>
    endif                                                    !CARBON CYCLE>>

  end subroutine user_deallocate_arrays


end module generic_BLING

