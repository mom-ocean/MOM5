
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Charles Stock 
! </CONTACT>
!
! <OVERVIEW>
! This module contains the generic version of the COBALT 1.0 model: "Carbon Ocean
! Biogeochemistry and Lower Trophics".  COBALT augments the foodweb dynamics 
! in TOPAZ to enable anaylisis of the energy flow through the planktonic
! foodweb and improve the mechanistic resolution of foodweb dynamics that
! influence biogeochemical processes.
! </OVERVIEW>
!<DESCRIPTION>
!       COBALT simulates the biogeochemical cycling of carbon, nitrogen,
!       phosphorous, iron, silica, calcium carbonate, and lithogenic
!       material in the ocean.  The code is built upon the TOPAZ code
!       developed by John Dunne.  The primary changes to TOPAZ are:
!
!          1) the addition of three zooplankton groups
!          2) The addition of bacteria
!          3) The expansion of the dissolved organic nitrogen and 
!             phosphorous groups to include three types each: labile,
!             semi-labile, and refractory
!          4) The division of small phytoplankton into low- and high-
!             light adapted varieties
!          5) The 1.0 version of the model is coded for constant P:N.  Code
!             related to the variable P:N formulation used in TOPAZ has
!             been retained, but phytoplankton phosphorous state variables
!             have been removed (commented out) for computational savings.
!       
!       Numerous other adjustments to TOPAZ have been made and are detailed in
!       the COBALT manual, which can be found at:
!
!
!       This manual provides the rationale and justification for the various
!       parameterizations used herein, as well as definitions for all variables
!       and parameters.  The 35 model state variables are:
!
!       alk: alkalinity
!       cadet_arag: calcium carbonate detritus (aragonite)                
!       cadet_calc: calcium carbonate detritus (calcite)                  
!       dic: dissolved inorganic carbon                                   
!       fed: dissolved iron                                               
!       fedi: diazotroph iron                                             
!       felg: large phytoplankton iron
!       fedet: iron detritus                                              
!       fesm: small phytoplankton iron
!       ldon: labile dissolved organic nitrogen                           
!       ldop: labile dissolved organic phosphorous
!       lith: lithogenic aluminosilicate particles                        
!       lithdet: lithogenic detritus                                      
!       nbact: bacteria
!       ndet: nitrogen detritus                                           
!       ndi: diazotroph nitrogen                                          
!       nlg: large phyto nitrogen
!       nsm: high-light adapted small phyto nitrogen
!       nh4: ammonia                                                      
!       no3: nitrate                                                      
!       o2: oxygen                                                        
!       pdet: phosphorous detritus                                        
!       po4: phosphate                                                    
!       srdon: semi-refractory dissolved organic nitrogen
!             (decays over years to decades)
!       srdop: semi-refractory dissolved organic phosphorous
!             (decays over years to decades)
!       sldon: semi-labile dissolved organic nitrogen 
!             (decays on monthly time scales)               
!       sldop: semi-labile dissolved organic phosphorous                
!             (decays on monthly time scales)
!       sidet: silica detritus                                            
!       silg: large phyto silica
!       sio4: silicate                                                    
!       nsmz: small zooplankton nitrogen
!       nmdz: medium zooplankton nitrogen
!       nlgz: large zooplankton nitrogen
!   
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! </REFERENCE>
! <DEVELOPER_NOTES>
! </DEVELOPER_NOTES>
! </INFO>
!----------------------------------------------------------------

module generic_COBALT

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len, fm_path_name_len
  use mpp_mod,           only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod,           only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE
  use time_manager_mod,  only: time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  
  use diag_manager_mod,  only: register_diag_field, send_data 
  use constants_mod,     only: WTMCO2, WTMO2
  use fms_mod,           only: write_version_number, FATAL, WARNING, stdout, stdlog

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer
  use g_tracer_utils, only : g_tracer_get_common,g_tracer_set_common 
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  
  use g_tracer_utils, only : g_diag_type, g_diag_field_add

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector

  implicit none ; private
!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: generic_COBALT.F90,v 20.0 2013/12/14 00:18:04 fms Exp $'
  character(len=128) :: tag = '$Name: tikal $'
!-----------------------------------------------------------------------

  character(len=fm_string_len), parameter :: mod_name       = 'generic_COBALT'
  character(len=fm_string_len), parameter :: package_name   = 'generic_cobalt'

  public do_generic_COBALT
  public generic_COBALT_register
  public generic_COBALT_init
  public generic_COBALT_register_diag
  public generic_COBALT_update_from_coupler
  public generic_COBALT_update_from_source
  public generic_COBALT_update_from_bottom
  public generic_COBALT_set_boundary_values
  public generic_COBALT_end

  !The following logical for using this module is overwritten 
  logical, save :: do_generic_COBALT = .false.

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30
  real,parameter :: missing_value1=-1.0e+10
  real, parameter :: missing_value_diag=-1.0e+10

  ! Declare phytoplankton, zooplankton and cobalt variable types, which contain
  ! the vast majority of all variables used in this module. 

  type phytoplankton
     real :: alpha,   &			
          fe_2_n_max,    &
          p_2_n_static,  &
          k_fe_2_n,      &
          k_fed,         &
          k_nh4,         &
          k_no3,         &
          k_po4,         &
          k_sio4,        &
          P_C_max,       &
          si_2_n_max,    &
          si_2_n_static, &
          thetamax,      &     
          bresp,         &
          agg,           &
          vir,           &            
          exu 
     real, ALLOCATABLE, dimension(:,:)  :: &
          jprod_n_100,      & 
          jprod_n_new_100,  & 
          jprod_n_n2_100,   &   
          jzloss_n_100,     &
          jaggloss_n_100,   &
          jvirloss_n_100,   &
          jexuloss_n_100,   &
          f_n_100       
     real, ALLOCATABLE, dimension(:,:,:)  :: &
          def_fe      , & 
          def_p       , & 
          f_fe        , & 
          f_n         , & 
          felim       , & 
          irrlim      , & 
          jzloss_fe   , & 
          jzloss_n    , & 
          jzloss_p    , & 
          jzloss_sio2 , & 
          jaggloss_fe , &  
          jaggloss_n  , & 
          jaggloss_p  , &
          jaggloss_sio2,& 
          jvirloss_fe , & 
          jvirloss_n  , & 
          jvirloss_p  , & 
          jvirloss_sio2,&
          jexuloss_fe , &
          jexuloss_n  , &
          jexuloss_p  , &
          jhploss_fe  , & 
          jhploss_n   , &
          jhploss_p   , & 
          jhploss_sio2, & 
          juptake_n2  , & 
          juptake_fe  , & 
          juptake_nh4 , & 
          juptake_no3 , & 
          juptake_po4 , & 
          juptake_sio4, &
          jprod_n     , & 
          liebig_lim  , & 
          mu          , & 
          nh4lim      , & 
          no3lim      , & 
          po4lim      , &
          o2lim       , & 
          q_fe_2_n    , & 
          q_p_2_n     , & 
          silim       , & 
          q_si_2_n    , & 
          theta            
     integer ::            &
          id_def_fe       = -1, & 
          id_def_p        = -1, &
          id_felim        = -1, &
          id_irrlim       = -1, &
          id_jzloss_fe    = -1, &
          id_jzloss_n     = -1, & 
          id_jzloss_p     = -1, & 
          id_jzloss_sio2  = -1, &
          id_jaggloss_fe  = -1, &
          id_jaggloss_n   = -1, &
          id_jaggloss_p   = -1, &
          id_jaggloss_sio2= -1, & 
          id_jvirloss_fe  = -1, & 
          id_jvirloss_n   = -1, &
          id_jvirloss_p   = -1, &
          id_jvirloss_sio2= -1, &
          id_jexuloss_n   = -1, &
          id_jexuloss_p   = -1, &
          id_jexuloss_fe  = -1, &
          id_jhploss_fe   = -1, & 
          id_jhploss_n    = -1, & 
          id_jhploss_p    = -1, &
          id_jhploss_sio2 = -1, &
          id_juptake_n2   = -1, &
          id_juptake_fe   = -1, &
          id_juptake_nh4  = -1, &
          id_juptake_no3  = -1, & 
          id_juptake_po4  = -1, &
          id_juptake_sio4 = -1, &
          id_jprod_n      = -1, & 
          id_liebig_lim   = -1, &
          id_mu           = -1, &
          id_nh4lim       = -1, &
          id_no3lim       = -1, &
          id_po4lim       = -1, &
          id_o2lim        = -1, &
          id_q_fe_2_n     = -1, &
          id_q_p_2_n      = -1, &
          id_silim        = -1, &
          id_q_si_2_n     = -1, & 
          id_theta        = -1, &
          id_jprod_n_100  = -1, &
          id_jprod_n_new_100  = -1, &     
          id_jprod_n_n2_100 = -1, &
          id_jzloss_n_100     = -1, &
          id_jaggloss_n_100   = -1, &
          id_jvirloss_n_100   = -1, &
          id_jexuloss_n_100   = -1, &
          id_f_n_100          = -1, &
          id_sfc_f_n          = -1, &
          id_sfc_chl          = -1, &
          id_sfc_def_fe       = -1, &
          id_sfc_felim        = -1, &
          id_sfc_q_fe_2_n     = -1, &
          id_sfc_nh4lim       = -1, &
          id_sfc_no3lim       = -1, &
          id_sfc_po4lim       = -1, &
          id_sfc_irrlim       = -1, &
          id_sfc_theta        = -1, &
          id_sfc_mu           = -1
  end type phytoplankton

  type zooplankton
    real ::  &
	  imax,             & ! maximum ingestion rate (sec-1)         
          ki,               & ! half-sat for ingestion (moles N m-3)
          gge_max,          & ! max gross growth efficiciency (approached as i >> bresp, dimensionless)
          nswitch,          & ! switching parameter (dimensionless)
          mswitch,          & ! switching parameter (dimensionless)
          bresp,            & ! basal respiration rate (sec-1)
          ktemp,            & ! temperature dependence of zooplankton rates (C-1)
          phi_det,          & ! fraction of ingested N to detritus
          phi_ldon,         & ! fraction of ingested N/P to labile don
          phi_sldon,        & ! fraction of ingested N/P to semi-labile don
          phi_srdon,        & ! fraction of ingested N/P to semi-refractory don
          phi_ldop,         & ! fraction of ingested N/P to labile dop
          phi_sldop,        & ! fraction of ingested N/P to semi-labile dop
          phi_srdop,        & ! fraction of ingested N/P to semi-refractory dop 
          phi_nh4,          & ! fraction of ingested N to nh4 due to ingestion-related metabolism
          phi_po4,	    & ! fraction of ingested N to po4 due to ingestion-related metabolism
          q_p_2_n,          & ! p:n ratio of zooplankton
          ipa_smp,          & ! innate prey availability of low-light adapt. small phytos 
          ipa_lgp,          & ! innate prey availability of large phytoplankton
          ipa_diaz,         & ! innate prey availability of diazotrophs 
          ipa_smz,          & ! innate prey availability of small zooplankton
          ipa_mdz,          & ! innate prey availability of large zooplankton
          ipa_lgz,          & ! innate prey availability of x-large zooplankton
          ipa_det,          & ! innate prey availability of detritus
          ipa_bact            ! innate prey availability for bacteria
    real, ALLOCATABLE, dimension(:,:)  :: &
          jprod_n_100,      &
          jingest_n_100,    &
          jzloss_n_100,     &
          jhploss_n_100,    &
          jprod_ndet_100,   &
          jprod_don_100,    &
          jremin_n_100,     &
          f_n_100          
    real, ALLOCATABLE, dimension(:,:,:) :: &
          f_n,              & ! zooplankton biomass
          jzloss_n,         & ! Losses of n due to consumption by other zooplankton groups
          jzloss_p,	    & ! Losses of p due to consumption by other zooplankton groups 
          jhploss_n,        & ! Losses of n due to consumption by unresolved higher preds
          jhploss_p,	    & ! Losses of p due to consumption by unresolved higher preds
          jingest_n,        & ! Total ingestion of n
          jingest_p,        & ! Total ingestion of p
          jingest_sio2,     & ! Total ingestion of silicate
          jingest_fe,	    & ! Total ingestion of iron
          jprod_ndet,       & ! production of nitrogen detritus by zooplankton group 
          jprod_pdet,       & ! production of phosphorous detritus by zooplankton group
          jprod_ldon,       & ! production of labile dissolved organic N by zooplankton group
          jprod_ldop,       & ! production of labile dissolved organic P by zooplankton group
          jprod_srdon,      & ! production of semi-refractory dissolved organic N by zooplankton group
          jprod_srdop,      & ! production of semi-refractory dissolved organic P by zooplankton group 
          jprod_sldon,      & ! production of semi-labile dissolved organic N by zooplankton group
          jprod_sldop,      & ! production of semi-labile dissolved organic P by zooplankton group
          jprod_fedet,      & ! production of iron detritus
          jprod_fed,	    & ! production of dissolved iron
          jprod_sidet,	    & ! production of silica detritus
          jprod_sio4,       & ! production of silicate via rapid dissolution at surface
          jprod_po4,        & ! phosphate production by zooplankton
          jprod_nh4,        & ! ammonia production by zooplankton
          jprod_n,          & ! zooplankton production
          temp_lim            ! Temperature limitation
    integer ::		    &
          id_jzloss_n       = -1, &
          id_jzloss_p       = -1, &
          id_jhploss_n      = -1, &
          id_jhploss_p      = -1, &
          id_jingest_n      = -1, &
          id_jingest_p      = -1, &
          id_jingest_sio2   = -1, &
          id_jingest_fe     = -1, &
          id_jprod_ndet     = -1, &
          id_jprod_pdet     = -1, &
          id_jprod_ldon     = -1, &
          id_jprod_ldop     = -1, &
          id_jprod_srdon    = -1, &
          id_jprod_srdop    = -1, &
          id_jprod_sldon    = -1, &
          id_jprod_sldop    = -1, &
          id_jprod_fedet    = -1, &
          id_jprod_fed      = -1, &
          id_jprod_sidet    = -1, &
          id_jprod_sio4     = -1, &
          id_jprod_po4      = -1, &
          id_jprod_nh4      = -1, &
          id_jprod_n        = -1, &
          id_temp_lim       = -1, &
          id_jprod_n_100    = -1, &
          id_jingest_n_100  = -1, &
          id_jzloss_n_100   = -1, &
          id_jhploss_n_100  = -1, &
          id_jprod_ndet_100 = -1, &
          id_jprod_don_100  = -1, &
          id_jremin_n_100   = -1, &
          id_f_n_100        = -1
  end type zooplankton

  type bacteria
    real ::  &
          mu_max,           & ! maximum bacterial growth rate (sec-1)
          k_ldon,           & ! half-sat for nitrogen-limited growth (mmoles N m-3)
          gge_max,          & ! max gross growth efficiciency (dimensionless)
          bresp,            & ! basal respiration rate (sec-1)
          ktemp,            & ! temperature dependence of bacterial rates (C-1)
          vir,              & ! virus-driven loss rate for bacteria (sec-1 mmole N m-3)
          q_p_2_n             ! p:n ratio for bacteria 
    real, ALLOCATABLE, dimension(:,:)  :: &
          jprod_n_100,      &
          jzloss_n_100,     &
          jvirloss_n_100,   &
          jremin_n_100,     &
          juptake_ldon_100, &
          f_n_100
    real, ALLOCATABLE, dimension(:,:,:) :: &
          f_n,              & ! bacteria biomass
          jzloss_n,         & ! Losses of n due to consumption by zooplankton 
          jzloss_p,         & ! Losses of p due to consumption by zooplankton
          jhploss_n,        & ! Losses of n due to consumption by unresolved higher preds
          jhploss_p,        & ! Losses of p due to consumption by unresolved higher preds
          jvirloss_n  ,     & ! nitrogen losses via viruses
          jvirloss_p  ,     & ! phosphorous losses via viruses
          juptake_ldon,     & ! Total uptake of ldon
          juptake_ldop,     & ! Total uptake of sldon
          jprod_nh4,        & ! production of ammonia bacteria  
          jprod_po4,        & ! production of phosphate by bacteria 
          jprod_n,          & ! bacterial production
          temp_lim            ! Temperature limitation
    integer ::              &
          id_jzloss_n       = -1, &
          id_jzloss_p       = -1, &
          id_jhploss_n      = -1, &
          id_jhploss_p      = -1, &
          id_jvirloss_n     = -1, &
          id_jvirloss_p     = -1, &
          id_juptake_ldon   = -1, &
          id_juptake_ldop   = -1, &
          id_jprod_nh4      = -1, &
          id_jprod_po4      = -1, &
          id_jprod_n        = -1, &
          id_temp_lim       = -1, &
          id_jprod_n_100    = -1, &
          id_jzloss_n_100   = -1, &
          id_jvirloss_n_100 = -1, &
          id_jremin_n_100   = -1, &
          id_juptake_ldon_100 = -1, &
          id_f_n_100
  end type bacteria

  integer, parameter :: NUM_PHYTO  = 3
  !
  ! Array allocations and flux calculations assume that phyto(1) is the
  ! only phytoplankton group cabable of nitrogen uptake by N2 fixation while phyto(2:NUM_PHYTO) 
  ! are only cabable of nitrgen uptake by NH4 and NO3 uptake
  !
  integer, parameter :: DIAZO      = 1
  integer, parameter :: LARGE      = 2
  integer, parameter :: SMALL      = 3
  type(phytoplankton), dimension(NUM_PHYTO), save :: phyto

  ! define three zooplankton types
  integer, parameter :: NUM_ZOO = 3
  type(zooplankton), dimension(NUM_ZOO), save :: zoo

  type(bacteria), dimension(1), save :: bact

  integer, parameter :: NUM_PREY = 8

  type generic_COBALT_type

     logical  ::       &
          init,             &                  ! If tracers should be initializated
          p_2_n_static,     &                  ! If P:N is fixed in phytoplankton
          tracer_debug

     real  ::          &
          atm_co2_flux,     &
          c_2_n,            & 
          ca_2_n_arag,      &
          ca_2_n_calc,      & 
          caco3_sat_max,    &
          fe_2_n_upt_fac,   &
          fe_2_n_sed,       &
          fe_coast,         &
          felig_2_don,      &
          felig_bkg ,       &
          gamma_cadet_arag, & 
          gamma_cadet_calc, & 
          gamma_irr_mem,    &
          gamma_ndet,       &
          gamma_nitrif,     &
          gamma_sidet,      &
          gamma_srdon,      &
          gamma_srdop,      &
          gamma_sldon,      &
          gamma_sldop,      &
          irr_inhibit,      &
          k_n_inhib_di,     &
          k_o2,             &
          kappa_eppley,     &
          kappa_remin,      &
          kfe_eq_lig_hl,    &
          kfe_eq_lig_ll,    &
          alpha_fescav,     &
          beta_fescav,      &
          gamma_fescav,     &
          ki_fescav,        & 
          io_fescav,        &
          remin_eff_fedet,  &
          k_lith,           &
          phi_lith,         &
          mass_2_n,         &
          alk_2_n_denit,    &
          n_2_n_denit,      &
          k_no3_denit,      &
          o2_min,           &
          o2_2_c,           &
          o2_2_nfix,        & 
          o2_2_nh4,         &
          o2_2_no3,         &
          o2_2_nitrif,      &
          o2_inhib_di_pow,  &
          o2_inhib_di_sat,  &
          P_C_max_assem,    &
          rpcaco3,          &
          rplith,           &
          rpsio2,           &
          thetamin,         &
          thetamin_nolim,   &
          vir_ktemp,        &
          lysis_phi_ldon,   &
          lysis_phi_srdon,  &
          lysis_phi_sldon,  & 
          lysis_phi_ldop,   & 
          lysis_phi_srdop,  &
          lysis_phi_sldop,  &
          wsink,            &
          z_sed,            &
          zeta,             &
          imax_hp,          & ! unresolved higher pred. max ingestion rate
          ki_hp,            & ! unresolved higher pred. half-sat
          ktemp_hp,         & ! temperature dependence for higher predators
          coef_hp,          & ! scaling between unresolved preds and available prey
          nswitch_hp,	    & ! higher predator switching behavior
          mswitch_hp,       & ! higher predator switching behavior
          hp_ipa_smp,       & ! innate prey availability of small phytos to hp's
          hp_ipa_lgp,       & ! "  "  "  "  "  "  "  "  "   large phytos to hp's
          hp_ipa_diaz,      & ! "  "  "  "  "  "  "  "  "   diazotrophs to hp's  
          hp_ipa_bact,      & ! "  "  "  "  "  "  "  "  "   bacteria to hp's
          hp_ipa_smz,       & ! "  "  "  "  "  "  "  "  "   small zooplankton to hp's
          hp_ipa_mdz,       & ! "  "  "  "  "  "  "  "  "   medium zooplankton to hp's
          hp_ipa_lgz,       & ! "  "  "  "  "  "  "  "  "   large zooplankton to hp's
          hp_ipa_det,       & ! "  "  "  "  "  "  "  "  "   detritus to hp's
          hp_phi_det,       & ! fraction of ingested N to detritus
          hp_phi_ldon,      & ! fraction of ingested N to labile don
          hp_phi_sldon,     & ! fraction of ingested N to semi-labile don
          hp_phi_srdon,     & ! fraction of ingested N to semi-refractory don
          hp_phi_ldop,      & ! fraction of ingested N to labile dop
          hp_phi_sldop,     & ! fraction of ingested N to semi-labile dop
          hp_phi_srdop,     & ! fraction of ingested N to semi-refractory dop
          hp_phi_nh4,       & ! fraction of ingested N to nh4 due to ingestion-related metabolism
          hp_phi_po4          ! fraction of ingested N to po4 due to ingestion-related metabolism

          
     real, dimension(3)                    :: total_atm_co2

     real    :: htotal_scale_lo, htotal_scale_hi, htotal_in
     real    :: Rho_0, a_0, a_1, a_2, a_3, a_4, a_5, b_0, b_1, b_2, b_3, c_0
     real    :: a1_co2, a2_co2, a3_co2, a4_co2, a1_o2, a2_o2, a3_o2, a4_o2

     logical, dimension(:,:), ALLOCATABLE ::  &
          mask_z_sat_arag,&
          mask_z_sat_calc

     real, dimension(:,:,:), ALLOCATABLE ::  &
          f_alk,&				! Other prognostic variables
          f_cadet_arag,&
          f_cadet_calc,&
          f_dic,&
          f_fed,&
          f_fedet,&
          f_ldon,&
          f_ldop,&
          f_lith,&
          f_lithdet,&
          f_ndet,&
          f_nh4,&
          f_no3,&
          f_o2,&
          f_pdet,&
          f_po4,&
          f_srdon,&
          f_srdop,&
          f_sldon,&
          f_sldop,&
          f_sidet,&
          f_silg,&
          f_sio4,&
          co3_sol_arag,&
          co3_sol_calc,&
          f_chl,&
          f_co3_ion,&
          f_htotal,&
          f_irr_mem,&       
          f_cased,&
          f_cadet_arag_btf,&
          f_cadet_calc_btf,&
          f_fedet_btf, &
          f_lithdet_btf, &
          f_ndet_btf,&
          f_pdet_btf,&
          f_sidet_btf,&
          jnbact,&
          jndi,&
          jnsm,&
          jnlg,&
          jnsmz,&
          jnmdz,&
          jnlgz,&
          jalk,&
          jcadet_arag,&
          jcadet_calc,&
          jdic,&
          jfed,&
          jfedi,&
          jfelg,&
          jfesm,&
          jfedet,&
          jldon,&
          jldop,&
          jlith,&
          jlithdet,&
          jndet,&
          jnh4,&
          jno3,&
          jo2,&
          jpdet,&
          jpo4,&
          jsrdon,&
          jsrdop,&
          jsldon,&
          jsldop,&
          jsidet,&
          jsilg,&
          jsio4,&                 
          jprod_ndet,&
          jprod_pdet,&
          jprod_ldon,&
          jprod_ldop,&
          jprod_sldon,&
          jprod_sldop,&
          jprod_srdon,&
          jprod_srdop,&
          jprod_fedet,&
          jprod_fed,&
          jprod_sidet,&
          jprod_sio4, &
          jprod_lithdet,&
          jprod_cadet_arag,&
          jprod_cadet_calc,&
          jprod_nh4,&
          jprod_po4,&
          det_jzloss_n,&
          det_jzloss_p,&
          det_jzloss_fe,&
          det_jzloss_si,&
          det_jhploss_n,&
          det_jhploss_p,&
          det_jhploss_fe,&
          det_jhploss_si,&
          jdiss_cadet_arag,&
          jdiss_cadet_calc,&
          jdiss_sidet,&
          jremin_ndet,&
          jremin_pdet,&
          jremin_fedet,&
          jfe_ads,&
          jfe_coast,&
          kfe_eq_lig,&
          expkT,&
          hp_temp_lim,&
          hp_jingest_n,&
          hp_jingest_p,&
          hp_jingest_fe,&
          hp_jingest_sio2,&
          irr_inst,&	
          irr_mix,&
          jno3denit_wc,&
          jnitrif,&
          omega_arag,&
          omega_calc,&                                                  
          tot_layer_int_c,&
          tot_layer_int_fe,&
          tot_layer_int_n,&
          tot_layer_int_p,&
          tot_layer_int_si,&
          total_filter_feeding,&
          net_prim_prod,&
          gross_prim_prod,&
          nlg_diatoms,&
          q_si_2_n_lg_diatoms,&
          zt, &
          zm

     real, dimension(:,:), ALLOCATABLE :: &
          b_alk,b_dic,b_fed,b_nh4,b_no3,b_o2,b_po4,b_sio4,&	! bottom flux terms
          co2_csurf,pco2_csurf,co2_alpha,&
          fcadet_arag_btm,&
          fcadet_calc_btm,&
          ffedet_btm,&
          flithdet_btm,&
          fpdet_btm,&
          fndet_btm,&
          fsidet_btm,&      
          fcased_burial,&
          fcased_input,&
          fcased_redis,&
          ffe_sed,&
          fnfeso4red_sed,&
          fno3denit_sed,&
          fnoxic_sed,&
          frac_burial,&
          fndet_burial,&
          fpdet_burial,&
          jprod_allphytos_100,&
          htotallo, htotalhi,&
          hp_jingest_n_100,&
          hp_jremin_n_100,&
          hp_jprod_ndet_100,&
          jprod_lithdet_100,&
          jprod_sidet_100,&
          jprod_cadet_calc_100,&
          jprod_cadet_arag_100,&
          jprod_mesozoo_200, &
          jremin_ndet_100, &
          f_ndet_100, &
          f_don_100, &
          f_silg_100, &
          f_mesozoo_200, &
          fndet_100, &
          fpdet_100, &
          fsidet_100, &
          fcadet_calc_100, &
          fcadet_arag_100, &
          ffedet_100, &
          flithdet_100, &
          btm_temp,     &
          btm_o2,       &
          o2min, & 
          z_o2min, & 
          z_sat_arag,&
          z_sat_calc

     real, dimension(:,:,:,:), pointer :: &
          p_alk,&
          p_cadet_arag,&
          p_cadet_calc,&
          p_dic,&
          p_fed,&
          p_fedi,&
          p_felg,&
          p_fedet,&
          p_fesm,&
          p_ldon,&
          p_ldop,&
          p_lith,&
          p_lithdet,&          
          p_nbact,&
          p_ndet,&
          p_ndi,&
          p_nlg,&
          p_nsm,&
          p_nh4,&
          p_no3,&
          p_o2,&
          p_pdet,&
          p_po4,&
          p_srdon,&
          p_srdop,&
          p_sldon,&
          p_sldop,&
          p_sidet,&
          p_silg,&
          p_sio4,&
          p_nsmz,&
          p_nmdz,&
          p_nlgz

      real, dimension (:,:), pointer :: &
          runoff_flux_alk,&
          runoff_flux_dic,&
          runoff_flux_lith,&
          runoff_flux_fed,&
          runoff_flux_no3,&
          runoff_flux_ldon,&
          runoff_flux_sldon,&
          runoff_flux_srdon,&
          runoff_flux_ndet,&
          runoff_flux_po4,&
          runoff_flux_ldop,&
          runoff_flux_sldop,&
          runoff_flux_srdop,&
          dry_fed, wet_fed,&
          dry_lith, wet_lith,&
          dry_no3, wet_no3,&
          dry_nh4, wet_nh4,&
          dry_po4, wet_po4

     integer :: nkml
     character(len=fm_string_len)          :: file
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file

     integer               ::          &
          id_ndi           = -1,       &
          id_nlg           = -1,       &
          id_nsm           = -1,       &
          id_nsmz          = -1,       &
          id_nmdz          = -1,       &
          id_nlgz          = -1,       & 
          id_nbact         = -1,       &
          id_alk           = -1,       &
          id_cadet_arag    = -1,       &
          id_cadet_calc    = -1,       &
          id_dic           = -1,       &
          id_fed           = -1,       &
          id_fedi          = -1,       &
          id_felg          = -1,       &
          id_fesm          = -1,       &
          id_fedet         = -1,       &
          id_ldon          = -1,       &
          id_ldop          = -1,       &
          id_lith          = -1,       &
          id_lithdet       = -1,       &
          id_ndet          = -1,       &
          id_nh4           = -1,       &
          id_no3           = -1,       &
          id_o2            = -1,       &
          id_pdet          = -1,       & 
          id_po4           = -1,       &
          id_srdop         = -1,       &
          id_srdon         = -1,       &
          id_sldon         = -1,       &
          id_sldop         = -1,       &
          id_sidet         = -1,       &
          id_silg          = -1,       &
          id_sio4          = -1,       &
          id_co3_sol_arag  = -1,       &
          id_co3_sol_calc  = -1,       &
          id_dep_dry_fed   = -1,       &
          id_dep_dry_nh4   = -1,       & 
          id_dep_dry_no3   = -1,       &
          id_dep_dry_po4   = -1,       & 
          id_dep_wet_fed   = -1,       & 
          id_dep_wet_nh4   = -1,       &
          id_dep_wet_no3   = -1,       &
          id_dep_wet_po4   = -1,       &
          id_dep_wet_lith  = -1,       &
          id_dep_dry_lith  = -1,       &
          id_omega_arag    = -1,       &
          id_omega_calc    = -1,       &
          id_chl           = -1,       &
          id_co3_ion       = -1,       &
          id_htotal        = -1,       &
          id_irr_mem       = -1,       &
          id_cased         = -1,       &
	  id_cadet_arag_btf = -1,      & 
          id_cadet_calc_btf = -1,      &
          id_fedet_btf     = -1,       & 
          id_lithdet_btf   = -1,       & 
          id_ndet_btf      = -1,       & 
          id_pdet_btf      = -1,       & 
          id_sidet_btf     = -1,       &
          id_jprod_ndet    = -1,       &
          id_jprod_pdet    = -1,       &
          id_jprod_sldon   = -1,       &
          id_jprod_ldon    = -1,       &
          id_jprod_srdon   = -1,       &
          id_jprod_sldop   = -1,       &
          id_jprod_ldop    = -1,       &
          id_jprod_srdop   = -1,       &
          id_jprod_fedet   = -1,       &
          id_jprod_fed     = -1,       &
          id_jprod_sidet   = -1,       &
          id_jprod_sio4    = -1,       &
          id_jprod_lithdet = -1,       &
          id_jprod_cadet_arag = -1,    &
          id_jprod_cadet_calc = -1,    & 
          id_jprod_po4     = -1,       &
          id_jprod_nh4     = -1,       &
          id_det_jzloss_n  = -1,       &
          id_det_jzloss_p  = -1,       &
          id_det_jzloss_fe = -1,       &
          id_det_jzloss_si = -1,       &
          id_det_jhploss_n = -1,       &
          id_det_jhploss_p = -1,       &
          id_det_jhploss_fe = -1,      &
          id_det_jhploss_si = -1,      &
          id_jdiss_sidet   = -1,       &
          id_jdiss_cadet_arag = -1,    &
          id_jdiss_cadet_calc = -1,    &
          id_jremin_ndet   = -1,       &
          id_jremin_pdet   = -1,       & 
          id_jremin_fedet  = -1,       &
          id_jfe_ads       = -1,       &
          id_jfe_coast     = -1,       &
          id_kfe_eq_lig    = -1,       &
          id_expkT         = -1,       &
          id_hp_temp_lim   = -1,       &
          id_hp_jingest_n  = -1,       &
          id_hp_jingest_p  = -1,       &
          id_hp_jingest_fe = -1,       &
          id_hp_jingest_sio2 = -1,     &                
          id_irr_inst      = -1,       &
          id_irr_mix       = -1,       &
          id_jno3denit_wc  = -1,       &
          id_jnitrif       = -1,       &
          id_co2_csurf     = -1,       & 
          id_pco2_csurf    = -1,       &
          id_co2_alpha     = -1,       &
          id_fcadet_arag   = -1,       &
          id_fcadet_calc   = -1,       &
          id_ffedet        = -1,       &
          id_fndet         = -1,       &
          id_fpdet         = -1,       &
          id_fsidet        = -1,       & 
          id_flithdet      = -1,       &
          id_fcadet_arag_btm = -1,     &
          id_fcadet_calc_btm = -1,     &
          id_ffedet_btm    = -1,       &
          id_flithdet_btm  = -1,       &
          id_fndet_btm     = -1,       &
          id_fpdet_btm     = -1,       &
          id_fsidet_btm    = -1,       &
          id_fcased_burial = -1,       &
          id_fcased_input  = -1,       &
          id_fcased_redis  = -1,       &
          id_ffe_sed       = -1,       &
          id_fnfeso4red_sed= -1,       &
          id_fno3denit_sed = -1,       &
          id_fnoxic_sed    = -1,       &
          id_frac_burial   = -1,       &
          id_fndet_burial  = -1,       &
          id_fpdet_burial  = -1,       &
          id_nphyto_tot    = -1,       &
          id_no3_in_source = -1,       &
          id_pco2surf      = -1,       &
          id_sfc_alk       = -1,       &
          id_sfc_cadet_arag= -1,       & 
          id_sfc_cadet_calc= -1,       & 
          id_sfc_dic       = -1,       & 
          id_sfc_fed       = -1,       & 
          id_sfc_ldon      = -1,       &
          id_sfc_sldon     = -1,       &
          id_sfc_srdon     = -1,       &
          id_sfc_no3       = -1,       &
          id_sfc_nh4       = -1,       &
          id_sfc_po4       = -1,       &
          id_sfc_sio4      = -1,       &
          id_sfc_htotal    = -1,       &
          id_sfc_o2        = -1,       &
          id_sfc_chl       = -1,       &
          id_sfc_irr       = -1,       &
          id_sfc_irr_mem   = -1,       &
          id_sfc_temp      = -1,       &
          id_btm_temp      = -1,       &
          id_btm_o2        = -1,       &
          id_sfc_co3_ion   = -1,       &
          id_sfc_co3_sol_arag = -1,    &
          id_sfc_co3_sol_calc = -1,    &
          id_runoff_flux_alk = -1,     &
          id_runoff_flux_dic = -1,     &
          id_runoff_flux_fed = -1,     &
          id_runoff_flux_lith = -1,    &
          id_runoff_flux_no3 = -1,     &
          id_runoff_flux_ldon = -1,    &
          id_runoff_flux_sldon = -1,   &
          id_runoff_flux_srdon = -1,   &
          id_runoff_flux_ndet = -1,    &
          id_runoff_flux_po4 = -1,     &
          id_runoff_flux_ldop = -1,    &
          id_runoff_flux_sldop = -1,   &
          id_runoff_flux_srdop = -1,   &
          id_tot_layer_int_c = -1,     & 
          id_tot_layer_int_fe = -1,    & 
          id_tot_layer_int_n = -1,     & 
          id_tot_layer_int_p = -1,     & 
          id_tot_layer_int_si = -1,    & 
          id_total_filter_feeding = -1,&
          id_net_prim_prod = -1,       &
          id_gross_prim_prod = -1,     &
          id_nlg_diatoms = -1,         &
          id_jprod_allphytos_100 = -1, &
          id_q_si_2_n_lg_diatoms = -1, &
          id_hp_jingest_n_100 = -1,    &
          id_hp_jremin_n_100 = -1,     &
          id_hp_jprod_ndet_100 = -1,   &
          id_jprod_lithdet_100 = -1,   &
          id_jprod_sidet_100 = -1,     &
          id_jprod_cadet_calc_100 = -1, &
          id_jprod_cadet_arag_100 = -1, &
          id_jprod_mesozoo_200 = -1,   &
          id_jremin_ndet_100 = -1,     &
          id_f_ndet_100 = -1,          &
          id_f_don_100 = -1,           &
          id_f_silg_100 = -1,          &
          id_f_mesozoo_200 = -1,       &
          id_fndet_100 = -1,           &
          id_fpdet_100 = -1,           &
          id_ffedet_100 = -1,          &
          id_fcadet_calc_100 = -1,     &
          id_fcadet_arag_100 = -1,     &
          id_flithdet_100 = -1,        &
          id_fsidet_100 = -1,          &
          id_o2min         = -1,       &
          id_z_o2min       = -1,       &
          id_z_sat_arag    = -1,       & ! Depth of Aragonite saturation
          id_z_sat_calc    = -1          ! Depth of Calcite saturation
  end type generic_COBALT_type

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

  type(generic_COBALT_type), save :: cobalt 
  
  type(CO2_dope_vector) :: CO2_dope_vec

  ! identification numbers for mpp clocks
  integer :: id_clock_carbon_calculations
  integer :: id_clock_phyto_growth
  integer :: id_clock_bacteria_growth
  integer :: id_clock_zooplankton_calculations
  integer :: id_clock_other_losses
  integer :: id_clock_production_loop
  integer :: id_clock_ballast_loops
  integer :: id_clock_source_sink_loop1
  integer :: id_clock_source_sink_loop2
  integer :: id_clock_source_sink_loop3
  integer :: id_clock_source_sink_loop4
  integer :: id_clock_source_sink_loop5
  integer :: id_clock_source_sink_loop6
  integer :: id_clock_cobalt_send_diagnostics
  integer :: id_clock_cobalt_calc_diagnostics

contains

  subroutine generic_COBALT_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_register'
  
    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)
    
  end subroutine generic_COBALT_register

  !  <SUBROUTINE NAME="generic_COBALT_init">
  !  <OVERVIEW>
  !   Initialize the generic COBALT module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the COBALT Tracers to the list of generic Tracers
  !       passed to it via utility subroutine g_tracer_add().
  !
  !       Adds all the parameters used by this module via utility
  !       subroutine g_tracer_add_param().
  !
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_COBALT_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate all the private work arrays used by this module.
    call user_allocate_arrays

    id_clock_carbon_calculations = mpp_clock_id('(Cobalt: carbon calcs)' ,grain=CLOCK_MODULE)
    id_clock_phyto_growth = mpp_clock_id('(Cobalt: phytoplankton growth calcs)',grain=CLOCK_MODULE)
    id_clock_bacteria_growth = mpp_clock_id('(Cobalt: bacteria growth calcs)',grain=CLOCK_MODULE)
    id_clock_zooplankton_calculations = mpp_clock_id('(Cobalt: zooplankton calculations)',grain=CLOCK_MODULE)
    id_clock_other_losses = mpp_clock_id('(Cobalt: other losses)',grain=CLOCK_MODULE)
    id_clock_production_loop = mpp_clock_id('(Cobalt: production loop)',grain=CLOCK_MODULE)
    id_clock_ballast_loops = mpp_clock_id('(Cobalt: ballasting loops)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop1 = mpp_clock_id('(Cobalt: source/sink loop 1)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop2 = mpp_clock_id('(Cobalt: source/sink loop 2)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop3 = mpp_clock_id('(Cobalt: source/sink loop 3)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop4 = mpp_clock_id('(Cobalt: source/sink loop 4)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop5 = mpp_clock_id('(Cobalt: source/sink loop 5)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop6 = mpp_clock_id('(Cobalt: source/sink loop 6)',grain=CLOCK_MODULE)
    id_clock_cobalt_send_diagnostics = mpp_clock_id('(Cobalt: send diagnostics)',grain=CLOCK_MODULE)
    id_clock_cobalt_calc_diagnostics = mpp_clock_id('(Cobalt: calculate diagnostics)',grain=CLOCK_MODULE)

  end subroutine generic_COBALT_init

  !   Register diagnostic fields to be used in this module. 
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !
  subroutine generic_COBALT_register_diag(diag_list)
    type(g_diag_type), pointer :: diag_list
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


    ! Register the diagnostics for the various phytoplankton 
    !
    ! Register Limitation Diagnostics
    !
    vardesc_temp = vardesc("def_fe_Di","Diaz. Phyto. Fe Deficiency",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("def_fe_Lg","Large Phyto. Fe Deficiency",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("def_fe_Sm","Small Phyto. Fe Deficiency",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Di","Diaz. Phyto. Fed uptake Limitation",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Lg","Large Phyto. Fed uptake Limitation",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Sm","Small Phyto. Fed uptake Limitation",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Di","Diaz. Phyto. Light Limitation",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Lg","Large Phyto. Light Limitation",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Sm","Small Phyto. Light Limitation",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("theta_Di","Diaz. Phyto. Chl:C",'h','L','s','g Chl (g C)-1','f')
    phyto(DIAZO)%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("theta_Lg","Large Phyto. Chl:C",'h','L','s','g Chl (g C)-1','f')
    phyto(LARGE)%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("theta_Sm","Small Phyto. Chl:C",'h','L','s','g Chl (g C)-1','f')
    phyto(SMALL)%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Di","Diaz. Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(DIAZO)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Lg","Large Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(LARGE)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Sm","Small Phyto. Growth Rate",'h','L','s','s-1','f')
    phyto(SMALL)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Lg","Ammonia Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Sm","Ammonia Limitation of Small Phyto",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Lg","Nitrate Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Sm","Nitrate Limitation of Small Phyto",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Di","Phosphate Limitation of Diaz. Phyto",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Lg","Phosphate Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Sm","Phosphate Limitation of Small Phyto",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2lim_Di","Oxygen Limitation of Diaz. Phyto",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Di","Fe:N ratio of Diaz. Phyto",'h','L','s','mol Fe/mol N','f')
    phyto(DIAZO)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Lg","Fe:N ratio of Large Phyto",'h','L','s','mol Fe/mol N','f')
    phyto(LARGE)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Sm","Fe:N ratio of Small Phyto",'h','L','s','mol Fe/mol N','f')
    phyto(SMALL)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("silim_Lg","SiO4 Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_silim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_si_2_n_Lg","Si:N ratio of Large Phyto",'h','L','s','mol Si/mol N','f')
    phyto(LARGE)%id_q_si_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for phytoplankton loss terms: zooplankton 
    ! CAS: loss diagnostics simplified to just N

    vardesc_temp = vardesc("jzloss_n_Di","Diazotroph nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Lg","Large phyto nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Sm","Small phyto nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    !  Register diagnostics for phytoplankton loss terms: aggregation 
    !

    vardesc_temp = vardesc("jaggloss_n_Di","Diazotroph nitrogen loss to aggregation layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jaggloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_n_Lg","Large phyto nitrogen loss to aggregation layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jaggloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_n_Sm","Small phyto nitrogen loss to aggregation layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jaggloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    !  Register diagnostics for phytoplankton loss terms: viruses
    !

    vardesc_temp = vardesc("jvirloss_n_Di","Diazotroph nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_n_Lg","Large phyto nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_n_Sm","Small phyto nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for phytoplankton exudation
    !
    vardesc_temp = vardesc("jexuloss_n_Di","Diazotroph nitrogen loss via exudation",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jexuloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_n_Lg","Large phyto nitrogen loss via exudation",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jexuloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_n_Sm","Small phyto nitrogen loss via exudation",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jexuloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register dynamic silicate diagnostics
    !
    vardesc_temp = vardesc("nlg_diatoms","Fraction of large phytos that are diatoms",&
                           'h','L','s','dimensionless','f')
    cobalt%id_nlg_diatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_si_2_n_lg_diatoms","Si:N ratio in large diatoms",&
                           'h','L','s','mol Si mol N','f')
    cobalt%id_q_si_2_n_lg_diatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    !  Register diagnostics for phytoplankton loss terms: higher predators 
    !

!    vardesc_temp = vardesc("jhploss_n_Di","Diazotroph nitrogen loss to higher predators layer integral",&
!                           'h','L','s','mol N m-2 s-1','f')
!    phyto(DIAZO)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
!         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
!
!    vardesc_temp = vardesc("jhploss_n_Lg","Large phyto nitrogen loss to higher predators layer integral",&
!                           'h','L','s','mol N m-2 s-1','f')
!    phyto(LARGE)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
!         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
!
!    vardesc_temp = vardesc("jhploss_n_Sm_hl","High light Sm. phyto nitrogen loss to higher preds layer integral",&
!                           'h','L','s','mol N m-2 s-1','f')
!    phyto(SMALL_HL)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
!         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
!
!    vardesc_temp = vardesc("jhploss_n_Sm_ll","Low light Sm. phyto nitrogen loss to higher preds layer integral",&
!                           'h','L','s','mol N m-2 s-1','f')
!    phyto(SMALL_LL)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
!         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)


    !
    ! Register Phytoplankton Production Diagnostics
    !

    vardesc_temp = vardesc("juptake_n2_Di","Nitrogen fixation layer integral",'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_juptake_n2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_fe_Di","Diaz. phyto. Fed uptake layer integral",'h','L','s','mol Fe m-2 s-1','f')
    phyto(DIAZO)%id_juptake_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_fe_Lg","Large phyto. Fed uptake layer integral",'h','L','s','mol Fe m-2 s-1','f')
    phyto(LARGE)%id_juptake_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_fe_Sm","Small phyto. Fed uptake layer integral",&
                           'h','L','s','mol Fe m-2 s-1','f')
    phyto(SMALL)%id_juptake_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_nh4_Di","Diaz. phyto. NH4 uptake layer integral",'h','L','s','mol NH4 m-2 s-1','f')
    phyto(DIAZO)%id_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_nh4_Lg","Large phyto. NH4 uptake layer integral",'h','L','s','mol NH4 m-2 s-1','f')
    phyto(LARGE)%id_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_nh4_Sm","Small phyto. NH4 uptake layer integral",&
                           'h','L','s','mol NH4 m-2 s-1','f')
    phyto(SMALL)%id_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_no3_Di","Diaz. phyto. NO3 uptake layer integral",'h','L','s','mol NO3 m-2 s-1','f')
    phyto(DIAZO)%id_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_no3_Lg","Large phyto. NO3 uptake layer integral",'h','L','s','mol NO3 m-2 s-1','f')
    phyto(LARGE)%id_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_no3_Sm","Small phyto. NO3 uptake layer integral",&
                           'h','L','s','mol NO3 m-2 s-1','f')
    phyto(SMALL)%id_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_po4_Di","Diaz. phyto. PO4 uptake layer integral",'h','L','s','mol PO4 m-2 s-1','f')
    phyto(DIAZO)%id_juptake_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_po4_Lg","Large phyto. PO4 uptake layer integral",'h','L','s','mol PO4 m-2 s-1','f')
    phyto(LARGE)%id_juptake_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_po4_Sm","Small phyto. PO4 uptake layer integral",&
                           'h','L','s','mol PO4 m-2 s-1','f')
    phyto(SMALL)%id_juptake_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_sio4_Lg","Large phyto. SiO4 uptake layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_juptake_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi","Diazotroph Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmp","Small phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgp","Large phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register zooplankton diagnostics, starting with losses of zooplankton to ingestion by zooplankton
    !

    vardesc_temp = vardesc("jzloss_n_Smz","Small zooplankton nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Mdz","Medium-sized zooplankton nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Lgz","Large zooplankton nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for zooplankton loss terms: higher predators
    !

    vardesc_temp = vardesc("jhploss_n_Smz","Small zooplankton nitrogen loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_n_Mdz","Medium-sized zooplankton nitrogen loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_n_Lgz","Large zooplankton nitrogen loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register zooplankton ingestion rates
    !

    vardesc_temp = vardesc("jingest_n_Smz","Ingestion of nitrogen by small zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_Mdz","Ingestion of nitrogen by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_Lgz","Ingestion of nitrogen by large zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_p_Smz","Ingestion of phosphorous by small zooplankton, layer integral", &
                           'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_p_Mdz","Ingestion of phosphorous by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_p_Lgz","Ingestion of phosphorous by large zooplankton, layer integral",&
                           'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_sio2_Smz","Ingestion of sio2 by small zooplankton, layer integral",&
                           'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(1)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_sio2_Mdz","Ingestion of sio2 by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(2)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_sio2_Lgz","Ingestion of sio2 by large zooplankton, layer integral",&
                           'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(3)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_fe_Smz","Ingestion of Fe by small zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(1)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_fe_Mdz","Ingestion of Fe by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol Fe m-2 s-1','f')
    zoo(2)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_fe_Lgz","Ingestion of Fe by large zooplankton, layer integral",&
                           'h','L','s','mol Fe m-2 s-1','f')
    zoo(3)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register detrital production terms for zooplankton
    !

    vardesc_temp = vardesc("jprod_ndet_Smz","Production of nitrogen detritus by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_Mdz","Production of nitrogen detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_Lgz","Production of nitrogen detritus by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet_Smz","Production of phosphorous detritus by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet_Mdz","Production of phosphorous detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet_Lgz","Production of phosphorous detritus by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_Smz","Production of opal detritus by small zooplankton, layer integral",&
                   'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(1)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_Mdz","Production of opal detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(2)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_Lgz","Production of opal detritus by large zooplankton, layer integral",&
                   'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(3)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Smz","Production of sio4 through grazing/dissolution, layer integral",&
                   'h','L','s','mol SiO4 m-2 s-1','f')
    zoo(1)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Mdz","Production of sio4 through grazing/dissolution, layer integral",&
                   'h','L','s','mol SiO4 m-2 s-1','f')
    zoo(2)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Lgz","Production of sio4 through grazing/dissolution, layer integral",&
                   'h','L','s','mol SiO4 m-2 s-1','f')
    zoo(3)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet_Smz","Production of iron detritus by small zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(1)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet_Mdz","Production of iron detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(2)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet_Lgz","Production of iron detritus by large zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(3)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register dissolved organic/inorganic production terms for zooplankton
    !
    ! Labile dissolved organic nitrogen 
    vardesc_temp = vardesc("jprod_ldon_Smz","Production of labile dissolved organic nitrogen by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldon_Mdz","Production of labile dissolved organic nitrogen by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldon_Lgz","Production of labile dissolved organic nitrogen by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Labile dissolved organic phosphorous
    vardesc_temp = vardesc("jprod_ldop_Smz","Production of labile dissolved organic phosphorous by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldop_Mdz","Production of labile dissolved organic phosphorous by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldop_Lgz","Production of labile dissolved organic phosphorous by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Refractory dissolved organic nitrogen
    vardesc_temp = vardesc("jprod_srdon_Smz","Production of semi-refractory dissolved organic nitrogen by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdon_Mdz","Production of semi-refractory dissolved organic nitrogen by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdon_Lgz","Production of semi-refractory dissolved organic nitrogen by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Labile dissolved organic phosphorous
    vardesc_temp = vardesc("jprod_srdop_Smz","Production of semi-refractory dissolved organic phosphorous by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdop_Mdz","Production of semi-refractory dissolved organic phosphorous by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdop_Lgz","Production of semi-refractory dissolved organic phosphorous by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! semi-labile dissolved organic nitrogen 
    vardesc_temp = vardesc("jprod_sldon_Smz","Production of semi-labile dissolved organic nitrogen by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldon_Mdz","Production of semi-labile dissolved organic nitrogen by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldon_Lgz","Production of semi-labile dissolved organic nitrogen by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! semi-labile dissolved organic phosphorous
    vardesc_temp = vardesc("jprod_sldop_Smz","Production of semi-labile dissolved organic phosphorous by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldop_Mdz","Production of semi-labile dissolved organic phosphorous by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldop_Lgz","Production of semi-labile dissolved organic phosphorous by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! dissolved iron
    vardesc_temp = vardesc("jprod_fed_Smz","Production of dissolved iron by small zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(1)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fed_Mdz","Production of dissolved iron by medium-sized zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(2)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fed_Lgz","Production of dissolved iron by large zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(3)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! phosphate
    vardesc_temp = vardesc("jprod_po4_Smz","Production of phosphate by small zooplankton, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    zoo(1)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Mdz","Production of phosphate by medium-sized zooplankton, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    zoo(2)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Lgz","Production of phosphate by large zooplankton, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    zoo(3)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! ammonia
    vardesc_temp = vardesc("jprod_nh4_Smz","Production of ammonia by small zooplankton, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    zoo(1)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Mdz","Production of ammonia by medium-sized zooplankton, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    zoo(2)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Lgz","Production of ammonia by large zooplankton, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    zoo(3)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register zooplankton production terms 
    !

    vardesc_temp = vardesc("jprod_nsmz","Production of new biomass (nitrogen) by small zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nmdz","Production of new biomass (nitrogen) by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgz","Production of new biomass (nitrogen) by large zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Smz","Temperature limitation of small zooplankton",'h','L','s','dimensionless','f')
    zoo(1)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Mdz","Temperature limitation of medium-sized zooplankton",'h','L','s','dimensionless','f')
    zoo(2)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Lgz","Temperature limitation of large zooplankton",'h','L','s','dimensionless','f')
    zoo(3)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register bacterial diagnostics, starting with losses of bacteria to ingestion by zooplankton
    ! CAS: limit loss terms to N

    vardesc_temp = vardesc("jzloss_n_Bact","Bacterial nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for bacteria loss terms: higher predators
    !

!    vardesc_temp = vardesc("jhploss_n_Bact","Bacterial nitrogen loss to higher predators layer integral",&
!                           'h','L','s','mol N m-2 s-1','f')
!    bact(1)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
!         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for bacteria loss terms: viruses 
    !

    vardesc_temp = vardesc("jvirloss_n_Bact","Bacterial nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register bacterial uptake terms
    !

    vardesc_temp = vardesc("juptake_ldon","Bacterial uptake of labile dissolved organic nitrogen",'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_juptake_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
         
    vardesc_temp = vardesc("juptake_ldop","Bacterial uptake of labile dissolved organic phosphorous",'h','L','s','mol P m-2 s-1','f')
    bact(1)%id_juptake_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register dissolved inorganic production terms for bacteria 
    !
    ! phosphate
    vardesc_temp = vardesc("jprod_po4_Bact","Production of phosphate by bacteria, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    bact(1)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! ammonia
    vardesc_temp = vardesc("jprod_nh4_Bact","Production of ammonia by bacteria, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    bact(1)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register bacterial production terms
    !

    vardesc_temp = vardesc("jprod_nbact","Production of new biomass (nitrogen) by bacteria, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Bact","Temperature limitation of bacteria",'h','L','s','dimensionless','f')
    bact(1)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register general COBALT diagnostics
    !

    vardesc_temp = vardesc("co3_sol_arag","Carbonate Ion Solubility for Aragonite",'h','L','s','mol kg-1','f')
    cobalt%id_co3_sol_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("co3_sol_calc","Carbonate Ion Solubility for Calcite",'h','L','s','mol kg-1','f')
    cobalt%id_co3_sol_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("omega_arag","Carbonate Ion Saturation State for Aragonite",'h','L','s','mol kg-1','f')
    cobalt%id_omega_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("omega_calc","Carbonate Ion Saturation State for Calcite",'h','L','s','mol kg-1','f')
    cobalt%id_omega_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! A few overall production diagnostics
    !

    vardesc_temp = vardesc("jprod_cadet_arag","Aragonite CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_calc","Calcite CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_lithdet","Lithogenic detritus production layer integral",'h','L','s','g m-2 s-1','f')
    cobalt%id_jprod_lithdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    !vardesc_temp = vardesc("jprod_sidet","opal detritus production layer integral",'h','L','s','mol SiO2 m-2 s-1','f')
    !cobalt%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_sio4","sio4 production layer integral",'h','L','s','mol SiO2 m-2 s-1','f')
    !cobalt%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_fedet","Detrital Fedet production layer integral",'h','L','s','mol Fe m-2 s-1','f')
    !cobalt%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_ndet","Detrital PON production layer integral",'h','L','s','mol N m-2 s-1','f')
    !cobalt%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_pdet","Detrital phosphorus production layer integral",'h','L','s','mol P m-2 s-1','f')
    !cobalt%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_ldon","labile dissolved organic nitrogen production layer integral",&
    !                        'h','L','s','mol N m-2 s-1','f')
    !cobalt%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    ! 
    !vardesc_temp = vardesc("jprod_ldop","labile dissolved organic phosphorous production layer integral",&
    !                        'h','L','s','mol P m-2 s-1','f')
    !cobalt%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_srdon","refractory dissolved organic nitrogen production layer integral",&
    !                        'h','L','s','mol N m-2 s-1','f')
    !cobalt%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_srdop","refractory dissolved organic phosphorous production layer integral",&
    !                        'h','L','s','mol P m-2 s-1','f')
    !cobalt%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_sldon","semi-labile dissolved organic nitrogen production layer integral",&
    !                        'h','L','s','mol N m-2 s-1','f')
    !cobalt%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_sldop","semi-labile dissolved organic phosphorous production layer integral",&
    !                        'h','L','s','mol P m-2 s-1','f')
    !cobalt%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_fed","dissolved iron production layer integral",&
    !                        'h','L','s','mol Fe m-2 s-1','f')
    !cobalt%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_po4","phosphate production layer integral",&
    !                        'h','L','s','mol PO4 m-2 s-1','f')
    !cobalt%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("jprod_nh4","NH4 production layer integral",'h','L','s','mol NH4 m-2 s-1','f')
    !cobalt%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !     init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !
    ! loss diagnostics: detrital loss terms
    !
    !
    !vardesc_temp = vardesc("det_jzloss_n","nitrogen detritus loss to zooplankton layer integral",&
    !                       'h','L','s','mol N m-2 s-1','f')
    !cobalt%id_det_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !   init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !vardesc_temp = vardesc("det_jhploss_n","nitrogen detritus loss to higher predators layer integral",&
    !                       'h','L','s','mol N m-2 s-1','f')
    !cobalt%id_det_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
    !   init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    !
    ! Loss diagnostics: dissolution and remineralization
    !

    vardesc_temp = vardesc("jdiss_sidet","SiO2 detritus dissolution, layer integral",&
                           'h','L','s','mol m-2 s-1','f')
    cobalt%id_jdiss_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdiss_cadet_arag","CaCO3 detritus dissolution, layer integral", &
                           'h','L','s','mol CaCO3 m-2 s-1','f')
    cobalt%id_jdiss_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdiss_cadet_calc","CaCO3 detritus dissolution, layer integral", &
                           'h','L','s','mol CaCO3 m-2 s-1','f')
    cobalt%id_jdiss_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_ndet","Nitrogen detritus remineralization, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    cobalt%id_jremin_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_pdet","Phosphorous detritus remineralization, layer integral",&
                           'h','L','s','mol P m-2 s-1','f')
    cobalt%id_jremin_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_fedet","Iron detritus remineralization, layer integral",&
                           'h','L','s','mol m-2 s-1','f')
    cobalt%id_jremin_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! iron cycling diagnostics 
    !

    vardesc_temp = vardesc("jfe_ads","Iron adsorption layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jfe_ads = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfe_coast","Coastal iron efflux layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jfe_coast = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("kfe_eq_lig","Effective ligand binding strength",'h','L','s','kg mol-1','f')
    cobalt%id_kfe_eq_lig = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Temperature limitation diagnostics
    !

    vardesc_temp = vardesc("expkT","Eppley temperature limitation factor",'h','L','s','dimensionless','f')
    cobalt%id_expkT = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("hp_temp_lim","Temperature limitation of higher predators",'h','L','s','dimensionless','f')
    cobalt%id_hp_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Some additional light field diagnostics
    !

    vardesc_temp = vardesc("irr_inst","Instantaneous Light",'h','L','s','W m-2','f')
    cobalt%id_irr_inst = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irr_mix","Light averaged over mixing layer",'h','L','s','W m-2','f')
    cobalt%id_irr_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Nitrification/Denitrification diagnostics
    !

    vardesc_temp = vardesc("jnitrif","Nitrification layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jnitrif = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jno3denit_wc","Water column Denitrification layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jno3denit_wc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)


    !
    ! Some useful total layer integrals
    !

    vardesc_temp = vardesc("nphyto_tot","Total NO3: Di+Lg+Sm",'h','L','s','mol m-2 s-1','f')
    cobalt%id_nphyto_tot = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_c","Total Carbon (DIC+OC+IC) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_fe","Total Iron (Fed_OFe) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_n","Total Nitrogen (NO3+NH4+ON) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_p","Total Phosphorus (PO4+OP) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_si","Total Silicon (SiO4+SiO2) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_si = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("total_filter_feeding","Total filter feeding by large organisms",'h','L','s','mol N m-2 s-1','f')
    cobalt%id_total_filter_feeding = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("net_prim_prod","net primary production by all phytoplankton",'h','L','s','mol C m-2 yr-1','f')
    cobalt%id_net_prim_prod = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("gross_prim_prod","gross primary production by all phytoplankton",'h','L','s','mol C m-2 yr-1','f')
    cobalt%id_gross_prim_prod = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    !  Save river, depositon and bulk elemental fluxes
    !

    vardesc_temp = vardesc("dep_dry_fed","Dry Deposition of Iron to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_lith","Dry Deposition of Lithogenic Material",'h','1','s','g m-2 s-1','f')
    cobalt%id_dep_dry_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_nh4","Dry Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_no3","Dry Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_po4","Dry Deposition of Phosphate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_fed","Wet Deposition of Iron to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_lith","Wet Deposition of Lithogenic Material",'h','1','s','g m-2 s-1','f')
    cobalt%id_dep_wet_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_nh4","Wet Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_no3","Wet Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_po4","Wet Deposition of Phosphate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_alk","Alkalinity runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_dic","Dissolved Inorganic Carbon runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_fed","Iron runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_lith","Lithogenic runoff flux to the ocean",'h','1','s','g m-2 s-1','f')
    cobalt%id_runoff_flux_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_no3","Nitrate runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_ldon","LDON runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_sldon","SLDON runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_srdon","SRDON runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_ndet","NDET runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_po4","PO4 runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_ldop","LDOP runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_sldop","SLDOP runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_srdop","SRDOP runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 3D sinking information 
    !

    vardesc_temp = vardesc("fcadet_arag","CaCO3 sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc","CaCO3 sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet","fedet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet","lithdet sinking flux",'h','1','s','g m-2 s-1','f')
    cobalt%id_flithdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet","ndet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet","pdet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet","sidet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 2D sinking, bottom source/sink and burial diagnostics
    !

    vardesc_temp = vardesc("fcadet_arag_btm","CaCO3 sinking flux at bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_arag_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_btm","CaCO3 sinking flux at bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_calc_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_burial","CaCO3 permanent burial flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcased_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_input","CaCO3 flux into sediment layer",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcased_input = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_redis","CaCO3 redissolution from sediments",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcased_redis = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_btm","fedet sinking flux burial",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffedet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffe_sed","Sediment iron efflux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffe_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_btm","Lithogenic detrital sinking flux burial",'h','1','s','g m-2 s-1','f')
    cobalt%id_flithdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_btm","ndet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fnfeso4red_sed","Sediment Ndet Fe and SO4 reduction flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fnfeso4red_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fno3denit_sed","Sediment denitrification flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fno3denit_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fnoxic_sed","Sediment oxic Ndet remineralization flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fnoxic_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_btm","pdet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_btm","sidet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsidet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("frac_burial","fraction of organic matter buried",'h','1','s','dimensionless','f')
    cobalt%id_frac_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_burial","ndet burial flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_burial","pdet burial flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Surface Diagnostics
    !

    vardesc_temp = vardesc("pco2surf","Oceanic pCO2",'h','1','s','uatm','f')
    cobalt%id_pco2surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_alk","Surface Alkalinity",'h','1','s','eq kg-1','f')
    cobalt%id_sfc_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_cadet_arag","Surface Detrital Aragonite",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_cadet_calc","Surface Detrital Calcite",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_dic","Surface Dissolved Inorganic Carbon",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_fed","Surface Dissolved Iron",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ldon","Surface Labile Dissolved Organic Nitrogen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_sldon","Surface semi-labile Dissolved Organic Nitrogen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_srdon","Surface semi-refractory Dissolved Organic Nitrogen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3","Surface NO3",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4","Surface NH4",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4","Surface PO4",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_sio4","Surface SiO4",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_htotal","Surface Htotal",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_o2","Surface Oxygen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl","Surface Chl",'h','1','s','ug kg-1','f')
    cobalt%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irr","Surface Irradiance",'h','1','s','W m-2','f')
    cobalt%id_sfc_irr = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irr_mem","Surface Irradiance memory",'h','1','s','W m-2','f')
    cobalt%id_sfc_irr_mem = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_temp","Surface Temperature",'h','1','s','deg C','f')
    cobalt%id_sfc_temp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_temp","Bottom Temperature",'h','1','s','deg C','f')
    cobalt%id_btm_temp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_o2","Bottom Oxygen",'h','1','s','mol kg-1','f')
    cobalt%id_btm_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_ion","Surface Carbonate Ion",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_co3_ion = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_sol_arag","Surface Carbonate Ion Solubility for Aragonite",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_co3_sol_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_sol_calc","Surface Carbonate Ion Solubility for Calcite ",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_co3_sol_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nsmp","Surface small phyto. nitrogen",'h','1','s','mol kg-1','f')
    phyto(SMALL)%id_sfc_f_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nlgp","Surface large phyto. nitrogen",'h','1','s','mol kg-1','f')
    phyto(LARGE)%id_sfc_f_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ndi","Surface diazotroph nitrogen",'h','1','s','mol kg-1','f')
    phyto(DIAZO)%id_sfc_f_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_smp","Surface small phyto. chlorophyll",'h','1','s','ug kg-1','f')
    phyto(SMALL)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_lgp","Surface large phyto. chlorophyll",'h','1','s','ug kg-1','f')
    phyto(LARGE)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_di","Surface diazotroph chlorophyll",'h','1','s','mol kg-1','f')
    phyto(DIAZO)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_def_fe_smp","Surface small phyto. iron deficiency",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_def_fe_lgp","Surface large phyto. iron deficiency",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_def_fe_di","Surface diazotroph iron deficiency",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_felim_smp","Surface small phyto. iron uptake limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_felim_lgp","Surface large phyto. iron uptake limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_felim_di","Surface diazotroph iron uptake limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_q_fe_2_n_di","Surface diazotroph iron:nitrogen",'h','1','s','moles Fe (moles N)-1','f')
    phyto(DIAZO)%id_sfc_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_q_fe_2_n_smp","Surface small phyto. iron:nitrogen",'h','1','s','moles Fe (moles N)-1','f')
    phyto(SMALL)%id_sfc_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_q_fe_2_n_lgp","Surface large phyto. iron:nitrogen",'h','1','s','moles Fe (moles N)-1','f')
    phyto(LARGE)%id_sfc_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irrlim_smp","Surface small phyto. light limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irrlim_lgp","Surface large phyto. light limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irrlim_di","Surface diazotroph light limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_theta_smp","Surface small phyto. Chl:C",'h','1','s','g Chl (g C)-1','f')
    phyto(SMALL)%id_sfc_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_theta_lgp","Surface large phyto. Chl:C",'h','1','s','g Chl (g C)-1','f')
    phyto(LARGE)%id_sfc_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_theta_di","Surface diazotroph Chl:C",'h','1','s','g Chl (g C)-1','f')
    phyto(DIAZO)%id_sfc_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

   vardesc_temp = vardesc("sfc_mu_smp","Surface small phyto. Chl:C",'h','1','s','sec-1','f')
    phyto(SMALL)%id_sfc_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_mu_lgp","Surface large phyto. Chl:C",'h','1','s','sec-1','f')
    phyto(LARGE)%id_sfc_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_mu_di","Surface diazotroph growth rate",'h','1','s','sec-1','f')
    phyto(DIAZO)%id_sfc_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3lim_smp","Surface small phyto. nitrate limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3lim_lgp","Surface large phyto. nitrate limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4lim_smp","Surface small phyto. ammonia limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4lim_lgp","Surface large phyto. ammonia limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4lim_smp","Surface small phyto. phosphate limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4lim_lgp","Surface large phyto. phosphate limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4lim_di","Surface diazotroph phosphate limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 100m integrated fluxes
    !

    vardesc_temp = vardesc("jprod_allphytos_100","Total Nitrogen prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_allphytos_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi_100","Diazotroph nitrogen prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmp_100","Small phyto. nitrogen  prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgp_100","Large phyto. nitrogen  prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi_new_100","Diazotroph new (NO3-based) prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_new_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmp_new_100","Small phyto. new (NO3-based) prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n_new_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgp_new_100","Large phyto. new (NO3-based) prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n_new_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi_n2_100","Diazotroph nitrogen fixation in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_n2_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_ndi_100","Diazotroph nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nsmp_100","Small phyto. nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nlgp_100","Large phyto. nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_nsmp_100","Small phyto. nitrogen aggregation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jaggloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_nlgp_100","Large phyto. nitrogen aggregation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jaggloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_nsmp_100","Small phyto. nitrogen virus loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jvirloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_ndi_100","Diazotroph nitrogen exudation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jexuloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_nsmp_100","Small phyto. nitrogen exudation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jexuloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_nlgp_100","Large phyto. nitrogen exudation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jexuloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmz_100","Small zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nmdz_100","Medium zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgz_100","Large zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_nsmz_100","Small zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_nmdz_100","Medium zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_nlgz_100","Large zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nsmz_100","Small zooplankton nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nmdz_100","Medium zooplankton nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_nmdz_100","Medium zooplankton nitrogen loss to higher preds. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jhploss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_nlgz_100","Large zooplankton nitrogen loss to higher preds. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jhploss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_nmdz_100","Medium zooplankton nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_nlgz_100","Large zooplankton nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_don_nsmz_100","Small zooplankton dissolved org. nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jprod_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_don_nmdz_100","Medium zooplankton dissolved org. nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jprod_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nsmz_100","Small zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nmdz_100","Medium zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nlgz_100","Large zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_hp_100","Higher predator nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_hp_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_hp_100","Higher predator ingestion of nitrogen integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_hp_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_hp_100","Higher predator nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_hp_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nbact_100","Bacteria nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nbact_100","Bacteria nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_nbact_100","Bacteria nitrogen loss to viruses integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jvirloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nbact_100","Bacteria nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_ldon_nbact_100","Bacterial uptake of labile dissolved org. nitrogen in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_juptake_ldon_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_lithdet_100","Lithogenic detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_lithdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_100","Silica detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_sidet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_calc_100","Calcite detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_calc_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_arag_100","Aragonite detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_arag_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_ndet_100","Remineralization of nitrogen detritus integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jremin_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_mesozoo_200","Mesozooplankton Production, 200m integration",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_mesozoo_200 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 100m integrated biomass
    !

    vardesc_temp = vardesc("nsmp_100","Small phytoplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    phyto(SMALL)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nlgp_100","Large phytoplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    phyto(LARGE)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ndi_100","Diazotroph nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    phyto(DIAZO)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nsmz_100","Small zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(1)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nmdz_100","Medium zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(2)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nlgz_100","Large zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(3)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nbact_100","Bacterial nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    bact(1)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("silgp_100","Large phytoplankton silicon biomass in upper 100m",'h','1','s','mol m-2','f')
    cobalt%id_f_silg_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ndet_100","Nitrogen detritus biomass in upper 100m",'h','1','s','mol m-2','f')
    cobalt%id_f_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("don_100","Dissolved organic nitrogen (sr+sl+l) in upper 100m",'h','1','s','mol m-2','f')
    cobalt%id_f_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mesozoo_200","Mesozooplankton biomass, 200m integral",'h','1','s','mol m-2','f')
    cobalt%id_f_mesozoo_200 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    ! sinking flux = 100m
    !

    vardesc_temp = vardesc("fndet_100","Nitrogen detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_100","Phosphorous detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_100","Iron detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffedet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_100","Silicon detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsidet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_100","Calcite detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_calc_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_arag_100","Aragonite detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_arag_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_100","Lithogenic detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_flithdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Oxygen minima (value and location
    !

    vardesc_temp = vardesc("o2min","Minimum Oxygen",'h','1','s','mol kg-1','f')
    cobalt%id_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("z_o2min","Depth of Oxygen minimum",'h','1','s','m','f')
    cobalt%id_z_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Calcite and aragonite saturation depths
    !

    vardesc_temp = vardesc("z_sat_arag","Depth of Aragonite Saturation",'h','1','s','m','f')
    cobalt%id_z_sat_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         mask_variant=.TRUE.)

    vardesc_temp = vardesc("z_sat_calc","Depth of Calcite Saturation",'h','1','s','m','f')
    cobalt%id_z_sat_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         mask_variant=.TRUE.)



  end subroutine generic_COBALT_register_diag

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !===============================@===============================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_COBALT_params type
    !==============================================================    

    !=============
    !Block Starts: g_tracer_add_param
    !=============
    !Add the known experimental parameters used for calculations
    !in this module.
    !All the g_tracer_add_param calls must happen between 
    !g_tracer_start_param_list and g_tracer_end_param_list  calls.
    !This implementation enables runtime overwrite via field_table.

    call g_tracer_start_param_list(package_name)
    call g_tracer_add_param('init', cobalt%init, .false. )

    call g_tracer_add_param('htotal_scale_lo', cobalt%htotal_scale_lo, 0.01)
    call g_tracer_add_param('htotal_scale_hi', cobalt%htotal_scale_hi, 100.0)

    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    call g_tracer_add_param('RHO_0', cobalt%Rho_0, 1035.0)
    call g_tracer_add_param('NKML' , cobalt%nkml, 1)
    !-----------------------------------------------------------------------
    !       coefficients for O2 saturation
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_0', cobalt%a_0, 2.00907)
    call g_tracer_add_param('a_1', cobalt%a_1, 3.22014)
    call g_tracer_add_param('a_2', cobalt%a_2, 4.05010)
    call g_tracer_add_param('a_3', cobalt%a_3, 4.94457)
    call g_tracer_add_param('a_4', cobalt%a_4, -2.56847e-01)
    call g_tracer_add_param('a_5', cobalt%a_5, 3.88767)
    call g_tracer_add_param('b_0', cobalt%b_0, -6.24523e-03)
    call g_tracer_add_param('b_1', cobalt%b_1, -7.37614e-03)
    call g_tracer_add_param('b_2', cobalt%b_2, -1.03410e-02 )
    call g_tracer_add_param('b_3', cobalt%b_3, -8.17083e-03)
    call g_tracer_add_param('c_0', cobalt%c_0, -4.88682e-07)
    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !-----------------------------------------------------------------------
    !
    !  Compute the Schmidt number of CO2 in seawater using the 
    !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    !  7373-7382).
    !-----------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_co2', cobalt%a1_co2,  2068.9)
    call g_tracer_add_param('a2_co2', cobalt%a2_co2, -118.63)
    call g_tracer_add_param('a3_co2', cobalt%a3_co2,  2.9311)
    call g_tracer_add_param('a4_co2', cobalt%a4_co2, -0.027)
    !---------------------------------------------------------------------
    !  Compute the Schmidt number of O2 in seawater using the 
    !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
    !  Cycles, 12, 141-163).
    !---------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_o2', cobalt%a1_o2, 1929.7)
    call g_tracer_add_param('a2_o2', cobalt%a2_o2, -117.46)
    call g_tracer_add_param('a3_o2', cobalt%a3_o2, 3.116)
    call g_tracer_add_param('a4_o2', cobalt%a4_o2, -0.0306)
    !
    !-----------------------------------------------------------------------
    ! Stoichiometry
    !-----------------------------------------------------------------------
    !
    ! Values taken from OCMIP-II Biotic protocols after Anderson
    ! and Sarmiento (1994)
    !
    call g_tracer_add_param('mass_2_n', cobalt%mass_2_n, 106.0 / 16.0 * 12.0 * 1.87)         ! g mol N-1
    call g_tracer_add_param('n_2_n_denit', cobalt%n_2_n_denit, 472.0/(5.0*16.0))             ! mol N NO3 mol N org-1
    call g_tracer_add_param('o2_2_c', cobalt%o2_2_c, 150.0 / 106)                            ! mol O2 mol C-1
    call g_tracer_add_param('o2_2_nfix', cobalt%o2_2_nfix, (118.0+3.0/(5.0+3.0)*(150.0-118.0))/16.0) ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_nh4', cobalt%o2_2_nh4, 118.0 / 16)                         ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_nitrif', cobalt%o2_2_nitrif, 2.0)                          ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_no3', cobalt%o2_2_no3, 150.0 / 16.0)                       ! mol O2 mol N-1
    !
    !-----------------------------------------------------------------------
    ! Nutrient Limitation Parameters (phytoplankton) 
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('k_fed_Di', phyto(DIAZO)%k_fed, 5.0e-10)                   ! mol Fed kg-1
    call g_tracer_add_param('k_fed_Lg', phyto(LARGE)%k_fed, 5.0e-10)                   ! mol Fed kg-1
    call g_tracer_add_param('k_fed_Sm', phyto(SMALL)%k_fed,  1.0e-10)                 ! mol Fed kg-1
    call g_tracer_add_param('k_nh4_Lg', phyto(LARGE)%k_nh4,  5.0e-7)                  ! mol NH4 kg-1
    call g_tracer_add_param('k_nh4_Sm', phyto(SMALL)%k_nh4,  1.0e-7)                  ! mol NH4 kg-1
    call g_tracer_add_param('k_nh4_Di', phyto(DIAZO)%k_nh4,  5.0e-7)                  ! mol NH4 kg-1
    call g_tracer_add_param('k_no3_Lg', phyto(LARGE)%k_no3,  2.5e-6)                  ! mol NO3 kg-1
    call g_tracer_add_param('k_no3_Sm', phyto(SMALL)%k_no3,  5.0e-7)                  ! mol NO3 kg-1
    call g_tracer_add_param('k_no3_Di', phyto(DIAZO)%k_no3,  2.5e-6)                  ! mol NO3 kg-1
    call g_tracer_add_param('k_po4_Di', phyto(DIAZO)%k_po4,  5.0e-8)                  ! mol PO4 kg-1
    call g_tracer_add_param('k_po4_Lg', phyto(LARGE)%k_po4,  5.0e-8)                  ! mol PO4 kg-1
    call g_tracer_add_param('k_po4_Sm', phyto(SMALL)%k_po4,  1.0e-8)                  ! mol PO4 kg-1
    call g_tracer_add_param('k_sio4_Lg',phyto(LARGE)%k_sio4, 2.0e-6)                        ! mol SiO4 kg-1
    call g_tracer_add_param('k_fe_2_n_Di', phyto(DIAZO)%k_fe_2_n, 25.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    call g_tracer_add_param('k_fe_2_n_Lg', phyto(LARGE)%k_fe_2_n, 6.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    call g_tracer_add_param('k_fe_2_n_Sm',phyto(SMALL)%k_fe_2_n, 3.0e-6*106.0/16.0)        ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Sm',phyto(SMALL)%fe_2_n_max, 50.e-6*106.0/16.0)     ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Lg', phyto(LARGE)%fe_2_n_max, 500.0e-6*106.0/16.0)  ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Di', phyto(DIAZO)%fe_2_n_max, 500.0e-6*106.0/16.0)  ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_upt_fac', cobalt%fe_2_n_upt_fac, 15.0e-6)               ! mol Fe mol N-1
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton light limitation/growth rate
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('alpha_Di', phyto(DIAZO)%alpha,  1.0e-5 * 2.77e18 / 6.022e17) ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('alpha_Lg', phyto(LARGE)%alpha,  1.0e-5 * 2.77e18 / 6.022e17) ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('alpha_Sm', phyto(SMALL)%alpha,2.0e-5*2.77e18/6.022e17)  ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('kappa_eppley', cobalt%kappa_eppley, 0.063)                    ! deg C-1
    call g_tracer_add_param('P_C_max_Di', phyto(DIAZO)%P_C_max, 0.50/sperd)                ! s-1
    call g_tracer_add_param('P_C_max_Lg', phyto(LARGE)%P_C_max, 1.25/sperd)                ! s-1
    call g_tracer_add_param('P_C_max_Sm', phyto(SMALL)%P_C_max, 1.125/sperd)                 ! s-1
    call g_tracer_add_param('thetamax_Di', phyto(DIAZO)%thetamax, 0.03)                    ! g Chl g C-1
    call g_tracer_add_param('thetamax_Lg', phyto(LARGE)%thetamax, 0.05)                    ! g Chl g C-1
    call g_tracer_add_param('thetamax_Sm', phyto(SMALL)%thetamax, 0.03)                    ! g Chl g C-1
    call g_tracer_add_param('bresp_Di', phyto(DIAZO)%bresp,0.025/sperd)                    ! sec-1 
    call g_tracer_add_param('bresp_Lg', phyto(LARGE)%bresp,0.025/sperd)                    ! sec-1 
    call g_tracer_add_param('bresp_Sm', phyto(SMALL)%bresp,0.0225/sperd)                     ! sec-1 
    call g_tracer_add_param('thetamin', cobalt%thetamin, 0.002)                            ! g Chl g C-1
    call g_tracer_add_param('thetamin_nolim', cobalt%thetamin_nolim, 0.0)                  ! g Chl g C-1
    call g_tracer_add_param('zeta', cobalt%zeta, 0.05)                                     ! dimensionless
    call g_tracer_add_param('gamma_irr_mem', cobalt%gamma_irr_mem, 1.0 / sperd)            ! s-1
    !
    !-----------------------------------------------------------------------
    ! Nitrogen fixation inhibition parameters
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('k_n_inhib_Di', cobalt%k_n_inhib_Di, 1.0e-6)                    ! mol NO3 kg-1
    call g_tracer_add_param('o2_inhib_Di_pow', cobalt%o2_inhib_Di_pow, 4.0)                 ! mol O2-1 m3
    call g_tracer_add_param('o2_inhib_Di_sat', cobalt%o2_inhib_Di_sat, 3.0e-4)              ! mol O2 kg-1
    !
    !-----------------------------------------------------------------------
    ! Other stoichiometry
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('p_2_n_static', cobalt%p_2_n_static, .true. )
    call g_tracer_add_param('c_2_n', cobalt%c_2_n, 106.0 / 16.0)
    call g_tracer_add_param('alk_2_n_denit', cobalt%alk_2_n_denit, 552.0/472.0)             ! eq. alk mol NO3-1
    call g_tracer_add_param('p_2_n_static_Di', phyto(DIAZO)%p_2_n_static,1.0/40.0 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_static_Lg', phyto(LARGE)%p_2_n_static,1.0/16.0 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_static_Sm', phyto(SMALL)%p_2_n_static,1.0/16.0 )         ! mol P mol N-1
    call g_tracer_add_param('si_2_n_static_Lg', phyto(LARGE)%si_2_n_static, 2.0)            ! mol Si mol N-1
    call g_tracer_add_param('si_2_n_max_Lg', phyto(LARGE)%si_2_n_max, 5.0)                  ! mol Si mol N-1
    call g_tracer_add_param('ca_2_n_arag', cobalt%ca_2_n_arag, 0.020 * 106.0 / 16.0)        ! mol Ca mol N-1
    call g_tracer_add_param('ca_2_n_calc', cobalt%ca_2_n_calc, 0.010 * 106.0 / 16.0)        ! mol Ca mol N-1
    call g_tracer_add_param('caco3_sat_max', cobalt%caco3_sat_max,10.0)                     ! dimensionless
    !
    !-----------------------------------------------------------------------
    ! Zooplankton Stoichiometry - presently static
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('q_p_2_n_smz',zoo(1)%q_p_2_n, 1.0/16.0)          ! mol P mol N-1 
    call g_tracer_add_param('q_p_2_n_mdz',zoo(2)%q_p_2_n, 1.0/16.0)          ! mol P mol N-1 
    call g_tracer_add_param('q_p_2_n_lgz',zoo(3)%q_p_2_n, 1.0/16.0)          ! mol P mol N-1 
    !
    !-----------------------------------------------------------------------
    ! Bacteria Stoichiometry - presently static
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('q_p_2_n_bact',bact(1)%q_p_2_n, 1.0/16.0)        ! mol P mol N-1
    !
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton aggregation
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('agg_Sm',phyto(SMALL)%agg,0.1*1e6 / sperd)          ! s-1 (mole N kg)-1
    call g_tracer_add_param('agg_Di',phyto(DIAZO)%agg,  0    / sperd)            ! s-1 (mole N kg)-1
    call g_tracer_add_param('agg_Lg',phyto(LARGE)%agg,0.3*1e6/ sperd)            ! s-1 (mole N kg)-1
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton and bacterial losses to viruses
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('vir_Sm',phyto(SMALL)%vir, 0.025*1e6/sperd )  ! s-1 (mole N kg)-1
    call g_tracer_add_param('vir_Di',phyto(DIAZO)%vir, 0.0 )        ! s-1 (mole N kg)-1
    call g_tracer_add_param('vir_Lg',phyto(LARGE)%vir, 0.0 )        ! s-1 (mole N kg)-1
    call g_tracer_add_param('vir_Bact',bact(1)%vir,   0.033*1e6/sperd)   ! s-1 (mole N kg)-1
    call g_tracer_add_param('ktemp_vir',cobalt%vir_ktemp, 0.063)       ! C-1
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton losses to exudation
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('exu_Sm',phyto(SMALL)%exu, 0.13) ! dimensionless (fraction of NPP)
    call g_tracer_add_param('exu_Di',phyto(DIAZO)%exu, 0.13)       ! dimensionless (fraction of NPP) 
    call g_tracer_add_param('exu_Lg',phyto(LARGE)%exu, 0.13)       ! dimensionless (fraction of NPP) 
    !
    !-----------------------------------------------------------------------
    ! Zooplankton ingestion parameterization and temperature dependence
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('imax_smz',zoo(1)%imax, 1.42 / sperd)          ! s-1
    call g_tracer_add_param('imax_mdz',zoo(2)%imax, 0.57 / sperd)              ! s-1
    call g_tracer_add_param('imax_lgz',zoo(3)%imax, 0.23 / sperd)              ! s-1
    call g_tracer_add_param('ki_smz',zoo(1)%ki, 1.25e-6)                        ! moles N kg-1
    call g_tracer_add_param('ki_mdz',zoo(2)%ki, 1.25e-6)                        ! moles N kg-1
    call g_tracer_add_param('ki_lgz',zoo(3)%ki, 1.25e-6)                        ! moles N kg-1
    call g_tracer_add_param('ktemp_smz',zoo(1)%ktemp, 0.063)                   ! C-1
    call g_tracer_add_param('ktemp_mdz',zoo(2)%ktemp, 0.063)                   ! C-1
    call g_tracer_add_param('ktemp_lgz',zoo(3)%ktemp, 0.063)                   ! C-1
    !
    !-----------------------------------------------------------------------
    ! Bacterial growth and uptake parameters
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('mu_max_bact',bact(1)%mu_max, 1.0/sperd )          ! s-1 
    call g_tracer_add_param('k_ldon_bact', bact(1)%k_ldon,  5.0e-7)            ! mol ldon kg-1
    call g_tracer_add_param('ktemp_bact', bact(1)%ktemp, 0.063)                ! C-1
    !
    !-----------------------------------------------------------------------
    ! Zooplankton switching and prey preference parameters
    !-----------------------------------------------------------------------
    !
    ! parameters controlling the extent of biomass-based switching between
    ! multiple prey options 
    call g_tracer_add_param('nswitch_smz',zoo(1)%nswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('nswitch_mdz',zoo(2)%nswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('nswitch_lgz',zoo(3)%nswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('mswitch_smz',zoo(1)%mswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('mswitch_mdz',zoo(2)%mswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('mswitch_lgz',zoo(3)%mswitch, 2.0)          ! dimensionless
    ! innate prey availability for small zooplankton 
    call g_tracer_add_param('smz_ipa_smp',zoo(1)%ipa_smp, 1.0)    ! dimensionless
    call g_tracer_add_param('smz_ipa_lgp',zoo(1)%ipa_lgp, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_diaz',zoo(1)%ipa_diaz,0.0)         ! dimensionless
    call g_tracer_add_param('smz_ipa_smz',zoo(1)%ipa_smz, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_mdz',zoo(1)%ipa_mdz, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_lgz',zoo(1)%ipa_lgz, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_bact',zoo(1)%ipa_bact,0.25)         ! dimensionless
    call g_tracer_add_param('smz_ipa_det',zoo(1)%ipa_det, 0.0)          ! dimensionless
    ! innate prey availability for large zooplankton 
    call g_tracer_add_param('mdz_ipa_smp',zoo(2)%ipa_smp, 0.0)    ! dimensionless
    call g_tracer_add_param('mdz_ipa_lgp',zoo(2)%ipa_lgp, 1.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_diaz',zoo(2)%ipa_diaz,1.0)         ! dimensionless
    call g_tracer_add_param('mdz_ipa_smz',zoo(2)%ipa_smz, 1.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_mdz',zoo(2)%ipa_mdz, 0.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_lgz',zoo(2)%ipa_lgz, 0.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_bact',zoo(2)%ipa_bact, 0.0)        ! dimensionless
    call g_tracer_add_param('mdz_ipa_det',zoo(2)%ipa_det, 0.0)          ! dimensionless
    ! innate prey availability large predatory zooplankton/krill
    call g_tracer_add_param('lgz_ipa_smp',zoo(3)%ipa_smp, 0.0)   ! dimensionless
    call g_tracer_add_param('lgz_ipa_lgp',zoo(3)%ipa_lgp, 1.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_diaz',zoo(3)%ipa_diaz, 1.0)       ! dimensionless
    call g_tracer_add_param('lgz_ipa_smz',zoo(3)%ipa_smz, 0.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_mdz',zoo(3)%ipa_mdz, 1.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_lgz',zoo(3)%ipa_lgz, 0.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_bact',zoo(3)%ipa_bact, 0.0)       ! dimensionless
    call g_tracer_add_param('lgz_ipa_det',zoo(3)%ipa_det, 0.0)         ! dimensionless
    !
    !----------------------------------------------------------------------
    ! Zooplankton bioenergetics
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('gge_max_smz',zoo(1)%gge_max, 0.4)                   ! dimensionless
    call g_tracer_add_param('gge_max_mdz',zoo(2)%gge_max, 0.4)                   ! dimensionless
    call g_tracer_add_param('gge_max_lgz',zoo(3)%gge_max, 0.4)                   ! dimensionless
    call g_tracer_add_param('bresp_smz',zoo(1)%bresp, 0.020 / sperd)        ! s-1
    call g_tracer_add_param('bresp_mdz',zoo(2)%bresp, 0.008 / sperd)        ! s-1
    call g_tracer_add_param('bresp_lgz',zoo(3)%bresp, 0.0032 / sperd)       ! s-1
    !
    !----------------------------------------------------------------------
    ! Bacterial bioenergetics
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('gge_max_bact',bact(1)%gge_max,0.4)              ! dimensionless
    call g_tracer_add_param('bresp_bact',bact(1)%bresp, 0.0075/sperd)         ! s-1
    !
    !----------------------------------------------------------------------
    ! Partitioning of zooplankton ingestion to other compartments
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('phi_det_smz',zoo(1)%phi_det, 0.00)            ! dimensionless
    call g_tracer_add_param('phi_det_mdz',zoo(2)%phi_det, 0.20)            ! dimensionless
    call g_tracer_add_param('phi_det_lgz',zoo(3)%phi_det, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_ldon_smz',zoo(1)%phi_ldon, 0.55*0.30)     ! dimensionless
    call g_tracer_add_param('phi_ldon_mdz',zoo(2)%phi_ldon, 0.55*0.10)     ! dimensionless
    call g_tracer_add_param('phi_ldon_lgz',zoo(3)%phi_ldon, 0.55*0.0)      ! dimensionless
    call g_tracer_add_param('phi_ldop_smz',zoo(1)%phi_ldop, 0.45*0.30)     ! dimensionless
    call g_tracer_add_param('phi_ldop_mdz',zoo(2)%phi_ldop, 0.45*0.10)     ! dimensionless
    call g_tracer_add_param('phi_ldop_lgz',zoo(3)%phi_ldop, 0.45*0.0)      ! dimensionless
    call g_tracer_add_param('phi_srdon_smz',zoo(1)%phi_srdon, 0.05*0.30)   ! dimensionless
    call g_tracer_add_param('phi_srdon_mdz',zoo(2)%phi_srdon, 0.05*0.10)   ! dimensionless
    call g_tracer_add_param('phi_srdon_lgz',zoo(3)%phi_srdon, 0.05*0.0)    ! dimensionless
    call g_tracer_add_param('phi_srdop_smz',zoo(1)%phi_srdop, 0.15*0.30)   ! dimensionless
    call g_tracer_add_param('phi_srdop_mdz',zoo(2)%phi_srdop, 0.15*0.10)   ! dimensionless
    call g_tracer_add_param('phi_srdop_lgz',zoo(3)%phi_srdop, 0.15*0.0)    ! dimensionless
    call g_tracer_add_param('phi_sldon_smz',zoo(1)%phi_sldon, 0.4*0.30)    ! dimensionless
    call g_tracer_add_param('phi_sldon_mdz',zoo(2)%phi_sldon, 0.4*0.10)    ! dimensionless
    call g_tracer_add_param('phi_sldon_lgz',zoo(3)%phi_sldon, 0.4*0.0)     ! dimensionless
    call g_tracer_add_param('phi_sldop_smz',zoo(1)%phi_sldop, 0.4*0.30)    ! dimensionless
    call g_tracer_add_param('phi_sldop_mdz',zoo(2)%phi_sldop, 0.4*0.10)    ! dimensionless
    call g_tracer_add_param('phi_sldop_lgz',zoo(3)%phi_sldop, 0.4*0.0)     ! dimensionless
    call g_tracer_add_param('phi_nh4_smz',zoo(1)%phi_nh4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_nh4_mdz',zoo(2)%phi_nh4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_nh4_lgz',zoo(3)%phi_nh4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_po4_smz',zoo(1)%phi_po4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_po4_mdz',zoo(2)%phi_po4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_po4_lgz',zoo(3)%phi_po4, 0.30)            ! dimensionless
    !
    !----------------------------------------------------------------------
    ! Partitioning of viral losses to various dissolved pools
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('phi_ldon_vir',cobalt%lysis_phi_ldon, 0.55)    ! dimensionless
    call g_tracer_add_param('phi_srdon_vir',cobalt%lysis_phi_srdon, 0.05)  ! dimensionless
    call g_tracer_add_param('phi_sldon_vir',cobalt%lysis_phi_sldon, 0.40)  ! dimensionless
    call g_tracer_add_param('phi_ldop_vir',cobalt%lysis_phi_ldop, 0.45)    ! dimensionless
    call g_tracer_add_param('phi_srdop_vir',cobalt%lysis_phi_srdop, 0.15)  ! dimensionless
    call g_tracer_add_param('phi_sldop_vir',cobalt%lysis_phi_sldop, 0.40)  ! dimensionless
    ! 
    !----------------------------------------------------------------------
    ! Parameters for unresolved higher predators
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('imax_hp',     cobalt%imax_hp, 0.09/sperd)     ! s-1 
    call g_tracer_add_param('ki_hp',       cobalt%ki_hp, 1.2e-6)           ! mol N kg-1
    call g_tracer_add_param('coef_hp',     cobalt%coef_hp, 2.0)            ! dimensionless
    call g_tracer_add_param('ktemp_hp',    cobalt%ktemp_hp, 0.063)         ! C-1 
    call g_tracer_add_param('nswitch_hp',  cobalt%nswitch_hp, 2.0)         ! dimensionless
    call g_tracer_add_param('mswitch_hp',  cobalt%mswitch_hp, 2.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_smp',  cobalt%hp_ipa_smp, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_lgp',  cobalt%hp_ipa_lgp, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_diaz', cobalt%hp_ipa_diaz, 0.0)        ! dimensionless
    call g_tracer_add_param('hp_ipa_smz',  cobalt%hp_ipa_smz, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_mdz',  cobalt%hp_ipa_mdz, 1.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_lgz',  cobalt%hp_ipa_lgz, 1.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_bact', cobalt%hp_ipa_bact,0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_det',  cobalt%hp_ipa_det, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_phi_det',  cobalt%hp_phi_det, 0.35)        ! dimensionless
    call g_tracer_add_param('hp_phi_ldon', cobalt%hp_phi_ldon, 0.0)        ! dimensionless
    call g_tracer_add_param('hp_phi_ldop', cobalt%hp_phi_ldop, 0.0)        ! dimensionless
    call g_tracer_add_param('hp_phi_srdon', cobalt%hp_phi_srdon, 0.0)      ! dimensionless
    call g_tracer_add_param('hp_phi_srdop', cobalt%hp_phi_srdop, 0.0)      ! dimensionless
    call g_tracer_add_param('hp_phi_sldon', cobalt%hp_phi_sldon, 0.0)      ! dimensionless
    call g_tracer_add_param('hp_phi_sldop', cobalt%hp_phi_sldop, 0.0)      ! dimensionless
    call g_tracer_add_param('hp_phi_nh4',  cobalt%hp_phi_nh4, 0.65)         ! dimensionless
    call g_tracer_add_param('hp_phi_po4',  cobalt%hp_phi_po4, 0.65)         ! dimensionless
    !
    !----------------------------------------------------------------------
    ! Iron chemistry
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('felig_bkg', cobalt%felig_bkg, 1.0e-9)                           ! mol Fe kg-1
    call g_tracer_add_param('felig_2_don', cobalt%felig_2_don, 0.0e-3 / 40.0 * 106.0 / 16.0) ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_sed', cobalt%fe_2_n_sed, 100.0e-5 * 106 / 16)            ! mol Fe mol N-1
    call g_tracer_add_param('fe_coast', cobalt%fe_coast,1.0e-11 )                            ! mol Fe m kg-1 s-1
    call g_tracer_add_param('alpha_fescav',cobalt%alpha_fescav, 15.0/spery)                  ! sec-1 
    call g_tracer_add_param('beta_fescav',cobalt%beta_fescav, 0.0/spery)		     ! mol N-1 sec-1  
    call g_tracer_add_param('remin_eff_fedet',cobalt%remin_eff_fedet, 0.1)                   ! unitless 
    call g_tracer_add_param('ki_fescav',cobalt%ki_fescav, 1.0 )                              ! watts m-2
    call g_tracer_add_param('io_fescav',cobalt%io_fescav, 10.0 )                             ! watts m-2
    call g_tracer_add_param('gamma_fescav',cobalt%gamma_fescav, 1.0 )                        ! watts m-2
    call g_tracer_add_param('kfe_eq_lig_ll',cobalt%kfe_eq_lig_ll, 1.0e12)                    ! mol lig-1 kg
    call g_tracer_add_param('kfe_eq_lig_hl',cobalt%kfe_eq_lig_hl, 1e8)                       ! mol lig-1 kg
    !
    !-------------------------------------------------------------------------
    ! Remineralization
    !-------------------------------------------------------------------------
    !
    call g_tracer_add_param('k_o2', cobalt%k_o2, 20.0e-6)                                    ! mol O2 kg-1
    call g_tracer_add_param('o2_min', cobalt%o2_min, 1.0 * 1.0e-06)                          ! mol O2 kg-1
    call g_tracer_add_param('rpcaco3', cobalt%rpcaco3, 0.070/12.0*16.0/106.0*100.0)          ! mol N mol Ca-1
    call g_tracer_add_param('rplith',  cobalt%rplith,  0.065/12.0*16.0/106.0)                ! mol N g lith-1
    call g_tracer_add_param('rpsio2',  cobalt%rpsio2,  0.026/12.0*16.0/106.0*60.0)           ! mol N mol Si-1
    call g_tracer_add_param('gamma_ndet',  cobalt%gamma_ndet, cobalt%wsink / 188.0 )         ! s-1
    call g_tracer_add_param('gamma_cadet_arag',cobalt%gamma_cadet_arag,cobalt%wsink/760.0)   ! s-1
    call g_tracer_add_param('gamma_cadet_calc',cobalt%gamma_cadet_calc,cobalt%wsink/1343.0)  ! s-1
    call g_tracer_add_param('gamma_sidet',  cobalt%gamma_sidet, cobalt%wsink / 2000.0 )      ! s-1
    call g_tracer_add_param('phi_lith' ,  cobalt%phi_lith, 0.002)                            ! kg mol-1
    call g_tracer_add_param('k_lith',  cobalt%k_lith, 1e-6/sperd )                           ! s-1
    call g_tracer_add_param('z_sed',  cobalt%z_sed, 0.1 )                                    ! m
    call g_tracer_add_param('k_no3_denit',cobalt%k_no3_denit,1.0e-6)                        ! mol NO3 kg-1
    !
    !-----------------------------------------------------------------------
    ! Dissolved Organic Material
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('gamma_srdon',  cobalt%gamma_srdon, 1.0 / (18.0 * spery))          ! s-1
    call g_tracer_add_param('gamma_srdop',  cobalt%gamma_srdop, 1.0 / (4.0 * spery))           ! s-1
    call g_tracer_add_param('gamma_sldon',  cobalt%gamma_sldon, 1.0 / (90 * sperd))           ! s-1
    call g_tracer_add_param('gamma_sldop',  cobalt%gamma_sldop, 1.0 / (90 * sperd))           ! s-1
    !
    !-----------------------------------------------------------------------
    ! Nitrification
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('gamma_nitrif',  cobalt%gamma_nitrif, 1.0 / (30.0 * sperd))      ! s-1
    call g_tracer_add_param('irr_inhibit',  cobalt%irr_inhibit, 0.1)                         ! m2 W-1
    !
    !-----------------------------------------------------------------------
    ! Miscellaneous
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('tracer_debug',  cobalt%tracer_debug, .false.)

    call g_tracer_end_param_list(package_name)
    !===========
    !Block Ends: g_tracer_add_param
    !===========

  end subroutine user_add_params

  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list


    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'

    !
    !Add here only the parameters that are required at the time of registeration 
    !(to make flux exchanging Ocean tracers known for all PE's) 
    !
    call g_tracer_start_param_list(package_name)
    !
    call g_tracer_add_param('htotal_in', cobalt%htotal_in, 1.0e-08)
    !
    ! Sinking velocity of detritus: a value of 20 m d-1 is consistent with a characteristic sinking
    ! velocity of 100 m d-1 of marine aggregates and a disaggregation rate constant
    ! of 5 d-1 in the surface ocean (Clegg and Whitfield, 1992; Dunne, 1999).  Alternatively, 100 m d-1 
    ! is more in line with the deep water synthesis of Berelson (2002; Particel settling rates increase
    ! with depth in the ocean, DSR-II, 49, 237-252).
    !
    call g_tracer_add_param('wsink',  cobalt%wsink, 100.0 / sperd)                             ! m s-1

    call g_tracer_add_param('ice_restart_file'   , cobalt%ice_restart_file   , 'ice_cobalt.res.nc')
    call g_tracer_add_param('ocean_restart_file' , cobalt%ocean_restart_file , 'ocean_cobalt.res.nc' )
    call g_tracer_add_param('IC_file'       , cobalt%IC_file       , '')
    !
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file    = cobalt%ice_restart_file,&
         ocean_restart_file  = cobalt%ocean_restart_file )

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
    !
    !       ALK (Total carbonate alkalinity)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'alk',         &
         longname   = 'Alkalinity',  &
         units      = 'mol/kg',      &
         prog       = .true.,        &
         flux_runoff= .true.,        &
         flux_param = (/ 1.0e-03 /), &
         flux_bottom= .true.         )
    !
    !       Aragonite (Sinking detrital/particulate CaCO3)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cadet_arag',     &
         longname   = 'Detrital CaCO3', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         sink_rate  = cobalt%wsink,     &
         btm_reservoir = .true.         )
    !
    !       Calcite 
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cadet_calc',     &
         longname   = 'Detrital CaCO3', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         sink_rate  = cobalt%wsink,     &
         btm_reservoir = .true.         )
    !
    !       DIC (Dissolved inorganic carbon)
    !
    call g_tracer_add(tracer_list,package_name,                       &
         name       = 'dic',                                           &
         longname   = 'Dissolved Inorganic Carbon',                    &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'co2_flux',                                  &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMCO2,                                      &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                  &
         flux_gas_restart_file  = 'ocean_cobalt_airsea_flux.res.nc',    &
         flux_runoff= .true.,                                          &
         flux_param = (/12.011e-03  /),                                &
         flux_bottom= .true.,                                          &
         init_value = 0.001                                            )
    !
    !       Dissolved Fe (assumed to be all available to phytoplankton)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fed',            &
         longname   = 'Dissolved Iron', & 
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_runoff= .true.,           &
         flux_wetdep= .true.,           &
         flux_drydep= .true.,           &
         flux_param = (/ 55.847e-03 /), &
         flux_bottom= .true.            )
    !
    !    Fedet (Sinking detrital/particulate iron)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fedet',         &
         longname   = 'Detrital Iron', &
         units      = 'mol/kg',        &
         prog       = .true.,          &
         sink_rate  = cobalt%wsink,     &
         btm_reservoir = .true.        )
    !
    !       Diazotroph Fe (Iron in N2-fixing phytoplankton for variable Fe:N ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fedi',            &
         longname   = 'Diazotroph Iron', &
         units      = 'mol/kg',          &
         prog       = .true.             )
    !
    !       Large Fe (Iron in large phytoplankton to allow for variable Fe:N ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'felg',       &
         longname   = 'Large Phytoplankton Iron', &
         units      = 'mol/kg',     &
         prog       = .true.        )
    !
    !       Small Fe
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fesm',       &
         longname   = 'Small Phytoplankton Iron', &
         units      = 'mol/kg',     &
         prog       = .true.        )
    !
    !       LDON (Labile dissolved organic nitrogen)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ldon',           &
         flux_runoff= .true.,         &
         longname   = 'labile DON',     &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)       ) 
    !
    !       LDOP (Labile dissolved organic phosphorous)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ldop',           &
         flux_runoff= .true.,         &
         longname   = 'labile DOP',     &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)       ) 
    !
    !       LITH (Lithogenic aluminosilicate particles)
    !
    call g_tracer_add(tracer_list,package_name,     &
         name       = 'lith',                       &
         longname   = 'Lithogenic Aluminosilicate', &
         units      = 'g/kg',                       &
         prog       = .true.,                       &
         flux_runoff= .true.,                       &
         flux_wetdep= .true.,                       &
         flux_drydep= .true.,                       &
         flux_param = (/ 1.0e-03 /)                 )
    !
    !     LITHdet (Detrital Lithogenic aluminosilicate particles)  
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'lithdet',   &
         longname   = 'lithdet',   &
         units      = 'g/kg',      &
         prog       = .true.,      &
         sink_rate  = cobalt%wsink, &
         btm_reservoir = .true.    )
    !
    !       NBact: Bacteria nitrogen
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nbact',          &
         longname   = 'bacterial',      &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !    Ndet (Sinking detrital/particulate Nitrogen)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndet',      &
         longname   = 'ndet',      &
         flux_runoff= .true.,      &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         sink_rate  = cobalt%wsink,&
         btm_reservoir = .true.,   &
         flux_param = (/ 1.0e-3 /) ) 
    !
    !    NDi (assumed to be facultative N2-fixers, with a variable N:P ratio
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndi',                 &
         longname   = 'Diazotroph Nitrogen', &
         units      = 'mol/kg',              &
         prog       = .true.                 )

    !
    !    NLg (assumed to be a dynamic combination of diatoms and other 
    !         eukaryotes all effectively greater than 5 um in diameter,
    !         and having a fixed C:N ratio)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nlg',            &
         longname   = 'Large Phytoplankton Nitrogen', &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !       NSm (Nitrogen in picoplankton and nanoplankton
    !            ~less than 5 um in diameter and having a fixed C:N:P ratio)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nsm',            &
         longname   = 'Small Phytoplankton Nitrogen', &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !       NH4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nh4',             &
         longname   = 'Ammonia',         &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_wetdep= .true.,            &
         flux_drydep= .true.,            &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )
    !
    !       NO3
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'no3',             &
         longname   = 'Nitrate',         &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_runoff= .true.,            &
         flux_wetdep= .true.,            &
         flux_drydep= .true.,            &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )
    !
    !       O2
    !
    call g_tracer_add(tracer_list,package_name,                        &
         name       = 'o2',                                            &
         longname   = 'Oxygen',                                        &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'o2_flux',                                   &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMO2,                                       &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                  &
         flux_gas_restart_file  = 'ocean_cobalt_airsea_flux.res.nc',    &
         flux_bottom= .true.             )
    !
    !    Pdet (Sinking detrital/particulate Phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,         &
         name       = 'pdet',                           &
         longname   = 'Detrital Phosphorus',            &
         units      = 'mol/kg',                         &
         prog       = .true.,                           &
         sink_rate  = cobalt%wsink,                      &
         btm_reservoir = .true.    )
    !
    !       PO4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'po4',       &
         longname   = 'Phosphate', &
         flux_runoff= .true.,      &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         flux_wetdep= .true.,      &
         flux_drydep= .true.,      &
         flux_bottom= .true.,      &
         flux_param = (/ 1.0e-3 /) )
    !
    !       SRDON (Semi-Refractory dissolved organic nitrogen)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'srdon',           &
         longname   = 'Semi-Refractory DON', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)  )
    !
    !       SRDOP (Semi-Refractory dissolved organic phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'srdop',           &
         longname   = 'Semi-Refractory DOP', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)      )
    !
    !       SLDON (Semilabile dissolved organic nitrogen)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sldon',           &
         longname   = 'Semilabile DON', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)      )
    !
    !       SLDOP (Semilabile dissolved organic phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sldop',           &
         longname   = 'Semilabile DOP', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)      )
    !
    !       Sidet (Sinking detrital/particulate Silicon)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sidet',            &
         longname   = 'Detrital Silicon', &
         units      = 'mol/kg',           &
         prog       = .true.,             &
         sink_rate  = cobalt%wsink,        &
         btm_reservoir = .true.    )
    !
    !    SiLg (Silicon in large phytoplankton for variable Si:N ratios
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'silg',          &
         longname   = 'Large Phytoplankton Silicon', &
         units      = 'mol/kg',        &
         prog       = .true.           )
    !
    !       SiO4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sio4',     &
         longname   = 'Silicate', &
         units      = 'mol/kg',   &
         prog       = .true.,     &
         flux_bottom= .true.      )

    !
    !     Small zooplankton N  
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nsmz',     &
         longname   = 'Small Zooplankton Nitrogen', &
         units      = 'mol/kg',   &
         prog       = .true.     ) 

    !
    !     Medium-sized zooplankton N
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nmdz',     &
         longname   = 'Medium-sized zooplankton Nitrogen', &
         units      = 'mol/kg',   &
         prog       = .true.     ) 

    !
    !     Large zooplankton N (Pred zoo + krill)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nlgz',     &
         longname   = 'large Zooplankton Nitrogen', &
         units      = 'mol/kg',   &
         prog       = .true.     ) 

    !===========================================================
    !Diagnostic Tracers
    !===========================================================
    !
    !    Cased (CaCO3 in sediments)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cased',          &
         longname   = 'Sediment CaCO3', &
         units      = 'mol m-3',       &
         prog       = .false.           )
    !
    !       Chl (Chlorophyll)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'chl',         &
         longname   = 'Chlorophyll', &
         units      = 'ug kg-1',     &
         prog       = .false.,       &
         init_value = 0.08           )
    !
    !       CO3_ion (Carbonate ion)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'co3_ion',       &
         longname   = 'Carbonate ion', &
         units      = 'mol/kg',        &
         prog       = .false.          )
    !
    !      cadet_arag_btf (Aragonite flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,  &
         name       = 'cadet_arag_btf',          &
         longname   = 'aragonite flux to Sediments', &
         units      = 'mol m-2 s-1',             &
         prog       = .false.                    )
    !
    !      cadet_calc_btf (Calcite flux to sediments) 
    !
    call g_tracer_add(tracer_list,package_name,  &
         name       = 'cadet_calc_btf',          &
         longname   = 'calcite flux to Sediments', &
         units      = 'mol m-2 s-1',             &
         prog       = .false.                    )
    !
    !      lithdet_btf (Lithogenic flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'lithdet_btf',            &
         longname   = 'Lith flux to Sediments', &
         units      = 'g m-2 s-1',              &
         prog       = .false.                   )
    !
    !      ndet_btf (N flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndet_btf',            &
         longname   = 'N flux to Sediments', &
         units      = 'mol m-2 s-1',         &
         prog       = .false.                )
    !
    !      pdet_btf (P flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'pdet_btf',            &
         longname   = 'P flux to Sediments', &
         units      = 'mol m-2 s-1',         &
         prog       = .false.                )
    !
    !      sidet_btf (SiO2 flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'sidet_btf',              &
         longname   = 'SiO2 flux to Sediments', &
         units      = 'mol m-2 s-1',            &
         prog       = .false.                   )
    !
    !       htotal (H+ ion concentration)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'htotal',               &
         longname   = 'H+ ion concentration', &
         units      = 'mol/kg',               &
         prog       = .false.,                &
         init_value = cobalt%htotal_in         )
    !
    !       Irr_mem (Irradiance Memory)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'irr_mem',           &
         longname   = 'Irradiance memory', &
         units      = 'Watts/m^2',         &
         prog       = .false.              )


  end subroutine user_add_tracers


  ! <SUBROUTINE NAME="generic_COBALT_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_COBALT_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_update_from_copler'

    real, dimension(:,:)  ,pointer    :: stf_alk,dry_no3,wet_no3

    !
    ! NO3 has deposition, river flux, and negative deposition contribution to alkalinity
    !
    call g_tracer_get_pointer(tracer_list,'no3','drydep',dry_no3)
    call g_tracer_get_pointer(tracer_list,'no3','wetdep',wet_no3)

    call g_tracer_get_pointer(tracer_list,'alk','stf',stf_alk)

    stf_alk = stf_alk - dry_no3 - wet_no3 ! update 'tracer%stf' thru pointer

    return
  end subroutine generic_COBALT_update_from_coupler

  ! <SUBROUTINE NAME="generic_COBALT_update_from_bottom">
  ! 
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !
  !   This routine calculates bottom fluxes for tracers with bottom reservoirs.
  !   It is called near the end of the time step, meaning that the fluxes 
  !   calculated pertain to the next time step.  
  ! 
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_bottom(tracer_list,dt, tau, model_time) 
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment 
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  !
  ! </SUBROUTINE>
  subroutine generic_COBALT_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    !
    ! The bottom reservoirs of aragonite and calcite are immediately redistributed to the
    ! water column as a bottom flux (btf) where they impact the alkalinity and DIC
    !
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',cobalt%fcadet_arag_btm,isd,jsd)
    cobalt%fcadet_arag_btm = cobalt%fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (cobalt%id_fcadet_arag_btm .gt. 0)           &
         used = send_data(cobalt%id_fcadet_arag_btm,cobalt%fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',cobalt%fcadet_calc_btm,isd,jsd)
    cobalt%fcadet_calc_btm = cobalt%fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (cobalt%id_fcadet_calc_btm .gt. 0)           &
         used = send_data(cobalt%id_fcadet_calc_btm, cobalt%fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',cobalt%ffedet_btm,isd,jsd)
    cobalt%ffedet_btm = cobalt%ffedet_btm/dt
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (cobalt%id_ffedet_btm .gt. 0)           &
         used = send_data(cobalt%id_ffedet_btm, cobalt%ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',cobalt%flithdet_btm,isd,jsd)
    cobalt%flithdet_btm = cobalt%flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (cobalt%id_flithdet_btm .gt. 0)           &
         used = send_data(cobalt%id_flithdet_btm, cobalt%flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of 
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',cobalt%fndet_btm,isd,jsd)
    cobalt%fndet_btm = cobalt%fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (cobalt%id_fndet_btm .gt. 0)           &
         used = send_data(cobalt%id_fndet_btm,cobalt%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',cobalt%fpdet_btm,isd,jsd)
    cobalt%fpdet_btm = cobalt%fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (cobalt%id_fpdet_btm .gt. 0)           &
         used = send_data(cobalt%id_fpdet_btm,cobalt%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',cobalt%fsidet_btm,isd,jsd)
    cobalt%fsidet_btm = cobalt%fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (cobalt%id_fsidet_btm .gt. 0)           &
         used = send_data(cobalt%id_fsidet_btm,    cobalt%fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine generic_COBALT_update_from_bottom

  ! <SUBROUTINE NAME="generic_COBALT_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This is the subroutine to contain most of the biogeochemistry for calculating the 
  !   interaction of tracers with each other and with outside forcings.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_source(tracer_list,Temp,Salt,dzt,hblt_depth,&
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
  subroutine generic_COBALT_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,hblt_depth,&
       ilb,jlb,tau,dt,grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band)

    type(g_tracer_type),            pointer    :: tracer_list
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

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_update_from_source'
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau, i, j, k , kblt, m, n, k_100, k_200 
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast,grid_kmt
    !
    !------------------------------------------------------------------------
    ! Local Variables
    !------------------------------------------------------------------------
    !
    logical :: used, first
    integer :: nb
    real :: r_dt
    real :: feprime
    real :: juptake_di_tot2nterm
    real :: log_btm_flx
    real :: P_C_m
    real :: p_lim_nhet
    real :: TK, PRESS, PKSPA, PKSPC
    real :: tmp_hblt, tmp_irrad, tmp_irrad_ML,tmp_opacity
    real :: drho_dzt
    real, dimension(:), Allocatable   :: tmp_irr_band
    real, dimension(:,:), Allocatable :: rho_dzt_100, rho_dzt_200
    real,dimension(1:NUM_ZOO,1:NUM_PREY) :: ipa_matrix,pa_matrix,ingest_matrix
    real,dimension(1:NUM_PREY) :: hp_ipa_vec,hp_pa_vec,hp_ingest_vec
    real,dimension(1:NUM_PREY) :: prey_vec,prey_p2n_vec,prey_fe2n_vec,prey_si2n_vec
    real,dimension(1:NUM_ZOO)  :: tot_prey
    real :: tot_prey_hp, sw_fac_denom, ingest_p2n, refuge_conc 
    real :: bact_ldon_lim, bact_uptake_ratio, vmax_bact
    real :: fpoc_btm, log_fpoc_btm

    r_dt = 1.0 / dt

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

    call mpp_clock_begin(id_clock_carbon_calculations)
    !Get necessary fields
    call g_tracer_get_values(tracer_list,'htotal','field', cobalt%f_htotal,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'po4'   ,'field', cobalt%f_po4,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'sio4'  ,'field', cobalt%f_sio4,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'alk'   ,'field', cobalt%f_alk,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dic'   ,'field', cobalt%f_dic  ,isd,jsd,ntau=tau)

    !---------------------------------------------------------------------
    !Calculate co3_ion
    !Also calculate co2 fluxes csurf and alpha for the next round of exchnage
    !---------------------------------------------------------------------
   
    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       cobalt%htotallo(i,j) = cobalt%htotal_scale_lo * cobalt%f_htotal(i,j,k)
       cobalt%htotalhi(i,j) = cobalt%htotal_scale_hi * cobalt%f_htotal(i,j,k)
    enddo; enddo ; !} i, j
 
    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         cobalt%f_dic(:,:,k),                          &
         cobalt%f_po4(:,:,k),                          &  
         cobalt%f_sio4(:,:,k),                         &
         cobalt%f_alk(:,:,k),                          &
         cobalt%htotallo, cobalt%htotalhi,&
                                !InOut
         cobalt%f_htotal(:,:,k),                       & 
                                !OUT
         co2star=cobalt%co2_csurf(:,:), alpha=cobalt%co2_alpha(:,:), &
         pCO2surf=cobalt%pco2_csurf(:,:), &
         co3_ion=cobalt%f_co3_ion(:,:,k))

    do k = 2, nk
       do j = jsc, jec ; do i = isc, iec  !{
          cobalt%htotallo(i,j) = cobalt%htotal_scale_lo * cobalt%f_htotal(i,j,k)
          cobalt%htotalhi(i,j) = cobalt%htotal_scale_hi * cobalt%f_htotal(i,j,k)
       enddo; enddo ; !} i, j
  
       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
            Temp(:,:,k), Salt(:,:,k),                    &
            cobalt%f_dic(:,:,k),                          &
            cobalt%f_po4(:,:,k),                          &  
            cobalt%f_sio4(:,:,k),                         &
            cobalt%f_alk(:,:,k),                          &
            cobalt%htotallo, cobalt%htotalhi,&
                                !InOut
            cobalt%f_htotal(:,:,k),                       & 
                                !OUT
            co3_ion=cobalt%f_co3_ion(:,:,k))
    enddo

    call g_tracer_set_values(tracer_list,'htotal','field',cobalt%f_htotal  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion','field',cobalt%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'dic','alpha',cobalt%co2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','csurf',cobalt%co2_csurf    ,isd,jsd)

    call mpp_clock_end(id_clock_carbon_calculations)


    !---------------------------------------------------------------------
    ! Get positive tracer concentrations
    !---------------------------------------------------------------------

    call mpp_clock_begin(id_clock_phyto_growth)

    call g_tracer_get_values(tracer_list,'cadet_arag','field',cobalt%f_cadet_arag ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'cadet_calc','field',cobalt%f_cadet_calc ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fed'    ,'field',cobalt%f_fed      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fedet'  ,'field',cobalt%f_fedet    ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ldon'   ,'field',cobalt%f_ldon     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ldop'   ,'field',cobalt%f_ldop     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'lith'   ,'field',cobalt%f_lith     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'lithdet','field',cobalt%f_lithdet  ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ndet'   ,'field',cobalt%f_ndet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nh4'    ,'field',cobalt%f_nh4      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'no3'    ,'field',cobalt%f_no3      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'o2'     ,'field',cobalt%f_o2       ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'pdet'   ,'field',cobalt%f_pdet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'po4'    ,'field',cobalt%f_po4      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'srdon'   ,'field',cobalt%f_srdon   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'srdop'   ,'field',cobalt%f_srdop   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sldon'   ,'field',cobalt%f_sldon   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sldop'   ,'field',cobalt%f_sldop   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sidet'  ,'field',cobalt%f_sidet    ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sio4'   ,'field',cobalt%f_sio4     ,isd,jsd,ntau=tau,positive=.true.)
    !
    ! phytoplankton fields
    !
    call g_tracer_get_values(tracer_list,'fedi'   ,'field',phyto(DIAZO)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'felg'   ,'field',phyto(LARGE)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fesm','field',phyto(SMALL)%f_fe(:,:,:),isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ndi'    ,'field',phyto(DIAZO)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nlg'    ,'field',phyto(LARGE)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nsm' ,'field',phyto(SMALL)%f_n(:,:,:),isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'silg'   ,'field',cobalt%f_silg     ,isd,jsd,ntau=tau,positive=.true.)
    !
    ! zooplankton fields
    !
    call g_tracer_get_values(tracer_list,'nsmz'    ,'field',zoo(1)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nmdz'    ,'field',zoo(2)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nlgz'    ,'field',zoo(3)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    !
    ! bacteria
    !
    call g_tracer_get_values(tracer_list,'nbact'   ,'field',bact(1)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    !
    ! diagnostic tracers that are passed between time steps (except chlorophyll)
    !
    call g_tracer_get_values(tracer_list,'cased'  ,'field',cobalt%f_cased    ,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'co3_ion','field',cobalt%f_co3_ion  ,isd,jsd,ntau=1,positive=.true.)
    call g_tracer_get_values(tracer_list,'cadet_arag_btf','field',cobalt%f_cadet_arag_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'cadet_calc_btf','field',cobalt%f_cadet_calc_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'lithdet_btf','field',cobalt%f_lithdet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'ndet_btf','field',cobalt%f_ndet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'pdet_btf','field',cobalt%f_pdet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'sidet_btf','field',cobalt%f_sidet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'irr_mem','field',cobalt%f_irr_mem ,isd,jsd,ntau=1)

    cobalt%zt = 0.0
    cobalt%zm = 0.0
    ! minimum concentration below which predation/basal respiration stops
    refuge_conc = 1.0e-9


!
!-----------------------------------------------------------------------------------
! 1: Phytoplankton growth and nutrient uptake calculations
!-----------------------------------------------------------------------------------
!
    !
    !-----------------------------------------------------------------------------------
    ! 1.1: Nutrient Limitation Calculations
    !-----------------------------------------------------------------------------------
    !
    ! Calculate iron cell quota
    !
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec
       do n = 1,NUM_PHYTO    !{
          phyto(n)%q_fe_2_n(i,j,k) = max(0.0, phyto(n)%f_fe(i,j,k)/ &
                 max(epsln,phyto(n)%f_n(i,j,k)))
          phyto(n)%q_p_2_n(i,j,k) = phyto(n)%p_2_n_static
       enddo  !} n
       !
       ! N limitation with NH4 inhibition after Frost and Franzen (1992)
       !
       do n= 2, NUM_PHYTO   !{
          phyto(n)%no3lim(i,j,k) = cobalt%f_no3(i,j,k) / &
             ( (phyto(n)%k_no3+cobalt%f_no3(i,j,k)) * (1.0 + cobalt%f_nh4(i,j,k)/phyto(n)%k_nh4) )
          phyto(n)%nh4lim(i,j,k) = cobalt%f_nh4(i,j,k) / (phyto(n)%k_nh4 + cobalt%f_nh4(i,j,k))
       enddo !} n
       !
       ! O2 inhibition term for diazotrophs
       !
       n = DIAZO
       phyto(n)%o2lim(i,j,k) = (1.0 - cobalt%f_o2(i,j,k)**cobalt%o2_inhib_Di_pow / &
         (cobalt%f_o2(i,j,k)**cobalt%o2_inhib_Di_pow+cobalt%o2_inhib_Di_sat**cobalt%o2_inhib_Di_pow))
       !
       ! SiO4, PO4 and Fe uptake limitation with Michaelis-Mentin 
       !
       phyto(LARGE)%silim(i,j,k) = cobalt%f_sio4(i,j,k) / (phyto(LARGE)%k_sio4 + cobalt%f_sio4(i,j,k))
       do n= 1, NUM_PHYTO   !{
          phyto(n)%po4lim(i,j,k) = cobalt%f_po4(i,j,k) / (phyto(n)%k_po4 + cobalt%f_po4(i,j,k))
          phyto(n)%felim(i,j,k)  = cobalt%f_fed(i,j,k) / (phyto(n)%k_fed + cobalt%f_fed(i,j,k))
          phyto(n)%def_fe(i,j,k) = phyto(n)%q_fe_2_n(i,j,k)**2.0 / (phyto(n)%k_fe_2_n**2.0 +  &
               phyto(n)%q_fe_2_n(i,j,k)**2.0)
       enddo !} n
    enddo;  enddo ;  enddo !} i,j,k
    !
    ! Calculate nutrient limitation based on the most limiting nutrient (liebig_lim)
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
       n=DIAZO
       phyto(n)%liebig_lim(i,j,k) = phyto(n)%o2lim(i,j,k)* &
          min(phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
       do n= 2, NUM_PHYTO   !{
          phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k),&
             phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
       enddo !} n
    enddo;  enddo ;  enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    ! 1.2: Light Limitation/Growth Calculations
    !-----------------------------------------------------------------------
    !
    ! Create relevant light fields based on incident radiation and opacity
    ! information passed from the ocean code
    !
    allocate(tmp_irr_band(nbands))
    do j = jsc, jec ; do i = isc, iec   !{

       do nb=1,nbands !{
          if (max_wavelength_band(nb) .lt. 710.0) then !{
             tmp_irr_band(nb) = sw_pen_band(nb,i,j)
          else
             tmp_irr_band(nb) = 0.0
          endif !}
       enddo !}

       kblt = 0 ; tmp_irrad_ML = 0.0 ; tmp_hblt = 0.0
       do k = 1, nk !{
          tmp_irrad = 0.0
          do nb=1,nbands !{
             tmp_opacity = opacity_band(nb,i,j,k)
             tmp_irrad = tmp_irrad + max(0.0,tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k) * 0.5))
             ! Change tmp_irr_band from being the value atop layer k to the value
             ! at the bottom of layer k.
             tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k))
          enddo !}
          cobalt%irr_inst(i,j,k) = tmp_irrad * grid_tmask(i,j,k)
          cobalt%irr_mix(i,j,k) = tmp_irrad * grid_tmask(i,j,k)
          if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then !{
             kblt = kblt+1
             tmp_irrad_ML = tmp_irrad_ML + cobalt%irr_mix(i,j,k) * dzt(i,j,k)
             tmp_hblt = tmp_hblt + dzt(i,j,k)
          endif !}
       enddo !} k-loop
       cobalt%irr_mix(i,j,1:kblt) = tmp_irrad_ML / max(1.0e-6,tmp_hblt)
    enddo;  enddo !} i,j

    deallocate(tmp_irr_band)
    !
    ! Calculate the temperature limitation (expkT) and the time integrated
    ! irradiance (f_irr_mem) to which the Chl:C ratio responds (~24 hours)
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
       cobalt%expkT(i,j,k) = exp(cobalt%kappa_eppley * Temp(i,j,k))
       cobalt%f_irr_mem(i,j,k) = (cobalt%f_irr_mem(i,j,k) + (cobalt%irr_mix(i,j,k) - &
          cobalt%f_irr_mem(i,j,k)) * min(1.0,cobalt%gamma_irr_mem * dt)) * grid_tmask(i,j,k)
    enddo; enddo ; enddo !} i,j,k
    !
    ! Phytoplankton growth rate calculation based on Geider et al. (1997)
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%gross_prim_prod(i,j,k) = 0.0
       cobalt%f_chl(i,j,k) = 0.0

       do n = 1, NUM_PHYTO   !{
          P_C_m = phyto(n)%liebig_lim(i,j,k)*phyto(n)%P_C_max*cobalt%expkT(i,j,k)+epsln
          phyto(n)%theta(i,j,k) = (phyto(n)%thetamax-cobalt%thetamin) / (1.0 +                   &
             phyto(n)%thetamax*phyto(n)%alpha*cobalt%f_irr_mem(i,j,k)*0.5 /  &
             P_C_m) + cobalt%thetamin
          cobalt%f_chl(i,j,k) = cobalt%f_chl(i,j,k)+cobalt%c_2_n*12.0e6*phyto(n)%theta(i,j,k)*   &
             phyto(n)%f_n(i,j,k)
          phyto(n)%irrlim(i,j,k) = (1.0-exp(-phyto(n)%alpha*cobalt%irr_inst(i,j,k)*              &
             phyto(n)%theta(i,j,k)/P_C_m))

          ! calculate the growth rate
          phyto(n)%mu(i,j,k) = P_C_m / (1.0 + cobalt%zeta) * phyto(n)%irrlim(i,j,k) - &
             cobalt%expkT(i,j,k)*phyto(n)%bresp*                                      &
             phyto(n)%f_n(i,j,k)/(refuge_conc + phyto(n)%f_n(i,j,k))

          cobalt%gross_prim_prod(i,j,k) = cobalt%gross_prim_prod(i,j,k) + P_C_m*phyto(n)%irrlim(i,j,k)* &
                                          phyto(n)%f_n(i,j,k)
          ! Negative growth assumed to go to cell death rather than respiration (see manual) 
          cobalt%net_prim_prod(i,j,k) = max(phyto(n)%mu(i,j,k),0.0)*phyto(n)%f_n(i,j,k)
       enddo !} n

       cobalt%gross_prim_prod(i,j,k) = cobalt%gross_prim_prod(i,j,k)*cobalt%c_2_n*spery
       cobalt%net_prim_prod(i,j,k) = cobalt%net_prim_prod(i,j,k)*cobalt%c_2_n*spery
    enddo;  enddo ; enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    ! 1.3: Nutrient uptake calculations 
    !-----------------------------------------------------------------------
    !
    ! Uptake of nitrate and ammonia
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       n = DIAZO
       !juptake_di_tot2nterm=max(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)* &
       !  (1.0/(cobalt%f_no3(i,j,k)+cobalt%f_nh4(i,j,k)+cobalt%k_n_inhib_Di)))
       phyto(n)%juptake_n2(i,j,k) =  max(0.0,(1.0 - phyto(LARGE)%no3lim(i,j,k) - phyto(LARGE)%nh4lim(i,j,k))* &
          phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
       phyto(n)%juptake_nh4(i,j,k) = max(0.0,phyto(LARGE)%nh4lim(i,j,k)* phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
       phyto(n)%juptake_no3(i,j,k) = max(0.0,phyto(LARGE)%no3lim(i,j,k)* phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)) 
       ! If growth is negative, net remineralization of organic material
       phyto(n)%juptake_nh4(i,j,k) = phyto(n)%juptake_nh4(i,j,k) + &
                                     min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
       phyto(n)%jprod_n(i,j,k) = phyto(n)%juptake_nh4(i,j,k) + phyto(n)%juptake_no3(i,j,k) + &
          phyto(n)%juptake_n2(i,j,k)
       do n = 2, NUM_PHYTO !{
          phyto(n)%juptake_no3(i,j,k) = max( 0.0, phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*   & 
             phyto(n)%no3lim(i,j,k)/(phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k)+epsln) )
          phyto(n)%juptake_nh4(i,j,k) = max( 0.0, phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*   & 
             phyto(n)%nh4lim(i,j,k)/(phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k)+epsln) )
          ! If growth is negative, net remineralization of organic material
          phyto(n)%juptake_nh4(i,j,k) = phyto(n)%juptake_nh4(i,j,k) + &
                                        min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
          phyto(n)%jprod_n(i,j,k) = phyto(n)%juptake_nh4(i,j,k) + phyto(n)%juptake_no3(i,j,k)
       enddo !} n
    enddo;  enddo ; enddo !} i,j,k
    !
    ! Phosphorous uptake
    ! 
    do k = 1, nk  ;    do j = jsc, jec ;      do i = isc, iec   !{
       n=DIAZO
       phyto(n)%juptake_po4(i,j,k) = (phyto(n)%juptake_n2(i,j,k)+phyto(n)%juptake_nh4(i,j,k) + &
          phyto(n)%juptake_no3(i,j,k))*phyto(n)%p_2_n_static
       do n = 2, NUM_PHYTO
          phyto(n)%juptake_po4(i,j,k) = (phyto(n)%juptake_no3(i,j,k)+   &
                  phyto(n)%juptake_nh4(i,j,k)) * phyto(n)%p_2_n_static
       enddo !} n
    enddo; enddo ; enddo !} i,j,k
    !
    ! Iron uptake
    ! 
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       do n = 1, NUM_PHYTO  !{
          if (phyto(n)%q_fe_2_n(i,j,k).lt.phyto(n)%fe_2_n_max) then
             phyto(n)%juptake_fe(i,j,k) = phyto(n)%P_C_max*cobalt%expkT(i,j,k)*phyto(n)%f_n(i,j,k)* &
                phyto(n)%felim(i,j,k)*cobalt%fe_2_n_upt_fac
          else 
             phyto(n)%juptake_fe(i,j,k) = 0.0
          endif
       enddo   !} n
    enddo; enddo ; enddo !} i,j,k
    !
    ! Silicate uptake
    !
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%nlg_diatoms(i,j,k)=phyto(LARGE)%f_n(i,j,k)*phyto(LARGE)%silim(i,j,k)
       cobalt%q_si_2_n_lg_diatoms(i,j,k)= cobalt%f_silg(i,j,k)/ &
             (cobalt%nlg_diatoms(i,j,k) + epsln)
       phyto(LARGE)%juptake_sio4(i,j,k) = &
             max(phyto(LARGE)%juptake_no3(i,j,k)+phyto(LARGE)%juptake_nh4(i,j,k),0.0)*phyto(LARGE)%silim(i,j,k)* &
             phyto(LARGE)%silim(i,j,k)*phyto(LARGE)%si_2_n_max 

       ! CAS: set q_si_2_n values for each of the phyto groups for consumption calculations
       ! Note that this is si_2_n in large phytoplankton pool, not in diatoms themselves 
       phyto(LARGE)%q_si_2_n(i,j,k) = cobalt%f_silg(i,j,k)/(phyto(LARGE)%f_n(i,j,k)+epsln)

    enddo; enddo ; enddo !} i,j,k
    call mpp_clock_end(id_clock_phyto_growth)
!
!-----------------------------------------------------------------------
! 2: Bacterial Growth and Uptake Calculations 
!-----------------------------------------------------------------------
!
    !
    ! calculate an effective maximum ldon uptake rate (at 0 deg. C) for bacteria
    ! from specified values of bact(1)%gge_max, bact(1)%mu_max and bact(1)%bresp
    !

    call mpp_clock_begin(id_clock_bacteria_growth)
    vmax_bact = (1/bact(1)%gge_max)*(bact(1)%mu_max + bact(1)%bresp)
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec   !{
       bact(1)%temp_lim(i,j,k) = exp(bact(1)%ktemp*Temp(i,j,k))
       bact_ldon_lim = cobalt%f_ldon(i,j,k)/(bact(1)%k_ldon + cobalt%f_ldon(i,j,k))
       bact(1)%juptake_ldon(i,j,k) = vmax_bact*bact(1)%temp_lim(i,j,k)*bact_ldon_lim* &
          bact(1)%f_n(i,j,k)
       bact_uptake_ratio = ( cobalt%f_ldop(i,j,k)/max(cobalt%f_ldon(i,j,k),epsln) )
       bact(1)%juptake_ldop(i,j,k) = bact(1)%juptake_ldon(i,j,k)*bact_uptake_ratio
       if (bact_uptake_ratio.lt.bact(1)%q_p_2_n) then
          bact(1)%jprod_n(i,j,k) = bact(1)%gge_max*bact(1)%juptake_ldop(i,j,k)*16.0 - &
             bact(1)%f_n(i,j,k)/(refuge_conc + bact(1)%f_n(i,j,k)) *                     &
             bact(1)%temp_lim(i,j,k)*bact(1)%bresp*bact(1)%f_n(i,j,k)
       else
          bact(1)%jprod_n(i,j,k) = bact(1)%gge_max*bact(1)%juptake_ldon(i,j,k) - &
             bact(1)%f_n(i,j,k)/(refuge_conc + bact(1)%f_n(i,j,k)) *                &
             bact(1)%temp_lim(i,j,k)*bact(1)%bresp*bact(1)%f_n(i,j,k)
       endif
    enddo; enddo ; enddo !} i,j,k
    call mpp_clock_end(id_clock_bacteria_growth)
!
!-----------------------------------------------------------------------
! 3: Plankton foodweb dynamics
!-----------------------------------------------------------------------
!
    !
    ! 3.1 Plankton foodweb dynamics: consumption by zooplankton and higher predators
    !

    call mpp_clock_begin(id_clock_zooplankton_calculations)

    !
    ! Set-up local matrices for calculating zooplankton ingestion of
    ! multiple prey types.  The rows are consumers (i.e., NUM_ZOO zooplankton
    ! groups), the columns are food sources (i.e., NUM_PREY potential food sources)
    !
    ! ipa_matrix = the innate prey availability matrix
    ! pa_matrix = prey availability matrix after accounting for switching 
    ! ingest_matrix = ingestion matrix
    ! tot_prey = total prey available to predator m 
    !
    ! The definition of predator-prey matrices is intended to allow for
    ! efficient experimentation with predator-prey interconnections.
    ! However, we are still working to reduce the runtime required to
    ! include this feature.  The matrix structures are thus included,
    ! but the standard COBALT interactions have been hard-coded such
    ! that changing linkages requires changing the prey availability
    ! values and adding additional code to handle the new linkages.
    !
    ! With regard to stoichiometry, the primary ingestion calculations
    ! (i.e., those within the i, j, k loops) are coded to allow for 
    ! variable stoichiometry.  Several sections of the code corresponding
    ! to predator-prey and other linkages not in included in the
    ! default COBALT parameterizations have been commented out to
    ! avoid unnecessary calculations.
    !

    do m = 1,NUM_ZOO !{
       ipa_matrix(m,1) = zoo(m)%ipa_diaz
       ipa_matrix(m,2) = zoo(m)%ipa_lgp
       ipa_matrix(m,3) = zoo(m)%ipa_smp
       ipa_matrix(m,4) = zoo(m)%ipa_bact
       ipa_matrix(m,5) = zoo(m)%ipa_smz
       ipa_matrix(m,6) = zoo(m)%ipa_mdz
       ipa_matrix(m,7) = zoo(m)%ipa_lgz
       ipa_matrix(m,8) = zoo(m)%ipa_det
       tot_prey(m) = 0.0
       do n = 1,NUM_PREY !{
           pa_matrix(m,n) = 0.0
           ingest_matrix(m,n) = 0.0
       enddo !} n
    enddo !} m

    !
    ! Set-up local matrices for calculating higher predator ingestion
    ! of multiple prey types
    !

    hp_ipa_vec(1) = cobalt%hp_ipa_diaz
    hp_ipa_vec(2) = cobalt%hp_ipa_lgp
    hp_ipa_vec(3) = cobalt%hp_ipa_smp
    hp_ipa_vec(4) = cobalt%hp_ipa_bact
    hp_ipa_vec(5) = cobalt%hp_ipa_smz
    hp_ipa_vec(6) = cobalt%hp_ipa_mdz
    hp_ipa_vec(7) = cobalt%hp_ipa_lgz
    hp_ipa_vec(8) = cobalt%hp_ipa_det
    tot_prey_hp = 0.0
    do n = 1,NUM_PREY  !{  
       hp_pa_vec(n) = 0.0                  
       hp_ingest_vec(n) = 0.0              
    enddo !} n

    ! 
    ! Set all static stoichiometric ratios outside k,j,i loop
    !

    prey_p2n_vec(1) = phyto(DIAZO)%p_2_n_static
    prey_p2n_vec(2) = phyto(LARGE)%p_2_n_static
    prey_p2n_vec(3) = phyto(SMALL)%p_2_n_static
    prey_p2n_vec(4) = bact(1)%q_p_2_n
    prey_p2n_vec(5) = zoo(1)%q_p_2_n
    prey_p2n_vec(6) = zoo(2)%q_p_2_n
    prey_p2n_vec(7) = zoo(3)%q_p_2_n

    prey_fe2n_vec(4) = 0.0
    prey_fe2n_vec(5) = 0.0
    prey_fe2n_vec(6) = 0.0
    prey_fe2n_vec(7) = 0.0

    prey_si2n_vec(1) = 0.0
    prey_si2n_vec(3) = 0.0
    prey_si2n_vec(4) = 0.0
    prey_si2n_vec(5) = 0.0
    prey_si2n_vec(6) = 0.0
    prey_si2n_vec(7) = 0.0

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; !{

       !
       ! 3.1.1: Calculate zooplankton ingestion fluxes
       !

       ! Prey vectors for ingestion and loss calculations 
       ! (note: ordering of phytoplankton must be consistent with
       !  DIAZO, LARGE, SMALL ordering inherited from TOPAZ)
       !
       prey_vec(1) = max(phyto(DIAZO)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(2) = max(phyto(LARGE)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(3) = max(phyto(SMALL)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(4) = max(bact(1)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(5) = max(zoo(1)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(6) = max(zoo(2)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(7) = max(zoo(3)%f_n(i,j,k) - refuge_conc,0.0)
       prey_vec(8) = max(cobalt%f_ndet(i,j,k) - refuge_conc,0.0)
       ! 
       ! Set dynamic stoichiometric rations inside k,j,i loop
       prey_p2n_vec(8) = cobalt%f_pdet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)
       prey_fe2n_vec(1) = phyto(DIAZO)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(2) = phyto(LARGE)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(3) = phyto(SMALL)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(8) = cobalt%f_fedet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)
       prey_si2n_vec(2) = phyto(LARGE)%q_si_2_n(i,j,k)
       prey_si2n_vec(8) = cobalt%f_sidet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)

       !
       ! Calculate zooplankton ingestion
       !
       ! Small zooplankton (m = 1) consuming small phytoplankton (3) and
       ! bacteria (4).  sw_fac_denom is the denominator of the abundance-
       ! based switching factor, tot_prey is the total available prey 
       ! after accounting for switching.
       !

       m = 1 
       zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*Temp(i,j,k)) 
       sw_fac_denom = (ipa_matrix(m,3)*prey_vec(3))**zoo(m)%nswitch + &
                      (ipa_matrix(m,4)*prey_vec(4))**zoo(m)%nswitch
       pa_matrix(m,3) = ipa_matrix(m,3)* &
                        ( (ipa_matrix(m,3)*prey_vec(3))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,4) = ipa_matrix(m,4)* &
                        ( (ipa_matrix(m,4)*prey_vec(4))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) = pa_matrix(m,3)*prey_vec(3) + pa_matrix(m,4)*prey_vec(4)
       ingest_matrix(m,3) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,3)* &
                            prey_vec(3)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,4) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,4)* &
                            prey_vec(4)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,3) + ingest_matrix(m,4)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,3)*prey_p2n_vec(3) + &
                                 ingest_matrix(m,4)*prey_p2n_vec(4)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,3)*prey_fe2n_vec(3)

       !
       ! Medium zooplankton (m = 2) consuming diazotrophs (1), large
       ! phytoplankton (2), and small zooplankton (5) 
       !

       m = 2 
       zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*Temp(i,j,k))
       sw_fac_denom = (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch + &
                      (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch + &
                      (ipa_matrix(m,5)*prey_vec(5))**zoo(m)%nswitch
       pa_matrix(m,1) = ipa_matrix(m,1)* &
                        ( (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,2) = ipa_matrix(m,2)* & 
                        ( (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,5) = ipa_matrix(m,5)* & 
                        ( (ipa_matrix(m,5)*prey_vec(5))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) = pa_matrix(m,1)*prey_vec(1) + pa_matrix(m,2)*prey_vec(2) + &
                     pa_matrix(m,5)*prey_vec(5)
       ingest_matrix(m,1) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,1)* &
                            prey_vec(1)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,2) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,2)* &
                            prey_vec(2)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,5) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,5)* &
                            prey_vec(5)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,1) + ingest_matrix(m,2) + &
                                 ingest_matrix(m,5)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,1)*prey_p2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_p2n_vec(2) + &
                                 ingest_matrix(m,5)*prey_p2n_vec(5)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,1)*prey_fe2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_fe2n_vec(2)
       zoo(m)%jingest_sio2(i,j,k) = ingest_matrix(m,2)*prey_si2n_vec(2)

       !
       ! Large zooplankton (m = 3) consuming diazotrophs (2), large phytoplankton (2)
       ! and medium zooplankton (6)
       !

       m = 3
       zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*Temp(i,j,k))
       sw_fac_denom = (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch + &
                      (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch + &
                      (ipa_matrix(m,6)*prey_vec(6))**zoo(m)%nswitch
       pa_matrix(m,1) = ipa_matrix(m,1)* &
                        ( (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,2) = ipa_matrix(m,2)* &
                        ( (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,6) = ipa_matrix(m,6)* &
                        ( (ipa_matrix(m,6)*prey_vec(6))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) = pa_matrix(m,1)*prey_vec(1) + pa_matrix(m,2)*prey_vec(2) + &
                     pa_matrix(m,6)*prey_vec(6)
       ingest_matrix(m,1) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,1)* &
                            prey_vec(1)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,2) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,2)* &
                            prey_vec(2)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,6) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*pa_matrix(m,6)* &
                            prey_vec(6)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,1) + ingest_matrix(m,2) + &
                                 ingest_matrix(m,6)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,1)*prey_p2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_p2n_vec(2) + &
                                 ingest_matrix(m,6)*prey_p2n_vec(6)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,1)*prey_fe2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_fe2n_vec(2)
       zoo(m)%jingest_sio2(i,j,k) = ingest_matrix(m,2)*prey_si2n_vec(2)

       cobalt%total_filter_feeding(i,j,k) = ingest_matrix(2,1) + ingest_matrix(2,2) + &
          ingest_matrix(2,3) + ingest_matrix(3,1) + ingest_matrix(3,2) + & 
          ingest_matrix(3,3) + hp_ingest_vec(1) + hp_ingest_vec(2) + hp_ingest_vec(3) 

       !
       ! Calculate losses to zooplankton
       !

       do n = 1,NUM_PHYTO
          phyto(n)%jzloss_n(i,j,k) = 0.0
       enddo

       do m = 1,NUM_ZOO !{
          phyto(DIAZO)%jzloss_n(i,j,k) = phyto(DIAZO)%jzloss_n(i,j,k) + ingest_matrix(m,DIAZO)
          phyto(LARGE)%jzloss_n(i,j,k) = phyto(LARGE)%jzloss_n(i,j,k) + ingest_matrix(m,LARGE)
          phyto(SMALL)%jzloss_n(i,j,k) = phyto(SMALL)%jzloss_n(i,j,k) + ingest_matrix(m,SMALL)
       enddo !} m

       do n = 1,NUM_PHYTO !{
          phyto(n)%jzloss_p(i,j,k) = phyto(n)%jzloss_n(i,j,k)*prey_p2n_vec(n)
          phyto(n)%jzloss_fe(i,j,k) = phyto(n)%jzloss_n(i,j,k)*prey_fe2n_vec(n)
          phyto(n)%jzloss_sio2(i,j,k) = phyto(n)%jzloss_n(i,j,k)*prey_si2n_vec(n)  
       enddo !} n

       !
       ! losses of bacteria to zooplankton 
       !

       bact(1)%jzloss_n(i,j,k) = 0.0
       do m = 1,NUM_ZOO !{
          bact(1)%jzloss_n(i,j,k) = bact(1)%jzloss_n(i,j,k) + ingest_matrix(m,4)
       enddo !} m
       bact(1)%jzloss_p(i,j,k) = bact(1)%jzloss_n(i,j,k)*prey_p2n_vec(4)

       !
       ! losses of zooplankton to zooplankton
       !

       do n = 1,NUM_ZOO !{
          zoo(n)%jzloss_n(i,j,k) = 0.0

          do m = 1,NUM_ZOO !{
             zoo(n)%jzloss_n(i,j,k) = zoo(n)%jzloss_n(i,j,k) + ingest_matrix(m,NUM_PHYTO+1+n)
          enddo !} m

          zoo(n)%jzloss_p(i,j,k) = zoo(n)%jzloss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+1+n)
       enddo !} n

       !
       ! losses of detritus to zooplankton (no detrivory in default settings) 
       !
       !cobalt%det_jzloss_n(i,j,k) = 0.0
       !
       !do m = 1,NUM_ZOO !{
       !   cobalt%det_jzloss_n(i,j,k) = cobalt%det_jzloss_n(i,j,k)+ingest_matrix(m,NUM_PHYTO+NUM_ZOO+2)
       !enddo !} m
       !
       !cobalt%det_jzloss_p(i,j,k) = cobalt%det_jzloss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+NUM_ZOO+2)
       !cobalt%det_jzloss_fe(i,j,k) = cobalt%det_jzloss_n(i,j,k)*prey_fe2n_vec(NUM_PHYTO+NUM_ZOO+2)
       !cobalt%det_jzloss_si(i,j,k) = cobalt%det_jzloss_si(i,j,k)*prey_si2n_vec(NUM_PHYTO+NUM_ZOO+2)

       !
       ! 3.1.2 Calculate ingestion by higher predators
       !

       ! The higher-predator ingestion calculations mirror those used for zooplankton
       !
       cobalt%hp_temp_lim(i,j,k) = exp(cobalt%ktemp_hp*Temp(i,j,k))
       sw_fac_denom = (hp_ipa_vec(6)*prey_vec(6))**cobalt%nswitch_hp + &
                      (hp_ipa_vec(7)*prey_vec(7))**cobalt%nswitch_hp
       hp_pa_vec(6) = hp_ipa_vec(6)* &
                      ( (hp_ipa_vec(6)*prey_vec(6))**cobalt%nswitch_hp / &
                        (sw_fac_denom+epsln) )**(1.0/cobalt%mswitch_hp)
       hp_pa_vec(7) = hp_ipa_vec(7)* &
                      ( (hp_ipa_vec(7)*prey_vec(7))**cobalt%nswitch_hp / &
                        (sw_fac_denom+epsln) )**(1.0/cobalt%mswitch_hp)
       tot_prey_hp = hp_pa_vec(6)*prey_vec(6) + hp_pa_vec(7)*prey_vec(7)
       hp_ingest_vec(6) = cobalt%hp_temp_lim(i,j,k)*cobalt%imax_hp*hp_pa_vec(6)* &
                            prey_vec(6)*tot_prey_hp**(cobalt%coef_hp-1)/ &
                            (cobalt%ki_hp+tot_prey_hp)
       hp_ingest_vec(7) = cobalt%hp_temp_lim(i,j,k)*cobalt%imax_hp*hp_pa_vec(7)* &
                            prey_vec(7)*tot_prey_hp**(cobalt%coef_hp-1)/ &
                            (cobalt%ki_hp+tot_prey_hp)
       cobalt%hp_jingest_n(i,j,k) = hp_ingest_vec(6) + hp_ingest_vec(7)
       cobalt%hp_jingest_p(i,j,k) = hp_ingest_vec(6)*prey_p2n_vec(6) + &
                                    hp_ingest_vec(7)*prey_p2n_vec(7)
       !
       ! No iron and sio2 ingestion by higher predators with default settings
       !
       !cobalt%hp_jingest_fe(i,j,k) = hp_ingest_vec(6)*prey_fe2n_vec(6) + &
       !                              hp_ingest_vec(7)*prey_fe2n_vec(7)
       !cobalt%hp_jingest_sio2(i,j,k) = hp_ingest_vec(6)*prey_si2n_vec(6) + &
       !                                hp_ingest_vec(7)*prey_si2n_vec(7)

       !
       ! Calculate losses to higher predators
       !

       ! losses of phytoplankton to higher predators (none with default settings)
       !
       !do n = 1,NUM_PHYTO !{
       !   phyto(n)%jhploss_n(i,j,k) = hp_ingest_vec(n)
       !   phyto(n)%jhploss_p(i,j,k) = phyto(n)%jhploss_n(i,j,k)*prey_p2n_vec(n)
       !   phyto(n)%jhploss_fe(i,j,k) = phyto(n)%jhploss_n(i,j,k)*prey_fe2n_vec(n)
       !   phyto(n)%jhploss_sio2(i,j,k) = phyto(n)%jhploss_n(i,j,k)*prey_si2n_vec(n)
       !enddo !} n
       !
       ! losses of bacteria to higher predators (none with default settings)
       !
       !   bact(1)%jhploss_n(i,j,k) = hp_ingest_vec(4)
       !   bact(1)%jhploss_p(i,j,k) = bact(1)%jhploss_n(i,j,k)*prey_p2n_vec(4)
       !
       ! losses of zooplankton to higher predators
       !
       do n = 1,NUM_ZOO !{
         zoo(n)%jhploss_n(i,j,k) = hp_ingest_vec(NUM_PHYTO+1+n)
         zoo(n)%jhploss_p(i,j,k) = zoo(n)%jhploss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+1+n)
       enddo !} n
       !
       ! losses of detritus to higher predators (none with default settings)
       !
       !cobalt%det_jhploss_n(i,j,k) = hp_ingest_vec(NUM_PHYTO+NUM_ZOO+2)
       !cobalt%det_jhploss_p(i,j,k) = cobalt%det_jhploss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+NUM_ZOO+2)
       !cobalt%det_jhploss_fe(i,j,k) = cobalt%det_jhploss_n(i,j,k)*prey_fe2n_vec(NUM_PHYTO+NUM_ZOO+2)
       !cobalt%det_jhploss_si(i,j,k) = cobalt%det_jhploss_si(i,j,k)*prey_si2n_vec(NUM_PHYTO+NUM_ZOO+2)

    enddo; enddo; enddo  !} i,j,k
    call mpp_clock_end(id_clock_zooplankton_calculations)

    !
    ! 3.2: Plankton foodweb dynamics: Other mortality and loss terms
    !

    call mpp_clock_begin(id_clock_other_losses)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; !{

       !  
       ! 3.2.1 Calculate losses of phytoplankton to aggregation 
       !

       do n = 1,NUM_PHYTO !{
            phyto(n)%jaggloss_n(i,j,k) = phyto(n)%agg*phyto(n)%f_n(i,j,k)**2.0 
            phyto(n)%jaggloss_p(i,j,k) = phyto(n)%jaggloss_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k)
            phyto(n)%jaggloss_fe(i,j,k) = phyto(n)%jaggloss_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k)
            phyto(n)%jaggloss_sio2(i,j,k) = phyto(n)%jaggloss_n(i,j,k)*phyto(n)%q_si_2_n(i,j,k)
       enddo !} n

       !
       ! 3.2.2 Calculate phytoplankton and bacterial losses to viruses
       !

       do n = 1,NUM_PHYTO !{
          phyto(n)%jvirloss_n(i,j,k) = bact(1)%temp_lim(i,j,k)*phyto(n)%vir*phyto(n)%f_n(i,j,k)**2.0 
          phyto(n)%jvirloss_p(i,j,k) = phyto(n)%jvirloss_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k)
          phyto(n)%jvirloss_fe(i,j,k) = phyto(n)%jvirloss_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k)
          phyto(n)%jvirloss_sio2(i,j,k) = phyto(n)%jvirloss_n(i,j,k)*phyto(n)%q_si_2_n(i,j,k)
       enddo !} n

       bact(1)%jvirloss_n(i,j,k) = bact(1)%temp_lim(i,j,k)*bact(1)%vir*bact(1)%f_n(i,j,k)**2.0
       bact(1)%jvirloss_p(i,j,k) = bact(1)%jvirloss_n(i,j,k)*bact(1)%q_p_2_n

       !
       ! 3.2.3 Calculate losses to exudation
       !

       n = DIAZO
       phyto(n)%jexuloss_n(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_no3(i,j,k)+ &
                                    phyto(n)%juptake_nh4(i,j,k)+phyto(n)%juptake_n2(i,j,k),0.0)
       phyto(n)%jexuloss_p(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_po4(i,j,k),0.0)
       phyto(n)%jexuloss_fe(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_fe(i,j,k),0.0)
       do n = 2,NUM_PHYTO !{
          phyto(n)%jexuloss_n(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_no3(i,j,k)+phyto(n)%juptake_nh4(i,j,k),0.0)
          phyto(n)%jexuloss_p(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_po4(i,j,k),0.0)
          phyto(n)%jexuloss_fe(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_fe(i,j,k),0.0)
       enddo
       ! Adjust silica uptake by large phytoplankton downward to maintain constant Si:N stoichimetry
       ! phyto(LARGE)%juptake_sio4(i,j,k) = (1-phyto(LARGE)%exu)*phyto(LARGE)%juptake_sio4(i,j,k)

    enddo; enddo; enddo  !} i,j,k
    call mpp_clock_end(id_clock_other_losses)

    !
    ! 3.3: Plankton foodweb dynamics: Production calculations
    !

    call mpp_clock_begin(id_clock_production_loop)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

       !
       ! 3.3.1: Calculate the production of detritus and dissolved organic material
       !

       ! initialize some cumulative COBALT-wide production diagnostics
       cobalt%jprod_ndet(i,j,k) = 0.0
       cobalt%jprod_pdet(i,j,k) = 0.0
       cobalt%jprod_sldon(i,j,k) = 0.0
       cobalt%jprod_ldon(i,j,k) = 0.0
       cobalt%jprod_srdon(i,j,k) = 0.0
       cobalt%jprod_sldop(i,j,k) = 0.0
       cobalt%jprod_ldop(i,j,k) = 0.0
       cobalt%jprod_srdop(i,j,k) = 0.0
       cobalt%jprod_fedet(i,j,k) = 0.0
       cobalt%jprod_fed(i,j,k) = 0.0
       cobalt%jprod_sidet(i,j,k) = 0.0
       cobalt%jprod_sio4(i,j,k) = 0.0
       cobalt%jprod_po4(i,j,k) = 0.0
       cobalt%jprod_nh4(i,j,k) = 0.0

       !
       ! Production of detritus and dissolved organic material from zooplankton egestion 
       !   

       do m = 1,NUM_ZOO
           zoo(m)%jprod_ndet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_pdet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_sldon(i,j,k) = zoo(m)%phi_sldon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_ldon(i,j,k) = zoo(m)%phi_ldon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_srdon(i,j,k) = zoo(m)%phi_srdon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_sldop(i,j,k) = zoo(m)%phi_sldop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_ldop(i,j,k) = zoo(m)%phi_ldop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_srdop(i,j,k) = zoo(m)%phi_srdop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_fedet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_fe(i,j,k)
           zoo(m)%jprod_sidet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_sio2(i,j,k)


           ! augment cumulative production with zooplankton terms
           cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) + zoo(m)%jprod_ndet(i,j,k)
           cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) + zoo(m)%jprod_pdet(i,j,k)
           cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + zoo(m)%jprod_sldon(i,j,k)
           cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + zoo(m)%jprod_ldon(i,j,k)
           cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + zoo(m)%jprod_srdon(i,j,k)
           cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + zoo(m)%jprod_sldop(i,j,k)
           cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + zoo(m)%jprod_ldop(i,j,k)
           cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + zoo(m)%jprod_srdop(i,j,k)
           cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + zoo(m)%jprod_fedet(i,j,k)
           cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + zoo(m)%jprod_sidet(i,j,k)
       enddo !} m

       !
       ! Production of detritus and dissolved organic material from higher predator egestion 
       ! (did not track individual terms, just add to cumulative total)
       !

       cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%hp_phi_sldon*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + cobalt%hp_phi_ldon*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%hp_phi_srdon*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%hp_phi_sldop*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + cobalt%hp_phi_ldop*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%hp_phi_srdop*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_fe(i,j,k)
       cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_sio2(i,j,k)
       
       !
       ! Sources from phytoplankton aggregation
       !

       do m = 1,NUM_PHYTO
           cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) + phyto(m)%jaggloss_n(i,j,k)
           cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) + phyto(m)%jaggloss_p(i,j,k)
           cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + phyto(m)%jaggloss_fe(i,j,k)
           cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + phyto(m)%jaggloss_sio2(i,j,k)
       enddo !} m

       !
       ! Sources due to phytoplankton mortality from adverse growth conditions (metabolic costs higher than
       ! photosynthetic capacity).  These conditions are assumed to lead to a source of detritus in large
       ! phytoplankton and diazotrophs. 
       !
       !n = DIAZO 
       !cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k),0.0)  
       !cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k),0.0)
       !cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k),0.0)
       !n = LARGE 
       !cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k),0.0)
       !cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k),0.0)
       !cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k),0.0)
       !cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) - min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_si_2_n(i,j,k),0.0)

       !
       ! Sources from viral lysis of phytoplankton (0 in default formulation) and exudation
       !

       do m = 1,NUM_PHYTO
           cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + cobalt%lysis_phi_ldon*phyto(m)%jvirloss_n(i,j,k) + &
                                      phyto(m)%jexuloss_n(i,j,k) 
           cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%lysis_phi_sldon*phyto(m)%jvirloss_n(i,j,k)
           cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%lysis_phi_srdon*phyto(m)%jvirloss_n(i,j,k)
           cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + cobalt%lysis_phi_ldop*phyto(m)%jvirloss_p(i,j,k) + &
                                      phyto(m)%jexuloss_p(i,j,k)
           cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%lysis_phi_sldop*phyto(m)%jvirloss_p(i,j,k)
           cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%lysis_phi_srdop*phyto(m)%jvirloss_p(i,j,k)
           cobalt%jprod_fed(i,j,k)   = cobalt%jprod_fed(i,j,k)   + phyto(m)%jvirloss_fe(i,j,k) + &
                                       phyto(m)%jexuloss_fe(i,j,k) 
           cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + phyto(m)%jvirloss_sio2(i,j,k)
       enddo !} m

       !
       ! Sources of dissolved organic material from small phytoplankton mortality (metabolic costs higher than photosynthetic
       ! capacity).  These conditions are assumed to lead to a lysis-like redistribution of small phyto organic matter.
       !

       !n = SMALL
       !cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) - cobalt%lysis_phi_ldon*min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k),0.0)
       !cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) - cobalt%lysis_phi_sldon*min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k),0.0)
       !cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) - cobalt%lysis_phi_srdon*min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k),0.0)
       !cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) - &
       !                           cobalt%lysis_phi_ldop*min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k),0.0)
       !cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) - &
       !                           cobalt%lysis_phi_sldop*min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k),0.0)
       !cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) - &
       !                           cobalt%lysis_phi_srdop*min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k),0.0)
       !cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) - &
       !                          min(phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k),0.0)

       !
       ! Sources of dissolved organic material from viral lysis due to bacteria 
       !

       cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + cobalt%lysis_phi_ldon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%lysis_phi_sldon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%lysis_phi_srdon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + cobalt%lysis_phi_ldop*bact(1)%jvirloss_p(i,j,k)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%lysis_phi_sldop*bact(1)%jvirloss_p(i,j,k)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%lysis_phi_srdop*bact(1)%jvirloss_p(i,j,k)

       !
       ! Sources of dissolved organic material from bacterial mortality (metabolic costs higher than food uptake).
       ! These conditions are assumed to lead to a lysis-like redistribution of bacteria organic matter.
       !

       cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) - cobalt%lysis_phi_ldon* &
                                  min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) - cobalt%lysis_phi_sldon* &
                                  min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) - cobalt%lysis_phi_srdon* &
                                  min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) - cobalt%lysis_phi_ldop* &
                                  min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) - cobalt%lysis_phi_sldop* &
                                  min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) - cobalt%lysis_phi_srdop* &
                                  min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       !
       ! 3.3.2: Calculate the remineralization of organic material by free-living bacteria
       !

       bact(1)%jprod_nh4(i,j,k) = bact(1)%juptake_ldon(i,j,k) - max(bact(1)%jprod_n(i,j,k),0.0)
       bact(1)%jprod_po4(i,j,k) = bact(1)%juptake_ldop(i,j,k) - &
                                  max(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + bact(1)%jprod_nh4(i,j,k)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + bact(1)%jprod_po4(i,j,k)

       !
       ! 3.3.3: Zooplankton production and excretion calculations
       !

       do m = 1,NUM_ZOO

          ingest_p2n = zoo(m)%jingest_p(i,j,k)/(zoo(m)%jingest_n(i,j,k)+epsln)

          if (ingest_p2n .lt. zoo(m)%q_p_2_n) then
             zoo(m)%jprod_n(i,j,k) = zoo(m)%gge_max*zoo(m)%jingest_p(i,j,k)*(1.0/zoo(m)%q_p_2_n)
          else
             zoo(m)%jprod_n(i,j,k) = zoo(m)%gge_max*zoo(m)%jingest_n(i,j,k)
          endif

          ! adjust production terms for basal respiration costs
          !if (zoo(m)%f_n(i,j,k).gt.refuge_conc) then
            zoo(m)%jprod_n(i,j,k) = zoo(m)%jprod_n(i,j,k) - &
                                     zoo(m)%f_n(i,j,k)/(refuge_conc + zoo(m)%f_n(i,j,k))* &
                                     zoo(m)%temp_lim(i,j,k)*zoo(m)%bresp*zoo(m)%f_n(i,j,k)
          !endif
          !
          ! Ingested material that does not go to zooplankton production, detrital production
          ! or production of dissolved organic material is excreted as nh4 or po4.  If production
          ! is negative, zooplankton are lost to large detritus 
          !
          if (zoo(m)%jprod_n(i,j,k) .gt. 0.0) then 
             zoo(m)%jprod_nh4(i,j,k) =  zoo(m)%jingest_n(i,j,k) - zoo(m)%jprod_ndet(i,j,k) -  &
                                        zoo(m)%jprod_n(i,j,k) - zoo(m)%jprod_ldon(i,j,k) - &
                                        zoo(m)%jprod_sldon(i,j,k) - zoo(m)%jprod_srdon(i,j,k)
             zoo(m)%jprod_po4(i,j,k) =  zoo(m)%jingest_p(i,j,k) - zoo(m)%jprod_pdet(i,j,k) - & 
                                        zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n - zoo(m)%jprod_ldop(i,j,k) -  &
                                        zoo(m)%jprod_sldop(i,j,k) - zoo(m)%jprod_srdop(i,j,k)
          else
             ! None of the ingestion material goes to zooplankton production
             zoo(m)%jprod_nh4(i,j,k) =  zoo(m)%jingest_n(i,j,k) - zoo(m)%jprod_ndet(i,j,k) - & 
                                        zoo(m)%jprod_ldon(i,j,k) - zoo(m)%jprod_sldon(i,j,k) - & 
                                        zoo(m)%jprod_srdon(i,j,k)
             zoo(m)%jprod_po4(i,j,k) =  zoo(m)%jingest_p(i,j,k) - zoo(m)%jprod_pdet(i,j,k) - &
                                        zoo(m)%jprod_ldop(i,j,k) - zoo(m)%jprod_sldop(i,j,k) - & 
                                        zoo(m)%jprod_srdop(i,j,k)

             ! The negative production (i.e., mortality) is lost to large detritus. Update values
             ! for zooplankton and for total.

             zoo(m)%jprod_ndet(i,j,k) = zoo(m)%jprod_ndet(i,j,k) - zoo(m)%jprod_n(i,j,k)
             zoo(m)%jprod_pdet(i,j,k) = zoo(m)%jprod_pdet(i,j,k) - zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n
             cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - zoo(m)%jprod_n(i,j,k)
             cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n 
          endif

          ! cumulative production of inorganic nutrients 
          cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + zoo(m)%jprod_nh4(i,j,k)
          cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + zoo(m)%jprod_po4(i,j,k)

          !
          ! Any ingested iron that is not allocated to detritus is routed back to the
          ! dissolved pool.       
          !
          zoo(m)%jprod_fed(i,j,k) = (1.0 - zoo(m)%phi_det)*zoo(m)%jingest_fe(i,j,k)
          cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + zoo(m)%jprod_fed(i,j,k)
          !
          ! Any ingested opal that is not allocated to detritus is assumed to undergo
          ! rapid dissolution to dissolved silica
          !
          zoo(m)%jprod_sio4(i,j,k) = (1.0 - zoo(m)%phi_det)*zoo(m)%jingest_sio2(i,j,k)
          cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + zoo(m)%jprod_sio4(i,j,k)
 
       enddo !} m

       !
       ! Excretion by higher predators
       !
       cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_fe(i,j,k)
       cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_sio2(i,j,k)
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + cobalt%hp_phi_nh4*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + cobalt%hp_phi_po4*cobalt%hp_jingest_p(i,j,k)

    enddo; enddo ; enddo !} i,j,k
    call mpp_clock_end(id_clock_production_loop)

    call mpp_clock_begin(id_clock_ballast_loops)
    do j = jsc, jec ; do i = isc, iec   !{
       cobalt%zt(i,j,1) = dzt(i,j,1)
       cobalt%zm(i,j,1) = 0.5*dzt(i,j,1)
    enddo; enddo !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%zt(i,j,k) = cobalt%zt(i,j,k-1) + dzt(i,j,k)
       cobalt%zm(i,j,k) = cobalt%zm(i,j,k-1) + dzt(i,j,k)
    enddo; enddo ; enddo !} i,j,k

!
!------------------------------------------------------------------------------------
! 4: Production of calcium carbonate (Calcite and Aragonite) and lithogenic material
!------------------------------------------------------------------------------------
!
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

    !
    ! 4.1: Calculate aragonite and calcite saturation states
    !

       TK = Temp(i,j,k) + 273.15
       PRESS = 0.1016 * cobalt%zt(i,j,k) + 1.013
       PKSPA = 171.945 + 0.077993 * TK - 2903.293 / TK - 71.595 * log10(TK) - (-0.068393 + 1.7276e-3 * &
          TK + 88.135 / TK) * sqrt(max(epsln, Salt(i,j,k))) + 0.10018 * max(epsln, Salt(i,j,k)) -      &
          5.9415e-3 * max(epsln, Salt(i,j,k))**(1.5) - 0.02 - (48.76 - 2.8 - 0.5304 * Temp(i,j,k)) *   &
          (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 - 0.3692 * Temp(i,j,k))) * (PRESS - 1.013) *&
          (PRESS - 1.013) / (382.92 * TK)
       cobalt%co3_sol_arag(i,j,k) = 10**(-PKSPA) / (2.937d-4 * max(5.0, Salt(i,j,k)))
       cobalt%omega_arag(i,j,k) = cobalt%f_co3_ion(i,j,k) / cobalt%co3_sol_arag(i,j,k)
       PKSPC = 171.9065 + 0.077993 * TK - 2839.319 / TK - 71.595 * log10(TK) - (-0.77712 + 2.8426e-3 * &
          TK + 178.34 / TK) * sqrt(max(epsln, Salt(i,j,k))) + 0.07711 * max(epsln, Salt(i,j,k)) -      &
          4.1249e-3 * max(epsln, Salt(i,j,k))**(1.5) - 0.02 - (48.76 - 0.5304 * Temp(i,j,k)) *         &
          (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 - 0.3692 * Temp(i,j,k))) * (PRESS - 1.013) *&
          (PRESS - 1.013) / (382.92 * TK)
       cobalt%co3_sol_calc(i,j,k) = 10**(-PKSPC) / (2.937d-4 * max(5.0, Salt(i,j,k)))
       cobalt%omega_calc(i,j,k) = cobalt%f_co3_ion(i,j,k) / cobalt%co3_sol_calc(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    !
    ! 4.2: Calculate the production rate of aragonite and calcite detritus 
    !

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
        cobalt%jprod_cadet_arag(i,j,k) = (zoo(2)%jzloss_n(i,j,k) + zoo(3)%jzloss_n(i,j,k) + &
                       zoo(2)%jhploss_n(i,j,k) + zoo(3)%jhploss_n(i,j,k))*cobalt%ca_2_n_arag* &
                       min(cobalt%caco3_sat_max, max(0.0,cobalt%omega_arag(i,j,k) - 1.0)) + epsln
        cobalt%jprod_cadet_calc(i,j,k) = (zoo(1)%jzloss_n(i,j,k) + phyto(SMALL)%jaggloss_n(i,j,k))*cobalt%ca_2_n_calc* &
                       min(cobalt%caco3_sat_max, max(0.0, cobalt%omega_calc(i,j,k) - 1.0)) + epsln
    enddo; enddo ; enddo !} i,j,k

    !
    ! 4.3: Lithogenic detritus production (repackaged from f_lith during filter feeding)
    !

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%jprod_lithdet(i,j,k)=( cobalt%total_filter_feeding(i,j,k)/ &
                                   ( phyto(LARGE)%f_n(i,j,k) + phyto(DIAZO)%f_n(i,j,k) + epsln ) * &  
                                    cobalt%phi_lith + cobalt%k_lith ) * cobalt%f_lith(i,j,k)
    enddo; enddo ; enddo !} i,j,k

!
!---------------------------------------------------------------------------------------------------------
! 5: Detrital dissolution and remineralization calculation
!---------------------------------------------------------------------------------------------------------
!

    !
    ! 5.1: Dissolution of aragonite, calcite and opal detrital particles
    !

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       cobalt%jdiss_cadet_arag(i,j,k) = cobalt%gamma_cadet_arag * & 
         max(0.0, 1.0 - cobalt%omega_arag(i,j,k)) * cobalt%f_cadet_arag(i,j,k)
       cobalt%jdiss_cadet_calc(i,j,k) = cobalt%gamma_cadet_calc * &
         max(0.0, 1.0 - cobalt%omega_calc(i,j,k)) * cobalt%f_cadet_calc(i,j,k)
       cobalt%jdiss_sidet(i,j,k) = cobalt%gamma_sidet * cobalt%f_sidet(i,j,k)
       cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + cobalt%jdiss_sidet(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    !
    ! 5.2: Remineralization of nitrogen, phosphorous and iron detritus accounting for oxygen 
    !      and mineral protection 
    !

    do k=1,nk ; do j=jsc,jec ; do i=isc,iec  !{
       cobalt%jno3denit_wc(i,j,k) = 0.0
       !
       !   Under oxic conditions
       !
       if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
          cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet * cobalt%f_o2(i,j,k) / & 
               ( cobalt%k_o2 + cobalt%f_o2(i,j,k) )*max( 0.0, cobalt%f_ndet(i,j,k) - &
               cobalt%rpcaco3*(cobalt%f_cadet_arag(i,j,k) + cobalt%f_cadet_calc(i,j,k)) - & 
               cobalt%rplith*cobalt%f_lithdet(i,j,k) - cobalt%rpsio2*cobalt%f_sidet(i,j,k) )
       !
       ! Under sub-oxic conditions
       !
       else !}{
          cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet * cobalt%o2_min / &
               (cobalt%k_o2 + cobalt%o2_min)* &
               cobalt%f_no3(i,j,k) / (phyto(SMALL)%k_no3 + cobalt%f_no3(i,j,k))* &
               max(0.0, cobalt%f_ndet(i,j,k) - &
               cobalt%rpcaco3*(cobalt%f_cadet_arag(i,j,k) + cobalt%f_cadet_calc(i,j,k)) - &
               cobalt%rplith*cobalt%f_lithdet(i,j,k) + cobalt%rpsio2*cobalt%f_sidet(i,j,k) )
          cobalt%jno3denit_wc(i,j,k) = cobalt%jremin_ndet(i,j,k) * cobalt%n_2_n_denit
       endif !}
       !
       ! P and Fe assumed to be protected similarly to N
       !
       cobalt%jremin_pdet(i,j,k) = cobalt%jremin_ndet(i,j,k)/(cobalt%f_ndet(i,j,k) + epsln)* &
         cobalt%f_pdet(i,j,k)
       cobalt%jremin_fedet(i,j,k) = cobalt%jremin_ndet(i,j,k) / (cobalt%f_ndet(i,j,k) + epsln) * &
         cobalt%remin_eff_fedet*cobalt%f_fedet(i,j,k)
    enddo; enddo; enddo  !} i,j,k

!
!--------------------------------------------------------------------------------------------
! 6: Miscellaneous sources and sinks: Nitrification, Iron Scavenging, Coastal Iron inputs
!--------------------------------------------------------------------------------------------
!

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

       !
       !  Nitrification
       !

       cobalt%jnitrif(i,j,k) = cobalt%gamma_nitrif * cobalt%expkT(i,j,k) * cobalt%f_nh4(i,j,k) * &
            phyto(SMALL)%nh4lim(i,j,k) * (1.0 - cobalt%f_irr_mem(i,j,k) / &
            (cobalt%irr_inhibit + cobalt%f_irr_mem(i,j,k)))

       !
       ! Solve for free iron
       !

       !cobalt%kfe_eq_lig(i,j,k) = 10**( log10(cobalt%kfe_eq_lig_ll) - &
       !     ( cobalt%irr_inst(i,j,k)/(cobalt%ki_fescav+cobalt%irr_inst(i,j,k)) ) * &
       !     (log10(cobalt%kfe_eq_lig_ll) - log10(cobalt%kfe_eq_lig_hl)) )
       cobalt%kfe_eq_lig(i,j,k) = min(cobalt%kfe_eq_lig_ll, 10**( log10(cobalt%kfe_eq_lig_hl) + &
            max(0.0,cobalt%gamma_fescav*log10(cobalt%io_fescav/cobalt%irr_inst(i,j,k))) ) ) 

       feprime = 1.0 + cobalt%kfe_eq_lig(i,j,k) * (cobalt%felig_bkg + cobalt%felig_2_don * &
            (cobalt%f_sldon(i,j,k) + cobalt%f_srdon(i,j,k)) - cobalt%f_fed(i,j,k))
       feprime = (-feprime + (feprime * feprime + 4.0 * cobalt%kfe_eq_lig(i,j,k) * &
            cobalt%f_fed(i,j,k))**(0.5)) / (2.0 * cobalt%kfe_eq_lig(i,j,k))

       !
       ! Iron adsorption to detrital particles
       !

       cobalt%jfe_ads(i,j,k) = min(r_dt,cobalt%alpha_fescav*feprime)
       if (cobalt%f_fed(i,j,k).gt.1.0e-9) then !{
          cobalt%jfe_ads(i,j,k) = min(r_dt,5.0*cobalt%alpha_fescav*cobalt%f_fed(i,j,k))
       endif !}
       !
       ! Coastal iron inputs (proxy for sediment inputs for areas with poorly resolved shelves)
       !

       cobalt%jfe_coast(i,j,k) = cobalt%fe_coast * mask_coast(i,j) * grid_tmask(i,j,k) / &
            sqrt(grid_dat(i,j))

    enddo; enddo; enddo  !} i,j,k

!
!-------------------------------------------------------------------------------------------------
! 7: Sedimentary fluxes/transformations
!-------------------------------------------------------------------------------------------------
!
    do j = jsc, jec; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{
          !
          ! Nitrogen flux from the sediments
          ! 
          if (cobalt%f_ndet_btf(i,j,1) .gt. 0.0) then !{
             ! fpoc_bottom in mmoles C m-2 day-1 for burial relationship
             fpoc_btm = (cobalt%f_ndet_btf(i,j,1)*cobalt%c_2_n*sperd*1000.0)
             cobalt%frac_burial(i,j) = (0.013 + 0.53*fpoc_btm**2.0)/((7.0+fpoc_btm)**2.0)
             cobalt%fndet_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%f_ndet_btf(i,j,1)
             cobalt%fpdet_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%f_pdet_btf(i,j,1)
             ! fpoc_bottom in micromoles C cm-2 day-1 for denitrification relationship, cap at 43
             ! to prevent anomalous extrapolation of the relationship
             log_fpoc_btm = log(min(43.0,0.1*fpoc_btm))
             cobalt%fno3denit_sed(i,j) = min(cobalt%f_no3(i,j,k)*cobalt%Rho_0*r_dt,  &      
                  min((cobalt%f_ndet_btf(i,j,1)-cobalt%fndet_burial(i,j))*cobalt%n_2_n_denit, & 
                  10.0**(-0.9543+0.7662*log_fpoc_btm - 0.235*log_fpoc_btm**2.0)/(cobalt%c_2_n*sperd*100.0)* &
                  cobalt%n_2_n_denit*cobalt%f_no3(i,j,k)/(cobalt%k_no3_denit + cobalt%f_no3(i,j,k))))
             if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
                cobalt%fnoxic_sed(i,j) = max(0.0, min(cobalt%f_o2(i,j,k)*cobalt%Rho_0*r_dt*(1.0/cobalt%o2_2_nh4), &
                                         cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j) - &
                                         cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit))
             else
                cobalt%fnoxic_sed(i,j) = 0.0
             endif !}
             cobalt%fno3denit_sed(i,j) = cobalt%fno3denit_sed(i,j) + &
                                         min(cobalt%f_no3(i,j,k)*cobalt%Rho_0*r_dt-cobalt%fno3denit_sed(i,j), &
                                         (cobalt%f_ndet_btf(i,j,1)-cobalt%fnoxic_sed(i,j)-cobalt%fndet_burial(i,j) - &
                                         cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)*cobalt%n_2_n_denit)
             cobalt%fnfeso4red_sed(i,j) = max(0.0, cobalt%f_ndet_btf(i,j,1)-cobalt%fnoxic_sed(i,j)- &
                                          cobalt%fndet_burial(i,j)-cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)
          else
             cobalt%fnfeso4red_sed(i,j) = 0.0
             cobalt%fno3denit_sed(i,j) = 0.0
             cobalt%fnoxic_sed(i,j) = 0.0
          endif !}

          ! iron from sediment 
          cobalt%ffe_sed(i,j) = cobalt%fe_2_n_sed * cobalt%f_ndet_btf(i,j,1)

          !
          ! Calcium carbonate flux and burial
          !
          cobalt%fcased_redis(i,j) = max(0.0, min(0.5 * cobalt%f_cased(i,j,1) * r_dt, min(0.5 *       &
             cobalt%f_cadet_calc_btf(i,j,1), 0.165 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n) +        &
             0.1244 / spery * max(0.0, 1.0 - cobalt%omega_calc(i,j,k) +   &
             4.38 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n * spery)**(2.91) *                        &
             max(1.0, cobalt%f_lithdet_btf(i,j,1) * spery + cobalt%f_cadet_calc_btf(i,j,1) * 100.0 *  &
             spery)**(-2.55) * cobalt%f_cased(i,j,1)))
          cobalt%fcased_burial(i,j) = max(0.0, cobalt%f_cadet_calc_btf(i,j,1) * cobalt%f_cased(i,j,1) /&
             8.1e3)
          cobalt%f_cased(i,j,1) = cobalt%f_cased(i,j,1) + (cobalt%f_cadet_calc_btf(i,j,1) -            &
             cobalt%fcased_redis(i,j) - cobalt%fcased_burial(i,j)) / cobalt%z_sed * dt *               &
             grid_tmask(i,j,k)
          !
          ! Bottom flux boundaries passed to the vertical mixing routine 
          !
          cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
             (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) + cobalt%alk_2_n_denit * cobalt%fno3denit_sed(i,j)
          cobalt%b_dic(i,j) =  - cobalt%fcased_redis(i,j) - cobalt%f_cadet_arag_btf(i,j,1) -            &
             (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) * cobalt%c_2_n
          cobalt%b_fed(i,j) = - cobalt%ffe_sed(i,j)
          cobalt%b_nh4(i,j) = - cobalt%f_ndet_btf(i,j,1) + cobalt%fndet_burial(i,j)
          cobalt%b_no3(i,j) = cobalt%fno3denit_sed(i,j)
          cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * (cobalt%fnoxic_sed(i,j) + cobalt%fnfeso4red_sed(i,j))
          cobalt%b_po4(i,j) = - cobalt%f_pdet_btf(i,j,1) + cobalt%fpdet_burial(i,j)
          cobalt%b_sio4(i,j)= - cobalt%f_sidet_btf(i,j,1)

       endif !}
    enddo; enddo  !} i, j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%f_cased(i,j,k) = 0.0
    enddo; enddo ; enddo  !} i,j,k

    call mpp_clock_end(id_clock_ballast_loops)

    call g_tracer_set_values(tracer_list,'alk',  'btf', cobalt%b_alk ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic',  'btf', cobalt%b_dic ,isd,jsd)
    call g_tracer_set_values(tracer_list,'fed',  'btf', cobalt%b_fed ,isd,jsd)
    call g_tracer_set_values(tracer_list,'nh4',  'btf', cobalt%b_nh4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'no3',  'btf', cobalt%b_no3 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2',   'btf', cobalt%b_o2  ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4',  'btf', cobalt%b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'sio4', 'btf', cobalt%b_sio4,isd,jsd)

    call mpp_clock_begin(id_clock_source_sink_loop1)

!
!-----------------------------------------------------------------------
! 8: Source/sink calculations 
!-----------------------------------------------------------------------
!
    !  
    !-------------------------------------------------------------------
    ! 4.1: Update the prognostics tracer fields via their pointers.
    !-------------------------------------------------------------------
    !
    call g_tracer_get_pointer(tracer_list,'alk'    ,'field',cobalt%p_alk    )
    call g_tracer_get_pointer(tracer_list,'cadet_arag','field',cobalt%p_cadet_arag)
    call g_tracer_get_pointer(tracer_list,'cadet_calc','field',cobalt%p_cadet_calc)
    call g_tracer_get_pointer(tracer_list,'dic'    ,'field',cobalt%p_dic    )
    call g_tracer_get_pointer(tracer_list,'fed'    ,'field',cobalt%p_fed    )
    call g_tracer_get_pointer(tracer_list,'fedi'   ,'field',cobalt%p_fedi   )
    call g_tracer_get_pointer(tracer_list,'felg'   ,'field',cobalt%p_felg   )
    call g_tracer_get_pointer(tracer_list,'fesm'   ,'field',cobalt%p_fesm )
    call g_tracer_get_pointer(tracer_list,'fedet'  ,'field',cobalt%p_fedet  )
    call g_tracer_get_pointer(tracer_list,'ldon'   ,'field',cobalt%p_ldon   )
    call g_tracer_get_pointer(tracer_list,'ldop'   ,'field',cobalt%p_ldop   )
    call g_tracer_get_pointer(tracer_list,'lith'   ,'field',cobalt%p_lith   )
    call g_tracer_get_pointer(tracer_list,'lithdet','field',cobalt%p_lithdet)
    call g_tracer_get_pointer(tracer_list,'nbact'  ,'field',cobalt%p_nbact  )
    call g_tracer_get_pointer(tracer_list,'ndet'   ,'field',cobalt%p_ndet   )
    call g_tracer_get_pointer(tracer_list,'ndi'    ,'field',cobalt%p_ndi    )
    call g_tracer_get_pointer(tracer_list,'nlg'    ,'field',cobalt%p_nlg    )
    call g_tracer_get_pointer(tracer_list,'nsm' ,'field',cobalt%p_nsm )
    call g_tracer_get_pointer(tracer_list,'nh4'    ,'field',cobalt%p_nh4    )
    call g_tracer_get_pointer(tracer_list,'no3'    ,'field',cobalt%p_no3    )
    call g_tracer_get_pointer(tracer_list,'o2'     ,'field',cobalt%p_o2     )
    call g_tracer_get_pointer(tracer_list,'pdet'   ,'field',cobalt%p_pdet   )
    call g_tracer_get_pointer(tracer_list,'po4'    ,'field',cobalt%p_po4    )
    call g_tracer_get_pointer(tracer_list,'srdon'   ,'field',cobalt%p_srdon   )
    call g_tracer_get_pointer(tracer_list,'srdop'   ,'field',cobalt%p_srdop   )
    call g_tracer_get_pointer(tracer_list,'sldon'   ,'field',cobalt%p_sldon   )
    call g_tracer_get_pointer(tracer_list,'sldop'   ,'field',cobalt%p_sldop   )
    call g_tracer_get_pointer(tracer_list,'sidet'  ,'field',cobalt%p_sidet  )
    call g_tracer_get_pointer(tracer_list,'silg'   ,'field',cobalt%p_silg   )
    call g_tracer_get_pointer(tracer_list,'sio4'   ,'field',cobalt%p_sio4   )
    call g_tracer_get_pointer(tracer_list,'nsmz'   ,'field',cobalt%p_nsmz   )
    call g_tracer_get_pointer(tracer_list,'nmdz'   ,'field',cobalt%p_nmdz   )
    call g_tracer_get_pointer(tracer_list,'nlgz'   ,'field',cobalt%p_nlgz   )


    if (cobalt%id_no3_in_source .gt. 0)                &
         used = send_data(cobalt%id_no3_in_source,         cobalt%f_no3,          &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    call mpp_clock_end(id_clock_source_sink_loop1)
    !
    !-----------------------------------------------------------------------
    ! 4.2: Source sink calculations
    !-----------------------------------------------------------------------
    !
    !     Phytoplankton Nitrogen and Phosphorus
    !
    call mpp_clock_begin(id_clock_source_sink_loop2)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Diazotrophic Phytoplankton Nitrogen
       !
       cobalt%jndi(i,j,k) = phyto(DIAZO)%mu(i,j,k)*phyto(DIAZO)%f_n(i,j,k) - &
                            phyto(DIAZO)%jzloss_n(i,j,k) -       &
                            phyto(DIAZO)%jhploss_n(i,j,k) - phyto(DIAZO)%jaggloss_n(i,j,k) -       &
                            phyto(DIAZO)%jvirloss_n(i,j,k) - phyto(DIAZO)%jexuloss_n(i,j,k)
       cobalt%p_ndi(i,j,k,tau) = cobalt%p_ndi(i,j,k,tau) + cobalt%jndi(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Nitrogen
       !
       cobalt%jnlg(i,j,k) = phyto(LARGE)%mu(i,j,k)*phyto(LARGE)%f_n(i,j,k) -    &
                            phyto(LARGE)%jzloss_n(i,j,k) - phyto(LARGE)%jhploss_n(i,j,k) -         &
                            phyto(LARGE)%jaggloss_n(i,j,k) - phyto(LARGE)%jvirloss_n(i,j,k) -      &
                            phyto(LARGE)%jexuloss_n(i,j,k)
       cobalt%p_nlg(i,j,k,tau) = cobalt%p_nlg(i,j,k,tau) + cobalt%jnlg(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Nitrogen
       !
       cobalt%jnsm(i,j,k) = phyto(SMALL)%mu(i,j,k)*phyto(SMALL)%f_n(i,j,k) -    &
                            phyto(SMALL)%jzloss_n(i,j,k) - phyto(SMALL)%jhploss_n(i,j,k) -         &
                            phyto(SMALL)%jaggloss_n(i,j,k) - phyto(SMALL)%jvirloss_n(i,j,k) -      &
                            phyto(SMALL)%jexuloss_n(i,j,k)                                         
       cobalt%p_nsm(i,j,k,tau) = cobalt%p_nsm(i,j,k,tau) + cobalt%jnsm(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    call mpp_clock_end(id_clock_source_sink_loop2)
    !
    !     Phytoplankton Silicon and Iron
    !
    call mpp_clock_begin(id_clock_source_sink_loop3)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Large Phytoplankton Silicon
       !
       cobalt%jsilg(i,j,k) = phyto(LARGE)%juptake_sio4(i,j,k) - & 
                             phyto(LARGE)%jzloss_sio2(i,j,k) - phyto(LARGE)%jhploss_sio2(i,j,k) - &
                             phyto(LARGE)%jaggloss_sio2(i,j,k) - phyto(LARGE)%jvirloss_sio2(i,j,k)
       cobalt%p_silg(i,j,k,tau) = cobalt%p_silg(i,j,k,tau) + cobalt%jsilg(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Diazotrophic Phytoplankton Iron
       !
       cobalt%jfedi(i,j,k) = phyto(DIAZO)%juptake_fe(i,j,k) - &
                             phyto(DIAZO)%jzloss_fe(i,j,k) - &
                             phyto(DIAZO)%jhploss_fe(i,j,k) - phyto(DIAZO)%jaggloss_fe(i,j,k) - &
                             phyto(DIAZO)%jvirloss_fe(i,j,k) - phyto(DIAZO)%jexuloss_fe(i,j,k)
       cobalt%p_fedi(i,j,k,tau) = cobalt%p_fedi(i,j,k,tau) + cobalt%jfedi(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Iron
       !
       cobalt%jfelg(i,j,k) = phyto(LARGE)%juptake_fe(i,j,k) - & 
                             phyto(LARGE)%jzloss_fe(i,j,k) - &
                             phyto(LARGE)%jhploss_fe(i,j,k) - phyto(LARGE)%jaggloss_fe(i,j,k) - &
                             phyto(LARGE)%jvirloss_fe(i,j,k) - phyto(LARGE)%jexuloss_fe(i,j,k)
       cobalt%p_felg(i,j,k,tau) = cobalt%p_felg(i,j,k,tau) + cobalt%jfelg(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Iron
       !
       cobalt%jfesm(i,j,k) = phyto(SMALL)%juptake_fe(i,j,k) - &
                                phyto(SMALL)%jzloss_fe(i,j,k) - &
                                phyto(SMALL)%jhploss_fe(i,j,k) - phyto(SMALL)%jaggloss_fe(i,j,k) - &
                                phyto(SMALL)%jvirloss_fe(i,j,k) - phyto(SMALL)%jexuloss_fe(i,j,k)
       cobalt%p_fesm(i,j,k,tau) = cobalt%p_fesm(i,j,k,tau) + cobalt%jfesm(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Bacteria
       !
       cobalt%jnbact(i,j,k) = bact(1)%jprod_n(i,j,k) - bact(1)%jzloss_n(i,j,k) - &
                              bact(1)%jvirloss_n(i,j,k) - bact(1)%jhploss_n(i,j,k)  
       cobalt%p_nbact(i,j,k,tau) = cobalt%p_nbact(i,j,k,tau) + cobalt%jnbact(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    call mpp_clock_end(id_clock_source_sink_loop3)
    !
    !    Zooplankton 
    !
    call mpp_clock_begin(id_clock_source_sink_loop4)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Small zooplankton
       !
       cobalt%jnsmz(i,j,k) = zoo(1)%jprod_n(i,j,k) - zoo(1)%jzloss_n(i,j,k) - &
                             zoo(1)%jhploss_n(i,j,k)
       cobalt%p_nsmz(i,j,k,tau) = cobalt%p_nsmz(i,j,k,tau) + cobalt%jnsmz(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Medium zooplankton
       !
       cobalt%jnmdz(i,j,k) = zoo(2)%jprod_n(i,j,k) - zoo(2)%jzloss_n(i,j,k) - &
                             zoo(2)%jhploss_n(i,j,k)
       cobalt%p_nmdz(i,j,k,tau) = cobalt%p_nmdz(i,j,k,tau) + cobalt%jnmdz(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Large zooplankton
       !
       cobalt%jnlgz(i,j,k) = zoo(3)%jprod_n(i,j,k) - zoo(3)%jzloss_n(i,j,k) - &
                             zoo(3)%jhploss_n(i,j,k)
       cobalt%p_nlgz(i,j,k,tau) = cobalt%p_nlgz(i,j,k,tau) + cobalt%jnlgz(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    call mpp_clock_end(id_clock_source_sink_loop4)
    !
    !     NO3
    !
    call mpp_clock_begin(id_clock_source_sink_loop5)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       cobalt%jno3(i,j,k) =  cobalt%jnitrif(i,j,k) - phyto(DIAZO)%juptake_no3(i,j,k) -  &
                             phyto(LARGE)%juptake_no3(i,j,k) - phyto(SMALL)%juptake_no3(i,j,k) - &
                             cobalt%jno3denit_wc(i,j,k)
       cobalt%p_no3(i,j,k,tau) = cobalt%p_no3(i,j,k,tau) + cobalt%jno3(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     Other nutrients
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! NH4
       !
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + cobalt%jremin_ndet(i,j,k)
       cobalt%jnh4(i,j,k) = cobalt%jprod_nh4(i,j,k) - phyto(DIAZO)%juptake_nh4(i,j,k) - &
                            phyto(LARGE)%juptake_nh4(i,j,k) - phyto(SMALL)%juptake_nh4(i,j,k) - &
                            cobalt%jnitrif(i,j,k)
       cobalt%p_nh4(i,j,k,tau) = cobalt%p_nh4(i,j,k,tau) + cobalt%jnh4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! PO4
       !
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + cobalt%jremin_pdet(i,j,k) 
       cobalt%jpo4(i,j,k) = cobalt%jprod_po4(i,j,k) - phyto(DIAZO)%juptake_po4(i,j,k) - &
                            phyto(LARGE)%juptake_po4(i,j,k) - phyto(SMALL)%juptake_po4(i,j,k)
       cobalt%p_po4(i,j,k,tau) = cobalt%p_po4(i,j,k,tau) + cobalt%jpo4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! SiO4
       !
       cobalt%jsio4(i,j,k) = cobalt%jprod_sio4(i,j,k) - phyto(LARGE)%juptake_sio4(i,j,k)
       cobalt%p_sio4(i,j,k,tau) = cobalt%p_sio4(i,j,k,tau) + cobalt%jsio4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! Fed
       !
       cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + &
                                 cobalt%jremin_fedet(i,j,k) + cobalt%jfe_coast(i,j,k)  
       cobalt%jfed(i,j,k) = cobalt%jprod_fed(i,j,k) - phyto(DIAZO)%juptake_fe(i,j,k) - &
                            phyto(LARGE)%juptake_fe(i,j,k) -  phyto(SMALL)%juptake_fe(i,j,k) - &
                            cobalt%jfe_ads(i,j,k)
       cobalt%p_fed(i,j,k,tau) = cobalt%p_fed(i,j,k,tau) + cobalt%jfed(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

    call mpp_clock_end(id_clock_source_sink_loop5)
    !
    !-----------------------------------------------------------------------
    !     Detrital Components
    !-----------------------------------------------------------------------
    !
    call mpp_clock_begin(id_clock_source_sink_loop6)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Cadet_arag
       !
       cobalt%jcadet_arag(i,j,k) = cobalt%jprod_cadet_arag(i,j,k) - cobalt%jdiss_cadet_arag(i,j,k) 
       cobalt%p_cadet_arag(i,j,k,tau) = cobalt%p_cadet_arag(i,j,k,tau) + cobalt%jcadet_arag(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Cadet_calc
       !
       cobalt%jcadet_calc(i,j,k) = cobalt%jprod_cadet_calc(i,j,k) - cobalt%jdiss_cadet_calc(i,j,k)
       cobalt%p_cadet_calc(i,j,k,tau) = cobalt%p_cadet_calc(i,j,k,tau) + cobalt%jcadet_calc(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Fedet
       !
       cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + cobalt%jfe_ads(i,j,k)
       cobalt%jfedet(i,j,k) = cobalt%jprod_fedet(i,j,k) - &
                              cobalt%jremin_fedet(i,j,k) - cobalt%det_jzloss_fe(i,j,k) - & 
                              cobalt%det_jhploss_fe(i,j,k)
       cobalt%p_fedet(i,j,k,tau) = cobalt%p_fedet(i,j,k,tau) + cobalt%jfedet(i,j,k)*dt*grid_tmask(i,j,k) 
       !
       ! Lithdet
       !
       cobalt%jlithdet(i,j,k) = cobalt%jprod_lithdet(i,j,k) 
       cobalt%p_lithdet(i,j,k,tau) = cobalt%p_lithdet(i,j,k,tau) + cobalt%jlithdet(i,j,k) * dt *  &
                                     grid_tmask(i,j,k)
       !
       ! Ndet
       !
       cobalt%jndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - cobalt%jremin_ndet(i,j,k) - &
                             cobalt%det_jzloss_n(i,j,k) - cobalt%det_jhploss_n(i,j,k)
       cobalt%p_ndet(i,j,k,tau) = cobalt%p_ndet(i,j,k,tau) + cobalt%jndet(i,j,k)*dt*grid_tmask(i,j,k)
       !cobalt%p_ndet(i,j,k,tau) = max(cobalt%p_ndet(i,j,k,tau),0.0)
       !
       ! Pdet
       !
       cobalt%jpdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - cobalt%jremin_pdet(i,j,k) - &
                             cobalt%det_jzloss_p(i,j,k) - cobalt%det_jhploss_p(i,j,k)         
       cobalt%p_pdet(i,j,k,tau) = cobalt%p_pdet(i,j,k,tau) + cobalt%jpdet(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Sidet
       !
       cobalt%jsidet(i,j,k) = cobalt%jprod_sidet(i,j,k) - & 
                              cobalt%jdiss_sidet(i,j,k) - cobalt%det_jzloss_si(i,j,k) - &
                              cobalt%det_jhploss_si(i,j,k)
       cobalt%p_sidet(i,j,k,tau) = cobalt%p_sidet(i,j,k,tau) + cobalt%jsidet(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     Dissolved Organic Matter
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Labile Dissolved Organic Nitrogen
       !
       cobalt%jldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + &
                             cobalt%gamma_sldon*cobalt%f_sldon(i,j,k) + &
                             cobalt%gamma_srdon*cobalt%f_srdon(i,j,k) - bact(1)%juptake_ldon(i,j,k)
       cobalt%p_ldon(i,j,k,tau) = cobalt%p_ldon(i,j,k,tau) +  cobalt%jldon(i,j,k)*dt*               &
            grid_tmask(i,j,k)
       !
       ! Labile Dissolved Organic Phosphorous
       !
       cobalt%jldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + &
                             cobalt%gamma_sldop*cobalt%f_sldop(i,j,k) + &
                             cobalt%gamma_srdop*cobalt%f_srdop(i,j,k) - bact(1)%juptake_ldop(i,j,k)
       cobalt%p_ldop(i,j,k,tau) = cobalt%p_ldop(i,j,k,tau) +  cobalt%jldop(i,j,k)*dt*               &
                             grid_tmask(i,j,k)
       !
       ! Semilabile Dissolved Organic Nitrogen
       !
       cobalt%jsldon(i,j,k) = cobalt%jprod_sldon(i,j,k) - &
                              cobalt%gamma_sldon*cobalt%f_sldon(i,j,k)
       cobalt%p_sldon(i,j,k,tau) = cobalt%p_sldon(i,j,k,tau) +  cobalt%jsldon(i,j,k) * dt *               &
            grid_tmask(i,j,k)
       !
       ! Semilabile dissolved organic phosphorous  
       !
       cobalt%jsldop(i,j,k) = cobalt%jprod_sldop(i,j,k) - &
                              cobalt%gamma_sldop*cobalt%f_sldop(i,j,k)
       cobalt%p_sldop(i,j,k,tau) = cobalt%p_sldop(i,j,k,tau) + cobalt%jsldop(i,j,k) * dt *                &
                                  grid_tmask(i,j,k)
       !
       ! Refractory Dissolved Organic Nitrogen
       ! 
       cobalt%jsrdon(i,j,k) = cobalt%jprod_srdon(i,j,k) -  cobalt%gamma_srdon * cobalt%f_srdon(i,j,k)
       cobalt%p_srdon(i,j,k,tau) = cobalt%p_srdon(i,j,k,tau) +  cobalt%jsrdon(i,j,k) * dt *               &
            grid_tmask(i,j,k)
       !
       ! Refractory dissolved organic phosphorous
       !
       cobalt%jsrdop(i,j,k) = cobalt%jprod_srdop(i,j,k) - cobalt%gamma_srdop * cobalt%f_srdop(i,j,k)
       cobalt%p_srdop(i,j,k,tau) = cobalt%p_srdop(i,j,k,tau) + cobalt%jsrdop(i,j,k) * dt *                &
                                  grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     O2
    !
    do k = 1, nk ; do j =jsc, jec ; do i = isc, iec  !{
       cobalt%jo2(i,j,k) = (cobalt%o2_2_no3 * (phyto(DIAZO)%juptake_no3(i,j,k) +   &
            phyto(LARGE)%juptake_no3(i,j,k) + phyto(SMALL)%juptake_no3(i,j,k)) + & 
            cobalt%o2_2_nh4 *       &
            (phyto(DIAZO)%juptake_nh4(i,j,k) + phyto(LARGE)%juptake_nh4(i,j,k) +      &
            phyto(SMALL)%juptake_nh4(i,j,k) + &  
            phyto(DIAZO)%juptake_n2(i,j,k))) * grid_tmask(i,j,k)
       if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
          cobalt%jo2(i,j,k) = cobalt%jo2(i,j,k) - cobalt%o2_2_nh4*cobalt%jprod_nh4(i,j,k) &
                              - cobalt%o2_2_nitrif*cobalt%jnitrif(i,j,k) 
       endif  !}
       cobalt%p_o2(i,j,k,tau) = cobalt%p_o2(i,j,k,tau) + cobalt%jo2(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     The Carbon system
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Alkalinity
       !
       cobalt%jalk(i,j,k) = (2.0 * (cobalt%jdiss_cadet_arag(i,j,k) +         &
          cobalt%jdiss_cadet_calc(i,j,k) - cobalt%jprod_cadet_arag(i,j,k) - &
          cobalt%jprod_cadet_calc(i,j,k)) + phyto(DIAZO)%juptake_no3(i,j,k) + &
          phyto(LARGE)%juptake_no3(i,j,k) + phyto(SMALL)%juptake_no3(i,j,k) + &
          cobalt%jprod_nh4(i,j,k) - phyto(DIAZO)%juptake_nh4(i,j,k) - & 
          phyto(LARGE)%juptake_nh4(i,j,k) - phyto(SMALL)%juptake_nh4(i,j,k) -  &
          2.0 * cobalt%jnitrif(i,j,k) + cobalt%alk_2_n_denit * cobalt%jno3denit_wc(i,j,k))
       cobalt%p_alk(i,j,k,tau) = cobalt%p_alk(i,j,k,tau) + cobalt%jalk(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! Dissolved Inorganic Carbon
       !
       cobalt%jdic(i,j,k) =(cobalt%c_2_n * (cobalt%jno3(i,j,k) + &
          cobalt%jnh4(i,j,k) + cobalt%jno3denit_wc(i,j,k) - phyto(DIAZO)%juptake_n2(i,j,k)) + &
          cobalt%jdiss_cadet_arag(i,j,k) + cobalt%jdiss_cadet_calc(i,j,k) - &
          cobalt%jprod_cadet_arag(i,j,k) - cobalt%jprod_cadet_calc(i,j,k))
       cobalt%p_dic(i,j,k,tau) = cobalt%p_dic(i,j,k,tau) + cobalt%jdic(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     Lithogenic aluminosilicate particulates
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       cobalt%p_lith(i,j,k,tau) = cobalt%p_lith(i,j,k,tau) - cobalt%jlithdet(i,j,k) * dt *        &
            grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    call mpp_clock_end(id_clock_source_sink_loop6)
    call mpp_clock_begin(id_clock_cobalt_calc_diagnostics)
    !
    !Set the diagnostics tracer fields.
    !
    call g_tracer_set_values(tracer_list,'cased',  'field',cobalt%f_cased    ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'chl',    'field',cobalt%f_chl      ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion','field',cobalt%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'irr_mem' ,'field',cobalt%f_irr_mem ,isd,jsd,ntau=1)
    !
    !-----------------------------------------------------------------------
    !       Save variables for diagnostics
    !-----------------------------------------------------------------------
    !

    do j = jsc, jec ; do i = isc, iec  !{
      if (grid_kmt(i,j) .gt. 0) then !{
        cobalt%o2min(i,j)=cobalt%p_o2(i,j,1,tau)
        cobalt%z_o2min(i,j)=cobalt%zt(i,j,1)
        cobalt%z_sat_arag(i,j)=missing_value1
        cobalt%z_sat_calc(i,j)=missing_value1
        cobalt%mask_z_sat_arag(i,j) = .FALSE.
        cobalt%mask_z_sat_calc(i,j) = .FALSE.
        if (cobalt%omega_arag(i,j,1) .le. 1.0) cobalt%z_sat_arag(i,j)=0.0
        if (cobalt%omega_calc(i,j,1) .le. 1.0) cobalt%z_sat_calc(i,j)=0.0
      endif !}
    enddo ; enddo  !} i,j,k
    do j = jsc, jec ; do i = isc, iec  !{
    first = .true.
      do k = 2, nk
         if (k .le. grid_kmt(i,j) .and. first) then !{
           if (cobalt%p_o2(i,j,k,tau) .lt. cobalt%p_o2(i,j,k-1,tau)) then
             cobalt%o2min(i,j)=cobalt%p_o2(i,j,k,tau)
             cobalt%z_o2min(i,j)=cobalt%zt(i,j,k)
           else
             first = .false.
           endif !}
         endif !}
      enddo;
    enddo ; enddo  !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec  !{
      if (k .le. grid_kmt(i,j)) then !{
        if (cobalt%omega_arag(i,j,k) .le. 1.0 .and. cobalt%z_sat_arag(i,j) .lt. 0.0) then
          cobalt%z_sat_arag(i,j)=cobalt%zt(i,j,k)
          cobalt%mask_z_sat_arag(i,j) = .TRUE.
        endif
        if (cobalt%omega_calc(i,j,k) .le. 1.0 .and. cobalt%z_sat_calc(i,j) .lt. 0.0) then
          cobalt%z_sat_calc(i,j)=cobalt%zt(i,j,k)
          cobalt%mask_z_sat_calc(i,j) = .TRUE.
        endif
      endif !}
    enddo; enddo ; enddo  !} i,j,k

    !
    !---------------------------------------------------------------------
    ! Calculate total carbon  = Dissolved Inorganic Carbon + Phytoplankton Carbon
    !   + Dissolved Organic Carbon (including refractory) + Heterotrophic Biomass
    !   + Detrital Orgainc and Inorganic Carbon
    ! For the oceanic carbon budget, a constant 42 uM of dissolved organic
    ! carbon is added to represent the refractory component.
    ! For the oceanic nitrogen budget, a constant 2 uM of dissolved organic
    ! nitrogen is added to represent the refractory component.
    !---------------------------------------------------------------------
    !
    cobalt%tot_layer_int_c(:,:,:) = (cobalt%p_dic(:,:,:,tau) + 4.2e-5 + cobalt%p_cadet_arag(:,:,:,tau) +&
         cobalt%p_cadet_calc(:,:,:,tau) + cobalt%c_2_n * (cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) +      &
         cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) + &
         cobalt%p_ldon(:,:,:,tau) + cobalt%p_sldon(:,:,:,tau) + cobalt%p_srdon(:,:,:,tau) +  &
         cobalt%p_ndet(:,:,:,tau) + cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + &
         cobalt%p_nlgz(:,:,:,tau))) * rho_dzt(:,:,:)

    cobalt%tot_layer_int_fe(:,:,:) = (cobalt%p_fed(:,:,:,tau) + cobalt%p_fedi(:,:,:,tau) + &
         cobalt%p_felg(:,:,:,tau) + cobalt%p_fesm(:,:,:,tau) + & 
         cobalt%p_fedet(:,:,:,tau)) * rho_dzt(:,:,:) 

    cobalt%tot_layer_int_n(:,:,:) = (cobalt%p_no3(:,:,:,tau) + &
         cobalt%p_nh4(:,:,:,tau) + cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) + &
         cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) + &
         cobalt%p_ldon(:,:,:,tau) + cobalt%p_sldon(:,:,:,tau) + cobalt%p_srdon(:,:,:,tau) +  cobalt%p_ndet(:,:,:,tau) + &
         cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + cobalt%p_nlgz(:,:,:,tau)) * & 
         rho_dzt(:,:,:)

    cobalt%tot_layer_int_p(:,:,:) = (cobalt%p_po4(:,:,:,tau) + &
         cobalt%p_ndi(:,:,:,tau)*phyto(1)%p_2_n_static + &
         cobalt%p_nlg(:,:,:,tau)*phyto(2)%p_2_n_static + &
         cobalt%p_nsm(:,:,:,tau)*phyto(3)%p_2_n_static + &
         cobalt%p_ldop(:,:,:,tau) + cobalt%p_sldop(:,:,:,tau) + &
         cobalt%p_srdop(:,:,:,tau) + cobalt%p_pdet(:,:,:,tau) + &
         bact(1)%q_p_2_n*cobalt%p_nbact(:,:,:,tau) + zoo(1)%q_p_2_n*cobalt%p_nsmz(:,:,:,tau) +  &
         zoo(2)%q_p_2_n*cobalt%p_nmdz(:,:,:,tau) + zoo(3)%q_p_2_n*cobalt%p_nlgz(:,:,:,tau))  &
         * rho_dzt(:,:,:)

    cobalt%tot_layer_int_si(:,:,:) = (cobalt%p_sio4(:,:,:,tau) + cobalt%p_silg(:,:,:,tau) +   &
         cobalt%p_sidet(:,:,:,tau)) * rho_dzt(:,:,:)

!****************************************************************************************************

    allocate(rho_dzt_100(isc:iec,jsc:jec))
    !
    !---------------------------------------------------------------------
    ! calculate upper 100 m vertical integrals
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec !{
       rho_dzt_100(i,j) = rho_dzt(i,j,1)
       do n = 1, NUM_PHYTO  !{
          phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jprod_n_new_100(i,j) = phyto(n)%juptake_no3(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jzloss_n_100(i,j) = phyto(n)%jzloss_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jexuloss_n_100(i,j) = phyto(n)%jexuloss_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%f_n_100(i,j) = phyto(n)%f_n(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n
       phyto(DIAZO)%jprod_n_n2_100(i,j) = phyto(DIAZO)%juptake_n2(i,j,1) * rho_dzt(i,j,1)
       phyto(SMALL)%jvirloss_n_100(i,j) = phyto(SMALL)%jvirloss_n(i,j,1) * rho_dzt(i,j,1)
       phyto(SMALL)%jaggloss_n_100(i,j) = phyto(SMALL)%jaggloss_n(i,j,1) * rho_dzt(i,j,1)
       phyto(LARGE)%jaggloss_n_100(i,j) = phyto(LARGE)%jaggloss_n(i,j,1) * rho_dzt(i,j,1)

       do n = 1, NUM_ZOO  !{
          zoo(n)%jprod_n_100(i,j) = zoo(n)%jprod_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jingest_n_100(i,j) = zoo(n)%jingest_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jremin_n_100(i,j) = zoo(n)%jprod_nh4(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%f_n_100(i,j) = zoo(n)%f_n(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n

       do n = 1,2  !{
          zoo(n)%jzloss_n_100(i,j) = zoo(n)%jzloss_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jprod_don_100(i,j) = (zoo(n)%jprod_ldon(i,j,1) + zoo(n)%jprod_sldon(i,j,1) + &
             zoo(n)%jprod_srdon(i,j,1))  * rho_dzt(i,j,1)
       enddo   !} n

       do n = 2,3  !{
          zoo(n)%jhploss_n_100(i,j) = zoo(n)%jhploss_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jprod_ndet_100(i,j) = zoo(n)%jprod_ndet(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n

       cobalt%hp_jingest_n_100(i,j) = cobalt%hp_jingest_n(i,j,1)*rho_dzt(i,j,1)
       cobalt%hp_jremin_n_100(i,j) =  cobalt%hp_jingest_n(i,j,1)*rho_dzt(i,j,1)*cobalt%hp_phi_nh4
       cobalt%hp_jprod_ndet_100(i,j) =  cobalt%hp_jingest_n(i,j,1)*rho_dzt(i,j,1)*cobalt%hp_phi_det

       bact(1)%jprod_n_100(i,j) = bact(1)%jprod_n(i,j,1) * rho_dzt(i,j,1)
       bact(1)%jzloss_n_100(i,j) = bact(1)%jzloss_n(i,j,1) * rho_dzt(i,j,1)
       bact(1)%jvirloss_n_100(i,j) = bact(1)%jvirloss_n(i,j,1) * rho_dzt(i,j,1)
       bact(1)%jremin_n_100(i,j) = bact(1)%jprod_nh4(i,j,1) * rho_dzt(i,j,1)
       bact(1)%juptake_ldon_100(i,j) = bact(1)%juptake_ldon(i,j,1) * rho_dzt(i,j,1)
       bact(1)%f_n_100(i,j) = bact(1)%f_n(i,j,1) * rho_dzt(i,j,1)

       cobalt%jprod_lithdet_100(i,j) = cobalt%jprod_lithdet(i,j,1) * rho_dzt(i,j,1)
       cobalt%jprod_sidet_100(i,j) = cobalt%jprod_sidet(i,j,1) * rho_dzt(i,j,1)
       cobalt%jprod_cadet_calc_100(i,j) = cobalt%jprod_cadet_calc(i,j,1) * rho_dzt(i,j,1)
       cobalt%jprod_cadet_arag_100(i,j) = cobalt%jprod_cadet_arag(i,j,1) * rho_dzt(i,j,1)
       cobalt%jremin_ndet_100(i,j) = cobalt%jremin_ndet(i,j,1) * rho_dzt(i,j,1)

       cobalt%f_ndet_100(i,j) = cobalt%f_ndet(i,j,1)*rho_dzt(i,j,1)
       cobalt%f_don_100(i,j) = (cobalt%f_ldon(i,j,1)+cobalt%f_sldon(i,j,1)+cobalt%f_srdon(i,j,1))* &
           rho_dzt(i,j,1)
       cobalt%f_silg_100(i,j) = cobalt%f_silg(i,j,1)*rho_dzt(i,j,1)

       cobalt%fndet_100(i,j) = cobalt%f_ndet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fpdet_100(i,j) = cobalt%f_pdet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%ffedet_100(i,j) = cobalt%f_fedet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%flithdet_100(i,j) = cobalt%f_lithdet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fsidet_100(i,j) = cobalt%f_sidet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fcadet_arag_100(i,j) = cobalt%f_cadet_arag(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fcadet_calc_100(i,j) = cobalt%f_cadet_calc(i,j,1) * cobalt%Rho_0 * cobalt%wsink
    enddo; enddo !} i,j

    do j = jsc, jec ; do i = isc, iec ; !{
       k_100 = 1
       do k = 2, grid_kmt(i,j)  !{
          if (rho_dzt_100(i,j) .lt. cobalt%Rho_0 * 100.0) then 
             k_100 = k
             rho_dzt_100(i,j) = rho_dzt_100(i,j) + rho_dzt(i,j,k)
             do n = 1, NUM_PHYTO !{
                phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n_100(i,j) + phyto(n)%jprod_n(i,j,k)* & 
                   rho_dzt(i,j,k)
                phyto(n)%jprod_n_new_100(i,j) = phyto(n)%jprod_n_new_100(i,j) + phyto(n)%juptake_no3(i,j,k)* &
                   rho_dzt(i,j,k)
                phyto(n)%jzloss_n_100(i,j) = phyto(n)%jzloss_n_100(i,j) + phyto(n)%jzloss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                phyto(n)%jexuloss_n_100(i,j) = phyto(n)%jexuloss_n_100(i,j) + phyto(n)%jexuloss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                phyto(n)%f_n_100(i,j) = phyto(n)%f_n_100(i,j) + phyto(n)%f_n(i,j,k)*rho_dzt(i,j,k) 
             enddo !} n
             phyto(DIAZO)%jprod_n_n2_100(i,j) = phyto(DIAZO)%jprod_n_n2_100(i,j) + &
                 phyto(DIAZO)%juptake_n2(i,j,k)*rho_dzt(i,j,k)
             phyto(SMALL)%jvirloss_n_100(i,j) = phyto(SMALL)%jvirloss_n_100(i,j) + &
                 phyto(SMALL)%jvirloss_n(i,j,k)*rho_dzt(i,j,k)
             phyto(SMALL)%jaggloss_n_100(i,j) = phyto(SMALL)%jaggloss_n_100(i,j) + &
                 phyto(SMALL)%jaggloss_n(i,j,k)*rho_dzt(i,j,k)
             phyto(LARGE)%jaggloss_n_100(i,j) = phyto(LARGE)%jaggloss_n_100(i,j) + &
                 phyto(LARGE)%jaggloss_n(i,j,k)*rho_dzt(i,j,k)

             do n = 1, NUM_ZOO !{
                zoo(n)%jprod_n_100(i,j) = zoo(n)%jprod_n_100(i,j) + zoo(n)%jprod_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jingest_n_100(i,j) = zoo(n)%jingest_n_100(i,j) + zoo(n)%jingest_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jremin_n_100(i,j) = zoo(n)%jremin_n_100(i,j) + zoo(n)%jprod_nh4(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%f_n_100(i,j) = zoo(n)%f_n_100(i,j) + zoo(n)%f_n(i,j,k)*rho_dzt(i,j,k)
             enddo !} n

             do n = 1,2 !{
                zoo(n)%jzloss_n_100(i,j) = zoo(n)%jzloss_n_100(i,j) + zoo(n)%jzloss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jprod_don_100(i,j) = zoo(n)%jprod_don_100(i,j) + (zoo(n)%jprod_ldon(i,j,k) + &
                   zoo(n)%jprod_sldon(i,j,k) + zoo(n)%jprod_srdon(i,j,k))*rho_dzt(i,j,k)
             enddo !} n

             do n = 2,3 !{
                zoo(n)%jhploss_n_100(i,j) = zoo(n)%jhploss_n_100(i,j) + zoo(n)%jhploss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jprod_ndet_100(i,j) = zoo(n)%jprod_ndet_100(i,j) + zoo(n)%jprod_ndet(i,j,k)* &
                   rho_dzt(i,j,k)
             enddo !} n

             cobalt%hp_jingest_n_100(i,j) = cobalt%hp_jingest_n_100(i,j) + cobalt%hp_jingest_n(i,j,k)* &
                 rho_dzt(i,j,k)
             cobalt%hp_jremin_n_100(i,j) = cobalt%hp_jremin_n_100(i,j) + cobalt%hp_jingest_n(i,j,k)* &
                 cobalt%hp_phi_nh4*rho_dzt(i,j,k)
             cobalt%hp_jprod_ndet_100(i,j) = cobalt%hp_jprod_ndet_100(i,j) + cobalt%hp_jingest_n(i,j,k)* &
                 cobalt%hp_phi_det*rho_dzt(i,j,k)

             bact(1)%jprod_n_100(i,j) = bact(1)%jprod_n_100(i,j) + bact(1)%jprod_n(i,j,k) * rho_dzt(i,j,k)
             bact(1)%jzloss_n_100(i,j) = bact(1)%jzloss_n_100(i,j) + bact(1)%jzloss_n(i,j,k) * rho_dzt(i,j,k)
             bact(1)%jvirloss_n_100(i,j) = bact(1)%jvirloss_n_100(i,j) + bact(1)%jvirloss_n(i,j,k) * rho_dzt(i,j,k)
             bact(1)%jremin_n_100(i,j) = bact(1)%jremin_n_100(i,j) + bact(1)%jprod_nh4(i,j,k) * rho_dzt(i,j,k)
             bact(1)%juptake_ldon_100(i,j) = bact(1)%juptake_ldon_100(i,j) + bact(1)%juptake_ldon(i,j,k) * rho_dzt(i,j,k)
             bact(1)%f_n_100(i,j) = bact(1)%f_n_100(i,j) + bact(1)%f_n(i,j,k)*rho_dzt(i,j,k)

             cobalt%jprod_lithdet_100(i,j) = cobalt%jprod_lithdet_100(i,j) + cobalt%jprod_lithdet(i,j,k) * rho_dzt(i,j,k)
             cobalt%jprod_sidet_100(i,j) = cobalt%jprod_sidet_100(i,j) + cobalt%jprod_sidet(i,j,k) * rho_dzt(i,j,k)
             cobalt%jprod_cadet_calc_100(i,j) = cobalt%jprod_cadet_calc_100(i,j) + cobalt%jprod_cadet_calc(i,j,k) * rho_dzt(i,j,k)
             cobalt%jprod_cadet_arag_100(i,j) = cobalt%jprod_cadet_arag_100(i,j) + cobalt%jprod_cadet_arag(i,j,k) * rho_dzt(i,j,k)
             cobalt%jremin_ndet_100(i,j) = cobalt%jremin_ndet_100(i,j) + cobalt%jremin_ndet(i,j,k) * rho_dzt(i,j,k)
             cobalt%f_ndet_100(i,j) = cobalt%f_ndet_100(i,j) + cobalt%f_ndet(i,j,k)*rho_dzt(i,j,k)
             cobalt%f_don_100(i,j) = cobalt%f_don_100(i,j) + (cobalt%f_ldon(i,j,k) + cobalt%f_sldon(i,j,k) + &
                cobalt%f_srdon(i,j,k))*rho_dzt(i,j,k)
             cobalt%f_silg_100(i,j) = cobalt%f_silg_100(i,j) + cobalt%f_silg(i,j,k)*rho_dzt(i,j,k)

             cobalt%fndet_100(i,j) = cobalt%f_ndet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fpdet_100(i,j) = cobalt%f_pdet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%ffedet_100(i,j) = cobalt%f_fedet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%flithdet_100(i,j) = cobalt%f_lithdet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fsidet_100(i,j) = cobalt%f_sidet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fcadet_arag_100(i,j) = cobalt%f_cadet_arag(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fcadet_calc_100(i,j) = cobalt%f_cadet_calc(i,j,k) * cobalt%Rho_0 * cobalt%wsink

          endif
       enddo  !} k

       if (k_100 .gt. 1 .and. k_100 .lt. grid_kmt(i,j)) then
          drho_dzt = cobalt%Rho_0 * 100.0 - rho_dzt_100(i,j)
          do n = 1, NUM_PHYTO !{
              phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n_100(i,j) + phyto(n)%jprod_n(i,j,k_100)* &
                 drho_dzt
              phyto(n)%jprod_n_new_100(i,j) = phyto(n)%jprod_n_new_100(i,j) + phyto(n)%juptake_no3(i,j,k_100)* &
                 drho_dzt
              phyto(n)%jzloss_n_100(i,j) = phyto(n)%jzloss_n_100(i,j) + phyto(n)%jzloss_n(i,j,k_100)* &
                 drho_dzt
             phyto(n)%jexuloss_n_100(i,j) = phyto(n)%jexuloss_n_100(i,j) + phyto(n)%jexuloss_n(i,j,k_100)* &
                 drho_dzt
              phyto(n)%f_n_100(i,j) = phyto(n)%f_n_100(i,j) + phyto(n)%f_n(i,j,k_100)*drho_dzt
           enddo !} n
           phyto(DIAZO)%jprod_n_n2_100(i,j) = phyto(DIAZO)%jprod_n_n2_100(i,j) + &
               phyto(DIAZO)%juptake_n2(i,j,k_100)*drho_dzt
           phyto(SMALL)%jvirloss_n_100(i,j) = phyto(SMALL)%jvirloss_n_100(i,j) + &
               phyto(SMALL)%jvirloss_n(i,j,k_100)*drho_dzt
           phyto(SMALL)%jaggloss_n_100(i,j) = phyto(SMALL)%jaggloss_n_100(i,j) + &
               phyto(SMALL)%jaggloss_n(i,j,k_100)*drho_dzt
           phyto(LARGE)%jaggloss_n_100(i,j) = phyto(LARGE)%jaggloss_n_100(i,j) + &
               phyto(LARGE)%jaggloss_n(i,j,k_100)*drho_dzt

           do n = 1, NUM_ZOO !{
               zoo(n)%jprod_n_100(i,j) = zoo(n)%jprod_n_100(i,j) + zoo(n)%jprod_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jingest_n_100(i,j) = zoo(n)%jingest_n_100(i,j) + zoo(n)%jingest_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jremin_n_100(i,j) = zoo(n)%jremin_n_100(i,j) + zoo(n)%jprod_nh4(i,j,k_100)* &
                 drho_dzt
               zoo(n)%f_n_100(i,j) = zoo(n)%f_n_100(i,j) + zoo(n)%f_n(i,j,k_100)*drho_dzt
           enddo !} n

           do n = 1,2 !{
               zoo(n)%jzloss_n_100(i,j) = zoo(n)%jzloss_n_100(i,j) + zoo(n)%jzloss_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jprod_don_100(i,j) = zoo(n)%jprod_don_100(i,j) + (zoo(n)%jprod_ldon(i,j,k_100) + &
                 zoo(n)%jprod_sldon(i,j,k_100) + zoo(n)%jprod_srdon(i,j,k_100))*drho_dzt
           enddo !} n

           do n = 2,3 !{
               zoo(n)%jhploss_n_100(i,j) = zoo(n)%jhploss_n_100(i,j) + zoo(n)%jhploss_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jprod_ndet_100(i,j) = zoo(n)%jprod_ndet_100(i,j) + zoo(n)%jprod_ndet(i,j,k_100)* &
                 drho_dzt
           enddo !} n

           cobalt%hp_jingest_n_100(i,j) = cobalt%hp_jingest_n_100(i,j) + cobalt%hp_jingest_n(i,j,k_100)* &
               drho_dzt
           cobalt%hp_jremin_n_100(i,j) = cobalt%hp_jremin_n_100(i,j) + cobalt%hp_jingest_n(i,j,k_100)* &
               cobalt%hp_phi_nh4*drho_dzt
           cobalt%hp_jprod_ndet_100(i,j) = cobalt%hp_jprod_ndet_100(i,j) + cobalt%hp_jingest_n(i,j,k_100)* &
               cobalt%hp_phi_det*drho_dzt

           bact(1)%jprod_n_100(i,j) = bact(1)%jprod_n_100(i,j) + bact(1)%jprod_n(i,j,k_100)* &
                drho_dzt
           bact(1)%jzloss_n_100(i,j) = bact(1)%jzloss_n_100(i,j) + bact(1)%jzloss_n(i,j,k_100)* & 
                drho_dzt
           bact(1)%jvirloss_n_100(i,j) = bact(1)%jvirloss_n_100(i,j) + bact(1)%jvirloss_n(i,j,k_100)* & 
                drho_dzt
           bact(1)%jremin_n_100(i,j) = bact(1)%jremin_n_100(i,j) + bact(1)%jprod_nh4(i,j,k_100)* & 
                drho_dzt
           bact(1)%juptake_ldon_100(i,j) = bact(1)%juptake_ldon_100(i,j) + bact(1)%juptake_ldon(i,j,k_100)* &
                drho_dzt
           bact(1)%f_n_100(i,j) = bact(1)%f_n_100(i,j) + bact(1)%f_n(i,j,k_100)*drho_dzt

           cobalt%jprod_lithdet_100(i,j) = cobalt%jprod_lithdet_100(i,j) + cobalt%jprod_lithdet(i,j,k_100)* &
                drho_dzt
           cobalt%jprod_sidet_100(i,j) = cobalt%jprod_sidet_100(i,j) + cobalt%jprod_sidet(i,j,k_100)* &
                drho_dzt
           cobalt%jprod_cadet_calc_100(i,j) = cobalt%jprod_cadet_calc_100(i,j) + cobalt%jprod_cadet_calc(i,j,k_100)* &
                drho_dzt
           cobalt%jprod_cadet_arag_100(i,j) = cobalt%jprod_cadet_arag_100(i,j) + cobalt%jprod_cadet_arag(i,j,k_100)* &
                drho_dzt
           cobalt%jremin_ndet_100(i,j) = cobalt%jremin_ndet_100(i,j) + cobalt%jremin_ndet(i,j,k_100)* &
                drho_dzt

           cobalt%f_ndet_100(i,j) = cobalt%f_ndet_100(i,j) + cobalt%f_ndet(i,j,k_100)*drho_dzt
           cobalt%f_don_100(i,j) = cobalt%f_don_100(i,j) + (cobalt%f_ldon(i,j,k_100) + cobalt%f_sldon(i,j,k_100) + &
              cobalt%f_srdon(i,j,k_100))*drho_dzt
           cobalt%f_silg_100(i,j) = cobalt%f_silg_100(i,j) + cobalt%f_silg(i,j,k_100)*drho_dzt

           cobalt%fndet_100(i,j) = cobalt%f_ndet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fpdet_100(i,j) = cobalt%f_pdet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%ffedet_100(i,j) = cobalt%f_fedet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%flithdet_100(i,j) = cobalt%f_lithdet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fsidet_100(i,j) = cobalt%f_sidet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fcadet_arag_100(i,j) = cobalt%f_cadet_arag(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fcadet_calc_100(i,j) = cobalt%f_cadet_calc(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
       endif

       cobalt%jprod_allphytos_100(i,j) = phyto(SMALL)%jprod_n_100(i,j) + phyto(LARGE)%jprod_n_100(i,j) + &
          phyto(DIAZO)%jprod_n_100(i,j) 
    enddo ; enddo  !} i,j
    deallocate(rho_dzt_100)

    do j = jsc, jec ; do i = isc, iec ; !{
      cobalt%btm_temp(i,j) = TEMP(i,j,grid_kmt(i,j))
      cobalt%btm_o2(i,j) = cobalt%f_o2(i,j,grid_kmt(i,j))      
    enddo; enddo  !} i, j

    !
    !---------------------------------------------------------------------
    ! calculate upper 200m vertical integrals for mesozooplankton
    ! quantities for comparison with COPEPOD database
    !---------------------------------------------------------------------
    !
    allocate(rho_dzt_200(isc:iec,jsc:jec))
    do j = jsc, jec ; do i = isc, iec !{
       rho_dzt_200(i,j) = rho_dzt(i,j,1)
       cobalt%jprod_mesozoo_200(i,j) = (zoo(2)%jprod_n(i,j,1) + zoo(3)%jprod_n(i,j,1))*rho_dzt(i,j,1)
       cobalt%f_mesozoo_200(i,j) = (zoo(2)%f_n(i,j,1)+zoo(3)%f_n(i,j,1))*rho_dzt(i,j,1)
    enddo; enddo !} i,j

    do j = jsc, jec ; do i = isc, iec ; !{
       k_200 = 1
       do k = 2, grid_kmt(i,j)  !{
          if (rho_dzt_200(i,j) .lt. cobalt%Rho_0 * 200.0) then
             k_200 = k
             rho_dzt_200(i,j) = rho_dzt_200(i,j) + rho_dzt(i,j,k)
             cobalt%jprod_mesozoo_200(i,j) = cobalt%jprod_mesozoo_200(i,j) + &
                (zoo(2)%jprod_n(i,j,k) + zoo(3)%jprod_n(i,j,k))*rho_dzt(i,j,k)
             cobalt%f_mesozoo_200(i,j) = cobalt%f_mesozoo_200(i,j) + &
                (zoo(2)%f_n(i,j,k)+zoo(3)%f_n(i,j,k))*rho_dzt(i,j,k)
          endif
       enddo  !} k

       if (k_200 .gt. 1 .and. k_200 .lt. grid_kmt(i,j)) then
          drho_dzt = cobalt%Rho_0 * 200.0 - rho_dzt_200(i,j)
          cobalt%jprod_mesozoo_200(i,j) = cobalt%jprod_mesozoo_200(i,j) + &
              (zoo(2)%jprod_n(i,j,k_200) + zoo(3)%jprod_n(i,j,k_200))*drho_dzt
          cobalt%f_mesozoo_200(i,j) = cobalt%f_mesozoo_200(i,j) + &
              (zoo(2)%f_n(i,j,k_200)+zoo(3)%f_n(i,j,k_200))*drho_dzt
       endif
    enddo ; enddo  !} i,j

    call g_tracer_get_pointer(tracer_list,'alk','runoff_tracer_flux',cobalt%runoff_flux_alk)
    call g_tracer_get_pointer(tracer_list,'dic','runoff_tracer_flux',cobalt%runoff_flux_dic)
    call g_tracer_get_pointer(tracer_list,'fed','runoff_tracer_flux',cobalt%runoff_flux_fed)
    call g_tracer_get_pointer(tracer_list,'fed','drydep',cobalt%dry_fed)
    call g_tracer_get_pointer(tracer_list,'fed','wetdep',cobalt%wet_fed)
    call g_tracer_get_pointer(tracer_list,'lith','drydep',cobalt%dry_lith)
    call g_tracer_get_pointer(tracer_list,'lith','wetdep',cobalt%wet_lith)
    call g_tracer_get_pointer(tracer_list,'lith','runoff_tracer_flux',cobalt%runoff_flux_lith)
    call g_tracer_get_pointer(tracer_list,'no3','runoff_tracer_flux',cobalt%runoff_flux_no3)
    call g_tracer_get_pointer(tracer_list,'no3','drydep',cobalt%dry_no3)
    call g_tracer_get_pointer(tracer_list,'no3','wetdep',cobalt%wet_no3)
    call g_tracer_get_pointer(tracer_list,'nh4','drydep',cobalt%dry_nh4)
    call g_tracer_get_pointer(tracer_list,'nh4','wetdep',cobalt%wet_nh4)
    call g_tracer_get_pointer(tracer_list,'po4','drydep',cobalt%dry_po4)
    call g_tracer_get_pointer(tracer_list,'po4','wetdep',cobalt%wet_po4)
    call g_tracer_get_pointer(tracer_list,'ldon','runoff_tracer_flux',cobalt%runoff_flux_ldon)
    call g_tracer_get_pointer(tracer_list,'sldon','runoff_tracer_flux',cobalt%runoff_flux_sldon)
    call g_tracer_get_pointer(tracer_list,'srdon','runoff_tracer_flux',cobalt%runoff_flux_srdon)
    call g_tracer_get_pointer(tracer_list,'ndet','runoff_tracer_flux',cobalt%runoff_flux_ndet)
    call g_tracer_get_pointer(tracer_list,'po4','runoff_tracer_flux',cobalt%runoff_flux_po4)
    call g_tracer_get_pointer(tracer_list,'ldop','runoff_tracer_flux',cobalt%runoff_flux_ldop)
    call g_tracer_get_pointer(tracer_list,'sldop','runoff_tracer_flux',cobalt%runoff_flux_sldop)
    call g_tracer_get_pointer(tracer_list,'srdop','runoff_tracer_flux',cobalt%runoff_flux_srdop)


!---------------------------------------------------------------------
! Add vertical integrals for diagnostics
!---------------------------------------------------------------------
!

    call mpp_clock_end(id_clock_cobalt_calc_diagnostics)
    call mpp_clock_begin(id_clock_cobalt_send_diagnostics)

!---------------------------------------------------------------------
!
! Send phytoplankton diagnostic data

    do n= 1, NUM_PHYTO
       if (phyto(n)%id_def_fe .gt. 0)          &
            used = send_data(phyto(n)%id_def_fe,     phyto(n)%def_fe,           &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_felim .gt. 0)           &
            used = send_data(phyto(n)%id_felim,      phyto(n)%felim,            &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_irrlim .gt. 0)          &
            used = send_data(phyto(n)%id_irrlim,     phyto(n)%irrlim,           &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jzloss_n .gt. 0)          &
            used = send_data(phyto(n)%id_jzloss_n, phyto(n)%jzloss_n*rho_dzt,      &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jaggloss_n .gt. 0)          &
            used = send_data(phyto(n)%id_jaggloss_n, phyto(n)%jaggloss_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jvirloss_n .gt. 0)          &
            used = send_data(phyto(n)%id_jvirloss_n, phyto(n)%jvirloss_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jexuloss_n .gt. 0)          &
            used = send_data(phyto(n)%id_jexuloss_n, phyto(n)%jexuloss_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
!       if (phyto(n)%id_jhploss_n .gt. 0)          &
!            used = send_data(phyto(n)%id_jhploss_n, phyto(n)%jhploss_n*rho_dzt,     &
!            model_time, rmask = grid_tmask,&
!            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_fe .gt. 0)          &
            used = send_data(phyto(n)%id_juptake_fe, phyto(n)%juptake_fe*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_nh4 .gt. 0)          &
            used = send_data(phyto(n)%id_juptake_nh4, phyto(n)%juptake_nh4*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_no3 .gt. 0)          &
            used = send_data(phyto(n)%id_juptake_no3, phyto(n)%juptake_no3*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_po4 .gt. 0)          &
            used = send_data(phyto(n)%id_juptake_po4, phyto(n)%juptake_po4*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_sio4 .gt. 0)          &
            used = send_data(phyto(n)%id_juptake_sio4, phyto(n)%juptake_sio4*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_n2 .gt. 0)          &
            used = send_data(phyto(n)%id_juptake_n2, phyto(n)%juptake_n2*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_n .gt. 0)          &
            used = send_data(phyto(n)%id_jprod_n, phyto(n)%jprod_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
!       if (phyto(n)%id_liebig_lim .gt. 0)      &
!            used = send_data(phyto(n)%id_liebig_lim,phyto(n)%liebig_lim,          &
!            model_time, rmask = grid_tmask,& 
!            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_mu .gt. 0)              &
            used = send_data(phyto(n)%id_mu,        phyto(n)%mu,                  &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_nh4lim .gt. 0)          &
            used = send_data(phyto(n)%id_nh4lim,     phyto(n)%nh4lim,             &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_no3lim .gt. 0)          &
            used = send_data(phyto(n)%id_no3lim,     phyto(n)%no3lim,             &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_po4lim .gt. 0)          &
            used = send_data(phyto(n)%id_po4lim,     phyto(n)%po4lim,             &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_o2lim .gt. 0)          &
            used = send_data(phyto(n)%id_o2lim,     phyto(n)%o2lim,             &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_fe_2_n .gt. 0)        &
            used = send_data(phyto(n)%id_q_fe_2_n,   phyto(n)%q_fe_2_n,           &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_silim .gt. 0)     &
            used = send_data(phyto(n)%id_silim, phyto(n)%silim,       &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_si_2_n .gt. 0)     &
            used = send_data(phyto(n)%id_q_si_2_n, phyto(n)%q_si_2_n,       &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_theta .gt. 0)           &
            used = send_data(phyto(n)%id_theta,      phyto(n)%theta,              &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo
    !--------------------------------------------------------------------------------------
    ! Send bacterial diagnostic data
    !

    if (bact(1)%id_jzloss_n .gt. 0)          &
       used = send_data(bact(1)%id_jzloss_n, bact(1)%jzloss_n*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
!    if (bact(1)%id_jhploss_n .gt. 0)          &
!       used = send_data(bact(1)%id_jhploss_n, bact(1)%jhploss_n*rho_dzt,           &
!       model_time, rmask = grid_tmask,&
!       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jvirloss_n .gt. 0)          &
       used = send_data(bact(1)%id_jvirloss_n, bact(1)%jvirloss_n*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_juptake_ldon .gt. 0)          &
       used = send_data(bact(1)%id_juptake_ldon, bact(1)%juptake_ldon*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_juptake_ldop .gt. 0)          &
       used = send_data(bact(1)%id_juptake_ldop, bact(1)%juptake_ldop*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jprod_nh4 .gt. 0)          &
       used = send_data(bact(1)%id_jprod_nh4, bact(1)%jprod_nh4*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jprod_po4 .gt. 0)          &
       used = send_data(bact(1)%id_jprod_po4, bact(1)%jprod_po4*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jprod_n .gt. 0)          &
       used = send_data(bact(1)%id_jprod_n, bact(1)%jprod_n*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_temp_lim .gt. 0)          &
       used = send_data(bact(1)%id_temp_lim, bact(1)%temp_lim,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    !--------------------------------------------------------------------------------------
    ! Send zooplankton diagnostic data
    !

    do n= 1, NUM_ZOO
       if (zoo(n)%id_jzloss_n .gt. 0)          &
            used = send_data(zoo(n)%id_jzloss_n, zoo(n)%jzloss_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jhploss_n .gt. 0)          &
            used = send_data(zoo(n)%id_jhploss_n, zoo(n)%jhploss_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_n .gt. 0)          &
            used = send_data(zoo(n)%id_jingest_n, zoo(n)%jingest_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_p .gt. 0)          &
            used = send_data(zoo(n)%id_jingest_p, zoo(n)%jingest_p*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_sio2 .gt. 0)          &
            used = send_data(zoo(n)%id_jingest_sio2, zoo(n)%jingest_sio2*rho_dzt,      &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_fe .gt. 0)          &
            used = send_data(zoo(n)%id_jingest_fe, zoo(n)%jingest_fe*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_ndet .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_ndet, zoo(n)%jprod_ndet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_pdet .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_pdet, zoo(n)%jprod_pdet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_ldon .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_ldon, zoo(n)%jprod_ldon*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_ldop .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_ldop, zoo(n)%jprod_ldop*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sldon .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_sldon, zoo(n)%jprod_sldon*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sldop .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_sldop, zoo(n)%jprod_sldop*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (zoo(n)%id_jprod_srdon .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_srdon, zoo(n)%jprod_srdon*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_srdop .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_srdop, zoo(n)%jprod_srdop*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_fedet .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_fedet, zoo(n)%jprod_fedet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_fed .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_fed,  zoo(n)%jprod_fed*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sidet .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_sidet, zoo(n)%jprod_sidet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sio4 .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_sio4, zoo(n)%jprod_sio4*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_po4 .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_po4,  zoo(n)%jprod_po4*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_nh4 .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_nh4,  zoo(n)%jprod_nh4*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_n .gt. 0)          &
            used = send_data(zoo(n)%id_jprod_n,   zoo(n)%jprod_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_temp_lim .gt. 0)          &
            used = send_data(zoo(n)%id_temp_lim, zoo(n)%temp_lim,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo
    !
    ! Production diagnostics
    !
    if (cobalt%id_jprod_cadet_arag .gt. 0)    &
       used = send_data(cobalt%id_jprod_cadet_arag, cobalt%jprod_cadet_arag * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_cadet_calc .gt. 0)    &
       used = send_data(cobalt%id_jprod_cadet_calc, cobalt%jprod_cadet_calc * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_ndet .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_ndet, cobalt%jprod_ndet*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_pdet .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_pdet, cobalt%jprod_pdet*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_srdon .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_srdon, cobalt%jprod_srdon*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_sldon .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_sldon, cobalt%jprod_sldon*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_ldon .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_ldon, cobalt%jprod_ldon*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_srdop .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_srdop, cobalt%jprod_srdop*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_sldop .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_sldop, cobalt%jprod_sldop*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_ldop .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_ldop, cobalt%jprod_ldop*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_nh4 .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_nh4, cobalt%jprod_nh4*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_po4 .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_po4, cobalt%jprod_po4*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_fedet .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_fedet,  cobalt%jprod_fedet*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_fed .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_fed, cobalt%jprod_fed*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_sidet .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_sidet, cobalt%jprod_sidet*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_sio4 .gt. 0)          &
    !    used = send_data(cobalt%id_jprod_sio4, cobalt%jprod_sio4*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_lithdet .gt. 0)          &
        used = send_data(cobalt%id_jprod_lithdet, cobalt%jprod_lithdet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_det_jzloss_n .gt. 0)          &
    !    used = send_data(cobalt%id_det_jzloss_n, cobalt%det_jzloss_n*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_det_jhploss_n .gt. 0)          &
    !    used = send_data(cobalt%id_det_jhploss_n, cobalt%det_jhploss_n*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdiss_cadet_arag .gt. 0)          &
        used = send_data(cobalt%id_jdiss_cadet_arag, cobalt%jdiss_cadet_arag*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdiss_cadet_calc .gt. 0)          &
        used = send_data(cobalt%id_jdiss_cadet_calc, cobalt%jdiss_cadet_calc*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdiss_sidet .gt. 0)          &
        used = send_data(cobalt%id_jdiss_sidet, cobalt%jdiss_sidet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jremin_ndet .gt. 0)          &
        used = send_data(cobalt%id_jremin_ndet, cobalt%jremin_ndet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jremin_pdet .gt. 0)          &
        used = send_data(cobalt%id_jremin_pdet, cobalt%jremin_pdet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jremin_fedet .gt. 0)          &
        used = send_data(cobalt%id_jremin_fedet, cobalt%jremin_fedet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jfe_ads .gt. 0)              &
         used = send_data(cobalt%id_jfe_ads,       cobalt%jfe_ads*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jfe_coast .gt. 0)            &
         used = send_data(cobalt%id_jfe_coast, cobalt%jfe_coast*rho_dzt,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_kfe_eq_lig .gt. 0)            &
         used = send_data(cobalt%id_kfe_eq_lig, log10(cobalt%kfe_eq_lig),         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_expkT .gt. 0)              &
         used = send_data(cobalt%id_expkT,       cobalt%expkT,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_hp_temp_lim .gt. 0)            &
         used = send_data(cobalt%id_hp_temp_lim, cobalt%hp_temp_lim,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_irr_inst .gt. 0)           &
         used = send_data(cobalt%id_irr_inst,      cobalt%irr_inst,              &
         model_time, rmask = grid_tmask(:,:,:),&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)             
    if (cobalt%id_irr_mix .gt. 0)           &
         used = send_data(cobalt%id_irr_mix,       cobalt%irr_mix,               &
         model_time, rmask = grid_tmask(:,:,:),&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jno3denit_wc .gt. 0)            &
         used = send_data(cobalt%id_jno3denit_wc,  cobalt%jno3denit_wc*rho_dzt,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jnitrif .gt. 0)              &
         used = send_data(cobalt%id_jnitrif,       cobalt%jnitrif*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_c .gt. 0)  &
         used = send_data(cobalt%id_tot_layer_int_c, cobalt%tot_layer_int_c,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_fe .gt. 0)  &
         used = send_data(cobalt%id_tot_layer_int_fe,cobalt%tot_layer_int_fe,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_n .gt. 0)  &
         used = send_data(cobalt%id_tot_layer_int_n,cobalt%tot_layer_int_n,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_p .gt. 0)  &
         used = send_data(cobalt%id_tot_layer_int_p,cobalt%tot_layer_int_p,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_si .gt. 0)  &
         used = send_data(cobalt%id_tot_layer_int_si,cobalt%tot_layer_int_si,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_total_filter_feeding .gt. 0)  &
         used = send_data(cobalt%id_total_filter_feeding,cobalt%total_filter_feeding,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_net_prim_prod .gt. 0)  &
         used = send_data(cobalt%id_net_prim_prod,cobalt%net_prim_prod*rho_dzt,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_gross_prim_prod.gt. 0)  &
         used = send_data(cobalt%id_gross_prim_prod,cobalt%gross_prim_prod*rho_dzt,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_nlg_diatoms.gt. 0)  &
         used = send_data(cobalt%id_nlg_diatoms,cobalt%nlg_diatoms,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_q_si_2_n_lg_diatoms.gt. 0)  &
         used = send_data(cobalt%id_q_si_2_n_lg_diatoms,cobalt%q_si_2_n_lg_diatoms,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_co2_csurf .gt. 0)             &
         used = send_data(cobalt%id_co2_csurf,      cobalt%co2_csurf,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_pco2_csurf .gt. 0)             &
         used = send_data(cobalt%id_pco2_csurf,      cobalt%pco2_csurf,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_co2_alpha .gt. 0)             &
         used = send_data(cobalt%id_co2_alpha,      cobalt%co2_alpha,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_arag_btm .gt. 0)           &  
         used = send_data(cobalt%id_fcadet_arag_btm,   cobalt%fcadet_arag_btm,      &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_calc_btm .gt. 0)           &
         used = send_data(cobalt%id_fcadet_calc_btm,   cobalt%fcadet_calc_btm,      &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_ffedet_btm .gt. 0)           &
         used = send_data(cobalt%id_ffedet_btm,   cobalt%ffedet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fndet_btm .gt. 0)            &
         used = send_data(cobalt%id_fndet_btm,    cobalt%fndet_btm,              &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fpdet_btm .gt. 0)            &
         used = send_data(cobalt%id_fpdet_btm,    cobalt%fpdet_btm,              &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fsidet_btm .gt. 0)           &
         used = send_data(cobalt%id_fsidet_btm,   cobalt%fsidet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_flithdet_btm .gt. 0)           &
         used = send_data(cobalt%id_flithdet_btm,   cobalt%flithdet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcased_burial .gt. 0)        &
         used = send_data(cobalt%id_fcased_burial, cobalt%fcased_burial,         &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fcased_input .gt. 0)           &
         used = send_data(cobalt%id_fcased_input,  cobalt%fcased_input,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fcased_redis .gt. 0)         &
         used = send_data(cobalt%id_fcased_redis,  cobalt%fcased_redis,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_ffe_sed .gt. 0)              &
         used = send_data(cobalt%id_ffe_sed,       cobalt%ffe_sed,               &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fnfeso4red_sed .gt. 0)           &
         used = send_data(cobalt%id_fnfeso4red_sed,cobalt%fnfeso4red_sed,        &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fno3denit_sed .gt. 0)           &
         used = send_data(cobalt%id_fno3denit_sed, cobalt%fno3denit_sed,         &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fnoxic_sed .gt. 0)           &
         used = send_data(cobalt%id_fnoxic_sed,    cobalt%fnoxic_sed,            &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_frac_burial .gt. 0)           &
         used = send_data(cobalt%id_frac_burial,    cobalt%frac_burial,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fndet_burial .gt. 0)           &
         used = send_data(cobalt%id_fndet_burial,    cobalt%fndet_burial,        &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fpdet_burial .gt. 0)           &
         used = send_data(cobalt%id_fpdet_burial,    cobalt%fpdet_burial,            &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_co3_sol_arag .gt. 0)        &
       used = send_data(cobalt%id_co3_sol_arag,  cobalt%co3_sol_arag,             &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_co3_sol_calc .gt. 0)        &
       used = send_data(cobalt%id_co3_sol_calc,  cobalt%co3_sol_calc,             &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_omega_arag .gt. 0)          &
       used = send_data(cobalt%id_omega_arag,  cobalt%omega_arag,                 &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_omega_calc .gt. 0)          &
       used = send_data(cobalt%id_omega_calc,  cobalt%omega_calc,                 &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fcadet_arag .gt. 0)               &
         used = send_data(cobalt%id_fcadet_arag, cobalt%p_cadet_arag(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:), &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fcadet_calc .gt. 0)               &
         used = send_data(cobalt%id_fcadet_calc, cobalt%p_cadet_calc(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink*grid_tmask(:,:,:), &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_ffedet .gt. 0)                &
         used = send_data(cobalt%id_ffedet,        cobalt%p_fedet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_flithdet .gt. 0)                &
         used = send_data(cobalt%id_flithdet,      cobalt%p_lithdet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fndet .gt. 0)                 &
         used = send_data(cobalt%id_fndet,         cobalt%p_ndet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fpdet .gt. 0)                 &
         used = send_data(cobalt%id_fpdet,         cobalt%p_pdet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fsidet .gt. 0)                &
         used = send_data(cobalt%id_fsidet,        cobalt%p_sidet(:,:,:,tau)  * cobalt%Rho_0 * &
         cobalt%wsink  *grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_nphyto_tot .gt. 0)                &
         used = send_data(cobalt%id_nphyto_tot,   (cobalt%p_ndi(:,:,:,tau) +  &
         cobalt%p_nlg(:,:,:,tau) + cobalt%p_nsm(:,:,:,tau)), &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !
    ! 2D COBALT fields
    !
   if (cobalt%id_pco2surf .gt. 0)              &
         used = send_data(cobalt%id_pco2surf,      cobalt%pco2_csurf,           &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_alk .gt. 0)             &
       used = send_data(cobalt%id_sfc_alk,       cobalt%p_alk(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_cadet_arag .gt. 0)      &
       used = send_data(cobalt%id_sfc_cadet_arag,cobalt%p_cadet_arag(:,:,1,tau),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_cadet_calc .gt. 0)      &
       used = send_data(cobalt%id_sfc_cadet_calc,cobalt%p_cadet_calc(:,:,1,tau),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_dic .gt. 0)             &
       used = send_data(cobalt%id_sfc_dic,       cobalt%p_dic(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_fed .gt. 0)             &
       used = send_data(cobalt%id_sfc_fed,       cobalt%p_fed(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_ldon .gt. 0)             &
       used = send_data(cobalt%id_sfc_ldon,       cobalt%p_ldon(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_sldon .gt. 0)             &
       used = send_data(cobalt%id_sfc_sldon,       cobalt%p_sldon(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_srdon .gt. 0)             &
       used = send_data(cobalt%id_sfc_srdon,       cobalt%p_srdon(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_no3 .gt. 0)             &
       used = send_data(cobalt%id_sfc_no3,       cobalt%p_no3(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_nh4 .gt. 0)             &
       used = send_data(cobalt%id_sfc_nh4,       cobalt%p_nh4(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_po4 .gt. 0)             &
       used = send_data(cobalt%id_sfc_po4,       cobalt%p_po4(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_sio4 .gt. 0)            &
       used = send_data(cobalt%id_sfc_sio4,      cobalt%p_sio4(:,:,1,tau),        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_htotal .gt. 0)          &
       used = send_data(cobalt%id_sfc_htotal,    cobalt%f_htotal(:,:,1),          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_o2 .gt. 0)             &
       used = send_data(cobalt%id_sfc_o2,       cobalt%p_o2(:,:,1,tau),           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_chl .gt. 0)             &
       used = send_data(cobalt%id_sfc_chl,       cobalt%f_chl(:,:,1),             &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_irr .gt. 0)              &
         used = send_data(cobalt%id_sfc_irr,      cobalt%irr_inst(:,:,1),        &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_irr_mem .gt. 0)          &
         used = send_data(cobalt%id_sfc_irr_mem,  cobalt%f_irr_mem(:,:,1),       &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_temp .gt. 0)            &
       used = send_data(cobalt%id_sfc_temp,      Temp(:,:,1),                    &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_temp .gt. 0)            &
       used = send_data(cobalt%id_btm_temp,      cobalt%btm_temp,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_o2 .gt. 0)            &
       used = send_data(cobalt%id_btm_o2,      cobalt%btm_o2,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_co3_ion .gt. 0)            &
       used = send_data(cobalt%id_sfc_co3_ion, cobalt%f_co3_ion(:,:,1),            &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_co3_sol_arag .gt. 0)            &
       used = send_data(cobalt%id_sfc_co3_sol_arag, cobalt%co3_sol_arag(:,:,1),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_co3_sol_calc .gt. 0)            &
       used = send_data(cobalt%id_sfc_co3_sol_calc, cobalt%co3_sol_calc(:,:,1),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    do n= 1, NUM_PHYTO
       if (phyto(n)%id_sfc_f_n .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_f_n, phyto(n)%f_n(:,:,1),            & 
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_sfc_chl .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_chl,    cobalt%c_2_n * 12.0e6 *      &
          phyto(n)%theta(:,:,1) * phyto(n)%f_n(:,:,1),                          &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_sfc_def_fe .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_def_fe, phyto(n)%def_fe(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_sfc_felim .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_felim, phyto(n)%felim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_sfc_irrlim .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_irrlim, phyto(n)%irrlim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (phyto(n)%id_sfc_theta .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_theta, phyto(n)%theta(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (phyto(n)%id_sfc_mu .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_mu, phyto(n)%mu(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(n)%id_sfc_po4lim .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_po4lim, phyto(n)%po4lim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(n)%id_sfc_q_fe_2_n .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_q_fe_2_n, phyto(n)%q_fe_2_n(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    enddo

   do n= 2,3
    if (phyto(n)%id_sfc_nh4lim .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_nh4lim, phyto(n)%nh4lim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (phyto(n)%id_sfc_no3lim .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_no3lim, phyto(n)%no3lim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   enddo


    ! 
    ! Save river, depositon and bulk elemental fluxes
    !
    if (cobalt%id_dep_dry_fed .gt. 0)     &
       used = send_data(cobalt%id_dep_dry_fed, cobalt%dry_fed,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_lith .gt. 0)     &
       used = send_data(cobalt%id_dep_dry_lith, cobalt%dry_lith,                        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_nh4 .gt. 0)     &
       used = send_data(cobalt%id_dep_dry_nh4, cobalt%dry_nh4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_no3 .gt. 0)     &
       used = send_data(cobalt%id_dep_dry_no3, cobalt%dry_no3,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_po4 .gt. 0)     &
       used = send_data(cobalt%id_dep_dry_po4, cobalt%dry_po4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_fed .gt. 0)     &
       used = send_data(cobalt%id_dep_wet_fed, cobalt%wet_fed,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_lith .gt. 0)     &
       used = send_data(cobalt%id_dep_wet_lith, cobalt%wet_lith,                        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_nh4 .gt. 0)     &
       used = send_data(cobalt%id_dep_wet_nh4, cobalt%wet_nh4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_no3 .gt. 0)     &
       used = send_data(cobalt%id_dep_wet_no3, cobalt%wet_no3,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_po4 .gt. 0)     &
       used = send_data(cobalt%id_dep_wet_po4, cobalt%wet_po4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_alk .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_alk, cobalt%runoff_flux_alk,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_dic .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_dic, cobalt%runoff_flux_dic,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_fed .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_fed, cobalt%runoff_flux_fed,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_lith .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_lith, cobalt%runoff_flux_lith,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_no3 .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_no3, cobalt%runoff_flux_no3,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_ldon .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_ldon, cobalt%runoff_flux_ldon,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_sldon .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_sldon, cobalt%runoff_flux_sldon,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_srdon .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_srdon, cobalt%runoff_flux_srdon,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_ndet .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_ndet, cobalt%runoff_flux_ndet,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_po4 .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_po4, cobalt%runoff_flux_po4,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_ldop .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_ldop, cobalt%runoff_flux_ldop,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_sldop .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_sldop, cobalt%runoff_flux_sldop,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_srdop .gt. 0)     &
       used = send_data(cobalt%id_runoff_flux_srdop, cobalt%runoff_flux_srdop,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    !
    ! Save 100m integral fluxes
    !
    if (cobalt%id_jprod_allphytos_100 .gt. 0)     &
       used = send_data(cobalt%id_jprod_allphytos_100, cobalt%jprod_allphytos_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    do n= 1, NUM_PHYTO  !{
       if (phyto(n)%id_jprod_n_100 .gt. 0)     &
          used = send_data(phyto(n)%id_jprod_n_100, phyto(n)%jprod_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_jprod_n_new_100 .gt. 0)     &
          used = send_data(phyto(n)%id_jprod_n_new_100, phyto(n)%jprod_n_new_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_jzloss_n_100 .gt. 0)     &
          used = send_data(phyto(n)%id_jzloss_n_100, phyto(n)%jzloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_jexuloss_n_100 .gt. 0)     &
          used = send_data(phyto(n)%id_jexuloss_n_100, phyto(n)%jexuloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_f_n_100 .gt. 0)     &
          used = send_data(phyto(n)%id_f_n_100, phyto(n)%f_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    enddo !} n
    if (phyto(DIAZO)%id_jprod_n_n2_100 .gt. 0)     &
       used = send_data(phyto(DIAZO)%id_jprod_n_n2_100, phyto(DIAZO)%jprod_n_n2_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(SMALL)%id_jvirloss_n_100 .gt. 0)     &
       used = send_data(phyto(SMALL)%id_jvirloss_n_100, phyto(SMALL)%jvirloss_n_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   if (phyto(SMALL)%id_jaggloss_n_100 .gt. 0)     &
       used = send_data(phyto(SMALL)%id_jaggloss_n_100, phyto(SMALL)%jaggloss_n_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   if (phyto(LARGE)%id_jaggloss_n_100 .gt. 0)     &
       used = send_data(phyto(LARGE)%id_jaggloss_n_100, phyto(LARGE)%jaggloss_n_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     do n= 1, NUM_ZOO  !{
       if (zoo(n)%id_jprod_n_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jprod_n_100, zoo(n)%jprod_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jingest_n_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jingest_n_100, zoo(n)%jingest_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jremin_n_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jremin_n_100, zoo(n)%jremin_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (zoo(n)%id_f_n_100 .gt. 0)     &
          used = send_data(zoo(n)%id_f_n_100, zoo(n)%f_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     enddo !} n

     do n= 1,2  !{
       if (zoo(n)%id_jzloss_n_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jzloss_n_100, zoo(n)%jzloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jprod_don_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jprod_don_100, zoo(n)%jprod_don_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     enddo !} n

     do n= 2,3  !{
       if (zoo(n)%id_jhploss_n_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jhploss_n_100, zoo(n)%jhploss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jprod_ndet_100 .gt. 0)     &
          used = send_data(zoo(n)%id_jprod_ndet_100, zoo(n)%jprod_ndet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     enddo !} n

     if (cobalt%id_hp_jingest_n_100 .gt. 0)     &
        used = send_data(cobalt%id_hp_jingest_n_100, cobalt%hp_jingest_n_100,         &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_hp_jremin_n_100 .gt. 0)     &
        used = send_data(cobalt%id_hp_jremin_n_100, cobalt%hp_jremin_n_100,         &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_hp_jprod_ndet_100 .gt. 0)     &
        used = send_data(cobalt%id_hp_jprod_ndet_100, cobalt%hp_jprod_ndet_100,         &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (bact(1)%id_jprod_n_100 .gt. 0)     &
          used = send_data(bact(1)%id_jprod_n_100, bact(1)%jprod_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_jzloss_n_100 .gt. 0)     &
          used = send_data(bact(1)%id_jzloss_n_100, bact(1)%jzloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_jvirloss_n_100 .gt. 0)     &
          used = send_data(bact(1)%id_jvirloss_n_100, bact(1)%jvirloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_jremin_n_100 .gt. 0)     &
          used = send_data(bact(1)%id_jremin_n_100, bact(1)%jremin_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_juptake_ldon_100 .gt. 0)     &
          used = send_data(bact(1)%id_juptake_ldon_100, bact(1)%juptake_ldon_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_f_n_100 .gt. 0)     &
          used = send_data(bact(1)%id_f_n_100, bact(1)%f_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (cobalt%id_jprod_lithdet_100 .gt. 0)     &
          used = send_data(cobalt%id_jprod_lithdet_100, cobalt%jprod_lithdet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_sidet_100 .gt. 0)     &
          used = send_data(cobalt%id_jprod_sidet_100, cobalt%jprod_sidet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_cadet_calc_100 .gt. 0)     &
          used = send_data(cobalt%id_jprod_cadet_calc_100, cobalt%jprod_cadet_calc_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_cadet_arag_100 .gt. 0)     &
          used = send_data(cobalt%id_jprod_cadet_arag_100, cobalt%jprod_cadet_arag_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jremin_ndet_100 .gt. 0)     &
          used = send_data(cobalt%id_jremin_ndet_100, cobalt%jremin_ndet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_mesozoo_200 .gt. 0)     &
          used = send_data(cobalt%id_jprod_mesozoo_200, cobalt%jprod_mesozoo_200,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (cobalt%id_f_ndet_100 .gt. 0)     &
          used = send_data(cobalt%id_f_ndet_100, cobalt%f_ndet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_f_don_100 .gt. 0)     &
          used = send_data(cobalt%id_f_don_100, cobalt%f_don_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_f_silg_100 .gt. 0)     &
          used = send_data(cobalt%id_f_silg_100, cobalt%f_silg_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_f_mesozoo_200 .gt. 0)     &
          used = send_data(cobalt%id_f_mesozoo_200, cobalt%f_mesozoo_200,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_fndet_100 .gt. 0)           &
       used = send_data(cobalt%id_fndet_100,     cobalt%fndet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fpdet_100 .gt. 0)           &
       used = send_data(cobalt%id_fpdet_100,     cobalt%fpdet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fsidet_100 .gt. 0)           &
       used = send_data(cobalt%id_fsidet_100,     cobalt%fsidet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_flithdet_100 .gt. 0)           &
       used = send_data(cobalt%id_flithdet_100,     cobalt%flithdet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_calc_100 .gt. 0)           &
       used = send_data(cobalt%id_fcadet_calc_100,     cobalt%fcadet_calc_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_arag_100 .gt. 0)           &
       used = send_data(cobalt%id_fcadet_arag_100,     cobalt%fcadet_arag_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_ffedet_100 .gt. 0)           &
       used = send_data(cobalt%id_ffedet_100,     cobalt%ffedet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    !
    !---------------------------------------------------------------------
    ! Save CaCO3 saturation and O2 minimum depths
    !---------------------------------------------------------------------
    !
    if (cobalt%id_o2min .gt. 0)               &
       used = send_data(cobalt%id_o2min,         cobalt%o2min,                    &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_z_o2min .gt. 0)             &
       used = send_data(cobalt%id_z_o2min,    cobalt%z_o2min,                     &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_z_sat_arag .gt. 0)          &
       used = send_data(cobalt%id_z_sat_arag,    cobalt%z_sat_arag,               &
       model_time, mask = cobalt%mask_z_sat_arag, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_z_sat_calc .gt. 0)          &
       used = send_data(cobalt%id_z_sat_calc,    cobalt%z_sat_calc,               &
       model_time, mask = cobalt%mask_z_sat_calc, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    call mpp_clock_end(id_clock_cobalt_send_diagnostics)

  end subroutine generic_COBALT_update_from_source

  ! <SUBROUTINE NAME="generic_COBALT_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_COBALT_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),   intent(in)   :: SST, SSS
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,tau

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: sal,ST,o2_saturation
    real    :: tt,tk,ts,ts2,ts3,ts4,ts5
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: o2_field,dic_field,po4_field,sio4_field,alk_field
    real, dimension(:,:,:), ALLOCATABLE :: htotal_field,co3_ion_field
    real, dimension(:,:), ALLOCATABLE :: co2_alpha,co2_csurf,co2_sc_no,o2_alpha,o2_csurf,o2_sc_no
    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_set_boundary_values'
    !
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    call g_tracer_get_pointer(tracer_list,'o2' ,'field',  o2_field)

    allocate(co2_alpha(isd:ied, jsd:jed)); co2_alpha=0.0
    allocate(co2_csurf(isd:ied, jsd:jed)); co2_csurf=0.0
    allocate(co2_sc_no(isd:ied, jsd:jed)); co2_sc_no=0.0
    allocate( o2_alpha(isd:ied, jsd:jed)); o2_alpha=0.0
    allocate( o2_csurf(isd:ied, jsd:jed)); o2_csurf=0.0
    allocate( o2_sc_no(isd:ied, jsd:jed)); o2_sc_no=0.0
    allocate(htotal_field(isd:ied,jsd:jed,nk),co3_ion_field(isd:ied,jsd:jed,nk))
    htotal_field=0.0 ; co3_ion_field=0.0

    !nnz: Since the generic_COBALT_update_from_source() subroutine is called by this time
    !     the following if block is not really necessary (since this calculation is already done in source).
    !     It is only neccessary if source routine is commented out for debugging.
    !Note: In order for this to work we should NOT zero out the coupler values for generic tracers
    !      This zeroing is done for non-generic TOPAZ by calling zero_ocean_sfc.
    !      Since the coupler values here are non-cumulative there is no need to zero them out anyway.

    if (cobalt%init ) then
       !Get necessary fields
       call g_tracer_get_pointer(tracer_list,'dic'   ,'field', dic_field)
       call g_tracer_get_pointer(tracer_list,'po4'   ,'field', po4_field)
       call g_tracer_get_pointer(tracer_list,'sio4'  ,'field', sio4_field)
       call g_tracer_get_pointer(tracer_list,'alk'   ,'field', alk_field)

       call g_tracer_get_values(tracer_list,'htotal' ,'field', htotal_field,isd,jsd,ntau=1)
       call g_tracer_get_values(tracer_list,'co3_ion','field',co3_ion_field,isd,jsd,ntau=1)

       do j = jsc, jec ; do i = isc, iec  !{
          cobalt%htotallo(i,j) = cobalt%htotal_scale_lo * htotal_field(i,j,1)
          cobalt%htotalhi(i,j) = cobalt%htotal_scale_hi * htotal_field(i,j,1)
       enddo; enddo ; !} i, j

       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                            &
            dic_field(:,:,1,tau),                          &
            po4_field(:,:,1,tau),                          &
            sio4_field(:,:,1,tau),                         &
            alk_field(:,:,1,tau),                          &
            cobalt%htotallo, cobalt%htotalhi,                &
                                !InOut
            htotal_field(:,:,1),                           &
                                !OUT
            co2star=co2_csurf(:,:), alpha=co2_alpha(:,:),  &
            pCO2surf=cobalt%pco2_csurf(:,:), &
            co3_ion=co3_ion_field(:,:,1))

       !Set fields !nnz: if These are pointers do I need to do this?
       call g_tracer_set_values(tracer_list,'htotal' ,'field',htotal_field ,isd,jsd,ntau=1)
       call g_tracer_set_values(tracer_list,'co3_ion','field',co3_ion_field,isd,jsd,ntau=1)

       call g_tracer_set_values(tracer_list,'dic','alpha',co2_alpha    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'dic','csurf',co2_csurf    ,isd,jsd)

       !!nnz: If source is called uncomment the following
       cobalt%init = .false. !nnz: This is necessary since the above two calls appear in source subroutine too.
    endif

    call g_tracer_get_values(tracer_list,'dic','alpha', co2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'dic','csurf', co2_csurf ,isd,jsd)

    call g_tracer_get_values(tracer_list,'o2','alpha', o2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'o2','csurf', o2_csurf ,isd,jsd)
    do j=jsc,jec ; do i=isc,iec
       !This calculation needs an input of SST and SSS
       sal = SSS(i,j) ; ST = SST(i,j)

       !nnz:
       !Note: In the following calculations in order to get results for co2 and o2
       !      identical with cobalt code in MOM cobalt%Rho_0 must be replaced with rho(i,j,1,tau)
       !      This is achieved by uncommenting the following if desired.
       !! cobalt%Rho_0 = rho(i,j,1,tau)
       !      But since %Rho_0 plays the role of a unit conversion factor in this module
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
       co2_sc_no(i,j) = cobalt%a1_co2 + ST * (cobalt%a2_co2 + ST * (cobalt%a3_co2 + ST * cobalt%a4_co2)) * &
            grid_tmask(i,j,1)
!       sc_no_term = sqrt(660.0 / (sc_co2 + epsln))
!
!       co2_alpha(i,j) = co2_alpha(i,j)* sc_no_term * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)
!       co2_csurf(i,j) = co2_csurf(i,j)* sc_no_term * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)
!
! in 'ocmip2_new' atmos_ocean_fluxes.F90 coupler formulation, the schmidt number is carried in explicitly
!
       co2_alpha(i,j) = co2_alpha(i,j) * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)
       co2_csurf(i,j) = co2_csurf(i,j) * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)

       !---------------------------------------------------------------------
       !     O2
       !---------------------------------------------------------------------
       !  Compute the oxygen saturation concentration at 1 atm total
       !  pressure in mol/kg given the temperature (t, in deg C) and
       !  the salinity (s, in permil)
       !
       !  From Garcia and Gosrdon (1992), Limnology and Oceonography.
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
            exp( cobalt%a_0 + cobalt%a_1*ts + cobalt%a_2*ts2 + cobalt%a_3*ts3 + cobalt%a_4*ts4 + cobalt%a_5*ts5 + &
            (cobalt%b_0 + cobalt%b_1*ts + cobalt%b_2*ts2 + cobalt%b_3*ts3 + cobalt%c_0*sal)*sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the
       !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
       !  Cycles, 12, 141-163).
       !---------------------------------------------------------------------
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       o2_sc_no(i,j)  = cobalt%a1_o2  + ST * (cobalt%a2_o2  + ST * (cobalt%a3_o2  + ST * cobalt%a4_o2 )) * &
            grid_tmask(i,j,1)
       !
       !      renormalize the alpha value for atm o2
       !      data table override for o2_flux_pcair_atm is now set to 0.21
       !
       o2_alpha(i,j) = (o2_saturation / 0.21)
       o2_csurf(i,j) = o2_field(i,j,1,tau) * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)


    enddo; enddo

    !
    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'dic','alpha',co2_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','csurf',co2_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','sc_no',co2_sc_no,isd,jsd)

    call g_tracer_set_values(tracer_list,'o2', 'alpha',o2_alpha, isd,jsd)
    call g_tracer_set_values(tracer_list,'o2', 'csurf',o2_csurf, isd,jsd)
    call g_tracer_set_values(tracer_list,'o2', 'sc_no',o2_sc_no, isd,jsd)

    deallocate(co2_alpha,co2_csurf,co2_sc_no,o2_alpha,o2_csurf,o2_sc_no)

  end subroutine generic_COBALT_set_boundary_values


  ! <SUBROUTINE NAME="generic_COBALT_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_end
  !  </TEMPLATE>
  ! </SUBROUTINE>


  subroutine generic_COBALT_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_end'
    call user_deallocate_arrays
  end subroutine generic_COBALT_end

  !
  !   This is an internal sub, not a public interface.
  !   Allocate all the work arrays to be used in this module.
  !
  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    !Allocate all the private arrays.

    !Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc ; CO2_dope_vec%iec = iec 
    CO2_dope_vec%jsc = jsc ; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd ; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd ; CO2_dope_vec%jed = jed    

    allocate(cobalt%htotallo(isd:ied,jsd:jed))
    allocate(cobalt%htotalhi(isd:ied,jsd:jed))

    !
    ! allocate and initialize array elements of all phytoplankton groups
    ! CAS: add fluxes for additional explicit phytoplankton loss terms

    do n = 1, NUM_PHYTO
       allocate(phyto(n)%def_fe(isd:ied,jsd:jed,nk))       ; phyto(n)%def_fe         = 0.0
       allocate(phyto(n)%def_p(isd:ied,jsd:jed,nk))        ; phyto(n)%def_p          = 0.0
       allocate(phyto(n)%f_fe(isd:ied,jsd:jed,nk))         ; phyto(n)%f_fe           = 0.0
       allocate(phyto(n)%f_n(isd:ied,jsd:jed,nk))          ; phyto(n)%f_n            = 0.0
       allocate(phyto(n)%felim(isd:ied,jsd:jed,nk))        ; phyto(n)%felim          = 0.0
       allocate(phyto(n)%irrlim(isd:ied,jsd:jed,nk))       ; phyto(n)%irrlim         = 0.0
       allocate(phyto(n)%jzloss_fe(isd:ied,jsd:jed,nk))    ; phyto(n)%jzloss_fe      = 0.0
       allocate(phyto(n)%jzloss_n(isd:ied,jsd:jed,nk))     ; phyto(n)%jzloss_n       = 0.0
       allocate(phyto(n)%jzloss_p(isd:ied,jsd:jed,nk))     ; phyto(n)%jzloss_p       = 0.0
       allocate(phyto(n)%jzloss_sio2(isd:ied,jsd:jed,nk))  ; phyto(n)%jzloss_sio2    = 0.0
       allocate(phyto(n)%jaggloss_fe(isd:ied,jsd:jed,nk))  ; phyto(n)%jaggloss_fe    = 0.0
       allocate(phyto(n)%jaggloss_n(isd:ied,jsd:jed,nk))   ; phyto(n)%jaggloss_n     = 0.0
       allocate(phyto(n)%jaggloss_p(isd:ied,jsd:jed,nk))   ; phyto(n)%jaggloss_p     = 0.0
       allocate(phyto(n)%jaggloss_sio2(isd:ied,jsd:jed,nk)); phyto(n)%jaggloss_sio2  = 0.0
       allocate(phyto(n)%jvirloss_fe(isd:ied,jsd:jed,nk))  ; phyto(n)%jvirloss_fe    = 0.0
       allocate(phyto(n)%jvirloss_n(isd:ied,jsd:jed,nk))   ; phyto(n)%jvirloss_n     = 0.0
       allocate(phyto(n)%jvirloss_p(isd:ied,jsd:jed,nk))   ; phyto(n)%jvirloss_p     = 0.0
       allocate(phyto(n)%jvirloss_sio2(isd:ied,jsd:jed,nk)); phyto(n)%jvirloss_sio2  = 0.0
       allocate(phyto(n)%jexuloss_fe(isd:ied,jsd:jed,nk))  ; phyto(n)%jexuloss_fe    = 0.0
       allocate(phyto(n)%jexuloss_n(isd:ied,jsd:jed,nk))   ; phyto(n)%jexuloss_n     = 0.0
       allocate(phyto(n)%jexuloss_p(isd:ied,jsd:jed,nk))   ; phyto(n)%jexuloss_p     = 0.0
       allocate(phyto(n)%jhploss_fe(isd:ied,jsd:jed,nk))   ; phyto(n)%jhploss_fe     = 0.0
       allocate(phyto(n)%jhploss_n(isd:ied,jsd:jed,nk))    ; phyto(n)%jhploss_n      = 0.0
       allocate(phyto(n)%jhploss_p(isd:ied,jsd:jed,nk))    ; phyto(n)%jhploss_p      = 0.0
       allocate(phyto(n)%jhploss_sio2(isd:ied,jsd:jed,nk)) ; phyto(n)%jhploss_sio2   = 0.0
       allocate(phyto(n)%juptake_fe(isd:ied,jsd:jed,nk))   ; phyto(n)%juptake_fe     = 0.0
       allocate(phyto(n)%juptake_nh4(isd:ied,jsd:jed,nk))  ; phyto(n)%juptake_nh4    = 0.0
       allocate(phyto(n)%juptake_no3(isd:ied,jsd:jed,nk))  ; phyto(n)%juptake_no3    = 0.0
       allocate(phyto(n)%juptake_po4(isd:ied,jsd:jed,nk))  ; phyto(n)%juptake_po4    = 0.0
       allocate(phyto(n)%jprod_n(isd:ied,jsd:jed,nk))      ; phyto(n)%jprod_n        = 0.0
       allocate(phyto(n)%liebig_lim(isd:ied,jsd:jed,nk))   ; phyto(n)%liebig_lim     = 0.0
       allocate(phyto(n)%mu(isd:ied,jsd:jed,nk))           ; phyto(n)%mu             = 0.0
       allocate(phyto(n)%po4lim(isd:ied,jsd:jed,nk))       ; phyto(n)%po4lim         = 0.0
       allocate(phyto(n)%q_fe_2_n(isd:ied,jsd:jed,nk))     ; phyto(n)%q_fe_2_n       = 0.0
       allocate(phyto(n)%q_p_2_n(isd:ied,jsd:jed,nk))      ; phyto(n)%q_p_2_n        = 0.0
       allocate(phyto(n)%q_si_2_n(isd:ied,jsd:jed,nk))     ; phyto(n)%q_si_2_n       = 0.0
       allocate(phyto(n)%theta(isd:ied,jsd:jed,nk))        ; phyto(n)%theta          = 0.0
    enddo
    !
    ! allocate and initialize array elements of only some phytoplankton groups
    !
    do n = 2, NUM_PHYTO
       allocate(phyto(n)%nh4lim(isd:ied,jsd:jed,nk))      ; phyto(n)%nh4lim          = 0.0
       allocate(phyto(n)%no3lim(isd:ied,jsd:jed,nk))      ; phyto(n)%no3lim          = 0.0
    enddo
    !
    ! allocate and initialize array elements of only one phytoplankton group
    !
    allocate(phyto(DIAZO)%juptake_n2(isd:ied,jsd:jed,nk))   ; phyto(DIAZO)%juptake_n2   = 0.0
    allocate(phyto(DIAZO)%o2lim(isd:ied,jsd:jed,nk))        ; phyto(DIAZO)%o2lim        = 0.0
    allocate(phyto(LARGE)%juptake_sio4(isd:ied,jsd:jed,nk)) ; phyto(LARGE)%juptake_sio4 = 0.0
    allocate(phyto(LARGE)%silim(isd:ied,jsd:jed,nk))        ; phyto(LARGE)%silim      = 0.0
    !
    ! allocate and initialize arrays for bacteria
    !
    allocate(bact(1)%f_n(isd:ied,jsd:jed,nk))              ; bact(1)%f_n             = 0.0
    allocate(bact(1)%jzloss_n(isd:ied,jsd:jed,nk))         ; bact(1)%jzloss_n        = 0.0
    allocate(bact(1)%jzloss_p(isd:ied,jsd:jed,nk))         ; bact(1)%jzloss_p        = 0.0
    allocate(bact(1)%jvirloss_n(isd:ied,jsd:jed,nk))       ; bact(1)%jvirloss_n      = 0.0
    allocate(bact(1)%jvirloss_p(isd:ied,jsd:jed,nk))       ; bact(1)%jvirloss_p      = 0.0
    allocate(bact(1)%jhploss_n(isd:ied,jsd:jed,nk))        ; bact(1)%jhploss_n       = 0.0
    allocate(bact(1)%jhploss_p(isd:ied,jsd:jed,nk))        ; bact(1)%jhploss_p       = 0.0
    allocate(bact(1)%juptake_ldon(isd:ied,jsd:jed,nk))     ; bact(1)%juptake_ldon    = 0.0
    allocate(bact(1)%juptake_ldop(isd:ied,jsd:jed,nk))     ; bact(1)%juptake_ldop    = 0.0
    allocate(bact(1)%jprod_nh4(isd:ied,jsd:jed,nk))        ; bact(1)%jprod_nh4       = 0.0
    allocate(bact(1)%jprod_po4(isd:ied,jsd:jed,nk))        ; bact(1)%jprod_po4       = 0.0
    allocate(bact(1)%jprod_n(isd:ied,jsd:jed,nk))          ; bact(1)%jprod_n      = 0.0
    allocate(bact(1)%temp_lim(isd:ied,jsd:jed,nk))         ; bact(1)%temp_lim        = 0.0
    !
    ! CAS: allocate and initialize array elements for all zooplankton groups
    !
    do n = 1, NUM_ZOO
       allocate(zoo(n)%f_n(isd:ied,jsd:jed,nk))           ; zoo(n)%f_n            = 0.0
       allocate(zoo(n)%jzloss_n(isd:ied,jsd:jed,nk))      ; zoo(n)%jzloss_n       = 0.0
       allocate(zoo(n)%jzloss_p(isd:ied,jsd:jed,nk))      ; zoo(n)%jzloss_p       = 0.0
       allocate(zoo(n)%jhploss_n(isd:ied,jsd:jed,nk))     ; zoo(n)%jhploss_n      = 0.0
       allocate(zoo(n)%jhploss_p(isd:ied,jsd:jed,nk))     ; zoo(n)%jhploss_p      = 0.0
       allocate(zoo(n)%jingest_n(isd:ied,jsd:jed,nk))     ; zoo(n)%jingest_n      = 0.0
       allocate(zoo(n)%jingest_p(isd:ied,jsd:jed,nk))     ; zoo(n)%jingest_p      = 0.0
       allocate(zoo(n)%jingest_sio2(isd:ied,jsd:jed,nk))  ; zoo(n)%jingest_sio2   = 0.0
       allocate(zoo(n)%jingest_fe(isd:ied,jsd:jed,nk))    ; zoo(n)%jingest_fe     = 0.0
       allocate(zoo(n)%jprod_ndet(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_ndet     = 0.0
       allocate(zoo(n)%jprod_pdet(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_pdet     = 0.0
       allocate(zoo(n)%jprod_ldon(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_ldon     = 0.0
       allocate(zoo(n)%jprod_ldop(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_ldop     = 0.0
       allocate(zoo(n)%jprod_srdon(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_srdon     = 0.0
       allocate(zoo(n)%jprod_srdop(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_srdop     = 0.0
       allocate(zoo(n)%jprod_sldon(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_sldon     = 0.0
       allocate(zoo(n)%jprod_sldop(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_sldop     = 0.0
       allocate(zoo(n)%jprod_fedet(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_fedet    = 0.0
       allocate(zoo(n)%jprod_fed(isd:ied,jsd:jed,nk))     ; zoo(n)%jprod_fed      = 0.0
       allocate(zoo(n)%jprod_sidet(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_sidet    = 0.0
       allocate(zoo(n)%jprod_sio4(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_sio4      = 0.0
       allocate(zoo(n)%jprod_po4(isd:ied,jsd:jed,nk))     ; zoo(n)%jprod_po4      = 0.0
       allocate(zoo(n)%jprod_nh4(isd:ied,jsd:jed,nk))     ; zoo(n)%jprod_nh4      = 0.0
       allocate(zoo(n)%jprod_n(isd:ied,jsd:jed,nk))      ; zoo(n)%jprod_n       = 0.0
       allocate(zoo(n)%temp_lim(isd:ied,jsd:jed,nk))      ; zoo(n)%temp_lim       = 0.0
    enddo

    ! higher predator ingestion
    allocate(cobalt%hp_jingest_n(isd:ied,jsd:jed,nk))     ; cobalt%hp_jingest_n      = 0.0
    allocate(cobalt%hp_jingest_p(isd:ied,jsd:jed,nk))     ; cobalt%hp_jingest_p      = 0.0
    allocate(cobalt%hp_jingest_sio2(isd:ied,jsd:jed,nk))  ; cobalt%hp_jingest_sio2   = 0.0
    allocate(cobalt%hp_jingest_fe(isd:ied,jsd:jed,nk))    ; cobalt%hp_jingest_fe     = 0.0 

    allocate(cobalt%f_alk(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_alk=0.0
    allocate(cobalt%f_cadet_arag(isd:ied, jsd:jed, 1:nk)) ; cobalt%f_cadet_arag=0.0
    allocate(cobalt%f_cadet_calc(isd:ied, jsd:jed, 1:nk)) ; cobalt%f_cadet_calc=0.0
    allocate(cobalt%f_dic(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_dic=0.0
    allocate(cobalt%f_fed(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_fed=0.0
    allocate(cobalt%f_fedet(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_fedet=0.0
    allocate(cobalt%f_ldon(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_ldon=0.0
    allocate(cobalt%f_ldop(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_ldop=0.0
    allocate(cobalt%f_lith(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_lith=0.0
    allocate(cobalt%f_lithdet(isd:ied, jsd:jed, 1:nk))    ; cobalt%f_lithdet=0.0
    allocate(cobalt%f_ndet(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_ndet=0.0
    allocate(cobalt%f_nh4(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_nh4=0.0
    allocate(cobalt%f_no3(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_no3=0.0
    allocate(cobalt%f_o2(isd:ied, jsd:jed, 1:nk))         ; cobalt%f_o2=0.0
    allocate(cobalt%f_pdet(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_pdet=0.0
    allocate(cobalt%f_po4(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_po4=0.0
    allocate(cobalt%f_srdon(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_srdon=0.0
    allocate(cobalt%f_srdop(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_srdop=0.0
    allocate(cobalt%f_sldon(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_sldon=0.0
    allocate(cobalt%f_sldop(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_sldop=0.0
    allocate(cobalt%f_sidet(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_sidet=0.0
    allocate(cobalt%f_silg(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_silg=0.0
    allocate(cobalt%f_sio4(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_sio4=0.0
    allocate(cobalt%co3_sol_arag(isd:ied, jsd:jed, 1:nk)) ; cobalt%co3_sol_arag=0.0
    allocate(cobalt%co3_sol_calc(isd:ied, jsd:jed, 1:nk)) ; cobalt%co3_sol_calc=0.0
    allocate(cobalt%omega_arag(isd:ied, jsd:jed, 1:nk))   ; cobalt%omega_arag=0.0
    allocate(cobalt%omega_calc(isd:ied, jsd:jed, 1:nk))   ; cobalt%omega_calc=0.0
    allocate(cobalt%f_chl(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_chl=0.0
    allocate(cobalt%f_co3_ion(isd:ied, jsd:jed, 1:nk))    ; cobalt%f_co3_ion=0.0
    allocate(cobalt%f_htotal(isd:ied, jsd:jed, 1:nk))     ; cobalt%f_htotal=0.0
    allocate(cobalt%f_irr_mem(isd:ied, jsd:jed, 1:nk))    ; cobalt%f_irr_mem=0.0
    allocate(cobalt%f_cased(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_cased=0.0
    allocate(cobalt%f_cadet_arag_btf(isd:ied, jsd:jed, 1:nk)); cobalt%f_cadet_arag_btf=0.0
    allocate(cobalt%f_cadet_calc_btf(isd:ied, jsd:jed, 1:nk)); cobalt%f_cadet_calc_btf=0.0
    allocate(cobalt%f_fedet_btf(isd:ied, jsd:jed, 1:nk))  ; cobalt%f_fedet_btf=0.0
    allocate(cobalt%f_lithdet_btf(isd:ied, jsd:jed, 1:nk)); cobalt%f_lithdet_btf=0.0
    allocate(cobalt%f_ndet_btf(isd:ied, jsd:jed, 1:nk))   ; cobalt%f_ndet_btf=0.0
    allocate(cobalt%f_pdet_btf(isd:ied, jsd:jed, 1:nk))   ; cobalt%f_pdet_btf=0.0
    allocate(cobalt%f_sidet_btf(isd:ied, jsd:jed, 1:nk))  ; cobalt%f_sidet_btf=0.0

    allocate(cobalt%jndi(isd:ied, jsd:jed, 1:nk))         ; cobalt%jndi=0.0
    allocate(cobalt%jnsm(isd:ied, jsd:jed, 1:nk))         ; cobalt%jnsm=0.0
    allocate(cobalt%jnlg(isd:ied, jsd:jed, 1:nk))         ; cobalt%jnlg=0.0
    allocate(cobalt%jnbact(isd:ied, jsd:jed, 1:nk))       ; cobalt%jnbact=0.0
    allocate(cobalt%jnsmz(isd:ied, jsd:jed, 1:nk))        ; cobalt%jnsmz=0.0
    allocate(cobalt%jnmdz(isd:ied, jsd:jed, 1:nk))        ; cobalt%jnmdz=0.0
    allocate(cobalt%jnlgz(isd:ied, jsd:jed, 1:nk))        ; cobalt%jnlgz=0.0
    allocate(cobalt%jalk(isd:ied, jsd:jed, 1:nk))         ; cobalt%jalk=0.0
    allocate(cobalt%jcadet_arag(isd:ied, jsd:jed, 1:nk))  ; cobalt%jcadet_arag=0.0
    allocate(cobalt%jcadet_calc(isd:ied, jsd:jed, 1:nk))  ; cobalt%jcadet_calc=0.0
    allocate(cobalt%jdic(isd:ied, jsd:jed, 1:nk))         ; cobalt%jdic=0.0
    allocate(cobalt%jfed(isd:ied, jsd:jed, 1:nk))         ; cobalt%jfed=0.0
    allocate(cobalt%jfedi(isd:ied, jsd:jed, 1:nk))        ; cobalt%jfedi=0.0
    allocate(cobalt%jfelg(isd:ied, jsd:jed, 1:nk))        ; cobalt%jfelg=0.0
    allocate(cobalt%jfesm(isd:ied, jsd:jed, 1:nk))     ; cobalt%jfesm=0.0
    allocate(cobalt%jfedet(isd:ied, jsd:jed, 1:nk))       ; cobalt%jfedet=0.0
    allocate(cobalt%jldon(isd:ied, jsd:jed, 1:nk))        ; cobalt%jldon=0.0
    allocate(cobalt%jldop(isd:ied, jsd:jed, 1:nk))        ; cobalt%jldop=0.0
    allocate(cobalt%jlith(isd:ied, jsd:jed, 1:nk))        ; cobalt%jlith=0.0
    allocate(cobalt%jlithdet(isd:ied, jsd:jed, 1:nk))     ; cobalt%jlithdet=0.0
    allocate(cobalt%jndet(isd:ied, jsd:jed, 1:nk))        ; cobalt%jndet=0.0
    allocate(cobalt%jnh4(isd:ied, jsd:jed, 1:nk))         ; cobalt%jnh4=0.0
    allocate(cobalt%jno3(isd:ied, jsd:jed, 1:nk))         ; cobalt%jno3=0.0
    allocate(cobalt%jo2(isd:ied, jsd:jed, 1:nk))          ; cobalt%jo2=0.0
    allocate(cobalt%jpdet(isd:ied, jsd:jed, 1:nk))        ; cobalt%jpdet=0.0
    allocate(cobalt%jpo4(isd:ied, jsd:jed, 1:nk))         ; cobalt%jpo4=0.0
    allocate(cobalt%jsrdon(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsrdon=0.0
    allocate(cobalt%jsrdop(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsrdop=0.0
    allocate(cobalt%jsldon(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsldon=0.0
    allocate(cobalt%jsldop(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsldop=0.0
    allocate(cobalt%jsidet(isd:ied, jsd:jed, 1:nk))       ; cobalt%jsidet=0.0
    allocate(cobalt%jsilg(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsilg=0.0
    allocate(cobalt%jsio4(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsio4=0.0
    allocate(cobalt%jprod_ndet(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_ndet=0.0
    allocate(cobalt%jprod_pdet(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_pdet=0.0
    allocate(cobalt%jprod_srdon(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_srdon=0.0
    allocate(cobalt%jprod_sldon(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_sldon=0.0
    allocate(cobalt%jprod_ldon(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_ldon=0.0
    allocate(cobalt%jprod_srdop(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_srdop=0.0
    allocate(cobalt%jprod_sldop(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_sldop=0.0
    allocate(cobalt%jprod_ldop(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_ldop=0.0
    allocate(cobalt%jprod_fedet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_fedet=0.0
    allocate(cobalt%jprod_fed(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_fed=0.0
    allocate(cobalt%jprod_sidet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_sidet=0.0
    allocate(cobalt%jprod_sio4(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_sio4=0.0
    allocate(cobalt%jprod_lithdet(isd:ied, jsd:jed, 1:nk)); cobalt%jprod_lithdet=0.0
    allocate(cobalt%jprod_cadet_arag(isd:ied, jsd:jed, 1:nk)); cobalt%jprod_cadet_arag=0.0
    allocate(cobalt%jprod_cadet_calc(isd:ied, jsd:jed, 1:nk)); cobalt%jprod_cadet_calc=0.0
    allocate(cobalt%jprod_nh4(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_nh4=0.0
    allocate(cobalt%jprod_po4(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_po4=0.0
    allocate(cobalt%det_jzloss_n(isd:ied, jsd:jed, 1:nk)) ; cobalt%det_jzloss_n=0.0
    allocate(cobalt%det_jzloss_p(isd:ied, jsd:jed, 1:nk)) ; cobalt%det_jzloss_p=0.0
    allocate(cobalt%det_jzloss_fe(isd:ied, jsd:jed, 1:nk)); cobalt%det_jzloss_fe=0.0
    allocate(cobalt%det_jzloss_si(isd:ied, jsd:jed, 1:nk)); cobalt%det_jzloss_si=0.0
    allocate(cobalt%det_jhploss_n(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_n=0.0
    allocate(cobalt%det_jhploss_p(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_p=0.0
    allocate(cobalt%det_jhploss_fe(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_fe=0.0
    allocate(cobalt%det_jhploss_si(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_si=0.0
    allocate(cobalt%jdiss_cadet_arag(isd:ied, jsd:jed, 1:nk)); cobalt%jdiss_cadet_arag=0.0
    allocate(cobalt%jdiss_cadet_calc(isd:ied, jsd:jed, 1:nk)); cobalt%jdiss_cadet_calc=0.0
    allocate(cobalt%jdiss_sidet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jdiss_sidet=0.0
    allocate(cobalt%jremin_ndet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jremin_ndet=0.0
    allocate(cobalt%jremin_pdet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jremin_pdet=0.0
    allocate(cobalt%jremin_fedet(isd:ied, jsd:jed, 1:nk)) ; cobalt%jremin_fedet=0.0
    allocate(cobalt%jfe_ads(isd:ied, jsd:jed, 1:nk))      ; cobalt%jfe_ads=0.0
    allocate(cobalt%jfe_coast(isd:ied, jsd:jed, 1:nk))    ; cobalt%jfe_coast=0.0
    allocate(cobalt%kfe_eq_lig(isd:ied, jsd:jed, 1:nk))   ; cobalt%kfe_eq_lig=0.0
    allocate(cobalt%expkT(isd:ied, jsd:jed, 1:nk))        ; cobalt%expkT=0.0
    allocate(cobalt%hp_temp_lim(isd:ied, jsd:jed, 1:nk))  ; cobalt%hp_temp_lim=0.0
    allocate(cobalt%irr_inst(isd:ied, jsd:jed, 1:nk))     ; cobalt%irr_inst=0.0
    allocate(cobalt%irr_mix(isd:ied, jsd:jed, 1:nk))      ; cobalt%irr_mix=0.0
    allocate(cobalt%jno3denit_wc(isd:ied, jsd:jed, 1:nk)) ; cobalt%jno3denit_wc=0.0
    allocate(cobalt%jnitrif(isd:ied, jsd:jed, 1:nk))      ; cobalt%jnitrif=0.0
    allocate(cobalt%total_filter_feeding(isd:ied,jsd:jed,1:nk)); cobalt%total_filter_feeding=0.0
    allocate(cobalt%zt(isd:ied, jsd:jed, 1:nk))           ; cobalt%zt=0.0
    allocate(cobalt%zm(isd:ied, jsd:jed, 1:nk))           ; cobalt%zm=0.0
    allocate(cobalt%tot_layer_int_c(isd:ied, jsd:jed,1:nk))  ; cobalt%tot_layer_int_c=0.0
    allocate(cobalt%tot_layer_int_fe(isd:ied, jsd:jed,1:nk)) ; cobalt%tot_layer_int_fe=0.0
    allocate(cobalt%tot_layer_int_n(isd:ied, jsd:jed, 1:nk)) ; cobalt%tot_layer_int_n=0.0
    allocate(cobalt%tot_layer_int_p(isd:ied, jsd:jed, 1:nk)) ; cobalt%tot_layer_int_p=0.0
    allocate(cobalt%tot_layer_int_si(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_si=0.0
    allocate(cobalt%net_prim_prod(isd:ied, jsd:jed, 1:nk)); cobalt%net_prim_prod=0.0
    allocate(cobalt%gross_prim_prod(isd:ied, jsd:jed, 1:nk)); cobalt%gross_prim_prod=0.0
    allocate(cobalt%nlg_diatoms(isd:ied, jsd:jed, 1:nk)); cobalt%nlg_diatoms=0.0
    allocate(cobalt%q_si_2_n_lg_diatoms(isd:ied, jsd:jed, 1:nk)); cobalt%q_si_2_n_lg_diatoms=0.0
    allocate(cobalt%b_alk(isd:ied, jsd:jed))              ; cobalt%b_alk=0.0
    allocate(cobalt%b_dic(isd:ied, jsd:jed))              ; cobalt%b_dic=0.0
    allocate(cobalt%b_fed(isd:ied, jsd:jed))              ; cobalt%b_fed=0.0
    allocate(cobalt%b_nh4(isd:ied, jsd:jed))              ; cobalt%b_nh4=0.0
    allocate(cobalt%b_no3(isd:ied, jsd:jed))              ; cobalt%b_no3=0.0
    allocate(cobalt%b_o2(isd:ied, jsd:jed))               ; cobalt%b_o2=0.0
    allocate(cobalt%b_po4(isd:ied, jsd:jed))              ; cobalt%b_po4=0.0
    allocate(cobalt%b_sio4(isd:ied, jsd:jed))             ; cobalt%b_sio4=0.0
    allocate(cobalt%pco2_csurf(isd:ied, jsd:jed))         ; cobalt%pco2_csurf=0.0
    allocate(cobalt%co2_csurf(isd:ied, jsd:jed))          ; cobalt%co2_csurf=0.0
    allocate(cobalt%co2_alpha(isd:ied, jsd:jed))          ; cobalt%co2_alpha=0.0
    allocate(cobalt%fcadet_arag_btm(isd:ied, jsd:jed))    ; cobalt%fcadet_arag_btm=0.0
    allocate(cobalt%fcadet_calc_btm(isd:ied, jsd:jed))    ; cobalt%fcadet_calc_btm=0.0
    allocate(cobalt%ffedet_btm(isd:ied, jsd:jed))         ; cobalt%ffedet_btm=0.0
    allocate(cobalt%flithdet_btm(isd:ied, jsd:jed))       ; cobalt%flithdet_btm=0.0
    allocate(cobalt%fpdet_btm(isd:ied, jsd:jed))          ; cobalt%fpdet_btm=0.0
    allocate(cobalt%fndet_btm(isd:ied, jsd:jed))          ; cobalt%fndet_btm=0.0
    allocate(cobalt%fsidet_btm(isd:ied, jsd:jed))         ; cobalt%fsidet_btm=0.0
    allocate(cobalt%fcased_burial(isd:ied, jsd:jed))      ; cobalt%fcased_burial=0.0
    allocate(cobalt%fcased_input(isd:ied, jsd:jed))       ; cobalt%fcased_input=0.0
    allocate(cobalt%fcased_redis(isd:ied, jsd:jed))       ; cobalt%fcased_redis=0.0
    allocate(cobalt%ffe_sed(isd:ied, jsd:jed))            ; cobalt%ffe_sed=0.0
    allocate(cobalt%fnfeso4red_sed(isd:ied, jsd:jed))     ; cobalt%fnfeso4red_sed=0.0
    allocate(cobalt%fno3denit_sed(isd:ied, jsd:jed))      ; cobalt%fno3denit_sed=0.0
    allocate(cobalt%fnoxic_sed(isd:ied, jsd:jed))         ; cobalt%fnoxic_sed=0.0
    allocate(cobalt%frac_burial(isd:ied, jsd:jed))        ; cobalt%frac_burial=0.0
    allocate(cobalt%fndet_burial(isd:ied, jsd:jed))       ; cobalt%fndet_burial=0.0
    allocate(cobalt%fpdet_burial(isd:ied, jsd:jed))       ; cobalt%fpdet_burial=0.0

    !
    ! allocate 100m integrated quantities
    !
    do n = 1, NUM_PHYTO
       allocate(phyto(n)%jprod_n_100(isd:ied,jsd:jed))      ; phyto(n)%jprod_n_100      = 0.0
       allocate(phyto(n)%jprod_n_new_100(isd:ied,jsd:jed))  ; phyto(n)%jprod_n_new_100  = 0.0
       allocate(phyto(n)%jzloss_n_100(isd:ied,jsd:jed))     ; phyto(n)%jzloss_n_100  = 0.0
       allocate(phyto(n)%jexuloss_n_100(isd:ied,jsd:jed))   ; phyto(n)%jexuloss_n_100  = 0.0
       allocate(phyto(n)%f_n_100(isd:ied,jsd:jed))          ; phyto(n)%f_n_100  = 0.0
    enddo
    allocate(phyto(DIAZO)%jprod_n_n2_100(isd:ied,jsd:jed)); phyto(DIAZO)%jprod_n_n2_100 = 0.0
    allocate(phyto(SMALL)%jvirloss_n_100(isd:ied,jsd:jed))  ; phyto(SMALL)%jvirloss_n_100 = 0.0
    allocate(phyto(SMALL)%jaggloss_n_100(isd:ied,jsd:jed))  ; phyto(SMALL)%jaggloss_n_100 = 0.0
    allocate(phyto(LARGE)%jaggloss_n_100(isd:ied,jsd:jed))  ; phyto(LARGE)%jaggloss_n_100 = 0.0
    allocate(cobalt%jprod_allphytos_100(isd:ied,jsd:jed))   ; cobalt%jprod_allphytos_100 = 0.0

   do n = 1, NUM_ZOO
       allocate(zoo(n)%jprod_n_100(isd:ied,jsd:jed))      ; zoo(n)%jprod_n_100      = 0.0
       allocate(zoo(n)%jingest_n_100(isd:ied,jsd:jed))    ; zoo(n)%jingest_n_100    = 0.0
       allocate(zoo(n)%jremin_n_100(isd:ied,jsd:jed))     ; zoo(n)%jremin_n_100     = 0.0
       allocate(zoo(n)%f_n_100(isd:ied,jsd:jed))          ; zoo(n)%f_n_100          = 0.0
    enddo

   do n = 1,2
       allocate(zoo(n)%jzloss_n_100(isd:ied,jsd:jed))     ; zoo(n)%jzloss_n_100     = 0.0
       allocate(zoo(n)%jprod_don_100(isd:ied,jsd:jed))    ; zoo(n)%jprod_don_100    = 0.0
   enddo

   do n = 2,3
       allocate(zoo(n)%jhploss_n_100(isd:ied,jsd:jed))    ; zoo(n)%jhploss_n_100     = 0.0
       allocate(zoo(n)%jprod_ndet_100(isd:ied,jsd:jed))   ; zoo(n)%jprod_ndet_100    = 0.0
   enddo

   allocate(cobalt%hp_jingest_n_100(isd:ied,jsd:jed))    ; cobalt%hp_jingest_n_100    = 0.0
   allocate(cobalt%hp_jremin_n_100(isd:ied,jsd:jed))     ; cobalt%hp_jremin_n_100     = 0.0
   allocate(cobalt%hp_jprod_ndet_100(isd:ied,jsd:jed))   ; cobalt%hp_jprod_ndet_100   = 0.0

   allocate(bact(1)%jprod_n_100(isd:ied,jsd:jed))   ; bact(1)%jprod_n_100 = 0.0
   allocate(bact(1)%jzloss_n_100(isd:ied,jsd:jed))  ; bact(1)%jzloss_n_100 = 0.0
   allocate(bact(1)%jvirloss_n_100(isd:ied,jsd:jed)); bact(1)%jvirloss_n_100 = 0.0
   allocate(bact(1)%jremin_n_100(isd:ied,jsd:jed))  ; bact(1)%jremin_n_100 = 0.0
   allocate(bact(1)%juptake_ldon_100(isd:ied,jsd:jed)) ; bact(1)%juptake_ldon_100 = 0.0
   allocate(bact(1)%f_n_100(isd:ied,jsd:jed))       ; bact(1)%f_n_100 = 0.0

   allocate(cobalt%jprod_lithdet_100(isd:ied,jsd:jed))      ; cobalt%jprod_lithdet_100 = 0.0
   allocate(cobalt%jprod_sidet_100(isd:ied,jsd:jed))        ; cobalt%jprod_sidet_100 = 0.0
   allocate(cobalt%jprod_cadet_calc_100(isd:ied,jsd:jed))   ; cobalt%jprod_cadet_calc_100 = 0.0
   allocate(cobalt%jprod_cadet_arag_100(isd:ied,jsd:jed))   ; cobalt%jprod_cadet_arag_100 = 0.0
   allocate(cobalt%jremin_ndet_100(isd:ied,jsd:jed))        ; cobalt%jremin_ndet_100 = 0.0
   allocate(cobalt%jprod_mesozoo_200(isd:ied,jsd:jed))      ; cobalt%jprod_mesozoo_200 = 0.0

   allocate(cobalt%f_ndet_100(isd:ied,jsd:jed))             ; cobalt%f_ndet_100 = 0.0
   allocate(cobalt%f_don_100(isd:ied,jsd:jed))              ; cobalt%f_don_100  = 0.0
   allocate(cobalt%f_silg_100(isd:ied,jsd:jed))             ; cobalt%f_silg_100 = 0.0
   allocate(cobalt%f_mesozoo_200(isd:ied,jsd:jed))          ; cobalt%f_mesozoo_200 = 0.0

   allocate(cobalt%fndet_100(isd:ied,jsd:jed))             ; cobalt%fndet_100 = 0.0
   allocate(cobalt%fpdet_100(isd:ied,jsd:jed))             ; cobalt%fpdet_100 = 0.0
   allocate(cobalt%fsidet_100(isd:ied,jsd:jed))            ; cobalt%fsidet_100 = 0.0
   allocate(cobalt%flithdet_100(isd:ied,jsd:jed))          ; cobalt%flithdet_100 = 0.0
   allocate(cobalt%fcadet_calc_100(isd:ied,jsd:jed))       ; cobalt%fcadet_calc_100 = 0.0
   allocate(cobalt%fcadet_arag_100(isd:ied,jsd:jed))       ; cobalt%fcadet_arag_100 = 0.0
   allocate(cobalt%ffedet_100(isd:ied,jsd:jed))            ; cobalt%ffedet_100 = 0.0

   allocate(cobalt%btm_temp(isd:ied,jsd:jed))            ; cobalt%btm_temp = 0.0
   allocate(cobalt%btm_o2(isd:ied,jsd:jed))            ; cobalt%btm_o2 = 0.0

   allocate(cobalt%o2min(isd:ied, jsd:jed)); cobalt%o2min=0.0
   allocate(cobalt%z_o2min(isd:ied, jsd:jed)); cobalt%z_o2min=0.0
   allocate(cobalt%z_sat_arag(isd:ied, jsd:jed)); cobalt%z_sat_arag=0.0
   allocate(cobalt%z_sat_calc(isd:ied, jsd:jed)); cobalt%z_sat_calc=0.0
   allocate(cobalt%mask_z_sat_arag(isd:ied, jsd:jed)); cobalt%mask_z_sat_arag = .FALSE.
   allocate(cobalt%mask_z_sat_calc(isd:ied, jsd:jed)); cobalt%mask_z_sat_calc = .FALSE.


  end subroutine user_allocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays
    integer n

    deallocate(cobalt%htotalhi,cobalt%htotallo)

    do n = 1, NUM_PHYTO
       deallocate(phyto(n)%def_fe)
       deallocate(phyto(n)%def_p)
       deallocate(phyto(n)%f_fe)
       deallocate(phyto(n)%f_n)
       deallocate(phyto(n)%felim)
       deallocate(phyto(n)%irrlim)
       deallocate(phyto(n)%jzloss_fe)
       deallocate(phyto(n)%jzloss_n)
       deallocate(phyto(n)%jzloss_p)
       deallocate(phyto(n)%jzloss_sio2)
       deallocate(phyto(n)%jaggloss_n) 
       deallocate(phyto(n)%jaggloss_p)
       deallocate(phyto(n)%jaggloss_fe)
       deallocate(phyto(n)%jaggloss_sio2)
       deallocate(phyto(n)%jvirloss_n)
       deallocate(phyto(n)%jvirloss_p)
       deallocate(phyto(n)%jvirloss_fe)
       deallocate(phyto(n)%jvirloss_sio2)
       deallocate(phyto(n)%jexuloss_n)
       deallocate(phyto(n)%jexuloss_p)
       deallocate(phyto(n)%jexuloss_fe)
       deallocate(phyto(n)%jhploss_fe)
       deallocate(phyto(n)%jhploss_n)
       deallocate(phyto(n)%jhploss_p)
       deallocate(phyto(n)%juptake_fe)
       deallocate(phyto(n)%juptake_nh4)
       deallocate(phyto(n)%juptake_no3)
       deallocate(phyto(n)%juptake_po4)
       deallocate(phyto(n)%jprod_n)
       deallocate(phyto(n)%liebig_lim)
       deallocate(phyto(n)%mu)
       deallocate(phyto(n)%po4lim)
       deallocate(phyto(n)%q_fe_2_n)
       deallocate(phyto(n)%q_p_2_n)
       deallocate(phyto(n)%q_si_2_n)
       deallocate(phyto(n)%theta)
    enddo
    do n = 2, NUM_PHYTO
       deallocate(phyto(n)%nh4lim)
       deallocate(phyto(n)%no3lim)
    enddo
    deallocate(phyto(DIAZO)%juptake_n2)
    deallocate(phyto(DIAZO)%o2lim)
    deallocate(phyto(LARGE)%juptake_sio4)
    deallocate(phyto(LARGE)%silim)

    ! bacteria
    deallocate(bact(1)%f_n)
    deallocate(bact(1)%jzloss_n)
    deallocate(bact(1)%jzloss_p)
    deallocate(bact(1)%jvirloss_n) 
    deallocate(bact(1)%jvirloss_p)   
    deallocate(bact(1)%jhploss_n)  
    deallocate(bact(1)%jhploss_p)
    deallocate(bact(1)%juptake_ldon)
    deallocate(bact(1)%juptake_ldop)
    deallocate(bact(1)%jprod_nh4)
    deallocate(bact(1)%jprod_po4)
    deallocate(bact(1)%jprod_n)
    deallocate(bact(1)%temp_lim)

    ! zooplankton
    do n = 1, NUM_ZOO
       deallocate(zoo(n)%f_n)
       deallocate(zoo(n)%jzloss_n)
       deallocate(zoo(n)%jzloss_p)
       deallocate(zoo(n)%jhploss_n)
       deallocate(zoo(n)%jhploss_p)
       deallocate(zoo(n)%jingest_n)
       deallocate(zoo(n)%jingest_p)
       deallocate(zoo(n)%jingest_sio2)
       deallocate(zoo(n)%jingest_fe)
       deallocate(zoo(n)%jprod_ndet)
       deallocate(zoo(n)%jprod_pdet)
       deallocate(zoo(n)%jprod_ldon)
       deallocate(zoo(n)%jprod_ldop)
       deallocate(zoo(n)%jprod_srdon)
       deallocate(zoo(n)%jprod_srdop)
       deallocate(zoo(n)%jprod_sldon)
       deallocate(zoo(n)%jprod_sldop)
       deallocate(zoo(n)%jprod_fedet)
       deallocate(zoo(n)%jprod_fed)
       deallocate(zoo(n)%jprod_sidet)
       deallocate(zoo(n)%jprod_sio4)
       deallocate(zoo(n)%jprod_po4)
       deallocate(zoo(n)%jprod_nh4)
       deallocate(zoo(n)%jprod_n)
       deallocate(zoo(n)%temp_lim)
    enddo

    deallocate(cobalt%hp_jingest_n)
    deallocate(cobalt%hp_jingest_p)
    deallocate(cobalt%hp_jingest_sio2)
    deallocate(cobalt%hp_jingest_fe)

    deallocate(&
    	cobalt%f_alk,&
    	cobalt%f_cadet_arag,&
        cobalt%f_cadet_calc,&
    	cobalt%f_dic,&
    	cobalt%f_fed,&
    	cobalt%f_fedet,&
    	cobalt%f_ldon,&
        cobalt%f_ldop,&
    	cobalt%f_lith,&
    	cobalt%f_lithdet,&
    	cobalt%f_ndet,&
    	cobalt%f_nh4,&
    	cobalt%f_no3,&
    	cobalt%f_o2,&
    	cobalt%f_pdet,&
    	cobalt%f_po4,&
        cobalt%f_srdon,&
        cobalt%f_srdop,&
    	cobalt%f_sldon,&
    	cobalt%f_sldop,&
    	cobalt%f_sidet,&
    	cobalt%f_silg,&
    	cobalt%f_sio4,&
    	cobalt%f_chl,&
    	cobalt%f_co3_ion,&
    	cobalt%f_htotal,&
    	cobalt%f_irr_mem,&
    	cobalt%f_cased,&
    	cobalt%f_cadet_arag_btf,&
        cobalt%f_cadet_calc_btf,&
    	cobalt%f_fedet_btf,&
    	cobalt%f_lithdet_btf,&
    	cobalt%f_ndet_btf,&
    	cobalt%f_pdet_btf,&
    	cobalt%f_sidet_btf,&
    	cobalt%jndi,&
    	cobalt%jnsm,&
    	cobalt%jnlg,&
    	cobalt%jnsmz,&
    	cobalt%jnmdz,&
    	cobalt%jnlgz,&
    	cobalt%jalk,&
    	cobalt%jcadet_arag,&
        cobalt%jcadet_calc,&
    	cobalt%jdic,&
    	cobalt%jfed,&
    	cobalt%jfedi,&
    	cobalt%jfelg,&
    	cobalt%jfesm,&
    	cobalt%jfedet,&
    	cobalt%jldon,&
        cobalt%jldop,&
    	cobalt%jlith,&
    	cobalt%jlithdet,&
    	cobalt%jndet,&
    	cobalt%jnh4,&
    	cobalt%jno3,&
    	cobalt%jo2,&
    	cobalt%jpdet,&
    	cobalt%jpo4,&
        cobalt%jsrdon,&
        cobalt%jsrdop,&
    	cobalt%jsldon,&
    	cobalt%jsldop,&
    	cobalt%jsidet,&
    	cobalt%jsilg,&
    	cobalt%jsio4,&
    	cobalt%jprod_ndet,&
    	cobalt%jprod_pdet,&
        cobalt%jprod_srdon,&
    	cobalt%jprod_sldon,&
    	cobalt%jprod_ldon,&
        cobalt%jprod_srdop,&
    	cobalt%jprod_sldop,&
    	cobalt%jprod_ldop,&
    	cobalt%jprod_fedet,&
    	cobalt%jprod_fed,&
    	cobalt%jprod_sidet,&
        cobalt%jprod_sio4,&
    	cobalt%jprod_lithdet,&
    	cobalt%jprod_cadet_arag,&
        cobalt%jprod_cadet_calc,&
    	cobalt%jprod_nh4,&
    	cobalt%jprod_po4,&
    	cobalt%det_jzloss_n,&
    	cobalt%det_jzloss_p,&
    	cobalt%det_jzloss_fe,&
        cobalt%det_jzloss_si,&
    	cobalt%det_jhploss_n,&
    	cobalt%det_jhploss_p,&
    	cobalt%det_jhploss_fe,&
        cobalt%det_jhploss_si,&
    	cobalt%jdiss_cadet_arag,&
        cobalt%jdiss_cadet_calc,&
    	cobalt%jdiss_sidet,&
    	cobalt%jremin_ndet,&
    	cobalt%jremin_pdet,&
    	cobalt%jremin_fedet,&
    	cobalt%jfe_ads,&
    	cobalt%jfe_coast,&
        cobalt%kfe_eq_lig,&
    	cobalt%expkT,&
    	cobalt%hp_temp_lim,&
    	cobalt%irr_inst,&
    	cobalt%irr_mix,&
    	cobalt%jno3denit_wc,&
    	cobalt%jnitrif,&
        cobalt%total_filter_feeding,&
        cobalt%net_prim_prod,&
        cobalt%gross_prim_prod,&
        cobalt%nlg_diatoms,&
        cobalt%q_si_2_n_lg_diatoms,&
    	cobalt%zt,&
        cobalt%zm,&
    	cobalt%tot_layer_int_c,&
    	cobalt%tot_layer_int_fe,&
    	cobalt%tot_layer_int_n,&
    	cobalt%tot_layer_int_p,&
    	cobalt%tot_layer_int_si)

    deallocate(&
        cobalt%b_alk,&
    	cobalt%b_dic,&
    	cobalt%b_fed,&
    	cobalt%b_nh4,&
    	cobalt%b_no3,&
    	cobalt%b_o2,&
    	cobalt%b_po4,&
    	cobalt%b_sio4,&
    	cobalt%pco2_csurf,&
    	cobalt%co2_csurf,&
    	cobalt%co2_alpha,&
    	cobalt%fcadet_arag_btm,&
        cobalt%fcadet_calc_btm,&
    	cobalt%ffedet_btm,&
    	cobalt%flithdet_btm,&
    	cobalt%fpdet_btm,&
    	cobalt%fndet_btm,&
    	cobalt%fsidet_btm,&
    	cobalt%fcased_burial,&
    	cobalt%fcased_input,&
    	cobalt%fcased_redis,&
    	cobalt%ffe_sed,&
    	cobalt%fnfeso4red_sed,&
    	cobalt%fno3denit_sed, &
        cobalt%fnoxic_sed, &
        cobalt%frac_burial,&
        cobalt%fndet_burial,&
        cobalt%fpdet_burial)

    do n = 1, NUM_PHYTO
       deallocate(phyto(n)%jprod_n_100)
       deallocate(phyto(n)%jprod_n_new_100)
       deallocate(phyto(n)%jzloss_n_100)
       deallocate(phyto(n)%jexuloss_n_100)
       deallocate(phyto(n)%f_n_100)
    enddo
    deallocate(phyto(DIAZO)%jprod_n_n2_100)
    deallocate(phyto(SMALL)%jvirloss_n_100)
    deallocate(phyto(SMALL)%jaggloss_n_100)
    deallocate(phyto(LARGE)%jaggloss_n_100)
    deallocate(cobalt%jprod_allphytos_100)

    do n = 1, NUM_ZOO
       deallocate(zoo(n)%jprod_n_100)
       deallocate(zoo(n)%jingest_n_100)
       deallocate(zoo(n)%jremin_n_100)
       deallocate(zoo(n)%f_n_100)
    enddo

    do n = 1,2
       deallocate(zoo(n)%jzloss_n_100)
       deallocate(zoo(n)%jprod_don_100)
    enddo

    do n = 2,3
       deallocate(zoo(n)%jhploss_n_100)
       deallocate(zoo(n)%jprod_ndet_100)
    enddo

    deallocate(cobalt%hp_jingest_n_100) 
    deallocate(cobalt%hp_jremin_n_100)  
    deallocate(cobalt%hp_jprod_ndet_100)  

    deallocate(bact(1)%jprod_n_100)
    deallocate(bact(1)%jzloss_n_100)
    deallocate(bact(1)%jvirloss_n_100)
    deallocate(bact(1)%jremin_n_100)
    deallocate(bact(1)%juptake_ldon_100)
    deallocate(bact(1)%f_n_100)

    deallocate(cobalt%jprod_lithdet_100)
    deallocate(cobalt%jprod_sidet_100)
    deallocate(cobalt%jprod_cadet_arag_100)
    deallocate(cobalt%jprod_cadet_calc_100)
    deallocate(cobalt%jremin_ndet_100)
    deallocate(cobalt%jprod_mesozoo_200)

    deallocate(cobalt%f_ndet_100)
    deallocate(cobalt%f_don_100)
    deallocate(cobalt%f_silg_100)
    deallocate(cobalt%f_mesozoo_200)

    deallocate(cobalt%fndet_100)
    deallocate(cobalt%fpdet_100)
    deallocate(cobalt%fsidet_100)
    deallocate(cobalt%flithdet_100)
    deallocate(cobalt%fcadet_calc_100)
    deallocate(cobalt%fcadet_arag_100)
    deallocate(cobalt%ffedet_100)
   
    deallocate(cobalt%btm_temp)
    deallocate(cobalt%btm_o2)

    deallocate(cobalt%o2min)
    deallocate(cobalt%z_o2min)
    deallocate(cobalt%z_sat_arag)
    deallocate(cobalt%z_sat_calc)
    deallocate(cobalt%mask_z_sat_arag)
    deallocate(cobalt%mask_z_sat_calc)

  end subroutine user_deallocate_arrays


end module generic_COBALT
