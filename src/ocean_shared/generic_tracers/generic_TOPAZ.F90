!----------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
! </CONTACT>
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
! This module contains the generic version of TOPAZ Tracers and their biogeochemistry.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
! The biogeochemistry calculations in this module are ported from GOLD_OCEAN_TOPAZ.F90
! released in omsk_2008_03 
! </OVERVIEW>
!<DESCRIPTION>
!       Phytoplankton Biogeochemistry: Includes an explicit ecological model
!       including three phytoplankton groups (small, large/diatoms and
!       diazotrophs), growth limitation by light, temperature and a suite of
!       nutrients including nitrate, ammonia, phosphate, iron and silicate,
!       dissolved inorganic carbon, alkalinity, two kinds of dissolved organic
!       material, O2, nitrogen fixation and denitrification. CO2 gas exchange
!       is function of the biologically and physically forced solubility.
!       Additionally, changes in the vertical distribution of phytoplankton
!       affect heat absorption with climate feedbacks. Food web processing in
!       the euphotic zone and remineralization/dissolution through the ocean
!       interior are handled as in Dunne et al. (in prep).  CO2 and O2
!       equilibria and gas exchange follow OCMIP2 protocols.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Biotic/HOWTO-Biotic.html
! </REFERENCE>
!
! <REFERENCE>
! Press, W. H., S. A. Teukosky, W. T. Vetterling, B. P. Flannery, 1992. 
! Numerical Recipes in FORTRAN, Second Edition, Cambridge University Press. 
! </REFERENCE>
!
! <REFERENCE>
! Enting, I.G., T. M. L. Wigley, M. Heimann, 1994. Future Emissions 
! and concentrations of carbon dioxide: key ocean / atmosphere / 
! land analyses, CSIRO Aust. Div. Atmos. Res. Tech. Pap. No. 31, 
! 118 pp.
! </REFERENCE>
!
! <DEVELOPER_NOTES>
! nnz: 
! A. Reproducing original GOLD TOPAZ results of GOLD_OCEAN_TOPAZ.F90
! If we bypass the update_from_source routine, then the tracers in this module reproduce the corresponding ones
! in non-generic module GOLD_OCEAN_TOPAZ.F90 with branch tag perth_gas_fluxes_nnz, 
! EXCEPT for tracers HTOTAL, CO3_ION  and DIC which still differ although only by roundoff errors. 
! This points to the slight differences between the FMS_ocmip2_co2calc subroutine in this module and 
! GOLD_ocmip2_co2calc used by GOLD_OCEAN_TOPAZ.F90 to calculate the HTOTAL and CO2 alpha and csurf.
! Particularly the difference is in CO2_FLUX_CSURF_OCN which depends on HTOTAL and not in CO2_FLUX_ALPHA_OCN.
! This difference in CO2_FLUX_CSURF_OCN gives rise to the difference in the final DIC concentration between
! the generic and non-generic TOPAZ.
!
! The reproducing of non-generic tracers in this case is sufficient evidence for the consistency of vertical diffusion
! (tracer_vertdiff) routine used in the two cases.
!
! If we include the update_from_source routine, tracers in this module no longer reproduce those 
! in GOLD_OCEAN_TOPAZ.F90 although they are within less than 1% of them. This is because the code
! path is sufficiently different in generic and non-generic cases if the update_from_source routine is included.
! 
! B. Reproducing original MOM TOPAZ results of ocean_topaz.F90
! Reproducing MOM TOPAZ results is generally less accurate since the whole vertical diffusion scheme is being
! done at a different point in the generic code path. 
! Preliminary tests with the simple "tester" prognostic tracer shows that the final results for generic vs. nongeneric
! code path for such simple non-intercating tracer is within less than 0.1% of each other.
!
! </DEVELOPER_NOTES>
! </INFO>
!----------------------------------------------------------------

module generic_TOPAZ

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len, fm_path_name_len
  use mpp_mod,           only: mpp_error, NOTE, WARNING, FATAL, stdout, mpp_chksum
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
  use g_tracer_utils, only : g_diag_type, g_diag_field_add

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, FMS_ocmip2_co2calc_old, CO2_dope_vector

  implicit none ; private
!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: generic_TOPAZ.F90,v 20.0 2013/12/14 00:18:07 fms Exp $'
  character(len=128) :: tag = '$Name: tikal $'
!-----------------------------------------------------------------------

  character(len=fm_string_len), parameter :: mod_name       = 'generic_TOPAZ'
  character(len=fm_string_len), parameter :: package_name   = 'generic_topaz'

  public do_generic_TOPAZ
  public generic_TOPAZ_register
  public generic_TOPAZ_init
  public generic_TOPAZ_register_diag
  public generic_TOPAZ_update_from_coupler
  public generic_TOPAZ_update_from_source
  public generic_TOPAZ_update_from_bottom
  public generic_TOPAZ_set_boundary_values
  public generic_TOPAZ_end

  !The following logical for using this module is overwritten 
  ! generic_tracer_nml namelist
  logical, save :: do_generic_TOPAZ = .false.

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30
  real, parameter :: missing_value1=-1.0e+20
  real, parameter :: missing_value_diag=-1.0e+10

  integer :: id_rho_dzt          = -1
  !
  !The following two types contain all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type phytoplankton
     real :: alpha,   &
          fdet0,         &
          fe_2_n_max,    &
          fe_2_n_norm,   &
          fe_2_n_static, &
          k_fe_2_n,      &
          k_fed,         &
          k_nh4,         &
          k_no3,         &
          k_p_2_n,       &
          k_po4,         &
          k_sio4,        &
          p_2_n_assem,   &
          p_2_n_static,  &
          P_C_max,       &
          plast_2_chl,   &
          r_nfix,        &
          r_other_min,   &
          r_uptake_max,  &
          si_2_n_max,    &
          si_2_n_static, &
          thetamax     
     real, ALLOCATABLE, dimension(:,:)  :: &
          jgraz_n_100,  & ! Nitrogen grazing integral in upper 100m
          jprod_n_100,  & ! Nitrogen production integral in upper 100m
          jprod_sio4_100  ! Silicon production integral in upper 100m

     real, pointer, dimension(:,:,:)  :: &
          chl         , & ! Phyto. Chl
          jgraz_n     , & ! Nitrogen grazing layer integral
          jprod_n     , & ! Total Nitrogen production layer integral
          jprod_sio4      ! Silicon production layer integral

     real, ALLOCATABLE, dimension(:,:,:)  :: &
          def_fe      , & ! Fe Deficiency
          def_p       , & ! P Deficiency
          f_fe        , & ! Phytoplankton Fe concentration
          f_n         , & ! Phytoplankton Nitrogen concentration
          f_p         , & ! Phytoplankton Phosphorus concentration
          felim       , & ! Fed Limitation
          irrlim      , & ! Light Limitation
          jgraz_fe    , & ! Fe grazing layer integral
          jgraz_sio2  , & ! Silicon grazing layer integral
          jprod_n2    , & ! Nitrogen fixation layer integral
          jprod_fe    , & ! Fe production layer integral
          jprod_nh4   , & ! NH4 production layer integral
          jprod_no3   , & ! NO3 production layer integral
          jprod_po4   , & ! PO4 production layer integral
          liebig_lim  , & ! Overall nutrient limitation
          mu          , & ! Overall growth rate
          nh4lim      , & ! Ammonia limitation
          no3lim      , & ! Nitrate limitation
          po4lim      , & ! Phosphate limitation
          q_fe_2_n    , & ! Fe:N ratio
          q_p_2_n     , & ! P:N ratio
          q_p_2_n_opt , & ! Optimal P:N ratio
          silim       , & ! SiO4 limitation
          theta           ! Chl:C ratio
     integer ::            &
          id_chl        = -1, & ! Diag id for Phyto. Chl
          id_def_fe     = -1, & ! Diag id for Phyto. Fe Deficiency
          id_def_p      = -1, & ! Diag id for Phyto. P Deficiency
          id_felim      = -1, & ! Diag id for Phyto. Fed Limitation
          id_irrlim     = -1, & ! Diag id for Phyto. Light Limitation
          id_jgraz_fe   = -1, & ! Diag id for iron grazing layer integral
          id_jgraz_n    = -1, & ! Diag id for nitrogen grazing layer integral
          id_jgraz_n_100= -1, & ! Diag id for nitrogen grazing integral in upper 100m
          id_jgraz_sio2 = -1, & ! Diag id for silicon grazing layer integral
          id_jprod_fe   = -1, & ! Diag id for phyto. Fed production layer integral
          id_jprod_n    = -1, & ! Diag id for total Nitrogen production layer integral
          id_jprod_n_100= -1, & ! Diag id for total Nitrogen production layer integral
          id_jprod_n2   = -1, & ! Diag id for Nitrogen fixation layer integral
          id_jprod_nh4  = -1, & ! Diag id for phyto. NH4 production layer integral
          id_jprod_no3  = -1, & ! Diag id for phyto. NO3 production layer integral
          id_jprod_po4  = -1, & ! Diag id for phyto. PO4 production layer integral
          id_jprod_sio4 = -1, & ! Diag id for phyto. SiO4 production layer integral
          id_jprod_sio4_100 = -1, & ! Diag id for phyto. SiO4 production integral in upper 100m
          id_liebig_lim = -1, & ! Diag id for Overall nutrient limitation
          id_mu         = -1, & ! Diag id for Overall growth rate
          id_nh4lim     = -1, & ! Diag id for Ammonia Limitation of Phyto
          id_no3lim     = -1, & ! Diag id for Nitrate Limitation of Phyto
          id_po4lim     = -1, & ! Diag id for Phosphate Limitation of Phyto
          id_q_fe_2_n   = -1, & ! Diag id for Fe:N ratio
          id_q_p_2_n    = -1, & ! Diag id for P:N ratio
          id_q_p_2_n_opt= -1, & ! Diag id for Optimal P:N ratio
          id_sfc_chl    = -1, & ! Diag id for Surface Phyto. Chl
          id_silim      = -1, & ! Diag id for SiO4 Limitation of Phyto
          id_theta      = -1
  end type phytoplankton

  integer, parameter :: NUM_PHYTO  = 3
  !
  ! Array allocations and flux calculations assume that phyto(1) is the
  ! only phytoplankton group cabable of nitrogen uptake by N2 fixation while phyto(2:NUM_PHYTO) 
  ! are only cabable of nitrgen uptake by NH4 and NO3 uptake
  !
  integer, parameter :: DIAZO      = 1
  integer, parameter :: LARGE      = 2
  integer, parameter :: SMALL      = 3
  type(phytoplankton), dimension(NUM_PHYTO), save  :: phyto


  type generic_TOPAZ_type

     logical  ::       &
          ca_2_n_static,    &                  ! If Ca:N is fixed in heterotrophs
          do_extra_diag_calcs,&                ! perform extra calculations for online diagnostics
          fe_2_n_static,    &                  ! If Fe:N is fixed in phytoplankton
          fe_ballast_assoc, &                  ! If iron scavenging is associated with ballast
          init,             &                  ! If tracers should be initializated
          p_2_n_static,     &                  ! If P:N is fixed in phytoplankton
          reproduce_esm2,   &                  ! Do not add changes post-CMIP5 version of ESM2
          si_2_n_static                        ! If Si:N is fixed in phytoplankton

     real  ::          &
          atm_co2_flux,     &
          bio_tau,          &                  ! Miscellaneous
          c_2_n,            &                  ! Stoichiometry
          ca_2_n_arag,      &                  ! Aragonite CaCO3 formation in heterotrophs
          ca_2_n_calc,      &                  ! Calcite CaCO3 formation in heterotrophs
          ca_2_n_het_static,&                  ! CaCO3 formation if fixed
          caco3_sat_max,    &                  ! Maximum saturation state effect
          fe_2_n_sed,       &                  ! Iron
          fe_coast,         &                  ! Iron
          felig_2_don,      &                  ! Iron
          felig_bkg ,       &                  ! Iron
          gamma_cadet_arag, &                  ! Aragonite Dissolution
          gamma_cadet_calc, &                  ! Calcite Dissolution
          gamma_irr_mem,    &                  ! Photosynthesis
          gamma_ldon,       &                  ! Dissolved Organic Material
          gamma_ndet,       &                  ! Remineralization
          gamma_nhet,       &                  ! Grazing
          gamma_nitrif,     &                  ! Miscellaneous
          gamma_sidet,      &                  ! Dissolution
          gamma_sdon,       &                  ! Dissolved Organic Material
          gamma_sdop,       &                  ! Dissolved Organic Material
          irr_inhibit,      &                  ! Miscellaneous
          k_diss_sio2,      &                  ! Grazing
          k_n_inhib_di,     &                  ! Photosynthesis
          k_o2,             &                  ! Stoichiometry
          kappa_eppley,     &                  ! Photosynthesis
          kappa_remin,      &                  ! Grazing
          kfe_2nd_order,    &                  ! Iron
          kfe_bal,          &
          kfe_des,          &                  ! Iron
          kfe_eq_lig,       &                  ! Iron
          kfe_org,          &
          k_lith,           &                  ! Miscellaneous
          lambda0,          &                  ! Grazing
          alk_2_n_denit,    &                  ! Stoichiometry
          mass_2_n,         &                  ! Stoichiometry
          n_2_n_denit,      &                  ! Stoichiometry
          o2_min,           &                  ! Grazing
          o2_2_c,           &                  ! Stoichiometry
          o2_2_nfix,        &                  ! Stoichiometry
          o2_2_nh4,         &                  ! Stoichiometry
          o2_2_no3,         &                  ! Stoichiometry
          o2_2_nitrif,      &                  ! Stoichiometry
          o2_inhib_di_pow,  &                  ! Photosynthesis
          o2_inhib_di_sat,  &                  ! Photosynthesis
          p_2_n_photo,      &                  ! P:N of photosynthesis machinery (Chloroplasts)
          p_2_n_RKR,        &                  ! P:N of standard Stoichiometry
          p_2_n_uptake,     &                  ! P:N of uptake machinery
          P_star,           &                  ! Grazing
          phi_nhet,         &                  ! Grazing
          phi_sdon,         &                  ! Dissolved Organic Material
          phi_sdop,         &                  ! Dissolved Organic Material
          phi_ldon,         &                  ! Dissolved Organic Material
          phi_lith,         &                  ! Miscellaneous
          phyto_min,        &                  ! Grazing
          plast_2_chl,      &                  ! Chloroplast:Chlorophyll conversion
          q_si_2_n_diss,    &                  ! Stoichiometric role in SiO2 dissolution
          r_bio_tau,        &                  ! Miscellaneous
          rpcaco3,          &                  ! Grazing
          rplith,           &                  ! Grazing
          rpsio2,           &                  ! Grazing
          thetamin,         &                  ! Photosynthesis
          thetamin_nolim,   &                  ! Photosynthesis
          wsink,            &                  ! Sinking
          z_sed,            &                  ! Sediment thickness
          zeta                                 ! Photosynthesis

     real, dimension(3)                    :: total_atm_co2

     real    :: htotal_scale_lo, htotal_scale_hi, htotal_in
     real    :: Rho_0, a_0, a_1, a_2, a_3, a_4, a_5, b_0, b_1, b_2, b_3, c_0
     real    :: a1_co2, a2_co2, a3_co2, a4_co2, a1_o2, a2_o2, a3_o2, a4_o2

     logical, dimension(:,:), ALLOCATABLE ::  &
          mask_z_sat_arag,&
          mask_z_sat_calc

     real, dimension(:,:,:), pointer ::  &
          chl_lg_diatoms,&
          chl_lg_nondiatoms,&
          co3_sol_arag,&
          co3_sol_calc,&
          f_cased,&
          f_chl,&
          f_co3_ion,&
          f_htotal,&
          f_irr_mem,&
          irr_inst,&
          jalk,&
          jalk_plus_btm,&
          jcadet_arag,&
          jcadet_calc,&
          jdic,&
          jdic_plus_btm,&
          jdin,&
          jdin_plus_btm,&
          jfe_ads,&
          jfe_des,&
          jfed,&
          jfed_plus_btm,&
          jfedet,&
          jgraz_ntot,&
          jndet,&
          jnh4_plus_btm,&
          jnitrif,&
          jno3_plus_btm,&
          jo2_plus_btm,&
          jpo4,&
          jpo4_plus_btm,&
          jprod_cadet_arag,&
          jprod_cadet_calc,&
          jprod_fetot,&
          jprod_ndet,&
          jprod_nlg_diatoms,&
          jprod_nlg_nondiatoms,&
          jprod_no3tot,&
          jprod_ntot,&
          jprod_ptot,&
          jsio4,&
          jsio4_plus_btm,&
          nLg_diatoms,&
          nLg_nondiatoms,&
          tot_bfe,&
          tot_bsi,&
          tot_doc,&
          tot_pon,&
          tot_pop,&
          tot_phyton,&
          tot_phytop,&
          tot_phytofe

     real, dimension(:,:,:), ALLOCATABLE ::  &
          f_alk,&
          f_cadet_arag,&
          f_cadet_arag_btf,&
          f_cadet_calc,&
          f_cadet_calc_btf,&
          f_dic,&
          f_fed,&
          f_fedet,&
          f_ldon,&
          f_lith,&
          f_lithdet,&
          f_lithdet_btf,&
          f_ndet,&
          f_ndet_btf,&
          f_nh4,&
          f_nhet,&
          f_no3,&
          f_o2,&
          f_pdet,&
          f_pdet_btf,&
          f_po4,&
          f_sdon,&
          f_sdop,&
          f_sidet,&
          f_sidet_btf,&
          f_silg,&
          f_sio4,&
          expkT,&
          frac_det_prod,&
          irr_mix,&
          jdiss_sio2,&
          jfe_graz,&
          jfe_coast,&
          jldon,&
          jnh4,&
          jnh4_graz,&
          jnhet,&
          jno3,&
          jno3denit_wc,&
          jo2,&
          jpdet,&
          jpo4_graz,&
          jprod_fedet,&
          jprod_lithdet,&
          jprod_nhet,&
          jprod_pdet,&
          jsdon,&
          jsdop,&
          jsidet,&
          omega_arag,&
          omega_calc,&
          q_si_2_n_Lg_diatoms,&
          tot_layer_int_c,&
          tot_layer_int_fe,&
          tot_layer_int_n,&
          tot_layer_int_p,&
          tot_layer_int_si,&
          zt

     real, dimension(:,:), ALLOCATABLE :: &
          b_alk,b_dic,b_fed,b_nh4,b_no3,b_o2,b_po4,b_sio4,&
          co2_alpha,&
          co2_csurf,&
          f_alk_int_100,&
          f_dic_int_100,&
          f_din_int_100,&
          f_fed_int_100,&
          f_po4_int_100,&
          f_sio4_int_100,&
          fca_sed_sfc,&
          fcadet_arag_100,&
          fcadet_arag_btm,&
          fcadet_calc_100,&
          fcadet_calc_btm,&
          fcased_burial,&
          fcased_redis,&
          ffe_sed,&
          ffedet_100,&
          ffedet_btm,&
          flithdet_100,&
          flithdet_btm,&
          fndet_100,&
          fndet_btm,&
          fnfeso4red_sed,&
          fno3denit_sed,&
          fno3denit_tot,&
          fnoxic_sed,&
          fpdet_100,&
          fpdet_btm,&
          fsidet_100,&
          fsidet_btm,&
          htotallo, htotalhi,&
          jalk_100,&
          jdic_100,&
          jdin_100,&
          jfed_100,&
          jnitrif_100,&
          jpo4_100,&
          jprod_cadet_arag_100,&
          jprod_cadet_calc_100,&
          jprod_fetot_100,&
          jprod_nlg_diatoms_100,&
          jprod_nlg_nondiatoms_100,&
          jprod_no3tot_100,&
          jprod_ntot_100,&
          jprod_ptot_100,&
          jsio4_100,&
          nongas_source_c,&
          nongas_source_fe,&
          nongas_source_n,&
          o2min,&
          pco2_csurf,&
          wc_vert_int_c,&
          wc_vert_int_dic,&
          wc_vert_int_jfe_coast,&
          wc_vert_int_jno3denit,&
          wc_vert_int_nfix,&
          z_o2min,&
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
          p_lith,&
          p_lithdet,&          
          p_ndet,&
          p_ndi,&
          p_nlg,&
          p_nsm,&
          p_nh4,&
          p_nhet,&
          p_no3,&
          p_o2,&
          p_pdet,&
          p_pdi,&
          p_plg,&
          p_po4,&
          p_psm,&
          p_sdon,&
          p_sdop,&
          p_sidet,&
          p_silg,&
          p_sio4

     integer :: nkml
     character(len=fm_string_len)          :: file
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file

     integer               ::          &
          id_cased_diag    = -1,       & ! Sediment CaCO3 integral via send diag
          id_chl_diag      = -1,       & ! Diag id for Phyto. Chl via send data
          id_chl_lg_diatoms= -1,       & ! large diatom Chlorophyll
          id_chl_lg_nondiatoms = -1,   & ! large non-diatom Chlorophyll
          id_co3_ion       = -1,       & ! Carbonate Ion
          id_co3_ion_diag  = -1,       & ! Carbonate Ion via send diag
          id_co3_sol_arag  = -1,       & ! Carbonate Ion Solubility for Aragonite
          id_co3_sol_calc  = -1,       & ! Carbonate Ion Solubility for Calcite
          id_dep_dry_fed   = -1,       & ! Iron dry deposition
          id_dep_dry_nh4   = -1,       & ! Ammonia dry deposition
          id_dep_dry_no3   = -1,       & ! Nitrate dry deposition
          id_dep_wet_fed   = -1,       & ! Iron wet deposition
          id_dep_wet_nh4   = -1,       & ! Ammonia wet deposition
          id_dep_wet_no3   = -1,       & ! Nitrate wet deposition
          id_f_alk_int_100 = -1,       & ! Inventory of Alkalinity in the upper 100 m
          id_f_dic_int_100 = -1,       & ! Inventory of Dissolved Inorganic Carbon in the upper 100 m
          id_f_din_int_100 = -1,       & ! Inventory of Inorganic nitrogen in the upper 100 m
          id_f_fed_int_100 = -1,       & ! Inventory of Dissolved Iron in the upper 100 m
          id_f_po4_int_100 = -1,       & ! Inventory of Dissolved Phosphate in the upper 100 m
          id_f_sio4_int_100= -1,       & ! Inventory of Dissolved Silicate in the upper 100 m
          id_fca_sed_sfc   = -1,       & ! Net Calcite CaCO3 flux across sediment surface (fcadet_calc_btm-redissolution)
          id_fcased_int    = -1,       & ! Sediment CaCO3 integral
          id_fcadet_arag   = -1,       & ! Aragonite CaCO3 sinking flux
          id_fcadet_arag_100 = -1,     & ! Aragonite CaCO3 sinking flux integral in upper 100m
          id_fcadet_arag_btm = -1,     & ! Aragonite CaCO3 sinking flux at bottom
          id_fcadet_calc   = -1,       & ! Calcite CaCO3 sinking flux
          id_fcadet_calc_100 = -1,     & ! Calcite CaCO3 sinking flux integral in upper 100m
          id_fcadet_calc_btm = -1,     & ! Calcite CaCO3 sinking flux at bottom
          id_fcased_burial = -1,       & ! Calcite CaCO3 sinking flux permanent burial
          id_fcased_redis  = -1,       & ! Calcite CaCO3 redissolution after sinking flux burial
          id_ffe_sed       = -1,       & ! Sediment iron efflux
          id_ffedet        = -1,       & ! POFe sinking flux
          id_ffedet_100    = -1,       & ! POFe sinking flux integral in upper 100m
          id_ffedet_btm    = -1,       & ! Iron sinking flux burial
          id_flithdet      = -1,       & ! Lith sinking flux
          id_flithdet_100  = -1,       & ! Lithogenic sinking flux at 100m
          id_flithdet_btm  = -1,       & ! Lithogenic sinking flux burial
          id_fndet         = -1,       & ! PON sinking flux
          id_fndet_100     = -1,       & ! PON sinking flux integral in upper 100m
          id_fndet_btm     = -1,       & ! Nitrogen sinking flux at bottom
          id_fnfeso4red_sed= -1,       & ! Sediment NDET Fe and SO4 reducing degradation flux
          id_fno3denit_sed = -1,       & ! Sediment NO3 Denitrification flux
          id_fno3denit_tot = -1,       & ! Total NO3 Denitrification flux
          id_fnoxic_sed    = -1,       & ! Sediment NDET oxic degradation flux
          id_fpdet         = -1,       & ! POP sinking flux
          id_fpdet_100     = -1,       & ! POP sinking flux integral in upper 100m
          id_fpdet_btm     = -1,       & ! Phosphorus sinking flux at bottom
          id_fsidet        = -1,       & ! Silicon sinking flux
          id_fsidet_100    = -1,       & ! Silicon sinking flux integral in upper 100m
          id_fsidet_btm    = -1,       & ! Silicon sinking flux at bottom
          id_htotal        = -1,       & ! H+ ion concentration
          id_htotal_diag   = -1,       & ! H+ ion concentration via send data
          id_irr_inst      = -1,       & ! Instantaneous Light
          id_irr_mem_diag  = -1,       & ! Irradiance Memory via send diag
          id_irr_mix       = -1,       & ! Light in mixing layer
          id_jalk          = -1,       & ! Alk source layer integral
          id_jalk_plus_btm = -1,       & ! Alk source layer integral
          id_jalk_100      = -1,       & ! Alk source integral in upper 100m
          id_jcadet_arag   = -1,       & ! Aragonite CaCO3 change layer integral
          id_jcadet_calc   = -1,       & ! Calcite CaCO3 change layer integral
          id_jdic          = -1,       & ! Dissolved Inorganic Carbon source layer integral
          id_jdic_plus_btm = -1,       & ! Dissolved Inorganic Carbon source layer integral
          id_jdic_100      = -1,       & ! Dissolved Inorganic Carbon source layer integral in upper 100m
          id_jdin          = -1,       & ! Dissolved Inorganic Nitrogen (NO3+NH4) source layer integral
          id_jdin_plus_btm = -1,       & ! Dissolved Inorganic Nitrogen (NO3+NH4) source layer integral
          id_jdin_100      = -1,       & ! Dissolved Inorganic Nitrogen (NO3+NH4) source integral in upper 100m
          id_jdiss_sio2    = -1,       & ! SiO2 Dissolution during grazing layer integral
          id_jfe_ads       = -1,       & ! Iron adsorption layer integral
          id_jfe_des       = -1,       & ! Iron desorption layer integral
          id_jfe_graz      = -1,       & ! Dissolved iron grazing source layer integral
          id_jgraz_ntot    = -1,       & ! Total Nitrogen grazing source layer integral
          id_jgraz_ntot_100= -1,       & ! Total Nitrogen grazing source integral in upper 100m
          id_jfe_coast     = -1,       & ! Coastal iron efflux layer integral
          id_jfed          = -1,       & ! Dissolved iron source layer integral
          id_jfed_plus_btm = -1,       & ! Dissolved iron source layer integral
          id_jfed_100      = -1,       & ! Dissolved iron source layer integral in upper 100m
          id_jfedet        = -1,       & ! Loss of sinking iron layer integral
          id_jldon         = -1,       & ! Labile DON source layer integral
          id_jndet         = -1,       & ! Loss of sinking nitrogen layer integral
          id_jnh4          = -1,       & ! NH4 source layer integral
          id_jnh4_graz     = -1,       & ! NH4 grazing source layer integral
          id_jnh4_plus_btm = -1,       & ! NH4 source layer integral
          id_jnhet         = -1,       & ! Heterotrophic Nitrogen remineralization layer integral
          id_jnitrif       = -1,       & ! Nitrification layer integral
          id_jnitrif_100   = -1,       & ! Nitrification integral in upper 100m
          id_jno3          = -1,       & ! NO3 source layer integral
          id_jno3_plus_btm = -1,       & ! NO3 source layer integral
          id_jno3denit_wc  = -1,       & ! Water column Denitrification layer integral
          id_jo2           = -1,       & ! O2 source layer integral
          id_jo2_plus_btm  = -1,       & ! O2 source layer integral
          id_jpdet         = -1,       & ! Loss of sinking phosphorus layer integral
          id_jpo4          = -1,       & ! Dissolved PO4 source layer integral
          id_jpo4_plus_btm = -1,       & ! Dissolved PO4 source layer integral
          id_jpo4_100      = -1,       & ! Dissolved PO4 source layer integral in upper 100m
          id_jpo4_graz     = -1,       & ! PO4 source from grazing layer integral
          id_jprod_cadet_arag = -1,    & ! Aragonite CaCO3 production layer integral  
          id_jprod_cadet_arag_100 = -1,& ! Aragonite CaCO3 production integral in upper 100m  
          id_jprod_cadet_calc = -1,    & ! Calcite CaCO3 production layer integral  
          id_jprod_cadet_calc_100 = -1,& ! Calcite CaCO3 production integral in upper 100m  
          id_jprod_fetot   = -1,       & ! Total Iron production layer integral
          id_jprod_fetot_100 = -1,     & ! Total Iron production integral in upper 100m
          id_jprod_fedet   = -1,       & ! Detrital Iron production layer integral
          id_jprod_lithdet = -1,       & ! Lithogenic removal to sinking layer integral
          id_jprod_nlg_diatoms = -1,   & ! Large Diatom Nitrogen production layer integral
          id_jprod_nlg_diatoms_100 = -1, & ! Large Diatom Nitrogen production integral in upper 100m
          id_jprod_nlg_nondiatoms = -1, & ! Large nonDiatom Nitrogen production layer integral
          id_jprod_nlg_nondiatoms_100 = -1, & ! Large nonDiatom Nitrogen production integral in upper 100m
          id_jprod_ntot    = -1,       & ! Total Nitrogen production layer integral
          id_jprod_ntot_100= -1,       & ! Total Nitrogen production integral in upper 100m
          id_jprod_nhet    = -1,       & ! Heterotrophic Nitrogen Production layer integral
          id_jprod_ndet    = -1,       & ! Detrital Nitrogen production layer integral
          id_jprod_no3tot  = -1,       & ! NO3 production layer integral
          id_jprod_no3tot_100= -1,     & ! NO3 production integral in upper 100m
          id_jprod_ptot    = -1,       & ! Total Phosphorus production layer integral
          id_jprod_ptot_100= -1,       & ! Total Phosphorus production integral in upper 100m
          id_jprod_pdet    = -1,       & ! Detrital Phosphorus production layer integral
          id_jsdon         = -1,       & ! Semilabile DON source layer integral
          id_jsdop         = -1,       & ! Semilabile DOP source layer integral
          id_jsidet        = -1,       & ! Loss of sinking silicon layer integral
          id_jsio4         = -1,       & ! Dissolved SiO4 source layer integral
          id_jsio4_plus_btm = -1,      & ! Dissolved SiO4 source layer integral
          id_jsio4_100     = -1,       & ! Dissolved SiO4 source layer integral in upper 100m
          id_nLg_diatoms   = -1,       & ! Large diatom nitrogen
          id_nLg_nondiatoms= -1,       & ! Large nondiatom nitrogen
          id_no3_in_source = -1,       & ! Nitrate Concentration used in source routine
          id_nongas_source_c= -1,      & ! Sum of all non-gas sources of Carbon to the ocean
          id_nongas_source_fe= -1,     & ! Sum of all non-gas sources of Iron to the ocean
          id_nongas_source_n= -1,      & ! Sum of all non-gas sources of Nitrogen to the ocean
          id_o2min         = -1,       & ! minimum oxygen
          id_omega_arag    = -1,       & ! Carbonate Ion Saturation state for Aragonite
          id_omega_calc    = -1,       & ! Carbonate Ion Saturation state for Calcite
          id_pco2surf      = -1,       & ! Oceanic pCO2
          id_sfc_alk       = -1,       & ! Surface Alkalinity
          id_sfc_cadet_arag= -1,       & ! Surface Detrital Aragonite
          id_sfc_cadet_calc= -1,       & ! Surface Detrital Calcite
          id_sfc_chl       = -1,       & ! Surface Chl
          id_sfc_chl_lg_diatoms = -1,  & ! Surface large diatom Chl
          id_sfc_chl_lg_nondiatoms = -1, & ! Surface large nondiatom Chl
          id_sfc_co3_ion       = -1,   & ! Surface Carbonate Ion
          id_sfc_co3_sol_arag  = -1,   & ! Surface Carbonate Ion Solubility for Aragonite
          id_sfc_co3_sol_calc  = -1,   & ! Surface Carbonate Ion Solubility for Calcite
          id_sfc_dic       = -1,       & ! Surface DIC 
          id_sfc_fed       = -1,       & ! Surface Fed 
          id_sfc_htotal    = -1,       & ! Surface Htotal
          id_sfc_irr       = -1,       & ! Surface Irr 
          id_sfc_irr_mem   = -1,       & ! Surface Irr_mem
          id_sfc_ndet      = -1,       & ! Surface Detrital Nitrogen
          id_sfc_ndi       = -1,       & ! Surface ndi
          id_sfc_nh4       = -1,       & ! Surface NH4
          id_sfc_nhet      = -1,       & ! Surface Heterotroph Nitrogen
          id_sfc_nlg       = -1,       & ! Surface nlg
          id_sfc_nlg_diatoms = -1,     & ! Surface nlg_diatoms
          id_sfc_nlg_nondiatoms = -1,  & ! Surface nlg_nondiatoms
          id_sfc_no3       = -1,       & ! Surface NO3
          id_sfc_nsm       = -1,       & ! Surface nsm
          id_sfc_o2        = -1,       & ! Surface O2
          id_sfc_po4       = -1,       & ! Surface PO4
          id_sfc_silg      = -1,       & ! Surface large phytoplankton Si
          id_sfc_sio4      = -1,       & ! Surface SiO4
          id_sfc_temp      = -1,       & ! Surface Temp
          id_sfc_tot_bfe   = -1,       & ! Surface total particulate biogenic Iron
          id_sfc_tot_bsi   = -1,       & ! Surface total particulate biogenic Silicon
          id_sfc_tot_doc   = -1,       & ! Surface total dissolved organic carbon
          id_sfc_tot_phytofe = -1,     & ! Surface total phytoplankton Iron
          id_sfc_tot_phyton  = -1,     & ! Surface total phytoplankton Nitrogen
          id_sfc_tot_phytop  = -1,     & ! Surface total phytoplankton Phosphorus
          id_sfc_tot_pon   = -1,       & ! Surface total particulate organic Nitrogen
          id_sfc_tot_pop   = -1,       & ! Surface total particulate organic Phosphorus
          id_tot_bfe = -1,             & ! Total Biogenic Iron (Phyto+Det)
          id_tot_bsi = -1,             & ! Total Biogenic Silica (Phyto+Det)
          id_tot_doc = -1,             & ! Total Dissolved Organic Carbon (refractory+labile+semilabile)
          id_tot_layer_int_c = -1,     & ! Total Carbon (DIC+OC+IC) boxwise layer integral
          id_tot_layer_int_fe = -1,    & ! Total Phosphorus (Fed+OFe) boxwise layer integral
          id_tot_layer_int_n = -1,     & ! Total Nitrogen (NO3+NH4+ON) boxwise layer integral
          id_tot_layer_int_p = -1,     & ! Total Phosphorus (PO4+OP) boxwise layer integral
          id_tot_layer_int_si = -1,    & ! Total Silicon (SiO4+SiO2) boxwise layer integral
          id_tot_pon = -1,             & ! Total Particulate Organic Nitrogen (Phyto+Nhet+Det)
          id_tot_pop = -1,             & ! Total Particulate Organic Phosphorus (Phyto+Nhet+Det)
          id_tot_phytofe   = -1,       & ! Total Iron in phytoplankton (Di+Lg+Sm)
          id_tot_phyton    = -1,       & ! Total Nitrogen in phytoplankton (Di+Lg+Sm)
          id_tot_phytop    = -1,       & ! Total Phosphorus in phytoplankton (Di+Lg+Sm)
          id_runoff_flux_alk = -1,     & ! Alkalinity runoff Flux to the ocean
          id_runoff_flux_dic = -1,     & ! DIC runoff Flux to the ocean
          id_runoff_flux_fed = -1,     & ! Iron runoff Flux to the ocean
          id_runoff_flux_lith = -1,    & ! Lithogenic runoff Flux to the ocean
          id_runoff_flux_nh4 = -1,     & ! Ammonia runoff Flux to the ocean
          id_runoff_flux_no3 = -1,     & ! Nitrate runoff Flux to the ocean
          id_wc_vert_int_c = -1,       & ! Total Carbon water column vertical integral
          id_wc_vert_int_dic = -1,     & ! DIC water column vertical integral
          id_wc_vert_int_jfe_coast=-1, & ! Coastal iron efflux vertical integral
          id_wc_vert_int_jno3denit=-1, & ! Water column Denitrification vertical integral
          id_wc_vert_int_nfix = -1,    & ! Nitrogen fixation integral in upper 100m
          id_z_o2min       = -1,       & ! Depth of oxygen minimum
          id_z_sat_arag    = -1,       & ! Depth of Aragonite saturation
          id_z_sat_calc    = -1          ! Depth of Calcite saturation

  end type generic_TOPAZ_type

  !An auxiliary type for storing variable names
  type, public :: vardesc
     character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
     character(len=fm_string_len) :: longname ! The long name of that variable.
     character(len=1)  :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
     character(len=1)  :: z_grid   ! The vert. grid:  L, i, or 1.
     character(len=1)  :: t_grid   ! The time description: s, a, m, or 1.
     character(len=fm_string_len) :: units    ! The dimensions of the variable.
     character(len=1)  :: mem_size ! The size in memory: d or f.
  end type vardesc

  type(generic_TOPAZ_type), save  :: topaz
  
  type(CO2_dope_vector) :: CO2_dope_vec

contains

  subroutine generic_TOPAZ_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_TOPAZ_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

    
  end subroutine generic_TOPAZ_register

  ! <SUBROUTINE NAME="generic_TOPAZ_init">
  !  <OVERVIEW>
  !   Initialize the generic TOPAZ module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the TOPAZ Tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_TOPAZ_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_TOPAZ_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_TOPAZ_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate all the private work arrays used by this module.
    call user_allocate_arrays

  end subroutine generic_TOPAZ_init

  !   Register diagnostic fields to be used in this module. 
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !
  subroutine generic_TOPAZ_register_diag(diag_list)
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

    !Below are the rate fields needed for IPCC5  as _z in GOLD:
    !set varlist  = (jprod_ntot jprod_no3tot jprod_ptot jprod_fetot jprod_sio4_Lg   \
    !                jprod_cadet_calc  jprod_cadet_arag jcadet_calc jcadet_arag     \
    !                jprod_nlg_diatoms jprod_nlg_nondiatoms jprod_ndi jprod_nsm     \
    !                jdic jdin jpo4 jfed jsio4 jalk                                 \
    !                jfe_ads jfe_des jfedet jgraz_ntot                              )
    !
    !In order to do this we register the diag fields via a subroutine call g_diag_field_add
    !that is basically a wrapper around the fms register_diag_field unless an optional argument
    !is passed to it which in that case delegates the registeration and subsequent send data
    !to GOLD.

    !Missing sends:jprod_sio4_Lg , jprod_ndi, jprod_nsm

    !JPD/JGJ: also need 16 tracer fields below for IPCC5 as _z in GOLD:
    !  tot_doc tot_phyton nlg_diatoms nlg_nondiatoms chl_lg_diatoms chl_lg_nondiatoms
    !  chl_di chl_sm tot_pon tot_pop tot_bfe tot_bsi tot_phytop tot_phytofe
    !  co3_sol_arag co3_sol_calc
    !2010/05/28 add 9 more fields for rates with btm fluxes added
    !  jalk_plus_btm, jdic_plus_btm, jdin_plus_btm, jfed_plus_btm
    !  jnh4_plus_btm, jno4_plus_btm, jo2_plus_btm, jpo4_plus_btm, jsio4_plus_btm
    !2012/01/09 JGJ add 4 more rate fields for JPD: 
    ! jgraz_n_lg_z, jnitrif_z, jprod_ndet_z, jndet_z
    !
    !2012/01/31 JGJ modify output diagnostics to get irr_inst_z

    ! Register rho_dzt

    vardesc_temp = vardesc("rho_dzt","ocean t-cell rho*thickness",'h','L','s','(kg/m^3)*m','f')
    id_rho_dzt = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Register the diagnostics for the various phytoplankton 
    !
    ! Register Limitation Diagnostics
    !
    vardesc_temp = vardesc("chl_Di","Chlorophyll in Diaz. Phyto",'h','L','s','ug kg-1','f')
    call g_diag_field_add(diag_list, phyto(DIAZO)%id_chl, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="chl_Di" ,Zlongname="Chlorophyll in Diaz. Phyto", &
         Zunits="ug kg-1", field_ptr=phyto(DIAZO)%chl)

    vardesc_temp = vardesc("chl_Lg","Chlorophyll in Large Phyto",'h','L','s','ug kg-1','f')
    phyto(LARGE)%id_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("chl_Sm","Chlorophyll in Small Phyto",'h','L','s','ug kg-1','f')
    call g_diag_field_add(diag_list, phyto(SMALL)%id_chl, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="chl_Sm" ,Zlongname="Chlorophyll in Small Phyto", &
         Zunits="ug kg-1", field_ptr=phyto(SMALL)%chl)

    vardesc_temp = vardesc("def_fe_Di","Diaz. Phyto. Iron Deficiency",'h','L','s','unknown units','f')
    phyto(DIAZO)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("def_fe_Lg","Large Phyto. Iron Deficiency",'h','L','s','unknown units','f')
    phyto(LARGE)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("def_fe_Sm","Small Phyto. Iron Deficiency",'h','L','s','unknown units','f')
    phyto(SMALL)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Di","Diaz. Phyto. Fed Limitation",'h','L','s','unknown units','f')
    phyto(DIAZO)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Lg","Large Phyto. Fed Limitation",'h','L','s','unknown units','f')
    phyto(LARGE)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Sm","Small Phyto. Fed Limitation",'h','L','s','unknown units','f')
    phyto(SMALL)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Di","Diaz. Phyto. Light Limitation",'h','L','s','unknown units','f')
    phyto(DIAZO)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Lg","Large Phyto. Light Limitation",'h','L','s','unknown units','f')
    phyto(LARGE)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Sm","Small Phyto. Light Limitation",'h','L','s','unknown units','f')
    phyto(SMALL)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("liebig_lim_Di","Diaz. Phyto. Overall Nutrient Limitation",'h','L','s','unknown units','f')
    phyto(DIAZO)%id_liebig_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("liebig_lim_Lg","Large Phyto. Overall Nutrient Limitation",'h','L','s','unknown units','f')
    phyto(LARGE)%id_liebig_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("liebig_lim_Sm","Small Phyto. Overall Nutrient Limitation",'h','L','s','unknown units','f')
    phyto(SMALL)%id_liebig_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Di","Diaz. Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(DIAZO)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Lg","Large Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(LARGE)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Sm","Small Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(SMALL)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Lg","Ammonia Limitation of Large Phyto",'h','L','s','unknown units','f')
    phyto(LARGE)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Sm","Ammonia Limitation of Small Phyto",'h','L','s','unknown units','f')
    phyto(SMALL)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Lg","Nitrate Limitation of Large Phyto",'h','L','s','unknown units','f')
    phyto(LARGE)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Sm","Nitrate Limitation of Small Phyto",'h','L','s','unknown units','f')
    phyto(SMALL)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Di","Phosphate Limitation of Diaz. Phyto",'h','L','s','unknown units','f')
    phyto(DIAZO)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Lg","Phosphate Limitation of Large Phyto",'h','L','s','unknown units','f')
    phyto(LARGE)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Sm","Phosphate Limitation of Small Phyto",'h','L','s','unknown units','f')
    phyto(SMALL)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Di","Fe:N ratio of Diaz. Phyto",'h','L','s','mol/mol','f')
    phyto(DIAZO)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Lg","Fe:N ratio of Large Phyto",'h','L','s','mol/mol','f')
    phyto(LARGE)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Sm","Fe:N ratio of Small Phyto",'h','L','s','mol/mol','f')
    phyto(SMALL)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_p_2_n_Di","P:N ratio of Diaz. Phyto",'h','L','s','mol/mol','f')
    phyto(DIAZO)%id_q_p_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_p_2_n_Lg","P:N ratio of Large Phyto",'h','L','s','mol/mol','f')
    phyto(LARGE)%id_q_p_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_p_2_n_Sm","P:N ratio of Small Phyto",'h','L','s','mol/mol','f')
    phyto(SMALL)%id_q_p_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_p_2_n_opt_Di","Optimal P:N ratio of Diaz. Phyto",'h','L','s','mol/mol','f')
    phyto(DIAZO)%id_q_p_2_n_opt = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_p_2_n_opt_Lg","Optimal P:N ratio of Large Phyto",'h','L','s','mol/mol','f')
    phyto(LARGE)%id_q_p_2_n_opt = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_p_2_n_opt_Sm","Optimal P:N ratio of Small Phyto",'h','L','s','mol/mol','f')
    phyto(SMALL)%id_q_p_2_n_opt = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_Di","Surface Chlorophyll in Diaz. Phyto",'h','L','s','ug kg-1','f')
    phyto(DIAZO)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_Lg","Surface Chlorophyll in Large Phyto",'h','L','s','ug kg-1','f')
    phyto(LARGE)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_Sm","Surface Chlorophyll in Small Phyto",'h','L','s','ug kg-1','f')
    phyto(SMALL)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("silim_Lg","SiO4 Limitation of Large Phyto",'h','L','s','unknown units','f')
    phyto(LARGE)%id_silim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register Grazing Diagnostics
    !

    vardesc_temp = vardesc("jgraz_n_Di_100","Diazotroph nitrogen grazing integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jgraz_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_n_Lg_100","Large nitrogen grazing integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jgraz_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_n_Sm_100","Small nitrogen grazing integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jgraz_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_fe_Di","Diazotroph iron grazing layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jgraz_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_fe_Lg","Large iron grazing layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jgraz_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_fe_Sm","Small iron grazing layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jgraz_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_n_Di","Diazotroph nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jgraz_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_n_Lg","Large nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, phyto(LARGE)%id_jgraz_n, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jgraz_n_Lg" ,Zlongname="Large nitrogen grazing", &
         Zunits="mol kg-1 s-1", field_ptr=phyto(LARGE)%jgraz_n)

    vardesc_temp = vardesc("jgraz_n_Sm","Small nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jgraz_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jgraz_sio2_Lg","Large silicon grazing layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jgraz_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    !
    ! Register Production Diagnostics
    !
    vardesc_temp = vardesc("jprod_n2_Di","Nitrogen fixation layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fe_Di","Diaz. phyto. Fed production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fe_Lg","Large phyto. Fed production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fe_Sm","Small phyto. Fed production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nDi","Diazotroph phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, phyto(DIAZO)%id_jprod_n, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_nDi" ,Zlongname="Diazotroph phyto. Nitrogen production", &
         Zunits="mol kg-1 s-1", field_ptr=phyto(DIAZO)%jprod_n)

    vardesc_temp = vardesc("jprod_nDi_100","Diazotroph phyto. Nitrogen production integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nLg","Large phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nLg_100","Large phyto. Nitrogen production integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nSm","Small phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, phyto(SMALL)%id_jprod_n, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_nSm" ,Zlongname="Small phyto. Nitrogen production", &
         Zunits="mol kg-1 s-1", field_ptr=phyto(SMALL)%jprod_n)

    vardesc_temp = vardesc("jprod_nSm_100","Small phyto. Nitrogen production integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Di","Diaz. phyto. NH4 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Lg","Large phyto. NH4 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Sm","Small phyto. NH4 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_no3_Di","Diaz. phyto. NO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_no3_Lg","Large phyto. NO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_no3_Sm","Small phyto. NO3 Production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Di","Diaz. phyto. PO4 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Lg","Large phyto. PO4 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Sm","Small phyto. PO4 production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Lg","Large phyto. SiO4 production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, phyto(LARGE)%id_jprod_sio4, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_sio4_Lg" ,Zlongname="Large phyto. SiO4 production", &
         Zunits="mol kg-1 s-1", field_ptr=phyto(LARGE)%jprod_sio4)

    vardesc_temp = vardesc("jprod_sio4_Lg_100","Large phyto. SiO4 production integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_sio4_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register non Phytoplankton Diagnostics
    !

    vardesc_temp = vardesc("f_alk_int_100","Inventory of alkalinity in upper 100m",'h','1','s','mol m-2','f')
    topaz%id_f_alk_int_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("f_dic_int_100","Inventory of Dissolved Inorganic Carbon in upper 100m",'h','1','s','mol m-2','f')
    topaz%id_f_dic_int_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("f_din_int_100","Inventory of Dissolved Inorganic Nitrogen in upper 100m",'h','1','s','mol m-2','f')
    topaz%id_f_din_int_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("f_fed_int_100","Inventory of Dissolved Iron in upper 100m",'h','1','s','mol m-2','f')
    topaz%id_f_fed_int_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("f_po4_int_100","Inventory of Phosphate in upper 100m",'h','1','s','mol m-2','f')
    topaz%id_f_po4_int_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("f_sio4_int_100","Inventory of Silicate in upper 100m",'h','1','s','mol m-2','f')
    topaz%id_f_sio4_int_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("cased_diag","Sediment CaCO3",'h','L','s','mol m-3','f')
    call g_diag_field_add(diag_list, topaz%id_cased_diag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="cased_diag" ,Zlongname="Sediment CaCO3", &
         Zunits="mol m-3", field_ptr= topaz%f_cased)

    vardesc_temp = vardesc("chl_diag","Chlorophyll",'h','L','s','ug kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_chl_diag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="chl_diag" ,Zlongname="Chlorophyll", &
         Zunits="ug kg-1", field_ptr= topaz%f_chl)

    vardesc_temp = vardesc("chl_lg_diatoms","large Diatom Chlorophyll",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_chl_lg_diatoms, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="chl_lg_diatoms" ,Zlongname="large Diatom Chlorophyll", &
         Zunits="mol kg-1", field_ptr=topaz%chl_lg_diatoms)

    vardesc_temp = vardesc("chl_lg_nondiatoms","large non-Diatom Chlorophyll",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_chl_lg_nondiatoms, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="chl_lg_nondiatoms" ,Zlongname="large non-Diatom Chlorophyll", &
         Zunits="mol kg-1", field_ptr=topaz%chl_lg_nondiatoms)

    vardesc_temp = vardesc("co3_ion_diag","Carbonate Ion",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_co3_ion_diag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="co3_ion_diag" ,Zlongname="Carbonate Ion", &
         Zunits="mol kg-1", field_ptr= topaz%f_co3_ion)

    vardesc_temp = vardesc("co3_sol_arag","Carbonate Ion Solubility for Aragonite",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_co3_sol_arag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="co3_sol_arag" ,Zlongname="Carbonate Ion Solubility for Aragonite", &
         Zunits="mol kg-1", field_ptr= topaz%co3_sol_arag)

    vardesc_temp = vardesc("co3_sol_calc","Carbonate Ion Solubility for Calcite",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_co3_sol_calc, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="co3_sol_calc" ,Zlongname="Carbonate Ion Solubility for Calcite", &
         Zunits="mol kg-1", field_ptr= topaz%co3_sol_calc)

    vardesc_temp = vardesc("fca_sed_sfc","Net CaCO3 flux across sediment surface (bottom sinking minus redissolution",'h','1','s','mol m-2 s-1','f')
    topaz%id_fca_sed_sfc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_int","Sediment CaCO3 integral",'h','1','s','mol m-3','f')
    topaz%id_fcased_int = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_arag","Aragonite CaCO3 detrital sinking flux",'h','L','s','mol m-2 s-1','f')
    topaz%id_fcadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_arag_100","Aragonite CaCO3 sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_fcadet_arag_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_arag_btm","Aragonite CaCO3 sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    topaz%id_fcadet_arag_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc","Calcite CaCO3 detrital sinking flux",'h','L','s','mol m-2 s-1','f')
    topaz%id_fcadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_100","Calcite CaCO3 sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_fcadet_calc_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_btm","Calcite CaCO3 sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    topaz%id_fcadet_calc_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_burial","CaCO3 permanent burial flux",'h','1','s','mol m-2 s-1','f')
    topaz%id_fcased_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_redis","CaCO3 redissolution from sediments",'h','1','s','mol m-2 s-1','f')
    topaz%id_fcased_redis = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet","Iron detrital sinking flux",'h','L','s','mol m-2 s-1','f')
    topaz%id_ffedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_100","Iron detrital sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_ffedet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_btm","Iron detrital sinking flux burial",'h','1','s','mol m-2 s-1','f')
    topaz%id_ffedet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffe_sed","Sediment iron efflux",'h','1','s','mol m-2 s-1','f')
    topaz%id_ffe_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet","Lithogenic detrital sinking flux",'h','L','s','g m-2 s-1','f')
    topaz%id_flithdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_100","Lithogenic detrital sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_flithdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_btm","Lithogenic detrital sinking flux burial",'h','1','s','g m-2 s-1','f')
    topaz%id_flithdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet","PON sinking flux",'h','L','s','mol m-2 s-1','f')
    topaz%id_fndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_100","PON sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_fndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_btm","PON sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    topaz%id_fndet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fnfeso4red_sed","Sediment Ndet, Iron and SO4 reduction flux",'h','1','s','mol m-2 s-1','f')
    topaz%id_fnfeso4red_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fno3denit_sed","Sediment denitrification flux",'h','1','s','mol m-2 s-1','f')
    topaz%id_fno3denit_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fno3denit_tot","Total denitrification flux",'h','1','s','mol m-2 s-1','f')
    topaz%id_fno3denit_tot = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fnoxic_sed","Sediment oxic Ndet remineralization flux",'h','1','s','mol m-2 s-1','f')
    topaz%id_fnoxic_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet","POP sinking flux",'h','L','s','mol m-2 s-1','f')
    topaz%id_fpdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_100","POP sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_fpdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_btm","POP sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    topaz%id_fpdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet","Si sinking flux",'h','L','s','mol m-2 s-1','f')
    topaz%id_fsidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_100","Si sinking flux at 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_fsidet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_btm","Si sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    topaz%id_fsidet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("htotal_diag","H+ ion concentration",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_htotal_diag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="htotal_diag" ,Zlongname="H+ ion concentration", &
         Zunits="mol kg-1", field_ptr= topaz%f_htotal)

    vardesc_temp = vardesc("irr_inst","Instantaneous Light",'h','L','s','W m-2','f')
    call g_diag_field_add(diag_list, topaz%id_irr_inst, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="irr_inst" ,Zlongname="Instantaneous Light", &
         Zunits="W m-2", field_ptr= topaz%irr_inst)

    vardesc_temp = vardesc("irr_mem_diag","Irradiance memory",'h','L','s','W m-2','f')
    call g_diag_field_add(diag_list, topaz%id_irr_mem_diag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="irr_mem_diag" ,Zlongname="Irradiance memory", &
         Zunits="W m-2", field_ptr= topaz%f_irr_mem)

    vardesc_temp = vardesc("irr_mix","Light averaged over mixing layer",'h','L','s','W m-2','f')
    topaz%id_irr_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jalk_100","Alkalinity source layer integral in upper 100m",'h','L','s','eq m-2 s-1','f')
    topaz%id_jalk_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jalk","Alkalinity source layer integral",'h','L','s','eq m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jalk, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jalk" ,Zlongname="Alkalinity source", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jalk)

    vardesc_temp = vardesc("jalk_plus_btm","Alkalinity source plus btm layer integral",'h','L','s','eq m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jalk_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jalk_plus_btm" ,Zlongname="Alkalinity source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jalk_plus_btm)

    vardesc_temp = vardesc("jcadet_arag","Aragonite CaCO3 change layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jcadet_arag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jcadet_arag" ,Zlongname="Aragonite CaCO3 change", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jcadet_arag)

    vardesc_temp = vardesc("jcadet_calc","Calcite CaCO3 change layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jcadet_calc, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jcadet_calc" ,Zlongname="Calcite CaCO3 change", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jcadet_calc)

    vardesc_temp = vardesc("jdic","Dissolved Inorganic Carbon source layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jdic, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jdic" ,Zlongname="Dissolved Inorganic Carbon source", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jdic)

    vardesc_temp = vardesc("jdic_plus_btm","Dissolved Inorganic Carbon source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jdic_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jdic_plus_btm" ,Zlongname="Dissolved Inorganic Carbon source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jdic_plus_btm)

    vardesc_temp = vardesc("jdic_100","Dissolved Inorganic Carbon source integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jdic_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdin","Dissolved Inorganic Nitrogen (NO3+NH4) source layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jdin, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jdin" ,Zlongname="Dissolved Inorganic Nitrogen (NO3+NH4) source", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jdin)

    vardesc_temp = vardesc("jdin_plus_btm","Dissolved Inorganic Nitrogen (NO3+NH4) source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jdin_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jdin_plus_btm" ,Zlongname="Dissolved Inorganic Nitrogen (NO3+NH4) source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jdin_plus_btm)

    vardesc_temp = vardesc("jdin_100","Dissolved Inorganic Nitrogen (NO3+NH4) source integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jdin_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdiss_sio2","SiO2 Dissolution during grazing layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jdiss_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfe_ads","Iron adsorption layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jfe_ads, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jfe_ads" ,Zlongname="Iron adsorption", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jfe_ads)

    vardesc_temp = vardesc("jfe_des","Iron desorption layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jfe_des, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jfe_des" ,Zlongname="Iron desorption", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jfe_des)

    vardesc_temp = vardesc("jfe_graz","Dissolved iron grazing source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jfe_graz = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfe_coast","Coastal iron efflux layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jfe_coast = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfed","Dissolved iron source layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jfed, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jfed" ,Zlongname="Dissolved iron source", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jfed)

    vardesc_temp = vardesc("jfed_plus_btm","Dissolved iron source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jfed_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jfed_plus_btm" ,Zlongname="Dissolved iron source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jfed_plus_btm)

    vardesc_temp = vardesc("jfed_100","Dissolved iron source integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jfed_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfedet","Loss of sinking iron layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jfedet, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jfedet" ,Zlongname="Loss of sinking iron", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jfedet)

    vardesc_temp = vardesc("jgraz_ntot","Total Nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jgraz_ntot, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jgraz_ntot" ,Zlongname="Total Nitrogen grazing", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jgraz_ntot)

    vardesc_temp = vardesc("jgraz_ntot_100","Total Nitrogen grazing integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jgraz_ntot_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jldon","Labile DON source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jndet","Loss of sinking nitrogen layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jndet, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jndet" ,Zlongname="Loss of sinking nitrogen", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jndet)

    vardesc_temp = vardesc("jnh4","NH4 source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jnh4_plus_btm","NH4 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jnh4_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jnh4_plus_btm" ,Zlongname="NH4 source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jnh4_plus_btm)

    vardesc_temp = vardesc("jnh4_graz","NH4 grazing source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jnh4_graz = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jnhet","heterotrophic Nitrogen remineralization layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jnhet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jnitrif","Nitrification layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jnitrif, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jnitrif" ,Zlongname="Nitrification", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jnitrif)

    vardesc_temp = vardesc("jnitrif_100","Nitrification integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jnitrif_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jno3","NO3 source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jno3_plus_btm","NO3 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jno3_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jno3_plus_btm" ,Zlongname="NO3 source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jno3_plus_btm)

    vardesc_temp = vardesc("jno3denit_wc","Water column Denitrification layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jno3denit_wc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jo2","O2 source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jo2_plus_btm","O2 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jo2_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jo2_plus_btm" ,Zlongname="O2 source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jo2_plus_btm)

    vardesc_temp = vardesc("jpdet","Loss of sinking phosphorus layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jpdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jpo4","PO4 source layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jpo4, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jpo4" ,Zlongname="PO4 source", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jpo4)

    vardesc_temp = vardesc("jpo4_plus_btm","PO4 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jpo4_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jpo4_plus_btm" ,Zlongname="PO4 source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jpo4_plus_btm)

    vardesc_temp = vardesc("jpo4_100","PO4 source integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jpo4_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jpo4_graz","PO4 source from grazing layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jpo4_graz = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_arag","Aragonite CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_cadet_arag, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_cadet_arag" ,Zlongname="Aragonite CaCO3 production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_cadet_arag)

    vardesc_temp = vardesc("jprod_cadet_arag_100","Aragonite CaCO3 production integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jprod_cadet_arag_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_calc","Calcite CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_cadet_calc, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_cadet_calc" ,Zlongname="Calcite CaCO3 production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_cadet_calc)

    vardesc_temp = vardesc("jprod_cadet_calc_100","Calcite CaCO3 production layer integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jprod_cadet_calc_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fetot","Total phytoplankton iron production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_fetot, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_fetot" ,Zlongname="Total phytoplankton iron production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_fetot)

    vardesc_temp = vardesc("jprod_fetot_100","Total Iron production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jprod_fetot_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet","Detrital Fedet production layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_lithdet","Lithogenic removal to sinking layer integral",'h','L','s','g m-2 s-1','f')
    topaz%id_jprod_lithdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet","Detrital PON production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_ndet, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_ndet" ,Zlongname="Detrital PON production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_ndet)

    vardesc_temp = vardesc("jprod_nhet","heterotrophic Nitrogen Production layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jprod_nhet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlg_diatoms","Large Diatom Nitrogen production layer integral",'h','1','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_nlg_diatoms, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_nlg_diatoms" ,Zlongname="Large Diatom Nitrogen production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_nlg_diatoms)

    vardesc_temp = vardesc("jprod_nlg_diatoms_100","Large Diatom Nitrogen production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jprod_nlg_diatoms_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlg_nondiatoms","Large non-Diatom Nitrogen production layer integral",'h','1','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_nlg_nondiatoms, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_nlg_nondiatoms" ,Zlongname="Large non-Diatom Nitrogen production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_nlg_nondiatoms)

    vardesc_temp = vardesc("jprod_nlg_nondiatoms_100","Large non-Diatom Nitrogen production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jprod_nlg_nondiatoms_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_no3tot","Total phytoplankton nitrate production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_no3tot, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_no3tot" ,Zlongname="Total phytoplankton nitrate production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_no3tot)

    vardesc_temp = vardesc("jprod_no3tot_100","Total NO3 production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jprod_no3tot_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ntot","Total phytoplankton nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_ntot, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_ntot" ,Zlongname="Total phytoplankton nitrogen production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_ntot)

    vardesc_temp = vardesc("jprod_ntot_100","Total Nitrogen production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jprod_ntot_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet","Detrital phosphorus production layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ptot","Total phytoplankton phosphorus production layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jprod_ptot, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jprod_ptot" ,Zlongname="Total phytoplankton phosphorus production", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jprod_ptot)

    vardesc_temp = vardesc("jprod_ptot_100","Total Phosphorus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    topaz%id_jprod_ptot_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jsdop","Semilabile DOP source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jsdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jsidet","Loss of sinking silicon layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jsidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jsio4","SiO4 source layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jsio4, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jsio4" ,Zlongname="SiO4 source", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jsio4)

    vardesc_temp = vardesc("jsio4_plus_btm","SiO4 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    call g_diag_field_add(diag_list, topaz%id_jsio4_plus_btm, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="jsio4_plus_btm" ,Zlongname="SiO4 source plus btm", &
         Zunits="mol kg-1 s-1", field_ptr=topaz%jsio4_plus_btm)

    vardesc_temp = vardesc("jsio4_100","SiO4 source integral in upper 100m",'h','L','s','mol m-2 s-1','f')
    topaz%id_jsio4_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jsdon","Semilabile DON source layer integral",'h','L','s','mol m-2 s-1','f')
    topaz%id_jsdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nLg_diatoms","Large diatom nitrogen",'h','L','s','unknown units','f')
    call g_diag_field_add(diag_list, topaz%id_nLg_diatoms, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="nLg_diatoms" ,Zlongname="Large diatom nitrogen", &
         Zunits="unknown units", field_ptr=topaz%nLg_diatoms)

    vardesc_temp = vardesc("nLg_nondiatoms","Large non-diatom nitrogen",'h','L','s','unknown units','f')
    call g_diag_field_add(diag_list, topaz%id_nLg_nondiatoms, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="nLg_nondiatoms" ,Zlongname="Large non-diatom nitrogen", &
         Zunits="unknown units", field_ptr=topaz%nLg_nondiatoms)

    vardesc_temp = vardesc("no3_in_source","Nitrate in source routine",'h','L','s','mol/kg','f')
    topaz%id_no3_in_source = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("omega_arag","Carbonate Ion Saturation State for Aragonite",'h','L','s','mol kg-1','f')
    topaz%id_omega_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("omega_calc","Carbonate Ion Saturation State for Calcite",'h','L','s','mol kg-1','f')
    topaz%id_omega_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_bfe","Total Particulate Biogenic Iron",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_bfe, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_bfe" ,Zlongname="Total Particulate Biogenic Iron", &
         Zunits="mol kg-1", field_ptr=topaz%tot_bfe)

    vardesc_temp = vardesc("tot_bsi","Total Particulate Biogenic Silica",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_bsi, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_bsi" ,Zlongname="Total Particulate Biogenic Silica", &
         Zunits="mol kg-1", field_ptr=topaz%tot_bsi)

    vardesc_temp = vardesc("tot_doc","Total Dissolved Organic Carbon (refractory+sdoc+ldoc)",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_doc, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_doc" ,Zlongname="Total Dissolved Organic Carbon (refractory+sdoc+ldoc)", &
         Zunits="mol kg-1", field_ptr=topaz%tot_doc)

    vardesc_temp = vardesc("tot_layer_int_c","Total Carbon (DIC+OC+IC) boxwise layer integral",'h','L','s','mol m-2','f')
    topaz%id_tot_layer_int_c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_fe","Total Iron (Fed_OFe) boxwise layer integral",'h','L','s','mol m-2','f')
    topaz%id_tot_layer_int_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_n","Total Nitrogen (NO3+NH4+ON) boxwise layer integral",'h','L','s','mol m-2','f')
    topaz%id_tot_layer_int_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_p","Total Phosphorus (PO4+OP) boxwise layer integral",'h','L','s','mol m-2','f')
    topaz%id_tot_layer_int_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_si","Total Silicon (SiO4+SiO2) boxwise layer integral",'h','L','s','mol m-2','f')
    topaz%id_tot_layer_int_si = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_phytofe","Total Phytoplankton Iron (Di+Lg+Sm)",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_phytofe, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_phytofe" ,Zlongname="Total Phytoplankton Iron (Di+Lg+Sm)", &
         Zunits="mol kg-1", field_ptr=topaz%tot_phytofe)

    vardesc_temp = vardesc("tot_phyton","Total Phytoplankton Nitrogen (Di+Lg+Sm)",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_phyton, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_phyton" ,Zlongname="Total Phytoplankton Nitrogen (Di+Lg+Sm)", &
         Zunits="mol kg-1", field_ptr=topaz%tot_phyton)

    vardesc_temp = vardesc("tot_phytop","Total Phytoplankton Phosphorus (Di+Lg+Sm)",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_phytop, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_phytop" ,Zlongname="Total Phytoplankton Phosphorus (Di+Lg+Sm)", &
         Zunits="mol kg-1", field_ptr=topaz%tot_phytop)

    vardesc_temp = vardesc("tot_pon","Total Particulate Organic Nitrogen",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_pon, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_pon" ,Zlongname="Total Particulate Organic Nitrogen", &
         Zunits="mol kg-1", field_ptr=topaz%tot_pon)

    vardesc_temp = vardesc("tot_pop","Total Particulate Organic Phosphorus",'h','L','s','mol kg-1','f')
    call g_diag_field_add(diag_list, topaz%id_tot_pop, package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         Z_diag=1, Zname="tot_pop" ,Zlongname="Total Particulate Organic Phosphorus", &
         Zunits="mol kg-1", field_ptr=topaz%tot_pop)

    !
    ! Surface Diagnostics
    !

    vardesc_temp = vardesc("sfc_alk","Surface Alkalinity",'h','1','s','eq kg-1','f')
    topaz%id_sfc_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_cadet_arag","Surface Detrital Aragonite",'h','1','s','mol kg-1','f')
    topaz%id_sfc_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_cadet_calc","Surface Detrital Calcite",'h','1','s','mol kg-1','f')
    topaz%id_sfc_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl","Surface Chl",'h','1','s','ug kg-1','f')
    topaz%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_lg_diatoms","Surface Large Diatom Chl",'h','1','s','ug kg-1','f')
    topaz%id_sfc_chl_lg_diatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_lg_nondiatoms","Surface Large nonDiatom Chl",'h','1','s','ug kg-1','f')
    topaz%id_sfc_chl_lg_nondiatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_ion","Surface Carbonate Ion",'h','1','s','mol kg-1','f')
    topaz%id_sfc_co3_ion = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_sol_arag","Surface Carbonate Ion Solubility for Aragonite",'h','1','s','mol kg-1','f')
    topaz%id_sfc_co3_sol_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_sol_calc","Surface Carbonate Ion Solubility for Calcite ",'h','1','s','mol kg-1','f')
    topaz%id_sfc_co3_sol_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_dic","Surface Dissolved Inorganic Carbon",'h','1','s','mol kg-1','f')
    topaz%id_sfc_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_fed","Surface Dissolved Iron",'h','1','s','mol kg-1','f')
    topaz%id_sfc_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_htotal","Surface Htotal",'h','1','s','mol kg-1','f')
    topaz%id_sfc_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irr","Surface Irradiance",'h','1','s','W m-2','f')
    topaz%id_sfc_irr = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irr_mem","Surface Irradiance memory",'h','1','s','W m-2','f')
    topaz%id_sfc_irr_mem = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ndet","Surface Detrital N",'h','1','s','mol kg-1','f')
    topaz%id_sfc_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ndi","Surface Diazotroph N",'h','1','s','mol kg-1','f')
    topaz%id_sfc_ndi = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4","Surface NH4",'h','1','s','mol kg-1','f')
    topaz%id_sfc_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nhet","Surface Heterotroph N",'h','1','s','mol kg-1','f')
    topaz%id_sfc_nhet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nlg","Surface Large Phytoplankton N",'h','1','s','mol kg-1','f')
    topaz%id_sfc_nlg = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nlg_diatoms","Surface Large Diatoms",'h','1','s','mol kg-1','f')
    topaz%id_sfc_nlg_diatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nlg_nondiatoms","Surface Large non-Diatoms",'h','1','s','mol kg-1','f')
    topaz%id_sfc_nlg_nondiatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3","Surface NO3",'h','1','s','mol kg-1','f')
    topaz%id_sfc_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nsm","Surface Small Phytoplankton Nitrogen",'h','1','s','mol kg-1','f')
    topaz%id_sfc_nsm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_o2","Surface Oxygen",'h','1','s','mol kg-1','f')
    topaz%id_sfc_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4","Surface PO4",'h','1','s','mol kg-1','f')
    topaz%id_sfc_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_silg","Surface Large Phytoplankton Silicon",'h','1','s','mol kg-1','f')
    topaz%id_sfc_silg = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_sio4","Surface SiO4",'h','1','s','mol kg-1','f')
    topaz%id_sfc_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_temp","Surface Temperature",'h','1','s','deg C','f')
    topaz%id_sfc_temp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_bfe","Surface Total Particulate Biogenic Iron",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_bfe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_bsi","Surface Total Particulate Biogenic Silicon",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_bsi = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_doc","Surface Total Dissolved Organic Carbon",'h','1','s','mol kg-1','f')
    topaz%id_sfc_tot_doc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_phytofe","Surface Total Phytoplankton Iron (Di+Lg+Sm)",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_phytofe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_phyton","Surface Total Phytoplankton Nitrogen (Di+Lg+Sm)",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_phyton = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_phytop","Surface Total Phytoplankton Phosphorus (Di+Lg+Sm)",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_phytop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_pon","Surface Total Particulate Organic Nitrogen",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_pon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_tot_pop","Surface Total Particulate Organic Phosphorus",'h','L','s','mol kg-1','f')
    topaz%id_sfc_tot_pop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_fed","Dry Deposition of Iron to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_dep_dry_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_nh4","Dry Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_dep_dry_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_no3","Dry Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_dep_dry_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_fed","Wet Deposition of Iron to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_dep_wet_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_nh4","Wet Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_dep_wet_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_no3","Wet Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_dep_wet_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_alk","Alkalinity runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_runoff_flux_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_dic","Dissolved Inorganic Carbon runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_runoff_flux_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_fed","Iron runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_runoff_flux_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_lith","Lithogenic runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_runoff_flux_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_nh4","Ammonia runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_runoff_flux_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_no3","Nitrate runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_runoff_flux_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nongas_source_c","Sum of all nongas sources of carbon to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_nongas_source_c = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nongas_source_fe","Sum of all nongas sources of iron to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_nongas_source_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nongas_source_n","Sum of all nongas sources of nitrogen to the ocean",'h','1','s','mol m-2 s-1','f')
    topaz%id_nongas_source_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2min","Minimum Oxygen",'h','1','s','mol kg-1','f')
    topaz%id_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("pco2surf","Oceanic pCO2",'h','1','s','uatm','f')
    topaz%id_pco2surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_c","Total Carbon (DIC+OC+IC) vertical integral",'h','1','s','mol m-2','f')
    topaz%id_wc_vert_int_c = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_dic","Dissolved Inorganic Carbon vertical integral",'h','1','s','mol m-2','f')
    topaz%id_wc_vert_int_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jfe_coast","Coastal iron efflux  vertical integral",'h','1','s','mol m-2 s-1','f')
    topaz%id_wc_vert_int_jfe_coast = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jno3denit","Water Column Denitrification vertical integral",'h','1','s','mol m-2 s-1','f')
    topaz%id_wc_vert_int_jno3denit = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_nfix","Nitrogen fixation vertical integral",'h','1','s','mol N m-2 s-1','f')
    topaz%id_wc_vert_int_nfix = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("z_o2min","Depth of Oxygen minimum",'h','1','s','m','f')
    topaz%id_z_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("z_sat_arag","Depth of Aragonite Saturation",'h','1','s','m','f')
    topaz%id_z_sat_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         mask_variant=.TRUE.)

    vardesc_temp = vardesc("z_sat_calc","Depth of Calcite Saturation",'h','1','s','m','f')
    topaz%id_z_sat_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         mask_variant=.TRUE.)

  end subroutine generic_TOPAZ_register_diag

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_TOPAZ_params type
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

    !     g_tracer_add_param(name   , variable   ,  default_value)
    !call g_tracer_add_param('', topaz%,  )

    call g_tracer_add_param('init', topaz%init, .false. )

    call g_tracer_add_param('htotal_scale_lo', topaz%htotal_scale_lo, 0.01)
    call g_tracer_add_param('htotal_scale_hi', topaz%htotal_scale_hi, 100.0)

    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    call g_tracer_add_param('RHO_0', topaz%Rho_0, 1035.0)
    call g_tracer_add_param('NKML' , topaz%nkml, 1)
    !-----------------------------------------------------------------------
    !       coefficients for O2 saturation
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_0', topaz%a_0, 2.00907)
    call g_tracer_add_param('a_1', topaz%a_1, 3.22014)
    call g_tracer_add_param('a_2', topaz%a_2, 4.05010)
    call g_tracer_add_param('a_3', topaz%a_3, 4.94457)
    call g_tracer_add_param('a_4', topaz%a_4, -2.56847e-01)
    call g_tracer_add_param('a_5', topaz%a_5, 3.88767)
    call g_tracer_add_param('b_0', topaz%b_0, -6.24523e-03)
    call g_tracer_add_param('b_1', topaz%b_1, -7.37614e-03)
    call g_tracer_add_param('b_2', topaz%b_2, -1.03410e-02 )
    call g_tracer_add_param('b_3', topaz%b_3, -8.17083e-03)
    call g_tracer_add_param('c_0', topaz%c_0, -4.88682e-07)
    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !-----------------------------------------------------------------------
    !
    !  Compute the Schmidt number of CO2 in seawater using the 
    !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    !  7373-7382).
    !-----------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_co2', topaz%a1_co2,  2068.9)
    call g_tracer_add_param('a2_co2', topaz%a2_co2, -118.63)
    call g_tracer_add_param('a3_co2', topaz%a3_co2,  2.9311)
    call g_tracer_add_param('a4_co2', topaz%a4_co2, -0.027)
    !---------------------------------------------------------------------
    !  Compute the Schmidt number of O2 in seawater using the 
    !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
    !  Cycles, 12, 141-163).
    !---------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_o2', topaz%a1_o2, 1929.7)
    call g_tracer_add_param('a2_o2', topaz%a2_o2, -117.46)
    call g_tracer_add_param('a3_o2', topaz%a3_o2, 3.116)
    call g_tracer_add_param('a4_o2', topaz%a4_o2, -0.0306)
    !
    !-----------------------------------------------------------------------
    ! Stoichiometry
    !-----------------------------------------------------------------------
    !
    ! Values taken from OCMIP-II Biotic protocols after Anderson (1995) for an
    ! organic material of C106H172O38N16(H3PO4) which gives an average oxidation state
    ! for carbon of (3*16+2*38-172)/106 = -0.4528.  These calculations ignore H3PO4.
    !
    ! Nitrate Production:
    !   16*H+ + 16*NO3- + 106*CO2 + 78*H2O <-> C106H172O38N16 + 150*O2
    !   Effect is to increase alkalinity by 16 NO3 equivalents.
    !
    ! Ammonia Production (and reverse for remineralization):
    !   16*NH4+ + 106*CO2 + 62*H2O <-> C106H172O38N16 + 118*O2 + 16*H+
    !   Effect is to decrease (increase for remineralization) alkalinity by 16 NH4 equivalents.
    !
    ! N2 Production:
    !   8*N2 + 106*CO2 + 86*H2O <-> C106H172O38N16 + 130*O2
    !   No effect on alkalinity.
    !
    ! Nitrification:
    !   NH4+ + 2*O2 <-> NO3- + H2O + 2*H+
    !   Effect is to decrease alkalinity by 2 NH4 equivalents.
    !
    ! Denitrification:
    !   C106H172O38N16 + 472/5*NO3- + 552/5*H+ <-> 106*CO2 + 16*NH4+ + 236/5*N2 + 546/5*H2O 
    !   Effect is to increase alkalinity by 552/472 = 1.169 NO3 equivalents.
    !   
    call g_tracer_add_param('alk_2_n_denit', topaz%alk_2_n_denit, 552.0 / 472.0)            ! eq. alk mol NO3-1
    call g_tracer_add_param('c_2_n', topaz%c_2_n, 106.0 / 16.0)                             ! mol C mol N-1
    call g_tracer_add_param('mass_2_n', topaz%mass_2_n, 106.0 / 16.0 * 12.0 * 1.87)         ! g mol N-1
    call g_tracer_add_param('n_2_n_denit', topaz%n_2_n_denit, 472.0 / (5.0 * 16.0))         ! mol N NO3 mol N org-1
    call g_tracer_add_param('o2_2_c', topaz%o2_2_c, 150.0 / 106.0)                          ! mol O2 mol C-1
    call g_tracer_add_param('o2_2_nfix', topaz%o2_2_nfix, (118.0+3.0/(5.0+3.0)*(150.0-118.0))/16.0) ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_nh4', topaz%o2_2_nh4, 118.0 / 16.0)                       ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_nitrif', topaz%o2_2_nitrif, 2.0)                          ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_no3', topaz%o2_2_no3, 150.0 / 16.0)                       ! mol O2 mol N-1
    !
    ! CaCO3 to nitrogen uptake ratio set in order to obtain a global CaCO3 export flux out of
    ! the euphotic zone of approximately 0.6 Pg C in CaCO3 after Sarmiento et al (2002)
    ! 2010/05/10 jgj: ca_2_n_calc changed to 0.005*106/16 (was 0.010*106/16)
    !
    call g_tracer_add_param('ca_2_n_arag', topaz%ca_2_n_arag, 0.010 * 106.0 / 16.0)           ! mol Ca mol N-1
    call g_tracer_add_param('ca_2_n_calc', topaz%ca_2_n_calc, 0.005 * 106.0 / 16.0)           ! mol Ca mol N-1
    !
    ! Upper limit with which CaCO3 supersaturation (co3/co3_sol - 1) can modulate the CaCO3
    ! to nitrogen uptake ratio
    !
    call g_tracer_add_param('caco3_sat_max', topaz%caco3_sat_max, 10.0)                     ! dimensionless
    !
    ! Static small phytoplankton/coccholith/foram-grazed/pteropod-grazed CaCO3 to
    ! nitrogen uptake ratio for a single, globally constant value instead of a
    ! dynamic value.
    !
    call g_tracer_add_param('ca_2_n_static', topaz%ca_2_n_static, .false.)
    call g_tracer_add_param('ca_2_n_het_static', topaz%ca_2_n_het_static, 0.012 * 106.0 / 16.0) ! mol Ca mol N-1
    !
    ! P:N limitation of phytoplankton growth, is taken from the fraction of the
    ! cell dedicated to assembly after Klausmeier et al (2004, Optimal nitrogen-to-phosphorus
    ! stoichiometry of phytoplankton, Nature, 429, 171-174), but also consistent with other
    ! work on the limitation of growth and optimal stoichiometry by Geider et al (1998, A 
    ! dynamic regulatory model of phytoplanktonic acclimation to light, nutrients,
    ! and temperature, Limnol. Oceanogr., 43, 679-694), Christian (2005, 
    ! Biogeochemical cycling in the oligotrophic ocean: Redfield and non-Redfield  
    ! models, Limnol. Oceanogr., 50, 646-657) and Hall et al. (2005,
    ! Constratints on primary producer N:P stoichiometry along N:P supply ratio gradients
    ! Ecology, 86, 1894-1904).  p_2_n_assem is based on the N:P of ribosomes of 7.2 for
    ! eukaryotes (large) and 6.0 for prokaryotes (small, diaz) from Elser et al. (1996, Organism
    ! size, life history, and N:P stoichiometry, BioScience, 46, 674-684).
    !
    call g_tracer_add_param('p_2_n_assem_Di', phyto(DIAZO)%p_2_n_assem, 1.0 / 6.0 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_assem_Lg', phyto(LARGE)%p_2_n_assem, 1.0 / 7.2 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_assem_Sm', phyto(SMALL)%p_2_n_assem, 1.0 / 6.0 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_photo', topaz%p_2_n_photo, 0.0128 )                      ! mol P mol N-1
    !
    ! Standard P:N value of 16 from Redfield, Ketchum and Richards (1963; The influence 
    ! of organisms on the composition of sea-water; The Sea, 2, 26-77).  Alternatively,
    ! some studies suggest a value of 15 from Goldman (1980) as 
    ! reprinted in Broecker and Peng (1982; Tracers in the Sea).
    !
    call g_tracer_add_param('p_2_n_RKR', topaz%p_2_n_RKR, 1.0 / 16.0)                       ! mol P mol N-1
    !
    ! In contrast to Klausmeier et al, a finite value of
    ! phosphorus in uptake machinery (p_2_n_uptake) is assumed in order to simulate
    ! the need for ATP in active transport.  This is set to p_2_n_photo somewhat arbitrarily as a
    ! low value.
    !
    call g_tracer_add_param('p_2_n_uptake', topaz%p_2_n_uptake, 0.0128 )                    ! mol P mol N-1
    !
    ! Calculation of the optimal P:N ratio is taken from Klausmeier et al (2004).  Allocation
    ! of phytoplankton cell structures to "other" is set to a minimum for large phytoplankton
    ! and small phytoplankton to be consistent with their different sizes, organelle structures and
    ! corresponding thetamax values.  Allocation to photosynthesis is taken from the realized
    ! Chlorophyll to Carbon ratio (theta) using the plastids to chlorophyll ratio.
    ! For diazotrophs, an additional allocation is specified for nitrogenase and other
    ! structures necessary for nitrogen fixation, and r_other_min is correspondingly lowered
    ! in order to match the observed range in N:P ratios of 14-182 in White, A. E., Y. H 
    ! Spitz, D. M. Karl and R. M. Letelier (2006, Flexible elemental stoichiometry in 
    ! Trichodesmium spp. and its ecological implications. Limnol Oceanogr, 51, 1777-1790.
    ! To avoid negative allocations, thetamax*plastid_2_chl+r_other_min+r_nfix must be < 1.
    !
    call g_tracer_add_param('r_nfix_Di',      phyto(DIAZO)%r_nfix, 0.5)                     ! dimensionless
    call g_tracer_add_param('r_nfix_Lg',      phyto(LARGE)%r_nfix, 0.0)                     ! dimensionless
    call g_tracer_add_param('r_nfix_Sm',      phyto(SMALL)%r_nfix, 0.0)                     ! dimensionless
    call g_tracer_add_param('r_other_min_Di', phyto(DIAZO)%r_other_min, 0.1)                ! dimensionless
    call g_tracer_add_param('r_other_min_Lg', phyto(LARGE)%r_other_min, 0.2)                ! dimensionless
    call g_tracer_add_param('r_other_min_Sm', phyto(SMALL)%r_other_min, 0.2)                ! dimensionless
    call g_tracer_add_param('r_uptake_max_Di',phyto(DIAZO)%r_uptake_max, 0.2)               ! dimensionless
    call g_tracer_add_param('r_uptake_max_Lg',phyto(LARGE)%r_uptake_max, 0.2)               ! dimensionless
    call g_tracer_add_param('r_uptake_max_Sm',phyto(SMALL)%r_uptake_max, 0.4)               ! dimensionless
    !
    ! The relationship between chlorophyll and carbon is taken from J. T. O. Kirk and R. A. E.
    ! Tilney-Bassett (1978, The Plastids: Their chemistry, structure, growth and inheritance.
    ! Revised second edition, Elsevier/North-Holland Biomedical Press, New York.)
    ! for eukaryotes (large) that have chloroplasts and 6.0 for prokaryotes (small, diaz) which
    ! only have the thylakoid.
    !
    call g_tracer_add_param('plast_2_chl_Di', phyto(DIAZO)%plast_2_chl, 1.0 / 1.87 / 0.08)  ! g C g Chl-1
    call g_tracer_add_param('plast_2_chl_Lg', phyto(LARGE)%plast_2_chl, 1.0 / 1.87 / 0.043) ! g C g Chl-1
    call g_tracer_add_param('plast_2_chl_Sm', phyto(SMALL)%plast_2_chl, 1.0 / 1.87 / 0.08)  ! g C g Chl-1
    !
    ! Static phosphorus to nitrogen uptake ratio for a single, globally constant value for each
    ! phytoplankton group instead of dynamic values for each.
    !
    call g_tracer_add_param('p_2_n_static', topaz%p_2_n_static, .false. )
    !
    ! If static, Diazotrophs are assumed to have low P:N after Letelier, R. M. & Karl, 
    ! D. M. Role of Trichodesmium spp. in the productivity of the subtropical North
    ! Pacific Ocean. Mar. Ecol. Prog. Ser. 133, 263-273 (1996).
    !
    call g_tracer_add_param('p_2_n_static_Di', phyto(DIAZO)%p_2_n_static, 1.0 / 40.0 )      ! mol P mol N-1
    !
    ! If static, Large and Small are assumed to have RKR P:N.
    !
    call g_tracer_add_param('p_2_n_static_Lg', phyto(LARGE)%p_2_n_static, 1.0 / 16.0 )      ! mol P mol N-1
    call g_tracer_add_param('p_2_n_static_Sm', phyto(SMALL)%p_2_n_static, 1.0 / 16.0 )      ! mol P mol N-1
    !
    ! Maximum diatom silicon to nitrogen uptake ratio after Mongin, M., D. M. Nelson,
    ! P. Pondaven, M. A. Brzezinski and P. Treguer (2003): Simulation of upper-ocean
    ! biogeochemistry with a flexible-composition phytoplankton model: C, N and Si
    ! cycling in the western Sargasso Sea. Deep-Sea Res. I, 50, 1445-1480.
    !
    call g_tracer_add_param('si_2_n_max_Lg', phyto(LARGE)%si_2_n_max, 5.0)                  ! mol Si mol N-1
    !
    ! Static large phytoplankton/diatom silicon to nitrogen uptake ratio for a single,
    ! globally constant value instead of a dynamic value.
    !
    call g_tracer_add_param('si_2_n_static', topaz%si_2_n_static, .false.)
    !
    ! Value of the static large phytoplankton/diatom silicon to nitrogen uptake ratio for
    ! a single, globally constant value instead of a dynamic value.  Rather than the
    ! canonical value of 1.0 after culture work of Brzezinski (1985), this 
    ! value is set at the iron-stressed value of 2.0 after Hutchins, D. A. and K. W.
    ! Bruland (1998): Iron-limited diatom growth and Si:N uptake ratios in a coastal
    ! upwelling regime. Nature, 393, 561-564 and Takeda, S. (1998): Influence of iron
    ! variability on nutrient consumption ratio of diatoms in oceanic waters. Nature, 393,
    ! 774-777.
    !
    call g_tracer_add_param('si_2_n_static_Lg', phyto(LARGE)%si_2_n_static, 2.0)            ! mol Si mol N-1
    !
    !----------------------------------------------------------------------------------
    ! Monod half saturation coefficients.  k_PO4 for small phytoplankton was taken
    ! from the minimum of 5 nM observed in Cotner JB, Ammerman JA, Peele ER, Bentzen 
    ! E (1997)  Phosphorus-limited bacterioplankton growth in the Sargasso Sea. 
    ! Aquat Microb Ecol 13:141.), and was assumed to follow the effective size
    ! relationship of Lg/Sm = 3.  k_nh4 was assumed to be equal to k_po4.  k_no3 was
    ! assumed to be a factor of five larger than k_nh4 to allow nh4 inhibition of
    ! no3 uptake after the formulation of Sharada and Yajnik (2005) and Frost and 
    ! Franzen (1992).  k_SiO4 taken from Dugdale, R.C. and Wilkerson, F.P. (1998, 
    ! Silicate regulation of new production in the equatorial Pacific upwelling.
    ! Nature, 391, 270-273).  
    ! k_fed for small phytoplankton was taken from the value for T. weissflogii
    ! of 3 nM in Hudson and Morel (1990; L & O 35, 1002-1020) and the value for large 
    ! assumed to follow the allometric relationship above.
    !----------------------------------------------------------------------------------
    !
    call g_tracer_add_param('k_fed_Di', phyto(DIAZO)%k_fed,  9.0e-9)                        ! mol Fed kg-1
    call g_tracer_add_param('k_fed_Lg', phyto(LARGE)%k_fed,  9.0e-9)                        ! mol Fed kg-1
    call g_tracer_add_param('k_fed_Sm', phyto(SMALL)%k_fed,  3.0e-9)                        ! mol Fed kg-1
    call g_tracer_add_param('k_nh4_Lg', phyto(LARGE)%k_nh4,  6.0e-7)                        ! mol NH4 kg-1
    call g_tracer_add_param('k_nh4_Sm', phyto(SMALL)%k_nh4,  2.0e-7)                        ! mol NH4 kg-1
    call g_tracer_add_param('k_no3_Lg', phyto(LARGE)%k_no3,  6.0e-6)                        ! mol NO3 kg-1
    call g_tracer_add_param('k_no3_Sm', phyto(SMALL)%k_no3,  2.0e-6)                        ! mol NO3 kg-1
    call g_tracer_add_param('k_po4_Di', phyto(DIAZO)%k_po4,  6.0e-7)                        ! mol PO4 kg-1
    call g_tracer_add_param('k_po4_Lg', phyto(LARGE)%k_po4,  6.0e-7)                        ! mol PO4 kg-1
    call g_tracer_add_param('k_po4_Sm', phyto(SMALL)%k_po4,  2.0e-7)                        ! mol PO4 kg-1
!    call g_tracer_add_param('k_sio4_Lg',phyto(LARGE)%k_sio4, 3.0e-6)                       ! mol SiO4 kg-1
!## jgj NOTE: k_sio4_Sm = 1.0e-6 in mom4p1
    call g_tracer_add_param('k_sio4_Lg',phyto(LARGE)%k_sio4, 1.0e-6)                        ! mol SiO4 kg-1

    !-----------------------------------------------------------------------
    ! Iron
    !-----------------------------------------------------------------------
    !
    ! Whether or not to allow mineral ballast dissolution to
    ! return iron to the dissolved phase - a "false" value
    ! assumes that all iron is associated with organic material.
    ! A true value assumes that iron is distributed between
    ! mineral and organic matter by mass leading to a deeper
    ! regeneration length scale.
    !
    call g_tracer_add_param('fe_ballast_assoc', topaz%fe_ballast_assoc, .true.)
    !
    ! Background concentration of iron ligand of 2.0e-9 mol kg-1 taken from Rue, E. L. and K. W.
    ! Bruland (1995) Mar. Chem., 50, 117-138.  Alternatively, a value of 1.0e-9 mol kg-1 can be
    ! taken from Parekh, P., M. J. Follows and E. A. Boyle (2005) Decoupling of iron
    ! and phosphate in the global ocean. Glob. Biogeochem. Cycles, 19, 
    ! doi: 10.1029/2004GB002280.
    !
    call g_tracer_add_param('felig_bkg', topaz%felig_bkg, 2.0e-9)                           ! mol Fe kg-1
    !
    ! Ratio of iron ligand to semilabile and labile don taken from the ratio of
    ! ligand to dissolved organic carbon in deep water.
    !
    call g_tracer_add_param('felig_2_don', topaz%felig_2_don, 2.0e-3 / 40.0 * 106.0 / 16.0) ! mol Fe mol N-1
    !
    ! Iron limitation of the Chl:C, through the def_fe factor, to allow
    ! iron to modulate diazotrophic phytoplankton light utilization
    ! efficiency.  This value is set to a relatively high value after Raven et al.
    !
    !call g_tracer_add_param('k_fe_2_n_Di', phyto(DIAZO)%k_fe_2_n, 24.0e-6 * 106.0 / 16.0)  ! mol Fe mol N-1
    !
    ! Iron limitation of the Chl:C, through the def_fe factor, to allow
    ! iron to modulate small phytoplankton light utilization efficiency.
    ! This value is set to a very low value after Sunda and Huntsman (1995, 
    ! Mar. Chem, 50, 189-206) for Small phytoplankton and is assumed to
    ! have an implicit allometric relationship of eff_size_Lg_2_Sm.
    !
    ! values raised to stimulate expression of iron limitation.
    !
    ! call g_tracer_add_param('k_fe_2_n_Di', phyto(DIAZO)%k_fe_2_n, 60.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    ! call g_tracer_add_param('k_fe_2_n_Lg', phyto(LARGE)%k_fe_2_n, 30.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    ! call g_tracer_add_param('k_fe_2_n_Sm', phyto(SMALL)%k_fe_2_n, 10.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1

    ! per JPD 2009/12/22 Values above were tuned for mom4p0/ESM2.1.  Try using original values from literature
    ! 2010/04/15 jgj: k_fe_2_n_Sm changed to 12e-6*106/16

    call g_tracer_add_param('k_fe_2_n_Di', phyto(DIAZO)%k_fe_2_n, 36.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    call g_tracer_add_param('k_fe_2_n_Lg', phyto(LARGE)%k_fe_2_n, 18.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    call g_tracer_add_param('k_fe_2_n_Sm', phyto(SMALL)%k_fe_2_n, 12.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    !
    ! Maximum Fe:N level where uptake ceases for Small Phytoplankton...
    ! that is, where the phytoplankton get "full" of iron.  This maximum
    ! value was set using the "fe replete" value of 23 umol Fe / mol C for
    ! Synechococcus observed by Porta et al (2003; J. Phycol, 39, 64-73).
    ! Since the Droop function diminishes to zero at half the maximum, the
    ! Droop maximum is set to twice the observed maximum.
    !
    call g_tracer_add_param('fe_2_n_max_Sm', phyto(SMALL)%fe_2_n_max, 46.e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    !
    ! Maximum Fe:N level where uptake ceases for Large Phytoplankton...
    ! that is, where the phytoplankton get "full" of iron.  This maximum
    ! value was set from an "fe replete" value of 333 umol Fe / mol C for the
    ! coast diatom T weissflogii observed by Sunda and Huntsman (1995, 
    ! Mar. Chem, 50, 189-206).
    !
    call g_tracer_add_param('fe_2_n_max_Lg', phyto(LARGE)%fe_2_n_max, 666.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Di', phyto(DIAZO)%fe_2_n_max, 666.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    !
    ! Normalization Fe:N level for phytoplankton set to the level of Fe:N for small phytoplankton under the
    ! assumption of that being the level excluding luxury uptake.
    ! 2010/05/02 jgj: fe_2_n_norm_Lg = fe_2_n_norm_Sm = 50e-6*106/16
    !
    call g_tracer_add_param('fe_2_n_norm_Sm', phyto(SMALL)%fe_2_n_norm, 50.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_norm_Lg', phyto(LARGE)%fe_2_n_norm, 50.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_norm_Di', phyto(DIAZO)%fe_2_n_norm, 100.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    
    ! Ratio of iron influx from bottom sediment boundaries based on nitrogen flux
    !
    call g_tracer_add_param('fe_2_n_sed', topaz%fe_2_n_sed, 100.0e-6 * 106.0 / 16.0)            ! mol Fe mol N-1
    !
    ! Static iron to nitrogen uptake ratio for a single, globally constant value for each
    ! phytoplankton group instead of dynamic values for each.
    !
    call g_tracer_add_param('fe_2_n_static', topaz%fe_2_n_static, .false.)
    !
    ! If static, Diazotrophs are assumed to have a Fe:N after the value of the Fe:C ratio
    ! necessary to achieve half the maximum growth rate of Trichodesmium in
    ! Berman-Frank, I., J. T. Cullen, Y. Shaked, R. M. Sherrell, and P. G. Falkowski
    ! (2001) Iron availability, cellular iron quotas, and nitrogen fixation in
    ! Trichodesmium. Limnol. Oceanogr., 46, 1249-1260.
    !
    !  call g_tracer_add_param('fe_2_n_static_Di', topaz%fe_2_n_static_Di, 30.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    !
    ! If static, large and small phytoplankton are assumed to have an Fe:N after the value of
    ! the Fe:C ratio necessary to achieve half the maximum growth rate of T. oceanica in 
    ! Sunda and Huntsman (1995) Iron uptake and growth limitation in oceanic and coastal 
    ! phytoplankton. Mar. Chem., 50, 189-196.
    !
    !  call g_tracer_add_param('fe_2_n_static_Lg', topaz%fe_2_n_static_Lg, 3.0e-6 * 106.0 / 16.0 ) ! mol Fe mol N-1
    !  call g_tracer_add_param('fe_2_n_static_Sm', topaz%fe_2_n_static_Sm, 3.0e-6 * 106.0 / 16.0 ) ! mol Fe mol N-1
    !
    ! values raised to stimulate expression of iron limitation.
    !
    call g_tracer_add_param('fe_2_n_static_Di', phyto(DIAZO)%fe_2_n_static, 60.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_static_Lg', phyto(LARGE)%fe_2_n_static, 45.0e-6 * 106.0 / 16.0 )! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_static_Sm', phyto(SMALL)%fe_2_n_static, 15.0e-6 * 106.0 / 16.0 )! mol Fe mol N-1
    !
    ! Rate kinetics of iron influx from side boundaries
    !
    call g_tracer_add_param('fe_coast', topaz%fe_coast,1.0e-11 )                            ! mol Fe m kg-1 s-1
    !
    ! Second-order iron scavenging in order to prevent high iron
    ! accumulations in high deposition regions (like the tropical
    ! Atlantic).
    !
    call g_tracer_add_param('kfe_2nd_order', topaz%kfe_2nd_order, 1.0e10/sperd)             ! mol Fe-1 kg s-1
    !
    ! Adsorption rate coefficient for ballast.  This was
    ! set to a low value to prevent iron from accumulating in
    ! the deep ocean and keep a no3-like profile instead of a sio4-like profile.
    !
    call g_tracer_add_param('kfe_bal', topaz%kfe_bal, 0.0/sperd)                            ! g ballast-1 m3 s-1
    !
    ! Desorption rate coefficient.  Set to 0.0068 d-1 after Bacon and Anderson (1982)
    ! J. Geophys. Res., 87, No. C3, 2045-2056. and Clegg and Whittfield (1993)
    ! Deep-Sea Res., 40, 1529-1545 from Thorium.
    !
    call g_tracer_add_param('kfe_des', topaz%kfe_des, 0.0068/sperd)                         ! s-1
    !
    ! Equilibrium constant for (free and inorganically bound) iron binding with organic
    ! ligands taken from Parekh, P., M. J. Follows and E. A. Boyle (2005) Decoupling of
    ! iron and phosphate in the global ocean. Glob. Biogeochem. Cycles, 19, 
    ! doi: 10.1029/2004GB002280.
    !
    call g_tracer_add_param('kfe_eq_lig', topaz%kfe_eq_lig, 1.0e11)                         ! mol lig-1 kg
    !
    ! Adsorption rate coefficient for detrital organic material.  This was
    ! set to a low value to prevent iron from accumulating in
    ! the deep ocean and keep a no3-like profile instead of a sio4-like profile.

    !
    call g_tracer_add_param('kfe_org', topaz%kfe_org, 0.0/sperd)                            ! g org-1 m3 s-1
    !
    !-----------------------------------------------------------------------
    ! Photosynthesis
    !-----------------------------------------------------------------------
    !
    ! Phytoplankton growth altered from Geider et al (1997)
    ! and Moore et al (2002).  Thetamax values
    ! are at the high end in order to account for the additional
    ! iron limitation term.  The factor of 6.022e17 is to convert
    ! from umol to quanta and 2.77e18 to convert from quanta/sec
    ! to Watts given the average energy spectrum for underwater
    ! PAR from the Seabird sensor.  Values of P_C_max are decreased
    ! relative to Geider et al. (1997) by a factor of 4 to account
    ! the difference in reference temperatures and, for Small and Large
    ! increased by a factor of 5 to account for the addition of po4 limitation.
    ! Values of thetamax are increased relative to Geider et al. (1997)
    ! by a factor of 2 to account for the def_fe factor.
    !
    ! alpha is assumed to not have a size relationship
    !
    ! theta_max is assumed to have an implicit allometric relationship of 
    ! eff_size**(2/3) from the surface area to volume relationship.
    !
!    call g_tracer_add_param('alpha_Di', phyto(DIAZO)%alpha,  1.0e-5 * 2.77e18 / 6.022e17)   ! g C g Chl-1 m2 W-1 s-1
!    call g_tracer_add_param('alpha_Lg', phyto(LARGE)%alpha,  2.0e-5 * 2.77e18 / 6.022e17)   ! g C g Chl-1 m2 W-1 s-1
!    call g_tracer_add_param('alpha_Sm', phyto(SMALL)%alpha,  2.0e-5 * 2.77e18 / 6.022e17)   ! g C g Chl-1 m2 W-1 s-1
!
! Values all reset 20% higher to reflect lower light regime in generic TOPAZ.
! High values of 2.4 equal to the maximum values seen in Geider et al., 1997 for
! skeletonema and Microcystis
!
! Value for Di set to half the total to reflect the relative assembly
! limit of r_assem_max_Di/r_assem_Lg = p_2_max_Di/p_2_n_max_Lg = (1-0.1-0.5)/(1-0.2)
!!
    call g_tracer_add_param('alpha_Di', phyto(DIAZO)%alpha,  1.2e-5 * 2.77e18 / 6.022e17)   ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('alpha_Lg', phyto(LARGE)%alpha,  2.4e-5 * 2.77e18 / 6.022e17)   ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('alpha_Sm', phyto(SMALL)%alpha,  2.4e-5 * 2.77e18 / 6.022e17)   ! g C g Chl-1 m2 W-1 s-1
    call g_tracer_add_param('kappa_eppley', topaz%kappa_eppley, 0.063)                      ! deg C-1
!    call g_tracer_add_param('P_C_max_Di', phyto(DIAZO)%P_C_max, 0.6e-5)                     ! s-1
!
! Value for Di set to half the total to reflect the relative assembly
! limit of r_assem_max_Di/r_assem_Lg = p_2_max_Di/p_2_n_max_Lg = (1-0.1-0.5)/(1-0.2)
!
    call g_tracer_add_param('P_C_max_Di', phyto(DIAZO)%P_C_max, 0.75e-5)                     ! s-1
    call g_tracer_add_param('P_C_max_Lg', phyto(LARGE)%P_C_max, 1.5e-5)                     ! s-1
    call g_tracer_add_param('P_C_max_Sm', phyto(SMALL)%P_C_max, 1.5e-5)                     ! s-1
    call g_tracer_add_param('thetamax_Di', phyto(DIAZO)%thetamax, 0.04)                     ! g Chl g C-1
    call g_tracer_add_param('thetamax_Lg', phyto(LARGE)%thetamax, 0.06)                     ! g Chl g C-1
    call g_tracer_add_param('thetamax_Sm', phyto(SMALL)%thetamax, 0.04)                     ! g Chl g C-1
    call g_tracer_add_param('thetamin', topaz%thetamin, 0.001)                              ! g Chl g C-1
    call g_tracer_add_param('thetamin_nolim', topaz%thetamin_nolim, 0.01)                   ! g Chl g C-1
    call g_tracer_add_param('zeta', topaz%zeta, 0.1)                                        ! dimensionless
    !
    ! Diazotrophs are assumed to be inhibited by nitrate after Holl and Montoya 
    ! (2005;  Interactions between nitrate uptake and nitrogen fixation in 
    ! continuous cultures of the marine diazotroph trichodesmium cyanobacteria); 
    ! J. Phycol, 41, 1178-1183).
    !
    call g_tracer_add_param('k_n_inhib_Di', topaz%k_n_inhib_Di, 7.0e-6)                     ! mol NO3 kg-1
    !
    ! Diazotrophs are also assumed to be inhibited by oxygen after Stewart and Pearson (1970;
    ! Effects of aerobic and anaerobic conditions on growth and metabolism of blue-green
    ! algae. Proc. Soc. Lond. B., 175, 293-311) and Berman-Frank et al (2005; Inhibition of
    ! nitrogenase by oxygen in marine cyanobacteria controls the global nitrogen and 
    ! oxygen cycles. Biogeosciences Discussions, 2, 261-273).
    !
    call g_tracer_add_param('o2_inhib_Di_pow', topaz%o2_inhib_Di_pow, 4.0)                  ! mol O2-1 m3
    call g_tracer_add_param('o2_inhib_Di_sat', topaz%o2_inhib_Di_sat, 3.0e-4)               ! mol O2 kg-1
    !
    ! Chl:C response rate constant for phytoplankton calibrated to 1 d-1
    ! after Owens et al (1980, Diel Periodicity in cellular Chlorophyll
    ! content of marine diatoms, Mar. Biol, 59, 71-77).
    !
    call g_tracer_add_param('gamma_irr_mem', topaz%gamma_irr_mem, 1.0 / sperd)              ! s-1
    !
    !-----------------------------------------------------------------------
    ! Grazing
    !-----------------------------------------------------------------------
    !
    ! Values of fractional detritus production from the global
    !  synthesis of Dunne et al. (submitted)
    !
    call g_tracer_add_param('fdet0_Lg', phyto(LARGE)%fdet0, 0.93)                           ! dimensionless
    call g_tracer_add_param('fdet0_Sm', phyto(SMALL)%fdet0, 0.18)                           ! dimensionless
    !
    ! Rate constant for remineralization of heterotrophic biomass
    ! 
    call g_tracer_add_param('gamma_nhet', topaz%gamma_nhet, 1.0 / (30.0 * sperd))           ! s-1
    !
    ! Dissolution of sio2 was set as a temperature-dependent
    ! fraction of grazed material to be roughly in line with
    ! the work of Kamatani (1982)
    !
    call g_tracer_add_param('k_diss_sio2', topaz%k_diss_sio2, 3.0 )     ! s-1
    !
    ! Half saturation oxygen concentration for oxic remineralization rate.
    !
    call g_tracer_add_param('k_o2', topaz%k_o2, 20.0e-6)                                    ! mol O2 kg-1
    !
    ! Temperature-dependence of fractional detritus production
    ! from the global synthesis of Dunne et al. (submitted)
    !
    call g_tracer_add_param('kappa_remin', topaz%kappa_remin, -0.032)                       ! deg C-1
    !
    ! T=0 phytoplankton specific grazing rate from the global
    ! synthesis of Dunne et al. (2005)
    !
    call g_tracer_add_param('lambda0', topaz%lambda0, 0.19 / sperd)                         ! s-1
    !
    ! Minimum oxygen concentration for oxic remineralization.
    ! this is necessary for both numerical stability and to
    ! qeue the switch to denitrification
    !
    call g_tracer_add_param('o2_min', topaz%o2_min, 1.0 * 1.0e-06)                          ! mol O2 kg-1
    !
    ! Pivot phytoplankton concentration for grazing-based
    ! variation in ecosystem structure from the global
    ! synthesis of Dunne et al. (submitted)
    !
    call g_tracer_add_param('P_star', topaz%P_star, 1.9e-6 * 16.0 / 106.0)                  ! mol N kg-1
    !
    ! Minimum phytoplankton concentration or grazing.  This is
    ! necessary for numerical stability.
    ! 
    call g_tracer_add_param('phyto_min', topaz%phyto_min, 1.0e-10)                          ! mol N kg-1
    !
    ! Fraction of phytoplankton grazing and dom consumption eventually going to NH4
    ! that is temporarily stored in heterotrophic biomass
    ! 
    call g_tracer_add_param('phi_nhet', topaz%phi_nhet, 0.50)                               ! dimensionless
    !
    ! SiO2 dissolution is set to globally dissolve 50% after Nelson et al. (1995) 
    ! through grazing.  The temperature functionality is set to a combination of
    ! stoichiometry and Eppley temperature formulation to give roughly the range of 
    ! observations in Kamatani (1982) with respect to frustrule thickness and
    ! temperature by utilizing the inverse of the Eppley temperature
    ! functionality and a normalization to stoichiometry ( q_si_2_n_diss).
    ! The value of q_si_2_n_diss was set so as to simultaneously reproduce the low
    ! silicon export efficiencies (~0.1) observed in the equatorial Pacific by
    ! Blain et al. (1997, DSR I; Dunne et al., 1999, GBC) and high export efficiencies
    ! of ~0.64 observed in the Southern Ocean by Bzrezinski et al., 2001, retaining
    ! a ~0.5 global average after Nelson et al. (1995).
    !
    call g_tracer_add_param('q_si_2_n_diss', topaz%q_si_2_n_diss, 3.0)                      ! mol mol Si mol N
    !
    ! Organic matter protection by mineral - after Klaas and
    ! Archer (2002)
    !
    call g_tracer_add_param('rpcaco3', topaz%rpcaco3, 0.070/12*16/106.0*100)                ! mol N mol Ca-1
    call g_tracer_add_param('rplith',  topaz%rplith,  0.065/12*16/106.0)                    ! mol N g lith-1
    call g_tracer_add_param('rpsio2',  topaz%rpsio2,  0.026/12*16/106.0*60)                 ! mol N mol Si-1
    !
    !-----------------------------------------------------------------------
    ! Remineralization length scales
    !-----------------------------------------------------------------------
    !

    !
    ! Value of gamma_ndet to approximate upper e-folding of the globally-tuned
    ! "Martin curve" used in the OCMIP-II Biotic configuration of (z/75)^-0.9
    ! that gives a value of exp(-1) at 228 m from 75 m for an e-folding scale
    ! of 188 m.
    !
    call g_tracer_add_param('gamma_ndet',  topaz%gamma_ndet, topaz%wsink / 188.0 )          ! s-1
    !
    ! Cadet_arag dissolution length scale fit to Betzer et al 1984 data from
    ! 20n - 46n, 170e attenuation of aragonite flux between 900-2200 m of 47.51%
    !
    call g_tracer_add_param('gamma_cadet_arag',  topaz%gamma_cadet_arag, topaz%wsink / 760.0)        ! s-1
    !
    ! Cadet_Calc dissolution rate constant: Most deep traps show little correlation
    ! with depth, suggesting little water column dissolution in much of the
    ! ocean.  Values were calibrated to get the observed 67% transfer efficiency at
    ! Ocean Weather Station PAPA in the North Pacific between 1000-3800 m
    ! equivalent to a (1-omega)-modulated dissolution length scale of 1343 m
    ! based on the Honjo data on the JGOFS database assuming an omega of 0.81 from
    ! GLODAP and converted to a dissolution rate constant assuming a sinking
    ! velocity of 20 m d-1.
    !
    call g_tracer_add_param('gamma_cadet_calc',  topaz%gamma_cadet_calc, topaz%wsink / 1343.0)        ! s-1
    !
    ! Sidet dissolution rate constant assuming a dissolution length scale of 2000 m
    ! consistent with Gnanadesikan (2000).
    !
    call g_tracer_add_param('gamma_sidet',  topaz%gamma_sidet, topaz%wsink / 2000.0 )       ! s-1
    !
    !-----------------------------------------------------------------------
    ! Dissolved Organic Material
    !-----------------------------------------------------------------------
    !
    !
    ! Dissolved Organic Material remineralization rate constants
    ! and fractional production ratios, all to be consistent
    ! with the work of Abell et al. (2000, Distributions of TOP, TON and TOC
    ! in the North pacific subtropical gyre: Implications for nutrient supply
    ! in the surface ocean and remineralization in the upper thermocline,
    ! J. Mar. Res., 58, 203-222) and DOC and DON data provided by
    ! Dennis Hansell (personal communication).  Assuming refractory/deep
    ! concentrations of DOC=42 uM, DON=1.8 uM, and DOP=0.0 uM, we allow the
    ! semilabile and labile pools to have RKR C:N and reproduce observed
    ! surface expression in north Pacific subtropical gyre (DOC=72 uM, 
    ! DON=6 uM, and DOP=0.2 uM) and depth penetration. 
    !
    ! Warning: phi_sdon + phi_ldon should be < 1.0.
    !
    !
!    call g_tracer_add_param('gamma_sdon',  topaz%gamma_sdon, 1.0 / (18.0 *365.0 * sperd))   ! s-1
!    call g_tracer_add_param('gamma_sdop',  topaz%gamma_sdop, 1.0 / (4.0 *365.0 * sperd))    ! s-1
!## jgj NOTE: to match mom4p1/gold
    call g_tracer_add_param('gamma_sdon',  topaz%gamma_sdon, 1.0 / (18.0 * spery))          ! s-1
    call g_tracer_add_param('gamma_sdop',  topaz%gamma_sdop, 1.0 / (4.0 * spery))           ! s-1
    ! 2010/04/15 jgj: phi_sdon/phi_sdop changed to 0.025/0.050 respectively 
    call g_tracer_add_param('phi_sdon'  ,  topaz%phi_sdon, 0.025)                            ! dimensionless
    call g_tracer_add_param('phi_sdop'  ,  topaz%phi_sdop, 0.050)                            ! dimensionless
    !
    ! The remineralization rate constant for labile DOP (bio_tau_ldon)
    ! was set to 3 months after Archer et al. (1997, GBC, 11, 435-452).
    ! The fraction going to labile DOC was inspired by data-model
    ! comparisons to Libby and Wheeler (1997, Deep-Sea Res. I, 44, 345-361)
    !
    call g_tracer_add_param('gamma_ldon',  topaz%gamma_ldon, 1.0 / (90.0 * sperd))          ! s-1
    call g_tracer_add_param('phi_ldon'  ,  topaz%phi_ldon, 0.06)                            ! dimensionless
    !
    !-----------------------------------------------------------------------
    ! Miscellaneous
    !-----------------------------------------------------------------------
    !
    ! Flag to perform extra calculations for online diagnostics
    !
    call g_tracer_add_param('do_extra_diag_calcs', topaz%do_extra_diag_calcs, .true.)
    !
    ! Prevent incorporation of changes since CMIP5 version of ESM2
    call g_tracer_add_param('reproduce_esm2', topaz%reproduce_esm2, .false.)
    !
    ! Nitrification rate constant assumed to be light-limited with an inhibition
    ! factor.  gamma_nitrif was tuned to reproduce the scaling observed in Ward et
    ! al. (1982; Microbial nitrification rates in the  primary nitrite maximum off
    ! southern California, Deep-Sea Res., 29, 247-255), and irr_inhibit was tuned to
    ! reproduce Olson (1981; Differential photoinhibition of marine nitrifying
    ! bacteria: a possible mechanism for the formulation of the primary nitrite
    ! maximum, J. Mar. Res., 39, 227-238).
    !
    call g_tracer_add_param('gamma_nitrif',  topaz%gamma_nitrif, 1.0 / (30.0 * sperd))      ! s-1
!    call g_tracer_add_param('irr_inhibit',  topaz%irr_inhibit, 2.0)                         ! m2 W-1
!   jpd - 12/04/09  irradiance functionality was lowered to reduce nitrification in surface waters.
    call g_tracer_add_param('irr_inhibit',  topaz%irr_inhibit, 0.1)                         ! m2 W-1
    !
    ! Scavenging rate coefficient for lithogenic material relative to large
    ! phytoplankton concentration via large phytoplankton grazing.
    !
    call g_tracer_add_param('phi_lith' ,  topaz%phi_lith, 0.002)                          ! dimensionless
    !
    ! Adsorption rate coefficient for lithogenic material onto sinking material.
    ! This was set to a small but non-zero to prevent lithogenic
    ! material from accumulating in the deep ocean.
    !
    !call g_tracer_add_param('k_lith',  topaz%k_lith, 10.0/spery )                           ! s-1
!## jgj NOTE: to match mom4p1/gold
    call g_tracer_add_param('k_lith',  topaz%k_lith, 1e-6/sperd )                           ! s-1
    call g_tracer_add_param('z_sed',  topaz%z_sed, 0.1 )                                    ! m

    !
    ! Constant for productivity mask
    !
    call g_tracer_add_param('bio_tau',  topaz%bio_tau, 30.0 * sperd)                        ! s

    call g_tracer_add_param('r_bio_tau',  topaz%r_bio_tau,1.0 / topaz%bio_tau)


    call g_tracer_end_param_list(package_name)
    !===========
    !Block Ends: g_tracer_add_param
    !===========



  end subroutine user_add_params

  !
  !   This is an internal sub, not a public interface.
  !   Add all the tracers to be used in this module. 
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list


    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'

    !
    !Add here only the parameters that are required at the time of registeration 
    !(to make flux exchanging Ocean tracers known for all PE's) 
    !
    call g_tracer_start_param_list(package_name)
    !
    call g_tracer_add_param('htotal_in', topaz%htotal_in, 1.0e-08)
    !
    ! Sinking velocity of detritus: a value of 20 m d-1 is consistent with a characteristic sinking
    ! velocity of 100 m d-1 of marine aggregates and a disaggregation rate constant
    ! of 5 d-1 in the surface ocean (Clegg and Whitfield, 1992; Dunne, 1999).  Alternatively, 100 m d-1 
    ! is more in line with the deep water synthesis of Berelson (2002; Particle settling rates increase
    ! with depth in the ocean, DSR-II, 49, 237-252).
    !
    call g_tracer_add_param('wsink',  topaz%wsink, 100.0 / sperd)                             ! m s-1

    call g_tracer_add_param('ice_restart_file'   , topaz%ice_restart_file   , 'ice_topaz.res.nc')
    call g_tracer_add_param('ocean_restart_file' , topaz%ocean_restart_file , 'ocean_topaz.res.nc' )
    call g_tracer_add_param('IC_file'       , topaz%IC_file       , '')
    !
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file    = topaz%ice_restart_file,&
         ocean_restart_file  = topaz%ocean_restart_file )

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
    !       Cadet_arag (Sinking detrital/particulate Aragonite CaCO3)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cadet_arag',         &
         longname   = 'Detrital Aragonite', &
         units      = 'mol/kg',             &
         prog       = .true.,               &
         sink_rate  = topaz%wsink,          &
         btm_reservoir = .true.             )
    !
    !       Cadet_Calc (Sinking detrital/particulate Calcite CaCO3)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cadet_calc',       &
         longname   = 'Detrital Calcite', &
         units      = 'mol/kg',           &
         prog       = .true.,             &
         sink_rate  = topaz%wsink,        &
         btm_reservoir = .true.           )
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
         flux_gas_restart_file  = 'ocean_topaz_airsea_flux.res.nc',    &
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
         sink_rate  = topaz%wsink,     &
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
         longname   = 'Large Iron', &
         units      = 'mol/kg',     &
         prog       = .true.        )
    !
    !       Small Fe(Iron in small phytoplankton to allow for variable Fe:N ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fesm',       &
         longname   = 'Small Iron', &
         units      = 'mol/kg',     &
         prog       = .true.        )
    !
    !       LDON (Labile organic nitrogen; assumed to have a fixed, Redfield, Ketchum
    !       and Richards (1963) C:N:P ratio)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ldon',           &
         longname   = 'labile DON',     &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_runoff= .true.,           &
         flux_param = (/ 14.0067e-03 /) )
    !
    !       LITH (Lithogenic aluminosilicate particles)
    !
    call g_tracer_add(tracer_list,package_name,     &
         name       = 'lith',                       &
         longname   = 'Lithogenic Aluminosilicate', &
         units      = 'mol/kg',                     &
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
         longname   = 'Detrital Lith',   &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         sink_rate  = topaz%wsink, &
         btm_reservoir = .true.    )
    !
    !    Ndet (Sinking detrital/particulate Nitrogen)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndet',      &
         longname   = 'Detrital Nitrogen',      &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         sink_rate  = topaz%wsink, &
         btm_reservoir = .true.    )
    !
    !    NDi (assumed to be facultative N2-fixers, with a variable N:P ratio
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndi',                 &
         longname   = 'Diazotroph Nitrogen', &
         units      = 'mol/kg',              &
         prog       = .true.                 )

    !
    !       Heterotrophic N (assumed to be the sum of heterotrphic bacteria,
    !            microzooplankton and mesozooplankton - effectively a biomass
    !            storage reservoir for nutrients and carbon; assumed to have 
    !            a fixed C:N:P ratio)
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'nhet',                   &
         longname   = 'Heterotrophic Nitrogen', &
         units      = 'mol/kg',                 &
         prog       = .true.                    )
    !
    !    NLg (assumed to be a dynamic combination of diatoms and other 
    !         eukaryotes all effectively greater than 5 um in diameter,
    !         and having a fixed C:N ratio)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nlg',            &
         longname   = 'Large Nitrogen', &
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
         flux_runoff= .true.,            &
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
    !       NSm (Nitrogen in picoplankton and nanoplankton - effectively less
    !            than 5 um in diameter and having a fixed C:N:P ratio))
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nsm',            &
         longname   = 'Small Nitrogen', &
         units      = 'mol/kg',         &
         prog       = .true.            )
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
         flux_gas_restart_file  = 'ocean_topaz_airsea_flux.res.nc',    &
         flux_bottom= .true.             )
    !
    !    Pdet (Sinking detrital/particulate Phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,         &
         name       = 'pdet',                           &
         longname   = 'Detrital Phosphorus', &
         units      = 'mol/kg',                         &
         prog       = .true.,                           &
         sink_rate  = topaz%wsink,                      &
         btm_reservoir = .true.    )
    !
    !    PDi (Phosphorus in diaz. phytoplankton for variable N:P ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'pdi',             &
         longname   = 'Diaz Phosphorus', &
         units      = 'mol/kg',          &
         prog       = .true.             )
    !
    !    PLg (Phosphorus in large phytoplankton for variable N:P ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'plg',              &
         longname   = 'Large Phosphorus', &
         units      = 'mol/kg',           &
         prog       = .true.              )
    !
    !       PO4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'po4',       &
         longname   = 'Phosphate', &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         flux_bottom= .true.       )
    !
    !    PSm (Phosphorus in small phytoplankton for variable N:P ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'psm',              &
         longname   = 'Small Phosphorus', &
         units      = 'mol/kg',           &
         prog       = .true.              )
    !
    !       SDON (Semilabile dissolved organic nitrogen for variable N:P ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sdon',           &
         longname   = 'Semilabile DON', &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !       SDOP (Semilabile dissolved organic phosphorus for variable N:P ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sdop',           &
         longname   = 'Semilabile DOP', &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !       Sidet (Sinking detrital/particulate Silicon)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sidet',            &
         longname   = 'Detrital Silicon', &
         units      = 'mol/kg',           &
         prog       = .true.,             &
         sink_rate  = topaz%wsink,        &
         btm_reservoir = .true.    )
    !
    !    SiLg (Silicon in large phytoplankton for variable Si:N ratios
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'silg',          &
         longname   = 'Large Silicon', &
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
!    !
!    !       Passive tracer for testing purposes
!    !
!    call g_tracer_add(tracer_list,package_name,&
!         name       = 'topaz_passive',        &
!         longname   = 'Topaz Passive Tracer', &
!         units      = 'mol/kg',               &
!         const_init_value = 1.0,              &
!         prog       = .true.                  )
!    call g_tracer_add(tracer_list,package_name,       &
!         name       = 'topaz_sinking_passive',        &
!         longname   = 'Topaz Passive Sinking Tracer', &
!         units      = 'mol/kg',                       &
!         const_init_value = 1.0,                      &
!         sink_rate  = topaz%wsink,                    &
!         prog       = .true.                          )

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
    call g_tracer_add(tracer_list,package_name,      &
         name       = 'cadet_arag_btf',              &
         longname   = 'Aragonite flux to Sediments', &
         units      = 'mol m-2 s-1',                 &
         prog       = .false.                        )
    !
    !      cadet_calc_btf (Calcite flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,    &
         name       = 'cadet_calc_btf',            &
         longname   = 'Calcite flux to Sediments', &
         units      = 'mol m-2 s-1',               &
         prog       = .false.                      )
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
         init_value = topaz%htotal_in         )
    !
    !       Irr_mem (Irradiance Memory)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'irr_mem',           &
         longname   = 'Irradiance memory', &
         units      = 'Watts/m^2',         &
         prog       = .false.              )


  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_TOPAZ_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_TOPAZ_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_TOPAZ_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_update_from_coupler'

    real, dimension(:,:)  ,pointer    :: stf_alk,dry_no3,wet_no3

    !
    ! NO3 has deposition, river flux, and negative deposition contribution to alkalinity
    !
    call g_tracer_get_pointer(tracer_list,'no3','drydep',dry_no3)
    call g_tracer_get_pointer(tracer_list,'no3','wetdep',wet_no3)
    call g_tracer_get_pointer(tracer_list,'alk','stf',stf_alk)

    stf_alk = stf_alk - dry_no3 - wet_no3 ! update 'tracer%stf' thru pointer

    return
  end subroutine generic_TOPAZ_update_from_coupler

  ! <SUBROUTINE NAME="generic_TOPAZ_update_from_bottom">
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracers have bottom fluxes and reservoirs. 
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_TOPAZ_update_from_bottom(tracer_list,dt, tau, model_time) 
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
  subroutine generic_TOPAZ_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    !---------------------------------------------------------------------
    ! Get bottom fluxes and reset bottom flux boundary condition - Fedet and Lithdet are lost through
    ! the bottom, while the rest are not.
    !---------------------------------------------------------------------

    call g_tracer_get_values(tracer_list,'cadet_arag'  ,'btm_reservoir', topaz%fcadet_arag_btm,isd,jsd)
    topaz%fcadet_arag_btm = topaz%fcadet_arag_btm /dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = topaz%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (topaz%id_fcadet_arag_btm .gt. 0)           &
         used = send_data(topaz%id_fcadet_arag_btm,    topaz%fcadet_arag_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc'  ,'btm_reservoir', topaz%fcadet_calc_btm,isd,jsd)
    topaz%fcadet_calc_btm = topaz%fcadet_calc_btm /dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = topaz%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (topaz%id_fcadet_calc_btm .gt. 0)           &
         used = send_data(topaz%id_fcadet_calc_btm,    topaz%fcadet_calc_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'fedet'  ,'btm_reservoir', topaz%ffedet_btm,isd,jsd)
    topaz%ffedet_btm = topaz%ffedet_btm /dt
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (topaz%id_ffedet_btm .gt. 0)           &
         used = send_data(topaz%id_ffedet_btm,    topaz%ffedet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir', topaz%flithdet_btm,isd,jsd)
    topaz%flithdet_btm = topaz%flithdet_btm /dt
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = topaz%flithdet_btm(:,:)
    if (topaz%id_flithdet_btm .gt. 0)           &
         used = send_data(topaz%id_flithdet_btm,    topaz%flithdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'ndet'  ,'btm_reservoir', topaz%fndet_btm,isd,jsd)
    topaz%fndet_btm = topaz%fndet_btm /dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = topaz%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (topaz%id_fndet_btm .gt. 0)           &
         used = send_data(topaz%id_fndet_btm,    topaz%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet'  ,'btm_reservoir', topaz%fpdet_btm,isd,jsd)
    topaz%fpdet_btm = topaz%fpdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = topaz%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (topaz%id_fpdet_btm .gt. 0)           &
         used = send_data(topaz%id_fpdet_btm,    topaz%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet'  ,'btm_reservoir', topaz%fsidet_btm,isd,jsd)
    topaz%fsidet_btm = topaz%fsidet_btm /dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = topaz%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (topaz%id_fsidet_btm .gt. 0)           &
         used = send_data(topaz%id_fsidet_btm,    topaz%fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine generic_TOPAZ_update_from_bottom

  ! <SUBROUTINE NAME="generic_TOPAZ_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This is the subroutine to contain most of the biogeochemistry for calculating the 
  !   interaction of tracers with each other and with outside forcings.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_TOPAZ_update_from_source(tracer_list,Temp,Salt,dzt,hblt_depth,&
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
  subroutine generic_TOPAZ_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,hblt_depth,&
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


    character(len=fm_string_len), parameter :: sub_name = 'generic_TOPAZ_update_from_source'
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,i,j,k,kblt,n,k_100
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast,grid_kmt

    logical :: used, first
    integer :: nb
    real :: r_dt
    real :: feprime
    real :: graz_lg_terms
    real :: jgraz_n
    real :: jgraz_p
    real :: jgraz_fe
    real :: jprod_di_tot2nterm
    real :: log_btm_flx
    real :: P_C_m
    real :: p_lim_nhet
    real :: r_assem,  r_photo,   r_uptake
    real :: thetamin_term
    real :: TK, PRESS, PKSPA, PKSPC
    real :: tmp_hblt, tmp_Irrad, tmp_irrad_ML,tmp_opacity
    real :: drho_dzt
    real, dimension(:), Allocatable   :: tmp_irr_band 
    real, dimension(:,:), Allocatable :: rho_dzt_100
    real, dimension(NUM_PHYTO)        :: p_2_n_max
    real, dimension(:,:)   ,pointer   :: runoff_flux_alk,runoff_flux_dic
    real, dimension(:,:)   ,pointer   :: runoff_flux_fed,wet_fed,dry_fed
    real, dimension(:,:)   ,pointer   :: runoff_flux_lith
    real, dimension(:,:)   ,pointer   :: runoff_flux_no3,wet_no3,dry_no3
    real, dimension(:,:)   ,pointer   :: runoff_flux_nh4,wet_nh4,dry_nh4

    r_dt = 1.0 / dt

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)


    !Get necessary fields
    call g_tracer_get_values(tracer_list,'htotal','field', topaz%f_htotal,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'po4'   ,'field', topaz%f_po4,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'sio4'  ,'field', topaz%f_sio4,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'alk'   ,'field', topaz%f_alk,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dic'   ,'field', topaz%f_dic  ,isd,jsd,ntau=tau)

    !---------------------------------------------------------------------
    !Calculate co3_ion
    !Also calculate co2 fluxes csurf and alpha for the next round of exchnage
    !---------------------------------------------------------------------
   
    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       topaz%htotallo(i,j) = topaz%htotal_scale_lo * topaz%f_htotal(i,j,k)
       topaz%htotalhi(i,j) = topaz%htotal_scale_hi * topaz%f_htotal(i,j,k)
    enddo; enddo ; !} i, j
 
    if (topaz%reproduce_esm2) then
       call FMS_ocmip2_co2calc_old(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         topaz%f_dic(:,:,k),                          &
         topaz%f_po4(:,:,k),                          &  
         topaz%f_sio4(:,:,k),                         &
         topaz%f_alk(:,:,k),                          &
         topaz%htotallo, topaz%htotalhi,&
                                !InOut
         topaz%f_htotal(:,:,k),                       & 
                                !OUT
         co2star=topaz%co2_csurf(:,:), alpha=topaz%co2_alpha(:,:), &
         pCO2surf=topaz%pco2_csurf(:,:), &
         co3_ion=topaz%f_co3_ion(:,:,k))
    else
    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         topaz%f_dic(:,:,k),                          &
         topaz%f_po4(:,:,k),                          &  
         topaz%f_sio4(:,:,k),                         &
         topaz%f_alk(:,:,k),                          &
         topaz%htotallo, topaz%htotalhi,&
                                !InOut
         topaz%f_htotal(:,:,k),                       & 
                                !OUT
         co2star=topaz%co2_csurf(:,:), alpha=topaz%co2_alpha(:,:), &
         pCO2surf=topaz%pco2_csurf(:,:), &
         co3_ion=topaz%f_co3_ion(:,:,k))
    endif

    if(topaz%reproduce_esm2)then
      do k = 2, nk
        do j = jsc, jec ; do i = isc, iec  !{
          topaz%htotallo(i,j) = topaz%htotal_scale_lo * topaz%f_htotal(i,j,k)
          topaz%htotalhi(i,j) = topaz%htotal_scale_hi * topaz%f_htotal(i,j,k)
        enddo; enddo ; !} i, j
        call FMS_ocmip2_co2calc_old(CO2_dope_vec,grid_tmask(:,:,k),&
            Temp(:,:,k), Salt(:,:,k),                    &
            topaz%f_dic(:,:,k),                          &
            topaz%f_po4(:,:,k),                          &  
            topaz%f_sio4(:,:,k),                         &
            topaz%f_alk(:,:,k),                          &
            topaz%htotallo, topaz%htotalhi,&
                                !InOut
            topaz%f_htotal(:,:,k),                       & 
                                !OUT
            co3_ion=topaz%f_co3_ion(:,:,k))
      enddo
    else
    do k = 2, nk
       do j = jsc, jec ; do i = isc, iec  !{
          topaz%htotallo(i,j) = topaz%htotal_scale_lo * topaz%f_htotal(i,j,k)
          topaz%htotalhi(i,j) = topaz%htotal_scale_hi * topaz%f_htotal(i,j,k)
       enddo; enddo ; !} i, j
       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
            Temp(:,:,k), Salt(:,:,k),                    &
            topaz%f_dic(:,:,k),                          &
            topaz%f_po4(:,:,k),                          &  
            topaz%f_sio4(:,:,k),                         &
            topaz%f_alk(:,:,k),                          &
            topaz%htotallo, topaz%htotalhi,&
                                !InOut
            topaz%f_htotal(:,:,k),                       & 
                                !OUT
            co3_ion=topaz%f_co3_ion(:,:,k))
    enddo
    endif

    call g_tracer_set_values(tracer_list,'htotal','field',topaz%f_htotal  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion','field',topaz%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'dic','alpha',topaz%co2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','csurf',topaz%co2_csurf    ,isd,jsd)



    !---------------------------------------------------------------------
    ! Get positive tracer concentrations
    !---------------------------------------------------------------------

    call g_tracer_get_values(tracer_list,'cadet_arag' ,'field',topaz%f_cadet_arag ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'cadet_calc' ,'field',topaz%f_cadet_calc ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fed'    ,'field',topaz%f_fed      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fedet'  ,'field',topaz%f_fedet    ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ldon'   ,'field',topaz%f_ldon     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'lith'   ,'field',topaz%f_lith     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'lithdet','field',topaz%f_lithdet  ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ndet'   ,'field',topaz%f_ndet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nh4'    ,'field',topaz%f_nh4      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nhet'   ,'field',topaz%f_nhet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'no3'    ,'field',topaz%f_no3      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'o2'     ,'field',topaz%f_o2       ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'pdet'   ,'field',topaz%f_pdet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'po4'    ,'field',topaz%f_po4      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sdon'   ,'field',topaz%f_sdon     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sdop'   ,'field',topaz%f_sdop     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sidet'  ,'field',topaz%f_sidet    ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sio4'   ,'field',topaz%f_sio4     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'silg'   ,'field',topaz%f_silg     ,isd,jsd,ntau=tau,positive=.true.)

    call g_tracer_get_values(tracer_list,'fedi'   ,'field',phyto(DIAZO)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'felg'   ,'field',phyto(LARGE)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fesm'   ,'field',phyto(SMALL)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ndi'    ,'field',phyto(DIAZO)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nlg'    ,'field',phyto(LARGE)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nsm'    ,'field',phyto(SMALL)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'pdi'    ,'field',phyto(DIAZO)%f_p(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'plg'    ,'field',phyto(LARGE)%f_p(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'psm'    ,'field',phyto(SMALL)%f_p(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)

    call g_tracer_get_values(tracer_list,'cased'  ,'field',topaz%f_cased    ,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'co3_ion','field',topaz%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'cadet_arag_btf','field',topaz%f_cadet_arag_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'cadet_calc_btf','field',topaz%f_cadet_calc_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'lithdet_btf','field',topaz%f_lithdet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'ndet_btf','field',topaz%f_ndet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'pdet_btf','field',topaz%f_pdet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'sidet_btf','field',topaz%f_sidet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'irr_mem' ,'field',topaz%f_irr_mem  ,isd,jsd,ntau=1)

    topaz%zt = 0.0

    n=DIAZO
    p_2_n_max(n) = phyto(DIAZO)%p_2_n_assem * (1.0 - phyto(n)%r_other_min - phyto(n)%r_nfix) +      &
       topaz%p_2_n_RKR * phyto(n)%r_other_min
    do n = 2, NUM_PHYTO   
       p_2_n_max(n) = phyto(n)%p_2_n_assem * (1.0 - phyto(n)%r_other_min) +                         &
          topaz%p_2_n_RKR * phyto(n)%r_other_min
    enddo


    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec  
       do n= 1, NUM_PHYTO   !{
          phyto(n)%q_fe_2_n(i,j,k) = min(phyto(n)%fe_2_n_max, max(0.0, phyto(n)%f_fe(i,j,k) /       &
             max(epsln,phyto(n)%f_n(i,j,k))))
          phyto(n)%q_p_2_n(i,j,k) = min(p_2_n_max(n), max(0.0, phyto(n)%f_p(i,j,k) /                &
             max(epsln,phyto(n)%f_n(i,j,k))))
       enddo !} n
       !
       !
       !---------------------------------------------------------------------
       !     Phytoplankton growth and grazing through the water column
       !---------------------------------------------------------------------
       !
       !-----------------------------------------------------------------------
       ! Calculate nutrient limitation terms
       ! 1d-30 added to avoid divide by zero where necessary
       !-----------------------------------------------------------------------
       !
       !-----------------------------------------------------------------------
       ! N limitation with NH4 inhibition after Frost and Franzen (1992)
       !-----------------------------------------------------------------------
       !
       do n= 2, NUM_PHYTO   !{
          phyto(n)%no3lim(i,j,k) = topaz%f_no3(i,j,k) / ((phyto(n)%k_no3 + topaz%f_no3(i,j,k)) *&
             (1.0 + topaz%f_nh4(i,j,k) / phyto(n)%k_nh4))
          phyto(n)%nh4lim(i,j,k) = topaz%f_nh4(i,j,k) / (phyto(n)%k_nh4 + topaz%f_nh4(i,j,k))
       enddo !} n
       !
       !-----------------------------------------------------------------------
       ! The rest are straight Michaelis Menten
       !-----------------------------------------------------------------------
       !
       ! phyto(LARGE)%k_sio4 was topaz(n)%k_sio4_Sm in MOM_TOPAZ code.
       phyto(LARGE)%silim(i,j,k) = topaz%f_sio4(i,j,k) / (phyto(LARGE)%k_sio4 + topaz%f_sio4(i,j,k))
       do n= 1, NUM_PHYTO   !{
          phyto(n)%po4lim(i,j,k) = topaz%f_po4(i,j,k) / (phyto(n)%k_po4 + topaz%f_po4(i,j,k))
          phyto(n)%felim(i,j,k)  = topaz%f_fed(i,j,k) / (phyto(n)%k_fed + topaz%f_fed(i,j,k))
          !
          !-----------------------------------------------------------------------
          ! Calculate phosphorus and iron deficiency.  Phosphorus deficiency is approximated
          ! from the deviation from the maximum based on the amount of phosphorus needed for
          ! assembly.  Iron deficiency is approximated based on a sigmoidal half-saturation
          ! based on the observed relationship of Sunda and Huntsman (1997).
          !-----------------------------------------------------------------------
          !
          phyto(n)%def_p(i,j,k) = phyto(n)%q_p_2_n(i,j,k) / p_2_n_max(n)
          phyto(n)%def_fe(i,j,k) = phyto(n)%q_fe_2_n(i,j,k) * phyto(n)%q_fe_2_n(i,j,k) /            &
             (phyto(n)%k_fe_2_n * phyto(n)%k_fe_2_n +                                               &
             phyto(n)%q_fe_2_n(i,j,k) * phyto(n)%q_fe_2_n(i,j,k))
       enddo !} n
    enddo;  enddo ;  enddo 

    !
    ! Calculate Leibig nutrient terms to for static and dynamic stoichiometry
    !
    if (topaz%p_2_n_static) then  !{ p:n is static
       if (topaz%fe_2_n_static) then  !{ p:n and fe:n are static
          do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
             n=DIAZO
             phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%po4lim(i,j,k),phyto(n)%felim(i,j,k))
             do n= 2, NUM_PHYTO   !{
                phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),   &
                   phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))
             enddo !} n
          enddo;  enddo ;  enddo !} i,j,k
       else  !{ p:n is static but fe:n is dynamic
          do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
             n=DIAZO
             phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
             do n= 2, NUM_PHYTO   !{
                phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),&
                     phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
             enddo !} n
          enddo;  enddo ;  enddo !} i,j,k
       endif
    else
       if (topaz%fe_2_n_static) then  !{ fe:n is static but p:n is dynamic
          do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
             n=DIAZO
             phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%felim(i,j,k), phyto(n)%def_p(i,j,k))
             do n= 2, NUM_PHYTO   !{
                phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),   &
                     phyto(n)%felim(i,j,k),phyto(n)%def_p(i,j,k))
             enddo !} n
          enddo;  enddo ;  enddo !} i,j,k
       else  !{ pe:n and fe:n are dynamic
          do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
             n=DIAZO
             phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%def_fe(i,j,k), phyto(n)%def_p(i,j,k))
             do n= 2, NUM_PHYTO   !{
                phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),   &
                     phyto(n)%def_fe(i,j,k),phyto(n)%def_p(i,j,k))
             enddo !} n
          enddo;  enddo ;  enddo !} i,j,k
       endif
    endif

    !
    !-----------------------------------------------------------------------
    ! Calculate general ancillary terms
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton in the actively mixed layer are assumed to be photoadapted to
    ! the mean light level in that layer as defined in the KPP routine plus one
    ! more vertical box to account for mixing directly below the boundary layer
    !-----------------------------------------------------------------------
    !

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
             tmp_Irrad = tmp_Irrad + max(0.0,tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k) * 0.5))
             !   Change tmp_irr_band from being the value atop layer k to the value
             ! at the bottom of layer k.
             tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k))
          enddo !}
          topaz%irr_inst(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
          topaz%irr_mix(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
          if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then !{
             kblt = kblt+1
             tmp_irrad_ML = tmp_irrad_ML + topaz%irr_mix(i,j,k) * dzt(i,j,k)
             tmp_hblt = tmp_hblt + dzt(i,j,k)
          endif !}
       enddo !} k-loop
       topaz%irr_mix(i,j,1:kblt) = tmp_irrad_ML / max(1.0e-6,tmp_hblt)

    enddo;  enddo !} i,j

    deallocate(tmp_irr_band)

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
       !
       !-----------------------------------------------------------------------
       ! Temperature functionality of growth and grazing
       !-----------------------------------------------------------------------
       !
       topaz%expkT(i,j,k) = exp(topaz%kappa_eppley * Temp(i,j,k))
       !
       !-----------------------------------------------------------------------
       ! Phytoplankton photoadaptation timescale
       !-----------------------------------------------------------------------
       !
       topaz%f_irr_mem(i,j,k) = (topaz%f_irr_mem(i,j,k) + (topaz%irr_mix(i,j,k) -                   &
          topaz%f_irr_mem(i,j,k)) * min(1.0,topaz%gamma_irr_mem * dt)) * grid_tmask(i,j,k)
    enddo; enddo ; enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    ! This model predicts the Chl:N ratio at each time-step as an equilibrated
    ! phytoplankton response to the combined pressures of light, major nutrient
    ! and iron limitation.  Phytoplankton uptake is generally modelled after
    ! Geider et al. (1997) of steady state N and CO2 uptake, but also includes
    ! the following important modifications:
    !
    ! 1) The temperature effect of Eppley (1972) is used instead
    !    of that in Geider et al (1997) for both simplicity and
    !    to incorporate combined effects on uptake, incorporation
    !    into organic matter and photorespiration.  Values of PCmax
    !    are normalized to 0C rather than 20C in Geider et al. (1997)
    ! 2) The Fe:N ratio is allowed to modulate the Chl:N ratio to be
    !    consistent with Sunda and Huntsman (1997) through the "def_fe"
    !    factor - the phytoplankton Fe:N ratio normalized to a saturated
    !    value (k_fe_2_n) necessary to synthesize chlorophyll.
    !    The def_fe calculation:
    !
    !    def_fe  = (Fe:N)**2 / [(Fe:N)irr**2 + (Fe:N)**2]
    !
    ! 3) Values of the maximum Chl:C ratio (thetamax) are increased and
    !    values of alpha decreased to account for the additional iron
    !    term in the theta equation.
    ! 4) A thetamin_nolim value is also incorporated to set a minimum level of
    !    chlorophyll per carbon at high light, high nutrients
    ! 4) A thetamin value is also incorporated to set an absolute minimum level of
    !    chlorophyll per carbon
    !
    ! While major nutrient limitation is handled through Michaelis Menten
    ! limitation of the phytoplankton specific growth prefactor (P_C_m), iron
    ! limitation is handled indirectly through modulatation of the Chl:N ratio. 
    ! This allows a compensatory relationship between irradiance and iron
    ! availability on phytoplankton specific growth, i.e. if they have a lot of
    ! light, they don't need a lot of iron and vice versa. Def_fe is assumed to be
    ! a quadratic function of the Fe:N ratio nomalized to vary between 0 and 1. 
    ! This relationship is a simple/crude representation of the complex
    ! physiological requirements and functionality of iron by separating
    ! phytoplankton iron into three components:
    !     1) a "basal" requirement of iron for phytoplankton respiration and protein
    !            synthesis (e.g. the electron transport chain)
    !     2) Chlorophyll synthesis for photosynthesis
    !     3) Luxury uptake
    ! While somewhat mathematically ad-hoc, this representation is grounded in the
    ! observed relationship between Chl:C, Fe:C, dissolved Fe and phytoplankton
    ! specific growth rates of Sunda and Huntsman (1997) as well our general
    ! understanding of the role of iron in phytoplankton physiology (e.g Geider and
    ! La Rocha; 1994; Photosynthesis Res., 39, 275-301).  P_C_m is calculated as a Liebig 
    ! nutrient- or stoichiometry- limited rate modulated by the Eppley-type Temperature
    ! function.
    !
    !   P_C_m = P_C_max * e**(k*T) * Liebig_lim
    !
    ! Chl:C calculation after Geider et al (1997) but with addition of a minimum Chl term
    ! that is also modulated by nutrient limitation after optimal allocation work:
    !
    !   Chl:C = (Chl:C_max-Chl:C_min) / (1+ (Chl:C_max-Chl:C_min) * a * I / 2 / P_C_m)
    !      + Chl:C_min * Liebig_lim
    !
    ! Growth rate calculation after Geider et al (1997):
    !
    ! mu = P_C_m / (1 + zeta) * (1 - e**( - alpha * I * Chl:C / P_C_m))
    !
    ! Total Chlorophyll is calculated for use in the short-wave absorption module,
    ! ocean_shortwave_pen.F90.  Note: For this chlorophyll concentration to have
    ! an active effect on irradiance and shortwave absorption, the variable
    ! read_chl=.false. must be set in ocean_shortwave_pen_nml.  Otherwise,
    ! that module will use the data override values.
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
       topaz%f_chl(i,j,k) = 0.0
       do n = 1, NUM_PHYTO   !{
          P_C_m = phyto(n)%liebig_lim(i,j,k) * phyto(n)%P_C_max * topaz%expkT(i,j,k) + epsln
          thetamin_term = max(0.0, topaz%thetamin_nolim - topaz%thetamin) *                         &
             phyto(n)%liebig_lim(i,j,k) + topaz%thetamin
          phyto(n)%theta(i,j,k) = (phyto(n)%thetamax - thetamin_term) / (1.0 +                      &
             (phyto(n)%thetamax - thetamin_term) * phyto(n)%alpha * topaz%f_irr_mem(i,j,k) * 0.5 /  &
             P_C_m) + thetamin_term
          topaz%f_chl(i,j,k) = topaz%f_chl(i,j,k) + topaz%c_2_n * 12.0e6 * phyto(n)%theta(i,j,k) *  &
             phyto(n)%f_n(i,j,k)
          phyto(n)%irrlim(i,j,k) = (1.0 - exp(-phyto(n)%alpha * topaz%irr_inst(i,j,k) *             &
             phyto(n)%theta(i,j,k) / P_C_m))
          phyto(n)%mu(i,j,k) = P_C_m / (1.0 + topaz%zeta) * phyto(n)%irrlim(i,j,k)
       enddo !} n
    enddo;  enddo ; enddo !} i,j,k

    !
    !-----------------------------------------------------------------------
    ! Apply production terms:
    !
    ! NO3 and NH4 uptake are calculated as fractions of total N
    ! uptake
    !
    ! Diazotrophs produce organic N from N2
    !
    ! PO4 production is assumed to be stoichiometric to N with
    ! same stoichiometry for Sm and Lg but higher for Di
    !
    ! Large and Small phytoplankton are allowed to always take up as much Iron
    ! as they can after Sunda and Huntman (1997).  Small phytoplankton are forced
    ! to diminish their uptake at saturated levels of the Fe:C ratio in small
    ! phytoplankton (to mimic their general lack of luxury storage capacity).
    !
    ! Si uptake is make to be consistent with the Si:N ratio
    ! synthesis of Martin-Jezequel et al (2000) and the Droop
    ! quota argument of Mongin et al. (submitted)
    !
    ! CaCO3 formation is set to go directly to detritus as a constant
    ! fraction of Sm production after Moore et al (2002)
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       !
       ! Nitrogen fixation is assumed to be inhibited by high oxygen after Stewart and
       ! Pearson (1970; Effects of aerobic and anaerobic conditions on growth and
       ! metabolism of blue-green algae. Proc. Soc. Lond. B., 175, 293-311).  Nitrogen
       ! fixation is also assumed to be limited by nitrate after Holl and Montoya
       ! (2005;  Interactions between fixed nitrogen uptake and nitrogen fixation in
       ! continuous cultures of the marine diazotroph trichodesmium cyanobacteria); J.
       ! Phycol, 41, 1178-1183).
       !
       n=DIAZO
       jprod_di_tot2nterm=phyto(n)%mu(i,j,k) * phyto(n)%f_n(i,j,k) * (1.0 -                         &
          topaz%f_o2(i,j,k)**topaz%o2_inhib_Di_pow / (topaz%f_o2(i,j,k)**topaz%o2_inhib_Di_pow +  &
          topaz%o2_inhib_Di_sat**topaz%o2_inhib_Di_pow)) / (topaz%f_no3(i,j,k) +                  &
          topaz%f_nh4(i,j,k) + topaz%k_n_inhib_Di)
       phyto(n)%jprod_n2(i,j,k) = jprod_di_tot2nterm * topaz%k_n_inhib_Di      
       phyto(n)%jprod_nh4(i,j,k) = jprod_di_tot2nterm * topaz%f_nh4(i,j,k)
       phyto(n)%jprod_no3(i,j,k) = jprod_di_tot2nterm * topaz%f_no3(i,j,k)
       phyto(n)%jprod_n(i,j,k) = phyto(n)%jprod_n2(i,j,k) + phyto(n)%jprod_nh4(i,j,k) + phyto(n)%jprod_no3(i,j,k)
       do n = 2, NUM_PHYTO
          phyto(n)%jprod_no3(i,j,k) = phyto(n)%mu(i,j,k) * phyto(n)%f_n(i,j,k) *                    &
             phyto(n)%no3lim(i,j,k) / (phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k) + epsln)
          phyto(n)%jprod_nh4(i,j,k) = phyto(n)%mu(i,j,k) * phyto(n)%f_n(i,j,k) *                    &
             phyto(n)%nh4lim(i,j,k) / (phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k) + epsln)
       phyto(n)%jprod_n(i,j,k) = phyto(n)%jprod_nh4(i,j,k) + phyto(n)%jprod_no3(i,j,k)
       enddo !} n
    enddo; enddo ; enddo !} i,j,k

    if (topaz%p_2_n_static) then  !{
       do k = 1, nk  ;    do j = jsc, jec ;      do i = isc, iec   !{
          n=DIAZO
          phyto(n)%jprod_po4(i,j,k) = (phyto(n)%jprod_n2(i,j,k) + phyto(n)%jprod_nh4(i,j,k) +       &
             phyto(n)%jprod_no3(i,j,k)) * phyto(n)%p_2_n_static
          do n = 2, NUM_PHYTO
             phyto(n)%jprod_po4(i,j,k) = (phyto(n)%jprod_no3(i,j,k) + phyto(n)%jprod_nh4(i,j,k)) *  &
                phyto(n)%p_2_n_static
          enddo !} n
       enddo; enddo ; enddo !} i,j,k
    else
       do k = 1, nk  ;    do j = jsc, jec ;      do i = isc, iec   !{
          !
          ! Phosphate uptake patterned after the optimal allocation theory of Klausmeier et al.
          ! 2004 in which N:P is determined by the allocation towards photosynthetic machinery,
          ! nutrient uptake, assemby, and other.  Allocation for "other" is held fixed at 0.2 
          ! after Klausmeier et al. (2004).  Allocation towards photosynthetic machinery
          ! is taken from the instantaneous C:Chl ratio which is renormalized to the ratio of
          ! chlorophyll to chloroplasts in the cell, allocation towards uptake is taken from
          ! the degree of nutrient limitation normalized to the remaining space, and assembly is
          ! then taken from the difference, i.e.:
          !    R_other=0.2
          !    R_photo=1/1.87/0.07*Chl2C
          !    R_uptake=(1-min(plim,nlim))*(1-R_other-R_photo)
          !    R_assembly=1-R_other-R_photo-R_uptake
          !    p2n_opt=p2n_assembly*R_assembly+p2n_photo*R_photo+p2n_other*R_other
          !
          ! For efficiency, this is all condensed into a single calculation.  Phytoplankton are
          ! then allowed to accelerate their approach towards this optimum by taking up more or
          ! less phosphate than the optimum N:P linearly with respect to their relative N:P offset.
          !
          ! In order to simulate the role of storage/vacuoles, Eukaryotes (Lg) are allowed to maximize
          ! r_assem while Prokaryotes (Sm, Di) are forced to maximize r_other
          !
          ! Here, the maximum assembly rate is assumed to be equal to the maximum photosynthetic
          ! rate multiplied by 1/(1 - r_other_min).  However, there is some indication that the
          ! maximum assembly rate should be double the maximum photosynthetic rate
          ! based on comparisons of the envelope for maximum growth rates of bacteria of
          ! ~0.72e-5 * exp(0.063*T) [s-1] for whole organism phytoplankton after Eppley 
          ! (1972, Fish. Bull., 70, 1063-1084) and ~1.6e-5 * exp(0.063*T) [s-1] for whole
          ! organism heterotrophic bacteria after Ducklow and Hill (1985, Limnol. 
          ! Oceanogr., 30, 239-259).  Assuming that the r_assem of heterotrophic bacteria is 0.8
          ! at the maximum growth rate (e.g. r_other=0.2 after Klausmeier et al., 2004) and 
          ! r_assem of phytoplankton is 0.64 at the maximum growth rate (Klausmeier et al., 2004)
          ! gives P_C_max_assem = 1.6e-5 / 0.72e-5 * 0.64 / 0.8 * P_C_max.  The conversion is highly
          ! uncertain however, and leave room for future tuning.

          n=DIAZO
          r_photo = phyto(n)%plast_2_chl * phyto(n)%theta(i,j,k)
          r_uptake = (1.0 - min(phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))) *                   &
             phyto(n)%r_uptake_max
          !
          ! diaz phytoplankton are assumed to take up only as much PO4 for assembly as they need
          ! based on their overall nutrient limited growth rate, thereby assuming that they cannot
          ! store extra PO4.
          !
          r_assem = min(1.0 - phyto(n)%r_nfix - phyto(n)%r_other_min - r_photo - r_uptake,          &
             (1.0 - phyto(n)%r_other_min - phyto(n)%r_nfix) * phyto(n)%liebig_lim(i,j,k))
          phyto(n)%q_p_2_n_opt(i,j,k) = phyto(n)%p_2_n_assem * r_assem + topaz%p_2_n_uptake *       &
             r_uptake + topaz%p_2_n_photo * r_photo + topaz%p_2_n_RKR * (1.0 - phyto(n)%r_nfix -    &
             r_photo - r_uptake - r_assem)
          phyto(n)%jprod_po4(i,j,k) = min(phyto(n)%po4lim(i,j,k) * phyto(n)%P_C_max *               &
             topaz%expkT(i,j,k) * phyto(n)%f_n(i,j,k), phyto(n)%jprod_n2(i,j,k) +                   &
             phyto(n)%jprod_no3(i,j,k) + phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%q_p_2_n_opt(i,j,k)
          n=LARGE
          r_photo = phyto(n)%plast_2_chl * phyto(n)%theta(i,j,k)
          r_uptake = (1.0 - min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),                    &
               phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))) * phyto(n)%r_uptake_max
          !
          ! large phytoplankton are assumed to take up only as much PO4 for assembly as they need
          ! based on their overall PO4 nutrient limited growth rate, thereby assuming that they can
          ! store extra PO4.
          !
          r_assem = min(1.0 - phyto(n)%r_other_min - r_photo - r_uptake,                            &
             (1.0 - phyto(n)%r_other_min) * phyto(n)%po4lim(i,j,k))
          phyto(n)%q_p_2_n_opt(i,j,k) = phyto(n)%p_2_n_assem * r_assem +                            &
             topaz%p_2_n_uptake * r_uptake + topaz%p_2_n_photo * r_photo + topaz%p_2_n_RKR *        &
             (1.0 -r_photo - r_uptake - r_assem)
          phyto(n)%jprod_po4(i,j,k) = min(phyto(n)%po4lim(i,j,k) * phyto(n)%P_C_max *               &
             topaz%expkT(i,j,k) * phyto(n)%f_n(i,j,k),phyto(n)%jprod_no3(i,j,k) +                   &
             phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%q_p_2_n_opt(i,j,k)
          n=SMALL
          r_photo = phyto(n)%plast_2_chl * phyto(n)%theta(i,j,k)
          r_uptake = (1.0 - min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),                    &
             phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))) * phyto(n)%r_uptake_max
          !
          ! small phytoplankton are assumed to take up only as much PO4 for assembly as they need
          ! based on their overall nutrient limited growth rate, thereby assuming that they cannot
          ! store extra PO4.
          !
          r_assem = min(1.0 - phyto(n)%r_other_min - r_photo - r_uptake,                            &
             (1.0 - phyto(n)%r_other_min) * phyto(n)%liebig_lim(i,j,k))
          phyto(n)%q_p_2_n_opt(i,j,k) = phyto(n)%p_2_n_assem * r_assem +                            &
             topaz%p_2_n_uptake * r_uptake + topaz%p_2_n_photo * r_photo + topaz%p_2_n_RKR *        &
             (1.0 - r_photo - r_uptake - r_assem)
          phyto(n)%jprod_po4(i,j,k) = min(phyto(n)%po4lim(i,j,k) * phyto(n)%P_C_max *               &
             topaz%expkT(i,j,k) * phyto(n)%f_n(i,j,k), phyto(n)%jprod_no3(i,j,k) +                  &
             phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%q_p_2_n_opt(i,j,k)
       enddo; enddo ; enddo !} i,j,k
    endif

    if (topaz%fe_2_n_static) then  !{
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
          phyto(DIAZO)%jprod_fe(i,j,k) = (phyto(DIAZO)%jprod_n2(i,j,k) +                            &
             phyto(DIAZO)%jprod_nh4(i,j,k) + phyto(DIAZO)%jprod_no3(i,j,k)) *                       &
             phyto(DIAZO)%fe_2_n_static
          do n = 2, NUM_PHYTO  !{
             phyto(n)%jprod_Fe(i,j,k) = (phyto(n)%jprod_no3(i,j,k) +                                &
                phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%fe_2_n_static
          enddo
       enddo; enddo ; enddo !} i,j,k
    else
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
          do n = 1, NUM_PHYTO  !{
             phyto(n)%jprod_Fe(i,j,k) = phyto(n)%P_C_max * topaz%expkT(i,j,k) *                     &
                phyto(n)%f_n(i,j,k) * phyto(n)%felim(i,j,k) * (phyto(n)%fe_2_n_max -                &
                phyto(n)%q_fe_2_n(i,j,k)) / phyto(n)%fe_2_n_max * phyto(n)%fe_2_n_norm
          enddo   !} n
       enddo; enddo ; enddo !} i,j,k
    endif

    if (topaz%si_2_n_static) then  !{
       do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec   !{
          topaz%nLg_diatoms(i,j,k) = phyto(LARGE)%f_n(i,j,k) * phyto(LARGE)%silim(i,j,k)
          topaz%q_si_2_n_Lg_diatoms(i,j,k) = phyto(LARGE)%si_2_n_static
          phyto(LARGE)%jprod_sio4(i,j,k) = (phyto(LARGE)%jprod_no3(i,j,k) +                         &
             phyto(LARGE)%jprod_nh4(i,j,k)) * phyto(LARGE)%silim(i,j,k) * phyto(LARGE)%si_2_n_static
       enddo; enddo ; enddo !} i,j,k
    else
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
          !
          ! The fraction of large phytoplankton that is diatoms is set equal to silim.  All
          ! other large phytoplankton are assumed to be non-diatom, e.g. green algae,
          ! dinoflagellates, phaeocystis, etc.  SiO4 uptake is limited by the maximum quota
          ! as for PO4 and Fed.
          !
          topaz%nLg_diatoms(i,j,k) = phyto(LARGE)%f_n(i,j,k) * phyto(LARGE)%silim(i,j,k)
          topaz%q_si_2_n_Lg_diatoms(i,j,k) = min(phyto(LARGE)%si_2_n_max, topaz%f_silg(i,j,k) /     &
             (topaz%nLg_diatoms(i,j,k) + epsln))
          phyto(LARGE)%jprod_sio4(i,j,k) = phyto(LARGE)%P_C_max * topaz%expkT(i,j,k) *              &
             topaz%nLg_diatoms(i,j,k) * phyto(LARGE)%silim(i,j,k) * max(0.0,                        &
             phyto(LARGE)%si_2_n_max - topaz%q_si_2_n_Lg_diatoms(i,j,k))
       enddo; enddo ; enddo !} i,j,k
    endif

    !
    !-----------------------------------------------------------------------
    !
    !     Food Web Processing
    !
    ! Phytoplankton loss:
    ! Sm loss is proportional to Sm**2,
    ! Lg loss is proportional to Lg**(4/3)
    ! after Dunne et al, 2005.
    ! for the purposes of grazing, Di and Lg are both considered part of the
    ! "Large" or traditional ecosystem pool. Prey switching/selectivity is assumed
    ! to factor into determining the overall grazing pressure on these two groups
    ! after Fasham et al. (1990; JMR, ) except that the effective size is used
    ! as a proxy for particle number for encounter rates with zooplankton predators.
    ! For this case, the sum of these species is used for the overall grazing
    ! pressure (nlg_tot=ndi+nlg) and the individual grazing pressure on each 
    ! component is determined by:
    !
    !  Pi / eff_size_i / sum((Pi / eff_size_i)^2)^(1/2)
    !
    ! Additionally two criteria for numerical stability are added:
    !    1) the absolute first order rate constant is
    ! never allowed to be greater than 1/dt for numerical stability.
    !    2) a Michaelis-Menton type of threshold using a half
    ! saturation value of phyto_min is set to prevent phytoplankton
    ! from going extinct at low concentrations.
    !
    ! Nitrification is set to be inhibited by light after Ward et al. (1982)
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       graz_Lg_terms = topaz%lambda0 * topaz%expkT(i,j,k) * ((phyto(DIAZO)%f_n(i,j,k) +             &
          phyto(LARGE)%f_n(i,j,k)) / topaz%P_star)**(1.0/3.0) * (phyto(DIAZO)%f_n(i,j,k) +          &
          phyto(LARGE)%f_n(i,j,k)) / (phyto(DIAZO)%f_n(i,j,k) + phyto(LARGE)%f_n(i,j,k) +           &
          topaz%phyto_min) / sqrt(phyto(DIAZO)%f_n(i,j,k) * phyto(DIAZO)%f_n(i,j,k) +               &
          phyto(LARGE)%f_n(i,j,k) * phyto(LARGE)%f_n(i,j,k) + epsln)
       phyto(DIAZO)%jgraz_n(i,j,k) = min( r_dt , graz_Lg_terms * phyto(DIAZO)%f_n(i,j,k)) *         &
          phyto(DIAZO)%f_n(i,j,k)
       phyto(LARGE)%jgraz_n(i,j,k) = min( r_dt , graz_Lg_terms * phyto(LARGE)%f_n(i,j,k)) *         &
          phyto(LARGE)%f_n(i,j,k)
       phyto(SMALL)%jgraz_n(i,j,k) = min( r_dt , topaz%lambda0 * topaz%expkT(i,j,k) *               &
          phyto(SMALL)%f_n(i,j,k) * phyto(SMALL)%f_n(i,j,k) / (topaz%P_star *                       &
          (phyto(SMALL)%f_n(i,j,k) + topaz%phyto_min))) * phyto(SMALL)%f_n(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       jgraz_n = 0.0
       jgraz_p = 0.0
       do n = 1, NUM_PHYTO  !{
          jgraz_n = jgraz_n + phyto(n)%jgraz_n(i,j,k)
          jgraz_p = jgraz_p + phyto(n)%q_p_2_n(i,j,k) * phyto(n)%jgraz_n(i,j,k)
       enddo   !} n
       topaz%jprod_ndet(i,j,k) = (phyto(SMALL)%fdet0 * (phyto(SMALL)%jgraz_n(i,j,k) +               &
          phyto(DIAZO)%jgraz_n(i,j,k)) + phyto(LARGE)%fdet0 * phyto(LARGE)%jgraz_n(i,j,k)) *        &
            (1.0 - topaz%phi_sdon - topaz%phi_ldon) * exp(topaz%kappa_remin * Temp(i,j,k))
       topaz%frac_det_prod(i,j,k) = topaz%jprod_ndet(i,j,k) / (jgraz_n + epsln)
       p_lim_nhet = min(1.0, jgraz_p / (jgraz_n + epsln) / topaz%p_2_n_RKR)
       topaz%jprod_pdet(i,j,k) = topaz%frac_det_prod(i,j,k) * jgraz_p * p_lim_nhet
       topaz%jsdon(i,j,k) = topaz%phi_sdon * jgraz_n
       topaz%jsdop(i,j,k) = topaz%phi_sdop * jgraz_p
       topaz%jldon(i,j,k) = topaz%phi_ldon * jgraz_n * p_lim_nhet
       topaz%jprod_nhet(i,j,k) = (jgraz_n - topaz%jprod_ndet(i,j,k) - topaz%jsdon(i,j,k) -          &
          topaz%jldon(i,j,k)) * topaz%phi_nhet * p_lim_nhet
       topaz%jnh4_graz(i,j,k) = jgraz_n - topaz%jprod_ndet(i,j,k) - topaz%jsdon(i,j,k) -            &
          topaz%jldon(i,j,k) - topaz%jprod_nhet(i,j,k)
       topaz%jpo4_graz(i,j,k) = jgraz_p - topaz%jprod_pdet(i,j,k) - topaz%jsdop(i,j,k) -            &
          (topaz%jldon(i,j,k) + topaz%jprod_nhet(i,j,k)) * topaz%p_2_n_RKR
       topaz%jnhet(i,j,k) = topaz%gamma_nhet *topaz%expkT(i,j,k) * topaz%f_nhet(i,j,k) 
       topaz%jnitrif(i,j,k) = topaz%gamma_nitrif * topaz%expkT(i,j,k) * topaz%f_nh4(i,j,k) *        &
          phyto(SMALL)%nh4lim(i,j,k) * (1.0 - topaz%f_irr_mem(i,j,k) /                              &
          (topaz%irr_inhibit + topaz%f_irr_mem(i,j,k)))
       !
       ! Lithogenic material is assumed to get converted into the sinking particulate phase
       ! through incorporation into mesozooplankton fecal pellets with an efficiency of phi_lith.
       !
       topaz%jprod_lithdet(i,j,k) = (phyto(LARGE)%jgraz_n(i,j,k) / (phyto(LARGE)%f_n(i,j,k) +       &
            epsln) * topaz%phi_lith + topaz%k_lith) * topaz%f_lith(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    do j = jsc, jec ;      do i = isc, iec   !{
       topaz%zt(i,j,1) = dzt(i,j,1)
    enddo; enddo !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       topaz%zt(i,j,k) = topaz%zt(i,j,k-1) + dzt(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       !
       !---------------------------------------------------------------------
       ! CaCO3 solubility taken from the UNESCO 87 recommendation (Mucci for 1-atm solubility with a
       ! fudge factor, times Millero's pressure dependence) gives a result
       ! that splits the difference between the low-end 
       ! solublities (Mucci/Morse) and the high-end solubilities (Ingle/
       ! Plath).  Following unesco, this calculation has a fudge factor built
       ! in that will give 4.5e-7 at 2c, 35 psu, 1-atm.  This is because 
       ! the UNESCO group couldn't decide between the low solubility 
       ! of Mucci (4.3 e-7) and the higher solubility of Plath (4.7e-7) and so 
       ! they just split the difference and called it good.
       ! The salinity normalization for Ca++ breaks down below S = 5.
       !---------------------------------------------------------------------
       !
       TK = Temp(i,j,k) + 273.15
       PRESS = 0.1016 * topaz%zt(i,j,k) + 1.013
       PKSPA = 171.945 + 0.077993 * TK - 2903.293 / TK - 71.595 * log10(TK) - (-0.068393 + 1.7276e-3 * &
          TK + 88.135 / TK) * sqrt(max(epsln, Salt(i,j,k))) + 0.10018 * max(epsln, Salt(i,j,k)) -      &
          5.9415e-3 * max(epsln, Salt(i,j,k))**(1.5) - 0.02 - (48.76 - 2.8 - 0.5304 * Temp(i,j,k)) *   &
          (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 - 0.3692 * Temp(i,j,k))) * (PRESS - 1.013) *&
          (PRESS - 1.013) / (382.92 * TK)
       topaz%co3_sol_arag(i,j,k) = 10**(-PKSPA) / (2.937d-4 * max(5.0, Salt(i,j,k)))
       topaz%omega_arag(i,j,k) = topaz%f_co3_ion(i,j,k) / topaz%co3_sol_arag(i,j,k)
       PKSPC = 171.9065 + 0.077993 * TK - 2839.319 / TK - 71.595 * log10(TK) - (-0.77712 + 2.8426e-3 * &
          TK + 178.34 / TK) * sqrt(max(epsln, Salt(i,j,k))) + 0.07711 * max(epsln, Salt(i,j,k)) -      &
          4.1249e-3 * max(epsln, Salt(i,j,k))**(1.5) - 0.02 - (48.76 - 0.5304 * Temp(i,j,k)) *         &
          (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 - 0.3692 * Temp(i,j,k))) * (PRESS - 1.013) *&
          (PRESS - 1.013) / (382.92 * TK)
       topaz%co3_sol_calc(i,j,k) = 10**(-PKSPC) / (2.937d-4 * max(5.0, Salt(i,j,k)))
       topaz%omega_calc(i,j,k) = topaz%f_co3_ion(i,j,k) / topaz%co3_sol_calc(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    if (topaz%ca_2_n_static) then  !{
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
          topaz%jprod_cadet_arag(i,j,k) = epsln
          topaz%jprod_cadet_calc(i,j,k) = (phyto(DIAZO)%jgraz_n(i,j,k) + phyto(LARGE)%jgraz_n(i,j,k) + &
             phyto(SMALL)%jgraz_n(i,j,k)) * topaz%ca_2_n_het_static + epsln
       enddo; enddo ; enddo !} i,j,k
    else
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
          ! Aragonite and Calcite CaCO3 production is assumed to be proportional to both calcite supersaturation and
          ! microzooplankton grazing after the calcite formulation of Dunne et al. (in prep)
          topaz%jprod_cadet_arag(i,j,k) = phyto(LARGE)%jgraz_n(i,j,k) * topaz%ca_2_n_arag *         &
             min(topaz%caco3_sat_max, max(0.0,topaz%omega_arag(i,j,k) - 1.0)) + epsln
          topaz%jprod_cadet_calc(i,j,k) = phyto(SMALL)%jgraz_n(i,j,k) * topaz%ca_2_n_calc *         &
             exp(-0.0539 * Temp(i,j,k)) * min(topaz%caco3_sat_max, max(0.0, topaz%omega_calc(i,j,k) - 1.0)) + epsln
       enddo; enddo ; enddo !} i,j,k
    endif
    !
    !-----------------------------------------------------------------------
    ! Iron and Silicon Processing:
    !
    ! Iron is recycled with the same efficiency as nitrogen.
    !
    ! SiO2 dissolution is set to globally dissolve 50% after Nelson et al. (1995) 
    ! through grazing.  The temperature functionality is set to a combination
    ! Michaelis  Menton and Eppley temperature formulation to give roughly the range
    ! of  observations in Nelson et al. (1995), Blain et al. (1999) and Bzrezenski's
    ! work... 
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       jgraz_fe = 0.0
       do n = 1, NUM_PHYTO  !{
          phyto(n)%jgraz_fe(i,j,k) = phyto(n)%jgraz_n(i,j,k) * phyto(n)%q_fe_2_n(i,j,k)
          jgraz_fe = jgraz_fe + phyto(n)%jgraz_fe(i,j,k)
       enddo   !} n
       topaz%jprod_fedet(i,j,k) = topaz%frac_det_prod(i,j,k) * jgraz_fe
       topaz%jfe_graz(i,j,k) = jgraz_fe - topaz%jprod_fedet(i,j,k)
       phyto(LARGE)%jgraz_sio2(i,j,k) = phyto(LARGE)%jgraz_n(i,j,k) * topaz%f_siLg(i,j,k) /         &
          (phyto(LARGE)%f_n(i,j,k) + epsln)
       topaz%jdiss_sio2(i,j,k) = phyto(LARGE)%jgraz_sio2(i,j,k) *                                   &
          exp(-topaz%q_si_2_n_Lg_diatoms(i,j,k) / (topaz%q_si_2_n_diss * topaz%expkT(i,j,k)))
    enddo; enddo ; enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    !   Ballast Protection Interior Remineralization Scheme and Iron scavenging
    !-----------------------------------------------------------------------
    !
    !
    do k=1,nk ; do j=jsc,jec ; do i=isc,iec  !{
       if (k .le. grid_kmt(i,j)) then !{
       topaz%jcadet_arag(i,j,k) = topaz%gamma_cadet_arag * max(0.0, 1.0 - topaz%omega_arag(i,j,k)) * topaz%f_cadet_arag(i,j,k)
       topaz%jcadet_calc(i,j,k) = topaz%gamma_cadet_calc * max(0.0, 1.0 - topaz%omega_calc(i,j,k)) * topaz%f_cadet_calc(i,j,k)
       topaz%jsidet(i,j,k) = topaz%gamma_sidet * topaz%f_sidet(i,j,k)
       topaz%jno3denit_wc(i,j,k) = 0.0
       !
       !---------------------------------------------------------------------
       ! Remineralization of unprotected organic material and
       ! previously protected particulate organic material
       !---------------------------------------------------------------------
       !
       !---------------------------------------------------------------------
       !   Under oxic conditions
       !---------------------------------------------------------------------
       !
       if (topaz%f_o2(i,j,k) .gt. topaz%o2_min) then  !{
          topaz%jndet(i,j,k) = topaz%gamma_ndet * topaz%f_o2(i,j,k) / (topaz%k_o2 +                 &
             topaz%f_o2(i,j,k)) * max(0.0, topaz%f_ndet(i,j,k) - (topaz%rpcaco3 *                   &
             (topaz%f_cadet_arag(i,j,k)+topaz%f_cadet_calc(i,j,k)) +                                &
             topaz%rplith * topaz%f_lithdet(i,j,k) + topaz%rpsio2 * topaz%f_sidet(i,j,k)))
       else !}{
          !
          !---------------------------------------------------------------------
          !   Under suboxic conditions, including NO3 limitation
          !---------------------------------------------------------------------
          !
          topaz%jndet(i,j,k) = topaz%gamma_ndet * topaz%o2_min / (topaz%k_o2 + topaz%o2_min) *      &
             topaz%f_no3(i,j,k) / (phyto(SMALL)%k_no3 + topaz%f_no3(i,j,k)) *                       &
             max(0.0, topaz%f_ndet(i,j,k) - (topaz%rpcaco3 * (topaz%f_cadet_arag(i,j,k) +           &
             topaz%f_cadet_calc(i,j,k)) + topaz%rplith * topaz%f_lithdet(i,j,k) + topaz%rpsio2 *    &
             topaz%f_sidet(i,j,k)))
          topaz%jno3denit_wc(i,j,k) = topaz%jndet(i,j,k) * topaz%n_2_n_denit
       endif !}
       !
       !---------------------------------------------------------------------
       ! Apply N change to P assuming equal partitioning between protected,
       ! previously protected and unprotected particulate organic material
       !---------------------------------------------------------------------
       !
       topaz%jpdet(i,j,k) = topaz%jndet(i,j,k) / (topaz%f_ndet(i,j,k) + epsln) * topaz%f_pdet(i,j,k)
       !
       !---------------------------------------------------------------------
       ! Apply N change to Fe incorporating adsorption and desorption
       !---------------------------------------------------------------------
       !
       !---------------------------------------------------------------------
       ! Calculate free and inorganically associated iron concentration for scavenging
       ! based on chemical equilibrium with ligand bound iron using the following
       ! three equations:
       ! kfe_eq_lig = [Feprime] * [L] / [FeL]
       ! [Fed] = [Feprime] + [FeL]
       ! [Ltotal] = [L] + [FeL] = felig_bkg + felig_2_don * (ldon+sdon) 
       !
       ! These equations are solved for Feprime using the quadratic formula:
       ! (-b +- sqrt(b^2 - 4ac))/(2a) where:
       ! a = kfe_eq_lig
       ! b = 1.0 + kfe_eq_lig * (felig_bkg + felig_2_don * (ldon + sdon) - fed)
       ! c = -1
       !---------------------------------------------------------------------
       !
       feprime = 1.0 + topaz%kfe_eq_lig * (topaz%felig_bkg + topaz%felig_2_don *                    &
          (topaz%f_ldon(i,j,k) + topaz%f_sdon(i,j,k)) - topaz%f_fed(i,j,k))
       feprime = (-feprime + sqrt(feprime * feprime + 4.0 * topaz%kfe_eq_lig *                      &
          topaz%f_fed(i,j,k))) / (2.0 * topaz%kfe_eq_lig)
       !
       !---------------------------------------------------------------------
       ! The absolute first order rate constant is never allowed to be greater than
       ! 1/dt for numerical stability.
       !---------------------------------------------------------------------
       !
       topaz%jfe_ads(i,j,k) = min(r_dt, topaz%kfe_org * (topaz%f_ndet(i,j,k) * topaz%mass_2_n +     &
          topaz%kfe_bal * (topaz%f_sidet(i,j,k) * 60.0 + (topaz%f_cadet_arag(i,j,k) +               &
          topaz%f_cadet_calc(i,j,k)) * 100 + topaz%f_lithdet(i,j,k))) * topaz%Rho_0 * topaz%wsink + &
          topaz%kfe_2nd_order * feprime) * feprime
       topaz%jfe_des(i,j,k)=topaz%kfe_des * topaz%f_fedet(i,j,k)
       !
       !---------------------------------------------------------------------
       ! Choose between associating particulate Fe with ballast and organic matter, 
       ! or just with organic matter
       !---------------------------------------------------------------------
       !
       if (topaz%fe_ballast_assoc) then  !{
          topaz%jfedet(i,j,k) = (topaz%jndet(i,j,k) * topaz%mass_2_n + topaz%jsidet(i,j,k) * 60.0 + &
             (topaz%jcadet_arag(i,j,k) + topaz%jcadet_calc(i,j,k)) * 100.0) /                       &
             (topaz%f_ndet(i,j,k) * topaz%mass_2_n + topaz%f_sidet(i,j,k) * 60.0 +                  &
             (topaz%f_cadet_arag(i,j,k) + topaz%f_cadet_calc(i,j,k)) * 100.0 +                      &
             topaz%f_lithdet(i,j,k) + epsln) * topaz%f_fedet(i,j,k)
       else !}{
          topaz%jfedet(i,j,k) = topaz%jndet(i,j,k) / (topaz%f_ndet(i,j,k) + epsln) *                &
             topaz%f_fedet(i,j,k)
       endif !}
       !
       !---------------------------------------------------------------------
       ! Calculate iron and lithogenic loss to sediments
       !---------------------------------------------------------------------
       !
       !---------------------------------------------------------------------
       ! Apply sediment flux to all ocean cells adjacent or corner to land
       !---------------------------------------------------------------------
       !
       topaz%jfe_coast(i,j,k) = topaz%fe_coast * mask_coast(i,j) * grid_tmask(i,j,k) /              &
            sqrt(grid_dat(i,j))
       endif !}
    enddo; enddo ; enddo  !} i,j,k

    !
    !---------------------------------------------------------------------
    ! Account for remineralization/dissolution of sinking flux, and
    ! sediment processed in bottom box
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{
          !
          !---------------------------------------------------------------------
          ! Subtract sedimentary denitrification after Middelburg et al. 1996,
          ! GBC, 10, 661-673 with additional limitation terms to disallow the consumption
          ! of both more ndet than the incoming flux (wsink*ndet) and more no3 that is available in the
          ! lowest 1 m (no3*rho*1m/dt).  The original algorithm was edited to add a maximum flux input
          ! value in order toavoid extrapolatin gthe polynomial into unobserved space which would 
          ! otherwise give low denitrification rates at very high bottom fluxes.
          !---------------------------------------------------------------------
          !
          if (topaz%f_ndet_btf(i,j,1) .gt. 0.0) then !{
          !
          ! Convert previous sediment N inventory into a total flux equivalent assuming a
          ! 10 cm thick sediment layer and  convert flux to umol C cm-2 d-1 and add max value.
             log_btm_flx = log10(min(43.0,topaz%f_ndet_btf(i,j,1) * topaz%c_2_n * sperd * 100.0))
             topaz%fno3denit_sed(i,j) = min(topaz%f_no3(i,j,k) * topaz%Rho_0 * r_dt,                &
                min(topaz%f_ndet_btf(i,j,1) * topaz%n_2_n_denit, 10**(-0.9543 + 0.7662 *            &
                log_btm_flx - 0.235 * log_btm_flx * log_btm_flx) / (topaz%c_2_n * sperd * 100.0) *  &
                topaz%n_2_n_denit * topaz%f_no3(i,j,k) / (phyto(SMALL)%k_no3 + topaz%f_no3(i,j,k))))
             if (topaz%f_o2(i,j,k) .gt. topaz%o2_min) then  !{
                topaz%fnoxic_sed(i,j) = max(0.0, min(topaz%f_o2(i,j,k) * topaz%Rho_0 * r_dt,        &
                  topaz%f_ndet_btf(i,j,1) - topaz%fno3denit_sed(i,j) / topaz%n_2_n_denit))
             else
                topaz%fnoxic_sed(i,j) = 0.0
             endif !}
             topaz%fno3denit_sed(i,j) = topaz%fno3denit_sed(i,j) + min(topaz%f_no3(i,j,k) *         &
               topaz%Rho_0 * r_dt, (topaz%f_ndet_btf(i,j,1) - topaz%fnoxic_sed(i,j) -                &
               topaz%fno3denit_sed(i,j) / topaz%n_2_n_denit) * topaz%n_2_n_denit)
             topaz%fnfeso4red_sed(i,j) = max(0.0, topaz%f_ndet_btf(i,j,1) - topaz%fnoxic_sed(i,j) - &
                topaz%fno3denit_sed(i,j) / topaz%n_2_n_denit)
          else
             topaz%fnfeso4red_sed(i,j) = 0.0
             topaz%fno3denit_sed(i,j) = 0.0
             topaz%fnoxic_sed(i,j) = 0.0
          endif !}
          !
          !---------------------------------------------------------------------
          ! Calculate iron addition from sediments as a function of organic matter
          ! supply
          !---------------------------------------------------------------------
          !
          topaz%ffe_sed(i,j) = topaz%fe_2_n_sed * topaz%f_ndet_btf(i,j,1)
          !
          !---------------------------------------------------------------------
          ! Determine the flux of CaCO3 retained in sediment using the Dunne et al (in prep)
          ! metamodel calibrated to the Hales (2003) steady state model of CaCO3 burial
          ! using bottom water omega and organic-based dissolution and ultimate burial of
          ! sediment CaCO3 assuming a 10 cm mixed layer advecting downward at lithogenic
          ! and CaCO3-based sediment accumulation rate assuming a density of 2.7 g cm-3,
          ! a porosity of 0.7 and molecular weight of 100 to give: 
          ! 2.7e6*(1-0.7)/100 = 8.1e3 mol m-3.  At steady state, burial
          ! efficiency (f = fcased_burial/f_cadet_calc_btf) reduces to:
          !
          !  f = min(f_cadet_calc_btf,0.165 * f_ndet_btf * topaz%c_2_n) / (0.1244 / spery * max(0.0,
          !      1.0 - topaz%omega_calc(i,j,k) + 4.38 * f_ndet_btf * c_2_n * spery)**2.91 *
          !      (f_lithdet * spery + f_cadet_calc_btf * 100.0 * spery)**-2.55 + f_cadet_calc_btf)
          !
          ! which can be used to generate a useful initial condition.
          ! For numerical stability, the redissolution rate is limited to only consume less than
          ! half of the sediment calcite in a time step.
          ! For bulk consistency with the mechanisms of Hales (2003), the ability of organic
          ! flux to instantaneously consume calcite is limited to half the calcite flux. 
          !---------------------------------------------------------------------
          !
          topaz%fcased_redis(i,j) = max(0.0, min(0.5 * topaz%f_cased(i,j,1) * r_dt, min(0.5 *       &
             topaz%f_cadet_calc_btf(i,j,1), 0.165 * topaz%f_ndet_btf(i,j,1) * topaz%c_2_n) +        &
             0.1244 / spery * max(0.0, 1.0 - topaz%omega_calc(i,j,k) +   &
             4.38 * topaz%f_ndet_btf(i,j,1) * topaz%c_2_n * spery)**(2.91) *                        &
             max(1.0, topaz%f_lithdet_btf(i,j,1) * spery + topaz%f_cadet_calc_btf(i,j,1) * 100.0 *  &
             spery)**(-2.55) * topaz%f_cased(i,j,1)))
          topaz%fcased_burial(i,j) = max(0.0, topaz%f_cadet_calc_btf(i,j,1) * topaz%f_cased(i,j,1) /&
             8.1e3)
          topaz%f_cased(i,j,1) = topaz%f_cased(i,j,1) + (topaz%f_cadet_calc_btf(i,j,1) -            &
             topaz%fcased_redis(i,j) - topaz%fcased_burial(i,j)) / topaz%z_sed * dt *               &
             grid_tmask(i,j,k)
          !
          !-----------------------------------------------------------------------
          !     Calculate external bottom fluxes for tracer_vertdiff
          !-----------------------------------------------------------------------
          !
          topaz%b_alk(i,j) = - 2.0 * (topaz%fcased_redis(i,j) + topaz%f_cadet_arag_btf(i,j,1)) -    &
             topaz%f_ndet_btf(i,j,1) - topaz%alk_2_n_denit * topaz%fno3denit_sed(i,j)
          topaz%b_dic(i,j) = - topaz%fcased_redis(i,j) - topaz%f_cadet_arag_btf(i,j,1) -            &
             topaz%f_ndet_btf(i,j,1) * topaz%c_2_n
          topaz%b_fed(i,j) = - topaz%ffe_sed(i,j)
          topaz%b_nh4(i,j) = - topaz%f_ndet_btf(i,j,1)
          topaz%b_no3(i,j) = topaz%fno3denit_sed(i,j)
          topaz%b_o2(i,j)  = topaz%o2_2_nh4 * (topaz%fnoxic_sed(i,j) + topaz%fnfeso4red_sed(i,j))
          topaz%b_po4(i,j) = - topaz%f_pdet_btf(i,j,1)
          topaz%b_sio4(i,j)= - topaz%f_sidet_btf(i,j,1)
       endif !}
    enddo; enddo  !} i, j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       topaz%f_cased(i,j,k) = 0.0
    enddo; enddo ; enddo  !} i,j,k

    call g_tracer_set_values(tracer_list,'alk',  'btf', topaz%b_alk ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic',  'btf', topaz%b_dic ,isd,jsd)
    call g_tracer_set_values(tracer_list,'fed',  'btf', topaz%b_fed ,isd,jsd)
    call g_tracer_set_values(tracer_list,'nh4',  'btf', topaz%b_nh4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'no3',  'btf', topaz%b_no3 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2',   'btf', topaz%b_o2  ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4',  'btf', topaz%b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'sio4', 'btf', topaz%b_sio4,isd,jsd)

    !
    !-----------------------------------------------------------------------
    !
    !     CALCULATE SOURCE/SINK TERMS FOR EACH TRACER
    !
    !-----------------------------------------------------------------------

    !
    !Update the prognostics tracer fields via their pointers.
    !

    call g_tracer_get_pointer(tracer_list,'alk'    ,'field',topaz%p_alk    )
    call g_tracer_get_pointer(tracer_list,'cadet_arag','field',topaz%p_cadet_arag)
    call g_tracer_get_pointer(tracer_list,'cadet_calc','field',topaz%p_cadet_calc)
    call g_tracer_get_pointer(tracer_list,'dic'    ,'field',topaz%p_dic    )
    call g_tracer_get_pointer(tracer_list,'fed'    ,'field',topaz%p_fed    )
    call g_tracer_get_pointer(tracer_list,'fedi'   ,'field',topaz%p_fedi   )
    call g_tracer_get_pointer(tracer_list,'felg'   ,'field',topaz%p_felg   )
    call g_tracer_get_pointer(tracer_list,'fesm'   ,'field',topaz%p_fesm   )
    call g_tracer_get_pointer(tracer_list,'fedet'  ,'field',topaz%p_fedet  )
    call g_tracer_get_pointer(tracer_list,'ldon'   ,'field',topaz%p_ldon   )
    call g_tracer_get_pointer(tracer_list,'lith'   ,'field',topaz%p_lith   )
    call g_tracer_get_pointer(tracer_list,'lithdet','field',topaz%p_lithdet)
    call g_tracer_get_pointer(tracer_list,'ndet'   ,'field',topaz%p_ndet   )
    call g_tracer_get_pointer(tracer_list,'ndi'    ,'field',topaz%p_ndi    )
    call g_tracer_get_pointer(tracer_list,'nlg'    ,'field',topaz%p_nlg    )
    call g_tracer_get_pointer(tracer_list,'nsm'    ,'field',topaz%p_nsm    )
    call g_tracer_get_pointer(tracer_list,'nh4'    ,'field',topaz%p_nh4    )
    call g_tracer_get_pointer(tracer_list,'nhet'   ,'field',topaz%p_nhet   )
    call g_tracer_get_pointer(tracer_list,'no3'    ,'field',topaz%p_no3    )
    call g_tracer_get_pointer(tracer_list,'o2'     ,'field',topaz%p_o2     )
    call g_tracer_get_pointer(tracer_list,'pdet'   ,'field',topaz%p_pdet   )
    call g_tracer_get_pointer(tracer_list,'pdi'    ,'field',topaz%p_pdi    )
    call g_tracer_get_pointer(tracer_list,'plg'    ,'field',topaz%p_plg    )
    call g_tracer_get_pointer(tracer_list,'psm'    ,'field',topaz%p_psm    )
    call g_tracer_get_pointer(tracer_list,'po4'    ,'field',topaz%p_po4    )
    call g_tracer_get_pointer(tracer_list,'sdon'   ,'field',topaz%p_sdon   )
    call g_tracer_get_pointer(tracer_list,'sdop'   ,'field',topaz%p_sdop   )
    call g_tracer_get_pointer(tracer_list,'sidet'  ,'field',topaz%p_sidet  )
    call g_tracer_get_pointer(tracer_list,'silg'   ,'field',topaz%p_silg   )
    call g_tracer_get_pointer(tracer_list,'sio4'   ,'field',topaz%p_sio4   )

    if (topaz%id_no3_in_source .gt. 0)       &
       used = send_data(topaz%id_no3_in_source,         topaz%f_no3,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !
    !-----------------------------------------------------------------------
    !     Phytoplankton Nitrogen and Phosphorus
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Diazotrophic Phytoplankton Nitrogen
       !
       topaz%p_ndi(i,j,k,tau) = topaz%p_ndi(i,j,k,tau) + (phyto(DIAZO)%jprod_n2(i,j,k) +            &
          phyto(DIAZO)%jprod_no3(i,j,k) + phyto(DIAZO)%jprod_nh4(i,j,k) -                           &
          phyto(DIAZO)%jgraz_n(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Nitrogen
       !
       topaz%p_nlg(i,j,k,tau) = topaz%p_nlg(i,j,k,tau) + (phyto(LARGE)%jprod_no3(i,j,k) +           &
          phyto(LARGE)%jprod_nh4(i,j,k) - phyto(LARGE)%jgraz_n(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Nitrogen
       !
       topaz%p_nsm(i,j,k,tau) = topaz%p_nsm(i,j,k,tau) + (phyto(SMALL)%jprod_no3(i,j,k) +           &
          phyto(SMALL)%jprod_nh4(i,j,k) - phyto(SMALL)%jgraz_n(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Diazotrophic Phytoplankton Phosphorus
       !
       topaz%p_pdi(i,j,k,tau) = topaz%p_pdi(i,j,k,tau) + (phyto(DIAZO)%jprod_po4(i,j,k) -           &
          phyto(DIAZO)%jgraz_n(i,j,k) * phyto(DIAZO)%q_p_2_n(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Phosphorus
       !
       topaz%p_plg(i,j,k,tau) = topaz%p_plg(i,j,k,tau) + (phyto(LARGE)%jprod_po4(i,j,k) -           &
          phyto(LARGE)%jgraz_n(i,j,k) * phyto(LARGE)%q_p_2_n(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Phosphorus
       !
       topaz%p_psm(i,j,k,tau) = topaz%p_psm(i,j,k,tau) + (phyto(SMALL)%jprod_po4(i,j,k) -           &
          phyto(SMALL)%jgraz_n(i,j,k) * phyto(SMALL)%q_p_2_n(i,j,k)) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     Phytoplankton Silicon and Iron
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Large Phytoplankton Silicon
       !
       topaz%p_silg(i,j,k,tau) = topaz%p_silg(i,j,k,tau) + (phyto(LARGE)%jprod_sio4(i,j,k) -        &
          phyto(LARGE)%jgraz_sio2(i,j,k)) * dt * grid_tmask(i,j,k)
       ! Diazotrophic Phytoplankton Iron
       !
       topaz%p_fedi(i,j,k,tau) = topaz%p_fedi(i,j,k,tau) + (phyto(DIAZO)%jprod_fe(i,j,k) -          &
          phyto(DIAZO)%jgraz_fe(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Iron
       !
       topaz%p_felg(i,j,k,tau) = topaz%p_felg(i,j,k,tau) + (phyto(LARGE)%jprod_fe(i,j,k) -          &
          phyto(LARGE)%jgraz_fe(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Iron
       !
       topaz%p_fesm(i,j,k,tau) = topaz%p_fesm(i,j,k,tau) + (phyto(SMALL)%jprod_fe(i,j,k) -          &
          phyto(SMALL)%jgraz_fe(i,j,k)) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     Heterotrophic N
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       topaz%p_nhet(i,j,k,tau) = topaz%p_nhet(i,j,k,tau) + (topaz%jprod_nhet(i,j,k) -               &
          topaz%jnhet(i,j,k)) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     NO3
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       topaz%jno3(i,j,k) =  topaz%jnitrif(i,j,k) - phyto(DIAZO)%jprod_no3(i,j,k) -                  &
          phyto(LARGE)%jprod_no3(i,j,k) - phyto(SMALL)%jprod_no3(i,j,k) - topaz%jno3denit_wc(i,j,k)
       topaz%p_no3(i,j,k,tau) = topaz%p_no3(i,j,k,tau) + topaz%jno3(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

    !
    !-----------------------------------------------------------------------
    !     Other nutrients
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! NH4
       !
       topaz%jnh4(i,j,k) = topaz%jnh4_graz(i,j,k) + topaz%jnhet(i,j,k) + topaz%gamma_ldon *         &
          topaz%f_ldon(i,j,k) + topaz%gamma_sdon * topaz%f_sdon(i,j,k) + topaz%jndet(i,j,k)-        &
          phyto(DIAZO)%jprod_nh4(i,j,k) - phyto(LARGE)%jprod_nh4(i,j,k) -                           &
          phyto(SMALL)%jprod_nh4(i,j,k) - topaz%jnitrif(i,j,k)
       topaz%p_nh4(i,j,k,tau) = topaz%p_nh4(i,j,k,tau) + topaz%jnh4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! PO4
       !
       topaz%jpo4(i,j,k) =  topaz%jpo4_graz(i,j,k) + (topaz%jnhet(i,j,k) + topaz%gamma_ldon *       &
          topaz%f_ldon(i,j,k)) * topaz%p_2_n_RKR + topaz%gamma_sdop * topaz%f_sdop(i,j,k) +         &
          topaz%jpdet(i,j,k) - phyto(DIAZO)%jprod_po4(i,j,k) - phyto(LARGE)%jprod_po4(i,j,k) -      &
          phyto(SMALL)%jprod_po4(i,j,k)
       topaz%p_po4(i,j,k,tau) = topaz%p_po4(i,j,k,tau) + topaz%jpo4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! SiO4
       !
       topaz%jsio4(i,j,k) = topaz%jsidet(i,j,k) - phyto(LARGE)%jprod_sio4(i,j,k) +                  &
          topaz%jdiss_sio2(i,j,k)
       topaz%p_sio4(i,j,k,tau) = topaz%p_sio4(i,j,k,tau) + topaz%jsio4(i,j,k) * dt *                &
          grid_tmask(i,j,k)
       !
       ! Fed
       !
       topaz%jfed(i,j,k) = topaz%jfe_graz(i,j,k) - phyto(DIAZO)%jprod_Fe(i,j,k) -                   &
          phyto(LARGE)%jprod_Fe(i,j,k) - phyto(SMALL)%jprod_Fe(i,j,k) - topaz%jfe_ads(i,j,k) +      &
          topaz%jfe_des(i,j,k) + topaz%jfedet(i,j,k) + topaz%jfe_coast(i,j,k)
       topaz%p_fed(i,j,k,tau) = topaz%p_fed(i,j,k,tau) + topaz%jfed(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     Detrital Components
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Cadet_Arag
       !
       topaz%p_cadet_arag(i,j,k,tau) = topaz%p_cadet_arag(i,j,k,tau) +                              &
          (topaz%jprod_cadet_arag(i,j,k) - topaz%jcadet_arag(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Cadet_Calc
       !
       topaz%p_cadet_calc(i,j,k,tau) = topaz%p_cadet_calc(i,j,k,tau) +                              &
          (topaz%jprod_cadet_calc(i,j,k) - topaz%jcadet_calc(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Fedet
       !
       topaz%p_fedet(i,j,k,tau) = topaz%p_fedet(i,j,k,tau) + (topaz%jprod_fedet(i,j,k) +            &
          topaz%jfe_ads(i,j,k) - topaz%jfe_des(i,j,k) - topaz%jfedet(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Lithdet
       !
       topaz%p_lithdet(i,j,k,tau) = topaz%p_lithdet(i,j,k,tau) + topaz%jprod_lithdet(i,j,k) * dt *  &
          grid_tmask(i,j,k)
       !
       ! Ndet
       !
       topaz%p_ndet(i,j,k,tau) = topaz%p_ndet(i,j,k,tau) + (topaz%jprod_ndet(i,j,k) -               &
          topaz%jndet(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Pdet
       !
       topaz%p_pdet(i,j,k,tau) = topaz%p_pdet(i,j,k,tau) + (topaz%jprod_pdet(i,j,k) -               &
          topaz%jpdet(i,j,k)) * dt * grid_tmask(i,j,k)
       !
       ! Sidet
       !
       topaz%p_sidet(i,j,k,tau) = topaz%p_sidet(i,j,k,tau) + (phyto(LARGE)%jgraz_sio2(i,j,k) -      &
          topaz%jdiss_sio2(i,j,k) - topaz%jsidet(i,j,k)) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     Dissolved Organic Matter
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Semilabile Dissolved Organic Nitrogen
       !
       topaz%jsdon(i,j,k) = topaz%jsdon(i,j,k) -  topaz%gamma_sdon * topaz%f_sdon(i,j,k)
       topaz%p_sdon(i,j,k,tau) = topaz%p_sdon(i,j,k,tau) +  topaz%jsdon(i,j,k) * dt *               &
          grid_tmask(i,j,k)
       !
       ! Semilabile Dissolved Organic Phosphorus
       !
       topaz%jsdop(i,j,k) = topaz%jsdop(i,j,k) - topaz%gamma_sdop * topaz%f_sdop(i,j,k)
       topaz%p_sdop(i,j,k,tau) = topaz%p_sdop(i,j,k,tau) + topaz%jsdop(i,j,k) * dt *                &
          grid_tmask(i,j,k)
       !
       ! Labile Dissolved Organic Nitrogen
       !
       topaz%jldon(i,j,k) = topaz%jldon(i,j,k) - topaz%gamma_ldon * topaz%f_ldon(i,j,k)
       topaz%p_ldon(i,j,k,tau) = topaz%p_ldon(i,j,k,tau) + topaz%jldon(i,j,k) * dt *                &
          grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     O2
    !
    ! O2 production from nitrate, ammonia and nitrogen fixation and
    ! O2 consumption from production of NH4 from non-sinking particles,
    ! sinking particles and DOM and O2 consumption from nitrification
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j =jsc, jec ; do i = isc, iec  !{
       topaz%jo2(i,j,k) = (topaz%o2_2_no3 * (phyto(DIAZO)%jprod_no3(i,j,k) +                        &
          phyto(LARGE)%jprod_no3(i,j,k) + phyto(SMALL)%jprod_no3(i,j,k)) + topaz%o2_2_nh4 *         &
          (phyto(DIAZO)%jprod_nh4(i,j,k) + phyto(LARGE)%jprod_nh4(i,j,k) +                          &
          phyto(SMALL)%jprod_nh4(i,j,k) + phyto(DIAZO)%jprod_n2(i,j,k))) * grid_tmask(i,j,k)
       !
       !-----------------------------------------------------------------------
       ! If O2 is present
       !-----------------------------------------------------------------------
       !
       if (topaz%f_o2(i,j,k) .gt. topaz%o2_min) then  !{
          topaz%jo2(i,j,k) = topaz%jo2(i,j,k) - topaz%o2_2_nh4 * (topaz%jnh4_graz(i,j,k) +          &
             topaz%jnhet(i,j,k) + topaz%jndet(i,j,k) + topaz%gamma_sdon * topaz%f_sdon(i,j,k) +     &
             topaz%gamma_ldon * topaz%f_ldon(i,j,k)) - topaz%o2_2_nitrif * topaz%jnitrif(i,j,k)
       endif  !}
       topaz%p_o2(i,j,k,tau) = topaz%p_o2(i,j,k,tau) + topaz%jo2(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     The Carbon system
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Alkalinity (see full N equations in parameter definitions section)
       ! CaCO3 cycling affects Alkalinity by a factor of 2.
       ! Nitrate-based Production increases alkalinity by NO3 equivalents.
       ! Ammonia-based Production decreases alkalinity (and reverse for remineralization) by NH4 equivalents.
       ! N2 Production has no effect on alkalinity.
       ! Nitrification decreases alkalinity by 2 NH4 equivalents.
       ! Denitrification decreases alkalinity by 552/472 = 1.169 NO3 equivalents.
       !
       topaz%jalk(i,j,k) = (2.0 * (topaz%jcadet_arag(i,j,k) +         &
          topaz%jcadet_calc(i,j,k) - topaz%jprod_cadet_arag(i,j,k) -                                &
          topaz%jprod_cadet_calc(i,j,k)) + phyto(DIAZO)%jprod_no3(i,j,k) +                          &
          phyto(LARGE)%jprod_no3(i,j,k) + phyto(SMALL)%jprod_no3(i,j,k) + topaz%jnh4_graz(i,j,k) +  &
          topaz%jnhet(i,j,k) + topaz%gamma_ldon * topaz%f_ldon(i,j,k) + topaz%gamma_sdon *          &
          topaz%f_sdon(i,j,k) + topaz%jndet(i,j,k) - phyto(DIAZO)%jprod_nh4(i,j,k) -                &
          phyto(LARGE)%jprod_nh4(i,j,k) - phyto(SMALL)%jprod_nh4(i,j,k) -                           &
          2.0 * topaz%jnitrif(i,j,k) + topaz%alk_2_n_denit * topaz%jno3denit_wc(i,j,k))
       topaz%p_alk(i,j,k,tau) = topaz%p_alk(i,j,k,tau) + topaz%jalk(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! Dissolved Inorganic Carbon
       !
       topaz%jdic(i,j,k) =(topaz%c_2_n * (topaz%jno3(i,j,k) +        &
          topaz%jnh4(i,j,k) + topaz%jno3denit_wc(i,j,k) - phyto(DIAZO)%jprod_n2(i,j,k)) +           &
          topaz%jcadet_arag(i,j,k) + topaz%jcadet_calc(i,j,k) - topaz%jprod_cadet_arag(i,j,k) -     &
          topaz%jprod_cadet_calc(i,j,k))
       topaz%p_dic(i,j,k,tau) = topaz%p_dic(i,j,k,tau) + topaz%jdic(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    !     Lithogenic aluminosilicate particulates
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       topaz%p_lith(i,j,k,tau) = topaz%p_lith(i,j,k,tau) - topaz%jprod_lithdet(i,j,k) * dt *        &
          grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k


    !
    !Set the diagnostics tracer fields.
    !
    call g_tracer_set_values(tracer_list,'cased',  'field',topaz%f_cased    ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'chl',    'field',topaz%f_chl      ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion','field',topaz%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'irr_mem' ,'field',topaz%f_irr_mem ,isd,jsd,ntau=1)

  if (topaz%do_extra_diag_calcs) then !{
    !
    !-----------------------------------------------------------------------
    !       Calculate variables for diagnostics
    !-----------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    ! Calculate CaCO3 saturation and O2 minimum depths
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec  !{
      if (grid_kmt(i,j) .gt. 0) then !{
        topaz%o2min(i,j)=topaz%p_o2(i,j,1,tau)
        topaz%z_o2min(i,j)=topaz%zt(i,j,1)
        topaz%z_sat_arag(i,j)=missing_value1
        topaz%z_sat_calc(i,j)=missing_value1
        topaz%mask_z_sat_arag(i,j) = .FALSE.
        topaz%mask_z_sat_calc(i,j) = .FALSE.
        if (topaz%omega_arag(i,j,1) .le. 1.0) topaz%z_sat_arag(i,j)=0.0
        if (topaz%omega_calc(i,j,1) .le. 1.0) topaz%z_sat_calc(i,j)=0.0
      endif !}
    enddo ; enddo  !} i,j,k
    do j = jsc, jec ; do i = isc, iec  !{
    first = .true.
      do k = 2, nk
         if (k .le. grid_kmt(i,j) .and. first) then !{
           if (topaz%p_o2(i,j,k,tau) .lt. topaz%p_o2(i,j,k-1,tau)) then
             topaz%o2min(i,j)=topaz%p_o2(i,j,k,tau)
             topaz%z_o2min(i,j)=topaz%zt(i,j,k)
           else
             first = .false.
           endif !}
         endif !}
      enddo; 
    enddo ; enddo  !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec  !{
      if (k .le. grid_kmt(i,j)) then !{
        if (topaz%omega_arag(i,j,k) .le. 1.0 .and. topaz%z_sat_arag(i,j) .lt. 0.0) then
          topaz%z_sat_arag(i,j)=topaz%zt(i,j,k)
          topaz%mask_z_sat_arag(i,j) = .TRUE.
        endif
        if (topaz%omega_calc(i,j,k) .le. 1.0 .and. topaz%z_sat_calc(i,j) .lt. 0.0) then
          topaz%z_sat_calc(i,j)=topaz%zt(i,j,k)
          topaz%mask_z_sat_calc(i,j) = .TRUE.
        endif
      endif !}
    enddo; enddo ; enddo  !} i,j,k
    !
    !---------------------------------------------------------------------
    ! Calculate total elemental layer integrals.
    ! Total carbon  = Dissolved Inorganic Carbon + Phytoplankton Carbon
    !   + Dissolved Organic Carbon (including refractory) + Heterotrophic Biomass
    !   + Detrital Orgainc and Inorganic Carbon
    ! For the oceanic carbon budget, a constant 42 uM of dissolved organic
    ! carbon is added to represent the refractory component.
    ! For the oceanic nitrogen budget, a constant 2 uM of dissolved organic
    ! nitrogen is added to represent the refractory component.
    !---------------------------------------------------------------------
    !
    ! 2010/05/02 jgj: changed  4.2e-5 to  4.0e-5 per JPD
    topaz%tot_layer_int_c(:,:,:) = (topaz%p_dic(:,:,:,tau) + 4.0e-5 +                               &
       topaz%p_cadet_arag(:,:,:,tau) + topaz%p_cadet_calc(:,:,:,tau) + topaz%c_2_n *                &
       (topaz%p_ndi(:,:,:,tau) + topaz%p_nlg(:,:,:,tau) + topaz%p_nsm(:,:,:,tau) +                  &
       topaz%p_ldon(:,:,:,tau) + topaz%p_sdon(:,:,:,tau) + topaz%p_nhet(:,:,:,tau) +                &
       topaz%p_ndet(:,:,:,tau))) * rho_dzt(:,:,:)

    topaz%tot_layer_int_fe(:,:,:) = (topaz%p_fed(:,:,:,tau) + topaz%p_fedi(:,:,:,tau) +             &
       topaz%p_felg(:,:,:,tau) + topaz%p_fesm(:,:,:,tau) + topaz%p_fedet(:,:,:,tau)) * rho_dzt(:,:,:)

    topaz%tot_layer_int_n(:,:,:) = (topaz%p_no3(:,:,:,tau) + 2.0e-6 + topaz%p_nh4(:,:,:,tau) +      &
       topaz%p_ndi(:,:,:,tau) + topaz%p_nlg(:,:,:,tau) + topaz%p_nsm(:,:,:,tau) +                   &
       topaz%p_nsm(:,:,:,tau) + topaz%p_ldon(:,:,:,tau) + topaz%p_sdon(:,:,:,tau) +                 &
       topaz%p_nhet(:,:,:,tau) + topaz%p_ndet(:,:,:,tau)) * rho_dzt(:,:,:)

    topaz%tot_layer_int_p(:,:,:) = (topaz%p_po4(:,:,:,tau) + topaz%p_pdi(:,:,:,tau) +               &
       topaz%p_plg(:,:,:,tau) + topaz%p_psm(:,:,:,tau) + topaz%p_sdop(:,:,:,tau) +                  &
       topaz%p_pdet(:,:,:,tau) + topaz%p_2_n_RKR * (topaz%p_nhet(:,:,:,tau) +                       &
       topaz%p_ldon(:,:,:,tau))) * rho_dzt(:,:,:)

    topaz%tot_layer_int_si(:,:,:) = (topaz%p_sio4(:,:,:,tau) + topaz%p_silg(:,:,:,tau) +            &
       topaz%p_sidet(:,:,:,tau)) * rho_dzt(:,:,:)
    !
    !---------------------------------------------------------------------
    ! calculate water column vertical integrals for diagnostics
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec !{
       topaz%wc_vert_int_jfe_coast(i,j) = 0.0
       topaz%wc_vert_int_jno3denit(i,j) = 0.0
       topaz%wc_vert_int_nfix(i,j) = 0.0
       topaz%wc_vert_int_c(i,j) = 0.0
       topaz%wc_vert_int_dic(i,j) = 0.0
    enddo; enddo !} i,j
    do j = jsc, jec ; do i = isc, iec ; do k = 1, nk  !{
          topaz%wc_vert_int_c(i,j) = topaz%wc_vert_int_c(i,j) + topaz%tot_layer_int_c(i,j,k)
          topaz%wc_vert_int_dic(i,j) = topaz%wc_vert_int_dic(i,j) + topaz%p_dic(i,j,k,tau) *        &
             rho_dzt(i,j,k) * grid_tmask(i,j,k)
          topaz%wc_vert_int_jfe_coast(i,j) = topaz%wc_vert_int_jfe_coast(i,j) +                     &
             topaz%jfe_coast(i,j,k) * rho_dzt(i,j,k) * grid_tmask(i,j,k)
          topaz%wc_vert_int_jno3denit(i,j) = topaz%wc_vert_int_jno3denit(i,j) +                     &
             topaz%jno3denit_wc(i,j,k) * rho_dzt(i,j,k) * grid_tmask(i,j,k)
          topaz%wc_vert_int_nfix(i,j) = topaz%wc_vert_int_nfix(i,j) + phyto(DIAZO)%jprod_n2(i,j,k) *&
             rho_dzt(i,j,k) * grid_tmask(i,j,k)
    enddo; enddo; enddo  !} i,j,k

    !
    !---------------------------------------------------------------------
    ! Add external bottom fluxes to specific rates
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec ; do k = 1, nk  !{
       topaz%jalk_plus_btm(i,j,k)  = topaz%jalk(i,j,k) 
       topaz%jdic_plus_btm(i,j,k)  = topaz%jdic(i,j,k) 
       topaz%jfed_plus_btm(i,j,k)  = topaz%jfed(i,j,k) 
       topaz%jnh4_plus_btm(i,j,k)  = topaz%jnh4(i,j,k) 
       topaz%jno3_plus_btm(i,j,k)  = topaz%jno3(i,j,k) 
       topaz%jo2_plus_btm(i,j,k)   = topaz%jo2(i,j,k) 
       topaz%jpo4_plus_btm(i,j,k)  = topaz%jpo4(i,j,k) 
       topaz%jsio4_plus_btm(i,j,k) = topaz%jsio4(i,j,k) 
       topaz%jdin_plus_btm(i,j,k)  = topaz%jno3(i,j,k) + topaz%jnh4(i,j,k)
    enddo; enddo; enddo  !} i,j,k

    do j = jsc, jec ; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{
          topaz%jalk_plus_btm(i,j,k)  = topaz%jalk(i,j,k) +                       &
            (2.0 * (topaz%fcased_redis(i,j) + topaz%f_cadet_arag_btf(i,j,1)) +    &
             topaz%f_ndet_btf(i,j,1) + topaz%alk_2_n_denit * topaz%fno3denit_sed(i,j)) / rho_dzt(i,j,k)
          topaz%jdic_plus_btm(i,j,k)  = topaz%jdic(i,j,k) +                       &
            (topaz%fcased_redis(i,j) + topaz%f_cadet_arag_btf(i,j,1) +            &
             topaz%f_ndet_btf(i,j,1) * topaz%c_2_n) / rho_dzt(i,j,k)
          topaz%jfed_plus_btm(i,j,k)  = topaz%jfed(i,j,k) + topaz%ffe_sed(i,j) / rho_dzt(i,j,k)
          topaz%jnh4_plus_btm(i,j,k)  = topaz%jnh4(i,j,k) + topaz%f_ndet_btf(i,j,1) / rho_dzt(i,j,k)     
          topaz%jno3_plus_btm(i,j,k)  = topaz%jno3(i,j,k) + topaz%fno3denit_sed(i,j) / rho_dzt(i,j,k)
          topaz%jo2_plus_btm(i,j,k)   = topaz%jo2(i,j,k) +                        &
            (topaz%o2_2_nh4 * (topaz%fnoxic_sed(i,j) + topaz%fnfeso4red_sed(i,j))) / rho_dzt(i,j,k)
          topaz%jpo4_plus_btm(i,j,k)  = topaz%jpo4(i,j,k) + topaz%f_pdet_btf(i,j,1) / rho_dzt(i,j,k)     
          topaz%jsio4_plus_btm(i,j,k) = topaz%jsio4(i,j,k) + topaz%f_sidet_btf(i,j,1) / rho_dzt(i,j,k)   
          topaz%jdin_plus_btm(i,j,k)  = topaz%jno3_plus_btm(i,j,k) + topaz%jnh4_plus_btm(i,j,k) 
       endif !}
    enddo; enddo  !} i, j


    allocate(rho_dzt_100(isc:iec,jsc:jec))
    !
    !---------------------------------------------------------------------
    ! calculate upper 100 m vertical integrals
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec !{
       rho_dzt_100(i,j) = rho_dzt(i,j,1)
       topaz%f_alk_int_100(i,j) = topaz%p_alk(i,j,1,tau) * rho_dzt(i,j,1)
       topaz%f_dic_int_100(i,j) = topaz%p_dic(i,j,1,tau) * rho_dzt(i,j,1)
       topaz%f_din_int_100(i,j) = (topaz%p_no3(i,j,1,tau) + topaz%p_nh4(i,j,1,tau)) * rho_dzt(i,j,1)
       topaz%f_fed_int_100(i,j) = topaz%p_fed(i,j,1,tau) * rho_dzt(i,j,1)
       topaz%f_po4_int_100(i,j) = topaz%p_po4(i,j,1,tau) * rho_dzt(i,j,1)
       topaz%f_sio4_int_100(i,j) = topaz%p_sio4(i,j,1,tau) * rho_dzt(i,j,1)
       topaz%fndet_100(i,j) = topaz%f_ndet(i,j,1) * topaz%Rho_0 * topaz%wsink
       topaz%fpdet_100(i,j) = topaz%f_pdet(i,j,1) * topaz%Rho_0 * topaz%wsink 
       topaz%ffedet_100(i,j) = topaz%f_fedet(i,j,1) * topaz%Rho_0 * topaz%wsink 
       topaz%flithdet_100(i,j) = topaz%f_lithdet(i,j,1) * topaz%Rho_0 * topaz%wsink 
       topaz%fsidet_100(i,j) = topaz%f_sidet(i,j,1) * topaz%Rho_0 * topaz%wsink
       topaz%fcadet_arag_100(i,j) = topaz%f_cadet_arag(i,j,1) * topaz%Rho_0 * topaz%wsink
       topaz%fcadet_calc_100(i,j) = topaz%f_cadet_calc(i,j,1) * topaz%Rho_0 * topaz%wsink
       topaz%jalk_100(i,j) = topaz%jalk(i,j,1) * rho_dzt(i,j,1)
       topaz%jdic_100(i,j) = topaz%jdic(i,j,1) * rho_dzt(i,j,1)
       topaz%jdin_100(i,j) = (topaz%jno3(i,j,1) +topaz%jnh4(i,j,1)) * rho_dzt(i,j,1)
       topaz%jfed_100(i,j) = topaz%jfed(i,j,1) * rho_dzt(i,j,1)
       topaz%jnitrif_100(i,j) = topaz%jnitrif(i,j,1) * rho_dzt(i,j,1)
       topaz%jpo4_100(i,j) = topaz%jpo4(i,j,1) * rho_dzt(i,j,1)
       topaz%jprod_cadet_arag_100(i,j) =  topaz%jprod_cadet_arag(i,j,1) * rho_dzt(i,j,1)
       topaz%jprod_cadet_calc_100(i,j) =  topaz%jprod_cadet_calc(i,j,1) * rho_dzt(i,j,1)
       topaz%jprod_fetot_100(i,j) = (phyto(DIAZO)%jprod_fe(i,j,1) + phyto(LARGE)%jprod_fe(i,j,1) +  &
          phyto(SMALL)%jprod_fe(i,j,1)) * rho_dzt(i,j,1)
       topaz%jprod_ntot_100(i,j) = (phyto(DIAZO)%jprod_n(i,j,1) + phyto(LARGE)%jprod_n(i,j,1) +     &
          phyto(SMALL)%jprod_n(i,j,1)) * rho_dzt(i,j,1)
       topaz%jprod_nlg_diatoms_100(i,j) = phyto(LARGE)%jprod_n(i,j,1) * min(1.0,                    &
          topaz%nLg_diatoms(i,j,1) / max(epsln, phyto(LARGE)%f_n(i,j,1))) * rho_dzt(i,j,1)
       topaz%jprod_nlg_nondiatoms_100(i,j) = max(0.0, phyto(LARGE)%jprod_n(i,j,1) * rho_dzt(i,j,1) -&
          topaz%jprod_nlg_diatoms_100(i,j))
       topaz%jprod_no3tot_100(i,j) = (phyto(DIAZO)%jprod_no3(i,j,1) + phyto(LARGE)%jprod_no3(i,j,1) &
          + phyto(SMALL)%jprod_no3(i,j,1)) * rho_dzt(i,j,1)
       topaz%jprod_ptot_100(i,j) = (phyto(DIAZO)%jprod_po4(i,j,1) + phyto(LARGE)%jprod_po4(i,j,1) + &
          phyto(SMALL)%jprod_po4(i,j,1)) * rho_dzt(i,j,1)
       topaz%jsio4_100(i,j) = topaz%jsio4(i,j,1) * rho_dzt(i,j,1)
       do n = 1, NUM_PHYTO  !{
          phyto(n)%jgraz_n_100(i,j) = phyto(n)%jgraz_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n
       phyto(LARGE)%jprod_sio4_100(i,j) = phyto(LARGE)%jprod_sio4(i,j,1) * rho_dzt(i,j,1)
    enddo; enddo !} i,j
    do j = jsc, jec ; do i = isc, iec ; !{
       k_100 = 1
       do k = 2, grid_kmt(i,j)  !{
          if (rho_dzt_100(i,j) .lt. topaz%Rho_0 * 100.0) then
             k_100 = k
             rho_dzt_100(i,j) = rho_dzt_100(i,j) + rho_dzt(i,j,k)
             topaz%f_alk_int_100(i,j) = topaz%f_alk_int_100(i,j) + topaz%p_alk(i,j,k,tau) *         &
                rho_dzt(i,j,k)
             topaz%f_dic_int_100(i,j) = topaz%f_dic_int_100(i,j) + topaz%p_dic(i,j,k,tau) *         &
                rho_dzt(i,j,k)
             topaz%f_din_int_100(i,j) = topaz%f_din_int_100(i,j) + (topaz%p_no3(i,j,k,tau) +        &
                topaz%p_nh4(i,j,k,tau)) * rho_dzt(i,j,k)
             topaz%f_fed_int_100(i,j) = topaz%f_fed_int_100(i,j) + topaz%p_fed(i,j,k,tau) *         &
                rho_dzt(i,j,k)
             topaz%f_po4_int_100(i,j) = topaz%f_po4_int_100(i,j) + topaz%p_po4(i,j,k,tau) *         &
                rho_dzt(i,j,k)
             topaz%f_sio4_int_100(i,j) = topaz%f_sio4_int_100(i,j) + topaz%p_sio4(i,j,k,tau) *      &
                rho_dzt(i,j,k)
             topaz%fndet_100(i,j) = topaz%f_ndet(i,j,k) * topaz%Rho_0 * topaz%wsink
             topaz%fpdet_100(i,j) = topaz%f_pdet(i,j,k) * topaz%Rho_0 * topaz%wsink
             topaz%ffedet_100(i,j) = topaz%f_fedet(i,j,k) * topaz%Rho_0 * topaz%wsink 
             topaz%flithdet_100(i,j) = topaz%f_lithdet(i,j,k) * topaz%Rho_0 * topaz%wsink 
             topaz%fsidet_100(i,j) = topaz%f_sidet(i,j,k)  *topaz%Rho_0 * topaz%wsink
             topaz%fcadet_arag_100(i,j) = topaz%f_cadet_arag(i,j,k) * topaz%Rho_0 * topaz%wsink
             topaz%fcadet_calc_100(i,j) = topaz%f_cadet_calc(i,j,k) * topaz%Rho_0 * topaz%wsink
             topaz%jalk_100(i,j) = topaz%jalk_100(i,j) + topaz%jalk(i,j,k) * rho_dzt(i,j,k)
             topaz%jdic_100(i,j) = topaz%jdic_100(i,j) + topaz%jdic(i,j,k) * rho_dzt(i,j,k)
             topaz%jdin_100(i,j) = topaz%jdin_100(i,j) + (topaz%jno3(i,j,k) +topaz%jnh4(i,j,k)) *   &
                rho_dzt(i,j,k)
             topaz%jfed_100(i,j) = topaz%jfed_100(i,j) + topaz%jfed(i,j,k) * rho_dzt(i,j,k)
             topaz%jnitrif_100(i,j) = topaz%jnitrif_100(i,j) + topaz%jnitrif(i,j,k) * rho_dzt(i,j,k)
             topaz%jpo4_100(i,j) = topaz%jpo4_100(i,j) + topaz%jpo4(i,j,k) * rho_dzt(i,j,k)
             topaz%jprod_cadet_arag_100(i,j) = topaz%jprod_cadet_arag_100(i,j) +                    &
                topaz%jprod_cadet_arag(i,j,k) * rho_dzt(i,j,k)
             topaz%jprod_cadet_calc_100(i,j) = topaz%jprod_cadet_calc_100(i,j) +                    &
                topaz%jprod_cadet_calc(i,j,k) * rho_dzt(i,j,k)
             topaz%jprod_fetot_100(i,j) = topaz%jprod_fetot_100(i,j) +                              &
                (phyto(DIAZO)%jprod_fe(i,j,k) + phyto(LARGE)%jprod_fe(i,j,k) +                      &
                phyto(SMALL)%jprod_fe(i,j,k)) * rho_dzt(i,j,k)
             topaz%jprod_ntot_100(i,j) = topaz%jprod_ntot_100(i,j) + (phyto(DIAZO)%jprod_n(i,j,k) + &
                phyto(LARGE)%jprod_n(i,j,k) + phyto(SMALL)%jprod_n(i,j,k)) * rho_dzt(i,j,k)
             topaz%jprod_nlg_diatoms_100(i,j) = topaz%jprod_nlg_diatoms_100(i,j) +                  &
                phyto(LARGE)%jprod_n(i,j,k) * min(1.0, topaz%nLg_diatoms(i,j,k) /                   &
                max(epsln, phyto(LARGE)%f_n(i,j,k))) * rho_dzt(i,j,k)
             topaz%jprod_nlg_nondiatoms_100(i,j) = topaz%jprod_nlg_nondiatoms_100(i,j) +            &
                max(0.0, phyto(LARGE)%jprod_n(i,j,k) * (1.0 - min(1.0, topaz%nLg_diatoms(i,j,k) /   &
                max(epsln, phyto(LARGE)%f_n(i,j,k))))) * rho_dzt(i,j,k) 
             topaz%jprod_no3tot_100(i,j) = topaz%jprod_no3tot_100(i,j) +                            &
                (phyto(DIAZO)%jprod_no3(i,j,k) + phyto(LARGE)%jprod_no3(i,j,k) +                    &
                phyto(SMALL)%jprod_no3(i,j,k)) * rho_dzt(i,j,k)
             topaz%jprod_ptot_100(i,j) = topaz%jprod_ptot_100(i,j) + (phyto(DIAZO)%jprod_po4(i,j,k) &
                + phyto(LARGE)%jprod_po4(i,j,k) + phyto(SMALL)%jprod_po4(i,j,k)) * rho_dzt(i,j,k)
             topaz%jsio4_100(i,j) = topaz%jsio4_100(i,j) + topaz%jsio4(i,j,k) * rho_dzt(i,j,k)
             do n = 1, NUM_PHYTO  !{
                phyto(n)%jgraz_n_100(i,j) = phyto(n)%jgraz_n_100(i,j) + phyto(n)%jgraz_n(i,j,k) *   &
                   rho_dzt(i,j,k)
                phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n_100(i,j) + phyto(n)%jprod_n(i,j,k) *   &
                   rho_dzt(i,j,k)
             enddo   !} n
             phyto(LARGE)%jprod_sio4_100(i,j) = phyto(LARGE)%jprod_sio4_100(i,j) +                  &
                phyto(LARGE)%jprod_sio4(i,j,k) * rho_dzt(i,j,k)
          endif
       enddo  !} k
       if (k_100 .gt. 1 .and. k_100 .lt. grid_kmt(i,j)) then
! 2013/05/22 CAS/JPD/JGJ remove line incrementing k_100
          drho_dzt = topaz%Rho_0 * 100.0 - rho_dzt_100(i,j)
          topaz%f_alk_int_100(i,j) = topaz%f_alk_int_100(i,j) + topaz%p_alk(i,j,k_100,tau) * drho_dzt
          topaz%f_dic_int_100(i,j) = topaz%f_dic_int_100(i,j) + topaz%p_dic(i,j,k_100,tau) * drho_dzt
          topaz%f_din_int_100(i,j) = topaz%f_din_int_100(i,j) + (topaz%p_no3(i,j,k_100,tau) +       &
             topaz%p_nh4(i,j,k_100,tau)) * drho_dzt
          topaz%f_fed_int_100(i,j) = topaz%f_fed_int_100(i,j) + topaz%p_fed(i,j,k_100,tau) * drho_dzt
          topaz%f_po4_int_100(i,j) = topaz%f_po4_int_100(i,j) + topaz%p_po4(i,j,k_100,tau) * drho_dzt
          topaz%f_sio4_int_100(i,j) = topaz%f_sio4_int_100(i,j) + topaz%p_sio4(i,j,k_100,tau) *     &
             drho_dzt
          topaz%fndet_100(i,j) = topaz%f_ndet(i,j,k_100) * topaz%Rho_0 * topaz%wsink
          topaz%fpdet_100(i,j) = topaz%f_pdet(i,j,k_100) * topaz%Rho_0 * topaz%wsink
          topaz%ffedet_100(i,j) = topaz%f_fedet(i,j,k_100) * topaz%Rho_0 * topaz%wsink 
          topaz%flithdet_100(i,j) = topaz%f_lithdet(i,j,k_100) * topaz%Rho_0 * topaz%wsink 
          topaz%fsidet_100(i,j) = topaz%f_sidet(i,j,k_100) * topaz%Rho_0 * topaz%wsink
          topaz%fcadet_arag_100(i,j) = topaz%f_cadet_arag(i,j,k_100) * topaz%Rho_0 * topaz%wsink
          topaz%fcadet_calc_100(i,j) = topaz%f_cadet_calc(i,j,k_100) * topaz%Rho_0 * topaz%wsink
          topaz%jalk_100(i,j) = topaz%jalk_100(i,j) + topaz%jalk(i,j,k_100) * drho_dzt
          topaz%jdic_100(i,j) = topaz%jdic_100(i,j) + topaz%jdic(i,j,k_100) * drho_dzt
          topaz%jdin_100(i,j) = topaz%jdin_100(i,j) + (topaz%jno3(i,j,k_100) +                      &
             topaz%jnh4(i,j,k_100)) * drho_dzt
          topaz%jfed_100(i,j) = topaz%jfed_100(i,j) + topaz%jfed(i,j,k_100) * drho_dzt
          topaz%jnitrif_100(i,j) = topaz%jnitrif_100(i,j) + topaz%jnitrif(i,j,k_100) * drho_dzt
          topaz%jpo4_100(i,j) = topaz%jpo4_100(i,j) + topaz%jpo4(i,j,k_100) * drho_dzt
          topaz%jprod_cadet_arag_100(i,j) = topaz%jprod_cadet_arag_100(i,j) +                       &
             topaz%jprod_cadet_arag(i,j,k_100) * drho_dzt
          topaz%jprod_cadet_calc_100(i,j) = topaz%jprod_cadet_calc_100(i,j) +                       &
             topaz%jprod_cadet_calc(i,j,k_100) * drho_dzt
          topaz%jprod_fetot_100(i,j) = topaz%jprod_fetot_100(i,j) +                                 &
             (phyto(DIAZO)%jprod_fe(i,j,k_100) + phyto(LARGE)%jprod_fe(i,j,k_100) +                 &
             phyto(SMALL)%jprod_fe(i,j,k_100)) * drho_dzt
          topaz%jprod_ntot_100(i,j) = topaz%jprod_ntot_100(i,j) + (phyto(DIAZO)%jprod_n(i,j,k_100) +&
             phyto(LARGE)%jprod_n(i,j,k_100) + phyto(SMALL)%jprod_n(i,j,k_100)) * drho_dzt
          topaz%jprod_nlg_diatoms_100(i,j) = topaz%jprod_nlg_diatoms_100(i,j) +                     &
             phyto(LARGE)%jprod_n(i,j,k_100) * min(1.0, topaz%nLg_diatoms(i,j,k_100) /              &
             max(epsln,phyto(LARGE)%f_n(i,j,k_100))) * rho_dzt(i,j,1)
          topaz%jprod_nlg_nondiatoms_100(i,j) = topaz%jprod_nlg_nondiatoms_100(i,j) +               &
             max(0.0, phyto(LARGE)%jprod_n(i,j,k_100) * drho_dzt - topaz%jprod_nlg_diatoms_100(i,j))
          topaz%jprod_no3tot_100(i,j) = topaz%jprod_no3tot_100(i,j) +                               &
             (phyto(DIAZO)%jprod_no3(i,j,k_100) + phyto(LARGE)%jprod_no3(i,j,k_100) +               &
             phyto(SMALL)%jprod_no3(i,j,k_100)) * drho_dzt
          topaz%jprod_ptot_100(i,j) = topaz%jprod_ptot_100(i,j) +                                   &
             (phyto(DIAZO)%jprod_po4(i,j,k_100) + phyto(LARGE)%jprod_po4(i,j,k_100) +               &
             phyto(SMALL)%jprod_po4(i,j,k_100)) * drho_dzt
          topaz%jsio4_100(i,j) = topaz%jsio4_100(i,j) + topaz%jsio4(i,j,k_100) * drho_dzt
          do n = 1, NUM_PHYTO  !{
             phyto(n)%jgraz_n_100(i,j) = phyto(n)%jgraz_n_100(i,j) + phyto(n)%jgraz_n(i,j,k_100) *  &
                drho_dzt
             phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n_100(i,j) + phyto(n)%jprod_n(i,j,k_100) *  &
                drho_dzt
          enddo   !} n
          phyto(LARGE)%jprod_sio4_100(i,j) = phyto(LARGE)%jprod_sio4_100(i,j) +                     &
             phyto(LARGE)%jprod_sio4(i,j,k_100) * drho_dzt
       endif
    enddo ; enddo  !} i,j
    deallocate(rho_dzt_100)
    !
    !---------------------------------------------------------------------
    ! Calculate bulk elemental fluxes
    !---------------------------------------------------------------------
    !
    ! Carbon sources from rivers, loss of calcite to sediments and sediment calcite dissolution
    call g_tracer_get_pointer(tracer_list,'alk','runoff_tracer_flux',runoff_flux_alk)
    call g_tracer_get_pointer(tracer_list,'dic','runoff_tracer_flux',runoff_flux_dic)
    call g_tracer_get_pointer(tracer_list,'fed','runoff_tracer_flux',runoff_flux_fed)
    call g_tracer_get_pointer(tracer_list,'fed','drydep',dry_fed)
    call g_tracer_get_pointer(tracer_list,'fed','wetdep',wet_fed)
    call g_tracer_get_pointer(tracer_list,'lith','runoff_tracer_flux',runoff_flux_lith)
    call g_tracer_get_pointer(tracer_list,'no3','runoff_tracer_flux',runoff_flux_no3)
    call g_tracer_get_pointer(tracer_list,'no3','drydep',dry_no3)
    call g_tracer_get_pointer(tracer_list,'no3','wetdep',wet_no3)
    call g_tracer_get_pointer(tracer_list,'nh4','runoff_tracer_flux',runoff_flux_nh4)
    call g_tracer_get_pointer(tracer_list,'nh4','drydep',dry_nh4)
    call g_tracer_get_pointer(tracer_list,'nh4','wetdep',wet_nh4)
    do j = jsc, jec ; do i = isc, iec  !{
      topaz%nongas_source_c(i,j) = runoff_flux_dic(i,j) + topaz%fcased_redis(i,j)
      topaz%nongas_source_fe(i,j) = runoff_flux_fed(i,j) + topaz%wc_vert_int_jfe_coast(i,j) + &
        dry_fed(i,j) + wet_fed(i,j) + topaz%ffe_sed(i,j)
      topaz%nongas_source_n(i,j) = runoff_flux_no3(i,j) + runoff_flux_nh4(i,j) +               &
        dry_no3(i,j) + wet_no3(i,j) + dry_nh4(i,j) + wet_nh4(i,j) + topaz%wc_vert_int_nfix(i,j)
    enddo ; enddo  !} i,j
    !
    !---------------------------------------------------------------------
    ! Save CaCO3 saturation and O2 minimum depths
    !---------------------------------------------------------------------
    !
    if (topaz%id_o2min .gt. 0)               &
       used = send_data(topaz%id_o2min,         topaz%o2min,                    &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_z_o2min .gt. 0)             &
       used = send_data(topaz%id_z_o2min,    topaz%z_o2min,                     &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_z_sat_arag .gt. 0)          &
       used = send_data(topaz%id_z_sat_arag,    topaz%z_sat_arag,               &
       model_time, mask = topaz%mask_z_sat_arag, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_z_sat_calc .gt. 0)          &
       used = send_data(topaz%id_z_sat_calc,    topaz%z_sat_calc,               &
       model_time, mask = topaz%mask_z_sat_calc, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    !
    !---------------------------------------------------------------------
    ! Save water column vertical integrals
    !---------------------------------------------------------------------
    !
    if (topaz%id_wc_vert_int_c .gt. 0)       &
       used = send_data(topaz%id_wc_vert_int_c,    topaz%wc_vert_int_c,         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_wc_vert_int_dic .gt. 0)     &
       used = send_data(topaz%id_wc_vert_int_dic,  topaz%wc_vert_int_dic,       &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_wc_vert_int_jfe_coast .gt. 0)&
       used = send_data(topaz%id_wc_vert_int_jfe_coast, topaz%wc_vert_int_jfe_coast, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_wc_vert_int_jno3denit .gt. 0)&
       used = send_data(topaz%id_wc_vert_int_jno3denit, topaz%wc_vert_int_jno3denit, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_wc_vert_int_nfix .gt. 0)    &
       used = send_data(topaz%id_wc_vert_int_nfix, topaz%wc_vert_int_nfix,      &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    !
    !---------------------------------------------------------------------
    ! Calculate total elemental concentrations and layer integrals.
    !

    topaz%tot_bfe = topaz%p_fedi(:,:,:,tau) + topaz%p_felg(:,:,:,tau)           &
       + topaz%p_fesm(:,:,:,tau) + topaz%p_fedet(:,:,:,tau)
    if (topaz%id_tot_bfe .gt. 0)             &
       used = send_data(topaz%id_tot_bfe,   topaz%tot_bfe,                      &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%tot_bsi = topaz%p_silg(:,:,:,tau) +  topaz%p_sidet(:,:,:,tau)     
    if (topaz%id_tot_bsi .gt. 0)             &
       used = send_data(topaz%id_tot_bsi,   topaz%tot_bsi,                      &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    ! 2010/05/02 jgj: changed  4.2e-5 to  4.0e-5 per JPD
    topaz%tot_doc = 4.0e-5 +                                                     &
       topaz%c_2_n * (topaz%p_ldon(:,:,:,tau) + topaz%p_sdon(:,:,:,tau))
    if (topaz%id_tot_doc .gt. 0)             &
       used = send_data(topaz%id_tot_doc, topaz%tot_doc,                        &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (topaz%id_tot_layer_int_c .gt. 0)     &
       used = send_data(topaz%id_tot_layer_int_c, topaz%tot_layer_int_c,        &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_tot_layer_int_fe .gt. 0)    &
       used = send_data(topaz%id_tot_layer_int_fe,topaz%tot_layer_int_fe,       &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_tot_layer_int_n .gt. 0)     &
       used = send_data(topaz%id_tot_layer_int_n,topaz%tot_layer_int_n,         &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_tot_layer_int_p .gt. 0)     &
       used = send_data(topaz%id_tot_layer_int_p,topaz%tot_layer_int_p,         &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_tot_layer_int_si .gt. 0)    &
       used = send_data(topaz%id_tot_layer_int_si,topaz%tot_layer_int_si,       &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%tot_phytofe =  topaz%p_fedi(:,:,:,tau) +                              &
       topaz%p_felg(:,:,:,tau) + topaz%p_fesm(:,:,:,tau)
    if (topaz%id_tot_phytofe .gt. 0)             &
       used = send_data(topaz%id_tot_phytofe, topaz%tot_phytofe,                &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%tot_phyton = topaz%p_ndi(:,:,:,tau) +                                 &
       topaz%p_nlg(:,:,:,tau) + topaz%p_nsm(:,:,:,tau)
    if (topaz%id_tot_phyton .gt. 0)             &
       used = send_data(topaz%id_tot_phyton, topaz%tot_phyton,                  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%tot_phytop = topaz%p_pdi(:,:,:,tau) +                                 &
       topaz%p_plg(:,:,:,tau) + topaz%p_psm(:,:,:,tau)
    if (topaz%id_tot_phytop .gt. 0)             &
       used = send_data(topaz%id_tot_phytop, topaz%tot_phytop,                  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%tot_pon = topaz%p_ndi(:,:,:,tau) + topaz%p_nlg(:,:,:,tau) +           &
       topaz%p_nsm(:,:,:,tau) + topaz%p_ndet(:,:,:,tau) + topaz%p_nhet(:,:,:,tau)
    if (topaz%id_tot_pon .gt. 0)             &
       used = send_data(topaz%id_tot_pon, topaz%tot_pon,                        & 
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%tot_pop = topaz%p_pdi(:,:,:,tau) + topaz%p_plg(:,:,:,tau) +           &
       topaz%p_psm(:,:,:,tau) + topaz%p_pdet(:,:,:,tau) +                       &
       topaz%p_2_n_RKR * topaz%p_nhet(:,:,:,tau)
    if (topaz%id_tot_pop .gt. 0)             &
       used = send_data(topaz%id_tot_pop, topaz%tot_pop,                        & 
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (topaz%id_sfc_nlg_diatoms .gt. 0)     &
       used = send_data(topaz%id_sfc_nlg_diatoms, topaz%nlg_diatoms(:,:,1),     &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_nlg_nondiatoms .gt. 0)  &
       used = send_data(topaz%id_sfc_nlg_nondiatoms, max(0.0,                   &
       phyto(LARGE)%f_n(:,:,1) - topaz%nLg_diatoms(:,:,1)),                     &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_bfe .gt. 0)         &
       used = send_data(topaz%id_sfc_tot_bfe,   topaz%p_fedi(:,:,1,tau) +       &
       topaz%p_felg(:,:,1,tau) + topaz%p_fesm(:,:,1,tau) + topaz%p_fedet(:,:,1,tau),&
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_bsi .gt. 0)         &
       used = send_data(topaz%id_sfc_tot_bsi,   topaz%p_silg(:,:,1,tau) +       &
       topaz%p_sidet(:,:,1,tau), &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    ! 2010/05/02 jgj: changed  4.2e-5 to  4.0e-5 per JPD
    if (topaz%id_sfc_tot_doc .gt. 0)         &
       used = send_data(topaz%id_sfc_tot_doc,   4.0e-5 +                        &
       topaz%c_2_n * (topaz%p_ldon(:,:,1,tau) + topaz%p_sdon(:,:,1,tau)),       &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_phytofe .gt. 0)     &
       used = send_data(topaz%id_sfc_tot_phytofe, topaz%p_fedi(:,:,1,tau) +     &
       topaz%p_felg(:,:,1,tau) + topaz%p_fesm(:,:,1,tau),                       &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_phyton .gt. 0)      &
       used = send_data(topaz%id_sfc_tot_phyton,       topaz%p_ndi(:,:,1,tau) + &
       topaz%p_nlg(:,:,1,tau) + topaz%p_nsm(:,:,1,tau), &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_phytop .gt. 0)             &
       used = send_data(topaz%id_sfc_tot_phytop,       topaz%p_pdi(:,:,1,tau) + &
       topaz%p_plg(:,:,1,tau) + topaz%p_psm(:,:,1,tau),                         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_pon .gt. 0)         &
       used = send_data(topaz%id_sfc_tot_pon,       topaz%p_ndi(:,:,1,tau) +    &
       topaz%p_nlg(:,:,1,tau) + topaz%p_nsm(:,:,1,tau) +                        &
       topaz%p_ndet(:,:,1,tau) + topaz%p_nhet(:,:,1,tau),                       &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_tot_pop .gt. 0)         &
       used = send_data(topaz%id_sfc_tot_pop,       topaz%p_pdi(:,:,1,tau) +    &
       topaz%p_plg(:,:,1,tau) + topaz%p_psm(:,:,1,tau) +                        &
       topaz%p_pdet(:,:,1,tau) + topaz%p_2_n_RKR * topaz%p_nhet(:,:,1,tau),     &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    !---------------------------------------------------------------------
    ! Save upper 100 m vertical integrals
    !---------------------------------------------------------------------
    !
    if (topaz%id_f_alk_int_100 .gt. 0)       &
       used = send_data(topaz%id_f_alk_int_100, topaz%f_alk_int_100,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_f_dic_int_100 .gt. 0)       &
       used = send_data(topaz%id_f_dic_int_100, topaz%f_dic_int_100,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_f_din_int_100 .gt. 0)       &
       used = send_data(topaz%id_f_din_int_100, topaz%f_din_int_100,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_f_fed_int_100 .gt. 0)       &
       used = send_data(topaz%id_f_fed_int_100, topaz%f_fed_int_100,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_f_po4_int_100 .gt. 0)       &
       used = send_data(topaz%id_f_po4_int_100, topaz%f_po4_int_100,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_f_sio4_int_100 .gt. 0)      &
       used = send_data(topaz%id_f_sio4_int_100, topaz%f_sio4_int_100,          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_fcadet_arag_100 .gt. 0)     &
       used = send_data(topaz%id_fcadet_arag_100, topaz%fcadet_arag_100,        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_fcadet_calc_100 .gt. 0)     &
       used = send_data(topaz%id_fcadet_calc_100, topaz%fcadet_calc_100,        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_ffedet_100 .gt. 0)          &
       used = send_data(topaz%id_ffedet_100,    topaz%ffedet_100,               &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_flithdet_100 .gt. 0)        &
       used = send_data(topaz%id_flithdet_100,    topaz%flithdet_100,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_fndet_100 .gt. 0)           &
       used = send_data(topaz%id_fndet_100,     topaz%fndet_100,                &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_fno3denit_tot .gt. 0)       &
       used = send_data(topaz%id_fno3denit_tot, topaz%fno3denit_sed + topaz%wc_vert_int_jno3denit, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fpdet_100 .gt. 0)           &
       used = send_data(topaz%id_fpdet_100,     topaz%fpdet_100,                &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_fsidet_100 .gt. 0)          &
       used = send_data(topaz%id_fsidet_100,    topaz%fsidet_100,               &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jalk_100 .gt. 0)            &
       used = send_data(topaz%id_jalk_100,      topaz%jalk_100,                 &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jdic_100 .gt. 0)            &
       used = send_data(topaz%id_jdic_100,      topaz%jdic_100,                 &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jdin_100 .gt. 0)            &
       used = send_data(topaz%id_jdin_100,      topaz%jdin_100,                 &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jfed_100 .gt. 0)            &
       used = send_data(topaz%id_jfed_100,      topaz%jfed_100,                 &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jgraz_ntot = phyto(DIAZO)%jgraz_n + phyto(LARGE)%jgraz_n + phyto(SMALL)%jgraz_n
    if (topaz%id_jgraz_ntot .gt. 0)          &
       used = send_data(topaz%id_jgraz_ntot, topaz%jgraz_ntot * rho_dzt,                  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jgraz_ntot_100 .gt. 0)      &
       used = send_data(topaz%id_jgraz_ntot_100,      (phyto(DIAZO)%jgraz_n_100 &
                + phyto(LARGE)%jgraz_n_100 + phyto(SMALL)%jgraz_n_100),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jnitrif_100 .gt. 0)         &
       used = send_data(topaz%id_jnitrif_100,   topaz%jnitrif_100,              &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jpo4_100 .gt. 0)            &
       used = send_data(topaz%id_jpo4_100,      topaz%jpo4_100,                 &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jprod_cadet_arag_100 .gt. 0)&
       used = send_data(topaz%id_jprod_cadet_arag_100, topaz%jprod_cadet_arag_100, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jprod_cadet_calc_100 .gt. 0)&
       used = send_data(topaz%id_jprod_cadet_calc_100, topaz%jprod_cadet_calc_100, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jprod_fetot=phyto(DIAZO)%jprod_fe + phyto(LARGE)%jprod_fe + phyto(SMALL)%jprod_fe
    if (topaz%id_jprod_fetot .gt. 0)         &
       used = send_data(topaz%id_jprod_fetot, topaz%jprod_fetot * rho_dzt,                &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_fetot_100 .gt. 0)     &
       used = send_data(topaz%id_jprod_fetot_100, topaz%jprod_fetot_100,        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jprod_nlg_diatoms=phyto(LARGE)%jprod_n * min(1.0, topaz%nLg_diatoms / max(epsln, phyto(LARGE)%f_n))
    if (topaz%id_jprod_nlg_diatoms .gt. 0)   &
       used = send_data(topaz%id_jprod_nlg_diatoms, topaz%jprod_nlg_diatoms * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_nlg_diatoms_100 .gt. 0)&
       used = send_data(topaz%id_jprod_nlg_diatoms_100, topaz%jprod_nlg_diatoms_100, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jprod_nlg_nondiatoms=max(0.0,phyto(LARGE)%jprod_n * (phyto(LARGE)%f_n - topaz%nLg_diatoms) /          &
       max(epsln, phyto(LARGE)%f_n))
    if (topaz%id_jprod_nlg_nondiatoms .gt. 0)&
       used = send_data(topaz%id_jprod_nlg_nondiatoms, topaz%jprod_nlg_nondiatoms * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_nlg_nondiatoms_100 .gt. 0) &
       used = send_data(topaz%id_jprod_nlg_nondiatoms_100, topaz%jprod_nlg_nondiatoms_100, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jprod_no3tot = phyto(DIAZO)%jprod_no3 + phyto(LARGE)%jprod_no3 + phyto(SMALL)%jprod_no3
    if (topaz%id_jprod_no3tot .gt. 0)        &
       used = send_data(topaz%id_jprod_no3tot,  topaz%jprod_no3tot * rho_dzt,     &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)


    if (topaz%id_jprod_no3tot_100 .gt. 0)    &
       used = send_data(topaz%id_jprod_no3tot_100, topaz%jprod_no3tot_100,      &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jprod_ntot = phyto(DIAZO)%jprod_n + phyto(LARGE)%jprod_n + phyto(SMALL)%jprod_n
    if (topaz%id_jprod_ntot .gt. 0)          &
       used = send_data(topaz%id_jprod_ntot,    topaz%jprod_ntot * rho_dzt,                  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (topaz%id_jprod_ntot_100 .gt. 0)      &
       used = send_data(topaz%id_jprod_ntot_100, topaz%jprod_ntot_100,          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    topaz%jprod_ptot=phyto(DIAZO)%jprod_po4 + phyto(LARGE)%jprod_po4 + phyto(SMALL)%jprod_po4
    if (topaz%id_jprod_ptot .gt. 0)          &
       used = send_data(topaz%id_jprod_ptot, topaz%jprod_ptot * rho_dzt,     &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_ptot_100 .gt. 0)      &
       used = send_data(topaz%id_jprod_ptot_100, topaz%jprod_ptot_100,          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_jsio4_100 .gt. 0)           &
       used = send_data(topaz%id_jsio4_100,     topaz%jsio4_100,                &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_nLg_diatoms .gt. 0)         &
       used = send_data(topaz%id_nLg_diatoms,   topaz%nLg_diatoms,              &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%nLg_nondiatoms = max(0.0, phyto(LARGE)%f_n - topaz%nLg_diatoms)
    if (topaz%id_nLg_nondiatoms .gt. 0)      &
       used = send_data(topaz%id_nLg_nondiatoms, topaz%nLg_nondiatoms,          &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !
    !---------------------------------------------------------------------
    ! Save river, depositon and bulk elemental fluxes
    !---------------------------------------------------------------------
    !
    if (topaz%id_dep_dry_fed .gt. 0)     &
       used = send_data(topaz%id_dep_dry_fed, dry_fed,                          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_dep_dry_nh4 .gt. 0)     &
       used = send_data(topaz%id_dep_dry_nh4, dry_nh4,                          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_dep_dry_no3 .gt. 0)     &
       used = send_data(topaz%id_dep_dry_no3, dry_no3,                          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_dep_wet_fed .gt. 0)     &
       used = send_data(topaz%id_dep_wet_fed, wet_fed,                          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_dep_wet_nh4 .gt. 0)     &
       used = send_data(topaz%id_dep_wet_nh4, wet_nh4,                          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_dep_wet_no3 .gt. 0)     &
       used = send_data(topaz%id_dep_wet_no3, wet_no3,                          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_runoff_flux_alk .gt. 0)     &
       used = send_data(topaz%id_runoff_flux_alk, runoff_flux_alk,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_runoff_flux_dic .gt. 0)     &
       used = send_data(topaz%id_runoff_flux_dic, runoff_flux_dic,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_runoff_flux_fed .gt. 0)     &
       used = send_data(topaz%id_runoff_flux_fed, runoff_flux_fed,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_runoff_flux_lith .gt. 0)     &
       used = send_data(topaz%id_runoff_flux_lith, runoff_flux_lith,         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_runoff_flux_nh4 .gt. 0)     &
       used = send_data(topaz%id_runoff_flux_nh4, runoff_flux_nh4,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_runoff_flux_no3 .gt. 0)     &
       used = send_data(topaz%id_runoff_flux_no3, runoff_flux_no3,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (topaz%id_nongas_source_c .gt. 0)     &
       used = send_data(topaz%id_nongas_source_c, topaz%nongas_source_c,        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_nongas_source_fe .gt. 0)    &
       used = send_data(topaz%id_nongas_source_fe, topaz%nongas_source_fe,      &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_nongas_source_n .gt. 0)     &
       used = send_data(topaz%id_nongas_source_n, topaz%nongas_source_n,        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif !}
!---------------------------------------------------------------------
! Save phyto type fields
!---------------------------------------------------------------------
    do n= 1, NUM_PHYTO
       !Beware that some of these might not be legitimate output variables
       ! e.g. The original code did not have phyto(DIAZO)%jprod_nh4
       ! However the calls to ocean_register_diag should have been done properly.

       phyto(n)%chl = topaz%c_2_n * 12.0e6 * phyto(n)%theta * phyto(n)%f_n
       if (phyto(n)%id_chl .gt. 0)              &
          used = send_data(phyto(n)%id_chl,  phyto(n)%chl,                      &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_sfc_chl .gt. 0)              &
          used = send_data(phyto(n)%id_sfc_chl,    topaz%c_2_n * 12.0e6 *       &
          phyto(n)%theta(:,:,1) * phyto(n)%f_n(:,:,1),                          &
          model_time, rmask = grid_tmask(:,:,1),& 
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_def_fe .gt. 0)           &
          used = send_data(phyto(n)%id_def_fe,     phyto(n)%def_fe,             &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_felim .gt. 0)            &
          used = send_data(phyto(n)%id_felim,      phyto(n)%felim,              &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_irrlim .gt. 0)           &
          used = send_data(phyto(n)%id_irrlim,     phyto(n)%irrlim,             &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jgraz_n .gt. 0)          &
          used = send_data(phyto(n)%id_jgraz_n,    phyto(n)%jgraz_n * rho_dzt,  &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jgraz_n_100 .gt. 0)      &
          used = send_data(phyto(n)%id_jgraz_n_100, phyto(n)%jgraz_n_100,       &
          model_time, rmask = grid_tmask(:,:,1),& 
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_jgraz_fe .gt. 0)         &
          used = send_data(phyto(n)%id_jgraz_fe,   phyto(n)%jgraz_fe * rho_dzt, &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_fe .gt. 0)         &
          used = send_data(phyto(n)%id_jprod_fe,   phyto(n)%jprod_fe * rho_dzt, &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_n .gt. 0)          &
          used = send_data(phyto(n)%id_jprod_n,    phyto(n)%jprod_n * rho_dzt,  &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_n_100 .gt. 0)      &
          used = send_data(phyto(n)%id_jprod_n_100, phyto(n)%jprod_n_100,       &
          model_time, rmask = grid_tmask(:,:,1),& 
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_jprod_nh4 .gt. 0)        &
          used = send_data(phyto(n)%id_jprod_nh4,  phyto(n)%jprod_nh4 * rho_dzt,&
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_no3 .gt. 0)        &
          used = send_data(phyto(n)%id_jprod_no3,  phyto(n)%jprod_no3 * rho_dzt,&
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_po4 .gt. 0)        &
          used = send_data(phyto(n)%id_jprod_po4,  phyto(n)%jprod_po4 * rho_dzt,&
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_liebig_lim .gt. 0)       &
          used = send_data(phyto(n)%id_liebig_lim,phyto(n)%liebig_lim,          &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_mu .gt. 0)               &
          used = send_data(phyto(n)%id_mu,         phyto(n)%mu,                 &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_nh4lim .gt. 0)           &
          used = send_data(phyto(n)%id_nh4lim,     phyto(n)%nh4lim,             &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_no3lim .gt. 0)           &
          used = send_data(phyto(n)%id_no3lim,     phyto(n)%no3lim,             &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_po4lim .gt. 0)           &
          used = send_data(phyto(n)%id_po4lim,     phyto(n)%po4lim,             &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_fe_2_n .gt. 0)         &
          used = send_data(phyto(n)%id_q_fe_2_n,   phyto(n)%q_fe_2_n,           &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_p_2_n .gt. 0)          &
          used = send_data(phyto(n)%id_q_p_2_n,    phyto(n)%q_p_2_n,            &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_p_2_n_opt .gt. 0)      &
          used = send_data(phyto(n)%id_q_p_2_n_opt, phyto(n)%q_p_2_n_opt,       &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_theta .gt. 0)            &
          used = send_data(phyto(n)%id_theta,      phyto(n)%theta,              &
          model_time, rmask = grid_tmask(:,:,:),& 
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo

    if (phyto(DIAZO)%id_jprod_n2 .gt. 0)     &
       used = send_data(phyto(DIAZO)%id_jprod_n2, phyto(DIAZO)%jprod_n2 * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (phyto(LARGE)%id_jgraz_sio2 .gt. 0)   &
       used = send_data(phyto(LARGE)%id_jgraz_sio2, phyto(LARGE)%jgraz_sio2 * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (phyto(LARGE)%id_jprod_sio4 .gt. 0)   &
       used = send_data(phyto(LARGE)%id_jprod_sio4, phyto(LARGE)%jprod_sio4 * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (phyto(LARGE)%id_jprod_sio4_100 .gt. 0) &
       used = send_data(phyto(LARGE)%id_jprod_sio4_100, phyto(LARGE)%jprod_sio4_100, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(LARGE)%id_silim .gt. 0)        &
       used = send_data(phyto(LARGE)%id_silim,  phyto(LARGE)%silim,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
!---------------------------------------------------------------------
! Save 3D topaz type fields
!---------------------------------------------------------------------

    topaz%chl_lg_diatoms = topaz%c_2_n * 12.0e6 * phyto(LARGE)%theta * topaz%nLg_diatoms
    if (topaz%id_chl_lg_diatoms .gt. 0)      &
       used = send_data(topaz%id_chl_lg_diatoms, topaz%chl_lg_diatoms,          &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%chl_lg_nondiatoms = topaz%c_2_n * 12.0e6 *                            &
          phyto(LARGE)%theta * (phyto(LARGE)%f_n - topaz%nLg_diatoms)
    if (topaz%id_chl_lg_nondiatoms .gt. 0)   &
       used = send_data(topaz%id_chl_lg_nondiatoms, topaz%chl_lg_nondiatoms,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (topaz%id_cased_diag .gt. 0)        &
       used = send_data(topaz%id_cased_diag,  topaz%f_cased,                    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_chl_diag .gt. 0)        &
       used = send_data(topaz%id_chl_diag,  topaz%f_chl,                        &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_co3_ion_diag .gt. 0)        &
       used = send_data(topaz%id_co3_ion_diag,  topaz%f_co3_ion,                &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_co3_sol_arag .gt. 0)        &
       used = send_data(topaz%id_co3_sol_arag,  topaz%co3_sol_arag,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_co3_sol_calc .gt. 0)        &
       used = send_data(topaz%id_co3_sol_calc,  topaz%co3_sol_calc,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_fca_sed_sfc .gt. 0)         &
       used = send_data(topaz%id_fca_sed_sfc,   topaz%fcadet_calc_btm - topaz%fcased_redis, &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fcased_int .gt. 0)       &
       used = send_data(topaz%id_fcased_int, topaz%f_cased(:,:,1),               &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fcadet_arag .gt. 0)         &
       used = send_data(topaz%id_fcadet_arag,   topaz%f_cadet_arag(:,:,:) * topaz%Rho_0 * topaz%wsink, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_fcadet_calc .gt. 0)         &
       used = send_data(topaz%id_fcadet_calc,   topaz%f_cadet_calc(:,:,:) * topaz%Rho_0 * topaz%wsink, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_fcased_burial .gt. 0)       &
       used = send_data(topaz%id_fcased_burial, topaz%fcased_burial,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fcased_redis .gt. 0)        &
       used = send_data(topaz%id_fcased_redis,  topaz%fcased_redis,             &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_ffe_sed .gt. 0)             &
       used = send_data(topaz%id_ffe_sed,       topaz%ffe_sed,                  &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_ffedet .gt. 0)              &
       used = send_data(topaz%id_ffedet,        topaz%f_fedet(:,:,:) *          &
       topaz%Rho_0 * topaz%wsink * grid_tmask(:,:,:),                           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_flithdet .gt. 0)            &
       used = send_data(topaz%id_flithdet,      topaz%f_lithdet(:,:,:) *        &
       topaz%Rho_0 * topaz%wsink * grid_tmask(:,:,:),                           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_fndet .gt. 0)               &
       used = send_data(topaz%id_fndet,         topaz%f_ndet(:,:,:) *           &
       topaz%Rho_0 * topaz%wsink * grid_tmask(:,:,:),                           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_fnfeso4red_sed .gt. 0)      &
       used = send_data(topaz%id_fnfeso4red_sed,topaz%fnfeso4red_sed,           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fno3denit_sed .gt. 0)       &
       used = send_data(topaz%id_fno3denit_sed, topaz%fno3denit_sed,            &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fnoxic_sed .gt. 0)          &
       used = send_data(topaz%id_fnoxic_sed,    topaz%fnoxic_sed,               &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_fpdet .gt. 0)               &
       used = send_data(topaz%id_fpdet,         topaz%f_pdet(:,:,:) *           &
       topaz%Rho_0 * topaz%wsink * grid_tmask(:,:,:),                           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_fsidet .gt. 0)              &
       used = send_data(topaz%id_fsidet,        topaz%f_sidet(:,:,:)  *         &
       topaz%Rho_0 * topaz%wsink * grid_tmask(:,:,:),                           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_htotal_diag .gt. 0)        &
       used = send_data(topaz%id_htotal_diag,  topaz%f_htotal,                  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_irr_inst .gt. 0)            &
       used = send_data(topaz%id_irr_inst,      topaz%irr_inst,                 &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_irr_mem_diag .gt. 0)        &
       used = send_data(topaz%id_irr_mem_diag,  topaz%f_irr_mem,                &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_irr_mix .gt. 0)             &
       used = send_data(topaz%id_irr_mix,       topaz%irr_mix,                  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_cadet_arag .gt. 0)    &
       used = send_data(topaz%id_jprod_cadet_arag, topaz%jprod_cadet_arag * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jalk .gt. 0)                &
       used = send_data(topaz%id_jalk,          topaz%jalk * rho_dzt,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jalk_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jalk_plus_btm, topaz%jalk_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jcadet_arag .gt. 0)         &
       used = send_data(topaz%id_jcadet_arag,   topaz%jcadet_arag * rho_dzt,      &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jcadet_calc .gt. 0)         &
       used = send_data(topaz%id_jcadet_calc,   topaz%jcadet_calc * rho_dzt,      &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jdic .gt. 0)                &
       used = send_data(topaz%id_jdic,          topaz%jdic * rho_dzt,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jdic_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jdic_plus_btm, topaz%jdic_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    topaz%jdin=topaz%jno3 +topaz%jnh4
    if (topaz%id_jdin .gt. 0)                &
       used = send_data(topaz%id_jdin, topaz%jdin * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jdin_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jdin_plus_btm, topaz%jdin_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jdiss_sio2 .gt. 0)          &
       used = send_data(topaz%id_jdiss_sio2,    topaz%jdiss_sio2*rho_dzt,         &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfe_ads .gt. 0)             &
       used = send_data(topaz%id_jfe_ads,       topaz%jfe_ads * rho_dzt,          &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfe_des .gt. 0)             &
       used = send_data(topaz%id_jfe_des,       topaz%jfe_des*rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfe_graz .gt. 0)            &
       used = send_data(topaz%id_jfe_graz,      topaz%jfe_graz*rho_dzt,           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfe_coast .gt. 0)           &
       used = send_data(topaz%id_jfe_coast, topaz%jfe_coast*rho_dzt,              &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfed .gt. 0)                &
       used = send_data(topaz%id_jfed,          topaz%jfed * rho_dzt,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfed_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jfed_plus_btm, topaz%jfed_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jfedet .gt. 0)              &
       used = send_data(topaz%id_jfedet,        topaz%jfedet * rho_dzt,           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jldon .gt. 0)               &
       used = send_data(topaz%id_jldon,         topaz%jldon*rho_dzt,              &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jndet .gt. 0)               &
       used = send_data(topaz%id_jndet,         topaz%jndet*rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jnh4 .gt. 0)                &
       used = send_data(topaz%id_jnh4,          topaz%jnh4*rho_dzt,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jnh4_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jnh4_plus_btm, topaz%jnh4_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jnh4_graz .gt. 0)           &
       used = send_data(topaz%id_jnh4_graz,     topaz%jnh4_graz*rho_dzt,        &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jnhet .gt. 0)               &
       used = send_data(topaz%id_jnhet,         topaz%jnhet*rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jnitrif .gt. 0)             &
       used = send_data(topaz%id_jnitrif,       topaz%jnitrif*rho_dzt,          &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jno3 .gt. 0)                &
       used = send_data(topaz%id_jno3,          topaz%jno3*rho_dzt,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jno3_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jno3_plus_btm, topaz%jno3_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jno3denit_wc .gt. 0)        &
       used = send_data(topaz%id_jno3denit_wc,  topaz%jno3denit_wc*rho_dzt,     &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jo2 .gt. 0)                 &
       used = send_data(topaz%id_jo2,           topaz%jo2*rho_dzt,              &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jo2_plus_btm .gt. 0)        &
       used = send_data(topaz%id_jo2_plus_btm,  topaz%jo2_plus_btm * rho_dzt,     &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jpdet .gt. 0)               &
       used = send_data(topaz%id_jpdet,         topaz%jpdet*rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jpo4 .gt. 0)                &
       used = send_data(topaz%id_jpo4,          topaz%jpo4 * rho_dzt,             &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jpo4_graz .gt. 0)           &
       used = send_data(topaz%id_jpo4_graz,     topaz%jpo4_graz*rho_dzt,        &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jpo4_plus_btm .gt. 0)       &
       used = send_data(topaz%id_jpo4_plus_btm, topaz%jpo4_plus_btm * rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_cadet_calc .gt. 0)    &
       used = send_data(topaz%id_jprod_cadet_calc, topaz%jprod_cadet_calc * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_fedet .gt. 0)         &
       used = send_data(topaz%id_jprod_fedet,   topaz%jprod_fedet*rho_dzt,      &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_lithdet .gt. 0)       &
       used = send_data(topaz%id_jprod_lithdet, topaz%jprod_lithdet*rho_dzt,    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_ndet .gt. 0)          &
       used = send_data(topaz%id_jprod_ndet,    topaz%jprod_ndet*rho_dzt,       &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_nhet .gt. 0)          &
       used = send_data(topaz%id_jprod_nhet,    topaz%jprod_nhet*rho_dzt,       &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jprod_pdet .gt. 0)          &
       used = send_data(topaz%id_jprod_pdet,    topaz%jprod_pdet*rho_dzt,       &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jsdon .gt. 0)               &
       used = send_data(topaz%id_jsdon,         topaz%jsdon*rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jsdop .gt. 0)               &
       used = send_data(topaz%id_jsdop,         topaz%jsdop*rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jsidet .gt. 0)              &
       used = send_data(topaz%id_jsidet,        topaz%jsidet*rho_dzt,           &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jsio4 .gt. 0)               &
       used = send_data(topaz%id_jsio4,         topaz%jsio4 * rho_dzt,            &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_jsio4_plus_btm .gt. 0)      &
       used = send_data(topaz%id_jsio4_plus_btm, topaz%jsio4_plus_btm * rho_dzt,  &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_omega_arag .gt. 0)          &
       used = send_data(topaz%id_omega_arag,  topaz%omega_arag,                 &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (topaz%id_omega_calc .gt. 0)          &
       used = send_data(topaz%id_omega_calc,  topaz%omega_calc,                 &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
!---------------------------------------------------------------------
! Save 2D topaz type fields
!---------------------------------------------------------------------
    if (topaz%id_omega_arag .gt. 0)          &
       used = send_data(topaz%id_omega_arag,    topaz%omega_arag,               &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_omega_calc .gt. 0)          &
       used = send_data(topaz%id_omega_calc,    topaz%omega_calc,               &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_pco2surf .gt. 0)            &
       used = send_data(topaz%id_pco2surf,      topaz%pco2_csurf,               &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_alk .gt. 0)             &
       used = send_data(topaz%id_sfc_alk,       topaz%p_alk(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_cadet_arag .gt. 0)      &
       used = send_data(topaz%id_sfc_cadet_arag,topaz%p_cadet_arag(:,:,1,tau),  &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_cadet_calc .gt. 0)      &
       used = send_data(topaz%id_sfc_cadet_calc,topaz%p_cadet_calc(:,:,1,tau),  &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_chl .gt. 0)             &
       used = send_data(topaz%id_sfc_chl,       topaz%f_chl(:,:,1),             &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_chl_lg_diatoms .gt. 0)  &
       used = send_data(topaz%id_sfc_chl_lg_diatoms, topaz%c_2_n * 12.0e6 *     &
          phyto(LARGE)%theta(:,:,1) * topaz%nLg_diatoms(:,:,1),                 &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_chl_lg_nondiatoms .gt. 0)  &
       used = send_data(topaz%id_sfc_chl_lg_nondiatoms, topaz%c_2_n * 12.0e6 *  &
          phyto(LARGE)%theta(:,:,1) * max(0.0, (phyto(LARGE)%f_n(:,:,1) -       &
          topaz%nLg_diatoms(:,:,1))),                                           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_co3_ion .gt. 0)        &
       used = send_data(topaz%id_sfc_co3_ion,  topaz%f_co3_ion(:,:,1),          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_sfc_co3_sol_arag .gt. 0)        &
       used = send_data(topaz%id_sfc_co3_sol_arag,  topaz%co3_sol_arag(:,:,1),  &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_sfc_co3_sol_calc .gt. 0)        &
       used = send_data(topaz%id_sfc_co3_sol_calc,  topaz%co3_sol_calc(:,:,1),  &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (topaz%id_sfc_dic .gt. 0)             &
       used = send_data(topaz%id_sfc_dic,       topaz%p_dic(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_fed .gt. 0)             &
       used = send_data(topaz%id_sfc_fed,       topaz%p_fed(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_htotal .gt. 0)          &
       used = send_data(topaz%id_sfc_htotal,    topaz%f_htotal(:,:,1),          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_irr .gt. 0)             &
       used = send_data(topaz%id_sfc_irr,       topaz%irr_inst(:,:,1),          &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_irr_mem .gt. 0)         &
       used = send_data(topaz%id_sfc_irr_mem,   topaz%f_irr_mem(:,:,1),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_ndet .gt. 0)            &
       used = send_data(topaz%id_sfc_ndet,       topaz%p_ndet(:,:,1,tau),       &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_ndi .gt. 0)             &
       used = send_data(topaz%id_sfc_ndi,       topaz%p_ndi(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_nh4 .gt. 0)             &
       used = send_data(topaz%id_sfc_nh4,       topaz%p_nh4(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_nhet .gt. 0)             &
       used = send_data(topaz%id_sfc_nhet,       topaz%p_nhet(:,:,1,tau),       &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_nlg .gt. 0)             &
       used = send_data(topaz%id_sfc_nlg,       topaz%p_nlg(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_no3 .gt. 0)             &
       used = send_data(topaz%id_sfc_no3,       topaz%p_no3(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_nsm .gt. 0)             &
       used = send_data(topaz%id_sfc_nsm,       topaz%p_nsm(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_o2 .gt. 0)             &
       used = send_data(topaz%id_sfc_o2,       topaz%p_o2(:,:,1,tau),           &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_po4 .gt. 0)             &
       used = send_data(topaz%id_sfc_po4,       topaz%p_po4(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_silg .gt. 0)            &
       used = send_data(topaz%id_sfc_silg,      topaz%p_silg(:,:,1,tau),        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_sio4 .gt. 0)            &
       used = send_data(topaz%id_sfc_sio4,      topaz%p_sio4(:,:,1,tau),        &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (topaz%id_sfc_temp .gt. 0)            &
       used = send_data(topaz%id_sfc_temp,      Temp(:,:,1),                    &
       model_time, rmask = grid_tmask(:,:,1),& 
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

!---------------------------------------------------------------------
! Save rho_dzt
!---------------------------------------------------------------------
    if (id_rho_dzt .gt. 0)                   &
       used = send_data(id_rho_dzt,  rho_dzt,                                    &
       model_time, rmask = grid_tmask(:,:,:),& 
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

  end subroutine generic_TOPAZ_update_from_source

  ! <SUBROUTINE NAME="generic_TOPAZ_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_TOPAZ_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_TOPAZ_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
    character(len=fm_string_len), parameter :: sub_name = 'generic_TOPAZ_set_boundary_values'


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

    !nnz: Since the generic_TOPAZ_update_from_source() subroutine is called by this time
    !     the following if block is not really necessary (since this calculation is already done in source).
    !     It is only neccessary if source routine is commented out for debugging.
    !Note: In order for this to work we should NOT zero out the coupler values for generic tracers
    !      This zeroing is done for non-generic TOPAZ by calling zero_ocean_sfc.  
    !      Since the coupler values here are non-cumulative there is no need to zero them out anyway.

    if (topaz%init ) then
       !Get necessary fields
       call g_tracer_get_pointer(tracer_list,'dic'   ,'field', dic_field)
       call g_tracer_get_pointer(tracer_list,'po4'   ,'field', po4_field)
       call g_tracer_get_pointer(tracer_list,'sio4'  ,'field', sio4_field)
       call g_tracer_get_pointer(tracer_list,'alk'   ,'field', alk_field)

       call g_tracer_get_values(tracer_list,'htotal' ,'field', htotal_field,isd,jsd,ntau=1)
       call g_tracer_get_values(tracer_list,'co3_ion','field',co3_ion_field,isd,jsd,ntau=1)

       do j = jsc, jec ; do i = isc, iec  !{
          topaz%htotallo(i,j) = topaz%htotal_scale_lo * htotal_field(i,j,1)
          topaz%htotalhi(i,j) = topaz%htotal_scale_hi * htotal_field(i,j,1)
       enddo; enddo ; !} i, j

       if(topaz%reproduce_esm2) then
       call FMS_ocmip2_co2calc_old(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                            &
            dic_field(:,:,1,tau),                          &
            po4_field(:,:,1,tau),                          &
            sio4_field(:,:,1,tau),                         &
            alk_field(:,:,1,tau),                          &
            topaz%htotallo, topaz%htotalhi,                &
                                !InOut
            htotal_field(:,:,1),                           & 
                                !OUT
            co2star=co2_csurf(:,:), alpha=co2_alpha(:,:),  &
            pCO2surf=topaz%pco2_csurf(:,:), &
            co3_ion=co3_ion_field(:,:,1))
      else
       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                            &
            dic_field(:,:,1,tau),                          &
            po4_field(:,:,1,tau),                          &
            sio4_field(:,:,1,tau),                         &
            alk_field(:,:,1,tau),                          &
            topaz%htotallo, topaz%htotalhi,                &
                                !InOut
            htotal_field(:,:,1),                           & 
                                !OUT
            co2star=co2_csurf(:,:), alpha=co2_alpha(:,:),  &
            pCO2surf=topaz%pco2_csurf(:,:), &
            co3_ion=co3_ion_field(:,:,1))
      endif

       !Set fields !nnz: if These are pointers do I need to do this?
       call g_tracer_set_values(tracer_list,'htotal' ,'field',htotal_field ,isd,jsd,ntau=1)
       call g_tracer_set_values(tracer_list,'co3_ion','field',co3_ion_field,isd,jsd,ntau=1)

       call g_tracer_set_values(tracer_list,'dic','alpha',co2_alpha    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'dic','csurf',co2_csurf    ,isd,jsd)

       !!nnz: If source is called uncomment the following
       topaz%init = .false. !nnz: This is necessary since the above two calls appear in source subroutine too.
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
       !      identical with topaz code in MOM topaz%Rho_0 must be replaced with rho(i,j,1,tau)
       !      This is achieved by uncommenting the following if desired.
       !! topaz%Rho_0 = rho(i,j,1,tau)
       !      But since topaz%Rho_0 plays the role of a unit conversion factor in this module
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
       co2_sc_no(i,j) = topaz%a1_co2 + ST * (topaz%a2_co2 + ST * (topaz%a3_co2 + ST * topaz%a4_co2)) * &
            grid_tmask(i,j,1)
!       sc_no_term = sqrt(660.0 / (sc_co2 + epsln)) 
!
!       co2_alpha(i,j) = co2_alpha(i,j)* sc_no_term * topaz%Rho_0 !nnz: MOM has rho(i,j,1,tau)  
!       co2_csurf(i,j) = co2_csurf(i,j)* sc_no_term * topaz%Rho_0 !nnz: MOM has rho(i,j,1,tau) 
!
! in 'ocmip2_new' atmos_ocean_fluxes.F90 coupler formulation, the schmidt number is carried in explicitly
!
       co2_alpha(i,j) = co2_alpha(i,j) * topaz%Rho_0 !nnz: MOM has rho(i,j,1,tau)  
       co2_csurf(i,j) = co2_csurf(i,j) * topaz%Rho_0 !nnz: MOM has rho(i,j,1,tau) 

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
            exp( topaz%a_0 + topaz%a_1*ts + topaz%a_2*ts2 + topaz%a_3*ts3 + topaz%a_4*ts4 + topaz%a_5*ts5 + &
            (topaz%b_0 + topaz%b_1*ts + topaz%b_2*ts2 + topaz%b_3*ts3 + topaz%c_0*sal)*sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the 
       !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
       !  Cycles, 12, 141-163).
       !---------------------------------------------------------------------
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       o2_sc_no(i,j)  = topaz%a1_o2  + ST * (topaz%a2_o2  + ST * (topaz%a3_o2  + ST * topaz%a4_o2 )) * &
            grid_tmask(i,j,1)
       !
       !      renormalize the alpha value for atm o2
       !      data table override for o2_flux_pcair_atm is now set to 0.21 
       !
       o2_alpha(i,j) = (o2_saturation / 0.21)
       o2_csurf(i,j) = o2_field(i,j,1,tau) * topaz%Rho_0 !nnz: MOM has rho(i,j,1,tau)


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

  end subroutine generic_TOPAZ_set_boundary_values

  ! <SUBROUTINE NAME="generic_TOPAZ_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_TOPAZ_end
  !  </TEMPLATE>
  ! </SUBROUTINE>


  subroutine generic_TOPAZ_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_TOPAZ_end'
    call user_deallocate_arrays
  end subroutine generic_TOPAZ_end

  !
  !   This is an internal sub, not a public interface.
  !   Allocate all the work arrays to be used in this module.
  !
  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    !Allocate all the private arrays.

    !Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc ; CO2_dope_vec%iec = iec 
    CO2_dope_vec%jsc = jsc ; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd ; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd ; CO2_dope_vec%jed = jed    

    allocate(topaz%htotallo(isd:ied,jsd:jed))
    allocate(topaz%htotalhi(isd:ied,jsd:jed))
    

    !
    ! allocate and initialize array elements of all phytoplankton groups
    !

    do n = 1, NUM_PHYTO
       allocate(phyto(n)%jprod_n_100(isd:ied,jsd:jed))     ; phyto(n)%jprod_n_100    = 0.0
       allocate(phyto(n)%jgraz_n_100(isd:ied,jsd:jed))     ; phyto(n)%jgraz_n_100    = 0.0
       allocate(phyto(n)%chl(isd:ied,jsd:jed,nk))          ; phyto(n)%chl            = 0.0
       allocate(phyto(n)%def_fe(isd:ied,jsd:jed,nk))       ; phyto(n)%def_fe         = 0.0
       allocate(phyto(n)%def_p(isd:ied,jsd:jed,nk))        ; phyto(n)%def_p          = 0.0
       allocate(phyto(n)%f_fe(isd:ied,jsd:jed,nk))         ; phyto(n)%f_fe           = 0.0
       allocate(phyto(n)%f_n(isd:ied,jsd:jed,nk))          ; phyto(n)%f_n            = 0.0
       allocate(phyto(n)%f_p(isd:ied,jsd:jed,nk))          ; phyto(n)%f_p            = 0.0
       allocate(phyto(n)%felim(isd:ied,jsd:jed,nk))        ; phyto(n)%felim          = 0.0
       allocate(phyto(n)%irrlim(isd:ied,jsd:jed,nk))       ; phyto(n)%irrlim         = 0.0
       allocate(phyto(n)%jgraz_fe(isd:ied,jsd:jed,nk))     ; phyto(n)%jgraz_fe       = 0.0
       allocate(phyto(n)%jgraz_n(isd:ied,jsd:jed,nk))      ; phyto(n)%jgraz_n        = 0.0
       allocate(phyto(n)%jprod_fe(isd:ied,jsd:jed,nk))     ; phyto(n)%jprod_fe       = 0.0
       allocate(phyto(n)%jprod_n(isd:ied,jsd:jed,nk))      ; phyto(n)%jprod_n        = 0.0
       allocate(phyto(n)%jprod_nh4(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_nh4      = 0.0
       allocate(phyto(n)%jprod_no3(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_no3      = 0.0
       allocate(phyto(n)%jprod_po4(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_po4      = 0.0
       allocate(phyto(n)%liebig_lim(isd:ied,jsd:jed,nk))   ; phyto(n)%liebig_lim     = 0.0
       allocate(phyto(n)%mu(isd:ied,jsd:jed,nk))           ; phyto(n)%mu             = 0.0
       allocate(phyto(n)%po4lim(isd:ied,jsd:jed,nk))       ; phyto(n)%po4lim         = 0.0
       allocate(phyto(n)%q_fe_2_n(isd:ied,jsd:jed,nk))     ; phyto(n)%q_fe_2_n       = 0.0
       allocate(phyto(n)%q_p_2_n(isd:ied,jsd:jed,nk))      ; phyto(n)%q_p_2_n        = 0.0
       allocate(phyto(n)%q_p_2_n_opt(isd:ied,jsd:jed,nk))  ; phyto(n)%q_p_2_n_opt    = 0.0
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
    allocate(phyto(DIAZO)%jprod_n2(isd:ied,jsd:jed,nk))   ; phyto(DIAZO)%jprod_n2   = 0.0
    allocate(phyto(LARGE)%jgraz_sio2(isd:ied,jsd:jed,nk)) ; phyto(LARGE)%jgraz_sio2 = 0.0
    allocate(phyto(LARGE)%jprod_sio4(isd:ied,jsd:jed,nk)) ; phyto(LARGE)%jprod_sio4 = 0.0
    allocate(phyto(LARGE)%silim(isd:ied,jsd:jed,nk))      ; phyto(LARGE)%silim      = 0.0
    allocate(phyto(LARGE)%jprod_sio4_100(isd:ied,jsd:jed)); phyto(LARGE)%jprod_sio4_100 = 0.0

    allocate(topaz%chl_lg_diatoms(isd:ied, jsd:jed, 1:nk))    ; topaz%chl_lg_diatoms=0.0
    allocate(topaz%chl_lg_nondiatoms(isd:ied, jsd:jed, 1:nk)); topaz%chl_lg_nondiatoms=0.0
    allocate(topaz%co3_sol_arag(isd:ied, jsd:jed, 1:nk)); topaz%co3_sol_arag=0.0
    allocate(topaz%co3_sol_calc(isd:ied, jsd:jed, 1:nk)); topaz%co3_sol_calc=0.0
    allocate(topaz%f_alk(isd:ied, jsd:jed, 1:nk)); topaz%f_alk=0.0
    allocate(topaz%f_cadet_arag(isd:ied, jsd:jed, 1:nk)); topaz%f_cadet_arag=0.0
    allocate(topaz%f_cadet_calc(isd:ied, jsd:jed, 1:nk)); topaz%f_cadet_calc=0.0
    allocate(topaz%f_cased(isd:ied, jsd:jed, 1:nk)); topaz%f_cased=0.0
    allocate(topaz%f_chl(isd:ied, jsd:jed, 1:nk)); topaz%f_chl=0.0
    allocate(topaz%f_co3_ion(isd:ied, jsd:jed, 1:nk)); topaz%f_co3_ion=0.0
    allocate(topaz%f_dic(isd:ied, jsd:jed, 1:nk)); topaz%f_dic=0.0
    allocate(topaz%f_cadet_arag_btf(isd:ied, jsd:jed, 1:nk)); topaz%f_cadet_arag_btf=0.0
    allocate(topaz%f_cadet_calc_btf(isd:ied, jsd:jed, 1:nk)); topaz%f_cadet_calc_btf=0.0
    allocate(topaz%f_fed(isd:ied, jsd:jed, 1:nk)); topaz%f_fed=0.0
    allocate(topaz%f_fedet(isd:ied, jsd:jed, 1:nk)); topaz%f_fedet=0.0
    allocate(topaz%f_lithdet_btf(isd:ied, jsd:jed, 1:nk)); topaz%f_lithdet_btf=0.0
    allocate(topaz%f_ndet_btf(isd:ied, jsd:jed, 1:nk)); topaz%f_ndet_btf=0.0
    allocate(topaz%f_pdet_btf(isd:ied, jsd:jed, 1:nk)); topaz%f_pdet_btf=0.0
    allocate(topaz%f_sidet_btf(isd:ied, jsd:jed, 1:nk)); topaz%f_sidet_btf=0.0
    allocate(topaz%f_htotal(isd:ied, jsd:jed, 1:nk)); topaz%f_htotal=0.0
    allocate(topaz%f_irr_mem(isd:ied, jsd:jed, 1:nk)); topaz%f_irr_mem=0.0
    allocate(topaz%f_ldon(isd:ied, jsd:jed, 1:nk)); topaz%f_ldon=0.0
    allocate(topaz%f_lith(isd:ied, jsd:jed, 1:nk)); topaz%f_lith=0.0
    allocate(topaz%f_lithdet(isd:ied, jsd:jed, 1:nk)); topaz%f_lithdet=0.0
    allocate(topaz%f_ndet(isd:ied, jsd:jed, 1:nk)); topaz%f_ndet=0.0
    allocate(topaz%f_nh4(isd:ied, jsd:jed, 1:nk)); topaz%f_nh4=0.0
    allocate(topaz%f_nhet(isd:ied, jsd:jed, 1:nk)); topaz%f_nhet=0.0
    allocate(topaz%f_no3(isd:ied, jsd:jed, 1:nk)); topaz%f_no3=0.0
    allocate(topaz%f_o2(isd:ied, jsd:jed, 1:nk)); topaz%f_o2=0.0
    allocate(topaz%f_pdet(isd:ied, jsd:jed, 1:nk)); topaz%f_pdet=0.0
    allocate(topaz%f_po4(isd:ied, jsd:jed, 1:nk)); topaz%f_po4=0.0
    allocate(topaz%f_sdon(isd:ied, jsd:jed, 1:nk)); topaz%f_sdon=0.0
    allocate(topaz%f_sdop(isd:ied, jsd:jed, 1:nk)); topaz%f_sdop=0.0
    allocate(topaz%f_sidet(isd:ied, jsd:jed, 1:nk)); topaz%f_sidet=0.0
    allocate(topaz%f_silg(isd:ied, jsd:jed, 1:nk)); topaz%f_silg=0.0
    allocate(topaz%f_sio4(isd:ied, jsd:jed, 1:nk)); topaz%f_sio4=0.0
    allocate(topaz%expkT(isd:ied, jsd:jed, 1:nk)); topaz%expkT=0.0
    allocate(topaz%frac_det_prod(isd:ied, jsd:jed, 1:nk)); topaz%frac_det_prod=0.0
    allocate(topaz%irr_inst(isd:ied, jsd:jed, 1:nk)); topaz%irr_inst=0.0
    allocate(topaz%irr_mix(isd:ied, jsd:jed, 1:nk)); topaz%irr_mix=0.0
    allocate(topaz%jalk(isd:ied, jsd:jed, 1:nk)); topaz%jalk=0.0
    allocate(topaz%jalk_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jalk_plus_btm=0.0
    allocate(topaz%jcadet_arag(isd:ied, jsd:jed, 1:nk)); topaz%jcadet_arag=0.0
    allocate(topaz%jcadet_calc(isd:ied, jsd:jed, 1:nk)); topaz%jcadet_calc=0.0
    allocate(topaz%jdic(isd:ied, jsd:jed, 1:nk)); topaz%jdic=0.0
    allocate(topaz%jdic_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jdic_plus_btm=0.0
    allocate(topaz%jdiss_sio2(isd:ied, jsd:jed, 1:nk)); topaz%jdiss_sio2=0.0
    allocate(topaz%jfe_ads(isd:ied, jsd:jed, 1:nk)); topaz%jfe_ads=0.0
    allocate(topaz%jfe_des(isd:ied, jsd:jed, 1:nk)); topaz%jfe_des=0.0
    allocate(topaz%jfe_graz(isd:ied, jsd:jed, 1:nk)); topaz%jfe_graz=0.0
    allocate(topaz%jfe_coast(isd:ied, jsd:jed, 1:nk)); topaz%jfe_coast=0.0
    allocate(topaz%jfed(isd:ied, jsd:jed, 1:nk)); topaz%jfed=0.0
    allocate(topaz%jfed_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jfed_plus_btm=0.0
    allocate(topaz%jfedet(isd:ied, jsd:jed, 1:nk)); topaz%jfedet=0.0
    allocate(topaz%jldon(isd:ied, jsd:jed, 1:nk)); topaz%jldon=0.0
    allocate(topaz%jndet(isd:ied, jsd:jed, 1:nk)); topaz%jndet=0.0
    allocate(topaz%jnh4(isd:ied, jsd:jed, 1:nk)); topaz%jnh4=0.0
    allocate(topaz%jnh4_graz(isd:ied, jsd:jed, 1:nk)); topaz%jnh4_graz=0.0
    allocate(topaz%jnh4_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jnh4_plus_btm=0.0
    allocate(topaz%jnhet(isd:ied, jsd:jed, 1:nk)); topaz%jnhet=0.0
    allocate(topaz%jnitrif(isd:ied, jsd:jed, 1:nk)); topaz%jnitrif=0.0
    allocate(topaz%jno3(isd:ied, jsd:jed, 1:nk)); topaz%jno3=0.0
    allocate(topaz%jno3_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jno3_plus_btm=0.0
    allocate(topaz%jno3denit_wc(isd:ied, jsd:jed, 1:nk)); topaz%jno3denit_wc=0.0
    allocate(topaz%jo2(isd:ied, jsd:jed, 1:nk)); topaz%jo2=0.0
    allocate(topaz%jo2_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jo2_plus_btm=0.0
    allocate(topaz%jpdet(isd:ied, jsd:jed, 1:nk)); topaz%jpdet=0.0
    allocate(topaz%jpo4(isd:ied, jsd:jed, 1:nk)); topaz%jpo4=0.0
    allocate(topaz%jpo4_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jpo4_plus_btm=0.0
    allocate(topaz%jpo4_graz(isd:ied, jsd:jed, 1:nk)); topaz%jpo4_graz=0.0
    allocate(topaz%jprod_cadet_arag(isd:ied, jsd:jed, 1:nk)); topaz%jprod_cadet_arag=0.0
    allocate(topaz%jprod_cadet_calc(isd:ied, jsd:jed, 1:nk)); topaz%jprod_cadet_calc=0.0
    allocate(topaz%jprod_lithdet(isd:ied, jsd:jed, 1:nk)); topaz%jprod_lithdet=0.0
    allocate(topaz%jprod_nhet(isd:ied, jsd:jed, 1:nk)); topaz%jprod_nhet=0.0
    allocate(topaz%jprod_fedet(isd:ied, jsd:jed, 1:nk)); topaz%jprod_fedet=0.0
    allocate(topaz%jprod_ndet(isd:ied, jsd:jed, 1:nk)); topaz%jprod_ndet=0.0
    allocate(topaz%jprod_pdet(isd:ied, jsd:jed, 1:nk)); topaz%jprod_pdet=0.0
    allocate(topaz%jsdon(isd:ied, jsd:jed, 1:nk)); topaz%jsdon=0.0
    allocate(topaz%jsdop(isd:ied, jsd:jed, 1:nk)); topaz%jsdop=0.0
    allocate(topaz%jsidet(isd:ied, jsd:jed, 1:nk)); topaz%jsidet=0.0
    allocate(topaz%jsio4(isd:ied, jsd:jed, 1:nk)); topaz%jsio4=0.0
    allocate(topaz%jsio4_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jsio4_plus_btm=0.0
    allocate(topaz%nLg_diatoms(isd:ied, jsd:jed, 1:nk)); topaz%nLg_diatoms=0.0
    allocate(topaz%nLg_nondiatoms(isd:ied, jsd:jed, 1:nk)); topaz%nLg_nondiatoms=0.0
    allocate(topaz%omega_arag(isd:ied, jsd:jed, 1:nk)); topaz%omega_arag=0.0
    allocate(topaz%omega_calc(isd:ied, jsd:jed, 1:nk)); topaz%omega_calc=0.0
    allocate(topaz%q_si_2_n_Lg_diatoms(isd:ied, jsd:jed, 1:nk)); topaz%q_si_2_n_Lg_diatoms=0.0
    allocate(topaz%tot_bfe(isd:ied, jsd:jed, 1:nk)); topaz%tot_bfe=0.0
    allocate(topaz%tot_bsi(isd:ied, jsd:jed, 1:nk)); topaz%tot_bsi=0.0
    allocate(topaz%tot_doc(isd:ied, jsd:jed, 1:nk)); topaz%tot_doc=0.0
    allocate(topaz%tot_pon(isd:ied, jsd:jed, 1:nk)); topaz%tot_pon=0.0
    allocate(topaz%tot_pop(isd:ied, jsd:jed, 1:nk)); topaz%tot_pop=0.0
    allocate(topaz%tot_phyton(isd:ied, jsd:jed, 1:nk)); topaz%tot_phyton=0.0
    allocate(topaz%tot_phytop(isd:ied, jsd:jed, 1:nk)); topaz%tot_phytop=0.0
    allocate(topaz%tot_phytofe(isd:ied, jsd:jed, 1:nk)); topaz%tot_phytofe=0.0
    allocate(topaz%zt(isd:ied, jsd:jed, 1:nk)); topaz%zt=0.0

    allocate(topaz%tot_layer_int_c(isd:ied, jsd:jed, 1:nk)); topaz%tot_layer_int_c=0.0
    allocate(topaz%tot_layer_int_fe(isd:ied, jsd:jed, 1:nk)); topaz%tot_layer_int_fe=0.0
    allocate(topaz%tot_layer_int_n(isd:ied, jsd:jed, 1:nk)); topaz%tot_layer_int_n=0.0
    allocate(topaz%tot_layer_int_p(isd:ied, jsd:jed, 1:nk)); topaz%tot_layer_int_p=0.0
    allocate(topaz%tot_layer_int_si(isd:ied, jsd:jed, 1:nk)); topaz%tot_layer_int_si=0.0

    allocate(topaz%jgraz_ntot(isd:ied, jsd:jed, 1:nk)); topaz%jgraz_ntot=0.0
    allocate(topaz%jprod_fetot(isd:ied, jsd:jed, 1:nk)); topaz%jprod_fetot=0.0
    allocate(topaz%jprod_nlg_diatoms(isd:ied, jsd:jed, 1:nk)); topaz%jprod_nlg_diatoms=0.0
    allocate(topaz%jprod_nlg_nondiatoms(isd:ied, jsd:jed, 1:nk)); topaz%jprod_nlg_nondiatoms=0.0
    allocate(topaz%jprod_no3tot(isd:ied, jsd:jed, 1:nk)); topaz%jprod_no3tot=0.0
    allocate(topaz%jprod_ntot(isd:ied, jsd:jed, 1:nk)); topaz%jprod_ntot=0.0
    allocate(topaz%jprod_ptot(isd:ied, jsd:jed, 1:nk)); topaz%jprod_ptot=0.0
    allocate(topaz%jdin(isd:ied, jsd:jed, 1:nk)); topaz%jdin=0.0
    allocate(topaz%jdin_plus_btm(isd:ied, jsd:jed, 1:nk)); topaz%jdin_plus_btm=0.0

    allocate(topaz%b_alk(isd:ied, jsd:jed)); topaz%b_alk=0.0
    allocate(topaz%b_dic(isd:ied, jsd:jed)); topaz%b_dic=0.0
    allocate(topaz%b_fed(isd:ied, jsd:jed)); topaz%b_fed=0.0
    allocate(topaz%b_nh4(isd:ied, jsd:jed)); topaz%b_nh4=0.0
    allocate(topaz%b_no3(isd:ied, jsd:jed)); topaz%b_no3=0.0
    allocate(topaz%b_o2(isd:ied, jsd:jed)); topaz%b_o2=0.0
    allocate(topaz%b_po4(isd:ied, jsd:jed)); topaz%b_po4=0.0
    allocate(topaz%b_sio4(isd:ied, jsd:jed)); topaz%b_sio4=0.0
    allocate(topaz%co2_csurf(isd:ied, jsd:jed)); topaz%co2_csurf=0.0
    allocate(topaz%co2_alpha(isd:ied, jsd:jed)); topaz%co2_alpha=0.0
    allocate(topaz%f_alk_int_100(isd:ied, jsd:jed)); topaz%f_alk_int_100=0.0
    allocate(topaz%f_dic_int_100(isd:ied, jsd:jed)); topaz%f_dic_int_100=0.0
    allocate(topaz%f_din_int_100(isd:ied, jsd:jed)); topaz%f_din_int_100=0.0
    allocate(topaz%f_fed_int_100(isd:ied, jsd:jed)); topaz%f_fed_int_100=0.0
    allocate(topaz%f_po4_int_100(isd:ied, jsd:jed)); topaz%f_po4_int_100=0.0
    allocate(topaz%f_sio4_int_100(isd:ied, jsd:jed)); topaz%f_sio4_int_100=0.0
    allocate(topaz%fcased_burial(isd:ied, jsd:jed)); topaz%fcased_burial=0.0
    allocate(topaz%fcased_redis(isd:ied, jsd:jed)); topaz%fcased_redis=0.0
    allocate(topaz%fcadet_arag_btm(isd:ied, jsd:jed)); topaz%fcadet_arag_btm=0.0
    allocate(topaz%fcadet_calc_btm(isd:ied, jsd:jed)); topaz%fcadet_calc_btm=0.0
    allocate(topaz%ffe_sed(isd:ied, jsd:jed)); topaz%ffe_sed=0.0
    allocate(topaz%ffedet_100(isd:ied, jsd:jed)); topaz%ffedet_100=0.0
    allocate(topaz%ffedet_btm(isd:ied, jsd:jed)); topaz%ffedet_btm=0.0
    allocate(topaz%flithdet_100(isd:ied, jsd:jed)); topaz%flithdet_100=0.0
    allocate(topaz%flithdet_btm(isd:ied, jsd:jed)); topaz%flithdet_btm=0.0
    allocate(topaz%fnfeso4red_sed(isd:ied, jsd:jed)); topaz%fnfeso4red_sed=0.0
    allocate(topaz%fndet_100(isd:ied, jsd:jed)); topaz%fndet_100=0.0
    allocate(topaz%fndet_btm(isd:ied, jsd:jed)); topaz%fndet_btm=0.0
    allocate(topaz%fno3denit_sed(isd:ied, jsd:jed)); topaz%fno3denit_sed=0.0
    allocate(topaz%fnoxic_sed(isd:ied, jsd:jed)); topaz%fnoxic_sed=0.0
    allocate(topaz%fpdet_100(isd:ied, jsd:jed)); topaz%fpdet_100=0.0
    allocate(topaz%fpdet_btm(isd:ied, jsd:jed)); topaz%fpdet_btm=0.0
    allocate(topaz%fcadet_arag_100(isd:ied, jsd:jed)); topaz%fcadet_arag_100=0.0
    allocate(topaz%fcadet_calc_100(isd:ied, jsd:jed)); topaz%fcadet_calc_100=0.0
    allocate(topaz%fsidet_100(isd:ied, jsd:jed)); topaz%fsidet_100=0.0
    allocate(topaz%fsidet_btm(isd:ied, jsd:jed)); topaz%fsidet_btm=0.0
    allocate(topaz%jalk_100(isd:ied, jsd:jed)); topaz%jalk_100=0.0
    allocate(topaz%jdic_100(isd:ied, jsd:jed)); topaz%jdic_100=0.0
    allocate(topaz%jdin_100(isd:ied, jsd:jed)); topaz%jdin_100=0.0
    allocate(topaz%jfed_100(isd:ied, jsd:jed)); topaz%jfed_100=0.0
    allocate(topaz%jnitrif_100(isd:ied, jsd:jed)); topaz%jnitrif_100=0.0
    allocate(topaz%jpo4_100(isd:ied, jsd:jed)); topaz%jpo4_100=0.0
    allocate(topaz%jprod_cadet_arag_100(isd:ied, jsd:jed)); topaz%jprod_cadet_arag_100=0.0
    allocate(topaz%jprod_cadet_calc_100(isd:ied, jsd:jed)); topaz%jprod_cadet_calc_100=0.0
    allocate(topaz%jprod_fetot_100(isd:ied, jsd:jed)); topaz%jprod_fetot_100=0.0
    allocate(topaz%jprod_nlg_diatoms_100(isd:ied, jsd:jed)); topaz%jprod_nlg_diatoms_100=0.0
    allocate(topaz%jprod_nlg_nondiatoms_100(isd:ied, jsd:jed)); topaz%jprod_nlg_nondiatoms_100=0.0
    allocate(topaz%jprod_ntot_100(isd:ied, jsd:jed)); topaz%jprod_ntot_100=0.0
    allocate(topaz%jprod_no3tot_100(isd:ied, jsd:jed)); topaz%jprod_no3tot_100=0.0
    allocate(topaz%jprod_ptot_100(isd:ied, jsd:jed)); topaz%jprod_ptot_100=0.0
    allocate(topaz%jsio4_100(isd:ied, jsd:jed)); topaz%jsio4_100=0.0
    allocate(topaz%nongas_source_c(isd:ied, jsd:jed)); topaz%nongas_source_c=0.0
    allocate(topaz%nongas_source_fe(isd:ied, jsd:jed)); topaz%nongas_source_fe=0.0
    allocate(topaz%nongas_source_n(isd:ied, jsd:jed)); topaz%nongas_source_n=0.0
    allocate(topaz%o2min(isd:ied, jsd:jed)); topaz%o2min=0.0
    allocate(topaz%pco2_csurf(isd:ied, jsd:jed)); topaz%pco2_csurf=0.0
    allocate(topaz%wc_vert_int_c(isd:ied, jsd:jed)); topaz%wc_vert_int_c=0.0
    allocate(topaz%wc_vert_int_dic(isd:ied, jsd:jed)); topaz%wc_vert_int_dic=0.0
    allocate(topaz%wc_vert_int_jfe_coast(isd:ied, jsd:jed)); topaz%wc_vert_int_jfe_coast=0.0
    allocate(topaz%wc_vert_int_jno3denit(isd:ied, jsd:jed)); topaz%wc_vert_int_jno3denit=0.0
    allocate(topaz%wc_vert_int_nfix(isd:ied, jsd:jed)); topaz%wc_vert_int_nfix=0.0
    allocate(topaz%z_o2min(isd:ied, jsd:jed)); topaz%z_o2min=0.0
    allocate(topaz%z_sat_arag(isd:ied, jsd:jed)); topaz%z_sat_arag=0.0
    allocate(topaz%z_sat_calc(isd:ied, jsd:jed)); topaz%z_sat_calc=0.0

    allocate(topaz%mask_z_sat_arag(isd:ied, jsd:jed)); topaz%mask_z_sat_arag = .FALSE.
    allocate(topaz%mask_z_sat_calc(isd:ied, jsd:jed)); topaz%mask_z_sat_calc = .FALSE.

  end subroutine user_allocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays
    integer n

    deallocate(topaz%htotalhi,topaz%htotallo)

    do n = 1, NUM_PHYTO
       deallocate(phyto(n)%jgraz_n_100)
       deallocate(phyto(n)%jprod_n_100)
       deallocate(phyto(n)%chl)
       deallocate(phyto(n)%f_fe)
       deallocate(phyto(n)%def_fe)
       deallocate(phyto(n)%def_p)
       deallocate(phyto(n)%felim)
       deallocate(phyto(n)%irrlim)
       deallocate(phyto(n)%jgraz_fe)
       deallocate(phyto(n)%jgraz_n)
       deallocate(phyto(n)%jprod_fe)
       deallocate(phyto(n)%jprod_n)
       deallocate(phyto(n)%jprod_nh4)
       deallocate(phyto(n)%jprod_no3)
       deallocate(phyto(n)%jprod_po4)
       deallocate(phyto(n)%liebig_lim)
       deallocate(phyto(n)%mu)
       deallocate(phyto(n)%f_n)
       deallocate(phyto(n)%f_p)
       deallocate(phyto(n)%po4lim)
       deallocate(phyto(n)%q_fe_2_n)
       deallocate(phyto(n)%q_p_2_n)
       deallocate(phyto(n)%q_p_2_n_opt)
       deallocate(phyto(n)%theta)
    enddo
    do n = 2, NUM_PHYTO
       deallocate(phyto(n)%nh4lim)
       deallocate(phyto(n)%no3lim)
    enddo
    deallocate(phyto(DIAZO)%jprod_n2)
    deallocate(phyto(LARGE)%jgraz_sio2)
    deallocate(phyto(LARGE)%jprod_sio4)
    deallocate(phyto(LARGE)%jprod_sio4_100)
    deallocate(phyto(LARGE)%silim)

    deallocate(&
         topaz%chl_lg_diatoms,&
         topaz%chl_lg_nondiatoms,&
         topaz%co3_sol_arag,&
         topaz%co3_sol_calc,&
         topaz%f_alk,&
         topaz%f_cadet_arag,&
         topaz%f_cadet_arag_btf,&
         topaz%f_cadet_calc,&
         topaz%f_cadet_calc_btf,&
         topaz%f_cased,&
         topaz%f_chl,&
         topaz%f_co3_ion,&
         topaz%f_dic,&
         topaz%f_fed,&
         topaz%f_fedet,&
         topaz%f_htotal,&
         topaz%f_irr_mem,&
         topaz%f_ldon,&
         topaz%f_lith,&
         topaz%f_lithdet,&
         topaz%f_lithdet_btf,&
         topaz%f_ndet,&
         topaz%f_ndet_btf,&
         topaz%f_nh4,&
         topaz%f_nhet,&
         topaz%f_no3,&
         topaz%f_o2,&
         topaz%f_pdet,&
         topaz%f_pdet_btf,&
         topaz%f_po4,&
         topaz%f_sdon,&
         topaz%f_sdop,&
         topaz%f_sidet,&
         topaz%f_sidet_btf,&
         topaz%f_silg,&
         topaz%f_sio4,&
         topaz%expkT,&
         topaz%frac_det_prod,&
         topaz%irr_inst,&
         topaz%irr_mix,&
         topaz%jalk,&
         topaz%jalk_plus_btm,&
         topaz%jcadet_arag,&
         topaz%jcadet_calc,&
         topaz%jdic,&
         topaz%jdic_plus_btm,&
         topaz%jdin,&
         topaz%jdin_plus_btm,&
         topaz%jdiss_sio2,&
         topaz%jfe_ads,&
         topaz%jfe_des,&
         topaz%jfe_graz,&
         topaz%jfe_coast,&
         topaz%jfed,&
         topaz%jfed_plus_btm,&
         topaz%jfedet,&
         topaz%jgraz_ntot,&
         topaz%jldon,&
         topaz%jndet,&
         topaz%jnh4,&
         topaz%jnh4_graz,&
         topaz%jnh4_plus_btm,&
         topaz%jnhet,&
         topaz%jnitrif,&
         topaz%jno3,&
         topaz%jno3_plus_btm,&
         topaz%jno3denit_wc,&
         topaz%jo2,&
         topaz%jo2_plus_btm,&
         topaz%jpdet,&
         topaz%jpo4,&
         topaz%jpo4_graz,&
         topaz%jpo4_plus_btm,&
         topaz%jprod_cadet_arag,&
         topaz%jprod_cadet_calc,&
         topaz%jprod_lithdet,&
         topaz%jprod_fedet,&
         topaz%jprod_fetot,&
         topaz%jprod_ndet,&
         topaz%jprod_nhet,&
         topaz%jprod_nlg_diatoms,&
         topaz%jprod_nlg_nondiatoms,&
         topaz%jprod_no3tot,&
         topaz%jprod_ntot,&
         topaz%jprod_pdet,&
         topaz%jprod_ptot,&
         topaz%jsdon,&
         topaz%jsdop,&
         topaz%jsidet,&
         topaz%jsio4,&
         topaz%jsio4_plus_btm,&
         topaz%nLg_diatoms,&
         topaz%nLg_nondiatoms,&
         topaz%omega_arag,&
         topaz%omega_calc,&
         topaz%q_si_2_n_Lg_diatoms,&
         topaz%tot_layer_int_c,&
         topaz%tot_layer_int_fe,&
         topaz%tot_layer_int_n,&
         topaz%tot_layer_int_p,&
         topaz%tot_layer_int_si,&
         topaz%tot_bfe,&
         topaz%tot_bsi,&
         topaz%tot_doc,&
         topaz%tot_pon,&
         topaz%tot_pop,&
         topaz%tot_phyton,&
         topaz%tot_phytop,&
         topaz%tot_phytofe,&
         topaz%zt )

    deallocate(&
         topaz%b_alk,&
         topaz%b_dic,&
         topaz%b_fed,&
         topaz%b_nh4,&
         topaz%b_no3,&
         topaz%b_o2,&
         topaz%b_po4,&
         topaz%b_sio4,&
         topaz%co2_csurf,topaz%co2_alpha,&
         topaz%f_alk_int_100,&
         topaz%f_dic_int_100,&
         topaz%f_din_int_100,&
         topaz%f_fed_int_100,&
         topaz%f_po4_int_100,&
         topaz%f_sio4_int_100,&
         topaz%fcadet_arag_100,&
         topaz%fcadet_arag_btm,&
         topaz%fcadet_calc_100,&
         topaz%fcadet_calc_btm,&
         topaz%fcased_burial,&
         topaz%fcased_redis,&
         topaz%ffe_sed,&
         topaz%ffedet_100,&
         topaz%ffedet_btm,&
         topaz%flithdet_100,&
         topaz%flithdet_btm,&
         topaz%fnfeso4red_sed,&
         topaz%fndet_100,&
         topaz%fndet_btm,&
         topaz%fno3denit_sed,&
         topaz%fnoxic_sed,&
         topaz%fpdet_100,&
         topaz%fpdet_btm,&
         topaz%fsidet_100,&
         topaz%fsidet_btm,&
         topaz%jalk_100,&
         topaz%jdic_100,&
         topaz%jdin_100,&
         topaz%jfed_100,&
         topaz%jnitrif_100,&
         topaz%jpo4_100,&
         topaz%jprod_cadet_arag_100,&
         topaz%jprod_cadet_calc_100,&
         topaz%jprod_fetot_100,&
         topaz%jprod_nlg_diatoms_100,&
         topaz%jprod_nlg_nondiatoms_100,&
         topaz%jprod_ntot_100,&
         topaz%jprod_no3tot_100,&
         topaz%jprod_ptot_100,&
         topaz%jsio4_100,&
         topaz%nongas_source_c,&
         topaz%nongas_source_fe,&
         topaz%nongas_source_n,&
         topaz%o2min,&
         topaz%pco2_csurf,&
         topaz%wc_vert_int_c,&
         topaz%wc_vert_int_dic,&
         topaz%wc_vert_int_nfix,&
         topaz%wc_vert_int_jfe_coast,&
         topaz%wc_vert_int_jno3denit,&
         topaz%z_o2min,&
         topaz%z_sat_arag,&
         topaz%z_sat_calc,& 
         topaz%mask_z_sat_arag,&
         topaz%mask_z_sat_calc )

  end subroutine user_deallocate_arrays


end module generic_TOPAZ

