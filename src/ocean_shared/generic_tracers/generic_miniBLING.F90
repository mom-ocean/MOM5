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
!   This module contains the generic version of miniBLING.
!   It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
!
!   WARNING: although the core components of the model (PO4, Fed, DOP, O2) 
!   have been reasonably well tested, the other components should be viewed as
!   developmental at this point. There may still be some bugs.
!
!   Also, the growth parameters have been tuned to produce a reasonable simulation
!   of PO4 and chl in a 3-degree CORE-forced version of MOM4p1. It is unlikely that 
!   these parameter choices will produce satisfactory simulations in other physical
!   model configurations, and may need to be adjusted. 
! </OVERVIEW>
!
!<DESCRIPTION>
!   Biogeochemistry with Light, Iron, Nutrient and Gas version zero (BLINGv0) 
!   includes an implicit ecological model of growth limitation by light,
!   temperature, phosphate and iron, along with dissolved organic
!   phosphorus and O2 pools.
!   Food web processing in the euphotic zone and remineralization/
!   dissolution through the ocean interior are handled as in Dunne et al. 
!   (2005).  O2 equilibria and gas exchange follow OCMIP2 protocols.
!   Additional functionality comes from an optional carbon cycle that is 
!   non-interactive, i.e. does not change the core miniBLING behaviour, as
!   well as tracers for radiocarbon (14c), a decomposition of carbon 
!   components by gas exchange and remineralization (carbon_pre), a
!   decomposition of oxygen as preformed and total (o2_pre), saturation and 
!   consumed, and a decomposition of phosphate as preformed and remineralized 
!   (po4_pre).
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! This model is available for public use. 
! The current version is BLINGv0. The version number refers to the core
! model behaviour; additional tracers exist in different iterations of the
! module. In publications it should be referenced as:
! Galbraith, E.D., Gnanadesikan, A., Dunne, J. and Hiscock, M. 2010.
! Regional impacts of iron-light colimitation in a global 
! biogeochemical model. Biogeosciences , 7, 1043-1064.
!
! All parameter values are as described in this paper.
! Note that this reference is only for the core model components, and
! does not include any of the additional functionalities, which remain
! undocumented. Please contact Eric Galbraith (eric.galbraith@mcgill.ca)
! for more information.
! </REFERENCE>
!
! <DEVELOPER_NOTES>
! This code was originally developed based on the template of Perth generic TOPAZ code.
! </DEVELOPER_NOTES>
! </INFO>
!
!<NAMELIST NAME="generic_miniBLING_nml">
!
!  <DATA NAME="do_14c" TYPE="logical">
!  If true, then simulate radiocarbon. Includes 2 prognostic tracers, DI14C
! and DO14C. Requires that do_carbon = .true. 
!  </DATA> 
!
!  <DATA NAME="do_carbon" TYPE="logical">
!  If true, then simulate the carbon cycle based on strict stoichiometry
! of C:P. Includes 1 prognostic tracer, DIC.
!  </DATA> 
!
!</NAMELIST>
!
!----------------------------------------------------------------

#include <fms_platform.h>

module generic_miniBLING_mod

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len, fm_path_name_len, fm_field_name_len
  use mpp_mod,           only: mpp_error, NOTE, FATAL
  use mpp_mod,           only: stdout
  use time_manager_mod,  only: time_type
  use diag_manager_mod,  only: register_diag_field, send_data 
  use constants_mod,     only: WTMCO2, WTMO2
  use data_override_mod, only: data_override

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer
  use g_tracer_utils, only : g_tracer_get_common
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_get_values, g_tracer_column_int, g_tracer_flux_at_depth

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector

  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_miniBLING'
  character(len=fm_string_len), parameter :: package_name   = 'generic_minibling'

  public do_generic_miniBLING
  public generic_miniBLING_register
  public generic_miniBLING_init
  public generic_miniBLING_register_diag
  public generic_miniBLING_update_from_coupler
  public generic_miniBLING_diag
  public generic_miniBLING_update_from_source
  public generic_miniBLING_update_from_bottom
  public generic_miniBLING_set_boundary_values
  public generic_miniBLING_end

  !The following logical for using this module is overwritten 
  ! generic_tracer_nml namelist
  logical, save :: do_generic_miniBLING = .false.
  logical, save :: module_is_initialized = .false.

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30


  !
  !The following two types contain all the parameters and arrays used in this module.

  type generic_miniBLING_type

     character(len=fm_string_len)               :: name           = '_'
     character(len=fm_field_name_len)           :: suffix         = ' '
     character(len=fm_field_name_len)           :: long_suffix    = ' '
     logical                                    :: prevent_neg_o2 = .true.

  ! Turn on additional complexity. Most relevant diagnostic variables and all  
  ! tracers are not activated unless the appropriate switch is set to true.
  
     logical                                    :: do_14c                  = .true.     ! Requires do_carbon = .true.
     logical                                    :: do_carbon               = .true.
     real                                       :: min_frac_pop  = 0.0        ! Set to 1 to turn off recycling

     character(len=fm_string_len)               :: alk_scheme    = 'normal'   ! Specify the scheme to use for calculating alkalinity
     character(len=fm_string_len)               :: biomass_type  = 'single'   ! Specify the scheme to use for calculating biomass
     real                                       :: alk_slope     = 32.0e-06   ! Slope of alk:salt equation
     real                                       :: alk_intercept = 1200.0e-06 ! Intercept of alk:salt equation
     real                                       :: alpha_photo          ! Quantum yield under low light
     real                                       :: c_2_p                ! Carbon to Phosphorus ratio
     real                                       :: chl_min              ! Minimum chl concentration allowed (for numerical stability)
     logical                                    :: fe_is_prognostic = .false. ! Set whether Fed is prognostic or diagnostic
     logical                                    :: fe_is_diagnostic = .false. ! Set whether Fed is diagnostic or data
     real                                       :: fe_restoring     = 10.0    ! Restoring time scale, in days, if Fed is diagnostic
     real                                       :: fe_coastal       = 2.0e-09 ! Coastal iron concentration, in mol/kg, if Fed is diagnostic
     real                                       :: fe_coastal_depth = 200.0   ! Coastal depth, in meters, if Fed is diagnostic
     real                                       :: fe_2_p_max           ! Iron to Phosphate uptake ratio scaling
     real                                       :: def_fe_min = 0.0     ! Minimum for iron deficiency
     real                                       :: fe_2_p_sed           ! Iron to Phosphorus ratio in sediments
     real                                       :: felig_bkg            ! Iron ligand concentration
     real                                       :: gamma_biomass        ! Biomass adjustment timescale
     real                                       :: gamma_irr_mem        ! Photoadaptation timescale
     real                                       :: gamma_pop            ! Patriculate Organic Phosphorus decay
     real                                       :: half_life_14c        ! Radiocarbon half-life
     real                                       :: k_fe_2_p             ! Fe:P half-saturation constant
     real                                       :: k_fe_uptake          ! Iron half-saturation concentration
     real                                       :: k_o2                 ! Oxygen half-saturation concentration
     real                                       :: k_po4                ! Phosphate half-saturation concentration
     real                                       :: k_po4_recycle        ! Phosphate half-saturation concentration
     real                                       :: kappa_eppley         ! Temperature dependence
     real                                       :: kappa_remin          ! Temperature dependence for particle fractionation
     real                                       :: kfe_inorg            ! Iron scavenging, 2nd order
     real                                       :: kfe_eq_lig_max       ! Maximum light-dependent iron ligand stability constant
     real                                       :: kfe_eq_lig_min       ! Minimum light-dependent iron ligand stability constant
     real                                       :: kfe_eq_lig_irr       ! Irradiance scaling for iron ligand stability constant
     real                                       :: kfe_eq_lig_femin     ! Low-iron threshold for ligand stability constant
     real                                       :: kfe_org              ! Iron scavenging, 1st order
     real                                       :: lambda0              ! Total mortality rate constant
     real                                       :: lambda_14c           ! Radiocarbon decay rate
     real                                       :: mass_2_p             ! Organic matter mass to Phosphorus ratio
     real                                       :: o2_2_p               ! Oxygen to Phosphorus ratio
     real                                       :: o2_min               ! Anaerobic respiration threshold
     real                                       :: P_star               ! Pivotal phytoplankton concentration
     real                                       :: pc_0                 ! Maximum carbon-specific growth rate
     real                                       :: phi_lg               ! Fraction of small phytoplankton converted to detritus
     real                                       :: phi_sm               ! Fraction of large phytoplankton converted to detritus
     real                                       :: po4_min              ! Minimum PO4 concentration
     real                                       :: remin_min            ! Minimum remineralization under low O2
     real                                       :: thetamax_hi          ! Maximum Chl:C ratio when iron-replete
     real                                       :: thetamax_lo          ! Maximum Chl:C ratio when iron-limited
     real                                       :: wsink_acc            ! Sinking rate acceleration with depth
     real                                       :: wsink0               ! Sinking rate at surface
     real                                       :: wsink0_z             ! Depth to which sinking rate remains constant

     real                                       :: htotal_scale_lo
     real                                       :: htotal_scale_hi
     real                                       :: htotal_in
     real                                       :: Rho_0
     real                                       :: a_0
     real                                       :: a_1
     real                                       :: a_2
     real                                       :: a_3
     real                                       :: a_4
     real                                       :: a_5
     real                                       :: b_0
     real                                       :: b_1
     real                                       :: b_2
     real                                       :: b_3
     real                                       :: c_0
     real                                       :: a1_co2
     real                                       :: a2_co2
     real                                       :: a3_co2
     real                                       :: a4_co2
     real                                       :: a1_o2
     real                                       :: a2_o2
     real                                       :: a3_o2
     real                                       :: a4_o2

!
!       the following arrays are used for calculation diagnostic integrals and fluxes at depth
!

     real,    dimension(:,:,:), _ALLOCATABLE       :: wrk_3d               _NULL
     real,    dimension(:,:),   _ALLOCATABLE       :: wrk_2d               _NULL
     integer, dimension(:,:),   _ALLOCATABLE       :: k_lev                _NULL
     real,    dimension(:,:),   _ALLOCATABLE       :: integral             _NULL
     real,    dimension(:,:),   _ALLOCATABLE       :: flux                 _NULL

! The prefix nomenclature is as follows:
! "f_t" = a "field", generally a working array for the concentration of tracer t
! "jt_process" = a source/sink term for tracer t due to a biogeochemical process.
!     * Note, j terms are in units of mol kg-1 in the code, but are saved to the diagnostic
!       file as layer integrals (i.e. multiplied by the layer thickness/density) 
! "b_t" = the flux of tracer t out of the ocean bottom
! "p_t" = a pointer, generally to the concentration of a tracer t

     real, dimension(:,:,:), _ALLOCATABLE       :: biomass_p_ts         _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: def_fe               _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: expkT                _NULL
     real, dimension(:,:,:), pointer            :: p_biomass_p          => NULL()
     real, dimension(:,:,:), _ALLOCATABLE       :: f_chl                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: f_fed                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: f_fed_data           _NULL
     real, dimension(:,:,:), pointer            :: p_phyto_lg           => NULL()
     real, dimension(:,:,:), pointer            :: p_phyto_sm           => NULL()
     real, dimension(:,:,:), pointer            :: p_htotal             => NULL()
     real, dimension(:,:,:), pointer            :: p_irr_mem            => NULL()
     real, dimension(:,:,:), _ALLOCATABLE       :: f_o2                 _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: f_po4                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: fe_2_p_uptake        _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: feprime              _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: fpofe                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: fpop                 _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: frac_lg              _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: frac_pop             _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: irr_inst             _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: irr_mix              _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: irrk                 _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jfe_ads_inorg        _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jfe_ads_org          _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jfe_recycle          _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jfe_reminp           _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jfe_uptake           _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jo2                  _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jp_recycle           _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jp_reminp            _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jp_uptake            _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jpo4                 _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jfeop                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jpop                 _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: kfe_eq_lig           _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: mu                   _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: pc_m                 _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: theta                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: thetamax_fe          _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: wsink                _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: zremin               _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: zbot                 _NULL

     real, dimension(:,:), _ALLOCATABLE         :: b_fed                _NULL
     real, dimension(:,:), _ALLOCATABLE         :: b_o2                 _NULL
     real, dimension(:,:), _ALLOCATABLE         :: b_po4                _NULL
     real, dimension(:,:), _ALLOCATABLE         :: fe_burial            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: ffe_sed              _NULL
     real, dimension(:,:), _ALLOCATABLE         :: o2_saturation        _NULL

     real, dimension(:,:), _ALLOCATABLE         :: b_dic                _NULL
     real, dimension(:,:), _ALLOCATABLE         :: co2_alpha            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: co2_csurf            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: htotallo             _NULL
     real, dimension(:,:), _ALLOCATABLE         :: htotalhi             _NULL
     real, dimension(:,:), _ALLOCATABLE         :: pco2_surf            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: surf_temp            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: surf_salt            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: surf_alk             _NULL
     real, dimension(:,:), _ALLOCATABLE         :: surf_po4             _NULL
     real, dimension(:,:), _ALLOCATABLE         :: surf_sio4            _NULL
     real, dimension(:,:), _ALLOCATABLE         :: surf_dic             _NULL

     real, dimension(:,:,:), _ALLOCATABLE       :: c14_2_p              _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: fpo14c               _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: j14c_decay_dic       _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: j14c_reminp          _NULL
     real, dimension(:,:,:), _ALLOCATABLE       :: jdi14c               _NULL
     real, dimension(:,:), _ALLOCATABLE         :: b_di14c              _NULL
     real, dimension(:,:), _ALLOCATABLE         :: c14o2_alpha          _NULL
     real, dimension(:,:), _ALLOCATABLE         :: c14o2_csurf          _NULL
     
     real, dimension(:,:,:,:), pointer          :: p_fed        => NULL()
     real, dimension(:,:,:), pointer            :: p_fed_diag   => NULL()
     real, dimension(:,:,:,:), pointer          :: p_o2         => NULL()
     real, dimension(:,:,:,:), pointer          :: p_po4        => NULL()

     real, dimension(:,:,:,:), pointer          :: p_di14c      => NULL()
     real, dimension(:,:,:,:), pointer          :: p_dic        => NULL()
     
     character(len=fm_string_len)               :: ice_restart_file
     character(len=fm_string_len)               :: ocean_restart_file
     character(len=fm_string_len)               :: IC_file

     real                                       :: diag_depth          = 100.0  ! Depth over which to integrate and at which to get flux
                                                                                ! for diagnostics
     character(len=16)                          :: diag_depth_str      = ' '    ! String to hold diag depth

     integer                                    :: id_b_dic            = -1     ! Bottom flux of DIC
     integer                                    :: id_b_di14c          = -1     ! Bottom flux of DI14C
     integer                                    :: id_b_fed            = -1     ! Bottom flux of Fe
     integer                                    :: id_b_o2             = -1     ! Bottom flux of O2
     integer                                    :: id_b_po4            = -1     ! Bottom flux of PO4
     integer                                    :: id_biomass_p_ts     = -1     ! Instantaneous P concentration in biomass
     integer                                    :: id_c14_2_p          = -1     ! DI14C to PO4 uptake ratio
     integer                                    :: id_c14o2_csurf      = -1     ! Surface water 14CO2*
     integer                                    :: id_c14o2_alpha      = -1     ! Surface water 14CO2* solubility 
     integer                                    :: id_co2_csurf        = -1     ! Surface water CO2*
     integer                                    :: id_co2_alpha        = -1     ! Surface water CO2* solubility
     integer                                    :: id_def_fe           = -1     ! Iron deficiency term 
     integer                                    :: id_expkT            = -1     ! Temperature dependence 
     integer                                    :: id_fe_2_p_uptake    = -1     ! Fed:PO4 of instantaneous uptake
     integer                                    :: id_feprime          = -1     ! Free (unbound) iron concentration
     integer                                    :: id_fe_burial        = -1     ! Flux of iron to sediment as particulate
     integer                                    :: id_ffe_sed          = -1     ! Sediment iron efflux
     integer                                    :: id_fpofe            = -1     ! POFe sinking flux
     integer                                    :: id_fpo14c           = -1     ! PO14C sinking flux
     integer                                    :: id_fpop             = -1     ! POP sinking flux
     integer                                    :: id_fpop_depth       = -1     ! POP sinking flux at depth
     integer                                    :: id_frac_lg          = -1     ! Fraction of production by large phytoplankton
     integer                                    :: id_frac_pop         = -1     ! Fraction of uptake converted to particulate
     integer                                    :: id_irr_inst         = -1     ! Instantaneous irradiance 
     integer                                    :: id_irr_mix          = -1     ! Mixed layer irradiance 
     integer                                    :: id_irrk             = -1     ! Effective susceptibility to light limitation 
     integer                                    :: id_j14c_decay_dic   = -1     ! Radioactive decay of DI14C
     integer                                    :: id_j14c_reminp      = -1     ! 14C particle remineralization layer integral
     integer                                    :: id_jdi14c           = -1     ! DI14C source layer integral
     integer                                    :: id_jfe_ads_inorg    = -1     ! Iron adsorption (2nd order) layer integral
     integer                                    :: id_jfe_ads_org      = -1     ! Iron adsorption to fpop layer integral
     integer                                    :: id_jfe_recycle      = -1     ! Iron fast recycling layer integral
     integer                                    :: id_jfe_reminp       = -1     ! Iron particle remineralization layer integral
     integer                                    :: id_jfe_uptake       = -1     ! Iron uptake layer integral
     integer                                    :: id_jo2              = -1     ! O2 source layer integral
     integer                                    :: id_jo2_depth        = -1     ! Depth integral of O2 source
     integer                                    :: id_jp_recycle       = -1     ! Phosphorus fast recycling layer integral
     integer                                    :: id_jp_recycle_depth = -1     ! Depth integral of Phosphorus fast recycling
     integer                                    :: id_jp_reminp        = -1     ! Phosphorus particle remineralization layer integral
     integer                                    :: id_jp_reminp_depth  = -1     ! Depth integral of Phosphorus particle remineralization
     integer                                    :: id_jp_uptake        = -1     ! Phosphorus uptake layer integral
     integer                                    :: id_jp_uptake_depth  = -1     ! Depth integral of Phosphorus uptake
     integer                                    :: id_jpo4             = -1     ! PO4 source layer integral
     integer                                    :: id_jpo4_depth       = -1     ! Depth integral of PO4 source layer integral
     integer                                    :: id_jfeop            = -1     ! Particulate organic iron source layer integral
     integer                                    :: id_jpop             = -1     ! Particulate organic phosphorus source layer integral
     integer                                    :: id_kfe_eq_lig       = -1     ! Iron-ligand stability constant
     integer                                    :: id_mu               = -1     ! Growth rate after respiratory loss(carbon specific)
     integer                                    :: id_o2_saturation    = -1     ! Surface water O2 saturation
     integer                                    :: id_pc_m             = -1     ! Light-saturated maximum photosynthesis rate (carbon specific)
     integer                                    :: id_pco2_surf        = -1     ! Surface water pCO2
     integer                                    :: id_temp_co2calc     = -1     ! Surface temp for co2calc
     integer                                    :: id_salt_co2calc     = -1     ! Surface salt for co2calc
     integer                                    :: id_alk_co2calc      = -1     ! Surface temp for co2calc
     integer                                    :: id_po4_co2calc      = -1     ! Surface temp for co2calc
     integer                                    :: id_sio4_co2calc     = -1     ! Surface temp for co2calc
     integer                                    :: id_dic_co2calc      = -1     ! Surface temp for co2calc
     integer                                    :: id_theta            = -1     ! Chl:C ratio
     integer                                    :: id_thetamax_fe      = -1     ! Iron-limited maximum Chl:C ratio
     integer                                    :: id_wsink            = -1     ! Sinking rate
     integer                                    :: id_zremin           = -1     ! Remineralization length scale
     integer                                    :: id_fed_data         = -1     ! Dissolved Iron data

     integer                                    :: id_di14c_surf       = -1     ! Surface dissolved inorganic radiocarbon Prognostic tracer
     integer                                    :: id_dic_surf         = -1     ! Surface dissolved inorganic carbon Prognostic tracer
     integer                                    :: id_fed_surf         = -1     ! Surface dissolved Iron Prognostic tracer
     integer                                    :: id_o2_surf          = -1     ! Surface oxygen Prognostic tracer
     integer                                    :: id_po4_surf         = -1     ! Surface phosphate Prognostic tracer
     integer                                    :: id_di14c_depth      = -1     ! Depth integral of dissolved inorganic radiocarbon Prognostic tracer
     integer                                    :: id_dic_depth        = -1     ! Depth integral of dissolved inorganic carbon Prognostic tracer
     integer                                    :: id_fed_depth        = -1     ! Depth integral of dissolved Iron Prognostic tracer
     integer                                    :: id_o2_depth         = -1     ! Depth integral of oxygen Prognostic tracer
     integer                                    :: id_po4_depth        = -1     ! Depth integral of phosphate Prognostic tracer

     integer                                    :: id_fed_data_surf    = -1     ! Surface dissolved Iron data
     integer                                    :: id_htotal_surf      = -1     ! Surface hydrogen ion Diagnostic tracer
     integer                                    :: id_chl_surf         = -1     ! Surface chlorophyll Diagnostic tracer
     integer                                    :: id_biomass_p_surf   = -1     ! Surface biomass Diagnostic tracer
     integer                                    :: id_phyto_lg_surf    = -1     ! Surface large phytoplankton
     integer                                    :: id_phyto_sm_surf    = -1     ! Surface small phytoplankton
     integer                                    :: id_irr_mem_surf     = -1     ! Surface irradiance Memory Diagnostic tracer
     integer                                    :: id_fed_data_depth   = -1     ! Depth integral of dissolved Iron data
     integer                                    :: id_chl_depth        = -1     ! Depth integral of chlorophyll Diagnostic tracer
     integer                                    :: id_biomass_p_depth  = -1     ! Depth integral of biomass Diagnostic tracer
     integer                                    :: id_phyto_lg_depth   = -1     ! Depth integral of large phytoplankton
     integer                                    :: id_phyto_sm_depth   = -1     ! Depth integral of small phytoplankton
     integer                                    :: id_irr_mem_depth    = -1     ! Depth integral of irradiance Memory Diagnostic tracer

     logical                                    :: override_surf_temp = .true. ! True if overriding surface properties
     logical                                    :: override_surf_salt = .true. !   Must be true for first try, and will then
     logical                                    :: override_surf_alk  = .true. !   be set accordingly by data_override
     logical                                    :: override_surf_po4  = .true.
     logical                                    :: override_surf_sio4 = .true.
     logical                                    :: override_surf_dic  = .true.

  end type generic_miniBLING_type

  !An auxiliary type for storing varible names and descriptions
  type, public :: vardesc
     character(len=fm_string_len)       :: name     ! The variable name in a NetCDF file.
     character(len=fm_string_len)       :: longname ! The long name of that variable.
     character(len=1)                   :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
     character(len=1)                   :: z_grid   ! The vert. grid:  L, i, or 1.
     character(len=1)                   :: t_grid   ! The time description: s, a, m, or 1.
     character(len=fm_string_len)       :: units    ! The dimensions of the variable.
     character(len=1)                   :: mem_size ! The size in memory: d or f.
  end type vardesc

  type(generic_miniBLING_type), save       :: bling
  integer, parameter                 :: num_instances = 1
  !type(generic_miniBLING_type), dimension(:), pointer       :: bling
  !integer                                               :: num_instances

  type(CO2_dope_vector) :: CO2_dope_vec

contains


!#######################################################################

  subroutine generic_miniBLING_register(tracer_list)

    type(g_tracer_type), pointer, intent(inout) :: tracer_list

!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!
    character(len=fm_string_len), parameter :: sub_name = 'generic_miniBLING_register'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter   :: warn_header =                                &
         '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

    integer                                             :: n
    integer                                             :: stdout_unit

    stdout_unit = stdout()

    !Add here only the parameters that are required at the time of registeration 
    !(to make flux exchanging Ocean tracers known for all PE's) 
    !
    call g_tracer_start_param_list(package_name)

    call g_tracer_add_param('name', bling%name, '_')

    ! Turn on additional complexity. Most relevant diagnostic variables and all  
    ! tracers are not activated unless the appropriate switch is set to true.
    call g_tracer_add_param('do_14c',                  bling%do_14c,                  .true.)
    call g_tracer_add_param('do_carbon',               bling%do_carbon,               .true.)

    call g_tracer_add_param('ice_restart_file'  , bling%ice_restart_file  , 'ice_minibling.res.nc')
    call g_tracer_add_param('ocean_restart_file', bling%ocean_restart_file, 'ocean_minibling.res.nc')
    call g_tracer_add_param('IC_file'           , bling%IC_file       , '')

    call g_tracer_end_param_list(package_name)

    !-----------------------------------------------------------------------

!       Set the suffixes for this instance

      if (bling%name(1:1) .eq. '_') then  
        bling%suffix = ' '
        bling%long_suffix = ' '
      else  !}{
        bling%suffix = '_' // bling%name
        bling%long_suffix = ' (' // trim(bling%name) // ')'
      endif  !}
 
    !  Check for some possible fatal problems in the namelist variables.

      if ((bling%do_14c) .and. (bling%do_carbon)) then
        write (stdout_unit,*) trim(note_header), 'Simulating radiocarbon for instance ' // trim(bling%name)
      else if ((bling%do_14c) .and. .not. (bling%do_carbon)) then
        call mpp_error(FATAL, trim(error_header) //        &
             ' Do_14c requires do_carbon for instance ' // trim(bling%name))
      endif

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file    = bling%ice_restart_file,&
         ocean_restart_file  = bling%ocean_restart_file )

    do n = 1, num_instances  

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
      if (bling%fe_is_prognostic) then
        call g_tracer_add(tracer_list, package_name,                                            &
             name       = 'fed' // bling%suffix,                                          &
             longname   = 'Dissolved Iron' // bling%long_suffix,                          &
             units      = 'mol/kg',                                                             &
             prog       = .true.,                                                               &
             flux_runoff    = .false.,                                                          &
             flux_wetdep    = .true.,                                                           &
             flux_drydep    = .true.,                                                           &
             flux_param     = (/ 55.847e-03 /),                                                 &
             flux_bottom    = .true. )
      elseif (bling%fe_is_diagnostic) then
        call g_tracer_add(tracer_list, package_name,                                            &
             name       = 'fed' // bling%suffix,                                          &
             longname   = 'Dissolved Iron' // bling%long_suffix,                          &
             units      = 'mol/kg',                                                             &
             prog       = .false.)
      else
        call mpp_error(NOTE, trim(note_header) // ' Fe is data overridden for instance ' // trim(bling%name))
      endif

    !       O2
    !
      call g_tracer_add(tracer_list, package_name,                                              &
           name       = 'o2' // bling%suffix,                                             &
           longname   = 'Oxygen' // bling%long_suffix,                                    &
           units      = 'mol/kg',                                                               &
           prog       = .true.,                                                                 &
           flux_gas       = .true.,                                                             &
           flux_gas_type  = 'air_sea_gas_flux_generic',                                          &
           flux_gas_name  = 'o2_flux' // trim(bling%suffix),                              &
           flux_gas_molwt = WTMO2,                                                              &
           flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                                         &
           flux_bottom    = .true.,                                                             &
           flux_gas_restart_file  = 'ocean_minibling_airsea_flux.res.nc' )

   
    !       PO4
    !
      call g_tracer_add(tracer_list, package_name,                                              &
           name       = 'po4' // bling%suffix,                                            &
           longname   = 'Phosphate' // bling%long_suffix,                                 &
           units      = 'mol/kg',                                                               &
           prog       = .true.,                                                                 &
           flux_bottom    = .true.     )

    
    !===========================================================
    !Diagnostic Tracers
    !===========================================================

    !       Chl (Chlorophyll)
    !
      call g_tracer_add(tracer_list, package_name,                                              &
           name       = 'chl' // bling%suffix,                                            &
           longname   = 'Chlorophyll' // bling%long_suffix,                               &
           units      = 'ug kg-1',                                                              &
           prog       = .false.,                                                                &
           init_value = 0.08          )

    !       Irr_mem (Irradiance Memory)
    !
      call g_tracer_add(tracer_list, package_name,                                              &
           name       = 'irr_mem' // bling%suffix,                                        &
           longname   = 'Irradiance memory' // bling%long_suffix,                         &
           units      = 'Watts/m^2',                                                            &
           prog       = .false.)

      if (bling%biomass_type .eq. 'single') then

      !       Biomass
      !
        call g_tracer_add(tracer_list, package_name,                                              &
             name       = 'biomass_p' // bling%suffix,                                      &
             longname   = 'Biomass in P units' // bling%long_suffix,                        &
             units      = 'mol P kg-1',                                                           &
             prog       = .false.)

      elseif (bling%biomass_type .eq. 'lg_sm_phyto') then

      !       Large phytoplankton biomass
      !
        call g_tracer_add(tracer_list, package_name,                                              &
             name       = 'phyto_lg' // bling%suffix,                                       &
             longname   = 'Large phytoplankton biomass in P units' // bling%long_suffix,    &
             units      = 'mol P kg-1',                                                           &
             prog       = .false.,                                                                &
             init_value = 4.e-07 )

      !       Small phytoplankton biomass
      !
        call g_tracer_add(tracer_list, package_name,                                              &
             name       = 'phyto_sm' // bling%suffix,                                       &
             longname   = 'Small phytoplankton biomass in P units' // bling%long_suffix,    &
             units      = 'mol P kg-1',                                                           &
             prog       = .false.,                                                                &
             init_value = 4.e-07 )

      else

        call mpp_error(FATAL, trim(error_header) // ' Unknown biomass type "' // trim(bling%biomass_type) // '"')

      endif

      if (bling%do_carbon) then                                    !<<CARBON CYCLE

    !       DIC (Dissolved inorganic carbon)
    !
        call g_tracer_add(tracer_list, package_name,                                            &
             name       = 'dic' // bling%suffix,                                          &
             longname   = 'Dissolved Inorganic Carbon' // bling%long_suffix,              &
             units      = 'mol/kg',                                                             &
             prog       = .true.,                                                               &
             flux_gas       = .true.,                                                           &
             flux_gas_type  = 'air_sea_gas_flux_generic',                                        &
             flux_gas_name  = 'co2_flux' // trim(bling%suffix),                           &
             flux_gas_molwt = WTMCO2,                                                           &
             flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                                       &
             flux_gas_restart_file  = 'ocean_minibling_airsea_flux.res.nc',  &
             flux_runoff    = .false.,                                                          &
             flux_param     = (/ 12.011e-03 /),                                                 &
             flux_bottom    = .true.,                                                           &
             init_value     = 0.001)

    !Diagnostic Tracers:

    !       Htotal (H+ ion concentration)
    !
        call g_tracer_add(tracer_list, package_name,                                            &
             name       = 'htotal' // bling%suffix,                                       &
             longname   = 'H+ ion concentration' // bling%long_suffix,                    &
             units      = 'mol/kg',                                                             &
             prog       = .false.,                                                              &
             init_value = bling%htotal_in)


        if (bling%do_14c) then                                      !<<RADIOCARBON
      !       DI14C (Dissolved inorganic radiocarbon)
      !
          call g_tracer_add(tracer_list, package_name,                                          &
             name       = 'di14c' // bling%suffix,                                        &
             longname   = 'Dissolved Inorganic Radiocarbon' // bling%long_suffix,         &
             units      = 'mol/kg',                                                             &
             prog       = .true.,                                                               &
             flux_gas       = .true.,                                                           &
             flux_gas_type  = 'air_sea_gas_flux_generic',                                        &
             flux_gas_name  = 'c14o2_flux' // trim(bling%suffix),                         &
             flux_gas_molwt = WTMCO2,                                                           &
             flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                                       &
             flux_gas_restart_file  = 'ocean_minibling_airsea_flux.res.nc',  &
             flux_param     = (/ 14.e-03 /),                                                    &
             flux_bottom    = .true.,                                                           &
             init_value     = 0.001)

        endif  !}                                               !RADIOCARBON>>

      endif  !}                                                !CARBON CYCLE>>

    enddo  !} n
    
  end subroutine generic_miniBLING_register


!#######################################################################
  ! <SUBROUTINE NAME="generic_miniBLING_init">
  !  <OVERVIEW>
  !   Initialize the generic miniBLING module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the miniBLING Tracers to the list of generic Tracers passed 
  !  to it via utility subroutine g_tracer_add(). Adds all the parameters 
  !  used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_miniBLING_init
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_miniBLING_init(tracer_list)

    type(g_tracer_type), pointer :: tracer_list

!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!
    character(len=64), parameter                              :: sub_name = 'generic_miniBLING_init'

    character(len=256)                                        :: caller_str
    character(len=256)                                        :: error_header
    character(len=256)                                        :: warn_header
    character(len=256)                                        :: note_header
    integer                                                   :: n
    character(len=fm_field_name_len)                          :: name
    integer                                                   :: nn
    character(len=fm_field_name_len), pointer, dimension(:)   :: names => NULL()
    integer                                                   :: stdout_unit
    character(len=fm_string_len)                              :: string
    integer                                                   :: package_index

    stdout_unit = stdout()
    
    !       Set up the field input

    caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[]'
    error_header = '==>Error from '   // trim(caller_str) // ':'
    warn_header =  '==>Warning from ' // trim(caller_str) // ':'
    note_header =  '==>Note from '    // trim(caller_str) // ':'

    write (stdout_unit,*)
    write (stdout_unit,*) trim(note_header), ' Processing generic tracer package miniBLING'

    do n = 1, num_instances  

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_miniBLING_params type
    !==============================================================    

    !Add the known experimental parameters used for calculations in this module.
    !All the g_tracer_add_param calls must happen between 
    !g_tracer_start_param_list and g_tracer_end_param_list  calls.
    !This implementation enables runtime overwrite via field_table.

      call g_tracer_start_param_list(package_name)

    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
      call g_tracer_add_param('RHO_0', bling%Rho_0, 1035.0)
 
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
      call g_tracer_add_param('alk_scheme', bling%alk_scheme, 'normal')
      call g_tracer_add_param('alk_slope', bling%alk_slope, 32.0e-06)
      call g_tracer_add_param('alk_intercept', bling%alk_intercept, 1200.0e-06)
      call g_tracer_add_param('biomass_type', bling%biomass_type, 'single')
      call g_tracer_add_param('alpha_photo', bling%alpha_photo, 1.e-5 * 2.77e+18 / 6.022e+17) ! g C g Chl-1 m2 W-1 s-1
      call g_tracer_add_param('kappa_eppley', bling%kappa_eppley, 0.063)              ! deg C-1
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

      call g_tracer_add_param('po4_min', bling%po4_min,  1.0e-8)                       ! mol PO4 kg-1

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

    !
    !   New paramter to help with non-prognostic iron
    !

      call g_tracer_add_param('def_fe_min', bling%def_fe_min, 0.0)   ! ?

    !
    ! If fe_is_prognostic is true, then Fed will be a prognostic variable, otherwise
    ! if fe_is_diagnostic is true, then it will be diagnostic, restoring to a 3-d field with
    ! a time-scale of fe_restoring (in days), otherwise fed will be data driven
    ! with any coastal increase.
    !
      call g_tracer_add_param('fe_is_prognostic', bling%fe_is_prognostic, .false.)
      call g_tracer_add_param('fe_is_diagnostic', bling%fe_is_diagnostic, .false.)
      call g_tracer_add_param('fe_restoring',     bling%fe_restoring,        10.0)  ! days
      call g_tracer_add_param('fe_coastal',       bling%fe_coastal,       2.0e-09)  ! mol/kg
      call g_tracer_add_param('fe_coastal_depth', bling%fe_coastal_depth,   200.0)  ! m

    ! The k_fe_2_p is the Fe:P at which the iron-limitation term has a 
    ! value of 0.5, chosen according to Sunda and Huntsman (Fig. 2, 
    ! Nature, 1997). Converted from Fe:C ratio.

      call g_tracer_add_param('k_fe_2_p', bling%k_fe_2_p,  7.e-6 * 106.)           ! mol Fe mol P-1

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

    ! Half saturation constant for fast recycling of P, very low to act only in nutrient-poor waters
      call g_tracer_add_param('k_po4_recycle', bling%k_po4_recycle,  2.0e-8)       ! mol PO4 kg-1

    !-----------------------------------------------------------------------
    ! Remineralization
    !-----------------------------------------------------------------------
    !
    ! Stoichiometric ratios taken from Anderson (1995) as discussed in
    ! Sarmiento and Gruber (2008), and Sarmiento et al. (2002) for Ca:P.
    !
      call g_tracer_add_param('c_2_p', bling%c_2_p, 106.0 )                        ! mol C mol P-1
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
    ! fully oxidized conditions. As this code is currently intended for short, high-resolution runs,
    ! this value is set to zero to cause a cessation of remineralization under suboxia/anoxia.
    ! This will allow P to sink past the OMZ, which lead lead to a downward expansion of the OMZ,
    ! but it hopefully won't be a huge problem on the timescale of 100-200 years.
    !
      call g_tracer_add_param('remin_min', bling%remin_min, 0.0)                   ! dimensionless
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
    ! Mimimum fraction of POP (for turning off recycling set to 1.0)
    !
      call g_tracer_add_param('min_frac_pop', bling%min_frac_pop, 0.0)
    !
    ! Depth for integral and flux diagnostics
    !
    
    !-----------------------------------------------------------------------
    ! Miscellaneous
    !-----------------------------------------------------------------------
    !
      call g_tracer_add_param('diag_depth', bling%diag_depth, 100.0)    ! use nearest integer

    !
      call g_tracer_end_param_list(package_name)

    !
    ! Check the diag depth and set a string for that depth
    !

      if (bling%diag_depth .gt. 0.0) then
        bling%diag_depth = nint(bling%diag_depth)
        write (bling%diag_depth_str, '(f10.0)') bling%diag_depth
        bling%diag_depth_str = adjustl(bling%diag_depth_str)
        bling%diag_depth_str = bling%diag_depth_str(1:len_trim(bling%diag_depth_str)-1)   ! remove trailing decimal point
      else
        call mpp_error(FATAL, trim(error_header) // ' diag_depth <= 0 for instance ' // trim(bling%name))
      endif

    enddo  !} n
    
    !   Allocate all the private work arrays used by this module.
    call user_allocate_arrays

  end subroutine generic_miniBLING_init



!#######################################################################
  !   Register diagnostic fields to be used in this module. 
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !

  subroutine generic_miniBLING_register_diag

    real,          parameter    :: missing_value1 = -1.0e+10
    type(vardesc)               :: vardesc_temp
    integer                     :: isc
    integer                     :: iec
    integer                     :: jsc
    integer                     :: jec
    integer                     :: isd
    integer                     :: ied
    integer                     :: jsd
    integer                     :: jed
    integer                     :: nk
    integer                     :: ntau
    integer                     :: n
    integer                     :: axes(3)
    type(time_type)             :: init_time 

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, axes = axes, init_time = init_time) 


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

    do n = 1, num_instances  

      if (bling%fe_is_prognostic) then
        vardesc_temp = vardesc&
        ("b_fed","Bottom flux of Fe into sediment",'h','1','s','mol m-2 s-1','f')
        bling%id_b_fed = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("b_o2","Bottom flux of O2 into sediment",'h','1','s','mol m-2 s-1','f')
      bling%id_b_o2 = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("b_po4","Bottom flux of PO4 into sediment",'h','1','s','mol m-2 s-1','f')
      bling%id_b_po4 = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                   &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      if (bling%biomass_type .eq. 'single') then
        vardesc_temp = vardesc&
        ("biomass_p_ts","Instantaneous P concentration in biomass",'h','L','s','mol kg-1','f')
        bling%id_biomass_p_ts = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,            &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("def_Fe","Iron deficiency term",'h','L','s','unitless','f')
      bling%id_def_fe = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                  &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("expkT","Temperature dependence",'h','L','s','unitless','f')
      bling%id_expkT = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                   &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("fe_2_p_uptake","Uptake ratio of Fed:PO4",'h','L','s','mol Fe mol P-1','f')
      bling%id_fe_2_p_uptake = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,         &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
           vardesc_temp%units, missing_value = missing_value1)
      if (bling%fe_is_prognostic) then
        vardesc_temp = vardesc&
        ("fe_burial","Sedimenting iron flux",'h','1','s','mol m-2 s-1','f')
        bling%id_fe_burial = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,               &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("feprime","Concentration of free, unbound iron",'h','L','s','mol kg-1','f')
        bling%id_feprime = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("ffe_sed","Sediment iron efflux",'h','1','s','mol m-2 s-1','f')
        bling%id_ffe_sed = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("fpofe","POFe sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
        bling%id_fpofe = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                   &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("fpop_" // trim(bling%diag_depth_str),"POP sinking flux at " // trim(bling%diag_depth_str) // " m",              &
           'h','L','s','mol m-2 s-1','f')
      bling%id_fpop_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("fpop","POP sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
      bling%id_fpop = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("frac_lg","Fraction of production by large phytoplankton",'h','L','s','unitless','f')
      bling%id_frac_lg = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("frac_pop","Particulate fraction of total uptake",'h','L','s','unitless','f')
      bling%id_frac_pop = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("irr_inst","Instantaneous light",'h','L','s','W m-2','f')
      bling%id_irr_inst = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("irr_mix","Mixed layer light",'h','L','s','W m-2','f')
      bling%id_irr_mix = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("irrk","Tendency to light limitation",'h','L','s','W m-2','f')
      bling%id_irrk = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      if (bling%fe_is_prognostic) then
        vardesc_temp = vardesc&
        ("jfe_ads_inorg","Iron adsorption (2nd order) layer integral",'h','L','s','mol m-2 s-1','f')
        bling%id_jfe_ads_inorg = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,           &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("jfe_ads_org","Iron adsorption to FPOP layer integral",'h','L','s','mol m-2 s-1','f')
        bling%id_jfe_ads_org = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("jfe_recycle","Fast recycling of iron layer integral",'h','L','s','mol m-2 s-1','f')
        bling%id_jfe_recycle = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      if (bling%fe_is_prognostic .or. bling%fe_is_diagnostic) then
        vardesc_temp = vardesc&
        ("jfe_reminp","Sinking particulate Fe decay layer integral",'h','L','s','mol m-2 s-1','f')
        bling%id_jfe_reminp = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("jfe_uptake","Iron production layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jfe_uptake = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jo2_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral of o2 source",           &
           'h','L','s','mol m-2 s-1','f')
      bling%id_jo2_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                     &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jo2","O2 source layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jo2 = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                     &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jp_recycle_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral of fast recycling of PO4",                &
           'h','L','s','mol m-2 s-1','f')
      bling%id_jp_recycle_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jp_recycle","Fast recycling of PO4 layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jp_recycle = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jp_reminp_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral of sinking particulate P decay",           &
           'h','L','s','mol m-2 s-1','f')
      bling%id_jp_reminp_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,               &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jp_reminp","Sinking particulate P decay layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jp_reminp = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,               &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jp_uptake_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral of pO4 uptake",            &
           'h','L','s','mol m-2 s-1','f')
      bling%id_jp_uptake_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,               &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jp_uptake","PO4 uptake layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jp_uptake = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,               &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jpo4_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral of PO4 source",         &
           'h','L','s','mol m-2 s-1','f')
      bling%id_jpo4_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                     &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jpo4","PO4 source layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jpo4 = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("jpop","Particulate P source layer integral",'h','L','s','mol m-2 s-1','f')
      bling%id_jpop = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      if (bling%fe_is_prognostic) then
        vardesc_temp = vardesc&
        ("jfeop","Particulate Fe source layer integral",'h','L','s','mol m-2 s-1','f')
        bling%id_jfeop = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                   &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("kfe_eq_lig","Iron ligand stability constant",'h','L','s','mol-1 kg','f')
        bling%id_kfe_eq_lig = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("mu","Net growth rate after respiratory loss",'h','L','s','s-1','f')
      bling%id_mu = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                      &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("o2_saturation","Saturation O2 concentration",'h','1','s','mol kg-1','f')
      bling%id_o2_saturation = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,           &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("pc_m","Light-saturated photosynthesis rate (carbon specific)",'h','L','s','s-1','f')
      bling%id_pc_m = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                    &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("theta","Chl:C ratio",'h','L','s','g Chl g C-1','f')
      bling%id_theta = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                   &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("thetamax_fe","Fe-limited max Chl:C",'h','L','s','g Chl g C-1','f')
      bling%id_thetamax_fe = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("wsink","Sinking rate",'h','L','s','m s-1','f')
      bling%id_wsink = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                   &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("zremin","Remineralization lengthscale",'h','L','s','m','f')
      bling%id_zremin = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                  &
           axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("o2_surf","Surface O2 concentration",'h','1','s','mol kg-1','f')
      bling%id_o2_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("o2_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral O2",              &
           'h','1','s','mol m-2','f')
      bling%id_o2_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("fed_surf","Surface Fed concentration",'h','1','s','mol kg-1','f')
      bling%id_fed_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("fed_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral Fed",            &
           'h','1','s','mol m-2','f')
      bling%id_fed_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      if (.not. bling%fe_is_prognostic) then
        vardesc_temp = vardesc&
        ("fed_data_surf","Surface Fed data concentration",'h','1','s','mol kg-1','f')
        bling%id_fed_data_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("fed_data_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral Fedconcentration",                &
             'h','1','s','mol m-2','f')
        bling%id_fed_data_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("po4_surf","Surface PO4 concentration",'h','1','s','mol kg-1','f')
      bling%id_po4_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("po4_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral PO4",            &
           'h','1','s','mol m-2','f')
      bling%id_po4_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("htotal_surf","Surface H+ concentration",'h','1','s','mol kg-1','f')
      bling%id_htotal_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      if (bling%biomass_type .eq. 'single') then
        vardesc_temp = vardesc&
        ("biomass_p_surf","Surface Biomass-P concentration",'h','1','s','mol kg-1','f')
        bling%id_biomass_p_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("biomass_p_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral BiomasP concentration",          &
             'h','1','s','mol m-2','f')
        bling%id_biomass_p_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      elseif (bling%biomass_type .eq. 'lg_sm_phyto') then
        vardesc_temp = vardesc&
        ("phyto_lg_surf","Surface large phytoplankton concentration",'h','1','s','mol kg-1','f')
        bling%id_phyto_lg_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("phyto_lg_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral largeconcentration",              &
             'h','1','s','mol m-2','f')
        bling%id_phyto_lg_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("phyto_sm_surf","Surface small phytoplankton concentration",'h','1','s','mol kg-1','f')
        bling%id_phyto_sm_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("phyto_sm_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral smallconcentration",              &
             'h','1','s','mol m-2','f')
        bling%id_phyto_sm_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif
      vardesc_temp = vardesc&
      ("chl_surf","Surface Chl concentration",'h','1','s','mol kg-1','f')
      bling%id_chl_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("chl_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral Chl",            &
           'h','1','s','mol m-2','f')
      bling%id_chl_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("irr_mem_surf","Surface IRR_mem concentration",'h','1','s','mol kg-1','f')
      bling%id_irr_mem_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,            &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc&
      ("irr_mem_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral IRR_mem",            &
           'h','1','s','mol m-2','f')
      bling%id_irr_mem_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,            &
           axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
           vardesc_temp%units, missing_value = missing_value1)
      if (.not. bling%fe_is_prognostic) then
        vardesc_temp = vardesc&
        ("fed_data","Fed data concentration",'h','1','s','mol kg-1','f')
        bling%id_fed_data = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                &
             axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                        &
             vardesc_temp%units, missing_value = missing_value1)
      endif

 
      if (bling%do_carbon) then                                    !<<CARBON CYCLE
        vardesc_temp = vardesc&
        ("dic_surf","Surface DIC concentration",'h','1','s','mol kg-1','f')
        bling%id_dic_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("dic_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral DIC",          &
             'h','1','s','mol m-2','f')
        bling%id_dic_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("b_dic","Bottom flux of DIC into sediment",'h','1','s','mol m-2 s-1','f')
        bling%id_b_dic = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,                 &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("co2_alpha","Saturation surface CO2* per uatm",'h','1','s','mol kg-1 atm-1','f')
        bling%id_co2_alpha = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("co2_csurf","CO2* concentration at surface",'h','1','s','mol kg-1','f')
        bling%id_co2_csurf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("pco2_surf","Seawater pCO2 in surface layer",'h','1','s','uatm','f')
        bling%id_pco2_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("temp_co2calc","Surface temperature used for co2calc",'h','1','s','deg C','f')
        bling%id_temp_co2calc = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,     &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("salt_co2calc","Surface salinity used for co2calc",'h','1','s','PSU','f')
        bling%id_salt_co2calc = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,     &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("alk_co2calc","Surface alkalinity used for co2calc",'h','1','s','eq kg-1','f')
        bling%id_alk_co2calc = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,      &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("po4_co2calc","Surface phosphate used for co2calc",'h','1','s','mol kg -1','f')
        bling%id_po4_co2calc = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,      &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("sio4_co2calc","Surface silicate used for co2calc",'h','1','s','mol kg -1','f')
        bling%id_sio4_co2calc = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,     &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc&
        ("dic_co2calc","Surface DIC used for co2calc",'h','1','s','mol kg -1','f')
        bling%id_dic_co2calc = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,      &
             axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                      &
             vardesc_temp%units, missing_value = missing_value1)

        if (bling%do_14c) then                                      !<<RADIOCARBON
          vardesc_temp = vardesc&
          ("di14c_surf","Surface DI14C concentration",'h','1','s','mol kg-1','f')
          bling%id_di14c_surf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
               axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("di14c_" // trim(bling%diag_depth_str),trim(bling%diag_depth_str) // " m integral DI14C",            &
               'h','1','s','mol m-2','f')
          bling%id_di14c_depth = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,          &
               axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("b_di14c","Bottom flux of DI14C into sediment",'h','1','s','mol m-2 s-1','f')
          bling%id_b_di14c = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
               axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("c14_2_p","Ratio of DI14C to PO4",'h','L','s','mol kg-1','f')
          bling%id_c14_2_p = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,             &
               axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("c14o2_alpha","Saturation surface 14CO2* per uatm",'h','1','s','mol kg-1 atm-1','f')
          bling%id_c14o2_alpha = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,         &
               axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("c14o2_csurf","14CO2* concentration at surface",'h','1','s','mol kg-1','f')
          bling%id_c14o2_csurf = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,         &
               axes(1:2), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("fpo14c","PO14C sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
          bling%id_fpo14c = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
               axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("j14c_decay_dic","DI14C radioactive decay layer integral",'h','L','s','mol m-2 s-1','f')
          bling%id_j14c_decay_dic = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,      &
               axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("j14c_reminp","Sinking PO14C remineralization layer integral",'h','L','s','mol m-2 s-1','f')
          bling%id_j14c_reminp = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,         &
               axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          vardesc_temp = vardesc&
          ("jdi14c","DI14C source layer integral",'h','L','s','mol m-2 s-1','f')
          bling%id_jdi14c = register_diag_field(package_name, trim(vardesc_temp%name) // bling%suffix,              &
               axes(1:3), init_time, trim(vardesc_temp%longname) // bling%long_suffix,                                    &
               vardesc_temp%units, missing_value = missing_value1)
          
        endif  !}                                               !RADIOCARBON>>

      endif  !}                                                !CARBON CYCLE>>

    enddo  !} n

  end subroutine generic_miniBLING_register_diag


!#######################################################################
! <SUBROUTINE NAME="generic_miniBLING_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracer fields could be modified after values are obtained from the 
  !  coupler. This subroutine is the place for specific tracer manipulations.
  !  miniBLING currently does not use this.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_miniBLING_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_miniBLING_update_from_coupler(tracer_list)

    type(g_tracer_type), pointer, intent(inout) :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_miniBLING_update_from_coupler'

  end subroutine generic_miniBLING_update_from_coupler


!#######################################################################
  ! <SUBROUTINE NAME="generic_miniBLING_update_from_bottom">
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracers could have bottom fluxes and reservoirs. 
  !   This subroutine is the place for specific tracer manipulations.
  !   miniBLING currently does not use this.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_miniBLING_update_from_bottom(tracer_list,dt, tau) 
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

  subroutine generic_miniBLING_update_from_bottom(tracer_list, dt, tau)

    type(g_tracer_type), pointer, intent(inout) :: tracer_list
    real,                         intent(in)    :: dt
    integer,                      intent(in)    :: tau

  end subroutine generic_miniBLING_update_from_bottom
  

!#######################################################################
  ! <SUBROUTINE NAME="generic_miniBLING_diag">
  !  <OVERVIEW>
  !   Do things which must be done after tronsports and sources have been applied
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine saves out surface diagnostic firlds for prognostic tracers
  !   after vertical transport has been calculated
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_miniBLING_diag(tracer_list,tau,model_time) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !  <IN NAME="model_time" TYPE="time_type">
  !   Model time
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_miniBLING_diag(tracer_list, ilb, jlb, tau, model_time, dzt, rho_dzt, caller)

    type(g_tracer_type),            pointer, intent(inout)              :: tracer_list
    integer,                                 intent(in)                 :: ilb
    integer,                                 intent(in)                 :: jlb
    integer,                                 intent(in)                 :: tau
    type(time_type),                         intent(in)                 :: model_time
    real, dimension(ilb:,jlb:,:),            intent(in)                 :: dzt
    real, dimension(ilb:,jlb:,:),            intent(in)                 :: rho_dzt
    character(len=*),                        intent(in),    optional    :: caller

!-----------------------------------------------------------------------
!     local parameters

    character(len=fm_string_len), parameter     :: sub_name = 'generic_miniBLING_diag'

    character(len=256)                          :: caller_str
    character(len=256)                          :: error_header
    character(len=256)                          :: warn_header
    character(len=256)                          :: note_header
    integer                                     :: isc
    integer                                     :: iec
    integer                                     :: jsc
    integer                                     :: jec
    integer                                     :: isd
    integer                                     :: ied
    integer                                     :: jsd
    integer                                     :: jed
    integer                                     :: nk
    integer                                     :: ntau
    integer                                     :: i
    integer                                     :: j
    integer                                     :: k 
    integer                                     :: n
    real,    dimension(:,:,:), pointer          :: grid_tmask
    logical                                     :: used
    integer                                     :: k_int
    logical                                     :: diag_initialized
    
    !  Set up the headers for stdout messages.

    if (present(caller)) then
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[' // trim(caller) // ']'
    else
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[]'
    endif
    error_header = '==> Error from '   // trim(caller_str) // ':'
    warn_header =  '==> Warning from ' // trim(caller_str) // ':'
    note_header =  '==> Note from '    // trim(caller_str) // ':'

    !  Set up the module if not already done

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau,  &
         grid_tmask = grid_tmask)

    !
    !-----------------------------------------------------------------------
    !       Save depth integrals and fluxes
    !-----------------------------------------------------------------------
    !

    k_int = 0
    diag_initialized = .false.

    do n = 1, num_instances

      if (bling%id_po4_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_po4(:,:,:,tau), dzt, rho_dzt,      &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_po4_depth, bling%integral,                                            &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_o2_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_o2(:,:,:,tau), dzt, rho_dzt,       &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_o2_depth, bling%integral,                                             &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_dic_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_dic(:,:,:,tau), dzt, rho_dzt,      &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_dic_depth, bling%integral,                                            &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_di14c_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_di14c(:,:,:,tau), dzt, rho_dzt,    &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_di14c_depth, bling%integral,                                          &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%fe_is_prognostic) then
        if (bling%id_fed_depth .gt. 0) then
          call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_fed(:,:,:,tau), dzt, rho_dzt,    &
               bling%wrk_3d, k_int, bling%integral)
          used = send_data(bling%id_fed_depth, bling%integral,                                          &
               model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        endif  !}
      elseif (bling%fe_is_diagnostic) then
        if (bling%id_fed_depth .gt. 0) then
          call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_fed_diag, dzt, rho_dzt,          &
               bling%wrk_3d, k_int, bling%integral)
          used = send_data(bling%id_fed_depth, bling%integral,                                          &
               model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        endif  !}
      else
        if (bling%id_fed_depth .gt. 0) then
          call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%f_fed, dzt, rho_dzt,               &
               bling%wrk_3d, k_int, bling%integral)
          used = send_data(bling%id_fed_depth, bling%integral,                                          &
               model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        endif  !}
      endif
      if (.not. bling%fe_is_prognostic) then
        if (bling%id_fed_data_depth .gt. 0) then
          call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%f_fed_data, dzt, rho_dzt,          &
               bling%wrk_3d, k_int, bling%integral)
          used = send_data(bling%id_fed_data_depth, bling%integral,                                     &
               model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        endif  !}
      endif  !}
      if (bling%id_chl_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%f_chl, dzt, rho_dzt,                 &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_chl_depth, bling%integral,                                            &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_biomass_p_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_biomass_p, dzt, rho_dzt,           &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_biomass_p_depth, bling%integral,                                      &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_phyto_lg_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_phyto_lg, dzt, rho_dzt,            &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_phyto_lg_depth, bling%integral,                                       &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_phyto_sm_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_phyto_sm, dzt, rho_dzt,            &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_phyto_sm_depth, bling%integral,                                       &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_irr_mem_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%p_irr_mem, dzt, rho_dzt,             &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_irr_mem_depth, bling%integral,                                        &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_jp_uptake_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%jp_uptake, dzt, rho_dzt,             &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_jp_uptake_depth, bling%integral,                                      &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_jp_recycle_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%jp_recycle, dzt, rho_dzt,            &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_jp_recycle_depth, bling%integral,                                     &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_jp_reminp_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%jp_reminp, dzt, rho_dzt,             &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_jp_reminp_depth, bling%integral,                                      &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_jpo4_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%jpo4, dzt, rho_dzt,                  &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_jpo4_depth, bling%integral,                                           &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_jo2_depth .gt. 0) then
        call g_tracer_column_int(bling%diag_depth, isd, jsd, bling%jo2, dzt, rho_dzt,                   &
             bling%wrk_3d, k_int, bling%integral)
        used = send_data(bling%id_jo2_depth, bling%integral,                                            &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}
      if (bling%id_fpop_depth .gt. 0) then
        call g_tracer_flux_at_depth(bling%diag_depth, isd, jsd, bling%fpop, dzt,                        &
             bling%k_lev, bling%wrk_2d, diag_initialized, bling%flux)
        used = send_data(bling%id_fpop_depth, bling%flux,                                               &
             model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif  !}

    enddo

    !
    !-----------------------------------------------------------------------
    !       Save surface prognostic variables for diagnostics, after vertical diffusion
    !-----------------------------------------------------------------------
    !

    do n = 1, num_instances

      if (bling%id_po4_surf .gt. 0)                                             &
           used = send_data(bling%id_po4_surf,    bling%p_po4(:,:,1,tau),       &
           model_time, rmask = grid_tmask(:,:,1),                               & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_o2_surf .gt. 0)                                              &
           used = send_data(bling%id_o2_surf,    bling%p_o2(:,:,1,tau),         &
           model_time, rmask = grid_tmask(:,:,1),                               & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_dic_surf .gt. 0)                                             &
           used = send_data(bling%id_dic_surf,    bling%p_dic(:,:,1,tau),       &
           model_time, rmask = grid_tmask(:,:,1),                               & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_di14c_surf .gt. 0)                                           &
           used = send_data(bling%id_di14c_surf,    bling%p_di14c(:,:,1,tau),   &
           model_time, rmask = grid_tmask(:,:,1),                               & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%fe_is_prognostic) then
        if (bling%id_fed_surf .gt. 0)                                           &
             used = send_data(bling%id_fed_surf,    bling%p_fed(:,:,1,tau),     &
             model_time, rmask = grid_tmask(:,:,1),                             & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      elseif (bling%fe_is_diagnostic) then
        if (bling%id_fed_surf .gt. 0)                                           &
             used = send_data(bling%id_fed_surf,    bling%p_fed_diag(:,:,1),    &
             model_time, rmask = grid_tmask(:,:,1),                             & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      else
        if (bling%id_fed_surf .gt. 0)                                           &
             used = send_data(bling%id_fed_surf,    bling%f_fed(:,:,1),         &
             model_time, rmask = grid_tmask(:,:,1),                             & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      endif

    enddo

    return

  end subroutine generic_miniBLING_diag


!#######################################################################
  ! <SUBROUTINE NAME="generic_miniBLING_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine contains most of the biogeochemistry for calculating the 
  !   interaction of the core set of tracers with each other and with outside forcings.
  !   Additional tracers (e.g. carbon, isotopes) are calculated in other subroutines.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_miniBLING_update_from_source(tracer_list,Temp,Salt,dzt,hblt_depth,&
  !                                         ilb,jlb,tau,dtts, grid_dat,sw_pen,opacity) 
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
  !  <IN NAME="dtts" TYPE="real">
  !   Time step increment
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_miniBLING_update_from_source(tracer_list, Temp, Salt,                  &
       rho_dzt, dzt, hblt_depth, ilb, jlb, tau, dtts, grid_dat, model_time, nbands,     &
       max_wavelength_band, sw_pen_band, opacity_band, grid_ht)

    type(g_tracer_type),            pointer, intent(inout)              :: tracer_list
    real, dimension(ilb:,jlb:,:),            intent(in)                 :: Temp
    real, dimension(ilb:,jlb:,:),            intent(in)                 :: Salt
    real, dimension(ilb:,jlb:,:),            intent(in)                 :: rho_dzt
    real, dimension(ilb:,jlb:,:),            intent(in)                 :: dzt
    real, dimension(ilb:,jlb:),              intent(in)                 :: hblt_depth
    real, dimension(ilb:,jlb:),              intent(in)                 :: grid_ht
    integer,                                 intent(in)                 :: ilb
    integer,                                 intent(in)                 :: jlb
    integer,                                 intent(in)                 :: tau
    real,                                    intent(in)                 :: dtts
    real, dimension(ilb:,jlb:),              intent(in)                 :: grid_dat
    type(time_type),                         intent(in)                 :: model_time
    integer,                                 intent(in)                 :: nbands
    real, dimension(:),                      intent(in)                 :: max_wavelength_band
    real, dimension(:,ilb:,jlb:),            intent(in)                 :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:),          intent(in)                 :: opacity_band

!-----------------------------------------------------------------------
!     local parameters

    character(len=fm_string_len), parameter     :: sub_name = 'generic_miniBLING_update_from_source'

    character(len=256)                          :: caller_str
    character(len=256)                          :: error_header
    character(len=256)                          :: warn_header
    character(len=256)                          :: note_header
    integer                                     :: isc
    integer                                     :: iec
    integer                                     :: jsc
    integer                                     :: jec
    integer                                     :: isd
    integer                                     :: ied
    integer                                     :: jsd
    integer                                     :: jed
    integer                                     :: nk
    integer                                     :: ntau
    integer                                     :: i
    integer                                     :: j
    integer                                     :: k 
    integer                                     :: kblt
    integer                                     :: n
    real,    dimension(:,:,:), pointer          :: grid_tmask
    integer, dimension(:,:),   pointer          :: grid_kmt
    logical                                     :: used
    integer                                     :: nb
    real                                        :: tmp_hblt
    real                                        :: tmp_Irrad
    real                                        :: tmp_irrad_ML
    real                                        :: tmp_phyto_lg_ML
    real                                        :: tmp_phyto_sm_ML
    real                                        :: tmp_opacity
    real,    dimension(:),     Allocatable      :: tmp_irr_band 
    real                                        :: s_over_p
    
    !  Set up the headers for stdout messages.

    caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[]'
    error_header = '==>Error from '   // trim(caller_str) // ':'
    warn_header =  '==>Warning from ' // trim(caller_str) // ':'
    note_header =  '==>Note from '    // trim(caller_str) // ':'

    !  Set up the module if not already done

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau,  &
         grid_tmask = grid_tmask, grid_kmt = grid_kmt)



  ! SURFACE GAS FLUXES
  ! 
  ! This subroutine coordinates the calculation of gas concentrations and solubilities 
  ! in the surface layer. The concentration of a gas is written as csurf, while the
  ! solubility (in mol kg-1 atm-1 or mol m-3 atm-1) is written as alpha. These two
  ! quantities are passed to the coupler, which multiplies their difference by the
  ! gas exchange piston velocity over the mixed layer depth to provide the gas
  ! exchange flux,
  !    Flux = Kw/dz * (alpha - csurf)
  !
  ! For CO2 and 14CO2, the carbon solubility and speciation are calculated by the
  ! subroutine co2calc, following the OCMIP2 protocol. These calculations are both made
  ! using total CO2, following which the surface CO2 concentration (CO2*, also known as
  ! H2CO3*) is scaled by the DI14C/DIC ratio to give the surface 14CO2 concentration.
  ! The speciation calculation uses in situ temperature, salinity, and PO4. 
  !
  !
  ! Oxygen solubility is calculated here, using in situ temperature and salinity.  

    !---------------------------------------------------------------------
    ! Get positive tracer concentrations for carbon calculation
    !---------------------------------------------------------------------

    allocate(tmp_irr_band(nbands))

    do n = 1, num_instances  

      bling%zbot = 0.0
      s_over_p = 0.0

    !---------------------------------------------------------------------
    ! Get positive concentrations for prognostic tracers
    !---------------------------------------------------------------------
      call g_tracer_get_values(tracer_list, 'po4' // bling%suffix, 'field', bling%f_po4, isd, jsd,          &
           ntau = tau, positive = .true.)
      if (bling%fe_is_prognostic) then
        call g_tracer_get_values(tracer_list, 'fed' // bling%suffix, 'field', bling%f_fed, isd, jsd,        &
             ntau = tau, positive = .true.)
      else
        call data_override('OCN', 'fed_data' // trim(bling%suffix), bling%f_fed_data, model_time)
        do k = 1, nk
          do j = jsc, jec
            do i = isc, iec
              bling%f_fed_data(i,j,k) =                                                                           &
                   max(bling%f_fed_data(i,j,k),                                                                   &
                       bling%fe_coastal * (1.0 - grid_ht(i,j)/bling%fe_coastal_depth)) * grid_tmask(i,j,k)
            enddo  !} i
          enddo  !} j
        enddo  !} k
        if (bling%fe_is_diagnostic) then
          call g_tracer_get_values(tracer_list, 'fed' // bling%suffix, 'field', bling%f_fed, isd, jsd,      &
               positive = .true.)
        else
          do k = 1, nk
            do j = jsc, jec
              do i = isc, iec
                bling%f_fed(i,j,k) = bling%f_fed_data(i,j,k)
              enddo  !} i
            enddo  !} j
          enddo  !} k
        endif
      endif
      call g_tracer_get_values(tracer_list, 'o2'  // bling%suffix, 'field', bling%f_o2,  isd, jsd,          &
           ntau = tau, positive = .true.)

    !---------------------------------------------------------------------
    ! Assign pointers for diagnostic tracers
    !---------------------------------------------------------------------
      if (bling%biomass_type .eq. 'single') then
        call g_tracer_get_pointer(tracer_list,'biomass_p' // bling%suffix,'field',bling%p_biomass_p)
      elseif (bling%biomass_type .eq. 'lg_sm_phyto') then
        call g_tracer_get_pointer(tracer_list,'phyto_lg'// bling%suffix,'field',bling%p_phyto_lg)
        call g_tracer_get_pointer(tracer_list,'phyto_sm'// bling%suffix,'field',bling%p_phyto_sm)
      endif
      call g_tracer_get_pointer(tracer_list,'irr_mem' // bling%suffix,'field',bling%p_irr_mem)


      if (bling%do_carbon) then                                    !<<CARBON CYCLE

        call g_tracer_get_pointer(tracer_list, 'htotal' // bling%suffix, 'field', bling%p_htotal)
        call g_tracer_get_pointer(tracer_list, 'dic'    // bling%suffix, 'field', bling%p_dic)

    !---------------------------------------------------------------------
    ! Calculate co2 fluxes csurf and alpha for the next round of exchange
    ! Note a scaled value of the PO4, rather than SiOH3, is used for all 
    ! calculations since there is no prognostic silica cycle. 
    ! Alkalinity is calculated from salinity, since there is no prognostic
    ! alkalinity cycle.
    !---------------------------------------------------------------------
   
        k=1
        do j = jsc, jec  
          do i = isc, iec  
            bling%htotallo(i,j) = bling%htotal_scale_lo * bling%p_htotal(i,j,k)
            bling%htotalhi(i,j) = bling%htotal_scale_hi * bling%p_htotal(i,j,k)
          enddo  !} i
        enddo  !}  j
     
        do j = jsc, jec  !{
          do i = isc, iec  !{
            bling%surf_temp(i,j) = Temp(i,j,k)
            bling%surf_salt(i,j) = Salt(i,j,k)
            bling%surf_po4(i,j)  = bling%f_po4(i,j,k)
            bling%surf_sio4(i,j) = bling%f_po4(i,j,k)
            bling%surf_dic(i,j)  = bling%p_dic(i,j,k,tau)
          enddo  !} i
        enddo  !} j


!       Optionally override surface values used in the gas exchange calculations
!       Note that the data_override routine wants the array to be over the computational grid,
!       so we need to pass only that part of the array (the arrays must be dimensioned on the data
!       domain, as that is what is needed for the co2calc routine). Also note that we only call the
!       data_override routine if we are actually overriding to avoid the implicit array copies implied
!       by passing a sub-array into and out of the subroutine.
!
!       There could be a problem with this scheme if the halo region values are actually used, as these may
!       be inconsistent with the overridden values on adjacent processors. I do not believe that this is
!       problem, however. -- Richard Slater (2012-02-07)

        if (bling%override_surf_temp) then
          call data_override('OCN', 'temp_co2_flux' // trim(bling%suffix),bling%surf_temp(isc:iec,jsc:jec),&
               model_time, override = bling%override_surf_temp)
        endif
        if (bling%override_surf_salt) then
          call data_override('OCN', 'salt_co2_flux' // trim(bling%suffix),bling%surf_salt(isc:iec,jsc:jec),&
               model_time, override = bling%override_surf_salt)
        endif
        if (bling%override_surf_alk) then
          call data_override('OCN', 'alk_co2_flux'  // trim(bling%suffix),bling%surf_alk(isc:iec,jsc:jec), &
               model_time, override = bling%override_surf_alk)
        endif
        if (bling%override_surf_po4) then
          call data_override('OCN', 'po4_co2_flux'  // trim(bling%suffix),bling%surf_po4(isc:iec,jsc:jec), &
               model_time, override = bling%override_surf_po4)
        endif
        if (bling%override_surf_sio4) then
          call data_override('OCN', 'sio4_co2_flux' // trim(bling%suffix),bling%surf_sio4(isc:iec,jsc:jec),&
               model_time, override = bling%override_surf_sio4)
        endif
        if (bling%override_surf_dic) then
          call data_override('OCN', 'dic_co2_flux'  // trim(bling%suffix),bling%surf_dic(isc:iec,jsc:jec), &
               model_time, override = bling%override_surf_dic)
        endif
     
        if (.not. bling%override_surf_alk) then
          !
          !     Calculate the surface alkalinity if not overridden above
          !
          if (bling%alk_scheme .eq. 'normal') then
            do j = jsc, jec  !{
              do i = isc, iec  !{
                ! This is an ad hoc regression, eyeballed from GLODAP vs WOA in ferret, to give ALK from salinity.
                ! Intercept is large, to keep the Southern Ocean alkalinity close to obs. Note this makes 
                ! the gyres low alkalinity, which will lead to outgassing there.
                ! Would probably be better to use an ALK/salinity map instead.
                bling%surf_alk(i,j)  = Salt(i,j,1) * bling%alk_slope + bling%alk_intercept
              enddo  !} i
            enddo  !} j
          elseif (bling%alk_scheme .eq. 'ratios') then
            call data_override('OCN', 'surf_alk' // trim(bling%suffix), bling%surf_alk(isc:iec,jsc:jec), model_time)
            do j = jsc, jec  !{
              do i = isc, iec  !{
                ! Use a map of the ratios of alkalinity to salinity and a constant intercept
                ! For this method, we use the data_override value for surf_alk to set the slopes
                bling%surf_alk(i,j)  = Salt(i,j,1) * bling%surf_alk(i,j) + bling%alk_intercept
              enddo  !} i
            enddo  !} j
          elseif (bling%alk_scheme .eq. 'intercepts') then
            call data_override('OCN', 'surf_alk' // trim(bling%suffix), bling%surf_alk(isc:iec,jsc:jec), model_time)
            do j = jsc, jec  !{
              do i = isc, iec  !{
                ! Use a map of the intercepts of alkalinity to salinity and a constant slope
                ! For this method, we use the data_override value for surf_alk to set the intercepts
                bling%surf_alk(i,j)  = Salt(i,j,1) * bling%alk_slope + bling%surf_alk(i,j)
              enddo  !} i
            enddo  !} j
          else
            call mpp_error(FATAL, trim(error_header) //        &
                 ' Illegal alk_scheme (' // trim(bling%alk_scheme) // ') for instance ' // trim(bling%name))
          endif
        endif

        call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k),        &
             bling%surf_temp(:,:), bling%surf_salt(:,:),    &
             bling%surf_dic(:,:),                                 &
             bling%surf_po4(:,:),                                 &
             bling%surf_sio4(:,:),                                &
             bling%surf_alk(:,:),                                 &
             bling%htotallo, bling%htotalhi,                &
                                    !InOut
             bling%p_htotal(:,:,k),                               &
                                    !OUT
             co2star=bling%co2_csurf(:,:),                        &
             alpha=bling%co2_alpha(:,:),                          &
             pCO2surf=bling%pco2_surf(:,:))                       

        call g_tracer_set_values(tracer_list,'dic' // bling%suffix,'alpha',bling%co2_alpha    ,isd,jsd)
        call g_tracer_set_values(tracer_list,'dic' // bling%suffix,'csurf',bling%co2_csurf    ,isd,jsd)

    
        if (bling%do_14c) then                                      !<<RADIOCARBON
      
          call g_tracer_get_pointer(tracer_list,'di14c' // bling%suffix ,'field', bling%p_di14c)
        
          do j = jsc, jec  
            do i = isc, iec  
            
              ! The surface p14CO2* concentration is calculated by scaling the total CO2* by the
              !   surface water 14C/12C.
              
              bling%c14o2_csurf(i,j) =  bling%co2_csurf(i,j) *                &
                   bling%p_di14c(i,j,1,tau) / (bling%p_dic(i,j,1,tau) + epsln)

              ! Alpha is here the same as co2. The air-sea flux depends on the atmospheric 
              ! p14CO2 given in the data table entry (which may vary over time, reflecting
              ! both changes in atmospheric pCO2 and D14CO2).

              bling%c14o2_alpha(i,j) =  bling%co2_alpha(i,j) 
              
            enddo  !} i
          enddo  !}  j

          call g_tracer_set_values(tracer_list,'di14c' // bling%suffix,'alpha',bling%c14o2_alpha ,isd,jsd)
          call g_tracer_set_values(tracer_list,'di14c' // bling%suffix,'csurf',bling%c14o2_csurf ,isd,jsd)

        endif                                                 !RADIOCARBON>>
      endif                                                  !CARBON CYCLE>>
    

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

      if (bling%biomass_type .eq. 'single') then

        do j = jsc, jec  
          do i = isc, iec  

            do nb = 1,nbands 
              if (max_wavelength_band(nb) .lt. 710) then 
                 tmp_irr_band(nb) = max(0.0,sw_pen_band(nb,i,j))
              else
                 tmp_irr_band(nb) = 0.0
              endif 
            enddo !} nbands

            kblt = 0
            tmp_irrad_ML = 0.0
            tmp_hblt = 0.0
            do k = 1, nk 
              tmp_Irrad = 0.0
              do nb = 1,nbands 
                tmp_opacity = opacity_band(nb,i,j,k)
                tmp_Irrad = tmp_Irrad + tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k) * 0.5)
                ! Change tmp_irr_band from being the value atop layer k to the value at the bottom of layer k.
                tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k))
              enddo !}  nbands
              bling%irr_inst(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
              bling%irr_mix(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
              if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then 
                kblt = kblt+1
                tmp_irrad_ML    = tmp_irrad_ML    + bling%irr_mix(i,j,k)    * dzt(i,j,k)
                tmp_hblt = tmp_hblt + dzt(i,j,k)
              endif 
            enddo !} k
            bling%irr_mix(i,j,1:kblt)    = tmp_irrad_ML / max(1.0e-6,tmp_hblt)

          enddo  !} i
        enddo  !} j

      elseif (bling%biomass_type .eq. 'lg_sm_phyto') then

        do j = jsc, jec  
          do i = isc, iec  

            do nb = 1,nbands 
              if (max_wavelength_band(nb) .lt. 710) then 
                 tmp_irr_band(nb) = max(0.0,sw_pen_band(nb,i,j))
              else
                 tmp_irr_band(nb) = 0.0
              endif 
            enddo !} nbands

            kblt = 0
            tmp_irrad_ML = 0.0
            tmp_phyto_lg_ML = 0.0
            tmp_phyto_sm_ML = 0.0
            tmp_hblt = 0.0
            do k = 1, nk 
              tmp_Irrad = 0.0
              do nb = 1,nbands 
                tmp_opacity = opacity_band(nb,i,j,k)
                tmp_Irrad = tmp_Irrad + tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k) * 0.5)
                ! Change tmp_irr_band from being the value atop layer k to the value at the bottom of layer k.
                tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k))
              enddo !}  nbands
              bling%irr_inst(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
              bling%irr_mix(i,j,k) = tmp_Irrad * grid_tmask(i,j,k)
              if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then 
                kblt = kblt+1
                tmp_irrad_ML    = tmp_irrad_ML    + bling%irr_mix(i,j,k)    * dzt(i,j,k)
                tmp_phyto_lg_ML = tmp_phyto_lg_ML + bling%p_phyto_lg(i,j,k) * dzt(i,j,k)
                tmp_phyto_sm_ML = tmp_phyto_sm_ML + bling%p_phyto_sm(i,j,k) * dzt(i,j,k)
                tmp_hblt = tmp_hblt + dzt(i,j,k)
              endif 
            enddo !} k
            bling%irr_mix(i,j,1:kblt)    = tmp_irrad_ML / max(1.0e-6,tmp_hblt)
            bling%p_phyto_lg(i,j,1:kblt) = tmp_phyto_lg_ML / max(1.0e-6,tmp_hblt)
            bling%p_phyto_lg(i,j,1:kblt) = tmp_phyto_sm_ML / max(1.0e-6,tmp_hblt)

          enddo  !} i
        enddo  !} j

      endif

      do k = 1, nk  
        do j = jsc, jec  
          do i = isc, iec         

       !--------------------------------------------------------------------
       ! Phytoplankton photoadaptation. This represents the fact that phytoplankton cells are
       ! adapted to the averaged light field to which they've been exposed over their lifetimes,
       ! rather than the instantaneous light. The timescale is set by gamma_irr_mem.
       
            bling%p_irr_mem(i,j,k) = (bling%p_irr_mem(i,j,k) +                   &
                 (bling%irr_mix(i,j,k) - bling%p_irr_mem(i,j,k)) * min( 1.0 ,   &
                 bling%gamma_irr_mem * dtts)) * grid_tmask(i,j,k)

       !--------------------------------------------------------------------
       ! Temperature functionality of growth and grazing
       ! NB The temperature effect of Eppley (1972) is used instead
       !    of that in Geider et al (1997) for both simplicity and
       !    to incorporate combined effects on uptake, incorporation
       !    into organic matter and photorespiration.  Values of PCmax
       !    are normalized to 0C rather than 20C in Geider et al. (1997)
        
            bling%expkT(i,j,k) = exp(bling%kappa_eppley * Temp(i,j,k))
     
          enddo  !} i
        enddo  !} j
      enddo  !} k

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
    ! of Michaelis-Menton PO4-limitation, or iron-limitation). The iron
    ! limitation term is scaled by (k_fe_2_p + fe_2_p_max) / fe_2_p_max
    ! so that it approaches 1 as fed approaches infinity. Thus, 
    ! it's of comparable magnitude to the PO4 limitation term.
    !
    ! Fe limitation acts by reducing the maximum achievable Chl:C ratio 
    ! (theta) below a prescribed, Fe-replete maximum value (thetamax), to 
    ! approach a prescribed minimum Chl:C (thetamin) under extreme
    ! Fe-limitation.
    !-----------------------------------------------------------------------
    
      do k = 1, nk  
        do j = jsc, jec  
          do i = isc, iec         
            bling%fe_2_p_uptake(i,j,k) = bling%fe_2_p_max *                         &
              bling%f_fed(i,j,k) / (bling%k_fe_uptake + bling%f_fed(i,j,k))
            bling%def_fe(i,j,k) = max(bling%def_fe_min,                             &
              (bling%fe_2_p_uptake(i,j,k) /                                               &
               (bling%k_fe_2_p + bling%fe_2_p_uptake(i,j,k)) *                      &
               (bling%k_fe_2_p + bling%fe_2_p_max) / bling%fe_2_p_max)) 
            bling%pc_m(i,j,k) = bling%pc_0 * bling%expkT(i,j,k) * min(        &
              max(0.,((bling%f_po4(i,j,k) - bling%po4_min) /                        & 
              (bling%k_po4 + bling%f_po4(i,j,k) - bling%po4_min))) ,          &
              bling%def_fe(i,j,k))
            bling%thetamax_fe(i,j,k) = bling%thetamax_lo +                          &
              (bling%thetamax_hi - bling%thetamax_lo) * bling%def_fe(i,j,k)
            
    !-----------------------------------------------------------------------
    ! Next, the nutrient-limited efficiency of algal photosystems, Irrk, is
    ! calculated. This requires a prescribed quantum yield, alpha.
    ! The iron deficiency term is included here as a multiplier of the 
    ! thetamax_fe to represent the importance of Fe in forming chlorophyll
    ! accessory antennae, which do not affect the Chl:C but still affect the
    ! phytoplankton ability to use light (eg Stzrepek & Harrison Nature 
    ! 2004).

            bling%irrk(i,j,k) = (bling%pc_m(i,j,k) / ( epsln +               &
              bling%alpha_photo * bling%thetamax_fe(i,j,k) )) +             &
              bling%p_irr_mem(i,j,k) * 0.5

    !-----------------------------------------------------------------------
    ! We also calculate the Chl:C ratio here, although it does not enter  
    ! into the uptake calculation and is only used for the diagnostic
    ! chlorophyll concentration, below.

            bling%theta(i,j,k) = bling%thetamax_fe(i,j,k) / (1. +              &
              bling%thetamax_fe(i,j,k) * bling%alpha_photo *                  &
              bling%p_irr_mem(i,j,k) / (epsln + 2. * bling%pc_m(i,j,k)))
            
    !-----------------------------------------------------------------------
    ! Now we can calculate the carbon-specific photosynthesis rate, mu. 
    
            bling%mu(i,j,k) = bling%pc_m(i,j,k) *                              &
              (1. - exp(-bling%irr_mix(i,j,k) / (epsln + bling%irrk(i,j,k))))

          enddo  !} i
        enddo  !} j
      enddo  !} k

    !-----------------------------------------------------------------------
    ! We now must convert this net carbon-specific growth rate to nutrient 
    ! uptake rates, the quantities we are interested in. Since we have no 
    ! explicit biomass tracer, we use the result of Dunne et al. (GBC, 2005) 
    ! to calculate an implicit biomass from the uptake rate through the  
    ! application of a simple idealized grazing law. This has the effect of 
    ! reducing uptake in low growth-rate regimes and increasing uptake in 
    ! high growth-rate regimes - essentially a non-linear amplification of 
    ! the growth rate variability. The result is:

      if (bling%biomass_type .eq. 'single') then
        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec         
      
              bling%biomass_p_ts(i,j,k) =                                                  &
                ((bling%mu(i,j,k)/(bling%lambda0 * bling%expkT(i,j,k)))**3     &
                + (bling%mu(i,j,k)/(bling%lambda0 * bling%expkT(i,j,k))))      &
                * bling%p_star

              bling%p_biomass_p(i,j,k) = bling%p_biomass_p(i,j,k) +            &
                (bling%biomass_p_ts(i,j,k) - bling%p_biomass_p(i,j,k)) *       &
                min(1.0, bling%gamma_biomass * dtts) * grid_tmask(i,j,k)
              
              bling%jp_uptake(i,j,k) = bling%p_biomass_p(i,j,k) *              &
                bling%mu(i,j,k)

    ! We can now use the diagnostic biomass to calculate the chlorophyll
    ! concentration:
    
              bling%f_chl(i,j,k) = max(bling%chl_min, bling%p_biomass_p(i,j,k) &
                * bling%c_2_p * 12.011e6 * bling%theta(i,j,k)) *                     &
                grid_tmask(i,j,k)
 
     ! As a helpful diagnostic, the implied fraction of production by large 
     ! phytoplankton is calculated, also following Dunne et al. 2005. This
     ! could be done more simply, but is done here in a complicated way as
     ! a sanity check. Note the calculation is made in P units, rather than C.

              s_over_p = ( -1. + ( 1. + 4. * bling%jp_uptake(i,j,k) /                         &
                (bling%expkT(i,j,k) * bling%lambda0 * bling%p_star))**0.5) * .5
              bling%frac_lg(i,j,k) = s_over_p / (1 + s_over_p)

            enddo  !} i
          enddo  !} j
        enddo  !} k

      elseif (bling%biomass_type .eq. 'lg_sm_phyto') then

        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec         
              bling%jp_uptake(i,j,k) = bling%mu(i,j,k)  *              &
                (bling%p_phyto_lg(i,j,k) + bling%p_phyto_sm(i,j,k))
            enddo  !} i
          enddo  !} j
        enddo  !} k

      endif

    !-----------------------------------------------------------------------
    ! Iron is then taken up as a function of PO4 uptake and iron limitation,
    ! with a maximum Fe:P uptake ratio of fe2p_max:
    
      do k = 1, nk  
        do j = jsc, jec  
          do i = isc, iec         
            bling%jfe_uptake(i,j,k) = bling%jp_uptake(i,j,k) *               &
              bling%fe_2_p_uptake(i,j,k)
          enddo  !} i
        enddo  !} j
      enddo  !} k
    
    
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
  ! The non-sinking fraction is recycled instantaneously to the inorganic 
  ! nutrient pool,
  ! representing the fast turnover of labile dissolved organic matter via
  ! the microbial loop, and the remainder is converted to semi-labile
  ! dissolved organic matter. Iron and phosphorus are treated identically 
  ! for the first step, but all iron is recycled instantaneously in the
  ! second step (i.e. there is no dissolved organic iron pool).
  !-------------------------------------------------------------------------
    
      do k = 1, nk  
        do j = jsc, jec  
          do i = isc, iec        

            bling%frac_pop(i,j,k) = max((bling%phi_sm + bling%phi_lg *             &
              (bling%mu(i,j,k)/(bling%lambda0*bling%expkT(i,j,k)))**2.)/       &
              (1. + (bling%mu(i,j,k)/(bling%lambda0*bling%expkT(i,j,k)))**2.)* &
              exp(bling%kappa_remin * Temp(i,j,k)) *                                       &
              ! Experimental! Reduce frac_pop under strong PO4 limitation
              bling%f_po4(i,j,k) / (bling%k_po4_recycle + bling%f_po4(i,j,k)),  &
              bling%min_frac_pop)

            bling%jpop(i,j,k) = bling%frac_pop(i,j,k) * bling%jp_uptake(i,j,k)
     
     ! Whatever isn't converted to sinking particulate is recycled to the dissolved pool.
     
            bling%jp_recycle(i,j,k) = bling%jp_uptake(i,j,k) -                  &
              bling%jpop(i,j,k)

          enddo  !] i
        enddo  !} j
      enddo  !} k
    
      if (bling%biomass_type .eq. 'lg_sm_phyto') then

        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec        

      ! Finally, update the biomass of total phytoplankton, and of diazotrophs.
      ! Use this to solve the Dunne et al. 2005 mortality term, with alpha=1/3 (eq. 5b). 
      ! Then, add this to the pre-exisiting phytoplankton biomass and the total uptake to give
                
              bling%p_phyto_lg(i,j,k) = bling%p_phyto_lg(i,j,k) +                              &
                bling%p_phyto_lg(i,j,k) * (bling%mu(i,j,k) -                                   &
                bling%lambda0 * bling%expkT(i,j,k) *                                           &
                (bling%p_phyto_lg(i,j,k) / bling%p_star)**(1./3.) ) *  dtts * grid_tmask(i,j,k)
                
              bling%p_phyto_sm(i,j,k) = bling%p_phyto_sm(i,j,k) +                         &
                bling%p_phyto_sm(i,j,k) * (bling%mu(i,j,k) -                              &
                bling%lambda0 * bling%expkT(i,j,k) *                                      &
                (bling%p_phyto_sm(i,j,k) / bling%p_star) ) * dtts * grid_tmask(i,j,k)

              bling%frac_lg(i,j,k) = bling%p_phyto_lg(i,j,k) / &
                (epsln + bling%p_phyto_lg(i,j,k)+bling%p_phyto_sm(i,j,k))

       ! Calculate the chlorophyll concentration:
      
               bling%f_chl(i,j,k) = max(bling%chl_min,                                  &
                 bling%c_2_p * 12.011e6 * bling%theta(i,j,k) *                          &
                 (bling%p_phyto_lg(i,j,k) + bling%p_phyto_sm(i,j,k))) * grid_tmask(i,j,k)

            enddo  !] i
          enddo  !} j
        enddo  !} k

      endif

      !
      ! perform recycling, as above, for the prognostic Fed tracer
      !
      if (bling%fe_is_prognostic) then
        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec        

              bling%jfeop(i,j,k) = bling%frac_pop(i,j,k)*bling%jfe_uptake(i,j,k)

              bling%jfe_recycle(i,j,k) = bling%jfe_uptake(i,j,k) -                &
                bling%jfeop(i,j,k)

            enddo  !] i
          enddo  !} j
        enddo  !} k
      endif
    

  !-------------------------------------------------------------------------
  ! SINKING AND REMINERALIZATION
  !-------------------------------------------------------------------------
    ! Calculate the depth of each grid cell (needs to be 3d for use with
    ! isopycnal co-ordinate model).

      do j = jsc, jec  
        do i = isc, iec  
          bling%zbot(i,j,1) = dzt(i,j,1)
        enddo  !} i
      enddo  !} j

      do k = 2, nk   
        do j = jsc, jec   
          do i = isc, iec   
            bling%zbot(i,j,k) = bling%zbot(i,j,k-1) + dzt(i,j,k)
          enddo  !} i
        enddo  !} j
      enddo  !} k

    !-----------------------------------------------------------------------
    ! Calculate the remineralization lengthscale matrix, zremin, a function 
    ! of z. Sinking rate (wsink) is constant over the upper wsink0_z metres,
    ! then  increases linearly with depth.
    ! The remineralization rate is a function of oxygen concentrations,
    ! to slow remineralization under suboxia/anoxia. The remineralization rate 
    ! approaches the remin_min as O2 approaches O2 min.
    
      do k = 1, nk  
        do j = jsc, jec  
          do i = isc, iec         
      
            if (bling%zbot(i,j,k) .lt. bling%wsink0_z) then  
              bling%wsink(i,j,k) = bling%wsink0
            else  
              bling%wsink(i,j,k) = (bling%wsink_acc * (bling%zbot(i,j,k) -     &
                bling%wsink0_z) + bling%wsink0)
            endif  

            bling%zremin(i,j,k) = bling%gamma_pop * (bling%f_o2(i,j,k) /    &
              (bling%k_o2 + bling%f_o2(i,j,k)) * (1. - bling%remin_min)+ &
              bling%remin_min) / (bling%wsink(i,j,k) + epsln)

          enddo  !} i
        enddo  !} j
      enddo  !} k

      if (bling%do_carbon) then                                    !<<CARBON CYCLE
 
        if (bling%do_14c) then                                      !<<RADIOCARBON

      ! Sinking particulate 14C is generated in the local ratio of 14C/12C
      ! to sinking 12C, which itself is strictly tied to P through a fixed
      ! C:P. Therefore, jpop can be used to calculate fpo14c.

          do j = jsc, jec  
            do i = isc, iec  
              bling%c14_2_p(i,j,1) = bling%c_2_p *                                   &
                bling%p_di14c(i,j,1,tau) / (epsln + bling%p_dic(i,j,1,tau))

              bling%fpo14c(i,j,1) = bling%jpop(i,j,1) * bling%c14_2_p(i,j,1) * &
                rho_dzt(i,j,1) / (1.0 + dzt(i,j,1) * bling%zremin(i,j,1))                                   

              bling%j14c_reminp(i,j,1) = (bling%jpop(i,j,1) *                        &
                bling%c14_2_p(i,j,1) * rho_dzt(i,j,1) - bling%fpo14c(i,j,1)) /       &
               (epsln + rho_dzt(i,j,1))
            enddo  !} i
          enddo  !} j

          do k = 2, nk  
            do j = jsc, jec  
              do i = isc, iec  
                bling%fpo14c(i,j,k) = (bling%fpo14c(i,j,k-1) +                       &
                  bling%jpop(i,j,k) * bling%c14_2_p(i,j,k) * rho_dzt(i,j,k)) /       &
                  (1.0 + dzt(i,j,k) * bling%zremin(i,j,k)) 

                bling%j14c_reminp(i,j,k) = (bling%fpo14c(i,j,k-1) +                  &
                 bling%jpop(i,j,k) * bling%c14_2_p(i,j,k) * rho_dzt(i,j,k) -         &
                 bling%fpo14c(i,j,k)) / (epsln + rho_dzt(i,j,k))
              enddo  !} i
            enddo  !} j
          enddo  !} k

     ! Decay the radiocarbon in DIC
     
          bling%lambda_14c = log(2.0) / (bling%half_life_14c * spery)

          do k = 1, nk  
            do j = jsc, jec  
              do i = isc, iec          

                bling%j14c_decay_dic(i,j,k) = bling%p_di14c(i,j,k,tau) *               &
                  bling%lambda_14c 

              enddo  !} i
            enddo  !} j
          enddo  !} k
          
          
        endif                                                 !RADIOCARBON>>
      endif                                                 !CARBON CYCLE>>

      if (bling%fe_is_prognostic) then
        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec          

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
     ! to be produced as a response to extreme iron stress.
     ! In anoxic waters, iron should be reduced, and therefore mostly 
     ! immune to scavenging. Easiest way to do this is to skip the feprime
     ! calculation if oxygen is less than 0.
     
              if (bling%f_o2(i,j,k) .gt. bling%o2_min) then  
                bling%kfe_eq_lig(i,j,k) = bling%kfe_eq_lig_max -                            &
                  (bling%kfe_eq_lig_max - bling%kfe_eq_lig_min) *                           &
                  (bling%irr_inst(i,j,k)**2. / (bling%irr_inst(i,j,k)**2. +                 &
                  bling%kfe_eq_lig_irr **2.)) * max(0., min(1., (bling%f_fed(i,j,k) -       &
                  bling%kfe_eq_lig_femin) / (epsln + bling%f_fed(i,j,k)) * 1.2))

                bling%feprime(i,j,k) = 1.0 + bling%kfe_eq_lig(i,j,k) *                      &
                  (bling%felig_bkg - bling%f_fed(i,j,k))
               
                bling%feprime(i,j,k) = (-bling%feprime(i,j,k) +(bling%feprime(i,j,k)* &
                  bling%feprime(i,j,k) + 4.0 * bling%kfe_eq_lig(i,j,k) *                    &
                  bling%f_fed(i,j,k))**(0.5)) /(2.0 * bling%kfe_eq_lig(i,j,k))
              else  !}{
                bling%feprime(i,j,k) = 0.
              endif  !}

              bling%jfe_ads_inorg(i,j,k) = min(0.5/dtts, bling%kfe_inorg *                  &
                 bling%feprime(i,j,k) ** 0.5) * bling%feprime(i,j,k) 
       
            enddo  !} i
          enddo  !} j
        enddo  !} k
      endif

    !---------------------------------------------------------------------
    ! In general, the flux at the bottom of a grid cell should equal
    ! Fb = (Ft + Prod*dz) / (1 + zremin*dz)
    ! where Ft is the flux at the top, and prod*dz is the integrated 
    ! production of new sinking particles within the layer.
    ! Since Ft=0 in the first layer,

      do j = jsc, jec  
        do i = isc, iec  

          bling%fpop(i,j,1) = bling%jpop(i,j,1) * rho_dzt(i,j,1) /                    &
            (1.0 + dzt(i,j,1) * bling%zremin(i,j,1)) 

    !-----------------------------------------------------------------------
    ! Calculate remineralization terms

          bling%jp_reminp(i,j,1) =                                                          &
           (bling%jpop(i,j,1) * rho_dzt(i,j,1) - bling%fpop(i,j,1)) /                 &
           (epsln + rho_dzt(i,j,1))
           
        enddo  !} i
      enddo  !} j
     

    !-----------------------------------------------------------------------
    ! Then, for the rest of water column, include flux from above:

      do k = 2, nk  
        do j = jsc, jec  
          do i = isc, iec  

            bling%fpop(i,j,k) = (bling%fpop(i,j,k-1) +                           &
              bling%jpop(i,j,k) * rho_dzt(i,j,k)) /                                    &
              (1.0 + dzt(i,j,k) * bling%zremin(i,j,k)) 

    !---------------------------------------------------------------------
    ! Calculate remineralization terms

            bling%jp_reminp(i,j,k) = (bling%fpop(i,j,k-1) +                      &
             bling%jpop(i,j,k) * rho_dzt(i,j,k) - bling%fpop(i,j,k)) /           &
             (epsln + rho_dzt(i,j,k))
             
          enddo  !} i
        enddo  !} j
      enddo  !} k
    

    !---------------------------------------------------------------------
    ! BOTTOM LAYER 
    ! Account for remineralization in bottom box, and bottom fluxes

      do j = jsc, jec  
        do i = isc, iec  
          k = grid_kmt(i,j)
          if (k .gt. 0) then 

      !---------------------------------------------------------------------
      ! Calculate external bottom fluxes for tracer_vertdiff. Positive fluxes
      ! are from the water column into the seafloor. For P, the bottom flux  
      ! puts the sinking flux reaching the bottom cell into the water column 
      ! through diffusion.
      ! For oxygen, the consumption of oxidant required to respire  
      ! the settling flux of organic matter (in support of the
      ! PO4 bottom flux) diffuses from the bottom water into the sediment.

            bling%b_po4(i,j) = - bling%fpop(i,j,k)

            if (bling%f_o2(i,j,k) .gt. bling%o2_min) then  
              bling%b_o2(i,j) = bling%o2_2_p * bling%fpop(i,j,k)
            else
              bling%b_o2(i,j) = 0.0
            endif 

          endif 
        enddo  !} i
      enddo  !} j

      if (bling%fe_is_prognostic) then

        do j = jsc, jec  
          do i = isc, iec  

    !-----------------------------------------------------------------------
    ! Now, calculate the Fe adsorption using this fpop:
    ! The absolute first order rate constant is calculated from the 
    ! concentration of organic particles, after Parekh et al. (2005). Never
    !  allowed to be greater than 1/2dt for numerical stability.

            bling%jfe_ads_org(i,j,1) = min (0.5/dtts,                                         &
              bling%kfe_org * (bling%fpop(i,j,1) / (epsln + bling%wsink(i,j,1)) * &
              bling%mass_2_p) ** 0.58) * bling%feprime(i,j,1)
             
            bling%fpofe(i,j,1) = (bling%jfeop(i,j,1) +bling%jfe_ads_inorg(i,j,1)  &
               + bling%jfe_ads_org(i,j,1)) * rho_dzt(i,j,1) /                                 &
               (1.0 + dzt(i,j,1) * bling%zremin(i,j,1)) 
      
    !-----------------------------------------------------------------------
    ! Calculate remineralization terms

            bling%jfe_reminp(i,j,1) =                                                         &
             ((bling%jfeop(i,j,1) + bling%jfe_ads_org(i,j,1) +                          &
             bling%jfe_ads_inorg(i,j,1)) * rho_dzt(i,j,1) -                                   &
             bling%fpofe(i,j,1)) / (epsln + rho_dzt(i,j,1))

          enddo  !} i
        enddo  !} j
     

    !-----------------------------------------------------------------------
    ! Then, for the rest of water column, include flux from above:

        do k = 2, nk  
          do j = jsc, jec  
            do i = isc, iec  

    !-----------------------------------------------------------------------
    ! Again, calculate the Fe adsorption using this fpop:

              bling%jfe_ads_org(i,j,k) = min (0.5/dtts, bling%kfe_org *            &
                (bling%fpop(i,j,k) / (epsln + bling%wsink(i,j,k)) *                &
                bling%mass_2_p) ** 0.58) * bling%feprime(i,j,k)
               
              bling%fpofe(i,j,k) = (bling%fpofe(i,j,k-1) +                         &
                (bling%jfe_ads_org(i,j,k) + bling%jfe_ads_inorg(i,j,k) +           &
                bling%jfeop(i,j,k)) *rho_dzt(i,j,k)) /                                   &
                (1.0 + dzt(i,j,k) * bling%zremin(i,j,k)) 

    !---------------------------------------------------------------------
    ! Calculate remineralization terms

              bling%jfe_reminp(i,j,k) = (bling%fpofe(i,j,k-1) +                    &
               (bling%jfe_ads_org(i,j,k) + bling%jfe_ads_inorg(i,j,k) +            & 
                bling%jfeop(i,j,k)) * rho_dzt(i,j,k) -                                  &
               bling%fpofe(i,j,k)) / (epsln + rho_dzt(i,j,k))

            enddo  !} i
          enddo  !} j
        enddo  !} k
    

    !---------------------------------------------------------------------
    ! BOTTOM LAYER 
    ! Account for remineralization in bottom box, and bottom fluxes

        do j = jsc, jec  
          do i = isc, iec  
            k = grid_kmt(i,j)
            if (k .gt. 0) then 

      !---------------------------------------------------------------------
      ! Calculate iron addition from sediments as a function of organic
      ! matter supply.

              bling%ffe_sed(i,j) = bling%fe_2_p_sed * bling%fpop(i,j,k)

        ! Added the burial flux of sinking particulate iron here as a 
        ! diagnostic, needed to calculate mass balance of iron.

              bling%fe_burial(i,j) = bling%fpofe(i,j,k)

      !---------------------------------------------------------------------
      ! Calculate external bottom fluxes for tracer_vertdiff. Positive fluxes
      ! are from the water column into the seafloor.  For iron, the sinking flux disappears into the 
      ! sediments if bottom waters are oxic (assumed adsorbed as oxides),
      ! while an efflux of dissolved iron occurs dependent on the supply of
      ! reducing organic matter (scaled by the org-P sedimentation rate).
      ! If bottom waters are anoxic, the sinking flux of Fe is returned to
      ! the water column. Note this is not appropriate for very long runs
      ! with an anoxic ocean (iron will keep accumulating forever).

              if (bling%f_o2(i,j,k) .gt. bling%o2_min) then  
                bling%b_fed(i,j) = - bling%ffe_sed(i,j) 
              else
                bling%b_fed(i,j) = - bling%ffe_sed(i,j) - bling%fpofe(i,j,k)
              endif 

            endif 
          enddo  !} i
        enddo  !} j
      endif

      if (bling%fe_is_prognostic) then
        call g_tracer_set_values(tracer_list,'fed' // bling%suffix, 'btf', bling%b_fed ,isd,jsd)
      endif
      call g_tracer_set_values(tracer_list,'po4' // bling%suffix, 'btf', bling%b_po4 ,isd,jsd)
      call g_tracer_set_values(tracer_list,'o2' // bling%suffix,  'btf', bling%b_o2 ,isd,jsd)

      if (bling%do_carbon) then                              !<<CARBON CYCLE
    ! Do bottom box calcs for carbon cycle
    
        do j = jsc, jec  
          do i = isc, iec  
            k = grid_kmt(i,j)
            if (k .gt. 0) then 
                ! Do not bury any C-org - all goes back to water column  
                bling%b_dic(i,j) = - bling%fpop(i,j,k) * bling%c_2_p 
            endif  !}
          enddo  !} i
        enddo  !} j
        call g_tracer_set_values(tracer_list,'dic' // bling%suffix, 'btf', bling%b_dic ,isd,jsd)

        if (bling%do_14c) then                                      !<<RADIOCARBON
          do j = jsc, jec  
            do i = isc, iec  
              k = grid_kmt(i,j)
              if (k .gt. 0) then 
                bling%b_di14c(i,j) = - bling%fpo14c(i,j,k)
              endif  !}
            enddo  !} i
          enddo  !} j
          call g_tracer_set_values(tracer_list,'di14c' // bling%suffix,'btf',bling%b_di14c,isd,jsd)
        endif  !}                                               !RADIOCARBON>>

      endif  !}                                                !CARBON CYCLE>>


  !-------------------------------------------------------------------------
  !     CALCULATE SOURCE/SINK TERMS FOR EACH TRACER
  !-------------------------------------------------------------------------

    !Update the prognostics tracer fields via their pointers.

      if (bling%fe_is_prognostic) then
        call g_tracer_get_pointer(tracer_list, 'fed' // bling%suffix, 'field', bling%p_fed)
      elseif (bling%fe_is_diagnostic) then
        call g_tracer_get_pointer(tracer_list, 'fed' // bling%suffix, 'field', bling%p_fed_diag)
      endif
      call g_tracer_get_pointer(tracer_list,'o2' // bling%suffix     ,'field',bling%p_o2     )
      call g_tracer_get_pointer(tracer_list,'po4' // bling%suffix    ,'field',bling%p_po4    )

      if (bling%do_carbon) then  
        call g_tracer_get_pointer(tracer_list,'dic' // bling%suffix,'field',bling%p_dic)
        if (bling%do_14c) then  
          call g_tracer_get_pointer(tracer_list,'di14c' // bling%suffix,'field',bling%p_di14c)
        endif   !}
      endif   !}

      do k = 1, nk  
        do j = jsc, jec  
          do i = isc, iec  
       
       !
       ! PO4
       ! Sum of fast recycling and decay of sinking POP, less uptake.
       !    
            bling%jpo4(i,j,k) = bling%jp_recycle(i,j,k) +                       &
              bling%jp_reminp(i,j,k) - bling%jp_uptake(i,j,k)
      
            bling%p_po4(i,j,k,tau) = bling%p_po4(i,j,k,tau) +                   &
               bling%jpo4(i,j,k) * dtts * grid_tmask(i,j,k)

    !-----------------------------------------------------------------------
    !     O2
    ! Assuming constant P:O ratio.
    ! Optional prevention of negative oxygen (does not conserve ocean 
    ! redox potential) or alternatively it can be allowed to go negative, 
    ! keeping track of an implicit nitrate deficit 
    ! plus sulfate reduction.
    !-----------------------------------------------------------------------

            if ( (bling%prevent_neg_o2) .and.                                               &
                 (bling%f_o2(i,j,k) .lt. bling%o2_min) ) then 
              bling%jo2(i,j,k) = 0. * grid_tmask(i,j,k)
            else
              bling%jo2(i,j,k) = - bling%o2_2_p * bling%jpo4(i,j,k)             &
                * grid_tmask(i,j,k)
            endif !}

            bling%p_o2(i,j,k,tau) = bling%p_o2(i,j,k,tau) + bling%jo2(i,j,k) *  &
               dtts * grid_tmask(i,j,k)
          enddo  !} i
        enddo  !} j
      enddo  !} k

      !
      ! Fed
      !

      if (bling%fe_is_prognostic) then
        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec  
              bling%p_fed(i,j,k,tau) = bling%p_fed(i,j,k,tau) +                                             &
                (bling%jfe_recycle(i,j,k) + bling%jfe_reminp(i,j,k) -                                       &
                bling%jfe_uptake(i,j,k) - bling%jfe_ads_org(i,j,k) -                                        &
                bling%jfe_ads_inorg(i,j,k) ) * dtts * grid_tmask(i,j,k)
            enddo  !} i
          enddo  !} j
        enddo  !} k
      elseif (bling%fe_is_diagnostic) then
        do k = 1, nk
          do j = jsc, jec
            do i = isc, iec
              bling%p_fed_diag(i,j,k) = bling%p_fed_diag(i,j,k) -                                           &
                   bling%jfe_uptake(i,j,k) * dtts * grid_tmask(i,j,k)
              bling%jfe_reminp(i,j,k) = (bling%f_fed_data(i,j,k) - bling%p_fed_diag(i,j,k)) *         &
                   (1.0 / (bling%fe_restoring * 86400.0)) * grid_tmask(i,j,k)
              bling%p_fed_diag(i,j,k) = bling%p_fed_diag(i,j,k) +                                           &
                   bling%jfe_reminp(i,j,k) * dtts
            enddo  !} i
          enddo  !} j
        enddo  !} k
      endif

      if (bling%do_carbon) then                                    !<<CARBON CYCLE
        do k = 1, nk  
          do j = jsc, jec  
            do i = isc, iec  
    
              bling%p_dic(i,j,k,tau) = bling%p_dic(i,j,k,tau) +                     &
                (bling%jpo4(i,j,k) * bling%c_2_p) * dtts * grid_tmask(i,j,k)
             
              if (bling%do_14c) then                                      !<<RADIOCARBON
                bling%jdi14c(i,j,k) = (bling%jp_recycle(i,j,k) -                    &
                  bling%jp_uptake(i,j,k)) * bling%c14_2_p(i,j,k) +                  &
                  bling%j14c_reminp(i,j,k) 
        
                bling%p_di14c(i,j,k,tau) = bling%p_di14c(i,j,k,tau) +               &
                  (bling%jdi14c(i,j,k) - bling%j14c_decay_dic(i,j,k)) * dtts        &
                  * grid_tmask(i,j,k)
              endif  !}                                               !RADIOCARBON>>
       
            enddo  !} i
          enddo  !} j
        enddo  !} k
      endif                                                    !CARBON CYCLE>>

    !
    !Set the diagnostics tracer fields.
    !
      call g_tracer_set_values(tracer_list,'chl' // bling%suffix,'field',bling%f_chl,isd,jsd, &
           ntau=1)

    !-----------------------------------------------------------------------
    !       Save variables for diagnostics
    !-----------------------------------------------------------------------
    !

      if (.not. bling%fe_is_prognostic) then
        if (bling%id_fed_data_surf .gt. 0)                                                 &
             used = send_data(bling%id_fed_data_surf,    bling%f_fed_data(:,:,1),    &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      endif
      if (bling%id_htotal_surf .gt. 0)                                                   &
           used = send_data(bling%id_htotal_surf,    bling%p_htotal(:,:,1),        &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_chl_surf .gt. 0)                                                      &
           used = send_data(bling%id_chl_surf,    bling%f_chl(:,:,1),              &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%biomass_type .eq. 'single') then
        if (bling%id_biomass_p_surf .gt. 0)                                                &
             used = send_data(bling%id_biomass_p_surf,    bling%p_biomass_p(:,:,1),  &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      elseif (bling%biomass_type .eq. 'lg_sm_phyto') then
        if (bling%id_phyto_lg_surf .gt. 0)                                                &
             used = send_data(bling%id_phyto_lg_surf,    bling%p_phyto_lg(:,:,1),  &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
        if (bling%id_phyto_sm_surf .gt. 0)                                                &
             used = send_data(bling%id_phyto_sm_surf,    bling%p_phyto_sm(:,:,1),  &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      endif
      if (bling%id_irr_mem_surf .gt. 0)                                                  &
           used = send_data(bling%id_irr_mem_surf,    bling%p_irr_mem(:,:,1),      &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_pco2_surf .gt. 0)                                             &
           used = send_data(bling%id_pco2_surf,    bling%pco2_surf,        &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_temp_co2calc .gt. 0)                                             &
           used = send_data(bling%id_temp_co2calc,    bling%surf_temp,        &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_salt_co2calc .gt. 0)                                             &
           used = send_data(bling%id_salt_co2calc,    bling%surf_salt,        &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_po4_co2calc .gt. 0)                                              &
           used = send_data(bling%id_po4_co2calc,    bling%surf_po4,          &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_alk_co2calc .gt. 0)                                              &
           used = send_data(bling%id_alk_co2calc,    bling%surf_alk,          &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_sio4_co2calc .gt. 0)                                             &
           used = send_data(bling%id_sio4_co2calc,    bling%surf_sio4,        &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_dic_co2calc .gt. 0)                                              &
           used = send_data(bling%id_dic_co2calc,    bling%surf_dic,          &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%fe_is_prognostic) then
        if (bling%id_b_fed .gt. 0)                                                         &
             used = send_data(bling%id_b_fed,          bling%b_fed,                  &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      endif
      if (bling%id_b_o2 .gt. 0)                                                          &
           used = send_data(bling%id_b_o2,           bling%b_o2,                   &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_b_po4 .gt. 0)                                                         &
           used = send_data(bling%id_b_po4,          bling%b_po4,                  &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%biomass_type .eq. 'single') then
        if (bling%id_biomass_p_ts .gt. 0)                                                  &
             used = send_data(bling%id_biomass_p_ts,   bling%biomass_p_ts,           &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (bling%id_def_fe .gt. 0)                                                        &
           used = send_data(bling%id_def_fe,         bling%def_fe,                 &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_expkT .gt. 0)                                                         &
           used = send_data(bling%id_expkT,          bling%expkT,                  &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_fe_2_p_uptake .gt. 0)                                                 &
           used = send_data(bling%id_fe_2_p_uptake,  bling%fe_2_p_uptake,          &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%fe_is_prognostic) then
        if (bling%id_feprime .gt. 0)                                                       &
             used = send_data(bling%id_feprime,        bling%feprime,                &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
        if (bling%id_fe_burial .gt. 0)                                                     &
             used = send_data(bling%id_fe_burial,      bling%fe_burial,              &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
        if (bling%id_ffe_sed .gt. 0)                                                       &
             used = send_data(bling%id_ffe_sed,        bling%ffe_sed,                &
             model_time, rmask = grid_tmask(:,:,1),                                              & 
             is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
        if (bling%id_fpofe .gt. 0)                                                         &
             used = send_data(bling%id_fpofe,          bling%fpofe,                  &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (bling%id_fpop .gt. 0)                                                          &
           used = send_data(bling%id_fpop,           bling%fpop,                   &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_frac_lg .gt. 0)                                                       &
           used = send_data(bling%id_frac_lg,        bling%frac_lg,                &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_frac_pop .gt. 0)                                                      &
           used = send_data(bling%id_frac_pop,       bling%frac_pop,               &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_irr_inst .gt. 0)                                                      &
           used = send_data(bling%id_irr_inst,       bling%irr_inst,               &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_irr_mix .gt. 0)                                                       &
           used = send_data(bling%id_irr_mix,        bling%irr_mix,                &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_irrk .gt. 0)                                                          &
           used = send_data(bling%id_irrk,           bling%irrk,                   &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%fe_is_prognostic) then
        if (bling%id_jfe_ads_inorg .gt. 0)                                                 &
             used = send_data(bling%id_jfe_ads_inorg,  bling%jfe_ads_inorg*rho_dzt,  &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
        if (bling%id_jfe_ads_org .gt. 0)                                                   &
             used = send_data(bling%id_jfe_ads_org,    bling%jfe_ads_org*rho_dzt,    &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
        if (bling%id_jfe_recycle .gt. 0)                                                   &
             used = send_data(bling%id_jfe_recycle,    bling%jfe_recycle*rho_dzt,    &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (bling%fe_is_prognostic .or. bling%fe_is_diagnostic) then
        if (bling%id_jfe_reminp .gt. 0)                                                    &
             used = send_data(bling%id_jfe_reminp,     bling%jfe_reminp*rho_dzt,     &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (bling%id_jfe_uptake .gt. 0)                                                    &
           used = send_data(bling%id_jfe_uptake,     bling%jfe_uptake*rho_dzt,     &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_jo2 .gt. 0)                                                           &
           used = send_data(bling%id_jo2,            bling%jo2*rho_dzt,            &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_jp_recycle .gt. 0)                                                    &
           used = send_data(bling%id_jp_recycle,     bling%jp_recycle*rho_dzt,     &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_jp_reminp .gt. 0)                                                     &
           used = send_data(bling%id_jp_reminp,      bling%jp_reminp*rho_dzt,      &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_jp_uptake .gt. 0)                                                     &
           used = send_data(bling%id_jp_uptake,      bling%jp_uptake*rho_dzt,      &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_jpo4 .gt. 0)                                                          &
           used = send_data(bling%id_jpo4,           bling%jpo4*rho_dzt,           &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_jpop .gt. 0)                                                          &
           used = send_data(bling%id_jpop,           bling%jpop*rho_dzt,           &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%fe_is_prognostic) then
        if (bling%id_jfeop .gt. 0)                                                         &
             used = send_data(bling%id_jfeop,          bling%jfeop*rho_dzt,          &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
        if (bling%id_kfe_eq_lig .gt. 0)                                                    &
             used = send_data(bling%id_kfe_eq_lig,     bling%kfe_eq_lig,             &
             model_time, rmask = grid_tmask,                                                     & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (bling%id_pc_m .gt. 0)                                                          &
           used = send_data(bling%id_pc_m,           bling%pc_m,                   &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_mu .gt. 0)                                                            &
           used = send_data(bling%id_mu,             bling%mu,                     &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_o2_saturation .gt. 0)                                                 &
           used = send_data(bling%id_o2_saturation,  bling%o2_saturation,          &
           model_time, rmask = grid_tmask(:,:,1),                                              & 
           is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      if (bling%id_theta .gt. 0)                                                         &
           used = send_data(bling%id_theta,          bling%theta,                  &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_thetamax_fe .gt. 0)                                                   &
           used = send_data(bling%id_thetamax_fe,    bling%thetamax_fe,            &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_wsink .gt. 0)                                                         &
           used = send_data(bling%id_wsink,          bling%wsink,                  &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (bling%id_zremin .gt. 0)                                                        &
           used = send_data(bling%id_zremin,         bling%zremin,                 &
           model_time, rmask = grid_tmask,                                                     & 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (.not. bling%fe_is_prognostic) then
        if (bling%id_fed_data .gt. 0)                                                    &
             used = send_data(bling%id_fed_data,         bling%f_fed_data,         &
             model_time, rmask = grid_tmask,                                                   & 
             is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      endif

    enddo  !} n
 
  deallocate(tmp_irr_band)

    return

  end subroutine generic_miniBLING_update_from_source




!#######################################################################
  ! <SUBROUTINE NAME="generic_miniBLING_set_boundary_values">
  !  <DESCRIPTION>
  !   Calculate and set coupler values at the surface / bottom of the ocean.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_miniBLING_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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

  subroutine generic_miniBLING_set_boundary_values(tracer_list, SST, SSS, rho, ilb, jlb, tau)

    type(g_tracer_type), pointer,   intent(inout)               :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in)                  :: SST
    real, dimension(ilb:,jlb:),     intent(in)                  :: SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in)                  :: rho
    integer,                        intent(in)                  :: ilb
    integer,                        intent(in)                  :: jlb
    integer,                        intent(in)                  :: tau

    integer                                     :: isc
    integer                                     :: iec
    integer                                     :: jsc
    integer                                     :: jec
    integer                                     :: isd
    integer                                     :: ied
    integer                                     :: jsd
    integer                                     :: jed
    integer                                     :: nk
    integer                                     :: ntau 
    integer                                     :: i
    integer                                     :: j
    integer                                     :: n
    real                                        :: sal
    real                                        :: ST
    real                                        :: sc_co2
    real                                        :: sc_o2
    !real                                        :: sc_no_term
    real                                        :: o2_saturation
    real                                        :: tt
    real                                        :: tk
    real                                        :: ts
    real                                        :: ts2
    real                                        :: ts3
    real                                        :: ts4
    real                                        :: ts5
    real, dimension(:,:,:),   pointer           :: grid_tmask
    real, dimension(:,:,:,:), pointer           :: o2_field
    real, dimension(:,:),     pointer           :: co2_alpha
    real, dimension(:,:),     pointer           :: co2_csurf
    real, dimension(:,:),     pointer           :: co2_schmidt
    real, dimension(:,:),     pointer           :: o2_alpha
    real, dimension(:,:),     pointer           :: o2_csurf
    real, dimension(:,:),     pointer           :: o2_schmidt
    real, dimension(:,:),     pointer           :: co2_sat_rate
    real, dimension(:,:),     pointer           :: c14o2_alpha
    real, dimension(:,:),     pointer           :: c14o2_csurf
    real, dimension(:,:),     pointer           :: c14o2_schmidt
    real                                        :: surface_rho

    character(len=fm_string_len), parameter :: sub_name = 'generic_miniBLING_set_boundary_values'

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

    !Get the necessary properties
    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, grid_tmask = grid_tmask) 

    do n = 1, num_instances  

      call g_tracer_get_pointer(tracer_list, 'o2' // bling%suffix ,'field',   o2_field)
      call g_tracer_get_pointer(tracer_list, 'o2' // bling%suffix, 'alpha',   o2_alpha)
      call g_tracer_get_pointer(tracer_list, 'o2' // bling%suffix, 'csurf',   o2_csurf)
      call g_tracer_get_pointer(tracer_list, 'o2' // bling%suffix, 'sc_no', o2_schmidt)

      do j = jsc, jec  
        do i = isc, iec  

          sal = SSS(i,j)
          ST = SST(i,j)
 
          surface_rho = bling%Rho_0

       !---------------------------------------------------------------------
       !     O2
       !---------------------------------------------------------------------
       !  Compute the oxygen saturation concentration at 1 atm total pressure in mol/kg 
       !  given the temperature (T, in deg C) and the salinity (S, in permil).
       !
       !  From Garcia and Gordon (1992), Limnology and Oceonography (page 1310, eq (8)).
       !  *** Note: the "a3*ts^2" term was erroneous, and not included here. ***
       !  Defined between T(freezing) <= T <= 40 deg C and 0 <= S <= 42 permil.
       !
       ! check value: T = 10 deg C, S = 35 permil, o2_saturation = 0.282015 mol m-3
       !---------------------------------------------------------------------

          tt = 298.15 - ST
          tk = 273.15 + ST
          ts = log(tt / tk)
          ts2 = ts  * ts
          ts3 = ts2 * ts
          ts4 = ts3 * ts
          ts5 = ts4 * ts

          o2_saturation = (1000.0/22391.6) * grid_tmask(i,j,1) *  &         !convert from ml/l to mol m-3
               exp(bling%a_0 + bling%a_1*ts  + bling%a_2*ts2 + bling%a_3*ts3 +  &
               bling%a_4*ts4 + bling%a_5*ts5 + (bling%b_0    + bling%b_1*ts  +  &
               bling%b_2*ts2 + bling%b_3*ts3 + bling%c_0 * sal) * sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the formulation proposed
       !  by Keeling et al. (1998, Global Biogeochem. Cycles, 12, 141-163).
       !---------------------------------------------------------------------

          sc_o2  = bling%a1_o2  + ST * (bling%a2_o2  + ST * (bling%a3_o2  + &
               ST * bling%a4_o2 )) * grid_tmask(i,j,1)
       
          o2_alpha(i,j)   = o2_saturation
          bling%o2_saturation(i,j) = o2_saturation / surface_rho
          o2_csurf(i,j)   = o2_field(i,j,1,tau) * surface_rho
          o2_schmidt(i,j) = sc_o2

        enddo  !} i
      enddo  !} j

      if (bling%do_carbon) then                              !<<CARBON CYCLE
    
        call g_tracer_get_pointer(tracer_list, 'dic' // bling%suffix, 'alpha',   co2_alpha)
        call g_tracer_get_pointer(tracer_list, 'dic' // bling%suffix, 'csurf',   co2_csurf)
        call g_tracer_get_pointer(tracer_list, 'dic' // bling%suffix, 'sc_no', co2_schmidt)

        do j = jsc, jec  
          do i = isc, iec  

            ST = SST(i,j)
            surface_rho = bling%Rho_0

       !---------------------------------------------------------------------
       !     CO2
       !---------------------------------------------------------------------
       !---------------------------------------------------------------------
       !  Compute the Schmidt number of CO2 in seawater using the formulation
       !   presented by Wanninkhof (1992, J. Geophys. Res., 97, 7373-7382).
       !---------------------------------------------------------------------
            sc_co2 = bling%a1_co2 + ST * (bling%a2_co2 + ST * &
                (bling%a3_co2 + ST * bling%a4_co2)) * grid_tmask(i,j,1)
         
            co2_alpha(i,j)   = co2_alpha(i,j) * surface_rho
            co2_csurf(i,j)   = co2_csurf(i,j) * surface_rho
            co2_schmidt(i,j) = sc_co2

          enddo  !} i
        enddo  !} j

        if (bling%do_14c) then  
    
          call g_tracer_get_pointer(tracer_list, 'di14c' // bling%suffix, 'alpha',   c14o2_alpha)
          call g_tracer_get_pointer(tracer_list, 'di14c' // bling%suffix, 'csurf',   c14o2_csurf)
          call g_tracer_get_pointer(tracer_list, 'di14c' // bling%suffix, 'sc_no', c14o2_schmidt)

          do j = jsc, jec  
            do i = isc, iec  

              ST = SST(i,j)
              surface_rho = bling%Rho_0
              sc_co2 = bling%a1_co2 + ST * (bling%a2_co2 + ST * (bling%a3_co2 +       &
                   ST * bling%a4_co2)) * grid_tmask(i,j,1)
         
              c14o2_alpha(i,j)   = c14o2_alpha(i,j) * surface_rho 
              c14o2_csurf(i,j)   = c14o2_csurf(i,j) * surface_rho
              c14o2_schmidt(i,j) = sc_co2

            enddo  
          enddo 
        endif 
      endif                                                              !CARBON CYCLE>>

    enddo  

    return

  end subroutine generic_miniBLING_set_boundary_values




!#######################################################################
  ! <SUBROUTINE NAME="generic_miniBLING_end">
  !  <DESCRIPTION>
  !   End the module. Deallocate all work arrays.
  !  </DESCRIPTION>
  ! </SUBROUTINE>

  subroutine generic_miniBLING_end

    character(len=fm_string_len), parameter :: sub_name = 'generic_miniBLING_end'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter   :: warn_header =                                &
         '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

    integer :: stdout_unit

    stdout_unit = stdout()

    call user_deallocate_arrays

    return
  end subroutine generic_miniBLING_end


!#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Allocate all the work arrays to be used in this module.
  !

  subroutine user_allocate_arrays

    integer     :: isc
    integer     :: iec
    integer     :: jsc
    integer     :: jec
    integer     :: isd
    integer     :: ied
    integer     :: jsd
    integer     :: jed
    integer     :: nk
    integer     :: ntau
    integer     :: n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    !Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc
    CO2_dope_vec%iec = iec 
    CO2_dope_vec%jsc = jsc
    CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd
    CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd
    CO2_dope_vec%jed = jed    

  do n = 1, num_instances  

    allocate(bling%wrk_3d          (isd:ied, jsd:jed, 1:nk)); bling%wrk_3d=0.0
    allocate(bling%wrk_2d          (isd:ied, jsd:jed)      ); bling%wrk_2d=0.0
    allocate(bling%flux            (isd:ied, jsd:jed)      ); bling%flux=0.0
    allocate(bling%integral        (isd:ied, jsd:jed)      ); bling%integral=0.0
    allocate(bling%k_lev           (isd:ied, jsd:jed)      ); bling%k_lev=0.0

    if (bling%biomass_type .eq. 'single') then
      allocate(bling%biomass_p_ts    (isd:ied, jsd:jed, 1:nk)); bling%biomass_p_ts=0.0
    endif
    allocate(bling%def_fe          (isd:ied, jsd:jed, 1:nk)); bling%def_fe=0.0
    allocate(bling%expkT           (isd:ied, jsd:jed, 1:nk)); bling%expkT=0.0
    allocate(bling%f_chl           (isd:ied, jsd:jed, 1:nk)); bling%f_chl=0.0
    allocate(bling%f_fed           (isd:ied, jsd:jed, 1:nk)); bling%f_fed=0.0
    if (.not. bling%fe_is_prognostic) then
      allocate(bling%f_fed_data      (isc:iec, jsc:jec, 1:nk)); bling%f_fed_data=0.0
    endif
    allocate(bling%f_o2            (isd:ied, jsd:jed, 1:nk)); bling%f_o2=0.0
    allocate(bling%f_po4           (isd:ied, jsd:jed, 1:nk)); bling%f_po4=0.0
    allocate(bling%fe_2_p_uptake   (isd:ied, jsd:jed, 1:nk)); bling%fe_2_p_uptake=0.0
    if (bling%fe_is_prognostic) then
      allocate(bling%feprime         (isd:ied, jsd:jed, 1:nk)); bling%feprime=0.0
      allocate(bling%fpofe           (isd:ied, jsd:jed, 1:nk)); bling%fpofe=0.0
    endif
    allocate(bling%fpop            (isd:ied, jsd:jed, 1:nk)); bling%fpop=0.0
    allocate(bling%frac_lg         (isd:ied, jsd:jed, 1:nk)); bling%frac_lg=0.0
    allocate(bling%frac_pop        (isd:ied, jsd:jed, 1:nk)); bling%frac_pop=0.0
    allocate(bling%irr_inst        (isd:ied, jsd:jed, 1:nk)); bling%irr_inst=0.0
    allocate(bling%irr_mix         (isd:ied, jsd:jed, 1:nk)); bling%irr_mix=0.0
    allocate(bling%irrk            (isd:ied, jsd:jed, 1:nk)); bling%irrk=0.0
    if (bling%fe_is_prognostic) then
      allocate(bling%jfe_ads_inorg   (isd:ied, jsd:jed, 1:nk)); bling%jfe_ads_inorg=0.0
      allocate(bling%jfe_ads_org     (isd:ied, jsd:jed, 1:nk)); bling%jfe_ads_org=0.0
      allocate(bling%jfe_recycle     (isd:ied, jsd:jed, 1:nk)); bling%jfe_recycle=0.0
    endif
    if (bling%fe_is_prognostic .or. bling%fe_is_diagnostic) then
      allocate(bling%jfe_reminp      (isd:ied, jsd:jed, 1:nk)); bling%jfe_reminp=0.0
    endif
    allocate(bling%jfe_uptake      (isd:ied, jsd:jed, 1:nk)); bling%jfe_uptake=0.0
    allocate(bling%jo2             (isd:ied, jsd:jed, 1:nk)); bling%jo2=0.0
    allocate(bling%jp_recycle      (isd:ied, jsd:jed, 1:nk)); bling%jp_recycle=0.0
    allocate(bling%jp_reminp       (isd:ied, jsd:jed, 1:nk)); bling%jp_reminp=0.0
    allocate(bling%jp_uptake       (isd:ied, jsd:jed, 1:nk)); bling%jp_uptake=0.0
    allocate(bling%jpo4            (isd:ied, jsd:jed, 1:nk)); bling%jpo4=0.0
    allocate(bling%jpop            (isd:ied, jsd:jed, 1:nk)); bling%jpop=0.0
    if (bling%fe_is_prognostic) then
      allocate(bling%jfeop           (isd:ied, jsd:jed, 1:nk)); bling%jfeop=0.0
      allocate(bling%kfe_eq_lig      (isd:ied, jsd:jed, 1:nk)); bling%kfe_eq_lig=0.0
    endif
    allocate(bling%mu              (isd:ied, jsd:jed, 1:nk)); bling%mu=0.0
    allocate(bling%pc_m            (isd:ied, jsd:jed, 1:nk)); bling%pc_m=0.0
    allocate(bling%theta           (isd:ied, jsd:jed, 1:nk)); bling%theta=0.0
    allocate(bling%thetamax_fe     (isd:ied, jsd:jed, 1:nk)); bling%thetamax_fe=0.0
    allocate(bling%wsink           (isd:ied, jsd:jed, 1:nk)); bling%wsink=0.0
    allocate(bling%zremin          (isd:ied, jsd:jed, 1:nk)); bling%zremin=0.0
    allocate(bling%zbot            (isd:ied, jsd:jed, 1:nk)); bling%zbot=0.0
    allocate(bling%b_o2            (isd:ied, jsd:jed));       bling%b_o2=0.0
    allocate(bling%b_po4           (isd:ied, jsd:jed));       bling%b_po4=0.0
    if (bling%fe_is_prognostic) then
      allocate(bling%b_fed           (isd:ied, jsd:jed));       bling%b_fed=0.0
      allocate(bling%fe_burial       (isd:ied, jsd:jed));       bling%fe_burial=0.0
      allocate(bling%ffe_sed         (isd:ied, jsd:jed));       bling%ffe_sed=0.0
    endif
    allocate(bling%o2_saturation   (isd:ied, jsd:jed));       bling%o2_saturation=0.0

  if (bling%do_carbon) then                                     !<<CARBON CYCLE
    allocate(bling%b_dic           (isd:ied, jsd:jed));       bling%b_dic=0.0
    allocate(bling%co2_alpha       (isd:ied, jsd:jed));       bling%co2_alpha=0.0
    allocate(bling%co2_csurf       (isd:ied, jsd:jed));       bling%co2_csurf=0.0
    allocate(bling%htotallo        (isd:ied, jsd:jed))
    allocate(bling%htotalhi        (isd:ied, jsd:jed))
    allocate(bling%pco2_surf       (isd:ied, jsd:jed));       bling%pco2_surf=0.0
    allocate(bling%surf_temp       (isd:ied, jsd:jed));       bling%surf_temp=0.0
    allocate(bling%surf_salt       (isd:ied, jsd:jed));       bling%surf_salt=0.0
    allocate(bling%surf_alk        (isd:ied, jsd:jed));       bling%surf_alk=0.0
    allocate(bling%surf_po4        (isd:ied, jsd:jed));       bling%surf_po4=0.0
    allocate(bling%surf_sio4       (isd:ied, jsd:jed));       bling%surf_sio4=0.0
    allocate(bling%surf_dic        (isd:ied, jsd:jed));       bling%surf_dic=0.0
  if (bling%do_14c) then                                        !<<RADIOCARBON
    allocate(bling%c14_2_p         (isd:ied, jsd:jed, 1:nk)); bling%c14_2_p=0.0
    allocate(bling%fpo14c          (isd:ied, jsd:jed, 1:nk)); bling%fpo14c=0.0
    allocate(bling%j14c_decay_dic  (isd:ied, jsd:jed, 1:nk)); bling%j14c_decay_dic=0.0
    allocate(bling%j14c_reminp     (isd:ied, jsd:jed, 1:nk)); bling%j14c_reminp=0.0
    allocate(bling%jdi14c          (isd:ied, jsd:jed, 1:nk)); bling%jdi14c=0.0
    allocate(bling%b_di14c         (isd:ied, jsd:jed));       bling%b_di14c=0.0
    allocate(bling%c14o2_alpha     (isd:ied, jsd:jed));       bling%c14o2_alpha=0.0
    allocate(bling%c14o2_csurf     (isd:ied, jsd:jed));       bling%c14o2_csurf=0.0
  endif                                                               !RADIOCARBON>>
  endif                                                               !CARBON CYCLE>>

  enddo  

    return

  end subroutine user_allocate_arrays



!#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays

    integer     :: n

    do n = 1, num_instances  

      deallocate(bling%wrk_3d)
      deallocate(bling%wrk_2d)
      deallocate(bling%flux)
      deallocate(bling%integral)
      deallocate(bling%k_lev)

      deallocate(bling%o2_saturation)
      if (bling%biomass_type .eq. 'single') then
        deallocate(bling%biomass_p_ts)
      endif
      deallocate(bling%def_fe)
      deallocate(bling%expkT)
      deallocate(bling%f_chl)
      deallocate(bling%f_fed)
      if (.not. bling%fe_is_prognostic) then
        deallocate(bling%f_fed_data)
      endif
      deallocate(bling%f_o2)
      deallocate(bling%f_po4)
      deallocate(bling%fe_2_p_uptake)
      if (bling%fe_is_prognostic) then
        deallocate(bling%feprime)
        deallocate(bling%fpofe)
      endif
      deallocate(bling%fpop)
      deallocate(bling%frac_lg)
      deallocate(bling%frac_pop)
      deallocate(bling%irr_inst)
      deallocate(bling%irr_mix)
      deallocate(bling%irrk)
      if (bling%fe_is_prognostic) then
        deallocate(bling%jfe_ads_inorg)
        deallocate(bling%jfe_ads_org)
        deallocate(bling%jfe_recycle)
      endif
      if (bling%fe_is_prognostic .or. bling%fe_is_diagnostic) then
        deallocate(bling%jfe_reminp)
      endif
      deallocate(bling%jfe_uptake)
      deallocate(bling%jo2)
      deallocate(bling%jp_recycle)
      deallocate(bling%jp_reminp)
      deallocate(bling%jp_uptake)
      deallocate(bling%jpo4)
      deallocate(bling%jpop)
      if (bling%fe_is_prognostic) then
        deallocate(bling%jfeop)
        deallocate(bling%kfe_eq_lig)
      endif
      deallocate(bling%pc_m)
      deallocate(bling%mu)
      deallocate(bling%theta)
      deallocate(bling%thetamax_fe)
      deallocate(bling%wsink)
      deallocate(bling%zremin)
      deallocate(bling%zbot)
      if (bling%fe_is_prognostic) then
        deallocate(bling%fe_burial)
        deallocate(bling%ffe_sed)
        deallocate(bling%b_fed)
      endif
      deallocate(bling%b_o2)
      deallocate(bling%b_po4)

      if (bling%do_carbon) then                                    !<<CARBON CYCLE
        deallocate(bling%surf_temp)
        deallocate(bling%surf_salt)
        deallocate(bling%surf_alk)
        deallocate(bling%surf_po4)
        deallocate(bling%surf_sio4)
        deallocate(bling%surf_dic)
        deallocate(bling%co2_csurf)
        deallocate(bling%pco2_surf)
        deallocate(bling%co2_alpha)
        deallocate(bling%htotallo)
        deallocate(bling%htotalhi)
        deallocate(bling%b_dic)
        if (bling%do_14c) then                                      !<<RADIOCARBON
          deallocate(bling%c14_2_p)
          deallocate(bling%fpo14c)
          deallocate(bling%j14c_decay_dic)
          deallocate(bling%j14c_reminp)
          deallocate(bling%jdi14c)
          deallocate(bling%c14o2_alpha)
          deallocate(bling%c14o2_csurf)
          deallocate(bling%b_di14c)
        endif   !}                                              !RADIOCARBON>>
      endif  !}                                                !CARBON CYCLE>>

    enddo  !} n

    return

  end subroutine user_deallocate_arrays


end module generic_miniBLING_mod
