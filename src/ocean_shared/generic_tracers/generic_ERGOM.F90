!----------------------------------------------------------------
! <CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> 
! </CONTACT>
! 
! <REVIEWER EMAIL="none yet"> 
! </REVIEWER>
!<OVERVIEW>
! This module contains the generic version of ERGOM modified for the project GENUS.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
! The genreal coding scheme follows that of the TOPAZ package
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the ERGOM equations
!       Fennel and Neumann
!       S. Schaefer
!       M. Schmidt, A Eggert, H. Radtke
!       Some code pieces are reused from the TOPAZ code - 
!       no need to invent the wheel twice. Thanks to John Dunne and Niki Zadeh.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.io-warnemuende.de
! </REFERENCE>
! <DEVELOPER_NOTES>
! very new
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_ERGOM

  use coupler_types_mod,   only: coupler_2d_bc_type
  use field_manager_mod,   only: fm_string_len
  use fms_mod,             only: write_version_number, open_namelist_file, close_file, check_nml_error  
  use mpp_mod,             only: mpp_error, NOTE, WARNING, FATAL, stdout, stdlog, input_nml_file
  use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
  use time_manager_mod,    only: time_type
  use fm_util_mod,         only: fm_util_start_namelist, fm_util_end_namelist  
  use diag_manager_mod,    only: register_diag_field, send_data 
  use constants_mod,       only: WTMO2

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  
  use g_tracer_utils, only : g_diag_type, g_diag_field_add


  implicit none ; private

  character(len=128) :: version = '$Id: generic_ERGOM.F90,v 20.0 2013/12/14 00:18:05 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  character(len=fm_string_len), parameter :: mod_name       = 'generic_ERGOM'
  character(len=fm_string_len), parameter :: package_name   = 'generic_ergom'

  public do_generic_ERGOM
  public generic_ERGOM_register
  public generic_ERGOM_init
  public generic_ERGOM_update_from_coupler
  public generic_ERGOM_update_from_source
  public generic_ERGOM_update_from_bottom
  public generic_ERGOM_set_boundary_values
  public generic_ERGOM_register_diag
  public generic_ERGOM_end

  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  !An auxiliary type for storing varible names
  type, public :: vardesc
     character(len=fm_string_len) :: name     ! variable name in a NetCDF file.
     character(len=fm_string_len) :: longname ! long name of that variable.
     character(len=1)  :: hor_grid            ! hor. grid:  u, v, h, q, or 1.
     character(len=1)  :: z_grid              ! vert. grid:  L, i, or 1.
     character(len=1)  :: t_grid              ! time description: s, a, m, or 1.
     character(len=fm_string_len) :: units    ! dimensions of the variable.
     character(len=1)  :: mem_size            ! size in memory: d or f.
  end type vardesc
  
  type generic_ERGOM_type
     real :: t1_nit, t2_nit, t3_nit, t4_nit
     real :: t0_o2, t1_o2, t2_o2, t3_o2, t4_o2, t5_o2  ! for oxygen solubility
     real :: b0_o2, b1_o2, b2_o2, b3_o2, c0_o2         ! for oxygen solubility
     real :: e1_nit
     ! dimensionless parameters
     real :: s1_nit, s2_nit, s3_nit
     ! dimension is 1/salinity
     real :: a1_nit, a2_nit, a3_nit, a4_nit  !  coefficients for Schmidt number
     real :: a1_o2 , a2_o2 , a3_o2 , a4_o2   !  coefficients for Schmidt number
     real :: Rho_0
     ! 
     real :: &
          q10_nit    = 0., &   ! q10 parameter for nitrification [1/Celsius]
          q10_h2s    = 0., &   ! q10 parameter for chemolithotrophs (so4 reduction) [1/Celsius]
          nf         = 0., &   ! nitrification rate [1/s]
          alpha_nit  = 0., &   ! half-saturation constant for nitrification [mol/kg]
          alp_o2     = 0., &   ! slope function for detritus recycling [kg/mol]
	  alp_no3    = 0., &   ! slope function for detritus recycling [kg/mol]
	  alp_h2s    = 0., &   ! slope function for detritus recycling [kg/mol]
	  alp_nh4    = 0., &   ! slope function for detritus recycling [kg/mol]
	  k_h2s_o2   = 0., &   ! reaction constant h2s oxidation with o2 [kg/mol/s]
          k_h2s_no3  = 0., &   ! reaction constant h2s oxidation with no3 [kg/mol/s]
          k_sul_o2   = 0., &   ! reaction constant sulfur oxidation with o2 [kg/mol/s]
          k_sul_no3  = 0., &   ! reaction constant sulfur oxidation with no3 [kg/mol/s]
          k_an0      = 0.      ! maximum anammox rate [1/s]
     real, dimension(:,:,:), ALLOCATABLE ::  &
          f_no3,   &           !no3 concentration [mol/kg]
          f_nh4,   &           !nh4 concentration [mol/kg]
          f_po4,   &           !po4 concentration [mol/kg]
          f_o2 ,   &           !o2 concentration [mol/kg]
          f_h2s,   &           !h2s concentration [mol/kg]
          f_sul,   &           !sulfur concentration [mol/kg]	       
	  f_chl,   &           !chlorophyll concentration [µg/kg] 
	  irr_inst,&           !instantaneous light [W/m2]
          jno3 ,   &           !time change of no3 concentration [mol/kg/s]
          jnh4 ,   &           !time change of nh4 concentration [mol/kg/s]
          jpo4 ,   &           !time change of po4 concentration [mol/kg/s]
          jo2  ,   &           !time change of o2 concentration [mol/kg/s]
          jh2s ,   &           !time change of h2s concentration [mol/kg/s]
          jsul ,   &           !time change of h2s concentration [mol/kg/s]
          jh2s_o2   , &        !time change of h2s concentration [mol/kg/s] by oxygen
          jh2s_no3  , &        !time change of h2s concentration [mol/kg/s] by nitrate
          jsul_o2   , &        !time change of sul concentration [mol/kg/s] by oxygen
          jsul_no3  , &        !time change of sul concentration [mol/kg/s] by nitrate
	  jrec_o2   , &        !nitrogen loss to nh4 by recycling with o2 [mol/kg/s]
	  jrec_no3  , &        !nitrogen loss to nh4 by recycling with no3 [mol/kg/s]
	  jrec_so4  , &        !nitrogen loss to nh4 by recycling with so4 [mol/kg/s]
	  jrec_ana  , &        !nitrogen loss to nh4 by recycling and subseqent anammox [mol/kg/s]
	  jdenit_wc,  &        !denitrification in water column [mol/kg/s]
	  jnitrif              !nitrification in water column [mol/kg/s]
     real, dimension(:,:), ALLOCATABLE :: &
          b_nh4,      &
	  b_no3,      &
	  b_o2,       &
	  b_po4,      &
	  b_nitrogen, &
	  b_h2s
     real, dimension(:,:,:,:), pointer :: &
          p_no3,      &        !no3 concentration [mol/kg]
          p_nh4,      &        !nh4 concentration [mol/kg]
          p_po4,      &        !po4 concentration [mol/kg] 
          p_o2,       &        !o2 concentration [mol/kg]
          p_h2s,      &        !h2s concentration [mol/kg]
          p_sul,      &        !sulfur concentration [mol/kg]
          p_nitrogen           !n2 [mol/kg]
     integer            ::    &
          id_no3        = -1, & ! no3 prognostic tracer
          id_nh4        = -1, & ! nh4 prognostic tracer
          id_po4        = -1, & ! po4 prognostic tracer
          id_o2         = -1, & ! o2 prognostic tracer
          id_h2s        = -1, & ! h2s prognostic tracer
          id_chl        = -1, & ! chlorophyll diagnostic tracer
          id_sul        = -1, & ! sulfur prognostic tracer
          id_dia        = -1, & ! nitrogen in diatoms prognostic tracer
          id_fla        = -1, & ! nitrogen in flagellates prognostic tracer
          id_cya        = -1, & ! nitrogen in cyanobacteria prognostic tracer
          id_zoo        = -1, & ! nitrogen in zooplankton prognostic tracer
          id_det        = -1, & ! nitrogen in detritus prognostic tracer
	  id_irr_inst   = -1, & ! instantaneous light 
          id_jno3       = -1, & ! no3 source layer integral [mol/m2/s]
          id_jnh4       = -1, & ! nh4 source layer integral [mol/m2/s]
          id_jpo4       = -1, & ! po4 source layer integral [mol/m2/s]
          id_jo2        = -1, & ! o2 source layer integral [mol/m2/s]
          id_jh2s_o2    = -1, & ! h2s loss with oxygen source layer integral [mol/m2/s]
          id_jh2s_no3   = -1, & ! h2s loss with nitrate source layer integral [mol/m2/s]
          id_jsul_o2    = -1, & ! sulfur loss with oxygen source layer integral [mol/m2/s]
          id_jsul_no3   = -1, & ! sulfur loss with nitrate source layer integral [mol/m2/s]
          id_jh2s       = -1, & ! h2s source layer integral [mol/m2/s]
          id_jsul       = -1, & ! sulfur source layer integral [mol/m2/s]
          id_jrec_o2    = -1, & ! id for layer integral nitrogen loss to nh4 by recycling with o2 [mol/m2/s]
          id_jrec_no3   = -1, & ! id for layer integral nitrogen loss to nh4 by recycling with no3 [mol/m2/s]
          id_jrec_so4   = -1, & ! id for layer integral nitrogen loss to nh4 by recycling with so4 [mol/m2/s]
          id_jrec_ana   = -1, & ! id for layer integral nitrogen loss to nh4 by recycling and subseqent anammox [mol/m2/s]
          id_jnitrif    = -1, & ! nitrification [mol/m2/s]
          id_jdenit_wc  = -1, & ! denitrification in water column [mol/m2/s]  
          id_dep_dry_po4  = -1,  & ! phosphate dry deposition
          id_dep_dry_nh4  = -1,  & ! Ammonia dry deposition
          id_dep_dry_no3  = -1,  & ! Nitrate dry deposition
          id_dep_wet_po4  = -1,  & ! phosphate wet deposition
          id_dep_wet_nh4  = -1,  & ! Ammonia wet deposition
          id_dep_wet_no3  = -1,  & ! Nitrate wet deposition
          id_runoff_flux_po4 = -1,  & ! phosphate runoff flux to the ocean
          id_runoff_flux_nh4 = -1,  & ! Ammonia runoff Flux to the ocean
          id_runoff_flux_no3 = -1     ! Nitrate runoff Flux to the ocean
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file, IC_file
  end type generic_ERGOM_type

  type phytoplankton
     real, ALLOCATABLE, dimension(:,:,:)  :: &
          f_n         , & ! nitrogen in phytoplankton, the intermediate value after physics [mol/kg]
          ilim        , & ! light limitation factor [dimensionless]
          jprod_no3   , & ! no3 uptake [mol/kg/s]
	  jprod_nh4   , & ! nh4 uptake [mol/kg/s]
          jprod_po4   , & ! po4 uptake [mol/kg/s]
          jprod_n2    , & ! n2 fixation [mol/kg/s]
	  jgraz_n     , & ! nitrogen loss by grazing [mol/kg/s]
	  jres_n      , & ! nitrogen loss by respiration [mol/kg/s]
	  jdet_n      , & ! nitrogen loss to detritus [mol/kg/s]
	  move            ! phytoplankton sinking velocity (<0 for sinking) [m/s]	  
     real, dimension(:,:,:,:), pointer :: &
          p_phyt          !nitrogen in phytoplankton [mol/kg]
     real :: &
          imin     = 0. , & ! minimum light [W/m2]
	  tmin     = 0. , & ! minimum temperature [Celsius]
	  smin     = 0. , & ! minimum phytoplankton salinity [g/kg]
	  smax     = 0. , & ! maximum phytoplankton salinity [g/kg]
          alpha    = 0. , & ! DIN half-saturation constant [mol/kg]
          talpha   = 0. , & ! Michaeles Menton-like temperature [Celsius]
          rp0      = 0. , & ! maximum uptake rate [1/s]
          p0       = 0. , & ! background concentration for initial growth [mol/kg]
	  pnr      = 0. , & ! P/N ratio
	  cnr      = 0. , & ! C/N ratio
	  lpd      = 0. , & ! phytoplankton loss to detritus [1/s]
	  lpr      = 0. , & ! phytoplankton loss by respiration [1/s]
	  wsink0   = 0.     ! surface phytoplankton sinking velocity [m/s]
     integer ::            &
          id_jprod_no3  = -1 , & ! Diag id for no3 production layer integral [mol/m2/s]
	  id_jprod_nh4  = -1 , & ! Diag id for nh4 production layer integral [mol/m2/s]
          id_jprod_po4  = -1 , & ! Diag id for po4 production layer integral [mol/m2/s]
          id_jprod_n2   = -1 , & ! Diag id for n2 fixation layer integral [mol/m2/s]
	  id_jgraz_n    = -1 , & ! Diag id for nitrogen grazing layer integral [mol/m2/s]
	  id_jres_n     = -1 , & ! Diag id for nitrogen respiration layer integral [mol/m2/s]
	  id_jdet_n     = -1 , & ! Diag id for nitrogen detritus layer integral [mol/m2/s]
          id_ilim       = -1 , & ! light limitation	     
          id_nlim       = -1 , & ! DIN limitation
	  id_plim       = -1     ! DIP limitation	   
     character(len=3) :: name	  
  end type phytoplankton 

  type zooplankton
     real, ALLOCATABLE, dimension(:)      :: pref_phy, pref_zoo, pref_det
     real, ALLOCATABLE, dimension(:,:,:)  :: &
          f_n       , & ! nitrogen in zooplankton, the intermediate value after physics [mol/kg]
	  jgraz_n   , & ! nitrogen loss by grazing [mol/kg/s]
          jgain_n   , & ! nitrogen gain by grazing [mol/kg/s]
	  jres_n    , & ! nitrogen loss by respiration [mol/kg/s]
	  jdet_n    , & ! nitrogen loss to detritus [mol/kg/s]
	  move      	! zooplankton movement [m/s]
     real, dimension(:,:,:,:), pointer :: &
          p_zoo           !nitrogen in zooplankton [mol/kg]
     real, dimension(:,:,:), pointer   :: &
          p_vmove     , & !vertical movement [m/s]
          p_vdiff         !vertical diffusion [m²/s]
     real :: &
	  pnr            = 0. , &  ! P/N ratio
	  cnr            = 0. , &  ! C/N ratio
	  t_opt          = 0. , &  ! optimal grazing temperature [Celsius]
          t_max          = 0. , &  ! maximal grazing temperature [Celsius]
          beta           = 0. , &  ! parameter for temperature dependence of grazing [dimensionless]
	  sigma_b        = 0. , &  ! zooplankton loss rate to detritus [1/s]
          oxy_sub        = 0. , &  ! oxygen level below which reduced respiration starts [mol/kg]
          oxy_min        = 0. , &  ! oxygen level below which no respiration takes place [mol/kg]
          resp_red       = 0. , &  ! reduction factor for respiration under suboxic conditions [dimensionless]
	  nue            = 0. , &  ! zooplankton loss rate to nh4 by respiration [1/s]
          food_to_nh4    = 0. , &  ! fraction of eaten food directly lost to respiration [dimensionless]
          food_to_det    = 0. , &  ! fraction of eaten food directly lost to detritus [dimensionless]
          food_to_nh4_2  = 0. , &  ! fraction of (food that could be eaten at optimal temperature) 
	                           ! directly lost to respiration [dimensionless]
          food_to_det_2  = 0. , &  ! fraction of (food that could be eaten at optimal temperature) 
	                           ! directly lost to detritus [dimensionless]
	  iv             = 0. , &  ! Ivlev constant [kg/mol]
	  zcl1           = 0. , &  ! closure parameter [kg/mol]
	  graz           = 0. , &  ! zooplankton maximum grazing rate [1/s]
          z0             = 0. , &  ! background concentration for initial growth [mol/kg]
	  Imax           = 0. , &  ! maximum light intensity [W/m^2]
          alpha          = 0. , &  ! light inhibition shape factor
          o2min          = 0. , &  ! minimum oxygen concentration where sinking stops [mol/kg]
          h2smax         = 0. , &  ! maximum h2s concentration where sinking stops [mol/kg]
          wtemp          = 0. , &  ! weight number for temperature sensitivity
          wo2            = 0. , &  ! weight number for o2 sensitivity
          wh2s           = 0. , &  ! weight number for h2s sensitivity
	  wsink0         = 0. , &  ! sink velocity (<0 for sinking) [m/s]
	  wrise0         = 0. , &  ! maximum rise velocity (>0 for rising) [m/s]
	  vdiff_max      = 0. , &  ! maximum enhanced diffusion
          dark_rise      = 0. , &  ! whether zooplankton rises in the dark independent off a food gradients [dimensionless]
	  wfood          = 0.      ! weight number for food gradiens 
     logical ::            &
          vertical_migration    = .false., & ! if true this special undergoes vertical migration
          blanchard_temperature = .false.    ! .false.: old ERGOM temperature dependence, 
	                                     ! .true. : Blanchard 1996 formula
     integer ::            &
	  graz_pref      = -1  , & ! flag to select grazing preferences
	  id_jgraz_n     = -1  , & ! Diag id for nitrogen grazing (loss) layer integral [mol/m2/s]
          id_jgain_n     = -1  , & ! Diag id for nitrogen grazing (gain) layer integral [mol/m2/s]
	  id_jres_n      = -1  , & ! Diag id for nitrogen respiration layer integral [mol/m2/s]
	  id_jdet_n      = -1  , & ! Diag id for nitrogen detritus layer integral [mol/m2/s]
	  id_vmove       = -1      ! Diag id for sink velocity 
     character(len=32) ::  &
          name             = 'none'  
  end type zooplankton

  type detritus
     ! The detritus variable is only a reference to a suspended matter 3d variable.
     ! It must have the same name.
     real, ALLOCATABLE, dimension(:,:,:)  :: &
          f_n      , &            ! nitrogen in detritus, the intermediate value after physics [mol/kg]
	  jgraz_n  , &            ! nitrogen loss to zooplankton by grazing [mol/kg/s]
	  jmort    	          ! nitrogen gain in detritus by mortality [mol/kg/s]
     real :: &
	  dn       = 0. , &       ! recycling rate [1/s]
	  q10_rec  = 0.	          ! q10 parameter for recycling of detritus [1/Celsius]
     integer ::        &
                                  ! detritus is suspended matter. The data field is in spm(i) 
          index_spm      = -1 , & ! stores the index index_spm in spm(index_spm)
          id_jgraz_n     = -1 , & ! Diag id for nitrogen loss by grazing layer integral [mol/m2/s]
          id_jmort       = -1     ! Diag id for nitrogen gain by mortality [mol/m2/s]
     character(len=32) ::  &
          name           = 'none' ! use the same name as a spm 3d variable
  end type detritus

  type sediment
     ! This type is for setting the biological rate parameters of the surface layer sediment.
     ! It contains no own data array.
     real, ALLOCATABLE, dimension(:,:)  :: &
          bioerosion  , &          ! local intensity of bioerosion (0.0 to 1.0)
	  jsed_n      , &          ! Nitrogen gain by sedimentation [mol/m2/s]
	  jrec_n      , &          ! Nitrogen loss to nh4 by recycling [mol/m2/s]
	  jdenit_sed  , &          ! Nitrogen loss to n2 by denitrification [mol/m2/s]
	  mode_sed                 ! support of bacterial matts 
     real :: &
	  dn              = 0. , & ! recycling rate [1/s]
          frac_dn_anoxic  = 0. , & ! fraction of recycling in shallow sediments for anoxic bottom water [dimensionless]
	  thio_bact_min   = 0. , & ! minimum amount of active sediment for thiomargarita [mol/m2]
	  q10_rec         = 0. , & ! q10 parameter for recycling [1/Celsius]
	  den_rate        = 0. , & ! proportion of denitrification at the redoxcline in sediment
	  pnr             = 0. , & ! P/N ratio
	  cnr             = 0. , & ! C/N ratio
	  wsed            = 0. , & ! sedimentation rate [m/s]
          po4_lib_rate    = 0. , & ! liberation rate for iron phosphate in the sediment [1/s]
          po4_retention   = 0. , & ! fraction of phosphorous retained in the sediment while recycled [dimensionless]
          po4_ret_plus_BB = 0. , & ! value added to po4_retention north of 60.75N 
                          	   ! (special treatment for the Bothnian Bay to supress cyanobacterial blooms)
          o2_bioerosion   = 0.     ! oxygen threshold for bioerosion  [mol/kg]
     character(len=32) ::  &
          name_redfield_sed   = 'sed' , &  ! name of the redfield-ratio sediment variable in sed
          name_iron_phosphate = 'ips'      ! name of the iron phosphate sediment variable in sed
     integer ::            &
          index_sed       = -1    , &  ! index of redfield-ratio sediment variable in sed
          index_ips       = -1    , &  ! index of iron phosphate sediment variable in sed
          id_jrec_n       = -1    , &  ! Diag id for recycling 
          id_jdenit_sed   = -1    , &  ! denitrification in water column [mol/m2/s]  
          id_mode_sed     = -1         ! support of bacterial matts 
  end type sediment

  type spm_type
     ! This is a type which describes a type of suspended particulate matter
     ! that is able to settle (such as detritus).
     ! You should specify "sediment_to" to allow sedimentation to a sed_type tracer.
     real, ALLOCATABLE, dimension(:,:,:)  :: &
	  move                          ! sinking velocity (<0 for sinking) [m/s]
     real, dimension(:,:,:,:), pointer :: &
          p_wat                         ! pointer to 3d variable for concentration in water column [mol/kg]
     real, ALLOCATABLE, dimension(:,:) :: &
          jsed                          ! sedimentation [mol/m2/s]
     real, ALLOCATABLE, dimension(:,:) :: &
          btf                           ! The total bottom flux [mol/m2/s]
     real :: &
	  wsink0          = 0. , &      ! sinking velocity (<0 for sinking) [m/d]
          wsed            = 0.          ! sedimentation rate [m/d]
     character(len=32) ::  &
          name          = 'none' , &    ! name of 3d tracer
          longname      = 'none' , &    ! long name for output
          sediment_to   = 'none'        ! to which sed_type 2d variable the sedimentation takes place, 
     integer ::               &
          index_sediment_to , &         ! Index i in the array sed(i) to which the sedimentation takes place
          id_jsed                       ! Diag id for sedimentation [mol/m2/s]
  end type spm_type
  
  type sed_type
     ! This is a type which describes a type of particles in the sediment 
     ! that is able to be resuspended (such as detritus).
     ! You should specify "suspend_to" to allow resuspension to a spm_type variable. 
     real, ALLOCATABLE, dimension(:,:,:) :: &
          f_sed                         ! f_sed(i, j, layer) -> 2d arrays for storing sediment concentration [mol/m2]
     real, ALLOCATABLE, dimension(:,:) :: &
          jgain_sed , &                 ! gain by transformation of other sediments [mol/m2/s]
          jloss_sed , &                 ! loss by transformation of sediment [mol/m2/s]
          jres      , &                 ! resuspension [mol/m2/s]
          jbiores                       ! bioerosion [mol/m2/s]
     real :: &
          erosion_rate    = 0. , &      ! erosion rate [1/d]
          bioerosion_rate = 0. , &      ! erosion rate by benthic animals [1/d]
          molar_volume    = 0. , &      ! volume of this tracer when deposited in the sediment [m3/mol]
          critical_stress = 160.        ! critical shear stress when erosion starts [N/m2]
     character(len=32) ::  &
          name          = 'none' , &    ! name of 2d tracer
          longname      = 'none' , &    ! long name for output
          suspend_to    = 'none'        ! to which spm_type 3d variable the resuspension takes place
     integer ::               &
          index_suspend_to  , &         ! Index i in the array spm(i) to which the resuspension takes place
          id_jgain_sed , &              ! Diag id for gain by transformation of other sediments [mol/m2/s]
          id_jloss_sed , &              ! Diag id for loss by transformation of sediment [mol/m2/s]
          id_jbiores   , &              ! Diag id for bio-resuspension [mol/m2/s]
          id_jres                       ! Diag id for resuspension [mol/m2/s]
  end type sed_type

  type sed_defs_type
    integer :: &
          NUM_LAYERS        = -1 , & ! Number of vertical layers
          layer_propagation = -1 , & ! Sediment layer propagation settings, 
	                             ! SLP_DOWNWARD=1, SLP_FULL_BOX=2, SLP_OLD_ERGOM=3
          erosion_mode      = -1     ! Sediment erosion mode, INDEPENDENT=1, MAXSTRESS=2, ORGANIC=3
    real, allocatable, dimension(:) :: layer_height  ! (maximum) height of vertical layers [m]. 
                                                     ! <0 means the layer may become infinitely thick.
  end type sed_defs_type

  type tracer_2d
    character(len=fm_string_len) ::  name               = 'none'  ! name of the 2d tracer
    character(len=fm_string_len) ::  longname           = 'none'  ! longname of the 2d tracer
    character(len=fm_string_len) ::  units              = 'none'  ! units of the 2d tracer
    character(len=fm_string_len) ::  name_of_3d_tracer  = 'none'  ! name of the 3d tracer
    integer                      ::  layer_in_3d_tracer = 0       ! z-level inside the 3d tracer
    integer                      ::  diag_field         = -1      ! handler returned from register_diag_field
    real, dimension(:,:), pointer::  p_field                      ! pointer to the 2d field stored
    logical                      ::  field_assigned     = .false. ! whether a 2d field has been assigned
  end type tracer_2d
  
  ! Zooplankton preference settings
  integer, parameter             :: ERGOM_PREFS=1, HUTSON=2, HUTSON_QUADRATIC=3, &
                                    GENUS_IVLEV=4, GENUS_IVLEV_QUADRATIC=5
  ! Sediment layer propagation settings
  integer, parameter             :: SLP_DOWNWARD=1, SLP_FULL_BOX=2, SLP_OLD_ERGOM=3
  ! Sediment erosion settings
  integer, parameter             :: INDEPENDENT=1, MAXSTRESS=2, ORGANIC=3

  type(generic_ERGOM_type), save   :: ergom
  integer :: id_def_nit

  !
  ! Array allocations and flux calculations assume that phyto(1) is the
  ! only phytoplankton group cabable of nitrogen uptake by N2 fixation 
  ! while phyto(2:NUM_PHYTO) 
  ! are only cabable of nitrogen uptake by NH4 and NO3 uptake
  !
  type(phytoplankton),        ALLOCATABLE, dimension(:), save :: phyto
  
  type(zooplankton),          ALLOCATABLE, dimension(:), save :: zoo
 
  type(detritus),             ALLOCATABLE, dimension(:), save :: det

  type(sediment), save                                        :: biosed

  type(sed_defs_type),                                   save :: sed_defs

  type(spm_type),         ALLOCATABLE, dimension(:), save :: spm
  type(sed_type),         ALLOCATABLE, dimension(:), save :: sed

  type(tracer_2d),            ALLOCATABLE, dimension(:), save :: tracers_2d
   
  real, parameter :: epsln=1.0e-30

! Most ecosystem parameters use 1/d as time unit. This is to convert  
  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: watt_2_my_einstein = 4.6
  real, parameter :: missing_value1=-1.0e+20
    
  !The following logical for using this module is overwritten 
  ! by generic_tracer_nml namelist
  logical, save :: do_generic_ERGOM = .false.
  
  integer, parameter :: maxphyt=3, maxzoo=4, maxdet=2 
  integer, parameter :: maxspm=2, maxsed=2, max_sediment_layers=2
  real, dimension(maxphyt) :: &
        imin   , &	  ! minimum light [W/m2]
        tmin   , &	    ! minimum phytoplankton growth temperature (for cyanobacteria only) [Celsius]
        smin   , &	    ! minimum phytoplankton salinity (for diatoms only) 
        smax   , &	    ! maximum phytoplankton salinity (for diatoms only) 
        alpha  , &	    ! half-saturation constants of nutrient uptake by phytoplankton [mol/kg]
        talpha , &	    ! Michaeles Menton-like temperature [Celsius]
        rp0    , &	    ! maximum growth rates of phytoplankton [1/d]
        p0     , &	    ! background concentration for initial phytoplankton growth [mol/kg]
	np     , &	    ! number of P atoms in uptake
	nn     , &	    ! number of N atoms in uptake
	nc     , &	    ! number of C atoms in uptake
	lpd    , &	    ! loss rate of phytoplankton to detritus [1/d]
	lpr    , &	    ! loss rate of phytoplankton by respiration [1/d]
	sinkp		    ! phytoplankton sinking velocity [m/d]
  character(len=32)        :: name_phyt(maxphyt)
                                  
  real, dimension(maxzoo)  :: &
        t_opt_zoo  , &      ! temperature optimum of grazing [Celsius]
        t_max_zoo  , &      ! maximal grazing temperature [Celsius]
        beta_zoo   , &      ! parameter for temperature dependence of grazing [dimensionless]
        oxy_sub_zoo , &     ! oxygen level below which reduced respiration starts [mol/kg]
        oxy_min_zoo, &      ! oxygen level below which no respiration takes place [mol/kg]
        resp_red_zoo, &     ! reduction factor for respiration under suboxic conditions [dimensionless]
	sigma_b, &	    ! loss rate of zooplankton to detritus [mol/kg/d]
	nue    , &	    ! loss rate of zooplankton to nh4 by respiration [mol/kg/d]
        food_to_nh4   , &   ! fraction of eaten food that is directly lost to respiration [dimensionless]
        food_to_det   , &   ! fraction of eaten food that is directly lost to detritus [dimensionless]
        food_to_nh4_2 , &   ! fraction of food eaten potentially at optimal temperature 
	                    ! directly lost to respiration [dimensionless]
        food_to_det_2 , &   ! fraction of food eaten potentially at optimal temperature 
	                    ! directly lost to detritus [dimensionless]
	iv     , &	    ! Ivlev constant [kg/mol]
	zcl1   , &	    ! closure parameter [kg2/mol2]
	graz   , &	    ! zooplankton grazing rate [1/d]
        z0     , &	    ! background concentration for initial zooplankton growth [mol/kg]
	Imax   , &	    ! zooplankton maximum light intensity [W/m^2]
        alpha_zoo, &	    ! light inhibition shape factor
	sinkz  , &	    ! zooplankton sink velocity (<0 for sinking) [m/d]
	risez  , &	    ! zooplankton rise velocity (>0 for rising) [m/d]
	vdiff_max  , &      ! zooplankton maximum enhanced diffusion[m^2/s]
        dark_rise,&	    ! whether zooplankton rises independent off a food gradient in the dark [dimensionless]
	wfood  , &	    ! weight for food gradients in zooplankton rise velocity 
        o2min  , &	    ! minimum oxygen concentration where sinking stops [mol/kg]
        h2smax , &	    ! maximum h2s concentration where sinking stops [mol/kg]
        wtemp  , &	    ! weight number for temperature sensitivity
        wo2    , &	    ! weight number for o2 sensitivity
        wh2s   , &	    ! weight number for h2s sensitivity
	np_zoo , &	    ! number of P atoms in the zooplankton, Redfield ratio
        nn_zoo , &	    ! number of N atoms in the zooplankton, Redfield ratio
        nc_zoo  	    ! number of C atoms in the zooplankton, Redfield ratio
  real, dimension(maxzoo,maxphyt) :: pref_phy    ! food preferences of zooplankton for phytoplankton
  real, dimension(maxzoo,maxzoo)  :: pref_zoo    ! food preferences of zooplankton for zooplankton
  real, dimension(maxzoo,maxdet)  :: pref_det    ! food preferences of zooplankton for detritus
  integer, dimension(maxzoo)      :: graz_pref   ! flag to select grazing preferences
  
  character(len=32)         :: name_zoo(maxzoo)
  logical, dimension(maxzoo):: &
        vertical_migration, & ! enables zooplankton migration				  
        blanchard_temperature ! .false.: old ERGOM temperature dependence, .true.: Blanchard 1996 formula

  real, dimension(maxdet)   :: &
        dn     , &          ! recycling rate
	q10_rec             ! q10 parameter for recycling [1/Celsius]
  character(len=32)         :: &
        name_det(maxdet), &
        name_redfield_sed   = 'sed', &
        name_iron_phosphate = 'ips'

  real, dimension(maxspm) :: &
        wsink0_spm         , &   ! sinking velocity (<0 for sinking) [m/d]
        wsed_spm                 ! sedimentation rate [m/d]
  real, dimension(maxsed) :: &
        erosion_rate_sed   , &   ! erosion rate [1/d]
        bioerosion_rate_sed, &   ! erosion rate by benthic animals [1/d]
        molar_volume_sed   , &   ! volume of this tracer when deposited in the sediment [m3/mol]
        critical_stress_sed      ! critical shear stress when erosion starts [N/m2]
  
  character(len=32) ::           &
        name_spm(maxspm)    , & ! name of this type of spm (suspended particulate matter)
        sediment_to(maxspm) , & ! name of the sed(:) tracer to which sedimentation takes place
        longname_spm(maxspm)	! long name for output

  character(len=32) ::           &
        name_sed(maxsed)     , & ! name of this type of sedimented matter
        suspend_to(maxsed)   , & ! name of spm(:) tracer to which the resuspension takes place
        longname_sed(maxsed)	 ! long name for output
                                  
  real, dimension(max_sediment_layers) :: &
        sed_layer_height         ! (maximum) height of vertical layers [m]. 
	                         ! < 0: the layer may become infinitely thick.
  real  :: &
        nf        = .1       , &  ! nitrification rate [1/d]
        q10_nit   = .11      , &  ! q10 parameter for nitrification [1/Celsius]
        alpha_nit = 3.75e-6  , &  ! half-saturation constant for nitrification [mol/kg]
        q10_h2s   = 0.0693   , &  ! q10 parameter for chemolithotrophs (h2s oxidation) [1/Celsius]
        k_h2s_o2  = 8.e5     , &  ! reaction constant h2s oxidation with o2  [kg/mol/d]
        k_h2s_no3 = 8.e5     , &  ! reaction constant h2s oxidation with no3 [kg/mol/d]
        k_sul_o2  = 2.e4     , &  ! reaction constant sulfur oxidation with o2 [kg/mol/d]
        k_sul_no3 = 2.e4     , &  ! reaction constant sulfur oxidation with no3 [kg/mol/d] 
        ldn_N	  = 5.3      , &  ! stochiometric ratio no3/det for detritus recycling (denitrification) [1]
        ldn_O	  = 6.625    , &  ! stochiometric ratio o2/det for detritus recycling (oxic) [1]
        ldn_S	  = 3.3125   , &  ! stochiometric ratio h2s/det for detritus recycling (sulfate reduction) [1]
	ldn_A	  = 13.25    , &  ! stochiometric ratio no3/det = nh4/det for detritus recycling (anammox) [1]
	alp_o2    =  5.e5    , &  ! smooth oxygen swithch function for detritus recycling [kg/mol]
	alp_no3   =  2.2e6   , &  ! smooth oxygen swithch function for detritus recycling [kg/mol]
	alp_h2s   =  5.e3    , &  ! smooth oxygen swithch function for detritus recycling [kg/mol]
	alp_nh4   =  2.2e6   , &  ! smooth oxygen swithch function for detritus recycling [kg/mol]
        k_an0	  =  .02     , &  ! maximum anammox rate [1/d]
	k_DN	  =  1.      , &  ! 
 	k_DS	  =  1.      , &  !
	den_rate  =  0.5     , &  ! proportion of denitrification at the sediment redoxcline		     
	dn_sed    =  0.003   , &  ! recycling rate of detritus in the sediment [1/d]		     
        frac_dn_anoxic = 0.3 , &  ! fraction of recycling rate in shallow sediments for anoxic bottom water [dimensionless]
        thio_bact_min  = 1.0 , &  ! minimum nitrogen content of active sediment for thiomargarita [mol/m2]
	q10_rec_sed = 0.0693 , &  ! q10 paramter for recycling of detritus in the sediment
	np_sed    =   1.     , &  ! number of P atoms in the sediment, Redfield ratio
        nn_sed    =  16.     , &  ! number of N atoms in the sediment, Redfield ratio
        nc_sed    = 106.     , &  ! number of C atoms in the sediment, Redfield ratio
        po4_lib_rate	 = 0., &  ! fraction of phosphorous retained in the sediment while recycled [1/d]
        po4_retention	 = 0., &  ! fraction of phosphorous retained in the sediment while recycled [dimensionless]
        po4_ret_plus_BB  = 0., &  ! value added to po4_retention north of 60.75N
        			  ! (special treatment for the Bothnian Bay to supress cyanobacterial blooms)
        o2_bioerosion	 = 6.5e-5 ! oxygen thresold to enable bioerosion [mol/kg]
   
  integer :: NUM_PHYTO   = 3 
  integer :: NUM_ZOO     = 1
  integer :: NUM_DET     = 1
  integer :: NUM_SPM     = 1  
  integer :: NUM_SED     = 1  
  integer :: NUM_SEDIMENT_LAYERS   = 2
  integer :: sed_layer_propagation = 1 ! Sediment layer propagation definitions, 
                                            ! SLP_DOWNWARD=1, SLP_FULL_BOX=2, SLP_OLD_ERGOM=3
  integer :: sed_erosion_mode      = 1 ! Sediment erosion mode, INDEPENDENT=1, MAXSTRESS=2, ORGANIC=3
          
  integer :: SLOW = 1
  integer :: FAST = 2
  integer :: n, m
  integer :: vlev_sed = 2                ! number of 2d tracers that may be stored in one diagnostic 3d tracer            
  
  !
  ! Array allocations and flux calculations assume that phyto(1) is the
  ! only phytoplankton group cabable of nitrogen uptake by N2 fixation while phyto(2:NUM_PHYTO) 
  ! are only cabable of nitrogen uptake by NH4 and NO3 uptake
  !
  integer :: DIA     ! = 1    ! diatoms
  integer :: FLA     ! = 2    ! flagellates
  integer :: CYA     ! = 3    ! cyanobacteria
  ! identification numbers for mpp clocks
  integer :: id_source, id_init, id_susp, id_alloc

  data (imin (n),    n=1, maxphyt) /25.,50.,50/                  ! minimum light saturation for phytoplankton growth [W/m^2]
  data (tmin(n),     n=1, maxphyt) /-20.,-20., 20./              ! minimum growth temperature for phytoplankton, 
                                                                 ! here only for diazotrophs [Celsius]
  data (smin(n),     n=1, maxphyt) /-20.,-20.,-20./              ! minimum growth salinity for phytoplankton, 
                                                                 ! here only for diazotrophs [psu]
  data (smax(n),     n=1, maxphyt) /60., 60., 60./               ! maximum growth salinity for phytoplankton, 
                                                                 ! here only for diazotrophs [psu]
  data (alpha(n),    n=1, maxphyt) /1.56e-6, 0.45e-6, 2.25e-6/   ! half-saturation constants of nutrient uptake 
                                                                 ! by phytoplankton [mol/kg]
  data (talpha(n),   n=1, maxphyt) /1.e10, 10., 0./              ! Michaeles Menton-like temperature [Celsius]
  data (rp0  (n),    n=1, maxphyt) /4.0, 1.6, 1.2/               ! growth rates of phytoplankton [1/d]
  data (p0   (n),    n=1, maxphyt) /1.0e-9, 1.0e-9, 1.0e-9/      ! background concentration for initial growth [mol/kg]
  data (np  (n),     n=1, maxphyt) /1.,  1.,  1./                ! number of P atoms in uptake, Redfield ratio
  data (nn  (n),     n=1, maxphyt) /16., 16., 16./               ! number of N atoms in uptake, Redfield ratio
  data (nc  (n),     n=1, maxphyt) /106. , 106., 106./           ! number of C atoms in uptake, Redfield ratio
  data (lpd  (n),    n=1, maxphyt) /0.02, 0.02, 0.02/            ! loss rate of phytoplankton to detritus [mol/kg/d]
  data (lpr  (n),    n=1, maxphyt) /0.01, 0.01, 0.01/            ! loss rate of phytoplankton by respiration [mol/kg/d]
  data (sinkp(n),    n=1, maxphyt) /1., 0., 0./                  ! phytoplankton sink velocity, here only for diatoms [m/d]
  data (name_phyt(n), n=1, maxphyt)/'dia','fla','cya'/           ! 
 

  data (nue    (n),  n=1, maxzoo)  /0.01, 0.01, 0.01, 0.01/      ! loss rate of zooplankton to nh4 by respiration [mol/kg/d]
  data (food_to_nh4  (n), n=1, maxzoo)  /0.0, 0.0, 0.0, 0.0/     ! fraction of eaten food that is directly lost to 
                                                                 ! respiration [dimensionless]
  data (food_to_det  (n), n=1, maxzoo)  /0.0, 0.0, 0.0, 0.0/     ! fraction of eaten food that is directly lost to 
                                                                 ! detritus [dimensionless]
  data (food_to_nh4_2(n), n=1, maxzoo)  /0.18, 0.18, 0.18, 0.18/ ! fraction of food eaten potentially at optimal 
                                                                 ! temperature directly lost to respiration [dimensionless]
  data (food_to_det_2(n), n=1, maxzoo)  /0.18, 0.18, 0.18, 0.18/ ! fraction of food eaten potentially at optimal
                                                                 ! temperature directly lost to detritus [dimensionless]
  data (graz   (n),  n=1, maxzoo)  /0.5, 0.5, 0.5, 0.5 /         ! zooplankton grazing rate [1/d]
  data (z0     (n),  n=1, maxzoo)  /1.0e-9, 1.0e-9, 1.0e-9, 1.0e-9/
                                                                 ! background concentration for initial growth [mol/kg]
  data (sigma_b(n),  n=1, maxzoo)  /0.03, 0.03, 0.03, 0.03/      ! loss rate of zooplankton to detritus [mol/kg/d]
  data (t_opt_zoo(n),   n=1, maxzoo)  /15.0, 10.0, 10.0, 10.0/	 ! temperature optimum of zooplankton for grazing [Celsius]
  data (t_max_zoo(n),   n=1, maxzoo)  /25.0, 15.0, 15.0, 15.0/	 ! maximal grazing temperature [Celsius]
  data (beta_zoo(n),    n=1, maxzoo)  /1.7, 1.7, 1.7, 1.7/       ! parameter for temperature dependence of grazing 
                                                                 ! [dimensionless]
  data (oxy_sub_zoo(n), n=1, maxzoo) /60.e-6, 60.e-6, 60.e-6, 60.e-6/	
                                                                 ! threshold oxygen level for reduced respiration [mol/kg]
  data (oxy_min_zoo(n), n=1, maxzoo)  /5.e-6, 5.e-6, 5.e-6, 5.e-6/ 
                                                                 ! threshold level for no respiration [mol/kg]
  data (resp_red_zoo(n),  n=1, maxzoo)  /1.0, 1.0, 1.0, 1.0 /    ! reduction factor for respiration 
                                                                 ! under suboxic conditions [dimensionless]
  data (iv     (n),  n=1, maxzoo)  /1.e6, 1.e6, 1.e6, 1.e6/      ! Ivlev paramter for zooplankton grazing [kg/mol]
  data (zcl1   (n),  n=1, maxzoo)  /6.67e6, 6.67e6, 6.67e6, 6.67e6/
                                                                 ! closure parameter for zooplankton [kg2/mol2]
  data (Imax   (n),  n=1, maxzoo)  /0.1, 0.1, 0.1, 0.1/          ! maximum light 10µE/m^2 / 4.6
  data (o2min  (n),  n=1, maxzoo)  /5.e-6, 5.e-6, 5.e-6, 5.e-6/  ! minimum oxygen concentration where sinking stops [mol/kg]
  data (h2smax (n),  n=1, maxzoo)  /1.e-6, 1.e-6, 1.e-6, 1.e-6/  ! maximum h2s concentration where sinking stops [mol/kg]
  data (wtemp  (n),  n=1, maxzoo)  /1.,	1., 1., 1. /		 ! weight number for temperature sensitivity
  data (wo2    (n),  n=1, maxzoo)  /1.,	1., 1., 1. /		 ! weight number for o2 sensitivity
  data (wh2s   (n),  n=1, maxzoo)  /1.,	1., 1., 1. /		 ! weight number for h2s sensitivity
  data (np_zoo (n),  n=1, maxzoo)  /1., 1., 1., 1. /		 ! number of P atoms in the zooplankton 
  data (nn_zoo (n),  n=1, maxzoo)  /16., 16., 16., 16. /         ! number of N atoms in the zooplankton 
  data (nc_zoo (n),  n=1, maxzoo)  /106., 106., 106., 106./      ! number of C atoms in the zooplankton 
  data (alpha_zoo(n),  n=1, maxzoo)  /0., 0.012, 0.012, 0.012/     ! light inhibitation parameter
  data (sinkz  (n),  n=1, maxzoo)  /-2000., -2000., -2000., -2000./! zooplankton sink velocity (<0 for sinking) [m/d]
  data (risez  (n),  n=1, maxzoo)  /4000., 4000., 4000., 4000./    ! zooplankton rise velocity  (>0 for rising) [m/d]
  data (vdiff_max  (n),  n=1, maxzoo)  /1., 1., 1., 1./            ! zooplankton maximum enhanced diffusion[m^2/s]
  data (dark_rise(n),n=1, maxzoo)  /1., 0., 0., 0./                ! zooplankton rise velocity  (>0 for rising) [m/d]
  data (wfood  (n),  n=1, maxzoo)  /1., 1., 1., 1./                ! weight for food gradients in zooplankton rise velocity 
  data (vertical_migration(n), n=1, maxzoo)  /.true.,.true.,.true.,.true./       ! logical for vertical migration
  data (blanchard_temperature(n), n=1, maxzoo)  /.true.,.true.,.true.,.true./    ! .false.: old ERGOM temperature dependence, 
                                                                   ! .true.: use Blanchard 1996 formula
  data ((pref_phy(n, m), n=1, maxzoo), m=1, maxphyt)&
        /0.39, 0.29, 0.  ,0.29, &
         0.39, 0.29, 0.  ,0.29, &
         0.02, 0.02, 0.  ,0.02  &
	/       ! 
  data ((pref_zoo(n, m), n=1, maxzoo), m=1, maxzoo)&
        /0. , 0.2, 0.5, 0., &
         0. , 0. , 0.2, 0., &
         0. , 0. , 0.1, 0., &
         0. , 0. ,  0., 0.  &
	/       ! 
  data ((pref_det(n, m), n=1, maxzoo), m=1, maxdet)&
        /0.2, 0.2, 0.2, 0.0, &
         0. , 0.0, 0.0, 0.0  &
	/       ! 
  data (name_zoo(n), n=1, maxzoo)  /'cop','kr1','kr2','sal'/ ! 
  data (graz_pref (n), n=1, maxzoo) /ERGOM_PREFS,GENUS_IVLEV,GENUS_IVLEV,0/ ! flag to select grazing preferences
  
  data (dn     (n),  n=1, maxdet)  /0.003,  0.003 /	  ! recycling rate of detritus [1/d]
  data (q10_rec(n),  n=1, maxdet)  /0.0693, 0.0693/	  ! q10 paramter for recycling of detritus [1/Celsius]
  data (name_det(n), n=1, maxdet)  /'det','fast'/         ! 

  data(wsink0_spm (n) , n=1, maxspm) /-3.0, -1.0/	  ! sinking velocity (<0 for sinking) [m/d]
  data(wsed_spm   (n) , n=1, maxspm) / 2.5,  0.5/	  ! sedimentation rate [m/d]
  data(name_spm(n)    , n=1, maxspm) /'det','ipw'/	  ! name of spm tracer
  data(sediment_to(n) , n=1, maxspm) /'none','none'/	  ! name of sed tracer to which sedimentation takes place
  
  data(name_sed(n)    , n=1, maxsed) /'sed','ips'/	  ! name of sed tracer
  data(suspend_to(n)  , n=1, maxsed) /'none','none'/	  ! name of spm tracer to which resuspension takes place
  data(longname_sed(n), n=1, maxsed) /'detritus','iron phosphate'/ 
                                                          ! long name for output
  data(erosion_rate_sed(n), n=1, maxsed) /6.0, 6.0/       ! erosion rate [1/d]
  data(bioerosion_rate_sed(n), n=1, maxsed) /0.0, 0.0/    ! erosion rate by benthic animals [1/d]
  data(molar_volume_sed(n),n=1, maxsed) /0.01111, 0.0/    ! molar tracer volume [m3/mol]
  data(critical_stress_sed(n),n=1, maxsed) /0.196, 0.196/ ! critical shear stress when erosion starts [N/m2]
  
  data(sed_layer_height(n),n=1,max_sediment_layers) /0.02222, -1.0/  
                                                          ! (maximum) height of vertical layers [m]. 
                                                          ! <0 mean the layer may become infinitely thick.
 
  namelist /ergom_nml/  &
! dimensions of tracer arrays
   NUM_PHYTO    , &
   NUM_ZOO      , &
   NUM_DET      , &
   NUM_SPM      , &
   NUM_SED      , &
   NUM_SEDIMENT_LAYERS , &
! phytoplankton parameters
   name_phyt    , & ! name of phytoplankton variable
   imin         , & ! minimum light saturation for phytoplankton growth [W/m²]
   tmin         , & ! minimum growth temperature for phytoplankton
   smin         , & ! minimum growth salinity for phytoplankton
   smax         , & ! maximum growth salinity for phytoplankton
   alpha        , & ! half-saturation constants of nutrient uptake by phytoplankton [mol/kg]
   talpha       , & ! Michaeles Menton-like temperature [Celsius]
   rp0          , & ! growth rates of phytoplankton [1/d]
   p0           , & ! background concentration for initial phytoplankton growth [mol/kg]
   np           , & ! number of P atoms in uptake, Redfield ratio
   nn           , & ! number of N atoms in uptake, Redfield ratio
   nc           , & ! number of C atoms in uptake, Redfield ratio
   lpd          , & ! phytoplankton loss rates to detritus [1/d]
   lpr          , & ! phytoplankton respiration rates to nh4 [1/d]
   sinkp        , & ! phytoplankton sinking velocities, i.e. here only for diatoms [m/d]
! zooplankton parameters
   name_zoo     , & ! name of zooplankton variable
   sinkz        , & ! zoplankton sinking velocities, [m/d]
   risez        , & ! zoplankton rise velocities, [m/d]
   vdiff_max    , & ! zoplankton maximum enhanced diffusion [m^2/s]
   wfood        , & ! weight for food gradients in zooplankton rise velocity 
   o2min        , & ! minimum oxygen concentration where sinking stops [mol/kg]
   h2smax       , & ! maximum h2s concentration where sinking stops [mol/kg]
   wtemp        , & ! weight number for temperature sensitivity
   wo2          , & ! weight number for o2 sensitivity
   wh2s         , & ! weight number for h2s sensitivity
   t_opt_zoo    , & ! temperature optimum of zooplankton for grazing [Celsius]
   t_max_zoo    , & ! temperature maximum of zooplankton for grazing [Celsius]
   beta_zoo     , & ! parameter for temperature dependence of grazing [dimensionless]
   oxy_sub_zoo ,  & ! oxygen level below which reduced respiration starts [mol/kg]
   oxy_min_zoo,   & ! oxygen level below which no respiration takes place [mol/kg]
   resp_red_zoo,  & ! reduction factor for respiration under suboxic conditions [dimensionless]
   sigma_b      , & ! zooplankton loss rate to detritus [1/d]
   nue          , & ! zooplankton respiration rate to nh4 [1/d]
   food_to_nh4  , & ! fraction of eaten food that is directly lost to respiration [dimensionless]
   food_to_det  , & ! fraction of eaten food that is directly lost to detritus [dimensionless]
   food_to_nh4_2, & ! fraction of food eaten potentially, directly lost to respiration [dimensionless]
   food_to_det_2, & ! fraction of food eaten potentially, directly lost to detritus [dimensionless]
   iv	        , & ! Ivlev paramter for zooplankton grazing [kg/mol]
   zcl1         , & ! closure parameter for zooplankton [kg/mol]
   graz         , & ! zooplankton maximum grazing rate [1/d] 
   z0           , & ! background concentration for initial zooplankton growth [mol/kg]
   Imax         , & ! maximum light for zooplankton activity [W/m²]
!   Iopt         , & ! optimum light [W/m²]
   vertical_migration, & ! zooplankton vertical migration
   blanchard_temperature, & ! .false.: old ERGOM temperature dependence, .true.: Blanchard 1996 formula
   pref_phy     , & ! preferences of zooplankton for phytoplankton(dimensionless)   
   pref_zoo     , & ! preferences of zooplankton for zooplankton (dimensionless)   
   pref_det     , & ! preferences of zooplankton for detritus (dimensionless)   
   graz_pref    , &  ! flag to select grazing preferences
! detritus parameters
   name_det     , & ! name of detritus, must be equal to a suspended particulate matter (spm) variable name
   dn	        , & ! recycling rate [1/d]
   q10_rec      , & ! q10 paramter for recycling of detritus [1/Celsius]
! generic_ERGOM_type parameters
   nf	        , & ! nitrification rate [1/d]
   q10_nit      , & ! q10 parameter for nitrification [1/Celsius]
   alpha_nit    , & ! half-saturation constant for nitrification [mol/kg]
   q10_h2s      , & ! q10 parameter for chemolithotrophs (h2s oxidation) [1/Celsius]
   k_h2s_o2     , & ! reaction constant h2s oxidation with o2  [kg/mol/d]
   k_h2s_no3    , & ! reaction constant h2s oxidation with no3 [kg/mol/d]
   k_sul_o2     , & ! reaction constant sulfur oxidation with o2 [kg/mol/d]
   k_sul_no3    , & ! reaction constant sulfur oxidation with no3 [kg/mol/d]
   k_an0        , & ! maximum anammox rate [1/d]
   k_DN         , & ! 
   k_DS         , & !
! sediment parameters
   name_redfield_sed    , &
   name_iron_phosphate  , &
   dn_sed       , & ! recycling rate of detritus in the sediment [1/d]
   frac_dn_anoxic, &! fraction of recycling rate in shallow sediments for anoxic bottom water [dimensionless]
   thio_bact_min, & ! minimum nitrogen content of active sediment for thiomargarita [mol/m2]
   q10_rec_sed  , & ! q10 paramter for recycling of detritus in the sediment [1/Celsius]
   den_rate     , & ! porportion of denitrification at the sediment redoxcline
   np_sed       , & ! number of P atoms in the sediment, Redfield ratio
   nn_sed       , & ! number of N atoms in the sediment, Redfield ratio
   nc_sed       , & ! number of C atoms in the sediment, Redfield ratio
   po4_lib_rate , & ! liberation rate of iron phosphate in the sediment [1/d]
   po4_retention, & ! fraction of phosphorous retained in the sediment while recycled [dimensionless]
   po4_ret_plus_BB, &! value added to po4_retention north of 60.75N 
                     !(special treatment for the Bothnian Bay to supress cyanobacterial blooms)
   o2_bioerosion  , &! oxygen thresold to enable bioerosion  [mol/kg]
! suspended particulate matter parameters
   name_spm       , & ! name of spm tracer
   longname_spm   , & ! long name for output
   wsink0_spm     , & ! sinking velocity (<0 for sinking) [m/d]
   wsed_spm       , & ! sedimentation rate [m/s]
   sediment_to    , & ! name of 2d tracer to which sedimentation takes place   
! settled matter parameters
   name_sed            , & ! name of sed tracer
   longname_sed        , & ! long name for output
   molar_volume_sed    , & ! specific volume of this sediment type [m3/mol]
   critical_stress_sed , & ! critical shear stress when erosion starts [N/m2]
   erosion_rate_sed    , & ! erosion rate [1/d]
   bioerosion_rate_sed , & ! erosion rate by benthic animals [1/d]
   suspend_to          , & ! name of 3d tracer to which resuspension  takes place
! sediment settings parameters
   sed_layer_height       , & ! maximum height of vertical layers [m]. <0 mean the layer may become infinitely thick.
   sed_layer_propagation  , & ! Sediment layer propagation settings, SLP_DOWNWARD=1, SLP_FULL_BOX=2, SLP_OLD_ERGOM=3
   sed_erosion_mode       , & ! Sediment erosion mode, INDEPENDENT=1, MAXSTRESS=2, ORGANIC=3
! other parameters
   vlev_sed                        ! number of 2d tracers that may be stored in one diagnostic 3d tracer
                                   ! set this parameter <= nk

contains

  subroutine generic_ERGOM_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_register'
    character(len=fm_string_len) :: errorstring
    integer :: ioun, io_status, ierr, i, j
    logical :: found
    integer :: stdoutunit,stdlogunit

    stdoutunit=stdout();stdlogunit=stdlog()

    call write_version_number(version, tagname)

    ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ergom_nml, iostat=io_status)
#else
    ioun = open_namelist_file()
    read  (ioun, ergom_nml,iostat=io_status)
    call close_file (ioun)
#endif
    ierr = check_nml_error(io_status,'ergom_nml')

    write (stdoutunit,'(/)')
    write (stdoutunit, ergom_nml)
    write (stdlogunit, ergom_nml)

    allocate(phyto  (NUM_PHYTO)  )
    allocate(zoo    (NUM_ZOO)    )
    allocate(det    (NUM_DET)    )
    allocate(spm    (NUM_SPM)    )
    allocate(sed    (NUM_SED)    )
    allocate(sed_defs%layer_height(NUM_SEDIMENT_LAYERS))

    ! now, set the indexes for spm variables referred to
    ! A ) Set the values in sed
    do i=1,NUM_SED
      ! search the index of the spm variables to which resuspension takes place
      found = .false.
      sed(i)%index_suspend_to = -1
      do j=1,NUM_SPM
	if (trim(adjustl(suspend_to(i))) .eq. trim(adjustl(name_spm(j)))) then
	  found = .true.
	  sed(i)%index_suspend_to = j
	endif
      enddo
      if (.not. found) then
	if ((trim(adjustl(suspend_to(i))) .eq. 'none') .and. (NUM_SPM .eq. NUM_SED)) then
	  sed(i)%index_suspend_to = i
	else
          write(errorstring, '(a)') &
		'Error: settled matter tracer '// &
		 trim(adjustl(name_sed(i)))	// &
		 ' shall be resuspended to tracer '	// & 
		 trim(adjustl(suspend_to(i))) // &
		 ', but that does not exist as an spm tracer.'
	  call  mpp_error(FATAL, errorstring)
	endif
      endif
    enddo
    ! B) Set the values in det
    do i=1,NUM_DET
      found = .false.
      do j=1,NUM_SPM
        if (trim(adjustl(name_det(i))) .eq. trim(adjustl(name_spm(j)))) then
          found = .true.
          det(i)%index_spm = j
        endif
      enddo
      if (.not. found) then
	write(errorstring, '(a)') &
                   'Error: detritus tracer '// &
                    trim(adjustl(name_det(i)))     // &
                   ' does not exist as a suspended particulate matter (spm) tracer.'
        call  mpp_error(FATAL, errorstring)
      endif
    enddo
    
    ! now, set the indexes for sed variables referred to
    ! A) Set the values in spm
    do i=1,NUM_SPM
      ! search the index of the sed variables to which sedimentation takes place
      found = .false.
      spm(i)%index_sediment_to = -1
      do j=1,NUM_SED
	if (trim(adjustl(sediment_to(i))) .eq. trim(adjustl(name_sed(j)))) then
	  found = .true.
	  spm(i)%index_sediment_to = j
	endif
      enddo
      if (.not. found) then
	if ((trim(adjustl(sediment_to(i))) .eq. 'none') .and. (NUM_SPM .eq. NUM_SED)) then
	  spm(i)%index_sediment_to = i
	else
          write(errorstring, '(a)') &
               'Error: suspended particulate matter (spm) tracer '// &
		trim(adjustl(name_spm(i)))     // &
		' shall be sedimented to tracer '     // & 
		trim(adjustl(sediment_to(i))) // &
		', but that does not exist as a settled matter (sed) tracer.'
	  call  mpp_error(FATAL, errorstring)
	endif
      endif
    enddo
    ! C) Set the values in biosed
    found = .false.
    do i=1,NUM_SED
      if (trim(adjustl(name_redfield_sed)) .eq. trim(adjustl(name_sed(i)))) then
        found=.true.
        biosed%index_sed=i
      endif
    enddo
    if (.not. found) then
      write(errorstring, '(a)') &
          'Error: redfield-ratio sediment tracer '// &
          trim(adjustl(name_redfield_sed))     // &
          ' does not exist as a settled matter (sed) tracer.'
      call  mpp_error(FATAL, errorstring)
    endif
    found = .false.
    do i=1,NUM_SED
      if (trim(adjustl(name_iron_phosphate)) .eq. trim(adjustl(name_sed(i)))) then
        found=.true.
        biosed%index_ips=i
      endif
    enddo
    if ((.not. found) .and. (po4_retention .gt. 0.0)) then
      write(errorstring, '(a)') &
          'Error: iron phosphate sediment tracer '// &
          trim(adjustl(name_iron_phosphate))     // &
          ' does not exist as a settled matter (sed) tracer, '// &
	  'but po4_retention is not zero.'
      call  mpp_error(FATAL, errorstring)
    endif

    DIA      = 1    ! diatoms
    FLA      = 2    ! flagellates
    CYA      = 3    ! cyanobacteria

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)
    
  end subroutine generic_ERGOM_register

  ! <SUBROUTINE NAME="generic_ERGOM_init">
  !  <OVERVIEW>
  !   Initialize the generic ERGOM module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the CFC Tracers to the list of generic Tracers passed to it 
  !       via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_ERGOM_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_init'
    integer :: stdoutunit,stdlogunit

    stdoutunit=stdout();stdlogunit=stdlog()

    id_init   =   mpp_clock_id('(ERGOM init) '           ,grain=CLOCK_ROUTINE)
    id_alloc  =   mpp_clock_id('(ERGOM allocate) '       ,grain=CLOCK_ROUTINE)
    id_source = mpp_clock_id('(ERGOM source terms) '     ,grain=CLOCK_ROUTINE)
    id_susp   =   mpp_clock_id('(ERGOM resuspension) '   ,grain=CLOCK_ROUTINE)
    call mpp_clock_begin(id_init)

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    call user_allocate_arrays 
    
    ! now print a summary of all parameters    
    write (stdoutunit,'(/)')
    write (stdoutunit,*) 'Summary of the ERGOM model setup'
    write (stdoutunit,'(a,I2)') 'Number of phytoplankton types : ', NUM_PHYTO
    write (stdoutunit,'(a,I2)') 'Number of zoooplankton types  : ', NUM_ZOO
    write (stdoutunit,'(a,I2)') 'Number of detritus types      : ', NUM_DET
    write (stdoutunit,'(a,I2)') 'Number of SPM types           : ', NUM_SPM
    write (stdoutunit,'(a,I2)') 'Number of settled matter types: ', NUM_SED
    do n=1, NUM_PHYTO
       write (stdoutunit,'(a)')          phyto(n)%name//':'
       write (stdoutunit,'(a)')         '  Parameters: '
       write (stdoutunit,'((a), e13.6)')'    P/N ratio, pnr			        : ', phyto(n)%pnr
       write (stdoutunit,'((a), e13.6)')'    C/N ratio, cnr			        : ', phyto(n)%cnr
       write (stdoutunit,'((a), e13.6)')'    seed concentration, p0, [mol/kg]	        : ', phyto(n)%p0
       write (stdoutunit,'(a)')         '  Parameters for growth: '
       write (stdoutunit,'((a), e13.6)')'    minimum light		  imin [W/m2]   : ', phyto(n)%imin  
       write (stdoutunit,'((a), e13.6)')'    minimum temperature  	  tmin [C]      : ', phyto(n)%tmin  
       write (stdoutunit,'((a), e13.6)')'    minimum phytoplankton salinity smin [g/kg]   : ', phyto(n)%smin  
       write (stdoutunit,'((a), e13.6)')'    maximum phytoplankton salinity smax [g/kg]   : ', phyto(n)%smax  
       write (stdoutunit,'((a), e13.6)')'    DIN half-sat constant, alpha [mol/kg]        : ', phyto(n)%alpha 
       write (stdoutunit,'((a), e13.6)')'    temperature dep. uptake, talpha  [Celsius]   : ', phyto(n)%talpha 
       write (stdoutunit,'((a), e13.6)')'    maximum uptake rate  	  rp0  [1/s]    : ', phyto(n)%rp0   
       write (stdoutunit,'(a)')         '  Parameters for losses: '
       write (stdoutunit,'((a), e13.6)')'    loss to detritus, lpd  [1/s] 	        : ', phyto(n)%lpd
       write (stdoutunit,'((a), e13.6)')'    loss by respiration, lpr [1/s]	        : ', phyto(n)%lpr
       write (stdoutunit,'((a), e13.6)')'    sinking velocity, wsink0 [m/s]	        : ', phyto(n)%wsink0
       write (stdoutunit,'(/)')
    enddo
    do n=1, NUM_ZOO
       write (stdoutunit,'(a)')             zoo(n)%name//':'
       write (stdoutunit,'(a,(5f7.4,2x))')'  Grazing preferences phytoplankton            : ', &
                                         (zoo(n)%pref_phy(m),m=1,NUM_PHYTO)
       write (stdoutunit,'(a,(5f7.4,2x))')'  Grazing preferences zooplankton              : ', &
                                         (zoo(n)%pref_zoo(m),m=1,NUM_ZOO)
       write (stdoutunit,'(a,(5f7.4,2x))')'  Grazing preferences detritus                 : ', &
                                         (zoo(n)%pref_det(m),m=1,NUM_DET)
       write (stdoutunit,'(a)')           '  Parameters: '
       write (stdoutunit,'((a), e13.6)')  '    P/N ratio, pnr				: ', zoo(n)%pnr
       write (stdoutunit,'((a), e13.6)')  '    C/N ratio, cnr				: ', zoo(n)%cnr
       write (stdoutunit,'((a), e13.6)')  '    seed concentration, z0, [mol/kg]		: ', zoo(n)%z0
       write (stdoutunit,'(a)')           '  Parameters for grazing: '
       write (stdoutunit,'((a), I2)')     '    grazing preference method                  : ', zoo(n)%graz_pref
       write (stdoutunit,'((a), e13.6)')  '    Ivlev constant, iv,         [kg/mol]       : ', zoo(n)%iv
       write (stdoutunit,'((a), e13.6)')  '    maximum grazing rate, graz, [1/s]          : ', zoo(n)%graz
       write (stdoutunit,'((a), e13.6)')  '    maximal grazing temperature, t_max [C]     : ', zoo(n)%t_max
       write (stdoutunit,'((a), e13.6)')  '    optimal grazing temperature, t_opt [C]     : ', zoo(n)%t_opt
       write (stdoutunit,'((a), e13.6)')  '    parameter for temperature dependence, beta : ', zoo(n)%beta
       write (stdoutunit,'((a), e13.6)')  '    fraction lost to respiration, food_to_nh4  : ', zoo(n)%food_to_nh4
       write (stdoutunit,'((a), e13.6)')  '    fraction lost to detritus                  : ', zoo(n)%food_to_det
       write (stdoutunit,'(a)')           '  Parameters for migration: '
       if (zoo(n)%vertical_migration) then
         write (stdoutunit,'(a)')         '    migrating '
         write (stdoutunit,'((a), e13.6)')'    maximum light intensity, Imax [W/m^2]	: ', zoo(n)%imax
         write (stdoutunit,'((a), e13.6)')'    minimum oxygen conc., o2min [mol/kg]	: ', zoo(n)%o2min
         write (stdoutunit,'((a), e13.6)')'    maximal temp for migration, t_max [C]	: ', zoo(n)%t_max
         write (stdoutunit,'((a), e13.6)')'    maximum rise velocity, wrise0 [m/s]	: ', zoo(n)%wrise0
         write (stdoutunit,'((a), e13.6)')'    maximum sink velocity, wsink0 [m/s]	: ', zoo(n)%wsink0
         write (stdoutunit,'((a), e13.6)')'    maximum enhanced diff., vdiff_max [m2/s]	: ', zoo(n)%vdiff_max
       else
         write (stdoutunit,'(a)')         '    not migrating '
       endif
       write (stdoutunit,'(a)')           '  Closure term: '
       write (stdoutunit,'((a), e13.6)')  '    closure parameter, zcl1 [kg/mol]           : ', zoo(n)%zcl1
       write (stdoutunit,'((a), e13.6)')  '    loss rate by respiration, nue [1/s]        : ', zoo(n)%nue
       write (stdoutunit,'((a), e13.6)')  '    loss rate to detritus, sigma_b [1/s]       : ', zoo(n)%sigma_b
       write (stdoutunit,'(a)')           '  Respiration: '
       write (stdoutunit,'((a), e13.6)')  '    reduction factor for respiration, resp_red : ', zoo(n)%resp_red
       write (stdoutunit,'((a), e13.6)')  '    max oxy. f. reduced resp., oxy_sub [mol/kg]: ', zoo(n)%oxy_sub
       write (stdoutunit,'((a), e13.6)')  '    min oxy. f. resp., oxy_min [mol/kg]        : ', zoo(n)%oxy_min

       write (stdoutunit,'(/)')
    enddo
    do n=1, NUM_DET
       write (stdoutunit,'(a)') det(n)%name//':'
       write (stdoutunit,'(/)')
    enddo
    do n=1, NUM_SPM
       write (stdoutunit,'(a)')            trim(spm(n)%name)//':'
       write (stdoutunit,'((a), e13.6)')  '  sinking velocity, wsink0 [m/d]               : ', spm(n)%wsink0  
       write (stdoutunit,'((a), e13.6)')  '  sedimentation rate, wsed [m/d]               : ', spm(n)%wsed    
       write (stdoutunit,'((a), (a))')    '  will sediment to tracer, (sediment_to)       : ', trim(sed(spm(n)%index_sediment_to)%name) 
    enddo
    do n=1, NUM_SED
       write (stdoutunit,'(a)')            trim(sed(n)%name)//':'
       write (stdoutunit,'((a), e13.6)')  '  critical shear stress, critical_stress [N/m2]: ', sed(n)%critical_stress  
       write (stdoutunit,'((a), (a))')    '  will be resuspended to tracer, (suspend_to)  : ', trim(spm(sed(n)%index_suspend_to)%name) 
    enddo
    write (stdoutunit,'(a)')          'Sediment parameters:'
    write (stdoutunit,'((a), e13.6)') '  recycling rate, dn [1/s]                                              : ', &
                                    biosed%dn	       
    write (stdoutunit,'((a), e13.6)') '  frac. rec.-rate in shallow sediments when anoxic, frac_dn_anoxic      : ', &
                                    biosed%frac_dn_anoxic   
    write (stdoutunit,'((a), e13.6)') '  minimum amount of active sed for thiomargarita, thio_bact_min [mol/m2]: ', & 
                                    biosed%thio_bact_min    
    write (stdoutunit,'((a), e13.6)') '  q10 parameter for recycling, q10_rec [1/C]                            : ', & 
                                    biosed%q10_rec	       
    write (stdoutunit,'((a), e13.6)') '  proportion of denit in sediment, den_rate                             : ', & 
                                    biosed%den_rate         
    write (stdoutunit,'((a), e13.6)') '  P/N ratio, pnr                                                        : ', &
                                    biosed%pnr	       
    write (stdoutunit,'((a), e13.6)') '  C/N ratio, cnr                                                        : ', &
                                    biosed%cnr	       
    write (stdoutunit,'((a), e13.6)') '  liberation rate for iron phosphate, po4_lib_rate  [1/s]               : ', & 
                                    biosed%po4_lib_rate     
    write (stdoutunit,'((a), e13.6)') '  fraction of phosphorous retained in the sediment, po4_retention       : ', & 
                                    biosed%po4_retention    
    write (stdoutunit,'((a), e13.6)') '  value added to po4_retention north of 60.75N, po4_ret_plus_BB         : ', & 
                                    biosed%po4_ret_plus_BB  
    write (stdoutunit,'(/)')
    write (stdoutunit,'(a)')          'Ergom parameters:'
    write (stdoutunit,'((a), e13.6)') '  q10 parameter for nitrification, q10_nit [1/C]                        : ', & 
                                    ergom%q10_nit        
    write (stdoutunit,'((a), e13.6)') '  q10 parameter for chemolithotrophs (so4 reduction), q10_h2s   [1/C]   : ', & 
                                    ergom%q10_h2s        
    write (stdoutunit,'((a), e13.6)') '  nitrification rate, nf [1/s]                                          : ', & 
                                    ergom%nf     
    write (stdoutunit,'((a), e13.6)') '  half-saturation constant for nitrification, alpha_nit [mol/kg]        : ', & 
                                    ergom%alpha_nit 
    write (stdoutunit,'((a), e13.6)') '  slope function for detritus recycling, alp_o2   [kg/mol]              : ', & 
                                    ergom%alp_o2 
    write (stdoutunit,'((a), e13.6)') '  slope function for detritus recycling, alp_no3  [kg/mol]              : ', & 
                                    ergom%alp_no3        
    write (stdoutunit,'((a), e13.6)') '  slope function for detritus recycling, alp_h2s  [kg/mol]              : ', & 
                                    ergom%alp_h2s        
    write (stdoutunit,'((a), e13.6)') '  slope function for detritus recycling, alp_nh4  [kg/mol]              : ', & 
                                    ergom%alp_nh4        
    write (stdoutunit,'((a), e13.6)') '  reaction constant h2s oxidation with o2, k_h2s_o2 [kg/mol/s]          : ', & 
                                    ergom%k_h2s_o2  
    write (stdoutunit,'((a), e13.6)') '  reaction constant h2s oxidation with no3, k_h2s_no3 [kg/mol/s]        : ', & 
                                    ergom%k_h2s_no3 
    write (stdoutunit,'((a), e13.6)') '  reaction constant sulfur oxidation with o2, k_sul_o2 [kg/mol/s]       : ', & 
                                    ergom%k_sul_o2  
    write (stdoutunit,'((a), e13.6)') '  reaction constant sulfur oxidation with no3, k_sul_no3 [kg/mol/s]     : ', & 
                                    ergom%k_sul_no3 
    write (stdoutunit,'((a), e13.6)') '  maximum anammox rate, k_an0 [1/s]                                     : ', & 
                                    ergom%k_an0	 
    write (stdoutunit,'(/)')

!! Do not delete, men at work
!!  food_to_nh4_2 ! fraction of food eaten potentially at optimal temperature directly lost to respiration [dimensionless]
!!  food_to_det_2 ! fraction of food eaten potentially at optimal temperature directly lost to detritus [dimensionless]
!!  alpha	  ! light inhibition shape factor
!!  h2smax	  ! maximum h2s concentration where sinking stops [mol/kg]
!!  wo2 	  ! weight number for o2 sensitivity
!!  wh2s	  ! weight number for h2s sensitivity
!!  dark_rise	  ! whether zooplankton rises independent from a food gradient if it is dark [dimensionless]
!!  wfood	  ! weight number for food gradiens 
!!

    call mpp_clock_end(id_init)


  end subroutine generic_ERGOM_init

  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, n, i
    character(len=fm_string_len) :: mystring

    call mpp_clock_begin(id_alloc)

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 
    !Allocate all the private arrays.
    allocate(ergom%f_no3(isd:ied, jsd:jed, 1:nk));         ergom%f_no3=0.0
    allocate(ergom%f_nh4(isd:ied, jsd:jed, 1:nk));         ergom%f_nh4=0.0
    allocate(ergom%f_po4(isd:ied, jsd:jed, 1:nk));         ergom%f_po4=0.0
    allocate(ergom%f_o2 (isd:ied, jsd:jed, 1:nk));         ergom%f_o2=0.0
    allocate(ergom%f_h2s (isd:ied, jsd:jed, 1:nk));        ergom%f_h2s=0.0
    allocate(ergom%f_sul (isd:ied, jsd:jed, 1:nk));        ergom%f_sul=0.0
    allocate(ergom%f_chl (isd:ied, jsd:jed, 1:nk));        ergom%f_chl=0.0
    allocate(ergom%irr_inst(isd:ied, jsd:jed, 1:nk));      ergom%irr_inst=0.0
    allocate(ergom%jno3 (isd:ied, jsd:jed, 1:nk));         ergom%jno3=0.0
    allocate(ergom%jnh4 (isd:ied, jsd:jed, 1:nk));         ergom%jnh4=0.0
    allocate(ergom%jpo4 (isd:ied, jsd:jed, 1:nk));         ergom%jnh4=0.0
    allocate(ergom%jdenit_wc (isd:ied, jsd:jed, 1:nk));    ergom%jdenit_wc=0.0
    allocate(ergom%jnitrif (isd:ied, jsd:jed, 1:nk));      ergom%jnitrif=0.0
    allocate(ergom%jo2  (isd:ied, jsd:jed, 1:nk));         ergom%jo2=0.0
    allocate(ergom%jh2s  (isd:ied, jsd:jed, 1:nk));        ergom%jh2s=0.0
    allocate(ergom%jsul  (isd:ied, jsd:jed, 1:nk));        ergom%jsul=0.0
    allocate(ergom%b_o2(isd:ied, jsd:jed));                ergom%b_o2=0.0
    allocate(ergom%b_no3(isd:ied, jsd:jed));               ergom%b_no3=0.0
    allocate(ergom%b_nh4(isd:ied, jsd:jed));               ergom%b_nh4=0.0
    allocate(ergom%b_po4(isd:ied, jsd:jed));               ergom%b_po4=0.0
    allocate(ergom%b_nitrogen(isd:ied, jsd:jed));          ergom%b_nitrogen=0.0
    allocate(ergom%b_h2s(isd:ied, jsd:jed));               ergom%b_h2s=0.0
    do n = 1, NUM_PHYTO
       allocate(phyto(n)%f_n(isd:ied,jsd:jed,nk));         ; phyto(n)%f_n        = 0.0
       allocate(phyto(n)%jgraz_n(isd:ied,jsd:jed,nk))      ; phyto(n)%jgraz_n    = 0.0
       allocate(phyto(n)%jprod_po4(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_po4  = 0.0
       allocate(phyto(n)%jres_n(isd:ied,jsd:jed,nk))       ; phyto(n)%jres_n     = 0.0
       allocate(phyto(n)%jdet_n(isd:ied,jsd:jed,nk))       ; phyto(n)%jdet_n     = 0.0
    enddo
    allocate(phyto(dia)%move  (isd:ied,jsd:jed,nk))        ; phyto(dia)%move     = 0.0
    allocate(phyto(cya)%move  (isd:ied,jsd:jed,nk))        ; phyto(cya)%move     = 0.0
    do n = 1, 2
       allocate(phyto(n)%jprod_nh4(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_nh4  = 0.0
       allocate(phyto(n)%jprod_no3(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_no3  = 0.0
    enddo
    allocate(phyto(CYA)%jprod_n2(isd:ied,jsd:jed,nk))      ; phyto(CYA)%jprod_n2   = 0.0
    do n = 1, NUM_ZOO
      allocate(zoo(n)%f_n (isd:ied, jsd:jed, 1:nk))	  ; zoo(n)%f_n      = 0.0
      allocate(zoo(n)%jgraz_n(isd:ied,jsd:jed,nk))	  ; zoo(n)%jgraz_n  = 0.0
      allocate(zoo(n)%jgain_n(isd:ied,jsd:jed,nk))	  ; zoo(n)%jgain_n  = 0.0
      allocate(zoo(n)%jres_n(isd:ied,jsd:jed,nk))	  ; zoo(n)%jres_n   = 0.0
      allocate(zoo(n)%jdet_n(isd:ied,jsd:jed,nk))	  ; zoo(n)%jdet_n   = 0.0
!       if(vertical_migration(n)) &
!       allocate(zoo(n)%move  (isd:ied,jsd:jed,nk))         ; zoo(n)%move         = 0.0
    enddo
    allocate(biosed%bioerosion(isd:ied,jsd:jed))            ; biosed%bioerosion       = 0.0
    do n = 1, NUM_SPM
      allocate(spm(n)%move     (isd:ied,jsd:jed,nk))  ; spm(n)%move = 0.0
      allocate(spm(n)%btf      (isd:ied,jsd:jed))     ; spm(n)%btf  = 0.0
      allocate(spm(n)%jsed     (isd:ied,jsd:jed))     ; spm(n)%jsed = 0.0
    enddo
    do n = 1, NUM_SED
      allocate(sed(n)%f_sed    (isd:ied,jsd:jed,NUM_SEDIMENT_LAYERS)); sed(n)%f_sed = 0.0
      allocate(sed(n)%jgain_sed(isd:ied,jsd:jed)); sed(n)%jgain_sed = 0.0
      allocate(sed(n)%jloss_sed(isd:ied,jsd:jed)); sed(n)%jloss_sed = 0.0
      allocate(sed(n)%jres     (isd:ied,jsd:jed)); sed(n)%jres      = 0.0
      allocate(sed(n)%jbiores  (isd:ied,jsd:jed)); sed(n)%jbiores   = 0.0
      do i=1,NUM_SEDIMENT_LAYERS
	write( mystring, '(i4)' )  i
	call user_2d_tracer_assign_array(trim(sed(n)%name)//'_'//trim(adjustl(mystring)), &
    					 sed(n)%f_sed(:,:,i))
      enddo
    enddo
    do n = 1, NUM_DET
       allocate(det(n)%f_n (isd:ied, jsd:jed, 1:nk))       ; det(n)%f_n            = 0.0
       allocate(det(n)%jgraz_n(isd:ied,jsd:jed,nk))        ; det(n)%jgraz_n        = 0.0
       allocate(det(n)%jmort(isd:ied,jsd:jed,nk))          ; det(n)%jmort          = 0.0
    enddo

    call mpp_clock_end(id_alloc)

  end subroutine user_allocate_arrays

  subroutine generic_ERGOM_register_diag(diag_list)
    type(g_diag_type), pointer :: diag_list
    type(vardesc)  :: vardesc_temp
    integer        :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, axes(3)
    type(time_type):: init_time 

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes=axes,init_time=init_time) 


    !   The following > types contain a package of metadata about each tracer,
    ! including, in order, the following elements: name; longname; horizontal
    ! staggering ('h') for collocation with thickness points ; vertical staggering
    ! ('L') for a layer variable ; temporal staggering ('s' for snapshot) ; units ;
    ! and precision in non-restart output files ('f' for 32-bit float or 'd' for
    ! 64-bit doubles). For most tracers, only the name, longname and units should
    ! be changed.  
    vardesc_temp = vardesc("dep_wet_nh4","Wet Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_dep_wet_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_no3","Wet Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_dep_wet_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_po4","Wet Deposition of Phosphate to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_dep_wet_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_nh4","Dry Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_dep_dry_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_no3","Dry Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_dep_dry_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_po4","Dry Deposition of Phosphate to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_dep_dry_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_po4","Phosphate runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_runoff_flux_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_nh4","Ammonia runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_runoff_flux_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_no3","Nitrate runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    ergom%id_runoff_flux_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)


    ! Register the diagnostics for the various phytoplankton 
    !
    ! Register Limitation Diagnostics
    !
    vardesc_temp = vardesc("def_nit","molecular nitrogen",'h','L','s','unknown units','f')
    id_def_nit = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("irr_inst","Instantaneous Light",'h','L','s','W m-2','f')
    ergom%id_irr_inst = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jnh4","NH4 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jno3","NO3 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jno3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jdenit_wc","Water column Denitrification layer integral",'h','L','s',&
                           'mol m-2 s-1','f')
    ergom%id_jdenit_wc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jo2","O2 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jh2s","H2S source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jh2s = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jh2s_o2","H2S with O2 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jh2s_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jh2s_no3","H2S with NO3 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jh2s_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jsul_o2","sulfur with O2 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jsul_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jsul_no3","sulfur with NO3 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jsul_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jsul","sulfur source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jsul = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jpo4","PO4 source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jpo4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jnitrif","Nitrification source layer integral",'h','L','s','mol m-2 s-1','f')
    ergom%id_jnitrif = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jrec_o2", "Nitrogen flux to NH4, recycling o2",'h','L','s','mol m-2 s-1','f')
    ergom%id_jrec_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jrec_no3","Nitrogen flux to NH4, recycling no3",'h','L','s','mol m-2 s-1','f')
    ergom%id_jrec_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jrec_so4","Nitrogen flux to NH4, recycling so4",'h','L','s','mol m-2 s-1','f')
    ergom%id_jrec_so4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jrec_ana","Nitrogen flux to NH4, recycling anamox",'h','L','s','mol m-2 s-1','f')
    ergom%id_jrec_ana = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    if(ergom%id_jh2s_o2 .gt. 0) then
      allocate(ergom%jh2s_o2(isd:ied,jsd:jed,nk)) 
      ergom%jh2s_o2  = 0.0
    endif
    if(ergom%id_jh2s_no3.gt. 0) then 
      allocate(ergom%jh2s_no3(isd:ied,jsd:jed,nk))
      ergom%jh2s_no3 = 0.0
    endif
    if(ergom%id_jsul_o2 .gt. 0) then 
      allocate(ergom%jsul_o2(isd:ied,jsd:jed,nk))
      ergom%jsul_o2  = 0.0
    endif
    if(ergom%id_jsul_no3.gt. 0) then 
      allocate(ergom%jsul_no3(isd:ied,jsd:jed,nk))
      ergom%jsul_no3 = 0.0
    endif
    if(ergom%id_jrec_o2 .gt. 0) then 
      allocate(ergom%jrec_o2(isd:ied,jsd:jed,nk))
      ergom%jrec_o2  = 0.0
    endif
    if(ergom%id_jrec_no3.gt. 0) then
      allocate(ergom%jrec_no3(isd:ied,jsd:jed,nk))
      ergom%jrec_no3 = 0.0
    endif
    if(ergom%id_jrec_so4.gt. 0) then
      allocate(ergom%jrec_so4(isd:ied,jsd:jed,nk))
      ergom%jrec_so4 = 0.0
    endif
    if(ergom%id_jrec_ana.gt. 0) then
      allocate(ergom%jrec_ana(isd:ied,jsd:jed,nk))
      ergom%jrec_ana = 0.0
    endif
    do n=1, NUM_SPM
        vardesc_temp = vardesc(trim(spm(n)%name)//"_jsed", &
               "Sedimentation of "//trim(spm(n)%longname),'h','L','s','mol m-2 s-1','f')
        spm(n)%id_jsed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    enddo

    do n=1, NUM_SED
        vardesc_temp = vardesc(trim(sed(n)%name)//"_jgain_sed", &
               "Gain of "//trim(sed(n)%longname)//" by transformation of other sediment classes", &
	               'h','L','s','mol m-2 s-1','f')
        sed(n)%id_jgain_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc(trim(sed(n)%name)//"_jloss_sed", &
               "Loss of "//trim(sed(n)%longname)//" by transformation of sediment",'h','L','s','mol m-2 s-1','f')
        sed(n)%id_jloss_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc(trim(sed(n)%name)//"_jres", &
               "Resuspension of "//trim(sed(n)%longname),'h','L','s','mol m-2 s-1','f')
        sed(n)%id_jres = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc(trim(sed(n)%name)//"_jbiores", &
               "resuspension of "//trim(sed(n)%longname)//" by benthic organisms",'h','L','s','mol m-2 s-1','f')
        sed(n)%id_jbiores = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    enddo
    do n=1, NUM_DET
      vardesc_temp = vardesc(trim(det(n)%name)//"_jgraz_n", "Nitrogen loss to zooplankton by grazing", &
                       'h','L','s','mol m-2 s-1','f')
      det(n)%id_jgraz_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(det(n)%name)//"_jmort","Detritus mort. source layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      det(n)%jmort = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    enddo	  
    do n=1, NUM_ZOO
      vardesc_temp = vardesc(trim(zoo(n)%name)//"_jgraz_n","Nitrogen loss to zooplankton by grazing", &
                       'h','L','s','mol m-2 s-1','f')
      zoo(n)%id_jgraz_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(zoo(n)%name)//"_jgain_n","Grazing nitrogen uptake layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      zoo(n)%id_jgain_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(zoo(n)%name)//"_jres_n","Respiration nitrogen loss layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      zoo(n)%id_jres_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(zoo(n)%name)//"_jdet_n","Zooplankton nitrogen loss to detritus layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      zoo(n)%id_jdet_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    enddo	  
    do n=1, NUM_PHYTO
      vardesc_temp = vardesc(trim(phyto(n)%name)//"_jgraz_n","Grazing nitrogen uptake layer integral",&
                       'h','L','s','mol m-2 s-1','f')
      phyto(n)%id_jgraz_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(phyto(n)%name)//"_jres_n","Respiration nitrogen loss layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      phyto(n)%id_jres_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(phyto(n)%name)//"_jdet_n","phytoplankton nitrogen loss to detritus layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      phyto(n)%id_jdet_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(phyto(n)%name)//"_jprod_po4","phytoplankton nitrate uptake layer integral", &
                       'h','L','s','mol m-2 s-1','f')
      phyto(n)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      vardesc_temp = vardesc(trim(phyto(n)%name)//"_ilim","Light limitation",'h','L','s','W m-2','f')
      phyto(n)%id_ilim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	 init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      if (phyto(n)%id_ilim .gt. 0) then 
         allocate(phyto(n)%ilim(isd:ied,jsd:jed,nk))   
	 phyto(n)%ilim   = 0.0	
      endif 
      if (n .ne. cya) then
        vardesc_temp = vardesc(trim(phyto(n)%name)//"_jprod_no3","phytoplankton nitrate uptake layer integral", &
	               'h','L','s','mol m-2 s-1','f')
        phyto(n)%id_jprod_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	   init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
        vardesc_temp = vardesc(trim(phyto(n)%name)//"_jprod_nh4","phytoplankton nitrate uptake layer integral", &
	               'h','L','s','mol m-2 s-1','f')
        phyto(n)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
	   init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      endif
    enddo	  
    vardesc_temp = vardesc("jrec_n_sed","nitrogen loss by mineralisation",'h','L','s','mol m-2 s-1','f')
    biosed%id_jrec_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jdenit_sed","nitrogen loss by denitrification in sediment",'h','L','s','mol m-2 s-1','f')
    biosed%id_jdenit_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("mode_sed","sulfur bacteria on sediment",'h','L','s','mol m-2 s-1','f')
    biosed%id_mode_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    if (biosed%id_jrec_n .gt. 0) then    
       allocate(biosed%jrec_n(isd:ied,jsd:jed))    
       biosed%jrec_n  = 0.0
    endif
    if (biosed%id_jdenit_sed .gt. 0) then 
       allocate(biosed%jdenit_sed(isd:ied,jsd:jed)) 
       biosed%jdenit_sed  = 0.0
    endif
    if (biosed%id_mode_sed .gt. 0) then  
       allocate(biosed%mode_sed(isd:ied,jsd:jed))  
       biosed%mode_sed  = 0.0
    endif
    call user_register_2d_tracers

  end subroutine generic_ERGOM_register_diag

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params
  
  !Specify all parameters used in this modules.
  !=============================================================================
  !User adds one call for each parameter below!
  !User also adds the definition of each parameter in generic_ERGOM_params type
  !=============================================================================    

  !=============================================================================
  !Block Starts: g_tracer_add_param
  !=============================================================================
    
  !Add the known experimental parameters used for calculations
  !in this module.
  !All the g_tracer_add_param calls must happen between 
  !g_tracer_start_param_list and g_tracer_end_param_list  calls.
  !This implementation enables runtime overwrite via field_table.
  
    integer                      :: n
    real                         :: pref_sum
    character(len=fm_string_len) :: mystring
    integer :: stdoutunit

    stdoutunit=stdout()
    call g_tracer_start_param_list(package_name)
    
    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !-----------------------------------------------------------------------
    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    !  Until we know better take coefficient for oxygen
    call g_tracer_add_param('a1_nit', ergom%a1_nit, 2206.1)
    call g_tracer_add_param('a2_nit', ergom%a2_nit, -144.86)
    call g_tracer_add_param('a3_nit', ergom%a3_nit, 4.5413)
    call g_tracer_add_param('a4_nit', ergom%a4_nit, -0.056988)
    !---------------------------------------------------------------------
    !  Compute the Schmidt number of O2 in seawater using the 
    !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
    !  Cycles, 12, 141-163).
    !---------------------------------------------------------------------
    !New Wanninkhof numbers
    call g_tracer_add_param('a1_o2',  ergom%a1_o2, 1929.7)
    call g_tracer_add_param('a2_o2',  ergom%a2_o2, -117.46)
    call g_tracer_add_param('a3_o2',  ergom%a3_o2, 3.116)
    call g_tracer_add_param('a4_o2',  ergom%a4_o2, -0.0306)
    !-----------------------------------------------------------------------
    !     Solubility coefficients for alpha in µmol/kg/atm 
    !     Hamme and Emerson, 2004
    !     Hamme, R.C. and Emerson, S.R., 2004. The solubility of neon, nitrogen and argon 
    !     in distilled water and seawater. 
    !     Deep Sea Research Part I: Oceanographic Research Papers, 51(11): 1517-1528.
    !-----------------------------------------------------------------------
    call g_tracer_add_param('RHO_0',  ergom%Rho_0, 1035.0)
    call g_tracer_add_param('t1_nit', ergom%t1_nit, 6.42931)
    call g_tracer_add_param('t2_nit', ergom%t2_nit, 2.92704)
    call g_tracer_add_param('t3_nit', ergom%t3_nit, 4.32531)
    call g_tracer_add_param('t4_nit', ergom%t4_nit, 4.69149)
    call g_tracer_add_param('e1_nit', ergom%e1_nit, 298.15)
    call g_tracer_add_param('s1_nit', ergom%s1_nit, 7.44129e-3)
    call g_tracer_add_param('s2_nit', ergom%s2_nit, 8.02566e-3)
    call g_tracer_add_param('s3_nit', ergom%s3_nit, 0.0146775)

    !-----------------------------------------------------------------------
    !     Solubility coefficients for alpha in µmol/kg/atm FOR OXYGEN
    call g_tracer_add_param('t0_o2',  ergom%t0_o2, 2.00907)
    call g_tracer_add_param('t1_o2',  ergom%t1_o2, 3.22014)
    call g_tracer_add_param('t2_o2',  ergom%t2_o2, 4.05010)
    call g_tracer_add_param('t3_o2',  ergom%t3_o2, 4.94457)
    call g_tracer_add_param('t4_o2',  ergom%t4_o2, -2.56847e-01)
    call g_tracer_add_param('t5_o2',  ergom%t5_o2, 3.88767)
    call g_tracer_add_param('b0_o2',  ergom%b0_o2, -6.24523e-03)
    call g_tracer_add_param('b1_o2',  ergom%b1_o2, -7.37614e-03)
    call g_tracer_add_param('b2_o2',  ergom%b2_o2, -1.03410e-02 )
    call g_tracer_add_param('b3_o2',  ergom%b3_o2, -8.17083e-03)
    call g_tracer_add_param('c0_o2',  ergom%c0_o2, -4.88682e-07)

    !-----------------------------------------------------------------------
       
    ! phytoplankton parameters
        
    do n=1, NUM_PHYTO
      call g_tracer_add_param(trim(name_phyt(n))//'imin',   phyto(n)%imin, imin(n))           
      call g_tracer_add_param(trim(name_phyt(n))//'tmin',   phyto(n)%tmin, tmin(n))           
      call g_tracer_add_param(trim(name_phyt(n))//'smin',   phyto(n)%smin, smin(n))           
      call g_tracer_add_param(trim(name_phyt(n))//'smax',   phyto(n)%smax, smax(n))           

      call g_tracer_add_param(trim(name_phyt(n))//'alpha',  phyto(n)%alpha, alpha(n))         
      call g_tracer_add_param(trim(name_phyt(n))//'talpha', phyto(n)%talpha, talpha(n))      
                                                                                             
      call g_tracer_add_param(trim(name_phyt(n))//'rp0',    phyto(n)%rp0, rp0(n)/sperd)       
      call g_tracer_add_param(trim(name_phyt(n))//'p0' ,    phyto(n)%p0,  p0(n))              
     
      call g_tracer_add_param(trim(name_phyt(n))//'lpd',    phyto(n)%lpd, lpd(n)/sperd)       
      call g_tracer_add_param(trim(name_phyt(n))//'lpr',    phyto(n)%lpr, lpr(n)/sperd)       

      call g_tracer_add_param(trim(name_phyt(n))//'pnr',    phyto(n)%pnr, np(n)/nn(n))        
      call g_tracer_add_param(trim(name_phyt(n))//'cnr',    phyto(n)%cnr, nc(n)/nn(n))        
     
      call g_tracer_add_param(trim(name_phyt(n))//'wsink0', phyto(n)%wsink0, sinkp(n)/sperd)  
      call g_tracer_add_param(trim(name_phyt(n))//'name' ,  phyto(n)%name, trim(name_phyt(n))) 
    enddo

    ! zooplankton parameters

    do n=1, NUM_ZOO
      call g_tracer_add_param(trim(name_zoo(n))//'nue'          , zoo(n)%nue          , nue(n)/sperd) 
      call g_tracer_add_param(trim(name_zoo(n))//'food_to_nh4'  , zoo(n)%food_to_nh4  , food_to_nh4(n))   
      call g_tracer_add_param(trim(name_zoo(n))//'food_to_det'  , zoo(n)%food_to_det  , food_to_det(n))   
      call g_tracer_add_param(trim(name_zoo(n))//'food_to_nh4_2', zoo(n)%food_to_nh4_2, food_to_nh4_2(n)) 
      call g_tracer_add_param(trim(name_zoo(n))//'food_to_det_2', zoo(n)%food_to_det_2, food_to_det_2(n)) 
      call g_tracer_add_param(trim(name_zoo(n))//'graz',          zoo(n)%graz         , graz(n)/sperd) 
      call g_tracer_add_param(trim(name_zoo(n))//'graz_pref',     zoo(n)%graz_pref    , graz_pref(n)) 
      call g_tracer_add_param(trim(name_zoo(n))//'z0' ,           zoo(n)%z0           , z0(n))            
      call g_tracer_add_param(trim(name_zoo(n))//'sigma_b',       zoo(n)%sigma_b      , sigma_b(n)/sperd) 
      call g_tracer_add_param(trim(name_zoo(n))//'t_opt',         zoo(n)%t_opt,         t_opt_zoo(n))     
      call g_tracer_add_param(trim(name_zoo(n))//'t_max',         zoo(n)%t_max,         t_max_zoo(n))     
      call g_tracer_add_param(trim(name_zoo(n))//'beta',          zoo(n)%beta,          beta_zoo(n))      
      call g_tracer_add_param(trim(name_zoo(n))//'oxy_sub',       zoo(n)%oxy_sub      , oxy_sub_zoo(n))   
      call g_tracer_add_param(trim(name_zoo(n))//'oxy_min',       zoo(n)%oxy_min      , oxy_min_zoo(n))   
      call g_tracer_add_param(trim(name_zoo(n))//'resp_red',      zoo(n)%resp_red     , resp_red_zoo(n))  
      call g_tracer_add_param(trim(name_zoo(n))//'iv',    zoo(n)%iv,       iv     (n))       
      call g_tracer_add_param(trim(name_zoo(n))//'zcl1',  zoo(n)%zcl1,     zcl1   (n))       
      call g_tracer_add_param(trim(name_zoo(n))//'Imax',  zoo(n)%Imax, Imax(n))              
      call g_tracer_add_param(trim(name_zoo(n))//'vertical_migration', zoo(n)%vertical_migration, vertical_migration(n))           
      call g_tracer_add_param(trim(name_zoo(n))//'blanchard_temperature', zoo(n)%blanchard_temperature, blanchard_temperature(n))  
      call g_tracer_add_param(trim(name_zoo(n))//'alpha', zoo(n)%alpha, alpha_zoo(n))        
      call g_tracer_add_param(trim(name_zoo(n))//'wsink0',zoo(n)%wsink0, sinkz(n)/sperd)     
      call g_tracer_add_param(trim(name_zoo(n))//'wrise0',zoo(n)%wrise0, risez(n)/sperd)     
      call g_tracer_add_param(trim(name_zoo(n))//'vdiff_max',zoo(n)%vdiff_max, vdiff_max(n)) 
      call g_tracer_add_param(trim(name_zoo(n))//'dark_rise',zoo(n)%dark_rise,dark_rise(n))  
      call g_tracer_add_param(trim(name_zoo(n))//'wfood' ,zoo(n)%wfood,  wfood(n))           
      call g_tracer_add_param(trim(name_zoo(n))//'o2min' ,zoo(n)%o2min,  o2min(n))	
      call g_tracer_add_param(trim(name_zoo(n))//'h2smax',zoo(n)%h2smax, h2smax(n))  	 
      call g_tracer_add_param(trim(name_zoo(n))//'wtemp', zoo(n)%wtemp,  wtemp(n))	     
      call g_tracer_add_param(trim(name_zoo(n))//'wo2' ,  zoo(n)%wo2,    wo2(n))	     
      call g_tracer_add_param(trim(name_zoo(n))//'wh2s' , zoo(n)%wh2s,   wh2s(n))	     
      call g_tracer_add_param(trim(name_zoo(n))//'name'  ,zoo(n)%name, trim(name_zoo(n)))    
      call g_tracer_add_param(trim(name_zoo(n))//'pnr',   zoo(n)%pnr, np_zoo(n)/nn_zoo(n))
      call g_tracer_add_param(trim(name_zoo(n))//'cnr',   zoo(n)%cnr, nc_zoo(n)/nn_zoo(n))
      allocate(zoo(n)%pref_phy(NUM_PHYTO)) ; zoo(n)%pref_phy = 0.0	  
      allocate(zoo(n)%pref_zoo(NUM_ZOO))   ; zoo(n)%pref_zoo = 0.0	  
      allocate(zoo(n)%pref_det(NUM_DET))   ; zoo(n)%pref_det = 0.0	  
!
!     The sum of all food weights should be 1      
      pref_sum = 0.
      do m=1, NUM_PHYTO
         pref_sum = pref_sum + pref_phy(n,m)
      enddo
      do m=1, NUM_ZOO
         pref_sum = pref_sum + pref_zoo(n,m)
      enddo
      do m=1, NUM_DET
         pref_sum = pref_sum + pref_zoo(n,m)
      enddo
      do m=1, NUM_PHYTO
         pref_phy(n,m) = pref_phy(n,m)/(pref_sum+epsln)
      enddo
      do m=1, NUM_ZOO
         pref_zoo(n,m) = pref_zoo(n,m)/(pref_sum+epsln)
      enddo
      do m=1, NUM_DET
         pref_zoo(n,m) = pref_zoo(n,m)/(pref_sum+epsln)
      enddo
      if (pref_sum .le. epsln) then 
         write (stdoutunit,'(a)') 'WARNING, all preferences of '//trim(zoo(n)%name)//' are zero.'
      endif
      do m=1, NUM_PHYTO
         call g_tracer_add_param(trim(zoo(n)%name)//'pref'//trim(phyto(m)%name), zoo(n)%pref_phy(m), pref_phy(n,m))
      enddo
      do m=1, NUM_ZOO
         call g_tracer_add_param(trim(zoo(n)%name)//'pref'//trim(zoo(m)%name),   zoo(n)%pref_zoo(m), pref_zoo(n,m))
      enddo
      do m=1, NUM_DET
         call g_tracer_add_param(trim(zoo(n)%name)//'pref'//trim(det(m)%name),   zoo(n)%pref_det(m), pref_det(n,m))
      enddo
    enddo
    
  
    ! detritus parameters    
    if (NUM_DET .eq. 1) FAST = SLOW
    do n=1, NUM_DET
    ! Detritus parameters for recycling and sinking
      call g_tracer_add_param(trim(name_det(n))//'dn',	   det(n)%dn,      dn   (n)/sperd)  
      call g_tracer_add_param(trim(name_det(n))//'q10_rec',det(n)%q10_rec, q10_rec(n)) 
      call g_tracer_add_param(trim(name_det(n))//'name'  , det(n)%name,    trim(name_det(n)))    
    enddo
      
    ! sediment parameters
     
    call g_tracer_add_param('dn'      ,  biosed%dn,       dn_sed/sperd)      
    call g_tracer_add_param('frac_dn_anoxic',biosed%frac_dn_anoxic,frac_dn_anoxic)
    call g_tracer_add_param('q10_rec',   biosed%q10_rec,  q10_rec_sed )      
    call g_tracer_add_param('thio_bact_min',  biosed%thio_bact_min,  thio_bact_min )  
    call g_tracer_add_param('pnr'     ,  biosed%pnr     , np_sed/nn_sed)     
    call g_tracer_add_param('cnr'     ,  biosed%cnr     , nc_sed/nn_sed)     
    call g_tracer_add_param('den_rate',  biosed%den_rate, den_rate)         
    call g_tracer_add_param('po4_lib_rate',    biosed%po4_lib_rate,    po4_lib_rate/sperd) 
    call g_tracer_add_param('po4_retention',   biosed%po4_retention,   po4_retention) 
    call g_tracer_add_param('po4_ret_plus_BB', biosed%po4_ret_plus_BB, po4_ret_plus_BB)
    call g_tracer_add_param('o2_bioerosion',   biosed%o2_bioerosion,   o2_bioerosion) 

    ! suspended particulate matter parameters
    do n=1, NUM_SPM
        call g_tracer_add_param(trim(name_spm(n))//'name',        spm(n)%name,        name_spm(n))      
        call g_tracer_add_param(trim(name_spm(n))//'wsink0',      spm(n)%wsink0,      wsink0_spm(n)/sperd) 
        call g_tracer_add_param(trim(name_spm(n))//'wsed',        spm(n)%wsed,        wsed_spm(n)/sperd)   
        call g_tracer_add_param(trim(name_spm(n))//'sediment_to', spm(n)%sediment_to, sediment_to(n))  
    enddo
    ! settled matter parameters
    do n=1, NUM_SED
        call g_tracer_add_param(trim(name_sed(n))//'name_2d',     sed(n)%name,        name_sed(n))      
        call g_tracer_add_param(trim(name_sed(n))//'suspend_to',  sed(n)%suspend_to,  suspend_to(n))   
        call g_tracer_add_param(trim(name_sed(n))//'erosion_rate',sed(n)%erosion_rate,& 
	                             erosion_rate_sed(n)/sperd) 
        call g_tracer_add_param(trim(name_sed(n))//'bioerosion_rate',sed(n)%bioerosion_rate, &
	                             bioerosion_rate_sed(n)/sperd) 
        call g_tracer_add_param(trim(name_sed(n))//'molar_volume',sed(n)%molar_volume,molar_volume_sed(n)) 
        call g_tracer_add_param(trim(name_sed(n))//'critical_stress',sed(n)%critical_stress, &
	                             critical_stress_sed(n))
      call g_tracer_add_param(trim(longname_sed(n))//'longname',   sed(n)%longname,    longname_sed(n))       
    enddo

    ! sediment settings parameters
    do n=1,NUM_SEDIMENT_LAYERS
      write( mystring, '(i4)' )  n
      ! (maximum) height of vertical layers [m]. Numbers <0 mean the layer may become infinitely thick.
      call g_tracer_add_param('sed_layer_height_'//trim(adjustl(mystring)), sed_defs%layer_height(n), & 
                                     sed_layer_height(n))
    enddo
    call g_tracer_add_param('sed_layer_propagation', sed_defs%layer_propagation, sed_layer_propagation) 
    call g_tracer_add_param('sed_erosion_mode', sed_defs%erosion_mode     , sed_erosion_mode     ) 
    call g_tracer_add_param('NUM_SEDIMENT_LAYERS'  , sed_defs%NUM_LAYERS       , NUM_SEDIMENT_LAYERS	   ) 

    call g_tracer_add_param('q10_nit',  ergom%q10_nit,   q10_nit)      
    call g_tracer_add_param('nf',       ergom%nf,        nf/sperd)     
    call g_tracer_add_param('alpha_nit',ergom%alpha_nit, alpha_nit)    

    ! parameter for sulfur cycle elements
        
    call g_tracer_add_param('q10_h2s'  , ergom%q10_h2s,   q10_h2s)             
    call g_tracer_add_param('k_h2s_o2' , ergom%k_h2s_o2,  k_h2s_o2/sperd)      
    call g_tracer_add_param('k_h2s_no3', ergom%k_h2s_no3, k_h2s_no3/sperd)     
    call g_tracer_add_param('k_sul_o2' , ergom%k_sul_o2,  k_sul_o2/sperd)      
    call g_tracer_add_param('k_sul_no3', ergom%k_sul_no3, k_sul_no3/sperd)     
    call g_tracer_add_param('alp_o2'   , ergom%alp_o2,    alp_o2)              
    call g_tracer_add_param('alp_h2s'  , ergom%alp_h2s,   alp_h2s)             
    call g_tracer_add_param('alp_no3'  , ergom%alp_no3,   alp_no3)             
    call g_tracer_add_param('alp_nh4'  , ergom%alp_nh4,   alp_nh4)             

    ! parameter for anammox

    call g_tracer_add_param('k_an0'    , ergom%k_an0,     k_an0/sperd)         

    

    call g_tracer_end_param_list(package_name)
    
    !=====================================================================
    !Block Ends: g_tracer_add_param
    !=====================================================================



  end subroutine user_add_params

  !
  !   This is an internal sub, not a public interface.
  !   Add all the tracers to be used in this module. 
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'
    integer                                 :: i, n
    character(len=fm_string_len)            :: mystring

    call g_tracer_start_param_list(package_name)
    call g_tracer_add_param('ice_restart_file'   , ergom%ice_restart_file   , 'ice_ergom.res.nc')
    call g_tracer_add_param('ocean_restart_file' , ergom%ocean_restart_file , 'ocean_ergom.res.nc' )
    call g_tracer_add_param('IC_file'            , ergom%IC_file       , '')
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file=ergom%ice_restart_file, ocean_restart_file=ergom%ocean_restart_file )

    !=====================================================
    !Specify all prognostic tracers of this modules.
    !=====================================================
    !User adds one call for each prognostic tracer below!
    !User should specify if fluxes must be extracted from boundary 
    !by passing one or more of the following methods as .true.  
    !and provide the corresponding parameters array
    !methods: flux_gas,flux_runoff,flux_wetdep,flux_drydep  
    !
    !prog_tracers: nitrogen 
    !diag_tracers: none
    !

    call user_init_2d_tracer_list

    do n=1, NUM_SPM
        call g_tracer_add(tracer_list,package_name,&
           name       = trim(name_spm(n)),             &
           longname   = trim(longname_spm(n))//' concentration in water', &
           units      = 'mol/kg',          &
	   prog       = .true. ,           &
           flux_bottom= .true.            )
    enddo
    do n=1, NUM_SED
        do i=1,NUM_SEDIMENT_LAYERS
          write( mystring, '(i4)' )  i
          call user_add_2d_tracer(tracer_list,    &
             name       = trim(name_sed(n))//'_'//trim(adjustl(mystring)),  &
             longname   = trim(longname_sed(n))//' concentration in sediment layer '//trim(adjustl(mystring)), &
             units      = 'mol/m^2')
        enddo
    enddo

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nitrogen',               &
         longname   = 'nitrogen Concentration', &
         units      = 'mol/kg',            &
         prog       = .true.,               &
         flux_gas       = .true.,           &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /), &
         flux_gas_restart_file  = 'ocean_ergom_airsea_flux.res.nc', &
         flux_bottom= .true.             )

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

    do n=1, NUM_PHYTO
       call g_tracer_add(tracer_list,package_name,&
         name       = trim(name_phyt(n)),        &
         longname   = trim(name_phyt(n))//' Nitrogen',&
         units      = 'mol/kg',          &
         prog       = .true.,            &
         btm_reservoir = .true.          )
    enddo

    call g_tracer_add(tracer_list,package_name,&
         name       = 'chl',         &
         longname   = 'Chlorophyll', &
         units      = 'ug kg-1',     &
         prog       = .false.,       &
         init_value = 0.08           )

    do n=1, NUM_ZOO
       call g_tracer_add(tracer_list,package_name,&
         name       = trim(name_zoo(n)),	  &
         longname   = trim(name_zoo(n))//' Nitrogen',	  &
         units      = 'mol/kg',                   &
         prog       = .true.  ,                   &
	 move_vertical = vertical_migration(n),   &
	 diff_vertical = vertical_migration(n))
    enddo

!!    call g_tracer_add(tracer_list,package_name,&
!!         name       = 'ndet_btf',            &
!!         longname   = 'N flux to Sediments', &
!!         units      = 'mol m-2 s-1',         &
!!         prog       = .false.                )

!    call g_tracer_add(tracer_list,package_name,&
!         name       = 'sed',             &
!         longname   = 'Sediment',        &
!         units      = 'mol/m^2',         &
!         prog       = .false.            )              
         

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nh4',             &
         longname   = 'Ammonium',         &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_runoff= .true.,            &
         flux_wetdep= .true.,            &
         flux_drydep= .true.,            &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )
    
    call g_tracer_add(tracer_list,package_name,&
         name       = 'po4',             &
         longname   = 'Phosphate',         &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_runoff= .true.,            &
         flux_wetdep= .true.,            &
         flux_drydep= .true.,            &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )
	 
    call g_tracer_add(tracer_list,package_name,&
         name       = 'h2s',             &
         longname   = 'Sulfide',         &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_runoff= .false.,           &
         flux_wetdep= .false.,           &
         flux_drydep= .false.,           &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )

    call g_tracer_add(tracer_list,package_name,                        &
         name       = 'o2',                                            &
         longname   = 'Oxygen',                                        &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
!         flux_gas_name  = 'o2_flux',                                   &
         flux_gas_molwt = WTMO2,                                       &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                  &
         flux_gas_restart_file  = 'ocean_ergom_airsea_flux.res.nc',    &
         flux_bottom= .true.             )
	 
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sul',             &
         longname   = 'Sulfur',          &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_runoff= .false.,           &
         flux_wetdep= .false.,           &
         flux_drydep= .false.,           &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .false.             )

  end subroutine user_add_tracers

  !
  !   This is an internal sub, not a public interface.
  !   initialize the list of 2d tracers with a length of zero
  !
  subroutine user_init_2d_tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'user_init_2d_tracer_list'
    allocate(tracers_2d(0))
  end subroutine user_init_2d_tracer_list
  !
  !   This is an internal sub, not a public interface.
  !   Insert a 2d tracer into the 2d tracer list
  !   2d tracers are stored  together in diagnostic 3d tracers.
  !   Those are registered for output and in the restart list.
  !
  subroutine user_add_2d_tracer(tracer_list,name,longname,units)
    type(g_tracer_type), pointer :: tracer_list
    character(len=*), intent(in) :: name, longname, units

    character(len=fm_string_len), parameter :: sub_name = 'user_add_2d_tracer'
    integer                      :: m, n, dummy
    character(len=fm_string_len) :: temp_string
    type(tracer_2d), ALLOCATABLE, dimension(:) :: temp_tracers_2d
    
    n  = size(tracers_2d)
    allocate(temp_tracers_2d(n))
    do m=1,n
      temp_tracers_2d(m)=tracers_2d(m)
    end do
    deallocate(tracers_2d)
    allocate(tracers_2d(n+1))
    do m=1,n
      tracers_2d(m)=temp_tracers_2d(m)
    end do
    deallocate(temp_tracers_2d)
    
    tracers_2d(n+1)%name     = trim(adjustl(name))
    tracers_2d(n+1)%longname = trim(adjustl(longname))
    tracers_2d(n+1)%units    = trim(adjustl(units))
    tracers_2d(n+1)%layer_in_3d_tracer = mod(n,vlev_sed)+1
    tracers_2d(n+1)%field_assigned     = .false.  ! whether a 2d field has been assigned
    
    m = n/vlev_sed+1
    write( temp_string, '(i4)' )  m; temp_string = adjustl(temp_string)
    tracers_2d(n+1)%name_of_3d_tracer = 'tracer_2d_'//trim(temp_string)
    if (mod(n,vlev_sed) .eq. 0) then
      call g_tracer_add(tracer_list,package_name,                              &
         name       = 'tracer_2d_'//trim(temp_string),                         &
         longname   = 'Tracer '//trim(temp_string)//' containing 2d variables',&
         units      = 'none',                                                  &
         prog       = .false.)
    endif
  end subroutine user_add_2d_tracer

  !   This is an internal sub, not a public interface.
  !   Register the 2d tracer for output via a 2d output field
  subroutine user_register_2d_tracers
    real,parameter :: missing_value1=-1.0e+20
    type(vardesc)  :: vardesc_temp
    integer        :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, axes(3)
    type(time_type):: init_time 
    integer        :: m,n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes=axes,init_time=init_time) 

    n=size(tracers_2d)
    do m=1,n
      vardesc_temp = vardesc(                                                                   &
      tracers_2d(m)%name, tracers_2d(m)%longname,'h','L','s',tracers_2d(m)%units,'f')
      tracers_2d(m)%diag_field = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
	      init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    end do
  end subroutine user_register_2d_tracers

  !   This is an internal sub, not a public interface.
  !   Let the p_field pointer of the 2d tracer point to the data array (for read and write)
  subroutine user_2d_tracer_assign_array(name,array)
    character(len=*), intent(in)           :: name
    real,dimension(:,:),target             :: array
    integer                                :: n
    logical                                :: found_tracer
    character(len=fm_string_len)           :: errorstring

    found_tracer=.false.
    do n=1,size(tracers_2d)
      if (trim(adjustl(name)) .eq. trim(adjustl(tracers_2d(n)%name))) then
        tracers_2d(n)%p_field => array
        tracers_2d(n)%field_assigned = .true.
        found_tracer=.true.
      end if
    end do
    if (.not. found_tracer) then
       write(errorstring, '(a)') &
       'array assigned to tracer '//trim(name)// &
       ', but that tracer was not added by user_add_2d_tracer'
       call  mpp_error(FATAL, errorstring)
    end if
  end subroutine user_2d_tracer_assign_array

  !   This is an internal sub, not a public interface.
  !   Load all data stored in the (3d tracers containing 2d tracer data) into their temporary 2d arrays
   subroutine user_get_2d_tracer_values(tracer_list,isd,ied,jsd,jed,nk)
    type(g_tracer_type),    pointer      :: tracer_list
    integer,                  intent(in) :: isd,ied,jsd,jed,nk
    real,dimension(:,:,:),allocatable    :: a3d
    integer                              :: n
    character(len=fm_string_len)         :: loaded_3d_tracer
    allocate(a3d(isd:ied,jsd:jed,1:nk))
    loaded_3d_tracer=''
    do n=1,size(tracers_2d)
      if (tracers_2d(n)%field_assigned) then
        if (trim(adjustl(tracers_2d(n)%name_of_3d_tracer)) .ne. trim(adjustl(loaded_3d_tracer))) then
          call g_tracer_get_values(     &
	  tracer_list,tracers_2d(n)%name_of_3d_tracer,'field',a3d,isd,jsd,ntau=1,positive=.true.)
          loaded_3d_tracer=tracers_2d(n)%name_of_3d_tracer
        end if
        tracers_2d(n)%p_field = a3d(:,:,tracers_2d(n)%layer_in_3d_tracer)
      end if
    end do
    deallocate(a3d)
  end subroutine user_get_2d_tracer_values

  !   This is an internal sub, not a public interface.
  !   Save all 2d tracer data from their temporary 2d arrays into the 3d tracers containing 2d tracer data
  subroutine user_set_2d_tracer_values(tracer_list ,model_time)
    type(g_tracer_type),      pointer         :: tracer_list
    type(time_type),          intent(in)      :: model_time
    real, dimension(:,:,:),   pointer         :: grid_tmask
    
    integer                                   :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau
    real,dimension(:,:,:),allocatable         :: a3d
    integer                                   :: n, layer
    character(len=fm_string_len)              :: loaded_3d_tracer
    logical                                   :: used

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    allocate(a3d(isd:ied,jsd:jed,1:nk))
    loaded_3d_tracer=''
    do n=1,size(tracers_2d)
      if (tracers_2d(n)%field_assigned) then
        if (trim(adjustl(tracers_2d(n)%name_of_3d_tracer)) .ne. trim(adjustl(loaded_3d_tracer))) then 
          ! a new 3d tracer starts -> first save the temporarily stored data into the old 3d tracer
          if (loaded_3d_tracer .ne. '') then
            call g_tracer_set_values(tracer_list,loaded_3d_tracer,'field',a3d,isd,jsd,ntau=1)
          end if
          ! then, nullify the temporary tracer
          a3d=0.0
          loaded_3d_tracer=tracers_2d(n)%name_of_3d_tracer
        end if
        ! save the data to the temporary tracer
        layer = tracers_2d(n)%layer_in_3d_tracer
        a3d(:,:,layer)  = tracers_2d(n)%p_field
        ! save the values to the diagnostic field
        used = send_data(tracers_2d(n)%diag_field, a3d(:,:,layer),      & 
                  model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      end if
    end do
    !now, save the last tracer
    if (loaded_3d_tracer .ne. '') then
      call g_tracer_set_values(tracer_list,loaded_3d_tracer,'field',a3d,isd,jsd,ntau=1)
    end if
    deallocate(a3d)
  end subroutine user_set_2d_tracer_values

  ! <SUBROUTINE NAME="generic_ERGOM_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for CFCs.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_ERGOM_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_update_from_coupler'
    !
    !Nothing specific to be done for CFC's
    !
    return
  end subroutine generic_ERGOM_update_from_coupler

  ! <SUBROUTINE NAME="sedimentation_and_resuspension">
  !  <OVERVIEW>
  !    Perform sedimentation and resuspension for all spm and sed tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    All tracers that are able to be sedimented are stored in the spm array.
  !    All tracers that are able to be resuspended are stored in the sed array.
  !    This subroutine performs the sedimentation and resuspension.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call sedimentation_and_resuspension(NUM_SPM, spm, NUM_SED, sed,    &
  !isc, iec, jsc, jec, isd, ied, jsd, jed, grid_kmt, dzt, rho_dzt, tau, dt, &
  !sed_defs, current_wave_stress, bioerosion)
  !  </TEMPLATE>
  ! </SUBROUTINE>
  subroutine sedimentation_and_resuspension(NUM_SPM, spm, NUM_SED, sed,    &
  isc, iec, jsc, jec, isd, ied, jsd, jed, grid_kmt, dzt, rho_dzt, tau, dt, &
  sed_defs, current_wave_stress, bioerosion)
    use mpp_mod,                only: mpp_pe
    integer, intent(in)                                     :: NUM_SPM,        &
                                                               NUM_SED,        &
                                                               isc, iec, jsc, jec, &
							       isd, ied, jsd, jed, tau
    type(spm_type), intent(inout), dimension(:)             :: spm
    type(sed_type), intent(inout), dimension(:)             :: sed
    integer, dimension(isd:,jsd:)                           :: grid_kmt
    real, intent(in), dimension(isd:,jsd:,:)                :: dzt, rho_dzt
    real, intent(in)                                        :: dt
    type(sed_defs_type), intent(inout)                      :: sed_defs
    real, dimension(isd:,jsd:),     intent(in)              :: current_wave_stress
    real, dimension(isd:,jsd:), intent(in)                  :: bioerosion


    
    integer                            :: i,j,k,l,n,n_sed
    real                               :: temp1, diff_stress
    real,  allocatable, dimension(:,:) :: work1
    real,  parameter                   :: lbi = 0.007/24/3600   ! 7 promille per day in [1/s]
    real,  parameter                   :: ironpo4norm = 0.0045  ! [mol/m/m]
    real,  parameter                   :: ironpo4th   = 0.1     ! [mol/m/m]
    real,  parameter                   :: smax = 4.5            ! maximum nitrogen in sediment
    
    call mpp_clock_begin(id_susp)
    
    do n=1,NUM_SPM
       spm(n)%btf = 0
    enddo !n
   
    !sedimentation
    do n=1,NUM_SPM
       if (spm(n)%index_sediment_to .gt. 0) then
          n_sed=spm(n)%index_sediment_to
          do j = jsc, jec; do i = isc, iec
             k=grid_kmt(i,j)
             if (k == 0) cycle
             temp1 = spm(n)%p_wat(i,j,k,tau) * spm(n)%wsed * ergom%rho_0             ! settling rate [mol m-2 s-1]
             spm(n)%btf(i,j) = spm(n)%btf(i,j) + temp1
             sed(n_sed)%f_sed(i,j,1)   = &
          		sed(n_sed)%f_sed(i,j,1) + temp1  * dt		      ! add suspended matter to sediment
             if (spm(n)%id_jsed .gt. 0) spm(n)%jsed(i,j) = temp1	      ! output of sedimentation
          enddo; enddo !i,j
       endif
    enddo !n

    !resuspension
    selectcase(sed_defs%erosion_mode)
    case(INDEPENDENT)
      do n=1,NUM_SED
          if (sed(n)%index_suspend_to .gt. 0) then
             n_sed=sed(n)%index_suspend_to
             do j = jsc, jec; do i = isc, iec
                k=grid_kmt(i,j)
                if (k == 0) cycle
                if (current_wave_stress(i,j) .gt. sed(n)%critical_stress) then
! constant erosion independent off stress. The more organic matter is in the
! sediment the more will be suspended. Inorganic matter is not considered here.
        	   temp1 = sed(n)%f_sed(i,j,1) * sed(n)%erosion_rate	   ! erosion rate [mol m-2 s-1]
        	   if (sed(n)%id_jres .gt. 0) sed(n)%jres(i,j) = temp1	   ! output of resuspension
        	   sed(n)%f_sed(i,j,1) = sed(n)%f_sed(i,j,1) - temp1  * dt ! remove sediment [mol m-2]
                   spm(n_sed)%btf(i,j) = spm(n_sed)%btf(i,j) - temp1
                else
        	   if (sed(n)%id_jres .gt. 0) sed(n)%jres(i,j) = 0.
                endif
             enddo; enddo !i,j
          endif
      enddo !n

    case(ORGANIC)
      do n=1,NUM_SED
         if (sed(n)%index_suspend_to .gt. 0) then
            n_sed=sed(n)%index_suspend_to
            do j = jsc, jec; do i = isc, iec
              k=grid_kmt(i,j)
              if (k == 0) cycle
	      diff_stress = current_wave_stress(i,j) - sed(n)%critical_stress
              if (diff_stress .gt. 0.0) then
        	 temp1 = diff_stress * sed(n)%erosion_rate	           ! erosion rate [mol m-2 s-1]
        	 if (sed(n)%id_jres .gt. 0) sed(n)%jres(i,j) = temp1	   ! output of resuspension
        	 sed(n)%f_sed(i,j,1) = sed(n)%f_sed(i,j,1) - temp1  * dt   ! remove sediment [mol m-2]
                 spm(n_sed)%btf(i,j) = spm(n_sed)%btf(i,j) - temp1
!        	 spm(n_sed)%p_wat(i,j,k,tau)   = &
!        		   spm(n_sed)%p_wat(i,j,k,tau) + temp1  * dt	   ! add concentration to spm
               else
        	  if (sed(n)%id_jres .gt. 0) sed(n)%jres(i,j) = 0.
               endif
            enddo; enddo !i,j
         endif
      enddo !n
    case(MAXSTRESS)
    case DEFAULT
       call  mpp_error(FATAL, &
      '==>Error from ERGOM, undefined erosion mode.')
    end select

    !bioresuspension
    do n=1,NUM_SED
       if (sed(n)%index_suspend_to .gt. 0) then
          n_sed=sed(n)%index_suspend_to
          do j = jsc, jec; do i = isc, iec
             k=grid_kmt(i,j)
             if (k == 0) cycle
             temp1 = sed(n)%f_sed(i,j,1) * sed(n)%bioerosion_rate * bioerosion(i,j) ! bioerosion [mol m-2 s-1]
             if (sed(n)%id_jbiores .gt. 0) sed(n)%jbiores(i,j) = temp1  	    ! output of bioerosion
             sed(n)%f_sed(i,j,1) = sed(n)%f_sed(i,j,1) - temp1  * dt		    ! remove sediment [mol m-2]
             spm(n_sed)%btf(i,j) = spm(n_sed)%btf(i,j) - temp1
          enddo; enddo !i,j
       endif
    enddo !n


    !layer propagation
    selectcase(sed_defs%layer_propagation)
    case(SLP_DOWNWARD)
      !if the surface box is overfull, the same volume propagates downward through all boxes
      allocate(work1(isc:iec,jsc:jec))
      do l=1, sed_defs%NUM_LAYERS
         !calculate the volume in the box
         work1 = 0.0
         do n=1, NUM_SED
            do j = jsc, jec; do i = isc, iec
               k = grid_kmt(i,j)
               if (k == 0) cycle
               work1(i,j) = work1(i,j) + sed(n)%f_sed(i,j,1)*sed(n)%molar_volume ![(mol/m2) * (m3/mol) = m]
            enddo; enddo !i,j
         enddo !n
         !calculate the ratio between current and maximal box height
         do j = jsc, jec; do i = isc, iec
            k=grid_kmt(i,j)
            if (k == 0) cycle
            work1(i,j) = max(work1(i,j)/sed_defs%layer_height(l),1.0)  ! if layer_height<0, nothing is transferred
         enddo; enddo !i,j
         !divide the amount in the current box by this ratio
         do n=1, NUM_SED
            do j = jsc, jec; do i = isc, iec
               k=grid_kmt(i,j)
               if (k == 0) cycle
               sed(n)%f_sed(i,j,1) = sed(n)%f_sed(i,j,1)/work1(i,j) 
            enddo; enddo !i,j
         enddo !n
         !calculate which fraction needs to be transferred downward
         do j = jsc, jec; do i= isc, iec
            k=grid_kmt(i,j)
            if (k == 0) cycle
            work1(i,j) = work1(i,j) - 1.0 !e.g. 0.5
         enddo; enddo !i,j
         !store everything in the box below
         if (l .lt. sed_defs%NUM_LAYERS) then
            do n=1, NUM_SED
               do j = jsc, jec; do i = isc, iec
         	  k=grid_kmt(i,j)
         	  if (k == 0) cycle
         	  sed(n)%f_sed(i,j,l+1) = sed(n)%f_sed(i,j,l+1) + sed(n)%f_sed(i,j,l)*work1(i,j)
               enddo; enddo !i,j
            enddo !n
         endif
      enddo !l
      deallocate(work1)
    case(SLP_FULL_BOX)
    case(SLP_OLD_ERGOM)
      !if the surface box is overfull, sediment detritus propagates down through all layers. 
      allocate(work1(isc:iec,jsc:jec))
      do l=1, NUM_SEDIMENT_LAYERS
         !calculate the volume in the box
         work1=0.0
         do n=1, 1 !NUM_SED
            do j = jsc, jec; do i = isc, iec
               k=grid_kmt(i,j)
               if (k == 0) cycle
               work1(i,j) = work1(i,j) + sed(n)%f_sed(i,j,1)*sed(n)%molar_volume ![(mol/m2) * (m3/mol) = m]
            enddo; enddo !i,j
         enddo !n
         !calculate the ratio between current and maximal box height
         do j = jsc, jec; do i = isc, iec
            k=grid_kmt(i,j)
            if (k == 0) cycle
            work1(i,j) = max(work1(i,j)/sed_defs%layer_height(l),1.0)  ! if layer_height<0, nothing is transferred
         enddo; enddo !i,j
         !divide the amount in the current box by this ratio
         do n=1, 1 !NUM_SED
            do j = jsc, jec; do i = isc, iec
               k=grid_kmt(i,j)
               if (k == 0) cycle
               sed(n)%f_sed(i,j,1) = sed(n)%f_sed(i,j,1)/work1(i,j) 
            enddo; enddo !i,j
         enddo !n
         !calculate which fraction needs to be transferred downward
         do j = jsc, jec; do i = isc, iec
            k=grid_kmt(i,j)
            if (k == 0) cycle
            work1(i,j) = work1(i,j) - 1.0 !e.g. 0.5
         enddo; enddo !i,j
         !store everything in the box below
         if (l .lt. sed_defs%NUM_LAYERS) then
            do n=1,1 !NUM_SED
               do j = jsc, jec; do i = isc, iec
         	  k=grid_kmt(i,j)
         	  if (k == 0) cycle
         	  sed(n)%f_sed(i,j,l+1) = sed(n)%f_sed(i,j,l+1) + sed(n)%f_sed(i,j,l)*work1(i,j)
               enddo; enddo !i,j
            enddo !n
         endif
      enddo !l
      ! now, store some iron phosphate from the uppermost box in the second one
      do j = jsc, jec; do i = isc, iec
         k=grid_kmt(i,j)
         if (k == 0) cycle
         temp1 = sed(2)%f_sed(i,j,1)  ! iron phosphate concentration [mol/m2]
         temp1 = max(temp1,temp1+temp1*(temp1-ironpo4th)/ironpo4norm)*lbi*temp1/smax ! burial of iron phosphate [mol/m2/s]
         sed(2)%f_sed(i,j,1) = sed(2)%f_sed(i,j,1) - temp1*dt                        ! remove from first layer [mol/m2]
         sed(2)%f_sed(i,j,2) = sed(2)%f_sed(i,j,2) + temp1*dt                        ! store in second layer [mol/m2]
      enddo; enddo
      deallocate(work1)
    
    case DEFAULT
       call  mpp_error(FATAL, &
      '==>Error from ERGOM, undefined sediment layer propagation.')
    end select
    call mpp_clock_end(id_susp)

  end subroutine sedimentation_and_resuspension

  ! <SUBROUTINE NAME="generic_ERGOM_update_from_bottom">
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracers have bottom fluxes and reservoirs. 
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_update_from_bottom(tracer_list,dt, tau, model_time) 
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
  subroutine generic_ERGOM_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:,:),ALLOCATABLE :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    !---------------------------------------------------------------------
    ! Get bottom flux of detritus and reset bottom flux boundary condition
    !---------------------------------------------------------------------
!!    call g_tracer_get_values(tracer_list,'det'  ,'btm_reservoir', ergom%fndet_btm,isd,jsd)
!!    ergom%fndet_btm = ergom%fndet_btm /dt
!!    call g_tracer_get_values(tracer_list,'det_btf'  ,'field', temp_field,isd,jsd)
!!    temp_field(:,:,1,tau) = ergom%fndet_btm(:,:)
!!    call g_tracer_set_values(tracer_list,'det_btf'  ,'field', temp_field,isd,jsd)
!!    call g_tracer_set_values(tracer_list,'det','btm_reservoir',0.0)
!!    if (ergom%id_fndet_btm .gt. 0)           &
!!         used = send_data(ergom%id_fndet_btm,    ergom%fndet_btm,          &
!!         model_time, rmask = grid_tmask(:,:,1),& 
!!         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

   end subroutine generic_ERGOM_update_from_bottom

  ! <SUBROUTINE NAME="generic_ERGOM_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for CFCs.
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_ERGOM_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,hblt_depth,&
       ilb,jlb,tau,dt,grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band, &
       current_wave_stress)

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
    real, dimension(ilb:,jlb:),     intent(in) :: current_wave_stress

    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_update_from_source'
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau, i, j, k , kblt, n, m
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast,grid_kmt
     
    logical                         :: used
    integer                         :: nb 
    real, parameter                 :: p5 = 0.5
    real                            :: ntemp, ptemp, ntemp_r, itemp, nlim, ilim, plim, tlim, Iopt, rp, jtemp,  &
                                       ttemp, stemp, temp_alpha, temp_talpha
    real                            :: toptz, temp1, temp2, temp3, tempq, teo2q, tominq, tosubq, graz_z, food, &
                                       food_pref, gg
    real                            :: lzn, lzd, lzd1, lzd2, zz0
    real                            :: o2_avail, no3_avail, h2s_avail, sul_avail, nh4_avail, q10
    real                            :: do2, dno3, dh2s, dsul, sca, o2_sul, o2_h2s, no3_sul, no3_h2s
    real                            :: det_recyc, det_left, rec_o2, rec_ana, rec_no3, rec_so4    
    real                            :: r_o2, r_no3, r_nh4, r_h2s, k_an         
    real                            :: den_rate, active_sed, pot_o2, dh2s_dt, do2_dt, dno3_dt, dsed
    real                            :: rhodztdt_r
    real, ALLOCATABLE, dimension(:,:,:) :: wrk1, wrk2, wrk3, wrk4

    real, dimension(:,:)   ,pointer :: runoff_flux_no3, wet_no3, dry_no3
    real, dimension(:,:)   ,pointer :: runoff_flux_nh4, wet_nh4, dry_nh4
    real, dimension(:,:)   ,pointer :: runoff_flux_po4, wet_po4, dry_po4

    
    call mpp_clock_begin(id_source)

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

    allocate(wrk1(isd:ied, jsd:jed, 1:nk));                wrk1=0.0
    allocate(wrk2(isd:ied, jsd:jed, 1:nk));                wrk2=0.0
    allocate(wrk3(isd:ied, jsd:jed, 1:nk));                wrk3=0.0
    allocate(wrk4(isd:ied, jsd:jed, 1:nk));                wrk4=0.0

    !---------------------------------------------------------------------
    ! Get positive tracer concentrations
    !---------------------------------------------------------------------
    call g_tracer_get_values(tracer_list,'no3'    ,'field',ergom%f_no3     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nh4'    ,'field',ergom%f_nh4     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'po4'    ,'field',ergom%f_po4     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'o2'     ,'field',ergom%f_o2      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'h2s'    ,'field',ergom%f_h2s     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sul'    ,'field',ergom%f_sul     ,isd,jsd,ntau=tau,positive=.true.)
    do n=1, NUM_PHYTO
       call g_tracer_get_values(tracer_list,trim(phyto(n)%name),'field',phyto(n)%f_n ,isd,jsd,ntau=tau,positive=.true.)
    enddo
    do n=1, NUM_ZOO
      call g_tracer_get_values(tracer_list,trim(zoo(n)%name)   ,'field',zoo(n)%f_n   ,isd,jsd,ntau=tau,positive=.true.)
    enddo
    do n=1, NUM_DET
      call g_tracer_get_values(tracer_list,trim(det(n)%name)   ,'field',det(n)%f_n   ,isd,jsd,ntau=tau,positive=.true.)
    enddo
    call user_get_2d_tracer_values(tracer_list,isd,ied,jsd,jed,nk) ! get all 2d tracer values
    !
    !
    !Get the pointers to the prognostic tracers after applying the physics.
    !
    call g_tracer_get_pointer(tracer_list,'no3'     ,'field',ergom%p_no3     )
    call g_tracer_get_pointer(tracer_list,'nh4'     ,'field',ergom%p_nh4     )
    call g_tracer_get_pointer(tracer_list,'po4'     ,'field',ergom%p_po4     )
    call g_tracer_get_pointer(tracer_list,'o2'      ,'field',ergom%p_o2      )
    call g_tracer_get_pointer(tracer_list,'h2s'     ,'field',ergom%p_h2s     )
    call g_tracer_get_pointer(tracer_list,'sul'     ,'field',ergom%p_sul     )
    do n=1, NUM_PHYTO
      call g_tracer_get_pointer(tracer_list,trim(phyto(n)%name) ,'field',phyto(n)%p_phyt     )
    enddo
    do n=1, NUM_ZOO
      if (zoo(n)%vertical_migration)  &
      call g_tracer_get_pointer(tracer_list,trim(zoo(n)%name)   ,'vmove',zoo(n)%p_vmove  )
      if (zoo(n)%vertical_migration)  &
      call g_tracer_get_pointer(tracer_list,trim(zoo(n)%name)   ,'vdiff',zoo(n)%p_vdiff  )
      call g_tracer_get_pointer(tracer_list,trim(zoo(n)%name)   ,'field',zoo(n)%p_zoo  )
    enddo
    do n=1, NUM_SPM
      call g_tracer_get_pointer(tracer_list,trim(spm(n)%name)   ,'field',spm(n)%p_wat     )
    enddo
    call g_tracer_get_pointer(tracer_list,'nitrogen','field',ergom%p_nitrogen)

!   Use wrk1 to store the switch between oxic and anoxic 
!   may be replced by a smoother function later

    wrk1 = 1.0
    where (ergom%f_o2 <= 5.e-6) wrk1 = 0.

! The first index is the light band, currently there is only one.
! This is about the ERGOM algorithm. Take only half of the radiation.
! The radiation part can be improved, e.g. including solar angle.

    nb = 1
    k = 1
    do j = jsc, jec ; do i = isc, iec  !{
      ergom%irr_inst(i,j,k) = p5 * sw_pen_band(nb,i,j) * exp(-opacity_band(nb,i,j,k)* dzt(i,j,k)*p5)       
    enddo; enddo  !} i,j,k
    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec  !{
      ergom%irr_inst(i,j,k) = ergom%irr_inst(i,j,k-1) * exp(-opacity_band(nb,i,j,k)* dzt(i,j,k))       
    enddo; enddo ; enddo  !} i,j,k

! Diatoms and Flagellates

    do n=1, NUM_PHYTO; do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!some temp-fields
      ttemp = temp(i,j,k)
      stemp = salt(i,j,k)
!Light
      itemp   = ergom%irr_inst(i,j,k)
! Iopt of photosynthesis defined as the maximum of Imin and Itemp/2 (??)      
! photosynthesis vs. light response curve (Steel 1962, Limnol Oceanogr), concept of Iopt      
      Iopt    = max(itemp/2.0, phyto(n)%imin)                           
      ilim    = itemp/Iopt * exp(1 - itemp/Iopt)                        
      if (phyto(n)%id_ilim .gt. 0) phyto(n)%ilim(i,j,k) = ilim
!DIN
      temp_alpha  = phyto(n)%alpha
      temp_talpha = phyto(n)%talpha
      if (n .ne. CYA) then
! available DIN defined as no3+nh4
         ntemp   = ergom%f_no3(i,j,k) + ergom%f_nh4(i,j,k)                 
         ntemp_r = 1./max(ntemp,epsln)
         ntemp   = ntemp*ntemp
! nutrient uptake response (Fennel & Neumann 2001, AMBIO)         
         nlim    = ntemp / (temp_alpha*temp_alpha+ntemp)
!increase of growth with temperature following Neumann
         tlim    = 1. + ttemp*ttemp/(temp_talpha*temp_talpha + ttemp*ttemp)
      else
! Minimum temperature and minimum and maximum salinity
! salinity limitation for diazotrophs not rlevant for GENUS area, but e.g. for the Baltic Sea
         tlim    = 1. / (1. + exp(phyto(n)%tmin - ttemp)) &
                      / (1. + exp(stemp - phyto(n)%smax)) &          
                      / (1. + exp(phyto(n)%smin - stemp))
      endif
!DIP
      ptemp      = ergom%f_po4(i,j,k)
      ptemp      = ptemp*ptemp
! calculation of po4-alpha via Redfield ratio
      temp_alpha = temp_alpha*phyto(n)%pnr                       
! nutrient uptake response (Fennel & Neumann 2001, AMBIO)      
      plim    = ptemp / (temp_alpha*temp_alpha + ptemp)          

            
      !What factor limits growth?      
      if (n .eq. CYA) then
        rp      = min(plim,ilim) * tlim * phyto(n)%rp0    ! no DIN limitation
        phyto(n)%jprod_n2(i,j,k) = rp * (phyto(n)%f_n(i,j,k) + phyto(n)%p0)  
        phyto(n)%jprod_po4(i,j,k) = phyto(n)%jprod_n2(i,j,k)*phyto(n)%pnr
      else
        ! Liebigs law of minimum: DIN, DIP, light limitation
        rp      = min(nlim,plim,ilim) * tlim * phyto(n)%rp0                      
        jtemp   = rp * (phyto(n)%f_n(i,j,k) + phyto(n)%p0) * ntemp_r
        phyto(n)%jprod_no3(i,j,k) = jtemp * ergom%f_no3(i,j,k)
        phyto(n)%jprod_nh4(i,j,k) = jtemp * ergom%f_nh4(i,j,k)
        phyto(n)%jprod_po4(i,j,k) = (phyto(n)%jprod_no3(i,j,k) &
	                            +phyto(n)%jprod_nh4(i,j,k))*phyto(n)%pnr
      endif
    enddo; enddo ; enddo; enddo   !} i,j,k,n

! Diatoms, flagellates, diazotrophs
    det(SLOW)%jmort = 0.
    det(FAST)%jmort = 0.
    ! Phytoplankon loss (remember! wrk1=0 for anoxic conditions, wrk1=1 for oxic conditions)
    do n=1, NUM_PHYTO; do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  
      ! respiration of phytoplankton only at oxic conditions (wrk1)
      phyto(n)%jres_n(i,j,k) = phyto(n)%lpr * wrk1(i,j,k)*phyto(n)%f_n(i,j,k)
      ! mortality of phytoplankton (loss to detritus) 10x higher at anoxic than oxic conditions   
      phyto(n)%jdet_n(i,j,k) = phyto(n)%lpd * (10.-9.*wrk1(i,j,k))* phyto(n)%f_n(i,j,k)    
    enddo; enddo ; enddo; enddo   !} i,j,k,n
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      det(FAST)%jmort(i,j,k) = det(FAST)%jmort(i,j,k) + 0.5 * phyto(DIA)%jdet_n(i,j,k) 
      det(SLOW)%jmort(i,j,k) = det(SLOW)%jmort(i,j,k) + 0.5 * phyto(DIA)%jdet_n(i,j,k) 
      det(SLOW)%jmort(i,j,k) = det(SLOW)%jmort(i,j,k) + phyto(FLA)%jdet_n(i,j,k) 
      det(SLOW)%jmort(i,j,k) = det(SLOW)%jmort(i,j,k) + phyto(CYA)%jdet_n(i,j,k) 
    enddo ; enddo; enddo   !} i,j,k

! Find the chlorophyll concentration for sea water optics use µg/kg
! This needs improvement for example after Geider et al (1997) to be done
! Add detritus. In the surface it carries chlorophyll - in deep layers the error does not matter
    ergom%f_chl = 0.0
    do n=1, NUM_PHYTO
      ergom%f_chl(:,:,:) = ergom%f_chl(:,:,:) + phyto(n)%cnr * 12.0e6 * 0.025 * phyto(n)%f_n(:,:,:)
    enddo
    do n=1, NUM_DET ! There is chlorophyll in detritus
      ergom%f_chl(:,:,:) = ergom%f_chl(:,:,:) + phyto(1)%cnr * 12.0e6 * 0.025 * det(n)%f_n(:,:,:)    
    enddo

! Zooplankton

    do m=1,NUM_PHYTO
      phyto(m)%jgraz_n = 0. 
    enddo
    do m=1,NUM_DET
      det (m)%jgraz_n = 0.    
    enddo
    do m=1,NUM_ZOO
      zoo(m)%jgraz_n = 0.    
      zoo(m)%jgain_n = 0.    
    enddo

    do n=1, NUM_ZOO 
      ! Store available food in wrk2
      wrk2 = 0. ; wrk3=0.
      
      selectcase(zoo(n)%graz_pref)
      case(ERGOM_PREFS)
      case(HUTSON)
      case(HUTSON_QUADRATIC)
      case(GENUS_IVLEV) 
      case(GENUS_IVLEV_QUADRATIC)
      case DEFAULT
         call  mpp_error(FATAL, &
        '==>Error from ERGOM, undefined grazing preference.')
      end select
      ! grazing on phytoplankton
      do m=1, NUM_PHYTO
        wrk2(:,:,:) = wrk2(:,:,:) + zoo(n)%pref_phy(m)*phyto(m)%f_n(:,:,:)   
        wrk3(:,:,:) = wrk3(:,:,:) + phyto(m)%f_n(:,:,:)                      
      enddo
      ! grazing on detritus
      do m=1, NUM_DET
        wrk2(:,:,:) = wrk2(:,:,:) + zoo(n)%pref_det(m)*det(m)%f_n(:,:,:)     
        wrk3(:,:,:) = wrk3(:,:,:) + det(m)%f_n(:,:,:)                        
      enddo
      ! grazing on zooplankton
      do m=1, NUM_ZOO
        wrk2(:,:,:) = wrk2(:,:,:) + zoo(n)%pref_zoo(m)*zoo(m)%f_n(:,:,:)     
        wrk3(:,:,:) = wrk3(:,:,:) + zoo(m)%f_n(:,:,:)                        
      enddo
      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
         wrk2(i,j,k) = max(wrk2(i,j,k),epsln)
         wrk3(i,j,k) = max(wrk3(i,j,k),epsln)
      enddo ; enddo; enddo   !} i,j,k

      ! Vertical movement of zooplankton
      wrk4 = 1.
      if(zoo(n)%vertical_migration) then
        call generic_ERGOM_find_vmove(zoo(n), temp, wrk2, ergom%f_o2, ergom%f_h2s, wrk4, dzt, isd, ied, jsd, jed)
        ! wrk4 is used to reduce grazing when sinking fast
        ! The movement and diffusion are carried out in the generic triagonal solver IOWtridiag in generic_tracer_utils
        ! Otherwise use
	! call generic_ERGOM_vmove(zoo(n)%p_vmove, dzt, zoo(n)%p_zoo(:,:,:,tau), dt, isd, ied, jsd, jed)
      endif

      ! temporarily store temperature dependence factor for zooplankton grazing in zoo%jdet_n
      if (zoo(n)%blanchard_temperature) then
        ! Zooplankton temperature dependence after Blanchard 1996
        ! P_max(T) = P_MAX * ((T_max-T)/(T_max-T_opt))^beta * exp(-beta*((T_max-T)/(T_max-T_opt)-1))
        do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
          toptz = max(0.0,(zoo(n)%t_max-temp(i,j,k))/(zoo(n)%t_max-zoo(n)%t_opt))
          zoo(n)%jdet_n(i,j,k) = toptz ** zoo(n)%beta * exp(-zoo(n)%beta*(toptz-1))
        enddo ; enddo; enddo   !} i,j,k
      else
        ! Zooplankton temperature dependence as used in earlier ERGOM versions (note: range 1.0..2.0)
        ! P_max(T) = P_MAX * ( 1 + (T/T_opt * exp(1-T/T_opt))^2 )
        do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!          toptz = zoo(n)%t_opt
!          temp1 = 2./toptz
!          temp2 = 2.7183/(toptz*toptz)
!          tempq = max(temp(i,j,k),0.)
!          tempq = tempq*tempq
!          zoo(n)%jdet_n(i,j,k) = (1.0 + temp2*tempq*exp(1.0 - temp(i,j,k)*temp1))
          tempq                = max(temp(i,j,k),0.) / zoo(n)%t_opt
          tempq                = tempq * exp( 1 - tempq )
          zoo(n)%jdet_n(i,j,k) = 1 + tempq*tempq 
        enddo ; enddo; enddo   !} i,j,k
      endif

      tominq = zoo(n)%oxy_min**2
      tosubq = zoo(n)%oxy_sub**2
      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  
 	food_pref  = wrk2(i,j,k)
	food       = wrk3(i,j,k)
        temp3 = zoo(n)%iv*food
	teo2q = ergom%f_o2(i,j,k)**2
	! no feeding if anoxic
        temp2 = teo2q / (tominq + teo2q) 
	! Ivlev^2-grazing model, grazing only at oxic conditions               
        gg    = wrk4(i,j,k) * zoo(n)%graz*(1. - exp(-temp3*temp3)) * temp2     
        
	! temperature-, oxygen- and food-dependent total grazing rate [1/s]
        ! (multiply by temporarily stored temperature dependence in zoo(n)%jdet_n)
        graz_z = gg * zoo(n)%jdet_n(i,j,k)  
        
	! mortality enhancement factor (10x if anoxic)
        ! mortality of zooplankton (loss to detritus), 10x higher at anoxic than oxic conditions
        ! loss to detritus by feeding, 10x higher at anoxic than oxic conditions
        temp1 = 1.0 + 9.0 * tominq / (tominq + teo2q)                                            
        lzd1  =  zoo(n)%sigma_b * temp1                                                          
        lzd2  = (zoo(n)%food_to_det_2*gg + zoo(n)%food_to_det*graz_z) * temp1                    
        lzd   = lzd1 + lzd2

        !respiration reduction factor (reduced if suboxic, none if anoxic)
        temp2 = zoo(n)%resp_red * temp2 + (1.0-zoo(n)%resp_red) * teo2q / (tosubq + teo2q)       
        ! respiration of zooplankton        
        lzn   = (zoo(n)%nue     + zoo(n)%food_to_nh4_2*gg  + zoo(n)%food_to_nh4*graz_z) * temp2  

        ! grazing is divided by total food to distribute it to different food types        
        graz_z = graz_z*(zoo(n)%f_n(i,j,k)+zoo(n)%z0)/(food + epsln)                             
        
        zoo(n)%jgain_n(i,j,k)= graz_z * food_pref
        do m=1,NUM_PHYTO
          phyto(m)%jgraz_n(i,j,k) = phyto(m)%jgraz_n(i,j,k) + graz_z * zoo(n)%pref_phy(m)*phyto(m)%f_n(i,j,k)
        enddo
        do m=1,NUM_DET
          det(m)%jgraz_n(i,j,k)   = det(m)%jgraz_n(i,j,k)   + graz_z * zoo(n)%pref_det(m)* det(m)%f_n(i,j,k)
        enddo
        do m=1,NUM_ZOO
          zoo(m)%jgraz_n(i,j,k)   = zoo(m)%jgraz_n(i,j,k)   + graz_z * zoo(n)%pref_zoo(m)* zoo(m)%f_n(i,j,k)
        enddo
        
        zz0 = zoo(n)%f_n(i,j,k)           
        zz0 = zz0*zz0*zoo(n)%zcl1        ! zcl1: closure parameter for zooplankton [kg2/mol2]
        zoo(n)%jdet_n(i,j,k) = lzd * zz0 
        zoo(n)%jres_n(i,j,k) = lzn * zz0
        
        det(FAST)%jmort(i,j,k) = det(FAST)%jmort(i,j,k) + lzd1 * zz0
        det(SLOW)%jmort(i,j,k) = det(SLOW)%jmort(i,j,k) + lzd2 * zz0
        
      enddo ; enddo; enddo   !} i,j,k,n
    enddo   !} n


!   Nitrification

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      temp1 = ergom%f_o2(i,j,k)
      q10   = ergom%q10_nit * temp(i,j,k)
      ! Michaelis-Menten-type reaction equation
      ergom%jnitrif(i,j,k) = ergom%nf*temp1/(temp1+ergom%alpha_nit) * exp(q10)*ergom%f_nh4(i,j,k)       
    enddo; enddo ; enddo  !} i,j,k
!
    !
    !-----------------------------------------------------------------------
    !
    !     CALCULATE SOURCE/SINK TERMS FOR EACH TRACER
    !
    !-----------------------------------------------------------------------
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       ergom%jno3(i,j,k) =  ergom%jnitrif(i,j,k)             &  ! nitrification releases no3
                           - phyto(DIA)%jprod_no3(i,j,k)     &  ! no3-uptake by diatoms and flagellates
			   - phyto(FLA)%jprod_no3(i,j,k)     
       ergom%jnh4(i,j,k) = - ergom%jnitrif(i,j,k)            &  ! nitrification consumes nh4
                           - phyto(DIA)%jprod_nh4(i,j,k)     &  ! nh4-uptake by diatoms and flagellates
			   - phyto(FLA)%jprod_nh4(i,j,k)     &
			   + phyto(DIA)%jres_n(i,j,k)        &  ! nh4-release by all phytoplankton and zooplankton
			   + phyto(FLA)%jres_n(i,j,k)        &
			   + phyto(CYA)%jres_n(i,j,k) 
			     ! O2-release by no3 assimilation (photosynthesis) of diatoms and flagellates       
       ergom%jo2(i,j,k)  =   8.625*(phyto(DIA)%jprod_no3(i,j,k) + phyto(FLA)%jprod_no3(i,j,k)) &        
                             ! O2-release by nh4 assimilation (photosynthesis) of diatoms and flagellates
			   + 6.625*(phyto(DIA)%jprod_nh4(i,j,k) + phyto(FLA)%jprod_nh4(i,j,k)  &        
                             ! O2-release by n2-fixation and subsequent nh4 assimilation (photosynthesis)of cyanobacteria
			          + phyto(CYA)%jprod_n2(i,j,k))                                &        
                             ! O2-consumption by respiration of all phytoplankton and zooplankton
	                   - 6.625*(phyto(DIA)%jres_n(i,j,k) + phyto(FLA)%jres_n(i,j,k)        &        
			          + phyto(CYA)%jres_n(i,j,k))                                  &
			     ! nitrification consumes 2 mol o2              
		           - 2.   * ergom%jnitrif(i,j,k)                                               
    enddo; enddo ; enddo  !} i,j,k

    ergom%jpo4 = 0.
    do n = 1, NUM_PHYTO;do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       ergom%jpo4(i,j,k) = ergom%jpo4(i,j,k)                               &             
			   - phyto(n)%jprod_po4(i,j,k)                     &
			   + phyto(n)%jres_n(i,j,k) * phyto(n)%pnr
    enddo; enddo; enddo ; enddo  !} i,j,k,n

    do n = 1, NUM_ZOO;do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       ergom%jnh4(i,j,k) =  ergom%jnh4(i,j,k) + zoo(n)%jres_n(i,j,k) 
       ergom%jo2 (i,j,k) =  ergom%jo2(i,j,k)  - 6.625 * zoo(n)%jres_n(i,j,k)                    
       ergom%jpo4(i,j,k) =  ergom%jpo4(i,j,k) + zoo(n)%jres_n(i,j,k) * zoo(n)%pnr         
    enddo; enddo; enddo ; enddo  !} i,j,k,n

    !-----------------------------------------------------------------------
    !
    !Update the prognostics tracer fields via their pointers.
    !
    !-----------------------------------------------------------------------
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       ergom%p_no3(i,j,k,tau) = ergom%p_no3(i,j,k,tau) + ergom%jno3(i,j,k) * dt * grid_tmask(i,j,k)
       ergom%p_nh4(i,j,k,tau) = ergom%p_nh4(i,j,k,tau) + ergom%jnh4(i,j,k) * dt * grid_tmask(i,j,k)
       ergom%p_po4(i,j,k,tau) = ergom%p_po4(i,j,k,tau) + ergom%jpo4(i,j,k) * dt * grid_tmask(i,j,k)
       phyto(DIA)%p_phyt(i,j,k,tau) = phyto(DIA)%p_phyt(i,j,k,tau) + (phyto(DIA)%jprod_no3(i,j,k)  &
                              + phyto(DIA)%jprod_nh4(i,j,k)                            &
                              - phyto(DIA)%jres_n(i,j,k) - phyto(DIA)%jdet_n(i,j,k)    &
			      - phyto(DIA)%jgraz_n(i,j,k)) * dt * grid_tmask(i,j,k)
       phyto(FLA)%p_phyt(i,j,k,tau) = phyto(FLA)%p_phyt(i,j,k,tau) + (phyto(FLA)%jprod_no3(i,j,k)  &
                              + phyto(FLA)%jprod_nh4(i,j,k)                            &
			      - phyto(FLA)%jres_n(i,j,k) - phyto(FLA)%jdet_n(i,j,k)    & 
			      - phyto(FLA)%jgraz_n(i,j,k)) * dt * grid_tmask(i,j,k)
       phyto(CYA)%p_phyt(i,j,k,tau) = phyto(CYA)%p_phyt(i,j,k,tau) + (phyto(CYA)%jprod_n2(i,j,k)   &
                              - phyto(CYA)%jres_n(i,j,k) - phyto(CYA)%jdet_n(i,j,k)    &  
			      - phyto(CYA)%jgraz_n(i,j,k)) * dt * grid_tmask(i,j,k)
       ergom%p_nitrogen(i,j,k,tau) = ergom%p_nitrogen(i,j,k,tau)                       &
                              - 0.5*phyto(CYA)%jprod_n2(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    do n = 1, NUM_DET
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       spm(det(n)%index_spm)%p_wat(i,j,k,tau) = spm(det(n)%index_spm)%p_wat(i,j,k,tau) + (&
       			         det(n)%jmort(i,j,k) - det(n)%jgraz_n(i,j,k)          &  
                                ) * dt * grid_tmask(i,j,k)
       enddo; enddo ; enddo  !} n,i,j,k
    enddo

    do n = 1, NUM_ZOO
       do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
          zoo(n)%p_zoo(i,j,k,tau) = zoo(n)%p_zoo(i,j,k,tau) + (zoo(n)%jgain_n(i,j,k)        & 
                                  - zoo(n)%jdet_n(i,j,k) - zoo(n)%jres_n(i,j,k)            & 
                                  - zoo(n)%jgraz_n(i,j,k) ) * dt * grid_tmask(i,j,k)
       enddo; enddo ; enddo  !} n,i,j,k
    enddo
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
       ergom%p_o2(i,j,k,tau) = ergom%p_o2(i,j,k,tau) + ergom%jo2(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

! Sulfur cycle elements

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      o2_avail  = max(ergom%p_o2 (i,j,k,tau), 0.0)                      ! semi-implicit approach
      no3_avail = max(ergom%p_no3(i,j,k,tau), 0.0) 
      h2s_avail = max(ergom%p_h2s(i,j,k,tau), 0.0) 
      sul_avail = max(ergom%p_sul(i,j,k,tau), 0.0) 
      q10       = exp(ergom%q10_h2s*temp(i,j,k))


      ! inhibition of nitrate reduction in the presence of oxygen      
      temp2    = 5.e-13/(5.e-13 + o2_avail*o2_avail)                    
      
      ! the following may be improved by MM shaped functions                                                                 
      no3_h2s  = q10 * no3_avail * ergom%k_h2s_no3 *temp2 * h2s_avail   ! mainly reacts at anoxic conditions
      no3_sul  = q10 * no3_avail * ergom%k_sul_no3 *temp2 * sul_avail   ! mainly reacts at anoxic conditions
      o2_h2s   = q10 * o2_avail  * ergom%k_h2s_o2	  * h2s_avail   ! mainly reacts at oxic conditions
      o2_sul   = q10 * o2_avail  * ergom%k_sul_o2	  * sul_avail   ! mainly reacts at oxic conditions

      ! make scheme positive, clip with available concentrations
      ! and rescale the reaction rates
      ! check availability of the reaction partners       
      dh2s    = dt * (o2_h2s + no3_h2s)	 
      ! scaling with sca                               
      sca     = min(dh2s, h2s_avail)/(dh2s + epsln)                     
      o2_h2s  = o2_h2s  * sca
      no3_h2s = no3_h2s * sca

      do2     = dt * ( 0.5 * o2_h2s + 1.5 * o2_sul)
      sca     = min(do2, o2_avail)/(do2 + epsln)
      o2_h2s  = o2_h2s  * sca
      o2_sul  = o2_sul  * sca
      
      dno3    = dt * ( 0.4 * no3_h2s + 1.2 * no3_sul)
      sca     = min(dno3, no3_avail)/(dno3 + epsln)
      no3_h2s = no3_h2s * sca
      no3_sul = no3_sul * sca

      dsul    = dt * (o2_sul + no3_sul)		
      sca     = min(dsul, sul_avail+dt*(o2_h2s+o2_h2s))/(dsul + epsln)
      o2_sul  = o2_sul  * sca
      no3_sul = no3_sul * sca

      ergom%p_no3(i,j,k,tau) = ergom%p_no3(i,j,k,tau) + dt * &
         ( - 0.4 * no3_h2s - 1.2 * no3_sul )           

      ergom%p_o2(i,j,k,tau)  = ergom%p_o2(i,j,k,tau)  + dt * &
         ( - 0.5 * o2_h2s  - 1.5 * o2_sul )                  

      ergom%p_h2s(i,j,k,tau) = ergom%p_h2s(i,j,k,tau) + dt * &
         ( - o2_h2s - no3_h2s )       

      ergom%p_sul(i,j,k,tau) = ergom%p_sul(i,j,k,tau) + dt * &
         ( o2_h2s - o2_sul + (no3_h2s - no3_sul))
      
      if(ergom%id_jh2s_o2 .gt. 0) ergom%jh2s_o2(i,j,k) = - o2_h2s 
      if(ergom%id_jh2s_no3.gt. 0) ergom%jh2s_no3(i,j,k)= - no3_h2s
      if(ergom%id_jsul_o2 .gt. 0) ergom%jsul_o2(i,j,k) = - o2_sul 
      if(ergom%id_jsul_no3.gt. 0) ergom%jsul_no3(i,j,k)= - no3_sul
      
      ergom%jno3(i,j,k) = ergom%jno3(i,j,k) - 0.4 * no3_h2s - 1.2 * no3_sul
      ergom%jo2(i,j,k)  = ergom%jo2(i,j,k)  - 0.5 * o2_h2s  - 1.5 * o2_sul
      ergom%jh2s(i,j,k) = ergom%jh2s(i,j,k) - o2_h2s - no3_h2s
      ergom%jsul(i,j,k) = ergom%jsul(i,j,k) + o2_h2s - o2_sul + no3_h2s - no3_sul

    enddo; enddo ; enddo  !} i,j,k

! Vertical movement of phytoplankton and suspended particulate matter

    phyto(dia)%move = - phyto(dia)%wsink0
    call generic_ERGOM_vmove(phyto(dia)%move, dzt, phyto(DIA)%p_phyt(:,:,:,tau), dt, isd, ied, jsd, jed)
    phyto(cya)%move = - phyto(cya)%wsink0
    call generic_ERGOM_vmove(phyto(cya)%move, dzt, phyto(CYA)%p_phyt(:,:,:,tau), dt, isd, ied, jsd, jed)
    do n = 1, NUM_SPM
      spm(n)%move     = - spm(n)%wsink0
      call generic_ERGOM_vmove(spm(n)%move, dzt, spm(n)%p_wat(:,:,:,tau), dt, isd, ied, jsd, jed)
    enddo

! Oxic and anoxic detritus recycling (denitrification, anammox)
 
    wrk1 = epsln
    wrk2 = 0.
    ! find a first guess for the recycling rate
    do n=1, NUM_DET; do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  
       q10  = det(n)%q10_rec * temp(i,j,k)              
       wrk1(i,j,k) = wrk1(i,j,k) + det(n)%f_n(i,j,k)    
       wrk2(i,j,k) = wrk2(i,j,k) + det(n)%dn * exp(q10)*det(n)%f_n(i,j,k)    
    enddo; enddo; enddo; enddo
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  
      ! temperature dependenc of anammox
      ! i.e. only between 7 and 30degC (Dalsgaard & Thamdrup 2002 in ApplEnvMicrobiol) 
      k_an = ergom%k_an0 * det(SLOW)%f_n(i,j,k) 	 &    ! only det(SLOW)
        	   / (1. + exp((temp(i,j,k) - 30.)*0.5)) &    
    		   / (1. + exp((7.  - temp(i,j,k))*0.5))      

      det_recyc = wrk2(i,j,k)				      ! first guess for the recycling rate
      o2_avail  = max(ergom%p_o2 (i,j,k,tau), 0.0)	      ! available amount of electron donors
      no3_avail = max(ergom%p_no3(i,j,k,tau), 0.0) 
      h2s_avail = max(ergom%p_h2s(i,j,k,tau), 0.0) 
      nh4_avail = max(ergom%p_nh4(i,j,k,tau), 0.0)	      ! nh4 for anammox
      
      ! functions tend to 1 for high concentrations and have a linear slope at X=0
      r_o2  = 2./(1.+exp(-2.*ergom%alp_o2  * o2_avail)) -1.   ! "soft switch"
      r_no3 = 2./(1.+exp(-2.*ergom%alp_no3 * no3_avail))-1.   
      r_h2s = 2./(1.+exp(-2.*ergom%alp_h2s * h2s_avail))-1.
      r_nh4 = 2./(1.+exp(-2.*ergom%alp_nh4 * nh4_avail))-1.
    									
      ! subdivide the recycling into the four possible electron donors
      ! oxygenrecycling = o (if oxygen is present)
      ! denitrification = (NOT o) AND n AND NOT (a AND (NOT h))
      ! sulfatreduction = (NOT o) AND (NOT n)
      ! anammox = (NOT o) AND n AND a AND (NOT h)									
      rec_o2  = det_recyc*r_o2  					
      rec_no3 = k_DN * det_recyc*(1-r_o2)*r_no3*(1-r_nh4*(1-r_h2s))	
      rec_so4 = k_DS * det_recyc*(1-r_o2)*(1-r_no3)			
      rec_ana = k_an * (1-r_o2)*r_no3*r_nh4*(1-r_h2s)			

      ! anammox is regarded as additional "bonus" recycling
      det_recyc = rec_o2+rec_no3+rec_so4+rec_ana			

      if(ergom%id_jrec_o2 .gt. 0) ergom%jrec_o2(i,j,k)  = rec_o2
      if(ergom%id_jrec_no3.gt. 0) ergom%jrec_no3(i,j,k) = rec_no3
      if(ergom%id_jrec_so4.gt. 0) ergom%jrec_so4(i,j,k) = rec_so4   
      if(ergom%id_jrec_ana.gt. 0) ergom%jrec_ana(i,j,k) = rec_ana 
      ergom%jdenit_wc(i,j,k) = ldn_N * rec_no3

      ergom%jno3(i,j,k)    = ergom%jno3(i,j,k) - ldn_A * rec_ana - ldn_N * rec_no3
      ergom%jnh4(i,j,k)    = ergom%jnh4(i,j,k) - ldn_A * rec_ana + det_recyc	   
      ergom%jo2(i,j,k)     = ergom%jo2(i,j,k)  - ldn_O * rec_o2 
      ergom%jh2s(i,j,k)    = ergom%jh2s(i,j,k) + ldn_S * rec_so4 
      ergom%jpo4(i,j,k)    = ergom%jpo4(i,j,k) + det_recyc * phyto(DIA)%pnr   

      ergom%p_o2(i,j,k,tau)  = ergom%p_o2 (i,j,k,tau) - dt * grid_tmask(i,j,k) * ldn_O * rec_o2
      ergom%p_nh4(i,j,k,tau) = ergom%p_nh4(i,j,k,tau) - dt * grid_tmask(i,j,k) * (ldn_A * rec_ana - det_recyc)
      ergom%p_no3(i,j,k,tau) = ergom%p_no3(i,j,k,tau) - dt * grid_tmask(i,j,k) * (ldn_A * rec_ana + ldn_N * rec_no3)
      ergom%p_h2s(i,j,k,tau) = ergom%p_h2s(i,j,k,tau) + dt * grid_tmask(i,j,k) * ldn_S * rec_so4
      ergom%p_po4(i,j,k,tau) = ergom%p_po4(i,j,k,tau) + dt * grid_tmask(i,j,k) * det_recyc * phyto(DIA)%pnr
      ergom%p_nitrogen(i,j,k,tau) = ergom%p_nitrogen(i,j,k,tau) &
    				             + dt * grid_tmask(i,j,k) * 0.5 * (2.*ldn_A * rec_ana + ldn_N * rec_no3)
      wrk2(i,j,k) = - dt * grid_tmask(i,j,k) * det_recyc
    enddo; enddo ; enddo  !} i,j,k
    do n=1, NUM_DET; do k = 1, nk ; do j = jsc, jec ; do i = isc, iec 
        spm(det(n)%index_spm)%p_wat(i,j,k,tau) = & 
           spm(det(n)%index_spm)%p_wat(i,j,k,tau) + wrk2(i,j,k) * det(n)%f_n(i,j,k)/wrk1(i,j,k)
    enddo; enddo; enddo; enddo

! Detritus recycling in the sediment by sulfur bacteria (e.g. Beggiatoa)

! Sulfur bacteria colonize anoxic sediment (no o2 or no3)
! i.e. organic matter is recycled with so4 reduction only and h2s is available.
! This h2s appears at the bottom of the lowest layer (mol/m2 -> mol/kg) 
! Assumption: The sediment is covered by bacteria, which prevent direct exchange with the water column.
    biosed%bioerosion=0.0
    if (biosed%id_jdenit_sed .gt. 0) biosed%jdenit_sed = 0.
    do j = jsc, jec ; do i = isc, iec  
      k=grid_kmt(i,j)
      if (k == 0) cycle
      rhodztdt_r = rho_dzt(i,j,k)/dt
      q10        = biosed%q10_rec * temp(i,j,k)
      ! assumption that only the upper layer (5 cm) is involved. The unit is unit = mol/m2   
      active_sed = max(sed(biosed%index_sed)%f_sed(i,j,1),0.)          
                                                                                                                                                                       
      det_recyc =  biosed%dn * exp(q10) * active_sed
            
      if (active_sed > biosed%thio_bact_min ) then   
         
	 ! the layer may support sulfur bacteria, this is anoxic sediment
	 det_recyc = biosed%frac_dn_anoxic*det_recyc                       
         
	 if (biosed%id_mode_sed .gt. 0) biosed%mode_sed  = 1.
	 dsed      =  dt*det_recyc
         dh2s_dt   = 0.5*ldn_O*det_recyc                                
         o2_avail  = max(ergom%p_o2(i,j,k,tau), 0.0)
         do2_dt    = min(o2_avail*rhodztdt_r, 2.*dh2s_dt)
	 
	 ! oxidation of h2s with o2 
         ergom%b_o2(i,j) = do2_dt                                       
         dh2s_dt   = dh2s_dt - 0.5 * do2_dt
         
	 ! oxidation of the remaining h2s with no3         
         no3_avail = max(ergom%p_no3(i,j,k,tau), 0.0)                   
         dno3_dt   = min(no3_avail*rhodztdt_r, dh2s_dt)
         dh2s_dt   = dh2s_dt - dno3_dt
	 
	 ! sulfur bacteria release nh4, i.e. no denitrification
	 ! release the remaining h2s into the water
	 ergom%b_no3(i,j) = dno3_dt
	 ergom%b_nh4(i,j) = - dno3_dt - det_recyc                       
         ergom%b_h2s(i,j) = - dh2s_dt                                   
         
	 ! po4 recycling
         ergom%b_po4(i,j) = - det_recyc*biosed%pnr                         

      else  ! sediment too thin for the growth of sulfur bacteria

         if (biosed%id_mode_sed .gt. 0) biosed%mode_sed  = 0.
         ! now unit back to mol/kg                                                                                                
         pot_o2 = ergom%p_o2(i,j,k,tau) - 2. * ergom%p_h2s(i,j,k,tau)
	 if (pot_o2 < 0.) det_recyc = biosed%frac_dn_anoxic*det_recyc
	 dsed      =  dt*det_recyc
                                                                
         ! use o2 for re-mineralisation
	 o2_avail  = max(ergom%p_o2(i,j,k,tau), 0.0)                    
         do2_dt    = min(o2_avail*rhodztdt_r, det_recyc*ldn_O)
         ergom%b_o2(i,j) = do2_dt
         det_left  =  det_recyc - do2_dt / ldn_O

         ! use no3 for re-mineralisation if o2-availability is too low
         no3_avail = max(ergom%p_no3(i,j,k,tau), 0.0)                   
         dno3_dt   = min(no3_avail*rhodztdt_r, det_left*ldn_N)
	 ergom%b_no3(i,j)      = dno3_dt
         ergom%b_nitrogen(i,j) = -0.5 * dno3_dt
         det_left = det_left - dno3_dt / ldn_N

         ! use so4 for re-mineralisation of the remaining part if o2 and no3-availability is too low
	 dh2s_dt   = 0.5*ldn_O*det_left                                 
         
         ergom%b_h2s(i,j) = - dh2s_dt    ! release of the remaining part into the water

         ! if the bottom water is oxic, denitrification happens in the sediment
         ! in this case less nh4 is released.
         ! oxygen used for nitrification in the sediment is 2.* biosed%den_rate * det_recyc	 
         temp2 =  dt * (2.*biosed%den_rate * det_recyc + do2_dt)                                   
         if (pot_o2 > temp2) then   !oxic sediment --> denitrification                                                            

	   ergom%b_o2(i,j) = ergom%b_o2(i,j) + 2.* biosed%den_rate * det_recyc                     
	   dsed = dsed + det_recyc*biosed%den_rate/ldn_N*dt
	   ergom%b_nh4(i,j) = - det_recyc *(1.+ biosed%den_rate*(1./ldn_N - 1.))
           ergom%b_nitrogen(i,j) = ergom%b_nitrogen(i,j) - det_recyc * 0.5* biosed%den_rate
           
	   ! po4 recycling
           ergom%b_po4(i,j) = - det_recyc*biosed%pnr*(1. + biosed%den_rate/ldn_N)                     
           if (biosed%id_jdenit_sed .gt. 0) biosed%jdenit_sed(i,j)  = det_recyc*biosed%den_rate
	 
	 else                        !unoxic sediment --> no denitrification  
	 
	   ergom%b_nh4(i,j) = - det_recyc 
           
	   ! po4 recycling
           ergom%b_po4(i,j) = - det_recyc*biosed%pnr                                               
	 endif 
      endif      
      if (biosed%id_jrec_n .gt. 0) biosed%jrec_n(i,j) = dsed/dt
      sed(biosed%index_sed)%f_sed(i,j,1)   = sed(biosed%index_sed)%f_sed(i,j,1) - dsed
    enddo; enddo  !} i,j

    if (biosed%po4_retention .gt. 0.0) then
      ! Apply retention of phosphate as iron phosphate in the sediment
      ! iron phosphate is stored in sed(biosed%index_ips)
      wrk1 = 1.0
      where (ergom%f_o2 <= 5.e-6) wrk1 = 0.
      do j = jsc, jec 
        temp1 = biosed%po4_retention
        ! to do: add biosed%po4_ret_plus_BB if y>60.75
        do i = isc, iec 
          k=grid_kmt(i,j)
          if (k == 0) cycle 
          ! retention of phosphate
          ! rate at which phosphorous is retained in the sediment [mol m-2 s-1], remember b_po4 is negative
          ! phosphate flux into the water column (b_po4, negative) is reduced
	  temp2 = -ergom%b_po4(i,j) * temp1 * wrk1(i,j,k)    
          ergom%b_po4(i,j) = ergom%b_po4(i,j)+temp2          
          sed(biosed%index_ips)%f_sed(i,j,1) = sed(biosed%index_ips)%f_sed(i,j,1) + temp2*dt
          if (sed(biosed%index_ips)%id_jgain_sed .gt. 0) &
	                                         sed(biosed%index_ips)%jgain_sed(i,j) = temp2
          ! liberation of phosphate under anoxic conditions
          ! rate at which phosphorous is liberated from sediment [mol m-2 s-1]
	  temp2 = sed(biosed%index_ips)%f_sed(i,j,1) * (1.0-wrk1(i,j,k)) * biosed%po4_lib_rate 
          ergom%b_po4(i,j) = ergom%b_po4(i,j) - temp2
          sed(biosed%index_ips)%f_sed(i,j,1) = sed(biosed%index_ips)%f_sed(i,j,1) - temp2 * dt
          if (sed(biosed%index_ips)%id_jloss_sed .gt. 0) & 
	                                         sed(biosed%index_ips)%jloss_sed(i,j) = temp2
        enddo !} i
      enddo !} j
    endif

    ! Do sedimentation and resuspension
    call sedimentation_and_resuspension(NUM_SPM, spm, NUM_SED, sed, &
    isc, iec, jsc, jec, isd, ied, jsd, jed, grid_kmt, dzt, rho_dzt, tau, dt, sed_defs, &
    current_wave_stress, biosed%bioerosion)

    call user_set_2d_tracer_values(tracer_list, model_time)  ! save all 2d tracers
    call g_tracer_set_values(tracer_list,'chl',     'field',ergom%f_chl, isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'nh4',      'btf', ergom%b_nh4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'no3',      'btf', ergom%b_no3 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2',       'btf', ergom%b_o2  ,isd,jsd)
    call g_tracer_set_values(tracer_list,'h2s',      'btf', ergom%b_h2s ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4',      'btf', ergom%b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'nitrogen', 'btf', ergom%b_nitrogen ,isd,jsd)
    
!   Diagnostics

    call g_tracer_get_pointer(tracer_list,'no3','runoff_tracer_flux',runoff_flux_no3)
    call g_tracer_get_pointer(tracer_list,'no3','drydep',dry_no3)
    call g_tracer_get_pointer(tracer_list,'no3','wetdep',wet_no3)
    call g_tracer_get_pointer(tracer_list,'nh4','runoff_tracer_flux',runoff_flux_nh4)
    call g_tracer_get_pointer(tracer_list,'nh4','drydep',dry_nh4)
    call g_tracer_get_pointer(tracer_list,'nh4','wetdep',wet_nh4)
    call g_tracer_get_pointer(tracer_list,'po4','runoff_tracer_flux',runoff_flux_po4)
    call g_tracer_get_pointer(tracer_list,'po4','drydep',dry_po4)
    call g_tracer_get_pointer(tracer_list,'po4','wetdep',wet_po4)
    do n=1,NUM_SED
       call g_tracer_set_values(tracer_list,trim(spm(n)%name),  'btf', spm(n)%btf ,isd,jsd)
    enddo !n


    if (ergom%id_runoff_flux_po4 .gt. 0)      &
       used = send_data(ergom%id_runoff_flux_po4, runoff_flux_po4,           &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_runoff_flux_nh4 .gt. 0)      &
       used = send_data(ergom%id_runoff_flux_nh4, runoff_flux_nh4,           &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_runoff_flux_no3 .gt. 0)      &
       used = send_data(ergom%id_runoff_flux_no3, runoff_flux_no3,           &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_dep_dry_po4 .gt. 0)          &
       used = send_data(ergom%id_dep_dry_po4, dry_po4,                       &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_dep_dry_nh4 .gt. 0)     &
       used = send_data(ergom%id_dep_dry_nh4, dry_nh4,                       &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_dep_dry_no3 .gt. 0)     &
       used = send_data(ergom%id_dep_dry_no3, dry_no3,                       &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_dep_wet_po4 .gt. 0)     &
       used = send_data(ergom%id_dep_wet_po4, wet_po4,                       &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_dep_wet_nh4 .gt. 0)     &
       used = send_data(ergom%id_dep_wet_nh4, wet_nh4,                       &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (ergom%id_dep_wet_no3 .gt. 0)     &
       used = send_data(ergom%id_dep_wet_no3, wet_no3,                       &
       model_time, rmask = grid_tmask(:,:,1), & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!   Diagnostics

    if (ergom%id_irr_inst .gt. 0)           &
         used = send_data(ergom%id_irr_inst,     ergom%irr_inst,             &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jnh4 .gt. 0)           &
         used = send_data(ergom%id_jnh4,         ergom%jnh4*rho_dzt,         &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jno3 .gt. 0)           &
         used = send_data(ergom%id_jno3,         ergom%jno3*rho_dzt,         &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jdenit_wc .gt. 0)           &
         used = send_data(ergom%id_jdenit_wc,    ergom%jdenit_wc*rho_dzt,    &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jo2 .gt. 0)           &
         used = send_data(ergom%id_jo2,          ergom%jo2*rho_dzt,          &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jpo4 .gt. 0)           &
         used = send_data(ergom%id_jpo4,         ergom%jpo4*rho_dzt,         &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jnitrif .gt. 0)           &
         used = send_data(ergom%id_jnitrif,      ergom%jnitrif*rho_dzt,      &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jh2s_o2 .gt. 0)           &
         used = send_data(ergom%id_jh2s_o2,      ergom%jh2s_o2*rho_dzt,      &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jh2s_no3 .gt. 0)           &
         used = send_data(ergom%id_jh2s_no3,     ergom%jh2s_no3*rho_dzt,     &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jsul_o2 .gt. 0)           &
         used = send_data(ergom%id_jsul_o2,      ergom%jsul_o2*rho_dzt,      &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jsul_no3 .gt. 0)           &
         used = send_data(ergom%id_jsul_no3,     ergom%jsul_no3*rho_dzt,     &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jrec_o2 .gt. 0)	   &
    	 used = send_data(ergom%id_jrec_o2,	 ergom%jrec_o2*rho_dzt,      &
    	 model_time, rmask = grid_tmask(:,:,:),& 
    	 is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jrec_no3 .gt. 0)	    &
    	 used = send_data(ergom%id_jrec_no3,	 ergom%jrec_no3*rho_dzt,     &
    	 model_time, rmask = grid_tmask(:,:,:),& 
    	 is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jrec_so4 .gt. 0)	    &
    	 used = send_data(ergom%id_jrec_so4,	 ergom%jrec_so4*rho_dzt,     &
    	 model_time, rmask = grid_tmask(:,:,:),& 
    	 is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (ergom%id_jrec_ana .gt. 0)	    &
    	 used = send_data(ergom%id_jrec_ana,	 ergom%jrec_ana*rho_dzt,     &
    	 model_time, rmask = grid_tmask(:,:,:),& 
    	 is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !
    do n=1, NUM_DET
      if (det(n)%id_jgraz_n .gt. 0)           &
           used = send_data(det(n)%id_jgraz_n,      det(n)%jgraz_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (det(n)%id_jmort .gt. 0)           &
           used = send_data(det(n)%id_jmort,         det(n)%jmort*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
   enddo
    !
    do n=1, NUM_PHYTO
      if (phyto(n)%id_jprod_no3 .gt. 0)           &
           used = send_data(phyto(n)%id_jprod_no3,      phyto(n)%jprod_no3*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (phyto(n)%id_jprod_nh4 .gt. 0)           &
           used = send_data(phyto(n)%id_jprod_nh4,      phyto(n)%jprod_nh4*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (phyto(n)%id_jprod_po4 .gt. 0)           &
           used = send_data(phyto(n)%id_jprod_po4,      phyto(n)%jprod_po4*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (phyto(n)%id_jgraz_n .gt. 0)           &
           used = send_data(phyto(n)%id_jgraz_n,      phyto(n)%jgraz_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (phyto(n)%id_jres_n .gt. 0)           &
           used = send_data(phyto(n)%id_jres_n,      phyto(n)%jres_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (phyto(n)%id_jdet_n .gt. 0)           &
           used = send_data(phyto(n)%id_jdet_n,      phyto(n)%jdet_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (phyto(n)%id_ilim .gt. 0)           &
           used = send_data(phyto(n)%id_ilim,      phyto(n)%ilim,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo
!
    do n=1, NUM_ZOO
      if (zoo(n)%id_jgraz_n .gt. 0)           &
           used = send_data(zoo(n)%id_jgraz_n,      zoo(n)%jgraz_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (zoo(n)%id_jgain_n .gt. 0)           &
           used = send_data(zoo(n)%id_jgain_n,      zoo(n)%jgain_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (zoo(n)%id_jres_n .gt. 0)           &
           used = send_data(zoo(n)%id_jres_n,      zoo(n)%jres_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (zoo(n)%id_jdet_n .gt. 0)           &
           used = send_data(zoo(n)%id_jdet_n,      zoo(n)%jdet_n*rho_dzt,              &
           model_time, rmask = grid_tmask(:,:,:),& 
           is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo
!
    if (biosed%id_jrec_n .gt. 0)           & 
           used = send_data(biosed%id_jrec_n,     biosed%jrec_n,      & 
           model_time, rmask = grid_tmask(:,:,1),                     & 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (biosed%id_jdenit_sed .gt. 0)           & 
           used = send_data(biosed%id_jdenit_sed,     biosed%jdenit_sed,      & 
           model_time, rmask = grid_tmask(:,:,1),                     & 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (biosed%id_mode_sed .gt. 0)           & 
           used = send_data(biosed%id_mode_sed,     biosed%mode_sed,      & 
           model_time, rmask = grid_tmask(:,:,1),                     & 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    do n=1, NUM_SPM
      if (spm(n)%id_jsed .gt. 0)           &
           used = send_data(spm(n)%id_jsed, spm(n)%jsed, &
           model_time, rmask = grid_tmask(:,:,1),& 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    enddo
    
    do n=1, NUM_SED
      if (sed(n)%id_jgain_sed .gt. 0)           &
           used = send_data(sed(n)%id_jgain_sed, sed(n)%jgain_sed, &
           model_time, rmask = grid_tmask(:,:,1),& 
           is_in=isc, js_in=jsc ,ie_in=iec, je_in=jec)
      if (sed(n)%id_jloss_sed .gt. 0)           &
           used = send_data(sed(n)%id_jloss_sed, sed(n)%jloss_sed, &
           model_time, rmask = grid_tmask(:,:,1),& 
           is_in=isc, js_in=jsc ,ie_in=iec, je_in=jec)
      if (sed(n)%id_jres .gt. 0)           &
           used = send_data(sed(n)%id_jres, sed(n)%jres, &
           model_time, rmask = grid_tmask(:,:,1),& 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (sed(n)%id_jbiores .gt. 0)           &
           used = send_data(sed(n)%id_jbiores, sed(n)%jbiores, &
           model_time, rmask = grid_tmask(:,:,1),& 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    enddo

    deallocate(wrk1)
    deallocate(wrk2)
    deallocate(wrk3)
    deallocate(wrk4)
    
    call mpp_clock_end(id_source)

    return
  end subroutine generic_ERGOM_update_from_source

  ! <SUBROUTINE NAME="generic_ERGOM_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_ERGOM_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: SST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac, sal, tk, ta, tr, ST, alpha_nit
    real    :: ts, ts2, ts3, ts4, ts5, o2_saturation
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: nitrogen_field, o2_field
    real, dimension(:,:), ALLOCATABLE :: nitrogen_alpha, nitrogen_csurf, nitrogen_sc_no
    real, dimension(:,:), ALLOCATABLE :: o2_alpha,       o2_csurf      , o2_sc_no

    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_set_boundary_values'


    !nnz: Can we treat these as source and move block to user_update_from_source?
    !
    !=============
    !Block Starts: Calculate the boundary values
    !=============
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    call g_tracer_get_pointer(tracer_list,'nitrogen','field',nitrogen_field)
    call g_tracer_get_pointer(tracer_list,'o2'      ,'field',o2_field)

    allocate(nitrogen_alpha(isd:ied, jsd:jed)); nitrogen_alpha=0.0
    allocate(nitrogen_csurf(isd:ied, jsd:jed)); nitrogen_csurf=0.0
    allocate(nitrogen_sc_no(isd:ied, jsd:jed)); nitrogen_sc_no=0.0
    allocate(o2_alpha(isd:ied, jsd:jed));       o2_alpha=0.0
    allocate(o2_csurf(isd:ied, jsd:jed));       o2_csurf=0.0
    allocate(o2_sc_no(isd:ied, jsd:jed));       o2_sc_no=0.0

    !The atmospheric code needs soluabilities in units of mol/m**3/atm
    !
    !MOM
    !       The factor 1.0e-06 Rho_0 is for the conversion  
    !       from µmol/(kg * atm) to mol/(m**3 * atm) 
    conv_fac = 1.0e-06 *  ergom%Rho_0 
    !
    !GOLD
    !conv_fac = 1.0e-09

    do j=jsc,jec ; do i=isc,iec

       !This calculation needs an input of SST and SSS
       sal = SSS(i,j)
       ST  = SST(i,j) 
       tk = ST + 273.15 
       ta = ergom%e1_nit - ST
       tr = LOG(ta/tk)

       nitrogen_sc_no(i,j)    = ergom%a1_nit  + ST * (ergom%a2_nit  + ST * (ergom%a3_nit  + ST * ergom%a4_nit )) * &
            grid_tmask(i,j,1)
       !---------------------------------------------------------------------
       !     Calculate solubility
       !     Hamme, R.C. and Emerson, S.R., 2004. The solubility of neon, nitrogen and argon in distilled water and seawater.
       !     Deep Sea Research Part I: Oceanographic Research Papers, 51(11): 1517-1528.
       !---------------------------------------------------------------------

       !nnz: MOM hmask=grid_tmask(i,j,1), GOLD hmask=G%hmask 
       ! units are µmol/kg, Hence multiply with 10**6 *rho0 to get mol/m**3
       alpha_nit = EXP(((ergom%t4_nit*tr + ergom%t3_nit-ergom%s3_nit*sal)*tr & 
                         + ergom%t2_nit-ergom%s2_nit*sal)*tr + ergom%t1_nit-ergom%s1_nit*sal)

       nitrogen_alpha(i,j) = alpha_nit * conv_fac

       nitrogen_csurf(i,j) = nitrogen_field(i,j,1,taum1)  * ergom%Rho_0

    enddo; enddo
    do j=jsc,jec ; do i=isc,iec
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
       sal = SSS(i,j)
       ST  = SST(i,j) 
       ta = 298.15 - ST
       tk = 273.15 + ST
       ts = log(ta / tk)
       ts2 = ts  * ts
       ts3 = ts2 * ts
       ts4 = ts3 * ts
       ts5 = ts4 * ts

       o2_saturation = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
            exp( ergom%t0_o2 + ergom%t1_o2*ts + ergom%t2_o2*ts2 + &
	         ergom%t3_o2*ts3 + ergom%t4_o2*ts4 + ergom%t5_o2*ts5 + &
            (ergom%b0_o2 + ergom%b1_o2*ts + ergom%b2_o2*ts2 + ergom%b3_o2*ts3 + ergom%c0_o2*sal)*sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the 
       !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
       !  Cycles, 12, 141-163).
       !---------------------------------------------------------------------

       o2_sc_no(i,j)  = ergom%a1_o2  + ST * (ergom%a2_o2  + ST * (ergom%a3_o2  + ST * ergom%a4_o2 )) * &
                  grid_tmask(i,j,1)
       o2_alpha(i,j) = o2_saturation 
       o2_csurf(i,j) = o2_field(i,j,1,taum1) * ergom%Rho_0 !nnz: MOM has rho(i,j,1,tau)
    enddo; enddo

    !=============
    !Block Ends: Calculate the boundary values
    !=============

    !
    !Set %csurf and %alpha for these tracers. This will mark them for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'nitrogen','alpha',nitrogen_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'nitrogen','csurf',nitrogen_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'nitrogen','sc_no',nitrogen_sc_no,isd,jsd)

    call g_tracer_set_values(tracer_list,'o2','alpha',o2_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2','csurf',o2_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2','sc_no',o2_sc_no,isd,jsd)

    deallocate(nitrogen_alpha,nitrogen_csurf,nitrogen_sc_no)
    deallocate(o2_alpha      ,o2_csurf      ,o2_sc_no)

  end subroutine generic_ERGOM_set_boundary_values

  ! <SUBROUTINE NAME="generic_ERGOM_find_vmove">
  !  <OVERVIEW>
  !   Calculates vertical movement of zooplankton 
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_find_vmove
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_ERGOM_find_vmove(zoo, temp, food, o2, h2s, graz, dzt, ilb, iub, jlb, jub)
    type(zooplankton)           , intent(inout) :: zoo
    real, dimension(ilb:,jlb:,:), intent(inout) :: graz
    real, dimension(ilb:,jlb:,:), intent(in)    :: temp, food, o2, h2s, dzt
    integer                     , intent(in)    :: ilb, iub, jlb, jub
    real, dimension(:,:,:)  , pointer  :: tmask
    integer, dimension(:,:) , pointer  :: kmt
    real                               :: Ieval, rand01, rand02, nue, uexp, Iz, Imax
    real                               :: oeval, o2z, o2min 
    real                               :: temp1, temp2, grad, feval 
    real                               :: tz, topt, tmax, tmig, teval
    integer                            :: i, j, k, km1
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau
    
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=tmask, grid_kmt=kmt)

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!  the light scheme
!  - always sinking
!  - light > Imax     ---> no activity, no grazing (outside this subroutine)
!  - light < Imax     ---> upward motion 
!  - light = Imax     ---> random motion
       Iz   = 2.*ergom%irr_inst(i,j,k) * watt_2_my_einstein !irr_inst is .5 * I_tot for phytopl. 
       Iz   = Iz*Iz
       Imax = zoo%Imax
       Imax = Imax*Imax
       Ieval = Imax*Imax/(Iz*Iz+Imax*Imax)
! no feeding when too much light (Iz>Imax). graz= 0..1
       graz(i,j,k) = Ieval   
! find a gradient following movement
! low  food ---> follow the gradient
! high food ---> random motion
       km1    = max(1,k-1)
       temp1  = zoo%iv*food(i,j,k)
       temp2  = exp(-temp1*temp1)
!       grad   = sign(1.,(food(i,j,km1) - food(i,j,k))) * zoo%iv
       grad   = (food(i,j,km1) - food(i,j,k)) * zoo%iv
       feval = zoo%wfood * temp2 * grad     
       
       o2z  = o2(i,j,k)
       o2min= zoo%o2min
       o2z  = o2z*o2z
       o2min= o2min*o2min
       oeval= zoo%wo2 * o2min/(o2z+o2min)

       tz = max(0.,temp(i,j,k)-zoo%t_opt)
       tz = tz*tz
       tmig = zoo%t_max-zoo%t_opt
       tmig = tmig*tmig
! to_do add zoo%wtemp
       teval = tmig/(tz+tmig)           ! disable rising if temperature is above the optimum temperature

       zoo%p_vmove(i,j,k) = min(Ieval + oeval,teval,1.) * (zoo%wrise0*(1. + feval)) + zoo%wsink0
       zoo%p_vdiff(i,j,k) = zoo%vdiff_max
    enddo; enddo ; enddo  !} i,j,k
  end subroutine generic_ERGOM_find_vmove

  ! <SUBROUTINE NAME="generic_ERGOM_vmove">
  !  <OVERVIEW>
  !   Performs vertical movement (up or down)
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Updates particulate tracer concentrations 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_vmove
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_ERGOM_vmove(move, dzt, field, dt, ilb, iub, jlb, jub)
    real, dimension(ilb:,jlb:,:) ,intent(in)    :: move, dzt
    real, dimension(ilb:,jlb:,:) ,intent(inout) :: field
    real,                   intent(in)          :: dt
    integer,                intent(in)          :: ilb, iub, jlb, jub
    
    real, dimension(:,:,:)  , pointer  :: tmask
    integer, dimension(:,:) , pointer  :: kmt
    real, dimension(ilb:iub,jlb:jub)   :: ft1, ft2
    real                               :: velocity, wpos, wneg
    integer                            :: k, i, j, kp1
    integer                            :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau
    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_vmove'
    
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=tmask, grid_kmt=kmt)

    ft1  = 0.0
    do k = 1, nk-1      
      kp1 = k+1
      do j = jsc, jec ; do i = isc, iec 
        velocity = 0.5*move(i,j,k)
        wpos     = velocity + abs(velocity)
        wneg     = velocity - abs(velocity) 
        ft2(i,j) = (wneg*field(i,j,k) + wpos*field(i,j,kp1)) &
                   *tmask(i,j,k)*tmask(i,j,kp1) 
        field(i,j,k) = field(i,j,k) - dt*tmask(i,j,k)*(ft1(i,j)-ft2(i,j))/dzt(i,j,k)
        ft1(i,j) = ft2(i,j)
      enddo; enddo 
    enddo
    k = nk
    do j = jsc, jec ; do i = isc, iec 
      field(i,j,k) = field(i,j,k) - dt*tmask(i,j,k)*ft1(i,j)/dzt(i,j,k)
    enddo; enddo 
    
  end subroutine generic_ERGOM_vmove
  

  ! <SUBROUTINE NAME="generic_ERGOM_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_ERGOM_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_ERGOM_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_ERGOM_end'
    call user_deallocate_arrays
  end subroutine generic_ERGOM_end
  
  subroutine user_deallocate_arrays
  integer :: n
    deallocate(ergom%f_no3)
    deallocate(ergom%f_nh4)
    deallocate(ergom%f_po4)
    deallocate(ergom%f_o2 )
    deallocate(ergom%f_h2s)
    deallocate(ergom%f_sul)
    deallocate(ergom%f_chl)
    deallocate(ergom%irr_inst)    
    deallocate(ergom%jno3 )
    deallocate(ergom%jnh4 )
    deallocate(ergom%jpo4 )
    deallocate(ergom%jo2  )
    deallocate(ergom%jh2s )
    deallocate(ergom%jsul )
    deallocate(ergom%jdenit_wc)
    deallocate(ergom%jnitrif )
    deallocate(ergom%b_o2 )
    deallocate(ergom%b_no3)
    deallocate(ergom%b_nh4)
    deallocate(ergom%b_po4)
    deallocate(ergom%b_h2s)
    deallocate(ergom%b_nitrogen)
    if(ergom%id_jh2s_o2 .gt. 0) deallocate(ergom%jh2s_o2) 
    if(ergom%id_jh2s_no3.gt. 0) deallocate(ergom%jh2s_no3)
    if(ergom%id_jsul_o2 .gt. 0) deallocate(ergom%jsul_o2) 
    if(ergom%id_jsul_no3.gt. 0) deallocate(ergom%jsul_no3)
    if(ergom%id_jrec_o2 .gt. 0) deallocate(ergom%jrec_o2)
    if(ergom%id_jrec_no3.gt. 0) deallocate(ergom%jrec_no3)
    if(ergom%id_jrec_so4.gt. 0) deallocate(ergom%jrec_so4)
    if(ergom%id_jrec_ana.gt. 0) deallocate(ergom%jrec_ana)
    do n = 1, 2
       deallocate(phyto(n)%f_n)
       if(phyto(n)%id_ilim .gt. 0) deallocate(phyto(n)%ilim)
       deallocate(phyto(n)%jgraz_n)
       deallocate(phyto(n)%jprod_nh4)
       deallocate(phyto(n)%jprod_no3)
       deallocate(phyto(n)%jprod_po4)
       deallocate(phyto(n)%jres_n)	 
       deallocate(phyto(n)%jdet_n)	 
    enddo
    n = DIA
    deallocate(phyto(n)%move)	 
    n = CYA
    deallocate(phyto(n)%move) 
    deallocate(phyto(n)%f_n)
    if(phyto(n)%id_ilim .gt. 0) deallocate(phyto(n)%ilim)
    deallocate(phyto(n)%jgraz_n)
    deallocate(phyto(n)%jprod_po4)
    deallocate(phyto(n)%jres_n)       
    deallocate(phyto(n)%jdet_n)       
    deallocate(phyto(n)%jprod_n2)
    do n = 1, NUM_ZOO
       deallocate(zoo(n)%f_n)	
       deallocate(zoo(n)%jgraz_n)	
       deallocate(zoo(n)%jgain_n)	
       deallocate(zoo(n)%jres_n)	
       deallocate(zoo(n)%jdet_n)	
!       deallocate(zoo(n)%move)	
    enddo
    do n = 1, NUM_DET
       deallocate(det(n)%f_n)
       deallocate(det(n)%jgraz_n)
       deallocate(det(n)%jmort)
    enddo
    deallocate(biosed%bioerosion)
    do n = 1, NUM_SPM
        deallocate(spm(n)%move)
        deallocate(spm(n)%btf)
        deallocate(spm(n)%jsed)
    enddo
    do n = 1, NUM_SED
        deallocate(sed(n)%f_sed)
        deallocate(sed(n)%jgain_sed)
        deallocate(sed(n)%jloss_sed)
        deallocate(sed(n)%jres)
        deallocate(sed(n)%jbiores)
    end do
    if (biosed%id_jrec_n .gt. 0)     deallocate(biosed%jrec_n)
    if (biosed%id_jdenit_sed .gt. 0) deallocate(biosed%jdenit_sed)
    if (biosed%id_mode_sed .gt. 0)   deallocate(biosed%mode_sed)
    deallocate(tracers_2d)
  end subroutine user_deallocate_arrays

end module generic_ERGOM
