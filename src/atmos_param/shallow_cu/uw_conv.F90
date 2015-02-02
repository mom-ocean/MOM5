#include <fms_platform.h>
MODULE UW_CONV_MOD

  use           mpp_mod, only : mpp_pe, mpp_root_pe, stdlog
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use   Time_Manager_Mod, ONLY: time_type, get_time 
  use           mpp_mod, only : input_nml_file
  use           fms_mod, only : write_version_number, open_namelist_file, check_nml_error,&
                                FILE_EXIST, ERROR_MESG,  &
                                lowercase, &
                                CLOSE_FILE, FATAL
  use  field_manager_mod, only: MODEL_ATMOS
  use  tracer_manager_mod, only: get_tracer_names, query_method, &
                                 get_tracer_index, NO_TRACER
  use  sat_vapor_pres_mod,only : sat_vapor_pres_init
  use atmos_tracer_utilities_mod, only : get_wetdep_param

  use  rad_utilities_mod, only : aerosol_type
  
  use  aer_ccn_act_mod, only :   aer_ccn_act_init
  use  conv_utilities_mod,only :   uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, &
                                   check_tracer_realizability, &
                                   qt_parcel_k, &
                                   adicloud, sounding, uw_params

  use  conv_plumes_k_mod,only    : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist, cwetdep_type

  use  conv_closures_mod,only    : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure

  use  deep_conv_mod,only        : deepc, cpn_copy, dpconv0, dpconv1, dpconv2

!---------------------------------------------------------------------
  implicit none
  private
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: uw_conv.F90,v 20.0 2013/12/13 23:21:41 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: uw_conv, uw_conv_init, uw_conv_end

  real, parameter :: aday = 1.
  real, parameter :: mv = -999.
  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'uw_conv'

  !namelist parameters for UW convection scheme
  integer :: iclosure = 0      ! 0: Bretherton UWShCu orginal / -CIN/TKE based
                               ! 1: Emanuel-Rayment: quasiequilibrium PBL
  real    :: rkm_sh   = 16.0   ! fractional lateral mixing rate for shallow
  real    :: cldhgt_max   = 4.e3
  real    :: landfact_m   = 0.0
  integer :: idpchoice = 0  
  logical :: do_deep = .false.
  logical :: do_relaxcape = .false.
  logical :: do_relaxwfn  = .false.
  logical :: do_coldT = .true.
  logical :: do_lands = .false.
  logical :: do_uwcmt = .false.   
  logical :: do_fast  = .false.
  logical :: do_ice   = .true.
  logical :: do_ppen  = .true.
  logical :: do_forcedlifting = .false.
  real    :: atopevap = 0.
  logical :: apply_tendency = .true.
  logical :: prevent_unreasonable = .true.
  real    :: aerol = 1.e-12
  real    :: gama     = 1.0    ! 
  real    :: tkemin   = 1.e-6
  real    :: wmin_ratio = 0.
  logical :: use_online_aerosol = .false.
  logical :: use_sub_seasalt = .true.
  logical :: do_auto_aero = .false.
  logical :: do_rescale   = .false.
  logical :: do_debug     = .false.
!miz
  logical :: do_imposing_cooling_drying = .false.
  real    :: tdt_rate = 0.0             
  real    :: qvdt_rate = 0.0
  real    :: pres_min = 0.0
  real    :: pres_max = 0.0
  logical :: do_imposing_rad_cooling = .false.
  logical :: do_no_uw_conv           = .false.
  real    :: cooling_rate = -1.5 !K/day
  real    :: t_thresh = 207.5    !K
  real    :: t_strato = 200.0    !K
  real    :: tau_rad  = 5.0      !day
!miz
  integer :: cush_choice  = 0
  real    :: pcp_min      = 3e-5
  real    :: pcp_max      = 1.5e-3
  real    :: rh0          = 0.8
  real    :: cush_ref     = 0.
  real    :: pblht0 = 500.
  real    :: tke0 = 1.
  real    :: lofactor0 = 1.
  integer :: lochoice  = 0
  real    :: wrel_min = 1.
  real    :: om_to_oc = 1.67
  real    :: sea_salt_scale = 0.1
  logical :: do_qctflx_zero = .false.
  logical :: do_detran_zero = .false.
  logical :: do_gust_cv     = .false.
  real    :: gustmax        = 3.! maximum gustiness wind (m/s)
  real    :: gustconst      = 10./86400.   ! constant in kg/m2/sec, default =
                                           ! 1 cm/day = 10 mm/day

  NAMELIST / uw_conv_nml / iclosure, rkm_sh, cldhgt_max, &
       do_deep, idpchoice, do_relaxcape, do_relaxwfn, do_coldT, do_lands, do_uwcmt,       &
       do_fast, do_ice, do_ppen, do_forcedlifting, &
       atopevap, apply_tendency, prevent_unreasonable, aerol, gama, tkemin,    &
       wmin_ratio, use_online_aerosol, use_sub_seasalt, landfact_m, pblht0, tke0, lofactor0, lochoice, &
       do_auto_aero, do_rescale, wrel_min, om_to_oc, sea_salt_scale,                  &
       do_debug, cush_choice, pcp_min, pcp_max, cush_ref,   &
       rh0, do_qctflx_zero, do_detran_zero, do_gust_cv, gustmax, gustconst, &
       do_imposing_cooling_drying, tdt_rate, qvdt_rate, pres_min, pres_max, &
       do_imposing_rad_cooling, do_no_uw_conv, cooling_rate, t_thresh, t_strato, tau_rad

  !namelist parameters for UW convective plume
  real    :: rle      = 0.10   ! for critical stopping distance for entrainment
  real    :: rpen     = 5.0    ! for entrainment efficiency
  real    :: rmaxfrac = 0.05   ! maximum allowable updraft fraction
  real    :: wmin     = 0.0    ! maximum allowable updraft fraction
  real    :: rbuoy    = 1.0    ! for nonhydrostatic pressure effects on updraft
  real    :: rdrag    = 1.0 
  real    :: frac_drs = 0.0    ! 
  real    :: bigc     = 0.7    ! for momentum transfer
  real    :: auto_th0 = 0.5e-3 ! threshold for precipitation
  real    :: auto_rate= 1.e-3
  real    :: tcrit    = -45.0  ! critical temperature 
  real    :: deltaqc0 = 0.5e-3 
  logical :: do_pdfpcp= .false.
  logical :: do_pmadjt= .false.
  logical :: do_emmax = .false.
  logical :: do_pnqv  = .false.
  real    :: rad_crit = 14.0   ! critical droplet radius
  real    :: emfrac_max = 1.0
  integer :: mixing_assumption = 0
  integer :: mp_choice = 1
  real    :: Nl_land   = 300.e6
  real    :: Nl_ocean  = 100.e6
  real    :: qi_thresh = 1.e-4
  real    :: r_thresh  = 12.e-6
  logical :: do_pevap = .false.
  real    :: cfrac     = 0.05
  real    :: hcevap    = 0.8
  logical :: do_weffect = .false.
  real    :: weffect    = 0.5
  real    :: peff_l     = 1.0
  real    :: peff_i     = 1.0
  real    :: t00        = 295

  NAMELIST / uw_plume_nml / rle, rpen, rmaxfrac, wmin, rbuoy, rdrag, frac_drs, bigc, &
       auto_th0, auto_rate, tcrit, deltaqc0, do_pdfpcp, do_pmadjt, do_emmax, do_pnqv, rad_crit, emfrac_max, &
       mixing_assumption, mp_choice, Nl_land, Nl_ocean, qi_thresh, r_thresh, do_pevap, cfrac, hcevap, &
       do_weffect, weffect, peff_l, peff_i, t00
  !namelist parameters for UW convective closure
  integer :: igauss   = 1      ! options for cloudbase massflux closure
                               ! 1: cin/gaussian closure, using TKE to compute CIN.
                               ! 2: cin/gaussian closure, using W* to compute CIN.
                               ! 0: cin/tke mapse-style closure; 
  real    :: rkfre    = 0.05   ! vertical velocity variance as fraction of tke
  real    :: tau_sh   = 7200.  ! 
  real    :: wcrit_min= 0.

  NAMELIST / uw_closure_nml / igauss, rkfre, tau_sh, wcrit_min


!========Option for deep convection=======================================
  real    :: rkm_dp1       = 10.
  real    :: rkm_dp2       = 1.
  real    :: cbmf_dp_frac1 = 0.
  real    :: cbmf_dp_frac2 = 1.
  real    :: crh_th_ocean  = 100.
  real    :: crh_th_land   = 100.
  real    :: cape_th       = 0.
  real    :: tau_dp        = 7200.
  real    :: rpen_d        = 5.0
  integer :: mixing_assumption_d = 0
  integer :: norder      = 1
  logical :: do_ppen_d   = .true.
  logical :: do_pevap_d  = .true.
  real    :: cfrac_d     = 0.05
  real    :: hcevap_d    = 0.8
  real    :: dcapedm_th  = 0
  real    :: frac_limit_d = 0.25
  real    :: lofactor_d   = 1.0
  real    :: auto_th0_d   = 1.0e-3
  real    :: tcrit_d      = -120
  logical :: do_forcedlifting_d = .false.
  logical :: do_lod_rkm   = .false.
  logical :: do_lod_cfrac = .false.
  logical :: do_lod_tcrit = .false.

  NAMELIST / deep_conv_nml / rkm_dp1, rkm_dp2, cbmf_dp_frac1, cbmf_dp_frac2, &
                 crh_th_ocean, crh_th_land, do_forcedlifting_d, frac_limit_d, &
                 cape_th, tau_dp, rpen_d, mixing_assumption_d, norder, &
                 do_ppen_d, do_pevap_d, cfrac_d, hcevap_d, lofactor_d, dcapedm_th, &
                 auto_th0_d, tcrit_d, do_lod_rkm, do_lod_cfrac, do_lod_tcrit
!========Option for deep convection=======================================

!------------------------------------------------------------------------

  integer :: nqv, nql, nqi, nqa ,nqn
  logical :: do_qn = .false.    ! use droplet number tracer field ?

  integer :: id_tdt_uwc, id_qdt_uwc, id_prec_uwc, id_snow_uwc,               &
       id_cin_uwc, id_cbmf_uwc, id_tke_uwc, id_plcl_uwc, id_zinv_uwc,  &
       id_cush_uwc, id_pct_uwc, id_pcb_uwc, id_plfc_uwc, id_enth_uwc,  &
       id_qldt_uwc, id_qidt_uwc, id_qadt_uwc, id_qndt_uwc, id_cmf_uwc, id_wu_uwc,   &
       id_fer_uwc,  id_fdr_uwc, id_fdrs_uwc, id_cqa_uwc, id_cql_uwc,   &
       id_cqi_uwc,  id_cqn_uwc, id_hlflx_uwc, id_qtflx_uwc,           &
       id_cape_uwc, id_dcin_uwc, id_dcape_uwc, id_dwfn_uwc, id_crh_uwc,&
       id_ocode_uwc, id_plnb_uwc, id_wrel_uwc, id_ufrc_uwc, id_qtmp_uwc,&
       id_tdt_pevap_uwc, id_qdt_pevap_uwc, id_xhlsrc_uwc, id_xqtsrc_uwc,&
       id_qldet_uwc, id_qidet_uwc, id_qadet_uwc, id_qtdt_uwc, id_dting_uwc, &
       id_cfq_uwc, id_fdp_uwc, id_hmo_uwc, id_hms_uwc, id_abu_uwc, id_peo_uwc, &
       id_tdt_rad_uwc


  integer, allocatable :: id_tracerdt_uwc(:), id_tracerdt_uwc_col(:), &
                          id_tracerdtwet_uwc(:), id_tracerdtwet_uwc_col(:)

!========Option for deep convection=======================================
  integer :: id_tdt_uwd, id_qdt_uwd, id_qtdt_uwd, id_prec_uwd, id_snow_uwd,   &
       id_cbmf_uwd, id_enth_uwd, id_qldt_uwd, id_qidt_uwd,             &
       id_qndt_uwd, id_qadt_uwd, id_cmf_uwd, id_wu_uwd, id_fer_uwd,    &
       id_fdr_uwd, id_fdrs_uwd, id_cqa_uwd, id_cql_uwd, id_cqi_uwd,    &
       id_cqn_uwd, id_hlflx_uwd, id_qtflx_uwd, id_dcin_uwd,            &
       id_dcapedm_uwd, id_dwfn_uwd, id_ocode_uwd,                      &
       id_tdt_pevap_uwd, id_qdt_pevap_uwd, id_rkm_uwd, id_cbu_uwd
!========Option for deep convection=======================================

  type(cwetdep_type), dimension(:), allocatable :: wetdep
  type(uw_params),  save  :: Uw_p
  character(len=32), dimension(:), allocatable   :: tracername 
  character(len=32), dimension(:), allocatable   :: tracer_units 

contains

!#####################################################################
!#####################################################################

  SUBROUTINE UW_CONV_INIT(do_strat, axes, Time, kd, tracers_in_uw)
    logical,         intent(in) :: do_strat
    integer,         intent(in) :: axes(4), kd
    type(time_type), intent(in) :: Time
    logical,         intent(in) :: tracers_in_uw(:)
    
!---------------------------------------------------------------------
!  intent(in) variables:
!
!      tracers_in_uw 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by UW convection
!
!-------------------------------------------------------------------

    integer   :: unit, io
    
    integer   :: ntracers, n, nn, ierr, logunit
    logical   :: flag
    character(len=200) :: text_in_scheme, control
     real :: frac_junk, frac_junk2
 
    ntracers = count(tracers_in_uw)

    call uw_params_init   (Uw_p)

!   Initialize lookup tables needed for findt and exn
!   sat_vapor_pres needs to be initialized if not already done
    call sat_vapor_pres_init     
    call exn_init_k (Uw_p)
    call findt_init_k (Uw_p)

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=uw_closure_nml, iostat=io)
      ierr = check_nml_error(io,'uw_closure_nml')
      read (input_nml_file, nml=uw_conv_nml, iostat=io)
      ierr = check_nml_error(io,'uw_conv_nml')
      read (input_nml_file, nml=uw_plume_nml, iostat=io)
      ierr = check_nml_error(io,'uw_plume_nml')
      read (input_nml_file, nml=deep_conv_nml, iostat=io)
      ierr = check_nml_error(io,'deep_conv_nml')
#else   
    if( FILE_EXIST( 'input.nml' ) ) then
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_closure_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'uw_closure_nml')
       end do
10     call close_file ( unit )

       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_conv_nml, iostat = io, end = 20 )
          ierr = check_nml_error(io,'uw_conv_nml')
       end do
20     call close_file ( unit )
       
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_plume_nml, iostat = io, end = 30 )
          ierr = check_nml_error(io,'uw_plume_nml')
       end do
30     call close_file ( unit )

!========Option for deep convection=======================================
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = deep_conv_nml, iostat = io, end = 40 )
          ierr = check_nml_error(io,'deep_conv_nml')
       end do
40     call close_file ( unit )
!========Option for deep convection=======================================
    end if
#endif
    call write_version_number (version, tagname)
    logunit = stdlog()
    WRITE( logunit, nml = uw_closure_nml )
    WRITE( logunit, nml = uw_conv_nml )
    WRITE( logunit, nml = uw_plume_nml )
    WRITE( logunit, nml = deep_conv_nml )

    if ( use_online_aerosol ) call aer_ccn_act_init

    nqv = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
    nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
    nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
    nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
    if (nqn /= NO_TRACER) do_qn = .true.
    if (ntracers > 0) then
      allocate ( tracername   (ntracers) )
      allocate ( tracer_units (ntracers) )
      allocate ( wetdep       (ntracers) )
      nn = 1
      do n=1,size(tracers_in_uw(:))
         if (tracers_in_uw(n)) then
             call get_tracer_names (MODEL_ATMOS, n,  &
                                    name = tracername(nn), &
                                    units = tracer_units(nn))
             flag = query_method( 'wet_deposition', MODEL_ATMOS, n, &
                                  text_in_scheme, control )
             call get_wetdep_param( text_in_scheme, control, &
                                    wetdep(nn)%scheme, &
                                    wetdep(nn)%Henry_constant, &
                                    wetdep(nn)%Henry_variable, &
                                    frac_junk, frac_junk2, &
                                    wetdep(nn)%alpha_r, &
                                    wetdep(nn)%alpha_s, &
                                    wetdep(nn)%Lwetdep, &
                                    wetdep(nn)%Lgas, &
                                    wetdep(nn)%Laerosol, &
                                    wetdep(nn)%Lice, &
                                    frac_in_cloud_uw=wetdep(nn)%frac_in_cloud )
             wetdep(nn)%scheme = lowercase( wetdep(nn)%scheme )
             nn = nn + 1
          endif
       end do
    endif

    id_xhlsrc_uwc = register_diag_field (mod_name,'xhlsrc_uwc', axes(1:2), Time, &
         'xhlsrc', 'J/kg' )
    id_xqtsrc_uwc = register_diag_field (mod_name,'xqtsrc_uwc', axes(1:2), Time, &
         'xqtsrc', 'kg/kg' )

    id_tdt_pevap_uwc = register_diag_field ( mod_name, 'tdt_pevap_uwc', axes(1:3), Time, &
         'Temperature tendency due to pevap from uw_conv', 'K/s', missing_value=mv )
    id_qdt_pevap_uwc = register_diag_field ( mod_name, 'qdt_pevap_uwc', axes(1:3), Time, &
         'Spec. humidity tendency due to pevap from uw_conv', 'kg/kg/s', missing_value=mv)

    id_tdt_uwc = register_diag_field ( mod_name, 'tdt_uwc', axes(1:3), Time, &
         'Temperature tendency from uw_conv', 'K/s', missing_value=mv )
    id_qdt_uwc = register_diag_field ( mod_name, 'qdt_uwc', axes(1:3), Time, &
         'Spec. humidity tendency from uw_conv', 'kg/kg/s', missing_value=mv)
    id_cmf_uwc = register_diag_field ( mod_name, 'cmf_uwc', axes(1:3), Time, &
         'Cloud vert. mass flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_cfq_uwc = register_diag_field ( mod_name, 'cfq_uwc', axes(1:3), Time,   &
         'Convective frequency', 'none', missing_value=mv)
    id_peo_uwc = register_diag_field ( mod_name, 'peo_uwc', axes(1:3), Time,   &
         'Convective precipitation efficiency', 'none', missing_value=mv)
    id_hmo_uwc = register_diag_field ( mod_name, 'hmo_uwc', axes(1:3), Time,   &
         'moist static energy', 'J/kg', missing_value=mv)
    id_hms_uwc = register_diag_field ( mod_name, 'hms_uwc', axes(1:3), Time,   &
         'moist static energy', 'J/kg', missing_value=mv)
    id_abu_uwc = register_diag_field ( mod_name, 'abu_uwc', axes(1:3), Time,   &
         'adiabatic buoyancy', 'K', missing_value=mv)
    id_wu_uwc = register_diag_field ( mod_name, 'wu_uwc', axes(1:3), Time,   &
         'Updraft vert. velocity from uw_conv', 'm/s', missing_value=mv)
    id_fer_uwc = register_diag_field ( mod_name, 'fer_uwc', axes(1:3), Time, &
         'Fractional entrainment rate from uw_conv', '1/Pa', missing_value=mv)
    id_fdr_uwc = register_diag_field ( mod_name, 'fdr_uwc', axes(1:3), Time, &
         'Fractional detrainment rate from uw_conv', '1/Pa', missing_value=mv)
    id_fdrs_uwc = register_diag_field (mod_name,'fdrs_uwc', axes(1:3), Time, &
         'Detrainment rate for sat. air from uw_conv', '1/Pa', missing_value=mv)
    id_cqa_uwc = register_diag_field ( mod_name, 'cqa_uwc', axes(1:3), Time, &
         'Updraft fraction from uw_conv', 'none', missing_value=mv)
    id_cql_uwc = register_diag_field ( mod_name, 'cql_uwc', axes(1:3), Time, &
         'Updraft liquid from uw_conv', 'kg/kg', missing_value=mv)
    id_cqi_uwc = register_diag_field ( mod_name, 'cqi_uwc', axes(1:3), Time, &
         'Updraft ice from uw_conv', 'kg/kg', missing_value=mv)
    id_cqn_uwc = register_diag_field ( mod_name, 'cqn_uwc', axes(1:3), Time, &
         'Updraft liquid drop from uw_conv', '/kg', missing_value=mv)
    id_hlflx_uwc=register_diag_field (mod_name,'hlflx_uwc',axes(1:3),Time, &
         'Liq.wat.pot.temp. flux from uw_conv', 'W/m2', missing_value=mv)
    id_qtflx_uwc = register_diag_field (mod_name,'qtflx_uwc',axes(1:3),Time, &
         'Total water flux from uw_conv', 'W/m2', missing_value=mv)
    id_prec_uwc = register_diag_field (mod_name,'prec_uwc', axes(1:2), Time, &
         'Precipitation rate from uw_conv', 'kg/m2/sec',                     &
         interp_method = "conserve_order1" )
    id_snow_uwc = register_diag_field (mod_name,'snow_uwc', axes(1:2), Time, &
         'Frozen precip. rate from uw_conv', 'kg/m2/sec',                       &
         interp_method = "conserve_order1" )
    id_cin_uwc = register_diag_field ( mod_name, 'cin_uwc', axes(1:2), Time, &
         'CIN from uw_conv', 'm2/s2' )
    id_cape_uwc= register_diag_field ( mod_name,'cape_uwc', axes(1:2), Time, &
         'CAPE from uw_conv', 'm2/s2' )
    id_crh_uwc= register_diag_field ( mod_name,'crh_uwc', axes(1:2), Time, &
         'Column RH from uw_conv', '%' )
    id_cbmf_uwc = register_diag_field (mod_name,'cbmf_uwc', axes(1:2), Time, &
         'Cloud-base mass flux from uw_conv', 'kg/m2/s' )
    id_wrel_uwc = register_diag_field (mod_name,'wrel_uwc', axes(1:2), Time, &
         'Release level vertical velocity from uw_conv', 'm/s' )
    id_ufrc_uwc = register_diag_field (mod_name,'ufrc_uwc', axes(1:2), Time, &
         'Release level updraft fraction from uw_conv', 'none' )
    id_tke_uwc = register_diag_field ( mod_name, 'tke_uwc', axes(1:2), Time, &
         'PBL mean TKE from uw_conv', 'm2/s2' )
    id_plcl_uwc = register_diag_field (mod_name,'plcl_uwc', axes(1:2), Time, &
         'LCL pressure from uw_conv', 'hPa' )
    id_plfc_uwc = register_diag_field (mod_name,'plfc_uwc', axes(1:2), Time, &
         'LFC pressure from uw_conv', 'hPa' )
    id_plnb_uwc = register_diag_field (mod_name,'plnb_uwc', axes(1:2), Time, &
         'LNB pressure from uw_conv', 'hPa' )
    id_zinv_uwc = register_diag_field (mod_name,'zinv_uwc', axes(1:2), Time, &
         'Inversion pressure from uw_conv', 'm' )
    id_pct_uwc = register_diag_field ( mod_name, 'pct_uwc', axes(1:2), Time, &
         'Cloud-top pressure from uw_conv', 'hPa' )
    id_pcb_uwc = register_diag_field ( mod_name, 'pcb_uwc', axes(1:2), Time, &
         'Cloud-base pressure from uw_conv', 'hPa' )
    id_cush_uwc = register_diag_field (mod_name,'cush_uwc', axes(1:2), Time, &
         'Convective scale height from uw_conv', 'm' )
    id_dcin_uwc = register_diag_field (mod_name, 'dcin_uwc', axes(1:2), Time, &
         'dCIN/cbmf from uw_conv', 'm2/s2/(kg/m2/s)' )
    id_dcape_uwc= register_diag_field (mod_name, 'dcape_uwc', axes(1:2), Time, &
         'dCAPE/cbmf from uw_conv', 'm2/s2/(kg/m2/s)' )
    id_dwfn_uwc = register_diag_field (mod_name, 'dwfn_uwc',  axes(1:2), Time, &
         'dwfn/cbmf from uw_conv', '(m2/s2)/(kg/m2/s)' )
    id_enth_uwc = register_diag_field (mod_name,'enth_uwc', axes(1:2), Time, &
         'Column-integrated enthalpy tendency from uw_conv', 'W/m2' )
    id_qtmp_uwc = register_diag_field (mod_name,'qtmp_uwc', axes(1:2), Time, &
         'Column-integrated water tendency from uw_conv', 'kg/m2/s' )
    id_dting_uwc = register_diag_field (mod_name,'dting_uwc', axes(1:2), Time, &
         'Column-integrated heating rate from uw_conv', 'W/m2' )
    id_ocode_uwc = register_diag_field (mod_name,'ocode_uwc', axes(1:2), Time, &
         'Out code from uw_conv', 'none' )
    id_fdp_uwc = register_diag_field ( mod_name, 'fdp_uwc',   axes(1:2), Time,   &
         'Deep convective frequency', 'none', missing_value=mv)
    if ( do_strat ) then
       id_qldt_uwc= register_diag_field (mod_name,'qldt_uwc',axes(1:3),Time, &
            'Liquid water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qidt_uwc= register_diag_field (mod_name,'qidt_uwc',axes(1:3),Time, &
            'Ice water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qadt_uwc= register_diag_field (mod_name,'qadt_uwc',axes(1:3),Time, &
            'CLD fraction tendency from uw_conv', '1/s', missing_value=mv )
       id_qndt_uwc= register_diag_field (mod_name,'qndt_uwc',axes(1:3),Time, &
            'Cloud droplet number fraction tendency from uw_conv', '#/kg/s', missing_value=mv )
       id_qldet_uwc = register_diag_field (mod_name,'qldet_uwc',axes(1:3),Time, &
            'ql detrainment', 'kg/kg/s', missing_value=mv)
       id_qidet_uwc = register_diag_field (mod_name,'qidet_uwc',axes(1:3),Time, &
            'qi detrainment', 'kg/kg/s', missing_value=mv)
       id_qadet_uwc = register_diag_field (mod_name,'qadet_uwc',axes(1:3),Time, &
            'qa detrainment', '1/s', missing_value=mv)
       id_qtdt_uwc= register_diag_field (mod_name,'qtdt_uwc',axes(1:3),Time, &
            'Total water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
    end if

    if (do_imposing_rad_cooling) then
       id_tdt_rad_uwc = register_diag_field ( mod_name, 'tdt_rad_uwc', axes(1:3), Time, &
         'Idealized radiative temperature tendency from uw_conv', 'K/s', missing_value=mv )
    end if

!========Option for deep convection=======================================
    if (do_deep) then
       id_tdt_pevap_uwd = register_diag_field ( mod_name, 'tdt_pevap_uwd', axes(1:3), Time, &
            'Temperature tendency due to pevap from deep_conv', 'K/s', missing_value=mv )
       id_qdt_pevap_uwd = register_diag_field ( mod_name, 'qdt_pevap_uwd', axes(1:3), Time, &
            'Spec. humidity tendency due to pevap from deep_conv', 'kg/kg/s', missing_value=mv)

       id_tdt_uwd = register_diag_field ( mod_name, 'tdt_uwd', axes(1:3), Time, &
            'Temperature tendency from deep_conv', 'K/s', missing_value=mv )
       id_qdt_uwd = register_diag_field ( mod_name, 'qdt_uwd', axes(1:3), Time, &
            'Spec. humidity tendency from deep_conv', 'kg/kg/s', missing_value=mv)
       id_qtdt_uwd= register_diag_field ( mod_name, 'qtdt_uwd', axes(1:3), Time, &
            'Total water spec. humidity tendency from deep_conv', 'kg/kg/s', missing_value=mv)
       id_cmf_uwd = register_diag_field ( mod_name, 'cmf_uwd', axes(1:3), Time, &
            'Cloud vert. mass flux from deep_conv', 'kg/m2/s', missing_value=mv)
       id_wu_uwd = register_diag_field ( mod_name, 'wu_uwd', axes(1:3), Time,   &
            'Updraft vert. velocity from deep_conv', 'm/s', missing_value=mv)
       id_cbu_uwd= register_diag_field ( mod_name, 'cbu_uwd', axes(1:3), Time,   &
            'deep plume buoyancy', 'K', missing_value=mv)
       id_fer_uwd = register_diag_field ( mod_name, 'fer_uwd', axes(1:3), Time, &
         'Fractional entrainment rate from deep_conv', '1/Pa', missing_value=mv)
       id_fdr_uwd = register_diag_field ( mod_name, 'fdr_uwd', axes(1:3), Time, &
            'Fractional detrainment rate from deep_conv', '1/Pa', missing_value=mv)
       id_fdrs_uwd = register_diag_field (mod_name,'fdrs_uwd', axes(1:3), Time, &
            'Detrainment rate for sat. air from deep_conv', '1/Pa', missing_value=mv)
       id_cqa_uwd = register_diag_field ( mod_name, 'cqa_uwd', axes(1:3), Time, &
            'Updraft fraction from deep_conv', 'none', missing_value=mv)
       id_cql_uwd = register_diag_field ( mod_name, 'cql_uwd', axes(1:3), Time, &
         'Updraft liquid from deep_conv', 'kg/kg', missing_value=mv)
       id_cqi_uwd = register_diag_field ( mod_name, 'cqi_uwd', axes(1:3), Time, &
            'Updraft ice from deep_conv', 'kg/kg', missing_value=mv)
       id_cqn_uwd = register_diag_field ( mod_name, 'cqn_uwd', axes(1:3), Time, &
            'Updraft liquid drop from deep_conv', '/kg', missing_value=mv)
       id_hlflx_uwd=register_diag_field (mod_name,'hlflx_uwd',axes(1:3),Time, &
            'Liq.wat.pot.temp. flux from deep_conv', 'W/m2', missing_value=mv)
       id_qtflx_uwd = register_diag_field (mod_name,'qtflx_uwd',axes(1:3),Time, &
            'Total water flux from deep_conv', 'W/m2', missing_value=mv)
       id_prec_uwd = register_diag_field (mod_name,'prec_uwd', axes(1:2), Time, &
            'Precipitation rate from deep_conv', 'kg/m2/sec' )
       id_snow_uwd = register_diag_field (mod_name,'snow_uwd', axes(1:2), Time, &
            'Frozen precip. rate from deep_conv', 'kg/m2/sec' )
       id_cbmf_uwd = register_diag_field (mod_name,'cbmf_uwd', axes(1:2), Time, &
            'Cloud-base mass flux from deep_conv', 'kg/m2/s' )
       id_dcapedm_uwd= register_diag_field (mod_name, 'dcapedm_uwd', axes(1:2), Time, &
            'dCAPE/cbmf from deep_conv', 'm2/s2/(kg/m2/s)' )
       id_dwfn_uwd = register_diag_field (mod_name, 'dwfn_uwd',  axes(1:2), Time, &
            'dwfn/cbmf from deep_conv', '(m2/s2)/(kg/m2/s)' )
       id_enth_uwd = register_diag_field (mod_name,'enth_uwd', axes(1:2), Time, &
            'Column-integrated enthalpy tendency from deep_conv', 'K/s' )
       id_ocode_uwd = register_diag_field (mod_name,'ocode_uwd', axes(1:2), Time, &
            'Out code from deep_conv', 'none' )
       id_rkm_uwd = register_diag_field (mod_name,'rkm_uwd', axes(1:2), Time, &
            'rkm for deep_conv', 'none' )
       if ( do_strat ) then
          id_qldt_uwd= register_diag_field (mod_name,'qldt_uwd',axes(1:3),Time, &
               'Liquid water tendency from deep_conv', 'kg/kg/s', missing_value=mv)
          id_qidt_uwd= register_diag_field (mod_name,'qidt_uwd',axes(1:3),Time, &
               'Ice water tendency from deep_conv', 'kg/kg/s', missing_value=mv)
          id_qadt_uwd= register_diag_field (mod_name,'qadt_uwd',axes(1:3),Time, &
               'CLD fraction tendency from deep_conv', '1/s', missing_value=mv )
       end if
    end if
!========Option for deep convection=======================================


    if ( ntracers>0 ) then
      allocate(id_tracerdt_uwc(ntracers), id_tracerdt_uwc_col(ntracers) )
       allocate(id_tracerdtwet_uwc(ntracers), id_tracerdtwet_uwc_col(ntracers))
      do nn = 1,ntracers
         id_tracerdt_uwc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdt_uwc_col(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_col', &
                                     axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
           id_tracerdtwet_uwc(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_wet', &
                                    axes(1:3), Time, &
                                   trim(tracername(nn)) //' tendency from uw_conv wetdep', &
                                   trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdtwet_uwc_col(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_wet_col', &
                                   axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv wetdep', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
        end do
     end if

    module_is_initialized = .true.

    
  end SUBROUTINE UW_CONV_INIT

!#####################################################################
!#####################################################################

  subroutine uw_conv_end
    call exn_end_k
    call findt_end_k
    module_is_initialized = .FALSE.
  end subroutine uw_conv_end

!#####################################################################
!#####################################################################

  SUBROUTINE uw_conv(is, js, Time, tb, qv, ub, vb, pmid, pint,zmid,  & !input
       zint, q, omega, delt, pblht, ustar, bstar, qstar, land, coldT,& !input
       asol,                                                         & !input
       cush, do_strat,  skip_calculation, max_available_cf,          & !input
       tten, qvten, qlten, qiten, qaten, qnten,                      & !output
       uten, vten, rain, snow,                                       & !output
       cmf, hlflx, qtflx, pflx, liq_pflx, ice_pflx, cldql, cldqi, cldqa,cldqn, cbmfo,  & !output
        tracers, trtend, uw_wetdep)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt 

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: q     !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    real, intent(in), dimension(:,:,:)   :: max_available_cf !  largest
                                     ! realizable value for uw cld frac
                                   ! after accounting for deep cld frac
    logical,intent(in), dimension(:,:)   :: skip_calculation ! do not
                                                 ! calculate where .true.
    logical,intent(in)                   :: do_strat !logical flag
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar !pbl height...
    real, intent(inout), dimension(:,:)  :: cush  ! convective scale height (m) 

    type(aerosol_type),  intent (in)     :: asol
   
    real, intent(out), dimension(:,:,:)  :: tten,qvten              ! T,qv tendencies
    real, intent(out), dimension(:,:,:)  :: qlten,qiten,qaten,qnten ! q tendencies
    real, intent(out), dimension(:,:,:)  :: uten,vten               ! u,v tendencies
   
    real, intent(out), dimension(:,:,:)  :: cldql,cldqi,cldqa, cldqn!in-updraft q
    real, intent(out), dimension(:,:,:)  :: cmf    ! mass flux at level above layer (kg/m2/s)
    real, intent(out), dimension(:,:,:)  :: pflx   ! precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: liq_pflx   ! liq precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: ice_pflx   ! solid precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: hlflx ! theta_l flux
    real, intent(out), dimension(:,:,:)  :: qtflx  ! qt  flux
    real, intent(out), dimension(:,:)    :: rain, snow
    real, intent(inout), dimension(:,:)  :: cbmfo  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers
    real, intent(out), dimension(:,:,:,:)  :: trtend          ! calculated tracer tendencies
    real, intent(out), dimension(:,:,:)  :: uw_wetdep       ! calculated wet depostion for tracers

    integer i, j, k, kl, klm, nk, naer, na, n

    real rhos0j
    real hlsrc, thcsrc, qctsrc, tmp, lofactor, crh_th
    real zsrc, psrc, cbmf_shallow, cbmf_old, cbmf_deep, rkm_sh1, rkm_dp, cbmf_dp_frac, dcrh, dcrh0, dpsum
    real, dimension(size(tb,1),size(tb,2)) :: &
         plcl,       &     ! pressure of lifting condensation level (Pa)
         plfc,       &     ! pressure of level of free convection (Pa)
         plnb,       &     ! pressure of level of neutral buoyancy (Pa)
         cino,       &     ! cin (m2/s2)
         capeo,      &     ! cape(m2/s2)
         tkeo,       &     ! tke (m2/s2)
         wrelo,      &     ! release level vertical velocity (m/s)
         ufrco,      &     ! cloud-base updraft fraction
         zinvo,      &     ! surface driven mixed-layer height
         denth,      &     
         dqtmp,      &
         dting,      &
         dcino,      &     ! dcin (m2/s2)
         dcapeo,     &     ! dcape(m2/s2)
         dwfno,      &     ! dwfn(m2/s2)
         ocode,      &
         xhlsrc,     &
         xqtsrc,     &
         crho,       &
         rkm_d,      &
         fdp

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: wuo,fero,fdro,fdrso, tten_pevap, qvten_pevap
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: qldet, qidet, qadet, cfq, peo, hmo, hms, abu

    real, dimension(size(tb,1),size(tb,2))            :: scale_uw
    real :: qtin, dqt, temp_1

!========Option for deep convection=======================================
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: uten_d, vten_d, tten_d, &
         qvten_d, qlten_d, qiten_d, qaten_d, qnten_d, cmf_d, cbu_d, pflx_d, hlflx_d, qtflx_d, qtten_d, &
         wuo_d, fero_d, fdro_d, fdrso_d, cldql_d, cldqi_d, cldqa_d, cldqn_d, tten_pevap_d, qvten_pevap_d
    real, dimension(size(tb,1),size(tb,2)) :: rain_d, snow_d, dcapedm_d, dwfn_d, denth_d, dting_d, dqtmp_d, cbmf_d
!========Option for deep convection=======================================

    real, dimension(size(tb,3)) :: am1, am2, am3, am4, am5, qntmp
    
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: pmass    ! layer mass (kg/m2)
    real, dimension(size(tb,1),size(tb,2))            :: tempdiag ! temporary diagnostic variable
    real, dimension(size(tracers,1),size(tracers,2),size(tracers,3),size(tracers,4))  :: trwet 
    ! calculated tracer wet deposition tendencies

    integer imax, jmax, kmax
    integer kd, ntracers
    integer ktop_tmp, kbot_tmp
    real :: tten_intg, qvten_intg
    real, dimension(size(tb,3)) :: tten_tmp, qvten_tmp
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tten_rad
    
    logical used
    type(sounding)          :: sd, sd1
    type(adicloud)          :: ac, ac1
    type(cclosure)          :: cc, cc1
    type(cplume)            :: cp, cp1
    type(ctend)             :: ct, ct1
    type(cpnlist)           :: cpn,dpn
    type(deepc)             :: dpc
    integer ::  ier
    character(len=256) :: ermesg

    kd = size(tracers,3)
    ntracers = size(tracers,4)
    call sd_init_k(kd,ntracers,sd);
    call sd_init_k(kd,ntracers,sd1);
    call ac_init_k(kd,ac);
    call ac_init_k(kd,ac1);
    call cp_init_k(kd,ntracers,cp)
    call cp_init_k(kd,ntracers,cp1)
    call ct_init_k(kd,ntracers,ct)
    call ct_init_k(kd,ntracers,ct1)
    !pack namelist parameters into plume and closure structure
    cpn % do_qctflx_zero = do_qctflx_zero
    cpn % do_detran_zero = do_detran_zero
    cpn % rle       = rle
    cpn % rpen      = rpen
    cpn % rmaxfrac  = rmaxfrac
    cpn % wmin      = wmin
    cpn % rbuoy     = rbuoy
    cpn % rdrag     = rdrag  
    cpn % frac_drs  = frac_drs
    cpn % bigc      = bigc    
    cpn % auto_th0  = auto_th0
    cpn % deltaqc0  = deltaqc0
    cpn % do_pdfpcp = do_pdfpcp
    cpn % do_pmadjt = do_pmadjt
    cpn % do_emmax  = do_emmax
    cpn % do_pnqv   = do_pnqv
    cpn % emfrac_max= emfrac_max
    cpn % auto_rate = auto_rate
    cpn % tcrit     = tcrit  
    cpn % cldhgt_max= cldhgt_max
    cpn % do_ice    = do_ice
    cpn % do_ppen   = do_ppen
    cpn % do_pevap  = do_pevap
    cpn % hcevap    = hcevap
    cpn % cfrac     = cfrac
    cpn % mixing_assumption= mixing_assumption
    cpn % mp_choice = mp_choice
    cpn % Nl_land   = Nl_land
    cpn % Nl_ocean  = Nl_ocean
    cpn % qi_thresh = qi_thresh
    cpn % r_thresh  = r_thresh
    cpn % peff_l    = peff_l
    cpn % peff_i    = peff_i
    cpn % t00       = t00
    cpn % rh0       = rh0
    cpn % do_forcedlifting= do_forcedlifting
    cpn % atopevap  = atopevap
    cpn % wtwmin_ratio = wmin_ratio*wmin_ratio
    cpn % do_auto_aero = do_auto_aero
    cpn % rad_crit = rad_crit
    cpn % wrel_min = wrel_min
    cpn % do_weffect = do_weffect
    cpn % weffect    = weffect
    cpn % use_online_aerosol = use_online_aerosol
    if (ntracers > 0) then
      allocate ( cpn%tracername   (ntracers) )
      allocate ( cpn%tracer_units (ntracers) )
      allocate ( cpn%wetdep       (ntracers) )
      cpn%tracername(:) = tracername(:)
      cpn%tracer_units(:) = tracer_units(:)
      cpn%wetdep(:)%scheme = wetdep(:)%scheme
      cpn%wetdep(:)%Henry_constant = wetdep(:)%Henry_constant
      cpn%wetdep(:)%Henry_variable = wetdep(:)%Henry_variable
      cpn%wetdep(:)%frac_in_cloud = wetdep(:)%frac_in_cloud
      cpn%wetdep(:)%alpha_r = wetdep(:)%alpha_r
      cpn%wetdep(:)%alpha_s = wetdep(:)%alpha_s
      cpn%wetdep(:)%Lwetdep = wetdep(:)%Lwetdep
      cpn%wetdep(:)%Lgas = wetdep(:)%Lgas
      cpn%wetdep(:)%Laerosol = wetdep(:)%Laerosol
      cpn%wetdep(:)%Lice = wetdep(:)%Lice
      allocate ( dpn%tracername   (ntracers) )
      allocate ( dpn%tracer_units (ntracers) )
      allocate ( dpn%wetdep       (ntracers) )
      dpn%tracername(:) = tracername(:)
      dpn%tracer_units(:) = tracer_units(:)
      dpn%wetdep(:)%scheme = wetdep(:)%scheme
      dpn%wetdep(:)%Henry_constant = wetdep(:)%Henry_constant
      dpn%wetdep(:)%Henry_variable = wetdep(:)%Henry_variable
      dpn%wetdep(:)%frac_in_cloud = wetdep(:)%frac_in_cloud
      dpn%wetdep(:)%alpha_r = wetdep(:)%alpha_r
      dpn%wetdep(:)%alpha_s = wetdep(:)%alpha_s
      dpn%wetdep(:)%Lwetdep = wetdep(:)%Lwetdep
      dpn%wetdep(:)%Lgas = wetdep(:)%Lgas
      dpn%wetdep(:)%Laerosol = wetdep(:)%Laerosol
      dpn%wetdep(:)%Lice = wetdep(:)%Lice
    endif
    call cpn_copy(cpn, dpn)

    cc  % igauss    = igauss
    cc  % rkfre     = rkfre
    cc  % rmaxfrac  = rmaxfrac
    cc  % wcrit_min = wcrit_min
    cc  % rbuoy     = rbuoy
    cc  % tau_sh    = tau_sh
!========Option for deep convection=======================================
    dpc % rkm_dp1             = rkm_dp1
    dpc % rkm_dp2             = rkm_dp2
    dpc % cbmf_dp_frac1       = cbmf_dp_frac1
    dpc % cbmf_dp_frac2       = cbmf_dp_frac2
    dpc % crh_th_ocean        = crh_th_ocean
    dpc % crh_th_land         = crh_th_land
    dpc % cape_th             = cape_th
    dpc % tau_dp              = tau_dp
    dpc % mixing_assumption_d = mixing_assumption_d
    dpc % do_ppen_d           = do_ppen_d
    dpc % rpen_d              = rpen_d
    dpc % do_pevap_d          = do_pevap_d
    dpc % cfrac_d             = cfrac_d
    dpc % hcevap_d            = hcevap_d
    dpc % frac_limit_d        = frac_limit_d
    dpc % dcapedm_th          = dcapedm_th
    dpc % do_forcedlifting_d  = do_forcedlifting_d
    dpc % lofactor_d          = lofactor_d
    dpc % auto_th0_d          = auto_th0_d
    dpc % tcrit_d             = tcrit_d
!========Option for deep convection=======================================
    imax  = size( tb, 1 )
    jmax  = size( tb, 2 )
    kmax  = size( tb, 3 )
    sd % kmax=kmax

    kl=kmax-1
    klm=kl-1

   !initialize 3D variables outside the loop

    tten=0.; qvten=0.; qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    uten=0.; vten =0.; rain =0.; snow =0.; plcl =0.; plfc=0.; plnb=0.;  
    cldqa=0.; cldql=0.; cldqi=0.; cldqn=0.;
    hlflx=0.; qtflx=0.; pflx=0.; am1=0.; am2=0.; am3=0.; am4=0.;
    tten_pevap=0.; qvten_pevap=0.;
    ice_pflx = 0. ; liq_pflx = 0.

    cino=0.; capeo=0.; tkeo=0.; wrelo=0.; ufrco=0.; zinvo=0.; wuo=0.; peo=0.; 
    fero=0.; fdro=0.; fdrso=0.; cmf=0.; denth=0.;  dqtmp=0.; ocode=0;
    dcapeo=0.; dcino=0.; dwfno=0.; xhlsrc=0.; xqtsrc=0.; fdp=0.;
    trtend=0.; qldet=0.; qidet=0.; qadet=0.; crho=0.; hmo=0.; hms=0.; abu=0.;
    trwet = 0.
    dting = 0.

    naer = size(asol%aerosol,4)

!========Option for deep convection=======================================
    if (do_deep) then
       tten_d=0.; qvten_d=0.; qlten_d=0.; qiten_d=0.; qaten_d=0.; qnten_d=0.;
       uten_d=0.; vten_d =0.; rain_d =0.; snow_d =0.; qtten_d=0.;
       cldqa_d=0.; cldql_d=0.; cldqi_d=0.; cldqn_d=0.;
       hlflx_d=0.; qtflx_d=0.; pflx_d=0.;
       wuo_d=0.; fero_d=0.; fdro_d=0.; fdrso_d=0.; 
       cmf_d=0.; cbu_d=0.;
       denth_d=0.; dting_d=0.; dqtmp_d=0.; cbmf_d=0.; dcapedm_d=0.; dcino=0.; dwfn_d=0.;
       tten_pevap_d=0.; qvten_pevap_d=0.; rkm_d=0.;
    end if
!========Option for deep convection=======================================

    do j = 1, jmax
       do i=1, imax

         do k=1,kmax
           pmass(i,j,k) = (pint(i,j,k+1) - pint(i,j,k))/GRAV
         enddo
    !relaxation TKE back to 0 with time-scale of disscale
    !tkeavg = ustar(i,j)*bstar(i,j)*disscale 
    !dissipate tke with length-scale of disscale
    !tkeavg=(ustar(i,j)*bstar(i,j)*disscale)**(2./3.)
    !below following Holtslag and Boville 1993

         if (pblht(i,j).lt.0.) then
           temp_1=0.0
         elseif (pblht(i,j).gt.5000.) then
           temp_1=5000.
         else
           temp_1=pblht(i,j)
         endif
         temp_1=ustar(i,j)**3.+0.6*ustar(i,j)*bstar(i,j)*temp_1
         if (temp_1 .gt. 0.) temp_1 = 0.5*temp_1**(2./3.)
         tkeo(i,j) = MAX (tkemin, temp_1)

         if (do_gust_cv) then
           if (cbmfo(i,j)>0) then
             tkeo(i,j) = tkeo(i,j)+(gustmax*sqrt(cbmfo(i,j)/(gustconst + cbmfo(i,j))))**2.
           endif
         endif

         if (skip_calculation(i,j)) then
           ocode(i,j) = 6
           go to 100
         endif
         call clearit(ac, cc, cp, ct, cp1, ct1);

! restrict grid-box area available to shallow convection to that which 
! is not involved with deep convection
          cp%maxcldfrac = minval(max_available_cf(i,j,:))
          cc%maxcldfrac = cp%maxcldfrac

          cc%scaleh = cush(i,j); 
          cush(i,j) = -1.;
          if(cc%scaleh.le.0.0) cc%scaleh=1000.

          am1(:) = 0.; am2(:) = 0.; am3(:) = 0.; am4(:) = 0.; am5(:) = 0.;

          do k=1,kmax
            tmp=1. / (zint(i,j,k)-zint(i,j,k+1)) * 1.0e9 * 1.0e-12
            if(use_online_aerosol) then
              do na = 1,naer
                if(asol%aerosol_names(na) == 'so4' .or. &
                   asol%aerosol_names(na) == 'so4_anthro' .or. &
                   asol%aerosol_names(na) == 'so4_natural') then
                           am1(k)=am1(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'omphilic' .or. &
                        asol%aerosol_names(na) == 'omphobic') then
                           am4(k)=am4(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'bcphilic' .or. &
                        asol%aerosol_names(na) == 'bcphobic' .or. &
                        asol%aerosol_names(na) == 'dust1' .or. &
                        asol%aerosol_names(na) == 'dust2' .or. &
                        asol%aerosol_names(na) == 'dust3' ) then
                           am2(k)=am2(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'seasalt1' .or. &
                        asol%aerosol_names(na) == 'seasalt2') then
                           am3(k)=am3(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'seasalt3' .or. &
                        asol%aerosol_names(na) == 'seasalt4' .or. &
                        asol%aerosol_names(na) == 'seasalt5' ) then
                           am5(k)=am5(k)+asol%aerosol(i,j,k,na)*tmp
                end if
              end do
              am2(k)=am2(k)+am3(k)+am4(k)
              if(.not. use_sub_seasalt) am3(k)=am3(k)+am5(k)
            else
              am1(k)= asol%aerosol(i,j,k,2)*tmp
              am2(k)= asol%aerosol(i,j,k,1)*tmp
              am3(k)= sea_salt_scale*asol%aerosol(i,j,k,5)*tmp
              am4(k)= om_to_oc*asol%aerosol(i,j,k,3)*tmp
            endif
          end do

!========Pack column properties into a sounding structure====================

          if (do_qn) then
             qntmp(:)=q(i,j,:,nqn)
          else
             qntmp(:)=0.
          end if
          call pack_sd_k(land(i,j), coldT(i,j), delt, pmid(i,j,:), pint(i,j,:),     &
               zmid(i,j,:), zint(i,j,:), ub(i,j,:), vb(i,j,:), tb(i,j,:),   &
               q(i,j,:,nqv), q(i,j,:,nql), q(i,j,:,nqi), q(i,j,:,nqa), qntmp,       &
               am1(:), am2(:), am3(:), am4(:), tracers(i,j,:,:), sd, Uw_p)

!========Finite volume intepolation==========================================

          call extend_sd_k(sd,  pblht(i,j),do_ice, Uw_p)
          sd%tke = tkeo(i,j)
          zinvo (i,j) = sd%zinv


!========Find source air, and do adiabatic cloud lifting======================

          zsrc  =sd%zs (1)
          psrc  =sd%ps (1)
          thcsrc=sd%thc(1)
          qctsrc=sd%qct(1)
          hlsrc =sd%hl (1)
          rkm_sh1=rkm_sh
          if (do_lands) then
            !wstar   = (ustar(i,j)*bstar(i,j)*pblht(i,j))**(1./3.)
             cpn % auto_th0 = auto_th0 * (1. + landfact_m * sd%land)
             call qt_parcel_k (sd%qs(1), qstar(i,j), pblht(i,j), sd%tke, sd%land, gama, &
                  pblht0, tke0, lofactor0, lochoice, qctsrc, lofactor)
             rkm_sh1 = rkm_sh   * lofactor
          end if

          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, do_fast, do_ice, ac)
          ac % usrc = sd%u(sd%ktoppbl)
          ac % vsrc = sd%v(sd%ktoppbl)
          if (ac%plfc.eq.0) ac%plfc=psrc
          if (ac%plnb.eq.0) ac%plnb=psrc
          plcl (i,j) = ac%plcl
          plfc (i,j) = ac%plfc
          plnb (i,j) = ac%plnb
          cino (i,j) = ac%cin
          capeo(i,j) = ac%cape
          xhlsrc(i,j)= ac%hlsrc; 
          xqtsrc(i,j)= ac%qctsrc; 
          crho(i,j)  = sd%crh;
          do k = 1,kmax
             nk = kmax+1-k
             hmo  (i,j,nk) = sd%hm(k);
             hms  (i,j,nk) = sd%hms(k);
             abu  (i,j,nk) = ac%buo(k);
          end do

          if (do_fast) then
!             if (ac%klcl.eq.0 .or. ac%plcl.gt.sd%ps(1) .or. ac%plcl.lt.20000.) then
             if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
                ocode(i,j)=1; goto 100 !cycle;
             end if
             if (ac%plfc.lt.500.) then
                ocode(i,j)=2; goto 100 !cycle;
             end if
          end if

!========Cumulus closure to determine cloud base mass flux===================

          cbmf_old=cbmfo(i,j); cc%cbmf=cbmf_old;

          if (iclosure.eq.0) then
             call cclosure_bretherton(sd%tke, cpn, sd, Uw_p, ac, cc)
          else if (iclosure.eq.1) then
             call cclosure_implicit(sd%tke, cpn, sd, Uw_p, ac, cc, delt, rkm_sh1, &
                  do_coldT, sd1, ac1, cc1, cp, ct, ier, ermesg) 
             if (ier /= 0) then
               call error_mesg ('subroutine uw_conv iclosure=1 ', ermesg, FATAL)
             endif
          else if (iclosure.eq.2) then
             call cclosure_relaxwfn(sd%tke, cpn, sd, Uw_p, ac, cc, cp, ct, delt,  &
                  rkm_sh1, do_coldT, sd1, ac1, cc1, cp1, ct1, ier, ermesg)
             if (ier /= 0) then
               call error_mesg ('subroutine uw_conv iclosure=2 ', ermesg, FATAL)
             endif
          end if

          cbmfo(i,j) = cc%cbmf
          wrelo(i,j) = cc%wrel
          ufrco(i,j) = cc%ufrc

          if (.not.do_fast) then
             if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
!             if (ac%klcl.eq.0 .or. ac%plcl.lt.20000.) then
                ocode(i,j)=1; goto 100 !cycle;
             end if
             if (ac%plfc.lt.500.) then
                ocode(i,j)=2; goto 100 !cycle;
             end if
          end if

          if(cc%cbmf.lt.1.e-6 .or. cc%wrel.eq.0.) then
             ocode(i,j)=3; goto 100 !cycle;
          end if

!========Do shallow cumulus plume calculation================================

          cpn%isdeep=.false.


          cbmf_shallow = cc%cbmf ! - cbmf_deep
          cpn%do_ppen=do_ppen
          cpn%rpen   =rpen
          call cumulus_plume_k(cpn, sd, ac, cp, rkm_sh1, cbmf_shallow, cc%wrel, cc%scaleh, Uw_p, ier, ermesg)
          if (ier /= 0) then
            call error_mesg ('subroutine uw_conv', ermesg, FATAL)
          endif
          if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
             ocode(i,j)=4; goto 100 !cycle;
          end if
          if(cp%cldhgt.ge.cldhgt_max) then
             ocode(i,j)=5; goto 100 !cycle;
          end if

          if (cpn%isdeep .EQV. .true.) then 
             fdp(i,j) = 1
          else
             fdp(i,j) = 0
          end if

          cush(i,j)=cp%cush


!========Calculate cumulus produced tendencies===============================

          call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)

!========Unpack convective tendencies========================================
          do k = 1,cp%ltop
             nk = kmax+1-k
             uten  (i,j,nk) = ct%uten (k)
             vten  (i,j,nk) = ct%vten (k)
             qlten (i,j,nk) = ct%qlten(k)
             qiten (i,j,nk) = ct%qiten(k)
             qaten (i,j,nk) = ct%qaten(k)
             qnten (i,j,nk) = ct%qnten(k)
             qldet (i,j,nk) = ct%qldet(k)
             qidet (i,j,nk) = ct%qidet(k)
             qadet (i,j,nk) = ct%qadet(k)
             qvten (i,j,nk) = ct%qvten(k)
             pflx  (i,j,nk) = ct%pflx (k)
             ice_pflx(i,j,nk) = cp%ppti(k)
             liq_pflx(i,j,nk) = cp%pptr(k)
             tten  (i,j,nk) = ct%tten (k)
             rhos0j = sd%ps(k)/(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*sd%exners(k))
             hlflx(i,j,nk) = ct%hlflx(k)
             qtflx (i,j,nk) = rhos0j*HLv*ct%qctflx(k)
             tten_pevap (i,j,nk) = ct%tevap (k)
             qvten_pevap(i,j,nk) = ct%qevap (k)
             
             cldqa (i,j,nk) = cp%ufrc(k)
             cldql (i,j,nk) = cp%qlu(k)
             cldqi (i,j,nk) = cp%qiu(k)
             cldqn (i,j,nk) = cp%qnu(k)
             cmf   (i,j,nk) = cp%umf(k)
             wuo   (i,j,nk) = cp%wu (k)
             peo   (i,j,nk) = cp%peff(k)
             fero  (i,j,nk) = cp%fer(k)
             fdro  (i,j,nk) = cp%fdr(k)
             fdrso (i,j,nk) = cp%fdrsat(k)*cp%umf(k)
          enddo

! make sure the predicted tracer tendencies do not produce negative
! tracers due to convective tendencies. if necessary, adjust the 
! tendencies.
          call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                           cp%tr, ct%trten, ct%trwet) 
          do k = 1,cp%ltop
             nk = kmax+1-k
             do n = 1, size(trtend,4)
               trtend(i,j,nk,n) = ct%trten(k,n) + ct%trwet(k,n)
               trwet(i,j,nk,n)  = ct%trwet(k,n)
             enddo
          enddo
          snow  (i,j)  = ct%snow
          rain  (i,j)  = ct%rain
          denth (i,j)  = ct%denth
          dqtmp (i,j)  = ct%dqtmp
          dting (i,j)  = ct%dting

!========Option for deep convection=======================================
100       if (do_deep) then
	     cbmf_deep = 0.
	     rkm_dp = 0.
             tmp   = max(min (sd%crh, 1.0), 0.0)
	     crh_th  = sd%land*dpc%crh_th_land+(1.-sd%land)*dpc%crh_th_ocean
	     dcrh  = tmp - crh_th
             dcrh0 = 1.0001-crh_th
	     if (dcrh .gt. 0) then
	        dcrh = dcrh/dcrh0;
	        dcrh = dcrh**norder
	        rkm_dp       = dpc%rkm_dp1      + dcrh * (dpc%rkm_dp2      -dpc%rkm_dp1)
	        cbmf_dp_frac = dpc%cbmf_dp_frac1+ dcrh * (dpc%cbmf_dp_frac2-dpc%cbmf_dp_frac1)
	        cbmf_deep    = 1. !%cbmf_dp_frac * cc%cbmf
                lofactor     = 1. - sd%land * (1. - dpc%lofactor_d)
	        if (do_lod_rkm) then
               	   rkm_dp       = rkm_dp * lofactor
	        elseif (do_lod_cfrac) then
		   dpc % cfrac_d= dpc % cfrac_d * lofactor
	        elseif (do_lod_tcrit) then
             	   dpc % tcrit_d= dpc % tcrit_d * lofactor
 	        end if
	     end if
             rkm_d(i,j) = rkm_dp;

             dpn % do_ppen  = dpc % do_ppen_d
             dpn % do_pevap = dpc % do_pevap_d
             dpn % cfrac    = dpc % cfrac_d
             dpn % hcevap   = dpc % hcevap_d
             dpn % tcrit    = dpc % tcrit_d
             dpn % auto_th0 = dpc % auto_th0_d
             dpn % mixing_assumption = dpc % mixing_assumption_d
             dpn % do_forcedlifting  = dpc % do_forcedlifting_d
             if (idpchoice.eq.0) then
                call  dpconv0(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, cp1, ct1, ocode(i,j), ier, ermesg)
             else if (idpchoice.eq.1) then
                call  dpconv1(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j), dcapedm_d(i,j), ier, ermesg)
             else if (idpchoice.eq.2) then
                call  dpconv2(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j), ier, ermesg)
             end if
             if (ier /= 0) then
                call error_mesg ('uw_conv calling dpconv', ermesg, FATAL)
             endif
             if(ocode(i,j).ge.6) cycle;
             do k = 1,kmax !cp1%ltop
                nk = kmax+1-k
                uten_d  (i,j,nk) = ct1%uten (k)
                vten_d  (i,j,nk) = ct1%vten (k)
                qlten_d (i,j,nk) = ct1%qlten(k)
                qiten_d (i,j,nk) = ct1%qiten(k)
                qaten_d (i,j,nk) = ct1%qaten(k) 
                qnten_d (i,j,nk) = ct1%qnten(k) 
                qvten_d (i,j,nk) = ct1%qvten(k)
                qtten_d (i,j,nk) = ct1%qctten(k)
                pflx_d  (i,j,nk) = ct1%pflx (k)
                tten_d  (i,j,nk) = ct1%tten (k)
                hlflx_d (i,j,nk) = ct1%hlflx(k) 
                qtflx_d (i,j,nk) = ct1%qctflx(k)
                cldqa_d (i,j,nk) = cp1%ufrc(k)
                cldql_d (i,j,nk) = cp1%qlu(k)
                cldqi_d (i,j,nk) = cp1%qiu(k)
                cldqn_d (i,j,nk) = cp1%qnu(k)
                cmf_d   (i,j,nk) = cp1%umf(k) + cp1%emf(k)
                cbu_d   (i,j,nk) = cp1%buo(k)
                tten_pevap_d (i,j,nk) = ct1%tevap (k)
                qvten_pevap_d(i,j,nk) = ct1%qevap (k)
                wuo_d   (i,j,nk) = cp1%wu (k)
                fero_d  (i,j,nk) = cp1%fer(k)
                fdro_d  (i,j,nk) = cp1%fdr(k) 
                fdrso_d (i,j,nk) = cp1%fdrsat(k)*cp1%fdr(k)*cp1%umf(k)
             enddo
             snow_d  (i,j)  = ct1%snow
             rain_d  (i,j)  = ct1%rain
             cbmf_d  (i,j)  = cbmf_deep
             denth_d (i,j)  = ct1%denth
             dting_d (i,j)  = ct1%dting
             dqtmp_d (i,j)  = ct1%dqtmp
             !dwfn_d  (i,j)  = cc%dwfn

!========Option for deep convection=======================================
            
             uten  (i,j,:) = uten  (i,j,:) + uten_d  (i,j,:)
             vten  (i,j,:) = vten  (i,j,:) + vten_d  (i,j,:)
             qlten (i,j,:) = qlten (i,j,:) + qlten_d (i,j,:)
             qiten (i,j,:) = qiten (i,j,:) + qiten_d (i,j,:)
             qaten (i,j,:) = qaten (i,j,:) + qaten_d (i,j,:) 
             qnten (i,j,:) = qnten (i,j,:) + qnten_d (i,j,:) 
             qvten (i,j,:) = qvten (i,j,:) + qvten_d (i,j,:)
             pflx  (i,j,:) = pflx  (i,j,:) + pflx_d  (i,j,:)
             tten  (i,j,:) = tten  (i,j,:) + tten_d  (i,j,:)
             hlflx (i,j,:) = hlflx (i,j,:) + hlflx_d (i,j,:) 
             qtflx (i,j,:) = qtflx (i,j,:) + qtflx_d (i,j,:)
             cmf   (i,j,:) = cmf   (i,j,:) + cmf_d   (i,j,:)
             tten_pevap (i,j,:)=tten_pevap (i,j,:) + tten_pevap_d (i,j,:) 
             qvten_pevap(i,j,:)=qvten_pevap(i,j,:) + qvten_pevap_d(i,j,:) 
             !wuo   (i,j,:) = wuo   (i,j,:) 
             !fero  (i,j,:) = fero  (i,j,:) 
             !fdro  (i,j,:) = fdro  (i,j,:)  
             !fdrso (i,j,:) = fdrso (i,j,:) 
             !do n = 1, size(trtend,4)
             !   trtend(i,j,:,n) = trtend(i,j,:,n) 
             !enddo
             snow  (i,j)  = snow  (i,j) + snow_d  (i,j)
             rain  (i,j)  = rain  (i,j) + rain_d  (i,j)
             denth (i,j)  = denth (i,j) + denth_d (i,j)
             dting (i,j)  = dting (i,j) + dting_d (i,j)
             dqtmp (i,j)  = dqtmp (i,j) + dqtmp_d (i,j)
             !cbmfo (i,j)  = cc%cbmf
             !dwfno (i,j)  = cc%dwfn
          end if
!========Option for deep convection=======================================
	if (do_no_uw_conv) then
          uten (i,j,:)=0.; vten (i,j,:)=0.; 
          tten (i,j,:)=0.; qvten(i,j,:)=0.; 
          cmf  (i,j,:)=0.; qlten(i,j,:)=0.;
          qiten(i,j,:)=0.; qaten(i,j,:)=0.; 
          qnten(i,j,:)=0.; rain (i,j)  =0.; snow(i,j)=0.;
	endif
	if (do_imposing_rad_cooling) then
           tten_rad (i,j,:)=0;
	   do k = 1,sd%kmax
            nk = kmax+1-k
            if (sd%t(k)>t_thresh) then
               tten_rad (i,j,nk) = cooling_rate/86400.
	    else
               tten_rad (i,j,nk) = (t_strato-sd%t(k))/(tau_rad*86400.)
            end if
           enddo
	end if
	if (do_imposing_cooling_drying) then
           kbot_tmp=1; 
           ktop_tmp=sd%kmax;
	   do k=1,sd%kmax
	      if (sd%p(k)>=7500) then
	      	 ktop_tmp=k
	      end if
	      if (sd%p(k)>=85000) then
	      	 kbot_tmp=k
	      end if
	   enddo
           tten_tmp (:)=0; 
           qvten_tmp(:)=0;
	   do k = kbot_tmp,ktop_tmp
            if (sd%p(k)>pres_min .and. sd%p(k)<=pres_max) then
               tten_tmp (k)=tdt_rate/86400.
	       qvten_tmp(k)=qvdt_rate/86400.
            end if
           enddo
!           tten_intg=0.
!           qvten_intg=0.
!           dpsum =0.
!           do k = kbot_tmp,ktop_tmp
!     	      tten_intg  = tten_intg  + tten_tmp (k)*sd%dp(k)
!	      qvten_intg = qvten_intg + qvten_tmp(k)*sd%dp(k)
!	      dpsum    = dpsum + sd%dp(k)
!	    end do
!           do k = kbot_tmp,ktop_tmp
!             tten_tmp(k) = tten_tmp(k)  - tten_intg / dpsum
!	      qvten_tmp(k)= qvten_tmp(k) - qvten_intg /dpsum
!           end do
           do k = 1,kmax
              nk = kmax+1-k
              tten_rad  (i,j,nk) = tten_rad (i,j,nk) + tten_tmp (k)
           end do
	end if

       enddo
    enddo

    call sd_end_k(sd)
    call sd_end_k(sd1)
    call ac_end_k(ac)
    call ac_end_k(ac1)
    call cp_end_k(cp)
    call cp_end_k(cp1)
    call ct_end_k(ct)
    call ct_end_k(ct1)
    if (_ALLOCATED ( cpn%tracername    ))  deallocate ( cpn%tracername    )
    if (_ALLOCATED ( cpn%tracer_units  ))  deallocate ( cpn%tracer_units  )
    if (_ALLOCATED ( cpn%wetdep        ))  deallocate ( cpn%wetdep        )
    if (_ALLOCATED ( dpn%tracername    ))  deallocate ( dpn%tracername    )
    if (_ALLOCATED ( dpn%tracer_units  ))  deallocate ( dpn%tracer_units  )
    if (_ALLOCATED ( dpn%wetdep        ))  deallocate ( dpn%wetdep        )
    if (.not.do_uwcmt) then
       uten=0.;
       vten=0.;
    end if

    if ( prevent_unreasonable ) then
      scale_uw=HUGE(1.0)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            if ((q(i,j,k,nqa) + qaten(i,j,k)*delt) .lt. 0. .and. (qaten(i,j,k).ne.0.)) then
              qaten(i,j,k) = -1.*q(i,j,k,nqa)/delt
            end if
            if ((q(i,j,k,nqa) + qaten(i,j,k)*delt) .gt. 1. .and. (qaten(i,j,k).ne.0.)) then
              qaten(i,j,k)= (1. - q(i,j,k,nqa))/delt
            end if
 
            if ((q(i,j,k,nql) + qlten(i,j,k)*delt) .lt. 0. .and. (qlten(i,j,k).ne.0.)) then
              tten (i,j,k) = tten(i,j,k) -(q(i,j,k,nql)/delt+qlten(i,j,k))*HLv/Cp_Air
              qvten(i,j,k) = qvten(i,j,k)+(q(i,j,k,nql)/delt+qlten(i,j,k))
              qlten(i,j,k) = qlten(i,j,k)-(q(i,j,k,nql)/delt+qlten(i,j,k))
            end if
 
            if ((q(i,j,k,nqi) + qiten(i,j,k)*delt) .lt. 0. .and. (qiten(i,j,k).ne.0.)) then
              tten (i,j,k) = tten(i,j,k) -(q(i,j,k,nqi)/delt+qiten(i,j,k))*HLs/Cp_Air
              qvten(i,j,k) = qvten(i,j,k)+(q(i,j,k,nqi)/delt+qiten(i,j,k))
    
              qiten(i,j,k) = qiten(i,j,k)-(q(i,j,k,nqi)/delt+qiten(i,j,k))
            end if

            if (do_qn) then
              if ((q(i,j,k,nqn) + qnten(i,j,k)*delt) .lt. 0. .and. (qnten(i,j,k).ne.0.)) then
                qnten(i,j,k) = qnten(i,j,k)-(q(i,j,k,nqn)/delt+qnten(i,j,k))
              end if
            endif
    !rescaling to prevent negative specific humidity for each grid point
            if (do_rescale) then
              qtin =  q(i,j,k,nqv)
              dqt  =  qvten(i,j,k) * delt
              if ( dqt.lt.0 .and. qtin+dqt.lt.1.e-10 ) then
                temp_1 = max( 0.0, -(qtin-1.e-10)/dqt )
              else
                temp_1 = 1.0
              endif
    !scaling factor for each column is the minimum value within that column
              scale_uw(i,j) = min( temp_1, scale_uw(i,j))
            endif
          enddo
        enddo
      enddo

!     where ((tracers(:,:,:,:) + trtend(:,:,:,:)*delt) .lt. 0.)
!        trtend(:,:,:,:) = -tracers(:,:,:,:)/delt
!     end where

      if (do_rescale) then
      !scale tendencies
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              uten (i,j,k)  = scale_uw(i,j) * uten (i,j,k)
              vten (i,j,k)  = scale_uw(i,j) * vten (i,j,k)
              tten (i,j,k)  = scale_uw(i,j) * tten (i,j,k)
              qvten(i,j,k)  = scale_uw(i,j) * qvten(i,j,k)
              qlten(i,j,k)  = scale_uw(i,j) * qlten(i,j,k)
              qiten(i,j,k)  = scale_uw(i,j) * qiten(i,j,k)
              qaten(i,j,k)  = scale_uw(i,j) * qaten(i,j,k)
              if (do_qn) qnten(i,j,k) = scale_uw(i,j) * qnten(i,j,k)
              if (k.eq.kmax) then
                rain(i,j) = scale_uw(i,j) * rain(i,j)
                snow(i,j) = scale_uw(i,j) * snow(i,j)
              endif
            end do
          end do
        end do
      end if
    endif


    do k=1,kmax
      do j=1,jmax
        do i=1,imax
          cfq(i,j,k) = 0
          if (wuo(i,j,k) .gt. 0.) then
            cfq(i,j,k) = 1
          endif
        enddo
      enddo
    enddo

    !diagnostic output
    used = send_data( id_xhlsrc_uwc,       xhlsrc,             Time, is, js)
    used = send_data( id_xqtsrc_uwc,       xqtsrc,             Time, is, js)
    used = send_data( id_tdt_pevap_uwc,    tten_pevap*aday , Time, is, js, 1)
    used = send_data( id_qdt_pevap_uwc,    qvten_pevap*aday, Time, is, js, 1)

    used = send_data( id_tdt_uwc,    tten*aday , Time, is, js, 1)
    used = send_data( id_qdt_uwc,    qvten*aday, Time, is, js, 1)
    used = send_data( id_cmf_uwc,    cmf,          Time, is, js, 1)
    used = send_data( id_cfq_uwc,    cfq,          Time, is, js, 1)
    used = send_data( id_wu_uwc,     wuo,          Time, is, js, 1)
    used = send_data( id_peo_uwc,    peo,          Time, is, js, 1)
    used = send_data( id_fer_uwc,    fero,         Time, is, js, 1)
    used = send_data( id_fdr_uwc,    fdro,         Time, is, js, 1)
    used = send_data( id_fdrs_uwc,   fdrso,        Time, is, js, 1)
    used = send_data( id_cqa_uwc,    cldqa,        Time, is, js, 1)
    used = send_data( id_cql_uwc,    cldql,        Time, is, js, 1)
    used = send_data( id_cqi_uwc,    cldqi,        Time, is, js, 1)
    used = send_data( id_cqn_uwc,    cldqn,        Time, is, js, 1)
    used = send_data( id_hlflx_uwc, hlflx,       Time, is, js, 1)
    used = send_data( id_qtflx_uwc,  qtflx,        Time, is, js, 1)
    used = send_data( id_hmo_uwc,    hmo,          Time, is, js, 1)
    used = send_data( id_hms_uwc,    hms,          Time, is, js, 1)
    used = send_data( id_abu_uwc,    abu,          Time, is, js, 1)
    used = send_data( id_tdt_rad_uwc,tten_rad*aday,Time, is, js, 1)!miz
  
    used = send_data( id_prec_uwc, (rain+snow)*aday, Time, is, js )
    used = send_data( id_snow_uwc, (snow)*aday,      Time, is, js )
    used = send_data( id_cin_uwc,  (cino),             Time, is, js )
    used = send_data( id_cape_uwc, (capeo),            Time, is, js )
    used = send_data( id_crh_uwc,  (crho),             Time, is, js )
    used = send_data( id_tke_uwc,  (tkeo),             Time, is, js )
    used = send_data( id_cbmf_uwc, (cbmfo),            Time, is, js )
    used = send_data( id_wrel_uwc, (wrelo),            Time, is, js )
    used = send_data( id_ufrc_uwc, (ufrco),            Time, is, js )
    used = send_data( id_plcl_uwc, (plcl*0.01),        Time, is, js )
    used = send_data( id_plfc_uwc, (plfc*0.01),        Time, is, js )
    used = send_data( id_plnb_uwc, (plnb*0.01),        Time, is, js )
    used = send_data( id_zinv_uwc, (zinvo),            Time, is, js )
    used = send_data( id_cush_uwc, (cush),             Time, is, js )
    used = send_data( id_dcin_uwc, (dcino),            Time, is, js )
    used = send_data( id_dcape_uwc,(dcapeo),           Time, is, js )
    used = send_data( id_dwfn_uwc, (dwfno),            Time, is, js )
    used = send_data( id_enth_uwc, (denth),            Time, is, js )
    used = send_data( id_qtmp_uwc, (dqtmp),            Time, is, js )
    used = send_data( id_dting_uwc,(dting),            Time, is, js )
    used = send_data( id_ocode_uwc,(ocode),            Time, is, js )
    used = send_data( id_fdp_uwc,  (fdp),              Time, is, js )
   
    if ( do_strat ) then
       used = send_data( id_qldt_uwc, qlten*aday,    Time, is, js, 1)
       used = send_data( id_qidt_uwc, qiten*aday,    Time, is, js, 1)
       used = send_data( id_qadt_uwc, qaten*aday,    Time, is, js, 1)
       used = send_data( id_qndt_uwc, qnten*aday,    Time, is, js, 1)
       used = send_data( id_qldet_uwc,  qldet*aday,  Time, is, js, 1)
       used = send_data( id_qidet_uwc,  qidet*aday,  Time, is, js, 1)
       used = send_data( id_qadet_uwc,  qadet*aday,  Time, is, js, 1)
       used = send_data( id_qtdt_uwc,(qvten+qlten+qiten)*aday,Time, is, js, 1)
    end if

    if ( allocated(id_tracerdt_uwc) ) then
       do n = 1,size(id_tracerdt_uwc)
          used = send_data( id_tracerdt_uwc(n), trtend(:,:,:,n), Time, is, js, 1)
       end do
    end if
    if ( allocated(id_tracerdt_uwc_col) ) then
       do n = 1,size(id_tracerdt_uwc_col)
          if ( id_tracerdt_uwc_col(n) > 0 ) then
            tempdiag = 0.
            do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trtend(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdt_uwc_col(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if
    if ( allocated(id_tracerdtwet_uwc) ) then
       do n = 1,size(id_tracerdtwet_uwc)
          used = send_data( id_tracerdtwet_uwc(n), trwet(:,:,:,n), Time, is, js, 1)
       end do
    end if
    if ( allocated(id_tracerdtwet_uwc_col) ) then
       uw_wetdep = 0.
       do n = 1,size(id_tracerdtwet_uwc_col)
          if ( id_tracerdtwet_uwc_col(n) > 0 ) then
             tempdiag = 0.
             do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trwet(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdtwet_uwc_col(n), tempdiag(:,:), Time, is, js)
            uw_wetdep(:,:,n) = tempdiag(:,:)
          end if
       end do
    end if

!========Option for deep convection=======================================
    if (do_deep) then
       used=send_data( id_tdt_pevap_uwd,    tten_pevap_d*aday , Time, is, js, 1)
       used=send_data( id_qdt_pevap_uwd,    qvten_pevap_d*aday, Time, is, js, 1)
       used=send_data( id_tdt_uwd,   tten_d*aday , Time, is, js, 1)
       used=send_data( id_qdt_uwd,   qvten_d*aday, Time, is, js, 1)
       used=send_data( id_qtdt_uwd,  qtten_d*aday, Time, is, js, 1)
       used=send_data( id_cmf_uwd,   cmf_d,          Time, is, js, 1)
       used=send_data( id_cbu_uwd,   cbu_d,          Time, is, js, 1)
       used=send_data( id_wu_uwd,    wuo_d,          Time, is, js, 1)
       used=send_data( id_fer_uwd,   fero_d,         Time, is, js, 1)
       used=send_data( id_fdr_uwd,   fdro_d,         Time, is, js, 1)
       used=send_data( id_fdrs_uwd,  fdrso_d,        Time, is, js, 1)
       used=send_data( id_cqa_uwd,   cldqa_d,        Time, is, js, 1)
       used=send_data( id_cql_uwd,   cldql_d,        Time, is, js, 1)
       used=send_data( id_cqi_uwd,   cldqi_d,        Time, is, js, 1)
       used=send_data( id_cqn_uwd,   cldqn_d,        Time, is, js, 1)
       used=send_data( id_hlflx_uwd, hlflx_d,        Time, is, js, 1)
       used=send_data( id_qtflx_uwd, qtflx_d,        Time, is, js, 1)
      
       used=send_data( id_prec_uwd, (rain_d+snow_d)*aday,Time, is, js )
       used=send_data( id_snow_uwd, (snow_d)*aday,       Time, is, js )
       used=send_data( id_cbmf_uwd, (cbmf_d),              Time, is, js )
       used=send_data( id_dcapedm_uwd,(dcapedm_d),         Time, is, js )
       used=send_data( id_dwfn_uwd, (dwfn_d),              Time, is, js )
       used=send_data( id_enth_uwd, (denth_d),             Time, is, js )
       used=send_data( id_rkm_uwd,  (rkm_d),               Time, is, js )
             
       if ( do_strat ) then
          used=send_data( id_qldt_uwd, qlten_d*aday,     Time, is, js, 1)
          used=send_data( id_qidt_uwd, qiten_d*aday,     Time, is, js, 1)
          used=send_data( id_qadt_uwd, qaten_d*aday,     Time, is, js, 1)
       end if
    end if
!========Option for deep convection=======================================


    if (.not.apply_tendency) then
       uten=0.; vten=0.; tten=0.; qvten=0.; cmf=0.; rain=0.; snow=0.;
       qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    end if
!miz
    if (do_imposing_rad_cooling) then
       	  do j = 1, jmax	   
       	     do i=1, imax	   
	     	tten  (i,j,:) =	tten (i,j,:) + tten_rad (i,j,:)
	     end do
	  end do
    end if	
!miz
    if (do_gust_cv) then
     	  do j = 1, jmax	   
       	     do i=1, imax	   
               cbmfo(i,j) = rain(i,j) + snow(i,j)
             end do
          end do
    end if

  END SUBROUTINE UW_CONV

!#####################################################################
!#####################################################################

  subroutine clearit(ac, cc, cp, ct, cp1, ct1)

    type(adicloud), intent(inout) :: ac
    type(cclosure), intent(inout) :: cc
    type(cplume),   intent(inout) :: cp,cp1
    type(ctend),    intent(inout) :: ct,ct1

    call ac_clear_k(ac); 
    ac%klcl =0;  ac%klfc =0;  ac%klnb =0; 

    cc%wrel=0.; cc%ufrc=0.; cc%scaleh=0.;

    call cp_clear_k(cp);  cp%maxcldfrac =1.;
    call ct_clear_k(ct);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);

  end subroutine clearit


!#####################################################################

end MODULE UW_CONV_MOD
