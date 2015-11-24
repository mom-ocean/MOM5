
                    module moist_processes_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             moist convective adjustment
!             relaxed arakawa-schubert
!             donner deep convection
!             large-scale condensation
!             stratiform prognostic cloud scheme 
!             rel humidity cloud scheme 
!             diagnostic cloud scheme 
!             lin cloud microphysics
!             betts-miller convective adjustment
!
!-----------------------------------------------------------------------

! fms modules
use sat_vapor_pres_mod,    only: compute_qs, lookup_es
use time_manager_mod,      only: time_type, get_time
use diag_manager_mod,      only: register_diag_field, send_data
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 file_exist, check_nml_error,    &
                                 open_namelist_file, close_file, &
                                 write_version_number,           &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 MPP_CLOCK_SYNC, read_data, write_data
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index,&
                                 get_number_tracers, &
                                 get_tracer_names, &
                                 query_method, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                 RDGAS, RVGAS, TFREEZE, WTMAIR, &
                                 SECONDS_PER_DAY, KAPPA
! atmos_param modules
use betts_miller_mod,      only: betts_miller, betts_miller_init
use bm_massflux_mod,       only: bm_massflux, bm_massflux_init
use bm_omp_mod,            only: bm_omp, bm_omp_init
use donner_deep_mod,       only: donner_deep_init,               &
                                 donner_deep_time_vary,  &
                                 donner_deep_endts,         &
                                 donner_deep, donner_deep_end,   &
                                 donner_deep_restart
use moist_conv_mod,        only: moist_conv, moist_conv_init
use lscale_cond_mod,       only: lscale_cond_init
use uw_conv_mod,           only: uw_conv_end, uw_conv_init
use lin_cld_microphys_mod, only: lin_cld_microphys_init, &
                                 lin_cld_microphys_end
use ras_mod,               only: ras_end, ras_init
use dry_adj_mod,           only: dry_adj, dry_adj_init
use strat_cloud_mod,       only: strat_cloud_init, strat_cloud_end, &
                                 strat_cloud_restart, strat_cloud_time_vary
use detr_ice_num_mod,      only: detr_ice_num, detr_ice_num_init,   &
                                 detr_ice_num_end

! ---> h1g
use mpp_mod,               only: mpp_chksum
use MG_microp_3D_mod,      only: MG_microp_3D_init, MG_microp_3D, &
                                            MG_microp_3D_end

use clubb_driver_mod,      only: clubb_init, clubb, clubb_end
! <--- h1g

use rh_clouds_mod,         only: rh_clouds_init, rh_clouds_end, &
                                 rh_clouds_sum
use diag_cloud_mod,        only: diag_cloud_init, diag_cloud_end, &
                                 diag_cloud_restart
use diag_integral_mod,     only: diag_integral_field_init, &
                                 sum_diag_integral_field
use cu_mo_trans_mod,       only: cu_mo_trans_init, cu_mo_trans, cu_mo_trans_end
use moz_hook_mod,          only: moz_hook
use rad_utilities_mod,     only: aerosol_type
use moist_proc_utils_mod,  only: capecalcnew, tempavg, column_diag, rh_calc, pmass

use moistproc_kernels_mod, only: moistproc_init, moistproc_end, moistproc_mca, &
                                 moistproc_ras, moistproc_lscale_cond,         &
                                 moistproc_strat_cloud, moistproc_cmt,         &
                                 moistproc_uw_conv, moistproc_scale_uw,        &
                                 moistproc_scale_donner,                       &
                                 rain_uw, snow_uw, ttnd_uw, qtnd_uw, utnd_uw,  &
                                 vtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw,        &
                                 qntnd_uw, qtruw, qlin, qiin, qain, delta_ql,  &
                                 delta_qi, delta_qa, qnitnd_uw
! atmos_shared modules
use atmos_tracer_utilities_mod, only : wet_deposition

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   moist_processes, moist_processes_init, moist_processes_end, &
            moist_alloc_init, moist_alloc_end,  set_cosp_precip_sources, &
            moist_processes_time_vary, moist_processes_endts, &
            doing_strat, moist_processes_restart
  

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

!--------------------- version number ----------------------------------
   character(len=128) :: &
   version = '$Id: moist_processes.F90,v 20.0 2013/12/13 23:18:25 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal $'

   character(len=5), private :: mod_name = 'moist'
   logical            :: moist_allocated = .false.
   logical            :: module_is_initialized = .false.

!-------------------- namelist data (private) --------------------------

!---------------- namelist variable definitions ------------------------
!
!   do_limit_donner = limit Donner deeo tendencies to prevent the
!                formation of grid points with negative water vapor,
!                liquid or ice.
!
!   do_limit_uw = limit UW shallow tendencies to prevent the formation
!                of grid points with negative total water specific 
!                humidities. This situation can occur because both
!                shallow and deep convection operate on the same
!                soundings without knowledge of what the other is doing
!
!   do_unified_convective_closure = use cloud base mass flux calculated
!                in uw_conv module as value for donner deep parameter-
!                ization; adjust cbmf available for uw shallow appropr-
!                iately. only available when uw shallow and donner deep
!                are the active convective schemes
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=true ]
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=true ]
!   do_ras   = switch to turn on/off relaxed arakawa shubert
!                [logical, default: do_ras=false ]
!   do_donner_deep = switch to turn on/off donner deep convection scheme
!                [logical, default: do_donner_deep=false ]
!   do_strat = switch to turn on/off stratiform cloud scheme
!                [logical, default: do_strat=false ]
!   do_rh_clouds = switch to turn on/off simple relative humidity cloud scheme
!                [logical, default: do_rh_clouds=false ]
!   do_diag_clouds = switch to turn on/off (Gordon's) diagnostic cloud scheme
!                [logical, default: do_diag_clouds=false ]
!   do_dryadj = switch to turn on/off dry adjustment scheme
!                [logical, default: do_dryadj=false ]
!   do_lin_cld_microphys = switch to turn on/off the Lin Cloud Micro-Physics scheme
!                [logical, default: do_lin_cld_microphys=false ]
!   do_liq_num = switch to turn on/off the prognostic droplet number scheme.
!                [logical, default: do_liq_num=false ]
!   use_tau  = switch to determine whether current time level (tau)
!                will be used or else future time level (tau+1).
!                if use_tau = true then the input values for t,q, and r
!                are used; if use_tau = false then input values
!                tm+tdt*dt, etc. are used.
!                [logical, default: use_tau=false ]
!
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!   tfreeze  = mean temperature used for snowfall determination (deg k)
!                [real, default: tfreeze=273.16]
!
!   do_gust_cv = switch to use convective gustiness (default = false)
!   gustmax    = maximum convective gustiness (m/s)
!   gustconst  = precip rate which defines precip rate which begins to
 !               matter for convective gustiness (kg/m2/sec)
!   cmt_mass_flux_source = parameterization(s) being used to supply the 
!                mass flux profiles seen by the cumulus momentum transport
!                module; currently either 'ras', 'donner', 'uw', 
!                'donner_and_ras', 'donner_and_uw', 'ras_and_uw', 
!                'donner_and_ras_and_uw' or 'all'
!
!   do_bm    = switch to turn on/off betts-miller scheme
!                [logical, default: do_bm=false ]
!   do_bmmass  = switch to turn on/off betts-miller massflux scheme
!                [logical, default: do_bmmass=false ]
!   do_bmomp  = switch to turn on/off olivier's version of the betts-miller 
!                scheme (with separated boundary layer)
!                [logical, default: do_bmomp=false ]
!   do_simple = switch to turn on alternative definition of specific humidity.
!                When true, specific humidity = (rdgas/rvgas)*esat/pressure
!
!   notes: 1) do_lsc and do_strat cannot both be true
!          2) pdepth and tfreeze are used to determine liquid vs. solid
!             precipitation for mca, lsc, and ras schemes, the 
!             stratiform scheme determines it's own precipitation type.
!          3) if do_strat=true then stratiform cloud tracers: liq_wat,
!             ice_wat, cld_amt must be present 
!          4) do_donner_deep and do_rh_clouds cannot both be true
!             (pending revision of code flow)
!
!-----------------------------------------------------------------------
! main convection/large-scale schemes
   logical :: do_bm=.false.
   logical :: do_bmmass =.false.
   logical :: do_bmomp  =.false.
   logical :: do_cmt=.false.
   logical :: do_diag_clouds=.false.
   logical :: do_donner_deep=.false.
   logical :: do_dryadj=.false.
   logical :: do_lin_cld_microphys=.false.
   logical :: do_lsc=.true.
   logical :: do_mca=.true. 
   logical :: do_ras=.false.
   logical :: do_rh_clouds=.false.
   logical :: do_strat=.false.
   logical :: do_uw_conv=.false.
! tracers 
   logical :: do_tracers_in_donner =.false.
   logical :: do_tracers_in_mca = .false.
   logical :: do_tracers_in_ras = .false.
   logical :: do_tracers_in_uw = .false.
! donner specific 
   logical :: do_donner_before_uw = .false.
   logical :: do_donner_mca=.true.
   logical :: do_donner_conservation_checks = .false.
   logical :: do_limit_donner = .false. ! .false. produces previous 
                                        ! behavior (cjg)
   logical :: force_donner_moist_conserv = .false.
! cmt specific
   logical :: cmt_uses_donner = .false.
   logical :: cmt_uses_ras = .false.
   logical :: cmt_uses_uw  = .false.
! others
   logical :: doing_diffusive
   logical :: use_updated_profiles_for_uw = .false.
   logical :: only_one_conv_scheme_per_column = .false.
   logical :: limit_conv_cloud_frac = .false.

! ---> h1g
   real    :: conv_frac_max = 0.99
   logical :: use_updated_profiles_for_clubb = .false.
   logical :: remain_detrain_bug = .false.
! <--- h1g

   logical :: include_donmca_in_cosp = .true.
   logical :: use_tau=.false.
   logical :: do_gust_cv = .false.
   logical :: do_liq_num = .false.
   logical :: do_simple =.false.
   logical :: do_unified_convective_closure = .false.
   logical :: do_limit_uw = .false.     ! .false. produces previous
                                        ! behavior (cjg )
   logical :: using_fms = .true.
   logical :: do_ice_num=.false.
   logical :: detrain_liq_num=.false.
   logical :: detrain_ice_num =.false.
   logical :: do_legacy_strat_cloud = .true.
   character(len=64)  :: cmt_mass_flux_source = 'ras'

   integer :: tau_sg = 0
   integer :: k_sg = 2

   real :: pdepth = 150.e2
   real :: gustmax = 3.                    ! maximum gustiness wind (m/s)
   real :: gustconst = 10./SECONDS_PER_DAY ! constant in kg/m2/sec, default =
                                           ! 1 cm/day = 10 mm/day

namelist /moist_processes_nml/ do_mca, do_lsc, do_ras, do_uw_conv, do_strat,     &
                               do_donner_before_uw, use_updated_profiles_for_uw, &
                               only_one_conv_scheme_per_column, do_diag_clouds,  &
                               limit_conv_cloud_frac, do_dryadj, pdepth,         &
                               include_donmca_in_cosp, &
                               do_unified_convective_closure, tau_sg, k_sg,      &
                               do_lin_cld_microphys, use_tau, do_rh_clouds,      &
                               cmt_mass_flux_source, do_donner_deep, do_cmt,     &
                               do_gust_cv, cmt_mass_flux_source, gustmax,        &
                               gustconst, do_liq_num, force_donner_moist_conserv,&
                               do_donner_conservation_checks, do_donner_mca,     &
                               do_limit_uw, do_limit_donner, using_fms,          &
                               do_bm, do_bmmass, do_bmomp, do_simple, &
                               do_ice_num, do_legacy_strat_cloud, &
                               detrain_liq_num, detrain_ice_num,  &
                               conv_frac_max, use_updated_profiles_for_clubb, remain_detrain_bug !h1g

!-------------------- clock definitions --------------------------------

integer :: convection_clock, largescale_clock, donner_clock, mca_clock, ras_clock, &
           donner_mca_clock, bm_clock, cmt_clock, closure_clock, lscalecond_clock, &
           stratcloud_clock, shallowcu_clock

!-------------------- diagnostics fields -------------------------------
! ---> h1g, dump cell and neso cloud fraction from donner-deep, 2011-08-08
integer :: id_cell_cld_frac,  id_meso_cld_frac, id_donner_humidity_area
! <--- h1g, dump cell and neso cloud fraction from donner-deep, 2011-08-08

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_snow_tot, id_tot_cld_amt, id_conv_freq, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_WVP, id_LWP, id_IWP, id_AWP, id_gust_conv, &

           id_tot_cloud_area,  id_tot_liq_amt,  id_tot_ice_amt,  &
           id_tot_h2o, id_tot_vapor, &
           id_lsc_cloud_area,  id_lsc_liq_amt,  id_lsc_ice_amt,  &
           id_conv_cloud_area, id_conv_liq_amt, id_conv_ice_amt, &
           id_LWP_all_clouds,  id_IWP_all_clouds, id_WP_all_clouds, &

           id_tdt_dadj, id_rh,  id_qs, id_mc, id_mc_donner, id_mc_full, &
           id_mc_donner_half, &
           id_rh_cmip, id_mc_conv_up, id_mc_half, &
           id_conv_cld_base, id_conv_cld_top, &
           id_tdt_deep_donner, id_qdt_deep_donner, &
           id_qadt_deep_donner, id_qldt_deep_donner, &
           id_qidt_deep_donner, &
           id_qndt_deep_donner,  id_qnidt_deep_donner, &
           id_tdt_mca_donner, id_qdt_mca_donner, &
           id_prec_deep_donner, id_prec_mca_donner,&
           id_tdt_uw, id_qdt_uw, &
           id_qadt_uw, id_qldt_uw, id_qidt_uw, id_qndt_uw, id_qnidt_uw, &
           id_prec1_deep_donner, &
           id_snow_deep_donner, id_snow_mca_donner, &
           id_qadt_ls, id_qldt_ls, id_qndt_ls, id_qidt_ls, id_qnidt_ls, &
           id_qadt_conv, id_qldt_conv, id_qndt_conv, id_qidt_conv, &
           id_qnidt_conv, &
           id_qa_ls_col, id_ql_ls_col, id_qn_ls_col, id_qi_ls_col, &
           id_qni_ls_col, &
           id_qa_conv_col, id_ql_conv_col, id_qn_conv_col,  &
           id_qni_conv_col, id_qi_conv_col, &
           id_bmflag, id_klzbs, id_invtaubmt, id_invtaubmq, &
           id_massflux, id_entrop_ls, &
           id_cape, id_cin, id_tref, id_qref, &
           id_q_conv_col, id_q_ls_col, id_t_conv_col, id_t_ls_col, &
           id_enth_moist_col, id_wat_moist_col, &
           id_enth_ls_col, id_wat_ls_col, &
           id_enth_conv_col, id_wat_conv_col, &
           id_enth_donner_col, id_wat_donner_col, &
           id_enth_donner_col2,  &
           id_enth_donner_col3,  &
           id_enth_donner_col4,  &
           id_enth_donner_col5,  &
           id_enth_donner_col6,  &
           id_enth_donner_col7,  &
           id_enth_mca_donner_col, id_wat_mca_donner_col, &
           id_enth_uw_col, id_wat_uw_col, &
           id_scale_donner, id_scale_uw, &
           id_ras_precip, id_ras_freq, id_don_precip, id_don_freq, &
           id_lsc_precip, id_lsc_freq, id_uw_precip, id_uw_snow, &
           id_uw_freq, &
           id_prod_no, id_m_cdet_donner, id_m_cellup, &
           id_conv_rain3d, id_conv_snow3d,   &
           id_lscale_rain3d, id_lscale_snow3d, id_lscale_precip3d
 
integer :: id_qvout, id_qaout, id_qlout, id_qiout
integer :: id_qnout, id_qniout

integer :: id_vaporint, id_condensint, id_precipint, id_diffint
integer :: id_vertmotion
integer :: id_max_enthalpy_imbal_don, id_max_water_imbal_don
integer :: id_max_enthalpy_imbal, id_max_water_imbal
integer :: id_enthint, id_lprcp, id_lcondensint, id_enthdiffint
integer :: id_wetdep_om, id_wetdep_SOA, id_wetdep_bc, &
           id_wetdep_so4, id_wetdep_so2, id_wetdep_DMS, &
           id_wetdep_NH4NO3, id_wetdep_salt, id_wetdep_dust
integer :: id_f_snow_berg, id_f_snow_berg_cond, id_f_snow_berg_wtd

integer, dimension(:), allocatable :: id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col, &
                                      id_tracerdt_mcadon, &
                                      id_tracerdt_mcadon_col, &
                                      id_wetdep, &
                                      id_wet_deposition
real :: missing_value = -999.

!-------------------- individual scheme tracers ------------------------
   logical, dimension(:), allocatable :: tracers_in_donner, tracers_in_uw, &
                                         tracers_in_mca, tracers_in_ras
   integer :: num_donner_tracers=0
   integer :: num_mca_tracers=0
   integer :: num_ras_tracers=0
   integer :: num_uw_tracers=0
   integer :: num_tracers=0

   integer :: nbcphobic =0
   integer :: nbcphilic =0
   integer :: nomphobic =0
   integer :: nomphilic =0
   integer :: nsalt1 =0
   integer :: nsalt2 =0
   integer :: nsalt3 =0
   integer :: nsalt4 =0
   integer :: nsalt5 =0
   integer :: ndust1    =0
   integer :: ndust2    =0
   integer :: ndust3    =0
   integer :: ndust4    =0
   integer :: ndust5    =0
   integer :: nDMS      =0
   integer :: nSO2      =0
   integer :: nSO4      =0
   integer :: nSOA      =0
   integer :: nNH4NO3   =0
   integer :: nNH4      =0
   

!------------------- other global variables and parameters -------------
   real, parameter :: epst=200.

   integer :: nsphum, nql, nqi, nqa, nqn   ! tracer indices for stratiform clouds
   integer :: nqni
   integer :: nqr, nqs, nqg                ! additional tracer indices for Lin Micro-Physics
   integer :: ktop                         ! top layer index for Lin Micro-Physics
   logical :: do_cosp, donner_meso_is_largescale
   real    :: strat_precip_in_cosp = 0.
   real    :: donner_precip_in_cosp = 0.
   real    :: uw_precip_in_cosp = 0.
!-->cjg
   integer :: do_clubb
!<--cjg


!------------------ allocatable moist processes variables --------------

   real, allocatable, dimension(:,:)   :: max_enthalpy_imbal, max_water_imbal, &
                                          max_enthalpy_imbal_don, max_water_imbal_don
   real, allocatable, dimension(:,:,:) :: tin, qin, rin, uin, vin, &
                                          ttnd, qtnd, rtnd, utnd, vtnd, ttnd_don, qtnd_don, &
                                          delta_temp, delta_vapor, delta_q, &
                                          donner_humidity_area, donner_humidity_factor
   real, allocatable, dimension(:,:,:) :: delta_qni, delta_qn
   real, allocatable, dimension(:,:,:) :: nllin, nilin
   real, allocatable, dimension(:,:,:) :: tin_orig, qin_orig, tdt_init, qdt_init
   real, allocatable, dimension(:,:,:) :: qtnd_wet,  &         ! specific humidity tendency (kg/kg/s)
                                          cloud_wet, &         ! cloud liquid+ice (kg/kg)
                                          cloud_frac           ! cloud area fraction
   real, allocatable, dimension(:,:,:) :: liquid_precip, frozen_precip
   real, allocatable, dimension(:,:,:) :: frz_meso, liq_meso, frz_cell
   real, allocatable, dimension(:,:,:) :: liq_cell, mca_frz, mca_liq
   real, allocatable, dimension(:,:,:) :: frz_mesoh, liq_mesoh, frz_cellh, &
                                          liq_precflx, ice_precflx, &
                                          liq_cellh, mca_frzh, mca_liqh,&
                                          ice_precflxh, liq_precflxh
   real, allocatable, dimension(:,:) ::   sumneg
   real, allocatable, dimension(:,:,:) :: ttnd_conv, qtnd_conv
   real, allocatable, dimension(:,:,:) :: qsat, det0, det_cmt       
   real, allocatable, dimension(:,:,:) :: mc_full, mc_donner, m_cdet_donner, massflux, mc_donner_up, &
                                          mc_half, mc_donner_half
   real, allocatable, dimension(:,:,:) :: RH, wetdeptnd, q_ref, t_ref
   real, allocatable, dimension(:,:,:) :: cf, cmf
   real, allocatable, dimension(:,:,:,:) :: tracer,tracer_orig, rdt_init, &
                                            qtr, q_tnd, donner_tracer

   real, allocatable, dimension(:,:)   :: prec_intgl  

! ---> h1g, save cloud condensate tendency due to convection (20120817) 
   real, allocatable, dimension(:,:,:) :: qldt_conv, qidt_conv, qadt_conv, qndt_conv, qnidt_conv
! <--- h1g
!-----------------------------------------------------------------------

                             contains

!#######################################################################
! used to allocate variables used throughout moist_processes
!--> cjg: code modification to allow diagnostic tracers in physics_up (20120508) 
!         lx is the number of prognostic tracers
!         mx is the total number of tracers (prognostic+diagnostic)

!subroutine moist_alloc_init (ix, jx, kx, lx)
!   integer, intent(in) :: ix,jx,kx,lx

subroutine moist_alloc_init (ix, jx, kx, lx, mx)
   integer, intent(in) :: ix,jx,kx,lx,mx
!<--cjg

   if (moist_allocated) return

   allocate( tin       (ix,jx,kx))                          !; tin                    = 0.0
   allocate( qin       (ix,jx,kx))                          !; qin                    = 0.0
   allocate( rin       (ix,jx,kx))                          !; rin                    = 0.0
   allocate( uin       (ix,jx,kx))                          !; uin                    = 0.0
   allocate( vin       (ix,jx,kx))                          !; vin                    = 0.0
   allocate( tin_orig  (ix,jx,kx))                          !; tin_orig               = 0.0
   allocate( qin_orig  (ix,jx,kx))                          !; qin_orig               = 0.0
   allocate( t_ref     (ix,jx,kx))                          ; t_ref                  = 0.0
   allocate( q_ref     (ix,jx,kx))                          ; q_ref                  = 0.0
   allocate( ttnd      (ix,jx,kx))                          ; ttnd                   = 0.0
   allocate( qtnd      (ix,jx,kx))                          ; qtnd                   = 0.0
   allocate( rtnd      (ix,jx,kx))                          ; rtnd                   = 0.0
   allocate( utnd      (ix,jx,kx))                          ; utnd                   = 0.0
   allocate( vtnd      (ix,jx,kx))                          ; vtnd                   = 0.0
   allocate( ttnd_don  (ix,jx,kx))                          ; ttnd_don               = 0.0
   allocate( qtnd_don  (ix,jx,kx))                          ; qtnd_don               = 0.0
   allocate( ttnd_conv (ix,jx,kx))                          ; ttnd_conv              = 0.0
   allocate( qtnd_conv (ix,jx,kx))                          ; qtnd_conv              = 0.0
   allocate( qtnd_wet  (ix,jx,kx))                          ; qtnd_wet               = 0.0
   allocate( tdt_init  (ix,jx,kx))                          ; tdt_init               = 0.0
   allocate( qdt_init  (ix,jx,kx))                          ; qdt_init               = 0.0
   allocate( cf        (ix,jx,kx))                          ; cf                     = 0.0
   allocate( cmf       (ix,jx,kx))                          ; cmf                    = 0.0
   allocate( delta_temp(ix,jx,kx))                          ; delta_temp             = 0.0
   allocate( delta_q   (ix,jx,kx))                          ; delta_q                = 0.0
   allocate( delta_vapor(ix,jx,kx))                         ; delta_vapor            = 0.0
   allocate( donner_humidity_area(ix,jx,kx))                ; donner_humidity_area   = 0.0
   allocate( donner_humidity_factor(ix,jx,kx))              ; donner_humidity_factor = 0.0
   allocate( cloud_wet  (ix,jx,kx))                         ; cloud_wet              = 0.0
   allocate( cloud_frac (ix,jx,kx))                         ; cloud_frac             = 0.0
   allocate( liquid_precip(ix,jx,kx))                       ; liquid_precip          = 0.0
   allocate( frozen_precip(ix,jx,kx))                       ; frozen_precip          = 0.0
   allocate( ice_precflx (ix,jx,kx))                        ; ice_precflx            = 0.0
   allocate( liq_precflx (ix,jx,kx))                        ; liq_precflx            = 0.0
   allocate( frz_meso  (ix,jx,kx))                          ; frz_meso               = 0.0
   allocate( liq_meso  (ix,jx,kx))                          ; liq_meso               = 0.0
   allocate( frz_cell  (ix,jx,kx))                          ; frz_cell               = 0.0
   allocate( liq_cell  (ix,jx,kx))                          ; liq_cell               = 0.0
   allocate( mca_frz   (ix,jx,kx))                          ; mca_frz                = 0.0
   allocate( mca_liq   (ix,jx,kx))                          ; mca_liq                = 0.0
   allocate( frz_mesoh (ix,jx,kx+1))                        ; frz_mesoh              = 0.0
   allocate( liq_mesoh (ix,jx,kx+1))                        ; liq_mesoh              = 0.0
   allocate( frz_cellh (ix,jx,kx+1))                        ; frz_cellh              = 0.0
   allocate( sumneg    (ix,jx))                        ; sumneg                 = 0.0
   allocate( liq_cellh (ix,jx,kx+1))                        ; liq_cellh              = 0.0
   allocate( mca_liqh  (ix,jx,kx+1))                        ; mca_liqh               = 0.0
   allocate( mca_frzh  (ix,jx,kx+1))                        ; mca_frzh               = 0.0
   allocate( ice_precflxh(ix,jx,kx+1))                      ; ice_precflxh           = 0.0
   allocate( liq_precflxh(ix,jx,kx+1))                      ; liq_precflxh           = 0.0
   allocate( qsat      (ix,jx,kx))                          ; qsat                   = 0.0
   allocate( det0      (ix,jx,kx))                          ; det0                   = 0.0
   allocate( det_cmt   (ix,jx,kx))                          ; det_cmt                = 0.0
   allocate( mc_full   (ix,jx,kx))                          ; mc_full                = 0.0
   allocate( mc_donner (ix,jx,kx))                          ; mc_donner              = 0.0
   allocate( mc_donner_up (ix,jx,kx))                       ; mc_donner_up           = 0.0
   allocate( mc_half      (ix,jx,kx+1))                     ; mc_half                = 0.0
   allocate( mc_donner_half (ix,jx,kx+1))                   ; mc_donner_half         = 0.0
   allocate( m_cdet_donner(ix,jx,kx))                       ; m_cdet_donner          = 0.0
   allocate( massflux  (ix,jx,kx))                          ; massflux               = 0.0
   allocate( RH        (ix,jx,kx))                          ; RH                     = 0.0
! pmass defined in moist_processes_utils
   allocate( pmass     (ix,jx,kx))                          ; pmass                  = 0.0
   allocate( wetdeptnd (ix,jx,kx))                          ; wetdeptnd              = 0.0
!--> cjg: code modification to allow diagnostic tracers in physics_up (20120508) 
!   allocate(tracer     (ix,jx,kx,lx))                       ; tracer                 = 0.0
!   allocate(tracer_orig(ix,jx,kx,lx))                       ; tracer_orig            = 0.0
   allocate(tracer     (ix,jx,kx,mx))                       ; tracer                 = 0.0
   allocate(tracer_orig(ix,jx,kx,mx))                       ; tracer_orig            = 0.0
!<--cjg
   allocate(q_tnd      (ix,jx,kx,lx))                       ; q_tnd                  = 0.0
   allocate(rdt_init   (ix,jx,kx,lx))                       ; rdt_init               = 0.0
   allocate(qtr          (ix,jx,kx,num_donner_tracers))     ; qtr                    = 0.0
   allocate(donner_tracer(ix,jx,kx,num_donner_tracers))     ; donner_tracer          = 0.0
   allocate(delta_qn   (ix,jx,kx))                          ; delta_qn               = 0.0
   allocate(delta_qni  (ix,jx,kx))                          ; delta_qni              = 0.0
   allocate(nllin      (ix,jx,kx))                          ; nllin                  = 0.0
   allocate(nilin      (ix,jx,kx))                          ; nilin                  = 0.0

! ---> h1g, allocate cloud condensate tendency due to convection (20120817) 
   allocate( qldt_conv (ix,jx,kx))
   allocate( qidt_conv (ix,jx,kx))
   allocate( qadt_conv (ix,jx,kx))
   if( do_liq_num ) allocate( qndt_conv (ix,jx,kx))
   if( do_ice_num ) allocate( qnidt_conv (ix,jx,kx))
! <--- h1g

   moist_allocated = .true.
  
end subroutine moist_alloc_init


!#######################################################################
! used to deallocate variables used throughout moist_processes
subroutine moist_alloc_end

   if (moist_allocated .eqv. .false. ) return
   deallocate( tin       )
   deallocate( qin       )
   deallocate( rin       )
   deallocate( uin       )
   deallocate( vin       )
   deallocate( tin_orig  )
   deallocate( qin_orig  )
   deallocate( t_ref     )
   deallocate( q_ref     )
   deallocate( ttnd      )
   deallocate( qtnd      )
   deallocate( rtnd      )
   deallocate( utnd      )
   deallocate( vtnd      )
   deallocate( ttnd_don  )
   deallocate( qtnd_don  )
   deallocate( ttnd_conv )
   deallocate( qtnd_conv )
   deallocate( qtnd_wet  )
   deallocate( tdt_init  )
   deallocate( qdt_init  )
   deallocate( cf        )
   deallocate( cmf       )
   deallocate( delta_temp)
   deallocate( delta_q   )
   deallocate( delta_vapor )
   deallocate( donner_humidity_area)
   deallocate( donner_humidity_factor)
   deallocate( cloud_wet  )
   deallocate( cloud_frac )
   deallocate( liquid_precip)
   deallocate( frozen_precip)
   deallocate( ice_precflx)
   deallocate( liq_precflx)
   deallocate( frz_meso  )
   deallocate( liq_meso  )
   deallocate( frz_cell  )
   deallocate( liq_cell  )
   deallocate( mca_frz   )
   deallocate( mca_liq   )
   deallocate( frz_mesoh )
   deallocate( liq_mesoh )
   deallocate( frz_cellh )
   deallocate( sumneg    )
   deallocate( liq_cellh )
   deallocate( mca_frzh  )
   deallocate( mca_liqh  )
   deallocate( ice_precflxh)
   deallocate( liq_precflxh)
   deallocate( qsat      )
   deallocate( det0      )
   deallocate( det_cmt   )
   deallocate( mc_full   )
   deallocate( mc_donner )
   deallocate( mc_donner_up )
   deallocate( mc_half      )
   deallocate( mc_donner_half      )
   deallocate( m_cdet_donner)
   deallocate( massflux  )
   deallocate( RH        )
   deallocate( pmass     )
   deallocate( wetdeptnd )
   deallocate(tracer     )
   deallocate(tracer_orig)
   deallocate(q_tnd      )
   deallocate(rdt_init   )
   deallocate(qtr        )
   deallocate(donner_tracer)
   deallocate(delta_qn   )
   deallocate(delta_qni  )
   deallocate(nllin      )
   deallocate(nilin      )

! ---> h1g, deallocate cloud condensate tendency due to convection (20120817) 
   deallocate( qldt_conv )
   deallocate( qidt_conv )
   deallocate( qadt_conv )
   if( do_liq_num ) deallocate( qndt_conv )
   if( do_ice_num ) deallocate( qnidt_conv )
! <--- h1g

   moist_allocated = .false.

end subroutine moist_alloc_end

!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
                            phalf, pfull, zhalf, zfull, omega, diff_t, &
                            radturbten, cush, cbmf,                    &
                            pblht, ustar, bstar, qstar,                &
                            t, q, r, u, v, tm, qm, rm, um, vm,         &
                            tdt, qdt, rdt, udt, vdt, diff_cu_mo,       &
                            convect, lprec, fprec, fl_lsrain,          &
                            fl_lssnow, fl_ccrain, fl_ccsnow, &
                            fl_donmca_rain, fl_donmca_snow, gust_cv,  &
                            area, lon, lat, lsc_cloud_area, lsc_liquid,     &
                            lsc_ice, lsc_droplet_number, &
                            lsc_ice_number, lsc_snow, lsc_rain,  &
                            lsc_snow_size, lsc_rain_size     , &
! ---> h1g
                            dcond_ls_liquid,     dcond_ls_ice,         &
                            Ndrop_act_CLUBB,     Icedrop_act_CLUBB,    &
                            ndust, rbar_dust,                          &
                            diff_t_clubb,                              &
                            tdt_shf,                                   &
                            qdt_lhf,                                   &
! <--- h1g
                            Aerosol, mask, kbot, &
                            shallow_cloud_area, shallow_liquid,  &
                            shallow_ice, shallow_droplet_number, &
                            shallow_ice_number, &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            cell_droplet_number, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size,  &
                            meso_droplet_number, nsum_out, &
                            hydrostatic, phys_hydrostatic)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         land       fraction of surface covered by land
!                      [real, dimension(nlon,nlat)]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         omega      omega (vertical velocity) at full levels
!                    in pascals per second
!                      [real, dimension(nlon,nlat,nlev)]
!
!         diff_t     vertical diffusion coefficient for temperature
!                    and tracer (m*m/sec) on half levels
!                      [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         r          tracer fields at full model levels,
!                    at the current time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         u, v,      zonal and meridional wind [m/s] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
! 
!         tm, qm     temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the previous time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rm         tracer fields at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         um, vm     zonal and meridional wind [m/s] at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev)]
!
!         area       grid box area (in m2)
!                      [real, dimension(nlon,nlat)]
!
!         lon        longitude in radians           ! h1g
!                      [real, dimension(nlon,nlat)] ! h1g
!
!         lat        latitude in radians
!                      [real, dimension(nlon,nlat)]
!  
! inout:  tdt, qdt   temperature (tdt) [deg k/sec] and specific
!                    humidity of water vapor (qdt) tendency [1/sec]
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rdt        tracer tendencies 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         udt, vdt   zonal and meridional wind tendencies [m/s/s]
! 
!   out:  convect    is moist convection occurring in this grid box?
!                      [logical, dimension(nlon,nlat)]
!
!         lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
! 
!         gust_cv    gustiness from convection  in m/s
!                      [real, dimension(nlon,nlat)]
!
!       optional
!  -----------------
! 
!    in:  mask       mask (1. or 0.) for grid boxes above or below
!                    the ground   [real, dimension(nlon,nlat,nlev)]
!
!         kbot       index of the lowest model level
!                      [integer, dimension(nlon,nlat)]
!
!
!-----------------------------------------------------------------------
   integer,         intent(in)           :: is,ie,js,je
   type(time_type), intent(in)           :: Time
   real, intent(in)                      :: dt
   real, intent(in) , dimension(:,:)     :: land, pblht, ustar, bstar, qstar
   real, intent(inout), dimension(:,:)   :: cush, cbmf
   real, intent(in) , dimension(:,:,:)   :: phalf, pfull, zhalf, zfull, omega, &
                                            diff_t, t, q, u, v, tm, qm, um, vm
   real, dimension(:,:,:), intent(in)    :: radturbten
   real, intent(inout), dimension(:,:,:,:) :: r, rm                      ! cjg: inout
   real, intent(inout),dimension(:,:,:)  :: tdt, qdt, udt, vdt
   real, intent(inout),dimension(:,:,:,:):: rdt
logical, intent(out), dimension(:,:)     :: convect
   real, intent(out), dimension(:,:)     :: lprec, fprec, gust_cv
   real, intent(out), dimension(:,:,:)   :: fl_lsrain, fl_lssnow, &
                                            fl_ccrain, fl_ccsnow, &
                                            fl_donmca_rain, fl_donmca_snow
   real, intent(out), dimension(:,:,:)   :: diff_cu_mo
   real, intent(in) , dimension(:,:)     :: area
   real, intent(in) , dimension(:,:)     :: lon
   real, intent(in) , dimension(:,:)     :: lat

! ---> h1g
    real, intent(inout), dimension(:,:,:), optional :: dcond_ls_liquid, dcond_ls_ice
    real, intent(inout), dimension(:,:,:), optional :: Ndrop_act_CLUBB,  Icedrop_act_CLUBB
    real, intent(inout), dimension(:,:,:), optional :: ndust, rbar_dust
    real, intent(inout), dimension(:,:,:), optional :: diff_t_clubb
    real, intent(inout), dimension(:,:),   optional :: tdt_shf,  qdt_lhf
! < --- h1g

   real, intent(out) , dimension(:,:,:)  ::   &
                       lsc_cloud_area, lsc_liquid, lsc_ice,   &
                       lsc_droplet_number, lsc_ice_number, lsc_snow, &
                       lsc_rain, lsc_snow_size, lsc_rain_size
   type(aerosol_type),intent(in),       optional :: Aerosol
   real, intent(in) , dimension(:,:,:), optional :: mask
   integer, intent(in), dimension(:,:), optional :: kbot

   logical, intent(in), optional :: hydrostatic, phys_hydrostatic
   integer, intent(inout), dimension(:,:), optional ::  nsum_out
   real, intent(inout), dimension(:,:,:), optional :: &      
                                  shallow_cloud_area, shallow_liquid,   &
                                  shallow_ice, shallow_droplet_number, &
                                  shallow_ice_number, &
                                  cell_cld_frac, cell_liq_amt, cell_liq_size, &
                                  cell_ice_amt, cell_ice_size, &
                                  cell_droplet_number, &
                                  meso_cld_frac, meso_liq_amt, meso_liq_size, &
                                  meso_ice_amt, meso_ice_size, &
                                  meso_droplet_number

!-----------------------------------------------------------------------
   integer :: secs, days
   integer :: n, nn, i, j, k, ix, jx, kx, nt, tr
   integer :: m, mm
   logical :: used, avgbl
   real    :: dtinv

   real, dimension(size(t,1),size(t,2)) :: cape, cin
   real, dimension(size(t,1),size(t,2)) :: precip, total_precip, lheat_precip, &
                                           precip_returned, precip_adjustment, &
                                           vert_motion
   real, dimension(size(t,1),size(t,2)) :: rain, snow, &
                                           rain_don, snow_don, &
                                           rain_ras, snow_ras, &
                                           rain_donmca, snow_donmca
   real, dimension(size(t,1),size(t,2)) :: bmflag, klzbs, invtaubmt, invtaubmq
   real, dimension(size(t,1),size(t,2)) :: scale
   real, dimension(size(t,1),size(t,2)) :: freq_count
   real, dimension(size(t,1),size(t,2)) :: enthint, lcondensint, enthdiffint,  &
                                           vaporint, condensint, precipint, diffint

   real, dimension(size(t,1),size(t,2),size(phalf,3)) :: rain3d, snow3d
   real, dimension(size(t,1),size(t,2),size(phalf,3)) :: snowclr3d
   real, dimension(size(t,1),size(t,2),size(t,3)+1) :: mc, m_cellup, mc_cmt
   real, dimension(size(t,1),size(t,2),size(pfull,3)) :: f_snow_berg


!     sfc_sh_flux      sensible heat flux across the surface
!                      [ watts / m**2 ]
!     sfc_vapor_flux   water vapor flux across the surface
!                      [ kg(h2o) / (m**2 sec) ]
!     tr_flux          tracer fux across the surface
!                      [ kg(tracer) / (m**2 sec) ]
   real, dimension(size(t,1),size(t,2)) :: sfc_sh_flux, sfc_vapor_flux
   real, dimension(size(t,1),size(t,2),num_donner_tracers) :: tr_flux  
   real, dimension(size(t,1),size(t,2),num_donner_tracers) :: &
                                                          donner_wetdep
   real, dimension(size(t,1),size(t,2),num_uw_tracers) :: &
                                                          uw_wetdep
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) :: total_wetdep
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) ::  &
                                                       total_wetdep_uw
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) ::   &
                                                     total_wetdep_donner
   real, dimension(size(t,1),size(t,2),size(rdt,4)   ) :: ls_wetdep
   real, dimension(size(t,1),size(t,2),size(t,3) ) :: total_conv_cloud,&
                           conv_cld_frac, tot_conv_liq, tot_conv_ice

!chemistry start
   real, parameter :: boltz = 1.38044e-16
   integer, dimension(size(rdt,1),size(rdt,2)) :: cldtop, cldbot
   real, dimension(size(rdt,1),size(rdt,2),size(rdt,3)) :: prod_no
   real, dimension(size(rdt,1),size(rdt,2),size(rdt,3),size(rdt,4)) :: wet_data
!chemistry end

   real, dimension(size(t,1),size(t,2))           ::  adjust_frac      
   real, dimension(size(t,1),size(t,2),size(t,3)) ::  ttnd_adjustment
   real, dimension(size(t,1),size(t,2),size(t,3)) ::  available_cf_for_uw

   logical, dimension(size(t,1),size(t,2)) :: conv_calc_completed
   logical, dimension(size(t,1),size(t,2)) :: coldT

!temporary variables
   real :: temp
   logical, dimension(size(t,1),size(t,2)) :: ltemp
   real, dimension(size(t,1),size(t,2)) :: temp_2d
   real, dimension(size(t,1),size(t,2)) :: tca2
   real, dimension(size(t,1),size(t,2),size(t,3)) :: total_cloud_area
   real, dimension(size(t,1),size(t,2),size(t,3)) :: temp_3d1, temp_3d2, temp_3d3

! ---> h1g, 2010-08-23
   real                                 ::       current_total_sec
   integer                              ::       current_sec, current_days

!  consider the donner-deep mass flux impacts on clubb
   real, dimension(size(omega,1),size(omega,2),size(omega,3))  :: conv_frac_clubb
   real, dimension(size(omega,1),size(omega,2),size(omega,3))  :: convective_humidity_ratio_clubb
   real                                                        :: qrf, env_fraction, env_qv
! <--- h1g, 2010-08-23


! ---> h1g, 2012-10-05
   real, dimension(size(omega,1),size(omega,2),size(omega,3))  :: qcvar_clubb
! <--- h1g, 2012-10-05      
!-------- input array size and position in global storage --------------
      ix=size(t,1); jx=size(t,2); kx=size(t,3); nt=size(rdt,4)

       
!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('moist_processes_mod',  &
                 'moist_processes_init has not been called.', FATAL)
      endif

      conv_calc_completed = .false.
      available_cf_for_uw = 1.0

!--------------------------------------------------------------------
!    define the inverse of the time step.
!--------------------------------------------------------------------
      dtinv = 1.0/dt

!--------------------------------------------------------------------
!    initialize the arrays which will be used in this subroutine.
!--------------------------------------------------------------------
      rain_don     = 0.0
      snow_don     = 0.0
      rain_donmca  = 0.0
      snow_donmca  = 0.0
      lprec        = 0.0  
      fprec        = 0.0
      fl_lsrain(:,:,:) = 0.
      fl_lssnow(:,:,:) = 0.
      fl_ccrain(:,:,:) = 0.
      fl_ccsnow(:,:,:) = 0.
      fl_donmca_rain(:,:,:) = 0.
      fl_donmca_snow(:,:,:) = 0.
      convect      = .false.
      gust_cv      = 0.0
      precip       = 0.0 
      rain3d       = 0.0
      snow3d       = 0.0

!---------------------------------------------------------------------
!    initialize local arrays which will hold sums.
!---------------------------------------------------------------------
      rdt_init(is:ie,js:je,:,:)  = rdt
      tdt_init(is:ie,js:je,:)  = tdt
      qdt_init(is:ie,js:je,:)  = qdt
!      ttnd_conv(is:ie,js:je,:) = 0.
!      qtnd_conv(is:ie,js:je,:) = 0.
!      qtnd(is:ie,js:je,:)      = 0.
!      q_tnd(is:ie,js:je,:,:)     = 0.

!---------------------------------------------------------------------
!    define input fields to be used, either the tau time level fields,
!    or the tau - 1 time level values updated with the time tendencies
!    thus far calculated on the current step. control is through nml
!    variable use_tau.
!---------------------------------------------------------------------
      if (use_tau) then
        tin(is:ie,js:je,:) = t
        qin(is:ie,js:je,:) = q
        uin(is:ie,js:je,:) = u
        vin(is:ie,js:je,:) = v
        do tr=1,size(r,4)
          tracer(is:ie,js:je,:,tr) = r(:,:,:,tr)
        end do  
      else
        tin(is:ie,js:je,:) = tm + tdt*dt
        qin(is:ie,js:je,:) = qm + qdt*dt
        uin(is:ie,js:je,:) = um + udt*dt
        vin(is:ie,js:je,:) = vm + vdt*dt
        do tr=1,size(rdt,4)
          tracer(is:ie,js:je,:,tr) = rm(:,:,:,tr) + rdt(:,:,:,tr)*dt
        end do  
        do tr=size(rdt,4) +1, size(r,4)
          tracer(is:ie,js:je,:,tr) = r(:,:,:,tr)
        end do  
      endif

!--------------------------------------------------------------------
!    if using eta vertical coordinate, define the appropriate values 
!    for any points located below the ground. values of 0.0 are given
!    to u, v and q, and a temperature value of EPST (=200. K) is given 
!    to sub-surface  points.
!--------------------------------------------------------------------
      if (present(mask) .and. present(kbot))  then
        tin(is:ie,js:je,:) = mask*tin(is:ie,js:je,:) + (1.0 - mask)*EPST 
        qin(is:ie,js:je,:) = mask*qin(is:ie,js:je,:)
        uin(is:ie,js:je,:) = mask*uin(is:ie,js:je,:)
        vin(is:ie,js:je,:) = mask*vin(is:ie,js:je,:)
        do tr=1,size(r,4)
          tracer(is:ie,js:je,:,tr) = mask(:,:,:)*tracer(is:ie,js:je,:,tr)
        end do  
      endif
   
!----------------------------------------------------------------------
!    compute the mass in each model layer.
!----------------------------------------------------------------------
      do k=1,kx
        pmass(is:ie,js:je,k) = (phalf(:,:,k+1) - phalf(:,:,k))/GRAV
      end do

!----------------------------------------------------------------------
!    output any requested convectively-transported tracer fields 
!    and / or their column sums before convective transport.
!----------------------------------------------------------------------
      do n=1,num_tracers
        used = send_data (id_conv_tracer(n), tracer(is:ie,js:je,:,n), Time, &
                          is, js, 1, rmask=mask)
        if (id_conv_tracer_col(n) > 0)  &
          call column_diag(id_conv_tracer_col(n), is, js, Time, &
                           tracer(is:ie,js:je,:,n), 1.0) 
      end do

!----------------------------------------------------------------------
!    compute the mean temperature in the lower atmosphere (the lowest
!    pdepth Pa), to be used to determine whether rain or snow reaches
!    the surface. define a logical variable coldT indicating whether
!    snow or rain falls in the column.
!    ????    SHOULD TIN BE USED RATHER THAN t ??
!----------------------------------------------------------------------
      call tempavg (pdepth, phalf, t, snow, mask)
      coldT = .false.
      where (snow(:,:) <= TFREEZE)
        coldT(:,:) = .true.
      endwhere
      
!---------------------------------------------------------------------
!    begin the clock timing the dry and moist convection parameter-
!    izations.
!---------------------------------------------------------------------
      call mpp_clock_begin (convection_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   DRY CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if dry adjustment is desired call subroutine dry_adj to obtain
!    the temperature tendencie swhich must be applied to adjust each
!    column to a non-superadiabatic lapse rate. 
!---------------------------------------------------------------------
      if (do_dryadj) then
        call dry_adj (tin(is:ie,js:je,:), pfull, phalf, delta_temp(is:ie,js:je,:), mask)

!-------------------------------------------------------------------
!    add the temperature change due to dry adjustment to the current
!    temperature. convert the temperature change to a heating rate and
!    add that to the temperature temndency array accumulating the ten-
!    dencies due to all physics processes.
!-------------------------------------------------------------------
        tin(is:ie,js:je,:)  = tin(is:ie,js:je,:) + delta_temp(is:ie,js:je,:)
        ttnd(is:ie,js:je,:) = delta_temp(is:ie,js:je,:)*dtinv
        tdt  = tdt + ttnd(is:ie,js:je,:)

!---------------------------------------------------------------------
!    output the temperature tendency from dry adjustment, if desired.
!---------------------------------------------------------------------
        used = send_data (id_tdt_dadj, ttnd(is:ie,js:je,:), Time, is, js, 1, rmask=mask )

!---------------------------------------------------------------------
!    add the temperature time tendency from dry adjustment to the array
!    accumulating the total temperature time tendency from convection.
!---------------------------------------------------------------------
        ttnd_conv(is:ie,js:je,:) = ttnd_conv(is:ie,js:je,:) + ttnd(is:ie,js:je,:)
      endif


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  MOIST CONVECTION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                0. UW SHALLOW CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  cmf(is:ie,js:je,:) = 0.
  tracer_orig(is:ie,js:je,:,:) = tracer(is:ie,js:je,:,:)
  if (.not. do_donner_before_uw) then
    call mpp_clock_begin (shallowcu_clock)
    if (do_uw_conv) then
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the uw_conv param-
!    eterization are present.
!---------------------------------------------------------------------
      if    &
         (present (shallow_cloud_area) .and.   &
          present (shallow_liquid) .and.   &
          present (shallow_ice) .and.  &
          present ( shallow_droplet_number)  .and. &
          present ( shallow_ice_number) ) then
      else
       call error_mesg ('moist_processes_mod', 'moist_processes: &
              &not all 4 optional arguments needed for uw_conv &
            &output are present', FATAL)
      endif

      call moistproc_uw_conv(Time, is, ie, js, je, dt, tin(is:ie,js:je,:), qin(is:ie,js:je,:), &
                             uin(is:ie,js:je,:), vin(is:ie,js:je,:), tracer(is:ie,js:je,:,:),    &
                             pfull, phalf, zfull, zhalf, omega, pblht,        &
                             ustar, bstar, qstar, land, coldT, Aerosol,       &
                             cush, cbmf, cmf(is:ie,js:je,:), conv_calc_completed,            &
                             available_cf_for_uw, tdt, qdt, udt, vdt, rdt,    &
                             ttnd_conv(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:), lprec, fprec, precip,      &
                             liq_precflx(is:ie,js:je,:),  &
                             ice_precflx(is:ie,js:je,:), &
                             do_strat, do_limit_uw, do_liq_num, num_tracers,  &
                             tracers_in_uw, num_uw_tracers, shallow_cloud_area,&
                             shallow_liquid, shallow_ice, shallow_droplet_number, uw_wetdep, &
                             do_ice_num, detrain_ice_num)                
    endif  !(do_uw_conv)
    call mpp_clock_end   (shallowcu_clock)
  else
    tin_orig(is:ie,js:je,:) = tin(is:ie,js:je,:)
    qin_orig(is:ie,js:je,:) = qin(is:ie,js:je,:)
!   tracer_orig = tracer
  endif  ! (.not do_donner_before_uw)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                A. DONNER DEEP CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if donner_deep convection is activated, execute the following code.
!---------------------------------------------------------------------
  if (do_donner_deep) then
    call mpp_clock_begin (donner_clock)
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the donner param-
!    eterization are present.
!---------------------------------------------------------------------
    if    &
      (present (cell_cld_frac) .and.   &
      present (cell_liq_amt) .and. present ( cell_liq_size) .and. &
      present (cell_ice_amt) .and. present ( cell_ice_size) .and. &
      present (cell_droplet_number) .and. &
      present (meso_cld_frac) .and.   &
      present (meso_liq_amt) .and. present ( meso_liq_size) .and. &
      present (meso_ice_amt) .and. present ( meso_ice_size) .and. &
      present (meso_droplet_number) .and. &
      present (nsum_out) ) then
    else
      call error_mesg ('moist_processes_mod', 'moist_processes: &
              &not all 13 optional arguments needed for donner_deep &
              &output are present', FATAL)
    endif

!--------------------------------------------------------------------
!    if strat_cloud_mod is activated, define the cloud liquid and 
!    cloud ice specific humidities and cloud area associated with 
!    strat_cloud_mod, so that they may be input to donner_deep_mod. 
!    if strat_cloud_mod is not activated, define these arrays to be 
!    zero. 
!--------------------------------------------------------------------
    if (do_strat) then
      qlin(is:ie,js:je,:) = tracer(is:ie,js:je,:,nql)
      qiin(is:ie,js:je,:) = tracer(is:ie,js:je,:,nqi)
      qain(is:ie,js:je,:) = tracer(is:ie,js:je,:,nqa)
      IF ( do_liq_num ) nllin(is:ie,js:je,:) =  tracer(is:ie,js:je,:,nqn)
      IF ( do_ice_num ) nilin(is:ie,js:je,:) =  tracer(is:ie,js:je,:,nqni)
    endif

!--------------------------------------------------------------------
!    convert vapor specific humidity to vapor mixing ratio so it may
!    be input to donner_deep_mod.
!--------------------------------------------------------------------
    rin(is:ie,js:je,:) = qin(is:ie,js:je,:)/(1.0 - qin(is:ie,js:je,:))

!---------------------------------------------------------------------
!    if any tracers are to be transported by donner convection, 
!    check each active tracer to find those to be transported and fill 
!    the donner_tracers array with these fields.
!---------------------------------------------------------------------
!    donner_tracer(is:ie,js:je,:,:) = 0.0
    nn = 1
    do n=1,num_tracers
      if (tracers_in_donner(n)) then
        donner_tracer(is:ie,js:je,:,nn) = tracer(is:ie,js:je,:,n)
        nn = nn + 1
      endif
    end do

!---------------------------------------------------------------------
!  NOTE 1: sfc_sh_flux, sfc_vapor_flux, tr_flux are the surface fluxes
!          that will have been obtained from the flux exchange module
!          and passed on to moist_processes and then to donner_deep.
!          FOR NOW, these values are defined herein, and given
!          values of 0.0
!---------------------------------------------------------------------
!   sfc_sh_flux    = INPUT_SFC_SH_FLUX_FROM_COUPLER
!   sfc_vapor_flux = INPUT_SFC_VAPOR_FLUX_FROM_COUPLER
    sfc_sh_flux    = 0.0
    sfc_vapor_flux = 0.0
    tr_flux        = 0.0
    nn = 1
    do n=1,num_tracers
      if (tracers_in_donner(n)) then
!       tr_flux(:,:,nn) = INPUT_SFC_FLUX_FROM_COUPLER(:,:,n)
        tr_flux(:,:,nn) = 0.0                                 
        nn = nn + 1
      endif
    end do
    temp_2d=pblht
    temp_2d=min(max(temp_2d, 0.0),5000.);
    temp_2d=ustar**3.+0.6*ustar*bstar*temp_2d
    where (temp_2d .gt. 0.)
      temp_2d = temp_2d**(2./3.)
    end where
    temp_2d = MAX (1.e-6, temp_2d)

!---------------------------------------------------------------------
!    call donner_deep to compute the effects of deep convection on the 
!    temperature, vapor mixing ratio, tracers, cloud liquid, cloud ice
!    cloud area and precipitation fields.
!---------------------------------------------------------------------
    call get_time (Time, secs, days)
    if (do_strat) then
      call donner_deep (is, ie, js, je, dt, tin(is:ie,js:je,:), rin(is:ie,js:je,:), pfull,        &
                        phalf, zfull, zhalf, omega, pblht, temp_2d, &
                        qstar, cush, coldT, land, sfc_sh_flux,      &!miz
                        sfc_vapor_flux, tr_flux,                    &
                        donner_tracer(is:ie,js:je,:,:), secs, days, cbmf,            &
                        cell_cld_frac, cell_liq_amt, cell_liq_size, &
                        cell_ice_amt, cell_ice_size,                &
                        cell_droplet_number,                        &
                        meso_cld_frac, meso_liq_amt, meso_liq_size, &
                        meso_ice_amt, meso_ice_size,                &
                        meso_droplet_number, nsum_out,              &
                        precip_returned, delta_temp(is:ie,js:je,:), delta_vapor(is:ie,js:je,:),   &
                        m_cdet_donner(is:ie,js:je,:), m_cellup, mc_donner(is:ie,js:je,:),         &
                        mc_donner_up(is:ie,js:je,:), mc_donner_half(is:ie,js:je,:), &
                        donner_humidity_area(is:ie,js:je,:), donner_humidity_factor(is:ie,js:je,:),&
                        qtr(is:ie,js:je,:,:),  &
                        donner_wetdep, &
                        lheat_precip, vert_motion,             &
                        total_precip, liquid_precip(is:ie,js:je,:), frozen_precip(is:ie,js:je,:), &
                     frz_meso(is:ie,js:je,:), liq_meso(is:ie,js:je,:), &
                  frz_cell(is:ie,js:je,:), liq_cell(is:ie,js:je,:), &
                        qlin(is:ie,js:je,:), qiin(is:ie,js:je,:), qain(is:ie,js:je,:), delta_ql(is:ie,js:je,:),                 &!optional
                        delta_qi(is:ie,js:je,:), delta_qa(is:ie,js:je,:))                          !optional  
    else
      call donner_deep (is, ie, js, je, dt, tin(is:ie,js:je,:), rin(is:ie,js:je,:), pfull,        &
                        phalf, zfull, zhalf, omega, pblht, temp_2d, &
                        qstar, cush, coldT, land, sfc_sh_flux,      &!miz
                        sfc_vapor_flux, tr_flux,                    &
                        donner_tracer(is:ie,js:je,:,:), secs, days,  cbmf,           &
                        cell_cld_frac, cell_liq_amt, cell_liq_size, &
                        cell_ice_amt, cell_ice_size,                &
                        cell_droplet_number,                        &
                        meso_cld_frac, meso_liq_amt, meso_liq_size, &
                        meso_ice_amt, meso_ice_size,                &
                        meso_droplet_number, nsum_out,              &
                        precip_returned, delta_temp(is:ie,js:je,:), delta_vapor(is:ie,js:je,:),   &
                        m_cdet_donner(is:ie,js:je,:), m_cellup, mc_donner(is:ie,js:je,:),         &
                        mc_donner_up(is:ie,js:je,:), mc_donner_half(is:ie,js:je,:), &
                        donner_humidity_area(is:ie,js:je,:), donner_humidity_factor(is:ie,js:je,:),&
                        qtr(is:ie,js:je,:,:), donner_wetdep, &
                        lheat_precip, vert_motion,             &
                        total_precip, liquid_precip(is:ie,js:je,:), &
                        frozen_precip(is:ie,js:je,:), &
                frz_meso(is:ie,js:je,:), liq_meso(is:ie,js:je,:),  &
                frz_cell(is:ie,js:je,:), liq_cell(is:Ie,js:je,:))
    endif

!---------------------------------------------------------------------
!    update the current timestep tracer changes with the contributions 
!    just obtained from donner transport.
!---------------------------------------------------------------------
    nn = 1
    do n=1, num_tracers
      if (tracers_in_donner(n)) then
        rdt(:,:,:,n) = rdt(:,:,:,n) + qtr(is:ie,js:je,:,nn)
        nn = nn + 1
      endif
    end do

    if (do_donner_conservation_checks) then
      vaporint = 0.
      lcondensint = 0.
      condensint = 0.
      diffint = 0.
      enthint = 0.
      enthdiffint = 0.
    
      do k=1,kx
        vaporint(:,:) = vaporint(:,:) + pmass (is:ie,js:je,k)*delta_vapor(is:ie,js:je,k)
        enthint(:,:) = enthint(:,:) + CP_AIR*pmass(is:ie,js:je,k)*delta_temp(is:ie,js:je,k)
        condensint(:,:) = condensint(:,:) + pmass(is:ie,js:je,k) *  &
                         (delta_ql(is:ie,js:je,k) + delta_qi(is:ie,js:je,k))
        lcondensint(:,:) = lcondensint(:,:) + pmass(is:ie,js:je,k) *  &
                         (HLV*delta_ql(is:ie,js:je,k) + HLS*delta_qi(is:ie,js:je,k))
      end do
      precipint = total_precip/seconds_per_day
      diffint = (vaporint + condensint)*dtinv  + precipint
      enthdiffint = (enthint - lcondensint)*dtinv -    &
                   lheat_precip/seconds_per_day - vert_motion/seconds_per_day 
      do j=1,size(enthdiffint,2)
       do i=1,size(enthdiffint,1)
         max_enthalpy_imbal_don(i,j) = max( abs(enthdiffint(i,j)), &
                                         max_enthalpy_imbal_don(i,j) )
         max_water_imbal_don(i,j) = max( abs(diffint(i,j)), &
                                         max_water_imbal_don(i,j) )
       end do
      end do

      used = send_data(id_max_enthalpy_imbal_don, max_enthalpy_imbal_don, Time, is, js)
      used = send_data(id_max_water_imbal_don, max_water_imbal_don, Time, is, js)
      used = send_data(id_vaporint, vaporint*dtinv, Time, is, js)
      used = send_data(id_condensint, condensint*dtinv, Time, is, js)
      used = send_data(id_vertmotion, vert_motion/seconds_per_day, Time, is, js)
      used = send_data(id_precipint, precipint, Time, is, js)
      used = send_data(id_diffint, diffint, Time, is, js)
      used = send_data(id_enthint, enthint*dtinv, Time, is, js)
      used = send_data(id_lcondensint, lcondensint*dtinv, Time, is, js)
      used = send_data(id_lprcp, lheat_precip/seconds_per_day, Time, is, js)
      used = send_data(id_enthdiffint, enthdiffint, Time, is, js)
    endif

!--------------------------------------------------------------------
!    obtain updated vapor specific humidity (qnew) resulting from deep 
!    convection. define the vapor specific humidity change due to deep 
!    convection (qtnd).
!--------------------------------------------------------------------
    do k=1,kx
     do j=js,je
      do i=is,ie
        if (delta_vapor(i,j,k) /= 0.0) then
!was qnew... now temp
          temp = (rin(i,j,k) + delta_vapor(i,j,k))/   &
                 (1.0 + (rin(i,j,k) + delta_vapor(i,j,k)))
          delta_q(i,j,k) = temp - qin(i,j,k)
        else
          delta_q(i,j,k) = 0.
        endif
      enddo
     enddo
    end do

!---------------------------------------------------------------------
!    scale Donner tendencies to prevent the formation of negative
!    total water specific humidities
!---------------------------------------------------------------------
    if (do_strat .and. do_limit_donner) then
      call moistproc_scale_donner(is,ie,js,je,qin(is:ie,js:je,:), delta_temp(is:ie,js:je,:), delta_q(is:ie,js:je,:), &
                                  precip_returned, total_precip, lheat_precip, liquid_precip(is:ie,js:je,:),    &
                                  frozen_precip(is:ie,js:je,:), num_tracers, tracers_in_donner,&
                                  qtr(is:ie,js:je,:,:), scale)
      used = send_data (id_scale_donner, scale, Time, is, js )
    else
      scale = 1.0
      used = send_data (id_scale_donner, scale, Time, is, js )
    end if ! (do_strat and do_limit_donner)

!---------------------------------------------------------------------
!    recalculate the precip using the delta specific humidity tenden-
!    cies. define precip_adjustment as the change in precipitation 
!    resulting from the recalculation.
!---------------------------------------------------------------------
    if (force_donner_moist_conserv) then
!---------------------------------------------------------------------
!    calculate the precipitation needed to balance the change in water
!    content in the column.
!---------------------------------------------------------------------
      temp_2d = 0.
      do k=1,kx
        temp_2d (:,:) = temp_2d (:,:) + (-delta_q(is:ie,js:je,k) -  &
                        delta_ql(is:ie,js:je,k) -delta_qi(is:ie,js:je,k))*  &
                        pmass(is:ie,js:je,k)
      end do
      precip_adjustment = (temp_2d - precip_returned)
      do j=1,jx
       do i=1,ix
         if (ABS(precip_adjustment(i,j)) < 1.0e-10) then
           precip_adjustment (i,j) = 0.0
         endif
       end do
      end do
!----------------------------------------------------------------------
!    now adjust the temperature change to balance the precip adjustment
!    and so conserve enthalpy in the column.
!--------------------------------------------------------------------- 
      do j=1,jx
       do i=1,ix
         if (precip_returned(i,j) > 0.0) then
           adjust_frac(i,j) = precip_adjustment(i,j)/precip_returned(i,j)
         else
           adjust_frac(i,j) = 0.
         endif
       end do
      end do
      do k=1,kx
        ttnd_adjustment(:,:,k) = &
                      ((HLV*liquid_precip(is:ie,js:je,k)*adjust_frac(:,:) + &
                        HLS*frozen_precip(is:ie,js:je,k)*adjust_frac(:,:))  &
                       *dt/seconds_per_day)/CP_AIR
        liquid_precip(is:ie,js:je,k) = liquid_precip(is:ie,js:je,k) * (1.0+adjust_frac(:,:))
        frozen_precip(is:ie,js:je,k) = frozen_precip(is:ie,js:je,k) * (1.0+adjust_frac(:,:))
      end do
    else ! (force_donner_moist_conserv)
      precip_adjustment(:,:) = 0.0
      adjust_frac      (:,:) = 0.0
      ttnd_adjustment(:,:,:) = 0.
    endif  ! (force_donner_moist_conserv)

    do k=1,kx
      rain_don(:,:) = rain_don(:,:) + liquid_precip(is:ie,js:je,k)* pmass(is:ie,js:je,k)/seconds_per_day
      snow_don(:,:) = snow_don(:,:) + frozen_precip(is:ie,js:je,k)* pmass(is:ie,js:je,k)/seconds_per_day
    end do

!----------------------------------------------------------------------
!   modify each of the 3d precip fluxes returned from donner_deep, as
!   needed.
!----------------------------------------------------------------------
    if (do_cosp) then
      do k=1, size(t,3)
        do j=js,je 
          do i=is,ie 
            frz_meso(i,j,k) = frz_meso(i,j,k)*pmass(i,j,k)* &
                              scale(i-is+1,j-js+1)*(1.0+adjust_frac(i-is+1,j-js+1))/ &
                                                       SECONDS_PER_DAY
            liq_meso(i,j,k) = liq_meso(i,j,k)*pmass(i,j,k)* &
                              scale(i-is+1,j-js+1)*(1.0+adjust_frac(i-is+1,j-js+1))/ &
                                                       SECONDS_PER_DAY
            frz_cell(i,j,k) = frz_cell(i,j,k)*pmass(i,j,k)* &
                              scale(i-is+1,j-js+1)*(1.0+adjust_frac(i-is+1,j-js+1))/ &
                                                        SECONDS_PER_DAY
            liq_cell(i,j,k) = liq_cell(i,j,k)*pmass(i,j,k)* &
                              scale(i-is+1,j-js+1)*(1.0+adjust_frac(i-is+1,j-js+1))/ &
                                                        SECONDS_PER_DAY
          end do
        end do
      end do
    endif
    
    if (only_one_conv_scheme_per_column) then
      conv_calc_completed = (rain_don + snow_don) > 0.0
    endif
!---------------------------------------------------------------------
!    convert the changes in temperature, vapor specific humidity and 
!    precipitation resulting from deep convection to time tendencies 
!    of these quantities.
!---------------------------------------------------------------------
    ttnd_don(is:ie,js:je,:) = delta_temp(is:ie,js:je,:)*dtinv 
    ttnd_don(is:ie,js:je,:) = ttnd_don(is:ie,js:je,:) + ttnd_adjustment*dtinv
    qtnd_don(is:ie,js:je,:) = delta_q(is:ie,js:je,:)*dtinv

!--------------------------------------------------------------------
!    save the tendencies of temperature and specific humidity resulting
!    from the deep convection component of the donner parameterization. 
!--------------------------------------------------------------------
    ttnd_conv(is:ie,js:je,:) = ttnd_conv(is:ie,js:je,:) + ttnd_don(is:ie,js:je,:)
    qtnd_conv(is:ie,js:je,:) = qtnd_conv(is:ie,js:je,:) + qtnd_don(is:ie,js:je,:)

!--------------------------------------------------------------------
!    add the contributions to the temperature and vapor specific 
!    humidity tendencies from donner_deep mod to the arrays accumulating
!    the total tendencies due to all physics processes.
!--------------------------------------------------------------------
    tdt = tdt + ttnd_don(is:ie,js:je,:) 
    qdt = qdt + qtnd_don(is:ie,js:je,:)

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    deep convection on this step to the arrays accumulating precip-
!    itation from all sources (lprec, fprec).
!--------------------------------------------------------------------
    lprec  = lprec + rain_don
    fprec  = fprec + snow_don

!--------------------------------------------------------------------
!    output the time tendencies of temperature, vapor specific humid-
!    ity, precipitation and mass flux due to deep convection.
!--------------------------------------------------------------------
    used = send_data (id_tdt_deep_donner, ttnd_don(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
    used = send_data (id_qdt_deep_donner, qtnd_don(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
    used = send_data (id_qadt_deep_donner, delta_qa(is:ie,js:je,:)*dtinv, Time, is, js, 1, rmask=mask )
    used = send_data (id_qldt_deep_donner, delta_ql(is:ie,js:je,:)*dtinv, Time, is, js, 1, rmask=mask )
    used = send_data (id_qidt_deep_donner, delta_qi(is:ie,js:je,:)*dtinv, Time, is, js, 1, rmask=mask )
    used = send_data (id_mc_donner, mc_donner(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
    used = send_data (id_mc_donner_half, mc_donner_half(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
    used = send_data (id_m_cdet_donner, m_cdet_donner(is:ie,js:je,:), Time,  is, js, 1, rmask=mask )
    used = send_data (id_m_cellup, m_cellup, Time, is, js, 1, rmask=mask )
    used = send_data (id_snow_deep_donner, snow_don, Time, is, js)
    used = send_data (id_prec_deep_donner, rain_don + snow_don, Time, is, js )
    used = send_data (id_prec1_deep_donner, precip_adjustment,  &
                         Time, is, js, mask = precip_returned > 0.0)

    if (do_donner_conservation_checks) then
      used = send_data (id_enth_donner_col2, -hlv*rain_don, Time, is, js)
      used = send_data (id_enth_donner_col3, -hls*snow_don, Time, is, js)
      if (id_enth_donner_col4 > 0) call column_diag(id_enth_donner_col4, is, js, Time, &
                                        ttnd_don(is:ie,js:je,:), CP_AIR)
      if (id_enth_donner_col5 > 0) call column_diag(id_enth_donner_col5, is, js, Time, &
                                        delta_ql(is:ie,js:je,:), -HLV*dtinv, delta_qi(is:ie,js:je,:), -HLS*dtinv)
      if (id_enth_donner_col6 > 0) call column_diag(id_enth_donner_col6, is, js, Time, &
                                        ttnd_adjustment, CP_AIR)
      used = send_data (id_enth_donner_col7, adjust_frac, Time, is, js)
       
      temp_2d = 0.
      do k=1,kx
        temp_2d(:,:) = temp_2d(:,:)  &
             + (-HLV*liquid_precip(is:ie,js:je,k)/seconds_per_day -  &
                hls*frozen_precip(is:ie,js:je,k)/seconds_per_day  + &
                CP_AIR*ttnd_don(is:ie,js:je,k)  &
             -  (HLV*delta_ql(is:ie,js:je,k)*dtinv + HLS*delta_qi(is:ie,js:je,k)*dtinv)  &
                )*pmass(is:ie,js:je,k)
      end do
      used = send_data (id_enth_donner_col, temp_2d, Time, is, js)

      if (id_wat_donner_col > 0) then
        temp_2d = rain_don + snow_don
        call column_diag(id_wat_donner_col, is, js, Time, qtnd_don(is:ie,js:je,:), 1.0, &
                         delta_ql(is:ie,js:je,:), dtinv, delta_qi(is:ie,js:je,:), dtinv, temp_2d)
      endif
    endif ! (donner_conservation_checks)

! ---> h1g, dump donner-deep  cell and meso cloud fraction, 2010-08-08
    used = send_data (id_cell_cld_frac,  cell_cld_frac(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
    used = send_data (id_meso_cld_frac,  meso_cld_frac(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
    used = send_data (id_donner_humidity_area,  donner_humidity_area(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
! <--- h1g, dump donner-deep  cell and meso cloud fraction, 2010-08-08

    call mpp_clock_end (donner_clock)

    if (do_donner_mca) then
      call mpp_clock_begin (donner_mca_clock)
!--------------------------------------------------------------------
!    call subroutine moist_conv to handle any shallow convection 
!    present in the grid. in this call do_strat is always set to .false.
!    so that no convective detrainment (and corresponding change in
!    large-scale cloud amount and area) from moist convective adjustment
!    is allowed, consistent with this call being constrained to handle
!    shallow convection.
!--------------------------------------------------------------------
      tin(is:ie,js:je,:) = tin(is:ie,js:je,:)+delta_temp(is:ie,js:je,:)
      qin(is:ie,js:je,:) = qin(is:ie,js:je,:)+delta_q(is:ie,js:je,:)
      call moist_conv (tin(is:ie,js:je,:), qin(is:ie,js:je,:), pfull, phalf, coldT, &
                       ttnd_don(is:ie,js:je,:), qtnd_don(is:ie,js:je,:), &
                       rain_donmca, snow_donmca, dtinv, Time, is, js,     &
                       donner_tracer(is:ie,js:je,:,:), qtr(is:ie,js:je,:,:), Lbot=kbot, mask=mask)

      if (include_donmca_in_cosp) then
        do j=js,je
          do i=is,ie
            if (coldT(i-is+1,j-js+1)) then
              do k=1,kx
                mca_frz(i,j,k) = -1.0*qtnd_don(i,j,k)*pmass(i,j,k)
                mca_liq(i,j,k) = 0.
              end do
            else
              do k=1,kx
                mca_frz(i,j,k) = 0.
                mca_liq(i,j,k) = -1.0*qtnd_don(i,j,k)*pmass(i,j,k)
              end do
            endif
          end do
        end do
      else
        mca_frz(is:ie,js:je,:) = 0.
        mca_liq(is:ie,js:je,:) = 0.
      endif
!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from moist convective adjustment. currently there
!    is no tracer transport by this process.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_donner(n)) then
          rdt(:,:,:,n) = rdt(:,:,:,n) + qtr(is:ie,js:je,:,nn)
          nn = nn + 1
        endif
      end do

!--------------------------------------------------------------------
!    define the heating, moistening and precipitation rates as the sum 
!    of the contributions from the deep convection pass and the moist 
!    convective adjustment pass of the donner parameterization. if 
!    ras_mod is also activated, store these values in temporary arrays
!    until the contributions from ras_mod is calculated.
!--------------------------------------------------------------------
      ttnd_conv(is:ie,js:je,:) = ttnd_conv(is:ie,js:je,:) + ttnd_don(is:ie,js:je,:)
      qtnd_conv(is:ie,js:je,:) = qtnd_conv(is:ie,js:je,:) + qtnd_don(is:ie,js:je,:)

!--------------------------------------------------------------------
!    add the contributions to the temperature and vapor specific 
!    humidity tendencies from the moist convective adjustment pass of
!    donner_deep_mod to the arrays accumulating the total tendencies 
!    due to all physics processes.
!--------------------------------------------------------------------
      tdt = tdt + ttnd_don(is:ie,js:je,:)
      qdt = qdt + qtnd_don(is:ie,js:je,:)

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    the moist convective adjustment pass of the donner parameterization
!    on this step to the arrays accumulating precipitation from all 
!    sources (lprec, fprec).
!--------------------------------------------------------------------
      lprec  = lprec + rain_donmca
      fprec  = fprec + snow_donmca

!--------------------------------------------------------------------
!    output the time tendencies of temperature, vapor specific humid-
!    ity, precipitation and snow due to the moist convective 
!    adjustment pass of the donner parameterization.
!--------------------------------------------------------------------
      used = send_data (id_tdt_mca_donner, ttnd_don(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
      used = send_data (id_qdt_mca_donner, qtnd_don(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
      used = send_data (id_prec_mca_donner, rain_donmca+snow_donmca, Time, is, js) 
      used = send_data (id_snow_mca_donner, snow_donmca, Time, is, js)

      if (id_enth_mca_donner_col > 0) then
        temp_2d = -HLV*rain_donmca -HLS*snow_donmca
        call column_diag(id_enth_mca_donner_col, is, js, Time, ttnd_don(is:ie,js:je,:), CP_AIR, temp_2d)
      endif

      if (id_wat_mca_donner_col > 0) then
        temp_2d = rain_donmca + snow_donmca
        call column_diag(id_wat_mca_donner_col, is, js, Time, qtnd_don(is:ie,js:je,:), 1.0, temp_2d)
      endif

!--------------------------------------------------------------------
!------- diagnostics for tracers from convection -------
!  allow any tracer to be activated here (allows control cases)
!--------------------------------------------------------------------
      do n=1,num_tracers
        used = send_data ( id_conv_tracer(n), tracer(is:ie,js:je,:,n), Time, is, js, 1, &
                           rmask=mask )
!------- diagnostics for tracers column integral tendency ------
        if ( id_conv_tracer_col(n) > 0 ) &
          call column_diag(id_conv_tracer_col(n), is, js, Time, tracer(is:ie,js:je,:,n), 1.0)
      enddo

!--------------------------------------------------------------------
!    output the time tendencies of tracer and of column tracer 
!    due to the moist convective adjustment pass of the donner 
!    parameterization. currently moist convective adjustment does not
!    affect the tracer fields, so these fields are always 0.0.
!--------------------------------------------------------------------
      do n = 1, num_donner_tracers
        if ( id_tracerdt_mcadon(n) > 0 ) &
          used = send_data(id_tracerdt_mcadon(n), qtr(is:ie,js:je,:,n), Time, is, js, 1, rmask=mask )
        if (id_tracerdt_mcadon_col(n) > 0 )  &
          call column_diag(id_tracerdt_mcadon_col(n), is, js, Time, qtr(is:ie,js:je,:,n), 1.0)
      enddo

      call mpp_clock_end (donner_mca_clock)
    endif !(do_donner_mca) 

!---------------------------------------------------------------------
!    if donner_deep_mod is not active, define input fields normally 
!    produced by donner_deep_mod and needed by strat_cloud_mod 
!    appropriately.
!---------------------------------------------------------------------
  else   ! (do_donner_deep)
!    mc_donner(is:ie,js:je,:) = 0.0
!    mc_donner_up(is:ie,js:je,:) = 0.0
!    mc_donner_half(is:ie,js:je, : ) = 0.0
!    m_cdet_donner(is:ie,js:je,:) = 0.0
    m_cellup = 0.0
!    donner_humidity_area(is:ie,js:je,:) = 0.
!    donner_humidity_factor(is:ie,js:je,:) = 0.
  endif  ! (do_donner_deep)
! ADD TENDENCIES HERE, IN SAME AORDER AS ORIGINAL:
  if (do_donner_deep) then
    if (limit_conv_cloud_frac) then
      ltemp = ANY(donner_humidity_area(is:ie,js:je,:) >= 0.999, dim = 3)
      where (ltemp(:,:)) conv_calc_completed(:,:) = .true.
      available_cf_for_uw = MAX(0.999 - donner_humidity_area(is:ie,js:je,:), 0.0)
    endif

    if (do_strat) then
      tracer(is:ie,js:je,:,nql) = qlin(is:ie,js:je,:) + delta_ql(is:ie,js:je,:)
      tracer(is:ie,js:je,:,nqi) = qiin(is:ie,js:je,:) + delta_qi(is:ie,js:je,:)
      tracer(is:ie,js:je,:,nqa) = qain(is:ie,js:je,:) + delta_qa(is:ie,js:je,:)
      rdt(:,:,:,nql) = rdt(:,:,:,nql) + delta_ql(is:ie,js:je,:)*dtinv
      rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + delta_qi(is:ie,js:je,:)*dtinv
      rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + delta_qa(is:ie,js:je,:)*dtinv
    endif

    IF (do_ice_num .AND. detrain_ice_num) THEN
      CALL detr_ice_num (tin(is:ie,js:je,:), delta_qi(is:ie,js:je,:), &
                                                delta_qni(is:ie,js:je,:))   
      tracer(is:ie,js:je,:,nqni) =  nilin(is:ie,js:je,:)  +   &
                                                 delta_qni(is:ie,js:je,:) 
      rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + delta_qni(is:ie,js:je,:)*dtinv

      used = send_data (id_qnidt_deep_donner,   &
               delta_qni(is:ie,js:je,:)*dtinv, Time, is, js, 1, rmask=mask)
    END IF

!-------------------------------------------------------------------------
!    detrain liquid droplets if desired. assume 10 micron mean volume 
!    radius for detrained droplets
!-------------------------------------------------------------------------
    IF (do_liq_num .AND. detrain_liq_num) THEN
      if ( remain_detrain_bug ) then
        delta_qn(is:ie,js:je,:) =  delta_ql(is:ie,js:je,:)/1000.*   &
                                                        3./(4.*3.14*10.e-15)
      else
        delta_qn(is:ie,js:je,:) =  delta_ql(is:ie,js:je,:)/1000.*   &
                                                        3./(4.*3.14e-15)
      endif !( remain_detrain_bug )
      tracer(is:ie,js:je,:,nqn) =  nllin (is:ie,js:je,:) +  &
                                                  delta_qn(is:ie,js:je,:) 
      rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + delta_qn(is:ie,js:je,:)*dtinv

      used = send_data (id_qndt_deep_donner,   &
              delta_qn(is:ie,js:je,:)*dtinv, Time, is, js, 1, rmask=mask)
    END IF

!---------------------------------------------------------------------
!    update the values of temperature and vapor specific humidity to
!    include the effects of deep convection.
!---------------------------------------------------------------------
    if (.not. do_donner_mca) then
      tin(is:ie,js:je,:) = tin(is:ie,js:je,:) + delta_temp(is:ie,js:je,:)
      qin(is:ie,js:je,:) = qin(is:ie,js:je,:) + delta_q(is:ie,js:je,:)
    endif
  endif !(do_donner_deep)

  if (do_donner_before_uw) then
    if (do_uw_conv) then
      call mpp_clock_begin (shallowcu_clock)
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the uw_conv param-
!    eterization are present.
!---------------------------------------------------------------------
      if    &
        (present (shallow_cloud_area) .and.   &
         present (shallow_liquid) .and.   &
         present (shallow_ice) .and.  &
         present ( shallow_droplet_number) ) then
      else
        call error_mesg ('moist_processes_mod', 'moist_processes: &
               &not all 4 optional arguments needed for uw_conv &
               &output are present', FATAL)
      endif

      if (use_updated_profiles_for_uw) then 
!---------------------------------------------------------------------
!    update tracer fields with tendencies due to donner convection and 
!    wet deposition by donner deep precipitation.
!---------------------------------------------------------------------
        do n=1,size(rdt,4)
          if (n /= nsphum) then
            if (.not. do_strat .or. ( n /= nql .and. n /= nqi .and.   &
                    n /= nqa .and. n /= nqn  .and.  n /= nqni  ) ) then
              tracer(is:ie,js:je,:,n) = tracer_orig(is:ie,js:je,:,n) +   &
                             (rdt(:,:,:,n) - rdt_init(is:ie,js:je,:,n)) *dt
            endif
          endif
        end do
        call moistproc_uw_conv(Time, is, ie, js, je, dt, tin(is:ie,js:je,:), qin(is:ie,js:je,:), &
                               uin(is:ie,js:je,:), vin(is:ie,js:je,:), tracer(is:ie,js:je,:,:),    &
                               pfull, phalf, zfull, zhalf, omega, pblht,        &
                               ustar, bstar, qstar, land, coldT, Aerosol,       &
                               cush, cbmf, cmf(is:ie,js:je,:), conv_calc_completed,            &
                               available_cf_for_uw, tdt, qdt, udt, vdt, rdt,    &
                               ttnd_conv(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:), lprec, fprec, precip,      &
                               liq_precflx(is:ie,js:je,:),  &
                               ice_precflx(is:ie,js:je,:), &
                               do_strat, do_limit_uw, do_liq_num, num_tracers,  &
                               tracers_in_uw, num_uw_tracers, shallow_cloud_area,&
                               shallow_liquid, shallow_ice, shallow_droplet_number, uw_wetdep, &
                               do_ice_num, detrain_ice_num)

      else ! (.not. use_updated_profiles_for_uw)
        call moistproc_uw_conv(Time, is, ie, js, je, dt, tin_orig(is:ie,js:je,:), qin_orig(is:ie,js:je,:), &
                               uin(is:ie,js:je,:), vin(is:ie,js:je,:), tracer_orig(is:ie,js:je,:,:),    &
                               pfull, phalf, zfull, zhalf, omega, pblht,        &
                               ustar, bstar, qstar, land, coldT, Aerosol,       &
                               cush, cbmf, cmf(is:ie,js:je,:), conv_calc_completed,            &
                               available_cf_for_uw, tdt, qdt, udt, vdt, rdt,    &
                               ttnd_conv(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:), lprec, fprec, precip,      &
                               liq_precflx(is:ie,js:je,:),  &
                               ice_precflx(is:ie,js:je,:), &
                               do_strat, do_limit_uw, do_liq_num, num_tracers,  &
                               tracers_in_uw, num_uw_tracers, shallow_cloud_area,&
                               shallow_liquid, shallow_ice, shallow_droplet_number, uw_wetdep, &
                               do_ice_num, detrain_ice_num)
      endif ! (use_updated_profiles_for_uw)
      call mpp_clock_end (shallowcu_clock)
    endif !(do_uw_conv)
  endif !(do_donner_before_uw)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                B. MOIST CONVECTIVE ADJUSTMENT             
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if (do_mca) then
     call mpp_clock_begin (mca_clock)
     call moistproc_mca(Time, is, js, tin(is:ie,js:je,:), qin(is:ie,js:je,:), tracer(is:ie,js:je,:,:), pfull, phalf, coldT, dtinv, &
                        tdt, qdt, rdt, q_tnd(is:ie,js:je,:,:), ttnd_conv(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:),                 &
                        lprec, fprec, do_strat, num_tracers, tracers_in_mca,        &
                        num_mca_tracers, kbot, mask)
     call mpp_clock_end (mca_clock)
   endif ! (do_mca)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!           X. BETTS-MILLER CONVECTION SCHEME 
!			
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   if ( any((/do_bm,do_bmmass,do_bmomp/)) ) then
     call mpp_clock_begin (bm_clock)

     if (do_bm) then ! betts-miller cumulus param scheme
       call betts_miller (dt,tin(is:ie,js:je,:),qin(is:ie,js:je,:),pfull,phalf,coldT,rain,snow,&
                         ttnd(is:ie,js:je,:),qtnd(is:ie,js:je,:), &
                         q_ref(is:ie,js:je,:),bmflag,klzbs,cape,cin,t_ref(is:ie,js:je,:),invtaubmt,       &
                         invtaubmq, mask=mask)
     endif

     if (do_bmmass) then ! betts-miller-style massflux cumulus param scheme
       call bm_massflux (dt,tin(is:ie,js:je,:),qin(is:ie,js:je,:),pfull,phalf,coldT,rain,snow,&
                         ttnd(is:ie,js:je,:),qtnd(is:ie,js:je,:),  &
                         q_ref(is:ie,js:je,:),bmflag,klzbs,t_ref(is:ie,js:je,:),massflux(is:ie,js:je,:), mask=mask)
     endif

     if (do_bmomp) then ! olivier's betts-miller cumulus param scheme
       call bm_omp (dt,tin(is:ie,js:je,:),qin(is:ie,js:je,:),pfull,phalf,coldT,rain,snow,&
                    ttnd(is:ie,js:je,:),qtnd(is:ie,js:je,:),       &
                    q_ref(is:ie,js:je,:),bmflag,klzbs,t_ref(is:ie,js:je,:), mask=mask)
     endif

!------- (update input values and) compute tendency -----
     tin(is:ie,js:je,:)=tin(is:ie,js:je,:)+ttnd(is:ie,js:je,:)
     qin(is:ie,js:je,:)=qin(is:ie,js:je,:)+qtnd(is:ie,js:je,:)
     ttnd(is:ie,js:je,:)=ttnd(is:ie,js:je,:)*dtinv
     qtnd(is:ie,js:je,:)=qtnd(is:ie,js:je,:)*dtinv
     rain=rain*dtinv
     snow=snow*dtinv
                                                                                   
!-------- add on tendency ----------
     tdt=tdt+ttnd(is:ie,js:je,:) 
     qdt=qdt+qtnd(is:ie,js:je,:)

!------- save total precip and snow ---------
     lprec=lprec+rain
     fprec=fprec+snow
     precip=precip+rain+snow
                                                                         
!------- compute rh clouds if desired ------
     if (do_rh_clouds) then
!calculate relative humidity
       call rh_calc(pfull,tin(is:ie,js:je,:),qin(is:ie,js:je,:),RH(is:ie,js:je,:),do_simple,mask)
!pass RH to rh_clouds_sum
       call rh_clouds_sum (is, js, RH(is:ie,js:je,:)) ! XXX  RH is not relative humidity when do_simple=.true.
     end if

! betts-miller diags
     used = send_data (id_tref, t_ref(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
     used = send_data (id_qref, q_ref(is:ie,js:je,:), Time, is, js, 1, rmask=mask )
     used = send_data (id_bmflag, bmflag, Time, is, js)
     used = send_data (id_klzbs, klzbs, Time, is, js)
     used = send_data (id_invtaubmt, invtaubmt, Time, is, js)
     used = send_data (id_invtaubmq, invtaubmq, Time, is, js)
     used = send_data (id_massflux, massflux(is:ie,js:je,:), Time, is, js, 1, rmask=mask)

     call mpp_clock_end (bm_clock)
   endif ! if ( any((/do_bm,do_bmmass,do_bmomp/)) )


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!           C. RELAXED ARAKAWA-SCHUBERT PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!    execute relaxed arakawa/schubert cumulus parameterization scheme,
!    if desired.
!-----------------------------------------------------------------------
   if (do_ras) then
     call mpp_clock_begin (ras_clock)
     call moistproc_ras(Time, is, js, dt, coldT, tin(is:ie,js:je,:), qin(is:ie,js:je,:), uin(is:ie,js:je,:), vin(is:ie,js:je,:), &
                        tracer(is:ie,js:je,:,:), pfull, phalf, zhalf, tdt, qdt, udt, vdt, rdt,       &
                        q_tnd(is:ie,js:je,:,:), ttnd(is:ie,js:je,:), qtnd(is:ie,js:je,:), &
                        ttnd_conv(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:), mc, det0(is:ie,js:je,:),  &
                        lprec, fprec, rain_ras, snow_ras, rain3d, snow3d,   &
                        Aerosol, do_strat, do_liq_num, num_tracers,         &
                        tracers_in_ras, num_ras_tracers, kbot, mask, &
                        do_ice_num, detrain_ice_num)
     call mpp_clock_end (ras_clock)
   else
!---------------------------------------------------------------------
!    if ras_mod is not activated, set the ras mass flux field to be 0.0.
!---------------------------------------------------------------------
     mc   = 0.0
!     det0(is:ie,js:je,:) = 0.0
     rain_ras = 0.0
     snow_ras = 0.0
   endif  ! (do_ras)

!---------------------------------------------------------------------
!    call subroutine cu_mo_trans if diffusive cumulus momentum 
!    transport is desired. 
!---------------------------------------------------------------------
   if (do_cmt) then
     call mpp_clock_begin (cmt_clock)
     diff_cu_mo(:,:,:)  = 0.0

!  if doing nonlocal cmt, call cu_mo_trans for each convective scheme
!  separately
     if (.not. doing_diffusive) then
       if (cmt_uses_ras) then
!        mc_cmt = mc
!        det_cmt = det0
         call moistproc_cmt ( Time, is, js, tin(is:ie,js:je,:), uin(is:ie,js:je,:), vin(is:ie,js:je,:), &
                              tracer(is:ie,js:je,:,:), pfull, phalf, &
                              zfull, zhalf, pmass(is:ie,js:je,:), tdt, udt, vdt, rdt,           &
                              ttnd_conv(is:ie,js:je,:), dt, mc, det0(is:ie,js:je,:), diff_cu_mo,               &
                              num_tracers)
       endif !(cmt_uses_ras)
!
       if (cmt_uses_donner) then
!        mc_cmt = m_cellup 
!        det_cmt = m_cdet_donner 
         call moistproc_cmt ( Time, is, js, tin(is:ie,js:je,:), uin(is:ie,js:je,:), vin(is:ie,js:je,:), &
                              tracer(is:ie,js:je,:,:), pfull, phalf, &
                              zfull, zhalf, pmass(is:ie,js:je,:), tdt, udt, vdt, rdt,           &
                              ttnd_conv(is:ie,js:je,:), dt, m_cellup, M_cdet_donner(is:ie,js:je,:), diff_cu_mo,&
                              num_tracers)
       endif
!
       if (cmt_uses_uw) then
         mc_cmt(:,:,1) = 0.
         mc_cmt(:,:,kx+1) = 0.
         do k=2,kx
           mc_cmt(:,:,k) = cmf(is:ie,js:je,k-1)
         end do
!   CURRENTLY no detrained mass flux provided from uw_conv; should only
!   use with 'diffusive' cmt scheme, not the non-local. (attempt to
!   use non-local will cause FATAL in _init routine.)
!         det_cmt(is:ie,js:je,:) = 0.0   
         call moistproc_cmt ( Time, is, js, tin(is:ie,js:je,:), uin(is:ie,js:je,:), vin(is:ie,js:je,:), &
                              tracer(is:ie,js:je,:,:), pfull, phalf, &
                              zfull, zhalf, pmass(is:ie,js:je,:), tdt, udt, vdt, rdt,           &
                              ttnd_conv(is:ie,js:je,:), dt, mc_cmt, det_cmt(is:ie,js:je,:), diff_cu_mo,        &
                              num_tracers)
       endif

     else ! (we are doing_diffusive)

!  if using diffusive cmt, call cu_mo_trans once with combined mass
!  fluxes from all desired convective schemes.
       mc_cmt = 0.
!       det_cmt(is:ie,js:je,:) = 0.
       if (cmt_uses_ras) then
         mc_cmt = mc_cmt + mc
       endif
       if (cmt_uses_donner) then
         mc_cmt = mc_cmt + m_cellup 
       endif
       if (cmt_uses_uw) then
         do k=2,kx
           mc_cmt(:,:,k) = mc_cmt(:,:,k) + cmf(is:ie,js:je,k-1)
         end do
       endif
       call moistproc_cmt ( Time, is, js, tin(is:ie,js:je,:), uin(is:ie,js:je,:), vin(is:ie,js:je,:), &
                            tracer(is:ie,js:je,:,:), pfull, phalf, &
                            zfull, zhalf, pmass(is:ie,js:je,:), tdt, udt, vdt, rdt,           &
                            ttnd_conv(is:ie,js:je,:), dt, mc_cmt, det_cmt(is:ie,js:je,:), diff_cu_mo,        &
                            num_tracers)
     endif ! (.not. doing_diffusive)
     call mpp_clock_end (cmt_clock)
   else  !(do_cmt)
     diff_cu_mo(:,:,:)  = 0.0
   endif  ! (do_cmt)

!---------------------------------------------------------------------
!    calculate the tracer tendency due to wet deposition (wetdeptnd)
!    caused by the convectively generated precipitation (rain, snow) for
!    any tracers for which wet deposition has been activated. add this 
!    tendency to the tracer tendency due to all physics (rdt). save it 
!    also in an array which will be combined with any wet deposition 
!    resulting from large-scale precip producing the total wet deposition
!    for the tracer (wet_data).
!---------------------------------------------------------------------
   f_snow_berg = 0.
   wet_data = 0.0
   qtnd_wet(is:ie,js:je,:) = qtnd(is:ie,js:je,:)
   if (do_strat) then
     qtnd_wet(is:ie,js:je,:) = qtnd_wet(is:ie,js:je,:) + q_tnd(is:ie,js:je,:,nql) + q_tnd(is:ie,js:je,:,nqi)
     cloud_wet(is:ie,js:je,:) = 1.e-3
   else
     cloud_wet(is:ie,js:je,:) = 1.e-3
   end if
   cloud_frac(is:ie,js:je,:) = 0.1
    do n=1,size(rdt,4)
     if ( n /= nsphum ) then
       if ( .not. do_strat .or. (n /= nql .and. n /= nqi .and.   &
                        n /= nqa .and. n /= nqn  .and. n /= nqni  ) ) then
         wetdeptnd(is:ie,js:je,:) = 0.0
         call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain_ras, snow_ras, &
                              qtnd_wet(is:ie,js:je,:), cloud_wet(is:ie,js:je,:), cloud_frac(is:ie,js:je,:),                      &
                              f_snow_berg, &
                              rain3d, snow3d, tracer(is:ie,js:je,:,n), wetdeptnd(is:ie,js:je,:),           &
                              Time, 'convect', is, js, dt )
         rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd(is:ie,js:je,:)
         wet_data(:,:,:,n) = wetdeptnd(is:ie,js:je,:)
       endif
     endif  
   end do

   mc_full(is:ie,js:je,1)=0.; 
   mc_half(is:ie,js:je,1)=0.; 
   do k=2,kx   
     mc_full(is:ie,js:je,k) = 0.5*(mc(:,:,k) + mc(:,:,k+1)) +   &
                      0.5*(cmf(is:ie,js:je,k)+cmf(is:ie,js:je,k-1)) +   &
                           mc_donner(is:ie,js:je,k)
   end do
   do k=2,kx+1   
     mc_half(is:ie,js:je,k) = mc(:,:,k) +    &
                      cmf(is:ie,js:je,k-1)+   &
                           mc_donner_half(is:ie,js:je,k)
   end do

   if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER &
       .or. id_conv_freq > 0 &
       .or. id_conv_cld_base > 0 &
       .or. id_conv_cld_top > 0 ) then

     cldbot = 0
     cldtop = 0
     do j = 1,jx
       do i = 1,ix
         do k = 1,kx
           if (mc_full(i+is-1,j+js-1,k) /= 0 ) then
             cldtop(i,j) = k
             exit
           endif
         enddo
         do k = size(r,3),1,-1
           if (mc_full(i+is-1,j+js-1,k) /= 0 ) then
             cldbot(i,j) = k
             exit
           endif
         enddo
       enddo
     enddo
   end if

   if ( id_conv_cld_base > 0 ) then
     temp_2d = missing_value
     do j = 1,jx
       do i = 1,ix
         if ( cldbot(i,j) > 0 ) temp_2d(i,j) = pfull(i,j,cldbot(i,j))
       end do
     end do
     used = send_data(id_conv_cld_base, temp_2d, Time, is_in=is,   &
                                           js_in=js,  mask = cldbot > 0)
   end if

   if ( id_conv_cld_top > 0 ) then
     temp_2d = missing_value
     do j = 1,jx
       do i = 1,ix
         if ( cldtop(i,j) > 0 ) temp_2d(i,j) = pfull(i,j,cldtop(i,j))
       end do
     end do
     used = send_data(id_conv_cld_top, temp_2d, Time, is_in=is, &
                                  js_in=js,  mask = cldtop > 0)
   end if

!-----------------------------------------------------------------------
! lightning NOx parameterization
!-----------------------------------------------------------------------
   if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER ) then
!     cldbot = 0
!     cldtop = 0
!     do i = 1,ix
!      do j = 1,jx
!       do k = 1,kx
!         if (mc_full(i+is-1,j+js-1,k) /= 0 ) then
!           cldtop(i,j) = k
!           exit
!         endif
!       enddo
!       do k = size(r,3),1,-1
!         if (mc_full(i+is-1,j+js-1,k) /= 0 ) then
!           cldbot(i,j) = k
!           exit
!         endif
!       enddo
!      enddo
!     enddo
     call moz_hook(cldtop, cldbot, land, zfull, zhalf, t, prod_no, area, lat, &
                   Time, is, js)
     rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) =  &
              rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) + &
              prod_no* ((boltz * t) / (10. * pfull))                     !  conc_air
     used = send_data(id_prod_no,prod_no, Time, is_in=is, js_in=js)
   endif

!-----------------------------------------------------------------------
!    define the total precipitation rate (precip).
!-----------------------------------------------------------------------
   precip = lprec + fprec

!-----------------------------------------------------------------------
!    calculate convective gustiness, if desired.
!-----------------------------------------------------------------------
   if (do_gust_cv) then
     where((precip) > 0.0)
       gust_cv = gustmax*sqrt( precip/(gustconst + precip) )
     end where
   end if

!---------------------------------------------------------------------
!    save a diagnostic indicating whether or not convection has occurred
!    within the column.
!---------------------------------------------------------------------
   where (precip > 0.) convect = .true.

!---------------------------------------------------------------------
!    apply changes resulting from uw_conv
!---------------------------------------------------------------------
   if (do_uw_conv) then
     if (do_limit_uw) then
       call moistproc_scale_uw(is,ie,js,je,dt, qin(is:ie,js:je,:), tracer(is:ie,js:je,:,:), tdt, qdt, udt, vdt, rdt,  &
                               ttnd_conv(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:), lprec, fprec, precip,&
                               do_strat, do_liq_num, num_tracers,         &
                               tracers_in_uw, scale, do_ice_num)
       used = send_data (id_scale_uw, scale, Time, is, js )
     else !(do_limit_uw) 
        scale = 1.0
        used = send_data (id_scale_uw, scale, Time, is, js )
     endif !(do_limit_uw)

!       update input fields with changes from uw_conv
     tin(is:ie,js:je,:) = tin(is:ie,js:je,:) + ttnd_uw(is:ie,js:je,:)*dt
     qin(is:ie,js:je,:) = qin(is:ie,js:je,:) + qtnd_uw(is:ie,js:je,:)*dt
     uin(is:ie,js:je,:) = uin(is:ie,js:je,:) + utnd_uw(is:ie,js:je,:)*dt
     vin(is:ie,js:je,:) = vin(is:ie,js:je,:) + vtnd_uw(is:ie,js:je,:)*dt
     tracer(is:ie,js:je,:,nql) = tracer(is:ie,js:je,:,nql) + qltnd_uw(is:ie,js:je,:)*dt
     tracer(is:ie,js:je,:,nqi) = tracer(is:ie,js:je,:,nqi) + qitnd_uw(is:ie,js:je,:)*dt
     tracer(is:ie,js:je,:,nqa) = tracer(is:ie,js:je,:,nqa) + qatnd_uw(is:ie,js:je,:)*dt
     if (do_liq_num) then
       tracer(is:ie,js:je,:,nqn) = tracer(is:ie,js:je,:,nqn) + qntnd_uw(is:ie,js:je,:)*dt
     endif
     if (do_ice_num) then
       tracer(is:ie,js:je,:,nqni) = tracer(is:ie,js:je,:,nqni) +   &
                                               qnitnd_uw(is:ie,js:je,:)*dt
     endif
   endif !(uw_conv)
 
!---------------------------------------------------------------------
!    update tracer fields with tendencies due to convection and wet 
!    deposition by convective precipitation.
!---------------------------------------------------------------------
   do n=1,size(rdt,4)
     if (n /= nsphum) then
       if (.not. do_strat .or. ( n /= nql .and. n /= nqi .and.   &
                        n /= nqa .and. n /= nqn .and. n /= nqni ) ) then
!        tracer(:,:,:,n) = tracer(:,:,:,n) +   &
         tracer(is:ie,js:je,:,n) = tracer_orig(is:ie,js:je,:,n) +   &
                           (rdt(:,:,:,n) - rdt_init(is:ie,js:je,:,n)) *dt
       endif
     endif
   end do

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   CONVECTION DIAGNOSTICS      
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   used = send_data (id_ras_precip, rain_ras+snow_ras, Time, is, js)
   used = send_data (id_don_precip, rain_don+snow_don+rain_donmca+snow_donmca, &
                     Time, is, js)
! uw_conv diags
   if ( id_uw_precip > 0 ) then
     used = send_data (id_uw_precip, rain_uw(is:ie,js:je) + snow_uw(is:ie,js:je), Time, is, js)
   endif
   used = send_data (id_uw_snow, snow_uw(is:ie,js:je), Time, is, js)
 if (do_uw_conv) then
   used = send_data (id_tdt_uw, ttnd_uw(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
   used = send_data (id_qdt_uw, qtnd_uw(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
   used = send_data (id_qadt_uw, qatnd_uw(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
   used = send_data (id_qldt_uw, qltnd_uw(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
   used = send_data (id_qidt_uw, qitnd_uw(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
   if (do_liq_num) then
   used = send_data (id_qndt_uw, qntnd_uw(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
   endif
   if (do_ice_num) then
     used = send_data (id_qnidt_uw, qnitnd_uw(is:ie,js:je,:),   &
                                              Time, is, js, 1, rmask=mask)
   end if
 endif
        
   if (id_ras_freq > 0) then
     ltemp = rain_ras > 0. .or. snow_ras > 0.0
     temp_2d = 0.
     where (ltemp) 
       temp_2d = 1.
     end where
     used = send_data (id_ras_freq, temp_2d,Time, is, js)
   endif

   if (id_don_freq > 0) then
     ltemp = rain_don > 0. .or. snow_don > 0.0 .or. &
                rain_donmca > 0. .or. snow_donmca > 0.0
     temp_2d = 0.
     where (ltemp) 
       temp_2d = 1.
     end where
     used = send_data (id_don_freq, temp_2d, Time, is, js)
   endif

   if (id_uw_freq > 0) then
     ltemp = rain_uw(is:ie,js:je) > 0. .or. snow_uw(is:ie,js:je) > 0.0
     temp_2d = 0.
     where (ltemp) 
       temp_2d = 1.
     end where
     used = send_data (id_uw_freq, temp_2d, Time, is, js)
   endif

   if (id_enth_uw_col > 0) then
     temp_2d = -HLV*rain_uw(is:ie,js:je) -HLS*snow_uw(is:ie,js:je)
     call column_diag(id_enth_uw_col, is, js, Time, ttnd_uw(is:ie,js:je,:), CP_AIR, qltnd_uw(is:ie,js:je,:), -HLV, &
                      qitnd_uw(is:ie,js:je,:), -HLS, temp_2d)
   endif

   if (id_wat_uw_col > 0) then
     temp_2d = rain_uw(is:ie,js:je) + snow_uw(is:ie,js:je)
     call column_diag(id_wat_uw_col, is, js, Time, qtnd_uw(is:ie,js:je,:), 1.0, qltnd_uw(is:ie,js:je,:), 1.0, &
                      qitnd_uw(is:ie,js:je,:), 1.0, temp_2d)
   endif
        
!---------------------------------------------------------------------
!    temperature change due to dry and moist convection:
!---------------------------------------------------------------------
   used = send_data (id_tdt_conv, ttnd_conv(is:ie,js:je,:), Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    vapor specific humidity change due to convection:
!---------------------------------------------------------------------
   used = send_data (id_qdt_conv, qtnd_conv(is:ie,js:je,:), Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    total precipitation due to convection:
!---------------------------------------------------------------------
   used = send_data (id_prec_conv, precip, Time, is, js)

!---------------------------------------------------------------------
!    frozen precipitation (snow) due to convection:
!---------------------------------------------------------------------
   used = send_data (id_snow_conv, fprec, Time, is, js)

!---------------------------------------------------------------------
!    convective frequency
!---------------------------------------------------------------------
   if (id_conv_freq > 0) then
     ltemp = precip > 0. .or. cldtop > 0
     where (ltemp)
       freq_count = 1.
     elsewhere
       freq_count = 0.
     end where
     used = send_data (id_conv_freq, freq_count, Time, is, js )
   endif

!---------------------------------------------------------------------
!------- diagnostics for 3D precip_conv -------
!---------------------------------------------------------------------
   used = send_data ( id_conv_rain3d, rain3d, Time, is, js, 1 )

!---------------------------------------------------------------------
!------- diagnostics for 3D snow_conv -------
!---------------------------------------------------------------------
   used = send_data ( id_conv_snow3d, snow3d, Time, is, js, 1 )

!---------------------------------------------------------------------
!    surface wind gustiness due to convection:
!---------------------------------------------------------------------
   used = send_data (id_gust_conv, gust_cv, Time, is, js)

!---------------------------------------------------------------------
!    water vapor path tendency due to convection:
!---------------------------------------------------------------------
   if (id_q_conv_col > 0) call column_diag(id_q_conv_col, is, js, Time, qtnd_conv(is:ie,js:je,:), 1.0)
   
!---------------------------------------------------------------------
!    dry static energy tendency due to dry and moist convection:
!---------------------------------------------------------------------
   if (id_t_conv_col > 0) call column_diag(id_t_conv_col, is, js, Time, ttnd_conv(is:ie,js:je,:), CP_AIR)
   
!---------------------------------------------------------------------
!    cloud liquid, ice and area tendencies due to convection:
!---------------------------------------------------------------------
   if (do_strat) then

!---------------------------------------------------------------------
!    if cloud liquid diagnostics requested:
!---------------------------------------------------------------------
     if (id_qldt_conv > 0 .or. id_ql_conv_col > 0) then
       temp_3d1 = rdt(:,:,:,nql) - rdt_init(is:ie,js:je,:,nql)

!---------------------------------------------------------------------
!    cloud liquid tendency due to convection:
!---------------------------------------------------------------------
       used = send_data (id_qldt_conv, temp_3d1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud liquid water path tendency due to convection:
!---------------------------------------------------------------------
       if (id_ql_conv_col > 0) call column_diag(id_ql_conv_col, is, js, Time, temp_3d1, 1.0)
     endif

!---------------------------------------------------------------------
!    if cloud drop diagnostics requested:
!---------------------------------------------------------------------
     if (id_qndt_conv > 0 .or. id_qn_conv_col > 0) then
       temp_3d1 = rdt(:,:,:,nqn) - rdt_init(is:ie,js:je,:,nqn)

!---------------------------------------------------------------------
!    cloud drop tendency due to convection:
!---------------------------------------------------------------------
       used = send_data (id_qndt_conv, temp_3d1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud drop water path tendency due to convection:
!---------------------------------------------------------------------
       if (id_qn_conv_col > 0) call column_diag(id_qn_conv_col, is, js, Time, temp_3d1, 1.0)
     endif

!---------------------------------------------------------------------
!    if cloud ice diagnostics requested:
!---------------------------------------------------------------------
     if (id_qidt_conv > 0 .or. id_qi_conv_col > 0) then
       temp_3d1 = rdt(:,:,:,nqi) - rdt_init(is:ie,js:je,:,nqi)

!---------------------------------------------------------------------
!    cloud ice tendency due to convection:
!---------------------------------------------------------------------
       used = send_data (id_qidt_conv, temp_3d1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud ice water path tendency due to convection:
!---------------------------------------------------------------------
       if (id_qi_conv_col > 0) call column_diag(id_qi_conv_col, is, js, Time, temp_3d1, 1.0)
     endif        

!---------------------------------------------------------------------
!    if cloud ice number diagnostics requested:
!---------------------------------------------------------------------
     if (do_ice_num .and.    &
                       (id_qnidt_conv > 0 .or. id_qni_conv_col > 0)) then
       temp_3d1 = rdt(:,:,:,nqni) - rdt_init(is:ie,js:je,:,nqni)

!---------------------------------------------------------------------
!    cloud ice number tendency due to convection:
!---------------------------------------------------------------------
       used = send_data (id_qnidt_conv, temp_3d1, Time,   &
                                                   is, js, 1, rmask=mask)
 
!---------------------------------------------------------------------
!    cloud column ice number tendency due to convection:
!---------------------------------------------------------------------
       if (id_qni_conv_col > 0)   &
           call column_diag(id_qni_conv_col, is, js, Time, temp_3d1, 1.0)
     endif

!---------------------------------------------------------------------
!    if cloud area diagnostics requested:
!---------------------------------------------------------------------
     if (id_qadt_conv > 0 .or.  id_qa_conv_col > 0 ) then
       temp_3d1 = rdt(:,:,:,nqa) - rdt_init(is:ie,js:je,:,nqa)

!---------------------------------------------------------------------
!    cloud area tendency due to convection:
!---------------------------------------------------------------------
       used = send_data (id_qadt_conv, temp_3d1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    column integrated cloud mass tendency due to convection:
!---------------------------------------------------------------------
       if (id_qa_conv_col > 0) call column_diag(id_qa_conv_col, is, js, Time, temp_3d1, 1.0)
     endif
   endif !(do_strat)
         
!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    convection parameterization:
!---------------------------------------------------------------------
   if (id_enth_conv_col > 0 .or. id_wat_conv_col > 0) then
     temp_3d1 = rdt(:,:,:,nql) - rdt_init(is:ie,js:je,:,nql)
     temp_3d2 = rdt(:,:,:,nqi) - rdt_init(is:ie,js:je,:,nqi)

     if (id_enth_conv_col > 0) then
       temp_2d = -HLV*precip -HLF*fprec
       call column_diag(id_enth_conv_col, is, js, Time, ttnd_conv(is:ie,js:je,:), CP_AIR, temp_3d1, -HLV, temp_3d2, -HLS, temp_2d)
     endif

     if (id_wat_conv_col > 0) then
       temp_2d = precip
       call column_diag(id_wat_conv_col, is, js, Time, qtnd_conv(is:ie,js:je,:), 1.0, temp_3d1, 1.0, temp_3d2, 1.0, temp_2d)
     endif
   endif

!---------------------------------------------------------------------
!    tracer tendencies due to convection:
!---------------------------------------------------------------------
   do n=1,size(rdt,4)
     if (tracers_in_donner(n) .or.  tracers_in_ras(n) .or.  &
         tracers_in_mca(n)    .or.  tracers_in_uw(n))    then

       if (id_tracerdt_conv(n) > 0 .or. id_tracerdt_conv_col(n) > 0) then
         temp_3d1 = rdt(:,:,:,n) - rdt_init(is:ie,js:je,:,n)
         used = send_data (id_tracerdt_conv(n), temp_3d1, Time, is, js, 1, rmask=mask )

!---------------------------------------------------------------------
!    tracer column tendencies due to convection:
!---------------------------------------------------------------------
         if (id_tracerdt_conv_col(n) > 0) &
           call column_diag(id_tracerdt_conv_col(n), is, js, Time, temp_3d1, 1.0)
       endif        
     endif
   end do

!---------------------------------------------------------------------
!    total convective updraft mass flux (uw + donner cell up + 
!    donner meso up
!---------------------------------------------------------------------
    if (id_mc_conv_up > 0 ) &
      used = send_data (id_mc_conv_up, cmf(is:ie,js:je,:) + &
                  mc_donner_up(is:ie,js:je,:), Time, is, js, 1, rmask=mask )

!---------------------------------------------------------------------
!    end the timing of the convection code section.
!---------------------------------------------------------------------
   call mpp_clock_end (convection_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!              LARGE-SCALE CONDENSATION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!---------------------------------------------------------------------
!    begin the timing of the large-scale condensation code section.
!---------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!         A. NON-PROGNOSTIC CONDENSATION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!    if a non-prognostic cloud scheme is active, then call lscale_cond 
!    to calculate the temperature and specific humidity tendencies 
!    related to the latent heat release associated with the large-scale 
!    supersaturation.
!-----------------------------------------------------------------------
    call mpp_clock_begin (largescale_clock)

! ---> h1g, if CLUBB is in moist-processes, CLUBB is called after donner-deep but before microphysics
    if( do_clubb == 2) then
! ---> modify tendencies to consider donner-deep effects
           call compute_qs (t, pfull, qsat)

           do k=1, kx
            do j=js, je
             do i=is, ie
              qrf = MAX (q(i,j,k), 0.0)

              if (do_uw_conv .and. do_donner_deep) then
                 conv_frac_clubb(i,j,k) = donner_humidity_area(i,j,k) +  &
                           shallow_cloud_area(i,j,k)
                 env_qv = qrf - qsat(i,j,k)*(cell_cld_frac(i,j,k) +   &
                           donner_humidity_factor(i,j,k) + shallow_cloud_area(i,j,k))
              elseif (do_donner_deep) then
                 conv_frac_clubb(i,j,k) = donner_humidity_area(i,j,k)
                 env_qv = qrf - qsat(i,j,k)*(cell_cld_frac(i,j,k) + donner_humidity_factor(i,j,k))
              elseif (do_uw_conv) then
                 conv_frac_clubb(i,j,k) = shallow_cloud_area(i,j,k)
                 env_qv = qrf -  shallow_cloud_area(i,j,k)*qsat(i,j,k)
              else
                 conv_frac_clubb(i,j,k) = 0.0
                 env_qv = qrf
              endif
              conv_frac_clubb(i,j,k) = min( conv_frac_clubb(i,j,k), conv_frac_max )
              env_fraction = 1.0 - conv_frac_clubb(i,j,k)
! <--- h1g, 2011-08-19

!---------------------------------------------------------------------
!    define the ratio of the grid-box relative humidity to the humidity
!    in the environment of the convective clouds.
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    grid box has vapor and there is vapor outside of the convective a
!    clouds available for condensation.
!----------------------------------------------------------------
              if (qrf /= 0.0 .and. env_qv > 0.0) then
 
!--------------------------------------------------------------------
!    there is grid box area not filled with convective clouds
!--------------------------------------------------------------------  
                 if (env_fraction > 0.0) then
                   convective_humidity_ratio_clubb(i,j,k) =    &
                      MAX (qrf*env_fraction/env_qv, 1.0)
 
!---------------------------------------------------------------------
!    grid box is filled with convective clouds.
!----------------------------------------------------------------------
                 else
                   convective_humidity_ratio_clubb(i,j,k) = -10.0
                 endif

!--------------------------------------------------------------------
!    either no vapor or all vapor taken up in convective clouds so 
!    none left for large-scale cd.
!---------------------------------------------------------------------
              else
                 convective_humidity_ratio_clubb(i,j,k) = 1.0
              endif
             end do
            end do
           end do

           if ( .not. use_updated_profiles_for_clubb ) then
             tdt = tdt - ttnd_conv
             qdt = qdt - qtnd_conv
             qldt_conv =  rdt(:,:,:,nql) - rdt_init(:,:,:,nql)
             qidt_conv =  rdt(:,:,:,nqi) - rdt_init(:,:,:,nqi)
             qadt_conv =  rdt(:,:,:,nqa) - rdt_init(:,:,:,nqa)
             IF (do_liq_num) qndt_conv =  rdt(:,:,:,nqn) - rdt_init(:,:,:,nqn)
             IF (do_ice_num) qnidt_conv=  rdt(:,:,:,nqni)- rdt_init(:,:,:,nqni)

             rdt(:,:,:,nql) = rdt(:,:,:,nql) - qldt_conv
             rdt(:,:,:,nqi) = rdt(:,:,:,nqi) - qidt_conv
             rdt(:,:,:,nqa) = rdt(:,:,:,nqa) - qadt_conv
             IF (do_liq_num) rdt(:,:,:,nqn)  = rdt(:,:,:,nqn)  - qndt_conv
             IF (do_ice_num) rdt(:,:,:,nqni) = rdt(:,:,:,nqni) - qnidt_conv
           endif

           call clubb( is, ie, js, je, lon, lat,                &
                       Time,                                    &
                       dt,                                      &
                       phalf, pfull, zhalf, zfull, omega,       &
                       t, q, r, u, v,                           &
                       ustar, bstar, qstar,                     &
                       tdt, qdt, rdt, udt, vdt,                 &
                       dcond_ls_liquid, dcond_ls_ice,           &
                       Ndrop_act_CLUBB, Icedrop_act_CLUBB,      &
                       ndust, rbar_dust,                        &
                       diff_t_clubb,                            &
                       qcvar_clubb=qcvar_clubb,                 &
                       tdt_shf = tdt_shf,                       &
                       qdt_lhf = qdt_lhf,                       &
                       Aerosol=Aerosol, mask=mask,              &
                       mc_full=mc_full,                         &
                       conv_frac_clubb=conv_frac_clubb,         &
                       convective_humidity_ratio_clubb=convective_humidity_ratio_clubb)

           if ( .not. use_updated_profiles_for_clubb ) then
             tdt = tdt + ttnd_conv
             qdt = qdt + qtnd_conv
             rdt(:,:,:,nql) = rdt(:,:,:,nql) + qldt_conv
             rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qidt_conv
             rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qadt_conv
             IF (do_liq_num) rdt(:,:,:,nqn)  = rdt(:,:,:,nqn)  + qndt_conv
             IF (do_ice_num) rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + qnidt_conv
           endif

! ---> updated profiles are used to call microphysics
           tin(is:ie,js:je,:) = t + tdt*dt
           qin(is:ie,js:je,:) = q + qdt*dt
           uin(is:ie,js:je,:) = u + udt*dt
           vin(is:ie,js:je,:) = v + vdt*dt
           do tr=1,size(rdt,4)
              tracer(is:ie,js:je,:,tr) = r(:,:,:,tr) + rdt(:,:,:,tr)*dt
           end do
           do tr=size(rdt,4) +1, size(r,4)
              tracer(is:ie,js:je,:,tr) = r(:,:,:,tr)
           end do
     endif  ! if do_clubb == 2
! <--- h1g, if CLUBB is in moist-processes, CLUBB is called after donner-deep but before microphysics

!zero out arrays for large scale precipitation
    rain   = 0.
    snow   = 0.
    rain3d = 0.
    snow3d = 0.
    snowclr3d = 0.
    ttnd(is:ie,js:je,:)   = 0.
    qtnd(is:ie,js:je,:)   = 0.

    if (do_lsc) then
      call mpp_clock_begin (lscalecond_clock)
      call moistproc_lscale_cond (is, js, tin(is:ie,js:je,:), qin(is:ie,js:je,:), pfull, phalf, tdt, qdt, &
                                  ttnd(is:ie,js:je,:), qtnd(is:ie,js:je,:), qtnd_conv(is:ie,js:je,:), lprec, fprec, precip,    &
                                  rain, snow, dtinv, omega, do_rh_clouds, do_simple,&
                                  do_diag_clouds, coldT, kbot=kbot, mask=mask)
      call mpp_clock_end (lscalecond_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!       B. TIEDTKE / ROTSTAYN / KLEIN PROGNOSTIC CLOUD SCHEME  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    else if (do_strat) then
      call mpp_clock_begin (stratcloud_clock)
      call moistproc_strat_cloud(Time, is, ie, js, je, lon, lat, ktop, dt, tm, tin(is:ie,js:je,:), qin(is:ie,js:je,:), &
                                 tracer(is:ie,js:je,:,:), pfull, phalf, zhalf, omega, radturbten, mc_full(is:ie,js:je,:), &
                                 diff_t, land, area, tdt, qdt, rdt, q_tnd(is:ie,js:je,:,:), ttnd(is:ie,js:je,:),  &
                                 qtnd(is:ie,js:je,:), lprec, fprec, &
                                 f_snow_berg, rain, snow, rain3d, snow3d,  &
                                  snowclr3d, &
                                 Aerosol, lsc_cloud_area, lsc_liquid, lsc_ice,    &
                                 lsc_droplet_number, donner_humidity_area(is:ie,js:je,:),        &
                                 donner_humidity_factor(is:ie,js:je,:), shallow_cloud_area,      &
                                 cell_cld_frac, meso_cld_frac,                    &
                                 do_uw_conv, do_donner_deep, do_liq_num,          &
                                 do_clubb,                                        &  ! cjg
                                 do_lin_cld_microphys, id_qvout, id_qlout,        &
                                 id_qaout, id_qiout, id_qnout, id_qniout, &
                                 limit_conv_cloud_frac, mask, &
                                 hydrostatic, phys_hydrostatic, &     
                                 zfull,                 &
                                 do_ice_num , lsc_ice_number,          &
                                lsc_snow, lsc_rain, lsc_snow_size,  &
                                 lsc_rain_size, do_legacy_strat_cloud, &
                                 dcond_ls_liquid, dcond_ls_ice,                   &           ! h1g, cjg
                                 Ndrop_act_CLUBB, Icedrop_act_CLUBB,              &           ! h1g, cjg
                                 ndust, rbar_dust, qcvar_clubb=qcvar_clubb)                   ! h1g, cjg
      call mpp_clock_end (stratcloud_clock)
    endif  ! (do_lsc)
!---------------------------------------------------------------------
!    calculate the wet deposition associated with the large scale 
!    condensation. 
!---------------------------------------------------------------------
    qtnd_wet(is:ie,js:je,:) = qtnd(is:ie,js:je,:)
    if (do_strat) then
      qtnd_wet(is:ie,js:je,:) = qtnd_wet(is:ie,js:je,:) + q_tnd(is:ie,js:je,:,nql) + q_tnd(is:ie,js:je,:,nqi)
! Count precipitation formed over timestep plus cloud amount at end of timestep
      if (do_lin_cld_microphys) then
        cloud_wet(is:ie,js:je,:) = tracer(is:ie,js:je,:,nqr) + tracer(is:ie,js:je,:,nqs) + tracer(is:ie,js:je,:,nqg)
      else
        cloud_wet(is:ie,js:je,:) = rain3d(:,:,2:kx+1) - rain3d(:,:,1:kx) &
                  + snow3d(:,:,2:kx+1) - snow3d(:,:,1:kx)
        cloud_wet(is:ie,js:je,:) = cloud_wet(is:ie,js:je,:) * dt / pmass(is:ie,js:je,:) ! convert from kg/m2/s to kg/kg
      endif
      cloud_wet(is:ie,js:je,:) = cloud_wet(is:ie,js:je,:) + tracer(is:ie,js:je,:,nql) + tracer(is:ie,js:je,:,nqi)
      cloud_frac(is:ie,js:je,:) = max( min( tracer(is:ie,js:je,:,nqa), 1. ), 0. )
      used = send_data( id_f_snow_berg, f_snow_berg(:,:,:), Time,  &
                           is,js, 1)
      used = send_data( id_f_snow_berg_cond, f_snow_berg(:,:,:), Time,  &
                           is,js, 1,  mask = f_snow_berg /= 0.0 )
      used = send_data( id_f_snow_berg_wtd,   &
               f_snow_berg(:,:,:)*(rain3d(:,:,2:)+snow3d(:,:,2:)), Time,  &
                                   is,js, 1,  mask = f_snow_berg /= 0.0  )
    else
!     cloud_wet = qtnd_wet * dt
      cloud_wet(is:ie,js:je,:) = 0.5e-3
      cloud_frac(is:ie,js:je,:) = 1.
    end if
     ls_wetdep = 0.
    do n=1,size(rdt,4)
      if ( n /= nsphum ) then
        if ( .not. do_strat .or. (n /= nql .and. n /= nqi .and.   &
                         n /= nqa .and. n /= nqn  .and. n /= nqni ) ) then
          wetdeptnd(is:ie,js:je,:) = 0.0
          call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain, snow,   &
                               qtnd_wet(is:ie,js:je,:), cloud_wet(is:ie,js:je,:), &
                               cloud_frac(is:ie,js:je,:), f_snow_berg, &
                               rain3d, snow3d, &
                               tracer(is:ie,js:je,:,n), wetdeptnd(is:ie,js:je,:), Time, 'lscale',     &
                  is, js, dt, sum_wdep_out=ls_wetdep(:,:,n) )
          rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd(is:ie,js:je,:)
          wet_data(:,:,:,n) = wet_data(:,:,:,n) + wetdeptnd(is:ie,js:je,:)

          used = send_data( id_wet_deposition(n), wet_data(:,:,:,n), &
                            Time,is_in=is,js_in=js )
        end if
      end if
    end do

!---------------------------------------------------------------------
!    output diagnostics associated with the large-scale condensation
!    scheme.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    temperature change due to large-scale condensation:
!---------------------------------------------------------------------
    used = send_data (id_tdt_ls, ttnd(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
!---------------------------------------------------------------------
!    dry static energy tendency due to large-scale condensation:
!---------------------------------------------------------------------
    if (id_t_ls_col > 0) &
      call column_diag(id_t_ls_col, is, js, Time, ttnd(is:ie,js:je,:), CP_AIR) 
!---------------------------------------------------------------------
!    water vapor path tendency due to large-scale condensation:
!---------------------------------------------------------------------
    if (id_q_ls_col > 0) &
      call column_diag(id_q_ls_col, is, js, Time, qtnd(is:ie,js:je,:), 1.0) 

!---------------------------------------------------------------------
!    specific humidity change due to large-scale condensation:
!---------------------------------------------------------------------
    used = send_data (id_qdt_ls, qtnd(is:ie,js:je,:), Time, is, js, 1, rmask=mask)

    used = send_data (id_lsc_precip, rain + snow, Time, is, js)
        
    if (id_lsc_freq > 0) then
      ltemp = rain > 0. .or. snow > 0.0
      temp_2d = 0.
      where (ltemp) 
        temp_2d = 1.
      end where
      used = send_data (id_lsc_freq, temp_2d, Time, is, js)
    endif

!---------------------------------------------------------------------
!    total precipitation rate due to large-scale condensation:
!---------------------------------------------------------------------
    used = send_data (id_prec_ls, rain+snow, Time, is, js)

!---------------------------------------------------------------------
!    snowfall rate due to large-scale condensation:
!---------------------------------------------------------------------
    used = send_data (id_snow_ls, snow, Time, is, js)

!---------------------------------------------------------------------
!    define diagnostics specific to the strat_cloud formulation:
!---------------------------------------------------------------------
    if (do_strat) then

!---------------------------------------------------------------------
!    total cumulus mass flux due to strat_cloud parameterization:
!---------------------------------------------------------------------
      used = send_data (id_mc_full, mc_full(is:ie,js:je,:), Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    total cumulus mass flux on half levels:
!---------------------------------------------------------------------
      used = send_data (id_mc_half, mc_half(is:ie,js:je,:), Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud liquid, ice and area tendencies due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
      used = send_data (id_qldt_ls, q_tnd(is:ie,js:je,:,nql), Time, is, js, 1, rmask=mask)
      if (do_liq_num) used = send_data (id_qndt_ls, q_tnd(is:ie,js:je,:,nqn ), Time, is, js, 1, rmask=mask)
      used = send_data (id_qidt_ls, q_tnd(is:ie,js:je,:,nqi), Time, is, js, 1, rmask=mask)
      used = send_data (id_qadt_ls, q_tnd(is:ie,js:je,:,nqa), Time, is, js, 1, rmask=mask)
      if (do_ice_num) then
        used = send_data (id_qnidt_ls, q_tnd(is:ie,js:je,:,nqni),   &
                                              Time, is, js, 1, rmask=mask)
      endif

!---------------------------------------------------------------------
!    cloud liquid and ice water path tendencies due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
      if (id_ql_ls_col > 0) &
        call column_diag(id_ql_ls_col, is, js, Time, q_tnd(is:ie,js:je,:,nql), 1.0) 
      if (id_qi_ls_col > 0) &
        call column_diag(id_qi_ls_col, is, js, Time, q_tnd(is:ie,js:je,:,nqi), 1.0) 
      if (do_liq_num .and. id_qn_ls_col > 0) &
        call column_diag(id_qn_ls_col, is, js, Time, q_tnd(is:ie,js:je,:,nqn), 1.0) 
      if (do_ice_num .and. id_qni_ls_col > 0) &
        call column_diag(id_qni_ls_col, is, js, Time,   &
                                        q_tnd(is:ie,js:je,:,nqni), 1.0)
      
!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    strat_cloud  parameterization:
!---------------------------------------------------------------------
      if (id_enth_ls_col > 0) then
        temp_2d = -HLV*rain -HLS*snow
        call column_diag(id_enth_ls_col, is, js, Time, ttnd(is:ie,js:je,:), CP_AIR, &
                q_tnd(is:ie,js:je,:,nql), -HLV, q_tnd(is:ie,js:je,:,nqi), -HLS, temp_2d) 
      endif
 
      if (id_wat_ls_col > 0) then
        temp_2d = rain+snow
        call column_diag(id_wat_ls_col, is, js, Time, qtnd(is:ie,js:je,:), 1.0, &
                q_tnd(is:ie,js:je,:,nql), 1.0, q_tnd(is:ie,js:je,:,nqi), 1.0, temp_2d) 
      endif

!---------------------------------------------------------------------
!    stratiform cloud volume tendency due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
      if (id_qa_ls_col > 0) &
        call column_diag(id_qa_ls_col, is, js, Time, q_tnd(is:ie,js:je,:,nqa), 1.0)

!---------------------------------------------------------------------
!---- diagnostics for large scale precip -----------
!---------------------------------------------------------------------
      used = send_data(id_lscale_rain3d, rain3d, Time, is, js, 1)

!---------------------------------------------------------------------
!---- diagnostics for large scale snow -------------
!---------------------------------------------------------------------
      used = send_data(id_lscale_snow3d, snow3d, Time, is, js, 1)

!---------------------------------------------------------------------
!---- diagnostics for large scale precip -------------
!---------------------------------------------------------------------
      used = send_data(id_lscale_precip3d, snow3d(:,:,2:)+rain3d(:,:,2:), &
                       Time, is, js, 1, mask = f_snow_berg /= 0.0 )

    endif ! (do_strat)

!---------------------------------------------------------------------
!    end the timing of the large-scale condensation code section.
!---------------------------------------------------------------------
    call mpp_clock_end (largescale_clock)
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  GENERAL MOISTURE DIAGNOSTICS 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!--------------------------------------------------------------------
!    output diagnostics obtained from the combination of convective and
!    large-scale parameterizations.  
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    output diagnostics obtained from the combination of convective and
!    large-scale parameterizations.  
!--------------------------------------------------------------------

     total_wetdep(:,:,:) = 0.
     total_wetdep_donner(:,:,:) = 0.
     total_wetdep_uw    (:,:,:) = 0.
     m=1
     mm=1
     do n=1, size(rdt,4)
       if (tracers_in_donner(n)) then
         total_wetdep(:,:,n) = total_wetdep(:,:,n) +  &
                                                donner_wetdep(:,:,m)
         total_wetdep_donner(:,:,n) = donner_wetdep(:,:,m)
         m=m+1
       endif
       if (tracers_in_uw(n)) then
         total_wetdep(:,:,n) = total_wetdep(:,:,n) +  &
                                                    uw_wetdep(:,:,mm)
         total_wetdep_uw    (:,:,n) = uw_wetdep(:,:,mm)
         mm=mm+1
       endif
     end do
     if (do_strat) then
       total_wetdep = total_wetdep + ls_wetdep
     endif
     do n=1, size(rdt,4)
       if (id_wetdep(n) > 0) then
         used = send_data (id_wetdep(n), total_wetdep(:,:,n),  &
                                                          Time, is, js)
       endif
     end do
     if (id_wetdep_om > 0) then
       used = send_data (id_wetdep_om,  &
               total_wetdep       (:,:,nomphilic) + &
               total_wetdep       (:,:,nomphobic) , &
                                                Time, is,js)
     endif
     if (id_wetdep_SOA > 0) then
       used = send_data (id_wetdep_SOA,  &
               total_wetdep(:,:,nSOA) , Time, is,js)
     endif
     if (id_wetdep_bc > 0) then
       used = send_data (id_wetdep_bc,  &
               total_wetdep       (:,:,nbcphilic) + &
               total_wetdep       (:,:,nbcphobic) , &
                                                Time, is,js)
     endif
     if (id_wetdep_so4 > 0) then
       used = send_data (id_wetdep_so4,  &
          (96.0/WTMAIR)*(total_wetdep_donner(:,:,nso4     ) + &
               total_wetdep_uw    (:,:,nso4     )) + &
               0.096*ls_wetdep          (:,:,nso4     ) , &
                                                Time, is,js)
     endif
     if (id_wetdep_so2 > 0) then
       used = send_data (id_wetdep_so2,  &
          (64.0/WTMAIR)*(total_wetdep_donner(:,:,nso2     ) + &
               total_wetdep_uw    (:,:,nso2     )) + &
               0.064*ls_wetdep          (:,:,nso2     ) , &
                                                Time, is,js)
     endif
     if (id_wetdep_DMS > 0) then
       used = send_data (id_wetdep_DMS,  &
          (62.0/WTMAIR)*(total_wetdep_donner(:,:,nDMS     ) + &
               total_wetdep_uw    (:,:,nDMS     )) + &
               0.062*ls_wetdep          (:,:,nDMS     ) , &
                                                Time, is,js)
     endif
     if (id_wetdep_NH4NO3 > 0) then
       used = send_data (id_wetdep_NH4NO3,  &
           (18.0/WTMAIR)*(total_wetdep_donner(:,:,nnH4NO3  ) + &
               total_wetdep_donner(:,:,nNH4     ) + &
               total_wetdep_uw    (:,:,nNH4NO3  ) + &
               total_wetdep_uw    (:,:,nNH4     )) + &
            0.018*(ls_wetdep          (:,:,nNH4NO3  ) + &
                   ls_wetdep          (:,:,nNH4     )) , &
                                                Time, is,js)
     endif
     if (id_wetdep_salt   > 0) then
       used = send_data (id_wetdep_salt  ,  &
               ( total_wetdep(:,:,nsalt1) + &
                 total_wetdep(:,:,nsalt2) + &
                 total_wetdep(:,:,nsalt3) + &
                 total_wetdep(:,:,nsalt4) + &
                 total_wetdep(:,:,nsalt5)),  Time, is,js)
     endif
     if (id_wetdep_dust   > 0) then
       used = send_data (id_wetdep_dust  ,  &
               ( total_wetdep(:,:,ndust1) + &
                 total_wetdep(:,:,ndust2) + &
                 total_wetdep(:,:,ndust3) + &
                 total_wetdep(:,:,ndust4) + &
                 total_wetdep(:,:,ndust5)),  Time, is,js)
     endif

!---------------------------------------------------------------------
!    total precipitation (all sources):
!---------------------------------------------------------------------
    precip = fprec + lprec
    if (id_precip > 0) then
      used = send_data (id_precip, precip, Time, is, js)
    endif

!---------------------------------------------------------------------
!    snowfall rate due to all sources:
!---------------------------------------------------------------------
    used = send_data (id_snow_tot, fprec, Time, is, js)

!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    moist processes:
!---------------------------------------------------------------------

    if (id_enth_moist_col > 0 .or. id_max_enthalpy_imbal > 0) then
      temp_3d1 = tdt - tdt_init(is:ie,js:je,:)
      temp_3d2 = rdt(:,:,:,nql) - rdt_init(is:ie,js:je,:,nql)
      temp_3d3 = rdt(:,:,:,nqi) - rdt_init(is:ie,js:je,:,nqi)
      temp_2d(:,:) = -HLV*precip -HLF*fprec
      call column_diag(id_enth_moist_col, is, js, Time, temp_3d1, CP_AIR, temp_3d2, -HLV, temp_3d3, -HLS, temp_2d)
      if (id_max_enthalpy_imbal > 0) then
        max_enthalpy_imbal = max( abs(temp_2d), max_enthalpy_imbal )
        used = send_data(id_max_enthalpy_imbal, max_enthalpy_imbal, Time, is, js)
      endif
    endif
  
    if (id_wat_moist_col > 0 .or. id_max_water_imbal > 0) then
      temp_3d1 = qdt - qdt_init(is:ie,js:je,:)
      temp_3d2 = rdt(:,:,:,nql) - rdt_init(is:ie,js:je,:,nql)
      temp_3d3 = rdt(:,:,:,nqi) - rdt_init(is:ie,js:je,:,nqi)
      temp_2d(:,:) = precip
      call column_diag(id_wat_moist_col, is, js, Time, temp_3d1, 1.0, temp_3d2, 1.0, temp_3d3, 1.0, temp_2d)
      if (id_max_water_imbal > 0) then
        max_water_imbal = max( abs(temp_2d), max_water_imbal )
        used = send_data(id_max_water_imbal, max_water_imbal, Time, is, js)
      endif
    endif

!---------------------------------------------------------------------
!    water vapor, liquid water and ice water column paths:
!---------------------------------------------------------------------
    if (id_WVP > 0) &
      call column_diag(id_WVP, is, js, Time, qin(is:ie,js:je,:), 1.0) 
    if (id_LWP > 0 .and. do_strat) &
      call column_diag(id_LWP, is, js, Time, tracer(is:ie,js:je,:,nql), 1.0) 
    if (id_IWP > 0 .and. do_strat) &
      call column_diag(id_IWP, is, js, Time, tracer(is:ie,js:je,:,nqi), 1.0) 

    if (id_lsc_cloud_area > 0 .and. do_strat ) then
      used = send_data (id_lsc_cloud_area, 100.*lsc_cloud_area, Time, is,  &
                        js, 1, rmask=mask)
    end if

!----------------------------------------------------------------------
!    define total convective cloud amount (grid-box mean).
!----------------------------------------------------------------------
    total_conv_cloud = 0.
    tot_conv_liq = 0.
    tot_conv_ice = 0.
    conv_cld_frac = 0.
    total_cloud_area = 0.
    if (do_strat) then
      total_cloud_area = total_cloud_area + lsc_cloud_area
    endif
    if (do_donner_deep) then
      total_conv_cloud = total_conv_cloud +   &
          cell_cld_frac*cell_ice_amt + meso_cld_frac*meso_ice_amt +  &
          cell_cld_frac*cell_liq_amt + meso_cld_frac*meso_liq_amt 
      conv_cld_frac = conv_cld_frac + cell_cld_frac + meso_cld_frac
      total_cloud_area = total_cloud_area + cell_cld_frac +  &
                                                         meso_cld_frac
      tot_conv_liq = tot_conv_liq + cell_cld_frac*cell_liq_amt + &
                                           meso_cld_frac*meso_liq_amt 
      tot_conv_ice = tot_conv_ice + cell_cld_frac*cell_ice_amt + &
                                           meso_cld_frac*meso_ice_amt 
    endif
    if (do_uw_conv) then
      total_conv_cloud = total_conv_cloud +   &
                          shallow_cloud_area*shallow_ice  +   &
                          shallow_cloud_area*shallow_liquid
      conv_cld_frac = conv_cld_frac + shallow_cloud_area
      total_cloud_area = total_cloud_area + shallow_cloud_area
      tot_conv_liq = tot_conv_liq + shallow_cloud_area*shallow_liquid
      tot_conv_ice = tot_conv_ice + shallow_cloud_area*shallow_ice     
    endif

    if (id_lsc_liq_amt > 0 .and. do_strat ) then
      used = send_data (id_lsc_liq_amt,  &
                        lsc_liquid/(1.0 + total_conv_cloud), &
                                          Time, is, js, 1, rmask=mask)
    end if

    if ( id_lsc_ice_amt  > 0 .and. do_strat ) then
       used = send_data (id_lsc_ice_amt,   &
                         lsc_ice/(1.0 + total_conv_cloud), &
                                         Time, is, js, 1, rmask=mask)
    end if


    used = send_data (id_tot_cloud_area, 100.*total_cloud_area,  &
                                          Time, is, js, 1, rmask=mask)

    used = send_data (id_conv_cloud_area, 100.*conv_cld_frac, &
                                           Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    define the total cloud area. 
!---------------------------------------------------------------------
    if (id_tot_cld_amt > 0 ) then
      tca2 = 1.0
      do k=1,kx
        tca2(:,:) = tca2(:,:)*(1.0 - total_cloud_area(:,:,k))
      end do
      tca2 = 100.*(1. - tca2)
      used = send_data (id_tot_cld_amt, tca2, Time, is, js)
    endif
 
!---------------------------------------------------------------------
!    define the total and convective liquid and liquid water path. 
!---------------------------------------------------------------------
    if (id_tot_liq_amt > 0 ) &
    used = send_data (id_tot_liq_amt, &
                 (lsc_liquid + tot_conv_liq)/(1.0 + total_conv_cloud), &
                                           Time, is, js, 1, rmask=mask)

    if (id_conv_liq_amt > 0 ) &
    used = send_data (id_conv_liq_amt, &
                   tot_conv_liq /(1.0 + total_conv_cloud), &
                                           Time, is, js, 1, rmask=mask)
 
    if (id_LWP_all_clouds > 0 ) &
    call column_diag (id_LWP_all_clouds, is, js, Time, &
                      lsc_liquid+tot_conv_liq, 1.0)

!---------------------------------------------------------------------
!    define the total and convective ice and ice water path. 
!---------------------------------------------------------------------

    if (id_tot_ice_amt > 0 ) &
     used = send_data (id_tot_ice_amt, &
                 (lsc_ice + tot_conv_ice)/(1.0 + total_conv_cloud), &
                                            Time, is, js, 1, rmask=mask)

    if (id_conv_ice_amt > 0 ) &
     used = send_data (id_conv_ice_amt, &
                  tot_conv_ice/(1.0 + total_conv_cloud), &
                                            Time, is, js, 1, rmask=mask)

    if (id_IWP_all_clouds > 0 ) &
     call column_diag (id_IWP_all_clouds, is, js, Time, &
                       lsc_ice+tot_conv_ice, 1.0)
 
!---------------------------------------------------------------------
!    define the total vapor, total water substance and condensate water
!    path. 
!---------------------------------------------------------------------
    used = send_data (id_tot_vapor, qin, Time, is, js, 1, rmask=mask)
    used = send_data (id_tot_h2o  , &
                     (qin(is:ie,js:je,:) + lsc_ice + tot_conv_ice + lsc_liquid +  &
                               tot_conv_liq)/(1.0 + total_conv_cloud), &
                                            Time, is, js, 1, rmask=mask)

    if (id_WP_all_clouds > 0 ) &
    call column_diag(id_WP_all_clouds, is, js, Time, &
              lsc_ice + tot_conv_ice +  lsc_liquid + tot_conv_liq, 1.0)

!---------------------------------------------------------------------
!    column integrated cloud mass:
!---------------------------------------------------------------------
    if (id_AWP > 0 .and. do_strat) &
      call column_diag(id_AWP, is, js, Time, tracer(is:ie,js:je,:,nqa), 1.0) 

!---------------------------------------------------------------------
!    relative humidity:         
!---------------------------------------------------------------------
    if (id_rh > 0) then
      if (.not. (do_rh_clouds .or. do_diag_clouds)) then 
        call rh_calc (pfull, tin(is:ie,js:je,:), qin(is:ie,js:je,:), RH(is:ie,js:je,:), do_simple, mask)
      endif
      used = send_data (id_rh, rh(is:ie,js:je,:)*100., Time, is, js, 1, rmask=mask)
    endif

!---------------------------------------------------------------------
!    relative humidity (CMIP formulation):         
!---------------------------------------------------------------------
    if (id_rh_cmip > 0) then
      if (.not. (do_rh_clouds .or. do_diag_clouds)) then
        call rh_calc (pfull, tin(is:ie,js:je,:), qin(is:ie,js:je,:), &
                       RH(is:ie,js:je,:), do_simple, do_cmip=.true., &
                                                            mask=mask)
      endif
      used = send_data (id_rh_cmip, rh(is:ie,js:je,:)*100.,  &
                                         Time, is, js, 1, rmask=mask)
    endif

!---------------------------------------------------------------------
!    saturation specific humidity:         
!---------------------------------------------------------------------
    if (id_qs > 0) then
      call compute_qs (tin(is:ie,js:je,:), pfull, qsat(is:ie,js:je,:), q = qin(is:ie,js:je,:))
      used = send_data (id_qs, qsat(is:ie,js:je,:), Time, is, js, 1, rmask=mask)
    endif

!------- diagnostics for CAPE and CIN, 

!!-- compute and write out CAPE and CIN--
    if ( id_cape > 0 .or. id_cin > 0) then
!! calculate r
      rin(is:ie,js:je,:) = qin(is:ie,js:je,:)/(1.0 - qin(is:ie,js:je,:)) ! XXX rin is not mixing ratio when do_simple=.true.
      avgbl = .false.
      do j=js,je
       do i=is,ie
         call capecalcnew( kx, pfull(i,j,:), phalf(i,j,:), CP_AIR, RDGAS, RVGAS, &
                   HLV, KAPPA, tin(i,j,:), rin(i,j,:), avgbl, cape(i,j), cin(i,j))
       end do
      end do
      if (id_cape > 0) used = send_data ( id_cape, cape, Time, is, js )
      if ( id_cin > 0 ) used = send_data ( id_cin, cin, Time, is, js )
    end if

!---------------------------------------------------------------------
!    output the global integral of precipitation in units of mm/day.
!---------------------------------------------------------------------
    prec_intgl(is:ie,js:je) = precip(:,:)*SECONDS_PER_DAY

!----------------------------------------------------------------------
!    define the precip fluxes needed for input to the COSP simulator 
!    package.
!---------------------------------------------------------------------
      if (do_cosp) then

!---------------------------------------------------------------------
!    define precip fluxes from convective schemes at each layer 
!    interface.  (index 1 is model lid)
!---------------------------------------------------------------------
        liq_mesoh(is:ie,js:je,1) = 0.
        frz_mesoh(is:ie,js:je,1) = 0.
        liq_cellh(is:ie,js:je,1) = 0.
        frz_cellh(is:ie,js:je,1) = 0.
        ice_precflxh(is:ie,js:je,1) = 0.
        liq_precflxh(is:ie,js:je,1) = 0.
        mca_liqh(is:ie,js:je,1) = 0.
        mca_frzh(is:ie,js:je,1) = 0.
        do k=2, size(t,3)+1
          liq_mesoh(is:ie,js:je,k) =  liq_mesoh (is:ie,js:je,k-1) + &
                                      liq_meso (is:ie,js:je,k-1)
          frz_mesoh(is:ie,js:je,k) =  frz_mesoh (is:ie,js:je,k-1) + &
                                      frz_meso (is:ie,js:je,k-1)
          liq_cellh(is:ie,js:je,k) =  liq_cellh (is:ie,js:je,k-1) + &
                                      liq_cell (is:ie,js:je,k-1)
          frz_cellh(is:ie,js:je,k) =  frz_cellh (is:ie,js:je,k-1) + &
                                      frz_cell (is:ie,js:je,k-1)
          ice_precflxh(is:ie,js:je,k) =                  &
                               ice_precflxh(is:ie,js:je,k-1) +  &
                                           ice_precflx(is:ie,js:je,k-1)
          liq_precflxh(is:ie,js:je,k) =        &
                               liq_precflxh(is:ie,js:je,k-1) +   &
                                           liq_precflx(is:ie,js:je,k-1)
          if (include_donmca_in_cosp) then
            mca_liqh (is:ie,js:je,k) = mca_liqh (is:ie,js:je,k-1) + &
                                        mca_liq(is:ie,js:je,k-1)
            mca_frzh (is:ie,js:je,k) = mca_frzh (is:ie,js:je,k-1) + &
                                        mca_frz(is:ie,js:je,k-1)
          endif
        end do

!--------------------------------------------------------------------
!    adjust precip fluxes to account for any negative values produced.
!    precip contribution is determined as the negative of the moisture
!    tendency, so at top of clouds a positive moisture tendency some-
!    times results in a negative precipitation contribution. 
!--------------------------------------------------------------------
        sumneg(is:ie,js:je) = 0.
        do k=2, size(t,3)+1
        do j=js,je        
        do i=is,ie          
          if (liq_mesoh(i,j,k) > 0.0) then
            if (liq_mesoh(i,j,k) > ABS(sumneg(i,j))) then
              liq_mesoh(i,j,k) = liq_mesoh(i,j,k) + sumneg(i,j)
              sumneg(i,j) = 0.
            else
              sumneg(i,j) = sumneg(i,j) + liq_mesoh(i,j,k)
              liq_mesoh(i,j,k) = 0.                        
            endif
          else
            sumneg(i,j) = sumneg(i,j) + liq_mesoh(i,j,k)
            liq_mesoh(i,j,k) = 0.
          endif
        end do
        end do
        end do
        sumneg(is:ie,js:je) = 0.
        do k=2, size(t,3)+1
        do j=js,je            
        do i=is,ie          
          if (frz_mesoh(i,j,k) > 0.0) then
            if (frz_mesoh(i,j,k) > ABS(sumneg(i,j))) then
              frz_mesoh(i,j,k) = frz_mesoh(i,j,k) + sumneg(i,j)
              sumneg(i,j) = 0.
            else
              sumneg(i,j) = sumneg(i,j) + frz_mesoh(i,j,k)
              frz_mesoh(i,j,k) = 0.                        
            endif
          else
            sumneg(i,j) = sumneg(i,j) + frz_mesoh(i,j,k)
            frz_mesoh(i,j,k) = 0.
          endif
        end do
        end do
        end do
        sumneg(is:ie,js:je) = 0.
        do k=2, size(t,3)+1
        do j=js,je          
        do i=is,ie            
          if (liq_cellh(i,j,k) > 0.0) then
            if (liq_cellh(i,j,k) > ABS(sumneg(i,j))) then
              liq_cellh(i,j,k) = liq_cellh(i,j,k) + sumneg(i,j)
              sumneg(i,j) = 0.
            else
              sumneg(i,j) = sumneg(i,j) + liq_cellh(i,j,k)
              liq_cellh(i,j,k) = 0.                        
            endif
          else
            sumneg(i,j) = sumneg(i,j) + liq_cellh(i,j,k)
            liq_cellh(i,j,k) = 0.
          endif
        end do
        end do
        end do
        sumneg(is:ie,js:je) = 0.
        do k=2, size(t,3)+1
        do j=js,je            
        do i=is,ie             
          if (frz_cellh(i,j,k) > 0.0) then
            if (frz_cellh(i,j,k) > ABS(sumneg(i,j))) then
              frz_cellh(i,j,k) = frz_cellh(i,j,k) + sumneg(i,j)
              sumneg(i,j) = 0.
            else
              sumneg(i,j) = sumneg(i,j) + frz_cellh(i,j,k)
              frz_cellh(i,j,k) = 0.                        
            endif
          else
            sumneg(i,j) = sumneg(i,j) + frz_cellh(i,j,k)
            frz_cellh(i,j,k) = 0.
          endif
        end do
        end do
        end do
        sumneg(is:ie,js:je) = 0.
        do k=2, size(t,3)+1
        do j=js,je           
        do i=is,ie           
          if (ice_precflxh(i,j,k) > 0.0) then
            if (ice_precflxh(i,j,k) > ABS(sumneg(i,j))) then
              ice_precflxh(i,j,k) = ice_precflxh(i,j,k) + sumneg(i,j)
              sumneg(i,j) = 0.
            else
              sumneg(i,j) = sumneg(i,j) + ice_precflxh(i,j,k)
              ice_precflxh(i,j,k) = 0.                        
            endif
          else
            sumneg(i,j) = sumneg(i,j) + ice_precflxh(i,j,k)
            ice_precflxh(i,j,k) = 0.
          endif
        end do
        end do
        end do
        sumneg(is:ie,js:je) = 0.
        do k=2, size(t,3)+1
        do j=js,je          
        do i=is,ie              
          if (liq_precflxh(i,j,k) > 0.0) then
            if (liq_precflxh(i,j,k) > ABS(sumneg(i,j))) then
              liq_precflxh(i,j,k) = liq_precflxh(i,j,k) + sumneg(i,j)
              sumneg(i,j) = 0.
            else
              sumneg(i,j) = sumneg(i,j) + liq_precflxh(i,j,k)
              liq_precflxh(i,j,k) = 0.                        
            endif
          else
            sumneg(i,j) = sumneg(i,j) + liq_precflxh(i,j,k)
            liq_precflxh(i,j,k) = 0.
          endif
        end do
        end do
        end do
        if (include_donmca_in_cosp) then
        sumneg(is:ie,js:je) = 0.
          do k=2, size(t,3)+1
          do j=js,je          
          do i=is,ie             
            if (mca_liqh(i,j,k) > 0.0) then
              if (mca_liqh(i,j,k) > ABS(sumneg(i,j))) then
                mca_liqh(i,j,k) = mca_liqh(i,j,k) + sumneg(i,j)
                sumneg(i,j) = 0.
              else
                sumneg(i,j) = sumneg(i,j) + mca_liqh(i,j,k)
                mca_liqh(i,j,k) = 0.                        
              endif
            else
              sumneg(i,j) = sumneg(i,j) + mca_liqh(i,j,k)
              mca_liqh(i,j,k) = 0.
            endif
          end do
          end do
          end do
        sumneg(is:ie,js:je) = 0.
          do k=2, size(t,3)+1
          do j=js,je          
          do i=is,ie            
            if (mca_frzh(i,j,k) > 0.0) then
              if (mca_frzh(i,j,k) > ABS(sumneg(i,j))) then
                mca_frzh(i,j,k) = mca_frzh(i,j,k) + sumneg(i,j)
                sumneg(i,j) = 0.
              else
                sumneg(i,j) = sumneg(i,j) + mca_frzh(i,j,k)
                mca_frzh(i,j,k) = 0.                        
              endif
            else
              sumneg(i,j) = sumneg(i,j) + mca_frzh(i,j,k)
              mca_frzh(i,j,k) = 0.
            endif
          end do
          end do
          end do
        endif

!----------------------------------------------------------------------
!     define the grid-box precip flux as the average of the interface 
!     fluxes.
!----------------------------------------------------------------------
        do k=1, size(t,3)
          do j=1, size(t,2)
            do i=1, size(t,1)
              if (donner_meso_is_largescale) then
                fl_lsrain(i,j,k) =  0.5*   &
                   ((strat_precip_in_cosp*  &
                                 (rain3d(i,j,k) + rain3d(i,j,k+1)) + &
                    donner_precip_in_cosp*  &
                                 (liq_mesoh(i+is-1,j+js-1,k) +    &
                                          liq_mesoh(i+is-1,j+js-1,k+1))))
                fl_lssnow(i,j,k) = 0.5*   &
                   ((strat_precip_in_cosp*   &
                          (snowclr3d(i,j,k) + snowclr3d(i,j,k+1)) + &
                    donner_precip_in_cosp*  &
                           (frz_mesoh(i+is-1,j+js-1,k) +  &
                                    frz_mesoh(i+is-1,j+js-1,k+1))))
                fl_ccrain(i,j,k) =  0.5*  &
                    ((donner_precip_in_cosp* &
                             (liq_cellh(i+is-1,j+js-1,k) +  &
                                     liq_cellh(i+is-1,j+js-1,k+1)) +  &
                      uw_precip_in_cosp*   &
                        (liq_precflxh(i+is-1,j+js-1,k) +  &
                                       liq_precflxh(i+is-1,j+js-1,k+1))))
                fl_ccsnow(i,j,k) =  0.5*  &
                    ((donner_precip_in_cosp*  &
                          (frz_cellh(i+is-1,j+js-1,k) +  &
                                     frz_cellh(i+is-1,j+js-1,k+1))  +  &
                       uw_precip_in_cosp*  &
                           (ice_precflxh(i+is-1,j+js-1,k) +   &
                                     ice_precflxh(i+is-1,j+js-1,k+1))))
              else
                fl_lsrain(i,j,k) =  0.5*   &
                   strat_precip_in_cosp*   &
                             (rain3d(i,j,k) + rain3d(i,j,k+1))
                fl_lssnow(i,j,k) =  0.5*  &
                   strat_precip_in_cosp*    &
                            (snowclr3d(i,j,k) + snowclr3d(i,j,k+1))
                fl_ccrain(i,j,k) =    0.5* &
                   ((donner_precip_in_cosp*  &
                        (liq_mesoh(i+is-1,j+js-1,k) +    &
                                      liq_mesoh(i+is-1,j+js-1,k+1) +  &
                        liq_cellh(i+is-1,j+js-1,k) +    &
                                      liq_cellh(i+is-1,j+js-1,k+1)) +  &
                    uw_precip_in_cosp*   &
                        (liq_precflxh(i+is-1,j+js-1,k) +    &
                                      liq_precflxh(i+is-1,j+js-1,k+1))))
                fl_ccsnow(i,j,k) =  0.5*  &
                   ((donner_precip_in_cosp*  &
                              (frz_mesoh(i+is-1,j+js-1,k) +    &
                                      frz_mesoh(i+is-1,j+js-1,k+1) +  &
                        frz_cellh(i+is-1,j+js-1,k) +  &
                                      frz_cellh(i+is-1,j+js-1,k+1)) +  &
                     uw_precip_in_cosp*  &
                        (ice_precflxh(i+is-1,j+js-1,k) +    &
                                      ice_precflxh(i+is-1,j+js-1,k+1))))
              endif
              if (include_donmca_in_cosp .and. &
                   donner_precip_in_cosp .eq. 1.0) then
                fl_donmca_snow(i,j,k) =                         0.5*  &
                                   (mca_frzh(i+is-1,j+js-1,k) +   &
                                           mca_frzh(i+is-1,j+js-1,k+1))
                fl_donmca_rain(i,j,k) =                         0.5*  &
                                   (mca_liqh(i+is-1,j+js-1,k) +   &
                                           mca_liqh(i+is-1,j+js-1,k+1))
              endif
            end do
          end do
        end do

      endif ! (do_cosp)

!-----------------------------------------------------------------------
end subroutine moist_processes


!#####################################################################
 
subroutine moist_processes_time_vary (dt)

real, intent(in) :: dt


      if (do_donner_deep) then
        call donner_deep_time_vary (dt)
      endif
      if (do_strat .and. .not. do_legacy_strat_cloud) then
        call strat_cloud_time_vary (dt, limit_conv_cloud_frac)
      endif

end subroutine moist_processes_time_vary


!#####################################################################

subroutine moist_processes_endts (is, js)
 
integer, intent(in) :: is,js

      if (do_donner_deep) then
        call donner_deep_endts
      endif 


      call sum_diag_integral_field ('prec', prec_intgl)
      prec_intgl = 0.0


end subroutine moist_processes_endts



!###################################################################

!#######################################################################
!---> h1g
!subroutine moist_processes_init ( id, jd, kd, lonb, latb, pref, &
subroutine moist_processes_init ( id, jd, kd, lonb, latb, lon, lat, phalf, pref, &
                                  axes, Time, doing_donner, &
                                  doing_uw_conv, num_uw_tracers_out,&
                                  do_strat_out,     &
                                  do_clubb_in,      &    ! cjg
                                  do_cosp_in,  &
!                                 doing_uw_conv, &
!                                 do_cosp_in,  &
                                  donner_meso_is_largescale_in, &
                                  include_donmca_in_cosp_out)
!<--- h1g

!-----------------------------------------------------------------------
integer,              intent(in)  :: id, jd, kd, axes(4)
real, dimension(:,:), intent(in)  :: lonb, latb
real,dimension(:,:),  intent(in)  :: lon,  lat    ! h1g
real,dimension(:,:,:),intent(in)  :: phalf        ! h1g
real, dimension(:),   intent(in)  :: pref
type(time_type),      intent(in)  :: Time
 logical,              intent(out) :: doing_donner, doing_uw_conv,   &
                                      do_strat_out
!logical,              intent(out) :: doing_donner, doing_uw_conv
 integer,              intent(out) :: num_uw_tracers_out
integer,              intent(in), optional :: do_clubb_in        ! cjg
logical,              intent(in), optional ::   &
                                     do_cosp_in, &
                                     donner_meso_is_largescale_in
logical,             intent(out), optional ::    &
                                     include_donmca_in_cosp_out
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!-----------------------------------------------------------------------

integer :: unit,io,ierr, n, logunit
character(len=80)  :: scheme
integer            :: secs, days
integer            :: k
!-----------------------------------------------------------------------

       if ( module_is_initialized ) return

!-->cjg
       if (present(do_clubb_in)) then
         do_clubb = do_clubb_in
       else
         do_clubb = 0
       endif
!<--cjg

       if (present(do_cosp_in)) then
         do_cosp = do_cosp_in
       else
         do_cosp = .false.
       endif
       if (present(donner_meso_is_largescale_in)) then
         donner_meso_is_largescale = donner_meso_is_largescale_in
       else
         donner_meso_is_largescale = .false.
       endif

       if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
         read (input_nml_file, nml=moist_processes_nml, iostat=io)
         ierr = check_nml_error(io,'moist_processes_nml')
#else

         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=moist_processes_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'moist_processes_nml')
         enddo
  10     call close_file (unit)
#endif

!--------- write version and namelist to standard log ------------

      call write_version_number(version, tagname)
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=moist_processes_nml )

       endif

      if (present(include_donmca_in_cosp_out)) then
        include_donmca_in_cosp_out = include_donmca_in_cosp
      endif

!------------------- dummy checks --------------------------------------
      if (do_ice_num .and. .not. do_liq_num) then
        call error_mesg ('moist_processes_mod',  &
                'do_ice_num can only be selected if do_liq_num is &
                                                      &selected', FATAL)
      endif

         if (include_donmca_in_cosp .and. (.not. do_donner_mca) ) &
           call error_mesg ('moist_processes_init', &
          'want to include donmca in COSP when donmca is inactive', &
                                                                  FATAL)

         if ( do_mca .and. do_ras ) call error_mesg   &
                   ('moist_processes_init',  &
                    'both do_mca and do_ras cannot be specified', FATAL)

         if ( do_mca .and. do_bm ) call error_mesg   &
                   ('moist_processes_init',  &
                    'both do_mca and do_bm cannot be specified', FATAL)
         if ( do_ras .and. do_bm ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bm and do_ras cannot be specified', FATAL)
         if ( do_bm .and. do_bmmass ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bm and do_bmmass cannot be specified', FATAL)
         if ( do_bm .and. do_bmomp ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bm and do_bmomp cannot be specified', FATAL)
         if ( do_bmomp .and. do_bmmass ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmomp and do_bmmass cannot be specified', FATAL)
         if ( do_bmmass .and. do_mca ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmmass and do_mca cannot be specified', FATAL)
         if ( do_bmmass .and. do_ras ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmmass and do_ras cannot be specified', FATAL)
         if ( do_bmomp .and. do_mca ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmomp and do_mca cannot be specified', FATAL)
         if ( do_bmomp .and. do_ras ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmomp and do_ras cannot be specified', FATAL)

         if ( do_lsc .and. do_strat ) call error_mesg   &
                 ('moist_processes_init',  &
                  'both do_lsc and do_strat cannot be specified', FATAL)
         if (.not. do_lsc .and. .not. do_strat) then
           call error_mesg ('moist_processes_mod', &
              'must activate either do_lsc or do_strat in order to &
                             &include large-scale condensation', FATAL)
         endif

         if ( (do_rh_clouds.or.do_diag_clouds) .and. do_strat .and. &
             mpp_pe() == mpp_root_pe() ) call error_mesg ('moist_processes_init', &
     'do_rh_clouds or do_diag_clouds + do_strat should not be specified', NOTE)

         if ( do_rh_clouds .and. do_diag_clouds .and. mpp_pe() == mpp_root_pe() ) &
            call error_mesg ('moist_processes_init',  &
       'do_rh_clouds and do_diag_clouds should not be specified', NOTE)

         if (do_mca .and. do_donner_deep) call error_mesg &
                 ('moist_processes_init',  &
            'both do_donner_deep and do_mca cannot be specified', FATAL)

         if (do_donner_deep .and. do_rh_clouds) then
           call error_mesg ('moist_processes_init',  &
            'Cannot currently activate donner_deep_mod with rh_clouds', FATAL)
         endif   
         
         if (force_donner_moist_conserv .and. &
               .not. do_donner_conservation_checks) then
           call error_mesg ('moist_processes', &
              'when force_donner_moist_conserv is .true., &
                &do_donner_conservation_checks must be .true.', FATAL)
         endif

         if (use_updated_profiles_for_uw .and.   &
             .not. (do_donner_before_uw) ) then
           call error_mesg ('moist_processes_init', &
            'use_updated_profiles_for_uw is only meaningful when &
                               &do_donner_before_uw is true', FATAL)
         endif

         if (only_one_conv_scheme_per_column .and.   &
             .not. (do_donner_before_uw) ) then
           call error_mesg ('moist_processes_init', &
            'only_one_conv_scheme_per_column is only meaningful when &
                               &do_donner_before_uw is true', FATAL)
         endif

         if (limit_conv_cloud_frac .and.   &
                 .not. do_donner_before_uw) then
           call error_mesg ('moist_processes', &
              'when limit_conv_cloud_frac is .true., &
                 &do_donner_before_uw must be .true.', FATAL)
         endif

         if (do_cosp .and. .not. do_donner_conservation_checks) then
           do_donner_conservation_checks = .true.
           call error_mesg ('moist_processes', &
              'setting do_donner_conservation_checks to true so that &
                 &needed fields for COSP are produced.', NOTE)
         endif

!RSH  endif

!---------------------------------------------------------------------
! --- Find the tracer indices 
!---------------------------------------------------------------------

      if (do_strat) then
        ! get tracer indices for stratiform cloud variables
        nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
        nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
        nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
        nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
        if (min(nql,nqi,nqa) <= 0) call error_mesg ('moist_processes', &
                                                    'stratiform cloud tracer(s) not found', FATAL)
        if (nql == nqi .or. nqa == nqi .or. nql == nqa) call error_mesg ('moist_processes',  &
                                 'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)
        if (mpp_pe() == mpp_root_pe()) &
            write (logunit,'(a,3i4)') 'Stratiform cloud tracer indices: nql,nqi,nqa =',nql,nqi,nqa
      endif

      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      if (nqn == NO_TRACER .and. do_liq_num ) &
        call error_mesg ('moist_processes', &
             'prognostic droplet number scheme requested but tracer not found', FATAL)
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
      if (nqni == NO_TRACER .and. do_ice_num ) &
         call error_mesg ('moist_processes', &
            'prognostic ice number scheme requested but &
                                              &tracer not found', FATAL) 

!------------ initialize various schemes ----------
      if (do_lsc) then
                     call lscale_cond_init ()
                     if (do_rh_clouds) call rh_clouds_init (id,jd,kd)
                     if (do_diag_clouds) call diag_cloud_init (id,jd,kd,ierr)
      endif

! ---> h1g, cjg
      if (do_strat) then
        if (do_clubb > 0) then
! --->h1g, if CLUBB is in moist-processes, CLUBB is initialized here.
          if( do_clubb == 2) then
             call clubb_init( id, jd, kd, lon, lat,    &
                         axes, Time,  phalf )
          endif
! <---h1g, if CLUBB is in moist-processes, CLUBB is initialized here.
          call MG_microp_3D_init(axes,Time,id,jd,kd)
        end if
        call strat_cloud_init (axes, Time, id, jd, kd,    &
                                 do_legacy_strat_cloud = do_legacy_strat_cloud)
       end if
! <--- h1g, cjg
      if (do_dryadj) call     dry_adj_init ()
      if (do_cmt)    call cu_mo_trans_init (axes,Time, doing_diffusive)
      if (do_bm)     call betts_miller_init () 
      if (do_bmmass) call bm_massflux_init()
      if (do_bmomp)  call bm_omp_init () 

      if (do_cmt) then
        if ( .not. do_ras .and. .not. do_donner_deep  .and. &
             .not. do_uw_conv) then
          call error_mesg ( 'moist_processes_mod', &
                'do_cmt specified but no cumulus schemes activated', &
                                                              FATAL)
        endif
        if (trim(cmt_mass_flux_source) == 'ras') then
          cmt_uses_ras = .true.
          cmt_uses_donner = .false.
          cmt_uses_uw = .false.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif

        else if (trim(cmt_mass_flux_source) == 'donner') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .true.
          cmt_uses_uw = .false.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif

        else if (trim(cmt_mass_flux_source) == 'uw') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .false.
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_uw then uw_conv_mod must&
                                               & be active', FATAL)
           endif

        else if (trim(cmt_mass_flux_source) == 'donner_and_ras') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_uw = .false.

        else if (trim(cmt_mass_flux_source) == 'donner_and_uw') then
          cmt_uses_uw = .true.
          if (.not. do_uw_conv) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_uw then uw_conv_mod must be active', FATAL)
           endif
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_ras = .false.

        else if (trim(cmt_mass_flux_source) == 'ras_and_uw') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_uw then uw_conv_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_donner = .false.

        else if (trim(cmt_mass_flux_source) == 'donner_and_ras_and_uw') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_uw then uw_conv_mod must&
                                               & be active', FATAL)
           endif
        else if (trim(cmt_mass_flux_source) == 'all') then
          if (do_ras) then
            cmt_uses_ras = .true.
          else
            cmt_uses_ras = .false.
          endif
          if (do_donner_deep)  then
            cmt_uses_donner = .true.
          else
            cmt_uses_donner = .false.
          endif
          if (do_uw_conv)  then
            cmt_uses_uw = .true.
          else
            cmt_uses_uw = .false.
          endif
        else
          call error_mesg ('moist_processes_mod', &
             'invalid specification of cmt_mass_flux_source', FATAL)
        endif

        if (cmt_uses_uw .and. .not. doing_diffusive) then
          call error_mesg ('moist_processes_mod', &
             'currently cannot do non-local cmt with uw as mass &
                                                &flux_source', FATAL)
        endif
          

      endif

  
!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')
   allocate (prec_intgl(id,jd))

!---------------------------------------------------------------------
!    define output variables indicating whether certain convection 
!    schemes have been activated.
!---------------------------------------------------------------------
     doing_donner = do_donner_deep
     doing_uw_conv = do_uw_conv

!----- initialize clocks -----
   convection_clock = mpp_clock_id( '   Physics_up: Moist Proc: Conv' , grain=CLOCK_MODULE )
   largescale_clock = mpp_clock_id( '   Physics_up: Moist Proc: LS'   , grain=CLOCK_MODULE )
   donner_clock     = mpp_clock_id( '   Moist Processes: Donner_deep' , grain=CLOCK_MODULE )
   mca_clock        = mpp_clock_id( '   Moist Processes: MCA'         , grain=CLOCK_MODULE )
   donner_mca_clock = mpp_clock_id( '   Moist Processes: Donner_MCA'  , grain=CLOCK_MODULE )
   ras_clock        = mpp_clock_id( '   Moist Processes: RAS'         , grain=CLOCK_MODULE )
   closure_clock    = mpp_clock_id( '   Moist Processes: conv_closure', grain=CLOCK_MODULE )
   shallowcu_clock  = mpp_clock_id( '   Moist Processes: Shallow_cu'  , grain=CLOCK_MODULE )
   cmt_clock        = mpp_clock_id( '   Moist Processes: CMT'         , grain=CLOCK_MODULE )
   lscalecond_clock = mpp_clock_id( '   Moist Processes: lscale_cond' , grain=CLOCK_MODULE )
   stratcloud_clock = mpp_clock_id( '   Moist Processes: Strat_cloud' , grain=CLOCK_MODULE )
   bm_clock         = mpp_clock_id( '   Moist Processes: Betts-Miller', grain=CLOCK_MODULE )

       
      nbcphobic = get_tracer_index(MODEL_ATMOS,'bcphob')
      nbcphilic = get_tracer_index(MODEL_ATMOS,'bcphil')
      nomphobic = get_tracer_index(MODEL_ATMOS,'omphob')
      nomphilic = get_tracer_index(MODEL_ATMOS,'omphil')
      ndust1    = get_tracer_index(MODEL_ATMOS,'dust1')
      ndust2    = get_tracer_index(MODEL_ATMOS,'dust2')
      ndust3    = get_tracer_index(MODEL_ATMOS,'dust3')
      ndust4    = get_tracer_index(MODEL_ATMOS,'dust4')
      ndust5    = get_tracer_index(MODEL_ATMOS,'dust5')
      nsalt1 = get_tracer_index(MODEL_ATMOS,'ssalt1')
      nsalt2 = get_tracer_index(MODEL_ATMOS,'ssalt2')
      nsalt3 = get_tracer_index(MODEL_ATMOS,'ssalt3')
      nsalt4 = get_tracer_index(MODEL_ATMOS,'ssalt4')
      nsalt5 = get_tracer_index(MODEL_ATMOS,'ssalt5')

      nDMS      = get_tracer_index(MODEL_ATMOS,'simpleDMS')
      if (nDMS == NO_TRACER) then
        nDMS      = get_tracer_index(MODEL_ATMOS,'dms')
      endif

      nSO2      = get_tracer_index(MODEL_ATMOS,'simpleSO2')
      if (nSO2 == NO_TRACER) then
        nSO2      = get_tracer_index(MODEL_ATMOS,'so2')
      endif

      nSO4      = get_tracer_index(MODEL_ATMOS,'simpleSO4')
      if (nSO4 == NO_TRACER) then
        nSO4      = get_tracer_index(MODEL_ATMOS,'so4')
      endif

      nSOA      = get_tracer_index(MODEL_ATMOS,'SOA')
      nNH4NO3   = get_tracer_index(MODEL_ATMOS,'nh4no3')
      nNH4      = get_tracer_index(MODEL_ATMOS,'nh4')
!---------------------------------------------------------------------
!    retrieve the number of registered tracers in order to determine 
!    which tracers are to be convectively transported.
!---------------------------------------------------------------------
      call get_number_tracers (MODEL_ATMOS, num_prog= num_tracers)
 
!---------------------------------------------------------------------
!    allocate logical arrays to indicate the tracers which are to be
!    transported by the various available convective schemes. 
!    initialize these arrays to .false..
!---------------------------------------------------------------------
      allocate (tracers_in_donner(num_tracers))
      allocate (tracers_in_mca(num_tracers))
      allocate (tracers_in_ras(num_tracers))
      allocate (tracers_in_uw(num_tracers))
      tracers_in_donner = .false.
      tracers_in_mca = .false.
      tracers_in_ras = .false.
      tracers_in_uw = .false.

!----------------------------------------------------------------------
!    for each tracer, determine if it is to be transported by convect-
!    ion, and the convection schemes that are to transport it. set a 
!    logical flag to .true. for each tracer that is to be transported by
!    each scheme and increment the count of tracers to be transported
!    by that scheme.
!----------------------------------------------------------------------
      do n=1, num_tracers
        if (query_method ('convection', MODEL_ATMOS, n, scheme)) then
          select case (scheme)
            case ("none")
            case ("donner")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
            case ("mca")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
            case ("ras")
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
            case ("uw")
               num_uw_tracers = num_uw_tracers + 1
               tracers_in_uw(n) = .true.
            case ("donner_and_ras")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
            case ("donner_and_mca")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
            case ("mca_and_ras")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
            case ("all")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
               num_uw_tracers = num_uw_tracers + 1
               tracers_in_uw(n) = .true.
            case ("all_nodonner")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
               num_uw_tracers = num_uw_tracers + 1
               tracers_in_uw(n) = .true.
            case default  ! corresponds to "none"
          end select
        endif
      end do

!--------------------------------------------------------------------
!    set a logical indicating if any tracers are to be transported by
!    each of the available convection parameterizations.
!--------------------------------------------------------------------
      if (num_donner_tracers > 0) then
        do_tracers_in_donner = .true.
      else
        do_tracers_in_donner = .false.
      endif
      if (num_mca_tracers > 0) then
        do_tracers_in_mca = .true.
      else
        do_tracers_in_mca = .false.
      endif
      if (num_ras_tracers > 0) then
        do_tracers_in_ras = .true.
      else
        do_tracers_in_ras = .false.
      endif
      if (num_uw_tracers > 0) then
        do_tracers_in_uw = .true.
      else
        do_tracers_in_uw = .false.
      endif     
     
!---------------------------------------------------------------------
!    check for proper use of do_unified_convective_closure.
!---------------------------------------------------------------------
      if (do_unified_convective_closure) then
        call error_mesg ('moist_processes_init', &
         'do_unified_convective_closure is currently not allowed &
               & - see rsh', FATAL)
      endif
      if (do_unified_convective_closure) then
        if (.not. (do_donner_deep) .or. .not. (do_uw_conv)   &
            .or. do_ras .or. do_mca ) then
          call error_mesg ('moist_processes_init',  &
             'must have only donner_deep and uw shallow activated &
                &when do_unified_convective_closure is .true.', FATAL)
         endif
      endif
        
!---------------------------------------------------------------------
!    allocate and initialize arrays to hold maximum enthalpy and water
!    imbalances in each column.
!---------------------------------------------------------------------
      allocate (max_enthalpy_imbal (id, jd))
      allocate (max_water_imbal (id, jd))
      max_enthalpy_imbal = 0.
      max_water_imbal = 0.


!--------------------------------------------------------------------
!    initialize the convection scheme modules.
!--------------------------------------------------------------------
      if (do_donner_deep) then
        call get_time (Time, secs, days)
        call donner_deep_init (lonb, latb, pref, axes, secs, days,  &
                               tracers_in_donner,  &
                               do_donner_conservation_checks, &
                               do_unified_convective_closure, using_fms)
        if (do_donner_conservation_checks) then
          allocate (max_enthalpy_imbal_don (id, jd))
          allocate (max_water_imbal_don (id, jd))
          max_enthalpy_imbal_don = 0.
          max_water_imbal_don = 0.
        endif
      endif ! (do_donner_deep)
 
      if (do_ras)  then
        call ras_init (do_strat, do_liq_num, axes,Time, tracers_in_ras)
      endif

      if (do_uw_conv) call uw_conv_init (do_strat, axes, Time, kd, &
                                          tracers_in_uw)

      if (do_mca .or. do_donner_deep)  then
        call  moist_conv_init (axes,Time, tracers_in_mca)
      endif
  
 
!----- initialize quantities for diagnostics output -----
 
      call diag_field_init ( axes, Time, num_tracers, num_donner_tracers )

      if (do_lin_cld_microphys) then
         if (.not. do_strat) call error_mesg ('moist_processes_init',  &
                    'must also activate do_strat when do_lin_cld_microphys is active', FATAL)
         if (do_liq_num) call error_mesg ('moist_processes_init',  &
                    'do_lin_cld_microphys cannot be active with prognostic droplet &
                   & scheme (do_liq_num)', FATAL)
         nqr = get_tracer_index (MODEL_ATMOS, 'rainwat')
         nqs = get_tracer_index (MODEL_ATMOS, 'snowwat')
         nqg = get_tracer_index (MODEL_ATMOS, 'graupel')
         call lin_cld_microphys_init (id, jd, kd, axes, Time)
         ktop = 1
         do k = 1, kd
            if (pref(k) > 10.E2) then
              ktop=k
              exit
            endif
         enddo
         if (mpp_pe() == mpp_root_pe()) &
               write(*,*) 'Top layer for lin_cld_microphys=', ktop, pref(ktop)
      endif

      num_uw_tracers_out = num_uw_tracers
      do_strat_out = do_strat

      call  detr_ice_num_init 

      module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine moist_processes_init

!#######################################################################
!---> h1g
!subroutine moist_processes_end
subroutine moist_processes_end( clubb_term_clock )
integer, intent (out), optional :: clubb_term_clock
!<--- h1g

      if( .not.module_is_initialized ) return


!----------------close various schemes-----------------

! ---> h1g, cjg
      if (do_strat) then
        if (do_clubb > 0) then
! ---> h1g, if CLUBB is in moist-process, CLUBB ends here
          if( do_clubb == 2) then
              call mpp_clock_begin ( clubb_term_clock )
              call clubb_end
              call mpp_clock_end ( clubb_term_clock )
          endif
! <--- h1g, if CLUBB is in moist-process, CLUBB ends here
          call MG_microp_3D_end
        else
          call strat_cloud_end
        end if
      end if
! <--- h1g, cjg

      call  detr_ice_num_end
      if (do_rh_clouds)   call   rh_clouds_end
      if (do_diag_clouds) call  diag_cloud_end
      if (do_donner_deep) call donner_deep_end
      if (do_cmt        ) call cu_mo_trans_end
      if (do_ras        ) call         ras_end
      if (do_uw_conv    ) call     uw_conv_end
      if (do_lin_cld_microphys) call lin_cld_microphys_end

      deallocate (max_water_imbal)
      deallocate (max_enthalpy_imbal)
      if (do_donner_deep .and. do_donner_conservation_checks) then
        deallocate (max_water_imbal_don)
        deallocate (max_enthalpy_imbal_don)
      endif

      module_is_initialized = .false.

!-----------------------------------------------------------------------

end subroutine moist_processes_end


!#######################################################################
! <SUBROUTINE NAME="moist_processes_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine moist_processes_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

  if (do_strat)       call strat_cloud_restart(timestamp)
  if (do_diag_clouds) call diag_cloud_restart(timestamp)
  if (do_donner_deep) call donner_deep_restart(timestamp)

end subroutine moist_processes_restart
! </SUBROUTINE> NAME="moist_processes_restart"


!#######################################################################

subroutine diag_field_init ( axes, Time, num_tracers, num_donner_tracers )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time
  integer, intent(in) :: num_donner_tracers
  integer, intent(in) :: num_tracers

  character(len=32) :: tracer_units, tracer_name
  character(len=128) :: diaglname
  integer, dimension(3) :: half = (/1,2,4/)
  integer   :: n, nn

!------------ initializes diagnostic fields in this module -------------

   if ( any((/do_bm,do_bmmass,do_bmomp/)) ) then
      id_qref = register_diag_field ( mod_name, &
        'qref', axes(1:3), Time, &
        'Adjustment reference specific humidity profile', &
        'kg/kg',  missing_value=missing_value               )

      id_tref = register_diag_field ( mod_name, &
        'tref', axes(1:3), Time, &
        'Adjustment reference temperature profile', &
        'K',  missing_value=missing_value                   )

      id_bmflag = register_diag_field (mod_name, &
         'bmflag', axes(1:2), Time, &
         'Betts-Miller flag', &
         'no units', missing_value=missing_value            )

      id_klzbs  = register_diag_field  (mod_name, &
         'klzbs', axes(1:2), Time, &
         'klzb', &
         'no units', missing_value=missing_value            )

      id_cape = register_diag_field ( mod_name, & 
        'cape', axes(1:2), Time, &
        'Convectively available potential energy',      'J/Kg')

      id_cin = register_diag_field ( mod_name, &
        'cin', axes(1:2), Time, &
        'Convective inhibition',                        'J/Kg')
   endif

   if ( do_bm ) then
      id_invtaubmt  = register_diag_field  (mod_name, &
         'invtaubmt', axes(1:2), Time, &
         'Inverse temperature relaxation time', &
         '1/s', missing_value=missing_value            )

      id_invtaubmq = register_diag_field  (mod_name, &
         'invtaubmq', axes(1:2), Time, &
         'Inverse humidity relaxation time', &
         '1/s', missing_value=missing_value            )
   end if  ! if ( do_bm )

   if (do_bmmass) then
      id_massflux = register_diag_field (mod_name, &
         'massflux', axes(1:3), Time, &
         'Massflux implied by temperature adjustment', &
         'm/s', missing_value=missing_value                 )
   end if  ! if ( do_bmmass )

   id_ras_precip = register_diag_field ( mod_name, &
     'ras_precip', axes(1:2), Time, &
    'Precipitation rate from ras ',       'kg/m2/s' )

   id_ras_freq = register_diag_field ( mod_name, &
     'ras_freq', axes(1:2), Time, &
    'frequency of precip from ras ',       'number' , &
         missing_value = missing_value                       )

   id_don_precip = register_diag_field ( mod_name, &
     'don_precip', axes(1:2), Time, &
    'Precipitation rate from donner ',       'kg/m2/s' )

   id_don_freq = register_diag_field ( mod_name, &
     'don_freq', axes(1:2), Time, &
    'frequency of precip from donner ',       'number', &
         missing_value = missing_value                       )

   id_lsc_precip = register_diag_field ( mod_name, &
     'lsc_precip', axes(1:2), Time, &
    'Precipitation rate from lsc ',       'kg/m2/s' )

   id_lsc_freq = register_diag_field ( mod_name, &
     'lsc_freq', axes(1:2), Time, &
    'frequency of precip from lsc ',       'number' , &
         missing_value = missing_value                       )

   id_uw_precip = register_diag_field ( mod_name, &
     'uw_precip', axes(1:2), Time, &
    'Precipitation rate from uw shallow',       'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_uw_snow = register_diag_field ( mod_name, &
     'uw_snow', axes(1:2), Time, &
    'Snow rate from uw shallow',       'kg/m2/s' , &
     interp_method = "conserve_order1" )

   id_uw_freq = register_diag_field ( mod_name, &
     'uw_freq', axes(1:2), Time, &
    'frequency of precip from uw shallow ',       'number' , &
         missing_value = missing_value                       )

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from convection ',    'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from convection ',  'kg/kg/s',  &
                        missing_value=missing_value               )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from convection ',   'kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from convection ','W/m2' )
   
   id_enth_conv_col = register_diag_field ( mod_name, &
     'enth_conv_col', axes(1:2), Time, &
     'Column enthalpy tendency from convection','W/m2' )
 
   id_wat_conv_col = register_diag_field ( mod_name, &
     'wat_conv_col', axes(1:2), Time, &
     'Column total water tendency from convection','kg(h2o)/m2/s' )

   id_enth_donner_col2 = register_diag_field ( mod_name, &
     'enth_donner_col2', axes(1:2), Time, &
     'column enthalpy tendency from Donner liq precip','W/m2' )
 
   id_enth_donner_col3 = register_diag_field ( mod_name, &
     'enth_donner_col3', axes(1:2), Time, &
      'Column enthalpy tendency from Donner frzn precip','W/m2' )
 
   id_enth_donner_col4 = register_diag_field ( mod_name, &
      'enth_donner_col4', axes(1:2), Time, &
     'Atmospheric column enthalpy tendency from Donner convection', &
                                                            'W/m2' )
 
   id_enth_donner_col5 = register_diag_field ( mod_name, &
      'enth_donner_col5', axes(1:2), Time, &
      'Column enthalpy tendency due to condensate xfer from Donner &
                                                &to lsc','W/m2' )

  id_enth_donner_col6 = register_diag_field ( mod_name, &
     'enth_donner_col6', axes(1:2), Time, &
      'Column enthalpy tendency from donner moisture  &
                     &conservation  adjustment','W/m2' )
 
   id_enth_donner_col7 = register_diag_field ( mod_name, &
      'enth_donner_col7', axes(1:2), Time, &
      'Precip adjustment needed to balance donner moisture  &
                                           &adjustment','kg(h2o)/m2/s' )

   id_enth_donner_col = register_diag_field ( mod_name, &
     'enth_donner_col', axes(1:2), Time, &
     'Column enthalpy imbalance from Donner convection','W/m2' )

   id_wat_donner_col = register_diag_field ( mod_name, &
     'wat_donner_col', axes(1:2), Time, &
  'Column total water tendency from Donner convection','kg(h2o)/m2/s' )

   id_enth_mca_donner_col = register_diag_field ( mod_name, &
     'enth_mca_donner_col', axes(1:2), Time, &
    'Column enthalpy imbalance from Donner MCA convection','W/m2' )

   id_wat_mca_donner_col = register_diag_field ( mod_name, &
     'wat_mca_donner_col', axes(1:2), Time, &
     'Column total water imbalance from Donner MCA convection', &
                                                'kg(h2o)/m2/s' )

   id_enth_uw_col = register_diag_field ( mod_name, &
     'enth_uw_col', axes(1:2), Time, &
     'Column enthalpy tendency from UW convection','W/m2' )
 
   id_wat_uw_col = register_diag_field ( mod_name, &
     'wat_uw_col', axes(1:2), Time, &
      'Column total water tendency from UW convection','kg(h2o)/m2/s' )

   id_scale_uw = register_diag_field ( mod_name, &
     'scale_uw', axes(1:2), Time, &
     'Scaling factor applied to UW convection tendencies','1' )
          
   id_scale_donner = register_diag_field ( mod_name, &
     'scale_donner', axes(1:2), Time, &
     'Scaling factor applied to UW convection tendencies','1' )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from convection ',       'kg(h2o)/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from convection ',       'kg(h2o)/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_tot  = register_diag_field ( mod_name, &
     'snow_tot ', axes(1:2), Time, &
     'Frozen precip rate from all sources',       'kg(h2o)/m2/s', &
      interp_method = "conserve_order1" )

   id_conv_freq = register_diag_field ( mod_name, &
     'conv_freq', axes(1:2), Time, &
    'frequency of convection ',       'number', &
     missing_value = missing_value                       )

   id_gust_conv = register_diag_field ( mod_name, &
     'gust_conv', axes(1:2), Time, &
    'Gustiness resulting from convection ',       'm/s' )

  id_conv_rain3d= register_diag_field ( mod_name, &
     'conv_rain3d', axes(half), Time, &
    'Rain fall rate from convection -3D ',       'kg(h2o)/m2/s' )

   id_conv_snow3d= register_diag_field ( mod_name, &
     'conv_snow3d', axes(half), Time, &
    'Snow fall rate from convection -3D',       'kg(h2o)/m2/s' )

   id_lscale_rain3d= register_diag_field ( mod_name, &
     'lscale_rain3d', axes(half), Time, &
    'Rain fall rate from lscale  -3D ',   'kg(h2o)/m2/s' )

   id_lscale_snow3d= register_diag_field ( mod_name, &
     'lscale_snow3d', axes(half), Time, &
    'Snow fall rate from lscale -3D',       'kg(h2o)/m2/s' )
   
   id_lscale_precip3d= register_diag_field ( mod_name, &
     'lscale_precip3d', axes(1:3), Time, &
     'LS Precip falling out of gridbox',       'kg(h2o)/m2/s' , &
      mask_variant = .true., missing_value = missing_value)

    id_max_enthalpy_imbal    = register_diag_field    &
       (mod_name, 'max_enth_imbal', axes(1:2), Time,  &
        'max enthalpy  imbalance from moist_processes  ', 'W/m2',   &
              missing_value=missing_value)
    id_max_water_imbal    = register_diag_field    &
         (mod_name, 'max_water_imbal', axes(1:2), Time,   &
      'max water  imbalance from moist_processes  ', 'kg(h2o)/m2/s',  &
              missing_value=missing_value)

    id_enth_moist_col = register_diag_field ( mod_name, &
     'enth_moist_col', axes(1:2), Time, &
     'Column enthalpy imbalance from moist processes','W/m2' )
  
    id_wat_moist_col = register_diag_field ( mod_name, &
      'wat_moist_col', axes(1:2), Time, &
      'Column total water imbalance from moist processes','kg/m2/s' )

    if (do_donner_conservation_checks) then
      id_enthint    = register_diag_field    &
            (mod_name, 'enthint_don', axes(1:2), Time,  &
          'atmospheric column enthalpy change from donner', 'W/m2',  &
          missing_value=missing_value)
     id_lcondensint    = register_diag_field    &
         (mod_name, 'lcondensint_don', axes(1:2), Time, &
         'enthalpy transferred by condensate from donner to lscale', &
            'W/m2',  missing_value=missing_value)
     id_lprcp    = register_diag_field    &
             (mod_name, 'lprcpint_don', axes(1:2),   &
              Time, 'enthalpy removed by donner precip', 'W/m2',   &
             missing_value=missing_value)
      id_vertmotion    = register_diag_field    &
             (mod_name, 'vertmotion_don', axes(1:2), Time,  &
           'enthalpy change due to cell and meso motion in donner',  &
             'W/m2', missing_value=missing_value)
     id_enthdiffint    = register_diag_field    &
            (mod_name, 'enthdiffint_don', axes(1:2),   &
             Time, 'enthalpy  imbalance due to donner', 'W/m2',   &
             missing_value=missing_value)
     id_vaporint    = register_diag_field    &
           (mod_name, 'vaporint_don', axes(1:2),   &
            Time, 'column water vapor change', 'kg(h2o)/m2/s',   &
            missing_value=missing_value)
     id_max_enthalpy_imbal_don    = register_diag_field    &
            (mod_name, 'max_enth_imbal_don', axes(1:2),   &
              Time, 'max enthalpy  imbalance from donner', 'W/m**2',  &
              missing_value=missing_value)
     id_max_water_imbal_don    = register_diag_field    &
            (mod_name, 'max_water_imbal_don', axes(1:2),   &
              Time, 'max water imbalance from donner', 'kg(h2o)/m2/s', &
         missing_value=missing_value)
     id_condensint    = register_diag_field    &
           (mod_name, 'condensint_don', axes(1:2), Time,  &
         'column condensate exported from donner to lscale', &
                         'kg(h2o)/m2/s',  missing_value=missing_value )
     id_precipint    = register_diag_field    &
            (mod_name, 'precipint_don', axes(1:2),   &
             Time, 'column precip from donner', 'kg(h2o)/m2/s',   &
              missing_value=missing_value)
     id_diffint    = register_diag_field    &
          (mod_name, 'diffint_don', axes(1:2),   &
            Time, 'water imbalance due to donner', 'kg(h2o)/m2/s',   &
              missing_value=missing_value)
  endif



if (do_strat ) then

   id_qldt_conv = register_diag_field ( mod_name, &
     'qldt_conv', axes(1:3), Time, &
     'Liquid water tendency from convection',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qndt_conv = register_diag_field ( mod_name, &
     'qndt_conv', axes(1:3), Time, &
     'Liquid drop tendency from convection',      '#/kg/s',  &
                         missing_value=missing_value               )

   id_qidt_conv = register_diag_field ( mod_name, &
     'qidt_conv', axes(1:3), Time, &
     'Ice water tendency from convection',         'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_conv = register_diag_field ( mod_name, &
     'qadt_conv', axes(1:3), Time, &
     'Cloud fraction tendency from convection',    '1/sec',    &
                        missing_value=missing_value               )

   id_ql_conv_col = register_diag_field ( mod_name, &
     'ql_conv_col', axes(1:2), Time, &
    'Liquid water path tendency from convection',  'kg/m2/s' )
   
   id_qn_conv_col = register_diag_field ( mod_name, &
     'qn_conv_col', axes(1:2), Time, &
     'Liquid drp tendency from convection',  'kg/m2/s' )
 
   id_qi_conv_col = register_diag_field ( mod_name, &
     'qi_conv_col', axes(1:2), Time, &
    'Ice water path tendency from convection',     'kg/m2/s' )
   
   id_qa_conv_col = register_diag_field ( mod_name, &
     'qa_conv_col', axes(1:2), Time, &
    'Cloud mass tendency from convection',         'kg/m2/s' )
      
   id_qnidt_conv = register_diag_field ( mod_name, &
     'qnidt_conv', axes(1:3), Time, &
     'Ice number tendency from convection',      '#/kg/s',  &
                         missing_value=missing_value               )

   id_qni_conv_col = register_diag_field ( mod_name, &
     'qni_conv_col', axes(1:2), Time, &
     'Ice number tendency from convection',  'kg/m2/s' )

endif

if ( do_lsc ) then

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
       'Temperature tendency from large-scale cond',   'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from large-scale cond', 'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from large-scale cond',     'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from large-scale cond',     'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from large-scale cond','kg/m2/s' )
   
   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from large-scale cond','W/m2' )
   
 endif

   id_conv_cld_base = register_diag_field ( mod_name, &
     'conv_cld_base', axes(1:2), Time, &
     'pressure at convective cloud base',   'Pa', &
                       mask_variant = .true., &
                       missing_value=missing_value               )

   id_conv_cld_top = register_diag_field ( mod_name, &
     'conv_cld_top', axes(1:2), Time, &
     'pressure at convective cloud top',   'Pa', &
                       mask_variant = .true., &
                       missing_value=missing_value               )

if ( do_strat ) then

   id_mc_full = register_diag_field ( mod_name, &
     'mc_full', axes(1:3), Time, &
     'Net Mass Flux from convection',   'kg/m2/s', &
                       missing_value=missing_value               )
   
   id_mc_half = register_diag_field ( mod_name, &
     'mc_half', axes(half), Time, &
     'Net Mass Flux from convection on half levs',   'kg/m2/s', &
                       missing_value=missing_value               )
   
   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
     'Temperature tendency from strat cloud',        'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from strat cloud',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from strat cloud',          'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from strat cloud',          'kg/m2/s' )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from strat cloud',   'kg/m2/s' )
   
   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from strat cloud','W/m2' )
   
   id_qldt_ls = register_diag_field ( mod_name, &
     'qldt_ls', axes(1:3), Time, &
     'Liquid water tendency from strat cloud',       'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qndt_ls = register_diag_field ( mod_name, &
     'qndt_ls', axes(1:3), Time, &
     'Drop number tendency from strat cloud',        '#/kg/s',  &
                         missing_value=missing_value               )
   id_qidt_ls = register_diag_field ( mod_name, &
     'qidt_ls', axes(1:3), Time, &
     'Ice water tendency from strat cloud',          'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qnidt_ls = register_diag_field ( mod_name, &
     'qnidt_ls', axes(1:3), Time, &
     'Ice number tendency from strat cloud',          '#/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_ls = register_diag_field ( mod_name, &
     'qadt_ls', axes(1:3), Time, &
     'Cloud fraction tendency from strat cloud',     '1/sec',    &
                        missing_value=missing_value               )

   id_ql_ls_col = register_diag_field ( mod_name, &
     'ql_ls_col', axes(1:2), Time, &
    'Liquid water path tendency from strat cloud',   'kg/m2/s' )
   
   id_qn_ls_col = register_diag_field ( mod_name, &
     'qn_ls_col', axes(1:2), Time, &
     'Column drop number tendency from strat cloud',  '#/m2/s' )

   id_qni_ls_col = register_diag_field ( mod_name, &
     'qni_ls_col', axes(1:2), Time, &
     'Column ice particle number tendency from strat cloud',  '#/m2/s' )

   id_qi_ls_col = register_diag_field ( mod_name, &
     'qi_ls_col', axes(1:2), Time, &
    'Ice water path tendency from strat cloud',      'kg/m2/s' )
   
   id_qa_ls_col = register_diag_field ( mod_name, &
     'qa_ls_col', axes(1:2), Time, &
    'Cloud mass tendency from strat cloud',          'kg/m2/s' )
      
   id_enth_ls_col = register_diag_field ( mod_name, &
     'enth_ls_col', axes(1:2), Time, &
     'Column enthalpy tendency from strat cloud','W/m2' )
 
   id_wat_ls_col = register_diag_field ( mod_name, &
     'wat_ls_col', axes(1:2), Time, &
     'Column total water tendency from strat cloud','kg/m2/s' )

endif

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s', &
      interp_method = "conserve_order1" )

   id_WVP = register_diag_field ( mod_name, &
     'WVP', axes(1:2), Time, &
        'Column integrated water vapor',                'kg/m2'  )

if ( do_strat ) then

   id_LWP = register_diag_field ( mod_name, &
     'LWP', axes(1:2), Time, &
        'Liquid water path',                            'kg/m2'   )

   id_IWP = register_diag_field ( mod_name, &
     'IWP', axes(1:2), Time, &
        'Ice water path',                               'kg/m2'   )

   id_AWP = register_diag_field ( mod_name, &
     'AWP', axes(1:2), Time, &
        'Column integrated cloud mass ',                'kg/m2'   )

    id_tot_cld_amt = register_diag_field    &
              (mod_name, 'cld_amt_2d', axes(1:2), Time, &
                'total cloud amount', 'percent')

    id_tot_cloud_area = register_diag_field ( mod_name, &
      'tot_cloud_area', axes(1:3), Time, &
      'Cloud area -- all clouds', 'percent', missing_value=missing_value )

    id_tot_h2o     = register_diag_field ( mod_name, &
      'tot_h2o', axes(1:3), Time, &
      'total h2o -- all phases', 'kg/kg', missing_value=missing_value)

    id_tot_vapor     = register_diag_field ( mod_name, &
       'tot_vapor', axes(1:3), Time, &
       'total vapor', 'kg/kg', missing_value=missing_value)

    id_tot_liq_amt = register_diag_field ( mod_name, &
      'tot_liq_amt', axes(1:3), Time, &
      'Liquid amount -- all clouds', 'kg/kg', missing_value=missing_value)

    id_tot_ice_amt = register_diag_field ( mod_name, &
      'tot_ice_amt', axes(1:3), Time, &
      'Ice amount -- all clouds', 'kg/kg', missing_value=missing_value )

    id_lsc_cloud_area = register_diag_field ( mod_name, &
      'lsc_cloud_area', axes(1:3), Time, &
      'Large-scale cloud area', 'percent', missing_value=missing_value )

    id_lsc_liq_amt = register_diag_field ( mod_name, &
      'lsc_liq_amt', axes(1:3), Time, &
      'Large-scale cloud liquid amount', 'kg/kg', missing_value=missing_value )

    id_lsc_ice_amt = register_diag_field ( mod_name, &
      'lsc_ice_amt', axes(1:3), Time, &
      'Large-scale cloud ice amount', 'kg/kg', missing_value=missing_value )

    id_conv_cloud_area = register_diag_field ( mod_name, &
      'conv_cloud_area', axes(1:3), Time, &
      'Convective cloud area', 'percent', missing_value=missing_value )

    id_conv_liq_amt = register_diag_field ( mod_name, &
      'conv_liq_amt', axes(1:3), Time, &
      'Convective cloud liquid amount', 'kg/kg', missing_value=missing_value )

    id_conv_ice_amt = register_diag_field ( mod_name, &
      'conv_ice_amt', axes(1:3), Time, &
      'Convective cloud ice amount', 'kg/kg', missing_value=missing_value)
 
    id_WP_all_clouds = register_diag_field ( mod_name, &
      'WP_all_clouds', axes(1:2), Time, &
      'Total  water path -- all clouds',              'kg/m2'   )

    id_LWP_all_clouds = register_diag_field ( mod_name, &
      'LWP_all_clouds', axes(1:2), Time, &
      'Liquid water path -- all clouds',              'kg/m2'   )

    id_IWP_all_clouds = register_diag_field ( mod_name, &
      'IWP_all_clouds', axes(1:2), Time, &
      'Ice water path -- all clouds',                 'kg/m2'   )

endif

   id_tdt_dadj = register_diag_field ( mod_name, &
     'tdt_dadj', axes(1:3), Time, &
   'Temperature tendency from dry conv adj',       'deg_K/s',  &
                        missing_value=missing_value               )

   id_rh = register_diag_field ( mod_name, &
     'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  & 
                        missing_value=missing_value               )

   id_rh_cmip = register_diag_field ( mod_name, &
     'rh_cmip', axes(1:3), Time, &
     'relative humidity',                            'percent',  &
                      missing_value=missing_value               )

   id_qs = register_diag_field ( mod_name, &
     'qs', axes(1:3), Time, &
         'saturation specific humidity',                 'kg/kg',    & 
                        missing_value=missing_value               )
   
if (do_donner_deep) then

   id_tdt_deep_donner= register_diag_field ( mod_name, &
           'tdt_deep_donner', axes(1:3), Time, &
           ' heating rate - deep portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_deep_donner = register_diag_field ( mod_name, &
           'qdt_deep_donner', axes(1:3), Time, &
           ' moistening rate - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qadt_deep_donner = register_diag_field ( mod_name, &
     'qadt_deep_donner', axes(1:3), Time, &
     ' cloud amount tendency - deep portion', '1/s', &
                        missing_value=missing_value               )

   id_qldt_deep_donner = register_diag_field ( mod_name, &
     'qldt_deep_donner', axes(1:3), Time, &
     ' cloud liquid tendency - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qidt_deep_donner = register_diag_field ( mod_name, &
     'qidt_deep_donner', axes(1:3), Time, &
     ' ice water tendency - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )
   if (do_liq_num) &
    id_qndt_deep_donner = register_diag_field ( mod_name, &
            'qndt_deep_donner', axes(1:3), Time, &
            'deep convection cloud drop tendency', '#/kg/s', &
                       missing_value=missing_value               )

   if (do_ice_num) &
     id_qnidt_deep_donner = register_diag_field ( mod_name, &
      'qnidt_deep_donner', axes(1:3), Time, &
     ' ice number tendency - deep portion', '#/kg/s', &
                         missing_value=missing_value               )

   id_tdt_mca_donner = register_diag_field ( mod_name, &
     'tdt_mca_donner', axes(1:3), Time, &
     ' heating rate - mca  portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_mca_donner = register_diag_field ( mod_name, &
           'qdt_mca_donner', axes(1:3), Time, &
           ' moistening rate - mca  portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_prec_deep_donner = register_diag_field ( mod_name, &
           'prc_deep_donner', axes(1:2), Time, &
           ' total precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_prec1_deep_donner = register_diag_field ( mod_name, &
           'prc1_deep_donner', axes(1:2), Time, &
           ' change in precip for conservation in donner', 'kg/m2/s ', &
              missing_value=missing_value, mask_variant = .true., &
             interp_method = "conserve_order1"  )

   id_prec_mca_donner = register_diag_field ( mod_name, &
           'prc_mca_donner', axes(1:2), Time, &
           ' total precip rate - mca  portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_snow_deep_donner = register_diag_field ( mod_name, &
           'snow_deep_donner', axes(1:2), Time, &
           ' frozen precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_snow_mca_donner = register_diag_field ( mod_name, &
           'snow_mca_donner', axes(1:2), Time, &
           ' frozen precip rate -  mca portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_mc_donner = register_diag_field ( mod_name, &
           'mc_donner', axes(1:3), Time, &
           'Net Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_mc_donner_half = register_diag_field ( mod_name, &
           'mc_donner_half', axes(half), Time, &
           'Net Mass Flux from donner at half levs',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_mc_conv_up = register_diag_field ( mod_name, &
           'mc_conv_up', axes(1:3), Time, &
          'Upward Mass Flux from convection',   'kg/m2/s', &
                       missing_value=missing_value               )

   id_m_cdet_donner = register_diag_field ( mod_name, &
           'm_cdet_donner', axes(1:3), Time, &
           'Detrained Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_m_cellup = register_diag_field ( mod_name, &
           'm_cellup', axes(half), Time, &
           'Upward Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

! ---> h1g, cell and meso-scale cloud fraction from donner deep, 2011-08-08
   id_cell_cld_frac = register_diag_field ( mod_name, &
           'cell_cld_frac', axes(1:3), Time, & 
           'cell cloud fraction from donner',   '', &
                        missing_value=missing_value               )

   id_meso_cld_frac = register_diag_field ( mod_name, &
           'meso_cld_frac', axes(1:3), Time, & 
           'meso-scale cloud fraction from donner',   '', &
                        missing_value=missing_value               )

   id_donner_humidity_area = register_diag_field ( mod_name, &
           'donner_humidity_area', axes(1:3), Time, &
           'donner humidity area',  '', &
                        missing_value=missing_value               )
! <--- h1g, cell and meso-scale cloud fraction from donner deep, 2011-08-08


endif


if (do_uw_conv) then

   id_tdt_uw = register_diag_field ( mod_name, &
           'tdt_uw', axes(1:3), Time, &
           'UW convection heating rate', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_uw = register_diag_field ( mod_name, &
           'qdt_uw', axes(1:3), Time, &
           'UW convection moistening rate', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qadt_uw = register_diag_field ( mod_name, &
           'qadt_uw', axes(1:3), Time, &
           'UW convection cloud amount tendency', '1/s', &
                        missing_value=missing_value               )

   id_qldt_uw = register_diag_field ( mod_name, &
           'qldt_uw', axes(1:3), Time, &
           'UW convection cloud liquid tendency', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qidt_uw = register_diag_field ( mod_name, &
           'qidt_uw', axes(1:3), Time, &
           'UW convection ice water tendency', 'kg/kg/s', &
                        missing_value=missing_value               )

   if (do_liq_num) &
    id_qndt_uw = register_diag_field ( mod_name, &
           'qndt_uw', axes(1:3), Time, &
           'UW convection cloud drop tendency', '#/kg/s', &
                        missing_value=missing_value               )

    if (do_ice_num) &
     id_qnidt_uw = register_diag_field ( mod_name, &
           'qnidt_uw', axes(1:3), Time, &
           'UW convection ice number tendency', '#/kg/s', &
                        missing_value=missing_value               )

endif

   id_qvout = register_diag_field ( mod_name, &
           'qvout', axes(1:3), Time, 'qv after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   id_qaout = register_diag_field ( mod_name, &
           'qaout', axes(1:3), Time, 'qa after strat_cloud', 'none', &
                        missing_value=missing_value               )

   id_qlout = register_diag_field ( mod_name, &
           'qlout', axes(1:3), Time, 'ql after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   id_qiout = register_diag_field ( mod_name, &
           'qiout', axes(1:3), Time, 'qi after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   if (do_liq_num) then
   id_qnout = register_diag_field ( mod_name, &
           'qnout', axes(1:3), Time, 'qn after strat_cloud', '#/kg', &
                        missing_value=missing_value               )
   endif

   if (do_ice_num) then
   id_qniout = register_diag_field ( mod_name, &
           'qniout', axes(1:3), Time, 'qni after strat_cloud', '#/kg', &
                        missing_value=missing_value               )
   endif

!---------------------------------------------------------------------
!    register diagnostics for lightning NOx
!---------------------------------------------------------------------

   if (get_tracer_index(MODEL_ATMOS,'no') > 0) &
     id_prod_no = register_diag_field ( 'tracers', &
             'hook_no', axes(1:3), Time, &
             'hook_no',   'molec/cm3/s')

!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!    register the diagnostics associated with convective tracer 
!    transport.
!---------------------------------------------------------------------
      allocate (id_tracerdt_conv    (num_tracers))
      allocate (id_tracerdt_conv_col(num_tracers))
      allocate (id_wet_deposition(num_tracers))
      allocate (id_wetdep       (num_tracers))
      allocate (id_conv_tracer           (num_tracers))
      allocate (id_conv_tracer_col(num_tracers))

      id_tracerdt_conv = -1
      id_tracerdt_conv_col = -1
      id_wet_deposition = -1
      id_wetdep = -1
      id_conv_tracer = -1
      id_conv_tracer_col = -1
      
 
      id_wetdep_om = &
                         register_diag_field ( mod_name, &
                         'om_wet_dep',  &
                         axes(1:2), Time,  &
                         'total om wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_SOA = &
                         register_diag_field ( mod_name, &
                         'SOA_wet_dep',  &
                         axes(1:2), Time,  &
                         'total SOA wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_bc = &
                         register_diag_field ( mod_name, &
                         'bc_wet_dep',  &
                         axes(1:2), Time,  &
                         'total bc wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_so4 = &
                         register_diag_field ( mod_name, &
                         'so4_wet_dep',  &
                         axes(1:2), Time,  &
                         'total so4 wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_so2 = &
                         register_diag_field ( mod_name, &
                         'so2_wet_dep',  &
                         axes(1:2), Time,  &
                         'total so2 wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_DMS = &
                         register_diag_field ( mod_name, &
                         'DMS_wet_dep',  &
                         axes(1:2), Time,  &
                         'total DMS wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_NH4NO3 =  &
                         register_diag_field ( mod_name, &
                         'totNH4_wet_dep',  &
                         axes(1:2), Time,  &
                         'total NH4 + NH3 wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_salt   =  &
                         register_diag_field ( mod_name, &
                         'ssalt_wet_dep',  &
                         axes(1:2), Time,  &
                         'total seasalt wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_wetdep_dust   =  &
                         register_diag_field ( mod_name, &
                         'dust_wet_dep',  &
                         axes(1:2), Time,  &
                         'total dust wet deposition', &
                         'kg/m2/s',  &
                         missing_value=missing_value)

      id_f_snow_berg   =  &
                         register_diag_field ( mod_name, &
                         'f_snow_berg',  &
                         axes(1:3), Time,  &
                         'fraction of snow/ice produced having IFN', &
                         'fraction',  &
                         missing_value=missing_value)

      id_f_snow_berg_cond   =  &
                         register_diag_field ( mod_name, &
                         'f_snow_berg_cond',  &
                         axes(1:3), Time,  &
                         'conditional fraction of snow/ice produced &
                         &having IFN', 'fraction',  &
                         mask_variant = .true., &
                         missing_value=missing_value)

      id_f_snow_berg_wtd   =  &
                         register_diag_field ( mod_name, &
                         'f_snow_berg_wtd',  &
                         axes(1:3), Time,  &
                         'product of snow/ice produced having IFN and &
                         &ls precip falling out of gridbox', &
                         'kg(h2o)/m2/s', mask_variant = .true.,   &
                         missing_value=missing_value)

      do n = 1,num_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        if (tracers_in_donner(n) .or. &
            tracers_in_ras(n)      .or.  &
            tracers_in_mca(n)      .or.  &
            tracers_in_uw(n)) then

          diaglname = trim(tracer_name)//  &
                        ' wet deposition from all precip'
          id_wetdep(n) = &
                       register_diag_field ( mod_name, &
                          TRIM(tracer_name)//'_wet_depo',  &
                          axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'/s',  &
                          missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                        ' total tendency from moist convection'
          id_tracerdt_conv(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                       ' total path tendency from moist convection'
          id_tracerdt_conv_col(n) =  &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv_col', &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'*(kg/m2)/s',   &
                         missing_value=missing_value)
         endif
 
         diaglname = trim(tracer_name)
         id_conv_tracer(n) =    &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
         diaglname =  ' column integrated' // trim(tracer_name)
         id_conv_tracer_col(n) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_col', &
                        axes(1:2), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,   &
                        missing_value=missing_value)
         id_wet_deposition(n) = register_diag_field( mod_name, &
           trim(tracer_name)//'_wetdep', axes(1:3), Time, &
           trim(tracer_name)//' tendency from wet deposition',TRIM(tracer_units)//'/sec', &
           missing_value=missing_value )
      end do

!------------------------------------------------------------------
!    register the variables associated with the mca component of 
!    donner_deep transport.
!------------------------------------------------------------------
     if (do_donner_deep) then
       allocate (id_tracerdt_mcadon  (num_donner_tracers))
       allocate (id_tracerdt_mcadon_col(num_donner_tracers))
 
       nn = 1
       do n = 1,num_tracers
         call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                units = tracer_units)
         if (tracers_in_donner(n) ) then
           diaglname = trim(tracer_name)//  &
                       ' tendency from donner-mca'
           id_tracerdt_mcadon(nn) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_donmca',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                        missing_value=missing_value)

           diaglname = trim(tracer_name)//  &
                       ' total path tendency from donner-mca'
           id_tracerdt_mcadon_col(nn) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_donmca_col', &
                        axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'*(kg/m2)/s',   &
                        missing_value=missing_value)
           nn = nn + 1
         endif
       end do
 

     endif

end subroutine diag_field_init


!#######################################################################
function doing_strat()
logical :: doing_strat

  if (.not. module_is_initialized) call error_mesg ('doing_strat',  &
                     'moist_processes_init has not been called.', FATAL)

  doing_strat = do_strat

end function doing_strat


!#######################################################################  

subroutine set_cosp_precip_sources (cosp_precip_sources)

character(len=16),        intent(in) :: cosp_precip_sources

     if (trim(cosp_precip_sources)  == 'stratdeepuw') then
       strat_precip_in_cosp = 1.
       donner_precip_in_cosp = 1.
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'stratdeep') then
       strat_precip_in_cosp = 1.
       donner_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'stratuw') then
       strat_precip_in_cosp = 1.
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'deepuw') then
       donner_precip_in_cosp = 1.
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'strat') then
       strat_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'deep') then
       donner_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'uw') then
       uw_precip_in_cosp = 1.
     else if (trim(cosp_precip_sources)  == 'noprecip') then
!     COSP run without any precip input     
     else
       call error_mesg ('moist_processes_mod:set_cosp_precip_sources', &
        'cosp_precip_sources does not match any currently allowed string',&
                                                                 FATAL)
     endif

end subroutine set_cosp_precip_sources


!#######################################################################




end module moist_processes_mod

  

