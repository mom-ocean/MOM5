module moistproc_kernels_mod

use sat_vapor_pres_mod,         only: compute_qs
use time_manager_mod,           only: time_type
use diag_manager_mod,           only: send_data
use constants_mod,              only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                      RDGAS, RVGAS, TFREEZE, &
                                      SECONDS_PER_DAY, KAPPA
use field_manager_mod,          only: MODEL_ATMOS
use tracer_manager_mod,         only: get_tracer_index
use betts_miller_mod,           only: betts_miller
use bm_massflux_mod,            only: bm_massflux
use bm_omp_mod,                 only: bm_omp
use diag_cloud_mod,             only: diag_cloud_sum
use donner_deep_mod,            only: donner_deep
use moist_conv_mod,             only: moist_conv
use lscale_cond_mod,            only: lscale_cond
use uw_conv_mod,                only: uw_conv
use lin_cld_microphys_mod,      only: lin_cld_microphys_driver
use ras_mod,                    only: ras
use strat_cloud_mod,            only: strat_cloud, strat_cloud_sum, &
                                      strat_cloud_new
use rh_clouds_mod,              only: rh_clouds_sum
use cu_mo_trans_mod,            only: cu_mo_trans
use atmos_tracer_utilities_mod, only: wet_deposition
use moz_hook_mod,               only: moz_hook
use rad_utilities_mod,          only: aerosol_type
use moist_proc_utils_mod,       only: rh_calc
use detr_ice_num_mod ,          only: detr_ice_num

!--->h1g
use  mpp_mod,                   only: mpp_chksum, mpp_pe, mpp_root_pe
use  MG_microp_3D_mod,          only: MG_microp_3D
! <--- h1g

implicit none
private
public  moistproc_init, moistproc_end, moistproc_mca, moistproc_ras, &
        moistproc_lscale_cond, moistproc_strat_cloud, moistproc_cmt, &
        moistproc_uw_conv, moistproc_scale_uw, moistproc_scale_donner


!--------------------- version number ----------------------------------
character(len=128) :: &
version = '$Id: moistproc_kernels.F90,v 20.0 2013/12/13 23:18:29 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!-----------------------------------------------------------------------
real, public, allocatable, dimension(:,:)     :: rain_uw, snow_uw
real, public, allocatable, dimension(:,:,:)   :: ttnd_uw, qtnd_uw,   &
                                                 utnd_uw, vtnd_uw,   &
                                                 qltnd_uw, qitnd_uw, &
                                                 qatnd_uw, qntnd_uw, &  
                                                 qnitnd_uw,          &
                                                 delta_ql, delta_qi, &
                                                 delta_qa,           &
                                                 qlin, qiin, qain
real, public, allocatable, dimension(:,:,:,:) :: qtruw

logical :: moistproc_initialized = .false.

contains


!#######################################################################
subroutine moistproc_init(ix, jx, kx, num_uw_tracers, do_strat)
  integer, intent(in) :: ix, jx, kx, num_uw_tracers
  logical, intent(in) :: do_strat

      if (moistproc_initialized) return

      allocate( rain_uw  (ix,jx) )                   ; rain_uw  = 0.0
      allocate( snow_uw  (ix,jx) )                   ; snow_uw  = 0.0
      allocate( ttnd_uw  (ix,jx,kx) )                ; ttnd_uw  = 0.0
      allocate( qtnd_uw  (ix,jx,kx) )                ; qtnd_uw  = 0.0
      allocate( utnd_uw  (ix,jx,kx) )                ; utnd_uw  = 0.0
      allocate( vtnd_uw  (ix,jx,kx) )                ; vtnd_uw  = 0.0
      allocate( qltnd_uw (ix,jx,kx) )                ; qltnd_uw = 0.0
      allocate( qitnd_uw (ix,jx,kx) )                ; qitnd_uw = 0.0
      allocate( qatnd_uw (ix,jx,kx) )                ; qatnd_uw = 0.0
      allocate( qntnd_uw (ix,jx,kx) )                ; qntnd_uw = 0.0
      allocate( qnitnd_uw (ix,jx,kx) )               ; qnitnd_uw= 0.0
      allocate( qtruw    (ix,jx,kx,num_uw_tracers) ) ; qtruw    = 0.0
      if (do_strat) then
        allocate( delta_ql (ix,jx,kx) )              ; delta_ql = 0.0
        allocate( delta_qi (ix,jx,kx) )              ; delta_qi = 0.0
        allocate( delta_qa (ix,jx,kx) )              ; delta_qa = 0.0
        allocate( qlin     (ix,jx,kx) )              ; qlin     = 0.0
        allocate( qiin     (ix,jx,kx) )              ; qiin     = 0.0
        allocate( qain     (ix,jx,kx) )              ; qain     = 0.0
      endif

      moistproc_initialized = .true.
end subroutine moistproc_init

!#######################################################################
subroutine moistproc_end(do_strat)
  logical, intent(in) :: do_strat

      if (moistproc_initialized .eqv. .false. ) return
      deallocate( rain_uw    )
      deallocate( snow_uw    )
      deallocate( ttnd_uw    )
      deallocate( qtnd_uw    )
      deallocate( utnd_uw    )
      deallocate( vtnd_uw    )
      deallocate( qltnd_uw   )
      deallocate( qitnd_uw   )
      deallocate( qatnd_uw   )
      deallocate( qntnd_uw   )
      deallocate( qnitnd_uw   )
      deallocate( qtruw      )
      if (do_strat) then
        deallocate( delta_ql )
        deallocate( delta_qi )
        deallocate( delta_qa )
        deallocate( qlin     )
        deallocate( qiin     )
        deallocate( qain     )
      endif

      moistproc_initialized = .false.
end subroutine moistproc_end


!#######################################################################
subroutine moistproc_cmt ( Time, is, js, t, u, v, tracer, pfull, phalf, &
                           zfull, zhalf, pmass, tdt, udt, vdt, rdt,     &
                           ttnd_conv, dt, mc_cmt, det_cmt, diff_cu_mo,  &
                           num_tracers)
  type(time_type), intent(in)   :: Time
  integer, intent(in)           :: is, js, num_tracers
  real, intent(in)              :: dt
  real, intent(in),    dimension(:,:,:) :: pfull, phalf, zfull, zhalf, pmass, &
                                           mc_cmt, det_cmt
  real, intent(inout), dimension(:,:,:) :: t, u, v, tdt, udt, vdt, ttnd_conv, &
                                           diff_cu_mo
  real, intent(inout), dimension(:,:,:,:) :: rdt, tracer

  integer :: n
  real, dimension(size(t,1), size(t,2), size(t,3)) :: ttnd, utnd, vtnd
  real, dimension(size(rdt,1), size(rdt,2), size(rdt,3),num_tracers) :: qtr

      call cu_mo_trans (is, js, Time, mc_cmt, t, phalf, pfull, &
                        zhalf, zfull, dt, u, v, tracer,        &
                        pmass, det_cmt, utnd, vtnd, ttnd,      &
                        qtr, diff_cu_mo  )

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from cu_mo_trans.
!---------------------------------------------------------------------
      do n=1, num_tracers
        rdt(:,:,:,n) = rdt(:,:,:,n) + qtr(:,:,:,n)
      end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
      tdt = tdt + ttnd 
      udt = udt + utnd
      vdt = vdt + vtnd
      ttnd_conv = ttnd_conv + ttnd


end subroutine moistproc_cmt


!#######################################################################
subroutine moistproc_lscale_cond (is, js, t, q, pfull, phalf, tdt, qdt, &
                                  ttnd, qtnd, qtnd_conv, lprec, fprec, precip,&
                                  rain, snow, dtinv, omega, do_rh_clouds, do_simple, &
                                  do_diag_clouds, coldT, kbot, mask)
  integer, intent(in)  :: is, js
  real, intent(in)     :: dtinv
  logical, intent(in)  :: do_rh_clouds, do_simple, do_diag_clouds
  logical, intent(in), dimension(:,:)   :: coldT
  real, intent(in),    dimension(:,:,:) :: pfull, phalf, omega
  real, intent(inout), dimension(:,:)   :: lprec, fprec, precip
  real, intent(inout), dimension(:,:,:) :: t, q, tdt, qdt, qtnd_conv, &
                                           ttnd, qtnd
  real, intent(out),   dimension(:,:)   :: rain, snow
  integer, intent(in) , dimension(:,:), optional :: kbot
  real, intent(in) , dimension(:,:,:),  optional :: mask

  real, dimension(size(t,1), size(t,2), size(t,3)) :: cnvcntq, rh

      call lscale_cond (t, q, pfull, phalf, coldT, rain, snow,  &
                        ttnd, qtnd, mask=mask)

!-----------------------------------------------------------------------
!    add the temperature and specific humidity increments to the updated
!    temperature and specific humidity fields (tin, qin). convert these
!    increments and the precipitation increments to rates and add to 
!    the arrays accumulating the total rates for all physical processes
!    (tdt, qdt, lprec, fprec).
!-----------------------------------------------------------------------
      t     = t   + ttnd 
      q     = q   + qtnd
      tdt   = tdt   + ttnd*dtinv 
      qdt   = qdt   + qtnd*dtinv
      lprec = lprec + rain*dtinv
      fprec = fprec + snow*dtinv

!--------------------------------------------------------------------
!    if rh_clouds is active, call rh_calc to determine the grid box
!    relative humidity. call rh_clouds_sum to pass this field to 
!    rh_clouds_mod so it may be used to determine the grid boxes which
!    will contain clouds for the radiation package.
!---------------------------------------------------------------------
      if (do_rh_clouds) then
        call rh_calc (pfull, t, q, rh, do_simple, mask)
        call rh_clouds_sum (is, js, rh)
      endif

!--------------------------------------------------------------------
!    if the gordon diagnostic cloud parameterization is active, set a 
!    flag to indicate those grid points where drying has resulted from 
!    convective activity (cnvcntq). call rh_calc to determine the grid 
!    box relative humidity. call diag_cloud_sum to define the cloud 
!    field that will be seen by the radiation package.
!---------------------------------------------------------------------
      if (do_diag_clouds) then
        cnvcntq (:,:,:) = 0.0
        where (qtnd_conv(:,:,:) < 0.0)
          cnvcntq (:,:,:) = 1.0
        end where
        call rh_calc (pfull, t, q, rh, do_simple, mask)
        call diag_cloud_sum (is, js, t, q, rh, omega, qtnd,  &
                             cnvcntq, precip, kbot)
      endif


end subroutine moistproc_lscale_cond


!#######################################################################
subroutine moistproc_mca( Time, is, js, t, q, tracer, pfull, phalf, coldT, dtinv, &
                          tdt, qdt, rdt, q_tnd, ttnd_conv, qtnd_conv,             &
                          lprec, fprec, do_strat, num_tracers, tracers_in_mca,    &
                          num_mca_tracers, kbot, mask)

  type(time_type), intent(in) :: Time
  integer, intent(in)         :: is, js, num_tracers, num_mca_tracers
  real, intent(in)            :: dtinv
  logical, intent(in)         :: do_strat
  logical, intent(in), dimension(:)       :: tracers_in_mca
  logical, intent(in), dimension(:,:)     :: coldT
  real, intent(in),    dimension(:,:,:)   :: pfull, phalf 
  real, intent(inout), dimension(:,:)     :: lprec, fprec
  real, intent(inout), dimension(:,:,:)   :: t, q, tdt, qdt, ttnd_conv, qtnd_conv
  real, intent(inout), dimension(:,:,:,:) :: rdt, q_tnd
  real, intent(out),   dimension(:,:,:,:) :: tracer
  integer, intent(in) , dimension(:,:), optional :: kbot
  real, intent(in) , dimension(:,:,:),  optional :: mask

  integer :: nn, n, nql, nqa, nqi, nqn
  real, dimension(size(t,1), size(t,2)) :: rain, snow
  real, dimension(size(t,1), size(t,2), size(t,3)) :: ttnd, qtnd
  real, dimension(size(rdt,1), size(rdt,2), size(rdt,3),num_mca_tracers) :: trcr, qtr


      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )

!---------------------------------------------------------------------
!    check each active tracer to find any that are to be transported 
!    by moist convective adjustment and fill the mca_tracers array with
!    these fields.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_mca(n)) then
          trcr(:,:,:,nn) = tracer(:,:,:,n)
          nn = nn + 1
        endif
      end do

!---------------------------------------------------------------------
!    call subroutine moist_conv to obtain the temperature, moisture
!    precipitation and tracer tendencies due to the moist convective
!    adjustment parameterization. currently there is no tracer tendency
!    due to this parameterization.
!---------------------------------------------------------------------
!++++yim Should also account for change in qn dut to moist convective adjustment.

      if (do_strat) then
        call moist_conv (t, q, pfull, phalf, coldT, ttnd, qtnd,          &
                         rain, snow, dtinv, Time, is, js,                &
                         trcr, qtr, Lbot= kbot, mask=mask,               &
                         ql=tracer(:,:,:,nql), qi=tracer(:,:,:,nqi),     &
                         cf=tracer(:,:,:,nqa), qldel=q_tnd(:,:,:,nql),   &
                         qidel=q_tnd(:,:,:,nqi), cfdel=q_tnd(:,:,:,nqa))
      else
        call moist_conv (t, q, pfull, phalf, coldT, ttnd, qtnd,          &
                         rain, snow, dtinv, Time, is, js,                &
                         trcr, qtr, Lbot=kbot, mask=mask)
      endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from moist convective adjustment. currently there
!    is no tracer transport by this process.
!    NOTE : the stratcloud tracers are updated within moist_conv.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_mca(n)) then
          rdt(:,:,:,n) = rdt(:,:,:,n) + qtr(:,:,:,nn)
          nn = nn + 1
        endif
      end do

!----------------------------------------------------------------------
!    add the temperature and specific humidity tendencies from moist
!    convective adjustment (ttnd, qtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt).
!----------------------------------------------------------------------
      tdt = tdt + ttnd 
      qdt = qdt + qtnd
      ttnd_conv = ttnd_conv + ttnd
      qtnd_conv = qtnd_conv + qtnd

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from moist convective adjustment.
!----------------------------------------------------------------------
      lprec  = lprec  + rain
      fprec  = fprec  + snow

!----------------------------------------------------------------------
!    if strat_cloud_mod is activated, add the cloud liquid, ice and area
!    tendencies from moist convective adjustment to the 
!    arrays accumulating these tendencies from all physics processes 
!    (rdt).
!----------------------------------------------------------------------
      if (do_strat) then
        rdt(:,:,:,nql) = rdt(:,:,:,nql) + q_tnd(:,:,:,nql)
        rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + q_tnd(:,:,:,nqi)
        rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + q_tnd(:,:,:,nqa)
      endif


end subroutine moistproc_mca


!#######################################################################
subroutine moistproc_ras(Time, is, js, dt, coldT, t, q, u, v, tracer,       &
                         pfull, phalf, zhalf, tdt, qdt, udt, vdt, rdt,      &
                         q_tnd, ttnd, qtnd, ttnd_conv, qtnd_conv, mc, det0, &
                         lprec, fprec, rain, snow, rain3d, snow3d,          &
                         Aerosol, do_strat, do_liq_num, num_tracers,        &
                         tracers_in_ras, num_ras_tracers, kbot, mask, & 
                         do_ice_num, detrain_ice_num)

  type(time_type), intent(in)   :: Time
  integer, intent(in)           :: is, js, num_tracers, num_ras_tracers
  logical, intent(in)           :: do_strat, do_liq_num
  real, intent(in)              :: dt
  logical, intent(in), dimension(:)     :: tracers_in_ras
  logical, intent(in), dimension(:,:)   :: coldT
  real, intent(in),    dimension(:,:,:) :: pfull, phalf, zhalf
  real, intent(inout), dimension(:,:)   :: lprec, fprec
  real, intent(inout), dimension(:,:,:) :: t, q, u, v, tdt, qdt, udt, vdt,   &
                                           ttnd, qtnd, ttnd_conv, qtnd_conv
  real, intent(inout), dimension(:,:,:,:) :: rdt, tracer, q_tnd
  real, intent(out),   dimension(:,:)     :: rain, snow
  real, intent(out),   dimension(:,:,:)   :: rain3d,  snow3d, mc, det0

  type(aerosol_type),intent(in), optional :: Aerosol
  integer, intent(in), dimension(:,:), optional :: kbot
  real, intent(in), dimension(:,:,:),  optional :: mask
  logical, intent(in )                     :: do_ice_num, detrain_ice_num

  integer :: nn, n, nql, nqa, nqi, nqn, nqni
  real, dimension(size(t,1), size(t,2), size(t,3)) :: utnd, vtnd
  real, dimension(size(rdt,1), size(rdt,2), size(rdt,3),num_ras_tracers) :: trcr, qtr

      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
!----------------------------------------------------------------------
!    if any tracers are to be transported by ras convection, check each
!    active tracer to find those to be transported and fill the 
!    ras_tracers array with these fields.
!---------------------------------------------------------------------
     nn = 1
     do n=1, num_tracers
       if (tracers_in_ras(n)) then
         trcr(:,:,:,nn) = tracer(:,:,:,n)
         nn = nn + 1
       endif
     end do

!----------------------------------------------------------------------
!    call subroutine ras to obtain the temperature, specific humidity,
!    velocity, precipitation and tracer tendencies and mass flux 
!    associated with the relaxed arakawa-schubert parameterization.
!----------------------------------------------------------------------
     if (do_strat .and. (.not.do_liq_num)) then
       call ras (is,   js,     Time,     t,   q,          &
                 u,  v,    pfull,    phalf, zhalf, coldT, &
                 dt,   ttnd,   qtnd,     utnd,  vtnd,     &
                 rain3d, snow3d, rain, snow,              &
                 trcr, qtr, mask,  kbot, mc, det0,        &
                 tracer(:,:,:,nql), tracer(:,:,:,nqi),    &
                 tracer(:,:,:,nqa), q_tnd(:,:,:,nql),     &
                 q_tnd(:,:,:,nqi), q_tnd(:,:,:,nqa))       

     elseif (do_strat .and. do_liq_num) then
       call ras (is,   js,     Time,     t,   q,          &
                 u,  v,    pfull,    phalf, zhalf, coldT, &
                 dt,   ttnd,   qtnd,     utnd,  vtnd,     &
                 rain3d, snow3d, rain, snow,              &
                 trcr, qtr, mask,  kbot, mc, det0,        &
                 tracer(:,:,:,nql), tracer(:,:,:,nqi),    &
                 tracer(:,:,:,nqa), q_tnd(:,:,:,nql),     &
                 q_tnd(:,:,:,nqi), q_tnd(:,:,:,nqa),      &
                 tracer(:,:,:,nqn), q_tnd(:,:,:,nqn),     &
                 do_strat, Aerosol)
     else
       call ras (is,   js,     Time,     t,   q,          &
                 u,  v,    pfull,    phalf, zhalf, coldT, &
                 dt,   ttnd,   qtnd,     utnd,  vtnd,     &
                 rain3d, snow3d, rain, snow,              &
                 trcr, qtr, mask,  kbot,  mc, det0)
     endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from ras transport.
!    NOTE : the stratcloud tracers are updated within ras.        
!---------------------------------------------------------------------
     nn = 1
     do n=1, num_tracers
       if (tracers_in_ras(n)) then
         rdt(:,:,:,n) = rdt(:,:,:,n) + qtr (:,:,:,nn)
         nn = nn + 1
       endif
     end do
!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
     tdt = tdt + ttnd 
     qdt = qdt + qtnd
     udt = udt + utnd
     vdt = vdt + vtnd
!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
     ttnd_conv = ttnd_conv + ttnd
     qtnd_conv = qtnd_conv + qtnd

!----------------------------------------------------------------------
!    if strat_cloud_mod is activated, add the cloud liquid, ice and area
!    tendencies from ras to the arrays accumulating these tendencies 
!    from all physics processes (rdt).
!----------------------------------------------------------------------
     if (do_strat) then
       rdt(:,:,:,nql) = rdt(:,:,:,nql) + q_tnd(:,:,:,nql)
       rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + q_tnd(:,:,:,nqi)
       rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + q_tnd(:,:,:,nqa)
       if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + q_tnd(:,:,:,nqn)

!------------------------------------------------------------------------
!    currently this is the ice number tendency due to detrainment 
!    proportional by ice mass
!------------------------------------------------------------------------
       IF (do_ice_num .AND. detrain_ice_num) THEN
         CALL detr_ice_num (t, q_tnd(:,:,:,nqi), q_tnd(:,:,:,nqni))   
         rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + q_tnd(:,:,:,nqni)
       END IF
     endif

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from ras.
!----------------------------------------------------------------------
     lprec  = lprec  + rain
     fprec  = fprec  + snow

end subroutine moistproc_ras


!#######################################################################
subroutine moistproc_strat_cloud(Time, is, ie, js, je, lon, lat, ktop, dt, tm, t, q, tracer,&  ! cjg
                                 pfull, phalf, zhalf, omega, radturbten, mc_full, &
                                 diff_t, land, area, tdt, qdt, rdt, q_tnd, ttnd,  &
                                 qtnd, lprec, fprec, f_snow_berg, rain, &
                                 snow, rain3d, snow3d,  &
                                 snowclr3d, &
                                 Aerosol, lsc_cloud_area, lsc_liquid, lsc_ice,    &
                                 lsc_droplet_number, donner_humidity_area,        &
                                 donner_humidity_factor, shallow_cloud_area,      &
                                 cell_cld_frac, meso_cld_frac,                    &
                                 do_uw_conv, do_donner_deep, do_liq_num,          &
                                 do_clubb,                                        &  ! cjg
                                 do_lin_cld_microphys, id_qvout, id_qlout,        &
                                 id_qaout, id_qiout, id_qnout, id_qniout, &
                                 limit_conv_cloud_frac, mask, &
                                 hydrostatic, phys_hydrostatic,           &
                                 zfull, do_ice_num,  lsc_ice_number,   &
                                 lsc_snow, lsc_rain, lsc_snow_size,   &
                                 lsc_rain_size, do_legacy_strat_cloud, &
! ---> h1g
                                 dcond_ls_liquid, dcond_ls_ice,                   &
                                 Ndrop_act_CLUBB, Icedrop_act_CLUBB,              &
                                 ndust, rbar_dust, qcvar_clubb )
! <--- h1g

  type(time_type), intent(in) :: Time
  integer, intent(in)         :: is, ie, js, je, ktop, id_qvout, id_qlout, &
                                 id_qaout, id_qiout, id_qnout, id_qniout
  real, intent(in), dimension(:,:)  :: lon, lat
  real, intent(in)            :: dt
  logical, intent(in)         :: do_uw_conv, do_donner_deep, do_liq_num, &
                                 do_lin_cld_microphys,  &
                                 limit_conv_cloud_frac, &
                                 do_ice_num, do_legacy_strat_cloud
  integer, intent(in)         :: do_clubb                                 ! cjg
  real, intent(in),    dimension(:,:)     :: land, area
  real, intent(in),    dimension(:,:,:)   :: tm, pfull, phalf, zhalf, omega,  &
                                             radturbten, mc_full, diff_t,     &
                                             donner_humidity_area, donner_humidity_factor
  real, intent(in),    dimension(:,:,:)   :: zfull
  real, intent(inout), dimension(:,:)     :: lprec, fprec
  real, intent(inout), dimension(:,:,:)   :: t, q, tdt, qdt, ttnd, qtnd
  real, intent(inout), dimension(:,:,:,:) :: rdt, tracer, q_tnd
  real, intent(out),   dimension(:,:)     :: rain, snow
  real, intent(out),   dimension(:,:,:)   :: f_snow_berg
  real, intent(out),   dimension(:,:,:)   ::                 &
                  rain3d, snow3d, snowclr3d, lsc_cloud_area, lsc_liquid,  &
                  lsc_ice, lsc_droplet_number, lsc_ice_number, lsc_snow,  &
                  lsc_rain, lsc_snow_size, lsc_rain_size
! ---> h1g
  real, intent(in) , dimension(:,:,:), optional :: dcond_ls_liquid, dcond_ls_ice
  real, intent(in) , dimension(:,:,:), optional :: Ndrop_act_CLUBB,  Icedrop_act_CLUBB
  real, intent(in) , dimension(:,:,:), optional :: ndust, rbar_dust
  real, intent(in) , dimension(:,:,:), optional :: qcvar_clubb
! <--- h1g

  type(aerosol_type),intent(in), optional :: Aerosol
  logical, intent(in), optional           :: hydrostatic, phys_hydrostatic
  real, intent(in) , dimension(:,:,:),  optional :: mask, cell_cld_frac, meso_cld_frac, &
                                                    shallow_cloud_area

  logical :: used
  integer :: i, j, k, ix, jx, kx, nql, nqi, nqa, nqn, nqg, nqr, nqs, nqni
  real :: qrf, env_fraction, env_qv, dtinv
  real, dimension(size(t,1), size(t,2)) :: ice_lin, graupel_lin
  real, dimension(size(t,1), size(t,2), size(t,3)) :: delp, delz, qsat, &
                                                      convective_humidity_area,     &
                                                      convective_humidity_ratio

      ix=size(t,1) 
      jx=size(t,2) 
      kx=size(t,3)
      dtinv = 1./dt

!----------------------------------------------------------------------
!    define the grid box specific humidity and saturation specific 
!    humidity.
!------------------------------------------------------------------
      call compute_qs (t, pfull, qsat)
 
!----------------------------------------------------------------------
!    define the grid box area whose humidity is affected by the 
!    convective clouds and the environmental fraction and environmental
!    rh.
!-------------------------------------------------------------------
      do k=1, kx
       do j=1, jx
        do i=1, ix
          qrf = MAX (q(i,j,k), 0.0)
          if (do_uw_conv .and. do_donner_deep) then
            convective_humidity_area(i,j,k) = donner_humidity_area(i,j,k) +  &
                           shallow_cloud_area(i,j,k)
            env_qv = qrf - qsat(i,j,k)*(cell_cld_frac(i,j,k) +   &
                           donner_humidity_factor(i,j,k) + shallow_cloud_area(i,j,k))
          else if (do_donner_deep) then
            convective_humidity_area(i,j,k) = donner_humidity_area(i,j,k)
            env_qv = qrf - qsat(i,j,k)*(cell_cld_frac(i,j,k) +   &
                           donner_humidity_factor(i,j,k))
          else if (do_uw_conv) then
            convective_humidity_area(i,j,k) = shallow_cloud_area(i,j,k)
            env_qv = qrf -  shallow_cloud_area(i,j,k)*qsat(i,j,k)
          else
            convective_humidity_area(i,j,k) = 0.0
            env_qv = qrf
          endif
          env_fraction = 1.0 - convective_humidity_area(i,j,k)

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
              convective_humidity_ratio(i,j,k) =    &
                      MAX (qrf*env_fraction/env_qv, 1.0)
 
!---------------------------------------------------------------------
!    grid box is filled with convective clouds.
!----------------------------------------------------------------------
            else
              convective_humidity_ratio(i,j,k) = -10.0
            endif

!--------------------------------------------------------------------
!    either no vapor or all vapor taken up in convective clouds so 
!    none left for large-scale cd.
!---------------------------------------------------------------------
          else
            convective_humidity_ratio(i,j,k) = 1.0
          endif
        end do
       end do
      end do
        
!-----------------------------------------------------------------------
!    call strat_cloud to integrate the prognostic cloud equations. 
!-----------------------------------------------------------------------
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
      if ( .not. do_lin_cld_microphys ) then
        if (do_clubb > 0 .and. do_liq_num) then
         call MG_microp_3D( Time, is, ie, js, je, lon, lat, dt,                                 &
                            pfull, phalf, zhalf, land,                                          &
                            t, q, tracer(:,:,:,nql), tracer(:,:,:,nqi), tracer(:,:,:,nqa),      &
                            tracer(:,:,:,nqn), tracer(:,:,:,nqni), convective_humidity_area,    &
                            dcond_ls_liquid, dcond_ls_ice,                                      &
                            Ndrop_act_CLUBB, Icedrop_act_CLUBB,                                 &
                            ndust, rbar_dust,                                                   &
                            ttnd, qtnd, q_tnd(:,:,:,nql), q_tnd(:,:,:,nqi), q_tnd(:,:,:,nqa),   &
                            q_tnd(:,:,:,nqn), q_tnd(:,:,:,nqni),                                &
                            rain3d, snow3d, rain, snow,                                         &
                            do_clubb=do_clubb, qcvar_clubb = qcvar_clubb,                       &
                            MASK3d=mask,  &
                            lsc_snow = lsc_snow,    &
                            lsc_rain = lsc_rain, &
                            lsc_snow_size = lsc_snow_size,   &
                            lsc_rain_size = lsc_rain_size )
        elseif ( do_legacy_strat_cloud ) then
          if (do_liq_num) then
            call strat_cloud (Time, is, ie, js, je, dt, pfull, phalf,    & 
                            radturbten, t, q, tracer(:,:,:,nql),         &
                            tracer(:,:,:,nqi), tracer(:,:,:,nqa),        &
                            omega, mc_full, diff_t, land, ttnd, qtnd,    &
                            q_tnd(:,:,:,nql), q_tnd(:,:,:,nqi),          &
                            q_tnd(:,:,:,nqa),  f_snow_berg,  &
                            rain3d, snow3d, snowclr3d, &
                            rain, snow, convective_humidity_ratio,   &
                            convective_humidity_area,     &  
                            limit_conv_cloud_frac, mask=mask,            &
                            qn=tracer(:,:,:,nqn), Aerosol=Aerosol,       &
                            SN=q_tnd(:,:,:,nqn))
          else
            call strat_cloud (Time, is, ie, js, je, dt, pfull, phalf,   & 
                          radturbten, t, q, tracer(:,:,:,nql),        &
                          tracer(:,:,:,nqi), tracer(:,:,:,nqa),      &
                          omega, mc_full, diff_t, land, ttnd, qtnd,    &
                          q_tnd(:,:,:,nql), q_tnd(:,:,:,nqi),   &
                          q_tnd(:,:,:,nqa), f_snow_berg,  &
                          rain3d, snow3d, snowclr3d,  &
                          rain, snow, convective_humidity_ratio,  &
                          convective_humidity_area, &
                          limit_conv_cloud_frac, mask=mask)
          endif  ! do_liq_num
        else  ! do_legacy_strat_cloud
          if (do_ice_num) then
            call strat_cloud_new (Time, is, ie, js, je, dt, pfull, phalf, &
                                  zhalf,  zfull, radturbten, t, q,  &
                                  tracer(:,:,:,nql), tracer(:,:,:,nqi), &
                                  tracer(:,:,:,nqa), omega, mc_full, &
                                  diff_t, land, ttnd, qtnd,           &
                                  q_tnd(:,:,:,nql), q_tnd(:,:,:,nqi),    &
                                  q_tnd(:,:,:,nqa),  f_snow_berg, &
                                  rain3d, snow3d,   &
                                  snowclr3d, rain, snow,       &
                                  convective_humidity_ratio,   &
                                  convective_humidity_area,    &
                                  limit_conv_cloud_frac, Aerosol, &       
                                  mask3d=mask, qn_in=tracer(:,:,:,nqn),  &
                                  SN_out = q_tnd(:,:,:,nqn),             &
                                  qni_in=tracer(:,:,:,nqni),   &
                                  SNi_out = q_tnd(:,:,:,nqni),  &
                                  lsc_snow = lsc_snow,    &
                                  lsc_rain = lsc_rain, &
                                  lsc_snow_size = lsc_snow_size,   &
                                  lsc_rain_size = lsc_rain_size )
          else
            if (do_liq_num) then
              call strat_cloud_new (Time, is, ie, js, je, dt, pfull, phalf, &
                                  zhalf, zfull, radturbten, t, q,   &
                                  tracer(:,:,:,nql), tracer(:,:,:,nqi),  &
                                  tracer(:,:,:,nqa), omega, mc_full, &
                                  diff_t, land, ttnd, qtnd,    &
                                  q_tnd(:,:,:,nql), q_tnd(:,:,:,nqi),   &
                                  q_tnd(:,:,:,nqa), f_snow_berg,   &
                                  rain3d, snow3d,   &
                                  snowclr3d,rain, snow,   &
                                  convective_humidity_ratio,   &
                                  convective_humidity_area,&
                                  limit_conv_cloud_frac, Aerosol,   &
                                  mask3d=mask, qn_in=tracer(:,:,:,nqn),  &
                                  SN_out= q_tnd(:,:,:,nqn),   &
                                  lsc_snow = lsc_snow,      &
                                  lsc_rain = lsc_rain, &
                                  lsc_snow_size = lsc_snow_size,   &
                                  lsc_rain_size = lsc_rain_size )
            else
              call strat_cloud_new (Time, is, ie, js, je, dt, pfull, phalf, &
                                  zhalf, zfull, radturbten, t, q,   &
                                  tracer(:,:,:,nql), tracer(:,:,:,nqi),  &
                                  tracer(:,:,:,nqa), omega, mc_full, &
                                  diff_t, land, ttnd, qtnd,    &
                                  q_tnd(:,:,:,nql), q_tnd(:,:,:,nqi),   &
                                  q_tnd(:,:,:,nqa), f_snow_berg,   &
                                  rain3d, snow3d,   &
                                  snowclr3d,rain, snow,   &
                                  convective_humidity_ratio,   &
                                  convective_humidity_area,&
                                  limit_conv_cloud_frac, Aerosol,   &
                                  mask3d=mask,   &
                                  lsc_snow = lsc_snow,      &
                                  lsc_rain = lsc_rain, &
                                  lsc_snow_size = lsc_snow_size,   &
                                  lsc_rain_size = lsc_rain_size )
            end if  ! do_ice_num
          end if  ! do_liq_num
        end if   ! do_legacy_strat_cloud
      else if ( do_lin_cld_microphys ) then
        nqr = get_tracer_index (MODEL_ATMOS, 'rainwat')
        nqs = get_tracer_index (MODEL_ATMOS, 'snowwat')
        nqg = get_tracer_index (MODEL_ATMOS, 'graupel')
        do k=1,kx
          delp(:,:,k) =  phalf(:,:,k+1) - phalf(:,:,k)
          delz(:,:,k) = (zhalf(:,:,k+1)-zhalf(:,:,k))*t(:,:,k)/tm(:,:,k)
        enddo

        call lin_cld_microphys_driver(q, tracer(:,:,:,nql), tracer(:,:,:,nqr), &
                        tracer(:,:,:,nqi), tracer(:,:,:,nqs), tracer(:,:,:,nqg), &
                        tracer(:,:,:,nqa), qtnd, q_tnd(:,:,:,nql),               &
                        q_tnd(:,:,:,nqr), q_tnd(:,:,:,nqi),                      &
                        q_tnd(:,:,:,nqs), q_tnd(:,:,:,nqg), q_tnd(:,:,:,nqa),    &
                        ttnd, t,   pfull, delz, delp, area,                    &
                        dt, land, rain, snow, ice_lin, graupel_lin,              &
                        hydrostatic, phys_hydrostatic,                           &
                        is, ie, js, je, 1, kx, ktop, kx, Time)

! Add all "solid" form of precipitation into surf_snow
        snow = (snow + ice_lin + graupel_lin) * dt/86400.
        rain =  rain * dt/86400.

! Update tendencies:
        rdt(:,:,:,nqr) = rdt(:,:,:,nqr) + q_tnd(:,:,:,nqr)
        rdt(:,:,:,nqs) = rdt(:,:,:,nqs) + q_tnd(:,:,:,nqs)
        rdt(:,:,:,nqg) = rdt(:,:,:,nqg) + q_tnd(:,:,:,nqg)

        ttnd =  ttnd * dt
        qtnd =  qtnd * dt
        q_tnd(:,:,:,nql) = q_tnd(:,:,:,nql) * dt
        q_tnd(:,:,:,nqi) = q_tnd(:,:,:,nqi) * dt
        q_tnd(:,:,:,nqa) = q_tnd(:,:,:,nqa) * dt

! Update rain_wat, snow_wat, graupel_wat
        tracer(:,:,:,nqr) = tracer(:,:,:,nqr) + q_tnd(:,:,:,nqr)*dt
        tracer(:,:,:,nqs) = tracer(:,:,:,nqs) + q_tnd(:,:,:,nqs)*dt
        tracer(:,:,:,nqg) = tracer(:,:,:,nqg) + q_tnd(:,:,:,nqg)*dt
      endif ! not do_lin_cld_microphys
    
!----------------------------------------------------------------------
!    upon return from strat_cloud, update the cloud liquid, ice and area.
!    update the temperature and specific humidity fields.
!----------------------------------------------------------------------
      tracer(:,:,:,nql) = tracer(:,:,:,nql) + q_tnd(:,:,:,nql)              
      tracer(:,:,:,nqi) = tracer(:,:,:,nqi) + q_tnd(:,:,:,nqi)
      tracer(:,:,:,nqa) = tracer(:,:,:,nqa) + q_tnd(:,:,:,nqa)
      if (do_liq_num) tracer(:,:,:,nqn) =    &
                                   tracer(:,:,:,nqn) + q_tnd(:,:,:,nqn)
      if (do_ice_num) tracer(:,:,:,nqni) =    &
                                   tracer(:,:,:,nqni) + q_tnd(:,:,:,nqni)

!   save the lsc fields for use in radiation package.
      lsc_cloud_area(:,:,:) = tracer(:,:,:,nqa)
      lsc_liquid(:,:,:)  =  tracer(:,:,:,nql)
      lsc_ice(:,:,:) =  tracer(:,:,:,nqi)
      if (do_liq_num) lsc_droplet_number(:,:,:) = tracer(:,:,:,nqn)
      if (do_ice_num) lsc_ice_number(:,:,:) = tracer(:,:,:,nqni)

      t = t + ttnd 
      q = q + qtnd
        
      used = send_data (id_qvout, q, Time, is, js, 1, rmask=mask)
      used = send_data (id_qaout, tracer(:,:,:,nqa), Time, is, js, 1, rmask=mask)
      used = send_data (id_qlout, tracer(:,:,:,nql), Time, is, js, 1, rmask=mask)
      used = send_data (id_qiout, tracer(:,:,:,nqi), Time, is, js, 1, rmask=mask)
      if (do_liq_num) then
      used = send_data (id_qnout, tracer(:,:,:,nqn), Time, is, js, 1, rmask=mask)
      endif
      if (do_ice_num) then
      used = send_data (id_qniout, tracer(:,:,:,nqni), Time, is, js, 1, rmask=mask)
      endif

!----------------------------------------------------------------------
!    call strat_cloud_sum to make the cloud variables available for 
!    access by the radiation package. NOTE: this is no longer necessary,
!    and can be judiciously removed (provided other affiliated code 
!    and options are nullified).
!----------------------------------------------------------------------
     call strat_cloud_sum (is, js, tracer(:,:,:,nql),  &
                             tracer(:,:,:,nqi), tracer(:,:,:,nqa))

!----------------------------------------------------------------------
!    convert increments to tendencies.
!----------------------------------------------------------------------
      ttnd = ttnd*dtinv 
      qtnd = qtnd*dtinv
      rain = rain*dtinv 
      snow = snow*dtinv
      q_tnd(:,:,:,nql) = q_tnd(:,:,:,nql)*dtinv
      q_tnd(:,:,:,nqi) = q_tnd(:,:,:,nqi)*dtinv
      q_tnd(:,:,:,nqa) = q_tnd(:,:,:,nqa)*dtinv
      if (do_liq_num) q_tnd(:,:,:,nqn) = q_tnd(:,:,:,nqn)*dtinv
      if (do_ice_num) q_tnd(:,:,:,nqni) = q_tnd(:,:,:,nqni)*dtinv
   
!----------------------------------------------------------------------
!    update the total tendency terms (temperature, vapor specific 
!    humidity, cloud liquid, cloud ice, cloud area, liquid precip,
!    frozen precip) with the contributions from the strat_cloud scheme.
!----------------------------------------------------------------------
      tdt = tdt + ttnd 
      qdt = qdt + qtnd
      rdt(:,:,:,nql) = rdt(:,:,:,nql) + q_tnd(:,:,:,nql)
      rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + q_tnd(:,:,:,nqi)
      rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + q_tnd(:,:,:,nqa)
      if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + q_tnd(:,:,:,nqn)
      if (do_ice_num)         rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + q_tnd(:,:,:,nqni)
      lprec = lprec + rain
      fprec = fprec + snow

end subroutine moistproc_strat_cloud


!#######################################################################
subroutine moistproc_uw_conv(Time, is, ie, js, je, dt, t, q, u, v, tracer,            &
                             pfull, phalf, zfull, zhalf, omega, pblht,        &
                             ustar, bstar, qstar, land, coldT, Aerosol,       &
                             cush, cbmf, cmf, conv_calc_completed,            &
                             available_cf_for_uw, tdt, qdt, udt, vdt, rdt,    &
                             ttnd_conv, qtnd_conv, lprec, fprec, precip,      &
                             liq_precflx, ice_precflx,    &
                             do_strat, do_limit_uw, do_liq_num, num_tracers,  &
                             tracers_in_uw, num_uw_tracers, shallow_cloud_area,&
                             shallow_liquid, shallow_ice,  &
                             shallow_droplet_number, uw_wetdep, &
                             do_ice_num, detrain_ice_num)

  type(time_type), intent(in)   :: Time
  type(aerosol_type),intent(in) :: Aerosol
  integer, intent(in)           :: is, ie,js, je, num_tracers, num_uw_tracers
  real, intent(in)              :: dt
  logical, intent(in)           :: do_strat, do_limit_uw, do_liq_num
  logical, intent(in), dimension(:)       :: tracers_in_uw
  logical, intent(in), dimension(:,:)     :: coldT, conv_calc_completed
  real, intent(in),    dimension(:,:)     :: land, ustar, bstar, qstar, pblht
  real, intent(in),    dimension(:,:,:)   :: pfull, phalf, zfull, zhalf, omega, &
                                             t, q, u, v, available_cf_for_uw
  real, intent(in),    dimension(:,:,:,:) :: tracer
  real, intent(inout), dimension(:,:)     :: lprec, fprec, precip, cush, cbmf
  real, intent(inout), dimension(:,:,:)   :: tdt, qdt, udt, vdt,   &
                                             ttnd_conv, qtnd_conv, cmf
  real, intent(inout), dimension(:,:,:,:) :: rdt
  real, intent(inout), dimension(:,:,:)   :: shallow_cloud_area, &
                                             shallow_liquid,     &
                                             shallow_ice,        &
                                             shallow_droplet_number
  logical, intent(in )                    :: do_ice_num, detrain_ice_num
  real, intent(out),   dimension(:,:,:)   :: liq_precflx, ice_precflx
  real, intent(out),   dimension(:,:,:)   :: uw_wetdep

  integer :: n, nn, nql, nqi, nqa, nqn, nqni
  real, dimension(size(t,1), size(t,2), size(t,3)) :: thlflx, qtflx, precflx
  real, dimension(size(rdt,1), size(rdt,2), size(rdt,3), num_uw_tracers) :: trcr

      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

!----------------------------------------------------------------------
!    if any tracers are to be transported by UW convection, check each
!    active tracer to find those to be transported and fill the 
!    ras_tracers array with these fields.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_uw(n)) then
          trcr(:,:,:,nn) = tracer(:,:,:,n)
          nn = nn + 1
        endif
      end do

      call uw_conv (is, js, Time, t, q, u, v, pfull, phalf, zfull, zhalf, &
                    tracer, omega, dt, pblht, ustar, bstar, qstar, land,  &
                    coldT, Aerosol, cush, do_strat,  conv_calc_completed, &
                    available_cf_for_uw, ttnd_uw(is:ie,js:je,:),          &
                    qtnd_uw(is:ie,js:je,:), qltnd_uw(is:ie,js:je,:),      &
                    qitnd_uw(is:ie,js:je,:), qatnd_uw(is:ie,js:je,:),     &
                    qntnd_uw(is:ie,js:je,:),                              &
                    utnd_uw(is:ie,js:je,:), vtnd_uw(is:ie,js:je,:),       &
                    rain_uw(is:ie,js:je), snow_uw(is:ie,js:je), cmf,      &
                    thlflx, qtflx, precflx, liq_precflx, ice_precflx,     &
                    shallow_liquid, shallow_ice, shallow_cloud_area,      &
                    shallow_droplet_number, cbmf, trcr,                   &
                    qtruw(is:ie,js:je,:,:), uw_wetdep)

!-------------------------------------------------------------------------
!    currently qnitnd_uw is the tendency due to detrainment proportional 
!    by ice mass
!-------------------------------------------------------------------------
      IF ( do_ice_num .AND. detrain_ice_num ) THEN
        CALL detr_ice_num (t, qitnd_uw(is:ie,js:je,:),  &
                                              qnitnd_uw(is:ie,js:je,:) )
      else
        qnitnd_uw(is:ie,js:je,:) = 0.
      END IF

      if (.not. do_limit_uw) then
        tdt=tdt+ttnd_uw(is:ie,js:je,:) 
        qdt=qdt+qtnd_uw(is:ie,js:je,:)
        udt=udt+utnd_uw(is:ie,js:je,:)
        vdt=vdt+vtnd_uw(is:ie,js:je,:)
        ttnd_conv = ttnd_conv + ttnd_uw(is:ie,js:je,:)
        qtnd_conv = qtnd_conv + qtnd_uw(is:ie,js:je,:)
        lprec=lprec+rain_uw(is:ie,js:je)
        fprec=fprec+snow_uw(is:ie,js:je)
        precip=precip+rain_uw(is:ie,js:je)+snow_uw(is:ie,js:je)

        if (do_strat) then
          rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd_uw(is:ie,js:je,:)
          rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd_uw(is:ie,js:je,:)
          rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd_uw(is:ie,js:je,:)
          if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) +   &
                                                 qntnd_uw(is:ie,js:je,:)
          if (do_ice_num) rdt(:,:,:,nqni) = rdt(:,:,:,nqni) +    &
                                                  qnitnd_uw(is:ie,js:je,:)
        endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from uw transport.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_uw(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw(is:ie,js:je,:,nn)
            nn = nn + 1
          endif
        end do
      endif  !(.not. do_limit_uw)

end subroutine moistproc_uw_conv


!#######################################################################
subroutine moistproc_scale_donner(is,ie,js,je,q, delta_temp, delta_q, precip_returned,      &
                                  total_precip, lheat_precip, liquid_precip,    &
                                  frozen_precip, num_tracers, tracers_in_donner,&
                                  qtr, scale)

  integer, intent(in) :: is, ie, js, je, num_tracers
  logical, intent(in), dimension(:)       :: tracers_in_donner
  real, intent(inout), dimension(:,:)     :: precip_returned, total_precip, &
                                             lheat_precip
  real, intent(inout), dimension(:,:,:)   :: q, delta_temp, delta_q,    &
                                             liquid_precip, frozen_precip
  real, intent(inout), dimension(:,:,:,:) :: qtr
  real, intent(out),   dimension(:,:)     :: scale

  integer :: n, nn, i, j, k, ix, jx, kx
  real    :: qvin, dqv
  real, dimension(size(q,1), size(q,2), size(q,3)) :: temp

      ix = size(q,1)
      jx = size(q,2)
      kx = size(q,3)

!     Tendencies coming out of Donner deep are adjusted to prevent
!     the formation of negative water vapor, liquid or ice.

!     (1) Prevent negative liquid and ice specific humidities after
!     tendencies are applied

      where ((qlin(is:ie,js:je,:)+delta_ql(is:ie,js:je,:)) .lt. 0.)
        delta_temp(:,:,:)  = delta_temp (:,:,:) - (qlin(is:ie,js:je,:)+delta_ql(is:ie,js:je,:))*HLV/CP_AIR
        delta_q(:,:,:)     = delta_q    (:,:,:) + (qlin(is:ie,js:je,:)+delta_ql(is:ie,js:je,:))
        delta_ql(is:ie,js:je,:)    = delta_ql   (is:ie,js:je,:) - (qlin(is:ie,js:je,:)+delta_ql(is:ie,js:je,:))
      end where

      where ((qiin(is:ie,js:je,:)+delta_qi(is:ie,js:je,:)) .lt. 0.)
        delta_temp(:,:,:)  = delta_temp (:,:,:) - (qiin(is:ie,js:je,:)+delta_qi(is:ie,js:je,:))*HLS/CP_AIR
        delta_q(:,:,:)     = delta_q    (:,:,:) + (qiin(is:ie,js:je,:)+delta_qi(is:ie,js:je,:))
        delta_qi(is:ie,js:je,:)    = delta_qi   (is:ie,js:je,:) - (qiin(is:ie,js:je,:)+delta_qi(is:ie,js:je,:))
      end where

      where (abs(delta_ql(is:ie,js:je,:) + delta_qi(is:ie,js:je,:)) .lt. 1.e-10 )
        delta_qa(is:ie,js:je,:) = 0.0
      end where

!     (2) Compute limit on Donner tendencies to prevent water vapor
!     from going below 1.e-10. The value of 1.e-10 is consistent with qmin
!      in strat_cloud.F90

!     scaling factor for each grid point
      temp = 1.0
      do k=1,kx
       do j=1,jx
        do i=1,ix
          qvin = q(i,j,k)
          dqv  = delta_q    (i,j,k)
          if ( dqv.lt.0 .and. qvin+dqv.lt.1.e-10 ) then
            temp(i,j,k) = max( 0.0, -(qvin-1.e-10)/dqv )
          endif
        end do
       end do
      end do

!     scaling factor for each column is the minimum value within that column
      scale = minval( temp, dim=3 )

!     scale tendencies
      do k=1,kx
        delta_temp(:,:,k)  = scale(:,:) * delta_temp(:,:,k)
        delta_q(:,:,k)     = scale(:,:) * delta_q    (:,:,k)
        delta_qa(is:ie,js:je,k)    = scale(:,:) * delta_qa(is:ie,js:je,k)
        delta_ql(is:ie,js:je,k)    = scale(:,:) * delta_ql(is:ie,js:je,k)
        delta_qi(is:ie,js:je,k)    = scale(:,:) * delta_qi(is:ie,js:je,k)
      end do

      nn = 1
      do n=1, num_tracers
        if (tracers_in_donner(n)) then
          do k=1,kx
            qtr(:,:,k,nn) = scale(:,:) * qtr(:,:,k,nn)
          end do
          nn = nn + 1
        endif
      end do

      precip_returned = scale*precip_returned

      total_precip = scale*total_precip
      lheat_precip = scale*lheat_precip
      do k=1, kx
        liquid_precip(:,:,k) = scale(:,:)*liquid_precip(:,:,k)
        frozen_precip(:,:,k) = scale(:,:)*frozen_precip(:,:,k)
      end do

end subroutine moistproc_scale_donner


!#######################################################################
subroutine moistproc_scale_uw(is,ie,js,je, dt, q, tracer, tdt, qdt, udt, vdt, rdt,    &
                              ttnd_conv, qtnd_conv, lprec, fprec, precip,&
                              do_strat, do_liq_num, num_tracers,         &
                              tracers_in_uw, scale, do_ice_num)

  integer, intent(in)           :: is, ie, js, je, num_tracers
  real, intent(in)              :: dt
  logical, intent(in)           :: do_strat, do_liq_num, do_ice_num
  logical, intent(in), dimension(:)       :: tracers_in_uw
  real, intent(in),    dimension(:,:,:)   :: q
  real, intent(in),    dimension(:,:,:,:) :: tracer
  real, intent(inout), dimension(:,:)     :: lprec, fprec, precip
  real, intent(inout), dimension(:,:,:)   :: tdt, qdt, udt, vdt,  &
                                             ttnd_conv, qtnd_conv
  real, intent(inout), dimension(:,:,:,:) :: rdt
  real, intent(out),   dimension(:,:)     :: scale

  integer :: n, nn, i, j, k, ix, jx, kx, nql, nqi, nqa, nqn, nqni
  real    :: qvin, dqv
  real, dimension(size(q,1), size(q,2), size(q,3)) :: temp

      ix = size(q,1)
      jx = size(q,2)
      kx = size(q,3)
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

!----------------------------------------------------------------------

!      Tendencies coming out of UW shallow are adjusted to prevent
!      the formation of negative water vapor, liquid or ice.
 
!      (1) Prevent negative liquid and ice specific humidities after tendencies are applied
       temp = tracer(:,:,:,nql)/dt + qltnd_uw(is:ie,js:je,:)
       where (temp(:,:,:) .lt. 0.)
         ttnd_uw(is:ie,js:je,:)  = ttnd_uw(is:ie,js:je,:)  - temp(:,:,:)*HLV/CP_AIR
         qtnd_uw(is:ie,js:je,:)  = qtnd_uw(is:ie,js:je,:)  + temp(:,:,:)
         qltnd_uw(is:ie,js:je,:) = qltnd_uw(is:ie,js:je,:) - temp(:,:,:)
       end where

       temp = tracer(:,:,:,nqi)/dt + qitnd_uw(is:ie,js:je,:)
       where (temp .lt. 0.)
         ttnd_uw(is:ie,js:je,:)  = ttnd_uw(is:ie,js:je,:)  - temp(:,:,:)*HLS/CP_AIR
         qtnd_uw(is:ie,js:je,:)  = qtnd_uw(is:ie,js:je,:)  + temp(:,:,:)
         qitnd_uw(is:ie,js:je,:) = qitnd_uw(is:ie,js:je,:) - temp(:,:,:)
       end where

       where (abs(qltnd_uw(is:ie,js:je,:)+qitnd_uw(is:ie,js:je,:))*dt .lt. 1.e-10 )
         qatnd_uw(is:ie,js:je,:) = 0.0
       end where

!      (2) Compute limit on UW tendencies to prevent water vapor
!      from going below 1.e-10. The value of 1.e-10 is consistent with qmin
!      in strat_cloud.F90

!      scaling factor for each grid point

       temp = 1.0
       do k=1,kx
        do j=1,jx
         do i=1,ix
           qvin = q(i,j,k) + tracer(i,j,k,nql) + tracer(i,j,k,nqi)
           dqv  = ( qtnd_uw(i+is-1,j+js-1,k) + qltnd_uw(i+is-1,j+js-1,k) + qitnd_uw(i+is-1,j+js-1,k) )*dt
           if ( dqv.lt.0 .and. qvin+dqv.lt.1.e-10 ) then
             temp(i,j,k) = max( 0.0, -(qvin-1.e-10)/dqv )
           endif
         end do
        end do
       end do

!      scaling factor for each column is the minimum value within that column
       scale = minval( temp, dim=3 )

!      scale tendencies
       do k=1,kx
         utnd_uw(is:ie,js:je,k)  = scale(:,:) * utnd_uw(is:ie,js:je,k)
         vtnd_uw(is:ie,js:je,k)  = scale(:,:) * vtnd_uw(is:ie,js:je,k)
         ttnd_uw(is:ie,js:je,k)  = scale(:,:) * ttnd_uw(is:ie,js:je,k)
         qtnd_uw(is:ie,js:je,k)  = scale(:,:) * qtnd_uw(is:ie,js:je,k)
         qltnd_uw(is:ie,js:je,k) = scale(:,:) * qltnd_uw(is:ie,js:je,k)
         qitnd_uw(is:ie,js:je,k) = scale(:,:) * qitnd_uw(is:ie,js:je,k)
         qatnd_uw(is:ie,js:je,k) = scale(:,:) * qatnd_uw(is:ie,js:je,k)
       end do

       if (do_liq_num) then
         do k=1,kx
          qntnd_uw(is:ie,js:je,k) = scale(:,:) * qntnd_uw(is:ie,js:je,k)
         end do
       end if

       if (do_ice_num) then
         do k=1,kx
           qnitnd_uw(is:ie,js:je,k) = scale(:,:) * qnitnd_uw(is:ie,js:je,k)
         end do
       end if
       rain_uw(is:ie,js:je) = scale(:,:) * rain_uw(is:ie,js:je)
       snow_uw(is:ie,js:je) = scale(:,:) * snow_uw(is:ie,js:je)

!      update tendencies
       tdt=tdt+ttnd_uw(is:ie,js:je,:)
       qdt=qdt+qtnd_uw(is:ie,js:je,:)
       udt=udt+utnd_uw(is:ie,js:je,:)
       vdt=vdt+vtnd_uw(is:ie,js:je,:)
       ttnd_conv = ttnd_conv + ttnd_uw(is:ie,js:je,:)
       qtnd_conv = qtnd_conv + qtnd_uw(is:ie,js:je,:)

!      update precipitation
       lprec=lprec+rain_uw(is:ie,js:je)
       fprec=fprec+snow_uw(is:ie,js:je)
       precip=precip+rain_uw(is:ie,js:je)+snow_uw(is:ie,js:je)

       if (do_strat) then
         rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd_uw(is:ie,js:je,:)
         rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd_uw(is:ie,js:je,:)
         rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd_uw(is:ie,js:je,:)
         if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) +    &
                                                 qntnd_uw(is:ie,js:je,:)
         if (do_ice_num) rdt(:,:,:,nqni) = rdt(:,:,:,nqni) +  &
                                                 qnitnd_uw(is:ie,js:je,:)
       endif

!------------------------------------------------------------------------
!      update the current tracer tendencies with the contributions 
!      obtained from uw transport.
!------------------------------------------------------------------------
       nn = 1
       do n=1, num_tracers
         if (tracers_in_uw(n)) then
           rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw (is:ie,js:je,:,nn)
           nn = nn + 1
         endif
       end do

!-------------------------------------------------------------------------

end subroutine moistproc_scale_uw



!#########################################################################



end module moistproc_kernels_mod
