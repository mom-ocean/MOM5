#include <fms_platform.h>
MODULE CONV_PLUMES_k_MOD

  use  aer_ccn_act_k_mod,   only: aer_ccn_act_k
  use  conv_utilities_k_mod,only: findt_k, exn_k, qsat_k, adicloud, sounding, uw_params
  use Sat_Vapor_Pres_k_Mod, ONLY: compute_qs_k

!---------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: conv_plumes_k.F90,v 20.0 2013/12/13 23:21:34 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: cp_init_k, cp_end_k, cp_clear_k, ct_init_k, ct_end_k,  &
             ct_clear_k, cumulus_plume_k, cumulus_tend_k

  character(len=11) :: mod_name = 'conv_plumes'

  public cwetdep_type
  type cwetdep_type
   character(len=200) :: scheme
   real :: Henry_constant
   real :: Henry_variable
   real :: frac_in_cloud
   real :: frac_in_cloud_snow
   real :: alpha_r
   real :: alpha_s
   logical :: Lwetdep, Lgas, Laerosol, Lice
  end type cwetdep_type

  public cpnlist
  type cpnlist
     integer :: mixing_assumption, mp_choice
     real :: rle, rpen, rmaxfrac, wmin, rbuoy, rdrag, frac_drs, bigc
     real :: auto_th0, auto_rate, tcrit, cldhgt_max, atopevap, rad_crit,  &
             wtwmin_ratio, deltaqc0, emfrac_max, wrel_min,                        &
             Nl_land, Nl_ocean, r_thresh, qi_thresh, peff_l, peff_i, peff, rh0, cfrac,hcevap, weffect,t00
     logical :: do_ice, do_ppen, do_forcedlifting, do_pevap, do_pdfpcp, isdeep, use_online_aerosol
     logical :: do_auto_aero, do_pmadjt, do_emmax, do_pnqv, do_weffect, do_qctflx_zero,do_detran_zero
     character(len=32), dimension(:), _ALLOCATABLE  :: tracername _NULL
     character(len=32), dimension(:), _ALLOCATABLE  :: tracer_units _NULL
     type(cwetdep_type), dimension(:), _ALLOCATABLE :: wetdep _NULL
  end type cpnlist

  public cplume
  type cplume
     integer :: ltop, let, krel
     real    :: cush, cldhgt, prel, zrel
     real    :: maxcldfrac
     real, _ALLOCATABLE :: thcu  (:) _NULL, qctu  (:) _NULL, uu    (:) _NULL
     real, _ALLOCATABLE :: vu    (:) _NULL, qlu   (:) _NULL, qiu   (:) _NULL
     real, _ALLOCATABLE :: pptr  (:) _NULL, ppti  (:) _NULL, wu    (:) _NULL
     real, _ALLOCATABLE :: umf   (:) _NULL, emf   (:) _NULL, thvu  (:) _NULL
     real, _ALLOCATABLE :: rei   (:) _NULL, fer   (:) _NULL, fdr   (:) _NULL
     real, _ALLOCATABLE :: dp    (:) _NULL, thc   (:) _NULL, qct   (:) _NULL
     real, _ALLOCATABLE :: ql    (:) _NULL, qi    (:) _NULL, qa    (:) _NULL
     real, _ALLOCATABLE :: u     (:) _NULL, v     (:) _NULL, p     (:) _NULL
     real, _ALLOCATABLE :: ps    (:) _NULL, ufrc  (:) _NULL, thvtop(:) _NULL
     real, _ALLOCATABLE :: thvbot(:) _NULL, fdrsat(:) _NULL, z     (:) _NULL
     real, _ALLOCATABLE :: qn    (:) _NULL, qnu   (:) _NULL, zs    (:) _NULL
     real, _ALLOCATABLE :: hlu   (:) _NULL, hl    (:) _NULL, clu   (:) _NULL
     real, _ALLOCATABLE :: ciu   (:) _NULL, buo   (:) _NULL, t     (:) _NULL
     real, _ALLOCATABLE :: crate (:) _NULL, prate (:) _NULL, peff  (:) _NULL
!++++yim
     real, _ALLOCATABLE :: tr  (:,:) _NULL, tru (:,:) _NULL, tru_dwet(:,:) _NULL
     real, _ALLOCATABLE :: pptn  (:) _NULL
  end type cplume

  public ctend
  type ctend
     integer :: botlev, toplev
     real    :: rain, snow, denth, uav, vav, conint, freint,  &
                dtint, dqint, dqtmp, dting
     real, _ALLOCATABLE :: uten  (:) _NULL, vten  (:) _NULL, tten  (:) _NULL
     real, _ALLOCATABLE :: qvten (:) _NULL, qlten (:) _NULL, qiten (:) _NULL
     real, _ALLOCATABLE :: qaten (:) _NULL, thcten(:) _NULL, qctten(:) _NULL
     real, _ALLOCATABLE :: qvdiv (:) _NULL, qldiv (:) _NULL, qidiv (:) _NULL
     real, _ALLOCATABLE :: thcflx(:) _NULL, qctflx(:) _NULL
     real, _ALLOCATABLE :: umflx (:) _NULL, vmflx (:) _NULL, qvflx (:) _NULL
     real, _ALLOCATABLE :: qlflx (:) _NULL, qiflx (:) _NULL, qaflx (:) _NULL
     real, _ALLOCATABLE :: qnflx (:) _NULL, qnten (:) _NULL, pflx  (:) _NULL
     real, _ALLOCATABLE :: hlflx (:) _NULL, hlten (:) _NULL
     real, _ALLOCATABLE :: tevap (:) _NULL, qevap (:) _NULL
     real, _ALLOCATABLE :: qldet (:) _NULL, qidet (:) _NULL, qadet (:) _NULL
     real, _ALLOCATABLE :: qndet (:) _NULL
!++++yim
     real, _ALLOCATABLE :: trflx(:,:) _NULL,trten (:,:) _NULL, trwet(:,:) _NULL
  end type ctend

contains

!#####################################################################
!#####################################################################

!++++yim
  subroutine cp_init_k(kd, num_tracers, cp)
!++++yim    
    integer, intent(in) :: kd, num_tracers
    type(cplume), intent(inout) :: cp
    
    allocate ( cp%hlu   (0:kd)); cp%hlu   =0.;
    allocate ( cp%thcu  (0:kd)); cp%thcu  =0.;
    allocate ( cp%qctu  (0:kd)); cp%qctu  =0.;
    allocate ( cp%uu    (0:kd)); cp%uu    =0.;
    allocate ( cp%vu    (0:kd)); cp%vu    =0.;
    allocate ( cp%qlu   (0:kd)); cp%qlu   =0.;
    allocate ( cp%qiu   (0:kd)); cp%qiu   =0.;
    allocate ( cp%clu   (0:kd)); cp%clu   =0.;
    allocate ( cp%ciu   (0:kd)); cp%ciu   =0.;
    allocate ( cp%buo   (0:kd)); cp%buo   =0.;
    allocate ( cp%t     (0:kd)); cp%t     =0.;
    allocate ( cp%crate (0:kd)); cp%crate =0.;
    allocate ( cp%prate (0:kd)); cp%prate =0.;
    allocate ( cp%qnu   (0:kd)); cp%qnu   =0.;
    allocate ( cp%pptn  (1:kd)); cp%pptn  =0.;
    allocate ( cp%pptr  (1:kd)); cp%pptr  =0.;
    allocate ( cp%ppti  (1:kd)); cp%ppti  =0.;
    allocate ( cp%wu    (0:kd)); cp%wu    =0.;
    allocate ( cp%umf   (0:kd)); cp%umf   =0.;
    allocate ( cp%emf   (0:kd)); cp%emf   =0.;
    allocate ( cp%thvu  (1:kd)); cp%thvu  =0.;
    allocate ( cp%rei   (1:kd)); cp%rei   =0.;
    allocate ( cp%fer   (1:kd)); cp%fer   =0.;
    allocate ( cp%fdr   (1:kd)); cp%fdr   =0.;
    allocate ( cp%dp    (1:kd)); cp%dp    =0.;
    allocate ( cp%hl    (1:kd)); cp%hl    =0.;
    allocate ( cp%thc   (1:kd)); cp%thc   =0.;
    allocate ( cp%qct   (1:kd)); cp%qct   =0.;
    allocate ( cp%u     (1:kd)); cp%u     =0.;
    allocate ( cp%v     (1:kd)); cp%v     =0.;
    allocate ( cp%ql    (1:kd)); cp%ql    =0.;
    allocate ( cp%qi    (1:kd)); cp%qi    =0.;
    allocate ( cp%qa    (1:kd)); cp%qa    =0.;
    allocate ( cp%qn    (1:kd)); cp%qn    =0.;
    allocate ( cp%p     (1:kd)); cp%p     =0.;
    allocate ( cp%ps    (0:kd)); cp%ps    =0.;
    allocate ( cp%z     (1:kd)); cp%z     =0.;
    allocate ( cp%zs    (0:kd)); cp%zs    =0.;
    allocate ( cp%ufrc  (1:kd)); cp%ufrc  =0.;
    allocate ( cp%thvbot(1:kd)); cp%thvbot=0.;
    allocate ( cp%thvtop(1:kd)); cp%thvtop=0.;
    allocate ( cp%fdrsat(1:kd)); cp%fdrsat=0.;
    allocate ( cp%peff  (1:kd)); cp%peff  =0.;
!++++yim
    allocate ( cp%tr(1:kd,1:num_tracers)); cp%tr=0.;
    allocate ( cp%tru(0:kd,1:num_tracers)); cp%tru=0.;
    allocate ( cp%tru_dwet(0:kd,1:num_tracers)); cp%tru_dwet=0.;
  end subroutine cp_init_k

!#####################################################################
!#####################################################################

  subroutine cp_end_k (cp)
    type(cplume), intent(inout) :: cp
    deallocate (cp%thcu, cp%qctu, cp%uu, cp%vu, cp%qlu, cp%qiu,  &
                cp%pptr, cp%ppti, cp%wu, cp%umf, cp%emf, cp%thvu,  &
                cp%rei, cp%fer, cp%fdr, cp%dp, cp%thc, cp%qct, cp%u, &
                cp%v, cp%p, cp%ps, cp%ufrc, cp%thvbot, cp%thvtop, &
                cp%fdrsat, cp%peff, cp%qnu, cp%ql, cp%qi, cp%qa, cp%qn, &
                cp%pptn, cp%z, cp%zs, cp%hl, cp%hlu, cp%clu, cp%ciu, &
                cp%buo, cp%t, cp%crate, cp%prate, cp%tr, cp%tru, &
                cp%tru_dwet)
  end subroutine cp_end_k

!#####################################################################
!#####################################################################

  subroutine cp_clear_k (cp)
    type(cplume), intent(inout) :: cp
    cp%thcu  =0.;    cp%qctu  =0.;    cp%uu    =0.;    cp%vu    =0.;
    cp%qlu   =0.;    cp%qiu   =0.;    cp%qnu   =0.;    cp%pptr  =0.;
    cp%ppti  =0.;    cp%wu    =0.;    cp%umf   =0.;    cp%emf   =0.;
    cp%thvu  =0.;    cp%rei   =0.;    cp%fer   =0.;    cp%fdr   =0.;
    cp%dp    =0.;    cp%thc   =0.;    cp%qct   =0.;    cp%u     =0.;
    cp%v     =0.;    cp%ql    =0.;    cp%qi    =0.;    cp%qa    =0.; 
    cp%qn    =0.;    cp%p     =0.;    cp%ps    =0.;
    cp%ufrc  =0.;    cp%thvbot=0.;    cp%thvtop=0.;    cp%hlu   =0.;
    cp%fdrsat=0.;    cp%z     =0.;    cp%zs    =0.;    cp%hl    =0.;
    cp%clu   =0.;    cp%ciu   =0.;    cp%buo   =0.;    cp%t     =0.;
    cp%crate =0.;    cp%prate =0.;    cp%peff  =0.;    !cp%maxcldfrac = 1.;
!++++yim
    cp%pptn  =0.;    cp%tr    =0.;    cp%tru   =0.;    cp%tru_dwet = 0.
  end subroutine cp_clear_k

!#####################################################################
!#####################################################################

!++++yim
  subroutine ct_init_k (kd, num_tracers, ct)
!++++yim
    integer, intent(in) :: kd, num_tracers
    type(ctend), intent(inout) :: ct
    allocate ( ct%uten  (1:kd)); ct%uten  =0.;
    allocate ( ct%vten  (1:kd)); ct%vten  =0.;
    allocate ( ct%tten  (1:kd)); ct%tten  =0.;
    allocate ( ct%qvten (1:kd)); ct%qvten =0.;
    allocate ( ct%qlten (1:kd)); ct%qlten =0.;
    allocate ( ct%qiten (1:kd)); ct%qiten =0.;
    allocate ( ct%qaten (1:kd)); ct%qaten =0.;
    allocate ( ct%qnten (1:kd)); ct%qnten =0.;
    allocate ( ct%qldet (1:kd)); ct%qldet =0.;
    allocate ( ct%qidet (1:kd)); ct%qidet =0.;
    allocate ( ct%qadet (1:kd)); ct%qadet =0.;
    allocate ( ct%qndet (1:kd)); ct%qndet =0.;
    allocate ( ct%hlten (1:kd)); ct%hlten =0.;
    allocate ( ct%thcten(1:kd)); ct%thcten=0.;
    allocate ( ct%qctten(1:kd)); ct%qctten=0.;
    allocate ( ct%qvdiv (1:kd)); ct%qvdiv =0.;
    allocate ( ct%qldiv (1:kd)); ct%qldiv =0.;
    allocate ( ct%qidiv (1:kd)); ct%qidiv =0.;
    allocate ( ct%hlflx (0:kd)); ct%hlflx =0.;
    allocate ( ct%thcflx(0:kd)); ct%thcflx=0.;
    allocate ( ct%qctflx(0:kd)); ct%qctflx=0.;
    allocate ( ct%qvflx (0:kd)); ct%qvflx =0.;
    allocate ( ct%qlflx (0:kd)); ct%qlflx =0.;
    allocate ( ct%qiflx (0:kd)); ct%qiflx =0.;
    allocate ( ct%qaflx (0:kd)); ct%qaflx =0.;
    allocate ( ct%qnflx (0:kd)); ct%qnflx =0.;
    allocate ( ct%umflx (0:kd)); ct%umflx =0.;
    allocate ( ct%vmflx (0:kd)); ct%vmflx =0.;
    allocate ( ct%pflx  (1:kd)); ct%pflx  =0.;
    allocate ( ct%tevap (1:kd)); ct%tevap =0.;
    allocate ( ct%qevap (1:kd)); ct%qevap =0.;
!++++yim
    allocate ( ct%trflx  (0:kd,1:num_tracers)); ct%trflx  =0.;
    allocate ( ct%trten  (1:kd,1:num_tracers)); ct%trten  =0.;
    allocate ( ct%trwet  (1:kd,1:num_tracers)); ct%trwet  =0.;
  end subroutine ct_init_k

!#####################################################################
!#####################################################################

  subroutine ct_end_k (ct)
    type(ctend), intent(inout) :: ct
    deallocate (ct%uten, ct%vten, ct%tten, ct%qvten, ct%qlten,  &
                ct%qiten, ct%qaten, ct%qnten, ct%hlten, ct%thcten, &
                ct%qldet, ct%qidet, ct%qadet, ct%qndet,            &
                ct%qctten, ct%qvdiv, ct%qldiv, ct%qidiv, ct%hlflx, &
                ct%thcflx, ct%qctflx, ct%qvflx, ct%qlflx, ct%qiflx, &
                ct%qaflx, ct%qnflx, ct%umflx, ct%vmflx, ct%pflx, &
                ct%tevap, ct%qevap, ct%trflx, ct%trten, ct%trwet)
  end subroutine ct_end_k

!#####################################################################
!#####################################################################

  subroutine ct_clear_k(ct)
    type(ctend), intent(inout) :: ct
    
    ct%uten  =0.;    ct%vten  =0.;    ct%tten  =0.;
    ct%qvten =0.;    ct%qlten =0.;    ct%qiten =0.;
    ct%qaten =0.;    ct%qnten =0.;    ct%thcten=0.;
    ct%qctten=0.;    
    ct%qldet =0.;    ct%qidet =0.;    ct%qadet =0.;
    ct%qndet =0.;    
    ct%qvdiv =0.;    ct%qldiv =0.;    ct%qidiv =0.;
    ct%thcflx=0.;    ct%qctflx=0.;    ct%qvflx =0.;
    ct%qlflx =0.;    ct%qiflx =0.;    ct%qaflx =0.;    ct%qnflx =0.;
    ct%umflx =0.;    ct%vmflx =0.;    ct%pflx  =0.;
    ct%hlflx =0.;    ct%hlten =0.;
    ct%denth =0.;    ct%rain  =0.;    ct%snow  =0.;    
    ct%tevap =0.;    ct%qevap =0.;
    ct%dting =0.;
!++++yim
    ct%trflx =0.;    ct%trten =0.;    ct%trwet =0.
  end subroutine ct_clear_k

!#####################################################################
!#####################################################################

  subroutine mixing_k (cpn, z0, p0, hl0, thc0, qct0, hlu, thcu, qctu, &
                       wu, scaleh, rei, fer, fdr, fdrsat, rho0j, rkm, &
                       Uw_p, umfkm1, dp, dt)       
  
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: z0, p0, hl0, thc0, qct0 !envirn. properties at level k
    real,           intent(in)    :: hlu, thcu, qctu, wu !updraft properties at level k-1
    real,           intent(in)    :: scaleh, rkm
    real,           intent(in)    :: umfkm1, dp, dt
    real,           intent(inout) :: rei, fer, fdr, fdrsat, rho0j
    
    real    :: excessu, excess0, hlfs, qtfs, thvfs,  &
               xbuo0, xsat, xs, xs1, xs2
    real    :: thj, qvj, qlj, qij, qse, thvj, thv0j
    real    :: aquad, bquad, cquad, ee2, ud2
    real    :: emmax

!-----A.  Entrainment and Detrainment
!     first, to determine fraction (xsat) of mixture that is to be detrained out 
!     of clouds, i.e., the mixture with negative buoyancy. We consider a thin 
!     layer between two interfaces, so using mid-point value to represent the 
!     mean value of the layer. The properties of updraft at midpoint is assumed
!     to be undiluted from the lower interface.  

!-----calculate fraction of mixture that is just saturated

    excessu = qctu - qsat_k((hlu-Uw_p%grav*z0)/Uw_p%cp_air, p0,Uw_p) 
              excessu = max(excessu,0.0)
    excess0 = qct0 - qsat_k((hl0-Uw_p%grav*z0)/Uw_p%cp_air, p0, Uw_p) 

    if(excessu*excess0.le.0)then
       xsat = -excessu/(excess0-excessu)
    else
       xsat = 1.0
    endif
    hlfs =(1.-xsat)*hlu  + xsat*hl0
    qtfs =(1.-xsat)*qctu + xsat*qct0
    call findt_k (z0,p0,hlfs,qtfs, thj, qvj, qlj, qij, qse, thvfs, &
                  cpn%do_ice, Uw_p)
    call findt_k (z0,p0,hlu, qctu, thj, qvj, qlj, qij, qse, thvj, &
                  cpn%do_ice, Uw_p)   
    call findt_k (z0,p0,hl0, qct0, thj, qvj, qlj, qij, qse, thv0j, &
                  cpn%do_ice, Uw_p)   
    rho0j = p0/(Uw_p%rdgas*thv0j*exn_k(p0,Uw_p))

!-----calculate fraction of mixture with zero buoyancy
    if(thvfs.ge.thv0j) then
       xbuo0=xsat
    else if(thvj.le.thv0j) then
       xbuo0=0.
    else
       xbuo0=xsat*(thvj-thv0j)/(thvj-thvfs)
    endif

    !-----calculate fraction of mixture with negative buoyancy but can 
    !     penetrate a critical distance lc=rle*scaleh
    if(thvfs.ge.thv0j.or.xsat.le.0.05) then
       xs=xsat !mixture has to be saturated
    else
       aquad = wu**2.
       bquad = -(2.*wu**2. + 2.*cpn%rbuoy*Uw_p%grav*cpn%rle*scaleh*&
                                              (thvj-thvfs)/thv0j/xsat)
       cquad = wu**2. - 2.*cpn%rbuoy*Uw_p%grav*cpn%rle*scaleh* &
                                                        (1-thvj/thv0j)
       call roots(aquad,bquad,cquad,xs1,xs2)
       xs=min(xs1,xs2)
    endif
    xs=min(xs,xsat)
    xs=max(xbuo0,xs)
    xs=min(1.0,xs)
    
    ee2     = xs**2.
    ud2     = 1. - 2.*xs + xs**2.
    rei     = rkm/scaleh/Uw_p%grav/rho0j  !make entrainment rate in unit of 1/Pa
    fer     = rei * ee2
    fdr     = rei * ud2

 if(cpn%do_emmax) then
    emmax = cpn%emfrac_max * dp / dt / Uw_p%GRAV
    if ((fer-fdr)*dp*umfkm1 .gt. emmax) then
       rei = emmax / dp / umfkm1 / (ee2-ud2)
       fer = rei * ee2
       fdr = rei * ud2
    end if
 end if

    fdrsat  = rei * (ud2-(1. - 2.*xsat + xsat**2.))

    if (fdr.ne.0) then
       fdrsat  = min(fdrsat/fdr, 1.)
    else
       fdrsat  = 0.
    end if
    fdrsat  = max(fdrsat, cpn%frac_drs)

  end subroutine mixing_k

!#####################################################################
!#####################################################################

  subroutine cumulus_plume_k (cpn, sd, ac, cp, rkm, cbmf, wrel, scaleh,&
                              Uw_p, ier, ermesg)       
  

    type(cpnlist),      intent(in)    :: cpn
    type(uw_params),    intent(inout) :: Uw_p
    type(sounding),     intent(in)    :: sd
    type(adicloud),     intent(in)    :: ac
    real,               intent(in)    :: rkm, cbmf, wrel, scaleh
    type(cplume),       intent(inout) :: cp
    integer,            intent(out)   :: ier
    character(len=*),   intent(out)   :: ermesg

    real, dimension(4)            :: totalmass
    integer                       :: tym
    real                          :: drop
   
    integer :: k, klm, km1, krel, let, ltop
    real    :: thv0rel, wtw, wtwtop
    real    :: thj, qvj, qlj, qij, qse, rhos0j, rho0j
    real    :: bogtop, bogbot, delbog, drage, expfac
    real    :: zrel, prel, nu, leff, qrj, qsj, temp
    real    :: qctu_new, hlu_new, qlu_new, qiu_new, clu_new, ciu_new
    real    :: scaleh1
    real    :: qct_env_k, hl_env_k
    real    :: t_mid, tv_mid, air_density, total_condensate,   &
               total_rain, total_snow, delta_tracer, delta_qn, wrel2, gamma
    real    :: cflim
    integer :: n
    logical :: kbelowlet

    ier = 0
    ermesg = ' '
    tym = size(totalmass,1)
    call cp_clear_k (cp)
    cp%p=sd%p; cp%ps=sd%ps; cp%dp=sd%dp; cp%u=sd%u; cp%v=sd%v;
    cp%hl=sd%hl; cp%thc=sd%thc; cp%qct=sd%qct; 
    cp%ql=sd%ql; cp%qi=sd%qi; cp%qa=sd%qa; cp%qn=sd%qn;
!++++yim
    cp%tr=sd%tr;
    cp%thvbot=sd%thvbot; cp%thvtop=sd%thvtop;
    cp%z=sd%z; cp%zs=sd%zs;

    wtw  = wrel*wrel
    wtwtop =  cpn%wtwmin_ratio * wtw

    delta_qn =0.

    !determine release height and parcel properties (krel, prel, thv0rel, thvurel)
    if(ac % plcl .gt. sd % pinv)then
       krel    = sd % kinv
       prel    = sd % pinv
       zrel    = sd % zinv
       thv0rel = sd % thvinv
    else
       krel    = ac % klcl
       prel    = ac % plcl
       zrel    = ac % zlcl
       thv0rel = ac % thv0lcl
    endif

    cp%krel=krel
    cp%prel=prel
    cp%zrel=zrel

    !(krel-1) represents the bottom of the updraft
    call findt_k (zrel,prel,ac%hlsrc,ac%qctsrc, thj, qvj, qlj,  &
                  qij, qse, cp%thvu(krel-1), cpn%do_ice, Uw_p)   
    cp%ps   (krel-1) = max(prel, cp%ps(krel)+1.) !prel
    cp%hlu  (krel-1) = ac % hlsrc
    cp%thcu (krel-1) = ac % thcsrc
    cp%qctu (krel-1) = ac % qctsrc
    cp%uu   (krel-1) = ac % usrc
    cp%vu   (krel-1) = ac % vsrc
    cp%umf  (krel-1) = cbmf*sd%rho(krel-1)/ac%rho0lcl
    cp%wu   (krel-1) = wrel
    cp%ufrc (krel-1) = cp%umf(krel-1)/(sd%rho(krel-1)*cp%wu(krel-1))
    cp%tru  (krel-1,:) = cp%tr(1,:)
    !==================================================
    !     yim's CONVECTIVE NUCLEATION
    !==================================================
    totalmass(1)=     sd%am1(krel-1); !totalmass(1)=aerol;
    totalmass(2)=     sd%am2(krel-1); !totalmass(2)=0.;
    totalmass(3)=     sd%am3(krel-1); !totalmass(3)=0.;
    totalmass(4)=     sd%am4(krel-1); !totalmass(4)=0.;
    if (SUM(totalmass(:)) /= 0.0 .and. cpn%use_online_aerosol) then
      wrel2 = wrel*cpn%wrel_min
      call aer_ccn_act_k(thj*exn_k(prel,Uw_p), prel, wrel2, totalmass, &
                         tym, drop, ier, ermesg)
      if (ier /= 0) then
        return
      endif
      cp%qnu(krel-1) = drop * 1.0e6 / (prel /    &
                      (Uw_p%rdgas*cp%thvu(krel-1)*exn_k(prel,Uw_p)))
    else
      drop = 0.
      cp%qnu(krel-1) = 0.0
    endif


    !(krel) represents the first partial updraft layer
    cp%z      (krel) = (cp%zs(krel-1) + cp%zs(krel))*0.5
    cp%p      (krel) = (cp%ps(krel-1) + cp%ps(krel))*0.5
    cp%dp     (krel) =  cp%ps(krel-1) - cp%ps(krel)   
    cp%thvbot (krel) = thv0rel
    if(krel.ne. sd % kinv) then
       cp%hl  (krel) = cp%hl (krel)+sd%sshl (krel)*  &
                                                (cp%p(krel)-sd%p(krel))
       cp%thc (krel) = cp%thc(krel)+sd%ssthc(krel)*   &
                                                (cp%p(krel)-sd%p(krel))
       cp%qct (krel) = cp%qct(krel)+sd%ssqct(krel)* &
                                               (cp%p(krel)-sd%p(krel))
       call findt_k (cp%z(krel),cp%p(krel),cp%hl(krel), cp%qct(krel), &
                     thj, qvj, qlj, qij, qse, cp%thvu(krel),  &
                     cpn%do_ice, Uw_p)
    endif

    !Compute updraft properties above the LCL
    kbelowlet=.true.
    let=krel
    klm=sd%ktopconv-1
    do k=krel,klm
       km1=k-1
       hl_env_k  = cp%hl(k)
       qct_env_k = cp%qct(k)

       !Calculation entrainment and detrainment rate
       if (cpn%mixing_assumption.eq.0) then
          scaleh1 = scaleh
          call mixing_k (cpn, cp%z(k), cp%p(k), hl_env_k, cp%thc(k), &
                         qct_env_k, cp%hlu(km1), cp%thcu(km1),  &
                         cp%qctu(km1), cp%wu(km1), scaleh1, cp%rei(k), &
                         cp%fer(k), cp%fdr(k), cp%fdrsat(k), rho0j, &
                         rkm, Uw_p, cp%umf(km1), cp%dp(k), sd%delt)      
       else if (cpn%mixing_assumption.eq.1) then
          temp         = sqrt(cp%ufrc(km1)) !scaleh for fixed length scale for donner_plumes
          rho0j        = sd%rho(k)
          cp%rei(k)    = rkm/temp/Uw_p%grav/rho0j
          cp%fer(k)    = cp%rei(k)
          cp%fdr(k)    = 0.
          cp%fdrsat(k) = 0.
       else if (cpn%mixing_assumption.eq.2) then
          scaleh1 = max(1000., cp%z(k))
          call mixing_k (cpn, cp%z(k), cp%p(k), hl_env_k, cp%thc(k), &
                         qct_env_k, cp%hlu(km1), cp%thcu(km1),  &
                         cp%qctu(km1), cp%wu(km1), scaleh1, cp%rei(k), &
                         cp%fer(k), cp%fdr(k), cp%fdrsat(k), rho0j, &
                         rkm, Uw_p, cp%umf(km1), cp%dp(k), sd%delt)      
       else if (cpn%mixing_assumption.eq.3) then
          scaleh1 = cpn%t00*Uw_p%rdgas/Uw_p%grav*log(cp%p(1)/cp%p(k))
          scaleh1 = max (1000., scaleh1)
          call mixing_k (cpn, cp%z(k), cp%p(k), hl_env_k, cp%thc(k), &
                         qct_env_k, cp%hlu(km1), cp%thcu(km1),      &
                         cp%qctu(km1), cp%wu(km1), scaleh1, cp%rei(k), &
                         cp%fer(k), cp%fdr(k), cp%fdrsat(k), rho0j, &
                         rkm, Uw_p, cp%umf(km1), cp%dp(k), sd%delt)
       else if (cpn%mixing_assumption.eq.4) then
          scaleh1 = 2000.
          call mixing_k (cpn, cp%z(k), cp%p(k), hl_env_k, cp%thc(k), &
                         qct_env_k, cp%hlu(km1), cp%thcu(km1),  &
                         cp%qctu(km1), cp%wu(km1), scaleh1, cp%rei(k), &
                         cp%fer(k), cp%fdr(k), cp%fdrsat(k), rho0j, &
                         rkm, Uw_p, cp%umf(km1), cp%dp(k), sd%delt)
       else
          scaleh1 = cpn%t00*Uw_p%rdgas/Uw_p%grav*log(cp%p(1)/cp%p(k))
          scaleh1 = max (1000., scaleh1)
          call mixing_k (cpn, cp%z(k), cp%p(k), hl_env_k, cp%thc(k), &
                         qct_env_k, cp%hlu(km1), cp%thcu(km1),      &
                         cp%qctu(km1), cp%wu(km1), scaleh1, cp%rei(k), &
                         cp%fer(k), cp%fdr(k), cp%fdrsat(k), rho0j, &
                         rkm, Uw_p, cp%umf(km1), cp%dp(k), sd%delt)
       end if

       !Calculate the mass flux
       cp%umf(k)=cp%umf(km1)*exp(cp%dp(k)*(cp%fer(k)-cp%fdr(k)))
       cp%emf(k)=0.0

       !Thermodynamics for the dilute plume
       cp%hlu (k)=hl_env_k -(hl_env_k -cp%hlu (km1))*  &
                                            exp(-cp%fer(k)*cp%dp(k))
       cp%qctu(k)=qct_env_k-(qct_env_k-cp%qctu(km1))*  &
                                            exp(-cp%fer(k)*cp%dp(k))
       cp%qnu (k)=cp%qn (k)-(cp%qn (k)-cp%qnu (km1))*  &
                                            exp(-cp%fer(k)*cp%dp(k))
       cp%tru (k,:)=cp%tr (k,:)-(cp%tr (k,:)-cp%tru (km1,:))*  &
                                            exp(-cp%fer(k)*cp%dp(k))
       if(cp%fer(k)*cp%dp(k).lt.1.e-4)then
          cp%uu(k)=cp%uu(km1) - sd%dudp(k)*cp%dp(k)
          cp%vu(k)=cp%vu(km1) - sd%dvdp(k)*cp%dp(k)
       else
          cp%uu(k)=cp%u(k)-cpn%bigc*sd%dudp(k)/cp%fer(k)-  &
                   exp(-cp%fer(k)*cp%dp(k))* &
                  (cp%u(k)-cpn%bigc*sd%dudp(k)/cp%fer(k) - cp%uu(km1))
          cp%vu(k)=cp%v(k)-cpn%bigc*sd%dvdp(k)/cp%fer(k)-  &
                   exp(-cp%fer(k)*cp%dp(k))* &
                  (cp%v(k)-cpn%bigc*sd%dvdp(k)/cp%fer(k) - cp%vu(km1))
       endif
       if (cpn%mp_choice.eq.0) then
          call micro_donner_k (cpn, cp%zs(k), cp%ps(k), cp%hlu(k), &
                               cp%qctu(k), cp%zs(km1), cp%qlu(km1), &
                               cp%clu(km1), cp%qiu(km1), cp%ciu(km1), &
                               cp%wu(km1), cp%crate(k), cp%prate(k),   &
                               qrj, qsj, qlu_new, clu_new, qiu_new, &
                               ciu_new, hlu_new, qctu_new, temp, &
                               cpn%do_ice, Uw_p)
       else if (cpn%mp_choice.eq.1) then
          call precipitation_k (cp%zs(k), cp%ps(k), cp%hlu(k), &
                                cp%qctu(k), cp%qnu(k), cpn, qrj, qsj, &
                                hlu_new, qctu_new, qlu_new, qiu_new, &
                                clu_new, ciu_new, temp, cpn%do_ice, &
                                delta_qn, Uw_p, kbelowlet)        
       else if (cpn%mp_choice.eq.2) then
          call precip_new_k    (cp%zs(k), cp%ps(k), cp%hlu(k), &
                                cp%qctu(k), cp%qnu(k), cpn, qrj, qsj, &
                                hlu_new, qctu_new, qlu_new, qiu_new, &
                                clu_new, ciu_new, temp, cpn%do_ice, &
                                delta_qn, Uw_p, kbelowlet, sd%land, sd%delt, cp%wu(km1))
       else if (cpn%mp_choice.eq.3) then
          call precip3_k       (cp%zs(k), cp%ps(k), cp%hlu(k), &
                                cp%qctu(k), cp%qnu(k), cpn, qrj, qsj, &
                                hlu_new, qctu_new, qlu_new, qiu_new, &
                                clu_new, ciu_new, temp, cpn%do_ice, &
                                delta_qn, Uw_p, kbelowlet)        
       else if (cpn%mp_choice.eq.4) then
          call precip4_k (cp%zs(k), cp%ps(k), cp%hlu(k), &
                                cp%qctu(k), cp%qnu(k), cpn, qrj, qsj, &
                                hlu_new, qctu_new, qlu_new, qiu_new, &
                                clu_new, ciu_new, temp, cpn%do_ice, &
                                delta_qn, Uw_p, kbelowlet)        
       else if (cpn%mp_choice.eq.5) then
          call precip5_k (cp%zs(k), cp%ps(k), cp%hlu(k), &
                                cp%qctu(k), cp%qnu(k), cpn, qrj, qsj, &
                                hlu_new, qctu_new, qlu_new, qiu_new, &
                                clu_new, ciu_new, temp, cpn%do_ice, &
                                delta_qn, Uw_p, kbelowlet, sd%dp(k))        
       end if

       cp%qctu(k)=qctu_new
       cp%hlu (k)=hlu_new
       cp%qlu (k)=qlu_new
       cp%qiu (k)=qiu_new
       cp%clu (k)=clu_new
       cp%ciu (k)=ciu_new
       cp%peff(k)=(qrj+qsj)/max(qlu_new+qiu_new+qrj+qsj,1.e-28);

       cp%thvu(k)=temp/exn_k(cp%ps(k),Uw_p)*(1.+Uw_p%zvir*(cp%qctu(k)-cp%qlu(k)-cp%qiu(k))-cp%qlu(k)-cp%qiu(k))
       cp%buo (k)=(cp%thvu(k)-cp%thvtop(k))!/cp%thvtop(k)*Uw_p%grav
       cp%t   (k)=temp
       nu = max(min((268. - temp)/20.,1.0),0.0)
       leff = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
       cp%thcu(k)=temp/exn_k(cp%ps(k),Uw_p)

       !Calculate vertical velocity
!!$       bogbot = (cp%thvu(km1)/cp%thvbot(k) - 1.)
!!$       if(bogbot.gt.0.)then
!!$          bogbot =  bogbot /cpn%rbuoy
!!$       else
!!$          bogbot =  bogbot*cpn%rbuoy
!!$       endif
!!$       bogtop = (cp%thvu(k)/cp%thvtop(k) - 1.)
!!$       if(bogtop.gt.0.)then
!!$          bogtop =  bogtop /cpn%rbuoy
!!$       else
!!$          bogtop =  bogtop*cpn%rbuoy
!!$       endif

       bogbot = (cp%thvu(km1)/cp%thvbot(k) - 1.)
       bogbot =  bogbot*cpn%rbuoy
       bogtop = (cp%thvu(k)/cp%thvtop(k) - 1.)
       bogtop =  bogtop*cpn%rbuoy

       if(bogbot.gt.0.and.bogtop.gt.0) then
          let = k
          kbelowlet = .true.
       else
          kbelowlet = .false.
       end if

       delbog = bogtop - bogbot
       drage = cp%fer(k) * ( 1. + cpn%rdrag )
       expfac = exp(-2.* drage * cp%dp(k))
       if(drage * cp%dp(k).gt.1.e-3) then
          wtw = wtw*expfac + (delbog + (1.-expfac) *               &
               (bogbot+delbog/(-2.*drage*cp%dp(k))))/(rho0j*drage)
       else
          wtw = wtw + cp%dp(k) * (bogbot+bogtop)/rho0j
       endif
       wtwtop = max( cpn%wtwmin_ratio * wtw, wtwtop )

       if (cpn%do_forcedlifting) then
          if (wtw.le.0. .and. k <= ac%klfc) then
             wtw=wrel*wrel
!         else if (wtw.le.0.) then
           else if (wtw.le.wtwtop) then
             exit
          end if
       else
!         if(wtw.le.0.) exit
          if(wtw.le.wtwtop) exit
       end if

       cp%wu(k) = sqrt(wtw)
       if(cp%wu(k).gt.100.)then
          print *, 'Very big wu in UW-ShCu',bogbot,bogtop,expfac,cp%fer(k)
          return
       endif
       
       rhos0j     = cp%ps(k)/(Uw_p%rdgas*0.5*(cp%thvbot(k+1)+  &
                                    cp%thvtop(k))*exn_k(cp%ps(k),Uw_p))
       cp%ufrc(k) = cp%umf(k)/(rhos0j*max(cp%wu(k), cpn%wmin))
       cflim = MIN(cpn%rmaxfrac, cp%maxcldfrac)
       if(cp%ufrc(k).gt.cflim        )then
          cp%ufrc(k) = cflim
          cp%umf (k) = cflim         * rhos0j * cp%wu(k)
          cp%fdr (k) = cp%fer(k)-log(cp%umf(k)/cp%umf(km1))/cp%dp(k)
       endif
       cp%pptr(k) = qrj*cp%umf(k)
       cp%ppti(k) = qsj*cp%umf(k)
       cp%pptn(k) = delta_qn*cp%umf(k)
!temperature (at full level)    
       t_mid = 0.5 * (cp%t(k) + cp%t(km1))
! virtual temperature (at full level)
       tv_mid = 0.5 * ( cp%t(k) * (1+Uw_p%zvir*(cp%qctu(k)-  &
                                               cp%qlu(k)-cp%qiu(k))) &
               + cp%t(km1) * (1+Uw_p%zvir*(cp%qctu(km1)-   &
                                             cp%qlu(km1)-cp%qiu(km1))) )
! air density (kg/m3)
       air_density = 0.5*(cp%ps(km1)+cp%ps(k)) / (Uw_p%rdgas*tv_mid)
       total_condensate = cp%qlu(k) + cp%qiu(k) + qrj + qsj ! kg/kg
       total_rain = qrj * air_density ! kg/m3
       total_snow = qsj * air_density ! kg/m3
       if (total_rain+total_snow > 0.) then
          do n=1,size(cp%tru,2)
             if (cpn%wetdep(n)%Lwetdep) then
                call wet_deposition_0D( cpn%wetdep(n)%Henry_constant, &
                                        cpn%wetdep(n)%Henry_variable, &
                                        cpn%wetdep(n)%frac_in_cloud, &
                                        cpn%wetdep(n)%alpha_r, &
                                        cpn%wetdep(n)%alpha_s, &
                                        t_mid, cp%ps(km1), cp%ps(k), &
                                        air_density, &
                                        total_condensate, total_rain, total_snow, &
                                        cp%tru(k,n), &
                                        cpn%wetdep(n)%Lgas, cpn%wetdep(n)%Laerosol, cpn%wetdep(n)%Lice, &
                                        delta_tracer )
                cp%tru_dwet(k,n) = -delta_tracer ! tracer source from wet deposition (negative=sink)
                cp%tru(k,n) = cp%tru(k,n) - delta_tracer ! adjust in-cloud concentration for wet dep
             end if
       end do
       end if
       
    enddo !End of Updraft Loop

    ltop = k 
    cp%umf (ltop) = 0.
    cp%pptr(ltop) = 0.
    cp%ppti(ltop) = 0.
    cp%pptn(ltop) = 0.
    cp%peff(ltop) = 0.

    cp%let    = let
    cp%ltop   = ltop
    cp%cldhgt = sd%z(ltop)-ac%zlcl
    !Restriction of convection too deep or too shallow
    if(cp%cldhgt.ge.cpn%cldhgt_max .or. ltop.lt.krel+2 .or.   &
                                                   let.le.krel+1) then
       return
    endif

    !convective scale height
    cp % cush=sd%z(ltop) - sd%zs(0)

    if (cpn%mixing_assumption.eq.2) then
    	 cp % cush = cpn%t00*Uw_p%rdgas/Uw_p%grav*log(cp%p(1)/cp%p(k))
    endif	 

    if (cpn%do_ppen) then !Calculate penetrative entrainment
       call penetrative_mixing_k(cpn, sd, Uw_p, cp) 
    else
       cp%fdr(ltop) = 1./sd%dp(ltop)
    end if

    if (cpn%mixing_assumption.eq.1) then
       cp%fdr(ltop)    = 1./sd%dp(ltop)
       cp%fdrsat(ltop) = 1.
    end if

  end subroutine cumulus_plume_k


!++++yim
  subroutine precipitation_k (zs, ps, hlu, qctu, qnu, cpn, qrj, qsj, &
                              hlu_new, qctu_new, qlu_new, qiu_new, &
                              clu_new, ciu_new, temp, doice, delta_qn, &
                              Uw_p, kbelowlet)       
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: zs, ps, hlu, qctu
    real,           intent(inout)    :: qnu, delta_qn
    real,           intent(inout) :: qrj, qsj, hlu_new, qctu_new,  &
                                     qlu_new, qiu_new, clu_new,  &
                                     ciu_new, temp
    logical,        intent(in)    :: doice, kbelowlet

    real    :: thj, qvj, qlj, qij, qse, thvj, nu, exnj,  &
               auto_th, leff, pcp, qctmp, deltaqc, auto_th2

    !Precip at the flux level
    call findt_k (zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice, &
                  Uw_p)  
    exnj=exn_k(ps,Uw_p)
    temp=thj*exnj-273.15
    if (temp.ge.0.0) then
       auto_th=cpn%auto_th0
    else
       auto_th=cpn%auto_th0*(1.0-temp/cpn%tcrit)
    end if
    auto_th=max(auto_th,0.0)

    if (.not.kbelowlet) auto_th=1.e10

    temp=temp+273.15

  if (.not.cpn%do_auto_aero) then
    qctmp   = qlj+qij;
    if (cpn%do_pdfpcp) then
       deltaqc = min(cpn%deltaqc0, auto_th)
       if (qctmp .lt. (auto_th - deltaqc)) then
          pcp = 0.
       else if (qctmp .lt. (auto_th + deltaqc)) then
          pcp = (qctmp + deltaqc - auto_th)**2./(4.*deltaqc)
       else
          pcp = qctmp-auto_th
       end if
    else
       pcp = max(qctmp-auto_th,0.)
    end if
    qctmp = 1./max(qctmp,1.e-28)
    qrj = pcp*qlj*qctmp
    qsj = pcp*qij*qctmp
    nu  = max(min((268. - temp)/20.,1.0),0.0)
 
    if (qlj.le.0) then
       delta_qn = -qnu
       qnu = 0
    else
       delta_qn = qnu * qrj * qctmp
       qnu      = qnu - delta_qn
    end if
 else
      auto_th2 = max (4.18667e-15*qnu*cpn%rad_crit**3., 0.0)
      if((qlj+qij).gt.auto_th)then
        qsj = (qlj+qij-auto_th)*qij/(qlj+qij)
        nu = max(min((268. - temp)/20.,1.0),0.0)
      else
        qsj = 0.0
        nu  = 0.0
      endif
 
      if(qlj .gt. auto_th2)then
        qrj = qlj-auto_th2
        nu = max(min((268. - temp)/20.,1.0),0.0)
!++++yim in-cloud removal of dropelts
        delta_qn = qnu
        if (qlj .gt. 0.) then
          qnu = qnu*auto_th2/qlj
        else
          qnu = 0.
        end if
        delta_qn = delta_qn - qnu
      else
        qrj = 0.0
        nu  = 0.0
        delta_qn = 0.
      endif
 end if

    leff     = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff
    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj
    clu_new  = qlu_new
    ciu_new  = qiu_new

    return
    
  end subroutine precipitation_k


  subroutine precip5_k (zs, ps, hlu, qctu, qnu, cpn, qrj, qsj, &
                              hlu_new, qctu_new, qlu_new, qiu_new, &
                              clu_new, ciu_new, temp, doice, delta_qn, &
                              Uw_p, kbelowlet, delp)       
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: zs, ps, hlu, qctu, delp
    real,           intent(inout)    :: qnu, delta_qn
    real,           intent(inout) :: qrj, qsj, hlu_new, qctu_new,  &
                                     qlu_new, qiu_new, clu_new,  &
                                     ciu_new, temp
    logical,        intent(in)    :: doice, kbelowlet

    real    :: thj, qvj, qlj, qij, qse, thvj, nu, exnj,  &
               auto_th, leff, pcp, qctmp, deltaqc, auto_th2, peff

    !Precip at the flux level
    call findt_k (zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice, &
                  Uw_p)  
    exnj=exn_k(ps,Uw_p)
    temp=thj*exnj-273.15
    if (temp.ge.0.0) then
       peff=cpn%peff_l*delp
    else
       peff=cpn%peff_i*delp
    end if
    peff=max(1.0-peff,0.0)

    if (.not.kbelowlet) peff=0.0

    temp=temp+273.15

    qctmp   = qlj+qij;
    pcp = max(qctmp*peff,0.)
    qctmp = 1./max(qctmp,1.e-28)
    qrj = pcp*qlj*qctmp
    qsj = pcp*qij*qctmp
    nu  = max(min((268. - temp)/20.,1.0),0.0)
 
    if (qlj.le.0) then
       delta_qn = -qnu
       qnu = 0
    else
       delta_qn = qnu * qrj * qctmp
       qnu      = qnu - delta_qn
    end if

    leff     = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff
    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj
    clu_new  = qlu_new
    ciu_new  = qiu_new

    return
    
  end subroutine precip5_k


  subroutine precip4_k (zs, ps, hlu, qctu, qnu, cpn, qrj, qsj, &
                              hlu_new, qctu_new, qlu_new, qiu_new, &
                              clu_new, ciu_new, temp, doice, delta_qn, &
                              Uw_p, kbelowlet)       
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: zs, ps, hlu, qctu
    real,           intent(inout)    :: qnu, delta_qn
    real,           intent(inout) :: qrj, qsj, hlu_new, qctu_new,  &
                                     qlu_new, qiu_new, clu_new,  &
                                     ciu_new, temp
    logical,        intent(in)    :: doice, kbelowlet

    real    :: thj, qvj, qlj, qij, qse, thvj, nu, exnj,  &
               auto_th, leff, pcp, qctmp, deltaqc, auto_th2, peff

    !Precip at the flux level
    call findt_k (zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice, &
                  Uw_p)  
    exnj=exn_k(ps,Uw_p)
    temp=thj*exnj-273.15
    if (temp.ge.0.0) then
       peff=cpn%peff_l
    else
!      peff=(1.0-cpn%peff)*min(temp/cpn%tcrit,1.0)
       peff=cpn%peff_i
    end if
!    peff=max(1.0-peff,0.0)

    if (.not.kbelowlet) peff=0.0

    temp=temp+273.15

    qctmp   = qlj+qij;
    pcp = max(qctmp*peff,0.)
    qctmp = 1./max(qctmp,1.e-28)
    qrj = pcp*qlj*qctmp
    qsj = pcp*qij*qctmp
    nu  = max(min((268. - temp)/20.,1.0),0.0)
 
    if (qlj.le.0) then
       delta_qn = -qnu
       qnu = 0
    else
       delta_qn = qnu * qrj * qctmp
       qnu      = qnu - delta_qn
    end if

    leff     = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff
    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj
    clu_new  = qlu_new
    ciu_new  = qiu_new

    return
    
  end subroutine precip4_k


  subroutine precip3_k (zs, ps, hlu, qctu, qnu, cpn, qrj, qsj, &
                              hlu_new, qctu_new, qlu_new, qiu_new, &
                              clu_new, ciu_new, temp, doice, delta_qn, &
                              Uw_p, kbelowlet)       
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: zs, ps, hlu, qctu
    real,           intent(inout)    :: qnu, delta_qn
    real,           intent(inout) :: qrj, qsj, hlu_new, qctu_new,  &
                                     qlu_new, qiu_new, clu_new,  &
                                     ciu_new, temp
    logical,        intent(in)    :: doice, kbelowlet

    real    :: thj, qvj, qlj, qij, qse, thvj, nu, exnj,  &
               auto_th, leff, pcp, qctmp, deltaqc, peff

    !Precip at the flux level
    call findt_k (zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice, &
                  Uw_p)  
    exnj=exn_k(ps,Uw_p)
    temp=thj*exnj-273.15
    if (temp.ge.0.0) then
       auto_th=cpn%auto_th0
    else
       auto_th=cpn%auto_th0*(1.0-temp/cpn%tcrit)
    end if
    auto_th=max(auto_th,0.0)

    if (.not.kbelowlet) auto_th=1.e10

    temp=temp+273.15

    qctmp   = qlj+qij; 
    if (cpn%do_pdfpcp) then
       deltaqc = min(cpn%deltaqc0, auto_th)
       if (qctmp .lt. (auto_th - deltaqc)) then
          pcp = 0.
       else if (qctmp .lt. (auto_th + deltaqc)) then
          pcp = (qctmp + deltaqc - auto_th)**2./(4.*deltaqc)
       else
          pcp = qctmp-auto_th
       end if
    else
       pcp = max(qctmp-auto_th,0.)
    end if

    if (temp.ge.0.0) then
       peff=cpn%peff_l
    else
       peff=cpn%peff_i
    end if
    pcp = qctmp*peff

    qctmp = 1./max(qctmp,1.e-28)
    qrj = pcp*qlj*qctmp
    qsj = pcp*qij*qctmp
    nu  = max(min((268. - temp)/20.,1.0),0.0)
    
    if (qlj.le.0) then
       delta_qn = -qnu
       qnu = 0
    else
       delta_qn = qnu * qrj * qctmp
       qnu      = qnu - delta_qn
    end if

    leff     = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff
    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj
    clu_new  = qlu_new
    ciu_new  = qiu_new

    return
    
  end subroutine precip3_k



  subroutine precip_new_k    (zs, ps, hlu, qctu, qnu, cpn, qrj, qsj, &
                              hlu_new, qctu_new, qlu_new, qiu_new, &
                              clu_new, ciu_new, temp, doice, delta_qn, &
                              Uw_p, kbelowlet, land, delt, wu)
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: zs, ps, hlu, qctu, land, delt, wu
    real,           intent(inout)    :: qnu, delta_qn
    real,           intent(inout) :: qrj, qsj, hlu_new, qctu_new,  &
                                     qlu_new, qiu_new, clu_new,  &
                                     ciu_new, temp
    logical,        intent(in)    :: doice, kbelowlet

    real    :: thj, qvj, qlj, qij, qse, thvj, nu, exnj,  &
               auto_th, leff, pcp, qctmp
    real    :: Nl, fliq, ql0, qi0

    !Precip at the flux level
    call findt_k (zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice, &
                  Uw_p)  
    exnj=exn_k(ps,Uw_p)
    temp=thj*exnj
    qctmp = qlj+qij; 
    if (qctmp > 0.0) then
       Nl = cpn%Nl_land*land + cpn%Nl_ocean*(1.-land) !Nl is mixing ratio=(1/m3)/(kg/m3)
       fliq = qlj/qctmp;

       !ql0=(4./3.*pi*1000) * (12.e-6)^3. * (100.e6) = 7.238e-4kg/kg
       ql0 =4188.79 * (cpn%r_thresh)**3. * Nl !4188.79=4/3*pi*1000; Nl=N/rho=N0/rho0;
       if (temp.ge.248.) then
          qi0 =cpn%qi_thresh
       else
          qi0 =cpn%qi_thresh*min(max(1.0+(248.-temp)/(cpn%tcrit+25.15),0.0),1.0)
       end if

       if (qij.eq.0.0) then
          auto_th=ql0
       else if (qlj.eq.0.0) then
          auto_th=qi0
       else
          auto_th=ql0*fliq+qi0*(1.-fliq)
       end if

       if (cpn%do_weffect) then
          auto_th=auto_th*(max(wu,1.))**cpn%weffect
       end if

       if (.not.kbelowlet) auto_th=1.e10

!       pcp = max(qctmp-auto_th,0.)*cpn%peff*delt/(1.+cpn%peff*delt)
!       pcp = max(qctmp-auto_th,0.)*cpn%peffdelt
       pcp = max(qctmp-auto_th,0.)

       qrj = pcp*fliq
       qsj = pcp*(1.-fliq)
       nu  = max(min((268. - temp)/20.,1.0),0.0)
    
       if (qlj.le.0) then
          delta_qn = -qnu
          qnu = 0
       else
          delta_qn = qnu * qrj / qctmp
          qnu      = qnu - delta_qn
       end if
       leff     = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
       qctu_new = qctu - (qrj + qsj)
       hlu_new  = hlu  + (qrj + qsj)*leff
       qlu_new  = qlj - qrj
       qiu_new  = qij - qsj
       clu_new  = qlu_new
       ciu_new  = qiu_new
    else
       qrj      = 0.
       qsj      = 0.
       qctu_new = qctu
       hlu_new  = hlu  
       qlu_new  = qlj 
       qiu_new  = qij 
       clu_new  = qlu_new
       ciu_new  = qiu_new
    end if

    return
    
  end subroutine precip_new_k



  subroutine micro_donner_k (cpn, zs, ps, hlu, qctu, zs1, qlu1, clu1, &
                             qiu1, ciu1, w1, cr12, pr12, qrj, qsj, &
                             qlu_new, clu_new, qiu_new, ciu_new, &
                             hlu_new, qctu_new, temp, doice, Uw_p)       
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    real,           intent(in)    :: zs, ps, hlu, qctu
    real,           intent(in)    :: zs1, qlu1, clu1, qiu1, ciu1, w1
    real,           intent(inout) :: qrj, qsj, cr12, pr12
    real,           intent(inout) :: qlu_new, clu_new, qiu_new, &
                                     ciu_new, hlu_new, qctu_new, temp
    logical,        intent(in)    :: doice

    real    :: thj, qvj, qlj, qij, qse, thvj, nu, leff
    real    :: dt_micro, rw1, cw1, drwa, drwb, flw, rw2, cw2, pw2, dcw

    call findt_k (zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice, &
                  Uw_p)   
    temp = thj*exn_k(ps,Uw_p)
    if (doice) then
      nu   = max(min((268. - temp)/20.,1.0),0.0)
    else
      nu = 0.
    endif

    leff = (1.-nu)*Uw_p%HLv + nu*Uw_p%HLs
    if (qlj+qij .gt. 0.0) then
       flw  = qlj/(qlj+qij)
    else
       qrj      = 0.
       qsj      = 0.
       qlu_new  = 0.
       qiu_new  = 0.
       qctu_new = qctu
       hlu_new  = hlu
       clu_new  = 0.
       ciu_new  = 0.
       return
    end if

    cw1 = clu1 + ciu1

    rw1 = qlu1 + qiu1 - cw1
    rw1 = max(rw1, 0.0)

    cw2 = qlj + qij - rw1

    dcw = cw2 - cw1; 

    dt_micro = (zs - zs1) / w1

    cr12 = dcw/dt_micro

    drwa = cpn%auto_rate * (cw2 - cpn%auto_th0) * dt_micro

    drwa=min(max(drwa, 0.0), cw2-cpn%auto_th0)

    cw2 = cw2 - drwa

    drwb = 5.26e-03 * cw2 * (rw1**0.875) * dt_micro

    drwb=min(max(drwb, 0.0), cw2)

    cw2 = cw2 - drwb

    rw2 = rw1 + drwa + drwb

    rw2 = max(rw2, 0.0)

    pw2 = 5.1*(rw2**1.125)*dt_micro/100.

    pw2 =min(pw2, rw2)

    pw2 = min(pw2, qlj + qij)

    pr12=pw2/dt_micro
   
    rw2 =rw2-pw2

    qrj = pw2*flw
    qsj = pw2*(1.-flw)

    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff

    cw2 = max(qlu_new + qiu_new -rw2, 0.0)
    clu_new = cw2*flw
    ciu_new  = cw2*(1. - flw)

    if (qlu_new .lt. 0. .or. qiu_new .lt. 0. .or. clu_new .lt. 0.0 .or.&
        ciu_new .lt. 0.0) then
       print*, qlu_new, qiu_new, clu_new, ciu_new,   &
                             qrj, qsj, qlj, qij, '??????????????????'
    end if
    return
    
  end subroutine micro_donner_k



!#####################################################################
!#####################################################################

  subroutine penetrative_mixing_k(cpn, sd, Uw_p, cp)
    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(uw_params),   intent(inout) :: Uw_p
    type(cplume),   intent(inout) :: cp

    integer :: k, ltop, let
    real    :: rhos0j, bogtop, bogbot
    real    :: aquad, bquad, cquad, xs1, xs2, ppen
    real    :: thj, qvj, qse, thvj, qctulet, hlulet, umflet
    real    :: dqct1, dqct2, qctflxkm1, tmp

    ltop=cp%ltop
    let =cp%let

    cp % emf (ltop)  = 0.0
    cp % fdr (ltop)  = 0.0
    cp % hlu (ltop)  = cp % hl (ltop)
    cp % qctu(ltop)  = cp % qct(ltop)
    cp % qlu (ltop)  = cp % ql (ltop)
    cp % qiu (ltop)  = cp % qi (ltop)
    cp % tru (ltop,:)= cp % tr (ltop,:)

    qctulet = cp%qctu(let)
    hlulet  = cp%hlu (let)
    umflet  = cp%umf (let)

    do k=ltop-1,let,-1

       cp%fdr(k)  = 0.
       cp%qlu(k)  = cp%ql(k)
       cp%qiu(k)  = cp%qi(k)

       rhos0j = cp%ps(k) /(Uw_p%rdgas*0.5*   &
                     (cp%thvbot(k+1)+cp%thvtop(k))*exn_k(cp%ps(k),Uw_p))
       if(k.eq.ltop-1)then
          !Calculate ppen
!!$          bogbot = (cp%thvu(k)/cp%thvbot(ltop) - 1.)/cpn%rbuoy
!!$          if(bogbot.gt.0.)then
!!$             bogbot =  bogbot /cpn%rbuoy
!!$          else
!!$             bogbot =  bogbot*cpn%rbuoy
!!$          endif
!!$          bogtop = (cp%thvu(ltop)/cp%thvtop(ltop) - 1.)/cpn%rbuoy
!!$          if(bogtop.gt.0.)then
!!$             bogtop =  bogtop /cpn%rbuoy
!!$          else
!!$             bogtop =  bogtop*cpn%rbuoy
!!$          endif

          bogbot = (cp%thvu(k)/cp%thvbot(ltop) - 1.)
          bogbot =  bogbot*cpn%rbuoy
          bogtop = (cp%thvu(ltop)/cp%thvtop(ltop) - 1.)
          bogtop =  bogtop*cpn%rbuoy

          aquad = (bogtop - bogbot) / (cp%ps(ltop)-cp%ps(k))
          bquad = 2*bogbot
          cquad = -cp%wu(k) * cp%ps(k) /   &
                      (Uw_p%rdgas*cp%thvbot(ltop)*exn_k(cp%ps(k),Uw_p))
          call roots(aquad,bquad,cquad,xs1,xs2)
          if(xs1.le.0..and.xs2.le.0.)then
             ppen = max(xs1,xs2)
          else
             ppen = min(xs1,xs2)
          endif
          ppen = min(0.,max(-cp%dp(k+1),ppen))
          if(xs1.eq.-9.99e33.or.xs2.eq.-9.99e33) ppen=0.
          !Calculate returning mass flux
          cp%emf (k)=max(cp%umf(k)*ppen*cp%rei(ltop)*  &
                                                 cpn%rpen,-0.1*rhos0j)
          cp%hlu (k)=cp%hl (ltop)+sd%sshl (ltop)*(cp%ps(k)-cp%p(ltop))
          cp%qctu(k)=cp%qct(ltop)+sd%ssqct(ltop)*(cp%ps(k)-cp%p(ltop))
          cp%tru(k,:)=cp%tr(ltop,:)+sd%sstr(ltop,:)*   &
                                                 (cp%ps(k)-cp%p(ltop))
       else
          cp%emf (k)=max(cp%emf(k+1)-cp%umf(k)*cp%dp(k+1)*  &
                                      cp%rei(k+1)*cpn%rpen,-0.1*rhos0j)
          if (cp%emf(k).ne.0) then
            cp%hlu (k)=(cp%hlu (k+1)*cp%emf(k+1)+cp%hl (k+1)* &
                                        (cp%emf(k)-cp%emf(k+1)))/cp%emf(k)
            cp%qctu(k)=(cp%qctu(k+1)*cp%emf(k+1)+cp%qct(k+1)*  &
                                        (cp%emf(k)-cp%emf(k+1)))/cp%emf(k)
            cp%tru(k,:)=(cp%tru(k+1,:)*cp%emf(k+1)+cp%tr(k+1,:)* &
                                        (cp%emf(k)-cp%emf(k+1)))/cp%emf(k)
          endif
       endif
       cp%umf(k)=0.0 
    enddo

    k=let
    cp%fdr (k) = 1./sd%dp(k)
    
 if (cpn%do_pmadjt) then
    dqct1=cp%qctu(k-1)-(cp%qct(k)  +sd%ssqct(k)  *(sd%ps(k-1)-sd%p(k)))
    dqct2=cp%qctu(k)  -(cp%qct(k)  +sd%ssqct(k)  *(sd%ps(k)-sd%p(k)))
    qctflxkm1=cp%umf(k-1)*dqct1
    if ((cp%emf(k)*dqct2.gt.qctflxkm1).and.(qctflxkm1.gt.0.).and.(dqct2.lt.0)) then
       tmp=qctflxkm1/dqct2/cp%emf(k)
       cp%emf=cp%emf*tmp
    end if

    tmp=umflet-cp%emf(k)
    if (tmp.gt.0) then
       qctulet = (umflet*qctulet - cp%emf(k)*cp%qctu(k))/tmp
       hlulet  = (umflet*hlulet  - cp%emf(k)*cp%hlu (k))/tmp
    else
       qctulet = cp%qctu(k)
       hlulet  = cp%hlu (k)
    end if
    call findt_k (cp%zs(k), cp%ps(k), hlulet, qctulet, thj, qvj,       &
                  cp%qlu(k), cp%qiu(k), qse, thvj, cpn%do_ice, Uw_p)
 end if

  end subroutine penetrative_mixing_k

!#####################################################################
!#####################################################################

  subroutine cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)
  
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    type(sounding), intent(in)    :: sd
    type(cplume),   intent(inout) :: cp
    type(ctend),    intent(inout) :: ct
    logical,        intent(in)    :: do_coldT

    integer :: k, krel, ltop, kp1, km1, ktop, i
    real    :: dpsum, qtdef, hldef, umftmp, qlutmp, qiutmp, qnutmp, fdrtmp
    real, dimension(size(cp%tr,2)) :: trdef
    real    :: dpevap, x1, x2, x3, xx1, xx2, xx3, q1, q2, emftmp
    real    :: dqt

    call ct_clear_k (ct);

    krel=cp%krel
    ltop=cp%ltop

    ! Calculate Fluxes of heat, moisture, momentum
    dpsum = 0.0
    do k = 1, krel-1
       dpsum = dpsum + sd%dp(k)
    enddo

    qtdef = max(0.,cp%umf(krel)*(cp%qctu(krel) - cp%qct(krel)))
    trdef(:) = max(0.,cp%umf(krel)*(cp%tru(krel,:) - cp%tr(krel,:)))
    !yy1  = min(0.,umf(krel)*(thcu(krel) - thc0(krel)))
    hldef = min(0.,cp%umf(krel)*(cp%hlu (krel) - cp%hl (krel)))
    do k=1,krel-1
!      ct%hlflx (k)=0.0; 
       ct%hlflx (k)=0.0; !ct%hlflx (k-1) + hldef*sd%dp(k)/dpsum;
!      ct%thcflx(k)=0.0; !thcflx(k)=thcflx(k-1) + yy1*dp(k)/dpsum
       ct%qctflx(k)=ct%qctflx(k-1) + qtdef*sd%dp(k)/dpsum;
       ct%thcflx(k)=0.0;
       ct%trflx(k,:)=ct%trflx(k-1,:) + trdef(:)*sd%dp(k)/dpsum
       ct%qlflx(k)=0.0;
       ct%qiflx(k)=0.0;
       ct%qnflx(k)=0.0;
       ct%umflx(k)=0.0;
       ct%vmflx(k)=0.0;
       cp%pptr (k)=0.0; 
       cp%ppti (k)=0.0;
       cp%pptn (k)=0.0;
    enddo

    if (cpn%do_qctflx_zero) then
       do k=1,krel-1
          ct%qctflx(k) =0.;
          ct%trflx (k,:)=0.;
       end do
    end if

    do k = krel,ltop-1 !pzhu do k = krel,ltop
       kp1 = k+1

       ct%hlflx (k)= cp%umf(k)*(cp%hlu (k)-(cp%hl (kp1)+  &
                       sd%sshl (kp1)*(sd%ps(k)-sd%p(kp1)))) + &
                       cp%emf(k) * (cp%hlu (k)-  &
                           (cp%hl (k)+sd%sshl (k)*(sd%ps(k)-sd%p(k))))
       ct%thcflx(k)= cp%umf(k)*(cp%thcu(k)-(cp%thc(kp1)+  &
                        sd%ssthc(kp1)*(sd%ps(k)-sd%p(kp1)))) + &
                       cp%emf(k) * (cp%thcu(k)-  &
                            (cp%thc(k)+sd%ssthc(k)*(sd%ps(k)-sd%p(k))))
       ct%qctflx(k)= cp%umf(k)*(cp%qctu(k)-(cp%qct(kp1)+  &
                        sd%ssqct(kp1)*(sd%ps(k)-sd%p(kp1)))) + &
                       cp%emf(k) * (cp%qctu(k)-  &
                            (cp%qct(k)+sd%ssqct(k)*(sd%ps(k)-sd%p(k))))
       
       ct%umflx(k) =cp%umf(k) * (cp%uu(k) - cp%u(kp1))  +   &
                                      cp%emf(k) * (cp%u(kp1) -cp%u(k))
       ct%vmflx(k) =cp%umf(k) * (cp%vu(k) - cp%v(kp1))  +  &
                                      cp%emf(k) * (cp%v(kp1) -cp%v(k))

       ct%qlflx(k) =cp%umf(k) * (cp%qlu(k)- cp%ql(kp1)) +   &
                                      cp%emf(k) * (cp%ql(kp1)-cp%ql(k))
       ct%qiflx(k) =cp%umf(k) * (cp%qiu(k)- cp%qi(kp1)) +   &
                                      cp%emf(k) * (cp%qi(kp1)-cp%qi(k))
       ct%qaflx(k) =cp%umf(k) * (1.       - cp%qa(kp1)) +   &
                                       cp%emf(k) * (cp%qa(kp1)-cp%qa(k))
       ct%qnflx(k) =cp%umf(k) * (cp%qnu(k)- cp%qn(kp1)) +  &
                                      cp%emf(k) * (cp%qn(kp1)-cp%qn(k))

!++++yim
       ct%trflx(k,:) =cp%umf(k) * (cp%tru(k,:)- (cp%tr(kp1,:)+  &
                             sd%sstr(kp1,:)*(sd%ps(k)-sd%p(kp1))) ) + &
                       cp%emf(k) * (cp%tru(kp1,:)-   &
                         (cp%tr(k,:)+sd%sstr(k,:)*(sd%ps(k)-sd%p(k))) )
       ct%qvflx(k) =ct%qctflx(k)-ct%qlflx(k)-ct%qiflx(k)
    enddo

    do k = 2,ltop
       km1 = k-1
       kp1 = k+1
       x1=sd%p(k)   -sd%ps(k)
       x2=sd%ps(km1)-sd%p (k)
       x3=sd%ps(km1)-sd%ps(k)
       umftmp = (cp%umf(km1)*x1+cp%umf(k)*x2)/x3
       emftmp = (cp%emf(km1)*x1+cp%emf(k)*x2)/x3
       qlutmp = (cp%qlu(km1)*x1+cp%qlu(k)*x2)/x3
       qiutmp = (cp%qiu(km1)*x1+cp%qiu(k)*x2)/x3
       qnutmp = (cp%qnu(km1)*x1+cp%qnu(k)*x2)/x3

       fdrtmp = cp%fdrsat(k)*cp%fdr(k)
       ct%qldet(k) =  Uw_p%grav*(umftmp       )*(fdrtmp*(qlutmp-sd%ql(k)))
       ct%qidet(k) =  Uw_p%grav*(umftmp       )*(fdrtmp*(qiutmp-sd%qi(k)))
       ct%qadet(k) =  Uw_p%grav*(umftmp-emftmp)*(fdrtmp*(1.    -sd%qa(k)))
       ct%qndet(k) =  Uw_p%grav*(umftmp       )*(fdrtmp*(qnutmp-sd%qn(k)))

       x1 =sd%ps(k)  -sd%p (kp1)
       x2 =sd%p (k)  -sd%ps(k)
       x3 =sd%p (k)  -sd%p (kp1)
       xx1=sd%ps(km1)-sd%p (k)
       xx2=sd%p (km1)-sd%ps(km1)
       xx3=sd%p (km1)-sd%p (k)

       q2=(sd%ql(k)  *x1  + sd%ql(kp1)*x2  )/x3
       q1=(sd%ql(km1)*xx1 + sd%ql(k)  *xx2 )/xx3
       ct%qlten(k) = - Uw_p%grav*(umftmp * (sd%ql(k)-q2      )/x2 + &
                                  emftmp * (q1      -sd%ql(k))/xx1 )

       q2=(sd%qi(k)  *x1  + sd%qi(kp1)*x2  )/x3
       q1=(sd%qi(km1)*xx1 + sd%qi(k)  *xx2 )/xx3
       ct%qiten(k) = - Uw_p%grav*(umftmp * (sd%qi(k)-q2      )/x2 + &
                                  emftmp * (q1      -sd%qi(k))/xx1 )

       q2=(sd%qa(k)  *x1  + sd%qa(kp1)*x2  )/x3
       q1=(sd%qa(km1)*xx1 + sd%qa(k)  *xx2 )/xx3
       ct%qaten(k) = - Uw_p%grav*(umftmp * (sd%qa(k)-q2      )/x2 + &
                                  emftmp * (q1      -sd%qa(k))/xx1 )

       q2=(sd%qn(k)  *x1  + sd%qn(kp1)*x2  )/x3
       q1=(sd%qn(km1)*xx1 + sd%qn(k)  *xx2 )/xx3
       ct%qnten(k) = - Uw_p%grav*(umftmp * (sd%qn(k)-q2      )/x2 + &
                                  emftmp * (q1      -sd%qn(k))/xx1 )
    end do

    ct%qlten = ct%qlten + ct%qldet 
    ct%qiten = ct%qiten + ct%qidet
    ct%qaten = ct%qaten + ct%qadet
    ct%qnten = ct%qnten + ct%qndet

    if (cpn%do_detran_zero) then
       ct%qlten = 0.
       ct%qiten = 0.
       ct%qaten = 0.
       ct%qnten = 0.
    end if

    ! Calculate model tendencies
    do k = 1,ltop
       km1 = k-1
       kp1 = k+1
       ct%uten (k) = (ct%umflx(km1)-ct%umflx(k))*Uw_p%grav/sd%dp(k)
       ct%vten (k) = (ct%vmflx(km1)-ct%vmflx(k))*Uw_p%grav/sd%dp(k)

       ct%thcten(k) = (ct%thcflx(km1) - ct%thcflx(k))*   &
                                                   Uw_p%grav/sd%dp(k)
!      ct%hlten (k) = (ct%hlflx (km1) - ct%hlflx (k)   &
!            +(cp%pptr(k) + cp%ppti(k)) *sd%leff(k))*Uw_p%grav/sd%dp(k)
       ct%hlten (k) = (ct%hlflx (km1) - ct%hlflx (k) +   &
                       Uw_p%HLv*cp%pptr(k) + Uw_p%HLs*cp%ppti(k))* &
                                                      Uw_p%grav/sd%dp(k)
       ct%qctten(k) = (ct%qctflx(km1) - ct%qctflx(k) -   &
                           cp%pptr(k) - cp%ppti(k)) *Uw_p%grav/sd%dp(k)
       ct%qvten (k) = ct%qctten(k)-ct%qlten(k)-ct%qiten(k)
       ct%pflx  (k) = cp%pptr(k) + cp%ppti(k)

!      ct%tten  (k) = (ct%hlten(k)+sd%leff(k)*  &
!                                (ct%qlten(k)+ct%qiten(k)))/Uw_p%cp_air
       ct%tten  (k) = (ct%hlten(k)+Uw_p%HLv*ct%qlten(k)+    &
                                      Uw_p%HLs*ct%qiten(k))/Uw_p%cp_air
       ct%trten(k,:) = (ct%trflx(km1,:)-ct%trflx(k,:))*   &
                                                    Uw_p%grav/sd%dp(k)
       ct%trwet(k,:) = cp%tru_dwet(k,:) *cp%umf(k)*Uw_p%grav/sd%dp(k)

    enddo

    do k = cp%let,ltop
       if (ct%qctten(k).gt.0 .and. ct%qvten(k).lt.0) then
          qlutmp     =(1.-sd%nu(k))*ct%qvten(k)
          qiutmp     =    sd%nu(k) *ct%qvten(k)
          ct%qlten(k)=ct%qlten(k)+qlutmp
          ct%qiten(k)=ct%qiten(k)+qiutmp
          ct%qvten(k)=0.
          ct%tten (k)=ct%tten(k)+(Uw_p%HLv*qlutmp+Uw_p%HLs*qiutmp)/Uw_p%cp_air
       end if
    end do

    ct%dtint=0.; ct%dqint=0.; ct%conint=0.; ct%freint=0.; !dpsum=0.;

    ct%qlflx(ltop)=(1.-cpn%atopevap)*ct%qlflx(ltop-1);
    ct%qiflx(ltop)=(1.-cpn%atopevap)*ct%qiflx(ltop-1);
    ktop=ltop

    do k = 1,ktop
       km1 = k-1
       ct%qldiv(k) = -(ct%qlflx(k) - ct%qlflx(km1) +   &
                                       cp%pptr(k))* Uw_p%grav/sd%dp(k)
       ct%qidiv(k) = -(ct%qiflx(k) - ct%qiflx(km1) +   &
                                       cp%ppti(k))* Uw_p%grav/sd%dp(k)
    end do

    do k = 1,ltop
       km1 = k-1
       ct%qvdiv(k) = -(ct%qvflx(k) - ct%qvflx(km1))* Uw_p%grav/sd%dp(k)
       ct%dtint    = ct%dtint    + ct%tten (k)*Uw_p%cp_air/sd%leff(k)* &
                                                    sd%dp(k)/Uw_p%grav
       ct%dqint    = ct%dqint    + ct%qvten(k)* sd%dp(k)/Uw_p%grav
       ct%conint   = ct%conint   + ct%qldiv(k)* sd%dp(k)/Uw_p%grav
       ct%freint   = ct%freint   + ct%qidiv(k)* sd%dp(k)/Uw_p%grav
    end do

    if (do_coldT) then
       do k = 1,ltop
          if (sd%coldT) then
             ct%tten(k)=ct%tten(k)+cp%pptr(k)*Uw_p%HLf*Uw_p%grav/  &
                                                   sd%dp(k)/Uw_p%cp_air
             cp%ppti(k)=cp%ppti(k)+cp%pptr(k)
             cp%pptr(k)=0.
          else
             ct%tten(k)=ct%tten(k)-cp%ppti(k)*Uw_p%HLf*Uw_p%grav/  &
                                                   sd%dp(k)/Uw_p%cp_air
             cp%pptr(k)=cp%pptr(k)+cp%ppti(k)
             cp%ppti(k)=0.
          end if
          ct%snow  = ct%snow  + cp%ppti(k)
          ct%rain  = ct%rain  + cp%pptr(k)
       end do
       if (cpn%do_pevap .and. ct%snow+ct%rain > 0.) then
          call precip_evap (sd, cp, cpn, ct, Uw_p, dpevap)
          ct%tten (:)=ct%tten (:)+ct%tevap(:)
          ct%qvten(:)=ct%qvten(:)+ct%qevap(:)
          if (sd%coldT) then
             ct%snow  = ct%snow - dpevap
          else
             ct%rain  = ct%rain - dpevap
          end if
       end if
    end if

    if(cpn%do_pnqv) then
       do k = sd%kmax,3,-1
          dqt  =  ct%qvten(k) * sd%delt
          if (dqt.lt.0 .and. sd%qv(k)+dqt.lt.1.e-10) then
             fdrtmp = -(sd%qv(k)-1.e-10)/sd%delt - ct%qvten(k)
             ct%qvten(k) = ct%qvten(k) + fdrtmp
             dpsum = 0.0
             do i = k-1,1,-1
                dpsum = dpsum + sd%dp(i)
             enddo
             do i = k-1,1,-1
                ct%qvten(i) = ct%qvten(i) - fdrtmp*sd%dp(k)/dpsum
             enddo
          end if
       end do
    end if

    ct%dting=0.;
    ct%denth=0.; ct%dqtmp=0.; ct%uav=0.;ct%vav=0.; dpsum=0.;
    do k = 1,sd%kmax! ltop
       ct%denth = ct%denth + (Uw_p%cp_air*ct%tten(k)-Uw_p%HLv*ct%qlten(k)-Uw_p%HLs*ct%qiten(k))*sd%dp(k)/Uw_p%grav
       ct%dqtmp = ct%dqtmp + (ct%qvten(k) + ct%qlten(k) + ct%qiten(k))*sd%dp(k)/Uw_p%grav
       ct%uav   = ct%uav   + ct%uten(k)*sd%dp(k)/Uw_p%grav
       ct%vav   = ct%vav   + ct%vten(k)*sd%dp(k)/Uw_p%grav
       ct%dting = ct%dting + Uw_p%cp_air*ct%tten(k)*sd%dp(k)/Uw_p%grav
!       dpsum    = dpsum    + sd%dp(k)
     end do
     ct%denth = ct%denth - Uw_p%HLv*ct%rain - Uw_p%HLs*ct%snow
     ct%dqtmp = ct%dqtmp + ct%rain + ct%snow

  end subroutine cumulus_tend_k

!#####################################################################
!#####################################################################

 subroutine roots(a,b,c,r1,r2)
   real a,b,c,r1,r2,q
   
   if(a.eq.0)then            ! form b*x + c = 0
      if(b.eq.0)then         ! failure: c = 0
         r1 = -9.99e33
      else                   ! b*x + c = 0
         r1 = -c / b
      endif
      r2 = r1
   else
      if(b.eq.0.)then        ! form a*x**2 + c = 0
         if(a*c.gt.0.)then   ! failure: x**2 = -c/a < 0
            r1 =  -9.99e33
         else                ! x**2 = -c/a
            r1 = sqrt(-c/a)
         endif
         r2 = -r1
      else 
         if((b**2. - 4.*a*c).lt.0.)then ! failure, no real(r8) roots
            r1 =  -9.99e33
            r2 = -r1
         else
            q = - 0.5 * ( b + sign(1.0,b) * sqrt(b**2. - 4.*a*c) )
            r1 = q/a
            r2 = c/q
         endif
      endif
   endif
   return
 end subroutine roots

!#####################################################################
!#####################################################################

  SUBROUTINE precip_evap (sd, cp, cpn, ct, Uw_p, dpevap)

    implicit none
    
    type(sounding), intent(in)      :: sd
    type(cplume),   intent(in)      :: cp
    type(cpnlist),  intent(in)      :: cpn
    type(ctend),    intent(inout)   :: ct
    type(uw_params),intent(inout)   :: Uw_p
    real,           intent(inout)   :: dpevap

    real, parameter :: cem    = 0.054
    real, parameter :: ceta   = -544.0E-6
    real, parameter :: d622   = 287.04/461.50
    real, parameter :: d378   = 1.0-d622

    real, dimension(size(sd%t)) :: mass, temp_new, qvap_new, pptp, pflx, pflx_evap
    real    :: prec, def, evef, prec_mmph, pfac, emx, dpcu, HL, dqs, qs
    real    :: hcevap, cfrac

    integer :: k, ier

    cfrac     = cpn%cfrac
    hcevap    = cpn%hcevap

    pflx      = 0.0
    pflx_evap = 0.0
    qvap_new  = sd%qv
    temp_new  = sd%t

    if (sd%coldT) then
       HL=Uw_p%HLs
    else
       HL=Uw_p%HLv
    end if

    !pptp (unit: kg/m2 or mm) - precipitated water in the layer dp
    pptp   = (cp%pptr(:)+cp%ppti(:))*sd%delt
    mass   = sd%dp/Uw_p%grav
    dpcu   = 0.0
    dpevap = 0.0
    do k = cp%ltop, 1, -1
       dpcu = dpcu + pptp(k)
       prec = MAX(dpcu - dpevap, 0.0 )

! --- Compute precipitation efficiency factor
       prec_mmph = prec * 3600.0 / sd%delt
       pfac      = SQRT( sd%p(k) / sd%ps(0) )
       emx       = SQRT( cem * cfrac * prec_mmph * pfac )   
       evef      = 1.0 - EXP( ceta * sd%delt * emx ) 

       def=0. !Evaporate precip where needed
       if ( sd%rh(k) <= hcevap .and. prec > 0.0 ) then
          call compute_qs_k (sd%t(k), sd%p(k), Uw_p%epsilo, Uw_p%zvir, &
                             qs, ier, dqsdT=dqs) 
          def=(hcevap*sd%qs(k) - sd%qv(k))/(1.+(HL*hcevap*dqs/Uw_p%Cp_Air ))
          def=evef*def
          def=MIN( def, prec/mass(k) )
          def=MAX( def, 0.0)
       else
          def=0.0
       end if
       pflx_evap(k)= def * mass(k)
       dpevap      = dpevap + pflx_evap(k)
       qvap_new(k) = sd%qv(k) + def
       temp_new(k) = sd%t (k) - (def * HL/Uw_p%Cp_Air)
       pflx    (k) = prec
    end do
    dpevap      = min(dpevap, dpcu) / sd%delt
    ct%tevap(:) = (temp_new(:) - sd%t (:))/sd%delt
    ct%qevap(:) = (qvap_new(:) - sd%qv(:))/sd%delt
  
  end SUBROUTINE PRECIP_EVAP

!#####################################################################
!#####################################################################

end MODULE CONV_PLUMES_k_MOD

