
MODULE DEEP_CONV_MOD

  use      fms_mod,         only : write_version_number
  use      Constants_Mod,   ONLY : tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use  conv_utilities_mod,  only : uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, &
                                   adicloud, sounding, uw_params

  use  conv_plumes_k_mod,   only : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist

  use  conv_closures_mod,   only : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure


!---------------------------------------------------------------------
  implicit none
  private
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: deep_conv.F90,v 19.0 2012/01/06 20:26:02 fms Exp $'
  character(len=128) :: tagname = '$Name: siena_201207 $'

!-------  interfaces --------

  public  :: dpconv0, dpconv1, dpconv2, dpconv3, DEEP_CONV_INIT, DEEP_CONV_END

  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'deep_conv'

  public deepc
  type deepc
     real, dimension(7)  :: rkm_dp
     real, dimension(7)  :: rat_dp
     real, dimension(50) :: rkm
     real, dimension(50) :: rat
     real, dimension(50) :: hgt
     real    :: omeg_th
     real    :: cape_th 
     real    :: tau_dp  
     real    :: cbmf_d  
     real    :: cwfn_d
     real    :: deepdepth
     integer :: ideep_closure
     integer :: mixing_assumption
     logical :: do_generation
     logical :: do_ppen
     logical :: do_pevap
  end type deepc

contains

!#####################################################################
!#####################################################################

  subroutine cpn_copy(cpn, dpn)
    type(cpnlist), intent(in)    :: cpn
    type(cpnlist), intent(inout) :: dpn

    dpn % rle              = cpn % rle
    dpn % rpen             = cpn % rpen
    dpn % rmaxfrac         = cpn % rmaxfrac
    dpn % wmin             = cpn % wmin
    dpn % rbuoy            = cpn % rbuoy
    dpn % rdrag            = cpn % rdrag  
    dpn % frac_drs         = cpn % frac_drs
    dpn % bigc             = cpn % bigc    
    dpn % auto_th0         = cpn % auto_th0
    dpn % auto_rate        = cpn % auto_rate
    dpn % tcrit            = cpn % tcrit  
    dpn % cldhgt_max       = cpn % cldhgt_max
    dpn % do_ice           = cpn % do_ice
    dpn % do_ppen          = cpn % do_ppen
    dpn % do_pevap         = cpn % do_pevap
    dpn % mixing_assumption= cpn % mixing_assumption
    dpn % mp_choice        = cpn % mp_choice
    dpn % do_forcedlifting = cpn % do_forcedlifting
    dpn % atopevap         = cpn % atopevap
    dpn % wtwmin_ratio     = cpn % wtwmin_ratio

  end subroutine cpn_copy

!#####################################################################
!#####################################################################

  subroutine dpconv0(dpc, cpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
       omeg_avg, rkm_sh, cp1, ct1, cbmf_deep, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(in)     :: cpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp, cp1
    type(ctend),     intent(inout)  :: ct, ct1
    real,            intent(inout)  :: cbmf_deep, ocode, rkm_sh, omeg_avg
    integer,            intent(out)   :: ier
    character(len=256), intent(out)   :: ermesg


    type(cpnlist) :: dpn
    real          :: rkm_dp, zcldtop

    if ( (ocode.ne.0) .or. (omeg_avg .gt.dpc%omeg_th)) then
       ocode=6; return
    end if

    call cpn_copy(cpn, dpn)
    dpn % do_ppen   = dpc % do_ppen
    dpn % do_pevap  = dpc % do_pevap
    rkm_dp  = dpc%rkm_dp(1) *  rkm_sh
    zcldtop = sd%z(cp%ltop)

    call cp_clear_k(cp1);
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=6; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
 
  end subroutine dpconv0


!#####################################################################
!#####################################################################

  subroutine dpconv1(dpc, cpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, sd1, ac1, &
       cc1, cp1, ct1, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(in)     :: cpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(inout)  :: ac
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cclosure),  intent(inout)  :: cc,cc1
    type(cplume),    intent(inout)  :: cp,cp1
    type(ctend),     intent(inout)  :: ct,ct1
    real,            intent(inout)  :: ocode
    integer,            intent(out)   :: ier
    character(len=256), intent(out)   :: ermesg

    integer :: i, ksrc
    real    :: cbmf0, cbmfs, cbmf_max, dcape, scaleh, wrel, tmp, cbmf, rkm
    real    :: zsrc, psrc, thcsrc, hlsrc, qctsrc, pdeet1, pdeet2

    type(cpnlist) :: dpn

    call cpn_copy(cpn, dpn)

    dpn % do_ppen          = .false.
    dpn % rmaxfrac         = 1000000.
    dpn % rbuoy            = 0.66666
    dpn % rdrag            = 3.0
    dpn % auto_th0         = 0.5e-3
    dpn % auto_rate        = 1.0e-3
    dpn % mixing_assumption= 1
    dpn % mp_choice        = 1
    dpn % do_forcedlifting = .true.

!!$    zsrc  =sd%zs (1);
!!$    psrc  =sd%ps (1); thcsrc=sd%thc(1)
!!$    hlsrc =cc%xhlsrc
!!$    qctsrc=cc%xqtsrc
!!$    call adi_cloud(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, .false., do_ice, ac1)
!!$    cc % dcape=(ac1%cape-cc%xcape)/sd%delt

    pdeet1=ac%plfc - ac%plnb
    pdeet2=ac%plfc - ac%plcl

    if  ((ac%cape  <= dpc%cape_th )                    .or.  &
!        (cc%dcape <= 0.  .and. dpc%do_dcape_closure ) .or.  &
         (pdeet1   <= 500.e02                        ) .or.  &
         (ac%cin   >= 100.))                           then
       dpc%cbmf_d=0.; 
       ocode=6; 
       return
    end if
 
    ksrc=2
    zsrc  =sd%zs (ksrc);
    psrc  =sd%ps (ksrc);    thcsrc=sd%thc(ksrc)
    qctsrc=sd%qct(ksrc)
    hlsrc =sd%hl (ksrc)
    call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, .false., do_ice, ac)

    cbmf0=0.0001; dpc%cbmf_d=0.; wrel=0.5; scaleh=1000.

    cbmf=1000000.*wrel;
    do i=1, size(dpc%rkm_dp(:))
       call ct_clear_k(ct)
       call cp_clear_k(cp)
       rkm = dpc%rkm_dp(i)
       call cumulus_plume_k(dpn, sd, ac, cp, rkm, cbmf, wrel, scaleh, Uw_p, ier, ermesg)
       if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
          dpc%cbmf_d=0.; ocode=6; return
       end if
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct, do_coldT)
       ct1%tten  = ct1%tten  + ct%tten  * dpc%rat_dp(i)
       ct1%qvten = ct1%qvten + ct%qvten * dpc%rat_dp(i)
       ct1%qlten = ct1%qlten + ct%qlten * dpc%rat_dp(i)
       ct1%qiten = ct1%qiten + ct%qiten * dpc%rat_dp(i)
       ct1%qaten = ct1%qaten + ct%qaten * dpc%rat_dp(i)
       ct1%qnten = ct1%qnten + ct%qnten * dpc%rat_dp(i)
       ct1%uten  = ct1%uten  + ct%uten  * dpc%rat_dp(i)
       ct1%vten  = ct1%vten  + ct%vten  * dpc%rat_dp(i)
       ct1%pflx  = ct1%pflx  + ct%pflx  * dpc%rat_dp(i)
       ct1%hlflx = ct1%hlflx + ct%hlflx * dpc%rat_dp(i)
       ct1%qctflx= ct1%qctflx+ ct%qctflx* dpc%rat_dp(i)
!       ct1%tevap = ct1%tevap + ct%tevap * dpc%rat_dp(i)
!       ct1%qevap = ct1%qevap + ct%qevap * dpc%rat_dp(i)
       ct1%rain  = ct1%rain  + ct%rain  * dpc%rat_dp(i)
       ct1%snow  = ct1%snow  + ct%snow  * dpc%rat_dp(i)
       ct1%denth = ct1%denth + ct%denth * dpc%rat_dp(i)

       cp1%ufrc  = cp1%ufrc  + cp%ufrc  * dpc%rat_dp(i)
       cp1%qlu   = cp1%qlu   + cp%qlu   * dpc%rat_dp(i)
       cp1%qiu   = cp1%qiu   + cp%qiu   * dpc%rat_dp(i)
       cp1%qnu   = cp1%qnu   + cp%qnu   * dpc%rat_dp(i)
       cp1%umf   = cp1%umf   + cp%umf   * dpc%rat_dp(i)
       cp1%wu    = cp1%wu    + cp%wu    * dpc%rat_dp(i)
       cp1%fdrsat= cp1%fdrsat+ cp%fdrsat* dpc%rat_dp(i)
       cp1%fdr   = cp1%fdr   + cp%fdr   * dpc%rat_dp(i)
    end do

    call sd_copy_k(sd, sd1)
    tmp      = cbmf0 / cbmf
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    dcape=(ac%cape-ac1%cape)/cbmf0
    if (dcape <= 0.) then
       dpc%cbmf_d=0.; ocode=6; return
    end if

    if (dpc%ideep_closure.eq.1) then
       cbmfs = (ac%cape - dpc%cape_th) / dcape / (dpc%tau_dp/sd%delt)
    else if (dpc%ideep_closure.eq.2) then
       cbmfs = cc%dcape / dcape 
    else
       cbmfs = 0.0
    end if
!!$       call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, cc%wrel, cc%scaleh)
!!$       cc%dwfn=0.; cc%wfn=0.; delp=0.;
!!$       do k=cp1%krel, cp1%let
!!$          cc % wfn  = cc % wfn  + 0.5*(cp %wu(k)*cp %wu(k)) * cp%dp(k)
!!$          cc % dwfn = cc % dwfn + 0.5*(cp1%wu(k)*cp1%wu(k) - cp%wu(k)*cp%wu(k)) * cp%dp(k)
!!$          delp      = delp + cp%dp(k)
!!$       end do
!!$       cc % wfn  = cc % wfn  / delp 
!!$       cc % dwfn = cc % dwfn / delp / cbmf0
!!$       if (do_cape_closure) then
!!$          cbmfs = - ac%cape / cc % dcape / (dpc%tau_dp/sd%delt)
!!$        elseif (do_relaxwfn) then
!!$          cbmfs = - cc%wfn  / cc % dwfn  / (dpc%tau_dp/sd%delt)
!!$       else
!!$          cbmfs = - cc%wfn  / cc % dwfn
!!$          tmp   = sd%delt/dpc%tau_dp
!!$          cbmfs = (cbmf_old+tmp*cbmfs)/(1.+tmp)
!!$       end if

    cbmf_max=(sd%ps(0) - sd%ps(cp%krel))*(0.25/sd%delt)/Grav
    dpc%cbmf_d = max(min(cbmfs, cbmf_max), 0.)
 
    if(dpc%cbmf_d.lt.1.e-10) then 
       dpc%cbmf_d=0.; ocode=6; return
    end if

    tmp       = dpc%cbmf_d/ cbmf
    ct1%tten  = ct1%tten  * tmp
    ct1%qvten = ct1%qvten * tmp
    ct1%qlten = ct1%qlten * tmp
    ct1%qiten = ct1%qiten * tmp
    ct1%qaten = ct1%qaten * tmp
    ct1%qnten = ct1%qnten * tmp
    ct1%uten  = ct1%uten  * tmp
    ct1%vten  = ct1%vten  * tmp
    ct1%pflx  = ct1%pflx  * tmp
    ct1%hlflx = ct1%hlflx * tmp
    ct1%qctflx= ct1%qctflx* tmp
!    ct1%tevap = ct1%tevap * tmp
!    ct1%qevap = ct1%qevap * tmp
    ct1%rain  = ct1%rain  * tmp
    ct1%snow  = ct1%snow  * tmp
    ct1%denth = ct1%denth * tmp

!    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, dpc%cbmf_d, cc%wrel, 10000.)
!    call cumulus_tend_k(dpn, sd, cp1, ct1, do_coldT)

  end subroutine dpconv1


!#####################################################################
!#####################################################################


  subroutine dpconv2(dpc, cpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, sd1, ac1, &
       cc1, cp1, ct1, cbmf_deep, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(in)     :: cpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(inout)  :: ac
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cclosure),  intent(inout)  :: cc,cc1
    type(cplume),    intent(inout)  :: cp,cp1
    type(ctend),     intent(inout)  :: ct,ct1
    real,            intent(inout)  :: ocode, cbmf_deep
    integer,            intent(out)   :: ier
    character(len=256), intent(out)   :: ermesg

    integer :: k, ksrc, n
    real    :: cbmf0, cbmfs, cbmf_max, dcape, scaleh, wrel, tmp, rkm
    real    :: cwfn_d, dcwfn, rat, ratsum, zcldtop
    type(cpnlist) :: dpn

    call cpn_copy(cpn, dpn)
    dpn % do_ppen   = dpc % do_ppen
    dpn % do_pevap  = dpc % do_pevap
    dpn % mixing_assumption = dpc % mixing_assumption

    if (ocode.ne.0) then
       return
    end if

    call sd_copy_k(sd, sd1)
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    dpc%cbmf_d  = 0.; 
    dpc%cwfn_d  = 0.;
    dpc%rat     = 0.; 
    dpc%rkm     = 0.; 
    dpc%hgt     = 0.; 
    scaleh      = 0.;
    ratsum      = 0.;
    cbmf0       = 0.0001; 
    wrel        = cc%wrel; 
    rkm         = dpc%rkm_dp(1);
    zcldtop     = sd%z(cp%ltop); 

    do n=1, sd%kmax
       scaleh        = zcldtop
       dpc%hgt(n)    = scaleh
       dpc%rkm(n)    = rkm/scaleh
       rat           = dpc%rkm(n)
       dpc%rat(n)    = rat
       
       if (dpc%do_generation) then
          do k=1,cp%ltop
             sd1 % qv(k) = sd1 % qs(k)
          end do
          call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
          sd1 % thvtop(:) = sd % thvtop(:)
          sd1 % thvbot(:) = sd % thvbot(:)
       end if

       call cumulus_plume_k(dpn, sd1, ac, cp, rkm, cbmf0, wrel, scaleh, Uw_p, ier, ermesg)
       if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
          dpc%cbmf_d=0.; ocode=6; return
       end if

       zcldtop   = sd%z(cp%ltop)
       if (zcldtop.le.scaleh) then
          exit
       end if

       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct, do_coldT)

       ct1%tten  = ct1%tten  + ct%tten  * rat
       ct1%qvten = ct1%qvten + ct%qvten * rat
       ct1%qlten = ct1%qlten + ct%qlten * rat
       ct1%qiten = ct1%qiten + ct%qiten * rat
       ct1%qaten = ct1%qaten + ct%qaten * rat
       ct1%qnten = ct1%qnten + ct%qnten * rat
       ct1%uten  = ct1%uten  + ct%uten  * rat
       ct1%vten  = ct1%vten  + ct%vten  * rat
       ct1%tevap = ct1%tevap + ct%tevap * rat
       ct1%qevap = ct1%qevap + ct%qevap * rat
       ct1%pflx  = ct1%pflx  + ct%pflx  * rat
       ct1%hlflx = ct1%hlflx + ct%hlflx * rat
       ct1%qctflx= ct1%qctflx+ ct%qctflx* rat
       ct1%rain  = ct1%rain  + ct%rain  * rat
       ct1%snow  = ct1%snow  + ct%snow  * rat
       ct1%denth = ct1%denth + ct%denth * rat

       cp1%ufrc  = cp1%ufrc  + cp%ufrc  * rat
       cp1%qlu   = cp1%qlu   + cp%qlu   * rat
       cp1%qiu   = cp1%qiu   + cp%qiu   * rat
       cp1%qnu   = cp1%qnu   + cp%qnu   * rat
       cp1%umf   = cp1%umf   + cp%umf   * rat
       cp1%wu    = cp1%wu    + cp%wu    * rat
       cp1%fdrsat= cp1%fdrsat+ cp%fdrsat* rat
       cp1%fdr   = cp1%fdr   + cp%fdr   * rat

       tmp       = maxval(cp%wu(:))
       dpc%cwfn_d= max(dpc%cwfn_d, tmp*tmp)

       ratsum    = ratsum + rat
    end do

    if (n > 1) then
       rat       = 1./ratsum
       dpc%rat   = dpc%rat*rat
       ct1%tten  = ct1%tten  * rat
       ct1%qvten = ct1%qvten * rat
       ct1%qlten = ct1%qlten * rat
       ct1%qiten = ct1%qiten * rat
       ct1%qaten = ct1%qaten * rat
       ct1%qnten = ct1%qnten * rat
       ct1%uten  = ct1%uten  * rat
       ct1%vten  = ct1%vten  * rat
       ct1%tevap = ct1%tevap * rat
       ct1%qevap = ct1%qevap * rat
       ct1%pflx  = ct1%pflx  * rat
       ct1%hlflx = ct1%hlflx * rat
       ct1%qctflx= ct1%qctflx* rat
       ct1%rain  = ct1%rain  * rat
       ct1%snow  = ct1%snow  * rat
       ct1%denth = ct1%denth * rat
       
       cp1%ufrc  = cp1%ufrc  * rat
       cp1%qlu   = cp1%qlu   * rat
       cp1%qiu   = cp1%qiu   * rat
       cp1%qnu   = cp1%qnu   * rat
       cp1%umf   = cp1%umf   * rat
       cp1%wu    = cp1%wu    * rat
       cp1%fdrsat= cp1%fdrsat* rat
       cp1%fdr   = cp1%fdr   * rat       
    end if


    if (dpc%ideep_closure.eq.0) then
       cbmfs = cbmf_deep
    else 
       call sd_copy_k(sd, sd1)
       sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
       sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
       sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
       sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
       sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
       sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
       sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
       sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
       call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
       ksrc=1
       call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
            sd1, Uw_p, .false., do_ice, ac1)
       if (dpc%ideep_closure.eq.1) then
          dcape = (ac%cape - ac1%cape)/cbmf0
          if (dcape <= 0.) then
             dpc%cbmf_d=0.; ocode=6; return
          end if
          cbmfs = (ac%cape - dpc%cape_th) / dcape / (dpc%tau_dp/sd%delt)
       else if (dpc%ideep_closure.eq.2) then
          call cumulus_plume_k(dpn, sd1, ac1, cp, rkm, cbmf0, wrel, scaleh, Uw_p, ier, ermesg)
          tmp    = maxval(cp%wu(:))
          cwfn_d = tmp*tmp
          dcwfn=(dpc%cwfn_d - cwfn_d)  /cbmf0
          cbmfs = (dpc%cwfn_d - 0.) / dcwfn / (dpc%tau_dp/sd%delt)
       end if
    end if

    cbmf_max=(sd%ps(0) - sd%ps(cp%krel))*(0.25/sd%delt)/Grav
    dpc%cbmf_d = max(min(cbmfs, cbmf_max), 0.)
 
    if(dpc%cbmf_d.lt.1.e-10) then 
       dpc%cbmf_d=0.; ocode=6; return
    end if

    tmp       = dpc%cbmf_d/ cbmf0
    ct1%tten  = ct1%tten  * tmp
    ct1%qvten = ct1%qvten * tmp
    ct1%qlten = ct1%qlten * tmp
    ct1%qiten = ct1%qiten * tmp
    ct1%qaten = ct1%qaten * tmp
    ct1%qnten = ct1%qnten * tmp
    ct1%uten  = ct1%uten  * tmp
    ct1%vten  = ct1%vten  * tmp
    ct1%tevap = ct1%tevap * tmp
    ct1%qevap = ct1%qevap * tmp
    ct1%pflx  = ct1%pflx  * tmp
    ct1%hlflx = ct1%hlflx * tmp
    ct1%qctflx= ct1%qctflx* tmp
    ct1%rain  = ct1%rain  * tmp
    ct1%snow  = ct1%snow  * tmp
    ct1%denth = ct1%denth * tmp

    cp1%ufrc  = cp1%ufrc  * tmp
    cp1%umf   = cp1%umf   * tmp

  end subroutine dpconv2

!#####################################################################
!#####################################################################

  subroutine dpconv3(dpc, cpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
       omeg_avg, rkm_sh, sd1, ac1, cp1, ct1, cbmf_deep, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(in)     :: cpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    real,            intent(in)     :: rkm_sh, omeg_avg
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp, cp1
    type(ctend),     intent(inout)  :: ct, ct1
    real,            intent(inout)  :: cbmf_deep, ocode
    integer,            intent(out)   :: ier
    character(len=256), intent(out)   :: ermesg

    type(cpnlist) :: dpn
    real          :: rkm_dp, zcldtop, cbmf0, dcapedm, cbmf_max, tmp
    integer       :: ksrc

    zcldtop = sd%z(cp%ltop)
    if ( (ocode.ne.0) .or. (ac%cape  <= dpc%cape_th) .or. zcldtop < dpc%deepdepth) then
       ocode=6; return
    end if

    call cpn_copy(cpn, dpn)
    dpn % do_ppen   = dpc % do_ppen
    dpn % do_pevap  = dpc % do_pevap
    rkm_dp  = dpc%rkm_dp(1) *  rkm_sh

    cbmf0 = 0.0001
    call cp_clear_k(cp1);
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf0, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=6; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

    call sd_copy_k(sd, sd1)
    tmp      = 1.
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    ksrc=1
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    dcapedm=(ac%cape-ac1%cape)/cbmf0

    if (dcapedm <= 0.) then
       cbmf_deep=0.; ocode=6; return
    else
       cbmf_deep= (ac%cape - dpc%cape_th) / dcapedm / (dpc%tau_dp/sd%delt)
    end if

    cbmf_max=(sd%ps(0) - sd%ps(cp%krel))*(0.25/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)
 
    if(cbmf_deep.lt.1.e-10) then 
       cbmf_deep=0.; ocode=6; return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=6; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
    
  end subroutine dpconv3


!#####################################################################
!#####################################################################

subroutine DEEP_CONV_INIT

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.
    


end subroutine DEEP_CONV_INIT
!#####################################################################
!#####################################################################
subroutine DEEP_CONV_END

      module_is_initialized = .false.


end subroutine DEEP_CONV_END

!#####################################################################
!#####################################################################

end MODULE DEEP_CONV_MOD
