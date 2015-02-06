
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

  character(len=128) :: version = '$Id: deep_conv.F90,v 20.0 2013/12/13 23:21:39 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

!-------  interfaces --------

  public  :: cpn_copy, dpconv0, dpconv1, dpconv2

  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'deep_conv'

  public deepc
  type deepc
     real    :: rkm_dp1
     real    :: rkm_dp2
     real    :: cbmf_dp_frac1
     real    :: cbmf_dp_frac2
     real    :: crh_th_land
     real    :: crh_th_ocean
     real    :: cape_th 
     real    :: tau_dp 
     real    :: rpen_d
     integer :: mixing_assumption_d
     logical :: do_ppen_d
     logical :: do_pevap_d
     real    :: cfrac_d
     real    :: hcevap_d
     real    :: dcapedm_th
     real    :: frac_limit_d
     real    :: lofactor_d
     real    :: tcrit_d
     real    :: auto_th0_d
     logical :: do_forcedlifting_d
  end type deepc

contains

!#####################################################################
!#####################################################################

  subroutine cpn_copy(cpn, dpn)
    type(cpnlist), intent(in)    :: cpn
    type(cpnlist), intent(inout) :: dpn

    dpn % do_qctflx_zero     = cpn % do_qctflx_zero
    dpn % do_detran_zero     = cpn % do_detran_zero
    dpn % rle                = cpn % rle
    dpn % rpen               = cpn % rpen
    dpn % rmaxfrac           = cpn % rmaxfrac
    dpn % wmin               = cpn % wmin
    dpn % rbuoy              = cpn % rbuoy
    dpn % rdrag              = cpn % rdrag  
    dpn % frac_drs           = cpn % frac_drs
    dpn % bigc               = cpn % bigc    
    dpn % auto_th0           = cpn % auto_th0
    dpn % deltaqc0           = cpn % deltaqc0
    dpn % do_pdfpcp          = cpn % do_pdfpcp
    dpn % do_pmadjt          = cpn % do_pmadjt
    dpn % do_emmax           = cpn % do_emmax
    dpn % do_pnqv            = cpn % do_pnqv
    dpn % emfrac_max         = cpn % emfrac_max
    dpn % auto_rate          = cpn % auto_rate
    dpn % tcrit              = cpn % tcrit  
    dpn % cldhgt_max         = cpn % cldhgt_max
    dpn % do_ice             = cpn % do_ice
    dpn % do_ppen            = cpn % do_ppen
    dpn % do_pevap           = cpn % do_pevap
    dpn % hcevap             = cpn % hcevap
    dpn % cfrac              = cpn % cfrac
    dpn % mixing_assumption  = cpn % mixing_assumption
    dpn % mp_choice          = cpn % mp_choice
    dpn % Nl_land            = cpn % Nl_land
    dpn % Nl_ocean           = cpn % Nl_ocean
    dpn % qi_thresh          = cpn % qi_thresh
    dpn % r_thresh           = cpn % r_thresh
    dpn % peff_l             = cpn % peff_l
    dpn % peff_i             = cpn % peff_i
    dpn % peff               = cpn % peff
    dpn % t00                = cpn % t00
    dpn % rh0                = cpn % rh0
    dpn % do_forcedlifting   = cpn % do_forcedlifting
    dpn % atopevap           = cpn % atopevap
    dpn % wtwmin_ratio       = cpn % wtwmin_ratio
    dpn % do_auto_aero       = cpn % do_auto_aero
    dpn % rad_crit           = cpn % rad_crit
    dpn % wrel_min           = cpn % wrel_min
    dpn % do_weffect         = cpn % do_weffect
    dpn % weffect            = cpn % weffect
    dpn % use_online_aerosol = cpn % use_online_aerosol
    dpn % isdeep             = cpn % isdeep

  end subroutine cpn_copy

!#####################################################################
!#####################################################################

  subroutine dpconv0(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, cp1, ct1, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, dcrh, cbmf_dp_frac

    ier = 0
    ermesg = ' '
    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
       ocode=6;
       return
    end if

    zcldtop = sd%z(cp%ltop)
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
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

  subroutine dpconv1(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, dcapedm, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, dcapedm
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, wrel, cbmf0, cbmf_max, tmp
    integer       :: ksrc

    ier = 0
    ermesg = ' '
    zcldtop = 2000 !sd%z(cp%ltop)
    wrel = max(cc%wrel, 0.1)

!    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
    if ( cbmf_deep.eq.0 ) then
       ocode=6;
       return
    end if
    if (ac%cape .lt. dpc%cape_th) then
       ocode=7; return
    end if

    cbmf0 = 0.0001
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=8; return
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

    call ac_clear_k(ac1);
    ksrc=1
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)

    dcapedm=(ac%cape-ac1%cape)/cbmf0

    if (dcapedm .lt. dpc%dcapedm_th) then
       cbmf_deep=0.; ocode=9; 
       call ct_clear_k(ct1);
       return
    else
       cbmf_deep= (ac%cape - dpc%cape_th) / dcapedm / (dpc%tau_dp/sd%delt)
    end if

    cbmf_max  = (sd%ps(0) - sd%ps(cp1%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)
 
    if(cbmf_deep.lt.1.e-10) then 
       cbmf_deep=0.; ocode=10; 
       call ct_clear_k(ct1);
       return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=11; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
    
  end subroutine dpconv1

!#####################################################################
!#####################################################################

  subroutine dpconv2(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, cbmf0, dcwfn, cwfn0, cwfn1, dpsum, cbmf_max, tmp
    integer       :: ksrc, k

    ier = 0
    ermesg = ' '
    zcldtop = sd%z(cp%ltop)

    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
       ocode=6;
       return
    end if

    cbmf0 = 0.0001
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf0, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=8; return
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

    call ac_clear_k(ac1);
    ksrc=1
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, cc%wrel, zcldtop, Uw_p, ier, ermesg)

    cwfn0=0.; dpsum=0.
    do k = cp%let,cp%ltop
       cwfn0 = cwfn0 + cp%buo(k)*sd%dp(k)
       dpsum = dpsum + sd%dp(k)
    end do
    cwfn0 = cwfn0 /dpsum

    cwfn1=0.; dpsum=0.
    do k = cp1%let,cp1%ltop
       cwfn1 = cwfn1 + cp%buo(k)*sd%dp(k)
       dpsum = dpsum + sd%dp(k)
    end do
    cwfn1 = cwfn1 /dpsum

    dcwfn =(cwfn0 - cwfn1)/cbmf0

    if (dcwfn <= 0.) then
       cbmf_deep=0.; ocode=9; return
    else
       cbmf_deep= (cwfn0-0) / dcwfn / (dpc%tau_dp/sd%delt)
    end if

    cbmf_max  = (sd%ps(0) - sd%ps(cp%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)
 
    if(cbmf_deep.lt.1.e-10) then 
       cbmf_deep=0.; ocode=10; 
       call ct_clear_k(ct1);
       return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=11; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
    
  end subroutine dpconv2
!#####################################################################
!#####################################################################

end MODULE DEEP_CONV_MOD
