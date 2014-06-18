
MODULE CONV_CLOSURES_MOD

! use Sat_Vapor_Pres_Mod, ONLY: ESCOMP, DESCOMP
! use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use fms_mod,              only: write_version_number
  use  conv_utilities_k_mod,only: sd_copy_k, adi_cloud_k, extend_sd_k,&
                                  adicloud, sounding, uw_params
  use  conv_plumes_k_mod,   only: cumulus_plume_k, cumulus_tend_k, &
                                  cplume, ctend, cpnlist

!---------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: conv_closures.F90,v 19.0 2012/01/06 20:25:26 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'
  logical            :: module_is_initialized=.false.  ! module initialized ?

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: cclosure_bretherton, cclosure_relaxcbmf, cclosure_emanuel, &
             cclosure_implicit, cclosure_relaxwfn, &
             conv_closures_init, conv_closures_end

  character(len=11) :: mod_name = 'conv_closures'

  public cclosure
  type cclosure
     real    :: cbmf, wrel, ufrc, scaleh, dcin, dcape, dwfn, wfn
     integer :: igauss
     real    :: rkfre, rmaxfrac, rbuoy, tau_sh, tau_dp, wcrit_min
     real    :: maxcldfrac
  end type cclosure
 
contains

!#####################################################################
!#####################################################################

  function erfccc(x)
    !--------------------------------------------------------------
    ! This numerical recipes routine calculates the complementary
    ! error function.
    !--------------------------------------------------------------
   
    real :: erfccc
    real, intent(in) :: x 
    real :: t,z
   
    z=abs(x)      
    t=1./(1.+0.5*z)
   
    erfccc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*      &
         (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*    &
         (1.48851587+t*(-.82215223+t*.17087277)))))))))
    
    if (x.lt.0.) erfccc=2.-erfccc
    
  end function erfccc

!#####################################################################
!#####################################################################

  subroutine solvecbmf(alpha, beta, x)
    !Newton iteration solving Eq. x = beta * exp (-alpha * x)
    implicit none
    real,    intent(in)    :: alpha, beta
    real,    intent(inout) :: x
    
    integer :: iteration, niteration=5, id_check
    real    :: dydx, x0, y0

    x0=1.
    do iteration = 1,niteration
       y0   = x0 - beta * exp(-alpha * x0)
       dydx = 1. + beta * alpha * exp(-alpha * x0)
       x0   = x0 - y0 / dydx
       if (abs(y0) < 0.0001) then
          x=x0; id_check=0
       else
          id_check=1
       end if
    end do
    if (id_check==1) then
      x=1.
      print*, 'ID_CHECK=1, in solvecbmfffffffffffff'
    endif
    
  end subroutine solvecbmf

!#####################################################################
!#####################################################################

  subroutine cclosure_bretherton(tkeavg, cpn, sd, Uw_p,  ac, cc, &
                                 cbmf_unmod)
  
    implicit none
    real,           intent(in)    :: tkeavg
    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(uw_params), intent(inout)    :: Uw_p
    type(adicloud), intent(in)    :: ac
    type(cclosure), intent(inout) :: cc
    real,   intent(out), optional :: cbmf_unmod
    
    real    :: sigmaw, wcrit, erfarg, cbmf, wexp, ufrc, wtw
    real    :: rmfk1=0.3, rmfk2=5.0, rmfk3=3.0

    cc%cbmf=0.; cc%wrel=0.; cc%ufrc=0.;
    if(cc%igauss.eq.0)then     !Use cin and pbl tke
       cbmf = rmfk1* ac % rho0lcl * sqrt(tkeavg) * exp(-rmfk2* ac % cin/tkeavg)
       wexp = rmfk3* sqrt(tkeavg) !Updraft vertical velocity at release height depends on tke

    elseif(cc%igauss.eq.1)then !Use cin and gaussian distribution of w
       wcrit  = sqrt(2. * ac % cin * cc%rbuoy)
       sigmaw = sqrt(cc%rkfre * tkeavg)
       wcrit = max(wcrit, cc%wcrit_min*sigmaw)
       cbmf   = ac % rho0lcl * sigmaw / 2.5066 * exp(-0.5*((wcrit/sigmaw)**2.))

      if (present (cbmf_unmod)) then
        cbmf_unmod = MAX(0.0, cbmf)
      endif
  
       !Diagnose updraft fraction sqrt(2.) = 1.4142
       erfarg=wcrit / (1.4142 * sigmaw)
       if(erfarg.lt.20.)then
          ufrc = min(cc%maxcldfrac, cc%rmaxfrac, 0.5*erfccc(erfarg))
       else
          ufrc = 0.
       endif

       if(ufrc.gt.0.0) then !Diagnose expected value of cloud base vertical velocity
           wexp = cbmf / ac % rho0lcl / ufrc
       else
          wexp = 0.
          cbmf = 0.
       endif
    endif

    wtw = wexp * wexp - 2 * ac % cin * cc%rbuoy !used for the runs of xx-hv1_amip and tropical storm 
    if(wtw.le.0.) then
       cc%wrel=0.; 
    else
       cc%wrel=sqrt(wtw)
    end if

    cc%cbmf=cbmf
    cc%wrel=min(cc%wrel, 50.)!cc%ufrc=min(cc%rmaxfrac, cc%ufrc)
    cbmf = (sd%ps(0) - ac%plcl ) * 0.25 / sd%delt / Uw_p%GRAV
    if (cc%cbmf .gt. cbmf) cc%cbmf = cbmf
    if (cc%wrel .gt. 0.) then
      cc%ufrc=cc%cbmf / cc%wrel /ac % rho0lcl
    else
      cc%ufrc=0.
    end if
    if (cc%ufrc > cc%maxcldfrac) then
       cc%ufrc = cc%maxcldfrac
       cc%cbmf = cc%wrel*ac%rho0lcl*cc%ufrc
    end if   

!    cc%cbmf=cbmf
!    cc%ufrc=ufrc
!
!    cbmf = (sd%ps(0) - ac%plcl ) * 0.25 / sd%delt / Uw_p%GRAV
!    if (cc%cbmf .gt. cbmf .and. cc%wrel .gt. 0) then
!       cc%cbmf = cbmf
!       cc%ufrc = cc%cbmf / wexp /ac % rho0lcl
!    end if
!    cc%wrel=min(cc%wrel, 50.)
!    cc%ufrc=min(cc%rmaxfrac, cc%ufrc)

    return

  end subroutine cclosure_bretherton

!#####################################################################
!#####################################################################

  subroutine cclosure_implicit(tkeavg, cpn, sd, Uw_p, ac, cc, delt, rkm, &
       do_coldT, sd1, ac1, cc1, cp1, ct1, ier, ermesg)
    implicit none
    real,           intent(in)    :: tkeavg, delt, rkm
    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(uw_params), intent(inout)    :: Uw_p
    type(adicloud), intent(in)    :: ac
    type(cclosure), intent(inout) :: cc, cc1
    type(sounding), intent(inout) :: sd1
    type(adicloud), intent(inout) :: ac1
    type(cplume),   intent(inout) :: cp1
    type(ctend),    intent(inout) :: ct1
    logical,        intent(in)    :: do_coldT
    integer,        intent(out)     :: ier
    character(len=256), intent(out) :: ermesg

    logical :: dofast=.false., doice=.true.

    real :: cbmf0=0.001, alpha, beta, phi

    call cclosure_bretherton(tkeavg, cpn, sd, Uw_p, ac, cc)
    if(cc%cbmf.eq.0.) then
      cc % dcin=0.
      return
    end if

    call cumulus_plume_k(cpn, sd, ac, cp1, rkm, cbmf0, cc%wrel, cc%scaleh, Uw_p, ier, ermesg)
    if (ier /= 0) then
      ermesg = 'Called from cclosure_implicit : '// trim(ermesg)
      return
    endif
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cc % dcin=0.
       return
    else
       call cumulus_tend_k(cpn, sd, Uw_p, cp1, ct1, do_coldT)
       call sd_copy_k(sd, sd1)
       sd1 % t  = sd1 % t  + ct1%tten  * delt
       sd1 % qv = sd1 % qv + ct1%qvten * delt
       sd1 % ql = sd1 % ql + ct1%qlten * delt
       sd1 % qi = sd1 % qi + ct1%qiten * delt
       sd1 % qa = sd1 % qa + ct1%qaten * delt
       sd1 % qn = sd1 % qn + ct1%qnten * delt
       sd1 % u  = sd1 % u  + ct1%uten  * delt
       sd1 % v  = sd1 % v  + ct1%vten  * delt

       call extend_sd_k(sd1, sd%pblht, doice, Uw_p)

       call adi_cloud_k(sd1%zs(1), sd1%ps(1), sd1%hl(1), sd1%thc(1), sd1%qct(1), sd1, Uw_p, dofast, doice, ac1)

       cc % dcin=(ac1%cin-ac%cin)/cbmf0

       alpha  = (2. * cc%rbuoy) / (2. * cc%rkfre * tkeavg) * cc % cbmf * cc % dcin
       beta   = 1.  ! ac % rho0lcl * sqrt(cc%rkfre * tkeavg) / 2.5066
       phi    = 1.
       if (alpha .gt. 0.) then
          call solvecbmf(alpha, beta, phi)
          cc % cbmf = phi * cc % cbmf
       end if
    end if

  end subroutine cclosure_implicit

!#####################################################################
!#####################################################################

  subroutine cclosure_relaxcbmf(tkeavg, cpn, sd, Uw_p, ac, cc, delt)
  
    implicit none
    real,           intent(in)    :: tkeavg, delt
    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(uw_params), intent(in)    :: Uw_p
    type(adicloud), intent(in)    :: ac
    type(cclosure), intent(inout) :: cc
    
    real    :: sigmaw, wcrit, erfarg, wexp, wtw
    real    :: cbmfs, tmp
    
    cc%wrel=0.; cc%ufrc=0.;   

    wcrit  = sqrt(2. * ac % cin * cc%rbuoy)
    sigmaw = sqrt(cc%rkfre * tkeavg)
    cbmfs  = ac % rho0lcl * sigmaw / 2.5066 * exp(-0.5*((wcrit/sigmaw)**2.))

    tmp    = delt/cc%tau_sh
    cc%cbmf= max((cc%cbmf+tmp*cbmfs)/(1.+tmp),0.0)

    !Diagnose updraft fraction
    erfarg=wcrit / (1.4142 * sigmaw)
    if(erfarg.lt.20.)then
       cc%ufrc = min(cc%maxcldfrac, cc%rmaxfrac, 0.5*erfccc(erfarg))
    else
       cc%ufrc = 0.
    endif

    if(cc%ufrc.gt.0.001)then !Diagnose expected value of cloud base vertical velocity
       wexp = cc%cbmf / ac % rho0lcl / cc%ufrc
    else
       wexp = 0.
       cc%cbmf = 0.
    endif

    wexp=min(wexp, 50.)
    wtw = wexp * wexp - 2 * ac % cin * cc%rbuoy
    if(wtw.le.0.) then
       cc%wrel=0.; 
    else
       cc%wrel=sqrt(wtw)
    end if

    return

  end subroutine cclosure_relaxcbmf

!#####################################################################
!#####################################################################

  subroutine cclosure_emanuel(tkeavg, cpn, sd, Uw_p, ac, cc, delt)
  
    implicit none
    real,           intent(in)    :: tkeavg, delt
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    type(sounding), intent(in)    :: sd
    type(adicloud), intent(in)    :: ac
    type(cclosure), intent(inout) :: cc
    
    integer :: k
    real    :: ufrc=0.01
    real    :: dtmin, dpsum, dtpbl, damps
    real    :: cbmf
    real    :: dtmax    = 0.9    ! MAXIMUM NEGATIVE TEMPERATURE PERTURBATION
                                 ! A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC
    real    :: damp     = 0.1    ! ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF
    real    :: alpha    = 0.1    ! APPROACH TO QUASI-EQUILIBRIUM

    cbmf=cc%cbmf; cc%cbmf=0.; cc%wrel=0.; cc%ufrc=0.;

    dpsum=0.; dtpbl=0.
    do k=1, ac%klcl-1
       dtpbl=dtpbl+(ac%thv(k)-sd%thv(k))*sd%exner(k)*sd%dp(k)
       dpsum=dpsum+sd%dp(k);
    end do
    dtpbl=dtpbl/dpsum
    dtmin=(ac%thvlcl-ac%thv0lcl)+dtpbl+dtmax

    damps=damp*delt/300.
    cbmf =(1.-damps)*cbmf+0.1*alpha*dtmin
    cc%cbmf=max(cbmf,0.0)
    cc%wrel=cbmf / ac % rho0lcl / ufrc
    cc%ufrc=ufrc

    return

  end subroutine cclosure_emanuel

!#####################################################################
!#####################################################################

  subroutine cclosure_relaxwfn(tkeavg, cpn, sd, Uw_p, ac, cc, cp, ct, delt, rkm, &
       do_coldT, sd1, ac1, cc1, cp1, ct1, ier, ermesg)
    implicit none
    real,           intent(in)    :: tkeavg, delt, rkm
    type(cpnlist),  intent(in)    :: cpn
    type(uw_params),  intent(inout)    :: Uw_p
    type(sounding), intent(in)    :: sd
    type(adicloud), intent(in)    :: ac
    type(cclosure), intent(inout) :: cc, cc1
    type(sounding), intent(inout) :: sd1
    type(adicloud), intent(inout) :: ac1
    type(cplume),   intent(inout) :: cp, cp1
    type(ctend),    intent(inout) :: ct, ct1
    logical,        intent(in)    :: do_coldT
    integer,        intent(out)   :: ier
    character(len=256), intent(out) :: ermesg
    logical :: dofast=.false., doice=.true.


    integer :: k
    real    :: cbmf0=0.0001, delp, cbmf_old, tmp, cbmfs

    cbmf_old= cc%cbmf
    call cclosure_bretherton(tkeavg, cpn, sd, Uw_p, ac, cc)

    call cumulus_plume_k(cpn, sd,  ac, cp, rkm, cbmf0, cc%wrel, cc%scaleh, Uw_p, ier, ermesg)
    if (ier /= 0) then
      ermesg = 'Called from cclosure_relaxwfn : '//trim(ermesg)
      return
    endif
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cc % dcin=0.
       return
    else
       call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)
       call sd_copy_k(sd, sd1)
       sd1 % t  = sd1 % t  + ct%tten  * delt
       sd1 % qv = sd1 % qv + ct%qvten * delt
       sd1 % ql = sd1 % ql + ct%qlten * delt
       sd1 % qi = sd1 % qi + ct%qiten * delt
       sd1 % qa = sd1 % qa + ct%qaten * delt
       sd1 % qn = sd1 % qn + ct%qnten * delt
       sd1 % u  = sd1 % u  + ct%uten  * delt
       sd1 % v  = sd1 % v  + ct%vten  * delt

       call extend_sd_k(sd1, sd%pblht, doice, Uw_p)

       call adi_cloud_k(sd1%zs(1), sd1%ps(1), sd1%hl(1), sd1%thc(1), sd1%qct(1), sd1, Uw_p, dofast, doice, ac1)
       cc % dcin=(ac1%cin-ac%cin)/cbmf0
       cc % dcape=(ac1%cape-ac%cape)/cbmf0

       call cumulus_plume_k(cpn, sd1, ac1, cp1, rkm, cbmf0, cc%wrel, cc%scaleh, Uw_p, ier, ermesg)
       if (ier /= 0) then
         ermesg = 'Called from cclosure_relaxwfn 2nd call : '//trim(ermesg)
       endif

       cc%dwfn=0.; cc%wfn=0.; delp=0.;
       do k=cp1%krel, cp1%let
          cc % wfn  = cc % wfn  + 0.5*(cp %wu(k)*cp %wu(k)) * cp%dp(k)
          cc % dwfn = cc % dwfn + 0.5*(cp1%wu(k)*cp1%wu(k) - cp%wu(k)*cp%wu(k)) * cp%dp(k)
          delp      = delp + cp%dp(k)
       end do
       cc % wfn  = cc % wfn  / delp 
       cc % dwfn = cc % dwfn / delp / cbmf0

       cbmfs = - cc%wfn  / cc % dwfn

       tmp    = delt/cc%tau_sh
       cc%cbmf= (cbmf_old+tmp*cbmfs)/(1.+tmp)

       cc % cbmf = max(cc%cbmf,0.)


    end if

  end subroutine cclosure_relaxwfn

!#####################################################################
!#####################################################################
subroutine conv_closures_init

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.
    


end subroutine conv_closures_init
!#####################################################################
!#####################################################################
subroutine conv_closures_end

      module_is_initialized = .false.


end subroutine conv_closures_end

!#####################################################################
!#####################################################################

end MODULE CONV_CLOSURES_MOD
