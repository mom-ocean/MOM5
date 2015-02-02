#include <fms_platform.h>

MODULE CONV_UTILITIES_k_MOD
  
  use Sat_Vapor_Pres_k_Mod, ONLY: compute_qs_k

!---------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: conv_utilities_k.F90,v 20.0 2013/12/13 23:21:37 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: sd_init_k, sd_copy_k, sd_end_k, ac_init_k, ac_clear_k,  &
             ac_end_k, qsat_k, qses_k, exn_k, exn_init_k, exn_end_k, &
             findt_k, findt_init_k, findt_end_k, uw_params_init_k, &
             conden_k, pack_sd_k, pack_sd_lsm_k, extend_sd_k,  &
             adi_cloud_k, check_tracer_realizability, qt_parcel_k





 public sounding
 type sounding
    logical  :: coldT
    integer  :: kmax, kinv, ktoppbl, ktopconv
    real     :: psfc, pinv, zinv, thvinv, land, pblht, qint, delt, crh, tke
    real, _ALLOCATABLE :: t     (:)_NULL, qv   (:)_NULL, u     (:)_NULL
    real, _ALLOCATABLE :: v     (:)_NULL, ql   (:)_NULL, qi    (:)_NULL
    real, _ALLOCATABLE :: qa    (:)_NULL, thc  (:)_NULL, qct   (:)_NULL
    real, _ALLOCATABLE :: thv   (:)_NULL, rh   (:)_NULL, p     (:)_NULL
    real, _ALLOCATABLE :: z     (:)_NULL, dp   (:)_NULL, dz    (:)_NULL
    real, _ALLOCATABLE :: rho   (:)_NULL, nu   (:)_NULL, leff  (:)_NULL
    real, _ALLOCATABLE :: exner (:)_NULL, ps   (:)_NULL, exners(:)_NULL
    real, _ALLOCATABLE :: zs    (:)_NULL, ssthc(:)_NULL, ssqct (:)_NULL
    real, _ALLOCATABLE :: dudp  (:)_NULL, dvdp (:)_NULL, thvbot(:)_NULL
    real, _ALLOCATABLE :: thvtop(:)_NULL, qn   (:)_NULL, qs    (:)_NULL
    real, _ALLOCATABLE :: am1   (:)_NULL, am2  (:)_NULL, am3   (:)_NULL
    real, _ALLOCATABLE :: am4   (:)_NULL
    real, _ALLOCATABLE :: hl    (:)_NULL, sshl (:)_NULL, hm    (:)_NULL
    real, _ALLOCATABLE :: hms   (:)_NULL
!++++yim     
    real, _ALLOCATABLE :: tr    (:,:)_NULL, sstr(:,:)_NULL
 end type sounding

 public adicloud
 type adicloud
   real     :: usrc, vsrc, hlsrc, thcsrc, qctsrc
    integer  :: klcl, klfc, klnb
    real     :: plcl, zlcl, thvlcl, thv0lcl, rho0lcl
    real     :: plfc, plnb, cape, cin
    real, _ALLOCATABLE :: t  (:)_NULL, qv  (:)_NULL, ql  (:)_NULL
    real, _ALLOCATABLE :: qi (:)_NULL, thc (:)_NULL, qct (:)_NULL
    real, _ALLOCATABLE :: thv(:)_NULL, nu  (:)_NULL, leff(:)_NULL
    real, _ALLOCATABLE :: hl (:)_NULL, buo (:)_NULL
 end type adicloud

 public uw_params
 type uw_params
   real  :: hlv, hls, hlf, cp_air, grav, kappa, rdgas, p00, epsilo,  &
            zvir, tkmin, tkmax
   integer :: me
   logical :: master
 end type uw_params

! Lookup table dimension and ranges for findt_k
! These values should not be changed without a careful evaluation
! of the accuracy implications (cjg).

    integer, parameter :: nta = 1700
    real, parameter :: tamin =  120.0
    real, parameter :: tamax =  700.0

    integer, parameter :: np1 = 300
    real, parameter :: p1min = 10.0e2
    real, parameter :: p1max = 100.0e2
    integer, parameter :: np2 = 1000
    real, parameter :: p2min = 100.0e2
    real, parameter :: p2max = 1100.0e2
 
    real dta, dp1, dp2
    real rdta, rdp1, rdp2
    real(kind=4), allocatable :: ta_lookup(:,:,:)
    logical :: ta_lookup_allocated = .false.
 
! Lookup table for exn_k function (cjg)

    integer, parameter :: npex = 100000
    real, parameter :: pexmin = 10.0e2
    real, parameter :: pexmax = 1100.0e2

    real dpex, rdpex
    real(kind=4), allocatable :: ex_lookup(:)
    logical :: ex_lookup_allocated = .false.

contains

!#####################################################################
!#####################################################################

  subroutine uw_params_init_k (hlv, hls, hlf, cp_air, grav, kappa,  &
                               rdgas, p00, epsilo, zvir, tkmin,  &
                               tkmax, me, root_pe, Uw_p)
    real, intent(in) :: hlv, hls, hlf, cp_air, grav, kappa, rdgas,  &
                        p00, epsilo, zvir, tkmin, tkmax
    integer, intent(in) :: me, root_pe
    type(uw_params), intent(inout) :: Uw_p

    Uw_p%hlv = hlv
    Uw_p%hls = hls
    Uw_p%hlf = hlf
    Uw_p%cp_air = cp_air
    Uw_p%grav   = grav  
    Uw_p%kappa  = kappa 
    Uw_p%rdgas  = rdgas 
    Uw_p%p00    = p00   
    Uw_p%epsilo = epsilo
    Uw_p%zvir   = zvir  
    Uw_p%tkmin  = tkmin
    Uw_p%tkmax  = tkmax 
    Uw_p%me     = me
    Uw_p%master = (me == root_pe)

  end subroutine uw_params_init_k

!#####################################################################
!#####################################################################

  subroutine sd_init_k(kd, num_tracers, sd)
    integer, intent(in) :: kd, num_tracers
    type(sounding), intent(inout) :: sd
    
    sd%coldT    = .false.
    sd%kmax     = kd
    sd%kinv     = 0
    sd%ktoppbl  = 0
    sd%ktopconv = 0
    sd%psfc     = 0.0
    sd%pinv     = 0.0
    sd%zinv     = 0.0
    sd%thvinv   = 0.0
    sd%land     = 0.0
    sd%pblht    = 0.0
    sd%qint     = 0.0
    sd%delt     = 0.0
    sd%crh      = 0.0
    sd%tke      = 0.0
    allocate ( sd%t     (1:kd)); sd%t     =0.;
    allocate ( sd%qv    (1:kd)); sd%qv    =0.;
    allocate ( sd%u     (1:kd)); sd%u     =0.;
    allocate ( sd%v     (1:kd)); sd%v     =0.;
    allocate ( sd%qs    (1:kd)); sd%qs    =0.;
    allocate ( sd%ql    (1:kd)); sd%ql    =0.;
    allocate ( sd%qi    (1:kd)); sd%qi    =0.;
    allocate ( sd%qa    (1:kd)); sd%qa    =0.;
    allocate ( sd%qn    (1:kd)); sd%qn    =0.;
    allocate ( sd%thc   (1:kd)); sd%thc   =0.;
    allocate ( sd%qct   (1:kd)); sd%qct   =0.;
    allocate ( sd%thv   (1:kd)); sd%thv   =0.;
    allocate ( sd%rh    (1:kd)); sd%rh    =0.;
    allocate ( sd%p     (1:kd)); sd%p     =0.;
    allocate ( sd%z     (1:kd)); sd%z     =0.;
    allocate ( sd%dp    (1:kd)); sd%dp    =0.;
    allocate ( sd%dz    (1:kd)); sd%dz    =0.;
    allocate ( sd%rho   (1:kd)); sd%rho   =0.;
    allocate ( sd%nu    (1:kd)); sd%nu    =0.;
    allocate ( sd%leff  (1:kd)); sd%leff  =0.;
    allocate ( sd%exner (1:kd)); sd%exner =0.;
    allocate ( sd%ps    (0:kd)); sd%ps    =0.;
    allocate ( sd%zs    (0:kd)); sd%zs    =0.;
    allocate ( sd%exners(0:kd)); sd%exners=0.;
    allocate ( sd%ssthc (1:kd)); sd%ssthc =0.;
    allocate ( sd%ssqct (1:kd)); sd%ssqct =0.;
    allocate ( sd%dudp  (1:kd)); sd%dudp  =0.;
    allocate ( sd%dvdp  (1:kd)); sd%dvdp  =0.;
    allocate ( sd%thvbot(1:kd)); sd%thvbot=0.;
    allocate ( sd%thvtop(1:kd)); sd%thvtop=0.;
    allocate ( sd%am1   (1:kd)); sd%am1   =0.;
    allocate ( sd%am2   (1:kd)); sd%am2   =0.;
    allocate ( sd%am3   (1:kd)); sd%am3   =0.;
    allocate ( sd%am4   (1:kd)); sd%am4   =0.;
    allocate ( sd%hl    (1:kd)); sd%hl    =0.;
    allocate ( sd%hm    (1:kd)); sd%hm    =0.;
    allocate ( sd%hms   (1:kd)); sd%hms   =0.;
    allocate ( sd%sshl  (1:kd)); sd%sshl  =0.;
!++++yim
    allocate ( sd%tr  (1:kd,1:num_tracers)); sd%tr  =0.;
    allocate ( sd%sstr  (1:kd,1:num_tracers)); sd%sstr  =0.;
    
  end subroutine sd_init_k

!#####################################################################
!#####################################################################

  subroutine sd_copy_k(sd, sd1)
    type(sounding), intent(in)    :: sd
    type(sounding), intent(inout) :: sd1
    
    sd1% ktopconv = sd % ktopconv
    sd1% kmax = sd % kmax
    sd1% land = sd % land
    sd1% coldT= sd % coldT
    sd1%p     = sd%p;    sd1%z     =sd%z;
    sd1%ps    = sd%ps;   sd1%zs    =sd%zs;
    sd1%t     = sd%t;    sd1%qv    =sd%qv;
    sd1%u     = sd%u;    sd1%v     =sd%v;
    sd1%ql    = sd%ql;   sd1%qi    =sd%qi;
    sd1%qa    = sd%qa;   sd1%qn    =sd%qn;    
    sd1%am1   = sd%am1;  sd1%am2   =sd%am2; 
    sd1%am3   = sd%am3;  sd1%am4   =sd%am4;
    sd1%hl    = sd%hl;   sd1%sshl  =sd%sshl;
    sd1%hm    = sd%hm;   sd1%hms   =sd%hms;
!++++yim
    sd1%tr  =sd%tr
  end subroutine sd_copy_k

!#####################################################################
!#####################################################################

  subroutine sd_end_k(sd)
    type(sounding), intent(inout) :: sd
    deallocate ( sd%t, sd%qv, sd%u, sd%v, sd%ql, sd%qi, sd%qa, sd%thc, sd%qct,&
         sd%thv, sd%rh, sd%p, sd%z, sd%dp, sd%dz, sd%rho, sd%nu, sd%leff,     &
         sd%exner, sd%ps, sd%exners, sd%zs, sd%ssthc, sd%ssqct, sd%dudp,      &
         sd%dvdp, sd%thvbot, sd%thvtop, sd%qn, sd%am1, sd%am2, sd%am3, sd%am4,&
         sd%qs, sd%hl, sd%hm, sd%hms, sd%sshl, sd%tr, sd%sstr)
  end subroutine sd_end_k

!#####################################################################
!#####################################################################

  subroutine ac_init_k(kd, ac)
    integer, intent(in) :: kd
    type(adicloud), intent(inout) :: ac
    
    ac%usrc    = 0.0
    ac%vsrc    = 0.0
    ac%hlsrc   = 0.0
    ac%thcsrc  = 0.0
    ac%qctsrc  = 0.0
    ac%klcl    = 0
    ac%klfc    = 0
    ac%klnb    = 0
    ac%plcl    = 0.0
    ac%zlcl    = 0.0
    ac%thvlcl  = 0.0
    ac%thv0lcl = 0.0
    ac%rho0lcl = 0.0
    ac%plfc    = 0.0
    ac%plnb    = 0.0
    ac%cape    = 0.0
    ac%cin     = 0.0
    allocate ( ac%t     (1:kd)); ac%t    =0.;
    allocate ( ac%qv    (1:kd)); ac%qv   =0.;
    allocate ( ac%ql    (1:kd)); ac%ql   =0.;
    allocate ( ac%qi    (1:kd)); ac%qi   =0.;
    allocate ( ac%thc   (1:kd)); ac%thc  =0.;
    allocate ( ac%qct   (1:kd)); ac%qct  =0.;
    allocate ( ac%thv   (1:kd)); ac%thv  =0.;
    allocate ( ac%nu    (1:kd)); ac%nu   =0.;
    allocate ( ac%leff  (1:kd)); ac%leff =0.;
    allocate ( ac%hl    (1:kd)); ac%hl   =0.;
    allocate ( ac%buo   (1:kd)); ac%buo  =0.;
  end subroutine ac_init_k

!#####################################################################
!#####################################################################

  subroutine ac_clear_k(ac)
    type(adicloud), intent(inout) :: ac
    ac%t    =0.;    ac%qv   =0.;    ac%ql   =0.;
    ac%qi   =0.;    ac%thc  =0.;    ac%qct  =0.;
    ac%thv  =0.;    ac%nu   =0.;    ac%leff =0.; ac%hl   =0.;
    ac%buo  =0.;
  end subroutine ac_clear_k

!#####################################################################
!#####################################################################

  subroutine ac_end_k(ac)
    type(adicloud), intent(inout) :: ac
    deallocate (ac%t, ac%qv, ac%ql, ac%qi, ac%thc, ac%qct,  &
                ac%thv, ac%nu, ac%leff, ac%hl, ac%buo )
  end subroutine ac_end_k

!#####################################################################
!#####################################################################

  subroutine pack_sd_k (land, coldT, delt, pmid, pint, zmid, zint, &
                      u, v, t, qv, ql, qi, qa, qn, am1, am2, am3, am4,&
                        tracers, sd, Uw_p)

    real,    intent(in)              :: land
    logical, intent(in)              :: coldT
    real,    intent(in)              :: delt
    real, intent(in), dimension(:)   :: pmid, zmid !pressure&height@mid level
    real, intent(in), dimension(:)   :: pint, zint !pressure&height@ interface level
    real, intent(in), dimension(:)   :: u, v       !wind profile (m/s)
    real, intent(in), dimension(:)   :: t, qv      !temperature and specific humidity
    real, intent(in), dimension(:)   :: ql, qi, qa, qn !cloud tracers
    real, intent(in), dimension(:)   :: am1, am2, am3, am4  ! aerosal species
    real, intent(in), dimension(:,:) :: tracers        !env. tracers    
    type(sounding), intent(inout)    :: sd
    type(uw_params), intent(inout)    :: Uw_p

!++++yim
    integer :: k, nk, m
    real, parameter :: ptopconv = 3000.

    !Pack environmental sounding; layers are numbered from bottom up!=
    sd % kmax   = size(t)
    sd % land   = land
    sd % coldT  = coldT
    sd % delt   = delt
    sd % ps(0)  = pint(sd%kmax+1);
    sd % zs(0)  = zint(sd%kmax+1);
    sd % ktopconv = 1

    do k=1, sd%kmax
       nk=sd%kmax-k+1
       sd % p     (k) = pmid(nk)
       sd % z     (k) = zmid(nk);
       sd % ps    (k) = pint(nk); 
       sd % zs    (k) = zint(nk); 
       sd % t     (k) = t   (nk)
       !prevent negative values for qv,ql,qi,qa,qn
       sd % qv    (k) = max(qv(nk), 4.e-10) 
       sd % ql    (k) = max(ql(nk), 0.)
       sd % qi    (k) = max(qi(nk), 0.)
       sd % qa    (k) = max(qa(nk), 0.)
       sd % qn    (k) = max(qn(nk), 0.)
       sd % u     (k) = u(nk)
       sd % v     (k) = v(nk)
       !yim's aerosol
       sd % am1  (k) = am1(nk)
       sd % am2  (k) = am2(nk)
       sd % am3  (k) = am3(nk)
       sd % am4  (k) = am4(nk)
!++++yim
       do m=1, size(tracers,2)
          sd % tr (k,m) = tracers (nk,m)
       end do
       if (sd % p (k) > ptopconv) sd % ktopconv = k
    end do
  end subroutine pack_sd_k

!#####################################################################
!#####################################################################

  subroutine extend_sd_k(sd, pblht, doice, Uw_p)
    type(sounding), intent(inout) :: sd
    real, intent(in)              :: pblht
    logical, intent(in)           :: doice
    type(uw_params), intent(inout)    :: Uw_p

    integer :: k, kl, ktoppbl
    real    :: sshl0a, sshl0b, ssthc0a, ssthc0b, ssqct0a, ssqct0b
    real    :: hl0bot, thc0bot, qct0bot, hl0top, thc0top, qct0top
    real    :: thj, qvj, qlj, qij, qse, qs_sum, qt_sum, dpsum
    real, dimension(size(sd%tr,2)) :: sstr0a, sstr0b

    sd % exners(0) = exn_k(sd%ps(0),Uw_p);
    if (doice) then
       sd%nu(:)= max(min((268. - sd % t(:))/20.,1.0),0.0);
    else
       sd%nu(:)=0.
    end if
    sd % leff(:) = (1-sd%nu(:))*Uw_p%HLv + sd%nu(:)*Uw_p%HLs
    sd % qct (:) = sd%qv(:)+sd%ql(:)+sd%qi(:)
    sd % hl  (:) = Uw_p%cp_air*sd%t(:)+Uw_p%grav*sd%z(:)-  &
                   sd%leff(:)*(sd%ql(:)+sd%qi(:))
    sd % qint = 0.
    sd % crh = 0.; qs_sum=0.; qt_sum=0.; !dpsum=0.;
    do k=1, sd%ktopconv !sd%kmax
       sd % dp    (k) = sd%ps(k-1)-sd%ps(k)
       sd % dz    (k) = sd%zs(k)  -sd%zs(k-1)
       sd % exner (k) = exn_k(sd%p (k), Uw_p)
       sd % exners(k) = exn_k(sd%ps(k),Uw_p)
       sd % thc   (k) = sd%t(k) / sd%exner(k) 
       sd % qs    (k) = qsat_k(sd%t(k), sd%p(k),Uw_p, qv=sd%qv(k))
       sd % rh    (k) = min(sd%qv(k)/sd%qs(k),1.)
       sd % thv   (k) = sd%t(k)/sd%exner(k) *   &
                        (1.+Uw_p%zvir*sd%qv(k)-sd%ql(k)-sd%qi(k))
       sd % rho   (k) = sd % p(k)/     &
                        (Uw_p%rdgas * sd % thv(k) * sd % exner(k))
       sd % qint      =sd % qint + sd%qct(k)*sd%dp(k)
       qs_sum = qs_sum + sd % qs(k)  * sd%dp(k)
       qt_sum = qt_sum + sd % qct(k) * sd%dp(k)
       !dpsum = dpsum + sd%dp(k)
    end do
    sd % qint = sd % qint / Uw_p%grav
    sd % crh  = qt_sum / qs_sum
    sd % hm  (:) = Uw_p%cp_air*sd%t(:)+Uw_p%grav*sd%z(:)+sd%leff(:)*sd%qv(:)
    sd % hms (:) = Uw_p%cp_air*sd%t(:)+Uw_p%grav*sd%z(:)+sd%leff(:)*sd%qs(:)

   !Finite-Volume intepolation
    kl=sd%ktopconv !sd%kmax-1
    sshl0b  = (sd%hl (2)-sd%hl (1))/(sd%p(2)-sd%p(1))
    ssthc0b = (sd%thc(2)-sd%thc(1))/(sd%p(2)-sd%p(1))
    ssqct0b = (sd%qct(2)-sd%qct(1))/(sd%p(2)-sd%p(1))
    sstr0b(:) = (sd%tr(2,:)-sd%tr(1,:))/(sd%p(2)-sd%p(1))

    do k=2,kl
       sshl0a  = (sd%hl (k)-sd%hl (k-1))/(sd%p(k)-sd%p(k-1))
       if(sshl0a.gt.0)then
          sd%sshl (k-1) = max(0.,min(sshl0a,sshl0b))
       else
          sd%sshl (k-1) = min(0.,max(sshl0a,sshl0b))
       endif
       sshl0b = sshl0a
       ssthc0a = (sd%thc(k)-sd%thc(k-1))/(sd%p(k)-sd%p(k-1))
       if(ssthc0a.gt.0)then
          sd%ssthc(k-1) = max(0.,min(ssthc0a,ssthc0b))
       else
          sd%ssthc(k-1) = min(0.,max(ssthc0a,ssthc0b))
       endif
       ssthc0b = ssthc0a
       ssqct0a = (sd%qct(k)-sd%qct(k-1))/(sd%p(k)-sd%p(k-1))
       if(ssqct0a.gt.0)then
          sd%ssqct(k-1) = max(0.,min(ssqct0a,ssqct0b))
       else
          sd%ssqct(k-1) = min(0.,max(ssqct0a,ssqct0b))
       endif
       ssqct0b = ssqct0a
       sstr0a(:) = (sd%tr(k,:)-sd%tr(k-1,:))/(sd%p(k)-sd%p(k-1))
       where (sstr0a(:) > 0)
          sd%sstr(k-1,:) = max(0.,min(sstr0a(:),sstr0b(:)))
       elsewhere
          sd%sstr(k-1,:) = min(0.,max(sstr0a(:),sstr0b(:)))
       end where
       sstr0b(:) = sstr0a(:)
    enddo
    do k = 2,kl-1 !wind shear
       sd%dudp(k) = (sd%u(k+1)-sd%u(k-1))/(sd%p(k+1)-sd%p(k-1))
       sd%dvdp(k) = (sd%v(k+1)-sd%v(k-1))/(sd%p(k+1)-sd%p(k-1))
    end do
    sd%sshl (kl)=sd%sshl (kl-1)
    sd%ssthc(kl)=sd%ssthc(kl-1)
    sd%ssqct(kl)=sd%ssqct(kl-1)
    sd%sstr(kl,:)=sd%sstr(kl-1,:)

    do k = 1,sd%ktopconv !kl cannot be pver since ps0(pver)=0
       hl0bot  = sd%hl (k)+sd%sshl (k)*(sd%ps(k-1)-sd%p(k))
       thc0bot = sd%thc(k)+sd%ssthc(k)*(sd%ps(k-1)-sd%p(k))
       qct0bot = sd%qct(k)+sd%ssqct(k)*(sd%ps(k-1)-sd%p(k))
       call findt_k(sd%zs(k-1),sd%ps(k-1),hl0bot,qct0bot,thj,  &
                    qvj,qlj,qij,qse,sd%thvbot(k),doice, Uw_p)
       hl0top  = sd%hl (k)+sd%sshl (k)*(sd%ps(k)-sd%p(k))
       thc0top = sd%thc(k)+sd%ssthc(k)*(sd%ps(k)-sd%p(k))
       qct0top = sd%qct(k)+sd%ssqct(k)*(sd%ps(k)-sd%p(k))
       call findt_k(sd%zs(k),sd%ps(k),hl0top,qct0top,thj,  &
                    qvj,qlj,qij,qse,sd%thvtop(k),doice, Uw_p)
    enddo

    ktoppbl=1;
    do k = 2, kl-1
       if ((pblht+sd%zs(0)+1.-sd%zs(k))*(pblht+sd%zs(0)+1.-sd%zs(k+1)).lt.0.) then
          ktoppbl=k; exit;
       endif
    end do
    !given a layer index k (here k=kinv); !its bottom interface 
    !level pressure is ps0(k-1) [here ps0(kinv-1)] and its bottom
    !interface level virt. pot. temperature is thv(k) [thv0bot(kinv)]
    sd % ktoppbl = ktoppbl
    sd % kinv    = ktoppbl+1 
    sd % pinv    = sd % ps    (sd % kinv-1) 
    sd % zinv    = sd % zs    (sd % kinv-1) 
    sd % thvinv  = sd % thvbot(sd % kinv)
    sd % pblht   = pblht
  end subroutine extend_sd_k

!#####################################################################
!#####################################################################

  subroutine adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd,   &
                          Uw_p, dofast, doice, ac, rmuz)
  
    real, intent(inout) :: zsrc, psrc, hlsrc, thcsrc, qctsrc
    type(sounding), intent(in)    :: sd 
    type(uw_params), intent(inout)    :: Uw_p
    logical,        intent(in)    :: dofast, doice
    type(adicloud), intent(inout) :: ac
    real, intent(in), optional :: rmuz
    
    integer :: k, kl, klcl
    real    :: qs
    real    :: hl0lcl, thc0lcl, qct0lcl, thv0lcl, rho0lcl
    real    :: cin, cinlcl, plfc, thvubot, thvutop
    real    :: thj, qvj, qlj, qij, qse, thvj
    real    :: cape, plnb, chi, tmp, rhtmp
    real    :: alpha

    call ac_clear_k(ac);
    ac%klcl=0; ac%klfc=0; ac%klnb=0; 
    ac%plcl=0.; ac%zlcl=0.; ac%thvlcl=0.; ac%thv0lcl=0; ac%rho0lcl=0.;
    ac%plfc=0.; ac%plnb=0.; ac%cape=0.; ac%cin=0.;

!!$!below is pzhu's version now commented out
!!$    esrc=qctsrc*psrc/100./(qctsrc+epsilo)             ! water vapor pressure
!!$    tdsrc=tfreeze/(1-tfreeze*rvgas*log(esrc/6.11)/HLv)! dew-point of source air
!!$    temsrc=thcsrc*exn(psrc)                           ! temperature of source air
!!$    zlcl=123.5*(temsrc-tdsrc)+zsrc                    ! from sea-level
!!$    tlcl=temsrc-0.0098*(zlcl-zsrc)
!!$    ac % zlcl =zlcl
!!$    ac % plcl =psrc*(tlcl/temsrc)**(1./Kappa)


    !calculate lifted condensation level of air at parcel origin level
    !(within 0.2% of formula of Bolton, mon. wea. rev.,1980)
    call findt_k(zsrc, psrc, hlsrc, qctsrc, thj, qvj,    &
                 qlj, qij, qse, thvj, doice, Uw_p)
    tmp=thj*exn_k(psrc,Uw_p)
    rhtmp=min(qctsrc/qse,1.)
    chi=tmp/(1669.0-122.0*rhtmp-tmp)
    ac%plcl=psrc*(rhtmp**chi); !Emanuel's calculation, results nearly identical to RAS

    klcl=0;  !klcl is the layer containing the LCL, i.e., ps0(klcl)<=plcl(i,j)
    do k=1,sd % ktopconv-1
       if(sd%ps(k).le.ac%plcl) then
          klcl=k; 
          ac%zlcl=sd%zs(k)-(ac%plcl-sd%ps(k))/sd%dp(k)*sd%dz(k);
          exit
       else
          klcl   =sd % ktopconv
          ac%zlcl=sd % zs(klcl)
       end if
    end do
    if (sd%ps(1).le.ac%plcl) then
       klcl=2; ac%plcl=sd%ps(1); ac%zlcl=sd%zs(1); 
    end if
    ac % klcl=klcl; 

    if (dofast.and.(ac%klcl.eq.0 .or. ac%plcl.gt.sd%ps(1) .or. ac%plcl.lt.20000.)) return;

    call findt_k(ac%zlcl, ac%plcl, hlsrc, qctsrc, thj,   &
                 qvj, qlj, qij, qse, ac%thvlcl, doice, Uw_p)
    ac % hlsrc  = hlsrc
    ac % thcsrc = thcsrc
    ac % qctsrc = qctsrc

    hl0lcl  = sd%hl (klcl)+sd%sshl (klcl)*(ac%plcl-sd%p(klcl))
    thc0lcl = sd%thc(klcl)+sd%ssthc(klcl)*(ac%plcl-sd%p(klcl))
    qct0lcl = sd%qct(klcl)+sd%ssqct(klcl)*(ac%plcl-sd%p(klcl))
    call findt_k(ac%zlcl,ac%plcl,hl0lcl,qct0lcl,thj,qvj,  &
                 qlj,qij,qse,thv0lcl, doice, Uw_p)
    rho0lcl = ac%plcl/(Uw_p%rdgas*thv0lcl*exn_k(ac%plcl,Uw_p))
    ac % thv0lcl= thv0lcl
    ac % rho0lcl= rho0lcl


    kl=sd % ktopconv-1

    if (present(rmuz)) then
      if (rmuz /= 0.0) then
        do k=1,kl
         call findt_k(sd%zs(k),sd%ps(k), hlsrc, qctsrc, thj, ac%qv(k), &
                      ac%ql(k), ac%qi(k), qs, ac%thv(k), doice, Uw_p)
         ac%t(k) = thj*exn_k(sd%ps(k),Uw_p)
         alpha = MIN(1.0, rmuz*(sd%zs(k+1) - sd%zs(k)))
         hlsrc = (1.0-alpha)*hlsrc + alpha*sd%hl(k+1)
         qctsrc =(1.0-alpha)*qctsrc + alpha*sd%qct(k+1)
       end do
    else
      do k=1,kl
        call findt_k(sd%zs(k),sd%ps(k), hlsrc, qctsrc, thj, ac%qv(k), &
                     ac%ql(k), ac%qi(k), qs, ac%thv(k), doice, Uw_p)
        ac%t(k) = thj*exn_k(sd%ps(k),Uw_p)
      end do
    endif
  else
    do k=1,kl
       call findt_k(sd%zs(k),sd%ps(k), hlsrc, qctsrc, thj, ac%qv(k), &
                    ac%ql(k), ac%qi(k), qs, ac%thv(k), doice, Uw_p)
       ac%t(k) = thj*exn_k(sd%ps(k),Uw_p)
       ac%buo(k) = ac%thv(k) - sd%thvtop(k)
    end do
  endif


    !Determine the convective inhibition (CIN)
    CIN = 0.
    cinlcl = 0.
    plfc   = 0.

    !define CIN based on LFC  
    do k = sd % kinv, kl-1
       if(k.eq.klcl-1) then !klcl-1 < layer < klcl
          thvubot=ac % thv (k); thvutop=ac % thvlcl
          call getcin_k(sd%ps(k),sd%thvtop(k),ac%plcl,  &
                        thv0lcl,thvubot,thvutop,plfc,cin,Uw_p)
          cinlcl = cin
          thvubot=thvutop; thvutop=ac % thv (k+1)
          call getcin_k(ac%plcl,thv0lcl,sd%ps(k+1),sd%thvtop(k+1),  &
                        thvubot,thvutop,plfc,cin,Uw_p)
          if(plfc.gt.0. .and. plfc.lt.ac%plcl) exit
       else
          thvubot=ac % thv (k); thvutop=ac % thv (k+1)
          call getcin_k(sd%ps(k),sd%thvtop(k),sd%ps(k+1),  &
                        sd%thvtop(k+1),thvubot,thvutop,plfc, cin, Uw_p)
          if(plfc.gt.0. .and. plfc.lt.ac%plcl) exit
       endif
    enddo
 
    ac % cin =cin; !CIN has been estimated
    ac % plfc=plfc;

    if (dofast .and. (ac%plfc .lt. 500.) ) return; !miz

    !calculate cape=================
    if (ac%plfc.eq.0.0) then
       ac%cape=0.0;
       ac%plnb=0.0;
    else
       ac % klfc=0; !klfc is the layer containing the plfc, i.e., ps0(klfc)<=plfc(i,j)
       do k=1,kl 
          if(sd%ps(k).le.ac%plfc) then
             ac % klfc=max(k,2); exit
          end if
       end do
       plnb = 0.; cape=0.0;
       do k = ac % klfc-1, kl !for m45l48: sd%kmax
          thvubot=ac % thv (k); thvutop=ac % thv (k+1)
          call getcape_k(sd%ps(k),sd%thvtop(k),sd%ps(k+1),  &
                         sd%thvtop(k+1),thvubot,thvutop,plnb,cape,Uw_p)
          if(plnb.gt.0.) exit
       enddo
       ac%cape=cape
       ac%plnb=plnb
    end if

  end subroutine adi_cloud_k

!#####################################################################
!#####################################################################

 function qsat_k(temp, p,Uw_p, qv)
   real, intent(in)    :: temp, p
   type(uw_params), intent(inout) :: Uw_p
   real, intent(in), optional :: qv
   real :: qsat_k, t
   integer :: ier
 
   t = min(max(temp,Uw_p%tkmin),Uw_p%tkmax)
   if (present(qv)) then
     call compute_qs_k (t, p, Uw_p%epsilo, Uw_p%zvir, qsat_k, ier, &
                                                               q = qv)
   else
     call compute_qs_k (t, p, Uw_p%epsilo, Uw_p%zvir, qsat_k, ier)
   endif
 
   return
   end function qsat_k


!#####################################################################
!#####################################################################

  subroutine qses_k(temp, p, qs, es,Uw_p                     )
   real, intent(in)    :: temp, p
   type(uw_params), intent(inout) :: Uw_p
   real, intent(inout) :: qs, es
   real :: t
   integer :: ier

   t = min(max(temp,Uw_p%tkmin),Uw_p%tkmax)
   call compute_qs_k (t, p, Uw_p%epsilo, Uw_p%zvir, qs, ier)
  
   return
   end subroutine qses_k


!#####################################################################
!#####################################################################

! Subroutine to initialize the lookup table used to speed up
! the computation of the exner function (cjg)
 
  subroutine exn_init_k(Uw_p)       
 
    type(uw_params), intent(in) :: Uw_p
    

    integer k
    real p

!   Initialize 1d lookup table for exner function

  if ( .not. ex_lookup_allocated) then
!   if ( allocated(ex_lookup) ) deallocate(ex_lookup)
    allocate( ex_lookup(npex) )

    dpex = (pexmax-pexmin)/(npex-1)
    rdpex = 1.0 / dpex
    do k=1,npex
      p = (k-1)*dpex + pexmin
      ex_lookup(k) = (p/Uw_p%p00) ** Uw_p%kappa
    end do
    ex_lookup_allocated = .true.
 endif

  end subroutine exn_init_k

!#####################################################################
!#####################################################################

  subroutine exn_end_k()
 
    if ( allocated(ex_lookup) ) deallocate(ex_lookup)
    ex_lookup_allocated = .false.
 
  end subroutine exn_end_k
 
!#####################################################################
!#####################################################################

! Subroutine to compute the exner function using a lookup
! table for better performance (cjg)

  function exn_k(p,Uw_p)      
 
    real :: exn_k
    real, intent(in)  :: p
    type(uw_params), intent(inout) :: Uw_p

    integer k, kp1
    real w

    k = 0
    if ( p-pexmin.gt.0.0 ) k = int( (p-pexmin)*rdpex ) + 1
 
    if ( k.ge.1 .and.  k.lt.npex ) then
      kp1 = k+1
      w = ( p  - (k-1)*dpex  - pexmin  )*rdpex
      exn_k = (1-w)*ex_lookup(k) + w*ex_lookup(kp1)
    else
      exn_k = (p/Uw_p%p00) ** Uw_p%kappa
    end if
 
  end function exn_k

! Old subroutine that doesn't use a lookup table

! function exn_k(p,Uw_p)
!   real :: exn_k
!   real, intent(in)  :: p
!  type(uw_params), intent(inout) :: Uw_p
!   exn_k = (p/Uw_p%p00) ** Uw_p%Kappa
! end function exn_k

!#####################################################################
!#####################################################################

  subroutine conden_k(p,thc,qt,th,qv,ql,qi,qs,thv, Uw_p)
    real,     intent(in)  :: p, thc, qt
   type(uw_params), intent(inout) :: Uw_p
    real,     intent(out) :: th, qv, ql, qi, qs, thv
    real      tc, exn, leff, nu, qc, temps, tc1
    integer   iteration, id_check
    integer :: niteration = 5

    exn = (p/Uw_p%p00)**Uw_p%Kappa
    tc = thc * exn
    nu = max(min((268.-tc)/20.,1.0),0.0);
    leff = (1-nu)*Uw_p%HLv + nu*Uw_p%HLs
  
    temps = tc
    qs=qsat_k(temps,p,Uw_p)
    
    if(qs.gt.qt) then
       id_check=0
    else
       do iteration = 1,niteration
          temps = temps + ((tc-temps)*Uw_p%cp_air/leff + (qt -qs))/        &
               (Uw_p%cp_air/leff + Uw_p%epsilo*leff*qs/Uw_p%rdgas/temps/temps)
          qs = qsat_k(temps,p,Uw_p)
       enddo
       tc1=temps-leff/Uw_p%cp_air*(qt-qs)
       if(abs(tc1-tc).lt.1.0) then
          id_check=0
       else
          id_check=1; print*,'ID_CHECK=11111111111111111111111111111111'
       endif
    endif
    qc = max(qt-qs, 0.)
    qv = qt - qc
    ql = (1-nu)*qc
    qi = nu*qc !temps=tc+leff/Uw_p%cp_air*qc
    th = temps/exn
    thv=th*(1.+Uw_p%zvir*qv-ql-qi)
  end subroutine conden_k

!#####################################################################
!#####################################################################
! Subroutine to initialize lookup tables used to speed up findt (cjg)

  subroutine findt_init_k (Uw_p)                               
 
    implicit none

!   real, intent(in) :: epsilo, hlv, hls, cp_air, tkmin, tkmax
    type(uw_params), intent(inout) :: Uw_p

    integer i, k
    real ta, p, t

!   Allocate memory to hold lookup table

  if (.not. ta_lookup_allocated) then
!   if ( allocated(ta_lookup) ) deallocate(ta_lookup)
    allocate( ta_lookup(nta,np1+np2,0:1) )

!   Initialize 2d look up tables for temperature conversion

    dta = (tamax-tamin)/(nta-1)
    rdta = 1.0 / dta
    dp1 = (p1max-p1min)/(np1-1)
    rdp1 = 1.0 / dp1
    dp2 = (p2max-p2min)/(np2-1)
    rdp2 = 1.0 / dp2

    do i=1,nta
      ta = (i-1)*dta + tamin
      do k=1,np1
        p = (k-1)*dp1 + p1min
         call solve_ta_k(ta,p,t,0,Uw_p)
         ta_lookup(i,k,0) = t
         call solve_ta_k(ta,p,t,1,Uw_p)
         ta_lookup(i,k,1) = t
       end do
       do k=np1+1,np1+np2
         p = (k-np1-1)*dp2 + p2min
         call solve_ta_k(ta,p,t,0,Uw_p)
         ta_lookup(i,k,0) = t
         call solve_ta_k(ta,p,t,1,Uw_p)
         ta_lookup(i,k,1) = t
       end do
     end do
     ta_lookup_allocated = .true.
  endif

  end subroutine findt_init_k
 
!#####################################################################
!#####################################################################

  subroutine findt_end_k()
 
    if ( allocated(ta_lookup) ) deallocate(ta_lookup)
    ta_lookup_allocated = .false.

  end subroutine findt_end_k

!#####################################################################
!#####################################################################

  subroutine getcin_k(pbot,thv0bot,ptop,thv0top,thvubot,  &
                      thvutop,plfc,cin,Uw_p)
    
    real,    intent(in)    :: pbot,thv0bot,ptop,thv0top,thvubot,thvutop
    real,    intent(inout) :: plfc,cin
    real                   :: frc, rhom, delp
   type(uw_params), intent(inout) :: Uw_p

    delp=(pbot-ptop)
    rhom = pbot/(Uw_p%rdgas*thv0bot*exn_k(pbot,Uw_p))+ptop/  &
           (Uw_p%rdgas*thv0top*exn_k(ptop,Uw_p))

    if(thvubot.gt.thv0bot.and.thvutop.gt.thv0top)then
       !Both top and bottom positively buoyant
       plfc = pbot
    elseif(thvubot.le.thv0bot.and.thvutop.le.thv0top)then 
       !Both top and bottom negatively buoyant
       cin  = cin  - ((thvubot/thv0bot-1.)+(thvutop/thv0top-1.)) * delp / rhom
    elseif(thvutop.le.thv0top.and.thvubot.gt.thv0bot)then
       !Top negatively buoyant; Bottom positively buoyant
       frc  = (thvutop/thv0top-1.)/((thvutop/thv0top-1.)-(thvubot/thv0bot-1.))
       delp = (ptop+frc*delp) - ptop
       cin  = cin - (thvutop/thv0top-1.) * delp / rhom
    else                                                  
       !Top positively buoyant; Bottom negatively buoyant
       frc = (thvubot/thv0bot-1.)/((thvubot/thv0bot-1.)-(thvutop/thv0top-1.))
       plfc = pbot - frc * (pbot-ptop)
       delp = pbot - plfc
       cin  = cin - (thvubot/thv0bot-1.) * delp / rhom
    endif
  end subroutine getcin_k

!#####################################################################
!#####################################################################

  subroutine getcape_k (pbot,thv0bot,ptop,thv0top,thvubot,  &
                        thvutop,plnb,cape,Uw_p)
   
    real,    intent(in)    :: pbot,thv0bot,ptop,thv0top,thvubot,thvutop
    real,    intent(inout) :: plnb,cape
   type(uw_params), intent(inout) :: Uw_p
    real                   :: frc, rhom, delp

    delp=(pbot-ptop)
    rhom = pbot/(Uw_p%rdgas*thv0bot*exn_k(pbot,Uw_p))+ptop/  &
           (Uw_p%rdgas*thv0top*exn_k(ptop,Uw_p))

    if(thvubot.gt.thv0bot.and.thvutop.gt.thv0top)then
       !Both top and bottom positively buoyant
       cape = cape + ((thvubot/thv0bot - 1.) + (thvutop/thv0top - 1.))*&
              delp/rhom
    elseif(thvubot.le.thv0bot.and.thvutop.le.thv0top)then 
       !Both top and bottom negatively buoyant
       plnb = pbot
    elseif(thvutop.le.thv0top.and.thvubot.gt.thv0bot)then
       !Top negatively buoyant; Bottom positively buoyant
       frc  = (thvubot/thv0bot-1.)/((thvubot/thv0bot-1.)-  &
              (thvutop/thv0top-1.))
       plnb = pbot - frc * (pbot-ptop)
       delp = pbot - plnb
       cape = cape + (thvubot/thv0bot-1.) * delp / rhom
    else                                                  
       !Top positively buoyant; Bottom negatively buoyant
       frc  = (thvutop/thv0top-1.)/((thvutop/thv0top-1.)- &
               (thvubot/thv0bot-1.))
       delp = (ptop+frc*delp) - ptop
       cape = cape + (thvutop/thv0top-1.) * delp / rhom
    endif

  end subroutine getcape_k


!###################################################################
!###################################################################

subroutine pack_sd_lsm_k (do_lands, land, coldT, dt, pf, ph, zf, zh, &
                          t, qv, tracers, sd)

  logical,            intent(in)    :: do_lands
  real,               intent(in)    :: land
  logical,            intent(in)    :: coldT
  real,               intent(in)    :: dt
  real, dimension(:), intent(in)    :: pf, ph, zf, zh, t, qv
!++++yim
  real, dimension(:,:), intent(in)    :: tracers
  type(sounding),     intent(inout) :: sd

!++++yim
  real, parameter :: ptopconv = 3000.
  integer :: k, nk, kmax, m

  kmax=size(t)
  sd % kmax   = kmax
  if (do_lands) then
    sd % land   = land  
    sd % coldT  = coldT    
  else
    sd % land   = 0.
    sd % coldT  = .false.
  endif 
  sd % delt   = dt   
  sd % ps(0)  = ph(kmax+1)
  sd % zs(0)  = zh(kmax+1)
  sd % ktopconv = 1

  do k=1, kmax
     nk=kmax-k+1
     sd % p (k) = pf(nk)
     sd % z (k) = zf(nk)
     sd % ps(k) = ph(nk)
     sd % zs(k) = zh(nk)
     sd % t (k) = t (nk)
     sd % qv(k) = max(qv(nk)/(1.+qv(nk)), 4.e-10) !for donner_deep where mixing-ratio passed in
     sd % ql(k) = 0.
     sd % qi(k) = 0.
     sd % qa(k) = 0.
     sd % qn(k) = 0.
     sd % u (k) = 0.
     sd % v (k) = 0.
     if (sd % p (k) > ptopconv) sd % ktopconv = k
!++++yim
       do m=1, size(tracers,2)
          sd % tr (k,m) = tracers (nk,m)
       end do
  end do

end subroutine pack_sd_lsm_k

!###################################################################
!###################################################################

! subroutine findt_k(z,p,hl,qt,th,qv,ql,qi,qs,thv,doice, Uw_p)
  subroutine findt_new_k(z,p,hl,qt,th,qv,ql,qi,qs,thv,doice, Uw_p)

      implicit none

      real, intent(in)  :: z, p, hl, qt
      real, intent(out) :: th, qv, ql, qi, qs, thv
      logical, intent(in) :: doice
      type(uw_params), intent(inout) :: Uw_p

      real hh, temp, temp_unsat, temp_sat, temp1, temp2, dtemp
      real tempmin, tempmax
      real f, f1, fmin, fmax
      real es, nu, leff, qc

      integer n, nmax
      integer hflag
      logical lbracket

      hflag = 1
      if (doice) hflag = 2

!     Definitely unsaturated case
!     The unsaturated temperature (temp_unsat) is always lower or equal 
!     to the actual temperature (temp). 
!     Therefore qs(temp_unsat) <= qs(temp), since qs(T) is monotically
!     increasing with T.

      temp = (hl-Uw_p%grav*z)/Uw_p%cp_air
      call qses_k(temp,p,qs,es,Uw_p)
      if ( qs.gt.qt ) then
        ql  = 0.
        qi  = 0.
        qv  = qt
        th  = temp/exn_k(p,Uw_p)
        thv = th*(1.+Uw_p%zvir*qv-ql-qi)
        return
      end if
      temp_unsat = temp

!     Possibly saturated case
!     Absolute bounds on the temperature

      hh = hl - Uw_p%grav*z
      tempmin = max( temp_unsat - 1.0, Uw_p%tkmin )
      tempmax = min( temp_unsat + Uw_p%hls*(qt-qs)/Uw_p%cp_air + 1.0,  &
                     Uw_p%tkmax )
      fmin = saturated_k(tempmin,hh,qt,p,hflag,Uw_p)
      fmax = saturated_k(tempmax,hh,qt,p,hflag,Uw_p)

!     The bounds on temperature are likely to be too large,
!     so we need to bracket the solution first. 
!     We search for the root closest to tempmin.

      lbracket = .false.
      dtemp = 10.0

!     Is the initial bracket good enough ?

      if ( tempmax-tempmin.le.dtemp ) then
        if ( fmin*fmax.le.0.0 ) then
          lbracket = .true.
          temp1 = tempmin
          temp2 = tempmax
        else
          dtemp = tempmax - tempmin
        endif
      endif

!     If not refine it

      do while (.not.lbracket .and. dtemp.ge.0.2)
      temp1 = tempmin
      temp2 = tempmax
      f1 = saturated_k(temp1,hh,qt,p,hflag,Uw_p)
      nmax = int( (temp2-temp1)/dtemp ) + 1
      temp = temp1
      do n=1,nmax
        temp = MIN(temp + dtemp, tempmax)
        f = saturated_k(temp,hh,qt,p,hflag,Uw_p)
        if (f1*f.le.0) then
          temp2 = temp
          lbracket = .true.
          exit
        endif
        temp1 = temp
        f1 = f
      enddo
      dtemp = 0.5 * dtemp
      enddo

!     Did we make it ?

      if (lbracket) then

!       Now find the root within the bracket

        temp_sat = zriddr_k(saturated_k,temp1,temp2,hh,qt,p,hflag,1.e-3,Uw_p)

!       Choose between one of the two choices, the highest value is temp

        if (temp_sat .gt. temp_unsat ) then
          temp = temp_sat
        else
          temp = temp_unsat
        endif

        call qses_k(temp,p,qs,es,Uw_p)
        if (doice) then
          nu = max(min((268.-temp)/20.,1.0),0.0)
        else
          nu = 0.0
        endif
        leff = (1.-nu)*Uw_p%hlv + nu*Uw_p%hls

        qc = max(qt-qs, 0.)
        qv = qt - qc
        ql = (1.-nu)*qc
        qi = nu*qc
        th = temp/exn_k(p,Uw_p)
        thv=th*(1.+Uw_p%zvir*qv-ql-qi)
        return

      else

!RSH2    write(*,*) 'WARNING findt_new_k: not bracketed'
!        write(*,*) 'Not bracketed i = ',i
!RSH :: need to properly process this error condition (if in fact it is
!       possible to occur)
!       write(*,*) 'Not bracketed temp1, temp2, dtemp, f1, f, tempmin&
!                    & tempmax, fmin, fmax = ', &
!               temp1, temp2, dtemp, f1, f, tempmin, tempmax, fmin, fmax
        temp = temp_unsat
        qv = qt
        ql = 0.0
        qi = 0.0
        th = temp/exn_k(p,Uw_p)
        thv=th*(1.+Uw_p%zvir*qv)
        return
        
      endif

      end subroutine findt_new_k
!     end subroutine findt_k

!     -----------------------------------------------------------------
!     To diagnose temp from hl and qt, we need to find the zero
!     of this function

      real function saturated_k(temp,hh,qt,p,hflag,Uw_p)
      implicit none

      real, intent(in) :: temp, hh, qt, p
      type(uw_params), intent(inout) :: Uw_p
  
      integer, intent(in) :: hflag

      real es, qs, leff, nu

      call qses_k(temp,p,qs,es,Uw_p)

      select case (hflag)
        case (0)
          leff = Uw_p%hls
        case (1)
          leff = Uw_p%hlv
        case (2)
          nu = max(min((268.0-temp)/20.,1.0),0.0)
          leff = (1.0-nu)*Uw_p%hlv + nu*Uw_p%hls
      end select
      saturated_k = hh - Uw_p%cp_air*temp + Leff*(qt-qs)

      return
      end function saturated_k

!###################################################################
!###################################################################
! Newest findt subroutine (cjg). Comparable in accuracy with
! the revised version (findt_new_k), but substantially
! faster thanks to the use of lookup tables.

! subroutine findt_fast_k(z,p,hl,qt,th,qv,ql,qi,qs,thv,doice, Uw_p)
  subroutine findt_k(z,p,hl,qt,th,qv,ql,qi,qs,thv,doice, Uw_p)

    implicit none

    real, intent(in)  :: z, p, hl, qt
    real, intent(out) :: th, qv, ql, qi, qs, thv
    logical, intent(in) :: doice
    type(uw_params), intent(inout) :: Uw_p

    integer i, il, ii, k
    integer ip1, kp1
    integer hflag

    real hh, tl, tal, tai, t0, t1, temp, temp_unsat, temp_sat
    real es
    real qc, nu, leff

    real u, w

    hflag = 1
    if (doice) hflag = 2

!   Definitely unsaturated case
!   The unsaturated temperature (temp_unsat) is always lower or equal 
!   to the actual temperature (temp). 
!   Therefore qs(temp_unsat) <= qs(temp), since qs(T) is monotically
!   increasing with T.

    tl = (hl-Uw_p%grav*z)/Uw_p%cp_air
    call qses_k(tl,p,qs,es,Uw_p)
    if ( qs.gt.qt .or. tl >= 372.9 ) then
      ql  = 0.
      qi  = 0.
      qv  = qt
      th  = tl/exn_k(p,Uw_p)          
      thv = th*(1.+Uw_p%zvir*qv-ql-qi)
      return
    end if
    temp_unsat = tl

!   Compute temperature assuming saturated air

    tai = ( hl - Uw_p%grav*z + Uw_p%hls*qt ) / Uw_p%cp_air
    tal = ( hl - Uw_p%grav*z + Uw_p%hlv*qt ) / Uw_p%cp_air

!   Are we within the lookup table range?

    il = 0
    if ( tal.gt.tamin ) il = int( (tal-tamin)*rdta ) + 1
    ii = 0
    if ( tai.gt.tamin ) ii = int( (tai-tamin)*rdta ) + 1
    k = 0
    if ( p.gt.p1min .and. p.le.p1max ) then
      k = int( (p-p1min)*rdp1 ) + 1
    else if ( p.gt.p2min ) then
      k = int( (p-p2min)*rdp2 ) + np1 + 1
    end if

    if (       il.ge.1 .and. il.lt.nta  &
         .and. ii.ge.1 .and. ii.lt.nta  &
         .and.  k.ge.1 .and.  k.lt.np1+np2 ) then

!     Inside lookup table range

!     Use bi-linear interpolation from lookup table values

      if (hflag.eq.2) then

!       Mixed phase case

!       Use bi-linear interpolation to find temperature values
!       from look up tabled

        kp1 = k+1
        if ( k.le.np1 ) then
          w = ( p  - (k-1)*dp1  - p1min  )*rdp1
        else
          w = ( p  - (k-np1-1)*dp2  - p2min  )*rdp2
        endif

!       t0 is the temperature assuming pure ice condensate

        i = ii
        ip1 = i+1
        u = ( tai - (i-1)*dta - tamin )*rdta
        t0 =   (1-u)*(1-w) * ta_lookup(i  ,k  ,0) &
             + (1-u)*w     * ta_lookup(i  ,kp1,0) &
             + u    *(1-w) * ta_lookup(ip1,k  ,0) &
             + u    *w     * ta_lookup(ip1,kp1,0)

!       t1 is the temperature assuming pure liquid condensate

        i = il
        ip1 = i+1
        u = ( tal - (i-1)*dta - tamin )*rdta
        t1 =   (1-u)*(1-w) * ta_lookup(i  ,k  ,1) &
             + (1-u)*w     * ta_lookup(i  ,kp1,1) &
             + u    *(1-w) * ta_lookup(ip1,k  ,1) &
             + u    *w     * ta_lookup(ip1,kp1,1)

!       Do either t0 or t1 fall inside the mixed phase temperature range?

        if (      (t0.ge.248.0 .and. t0.le.268.0)  &
             .or. (t1.ge.248.0 .and. t1.le.268.0) ) then

!         Yes, use t0 and t1 as initial brackets for the root. Because
!         t0 and t1 are derived from lookup tables, the actual bracket
!         must be increased slightly

          hh = hl - Uw_p%grav*z
          temp_sat = zriddr_k(sat1_k,t1-0.01,t0+0.01,hh,qt,p,hflag,1.e-3, &
                              Uw_p)

        elseif ( t0.lt.248.0 .and. t1.lt.248.0 ) then

!         No, the condensate is definitely pure ice

          temp_sat = t0

        elseif ( t0.gt.268.0 .and. t1.gt.268.0 ) then

!         No, the condensate is definitely pure liquid

          temp_sat = t1

        else

         write(*,*) 'WARNING findt_fast_k: never get there'

        endif

      else  ! hflag.ne.2

!     Not mixed phase case: hflag = 0 : ice only
!                           hflag = 1 : liquid only (doice = .false.)

      i = il
      ip1 = i+1
      kp1 = k+1
      u = ( tal - (i-1)*dta - tamin )*rdta
      if ( k.le.np1 ) then
        w = ( p  - (k-1)*dp1  - p1min  )*rdp1
      else
        w = ( p  - (k-np1-1)*dp2  - p2min  )*rdp2
      endif

      temp_sat =   (1-u)*(1-w) * ta_lookup(i  ,k  ,hflag) &
                 + (1-u)*w     * ta_lookup(i  ,kp1,hflag) &
                 + u    *(1-w) * ta_lookup(ip1,k  ,hflag) &
                 + u    *w     * ta_lookup(ip1,kp1,hflag)

      endif

    else

!     Outside lookup table range: use root finding algorithm

      call solve_hl_k(z,p,hl,qt,temp_sat,hflag,Uw_p)

    endif

!   Choose between one of the two choices (temp_unsat, temp_sat):
!   the highest value is temp

    temp = max( temp_unsat, temp_sat )

!   Compute output variables

    call qses_k(temp,p,qs,es,Uw_p)
    if (doice) then
      nu = max(min((268.-temp)/20.,1.0),0.0)
    else
      nu = 0.0
    endif
    leff = (1.-nu)*Uw_p%hlv + nu*Uw_p%hls

    qc = max(qt-qs, 0.)
    qv = qt - qc
    ql = (1.-nu)*qc
    qi = nu*qc
    th = temp/exn_k(p,Uw_p)
    thv=th*(1.+Uw_p%zvir*qv-ql-qi)

    return

!  end subroutine findt_fast_k
  end subroutine findt_k

!     -----------------------------------------------------------------
!     Subroutine to solve for temp as a function of z, p, hl, qt

      subroutine solve_hl_k(z,p,hl,qt,temp,hflag,Uw_p)
      implicit none

      real, intent(in)  :: z, p, hl, qt
      real, intent(out) :: temp
      integer, intent(in) :: hflag
      type(uw_params), intent(inout) :: Uw_p

      real hh, temp_unsat, temp_sat, temp1, temp2, dtemp
      real tempmin, tempmax
      real f, f1, fmin, fmax
      real qs, es

      integer n, nmax
      logical lbracket

!     Definitely unsaturated case
!     The unsaturated temperature (temp_unsat) is always lower or equal 
!     to the actual temperature (temp). 
!     Therefore qs(temp_unsat) <= qs(temp), since qs(T) is monotically
!     increasing with T.

      temp = (hl-Uw_p%grav*z)/Uw_p%cp_air
      call qses_k(temp,p,qs,es,Uw_p)        
      if ( qs.gt.qt ) return

!     Absolute bounds on the temperature. Return immediately
!     if the bounds are outside the range of the saturation vapor 
!     pressure lookup tables.

      tempmin = temp - 1.0
      tempmax = temp + Uw_p%hls*(qt-qs)/Uw_p%cp_air + 1.0

      if ( tempmin.lt.Uw_p%tkmin .and. tempmax.lt.Uw_p%tkmin ) then
        temp = (hl-Uw_p%grav*z+Uw_p%hls*max(qt-qs,0.))/Uw_p%cp_air
        return
      endif

      if ( tempmin.gt.Uw_p%tkmax .and. tempmax.gt.Uw_p%tkmax ) then
        temp = (hl-Uw_p%grav*z+Uw_p%hlv*max(qt-qs,0.))/Uw_p%cp_air
        return
      endif

      temp_unsat = temp
      hh = hl - Uw_p%grav*z
      fmin = sat1_k(tempmin,hh,qt,p,hflag,Uw_p)
      fmax = sat1_k(tempmax,hh,qt,p,hflag,Uw_p)  

!     The bounds on temperature are likely to be too large,
!     so we need to bracket the solution first. 
!     We search for the root closest to tempmin.

      lbracket = .false.
      dtemp = 10.0

!     Is the initial bracket good enough ?

      if ( tempmax-tempmin.le.dtemp ) then
        if ( fmin*fmax.le.0.0 ) then
          lbracket = .true.
          temp1 = tempmin
          temp2 = tempmax
        else
          dtemp=tempmax-tempmin
        endif
      endif

!     If not refine it

      do while (.not.lbracket .and. dtemp.ge.0.2)
      temp1 = tempmin
      temp2 = tempmax
      f1 = sat1_k(temp1,hh,qt,p,hflag,Uw_p)  
      nmax = int( (temp2-temp1)/dtemp ) + 1
      temp = temp1
      do n=1,nmax
        temp = min(temp + dtemp, tempmax)
        f = sat1_k(temp,hh,qt,p,hflag,Uw_p)
        if (f1*f.le.0) then
          temp2 = temp
          lbracket = .true.
          exit
        endif
        temp1 = temp
        f1 = f
      enddo
      dtemp = 0.5 * dtemp
      enddo

!     Did we manage to bracket the root ?

      if (lbracket) then

!       Yes, now find the root within the bracket

        temp_sat = zriddr_k(sat1_k,temp1,temp2,hh,qt,p,hflag,1.e-3, &
                            Uw_p)   

!       Choose between one of the two choices (temp_unsat, temp_sat):
!       the highest value is temp

        temp = max( temp_unsat, temp_sat )

      else

        if (Uw_p%me == 0) then
        write(*,'(a,4e20.12)') 'WARNING solve_hl_k: not bracketed',z,p,hl,qt
        endif
        temp = temp_unsat
        
      endif

      return
      end subroutine solve_hl_k

!     -----------------------------------------------------------------
!     Subroutine to solve for temp as a function of ta, p

      subroutine solve_ta_k(ta,p,t,hflag,Uw_p)
      implicit none

      real, intent(in)  :: ta, p
      real, intent(out) :: t
      integer, intent(in) :: hflag
      type(uw_params), intent(inout) :: Uw_p

      real temp, temp1, temp2, dtemp
      real tempmin, tempmax
      real f, f1, fmin, fmax

      integer n, nmax
      logical lbracket

!     Absolute bounds on the temperature

      tempmin = Uw_p%tkmin
      tempmax = ta + 1.0
      fmin = sat2_k(tempmin,ta,p,0.0,hflag,Uw_p)  
      fmax = sat2_k(tempmax,ta,p,0.0,hflag,Uw_p)

!     The bounds on temperature are likely to be too large,
!     so we need to bracket the solution first. 
!     We search for the root closest to tempmin.

      lbracket = .false.
      dtemp = 50.0

!     Is the initial bracket good enough ?

      if ( tempmax-tempmin.le.dtemp ) then
        if ( fmin*fmax.le.0.0 ) then
          lbracket = .true.
          temp1 = tempmin
          temp2 = tempmax
        else
          dtemp=tempmax-tempmin
        endif
      endif

!     If not refine it

      do while (.not.lbracket .and. dtemp.ge.0.1)
      temp1 = tempmin
      temp2 = tempmax
      f1 = sat2_k(temp1,ta,p,0.0,hflag,Uw_p)
      nmax = int( (temp2-temp1)/dtemp ) + 1
      temp = temp1
      do n=1,nmax
        temp = min(temp + dtemp, tempmax)
        f = sat2_k(temp,ta,p,0.0,hflag,Uw_p)
        if (f1*f.le.0) then
          temp2 = temp
          lbracket = .true.
          exit
        endif
        temp1 = temp
        f1 = f
      enddo
      dtemp = 0.5 * dtemp
      enddo

!     Did we manage to bracket the root ?

      if (lbracket) then

!       Yes, now find the root within the bracket
        t = zriddr_k(sat2_k,temp1,temp2,ta,p,0.0,hflag,1.e-10,Uw_p)
!                    epsilo,hlv,hls,cp_air,tkmin,tkmax)
        return

      else

!       No, bracketing failed
        write(*,*) 'WARNING solve_ta_k: not bracketed'
        t = -999.0
        return
        
      endif

      end subroutine solve_ta_k

!     -----------------------------------------------------------------
!     To diagnose temp from hl and qt, we need to find the zero
!     of this function

      real function sat1_k(temp,hh,qt,p,hflag,Uw_p)  
      implicit none

      real, intent(in) :: temp, hh, qt, p
      integer, intent(in) :: hflag
      type(uw_params), intent(inout) :: Uw_p

      integer ier
      real qs, leff, nu
      real t

!     In-line qses computation here for better performance
      t = min(max(temp,Uw_p%tkmin),Uw_p%tkmax)
!   INLINING ????
      call compute_qs_k (t, p, Uw_p%epsilo, Uw_p%zvir, qs, ier)

      select case (hflag)
        case (0)
          leff = Uw_p%hls
        case (1)
          leff = Uw_p%hlv
        case (2)
          nu = max(min((268.0-temp)/20.,1.0),0.0)
          leff = (1.0-nu)*Uw_p%hlv + nu*Uw_p%hls
      end select
      sat1_k = hh - Uw_p%cp_air*temp + leff*(qt-qs)

      return
      end function sat1_k

!     -----------------------------------------------------------------
!     To diagnose temp from ta and p, we need to find the zero
!     of this function

      real function sat2_k(temp,ta,p,tmp,hflag, Uw_p)
      implicit none

      real, intent(in) :: temp, ta, p, tmp
      integer, intent(in) :: hflag
      type(uw_params), intent(inout) :: Uw_p

      integer ier
      real qs, leff, nu
      real t

!     In-line qses computation here for better performance
      t = min(max(temp,Uw_p%tkmin),Uw_p%tkmax)
!   INLINING ????
      call compute_qs_k (t, p, Uw_p%epsilo, Uw_p%zvir, qs, ier)

      select case (hflag)
        case (0)
          leff = Uw_p%hls
        case (1)
          leff = Uw_p%hlv
        case (2)
          nu = max(min((268.0-temp)/20.,1.0),0.0)
          leff = (1.0-nu)*Uw_p%hlv + nu*Uw_p%hls
      end select
      sat2_k = ta - temp - leff/Uw_p%cp_air*qs

      return
      end function sat2_k

!     -----------------------------------------------------------------
!     Function to find zero of function 'saturated' using Ridders'
!     method

      real function zriddr_k(func, x1,x2,ya,yb,yc,ld,xacc,Uw_p)
      implicit none

      integer maxit
      real unused
      type(uw_params) Uw_p
      parameter (maxit=60,unused=-1.11e30)

      real func
      external func

      real x1,x2,xacc
      real ya,yb,yc
      integer ld
      
      integer j
      real fh,fl,fm,fnew,s,xh,xl,xm,xnew

      fl=func       (x1,ya,yb,yc,ld,Uw_p)
      fh=func       (x2,ya,yb,yc,ld,Uw_p)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr_k=unused
        do j=1,maxit
          xm=0.5*(xl+xh)
          fm=func       (xm,ya,yb,yc,ld,Uw_p)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
          if (abs(xnew-zriddr_k).le.xacc) return
          zriddr_k=xnew
          fnew=func       (zriddr_k,ya,yb,yc,ld,Uw_p)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr_k
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr_k
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr_k
            fl=fnew
          else
!            write(*,*) 'WARNING: never get here in zriddr'
!            write(*,*) 'WARNING zriddr_k: never get there'
          endif
          if(abs(xh-xl).le.xacc) return
        enddo
!       write(*,*) 'WARNING zriddr_k: exceeded maximum iterations'
!        write(*,*) 'WARNING: zriddr exceed maximum iterations'
!  need to properly handle error condition
!        write(*,*) 'WARNING: zriddr exceed maximum iterations', &
!         x1, x2, xl, xh
      else if (fl.eq.0.) then
        zriddr_k=x1
      else if (fh.eq.0.) then
        zriddr_k=x2
      else
!        write(*,*) 'WARNING zriddr_k: root must be bracketed'
!        write(*,*) 'WARNING: root must be bracketed in zriddr'
      endif

      return
      end function zriddr_k

!######################################################################

!++lwh
subroutine check_tracer_realizability(kmax, ntr, dt, &
                                             tracers, trten, trwet)
!---------------------------------------------------------------------
!  Check for tracer realizability. If convective tendencies would
!  produce negative tracer mixing ratios, scale down tracer tendency
!  terms uniformly for this tracer throughout convective column. This is
!  equivalent to limiting the cell areas.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  Dummy arguments
!---------------------------------------------------------------------
integer,                 intent(in)     :: kmax, ntr
real,                    intent(in)     :: dt 
real, dimension(kmax,ntr), &
                         intent(in)     :: tracers        
real,dimension(kmax,ntr),intent(inout)  :: trten, trwet

!---------------------------------------------------------------------
!   intent(in) variables:
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     kmax           number of model layers in large-scale model
!     dt             physics time step [ sec ]
!
!   intent(inout) variables:
!     trten          tracer tendency
!     trwet          tracer wet deposition tendency
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  Local variables
!---------------------------------------------------------------------

   integer :: n,k
   real, dimension(kmax) :: tracer0, trtend, trtendw, tracer1, tracer1w
   real :: ratio, tracer_max, tracer_min

!---------------------------------------------------------------------
!   local variables:
!
!     tracers        tracer mixing ratios of tracers transported by the
!                    donner deep convection parameterization
!                    [ tracer units, e.g., kg(tracer) / kg (dry air) ]
!     tracer0        column tracer mixing ratios before convection
!     trtend         column tracer mixing ratio tendencies due to convection [ (tracer units) / s ]
!     tracer1        column tracer mixing ratios after convection
!     k, n     do-loop indices
!     ratio          ratio by which tracer convective tendencies need to 
!                    be reduced to permit realizability (i.e., to prevent
!                    negative tracer mixing ratios)
!
!---------------------------------------------------------------------

   do n = 1,ntr
      
      tracer0(:)  = tracers(:,n)
      trtend(:)   = trten(:,n)
      trtendw(:)  = trtend(:) + trwet(:,n)
      tracer1(:)  = tracer0 + dt * trtend(:)
      tracer1w(:) = tracer0 + dt * trtendw(:)

      tracer_min = 1.e20
      tracer_max = -1.e20

      do k = 1,kmax
         if (trtend(k) /= 0.) then
            tracer_max = max(tracer0(k),tracer_max)
            tracer_min = min(tracer0(k),tracer_min)
         end if
      end do
 
      ratio = 1.
      do k = 1,kmax
         if (tracer0(k)>0. .and. tracer1w(k)<0. ) then
            ratio = MIN( ratio,tracer0(k)/(-trtendw(k)*dt) )
         end if
         if (tracer1(k)<tracer_min .and. trtend(k) /= 0.0 ) then
            ratio = MIN( ratio,(tracer0(k)-tracer_min)/(-trtend(k)*dt) )
         end if
         if (tracer1(k)>tracer_max  .and. trtend(k) /= 0.0 ) then
            ratio = MIN( ratio,(tracer_max-tracer0(k))/(trtend(k)*dt) )
         end if
      end do
      ratio = MAX(0.,MIN(1.,ratio))
      if (ratio /= 1.) then
         trten(:,n) =  trten(:,n)*ratio
         trwet(:,n) =  trwet(:,n)*ratio
      end if
   end do


end subroutine check_tracer_realizability
!--lwh

!#####################################################################
subroutine qt_parcel_k (qs, qstar, pblht, tke, land, gama, pblht0, tke0, lofactor0, &
     lochoice, qt, lofactor)
    real,    intent(in)    :: qs, qstar, pblht, tke, land, gama, pblht0, tke0, lofactor0
    integer, intent(in)    :: lochoice
    real,    intent(inout) :: qt, lofactor

    real :: qttmp

    if (lochoice .eq. 0) then
       lofactor = 1. - land * (1. - lofactor0)
    elseif (lochoice .eq. 1) then
       lofactor = pblht0 / max(pblht,  pblht0)
    elseif (lochoice .eq. 2) then
       lofactor = tke0   / max(tke, tke0  )
    elseif (lochoice .eq. 3) then
       lofactor = tke0   / max(tke, tke0  )
       lofactor = sqrt(lofactor)
    else
       lofactor = 1.
    end if

    qttmp = qt*(1. + gama * land)
    qt    = max(qt, min(qttmp, qs))
 
  end subroutine qt_parcel_k


end MODULE CONV_UTILITIES_k_MOD
