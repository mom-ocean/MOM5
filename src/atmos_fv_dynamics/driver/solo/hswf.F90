!****6***0*********0*********0*********0*********0*********0**********72
      subroutine hswf(im, jm, km, jfirst, jlast,                    &
                      u,v,pt,pe, delp, peln, pkz,pdt,akap0,grav,rg,  &
                      strat, rayf, master, sinp,cosp,sine,cose,    &
                      ng_s, ng_d)
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:

#if defined( SPMD )
      use mod_comm, only :  gid, mp_send3d, mp_recv3d
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km
      integer jfirst, jlast
      integer ng_s, ng_d
      integer pdt
      real akap0, grav, rg
      logical strat
      logical rayf
      logical master
      real cosp(jm),sinp(jm),cose(jm),sine(jm)

! !INPUT/OUTPUT PARAMETERS:

      real :: u(im,jfirst-ng_d:jlast+ng_s,km)
      real :: v(im,jfirst-ng_d:jlast+ng_d,km)
      real :: pt(im,jfirst-ng_d:jlast+ng_d,km)
      real :: delp(im,jfirst:jlast,km)
      real ::   pe(im,km+1,jfirst:jlast)
      real :: peln(im,km+1,jfirst:jlast)
      real :: pkz(im,jfirst:jlast,km)


! !DESCRIPTION:
!    Author: Shian-Jiann Lin, JCET, UMBC and NASA/GSFC
!
! !REVISION HISTORY:
!   SJL 99.09.30:  Delivery
!   WS  99.10.28:  jfirst:jlast; documentation
!   WS  99.11.07:  pruned arrays
!   WS  00.07.10:  Incorporated simplfications to PILGRIM
!
!EOP
!-----------------------------------------------------------------------
!BOC
!

      real  p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real  tmp
      real  ap0k, algpk, rfc
      real  tey, tez, fac, pw, sigl
      real  pl(im,km),frac(im,jm)
      real  teq(im,km)
      real  h0, dz
      real  dt_tropic
      real  rmr, rms
      real  relx, tau
      real  t_st, t_ms
      real  f1
      real  ptop
      real  ak(km+1), bk(km+1)
      real  pc, c1
      real  t2(im, km)
      real  dp(im, km)

      real  ua(im,jfirst:jlast,km)
      real  va(im,jfirst:jlast,km)
      real  u2(im,km)
      real  v2(im,km)
      real  fu(im,km)
      real  fv(im,km)

#if defined( SPMD )
      real  , allocatable :: pesouth(:,:), uasouth(:,:)    !will be removed later
!      real  :: pesouth(im, km+1), uasouth(im, km)    !will be removed later
#endif
      real rdt
      real tz
      real akap
      integer ks
      save ks

      integer i, j, k, js2g0, js2gm1, jn2g0, count, ierror

      logical initialized
      data initialized /.false./

      real, allocatable, save :: sinp2(:),cosp2(:),cosp4(:), rf(:)

      js2g0  = max(2,jfirst)
      js2gm1 = max(2,jfirst+1)
      jn2g0  = min(jm-1,jlast)

      tz = 10.             ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24*3600
      rkv = 0.5*pdt/sday
      rka = pdt/ (40.*sday)      ! was 40 days
      rfc = 1./(1.+rka)
      rks = pdt/ (4.*sday)       ! was 4 days

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

      if ( .not. initialized ) then
        allocate( sinp2(jm), cosp2(jm), cosp4(jm), rf(km) )
        do j=2,jm-1
          sinp2(j) = sinp(j)**2
          cosp2(j) = cosp(j)**2
        enddo

        sinp2(1) = ( 0.5*(-1.+sine(2)) )**2
        sinp2(jm) = sinp2(1)
        cosp2(1) = ( 0.5*cose(2) ) **2
        cosp2(jm) = cosp2(1)

        do j=1,jm
          cosp4(j) = cosp2(j)**2
        enddo

        call set_eta(km, ks, ptop, ak, bk)

        if ( rayf ) then
          c1 = 1. / (36.*3600)
          pc = 1.
          do k=1,ks
            tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
            rf(k) = c1*(1.+tanh(log10(pc/tmp)))
            if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
	    rf(k) = 1./(1.+pdt*rf(k))
          enddo
        endif

	initialized = .true.
      endif

!$omp  parallel do                                      &
!$omp  default(shared)                                  &
!$omp  private(i,j,k,pl,tey,tez,dz,relx,dt_tropic, rdt) &
!$omp  private(teq,sigl,f1,rkt, tmp, u2, v2, fu, fv, t2, dp)

      do 1000 j=jfirst,jlast

        tey = ap0k*(315.-60.*sinp2(j))
!       tez = ap0k*10./akap*cosp2(j)
        tez = tz*(ap0k/akap)*cosp2(j)

        do k=1,km
          do i=1,im
            pl(i,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
          enddo
        enddo

        do k=km,1,-1
          do i=1,im
            if (strat .and. pl(i,k) < 10000.    &
                      .and. pl(i,k) > 100.  )  then
              dz = h0 * log(pl(i,k+1)/pl(i,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
              relx =  t_ms + tau*log(0.01*pl(i,k))
              relx = pdt/(relx*sday)
              dt_tropic = 2.25*cosp(j) * dz
              teq(i,k)  =  teq(i,k+1) + dt_tropic
              pt(i,j,k) = (pt(i,j,k)+relx*teq(i,k))/(1.+relx)
            elseif (strat .and. pl(i,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
              dz = h0 * log(pl(i,k+1)/pl(i,k))
              dt_tropic = -2.25*cosp(j) * dz
              teq(i,k) = teq(i,k+1) + dt_tropic
              pt(i,j,k) = (pt(i,j,k)+rms*teq(i,k))*rmr
            else
!
! Trop:  strictly Held-Suarez
!
              sigl = pl(i,k)/pe(i,km+1,j)
              f1 = max(0., (sigl-sigb) * rsgb )
              teq(i,k) = tey - tez*(log(pkz(i,j,k))+algpk)
              teq(i,k) = max(t0, teq(i,k)*pkz(i,j,k))
              rkt = rka + (rks-rka)*f1*cosp4(j)
              pt(i,j,k) = (pt(i,j,k)+rkt*teq(i,k))/(1.+rkt)
            endif
          enddo     !i-loop
        enddo     !k-loop

1000  continue

#if defined( SPMD )
!
! Communication might include ua and/or pe on the south only

      allocate (uasouth(im, km)) 
      allocate (pesouth(im, km+1)) 

1001  continue
      call mp_send3d(gid+1, gid-1, im, km+1, jm, 1, im, 1, km+1, jfirst, jlast,&
                            1, im, 1, km+1, jlast, jlast, pe)
      call mp_recv3d(gid-1, im, km+1, jm, 1, im, 1, km+1, jfirst-1, jfirst-1, &
                            1, im, 1, km+1, jfirst-1, jfirst-1, pesouth)
#endif

!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k,sigl,fac,frac)


      do 2000 k=1,km

        if (rayf .and. k <= ks) then
! Apply Rayleigh friction
          do j=js2g0,jlast
            do i=1,im
              u(i,j,k) = u(i,j,k)*rf(k)
            enddo
          enddo

          do j=js2g0,jn2g0
            do i=1,im
              v(i,j,k) = v(i,j,k)*rf(k)
            enddo
          enddo
        else
! Surface Rayleigh friction according to Held-Suarez
          do j=jfirst,jlast
            do i=1,im
              sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,km+1,j)
              frac(i,j) = max(0., (sigl-sigb)*rsgb )
            enddo
          enddo
#if defined( SPMD )
          if ( jfirst > 1 ) then
            do i=1,im
              sigl = 0.5*(pesouth(i,k)+pesouth(i,k+1)) / pesouth(i,km+1)
              frac(i,jfirst-1) = max(0., (sigl-sigb)*rsgb )
            enddo
          endif
#endif

! Backward adjustment
          do j=js2g0,jlast
            do i=1,im
              fac = frac(i,j)+frac(i,j-1)
              if (fac .gt. 0.) then
                u(i,j,k) = u(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo

          do j=js2g0,jn2g0
            do i=2,im
              fac = frac(i,j)+frac(i-1,j)
              if (fac .gt. 0.) then
                 v(i,j,k) = v(i,j,k)/(1.+rkv*fac)
              endif
            enddo
          enddo

          do j=js2g0,jn2g0
             fac = frac(1,j)+frac(im,j)
             if (fac .gt. 0.) then
                 v(1,j,k) = v(1,j,k)/(1.+rkv*fac)
             endif
          enddo
        endif
2000  continue


#if defined (SPMD)
      deallocate (pesouth)
      deallocate (uasouth)
#endif

      return
      end
