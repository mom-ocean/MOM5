module sw_core
 
! This module contains vertical independent part of the Lagrangian dynamics
! This is esentially shallow-water dynamics

  use tp_core,      only: tp2c, tpcc, xtp, ytp, ycc
#ifndef USE_LIMA
  use tp_core,      only: tp2d
#else
  use tp_core,      only: tp2d_lima
  use shr_kind_mod, only: r8 => shr_kind_r8
#endif
  implicit none

contains

#ifndef USE_LIMA
 subroutine c_sw(u,      v,       pt,     delp,                       &
                 uc,     vc,      ptc,    delpf,     ptk,             &
                 cosp,   acosp,   cose,   coslon,   sinlon,           &
                 dxdt,   dxe,     dtdx2,                              &
                 dtdx4,  dtxe5,   rdxe,   dycp,     dydt,             & 
                 dtdy5,  cye,     fc, zt_c,     tiny,     rcap,       &
                 im,     jm,      jfirst, jlast,    ng_c,             &
                 ng_d,   ng_s,    js2g0,  jn2g0,    js2gc,            &
                 jn1gc,  iord,    jord,   cosl5,    sinl5, dyn_step)
#else
 subroutine c_sw(u,      v,       pt,     delp,                       &
                 uc,     vc,      ptc,    delpf,     ptk,             &
                 cosp,   acosp,   cose,   coslon,   sinlon,           &
                 dxdt,   dxe,     dtdx2,                              &
                 dtdx4,  dtxe5,   rdxe,   dycp,     dydt,             &
                 dtdy5,  cye,     fc,     ifax,     trigs,            &
                 dc,     sc,      zt_c,   tiny,     rcap,             &
                 im,     jm,      jfirst, jlast,    ng_c,             &
                 ng_d,   ng_s,    js2g0,  jn2g0,    js2gc,            &
                 jn1gc,  iord,    jord,   cosl5,    sinl5  )
#endif
! Routine for shallow water dynamics on the C-grid

! !USES:

  use pft_module, only : pft2d


! INPUT:
  integer, intent(in):: im
  integer, intent(in):: jm
  integer, intent(in):: jfirst
  integer, intent(in):: jlast
  integer, intent(in):: js2g0
  integer, intent(in):: jn2g0
  integer, intent(in):: js2gc
  integer, intent(in):: jn1gc
  integer, intent(in):: iord
  integer, intent(in):: jord
  integer, intent(in):: ng_c
  integer, intent(in):: ng_s
  integer, intent(in):: ng_d
#ifndef USE_LIMA
  integer, intent(in):: dyn_step
#else
! polar filter related input arrays:
  integer, intent(in)::  ifax(13)                      !ECMWF fft
  real(r8), intent(in):: trigs(3*im/2+1)
  real(r8), intent(in):: dc(im,js2g0:jn2g0)
  real(r8), intent(in):: sc(js2g0:jn2g0)
#endif

! Prognostic variables:
  real, intent(in)   :: u(im,jfirst-ng_d:jlast+ng_s)
  real, intent(inout):: v(im,jfirst-ng_d:jlast+ng_d)  ! Wind in Y
  real, intent(in):: pt(im,jfirst-ng_d:jlast+ng_d) ! Wind in Y
  real, intent(in):: delp(im,jfirst:jlast)         ! Delta pressure
  real, intent(in):: delpf(im,jfirst-ng_d:jlast+ng_d)

#ifndef USE_LIMA
  real, intent(in):: fc(im,js2gc:jn1gc)
#else               
  real, intent(in):: fc(js2gc:jn1gc)
#endif

  real, intent(in):: cosp(jm)
  real, intent(in):: acosp(jm)
  real, intent(in):: cose(jm)

  real, intent(in):: dxdt(jm)
  real, intent(in):: dxe(jm)
  real, intent(in):: rdxe(jm)
  real, intent(in):: dtdx2(jm)
  real, intent(in):: dtdx4(jm)
  real, intent(in):: dtxe5(jm)
  real, intent(in):: dycp(jm)
  real, intent(in)::  cye(jm)

  real, intent(in):: sinlon(im)
  real, intent(in):: coslon(im)
  real, intent(in):: sinl5(im)
  real, intent(in):: cosl5(im)

  real, intent(in):: zt_c
  real, intent(in):: rcap
  real, intent(in):: tiny
  real, intent(in):: dydt
  real, intent(in):: dtdy5

! Output:
  real, intent(out):: uc(im,jfirst-ng_d:jlast+ng_d)
  real, intent(out):: vc(im,jfirst-2   :jlast+2 )
  real, intent(out):: ptc(im,jfirst:jlast)
  real, intent(out):: ptk(im,jfirst:jlast)

!--------------------------------------------------------------
! Local 

    real    fx(im,jfirst:jlast)
    real   xfx(im,jfirst:jlast)
    real   tm2(im,jfirst:jlast)

    real    va(im,jfirst-1:jlast)

    real   wk1(im,jfirst-1:jlast+1)
    real   cry(im,jfirst-1:jlast+1)
    real    fy(im,jfirst-1:jlast+1)

    real  ymass(im,jfirst: jlast+1) 
    real    yfx(im,jfirst: jlast+1)

    real   crx(im,jfirst-ng_c:jlast+ng_c)
    real    u2(im,jfirst-ng_d:jlast+ng_d)
    real    v2(im,jfirst-ng_d:jlast+ng_d)

    real  fxj(im)
    real  p1d(im)
    real  cx1(im)

    real  qtmp(-im/3:im+im/3)
    real  slope(-im/3:im+im/3)
    real  al(-im/3:im+im/3)
    real  ar(-im/3:im+im/3)
    real  a6(-im/3:im+im/3)

    real  us, vs, un, vn
    real  p1ke, p2ke

    logical ffsl(jm)
    logical sld

    integer i, j, im2
    integer js1g1
    integer js2g1
    integer js2gc1
    integer js2gcp1
    integer jn2gc
    integer jn1g1

    im2 = im/2
#ifdef USE_LIMA
    ffsl = .true.
#endif

! Set loop limits

    js1g1 = max(1,jfirst-1)
    js2g1 = max(2,jfirst-1)
    js2gcp1 = max(2,jfirst-ng_c-1)   ! NG-1 latitudes on S (starting at 2)
    jn1g1 = min(jm,jlast+1)
    jn2gc = min(jm-1,jlast+ng_c)   ! NG latitudes on N (ending at jm-1)
 
!
! Treat the special case of ng_c = 1
!
    if ( ng_c == 1 .AND. ng_d > 1 ) THEN
        js2gc1 = js2gc
    else
        js2gc1 = max(2,jfirst-ng_c+1)   ! NG-1 latitudes on S (starting at 2)
    endif

! Get D-grid V-wind at the poles.
    call vpol5(u, v, im, jm,            &
               coslon, sinlon, cosl5, sinl5, ng_d, ng_s, jfirst, jlast )

    do j=max(2,jfirst-ng_d),  min(jm-1,jlast+ng_d)
       do i=1,im-1
          v2(i,j) = v(i,j) + v(i+1,j)
       enddo
          v2(im,j) = v(im,j) + v(1,j)
    enddo

    do j=max(2,jfirst-ng_d), min(jm-1,jlast+ng_s-1)
       do i=1,im
          u2(i,j) = u(i,j) + u(i,j+1)
       enddo
    enddo

       if ( jfirst-ng_d <= 1 ) then
! Projection at SP
          us = 0.
          vs = 0.

          do i=1,im2
            us = us + (u2(i+im2,2)-u2(i,2))*sinlon(i)         &
                    + (v2(i,2)-v2(i+im2,2))*coslon(i)
            vs = vs + (u2(i+im2,2)-u2(i,2))*coslon(i)         &
                    + (v2(i+im2,2)-v2(i,2))*sinlon(i)
          enddo

          us = us/im
          vs = vs/im

! SP
          do i=1,im2
            u2(i,1)  = -us*sinlon(i) - vs*coslon(i)
            v2(i,1)  =  us*coslon(i) - vs*sinlon(i)
            u2(i+im2,1)  = -u2(i,1)
            v2(i+im2,1)  = -v2(i,1)
          enddo

          p1ke = 0.125*(u2(1, 1)**2 + v2(1, 1)**2)

        endif

        if ( jlast+ng_d >= jm ) then

! Projection at NP
          un = 0.
          vn = 0.

          j = jm-1
          do i=1,im2
            un = un + (u2(i+im2,j)-u2(i,j))*sinlon(i)        &
                    + (v2(i+im2,j)-v2(i,j))*coslon(i)
            vn = vn + (u2(i,j)-u2(i+im2,j))*coslon(i)        &
                    + (v2(i+im2,j)-v2(i,j))*sinlon(i)
          enddo

          un = un/im
          vn = vn/im

! NP
          do i=1,im2
            u2(i,jm) = -un*sinlon(i) + vn*coslon(i)
            v2(i,jm) = -un*coslon(i) - vn*sinlon(i)
            u2(i+im2,jm) = -u2(i,jm)
            v2(i+im2,jm) = -v2(i,jm)
          enddo

          p2ke = 0.125*(u2(1,jm)**2 + v2(1,jm)**2)

        endif

! A -> C
        do j=js2gc,jn2gc                ! uc needed N*ng S*ng
!       do j=max(2,jfirst-2),min(jm-1,jlast+2)
! i=1
            uc(1,j) = 0.25*(u2(1,j)+u2(im,j))
          do i=2,im
            uc(i,j) = 0.25*(u2(i,j)+u2(i-1,j))
          enddo
        enddo

        do j=max(2,jfirst-2), min(jm,jlast+2) ! ng_d must be 3
          do i=1,im
            vc(i,j) = 0.25*(v2(i,j)+v2(i,j-1))  ! v2 needed (-3:+2)
          enddo
        enddo

        if ( jfirst /= 1 ) then
          do i=1,im
            cry(i,jfirst-1) = dtdy5*vc(i,jfirst-1)
          enddo
        endif

        do j=js2g0,jn1g1                     ! ymass needed on NS
          do i=1,im
               cry(i,j) = dtdy5*vc(i,j)
             ymass(i,j) = cry(i,j)*cose(j)
          enddo
        enddo

! New va definition
        do j=js2g1,jn2g0                     ! va needed on S (for YCC, iv==1)
          do i=1,im
            va(i,j) = 0.5*(cry(i,j)+cry(i,j+1))
          enddo
        enddo

! SJL: Check if FFSL integer fluxes need to be computed

        do 2222 j=js2gc,jn2gc                ! ffsl needed on N*sg S*sg
          do i=1,im
            crx(i,j) = uc(i,j)*dtdx2(j)
          enddo
          ffsl(j) = .false.
          if( cosp(j) < zt_c ) then
            do i=1,im
              if( abs(crx(i,j)) > 1. ) then
                ffsl(j) = .true. 
                go to 2222
              endif
            enddo
          endif
2222    continue

! 2D transport of polar filtered delp (for computing fluxes!)
! Update is done on the unfiltered delp

#ifndef USE_LIMA
      call tp2c( ptk,  delpf(1,jfirst-ng_c),     &
              crx(1,jfirst-ng_c), cry(1,jfirst),             &
              im, jm, iord, jord, ng_c, xfx,                 &
              yfx, ffsl, rcap, acosp,                        &
              crx(1,jfirst), ymass, cosp,                    &
              0, jfirst, jlast, delp, dyn_step)
#else
      call tp2c( ptk,  va(1,jfirst),  delpf(1,jfirst-ng_c),     &
              crx(1,jfirst-ng_c), cry(1,jfirst),             &
              im, jm, iord, jord, ng_c, xfx,                 &
              yfx, ffsl, rcap, acosp,                        &
              crx(1,jfirst), ymass, cosp,                    &
              0, jfirst, jlast)   
#endif

   do j=js2g0,jn2g0                      ! xfx not ghosted
      if( ffsl(j) ) then
         do i=1,im
           xfx(i,j) = xfx(i,j)/sign(max(abs(crx(i,j)),tiny),crx(i,j))
         enddo
      endif
   enddo

! pt-advection using pre-computed mass fluxes
! tm2 below as the storage for pt increment
! WS 99.09.20 : pt, crx need on N*ng S*ng, yfx on N

#if !defined (SW_DYN)

#ifndef USE_LIMA
    call tp2c(tm2, pt(1,jfirst-ng_c),       &
              crx(1,jfirst-ng_c), cry(1,jfirst),          &
              im, jm,  iord, jord, ng_c, fx,              &
              fy(1,jfirst), ffsl, rcap, acosp,            &
              xfx, yfx, cosp, 1, jfirst, jlast, delp, dyn_step)
#else
    call tp2c(tm2 ,va(1,jfirst), pt(1,jfirst-ng_c),       &
              crx(1,jfirst-ng_c), cry(1,jfirst),          &
              im, jm,  iord, jord, ng_c, fx,              &
              fy(1,jfirst), ffsl, rcap, acosp,            &
              xfx, yfx, cosp, 1, jfirst, jlast)

#endif

#endif

#if ( !defined ALT_PFT )
! v2, crx as work arrays
#ifndef USE_LIMA
         call pft2d(ptk(1,js2g0), im, jn2g0-js2g0+1, v2, crx, 1)
#else
         call pft2d(ptk(1,js2g0), sc(js2g0), dc(1,js2g0), im,   &
                jn2g0-js2g0+1,  ifax, trigs, v2, crx )
#endif

#if !defined (SW_DYN)

#ifndef USE_LIMA
         call pft2d(tm2(1,js2g0), im, jn2g0-js2g0+1, v2, crx, 1)
#else
         call pft2d(tm2(1,js2g0), sc(js2g0), dc(1,js2g0),  im,  &
                jn2g0-js2g0+1,  ifax, trigs, v2, crx )
#endif

#endif

#endif

    do j=jfirst,jlast
       do i=1,im
          ptk(i,j) = delp(i,j) + ptk(i,j)
#ifdef SW_DYN
          ptc(i,j) = pt(i,j)
#else
          ptc(i,j) = (pt(i,j)*delp(i,j) + tm2(i,j))/ptk(i,j)
#endif
       enddo
    enddo

!------------------
! Momentum equation
!------------------
     call ycc(im, jm, fy, vc(1,jfirst-2), va(1,jfirst-1),   &
              va(1,jfirst-1), jord, 1, jfirst, jlast)

     do j=js2g1,jn2g0

          do i=1,im
            cx1(i) = dtdx4(j)*u2(i,j)
          enddo

          sld = .false.
          if( cosp(j) < zt_c ) then
            do i=1,im
              if( abs(cx1(i)) > 1. ) then
                sld = .true. 
                go to 3333
              endif
            enddo
          endif
3333      continue

          p1d(im) = uc(1,j)
          do i=1,im-1
            p1d(i) = uc(i+1,j)
          enddo

          call xtp(im,   sld, fxj, p1d, cx1, iord,   &
                   cx1,  cosp(j),  0,   slope,       &
                   qtmp, al,       ar,  a6)
 
          do i=1,im
            wk1(i,j) = dxdt(j)*fxj(i) + dydt*fy(i,j)
          enddo
     enddo


     if ( jfirst-1 <= 1 ) then
          do i=1,im
            wk1(i,1) = p1ke
          enddo
     endif

     if ( jlast+1 >= jm ) then
          do i=1,im
            wk1(i,jm) = p2ke
          enddo
     endif

! crx redefined
     do j=js2gc1,jn1gc
            crx(1,j) = dtxe5(j)*u(im,j)
          do i=2,im
            crx(i,j) = dtxe5(j)*u(i-1,j)
          enddo
     enddo

     if ( jfirst /=1 ) then 
          do i=1,im
             cry(i,jfirst-1) = dtdy5*v(i,jfirst-1)
          enddo
     endif

     do j=jfirst,jlast
        do i=1,im
             cry(i,j) = dtdy5*v(i,j)
           ymass(i,j) = cry(i,j)*cosp(j)       ! ymass actually unghosted
        enddo
     enddo

     do j=js2g0,jlast
          do i=1,im
            tm2(i,j) = 0.5*(cry(i,j)+cry(i,j-1)) ! cry ghosted on S 
          enddo
     enddo

!    Compute absolute vorticity on the C-grid.

     if ( jfirst-ng_d <= 1 ) then
          do i=1,im
            u2(i,1) = 0.
          enddo
     endif

     do j=js2gc,jn2gc
         do i=1,im
            u2(i,j) = uc(i,j)*cosp(j)
         enddo
     enddo

     if ( jlast+ng_d >= jm ) then
          do i=1,im
            u2(i,jm) = 0.
          enddo
     endif


#ifndef USE_LIMA
     do j=js2gc1,jn1gc
! The computed absolute vorticity on C-Grid is assigned to v2
          v2(1,j) = fc(1,j) + (u2(1,j-1)-u2(1,j))*cye(j) +     &
                    (vc(1,j) - vc(im,j))*rdxe(j)

          do i=2,im
             v2(i,j) = fc(i,j) + (u2(i,j-1)-u2(i,j))*cye(j) +  &
                       (vc(i,j) - vc(i-1,j))*rdxe(j)
          enddo
     enddo
#else               
     do j=js2gc1,jn1gc
! The computed absolute vorticity on C-Grid is assigned to v2
          v2(1,j) = fc(j) + (u2(1,j-1)-u2(1,j))*cye(j) +     &
                    (vc(1,j) - vc(im,j))*rdxe(j)

          do i=2,im
             v2(i,j) = fc(j) + (u2(i,j-1)-u2(i,j))*cye(j) +  &
                       (vc(i,j) - vc(i-1,j))*rdxe(j)
          enddo
     enddo
#endif

!         ffsl = .true.
     do 2233 j=js2gc1,jn1gc          ! ffsl needed on N*ng S*(ng-1)
          ffsl(j) = .false.
          if( cose(j) < zt_c ) then
            do i=1,im
              if( abs(crx(i,j)) > 1. ) then
                ffsl(j) = .true. 
                go to 2233
              endif
            enddo
          endif
2233    continue

   call tpcc( tm2, ymass, v2(1,jfirst-ng_d), crx(1,jfirst-ng_c),  &
              cry(1,jfirst), im, jm, ng_c, ng_d,                  &
              iord, jord, fx, fy(1,jfirst), ffsl, cose,           &
              jfirst, jlast, slope, qtmp, al, ar, a6 )

   do j=js2g0,jn2g0
         uc(1,j) = uc(1,j) + dtdx2(j)*(wk1(im,j)-wk1(1,j)) + dycp(j)*fy(1,j)
      do i=2,im
         uc(i,j) = uc(i,j) + dtdx2(j)*(wk1(i-1,j)-wk1(i,j)) + dycp(j)*fy(i,j)
      enddo
   enddo

   do j=js2g0,jlast
        do i=1,im-1
           vc(i,j) = vc(i,j) + dtdy5*(wk1(i,j-1)-wk1(i,j))-dxe(j)*fx(i+1,j)
        enddo
           vc(im,j) = vc(im,j) + dtdy5*(wk1(im,j-1)-wk1(im,j))-dxe(j)*fx(1,j)
   enddo

 end subroutine c_sw


!--------------------------------------------------------------------------
 subroutine d_sw( u,     v,      uc,       vc,                       &   
                  pt,    delp,   delpf,    cx3,                      &
                  cy3,    mfx,    mfy,     cdx,    cdy,              &
                  dtdx,  dtdxe,  dtxe5,  txe5,  dyce,   rdx,         &
                  cy,    dx,    f0,     js2g0,  jn1g1,               &  
                  im,    jm,     jfirst, jlast, ng_d,   ng_s,        &
                  nq,    iord,   jord,   sord,  zt_d,  rcap,         &
                  tiny,  dtdy,  dtdy5,  tdy5,   rdy,   cosp,         &
#ifndef USE_LIMA
                  acosp, cose,  coslon, sinlon, cosl5, sinl5, dyn_step )
#else
                  acosp, cose,  coslon, sinlon, cosl5, sinl5 )
#endif
!--------------------------------------------------------------------------
! Routine for shallow water dynamics on the D-grid

! INPUT:
  integer, intent(in):: im
  integer, intent(in):: jm
  integer, intent(in):: jfirst
  integer, intent(in):: jlast
  integer, intent(in):: iord
  integer, intent(in):: jord
  integer, intent(in):: sord               ! for pt and delp
  integer, intent(in):: js2g0
  integer, intent(in):: jn1g1
  integer, intent(in):: ng_d
  integer, intent(in):: ng_s
  integer, intent(in):: nq
#ifndef USE_LIMA
  integer, intent(in):: dyn_step
#endif

! Prognostic variables:
  real, intent(in):: u(im,jfirst-ng_d:jlast+ng_s)          ! Wind in X
  real, intent(in):: v(im,jfirst-ng_d:jlast+ng_d)          ! Wind in Y
  real, intent(inout):: delp(im,jfirst:jlast)         ! Delta pressure
  real, intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d) ! (Virtual) Potential temperature
  real, intent(inout):: delpf(im,jfirst-ng_d:jlast+ng_d)

  real, intent(in):: cosp(jm)
  real, intent(in):: acosp(jm)
  real, intent(in):: cose(jm)

  real, intent(in):: sinlon(im)
  real, intent(in):: coslon(im)
  real, intent(in):: sinl5(im)
  real, intent(in):: cosl5(im)

  real, intent(in):: dtdx(jm)
  real, intent(in):: dtdxe(jm)
  real, intent(in):: dx(jm)
  real, intent(in):: rdx(jm)
  real, intent(in):: cy(jm)
  real, intent(in):: dyce(jm)
  real, intent(in):: dtxe5(jm)
  real, intent(in):: txe5(jm)

  real, intent(in):: cdx(js2g0:jn1g1)
  real, intent(in):: cdy(js2g0:jn1g1)
#ifndef USE_LIMA
  real, intent(in):: f0(im,jfirst-ng_d:jlast+ng_d)
#else               
  real, intent(in):: f0(jfirst-ng_d:jlast+ng_d)
#endif

  real, intent(in):: zt_d
  real, intent(in):: rcap
  real, intent(in):: tiny
  real, intent(in):: tdy5
  real, intent(in):: rdy
  real, intent(in):: dtdy
  real, intent(in):: dtdy5

! INPUT/OUTPUT:
  real, intent(inout):: uc(im,jfirst-ng_d:jlast+ng_d)
  real, intent(inout):: vc(im,jfirst-1   :jlast+2 )
  real, intent(inout):: cx3(im,jfirst-ng_d:jlast+ng_d)! Accumulated Courant number in X
  real, intent(inout):: cy3(im,jfirst:jlast+1)        ! Accumulated Courant number in Y
  real, intent(inout):: mfx(im,jfirst:jlast)          ! Mass flux in X  (unghosted)
  real, intent(inout):: mfy(im,jfirst:jlast+1)        ! Mass flux in Y

! Local 
    real    fx(im,jfirst:jlast)
    real   xfx(im,jfirst:jlast)

    real   wk1(im,jfirst-1:jlast+1)
    real   cry(im,jfirst-1:jlast+1)
    real    fy(im,jfirst-1:jlast+1)

    real  ymass(im,jfirst: jlast+1) 
    real    yfx(im,jfirst: jlast+1)

    real    ub(im,jfirst:  jlast+1)
#ifdef USE_LIMA
    real(r8)   va(im,jfirst-1:jlast)
#endif

    real   crx(im,jfirst-ng_d:jlast+ng_d)

    real  fxj(im)
    real  qtmp(-im/3:im+im/3)
    real  slope(-im/3:im+im/3)
    real  al(-im/3:im+im/3)
    real  ar(-im/3:im+im/3)
    real  a6(-im/3:im+im/3)

    real   c1, c2

    logical ffsl(jm)
    logical sld

    integer i, j
    integer js2gd, jn2g0, jn2g1, jn2gd, jn1gd

! Set loop limits

  jn2g0 = min(jm-1,jlast)
  jn2g1 = min(jm-1,jlast+1)
  js2gd = max(2,jfirst-ng_d)     ! NG latitudes on S (starting at 1)
  jn2gd = min(jm-1,jlast+ng_d)   ! NG latitudes on S (ending at jm-1)
  jn1gd = min(jm,jlast+ng_d)     ! NG latitudes on N (ending at jm)


! Get C-grid U-wind at poles.
! call upol5(uc,   vc,   im,  jm,  coslon,  sinlon,  &
!            cosl5, sinl5, jfirst, jlast, ng_d)

  do j=js2gd,jn2gd                     ! crx needed on N*ng S*ng
     do i=1,im
        crx(i,j) = dtdx(j)*uc(i,j)
     enddo
  enddo

  do 2225 j=js2gd,jn2gd                ! ffsl needed on N*ng S*ng
     ffsl(j) = .false.
     if( cosp(j) < zt_d ) then
         do i=1,im
            if( abs(crx(i,j)) > 1. ) then
               ffsl(j) = .true. 
               go to 2225
            endif
         enddo
      endif
2225    continue

  do j=js2g0,jn1g1                       ! cry, ymass needed on N
     do i=1,im
        cry(i,j) = dtdy*vc(i,j)
        ymass(i,j) = cry(i,j)*cose(j)
     enddo
  enddo

#ifdef USE_LIMA
  do j=js2g0,jn2g0                         ! No ghosting
     do i=1,im
        if( cry(i,j)*cry(i,j+1) > 0. ) then
           if( cry(i,j) > 0. ) then
              va(i,j) = cry(i,j)
           else
              va(i,j) = cry(i,j+1)         ! cry ghosted on N
           endif
        else
           va(i,j) = 0.
        endif
     enddo
  enddo
#endif


! transport polar filtered delp
#ifndef USE_LIMA
      call tp2c(ub(1,jfirst), delpf(1,jfirst-ng_d),   &
                crx(1,jfirst-ng_d),cry(1,jfirst),im,jm,sord,sord,   &
                ng_d, xfx, yfx, ffsl,                               &
                rcap, acosp,crx(1,jfirst), ymass,                   &
                cosp, 0, jfirst, jlast, delp, dyn_step)
#else
      call tp2c(ub(1,jfirst), va(1,jfirst), delpf(1,jfirst-ng_d),   &
                crx(1,jfirst-ng_d),cry(1,jfirst),im,jm,sord,sord,   &
                ng_d, xfx, yfx, ffsl,                               &
                rcap, acosp,crx(1,jfirst), ymass,                   &
                cosp, 0, jfirst, jlast)
#endif


! <<< Save necessary data for large time step tracer transport >>>
      if( nq > 0 ) then
          do j=js2g0,jn2g0                       ! No ghosting needed
  do i=1,im
              cx3(i,j) = cx3(i,j) + crx(i,j)
              mfx(i,j) = mfx(i,j) + xfx(i,j)
            enddo
          enddo

          do j=js2g0,jlast                      ! No ghosting needed
            do i=1,im
              cy3(i,j) = cy3(i,j) + cry(i,j)
              mfy(i,j) = mfy(i,j) + yfx(i,j)
            enddo
          enddo
      endif

     do j=js2g0,jn2g0                         ! No ghosting needed
        if( ffsl(j) ) then
          do i=1,im
             xfx(i,j) = xfx(i,j)/sign(max(abs(crx(i,j)),tiny),crx(i,j))
          enddo
        endif
     enddo

! Update delp
        do j=jfirst,jlast
          do i=1,im
! SAVE old delp: pressure thickness ~ "air density"
            wk1(i,j) = delp(i,j)
            delp(i,j) = wk1(i,j) + ub(i,j)
          enddo
        enddo

#if !defined ( SW_DYN )
! pt Advection
#ifndef USE_LIMA
  call tp2c(ub(1,jfirst), pt(1,jfirst-ng_d),    &
            crx(1,jfirst-ng_d),cry(1,jfirst),               &
            im,jm, sord, sord, ng_d, fx, fy(1,jfirst),      &
            ffsl, rcap, acosp, xfx, yfx(1,jfirst), cosp,    &
            1, jfirst,jlast, wk1(1,jfirst), dyn_step)
#else
  call tp2c(ub(1,jfirst),va(1,jfirst),pt(1,jfirst-ng_d),    &
            crx(1,jfirst-ng_d),cry(1,jfirst),               &
            im,jm, sord, sord, ng_d, fx, fy(1,jfirst),      &
            ffsl, rcap, acosp,                              &
            xfx, yfx(1,jfirst), cosp, 1, jfirst,jlast)
#endif

! Update pt.
      do j=jfirst,jlast
         do i=1,im
            pt(i,j) = (pt(i,j)*wk1(i,j)+ub(i,j)) / delp(i,j)
         enddo
      enddo
#endif

#ifdef ALT_KE
!----------------------------------------------------
! Computing KE directly on the corners using C-winds
!----------------------------------------------------
  do j=max(1,jfirst-1),jn1g1
     do i=1,im
        fy(i,j) = uc(i,j)**2 
     enddo
  enddo
 
  do j=js2g0,jn1g1 
     do i=1,im
        fxj(i) = vc(i,j)**2 
     enddo
        wk1(1,j) = fxj(im) + fxj(1)
     do i=2,im
        wk1(i,j) = fxj(i-1) + fxj(i)
     enddo
     do i=1,im
        wk1(i,j) = 0.25*(wk1(i,j) + fy(i,j) + fy(i,j-1))
     enddo
  enddo
#else
! Start using ub as v (CFL) on B-grid (cell corners)
      do j=js2g0,jn1g1                          ! ub needed on N
           ub(1,j) = dtdy5*(vc(1,j) + vc(im,j))  
         do i=2,im
            ub(i,j) = dtdy5*(vc(i,j) + vc(i-1,j))
         enddo
      enddo

      call ytp(im, jm, fy(1,jfirst), v(1,jfirst-ng_d), ub(1,jfirst),  &
               ub(1,jfirst), ng_d, jord, 1, jfirst, jlast)
! End using ub as v (CFL) on B-grid

   do j=js2g0,jn1g1                 ! ub needed on N
       do i=1,im                
          ub(i,j) = dtxe5(j)*(uc(i,j) + uc(i,j-1))
!                        uc will be used as wrok array after this point
       enddo
   enddo

  do j=js2g0,jn1g1                       ! wk1 needed on N
          sld = .false.
          if( cose(j) < zt_d ) then
            do i=1,im
              if( abs(ub(i,j)) > 1. ) then    ! ub ghosted on N
                sld = .true. 
                go to 2235
              endif
            enddo
          endif
2235      continue

     call xtp(im,  sld, fxj, u(1,j), ub(1,j),          &
              iord, ub(1,j), cose(j), 0,               &
              slope, qtmp, al, ar, a6 )

     do i=1,im
        wk1(i,j) =  txe5(j)*fxj(i) + tdy5*fy(i,j)  ! fy ghosted on N
     enddo
  enddo
#endif

! Add divergence damping to vector invariant form of the momentum eqn
! (absolute vorticity is damped by ffsl scheme, therefore divergence damping
! provides more consistent dissipation to divergent part of the flow)

!--------------------------
! Perform divergence damping 
!--------------------------

        do j=max(2,jfirst-1), jn2g1                   ! fy need on NS (below)
            do i=1,im
              fy(i,j) = v(i,j)*cosp(j)      ! v ghosted on NS at least
            enddo
        enddo

        do j=js2g0,jn1g1
! i=1
              uc(1,j) = u(im,j) - u(1,j)    ! u ghosted on N at least
            do i=2,im
              uc(i,j) = u(i-1,j) - u(i,j)
            enddo
        enddo

      if ( jfirst-1 <= 2 ) then
! j=2
           do i=1,im
              wk1(i,2) = wk1(i,2) - cdy(2)*fy(i, 2) + cdx(2)*uc(i,2)
           enddo
      endif

        do j=max(3,jfirst),jn2g1
            do i=1,im
              wk1(i,j) = wk1(i,j) + cdy(j)*(fy(i,j-1) - fy(i,j))  &
                                  + cdx(j)*uc(i,j)
            enddo
        enddo

      if ( jlast+1 >= jm ) then
           do i=1,im
              wk1(i,jm) = wk1(i,jm) + cdy(jm)*fy(i,jm-1) + cdx(jm)*uc(i,jm)
           enddo
      endif
!------------------------------------
! End divergence damping computation
!------------------------------------


! Compute Vorticity on the D grid
! delpf used as work array

      do j=js2gd,jn1gd
         do i=1,im
            delpf(i,j) = u(i,j)*cose(j)   ! u ghosted on N*ng S*ng
         enddo
      enddo

      if ( jfirst-ng_d <= 1 ) then
          c1 = 0.
          do i=1,im
            c1 = c1 + delpf(i,2)
          end do
          c1 = -c1*rdy*rcap

          do i=1,im
            uc(i,1) = c1
          enddo
      endif

      if ( jlast+ng_d >= jm ) then
          c2 = 0.
          do i=1,im
            c2 = c2 + delpf(i,jm)
          end do
          c2 = c2*rdy*rcap

          do i=1,im
            uc(i,jm) = c2
          enddo
      else

! This is an attempt to avoid ghosting u on N*(ng+1)
          do i=1,im
             uc(i,jn2gd) = 1.e40
          enddo
      endif

      do j=js2gd, min(jm-1,jlast+ng_d-1)
          do i=1,im-1
             uc(i,j) = ( delpf(i,j) - delpf(i,j+1)) * cy(j)  +         &
                        (v(i+1,j) - v(i,j))    * rdx(j)
          enddo
            uc(im,j) = (delpf(im,j) - delpf(im,j+1)) *  cy(j) +        &
                       (v(1,j) - v(im,j)) * rdx(j)
      enddo

! uc is relative vorticity at this point

#ifndef USE_LIMA
      do j=max(1,jfirst-ng_d), jn1gd
          do i=1,im
             uc(i,j) = uc(i,j) + f0(i,j)
! uc is absolute vorticity
          enddo
      enddo
#else               
      do j=max(1,jfirst-ng_d), jn1gd
          do i=1,im
             uc(i,j) = uc(i,j) + f0(j)
! uc is absolute vorticity
          enddo
      enddo
#endif

!-----------------------------
! Transport absolute vorticity
!-----------------------------
#ifndef USE_LIMA
      call tp2d(uc(1,jfirst-ng_d), crx(1,jfirst-ng_d),   &
                cry(1,jfirst), im, jm, iord, jord, ng_d, fx,           &
                fy(1,jfirst), ffsl, crx(1,jfirst),                     &
                ymass, cosp, acosp, 0, jfirst, jlast, delp, dyn_step)

#else
      call tp2d_lima(va(1,jfirst), uc(1,jfirst-ng_d), crx(1,jfirst-ng_d),   &
                cry(1,jfirst), im, jm, iord, jord, ng_d, fx,           &
                fy(1,jfirst), ffsl, crx(1,jfirst),                     &
                ymass, cosp, 0, jfirst, jlast)

#endif

      do j=js2g0,jlast
          do i=1,im-1
            uc(i,j) = dtdxe(j)*(wk1(i,j)-wk1(i+1,j)) + dyce(j)*fy(i,j)
          enddo
           uc(im,j) = dtdxe(j)*(wk1(im,j)-wk1(1,j)) + dyce(j)*fy(im,j)
      enddo

      do j=js2g0,jn2g0
          do i=1,im
            vc(i,j) = dtdy*(wk1(i,j)-wk1(i,j+1)) - dx(j)*fx(i,j)
          enddo
      enddo

 end subroutine d_sw


 subroutine upol5(u, v, im, jm, coslon, sinlon, cosl5, sinl5,    &
                  jfirst, jlast, ng)

! !INPUT PARAMETERS:
      integer im                       ! Total longitudes
      integer jm                       ! Total latitudes
      integer jfirst                   ! First PE latitude (no ghosting)
      integer jlast                    ! Last  PE latitude (no ghosting)
      integer ng

      real, intent(in):: coslon(im), sinlon(im)
      real, intent(in):: cosl5(im),  sinl5(im)
      real, intent(in):: v(im,jfirst:jlast)   ! Winds in Y (C-grid)
      real, intent(inout):: u(im,jfirst:jlast)      ! Winds in X (C-grid)

! LOCAL VARIABLES:

      integer i, imh
      real  uanp(im), uasp(im), vanp(im), vasp(im)
      real  un, vn, us, vs, r2im

      imh = im/2
      r2im = 0.5d0/dble(im)

      if ( jfirst == 1 ) then
!
! Treat SP
!
      do i=1,im-1
         uasp(i) = u(i,  2) + u(i+1,2)
      enddo
         uasp(im) = u(im,  2) + u(1,2)

      do i=1,im
         vasp(i) = v(i,  2) + v(i,  3)
      enddo

! Projection at SP

      us = 0.; vs = 0.
      do i=1,imh
         us = us + (uasp(i+imh)-uasp(i))*sinlon(i)     &
                 + (vasp(i)-vasp(i+imh))*coslon(i)
         vs = vs + (uasp(i+imh)-uasp(i))*coslon(i)     &
                 + (vasp(i+imh)-vasp(i))*sinlon(i)
      enddo
      us = us*r2im
      vs = vs*r2im

! get U-wind at SP

      do i=1,imh
         u(i,  1) = -us*sinl5(i) - vs*cosl5(i)
         u(i+imh,  1) = -u(i,  1)
      enddo

      endif

      if ( jlast == jm ) then
!
! Treat NP
!
      do i=1,im-1
         uanp(i) = u(i,jm-1) + u(i+1,jm-1)
      enddo
        uanp(im) = u(im,jm-1) + u(1,jm-1)

      do i=1,im
         vanp(i) = v(i,jm-1) + v(i,jm)
      enddo

! Projection at NP

      un = 0.;  vn = 0.
      do i=1,imh
         un = un + (uanp(i+imh)-uanp(i))*sinlon(i)  &
                 + (vanp(i+imh)-vanp(i))*coslon(i)
         vn = vn + (uanp(i)-uanp(i+imh))*coslon(i)  &
                 + (vanp(i+imh)-vanp(i))*sinlon(i)
      enddo
      un = un*r2im
      vn = vn*r2im

! get U-wind at NP
      do i=1,imh
         u(i,jm) = -un*sinl5(i) + vn*cosl5(i)
         u(i+imh,jm) = -u(i,jm)
      enddo

      endif

!EOC
 end subroutine upol5
!----------------------------------------------------------------------- 

!----------------------------------------------------------------------- 
!BOP
!
 subroutine vpol5(u, v, im, jm, coslon, sinlon, cosl5, sinl5,    &
                  ng_d,  ng_s,  jfirst, jlast)

! !INPUT PARAMETERS:
      integer im                       ! Total longitudes
      integer jm                       ! Total latitudes
      integer jfirst                   ! First PE latitude (no ghosting)
      integer jlast                    ! Last  PE latitude (no ghosting)
      integer, intent(in):: ng_s, ng_d
      real, intent(in):: coslon(im), sinlon(im)
      real, intent(in):: cosl5(im),sinl5(im)
      real, intent(in):: u(im,jfirst-ng_d:jlast+ng_s)

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: v(im,jfirst-ng_d:jlast+ng_d)

! !DESCRIPTION:
!
!   Treat the V winds at the poles.  This requires an average 
!   of the U- and V-winds, weighted by their angles of incidence
!   at the pole points.     
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer i, imh
      real  uanp(im), uasp(im), vanp(im), vasp(im)
      real  un, vn, us, vs, r2im

! WS 99.05.25 :  Replaced conversions of IMR with IM
      r2im = 0.5d0/dble(im)
      imh  = im / 2

! WS 990726 :  Added condition to decide if poles are on this processor

   if ( jfirst-ng_d <= 1 ) then
         do i=1,im
            uasp(i) = u(i,  2) + u(i,3)
         enddo

         do i=1,im-1
            vasp(i)  = v(i,  2) + v(i+1,2)
         enddo
            vasp(im) = v(im,2) + v(1,2)

! Projection at SP
      us = 0.; vs = 0.

      do i=1,imh
         us = us + (uasp(i+imh)-uasp(i))*sinlon(i)    &
                 + (vasp(i)-vasp(i+imh))*coslon(i)
         vs = vs + (uasp(i+imh)-uasp(i))*coslon(i)    &
                 + (vasp(i+imh)-vasp(i))*sinlon(i)
      enddo
      us = us*r2im
      vs = vs*r2im

! get V-wind at SP

      do i=1,imh
         v(i,    1) =  us*cosl5(i) - vs*sinl5(i)
         v(i+imh,1) = -v(i,1)
      enddo

   endif

   if ( jlast+ng_d >= jm ) then

      do i=1,im
         uanp(i) = u(i,jm-1) + u(i,jm)
      enddo

      do i=1,im-1
         vanp(i) = v(i,jm-1) + v(i+1,jm-1)
      enddo
         vanp(im) = v(im,jm-1) + v(1,jm-1)

! Projection at NP

      un = 0.
      vn = 0.
      do i=1,imh
         un = un + (uanp(i+imh)-uanp(i))*sinlon(i)   &
                 + (vanp(i+imh)-vanp(i))*coslon(i)
         vn = vn + (uanp(i)-uanp(i+imh))*coslon(i)   &
                 + (vanp(i+imh)-vanp(i))*sinlon(i)
      enddo
      un = un*r2im
      vn = vn*r2im

! get V-wind at NP

      do i=1,imh
         v(i,    jm) = -un*cosl5(i) - vn*sinl5(i)
         v(i+imh,jm) = -v(i,jm)
      enddo

   endif

 end subroutine vpol5

!-----------------------------------------------------------------------
 end module sw_core
