module dyn_core

! !USES:
#ifndef USE_LIMA
   use fv_pack,      only: get_eta_level, pi, n_spong, n_diffu,  &
                           a_div, b_div, tra_step
#else
   use shr_kind_mod, only: r8 => shr_kind_r8
   use fv_pack,      only: get_eta_level, sc, se, dc, de, pi, ipft1, ipft2,  &
                           ifax, trigs, n_spong, n_diffu, too_big,           &
                           a_div, b_div
#endif

#ifdef MARS_GCM
   use fv_pack,      only: p_ref, fv_sponge_damp, fv_sponge_lev
#endif MARS_GCM


   use pft_module,   only: pft2d
   use sw_core,      only: c_sw, d_sw, upol5
   use timingModule, only: timing_on, timing_off

#ifdef SPMD
   use mod_comm,     only: mp_send4d_ns, mp_recv4d_ns,      &
                           mp_send3d_2,  mp_recv3d_2,       &
                           mp_send3d,    mp_recv3d,         &
                           mp_send2_ns,  mp_recv2_ns,  gid
#endif
   use fv_arrays_mod, only: fv_array_check, fv_array_sync, fv_array_limits, &
        fv_print_chksum, fv_print_chksums, ksp, kep
   use pmaxmin_mod, only: prt_maxmin_local

   implicit none

   private
   public cd_core


CONTAINS 

  subroutine cd_core(im,  jm,  km,   nq,  nx, jfirst, jlast,      &
       u,    v,  pt, delp,  pe,   pk,               &
       ns,  dt, ptop, umax, ae,  rcap, cp, akap,    &
       iord_mt, iord_tm, jord_mt, jord_tm,          &
       ng_d, ng_s, it, n_split, om, hs, sinp,       &
       cosp, cose, acosp, sinlon, coslon,           &
       cosl5, sinl5, f0, dx, rdx,  rdy, cx3, cy3,   &
       mfx, mfy, delpf, uc, vc, dpt, ptc,           &
#ifdef USE_LIMA
    wz3, pkc, wz, delpc, master, ig )
#else
    wz3, pkc, wz, delpc, master, ig, dyn_step )
#endif


! !INPUT PARAMETERS:

! Input paraterers:
    integer, intent(in):: im, jm, km
    integer, intent(in):: ig
    integer, intent(in):: nq
    integer, intent(in):: ns           ! # of time splits
    integer, intent(in):: nx           ! # of split pieces in longitude direction
    integer, intent(in):: jfirst       ! first latitude of the subdomain
    integer, intent(in):: jlast        ! last latitude of the subdomain
    integer, intent(in):: it
    integer, intent(in):: n_split

    real, intent(in):: ae          ! Radius of the Earth (m)
    real, intent(in):: om          ! rotation rate of the earth
    real, intent(in):: ptop        ! Model top pressure (pa)
    real, intent(in):: umax        ! Estimated upper bound of u-wind (m/s)
    real, intent(in):: dt          ! small time step in seconds
    real, intent(in):: rcap
    real, intent(in):: cp
    real, intent(in):: akap

!Grid Geometry
    real, intent(in):: rdy
    real, intent(in):: dx(jm)
    real, intent(in):: rdx(jm)

    integer, intent(in):: iord_mt, jord_mt
    integer, intent(in):: iord_tm, jord_tm
    integer, intent(in):: ng_d          ! Max NS dependencies
    integer, intent(in):: ng_s
#ifndef USE_LIMA
    integer, intent(in):: dyn_step
#endif
    logical, intent(in):: master   ! if master process (for printing)

! Input time independent arrays:
    real, intent(in):: hs(im,jfirst-ig:jlast+ig) ! surface geopotential
    real, intent(in)::  sinp(jm)
    real, intent(in)::  cosp(jm)
    real, intent(in):: acosp(jm)
    real, intent(in)::  cose(jm)

    real, intent(in):: sinlon(im)
    real, intent(in):: coslon(im)
    real, intent(in)::  sinl5(im)
    real, intent(in)::  cosl5(im)

#ifndef USE_LIMA
    real, intent(in):: f0(im,jfirst-ng_d:jlast+ng_d)
#else
    real, intent(in):: f0(jfirst-ng_d:jlast+ng_d)
#endif

! !INPUT/OUTPUT PARAMETERS:
    real, intent(inout)::  u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-Wind (m/s)
    real, intent(inout)::  v(im,jfirst-ng_d:jlast+ng_d,km)  ! v-wind (m/s)
    real, intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d,km)  ! Scaled-Potential temperature 
    real, intent(inout):: delp(im,jfirst:jlast,km)            ! Delta pressure (pa)
    real, intent(inout):: delpf(im,jfirst-ng_d:jlast+ng_d,km)    ! filtered delp

! Input/output: accumulated winds & mass fluxes on c-grid for large-
!               time-step transport
    real, intent(inout):: cx3(im,jfirst-ng_d:jlast+ng_d,km)!Accumulated Courant number in X
    real, intent(inout):: cy3(im,jfirst:jlast+1,km)        !Accumulated Courant number in Y
    real, intent(inout):: mfx(im,jfirst:jlast,km)          !Mass flux in X  (unghosted)
    real, intent(inout):: mfy(im,jfirst:jlast+1,km)        !Mass flux in Y

! !OUTPUT PARAMETERS:
    real, intent(out):: pe(im,km+1,jfirst:jlast)         ! Edge pressure (pascal)
    real, intent(out):: pk(im,jfirst:jlast,km+1)         ! pe ** kappa

! work arrays (useless output)
    real, intent(out)::    uc(im,jfirst-ng_d:jlast+ng_d,km)
    real, intent(out)::    vc(im,jfirst-2:   jlast+2,   km)
    real, intent(out):: delpc(im,jfirst:jlast,km)

    real, intent(out):: ptc(im,jfirst:jlast,km)
    real, intent(out):: dpt(im,jfirst-1:jlast+1,km)

    real, intent(out):: wz3(im,jfirst-1:jlast  ,km+1)

    real, intent(out):: pkc(im,jfirst-1:jlast+1,km+1) 
    real, intent(out)::  wz(im,jfirst-1:jlast+1,km+1)

! ! !DESCRIPTION:
!    Perform a dynamical update for one small time step; the small
!    time step is limitted by the fastest wave within the Lagrangian control-
!    volume 
!
! !REVISION HISTORY:
!     SJL  99.01.01:   Original SMP version
!     WS   99.04.13:   MPI-1; added jfirst:jlast concept
!     SJL  99.07.15:   Merged c_core and d_core to this routine
!     SJL  99.12.23:   More comments; general optimization; reduction
!                      of redundant computation & communication
!     SJL  03.11.15:   Restructure c_sw to reduce communication (added ig)
!     SJL  04.01.30:   Fix bug(s) at poles; revert to old sw_core codes
!
!EOP
!---------------------------------------------------------------------
!BOC
! Local 2D arrays:
    integer, parameter:: pft_buff = 1  ! increase buffer size for pft2d
    real   wk(im,jfirst:jlast+2+pft_buff)
    real  wk2(im,jfirst:jlast+1+pft_buff)
    real  wk1(im,jfirst-1:jlast+1+pft_buff)
    real  wk3(im,jfirst-1:jlast+1)
    real  p1d(im)
    real   phalf(km+1)

!Local scalars
    real  dt0, dt5
    real  dl, dp
    real  dt_divg

    integer i, j, k
    integer js1g1, js2g0
    integer jn2g0, jn1g1, jsfc, jnfc
    integer iord , jord, sord

    real  dtdy, dydt, dtdy5, tdy5
    real  zt_c   
    real  zt_d  
    real  tau, fac, pk4
    real, parameter:: esl = 1.e-20

! Declare permanent local arrays
    real, allocatable, save :: pfull(:)
#ifndef USE_LIMA
    real, allocatable, save :: fc(:,:)
#else
    real, allocatable, save :: fc(:)
#endif
    real, allocatable, save :: cdx(:,:), cdy(:,:)
    real, allocatable, save :: dtdx(:), dtdxe(:), txe5(:), dtxe5(:)
    real, allocatable, save :: dyce(:), cy(:)
    real, allocatable, save :: dtdx2(:), dtdx4(:),  dxdt(:), dxe(:)
    real, allocatable, save :: cye(:),    dycp(:),  rdxe(:)
    integer ng_c, js2gc, jn1gc
    integer :: kepp !use to do km+1: Balaji

    data dt0 / 0./
    save dtdy, dydt, dtdy5, tdy5
    save zt_c, zt_d

    call fv_print_chksums( 'Entering  cd_core' )
    call fv_array_check( LOC(u) )
    call fv_array_check( LOC(v) )
    call fv_array_check( LOC(pt) )
    call fv_array_check( LOC(delp) )
    call fv_array_check( LOC(delpf) )
    call fv_array_check( LOC(cx3) )
    call fv_array_check( LOC(cy3) )
    call fv_array_check( LOC(mfx) )
    call fv_array_check( LOC(mfy) )
    call fv_array_check( LOC(pe) )
    call fv_array_check( LOC(pk) )
    call fv_array_check( LOC(uc) )
    call fv_array_check( LOC(vc) )
    call fv_array_check( LOC(delpc) )
    call fv_array_check( LOC(ptc) )
    call fv_array_check( LOC(dpt) )
    call fv_array_check( LOC(wz3) )
    call fv_array_check( LOC(pkc) )
    call fv_array_check( LOC(wz) )
#ifdef SPMD
    ng_c = 2
#else
    ng_c = 0
#endif
    js2gc  = max(2,jfirst-ng_c)
    jn1gc  = min(jm,jlast+ng_c)

! Set loop limits
    js1g1 = max(1,jfirst-1)
    js2g0 = max(2,jfirst)
    jsfc  = max(2,jfirst-2) 

    jn2g0 = min(jm-1,jlast)
    jn1g1 = min(jm,  jlast+1)
    jnfc  = min(jm,  jlast+2)
    kepp = kep; if( kep.EQ.km )kepp = kep + 1

    if( abs(dt0-dt) > 0.1 ) then
        if ( .not. allocated( dtdx ) ) then
            allocate( dtdx(jm), dtdx2(jm), dtdx4(jm), dtdxe(jm), dxdt(jm),    &
                 dxe(jm),   cye(jm),  dycp(jm),  rdxe(jm), txe5(jm),    &
                 dtxe5(jm),  dyce(jm),    cy(jm) )

            allocate( cdx(js2g0:jn1g1,km) )
            allocate( cdy(js2g0:jn1g1,km) )
#ifndef USE_LIMA
            allocate( fc(im,jsfc:jnfc) )
#else
            allocate( fc(jsfc:jnfc) )
#endif 


#ifndef USE_LIMA

! Compute coriolis parameter at cell corners.
!         do j=2,jm
            do j=jsfc, jnfc
               do i=1,im
                  p1d(i) = 0.25*(f0(i,j) + f0(i,j-1))
               enddo
               fc(1,j) = p1d(im) + p1d(1)
               do i=2,im
                  fc(i,j) = p1d(i-1) + p1d(i)
               enddo
            enddo
#else               

! Compute coriolis parameter at cell corners.
!         do j=2,jm
            do j=jsfc, jnfc
               fc(j) = 0.5*(f0(j) + f0(j-1))
            enddo

#endif

        endif

        dt0 = dt
        dt5 = 0.5*dt

        dl = (pi+pi)/im
        dp = pi/(jm-1)

        dtdy  = dt *rdy
        dtdy5 = dt5*rdy
        dydt  = (ae*dp) / dt
        tdy5  = 0.5/dtdy

        do j=2,jm-1
           dtdx(j)  = dt / dx(j)
           dxdt(j)  = dx(j) / dt
           dtdx2(j) = 0.5*dtdx(j)
           dtdx4(j) = 0.5*dtdx2(j)
           dycp(j)  = ae*dp/cosp(j)
           cy(j)    =  rdy * acosp(j)
        enddo

        do j=2,jm
           dxe(j)   = ae*dl*cose(j)
           rdxe(j)  = 1. / dxe(j)
           dtdxe(j) = dt / dxe(j)
           dtxe5(j) = 0.5*dtdxe(j)
           txe5(j)  = 0.5/dtdxe(j)
           cye(j)  =  1. / (ae*cose(j)*dp)
           dyce(j)  = ae*dp/cose(j)
        enddo

! Checking if Lagrangian (integer part) flux is needed or not
        zt_c = abs(umax*dt5) / (dl*ae)             ! for C-grid
        zt_d = abs(umax*dt)  / (dl*ae)             ! for D-grid

!--------------------------------------------------------
! Divergence damping coeff: Cd = dx*dy / dt_divg [m**2/s]
!--------------------------------------------------------
        if(master) write(6,*) 'Divergence damping coefficient (m**2/sec x E6):'
        dt_divg = 450.*180./max(im, 2*(jm-1))  ! small-time-step as determined
! by the default setup (n_split=0)
        allocate ( pfull(km) )
#ifdef MARS_GCM
        call get_eta_level(km, p_ref, pfull, phalf, 0.01) 
#else
        call get_eta_level(km, 1.E5, pfull, phalf, 0.01) 
#endif MARS_GCM 


        do k=1,km
           if ( km <= 2 ) then 
! Shallow-water mode
               tau = 1.
           else 
#ifdef OLD_DIVD
               if ( n_spong == 0 ) then
                   tau = 1.0              ! weak damping
               else
                   if(k == 1) then
                       tau = 4.0
                   elseif(k == 2) then
                       tau = 2.8
                   else
                       tau = 1.8 + tanh(0.5*log(50./pfull(k))) 
                   endif
               endif
#elif MARS_GCM
!           Set break-point at fv_sponge_lev  ( 0.5 mb is the default) 
             tau =  1.0 +  0.5*fv_sponge_damp*(1.0 + tanh(0.9*log(fv_sponge_lev/pfull(k)) ) )
#else
               if ( n_spong == 0 ) then
                   if(k == 1) then
                       tau = 2.0
                   else
                       tau = 1.
                   endif
               else
                   if(k == 1) then
                       tau = 3.0
                   elseif( k == 2 .and. pfull(k) < 10. ) then
                       tau = 1.25
                   else
                       tau = 1.
                   endif
               endif
#endif
           endif
           tau = max(tau, 1.) / (128.*dt_divg)
           do j=js2g0,jn1g1
              fac = tau * ae * ( a_div/cose(j) + b_div )
              cdx(j,k) = fac*dp
              cdy(j,k) = fac*dl
           enddo
           if(master) write(6,*) k, pfull(k), cdx(js2g0,k)*ae*dl*cose(js2g0)*1.E-6
        enddo
    endif ! end init

#ifdef SPMD
    if( it /= 1 ) then             ! not the first call
        call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_s, u)
        call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, v)
    endif
    call fv_array_sync()
#endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the C-grid
!-----------------------------------------------------------------

!$omp parallel do default(shared) private(k, iord, jord)
    do k=ksp,kep !1,km
       iord = 1
       jord = 1
#ifndef USE_LIMA

       call c_sw(  u(1,jfirst-ng_d,k),   v(1,jfirst-ng_d,k),        &
            pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),          &
            uc(1,jfirst-ng_d,k),   vc(1,jfirst-2,k),          &
            ptc(1,jfirst,k),        delpf(1,jfirst-ng_d,k),    &
            delpc(1,jfirst,k),                                 &
            cosp,   acosp,   cose,   coslon,   sinlon,         &
            dxdt,   dxe,     dtdx2,  dtdx4,    dtxe5,  rdxe,   &
            dycp,   dydt,    dtdy5,  cye,      fc,             &
            zt_c,   esl,     rcap,   im,       jm,             &
            jfirst, jlast,   ng_c,    ng_d,    ng_s,           &
            js2g0,  jn2g0,   js2gc,   jn1gc,                   &
            iord,   jord,    cosl5,   sinl5, 0 )
#else               
       call c_sw(  u(1,jfirst-ng_d,k),   v(1,jfirst-ng_d,k),        &
            pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),          &
            uc(1,jfirst-ng_d,k),   vc(1,jfirst-2,k),          &
            ptc(1,jfirst,k),        delpf(1,jfirst-ng_d,k),    &
            delpc(1,jfirst,k),                                 &
            cosp,   acosp,   cose,   coslon,   sinlon,         &
            dxdt,   dxe,     dtdx2,  dtdx4,    dtxe5,  rdxe,   &
            dycp,   dydt,    dtdy5,  cye,      fc,             &
            ifax,   trigs,   dc(1,js2g0),      sc,             &
            zt_c,   esl,     rcap,   im,       jm,             &
            jfirst, jlast,   ng_c,    ng_d,    ng_s,           &
            js2g0,  jn2g0,   js2gc,   jn1gc,                   &
            iord,   jord,    cosl5,   sinl5 )
#endif
    enddo
    call fv_array_sync

    call geopk(ptop, pe, delpc, pkc, wz, hs(1,jfirst), ptc, im, jm, km,   &
         jfirst, jlast, 0, akap, nx, 0, .false.)

#ifdef SPMD
    call mp_send3d_2(gid+1, gid-1, im, jm, km+1, &
         1, im, jfirst-1, jlast+1, 1, km+1, &
         1, im, jlast, jlast, 1, km+1, pkc, wz)
#endif


!$omp parallel do private(i, j, k, p1d, wk, wk2)
    do k=ksp,kep !1,km
       do j=js2g0,jn2g0
          do i=1,im
             p1d(i) = pkc(i,j,k+1) - pkc(i,j,k)
          enddo

          uc(1,j,k) = uc(1,j,k) + dtdx2(j) * (                      &
               (wz(im,j,k+1)-wz(1,j,k))*(pkc(1,j,k+1)-pkc(im,j,k))   &
               + (wz(im,j,k)-wz(1,j,k+1))*(pkc(im,j,k+1)-pkc(1,j,k)))  &
               / (p1d(1)+p1d(im))
          do i=2,im
             uc(i,j,k) = uc(i,j,k) + dtdx2(j) * (                     &
                  (wz(i-1,j,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i-1,j,k))    &
                  + (wz(i-1,j,k)-wz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k)) )  &
                  / (p1d(i)+p1d(i-1))
          enddo
       enddo

#ifndef USE_LIMA
       call pft2d(uc(1,js2g0,k), im, jn2g0-js2g0+1,  wk, wk2, 1)
#else               
       call pft2d(uc(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,        &
            jn2g0-js2g0+1,  ifax,   trigs, wk, wk2)
#endif


       call upol5(uc(1,jfirst,k), vc(1,jfirst,k), im, jm,    &
            coslon,  sinlon,  &
            cosl5, sinl5, jfirst, jlast, ng_d)
    enddo                 ! end k paralle loop
    call fv_array_sync

#ifdef SPMD
    call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, uc)

    call mp_recv3d_2(gid-1, im, jm, km+1, &
         1, im, jfirst-1, jlast+1, 1, km+1, &
         1, im, jfirst-1, jfirst-1, 1, km+1, pkc, wz)
#endif
    call fv_array_sync

!$omp parallel do private(i, j, k, wk, wk1)
    do k=ksp,kep !1,km
       do j=js1g1,jlast
          do i=1,im
             wk1(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
          enddo
       enddo

       do j=js2g0,jlast
          do i=1,im
             vc(i,j,k) = vc(i,j,k) + dtdy5/(wk1(i,j)+wk1(i,j-1)) *   &
                  ((wz(i,j-1,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j-1,k))    &
                  + (wz(i,j-1,k)-wz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)))
          enddo
       enddo

#ifndef USE_LIMA
       call pft2d(vc(1,js2g0,k), im, jlast-js2g0+1, wk, wk1, 2)
#else               
       call pft2d(vc(1,js2g0,k), se(js2g0), de(1,js2g0), im,      &
            jlast-js2g0+1,  ifax, trigs, wk, wk1 )
#endif

    enddo
    call fv_array_sync

#ifdef SPMD
!------------
! Receive uc
!------------
    call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, uc)
!------------
! Send/recv vc
!------------
    call mp_send3d(gid-1, gid+1, im, jm, km,               &
         1, im, jfirst-2, jlast+2, 1, km,        &
         1, im, jfirst, jfirst, 1, km, vc)
    call mp_recv3d(gid+1, im, jm, km,                       &
         1, im, jfirst-2, jlast+2, 1, km, &
         1, im, jlast+1, jlast+1, 1, km, vc)
#endif
!the vc array shouldn't be used until receive complete?
    call fv_array_sync
#ifdef DEBUG
    call prt_maxmin_local(gid, 'uc', uc, im, jm, km, jfirst-ng_d, jlast+ng_d)
    call prt_maxmin_local(gid, 'vc', vc, im, jm, km, jfirst-2, jlast+2)
#endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the D-grid
!-----------------------------------------------------------------

!    call timing_on(' D_SW')
!$omp parallel do default(shared) private(k, iord, jord, sord)
    do k=ksp,kep !1,km
       sord = min(iord_tm, jord_tm)
       if( k <= n_spong ) then
! Apply first order scheme for damping the sponge layer
           iord = 1
           jord = 1
           sord = 1                          ! scheme for pt & delp
       elseif( k == n_spong+1 .and. pfull(k) < 10.0  ) then
           iord = min(2, iord_mt)
           jord = min(2, jord_mt)
           sord = min(2, iord_tm, jord_mt)
       elseif( k >= km-n_diffu+1 .and. pfull(k) > 800.0  ) then
           iord = 2
           jord = 2
           sord = 2
       else
! Apply the chosen high-order scheme
           iord = iord_mt
           jord = jord_mt
       endif

#ifndef USE_LIMA
       call d_sw( u(1,jfirst-ng_d,k),      v(1,jfirst-ng_d,k),     &
            uc(1,jfirst-ng_d,k),    vc(1,jfirst-1,k),        &
            pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),         &
            delpf(1,jfirst-ng_d,k), cx3(1,jfirst-ng_d,k),    &
            cy3(1,jfirst,k),        mfx(1,jfirst,k),         &
            mfy(1,jfirst,k), cdx(js2g0,k),  cdy(js2g0,k),    &
            dtdx,   dtdxe,  dtxe5,  txe5,  dyce,  rdx,  cy,  &
            dx,  f0(1,jfirst-ng_d), js2g0,  jn1g1, im,  jm,    &
            jfirst, jlast,  ng_d,  ng_s,   nq,    iord,      &
            jord,   sord,   zt_d,   rcap,  esl,   dtdy,      &
            dtdy5,  tdy5,   rdy,    cosp,  acosp, cose,      &
            coslon, sinlon, cosl5, sinl5, dyn_step )
#else               
       call d_sw( u(1,jfirst-ng_d,k),      v(1,jfirst-ng_d,k),     &
            uc(1,jfirst-ng_d,k),    vc(1,jfirst-1,k),        &
            pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),         &
            delpf(1,jfirst-ng_d,k), cx3(1,jfirst-ng_d,k),    &
            cy3(1,jfirst,k),        mfx(1,jfirst,k),         &
            mfy(1,jfirst,k), cdx(js2g0,k),  cdy(js2g0,k),    &
            dtdx,   dtdxe,  dtxe5,  txe5,  dyce,  rdx,  cy,  &
            dx,  f0(jfirst-ng_d), js2g0,  jn1g1, im,  jm,    &
            jfirst, jlast,  ng_d,  ng_s,   nq,    iord,      &
            jord,   sord,   zt_d,   rcap,  esl,   dtdy,      &
            dtdy5,  tdy5,   rdy,    cosp,  acosp, cose,      &
            coslon, sinlon, cosl5, sinl5 )
#endif

    enddo
    call fv_array_sync

!    call timing_off(' D_SW')

!      call timing_on(' geop_d')
    call geop_d(ptop, pe, delp, pkc, wz, hs(1,jfirst), pt, im, jm,  &
         km,   jfirst, jlast, ng_d, cp, akap, nx, it==n_split, .true.)
!      call timing_off(' geop_d')

#ifdef DEBUG
    call prt_maxmin_local(gid, 'pt_2', pt, im, jm, km, jfirst-ng_d, jlast+ng_d)
#endif

#ifdef SPMD
!-------------------------------------
!  Send pkc and wz NS boundary regions
!-------------------------------------
    call mp_send2_ns(im, jm, km+1, jfirst, jlast, 1, km+1, 1, pkc, wz)
#endif

    if ( it /= n_split ) then          !  not the last call
!$omp parallel do private(i, j, k, wk, wk2)
        do k=ksp,kep !1,km
           do j=jfirst,jlast
              do i=1,im
                 delpf(i,j,k) = delp(i,j,k)
              enddo
           enddo

#ifndef USE_LIMA
           call pft2d(delpf(1,js2g0,k), im, jn2g0-js2g0+1, wk, wk2, 1)
#else               
           call pft2d(delpf(1,js2g0,k), sc(js2g0), dc(1,js2g0),   &
                im, jn2g0-js2g0+1, ifax, trigs, wk, wk2)
#endif


        enddo
    else
! Last call
!$omp parallel do private(i, j, k)
        do k=ksp,kepp !1,km+1 Balaji
           do j=jfirst,jlast
              do i=1,im
                 pk(i,j,k) = pkc(i,j,k)
              enddo
           enddo
        enddo
    endif
    call fv_array_sync

#ifdef SPMD
    call mp_recv2_ns(im, jm, km+1, jfirst, jlast, 1, km+1, 1, pkc, wz)

    if ( it == n_split ) then          !  last call; send data for tracer transport
#ifndef USE_LIMA
        if ( tra_step /= 1 )       &
             call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, cx3)
#else
        call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, cx3)
#endif
        call mp_send3d_2(gid-1, gid+1, im, jm, km,                            &
             1, im, jfirst, jlast+1, 1, km,                   &
             1, im, jfirst, jfirst,  1, km, cy3, mfy)
    else 

#if !defined (SW_DYN)
        call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, pt)
#endif
        call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, delpf)

    endif
    call fv_array_sync
#endif

!$omp parallel do private(i, j, k)
    do k=ksp,kep !1,km
       do j=js1g1,jn1g1                  ! dpt needed NS
          do i=1,im                      ! wz, pkc ghosted NS
             dpt(i,j,k)=(wz(i,j,k+1)+wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j,k))
          enddo
       enddo
    enddo
    call fv_array_sync

    pk4 = 4.*ptop**akap

!$omp parallel do private(i, j, k, wk1, wk3)

    do 4500 k=ksp,kepp !1,km+1 Balaji
       if ( k == 1 ) then
           do j=js2g0,jlast
              do i=1,im
                 wz3(i,j,1) = 0.
                 wz(i,j,1) = 0.
              enddo
           enddo
           do j=js2g0,jn1g1
              do i=1,im
                 pkc(i,j,1) = pk4
              enddo
           enddo

       else

           do j=max(2,jfirst-1), jn2g0
              wk3(1,j) = (wz(1,j,k)+wz(im,j,k)) *        &
                   (pkc(1,j,k)-pkc(im,j,k))
              do i=2,im
                 wk3(i,j) = (wz(i,j,k)+wz(i-1,j,k)) *       &
                      (pkc(i,j,k)-pkc(i-1,j,k))
              enddo
           enddo

           do j=max(2,jfirst-1), jn2g0
              do i=1,im-1
                 wk1(i ,j) = wk3(i, j) + wk3(i+1,j)
              enddo
              wk1(im,j) = wk3(im,j) + wk3(1,j)
           enddo

           if ( jfirst == 1 ) then
               do i=1,im
                  wk1(i,1) = 0.
               enddo
           endif

           if ( jlast == jm ) then
               do i=1,im
                  wk1(i,jm) = 0.
               enddo
           endif

           do j=js2g0,jlast                          ! wk1 ghosted S
              do i=1,im
                 wz3(i,j,k) = wk1(i,j) + wk1(i,j-1)
              enddo
           enddo

! N-S walls
           do j=js2g0,jn1g1                        ! wk1 needed N
              do i=1,im                            ! wz, pkc ghosted NS
                 wk1(i,j) = (wz(i,j,k)+wz(i,j-1,k))*(pkc(i,j,k)-pkc(i,j-1,k))
              enddo
           enddo

           do j=js2g0,jn1g1                         ! wk3 needed N
              wk3(1,j) = wk1(1,j) + wk1(im,j)    ! wk1 ghosted N
              do i=2,im
                 wk3(i,j) = wk1(i,j) + wk1(i-1,j)   ! wk1 ghosted N
              enddo
           enddo

           do j=js2g0,jn2g0
              do i=1,im
                 wz(i,j,k) = wk3(i,j) + wk3(i,j+1)  ! wk3 ghosted N
              enddo
           enddo

           do j=js1g1,jn1g1
              wk1(1,j) = pkc(1,j,k) + pkc(im,j,k)
              do i=2,im
                 wk1(i,j) = pkc(i,j,k) + pkc(i-1,j,k)
              enddo
           enddo

           do j=js2g0,jn1g1
              do i=1,im
                 pkc(i,j,k) = wk1(i,j) + wk1(i,j-1)
              enddo
           enddo
       endif
4500 end do !continue
    call fv_array_sync
!$omp parallel do private(i, j, k, wk, wk1, wk2, wk3)
    do 6000 k=ksp,kep !1,km

       do j=js1g1,jn1g1
          wk1(1,j) = dpt(1,j,k) + dpt(im,j,k)
          do i=2,im
             wk1(i,j) = dpt(i,j,k) + dpt(i-1,j,k)
          enddo
       enddo

       do j=js2g0,jn1g1
          do i=1,im
             wk2(i,j) = wk1(i,j)     + wk1(i,j-1)
             wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
          enddo
       enddo

       do j=js2g0,jlast
          do i=1,im-1
             wk3(i,j) = uc(i,j,k) + dtdxe(j)/(wk(i,j) + wk(i+1,j))  &
                  * (wk2(i,j)-wk2(i+1,j)+wz3(i,j,k+1)-wz3(i,j,k))
          enddo
          wk3(im,j) = uc(im,j,k) + dtdxe(j)/(wk(im,j) + wk(1,j))   &
               * (wk2(im,j)-wk2(1,j)+wz3(im,j,k+1)-wz3(im,j,k))
       enddo

       do j=js2g0,jn2g0
          do i=1,im
             wk1(i,j) = vc(i,j,k) + dtdy/(wk(i,j)+wk(i,j+1)) *     &
                  (wk2(i,j)-wk2(i,j+1)+wz(i,j,k+1)-wz(i,j,k))
          enddo
       enddo

#if (!defined ALT_PFT)

#ifndef USE_LIMA
       call pft2d(wk3(1,js2g0), im, jlast-js2g0+1, wk, wk2, 2)
       call pft2d(wk1(1,js2g0), im, jn2g0-js2g0+1, wk, wk2, 1)
#else               
       call pft2d(wk3(1,js2g0), se(js2g0), de(1,js2g0), im,     &
            jlast-js2g0+1,  ifax, trigs, wk, wk2 )
       call pft2d(wk1(1,js2g0), sc(js2g0), dc(1,js2g0), im,     &
            jn2g0-js2g0+1,  ifax, trigs, wk, wk2 )

#endif

#endif

       do j=js2g0,jn2g0
          do i=1,im
             v(i,j,k) = v(i,j,k) + wk1(i,j)
             u(i,j,k) = u(i,j,k) + wk3(i,j)
          enddo
       enddo

       if ( jlast == jm ) then
           do i=1,im
              u(i,jlast,k) = u(i,jlast,k) + wk3(i,jlast)
           enddo
       endif

#ifdef ALT_PFT

#ifndef USE_LIMA
       call pft2d(u(1,js2g0,k), im, jlast-js2g0+1, wk, wk2, 2)
       call pft2d(v(1,js2g0,k), im, jn2g0-js2g0+1, wk, wk2, 1)
#else               
       call pft2d(u(1,js2g0,k), se(js2g0), de(1,js2g0), im,       &
            jlast-js2g0+1,  ifax, trigs, wk, wk2)
       call pft2d(v(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,       &
            jn2g0-js2g0+1,  ifax, trigs, wk, wk2)
#endif

#endif
6000 end do !continue
    call fv_array_sync
#ifdef SPMD
    if ( it /= n_split ) then            ! not the last call
#if !defined (SW_DYN)
        call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, pt)
#endif
        call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, delpf)
        call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_s, u)
        call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, v)

    endif
#endif
!EOC
    call fv_array_sync
    call fv_print_chksums( 'Exiting  cd_core' )
  end subroutine cd_core
!

 subroutine geop_d(ptop, pe, delp, pk, wz, hs, pt, im, jm, km,  &
                   jfirst, jlast, nd, cp, akap, nx, last_call, dp_check)

! !INPUT PARAMETERS:

   integer, intent(in):: im, jm, km, jfirst, jlast
   integer, intent(in):: nx                      ! # of pieces in longitude direction
   integer, intent(in):: nd
   real, intent(in):: akap, cp, ptop
   real, intent(in):: hs(im,jfirst:jlast)
   logical, intent(in):: dp_check
   logical, intent(in):: last_call

! !INPUT/OUTPUT PARAMETERS:
! pt and delp could be adjusted if dp_check is true
   real, intent(inout)::  pt(im,jfirst-nd:jlast+nd,km) 
   real, intent(inout):: delp(im,jfirst:jlast,km)     

! !OUTPUT PARAMETERS
   real, intent(out):: wz(im,jfirst-1:jlast+1,km+1)  ! space N*1 S*1
   real, intent(out):: pk(im,jfirst-1:jlast+1,km+1)  ! space N*1 S*1
   real, intent(out):: pe(im,km+1,jfirst:jlast)      ! only if last_call

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!
! !REVISION HISTORY:
!
!  WS  99.10.22: MPIed SJ's original SMP version
!  SJL 00.01.01: Merged C-core and D-core computation
!                SMP "decmposition" in E-W by combining i and j loops
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      integer i, j, k
      integer ixj, jp, it, i1, i2

      real p1d(im)
      real ptk
      real dpmin
      real dp
#ifndef MARS_GCM
      parameter (dpmin = 0.02)              ! unit: pascal
#endif 
      integer :: ixjs, ixje


#ifdef MARS_GCM
!      this valuation has been adopted in the cubed-sphere code
      dpmin = 0.01*ptop
#endif

      call fv_print_chksums( 'Entering  geop_d' )
      call fv_array_check( LOC(hs) )
      call fv_array_check( LOC(pt) )
      call fv_array_check( LOC(delp) )
      call fv_array_check( LOC(wz) )
      call fv_array_check( LOC(pk) )
      call fv_array_check( LOC(pe) )

      it = im / nx
      jp = nx * ( jlast - jfirst + 1 )
      call fv_array_limits( 1, jp, ixjs, ixje )

!$omp parallel do default(shared) private(i1, i2, ixj, i, j, k, p1d, ptk, dp)
!     do 2000 j=jfirst,jlast
      do 2000 ixj=ixjs, ixje !1, jp
!        
         j  = jfirst + (ixj-1)/nx
         i1 = 1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1
        
        ptk  = ptop ** akap
        do i=i1,i2
           p1d(i) = ptop
           pk(i,j,1) = ptk
           wz(i,j,km+1) = hs(i,j)
        enddo

        if( last_call ) then
            do i=i1,i2
               pe(i,1,j) = ptop
            enddo
        endif

#if !defined(SW_DYN)
        if( dp_check ) then
            do k=1, km-1
               do i=i1,i2
                  if(delp(i,j,k) < dpmin) then
! Remap from below and mix pt
                      dp = dpmin - delp(i,j,k)
                      pt(i,j,k) = (pt(i,j,k)*delp(i,j,k) + pt(i,j,k+1)*dp) / dpmin
                      delp(i,j,k) = dpmin
                      delp(i,j,k+1) = delp(i,j,k+1) - dp
                  endif
               enddo
            enddo

! Bottom (k=km):
            do i=i1,i2
               if(delp(i,j,km) < dpmin) then
! Remap from above and mix pt
                   dp = dpmin - delp(i,j,km)
                   pt(i,j,km) = (pt(i,j,km)*delp(i,j,km) + pt(i,j,km-1)*dp)/dpmin
                   delp(i,j,km) = dpmin
                   delp(i,j,km-1) = delp(i,j,km-1) - dp
               endif
            enddo
        endif
#endif

! Top down
        do k=2,km+1
           do i=i1,i2
              p1d(i)  = p1d(i) + delp(i,j,k-1)
              pk(i,j,k) = p1d(i) ** akap
           enddo
           if( last_call ) then
               do i=i1,i2
                  pe(i,k,j) = p1d(i)
               enddo
           endif
        enddo

! Bottom up
        do k=km,1,-1
           do i=i1,i2
              wz(i,j,k) = wz(i,j,k+1) + pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
           enddo
        enddo
2000    continue
        call fv_array_sync()
        call fv_print_chksums( 'Exiting  geop_d' )
 end subroutine geop_d


 subroutine geopk(ptop, pe, delp, pk, wz, hs, pt, im, jm, km,  &
      jfirst, jlast, nd, akap, nx, id, dp_check)

   integer im, jm, km, jfirst, jlast, id
   integer nx                        ! # of pieces in longitude direction
   integer nd
   real    akap, ptop
   logical dp_check
   real hs(im,jfirst:jlast)

! !INPUT/OUTPUT PARAMETERS:
   real  pt(im,jfirst-nd:jlast+nd,km)  ! only altered if dp_check
   real delp(im,jfirst:jlast,km)      ! only altered if dp_check

! !OUTPUT PARAMETERS
   real wz(im,jfirst-1:jlast+1,km+1)  ! space N*1 S*1
   real pk(im,jfirst-1:jlast+1,km+1)  ! space N*1 S*1
   real pe(im,km+1,jfirst:jlast)      ! only if id == 1

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!
! !REVISION HISTORY:
!
!  WS  99.10.22: MPIed SJ's original SMP version
!  SJL 00.01.01: Merged C-core and D-core computation
!                SMP "decmposition" in E-W by combining i and j loops
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
   integer i, j, k
   real    p1d(im)
   integer ixj, jp, it, i1, i2

   integer :: ixjs, ixje

   real ptk
   real dpmin
   real dp
#ifdef MARS_GCM
!      this valuation has been adopted in the cubed-sphere code
     dpmin = 0.01*ptop
#else
!    rjw    Note that dpmin was 0.02 in the memphis code  
   parameter (dpmin = 0.1)              ! unit: pascal
#endif MARS_GCM

   call fv_print_chksums( 'Entering  geopk' )
   call fv_array_check( LOC(hs) )
   call fv_array_check( LOC(pt) )
   call fv_array_check( LOC(delp) )
   call fv_array_check( LOC(wz) )
   call fv_array_check( LOC(pk) )
   call fv_array_check( LOC(pe) )

   it = im / nx
   jp = nx * ( jlast - jfirst + 1 )
   call fv_array_limits( 1, jp, ixjs, ixje )

!$omp parallel do default(shared) private(i1, i2, ixj, i, j, k, p1d, ptk, dp)

!     do 2000 j=jfirst,jlast
   do 2000 ixj=ixjs, ixje !1, jp

      j  = jfirst + (ixj-1)/nx
      i1 = 1 + it * mod(ixj-1, nx)
      i2 = i1 + it - 1

      ptk  = ptop ** akap
      do i=i1,i2
         p1d(i) = ptop
         pk(i,j,1) = ptk
         wz(i,j,km+1) = hs(i,j)
      enddo

      if(id .eq. 1) then
          do i=i1,i2
             pe(i,1,j) = ptop
          enddo
      endif

      if( dp_check ) then

          do k=1, km-1
             do i=i1,i2
                if(delp(i,j,k) < dpmin) then
! Remap from below and mix pt
                    dp = dpmin - delp(i,j,k)
                    pt(i,j,k) = (pt(i,j,k)*delp(i,j,k) + pt(i,j,k+1)*dp) / dpmin
                    delp(i,j,k) = dpmin
                    delp(i,j,k+1) = delp(i,j,k+1) - dp
                endif
             enddo
          enddo

! Bottom (k=km):
          do i=i1,i2
             if(delp(i,j,km) < dpmin) then
! Remap from above and mix pt
                 dp = dpmin - delp(i,j,km)
                 pt(i,j,km) = (pt(i,j,km)*delp(i,j,km) + pt(i,j,km-1)*dp)/dpmin
                 delp(i,j,km) = dpmin
                 delp(i,j,km-1) = delp(i,j,km-1) - dp
             endif
          enddo
      endif

! Top down
      do k=2,km+1
         do i=i1,i2
            p1d(i)  = p1d(i) + delp(i,j,k-1)
            pk(i,j,k) = p1d(i) ** akap
         enddo
         if(id == 1) then
             do i=i1,i2
                pe(i,k,j) = p1d(i)
             enddo
         endif
      enddo

! Bottom up
      do k=km,1,-1
         do i=i1,i2
            wz(i,j,k) = wz(i,j,k+1)+pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
         enddo
      enddo
2000 end do !continue
   call fv_array_sync()
   call fv_print_chksums( 'Exiting  geopk' )
 end subroutine geopk

end module dyn_core
