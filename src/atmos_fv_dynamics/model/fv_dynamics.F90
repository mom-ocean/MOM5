module fv_dynamics_mod
  implicit none
  private
  public :: fv_dynamics
contains
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: driver for the finite-volume dynamical core,  
!           MPI in latitudinal direction SMP in (mostly) vertical direction
!
! !INTERFACE:
#ifndef USE_LIMA
  subroutine fv_dynamics(im,   jm,    km,    jfirst, jlast,     &
       nq,   pnats, p_map, consv,             &
       u,    v,     delp,  pt,     q,         &
       ps,   pe,    pk,    pkz,    phis,      &
       omga, peln,  ptop,  om,     ndt,       &
       r_vir, cp,    rg,    cappa,  ae, ua, va, Time_next )

! !USES:
    use fv_pack, only: acap, cosp, sinp, cose, sine, coslon, sinlon,    &
         rcap, acosp, ng_d, ng_s, cosl5, sinl5, f_d,      &
         one_by_dy,  dx,  one_by_dx, u_ghosted, n_split,  &
         iord_mt,  iord_tm, iord_tr,                      &
         jord_mt,  jord_tm, jord_tr,                      &
         kord_mt,  kord_tm, kord_tr,                      &
         tra_step, dyn_step,  fill, master, n_zonal,  &
         use_tendency, &
         ighost, ks, ak, bk, cp_vir, icase

    use mapz_module,    only: te_map, geo_map, benergy, p_energy
    use dyn_core,       only: cd_core

    use pft_module,     only: pft2d
    use timingModule,   only: timing_on, timing_off

#ifdef SPMD
    use mod_comm,       only: gid, mp_send4d_ns, mp_recv4d_ns
#endif
    use time_manager_mod, only: time_type

#else

    subroutine fv_dynamics(im,   jm,    km,    jfirst, jlast,     &
         nq,   pnats, p_map, consv,             &
         u,    v,     delp,  pt,     q,         &
         ps,   pe,    pk,    pkz,    phis,      &
         omga, peln,  ptop,  om,     ndt,       &
         r_vir, cp,    rg,    cappa,  ae, ua, va )

! !USES:
      use fv_pack,        only: acap, cosp, sinp, cose, sine, coslon, sinlon,    &
           rcap, acosp, ng_d, ng_s, cosl5, sinl5, f_d,      &
           one_by_dy,  dx,  one_by_dx, u_ghosted, n_split,  &
           iord_mt,  iord_tm, iord_tr,                      &
           jord_mt,  jord_tm, jord_tr,                      &
           kord_mt,  kord_tm, kord_tr,                      &
           fill, master, ifax, trigs, sc, dc, n_zonal,      &
           use_tendency,    &
           ighost, ks, ak, bk, cp_vir
      use fv_arrays_mod, only: ps_bp, u_phys, v_phys, t_phys, q_phys
      use mapz_module,    only: te_map, geo_map, benergy, p_energy
      use dyn_core,       only: cd_core

      use pft_module,     only: fftfax, pft2d
      use timingModule,   only: timing_on, timing_off

#ifdef SPMD
      use mod_comm,       only: gid, mp_send4d_ns, mp_recv4d_ns
#endif

      use shr_kind_mod,   only: r8 => shr_kind_r8

#endif
!  implicit none
      use fv_arrays_mod, only: fv_array_check, fv_array_sync, fv_stack_push, &
           fv_stack_reset, fv_print_chksums, jsp, jep, ksp, kep
#ifdef use_shared_pointers
      use fv_arrays_mod, only: is,ie,js,je, isd,ied,jsd,jed, isg,ieg,jsg,jeg, &
           nlev, ncnst, &
           ptr_ps_bp, ptr_u_phys, ptr_v_phys, ptr_t_phys, ptr_q_phys
      implicit none
      real :: ps_bp (isg:ieg, js:je)
      real :: u_phys(isg:ieg, js:je, nlev)
      real :: v_phys(isg:ieg, js:je, nlev)
      real :: t_phys(isg:ieg, js:je, nlev)
      real :: q_phys(isg:ieg, js:je, nlev, ncnst)
      pointer( p_ps_bp, ps_bp )
      pointer( p_u_phys, u_phys )
      pointer( p_v_phys, v_phys )
      pointer( p_t_phys, t_phys )
      pointer( p_q_phys, q_phys )
#else
      use fv_arrays_mod, only: ps_bp, u_phys, v_phys, t_phys, q_phys
      implicit none
#endif

! !INPUT PARAMETERS:
      integer, intent(in):: im       ! dimension in east-west
      integer, intent(in):: jm       ! dimension in North-South
      integer, intent(in):: km       ! number of Lagrangian layers
      integer, intent(in):: jfirst   ! starting latitude index for MPI
      integer, intent(in):: jlast    ! ending latitude index for MPI
      integer, intent(in):: nq       ! total # of tracers to be advected
      integer, intent(in):: pnats    ! number of non-advected tracers
      integer, intent(in):: ndt      ! the large time step in seconds

      real, intent(in):: ptop    ! constant pressure at model top (pascal)
      real, intent(in):: om      ! angular velocity of earth's rotation  
      real, intent(in):: r_vir    ! constant for the virtual effect
      real, intent(in):: cp      ! heat capacity of air at constant pressure
      real, intent(in):: rg
      real, intent(in):: cappa   ! rg/cp
      real, intent(in):: ae      ! radius of the earth (m)

      logical, intent(in):: p_map    ! Partial remapping 
      real,intent(in):: consv    ! energy conserving correcton factor: 0 = correction
!                                     1 = full correction
#ifndef USE_LIMA
      type(time_type), optional, intent(in) :: Time_next
#endif

      real, intent(inout):: phis(im,jfirst-ighost:jlast+ighost) ! surface geopotential (grav*zs)

! !INPUT/OUTPUT
      real, intent(inout):: ps(im,jfirst:jlast)      ! surface pressure (pascal)
      real, intent(inout):: delp(im,jfirst:jlast,km) ! pressure thickness (pascal)

! Ghosted prog arrays:
      real, intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-wind (m/s)
      real, intent(inout):: v(im,jfirst-ng_d:jlast+ng_d,km)  ! v-wind (m/s)
      real, intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d,km) ! temperature (K)
      real, intent(inout):: q(im,jfirst-ng_d:jlast+ng_d,km,nq)
! tracers (e.g., specific humidity)
! tracer mass / moist_air_mass

      real, intent(inout):: pkz(im,jfirst:jlast,km)   ! finite-volume mean of pk

! peln are not needed as input if consv = 0.
      real, intent(inout):: pe(im,km+1,jfirst:jlast)  ! pressure (pascal) at layer edges
      real, intent(inout):: peln(im,km+1,jfirst:jlast)  ! log pressure (pe) at layer edges
      real, intent(inout):: ua(im,jfirst:jlast,km)    ! A grid u-wind (m/s)
      real, intent(inout):: va(im,jfirst:jlast,km)    ! A grid v-wind (m/s)

! !OUTPUT (input values are not used):
      real, intent(inout):: pk(im,jfirst:jlast,km+1)  ! pe**cappa (cappa = rg/cp)
! outputed values can be used by physdrv
      real, intent(out):: omga(im,jfirst:jlast,km)    ! vertical pressure velocity (pa/sec)
! This is the rate of change of the Lagrangian
! interface in pascal/sec unit.

! !DESCRIPTION:
!
! Developer: Shian-Jiann (SJ) Lin, NOAA/GFDL; email: Shian-Jiann.Lin
!
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!               u(i,j+1)
!                 |
!      v(i,j)---delp(i,j)---v(i+1,j)
!                 |
!               u(i,j)
!
! The C grid component is hidden from the user

! External routine required: the user needs to supply a subroutine to set up
!                            "Eulerian vertical coordinate" for remapping purpose.
!                            Currently this routine is named as set_eta()
!                            In principle any terrian following vertical
!                            coordinate can be used. The input to fv_dynamics
!                            need not be on the same vertical coordinate
!                            as the output.
!
! Remarks: values at poles for both u and v need not be defined; but values for
!          all other scalars needed to be defined at both poles (as polar cap mean
!          quantities). Tracer advection is done "off-line" using the
!          large time step. Consistency is maintained by using the time accumulated
!          Courant numbers and horizontal mass fluxes for the FFSL algorithm.
!
!          The user may set the value of n_zonal to optimize the SMP performance
!          The optimal valuse of n_zonal depends on the total number of available
!          shared memory CPUs per node.
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Initial SMP version delivered to Will Sawyer
!   WS  99.10.03:  1D MPI completed and tested; 
!   WS  99.10.11:  Additional documentation
!   WS  99.10.28:  benergy and te_map added; arrays pruned
!   SJL 00.01.01:  SMP and MPI enhancements; documentation
!   WS  00.07.13:  Changed PILGRIM API
!   SJL 03.11.21:  3rd major revision
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Local variables:

      integer i, j, k
      integer it, padv, iq

      real :: umax = 300.   ! estimated upper bound of the maximum u-wind (m/s)
      real  dt
      real  bdt
      real  kappa, cp_dyn
      real frac       ! tracer time split fraction

! Local auto arrays
      real  wk1(im+2, jfirst:jlast+2)
      real  wk2(im+2, jfirst:jlast+2)
      real  tte(jfirst:jlast)
!temporary work arrays
      real :: cx(im,jfirst-ng_d:jlast+ng_d,km)
      real :: mfx(im,jfirst:jlast,km)
      real :: cy(im,jfirst:     jlast+1,   km)
      real :: mfy(im,jfirst:     jlast+1,   km)
      real :: dpt(im,jfirst-1:   jlast+1,   km)
      real :: vc(im,jfirst-2:   jlast+2,   km)
      real :: delpf(im,jfirst-ng_d:jlast+ng_d,km)
      real :: uc(im,jfirst-ng_d:jlast+ng_d,km)
      real :: dwz(im,jfirst-1:   jlast,     km+1)
      real :: pkc(im,jfirst-1:   jlast+1,   km+1) 
      real :: wz(im,jfirst-1:   jlast+1,   km+1)
! Save a copy of pe for the computation of vertical velocity

      real :: pem(im,km+1,jfirst:jlast)
#ifdef use_shared_pointers
      pointer( p_cx, cx )
      pointer( p_mfx, mfx )
      pointer( p_cy, cy )
      pointer( p_mfy, mfy )
      pointer( p_dpt, dpt )
      pointer( p_vc, vc )
      pointer( p_delpf, delpf )
      pointer( p_uc, uc )
      pointer( p_dwz, dwz )
      pointer( p_pkc, pkc )
      pointer( p_pem, pem )
      pointer( p_wz, wz )
      pointer( p_tte, tte )
#endif

      integer js2g0, jn2g0

#ifdef SW_DYN
! The following code segment is for shallow water test cases
      real yy, tday, aoft, pi, gv, ssec
      data ssec /0./

      kappa  = 1.
      cp_dyn = 1.
#else
      kappa  = cappa
      cp_dyn = cp
#endif

      js2g0 = max(2,jfirst)
      jn2g0 = min(jm-1,jlast) 
      padv = nq - pnats
#ifdef use_shared_pointers
      call fv_print_chksums( 'Entering  fv_dynamics' )
      p_ps_bp =  ptr_ps_bp
      p_u_phys =  ptr_u_phys
      p_v_phys =  ptr_v_phys
      p_t_phys =  ptr_t_phys
      p_q_phys =  ptr_q_phys
      call fv_stack_reset !zeroes stack from previous invocation
      call fv_stack_push( p_cx,    im*km*(jlast-jfirst+2*ng_d+1) )
      call fv_stack_push( p_mfx,   im*km*(jlast-jfirst+1)        )
      call fv_stack_push( p_cy,    im*km*(jlast-jfirst+2)        )
      call fv_stack_push( p_mfy,   im*km*(jlast-jfirst+2)        )
      call fv_stack_push( p_dpt,   im*km*(jlast-jfirst+3)        )
      call fv_stack_push( p_vc,    im*km*(jlast-jfirst+5)        )
      call fv_stack_push( p_delpf, im*km*(jlast-jfirst+2*ng_d+1) )
      call fv_stack_push( p_uc,    im*km*(jlast-jfirst+2*ng_d+1) )
      call fv_stack_push( p_dwz,   im*(km+1)*(jlast-jfirst+2)    )
      call fv_stack_push( p_pkc,   im*(km+1)*(jlast-jfirst+3)    )
      call fv_stack_push( p_wz,    im*(km+1)*(jlast-jfirst+3)    )
      call fv_stack_push( p_pem,   im*(km+1)*(jlast-jfirst+1)    )
      call fv_stack_push( p_tte,              jlast-jfirst+1     )
#endif
      call fv_array_check( LOC(phis) )
      call fv_array_check( LOC(ps) )
      call fv_array_check( LOC(delp) )
      call fv_array_check( LOC(u) )
      call fv_array_check( LOC(v) )
      call fv_array_check( LOC(pt) )
      call fv_array_check( LOC(q) )
      call fv_array_check( LOC(pkz) )
      call fv_array_check( LOC(pe) )
      call fv_array_check( LOC(peln) )
      call fv_array_check( LOC(ua) )
      call fv_array_check( LOC(va) )
      call fv_array_check( LOC(pk) )
      call fv_array_check( LOC(omga) )
      call fv_array_check( LOC(ps_bp) )
      call fv_array_check( LOC(u_phys) )
      call fv_array_check( LOC(v_phys) )
      call fv_array_check( LOC(t_phys) )
      call fv_array_check( LOC(q_phys) )

#ifdef SPMD
      call fv_array_sync()
      call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_s, u)
      call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, v)
#endif

      if ( use_tendency ) then
!$omp parallel do private(i,j,k,iq)
          do k=ksp,kep !1,km
             do j=jfirst,jlast
                do i=1,im
                   u_phys(i,j,k) = ua(i,j,k)
                   v_phys(i,j,k) = va(i,j,k)
                   t_phys(i,j,k) = pt(i,j,k)
                enddo
             enddo

             do iq=1,nq
                do j=jfirst,jlast
                   do i=1,im
                      q_phys(i,j,k,iq) = q(i,j,k,iq)
                   enddo
                enddo
             enddo
          enddo
      endif

! Mass conservation for non-advected tracers
      if(pnats /= 0) then
          do iq=padv+1,nq
!$omp parallel do private(i,j,k)
             do k=ksp,kep !1,km
                do j=jfirst,jlast
                   do i=1,im
                      q(i,j,k,iq) = q(i,j,k,iq)*delp(i,j,k)
                   enddo
                enddo
             enddo
          enddo
      endif

!$omp parallel do private(i,j,k)
      do k=ksp,kep !1,km
! Initialize the CFL number accumulators: (cx, cy)
! Initialize total mass fluxes: (mfx, mfy)
         if( padv > 0 ) then
             do j=jfirst-ng_d,jlast+ng_d
                do i=1,im
                   cx(i,j,k) = 0.
                enddo
             enddo
             do j=jfirst,jlast
                do i=1,im
                   cy(i,j,k) = 0.
                   mfx(i,j,k) = 0.
                   mfy(i,j,k) = 0.
                enddo
             enddo
         endif
      enddo        ! end parallel k-loop
      call fv_array_sync()
!      call fv_print_chksums( 'after cx cy mfx mfy' )

      if ( kord_tm > 0 ) then

!$omp parallel do private (i, j, k)
          do j=jsp,jep !jfirst,jlast
!    do k=1,km+1
!       do i=1,im
!          pem(i,k,j) = pe(i,k,j)
!       enddo
!    enddo

             do k=1,ks+1
                do i=1,im
                   pem(i,k,j) = ak(k)
                enddo
             enddo

             do k=ks+2,km+1
                do i=1,im
                   pem(i,k,j) = ak(k) + bk(k) * ps_bp(i,j)
                enddo
             enddo
          enddo
          call fv_array_sync()
      endif

#ifdef SPMD
      call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_s, u)
      u_ghosted = .true.
      call fv_array_sync()
#endif

      if( km > 1 .and. consv > 1.E-5 ) then

!---------------------------------------------------------------------
! Compute globally integrated Total Energy
! Assuming the 1st tracer is specific humidity
!---------------------------------------------------------------------

          if ( kord_tm > 0 ) then
              call benergy(im, jm, km, u, v, pt, delp, q(1,jfirst-ng_d,1,1),  &
                   pe, peln, phis(1,jfirst),  u_ghosted,              &
                   ng_d, ng_s, r_vir, cp_vir, cp, rg,  tte,           &
                   jfirst, jlast, acap, cosp, ua, va)
          else
              call p_energy(im, jm, km, u, v, pt, delp, q(1,jfirst-ng_d,1,1),  &
                   pe, peln, phis(1,jfirst),                           &
                   ng_d, ng_s, r_vir, cp_vir, cp, rg,  tte,            &
                   jfirst, jlast, acap, cosp, ua, va)
          endif
      endif
!      call fv_print_chksums( 'after benergy penergy' )

!$omp parallel do private(i, j, k, wk1, wk2)
      do k=ksp,kep !1,km
         do j=jfirst,jlast
            do i=1,im
!               Save initial delp field before the small-time-step
               va(i,j,k) = delp(i,j,k)
               delpf(i,j,k) = delp(i,j,k)
#ifdef SW_DYN
               pt(i,j,k) = 1.
#else
! Convert pt to virtual potential temperature * CP
               pt(i,j,k) = cp*pt(i,j,k)/pkz(i,j,k)*(1.+r_vir*q(i,j,k,1))
#endif
            enddo
         enddo

#ifndef USE_LIMA
         call pft2d( delpf(1,js2g0,k), im, jn2g0-js2g0+1, wk1, wk2, 1) 
#else               
         call pft2d( delpf(1,js2g0,k),  sc(js2g0), dc(1,js2g0),   &
              im, jn2g0-js2g0+1, ifax, trigs, wk1, wk2 )
#endif

      enddo !k-loop
      call fv_array_sync()
!      call fv_print_chksums( 'after delpf' )


#ifdef SPMD
      call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, v)

#if !defined (SW_DYN)                                                                 
      call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, pt)      
#endif                                                                                
      call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, delpf)   
#if !defined (SW_DYN)                                                                 
      call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, pt)      
#endif                                                                                
      call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng_d, ng_d, delpf)   
      call fv_array_sync()
#endif

      bdt = ndt
      dt  = bdt / real(n_split)

      do it=1, n_split

#ifdef SW_DYN                                                                 
         pi = 4.*atan(1.)
         gv = 9.8*720.

         if (icase == 8 ) then
             tday = ssec/(24.*3600.)
             if(tday >= 20.) then
                 ssec = 0                  ! reset to 0
                 aoft = 0.
             elseif(tday <= 4.) then
                 aoft = 0.5*(1.-cos(0.25*pi*tday))
             elseif(tday <= 16.) then
                 aoft = 1.
             else
                 aoft = 0.5*(1.+cos(0.25*pi*(tday-16.)))
             endif

             aoft = aoft * gv

! Time invariant part of the terrain BC
             do j=jfirst, jlast
                if(sinp(j) > 0.) then
                    yy = (cosp(j)/sinp(j))**2
                    do i=1,im
                       phis(i,j) =  aoft*yy*exp(1.-yy)*sinlon(i)
                    enddo
                else
                    do i=1,im
                       phis(i,j) = 0.
                    enddo
                endif
             enddo
             ssec = ssec + dt
         endif
#endif

!--------------------------------------------------------
! Call the Lagrangian dynamical core using small tme step
!--------------------------------------------------------

         call timing_on('CD_CORE')


#ifndef USE_LIMA
         call cd_core(im,  jm,  km,  padv, n_zonal, jfirst, jlast, &
              u,   v,  pt,  delp, pe,  pk, n_split, dt,    &
              ptop, umax,   ae, rcap, cp_dyn, kappa,       &
              iord_mt, iord_tm, jord_mt, jord_tm,          &
              ng_d,   ng_s,  it, n_split,                  &
              om, phis, sinp, cosp, cose, acosp,           &
              sinlon, coslon, cosl5, sinl5, f_d, dx,       &
              one_by_dx, one_by_dy, cx  , cy , mfx, mfy,   &
              delpf, uc, vc, dpt, pkz, dwz, pkc, wz,       &
              omga,  master, ighost, dyn_step)

         if (dyn_step == 1 ) then
             dyn_step = 2
         elseif(dyn_step == 2 ) then
             dyn_step = 1
         endif

#else               
         call cd_core(im,  jm,  km,  padv, n_zonal, jfirst, jlast, &
              u,   v,  pt,  delp, pe,  pk, n_split, dt,    &
              ptop, umax,   ae, rcap, cp_dyn, kappa,       &
              iord_mt, iord_tm, jord_mt, jord_tm,          &
              ng_d,   ng_s,  it, n_split,                  &
              om, phis, sinp, cosp, cose, acosp,           &
              sinlon, coslon, cosl5, sinl5, f_d, dx,       &
              one_by_dx, one_by_dy, cx  , cy , mfx, mfy,   &
              delpf, uc, vc, dpt,                          &
              pkz, dwz, pkc, wz, omga,  master, ighost )

#endif

         call timing_off('CD_CORE')

      enddo

      if(padv /= 0) then

!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
!--------------------------------------------------------

          call timing_on('tracer_2d')
#ifndef USE_LIMA
          call tracer_2d(q, padv, va, cx, cy, mfx, mfy, iord_tr, jord_tr,    &
               ng_d, cose, cosp, acosp, acap, rcap, tra_step,      &
               fill, im, jm, km, jfirst, jlast, pkz, bdt, frac)
#else
          call tracer_2d_lima(q, padv, va, cx, cy, mfx, mfy, iord_tr, jord_tr,  &
               ng_d, cose, cosp, acosp, acap, rcap, fill,        &
               im, jm, km, jfirst, jlast, pkz, ua, bdt, frac)
#endif
          call timing_off('tracer_2d')


      endif
#ifdef SW_DYN
      do j=jfirst,jlast
         do i=1,im
            ps(i,j) = delp(i,j,1)
         enddo
      enddo
#endif

      if ( km > 2 ) then

!------------------------------------------------------------------------
! Peroform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------

          call timing_on('Remapping')

          if ( kord_tm < 0 ) then
              call geo_map(ps, omga,  pe, delp, pkz, pk, ndt,             &
                   im, jm, km, n_zonal, jfirst, jlast, padv,      &
                   u,   v,  pt,  q, r_vir, cp_vir, cp,  cappa,    &
                   abs(kord_mt), abs(kord_tr), abs(kord_tm),      &
                   peln, phis(1,jfirst),  ng_d, ng_s,  &
                   rg, acap, cosp, ua, va, consv, tte, fill)
          else
#ifndef USE_LIMA
              call te_map(consv, ps, omga,  pe, pem, delp, pkz, pk, ndt,    &
                   im, jm, km, p_map, n_zonal, jfirst, jlast, padv,  &
                   u,   v,  pt,  q, phis(1,jfirst), r_vir, cp_vir,   &
                   cp, cappa, kord_mt, kord_tr, kord_tm, peln, tte,  &
                   ng_d, ng_s, ua, va, fill, Time_next)
#else
              call te_map(consv, ps, omga,  pe, pem, delp, pkz, pk, ndt,    &
                   im, jm, km, p_map, n_zonal, jfirst, jlast, padv,  &
                   u,   v,  pt,  q, phis(1,jfirst), r_vir, cp_vir,   &
                   cp, cappa, kord_mt, kord_tr, kord_tm, peln, tte,  &
                   ng_d, ng_s, ua, va, fill)
#endif
          endif

          call timing_off('Remapping')
      endif

!-----------------------------------------------------
! Add the v*grad(p) term to omega (dp/dt) for physics
!-----------------------------------------------------
      call compute_vdot_gradp(im, jm, km, jfirst, jlast, ng_d,       &
           bdt, frac, cx, cy, pe, omga)
#ifdef OMGA_PHYS
#endif
      if(pnats /= 0) then
! Mass conservation for non-advected tracers
          do iq=padv+1,nq
!$omp parallel do private(i,j,k)
             do k=ksp,kep !1,km
                do j=jfirst,jlast
                   do i=1,im
                      q(i,j,k,iq) = q(i,j,k,iq) / delp(i,j,k)
                   enddo
                enddo
             enddo
          enddo
      endif

!EOC
      call fv_print_chksums( 'Exiting  fv_dynamics' )
    end subroutine fv_dynamics

!-----------------------------------------------------------------------

    subroutine compute_vdot_gradp(im, jm, km, jfirst, jlast, ng_d,       &
         dt, frac, cx, cy, pe, omga)

      use pft_module,     only: pft2d
#ifdef SPMD
      use mod_comm, only: gid,  mp_send3d, mp_recv3d
#endif
      use fv_arrays_mod, only: fv_array_check, fv_array_sync, fv_stack_push, &
           fv_print_chksums, ksp, kep
#ifndef USE_LIMA
#ifdef use_shared_pointers
      use fv_arrays_mod, only: isg,ieg,nlev
      use fv_arrays_mod, only: ptr_pesouth
      implicit none
      real :: pesouth(isg:ieg, nlev+1)
      pointer( p_pesouth, pesouth )
#else
      use fv_arrays_mod, only: pesouth
      implicit none
#endif
#else               
      use fv_pack,        only: ifax, trigs, sc, dc
      use fv_arrays_mod, only: pesouth
      implicit none
#endif

! !INPUT PARAMETERS:
      integer, intent(in):: im       ! dimension in east-west
      integer, intent(in):: jm       ! dimension in North-South
      integer, intent(in):: km       ! number of Lagrangian layers
      integer, intent(in):: jfirst   ! starting latitude index for MPI
      integer, intent(in):: jlast    ! ending latitude index for MPI
      integer, intent(in):: ng_d

      real, intent(in):: dt
      real, intent(in):: frac

      real, intent(in):: cx(im,jfirst-ng_d:jlast+ng_d,km) 
      real, intent(in):: cy(im,jfirst:     jlast+1,   km)
      real, intent(in):: pe(im,km+1,jfirst:jlast)     ! pressure (pascal) at layer edges
      real, intent(inout):: omga(im,jfirst:jlast,km)  ! vertical pressure velocity (pa/sec)

! Local 
      real pm(im, jfirst-1:jlast+4)
      real penorth(im, km+1)
      real grad(im, jfirst:jlast+4)
      real fac, sum1

      integer i,j,k, js2g0, jn2g0
#ifdef use_shared_pointers
      pointer( p_penorth, penorth )
      call fv_stack_push( p_penorth, im*(km+1) )
      p_pesouth = ptr_pesouth
#endif

      call fv_print_chksums( 'Entering  compute_vdot_gradp' )
      js2g0 = max(2,jfirst)
      jn2g0 = min(jm-1,jlast) 

      fac = 0.5 / (dt * frac)

      call fv_array_check( LOC(cx) )
      call fv_array_check( LOC(cy) )
      call fv_array_check( LOC(pe) )
      call fv_array_check( LOC(omga) )

#ifdef SPMD
      call mp_send3d(gid+1, gid-1, im, km+1, jm, 1, im, 1, km+1,               &
           jfirst, jlast, 1, im, 1, km+1, jlast, jlast, pe)
      call mp_send3d(gid-1, gid+1, im, km+1, jm, 1, im, 1, km+1,               &
           jfirst, jlast, 1, im, 1, km+1, jfirst, jfirst, pe)
      call mp_recv3d(gid-1, im, km+1, jm, 1, im, 1, km+1, 1, 1, &
           1, im,    1, km+1, 1, 1, pesouth)
      call mp_recv3d(gid+1, im, km+1, jm, 1, im, 1, km+1, 1, 1, &
           1, im,    1, km+1, 1, 1, penorth)
      call fv_array_sync()
#endif

!$omp parallel do private(i,j,k,pm,grad, sum1)
      do k=ksp,kep !1,km

! Compute layer mean p
         do j=jfirst,jlast
            do i=1,im
               pm(i,j) = 0.5 * ( pe(i,k,j) + pe(i,k+1,j) )
            enddo
         enddo

         if ( jfirst/=1 ) then
             do i=1,im
                pm(i,jfirst-1) = 0.5 * ( pesouth(i,k) + pesouth(i,k+1))
             enddo
         endif

         if ( jlast/=jm ) then
             do i=1,im
                pm(i,jlast+1) = 0.5 * ( penorth(i,k) + penorth(i,k+1))
             enddo
         endif

         do j=js2g0,jn2g0
            i=1
            grad(i,j) = fac * cx(i,j,k) * (pm(i,j)-pm(im,j)) 
            do i=2,im
               grad(i,j) = fac * cx(i,j,k) * (pm(i,j)-pm(i-1,j)) 
            enddo
         enddo

         do j=js2g0,jn2g0
            do i=1,im-1
               omga(i,j,k) = omga(i,j,k) + grad(i,j) + grad(i+1,j)
            enddo
            i=im
            omga(i,j,k) = omga(i,j,k) + grad(i,j) + grad(1,j)
         enddo

         do j=js2g0,min(jm,jlast+1)
            do i=1,im
               grad(i,j) = fac * cy(i,j,k) * (pm(i,j)-pm(i,j-1)) 
            enddo
         enddo

         do j=js2g0,jn2g0
            do i=1,im
               omga(i,j,k) = omga(i,j,k) + grad(i,j) + grad(i,j+1)
            enddo
         enddo

!-------------------
! Apply polar filter
!-------------------

#ifndef USE_LIMA
         call pft2d( omga(1,js2g0,k), im, jn2g0-js2g0+1, pm, grad, 1) 
#else               
         call pft2d( omga(1,js2g0,k),  sc(js2g0), dc(1,js2g0),   &
              im, jn2g0-js2g0+1, ifax, trigs, pm, grad)
#endif


! Note: Since V*grad(P) at poles are harder to compute accurately we use the average of sourding points
!       to be used as input to physics.

         if ( jfirst==1 ) then
             sum1 = 0.
             do i=1,im
                sum1 = sum1 + omga(i,2,k)
             enddo
             sum1 = sum1 / real(im)
             do i=1,im
                omga(i,1,k) = sum1
             enddo
         endif

         if ( jlast==jm ) then
             sum1 = 0.
             do i=1,im
                sum1 = sum1 + omga(i,jm-1,k)
             enddo
             sum1 = sum1 / real(im)
             do i=1,im
                omga(i,jm,k) = sum1
             enddo
         endif
      enddo
      call fv_array_sync()
      call fv_print_chksums( 'Exiting  compute_vdot_gradp' )
    end subroutine compute_vdot_gradp

end module fv_dynamics_mod
