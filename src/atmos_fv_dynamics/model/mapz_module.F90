module mapz_module

  use fv_pack,      only: ak, bk, ks, ptop, cosp, acap,     &
                         u_ghosted, coslon, sinlon,          &
                         cose,  acosu, ptop_min, d2a3d, master &
                        ,  beglon, endlon, beglat, endlat

  use fill_module,  only: fillz

#ifdef VERT_FLUX
  use fv_diagnostics,   only: id_dnflux, id_upflux
  use diag_manager_mod, only: send_data
#endif

#ifdef SPMD
      use mod_comm, only: gid, mp_send3d, mp_recv3d
#endif
      use time_manager_mod, only: time_type      ! This is not needed if not under FMS
      use fv_arrays_mod, only: fv_array_check, fv_array_limits, &
           fv_print_chksums, fv_stack_push, fv_array_sync, fv_print_chksum, &
           isp, iep, jsp, jep, ksp, kep

      implicit none

      private
! !PUBLIC MEMBER FUNCTIONS:
      public benergy, geo_map, te_map, p_energy

#ifndef SPMD
      integer :: gid=0
#endif

CONTAINS

  subroutine te_map(consv, ps, omga, pe, pem, delp, pkz, pk, mdt,        &
       im, jm, km, p_map, nx, jfirst, jlast, nq, u, v,   &
       pt,  q, hs, r_vir, cp_vir, cp, akap, kord_mt,     &
       kord_tr, kord_tm,  peln,                          &
       tte0, ng_d,  ng_s, ua, va, fill, Time_next)
#ifdef use_shared_pointers
    use fv_arrays_mod, only: is,ie,js,je, isd,ied,jsd,jed, isg,ieg,jsg,jeg
    use fv_arrays_mod, only: nlev
    use fv_arrays_mod, only: ptr_ps_bp, ptr_pesouth
    real :: ps_bp (isg:ieg, js:je)
    pointer( p_ps_bp, ps_bp )
    real :: pesouth(isg:ieg, nlev+1)
    pointer( p_pesouth, pesouth )
#else
    use fv_arrays_mod, only: ps_bp, pesouth
#endif


! !ROUTINE: te_map --- Perform (optionally partial) remapping from surface to
!                      a chosen level


! !INPUT PARAMETERS:
    integer, intent(in):: mdt                   ! mapping time step (same as phys)
    integer, intent(in):: im, jm, km            ! x, y, z dimensions
    integer, intent(in):: nq                    ! number of tracers (including h2o)
    integer, intent(in):: nx                    ! number of SMP "decomposition" in x
    integer, intent(in):: ng_d
    integer, intent(in):: ng_s
    integer, intent(in):: jfirst, jlast         ! starting & ending latitude index
    integer, intent(in):: kord_mt                  ! Mapping oder for the vector winds
    integer, intent(in):: kord_tr                  ! Mapping oder for tracers
    integer, intent(in):: kord_tm                  ! Mapping oder for thermodynamics

    real, intent(in):: consv                 ! factor for TE conservation
    real, intent(in)::  r_vir
    real, intent(in):: cp_vir
    real, intent(in):: cp
    real, intent(in):: tte0(jfirst:jlast)
    real, intent(in):: hs(im,jfirst:jlast)  ! surface geopotential

    logical, intent(in):: p_map                 ! Do partial remapping
! layers between [1,ks] would stay Lagrangian
    logical, intent(in):: fill                  ! fill negative tracers
    type(time_type), optional, intent(in) :: Time_next

! !INPUT/OUTPUT
    real, intent(inout):: pk(im,jfirst:jlast,km+1) ! pe to the kappa
!Balaji: no assumed size?  real, intent(inout):: q(im,jfirst-ng_d:jlast+ng_d,km,*)
    real, intent(inout):: q(im,jfirst-ng_d:jlast+ng_d,km,nq)
    real, intent(inout):: delp(im,jfirst:jlast,km) ! pressure thickness
    real, intent(inout)::  pe(im,km+1,jfirst:jlast) ! pressure at layer edges
    real, intent(inout):: pem(im,km+1,jfirst:jlast) ! (old) pressure at layer edges
    real, intent(inout):: ps(im,jfirst:jlast)      ! surface pressure

! u-wind will be ghosted one latitude to the north upon exit
    real, intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)   ! u-wind (m/s)
    real, intent(inout):: v(im,jfirst-ng_d:jlast+ng_d,km)   ! v-wind (m/s)
    real, intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d,km)  ! virtual potential temperature 
! as input; output: temperature

    real, intent(inout):: ua(im,jfirst:jlast,km)   ! u-wind (m/s) on physics grid
    real, intent(inout):: va(im,jfirst:jlast,km)   ! v-wind (m/s) on physics grid
    real, intent(out):: omga(im,jfirst:jlast,km)    ! vertical press. velocity (pascal/sec)
    real, intent(out):: peln(im,km+1,jfirst:jlast)  ! log(pe)
    real, intent(out):: pkz(im,jfirst:jlast,km)     ! layer-mean pk for converting t to pt

! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!EOP
!-----------------------------------------------------------------------

#ifdef VERT_FLUX
    real :: mflux(im,jfirst:jlast,km)
    logical used
#endif
    real rmin(nx*jm), rmax(nx*jm)
    real tte(jfirst:jlast)
! x-y
    real  u2(im,jfirst:jlast+1)
    real  v2(im,jfirst:jlast)
! y-z
    real   q2(im/nx,km)
    real  dp2(im/nx,km)
    real  pe0(im,km+1)
    real  pe1(im,km+1)
    real  pe2(im,km+1)
    real  pe3(im,km+1)
    real phis(im,km+1)
! x
    real    gz(im)
    real ratio(im)
    real   bte(im)

    integer i, j, k, js2g0, jn2g0
    integer kmap

    real akap, dak, bkh, qmax, qmin
    real te_sp,  te_np
    real xsum, ysum, tsum
    real zsum(im)
    real dtmp
    real rdt5
    real rg
    real tm
    real te1
    real dlnp

    integer ixj, jp, it, i1, i2
    integer iq
    logical:: diag = .false.

    integer :: ixjs, ixje
#ifdef use_shared_pointers
    pointer( p_tte, tte )
#ifdef VERT_FLUX
    pointer( p_mflux, mflux )
    call fv_stack_push( p_mflux, im*km*(jlast-jfirst+1) )
#endif
    call fv_stack_push( p_tte, jlast-jfirst+1 )
    p_ps_bp = ptr_ps_bp
    p_pesouth = ptr_pesouth
#endif

    call fv_print_chksums( 'Entering  te_map' )
!     diag = .true.
    call fv_array_check( LOC(hs) )
    call fv_array_check( LOC(pk) )
    call fv_array_check( LOC(q) )
    call fv_array_check( LOC(delp) )
    call fv_array_check( LOC(pe) )
    call fv_array_check( LOC(pem) )
    call fv_array_check( LOC(ps) )
    call fv_array_check( LOC(u) )
    call fv_array_check( LOC(v) )
    call fv_array_check( LOC(pt) )
    call fv_array_check( LOC(ua) )
    call fv_array_check( LOC(va) )
    call fv_array_check( LOC(omga) )
    call fv_array_check( LOC(peln) )
    call fv_array_check( LOC(pkz) )

    if ( p_map ) then
        kmap = ks+1
    else
        kmap = 1
    endif

#ifdef SPMD
    call mp_send3d(gid-1, gid+1, im, jm, km, 1, im, jfirst-ng_d, jlast+ng_s, &
         1, km, 1, im, jfirst, jfirst, 1, km, u)
    call mp_send3d(gid+1, gid-1, im, km+1, jm, 1, im, 1, km+1,               &
         jfirst, jlast, 1, im, 1, km+1, jlast, jlast, pe)
#endif

    js2g0 = max(2,   jfirst)
    jn2g0 = min(jm-1,jlast)

    call pkez(nx, im, km, jfirst, jlast,             &
         pe, pk, akap, ks, peln, pkz, .false.)

#ifdef SPMD
    call mp_recv3d(gid+1, im, jm, km, 1, im, jfirst-ng_d, jlast+ng_s, &
         1, km, 1, im, jlast+1, jlast+1, 1, km, u)
    call fv_array_sync()
#endif

!$omp parallel do private(i, j, k, u2, v2, te_sp, te_np)
    do 1000 k=ksp,kep !1,km
! Compute cp*T + KE

       do j=js2g0,min(jlast+1,jm)    ! u ghosted one lat to the North
          do i=1,im
             u2(i,j) = cose(j) * u(i,j,k)**2
          enddo
       enddo

       do j=js2g0,jn2g0
          do i=1,im
             v2(i,j) = v(i,j,k)**2
          enddo

          do i=1,im-1
             ua(i,j,k) = 0.25 * ( (u2(i,j)+u2(i,j+1))*acosu(j) +   &
                  v2(i,j) + v2(i+1,j)  ) +         &
                  pt(i,j,k)*pkz(i,j,k)
!                        pt(i,j,k)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,1))  &
!                       * (1.+cp_vir*q(i,j,k,1))
          enddo
! i=im
          ua(im,j,k) = 0.25 * ( (u2(im,j)+u2(im,j+1))*acosu(j) +   &
               v2(im,j) + v2(1,j)  ) +            &
               pt(im,j,k)*pkz(im,j,k)
!                       * (1.+cp_vir*q(im,j,k,1))
!                       * (1.+cp_vir*q(im,j,k,1))
       enddo

       if( jfirst == 1 ) then
           te_sp = 0.
           do i=1,im
!            te_sp = te_sp + u2(i,2) + v2(i,2)
! SJL: 01/26/04
              te_sp = te_sp + 0.5*u2(i,2)/cose(2)
           enddo
           te_sp = 0.5*te_sp/float(im) +  pt(1,1,k)*pkz(1,1,k)
!                    * (1.+cp_vir*q(1,1,k,1)) / (1.+r_vir*q(1,1,k,1))
           do i=1,im
              ua(i,1,k) = te_sp
           enddo
       endif

       if ( jlast == jm ) then
           te_np = 0.
           do i=1,im
!            te_np = te_np + u2(i,jm) + v2(i,jm-1)
! SJL: 01/26/04
              te_np = te_np + 0.5*u2(i,jm)/cose(jm)
           enddo
           te_np = 0.5*te_np/float(im) + pt(1,jm,k)*pkz(1,jm,k)
!              * (1.+cp_vir*q(1,jm,k,1)) / (1.+r_vir*q(1,jm,k,1))

           do i=1,im
              ua(i,jm,k) = te_np
           enddo
       endif

! Compute va; geopotential increments
       do j=jfirst,jlast
          do i=1,im
             va(i,j,k) = pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
          enddo
       enddo
1000 end do !continue
    call fv_array_sync()

#ifdef SPMD
    call mp_recv3d(gid-1, im, km+1, jm, 1, im, 1, km+1, jfirst-1, jfirst-1, &
         1, im,    1, km+1, jfirst-1, jfirst-1, pesouth)
    call fv_array_sync()
#endif

    it = im / nx
    jp = nx * ( jlast - jfirst + 1 )

    call fv_array_limits( 1, jp, ixjs, ixje )
!$omp parallel do default(shared)                 &
!$omp private(i, j, k, pe0, pe1, pe2, pe3, ratio) &
!$omp private(dak,bkh,rdt5,phis,ixj,i1,i2, iq, q2, dp2)
    do 2000 ixj=ixjs, ixje !1, jp

       j  = jfirst + (ixj-1) / nx
       i1 = 1 + it * mod(ixj-1, nx)
       i2 = i1 + it - 1

! Copy data to local 2D arrays.
       do k=1,km+1
          do i=i1,i2
             pe1(i,k) = pe(i,k,j)
          enddo
       enddo

       rdt5 = 0.5 / float(mdt)
#ifdef VERT_FLUX
! use pe0 for net flux; pkz for downward (positive) flux; mflux for upward flux
       do k=1,km+1
          do i=i1,i2
             pe0(i,k) = pe1(i,k) - (ak(k)+bk(k)*pe1(i,km+1))
          enddo
       enddo
       do k=1,km
          do i=i1,i2
             dak = rdt5*(pe0(i,k)+pe0(i,k+1))
             if ( dak < 0. ) then
                 mflux(i,j,k) = dak
                 pkz(i,j,k) = 0.
             else
                 mflux(i,j,k) = 0.
                 pkz(i,j,k) = dak
             endif
          enddo
       enddo
#endif

       if ( p_map ) then
           do k=1,kmap
              do i=i1,i2
                 pe2(i,k) = pe1(i,k)     ! Stay Lagrangian
              enddo
           enddo
           do k=kmap+1,km
              do i=i1,i2
                 pe2(i,k) = (1.-bk(k))*pe1(i,k) + bk(k)*(ak(k)+bk(k)*pe1(i,km+1))
              enddo
           enddo
       else
           do k=1,ks+1
              do i=i1,i2
                 pe2(i,k) = ak(k)
              enddo
           enddo
           do k=ks+2,km
              do i=i1,i2
                 pe2(i,k) = ak(k) + bk(k)*pe1(i,km+1)
              enddo
           enddo
       endif

       do i=i1,i2
          pe2(i,km+1) = pe1(i,km+1)
       enddo

! Compute omga (dp/dt)
       do k=2,km+1
          do i=i1,i2
             pe0(i,k) = pe1(i,k) - pem(i,k,j)
          enddo
       enddo

       do i=i1,i2
! update ps
          ps(i,j)    = pe1(i,km+1)
          ps_bp(i,j) = ps(i,j)
          omga(i,j,1)  = rdt5 * pe0(i,2)
       enddo

       do k=2,km
          do i=i1,i2
             omga(i,j,k) = rdt5 * ( pe0(i,k) + pe0(i,k+1) )
          enddo
       enddo

       if( ks /= 0 ) then
           if ( .not.  p_map ) then
               do k=1,ks
                  dak = ak(k+1) - ak(k)
                  do i=i1,i2
                     delp(i,j,k) = dak
                  enddo
               enddo
           endif
       endif

       do k=ks+1,km
          do i=i1,i2
             delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
          enddo
       enddo

!---------------------
! Compute Total Energy
!---------------------
       do i=i1,i2
          phis(i,km+1) = hs(i,j)      
       enddo

       do k=km,1,-1
          do i=i1,i2
             phis(i,k) = phis(i,k+1) + va(i,j,k)   
          enddo
       enddo

       do k=1,km+1
          do i=i1,i2
             phis(i,k) = phis(i,k) * pe1(i,k)
          enddo
       enddo

       do k=1,km
          do i=i1,i2
! ua is used as the total energy
             ua(i,j,k) = ua(i,j,k)+(phis(i,k+1)-phis(i,k))/(pe1(i,k+1)-pe1(i,k))
          enddo
       enddo

!----------------
! Map Total Energy
!----------------

! kord=7 combined with high vertical resolutiona in the stratosphere
! may produce insufficient vertical dampping; inertial instability like oscillations in
! v-wind  may occur
! Therefore, use kord=7 only for low vertical resolution cases.

       call map1_ppm (km,   pe1,   ua,  kmap,                &
            km,   pe2,   ua,  0,     0,            &
            im,   i1, i2, j, jfirst, jlast, 1, kord_tm)

!----------------
! Map constituents
!----------------
       if(nq /= 0) then

           do k=1,km
              do i=i1,i2
                 dp2(i-i1+1,k) = pe2(i,k+1) - pe2(i,k)
              enddo
           enddo

!------------------------------------------------------------------
! Do remapping one tracer at a time; seems to be faster on the SGI
! It requires less memory than mapn_ppm

           do iq=1,nq

              call map1_q2(km, pe1, q(1,jfirst-ng_d,1,iq), kmap,  &
                   km, pe2, q2, im,    &
                   i1, i2, 0, kord_tr, j, jfirst, jlast, ng_d)

              if (fill) call fillz(i2-i1+1, km, 1, q2, dp2)
              do k=1,km
                 do i=i1,i2
                    q(i,j,k,iq) = q2(i-i1+1,k)
                 enddo
              enddo
           enddo
       endif

!------
! map u
!------
       if(j /= 1) then

           do i=i1,i2
              pe0(i,1) = pe1(i,1)
           enddo

           if (j > jfirst) then

               do k=2,km+1
                  do i=i1,i2
                     pe0(i,k) = 0.5*(pe1(i,k)+pe(i,k,j-1))
                  enddo
               enddo

               do k=ks+2,km+1
                  bkh = 0.5*bk(k)
                  do i=i1,i2
                     pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pe(i,km+1,j-1))
                  enddo
               enddo

#ifdef SPMD
           else
               do k=2,km+1
                  do i=i1,i2
                     pe0(i,k) = 0.5*(pe1(i,k)+pesouth(i,k))
                  enddo
               enddo

               do k=ks+2,km+1
                  bkh = 0.5*bk(k)
                  do i=i1,i2
                     pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pesouth(i,km+1))
                  enddo
               enddo
#endif
           endif

           if ( p_map ) then
               do k=1,kmap
                  do i=i1,i2
                     pe3(i,k) = pe0(i,k) 
                  enddo
               enddo
               do k=kmap+1,km
                  do i=i1,i2
                     pe3(i,k) =  (1.-bk(k))*pe0(i,k) + bk(k)*pe3(i,k)
                  enddo
               enddo
           else
               do k=1,ks+1
                  do i=i1,i2
                     pe3(i,k) = ak(k)
                  enddo
               enddo
           endif

           call map1_ppm( km,   pe0,   u, kmap,                  &
                km,   pe3,   u, ng_d, ng_s,            &
                im, i1, i2, j, jfirst, jlast, -1, kord_mt)
       endif

!------
! map v
!------
       if(j /= 1 .and. j /= jm) then
           do k=2,km+1

!-----------------------------------------------------------------------------
              if( i1 > 1 ) then
                  pe0(i1,k) = 0.5*(pe1(i1,k)+pe(i1-1,k,j))
              else
                  pe0(1,k) = 0.5*(pe1(1,k)+pe(im,k,j))
              endif
!-----------------------------------------------------------------------------

              do i=i1+1,i2
                 pe0(i ,k) = 0.5*(pe1(i,k)+pe1(i-1,k))
              enddo
           enddo

           do k=ks+2,km+1
!-----------------------------------------------------------------------------
              if ( i1 > 1 ) then
                  pe3(i1,k) = 0.5 * ( pe2(i1,k) +            &
                       (ak(k) + bk(k)*pe(i1-1,km+1,j)))
              else
                  pe3(1,k) = 0.5 * ( pe2(1,k) +            &
                       (ak(k) + bk(k)*pe(im,km+1,j)))
              endif
!-----------------------------------------------------------------------------

              do i=i1+1,i2
                 pe3(i,k) = 0.5*(pe2(i,k)+pe2(i-1,k))
              enddo
           enddo

           if ( p_map ) then
               do k=1,kmap
                  do i=i1,i2
                     pe3(i,k) = pe0(i,k) 
                  enddo
               enddo
               do k=kmap+1,km
                  do i=i1,i2
                     pe3(i,k) = (1.-bk(k))*pe0(i,k)+ bk(k)*pe3(i,k)
                  enddo
               enddo
           endif

           call map1_ppm ( km,   pe0,   v, kmap,                &
                km,   pe3,   v, ng_d, ng_d,          &
                im, i1, i2, j, jfirst, jlast, -1, kord_mt)
       endif

! Save new PE to pem
       do k=2,km+1
          do i=i1,i2
             pem(i,k,j) = pe2(i,k)
          enddo
       enddo

! Check deformation.
       if( diag .and. .not. p_map ) then
           rmax(ixj) = 0.
           rmin(ixj) = 1.
           do k=1,km
              do i=i1,i2
                 ratio(i) = (pe1(i,k+1)-pe1(i,k)) / (pe2(i,k+1)-pe2(i,k))
              enddo

              do i=i1,i2
                 if(ratio(i) > rmax(ixj)) then
                     rmax(ixj) = ratio(i)
                 elseif(ratio(i) < rmin(ixj)) then
                     rmin(ixj) = ratio(i)
                 endif
              enddo
           enddo
       endif
2000 end do !continue
    call fv_array_sync()
#ifdef SPMD
    call mp_send3d(gid-1, gid+1, im, jm, km, 1, im, jfirst-ng_d, jlast+ng_s, &
         1, km, 1, im, jfirst, jfirst, 1, km, u)
#endif

#ifdef VERT_FLUX
! use pkz for downward flux
    if(id_dnflux > 0) used = send_data ( id_dnflux,   pkz(beglon:endlon,beglat:endlat), Time_next )
    if(id_upflux > 0) used = send_data ( id_upflux, mflux(beglon:endlon,beglat:endlat), Time_next )
!      deallocate (mflux)
#endif

    if( diag .and. .not. p_map ) then
        qmin = rmin(1)
        do ixj=2, jp
           if(rmin(ixj) < qmin) then
               qmin = rmin(ixj)
           endif
        enddo

        qmax = rmax(1)
        do ixj=2, jp
           if(rmax(ixj) > qmax) then
               qmax = rmax(ixj)
           endif
        enddo
        if(qmin < 0.5 .or. qmax > 1.5) then
            write(6,*) 'PE=',gid, '[Rmin, Rmax] = [', qmin, ', ', qmax,']'
        endif
    endif

!$omp parallel do private(i,j,k)
    do j=jsp,jep !jfirst,jlast
       if ( p_map ) then
           do k=kmap+1,km
              do i=1,im
                 pe(i,k,j) = pem(i,k,j)
                 pk(i,j,k) = pe(i,k,j) ** akap
              enddo
           enddo
       else
           do k=2,km
              do i=1,im
                 pe(i,k,j) = pem(i,k,j)
              enddo
           enddo
       endif
    enddo
    call fv_array_sync()

    call pkez( nx, im, km, jfirst, jlast,               &
         pe, pk, akap, ks, peln, pkz, .not. p_map )

#ifdef SPMD
    call mp_recv3d(gid+1, im, jm, km, 1, im, jfirst-ng_d, jlast+ng_s, &
         1, km, 1, im, jlast+1, jlast+1, 1, km, u)
    u_ghosted = .true.
    call fv_array_sync()
#endif

!----- compute globally integrated TE ---------------
    if( consv > 1.E-5 ) then
!$omp parallel do private(i,j,k)
        do k=ksp,kep !1,km
           do j=jfirst,jlast
              do i=1,im
                 va(i,j,k) = ua(i,j,k) * delp(i,j,k)
              enddo
           enddo
        enddo
        call fv_array_sync()

!$omp parallel do private(i, j, k, bte, xsum)
        do j=jsp,jep !jfirst,jlast
! Perform vertical integration
           if ( j == 1 ) then
! SP
               tte(1) = 0.

               do k=1,km
                  tte(1) = tte(1) + va(1,1,k)
               enddo
               tte(1)  = acap * tte(1)

           elseif ( j == jm) then
! NP
               tte(jm) = 0.

               do k=1,km
                  tte(jm) = tte(jm) + va(1,jm,k)
               enddo
               tte(jm) = acap * tte(jm)

           else
! Interior
               do i=1, im
                  bte(i) = 0.
               enddo

               do k=1,km
                  do i=1,im
                     bte(i) = bte(i) + va(i,j,k)
                  enddo
               enddo

               xsum = 0.
               do i=1,im
                  xsum = xsum + bte(i)
               enddo
               tte(j) = xsum*cosp(j)

           endif

           tte(j) = tte0(j) - tte(j)
        enddo           ! end openMP j-loop
        call fv_array_sync()

        call par_vecsum(jm, jfirst, jlast, tte, te1)

!$omp parallel do private(i, j, xsum, ysum, zsum)
        do j=jsp,jep !jfirst,jlast

#ifdef FIX_TEMP
           if( j == 1 ) then
               tte(1) = acap*cp * (ps(1,1) - ptop -                  &
                    akap*ptop*(peln(1,km+1,1) - peln(1,1,1) ) )
           elseif( j == jm ) then
               tte(jm)= acap*cp * (ps(1,jm) - ptop -                 &
                    akap*ptop*(peln(1,km+1,jm) - peln(1,1,jm) ) )
           else
               xsum = 0.
               ysum = 0.
               do i=1,im
                  xsum = xsum + ps(i,j)
                  ysum = ysum + peln(i,km+1,j)
               enddo
               tte(j) = cp*cosp(j)*(xsum - ptop*im -                   &
                    akap*ptop*(ysum - peln(1,1,j)*im) )
           endif
#else
           if( j==1 .or. j==jm ) then
               ysum = ptop*(pk(1,j,1)-pk(1,j,km+1))
               do k=1,km
                  ysum = ysum + pkz(1,j,k)*delp(1,j,k)
               enddo

               tte(j) = cp*acap*ysum
           else

               do i=1,im
                  zsum(i) = ptop*(pk(i,j,1)-pk(i,j,km+1))
               enddo

               do k=1,km
                  do i=1,im
                     zsum(i) = zsum(i) + pkz(i,j,k)*delp(i,j,k)
                  enddo
               enddo

               xsum = 0.
               do i=1,im
                  xsum = xsum + zsum(i)
               enddo

               tte(j) = cp*cosp(j)*xsum
           endif
#endif
        enddo
        call fv_array_sync()

        call par_vecsum(jm, jfirst, jlast, tte, tsum)

        dtmp = min(1., consv) * te1 / tsum        ! do not allow over-correction
!       if(master) write(6,*) 'Energy correction in T/PT (deg/yr) =',dtmp*86400*365/mdt
    else
        dtmp = 0.
    endif        ! end consv check


!$omp parallel do private(i, j, k, u2, v2, te_sp, te_np)
    do 8000 k=max(kmap,ksp),kep !kmap,km
! Compute KE
       do j=js2g0, min(jlast+1,jm)
          do i=1,im
             u2(i,j) = cose(j) * u(i,j,k)**2
          enddo
       enddo

       do j=js2g0, jn2g0
          do i=1,im
             v2(i,j) = v(i,j,k)**2
          enddo
       enddo

       do j=js2g0, jn2g0
          do i=1,im-1
             ua(i,j,k) = ua(i,j,k) - 0.25 * ( (u2(i,j)+u2(i,j+1))*acosu(j)  &
                  +v2(i,j) + v2(i+1,j) )
          enddo
          ua(im,j,k) = ua(im,j,k) - 0.25*( (u2(im,j)+u2(im,j+1))*acosu(j) &
               +v2(im,j) + v2(1,j) )
       enddo

! poles
       if ( jfirst == 1 ) then
           te_sp = 0.
           do i=1,im
!           te_sp = te_sp + u2(i,2) + v2(i,2)
! SJL: 01/26/04
              te_sp = te_sp + 0.5*u2(i,2)/cose(2)
           enddo
           te_sp = ua(1,1,k) - 0.5*te_sp/float(im)

           do i=1,im
              ua(i,1,k) = te_sp
           enddo
       endif

       if ( jlast == jm ) then
           te_np = 0.
           do i=1,im
!           te_np = te_np + u2(i,jm) + v2(i,jm-1)
! SJL: 01/26/04
              te_np = te_np + 0.5*u2(i,jm)/cose(jm)
           enddo

           te_np = ua(1,jm,k) - 0.5*te_np/float(im)
           do i=1,im
              ua(i,jm,k) = te_np
           enddo
       endif
8000 end do !continue
    call fv_array_sync()

! Recover temperature
!$omp parallel do private(ixj, i1, i2, i, j, k, rg, gz, tm, dlnp)

    do ixj=ixjs,ixje !1,jp

       j  = jfirst + (ixj-1) / nx
       i1 = 1 + it * mod(ixj-1, nx)
       i2 = i1 + it - 1

       rg = akap * cp
       do i=i1,i2
          gz(i) = hs(i,j)      
       enddo

       if ( p_map ) then
           do k=1,kmap-1
              do i=i1,i2
! Recover temperature from Cp * virt_potential temp
                 pt(i,j,k) = pt(i,j,k)*pkz(i,j,k) / (cp*(1.+r_vir*q(i,j,k,1)))
              enddo
           enddo
       endif

       do k=km,kmap,-1
          do i=i1,i2
             dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
!             tm = (ua(i,j,k)-gz(i)) / ( cp*(1.+cp_vir*q(i,j,k,1)) -      &
!                        (1.+r_vir*q(i,j,k,1))*pe(i,k,j)*dlnp/delp(i,j,k) )
             tm = (ua(i,j,k)-gz(i)) / ( (cp -        &
                  pe(i,k,j)*dlnp/delp(i,j,k))*(1.+r_vir*q(i,j,k,1)) )
#ifdef FIX_TEMP
             pt(i,j,k) = tm + dtmp        ! pt is now temperature
#else
             pt(i,j,k) = tm + dtmp*pkz(i,j,k)/(1.+r_vir*q(i,j,k,1))
#endif
             gz(i) = gz(i) + dlnp*tm*(1.+r_vir*q(i,j,k,1))
          enddo
       enddo           ! end k-loop

    enddo
    call fv_array_sync()

!----------------------
! Output A grid winds:
!----------------------
    call d2a3d(u, v, ua, va, im, jm, km, jfirst, jlast,    &
         ng_d, ng_s, coslon, sinlon)


    call fv_print_chksums( 'Exiting  te_map' )
  end subroutine te_map

  subroutine geo_map(ps, omga, pe, delp, pkz, pk, mdt,            &
       im, jm, km, nx, jfirst, jlast, nq,  u,  v,   &
       pt, q, r_vir, cp_vir,  cp, akap,             &
       kord_mt, kord_tr, kord_tm,  peln,            &
       phis, ng_d, ng_s, rg, acap, cosp, ua, va,    &
       consv, tte0, fill)
! Alternative remapping algorithm based on preservation of
! geopotential; total energy is not conserved. This algorithm has the
! potential of being more accurate by NOT perturbing the pressure gradient
! with remapping
!
#ifdef use_shared_pointers
    use fv_arrays_mod, only: is,ie,js,je, isd,ied,jsd,jed, isg,ieg,jsg,jeg
    use fv_arrays_mod, only: nlev
    use fv_arrays_mod, only: ptr_ps_bp, ptr_pesouth
    real :: ps_bp (isg:ieg, js:je)
    pointer( p_ps_bp, ps_bp )
    real :: pesouth(isg:ieg, nlev+1)
    pointer( p_pesouth, pesouth )
#else
    use fv_arrays_mod, only: ps_bp, pesouth
#endif

! !INPUT PARAMETERS:
    integer, intent(in):: mdt        ! mapping time step (same as phys)
    integer, intent(in):: im         ! E-W dimension
    integer, intent(in):: jm         ! N-S dimension
    integer, intent(in):: km         ! Vertical dimension
    integer, intent(in):: nq         ! number of tracers (including h2o)
    integer, intent(in):: ng_d
    integer, intent(in):: ng_s
    integer, intent(in):: nx            ! number of SMP "decomposition" in x
    integer, intent(in):: jfirst, jlast ! starting & ending latitude index
    integer, intent(in):: kord_mt, kord_tr, kord_tm

    real, intent(in)::  r_vir
    real, intent(in):: cp_vir
    real, intent(in):: cp, rg, acap
    real, intent(in):: akap
    real, intent(in):: phis(im,jfirst:jlast)
    real, intent(in):: cosp(jm)
    real, intent(in):: consv
    logical,  intent(in):: fill

! !INPUT/OUTPUT PARAMETERS:
    real, intent(inout):: tte0(jfirst:jlast)
    real, intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)   ! u-wind (m/s)
    real, intent(inout):: v(im,jfirst-ng_d:jlast+ng_d,km)   ! v-wind (m/s)
    real, intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d,km)  ! virtual potential temperature
    real, intent(inout):: delp(im,jfirst:jlast,km) ! pressure thickness
    real, intent(inout):: q(im,jfirst-ng_d:jlast+ng_d,km,*)
    real, intent(inout):: pe(im,km+1,jfirst:jlast) ! pressure at layer edges
    real, intent(inout):: pk(im,jfirst:jlast,km+1) ! pe to the kappa
    real, intent(inout):: ps(im,jfirst:jlast)      ! surface pressure

! !OUTPUT PARAMETERS:
    real, intent(out):: omga(im,jfirst:jlast,km)    ! vertical press. velocity (pascal/sec)
    real, intent(out):: peln(im,km+1,jfirst:jlast)  ! log(pe)
    real, intent(out):: pkz(im,jfirst:jlast,km)     ! layer-mean pk for converting t to pt

! Work
    real, intent(out):: ua(im,jfirst:jlast,km)
    real, intent(out):: va(im,jfirst:jlast,km)
! !DESCRIPTION:
! This routine is based on te_map. However, total energy is no longer conserved.
! Instead, the geopotential is conserved.
!
! !REVISION HISTORY:
!
! SJL 02.04.07: Geopotential conserving mapping initial version based on te_map
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Local arrays:
    real tte(jfirst:jlast)
    real rmin(nx*jm), rmax(nx*jm)
! y-z
    real  q2(im/nx,km,nq)
    real dp2(im/nx,km)
    real pe0(im,km+1)
    real pe1(im,km+1)
    real pe2(im,km+1)
    real pe3(im,km+1)
    real ratio(im)
    real pe1w(km+1)
    real pe2w(km+1)

    real dak, bkh, qmax, qmin
    real rdt5
    real pktmp, tmp, ak1
    real xsum, ysum, tsum, dtmp, te1

    integer i, j, k, iq
    integer ixj, jp, it, i1, i2

    integer :: ixjs, ixje
    logical diag
    data diag    /.false./
#ifdef use_shared_pointers
    pointer( p_tte, tte )
    call fv_stack_push( p_tte, jlast-jfirst+1 )
#endif

!     diag = .true.
#ifdef use_shared_pointers
    p_ps_bp = ptr_ps_bp
    p_pesouth = ptr_pesouth
#endif

    call fv_print_chksums( 'Entering  geo_map' )
    call fv_array_check( LOC(phis) )
    call fv_array_check( LOC(u) )
    call fv_array_check( LOC(v) )
    call fv_array_check( LOC(pt) )
    call fv_array_check( LOC(delp) )
    call fv_array_check( LOC(q) )
    call fv_array_check( LOC(pe) )
    call fv_array_check( LOC(pk) )
    call fv_array_check( LOC(ps) )
    call fv_array_check( LOC(omga) )
    call fv_array_check( LOC(peln) )
    call fv_array_check( LOC(pkz) )
    call fv_array_check( LOC(ua) )
    call fv_array_check( LOC(va) )

#ifdef SPMD
    call mp_send3d(gid+1, gid-1, im, km+1, jm, 1, im, 1, km+1, jfirst, jlast,&
         1, im, 1, km+1, jlast, jlast, pe)
    call mp_recv3d(gid-1, im, km+1, jm, 1, im, 1, km+1, jfirst-1, jfirst-1, &
         1, im, 1, km+1, jfirst-1, jfirst-1, pesouth)
    u_ghosted = .false.
    call fv_array_sync()
#endif

    it = im / nx
    jp = nx * ( jlast - jfirst + 1 )

    call fv_array_limits( 1, jp, ixjs, ixje )
!$omp parallel do default(shared)                 &
!$omp private(i,j,k,pe0,pe1,pe2,pe3,ratio)        &
!$omp private(dak, bkh, rdt5, ixj, i1, i2, iq)   &
!$omp private(pe1w, pe2w, pktmp, q2, dp2)

    do 2000 ixj=ixjs, ixje !1, jp

       j  = jfirst + (ixj-1) / nx
       i1 = 1 + it * mod(ixj-1, nx)
       i2 = i1 + it - 1

! Copy data to local 2D arrays.
       do k=1,km+1
          do i=i1,i2
             pe1(i,k) = pe(i,k,j)
          enddo
       enddo

       do k=1,ks+1
          do i=i1,i2
             pe0(i,k) = ak(k)
             pe2(i,k) = ak(k)
             pe3(i,k) = ak(k)
          enddo
       enddo

       do k=ks+2,km
          do i=i1,i2
             pe0(i,k) = ak(k) + bk(k)* ps_bp(i,j)
             pe2(i,k) = ak(k) + bk(k)*pe1(i,km+1)
          enddo
       enddo

       do i=i1,i2
          pe0(i,km+1) =  ps_bp(i,j)
          pe2(i,km+1) = pe1(i,km+1)
       enddo

! Compute omga (dp/dt)
       do k=2,km+1
          do i=i1,i2
             pe0(i,k) = pe1(i,k) - pe0(i,k)
          enddo
       enddo

       rdt5 = 0.5 / float(mdt)
       do i=i1,i2
! update ps
          ps(i,j) = pe1(i,km+1)
          ps_bp(i,j) = ps(i,j)
          omga(i,j,1) = rdt5 * pe0(i,2)
       enddo

       do k=2,km
          do i=i1,i2
             omga(i,j,k) = rdt5 * ( pe0(i,k) + pe0(i,k+1) )
          enddo
       enddo

       if(ks /= 0) then
           do k=1,ks
              dak = ak(k+1) - ak(k)
              do i=i1,i2
                 delp(i,j,k) = dak
              enddo
           enddo
       endif

       do k=ks+1,km
          do i=i1,i2
             delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
          enddo
       enddo

! Check deformation.
       if( diag ) then
           rmax(ixj) = 0.
           rmin(ixj) = 1.
           do k=1,km
              do i=i1,i2
                 ratio(i) = (pe1(i,k+1)-pe1(i,k)) / (pe2(i,k+1)-pe2(i,k))
              enddo

              do i=i1,i2
                 if(ratio(i) > rmax(ixj)) then
                     rmax(ixj) = ratio(i)
                 elseif(ratio(i) .lt. rmin(ixj)) then
                     rmin(ixj) = ratio(i)
                 endif
              enddo
           enddo
       endif

!-----------------
! Map constituents
!-----------------

       if(nq /= 0) then

           do k=1,km
              do i=i1,i2
                 dp2(i-i1+1,k) = pe2(i,k+1) - pe2(i,k)
              enddo
           enddo


           call mapn_ppm ( km,   pe1,   q, nq,   1,                &
                km,   pe2,  q2, ng_d, ng_d,             &
                im,   i1, i2, j, jfirst, jlast, 0, kord_tr)

           if (fill) call fillz(i2-i1+1, km, nq, q2, dp2)
           do iq=1,nq
              do k=1,km
                 do i=i1,i2
                    q(i,j,k,iq) = q2(i-i1+1,k,iq)
                 enddo
              enddo
           enddo
       endif

!-------
! map u
!-------
       if(j /= 1) then

! WS 99.07.29 : protect j==jfirst case
           if (j > jfirst) then
               do k=2,km+1
                  do i=i1,i2
                     pe0(i,k) = 0.5*(pe1(i,k)+pe(i,k,j-1))
                  enddo
               enddo

               do k=ks+2,km+1
                  bkh = 0.5*bk(k)
                  do i=i1,i2
                     pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pe(i,km+1,j-1))
                  enddo
               enddo

#ifdef SPMD
           else
               do k=2,km+1
                  do i=i1,i2
                     pe0(i,k) = 0.5*(pe1(i,k)+pesouth(i,k))
                  enddo
               enddo

               do k=ks+2,km+1
                  bkh = 0.5*bk(k)
                  do i=i1,i2
                     pe3(i,k) = ak(k) + bkh*(pe1(i,km+1)+pesouth(i,km+1))
                  enddo
               enddo
#endif
           endif

           call map1_ppm( km,   pe0,   u, 1,                      &
                km,   pe3,   u, ng_d, ng_s,             &
                im, i1, i2, j, jfirst, jlast, -1, kord_mt)
       endif

! map v
       if(j .ne. 1 .and. j .ne. jm) then

           do k=1,km+1
              if( i1 == 1 ) then
                  pe1w(k) = pe(im,k,j)
              else
                  pe1w(k) = pe(i1-1,k,j)
              endif
           enddo

           do k=ks+2,km
              pe2w(k) = ak(k) + bk(k)*pe1w(km+1)
           enddo
           pe2w(km+1) = pe1w(km+1)

           do k=2,km+1
              pe0(i1,k) = 0.5*(pe1(i1,k)+pe1w(k))
              do i=i1+1,i2
                 pe0(i ,k) = 0.5*(pe1(i,k)+pe1(i-1,k))
              enddo
           enddo

           do k=ks+2,km+1
              pe3(i1,k) = 0.5*(pe2(i1,k)+pe2w(k))
              do i=i1+1,i2
                 pe3(i,k) = 0.5*(pe2(i,k)+pe2(i-1,k))
              enddo
           enddo

           call map1_ppm ( km,   pe0,   v, 1,                    &
                km,   pe3,   v, ng_d, ng_d,           &
                im, i1, i2, j, jfirst, jlast, -1, kord_mt)
       endif

! Save new PE to peln
       do k=2,km
          do i=i1,i2
             peln(i,k,j) = pe2(i,k)
          enddo
       enddo

!-------
! Map pt
!-------
       do k=1,km+1
          do i=i1,i2
             pe1(i,k) = pk(i,j,k)
          enddo
       enddo

       do k=1,ks+1
          pktmp = ak(k) ** akap
          do i=i1,i2
             pe2(i,k)  = pktmp
             pk(i,j,k) = pktmp
          enddo
       enddo

       do k=ks+2,km+1
          do i=i1,i2
             pe2(i,k)  = pe2(i,k) ** akap
             pk(i,j,k) = pe2(i,k)
          enddo
       enddo

       call map1_ppm ( km,   pe1,   pt,  1,                   &
            km,   pe2,   pt, ng_d, ng_d,           &
            im, i1, i2, j, jfirst, jlast, 1, kord_tr)
2000 end do !continue
    call fv_array_sync()
    if( diag ) then
        qmin = rmin(1)
        do ixj=2, jp
           if(rmin(ixj) < qmin) then
               qmin = rmin(ixj)
           endif
        enddo

        qmax = rmax(1)
        do ixj=2, jp
           if(rmax(ixj) > qmax) then
               qmax = rmax(ixj)
           endif
        enddo
        if(master) write(6,*) 'Rmin=', qmin, 'Rmax=', qmax
    endif

    ak1 = (akap + 1.) / akap

!$omp parallel do private(ixj, i1, i2, i, j, k, tmp)
    do 4000 ixj=ixjs,ixje !1,jp

       j  = jfirst + (ixj-1) / nx
       i1 = 1 + it * mod(ixj-1, nx)
       i2 = i1 + it - 1

       do k=2,km
          do i=i1,i2
             pe(i,k,j) = peln(i,k,j)
          enddo
       enddo

! Compute updated peln

       do k=2,ks+1
          tmp = log(ak(k))
          do i=i1,i2
             peln(i,k,j) = tmp
          enddo
       enddo

       do k=ks+2,km+1
          do i=i1,i2
             peln(i,k,j) = log(pe(i,k,j))
          enddo
       enddo

       if( ptop < ptop_min ) then
           do i=i1,i2
              peln(i,1,j) = peln(i,2,j) - ak1
           enddo
       else
           tmp = log( ptop )
           do i=i1,i2
              peln(i,1,j) = tmp
           enddo
       endif

! Compute pkz
       do k=1,km
          do i=i1,i2
             pkz(i,j,k) = ( pk(i,j,k+1) - pk(i,j,k) )  /         &
                  (akap*(peln(i,k+1,j) - peln(i,k,j)) )
          enddo
       enddo

! Recover temperature
       do k=1,km
          do i=i1,i2
             pt(i,j,k) = pt(i,j,k) * pkz(i,j,k) / ( cp*(1.+r_vir*q(i,j,k,1)) )
          enddo
       enddo
4000 end do !continue
    call fv_array_sync()

    if ( consv > 1.E-5 ) then
        call p_energy(im, jm, km, u, v, pt, delp, q(1,jfirst-ng_d,1,1),  &
             pe, peln, phis, ng_d, ng_s,                         &
             r_vir, cp_vir,  cp, rg,  tte,                    &
             jfirst, jlast, acap, cosp, ua, va)

!$omp parallel do private(i, j, xsum, ysum)
        do j=jsp,jep !jfirst,jlast

           tte0(j) = tte0(j) - tte(j)

           if( j == 1 ) then
               tte(1) = acap*cp * (ps(1,1) - ptop -                  &
                    akap*ptop*(peln(1,km+1,1) - peln(1,1,1) ) )
           elseif( j == jm ) then
               tte(jm)= acap*cp * (ps(1,jm) - ptop -                 &
                    akap*ptop*(peln(1,km+1,jm) - peln(1,1,jm) ) )
           else
               xsum = 0.
               ysum = 0.
               do i=1,im
                  xsum = xsum + ps(i,j)
                  ysum = ysum + peln(i,km+1,j)
               enddo
               tte(j) = cp*cosp(j)*(xsum - ptop*im -                   &
                    akap*ptop*(ysum - peln(1,1,j)*im) )
           endif
        enddo
        call fv_array_sync()

        call par_vecsum(jm, jfirst, jlast, tte0, te1)
        call par_vecsum(jm, jfirst, jlast, tte,  tsum)

        dtmp = min(1., consv) * te1 / tsum

!     if(master) write(6,*) 'dtmp=', dtmp*86400*365/mdt

!$omp parallel do private(i, j, k)
        do k=ksp,kep !1,km
           do j=jfirst,jlast
              do i=1,im
                 pt(i,j,k) = pt(i,j,k) + dtmp
              enddo
           enddo
        enddo
        call fv_array_sync()

    else
!----------------------
! Output A grid winds:
!----------------------
        call d2a3d(u, v, ua,  va, im, jm, km, jfirst, jlast,      &
             ng_d, ng_s, coslon, sinlon)

    endif

    call fv_print_chksums( 'Exiting  geo_map' )
  end subroutine geo_map


  subroutine p_energy(im, jm, km, u, v, pt, delp, q, pe, peln, phis,   &
       ng_d, ng_s, r_vir, cp_vir,  cp, rg, tte,       &
       jfirst, jlast, acap, cosp, ua, va )
! !INPUT PARAMETERS:
    integer,  intent(in):: im, jm, km, jfirst, jlast          ! Dimensions
    integer,  intent(in):: ng_d, ng_s
    real, intent(in):: r_vir
    real, intent(in):: cp_vir
    real, intent(in):: cp
    real, intent(in):: rg
    real, intent(in):: acap

    real, intent(in):: cosp(jm)
    real, intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-wind
    real, intent(in):: v(im,jfirst-ng_d:jlast+ng_d,km)  ! Winds y
    real, intent(in):: pt(im,jfirst-ng_d:jlast+ng_d,km)  ! temperature
    real, intent(in):: delp(im,jfirst:jlast,km)          ! Delta pressure
    real, intent(in):: q(im,jfirst-ng_d:jlast+ng_d,km)   ! specific humidity

    real, intent(in):: pe(im, km+1, jfirst:jlast)        ! Edge pressure
    real, intent(in):: peln(im, km+1, jfirst:jlast)      ! log(pe)
    real, intent(in):: phis(im,jfirst:jlast)

    real, intent(out):: ua(im, jfirst:jlast, km)
    real, intent(out):: va(im, jfirst:jlast, km)
    real, intent(out):: tte(jfirst:jlast)      ! column integrated Total Energy

! Local
    real bte(im)
    real gztop(im)

    real te_sp
    real te_np
    real xsum

    integer i, j, k, js2g0, jn2g0


    call fv_print_chksums( 'Entering  p_energy' )
    call fv_array_check( LOC(u) )
    call fv_array_check( LOC(v) )
    call fv_array_check( LOC(pt) )
    call fv_array_check( LOC(delp) )
    call fv_array_check( LOC(q) )
    call fv_array_check( LOC(pe) )
    call fv_array_check( LOC(peln) )
    call fv_array_check( LOC(phis) )
    call fv_array_check( LOC(ua) )
    call fv_array_check( LOC(va) )
    call fv_array_check( LOC(tte) )
    js2g0  = max(2,jfirst)
    jn2g0  = min(jm-1,jlast)

!----------------------
! Output A grid winds:
!----------------------
    call d2a3d(u, v, ua,  va, im, jm, km, jfirst, jlast,      &
         ng_d, ng_s, coslon, sinlon)

!$omp parallel do private(i,j,k,bte, xsum, gztop)
    do j=jsp,jep !jfirst,jlast
! Perform vertical integration
       do i=1,im
          gztop(i) = phis(i,j)
          do k=1,km
             gztop(i) = gztop(i) + (peln(i,k+1,j)-peln(i,k,j)) *   &
                  rg*pt(i,j,k)*(1.+r_vir*q(i,j,k))
          enddo
       enddo

       if (j == 1) then
! SP
           tte(1) = pe(1,km+1,1)*phis(1,1) - pe(1,1,1)*gztop(1)
           do k=1,km
              tte(1) = tte(1) + delp(1,1,k) * (0.5*(ua(1,1,k)**2+va(1,1,k)**2) + &
                   cp*pt(1,1,k)*(1.+r_vir*q(1,1,k)))
           enddo
           tte(1)  = acap * tte(1)
       elseif (j == jm) then
! NP
           tte(jm) = pe(1,km+1,jm)*phis(1,jm) - pe(1,1,jm)*gztop(1)
           do k=1,km
              tte(jm) = tte(jm) + delp(1,jm,k) * (0.5*(ua(1,jm,k)**2+va(1,jm,k)**2) + &
                   cp*pt(1,jm,k)*(1.+r_vir*q(1,jm,k)))
           enddo
           tte(jm) = acap * tte(jm)
       else
           do i=1,im
              bte(i) = pe(i,km+1,j)*phis(i,j) - pe(i,1,j)*gztop(i)
           enddo

           do k=1,km
              do i=1,im
                 bte(i) = bte(i) + delp(i,j,k) * (0.5*(ua(i,j,k)**2+va(i,j,k)**2)  + &
                      cp*pt(i,j,k)*(1.+r_vir*q(i,j,k)))
              enddo
           enddo

           xsum = 0.
           do i=1,im
              xsum = xsum + bte(i)
           enddo

           tte(j) = xsum*cosp(j)

       endif
    enddo
    call fv_array_sync()

    call fv_print_chksums( 'Exiting  p_energy' )
  end subroutine p_energy

  subroutine map1_q2( km,   pe1,   q1,   kmap,         &
                      kn,   pe2,   q2,   itot,         &
                      i1,   i2,    iv,   kord, j, jfirst, jlast, ng)


! !INPUT PARAMETERS:
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng                ! Ghosted latitudes south
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kmap              ! starting remapping k-index
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(itot,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real, intent(in) ::  pe2(itot,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real, intent(in) ::  q1(itot,jfirst-ng:jlast+ng,km) ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i2-i1+1,kn) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    incorporated latest FVGCM version
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real, parameter:: r3 = 1./3., r23 = 2./3.
      real   dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, delp, esl

      integer i, k, l, ll, k0

      do k=max(1,kmap-3),km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, kmap, i1, i2, iv, kord )

! Mapping
      do 1000 i=i1,i2
         k0 = kmap
      do 555 k=kmap,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i-i1+1,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                 delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                 qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i-i1+1,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

!EOC
 end subroutine map1_q2

 subroutine map1_ppm( km,   pe1,    q1,   kmap,         &
                      kn,   pe2,    q2,                 &
                      ng_s, ng_n,   itot,  i1,  i2,     &
                      j,    jfirst, jlast, iv,  kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kmap              ! starting remapping k-index
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(itot,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real, intent(in) ::  pe2(itot,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real, intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km) ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    incorporated latest FVGCM version
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real, parameter:: r3 = 1./3., r23 = 2./3.
      real    dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real    pl, pr, qsum, delp, esl

      integer i, k, l, ll, k0

      do k=max(1,kmap-3),km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, kmap, i1, i2, iv, kord )

! Mapping
      do 1000 i=i1,i2
         k0 = kmap
      do 555 k=kmap,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                 delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                 qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

!EOC
 end subroutine map1_ppm
!--------------------------------------------------------- 

 subroutine mapn_ppm(km,     pe1,   q1,  nq,   kmap,     &
                     kn,     pe2,   q2,  ng_s, ng_n,     &
                     itot,   i1,    i2,  j,              &
                     jfirst, jlast, iv, kord)

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == other scalars
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kmap              ! partial remapping
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: nq                ! Number of tracers

      real , intent(in) :: pe1(itot,km+1)   ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real , intent(in) :: pe2(itot,kn+1)   ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real , intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km,nq) ! Field input
! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout)::  q2(i2-i1+1,kn,nq) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

   real , parameter:: r3= 1./3., r23 = 2./3.
   real  dp1(i1:i2,km)
   real   q4(4,i1:i2,km,nq)
   real   qsum(nq)
   real   pl, pr, delp, esl
   real   tmp1, tmp2
   integer i, k, l, m, mm, k0, iq

   do k=max(1,kmap-3),km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
      enddo
   enddo

   do iq=1,nq
      do k=max(1,kmap-3),km
         do i=i1,i2
            q4(1,i,k,iq) = q1(i,j,k,iq)
         enddo
      enddo
! Compute vertical subgrid distribution
      call ppm2m( q4(1,i1,1,iq), dp1, km, kmap, i1, i2, iv, kord )
   enddo

! Mapping
      do 1000 i=i1,i2
         k0 = kmap
      do 555 k=kmap,kn
      do 100 m=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,m) .and. pe2(i,k) <= pe1(i,m+1)) then
         pl = (pe2(i,k)-pe1(i,m)) / dp1(i,m)
         if(pe2(i,k+1) <= pe1(i,m+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,m)) / dp1(i,m)
            tmp2 = pr + pl
            tmp1 = r3 * ( pr*tmp2 + pl**2 ) 
            do iq=1,nq
               q2(i-i1+1,k,iq) = q4(2,i,m,iq) + 0.5*(q4(4,i,m,iq)+q4(3,i,m,iq) -    &
                              q4(2,i,m,iq))*tmp2 - q4(4,i,m,iq)*tmp1
            enddo
            k0 = m
            goto 555
         else
! Fractional area...
            tmp1 = r3 * ( 1. + pl*(1.+pl) )
            do iq=1,nq
               qsum(iq) = (pe1(i,m+1)-pe2(i,k))*(q4(2,i,m,iq)+0.5*(q4(4,i,m,iq)+     &
                           q4(3,i,m,iq)-q4(2,i,m,iq))*(1.+pl)-q4(4,i,m,iq)*tmp1 )
            enddo
            do mm=m+1,km
! locate the bottom edge: pe2(i,k+1)
               if(pe2(i,k+1) > pe1(i,mm+1) ) then
! Whole layer..
                 do iq=1,nq
                    qsum(iq) = qsum(iq) + dp1(i,mm)*q4(1,i,mm,iq)
                 enddo
               else
                 delp = pe2(i,k+1) - pe1(i,mm)
                  esl = delp / dp1(i,mm)
                 tmp1 = 1. - r23*esl
                 do iq=1,nq
                    qsum(iq) = qsum(iq) + delp*(q4(2,i,mm,iq)+0.5*esl*            &
                       (q4(3,i,mm,iq)-q4(2,i,mm,iq)+q4(4,i,mm,iq)*tmp1 ))
                 enddo
                 k0 = mm
                 goto 123
               endif
            enddo
            goto 123
         endif
      endif
100   continue
123   do iq=1,nq 
         q2(i-i1+1,k,iq) = qsum(iq) / ( pe2(i,k+1) - pe2(i,k))
      enddo
555   continue
1000  continue

 end subroutine mapn_ppm

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m(a4, delp, km, kmap, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kmap    ! partial remap to start
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real , intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real , intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real    dc(i1:i2,km)
      real    h2(i1:i2,km)
      real  delq(i1:i2,km)
      real   df2(i1:i2,km)
      real    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt
      integer it
      real  fac
      real  a1, a2, c1, c2, c3, d1, d2
      real  qmax, qmin, cmax, cmin
      real  qm, dq, tmp
      real  qmp, pmp
      real  lac

      km1 = km - 1
       it = i2 - i1 + 1

      do k=max(2,kmap-2),km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo
 
      do k=max(2,kmap-2),km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=max(3,kmap), km1
      do i=i1,i2
        c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                  ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      if(km>8 .and. kord>3) call steepz(i1, i2, km, kmap, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      if ( kmap <= 2 ) then
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
         dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
         cmax = max(a4(1,i,1), a4(1,i,2))
         cmin = min(a4(1,i,1), a4(1,i,2))
         a4(2,i,2) = max(cmin,a4(2,i,2))
         a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

      do k=max(1,kmap),km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo

! Enforce monotonicity of the "slope" within the top layer
      if ( kmap <= 2 ) then
      do i=i1,i2
         if ( a4(2,i,1) * a4(1,i,1) <= 0. ) then 
              a4(2,i,1) = 0.
                dc(i,1) = a4(1,i,1)
         endif
         if ( dc(i,1) * (a4(2,i,2) - a4(1,i,1)) <= 0. ) then
! Setting DC==0 will force piecewise constant distribution after
! calling kmppm
              dc(i,1) = 0.
         endif
      enddo
      endif

! Enforce constraint on the "slope" at the surface

      do i=i1,i2
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) then
               dc(i,km) = 0.
         endif
         if( dc(i,km) * (a4(1,i,km) - a4(2,i,km)) <= 0. ) then
             dc(i,km) = 0.
         endif
      enddo
 
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      if ( kmap <= 2 ) then
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
            call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
      endif

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=max(2,kmap-1), km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=max(3,kmap), km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=max(3,kmap), km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
!EOC
 end subroutine ppm2m
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm(dm, a4, itot, lmt)

! !INPUT PARAMETERS:
      real , intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real , parameter:: r12 = 1./12.
      real  qmp
      real  da1, da2, a6da
      real  fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2003)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

!EOC
 end subroutine kmppm
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz(i1, i2, km, kmap, a4, df2, dm, dq, dp, d4)

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: kmap                 ! 
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real , intent(in) ::  dp(i1:i2,km)       ! grid size
      real , intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real , intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real , intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real , intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real , intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real  alfa(i1:i2,km)
      real     f(i1:i2,km)
      real   rat(i1:i2,km)
      real   dg2

! Compute ratio of dq/dp
      do k=max(2,kmap-1),km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=max(2,kmap-1),km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=max(3,kmap),km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=max(4,kmap+1),km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

!EOC
 end subroutine steepz


! !ROUTINE: benergy --- Calculate the total energy Before the Lagrangian dynamics
 
 subroutine benergy(im, jm, km, u, v, pt, delp, q, pe, peln, phis,    &
      u_ghosted,  ng_d, ng_s, r_vir, cp_vir, cp, rg,    &
      tte, jfirst, jlast, acap, cosp, te, dz )
! !USES:


! !INPUT PARAMETERS:
   integer,  intent(in):: im, jm, km, jfirst, jlast          ! Dimensions
   integer,  intent(in):: ng_d, ng_s
   real, intent(in):: r_vir
   real, intent(in):: cp_vir
   real, intent(in):: cp
   real, intent(in):: rg
   real, intent(in):: acap

   real, intent(in):: cosp(jm)
   real, intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-wind be ghosted one lat N
   real, intent(in)::  v(im,jfirst-ng_d:jlast+ng_d,km)  ! Winds y
   real, intent(in):: pt(im,jfirst-ng_d:jlast+ng_d,km)  ! temperature
   real, intent(in):: delp(im,jfirst:jlast,km)          ! Delta pressure
   real, intent(in):: q(im,jfirst-ng_d:jlast+ng_d,km)   ! specific humidity

   real, intent(in):: pe(im, km+1, jfirst:jlast)        ! Edge pressure
   real, intent(in):: peln(im, km+1, jfirst:jlast)      ! log(pe)
   real, intent(in):: phis(im,jfirst:jlast)
   real, intent(out):: tte(jfirst:jlast)                ! column integrated Total Energy

   logical, intent(inout):: u_ghosted


! !REVISION HISTORY:
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
   real u2(im,jfirst:jlast+1)
   real v2(im,jfirst:jlast)
   real tm(im,jfirst:jlast)

   real bte(im)
   real gztop(im)

   real te_sp
   real te_np
   real xsum

   integer i, j, k, js2g0, jn2g0
!Balaji: te/dz are listed as local work arrays but are actually
! reusing ua/va which are shared arrays!
! work arrays:
   real te(im, jfirst:jlast, km)
   real dz(im, jfirst:jlast, km)


   call fv_print_chksums( 'Entering  benergy' )
   call fv_array_check( LOC(u) )
   call fv_array_check( LOC(v) )
   call fv_array_check( LOC(pt) )
   call fv_array_check( LOC(delp) )
   call fv_array_check( LOC(q) )
   call fv_array_check( LOC(pe) )
   call fv_array_check( LOC(peln) )
   call fv_array_check( LOC(phis) )
   call fv_array_check( LOC(te) )
   call fv_array_check( LOC(dz) )
   call fv_array_check( LOC(tte) )

   js2g0  = max(2,jfirst)
   jn2g0  = min(jm-1,jlast)

!     if ( .not. u_ghosted ) then
!!        call error_mesg('Benergy:', 'U not ghosted',FATAL)
!     call mp_send3d(gid-1, gid+1, im, jm, km, 1, im, jfirst-ng_d, jlast+ng_s, &
!                        1, km, 1, im, jfirst, jfirst, 1, km, u)
!     call mp_recv3d(gid+1, im, jm, km, 1, im, jfirst-ng_d, jlast+ng_s, &
!                        1, km, 1, im, jlast+1, jlast+1, 1, km, u)
!           u_ghosted = .true.
!     endif

!$omp parallel do private(i, j, k, u2, v2, tm, te_sp, te_np)
   do k=ksp,kep !1,km
      do j=js2g0,min(jlast+1,jm)
         do i=1,im
            u2(i,j) = cose(j) * u(i,j,k)**2
         enddo
      enddo

      do j=js2g0,jn2g0
         do i=1,im
            v2(i,j) = v(i,j,k)**2
         enddo
      enddo

      do j=js2g0,jn2g0
         do i=1,im-1
            te(i,j,k) = 0.25*((u2(i,j)+u2(i,j+1))*acosu(j) + v2(i,j) + v2(i+1,j))
         enddo
         te(im,j,k) = 0.25*((u2(im,j)+u2(im,j+1))*acosu(j) + v2(im,j) + v2(1,j))
      enddo

      do j=jfirst,jlast
         do i=1,im
            tm(i,j) = pt(i,j,k)*(1.+r_vir*q(i,j,k))
         enddo
      enddo


      do j=js2g0,jn2g0
         do i=1,im
            te(i,j,k) = delp(i,j,k) * ( te(i,j,k)    +            &
                 cp*tm(i,j) )
         enddo
      enddo

      if ( jfirst == 1 ) then
          te_sp = 0.
          do i=1,im
!          te_sp = te_sp + u2(i,2) + v2(i,2)
! SJL 01/26/04
             te_sp = te_sp + 0.5*u2(i,2)/cose(2)
          enddo

          te_sp =  delp(1,1,k) * (0.5*te_sp/float(im) +           &
               cp*tm(1,1) )
          do i=1,im
             te(i,1,k) = te_sp
          enddo
      endif

      if ( jlast == jm ) then
          te_np = 0.
          do i=1,im
!          te_np = te_np + u2(i,jm) + v2(i,jm-1)
             te_np = te_np + 0.5*u2(i,jm)/cose(jm)
          enddo
          te_np = delp(1,jm,k) * (0.5*te_np/float(im) +          &
               cp*tm(1,jm) )
          do i=1,im
             te(i,jm,k) = te_np
          enddo
      endif

      do j=jfirst,jlast
         do i=1,im
            dz(i,j,k) = rg*tm(i,j)
         enddo
      enddo
   enddo
   call fv_array_sync()

!$omp parallel do private(i,j,k,bte, xsum, gztop)
   do j=jsp,jep !jfirst,jlast
! Perform vertical integration
      do i=1,im
         gztop(i) = phis(i,j)
         do k=1,km
            gztop(i) = gztop(i) + dz(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
         enddo
      enddo

      if (j == 1) then
! SP
          tte(1) = pe(1,km+1,1)*phis(1,1) - pe(1,1,1)*gztop(1)
          do k=1,km
             tte(1) = tte(1) + te(1,1,k)
          enddo
          tte(1)  = acap * tte(1)

      elseif (j == jm) then
! NP
          tte(jm) = pe(1,km+1,jm)*phis(1,jm) - pe(1,1,jm)*gztop(1)
          do k=1,km
             tte(jm) = tte(jm) + te(1,jm,k)
          enddo
          tte(jm) = acap * tte(jm)

      else
! Interior

          do i=1, im
             bte(i) = pe(i,km+1,j)*phis(i,j) - pe(i,1,j)*gztop(i)
          enddo

          do k=1,km
             do i=1,im
                bte(i) = bte(i) + te(i,j,k)
             enddo
          enddo

          xsum = 0.
          do i=1,im
             xsum = xsum + bte(i)
          enddo
          tte(j) = xsum*cosp(j)

      endif
   enddo
   call fv_array_sync()

! call par_vecsum(jm, jfirst, jlast, tte, te0)

   call fv_print_chksums( 'Exiting  benergy' )
 end subroutine benergy


 subroutine pkez(nx, im, km, jfirst, jlast,      &
      pe, pk, akap, ks, peln, pkz, eta)

! !INPUT PARAMETERS:
   integer, intent(in):: nx                   ! SMP decomposition in x
   integer, intent(in):: im, km               ! Dimensions
   integer, intent(in):: jfirst, jlast        ! Latitude strip
   integer, intent(in):: ks
   logical, intent(in):: eta     ! Is data on ETA coordinate?
! True:  input pe;     output pk, pkz, peln
! False: input pe, pk; output pkz and peln
   real, intent(in):: pe(im, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(in):: akap

! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout)::  pk(im,jfirst:jlast,km+1)
   real, intent(inout):: pkz(im,jfirst:jlast,km)

! !OUTPUT
   real, intent(out):: peln(im, km+1, jfirst:jlast)   ! log (pe)

! !DESCRIPTION:
!
!
! !CALLED FROM:
!     te_map and fvgcm
!
! !REVISION HISTORY:
!
!     WS  99.07.27 : Limited region to jfirst:jlast
!     WS  99.10.22 : Deleted cp as argument (was not used)
!     WS  99.11.05 : Documentation; pruning of arguments
!     SJL 00.01.02: SMP decomposition in i
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
   real pk2(im, km+1)
   real pek
   real lnp
   real ak1
   integer i, j, k
   integer ixj, jp, it, i1, i2

   integer :: ixjs, ixje


   call fv_print_chksums( 'Entering  pkez' )
   call fv_array_check( LOC(pe) )
   call fv_array_check( LOC(pk) )
   call fv_array_check( LOC(pkz) )
   call fv_array_check( LOC(peln) )

   it = im / nx
   jp = nx * ( jlast - jfirst + 1 )

   ak1 = (akap + 1.) / akap

   call fv_array_limits( 1, jp, ixjs, ixje )

!$omp  parallel do default(shared) private(ixj, i1, i2, i, j, k, pek, lnp, pk2)
!     do 1000 j=jfirst, jlast
!        i1 = 1
!        i2 = im
   do 1000 ixj=ixjs, ixje !1,jp

      j  = jfirst + (ixj-1) / nx
      i1 =  1 + it * mod(ixj-1, nx)
      i2 = i1 + it - 1

      if ( eta ) then
!---- Eta cordinate Coordinate -----------
          pek = pe(i1,1,j) ** akap

          do i=i1,i2
             pk2(i,1)   = pek
          enddo

          if(ks /= 0) then
              do k=2, ks+1
                 pek = pe(i1,k,j)**akap
                 lnp = log(pe(i1,k,j))
                 do i=i1,i2
                    pk2(i,k)   = pek
                    peln(i,k,j) =  lnp
                 enddo
              enddo

              if( ptop < ptop_min ) then
                  do i=i1,i2
                     peln(i,1,j) = peln(i,2,j) - ak1
                  enddo
              else
                  lnp = log( ptop )
                  do i=i1,i2
                     peln(i,1,j) = lnp
                  enddo
              endif

              do k=1,ks
                 pek = (pk2(i1,k+1) - pk2(i1,k)) / (akap*(peln(i1,k+1,j) - peln(i1,k,j)) )
                 do i=i1,i2
                    pkz(i,j,k) = pek
                 enddo
              enddo
          endif

          do k=ks+2,km
             do i=i1,i2
                pk2(i,k) = pe(i,k,j)**akap
             enddo
          enddo

          do i=i1,i2
             pk2(i,km+1) = pk(i,j,km+1)
          enddo

          do k=ks+2,km+1
             do i=i1,i2
                peln(i,k,j) =  log(pe(i,k,j))
             enddo
          enddo

          if ( ks == 0 ) then
              if( ptop < ptop_min ) then
                  do i=i1,i2
                     peln(i,1,j) = peln(i,2,j) - ak1
                  enddo
              else
                  lnp = log( ptop )
                  do i=i1,i2
                     peln(i,1,j) = lnp
                  enddo
              endif
          endif

          do k=ks+1,km
             do i=i1,i2
                pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k)) /           &
                     (akap*(peln(i,k+1,j) - peln(i,k,j)) )
             enddo
          enddo

          do k=2,km
             do i=i1,i2
                pk(i,j,k) = pk2(i,k)
             enddo
          enddo
      else
! ---- General Coordinate -------

          pek = pk(i1,j,1)
          do i=i1,i2
             pk2(i,1) = pek
          enddo

          do k=2,km+1
             do i=i1,i2
                peln(i,k,j) =  log(pe(i,k,j))
                pk2(i,k) =  pk(i,j,k)
             enddo
          enddo

!---- GFDL modification
          if( ptop < ptop_min ) then
              do i=i1,i2
                 peln(i,1,j) = peln(i,2,j) - ak1
              enddo
          else
              lnp = log( ptop )
              do i=i1,i2
                 peln(i,1,j) = lnp
              enddo
          endif
!---- GFDL modification
          do k=1,km
             do i=i1,i2
                pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                     (akap*(peln(i,k+1,j) - peln(i,k,j)) )
             enddo
          enddo

      endif
1000 end do !continue
   call fv_array_sync()
   call fv_print_chksums( 'Exiting  pkez' )
  end subroutine pkez

end module mapz_module
