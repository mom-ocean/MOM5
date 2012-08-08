module update_fv_phys_mod
  public :: update_fv_phys
contains
  subroutine update_fv_phys ( dt, nq, full_phys, nudge, Time )

    use fv_pack, only: nlon, mlat, beglon, endlon, beglat, endlat, &
         age_tracer,  age_time, do_ch4_chem,   &
         u_ghosted, ng_d, ng_s, ak, bk,        &
         use_tendency, get_eta_level, pft_phys, polavg
#ifdef USE_LIMA
    use fv_pack, only : sc
#endif

    use pft_module, only: pft2d_phys
    use constants_mod,      only: kappa
    use field_manager_mod,  only: MODEL_ATMOS
    use tracer_manager_mod, only: get_tracer_index
    use atmos_nudge_mod,    only: get_atmos_nudge, do_ps
    use time_manager_mod,   only: time_type

#ifdef SPMD
    use mod_comm,     only: gid, mp_send3d, mp_recv3d
#endif
    use mpp_mod, only: mpp_sync
    use fv_arrays_mod, only: fv_stack_push, fv_array_sync, fv_print_chksums

! implicit none
#include "fv_arrays.h"

    integer, intent(in):: nq         ! tracers modified by physics 
    logical, intent(in):: full_phys   ! full physics
    logical, intent(in):: nudge
    type (time_type), intent(in) :: Time


    real, intent(in)   :: dt

! Local arrays:
    real du_s(nlon,nlev)
! real  pk2(nlon,nlev+1)
    real  pk2(nlon,nlev+1) !private, only used in j-parallel loop
    real  dt5
    integer i, j, k, m
    integer cld_amt

!***********
! Haloe Data for CH4
!***********

    real, parameter::    q1_h2o = 2.2E-6
    real, parameter::    q7_h2o = 3.8E-6
    real, parameter::  q100_h2o = 3.8E-6
    real, parameter:: q1000_h2o = 3.1E-6
    real, parameter:: q2000_h2o = 2.8E-6
    real, parameter::   tau_h2o = 120.*86400.

    integer js2g0, jn2g0

    real  phalf(nlev+1), pfull(nlev)
    real  qstar
    real  rdt
    real  dbk
! real, allocatable, dimension(:,:) :: ps_dt
    real :: ps_dt(nlon, beglat:endlat)
#ifdef use_shared_pointers
    pointer( ptr_du_s, du_s )
    pointer( ptr_ps_dt, ps_dt )

#include "fv_point.inc"
    call fv_stack_push( ptr_du_s, nlon*nlev )
    call fv_stack_push( ptr_ps_dt, nlon*(endlat-beglat+1) )
#endif

    call fv_print_chksums( 'Entering  update_fv_phys' )
    dt5 = 0.5 * dt
    rdt = 1./ dt

    if ( use_tendency ) then
!----------------------------------------
! Compute wind tendencies due to physics:
!----------------------------------------
!$omp parallel do private (i,j,k)
        do k=ksp,kep !1,nlev
           do j=beglat,endlat
              do i=1,nlon !beglon,endlon
                 u_dt(i,j,k) = u_dt(i,j,k) - (ua(i,j,k)-u_phys(i,j,k))*rdt
                 v_dt(i,j,k) = v_dt(i,j,k) - (va(i,j,k)-v_phys(i,j,k))*rdt
              enddo
           enddo
        enddo
    endif

    if ( pft_phys .and. full_phys ) then

        js2g0 = max(2,beglat)
        jn2g0 = min(mlat-1,endlat)

! 3-point physics tendency filter: u_dt, v_dt, t_dt
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Need some mods for 2D to work
!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private (k)
        do k=ksp,kep !1,nlev
#if ( !defined USE_LIMA)
           call pft2d_phys(u_dt(1,js2g0,k),  nlon, jn2g0-js2g0+1)
           call pft2d_phys(v_dt(1,js2g0,k),  nlon, jn2g0-js2g0+1)
           call pft2d_phys(t_dt(1,js2g0,k),  nlon, jn2g0-js2g0+1)
#else               
           call pft2d_phys(u_dt(1,js2g0,k),  sc(js2g0), nlon, jn2g0-js2g0+1)
           call pft2d_phys(v_dt(1,js2g0,k),  sc(js2g0), nlon, jn2g0-js2g0+1)
           call pft2d_phys(t_dt(1,js2g0,k),  sc(js2g0), nlon, jn2g0-js2g0+1)
#endif
        enddo

    endif

#if ( !defined MARS_GCM )
    cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')

    call get_eta_level(nlev, 1.0E5, pfull, phalf)
#endif


!$omp parallel do private (i,j,k,m, qstar)
    do k=ksp, kep !1,nlev

#if ( !defined MARS_GCM)
! Do idealized Ch4 chemistry
       if ( do_ch4_chem .and. pfull(k) < 20.E2 ) then

           if ( pfull(k) < 1. ) then
               qstar = q1_h2o
           elseif ( pfull(k) <   7. .and. pfull(k) >=    1. ) then
               qstar = q1_h2o + (q7_h2o-q1_h2o)*log(pfull(k)/1.)/log(7.)
           elseif ( pfull(k) <  100. .and. pfull(k) >=    7. ) then
               qstar = q7_h2o + (q100_h2o-q7_h2o)*log(pfull(k)/7.)/log(100./7.)
           elseif ( pfull(k) < 1000. .and. pfull(k) >=  100. ) then
               qstar = q100_h2o + (q1000_h2o-q100_h2o)*log(pfull(k)/1.E2)/log(10.)
           elseif ( pfull(k) < 2000. .and. pfull(k) >= 1000. ) then
               qstar = q1000_h2o + (q2000_h2o-q1000_h2o)*log(pfull(k)/1.E3)/log(2.)
           endif

           do j=beglat,endlat
              do i=1,nlon !beglon,endlon
                 q_dt(i,j,k,1) = q_dt(i,j,k,1) + (qstar-q(i,j,k,1))/(tau_h2o+dt)
              enddo
           enddo

       endif
#endif 

! Average tendencies for scalars at poles (make it single-valued)
! This is necessary because the the physics/land/ocean use different grid
! than the dynamics at the poles (which has polar caps)
!----------------------------------------
! Needs some mods for phys 2D decomp
!----------------------------------------
       call polavg(t_dt(1,beglat,k), nlon, mlat, beglat, endlat)

#ifdef MARS_GCM
!!!!       call polavg( delp_dt(1,beglat,k), nlon, mlat, beglat, endlat)
!!!!       This is carried out in fv_phys.F90
#endif MARS_GCM

       if ( use_tendency ) then
           do j=beglat,endlat
              do i=1,nlon !beglon,endlon
                 ua(i,j,k) = ua(i,j,k) + dt*u_dt(i,j,k)
                 va(i,j,k) = va(i,j,k) + dt*v_dt(i,j,k)
                 pt(i,j,k) = t_phys(i,j,k) + dt*t_dt(i,j,k)
              enddo
           enddo
           do m=1,nq
              do j=beglat,endlat
                 do i=1,nlon !beglon,endlon
                    q_dt(i,j,k,m) = q_dt(i,j,k,m) - (q(i,j,k,m)-q_phys(i,j,k,m))*rdt
                 enddo
              enddo
           enddo
       else
           do j=beglat,endlat
              do i=1,nlon !beglon,endlon
                 ua(i,j,k) = ua(i,j,k) + dt*u_dt(i,j,k)
                 va(i,j,k) = va(i,j,k) + dt*v_dt(i,j,k)
                 pt(i,j,k) = pt(i,j,k) + dt*t_dt(i,j,k)
              enddo
           enddo
       endif


!----------------
! Update tracers:
!----------------
       do m=1,nq
!----------------------------------------
! Needs some mods for phys 2D decomp
!----------------------------------------
          call polavg(q_dt(1,beglat,k,m), nlon, mlat, beglat, endlat)
          do j=beglat,endlat
             do i=1,nlon !beglon,endlon
                q(i,j,k,m) = q(i,j,k,m) + dt*q_dt(i,j,k,m)
             enddo
          enddo
       enddo


       if ( full_phys ) then

#ifdef MARS_GCM
!   The adjustment of tracer mixing rations to reflect the 
!    change in atmospheric mass is carried out in fv_phys.F90 

#else
           do j=beglat,endlat
              do i=1,nlon !beglon,endlon
!----------------------------------------------------
! Adjust total air mass due to changes in water vaper
! (cloud liquid/ice effects ignored)
!----------------------------------------------------
                 t_dt(i,j,k) = 1. + dt*q_dt(i,j,k,1)
                 delp(i,j,k) = delp(i,j,k) * t_dt(i,j,k)
              enddo
           enddo

           do m=1,ncnst
!      if( m /= cld_amt ) then  ! cloud fraction in GFDL physics
!-----------------------------------------
! Adjust mass mixing ratio of all tracers 
!-----------------------------------------
              do j=beglat,endlat
                 do i=1,nlon !beglon,endlon
                    q(i,j,k,m) = q(i,j,k,m) / t_dt(i,j,k)
                 enddo
              enddo
!      endif
           enddo
#endif 
       endif

    enddo         ! -----------------------k-loop------------------
    call fv_array_sync()

! [delp, (ua, va), pt, q] updated. Perform nudging if requested

!------- nudging of atmospheric variables toward specified data --------
    if (nudge) then
!      allocate ( ps_dt(beglon:endlon, beglat:endlat) )
        ps_dt(:,:) = 0.
!Balaji: get_atmos_nudge is parallelized over i (beglon:endlon of shared arrays)
!--------------------------------------------
! All fields will be updated; tendencies added
!--------------------------------------------
        call get_atmos_nudge ( Time, dt, beglon, endlon, beglat, endlat,    &
             nlev, ng_d, ps(beglon:endlon,:), ua(beglon:endlon,:,:), &
             va(beglon:endlon,:,:), pt(beglon:endlon,:,:), &
             q(beglon:endlon,:,:,:), ps_dt(beglon:endlon,:), u_dt(beglon:endlon,:,:),  & 
             v_dt(beglon:endlon,:,:), t_dt(beglon:endlon,:,:), &
             q_dt(beglon:endlon,:,:,:) )
        call fv_array_sync()
        if (do_ps) then
!$omp parallel do private (i,j,k, dbk)
            do k=ksp,kep !1,nlev
!--------------
! Update delp
!--------------
               dbk = dt * (bk(k+1) - bk(k))
               do j=beglat,endlat
                  do i=1,nlon !beglon,endlon
                     delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
                  enddo
               enddo
            enddo
        endif
!    deallocate ( ps_dt )
    endif
!-----------------------------------------------------------------------
    call fv_array_sync()
#ifdef SPMD
    call mp_send3d(gid+1, gid-1, nlon, mlat, nlev, 1, nlon,  beglat, endlat, &
         1, nlev, 1, nlon, endlat, endlat, 1, nlev, u_dt)
#endif

!$omp parallel do private (i,j,k)
    do k=ksp,kep !1,nlev
!--------------
! Update u-wind
!--------------
       do j=beglat+1,endlat
          do i=1,nlon !beglon,endlon
             u(i,j,k) = u(i,j,k) + dt5*(u_dt(i,j-1,k)+u_dt(i,j,k))
          enddo
       enddo

!---------
! Update v
!---------
       do j=max(2,beglat), min(mlat-1,endlat)
!----------------------------------------
! Needs some mods for phys 2D decomp
!----------------------------------------
          v(1,j,k) = v(1,j,k) + dt5*(v_dt(nlon,j,k)+v_dt(1,j,k))
          do i=2,nlon
             v(i,j,k) = v(i,j,k) + dt5*(v_dt(i-1,j,k)+v_dt(i,j,k))
          enddo
       enddo
    enddo         ! k-loop
    call fv_array_sync()
    if ( full_phys ) then
!----------------------------------------
! Update pe, peln, pkz, and surface winds
!----------------------------------------
!$omp parallel do private(i, j, k, pk2)
        do j=jsp,jep !beglat,endlat

           do i=1,nlon !beglon,endlon
              u_srf(i,j) = ua(i,j,nlev)
              v_srf(i,j) = va(i,j,nlev)
              pk2(i,1) = pk(i,j,1)
           enddo

           do k=2,nlev+1                                                                             
              do i=1,nlon !beglon,endlon
                 pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
                 pk2(i,k)  = pe(i,k,j) ** kappa
                 peln(i,k,j) = log(pe(i,k,j))
              enddo
           enddo

           do i=1,nlon
              ps(i,j) = pe(i,nlev+1,j)
           enddo

           do k=1,nlev                                                                             
              do i=1,nlon !beglon,endlon
                 pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                      (kappa*(peln(i,k+1,j) - peln(i,k,j)) )
              enddo
           enddo
        enddo      ! j-loop
        call fv_array_sync()
    endif

#ifdef SPMD 
    call mp_recv3d(gid-1, nlon, mlat, nlev, 1, nlon, beglat-1, beglat-1, 1, &
         nlev, 1, nlon, beglat-1, beglat-1, 1, nlev, du_s)
    call fv_array_sync()
!----------------
! Finish u update
!----------------
!Balaji: needs work!
    if ( beglat > 1 ) then 
!$omp parallel do private(i,k)
        do k=ksp,kep !1,nlev
           do i=1,nlon !beglon,endlon
              u(i,beglat,k) = u(i,beglat,k) + dt5*(u_dt(i,beglat,k)+du_s(i,k))
           enddo
        enddo
    endif
    call fv_array_sync()
    u_ghosted = .false.
#endif  SPMD 
    call fv_print_chksums( 'Exiting update_fv_phys' )
  end subroutine update_fv_phys

end module update_fv_phys_mod
