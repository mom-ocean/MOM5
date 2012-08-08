#ifndef USE_LIMA
  subroutine init_dry_atm(mountain, kappa, grav, rg)
  use fv_pack,      only: nlon, mlat, beglat, endlat,  &
                          ak, bk, ptop,      &
                          ng_d, ng_s,  &         
                          full_phys,   &
                          d2a3d, p_var, coslon, sinlon, ighost,     &
                          cosp, sinp, cose, sine, dry_mass,  master
#ifdef MARS_GCM
    use fv_pack,   only: p_ref
#endif MARS_GCM


  use fms_mod,      only: file_exist,  error_mesg,  FATAL

#ifdef SPMD
   use mod_comm, only: mp_scatter3d
#endif
   use mpp_io_mod, only: mpp_open, mpp_close
   use mpp_io_mod, only: MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_NATIVE
   use fv_arrays_mod, only: fv_stack_push, fv_array_sync

!  implicit none
#include "fv_arrays.h"

  logical,  intent(in):: mountain
  real, intent(in):: kappa
  real, intent(in):: grav
  real, intent(in):: rg

  integer i,j,k,m
  integer im, jm, km
  real gmean, ftop 
  real:: slp0 = 1.E5

#include "fv_point.inc"

#ifdef POLVANI
 
! This is Lorenzo Polvani's idealized test case
 
  real, allocatable:: u_g(:,:,:), pt_g(:,:,:)
  integer :: unit

!$omp parallel do private (i, j)
  do j=jsp,jep
     do i=1,nlon
        ps(i,j) = slp0
     enddo
  enddo
  call fv_array_sync()

!$omp parallel do private (i, j, k, m)
  do k=ksp,kep
      do j=beglat,endlat
         do i=1,nlon
            v(i,j,k) = 0.
            delp(i,j,k) = ak(k+1) - ak(k) + ps(i,j) * (bk(k+1) - bk(k))
         enddo
      enddo
   enddo
   call fv_array_sync()

!Balaji: should not use unit "61"... convert to use mpp_open
   call mpp_open( unit, action=MPP_RDONLY, form=MPP_NATIVE, threading=MPP_MULTI )
   read(unit) im, jm, km
   write(*,*) im, jm, km
   if( im /= nlon .or. jm /= mlat .or. km /= nlev ) then
       call error_mesg('init_dry_atm','Dimension inconsistent', FATAL)
   endif
   allocate (  u_g(im,jm,km) )
   allocate ( pt_g(im,jm,km) )
   read(unit) u_g
   read(unit) pt_g
   call mpp_close(unit)
   do k = ksp,kep
      do j = jsd,jed
         do i =1,nlon
            u(i,j,k) = u_g(i,j,k)
            pt(i,j,k) = pt_g(i,j,k)
         end do
      end do
   end do
   call fv_array_sync
   deallocate (  u_g )
   deallocate ( pt_g )
   
! Read (u,PT) data from Palvani
!      if ( master ) then
!           if ( mlat == 91 ) then
!               open(unit=61,file='/home/sjl/ICs/ut_Polvani.N45L20',     &
!                    form='unformatted', status='old')
!           elseif ( mlat == 181 ) then
!               open(unit=61,file='/home/sjl/ICs/ut_Polvani.N90L20',     &
!                    form='unformatted', status='old')
!           elseif ( mlat == 361 ) then
!               open(unit=61,file='/home/sjl/ICs/ut_Polvani.N180L20',     &
!                    form='unformatted', status='old')
!           endif
!           read(61) im,jm,km
!           write(*,*) im, jm, km
!           if( im /= nlon .or. jm /= mlat .or. km /= nlev ) then
!               call error_mesg('init_dry_atm','Dimension inconsistent', FATAL)
!           endif
!      endif
!
!#ifdef SPMD
!      if ( master ) then
!           allocate (  u_g(im,jm,km) )
!           allocate ( pt_g(im,jm,km) )
!           read(61) u_g
!           read(61) pt_g
!           close(61)
!      endif
!! Scatter data
!        call mp_scatter3d(u_g, u, nlon, mlat, nlev, beglat, endlat,  &
!                          1, nlev, ng_d, ng_s, 0)
!        call mp_scatter3d(pt_g, pt, nlon, mlat, nlev, beglat, endlat,  &
!                          1, nlev, ng_d, ng_d, 0)
!
!      if ( master ) then
!           deallocate (  u_g )
!           deallocate ( pt_g )
!      endif
!#else
!           read(61) u
!           read(61) pt
!           close(61)
!#endif

! Initialize tracers
           q  = 0.
     do m=1,ncnst
        do j=jsp,jep
           do i=1,nlon
              q(i,j,nlev,m) = 1.
           enddo
         enddo
      enddo
      call fv_array_sync()

  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, .false. )

#else

#ifdef LIN2004
  real dpt, pi, dl
  real xc, yc, br, xo, yo, zo, xg, yg, zg
  real zam, prc, dist
  real u2(mlat,nlev), t2(mlat,nlev)
#ifdef use_shared_pointers
  pointer( p_u2, u2 )
  pointer( p_t2, t2 )
  call fv_stack_push( p_u2, mlat*nlev )
  call fv_stack_push( p_t2, mlat*nlev )
#endif
!$omp parallel do private (i, j)
  do j=jsp,jep
     do i=1,nlon
        ps(i,j) = slp0
     enddo
  enddo
  call fv_array_sync()

!$omp parallel do private (i, j, k, m)
  do k=ksp,kep
      do j=beglat,endlat
         do i=1,nlon
            v(i,j,k) = 0.
            delp(i,j,k) = ak(k+1) - ak(k) + ps(i,j)*(bk(k+1) - bk(k))
         enddo
      enddo
   enddo
   call fv_array_sync()

  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, .false. )

! Define jet and temperature as in Lin 2004, MWR
  call jet2d_symm(mlat, nlev, u2, t2, ak, bk, slp0, cosp, sinp, cose, sine)

!$omp parallel do private (i, j, k, m)
  do k=ksp,kep
      do j=beglat,endlat
         do i=1,nlon
            u(i,j,k) = u2(j,k)
         enddo
      enddo
   enddo
   call fv_array_sync()

!----------------------------------
! Define temperature perturbation 
!----------------------------------
    dpt = 2. / (slp0**(2./7.))           ! 2 deg (k) thermal bubble
    pi = 4.*atan(1.)
    yc = pi/4.            ! center in lat.
    xc = .5*pi            ! center: longitude
    br = 2.*pi / 36       ! 10 degree width
    dl = 2.*pi / real(nlon)

!
! Perturbation center location in (x, y, z) Cartensian coordinate: unit sphere
!
    xo = cos(yc) * cos(xc)
    yo = cos(yc) * sin(xc)
    zo = sin(yc)

   do k=1,nlev
      do j=beglat,endlat
         do i=1,nlon
            zam  = dl*real(i-1)
            prc = 0.5*(pe(i,k,j)+pe(i,k+1,j))
!
! Great circle distance on unit sphere
!
            xg = cosp(j) * coslon(i)
            yg = cosp(j) * sinlon(i)
            zg = sinp(j)
            dist = acos(xg*xo + yg*yo + zg*zo)
            pt(i,j,k) = t2(j,k) + pkz(i,j,k)*dpt*exp(-(dist/br)**2 )   &
                        * max(0., (prc-800.E2)/200.E2)
         enddo
      enddo
   enddo

! Initialize tracers
     q = 0.
     do m=1,ncnst
        do j=beglat, endlat
           do i=1,nlon
              q(i,j,nlev,m) = 1.
           enddo
         enddo
     enddo
#else

     ftop = gmean(nlon, mlat, beglat, endlat, phis(1,beglat))
     if(master) write(6,*) 'mean terrain height (m)=', ftop/grav

!$omp parallel do private (i, j, k, m)
     do k=ksp,kep
        do j=beglat, endlat
           do i=1,nlon
              q(i,j,k,1) = 3.E-6      ! specific humidity in the strat
           enddo
        enddo

        if( ncnst > 1 ) then
            do m=2,ncnst
               do j=beglat, endlat
                  do i=1,nlon
                     q(i,j,k,m) = 0.
                  enddo
               enddo
            enddo
        endif
     enddo

! -- TEST code for tracer mass conservation -------------
#ifdef TEST_TRACER
        q = 0.
       do j=beglat, endlat
          do i=1,nlon
             q(i,j,nlev,1) = 1.e-6
             q(i,j,nlev,2) = 1.e-12
             q(i,j,nlev,3) = 1.e-24
             q(i,j,nlev,4) = 1.
          enddo
       enddo
#endif

#ifdef MARS_GCM
  call hydro_eq(nlon, mlat, nlev, beglat, endlat, ps, phis, dry_mass,  &
                delp, ak, bk, u, v, pt, ng_d, ng_s, 0, grav, rg,   &
                p_ref, mountain, master)

!!!  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
!!!               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, adjust_dry_mass)

  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, .false. )


#else
  call hydro_eq(nlon, mlat, nlev, beglat, endlat, ps, phis, dry_mass,  &
                delp, ak, bk, u, v, pt, ng_d, ng_s, 0, grav, rg,   &
                mountain, master)

  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, full_phys)
#endif MARS_GCM

  call d2a3d(u(1,beglat-ng_d,nlev), v(1,beglat-ng_d,nlev), u_srf, v_srf, &
             nlon, mlat, 1, beglat, endlat, ng_d, ng_s, coslon, sinlon)

  do j=beglat,endlat
     do i=1,nlon
        ps_bp(i,j) = ps(i,j)
     enddo
  enddo

#endif  LIN2004

#endif  POLVANI

 end subroutine init_dry_atm

#ifdef MARS_GCM
 subroutine hydro_eq(im, jm, km, beglat, endlat, ps, hs, drym,  &
                    delp, ak, bk, u, v, pt, ng_d, ng_s, ig, grav,  &
                    rg, pref, mountain, master)
#else
 subroutine hydro_eq(im, jm, km, beglat, endlat, ps, hs, drym,  &
                    delp, ak, bk, u, v, pt, ng_d, ng_s, ig, grav,  &
                    rg, mountain, master)
#endif MARS_GCM
!============================================================

   use fv_arrays_mod, only: fv_array_check, fv_array_sync
   use fv_arrays_mod, only: isp, iep, jsp, jep, ksp, kep
   use pmaxmin_mod, only: pmaxming

  implicit none
! Input:
  integer, intent(in):: im, jm, km
  integer, intent(in):: ng_d, ng_s
  integer, intent(in):: beglat, endlat, ig
  real, intent(in):: ak(km+1), bk(km+1)
  real, intent(in):: hs(im,beglat-ig:endlat+ig)
  real, intent(in):: grav, rg, drym
  logical, intent(in):: mountain
  logical, intent(in):: master
#ifdef MARS_GCM
  real, intent(in):: pref
#endif MARS_GCM

! Output
  real, intent(out)::  ps(im,beglat:endlat)
  real, intent(out)::  u(im,beglat-ng_d:endlat+ng_s,km)
  real, intent(out)::  v(im,beglat-ng_d:endlat+ng_d,km)
  real, intent(out):: pt(im,beglat-ng_d:endlat+ng_d,km)
  real, intent(out):: delp(im,beglat:endlat,km)

! Local
  integer  i,j,k
  real z1, t1, p1, t0, a0, psm, dps
  real gmean, ztop, p_top
  real gz(im,km+1)
  real ph(im,km+1)
  integer kt
  real:: mslp = 1010.6*100.

  call fv_array_check( LOC(hs) )
  call fv_array_check( LOC(ps) )
  call fv_array_check( LOC(u) )
  call fv_array_check( LOC(v) )
  call fv_array_check( LOC(pt) )
  call fv_array_check( LOC(delp) )

  t1 = 215.
  if ( mountain ) then
! First guess: Mean sel-level pressure 1010.6
#ifdef MARS_GCM 
     print *, 'Mars topo -> ps calc:  ',  pref 
     do j=beglat, endlat
         do i=1,im
            ps(i,j) = pref*exp( -hs(i,j)/( rg*190.0 ) )
         enddo
     enddo

#else

! Given p1 and z1 (250mb, 10km)
        p_top = 10.*100.
        p1 = 25000.
        z1 = 10.E3 * grav
        t0 = 280.            ! sea-level temp.
        a0 = (t1-t0)/z1
        ztop = z1 + (rg*t1)*log(p1/p_top)

     do j=beglat, endlat
        do i=1,im
!          ps(i,j) = p1*((z1+t0/a0)/(hs(i,j)+t0/a0))**(1./(a0*rg))  
           ps(i,j) = mslp*((t0/a0)/(hs(i,j)+t0/a0))**(1./(a0*rg))  
        enddo
     enddo
#endif MARS_GCM 

     psm = gmean(im, jm, beglat, endlat, ps(1,beglat))
     dps = drym - psm
     if(master) write(6,*) 'Computed mean ps=', psm
     if(master) write(6,*) 'Correction delta-ps=', dps

!$omp parallel do private (i, j)
        do j=jsp,jep
           do i=1,im
              ps(i,j) = ps(i,j) + dps
           enddo
        enddo

! Compute surface pressure
!$omp parallel do private (i, j, k, ph, gz, kt)
     do j=jsp,jep
#ifdef MARS_GCM
        do k= 1, km
           do i=1,im
              pt(i,j,k) = 190.0
           enddo
        enddo
#else
        do i=1,im
           gz(i,km+1) = hs(i,j)
        enddo

        do k=1,km+1
           do i=1,im
              ph(i,k) = ak(k) + bk(k)*ps(i,j)
           enddo
        enddo

        do k=1,km
           do i=1,im
                 pt(i,j,k) = t1
              if ( ph(i,k+1) <= p1 ) then
                 kt = k
                 gz(i,k) = ztop + rg*t1*log(p_top/ph(i,k))
              else
! Constant lapse rate region (troposphere)
                 gz(i,k) = (hs(i,j)+t0/a0)/(ph(i,k)/ph(i,km+1))**(a0*rg)  &
                         - t0/a0
              endif
           enddo
        enddo

        do k=kt+1,km
           do i=1,im
! Convert geopotential to Temperature
              pt(i,j,k) = (gz(i,k)-gz(i,k+1))/(rg*(log(ph(i,k+1)/ph(i,k))))
              pt(i,j,k) = max(min(t0, pt(i,j,k)), t1) 
           enddo
        enddo
#endif  MARS_GCM

     enddo

  else
!$omp parallel do private (i, j)
      do j=jsp,jep

#ifdef MARS
         do i=1,im
            ps(i,j) = pref
         enddo
         do k=1,km
            do i=1,im
               pt(i,j,k) = 190.0
            enddo
         enddo

!       Add a small random perturbation to the temp field at the surface
         DO i= 1, im
           pt(i,j,km)= pt(i,j,km) + 1.E-5 * cos( (i+j)*1.0 ) 
         ENDDO 
#else 
         do i=1,im
            ps(i,j) = 1.E5
         enddo
         do k=1,km
            do i=1,im
               pt(i,j,k) = t1
            enddo
         enddo
#endif MARS

      end do
    endif  !  ------- end of no-topography case ------------


    call fv_array_sync()

!$omp parallel do private (i, j, k)
      do k=ksp,kep
         do j=beglat,endlat
            do i=1,im
                  u(i,j,k) = 0.
                  v(i,j,k) = 0.
               delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
            enddo
         enddo
      enddo
      call fv_array_sync()

  call pmaxming('Computed T ', pt, im, jm, km, beglat, endlat, ng_d, ng_d, 1.)

 end subroutine hydro_eq


#ifdef LIN2004
 subroutine jet2d_symm(jm, km, u2, t2, ak, bk, p0, cosp, sinp, cose, sine)
!=========================================================================
!
! Compute 2D zonally balanced flow (U,T)
!
 use constants_mod,  only: radius, omega, grav, rdgas
 implicit none

  
 integer, intent(in):: jm, km
 real, intent(in):: p0
 real, intent(in):: ak(km+1), bk(km+1)
 real, intent(in):: cosp(jm), sinp(jm)
 real, intent(in):: cose(jm), sine(jm)

 real, intent(out):: u2(jm,km)
 real, intent(out):: t2(jm,km)

!------
! local:
!------
 real pe1(km+1)
 real up(jm,km+1)
 real ue(jm,km+1)
 real gz(jm,km+1)
 integer j, k
 real:: u0 = 35.
 real z, uz, sin2p, sin2e
 real p1, t1, t0, a0, gz1, vort
 real dp

  do k=1,km+1
     pe1(k) = ak(k) + bk(k)*p0
     z = log(p0/pe1(k))
     uz = u0*z*exp( -z**2/4. )
     do j=1,jm
!
! Y-functional form: u = u(z) * sin(theta)**2
!
        sin2p   = 2.*sinp(j)*cosp(j)
        up(j,k) = uz*sin2p**4 + 0.1*z*(uz-15.)*cosp(j)**8
     enddo

     do j=2,jm
        sin2e = 2.*sine(j)*cose(j)
        ue(j,k) = uz*sin2e**4 + 0.1*z*(uz-15.)*cosp(j)**8
     enddo
  enddo

! Compute gz at p1 (250mb)
     p1 = 250.E2
     t1 = 210.
     t0 = 260.
     a0 = -5.E-3/grav
! height at 250-mb
     gz1 = (t0/a0) / (p1/p0)**(a0*rdgas) - t0/a0

! gz at ptop = gz1 + rdgas*t1*log(p1/ptop)

     do j=1,jm
        gz(j,km+1) = 0.        ! assuming zs=0
     enddo

! j=1
     do k=1,km
        if(pe1(k) <= p1) then
           gz(1,k) = gz1 + rdgas*t1*log(p1/pe1(k))
        else
           gz(1,k) = (t0/a0)/(pe1(k)/p0)**(a0*rdgas) - t0/a0
        endif
     enddo

! Compute gz from SP to NP

     dp = 4.*atan(1.) / (jm-1)

     do k=1,km+1
        do j=2,jm
           vort  = 2.*omega*sine(j) - (up(j,k)*cosp(j) -    &
                   up(j-1,k)*cosp(j-1)) / (radius*cose(j)*dp)
           gz(j,k) = gz(j-1,k) - ue(j,k)*(radius*vort*dp +   &
                                       up(j,k)-up(j-1,k))
        enddo
     enddo

! Compute pt (temperature)
     do k=1,km
        do j=2,jm
           u2(j,k) = 0.5*(ue(j,k)+ue(j,k+1))
        enddo
        do j=1,jm
           t2(j,k) = (gz(j,k)-gz(j,k+1)) / (rdgas*(log(pe1(k+1)/pe1(k))))
        enddo
     enddo

 end subroutine jet2d_symm
#endif

#else

  subroutine init_dry_atm(mountain, kappa, grav, rg)

  use shr_kind_mod, only: r8 => shr_kind_r8
!  use fv_pack,      only: nlon, mlat, nlev, ncnst, beglat, endlat,  &
!                          ps, u, v, delp, pt, q, ak, bk, ptop,      &
!                          pe, peln, pk, pkz, ng_d, ng_s, phis,      &         
!                          ps_bp, ua, va, u_srf, v_srf, full_phys,   &
!                          d2a3d, p_var, coslon, sinlon, ighost,     &
!                          cosp, sinp, cose, sine, dry_mass,  master
  use fv_pack,      only: nlon, mlat, beglat, endlat,  &
                          ak, bk, ptop,      &
                          ng_d, ng_s,  &         
                          full_phys,  master,   &
                          p_var, ighost
  use fms_mod,      only: file_exist,  error_mesg,  FATAL

#ifdef SPMD
   use mod_comm, only: mp_scatter3d
#endif
   use mpp_io_mod, only: mpp_open, mpp_close
   use mpp_io_mod, only: MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_NATIVE
   use fv_arrays_mod, only: fv_stack_push, fv_array_sync

!  implicit none
#include "fv_arrays.h"

  logical,  intent(in):: mountain
  real(r8), intent(in):: kappa
  real(r8), intent(in):: grav
  real(r8), intent(in):: rg

  real(r8):: slp0 = 1.E5
  real(r8):: tref = 273.16
  real(r8):: drym = 98288.        ! Dry air mass correct for USGS
! real(r8):: drym = 98290.        ! Revised dry air mass correct for USGS
  integer i,j,k,m
  integer im, jm, km
  real(r8) gmean, ftop 

#ifdef POLVANI
 
! This is Lorenzo Polvani's idealized test case
!
  real(r8), allocatable:: u_g(:,:,:), pt_g(:,:,:)

!$omp parallel do private (i, j)
  do j=beglat, endlat
     do i=1,nlon
        ps(i,j) = slp0
     enddo
  enddo

!$omp parallel do private (i, j, k, m)
  do k=1,nlev
      do j=beglat,endlat
         do i=1,nlon
            v(i,j,k) = 0.
            delp(i,j,k) = ak(k+1) - ak(k) + ps(i,j) * (bk(k+1) - bk(k))
         enddo
      enddo
  enddo

! Read (u,PT) data from Palvani
      if ( master ) then
           if ( mlat == 91 ) then
               open(unit=61,file='/home/sjl/ICs/ut_Polvani.N45L20',     &
                    form='unformatted', status='old')
           elseif ( mlat == 181 ) then
               open(unit=61,file='/home/sjl/ICs/ut_Polvani.N90L20',     &
                    form='unformatted', status='old')
           elseif ( mlat == 361 ) then
               open(unit=61,file='/home/sjl/ICs/ut_Polvani.N180L20',     &
                    form='unformatted', status='old')
           endif
           read(61) im,jm,km
           write(*,*) im, jm, km
           if( im /= nlon .or. jm /= mlat .or. km /= nlev ) then
               call error_mesg('init_dry_atm','Dimension inconsistent', FATAL)
           endif
      endif

#ifdef SPMD
      if ( master ) then
           allocate (  u_g(im,jm,km) )
           allocate ( pt_g(im,jm,km) )
           read(61) u_g
           read(61) pt_g
           close(61)
      endif
! Scatter data
        call mp_scatter3d(u_g, u, nlon, mlat, nlev, beglat, endlat,  &
                          1, nlev, ng_d, ng_s, 0)
        call mp_scatter3d(pt_g, pt, nlon, mlat, nlev, beglat, endlat,  &
                          1, nlev, ng_d, ng_d, 0)

      if ( master ) then
           deallocate (  u_g )
           deallocate ( pt_g )
      endif
#else
           read(61) u
           read(61) pt
           close(61)
#endif

! Initialize tracers
           q  = 0.
     do m=1,ncnst
        do j=beglat, endlat
           do i=1,nlon
              q(i,j,nlev,m) = 1.
           enddo
         enddo
     enddo

  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, .false. )

#else

     ftop = gmean(nlon, mlat, beglat, endlat, phis(1,beglat))
     if(master) write(6,*) 'mean terrain height (m)=', ftop/grav

     call hydro_eq(nlon, mlat, nlev, beglat, endlat, ps, phis, drym,  &
                   delp, ak, bk, ighost, grav, rg, mountain, master)
     u_srf = 0.; v_srf = 0.

!$omp parallel do private (i, j, k, m)
  do k=1,nlev
      do j=beglat,endlat
         do i=1,nlon
            u(i,j,k) = 0.
            v(i,j,k) = 0.
         enddo
      enddo

! Initialize tracers
        do j=beglat, endlat
           do i=1,nlon
              q(i,j,k,1) = 3.E-6      ! specific humidity in the strat
           enddo
         enddo

     if( ncnst > 1 ) then
     do m=2,ncnst
        do j=beglat, endlat
           do i=1,nlon
              q(i,j,k,m) = 0.
           enddo
         enddo
     enddo
     endif
  enddo

! -- TEST code for tracer mass conservation -------------
#ifdef TEST_TRACER
        q = 0.
       do j=beglat, endlat
          do i=1,nlon
             q(i,j,nlev,1) = 1.
!            q(i,j,nlev,1) = 1.e-24
!            q(i,j,nlev,2) = 1.e-12
!            q(i,j,nlev,3) = 1.e-6
!            q(i,j,nlev,4) = 1.
          enddo
       enddo
#endif

  call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,   &
               pe, peln, pk, pkz, kappa, q, ng_d, ncnst, full_phys )

!$omp parallel do private (i, j, k, m)
  do k=1,nlev
      do j=beglat,endlat
         do i=1,nlon
            pt(i,j,k) = tref
         enddo
      enddo
  enddo
#endif

end subroutine init_dry_atm


subroutine hydro_eq(im, jm, km, beglat, endlat, ps, hs, drym,  &
                    delp, ak, bk, ig, grav, rg, mountain, master)
  implicit none
! Input:
  integer, intent(in):: im, jm, km
  integer, intent(in):: beglat, endlat, ig
  real, intent(in):: ak(km+1), bk(km+1)
  real, intent(in):: hs(im,beglat-ig:endlat+ig)
  real, intent(in):: grav, rg, drym
  logical, intent(in):: mountain
  logical, intent(in):: master

! Output
  real, intent(out)::  ps(im,beglat:endlat)
  real, intent(out):: delp(im,beglat:endlat,km)

! Local
      integer  i,j,k

      real z1, t1, p1, t0, a0, psm, dps
      real gmean

      if ( mountain ) then

! Given p1 and z1 (250mb, 10km)
        p1 = 25000.
        z1 = 10.E3 * grav
        t1 = 212.
        t0 = 280.            ! sea-level temp.
        a0 = (t1-t0)/z1
! Compute surface pressure
!$omp parallel do private (i, j)
        do j=beglat, endlat
           do i=1,im
              ps(i,j) = p1*((z1+t0/a0)/(hs(i,j)+t0/a0))**(1./(a0*rg))  
           enddo
        enddo

!---------
! Adjust ps
!---------
     psm = gmean(im, jm, beglat, endlat, ps(1,beglat))
     if(master) write(6,*) 'Computed mean ps=', psm
     dps = drym - psm

!$omp parallel do private (i, j)
        do j=beglat, endlat
           do i=1,im
              ps(i,j) = ps(i,j) + dps
           enddo
        enddo

      else
!$omp parallel do private (i, j)
        do j=beglat, endlat
           do i=1,im
              ps(i,j) = 1.E5
           enddo
        enddo
      endif

!$omp parallel do private (i, j, k)
      do k=1,km
         do j=beglat,endlat
            do i=1,im
               delp(i,j,k) = ak(k+1) - ak(k) + ps(i,j) * (bk(k+1) - bk(k))
            enddo
         enddo
      enddo

end subroutine hydro_eq
#endif
