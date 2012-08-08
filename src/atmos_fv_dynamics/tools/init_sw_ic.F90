#ifndef USE_LIMA
  subroutine init_sw_ic( icase, im, jm, km, nq, jfirst, jlast,   &
                         ng_d,  ng_s,  phis, u, v, pt, delp, q,  &
                         grav,  ae,   omega, f0, ig )

!-----------------------------------
! Shallow water equivalent settings:
!-----------------------------------

! del(WZ) = CP * PT * del(PK)  = delp
! CP = kappa = 1; pt = 1 ---> Pk == pe;  del(g*Z) == delp
! Example: del(WZ) = 9.81 * 10 (km) ==> delp ~ 981 (mb)


  use fv_pack,      only: sinp, cosp, cose, sine, sinlon, coslon, master
  use pv_module,    only: vort_d

#ifdef SPMD
  use mod_comm,     only: mp_send4d_ns, mp_recv4d_ns
#endif
  implicit none

     integer, intent(in):: ig
     integer, intent(in):: icase           ! Case-2:
                                           ! Case-5: geostrophic flow over isolated mountain
                                           ! Case-6: Rossby-Haurwitz wave#4
     integer, intent(in):: im, jm, km, nq
     integer, intent(in):: jfirst, jlast
     integer, intent(in):: ng_d, ng_s

     real, intent(in):: grav
     real, intent(in):: ae
     real, intent(in):: omega

     real, intent(inout):: f0(im,jfirst-ng_d:jlast+ng_d)
     real, intent(out):: phis(im,jfirst-ig:jlast+ig)   ! surface geopotential
     real, intent(out):: delp(im,jfirst:jlast,km) ! pressure thickness (pascal)

     real, intent(out):: u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-wind (m/s)
     real, intent(out):: v(im,jfirst-ng_d:jlast+ng_d,km)  ! v-wind (m/s)
     real, intent(out):: pt(im,jfirst-ng_d:jlast+ng_d,km)
     real, intent(out):: q(im,jfirst-ng_d:jlast+ng_d,km,nq)

! local:
     integer i,j
     real  vort(im,jfirst:jlast)
     real  pi, dp, dl
     real  alpha, u0
     real  v0, h0, hs0, br, zc, pc, phi, zam, tmp, sr
     real  omg, k, R, cosz, zamc, a, b, c
     real hsp, deg, ddeg
     real h1(jm), u1(jm)

     pi = 4. * atan(1.)
     dl = 2.*pi / im
     dp = pi / (jm-1)

     select case (icase)

     case(2)

      alpha = 0.5*pi
      u0 = 2.*pi*ae / (12.*24.*3600.)
      h0 = 2.94e4

      do j = max(1,jfirst-ng_d), min(jm,jlast+ng_d)
         do i=1,im
            f0(i,j) = 2.*omega*(sinp(j)*cos(alpha)-coslon(i)*cosp(j)*sin(alpha)) 
         enddo
      enddo

      do j=jfirst, jlast
         do i=1,im
            u(i,j,1) = u0*(cose(j)*cos(alpha)+sine(j)*coslon(i)*sin(alpha))
         enddo
      enddo

      DO j=jfirst,jlast
         DO i=1,im
            zamc = dl*real(i-1)
             v(i,j,1) = -u0*sin(zamc)*sin(alpha)
          delp(i,j,1) = h0 - (ae*omega + 0.5*u0)*u0 *        &
                  (sinp(j)*cos(alpha)-coslon(i)*cosp(j)*sin(alpha))**2
            phis(i,j) = 0.
         enddo
      enddo

     case (5)                 ! Williamson et al's test case #5
!----------------------------
! Flow over isolated mountain
!----------------------------
      if(master) write(6,*) 'Generating IC for Flow over isolated mountain case'
      v0  = 20.
      h0  = 5960.*grav
      hs0 = 2000.*grav

      br  = pi/9.
      zc  = .5*pi - pi
      pc  = pi/6.

      do j=jfirst,jlast
         phi = - pi*0.5 + (j-1)*dp
         do i=1,im
            zam = dl*(real(i)-0.5) - pi
            tmp = (zam-zc)**2 + (phi-pc)**2
            sr = min(br**2, tmp)
            phis(i,j) = hs0*(1.-sqrt(sr)/br)
             u(i,j,1)  =  v0*cose(j)
             v(i,j,1)  =  0.
! Total geopotential
          delp(i,j,1)  =  h0 - (ae*omega+0.5*v0)*v0*sinp(j)**2
         enddo
      enddo

     case(6)                  ! Williamson et al's test case #6
!-----------------------
! Rossby-Haurwitz wave#4
!-----------------------

      if(master) write(6,*) 'Generating IC for Rossby-Haurwitz wave#4 case'

      omg = 7.848E-6
      k = omg
      R = 4.
      h0 = 8.E3 * grav

      do j=jfirst, jlast
         do i=1,im
                zam  = dl*(real(i)-0.5)
            u(i,j,1) = AE*omg*cose(j) +        &
             AE*k*(cose(j)**(R-1))*(R*sine(j)**2-cose(j)**2)*cos(R*zam)
         enddo
      enddo

      do j=jfirst, jlast
         cosz = cosp(j)
         do i=1,im
            zam  = DL*(real(i)-0.5)
            zamc = DL*real(i-1)
            v(i,j,1) = -AE*k*R*sinp(j)*sin(R*zamc)*cosz**(R-1)
            A = 0.5*omg*(2.*OMEGA+omg)*cosz**2 + 0.25*k*k*cosz**(R+R) *  &
                ((R+1)*cosz**2 + (2.*R*R-R-2.))-0.5*(k*R)**2*cosz**(R+R-2)
            B = 2.*(OMEGA+omg)*k / ((R+1)*(R+2)) *(cosz**R) *                  &
                ( (R*R+2.*R+2.) - ((R+1.)*cosz)**2 )
            C = 0.25*k*k* ((R+1)*cosz**2-(R+2)) * cosz**(2.*R)
           delp(i,j,1) = h0 + AE*AE*(A+B*cos(R*zam)+C*cos(2.*R*zam))
           phis(i,j) = 0.
         enddo
      enddo

     case(8)
!-----------------------------------------
! Vortex breaking (Lin and Rood 1997, QJ)
!-----------------------------------------

      if(master) write(6,*) 'Generating IC for Rossby wave breaking case'


      hsp = 6000. * grav
      ddeg = 180./real(jm-1)

      do j=1,jlast
         deg = -90. + (real(j-1)-0.5)*ddeg
         if(deg <= 0.) then
            u1(j) = -10.*(deg+90.)/90.
         elseif(deg <= 60.) then
            u1(j) = -10. +  deg
         else
            u1(j) = 50. - (50./30.)* (deg - 60.)
         endif
      enddo

      h1(1) = hsp
      do j=2,jlast
         h1(j) = h1(j-1) - dp*sine(j)*(AE*2.*omega+u1(j)/cose(j))*u1(j)
      enddo
      
      do j=jfirst, jlast
         do i=1,im
            u(i,j,1) = u1(j)
            v(i,j,1) = 0.
            delp(i,j,1) = h1(j)
            phis(i,j) = 0.
         enddo
      enddo

     end select

! For all cases:

#ifdef SPMD
     call mp_send4d_ns(im, jm, 1, 1, jfirst, jlast, 1, 1, ig, ig, phis)
     call mp_recv4d_ns(im, jm, 1, 1, jfirst, jlast, 1, 1, ig, ig, phis)
#endif

!--------------------------------------
! Initialize the passive tracer with PV
!--------------------------------------
     call vort_d(im, jm, 1, jfirst, jlast, u, v, vort, ng_s, ng_d)

     do j=jfirst, jlast
        do i=1,im
! Substract out surface geopotential to get Layer thickness
           delp(i,j,1)   =  delp(i,j,1) - phis(i,j)
              q(i,j,1,1) = (f0(i,j) + vort(i,j)) / delp(i,j,1)
             pt(i,j,1)   =  1.
        enddo
     enddo

  end subroutine init_sw_ic
#else  
  subroutine init_sw_ic( icase, im, jm, km, nq, jfirst, jlast,   &
                         ng_d,  ng_s,  phis, u, v, pt, delp, q,  &
                         grav,  ae,   omega, f0, ig )

!-----------------------------------
! Shallow water equivalent settings:
!-----------------------------------

! del(WZ) = CP * PT * del(PK)  = delp
! CP = kappa = 1; pt = 1 ---> Pk == pe;  del(g*Z) == delp
! Example: del(WZ) = 9.81 * 10 (km) ==> delp ~ 981 (mb)


  use shr_kind_mod, only: r8 => shr_kind_r8
  use fv_pack,      only: sinp, cosp, cose, sine, master
  use pv_module,    only: vort_d

#ifdef SPMD
  use mod_comm,     only: mp_send4d_ns, mp_recv4d_ns
#endif
  implicit none

     integer, intent(in):: ig
     integer, intent(in):: icase           ! Case-2:
                                           ! Case-5: geostrophic flow over isolated mountain
                                           ! Case-6: Rossby-Haurwitz wave#4
     integer, intent(in):: im, jm, km, nq
     integer, intent(in):: jfirst, jlast
     integer, intent(in):: ng_d, ng_s

     real(r8), intent(in):: grav
     real(r8), intent(in):: ae
     real(r8), intent(in):: omega
     real(r8), intent(in):: f0(jfirst-ng_d:jlast+ng_d)

     real(r8), intent(out):: phis(im,jfirst-ig:jlast+ig)   ! surface geopotential
     real(r8), intent(out):: delp(im,jfirst:jlast,km) ! pressure thickness (pascal)

     real(r8), intent(out):: u(im,jfirst-ng_d:jlast+ng_s,km)  ! u-wind (m/s)
     real(r8), intent(out):: v(im,jfirst-ng_d:jlast+ng_d,km)  ! v-wind (m/s)
     real(r8), intent(out):: pt(im,jfirst-ng_d:jlast+ng_d,km)
     real(r8), intent(out):: q(im,jfirst-ng_d:jlast+ng_d,km,nq)

! local:
     integer i,j
     real(r8)  vort(im,jfirst:jlast)
     real(r8)  pi, dp, dl
     real(r8)  alpha, u0
     real(r8)  v0, h0, hs0, br, zc, pc, phi, zam, tmp, sr
     real(r8)  omg, k, R, cosz, zamc, a, b, c

     pi = 4. * atan(1.)
     dl = 2.*pi / im
     dp = pi / (jm-1)

     select case (icase)

     case(2)
      alpha = 0.
      u0 = 2.*pi*ae / (12.*24.*3600.)
      h0 = 2.94e4

      do j=jfirst, jlast
         do i=1,im
                  zam = dl*(i-1) - PI
            u(i,j,1) = u0*(cose(j)*cos(alpha)+sine(j)*cos(zam)*sin(alpha))
         enddo
      enddo

      DO j=jfirst,jlast
         DO i=1,im
            zam  = DL*(i-1) - PI
            zamc = zam - 0.5*DL
             v(i,j,1) = -u0*sin(zamc)*sin(alpha)
          delp(i,j,1) = h0 - (ae*omega + 0.5*u0)*u0 *        &
                  (sinp(j)*cos(alpha)-cos(zam)*cosp(j)*sin(alpha))**2
            phis(i,j) = 0.
         enddo
      enddo

     case (5)                 ! Williamson et al's test case #5
!----------------------------
! Flow over isolated mountain
!----------------------------
      if(master) write(6,*) 'Generating IC for Flow over isolated mountain case'
      v0  = 20.
      h0  = 5960.*grav
      hs0 = 2000.*grav

      br  = pi/9.
      zc  = .5*pi - pi
      pc  = pi/6.

      do j=jfirst,jlast
         phi = - pi*0.5 + (j-1)*dp
         do i=1,im
            zam = dl*(i-1) - pi
            tmp = (zam-zc)**2 + (phi-pc)**2
            sr = min(br**2, tmp)
            phis(i,j) = hs0*(1.-sqrt(sr)/br)
             u(i,j,1)  =  v0*cose(j)
             v(i,j,1)  =  0.
! Total geopotential
          delp(i,j,1)  =  h0 - (ae*omega+0.5*v0)*v0*sinp(j)**2
         enddo
      enddo

     case(6)                  ! Williamson et al's test case #6
!-----------------------
! Rossby-Haurwitz wave#4
!-----------------------

      if(master) write(6,*) 'Generating IC for Rossby-Haurwitz wave#4 case'

      omg = 7.848E-6
      k = omg
      R = 4.
      h0 = 8.E3 * grav

      do j=jfirst, jlast
         do i=1,im
                zam  = DL*(i-1) - PI
            u(i,j,1) = AE*omg*cose(j) +        &
             AE*k*(cose(j)**(R-1))*(R*sine(j)**2-cose(j)**2)*cos(R*zam)
         enddo
      enddo

      do j=jfirst, jlast
         cosz = cosp(j)
         do i=1,im
            zam  = DL*(i-1) - PI
            zamc = zam - 0.5*DL
            v(i,j,1) = -AE*k*R*sinp(j)*sin(R*zamc)*cosz**(R-1)
            A = 0.5*omg*(2.*OMEGA+omg)*cosz**2 + 0.25*k*k*cosz**(R+R) *  &
                ((R+1)*cosz**2 + (2.*R*R-R-2.))-0.5*(k*R)**2*cosz**(R+R-2)
            B = 2.*(OMEGA+omg)*k / ((R+1)*(R+2)) *(cosz**R) *                  &
                ( (R*R+2.*R+2.) - ((R+1.)*cosz)**2 )
            C = 0.25*k*k* ((R+1)*cosz**2-(R+2)) * cosz**(2.*R)
           delp(i,j,1) = h0 + AE*AE*(A+B*cos(R*zam)+C*cos(2.*R*zam))
           phis(i,j) = 0.
         enddo
      enddo

     end select

! For all cases:

#ifdef SPMD
     call mp_send4d_ns(im, jm, 1, 1, jfirst, jlast, 1, 1, ig, ig, phis)
     call mp_recv4d_ns(im, jm, 1, 1, jfirst, jlast, 1, 1, ig, ig, phis)
#endif

!--------------------------------------
! Initialize the passive tracer with PV
!--------------------------------------
     call vort_d(im, jm, 1, jfirst, jlast, u, v, vort, ng_s, ng_d)

     do j=jfirst, jlast
        do i=1,im
! Substract out surface geopotential to get Layer thickness
           delp(i,j,1)   =  delp(i,j,1) - phis(i,j)
              q(i,j,1,1) = (f0(j) + vort(i,j)) / delp(i,j,1)
             pt(i,j,1)   =  1.
        enddo
     enddo

  end subroutine init_sw_ic
#endif
