module pv_module

  use fv_pack,      only: u_ghosted, cose, rcap, one_by_dy,   &
                          one_by_dx, f_d, acosp

#ifdef SPMD
  use mod_comm, only: gid, mp_send3d, mp_recv3d
#endif

  use fv_arrays_mod, only: fv_array_check, fv_array_sync, jsp, jep, ksp, kep
  implicit none

public vort_d,  pv_entropy

contains

  subroutine vort_d(im, jm, km, jfirst, jlast, u, v, vort, ng_s, ng_d)

! Computes relative vorticity on the D grid using circulation theorem

! !INPUT PARAMETERS:
    integer im, jm, km                   ! Horizontal dimensions
    integer jfirst, jlast                ! Latitude strip
    integer ng_s, ng_d

    real:: u(im,jfirst-ng_d:jlast+ng_s,km) 
    real:: v(im,jfirst-ng_d:jlast+ng_d,km) 

! !OUTPUT PARAMETERS:
    real, intent(out):: vort(im,jfirst:jlast,km)

! local
    real fx(im,jfirst-ng_d:jlast+ng_d)
! Geometric arrays
    real cy(jfirst:jlast)

    integer i, j, k,  js2g0, jn2g0
    real c1, c2

    call fv_array_check( LOC(u) )
    call fv_array_check( LOC(v) )
    call fv_array_check( LOC(vort) )

#ifdef SPMD
    call fv_array_sync()
    if ( .not. u_ghosted ) then
        call mp_send3d(gid-1, gid+1, im, jm, km, &
             1, im, jfirst-ng_d, jlast+ng_s, 1, km, &
             1, im, jfirst, jfirst, 1, km, u)
        call mp_recv3d(gid+1, im, jm, km, &
             1, im, jfirst-ng_d, jlast+ng_s, 1, km, &
             1, im, jlast+1, jlast+1, 1, km, u)
        u_ghosted = .true.
    endif
    call fv_array_sync()
#endif

    js2g0 = max(2,jfirst)
    jn2g0 = min(jm-1,jlast)

!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k, fx, c1, c2, cy)

    do k=ksp,kep

       do j=js2g0,min(jlast+1,jm)
          do i=1,im
             fx(i,j) = u(i,j,k)*cose(j)
          enddo
       enddo

       do j=js2g0,jn2g0
          cy(j) =  one_by_dy * acosp(j)
          do i=1,im-1
             vort(i,j,k) = ( fx(i,j) - fx(i,j+1) ) * cy(j) +      &
                  ( v(i+1,j,k) - v(i,j,k) ) * one_by_dx(j)
          enddo
          i = im
          vort(i,j,k) = ( fx(i,j) - fx(i,j+1) ) * cy(j)    +      &
               ( v(1,j,k) - v(i,j,k) ) * one_by_dx(j)
       enddo

! Vort at poles computed by circulation theorem

       if ( jfirst == 1 ) then
           c1 = -SUM(fx(1:im,2))*one_by_dy*rcap
           do i=1,im
              vort(i,1,k) = c1
           enddo
       endif

       if ( jlast == jm )  then
           c2 = SUM(fx(1:im,jm))*one_by_dy*rcap
           do i=1,im
              vort(i,jm,k) = c2
           enddo
       endif
    enddo
    call fv_array_sync()

  end subroutine vort_d


  subroutine pv_entropy(im, jm, km, jfirst, jlast, vort, pt, pkz, delp, ng, grav)
#ifdef use_shared_pointers
   use fv_arrays_mod, only: fv_stack_push
#endif

! !INPUT PARAMETERS:
   integer, intent(in):: im, jm, km               ! Horizontal dimensions
   integer, intent(in):: jfirst, jlast            ! Latitude strip
   integer, intent(in):: ng
   real, intent(in):: grav
   real, intent(in):: pt(im,jfirst-ng:jlast+ng,km) 
   real, intent(in):: pkz(im,jfirst:jlast,km) 
   real, intent(in):: delp(im,jfirst:jlast,km)

! vort is relative vorticity as input. Becomes PV on output
      real, intent(inout):: vort(im,jfirst:jlast,km) 

! !DESCRIPTION:
!        EPV = 1/r * (vort+f_d) * d(S)/dz; where S is a conservative scalar
!        r the fluid density, and S is chosen to be the entropy here: S = log(pt)
!        pt == potential temperature.
! Computation are performed assuming the input is on "x-y-z" Cartesian coordinate.
! The approximation made here is that the relative vorticity computed on constant
! z-surface is not that different from the hybrid sigma-p coordinate.
! See page 39, Pedlosky 1979: Geophysical Fluid Dynamics
!
! The follwoing simplified form is strictly correct only if vort is computed on 
! constant z surfaces. In addition hydrostatic approximation is made.
!        EPV = - GRAV * (vort+f_d) / del(p) * del(pt) / pt 
! where del() is the vertical difference operator.
!
! programmer: S.-J. Lin, sjl
!
!EOP
!---------------------------------------------------------------------
!BOC
      real w3d(im,jfirst:jlast,km) 
      real te(im,jm,km+1), t2(im,km), delp2(im,km)
      real te2(im,km+1)
      integer i, j, k
#ifdef use_shared_pointers
      pointer( p_w3d, w3d )
      pointer( p_te, te )
      call fv_stack_push( p_w3d, im*km*(jlast-jfirst+1) )
      call fv_stack_push( p_te,  im*jm*(km+1)           )
#endif
      call fv_array_check( LOC(pt) )
      call fv_array_check( LOC(pkz) )
      call fv_array_check( LOC(delp) )
      call fv_array_check( LOC(vort) )

#ifdef SW_DYN
        do j=jsp,jep !jfirst,jlast
          do i=1,im
#ifndef USE_LIMA
            vort(i,j,1) = grav * (vort(i,j,1)+f_d(i,j)) / delp(i,j,1)
#else
            vort(i,j,1) = grav * (vort(i,j,1)+f_d(j)) / delp(i,j,1)
#endif
          enddo
        enddo
        call fv_array_sync()
#else
! Compute PT at layer edges.

!$omp  parallel do                 &
!$omp  default(shared)             &
!$omp  private(i,j,k,t2,delp2,te2)

     do j=jsp,jep !jfirst,jlast

        do k=1,km
          do i=1,im
               t2(i,k) = pt(i,j,k) / pkz(i,j,k)
              w3d(i,j,k) = t2(i,k)
            delp2(i,k) = delp(i,j,k)
          enddo
        enddo

        call ppme(t2, te2, delp2, im, km)

        do k=1,km+1
           do i=1,im
              te(i,j,k) = te2(i,k)
           enddo
        enddo
     enddo
     call fv_array_sync()

!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k)

     do k=ksp,kep !1,km
        do j=jfirst,jlast
          do i=1,im
! Entropy is the thermodynamic variable in the following form
#ifndef USE_LIMA
            vort(i,j,k) = (vort(i,j,k)+f_d(i,j)) * ( te(i,j,k)-te(i,j,k+1) )  &
                          / ( w3d(i,j,k)*delp(i,j,k) )  * grav
#else
            vort(i,j,k) = (vort(i,j,k)+f_d(j)) * ( te(i,j,k)-te(i,j,k+1) )  &
                          / ( w3d(i,j,k)*delp(i,j,k) )  * grav
#endif
          enddo
        enddo
     enddo
     call fv_array_sync()
#endif

  end subroutine pv_entropy


  subroutine ppme(p,qe,delp,im,km)

  integer, intent(in):: im, km
  real, intent(in)::    p(im,km)
  real, intent(in):: delp(im,km)
  real, intent(out)::qe(im,km+1)

! local arrays.
      real dc(im,km),delq(im,km), a6(im,km)
      real c1, c2, c3, tmp, qmax, qmin
      real a1, a2, s1, s2, s3, s4, ss3, s32, s34, s42
      real a3, b2, sc, dm, d1, d2, f1, f2, f3, f4
      real qm, dq
      integer i, k, km1

      km1 = km - 1

      do 500 k=2,km
      do 500 i=1,im
500   a6(i,k) = delp(i,k-1) + delp(i,k)

      do 1000 k=1,km1
      do 1000 i=1,im
      delq(i,k) = p(i,k+1) - p(i,k)
1000  continue

      do 1220 k=2,km1
      do 1220 i=1,im
      c1 = (delp(i,k-1)+0.5*delp(i,k))/a6(i,k+1)
      c2 = (delp(i,k+1)+0.5*delp(i,k))/a6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /    &
                                    (a6(i,k)+delp(i,k+1))
      qmax = max(p(i,k-1),p(i,k),p(i,k+1)) - p(i,k)
      qmin = p(i,k) - min(p(i,k-1),p(i,k),p(i,k+1))
      dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
1220  continue

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

   do k=3,km1
      do i=1,im
         c1 = delq(i,k-1)*delp(i,k-1) / a6(i,k)
         a1 = a6(i,k-1) / (a6(i,k) + delp(i,k-1))
         a2 = a6(i,k+1) / (a6(i,k) + delp(i,k))
         qe(i,k) = p(i,k-1) + c1 + 2./(a6(i,k-1)+a6(i,k+1)) *        &
                   ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -         &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
   enddo

! three-cell parabolic subgrid distribution at model top

   do i=1,im
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2
      s1 = delp(i,1)
      s2 = delp(i,2) + s1
!
      s3 = delp(i,2) + delp(i,3)
      s4 = s3 + delp(i,4)
      ss3 =  s3 + s1
      s32 = s3*s3
      s42 = s4*s4
      s34 = s3*s4
! model top
      a3 = (delq(i,2) - delq(i,1)*s3/s2) / (s3*ss3)
!
      if(abs(a3) .gt. 1.e-14) then
         b2 =  delq(i,1)/s2 - a3*(s1+s2)
         sc = -b2/(3.*a3)
         if(sc .lt. 0. .or. sc .gt. s1) then
             qe(i,1) = p(i,1) - s1*(a3*s1 + b2)
         else
             qe(i,1) = p(i,1) - delq(i,1)*s1/s2
         endif
      else
! Linear
         qe(i,1) = p(i,1) - delq(i,1)*s1/s2
      endif
      dc(i,1) = p(i,1) - qe(i,1)
! compute coef. for the off-centered area preserving cubic poly.
      dm = delp(i,1) / (s34*ss3*(delp(i,2)+s3)*(s4+delp(i,1)))
      f1 = delp(i,2)*s34 / ( s2*ss3*(s4+delp(i,1)) )
      f2 = (delp(i,2)+s3) * (ss3*(delp(i,2)*s3+s34+delp(i,2)*s4)   &
            + s42*(delp(i,2)+s3+s32/s2))
      f3 = -delp(i,2)*( ss3*(s32*(s3+s4)/(s4-delp(i,2))            &
            + (delp(i,2)*s3+s34+delp(i,2)*s4))                     &
            + s42*(delp(i,2)+s3) )
      f4 = ss3*delp(i,2)*s32*(delp(i,2)+s3) / (s4-delp(i,2))
      qe(i,2) = f1*p(i,1)+(f2*p(i,2)+f3*p(i,3)+f4*p(i,4))*dm
   enddo

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
   do i=1,im
      d1 = delp(i,km)
      d2 = delp(i,km1)
      qm = (d2*p(i,km)+d1*p(i,km1)) / (d1+d2)
      dq = 2.*(p(i,km1)-p(i,km)) / (d1+d2)
      c1 = (qe(i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
      c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
      qe(i,km  ) = qm - c1*d1*d2*(d2+3.*d1)
      qe(i,km+1) = d1*(8.*c1*d1**2-c3) + qe(i,km)
   enddo

  end subroutine ppme

end module pv_module
