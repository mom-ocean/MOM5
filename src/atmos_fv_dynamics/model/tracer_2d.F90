!BOP
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!
! !INTERFACE:
#ifndef USE_LIMA
subroutine tracer_2d(q, nq, dp1, cx, cy, mfx, mfy, iord, jord,   &
                     ng, cose, cosp, acosp, acap, rcap, tra_step,  &
                     fill, im, jm, km, jfirst, jlast, flx, dt, frac)
#else
subroutine tracer_2d_lima(q, nq, dp1, cx, cy, mfx, mfy, iord, jord,   &
                     ng, cose, cosp, acosp, acap, rcap, fill,    &
                     im, jm, km, jfirst, jlast, va, flx, dt, frac)
#endif
 
! !USES:
use fv_pack, only : beglon, endlon, beglat, endlat
#ifndef USE_LIMA
   use tp_core,          only: split_trac
   use diag_manager_mod, only: send_data
   use fv_diagnostics,   only: id_divg, fv_time
#else
   use shr_kind_mod,     only: r8 => shr_kind_r8
   use tp_core,          only: tp2c
   use fill_module,      only: fillxy
   use diag_manager_mod, only: send_data
   use fv_diagnostics,   only: id_divg, fv_time
#endif

#ifdef SPMD
   use mod_comm,         only: mp_send4d_ns, mp_recv4d_ns, mp_send3d_2,   &
                               mp_recv3d_2,  mp_reduce_max, gid
#endif
   use fv_arrays_mod, only: fv_array_check, fv_array_sync, jsp, jep, ksp, kep
   use mpp_mod, only: mpp_sync
   use fv_arrays_mod, only: fv_stack_push

   implicit none

! !INPUT PARAMETERS:

   integer, intent(in):: im, jm, km
   integer, intent(in):: jfirst, jlast
   integer, intent(in):: ng          ! Max number of ghost latitudes
   integer, intent(in):: nq
   integer, intent(in):: iord,  jord
#ifndef USE_LIMA
   integer, intent(inout):: tra_step
#endif

   real, intent(in)::  cose(jm)
   real, intent(in)::  cosp(jm)
   real, intent(in):: acosp(jm)
   real, intent(in):: acap, rcap
   real, intent(in):: dt

   logical, intent(in):: fill

! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout)::  cx(im,jfirst-ng:jlast+ng,km)
   real, intent(inout):: mfx(im,jfirst:jlast,km)

   real, intent(inout)::  cy(im,jfirst:jlast+1,km)
   real, intent(inout):: mfy(im,jfirst:jlast+1,km)

   real, intent(inout):: dp1(im,jfirst:jlast,km)
   real, intent(inout)::   q(im,jfirst-ng:jlast+ng,km,nq)
   real, intent(out):: frac

! Input work arrays (useless output)
   real, intent(out):: flx(im,jfirst:jlast,km)

#ifdef USE_LIMA 
   real, intent(out)::  va(im,jfirst:jlast,km)
#endif
! !DESCRIPTION:
!
!  Perform large-time-step tracer transport using accumulated Courant
!  numbers (cx, cy) and the mass fluxes (mfx, mfy) within the Lagrangian
!  layers.  This routine is 100\% parallel in the vertical direction
!  (with SMP).  Merdional Courant number will be further split, if
!  necessary, to ensure stability.  Cy <= 1 away from poles; Cy $\le$
!  1/2 at the latitudes closest to the poles.
!
! !CALLED FROM:
!     fvcore
!
!EOP
!---------------------------------------------------------------------
!BOC
      real, parameter:: esl = 1.e-10

! Local variables:
      real  cymax(km)
      real  cy_global
      real  cmax
      real  sum1, sum2
      real  rdt
      real  dp2(im,jfirst:jlast,km)


! Local 2d arrays
#ifdef USE_LIMA
      real(r8) a2(im,jfirst:jlast)
#endif
      real  fx(im,jfirst:jlast)
      real  fy(im,jfirst:jlast+1)
      real  dh(im,jfirst:jlast)

      logical ffsl(jm,km)
      logical used
      integer i, j, k
      integer it, iq, nsplt
      integer js2g0, js2gd, jn2g0,jn2gd


#ifdef use_shared_pointers
      pointer( p_dp2, dp2 )
      pointer( p_cymax, cymax )
      call fv_stack_push( p_dp2, im*km*(jlast-jfirst+1) )
      call fv_stack_push( p_cymax, km )
#endif
      js2g0 = max(2,jfirst)
      jn2g0 = min(jm-1,jlast)
      js2gd = max(2,jfirst-ng)     ! NG latitudes on S (starting at 2)
      jn2gd = min(jm-1,jlast+ng)   ! NG latitudes on S (ending at jm-1)

      call fv_array_check( LOC(cx) )
      call fv_array_check( LOC(mfx) )
      call fv_array_check( LOC(cy) )
      call fv_array_check( LOC(mfy) )
      call fv_array_check( LOC(dp1) )
      call fv_array_check( LOC(q) )
      call fv_array_check( LOC(flx) )
#ifdef USE_LIMA 
      call fv_array_check( LOC(va) )
#endif

!$omp parallel do private(i,j,k,cmax)
      do k=ksp,kep !1,km
         cymax(k) = 0.
         do j=js2g0,jlast
            cmax = 0.
            do i=1,im
               cmax = max( abs(cy(i,j,k)), cmax)
            enddo
            if( j==2 .or. j==jm ) then
! Do not allow N-S CFL number at poles to be greater than 0.5
! The limit is 1.0 elsewhere
                cymax(k) = max(cymax(k), 2.*cmax)
            else
                cymax(k) = max(cymax(k), cmax)
            endif
         enddo
      enddo
      call fv_array_sync()

#ifdef SPMD
#ifndef USE_LIMA
      if ( tra_step/=1 )  &
           call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng, ng, cx)
#else
      call mp_recv4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng, ng, cx)
#endif
      call mp_recv3d_2(gid+1,  im,  jm,      km,                       &
           1,  im,  jfirst,  jlast+1, 1, km,           &
           1,  im,  jlast+1, jlast+1, 1, km, cy, mfy )
      call mp_reduce_max(km, cymax)
! Send data on the way for ghosting the first tracer
      call mp_send4d_ns( im, jm, km, 1, jfirst, jlast, 1, km,       &
                         ng, ng, q(1,jfirst-ng,1,1) )
#endif
      call fv_array_sync()

! Start FMS diagnostics ----------

      if ( id_divg > 0 ) then
!------------------------------------------------------------
! (cx, cy) now ghosted; Compute mean divergence on the C grid
!------------------------------------------------------------
          rdt = 1./ dt
!$omp parallel do private(i,j,k, sum1)
          do k=ksp,kep !1,km
             do j=js2g0,jn2g0
                do i=1,im-1
                   flx(i,j,k) = cx(i+1,j,k) - cx(i,j,k) + (cy(i,j+1,k)*cose(j+1) -  &
                        cy(i,j,k)*cose(j))*acosp(j) 
                enddo
!              i=im
                flx(im,j,k) = cx(1,j,k) - cx(im,j,k) + (cy(im,j+1,k)*cose(j+1) -    &
                     cy(im,j,k)*cose(j))*acosp(j) 
             enddo
! Poles:
             if ( jfirst == 1 ) then
                 sum1 = 0.
                 do i=1,im
                    sum1 = sum1 + cy(i,1,k)
                 enddo
                 sum1 =-sum1*rcap*cose(2)
                 do i=1,im
                    flx(i,1,k) = sum1
                 enddo
             endif
             if ( jlast == jm ) then
                 sum1 = 0.
                 do i=1,im
                    sum1 = sum1 + cy(i,jm,k)
                 enddo
                 sum1 =-sum1*rcap*cose(jm)
                 do i=1,im
                    flx(i,jm,k) = sum1
                 enddo
             endif
             do j=jfirst, jlast
                do i=1,im
                   flx(i,j,k) = flx(i,j,k) * rdt
                enddo
             enddo
          enddo

          call fv_array_sync()
          used = send_data( id_divg, flx(beglon:endlon,beglat:endlat,:), fv_time ) 
      endif
! End FMS diagnostics ----------

! find global max cymax
      cy_global = cymax(1)
      if ( km /= 1 ) then                ! if NOT shallow water test case
         do k=2,km
            cy_global = max(cymax(k), cy_global)
         enddo
      endif

      nsplt = int(1. + cy_global)
      frac  = 1. / float(nsplt)

!     if ( gid == 0 .and. nsplt > 1 )  write(6,*) 'Tracer_2d_split=', nsplt, cy_global

      
      call fv_array_sync() ! don't modify flx before send_data completed
!$omp parallel do private(i,j,k)
      do 4000 k=ksp,kep !1,km

         if( nsplt /= 1 ) then
!!!     do j=2,jm-1
             do j=js2gd,jn2gd                  
                do i=1,im
                   cx(i,j,k) =  cx(i,j,k) * frac      ! cx ghosted on N*ng S*ng
                enddo
             enddo

!!!     do j=2,jm-1
             do j=js2g0,jn2g0
                do i=1,im
                   mfx(i,j,k) = mfx(i,j,k) * frac
                enddo
             enddo

!!!     do j=2,jm
             do j=js2g0,min(jm,jlast+1)
                do i=1,im
                   cy(i,j,k) =  cy(i,j,k) * frac    ! cy ghosted on N
                   mfy(i,j,k) = mfy(i,j,k) * frac    ! mfy ghosted on N
                enddo
             enddo
         endif

#ifdef USE_LIMA  
!!!     do j=2,jm-1
         do j=js2g0,jn2g0
            do i=1,im
               if(cy(i,j,k)*cy(i,j+1,k) > 0.) then
                   if( cy(i,j,k) > 0.) then
                       va(i,j,k) = cy(i,j,k)
                   else
                       va(i,j,k) = cy(i,j+1,k)      ! cy ghosted on N
                   endif
               else
                   va(i,j,k) = 0.
               endif
            enddo
         enddo
#endif

! Check if FFSL extension is needed.

         do 2222 j=js2gd,jn2gd             ! flux needed on N*ng S*ng
            ffsl(j,k) = .false.
            do i=1,im
               if(abs(cx(i,j,k)) > 1.) then  ! cx ghosted on N*ng S*ng
                   ffsl(j,k) = .true.
                   go to 2222
               endif
            enddo
2222     continue


! Scale E-W mass fluxes
         do j=js2g0,jn2g0
            if( ffsl(j,k) ) then
                do i=1,im
                   flx(i,j,k) = mfx(i,j,k)/sign(max(abs(cx(i,j,k)), esl),cx(i,j,k))
                enddo
            else
                do i=1,im
                   flx(i,j,k) = mfx(i,j,k)
                enddo
            endif
         enddo
4000  continue
!        call fv_array_sync()
        
      do 6000 it=1, nsplt

!$omp parallel do private(i, j, k, sum1, sum2)
         do k=ksp,kep !1,km
            do j=js2g0,jn2g0
               do i=1,im-1
                  dp2(i,j,k) =  dp1(i,j,k) + mfx(i,j,k) - mfx(i+1,j,k) +  &
                       (mfy(i,j,k) - mfy(i,j+1,k)) * acosp(j)
               enddo
               dp2(im,j,k) = dp1(im,j,k) + mfx(im,j,k) - mfx(1,j,k) +     &
                    (mfy(im,j,k) - mfy(im,j+1,k)) * acosp(j)
            enddo

! Poles
            if ( jfirst == 1  ) then
                sum1 = 0.
                do i=1,im
                   sum1 = sum1 + mfy(i,2,k)
                enddo

                sum1 = - sum1 * rcap
                do i=1,im
                   dp2(i,1,k) = dp1(i, 1,k) +  sum1
                enddo
            endif

            if ( jlast  == jm ) then
                sum2 = 0.
                do i=1,im
                   sum2 = sum2 + mfy(i,jm,k)
                enddo

                sum2 = sum2 * rcap
                do i=1,im
                   dp2(i,jm,k) = dp1(i,jm,k) +  sum2
                enddo
            endif
         enddo
         call fv_array_sync()

         do 5000 iq=1,nq

#ifdef SPMD
!            if( iq.EQ.1 ) &
!                 call mp_send4d_ns( im, jm, km, 1, jfirst, jlast, 1, km, &
!                 ng, ng, q(1,jfirst-ng,1,1) )
!------------------------------------
! Receive data for the current tracer
!------------------------------------
            call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast, 1, km,    &
                 ng, ng, q(1,jfirst-ng,1,iq))
!-------------------------------
! Send the next tracer if needed
!-------------------------------
            if ( iq /= nq )  then
                call mp_send4d_ns( im, jm, km, 1, jfirst, jlast, 1, km, &
                     ng, ng, q(1,jfirst-ng,1, iq+1) )
            else
!Last tracer is reached for this sub-step; send the next (first) tracer
                if( it /= nsplt .and. nq /=1 )                            &
                     call mp_send4d_ns( im, jm, km, 1, jfirst, jlast, 1, km, &
                     ng, ng, q(1,jfirst-ng,1, 1) )
            endif
            call fv_array_sync()
#endif

!$omp parallel do private(i, j, k, fx, fy, dh)
            do k=ksp,kep !1,km
#ifndef USE_LIMA
               call split_trac(dp1(1,jfirst,k), dp2(1,jfirst,k), dh,     &
                    q(1,jfirst-ng,k,iq), cx(1,jfirst-ng,k),   &
                    cy(1,jfirst,k), im, jm, iord, jord, ng,   &
                    fx, fy, ffsl(1,k), tra_step, fill,   &
                    acap, acosp, flx(1,jfirst,k), mfy(1,jfirst,k),  & 
                    cosp, 1, jfirst, jlast )
#else
               call tp2c(a2, va(1,jfirst,k),     q(1,jfirst-ng,k,iq),  &
                    cx(1,jfirst-ng,k) , cy(1,jfirst,k),       &
                    im, jm,  iord,       jord,   ng,              &
                    fx, fy,  ffsl(1,k),  rcap,   acosp,           &
                    flx(1,jfirst,k),     mfy(1,jfirst,k),         &
                    cosp, 1, jfirst, jlast )

               do j=jfirst,jlast
                  do i=1,im
                     a2(i,j) = q(i,j,k,iq)*dp1(i,j,k) + a2(i,j)
                  enddo
               enddo

               if (fill) call fillxy ( a2,  im, jm, jfirst, jlast,  &
                    acap, cosp, acosp )

               do j=jfirst,jlast
                  do i=1,im
                     q(i,j,k,iq) = a2(i,j) / dp2(i,j,k)
                  enddo
               enddo

#endif
            enddo
            call fv_array_sync()
5000     end do !continue

         if( it /= nsplt ) then           ! not the last

#ifdef SPMD
!------------------
! Handle nq=1 case:
!------------------
             if( nq==1 )       &
                  call mp_send4d_ns(im, jm, km, 1, jfirst, jlast, 1, km, ng, ng, q)
#endif

!$omp parallel do private(i, j, k)
             do k=ksp,kep !1,km
                do j=jfirst,jlast
                   do i=1,im
                      dp1(i,j,k) = dp2(i,j,k)
                   enddo
                enddo
             enddo
             call fv_array_sync()
         endif
6000  end do !continue

#ifndef USE_LIMA
      if (tra_step == 1) then
             tra_step = 2
      elseif(tra_step == 2) then
             tra_step = 1
      endif
#endif

#ifndef USE_LIMA
end subroutine tracer_2d
#else               
end subroutine tracer_2d_lima
#endif
