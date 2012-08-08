module fill_module
!-----------------------------------------------------------------------
!BOP
!
     implicit none
! !MODULE: fill_module --- utilities for filling in nagative values
!
! !PUBLIC MEMBER FUNCTIONS:
      public filew, fillxy, fillz, pfix

!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!   99.03.01   Lin        Creation
!   01.02.14   Lin        Routines coalesced into this module
!   01.03.26   Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fillxy --- Fill from east, west, north and south neighbors
!
! !INTERFACE: 
 subroutine fillxy(q, im, jm, jfirst, jlast, acap, cosp, acosp)

 integer, intent(in):: im           ! Longitudes
 integer, intent(in):: jm           ! Total latitudes
 integer, intent(in):: jfirst       ! Starting latitude
 integer, intent(in):: jlast        ! Finishing latitude

 real , intent(in):: acap 
 real , intent(in):: cosp(jm)
 real , intent(in):: acosp(jm) 
 
 real , intent(inout):: q(im,jfirst:jlast) ! Field to adjust

!
! LOCAL VARIABLES:
  integer ipx, i,j
  real :: tiny = 1.e-40

  call filew(q,im,jm,jfirst,jlast,acap,ipx,tiny, cosp(2), cosp(3))

!  if(ipx.ne.0) then
! Fill zonally
!     call xfix(q,IM,JM,tiny,qt)
!  endif

!EOC
 end subroutine fillxy

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: filew --- Fill from east and west neighbors; essentially
!                      performing local flux adjustment
!
! !INTERFACE: 
 subroutine filew(q, im, jm, jfirst, jlast, acap, ipx, tiny, cosp2, cosp3)

! !INPUT PARAMETERS:
 integer im                  ! Longitudes
 integer jm                  ! Total latitudes
 integer jfirst              ! Starting latitude
 integer jlast               ! Finishing latitude

 real  tiny               ! A small number to pump up value
 real  acap               ! 1/(polar cap area)
 real  cosp2              ! cosine(lat) at j=2
 real  cosp3              ! cosine(lat) at j=3

! !INPUT/OUTPUT PARAMETERS:
 real  q(im,jfirst:jlast) ! Field to adjust

! !OUTPUT PARAMETERS:
 integer ipx                 ! Flag:  0 if Q not change, 1 if changed

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !REVISION HISTORY:
!   99.10.01   Lin        Creation
!   07.30.01   Lin        Improvement
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
 real  d0, d1, d2
 real  qtmp(jfirst:jlast,im)

 integer i, j, jm1, ip2
 integer j1, j2

 jm1 = jm-1
 j1 = max( jfirst,  2 )
 j2 = min( jlast, jm1 )

  ipx = 0
  do j=jfirst,jlast
     do i=1,im
        if( q(i,j) < 0. ) then
            ipx = 1
            go to 444
        endif
     enddo
  enddo

  if ( ipx == 0 ) return

! Copy & swap direction for vectorization.
444   do i=1,im
         do j=j1,j2
            qtmp(j,i) = q(i,j)
         enddo
      enddo

 
  do i=2,im-1
     do j=j1,j2
        if(qtmp(j,i) < 0.) then
! west
           d0 = max(0.,qtmp(j,i-1))
           d1 = min(-qtmp(j,i),d0)
           qtmp(j,i-1) = qtmp(j,i-1) - d1
           qtmp(j,i) = qtmp(j,i) + d1
! east
           d0 = max(0.,qtmp(j,i+1))
           d2 = min(-qtmp(j,i),d0)
           qtmp(j,i+1) = qtmp(j,i+1) - d2
           qtmp(j,i) = qtmp(j,i) + d2 + tiny
        endif
    enddo
  enddo
 
     i=1
  do j=j1,j2
     if(qtmp(j,i) < 0.) then
! west
        d0 = max(0.,qtmp(j,im))
        d1 = min(-qtmp(j,i),d0)
        qtmp(j,im) = qtmp(j,im) - d1
        qtmp(j,i) = qtmp(j,i) + d1
! east
        d0 = max(0.,qtmp(j,i+1))
        d2 = min(-qtmp(j,i),d0)
        qtmp(j,i+1) = qtmp(j,i+1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
  enddo

     i=im
  do j=j1,j2
     if(qtmp(j,i) < 0.) then
! west
        d0 = max(0.,qtmp(j,i-1))
        d1 = min(-qtmp(j,i),d0)
        qtmp(j,i-1) = qtmp(j,i-1) - d1
        qtmp(j,i) = qtmp(j,i) + d1
! east
        d0 = max(0.,qtmp(j,1))
        d2 = min(-qtmp(j,i),d0)
        qtmp(j,1) = qtmp(j,1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
     endif
  enddo

!-----------
! Final pass
!-----------
    do i=1,im-1
       do j=j1,j2
          if (qtmp(j,i) < 0. ) then
! Take mass from east (essentially adjusting fx(i+1,j))
              qtmp(j,i+1) = qtmp(j,i+1) + qtmp(j,i)
              qtmp(j,i) = 0.
          endif
       enddo
    enddo

    do i=im,2,-1
       do j=j1,j2
          if (qtmp(j,i) < 0. ) then
! Take mass from west (essentially adjusting fx(i,j))
              qtmp(j,i-1) = qtmp(j,i-1) + qtmp(j,i)
              qtmp(j,i) = 0.
          endif
       enddo
    enddo

    do j=j1,j2
       do i=1,im
          q(i,j) = qtmp(j,i)
       enddo
    enddo
 
! Check Poles.

 if ( jfirst == 1 ) then
      ip2 = 0
      if(q(1,1) < 0.) then
         call pfix(q(1,2),q(1,1),im,ipx,acap,cosp2)
      else
! Check the latitude next to the SP
         do i=1,im
            if(q(i,2) < 0.) then
               ip2 = 1
               go to 322
            endif
         enddo
322      continue
         if(ip2.ne.0) call pfix(q(1,2),q(1,1),im,ipx,acap,cosp2)
      endif
 endif
 
 if ( jlast == jm ) then
      if(q(1,jm) < 0.) then
         call pfix(q(1,jm1),q(1,jm),im,ipx,acap,cosp2)
      else
! Check the latitude next to the NP
              ip2 = 0
         do i=1,im
            if(q(i,jm1) < 0.) then
               ip2 = 1
               go to 323
            endif
         enddo
323      continue
         if(ip2.ne.0) call pfix(q(1,jm1),q(1,jm),im,ipx,acap,cosp2)
      endif
 endif

!EOC
 end subroutine filew
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fillz --- Fill from neighbors below and above
!
! !INTERFACE: 
 subroutine fillz(im, km, nq, q, dp)

! !INPUT PARAMETERS:
   integer,  intent(in):: im                ! No. of longitudes
   integer,  intent(in):: km                ! No. of levels
   integer,  intent(in):: nq                ! Total number of tracers

   real , intent(in)::  dp(im,km)       ! pressure thickness
! !INPUT/OUTPUT PARAMETERS:
   real , intent(inout) :: q(im,km,nq)   ! tracer mixing ratio

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !BUGS:
!   Currently this routine only performs the east-west fill algorithm.
!   This is because the N-S fill is very hard to do in a reproducible
!   fashion when the problem is decomposed by latitudes.
!
! !REVISION HISTORY:
!   00.04.01   Lin        Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer i, k, ic
   real  qup, qly, dup

   do ic=1,nq
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0.) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      do k=2,km-1
         do i=1,im
         if( q(i,k,ic) < 0. ) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( 0.75*qly, qup )        !borrow no more than 75% from top
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
! Borrow from below: q(i,k,ic) is still negative at this stage
             q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1) 
             q(i,k  ,ic) = 0.
          endif
          enddo
      enddo
 
! Bottom layer
      k = km
      do i=1,im
         if( q(i,k,ic) < 0.) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( qly, qup )
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,ic) = 0.
          endif
      enddo
   enddo
end subroutine fillz
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
 subroutine pfix(q, qp, im, ipx, acap, cosp2)

 integer im                  ! Longitudes
 real  acap               ! ???
 real  cosp2              ! ???

 real  q(im)              ! Latitude-level field to adjust
 real  qp(im)             ! Second latitude-level field to adjust (usually pole)

! !OUTPUT PARAMETERS:
 integer ipx                 ! Flag:  0 if Q not change, 1 if changed


! !LOCAL VARIABLES:
 integer i
 real  summ, sump, pmean
 
   summ = 0.
   sump = 0.
   do i=1,im
     summ = summ + q(i)
     sump = sump + qp(i)
   enddo
 
   sump = sump/im
   pmean = (sump*acap + summ*cosp2) / (acap + cosp2*im)
 
   do i=1,im
      q(i) = pmean
      qp(i) = pmean
   enddo
 
   if( qp(1) < 0. ) then
      ipx = 1
   endif

!EOC
 end subroutine pfix
!-----------------------------------------------------------------------

end module fill_module
