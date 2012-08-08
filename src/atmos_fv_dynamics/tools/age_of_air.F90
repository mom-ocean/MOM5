      subroutine age_of_air(im, jm, km, jfirst, jlast, ng, time, pfull, q)

      implicit none
      integer im
      integer jm
      integer km
      integer jfirst
      integer jlast  
      integer ng

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

      real, intent(in):: pfull(im,jfirst:jlast,km)   ! model full level
      real, intent(in):: time        ! accumulated time since init
      real, intent(inout):: q(im,jfirst-ng:jlast+ng,km)

! Local
      integer i, j, k
      real p_source      ! source level (pa)
      real ascale
      real tiny
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( ascale = 5.e-6 / 60. )

!$omp parallel do private(i, j, k)
      do k=1,km
        do j=jfirst, jlast
          do i=1,im
            if( time < tiny ) then
                q(i,j,k) = 0.
            elseif( pfull(i,j,k) >= p_source ) then
                q(i,j,k) = ascale * time
            endif
          enddo
        enddo          ! j-loop
      enddo             ! k-loop

      return
      end
