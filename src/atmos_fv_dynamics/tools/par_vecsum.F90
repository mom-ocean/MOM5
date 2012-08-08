 subroutine par_vecsum(jm, jfirst, jlast, InVector, te0)

! !ROUTINE: par_vecsum --- Calculate vector sum bit-wise consistently

#ifdef SPMD
      use mod_comm, only : mp_sum1d
#ifdef use_shared_pointers
      use fv_arrays_mod, only: fv_thread_bcast
#endif
#endif
      implicit none
! Input:
      integer jm                   ! global latitudes
      integer jfirst               ! first latitude on this PE
      integer jlast                ! last latitude on this PE
      real, intent(in):: InVector(jfirst:jlast)  ! input vector to be summed

! OUTPUT:
      real, intent(out):: te0   ! sum of all vector entries
! Local
      integer j

#ifdef SPMD 
      call mp_sum1d(jm, jfirst, jlast, InVector, te0)
#ifdef use_shared_pointers
      call fv_thread_bcast(te0)
#endif
#else
      te0 = 0.0
      do j=1,jm
        te0 = te0 + InVector(j)
      enddo
#endif

 end subroutine par_vecsum
