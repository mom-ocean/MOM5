! Parallelized utility routine for computing
! max of an input scalar
!
      subroutine getmax( psize, pmax )

#ifdef SPMD
      use mod_comm, only : mp_reduce_max
#endif

      implicit none

      integer psize
      real pmax(psize)

#ifdef SPMD
      call mp_reduce_max(psize, pmax)
#endif

      return
      end
