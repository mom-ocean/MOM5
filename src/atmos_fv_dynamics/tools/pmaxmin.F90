module pmaxmin_mod

  implicit none
  private
  public :: pmaxming, pmaxmin, prt_maxmin_local

contains
! Parallelized utility routines for computing/printing
! max/min of an input array
!

  subroutine pmaxmin( qname, a, pmin, pmax, im, jt, fac )

#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
    use mod_comm, only : gid, mp_reduce_max
#else
#define CPP_PRT_PREFIX
#endif
    use fv_arrays_mod, only: fv_array_check, fv_array_limits, fv_stack_push, fv_array_sync
    use mpp_mod, only: mpp_sync

    implicit none

    character*(*)  qname
    integer im, jt
    integer i, j
    real a(im,jt)

    real qmin(jt), qmax(jt)
    real pmax, pmin
    real fac                     ! multiplication factor
    real pm1(2)

    integer :: js, je
#ifdef use_shared_pointers
    pointer( p_qmin, qmin )
    pointer( p_qmax, qmax )
    call fv_stack_push( p_qmin, jt )
    call fv_stack_push( p_qmax, jt )
#endif

    call fv_array_check( LOC(a) )
    call fv_array_limits( 1, jt, js, je )
!$omp parallel do private(i, j, pmax, pmin)

    do j=js,je
       pmax = a(1,j)
       pmin = a(1,j)
       do i=2,im
          pmax = max(pmax, a(i,j))
          pmin = min(pmin, a(i,j))
       enddo
       qmax(j) = pmax
       qmin(j) = pmin
    enddo
    call fv_array_sync()
!
! Now find max/min of amax/amin
!
    pmax = qmax(1)
    pmin = qmin(1)
    do j=2,jt
       pmax = max(pmax, qmax(j))
       pmin = min(pmin, qmin(j))
    enddo


#if defined( SPMD )
    pm1(1) = pmax
    call mp_reduce_max(1, pm1)
    pmax=pm1(1)
    pm1(1) = -pmin
    call mp_reduce_max(1, pm1)
    pmin=-pm1(1)
#endif

    CPP_PRT_PREFIX write(6,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

  end subroutine pmaxmin

  subroutine pmaxming(qname, a, im, jm, km,     &
       jfirst, jlast, ng_d, ng_s, fac)

    use fv_arrays_mod, only: fv_array_check, fv_array_sync, jsp, jep
    use mpp_mod, only: mpp_sync
#ifdef use_shared_pointers
   use fv_arrays_mod, only: fv_stack_push
#endif
    implicit none

    character*(*)  qname
    integer im, jm, km, jfirst, jlast, ng_d, ng_s
    real a(im,jfirst-ng_d:jlast+ng_s,km)
    real tmp(im, jfirst:jlast,km)
    real pmax, pmin, fac
    integer i, j, k
#ifdef use_shared_pointers
    pointer( p_tmp, tmp )
    call fv_stack_push( p_tmp, im*km*(jlast-jfirst+1) )
#endif
    call fv_array_check( LOC(a) )

!$omp parallel do private(i, j, k)
!balaji: note, parallelizing j not k... km may be 1.
    do k=1, km
       do j=jsp,jep
          do i=1, im
             tmp(i,j,k)=a(i,j,k)
          enddo
       enddo
    enddo
    call fv_array_sync()

    call pmaxmin(qname, tmp, pmin, pmax, im*(jlast-jfirst+1),km, fac)

  end subroutine pmaxming


! Routine for debug domain decomposed data

  subroutine prt_maxmin_local(gid, qname, q, im, jm, km, j1, j2)

    implicit none
    character*(*)  qname
    integer im, jm, km, j1, j2, gid
    real q(im,j1:j2,km)
    integer i,j,k,jp
    real qmax, qmin

    jp = max(1,j1)

    qmax = q(1,jp,1)
    qmin = qmax

    do k=1,km
       do j=max(2,j1), min(jm-1,j2)
          do i=1,im
             qmax = max(qmax, q(i,j,k))
             qmin = min(qmin, q(i,j,k))
          enddo
       enddo
    enddo
    write(*,*) 'GID=', gid, qname, ': Max=', qmax, ' Min=', qmin

  end subroutine prt_maxmin_local

end module pmaxmin_mod
