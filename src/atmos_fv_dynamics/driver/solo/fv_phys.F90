module fv_phys_mod

  public :: fv_phys

contains

  subroutine fv_phys( Time , ndt )

    use time_manager_mod, only: time_type
    use hs_forcing_mod,   only: hs_forcing
    use constants_mod,    only: grav, kappa, rdgas
    use fv_pack, only: nt_phys, beglat, endlat, nlon, nlev, rlon, rlat
    use shr_kind_mod,     only : r8 => shr_kind_r8
    use update_fv_phys_mod, only: update_fv_phys
#include <fv_arrays.h>
!  implicit none

    type (time_type), intent(in):: Time
    integer,          intent(in):: ndt         ! physics time step (sec)

! Local arrays:
!    real(r8)  dt_u(nlon, beglat:endlat, nlev)
!    real(r8)  dt_v(nlon, beglat:endlat, nlev)
!    real(r8)  dt_t(nlon, beglat:endlat, nlev)
    real(r8)  dt_q(nlon, beglat:endlat, nlev, ncnst)

    real(r8) p_full(nlon, beglat:beglat, nlev)
    real(r8) p_half(nlon, beglat:beglat, nlev+1)

    real(r8) dt, dt5
    integer i, j, k, m
    integer ij, nx, tsiz, isiz
#ifndef use_shared_pointers
    integer is, ie
#endif
    logical rayf
    logical strat
#include <fv_point.inc>

    dt = ndt

    dt5 = 0.5 * dt

#ifdef HS_TSM

    rayf  = .false.
    strat = .false.

    call hswf(nlon,  mlat,  nlev,  beglat,  endlat,       &
         u, v, pt, pe, delp, peln, pkz, ndt,         &
         kappa, grav, rdgas,  strat, rayf, master,   &
         sinp, cosp, sine, cose,  ng_s, ng_d)
#else
!    nt_phys = 1

!   isiz = window(1)
    isiz = nlon/4
    nx = nlon/isiz
    tsiz = nx * ( endlat - beglat + 1 )         ! Total loop length


! Calling phys one latitude at a time using a big/fat OpenMP loop
! For cache performance, this could be changed to finer decomposition if needed

!$omp parallel do private(ij,i,j,k,m,is,ie,p_full,p_half)
    do ij=1,tsiz              

       j  = beglat + (ij-1) / nx

       is = 1 + isiz * mod(ij-1, nx)
       ie = is + isiz - 1

       do k=1,nlev
          do i=is,ie
             u_dt(i,j,k) = 0.
             v_dt(i,j,k) = 0.
             t_dt(i,j,k) = 0.
          enddo
       enddo

       do m=1,ncnst
          do k=1,nlev
             do i=is,ie
                q_dt(i,j,k,m) = 0.
             enddo
          enddo
       enddo


       do k=1,nlev+1
          do i=is,ie
             p_half(i,beglat,k) = pe(i,k,j)
          enddo
       enddo

       do k=1,nlev
          do i=is,ie
             p_full(i,beglat,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo

       call hs_forcing(is, ie, j-beglat+1, j-beglat+1, dt, Time,   &
            rlon(is:ie,j:j), rlat(is:ie,j:j),                      &
            p_half(is:ie,beglat:beglat,:),                         &
            p_full(is:ie,beglat:beglat,:),                         &
            ua(is:ie,j:j,1:nlev), va(is:ie,j:j,1:nlev),            &
            pt(is:ie,j:j,1:nlev),  q(is:ie,j:j,1:nlev,1:ncnst),    &
            ua(is:ie,j:j,1:nlev), va(is:ie,j:j,1:nlev),            &
            pt(is:ie,j:j,1:nlev),  q(is:ie,j:j,1:nlev,1:nt_phys),  &
            u_dt(is:ie,j:j,:), v_dt(is:ie,j:j,:),                  &
            t_dt(is:ie,j:j,:), q_dt(is:ie,j:j,:,1:nt_phys) )
    enddo

    call update_fv_phys ( dt, nt_phys, .false., .false., Time )
#endif

  end subroutine fv_phys

end module fv_phys_mod
