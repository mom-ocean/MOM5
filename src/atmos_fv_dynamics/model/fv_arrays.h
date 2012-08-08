! $Id: fv_arrays.h,v 17.0 2009/07/21 02:53:09 fms Exp $
! this include file associates the shared arrays
! since this file contains both use statements and declarations it
! should be inserted in subroutines after all their own use
! statements and before any other declarations
  use fv_arrays_mod, only: isp, iep, jsp, jep, ksp, kep
  use fv_arrays_mod, only: fv_array_check
  use fv_arrays_mod, only: nlev, ncnst

#ifdef use_shared_pointers
! get the pointers and sizes from fv_arrays_mod
  use fv_arrays_mod, only: ptr_u, ptr_v, ptr_delp, ptr_pt, ptr_q, &
       ptr_u_phys, ptr_v_phys, ptr_t_phys, ptr_q_phys, &
       ptr_phis, ptr_ps, ptr_omga, ptr_pkz, &
       ptr_pk, ptr_pe, ptr_peln, ptr_pesouth, ptr_ua, ptr_va, ptr_ps_bp, &
       ptr_u_srf, ptr_v_srf, ptr_u_dt, ptr_v_dt, ptr_t_dt, ptr_q_dt
  use fv_arrays_mod, only: is,ie,js,je, isd,ied,jsd,jed, isg,ieg,jsg,jeg
#  ifdef MARS_GCM
  use fv_arrays_mod, only: ptr_delp_dt, ptr_mars_sfc_budg
#  endif MARS_GCM

#else
!get arrays from fv_arrays_mod
  use fv_arrays_mod, only: u, v, delp, pt, q, u_phys, v_phys, &
       t_phys, q_phys, phis, ps, omga, pkz, &
       pk, pe, peln, pesouth, ua, va, ps_bp, &
       u_srf, v_srf, u_dt, v_dt, t_dt, q_dt

#  ifdef MARS_GCM
  use fv_arrays_mod, only: delp_dt, mars_sfc_budg
#  endif MARS_GCM

#endif  use_shared_pointers

  implicit none

