! $Id: fv_point.inc,v 17.0 2009/07/21 02:53:14 fms Exp $
#ifdef use_shared_pointers
! declare pointee arrays
  real ::    u(isg:ieg, jsd:jed, nlev)
  real ::    v(isg:ieg, jsd:jed, nlev)
  real :: delp(isg:ieg, js:je,   nlev)
  real ::   pt(isg:ieg, jsd:jed, nlev)
  real ::    q(isg:ieg, jsd:jed, nlev, ncnst)
  real ::  u_phys(isg:ieg, js:je, nlev)
  real ::  v_phys(isg:ieg, js:je, nlev)
  real ::  t_phys(isg:ieg, js:je, nlev)
  real ::  q_phys(isg:ieg, js:je, nlev, ncnst)
  real :: phis(isg:ieg, js:je)
  real :: ps  (isg:ieg, js:je)
  real :: omga(isg:ieg, js:je,nlev)    
  real :: pkz (isg:ieg, js:je, nlev)
  real :: pk  (isg:ieg, js:je, nlev+1)
  real :: pe  (isg:ieg, nlev+1, js:je) 
  real :: peln(isg:ieg, nlev+1, js:je) 
  real :: pesouth(isg:ieg, nlev+1)
  real ::  ua(isg:ieg, js:je, nlev)
  real ::  va(isg:ieg, js:je, nlev)
  real :: ps_bp (isg:ieg, js:je)
  real ::  u_srf(isg:ieg, js:je)
  real ::  v_srf(isg:ieg, js:je)
  real :: u_dt(isg:ieg, js:je, nlev)
  real :: v_dt(isg:ieg, js:je, nlev)
  real :: t_dt(isg:ieg, js:je, nlev)
  real :: q_dt(isg:ieg, js:je, nlev, ncnst)
#ifdef MARS_GCM
  real :: delp_dt(isg:ieg, js:je, nlev)
  real :: mars_sfc_budg(isg:ieg, js:je, 5)
#endif MARS_GCM
  pointer( p_u, u )
  pointer( p_v, v )
  pointer( p_delp, delp )
  pointer( p_pt, pt )
  pointer( p_q, q )
  pointer( p_u_phys, u_phys )
  pointer( p_v_phys, v_phys )
  pointer( p_t_phys, t_phys )
  pointer( p_q_phys, q_phys )
  pointer( p_phis, phis )
  pointer( p_ps, ps )
  pointer( p_omga, omga )
  pointer( p_pkz, pkz )
  pointer( p_pk, pk )
  pointer( p_pe, pe )
  pointer( p_peln, peln )
  pointer( p_pesouth, pesouth )
  pointer( p_ua, ua )
  pointer( p_va, va )
  pointer( p_ps_bp, ps_bp )
  pointer( p_u_srf, u_srf )
  pointer( p_v_srf, v_srf )
  pointer( p_u_dt, u_dt )
  pointer( p_v_dt, v_dt )
  pointer( p_t_dt, t_dt )
  pointer( p_q_dt, q_dt )
#ifdef MARS_GCM
  pointer( p_delp_dt, delp_dt )
  pointer( p_mars_sfc_budg, mars_sfc_budg )
#endif MARS_GCM
!point all the shared pointers
  p_u = ptr_u
  p_v = ptr_v
  p_delp = ptr_delp
  p_pt = ptr_pt
  p_q = ptr_q
  p_u_phys = ptr_u_phys
  p_v_phys = ptr_v_phys
  p_t_phys = ptr_t_phys
  p_q_phys = ptr_q_phys
  p_phis = ptr_phis
  p_ps = ptr_ps
  p_omga = ptr_omga
  p_pkz = ptr_pkz
  p_pk = ptr_pk
  p_pe = ptr_pe
  p_peln = ptr_peln
  p_pesouth = ptr_pesouth
  p_ua = ptr_ua
  p_va = ptr_va
  p_ps_bp = ptr_ps_bp
  p_u_srf = ptr_u_srf
  p_v_srf = ptr_v_srf
  p_u_dt = ptr_u_dt
  p_v_dt = ptr_v_dt
  p_t_dt = ptr_t_dt
  p_q_dt = ptr_q_dt
#ifdef MARS_GCM
  p_delp_dt = ptr_delp_dt
  p_mars_sfc_budg = ptr_mars_sfc_budg
#endif MARS_GCM

#endif  use_shared_pointers
