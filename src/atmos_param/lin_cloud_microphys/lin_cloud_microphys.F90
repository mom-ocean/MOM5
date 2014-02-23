!
! Cloud micro-physics package for GFDL global cloud resolving model
! The algorithms are originally based on Lin et al 1983. Many key
! elements have been changed/improved based on several other publications
! Developer: Shian-Jiann Lin
!
module lin_cld_microphys_mod
 use time_manager_mod,  only: time_type
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file,  &
                              error_mesg, FATAL

 implicit none
 private

 public  lin_cld_microphys_driver, lin_cld_microphys_init, lin_cld_microphys_end
 public  qsmith_init, qsmith, g_sum, wqsat_moist, wqsat2_moist, sat_adj2

!---- version number -----
 character(len=128) :: version = '$Id: lin_cloud_microphys.F90,v 20.0.2.1 2013/12/17 19:46:01 Niki.Zadeh Exp $'
 character(len=128) :: tagname = '$Name: nullify_rab_nnz $'

 contains


  subroutine lin_cld_microphys_driver(qv,    ql,    qr,    qi,    qs,    qg,    qa,  &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
                               pt_dt, pt, p3, dz,  delp, area, dt_in,                &
                               land,  rain, snow, ice, graupel,                      &
                               hydrostatic, phys_hydrostatic,                        &
                               iis,iie, jjs,jje, kks,kke, ktop, kbot, time)

  type(time_type), intent(in):: time
  logical,         intent(in):: hydrostatic, phys_hydrostatic
  integer,         intent(in):: iis,iie, jjs,jje  ! physics window
  integer,         intent(in):: kks,kke           ! vertical dimension
  integer,         intent(in):: ktop, kbot        ! vertical compute domain
  real,            intent(in):: dt_in

  real, intent(in   ), dimension(:,:)  :: area
  real, intent(in   ), dimension(:,:)  :: land  !land fraction
  real, intent(out  ), dimension(:,:)  :: rain, snow, ice, graupel
  real, intent(in   ), dimension(:,:,:):: p3, delp, dz    ! p3 not used
  real, intent(in   ), dimension(:,:,:):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(:,:,:):: pt_dt,  qa_dt
  real, intent(inout), dimension(:,:,:):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                          qs_dt, qg_dt

  call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

  end subroutine lin_cld_microphys_driver


 subroutine sat_adj2(mdt, is, ie, js, je, ng, km, k, hydrostatic, consv_te, &
                     te0, qv, ql, qi, qr, qs, qa, area, peln, delz, pt, dp, last_step)
! This is designed for 6-class micro-physics schemes
! input pt is T_vir
 real, intent(in):: mdt
 integer, intent(in):: is, ie, js, je, km, ng, k
 logical, intent(in):: hydrostatic, last_step
 logical, intent(in):: consv_te
 real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: dp, area
 real, intent(in):: delz(is:ie,js:je)      ! Delta p at each model level
 real, intent(in):: peln(is:ie,km+1,js:je)           ! ln(pe)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng):: pt, qv, ql, qi, qr, qs, qa
 real, intent(inout):: te0(is:ie,js:je)

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end subroutine sat_adj2


 subroutine lin_cld_microphys_init(id, jd, kd, axes, time)

    integer,         intent(in) :: id, jd, kd
    integer,         intent(in) :: axes(4)
    type(time_type), intent(in) :: time

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end subroutine lin_cld_microphys_init


 subroutine lin_cld_microphys_end

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end subroutine lin_cld_microphys_end


 subroutine qsmith_init

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end subroutine qsmith_init


 real function wqsat2_moist(ta, qv, pa, dqdt)
  real, intent(in):: ta, pa, qv
  real, intent(out):: dqdt

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end function wqsat2_moist

 real function wqsat_moist(ta, qv, pa)
  real, intent(in):: ta, pa, qv

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end function wqsat_moist


 subroutine qsmith(im, km, ks, t, p, q, qs, dqdt)
  integer, intent(in):: im, km, ks
  real, intent(in),dimension(im,km):: t, p, q
  real, intent(out),dimension(im,km):: qs
  real, intent(out), optional:: dqdt(im,km)

  call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

  end subroutine qsmith

 real function g_sum(p, ifirst, ilast, jfirst, jlast, area, mode)
 use mpp_mod,           only: mpp_sum
 integer, intent(IN) :: ifirst, ilast
 integer, intent(IN) :: jfirst, jlast
 integer, intent(IN) :: mode  ! if ==1 divided by area
 real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
 real, intent(IN) :: area(ifirst:ilast,jfirst:jlast)

 call error_mesg ('lin_cloud_microphys_mod', 'lin_cloud_microphysics should not be active', FATAL)

 end function g_sum

end module lin_cld_microphys_mod
