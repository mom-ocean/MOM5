module clubb_driver_mod

! --- Imported from FMS modules ---

use       constants_mod, only: RAD_TO_DEG
use             mpp_mod, only: mpp_pe, mpp_root_pe, stdlog, mpp_chksum,        &
                               mpp_clock_id, mpp_clock_begin, mpp_clock_end,   &
                               CLOCK_MODULE_DRIVER
use    diag_manager_mod, only: register_diag_field, send_data
use    time_manager_mod, only: time_type, get_time, set_time, get_date,        &
                               operator(+), operator(-)
use             fms_mod, only: write_version_number, open_file,                &
                               open_namelist_file, check_nml_error,            &
                               file_exist, error_mesg, close_file,             &
                               read_data, write_data,                          &
                               mpp_error, FATAL, NOTE
use   field_manager_mod, only: MODEL_ATMOS
use  tracer_manager_mod, only: get_number_tracers, get_tracer_index,           &
                               get_tracer_names
use   rad_utilities_mod, only: aerosol_type
use     aer_ccn_act_mod, only: aer_ccn_act_init, aer_ccn_act_end
use   aer_ccn_act_k_mod, only: aer_ccn_act_k
use        ice_nucl_mod, only: ice_nucl_wpdf_init, ice_nucl_wpdf_end
implicit none 
public :: clubb_setup, & 
          clubb_init,  & 
          clubb,       & 
          clubb_end

!--------------------- version number ----------------------------------
character(len=128)   :: version = '$Id: CLUBB_driver_SCM.F90,v 1.1.6.2.2.2.2.1 2013/12/17 19:46:01 Niki.Zadeh Exp $'
character(len=128)   :: tagname = '$Name: nullify_rab_nnz $'

contains

! NULL routines return error if called but not compiled for clubb
  subroutine clubb(is, ie, js, je, lon, lat,                  &
                   Time_next,                                 &
                   dtmain,                                    &
                   phalf, pfull, zhalf, zfull, omega_avg,     &
                   t, q, r, u, v,                             &
                   u_star, b_star, q_star,                    &
                   tdt, qdt, rdt, udt, vdt,                   &
                   dcond_ls_liquid, dcond_ls_ice,             &
                   Ndrop_act_clubb, Icedrop_act_clubb,        &
                   ndust, rbar_dust,                          &
                   diff_t_clubb,                              &
                   qcvar_clubb,                               &
                   tdt_shf,  qdt_lhf ,                        &                   
                   Aerosol, mask,                             &
                   mc_full,                                   &
                   conv_frac_clubb,                           &
                   convective_humidity_ratio_clubb)

  integer, intent(in)                           ::  is, ie, js, je
  real, intent(in), dimension(:,:)              ::  lon, lat
  type(time_type), intent(in)                   ::  Time_next
  real, intent(in)                              ::  dtmain
  real, intent(in), dimension(:,:,:)            ::  phalf, pfull, zhalf, zfull, omega_avg
  real, intent(in), dimension(:,:,:)            ::  t, q, u, v
  real, intent(inout), dimension(:,:,:,:)       ::  r
  real, intent(in), dimension(:,:)              ::  u_star, b_star, q_star
  real, intent(inout), dimension(:,:,:)         ::  tdt, qdt, udt, vdt
  real, intent(inout), dimension(:,:,:,:)       ::  rdt
  real, intent(out), dimension(:,:,:)           ::  dcond_ls_liquid
  real, intent(out), dimension(:,:,:)           ::  dcond_ls_ice
  real, intent(out), dimension(:,:,:)           ::  Ndrop_act_clubb
  real, intent(out), dimension(:,:,:)           ::  Icedrop_act_clubb
  real, intent(out), dimension(:,:,:)           ::  ndust, rbar_dust
  real, intent(out), dimension(:,:,:)           ::  diff_t_clubb
  real, intent(out), optional, dimension(:,:,:) ::  qcvar_clubb   
  type(aerosol_type), intent(in), optional      ::  Aerosol
  real, intent(in), optional, dimension(:,:,:)  ::  mask
  real, intent(in), optional, dimension(:,:,:)  ::  mc_full
  real, intent(in), optional, dimension(:,:,:)  ::  conv_frac_clubb
  real, intent(in), optional, dimension(:,:,:)  ::  convective_humidity_ratio_clubb
  real, intent(in), optional, dimension(:,:)    ::  tdt_shf,  qdt_lhf

  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)

  end subroutine clubb

  !=====================================================================
  !=====================================================================

  subroutine clubb_init(id, jd, kd, lon, lat, axes, Time, phalf )

  integer, intent(in)                  :: id, jd, kd
  real, dimension(:,:), intent(in)     :: lon, lat
  integer, dimension(4), intent(in)    :: axes
  type(time_type), intent(in)          :: Time
  real, dimension(:,:,:), intent(in)   :: phalf
 
  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)

  end subroutine clubb_init
  !=====================================================================

  !=====================================================================
  subroutine clubb_setup(id, jd, phalf)

  !---------------------------------------------------------------------
  !  id, jd                             input
  !    subdomain dimensions
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !    [real, dimension(nlon,nlat,nlev+1)]
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                 :: id, jd
  real, dimension(:,:,:), intent(in)  ::  phalf

  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)

  end subroutine clubb_setup
  !=====================================================================

  !=====================================================================
  subroutine clubb_end

  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)

  end subroutine clubb_end
  !=====================================================================
end module  clubb_driver_mod
