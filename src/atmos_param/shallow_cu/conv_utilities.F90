
MODULE CONV_UTILITIES_MOD
  
use      Constants_Mod, ONLY:   tfreeze, HLv, HLf, HLs, CP_AIR, GRAV, &
                                Kappa,rdgas,rvgas
use fms_mod,            only:   mpp_pe, mpp_root_pe
use conv_utilities_k_mod, only: uw_params_init_k, uw_params
! use conv_utilities_k_mod, only: sd_init_k, sd_copy_k, sd_end_k, &
!                                 ac_init_k, ac_clear_k, ac_end_k, &
!                                 uw_params_init_k, uw_params, &
!                                 pack_sd_k, extend_sd_k, adi_cloud_k,&
!                                 qsat_k, qses_k, exn_k, conden_k, &
!                                 findt_k,                      &
!                                 pack_sd_lsm_k, sounding, adicloud

!---------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: conv_utilities.F90,v 17.0 2009/07/21 02:58:03 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

! public  :: sd_init, sd_copy, sd_end, ac_init, ac_clear, ac_end, qsat, qses, exn, &
  public  ::   &
       uw_params_init
!      uw_params_init, &
!      conden, findt, pack_sd, pack_sd_lsm, extend_sd, adi_cloud

  real, parameter :: p00   = 1.E5
  real, parameter :: epsilo= rdgas/rvgas      !ratio of h2o to dry air molecular weights 
  real, parameter :: zvir  = rvgas/rdgas - 1. !rh2o/rair - 1
  real, parameter :: tkmin = -160 + tfreeze   ! tcmin from sat_vapor_pres.f90
  real, parameter :: tkmax =  100 + tfreeze   ! tcmax from sat_vapor_pres.f90

  character(len=7) :: mod_name = 'conv_utilities'

contains

!#####################################################################
!#####################################################################

  subroutine uw_params_init (Uw_p)

  type(uw_params), intent(inout) :: Uw_p
  
    integer :: me, root_pe
    
    me = mpp_pe()
    root_pe = mpp_root_pe()
    call uw_params_init_k (hlv, hls, hlf, cp_air, grav, kappa, rdgas, &
                           p00, epsilo, zvir, tkmin, tkmax, me,  &
                                                           root_pe,Uw_p)
    
  end subroutine uw_params_init

!#####################################################################
!#####################################################################



end MODULE CONV_UTILITIES_MOD
