module MG_microp_3D_mod

  ! --- external modules ---

  use              fms_mod, only :  file_exist, open_namelist_file, close_file,    &
                                    error_mesg, FATAL, check_nml_error
  use     time_manager_mod, only :  time_type, get_time, set_date

  implicit none


  ! --- available public interfaces ---
  public  MG_microp_3D_init, MG_microp_3D, MG_microp_3D_end
        
contains

! NULL routines return error if called but not compiled for clubb
subroutine MG_microp_3D_init(axes,Time,idim,jdim,kdim)

  ! --- calling arguments ---
  integer, intent (in)        :: axes(4)          ! x,y,z,z_half axes types
  integer, intent (in)        :: idim,jdim,kdim   ! dimensions
  type(time_type), intent(in) :: Time             ! time
    call error_mesg('MG_microp_3D_init','Not compiled with -DCLUBB',FATAL)
end subroutine MG_microp_3D_init
!##############################################################################


!##############################################################################
subroutine MG_microp_3D( Time, is, ie, js, je, lon, lat, dtcloud,              &
                         pfull3d, phalf3d, zhalf3d, LAND,                      &
                         T3d, qv3d, ql3d, qi3d, qa3d, qn3d, qni3d, ahuco3d,    &
                         dcond_ls_liquid, dcond_ls_ice,                        &
                         Ndrop_act_CLUBB, Icedrop_act_CLUBB,                   &
                         ndust, rbar_dust,                                     &
                         ST3d, SQ3d, SL3d, SI3d, SA3d, SN3d, SNi3d,            &
                         rain3d, snow3d, surfrain, surfsnow,                   &
                         do_clubb,  qcvar_clubb, MASK3d,                       &
                         lsc_snow, lsc_rain, lsc_snow_size, lsc_rain_size )

  ! --- calling arguments ---
  type(time_type), intent (in)                         :: Time
  integer, intent (in)                                 :: is,ie,js,je
  real, intent (in),    dimension(:,:)                 :: lon,lat
  real, intent (in)                                    :: dtcloud
  real, intent (in),    dimension(:,:,:)               :: pfull3d,phalf3d
  real, intent (in),    dimension(:,:,:)               :: zhalf3d
  real, intent (in),    dimension(:,:)                 :: LAND
  real, intent (in),    dimension(:,:,:)               :: T3d,qv3d,ql3d,qi3d,qa3d,qn3d,qni3d
  real, intent (in),    dimension(:,:,:)               :: ahuco3d
  real, intent (in),    dimension(:,:,:)               :: dcond_ls_liquid,dcond_ls_ice
  real, intent (in),    dimension(:,:,:)               :: Ndrop_act_CLUBB,Icedrop_act_CLUBB
  real, intent (in),    dimension(:,:,:)               :: ndust, rbar_dust
  real, intent (out),   dimension(:,:,:)               :: ST3d,SQ3d,SL3d,SI3d,SA3d,SN3d,SNi3d
  real, intent (out),   dimension(:,:,:)               :: rain3d,snow3d
  real, intent (out),   dimension(:,:)                 :: surfrain,surfsnow

  integer, intent (in),  optional                      :: do_clubb
  real, intent (in),  optional, dimension(:,:,:)       :: qcvar_clubb
  real, intent (in),  optional, dimension(:,:,:)       :: MASK3d
  real, intent (out), optional, dimension(:,:,:)       :: lsc_snow,      &
                                                          lsc_rain,      &
                                                          lsc_snow_size, &
                                                          lsc_rain_size
    call error_mesg('MG_microp_3D','Not compiled with -DCLUBB',FATAL)
end subroutine MG_microp_3D
!##############################################################################


!##############################################################################
subroutine MG_microp_3D_end()
    call error_mesg('MG_microp_3D_end','Not compiled with -DCLUBB',FATAL)
end subroutine MG_microp_3D_end
!##############################################################################

end module MG_microp_3D_mod
