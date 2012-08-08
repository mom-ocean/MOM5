module cosp_driver_mod

use time_manager_mod,       only: time_type

use fms_mod,                only: error_mesg, FATAL, WARNING

implicit none
private

character(len=128)  :: version =  '$Id: cosp_driver.F90,v 19.0 2012/01/06 20:04:34 fms Exp $'
character(len=128)  :: tagname =  '$Name: siena_201207 $'

public cosp_driver, cosp_driver_init, cosp_driver_end

contains

!######################################################################

subroutine cosp_driver_init (lonb, latb, Time_diag, axes,kd_in, ncol_in)

   real, dimension(:,:), intent(in) :: lonb, latb
   type(time_type), intent(in) :: Time_diag
   integer, dimension(4), intent(in) :: axes
   integer,               intent(in) :: kd_in, ncol_in

call error_mesg('cosp_driver_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cosp_driver_init



!#####################################################################


subroutine cosp_driver   &
        (lat_in, lon_in, daytime_in, phalf_plus, p_full_in, zhalf_plus,&
         z_full_in, u_wind_in, v_wind_in, mr_ozone_in, &
         T_in, sh_in, tca_in, cca_in, lsliq_in, lsice_in, ccliq_in, &
         ccice_in, fl_lsrain_in, fl_lssnow_in, fl_lsgrpl_in, &
         fl_ccrain_in,  &
         fl_ccsnow_in, reff_lsclliq_in, reff_lsclice_in,   &
         reff_lsprliq_in, reff_lsprice_in, reff_ccclliq_in,  &
         reff_ccclice_in, reff_ccprliq_in, reff_ccprice_in,  &
         skt_in, land_in, Time_diag, is, js, stoch_mr_liq_in, &
         stoch_mr_ice_in, stoch_size_liq_in, stoch_size_frz_in, &
         tau_stoch_in, lwem_stoch_in, stoch_cloud_type_in)
!--------------------------------------------------------------------
!    subroutine cosp_driver is the interface between the cosp simulator 
!    code and the AM model.
!--------------------------------------------------------------------
real, dimension(:,:),   intent(in) :: lat_in, lon_in, skt_in, land_in, &
                                      u_wind_in, v_wind_in
real, dimension(:,:), intent(in) :: daytime_in
real, dimension(:,:,:), intent(in) :: phalf_plus, p_full_in, &
        zhalf_plus, z_full_in, T_in, sh_in, &
        tca_in, cca_in, lsliq_in, lsice_in, ccliq_in, ccice_in, &
        fl_lsrain_in, fl_lssnow_in, fl_lsgrpl_in, fl_ccrain_in, &
        fl_ccsnow_in, mr_ozone_in, &
        reff_lsclliq_in, reff_lsclice_in, reff_lsprliq_in, &
        reff_lsprice_in, reff_ccclliq_in, reff_ccclice_in, &
        reff_ccprliq_in, reff_ccprice_in
real, dimension(:,:,:,:), intent(in), optional ::  &
               tau_stoch_in, lwem_stoch_in, stoch_cloud_type_in, &
               stoch_mr_liq_in, stoch_mr_ice_in, stoch_size_liq_in, &
               stoch_size_frz_in
type(time_type), intent(in) :: Time_diag
integer, intent(in) :: is, js

call error_mesg('cosp_driver', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cosp_driver



!#####################################################################

subroutine cosp_driver_end 

call error_mesg('cosp_driver_end', &
      'This module is not supported as part of the public release', FATAL)


end subroutine cosp_driver_end

!####################################################################

end module cosp_driver_mod
