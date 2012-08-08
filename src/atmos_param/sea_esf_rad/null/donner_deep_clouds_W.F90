                 module donner_deep_clouds_W_mod

use time_manager_mod,       only: time_type
use       fms_mod,          only: error_mesg, FATAL, WARNING
use rad_utilities_mod,      only: microphysics_type

implicit none
private

   character(len=128)  :: version =  '$Id: donner_deep_clouds_W.F90,v 15.0 2007/08/14 03:55:08 fms Exp $'
   character(len=128)  :: tagname =  '$Name: siena_201207 $'

public :: donner_deep_clouds_W_init,  &
          donner_deep_clouds_W_end , donner_deep_clouds_amt

contains 

subroutine donner_deep_clouds_W_init  (pref, lonb, latb, axes, Time)

real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: lonb, latb
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time

      call error_mesg('donner_deep_clouds_W_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_clouds_W_init

!#################################################################

subroutine donner_deep_clouds_W_end
       
      call error_mesg('donner_deep_clouds_W_end', &
      'This module is not supported as part of the public release', WARNING)

end subroutine donner_deep_clouds_W_end

!#################################################################

subroutine donner_deep_clouds_amt (is, ie, js, je,   &
                   cell_cloud_frac, cell_liquid_amt, cell_liquid_size, &
                   cell_ice_amt, cell_ice_size, &
                   cell_droplet_number, &
                   meso_cloud_frac, meso_liquid_amt, meso_liquid_size, &
                   meso_ice_amt, meso_ice_size, &
                   meso_droplet_number, nsum_out, &
                   Cell_microphys,  Meso_microphys)

integer,                intent(in)    :: is,ie,js,je
real, dimension(:,:,:), intent(inout) ::   &
                   cell_cloud_frac, cell_liquid_amt, cell_liquid_size, &
                   cell_ice_amt, cell_ice_size, &
                   cell_droplet_number, &
                   meso_cloud_frac, meso_liquid_amt, meso_liquid_size, &
                   meso_ice_amt, meso_ice_size, &
                   meso_droplet_number
integer, dimension(:,:), intent(inout) ::  nsum_out
type(microphysics_type), intent(inout) :: Cell_microphys, Meso_microphys

      call error_mesg('donner_deep_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_clouds_amt  

!####################################################################

       end module donner_deep_clouds_W_mod
