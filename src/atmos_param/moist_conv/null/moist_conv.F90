module moist_conv_mod

use   time_manager_mod, only: time_type
use            fms_mod, only: error_mesg, FATAL, WARNING

implicit none
private

public :: moist_conv, moist_conv_init, moist_conv_end

CONTAINS

!#######################################################################

 subroutine moist_conv ( Tin, Qin, Pfull, Phalf, coldT,        &
                         Tdel, Qdel, Rain, Snow,               &
                         dtinv, Time, is, js, tracers, qtrmca, &
                         Lbot, mask, Conv,                     &
                         ql, qi, cf, qldel, qidel, cfdel)

    real, intent(INOUT), dimension(:,:,:)           :: Tin, Qin
    real, intent(IN) ,   dimension(:,:,:)           :: Pfull, Phalf
 logical, intent(IN) ,   dimension(:,:)             :: coldT
    real, intent(OUT),   dimension(:,:,:)           :: Tdel, Qdel
    real, intent(OUT),   dimension(:,:)             :: Rain, Snow
    real, intent(IN)                                :: dtinv
type(time_type), intent(in)                         :: Time
integer, intent(IN)                                :: is, js
    real, dimension(:,:,:,:), intent(in)            :: tracers
    real, dimension(:,:,:,:), intent(out)           :: qtrmca
 integer, intent(IN) ,   dimension(:,:),   optional :: Lbot
    real, intent(IN) ,   dimension(:,:,:), optional :: mask
 logical, intent(OUT),   dimension(:,:,:), optional :: Conv
    real, intent(INOUT), dimension(:,:,:), optional :: ql, qi, cf
    real, intent(OUT),   dimension(:,:,:), optional :: qldel, qidel, cfdel

         
      call error_mesg('moist_conv', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine moist_conv

!#######################################################################

 subroutine moist_conv_init (axes, Time, tracers_in_mca)

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 logical, dimension(:), intent(in), optional :: tracers_in_mca

      call error_mesg('moist_conv_init', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine moist_conv_init


!#######################################################################
 subroutine moist_conv_end

 call error_mesg('donner_deep', &
      'This module is not supported as part of the public release', WARNING)

 end subroutine moist_conv_end

!#######################################################################

end module moist_conv_mod
