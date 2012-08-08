
module lscale_cond_mod

!-----------------------------------------------------------------------
use            fms_mod, only:  file_exist, error_mesg, open_namelist_file,  &
                               check_nml_error, mpp_pe, mpp_root_pe, FATAL,  &
                               close_file, write_version_number, stdlog
use sat_vapor_pres_mod, only:  escomp, descomp
use      constants_mod, only:  HLv,HLs,Cp_Air,Grav,rdgas,rvgas

implicit none
private
!-----------------------------------------------------------------------
!  ---- public interfaces ----

   public  lscale_cond, lscale_cond_init, lscale_cond_end

!-----------------------------------------------------------------------
!   ---- version number ----

 character(len=128) :: version = '$Id: lscale_cond.F90,v 10.0 2003/10/24 22:00:34 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
 logical            :: module_is_initialized=.false.


contains

!#######################################################################

   subroutine lscale_cond (tin, qin, pfull, phalf, coldT, &
                           rain, snow, tdel, qdel, mask, conv)

!-----------------------------------------------------------------------
!
!                      large scale condensation
!
!-----------------------------------------------------------------------
!
!   input:  tin      temperature at full model levels
!           qin      specific humidity of water vapor at full
!                      model levels
!           pfull    pressure at full model levels
!           phalf    pressure at half (interface) model levels
!           coldT    should precipitation be snow at this point?
!   optional:
!           mask     optional mask (0 or 1.) 
!           conv     logical flag; if true then no large-scale
!                       adjustment is performed at that grid-point or
!                       model level
!
!  output:  rain     liquid precipitation (kg/m2)
!           snow     frozen precipitation (kg/m2)
!           tdel     temperature tendency at full model levels
!           qdel     specific humidity tendency (of water vapor) at
!                      full model levels
!
!-----------------------------------------------------------------------
!--------------------- interface arguments -----------------------------

   real   , intent(in) , dimension(:,:,:) :: tin, qin, pfull, phalf
   logical   , intent(in) , dimension(:,:):: coldT
   real   , intent(out), dimension(:,:)   :: rain,snow
   real   , intent(out), dimension(:,:,:) :: tdel, qdel
   real   , intent(in) , dimension(:,:,:), optional :: mask
   logical, intent(in) , dimension(:,:,:), optional :: conv
!-----------------------------------------------------------------------

      call error_mesg('lscale_cond', &
      'This module is not supported as part of the public release', FATAL)

   end subroutine lscale_cond

!#######################################################################

   subroutine lscale_cond_init ()

!-----------------------------------------------------------------------
!
!        initialization for large scale condensation
!
!-----------------------------------------------------------------------

!---------- output namelist --------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized=.true.

!---------------------------------------------------------------------

      call error_mesg('lscale_cond_init', &
      'This module is not supported as part of the public release', FATAL)

   end subroutine lscale_cond_init

!#######################################################################

   subroutine lscale_cond_end

      module_is_initialized=.false.

!---------------------------------------------------------------------

      call error_mesg('lscale_cond_end', &
      'This module is not supported as part of the public release', FATAL)

   end subroutine lscale_cond_end

!#######################################################################

end module lscale_cond_mod

