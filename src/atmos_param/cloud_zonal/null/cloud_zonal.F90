
module cloud_zonal_mod

!=======================================================================
!
!       determines zonal cloud amounts and model levels.
!
!=======================================================================

use time_manager_mod, only:  time_type
use  time_interp_mod, only:  fraction_of_year
use          fms_mod, only:  error_mesg, FATAL, open_namelist_file, &
                             close_file, mpp_pe, mpp_root_pe, &
                             write_version_number

implicit none
private

public   cloud_zonal, cloud_zonal_init, cloud_zonal_end, getcld

!------------------- private data used by this module ------------------

   character(len=128) :: version = '$Id: cloud_zonal.F90,v 10.0 2003/10/24 22:00:25 fms Exp $'
   character(len=128) :: tagname = '$Name: siena_201207 $'
   logical            :: module_is_initialized=.false.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine cloud_zonal_init (season)

!-----------------------------------------------------------------------
!
!             initialization routine for retrieval of 
!             zonal cloud amounts and level indices.
!
!   input argument
!   --------------
!
!      season     scalar integer between 1-5
!                 where 1-4 uses fixed data (1=winter, 2=spring, etc.)
!                 season=5 is seasonal varying clouds
!

      integer, intent(in) :: season



!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)


!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('cloud_zonal_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cloud_zonal_init

!#######################################################################

subroutine getcld (time, lat, phalf, ktopsw, kbtmsw, cldamt)

!-----------------------------------------------------------------------
!
!  routine for retrieval of zonal cloud amounts and level indices.
!
!   input arguments
!   --------------
!
!      time       time of year (time_type)
!      lat        latitudes in radians, dimensioned by ncol   
!      phalf      pressure at model layer interfaces,
!                    dimensioned mxcol x nlev, although only
!                    the first ncol points of the first dimension
!                    are processed
!
!   output arguments
!   ----------------
!
!   (all output arguments are dimensioned ncol x 3; the second
!    dimension represents high, middle, and low clouds)
!
!      ktopsw     model layer interface indices for cloud tops
!      kbtmsw     model layer interface indices for cloud bottoms
!      cldamt     fractional cloud amounts
!

type(time_type), intent(in)  :: time
real,            intent(in)  :: lat(:,:), phalf(:,:,:)
integer,         intent(out) :: ktopsw(:,:,:),kbtmsw(:,:,:)
real , intent(out), optional :: cldamt(:,:,:)

!-----------------------------------------------------------------------

      call error_mesg('getcld', &
      'This module is not supported as part of the public release', FATAL)

!-----------------------------------------------------------------------

end subroutine getcld

!#######################################################################


subroutine cloud_zonal (time, lat, phalf,  &
                        nclds, ktopsw, kbtmsw, ktoplw, kbtmlw,  &
                        cldamt, cuvrf, cirrf, cirab, emcld)

!-----------------------------------------------------------------------
type(time_type), intent(in) :: time
           real, intent(in) :: lat(:,:), phalf(:,:,:)
integer, intent(out), dimension(:,:)   :: nclds
integer, intent(out), dimension(:,:,:) :: ktopsw,kbtmsw,ktoplw,kbtmlw
   real, intent(out), dimension(:,:,:) :: cldamt,cuvrf,cirrf,cirab,emcld
!-----------------------------------------------------------------------

      call error_mesg('cloud_zonal', &
      'This module is not supported as part of the public release', FATAL)

!-----------------------------------------------------------------------

end subroutine cloud_zonal

!#######################################################################


subroutine cloud_zonal_end

    module_is_initialized=.false.

!-----------------------------------------------------------------------

      call error_mesg('cloud_zonal_end', &
      'This module is not supported as part of the public release', FATAL)


end subroutine cloud_zonal_end

!#######################################################################

end module cloud_zonal_mod

