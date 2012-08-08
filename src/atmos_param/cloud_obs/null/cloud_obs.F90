
                    module cloud_obs_mod

!-----------------------------------------------------------------------
!
!           sets up observed (climatological) clouds
!
!-----------------------------------------------------------------------

use horiz_interp_mod, only: horiz_interp
use          fms_mod, only: file_exist, error_mesg, FATAL,  &
                            open_namelist_file, close_file,          &
                            check_nml_error, mpp_pe, mpp_root_pe, &
                            write_version_number, stdlog
use time_manager_mod, only: time_type, get_date
use  time_interp_mod, only: time_interp

implicit none
private

!---------- public interfaces ----------

public  cloud_obs, cloud_obs_init, cloud_obs_end

!-----------------------------------------------------------------------
!   ---------- private data ------------

   character(len=128) :: version = '$Id: cloud_obs.F90,v 15.0 2007/08/14 03:52:46 fms Exp $'
   character(len=128) :: tagname = '$Name: siena_201207 $'

   logical :: module_is_initialized = .false.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine cloud_obs ( is, js, Time, cldamt )

!-----------------------------------------------------------------------
!    routine that reads monthly records of climatological
!    isccp cloud amount and then linearly interpolates between them
!-----------------------------------------------------------------------
!     input
!     -----
!     is, js   starting i,j indices (dimension(2))
!     Time     current time (time_type)
!
!     output
!     ------
!     cldamt    cloud amount data on horizontal grid,
!               dimensioned ix x jx x 3, for high,med, & low clouds.
!-----------------------------------------------------------------------
        integer, intent(in)                    :: is, js
type(time_type), intent(in)                    :: Time
           real, intent(out), dimension(:,:,:) :: cldamt
!-----------------------------------------------------------------------

      call error_mesg('cloud_obs', &
      'This module is not supported as part of the public release', FATAL)


!-----------------------------------------------------------------------

 end subroutine cloud_obs

!#######################################################################

 subroutine cloud_obs_init (lonb,latb)

!-----------------------------------------------------------------------
!  lonb  =   longitude in radians of the grid box edges
!  latb  =   longitude in radians of the grid box edges
!-----------------------------------------------------------------------
   real, intent(in), dimension(:,:) :: lonb,latb
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)


!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('cloud_obs_init', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine cloud_obs_init

!#######################################################################

 subroutine cloud_obs_end
 
      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('cloud_obs_end', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine cloud_obs_end

!#######################################################################

end module cloud_obs_mod

