
                      module clouds_mod

!=======================================================================
!
!            determines cloud properties necessary for 
!                    fels-schwartzkopf radiation
!
!=======================================================================

use    cloud_rad_mod, only:  cloud_rad_init, cloud_summary
use  cloud_zonal_mod, only:  cloud_zonal
use    cloud_obs_mod, only:  cloud_obs, cloud_obs_init
use time_manager_mod, only:  time_type
use          fms_mod, only:  error_mesg, FATAL, file_exist,   &
                             check_nml_error, open_namelist_file,      &
                             mpp_pe, mpp_root_pe, close_file, &
                             write_version_number, stdlog
use    rh_clouds_mod, only:  do_rh_clouds, rh_clouds, rh_clouds_avg
use  strat_cloud_mod, only:  do_strat_cloud, strat_cloud_avg
use   diag_cloud_mod, only:  do_diag_cloud, diag_cloud_driver, &
                             diag_cloud_avg
use diag_manager_mod, only:  register_diag_field, send_data
use isccp_clouds_mod, only:  isccp_clouds_init

implicit none
private

!------------------- public interfaces ---------------------------------

public   clouds, clouds_init, clouds_end

!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------
 character(len=128) :: version = '$Id: clouds.F90,v 15.0 2007/08/14 03:52:54 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
!-----------------------------------------------------------------------

      logical :: module_is_initialized=.false.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine clouds  (is, js, clear_sky, Time, Time_diag, lat, &
                    land, tsfc, pfull, phalf, t, q, cosz,    &
                    nclds, ktopsw, kbtmsw, ktoplw, kbtmlw,   &
                    cldamt, cuvrf, cirrf, cirab, emcld, mask, kbot)

!-----------------------------------------------------------------------
        integer, intent(in)                    :: is, js
        logical, intent(in)                    :: clear_sky
type(time_type), intent(in)                    :: Time, Time_diag

   real, intent(in), dimension(:,:)    :: lat
   real, intent(in), dimension(:,:)    :: land,tsfc
   real, intent(in), dimension(:,:,:)  :: pfull,phalf,t,q
   real, intent(in), dimension(:,:)    :: cosz
integer, intent(out), dimension(:,:)   :: nclds
integer, intent(out), dimension(:,:,:) :: ktopsw,kbtmsw,ktoplw,kbtmlw
   real, intent(out), dimension(:,:,:) :: cldamt,cuvrf,cirrf,cirab,emcld
   real, intent(in),  dimension(:,:,:),optional :: mask
integer, intent(in),  dimension(:,:),  optional :: kbot
!-----------------------------------------------------------------------

      call error_mesg('clouds', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine clouds

!#######################################################################

      subroutine clouds_init ( lonb, latb, axes, Time )

!-----------------------------------------------------------------------
           real, intent(in), dimension(:,:) :: lonb, latb
        integer, intent(in), dimension(4)   :: axes
type(time_type), intent(in)                 :: Time

!-----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)


!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('clouds_init', &
      'This module is not supported as part of the public release', FATAL)


      end subroutine clouds_init

!#######################################################################

      subroutine clouds_end

!-----------------------------------------------------------------------
      module_is_initialized=.false.
!-----------------------------------------------------------------------

      call error_mesg('clouds_end', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine clouds_end

!#######################################################################

end module clouds_mod

