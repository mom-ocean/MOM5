
module rh_clouds_mod

!=======================================================================
!
!                          RH_CLOUDS MODULE
!
!=======================================================================

use       fms_mod, only:  error_mesg, FATAL, file_exist,    &
                          check_nml_error, open_namelist_file,       &
                          close_file, &
                          read_data, write_data, mpp_pe, mpp_root_pe, &
                          write_version_number, stdlog

!=======================================================================

implicit none
private

public  rh_clouds, rh_clouds_init, rh_clouds_end,  &
        rh_clouds_sum, rh_clouds_avg, do_rh_clouds


interface rh_clouds
    module procedure  rh_clouds_3d, rh_clouds_2d, rh_clouds_1d
end interface

!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: rh_clouds.F90,v 10.0 2003/10/24 22:00:39 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'
logical            :: module_is_initialized = .false.

contains

!#######################################################################

subroutine rh_clouds_init (nlon, nlat, nlev)

integer, intent(in) :: nlon, nlat, nlev

!---------- output namelist to log-------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.
!---------------------------------------------------------------------

      call error_mesg('rh_clouds_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_clouds_init

!#######################################################################

subroutine rh_clouds_end

    module_is_initialized = .true.
!---------------------------------------------------------------------

      call error_mesg('rh_clouds_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_clouds_end

!#######################################################################

 function do_rh_clouds ( ) result (answer)
   logical :: answer

!  returns logical value for whether rh_clouds has been initialized
!  presumably if initialized then rh_cloud will be used

   answer = module_is_initialized

!---------------------------------------------------------------------

   call error_mesg('do_rh_clouds', &
      'This module is not supported as part of the public release', FATAL)

 end function do_rh_clouds

!#######################################################################

 subroutine rh_clouds_sum (is, js, rh)

!-----------------------------------------------------------------------
   integer, intent(in)                   :: is, js
      real, intent(in), dimension(:,:,:) :: rh
!-----------------------------------------------------------------------
!---------------------------------------------------------------------

      call error_mesg('rh_clouds_sum', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine rh_clouds_sum

!#######################################################################

 subroutine rh_clouds_avg (is, js, rh, ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(out), dimension(:,:,:) :: rh
   integer, intent(out)                   :: ierr
!---------------------------------------------------------------------

      call error_mesg('rh_clouds_avg', &
      'This module is not supported as part of the public release', FATAL)

 end subroutine rh_clouds_avg

!#######################################################################

subroutine rh_clouds_3d(rh, p_full, p_surf, zenith, deg_lat,&
            n_cloud,top,bot,cldamt,alb_uv,alb_nir,abs_uv,abs_nir,emiss)

real   , intent(in) , dimension(:,:,:)   :: rh, p_full
real   , intent(in) , dimension(:,:)     :: p_surf, zenith, deg_lat
integer, intent(out), dimension(:,:,:)   :: top, bot
integer, intent(out), dimension(:,:)     :: n_cloud
real   , intent(out), dimension(:,:,:)   :: cldamt,emiss
real   , intent(out), dimension(:,:,:)   :: alb_uv,alb_nir,abs_uv,abs_nir

!---------------------------------------------------------------------

      call error_mesg('rh_clouds_3d', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_clouds_3d

!#######################################################################
!  THE FOLLOWING CODE ALLOWS RH_CLOUDS TO BE USED IN 2D AND 1D MODELS
!#######################################################################

subroutine rh_clouds_2d(rh, p_full, p_surf, zenith, deg_lat,&
            n_cloud,top,bot,cldamt,alb_uv,alb_nir,abs_uv,abs_nir,emiss)

real   , intent(in) , dimension(:,:)   :: rh, p_full
real   , intent(in) , dimension(:)     :: p_surf,zenith,deg_lat
integer, intent(out), dimension(:,:)   :: top, bot
integer, intent(out), dimension(:)     :: n_cloud
real   , intent(out), dimension(:,:)   :: cldamt,alb_uv,alb_nir,abs_uv
real   , intent(out), dimension(:,:)   :: abs_nir,emiss

!---------------------------------------------------------------------

      call error_mesg('rh_clouds_2d', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_clouds_2d

!#######################################################################

subroutine rh_clouds_1d(rh, p_full, p_surf, zenith, deg_lat,&
            n_cloud,top,bot,cldamt,alb_uv,alb_nir,abs_uv,abs_nir,emiss)

real   , intent(in) , dimension(:)   :: rh, p_full
real   , intent(in)                  :: p_surf,zenith,deg_lat
integer, intent(out), dimension(:)   :: top, bot
integer, intent(out)                 :: n_cloud
real   , intent(out), dimension(:)   :: cldamt,alb_uv,alb_nir,abs_uv
real   , intent(out), dimension(:)   :: abs_nir,emiss

!---------------------------------------------------------------------

      call error_mesg('rh_clouds_1d', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_clouds_1d


!#######################################################################

end module rh_clouds_mod

