                module original_fms_rad_mod


!-----------------------------------------------------------------------
!                 radiation interface module 
!-----------------------------------------------------------------------

use     time_manager_mod, only: time_type
use              fms_mod, only: FATAL,&
                                error_mesg, &
                                mpp_pe, mpp_root_pe,&
                                write_version_number

use    rad_utilities_mod, only:  radiative_gases_type, &
                                cldrad_properties_type, &
                                cld_specification_type, &
                                astronomy_type, &
                                atmos_input_type, &
                                surface_type, &
                                fsrad_output_type

implicit none 
private 

!----------- public interfaces in this module -----------------------

public    original_fms_rad_init, original_fms_rad_end, original_fms_rad

!-----------------------------------------------------------------------
!------------ version number for this module ---------------------------
character(len=128) :: version = '$Id: original_fms_rad.F90,v 15.0 2007/08/14 03:55:13 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'

     logical :: module_is_initialized = .false.

!---------------------------------------------------------------------
!---------------------------------------------------------------------

                         contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!####################################################################

   subroutine original_fms_rad_init ( lonb, latb, pref, axes, Time , &
                      kmax)

!-----------------------------------------------------------------------
           integer, intent(in)  :: kmax
           real, intent(in), dimension(:,:) :: lonb, latb
           real, intent(in), dimension(:,:) :: pref
        integer, intent(in), dimension(4)   :: axes
type(time_type), intent(in)                 :: Time

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('original_fms_rad_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine original_fms_rad_init

!#######################################################################

subroutine original_fms_rad_end

!-----------------------------------------------------------------------
module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('original_fms_rad_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine original_fms_rad_end
      
!###################################################################

subroutine original_fms_rad (is, ie, js, je, phalf, lat_in, lon_in, &
                             do_clear_sky_pass, &
                             Rad_time, Time_diag, Atmos_input, &
                             Surface, &
                             Astro, Rad_gases, Cldrad_props, Cld_spec, &
                             Fsrad_output, mask, kbot) 

integer,                      intent(in)           :: is, ie, js, je
real, dimension(:,:,:),       intent(in)           :: phalf
real, dimension(:,:),         intent(in)           :: lat_in, lon_in
type(time_type),              intent(in)           :: Rad_time, Time_diag
logical,                      intent(in)           :: do_clear_sky_pass
type(atmos_input_type),       intent(in)           :: Atmos_input
type(surface_type),           intent(in)           :: Surface
type(astronomy_type),         intent(in)           :: Astro        
type(radiative_gases_type),   intent(in)           :: Rad_gases
type(cldrad_properties_type), intent(in)           :: Cldrad_props
type(cld_specification_type), intent(in)           :: Cld_spec
type(fsrad_output_type),      intent(inout)        :: Fsrad_output
real, dimension(:,:,:),       intent(in), optional :: mask
integer, dimension(:,:),      intent(in), optional :: kbot
!--------------------------------------------------------------------

      call error_mesg('original_fms_rad', &
      'This module is not supported as part of the public release', FATAL)

end subroutine original_fms_rad

!#####################################################################

                 end module original_fms_rad_mod
