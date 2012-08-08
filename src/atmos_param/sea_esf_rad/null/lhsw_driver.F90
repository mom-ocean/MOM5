      module lhsw_driver_mod

use rad_utilities_mod,     only: astronomy_type, &
                                 atmos_input_type, &
                                 surface_type, &
                                 sw_output_type, &
                                 cld_space_properties_type, &
                                 radiative_gases_type, &
                                 cld_specification_type, &
                                 cldrad_properties_type
use        fms_mod,        only: error_mesg, &  
                                 FATAL, &
                                 mpp_pe, mpp_root_pe, &
                                 write_version_number

!--------------------------------------------------------------------

implicit none
private


!--------------------------------------------------------------------
!            lacis-hansen shortwave parameterization
!
!------------------------------------------------------------------


!--------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

    character(len=128)  :: version =  '$Id: lhsw_driver.F90,v 12.0 2005/04/14 15:49:39 fms Exp $'
    character(len=128)  :: tagname =  '$Name: siena_201207 $'
    logical             :: module_is_initialized = .false.



!---------------------------------------------------------------------
!-------  interfaces --------
 
public  lhsw_driver_init, lhsw_driver_end, swrad

!------------------------------------------------------------------
!------------------------------------------------------------------

contains


subroutine lhsw_driver_init (          pref )

real, dimension(:,:), intent(in) :: pref

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

!------------------------------------------------------------------- 
      module_is_initialized = .true.
!---------------------------------------------------------------------

      call error_mesg('lhsw_driver_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine lhsw_driver_init





!######################################################################
subroutine lhsw_driver_end

!------------------------------------------------------------------- 
      module_is_initialized = .true.
!---------------------------------------------------------------------

      call error_mesg('lhsw_driver_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine lhsw_driver_end





!######################################################################
 
subroutine swrad ( is, ie, js, je,                        &
                   Astro,  with_clouds,    Atmos_input,   &
                   Surface,                               &
                   Rad_gases,                             &
                   Cldrad_props, Cld_spec, Sw_output, Cldspace_rad, gwt)
       
!-----------------------------------------------------------------------
!
!     Swrad solves for shortwave radiation.
!
!     references:
!
!     (1)  lacis, a. a. and j. e. hansen, "a parameterization for the
!          absorption of solar radiation in the earth's atmosphere," 
!          journal of the atmospheric sciences, 31 (1974), 118-133.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
!     intent in:
!
!     fracday =  fraction of day (or timestep) that sun is above 
!                horizon.
! 
!     press   =  pressure at data levels of model.
!
!     qo3     =  mass mixing ratio of o3 at model data levels.
!
!     rh2o    =  mass mixing ratio of h2o at model data levels.
!
!     ssolar  =  solar constant (may vary over one year). units: Wm-2.
!
! cosangsolar =  zenith angle at grid point.
!-----------------------------------------------------------------------

integer,                         intent(in)    :: is, ie, js, je
logical,                         intent(in)    :: with_clouds
type(cldrad_properties_type),    intent(in)    :: Cldrad_props
type(cld_specification_type),    intent(in)    :: Cld_spec       
real, dimension(:), optional,    intent(in)    :: gwt
type(astronomy_type),            intent(in)    :: Astro
type(atmos_input_type),          intent(in)    :: Atmos_input
type(surface_type),              intent(in)    :: Surface
type(radiative_gases_type),      intent(in)    :: Rad_gases  

type(sw_output_type),            intent(inout) :: Sw_output
type(cld_space_properties_type), intent(inout) :: Cldspace_rad

!-----------------------------------------------------------------------
!     intent out:
!
!     dfsw    =  downward radiation at all pressure levels.
!
!     fsw     =  net radiation (up-down) at all pressure levels.
!
!     hsw     =  radiation heating rates at all pressure layers.
!
!     ufsw    =  upward radiation at all pressure levels.
!-----------------------------------------------------------------------


!---------------------------------------------------------------------

      call error_mesg('swrad', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine swrad 


             end module lhsw_driver_mod

