               module bulkphys_rad_mod
 
!    shared modules:

use fms_mod,                only: write_version_number,&
                                  error_mesg,   &
                                  FATAL
!    shared radiation package modules:

use rad_utilities_mod,      only: cldrad_properties_type, &
                                  cld_specification_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    bulkphys_rad_mod defines cloud radiative properties based on
!    bulk cloud physics values in contrast to microphysically-based
!    properties.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: bulkphys_rad.F90,v 15.0 2007/08/14 03:55:05 fms Exp $'
character(len=128)  :: tagname =  '$Name: siena_201207 $'



!---------------------------------------------------------------------
!-------  interfaces --------

public                                            &
          bulkphys_rad_init, bulkphys_lw_driver,  &
          bulkphys_sw_driver, bulkphys_rad_end

!-------------------------------------------------------------------
!    logical flag.
!-------------------------------------------------------------------
logical  :: module_is_initialized = .false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------



                           contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!#####################################################################

subroutine bulkphys_rad_init (min_cld_drop_rad_in, max_cld_drop_rad_in,&
                              min_cld_ice_size_in, max_cld_ice_size_in,&
                              pref, lonb, latb)

!---------------------------------------------------------------------
!    subroutine bulkphys_rad_init is the constructor for 
!    bulkphys_rad_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,                 intent(in) :: min_cld_drop_rad_in, &
                                    max_cld_drop_rad_in,&
                                    min_cld_ice_size_in,  &
                                    max_cld_ice_size_in
real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: lonb, latb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      array of model longitudes on cell boundaries 
!                 [ radians ]
!       latb      array of model latitudes at cell boundaries [radians]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('bulkphys_rad_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine bulkphys_rad_init



!#################################################################

subroutine bulkphys_sw_driver (is, ie, js, je, cosz, Cld_spec,   &
                              Cldrad_props)

!---------------------------------------------------------------------
!    bulkphys_sw_driver obtains bulk shortwave cloud radiative 
!    properties for the active cloud scheme.
!---------------------------------------------------------------------
 
integer,                      intent(in)    :: is, ie, js, je
real,    dimension(:,:),      intent(in)    :: cosz
type(cld_specification_type), intent(in)    :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle [ dimensionless ]
!      Cld_spec     cloud specification arrays defining the 
!                   location, amount and type (hi, middle, lo)
!                   of clouds that are present, provides input 
!                   to this subroutine
!                   [ cld_specification_type ]
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
!
!               the following components of this variable are output 
!               from this routine:
!
!                    %cirabsw   absorptivity of clouds in the 
!                               infrared frequency band
!                               [ dimensionless ]
!                    %cirrfsw   reflectivity of clouds in the 
!                               infrared frequency band
!                               [ dimensionless ]
!                    %cvisrfsw  reflectivity of clouds in the 
!                               visible frequency band
!                               [ dimensionless ]
!
!---------------------------------------------------------------------

      call error_mesg('bulkphys_sw_driver', &
      'This module is not supported as part of the public release', FATAL)

end subroutine bulkphys_sw_driver



!####################################################################

subroutine bulkphys_lw_driver (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    bulkphys_lw_driver defines bulk longwave cloud radiative 
!    properties for the active cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(cld_specification_type), intent(in)    :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Cld_spec          cloud specification arrays defining the 
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input 
!                        to this subroutine
!                        [ cld_specification_type ]
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
!
!               the following components of this variable are output 
!               from this routine:
!
!                    %emrndlw   longwave cloud emissivity for 
!                               randomly overlapped clouds
!                               in each of the longwave 
!                               frequency bands  [ dimensionless ]
!                    %emmxolw   longwave cloud emissivity for 
!                               maximally overlapped clouds
!                               in each of the longwave 
!                               frequency bands  [ dimensionless ]
!
!---------------------------------------------------------------------

      call error_mesg('bulkphys_lw_driver', &
      'This module is not supported as part of the public release', FATAL)

end subroutine bulkphys_lw_driver



!###################################################################
 
subroutine bulkphys_rad_end

!-------------------------------------------------------------------
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
     module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('bulkphys_rad_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine bulkphys_rad_end



                     end module bulkphys_rad_mod

 
