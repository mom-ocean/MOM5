 
                 module standalone_clouds_mod

use fms_mod,                    only: mpp_pe, mpp_root_pe, &
                                      error_mesg, FATAL, &
                                      write_version_number
use rad_utilities_mod,          only: cld_specification_type, &
                                      cldrad_properties_type,  &
                                      microphysics_type,  &
                                      microrad_properties_type
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!   standalone cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: standalone_clouds.F90,v 15.0 2007/08/14 03:55:18 fms Exp $'
  character(len=128)  :: tagname =  '$Name: siena_201207 $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          standalone_clouds_init,                           &
          standalone_clouds_end,                           &
          define_column_properties, &
          standalone_clouds_amt, obtain_micro_lw_sa, obtain_micro_sw_sa,  &
          obtain_bulk_lw_sa, obtain_bulk_sw_sa

logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
                        contains 


!####################################################################

subroutine standalone_clouds_init (pref, lonb, latb)

!--------------------------------------------------------------------
!    subroutine standalone_clouds_init is the constructor for the
!    standalone_clouds_mod.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in)    ::  pref        
real, dimension(:,:), intent(in)    ::  lonb, latb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      array of model longitudes on cell boundaries 
!                 [ radians ]
!       latb      array of model latitudes at cell boundaries [radians]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('standalone_clouds_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine standalone_clouds_init



!####################################################################

subroutine define_column_properties (pref, lonb, latb)

!---------------------------------------------------------------------
!    subroutine define_column_properties defines values for lw emiss-
!    ivity, visible and nir reflectivity and nir absorption to be used
!    with standalone clouds.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
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

      call error_mesg('define_column_properties', &
      'This module is not supported as part of the public release', FATAL)

end subroutine define_column_properties

!######################################################################

subroutine standalone_clouds_end
        
!----------------------------------------------------------------------
!    standalone_clouds_end is the destructor for standalone_clouds_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
        
      call error_mesg('standalone_clouds_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine standalone_clouds_end

!#################################################################

subroutine standalone_clouds_amt (is, ie, js, je, lat, press_mks,  &
                                  Cld_spec)

!---------------------------------------------------------------------
!    standalone_clouds_amt defines the number, amount (cloud fraction), 
!    and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------

integer,                      intent(in)     ::  is, ie, js, je
real,    dimension(:,:),      intent(in)     ::  lat  
real,    dimension(:,:,:),    intent(in)     ::  press_mks
type(cld_specification_type), intent(inout)  ::  Cld_spec

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      press_mks    pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!
!   intent(inout) variables:
!
!      Cld_spec     cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!
!               the following elements of Cld_spec are defined here:
!
!                  %cmxolw  fraction of maximally overlapped clouds
!                           seen by the longwave radiation 
!                           [ dimensionless ]
!                  %crndlw  fraction of randomly overlapped clouds
!                           seen by the longwave radiation 
!                           [ dimensionless ]
!                  %camtsw  cloud fraction seen by the shortwave
!                           radiation; the sum of the maximally
!                           overlapped and randomly overlapped 
!                           longwave cloud fractions  [ dimensionless ]
!                  %nmxolw  number of maximally overlapped longwave 
!                           clouds in each grid column.
!                  %nrndlw  number of randomly overlapped longwave 
!                           clouds in each grid column.
!                  %ncldsw  number of clouds seen by he shortwave
!                           radiation in each grid column.
!                  %hi_cld  logical flag indicating the presence of 
!                           high clouds in a grid box
!                 %mid_cld  logical flag indicating the presence of 
!                           middle clouds in a grid box
!                 %low_cld  logical flag indicating the presence of 
!                           low clouds in a grid box
!                                                                  
!---------------------------------------------------------------------

      call error_mesg('standalone_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine standalone_clouds_amt   



!#####################################################################

subroutine obtain_micro_lw_sa (is, ie, js, je, Lsc_microphys, &
                               Meso_microphys, Cell_microphys, &
                               Lscrad_props,  Mesorad_props, &
                               Cellrad_props)

!---------------------------------------------------------------------
!    obtain_micro_lw_sa defines microphysically-based longwave cloud 
!    radiative properties when the code is executed in standalone 
!    columns mode.
!---------------------------------------------------------------------

integer,                        intent(in)    :: is, ie, js, je
type(microphysics_type),        intent(inout) :: Lsc_microphys, &
                                                 Meso_microphys, &
                                                 Cell_microphys
type(microrad_properties_type), intent(inout) :: Lscrad_props, &
                                                 Mesorad_props, &
                                                 Cellrad_props
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Lsc_microphys     microphysical specification for large-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Meso_microphys    microphysical specification for meso-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Cell_microphys    microphysical specification for cell-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Lscrad_props      cloud radiative properties on model grid,
!                        [ microrad_properties_type ]
!      Mesorad_props     meso-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!      Cellrad_props     cell-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!
!               the following component of the **_props variables is 
!               output from this routine:
!
!                    %abscoeff  absorption coefficient for  
!                               clouds in each of the longwave 
!                               frequency bands  [ km **(-1) ]
!
!---------------------------------------------------------------------

      call error_mesg('obtain_micro_lw_sa', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_micro_lw_sa     

!#####################################################################

subroutine obtain_micro_sw_sa (is, ie, js, je, Lsc_microphys,   &
                               Meso_microphys, Cell_microphys,   &
                               Lscrad_props, Mesorad_props,   &
                               Cellrad_props)

!--------------------------------------------------------------------
!    obtain_micro_sw_sa defines microphysically-based shortwave cloud 
!    radiative properties for the standalone cloud scheme when run in 
!    columns mode.
!---------------------------------------------------------------------

integer,                         intent(in)    ::  is, ie, js, je
type(microphysics_type),         intent(inout) ::  Lsc_microphys, &
                                                   Meso_microphys,   &
                                                   Cell_microphys
type(microrad_properties_type),  intent(inout) ::  Lscrad_props,   &
                                                   Mesorad_props,  &
                                                   Cellrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Lsc_microphys     microphysical specification for large-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Meso_microphys    microphysical specification for meso-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Cell_microphys    microphysical specification for cell-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Lscrad_props      large-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!      Mesorad_props     meso-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!      Cellrad_props     cell-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!
!               the following components of the microrad_properties
!               variables are output from this routine:
!
!                   %cldext    sw extinction coefficient for  
!                              clouds in each of the shortwave 
!                              frequency bands  [ km **(-1) ]
!                   %cldsct    sw scattering coefficient for
!                              clouds in each of the shortwave
!                              frequency bands  [ km **(-1) ]
!                   %cldasymm  sw asymmetry factor for
!                              clouds in each of the shortwave 
!                              frequency bands  [ dimensionless ]
!
!-----------------------------------------------------------------

      call error_mesg('obtain_micro_sw_sa', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_micro_sw_sa     

!#####################################################################

subroutine obtain_bulk_lw_sa (is, ie, js, je, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_sa defines bulk longwave cloud radiative properties 
!    when using specified clouds in the standalone columns mode.
!---------------------------------------------------------------------
 
integer,                      intent(in)    :: is, ie, js, je
type(cldrad_properties_type), intent(inout) :: Cldrad_props
 
!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
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

      call error_mesg('obtain_bulk_lw_sa', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_lw_sa     


!#####################################################################

subroutine obtain_bulk_sw_sa (is, ie, js, je, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_sa defines bulk shortwave cloud radiative 
!    properties for the specified cloud scheme when running in 
!    standalone columns mode.
!---------------------------------------------------------------------

integer,                      intent(in)    ::   is, ie, js, je
type(cldrad_properties_type), intent(inout) ::   Cldrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
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

      call error_mesg('obtain_bulk_sw_sa', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_sw_sa     

!#####################################################################

       end module standalone_clouds_mod




 
