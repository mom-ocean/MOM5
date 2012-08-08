
                 module diag_clouds_W_mod

use time_manager_mod,       only: time_type
use       fms_mod,          only: error_mesg, FATAL,&
                                  mpp_pe, mpp_root_pe, &
                                  write_version_number
use rad_utilities_mod,      only: microphysics_type, &
                                  cld_specification_type, &
                                  cldrad_properties_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!           diag cloud radiative properties module
!            currently a wrapper until SKYHI goes away and this
!            module can be consolidated with diag_cloud_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: diag_clouds_W.F90,v 12.0 2005/04/14 15:49:30 fms Exp $'
   character(len=128)  :: tagname =  '$Name: siena_201207 $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          diag_clouds_W_init,    &
          diag_clouds_W_end,    &
          diag_clouds_amt,  &
          obtain_bulk_lw_diag, &
          obtain_bulk_sw_diag

!----------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





subroutine diag_clouds_W_init  (num_slingo_bands_out)


integer, intent(out) :: num_slingo_bands_out

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('diag_clouds_W_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine diag_clouds_W_init


!#####################################################################

subroutine diag_clouds_W_end
 
!----------------------------------------------------------------------
!    diag_clouds_end is the destructor for diag_clouds_W_mod.
!----------------------------------------------------------------------
 
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------

      call error_mesg('diag_clouds_W_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine diag_clouds_W_end



!#################################################################

subroutine diag_clouds_amt (is, ie, js, je, lat, pflux, press,   &
                            Rad_time, Cld_spec, Lsc_microphys) 

!----------------------------------------------------------------------
!    diag_clouds_amt defines the location, amount (cloud fraction), 
!    number, optical depth, thickness and liquid percentage of clouds 
!    present on the model grid.
!----------------------------------------------------------------------

integer,                      intent(in)     ::  is, ie, js, je
real,    dimension(:,:),      intent(in)     ::  lat
real,    dimension(:,:,:),    intent(in)     ::  pflux, press
type(time_type),              intent(in)     ::  Rad_time     
type(cld_specification_type), intent(inout)  ::  Cld_spec
type(microphysics_type),      intent(inout)  ::  Lsc_microphys

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ] 
!      press        pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      Rad_time     time at which the climatologically-determined, 
!                   time-varying zonal cloud fields should apply
!                   [ time_type, days and seconds]
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
!                  %liq_frac 
!                           percentage of cloud condensate in a grid 
!                           box which is liquid  [ dimensionless ]
!                  %tau     cloud optical depth  [ dimensionless ]
!                  %cloud_thickness
!                           number of model layers over which the cloud
!                           in this grid box extends
!                  %ice_cloud  
!                           logical variable, which if true, indicates 
!                           that the grid box will contain ice cloud; 
!                           if false, the box will contain liquid cloud
!
!---------------------------------------------------------------------

      call error_mesg('diag_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine diag_clouds_amt 


!#####################################################################

subroutine obtain_bulk_lw_diag (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_diag defines bulk longwave cloud radiative 
!    properties for the gordon diag cloud scheme.
!---------------------------------------------------------------------
 
integer,                     intent(in)     :: is, ie, js, je
type(cld_specification_type), intent(in   ) :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Cld_spec          cloud specification arrays defining the 
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input 
!                        to this subroutine
!                        [ cld_specification_type ]
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

      call error_mesg('obtain_bulk_lw_diag', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_lw_diag




!#####################################################################

subroutine obtain_bulk_sw_diag (is, ie, js, je, cosz, Cld_spec,  &   
                                Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_diag defines bulk shortwave cloud radiative 
!    properties for the gordon diag cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    ::  is, ie, js, je
real, dimension(:,:),         intent(in)    ::  cosz
type(cld_specification_type), intent(in   ) ::  Cld_spec
type(cldrad_properties_type), intent(inout) ::  Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle  [ dimensionless ]
!
!   intent(inout) variables:
!
!      Cld_spec          cloud specification arrays defining the 
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input 
!                        to this subroutine
!                        [ cld_specification_type ]
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

      call error_mesg('obtain_bulk_sw_diag', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_sw_diag



!####################################################################


       end module diag_clouds_W_mod
       
