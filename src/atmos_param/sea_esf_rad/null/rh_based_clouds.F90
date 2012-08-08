
                 module rh_based_clouds_mod

use fms_mod,                only: mpp_pe, &
                                  mpp_root_pe, &
                                  write_version_number, &
                                  error_mesg,   &
                                  FATAL      
use rad_utilities_mod,      only: cldrad_properties_type, &
                                  cld_specification_type
                                 

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!           module which defines cloud locations
!                     based on model relative humidity
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: rh_based_clouds.F90,v 12.0 2005/04/14 15:49:51 fms Exp $'
  character(len=128)  :: tagname =  '$Name: siena_201207 $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          rh_based_clouds_init,  &
          rh_clouds_amt,  &
          obtain_bulk_lw_rh, obtain_bulk_sw_rh, &
          rh_based_clouds_end, &
          cldalb, albcld_lw, albcld_sw

!----------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.

!----------------------------------------------------------------------
!----------------------------------------------------------------------

                           contains 


subroutine rh_based_clouds_init 



!--------------------------------------------------------------------
     integer                           :: unit, ierr, io


!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('rh_based_clouds_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_based_clouds_init

!####################################################################

subroutine rh_based_clouds_end

!----------------------------------------------------------------------
!    rh_clouds_end is the destructor for rh_based_cloouds_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('rh_based_clouds_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_based_clouds_end



!######################################################################

subroutine rh_clouds_amt (is, ie, js, je, press, lat, Cld_spec)

!----------------------------------------------------------------------
!    rh_clouds_amt defines the location, amount (cloud fraction), number
!    and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------

integer,                      intent(in)    ::  is, ie, js, je
real,    dimension(:,:,:),    intent(in)    ::  press
real,    dimension(:,:),      intent(in)    ::  lat                    
type(cld_specification_type), intent(inout) ::  Cld_spec       

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      press        pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      lat          latitude of model points  [ radians ]
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
 
      call error_mesg('rh_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine rh_clouds_amt 

!####################################################################

subroutine obtain_bulk_lw_rh (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_rh defines bulk longwave cloud radiative 
!    properties for the rh cloud scheme.
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


      call error_mesg('obtain_bulk_lw_rh', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_lw_rh

!######################################################################

subroutine obtain_bulk_sw_rh (is, ie, js, je, cosz, Cld_spec,   &
                              Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_rh defines bulk shortwave cloud radiative 
!    properties for the rh cloud scheme.
!---------------------------------------------------------------------
 
integer,                      intent(in)    :: is, ie, js, je
real,    dimension(:,:),      intent(in)    :: cosz
type(cld_specification_type), intent(in   ) :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle [ dimensionless ]
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

      call error_mesg('obtain_bulk_sw_rh', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_sw_rh

!####################################################################

subroutine cldalb (zenith)

!---------------------------------------------------------------------
!     cldalb calculates a zenith angle dependency for the cloud albedos.
!     the cloud albedos are interpolated using data adapted from fritz 
!     (1954).  the solar zenith angle is the only input required.
!-----------------------------------------------------------------------

real, intent(in)           ::  zenith

!---------------------------------------------------------------------

      call error_mesg('cldalb', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cldalb

!##################################################################

subroutine albcld_lw(hi_cloud, mid_cloud, low_cloud,       &
             cmxolw, crndlw, emmxolw, emrndlw)

!-----------------------------------------------------------------------
!     albcld_lw computes the lw cloud emissivities. This calculation is 
!     based on sigma and cloud thickness in the old scheme (cldht60) 
!     and sigma, cloud thickness and latitude in the new scheme 
!     (cldht93).
!-----------------------------------------------------------------------

real, dimension(:,:,:),    intent(in)    :: cmxolw, crndlw
real, dimension(:,:,:,:),  intent(inout) :: emmxolw, emrndlw
logical, dimension(:,:,:), intent(in)    :: hi_cloud, mid_cloud,   &
                                           low_cloud
!---------------------------------------------------------------------

      call error_mesg('albcld_lw', &
      'This module is not supported as part of the public release', FATAL)

end subroutine albcld_lw

!####################################################################

subroutine albcld_sw(i,j, hi_cloud, mid_cloud, low_cloud,         &
     camtsw, cmxolw, crndlw, cvisrfsw, cirrfsw, cirabsw)

!-----------------------------------------------------------------------
!     albcld_sw computes the cloud albedos. This calculation is based on
!     sigma and cloud thickness in the old scheme (cldht60) and sigma, 
!     cloud thickness  and latitude in the new scheme (cldht93).
!-----------------------------------------------------------------------

real, dimension(:,:,:),    intent(in)    :: camtsw, cmxolw, crndlw
real, dimension(:,:,:),    intent(inout) :: cvisrfsw, cirrfsw, cirabsw
logical, dimension(:,:,:), intent(in)    :: hi_cloud, mid_cloud,   &
                                           low_cloud
integer,                   intent(in)    :: i, j
!---------------------------------------------------------------------

      call error_mesg('albcld_sw', &
      'This module is not supported as part of the public release', FATAL)

end subroutine albcld_sw



       end module rh_based_clouds_mod


