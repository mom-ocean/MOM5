
                 module mgrp_prscr_clds_mod

use fms_mod,                only: error_mesg,   &
                                  FATAL, &
                                  mpp_pe, mpp_root_pe, &
                                  write_version_number

use rad_utilities_mod,      only: cldrad_properties_type, &
                                  cld_specification_type


!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!       mgroup prescribed cloud properties module
!               (this module runnable in SKYHI and FMS; 
!                zonal_clouds_mod is FMS native equivalent)
!
!!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: mgrp_prscr_clds.F90,v 15.0 2007/08/14 03:55:11 fms Exp $'
  character(len=128)  :: tagname =  '$Name: siena_201207 $'



!---------------------------------------------------------------------
!-------  interfaces --------

public    mgrp_prscr_clds_init, &
          mgrp_prscr_clds_end,  &
          prscr_clds_amt,       &
          obtain_bulk_lw_prscr, &
          obtain_bulk_sw_prscr 


logical :: module_is_initialized = .false.

!----------------------------------------------------------------------
!----------------------------------------------------------------------

 contains 


subroutine mgrp_prscr_clds_init (    pref, latb      )

!------------------------------------------------------------------
real, dimension(:),   intent(in)             ::  latb      
real, dimension(:,:), intent(in)             :: pref          
!--------------------------------------------------------------------

      call write_version_number (version, tagname)

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('mgrp_prscr_clds_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine mgrp_prscr_clds_init

!######################################################################

subroutine mgrp_prscr_clds_end
        
!----------------------------------------------------------------------
!    mgrp_prscr_clds_end is the destructor for mgrp_prscr_clds_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
        
!--------------------------------------------------------------------
 
 
!---------------------------------------------------------------------

      call error_mesg('mgrp_prscr_clds_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine mgrp_prscr_clds_end




!#####################################################################

subroutine prscr_clds_amt (is, ie, js, je, Cld_spec)

!---------------------------------------------------------------------
!    prscr_clds_amt defines the location, amount (cloud fraction), 
!    number and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------

integer, intent(in)                          :: is, ie, js, je
type(cld_specification_type), intent(inout)  :: Cld_spec       

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
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

      call error_mesg('prscr_clds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine prscr_clds_amt


!######################################################################

subroutine obtain_bulk_lw_prscr (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_prscr defines bulk longwave cloud radiative 
!    properties for the mgrp_prscr_clds cloud scheme.
!---------------------------------------------------------------------

integer,                     intent(in)     :: is, ie, js, je
type(cld_specification_type), intent(inout) :: Cld_spec
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

!---------------------------------------------------------------------

      call error_mesg('obtain_bulk_lw_prscr', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_lw_prscr 



!#####################################################################

subroutine obtain_bulk_sw_prscr (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_zonal defines bulk shortwave cloud radiative 
!    properties for the zonal cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(cld_specification_type), intent(inout) :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props
!-------------------------------------------------------------------

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

!---------------------------------------------------------------------

      call error_mesg('obtain_bulk_sw_prscr', &
      'This module is not supported as part of the public release', FATAL)

end subroutine obtain_bulk_sw_prscr 

!####################################################################

       end module mgrp_prscr_clds_mod




