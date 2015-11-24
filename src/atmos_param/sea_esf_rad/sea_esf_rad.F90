                    module sea_esf_rad_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   smf
! </REVIEWER>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!   Code to initialize, commpute, and clean up radiation calculation. 
! </OVERVIEW>
! <DESCRIPTION>
!   The radiation component that initializes, deployes, and ends longwave,
!   shortwave, and diagnostics calculation in the FMS model.
! </DESCRIPTION>

!  shared modules:

use mpp_mod,              only: input_nml_file
use fms_mod,              only: open_namelist_file, fms_init, &
                                mpp_pe, mpp_root_pe, stdlog, &
                                file_exist, write_version_number, &
                                check_nml_error, error_mesg, &
                                FATAL, close_file, &
                                mpp_clock_id, mpp_clock_begin, &
                                mpp_clock_end, CLOCK_ROUTINE, &
                                CLOCK_MODULE
use time_manager_mod,     only: time_manager_init, time_type

!  shared radiation package modules:

use rad_utilities_mod,    only: rad_utilities_init, Rad_control, &
                                radiative_gases_type, & 
                                cldrad_properties_type, &
                                cld_specification_type,  &
                                astronomy_type, atmos_input_type, &
                                surface_type, lw_diagnostics_type, &
                                aerosol_diagnostics_type, &
                                cld_space_properties_type, &
                                lw_table_type, &
                                aerosol_type, aerosol_properties_type,&
                                sw_output_type, lw_output_type, &
                                Sw_control, Lw_parameters

!   radiation package modules:

use radiation_diag_mod,   only: radiation_diag_init,   &
                                radiation_diag_driver, &
                                radiation_diag_end
use longwave_driver_mod,  only: longwave_driver_init,   &
                                longwave_driver_time_vary, &
                                longwave_driver, &
                                longwave_driver_endts, &
                                longwave_driver_end
use shortwave_driver_mod, only: shortwave_driver_init,  &
                                shortwave_driver,  &
                                shortwave_driver_end

!----------------------------------------------------------------------

implicit none 
private 

!-----------------------------------------------------------------------
!    sea_esf_rad_mod is the driver for the sea_esf_rad radiation 
!    package.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!------------ version number for this module ---------------------------

character(len=128) :: version = '$Id: sea_esf_rad.F90,v 19.0 2012/01/06 20:23:33 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


!--------------------------------------------------------------------
!-- interfaces -----

public       &
            sea_esf_rad_init, sea_esf_rad, sea_esf_rad_time_vary,  &
            sea_esf_rad_endts,  sea_esf_rad_end


private      &

! called from sea_esf_rad:
             deallocate_arrays


!---------------------------------------------------------------------
!--- namelist ---

logical :: dummy


namelist /sea_esf_rad_nml/   &
                            dummy

!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----


logical :: module_is_initialized = .false.    ! module initialized ?
integer :: longwave_clock, shortwave_clock    ! timing clocks


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################
! <SUBROUTINE NAME="sea_esf_rad_init">
!   <OVERVIEW>
!     Routine to initialize the radiation calculation
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine initializes the utilities and radiation utilities
!     modules. Then it reads in the radiation namelist from the input
!     namelist file and log the namelist in an output log file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     CALL sea_esf_rad_init (lonb, latb, pref_r)
!   </TEMPLATE>
!
!   <IN NAME="lonb" TYPE="real">
!     2d array of model longitudes at cell corners in [radians]
!   </IN>
!   <IN NAME="latb" TYPE="real">
!     2d array of model latitudes at cell corners in [radians]
!   </IN>
!   <IN NAME="pref_r" TYPE="real">
!     Array containing two reference pressure profiles 
!     on the radiation grid for use in defining 
!     transmission functions in [pascals]
!   </IN>
! </SUBROUTINE>
subroutine sea_esf_rad_init (lonb, latb, pref_r)

!---------------------------------------------------------------------
!   sea_esf_rad_init is the constructor for sea_esf_rad_mod.
!---------------------------------------------------------------------

real, dimension(:,:),    intent(in)  :: lonb, latb
real, dimension(:,:),    intent(in)  :: pref_r

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb      2d array of model longitudes at cell corners 
!                 [radians]
!       latb      2d array of model latitudes at cell corners 
!                 [radians]
!       pref_r    array containing two reference pressure profiles 
!                 on the radiation grid for use in defining 
!                 transmission functions 
!                 [pascals]
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables

      integer                           :: unit, io, ierr, logunit
      type(lw_table_type)               :: Lw_tables

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!        end
!        Lw_tables
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call time_manager_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=sea_esf_rad_nml, iostat=io)
      ierr = check_nml_error(io,'sea_esf_rad_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=sea_esf_rad_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'sea_esf_rad_nml')
        end do
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                           write (logunit, nml=sea_esf_rad_nml)

!---------------------------------------------------------------------
!    initialize the modules called by this module.
!---------------------------------------------------------------------
      call longwave_driver_init  (latb, lonb, pref_r, Lw_tables)
      call shortwave_driver_init (latb, pref_r)
      call radiation_diag_init   (latb, lonb, Lw_tables)

!---------------------------------------------------------------------
!    initialize clocks to time various modules called by this module.
!---------------------------------------------------------------------
      longwave_clock =      &
                  mpp_clock_id ('   Physics_down: Radiation: lw', &
                        grain=CLOCK_ROUTINE)
      shortwave_clock =     &
                  mpp_clock_id ('   Physics_down: Radiation: sw', &
                        grain=CLOCK_ROUTINE)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------



end subroutine sea_esf_rad_init


!####################################################################
 
subroutine sea_esf_rad_time_vary (Time, Rad_gases_tv)

!----------------------------------------------------------------------
type(time_type), intent(in)  :: Time
type(radiative_gases_type), intent(inout) :: Rad_gases_tv

 
      call longwave_driver_time_vary (Time, Rad_gases_tv)
 

end subroutine sea_esf_rad_time_vary
 

!#######################################################################        ######

subroutine sea_esf_rad_endts (Rad_gases_tv)
 
type(radiative_gases_type), intent(in) :: Rad_gases_tv
 
    call longwave_driver_endts (Rad_gases_tv)

end subroutine sea_esf_rad_endts 



!#####################################################################
! <SUBROUTINE NAME="sea_esf_rad">
!   
!   <OVERVIEW>
!     The radiation component interface of the climate model
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calls longwave radiation computation subroutine, 
!     shortwave radiation computation subroutine, radiation diagnostics
!     computation routine, and finally it deallocates all previously
!     allocated memory spaces of temporary arrays.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call sea_esf_rad (is, ie, js, je, Atmos_input, Surface, Astro, Rad_gases, &
!    Aerosol, Cldrad_props, Cld_spec, Cld_diagnostics, Lw_output, Sw_output)
!   </TEMPLATE>
!
!   <IN NAME="is" TYPE="integer">
!     Starting subdomain i indice of data in the physics window being
!     modeled (longitudinal)
!   </IN>
!   <IN NAME="js" TYPE="integer">
!     Starting subdomain j indice of data in the physics window being
!     modeled (latitudinal)
!   </IN>
!   <IN NAME="ie" TYPE="integer">
!     Ending subdomain i indice of data in the physics window being
!     modeled  (longitudinal)
!   </IN>
!   <IN NAME="je" TYPE="integer">
!     Ending subdomain j indice of data in the physics window being
!     modeled (latitudinal)
!   </IN>
!   <IN NAME="Atmos_input" TYPE="atmos_input_type">
!     Atmos_input_type variable containing the atmospheric
!     input fields on the radiation grid 
!   </IN>
!   <IN NAME="Astro" TYPE="astronomy_type">
!     Astronomy_type variable containing the astronomical
!     input fields on the radiation grid  
!   </IN>
!   <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!     Radiative_gases_type variable containing the radiative 
!     gas input fields on the radiation grid 
!   </IN>
!   <IN NAME="Aerosol" TYPE="aerosol_type">
!     Aerosol input data to the shortwave radiation calculation
!   </IN>
!   <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!     The cloud radiative property input fields on the
!     radiation grid
!   </IN>
!   <IN NAME="Cld_diagnostics" TYPE="cld_diagnostics_type">
!     The cloud diagnostics input fields on the radiation grid
!   </IN>
!   <INOUT NAME="Lw_output" TYPE="lw_output_type">
!     The longwave radiation calculation result
!   </INOUT>
!   <INOUT NAME="Sw_output" TYPE="sw_output_type">
!     The shortwave radiation calculation result
!   </INOUT>
!   <IN NAME="Surface" TYPE="surface_type">
!    Surface data as boundary condition to radiation
!   </IN>
!   <IN NAME="Cld_spec" TYPE="cld_specification_type">
!    Cloud specification data as initial condition to radiation
!   </IN>
! </SUBROUTINE>

subroutine sea_esf_rad (is, ie, js, je, Rad_time, Atmos_input, Surface,&
                        Astro, Rad_gases, Aerosol, Aerosol_props,    &
                        Cldrad_props, Cld_spec, Lw_output, Sw_output, &
                        Aerosol_diags, r)

!-----------------------------------------------------------------------
!     sea_esf_rad calls the modules which calculate the long- and short-
!     wave radiational heating terms and fluxes and the radiation diag-
!     nostics module which provides radiation package diagnostics.
!-----------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
type(time_type),              intent(in)     :: Rad_time
type(atmos_input_type),       intent(in)     :: Atmos_input
type(surface_type),           intent(in)     :: Surface     
type(astronomy_type),         intent(in)     :: Astro
type(radiative_gases_type),   intent(inout)  :: Rad_gases
type(aerosol_type),           intent(in)     :: Aerosol      
type(aerosol_properties_type),intent(inout)  :: Aerosol_props
type(cldrad_properties_type), intent(in)     :: Cldrad_props
type(cld_specification_type), intent(in)     :: Cld_spec       
type(lw_output_type), dimension(:), intent(inout)  :: Lw_output
type(sw_output_type), dimension(:), intent(inout)  :: Sw_output 
type(aerosol_diagnostics_type), intent(inout)  :: Aerosol_diags
real, dimension(:,:,:,:),     intent(inout)  :: r
!---------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je   starting/ending subdomain i,j indices of data in 
!                    the physics_window being integrated
!      Rad_time      time at which the climatologically-determined, 
!                    time-varying input fields to radiation should 
!                    apply    
!                    [ time_type, days and seconds]
!      Atmos_input   atmospheric input fields          
!                    [ atmos_input_type ]
!      Surface       surface variables 
!                    [ surface_type ]
!      Astro         astronomical input fields            
!                    [ astronomy_type ]
!      Rad_gases     radiative gas input fields   
!                    [ radiative_gases_type ]
!      Aerosol       aerosol input fields 
!                    [ aerosol_type ]
!      Cldrad_props  cloud radiative property input fields
!                    [ cldrad_properties_type ]
!      Cld_spec      cloud specification input fields 
!                    [ cld_specification_type ]
!
!  intent(out) variables:
!
!      Aerosol_props aerosol radiative properties
!                    [ aerosol_properties_type ]
!      Lw_output     longwave radiation output data from the sea_esf_rad
!                    radiation package
!                    [ lw_output_type ]
!      Sw_output     shortwave radiation output data from the
!                    sea_esf_rad radiation package 
!                    [ sw_output_type ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      type(lw_diagnostics_type)         :: Lw_diagnostics
      type(cld_space_properties_type)   :: Cldspace_rad

!---------------------------------------------------------------------
!   local variables
!
!         Lw_diagnostics      used to hold desired diagnostics from 
!                             longwave_driver_mod so they may be passed
!                             to radiation_diag_mod
!                             [ lw_diagnostics_type ]
!         Cldspace_rad        used to hold output from bulk sw routine
!                             so that it may be passed to the 
!                             radiation_diag_mod    
!                             [ cld_space_properties_type ]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('sea_esf_rad_mod',   &
              'module has not been initialized', FATAL )
      endif

!----------------------------------------------------------------------
!    compute longwave radiation.
!----------------------------------------------------------------------
    if (Rad_control%do_lw_rad) then
      call mpp_clock_begin (longwave_clock)
      call longwave_driver (is, ie, js, je, Rad_time, Atmos_input,  &
                            Rad_gases, Aerosol, Aerosol_props,   &
                            Cldrad_props, Cld_spec, Aerosol_diags, &
                            Lw_output, Lw_diagnostics)
      call mpp_clock_end (longwave_clock)
    endif

!----------------------------------------------------------------------
!    compute shortwave radiation.
!----------------------------------------------------------------------
    if (Rad_control%do_sw_rad) then
      call mpp_clock_begin (shortwave_clock)
      call shortwave_driver (is, ie, js, je, Atmos_input, Surface,  &
                             Astro, Aerosol, Aerosol_props, Rad_gases, &
                             Cldrad_props, Cld_spec, Sw_output,   &
                             Cldspace_rad, Aerosol_diags, r)
      call mpp_clock_end (shortwave_clock)
    endif

!--------------------------------------------------------------------
!    call radiation_diag_driver to compute radiation diagnostics at 
!    desired points.
!--------------------------------------------------------------------
    if (Rad_control%do_sw_rad .and. Rad_control%do_lw_rad) then
      call radiation_diag_driver (is, ie, js, je, Atmos_input, Surface,&
                                  Astro, Rad_gases, Cldrad_props,   &
                                  Cld_spec, Sw_output, Lw_output, &
                                  Lw_diagnostics, Cldspace_rad)
    endif

!---------------------------------------------------------------------
!    call deallocate_arrays to deallocate the array components of the 
!    local derived-type variables.
!---------------------------------------------------------------------
      call deallocate_arrays (Lw_diagnostics, Cldspace_rad)

!--------------------------------------------------------------------

end subroutine sea_esf_rad




!###################################################################
! <SUBROUTINE NAME="sea_esf_rad_end">
! 
!   <OVERVIEW>
!     Ends radiation calculation.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine ends longwave, shortwave, and radiation
!     diagnostics calculation.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call sea_esf_rad_end
!   </TEMPLATE>
! </SUBROUTINE>

subroutine sea_esf_rad_end
 
!-------------------------------------------------------------------
!    sea_esf_rad_end is the destructor for the sea_esf_rad module.
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('sea_esf_rad_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    close out the modules initialized by this module.
!--------------------------------------------------------------------
      call longwave_driver_end
      call shortwave_driver_end
      call radiation_diag_end
 
!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine sea_esf_rad_end


!####################################################################
      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="deallocate_arrays">
!  
!   <OVERVIEW>
!     A routine to deallocate arrays allocated temporarily during
!     radiation calculation.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine deallocates arrays used in longwave 
!     diagnostics and cloud space parameters used in the
!     lacis-hansen formulation.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call deallocate_arrays (Lw_diagnostics, Cldspace_rad)
!   </TEMPLATE>
!
!   <IN NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!     Desired diagnostics from longwave_driver
!     so they may be passed to radiation_diag_mod
!   </IN>
!   <IN NAME="Cldspace_rad" TYPE="cld_space_properties_type">
!     Cld_space_properties_type variable which
!     holds lacis-hansen sw cloud-radiation
!     variables in cloud-space, rather than 
!     k-space, as the third dimension.
!   </IN>
! </SUBROUTINE>

subroutine deallocate_arrays (Lw_diagnostics, Cldspace_rad)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array cpomponents of local
!    derived-type variables.
!---------------------------------------------------------------------

type(lw_diagnostics_type),       intent(inout)   :: Lw_diagnostics
type(cld_space_properties_type), intent(inout)   :: Cldspace_rad

!---------------------------------------------------------------------
!  intent(inout) variables:
!
!         Lw_diagnostics      lw_diagnostics_type variable to hold
!                             desired diagnostics from longwave_driver
!                             so they may be passed to 
!                             radiation_diag_mod
!         Cldspace_rad        cld_space_properties_type variable which
!                             holds lacis-hansen sw cloud-radiation
!                             variables in cloud-space, rather than 
!                             k-space, as the third dimension.
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate the components of Lw_diagnostics.
!--------------------------------------------------------------------
    if (Rad_control%do_lw_rad) then
      deallocate (Lw_diagnostics%flx1e1)
      deallocate (Lw_diagnostics%fluxn )
      deallocate (Lw_diagnostics%cts_out)
      deallocate (Lw_diagnostics%cts_outcf)
      deallocate (Lw_diagnostics%gxcts )
      deallocate (Lw_diagnostics%excts )
      deallocate (Lw_diagnostics%exctsn)
      deallocate (Lw_diagnostics%fctsg )
      if (Lw_parameters%nbtrge > 0) then
        deallocate (Lw_diagnostics%flx1e1f)
      endif
      if (Rad_control%do_totcld_forcing) then
        deallocate (Lw_diagnostics%fluxncf)
      endif
   endif

!--------------------------------------------------------------------
!    deallocate the components of Cldspace_rad. these arrays are only
!    allocated when the lh sw code is called with clouds present; 
!    therefore one must test for pointer association before deallo-
!    cating the memory.
!--------------------------------------------------------------------
   if (Rad_control%do_sw_rad) then
      if (Sw_control%do_lhsw) then
        if (associated ( Cldspace_rad%camtswkc) ) then
          deallocate (Cldspace_rad%camtswkc )
          deallocate (Cldspace_rad%cirabswkc )
          deallocate (Cldspace_rad%cirrfswkc )
          deallocate (Cldspace_rad%cvisrfswkc )
          deallocate (Cldspace_rad%ktopswkc )
          deallocate (Cldspace_rad%kbtmswkc )
        endif
      endif
    endif

!--------------------------------------------------------------------


end subroutine deallocate_arrays 


!####################################################################



                 end module sea_esf_rad_mod

