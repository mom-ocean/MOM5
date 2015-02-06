                 module cloudrad_package_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  Module that supplies cloud radiative properties
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

! shared modules:

use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: fms_init, open_namelist_file, &
                                    write_version_number, mpp_pe, &
                                    mpp_root_pe, stdlog, file_exist,  &
                                    check_nml_error, error_mesg,   &
                                    FATAL, close_file
use time_manager_mod,         only: time_type, time_manager_init

! shared radiation package modules:

use rad_utilities_mod,        only: rad_utilities_init, Lw_control, &
                                    Sw_control, cldrad_properties_type,&
                                    cld_specification_type, &
                                    microrad_properties_type, &
                                    microphysics_type, astronomy_type, &
                                    atmos_input_type, Cldrad_control
use esfsw_parameters_mod,     only: esfsw_parameters_init, Solar_spect

! radiation package modules:

use cloudrad_diagnostics_mod, only: cloudrad_diagnostics_init, &
                                    cloudrad_netcdf, &
                                    cloudrad_diagnostics_end
use bulkphys_rad_mod,         only: bulkphys_rad_init, &
                                    bulkphys_rad_end, &
                                    bulkphys_lw_driver, &
                                    bulkphys_sw_driver
use microphys_rad_mod,        only: lwemiss_calc, comb_cldprops_calc, &
                                    microphys_rad_init, &
                                    microphys_rad_end, &
                                    microphys_lw_driver, &
                                    microphys_sw_driver

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloudrad_package_mod computes cloud radiative properties consistent
!    with the activated radiation package options and returns them to
!    radiation_driver_mod.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloudrad_package.F90,v 20.0 2013/12/13 23:19:07 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloudrad_package_init, cloud_radiative_properties, &
         cldrad_props_dealloc, cloudrad_package_end


private          &
!  called from cloud_radiative_properties:
         initialize_cldrad_props, combine_cloud_properties,  &
         cloudrad_package_dealloc

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  :: microphys_form =  '     ' ! level of microphysics
                                                ! being used; either
                                                ! 'none', 'prescribed',
                                                ! or 'predicted'
real               :: min_cld_drop_rad = 4.2    ! smallest allowable 
                                                ! cloud drop radius 
                                                ! (microns) allowed in
                                                ! slingo scheme
real               :: max_cld_drop_rad = 16.6   ! largest allowable 
                                                ! cloud drop radius 
                                                ! (microns) allowed in
                                                ! slingo scheme
real               :: min_cld_ice_size = 18.6   ! smallest allowable 
                                                ! cloud ice size    
                                                ! (microns) allowed in
                                                ! fu-liou scheme
real               :: max_cld_ice_size = 130.2  ! largest allowable 
                                                ! cloud ice size    
                                                ! (microns) allowed in
                                                ! fu-liou scheme

namelist /cloudrad_package_nml /     &
                               min_cld_drop_rad, max_cld_drop_rad, &
                               min_cld_ice_size, max_cld_ice_size, &
                               microphys_form

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

logical      :: module_is_initialized = .false.  ! module initialized?



!----------------------------------------------------------------------
!----------------------------------------------------------------------



                          contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="cloudrad_package_init">
!  <OVERVIEW>
!   Contructor of cloudrad_package module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Contructor of cloudrad_package module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_package_init ( pref, lonb, latb, axes, Time)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference pressure levels
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   the longitude array of the model grid box corners
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   the latitude array of the model grid box corners
!  </IN>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_package_init (pref, lonb, latb, axes, Time, &
                                  donner_meso_is_largescale)

!---------------------------------------------------------------------
!    cloudrad_package_init is the constructor for cloudrad_package_mod.
!----------------------------------------------------------------------

real,    dimension(:,:), intent(in)    ::   pref
real,    dimension(:,:), intent(in)    ::   lonb, latb
integer, dimension(4),   intent(in)    ::   axes
type(time_type),         intent(in)    ::   Time
logical,                 intent(in)    ::   donner_meso_is_largescale

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      2d array of model longitudes on cell corners[ radians ]
!       latb      2d array of model latitudes at cell corners [radians]
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer         :: unit, io, ierr, logunit

!---------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file 
!      io       error status returned from io operation  
!      ierr     error code
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call rad_utilities_init

!---------------------------------------------------------------------
!    read namelist.
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cloudrad_package_nml, iostat=io)
      ierr = check_nml_error(io,"cloudrad_package_nml")
#else
!---------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=cloudrad_package_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'cloudrad_package_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                       write (logunit, nml=cloudrad_package_nml)
 
!-------------------------------------------------------------------
!    verify that Lw_control%do_lwcldemiss has been defined.
!-------------------------------------------------------------------
      if (Lw_control%do_lwcldemiss_iz) then
      else
        call error_mesg ('cloudrad_package_mod', &
         'Lw_control%do_lwcldemiss has not yet been defined', FATAL)
      endif

!-------------------------------------------------------------------
!    verify that Sw_control%do_esfsw has been defined.
!-------------------------------------------------------------------
      if (Sw_control%do_esfsw_iz) then
      else
        call error_mesg ('cloudrad_package_mod', &
         'Sw_control%do_esfsw has not yet been defined', FATAL)
      endif
!-------------------------------------------------------------------
!    verify that the component cloud scheme elements of Cldrad_control 
!    have been defined.
!-------------------------------------------------------------------
      if (Cldrad_control%do_rh_clouds_iz .and.  &
          Cldrad_control%do_strat_clouds_iz .and.  &
          Cldrad_control%do_zonal_clouds_iz .and.  &
          Cldrad_control%do_mgroup_prescribed_iz .and.  &
          Cldrad_control%do_obs_clouds_iz .and.  &
          Cldrad_control%do_no_clouds_iz .and.  &
          Cldrad_control%do_diag_clouds_iz .and.  &
          Cldrad_control%do_specified_clouds_iz .and.  &
          Cldrad_control%do_uw_clouds_iz .and.  &
          Cldrad_control%do_donner_deep_clouds_iz ) then  
      else
        call error_mesg ('cloudrad_package_mod', &
         'Cldrad_control%do_{some cloud type} not yet been defined', &
                                                                 FATAL)
      endif


!-------------------------------------------------------------------
!    define variables which denote whether microphysically-based cloud
!    radiative properties will be defined for the longwave and shortwave
!    radiation calculations. if esf sw is being used, then do_sw_micro
!    will be .true.; if do_lwcldemiss is .true. or if esf sw is being 
!    used along with strat clouds or diag clouds or specified clouds, 
!    then do_lw_micro will be .true..
!----------------------------------------------------------------------
      if (Lw_control%do_lwcldemiss) then
        Cldrad_control%do_lw_micro = .true.
      endif
      if (Sw_control%do_esfsw) then
        call esfsw_parameters_init
        Cldrad_control%do_sw_micro = .true.
        if (Cldrad_control%do_strat_clouds .or. &
            Cldrad_control%do_specified_clouds .or. &
            Cldrad_control%do_diag_clouds) then
          Cldrad_control%do_lw_micro = .true.
        endif
      endif

!--------------------------------------------------------------------
!    mark the logical controls as initialized.
!--------------------------------------------------------------------
      Cldrad_control%do_lw_micro_iz = .true.
      Cldrad_control%do_sw_micro_iz = .true.

!----------------------------------------------------------------------
!    define the microphysical use desired for this experiment. different
!    levels may be used for the lw and sw parameterizations. 
!----------------------------------------------------------------------
      if (trim(microphys_form) == 'predicted') then 

!---------------------------------------------------------------------
!    if microphys_form asks for predicted microphysics, then either
!    strat or donner deep clouds must be activated, and either one or 
!    both of do_lw_micro and do_sw_micro must be .true..  if these
!    conditions are met, set do_pred_cld_microphys to .true.. if only
!    one of do_sw_micro and do_lw_micro are true, then also set 
!    do_bulk_microphys to .true. so that the bulk scheme initialization
!    may be completed. if  neither do_lw_micro or do_sw_micro are .true.
!    or if a different cloud scheme has been activated, stop execution 
!    with an error message.
!---------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds .or. &
            Cldrad_control%do_donner_deep_clouds .or. &
            Cldrad_control%do_uw_clouds .or. &
            Cldrad_control%do_zetac_clouds) then
          if (Cldrad_control%do_sw_micro .and.   &
              Cldrad_control%do_lw_micro) then
            Cldrad_control%do_pred_cld_microphys = .true.
          else if (Cldrad_control%do_sw_micro .or.    &
                   Cldrad_control%do_lw_micro) then
            Cldrad_control%do_pred_cld_microphys = .true.
            Cldrad_control%do_bulk_microphys = .true.
          else
            call error_mesg( 'cloudrad_package_mod',  &
             ' not using microphysics -- set microphys_form '//&
                                                    'to none.', FATAL)
          endif
        else
          call error_mesg( 'cloudrad_package_mod',  &
                    ' predicted microphys not available with this '//&
                                                 'cloud scheme.', FATAL)
        endif

!---------------------------------------------------------------------
!    if prescribed microphysics are requested, make sure the cloud 
!    scheme requested has the capability of using the microphysical
!    properties, and that either the sw or lw scheme requested is 
!    microphysically based. if only one of do_sw_micro and do_lw_micro
!    is .true., then set do_bulk_microphys to .true., so that the bulk
!    scheme may be initialized. if neither is .true. or if a cloud
!    scheme has been requested that cannot use prescribed microphysics,
!    stop execution with an error message.
!---------------------------------------------------------------------
      else if (trim(microphys_form) == 'prescribed') then
        if (Cldrad_control%do_rh_clouds .or.  &
            Cldrad_control%do_mgroup_prescribed .or.  &
            Cldrad_control%do_specified_clouds  .or.  &
            Cldrad_control%do_diag_clouds       .or.  &
            Cldrad_control%do_no_clouds)    then
          if (Cldrad_control%do_sw_micro .and.    &
              Cldrad_control%do_lw_micro) then
            Cldrad_control%do_presc_cld_microphys = .true.
          else if (Cldrad_control%do_sw_micro .or.     &
                   Cldrad_control%do_lw_micro) then
            Cldrad_control%do_presc_cld_microphys = .true.
            Cldrad_control%do_bulk_microphys = .true.
          else
            call error_mesg( 'cloudrad_package_mod',  &
                ' not using microphysics -- set microphys_form '//&
                  'to none.', FATAL)
          endif
        else
          call error_mesg( 'cloudrad_package_mod',  &
             ' prescribed microphys not allowed with this cloud '//&
                                     'scheme.',  FATAL) 
        endif

!---------------------------------------------------------------------
!    if no microphysics is requested, make sure that donner_deep clouds
!    has not been requested (must use predicted cloud microphysics
!    for that scheme -- all others can be run without microphysics).  
!    also verify that the lw and sw schemes requested are not micro-
!    physically_based. if all is ok, set do_bulk_microphys to .true.; 
!    if not ok, write an error message and stop.
!---------------------------------------------------------------------
      else if (trim(microphys_form) == 'none') then
        if (Cldrad_control%do_donner_deep_clouds .or.  &
            Cldrad_control%do_uw_clouds) then        
          call error_mesg( 'cloudrad_package_mod',  &
            ' use predicted microphys with donner or uw clouds.', FATAL)
        else     
          if (Cldrad_control%do_sw_micro .or.    &
              Cldrad_control%do_lw_micro) then
            call error_mesg ('cloudrad_package_mod', &
               'must specify microphys_form when using microphysica'//&
                'lly-based cld rad scheme', FATAL)
          else
              Cldrad_control%do_bulk_microphys = .true.
          endif
        endif

!----------------------------------------------------------------------
!    error condition.
!----------------------------------------------------------------------
      else
        call error_mesg( 'cloudrad_package_mod',  &
           ' microphys_form is not an acceptable value.', FATAL)
      endif

!---------------------------------------------------------------------
!    define variables indicating that the desired cloud microphysics
!    type control variables have been defined.
!---------------------------------------------------------------------
      Cldrad_control%do_bulk_microphys_iz = .true.
      Cldrad_control%do_pred_cld_microphys_iz = .true.
      Cldrad_control%do_presc_cld_microphys_iz = .true.

!---------------------------------------------------------------------
!    if a microphysically-based scheme is being used, initialize the 
!    microphys_rad module.
!---------------------------------------------------------------------
      if (Cldrad_control%do_presc_cld_microphys  .or.  &
          Cldrad_control%do_pred_cld_microphys) then
        call microphys_rad_init ( min_cld_drop_rad, max_cld_drop_rad, &
                                  min_cld_ice_size, max_cld_ice_size, &
                                  axes, Time, lonb, latb)
      endif

!---------------------------------------------------------------------
!    if a bulk physics scheme is to be used, call bulkphys_rad_init. 
!---------------------------------------------------------------------
      if (Cldrad_control%do_bulk_microphys) then
        call bulkphys_rad_init (min_cld_drop_rad, max_cld_drop_rad, &
                                min_cld_ice_size, max_cld_ice_size, &
                                pref, lonb, latb)
      endif

!-------------------------------------------------------------------
!    when running the gcm or the standalone model with a cloud scheme,
!    call cloudrad_diagnostics_init to initialize the netcdf diagnostics
!    associated with the cloudrad package.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
        call cloudrad_diagnostics_init (min_cld_drop_rad,   &
                                        max_cld_drop_rad, &
                                  min_cld_ice_size, max_cld_ice_size, &
                                        axes, Time,    &
                                                donner_meso_is_largescale)
      endif

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------



end subroutine cloudrad_package_init



!####################################################################
! <SUBROUTINE NAME="cloud_radiative_properties">
!  <OVERVIEW>
!   Subroutine to calculate cloud radiative properties
!    appropriate for the radiation options that are active.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to calculate cloud radiative properties
!    appropriate for the radiation options that are active.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_radiative_properties (is, ie, js, je, Rad_time, Time, &
!                                       Time_next,  &
!                                       Astro, Atmos_input, Cld_spec,  &
!                                       Lsc_microphys, Meso_microphys, &
!                                       Cell_microphys, 
!                                     Shallow_microphys, Cldrad_props, &
!                                       kbot, mask)
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    Time      The current time.  [time_type]
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical properties needed by radiation
!                        package
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </IN>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties on model grid
!  </INOUT>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine cloud_radiative_properties (is, ie, js, je, Rad_time,   &
                                       Time, Time_next,  &
                                       Astro, Atmos_input, Cld_spec,  &
                                       Lsc_microphys, Meso_microphys, &
                                       Cell_microphys,   &
                                       Shallow_microphys,Cldrad_props, &
                                       Model_microphys, &
                                       kbot, mask)

!----------------------------------------------------------------------
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
type(time_type),              intent(in)             :: Rad_time, Time, &
                                                        Time_next
type(astronomy_type),         intent(in)             :: Astro
type(atmos_input_type),       intent(in)             :: Atmos_input
type(cld_specification_type), intent(inout)          :: Cld_spec    
type(microphysics_type),      intent(in)             :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys,&
                                                     Shallow_microphys
type(cldrad_properties_type), intent(inout)          :: Cldrad_props
type(microphysics_type),      intent(inout)          :: Model_microphys
integer, dimension(:,:),      intent(in),   optional :: kbot
real, dimension(:,:,:),       intent(in),   optional :: mask
!-------------------------------------------------------------------
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Time_next         time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!      Time              The current time.  [time_type]
!      Astro             astronomical properties needed by radiation
!                        package
!                        [ astronomy_type ]
!      Atmos_input       atmospheric input fields on model grid,
!                        [ atmos_input_type ] 
!      Cld_spec          cloud specification properties on model grid,
!                        [ cld_specification_type ]
!      Lsc_microphys     microphysical specification for large-scale 
!                        clouds
!                        [ microphysics_type ]
!      Meso_microphys    microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!                        [ microphysics_type ]
!      Cell_microphys    microphysical specification for convective cell
!                        clouds associated with donner convection
!                        [ microphysics_type ]
!      Shallow_microphys   
!                        microphysical specification for 
!                        clouds associated with uw shallow convection
!                        [ microphysics_type ]
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
!
!   intent(in), optional variables:
!
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!      mask              present when running eta vertical coordinate,
!                        mask to remove points below ground
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      type(microrad_properties_type) :: Lscrad_props, Cellrad_props, &
                                        Mesorad_props, Shallowrad_props
      integer  ::   ix, jx, kx  
      logical  ::   donner_flag = .true.
      logical  ::   donner_flag_uw = .false.

!---------------------------------------------------------------------
!   local variables:
!
!       Lscrad_props   cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!       Mesorad_props  cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!       Cellrad_props  cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!       Shallowrad_props   
!                      cloud radiative properties for 
!                      clouds associated with uw shallow convection  
!                      [ microrad_properties_type ]
!       ix             x dimension of current physics window
!       jx             y dimension of current physics window
!       kx             z dimension of current physics window
!       donner_flag    optional argument to microphys_rad module
!                      indicating that the meso or cell cloud compon-
!                      ent is currently being processed. needed because
!                      a different fu ice parameterization is used in
!                      these cases than is used for large-scale clouds.
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the dimensions of the current physics window.
!---------------------------------------------------------------------
      ix = size(Cld_spec%camtsw,1)
      jx = size(Cld_spec%camtsw,2)
      kx = size(Cld_spec%camtsw,3)

!--------------------------------------------------------------------
!    call initialize_cldrad_props to allocate and initialize the cloud 
!    radiative property arrays.
!---------------------------------------------------------------------
      call initialize_cldrad_props (ix, jx, kx, Lscrad_props,   &
                                    Mesorad_props, Cellrad_props, &
                                    Shallowrad_props, Cldrad_props)

!--------------------------------------------------------------------
!    if bulkphys_rad routines are needed, limit the condensate sizes
!    to that range acceptable for the radiative parameterizations.
!--------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro .and. &
        Cldrad_control%do_sw_micro ) then 
      else
        Cld_spec%reff_liq_lim = MAX(MIN(Cld_spec%reff_liq,  &
                                  max_cld_drop_rad), min_cld_drop_rad)
        Cld_spec%reff_ice_lim = MAX(MIN(Cld_spec%reff_ice,  &
                                  max_cld_ice_size), min_cld_ice_size)
      endif

!----------------------------------------------------------------------
!    if a cloud scheme is activated (in contrast to running without any
!    clouds), call either the microphysically-based or bulkphysics-based
!    modules to define the cloud lw and sw radiative properties. if the
!    model is being run with do_no_clouds = .true., exit from this 
!    routine, leaving the cloud radiative property variables as they 
!    were initialized (to values compatible to the non-existence of 
!    cloudiness).
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
        if (Cldrad_control%do_lw_micro) then
          if (Cldrad_control%do_ica_calcs) then
            call microphys_lw_driver (is, ie, js, je, Lsc_microphys,  &
                                      Cloud_rad_props=Cldrad_props)
          else
            call microphys_lw_driver (is, ie, js, je, Lsc_microphys,  &
                                      Micro_rad_props=Lscrad_props)
          endif
        else
          call bulkphys_lw_driver (is, ie, js, je, Cld_spec,    &
                                   Cldrad_props)
        endif
        if (Cldrad_control%do_sw_micro) then
          if (Cldrad_control%do_ica_calcs) then
            call microphys_sw_driver (is, ie, js, je, Lsc_microphys,  &
                                      Cloud_rad_props=Cldrad_props)
          else
            call microphys_sw_driver (is, ie, js, je, Lsc_microphys,  &
                                      Micro_rad_props=Lscrad_props)
          endif
        else
          call bulkphys_sw_driver (is, ie, js, je, Astro%cosz,   &
                                   Cld_spec, Cldrad_props)
        endif

!--------------------------------------------------------------------
!    if donner_deep_clouds is active, obtain the cloud radiative prop-
!    erties associated with the mesoscale and cell-scale convective
!    components. only microphysically-based properties are available.
!    the optional argument  donner_flag is used to indicate that prop-
!    erties associated with the clouds produced by the donner_deep_mod 
!    are being processed, since a different ice parameterization is
!    used for donner_deep relative to large-scale clouds.
!----------------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds) then
          donner_flag = .true.
          call microphys_lw_driver (is, ie, js, je, Meso_microphys,  &
                                    Micro_rad_props=Mesorad_props,   &
                                    donner_flag=donner_flag)
          call microphys_lw_driver (is, ie, js, je, Cell_microphys,  &
                                    Micro_rad_props=Cellrad_props, &
                                    donner_flag=donner_flag)
          call microphys_sw_driver (is, ie, js, je, Meso_microphys,  &
                                    Micro_rad_props=Mesorad_props, &
                                    donner_flag=donner_flag)
          call microphys_sw_driver (is, ie, js, je, Cell_microphys,  &
                                    Micro_rad_props=Cellrad_props, &
                                    donner_flag=donner_flag)
        endif

!--------------------------------------------------------------------
!    if the uw shallow convection scheme is active, obtain the cloud 
!    radiative properties associated with its clouds. only micro-
!    physically-based properties are available.
!    NOTE FOR NOW:
!   the optional argument  donner_flag is set to .false. when processing
!   shallow clouds. the ice cloud radiative properties are obtained from
!    the parameterization used by strat_cloud (effective size), rather 
!    than that used by donner_deep (generalized effective size).
!----------------------------------------------------------------------
       if (Cldrad_control%do_uw_clouds) then
         donner_flag_uw = .false.
         call microphys_lw_driver (is, ie, js, je, Shallow_microphys, &
                                   Micro_rad_props=Shallowrad_props, &
                                   donner_flag=donner_flag_uw)
         call microphys_sw_driver (is, ie, js, je, Shallow_microphys, &
                                   Micro_rad_props=Shallowrad_props,&
                                   donner_flag=donner_flag_uw)
        endif
      endif ! ( .not. do_no_clouds)

!---------------------------------------------------------------------
!    call combine_cloud_properties to define a set of total-cloud cloud
!    radiative properties. if donner deep and / or uw shallow clouds 
!    are active, this requires the combination of the cloud properties 
!    associated with the different types of cloud present (large-scale, 
!    donner meso and  cell, uw shallow). for other schemes the 
!    total-cloud values are simply the large-scale cloud values. 
!    this procedure is only needed when microphysically-based properties
!    are being used.
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_ica_calcs) then
        if (Cldrad_control%do_sw_micro  .or. &
            Cldrad_control%do_lw_micro) then  
          call combine_cloud_properties (is, js, Rad_time, Time_next, &
                                         Atmos_input%deltaz, &
                                         Cld_spec, &
                                         Lsc_microphys, Meso_microphys,&
                                         Cell_microphys,   &
                                         Shallow_microphys,  &
                                         Lscrad_props, &
                                         Mesorad_props, Cellrad_props, &
                                         Shallowrad_props, &
                                         Cldrad_props)
        endif  
      endif  

!----------------------------------------------------------------------
!    call lwemiss_calc to compute lw emissivity from the absorption 
!    coefficient when a microphysically-based lw emissivity scheme 
!    is being used.
!----------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        call lwemiss_calc (Atmos_input%clouddeltaz,   &
                           Cldrad_props%abscoeff, Cldrad_props%cldemiss)
        Cldrad_props%emmxolw = Cldrad_props%cldemiss
        Cldrad_props%emrndlw = Cldrad_props%cldemiss
      endif

!-------------------------------------------------------------------
!    if running the gcm or feedback program with a cloud scheme active,
!    call cloudrad_netcdf to generate netcdf output fields.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then 
          call cloudrad_netcdf (is, js, Time, Time_next, Atmos_input,&
                                Astro%cosz, Lsc_microphys, &
                                Meso_microphys, Cell_microphys,   &
                                Shallow_microphys, &
                                Lscrad_props, Mesorad_props,    &
                                Cellrad_props, Shallowrad_props, &
                                Cldrad_props, Cld_spec, Model_microphys)
      endif   ! (do_no_clouds)

!--------------------------------------------------------------------
!    call cloudrad_package_dealloc to deallocate the local derived type
!    variable arrays.
!--------------------------------------------------------------------
      call cloudrad_package_dealloc (Lscrad_props, Mesorad_props,   &
                                     Cellrad_props, Shallowrad_props)

!---------------------------------------------------------------------


end subroutine cloud_radiative_properties     



!#####################################################################
! <SUBROUTINE NAME="cldrad_properties_dealloc">
!  <OVERVIEW>
!    Subroutine to deallocate the array elements of the
!    cldrad_properties_type variable that is input.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Subroutine to deallocate the array elements of the
!    cldrad_properties_type variable that is input.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cldrad_props_dealloc (Cldrad_props)
!  </TEMPLATE>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cldrad_properties_type variable containing cloud 
!   radiative properties
!  </INOUT>
! </SUBROUTINE>
!
subroutine cldrad_props_dealloc (Cldrad_props)

!------------------------------------------------------------------
!    cldrad_props_dealloc deallocates the array elements of the
!    cldrad_properties_type variable that is input.
!------------------------------------------------------------------

type(cldrad_properties_type), intent(inout) :: Cldrad_props

!-------------------------------------------------------------------
!  intent(inout) variables:
!
!    Cldrad_props    cldrad_properties_type variable containing
!                    cloud radiative properties
!
!--------------------------------------------------------------------

!------------------------------------------------------------------
!    deallocate the array elements of Cldrad_props. different variables
!    exist dependent on the sw parameterization being used.
!-------------------------------------------------------------------
      deallocate (Cldrad_props%emmxolw   )
      deallocate (Cldrad_props%emrndlw   )
      deallocate (Cldrad_props%abscoeff  )
      deallocate (Cldrad_props%cldemiss  )

      if ( Cldrad_control%do_sw_micro ) then
        deallocate (Cldrad_props%cldext    )
        deallocate (Cldrad_props%cldasymm  )
        deallocate (Cldrad_props%cldsct    )
      else
        deallocate (Cldrad_props%cvisrfsw  )
        deallocate (Cldrad_props%cirabsw   )
        deallocate (Cldrad_props%cirrfsw   )
      endif

!-------------------------------------------------------------------


end subroutine cldrad_props_dealloc




!####################################################################
! <SUBROUTINE NAME="cloudrad_package_end">
!  <OVERVIEW>
!   Destructor of the cloudrad_package module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Destructor of the cloudrad_package module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_package_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cloudrad_package_end

!--------------------------------------------------------------------
!    cloudrad_package_end is the destructor for cloudrad_package_mod.
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    deactivate the modules which are component modules of
!    cloudrad_package_mod.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
        call cloudrad_diagnostics_end
      endif
      if (Cldrad_control%do_presc_cld_microphys  .or.  &
          Cldrad_control%do_pred_cld_microphys) then
        call microphys_rad_end
      endif
      if (Cldrad_control%do_bulk_microphys) then
        call bulkphys_rad_end
      endif

!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------



end subroutine cloudrad_package_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!####################################################################
! <SUBROUTINE NAME="initialize_cldrad_props">
!  <OVERVIEW>
!   initialize_cldrad_props allocates and initializes those fields
!    which define the cloud radiative properties needed by the
!    radiation package.
!  </OVERVIEW>
!  <DESCRIPTION>
!   initialize_cldrad_props allocates and initializes those fields
!    which define the cloud radiative properties needed by the
!    radiation package.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call initialize_cldrad_props (ix, jx, kx, Lscrad_props,    &
!                                    Mesorad_props, Cellrad_props, &
!                                    Shallowrad_props, Cldrad_props )
!  </TEMPLATE>
!  <IN NAME="ix, jx, kx" TYPE="integer">
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!  </IN>
!  <INOUT NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds   
!  </INOUT>
!  <INOUT NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale
!                      clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cellrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the convective cell
!                      clouds associated with donner convection 
!  </INOUT>
!  <INOUT NAME="Shallowrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the 
!                      clouds associated with uw shallow convection 
!  </INOUT>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties on model grid
!  </INOUT>
! </SUBROUTINE>
!
subroutine initialize_cldrad_props (ix, jx, kx, Lscrad_props,    &
                                    Mesorad_props, Cellrad_props, &
                                    Shallowrad_props, Cldrad_props )

!--------------------------------------------------------------------
!    initialize_cldrad_props allocates and initializes those fields
!    which define the cloud radiative properties needed by the
!    radiation package.
!---------------------------------------------------------------------

integer,                        intent(in)    :: ix, jx, kx
type(microrad_properties_type), intent(inout) :: Lscrad_props, &
                                                 Mesorad_props, &
                                                 Cellrad_props, &
                                                 Shallowrad_props
type(cldrad_properties_type),   intent(inout) :: Cldrad_props

!----------------------------------------------------------------------
!    intent(in) variables:
! 
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!
!   intent(inout) variables:
!
!       Lscrad_props   cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!       Mesorad_props  cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!       Cellrad_props  cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!       Shallowrad_props 
!                      cloud radiative properties for
!                      clouds associated with uw shallow convection  
!                      [ microrad_properties_type ]
!          the components of a microrad_structure are:
!            %cldext   parameterization band values of the cloud      
!                      extinction coefficient [ km**(-1) ]   
!            %cldsct   parameterization band values of the cloud      
!                      scattering coefficient [ km**(-1) ]
!            %cldasymm parameterization band values of the asymmetry  
!                      factor [ dimensionless ]
!            %abscoeff combined absorption coefficient for clouds in 
!                      each of the longwave frequency bands [ km**(-1) ]
!
!       Cldrad_props   cloud radiative properties on model grid,
!                      [ cldrad_properties_type ]
!          the components of a cldrad_properties_type strucure are:
!            %emmxolw  longwave cloud emissivity for maximally over-
!                      lapped clouds [ dimensionless ] 
!            %emrndlw  longwave cloud emissivity for randomly overlapped
!                      clouds  [ dimensionless ]
!            %cldext   parameterization band values of the cloud      
!                      extinction coefficient [ km**(-1) ]   
!            %cldsct   parameterization band values of the cloud      
!                      scattering coefficient [ km**(-1) ]
!            %cldasymm parameterization band values of the asymmetry  
!                      factor [ dimensionless ]
!            %abscoeff combined absorption coefficient for clouds in 
!                      each of the longwave frequency bands [ km**(-1) ]
!            %cldemiss longwave emissivity calculated using abscoeff
!                      [ dimensionless ]
!            %cirabsw  absorptivity of clouds in the infrared frequency 
!                      band. may be zenith angle dependent. 
!                      [ dimensionless ]
!            %cirrfsw  reflectivity of clouds in the infrared frequency
!                      band. may be zenith angle dependent.
!                      [ dimensionless ]
!            %cvisrfsw reflectivity of clouds in the visible frequency 
!                      band. may be zenith angle dependent.
!                      [ dimensionless ]
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variable:

      integer :: n_esfsw_bands
      integer :: nlwcldb

!-------------------------------------------------------------------
!   local variable:
!
!          n_esfsw_bands   number of spectral bands resolved by the 
!                          sw radiation package; set to zero when
!                          bulk-based sw (lhsw) is active
!          nlwcldb         number of frequency bands for which longwave
!                          emissivities are defined
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    allocate the arrays used to define the longwave cloud radiative
!    properties. initialize to appropriate non-cloudy values.
!--------------------------------------------------------------------
      nlwcldb = Cldrad_control%nlwcldb
      if (Cldrad_control%do_ica_calcs) then
        allocate (Cldrad_props%emmxolw  (ix, jx, kx, nlwcldb,nlwcldb) )
        allocate (Cldrad_props%emrndlw  (ix, jx, kx, nlwcldb,nlwcldb) )
        allocate (Cldrad_props%abscoeff (ix, jx, kx, nlwcldb,nlwcldb) )
        allocate (Cldrad_props%cldemiss (ix, jx, kx, nlwcldb,nlwcldb) )
      else
        allocate (Cldrad_props%emmxolw  (ix, jx, kx, nlwcldb,1) )
        allocate (Cldrad_props%emrndlw  (ix, jx, kx, nlwcldb,1) )
        allocate (Cldrad_props%abscoeff (ix, jx, kx, nlwcldb,1) )
        allocate (Cldrad_props%cldemiss (ix, jx, kx, nlwcldb,1) )
      endif
      Cldrad_props%emmxolw           = 1.0E+00
      Cldrad_props%emrndlw           = 1.0E+00
      Cldrad_props%abscoeff          = 0.0E+00
      Cldrad_props%cldemiss          = 0.0E+00

!---------------------------------------------------------------------
!    allocate and initialize the microphysically-based shortwave cloud 
!    radiative properties.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        n_esfsw_bands = Solar_spect%nbands 
        if (Cldrad_control%do_ica_calcs) then
          allocate (Cldrad_props%cldext  (ix, jx, kx, n_esfsw_bands, &
                    n_esfsw_bands) )
          allocate (Cldrad_props%cldsct  (ix, jx, kx, n_esfsw_bands, &
                    n_esfsw_bands) )
          allocate (Cldrad_props%cldasymm(ix, jx, kx, n_esfsw_bands, &
                    n_esfsw_bands) )
        else
          allocate (Cldrad_props%cldext  (ix, jx, kx, n_esfsw_bands, 1))
          allocate (Cldrad_props%cldsct  (ix, jx, kx, n_esfsw_bands, 1))
          allocate (Cldrad_props%cldasymm(ix, jx, kx, n_esfsw_bands, 1))
        endif
        Cldrad_props%cldsct            = 0.0E+00
        Cldrad_props%cldext            = 0.0E+00
        Cldrad_props%cldasymm          = 1.0E+00

!---------------------------------------------------------------------
!    allocate and initialize the bulk-based shortwave cloud 
!    radiative properties.
!---------------------------------------------------------------------
      else
        allocate (Cldrad_props%cirabsw (ix, jx, kx) )
        allocate (Cldrad_props%cirrfsw (ix, jx, kx) )
        allocate (Cldrad_props%cvisrfsw(ix, jx, kx) )
        Cldrad_props%cirrfsw (:,:,:) = 0.0E+00
        Cldrad_props%cvisrfsw(:,:,:) = 0.0E+00
        Cldrad_props%cirabsw (:,:,:) = 0.0E+00
      endif

!---------------------------------------------------------------------
!    allocate and initialize the cloud radiative properties associated
!    with large-scale clouds for those cases where mesoscale and cell-
!    scale clouds may also be present.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        allocate (Lscrad_props%cldext(ix, jx, kx, n_esfsw_bands) )
        allocate (Lscrad_props%cldsct(ix, jx, kx, n_esfsw_bands) )
        allocate (Lscrad_props%cldasymm(ix, jx, kx, n_esfsw_bands) )
        Lscrad_props%cldext   = 0.
        Lscrad_props%cldsct   = 0.
        Lscrad_props%cldasymm = 1.
      endif
      allocate (Lscrad_props%abscoeff (ix, jx, kx, nlwcldb) )
      Lscrad_props%abscoeff = 0.

!---------------------------------------------------------------------
!    allocate and initialize the cloud radiative properties associated
!    with mesoscale and cellscale clouds when they are present.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
        allocate (Cellrad_props%cldext(ix, jx, kx, n_esfsw_bands) )
        allocate (Cellrad_props%cldsct(ix, jx, kx, n_esfsw_bands) )
        allocate (Cellrad_props%cldasymm(ix, jx, kx, n_esfsw_bands) )
        allocate (Cellrad_props%abscoeff (ix, jx, kx, nlwcldb))
        allocate (Mesorad_props%cldext(ix, jx, kx, n_esfsw_bands) )
        allocate (Mesorad_props%cldsct(ix, jx, kx, n_esfsw_bands) )
        allocate (Mesorad_props%cldasymm(ix, jx, kx, n_esfsw_bands) )
        allocate (Mesorad_props%abscoeff (ix, jx, kx, nlwcldb) )
        Cellrad_props%cldext   = 0.
        Cellrad_props%cldsct   = 0.
        Cellrad_props%cldasymm = 1.
        Cellrad_props%abscoeff = 0.
        Mesorad_props%cldext   = 0.
        Mesorad_props%cldsct   = 0.
        Mesorad_props%cldasymm = 1.
        Mesorad_props%abscoeff = 0.
      endif

!---------------------------------------------------------------------
!    allocate and initialize the cloud radiative properties associated
!    with the clouds from the  uw shallow convection parameterization
!    when they are present.
!---------------------------------------------------------------------
       if (Cldrad_control%do_uw_clouds) then
         allocate (Shallowrad_props%cldext(ix, jx, kx, n_esfsw_bands) )
        allocate (Shallowrad_props%cldsct(ix, jx, kx, n_esfsw_bands) )
        allocate (Shallowrad_props%cldasymm(ix, jx, kx, n_esfsw_bands) )
        allocate (Shallowrad_props%abscoeff (ix, jx, kx, nlwcldb))
        Shallowrad_props%cldext   = 0.
        Shallowrad_props%cldsct   = 0.
        Shallowrad_props%cldasymm = 1.
        Shallowrad_props%abscoeff = 0.
      endif

!----------------------------------------------------------------------


end subroutine initialize_cldrad_props         



!#####################################################################
! <SUBROUTINE NAME="combine_cloud_properties">
!  <OVERVIEW>
!   combine_cloud_properties produces cloud-radiative properties fields
!    for the total-cloud field in each grid box.
!  </OVERVIEW>
!  <DESCRIPTION>
!   combine_cloud_properties produces cloud-radiative properties fields
!    for the total-cloud field in each grid box, using as input the 
!    properties and characteristics of the various cloud types that may 
!    be present (large-scale, donner mesoscale and cell-scale, uw 
!    shallow).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call combine_cloud_properties (is, js, Rad_time, deltaz,    &
!                                     Lsc_microphys, Meso_microphys,  &
!                                     Cell_microphys,   &
!                                     Shallow_microphys, &
!                                     Lscrad_props,   &
!                                     Mesorad_props,  Cellrad_props,  &
!                                     Shallowrad_props, &
!                                     Cldrad_props)
!  </TEMPLATE>
!  <IN NAME="ix, jx, kx" TYPE="integer">
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!    microphysical specification for large-scale 
!                      clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!    microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!    microphysical specification for  convective cell
!                      clouds associated with donner convection
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!    microphysical specification for 
!                      clouds associated with uw shallow convection
!  </IN>
!  <IN NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds   
!  </IN>
!  <IN NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale
!                      clouds   
!  </IN>
!  <IN NAME="Cellrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the convective cell
!                      clouds associated with donner convection 
!  </IN>
!  <IN NAME="Shallowrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the 
!                      clouds associated with uw shallow convection 
!  </IN>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties on model grid
!  </INOUT>
! </SUBROUTINE>
!
subroutine combine_cloud_properties (is, js, Rad_time, Time_next,  &
                                     deltaz,  Cld_spec,   &
                                     Lsc_microphys, Meso_microphys,  &
                                     Cell_microphys,   &
                                     Shallow_microphys, &
                                     Lscrad_props,   &
                                     Mesorad_props,  Cellrad_props,  &
                                     Shallowrad_props, &
                                     Cldrad_props)

!----------------------------------------------------------------------
!    combine_cloud_properties produces cloud-radiative properties fields
!    for the total-cloud field in each grid box, using as input the 
!    properties and characteristics of the various cloud types that may 
!    be present (large-scale, donner mesoscale and cell-scale, uw 
!    shallow).
!----------------------------------------------------------------------

integer,                        intent(in)    :: is, js
type(time_type),                intent(in)    :: Rad_time, Time_next
real, dimension(:,:,:),         intent(in)    :: deltaz
type(cld_specification_type),   intent(inout) :: Cld_spec
type(microphysics_type),        intent(in)    :: Lsc_microphys, &
                                                 Meso_microphys, &
                                                 Cell_microphys, &
                                                 Shallow_microphys
type(microrad_properties_type), intent(in)    :: Lscrad_props,  &
                                                 Mesorad_props,  &
                                                 Cellrad_props,&
                                                 Shallowrad_props
type(cldrad_properties_type), intent(inout)   :: Cldrad_props

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       Lsc_microphys  microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!       Meso_microphys microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!       Cell_microphys microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!       Shallow_microphys 
!                      microphysical specification for 
!                      clouds associated with uw shallow convection
!                      [ microphysics_type ]
!       Lscrad_props   cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!       Mesorad_props  cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!       Cellrad_props  cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!       Shallowrad_props  
!                      cloud radiative properties for
!                      clouds associated with uw shallow convection  
!                      [ microrad_properties_type ]
!
!    intent(inout) variables:
!
!      Cldrad_props    cloud radiative properties on model grid,
!                      [ cldrad_properties_type ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    call comb_cldprops_calc with the appropriate arguments dependent 
!    upon the cloud / convection scheme which is active.
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds .and.    &
          Cldrad_control%do_uw_clouds .and. &
          Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    if strat_cloud, donner_deep and uw shallow are all active, then 
!    lw and sw cloud radiative properties are microphysically based, 
!    and large-scale, donner mesoscale and cell-scale, and uw shallow
!    cloud properties are available. call comb_cldprops_calc to combine
!    these into a single set of cloud radiative properties to be used 
!    by the radiation package. 
!---------------------------------------------------------------------
        call comb_cldprops_calc (is, js, Rad_time, Time_next, deltaz,  &
                                 Cld_spec%stoch_cloud_type, &
                                 Cldrad_props%cldext,   &
                                 Cldrad_props%cldsct,   &
                                 Cldrad_props%cldasymm,  &
                                 Cldrad_props%abscoeff,   &
                                 Lsc_microphys = Lsc_microphys, &
                                 Meso_microphys = Meso_microphys, &
                                 Cell_microphys = Cell_microphys, &
                                 Shallow_microphys = Shallow_microphys,&
                                 Lscrad_props=Lscrad_props,  &
                                 Mesorad_props = Mesorad_props, &
                                 Cellrad_props = Cellrad_props, &
                                Shallowrad_props = Shallowrad_props)
     else if (Cldrad_control%do_strat_clouds .and.    &
             Cldrad_control%do_donner_deep_clouds) then
 

!----------------------------------------------------------------------
!    if strat_cloud and donner_deep are both active, then both lw and 
!    sw cloud radiative properties are microphysically based, and large-
!    scale, mesoscale and cell-scale properties are available. call
!    comb_cldprops_calc to combine these into a single set of cloud
!    radiative properties to be used by the radiation package. 
!---------------------------------------------------------------------
        call comb_cldprops_calc (is, js, Rad_time, Time_next, deltaz,   &
                                 Cld_spec%stoch_cloud_type, &
                                 Cldrad_props%cldext,   &
                                 Cldrad_props%cldsct,   &
                                 Cldrad_props%cldasymm,  &
                                 Cldrad_props%abscoeff,   &
                                 Lsc_microphys = Lsc_microphys, &
                                 Meso_microphys = Meso_microphys, &
                                 Cell_microphys = Cell_microphys, &
                                 Lscrad_props=Lscrad_props,  &
                                 Mesorad_props = Mesorad_props, &
                                 Cellrad_props = Cellrad_props)

      else if (Cldrad_control%do_strat_clouds .and.    &
               Cldrad_control%do_uw_clouds) then

!----------------------------------------------------------------------
!    if strat_cloud and uw shallow are both active, then both lw and 
!    sw cloud radiative properties are microphysically based, and large-
!    scale and uw shallow cloud properties are available. call
!    comb_cldprops_calc to combine these into a single set of cloud
!    radiative properties to be used by the radiation package. 
!---------------------------------------------------------------------
        call comb_cldprops_calc (is, js, Rad_time, Time_next, deltaz,  &
                                 Cld_spec%stoch_cloud_type, &
                                  Cldrad_props%cldext,   &
                               Cldrad_props%cldsct,   &
                                 Cldrad_props%cldasymm,  &
                                  Cldrad_props%abscoeff,   &
                                  Lsc_microphys = Lsc_microphys, &
                                 Shallow_microphys = Shallow_microphys,&
                                  Lscrad_props=Lscrad_props,  &
                                  Shallowrad_props = Shallowrad_props)

   else if (Cldrad_control%do_uw_clouds .and.    &
            Cldrad_control%do_donner_deep_clouds) then
 
!----------------------------------------------------------------------
!    if uw shallow and donner_deep clouds are both active, then both 
!    lw and sw cloud radiative properties are microphysically based, 
!    and donner mesoscale and cell-scale and uw shallow cloud properties
!    are available. call comb_cldprops_calc to combine these into a 
!    single set of cloud radiative properties to be used by the 
!    radiation package. 
!---------------------------------------------------------------------
      call comb_cldprops_calc (is, js, Rad_time, Time_next, deltaz,  &
                                 Cld_spec%stoch_cloud_type, &
                                Cldrad_props%cldext,   &
                                 Cldrad_props%cldsct,   &
                                Cldrad_props%cldasymm,  &
                                  Cldrad_props%abscoeff,   &
                                Meso_microphys = Meso_microphys, &
                                 Cell_microphys = Cell_microphys, &
                                 Shallow_microphys = Shallow_microphys,&
                                 Mesorad_props = Mesorad_props, &
                                 Cellrad_props = Cellrad_props, &
                                Shallowrad_props=Shallowrad_props)
      

!---------------------------------------------------------------------
!    if donner_deep alone is active, then the mesoscale
!    and cell-scale properties must be combined. 
!----------------------------------------------------------------------
      else if (Cldrad_control%do_donner_deep_clouds) then
        call comb_cldprops_calc (is, js, Rad_time, Time_next, deltaz,  &
                                 Cld_spec%stoch_cloud_type, &
                                 Cldrad_props%cldext,   &
                                 Cldrad_props%cldsct,   &
                                 Cldrad_props%cldasymm, &
                                 Cldrad_props%abscoeff,   &
                                 Meso_microphys = Meso_microphys, &
                                 Cell_microphys = Cell_microphys, &
                                 Mesorad_props = Mesorad_props, &
                                 Cellrad_props = Cellrad_props)

!---------------------------------------------------------------------
!    if uw shallow alone is active, then total cloud values are 
!    defined as the uw shallow values.
!----------------------------------------------------------------------
     else if (Cldrad_control%do_uw_clouds) then
       if (Cldrad_control%do_sw_micro) then
         Cldrad_props%cldsct(:,:,:,:,1) =    &
                                      Shallowrad_props%cldsct(:,:,:,:)
         Cldrad_props%cldext(:,:,:,:,1) =    &
                                      Shallowrad_props%cldext(:,:,:,:)
         Cldrad_props%cldasymm(:,:,:,:,1) =  &
                                     Shallowrad_props%cldasymm(:,:,:,:)
      endif
      if (Cldrad_control%do_lw_micro) then
        Cldrad_props%abscoeff(:,:,:,:,1) =   &
                                  Shallowrad_props%abscoeff(:,:,:,:)
     endif
      
!----------------------------------------------------------------------
!    if microphysically-based properties have been generated without
!    donner_deep being active, total-cloud values are defined as the
!    large-scale cloud values.
!----------------------------------------------------------------------
      else if (Cldrad_control%do_strat_clouds) then
        if (Cldrad_control%do_sw_micro) then
          Cldrad_props%cldsct(:,:,:,:,1) = Lscrad_props%cldsct(:,:,:,:)
          Cldrad_props%cldext(:,:,:,:,1) = Lscrad_props%cldext(:,:,:,:)
          Cldrad_props%cldasymm(:,:,:,:,1) =  &
                                         Lscrad_props%cldasymm(:,:,:,:)
        endif
        if (Cldrad_control%do_lw_micro) then
          Cldrad_props%abscoeff(:,:,:,:,1) =   &
                                      Lscrad_props%abscoeff(:,:,:,:)
        endif
      endif

!--------------------------------------------------------------------



end subroutine  combine_cloud_properties 



!####################################################################
! <SUBROUTINE NAME="cloudrad_package_dealloc">
!  <OVERVIEW>
!   Subroutine to deallocate the space cloud radiative properties use
!   in the model
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to deallocate the space cloud radiative properties use
!   in the model
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  cloudrad_package_dealloc (Lscrad_props, Mesorad_props,  &
!                                   Cellrad_props)
!  </TEMPLATE>
!  <IN NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds   
!  </IN>
!  <IN NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale
!                      clouds   
!  </IN>
!  <IN NAME="Cellrad_props" TYPE="cldrad_prperties_type">
!   cldrad_prperties_type variable containing the cloud radiative
!   properties for the donner cell clouds
!  </IN>
!  <IN NAME="Shallowrad_props" TYPE="cldrad_properties_type">
!   cldrad_prperties_type variable containing the cloud radiative
!   properties for the uw shallow clouds
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_package_dealloc (Lscrad_props, Mesorad_props,   &
                                     Cellrad_props, Shallowrad_props)

!---------------------------------------------------------------------
!    cloudrad_package_dealloc deallocates the components of the local
!    derived-type variables.
!---------------------------------------------------------------------
        
type(microrad_properties_type), intent(inout) :: Lscrad_props,  &
                                                 Mesorad_props, &
                                                 Cellrad_props, &
                                                 Shallowrad_props
!---------------------------------------------------------------------
!   intent(in) variables:
!
!       Lscrad_props   cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!       Mesorad_props  cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!       Cellrad_props  cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!     Shallowrad_props  cloud radiative properties for 
!                      clouds associated with uw shallow convection  
!                      [ microrad_properties_type ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate the elements of Lscrad_props.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        deallocate (Lscrad_props%cldext    )
        deallocate (Lscrad_props%cldsct    )
        deallocate (Lscrad_props%cldasymm  )
      endif
      deallocate (Lscrad_props%abscoeff    )

!--------------------------------------------------------------------
!    deallocate the elements of Cellrad_props and Mesorad_props.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
        deallocate (Cellrad_props%cldext   )
        deallocate (Cellrad_props%cldsct   )
        deallocate (Cellrad_props%cldasymm )
        deallocate (Cellrad_props%abscoeff )
        deallocate (Mesorad_props%cldext   )
        deallocate (Mesorad_props%cldsct   )
        deallocate (Mesorad_props%cldasymm )
        deallocate (Mesorad_props%abscoeff )
      endif

!--------------------------------------------------------------------
!    deallocate the elements of Shallowrad_props.
!---------------------------------------------------------------------
      if (Cldrad_control%do_uw_clouds) then
        deallocate (Shallowrad_props%cldext   )
        deallocate (Shallowrad_props%cldsct   )
        deallocate (Shallowrad_props%cldasymm )
        deallocate (Shallowrad_props%abscoeff )
      endif

!---------------------------------------------------------------------



end subroutine cloudrad_package_dealloc



!###################################################################





                   end module cloudrad_package_mod

