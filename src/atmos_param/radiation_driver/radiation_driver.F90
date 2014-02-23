                module radiation_driver_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
!  
! </REVIEWER>
! <OVERVIEW>
!    radiation_driver_mod is the interface between physics_driver_mod
!    and a specific radiation parameterization, currently either the
!    original_fms_rad or sea_esf_rad radiation package. it provides 
!    radiative heating rates, boundary radiative fluxes, and any other 
!    radiation package output fields to other component models of the
!    modeling system.
! </OVERVIEW>
! <DESCRIPTION>
! The following modules are called from this driver module:
!
!   1) astronomy
!
!   2) cloud properties
!
!   3) prescribed zonal ozone
!
!   4) longwave and shortwave radiation driver
! </DESCRIPTION>
!  <DIAGFIELDS>
!  Diagnostic fields may be output to a netcdf file by specifying the
!  module name radiation and the desired field names (given below)
!  in file diag_table. See the documentation for diag_manager.
!  
!  Diagnostic fields for module name: radiation
!  
!     field name      field description
!     ----------      -----------------
!  
!     alb_sfc         surface albedo (percent)
!     coszen          cosine of the solar zenith angle
!  
!     tdt_sw          temperature tendency for SW radiation (deg_K/sec)
!     tdt_lw          Temperature tendency for LW radiation (deg_K/sec)
!     swdn_toa        SW flux down at TOA (watts/m2)
!     swup_toa        SW flux up at TOA (watts/m2)
!     olr             outgoing longwave radiation (watts/m2)
!     swup_sfc        SW flux up at surface (watts/m2)
!     swdn_sfc        SW flux down at surface (watts/m2)
!     lwup_sfc        LW flux up at surface  (watts/m2)
!     lwdn_sfc        LW flux down at surface (watts/m2)
!  
!  NOTE: When namelist variable do_clear_sky_pass = .true. an additional clear sky
!        diagnostic fields may be saved.
!  
!     tdt_sw_clr      clear sky temperature tendency for SW radiation (deg_K/sec)
!     tdt_lw_clr      clear sky Temperature tendency for LW radiation (deg_K/sec)
!     swdn_toa_clr    clear sky SW flux down at TOA (watts/m2)
!     swup_toa_clr    clear sky SW flux up at TOA (watts/m2)
!     olr_clr         clear sky outgoing longwave radiation (watts/m2)
!     swup_sfc_clr    clear sky SW flux up at surface (watts/m2)
!     swdn_sfc_clr    clear sky SW flux down at surface (watts/m2)
!     lwup_sfc_clr    clear sky LW flux up at surface  (watts/m2)
!     lwdn_sfc_clr    clear sky LW flux down at surface (watts/m2)
!  </DIAGFIELDS>

! <INFO>

!   <REFERENCE>  For a specific list of radiation references see the
!     longwave and shortwave documentation.          </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>
!For some of the diagnostics fields that represent fractional amounts,
!    such as reflectivity and absorptivity, the units are incorrectly
!    given as percent.
!</BUG>
!   <NOTE> 
!CHANGE HISTORY
!changes prior to 1/24/2000
!
!  * Modified the radiation alarm. 
!    The module can now be stopped/started on a time step that is not the
!    radiation time step.
!
!  * Modified the radiation restart format.
!    Added a version number and radiation alarm information.
!
!  * Fixed a bug that occurred when namelist variable do_average = true.
!    An addition averaging variable was added for array "solar". 
!    This averaging information was also added to the restart file.
!    ***NOTE: As of this code, this namelist variable has been removed.***
!
!  * Removed the initialization for the astronomy package. This is now done
!    by the astronomy namelist.
!
!changes prior to 10/4/1999
!
!  * MPP version created. Changes to open_file and error_mesg arguments, 
!    Fortran write statements to standard output only on PE 0, Fortran close
!    statement changed to call close_file, and Fortran read/write statements
!    for restart files changed to call read_data/write_data.
!
!  * Implementation of the new MPP diagnostics package. This required major
!    changes to the diagnostic interface and the manner in which diagnostics
!    quantities are selected.
!
!  * There were no changes made that would cause answers to changes.
!
!changes prior to 5/26/1999
!
!  * added namelist variables for modifying the co2 mixing ratio.
!
!  * changed the units of namelist variable solar_constant from ly/min to watts/m2.
!
!   </NOTE>
!   <FUTURE>               </FUTURE>

! </INFO>
!   shared modules:

use mpp_mod,               only: input_nml_file
use fms_mod,               only: fms_init, mpp_clock_id, &
                                 mpp_clock_begin, mpp_clock_end, &
                                 CLOCK_MODULE,  field_exist, &
                                 field_size, &
                                 mpp_pe, mpp_root_pe, &
                                 open_namelist_file, stdlog, &
                                 file_exist, FATAL, WARNING, NOTE, &
                                 close_file, read_data, write_data, &
                                 write_version_number, check_nml_error,&
                                 error_mesg, &
                                 read_data, mpp_error
use fms_io_mod,            only: restore_state, &
                                 register_restart_field, restart_file_type, &
                                 save_restart, get_mosaic_tile_file
use diag_manager_mod,      only: register_diag_field, send_data, &
                                 diag_manager_init, get_base_time
use time_manager_mod,      only: time_type, set_date, set_time,  &
                                 get_time,    operator(+),       &
                                 print_date, time_manager_init, &
                                 assignment(=), &
                                 operator(-), operator(/=), get_date,&
                                 operator(<), operator(>=), operator(>)
use sat_vapor_pres_mod,    only: sat_vapor_pres_init, compute_qs
use constants_mod,         only: constants_init, RDGAS, RVGAS,   &
                                 STEFAN, GRAV, SECONDS_PER_DAY,  &
                                 RADIAN, diffac
use data_override_mod,     only: data_override

! shared radiation package modules:

use rad_utilities_mod,     only: radiation_control_type, Rad_control, &
                                 radiative_gases_type, &
                                 check_derived_types, &
                                 cldrad_properties_type, &
                                 astronomy_type, surface_type, &
                                 cld_specification_type, &
                                 aerosol_diagnostics_type, &
                                 atmos_input_type, rad_utilities_init,&
                                 aerosol_properties_type, aerosol_type,&
                                 sw_output_type, lw_output_type, &
                                 rad_output_type, microphysics_type, &
                                 shortwave_control_type, Sw_control, &
                                 Lw_control, &
                                 fsrad_output_type, &
                                 astronomy_inp_type, &
                                 cloudrad_control_type, Cldrad_control,&
                                 rad_utilities_end
use esfsw_parameters_mod,  only: Solar_spect, esfsw_parameters_init  

!  physics support modules:

use diag_integral_mod,     only: diag_integral_init, &
                                 diag_integral_field_init, &
                                 sum_diag_integral_field
use astronomy_mod,         only: astronomy_init, annual_mean_solar, &
                                 daily_mean_solar, diurnal_solar, &
                                 astronomy_end

!  component modules:

use original_fms_rad_mod,  only: original_fms_rad_init,  &
                                 original_fms_rad, &
                                 original_fms_rad_end
use sea_esf_rad_mod,       only: sea_esf_rad_init, sea_esf_rad, &
                                 sea_esf_rad_time_vary,  &
                                 sea_esf_rad_endts, & 
                                 sea_esf_rad_end
use rad_output_file_mod,   only: rad_output_file_init, &
                                 write_rad_output_file,    &
                                 rad_output_file_end
use cloudrad_package_mod,  only: cloudrad_package_init, &
                                 cloud_radiative_properties, &
                                 cldrad_props_dealloc, &
                                 cloudrad_package_end
use cloudrad_diagnostics_mod,      &
                           only: model_micro_dealloc, &
                                 obtain_cloud_tau_and_em, &
                                 modis_yim, modis_cmip
use microphys_rad_mod,     only: isccp_microphys_sw_driver, &
                                 isccp_microphys_lw_driver
use aerosolrad_package_mod, only: aerosolrad_package_init,    &
                                  aerosolrad_package_alloc, &
                                  aerosolrad_package_endts, &
                                  aerosolrad_package_time_vary, &
                                  aerosol_radiative_properties, &
                                  aerosolrad_package_end
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index, NO_TRACER

!--------------------------------------------------------------------

implicit none 
private 

!----------------------------------------------------------------------
!    radiation_driver_mod is the interface between physics_driver_mod
!    and a specific radiation parameterization, currently either the
!    original_fms_rad or sea_esf_rad radiation package. it provides 
!    radiative heating rates, boundary radiative fluxes, and any other 
!    radiation package output fields to other component models of the
!    modeling system.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = '$Id: radiation_driver.F90,v 20.0 2013/12/13 23:18:44 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


!---------------------------------------------------------------------
!------ interfaces -----
! <PUBLIC>
!use radiation_driver_mod [,only: radiation_driver_init,
                                   
!                                 radiation_driver,
!                                 radiation_driver_end]
!   radiation_driver_init
!       Must be called once before subroutine radiation_driver to
!       initialize the module (read namelist input and restart file).
!       Also calls the initialization routines for other modules used.
!   radiation_driver
!       Called every time step (not on the radiation time step)
!       to compute the longwave and shortwave radiative tendencies.
!   radiation_driver_end
!       Called once at the end of a model run to terminate the module (write
!       a restart file). Also calls the termination routines for other
!       modules used.
!Notes:
! 1) A namelist interface controls runtime options.
! 3) A restart file radiation_driver.res is generated by this module.
!</PUBLIC>

public    radiation_driver_init, radiation_driver, return_cosp_inputs, &
          radiation_driver_time_vary, radiation_driver_endts, &
          define_rad_times, define_atmos_input_fields,  &
          define_surface, surface_dealloc, atmos_input_dealloc, &
          microphys_dealloc, &
          radiation_driver_end, radiation_driver_restart

private  & 

! called from radiation_driver_init:
          initialize_diagnostic_integrals,   &
          diag_field_init, read_restart_nc, &

! called from radiation_driver_end:
          write_restart_nc, &

! called from radiation_driver:
          obtain_astronomy_variables, radiation_calc,    &
          update_rad_fields, produce_radiation_diagnostics,  &
          deallocate_arrays, &
          flux_trop_calc, &

! called from define_atmos_input_fields:
          calculate_auxiliary_variables


!-----------------------------------------------------------------------
!------- namelist ---------
logical :: using_restart_file = .true. ! if set to .false, restart file
                                       ! will NOT be written by this 
                                       ! module; this will not affect
                                       ! answers as long as job is 
                                       ! restarted on a radiation
                                       ! timestep
integer ::  rad_time_step = 0         !  radiative time step in seconds


integer ::  sw_rad_time_step = 0      !  radiative time step in seconds
logical :: use_single_lw_sw_ts = .true. ! lw and sw are integrated
                                        ! using rad_time_step ? if 
                                       ! false, then lw uses 
                                       ! rad_time_step, sw uses 
                                       ! sw_rad_time_step
logical ::  use_hires_coszen = .false. ! calculate for multiple zen angs
                                       ! within sw calc?
integer :: nzens_per_sw_rad_timestep = 1  !  number of cloudy
                                          ! sw calcs done on a sw rad 
                                          ! timestep
logical :: allow_nonrepro_across_restarts = .false.
                                      !  when set true, allows the 
                                      !  use_hires_coszen case to
                                      !  restart on non-radiation steps,
                                      ! with solution dependent on 
                                      ! restart interval
                                      ! (temporary until needed vari-
                                      ! ables added to restart file)
logical ::  do_clear_sky_pass= .false.!  are the clear-sky radiation
                                      !  diagnostics to be calculated ?
character(len=24) ::    &
            zenith_spec = '      '    !  string defining how zenith 
                                      !  angle is computed. acceptable
                                      !  values: 'daily_mean', 'annual_
                                      !  mean', 'diurnally_varying'
character(len=16) ::   &
                rad_package='sea_esf' !  string defining the radiation
                                      !  package being used. acceptable
                                      !  values : 'sea_esf', 
                                      !  'original_fms'     
logical ::    &
         calc_hemi_integrals = .false.!  are hemispheric integrals 
                                      !  desired ? 
logical ::     &
        all_step_diagnostics = .false.!  are lw and sw radiative bdy
                                      !  fluxes and atmospheric heating 
                                      !  rates to be output on physics 
                                      !  steps ?
logical ::     &
         renormalize_sw_fluxes=.false.!  should sw fluxes and the zenith
                                      !  angle be renormalized on each 
                                      !  timestep because of the 
                                      !  movement of earth wrt the sun ?
integer, dimension(6) ::    &
    rad_date = (/ 0, 0, 0, 0, 0, 0 /) !  fixed date for which radiation
                                      !  is to be valid (applies to
                                      !  solar info, ozone, clouds)
                                      !  [yr, mo, day, hr, min, sec]
logical  ::  &
         all_level_radiation = .true. !  is radiation to be calculated 
                                      !  at all model levels ?
integer ::    &
          topmost_radiation_level=-99 !  if all_level_radiation is 
                                      !  false., this is the lowest
                                      !  model index at which radiation
                                      !  is calculated
logical ::    &
          drop_upper_levels = .false. !  if all_level_radiation is false
                                      !  and drop_upper_levels is true,
                                      !  radiation will be calculated
                                      !  at all model levels from
                                      !  topmost_radiation_level to the
                                      !  surface
logical ::  &
         all_column_radiation = .true.!  is radiation to be calculated
                                      !  in all model columns ?
logical :: rsd=.false.                !  (repeat same day) - call 
                                      !  radiation for the specified 
                                      !  rad_date (yr,mo,day), but run 
                                      !  through the diurnal cycle (hr,
                                      !  min,sec)

logical :: use_mixing_ratio = .false. !  assumes q is mixing ratio
                                      !  rather than specific humidity
real    :: solar_constant = 1365.0    !  annual mean solar flux at top 
                                      !  of atmosphere [ W/(m**2) ]
logical :: doing_data_override = .false.  
                                      !  input fields to the radiation
                                      !  package are being overriden
                                      !  using data_override_mod ?
logical :: overriding_temps = .false. !  temperature and ts fields are
                                      !  overriden ?
logical :: overriding_sphum = .false. !  specific humidity field is
                                      !  overriden ?
logical :: overriding_clouds = .false.!  cloud specification fields are
                                      !  overriden ?
logical :: overriding_albedo = .false.!  surface albedo field is
                                      !  overriden ?
logical :: overriding_aerosol = .false.
                                      !  aerosol fields are overriden ?
logical :: use_co2_tracer_field = .false.
                                      !  obtain co2 field for use by 
                                      !  radiation package from co2
                                      !  tracer field ?
logical :: do_swaerosol_forcing = .false.
                                      !  calculating aerosol forcing in
                                      !  shortwave ?
logical :: do_lwaerosol_forcing = .false.
                                      !  calculating aerosol forcing in
                                      !  longwave ?
real    :: trop_ht_at_poles = 30000.  !  assumed height of tropoause at
                                      !  poles for case of tropause
                                      !  linearly varying with latitude
                                      !  [ Pa ]
real    :: trop_ht_at_eq    = 10000.  !  assumed height of tropoause at
                                      !  equator for case of tropause
                                      !  linearly varying with latitude
                                      !  [ Pa ]
real    :: trop_ht_constant = 20000.  !  assumed height of tropoause   
                                      !  when assumed constant       
                                      !  [ Pa ]
logical :: constant_tropo = .true.    !  generate tropopause fluxes when
                                      !  tropopause ht assumed constant?
logical :: linear_tropo   = .true.    !  generate tropopause fluxes when
                                      !  tropopause assumed to vary
                                      !  linearly with latitude?
logical :: thermo_tropo   = .false.   !  generate tropopause fluxes when
                                      !  tropopause determined thermo-
                                      !  dynamically ?
logical :: time_varying_solar_constant = .false. 
                                      !  solar_constant is to vary with
                                      !  time ?
logical :: use_uniform_solar_input = .false.
                                      !  the (lat,lon) values used to
                                      !  calculate zenith angle are
                                      !  uniform across the grid ?
real    :: lat_for_solar_input = 100. !  latitude to be used when uni-
                                      !  form solar input is activated
                                      !  [ degrees ]
real    :: lon_for_solar_input = 500. !  longitude to be used when uni-
                                      !  form solar input is activated
                                      !  [ degrees ]

logical :: always_calculate = .false. !  radiation calculation is done
                                      !  on every call to 
                                      !  radiation_driver ?
logical :: do_h2o         = .true.    !  h2o radiative effects are 
                                      !  included in the radiation 
                                      !  calculation ? 
logical :: do_o3          = .true.    !  o3 radiative effects are 
                                      !  included in the radiation 
                                      !  calculation ? 
integer, dimension(6) :: solar_dataset_entry = (/ 1, 1, 1, 0, 0, 0 /)
                                      ! time in solar data set corresp-
                                      ! onding to model initial time
                                      ! (yr, mo, dy, hr, mn, sc)
logical               :: treat_sfc_refl_dir_as_dif = .true.
                                      ! when true, solar direct  beam
                                      ! radiation reflected from the
                                      ! surface is seen as diffuse by
                                      ! the exchange grid. when false, it 
                                      ! is seen as direct, changing solar
                                      ! input to the sfc, and eliminating
                                      ! negative diffuse sw fluxes at the
                                      ! sfc, which cause problems in ESM.
real    :: solar_scale_factor = -1.0  ! factor to multiply incoming solar 
                                      ! spectral irradiances. default is to
                                      ! perform no computation. used to 
                                      ! change "solar constant"

! <NAMELIST NAME="radiation_driver_nml">
!  <DATA NAME="rad_time_step" UNITS="" TYPE="integer" DIM="" DEFAULT="14400">
!The radiative time step in seconds.
!  </DATA>
!  <DATA NAME="do_clear_sky_pass" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! are the clear-sky radiation
!  diagnostics to be calculated ?
!  </DATA>
!  <DATA NAME="zenith_spec" UNITS="" TYPE="character" DIM="" DEFAULT="">
!string defining how zenith 
!  angle is computed. acceptable
!  values: 'daily_mean', 'annual_
!  mean', 'diurnally_varying'
!  </DATA>
!  <DATA NAME="rad_package" UNITS="" TYPE="character" DIM="" DEFAULT="">
!string defining the radiation
!  package being used. acceptable
!  values : 'sea_esf', 
!  'original_fms'
!  </DATA>
!  <DATA NAME="calc_hemi_integrals" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!are hemispheric integrals 
!  desired ?
!  </DATA>
!  <DATA NAME="all_step_diagnostics" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!are lw and sw radiative bdy
!  fluxes and atmospheric heating 
!  rates to be output on physics 
!  steps ?
!  </DATA>
!  <DATA NAME="renormalize_sw_fluxes" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!should sw fluxes and the zenith
!  angle be renormalized on each 
!  timestep because of the 
!  movement of earth wrt the sun ?
!  </DATA>
!  <DATA NAME="rad_date" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!fixed date for which radiation
!  is to be valid (applies to
!  solar info, ozone, clouds)
!  [yr, mo, day, hr, min, sec]
!  </DATA>
!  <DATA NAME="all_level_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!is radiation to be calculated 
!  at all model levels ?
!  </DATA>
!  <DATA NAME="topmost_radiation_level" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!if all_level_radiation is 
!  false., this is the lowest
!  model index at which radiation
!  is calculated
!  </DATA>
!  <DATA NAME="drop_upper_levels" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!if all_level_radiation is false
!  and drop_upper_levels is true,
!  radiation will be calculated
!  at all model levels from
!  topmost_radiation_level to the
!  surface
!  </DATA>
!  <DATA NAME="all_column_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!is radiation to be calculated
!  in all model columns ?
!  </DATA>
!  <DATA NAME="rsd" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!(repeat same day) - call 
!  radiation for the specified 
!  rad_date (yr,mo,day), but run 
!  through the diurnal cycle (hr,
!  min,sec)
!  </DATA>
!  <DATA NAME="use_mixing_ratio" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!assumes q is mixing ratio
!  rather than specific humidity
!  </DATA>
!  <DATA NAME="solar_constant" UNITS="" TYPE="real" DIM="" DEFAULT="">
!annual mean solar flux at top 
!  of atmosphere [ W/(m**2) ]
!  </DATA>
!  <DATA NAME="doing_data_override" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!input fields to the radiation
!  package are being overriden
!  using data_override_mod ?
!  </DATA>
!  <DATA NAME="overriding_temps" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!temperature and ts fields are
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_sphum" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! specific humidity field is
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_clouds" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! cloud specification fields are
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_albedo" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! surface albedo field is
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_aerosol" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!aerosol fields are overriden ?
!  </DATA>
!  <DATA NAME="use_co2_tracer_field" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!use co2 value from co2 tracer field?
!  </DATA>
!  <DATA NAME="trop_ht_at_poles" UNITS="" TYPE="" DIM="" DEFAULT="">
!assumed height of tropoause at
!  poles for case of tropause
!  linearly varying with latitude
!  [ Pa ]
!  </DATA>
!  <DATA NAME="trop_ht_at_eq" UNITS="" TYPE="real" DIM="" DEFAULT="">
!assumed height of tropoause at
!  equator for case of tropause
!  linearly varying with latitude
!  [ Pa ]
!  </DATA>
!  <DATA NAME="trop_ht_constant" UNITS="" TYPE="real" DIM="" DEFAULT="">
!assumed height of tropoause   
!  when assumed constant       
!  [ Pa ]
!  </DATA>
!  <DATA NAME="constant_tropo" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! generate tropopause fluxes when
!  tropopause ht assumed constant?
!  </DATA>
!  <DATA NAME="linear_tropo" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! generate tropopause fluxes when
!  tropopause assumed to vary
!  linearly with latitude?
!  </DATA>
!  <DATA NAME="thermo_tropo" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! generate tropopause fluxes when
!  tropopause determined thermo-
!  dynamically ?
!  </DATA>
!  <DATA NAME="time_varying_solar_constant" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!solar_constant is to vary with
!  time ?
!  </DATA>
!  <DATA NAME="solar_scale_factor" UNITS="" TYPE="real" DIM="" DEFAULT="-1.">
! factor to multiply solar irradiance. default is not to perform calculation.
!  </DATA>
!  <DATA NAME="always_calculate" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  calculate radiative fluxes and heating rates on every call to 
!  radiation_driver ?
!  </DATA>
!  <DATA NAME="do_h2o" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  include h2o effects in radiation calculation ?
!  </DATA>
!  <DATA NAME="do_o3" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  include o3 effects in radiation calculation ?
!  </DATA>
!  <DATA NAME="solar_dataset_entry" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!time in solar data set corresp-
! onding to model initial time
! (yr, mo, dy, hr, mn, sc)
!  </DATA>
!  <DATA NAME="always_calculate" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!fluxes and heating rates should
! be calculatd on each call to
! radiation_driver ? (true for
! standalone applications)
!  </DATA>
!  <DATA NAME="use_uniform_solar_input" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  the (lat,lon) values used to
!  calculate zenith angle are
!  uniform across the grid ?
!  </DATA>
!  <DATA NAME="lat_for_solar_input" UNITS="" TYPE="real" DIM="" DEFAULT="">
!  latitude to be used when uni-
!  form solar input is activated
!  [ degrees ]
!  </DATA>
!  <DATA NAME="lon_for_solar_input" UNITS="" TYPE="real" DIM="" DEFAULT="">
!  longitude to be used when uni-
!  form solar input is activated
!  [ degrees ]
!  </DATA>
! </NAMELIST>
!
namelist /radiation_driver_nml/ rad_time_step, do_clear_sky_pass, &
                                using_restart_file, &
                                sw_rad_time_step,  &
                                use_single_lw_sw_ts, &
                                use_hires_coszen, &
                                allow_nonrepro_across_restarts, &
                                nzens_per_sw_rad_timestep, &
                                zenith_spec, rad_package,    &
                                calc_hemi_integrals,     &
                                all_step_diagnostics, &
                                renormalize_sw_fluxes, &
                                rad_date, all_level_radiation, &
                                topmost_radiation_level,   &
                                drop_upper_levels,  &
                                all_column_radiation, rsd,    &
                                use_mixing_ratio, solar_constant, &
                                doing_data_override, &
                                overriding_temps, overriding_sphum, &
                                overriding_clouds, overriding_albedo, &
                                overriding_aerosol, &
                                use_co2_tracer_field, &
                                do_swaerosol_forcing,   &
                                do_lwaerosol_forcing, &
                                trop_ht_at_poles, trop_ht_at_eq, &
                                trop_ht_constant, constant_tropo, &
                                linear_tropo, thermo_tropo, &
                                time_varying_solar_constant, &
                                solar_dataset_entry, &
                                solar_scale_factor,  &
                                treat_sfc_refl_dir_as_dif, &
                                always_calculate,  do_h2o, do_o3, &
                                use_uniform_solar_input, &
                                lat_for_solar_input, lon_for_solar_input
!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----
!-- for netcdf restart
type(restart_file_type), pointer, save :: Rad_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
integer                                :: int_renormalize_sw_fluxes
integer                                :: int_do_clear_sky_pass


!---------------------------------------------------------------------
!    logical  flags.

logical ::  module_is_initialized = .false. ! module initialized?
logical ::  do_rad                          ! is this a radiation step ?
logical ::  do_lw_rad, do_sw_rad            ! is this a radiation step ?
logical ::  use_rad_date                    ! specify time of radiation
                                            ! independent of model time?
logical ::  do_sea_esf_rad                  ! using sea_esf_rad package?

!---------------------------------------------------------------------
!    list of restart files readable by this module.
!
!                 sea_esf_rad.res:
!
!     version 1:  sea_esf_rad.res file version used initially in 
!                 AM2 model series (through galway code, AM2p8). this
!                 is the only version of sea_esf_rad.res ever produced.
!
!                 radiation_driver.res:
!
!     version 1:  not readable by this module.
!     version 2:  added cosine of zenith angle as an output to
!                 radiation_driver.res  (6/27/00)
!     version 3:  added restart variables needed when sw renormalization
!                 is active. (3/21/02)
!     version 4:  added longwave heating rate as separate output 
!                 variable, since it is needed as input to edt_mod
!                 and entrain_mod. (7/17/02)
!     version 5:  removed variables associated with the former 
!                 do_average namelist option (7/23/03)
!     version 6:  added writing of sw tropospheric fluxes (up and
!                 down) so that they are available for the renormal-
!                 ization case (developed by ds, 10/03; added to 
!                 trunk code 01/14/04).
!     version 7:  added swdn to saved variables (developed by slm 
!                 11/23/03, added to trunk code 01/14/04).
!     version 8:  includes additional sw fluxes at sfc, used with
!                 land model (11/13/03).
!     version 9:  consolidation of version 6 and version 8. (version 7
!                 replaced by version 8.)
!     version 10: adds 2 clr sky sw down diffuse and direct sfc flux
!                 diagnostic variables (10/18/04)
!     version 11: adds flux_sw_down_vis_clr diagnostic variable for use
!                 in assessing polar ice maintainability (6/19/07)
!---------------------------------------------------------------------
integer, dimension(10) :: restart_versions     = (/ 2, 3, 4, 5, 6,  &
                                                   7, 8, 9, 10, 11 /)
integer                :: vers ! version number of the restart file being read

!-----------------------------------------------------------------------
!    these arrays must be preserved across timesteps:
!
!    Rad_output is a rad_output_type variable with the following 
!    components:
!          tdt_rad        radiative (sw + lw) heating rate
!          flux_sw_surf   net (down-up) sw flux at surface
!          flux_sw_surf_dir   net (down-up) sw flux at surface
!          flux_sw_surf_refl_dir   dir sw flux reflected at surface
!          flux_sw_surf_dif   net (down-up) sw flux at surface
!          flux_sw_down_vis_dir  downward visible sw flux at surface
!          flux_sw_down_vis_dif  downward visible sw flux at surface
!          flux_sw_down_total_dir  downward total sw flux at surface
!          flux_sw_down_total_dif  downward total sw flux at surface
!          flux_sw_down_total_dir_clr  downward total direct sw flux at 
!                                      surface  (clear sky)
!          flux_sw_down_total_dif_clr  downward total diffuse sw flux 
!                                      at surface   (clear sky)
!          flux_sw_down_vis_clr  downward visible sw flux at surface
!                                       (clear sky)
!          flux_sw_vis    net visible sw flux at surface
!          flux_sw_vis_dir    net visible sw flux at surface
!          flux_sw_refl_vis_dir reflected direct visible sw flux at surface
!          flux_sw_vis_dif net visible sw flux at surface
!          flux_lw_surf   downward lw flux at surface
!          coszen_angle   cosine of the zenith angle (used for the 
!                         last radiation calculation)
!          tdt_rad_clr    net radiative heating rate in the absence of
!                         cloud
!          tdtsw          shortwave heating rate
!          tdtsw_clr      shortwave heating rate in he absence of cloud
!          tdtlw_clr       longwave heating rate in he absence of cloud
!          tdtlw          longwave heating rate
!          ufsw          upward sw flux
!          dfsw          downward sw flux
!          ufsw_clr      upward sw flux
!          dfsw_clr      downward sw flux
!          flxnet        net lw flux
!          flxnetcf      net lw flux, cloud free

!    solar_save is used when renormalize_sw_fluxes is active, to save
!    the solar factor (fracday*cosz/r**2) from the previous radiation
!    step so that the radiative forcing terms may be adjusted on each
!    timestep to reflect the current solar forcing.
!
!    sw_heating_clr, tot_heating_clr_save, sw_heating_save, 
!    tot_heating_save, flux_sw_surf_save, flux_sw_surf_dir_save,
!    flux_sw_surf_refl_dir_save,  flux_sw_refl_vis_dir_save, &
!    flux_sw_surf_dif_save, flux_sw_down_vis_dir_save, 
!    flux_sw_down_vis_dif_save, flux_sw_down_vis_clr_save,
!    flux_sw_down_total_dir_clr_save, flux_sw_down_total_dif_clr_save,
!    flux_sw_down_total_dir_save, flux_sw_down_total_dif_save and 
!    flux_sw_vis_save, flux_sw_vis_dir_save, flux_sw_vis_dif_save are 
!    the radiative forcing terms on radiation steps which also must be 
!    saved when renormalization is activated.

!    swdn_special_save, swup_special_save, swdn_special_clr_save,
!    swup_special_clr_save are also saved.
!
!    the ***sw_save arrays are currently saved so that their values may
!    be adjusted during sw renormalization for diagnostic purposes.
!                               
!    the **lw_save arrays are currently saved so that they may be output
!    in the diagnostics file on every physics step, if desired, so that
!    when renormalize_sw_fluxes is active, total radiative terms may be
!    easily generated.
!-----------------------------------------------------------------------

type(rad_output_type),save          ::  Rad_output
real, allocatable, dimension(:,:)   ::  solar_save, &
                                        dum_idjd
real, allocatable, dimension(:,:,:)   ::  &
                                    flux_sw_down_total_dir_clr_save, &
                                    flux_sw_down_total_dif_clr_save, &
                                        flux_sw_down_vis_clr_save
real, allocatable, dimension(:,:,:)   ::   flux_sw_surf_save, &
                                        flux_sw_surf_dir_save, &
                                        flux_sw_surf_refl_dir_save, &
                                        flux_sw_surf_dif_save, &
                                        flux_sw_down_vis_dir_save, &
                                        flux_sw_down_vis_dif_save, &
                                        flux_sw_down_total_dir_save, &
                                        flux_sw_down_total_dif_save, &
                                        flux_sw_vis_save, &
                                        flux_sw_vis_dir_save, &
                                        flux_sw_refl_vis_dir_save, &
                                        flux_sw_vis_dif_save
real, allocatable, dimension(:,:,:,:) ::  sw_heating_save,    &
                                        tot_heating_save, &
                                        dfsw_save, ufsw_save, fsw_save,&
                                        hsw_save
real, allocatable, dimension(:,:,:,:) ::  sw_heating_clr_save, &
                                        tot_heating_clr_save, &
                                        dfswcf_save,   &
                                        ufswcf_save, fswcf_save, &
                                        hswcf_save
real, allocatable, dimension(:,:,:) ::  tdtlw_save, tdtlw_clr_save
 real, allocatable, dimension(:,:,:) ::  flxnet_save, flxnetcf_save
real, allocatable, dimension(:,:)   ::  olr_save, lwups_save, &
                                        lwdns_save, olr_clr_save, &
                                        lwups_clr_save, lwdns_clr_save
real, allocatable, dimension(:,:,:,:) ::  swdn_special_save, &
                                        swdn_special_clr_save, &  
                                        swup_special_save,&
                                        swup_special_clr_save
real, allocatable, dimension(:,:,:) ::  netlw_special_save, &
                                        netlw_special_clr_save
real, allocatable, dimension(:,:,:,:) ::  dfsw_ad_save, ufsw_ad_save
real, allocatable, dimension(:,:,:,:) ::  dfswcf_ad_save, ufswcf_ad_save
real, allocatable, dimension(:,:)   ::  olr_ad_save, lwups_ad_save, &
                                        lwdns_ad_save, olr_ad_clr_save, &
                                    lwups_ad_clr_save, lwdns_ad_clr_save

!-----------------------------------------------------------------------
!    time-step-related constants
 
integer    :: lwrad_alarm    !  time interval until the next radiation 
                             !  calculation (seconds)
integer    :: swrad_alarm    !  time interval until the next radiation 
                             !  calculation (seconds)
integer    :: current_sw_zenith_step = 1  
                             !  current zenith angle index being used  
                             !  for cloudy sw calculations when 
                             !  use_hires_coszen is .true.
integer    :: num_pts=0      !  counter for current number of grid 
                             !  columns processed (when num_pts=0 or 
                             !  num_pts=total_pts certain things happen)
integer    :: total_pts      !  number of grid columns to be processed 
                             !  every time step (note: all grid columns
                             !  must be processed every time step)
type(time_type) :: Rad_time  !  time at which the climatologically-
                             !  determined, time-varying input fields to
                             !  radiation should apply 
                             !  [ time_type (days, seconds)]
integer    :: dt             !  physics time step (frequency of calling 
                             !  radiation_driver)  [ seconds ]
integer   :: lw_rad_time_step

!-----------------------------------------------------------------------
!    diagnostics variables
integer, parameter :: MX_SPEC_LEVS = 4 
                             ! number of special levels at
                             ! which radiative fluxes are to be 
                             ! calculated for diagnostic purposes

character(len=16)            :: mod_name = 'radiation'
integer                      :: id_alb_sfc, id_cosz, id_fracday, &
                                id_alb_sfc_avg, &
                                id_alb_sfc_vis_dir, id_alb_sfc_nir_dir,&
                                id_alb_sfc_vis_dif, id_alb_sfc_nir_dif
integer                      :: id_flux_sw_dir, id_flux_sw_dif, &
                                id_flux_sw_refl_dir,  &
                                id_flux_sw_refl_vis_dir, id_flux_sw, &
                                id_flux_sw_down_vis_dir, &
                                id_flux_sw_down_vis_dif, &
                                id_flux_sw_down_total_dir, &
                                id_flux_sw_down_total_dif, &
                                id_flux_sw_down_total_dir_clr, &
                                id_flux_sw_down_total_dif_clr, &
                                id_flux_sw_down_vis_clr, &
                                id_flux_sw_vis, &
                                id_flux_sw_vis_dir, &
                                id_flux_sw_vis_dif, &
                                id_rrvco2, id_rrvf11, id_rrvf12, &
                                id_rrvf113, id_rrvf22, id_rrvch4, &
                                id_rrvn2o, id_co2_tf, id_ch4_tf, &
                                id_n2o_tf, id_sol_con
integer                      :: id_conc_drop, id_conc_ice

integer                      :: id_allradp
integer, dimension(2)        :: id_tdt_sw,   id_tdt_lw,  &
                                id_ufsw, id_dfsw,  &
                                id_flxnet, &
                                id_swdn_toa, id_swup_toa, id_olr, &
                                id_netrad_toa,  id_netrad_1_Pa,  &
                                id_swup_sfc, id_swdn_sfc,         &
                                id_lwup_sfc, id_lwdn_sfc
integer, dimension(MX_SPEC_LEVS,2)   :: id_swdn_special,   &
                                        id_swup_special,  &
                                        id_netlw_special
integer, dimension(2)        :: id_swtoa, id_swsfc,               &
                                id_lwsfc,                         &
                                id_swtoa_ad, id_swsfc_ad,         &
                                id_swdn_sfc_ad,                   &
                                id_swup_sfc_ad,                   &
                                id_swup_toa_ad,                   &
                                id_olr_ad, id_lwsfc_ad



real                         :: missing_value = -999.
character(len=8)             :: std_digits   = 'f8.3'
character(len=8)             :: extra_digits = 'f16.11'

!-----------------------------------------------------------------------
!    timing clocks       

integer                      :: misc_clock, clouds_clock, calc_clock

!--------------------------------------------------------------------
! miscellaneous variables and indices

integer        ::  ks         !  model grid coordinate of top level
                              !  at which radiation is calculated 
                              !  (topmost_radiation_level)
integer        ::  ke         !  model grid coordinate of bottommost
                              !  level at which radiation is calculated

integer        ::  ksrad=1    !  always set to 1
integer        ::  kerad      !  number of layers in radiation grid

real           ::  rh2o_lower_limit_orig=3.0E-06
                              !  smallest value of h2o mixing ratio 
                              !  allowed with original_fms_rad package
real           ::  rh2o_lower_limit_seaesf=2.0E-07
                              !  smallest value of h2o mixing ratio 
                              !  allowed with sea_esf_rad package
real           ::  rh2o_lower_limit
                              !  smallest value of h2o mixing ratio 
                              !  allowed in the current experiment
real           ::  temp_lower_limit=100.0  ! [ K ]
                              !  smallest value of temperature      
                              !  allowed in the current experiment
real           ::  temp_upper_limit=370.00  ! [ K ]
                              !  largest value of temperature 
                              !  allowed in the current experiment

real           ::  surf_flx_init=50.0  ! [w / m^2 ]
                              !  value to which surface lw and sw fluxes
                              !  are set in the absence of a .res file
                              !  containing them

real           ::  coszen_angle_init=0.50
                              !  value to which cosine of zenith angle  
                              !  is set in the absence of a .res file
                              !  containing it

real           ::  log_p_at_top=2.0
                              !  assumed value of ln of ratio of pres-
                              !  sure at flux level 2 to that at model
                              !  top (needed for deltaz calculation,
                              !  is infinite for model top at p = 0.0,
                              !  this value is used to give a reasonable
                              !  deltaz)
real,parameter ::  D608 = (RVGAS-RDGAS)/RDGAS
                              !  virtual temperature factor  
real,parameter ::  D622 = RDGAS/RVGAS
                              ! ratio of gas constants - dry air to 
                              ! water vapor
real,parameter ::  D378 = 1.0 - D622  
                              ! 1 - gas constant ratio
integer :: id, jd
integer        ::  size_of_lwoutput = 1
integer        ::  size_of_swoutput = 1
integer        ::  indx_lwaf = 0
integer        ::  indx_swaf = 0

real, dimension(:,:), allocatable :: solflxtot_lean
real            :: solflxtot_lean_ann_1882, solflxtot_lean_ann_2000
integer         ::   first_yr_lean, last_yr_lean,   &
                     nvalues_per_year_lean, numbands_lean
integer         ::   years_of_data_lean 

type(time_type) :: Model_init_time, Solar_offset, &
                   Solar_entry
logical         :: negative_offset = .false.
real, dimension(:,:), allocatable :: swups_acc, swdns_acc
real, dimension(:,:), allocatable :: olr_intgl, swabs_intgl

! <DATASET NAME="CO2 transmission functions">
! Several ascii files are required that can be easily setup by a
!     script for getting physics data sets.
! </DATASET>
! <DATASET NAME="Restart file">
! A restart data set called radiation_driver.res(.nc) saves the
!     global fields for the current radiative tendency, net shortwave
!     surface flux, downward longwave surface flux, and cosine of the
!     zenith angle. If the namelist variable do_average=true,
!     then additional time averaged global data is written.
!     If the restart file is not present when initializing then the
!     radiative tendency is set to zero, the SW and LW surface fluxes
!     to 50 watts/m2, and the cosine of the zenith angle to 0.50.
!     Since radiation is usually computed on the first time step when
!     restarting, these values may have little or no effect.  If the
!     restart file is not present time average data is also set to zero.
! </DATASET>
!<REFERENCE> For a specific list of radiation references see the
!     longwave and shortwave documentation.</REFERENCE>
!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################
! <SUBROUTINE NAME="radiation_driver_init">
!  <OVERVIEW>
!   radiation_driver_init is the constructor for radiation_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_driver_init is the constructor for radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver_init (lonb, latb, pref, axes, Time, &
!                                  aerosol_names)
!
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
!    lonb      Longitude in radians for all (i.e., the global size)
!              grid box corners, the size of lonb should be one more
!              than the number of points along the x-axis and y-axis.
!                 [real, dimension(:,:)]
!  </IN>
!  <IN NAME="latb" TYPE="real">
!    latb      Latitude in radians for all (i.e., the global size)
!              grid box corners, the size of latb should be one more
!              than the number of latitude points along the x-axis and y-axis.
!                 [real, dimension(:,:)]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!    pref      Two reference profiles of pressure at full model levels
!              plus the surface (nlev+1). The first profile assumes a surface
!              pressure of 101325 pa, and the second profile assumes 
!              81060 pa.  [real, dimension(nlev+1,2)]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!    axes      The axis indices that are returned by previous calls to
!              diag_axis_init. The values of this array correspond to the
!              x, y, full (p)level, and half (p)level axes. These are the
!              axes that diagnostic fields are output on.
!                 [integer, dimension(4)]
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    Time      The current time.  [time_type]
!  </IN>
!  <IN NAME="aerosol_names" TYPE="character">
!   Aerosol names
!  </IN>
!   <ERROR MSG="must have two reference pressure profile" STATUS="FATAL">
!     The input argument pref must have a second dimension size of 2.
!   </ERROR>
!   <ERROR MSG="restart version ## cannot be read by this module version" STATUS="FATAL">
!     You have attempted to read a radiation_driver.res file with either
!       no restart version number or an incorrect restart version number.
!   </ERROR>
!
!   <NOTE>
!    radiation time step has changed, next radiation time also changed
!       The radiation time step from the namelist input did not match
!       the radiation time step from the radiation restart file.
!       The next time for radiation will be adjusted for the new  namelist
!       input) value.
!   </NOTE>
! </SUBROUTINE>
!
subroutine radiation_driver_init (lonb, latb, pref, axes, Time, &
                                  donner_meso_is_largescale, &
                                  aerosol_names, aerosol_family_names,&
                                  do_cosp, ncol)

!---------------------------------------------------------------------
!   radiation_driver_init is the constructor for radiation_driver_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real, dimension(:,:),            intent(in)  :: lonb, latb
real, dimension(:,:),            intent(in)  :: pref
integer, dimension(4),           intent(in)  :: axes
type(time_type),                 intent(in)  :: Time
logical,                         intent(in)  :: donner_meso_is_largescale
character(len=*), dimension(:), intent(in)   :: aerosol_names
character(len=*), dimension(:), intent(in)   :: aerosol_family_names
logical,                         intent(in)  :: do_cosp
integer,                         intent(out) :: ncol
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb           2d array of model longitudes on cell corners 
!                      [ radians ]
!       latb           2d array of model latitudes at cell corners 
!                      [ radians ]
!       pref           array containing two reference pressure profiles 
!                      for use in defining transmission functions
!                      [ pascals ]
!       axes           diagnostic variable axes
!       Time           current time [time_type(days, seconds)]
!       aerosol_names  names associated with the activated aerosol
!                      species
!       aerosol_family_names  
!                      names associated with the activated aerosol
!                      families
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      integer           ::   unit, io, ierr, logunit
      integer           ::   kmax 
      integer           ::   nyr, nv, nband
      integer           ::   yr, month, year, dum
      integer           ::   ico2
      integer           ::   nzens

!---------------------------------------------------------------------
!   local variables
! 
!        unit    io unit number for namelist file
!        io      error status returned from io operation
!        ierr    error code
!        id      number of grid points in x direction (on processor)
!        jd      number of grid points in y direction (on processor)
!        kmax    number of model layers
!                
!---------------------------------------------------------------------

      
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized. note that data_override_init cannot
!    be called successfully from here (a data_override_mod feature);
!    instead it relies upon a check for previous initialization when
!    subroutine data_override is called.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call diag_manager_init
      call time_manager_init
      call sat_vapor_pres_init
      call constants_init
      call diag_integral_init
      call esfsw_parameters_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=radiation_driver_nml, iostat=io)
      ierr = check_nml_error(io,'radiation_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=radiation_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'radiation_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!--------------------------------------------------------------------
!    make sure other namelist variables are consistent with 
!    doing_data_override. Validate here to prevent potentially mis-
!    leading values from going into the stdlog file.
!--------------------------------------------------------------------
      if (.not. doing_data_override) then
        overriding_temps   = .false.
        overriding_sphum   = .false.
        overriding_albedo  = .false.
        overriding_clouds  = .false.
        overriding_aerosol = .false.
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=radiation_driver_nml)

!---------------------------------------------------------------------
!    set logical variable defining the radiation scheme desired from the
!    namelist-input character string. set lower limit to water vapor 
!    mixing ratio that the radiation code will see, to assure keeping 
!    within radiation lookup tables. exit if value is invalid.
!---------------------------------------------------------------------
      if (rad_package == 'original_fms') then
        do_sea_esf_rad = .false.
        rh2o_lower_limit = rh2o_lower_limit_orig
      else if (rad_package == 'sea_esf') then
        do_sea_esf_rad = .true.
        rh2o_lower_limit = rh2o_lower_limit_seaesf
      else
        call error_mesg ('radiation_driver_mod',  &
           'string provided for rad_package is not valid', FATAL)
      endif

!---------------------------------------------------------------------
!    set control variable indicating whether water vapor effects are to
!    be included in the radiative calculation. if h2o effects are not 
!    to be included in the radiative calculations, set the lower limit 
!    for h2o to zero. set flag to indicate do_h2o has been initialized.
!---------------------------------------------------------------------
      Lw_control%do_h2o = do_h2o 
      if (.not. do_h2o) then
        rh2o_lower_limit = 0.0
      endif
      Lw_control%do_h2o_iz = .true.

!---------------------------------------------------------------------
!    set control variable indicating whether ozone effects are to be
!    included in the radiative calculation. set flag to indicate the
!    control variable has been initialized.
!---------------------------------------------------------------------
      Lw_control%do_o3 = do_o3 
      Lw_control%do_o3_iz = .true.
      
!---------------------------------------------------------------------
!    stop execution if overriding of aerosol data has been requested.
!    code to do so has not yet been written.
!---------------------------------------------------------------------
      if (overriding_aerosol) then
        call error_mesg ('radiation_driver_mod', &
                'overriding of aerosol data not yet implemented', FATAL)
      endif

!RSH:
!RSH    if use_co2_tracer_field is .true., verify here that there is
!RSH   in fact a co2 field included within the tracer array. if not,
!RSH   call error_mesg and abort execution.
!RSH
      if(use_co2_tracer_field) then
         ico2 = get_tracer_index(MODEL_ATMOS, 'co2')
         if(ico2 == NO_TRACER) then
            call error_mesg('radiation_driver_mod', &
                 'co2 must be present as a tracer when use_co2_tracer_field is .true.', FATAL)
         endif
      endif

!--------------------------------------------------------------------
!    set logical variables defining how the solar zenith angle is to
!    be  defined from the namelist-input character string.  exit if the
!    character string is invalid.
!--------------------------------------------------------------------
      if (zenith_spec == 'diurnally_varying') then
        Sw_control%do_diurnal = .true.
        Sw_control%do_annual = .false.
        Sw_control%do_daily_mean = .false.
      else if (zenith_spec == 'daily_mean') then
        Sw_control%do_diurnal = .false.
        Sw_control%do_annual = .false.
        Sw_control%do_daily_mean = .true.
      else if (zenith_spec == 'annual_mean') then
        Sw_control%do_diurnal = .false.
        Sw_control%do_annual = .true.
        Sw_control%do_daily_mean = .false.
      else
        call error_mesg ('radiation_driver_mod', &    
            'string provided for zenith_spec is invalid', FATAL)
      endif

!--------------------------------------------------------------------
!    check if spacially-uniform solar input has been requested. if it
!    has, verify that the requested lat and lon are valid, and convert
!    them to radians.
!--------------------------------------------------------------------
      if (use_uniform_solar_input) then
        if (lat_for_solar_input < -90. .or. &
            lat_for_solar_input >  90. ) then
          call error_mesg ('radiation_driver_mod', &
            'specified latitude for uniform solar input is invalid', &
                                                            FATAL)
        else
          lat_for_solar_input = lat_for_solar_input/RADIAN
        endif
        if (lon_for_solar_input < 0. .or. &
            lon_for_solar_input > 360. ) then
          call error_mesg ('radiation_driver_mod', &
            'specified longitude for uniform solar input is invalid', &
                                                             FATAL)
        else
          lon_for_solar_input = lon_for_solar_input/RADIAN
        endif
      endif

!--------------------------------------------------------------------
!     code to handle time-varying solar input
!--------------------------------------------------------------------
        if (file_exist('INPUT/lean_solar_spectral_data.dat')) then
          unit = open_namelist_file   &
                                 ('INPUT/lean_solar_spectral_data.dat')
          read (unit, FMT = '(4i8)') first_yr_lean, last_yr_lean,  &
                                   nvalues_per_year_lean, numbands_lean
          if (numbands_lean /= Solar_spect%nbands) then
            call error_mesg ('radiation_driver_mod', &
            ' number of sw parameterization bands in solar_spectral &
            &data file differs from that defined in esfsw_parameters',&
                                                           FATAL)
          endif
          years_of_data_lean = last_yr_lean - first_yr_lean + 1
          allocate (solflxtot_lean   &
                           (years_of_data_lean, nvalues_per_year_lean))
          allocate (Solar_spect%solflxband_lean   &
             (years_of_data_lean, nvalues_per_year_lean, numbands_lean))
          allocate (Solar_spect%solflxband_lean_ann_1882(numbands_lean))
          read (unit, FMT = '(2i6,f17.4)') yr, month, &
                                          solflxtot_lean_ann_1882
          read (unit, FMT = '(6e12.5 )')   &
                 (Solar_spect%solflxband_lean_ann_1882 &
                                 (nband), nband =1,numbands_lean)
          if (solar_scale_factor >= 0.0) then
            solflxtot_lean_ann_1882 = solflxtot_lean_ann_1882*  &
                                                        solar_scale_factor
            do nband = 1,numbands_lean
              Solar_spect%solflxband_lean_ann_1882(nband) =    &
                          Solar_spect%solflxband_lean_ann_1882(nband)* &
                                                         solar_scale_factor
            enddo
          endif
          do nyr=1,years_of_data_lean
            do nv=1,nvalues_per_year_lean
              read (unit, FMT = '(2i6,f17.4)') yr, month, &
                                       solflxtot_lean(nyr,nv)
              read (unit, FMT = '(6e12.5 )')   &
                 (Solar_spect%solflxband_lean  &
                                (nyr,nv,nband), nband =1,numbands_lean)
            end do
          end do
          if (solar_scale_factor >= 0.0) then
            do nv=1,nvalues_per_year_lean
              do nyr=1,years_of_data_lean
                solflxtot_lean(nyr,nv) =    &
                         solflxtot_lean(nyr,nv)*solar_scale_factor
              enddo
            enddo
            do nband = 1,numbands_lean
              do nv=1,nvalues_per_year_lean
                do nyr=1,years_of_data_lean
                  Solar_spect%solflxband_lean(nyr,nv,nband) =    &
                           Solar_spect%solflxband_lean(nyr,nv,nband)*   &
                                                        solar_scale_factor
                enddo
              enddo
            enddo
          endif
          allocate (Solar_spect%solflxband_lean_ann_2000(numbands_lean))
          read (unit, FMT = '(2i6,f17.4)') yr, month, &
                                           solflxtot_lean_ann_2000
          read (unit, FMT = '(6e12.5 )')   &
              (Solar_spect%solflxband_lean_ann_2000 &
                              (nband), nband =1,numbands_lean)
          call close_file (unit) 
          if (solar_scale_factor >= 0.0) then
            solflxtot_lean_ann_2000 = solflxtot_lean_ann_2000*  &
                                                       solar_scale_factor
            do nband = 1,numbands_lean
              Solar_spect%solflxband_lean_ann_2000(nband) =    &
                        Solar_spect%solflxband_lean_ann_2000(nband)* &
                                                        solar_scale_factor
            enddo
          endif
        else
          if (time_varying_solar_constant) then
            call error_mesg ('radiation_driver_mod', &
             'desired solar_spectral_data input file is not present', &
                                                             FATAL)
          endif
        endif
        
        if (time_varying_solar_constant) then
!----------------------------------------------------------------------
!    define the model base time.
!----------------------------------------------------------------------
          Model_init_time = get_base_time()

!----------------------------------------------------------------------
!    if no solar_dataset_entry is supplied, use the model base time,
!    meaning that the timeseries data will be mapped to the model time
!    without any offset.
!----------------------------------------------------------------------
        if (solar_dataset_entry(1) == 1 .and. &
            solar_dataset_entry(2) == 1 .and. &
            solar_dataset_entry(3) == 1 .and. &
            solar_dataset_entry(4) == 0 .and. &
            solar_dataset_entry(5) == 0 .and. &
            solar_dataset_entry(6) == 0 ) then
          Solar_entry = Model_init_time
 
!----------------------------------------------------------------------
!    if a solar_dataset_entry is supplied, define a corresponding
!    time-type variable.
!----------------------------------------------------------------------
        else
          Solar_entry  = set_date (solar_dataset_entry(1), &
                                   solar_dataset_entry(2), &
                                   solar_dataset_entry(3), &
                                   solar_dataset_entry(4), &
                                   solar_dataset_entry(5), &
                                   solar_dataset_entry(6))
        endif

        call error_mesg ('radiation_driver_mod', &
             'Solar data is varying in time', NOTE)
        call print_date (Solar_entry , str='Data from solar timeseries &
                                            &at time:')
        call print_date (Model_init_time , str='This data is mapped to &
                                             &model time:')
        Solar_offset = Solar_entry - Model_init_time
 
        if (Model_init_time > Solar_entry) then
          negative_offset = .true.
        else
          negative_offset = .false.
        endif
        Rad_control%using_solar_timeseries_data = .true.
        Rad_control%using_solar_timeseries_data_iz = .true.

!---------------------------------------------------------------------
!    if solar input not time-varying, define solar constant and set 
!    offset to 0.0.
!---------------------------------------------------------------------
      else
        if (solar_dataset_entry(1) == 1 .and. &
          solar_dataset_entry(2) == 1 .and. &
          solar_dataset_entry(3) == 1 .and. &
          solar_dataset_entry(4) == 0 .and. &
          solar_dataset_entry(5) == 0 .and. &
          solar_dataset_entry(6) == 0 ) then
          Sw_control%solar_constant = solar_constant
          Solar_offset = set_time(0,0)
          call error_mesg ('radiation_driver_mod', &
                   'Solar data is fixed in time at nml value', NOTE)
          Rad_control%using_solar_timeseries_data = .false.
          Rad_control%using_solar_timeseries_data_iz = .true.
        else
 
!----------------------------------------------------------------------
!    convert solar_dataset_entry to a time_type variable.
!----------------------------------------------------------------------
          Solar_entry  = set_date (solar_dataset_entry(1), &
                                   solar_dataset_entry(2), &
                                   solar_dataset_entry(3), &
                                   solar_dataset_entry(4), &
                                   solar_dataset_entry(5), &
                                   solar_dataset_entry(6))
          call error_mesg ('radiation_driver_mod', &
                                'Solar data is fixed in time', NOTE)
          call print_date (Solar_entry ,    &
             str='Data used in this experiment is from solar &
                  &timeseries at time:')
          if (size(Solar_spect%solflxband(:)) /= numbands_lean) then
            call error_mesg ('radiation_driver_mod', &
             'bands present in solar constant time data differs from &
               &model parameterization band number', FATAL)
          endif

!--------------------------------------------------------------------
!    define time to be used for solar input data.
!--------------------------------------------------------------------
          call get_date (Solar_entry, year, month, dum, dum, dum, dum)

!--------------------------------------------------------------------
!    define input value based on year and month of Solar_time.
!--------------------------------------------------------------------
          if (year < first_yr_lean) then
            Sw_control%solar_constant = solflxtot_lean_ann_1882
            do nband=1,numbands_lean
              Solar_spect%solflxband(nband) =  &
                       Solar_spect%solflxband_lean_ann_1882(nband)
            end do
          else if (year > last_yr_lean) then
            Sw_control%solar_constant = solflxtot_lean_ann_2000    
            do nband=1,numbands_lean
              Solar_spect%solflxband(nband) =  &
                  Solar_spect%solflxband_lean_ann_2000(nband)
            end do
          else
            Sw_control%solar_constant =   & 
                            solflxtot_lean(year-first_yr_lean+1, month)
            do nband=1,numbands_lean
              Solar_spect%solflxband(nband) =  &
         Solar_spect%solflxband_lean(year-first_yr_lean+1, month, nband)
            end do
          endif
          Rad_control%using_solar_timeseries_data = .true.
          Rad_control%using_solar_timeseries_data_iz = .true.
        endif
      endif

!---------------------------------------------------------------------
!     include logical control in Rad_control derived-type variable.
!---------------------------------------------------------------------
      Rad_control%time_varying_solar_constant =  &
                                        time_varying_solar_constant
      Rad_control%time_varying_solar_constant_iz = .true.

!---------------------------------------------------------------------
!    set flags indicating that the Sw_control variables have been 
!    defined.
!---------------------------------------------------------------------
      Sw_control%do_diurnal_iz = .true.
      Sw_control%do_annual_iz = .true.
      Sw_control%do_daily_mean_iz = .true.

!---------------------------------------------------------------------
!    be sure that sw renormalization and hi-res zenith angle are not
!    both selected as options. they are mutually exclusive.
!---------------------------------------------------------------------
      if (renormalize_sw_fluxes .and. use_hires_coszen) then
        call error_mesg ('radiation_driver_init', &
         ' cannot select both hi-res zenith angle and sw &
              &renormalization at same time -- choose only one', FATAL)
      endif

!---------------------------------------------------------------------
!    verify that radiation has been requested at all model levels and in
!    all model columns when the original fms radiation is activated.    
!    verify that renormalize_sw_fluxes has not been requested along
!    with the original fms radiation package. verify that all_step_diag-
!    nostics has not been requested with the original fms radiation
!    package.
!---------------------------------------------------------------------
      if (.not. do_sea_esf_rad) then
        if (.not. all_level_radiation .or. &
            .not. all_column_radiation) then
          call error_mesg ( 'radiation_driver_mod', &
        ' must specify all_level_radiation and all_column_radiation'//&
            ' as true when using original fms radiation', FATAL)
        endif
        if (renormalize_sw_fluxes) then
          call error_mesg ( 'radiation_driver_mod', &
           ' cannot renormalize shortwave fluxes with original_fms '//&
                 'radiation package.', FATAL)
        endif
        if (all_step_diagnostics) then
          call error_mesg ( 'radiation_driver_mod', &
            ' cannot request all_step_diagnostics with original_fms '//&
              'radiation package.', FATAL)
        endif
      endif

!---------------------------------------------------------------------
!    can only renormalize shortwave fluxes when diurnally_varying
!    radiation is used.
!---------------------------------------------------------------------
     if (renormalize_sw_fluxes .and. .not. Sw_control%do_diurnal) then
       call error_mesg ('radiation_driver_mod',  &
       ' can only renormalize sw fluxes when using diurnally-varying'//&
                       ' solar radiation', FATAL)
     endif

!---------------------------------------------------------------------
!    verify that a valid radiation time step has been specified.
!---------------------------------------------------------------------
      if (rad_time_step <= 0) then
        call error_mesg ('radiation_driver_mod', &
            ' radiation timestep must be set to a positive integer', &
              FATAL)
      endif
      if (.not. use_single_lw_sw_ts) then
        if (sw_rad_time_step <= 0) then
          call error_mesg ('radiation_driver_mod', &
           ' sw radiation timestep must be set to a positive integer', &
              FATAL)
        endif
      endif

      if (use_single_lw_sw_ts .and. (sw_rad_time_step /= 0.0 .and. &
          sw_rad_time_step /= rad_time_step) ) then
        call error_mesg ('radiation_driver', &
         'to avoid confusion, sw_rad_time_step must either remain at &
                  &default value of 0.0, or be same as rad_time_step &
                        &when use_single_lw_sw_ts is .true.', FATAL)
      endif
      if (use_single_lw_sw_ts) then
        sw_rad_time_step = rad_time_step
      endif
      lw_rad_time_step = rad_time_step
      Rad_control%rad_time_step = rad_time_step
      Rad_control%rad_time_step_iz  = .true.         
      Rad_control%lw_rad_time_step = lw_rad_time_step
      Rad_control%lw_rad_time_step_iz  = .true.         
      Rad_control%sw_rad_time_step = sw_rad_time_step
      Rad_control%sw_rad_time_step_iz  = .true.         

      if (MOD(INT(SECONDS_PER_DAY), lw_rad_time_step) /= 0) then
        call error_mesg ('radiation_driver_mod', &
             'lw radiation timestep currently restricted to be an &
                       &integral factor of seconds in a day', FATAL)
      endif
      if (MOD(INT(SECONDS_PER_DAY), sw_rad_time_step) /= 0) then
        call error_mesg ('radiation_driver_mod', &
             'sw radiation timestep currently restricted to be an &
                       &integral factor of seconds in a day', FATAL)
      endif

!----------------------------------------------------------------------
!    store the radiation time step in a derived-type variable for 
!    transfer to other modules.
!----------------------------------------------------------------------
      Rad_control%rad_time_step = rad_time_step
      Rad_control%rad_time_step_iz = .true.

!----------------------------------------------------------------------
!    store the controls for hires cloudy coszen calculations.
!----------------------------------------------------------------------
      if (use_hires_coszen) then
        Rad_control%hires_coszen = .true.
      else
        Rad_control%hires_coszen = .false.
      endif
      Rad_control%hires_coszen_iz = .true.
      
      if (nzens_per_sw_rad_timestep > 1 .and. &
          .not. (use_hires_coszen) ) then
        call error_mesg ('radiation_driver_init', &
           'uncertainty in what is desired wrt nzens; if &
           &nzens_per_sw_rad_timestep is not default, &
           &use_hires_coszen must be set to .true.' , FATAL)
      endif
      if (use_hires_coszen)  then
        Rad_control%nzens = nzens_per_sw_rad_timestep
      else
        Rad_control%nzens = 1
      endif
      Rad_control%nzens_iz = .true.

!---------------------------------------------------------------------
!    define the dimensions of the local processors portion of the grid.
!---------------------------------------------------------------------
      id    = size(lonb,1) - 1 
      jd    = size(latb,2) - 1
      kmax  = size(pref,1) - 1 

!---------------------------------------------------------------------
!    save the number of special levels at which fluxes may be defined
!    for diagnostic purposes. 
!---------------------------------------------------------------------
      Rad_control%mx_spec_levs = MX_SPEC_LEVS
      Rad_control%mx_spec_levs_iz = .true.        

!---------------------------------------------------------------------
!    check for consistency if drop_upper_levels is activated.
!----------------------------------------------------------------------
      if (drop_upper_levels .and. all_level_radiation) then
          call error_mesg ( 'radiation_driver_mod',  &
            ' drop_upper_levels and all_level_radiation are '//&
                                         'incompatible', FATAL)
      endif

!---------------------------------------------------------------------
!    define the starting and ending vertical indices of the radiation
!    grid. if all_level_radiation is .true., then radiation is done
!    at all model levels. ks, ke are model-based coordinates, while
!    ksrad and kerad are radiation-grid based coordinates (ksrad always
!    is equal to 1). 
!---------------------------------------------------------------------
      if (all_level_radiation) then
        ks = 1
        ke = kmax
        kerad = kmax
        topmost_radiation_level = 1
      else
        if (topmost_radiation_level <= 0) then
          call error_mesg ('radiation_driver_mod', &
          ' when all_level_radiation is .false., topmost_radiation'//&
              '_level must be specified as a positive integer.', FATAL)
        endif
        if (drop_upper_levels) then
          ks = topmost_radiation_level
          ke = kmax
          kerad = ke - ks + 1
          call error_mesg ( ' radiation_driver_mod', &
            ' code has not been validated for all_level_radiation = '//&
               'false. DO NOT USE!', FATAL)
        else
          call error_mesg ( ' radiation_driver_mod', &
         ' currently only drop_upper_levels is available as option '//&
                           'when all_level_radiation = false.', FATAL)
        endif
      endif
       
!---------------------------------------------------------------------
!    exit if all_column_radiation is not .true. -- this option is not
!    yet certified.
!---------------------------------------------------------------------
      if (.not. all_column_radiation) then
        call error_mesg ('radiation_driver_mod',  &
          ' code currently not validated for all_column_radiation = '//&
                                  'false. DO NOT USE!', FATAL)
      endif

!----------------------------------------------------------------------
!    be sure both reference pressure profiles have been provided.
!----------------------------------------------------------------------
      if (size(pref,2) /= 2)    &
        call error_mesg ('radiation_driver_mod', &
         'must provide two reference pressure profiles (pref).', FATAL)

!---------------------------------------------------------------------
!    allocate space for variables which must be saved when sw fluxes
!    are renormalized or diagnostics are desired to be output on every
!    physics step.
!---------------------------------------------------------------------
        nzens = Rad_control%nzens
      if (renormalize_sw_fluxes .or. all_step_diagnostics) then
        allocate (solar_save             (id,jd))
        allocate (dum_idjd               (id,jd))
        allocate (flux_sw_surf_save      (id,jd,nzens))
        allocate (flux_sw_surf_dir_save      (id,jd,nzens))
        allocate (flux_sw_surf_refl_dir_save      (id,jd,nzens))
        allocate (flux_sw_surf_dif_save      (id,jd,nzens))
        allocate (flux_sw_down_vis_dir_save      (id,jd,nzens))
        allocate (flux_sw_down_vis_dif_save      (id,jd,nzens))
        allocate (flux_sw_down_total_dir_save      (id,jd,nzens))
        allocate (flux_sw_down_total_dif_save      (id,jd,nzens))
        allocate (flux_sw_vis_save      (id,jd,nzens))
        allocate (flux_sw_vis_dir_save      (id,jd,nzens))
        allocate (flux_sw_refl_vis_dir_save      (id,jd,nzens))
        allocate (flux_sw_vis_dif_save      (id,jd,nzens))
        allocate (sw_heating_save        (id,jd,kmax,nzens))
        allocate (tot_heating_save       (id,jd,kmax,nzens))
        allocate (dfsw_save              (id,jd,kmax+1,nzens))
        allocate (ufsw_save              (id,jd,kmax+1,nzens))
        allocate ( fsw_save              (id,jd,kmax+1,nzens))
        allocate ( hsw_save              (id,jd,kmax,nzens))
        allocate (swdn_special_save      (id,jd,MX_SPEC_LEVS,nzens))
        allocate (swup_special_save      (id,jd,MX_SPEC_LEVS,nzens))
        if (do_swaerosol_forcing) then
          allocate (dfsw_ad_save              (id,jd,kmax+1,nzens))
          allocate (ufsw_ad_save              (id,jd,kmax+1,nzens))
        endif
        if (do_clear_sky_pass) then 
          allocate (sw_heating_clr_save  (id,jd,kmax,nzens))
          allocate (tot_heating_clr_save (id,jd,kmax,nzens))
          allocate (dfswcf_save          (id,jd,kmax+1,nzens))
          allocate (ufswcf_save          (id,jd,kmax+1,nzens))
          allocate ( fswcf_save          (id,jd,kmax+1,nzens))
          allocate ( hswcf_save          (id,jd,kmax,nzens))
          allocate (flux_sw_down_total_dir_clr_save  (id,jd,nzens))
          allocate (flux_sw_down_total_dif_clr_save  (id,jd,nzens))
          allocate (flux_sw_down_vis_clr_save (id,jd,nzens)) 
          allocate (swdn_special_clr_save(id,jd, MX_SPEC_LEVS,nzens))
          allocate (swup_special_clr_save(id,jd, MX_SPEC_LEVS,nzens))
          if (do_swaerosol_forcing) then
            allocate (dfswcf_ad_save          (id,jd,kmax+1,nzens))
            allocate (ufswcf_ad_save          (id,jd,kmax+1,nzens))
          endif
        endif
      endif

!---------------------------------------------------------------------
!    allocate space for variables which must be saved when lw fluxes
!    are to be output on every physics step.
!---------------------------------------------------------------------
      if (all_step_diagnostics) then
        allocate (olr_save             (id,jd))
        allocate (lwups_save           (id,jd))
        allocate (lwdns_save           (id,jd))
        allocate (tdtlw_save           (id,jd,kmax))
        allocate (flxnet_save           (id,jd,kmax+1))
        allocate (netlw_special_save   (id,jd,MX_SPEC_LEVS))
        if (do_lwaerosol_forcing) then
          allocate (olr_ad_save             (id,jd))
          allocate (lwups_ad_save           (id,jd))
          allocate (lwdns_ad_save           (id,jd))
        endif
        if (do_clear_sky_pass) then
          allocate (olr_clr_save           (id,jd))
          allocate (lwups_clr_save         (id,jd))
          allocate (lwdns_clr_save         (id,jd))
          allocate (tdtlw_clr_save         (id,jd,kmax))
          allocate (flxnetcf_save         (id,jd,kmax+1))
          allocate (netlw_special_clr_save (id,jd,MX_SPEC_LEVS))
          if (do_lwaerosol_forcing) then
            allocate (olr_ad_clr_save           (id,jd))
            allocate (lwups_ad_clr_save         (id,jd))
            allocate (lwdns_ad_clr_save         (id,jd))
          endif
        endif
      endif

!---------------------------------------------------------------------
!    allocate space for the global integrals being accumulated in 
!    this module.
!---------------------------------------------------------------------
      allocate (olr_intgl(id,jd))
      allocate (swabs_intgl(id,jd))

!---------------------------------------------------------------------
!    allocate space for module variables to contain values which must
!    be saved between timesteps (these are used on every timestep,
!    but only calculated on radiation steps).
!---------------------------------------------------------------------
        allocate (Rad_output%tdt_rad     (id,jd,kmax,nzens))
        allocate (Rad_output%tdt_rad_clr (id,jd,kmax,nzens))
        allocate (Rad_output%tdtsw       (id,jd,kmax,nzens))
        allocate (Rad_output%tdtsw_clr   (id,jd,kmax,nzens))
        allocate (Rad_output%ufsw        (id,jd,kmax+1,nzens))
        allocate (Rad_output%dfsw        (id,jd,kmax+1,nzens))
        allocate (Rad_output%ufsw_clr    (id,jd,kmax+1,nzens))
        allocate (Rad_output%dfsw_clr    (id,jd,kmax+1,nzens))
        allocate (Rad_output%flxnet      (id,jd,kmax+1))
        allocate (Rad_output%flxnetcf    (id,jd,kmax+1))
        allocate (Rad_output%tdtlw       (id,jd,kmax))
        allocate (Rad_output%tdtlw_clr   (id,jd,kmax))
        allocate (Rad_output%flux_sw_surf_dir(id,jd,nzens))
        allocate (Rad_output%flux_sw_surf_refl_dir(id,jd,nzens))
        allocate (Rad_output%flux_sw_surf_dif(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_vis_dir(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_vis_dif(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_total_dir(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_total_dif(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_total_dir_clr(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_total_dif_clr(id,jd,nzens))
        allocate (Rad_output%flux_sw_down_vis_clr(id,jd,nzens))
        allocate (Rad_output%flux_sw_vis(id,jd,nzens))
        allocate (Rad_output%flux_sw_vis_dir(id,jd,nzens))
        allocate (Rad_output%flux_sw_refl_vis_dir(id,jd,nzens))
        allocate (Rad_output%flux_sw_vis_dif(id,jd,nzens))
        allocate (Rad_output%flux_sw_surf(id,jd,nzens))
        allocate (Rad_output%flux_lw_surf(id,jd))
        allocate (Rad_output%coszen_angle(id,jd))
        Rad_output%tdtsw     = 0.0
        Rad_output%tdtsw_clr = 0.0
        Rad_output%ufsw      = 0.0
        Rad_output%dfsw      = 0.0
        Rad_output%ufsw_clr  = 0.0
        Rad_output%dfsw_clr  = 0.0
        Rad_output%flxnet    = 0.0
        Rad_output%flxnetcf  = 0.0
        Rad_output%tdtlw     = 0.0
        Rad_output%tdtlw_clr = 0.0

!-----------------------------------------------------------------------
!    if two radiation restart files exist, exit.
!-----------------------------------------------------------------------
        if ( file_exist('INPUT/sea_esf_rad.res')  .and.     &
             file_exist('INPUT/radiation_driver.res') ) then 
          call error_mesg ('radiation_driver_mod',  &
         ' both sea_esf_rad.res and radiation_driver.res files are'//&
               ' present in INPUT directory. which one to use ?', FATAL)
        endif

   if  (using_restart_file) then

!----------------------------------------------------------------------
!    Register fields to be written out to restart file.
     call rad_driver_register_restart('radiation_driver.res.nc')

!-----------------------------------------------------------------------
!    if a valid restart file exists, call read_restart_file to read it.
!-----------------------------------------------------------------------
        if ( file_exist('INPUT/radiation_driver.res.nc')) then
          call restore_state(Rad_restart)
          if(in_different_file) call restore_state(Til_restart)
        else if ( (do_sea_esf_rad .and.   &
             (file_exist('INPUT/sea_esf_rad.res')  .or. &
              file_exist('INPUT/radiation_driver.res') )  ) .or. &
             (.not. do_sea_esf_rad .and.   &
               file_exist('INPUT/radiation_driver.res') )  ) then
          call error_mesg ('radiation_driver_mod', &
              'Native restarts no longer supported', FATAL)
!----------------------------------------------------------------------
!    if no restart file is present, initialize the needed fields until
!    the radiation package may be called. initial surface flux is set 
!    to 100 wm-2, and is only used for initial guess of sea ice temp.
!    set rad_alarm to be 1 second from now, ie., on the first step of 
!    the job.
!-----------------------------------------------------------------------
        else
          lwrad_alarm                = 1
          swrad_alarm                = 1
          if (mpp_pe() == mpp_root_pe() ) then
          call error_mesg ('radiation_driver_mod', &
           'radiation to be calculated on first step: no restart file&
                                                 & present', NOTE)
          endif
          Rad_output%tdt_rad       = 0.0
          Rad_output%tdt_rad_clr   = 0.0
          Rad_output%tdtlw         = 0.0
          Rad_output%flux_sw_surf  = surf_flx_init
!!! BETTER INITIAL VALUES FOR THESE ARRAYS NEEDED ??
          Rad_output%flux_sw_surf_dir  = surf_flx_init
          Rad_output%flux_sw_surf_refl_dir  = surf_flx_init
          Rad_output%flux_sw_surf_dif  = surf_flx_init
!!! BETTER INITIAL VALUES FOR THESE ARRAYS NEEDED ??
          Rad_output%flux_sw_down_vis_dir  = 0.0
          Rad_output%flux_sw_down_vis_dif  = 0.0
          Rad_output%flux_sw_down_total_dir  = 0.0
          Rad_output%flux_sw_down_total_dif  = 0.0
          Rad_output%flux_sw_down_total_dir_clr  = 0.0
          Rad_output%flux_sw_down_total_dif_clr  = 0.0
          Rad_output%flux_sw_down_vis_clr  = 0.0
          Rad_output%flux_sw_vis  = 0.0
          Rad_output%flux_sw_vis_dir  = 0.0
          Rad_output%flux_sw_refl_vis_dir  = 0.0
          Rad_output%flux_sw_vis_dif  = 0.0
          Rad_output%flux_lw_surf  = surf_flx_init
          Rad_output%coszen_angle  = coszen_angle_init
          if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg ('radiation_driver_mod', &
           'no acceptable radiation restart file present; therefore'//&
           ' will initialize input fields', NOTE)
          endif
        endif

!---------------------------------------------------------------------
!    if not using restart file, then initialize fields it would contain.
!    it is the responsibility of the user to assure restart is on a
!    radiation timestep so that restart seamlessness is maintained. if
!    restart is done on a non-radiation step, restart seamlessness will 
!    be lost if a restart file is not available.
!---------------------------------------------------------------------
   else  ! (using_restart_file)
     lwrad_alarm                = 1
     swrad_alarm                = 1
     if (mpp_pe() == mpp_root_pe() ) then
       call error_mesg ('radiation_driver_mod', &
          'radiation to be calculated on first step: user asserts that&
           & this is a scheduled radiation step;  if it is not, &
                           &restart seamlessness will be lost ', NOTE)
     endif
     Rad_output%tdt_rad       = 0.0
     Rad_output%tdt_rad_clr   = 0.0
     Rad_output%tdtlw         = 0.0
     Rad_output%flux_sw_surf  = surf_flx_init
     Rad_output%flux_sw_surf_dir  = surf_flx_init
     Rad_output%flux_sw_surf_refl_dir  = surf_flx_init
     Rad_output%flux_sw_surf_dif  = surf_flx_init
     Rad_output%flux_sw_down_vis_dir  = 0.0
     Rad_output%flux_sw_down_vis_dif  = 0.0
     Rad_output%flux_sw_down_total_dir  = 0.0
     Rad_output%flux_sw_down_total_dif  = 0.0
     Rad_output%flux_sw_down_total_dir_clr  = 0.0
     Rad_output%flux_sw_down_total_dif_clr  = 0.0
     Rad_output%flux_sw_down_vis_clr  = 0.0
     Rad_output%flux_sw_vis  = 0.0
     Rad_output%flux_sw_vis_dir  = 0.0
     Rad_output%flux_sw_refl_vis_dir  = 0.0
     Rad_output%flux_sw_vis_dif  = 0.0
     Rad_output%flux_lw_surf  = surf_flx_init
     Rad_output%coszen_angle  = coszen_angle_init
   endif ! (using_restart_file)

!--------------------------------------------------------------------
!    do the initialization specific to the sea_esf_rad radiation
!    package.
!--------------------------------------------------------------------
      if (do_sea_esf_rad) then 

!---------------------------------------------------------------------
!    define control variables indicating whether the clear-sky forcing
!    should be calculated. set a flag to indicate that the variable
!    has been defined.
!---------------------------------------------------------------------
        Rad_control%do_totcld_forcing = do_clear_sky_pass
        Rad_control%do_totcld_forcing_iz = .true.

!---------------------------------------------------------------------
!    define control variables indicating whether the aerosol forcings
!    should be calculated. set a flag to indicate that the variables
!    have been defined.
!---------------------------------------------------------------------
        Rad_control%do_lwaerosol_forcing = do_lwaerosol_forcing
        Rad_control%do_lwaerosol_forcing_iz = .true.
        Rad_control%do_swaerosol_forcing = do_swaerosol_forcing
        Rad_control%do_swaerosol_forcing_iz = .true.
        if (do_lwaerosol_forcing) then
          size_of_lwoutput = size_of_lwoutput + 1
          indx_lwaf = size_of_lwoutput 
          Rad_control%indx_lwaf = indx_lwaf
        endif
        if (do_swaerosol_forcing) then
          size_of_swoutput = size_of_swoutput + 1
          indx_swaf = size_of_swoutput 
          Rad_control%indx_swaf = indx_swaf
        endif
        Rad_control%indx_lwaf_iz = .true.     
        Rad_control%indx_swaf_iz = .true.

!---------------------------------------------------------------------
!    initialize the modules that are accessed from radiation_driver_mod.
!---------------------------------------------------------------------
        call sea_esf_rad_init        (lonb, latb, pref(ks:ke+1,:))
        call cloudrad_package_init   (pref(ks:ke+1,:), lonb, latb,  &
                                      axes, Time,   &
                                      donner_meso_is_largescale)
        call aerosolrad_package_init (kmax, aerosol_names, lonb, latb)
        call rad_output_file_init    (axes, Time, aerosol_names, &
                                      aerosol_family_names)

!---------------------------------------------------------------------
!    do the initialization specific to the original fms radiation
!    package. 
!---------------------------------------------------------------------
      else
        call original_fms_rad_init (lonb, latb, pref, axes, Time, kmax)
      endif

!--------------------------------------------------------------------
!    initialize the astronomy_package.
!--------------------------------------------------------------------
      if (Sw_control%do_annual) then
        call astronomy_init (latb, lonb)
      else
        call astronomy_init
      endif

!---------------------------------------------------------------------
!    initialize the total number of columns in the processor's domain.
!---------------------------------------------------------------------
        total_pts = id*jd 

!-----------------------------------------------------------------------
!    check if optional radiative date should be used.
!-----------------------------------------------------------------------
        if (rad_date(1) > 1900 .and.                        &
            rad_date(2) >   0  .and. rad_date(2) < 13 .and. &
            rad_date(3) >   0  .and. rad_date(3) < 32 ) then
          use_rad_date = .true.
        else
          use_rad_date = .false.
        endif

!----------------------------------------------------------------------
!    define characteristics of desired diagnostic integrals. 
!----------------------------------------------------------------------
        call initialize_diagnostic_integrals

!----------------------------------------------------------------------
!    register the desired netcdf output variables with the 
!    diagnostics_manager.
!----------------------------------------------------------------------
        call diag_field_init (Time, axes)

!--------------------------------------------------------------------
!    initialize clocks to time portions of the code called from 
!    radiation_driver.
!--------------------------------------------------------------------
      misc_clock =    &
            mpp_clock_id ('   Physics_down: Radiation: misc', &
                grain = CLOCK_MODULE)
      clouds_clock =   &
            mpp_clock_id ('   Physics_down: Radiation: clds', &
               grain = CLOCK_MODULE)
      calc_clock =    &
            mpp_clock_id ('   Physics_down: Radiation: calc', &
                grain = CLOCK_MODULE)

!---------------------------------------------------------------------
!     if Rad_time is unchanging between timesteps, or the same day is being
!     repeated, switch to the alternative seed generation procedure to
!     assure unique temporal and spatial seeds for the stochastic cloud
!     parameterization.
!---------------------------------------------------------------------
      if ( (rsd .or. use_rad_date)) then  
        Cldrad_control%use_temp_for_seed = .true.
        call error_mesg ('cloud_spec_init', &
             'Will use temp as basis for stochastic cloud seed; &
                                &Rad_time is not monotonic', NOTE)
      endif
      Cldrad_control%use_temp_for_seed_iz = .true.

!---------------------------------------------------------------------
!    call check_derived_types to verify that all logical elements of
!    public derived-type variables stored in rad_utilities_mod but
!    initialized elsewhere have been initialized.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        call check_derived_types
      endif

!---------------------------------------------------------------------
!    verify that stochastic clouds have been activated if the COSP 
!    simulator output has been requested.
!---------------------------------------------------------------------
      if (do_cosp .and.   &
          (.not. Cldrad_control%do_stochastic_clouds) ) then
        call error_mesg ('radiation_driver_init', &
         'cannot call COSP simulator unless stochastic clouds are &
           &activated (do_stochastic_clouds in strat_clouds_W_nml)', &
                                                                  FATAL)
      endif

!--------------------------------------------------------------------
!    return the potential number of stochastic columns.
!--------------------------------------------------------------------
      ncol = Solar_spect%nbands + Cldrad_control%nlwcldb
 
!---------------------------------------------------------------------
!    set flag to indicate that module has been successfully initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine radiation_driver_init


!######################################################################
 
subroutine radiation_driver_time_vary (Time, Rad_gases_tv)
 
!---------------------------------------------------------------------
!    radiation_driver_time_vary calculates time-dependent, 
!    space-independent quantities needed within the modules of the 
!    radiation package.
!---------------------------------------------------------------------
 
type(time_type),               intent(in)     :: Time
type(radiative_gases_type),    intent(inout)  ::  Rad_gases_tv
 

      call aerosolrad_package_time_vary (Time)
      call sea_esf_rad_time_vary (Time, Rad_gases_tv)
 
end subroutine radiation_driver_time_vary
 
 
!####################################################################

subroutine radiation_driver_endts (is, js, Rad_gases_tv)

integer,                       intent(in)  :: is,js
type(radiative_gases_type),    intent(in)  ::  Rad_gases_tv

!---------------------------------------------------------------------

      call sum_diag_integral_field ('olr',    olr_intgl)
      call sum_diag_integral_field ('abs_sw', swabs_intgl )

      call aerosolrad_package_endts

      call sea_esf_rad_endts (Rad_gases_tv)

!---------------------------------------------------------------------
!    complete radiation step. if this was a radiation step, set the 
!    radiation alarm to go off rad_time_step seconds from now, and
!    set do_rad to false, so that radiation will not be calculated 
!    again until the alarm goes off.
!--------------------------------------------------------------------
      if (.not. always_calculate) then
        if (do_lw_rad) then
          lwrad_alarm = lwrad_alarm + lw_rad_time_step
          do_lw_rad = .false.
        endif
        if (do_sw_rad) then
          swrad_alarm = swrad_alarm + sw_rad_time_step
          do_sw_rad = .false.
        endif

        if (.not. do_lw_rad .and. .not. do_sw_rad)  then
          do_rad = .false.
        else
          do_rad = .true.
        endif

      endif  ! (always_calculate)

      Rad_control%do_lw_rad = do_lw_rad
      Rad_control%do_sw_rad = do_sw_rad

 
end subroutine radiation_driver_endts


!#####################################################################
! <SUBROUTINE NAME="radiation_driver">
!  <OVERVIEW>
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver (is, ie, js, je, Time, Time_next,  &
!                             lat, lon, Surface, Atmos_input, &
!                             Aerosol, Cld_spec, Rad_gases, &
!                             Lsc_microphys, Meso_microphys,    &
!                             Cell_microphys, Radiation, mask, kbot)
!
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <IN NAME="Time" TYPE="time_type">
!   current model time 
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   The time used for diagnostic output
!  </IN>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </INOUT>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </INOUT>
!  <INOUT NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data to radiation package
!  </INOUT>
!  <INOUT NAME="r" TYPE="real">
!   4 dimensional tracer array, last index is the number of all tracers
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud microphysical and physical parameters to radiation package, 
!                     contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!  </INOUT>
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties to radiation package, , contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale
!                      clouds
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                      clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale
!                      clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Radiation" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! <ERROR MSG="radiation_driver_init must first be called" STATUS="FALTA">
! You have not called radiation_driver_init before calling
!       radiation_driver.
! </ERROR>
! <ERROR MSG="Time_next <= Time" STATUS="FALTA">
! Time arguments to radiation_driver are producing a time step <= 0.
!       Check that the time argumnets passed to the physics_driver are
!       correct.
! </ERROR>
! </SUBROUTINE>
!
subroutine radiation_driver (is, ie, js, je, Time, Time_next,  &
                             lat, lon, Surface, Atmos_input, &
                             Aerosol, r, Cld_spec, Rad_gases, &
                             Lsc_microphys, Meso_microphys,    &
                             Cell_microphys, Shallow_microphys, &
                             Model_microphys, &
                             Radiation, Astronomy_inp, &
                             mask, kbot)

!---------------------------------------------------------------------
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!---------------------------------------------------------------------
 
!--------------------------------------------------------------------
integer,                      intent(in)           :: is, ie, js, je
type(time_type),              intent(in)           :: Time, Time_next
real, dimension(:,:),         intent(in)           :: lat, lon
type(surface_type),           intent(inout)        :: Surface
type(atmos_input_type),       intent(inout)        :: Atmos_input
type(aerosol_type),           intent(inout)        :: Aerosol  
real, dimension(:,:,:,:),     intent(inout)        :: r
type(cld_specification_type), intent(inout)        :: Cld_spec
type(radiative_gases_type),   intent(inout)        :: Rad_gases
type(microphysics_type),      intent(inout)        :: Lsc_microphys,&
                                                      Meso_microphys,&
                                                      Cell_microphys, &
                                                    Shallow_microphys, &
                                                      Model_microphys
type(rad_output_type),     intent(inout), optional :: Radiation
type(astronomy_inp_type),  intent(inout), optional :: Astronomy_inp
real, dimension(:,:,:),    intent(in),    optional :: mask
integer, dimension(:,:),   intent(in),    optional :: kbot
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time           current model time [ time_type (days, seconds) ] 
!      Time_next      time on next timestep, used as stamp for diagnos-
!                     tic output  [ time_type  (days, seconds) ]  
!      lat            latitude of model points  [ radians ]
!      lon            longitude of model points [ radians ]
!
!   intent(inout) variables:
!
!      Surface        surface_type structure, contains variables 
!                     defining the surface characteristics, including
!                     the following component referenced in this 
!                     routine:
!
!         asfc          surface albedo  [ dimensionless ]
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state, including the follow-
!                     ing component referenced in this routine
!
!         tsfc          surface temperature [ deg K ]
!
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds, passed through to lower
!                     level routines
!      Meso_microphys microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     meso-scale clouds, passed through to lower
!                     level routines
!      Cell_microphys microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     convective cell-scale clouds, passed through to 
!                     lower level routines
!
!   intent(inout), optional variables:
!
!      Radiation      rad_output_type structure, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.  present when running gcm,
!                     not present when running sa_gcm or standalone
!                     columns mode. variables defined here are:
!
!        tdt_rad         radiative (sw + lw) heating rate
!                        [ deg K / sec ]
!        flux_sw_surf    net (down-up) sw surface flux 
!                        [ watts / m^^2 ]
!        flux_lw_surf    downward lw surface flux 
!                        [ watts / m^^2 ]
!        coszen_angle    cosine of the zenith angle which will be used 
!                        for the next ocean_albedo calculation 
!                        [ dimensionless ]
!        tdtlw           longwave heating rate
!                        [ deg K / sec ]
!      Astronomy_inp  astronomy_input_type structure, optionally used
!                     to input astronomical forcings, when it is desired
!                     to specify them rather than use astronomy_mod.
!                     Used in various standalone applications.
!
!   intent(in), optional variables:
!
!        mask            present when running eta vertical coordinate,
!                        mask to remove points below ground
!        kbot            present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      type(cldrad_properties_type)       :: Cldrad_props
      type(astronomy_type)               :: Astro, Astro2
      type(lw_output_type), dimension(size_of_lwoutput) :: Lw_output
      type(sw_output_type), dimension(size_of_swoutput) :: Sw_output
      type(fsrad_output_type)            :: Fsrad_output
      type(aerosol_properties_type)      :: Aerosol_props
      type(aerosol_diagnostics_type)     :: Aerosol_diags

      real, dimension (ie-is+1, je-js+1) :: flux_ratio, &
                                            lat_uniform, lon_uniform
      integer :: nz
 
!-------------------------------------------------------------------
!   local variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [cldrad_properties_type]
!      Astro             astronomical properties on model grid, usually
!                        valid over radiation timestep
!                        [astronomy_type]
!      Astro2            astronomical properties on model grid, valid 
!                        over current physics timestep
!                        [astronomy_type]
!      Lw_output         sea longwave output fields on model grid,
!                        [lw_output_type]
!      Sw_output         esf shortwave output fields on model grid,
!                        [sw_output_type]
!      Fsrad_output      original fms radiation output fields on model
!                        grid, [fsrad_output_type]
!      flux_ratio        value  used to renormalize sw fluxes and 
!                        heating rates to account for earth-sun motion
!                        during the radiation timestep
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)
     
!---------------------------------------------------------------------
!    if this is a radiation step, or if the astronomical inputs to
!    radiation (solar, cosz, fracday, rrsun) need to be obtained 
!    because of time averaging or renormalization, call 
!    obtain_astronomy_variables to do so.
!---------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      if (do_rad .or. renormalize_sw_fluxes .or.   &
          present(Astronomy_inp)) then 
        if (use_uniform_solar_input) then
          if (present (Astronomy_inp)) then
            call error_mesg ('radiation_driver_mod', &
              'cannot specify both use_uniform_solar_input AND use&
              & Astronomy_inp to specify astronomical variables', &
                                                               FATAL)
          endif
          lat_uniform(:,:) = lat_for_solar_input
          lon_uniform(:,:) = lon_for_solar_input
          call obtain_astronomy_variables (is, ie, js, je,  &
                                           lat_uniform, lon_uniform,  &
                                           Astro, Astro2)
        else
          if (present (Astronomy_inp)) then
            Sw_control%do_diurnal = .false.
            Sw_control%do_annual = .false.
            Sw_control%do_daily_mean = .false.
          endif
          call obtain_astronomy_variables (is, ie, js, je, lat, lon,  &
                                           Astro, Astro2, &
                                           Astronomy_inp =  &
                                                          Astronomy_inp)
        endif
      endif

!     print *, 'before aerosol  ', mpp_pe()
      if (do_rad) then
        if (Rad_control%do_aerosol) then
          call aerosolrad_package_alloc (ie-is+1, je-js+1,  &
                              size(Aerosol%aerosol,3), Aerosol_props)
          call aerosol_radiative_properties (is, ie, js, je, &
                                             Rad_time,   &
                                             Atmos_input%pflux, &
                                             Aerosol_diags, &
                                             Aerosol, Aerosol_props)
!         allocate (Aerosol_diags%extopdep (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                           size(Aerosol%aerosol,3), &
!                                           size(Aerosol%aerosol,4) ))
!         Aerosol_diags%extopdep = 0.0
!         allocate (Aerosol_diags%absopdep (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                           size(Aerosol%aerosol,3), &
!                                           size(Aerosol%aerosol,4) ))
!         Aerosol_diags%absopdep = 0.0
!         allocate (Aerosol_diags%extopdep_vlcno    &
!                                          (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                           size(Aerosol%aerosol,3),3)) 
!         Aerosol_diags%extopdep_vlcno = 0.0
!         allocate (Aerosol_diags%absopdep_vlcno  &
!                                          (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                           size(Aerosol%aerosol,3),3))
!         Aerosol_diags%absopdep_vlcno = 0.0
!         allocate (Aerosol_diags%sw_heating_vlcno  &
!                                          (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                           size(Aerosol%aerosol,3)))
!         Aerosol_diags%sw_heating_vlcno = 0.0
!         allocate (Aerosol_diags%lw_extopdep_vlcno    &
!                                          (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                         size(Aerosol%aerosol,3)+1,3)) 
!         Aerosol_diags%lw_extopdep_vlcno = 0.0
!         allocate (Aerosol_diags%lw_absopdep_vlcno  &
!                                          (size(Aerosol%aerosol,1), &
!                                           size(Aerosol%aerosol,2), &
!                                         size(Aerosol%aerosol,3)+1,3))
!         Aerosol_diags%lw_absopdep_vlcno = 0.0
        
        endif
      endif
      call mpp_clock_end (misc_clock)

!--------------------------------------------------------------------
!    when using the sea-esf radiation, call cloud_radiative_properties
!    to obtain the cloud-radiative properties needed for the radiation 
!    calculation. (these properties are obtained within radiation_calc
!    when executing the original fms radiation code). if these fields 
!    are to be time-averaged, this call is made on all steps; otherwise
!    just on radiation steps.
!--------------------------------------------------------------------
!     print *, 'before cloud_rad', mpp_pe()
      call mpp_clock_begin (clouds_clock)
      if (do_rad) then
        if (do_sea_esf_rad) then
          if (present(kbot) ) then
            call cloud_radiative_properties (     &
                         is, ie, js, je, Rad_time, Time, Time_next, Astro,&
                         Atmos_input, Cld_spec, Lsc_microphys,  &
                         Meso_microphys, Cell_microphys,    &
                         Shallow_microphys, Cldrad_props,  &
                         Model_microphys, kbot=kbot, mask=mask)
          else    

            call cloud_radiative_properties (      &
                         is, ie, js, je, Rad_time, Time, Time_next, Astro,&
                         Atmos_input, Cld_spec, Lsc_microphys,   &
                         Meso_microphys, Cell_microphys,    &
                       Shallow_microphys, Cldrad_props, Model_microphys)
          endif
        endif
      endif
      call mpp_clock_end (clouds_clock)

!---------------------------------------------------------------------
!    on radiation timesteps, call radiation_calc to determine new radia-
!    tive fluxes and heating rates.
!---------------------------------------------------------------------
!     print *, 'before _calc    ', mpp_pe()
      call mpp_clock_begin (calc_clock)
      if (do_rad) then
        call radiation_calc (is, ie, js, je, Rad_time, Time_next, lat, &
                             lon, Atmos_input, Surface, Rad_gases,  &
                             Aerosol_props, Aerosol, r, Cldrad_props, &
                             Cld_spec, Astro, Rad_output, Lw_output, &
                             Sw_output, Fsrad_output, Aerosol_diags, &
                             mask=mask,   &
                             kbot=kbot)       
      endif

      call mpp_clock_end (calc_clock)
!-------------------------------------------------------------------
!    on all timesteps, call update_rad_fields to update the temperature 
!    tendency and define the fluxes needed by other component models.
!    if the shortwave fluxes are to be renormalized because of the 
!    change in zenith angle since the last radiation timestep, that also
!    is done in this subroutine. 
!-------------------------------------------------------------------
!     print *, 'before update   ', mpp_pe()
      call mpp_clock_begin (misc_clock)
!     if (Environment%running_gcm .or.  &
!         Environment%running_sa_model .or. &
!         (Environment%running_standalone .and. &
!          Environment%column_type == 'fms')) then
        call update_rad_fields (is, ie, js, je, Time_next, Astro2, &
                                Sw_output, Astro, Rad_output,    &
                                flux_ratio)

!-------------------------------------------------------------------
!    call produce_radiation_diagnostics to produce radiation 
!    diagnostics, both fields and integrals.
!-------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call produce_radiation_diagnostics        &
                            (is, ie, js, je, Time_next, Time, lat, &
                             Atmos_input%tsfc, Surface,  &
                             flux_ratio,  Astro, Rad_output,  &
                             Rad_gases, Lw_output=Lw_output,&
                             Sw_output=Sw_output,  &
                             Cld_spec=Cld_spec,  &
                             Lsc_microphys=Lsc_microphys)
        else
          call produce_radiation_diagnostics        &
                            (is, ie, js, je, Time_next, Time, lat, &
                             Atmos_input%tsfc, Surface,  &
                             flux_ratio,  Astro, Rad_output,  &
                             Rad_gases, Fsrad_output=Fsrad_output, &
                             mask=mask)
        endif 

!---------------------------------------------------------------------
!    call write_rad_output_file to produce a netcdf output file of 
!    radiation-package-relevant variables. note that this is called
!    only on radiation steps, so that the effects of sw renormalization
!    will not be seen in the variables of the data file written by
!    write_rad_output_file.
!---------------------------------------------------------------------
        if (do_lw_rad .and. do_sw_rad .and. do_sea_esf_rad) then
          if (Rad_control%do_aerosol) then
            call write_rad_output_file (is, ie, js, je,  &
                                        Atmos_input, Surface,   &
                                        Rad_output, Sw_output(1),  &
                                        Lw_output(1), Rad_gases,   & 
                                        Cldrad_props, Cld_spec, & 
                                        Time_next, Time, &
                                        Aerosol=Aerosol, &
                                        Aerosol_props=Aerosol_props, &
                                        Aerosol_diags=Aerosol_diags)
          else
            call write_rad_output_file (is, ie, js, je,  &
                                        Atmos_input,Surface, &
                                        Rad_output, Sw_output(1),   &
                                        Lw_output(1), Rad_gases,   &
                                        Cldrad_props, Cld_spec, &
                                        Time_next, Time)
          endif
        endif ! (do_rad and do_sea_esf_rad)
!     endif  ! (running_gcm)

!---------------------------------------------------------------------
!    call deallocate_arrays to deallocate the array space associated 
!    with stack-resident derived-type variables.
!---------------------------------------------------------------------
        call deallocate_arrays (Cldrad_props, Astro, Astro2,    &
                                Aerosol_props, &
                                Lw_output, Fsrad_output, Sw_output, &
                                Aerosol_diags)

!--------------------------------------------------------------------
!    define the elements of the rad_output_type variable which will
!    return the needed radiation package output to the calling routine.
!    Radiation is currently present when running within a gcm, but
!    not present for other applications.
!--------------------------------------------------------------------
      if (present (Radiation)) then 
        nz = current_sw_zenith_step
        Radiation%coszen_angle(:,:) =      &
                                  Rad_output%coszen_angle(is:ie,js:je)
        Radiation%tdt_rad(:,:,:,1) =   &
                                  Rad_output%tdt_rad(is:ie,js:je,:,nz)
        Radiation%flux_sw_surf(:,:,1) =    &
                                Rad_output%flux_sw_surf(is:ie,js:je,nz)

        if (treat_sfc_refl_dir_as_dif) then
        Radiation%flux_sw_surf_dir(:,:,1) =   &
                            Rad_output%flux_sw_surf_dir(is:ie,js:je,nz)
        Radiation%flux_sw_surf_dif(:,:,1) =   &
                            Rad_output%flux_sw_surf_dif(is:ie,js:je,nz)
        Radiation%flux_sw_vis_dir (:,:,1) =   &
                            Rad_output%flux_sw_vis_dir (is:ie,js:je,nz)
        Radiation%flux_sw_vis_dif (:,:,1) =   &
                            Rad_output%flux_sw_vis_dif (is:ie,js:je,nz)
        else
        Radiation%flux_sw_surf_dir(:,:,1) =   &
                            Rad_output%flux_sw_surf_dir(is:ie,js:je,nz) - &
                      Rad_output%flux_sw_surf_refl_dir(is:ie,js:je,nz)
        Radiation%flux_sw_surf_dif(:,:,1) =   &
                            Rad_output%flux_sw_surf_dif(is:ie,js:je,nz) + &
                      Rad_output%flux_sw_surf_refl_dir(is:ie,js:je,nz)
        Radiation%flux_sw_vis_dir (:,:,1) =   &
                           Rad_output%flux_sw_vis_dir (is:ie,js:je,nz) - &
                       Rad_output%flux_sw_refl_vis_dir (is:ie,js:je,nz)
        Radiation%flux_sw_vis_dif (:,:,1) =   &
                          Rad_output%flux_sw_vis_dif (is:ie,js:je,nz) + &
                       Rad_output%flux_sw_refl_vis_dir (is:ie,js:je,nz)
        endif

        Radiation%flux_sw_down_vis_dir(:,:,1) =   &
                         Rad_output%flux_sw_down_vis_dir(is:ie,js:je,nz)
        Radiation%flux_sw_down_vis_dif(:,:,1) =   &
                         Rad_output%flux_sw_down_vis_dif(is:ie,js:je,nz)
        Radiation%flux_sw_down_total_dir(:,:,1) =   &
                       Rad_output%flux_sw_down_total_dir(is:ie,js:je,nz)
        Radiation%flux_sw_down_total_dif(:,:,1) =   &
                      Rad_output%flux_sw_down_total_dif(is:ie,js:je,nz)
        Radiation%flux_sw_vis (:,:,1) =   &
                               Rad_output%flux_sw_vis (is:ie,js:je,nz)
        Radiation%flux_lw_surf(:,:)    =   &
                                  Rad_output%flux_lw_surf(is:ie,js:je)
        Radiation%flxnet(:,:,:)         =     &
                                  Rad_output%flxnet(is:ie,js:je,:)
        Radiation%tdtlw(:,:,:)         =     &
                                  Rad_output%tdtlw(is:ie,js:je,:)   
        Radiation%ufsw(:,:,:,1) = Rad_output%ufsw(is:ie,js:je,:,nz)
        Radiation%dfsw(:,:,:,1) = Rad_output%dfsw(is:ie,js:je,:,nz)
        Radiation%flxnetcf(:,:,:)         =     &
                                  Rad_output%flxnetcf(is:ie,js:je,:)
        Radiation%ufsw_clr(:,:,:,1) = Rad_output%ufsw_clr(is:ie,js:je,:,nz)
        Radiation%dfsw_clr(:,:,:,1) = Rad_output%dfsw_clr(is:ie,js:je,:,nz)
      endif
      call mpp_clock_end (misc_clock)

!---------------------------------------------------------------------




end subroutine radiation_driver



!#####################################################################
! <SUBROUTINE NAME="define_rad_times">
!  <OVERVIEW>
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  define_rad_times (Time, Time_next, Rad_time_out,    &
!                             need_aerosols, need_clouds, need_gases,  &
!                             need_basic)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current model time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   The time used for diagnostic output
!  </IN>
!  <INOUT NAME="Rad_time_out" TYPE="time_type">
!   time at which the climatologically-determined,
!                     time-varying input fields to radiation should 
!                     apply    
!  </INOUT>
!  <OUT NAME="need_aerosols" TYPE="logical">
!   aersosol input data is needed on this step ?
!  </OUT>
!  <OUT NAME="need_clouds" TYPE="logical">
!   cloud input data is needed on this step ?
!  </OUT>
!  <OUT NAME="need_gases" TYPE="logical">
!   radiative gas input data is needed on this step ?
!  </OUT>
!  <OUT NAME="need_basic" TYPE="logical">
!   atmospheric input fields are needed on this step ?
!  </OUT>
! </SUBROUTINE>
!
subroutine define_rad_times (Time, Time_next, Rad_time_out,    &
                             need_aerosols, need_clouds, need_gases,  &
                             need_basic)

!--------------------------------------------------------------------
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!-------------------------------------------------------------------- 

!---------------------------------------------------------------------
type(time_type), intent(in)     ::  Time, Time_next
type(time_type), intent(inout)  ::  Rad_time_out
logical,         intent(out)    ::  need_aerosols, need_clouds,   &
                                    need_gases, need_basic
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     Time            current model time  
!                     [ time_type, days and seconds]
!     Time_next       model time on the next atmospheric timestep
!                     [ time_type, days and seconds]
!     
!   intent(inout) variables:
!
!     Rad_time_out    time at which the climatologically-determined, 
!                     time-varying input fields to radiation should 
!                     apply    
!                     [ time_type, days and seconds]
!
!   intent(out) variables:
!
!     need_aerosols   aersosol input data is needed on this step ?
!     need_clouds     cloud input data is needed on this step ?
!     need_gases      radiative gas input data is needed on this step ?
!     need_basic      atmospheric input fields are needed on this step ?
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer        :: year, month, day, sec
      integer        :: dum, tod(3)
      integer        :: nband
      type(time_type) :: Solar_time

!---------------------------------------------------------------------
!   local variables:
!
!      day            day component of atmospheric timestep
!                     [ days ]
!      sec            seconds component of atmospheric timestep
!                     [ seconds ]
!      dum            dummy variable
!      tod            hours, minutes and seconds components of current
!                     time
!                     [ hours, minutes, seconds ]
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    store the atmospheric timestep into a module variable for later
!    use.
!--------------------------------------------------------------------
      call get_time (Time_next-Time, sec, day)    
      dt = day*SECONDS_PER_DAY + sec

!--------------------------------------------------------------------
!    verify that the radiation timestep is an even multiple of the 
!    physics timestep.
!---------------------------------------------------------------------
      if (MOD(lw_rad_time_step, dt) /= 0) then
        call error_mesg ('radiation_driver_mod',  &
    ' lw radiation timestep is not integral multiple of physics step', &
                                                           FATAL)
      endif
      if (MOD(sw_rad_time_step, dt) /= 0) then
        call error_mesg ('radiation_driver_mod',  &
       ' sw radiation timestep is not integral multiple of physics step', &
                                                           FATAL)
      endif

      if (MOD(sw_rad_time_step/nzens_per_sw_rad_timestep, dt) /= 0) then
        call error_mesg ( 'radiation_driver_mod', &
         'requested nzens per sw timestep incompatible with physics &
                                                    &timestep', FATAL)
      endif

!-------------------------------------------------------------------
!    for the standalone case, new radiation outputs are calculated on 
!    every step, using climatological variable values at the time spec-
!    ified by the input argument Time. 
!-------------------------------------------------------------------
      if (always_calculate) then
        do_rad = .true.
        do_sw_rad = .true.
        do_lw_rad = .true.
        Rad_time = Time
        current_sw_zenith_step = 1
        Rad_control%do_lw_rad = do_lw_rad
        Rad_control%do_sw_rad = do_sw_rad

!--------------------------------------------------------------------
!    if running a gcm aplication, if this is the first call by this
!    processor on this time step to radiation_driver (i.e. num_pts = 0),
!    determine if this is a radiation time step by decrementing the time
!    to alarm by the current model timestep.  if the alarm "goes off", 
!    i.e., is .le. 0, set do_rad to true, indicating this is a radiation
!    step. otherwise set it to .false. . 
!--------------------------------------------------------------------
      else
        if (num_pts == 0)  then
          lwrad_alarm = lwrad_alarm -  dt
          swrad_alarm = swrad_alarm -  dt
        endif
        if (lwrad_alarm <= 0) then
          do_lw_rad = .true.
        else
          do_lw_rad = .false.
        endif
        if (swrad_alarm <= 0) then
          do_sw_rad = .true.
          current_sw_zenith_step = 1
        else
          do_sw_rad = .false.
          if (use_hires_coszen) then
            current_sw_zenith_step = current_sw_zenith_step + 1
          endif
        endif
        if (do_sw_rad .or. do_lw_rad) then
           do_rad = .true.
        else
          do_rad = .false.
        endif
      Rad_control%do_lw_rad = do_lw_rad
      Rad_control%do_sw_rad = do_sw_rad

!-------------------------------------------------------------------
!    define the time to be used in defining the time-varying input 
!    fields for the radiation calculation (Rad_time). 
!-------------------------------------------------------------------
        if (rsd) then

!--------------------------------------------------------------------
!    if this is a repeat-same-day (rsd) experiment, define Rad_time
!    as the specified year-month-day (rad_date(1:3)), and the 
!    hr-min-sec of the current time (Time).
!---------------------------------------------------------------------
          if (.not. use_rad_date)   &
            call error_mesg ('radiation_driver_mod', &  
              'if (rsd), must set rad_date(1:3) to valid date', FATAL)
            call get_date (Time, dum, dum, dum, tod(1), tod(2), tod(3))
            Rad_time = set_date (rad_date(1), rad_date(2),& 
                                 rad_date(3), tod(1), tod(2), &
                                 tod(3))

!---------------------------------------------------------------------
!    if the specified date option is active, define Rad_time to be that
!    date and time.
!----------------------------------------------------------------------
        else if (use_rad_date) then
          Rad_time = set_date (rad_date(1), rad_date(2), rad_date(3),  &
                               rad_date(4), rad_date(5), rad_date(6))

!---------------------------------------------------------------------
!    if neither of these special cases is active, define Rad_time as
!    the current time (Time).
!---------------------------------------------------------------------
        else
          Rad_time = Time
        endif  ! (rsd)
      endif  ! (always_calculate)

!---------------------------------------------------------------------
!    define the solar_constant appropriate at Rad_time, including any
!    offset defined via the namelist.
!---------------------------------------------------------------------
      if (Rad_control%time_varying_solar_constant) then
        if (size(Solar_spect%solflxband(:)) /= numbands_lean) then
          call error_mesg ('radiation_driver_mod', &
             'bands present in solar constant time data differs from &
               &model parameterization band number', FATAL)
        endif

!--------------------------------------------------------------------
!    define time to be used for solar input data.
!--------------------------------------------------------------------
        if (negative_offset) then
          Solar_time = Rad_time - Solar_offset
        else
          Solar_time = Rad_time + Solar_offset
        endif
        call get_date (Solar_time, year, month, dum, dum, dum, dum)

!--------------------------------------------------------------------
!    define input value based on year and month of Solar_time.
!--------------------------------------------------------------------
        if (year < first_yr_lean) then
          Sw_control%solar_constant = solflxtot_lean_ann_1882
          do nband=1,numbands_lean
            Solar_spect%solflxband(nband) =  &
                       Solar_spect%solflxband_lean_ann_1882(nband)
          end do
        else if (year > last_yr_lean) then
          Sw_control%solar_constant = solflxtot_lean_ann_2000           
          do nband=1,numbands_lean
            Solar_spect%solflxband(nband) =  &
                  Solar_spect%solflxband_lean_ann_2000(nband)
          end do
        else
          Sw_control%solar_constant =   &
                            solflxtot_lean(year-first_yr_lean+1, month)
          do nband=1,numbands_lean
            Solar_spect%solflxband(nband) =  &
            Solar_spect%solflxband_lean(year-first_yr_lean+1, month, nband)
          end do
        endif
      endif
          
!--------------------------------------------------------------------
!    set a logical variable indicating whether radiative gas input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (do_rad) then
        need_gases = .true.
      else
        need_gases = .false.
      endif

!--------------------------------------------------------------------
!    set a logical variable indicating whether aerosol input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (do_rad  .and. Rad_control%do_aerosol) then
        need_aerosols = .true.
      else
        need_aerosols = .false.
      endif

!--------------------------------------------------------------------
!    set a logical variable indicating whether cloud input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (do_sea_esf_rad .and. do_rad) then
        need_clouds = .true.
      else
        need_clouds = .false.
      endif
      
!--------------------------------------------------------------------
!    set a logical variable indicating whether atmospheric input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (need_clouds .or. need_aerosols .or. need_gases) then
        need_basic = .true.
      else
        need_basic = .false.
      endif
   
!---------------------------------------------------------------------
!    place the time at which radiation is to be applied into an output
!    variable.
!---------------------------------------------------------------------
      Rad_time_out = Rad_time

!---------------------------------------------------------------------



end subroutine define_rad_times


!######################################################################
! <SUBROUTINE NAME="define_atmos_input_fields">
!  <OVERVIEW>
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
!                                      t, q, ts, r, gavg_rrv, Atmos_input, &
!                                      cloudtemp, cloudvapor, &
!                                      aerosoltemp, aerosolvapor, &
!                                      aerosolpress, kbot)  
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure at half levels
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at full levels
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity of water vapor at full levels
!  </IN>
!  <IN NAME="ts" TYPE="real">
!   surface temperature
!  </IN>
!  <IN NAME="r" TYPE="real">
!   tracer array
!  </IN>
!  <IN NAME="gavg_rrv" TYPE="real">
!   global average array of tracer volume mixxing ratio
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!  </INOUT>
!  <IN NAME="cloudtemp" TYPE="real">
!    temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!  </IN>
!  <IN NAME="cloudvapor" TYPE="real">
!   water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!  </IN>
!  <IN NAME="aerosoltemp" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by aerosol param-
!                         eterization 
!  </IN>
!  <IN NAME="aerosolvapor" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by aerosol param-
!                         eterization 
!  </IN>
!  <IN NAME="aerosolpress" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         pressure field to be used by aerosol param-
!                         eterization
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!   present when running eta vertical coordinate,
!                         index of lowest model level above ground
!  </IN>
! </SUBROUTINE>
!
subroutine define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
                                      t, q, ts, r, gavg_rrv, Atmos_input, &
                                      cloudtemp, cloudvapor, &
                                      aerosoltemp, aerosolvapor, &
                                      aerosolpress, kbot)     

!---------------------------------------------------------------------
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:,:),  intent(in)              :: pfull, phalf, t, q
real, dimension(:,:),    intent(in)              :: ts
real, dimension(:),      intent(in)              :: gavg_rrv
real, dimension(:,:,:,:),intent(in)              :: r
type(atmos_input_type),  intent(inout)           :: Atmos_input
integer, dimension(:,:), intent(in), optional    :: kbot
real, dimension(:,:,:),  intent(in), optional    :: cloudtemp,    &
                                                    cloudvapor, &
                                                    aerosoltemp, &
                                                    aerosolvapor, &
                                                    aerosolpress

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      pfull        pressure at full levels [ kg / (m s^2) ]
!      phalf        pressure at half levels [ kg / (m s^2) ]
!      t            temperature at full levels [ deg K]
!      q            specific humidity of water vapor at full levels
!                   [ dimensionless ]
!      ts           surface temperature  [ deg K ]
!
!   intent(out) variables:
!
!      Atmos_input   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!         psfc          surface pressure 
!                       [ (kg /( m s^2) ] 
!         tsfc          surface temperature
!                       [ deg K ]
!         temp          temperature at model levels (1:nlev), surface
!                       temperature is stored at value nlev+1; if eta
!                       coordinates, surface value stored in below 
!                       ground points
!                       [ deg K ]
!         press         pressure at model levels (1:nlev), surface 
!                       pressure is stored at index value nlev+1
!                       [ (kg /( m s^2) ] 
!         rh2o          mixing ratio of water vapor at model full levels
!                       [ non-dimensional ]
!         deltaz        model vertical grid separation
!                       [meters]
!         pflux         average of pressure at adjacent model levels
!                       [ (kg /( m s^2) ] 
!         tflux         average of temperature at adjacent model levels
!                       [ deg K ]
!         rel_hum       relative humidity
!                       [ dimensionless ]
!         cloudtemp     temperature to be seen by clouds (used in 
!                       sa_gcm feedback studies) 
!                       [ degrees K ]
!         cloudvapor    water vapor to be seen by clouds (used in 
!                       sa_gcm feedback studies) 
!                       [ nondimensional ]
!         clouddeltaz   deltaz to be used in defining cloud paths (used
!                       in sa_gcm feedback studies)
!                       [ meters ]
!         aerosoltemp   temperature to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ degrees K ]
!         aerosolvapor  water vapor to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ nondimensional ]
!         aerosolpress  pressure field to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ Pa ]
!         aerosolrelhum relative humidity seen by aerosol package,
!                       used in sa_gcm feedback studies
!                       [ dimensionless ]
!
!   intent(in), optional variables:
!
!      kbot               present when running eta vertical coordinate,
!                         index of lowest model level above ground (???)
!      cloudtemp          temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ degrees K ]
!      cloudvapor         water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ nondimensional ]
!      aerosoltemp        required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by aerosol param-
!                         eterization 
!      aerosolvapor       required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by aerosol param-
!                         eterization 
!      aerosolpress       required in sa_gcm mode, absent otherwise:
!                         pressure field to be used by aerosol param-
!                         eterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
 
      integer :: i, j, k, kb
      integer :: kmax
      logical :: override
      type(time_type)  :: Data_time
      real, dimension (size(q,1), size(q,2), size(q,3)) :: q2
      real, dimension (size(t,1), size(t,2), size(t,3)) :: t2, pfull2
      real, dimension (size(t,1), size(t,2), size(t,3)+1) ::  phalf2
      real, dimension (size(ts,1), size(ts,2)) ::  ts2
      integer                  :: ico2

!---------------------------------------------------------------------
!  local variables
!
!     i, j, k      do loop indices
!     kb           vertical index of lowest atmospheric level (when
!                  using eta coordinates)
!     kmax         number of model layers
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!----------------------------------------------------------------------
!    define the number of model layers.
!----------------------------------------------------------------------
      kmax = size(t,3)

!---------------------------------------------------------------------
!    if the temperature, cloud, or aerosol input data is to be over-
!    riden, define the time slice of data which is to be used. allocate
!    storage for the temperature data which will be needed for these
!    cases.
!---------------------------------------------------------------------
      if (doing_data_override) then
        Data_time = Rad_time                                  
        if (overriding_temps .or. overriding_aerosol .or. &
            overriding_clouds) then

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's temper-
!    ature data from the override file. if the process fails, write
!    an error message; if it succeeds move the data fro the current
!    window into array t2.
!---------------------------------------------------------------------
          call data_override ('ATM', 'tnew', t2(:,:,1:kmax), Data_time ,  &
                              override=override,                  &
                              is_in=is, ie_in=ie, js_in=js, je_in=je)
          if ( .not. override) then
            call error_mesg ('radiation_driver_mod', &
                      'temp => t not overridden successfully', FATAL)
          endif
        else
          t2 = t
        endif

!---------------------------------------------------------------------
!    if the temperature data is to be overriden, allocate storage for 
!    the surface temperature data which will be needed in this cases.
!---------------------------------------------------------------------
        if (overriding_temps) then

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's surface
!    temperature data from the override file. if the process fails,
!    write an error message; if it succeeds move the data from the 
!    current window into array ts2, and also into array ts2.
!---------------------------------------------------------------------
          call data_override ('ATM', 'ts', ts2, Data_time ,  &
                              override=override,             &
                              is_in=is, ie_in=ie, js_in=js, je_in=je)
          if ( .not. override) then
            call error_mesg ('radiation_driver_mod', &
              't_surf => ts not overridden successfully', FATAL)
          else
            t2(:,:,kmax+1) = ts2(:,:)
          endif
        else
          ts2 = ts
        endif

!---------------------------------------------------------------------
!    if the humidity, cloud, or aerosol input data is to be over-
!    riden, define the time slice of data which is to be used. allocate
!    storage for the humidity data which will be needed for these
!    cases.
!---------------------------------------------------------------------
        if (overriding_sphum .or. overriding_aerosol .or. &
            overriding_clouds) then

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's surface
!    humidity data from the override file. if the process fails,
!    write an error message; if it succeeds move the data from the 
!    current window into array q2.
!---------------------------------------------------------------------
          call data_override ('ATM', 'q', q2, Data_time ,  &
                              override=override,           &
                              is_in=is, ie_in=ie, js_in=js, je_in=je)
          if ( .not. override) then
            call error_mesg ('radiation_driver_mod', &
                 'sphum => q not overridden successfully', FATAL)
          endif
        else
          q2 = q
        endif

!---------------------------------------------------------------------
!    if the aerosol input data is to be overriden, allocate storage 
!    for the pressure data which will be needed in this case.
!---------------------------------------------------------------------
        if (overriding_aerosol) then

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's pressure
!    data from the override file. if the process fails, write an error
!    message; if it succeeds move the data from the current window into
!    array pfull2 and phalf2.
!---------------------------------------------------------------------
          call data_override ('ATM', 'pfull2', pfull2,  &
                              Data_time , override=override, &
                              is_in=is, ie_in=ie, js_in=js, je_in=je)
          if ( .not. override) then
            call error_mesg ('radiation_driver_mod', &
                 'pressm => pfull2 not overridden successfully', FATAL)
          endif
          call data_override ('ATM', 'phalf2', phalf2,  &
                              Data_time, override=override,&
                              is_in=is, ie_in=ie, js_in=js, je_in=je)
          if ( .not. override) then
            call error_mesg ('radiation_driver_mod', &
                 'phalfm => phalf2 not overridden successfully', FATAL)
          endif
        else
          pfull2 = pfull
          phalf2(:,:,kmax+1) = phalf(:,:,kmax+1)
        endif
        
!---------------------------------------------------------------------
!    if not doing data_override, define the arrays which will be 
!    used to define the components of Atmos_input%.
!---------------------------------------------------------------------
      else
        t2 = t
        ts2 = ts
        q2 = q
        pfull2 = pfull
        phalf2(:,:,kmax+1) = phalf(:,:,kmax+1)
      endif

!---------------------------------------------------------------------
!    allocate space for the components of the derived type variable
!    Atmos_input.
!---------------------------------------------------------------------
      allocate ( Atmos_input%press(size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%phalf(size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%temp (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%rh2o (size(t,1), size(t,2), size(t,3)  ) )
      allocate ( Atmos_input%rel_hum(size(t,1), size(t,2),    &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%cloudtemp(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%cloudvapor(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%clouddeltaz(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%aerosoltemp(size(t,1), size(t,2),   &
                                                    size(t,3)  ) )
      allocate ( Atmos_input%aerosolpress(size(t,1), size(t,2),    &
                                                     size(t,3)+1) )
      allocate ( Atmos_input%aerosolvapor(size(t,1), size(t,2),   &
                                                     size(t,3)  ) )
      allocate ( Atmos_input%aerosolrelhum(size(t,1), size(t,2),   &
                                                      size(t,3)  ) )
      allocate ( Atmos_input%deltaz(size(t,1), size(t,2), size(t,3) ) )
      allocate ( Atmos_input%pflux (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%tflux (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%psfc (size(t,1), size(t,2)             ) )
      allocate ( Atmos_input%tsfc (size(t,1), size(t,2)             ) )

      if (use_co2_tracer_field) then
         allocate ( Atmos_input%tracer_co2(size(t,1), size(t,2), size(t,3) ) )
      endif
!---------------------------------------------------------------------
!    define the cloudtemp component of Atmos_input. 
!---------------------------------------------------------------------
      if (present (cloudtemp) ) then
        Atmos_input%cloudtemp(:,:,:)   = cloudtemp(:,:,:)
      else
        if (overriding_clouds) then
          Atmos_input%cloudtemp(:,:,:)   = t2(:,:,:)
        else
          Atmos_input%cloudtemp(:,:,:)   = t(:,:,:)
        endif
      endif

!---------------------------------------------------------------------
!    define the cloudvapor component of Atmos_input.
!---------------------------------------------------------------------
      if (present (cloudvapor) ) then
        Atmos_input%cloudvapor(:,:,:)   = cloudvapor(:,:,:)
      else
        if (overriding_clouds) then
          Atmos_input%cloudvapor(:,:,:)   = q2(:,:,:)
        else
          Atmos_input%cloudvapor(:,:,:)   = q(:,:,:)
        endif
      endif

!---------------------------------------------------------------------
!    define the aerosoltemp component of Atmos_input.
!---------------------------------------------------------------------
      if (present (aerosoltemp) ) then
        Atmos_input%aerosoltemp(:,:,:)   = aerosoltemp(:,:,:)
      else
        if (overriding_aerosol) then
          Atmos_input%aerosoltemp(:,:,:)   = t2(:,:,:)
        else
          Atmos_input%aerosoltemp(:,:,:)   = t(:,:,:)
        endif
      endif
 
!---------------------------------------------------------------------
!    define the aerosolvapor component of Atmos_input.
!---------------------------------------------------------------------
      if (present (aerosolvapor) ) then
        Atmos_input%aerosolvapor(:,:,:)   = aerosolvapor(:,:,:)
      else
        if (overriding_aerosol) then
          Atmos_input%aerosolvapor(:,:,:)   = q2(:,:,:)
        else
          Atmos_input%aerosolvapor(:,:,:)   = q(:,:,:)
        endif
      endif

!---------------------------------------------------------------------
!    define values of surface pressure and temperature.
!--------------------------------------------------------------------
      if (present(kbot)) then
        do j=1,je-js+1
          do i=1,ie-is+1
            kb = kbot(i,j)
            Atmos_input%psfc(i,j) = phalf2(i,j,kb+1)
          end do
        end do
      else
        Atmos_input%psfc(:,:) = phalf2(:,:,kmax+1)
      endif

      Atmos_input%tsfc(:,:) = ts2(:,:)

!------------------------------------------------------------------
!    define the atmospheric pressure and temperature arrays.
!------------------------------------------------------------------
      do k=1,kmax 
        Atmos_input%press(:,:,k) = pfull2(:,:,k)
        Atmos_input%phalf(:,:,k) = phalf(:,:,k)
        Atmos_input%temp (:,:,k) = t2(:,:,k)
      end do
      Atmos_input%press(:,:,kmax+1) = phalf2(:,:,kmax+1)
      Atmos_input%phalf(:,:,kmax+1) = phalf2(:,:,kmax+1)
      Atmos_input%temp (:,:,kmax+1) = ts2 (:,:)

!---------------------------------------------------------------------
!    define the aerosolpress component of Atmos_input.
!---------------------------------------------------------------------
      if (present (aerosolpress) ) then
        do k=1,kmax
          Atmos_input%aerosolpress(:,:,k)   = aerosolpress(:,:,k)
        end do
      else
        if (overriding_aerosol) then
          do k=1,kmax
            Atmos_input%aerosolpress(:,:,k)   = pfull2(:,:,k)
          end do
          Atmos_input%aerosolpress(:,:,kmax+1)   = phalf2(:,:,kmax+1)
        else
          do k=1,kmax
            Atmos_input%aerosolpress(:,:,k)   = pfull(:,:,k)
          end do
          Atmos_input%aerosolpress(:,:,kmax+1)   = phalf(:,:,kmax+1)
        endif
      endif
 
!------------------------------------------------------------------
!    if in eta coordinates, fill in underground temperatures with 
!    surface value.
!------------------------------------------------------------------
      if (present(kbot)) then
        do j=1,je-js+1
          do i=1,ie-is+1
            kb = kbot(i,j)
            if (kb < kmax) then
              do k=kb+1,kmax
                Atmos_input%temp(i,j,k) = Atmos_input%temp(i,j,kmax+1)
              end do
            endif
          end do
        end do
      endif

!------------------------------------------------------------------
!    when running the gcm, convert the input water vapor specific 
!    humidity field to mixing ratio. it is assumed that water vapor 
!    mixing ratio is the input in the standalone case.
!------------------------------------------------------------------
        if (use_mixing_ratio) then
          Atmos_input%rh2o (:,:,:) = q2(:,:,:)
        else
          if (.not. overriding_sphum .and. &
              .not. overriding_clouds .and.  &
              .not. overriding_aerosol) then
            Atmos_input%rh2o (:,:,:) = q2(:,:,:)/(1.0 - q2(:,:,:))
          else ! for override, values are already mixing ratio
            Atmos_input%rh2o (:,:,:) = q2(:,:,:)
          endif
          if (.not. overriding_clouds) then
            Atmos_input%cloudvapor(:,:,:) =    &
                                 Atmos_input%cloudvapor(:,:,:)/  &
                            (1.0 - Atmos_input%cloudvapor(:,:,:))
          endif
          if (.not. overriding_aerosol) then
            Atmos_input%aerosolvapor(:,:,:) =    &
                                 Atmos_input%aerosolvapor(:,:,:)/  &
                            (1.0 - Atmos_input%aerosolvapor(:,:,:))
          endif
        endif
 
!------------------------------------------------------------------
!    be sure that the magnitude of the water vapor mixing ratio field 
!    to be input to the radiation code is no smaller than the value of 
!    rh2o_lower_limit, which is 2.0E-07 when running the sea_esf
!    radiation code and 3.0e-06 when running the original radiation
!    code. Likewise, the temperature that the radiation code sees is
!    constrained to lie between 100K and 370K. these are the limits of
!    the tables referenced within the radiation package.
!      exception:
!    if do_h2o is false, the lower limit of h2o is zero, and radiation
!    tables will not be called.
!-----------------------------------------------------------------------
      if (do_rad) then
        Atmos_input%rh2o(:,:,ks:ke) =    &
            MAX(Atmos_input%rh2o(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%cloudvapor(:,:,ks:ke) =    &
            MAX(Atmos_input%cloudvapor(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%aerosolvapor(:,:,ks:ke) =    &
                    MAX(Atmos_input%aerosolvapor(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%temp(:,:,ks:ke) =     &
                    MAX(Atmos_input%temp(:,:,ks:ke), temp_lower_limit)
        Atmos_input%temp(:,:,ks:ke) =     &
                    MIN(Atmos_input%temp(:,:,ks:ke), temp_upper_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) =     &
                    MAX(Atmos_input%cloudtemp(:,:,ks:ke), temp_lower_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) =     &
                    MIN(Atmos_input%cloudtemp(:,:,ks:ke), temp_upper_limit)
        Atmos_input%aerosoltemp(:,:,ks:ke) =     &
                    MAX(Atmos_input%aerosoltemp(:,:,ks:ke), temp_lower_limit)
     Atmos_input%aerosoltemp(:,:,ks:ke) =     &
                    MIN(Atmos_input%aerosoltemp(:,:,ks:ke), temp_upper_limit)
      endif

!--------------------------------------------------------------------
!    call calculate_aulixiary_variables to compute pressure and 
!    temperature arrays at flux levels and an array of model deltaz.
!--------------------------------------------------------------------
      if (do_rad) then
        call calculate_auxiliary_variables (Atmos_input)
      endif

!RSH
!RSH   define here the values for Atmos_input%tracer_co2.
!RSH
!fil   the error message should never be printed as that code should never
!      be executed, it's an extra guard against user error.
      if (use_co2_tracer_field ) then
         ico2 = get_tracer_index(MODEL_ATMOS, 'co2')
         if(ico2 /= NO_TRACER) then
            Atmos_input%tracer_co2(:,:,:) = r(:,:,:,ico2)
            Atmos_input%g_rrvco2 = gavg_rrv(ico2)
         else
            call error_mesg('radiation_driver', &
              'ico2 cannot be NO_TRACER when use_co2_tracer_field is .true.', FATAL)
         endif
      endif


!----------------------------------------------------------------------


end subroutine define_atmos_input_fields 




!#####################################################################
! <SUBROUTINE NAME="define_surface">
!  <OVERVIEW>
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_surface (is, ie, js, je, albedo, land, Surface)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="land" TYPE="real">
!   fraction of grid box which is land 
!  </IN>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   surface_type structure to be valued
!  </INOUT>
! </SUBROUTINE>
!
subroutine define_surface (is, ie, js, je, albedo, albedo_vis_dir,   &
                           albedo_nir_dir, albedo_vis_dif, &
                           albedo_nir_dif, land, Surface)

!---------------------------------------------------------------------
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface.  
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:),    intent(in)              :: albedo, land, &
                                                    albedo_vis_dir,    &
                                                    albedo_nir_dir, &
                                                    albedo_vis_dif,    &
                                                    albedo_nir_dif
type(surface_type),      intent(inout)           :: Surface     

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      albedo       surface albedo  [ dimensionless ]
!      albedo_vis_dir surface visible direct albedo  [ dimensionless ]
!      albedo_nir_dir surface nir direct albedo  [ dimensionless ]
!      albedo_vis_dif surface visible diffuse albedo  [ dimensionless ]
!      albedo_nir_dif surface nir diffuse albedo  [ dimensionless ]
!      land         fraction of grid box which is land [ dimensionless ]
!
!   intent(out) variables:
!
!      Surface       surface_type structure, contains the 
!                    following components defined in this subroutine
!         asfc          surface albedo
!                       [ non-dimensional ]
!         asfc_vis_dir  surface direct visible albedo
!                       [ non-dimensional ]
!         asfc_nir_dir  surface direct nir albedo
!                       [ non-dimensional ]
!         asfc_vis_dif  surface diffuse visible albedo
!                       [ non-dimensional ]
!         asfc_nir_dif  surface diffuse nir albedo
!                       [ non-dimensional ]
!         land          fraction of grid box covered by land
!                       [ non-dimensional ]
!
!---------------------------------------------------------------------

     logical :: override
     type(time_type)  :: Data_time
     real, dimension (size(albedo,1), size(albedo,2)) :: albedo_vis_dir2,   &
                                                         albedo_nir_dir2, &
                                                         albedo_vis_dif2,  &
                                                         albedo_nir_dif2

!-------------------------------------------------------------------
!    verify that the module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

      if (do_rad) then
        if (doing_data_override) then        
!---------------------------------------------------------------------
!    if the albedo data is to be overriden, define the time from which
!    the data is to be retrieved.
!---------------------------------------------------------------------
          if (overriding_albedo) then
            Data_time = Rad_time 


!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's surface
!    albedo data from the override file. if the process fails,
!    write an error message; if it succeeds move the data from the 
!    current window into array albedo2.
!---------------------------------------------------------------------
!           call data_override ('ATM', 'albedonew', albedo_proc,   &
!                             Data_time, override=override)
!           if ( .not. override) then
!             call error_mesg ('radiation_driver_mod', &
!             'cvisrfgd => albedo not overridden successfully', FATAL)
!           else
!             albedo2(:,:) =      albedo_proc(is:ie,js:je)
!           endif
            
            call data_override ('ATM', 'albedo_nir_dir_new',   &
                                albedo_nir_dir2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'nirdir => albedo not overridden successfully', FATAL)
            endif

            call data_override ('ATM', 'albedo_nir_dif_new',   &
                                albedo_nir_dif2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'nirdif => albedo not overridden successfully', FATAL)
            endif

            call data_override ('ATM', 'albedo_vis_dir_new',   &
                                albedo_vis_dir2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
               'visdir => albedo not overridden successfully', FATAL)
            endif

            call data_override ('ATM', 'albedo_vis_dif_new',   &
                                albedo_vis_dif2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
              'visdif => albedo not overridden successfully', FATAL)
           endif

!--------------------------------------------------------------------
!    if albedo data is not being overriden, define albedo2 to be the 
!    model value of albedo.
!--------------------------------------------------------------------
          else
!           albedo2 = albedo
            albedo_vis_dir2 = albedo_vis_dir
            albedo_nir_dir2 = albedo_nir_dir
            albedo_vis_dif2 = albedo_vis_dif
            albedo_nir_dif2 = albedo_nir_dif
          endif
        else ! (doing data_override)       
!         albedo2 = albedo
          albedo_vis_dir2 = albedo_vis_dir
          albedo_nir_dir2 = albedo_nir_dir
          albedo_vis_dif2 = albedo_vis_dif
          albedo_nir_dif2 = albedo_nir_dif
        endif
      else ! (do_rad)
!       albedo2 = albedo
        albedo_vis_dir2 = albedo_vis_dir
        albedo_nir_dir2 = albedo_nir_dir
        albedo_vis_dif2 = albedo_vis_dif
        albedo_nir_dif2 = albedo_nir_dif
      endif ! (do_rad)

!---------------------------------------------------------------------
!    allocate space for the components of the derived type variable
!    Surface.     
!---------------------------------------------------------------------
      allocate (Surface%asfc (size(albedo,1), size(albedo,2)) )
      allocate (Surface%asfc_vis_dir (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%asfc_nir_dir (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%asfc_vis_dif (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%asfc_nir_dif (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%land (size(albedo,1), size(albedo,2)) )

 
!------------------------------------------------------------------
!    define the fractional land area of each grid box and the surface
!    albedo from the input argument values.
!------------------------------------------------------------------
      Surface%land(:,:) = land(:,:)
      Surface%asfc(:,:) = albedo (:,:)

!pjp  Should the albedos below all be set to albedo2,
!pjp  or should they be included in the override data,
!pjp  or should it not be changed?

      Surface%asfc_vis_dir(:,:) = albedo_vis_dir2(:,:)
      Surface%asfc_nir_dir(:,:) = albedo_nir_dir2(:,:)
      Surface%asfc_vis_dif(:,:) = albedo_vis_dif2(:,:)
      Surface%asfc_nir_dif(:,:) = albedo_nir_dif2(:,:)
     

!----------------------------------------------------------------------


end subroutine define_surface    



!#####################################################################
! <SUBROUTINE NAME="surface_dealloc">
!  <OVERVIEW>
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!  </OVERVIEW>
!  <DESCRIPTION>
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call surface_dealloc (Surface)
!  </TEMPLATE>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   surface_type structure to be deallocated
!  </INOUT>
! </SUBROUTINE>
!
subroutine surface_dealloc (Surface)

!----------------------------------------------------------------------
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!----------------------------------------------------------------------

type(surface_type), intent(inout) :: Surface

!--------------------------------------------------------------------
!   intent(inout) variable:
!
!      Surface        surface_type structure, contains variables 
!                     defining the surface albedo and land fraction
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    deallocate components of surface_type structure.
!-------------------------------------------------------------------
      deallocate (Surface%asfc)
      deallocate (Surface%asfc_vis_dir )
      deallocate (Surface%asfc_nir_dir )
      deallocate (Surface%asfc_vis_dif )
      deallocate (Surface%asfc_nir_dif )
      deallocate (Surface%land)

!--------------------------------------------------------------------


end subroutine surface_dealloc 



!#####################################################################
! <SUBROUTINE NAME="atmos_input_dealloc">
!  <OVERVIEW>
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!  </OVERVIEW>
!  <DESCRIPTION>
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call atmos_input_dealloc (Atmos_input)
!  </TEMPLATE>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!      atmos_input_type structure, contains variables 
!                     defining the atmospheric pressure, temperature
!                     and moisture distribution.
!  </INOUT>
! </SUBROUTINE>
!
subroutine atmos_input_dealloc (Atmos_input)

!----------------------------------------------------------------------
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!----------------------------------------------------------------------

type(atmos_input_type), intent(inout) :: Atmos_input

!--------------------------------------------------------------------
!   intent(inout) variable:
!
!      Atmos_input    atmos_input_type structure, contains variables 
!                     defining the atmospheric pressure, temperature
!                     and moisture distribution.
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    deallocate components of atmos_input_type structure.
!---------------------------------------------------------------------
      deallocate (Atmos_input%press      )
      deallocate (Atmos_input%phalf      )
      deallocate (Atmos_input%temp       )
      deallocate (Atmos_input%rh2o       )
      deallocate (Atmos_input%rel_hum    )
      deallocate (Atmos_input%pflux      )
      deallocate (Atmos_input%tflux      )
      deallocate (Atmos_input%deltaz     )
      deallocate (Atmos_input%psfc       )
      deallocate (Atmos_input%tsfc       )
      deallocate (Atmos_input%cloudtemp  )
      deallocate (Atmos_input%cloudvapor )
      deallocate (Atmos_input%clouddeltaz)
      deallocate (Atmos_input%aerosoltemp)
      deallocate (Atmos_input%aerosolvapor )
      deallocate (Atmos_input%aerosolpress )
      deallocate (Atmos_input%aerosolrelhum )
      if(ASSOCIATED(Atmos_input%tracer_co2)) deallocate(Atmos_input%tracer_co2)
!--------------------------------------------------------------------


end subroutine atmos_input_dealloc 


!#####################################################################

 
subroutine microphys_dealloc (Model_microphys)
 
type(microphysics_type), intent(inout) :: Model_microphys
 
!----------------------------------------------------------------------
!   microphys_dealloc calls model_micro_dealloc to deallocate the 
!   array components of the microphysics_type structure Model_microphys.
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
                 'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    deallocate the components of module variable Model_microphys.
!---------------------------------------------------------------------
      call model_micro_dealloc (Model_microphys)

!--------------------------------------------------------------------
 

end subroutine microphys_dealloc


!#####################################################################
! <SUBROUTINE NAME="radiation_driver_end">
!  <OVERVIEW>
!   radiation_driver_end is the destructor for radiation_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_driver_end is the destructor for radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine radiation_driver_end

!----------------------------------------------------------------------
!    radiation_driver_end is the destructor for radiation_driver_mod.
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    write restart file if desired; the file is not necessary if job 
!    ends on step prior to radiation ts, or if restart seamlessness 
!    is not required.
!---------------------------------------------------------------------
    if (using_restart_file) then
! Make sure that the restart_versions variable is up to date.
      vers = restart_versions(size(restart_versions(:)))
      call radiation_driver_restart
    endif

!---------------------------------------------------------------------
!    wrap up modules initialized by this module.
!---------------------------------------------------------------------
      call astronomy_end

!---------------------------------------------------------------------
!    wrap up modules specific to the radiation package in use.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        call cloudrad_package_end
        call aerosolrad_package_end
        call rad_output_file_end
        call sea_esf_rad_end
      else
        call original_fms_rad_end
      endif


!---------------------------------------------------------------------
!    release space for renormalization arrays, if that option is active.
!---------------------------------------------------------------------
!     if (renormalize_sw_fluxes .or. use_hires_coszen .or. &
      if (renormalize_sw_fluxes .or.                       &
          all_step_diagnostics)  then
        deallocate (solar_save, flux_sw_surf_save, sw_heating_save, &
                    dum_idjd,   &
                    flux_sw_surf_dir_save,   &
                    flux_sw_surf_refl_dir_save,   &
                    flux_sw_surf_dif_save,   &
                    flux_sw_down_vis_dir_save,   &
                    flux_sw_down_vis_dif_save,   &
                    flux_sw_down_total_dir_save, &
                    flux_sw_down_total_dif_save, &
                    flux_sw_vis_save, &
                    flux_sw_vis_dir_save, &
                    flux_sw_refl_vis_dir_save, &
                    flux_sw_vis_dif_save, &
                    tot_heating_save, dfsw_save, ufsw_save,   &
                    swdn_special_save, swup_special_save,          &
                    fsw_save, hsw_save)
        if (do_swaerosol_forcing) then
          deallocate (dfsw_ad_save, ufsw_ad_save)
        endif
        if (do_clear_sky_pass) then
          deallocate (sw_heating_clr_save, tot_heating_clr_save,  &
                      dfswcf_save, ufswcf_save, fswcf_save,   &
                      swdn_special_clr_save, swup_special_clr_save,   &
                     flux_sw_down_total_dir_clr_save, &
                     flux_sw_down_total_dif_clr_save, &
                      flux_sw_down_vis_clr_save,   &
                      hswcf_save)
          if (do_swaerosol_forcing) then
            deallocate (dfswcf_ad_save, ufswcf_ad_save)
          endif
        endif
      endif

!---------------------------------------------------------------------
!    release space needed when all_step_diagnostics is active.
!---------------------------------------------------------------------
      if (all_step_diagnostics)  then
        deallocate (olr_save, lwups_save, lwdns_save, flxnet_save, &
                                                           tdtlw_save)
        deallocate (netlw_special_save)
        if (do_lwaerosol_forcing) then
          deallocate (olr_ad_save, lwups_ad_save, lwdns_ad_save)
        endif
        if (do_clear_sky_pass) then
          deallocate (olr_clr_save, lwups_clr_save, lwdns_clr_save, &
                      flxnetcf_save, tdtlw_clr_save)
          deallocate (netlw_special_clr_save)
          if (do_lwaerosol_forcing) then
            deallocate (olr_ad_clr_save, lwups_ad_clr_save,  &
                        lwdns_ad_clr_save)
          endif
        endif
      endif

!---------------------------------------------------------------------
!    release space used for module variables that hold data between
!    timesteps.
!---------------------------------------------------------------------
        deallocate (Rad_output%tdt_rad, Rad_output%tdt_rad_clr,  & 
                    Rad_output%tdtsw, Rad_output%tdtsw_clr, &
                    Rad_output%ufsw, Rad_output%dfsw, &
                    Rad_output%ufsw_clr, Rad_output%dfsw_clr, &
                    Rad_output%tdtlw_clr, &
                    Rad_output%flxnet, Rad_output%flxnetcf, &
                    Rad_output%tdtlw, Rad_output%flux_sw_surf,  &
                    Rad_output%flux_sw_surf_dir,  &
                    Rad_output%flux_sw_surf_refl_dir,  &
                    Rad_output%flux_sw_surf_dif,  &
                    Rad_output%flux_sw_down_vis_dir,  &
                    Rad_output%flux_sw_down_vis_dif,  &
                    Rad_output%flux_sw_down_total_dir,  &
                    Rad_output%flux_sw_down_total_dif,  &
                   Rad_output%flux_sw_down_total_dir_clr,  &
                   Rad_output%flux_sw_down_total_dif_clr,  &
                    Rad_output%flux_sw_down_vis_clr,  &
                    Rad_output%flux_sw_vis,  &
                    Rad_output%flux_sw_vis_dir,  &
                    Rad_output%flux_sw_refl_vis_dir,  &
                    Rad_output%flux_sw_vis_dif,  &
                    Rad_output%flux_lw_surf, Rad_output%coszen_angle)

!----------------------------------------------------------------------
!    deallocate arrays related to the time_varying solar constant.
!----------------------------------------------------------------------
      if (time_varying_solar_constant) then
        deallocate (solflxtot_lean, Solar_spect%solflxband_lean, &
                    Solar_spect%solflxband_lean_ann_1882, &
                    Solar_spect%solflxband_lean_ann_2000)  
      endif

!---------------------------------------------------------------------
!    call rad_utilities_end to uninitialize that module.
!---------------------------------------------------------------------
        call rad_utilities_end

!----------------------------------------------------------------------
!    set initialization status flag.
!----------------------------------------------------------------------
      module_is_initialized = .false.



end subroutine radiation_driver_end



!######################################################################

subroutine return_cosp_inputs (  &
                      is, ie, js, je, donner_meso_is_largescale,  &
                            Time_diag, Atmos_input, stoch_cloud_type, &
                      stoch_conc_drop, stoch_conc_ice, stoch_size_drop,&
                      stoch_size_ice, tau_stoch, lwem_stoch, &
                      Model_microphys, &
                      do_cosp, do_modis_yim, Lsc_microphys)

!---------------------------------------------------------------------
!    subroutine return_cosp_inputs calculates and returns the fields 
!    needed as input by the COSP simulator.
!---------------------------------------------------------------------

integer, intent(in)                      :: is,ie, js, je
logical, intent(in)                      :: donner_meso_is_largescale
logical, intent(in)                      :: do_cosp, do_modis_yim
type(time_type), intent(in)              :: Time_diag
type(atmos_input_type), intent(inout)    :: Atmos_input     
type(microphysics_type), intent(inout)   :: Model_microphys
type(microphysics_type), intent(in)      :: Lsc_microphys
real, dimension(:,:,:,:), intent(inout)  ::    &
                                stoch_cloud_type, stoch_conc_drop, &
                                stoch_conc_ice, stoch_size_drop,  &
                                stoch_size_ice, tau_stoch, lwem_stoch

!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      call obtain_cloud_tau_and_em (is, js, Model_microphys, &
                                    Atmos_input, &
                                    tau_stoch(is:ie,js:je,:,:),  &
                                    lwem_stoch(is:ie,js:je,:,:) )

!-------------------------------------------------------------------
      if (do_cosp) then

!---------------------------------------------------------------------
!    save the stochastic cloud type in each subcolumn.
!    output values of 0 --> no cloud
!           values of 1 --> stratiform cloud
!           values of 2 --> convective cloud
!    input values are 0(none), 1(strat), 2(donnermeso), 3(donnercell), 
!    4(uw)
!---------------------------------------------------------------------
        stoch_cloud_type(is:ie,js:je,:,:) =   &
                         Model_microphys%stoch_cloud_type(:,:,:,:)
         
!---------------------------------------------------------------------
!    donner meso clouds may be treated either as large-scale or
!    convective clouds, dependent on donner_meso_is_largescale.
!---------------------------------------------------------------------
        if (donner_meso_is_largescale) then
          where (stoch_cloud_type(is:ie,js:je,:,:) == 2)
            stoch_cloud_type(is:ie,js:je,:,:) = 1
          end where
          where (stoch_cloud_type(is:ie,js:je,:,:) >= 3)
            stoch_cloud_type(is:ie,js:je,:,:) = 2
          end where
        else
          where (stoch_cloud_type(is:ie,js:je,:,:) >= 2)
            stoch_cloud_type(is:ie,js:je,:,:) = 2
          end where
        endif    

!---------------------------------------------------------------------
!    save the particle concentrations and sizes seen by the radiation
!    package in each stochastic column.
!---------------------------------------------------------------------
        stoch_conc_drop(is:ie,js:je,:,:) =  &
                             Model_microphys%stoch_conc_drop(:,:,:,:)
        stoch_conc_ice (is:ie,js:je,:,:) =  &
                             Model_microphys%stoch_conc_ice (:,:,:,:)
        stoch_size_drop(is:ie,js:je,:,:) =  &
                             Model_microphys%stoch_size_drop(:,:,:,:)
        stoch_size_ice (is:ie,js:je,:,:) =  &
                             Model_microphys%stoch_size_ice (:,:,:,:)

      endif

!-------------------------------------------------------------------
      if (do_modis_yim) then
        call  modis_yim (is, js, Time_diag, Tau_stoch(is:ie,js:je,:,:),&
                         Model_microphys, Atmos_input)
      endif
      call modis_cmip (is, js, Time_diag, Lsc_microphys, &
                       Atmos_input)

!-------------------------------------------------------------------


end subroutine return_cosp_inputs 



!#######################################################################
!#######################################################################
! <SUBROUTINE NAME="radiation_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine radiation_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

! Make sure that the restart_versions variable is up to date.
  vers = restart_versions(size(restart_versions(:)))
  call write_restart_nc(timestamp)

end subroutine radiation_driver_restart
! </SUBROUTINE> NAME="radiation_driver_restart"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------

subroutine write_restart_nc(timestamp)
  character(len=*), intent(in), optional :: timestamp


  if( .not. using_restart_file ) return
!---------------------------------------------------------------------
!    only the root pe will write control information -- the last value 
!    in the list of restart versions and the alarm information.
!---------------------------------------------------------------------
        if (mpp_pe() == mpp_root_pe() ) then
         call error_mesg('radiation_driver_mod', 'Writing netCDF formatted restart file: RESTART/radiation_driver.res.nc', NOTE)
        endif

!---------------------------------------------------------------------
!    write out the optional time average restart data. note that 
!    do_average and renormalize_sw_fluxes may not both be true.
!---------------------------------------------------------------------
        int_renormalize_sw_fluxes = 0
        int_do_clear_sky_pass = 0
        if(renormalize_sw_fluxes) then
          int_renormalize_sw_fluxes = 1
        else if (use_hires_coszen) then
          if (current_sw_zenith_step == nzens_per_sw_rad_timestep) then
            int_renormalize_sw_fluxes = 2
          else
            int_renormalize_sw_fluxes = -2
            call error_mesg ('radiation_driver/write_restart_nc', &
             ' you are writing restart file on a non-radiation &
               &timestep. As a consequence, model results will be &
               &different if the model is run with different restart &
               & intervals. To correct, make sure rad_time_step &
               & is an integral factor of the requested run length.', &
                                                                  NOTE)
          endif
        endif
        if(do_clear_sky_pass) int_do_clear_sky_pass = 1

! Make sure that the restart_versions variable is up to date.
        vers = restart_versions(size(restart_versions(:)))
        call save_restart(Rad_restart, timestamp)
        if(in_different_file) call save_restart(Til_restart, timestamp)

end subroutine write_restart_nc


!#####################################################################
subroutine rad_driver_register_restart(fname)
  character(len=*), intent(in) :: fname
  character(len=64)            :: fname2
  integer                      :: id_restart

   call get_mosaic_tile_file(fname, fname2, .false. ) 
   allocate(Rad_restart)
   if(trim(fname2) == trim(fname)) then
      Til_restart => Rad_restart
      in_different_file = .false.
   else
      in_different_file = .true.
      allocate(Til_restart)
   endif

  id_restart = register_restart_field(Rad_restart, fname, 'vers', vers, no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'lwrad_alarm', lwrad_alarm, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'swrad_alarm', swrad_alarm, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'lw_rad_time_step', lw_rad_time_step, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'sw_rad_time_step', sw_rad_time_step, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'renormalize_sw_fluxes', int_renormalize_sw_fluxes,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'do_clear_sky_pass', int_do_clear_sky_pass,no_domain=.true.)
  id_restart = register_restart_field(Til_restart, fname, 'tdt_rad', Rad_output%tdt_rad(:,:,:,1) )
  id_restart = register_restart_field(Til_restart, fname, 'tdtlw', Rad_output%tdtlw)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf', Rad_output%flux_sw_surf)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_dir', Rad_output%flux_sw_surf_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_refl_dir', Rad_output%flux_sw_surf_refl_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_dif', Rad_output%flux_sw_surf_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_dir', Rad_output%flux_sw_down_vis_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_dif', Rad_output%flux_sw_down_vis_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dir', Rad_output%flux_sw_down_total_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dif', Rad_output%flux_sw_down_total_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis', Rad_output%flux_sw_vis)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_dir', Rad_output%flux_sw_vis_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_refl_vis_dir', Rad_output%flux_sw_refl_vis_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_dif', Rad_output%flux_sw_vis_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_lw_surf', Rad_output%flux_lw_surf)
  id_restart = register_restart_field(Til_restart, fname, 'coszen_angle', Rad_output%coszen_angle)
  if (renormalize_sw_fluxes ) then   
     id_restart = register_restart_field(Til_restart, fname, 'solar_save', solar_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_save', flux_sw_surf_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_dir_save', flux_sw_surf_dir_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_refl_dir_save', flux_sw_surf_refl_dir_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_dif_save', flux_sw_surf_dif_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_dir_save', flux_sw_down_vis_dir_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_dif_save', flux_sw_down_vis_dif_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dir_save', flux_sw_down_total_dir_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dif_save', flux_sw_down_total_dif_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_save', flux_sw_vis_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_dir_save', flux_sw_vis_dir_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_refl_vis_dir_save', flux_sw_refl_vis_dir_save)
     id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_dif_save', flux_sw_vis_dif_save)
     id_restart = register_restart_field(Til_restart, fname, 'sw_heating_save', sw_heating_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'tot_heating_save', tot_heating_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'dfsw_save', dfsw_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'ufsw_save', ufsw_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'fsw_save', fsw_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'hsw_save', hsw_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'swdn_special_save', swdn_special_save(:,:,:,1))
     id_restart = register_restart_field(Til_restart, fname, 'swup_special_save', swup_special_save(:,:,:,1))
     if (do_clear_sky_pass) then
        id_restart = register_restart_field(Til_restart, fname, 'sw_heating_clr_save', sw_heating_clr_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'tot_heating_clr_save', tot_heating_clr_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'dfswcf_save', dfswcf_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'ufswcf_save', ufswcf_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'fswcf_save', fswcf_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'hswcf_save', hswcf_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dir_clr_save', &
             flux_sw_down_total_dir_clr_save(:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dif_clr_save', &
             flux_sw_down_total_dif_clr_save(:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_clr_save', &
             flux_sw_down_vis_clr_save(:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'swdn_special_clr_save', swdn_special_clr_save(:,:,:,1))
        id_restart = register_restart_field(Til_restart, fname, 'swup_special_clr_save', swup_special_clr_save(:,:,:,1))
     endif
  endif

end subroutine rad_driver_register_restart


!#####################################################################
! <SUBROUTINE NAME="read_restart_nc">
!  <OVERVIEW>
!    read_restart_nc reads a netcdf restart file containing radiation
!    restart information.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_restart_nc reads a netcdf restart file containing radiation
!    restart information.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_restart_nc
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_restart_nc

  character(len=64) :: fname='INPUT/radiation_driver.res.nc'
  real              :: flag1, flag2
  logical           :: renorm_present, cldfree_present
  integer           :: new_rad_time
  integer           :: lw_old_time_step, sw_old_time_step
  integer           :: siz(4)
!----------------------------------------------------------------------
!    when running in gcm, read a restart file. this is not done in the
!    standalone case.
!---------------------------------------------------------------------
  if (mpp_pe() == mpp_root_pe() ) then
    call error_mesg('radiation_driver_mod', 'Reading netCDF formatted restart file: INPUT/radiation_driver.res.nc', NOTE)
  endif
  call read_data(fname, 'vers', vers, no_domain=.true.)

!--------------------------------------------------------------------
  if (field_exist (fname, 'rad_alarm')) then
    call read_data(fname, 'rad_alarm', lwrad_alarm,no_domain=.true.)
    call read_data(fname, 'rad_alarm', swrad_alarm,no_domain=.true.)
  else
    call read_data(fname, 'lwrad_alarm', lwrad_alarm,no_domain=.true.)
    call read_data(fname, 'swrad_alarm', swrad_alarm,no_domain=.true.)
  endif

!--------------------------------------------------------------------
  if (field_exist (fname, 'rad_time_step')) then
    call read_data(fname, 'rad_time_step', lw_old_time_step,no_domain=.true.)
    call read_data(fname, 'rad_time_step', sw_old_time_step,no_domain=.true.)
  else
  call read_data(fname, 'sw_rad_time_step', sw_old_time_step,no_domain=.true.)
  call read_data(fname, 'lw_rad_time_step', lw_old_time_step,no_domain=.true.)
  endif


  call read_data(fname, 'renormalize_sw_fluxes', flag1,no_domain=.true.)
  call read_data(fname, 'do_clear_sky_pass', flag2,no_domain=.true.)
  renorm_present = .false.
  cldfree_present = .false.
  if (flag1 ==   1.0)  then
    renorm_present = .true.
  else if (flag1 == 2.0) then
  else if (flag1 == -2.0) then
    if (.not. allow_nonrepro_across_restarts) then
      call error_mesg ( 'radiation_driver/read_restart_nc', &
       'the restart was written on a non-radiation step, so model&
        & solution will NOT be independent of restart interval. If &
        & you dont care about this, set nml variable &
        &allow_nonrepro_across_restarts to .true. and resubmit; &
        &if you do, contact developer so additional code may be &
        &added to allow seamless restart, OR rerun last job segment &
        & so that it is an integral number of rad_time_steps  &
                                                       &long.', FATAL)
    else
      call error_mesg ( 'radiation_driver/read_restart_nc', &
       'the restart was written on a non-radiation step, so model&
        & solution will NOT be independent of restart interval. You &
        & have chosen to proceed anyway by setting nml variable &
        &allow_nonrepro_across_restarts to .true.', NOTE )
      swrad_alarm = 1
      lwrad_alarm = 1
    endif

  endif
    
  if(flag2 .EQ. 1.0) cldfree_present = .true.

!---------------------------------------------------------------------
!    read the restart data.
!    currently this need not be done when hires_coszen = .true.
!---------------------------------------------------------------------
      if (flag1 == 0.0 .or. flag1 == 1.0) then
        call read_data (fname, 'tdt_rad',                Rad_output%tdt_rad(:,:,:,1))
        call read_data (fname, 'tdtlw',                  Rad_output%tdtlw)
        call read_data (fname, 'flux_sw_surf',           Rad_output%flux_sw_surf(:,:,1))
        call read_data (fname, 'flux_sw_surf_dir',       Rad_output%flux_sw_surf_dir(:,:,1))
        if (field_exist (fname, 'flux_sw_surf_refl_dir')) then
          call read_data (fname, 'flux_sw_surf_refl_dir',       Rad_output%flux_sw_surf_refl_dir(:,:,1))
        else
          Rad_output%flux_sw_surf_refl_dir(:,:,1) = 0.0
        endif
        call read_data (fname, 'flux_sw_surf_dif',       Rad_output%flux_sw_surf_dif(:,:,1))
        call read_data (fname, 'flux_sw_down_vis_dir',   Rad_output%flux_sw_down_vis_dir(:,:,1))
        call read_data (fname, 'flux_sw_down_vis_dif',   Rad_output%flux_sw_down_vis_dif(:,:,1))
        call read_data (fname, 'flux_sw_down_total_dir', Rad_output%flux_sw_down_total_dir(:,:,1))
        call read_data (fname, 'flux_sw_down_total_dif', Rad_output%flux_sw_down_total_dif(:,:,1))
        call read_data (fname, 'flux_sw_vis',            Rad_output%flux_sw_vis(:,:,1))
        call read_data (fname, 'flux_sw_vis_dir',        Rad_output%flux_sw_vis_dir(:,:,1))
        if (field_exist (fname, 'flux_sw_refl_vis_dir')) then
          call read_data (fname, 'flux_sw_refl_vis_dir',        Rad_output%flux_sw_refl_vis_dir(:,:,1))
        else
          Rad_output%flux_sw_refl_vis_dir(:,:,1) = 0.0
       endif
        call read_data (fname, 'flux_sw_vis_dif',        Rad_output%flux_sw_vis_dif(:,:,1))
        call read_data (fname, 'flux_lw_surf',           Rad_output%flux_lw_surf)
        call read_data (fname, 'coszen_angle',           Rad_output%coszen_angle)
   endif

!---------------------------------------------------------------------
!    read the optional shortwave renormalization data. 
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes ) then   
         if(renorm_present) then
          call read_data (fname, 'solar_save', solar_save)
          call read_data (fname, 'flux_sw_surf_save', flux_sw_surf_save(:,:,1))
          call read_data (fname, 'flux_sw_surf_dir_save', flux_sw_surf_dir_save(:,:,1))
          if (field_exist (fname, 'flux_sw_surf_refl_dir_save')) then
            call read_data (fname, 'flux_sw_surf_refl_dir_save', flux_sw_surf_refl_dir_save(:,:,1))
          else
            flux_sw_surf_refl_dir_save(:,:,1) = 0.
          endif
          call read_data (fname, 'flux_sw_surf_dif_save', flux_sw_surf_dif_save(:,:,1))
          call read_data (fname, 'flux_sw_down_vis_dir_save', flux_sw_down_vis_dir_save(:,:,1))
          call read_data (fname, 'flux_sw_down_vis_dif_save', flux_sw_down_vis_dif_save(:,:,1))
          call read_data (fname, 'flux_sw_down_total_dir_save', flux_sw_down_total_dir_save(:,:,1))
          call read_data (fname, 'flux_sw_down_total_dif_save', flux_sw_down_total_dif_save(:,:,1))
          call read_data (fname, 'flux_sw_vis_save', flux_sw_vis_save(:,:,1))
          call read_data (fname, 'flux_sw_vis_dir_save', flux_sw_vis_dir_save(:,:,1))
          if (field_exist (fname, 'flux_sw_refl_vis_dir_save')) then
            call read_data (fname, 'flux_sw_refl_vis_dir_save', flux_sw_refl_vis_dir_save(:,:,1))
          else
            flux_sw_refl_vis_dir_save(:,:,1) = 0.0
          endif
          call read_data (fname, 'flux_sw_vis_dif_save', flux_sw_vis_dif_save(:,:,1))
          call read_data (fname, 'sw_heating_save', sw_heating_save(:,:,:,1))
          call read_data (fname, 'tot_heating_save', tot_heating_save(:,:,:,1))
          call read_data (fname, 'dfsw_save', dfsw_save(:,:,:,1)) 
          call read_data (fname, 'ufsw_save', ufsw_save(:,:,:,1)) 
          call read_data (fname, 'fsw_save', fsw_save(:,:,:,1))  
          call read_data (fname, 'hsw_save', hsw_save(:,:,:,1))
          call field_size (fname, 'swdn_special_save', siz)
          call read_data (fname, 'swdn_special_save', swdn_special_save(:,:,1:siz(3),1))
          call read_data (fname, 'swup_special_save', swup_special_save(:,:,1:siz(3),1))
          if (MX_SPEC_LEVS > siz(3)) then
            swdn_special_save(:,:,siz(3)+1:MX_SPEC_LEVS,1) = 0.
            swup_special_save(:,:,siz(3)+1:MX_SPEC_LEVS,1) = 0.
          endif
          if (do_clear_sky_pass) then
           if(cldfree_present) then
            call read_data (fname, 'sw_heating_clr_save', sw_heating_clr_save(:,:,:,1))
            call read_data (fname, 'tot_heating_clr_save', tot_heating_clr_save(:,:,:,1))
            call read_data (fname, 'dfswcf_save', dfswcf_save(:,:,:,1)) 
            call read_data (fname, 'ufswcf_save', ufswcf_save(:,:,:,1)) 
            call read_data (fname, 'fswcf_save', fswcf_save(:,:,:,1))  
            call read_data (fname, 'hswcf_save', hswcf_save(:,:,:,1))
            if (vers >= 10) then
              call read_data (fname, 'flux_sw_down_total_dir_clr_save', flux_sw_down_total_dir_clr_save(:,:,1))
              call read_data (fname, 'flux_sw_down_total_dif_clr_save', flux_sw_down_total_dif_clr_save(:,:,1))
            else
              flux_sw_down_total_dir_clr_save = 0.0
              flux_sw_down_total_dif_clr_save = 0.0
            endif
            if (vers >= 11) then
              call read_data (fname, 'flux_sw_down_vis_clr_save', flux_sw_down_vis_clr_save(:,:,1))
            else
              flux_sw_down_vis_clr_save = 0.0
            endif
            call field_size (fname, 'swdn_special_clr_save', siz)
            call read_data (fname, 'swdn_special_clr_save', swdn_special_clr_save(:,:,1:siz(3),1))
            call read_data (fname, 'swup_special_clr_save', swup_special_clr_save(:,:,1:siz(3),1))
            if (MX_SPEC_LEVS > siz(3)) then
              swdn_special_clr_save(:,:,siz(3)+1:MX_SPEC_LEVS,1) = 0.
              swup_special_clr_save(:,:,siz(3)+1:MX_SPEC_LEVS,1) = 0.
            endif
            endif
         endif
      endif  ! (do_clear_sky_pass)
     endif  ! (renormalize_sw_fluxes)
!----------------------------------------------------------------------
!    if all_step_diagnostics is active and rad_alarm is not 1, abort 
!    job with error message. all_step_diagnostics may only be activated
!    when radiation is to be calculated on the first step of a job, 
!    unless additional arrays are added to the radiation restart file.
!----------------------------------------------------------------------
   if (lwrad_alarm /= 1 .and. all_step_diagnostics) then
     if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('radiation_driver_mod', &
       'cannot set all_step_diagnostics to be .true. unless &
         & starting job on step just prior to radiation call; &
         &doing so will lead to non-reproducibility of restarts', &
                                                                 FATAL)
     endif
  endif
   if (swrad_alarm /= 1 .and. all_step_diagnostics) then
     if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('radiation_driver_mod', &
       'cannot set all_step_diagnostics to be .true. unless &
         & starting job on step just prior to radiation call; &
         &doing so will lead to non-reproducibility of restarts', &
                                                                 FATAL)
     endif
  endif

   if (lwrad_alarm /= 1 .and.    &
                (do_lwaerosol_forcing  .or. do_swaerosol_forcing)) then
     if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('radiation_driver_mod', &
       'aerosol forcing diagnostics will only be strictly valid &
      &when restarting a job on the step just prior to radiation&
      &call; not doing so will lead to invalid diagnostics between time&
      & of restart and next radiation calculation, since these fields &
      &are not saved in the restart file', FATAL)
     endif
  endif
   if (swrad_alarm /= 1 .and.    &
                (do_lwaerosol_forcing  .or. do_swaerosol_forcing)) then
     if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('radiation_driver_mod', &
       'aerosol forcing diagnostics will only be strictly valid &
      &when restarting a job on the step just prior to radiation&
      &call; not doing so will lead to invalid diagnostics between time&
      & of restart and next radiation calculation, since these fields &
      &are not saved in the restart file', FATAL)
     endif
  endif

  !----------------------------------------------------------------------
  !    adjust radiation alarm if radiation step has changed from restart 
  !    file value, if it has not already been set to the first step.
  !----------------------------------------------------------------------
  if (lwrad_alarm /= 1) then
     !     if (rad_alarm == 1) then
     !       if (mpp_pe() == mpp_root_pe() ) then
     !         call error_mesg ('radiation_driver_mod',          &
     !              'radiation will be called on first step of run', NOTE)
     !       endif
     !     else
     if (rad_time_step /= lw_old_time_step ) then
        new_rad_time = lwrad_alarm - lw_old_time_step + lw_rad_time_step
        if ( new_rad_time > 0 ) then
           if (mpp_pe() == mpp_root_pe() ) then
              print *, 'radiation time step has changed, therefore '//&
                   'next time to next do lw radiation also changed;  &
                   &new lwrad_alarm is', new_rad_time
           endif
           lwrad_alarm = new_rad_time
        else
           lwrad_alarm = 1
           if (mpp_pe() == mpp_root_pe() ) then
              call error_mesg ('radiation_driver_mod', &
                   ' radiation to be calculated on first step: lw radiation &
                   &timestep has gotten shorter and is past due', NOTE)
           endif
        endif
     endif
  endif   ! (lwrad_alarm == 1)
  if (swrad_alarm /= 1) then
     !     if (rad_alarm == 1) then
     !       if (mpp_pe() == mpp_root_pe() ) then
     !         call error_mesg ('radiation_driver_mod',          &
     !              'radiation will be called on first step of run', NOTE)
     !       endif
     !     else
     if (sw_rad_time_step /= sw_old_time_step ) then
        new_rad_time = swrad_alarm - sw_old_time_step + sw_rad_time_step
        if ( new_rad_time > 0 ) then
           if (mpp_pe() == mpp_root_pe() ) then
              print *, 'radiation time step has changed, therefore '//&
                   'next time to next do sw radiation also changed;  &
                   &new swrad_alarm is', new_rad_time
           endif
           swrad_alarm = new_rad_time
        else
           swrad_alarm = 1
           if (mpp_pe() == mpp_root_pe() ) then
              call error_mesg ('radiation_driver_mod', &
                   ' radiation to be calculated on first step: sw radiation &
                   &timestep has gotten shorter and is past due', NOTE)
           endif
        endif
     endif
  endif   ! (swrad_alarm == 1)

  vers = restart_versions(size(restart_versions(:)))

end subroutine read_restart_nc


!#####################################################################
! <SUBROUTINE NAME="initialize_diagnostic_integrals">
!  <OVERVIEW>
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call initialize_diagnostic_integrals
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine initialize_diagnostic_integrals

!---------------------------------------------------------------------
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    initialize standard global quantities for integral package. 
!----------------------------------------------------------------------
      call diag_integral_field_init ('olr',    std_digits)
      call diag_integral_field_init ('abs_sw', std_digits)
!     call diag_integral_field_init ('olr_clr',    std_digits)
!     call diag_integral_field_init ('abs_sw_clr', std_digits)

!----------------------------------------------------------------------
!    if hemispheric integrals and global integrals with extended signif-
!    icance are desired, inform diag_integrals_mod.
!----------------------------------------------------------------------
      if (calc_hemi_integrals) then
        call diag_integral_field_init ('sntop_tot_sh', extra_digits)
        call diag_integral_field_init ('lwtop_tot_sh', extra_digits)
        call diag_integral_field_init ('sngrd_tot_sh', extra_digits)
        call diag_integral_field_init ('lwgrd_tot_sh', extra_digits)
        call diag_integral_field_init ('sntop_tot_nh', extra_digits)
        call diag_integral_field_init ('lwtop_tot_nh', extra_digits)
        call diag_integral_field_init ('sngrd_tot_nh', extra_digits)
        call diag_integral_field_init ('lwgrd_tot_nh', extra_digits)
        call diag_integral_field_init ('sntop_tot_gl', extra_digits)
        call diag_integral_field_init ('lwtop_tot_gl', extra_digits)
        call diag_integral_field_init ('sngrd_tot_gl', extra_digits)
        call diag_integral_field_init ('lwgrd_tot_gl', extra_digits)

!---------------------------------------------------------------------
!    if clear-sky integrals are desired, include them.
!---------------------------------------------------------------------
        if (do_clear_sky_pass) then
          call diag_integral_field_init ('sntop_clr_sh', extra_digits)
          call diag_integral_field_init ('lwtop_clr_sh', extra_digits)
          call diag_integral_field_init ('sngrd_clr_sh', extra_digits)
          call diag_integral_field_init ('lwgrd_clr_sh', extra_digits)
          call diag_integral_field_init ('sntop_clr_nh', extra_digits)
          call diag_integral_field_init ('lwtop_clr_nh', extra_digits)
          call diag_integral_field_init ('sngrd_clr_nh', extra_digits)
          call diag_integral_field_init ('lwgrd_clr_nh', extra_digits)
          call diag_integral_field_init ('sntop_clr_gl', extra_digits)
          call diag_integral_field_init ('lwtop_clr_gl', extra_digits)
          call diag_integral_field_init ('sngrd_clr_gl', extra_digits)
          call diag_integral_field_init ('lwgrd_clr_gl', extra_digits)
        endif
      endif

!--------------------------------------------------------------------


end subroutine initialize_diagnostic_integrals



!#######################################################################
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init ( Time, axes )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   Current time
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
! </SUBROUTINE>
!
subroutine diag_field_init ( Time, axes )

!---------------------------------------------------------------------
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)

!--------------------------------------------------------------------
!  intent(in) variables
!
!      Time        current time
!      axes        data axes for use with diagnostic fields
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      character(len=8)  ::   clr
      character(len=16) ::   clr2, lwaer_prep, swaer_prep
      integer           ::   bxes(4)
      integer           ::   i, n

!--------------------------------------------------------------------
!  local variables:
!
!       clr          character string used in netcdf variable short name
!       clr2         character string used in netcdf variable long name
!       n            number of passes through name generation loop
!       i            do-loop index
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    define variable axis array with elements (1:3) valid for variables
!    defined at flux levels.
!-------------------------------------------------------------------
      bxes(1:2) = axes(1:2)
      bxes(3) = axes(4)
      bxes(4) = axes(4)

!---------------------------------------------------------------------
!    determine how many passes are needed through the name generation 
!    loop. 
!---------------------------------------------------------------------
      if (do_clear_sky_pass) then
        n= 2
      else
        n= 1
      endif

      if (Sw_control%do_swaerosol ) then
        swaer_prep = 'without'
      else
        swaer_prep = 'with'
      endif
      if (Lw_control%do_lwaerosol ) then
        lwaer_prep = 'without'
      else
        lwaer_prep = 'with'
      endif

!---------------------------------------------------------------------
!    generate names for standard and clear sky diagnostic fields. if 
!    clear sky values being generated, generate the clear sky names
!    on pass 1, followed by the standard names.
!---------------------------------------------------------------------
      do i = 1, n
        if ( i == n) then
          clr  = "    "
          clr2 = "          "
        else
          clr  = "_clr"
          clr2 = "clear sky "
        endif

        id_swdn_special(1,i) = register_diag_field (mod_name,   &
                'swdn_200hPa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at 200 hPa', &
                'watts/m2', missing_value=missing_value)
        id_swdn_special(2,i) = register_diag_field (mod_name,   &
                'swdn_lin_trop'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at linear tropopause', &
                'watts/m2', missing_value=missing_value)
        id_swdn_special(3,i) = register_diag_field (mod_name,   &
                'swdn_therm_trop'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at thermo tropopause', &
                'watts/m2', missing_value=missing_value)
        id_swdn_special(4,i) = register_diag_field (mod_name,   &
                'swdn_1_Pa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at 1 Pa', &
                'watts/m2', missing_value=missing_value)

        id_swup_special(1,i) = register_diag_field (mod_name,   &
                'swup_200hPa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at 200 hPa', &
                'watts/m2', missing_value=missing_value)
        id_swup_special(2,i) = register_diag_field (mod_name,   &
                'swup_lin_trop'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at linear tropopause', &
                'watts/m2', missing_value=missing_value)
        id_swup_special(3,i) = register_diag_field (mod_name,   &
                'swup_therm_trop'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at thermo tropopause', &
                'watts/m2', missing_value=missing_value)
        id_swup_special(4,i) = register_diag_field (mod_name,   &
                'swup_1_Pa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at 1 Pa', &
                'watts/m2', missing_value=missing_value)

        id_netlw_special(1,i) = register_diag_field (mod_name,   &
                'netlw_200hPa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net LW flux at 200 hPa', &
                'watts/m2', missing_value=missing_value)
        id_netlw_special(2,i) = register_diag_field (mod_name,   &
                'netlw_lin_trop'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net LW flux at linear tropopause', &
                'watts/m2', missing_value=missing_value)
        id_netlw_special(3,i) = register_diag_field (mod_name,   &
                'netlw_therm_trop'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net LW flux at thermo tropopause', &
                'watts/m2', missing_value=missing_value)
        id_netlw_special(4,i) = register_diag_field (mod_name,   &
                'netlw_1_Pa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net LW flux at 1 Pa', &
                'watts/m2', missing_value=missing_value)

        id_tdt_sw(i) = register_diag_field (mod_name,   &
                'tdt_sw'//trim(clr), axes(1:3), Time, & 
                trim(clr2)//'temperature tendency for SW radiation', &
                'deg_K/sec', missing_value=missing_value) 

        id_ufsw(i) = register_diag_field (mod_name,   &
               'allufsw'//trim(clr), bxes(1:3), Time, &
               trim(clr2)//'upward sw flux', &
               'watts/m2', missing_value=missing_value)

        id_dfsw(i) = register_diag_field (mod_name,   &
               'alldfsw'//trim(clr), bxes(1:3), Time, &
               trim(clr2)//'downward sw flux', &
               'watts/m2', missing_value=missing_value)

        id_flxnet(i) = register_diag_field (mod_name,   &
               'allnetlw'//trim(clr), bxes(1:3), Time, &
               trim(clr2)//'net lw flux', &
               'watts/m2', missing_value=missing_value)

        id_tdt_lw(i) = register_diag_field (mod_name,    &
                'tdt_lw'//trim(clr), axes(1:3), Time, &
                trim(clr2)//'temperature tendency for LW radiation', &
                'deg_K/sec', missing_value=missing_value)

        id_swdn_toa(i) = register_diag_field (mod_name,   &
                'swdn_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at TOA', &
                'watts/m2', missing_value=missing_value)

        id_swup_toa(i) = register_diag_field (mod_name,    &
                'swup_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at TOA', &
                'watts/m2', missing_value=missing_value)

        id_olr(i) = register_diag_field (mod_name,   &
                'olr'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'outgoing longwave radiation', &
                'watts/m2', missing_value=missing_value)

        id_netrad_toa(i) = register_diag_field (mod_name,   &
                'netrad_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net radiation (lw + sw) at toa', &
                'watts/m2', missing_value=missing_value)

        id_netrad_1_Pa(i) = register_diag_field (mod_name,   &
                'netrad_1_Pa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net radiation (lw + sw) at 1 Pa', &
                'watts/m2', missing_value=missing_value)

        id_swup_sfc(i) = register_diag_field (mod_name,    &
                'swup_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at surface', &
                'watts/m2', missing_value=missing_value)

        id_swdn_sfc(i) = register_diag_field (mod_name,     &
                'swdn_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at surface', &
                'watts/m2', missing_value=missing_value)

        id_lwup_sfc(i) = register_diag_field (mod_name,   &
                'lwup_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'LW flux up at surface', &
                'watts/m2', missing_value=missing_value)

        id_lwdn_sfc(i) = register_diag_field (mod_name,    &
                'lwdn_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'LW flux down at surface', &
                'watts/m2', missing_value=missing_value)

        id_swtoa(i) = register_diag_field (mod_name,    &
                 'swtoa'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net SW flux at TOA ', &
                  'watts/m2', missing_value=missing_value)

        id_swsfc(i) = register_diag_field (mod_name,    &
                  'swsfc'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net SW flux at surface', &
                  'watts/m2', missing_value=missing_value)
 
        id_lwsfc(i) = register_diag_field (mod_name,    &
                  'lwsfc'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' Net LW flux at surface', &
                  'watts/m2', missing_value=missing_value)
 
        id_swtoa_ad(i) = register_diag_field (mod_name,    &
                  'swtoa_ad'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net SW flux at TOA '// trim(swaer_prep) &
                                 // ' aerosol', &
                  'watts/m2', missing_value=missing_value)
 
        id_swsfc_ad(i) = register_diag_field (mod_name,    &
                 'swsfc_ad'//trim(clr), axes(1:2), Time, &
            trim(clr2)//' Net SW flux at surface '// trim(swaer_prep) &
           // ' aerosol', &
                'watts/m2', missing_value=missing_value)

       id_swdn_sfc_ad(i) = register_diag_field (mod_name,    &
                 'swdn_sfc_ad'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' SW flux down at surface '// &
                   trim(swaer_prep) // ' aerosol', &
                 'watts/m2', missing_value=missing_value)
 
        id_swup_sfc_ad(i) = register_diag_field (mod_name,    &
                 'swup_sfc_ad'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' SW flux up at surface ' //   &
                   trim(swaer_prep) // ' aerosol', &
                  'watts/m2', missing_value=missing_value)

        id_swup_toa_ad(i) = register_diag_field (mod_name,    &
                 'swup_toa_ad'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' SW flux up at TOA '  //  &
                   trim(swaer_prep) // ' aerosol', &
                 'watts/m2', missing_value=missing_value)
 
         id_olr_ad(i) = register_diag_field (mod_name,    &
                  'lwtoa_ad'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net LW flux at TOA (olr) ' //  &
                   trim(lwaer_prep) // ' aerosol', &
                  'watts/m2', missing_value=missing_value)

         id_lwsfc_ad(i) = register_diag_field (mod_name,    &
                  'lwsfc_ad'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net LW flux at surface  ' //   &
                   trim(lwaer_prep) // ' aerosol', &
                 'watts/m2', missing_value=missing_value)
 
       end do

         id_allradp   = register_diag_field (mod_name,   &
                 'allradp', axes(1:3), Time, &
                 'temperature tendency for SW + LW radiation', &
                 'deg_K/sec', missing_value=missing_value)

!----------------------------------------------------------------------
!    register fields that are not clear-sky depedent.
!----------------------------------------------------------------------
        id_conc_drop = register_diag_field (mod_name,   &
                   'conc_drop', axes(1:3), Time, & 
                   'drop concentration ', &
                   'g/m^3', missing_value=missing_value) 
    
        id_conc_ice = register_diag_field (mod_name,   &
                   'conc_ice', axes(1:3), Time, & 
                   'ice concentration ', &
                   'g/m^3', missing_value=missing_value) 
 
      
      id_flux_sw_dir = register_diag_field (mod_name,    &
                'flux_sw_dir', axes(1:2), Time, &
                'net direct sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_refl_dir = register_diag_field (mod_name,    &
                'flux_sw_refl_dir', axes(1:2), Time, &
                'refl sw dir from sfc', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_refl_vis_dir = register_diag_field (mod_name,    &
                'flux_sw_refl_vis_dir', axes(1:2), Time, &
                'refl sw vis dir from sfc', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_dif = register_diag_field (mod_name,    &
                'flux_sw_dif', axes(1:2), Time, &
                'net diffuse sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw     = register_diag_field (mod_name,    &
                'flux_sw', axes(1:2), Time, &
                'net sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_down_vis_dir = register_diag_field (mod_name,    &
                'flux_sw_down_vis_dir', axes(1:2), Time, &
                'downward direct visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_down_vis_dif = register_diag_field (mod_name,    &
                'flux_sw_down_vis_dif', axes(1:2), Time, &
                'downward diffuse visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_down_total_dir = register_diag_field (mod_name,    &
                'flux_sw_down_total_dir', axes(1:2), Time, &
                'downward direct total sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_down_total_dif = register_diag_field (mod_name,    &
               'flux_sw_down_total_dif', axes(1:2), Time, &
               'downward diffuse total sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

    if (do_clear_sky_pass) then

      id_flux_sw_down_total_dir_clr = register_diag_field (mod_name,  &
               'flux_sw_down_total_dir_clr', axes(1:2), Time, &
               'downward clearsky direct total sfc sw flux',  &
               'watts/m2',  missing_value=missing_value)
  
      id_flux_sw_down_total_dif_clr = register_diag_field (mod_name,  &
               'flux_sw_down_total_dif_clr', axes(1:2), Time, &
               'downward clearsky diffuse total sfc sw flux',  &
               'watts/m2', missing_value=missing_value)

      id_flux_sw_down_vis_clr = register_diag_field (mod_name,    &
                'flux_sw_down_vis_clr', axes(1:2), Time, &
                'downward visible sfc sw flux clear sky', 'watts/m2', &
                 missing_value=missing_value)

    endif 

      id_flux_sw_vis = register_diag_field (mod_name,    &
               'flux_sw_vis', axes(1:2), Time, &
               'net visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_vis_dir = register_diag_field (mod_name,    &
               'flux_sw_vis_dir', axes(1:2), Time, &
               'net direct visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_vis_dif = register_diag_field (mod_name,    &
                'flux_sw_vis_dif', axes(1:2), Time, &
                'net diffuse visible sfc sw flux', 'watts/m2', &
                  missing_value=missing_value)


      id_sol_con = register_diag_field (mod_name,    &
                  'solar_constant', Time, &
                  'solar constant', 'watts/m2', &
                  missing_value=missing_value)      
                           
      id_co2_tf = register_diag_field (mod_name,    &
                  'co2_tf', Time, &
                  'co2 mixing ratio used for tf calculation', 'ppmv', &
                  missing_value=missing_value)      
                           
      id_ch4_tf = register_diag_field (mod_name,    &
                  'ch4_tf', Time, &
                  'ch4 mixing ratio used for tf calculation', 'ppbv', &
                  missing_value=missing_value)      
                           
      id_n2o_tf = register_diag_field (mod_name,    &
                  'n2o_tf', Time, &
                  'n2o mixing ratio used for tf calculation', 'ppbv', &
                  missing_value=missing_value)      
                           
      id_rrvco2 = register_diag_field (mod_name,    &
                  'rrvco2', Time, &
                  'co2 mixing ratio', 'ppmv', &
                  missing_value=missing_value)      
                           
      id_rrvf11 = register_diag_field (mod_name,    &
                  'rrvf11', Time, &
                  'f11 mixing ratio', 'pptv', &
                  missing_value=missing_value)
        
      id_rrvf12 = register_diag_field (mod_name,    &
                  'rrvf12', Time, &
                  'f12 mixing ratio', 'pptv', &
                  missing_value=missing_value)

      id_rrvf113 = register_diag_field (mod_name,    &
                   'rrvf113', Time, &
                   'f113 mixing ratio', 'pptv', &
                   missing_value=missing_value)
 
       id_rrvf22 = register_diag_field (mod_name,    &
                   'rrvf22', Time, &
                   'f22 mixing ratio', 'pptv', &
                   missing_value=missing_value)

       id_rrvch4 = register_diag_field (mod_name,    &
                   'rrvch4', Time, &
                   'ch4 mixing ratio', 'ppbv', &
                   missing_value=missing_value)

       id_rrvn2o = register_diag_field (mod_name,    &
                   'rrvn2o', Time, &
                   'n2o mixing ratio', 'ppbv', &
                   missing_value=missing_value)

         id_alb_sfc_avg = register_diag_field (mod_name,    &
                 'averaged_alb_sfc', axes(1:2), Time, &
                 'surface albedo', 'percent', &
                 missing_value=missing_value)
         if (id_alb_sfc_avg > 0) then
           allocate (swdns_acc(id,jd))
           allocate (swups_acc(id,jd))
           swups_acc = 0.0
           swdns_acc = 1.0e-35
         endif
      id_alb_sfc = register_diag_field (mod_name,    &
                'alb_sfc', axes(1:2), Time, &
                'surface albedo', 'percent', &
                  missing_value=missing_value) 

      id_alb_sfc_vis_dir = register_diag_field (mod_name,    &
                'alb_sfc_vis_dir', axes(1:2), Time, &
!               'surface albedo_vis_dir', 'percent')
! BUGFIX
                'surface albedo_vis_dir', 'percent', &
                 missing_value=missing_value)
      id_alb_sfc_nir_dir = register_diag_field (mod_name,    &
                'alb_sfc_nir_dir', axes(1:2), Time, &
!               'surface albedo_nir', 'percent')
! BUGFIX
                'surface albedo_nir_dir', 'percent', &
                  missing_value=missing_value)
 
      id_alb_sfc_vis_dif = register_diag_field (mod_name,    &
                 'alb_sfc_vis_dif', axes(1:2), Time, &
!               'surface albedo_vis', 'percent')
! BUGFIX
                 'surface albedo_vis_dif', 'percent', &
                  missing_value=missing_value)
      id_alb_sfc_nir_dif = register_diag_field (mod_name,    &
                 'alb_sfc_nir_dif', axes(1:2), Time, &
!               'surface albedo_nir', 'percent')
! BUGFIX
                 'surface albedo_nir_dif', 'percent', &
                   missing_value=missing_value)
      id_cosz = register_diag_field (mod_name,    &
                'cosz',axes(1:2),  Time,    &
                'cosine of zenith angle',    &
                'none', missing_value=missing_value)

      id_fracday = register_diag_field (mod_name,   &
                'fracday',axes(1:2), Time,   &
                'daylight fraction of radiation timestep',   &
                'percent', missing_value=missing_value)

!-----------------------------------------------------------------------


end subroutine diag_field_init



!######################################################################
! <SUBROUTINE NAME="obtain_astronomy_variables">
!  <OVERVIEW>
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_astronomy_variables (is, ie, js, je, lat, lon,     &
!                                       Astro, Astro2)  
!
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!     astronomy_type structure; It will
!                    be used to determine the insolation at toa seen
!                    by the shortwave radiation code
!  </INOUT>
!  <INOUT NAME="Astro2" TYPE="astronomy_type">
!     astronomy_type structure, defined when renormal-
!                    ization is active. the same components are defined
!                    as for Astro, but they are valid over the current
!                    physics timestep.
!  </INOUT>
!  <INOUT NAME="Astronomy_inp" TYPE="astronomy_inp_type">
!     astronomy_inp_type structure, optionally used to input astronom-
!     ical forcings, when it is desired to specify them rather than use
!     astronomy_mod. Used in various standalone applications.
!  </INOUT>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
! </SUBROUTINE>
!
subroutine obtain_astronomy_variables (is, ie, js, je, lat, lon,     &
                                       Astro, Astro2, Astronomy_inp)  

!---------------------------------------------------------------------
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!---------------------------------------------------------------------
integer,                     intent(in)    ::  is, ie, js, je
real, dimension(:,:),        intent(in)    ::  lat, lon
type(astronomy_type),        intent(inout) ::  Astro, Astro2
type(astronomy_inp_type),   intent(inout), optional ::  &
                                               Astronomy_inp

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  
!                   [ radians ]
!      lon          longitude of model points 
!                   [ radians ]
!
!   intent(inout) variables:
!
!      Astro         astronomy_type structure; contains the following
!                    components defined in this subroutine that will
!                    be used to determine the insolation at toa seen
!                    by the shortwave radiation code 
!         solar         shortwave flux factor: cosine of zenith angle *
!                       daylight fraction / (earth-sun distance squared)
!                       [ non-dimensional ]
!         cosz          cosine of zenith angle --  mean value over
!                       appropriate averaging interval
!                       [ non-dimensional ]
!         fracday       fraction of timestep during which the sun is 
!                       shining
!                       [ non-dimensional ]
!         rrsun         inverse of square of earth-sun distance, 
!                       relative to the mean square of earth-sun 
!                       distance
!                       [ non-dimensional ]
!
!      Astro2        astronomy_type structure, defined when renormal-
!                    ization is active. the same components are defined
!                    as for Astro, but they are valid over the current
!                    physics timestep.
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      type(time_type)                   :: Dt_zen, Dt_zen2
      type(time_type)                   :: Rad1    
      real, dimension(ie-is+1, je-js+1) ::                            &
                                           cosz_r, solar_r, fracday_r, &
                                           cosz_p, solar_p, fracday_p, &
                                           cosz_a, fracday_a
      real                              :: rrsun_r, rrsun_p, rrsun_a
      integer                           :: nz
      

!--------------------------------------------------------------------
!  local variables:
!
!     Dt_zen        time-type variable containing the components of the
!                   radiation time step, needed unless do_average is
!                   true or this is not a radiation step and renormal-
!                   ize_sw_fluxes is true
!     Dt_zen2       time-type variable containing the components of the
!                   physics time step, needed when renormalize_sw_fluxes
!                   or do_average is true
!     cosz_r        cosine of zenith angle --  mean value over
!                   radiation time step            
!                   [ non-dimensional ]
!     solar_r       shortwave flux factor relevant over radiation time
!                   step: cosine of zenith angle * daylight fraction / 
!                   (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_r     fraction of timestep during which the sun is 
!                   shining over radiation time step
!                   [ non-dimensional ]
!     cosz_p        cosine of zenith angle --  mean value over
!                   physics time step            
!                   [ non-dimensional ]
!     solar_p       shortwave flux factor relevant over physics time
!                   step: cosine of zenith angle * daylight fraction / 
!                   (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_p     fraction of timestep during which the sun is 
!                   shining over physics time step
!                   [ non-dimensional ]
!     cosz_a        cosine of zenith angle --  mean value over
!                   next radiation time step            
!                   [ non-dimensional ]
!     solar_a       shortwave flux factor relevant over next radiation 
!                   time step: cosine of zenith angle * daylight 
!                   fraction / (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_a     fraction of timestep during which the sun is 
!                   shining over next radiation time step
!                   [ non-dimensional ]
!     rrsun_r       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over radiation time step
!                   [ non-dimensional ]
!     rrsun_p       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over physics time step
!                   [ non-dimensional ]
!     rrsun_a       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over next radiation time step
!                   [ non-dimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    allocate the components of the astronomy_type structure which will
!    return the astronomical inputs to radiation (cosine of zenith 
!    angle, daylight fraction, solar flux factor and earth-sun distance)
!    that are to be used on the current step.
!---------------------------------------------------------------------
      allocate ( Astro%cosz   (size(lat,1), size(lat,2) ) )
      allocate ( Astro%fracday(size(lat,1), size(lat,2) ) )
      allocate ( Astro%solar  (size(lat,1), size(lat,2) ) )
      allocate ( Astro%cosz_p   (size(lat,1), size(lat,2),   &
                                                 Rad_control%nzens)  )
      allocate ( Astro%fracday_p(size(lat,1), size(lat,2),   &
                                                 Rad_control%nzens)  )
      allocate ( Astro%solar_p  (size(lat,1), size(lat,2),   &
                                                 Rad_control%nzens)  )

!---------------------------------------------------------------------
!    case 0: input parameters.
!---------------------------------------------------------------------
      if (present (Astronomy_inp)) then
        Astro%rrsun = Astronomy_inp%rrsun
        Astro%fracday(:,:) = Astronomy_inp%fracday(is:ie,js:je)
        Astro%cosz (:,:) = cos(   &
                        Astronomy_inp%zenith_angle(is:ie,js:je)/RADIAN)
        Astro%solar(:,:) = Astro%cosz(:,:)*Astro%fracday(:,:)* &
                           Astro%rrsun
        Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)
        do nz = 1, Rad_control%nzens
          Astro%fracday_p(:,:,nz) = Astro%fracday(:,:)
          Astro%cosz_p(:,:,nz) = Astro%cosz(:,:)
          Astro%solar_p(:,:,nz) = Astro%solar(:,:)
        end do

!---------------------------------------------------------------------
!    case 1: diurnally-varying shortwave radiation.
!---------------------------------------------------------------------
      else if (Sw_control%do_diurnal) then

!-------------------------------------------------------------------
!    convert the radiation timestep and the model physics timestep
!    to time_type variables.
!-------------------------------------------------------------------
        Dt_zen  = set_time (sw_rad_time_step, 0)
        Dt_zen2 = set_time (dt, 0)
        
!---------------------------------------------------------------------
!    calculate the astronomical factors averaged over the radiation time
!    step between Rad_time and Rad_time + Dt_zen. these values are 
!    needed on radiation steps. output is stored in Astro_rad.
!---------------------------------------------------------------------
        if (do_sw_rad) then
          if (Rad_control%hires_coszen) then
            Rad1 = Rad_time
            do nz=1,Rad_control%nzens
              call diurnal_solar (lat, lon, Rad1, cosz_r,  &
                                   fracday_r, rrsun_r, dt_time=Dt_zen2)
              fracday_r = MIN (fracday_r, 1.00)
              solar_r = cosz_r*fracday_r*rrsun_r
              Astro%cosz_p(:,:,nz)    = cosz_r
              Astro%fracday_p(:,:,nz) = fracday_r
              Astro%solar_p(:,:,nz)   = solar_r
              Rad1 = Rad1 + Dt_zen2
            end do
          endif
!  calculation for full radiation step:
          call diurnal_solar (lat, lon, Rad_time, cosz_r, fracday_r, &
                              rrsun_r, dt_time=Dt_zen)
          fracday_r = MIN (fracday_r, 1.00)
          solar_r = cosz_r*fracday_r*rrsun_r
        endif

!---------------------------------------------------------------------
!    calculate the astronomical factors averaged over the physics time
!    step between Rad_time and Rad_time + Dt_zen2. these values are
!    needed if either renormalization or time-averaging is active. store
!    the astronomical outputs in Astro_phys.
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then
          call diurnal_solar (lat, lon, Rad_time, cosz_p, fracday_p, &
                              rrsun_p, dt_time=Dt_zen2)
          fracday_p = MIN (fracday_p, 1.00)
          solar_p = cosz_p*fracday_p*rrsun_p
        endif

!--------------------------------------------------------------------
!    define the astronomy_type variable(s) to be returned and used in 
!    the radiation calculation. Astro contains the values to be used
!    in the radiation calculation, Astro2 contains values relevant 
!    over the current physics timestep and is used for renormalization.
!    when renormalization is active, the physics step set is always 
!    needed, and in addition on radiation steps, the radiation step
!    values are needed. 
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then
          if (.not. do_sw_rad) then
            Astro%cosz    = cosz_p
            Astro%fracday = fracday_p
            Astro%solar   = solar_p
            Astro%rrsun   = rrsun_p
          else 
            Astro%cosz    = cosz_r
            Astro%fracday = fracday_r
            Astro%solar   = solar_r
            Astro%rrsun   = rrsun_r
            allocate ( Astro2%fracday(size(lat,1), size(lat,2) ) )
            allocate ( Astro2%cosz   (size(lat,1), size(lat,2) ) )
            allocate ( Astro2%solar  (size(lat,1), size(lat,2) ) )
            Astro2%cosz    = cosz_p
            Astro2%fracday = fracday_p
            Astro2%solar   = solar_p
            Astro2%rrsun   = rrsun_p
          endif

!---------------------------------------------------------------------
!    if renormalization is active, then only the values applicable over
!    radiation steps are needed. 
!---------------------------------------------------------------------
        else                 
          Astro%cosz    = cosz_r
          Astro%fracday = fracday_r
          Astro%solar   = solar_r
          Astro%rrsun   = rrsun_r
        endif

!---------------------------------------------------------------------
!    when in the gcm and on a radiation calculation step, define cosine
!    of zenith angle valid over the next radiation step. this is needed 
!    so that the ocean albedo (function of zenith angle) may be properly
!    defined and provided as input to the radiation package on the next
!    timestep.
!----------------------------------------------------------------------
        if (do_sw_rad) then
          call diurnal_solar (lat, lon, Rad_time+Dt_zen, cosz_a,   &
                              fracday_a, rrsun_a, dt_time=Dt_zen)
          Rad_output%coszen_angle(is:ie,js:je) = cosz_a(:,:)
        endif  ! (do_sw_rad)

!---------------------------------------------------------------------
!    case 2: annual-mean shortwave radiation.
!---------------------------------------------------------------------
      else if (Sw_control%do_annual) then
        call annual_mean_solar (js, je, lat, Astro%cosz, Astro%solar,&
                                Astro%fracday, Astro%rrsun)

!---------------------------------------------------------------------
!    save the cosine of zenith angle on the current step to be used to 
!    calculate ocean albedo for use on the next radiation timestep.
!---------------------------------------------------------------------
        Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)

!---------------------------------------------------------------------
!    case 3: daily-mean shortwave radiation.
!---------------------------------------------------------------------
      else if (Sw_control%do_daily_mean) then
        call daily_mean_solar (lat, Rad_time, Astro%cosz,  &
                               Astro%fracday, Astro%rrsun)
        Astro%solar = Astro%cosz*Astro%rrsun*Astro%fracday

!---------------------------------------------------------------------
!    save the cosine of zenith angle on the current step to be used to 
!    calculate ocean albedo for use on the next radiation timestep.
!---------------------------------------------------------------------
        Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)

!----------------------------------------------------------------------
!    if none of the above options are active, write an error message and
!    stop execution.
!----------------------------------------------------------------------
      else
        call error_mesg('radiation_driver_mod', &
             ' no valid zenith angle specification', FATAL)
      endif

!-------------------------------------------------------------------


end subroutine obtain_astronomy_variables 



!####################################################################
! <SUBROUTINE NAME="radiation_calc">
!  <OVERVIEW>
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_calc (is, ie, js, je, Rad_time, Time_diag,  &
!                           lat, lon, Atmos_input, Surface, Rad_gases, &
!                           Aerosol_props, Aerosol, r, Cldrad_props,   &
!                           Cld_spec, Astro, Rad_output, Lw_output,   &
!                           Sw_output, Fsrad_output, mask, kbot)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <IN NAME="Rad_time" TYPE="time_type">
!      Rad_time          time at which the radiative fluxes are to apply
!                        [ time_type (days, seconds) ] 
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data to radiation package
!  </IN>
!  <INOUT NAME="Aerosol_props" TYPE="aerosol_properties_type">
!   Aerosol radiative properties
!  </INOUT>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   Cloud radiative properties
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud microphysical and physical parameters to radiation package, 
!                     contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical input data for the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties to radiation package, , contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!  </IN>
!  <INOUT NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Fsrad_output" TYPE="Fsrad_output_type">
!   radiation output data from the original_fms_rad
!                        radiation package, when that package 
!                        is active
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
subroutine radiation_calc (is, ie, js, je, Rad_time, Time_diag,  &
                           lat, lon, Atmos_input, Surface, Rad_gases, &
                           Aerosol_props, Aerosol, r, Cldrad_props,   &
                           Cld_spec, Astro, Rad_output, Lw_output,   &
                           Sw_output, Fsrad_output, Aerosol_diags, &
                           mask, kbot)

!--------------------------------------------------------------------
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
type(time_type),              intent(in)             :: Rad_time,   &
                                                        Time_diag
real, dimension(:,:),         intent(in)             :: lat, lon
type(atmos_input_type),       intent(in)             :: Atmos_input 
type(surface_type),           intent(in)             :: Surface
type(radiative_gases_type),   intent(inout)          :: Rad_gases
type(aerosol_type),           intent(in)             :: Aerosol
real, dimension(:,:,:,:),     intent(inout)          :: r   
type(aerosol_properties_type),intent(inout)          :: Aerosol_props
type(cldrad_properties_type), intent(in)             :: Cldrad_props
type(cld_specification_type), intent(in)             :: Cld_spec
type(astronomy_type),         intent(in)             :: Astro
type(rad_output_type),        intent(inout)          :: Rad_output
type(lw_output_type), dimension(:),  intent(inout)   :: Lw_output
type(sw_output_type), dimension(:), intent(inout)     :: Sw_output
type(fsrad_output_type),      intent(inout)          :: Fsrad_output
type(aerosol_diagnostics_type), intent(inout)        :: Aerosol_diags
real, dimension(:,:,:),       intent(in),   optional :: mask
integer, dimension(:,:),      intent(in),   optional :: kbot

!-----------------------------------------------------------------------
!    intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Rad_time          time at which the radiative fluxes are to apply
!                        [ time_type (days, seconds) ] 
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]  
!      lat               latitude of model points on model grid 
!                        [ radians ]
!      lon               longitude of model points on model grid 
!                        [ radians ]
!      Atmos_input       atmospheric input data for the radiation 
!                        package 
!                        [ atmos_input_type ]
!      Surface           surface input data to the radiation package
!                        [ surface_type ]
!      Rad_gases         radiative gas input data for the radiation 
!                        package
!                        [ radiative_gases_type ]
!      Aerosol_props     aerosol radiative property input data for the 
!                        radiation package
!                        [ aerosol_properties_type ]
!      Aerosol           aerosol input data to the radiation package
!                        [ aerosol_type ]
!      Cldrad_props      cloud radiative property input data for the 
!                        radiation package 
!                        [ cldrad_properties_type ]
!      Cld_spec          cloud specification input data for the 
!                        radiation package
!                        [ cld_specification_type ]
!      Astro             astronomical input data for the radiation 
!                        package 
!                        [ astronomy_type ]
!      Aerosol_diags     aerosol diagnostic output                  
!                        [ aerosol_diagnostics_type ]
!
!
!    intent(out) variables:
!
!      Rad_output        radiation output data needed by other modules
!                        [ rad_output_type ]
!      Lw_output         longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!                        [ lw_output_type ]
!          The following are the components of Lw_output:
!                 flxnet    net longwave flux at model flux levels 
!                           (including the ground and the top of the 
!                           atmosphere).
!                 heatra    longwave heating rates in model layers.
!                 flxnetcf  net longwave flux at model flux levels 
!                           (including the ground and the top of the 
!                           atmosphere) computed for cloud-free case.
!                 heatra    longwave heating rates in model layers 
!                           computed for cloud-free case.
!      Sw_output         shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!                        [ sw_output_type ]
!      Fsrad_output      radiation output data from the original_fms_rad
!                        radiation package, when that package 
!                        is active
!                        [ fsrad_output_type ]
!
!    intent(in), optional variables:
!
!      mask              present when running eta vertical coordinate,
!                        mask to define values at points below ground   
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!
!----------------------------------------------------------------------

      integer :: kmax, nz

!---------------------------------------------------------------------
!    all_column_radiation and all_level_radiation are included as 
!    future controls which may be utiliized to execute the radiation
!    code on a grid other than the model grid. in the current release
!    however, both must be .true.. 
!---------------------------------------------------------------------
      if (all_column_radiation .and. all_level_radiation) then

!--------------------------------------------------------------------
!    call routines to perform radiation calculations, either using the
!    sea_esf_rad or original_fms_rad radiation package.
!---------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call sea_esf_rad (is, ie, js, je, Rad_time, Atmos_input, &
                            Surface, Astro, Rad_gases, Aerosol,  &
                            Aerosol_props, Cldrad_props, Cld_spec,  &
                            Lw_output, Sw_output, Aerosol_diags, r)

!--------------------------------------------------------------------
!    define tropopause fluxes for diagnostic use later.
!--------------------------------------------------------------------
        if (do_lw_rad .or. do_sw_rad) then
          call flux_trop_calc (is, ie, js, je, lat,  &
                               Atmos_input, Lw_output(1), Sw_output(1) )
        endif

        else  
          call original_fms_rad (is, ie, js, je, Atmos_input%phalf,  &
                                 lat, lon, do_clear_sky_pass,   &
                                 Rad_time, Time_diag, Atmos_input,  &
                                 Surface, Astro, Rad_gases,   &
                                 Cldrad_props, Cld_spec,    &
                                 Fsrad_output, mask=mask, kbot=kbot) 
        endif
      else  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    when this option is coded, replace this error_mesg code with
!    code which will map the input fields from the model grid to
!    the desired radiation grid. A preliminary version of code to per-
!    form this task (at least some of it) is found with the inchon
!    tagged version of this module. it is removed here, since it has
!    not been tested or validated and is considered undesirable in
!    a code being prepared for public release. no immediate need for
!    it is seen at this time, but it will be added back when such need
!    arises.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call error_mesg ('radiation_driver_mod', &
               ' ability to calculate radiation on subset of columns'//&
              ' and/or levels not yet implemented',  FATAL)

      endif  ! (all_column .and. all_level)

!---------------------------------------------------------------------
!    define the components of Rad_output to be passed back to 
!    radiation_driver --  total and shortwave radiative heating rates 
!    for standard and clear-sky case (if desired), and surface long- 
!    and short-wave fluxes.  mask out any below ground values if 
!    necessary.
!---------------------------------------------------------------------
        if (do_sea_esf_rad) then
          if (do_sw_rad) then
            Rad_output%tdtsw(is:ie,js:je,:,:) =    &
                              Sw_output(1)%hsw(:,:,:,:)/SECONDS_PER_DAY
            Rad_output%ufsw(is:ie,js:je,:,:) =    &
                              Sw_output(1)%ufsw(:,:,:,:)
            Rad_output%dfsw(is:ie,js:je,:,:) =    &
                              Sw_output(1)%dfsw(:,:,:,:)
          endif
          if (present(mask)) then
            if (do_lw_rad) then
              Rad_output%tdtlw(is:ie,js:je,:) =   &
                        (Lw_output(1)%heatra(:,:,:)/SECONDS_PER_DAY)*  &
                                                            mask(:,:,:)
              Rad_output%flxnet(is:ie,js:je,:) =  &
                         Lw_output(1)%flxnet(:,:,:)*mask(:,:,:)
            endif
            do nz = 1, Rad_control%nzens
              Rad_output%tdt_rad (is:ie,js:je,:,nz) =   &
                         (Rad_output%tdtsw(is:ie,js:je,:,nz) + &
                            Rad_output%tdtlw(is:ie,js:je,:))*mask(:,:,:)
            end do
          else
            if (do_lw_rad) then
               Rad_output%tdtlw(is:ie,js:je,:) =   &
                             Lw_output(1)%heatra(:,:,:)/SECONDS_PER_DAY
               Rad_output%flxnet(is:ie,js:je,:) =  &
                          Lw_output(1)%flxnet(:,:,:)
            endif
            do nz = 1, Rad_control%nzens
              Rad_output%tdt_rad (is:ie,js:je,:,nz) =  &
                             (Rad_output%tdtsw(is:ie,js:je,:,nz) +   &
                                       Rad_output%tdtlw(is:ie,js:je,:))
            end do
          endif
          if (do_clear_sky_pass) then
            do nz = 1, Rad_control%nzens
              if (do_sw_rad) then
                Rad_output%tdtsw_clr(is:ie,js:je,:,nz) =   &
                           Sw_output(1)%hswcf(:,:,:,nz)/SECONDS_PER_DAY
                Rad_output%ufsw_clr(is:ie,js:je,:,nz) =   &
                           Sw_output(1)%ufswcf(:,:,:,nz)
                Rad_output%dfsw_clr(is:ie,js:je,:,nz) =   &
                           Sw_output(1)%dfswcf(:,:,:,nz)
                Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je,nz) =&
                               Sw_output(1)%dfsw_dir_sfc_clr(:,:,nz)
                Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je,nz) =&
                                Sw_output(1)%dfsw_dif_sfc_clr(:,:,nz)
                Rad_output%flux_sw_down_vis_clr(is:ie,js:je,nz) =   &
                                 Sw_output(1)%dfsw_vis_sfc_clr(:,:,nz)
              endif
              if (do_lw_rad) then
                Rad_output%tdtlw_clr(is:ie,js:je,:) =   &
                        Lw_output(1)%heatracf(:,:,:)/SECONDS_PER_DAY
                Rad_output%flxnetcf(is:ie,js:je,:) =  &
                              Lw_output(1)%flxnet(:,:,:)
              endif
              if (present(mask)) then
                Rad_output%tdt_rad_clr(is:ie,js:je,:,nz) =    &
                    (Rad_output%tdtsw_clr(is:ie,js:je,:,nz) +  &
                       Rad_output%tdtlw_clr(is:ie,js:je,:))*mask(:,:,:)
              else
                Rad_output%tdt_rad_clr(is:ie,js:je,:,nz) =    &
                           (Rad_output%tdtsw_clr(is:ie,js:je,:,nz) +  &
                                  Rad_output%tdtlw_clr(is:ie,js:je,:))
              endif
            end do
          endif

          kmax = size (Rad_output%tdtsw,3)
          if (do_sw_rad) then
            Rad_output%flux_sw_surf(is:ie,js:je,:) =   &
                            Sw_output(1)%dfsw(:,:,kmax+1,:) - &
                                        Sw_output(1)%ufsw(:,:,kmax+1,:)
            Rad_output%flux_sw_surf_dir(is:ie,js:je,:) =   &
                                    Sw_output(1)%dfsw_dir_sfc(:,:,:)
            Rad_output%flux_sw_surf_refl_dir(is:ie,js:je,:) =   &
                                    Sw_output(1)%ufsw_dir_sfc(:,:,:)
            Rad_output%flux_sw_surf_dif(is:ie,js:je,:) =   &
                                 Sw_output(1)%dfsw_dif_sfc(:,:,:) - &
                                      Sw_output(1)%ufsw_dif_sfc(:,:,:)
            Rad_output%flux_sw_down_vis_dir(is:ie,js:je,:) =   &
                                   Sw_output(1)%dfsw_vis_sfc_dir(:,:,:)
            Rad_output%flux_sw_down_vis_dif(is:ie,js:je,:) =   &
                                   Sw_output(1)%dfsw_vis_sfc_dif(:,:,:)
            Rad_output%flux_sw_down_total_dir(is:ie,js:je,:) =   &
                                       Sw_output(1)%dfsw_dir_sfc(:,:,:)
            Rad_output%flux_sw_down_total_dif(is:ie,js:je,:) =   &
                                      Sw_output(1)%dfsw_dif_sfc(:,:,:)
            Rad_output%flux_sw_vis (is:ie,js:je,:) =   &
                               Sw_output(1)%dfsw_vis_sfc(:,:,:) - &
                                        Sw_output(1)%ufsw_vis_sfc(:,:,:)
            Rad_output%flux_sw_vis_dir (is:ie,js:je,:) =   &
                                   Sw_output(1)%dfsw_vis_sfc_dir(:,:,:)
            Rad_output%flux_sw_refl_vis_dir (is:ie,js:je,:) =   &
                                   Sw_output(1)%ufsw_vis_sfc_dir(:,:,:)
            Rad_output%flux_sw_vis_dif (is:ie,js:je,:) =   &
                            Sw_output(1)%dfsw_vis_sfc_dif(:,:,:) - &
                                   Sw_output(1)%ufsw_vis_sfc_dif(:,:,:)
          endif
          if (do_lw_rad) then
            Rad_output%flux_lw_surf(is:ie,js:je) =    &
                         STEFAN*Atmos_input%temp(:,:,kmax+1)**4 -   &
                                       Lw_output(1)%flxnet(:,:,kmax+1)
          endif
        else
          Rad_output%tdtsw(is:ie,js:je,:,1) = Fsrad_output%tdtsw(:,:,:)
          if (present(mask)) then
            Rad_output%tdt_rad (is:ie,js:je,:,1) =   &
                               (Rad_output%tdtsw(is:ie,js:je,:,1) + &
                                Fsrad_output%tdtlw (:,:,:))*mask(:,:,:)
          else
            Rad_output%tdt_rad (is:ie,js:je,:,1) =   &
                               (Rad_output%tdtsw(is:ie,js:je,:,1) +   &
                                Fsrad_output%tdtlw (:,:,:))
          endif
          if (do_clear_sky_pass) then
            Rad_output%tdtsw_clr(is:ie,js:je,:,1) =    &
                                          Fsrad_output%tdtsw_clr(:,:,:)
            if (present(mask)) then
              Rad_output%tdt_rad_clr(is:ie,js:je,:,1) =    &
                         (Rad_output%tdtsw_clr(is:ie,js:je,:,1) +  &
                          Fsrad_output%tdtlw_clr(:,:,:))*mask(:,:,:)
            else
              Rad_output%tdt_rad_clr(is:ie,js:je,:,1) =    &
                         (Rad_output%tdtsw_clr(is:ie,js:je,:,1) +    &
                          Fsrad_output%tdtlw_clr(:,:,:))
            endif
          endif
          Rad_output%flux_sw_surf(is:ie,js:je,1) =    &
                                             Fsrad_output%swdns(:,:) - &
                                             Fsrad_output%swups(:,:)
          Rad_output%flux_lw_surf(is:ie,js:je) = Fsrad_output%lwdns(:,:)
        endif ! (do_sea_esf_rad)


!---------------------------------------------------------------------




end subroutine radiation_calc




!######################################################################
! <SUBROUTINE NAME="update_rad_fields">
!  <OVERVIEW>
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!  </OVERVIEW>
!  <DESCRIPTION>
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call update_rad_fields (is, ie, js, je, Time_diag, Astro2,   &
!                              Sw_output, Astro, Rad_output, flux_ratio)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!      Time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!   astronomical properties on model grid, usually
!                   valid over radiation timestep on entry, on exit are 
!                   valid over model timestep when renormalizing
!  </INOUT>
!  <IN NAME="Astro2" TYPE="astronomy_type">
!   astronomical properties on model grid, valid over 
!                   physics timestep, used when renormalizing sw fluxes
!  </IN>
!  <INOUT NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </IN>
!  <OUT NAME="flux_ratio" TYPE="real">
!   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!  </OUT>
! </SUBROUTINE>
!
subroutine update_rad_fields (is, ie, js, je, Time_diag, Astro2,   &
                              Sw_output, Astro, Rad_output, flux_ratio)

!---------------------------------------------------------------------
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!--------------------------------------------------------------------

integer,                 intent(in)    ::  is, ie, js, je
type(time_type),         intent(in)    ::  Time_diag
type(astronomy_type),    intent(in)    ::  Astro2
type(sw_output_type), dimension(:),   intent(inout)    ::  Sw_output
type(astronomy_type),    intent(inout) ::  Astro
type(rad_output_type),   intent(inout) ::  Rad_output
real,  dimension(:,:),   intent(out)   ::  flux_ratio

!-------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      Time_diag    time on next timestep, used as stamp for diag-
!                   nostic output  [ time_type  (days, seconds) ]  
!      Astro2       astronomical properties on model grid, valid over 
!                   physics timestep, used when renormalizing sw fluxes
!                   [astronomy_type]
!      Sw_output    shortwave output variables on model grid,
!                   [sw_output_type]     
!
!  intent(inout) variables:
!
!      Astro        astronomical properties on model grid, usually
!                   valid over radiation timestep on entry, on exit are 
!                   valid over model timestep when renormalizing
!                   [astronomy_type]
!      Rad_output   radiation output variables on model grid, valid
!                   on entry over either physics or radiation timestep, 
!                   on exit are valid over physics step when renormal-
!                   izing sw fluxes
!                   [rad_output_type]     
!
!  intent(out) variables:
!
!      flux_ratio   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
      real, dimension (is:ie, js:je,                                  &
                      size(Rad_output%tdt_rad,3))  ::  tdtlw, tdtlw_clr
      integer   :: i, j, k
      integer   :: nz

!---------------------------------------------------------------------
!  local variables:
!
!     tdtlw              longwave heating rate
!                        [ deg K sec(-1) ]
!     tdtlw_clr          longwave heating rate under clear sky 
!                        conditions
!                        [ deg K sec(-1) ]
!     i,j,k              do-loop indices
!
!---------------------------------------------------------------------

      if (renormalize_sw_fluxes) then

!----------------------------------------------------------------------
!    if sw fluxes are to be renormalized, save the heating rates, fluxes
!    and solar factor calculated on radiation steps.
!---------------------------------------------------------------------
        if (do_sw_rad) then
          solar_save(is:ie,js:je)  = Astro%solar(:,:)
          dfsw_save(is:ie,js:je,:,:) = Sw_output(1)%dfsw(:, :,:,:)
          ufsw_save(is:ie,js:je,:,:) = Sw_output(1)%ufsw(:, :,:,:)
          if (do_swaerosol_forcing) then
            dfsw_ad_save(is:ie,js:je,:,:) =   &
                                    Sw_output(indx_swaf)%dfsw(:, :,:,:)
            ufsw_ad_save(is:ie,js:je,:,:) =    &
                                    Sw_output(indx_swaf)%ufsw(:, :,:,:)
          endif
          fsw_save(is:ie,js:je,:,:)  = Sw_output(1)%fsw(:, :,:,:)
          hsw_save(is:ie,js:je,:,:)  = Sw_output(1)%hsw(:, :,:,:)
          flux_sw_surf_save(is:ie,js:je,:) =    &
                                 Rad_output%flux_sw_surf(is:ie,js:je,:)
          flux_sw_surf_dir_save(is:ie,js:je,:) =    &
                             Rad_output%flux_sw_surf_dir(is:ie,js:je,:)
          flux_sw_surf_refl_dir_save(is:ie,js:je,:) =    &
                             Rad_output%flux_sw_surf_refl_dir(is:ie,js:je,:)
          flux_sw_surf_dif_save(is:ie,js:je,:) =    &
                             Rad_output%flux_sw_surf_dif(is:ie,js:je,:)
          flux_sw_down_vis_dir_save(is:ie,js:je,:) =    &
                         Rad_output%flux_sw_down_vis_dir(is:ie,js:je,:)
          flux_sw_down_vis_dif_save(is:ie,js:je,:) =    &
                         Rad_output%flux_sw_down_vis_dif(is:ie,js:je,:)
          flux_sw_down_total_dir_save(is:ie,js:je,:) =    &
                      Rad_output%flux_sw_down_total_dir(is:ie,js:je,:)
          flux_sw_down_total_dif_save(is:ie,js:je,:) =    &
                      Rad_output%flux_sw_down_total_dif(is:ie,js:je,:)
          flux_sw_vis_save(is:ie,js:je,:) =    &
                              Rad_output%flux_sw_vis(is:ie,js:je,:)
          flux_sw_vis_dir_save(is:ie,js:je,:) =    &
                   Rad_output%flux_sw_vis_dir(is:ie,js:je,:)
          flux_sw_refl_vis_dir_save(is:ie,js:je,:) =    &
                   Rad_output%flux_sw_refl_vis_dir(is:ie,js:je,:)
          flux_sw_vis_dif_save(is:ie,js:je,:) =    &
                               Rad_output%flux_sw_vis_dif(is:ie,js:je,:)
          sw_heating_save(is:ie,js:je,:,:) =    &
                              Rad_output%tdtsw(is:ie,js:je,:,:)
          tot_heating_save(is:ie,js:je,:,:) =    &
                              Rad_output%tdt_rad(is:ie,js:je,:,:)
          swdn_special_save(is:ie,js:je,:,:) =   &
                                     Sw_output(1)%swdn_special(:,:,:,:)
          swup_special_save(is:ie,js:je,:,:) =   &
                                     Sw_output(1)%swup_special(:,:,:,:)
          if (do_clear_sky_pass) then
            sw_heating_clr_save(is:ie,js:je,:,:) =    &
                              Rad_output%tdtsw_clr(is:ie,js:je,:,:)
            tot_heating_clr_save(is:ie,js:je,:,:) =    &
                              Rad_output%tdt_rad_clr(is:ie,js:je,:,:)
            dfswcf_save(is:ie,js:je,:,:) = Sw_output(1)%dfswcf(:, :,:,:)
            ufswcf_save(is:ie,js:je,:,:) = Sw_output(1)%ufswcf(:, :,:,:)
            if (do_swaerosol_forcing) then
              dfswcf_ad_save(is:ie,js:je,:,:) =    &
                                  Sw_output(indx_swaf)%dfswcf(:, :,:,:)
              ufswcf_ad_save(is:ie,js:je,:,:) =   &
                                  Sw_output(indx_swaf)%ufswcf(:, :,:,:)
            endif
            fswcf_save(is:ie,js:je,:,:)  = Sw_output(1)%fswcf(:, :,:,:)
            hswcf_save(is:ie,js:je,:,:)  = Sw_output(1)%hswcf(:, :,:,:)
            flux_sw_down_total_dir_clr_save(is:ie,js:je,:) =    &
                   Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je,:)
            flux_sw_down_total_dif_clr_save(is:ie,js:je,:) =    &
                    Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je,:)
            flux_sw_down_vis_clr_save(is:ie,js:je,:) =    &
                    Rad_output%flux_sw_down_vis_clr(is:ie,js:je,:)
            swdn_special_clr_save(is:ie,js:je,:,:) =  &
                                Sw_output(1)%swdn_special_clr(:,:,:,:)
            swup_special_clr_save(is:ie,js:je,:,:) =  &
                                Sw_output(1)%swup_special_clr(:,:,:,:)
          endif 

!---------------------------------------------------------------------
!    define the ratio of the solar factor valid over this physics step
!    to that valid over the current radiation timestep.
!---------------------------------------------------------------------
          do j=1,je-js+1
            do i=1,ie-is+1
              if (solar_save(i+is-1,j+js-1) /= 0.0) then
                flux_ratio(i, j) = Astro2%solar(i,j)/   &
                                            solar_save(i+is-1,j+js-1)
              else
                flux_ratio(i,j) = 0.0
              endif

!---------------------------------------------------------------------
!    move the physics-step values(Astro2) to Astro, which will be used 
!    to calculate diagnostics. the radiation_step values (Astro) are no
!    longer needed.
!---------------------------------------------------------------------
              Astro%cosz(i,j) = Astro2%cosz(i,j)
              Astro%fracday(i,j) = Astro2%fracday(i,j)
              Astro%solar(i,j) = Astro2%solar(i,j)
              Astro%rrsun = Astro2%rrsun
            end do
          end do

!----------------------------------------------------------------------
!    on non-radiation steps define the ratio of the current solar factor
!    valid for this physics step to that valid for the last radiation 
!    step. 
!----------------------------------------------------------------------
        else 
          do j=1,je-js+1
            do i=1,ie-is+1
              if (solar_save(i+is-1,j+js-1) /= 0.0) then
                flux_ratio(i, j) = Astro%solar(i,j)/   &
                                            solar_save(i+is-1,j+js-1)
              else
                flux_ratio(i,j) = 0.0
              endif
            end do
          end do
        endif  ! (do_sw_rad)

!---------------------------------------------------------------------
!    redefine the total and shortwave heating rates, along with surface
!    sw fluxes, as a result of the difference in solar factor (the 
!    relative earth-sun motion) between the current physics and current
!    radiation timesteps.
!---------------------------------------------------------------------
        nz = current_sw_zenith_step
        tdtlw(:,:,:) = tot_heating_save(is:ie,js:je,:,nz) -    &
                       sw_heating_save(is:ie,js:je,:,nz)
        do k=1, size(Rad_output%tdt_rad,3)
          Rad_output%tdtsw(is:ie,js:je,k,nz) =    &
                       sw_heating_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
        end do
        do k=1, size(Rad_output%tdt_rad,3)+1
          Rad_output%ufsw(is:ie,js:je,k,nz) =    &
                       ufsw_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          Rad_output%dfsw(is:ie,js:je,k,nz) =    &
                       dfsw_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
        end do
        Rad_output%tdt_rad(is:ie,js:je,:,nz) = tdtlw(:,:,:) +    &
                                     Rad_output%tdtsw(is:ie,js:je,:,nz)
        Rad_output%flux_sw_surf(is:ie,js:je,nz) = flux_ratio(:,:)*    &
                                      flux_sw_surf_save(is:ie,js:je,nz)
        Rad_output%flux_sw_surf_dir(is:ie,js:je,nz) = flux_ratio(:,:)* &
                                  flux_sw_surf_dir_save(is:ie,js:je,nz)
        Rad_output%flux_sw_surf_refl_dir(is:ie,js:je,nz) =flux_ratio(:,:)*&
                               flux_sw_surf_refl_dir_save(is:ie,js:je,nz)
        Rad_output%flux_sw_surf_dif(is:ie,js:je,nz) = flux_ratio(:,:)* &
                           flux_sw_surf_dif_save(is:ie,js:je,nz)
        Rad_output%flux_sw_down_vis_dir(is:ie,js:je,nz) =    &
               flux_ratio(:,:)*flux_sw_down_vis_dir_save(is:ie,js:je,nz)
        Rad_output%flux_sw_down_vis_dif(is:ie,js:je,nz) =  &
              flux_ratio(:,:)*flux_sw_down_vis_dif_save(is:ie,js:je,nz)
        Rad_output%flux_sw_down_total_dir(is:ie,js:je,nz) =   &
             flux_ratio(:,:)*flux_sw_down_total_dir_save(is:ie,js:je,nz)
        Rad_output%flux_sw_down_total_dif(is:ie,js:je,nz) =   &
             flux_ratio(:,:)*flux_sw_down_total_dif_save(is:ie,js:je,nz)
        Rad_output%flux_sw_vis(is:ie,js:je,nz) = flux_ratio(:,:)*    &
                                    flux_sw_vis_save(is:ie,js:je,nz)
        Rad_output%flux_sw_vis_dir(is:ie,js:je,nz) = flux_ratio(:,:)*  &
                                   flux_sw_vis_dir_save(is:ie,js:je,nz)
        Rad_output%flux_sw_refl_vis_dir(is:ie,js:je,nz) = flux_ratio(:,:)*&
                                flux_sw_refl_vis_dir_save(is:ie,js:je,nz)
        Rad_output%flux_sw_vis_dif(is:ie,js:je,nz) = flux_ratio(:,:)* &
                                   flux_sw_vis_dif_save(is:ie,js:je,nz)
        if (do_clear_sky_pass) then
          tdtlw_clr(:,:,:) = tot_heating_clr_save(is:ie,js:je,:,nz) -  &
                             sw_heating_clr_save (is:ie,js:je,:,nz)
          Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je,nz) =   &
                  flux_ratio(:,:)*     &
                        flux_sw_down_total_dir_clr_save(is:ie,js:je,nz)
          Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je,nz) =   &
                  flux_ratio(:,:) * &
                        flux_sw_down_total_dif_clr_save(is:ie,js:je,nz)
          Rad_output%flux_sw_down_vis_clr(is:ie,js:je,nz) =   &
               flux_ratio(:,:)*flux_sw_down_vis_clr_save(is:ie,js:je,nz)
          do k=1, size(Rad_output%tdt_rad,3)
            Rad_output%tdtsw_clr(is:ie,js:je,k,nz) =   &
                        sw_heating_clr_save (is:ie,js:je,k,nz)*   &
                                                        flux_ratio(:,:)
          end do
          do k=1, size(Rad_output%tdt_rad,3)+1
            Rad_output%ufsw_clr(is:ie,js:je,k,nz) =    &
                        ufswcf_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
            Rad_output%dfsw_clr(is:ie,js:je,k,nz) =    &
                        dfswcf_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          end do
          Rad_output%tdt_rad_clr(is:ie,js:je,:,nz) = tdtlw_clr(:,:,:) +&
                           Rad_output%tdtsw_clr(is:ie,js:je,:,nz)
        endif
      else if (all_step_diagnostics) then

!----------------------------------------------------------------------
!    if sw fluxes are to be output on every physics step, save the 
!    heating rates and fluxes calculated on radiation steps.
!---------------------------------------------------------------------
        if (do_sw_rad) then
          if (do_swaerosol_forcing) then
            dfsw_ad_save(is:ie,js:je,:,:) =    &
                                   Sw_output(indx_swaf)%dfsw(:, :,:,:)
            ufsw_ad_save(is:ie,js:je,:,:) =    &
                                   Sw_output(indx_swaf)%ufsw(:, :,:,:)
          endif
          dfsw_save(is:ie,js:je,:,:) = Sw_output(1)%dfsw(:, :,:,:)
          ufsw_save(is:ie,js:je,:,:) = Sw_output(1)%ufsw(:, :,:,:)
          fsw_save(is:ie,js:je,:,:)  = Sw_output(1)%fsw(:, :,:,:)
          hsw_save(is:ie,js:je,:,:)  = Sw_output(1)%hsw(:, :,:,:)
          flux_sw_surf_save(is:ie,js:je,:) =    &
                              Rad_output%flux_sw_surf(is:ie,js:je,:)
          flux_sw_surf_dir_save(is:ie,js:je,:) =    &
                              Rad_output%flux_sw_surf_dir(is:ie,js:je,:)
          flux_sw_surf_refl_dir_save(is:ie,js:je,:) =    &
                           Rad_output%flux_sw_surf_refl_dir(is:ie,js:je,:)
          flux_sw_surf_dif_save(is:ie,js:je,:) =    &
                              Rad_output%flux_sw_surf_dif(is:ie,js:je,:)
          flux_sw_down_vis_dir_save(is:ie,js:je,:) =    &
                          Rad_output%flux_sw_down_vis_dir(is:ie,js:je,:)
          flux_sw_refl_vis_dir_save(is:ie,js:je,:) =    &
                            Rad_output%flux_sw_refl_vis_dir(is:ie,js:je,:)
          flux_sw_down_vis_dif_save(is:ie,js:je,:) =    &
                         Rad_output%flux_sw_down_vis_dif(is:ie,js:je,:)
          flux_sw_down_total_dir_save(is:ie,js:je,:) =    &
                       Rad_output%flux_sw_down_total_dir(is:ie,js:je,:)
          flux_sw_down_total_dif_save(is:ie,js:je,:) =    &
                       Rad_output%flux_sw_down_total_dif(is:ie,js:je,:)
          flux_sw_vis_save(is:ie,js:je,:) =    &
                               Rad_output%flux_sw_vis(is:ie,js:je,:)
          flux_sw_vis_dir_save(is:ie,js:je,:) =    &
                            Rad_output%flux_sw_vis_dir(is:ie,js:je,:)
          flux_sw_vis_dif_save(is:ie,js:je,:) =    &
                               Rad_output%flux_sw_vis_dif(is:ie,js:je,:)
          sw_heating_save(is:ie,js:je,:,:) =    &
                             Rad_output%tdtsw(is:ie,js:je,:,:)
          tot_heating_save(is:ie,js:je,:,:) =    &
                               Rad_output%tdt_rad(is:ie,js:je,:,:)
          swdn_special_save(is:ie,js:je,:,:) =   &
                                 Sw_output(1)%swdn_special(:,:,:,:)
          swup_special_save(is:ie,js:je,:,:) =   &
                                 Sw_output(1)%swup_special(:,:,:,:)
          if (do_clear_sky_pass) then
            sw_heating_clr_save(is:ie,js:je,:,:) =    &
                        Rad_output%tdtsw_clr(is:ie,js:je,:,:)
            tot_heating_clr_save(is:ie,js:je,:,:) =    &
                          Rad_output%tdt_rad_clr(is:ie,js:je,:,:)
            dfswcf_save(is:ie,js:je,:,:) = Sw_output(1)%dfswcf(:, :,:,:)
            ufswcf_save(is:ie,js:je,:,:) = Sw_output(1)%ufswcf(:, :,:,:)
            if (do_swaerosol_forcing) then
              dfswcf_ad_save(is:ie,js:je,:,:) =  &
                                 Sw_output(indx_swaf)%dfswcf(:, :,:,:)
              ufswcf_ad_save(is:ie,js:je,:,:) =   &
                                 Sw_output(indx_swaf)%ufswcf(:, :,:,:)
            endif
            fswcf_save(is:ie,js:je,:,:)  = Sw_output(1)%fswcf(:, :,:,:)
            hswcf_save(is:ie,js:je,:,:)  = Sw_output(1)%hswcf(:, :,:,:)
            flux_sw_down_total_dir_clr_save(is:ie,js:je,:) =    &
                  Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je,:)
            flux_sw_down_total_dif_clr_save(is:ie,js:je,:) =    &
                   Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je,:)
            flux_sw_down_vis_clr_save(is:ie,js:je,:) =    &
                         Rad_output%flux_sw_down_vis_clr(is:ie,js:je,:)
            swdn_special_clr_save(is:ie,js:je,:,:) =   &
                                Sw_output(1)%swdn_special_clr(:,:,:,:)
            swup_special_clr_save(is:ie,js:je,:,:) =   &
                                Sw_output(1)%swup_special_clr(:,:,:,:)
          endif
        endif
      else
        flux_ratio(:,:) = 1.0
      endif  ! (renormalize_sw_fluxes)

!--------------------------------------------------------------------



end subroutine update_rad_fields 



!####################################################################

! <SUBROUTINE NAME="flux_trop_calc">
!  <OVERVIEW>
!    flux_trop_calc defines the shortwave and longwave fluxes at the
!    tropopause immediately after the computation of fluxes at model
!    levels by the radiation algorithms (invoked by radiation_calc).
!  </OVERVIEW>
!  <DESCRIPTION>
!    flux_trop_calc defines the shortwave and longwave fluxes at the
!    tropopause immediately after the computation of fluxes at model
!    levels by the radiation algorithms (invoked by radiation_calc).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_trop_calc          (is, ie, js, je, lat,      &
!                                 Atmos_input,              &
!                                 Lw_output, Sw_output )
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    mean latitude (in radians) of all grid boxes processed by this
!    call to flux_trop_calc   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </IN>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </INOUT>
! </SUBROUTINE>

subroutine flux_trop_calc    (is, ie, js, je, lat, Atmos_input, &
                              Lw_output, Sw_output )

integer,                 intent(in)             :: is, ie, js, je
real,dimension(:,:),     intent(in)             :: lat
type(atmos_input_type),  intent(in)             :: Atmos_input
type(lw_output_type),    intent(inout)          :: Lw_output
type(sw_output_type),    intent(inout)          :: Sw_output

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      Atmos_input  component pflux (pressure at layer boundaries [ Pa ]
!                   is used
!
!   intent(inout) variables:
!      Lw_output    lw_output_type variable containing output from 
!                   the longwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Sw_output    sw_output_type variable containing output from 
!                   the shortwave radiation code of the
!                   sea_esf_rad package, on the model grid

      real, dimension (ie-is+1,je-js+1) ::      lat_deg, tropo_ht
      real, dimension (ie-is+1,je-js+1) ::      netlw_trop,  &
                                                swdn_trop, swup_trop, &
                                                netlw_trop_clr, &
                                                swdn_trop_clr,  &
                                                swup_trop_clr
      integer           :: j, k, nz
      integer           :: ki, i
      integer           :: kmax
      real              :: wtlo, wthi

      kmax = size(Atmos_input%pflux,3) - 1

!---------------------------------------------------------------------
!    compute net downward flux at 1 Pa (top of dynmamical model)
!    here dynamical pressure top is hard-wired to 1 Pa.
!---------------------------------------------------------------------
 
      do j = 1,je-js+1
        do i = 1,ie-is+1
          wtlo = (1.0           - Atmos_input%pflux(i,j,1))/ &
               (Atmos_input%pflux(i,j,2) - Atmos_input%pflux(i,j,1))
          wthi = 1.0 - wtlo
          if (Rad_control%do_lw_rad) then
            netlw_trop(i,j) = wthi*Lw_output%flxnet(i,j,1) + &
                              wtlo*Lw_output%flxnet(i,j,2)
            Lw_output%netlw_special(i,j,4) = netlw_trop(i,j)
            if (do_clear_sky_pass) then
              netlw_trop_clr(i,j) = wthi*Lw_output%flxnetcf(i,j,1) + &
                                    wtlo*Lw_output%flxnetcf(i,j,2)
              Lw_output%netlw_special_clr(i,j,4) = netlw_trop_clr(i,j)
            endif
          endif
          if (Rad_control%do_sw_rad) then
            do nz = 1,Rad_control%nzens
              swdn_trop(i,j) = wthi*Sw_output%dfsw(i,j,1,nz) + &
                               wtlo*Sw_output%dfsw(i,j,2,nz)
              swup_trop(i,j) = wthi*Sw_output%ufsw(i,j,1,nz) + &
                               wtlo*Sw_output%ufsw(i,j,2,nz)
              Sw_output%swdn_special(i,j,4,nz) = swdn_trop(i,j)
              Sw_output%swup_special(i,j,4,nz) = swup_trop(i,j)
              if (do_clear_sky_pass) then
                swdn_trop_clr(i,j) = wthi*Sw_output%dfswcf(i,j,1,nz) +&
                                     wtlo*Sw_output%dfswcf(i,j,2,nz)
                swup_trop_clr(i,j) = wthi*Sw_output%ufswcf(i,j,1,nz) +&
                                     wtlo*Sw_output%ufswcf(i,j,2,nz)
                Sw_output%swdn_special_clr(i,j,4,nz) =   &
                                                  swdn_trop_clr(i,j)
                Sw_output%swup_special_clr(i,j,4,nz) =    &
                                                swup_trop_clr(i,j)
              endif
            end do
          endif
        enddo
      enddo


      if (constant_tropo) then
        tropo_ht(:,:) = trop_ht_constant
! interpolate the fluxes between the appropriate pressures bracketing
! (trop) 

        do j = 1,je-js+1
          do i = 1,ie-is+1
            do k = kmax+1,2,-1
              if (Atmos_input%pflux(i,j,k) >= tropo_ht(i,j) .and.     &
                  Atmos_input%pflux(i,j,k-1) < tropo_ht(i,j))      then
                ki = k
!   the indices for high,low pressure bracketing "tropo_ht" are ki, ki-1
            wtlo = (tropo_ht(i,j) - Atmos_input%pflux(i,j,ki-1))/ &
              (Atmos_input%pflux(i,j,ki) - Atmos_input%pflux(i,j,ki-1))
            wthi = 1.0 - wtlo
            if (Rad_control%do_lw_rad) then
              netlw_trop(i,j) = wtlo*Lw_output%flxnet(i,j,ki) + &
                              wthi*Lw_output%flxnet(i,j,ki-1)
              Lw_output%netlw_special(i,j,1) = netlw_trop(i,j)
              if (do_clear_sky_pass) then
               netlw_trop_clr(i,j) = wtlo*Lw_output%flxnetcf(i,j,ki) + &
                                    wthi*Lw_output%flxnetcf(i,j,ki-1)
               Lw_output%netlw_special_clr(i,j,1) = netlw_trop_clr(i,j)
              endif
            endif
            if (Rad_control%do_sw_rad) then
              do nz = 1,Rad_control%nzens
              swdn_trop(i,j) = wtlo*Sw_output%dfsw(i,j,ki,nz) + &
                               wthi*Sw_output%dfsw(i,j,ki-1,nz)
              swup_trop(i,j) = wtlo*Sw_output%ufsw(i,j,ki,nz) + &
                               wthi*Sw_output%ufsw(i,j,ki-1,nz)
              Sw_output%swdn_special(i,j,1,nz) = swdn_trop(i,j)
              Sw_output%swup_special(i,j,1,nz) = swup_trop(i,j)
              if (do_clear_sky_pass) then
                swdn_trop_clr(i,j) = wtlo*Sw_output%dfswcf(i,j,ki,nz) +&
                                     wthi*Sw_output%dfswcf(i,j,ki-1,nz)
                swup_trop_clr(i,j) = wtlo*Sw_output%ufswcf(i,j,ki,nz) +&
                                     wthi*Sw_output%ufswcf(i,j,ki-1,nz)
                Sw_output%swdn_special_clr(i,j,1,nz) =   &
                                               swdn_trop_clr(i,j)
                Sw_output%swup_special_clr(i,j,1,nz) =    &
                                               swup_trop_clr(i,j)
              endif
              end do
            endif
            exit
          endif
          enddo
        enddo
        enddo
      endif
      if (linear_tropo) then
        lat_deg(:,:) = lat(:,:)*RADIAN
        tropo_ht(:,:) = trop_ht_at_eq + ABS(lat_deg(:,:))*  &
                        (trop_ht_at_poles - trop_ht_at_eq)/90.
! interpolate the fluxes between the appropriate pressures bracketing
! (trop) 

        do i = 1,ie-is+1
         do j = 1,je-js+1
          do k = kmax+1,2,-1
          if (Atmos_input%pflux(i,j,k) >= tropo_ht(i,j) .and.     &
              Atmos_input%pflux(i,j,k-1) < tropo_ht(i,j))      then
            ki = k
!   the indices for high,low pressure bracketing "tropo_ht" are ki, ki-1
            wtlo = (tropo_ht(i,j) - Atmos_input%pflux(i,j,ki-1))/ &
              (Atmos_input%pflux(i,j,ki) - Atmos_input%pflux(i,j,ki-1))
            wthi = 1.0 - wtlo
            if (Rad_control%do_lw_rad) then
            netlw_trop(i,j) = wtlo*Lw_output%flxnet(i,j,ki) + &
                              wthi*Lw_output%flxnet(i,j,ki-1)
            Lw_output%netlw_special(i,j,2) = netlw_trop(i,j)
            if (do_clear_sky_pass) then
              netlw_trop_clr(i,j) = wtlo*Lw_output%flxnetcf(i,j,ki) + &
                                    wthi*Lw_output%flxnetcf(i,j,ki-1)
              Lw_output%netlw_special_clr(i,j,2) = netlw_trop_clr(i,j)
            endif
            endif
            if (Rad_control%do_sw_rad) then
              do nz = 1,Rad_control%nzens
            swdn_trop(i,j) = wtlo*Sw_output%dfsw(i,j,ki,nz) + &
                             wthi*Sw_output%dfsw(i,j,ki-1,nz)
            swup_trop(i,j) = wtlo*Sw_output%ufsw(i,j,ki,nz) + &
                             wthi*Sw_output%ufsw(i,j,ki-1,nz)
            Sw_output%swdn_special(i,j,2,nz) = swdn_trop(i,j)
            Sw_output%swup_special(i,j,2,nz) = swup_trop(i,j)
            if (do_clear_sky_pass) then
              swdn_trop_clr(i,j) = wtlo*Sw_output%dfswcf(i,j,ki,nz) + &
                                   wthi*Sw_output%dfswcf(i,j,ki-1,nz)
              swup_trop_clr(i,j) = wtlo*Sw_output%ufswcf(i,j,ki,nz) + &
                                   wthi*Sw_output%ufswcf(i,j,ki-1,nz)
              Sw_output%swdn_special_clr(i,j,2,nz) = swdn_trop_clr(i,j)
              Sw_output%swup_special_clr(i,j,2,nz) = swup_trop_clr(i,j)
            endif
            end do
            endif
            exit
          endif
          enddo
        enddo
        enddo
      endif
      if (thermo_tropo) then
        call error_mesg ( 'radiation_driver_mod', &
              'thermo_tropo option not yet available', FATAL)
! interpolate the fluxes between the appropriate pressures bracketing
! (trop) 

        do i = 1,ie-is+1
         do j = 1,je-js+1
          do k = kmax+1,2,-1
          if (Atmos_input%pflux(i,j,k) >= tropo_ht(i,j) .and.     &
              Atmos_input%pflux(i,j,k-1) < tropo_ht(i,j))      then
            ki = k
!   the indices for high,low pressure bracketing "tropo_ht" are ki, ki-1
            wtlo = (tropo_ht(i,j) - Atmos_input%pflux(i,j,ki-1))/ &
              (Atmos_input%pflux(i,j,ki) - Atmos_input%pflux(i,j,ki-1))
            wthi = 1.0 - wtlo
            if (Rad_control%do_lw_rad) then
            netlw_trop(i,j) = wtlo*Lw_output%flxnet(i,j,ki) + &
                              wthi*Lw_output%flxnet(i,j,ki-1)
            Lw_output%netlw_special(i,j,3) = netlw_trop(i,j)
            if (do_clear_sky_pass) then
              netlw_trop_clr(i,j) = wtlo*Lw_output%flxnetcf(i,j,ki) + &
                                    wthi*Lw_output%flxnetcf(i,j,ki-1)
              Lw_output%netlw_special_clr(i,j,3) = netlw_trop_clr(i,j)
             endif
             endif
            if (Rad_control%do_sw_rad) then
              do nz = 1,Rad_control%nzens
            swdn_trop(i,j) = wtlo*Sw_output%dfsw(i,j,ki,nz) + &
                             wthi*Sw_output%dfsw(i,j,ki-1,nz)
            swup_trop(i,j) = wtlo*Sw_output%ufsw(i,j,ki,nz) + &
                             wthi*Sw_output%ufsw(i,j,ki-1,nz)
            Sw_output%swdn_special(i,j,3,nz) = swdn_trop(i,j)
            Sw_output%swup_special(i,j,3,nz) = swup_trop(i,j)
            if (do_clear_sky_pass) then
              swdn_trop_clr(i,j) = wtlo*Sw_output%dfswcf(i,j,ki,nz) + &
                                   wthi*Sw_output%dfswcf(i,j,ki-1,nz)
              swup_trop_clr(i,j) = wtlo*Sw_output%ufswcf(i,j,ki,nz) + &
                                   wthi*Sw_output%ufswcf(i,j,ki-1,nz)
              Sw_output%swdn_special_clr(i,j,3,nz) = swdn_trop_clr(i,j)
              Sw_output%swup_special_clr(i,j,3,nz) = swup_trop_clr(i,j)
            endif
            end do
            endif
            exit
          endif
          enddo
        enddo
        enddo
      endif

end subroutine flux_trop_calc 


!####################################################################
! <SUBROUTINE NAME="produce_radiation_diagnostics">
!  <OVERVIEW>
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call produce_radiation_diagnostics          &
!                            (is, ie, js, je, Time_diag, lat, ts, asfc, &
!                             flux_ratio, Astro, Rad_output, Lw_output, &
!                             Sw_output, Cld_spec, Lsc_microphys,    &
!                             Fsrad_output, mask)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 

!  <IN NAME="Time_diag" TYPE="time_type">
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>

!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="ts" TYPE="real">
!   Surface skin temperature
!  </IN>
!  <IN NAME="asfc" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="flux_ratio" TYPE="real">
!   renormalization factor for sw fluxes and heating 
!                   rates
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical input data for the radiation package
!  </IN>
!  <IN NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud microphysical and physical parameters to radiation package, 
!                        when the microphysical package is active
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale clouds,
!                        when the microphysical package is active
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine produce_radiation_diagnostics          &
                 (is, ie, js, je, Time_diag, Time, lat, ts, Surface, &
                  flux_ratio, Astro, Rad_output, Rad_gases,&
                  Lw_output, Sw_output, Cld_spec,   &
                  Lsc_microphys, Fsrad_output, mask)

!--------------------------------------------------------------------
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                 intent(in)             :: is, ie, js, je
type(time_type),         intent(in)             :: Time_diag
type(time_type),         intent(in)             :: Time
real,dimension(:,:),     intent(in)             :: lat, ts
type(surface_type),      intent(in)             :: Surface
real,dimension(:,:),     intent(in)             :: flux_ratio
type(astronomy_type),    intent(in)             :: Astro
type(rad_output_type),   intent(in)             :: Rad_output
type(radiative_gases_type), intent(in)          :: Rad_gases
type(lw_output_type), dimension(:), intent(in), optional   :: Lw_output
type(fsrad_output_type), intent(in), optional   :: Fsrad_output
type(sw_output_type),  dimension(:),  intent(in), optional :: Sw_output
type(cld_specification_type), intent(in), optional   :: Cld_spec
type(microphysics_type), intent(in), optional   :: Lsc_microphys
real,dimension(:,:,:),   intent(in), optional   :: mask
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      Time_diag    time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]  
!      lat          latitude of model points  [ radians ]
!      ts           surface temperature  [ deg K ]
!      asfc         surface albedo  [ dimensionless ]
!      flux_ratio   renormalization factor for sw fluxes and heating 
!                   rates [ dimensionless ]
!      Astro        astronomical  variables input to the radiation
!                   package [ dimensionless ]
!      Rad_output   rad_output_type variable containing radiation 
!                   output fields
!      Rad_gases    radiative_gases_type variable containing co2 mixing
!                   ratio
!
!
!    intent(in) optional variables:
!
!      Lw_output    lw_output_type variable containing output from 
!                   the longwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Sw_output    sw_output_type variable containing output from 
!                   the shortwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Cld_spec     cloud specification input data for the 
!                   radiation package
!                   [ cld_specification_type ]
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds, passed through to lower
!                     level routines
!      Fsrad_output fsrad_output_type variable containing 
!                   output from the original_fms_rad radiation
!                   package, on the model grid
!      mask         present when running eta vertical coordinate,
!                   mask to remove points below ground
!        
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables

      real, dimension (ie-is+1,je-js+1) ::           & 
                                                swin, swout, olr, &
                                                swups, swdns, lwups, &
                                                lwdns, swin_clr,   &
                                                swout_clr, olr_clr, &
                                                swups_clr, swdns_clr,&
                                                lwups_clr, lwdns_clr   

      real, dimension (ie-is+1,je-js+1, MX_SPEC_LEVS) ::           & 
                                                swdn_trop,  &
                                                swdn_trop_clr, &
                                                swup_trop, &
                                                swup_trop_clr, &
                                                netlw_trop, &
                                                netlw_trop_clr

      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)) ::  &
                                                tdtlw, tdtlw_clr,&
                                                hsw, hswcf

      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)+1) :: &
                                                dfsw, ufsw,  &
                                                dfswcf, ufswcf,&
                                                flxnet, flxnetcf, &
                                                fsw, fswcf
      real, dimension (ie-is+1,je-js+1) ::      &
                                         swin_ad,     swout_ad, olr_ad,&
                              swups_ad,    swdns_ad, lwups_ad,lwdns_ad,&
                                 swin_ad_clr, swout_ad_clr, olr_ad_clr,&
                  swups_ad_clr, swdns_ad_clr, lwups_ad_clr, lwdns_ad_clr
      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)+1) :: &
                                               dfsw_ad, ufsw_ad,  &
                                               dfswcf_ad, ufswcf_ad

      integer           :: j, k
      integer           :: ipass
      logical           :: used
      integer           :: iind, jind
      integer           :: kmax
      integer           :: nz

!      asfc         surface albedo  [ dimensionless ]
!      asfc_vis_dir surface visible albedo  [ dimensionless ]
!      asfc_nir_dir surface nir albedo  [ dimensionless ]
!      asfc_vis_dif surface visible albedo  [ dimensionless ]
!      asfc_nir_dif surface nir albedo  [ dimensionless ]

!---------------------------------------------------------------------
!    if sw flux renormalization is active, modify the fluxes calculated
!    on the last radiation step by the normalization factor based on
!    the difference in solar factor between the current model step and
!    the current radiation step.
!----------------------------------------------------------------------
      nz = current_sw_zenith_step
      kmax = size (Rad_output%tdtsw,3)
      if (renormalize_sw_fluxes) then
        do k=1, kmax         
          hsw(:,:,k) = hsw_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
        end do
        do k=1, kmax+1             
          if (do_swaerosol_forcing) then
            dfsw_ad(:,:,k) =   &
                        dfsw_ad_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
            ufsw_ad(:,:,k) =   &
                        ufsw_ad_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          endif
          dfsw(:,:,k) = dfsw_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          ufsw(:,:,k) = ufsw_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          fsw(:,:,k) = fsw_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
        end do
        do k=1,Rad_control%mx_spec_levs
          swdn_trop(:,:,k) = swdn_special_save(is:ie,js:je,k,nz)*  &
                             flux_ratio(:,:)
          swup_trop(:,:,k) = swup_special_save(is:ie,js:je,k,nz)*   &
                             flux_ratio(:,:)
        end do
        if (do_clear_sky_pass) then
          do k=1, kmax            
            hswcf(:,:,k) = hswcf_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          end do
          do k=1, kmax+1            
            if (do_swaerosol_forcing) then
              dfswcf_ad(:,:,k) = dfswcf_ad_save(is:ie,js:je,k,nz)*  &
                                 flux_ratio(:,:)
              ufswcf_ad(:,:,k) = ufswcf_ad_save(is:ie,js:je,k,nz)*  &
                                 flux_ratio(:,:)
            endif
            dfswcf(:,:,k) =    &
                          dfswcf_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
            ufswcf(:,:,k) =    &
                          ufswcf_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
            fswcf(:,:,k) = fswcf_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          end do
          do k=1,Rad_control%mx_spec_levs
            swdn_trop_clr(:,:,k) =    &
                  swdn_special_clr_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
            swup_trop_clr(:,:,k) =    &
                  swup_special_clr_save(is:ie,js:je,k,nz)*flux_ratio(:,:)
          end do
        endif

!----------------------------------------------------------------------
!    if renormalization is not active and this is a radiation step
!    (i.e., diagnostics desired), define the variables to be output as
!    the values present in Sw_output.
!---------------------------------------------------------------------
      else if (do_sw_rad .and. do_sea_esf_rad) then
        do k=1, kmax            
          hsw(:,:,k) = Sw_output(1)%hsw(:,:,k,nz)
        end do
        do k=1, kmax+1             
          if (do_swaerosol_forcing) then
            dfsw_ad(:,:,k) = Sw_output(indx_swaf)%dfsw(:,:,k,nz)
            ufsw_ad(:,:,k) = Sw_output(indx_swaf)%ufsw(:,:,k,nz)
          endif
          dfsw(:,:,k) = Sw_output(1)%dfsw(:,:,k,nz)
          ufsw(:,:,k) = Sw_output(1)%ufsw(:,:,k,nz)
          fsw(:,:,k) = Sw_output(1)%fsw(:,:,k,nz)
        end do
        swdn_trop(:,:,:) = Sw_output(1)%swdn_special(:,:,:,nz)
        swup_trop(:,:,:) = Sw_output(1)%swup_special(:,:,:,nz)
        if (do_clear_sky_pass) then
          do k=1, kmax             
            hswcf(:,:,k) = Sw_output(1)%hswcf(:,:,k,nz)
          end do
          do k=1, kmax+1            
            if (do_swaerosol_forcing) then
              dfswcf_ad(:,:,k) = Sw_output(indx_swaf)%dfswcf(:,:,k,nz)
              ufswcf_ad(:,:,k) = Sw_output(indx_swaf)%ufswcf(:,:,k,nz)  
            endif
            dfswcf(:,:,k) = Sw_output(1)%dfswcf(:,:,k,nz)
            ufswcf(:,:,k) = Sw_output(1)%ufswcf(:,:,k,nz)
            fswcf(:,:,k) = Sw_output(1)%fswcf(:,:,k,nz)
          end do
          swdn_trop_clr(:,:,:) = Sw_output(1)%swdn_special_clr(:,:,:,nz)
          swup_trop_clr(:,:,:) = Sw_output(1)%swup_special_clr(:,:,:,nz)
        endif

!----------------------------------------------------------------------
!    if renormalization is not active and this is not a radiation step
!    but all_step_diagnostics is activated (i.e., diagnostics desired),
!    define the variables to be output as the values previously saved
!    in the xxx_save variables.
!---------------------------------------------------------------------
      else if (do_sea_esf_rad .and. all_step_diagnostics) then
        do k=1, kmax
          hsw(:,:,k) = hsw_save(is:ie,js:je,k,nz)
        end do
        do k=1, kmax+1
          if (do_swaerosol_forcing) then
            dfsw_ad(:,:,k) = dfsw_ad_save(is:ie,js:je,k,nz)
            ufsw_ad(:,:,k) = ufsw_ad_save(is:ie,js:je,k,nz)
          endif
          dfsw(:,:,k) = dfsw_save(is:ie,js:je,k,nz)
          ufsw(:,:,k) = ufsw_save(is:ie,js:je,k,nz)
          fsw(:,:,k) = fsw_save(is:ie,js:je,k,nz)
        end do
        swdn_trop(:,:,:) = swdn_special_save(is:ie,js:je,:,nz)
        swup_trop(:,:,:) = swup_special_save(is:ie,js:je,:,nz)
        if (do_clear_sky_pass) then
          do k=1, kmax
            hswcf(:,:,k) = hswcf_save(is:ie,js:je,k,nz)
          end do
          do k=1, kmax+1
            if (do_swaerosol_forcing) then
              dfswcf_ad(:,:,k) = dfswcf_ad_save(is:ie,js:je,k,nz)
              ufswcf_ad(:,:,k) = ufswcf_ad_save(is:ie,js:je,k,nz)
            endif
            dfswcf(:,:,k) = dfswcf_save(is:ie,js:je,k,nz)
            ufswcf(:,:,k) = ufswcf_save(is:ie,js:je,k,nz)
            fswcf(:,:,k) = fswcf_save(is:ie,js:je,k,nz)
          end do
          swdn_trop_clr(:,:,:) = swdn_special_clr_save(is:ie,js:je,:,nz)
          swup_trop_clr(:,:,:) = swup_special_clr_save(is:ie,js:je,:,nz)
        endif
      endif

!---------------------------------------------------------------------
!    define the sw diagnostic arrays.
!---------------------------------------------------------------------
      if (renormalize_sw_fluxes .or. do_sw_rad .or.    &
           use_hires_coszen .or.   all_step_diagnostics) then
        if (do_sea_esf_rad) then
          if (do_swaerosol_forcing) then
            swin_ad (:,:) = dfsw_ad(:,:,1)
            swout_ad(:,:) = ufsw_ad(:,:,1)
            swups_ad(:,:) = ufsw_ad(:,:,kmax+1)
            swdns_ad(:,:) = dfsw_ad(:,:,kmax+1)
          endif
          swin (:,:) = dfsw(:,:,1)
          swout(:,:) = ufsw(:,:,1)
          swups(:,:) = ufsw(:,:,kmax+1)
          swdns(:,:) = dfsw(:,:,kmax+1)
          if (do_clear_sky_pass) then
            if (do_swaerosol_forcing) then
              swin_ad_clr (:,:) = dfswcf_ad(:,:,1)
              swout_ad_clr(:,:) = ufswcf_ad(:,:,1)
              swups_ad_clr(:,:) = ufswcf_ad(:,:,kmax+1)
              swdns_ad_clr(:,:) = dfswcf_ad(:,:,kmax+1)
            endif
            swin_clr (:,:) = dfswcf(:,:,1)
            swout_clr(:,:) = ufswcf(:,:,1)
            swups_clr(:,:) = ufswcf(:,:,kmax+1)
            swdns_clr(:,:) = dfswcf(:,:,kmax+1)
          endif
        else   ! original fms rad
          swin (:,:) = Fsrad_output%swin(:,:)               
          swout(:,:) = Fsrad_output%swout(:,:)         
          swups(:,:) = Fsrad_output%swups(:,:)
          swdns(:,:) = Fsrad_output%swdns(:,:)
          if (do_clear_sky_pass) then
            swin_clr (:,:) = Fsrad_output%swin_clr(:,:)               
            swout_clr(:,:) = Fsrad_output%swout_clr(:,:)         
            swups_clr(:,:) = Fsrad_output%swups_clr(:,:)
            swdns_clr(:,:) = Fsrad_output%swdns_clr(:,:)
          endif
        endif  ! do_sea_esf_rad

        if (id_alb_sfc_avg > 0) then
          swups_acc(is:ie,js:je) = swups_acc(is:ie, js:je) + swups(:,:)
          swdns_acc(is:ie,js:je) = swdns_acc(is:ie, js:je) + swdns(:,:)
        endif
 
!---------------------------------------------------------------------
!   send standard sw diagnostics to diag_manager.
!---------------------------------------------------------------------
      if (Time_diag > Time) then
        if (do_clear_sky_pass) then
          ipass = 2
        else
          ipass = 1
        endif

!------- sw tendency -----------
        if (id_tdt_sw(ipass) > 0 ) then
          used = send_data (id_tdt_sw(ipass),    &
                            Rad_output%tdtsw(is:ie,js:je,:,nz),   &
                            Time_diag, is, js, 1, rmask=mask )
        endif

!---- 3d upward sw flux ---------
        if (id_ufsw(ipass) > 0 ) then
          used = send_data (id_ufsw(ipass),    &
                            Rad_output%ufsw(is:ie,js:je,:,nz),   &
                            Time_diag, is, js, 1, rmask=mask )
        endif

!---- 3d downward sw flux ---------
        if (id_dfsw(ipass) > 0 ) then
          used = send_data (id_dfsw(ipass),    &
                            Rad_output%dfsw(is:ie,js:je,:,nz),   &
                            Time_diag, is, js, 1, rmask=mask )
        endif


!------- incoming sw flux toa -------
        if (id_swdn_toa(ipass) > 0 ) then
          used = send_data (id_swdn_toa(ipass), swin,   &
                            Time_diag, is, js )
        endif

!------- outgoing sw flux toa -------
        if (id_swup_toa(ipass) > 0 ) then
          used = send_data (id_swup_toa(ipass), swout,    &
                            Time_diag, is, js )
        endif

!------- incoming sw flux trop -------
        if (id_swdn_special(1,ipass) > 0 ) then
          used = send_data (id_swdn_special(1,ipass),   &
                            swdn_trop(:,:,1), &
                            Time_diag, is, js )
        endif
!------- incoming sw flux trop -------
        if (id_swdn_special(2,ipass) > 0 ) then
          used = send_data (id_swdn_special(2,ipass),   &
                            swdn_trop(:,:,2), &
                            Time_diag, is, js )
        endif
!------- incoming sw flux trop -------
        if (id_swdn_special(3, ipass) > 0 ) then
          used = send_data (id_swdn_special(3,ipass),   &
                            swdn_trop(:,:,3), &
                            Time_diag, is, js )
        endif

!------- net sw downward flux at model dynamics top (1 Pa) ----
        if (id_swdn_special(4, ipass) > 0 ) then
          used = send_data (id_swdn_special(4,ipass),   &
                            swdn_trop(:,:,4), &
                            Time_diag, is, js )
        endif


!------- outgoing sw flux trop -------
        if (id_swup_special(1,ipass) > 0 ) then
          used = send_data (id_swup_special(1,ipass),   &
                            swup_trop(:,:,1), &
                            Time_diag, is, js )
        endif
!------- outgoing sw flux trop -------
        if (id_swup_special(2,ipass) > 0 ) then
          used = send_data (id_swup_special(2,ipass),   &
                            swup_trop(:,:,2), &
                            Time_diag, is, js )
        endif
!------- outgoing sw flux trop -------
        if (id_swup_special(3,ipass) > 0 ) then
          used = send_data (id_swup_special(3,ipass),   &
                            swup_trop(:,:,3), &
                            Time_diag, is, js )
        endif
 
!------- net sw upward flux at model dynamics top (1 Pa) ----
        if (id_swdn_special(4, ipass) > 0 ) then
          used = send_data (id_swup_special(4,ipass),   &
                            swup_trop(:,:,4), &
                            Time_diag, is, js )
        endif

!------- upward sw flux surface -------
        if (id_swup_sfc(ipass) > 0 ) then
          used = send_data (id_swup_sfc(ipass), swups,    &
                            Time_diag, is, js )
        endif

!------- downward sw flux surface -------
        if (id_swdn_sfc(ipass) > 0 ) then
          used = send_data (id_swdn_sfc(ipass), swdns,   &
                            Time_diag, is, js )
        endif
        
!------- net sw flux at toa -------
       if (id_swtoa(ipass) > 0 ) then
          used = send_data (id_swtoa(ipass), swin-swout,   &
                             Time_diag, is, js )
         endif

!------- net sw flux at surface -------
         if (id_swsfc(ipass) > 0 ) then
           used = send_data (id_swsfc(ipass), swdns-swups,   &
                             Time_diag, is, js )
         endif

      if (do_swaerosol_forcing) then

!------- net sw flux at toa -------
         if (id_swtoa_ad(ipass) > 0 ) then
           used = send_data (id_swtoa_ad(ipass), swin_ad-swout_ad,   &
                             Time_diag, is, js )
        endif
 
!------- net sw flux at surface -------
         if (id_swsfc_ad(ipass) > 0 ) then
           used = send_data (id_swsfc_ad(ipass), swdns_ad-swups_ad,   &
                           Time_diag, is, js )
         endif
 
!------- sw flux down at surface -------
         if (id_swdn_sfc_ad(ipass) > 0 ) then
           used = send_data (id_swdn_sfc_ad(ipass), swdns_ad,   &
                             Time_diag, is, js )
         endif

!------- sw flux up at surface -------
         if (id_swup_sfc_ad(ipass) > 0 ) then
           used = send_data (id_swup_sfc_ad(ipass), swups_ad,   &
                            Time_diag, is, js )
         endif

!------- outgoing sw flux toa -------
         if (id_swup_toa_ad(ipass) > 0 ) then
           used = send_data (id_swup_toa_ad(ipass), swout_ad,    &
                            Time_diag, is, js )
        endif
     endif

!----------------------------------------------------------------------
!    now pass clear-sky diagnostics, if they have been calculated.
!----------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 1

!------- sw tendency -----------
          if (id_tdt_sw(ipass) > 0 ) then
            used = send_data (id_tdt_sw(ipass),   &
                              Rad_output%tdtsw_clr(is:ie,js:je,:,nz),  &
                              Time_diag, is, js, 1, rmask=mask )
          endif

!---- 3d upward sw flux ---------
         if (id_ufsw(ipass) > 0 ) then
           used = send_data (id_ufsw(ipass),    &
                             Rad_output%ufsw_clr(is:ie,js:je,:,nz),   &
                             Time_diag, is, js, 1, rmask=mask )
         endif
 
!---- 3d downward sw flux ---------
         if (id_dfsw(ipass) > 0 ) then
           used = send_data (id_dfsw(ipass),    &
                             Rad_output%dfsw_clr(is:ie,js:je,:,nz),   &
                             Time_diag, is, js, 1, rmask=mask )
         endif

!------- incoming sw flux toa -------
          if (id_swdn_toa(ipass) > 0 ) then
            used = send_data (id_swdn_toa(ipass), swin_clr,    &
                              Time_diag, is, js )
          endif

!------- outgoing sw flux toa -------
          if (id_swup_toa(ipass) > 0 ) then
            used = send_data (id_swup_toa(ipass), swout_clr,  &
                              Time_diag, is, js )
          endif

!------- incoming sw flux trop -------
          if (id_swdn_special(1,ipass) > 0 ) then
            used = send_data (id_swdn_special(1, ipass),    &
                              swdn_trop_clr(:,:,1), &
                              Time_diag, is, js )
          endif
!------- incoming sw flux trop -------
          if (id_swdn_special(2,ipass) > 0 ) then
            used = send_data (id_swdn_special(2, ipass),    &
                              swdn_trop_clr(:,:,2), &
                              Time_diag, is, js )
          endif
!------- incoming sw flux trop -------
          if (id_swdn_special(3,ipass) > 0 ) then
            used = send_data (id_swdn_special(3, ipass),    &
                              swdn_trop_clr(:,:,3), &
                              Time_diag, is, js )
          endif

!------- outgoing sw flux trop -------
          if (id_swup_special(1,ipass) > 0 ) then
            used = send_data (id_swup_special(1,ipass),   &
                              swup_trop_clr(:,:,1),  &
                              Time_diag, is, js )
          endif
!------- outgoing sw flux trop -------
          if (id_swup_special(2,ipass) > 0 ) then
            used = send_data (id_swup_special(2,ipass),   &
                              swup_trop_clr(:,:,2),  &
                              Time_diag, is, js )
          endif
!------- outgoing sw flux trop -------
          if (id_swup_special(3,ipass) > 0 ) then
            used = send_data (id_swup_special(3,ipass),   &
                              swup_trop_clr(:,:,3),  &
                              Time_diag, is, js )
          endif

!------- upward sw flux surface -------
          if (id_swup_sfc(ipass) > 0 ) then
            used = send_data (id_swup_sfc(ipass), swups_clr,   &
                              Time_diag, is, js )
          endif

!------- downward sw flux surface -------
          if (id_swdn_sfc(ipass) > 0 ) then
            used = send_data (id_swdn_sfc(ipass), swdns_clr,    &
                              Time_diag, is, js )
          endif

!------- net sw flux at toa -------
        if (id_swtoa(ipass) > 0 ) then
           used = send_data (id_swtoa(ipass), swin_clr-swout_clr,   &
                             Time_diag, is, js )
         endif
 
!------- net sw flux at surface -------
         if (id_swsfc(ipass) > 0 ) then
           used = send_data (id_swsfc(ipass), swdns_clr-swups_clr,   &
                            Time_diag, is, js )
         endif
     if (do_swaerosol_forcing) then

!------- net sw flux at toa -------
        if (id_swtoa_ad(ipass) > 0 ) then
          used = send_data (id_swtoa_ad(ipass), swin_ad_clr-swout_ad_clr,   &
                            Time_diag, is, js )
        endif
 
!------- net sw flux at surface -------
         if (id_swsfc_ad(ipass) > 0 ) then
           used = send_data (id_swsfc_ad(ipass), swdns_ad_clr-swups_ad_clr,   &
                             Time_diag, is, js )
         endif
 
!------- sw flux down at surface -------
         if (id_swdn_sfc_ad(ipass) > 0 ) then
           used = send_data (id_swdn_sfc_ad(ipass), swdns_ad_clr,   &
                             Time_diag, is, js )
        endif
 
!------- sw flux up at surface -------
        if (id_swup_sfc_ad(ipass) > 0 ) then
          used = send_data (id_swup_sfc_ad(ipass), swups_ad_clr,   &
                            Time_diag, is, js )
        endif
 
!------- outgoing sw flux toa -------
        if (id_swup_toa_ad(ipass) > 0 ) then
           used = send_data (id_swup_toa_ad(ipass), swout_ad_clr,    &
                            Time_diag, is, js )
        endif

     endif
         endif  ! (do_clear_sky_pass)

!-----------------------------------------------------------------------
!    send cloud-forcing-independent diagnostics to diagnostics manager.
!-----------------------------------------------------------------------

!---- 3d total radiative heating ---------
        if (id_allradp > 0 ) then
          used = send_data (id_allradp    ,    &
                            Rad_output%tdt_rad(is:ie,js:je,:,nz),   &
                            Time_diag, is, js, 1, rmask=mask )
        endif

!------- conc_drop  -------------------------
          if (do_rad) then
            if ( id_conc_drop > 0 ) then
              used = send_data (id_conc_drop, Lsc_microphys%conc_drop, &
                                Time_diag, is, js, 1, rmask=mask )
            endif
          endif
  
!------- conc_ice  -------------------------
          if (do_rad) then
            if (id_conc_ice > 0 ) then
              used = send_data (id_conc_ice, Lsc_microphys%conc_ice, &
                                Time_diag, is, js, 1, rmask=mask )
            endif
          endif

!------- solar constant  -------------------------
        if (do_rad) then
          if ( id_sol_con > 0 ) then
            used = send_data ( id_sol_con, Sw_control%solar_constant,  &
                               Time_diag )
          endif
        endif

!------- co2 mixing ratio used for tf calculation  -------------------
        if (do_rad) then
          if ( id_co2_tf > 0 ) then
            used = send_data ( id_co2_tf,   &
                               1.0E6*Rad_gases%co2_for_last_tf_calc,  &
                               Time_diag )
          endif
        endif
 
!------- ch4 mixing ratio used for tf calculation   ---------------
        if (do_rad) then
          if ( id_ch4_tf > 0 ) then
            used = send_data ( id_ch4_tf,  &
                               1.0E9*Rad_gases%ch4_for_last_tf_calc,  &
                               Time_diag )
          endif
        endif
 
!------- n2o mixing ratio used for tf calculation  ---------------
        if (do_rad) then
          if ( id_n2o_tf > 0 ) then
            used = send_data ( id_n2o_tf,   &
                               1.0E9*Rad_gases%n2o_for_last_tf_calc,  &
                               Time_diag )
          endif
        endif
 
!------- co2 mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvco2 > 0 ) then
            used = send_data ( id_rrvco2, 1.0E6*Rad_gases%rrvco2,  &
                               Time_diag )
          endif
        endif
 
!------- f11 mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvf11 > 0 ) then
            used = send_data ( id_rrvf11, 1.0E12*Rad_gases%rrvf11,  &
                               Time_diag )
          endif
        endif
 
!------- f12 mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvf12 > 0 ) then
            used = send_data ( id_rrvf12, 1.0E12*Rad_gases%rrvf12,  &
                               Time_diag )
          endif
        endif
 
!------- f113 mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvf113 > 0 ) then
            used = send_data ( id_rrvf113, 1.0E12*Rad_gases%rrvf113,  &
                               Time_diag )
          endif
        endif

!------- f22 mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvf22 > 0 ) then
            used = send_data ( id_rrvf22, 1.0E12*Rad_gases%rrvf22,  &
                               Time_diag )
          endif
        endif

!------- ch4 mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvch4 > 0 ) then
            used = send_data ( id_rrvch4, 1.0E9*Rad_gases%rrvch4,  &
                               Time_diag )
          endif
        endif

!------- n2o mixing ratio  -------------------------
        if (do_rad) then
          if ( id_rrvn2o > 0 ) then
            used = send_data ( id_rrvn2o, 1.0E9*Rad_gases%rrvn2o,  &
                               Time_diag )
          endif
        endif

!------- surface albedo  -------------------------
        if ( id_alb_sfc_avg > 0 ) then
!         used = send_data ( id_alb_sfc, 100.*Surface%asfc, &
          used = send_data ( id_alb_sfc_avg,  &
                  100.*swups_acc(is:ie,js:je)/swdns_acc(is:ie,js:je), &
                     Time_diag, is, js )
        endif
        if ( id_alb_sfc > 0 ) then
!         used = send_data ( id_alb_sfc, 100.*Surface%asfc, &
          used = send_data ( id_alb_sfc, 100.*swups/(1.0e-35 + swdns), &
                     Time_diag, is, js )
        endif

!------- surface visible albedo  -------------------------
        if ( id_alb_sfc_vis_dir > 0 ) then
          used = send_data ( id_alb_sfc_vis_dir, &
                         100.*Surface%asfc_vis_dir, Time_diag, is, js )
        endif
        if ( id_alb_sfc_vis_dif > 0 ) then
          used = send_data ( id_alb_sfc_vis_dif, &
                         100.*Surface%asfc_vis_dif, Time_diag, is, js )
        endif
 
!------- surface nir albedo  -------------------------
        if ( id_alb_sfc_nir_dir > 0 ) then
          used = send_data ( id_alb_sfc_nir_dir, &
                         100.*Surface%asfc_nir_dir, Time_diag, is, js )
        endif
        if ( id_alb_sfc_nir_dif > 0 ) then
           used = send_data ( id_alb_sfc_nir_dif, &
                         100.*Surface%asfc_nir_dif, Time_diag, is, js )
        endif
 
!------- surface net sw flux, direct and diffuse  --------------------
        if ( id_flux_sw > 0 ) then
         used = send_data ( id_flux_sw, &
          (Rad_output%flux_sw_surf_dir( is:ie,js:je,nz) + &
          Rad_output%flux_sw_surf_dif( is:ie,js:je,nz)),   &
                                         Time_diag, is, js )
        endif

        if ( id_flux_sw_dir > 0 ) then
         used = send_data ( id_flux_sw_dir, &
          Rad_output%flux_sw_surf_dir( is:ie,js:je,nz), Time_diag,  &
                                                              is, js )
        endif
        
        if ( id_flux_sw_refl_dir > 0 ) then
         used = send_data ( id_flux_sw_refl_dir, &
          Rad_output%flux_sw_surf_refl_dir( is:ie,js:je,nz), Time_diag,  &
                                                              is, js )
        endif

        if ( id_flux_sw_refl_vis_dir > 0 ) then
         used = send_data ( id_flux_sw_refl_vis_dir, &
          Rad_output%flux_sw_refl_vis_dir( is:ie,js:je,nz), Time_diag,  &
                                                              is, js )
        endif
        if ( id_flux_sw_dif > 0 ) then
          used = send_data ( id_flux_sw_dif, &
           Rad_output%flux_sw_surf_dif(is:ie,js:je,nz), Time_diag, &
                                                              is, js )
        endif

!------- surface downward visible sw flux, direct and diffuse ----------
        if ( id_flux_sw_down_vis_dir > 0 ) then
          used = send_data ( id_flux_sw_down_vis_dir, &
                     Rad_output%flux_sw_down_vis_dir(is:ie,js:je,nz), &
                     Time_diag, is, js )
        endif
        if ( id_flux_sw_down_vis_dif > 0 ) then
          used = send_data ( id_flux_sw_down_vis_dif, &
                     Rad_output%flux_sw_down_vis_dif(is:ie, js:je,nz), &
                     Time_diag, is, js )
        endif
 
!------- surface downward total sw flux, direct and diffuse  ----------
        if ( id_flux_sw_down_total_dir > 0 ) then
          used = send_data ( id_flux_sw_down_total_dir,  &
                     Rad_output%flux_sw_down_total_dir(is:ie,js:je,nz),  &
                     Time_diag, is, js )
        endif
        if ( id_flux_sw_down_total_dif > 0 ) then
         used = send_data ( id_flux_sw_down_total_dif,  &
                    Rad_output%flux_sw_down_total_dif(is:ie,js:je,nz),  &
                    Time_diag, is, js )
        endif

      if (do_clear_sky_pass) then
 
!------- surface downward total sw flux, direct and diffuse  ----------
        if ( id_flux_sw_down_total_dir_clr > 0 ) then
          used = send_data ( id_flux_sw_down_total_dir_clr,  &
                 Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je,nz),&
                 Time_diag, is, js )
        endif
        if ( id_flux_sw_down_total_dif_clr > 0 ) then
          used = send_data ( id_flux_sw_down_total_dif_clr,  &
                  Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je,nz),  &
                  Time_diag, is, js )
        endif
        if ( id_flux_sw_down_vis_clr > 0 ) then
          used = send_data ( id_flux_sw_down_vis_clr, &
                     Rad_output%flux_sw_down_vis_clr(is:ie, js:je,nz), &
                     Time_diag, is, js )
        endif
      endif
 
!------- surface net visible sw flux, total, direct and diffuse -------
        if ( id_flux_sw_vis > 0 ) then
          used = send_data ( id_flux_sw_vis,   &
                Rad_output%flux_sw_vis(is:ie,js:je,nz), Time_diag, is, js )
        endif
        if ( id_flux_sw_vis_dir > 0 ) then
          used = send_data ( id_flux_sw_vis_dir,   &
            Rad_output%flux_sw_vis_dir(is:ie,js:je,nz), Time_diag, is, js )
        endif
        if ( id_flux_sw_vis_dif > 0 ) then
          used = send_data ( id_flux_sw_vis_dif,  &
            Rad_output%flux_sw_vis_dif(is:ie,js:je,nz), Time_diag, is, js )
        endif

!------- cosine of zenith angle ----------------
        if ( id_cosz > 0 ) then
          used = send_data ( id_cosz, Astro%cosz, Time_diag, is, js )
        endif

!------- daylight fraction  --------------
        if ( id_fracday > 0 ) then
          used = send_data (id_fracday, Astro%fracday, Time_diag,   &
                            is, js )
        end if
      endif
      endif   ! (renormalize_sw_fluxes .or. do_rad .or.   
              !  all_step_diagnostics)

!---------------------------------------------------------------------
!    define the longwave diagnostic arrays for the sea-esf radiation 
!    package.  convert to mks units.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
!       if (do_rad) then
        if (do_lw_rad) then
          olr  (:,:)   = Lw_output(1)%flxnet(:,:,1)
          lwups(:,:)   =   STEFAN*ts(:,:  )**4
          lwdns(:,:)   = lwups(:,:) - Lw_output(1)%flxnet(:,:,kmax+1)
          tdtlw(:,:,:) = Lw_output(1)%heatra(:,:,:)/ SECONDS_PER_DAY
          netlw_trop(:,:,:) = Lw_output(1)%netlw_special(:,:,:)
          flxnet(:,:,:) = Lw_output(1)%flxnet(:,:,:)
          if (do_lwaerosol_forcing) then
            olr_ad  (:,:)   = Lw_output(indx_lwaf)%flxnet(:,:,1)
            lwups_ad(:,:)   = STEFAN*ts(:,:  )**4
            lwdns_ad(:,:)   = lwups_ad(:,:) -    &
                                  Lw_output(indx_lwaf)%flxnet(:,:,kmax+1)
          endif

          if (do_clear_sky_pass) then
            olr_clr  (:,:)   = Lw_output(1)%flxnetcf(:,:,1)
            lwups_clr(:,:)   =              STEFAN*ts(:,:  )**4
            lwdns_clr(:,:)   = lwups_clr(:,:) -    & 
                               Lw_output(1)%flxnetcf(:,:,kmax+1)
            tdtlw_clr(:,:,:) = Lw_output(1)%heatracf(:,:,:)/SECONDS_PER_DAY
            netlw_trop_clr(:,:,:) = Lw_output(1)%netlw_special_clr(:,:,:)
            flxnetcf(:,:,:) = Lw_output(1)%flxnetcf(:,:,:)
            if (do_lwaerosol_forcing) then
              olr_ad_clr  (:,:)   = Lw_output(indx_lwaf)%flxnetcf(:,:,1)
              lwups_ad_clr(:,:)   = STEFAN*ts(:,:  )**4
              lwdns_ad_clr(:,:)   = lwups_ad_clr(:,:) -    &
                               Lw_output(indx_lwaf)%flxnetcf(:,:,kmax+1)
            endif
          endif

!---------------------------------------------------------------------
!    if diagnostics are desired on all physics steps, save the arrays 
!    for later use.
!---------------------------------------------------------------------
          if (all_step_diagnostics) then
            if (do_lwaerosol_forcing) then
              olr_ad_save  (is:ie,js:je)   = olr_ad(:,:)
              lwups_ad_save(is:ie,js:je)   = lwups_ad(:,:)
              lwdns_ad_save(is:ie,js:je)   = lwdns_ad(:,:)
            endif
            olr_save  (is:ie,js:je)   = olr(:,:)
            lwups_save(is:ie,js:je)   = lwups(:,:)
            lwdns_save(is:ie,js:je)   = lwdns(:,:)
            tdtlw_save(is:ie,js:je,:) = tdtlw(:,:,:)
            flxnet_save(is:ie,js:je,:) = Lw_output(1)%flxnet(:,:,:)
            netlw_special_save(is:ie,js:je,:) = netlw_trop(:,:,:)

            if (do_clear_sky_pass) then
              if (do_lwaerosol_forcing) then 
                olr_ad_clr_save  (is:ie,js:je)   = olr_ad_clr(:,:)
                lwups_ad_clr_save(is:ie,js:je)   = lwups_ad_clr(:,:)
                lwdns_ad_clr_save(is:ie,js:je)   = lwdns_ad_clr(:,:)
              endif
              olr_clr_save  (is:ie,js:je)   = olr_clr(:,:)
              flxnetcf_save(is:ie,js:je,:) = Lw_output(1)%flxnetcf(:,:,:)
              lwups_clr_save(is:ie,js:je)   = lwups_clr(:,:)
              lwdns_clr_save(is:ie,js:je)   = lwdns_clr(:,:)
              tdtlw_clr_save(is:ie,js:je,:) = tdtlw_clr(:,:,:)
              netlw_special_clr_save(is:ie,js:je,:) =   &
                                         netlw_trop_clr(:,:,:)
            endif
           endif

!---------------------------------------------------------------------
!    if this is not a radiation step, but diagnostics are desired,
!    define the fields from the xxx_save variables.
!---------------------------------------------------------------------
!        else if (all_step_diagnostics) then  ! (do_rad)
         else if (all_step_diagnostics) then  ! (do_lw_rad)
           if (do_lwaerosol_forcing) then 
             olr_ad(:,:)     = olr_ad_save  (is:ie,js:je)
             lwups_ad(:,:)   = lwups_ad_save(is:ie,js:je)
             lwdns_ad(:,:)   = lwdns_ad_save(is:ie,js:je) 
           endif
           olr(:,:)     = olr_save  (is:ie,js:je)
           lwups(:,:)   = lwups_save(is:ie,js:je)
           lwdns(:,:)   = lwdns_save(is:ie,js:je)
           tdtlw(:,:,:) = tdtlw_save(is:ie,js:je,:)
           flxnet(:,:,:) = flxnet_save(is:ie,js:je,:)
           netlw_trop(:,:,:) = netlw_special_save(is:ie,js:je,:)

           if (do_clear_sky_pass) then
             if (do_lwaerosol_forcing) then
               olr_ad_clr(:,:)     = olr_ad_clr_save  (is:ie,js:je)
               lwups_ad_clr(:,:)   = lwups_ad_clr_save(is:ie,js:je)
               lwdns_ad_clr(:,:)   = lwdns_ad_clr_save(is:ie,js:je)
             endif
             olr_clr(:,:)     = olr_clr_save  (is:ie,js:je)
             lwups_clr(:,:)   = lwups_clr_save(is:ie,js:je)
             lwdns_clr(:,:)   = lwdns_clr_save(is:ie,js:je)
             tdtlw_clr(:,:,:) = tdtlw_clr_save(is:ie,js:je,:)
             flxnetcf (:,:,:) = flxnetcf_save(is:ie,js:je,:)
             netlw_trop_clr(:,:,:) =   &
                               netlw_special_clr_save(is:ie,js:je,:)
           endif
         endif

!---------------------------------------------------------------------
!    on radiation steps, define the longwave diagnostic arrays for the
!    original_fms_rad package.        
!---------------------------------------------------------------------
      else   ! original fms rad
        if (do_lw_rad) then
          olr  (:,:)   = Fsrad_output%olr(:,:)
          lwups(:,:)   = Fsrad_output%lwups(:,:)
          lwdns(:,:)   = Fsrad_output%lwdns(:,:)
          tdtlw(:,:,:) = Fsrad_output%tdtlw(:,:,:)

          if (do_clear_sky_pass) then
            olr_clr  (:,:)   = Fsrad_output%olr_clr(:,:)
            lwups_clr(:,:)   = Fsrad_output%lwups_clr(:,:)
            lwdns_clr(:,:)   = Fsrad_output%lwdns_clr(:,:)
            tdtlw_clr(:,:,:) = Fsrad_output%tdtlw_clr(:,:,:)
          endif
        endif
      endif  ! do_sea_esf_rad

      if (do_lw_rad .or. all_step_diagnostics) then
      if (Time_diag > Time) then
!---------------------------------------------------------------------
!   send standard lw diagnostics to diag_manager.
!---------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 2
        else
          ipass = 1
        endif

!---- net lw flux ---------
        if (id_flxnet(ipass) > 0 ) then
          used = send_data (id_flxnet(ipass),    &
                           flxnet,   &
                           Time_diag, is, js, 1, rmask=mask )
        endif

!------- lw tendency -----------
        if (id_tdt_lw(ipass) > 0 ) then
          used = send_data (id_tdt_lw(ipass), tdtlw,    &
                            Time_diag, is, js, 1, rmask=mask )
        endif

!------- outgoing lw flux toa (olr) -------
        if (id_olr(ipass) > 0 ) then
          used = send_data (id_olr(ipass), olr,    &
                            Time_diag, is, js )
        endif

!------- net radiation (lw + sw) at toa -------
        if (id_netrad_toa(ipass) > 0 ) then
          used = send_data (id_netrad_toa(ipass),   &
                            swin - swout - olr, &
                            Time_diag, is, js )
        endif

!------- net radiation (lw + sw) at 1 Pa-------
        if (id_netrad_1_Pa(ipass) > 0 ) then
          used = send_data (id_netrad_1_Pa(ipass),   &
               swdn_trop(:,:,4) -swup_trop(:,:,4) -netlw_trop(:,:,4), &
                             Time_diag, is, js )
        endif

!------- net lw flux trop (netlw_trop) -------
        if (id_netlw_special(1,ipass) > 0 ) then
          used = send_data (id_netlw_special(1, ipass),   &
                            netlw_trop(:,:,1),  &
                            Time_diag, is, js )
        endif
!------- net lw flux trop (netlw_trop) -------
        if (id_netlw_special(2,ipass) > 0 ) then
          used = send_data (id_netlw_special(2, ipass),   &
                            netlw_trop(:,:,2),  &
                            Time_diag, is, js )
        endif
!------- net lw flux trop (netlw_trop) -------
        if (id_netlw_special(3,ipass) > 0 ) then
          used = send_data (id_netlw_special(3, ipass),   &
                            netlw_trop(:,:,3),  &
                            Time_diag, is, js )
        endif
!------- net lw flux 1 Pa (netlw_trop) -------
        if (id_netlw_special(4,ipass) > 0 ) then
          used = send_data (id_netlw_special(4, ipass),   &
                            netlw_trop(:,:,4),  &
                            Time_diag, is, js )
        endif

!------- upward lw flux surface -------
        if ( id_lwup_sfc(ipass) > 0 ) then
          used = send_data (id_lwup_sfc(ipass), lwups,    &
                            Time_diag, is, js )
        endif

!------- downward lw flux surface -------
        if (id_lwdn_sfc(ipass) > 0 ) then
          used = send_data (id_lwdn_sfc(ipass), lwdns,    &
                            Time_diag, is, js )
        endif

!------- net lw flux surface -------
         if ( id_lwsfc(ipass) > 0 ) then
           used = send_data (id_lwsfc(ipass), lwups-lwdns,    &
                             Time_diag, is, js )
         endif
 
     if (do_lwaerosol_forcing) then

!------- outgoing lw flux toa (olr) with aerosols-------
        if (id_olr_ad(ipass) > 0 ) then
          used = send_data (id_olr_ad(ipass), olr_ad,    &
                             Time_diag, is, js )
        endif

!------- net lw flux surface -------
        if ( id_lwsfc_ad(ipass) > 0 ) then
           used = send_data (id_lwsfc_ad(ipass), lwups_ad-lwdns_ad,    &
                             Time_diag, is, js )
        endif
     endif

!----------------------------------------------------------------------
!    now pass clear-sky diagnostics, if they have been calculated.
!----------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 1

!---- net lw flux ---------
        if (id_flxnet(ipass) > 0 ) then
          used = send_data (id_flxnet(ipass),    &
                            flxnetcf,   &
                            Time_diag, is, js, 1, rmask=mask )
        endif

!------- lw tendency -----------
          if (id_tdt_lw(ipass) > 0 ) then
            used = send_data (id_tdt_lw(ipass), tdtlw_clr,    &
                              Time_diag, is, js, 1, rmask=mask )
          endif

!------- outgoing lw flux toa (olr) -------
          if (id_olr(ipass) > 0 ) then
            used = send_data (id_olr(ipass), olr_clr,   &
                              Time_diag, is, js )
          endif

!------- net radiation (lw + sw) toa -------
          if (id_netrad_toa(ipass) > 0 ) then
            used = send_data (id_netrad_toa(ipass),   &
                              swin_clr - swout_clr - olr_clr,   &
                              Time_diag, is, js )
          endif

!------- net lw flux trop (netlw_trop) -------
          if (id_netlw_special(1,ipass) > 0 ) then
            used = send_data (id_netlw_special(1, ipass),    &
                              netlw_trop_clr(:,:,1),    &
                              Time_diag, is, js )
          endif

!------- net lw flux trop (netlw_trop) -------
          if (id_netlw_special(2,ipass) > 0 ) then
            used = send_data (id_netlw_special(2, ipass),    &
                              netlw_trop_clr(:,:,2),    &
                              Time_diag, is, js )
          endif
!------- net lw flux trop (netlw_trop) -------
          if (id_netlw_special(3,ipass) > 0 ) then
            used = send_data (id_netlw_special(3, ipass),    &
                              netlw_trop_clr(:,:,3),    &
                              Time_diag, is, js )
          endif

!------- upward lw flux surface -------
          if (id_lwup_sfc(ipass) > 0 ) then
            used = send_data (id_lwup_sfc(ipass), lwups_clr,   &
                              Time_diag, is, js )
          endif

!------- downward lw flux surface -------
          if (id_lwdn_sfc(ipass) > 0 ) then
            used = send_data (id_lwdn_sfc(ipass), lwdns_clr,   &
                              Time_diag, is, js )
          endif

!------- net lw flux surface -------
         if ( id_lwsfc(ipass) > 0 ) then
           used = send_data (id_lwsfc(ipass), lwups_clr-lwdns_clr,    &
                             Time_diag, is, js )
        endif   
   
     if (do_lwaerosol_forcing) then

!------- outgoing lw flux toa (olr) with aerosols-------
         if (id_olr_ad(ipass) > 0 ) then
          used = send_data (id_olr_ad(ipass), olr_ad_clr,    &
                            Time_diag, is, js )
         endif   
   
!------- net lw flux surface -------
         if ( id_lwsfc_ad(ipass) > 0 ) then
           used = send_data (id_lwsfc_ad(ipass), lwups_ad_clr-lwdns_ad_clr,    &
                           Time_diag, is, js )
        endif
      endif

        endif  ! (do_clear_sky_pass)
        endif
      endif  ! (do_lw_rad .or. all_step_diagnostics)

!--------------------------------------------------------------------
!    now define various diagnostic integrals.
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!    accumulate global integral quantities 
!--------------------------------------------------------------------
        olr_intgl(is:ie,js:je) = olr(:,:)
        swabs_intgl(is:ie,js:je) = swin(:,:) - swout(:,:)
!       call sum_diag_integral_field ('olr_clr',    olr_clr, is, js)
!       call sum_diag_integral_field ('abs_sw_clr',    &
!                                           swin_clr-swout_clr, is, js)

!--------------------------------------------------------------------
!    accumulate hemispheric integral quantities, if desired. 
!--------------------------------------------------------------------
        if (calc_hemi_integrals) then
          do j=js,je        
            jind = j - js + 1
            iind = 1  ! are assuming all i points are at same latitude

!---------------------------------------------------------------------
!    calculate southern hemisphere integrals.
!---------------------------------------------------------------------
            if (lat(iind,jind) <= 0.0) then
              call sum_diag_integral_field ('sntop_tot_sh ',   &
                                            swin-swout, is, ie, j, j)
              call sum_diag_integral_field ('lwtop_tot_sh ', olr,     &
                                            is, ie, j, j)
              call sum_diag_integral_field ('sngrd_tot_sh ',   &
                                            swdns-swups, is, ie, j, j)
               call sum_diag_integral_field ('lwgrd_tot_sh ',   &
                                          Lw_output(1)%flxnet(:,:,kmax+1),&
                                          is, ie,  j, j)
              if (do_clear_sky_pass) then
                call sum_diag_integral_field ('sntop_clr_sh ',   &
                                              swin_clr-swout_clr, &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwtop_clr_sh ', olr_clr,&
                                              is, ie, j, j)
                call sum_diag_integral_field ('sngrd_clr_sh ',   &
                                              swdns_clr-swups_clr, &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwgrd_clr_sh ',    &
                                       Lw_output(1)%flxnetcf(:,:,kmax+1),&
                                       is, ie, j, j)
              endif

!---------------------------------------------------------------------
!    calculate northern hemisphere integrals.
!---------------------------------------------------------------------
            else
              call sum_diag_integral_field ('sntop_tot_nh ',    &
                                            swin-swout, is, ie, j, j)
              call sum_diag_integral_field ('lwtop_tot_nh ', olr,     &
                                            is, ie, j, j)
              call sum_diag_integral_field ('sngrd_tot_nh ',   &
                                            swdns-swups, is, ie, j, j)
              call sum_diag_integral_field ('lwgrd_tot_nh ',   &
                                          Lw_output(1)%flxnet(:,:,kmax+1),&
                                          is, ie, j, j)        
              if (do_clear_sky_pass) then
                call sum_diag_integral_field ('sntop_clr_nh ',   &
                                              swin_clr-swout_clr,  &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwtop_clr_nh ', olr_clr,&
                                              is, ie, j, j)
                call sum_diag_integral_field ('sngrd_clr_nh ',   &
                                              swdns_clr-swups_clr, &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwgrd_clr_nh ',   &
                                       Lw_output(1)%flxnetcf(:,:,kmax+1),&
                                       is, ie, j, j)
              endif
            endif
          end do

!--------------------------------------------------------------------
!    accumulate global integral quantities 
!--------------------------------------------------------------------
          call sum_diag_integral_field ('sntop_tot_gl ', swin-swout,  &
                                        is, js)
          call sum_diag_integral_field ('lwtop_tot_gl ', olr, is, js)
          call sum_diag_integral_field ('sngrd_tot_gl ', swdns-swups, &
                                        is, js)
          call sum_diag_integral_field ('lwgrd_tot_gl ',  &
                                  Lw_output(1)%flxnet(:,:,kmax+1), is, js)
          if (do_clear_sky_pass) then
            call sum_diag_integral_field ('sntop_clr_gl ',   &
                                          swin_clr-swout_clr, is, js)
            call sum_diag_integral_field ('lwtop_clr_gl ', olr_clr,   &
                                          is, js)
            call sum_diag_integral_field ('sngrd_clr_gl ',   &
                                          swdns_clr-swups_clr, is, js)
            call sum_diag_integral_field ('lwgrd_clr_gl ',   &
                                        Lw_output(1)%flxnetcf(:,:,kmax+1),&
                                        is, js)
          endif
        endif   ! (calc_hemi_integrals)

!---------------------------------------------------------------------



end subroutine produce_radiation_diagnostics



!###################################################################
! <SUBROUTINE NAME="deallocate_arrays">
!  <OVERVIEW>
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call deallocate_arrays (Cldrad_props, Astro, Astro2, Lw_output, &
!                              Fsrad_output, Sw_output)
!  </TEMPLATE>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   Cloud radiative properties
!  </INOUT>
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!   astronomical data for the radiation package
!  </INOUT>
!  <INOUT NAME="Astro2" TYPE="astronomy_type">
!   astronomical data for the radiation package
!  </INOUT>
!  <INOUT NAME="Fsrad_output" TYPE="rad_output_type">
!   radiation output data from the 
!                        original_fms_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </INOUT>
! </SUBROUTINE>
!
subroutine deallocate_arrays (Cldrad_props, Astro, Astro2,  &
                              Aerosol_props, Lw_output, &
                              Fsrad_output, Sw_output, Aerosol_diags)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!---------------------------------------------------------------------

type(cldrad_properties_type), intent(inout)   :: Cldrad_props
type(astronomy_type)        , intent(inout)   :: Astro, Astro2
type(aerosol_properties_type), intent(inout)  :: Aerosol_props
type(lw_output_type),dimension(:), intent(inout)   :: Lw_output
type(fsrad_output_type)     , intent(inout)   :: Fsrad_output
type(sw_output_type),dimension(:), intent(inout)   :: Sw_output
type(aerosol_diagnostics_type), intent(inout)  :: Aerosol_diags

      integer  ::  n

!--------------------------------------------------------------------
!    deallocate the variables in Aerosol_props.
!--------------------------------------------------------------------
      if ( do_rad .and. Rad_control%do_aerosol) then 
        if (Rad_control%volcanic_sw_aerosols) then
          deallocate (Aerosol_props%sw_ext)
          deallocate (Aerosol_props%sw_ssa)
          deallocate (Aerosol_props%sw_asy)
        endif
        if (Rad_control%volcanic_lw_aerosols) then
          deallocate (Aerosol_props%lw_ext)
          deallocate (Aerosol_props%lw_ssa)
          deallocate (Aerosol_props%lw_asy)
        endif
        deallocate (Aerosol_props%ivol)
        if (Sw_control%do_swaerosol .or. &
                      Rad_control%do_swaerosol_forcing) then
          deallocate (Aerosol_props%aerextband)
          deallocate (Aerosol_props%aerssalbband)
          deallocate (Aerosol_props%aerasymmband)
        endif
        if (Lw_control%do_lwaerosol .or. &
                      Rad_control%do_lwaerosol_forcing) then
          deallocate (Aerosol_props%aerextbandlw)
          deallocate (Aerosol_props%aerssalbbandlw)
          deallocate (Aerosol_props%aerextbandlw_cn)
          deallocate (Aerosol_props%aerssalbbandlw_cn)
        endif
        deallocate (Aerosol_props%sulfate_index)
        deallocate (Aerosol_props%optical_index)
        deallocate (Aerosol_props%omphilic_index)
        deallocate (Aerosol_props%bcphilic_index)
        deallocate (Aerosol_props%seasalt1_index)
        deallocate (Aerosol_props%seasalt2_index)
        deallocate (Aerosol_props%seasalt3_index)
        deallocate (Aerosol_props%seasalt4_index)
        deallocate (Aerosol_props%seasalt5_index)
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Astro and Astro2.
!--------------------------------------------------------------------
      if ( do_rad .or. renormalize_sw_fluxes ) then 
        deallocate (Astro%solar)
        deallocate (Astro%cosz )
        deallocate (Astro%fracday)
        deallocate (Astro%solar_p)
        deallocate (Astro%cosz_p )
        deallocate (Astro%fracday_p)
          if ( do_sw_rad .and. renormalize_sw_fluxes   &
              .and. Sw_control%do_diurnal ) then 
            deallocate (Astro2%solar)
            deallocate (Astro2%cosz )
            deallocate (Astro2%fracday)
        endif
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Lw_output.
!--------------------------------------------------------------------
      if (do_sea_esf_rad) then
        if (do_lw_rad) then
         do n=1,size_of_lwoutput
          deallocate (Lw_output(n)%heatra    )
          deallocate (Lw_output(n)%flxnet    )
          deallocate (Lw_output(n)%netlw_special)
          deallocate (Lw_output(n)%bdy_flx)
          if (Rad_control%do_totcld_forcing) then
            deallocate (Lw_output(n)%heatracf  )
            deallocate (Lw_output(n)%flxnetcf  )
            deallocate (Lw_output(n)%netlw_special_clr)
            deallocate (Lw_output(n)%bdy_flx_clr)
          endif
        end do
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Sw_output.
!--------------------------------------------------------------------
        if (do_sw_rad) then
        do n=1,size_of_swoutput
          deallocate (Sw_output(n)%dfsw     )
          deallocate (Sw_output(n)%ufsw     )
          deallocate (Sw_output(n)%dfsw_dir_sfc )
          deallocate (Sw_output(n)%ufsw_dir_sfc )
          deallocate (Sw_output(n)%dfsw_dif_sfc )
          deallocate (Sw_output(n)%ufsw_dif_sfc )
          deallocate (Sw_output(n)%fsw     )
          deallocate (Sw_output(n)%hsw     )
          deallocate (Sw_output(n)%dfsw_vis_sfc    )
          deallocate (Sw_output(n)%ufsw_vis_sfc    )
          deallocate (Sw_output(n)%ufsw_vis_sfc_dir)
          deallocate (Sw_output(n)%dfsw_vis_sfc_dir    )
          deallocate (Sw_output(n)%dfsw_vis_sfc_dif    )
          deallocate (Sw_output(n)%ufsw_vis_sfc_dif    )
          deallocate (Sw_output(n)%swdn_special)
          deallocate (Sw_output(n)%swup_special)
          deallocate (Sw_output(n)%bdy_flx)
          if (Rad_control%do_totcld_forcing) then
            deallocate (Sw_output(n)%dfswcf   )
            deallocate (Sw_output(n)%ufswcf   )
            deallocate (Sw_output(n)%fswcf   )
            deallocate (Sw_output(n)%hswcf   )
            deallocate (Sw_output(n)%dfsw_dir_sfc_clr)
            deallocate (Sw_output(n)%dfsw_dif_sfc_clr)
            deallocate (Sw_output(n)%dfsw_vis_sfc_clr)
            deallocate (Sw_output(n)%swdn_special_clr)
            deallocate (Sw_output(n)%swup_special_clr)
            deallocate (Sw_output(n)%bdy_flx_clr)
          endif
       end do
        endif
      endif

!--------------------------------------------------------------------
!    call cldrad_props_dealloc to deallocate the variables in 
!    Cldrad_props. 
!--------------------------------------------------------------------
      if (do_rad .and. do_sea_esf_rad) then
        call cldrad_props_dealloc (Cldrad_props)
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Fsrad_output.
!--------------------------------------------------------------------

      if (.not. do_sea_esf_rad .and. do_rad) then
         deallocate (Fsrad_output%tdtsw   )
         deallocate (Fsrad_output%tdtlw   )
         deallocate (Fsrad_output%swdns   )
         deallocate (Fsrad_output%swups   )
         deallocate (Fsrad_output%lwdns   )
         deallocate (Fsrad_output%lwups   )
         deallocate (Fsrad_output%swin    )
         deallocate (Fsrad_output%swout   )
         deallocate (Fsrad_output%olr     )
        if (do_clear_sky_pass) then
           deallocate (Fsrad_output%tdtsw_clr )
           deallocate (Fsrad_output%tdtlw_clr )
           deallocate (Fsrad_output%swdns_clr )
           deallocate (Fsrad_output%swups_clr )
           deallocate (Fsrad_output%lwdns_clr )
           deallocate (Fsrad_output%lwups_clr )
           deallocate (Fsrad_output%swin_clr  )
           deallocate (Fsrad_output%swout_clr )
           deallocate (Fsrad_output%olr_clr   )
        endif 
      endif 

!--------------------------------------------------------------------
!    deallocate the window-resident variables in Aerosol_props.
!--------------------------------------------------------------------
      if (do_rad .and. Rad_control%do_aerosol) then
        deallocate (Aerosol_diags%extopdep)
        deallocate (Aerosol_diags%absopdep)
        deallocate (Aerosol_diags%asymdep)
        deallocate (Aerosol_diags%extopdep_vlcno)
        deallocate (Aerosol_diags%absopdep_vlcno)
        deallocate (Aerosol_diags%sw_heating_vlcno)
        deallocate (Aerosol_diags%lw_extopdep_vlcno)
        deallocate (Aerosol_diags%lw_absopdep_vlcno)
      endif

!---------------------------------------------------------------------



end subroutine deallocate_arrays 



!#####################################################################
! <SUBROUTINE NAME="calculate_auxiliary_variables">
!  <OVERVIEW>
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!  </OVERVIEW>
!  <DESCRIPTION>
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call calculate_auxiliary_variables (Atmos_input)
!  </TEMPLATE>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable, its press and temp
!                   components are input, and its deltaz, rel_hum, 
!                   pflux, tflux and aerosolrelhum components are 
!                   calculated here and output.
!  </INOUT>
! </SUBROUTINE>
!
subroutine calculate_auxiliary_variables (Atmos_input)

!----------------------------------------------------------------------
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!---------------------------------------------------------------------

type(atmos_input_type), intent(inout)  :: Atmos_input

!--------------------------------------------------------------------
!   intent(inout) variables
!
!      Atmos_input  atmos_input_type variable, its press and temp
!                   components are input, and its deltaz, rel_hum, 
!                   pflux, tflux and aerosolrelhum components are 
!                   calculated here and output.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables

      real, dimension (size(Atmos_input%temp, 1), &
                       size(Atmos_input%temp, 2), &
                       size(Atmos_input%temp, 3) - 1) :: &
                                                     qsat, qv, tv
      integer   ::  k
      integer   ::  kmax


!--------------------------------------------------------------------
!    define flux level pressures (pflux) as midway between data level
!    (layer-mean) pressures. specify temperatures at flux levels
!    (tflux).
!--------------------------------------------------------------------
      do k=ks+1,ke
        Atmos_input%pflux(:,:,k) = 0.5E+00*  &
                (Atmos_input%press(:,:,k-1) + Atmos_input%press(:,:,k))
        Atmos_input%tflux(:,:,k) = 0.5E+00*  &
                (Atmos_input%temp (:,:,k-1) + Atmos_input%temp (:,:,k))
      end do
      Atmos_input%pflux(:,:,ks  ) = 0.0E+00
      Atmos_input%pflux(:,:,ke+1) = Atmos_input%press(:,:,ke+1)
      Atmos_input%tflux(:,:,ks  ) = Atmos_input%temp (:,:,ks  )
      Atmos_input%tflux(:,:,ke+1) = Atmos_input%temp (:,:,ke+1)

!-------------------------------------------------------------------
!    define deltaz in meters.
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%temp(:,:,ks:ke)*    &
                  (1.0 + D608*Atmos_input%rh2o(:,:,:))
      Atmos_input%deltaz(:,:,ks) = log_p_at_top*RDGAS*tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%deltaz(:,:,k) = alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      end do

!-------------------------------------------------------------------
!    define deltaz in meters to be used in cloud feedback analysis.
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%cloudtemp(:,:,ks:ke)*    &
                  (1.0 + D608*Atmos_input%cloudvapor(:,:,:))
      Atmos_input%clouddeltaz(:,:,ks) = log_p_at_top*RDGAS*  &
                                        tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%clouddeltaz(:,:,k) =    &
                            alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      end do

!------------------------------------------------------------------
!    define the relative humidity.
!------------------------------------------------------------------
      kmax = size(Atmos_input%temp,3) - 1
      qv(:,:,1:kmax) = Atmos_input%rh2o(:,:,1:kmax) /    &
                                   (1.0 + Atmos_input%rh2o(:,:,1:kmax))
      call compute_qs (Atmos_input%temp(:,:,1:kmax),  &
                       Atmos_input%press(:,:,1:kmax),  &
                       qsat(:,:,1:kmax), q = qv(:,:,1:kmax))
      do k=1,kmax
        Atmos_input%rel_hum(:,:,k) = qv(:,:,k) / qsat(:,:,k)
        Atmos_input%rel_hum(:,:,k) =    &
                                  MIN (Atmos_input%rel_hum(:,:,k), 1.0)
      end do

!------------------------------------------------------------------
!    define the relative humidity seen by the aerosol code.
!------------------------------------------------------------------
        qv(:,:,1:kmax) = Atmos_input%aerosolvapor(:,:,1:kmax) /    &
                         (1.0 + Atmos_input%aerosolvapor(:,:,1:kmax))
        call compute_qs (Atmos_input%aerosoltemp(:,:,1:kmax), &
                         Atmos_input%aerosolpress(:,:,1:kmax),  &
                         qsat(:,:,1:kmax), q = qv(:,:,1:kmax))
      do k=1,kmax
         Atmos_input%aerosolrelhum(:,:,k) = qv(:,:,k) / qsat(:,:,k)
         Atmos_input%aerosolrelhum(:,:,k) =    &
                            MIN (Atmos_input%aerosolrelhum(:,:,k), 1.0)
      end do
 
!----------------------------------------------------------------------


end subroutine calculate_auxiliary_variables


!#######################################################################


                 end module radiation_driver_mod
