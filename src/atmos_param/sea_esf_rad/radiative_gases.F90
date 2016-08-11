                  module radiative_gases_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  Module that defines mixing ratios of radiatively-active 
!    gases to be used in the calculation of longwave and shortwave
!    radiative fluxes and heating rates in the sea_esf_rad radiation
!    package.
! </OVERVIEW>
! <DESCRIPTION>
!  Module that defines mixing ratios of radiatively-active 
!    gases to be used in the calculation of longwave and shortwave
!    radiative fluxes and heating rates in the sea_esf_rad radiation
!    package.
! </DESCRIPTION>

!  shared modules:

use time_manager_mod,    only: time_manager_init, time_type, set_date, &
                               get_calendar_type, GREGORIAN, &
                               operator(>=), operator(-), operator(<=),&
                               operator(>),  operator (<), get_date,  &
                               set_time, operator(+), print_date,  &
                               days_in_year, get_time, length_of_year
use diag_manager_mod,    only: diag_manager_init, get_base_time
use mpp_mod,             only: input_nml_file
use fms_mod,             only: open_namelist_file, fms_init, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               file_exist, write_version_number, &
                               check_nml_error, error_mesg, &
                               FATAL, NOTE, close_file, &
                               open_restart_file, read_data
use fms_io_mod,          only: get_restart_io_mode, &
                               register_restart_field, restart_file_type, &
                               save_restart, restore_state, query_initialized
use time_interp_mod,     only: time_interp_init, time_interp
use tracer_manager_mod,  only: get_tracer_index, NO_TRACER
use field_manager_mod,   only: MODEL_ATMOS

!  shared radiation package modules:

use rad_utilities_mod,   only: rad_utilities_init, Lw_control, &
                               atmos_input_type, &
                               radiative_gases_type, Rad_control

! component modules:

use ozone_mod,           only: ozone_driver, ozone_init,  &
                               ozone_time_vary, ozone_endts, ozone_end

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    radiative_gases_mod defines mixing ratios of radiatively-active 
!    gases to be used in the calculation of longwave and shortwave
!    radiative fluxes and heating rates in the sea_esf_rad radiation
!    package.
!---------------------------------------------------------------------
 

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  &
'$Id: radiative_gases.F90,v 20.0 2013/12/13 23:20:28 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

public     &
         radiative_gases_init, define_radiative_gases,  &
         radiative_gases_time_vary, radiative_gases_endts,  &
         radiative_gases_end, radiative_gases_dealloc,  &
         radiative_gases_restart

private    &
! called from radiative_gases_init:
         validate_time_varying_inputs, &
         define_ch4, define_n2o, define_f11, &
         define_f12, define_f113, define_f22, &
         define_co2, read_gas_timeseries,  &

! called from define_radiative_gases:
         define_gas_amount,     &

! called from radiative_gases_end:
         write_restart_radiative_gases


!---------------------------------------------------------------------
!-------- namelist  ---------

integer          :: verbose = 0   ! verbosity levels, values run from
                                  ! 0 (least output) to 5 (most output)
integer          :: gas_printout_freq = 32*24
                                  ! frequency of outputting gas concen-
                                  ! trations [ hours ] default: 32 days

!--------------------------------------------------------------------
!    time_varying_xxx  logical flag indicating whether the vol. mixing 
!                      ratio of gas xxx varies in time
!    xxx_data_source   source of input data for initial value of gas xxx
!    xxx_floor         smallest value allowed for gas xxx vol. mixing 
!                      ratio [ no. / no. ]
!    xxx_ceiling       largest value allowed for gas xxx vol. mixing 
!                      ratio [ no. / no. ]
!    xxx_specification_type
!                      indicator as to the form of time variation of
!                      vol. mixing ratio for gas xxx; either 
!                      'base_and_trend' or 'time_series'
!    xxx_variation_type
!                      indicator as to the form of time variation of
!                      the vol. mixing ratio of gas xxx; either 'linear'
!                      or 'logarithmic'. Must be 'linear' for 
!                      'time_series'.
!      
!         The following variables only have relevance when 
!         specification_type = 'base_and_trend':
!
!    xxx_base_value    initial value of gas xxx vol. mixing ratio when
!                      xxx_data_source is 'namelist' [ no. / no. ]
!    xxx_base_time     time at which xxx_base_value is relevant spec-
!                      ified as (year, month, day, 0, 0, 0). (Can only 
!                      be specified as 00Z on the particular day).
!
!    xxx_change_rate   time rate of change of gas xxx vol. mixing ratio.
!                      [  1 +/- % per year ]
!
!--------------------------------------------------------------------

logical              :: time_varying_co2 = .false.
character(len=16)    :: co2_data_source  = '   '
real                 :: co2_base_value    = 0.0
integer,dimension(6) :: co2_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: co2_change_rate   = 0.0
real                 :: co2_floor = 0.0
real                 :: co2_ceiling = 1.0E6
character(len=16)    :: co2_specification_type = '              '
character(len=16)    :: co2_variation_type = '           '

logical              :: time_varying_ch4 = .false.
character(len=16)    :: ch4_data_source  = '   '
real                 :: ch4_base_value    = 0.0
integer,dimension(6) :: ch4_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: ch4_change_rate   = 0.0
real                 :: ch4_floor = 0.0
real                 :: ch4_ceiling = 1.0E6
character(len=16)    :: ch4_specification_type = '              '
character(len=16)    :: ch4_variation_type = '           '

logical              :: time_varying_n2o = .false.
character(len=16)    :: n2o_data_source  = '   '
real                 :: n2o_base_value    = 0.0
integer,dimension(6) :: n2o_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: n2o_change_rate   = 0.0
real                 :: n2o_floor = 0.0
real                 :: n2o_ceiling = 1.0E6
character(len=16)    :: n2o_specification_type = '              '
character(len=16)    :: n2o_variation_type = '           '

logical              :: time_varying_f11 = .false.
character(len=16)    :: f11_data_source  = '   '
real                 :: f11_base_value    = 0.0
integer,dimension(6) :: f11_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: f11_change_rate   = 0.0
real                 :: f11_floor = 0.0
real                 :: f11_ceiling = 1.0E6
character(len=16)    :: f11_specification_type = '              '
character(len=16)    :: f11_variation_type = '           '

logical              :: time_varying_f12 = .false.
character(len=16)    :: f12_data_source  = '   '
real                 :: f12_base_value    = 0.0
integer,dimension(6) :: f12_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: f12_change_rate   = 0.0
real                 :: f12_floor = 0.0
real                 :: f12_ceiling = 1.0E6
character(len=16)    :: f12_specification_type = '              '
character(len=16)    :: f12_variation_type = '           '

logical              :: time_varying_f113 = .false.
character(len=16)    :: f113_data_source  = '   '
real                 :: f113_base_value    = 0.0
integer,dimension(6) :: f113_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: f113_change_rate   = 0.0
real                 :: f113_floor = 0.0
real                 :: f113_ceiling = 1.0E6
character(len=16)    :: f113_specification_type = '              '
character(len=16)    :: f113_variation_type = '           '

logical              :: time_varying_f22 = .false.
character(len=16)    :: f22_data_source  = '   '
real                 :: f22_base_value    = 0.0
integer,dimension(6) :: f22_base_time     = (/ 0,0,0,0,0,0 /)
real                 :: f22_change_rate   = 0.0
real                 :: f22_floor = 0.0
real                 :: f22_ceiling = 1.0E6
character(len=16)    :: f22_specification_type = '              '
character(len=16)    :: f22_variation_type = '           '

integer, dimension(6) ::       &
                         co2_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in co2  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
integer, dimension(6) ::       &
                         ch4_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in ch4  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
integer, dimension(6) ::       &
                         n2o_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in n2o  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
integer, dimension(6) ::       &
                         f11_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in f11  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
integer, dimension(6) ::       &
                         f12_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in f12  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
integer, dimension(6) ::       &
                         f113_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in f113  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
integer, dimension(6) ::       &
                         f22_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in f22  data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
logical              :: time_varying_restart_bug = .false.
logical              :: use_globally_uniform_co2 = .true.
namelist /radiative_gases_nml/ verbose, &
        use_globally_uniform_co2, &
        gas_printout_freq, time_varying_restart_bug, &
        co2_dataset_entry, ch4_dataset_entry, n2o_dataset_entry,  &
        f11_dataset_entry, f12_dataset_entry, f113_dataset_entry, &
        f22_dataset_entry, &

        time_varying_co2, co2_data_source, co2_base_value,  &
        co2_base_time, co2_change_rate, co2_floor, co2_ceiling,  &
        co2_specification_type, co2_variation_type, &

                                time_varying_ch4, ch4_data_source, &
                                ch4_base_value, ch4_base_time, &
                                ch4_change_rate, ch4_floor,   &
                                ch4_ceiling, ch4_specification_type, &
                                ch4_variation_type, &

                                time_varying_n2o, n2o_data_source, &
                                n2o_base_value, n2o_base_time, &
                                n2o_change_rate, n2o_floor,   &
                                n2o_ceiling, n2o_specification_type, &
                                n2o_variation_type, &

                                time_varying_f11, f11_data_source, &
                                f11_base_value, f11_base_time, &
                                f11_change_rate, f11_floor,   &
                                f11_ceiling, f11_specification_type, &
                                f11_variation_type, &

                                time_varying_f12, f12_data_source, &
                                f12_base_value, f12_base_time, &
                                f12_change_rate, f12_floor,   &
                                f12_ceiling, f12_specification_type, &
                                f12_variation_type, &

                                time_varying_f113, f113_data_source, &
                                f113_base_value, f113_base_time, &
                                f113_change_rate, f113_floor,   &
                                f113_ceiling, f113_specification_type, &
                                f113_variation_type, &

                                time_varying_f22, f22_data_source, &
                                f22_base_value, f22_base_time, &
                                f22_change_rate, f22_floor,   &
                                f22_ceiling, f22_specification_type, &
                                f22_variation_type

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

!--- for netcdf restart
type(restart_file_type), save :: Rad_restart
logical                       :: do_netcdf_restart= .true.
integer                       ::  vers   ! version number of restart file 
!---------------------------------------------------------------------
!    list of restart versions of radiation_driver.res readable by this 
!    module.
!---------------------------------------------------------------------
integer, dimension(3)    ::  restart_versions = (/ 1, 2, 3 /)

!--------------------------------------------------------------------
!    initial mixing ratios of the various radiative gases. if a gas 
!    is not active, its mixing ratio is set to zero.
!--------------------------------------------------------------------
real         ::  rch4, rn2o, rf11, rf12, rf113, rf22, rco2

!--------------------------------------------------------------------
!    this is the mixing ratio of gas used when the transmission 
!    functions were last calculated.
!--------------------------------------------------------------------
real         ::  co2_for_last_tf_calc
real         ::  ch4_for_last_tf_calc
real         ::  n2o_for_last_tf_calc

!RSH
!  Need these as module variables rather than components of derived type
!   since they are needed in region executed by master thread only:
real         ::  co2_for_next_tf_calc
real         ::  ch4_for_next_tf_calc
real         ::  n2o_for_next_tf_calc
real         ::  ch4_tf_offset
real         ::  n2o_tf_offset
real         ::  co2_tf_offset
!--------------------------------------------------------------------
!    these variables are .true. if transmission functions are cal-
!    culated for the referenced gas.
!--------------------------------------------------------------------
logical, parameter   :: co2_uses_tfs  = .true.
logical, parameter   :: ch4_uses_tfs  = .true.
logical, parameter   :: n2o_uses_tfs  = .true.
logical, parameter   :: f11_uses_tfs  = .false.
logical, parameter   :: f12_uses_tfs  = .false.
logical, parameter   :: f113_uses_tfs = .false.
logical, parameter   :: f22_uses_tfs  = .false.

!--------------------------------------------------------------------
!    these variables contain the mixing ratios of the radiative gases
!    at the current time. if the gases are fixed, this is the same as
!    the initial value; otherwise it will be time-varying.
!--------------------------------------------------------------------
real   :: rrvco2, rrvf11, rrvf12, rrvf113, rrvf22, rrvch4, rrvn2o 

!---------------------------------------------------------------------
!    variables to specify data for gas xxx, when xxx_specification_type
!    is 'time_series':
!
!    xxx_time_list          list of times (time_type variable) for 
!                           which values for gas xxx are specified 
!                           for 'time_series', these define the data 
!                           points, for 'base_and_trend', the single
!                           entry is the xxx_base_time  [ time_type ]
!    xxx_value              values (vol. mixing ratio) of gas xxx at
!                           times given by xxx_time_list. data may have
!                           to be converted from other units (eg, ppmv 
!                           or ppbv or pptv). [real]
!---------------------------------------------------------------------
type(time_type), dimension(:), pointer :: Co2_time_list, N2o_time_list,&
                                          Ch4_time_list, F11_time_list,&
                                          F12_time_list, F22_time_list,&
                                          F113_time_list
real,            dimension(:), pointer :: co2_value, ch4_value,  &
                                          n2o_value, f11_value,  &
                                          f12_value, f113_value, &
                                          f22_value

!---------------------------------------------------------------------
!    miscellaneous variables:
!
!    restart_present        restart file present ?
!    module_is_initialized  module is initialized ?
!    pts_processed          number of processor's columns that have 
!                           been processed on the current time step 
!    total_points           number of model columns on the processor
!    define_xxx_for_last_tf_calc
!                           the tf's that were used for gas xxx on the 
!                           last step of the previous job must be 
!                           recalculated ? this can only be true when
!                           a restart file version earlier than version
!                           3 is being read. 

!---------------------------------------------------------------------
logical      ::  restart_present =  .false.   
logical      ::  module_is_initialized =  .false. 
integer      ::  ico2
logical      ::  define_co2_for_last_tf_calc = .false.
logical      ::  define_ch4_for_last_tf_calc = .false.
logical      ::  define_n2o_for_last_tf_calc = .false.
logical      ::  printed_current_floor_msg = .false.
logical      ::  printed_current_ceiling_msg = .false.
logical      ::  printed_next_floor_msg = .false.
logical      ::  printed_next_ceiling_msg = .false.
integer      ::  print_alarm = 0
 
type(time_type) :: Model_init_time  ! initial calendar time for model  
                                    ! [ time_type ]

!-------------------------------------------------------------------
!   xxx_offset  ! difference between model initial time and gas time-
                ! series mapped to model initial time
                ! [ time_type ]
!   xxx_entry   ! time in gas timeseries which is mapped to model 
                ! initial time
                ! [ time_type ]
!   negative_offset_xxx 
                !  the model initial time is later than the gas 
                !  xxx_dataset_entry time  ?
!-------------------------------------------------------------------
type(time_type)    :: Co2_offset,  Co2_entry
type(time_type)    :: Ch4_offset,  Ch4_entry
type(time_type)    :: N2o_offset,  N2o_entry
type(time_type)    :: F11_offset,  F11_entry
type(time_type)    :: F12_offset,  F12_entry
type(time_type)    :: F113_offset,  F113_entry
type(time_type)    :: F22_offset,  F22_entry

logical    :: negative_offset_co2 = .false.
logical    :: negative_offset_ch4 = .false.
logical    :: negative_offset_n2o = .false.
logical    :: negative_offset_f11 = .false.
logical    :: negative_offset_f12 = .false.
logical    :: negative_offset_f113 = .false.
logical    :: negative_offset_f22 = .false.

logical    :: co2_tfs_needed = .true.
logical    :: ch4_tfs_needed = .true.
logical    :: n2o_tfs_needed = .true.

!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!####################################################################
! <SUBROUTINE NAME="radiative_gases_init">
!  <OVERVIEW>
!   Subroutine to initialize radiative_gases module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize radiative_gases module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiative_gases_init (pref, latb, lonb)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference prssure profiles
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   2d array of model longitudes at cell corners [radians]
!  </IN>
! </SUBROUTINE>
!
subroutine radiative_gases_init (pref, latb, lonb)

!---------------------------------------------------------------------
!    radiative_gases_init is the constructor for radiative_gases_mod.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: latb, lonb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [pascals]
!       latb      2d array of model latitudes at cell corners [radians]
!       lonb      2d array of model longitudes at cell corners [radians]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer              :: unit     ! unit number for i/o operation
      integer              :: ierr     ! error code 
      integer              :: io       ! io status upon completion 
      integer              :: calendar ! calendar type used in model
      character(len=8)     :: gas_name ! name associated with current
                                       ! gas being processed
      character(len=32)    :: restart_file
      integer              :: id_restart
      integer              :: logunit  ! unit number for writing to logfile.

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
      call diag_manager_init
      call time_interp_init

!-----------------------------------------------------------------------
!    read namelist.              
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=radiative_gases_nml, iostat=io)
      ierr = check_nml_error(io,'radiative_gases_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=radiative_gases_nml, iostat=io, end=10) 
        ierr = check_nml_error(io,'radiative_gases_nml')
        end do                   
10      call close_file (unit)   
      endif                      
#endif
      call get_restart_io_mode(do_netcdf_restart)

                                  
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                       write (logunit, nml=radiative_gases_nml)

!---------------------------------------------------------------------
!    force xxx_data_source to be 'input' when gas is time-varying and
!    that variation is specified via a time-series file.
!----------------------------------------------------------------------
      if (time_varying_ch4 .and.   &
          trim(ch4_specification_type) == 'time_series' .and. &
          trim (ch4_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when ch4 is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

      if (time_varying_n2o .and.   &
          trim(n2o_specification_type) == 'time_series' .and. &
          trim (n2o_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when n2o is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

      if (time_varying_co2 .and.   &
          trim(co2_specification_type) == 'time_series' .and. &
          trim (co2_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when co2 is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

      if (time_varying_f11 .and.   &
          trim(f11_specification_type) == 'time_series' .and. &
          trim (f11_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when f11 is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

      if (time_varying_f12 .and.   &
          trim(f12_specification_type) == 'time_series' .and. &
          trim (f12_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when f12 is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

      if (time_varying_f113 .and.   &
          trim(f113_specification_type) == 'time_series' .and. &
          trim (f113_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when f113 is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

      if (time_varying_f22 .and.   &
          trim(f22_specification_type) == 'time_series' .and. &
          trim (f22_data_source) /= 'input') then
        call error_mesg ('radiative_gases_mod', &
          ' when f22 is time-varying and comes from a timeseries,&
          & data source must be specified as "input" ', FATAL)  
      endif

!--------------------------------------------------------------------
!    time variation of radiative gases is not currently available if the
!    gregorian calendar is being employed. NOTE : gregorian calendar
!    not available currently in FMS; if it becomes available, then
!    code to properly handle the time variation of gases with that
!    calendar may be developed.
!--------------------------------------------------------------------
      if (time_varying_co2 .or. &
          time_varying_ch4 .or. &
          time_varying_n2o .or. &
          time_varying_f11 .or. &
          time_varying_f12 .or. &
          time_varying_f113 .or. &
          time_varying_f22 )  then
        calendar = get_calendar_type()
        if (calendar == GREGORIAN ) then 
          call error_mesg ('radiative_gases_mod', &
               'code not available to handle time-varying radiative &
                &gases with gregorian calendar', FATAL)
        endif
      endif

!---------------------------------------------------------------------
!    if present, read the radiative gases restart file. set a flag 
!    indicating the presence of the file.
!---------------------------------------------------------------------
     restart_file = 'radiative_gases.res.nc'
      if(do_netcdf_restart) then
         id_restart = register_restart_field(Rad_restart, restart_file, 'vers', vers, no_domain = .true. )
         id_restart = register_restart_field(Rad_restart, restart_file, 'rco2', rco2, no_domain = .true. )          
         id_restart = register_restart_field(Rad_restart, restart_file, 'rf11', rf11, no_domain = .true. )
         id_restart = register_restart_field(Rad_restart, restart_file, 'rf12', rf12, no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'rf113', rf113, no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'rf22', rf22, no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'rch4', rch4, no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'rn2o', rn2o, no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'co2_for_last_tf_calc', &
                                             co2_for_last_tf_calc, mandatory=.false., no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'ch4_for_last_tf_calc', &
                                             ch4_for_last_tf_calc, mandatory=.false., no_domain = .true. ) 
         id_restart = register_restart_field(Rad_restart, restart_file, 'n2o_for_last_tf_calc', &
                                             n2o_for_last_tf_calc, mandatory=.false., no_domain = .true. )     
      endif

      restart_present = .false.
      if (file_exist('INPUT/radiative_gases.res.nc')) then
         if (mpp_pe() == mpp_root_pe()) call error_mesg ('radiative_gases_mod', &
              'Reading NetCDF formatted restart file: INPUT/radiative_gases.res.nc', NOTE)
         if(.not. do_netcdf_restart) call error_mesg ('radiative_gases_mod', &
              'netcdf format restart file INPUT/radiative_gases.res.nc exist, but do_netcdf_restart is false.', FATAL)
         call restore_state(Rad_restart)
         restart_present = .true.
         if(vers >= 3) then
            if(.NOT. query_initialized(Rad_restart, id_restart) ) call error_mesg('radiative_gases_mod', &
                'vers >=3 and INPUT/radiative_gases.res.nc exist, but field n2o_for_last_tf_calc does not in that file', FATAL)
         else
            define_co2_for_last_tf_calc = .true.
            define_ch4_for_last_tf_calc = .true.
            define_n2o_for_last_tf_calc = .true.
         endif       
         vers = restart_versions(size(restart_versions(:)))     
      else
         if (file_exist ('INPUT/radiative_gases.res')) then
           call error_mesg ('radiative_gases_mod', &
                 'Native formatted restart file no longer supported.', FATAL)
         endif
      endif

!---------------------------------------------------------------------
!    call a routine for each gas to initialize its mixing ratio
!    and set a flag indicating whether it is fixed in time or time-
!    varying.  fixed-in-time gases will be defined from the 
!    source specified in the namelist.
!---------------------------------------------------------------------
      call define_ch4 (ch4_data_source)
      call define_n2o (n2o_data_source)
      call define_f11 (f11_data_source)
      call define_f12 (f12_data_source)
      call define_f113(f113_data_source)
      call define_f22 (f22_data_source)
      call define_co2 (co2_data_source)

!---------------------------------------------------------------------
!    define logical variable indicating whether ch4 is active.
!---------------------------------------------------------------------
      if ((.not. time_varying_ch4) .and. rch4 == 0.0) then
        Lw_control%do_ch4 = .false.
      else
        Lw_control%do_ch4 = .true.
      endif

!---------------------------------------------------------------------
!    define logical variable indicating whether n2o is active.
!---------------------------------------------------------------------
      if ((.not. time_varying_n2o) .and. rn2o == 0.0) then
        Lw_control%do_n2o = .false.
      else
        Lw_control%do_n2o = .true.
      endif

!--------------------------------------------------------------------
!    set flag to indicate variable has been initialized.
!--------------------------------------------------------------------
      Lw_control%do_ch4_iz = .true.
      Lw_control%do_n2o_iz = .true.

!---------------------------------------------------------------------
!    if any of the cfcs are activated, set a flag indicating that cfcs
!    are active.
!---------------------------------------------------------------------
      if ((.not. time_varying_f11) .and. rf11 == 0.0 .and. &
          (.not. time_varying_f12) .and. rf12 == 0.0 .and. &
          (.not. time_varying_f113) .and. rf113 == 0.0 .and. &
          (.not. time_varying_f22) .and. rf22 == 0.0 )  then 
        Lw_control%do_cfc = .false.
      else
        Lw_control%do_cfc = .true.
      endif

!---------------------------------------------------------------------
!    set flag to indicate variable has been initialized.
!---------------------------------------------------------------------
      Lw_control%do_cfc_iz = .true.

!---------------------------------------------------------------------
!    define a logical variable indicating whether co2 is to be 
!    activated. currently co2 must be activated. 
!---------------------------------------------------------------------
      if ((.not. time_varying_co2) .and. rco2 == 0.0) then
        Lw_control%do_co2 = .false.
      else
        Lw_control%do_co2 = .true.
      endif

!---------------------------------------------------------------------
!    set flag to indicate variable has been initialized.
!---------------------------------------------------------------------
      Lw_control%do_co2_iz = .true.

!--------------------------------------------------------------------
!    define module variable which will be contain mixing ratio of each 
!    gas as model is integrated in time.
!--------------------------------------------------------------------
      rrvch4  = rch4
      rrvn2o  = rn2o
      rrvf11  = rf11
      rrvf12  = rf12
      rrvf113 = rf113
      rrvf22  = rf22
      rrvco2  = rco2

!--------------------------------------------------------------------
!    verify that the nml variables used when co2 is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_co2) then
        gas_name = 'co2 '
        call validate_time_varying_inputs   &
                        (gas_name, co2_base_time, co2_base_value,  &
                         co2_specification_type, co2_change_rate, &
                         co2_dataset_entry,  negative_offset_co2, &
                         Co2_offset, Co2_entry, &
                         co2_variation_type, Co2_time_list, co2_value)
      endif

!--------------------------------------------------------------------
!    verify that the nml variables used when ch4 is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_ch4) then
        gas_name = 'ch4 '
        call validate_time_varying_inputs   &
                        (gas_name, ch4_base_time, ch4_base_value,  &
                         ch4_specification_type, ch4_change_rate, &
                         ch4_dataset_entry,  negative_offset_ch4, &
                         Ch4_offset, Ch4_entry, &
                         ch4_variation_type, Ch4_time_list, ch4_value)
      endif

!--------------------------------------------------------------------
!    verify that the nml variables used when n2o is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_n2o) then
        gas_name = 'n2o '
        call validate_time_varying_inputs   &
                        (gas_name, n2o_base_time, n2o_base_value,  &
                         n2o_specification_type, n2o_change_rate, &
                         n2o_dataset_entry,  negative_offset_n2o, &
                         N2o_offset, N2o_entry, &
                         n2o_variation_type, N2o_time_list, n2o_value)
      endif

!--------------------------------------------------------------------
!    verify that the nml variables used when f11 is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_f11) then
        gas_name = 'f11 '
        call validate_time_varying_inputs  &
                        (gas_name, f11_base_time, f11_base_value,  &
                         f11_specification_type, f11_change_rate, &
                         f11_dataset_entry,  negative_offset_f11, &
                         F11_offset, F11_entry, &
                         f11_variation_type, F11_time_list, f11_value)
      endif

!--------------------------------------------------------------------
!    verify that the nml variables used when f12 is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_f12) then
        gas_name = 'f12 '
        call validate_time_varying_inputs   &
                        (gas_name, f12_base_time, f12_base_value, &
                         f12_specification_type, f12_change_rate, &
                         f12_dataset_entry,  negative_offset_f12, &
                         F12_offset, F12_entry, &
                         f12_variation_type, F12_time_list, f12_value)
      endif

!--------------------------------------------------------------------
!    verify that the nml variables used when f113 is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_f113) then
        gas_name = 'f113'
        call validate_time_varying_inputs   &
                        (gas_name, f113_base_time, f113_base_value, &
                         f113_specification_type, f113_change_rate, &
                         f113_dataset_entry,  negative_offset_f113, &
                         F113_offset, F113_entry, &
                         f113_variation_type, F113_time_list,   &
                         f113_value)
      endif

!--------------------------------------------------------------------
!    verify that the nml variables used when f22 is varying with time
!    have been given acceptable values.
!--------------------------------------------------------------------
      if (time_varying_f22) then
        gas_name = 'f22 '
        call validate_time_varying_inputs  &
                        (gas_name, f22_base_time, f22_base_value, &
                         f22_specification_type, f22_change_rate, &
                         f22_dataset_entry,  negative_offset_f22, &
                         F22_offset, F22_entry, &
                         f22_variation_type, F22_time_list, f22_value)
      endif

      ico2 = get_tracer_index(MODEL_ATMOS, 'co2')
      if (ico2 == NO_TRACER .and. trim(co2_data_source) == 'predicted') then
        call error_mesg('radiation_driver_mod', &
        'co2 must be a tracer when predicted co2 desired for radiation.', FATAL)
      endif

!--------------------------------------------------------------------- 
!    call ozone_init to initialize the ozone field.
!---------------------------------------------------------------------
      call ozone_init (latb, lonb)

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine radiative_gases_init



!####################################################################


! <SUBROUTINE NAME="define_radiative_gases">
!  <OVERVIEW>
!   Subroutine that returns the current values of the radiative 
!    gas mixing ratios to radiation_driver in Rad_gases.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that returns the current values of the radiative 
!    gas mixing ratios to radiation_driver in Rad_gases.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_radiative_gases (is, ie, js, je, Rad_time, lat, &
!                                Atmos_input, Time_next, Rad_gases)
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which radiation is to be calculated
!                   [ time_type (days, seconds) ] 
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   latitude of model points [ radians ] 
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]
!  </IN>
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!  </INOUT>
! </SUBROUTINE>
!
subroutine define_radiative_gases (is, ie, js, je, Rad_time, lat, &
                                   Atmos_input, r, Time_next, Rad_gases)

!-------------------------------------------------------------------
!    define_radiative_gases returns the current values of the radiative 
!    gas mixing ratios to radiation_driver in Rad_gases.
!-------------------------------------------------------------------

integer,                    intent(in)    :: is, ie, js, je
type(time_type),            intent(in)    :: Rad_time, Time_next
real, dimension(:,:),       intent(in)    :: lat
type(atmos_input_type),     intent(in)    :: Atmos_input
real, dimension(:,:,:,:),   intent(in)    :: r
type(radiative_gases_type), intent(inout) :: Rad_gases

!---------------------------------------------------------------------
!
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Rad_Time     time at which radiation is to be calculated
!                   [ time_type (days, seconds) ] 
!      Time_next    time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]  
!      lat          latitude of model points  [ radians ]
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!      r            array of tracers, some of which may represent
!                   evolving, radiatively active gases (e.g. ozone) 
!
!
!  intent(inout) variables:
!
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
      character(len=8)   :: gas_name  ! name associated with the 
                                      ! radiative gas
      integer            :: yr, mo, dy, hr, mn, sc 
                                      ! components of Rad_time
      type(time_type)    :: Gas_time

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ( 'radiative_gases_mod', &
             'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    if time-varying co2 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of co2 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_co2) then 
        else
           if (trim(co2_data_source) == 'predicted') then
              if (.not. use_globally_uniform_co2) then    
! if predicted, not globally uniform co2 for radiation is desired, then
! rrvco2 will have become an array and will be filled here:
!               rrvco2(:,:,:) = Atmos_input%tracer_co2(is:ie,js:je,:)  
              endif
           endif
        endif

!--------------------------------------------------------------------
!    fill the contents of the radiative_gases_type variable which
!    will be returned to the calling routine.  values of the gas mixing
!    ratio at the current time and a flag indicating if the gas is time-
!    varying are returned for all gases, and for gases for which tfs are
!    calculated, a variable indicating how long they have been varying,
!    and the value of gas mixing ratio used for the last tf calculation
!    are returned.
!---------------------------------------------------------------------
!   these values must now be filled from the module variables:
!      Rad_gases%ch4_tf_offset = ch4_tf_offset
!      Rad_gases%n2o_tf_offset = n2o_tf_offset
!      Rad_gases%co2_tf_offset = co2_tf_offset
!      Rad_gases%ch4_for_next_tf_calc = ch4_for_next_tf_calc
!      Rad_gases%n2o_for_next_tf_calc = n2o_for_next_tf_calc
!      Rad_gases%co2_for_next_tf_calc = co2_for_next_tf_calc

!      Rad_gases%rrvch4  = rrvch4
!      Rad_gases%rrvn2o  = rrvn2o
!      Rad_gases%rrvf11  = rrvf11
!      Rad_gases%rrvf12  = rrvf12
!      Rad_gases%rrvf113 = rrvf113
!      Rad_gases%rrvf22  = rrvf22
!      Rad_gases%rrvco2  = rrvco2
!      Rad_gases%time_varying_co2  = time_varying_co2
!      Rad_gases%time_varying_ch4  = time_varying_ch4
!      Rad_gases%time_varying_n2o  = time_varying_n2o
!      Rad_gases%time_varying_f11  = time_varying_f11
!      Rad_gases%time_varying_f12  = time_varying_f12
!      Rad_gases%time_varying_f113 = time_varying_f113
!      Rad_gases%time_varying_f22  = time_varying_f22
!      if (time_varying_co2) then
!        Rad_gases%Co2_time = Co2_time_list(1)
!      endif
!      if (time_varying_ch4) then
!        Rad_gases%Ch4_time = Ch4_time_list(1)
!      endif
!      if (time_varying_n2o) then
!        Rad_gases%N2o_time = N2o_time_list(1)
!      endif
!RSH    define value for the new variable 
!RSH                     Rad_gases%use_model_supplied_co2, .true. for
!RSH   co2_data_source = 'predicted', .false. otherwise.
!      if (trim(co2_data_source) == 'predicted') then
!         Rad_gases%use_model_supplied_co2 = .true.
!      else
!         Rad_gases%use_model_supplied_co2 = .false.
!      endif

!      Rad_gases%co2_for_last_tf_calc = co2_for_last_tf_calc
!      Rad_gases%ch4_for_last_tf_calc = ch4_for_last_tf_calc
!      Rad_gases%n2o_for_last_tf_calc = n2o_for_last_tf_calc

!--------------------------------------------------------------------
!    allocate an array in a radiative_gases_type variable to hold the
!    model ozone field at the current time. call ozone_driver to define
!    this field for use in the radiation calculation.
!--------------------------------------------------------------------
      allocate (Rad_gases%qo3(ie-is+1, je-js+1,    &
                              size(Atmos_input%press,3) - 1))
      Rad_gases%qo3 = 0.
      call ozone_driver (is, ie, js, je, lat, Rad_time, Atmos_input, &
                         r, Rad_gases)


!---------------------------------------------------------------------


end subroutine define_radiative_gases


!#####################################################################

subroutine radiative_gases_time_vary (Rad_time, gavg_rrv, Rad_gases_tv)

!---------------------------------------------------------------------
!     subroutine radiative_gases_time_vary calculates time-dependent, space-
!     independent quantities needed by this module
!---------------------------------------------------------------------

type(time_type),    intent(in)   :: Rad_time
real, dimension(:), intent(in)   :: gavg_rrv
type(radiative_gases_type), intent(inout)  :: Rad_gases_tv

!---------------------------------------------------------------------
!  local variables:
!
      character(len=8)   :: gas_name  ! name associated with the 
                                      ! radiative gas
      integer            :: yr, mo, dy, hr, mn, sc 
                                      ! components of Rad_time
      type(time_type)    :: Gas_time

!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    if time-varying ch4 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of ch4 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_ch4) then 
            if (negative_offset_ch4) then
              Gas_time = Rad_time - Ch4_offset
            else
              Gas_time = Rad_time + Ch4_offset
            endif
          if (Gas_time >= Ch4_time_list(1)) then
!         if (Rad_time >= Ch4_time_list(1)) then
            gas_name = 'ch4 '
            call define_gas_amount      &
                (gas_name, Gas_time, ch4_specification_type,  &
!               (gas_name, Rad_time, ch4_specification_type,  &
!                negative_offset_ch4, Ch4_offset, &
                 ch4_variation_type, ch4_floor, ch4_ceiling, rch4,   &
                 ch4_uses_tfs, ch4_change_rate, rrvch4, Ch4_time_list, &
                 ch4_value, &
                 gas_tf_calc_intrvl =        &
                            Rad_control%ch4_tf_calc_intrvl,         &
                 gas_tf_time_displacement =  &
                            Rad_control%ch4_tf_time_displacement,   &
                 calc_gas_tfs_on_first_step =  &
                            Rad_control%calc_ch4_tfs_on_first_step, &
                 calc_gas_tfs_monthly       =  &
                            Rad_control%calc_ch4_tfs_monthly,       &
                 use_current_gas_for_tf = &
                            Rad_control%use_current_ch4_for_tf,  &
                 gas_tf_offset = &
                            ch4_tf_offset,  &
                 gas_for_last_tf_calc =   &
                            ch4_for_last_tf_calc,    &
                 gas_for_next_tf_calc = &
                            ch4_for_next_tf_calc, &
                 gas_tfs_needed = ch4_tfs_needed, &
                 define_gas_for_last_tf_calc = &
                            define_ch4_for_last_tf_calc)
                 if (Rad_control%calc_ch4_tfs_on_first_step) then
                   ch4_tfs_needed = .false.
                 endif

!---------------------------------------------------------------------
!    if time-variation is desired, but it is not yet time to begin
!    variation, define the ch4 mixing ratio that was used for the last
!    transmission function calculation, so that the tfs may be calcul-
!    ated as previously. (this is only done on initial timestep of a 
!    job).
!---------------------------------------------------------------------
          else  ! (Rad_time > Ch4_time)
            ch4_for_last_tf_calc = rrvch4
            ch4_tf_offset = 0.0
          endif   ! (Rad_time > Ch4_time)
        else
          ch4_for_last_tf_calc = rrvch4
        endif  ! (time_varying_ch4)

!--------------------------------------------------------------------
!    if time-varying n2o is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of n2o needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_n2o) then 
            if (negative_offset_n2o) then
              Gas_time = Rad_time - N2o_offset
            else
              Gas_time = Rad_time + N2o_offset
            endif
          if (Gas_time >= N2o_time_list(1)) then
            gas_name = 'n2o '
            call define_gas_amount      &
                (gas_name, Gas_time, n2o_specification_type,  &
!                negative_offset_n2o, N2o_offset, &
                 n2o_variation_type, n2o_floor, n2o_ceiling, rn2o,  &
                 n2o_uses_tfs, n2o_change_rate, rrvn2o, N2o_time_list, &
                 n2o_value, &
                 gas_tf_calc_intrvl =        &
                               Rad_control%n2o_tf_calc_intrvl,         &
                 gas_tf_time_displacement =  &
                               Rad_control%n2o_tf_time_displacement,   &
                 calc_gas_tfs_on_first_step =  &
                               Rad_control%calc_n2o_tfs_on_first_step, &
                 calc_gas_tfs_monthly       =  &
                            Rad_control%calc_n2o_tfs_monthly,       &
                 use_current_gas_for_tf = &
                               Rad_control%use_current_n2o_for_tf,  &
                 gas_tf_offset = &
                               n2o_tf_offset,  &
                 gas_for_last_tf_calc =   &
                               n2o_for_last_tf_calc,    &
                 gas_for_next_tf_calc = &
                               n2o_for_next_tf_calc, &
                 gas_tfs_needed = n2o_tfs_needed, &
                 define_gas_for_last_tf_calc = &
                               define_n2o_for_last_tf_calc)
                 if (Rad_control%calc_n2o_tfs_on_first_step) then
                   n2o_tfs_needed = .false.
                 endif

!---------------------------------------------------------------------
!    if time-variation is desired, but it is not yet time to begin
!    variation, define the n2o mixing ratio that was used for the last
!    transmission function calculation, so that the tfs may be calcul-
!    ated as previously. (this is only done on initial timestep of a 
!    job).
!---------------------------------------------------------------------
          else  ! (Rad_time > N2o_time)
            n2o_for_last_tf_calc = rrvn2o
            n2o_tf_offset = 0.0
          endif   ! (Rad_time > N2o_time)
        else
          n2o_for_last_tf_calc = rrvn2o
        endif  ! (time_varying_n2o)

!--------------------------------------------------------------------
!    if time-varying f11 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of f11 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_f11) then 
            if (negative_offset_f11) then
              Gas_time = Rad_time - F11_offset
            else
              Gas_time = Rad_time + F11_offset
            endif
          if (Gas_time >= F11_time_list(1)) then
            gas_name = 'f11 '
            call define_gas_amount      &
                (gas_name, Gas_time, f11_specification_type,  &
!                negative_offset_f11, F11_offset, &
                 f11_variation_type, f11_floor, f11_ceiling, rf11,   &
                 f11_uses_tfs, f11_change_rate, rrvf11, F11_time_list, &
                 f11_value)
          endif   ! (Rad_time > F11_time)
        endif  ! (time_varying_f11)

!--------------------------------------------------------------------
!    if time-varying f12 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of f12 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_f12) then 
            if (negative_offset_f12) then
              Gas_time = Rad_time - F12_offset
            else
              Gas_time = Rad_time + F12_offset
            endif
          if (Gas_time >= F12_time_list(1)) then
            gas_name = 'f12 '
            call define_gas_amount      &
                (gas_name, Gas_time, f12_specification_type,   &
!                negative_offset_f12, F12_offset, &
                 f12_variation_type, f12_floor, f12_ceiling, rf12,   &
                 f12_uses_tfs, f12_change_rate, rrvf12, F12_time_list, &
                 f12_value)
          endif   ! (Rad_time > F12_time)
        endif  ! (time_varying_f12)

!--------------------------------------------------------------------
!    if time-varying f113 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of f113 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_f113) then 
            if (negative_offset_f113) then
              Gas_time = Rad_time - F113_offset
            else
              Gas_time = Rad_time + F113_offset
            endif
          if (Gas_time >= F113_time_list(1)) then
            gas_name = 'f113'
            call define_gas_amount      &
                (gas_name, Gas_time, f113_specification_type,  &
!                negative_offset_f113, F113_offset, &
                 f113_variation_type, f113_floor, f113_ceiling, rf113, &
                 f113_uses_tfs, f113_change_rate, rrvf113,  &
                 F113_time_list, f113_value)
          endif   ! (Rad_time > F113_time)
        endif  ! (time_varying_f113)

!--------------------------------------------------------------------
!    if time-varying f22 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of f22 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_f22) then 
            if (negative_offset_f22) then
              Gas_time = Rad_time - F22_offset
            else
              Gas_time = Rad_time + F22_offset
            endif
          if (Gas_time >= F22_time_list(1)) then
            gas_name = 'f22 '
            call define_gas_amount      &
                (gas_name, Gas_time, f22_specification_type,  &
!                negative_offset_f22, F22_offset, &
                 f22_variation_type, f22_floor, f22_ceiling, rf22, &
                 f22_uses_tfs, f22_change_rate, rrvf22, F22_time_list, &
                 f22_value)
          endif   ! (Rad_time > F22_time)
        endif  ! (time_varying_f22)

!--------------------------------------------------------------------
!    if time-varying co2 is desired, and the time at which variation 
!    was to begin has been exceeded, define the gas_name variable and
!    call define_gas_amount to return the values of co2 needed on this
!    timestep.
!--------------------------------------------------------------------
        if (time_varying_co2) then 
            if (negative_offset_co2) then
              Gas_time = Rad_time - Co2_offset
            else
              Gas_time = Rad_time + Co2_offset
            endif
          if (Gas_time >= Co2_time_list(1)) then
            gas_name = 'co2 '
            call define_gas_amount      &
                (gas_name, Gas_time, co2_specification_type,  &
!                negative_offset_co2, Co2_offset, &
                 co2_variation_type, co2_floor, co2_ceiling, rco2, &
                 co2_uses_tfs, co2_change_rate, rrvco2, Co2_time_list, &
                 co2_value,  &
                 gas_tf_calc_intrvl =        &
                               Rad_control%co2_tf_calc_intrvl,         &
                 gas_tf_time_displacement =  &
                               Rad_control%co2_tf_time_displacement,   &
                 calc_gas_tfs_on_first_step =  &
                               Rad_control%calc_co2_tfs_on_first_step, &
                 calc_gas_tfs_monthly       =  &
                            Rad_control%calc_co2_tfs_monthly,       &
                 use_current_gas_for_tf = &
                               Rad_control%use_current_co2_for_tf,  &
                 gas_tf_offset = &
                               co2_tf_offset,  &
                 gas_for_last_tf_calc =   &
                               co2_for_last_tf_calc,    &
                 gas_for_next_tf_calc = &
                               co2_for_next_tf_calc, &
                 gas_tfs_needed = co2_tfs_needed, &
                 define_gas_for_last_tf_calc = &
                               define_co2_for_last_tf_calc)
                 if (Rad_control%calc_co2_tfs_on_first_step) then
                   co2_tfs_needed = .false.
                 endif

!---------------------------------------------------------------------
!    if time-variation is desired, but it is not yet time to begin
!    variation, define the co2 mixing ratio that was used for the last
!    transmission function calculation, so that the tfs may be calcul-
!    ated as previously. (this is only needed on initial timestep of a 
!    job).
!---------------------------------------------------------------------
          else  ! (Rad_time > Co2_time)
            co2_for_last_tf_calc = rrvco2
            co2_tf_offset = 0.0
          endif   ! (Rad_time > Co2_time)
        else
           if (trim(co2_data_source) == 'predicted') then
              if (use_globally_uniform_co2) then    
                 rrvco2 = gavg_rrv(ico2)
              else
!    if 3d co2 distribution desired for radiation, it will be defined in
!    define_radiatiove_gases.
              endif
           else !trim(co2_data_source) == 'predicted')
              co2_for_last_tf_calc = rrvco2
           endif  !(trim(co2_data_source) == 'predicted')
        endif  ! (time_varying_co2)


!---------------------------------------------------------------------
!    print out the current gas mixing ratios, if desired.
!---------------------------------------------------------------------
      if ((mpp_pe() == mpp_root_pe()) .and. verbose >= 3) then
        if (Rad_control%do_lw_rad) then
          print_alarm = print_alarm - Rad_control%lw_rad_time_step
        endif
        if (print_alarm <= 0) then
          call get_date (Rad_time, yr, mo, dy, hr, mn, sc)
          write (*, FMT ='(a, i5, i3, i3, i3, i3, i3,   &
                           & a, f11.6, a, f11.6, a, /,  &
                           & 7x, a, f11.6, a, f11.6,  &
                           & a, f11.6, a , /, 7x, a,  f11.6, &
                           & a, f11.6, a)')  &
                  'Time =',  yr, mo, dy, hr, mn, sc, &
                  '    ch4 = ', rrvch4*1.0e09, 'ppb  n2o = ' ,  &
                  rrvn2o*1.0e9, 'ppb', 'co2  = ', rrvco2*1.0e6,  &
                  'ppm  f11 = ' , rrvf11*1.0e12,  'ppt  f12 = ', &
                  rrvf12*1.0e12, 'ppt', 'f113 = ', rrvf113*1.0e12, &
                  'ppt  f22 = ', rrvf22*1.0e12, 'ppt'
          print_alarm = gas_printout_freq*3600
        endif
      endif

!--------------------------------------------------------------------
!    fill the contents of the radiative_gases_type variable which
!    will be returned to the calling routine.  values of the gas mixing
!    ratio at the current time and a flag indicating if the gas is time-
!    varying are returned for all gases, and for gases for which tfs are
!    calculated, a variable indicating how long they have been varying,
!    and the value of gas mixing ratio used for the last tf calculation
!    are returned.
!---------------------------------------------------------------------
!   these values must now be filled from the module variables:
      Rad_gases_tv%ch4_tf_offset = ch4_tf_offset
      Rad_gases_tv%n2o_tf_offset = n2o_tf_offset
      Rad_gases_tv%co2_tf_offset = co2_tf_offset
      Rad_gases_tv%ch4_for_next_tf_calc = ch4_for_next_tf_calc
      Rad_gases_tv%n2o_for_next_tf_calc = n2o_for_next_tf_calc
      Rad_gases_tv%co2_for_next_tf_calc = co2_for_next_tf_calc

      Rad_gases_tv%rrvch4  = rrvch4
      Rad_gases_tv%rrvn2o  = rrvn2o
      Rad_gases_tv%rrvf11  = rrvf11
      Rad_gases_tv%rrvf12  = rrvf12
      Rad_gases_tv%rrvf113 = rrvf113
      Rad_gases_tv%rrvf22  = rrvf22
      Rad_gases_tv%rrvco2  = rrvco2
      Rad_gases_tv%time_varying_co2  = time_varying_co2
      Rad_gases_tv%time_varying_ch4  = time_varying_ch4
      Rad_gases_tv%time_varying_n2o  = time_varying_n2o
      Rad_gases_tv%time_varying_f11  = time_varying_f11
      Rad_gases_tv%time_varying_f12  = time_varying_f12
      Rad_gases_tv%time_varying_f113 = time_varying_f113
      Rad_gases_tv%time_varying_f22  = time_varying_f22
      if (time_varying_co2) then
        Rad_gases_tv%Co2_time = Co2_time_list(1)
      endif
      if (time_varying_ch4) then
        Rad_gases_tv%Ch4_time = Ch4_time_list(1)
      endif
      if (time_varying_n2o) then
        Rad_gases_tv%N2o_time = N2o_time_list(1)
      endif
!    define value for the new variable 
!    Rad_gases%use_model_supplied_co2, .true. for
!                    co2_data_source = 'predicted', .false. otherwise.
      if (trim(co2_data_source) == 'predicted') then
         Rad_gases_tv%use_model_supplied_co2 = .true.
      else
         Rad_gases_tv%use_model_supplied_co2 = .false.
      endif

      Rad_gases_tv%co2_for_last_tf_calc = co2_for_last_tf_calc
      Rad_gases_tv%ch4_for_last_tf_calc = ch4_for_last_tf_calc
      Rad_gases_tv%n2o_for_last_tf_calc = n2o_for_last_tf_calc


!----------------------------------------------------------------------

            call ozone_time_vary (Rad_time)

!--------------------------------------------------------------------


end subroutine radiative_gases_time_vary



!##################################################################

subroutine radiative_gases_endts

     call ozone_endts


end subroutine radiative_gases_endts




!##################################################################

!####################################################################
! <SUBROUTINE NAME="radiative_gases_end">
!  <OVERVIEW>
!    radiative_gases_end is the destructor for radiative_gases_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiative_gases_end is the destructor for radiative_gases_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiative_gases_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine radiative_gases_end

!---------------------------------------------------------------------
!    radiative_gases_end is the destructor for radiative_gases_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ( 'radiative_gases_mod', &
               'module has not been initialized', FATAL )
      endif

      call radiative_gases_restart

!--------------------------------------------------------------------
!    deallocate the timeseries arrays.
!--------------------------------------------------------------------
      deallocate (co2_value, ch4_value, n2o_value, f11_value, &
                  f12_value, f113_value, f22_value)    
      deallocate (Co2_time_list, Ch4_time_list, N2o_time_list, &
                  F11_time_list, F12_time_List, &
                  F113_time_list, F22_time_list)

!---------------------------------------------------------------------
!    call the destructors for the component gas module(s).
!---------------------------------------------------------------------
      call ozone_end

!--------------------------------------------------------------------
      module_is_initialized = .false.


end subroutine radiative_gases_end

!#######################################################################
! <SUBROUTINE NAME="radiative_gases_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine radiative_gases_restart(timestamp)
   character(len=*), intent(in), optional :: timestamp

! Make sure that the restart_versions variable is up to date.
   vers = restart_versions(size(restart_versions(:)))     
   if( do_netcdf_restart ) then
      if(mpp_pe() == mpp_root_pe() ) then
         call error_mesg ('radiative_gases_mod', 'Writing NetCDF formatted restart file: RESTART/radiative_gases.res.nc', NOTE)
      endif
      call save_restart(Rad_restart, timestamp)
   else
      call error_mesg ('radiative_gases_mod', &
         'Native intermediate restart files are not supported.', FATAL)
   endif

end subroutine radiative_gases_restart
! </SUBROUTINE>

!####################################################################
! <SUBROUTINE NAME="radiative_gases_dealloc">
!  <OVERVIEW>
!    radiative_gases_end is the destructor for radiative_gases_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiative_gases_end is the destructor for radiative_gases_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiative_gases_dealloc (Rad_gases)
!  </TEMPLATE>
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!              ative gas input fields needed by the radiation package
!  </INOUT>
! </SUBROUTINE>
!
subroutine radiative_gases_dealloc (Rad_gases)

!---------------------------------------------------------------------
!
!--------------------------------------------------------------------

type(radiative_gases_type), intent(inout)  :: Rad_gases

!---------------------------------------------------------------------
!  intent(inout) variable:
!
!     Rad_gases
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ( 'radiative_gases_mod', &
             'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Rad_gases.
!--------------------------------------------------------------------
      deallocate (Rad_gases%qo3)


!--------------------------------------------------------------------
!    save the value of the gas mixing ratio used to calculate the gas
!    transmission functions (it may have been updated on this step)
!    so that it is available to write to the restart file.
!--------------------------------------------------------------------
      co2_for_last_tf_calc = Rad_gases%co2_for_last_tf_calc
      ch4_for_last_tf_calc = Rad_gases%ch4_for_last_tf_calc
      n2o_for_last_tf_calc = Rad_gases%n2o_for_last_tf_calc

!---------------------------------------------------------------------


end subroutine radiative_gases_dealloc 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="validate_time_varying_inputs">
!  <OVERVIEW>
!   validate_time_varying_inputs checks for consistency among the 
!   namelist parameters defining the time variation of the input gas.
!   NOTE: THIS IS A PRIVATE SUBROUTINE.
!
!  </OVERVIEW>
!
!  <DESCRIPTION>
!   Subroutine validate_time_varying_inputs performs the following 
!   checks of the namelist variables:
!     1) Verifies that base_time has a non-zero month and day number,
!        indicating a change from the default;
!     2) Verifies that the base value of gas mixing ratio is non-zero;
!     3) Verifies that specification_type is valid, either 
!        'base_and_trend' or 'time_series';
!     4) Verifies that variation_type is valid, either 'linear' or
!        'logarithmic';
!     5) When specification_type is 'base_and_trend', verifies that 
!        change_rate is non-zero;
!  </DESCRIPTION>
!  <TEMPLATE>
!   call validate_time_varying_inputs   &
!                     (gas, base_time, base_value, specification_type, &
!                      change_rate, variation_type, Gas_time_list,  &
!                      gas_value)
!  </TEMPLATE>
!  <IN NAME="gas">
!   name associated with the current gas
!  </IN>
!  <IN NAME="base_time">
!   time at which the base_value is applicable
!   [ year, month, day, hour, minute, second ]
!  </IN>
!  <IN NAME="base_value">
!   base value for vol. mixing ratio of gas 
!   [ number / number ]
!  </IN>
!  <IN NAME="specification_type">
!   specification of form of time variation of gas 
!  </IN>
!  <IN NAME="change_rate">
!   rate of change of gas; 1.00 + [ percent change / year ]
!  </IN>
!  <IN NAME="variation_type">
!   form of the temporal behavior of gas; either 'linear' or 
!   'logarithmic'
!  </IN>
!  <INOUT NAME="Gas_time_list">
!   array of time_type variables defining the data times for the gas;
!   if specification_type is timeseries, then it is the set of times
!   in the daa set, if specification type is base_and_trend, then it
!   is an array of dimension 1 containing the xxx_base_time. 
!   [ time_type ]
!  </INOUT>
!  <IN NAME="gas_value">
!   array of values of gas concentrations corresponding to the times
!   in Gas_time_list [ number / number ]
!  </IN>
!

!<PUBLICROUTINE>
!
!NOTE: THIS IS A PRIVATE SUBROUTINE.
!

subroutine validate_time_varying_inputs   &
                      (gas, base_time, base_value, specification_type, &
                       change_rate,  gas_dataset_entry,  &
                       negative_offset_gas, Gas_offset, Gas_entry, &
                       variation_type, Gas_time_list, gas_value)

!---------------------------------------------------------------------
!    validate_time_varying_inputs checks for consistency among the 
!    namelist parameters defining the time variation of the input gas.
!--------------------------------------------------------------------

character(len=*),      intent(in)  :: gas
integer, dimension(6), intent(in)  :: base_time, gas_dataset_entry
real,                  intent(in)  :: base_value,            &
                                      change_rate
logical,               intent(inout)  :: negative_offset_gas
character(len=*),     intent(in)  :: specification_type, variation_type
type(time_type),dimension(:), intent(inout)  :: Gas_time_list
type(time_type), intent(inout)  :: Gas_offset, Gas_entry
real, dimension(:), intent(in) :: gas_value

!--------------------------------------------------------------------
!   local variables:

      integer :: n

!</PUBLICROUTINE>

!------------------------------------------------------------------
!    perform checks for the base_and_trend specification_type.
!------------------------------------------------------------------
      if (trim(specification_type) == 'base_and_trend') then

!------------------------------------------------------------------
!    verify that base_time has a valid value (month and day /= 0 , 
!    hr, min and sec == 0). 
!------------------------------------------------------------------
        if (base_time(4) /= 0 .or. base_time(5) /= 0  .or. &
            base_time(6) /= 0)  then  
          call error_mesg ('radiative_gases_mod', &
           'base_time for gas' // trim(gas)//'must be specified as 00Z &
            &of desired day', FATAL)
        endif

!------------------------------------------------------------------
!    verify that base_value is non-zero. 
!------------------------------------------------------------------
        if (base_value < 0.0) then
          call error_mesg ('radiative_gases_mod', &
                        trim(gas)//'_base_value must be >= 0.0', FATAL)
        endif

!---------------------------------------------------------------------
!    convert the base_time to a time_type and store in Gas_time_list(1).
!---------------------------------------------------------------------
        if (base_time(2) /= 0 .and. base_time(3) /= 0 ) then  
          Gas_time_list(1) = set_date (base_time(1), base_time(2), &
                                       base_time(3), base_time(4), &
                                       base_time(5), base_time(6))  
        else
          call error_mesg ('radiative_gases_mod', &
           'must supply valid date for '//trim(gas)//'_base_time when&
                              & using time_varying '//trim(gas), FATAL)
        endif

!------------------------------------------------------------------
!    make sure that change_rate has been specified as non-zero; a value
!    of 0.0 corresponds to non-time-varying, and is not allowed.
!------------------------------------------------------------------
        if (change_rate == 0.0) then
          call error_mesg ('radiative_gases_mod', &
            ' have specified base_and_trend '//gas//' variation but ' &
                            //trim(gas)//'_change_rate is zero', FATAL)
        endif
        
!------------------------------------------------------------------
!    set Gas_offset to 0 when timeseries is not being used.
!------------------------------------------------------------------
        Gas_offset = set_time(0,0)

!------------------------------------------------------------------
!    perform checks for the time_series specification_type.
!------------------------------------------------------------------
      else if (trim(specification_type) ==  'time_series') then

!------------------------------------------------------------------
!    make sure that the entries are in chronological order.
!------------------------------------------------------------------
        do n = 2, size(Gas_time_list(:)) 
          if (Gas_time_list(n) < Gas_time_list(n-1)) then
            call error_mesg ('radiative_gases_mod', &
                 'times  for ' // trim(gas) //   &
                 ' in Gas_time_list are not in sequence', FATAL) 
          endif
        end do
      
!------------------------------------------------------------------
!    make sure that all gas concentrations are acceptable in magnitude.
!------------------------------------------------------------------
        do n = 1, size(Gas_time_list(:)) 
          if (gas_value(n) < 0) then
          call error_mesg ('radiative_gases_mod', &
             trim(gas)//'_value must be >= 0.0', FATAL)
          endif
        end do
        
!------------------------------------------------------------------
!    if 'time_series' is specified, variation_type
!    must be 'linear'.
!------------------------------------------------------------------
        if (trim(variation_type) == 'logarithmic') then
          call error_mesg ('radiative_gases_mod', &
            'logarithmic variation not allowed with time_series &
                                   &specification', FATAL)
        endif

!---------------------------------------------------------------------
!    define model initial time (from diag_table).
!---------------------------------------------------------------------
        Model_init_time = get_base_time()
 
!---------------------------------------------------------------------
!    if an entry into the gas timeseries has not been specified, use
!    the model base time as the entry point.
!---------------------------------------------------------------------
        if (gas_dataset_entry(1) == 1 .and. &
            gas_dataset_entry(2) == 1 .and. &
            gas_dataset_entry(3) == 1 .and. &
            gas_dataset_entry(4) == 0 .and. &
            gas_dataset_entry(5) == 0 .and. &
            gas_dataset_entry(6) == 0 ) then
          Gas_entry = Model_init_time

!---------------------------------------------------------------------
!    define time for which gas data is desired.
!---------------------------------------------------------------------
        else
          Gas_entry  = set_date (gas_dataset_entry(1), &
                                 gas_dataset_entry(2), &
                                 gas_dataset_entry(3), &
                                 gas_dataset_entry(4), &
                                 gas_dataset_entry(5), &
                                 gas_dataset_entry(6))
        endif

        call error_mesg ( 'radiative_gases_mod', &
              'PROCESSING TIMESERIES FOR ' // trim(gas), NOTE)
        call print_date (Gas_entry , str='Data from timeseries &
                                           &at time:')
        call print_date(Model_init_time , str='This data is mapped to&
                                                  & model time:')

!----------------------------------------------------------------------
!    define the offset from model base time to gas_dataset_entry
!    as a time_type variable.
!----------------------------------------------------------------------
        Gas_offset = Gas_entry - Model_init_time

        if (Model_init_time > Gas_entry) then
          negative_offset_gas = .true.
        else
          negative_offset_gas = .false.
        endif

!------------------------------------------------------------------
!    if specification_type is invalid, issue an error message.
!------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
               ' invalid '//trim(gas)//'_specification_type', FATAL)
      endif

!------------------------------------------------------------------
!    verify that variation_type is valid.
!------------------------------------------------------------------
      if (trim(variation_type) == 'linear' .or. &
          trim(variation_type) == 'logarithmic' ) then
      else
        call  error_mesg ('radiative_gases_mod', &
          trim(gas)//'_variation_type must be "linear" or &
                                             &"logarithmic" ', FATAL)
      endif

!---------------------------------------------------------------------



end subroutine validate_time_varying_inputs


! </SUBROUTINE>

!###################################################################
! <SUBROUTINE NAME="define_ch4">
!  <OVERVIEW>
!   Subroutine that provides initial values for ch4 mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that provides initial values for ch4 mixing ratio.if ch4
!    is fixed in time, the value is given by the namelist specification.
!    if ch4 is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_ch4 (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define ch4 initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_ch4 (data_source) 

!---------------------------------------------------------------------
!    define_ch4 provides initial values for ch4 mixing ratio. if ch4
!    is fixed in time, the value is given by the namelist specification.
!    if ch4 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define ch4 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rch4_ipcc_80  = 1.56900E-06
      real       ::  rch4_ipcc_92  = 1.71400E-06
      real       ::  rch4_icrccm   = 1.75000E-06
      real       ::  rch4_ipcc_98  = 1.82120E-06




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial ch4 mixing ratios to be used.
!    'icrccm'     --> rch4_icrccm
!    'ipcc80'     --> rch4_ipcc_80
!    'ipcc92'     --> rch4_ipcc_92
!    'ipcc98'     --> rch4_ipcc_98
!    'input'      --> file INPUT/id1ch4n2o, record 1
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(ch4_specification_type) /= 'time_series') then
        allocate (Ch4_time_list(1))
        allocate (ch4_value(1))
      endif

      if (trim(data_source)      == 'icrccm') then
        rch4   = rch4_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rch4   = rch4_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rch4   = rch4_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rch4   = rch4_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(ch4_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_namelist_file ('INPUT/id1ch4n2o')
            read (inrad, FMT = '(5e18.10)')  rch4
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired ch4_n2o input file is not present', FATAL)
          endif
        else
          gas_name = 'ch4 '
          call read_gas_timeseries (gas_name, ch4_value, &
                                    Ch4_time_list, time_varying_ch4, &
                                    ch4_dataset_entry, rch4)   
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart ch4 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_namelist_file ('INPUT/id1ch4n2o')
            read (inrad, FMT = '(5e18.10)') rch4
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor ch4 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_namelist_file ('INPUT/id1ch4n2o')
            read (inrad, FMT = '(5e18.10)') rch4
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor ch4 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of ch4 is obtained
!    from the namelist variable ch4_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rch4 = ch4_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for ch4' , FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rch4 < ch4_floor) then
        call error_mesg ('radiative_gases_mod', &
              'base ch4 mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rch4 > ch4_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base ch4 mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------



end subroutine define_ch4



!#####################################################################
! <SUBROUTINE NAME="define_n2o">
!  <OVERVIEW>
!   Subroutine that provides initial values for n2o mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that provides initial values for n2o mixing ratio.if n2o
!    is fixed in time, the value is given by the namelist specification.
!    if n2o is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_n2o (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define n2o initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_n2o (data_source) 

!---------------------------------------------------------------------
!    define_n2o provides initial values for n2o mixing ratio. if n2o
!    is fixed in time, the value is given by the namelist specification.
!    if n2o is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define n2o initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rn2o_ipcc_80  = 3.02620E-07
      real       ::  rn2o_ipcc_92  = 3.11000E-07
      real       ::  rn2o_icrccm   = 2.80000E-07
      real       ::  rn2o_ipcc_98  = 3.16000E-07




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial n2o mixing ratios to be used.
!    'icrccm'     --> rn2o_icrccm
!    'ipcc80'     --> rn2o_ipcc_80
!    'ipcc92'     --> rn2o_ipcc_92
!    'ipcc98'     --> rn2o_ipcc_98
!    'input'      --> file INPUT/id1ch4n2o, record 2
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(n2o_specification_type) /= 'time_series') then
        allocate (N2o_time_list(1))
        allocate (n2o_value(1))
      endif

      if (trim(data_source)      == 'icrccm') then
        rn2o   = rn2o_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rn2o   = rn2o_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rn2o   = rn2o_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rn2o   = rn2o_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(n2o_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_namelist_file ('INPUT/id1ch4n2o')
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  rn2o
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired ch4_n2o input file is not present', FATAL)
          endif
        else
          gas_name = 'n2o '
          call read_gas_timeseries (gas_name, n2o_value,   &
                                    N2o_time_list, time_varying_n2o, &
                                    n2o_dataset_entry, rn2o)   
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart n2o values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_namelist_file ('INPUT/id1ch4n2o')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rn2o
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor n2o input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_namelist_file ('INPUT/id1ch4n2o')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rn2o
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor n2o input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of n2o is obtained
!    from the namelist variable n2o_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rn2o = n2o_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for n2o ', FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rn2o < n2o_floor) then
        print *, 'rn2o, n2o_floor', rn2o, n2o_floor, mpp_pe()
        call error_mesg ('radiative_gases_mod', &
              'base n2o mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rn2o > n2o_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base n2o mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------




end subroutine define_n2o



!#####################################################################
! <SUBROUTINE NAME="define_f11">
!  <OVERVIEW>
!   Subroutine that provides initial values for f11 mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_f11 provides initial values for f11 mixing ratio. if f11
!    is fixed in time, the value is given by the namelist specification.
!    if f11 is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_f11 (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define f11 initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_f11 (data_source) 

!---------------------------------------------------------------------
!    define_f11 provides initial values for f11 mixing ratio. if f11
!    is fixed in time, the value is given by the namelist specification.
!    if f11 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f11 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf11_icrccm   = 1.00000E-09
      real       ::  rf11_ipcc_80  = 1.57500E-10
      real       ::  rf11_ipcc_92  = 2.68000E-10
      real       ::  rf11_ipcc_98  = 2.68960E-10




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial f11 mixing ratios to be used.
!    'icrccm'     --> rf11_icrccm
!    'ipcc80'     --> rf11_ipcc_80
!    'ipcc92'     --> rf11_ipcc_92
!    'ipcc98'     --> rf11_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 1
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(f11_specification_type) /= 'time_series') then
        allocate (F11_time_list(1))
        allocate (f11_value(1))
     endif

      if (trim(data_source)      == 'icrccm') then
        rf11   = rf11_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf11   = rf11_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf11   = rf11_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf11   = rf11_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(f11_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)')  rf11
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired cfc input file is not present', FATAL)
          endif
        else
          gas_name = 'f11 '
          call read_gas_timeseries (gas_name, f11_value,   &
                                    F11_time_list, time_varying_f11, &
                                    f11_dataset_entry, rf11)   
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart f11 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') rf11
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f11 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') rf11
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f11 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of f11 is obtained
!    from the namelist variable f11_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rf11 = f11_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f11 ', FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rf11 < f11_floor) then
        call error_mesg ('radiative_gases_mod', &
              'base f11 mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rf11 > f11_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base f11 mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------




end subroutine define_f11




!#####################################################################
! <SUBROUTINE NAME="define_f12">
!  <OVERVIEW>
!   Subroutine that provides initial values for f12 mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_f12 provides initial values for f12 mixing ratio. if f12
!    is fixed in time, the value is given by the namelist specification.
!    if f12 is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_f12 (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define f12 initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_f12 (data_source) 

!---------------------------------------------------------------------
!    define_f12 provides initial values for f12 mixing ratio. if f12
!    is fixed in time, the value is given by the namelist specification.
!    if f12 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f12 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf12_icrccm   = 1.00000E-09
      real       ::  rf12_ipcc_80  = 2.72500E-10
      real       ::  rf12_ipcc_92  = 5.03000E-10
      real       ::  rf12_ipcc_98  = 5.31510E-10




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial f12 mixing ratios to be used.
!    'icrccm'     --> rf12_icrccm
!    'ipcc80'     --> rf12_ipcc_80
!    'ipcc92'     --> rf12_ipcc_92
!    'ipcc98'     --> rf12_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 2
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(f12_specification_type) /= 'time_series') then
        allocate (F12_time_list(1))
        allocate (f12_value(1))
      endif

      if (trim(data_source)      == 'icrccm') then
        rf12   = rf12_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf12   = rf12_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf12   = rf12_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf12   = rf12_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(f12_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  rf12
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired cfc input file is not present', FATAL)
          endif
        else
          gas_name = 'f12 '
          call read_gas_timeseries (gas_name, f12_value,   &
                                    F12_time_list, time_varying_f12, &
                                    f12_dataset_entry, rf12)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart f12 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf12
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f12 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf12
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f12 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of f12 is obtained
!    from the namelist variable f12_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rf12 = f12_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f12', FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rf12 < f12_floor) then
        call error_mesg ('radiative_gases_mod', &
              'base f12 mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rf12 > f12_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base f12 mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------




end subroutine define_f12




!#####################################################################
! <SUBROUTINE NAME="define_f113">
!  <OVERVIEW>
!   Subroutine that provides initial values for f113 mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_f113 provides initial values for f113 mixing ratio. if f113
!    is fixed in time, the value is given by the namelist specification.
!    if f113 is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_f113 (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define f113 initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_f113 (data_source) 

!---------------------------------------------------------------------
!    define_f113 provides initial values for f113 mixing ratio. if f113
!    is fixed in time, the value is given by the namelist specification.
!    if f113 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f113 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf113_icrccm  = 1.00000E-09
      real       ::  rf113_ipcc_80 = 2.31400E-11
      real       ::  rf113_ipcc_92 = 8.20000E-11
      real       ::  rf113_ipcc_98 = 8.58100E-11




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial f113 mixing ratios to be used.
!    'icrccm'     --> rf113_icrccm
!    'ipcc80'     --> rf113_ipcc_80
!    'ipcc92'     --> rf113_ipcc_92
!    'ipcc98'     --> rf113_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 3
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(f113_specification_type) /= 'time_series') then
        allocate (F113_time_list(1))
        allocate (f113_value(1))
      endif

      if (trim(data_source)      == 'icrccm') then
        rf113   = rf113_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf113   = rf113_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf113   = rf113_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf113   = rf113_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(f113_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  rf113
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired f113 input file is not present', FATAL)
          endif
        else
          gas_name = 'f113'
          call read_gas_timeseries (gas_name, f113_value,   &
                                  F113_time_list, time_varying_f113, &
                                  f113_dataset_entry, rf113 )
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart f113 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf113
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f113 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf113
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f113 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of f113 is obtained
!    from the namelist variable f113_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rf113 = f113_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f113', FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rf113 < f113_floor) then
        call error_mesg ('radiative_gases_mod', &
              'base f113 mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rf113 > f113_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base f113 mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------




end subroutine define_f113



!#####################################################################
! <SUBROUTINE NAME="define_f22">
!  <OVERVIEW>
!   Subroutine that provides initial values for f22 mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_f22 provides initial values for f22 mixing ratio. if f22
!    is fixed in time, the value is given by the namelist specification.
!    if f22 is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_f22 (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define f22 initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_f22 (data_source) 

!---------------------------------------------------------------------
!    define_f22 provides initial values for f22 mixing ratio. if f22
!    is fixed in time, the value is given by the namelist specification.
!    if f22 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f22 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf22_icrccm   = 1.00000E-09
      real       ::  rf22_ipcc_80  = 6.20200E-11
      real       ::  rf22_ipcc_92  = 1.05000E-10
      real       ::  rf22_ipcc_98  = 1.26520E-10




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial f22 mixing ratios to be used.
!    'icrccm'     --> rf22_icrccm
!    'ipcc80'     --> rf22_ipcc_80
!    'ipcc92'     --> rf22_ipcc_92
!    'ipcc98'     --> rf22_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 4
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(f22_specification_type) /= 'time_series') then
        allocate (F22_time_list(1))
        allocate (f22_value(1))
      endif

      if (trim(data_source)      == 'icrccm') then
        rf22   = rf22_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf22   = rf22_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf22   = rf22_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf22   = rf22_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(f22_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  
            read (inrad, FMT = '(5e18.10)')  rf22
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired cfc input file is not present', FATAL)
          endif
        else
          gas_name = 'f22 '
          call read_gas_timeseries (gas_name, f22_value,  &
                                    F22_time_list, time_varying_f22, &
                                    f22_dataset_entry, rf22)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart f22 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf22
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f22 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_namelist_file ('INPUT/id1cfc')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf22
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor f22 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of f22 is obtained
!    from the namelist variable f22_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rf22 = f22_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f22', FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rf22 < f22_floor) then
        call error_mesg ('radiative_gases_mod', &
              'base f22 mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rf22 > f22_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base f22 mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------




end subroutine define_f22



!#####################################################################
! <SUBROUTINE NAME="define_co2">
!  <OVERVIEW>
!   Subroutine that provides initial values for co2 mixing ratio.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_co2 provides initial values for co2 mixing ratio. if co2
!    is fixed in time, the value is given by the namelist specification.
!    if co2 is time-varying, the values are obtained from either a
!    restart or input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_co2 (data_source) 
!  </TEMPLATE>
!  <IN NAME="data_source" TYPE="character">
!   character string defining source to use to define co2 initial values
!  </IN>
! </SUBROUTINE>
!
subroutine define_co2 (data_source) 

!---------------------------------------------------------------------
!    define_co2 provides initial values for co2 mixing ratio. if co2
!    is fixed in time, the value is given by the namelist specification.
!    if co2 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define co2 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rco2_icrccm   = 3.00000E-04
      real       ::  rco2_ipcc_92  = 3.56000E-04
      real       ::  rco2_ipcc_80  = 3.37320E-04
      real       ::  rco2_ipcc_98  = 3.69400E-04
      real       ::  rco2_330ppm   = 3.30000E-04
      real       ::  rco2_660ppm   = 6.60000E-04
      real       ::  rco2_720ppm   = 7.20000E-04




      integer    :: inrad   ! unit number for i/o

character(len=8)     :: gas_name ! name associated with current
                                 ! gas being processed

!---------------------------------------------------------------------
!    define initial co2 mixing ratios to be used.
!    'icrccm'     --> rco2_icrccm
!    'ipcc80'     --> rco2_ipcc_80
!    'ipcc92'     --> rco2_ipcc_92
!    'ipcc98'     --> rco2_ipcc_98
!    '330ppm'     --> rco2_330ppm 
!    '660ppm'     --> rco2_660ppm 
!    'input'      --> file INPUT/id1co2, record 2
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if specification_type is 'base_and_trend', allocate length (1)
!    arrays to hold the base_time and base_value.
!---------------------------------------------------------------------
      if (trim(co2_specification_type) /= 'time_series') then
        allocate (Co2_time_list(1))
        allocate (co2_value(1))
      endif

      if (trim(data_source)      == 'icrccm') then
        rco2   = rco2_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rco2   = rco2_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rco2   = rco2_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rco2   = rco2_ipcc_98

      else if (trim(data_source) == '330ppm') then
        rco2   = rco2_330ppm

      else if (trim(data_source) == '660ppm') then
        rco2   = rco2_660ppm

      else if (trim(data_source) == '720ppm') then
        rco2   = rco2_720ppm

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (trim(co2_specification_type) /= 'time_series') then
          if (file_exist ('INPUT/id1co2') ) then
            inrad = open_namelist_file ('INPUT/id1co2')
            read (inrad, FMT = '(5e18.10)')  rco2
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
                   'desired co2 input file is not present', FATAL)
          endif
        else
          gas_name = 'co2 '
          call read_gas_timeseries (gas_name, co2_value, &
                                    Co2_time_list, time_varying_co2, &
                                    co2_dataset_entry, rco2)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ( 'radiative_gases_mod', &
           'cannot use restart co2 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1co2') ) then
            inrad = open_namelist_file ('INPUT/id1co2')
            read (inrad, FMT = '(5e18.10)') rco2
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor co2 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"prescribed" ', FATAL)
          endif
        endif

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1co2') ) then
            inrad = open_namelist_file ('INPUT/id1co2')
            read (inrad, FMT = '(5e18.10)') rco2
            call close_file (inrad)
            co2_for_last_tf_calc = rco2
          else
            call error_mesg ( 'radiative_gases_mod', &
              'neither restart nor co2 input file is present. one '//&
          'of these is required when the data_source is  '//&
          '"predicted" ', FATAL)
          endif
        endif

!-------------------------------------------------------------------
!    when the data_source is 'namelist', the value of co2 is obtained
!    from the namelist variable co2_base_value.
!-------------------------------------------------------------------
      else if (trim(data_source) == 'namelist') then
        rco2 = co2_base_value

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for co2 input', FATAL)
      endif

!---------------------------------------------------------------------
!    be sure value is within range of acceptable values. a later check
!    in lw_gases_stdtf_mod will further limit values to be those for
!    which tfs may be determined.
!---------------------------------------------------------------------
      if (rco2 < co2_floor) then
        call error_mesg ('radiative_gases_mod', &
              'base co2 mixing ratio LOWER THAN FLOOR value', FATAL)
      endif
      if (rco2 > co2_ceiling) then
        call error_mesg ('radiative_gases_mod', &
              ' base co2 mixing ratio HIGHER THAN CEILING value', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_co2


!#####################################################################
! <SUBROUTINE NAME="read_gas_timeseries">
!  <OVERVIEW>
!   Subroutine that reads in data values for well-mixed greenhouse
!   gases at specified times.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_gas_timeseries obtains global values for well-mixed
!    greenhouse gases from observed data sources for gas xxx.
!     the data are obtained from an input file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_gas_timeseries (gas, gas_value, Gas_time_list, rgas)
!  </TEMPLATE>
!  <IN NAME="gas">
!   name associated with the current gas
!  </IN>
!  <INOUT NAME="gas_value">
!   array of volume mixing ratio of gas 'gas' for each time in 
!   Gas_time_list gas_year. [no. /no. ]
!  </INOUT>
!  <INOUT NAME="Gas_time_list">
!   list of times (time_type) associated with the gas_value data 
!  </INOUT>
!  <OUT NAME="rgas">
!   gas volume mixing ratio at the start of the timeseries
!  </OUT>
! </SUBROUTINE>
!
subroutine read_gas_timeseries (gas, gas_value, Gas_time_list,   &
                                time_varying_gas, gas_dataset_entry, &
                                rgas)
 
!--------------------------------------------------------------------
!    read_gas_timeseries obtains global values for well-mixed
!    greenhouse gases from observed data sources for gas xxx.
!     the data are obtained from an input file.
!--------------------------------------------------------------------

character(len=*),                   intent(in)    :: gas
real, dimension(:),                 pointer       :: gas_value
type(time_type), dimension(:),      pointer       :: Gas_time_list
logical,                            intent(in)    :: time_varying_gas
integer,dimension(:),               intent(in)    :: gas_dataset_entry
real,                               intent(out)   :: rgas

!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!      gas     name associated with the current gas
!
!   pointer variables:
!
!      gas_value
!      Gas_time_list
!
!   intent(out) variables:
!
!      rgas
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variables :

      real,dimension(:), allocatable :: input_time
      type(time_type)                :: Year_t, Extra_time
      character(len=64)              :: file_name

      real       :: extra_seconds
      integer    :: inrad   ! unit number for i/o
      integer    :: i       ! do loop index
      integer    :: year, diy, yr, mo, dy, hr, mn, sc
      integer    :: series_length
      integer    :: index1, index2
      real       :: percent_of_period
      type(time_type) :: Gas_entry
      character(len=256) :: err_msg

!-------------------------------------------------------------------
!    define the gas_name which is currently being processed.
!-------------------------------------------------------------------
      file_name = 'INPUT/' // trim(gas) //'_gblannualdata'

!--------------------------------------------------------------------
!    process the gas timeseries file.
!--------------------------------------------------------------------
      if (file_exist (file_name) ) then
        inrad = open_namelist_file (file_name)

!--------------------------------------------------------------------
!    read the number of data points in the timeseries.
!--------------------------------------------------------------------
        read (inrad, FMT = '(i12)')  series_length

!----------------------------------------------------------------------
!    allocate gas_value, input_time and Gas_time_list arrays.
!----------------------------------------------------------------------
        allocate (gas_value(series_length))
        allocate (input_time(series_length))
        allocate (Gas_time_list(series_length))

!---------------------------------------------------------------------
!    read the timeseries data, time and then gas value.
!---------------------------------------------------------------------
        do i=1,series_length
          read (inrad, FMT = '(2f12.4)') input_time(i), &
                                         gas_value(i)
        end do

!---------------------------------------------------------------------
!    convert from input units  to [ no. / no. ].
!---------------------------------------------------------------------
        if ( trim(gas) == 'co2') then
          gas_value = 1.0e-6*gas_value
        else if (trim(gas) == 'ch4' .or. &
                 trim(gas) == 'n2o' ) then   
          gas_value = 1.0e-9*gas_value
        else if (trim(gas) == 'f11' .or. &
                 trim(gas) == 'f12' .or. &   
                 trim(gas) == 'f113' .or. &   
                 trim(gas) == 'f22' ) then   
          gas_value = 1.0e-12*gas_value
        endif

!---------------------------------------------------------------------
!    close the input file.
!---------------------------------------------------------------------
        call close_file (inrad)

!---------------------------------------------------------------------
!    convert the time stamps of the series to time_type variables.     
!---------------------------------------------------------------------
        if (verbose > 3) then
          if ( mpp_pe() == mpp_root_pe() ) then
            print *, 'time series entries for ' // trim (gas)
          endif
        endif
        do i=1,series_length
          year = INT(input_time(i))
          Year_t = set_date (year,1,1,0,0,0)
          diy = days_in_year (Year_t)
          extra_seconds = (input_time(i) - year)*diy*86400. 
          Extra_time=    set_time(NINT(extra_seconds), 0)
          Gas_time_list(i)    = Year_t + Extra_time
          call get_date (Gas_time_list(i), yr, mo, dy, hr, mn, sc)
          if (verbose > 3) then
            if ( mpp_pe() == mpp_root_pe() ) then
              print *, i, yr, mo, dy, hr, mn, sc, gas_value(i)
            endif
          endif
        end do

!--------------------------------------------------------------------
!    if the gas is not time-varying, its value must be defined here.
!--------------------------------------------------------------------
        if (.not. time_varying_gas) then

!--------------------------------------------------------------------
!    if a dataset entry time has not been specified, send an error
!    message.
!--------------------------------------------------------------------
          if (gas_dataset_entry(1) == 1 .and. &
              gas_dataset_entry(2) == 1 .and. &
              gas_dataset_entry(3) == 1 .and. &
              gas_dataset_entry(4) == 0 .and. &
              gas_dataset_entry(5) == 0 .and. &
              gas_dataset_entry(6) == 0 ) then   
            call error_mesg ('radiative_gases_mod', &
             'timeseries selected but no valid datset_entry time &
                             &provided for ' // trim(gas) , FATAL)

!---------------------------------------------------------------------
!    convert the dataset entry that is provided to a time_type variable,
!    and determine the timeseries value that corresponds.
!---------------------------------------------------------------------
          else
            Gas_entry = set_date (gas_dataset_entry(1), &
                                  gas_dataset_entry(2), &
                                  gas_dataset_entry(3), &
                                  gas_dataset_entry(4), &
                                  gas_dataset_entry(5), &
                                  gas_dataset_entry(6))     
            call time_interp (Gas_entry, Gas_time_list,  &
                              percent_of_period, index1, index2, err_msg=err_msg)
            if(err_msg /= '') then
               call error_mesg('radiative_gases_mod ',trim(err_msg)//' file='//trim(file_name), FATAL)
            endif
            rgas = gas_value(index1) + percent_of_period*  &
                   (gas_value(index2) - gas_value(index1))
            call error_mesg ( 'radiative_gases_mod', &
                   'PROCESSING TIMESERIES FOR ' // trim(gas), NOTE)
            call print_date (Gas_entry , str='Gas value is taken from &
                                     &timeseries at time:')
            if (mpp_pe() == mpp_root_pe() ) then
              print *, trim(gas) // ' value is ', rgas
            endif
          endif

!---------------------------------------------------------------------
!    if gas is time_varying, define the initial gas mixing ratio to be 
!    the first value in the timeseries. this value will be replaced 
!    before it is used.
!---------------------------------------------------------------------
        else
          rgas = gas_value(1)
        endif

!---------------------------------------------------------------------
!    if the requested input file is not present, write an error message.
!---------------------------------------------------------------------
      else
        call error_mesg ( 'radiative_gases_mod', &
             'desired ' // file_name // ' input file is not present', &
                                                                FATAL)
      endif

!--------------------------------------------------------------------
 


end subroutine read_gas_timeseries
 

!#####################################################################
! <SUBROUTINE NAME="define_gas_amount">
!  <OVERVIEW>
!   define_gas_amount defines the values of the gas mixing ratio needed !   at the current time, when the gas is varying with time.
!  </OVERVIEW>
!  <DESCRIPTION>
!   define_gas_amount performs the following actions:
!     1) checks for the presence of needed optional arguments;
!     2) determines how long the gas has been varying with time;
!     3) calculates values for the gas mixing ratio at the current time;
!     4) constrains calculated gas values to lie between the
!        specified floor and ceiling;
!
!   if transmission functions are calculated for the gas, then:
!     1) it is determined if the gas value used when the tfs were last
!        calculated is needed;
!     2) if the gas does not use its current value at the time when tfs
!        are calculated, the offset from the current time to the time
!        used for tfs is obtained;
!     3) if the gas value used when the tfs were last calculated is 
!        needed, it is calculated along with the time offset of that
!        time from the present time;
!     4) if the gas value at the time when the tfs are next to be cal-
!        culated is needed, it is calculated;
!     5) gas values relevant at the time when tfs are next to be 
!        calculated are constrained to lie between the specified floor 
!        and ceiling;
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_gas_amount      &
!        (gas, Rad_time, gas_specification_type, gas_variation_type, &
!         gas_floor, gas_ceiling, rgas, gas_uses_tfs, gas_change_rate, &
!         rrvgas, Gas_time_list, gas_value, gas_tf_calc_intrvl,  &
!         gas_tf_time_displacement, calc_gas_tfs_on_first_step,  &
!         use_current_gas_for_tf, gas_tf_offset, gas_for_last_tf_calc, &
!         gas_for_next_tf_calc, define_gas_for_last_tf_calc)
!  </TEMPLATE>
!  <IN NAME="gas">
!   character string associated with the gas being processed
!  </IN>
!  <IN NAME="Rad_time">
!   time at which radiation calculation is to apply
!  </IN>
!  <IN NAME="gas_specification_type">
!   indicator as to the form of time variation of vol. mixing ratio;
!   either 'base_and_trend' or 'time_series'.
!  </IN>
!  <IN NAME="gas_variation_type">
!   indicator as to the form of time variation of the vol. mixing ratio
!   of gas; either 'linear' or 'logarithmic'
!  </IN>
!  <IN NAME="gas_floor">
!   smallest value allowed for gas xxx vol. mixing ratio [ no. / no. ]
!  </IN>
!  <IN NAME="gas_ceiling">
!   largest value allowed for gas xxx vol. mixing ratio [ no. / no. ]
!  </IN>
!  <IN NAME="rgas">
!   initially specified gas mixing ratio [ no. / no. ]
!  </IN>
!  <IN NAME="gas_uses_tfs">
!   this gas has transmission functions associated with it ?
!  </IN>
!  <INOUT NAME="gas_change_rate">
!   time rate of change of gas xxx vol. mixing ratio
!   [  1 +/- % per year ]
!  </INOUT>
!  <INOUT NAME="rrvgas">
!   gas mixing ratio at current time [ no. / no. ]
!  </INOUT>
!  <IN NAME="Gas_time_list">
!   list of times in gas timeseries [ time_type ]
!  </IN>
!  <IN NAME="gas_value">
!   gas concentrations [ no. / no. ] associated with the times 
!   in Gas_time_list
!  </IN>
!  <IN NAME="gas_tf_calc_intrvl">
!   time interval between calculating gas tfs  [ hours ]
!   OPTIONAL: present only when the gas has tfs associated with it
!  </IN>
!  <IN NAME="gas_tf_time_displacement">
!   time displacement from present to the time at which gas values are
!   to be used in the calculation of tfs. may be <0, ==0, or > 0.
!   [ hours ]
!   OPTIONAL: present only when the gas has tfs associated with it, only
!   used when calc_gas_tfs_on_first_step is .true.
!  </IN>
!  <IN NAME="calc_gas_tfs_on_first_step">
!   if true, tfs are calculated ONLY on the first time step of a run,
!   using gas mixing ratios valid gas_tf_time_displacement hours from 
!   the start time
!   OPTIONAL: present only when the gas has tfs associated with it
!  </IN>
!  <IN NAME="use_current_gas_for_tf">
!   if true, the gas  mixing ratio at the current time is used to cal-
!   culate the gas tfs
!   OPTIONAL: present only when the gas has tfs associated with it
!  </IN>
!  <OUT NAME = "gas_tf_offset">
!   time between last tf calculation and present [ hours ]
!   OPTIONAL: present only when the gas has tfs associated with it
!  </OUT>
!  <OUT NAME = "gas_for_last_tf_calc">
!   value of gas mixing ratio used in last tf calculation [ no. / no. ]
!   OPTIONAL: present only when the gas has tfs associated with it
!  </OUT>
!  <OUT NAME = "gas_for_next_tf_calc">
!   value of gas mixing ratio to be used in next tf calculation 
!   OPTIONAL: present only when the gas has tfs associated with it
!   [ no. / no. ]
!  </OUT>
!  <INOUT NAME = "define_gas_for_last_tf_calc">
!   logical indicating if the gas value used for the last tf calculation
!   must be obtained
!  </INOUT>
!

!<PUBLICROUTINE>
!
!NOTE: THIS IS A PRIVATE SUBROUTINE>
!
subroutine define_gas_amount      &
         (gas, Rad_time, gas_specification_type,  &
!         negative_offset_gas, Gas_offset, gas_variation_type, &
                                           gas_variation_type, &
          gas_floor, gas_ceiling, rgas, gas_uses_tfs, gas_change_rate, &
          rrvgas, Gas_time_list, gas_value, gas_tf_calc_intrvl,  &
          gas_tf_time_displacement, calc_gas_tfs_on_first_step,  &
          calc_gas_tfs_monthly, &
          use_current_gas_for_tf, gas_tf_offset, gas_for_last_tf_calc, &
          gas_for_next_tf_calc,  gas_tfs_needed,  &
          define_gas_for_last_tf_calc)

!--------------------------------------------------------------------
character(len=*),              intent(in)    :: gas
type(time_type),               intent(in)    :: Rad_time
character(len=*),              intent(in)    :: gas_specification_type,&
                                                gas_variation_type
real,                          intent(in)    :: gas_floor, gas_ceiling,&
                                                rgas
logical,                       intent(in)    :: gas_uses_tfs  
!logical,                       intent(in)    :: negative_offset_gas  
real,                          intent(inout) :: gas_change_rate, rrvgas

type(time_type), dimension(:), intent(in)    :: Gas_time_list
!type(time_type),               intent(inout) :: Gas_offset      
real, dimension(:),            intent(in)    :: gas_value
real,               intent(in),    optional  :: gas_tf_calc_intrvl,    &
                                                gas_tf_time_displacement
logical,            intent(in),    optional  ::   &
                                           gas_tfs_needed, &
                                          calc_gas_tfs_on_first_step, &
                                          calc_gas_tfs_monthly,  &
                                          use_current_gas_for_tf 
real,               intent(out),   optional  :: gas_tf_offset, &
                                                gas_for_last_tf_calc, &
                                                gas_for_next_tf_calc
logical,            intent(inout), optional  :: &  
                                          define_gas_for_last_tf_calc

!---------------------------------------------------------------------
!  local variables:

     type(time_type)    :: Gas_yrs   
     integer            :: days, seconds
     real               :: years_of_gas, years_of_gas_till_next
     integer            :: days2, seconds2
     integer            :: days3, seconds3
     real               :: mean_days, calc_time
     character(len=16)  :: chvers7, chvers8, chvers9
     integer            :: alarm, minutes_from_start

     real               :: percent_of_period
     type(time_type)    :: Tf_offset, Tf_calc_intrvl 
     real               :: rseconds3
     integer            :: index1, index2
     integer            :: yr, mo, dy, hr, mn, sc
     integer            :: days7, seconds7
     type(time_type)    :: Tf_displ, First_of_month, Gas_tf_next, &
                           Time_left
     character(len=256) :: err_msg
!---------------------------------------------------------------------
!  local variables:
!    
!     Gas_yrs                 time interval from start of time variation
!                             until current time [ time_type ]
!     days                    days component of Gas_yrs  [ days ]   
!     seconds                 seconds component of Gas_yrs  [ seconds ]
!     minutes_from_start      time interval from start of time variation
!                             until current time [ minutes ]
!     years_of_gas            time interval from start of time variation
!                             until current time [ years ]
!     years_of_gas_till_next  time interval from start of time variation
!                             until next tf calculation [ years ]
!     days2                   days component of the mean length of year
!                             time_type variable [ days ]
!     seconds2                seconds component of the mean length of 
!                             year time_type variable [ seconds ]
!     mean_days               average number of days in a year [ days ]
!     calc_time               time at which tfs were last calculated
!                             [ years from start of gas time variation ]
!     chvers7                 character variable used to output model
!                             variables through error_mesg interface
!     chvers8                 character variable used to output model
!                             variables through error_mesg interface
!     chvers9                 character variable used to output model
!                             variables through error_mesg interface
!     chvers11                character variable used to output model
!                             variables through error_mesg interface
!     alarm                   time since last tf calculation until
!                             current time [ minutes ]              
!
!--------------------------------------------------------------------

!     type(time_type)  :: Gas_time  ! time for which gas data is desired

!</PUBLICROUTINE>

!--------------------------------------------------------------------
!    define the time for which gas data is desired.
!--------------------------------------------------------------------
!     if (negative_offset_gas) then
!       Gas_time = Rad_time - Gas_offset
!     else
!       Gas_time = Rad_time + Gas_offset
!     endif

!---------------------------------------------------------------------
!    if this gas calculates transmission functions, make sure all 
!    optional arguments needed when tfs are in use are present.
!---------------------------------------------------------------------
      if (gas_uses_tfs) then
        if (present( gas_tf_calc_intrvl) .and. &
            present( gas_tf_time_displacement) .and. &
            present( gas_tf_offset) .and. &
            present( gas_for_last_tf_calc) .and. &
            present( gas_for_next_tf_calc) .and. &
            present( calc_gas_tfs_on_first_step) .and. &
            present( calc_gas_tfs_monthly) .and. &
            present( define_gas_for_last_tf_calc) .and. &
            present( use_current_gas_for_tf) ) then 
        else
          call error_mesg ('radiative_gases_mod', &
          'necessary optional arguments for '//trim(gas)//' call to&
           & subroutine define_gas_amount are not present', FATAL)
        endif
      endif

!--------------------------------------------------------------------
!    define the mean length of the year in units of days. this will
!    be a function of the calendar being used. this will be the unit
!    used to define how long gas variation has been occurring 
!    (years_of_gas), not exactly equivalent to calendar years.
!--------------------------------------------------------------------
      call get_time (length_of_year(), seconds2, days2)
      mean_days = days2 + seconds2/86400.

!---------------------------------------------------------------------
!    define how long the gas variation has been occurring, expressed
!    as a time_type, the components of the time_type and then in units
!    of gas-years.
!---------------------------------------------------------------------
      Gas_yrs = Rad_time - Gas_time_list(1)
!     Gas_yrs = Gas_time - Gas_time_list(1)
      call get_time (Gas_yrs, seconds, days)
      years_of_gas = (days + seconds/86400.)/mean_days

!---------------------------------------------------------------------
!    define the current value of gas. this value will be used in sw
!    calculations. the following expressions are available:
!      base_and_trend, logarithmic: 
!                    g(t) = g(t0)*(gas_change_rate)**(t-t0)
!      base_and_trend, linear:   
!                    g(t) = g(t0)*(1.0 + (t-t0)*(gas_change_rate - 1.0)
!    where t0 is the base time.
!---------------------------------------------------------------------
      if (trim(gas_specification_type) =='base_and_trend') then
        if (trim(gas_variation_type) == 'logarithmic') then
          if (gas_change_rate /=0.0) then
            rrvgas = rgas* exp(alog(gas_change_rate)*years_of_gas)
          else
            rrvgas = rgas
          endif
        else 
          rrvgas = rgas*(1.0 + years_of_gas*(gas_change_rate-1.0))
        endif

      else if (trim(gas_specification_type) == 'time_series') then
        call time_interp (Rad_time, Gas_time_list,   &
!       call time_interp (Gas_time, Gas_time_list,   &
                          percent_of_period, index1, index2, err_msg=err_msg)
        if(err_msg /= '') then
           call error_mesg('radiative_gases_mod 1',trim(err_msg), FATAL)
        endif
        rrvgas   = gas_value(index1) + percent_of_period*  &
                   (gas_value(index2) - gas_value(index1))
      endif

!---------------------------------------------------------------------
!    be sure that newly calculated current gas mixing ratio remains 
!    within the floor / ceiling range. if either is exceeded, reset the
!    gas amount to the floor/ ceiling value.
!---------------------------------------------------------------------
      if (rrvgas < gas_floor) then
        if (verbose >= 1) then
          if (.not. printed_current_floor_msg) then
            write (chvers7, '(3pe15.7)') rrvgas
            write (chvers8, '(3pe15.7)') gas_floor
            write (chvers9, '(f9.5)') years_of_gas
            call error_mesg ('radiative_gases_mod', &
           'calculated '//trim(gas)//' mixing ratio ('//chvers7//  &
           ') LOWER THAN FLOOR ('//chvers8//') after'//chvers9//'years&
           & of '//trim(gas)// ' variation; reset to floor value ',  &
                                                                  NOTE)
            printed_current_floor_msg = .true.
          endif
        endif ! (verbose)
        rrvgas = gas_floor
      endif
      if (rrvgas > gas_ceiling) then
        if (verbose >= 1) then
          if (.not. printed_current_ceiling_msg) then
            write (chvers7, '(3pe15.7)') rrvgas
            write (chvers8, '(3pe15.7)') gas_ceiling
            write (chvers9, '(f 9.5)') years_of_gas
            call error_mesg ('radiative_gases_mod', &
           'calculated '//trim(gas)// ' mixing ratio ('//chvers7// &
           ') HIGHER THAN CEILING ('//chvers8//') after'//chvers9// &
           'years of '//trim(gas)// ' variation; reset to ceiling &
                                                        &value ', NOTE)
            printed_current_ceiling_msg = .true.
          endif
        endif ! (verbose)
        rrvgas = gas_ceiling
      endif

!---------------------------------------------------------------------
!    execute the following code if tfs are to be calculated for this 
!    gas.
!---------------------------------------------------------------------
      if (gas_uses_tfs) then
!     if (gas_uses_tfs .and. gas_tfs_needed) then
        if (gas_tfs_needed) then

!---------------------------------------------------------------------
!    if this gas uses tfs, determine if the gas mixing ratio used the
!    last time tfs were calculated is needed. in general it is  
!    needed if it was not read in from the restart file (a non-
!    current version), and gas time variation has already begun. 
!    however if time variation has not yet begun (the usual case if 
!    using an old restart) or if the tfs are to be calculated on the 
!    first step of the job, it will not be needed, and the flag is 
!    reset to so indicate. 
!---------------------------------------------------------------------
        if (years_of_gas == 0.0 .or. calc_gas_tfs_on_first_step .or. &
             calc_gas_tfs_monthly) then
          define_gas_for_last_tf_calc =.false.
        endif

!--------------------------------------------------------------------
!    define the time to which the gas mixing ratio used in the last
!    calculation corresponded. 
!--------------------------------------------------------------------
        if (define_gas_for_last_tf_calc) then
          calc_time = years_of_gas - MOD (years_of_gas, &
                      gas_tf_calc_intrvl/(24.0*mean_days))
          if (verbose >= 4 .and. mpp_pe() == mpp_root_pe() ) then
            print *, 'last tf:calc_time, years_of_gas', calc_time, &
                     years_of_gas
          endif

!--------------------------------------------------------------------
!    if the value of the gas mixing ratio used for the last tf calcul-
!    ation is needed, define it and the time to which it applies. after
!    calculation, set the flag indicating its need to .false.
!--------------------------------------------------------------------
          if (trim(gas_specification_type) =='base_and_trend') then
            if (trim(gas_variation_type) == 'logarithmic') then
              gas_for_last_tf_calc = rgas*exp(alog(gas_change_rate)* &
                                              calc_time)
              gas_tf_offset = (calc_time - years_of_gas)*24.0*mean_days
            else 
              gas_for_last_tf_calc  = rgas*(1.0 + calc_time*  &
                                      (gas_change_rate-1.0))
              gas_tf_offset = (calc_time - years_of_gas)*24.0*mean_days
            endif
          else if (trim(gas_specification_type) == 'time_series') then
            Tf_calc_intrvl = set_time (NINT(calc_time*24.0*  &
                                            mean_days*3600.), 0)
            call time_interp (Rad_time - Tf_calc_intrvl, Gas_time_list, &
!           call time_interp (Gas_time - Tf_calc_intrvl, Gas_time_list, &
                              percent_of_period, index1, index2, err_msg=err_msg)
            if(err_msg /= '') then
               call error_mesg('radiative_gases_mod 2',trim(err_msg), FATAL)
            endif
            gas_for_last_tf_calc   = gas_value(index1) +    &
                                     percent_of_period*  &
                                     (gas_value(index2) -   &
                                      gas_value(index1))
            if (mpp_pe() == mpp_root_pe() )then
              print *, 'gas_for_last_tf_calc, , *3,  &
                       &index1, index2, days3, seconds3, % ',   &
                         gas_for_last_tf_calc,   &
                         index1, index2, percent_of_period
            endif
            gas_tf_offset = (calc_time - years_of_gas)*24.0* mean_days
          endif
          define_gas_for_last_tf_calc = .false.
        endif !(define_gas_for_last_tf_calc) 

!---------------------------------------------------------------------
!    if tfs are calculated using other than the gas values at the 
!    current time, define the time that has elapsed since the beginning
!    of variation, and determine if this is a time at which tfs are
!    due to be calculated (alarm = 0).
!---------------------------------------------------------------------
        if (.not. use_current_gas_for_tf) then
          minutes_from_start = INT(days*1440.0 + real(seconds)/60.)
          if (gas_tf_calc_intrvl /= 0.0) then
            alarm = MOD(minutes_from_start,   &
                        INT(gas_tf_calc_intrvl*60.0))
          else
            alarm = 0
          endif

!---------------------------------------------------------------------
!    if alarm is 0 (indicating this is a step on which to calculate the
!    tfs), or if the option has been chosen to always and only calculate
!    tfs on the first step of a job, define the time to use to obtain 
!    the gas mixing ratio used to calculate those tfs.
!---------------------------------------------------------------------
          if (alarm == 0 .or. calc_gas_tfs_on_first_step .or. &
             calc_gas_tfs_monthly) then

!--------------------------------------------------------------------
!    if calc_gas_tfs_on_first_step is true, the gas mixing ratio at a
!    time gas_tf_time_displacement hours from now (plus, zero or minus 
!    allowed) will be used.
!--------------------------------------------------------------------
            if (calc_gas_tfs_on_first_step  ) then
              years_of_gas_till_next = years_of_gas +   &
                              gas_tf_time_displacement/(24.0*mean_days)
           else if (calc_gas_tfs_monthly) then
              call get_date (Rad_time, yr, mo, dy, hr, mn, sc)
              Tf_displ =   &
                      set_time(NINT(gas_tf_time_displacement*60*60), 0)
              First_of_month = set_date (yr, mo, 1, 0, 0, 0)
              Gas_tf_next =  First_of_month + Tf_displ
              if (Gas_tf_next > Rad_time) then
                Time_left = Gas_tf_next - Rad_time
                call get_time (Time_left, seconds7, days7)
                years_of_gas_till_next = years_of_gas + (days7 +   &
                                         seconds7/86400.)/mean_days
              else
                Time_left = Rad_time - Gas_tf_next
                call get_time (Time_left, seconds7, days7)
                years_of_gas_till_next = years_of_gas - (days7 +   &
                                         seconds7/86400.)/mean_days
!               years_of_gas_till_next = years_of_gas            
              endif
!--------------------------------------------------------------------
!    if alarm is 0, the gas mixing ratio at the mid point of the time
!    interval between this calculation time and the next is used.
!--------------------------------------------------------------------
            else if (alarm == 0 ) then   
              years_of_gas_till_next = years_of_gas +   &
                              0.5*(gas_tf_calc_intrvl)/(24.0*mean_days)
            endif

!--------------------------------------------------------------------
!    calculate the difference in time (hours) between current time and 
!    time used to define the gas mixing ratio used for the next tf
!    calculation.
!--------------------------------------------------------------------
            gas_tf_offset = (years_of_gas_till_next - years_of_gas)* &
                             24.0*mean_days

!--------------------------------------------------------------------
!    if the value of the gas mixing ratio to be used for the next tf 
!    calculation is needed, define it here. it will be needed if the
!    gas mixing ratio to be used in defining the tfs is not the
!    current value and either it is time to do the calculation or the
!    calculation is desired on the first step of the job.
!--------------------------------------------------------------------
            if (trim(gas_specification_type) =='base_and_trend') then
              if (trim(gas_variation_type) == 'logarithmic') then
                gas_for_next_tf_calc = rgas*exp(alog(gas_change_rate)*&
                                       years_of_gas_till_next)
              else 
                gas_for_next_tf_calc = rgas*(1.0 +   &
                                       years_of_gas_till_next*  &
                                       (gas_change_rate-1.0))
              endif
            else if (trim(gas_specification_type) == 'time_series') then
             if (calc_gas_tfs_monthly) then
               Tf_offset = Time_left
               if (gas_tf_offset > 0) then
                 call time_interp (Rad_time + Tf_offset, Gas_time_list,&
                                percent_of_period, index1, index2, err_msg=err_msg)
               else
                 call time_interp (Rad_time - Tf_offset, Gas_time_list,&
                                percent_of_period, index1, index2, err_msg=err_msg)
               endif
               if(err_msg /= '') then
                  call error_mesg('radiative_gases_mod 3',trim(err_msg), FATAL)
               endif
             else
               days3 = NINT(gas_tf_offset/24.0)
               rseconds3 = (gas_tf_offset - days3*24)*3600.0
               seconds3 = NINT(rseconds3)
               Tf_offset = set_time (seconds3, days3)
               call time_interp (Rad_time + Tf_offset, Gas_time_list,  &
                                percent_of_period, index1, index2,err_msg=err_msg)
               if(err_msg /= '') then
                  call error_mesg('radiative_gases_mod 4',trim(err_msg), FATAL)
               endif
              endif
              gas_for_next_tf_calc   = gas_value(index1) +    &
                                       percent_of_period*  &
                                       (gas_value(index2) -   &
                                        gas_value(index1))
            endif

!---------------------------------------------------------------------
!    if the current value is not being used, be sure that the gas mixing
!    ratio calculated for use when tfs are next calculated is within 
!    the floor / ceiling range. if either is exceeded, reset the gas 
!    amount to the floor/ ceiling value.
!---------------------------------------------------------------------
            if (gas_for_next_tf_calc < gas_floor) then
              if (verbose >= 1) then
                if (.not. printed_next_floor_msg) then
                  write (chvers7, '(3pe15.7)')  gas_for_next_tf_calc
                  write (chvers8, '(3pe15.7)') gas_floor
                  write (chvers9, '(f9.5)') years_of_gas_till_next
                  call error_mesg ('radiative_gases_mod', &    
                  'calculated '//trim(gas)// ' mixing ratio to be used&
                  & for tf calcs ('//chvers7//') LOWER  &
                  &THAN FLOOR ('//chvers8//') after'//chvers9//'years&
                  & of '//trim(gas)// ' variation; reset to floor  &
                                                         &value ', NOTE)
                  printed_next_floor_msg = .true.
                endif
              endif ! (verbose)
              gas_for_next_tf_calc = gas_floor
            endif
            if (gas_for_next_tf_calc > gas_ceiling) then
              if (verbose >= 1) then
                if (.not. printed_next_ceiling_msg) then
                  write (chvers7, '(3pe15.7)')   gas_for_next_tf_calc
                  write (chvers8, '(3pe15.7)') gas_ceiling
                  write (chvers9, '(f 9.5)') years_of_gas_till_next
                  call error_mesg ('radiative_gases_mod', &
                 'calculated '//trim(gas)// ' mixing ratio to be used&
                  & for tf calcs (' //chvers7//') HIGHER  &
                  &THAN CEILING ('//chvers8//') after'//chvers9//'years&
                  & of '//trim(gas)// ' variation; reset to ceiling &
                                                        &value ', NOTE)
                  printed_next_ceiling_msg = .true.
                endif
              endif ! (verbose)
              gas_for_next_tf_calc        = gas_ceiling
            endif
          else

!--------------------------------------------------------------------
!    set the value to be the current value. in this case it will not 
!    be used.
!--------------------------------------------------------------------
            gas_for_next_tf_calc = rrvgas
          endif

!---------------------------------------------------------------------
!    if the current value of the gas mixing ratio is to be used for the
!    next tf calculation, reset the values to the just-adjusted current
!    value rrvgas and set the time offset to 0.0.
!---------------------------------------------------------------------
        else ! (.not. use_current_gas_for_tf)
          gas_for_next_tf_calc = rrvgas
          gas_tf_offset = 0.0
        endif
 else ! (gas_tfs_needed)
          define_gas_for_last_tf_calc =.false.
          gas_for_next_tf_calc = rrvgas
          gas_tf_offset = 0.0
 endif ! (gas_tfs_needed)
      endif ! (gas_uses_tfs) 

!---------------------------------------------------------------------



end subroutine define_gas_amount 


! </SUBROUTINE>

!####################################################################
! <SUBROUTINE NAME="write_restart_radiative_gases">
!  <OVERVIEW>
!   Subroutine to write the radiative restart files
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to write the radiative restart files
!  </DESCRIPTION>
!  <TEMPLATE>
!   call write_restart_radiative_gases
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine write_restart_radiative_gases

!---------------------------------------------------------------------
!    write_restart_radiative_gases writes the radiative_gases.res file.
!---------------------------------------------------------------------

      integer    :: unit    ! unit number for i/o



!---------------------------------------------------------------------
!    open unit and write radiative gas restart file.
!---------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
         call error_mesg ('radiative_gases_mod', 'Writing native formatted restart file: RESTART/radiative_gases.res', NOTE)
        unit = open_restart_file ('RESTART/radiative_gases.res',   &
                                   action= 'write')
        write (unit) restart_versions(size(restart_versions(:)))
        write (unit) rrvco2
        write (unit) rrvf11, rrvf12, rrvf113, rrvf22
        write (unit) rrvch4, rrvn2o
        write (unit) co2_for_last_tf_calc
        write (unit) ch4_for_last_tf_calc
        write (unit) n2o_for_last_tf_calc
        call close_file (unit)
      endif

!----------------------------------------------------------------------
      


end subroutine write_restart_radiative_gases

!####################################################################


                  end module radiative_gases_mod

