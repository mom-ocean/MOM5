!
!  coupler_main couples component models and controls the time integration
!
program coupler_main
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> V. Balaji </CONTACT>


! <OVERVIEW>
!  A main program that couples component models for atmosphere, ocean, land, 
!  and sea ice on independent grids. 
! </OVERVIEW>

! <DESCRIPTION>
!  This version couples model components representing atmosphere, ocean, land 
!  and sea ice on independent grids. Each model component is represented by a 
!  data type giving the instantaneous model state.
!
!  The component models are coupled to allow implicit vertical diffusion of 
!  heat and moisture at the interfaces of the atmosphere, land, and ice models. 
!  As a result, the atmosphere, land, and ice models all use the same time step. 
!  The atmospheric model has been separated into down and up calls that 
!  correspond to the down and up sweeps of the standard tridiagonal elimination.
!
!  The ocean interface uses explicit mixing. Fluxes to and from the ocean must
!  be passed through the ice model. This includes atmospheric fluxes as well as 
!  fluxes from the land to the ocean (runoff).
!
!  This program contains the model's main time loop. Each iteration of the 
!  main time loop is one coupled (slow) time step. Within this slow time step 
!  loop is a fast time step loop, using the atmospheric time step, where the
!  tridiagonal vertical diffusion equations are solved. Exchange between sea 
!  ice and ocean occurs once every slow timestep.
!
! <PRE>
!      MAIN PROGRAM EXAMPLE
!      --------------------
!
!         DO slow time steps (ocean)
!
!              call flux_ocean_to_ice
!
!              call ICE_SLOW_UP
!
!              DO fast time steps (atmos)
!
!                   call flux_calculation
!
!                   call ATMOS_DOWN
!
!                   call flux_down_from_atmos
!
!                   call LAND_FAST
!
!                   call ICE_FAST
!
!                   call flux_up_to_atmos
!
!                   call ATMOS_UP
!
!              END DO
!
!              call ICE_SLOW_DN
!
!              call flux_ice_to_ocean
!
!              call OCEAN
!
!         END DO

!  </PRE>

! </DESCRIPTION>
! <INFO>
!   <NOTE>
!     <PRE>
!   1.If no value is set for current_date, start_date, or calendar (or default value 
!     specified) then the value from restart file "INPUT/coupler.res" will be used. 
!     If neither a namelist value or restart file value exist the program will fail. 
!   2.The actual run length will be the sum of months, days, hours, minutes, and 
!     seconds. A run length of zero is not a valid option. 
!   3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>

!   <ERROR MSG="no namelist value for current_date " STATUS="FATAL">
!     A namelist value for current_date must be given if no restart file for
!     coupler_main (INPUT/coupler.res) is found. 
!   </ERROR>
!   <ERROR MSG="invalid namelist value for calendar" STATUS="FATAL">
!     The value of calendar must be 'julian', 'noleap', or 'thirty_day'. 
!     See the namelist documentation. 
!   </ERROR>
!   <ERROR MSG="no namelist value for calendar" STATUS="FATAL">
!     If no restart file is present, then a namelist value for calendar 
!     must be specified. 
!   </ERROR>
!   <ERROR MSG="initial time is greater than current time" STATUS="FATAL">
!     If a restart file is present, then the namelist value for either 
!     current_date or start_date was incorrectly set. 
!   </ERROR>
!   <ERROR MSG="run length must be multiple of ocean time step " STATUS="FATAL">
!     There must be an even number of ocean time steps for the requested run length. 
!   </ERROR>
!   <ERROR MSG="final time does not match expected ending time " STATUS="WARNING">
!     This error should probably not occur because of checks done at initialization time. 
!   </ERROR>

! </INFO>

  use constants_mod,           only: constants_init

  use time_manager_mod,        only: time_type, set_calendar_type, set_time
  use time_manager_mod,        only: set_date, get_date, days_in_month, month_name
  use time_manager_mod,        only: operator(+), operator(-), operator (<)
  use time_manager_mod,        only: operator (>), operator ( /= ), operator ( / )
  use time_manager_mod,        only: operator (*), THIRTY_DAY_MONTHS, JULIAN
  use time_manager_mod,        only: NOLEAP, NO_CALENDAR, INVALID_CALENDAR
  use time_manager_mod,        only: date_to_string, increment_date
  use time_manager_mod,        only: operator(>=), operator(<=), operator(==)

  use fms_mod,                 only: open_namelist_file, field_exist, file_exist, check_nml_error
  use fms_mod,                 only: uppercase, error_mesg, write_version_number
  use fms_mod,                 only: fms_init, fms_end, stdout
  use fms_mod,                 only: read_data, write_data

  use fms_io_mod,              only: fms_io_exit
  use fms_io_mod,              only: restart_file_type, register_restart_field, save_restart

  use diag_manager_mod,        only: diag_manager_init, diag_manager_end, diag_grid_end
  use diag_manager_mod,        only: DIAG_OCEAN, DIAG_OTHER, DIAG_ALL, get_base_date
  use diag_manager_mod,        only: diag_manager_set_time_end

  use field_manager_mod,       only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE

  use tracer_manager_mod,      only: tracer_manager_init, get_tracer_index
  use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, NO_TRACER

  use coupler_types_mod,       only: coupler_types_init

  use data_override_mod,       only: data_override_init

!
! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!

  use atmos_model_mod,         only: atmos_model_init, atmos_model_end
  use atmos_model_mod,         only: update_atmos_model_down
  use atmos_model_mod,         only: update_atmos_model_up
  use atmos_model_mod,         only: atmos_data_type
  use atmos_model_mod,         only: land_ice_atmos_boundary_type
  use atmos_model_mod,         only: atmos_data_type_chksum
  use atmos_model_mod,         only: lnd_ice_atm_bnd_type_chksum
  use atmos_model_mod,         only: lnd_atm_bnd_type_chksum
  use atmos_model_mod,         only: ice_atm_bnd_type_chksum
  use atmos_model_mod,         only: atmos_model_restart

  use land_model_mod,          only: land_model_init, land_model_end
  use land_model_mod,          only: land_data_type, atmos_land_boundary_type
  use land_model_mod,          only: update_land_model_fast, update_land_model_slow
  use land_model_mod,          only: atm_lnd_bnd_type_chksum
  use land_model_mod,          only: land_data_type_chksum
  use land_model_mod,          only: land_model_restart

  use ice_model_mod,           only: ice_model_init, ice_model_end
  use ice_model_mod,           only: update_ice_model_slow_up
  use ice_model_mod,           only: update_ice_model_fast
  use ice_model_mod,           only: update_ice_model_slow_dn
  use ice_model_mod,           only: ice_data_type, land_ice_boundary_type
  use ice_model_mod,           only: ocean_ice_boundary_type, atmos_ice_boundary_type
  use ice_model_mod,           only: ice_model_restart
  use ice_model_mod,           only: ice_data_type_chksum, ocn_ice_bnd_type_chksum
  use ice_model_mod,           only: atm_ice_bnd_type_chksum, lnd_ice_bnd_type_chksum

  use ocean_model_mod,         only: update_ocean_model, ocean_model_init
  use ocean_model_mod,         only: ocean_model_end, ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ocean_model_mod,         only: ocean_model_restart
  use ocean_model_mod,         only: ocean_public_type_chksum, ice_ocn_bnd_type_chksum
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!

  use flux_exchange_mod,       only: flux_exchange_init
  use flux_exchange_mod,       only: sfc_boundary_layer
  use flux_exchange_mod,       only: generate_sfc_xgrid
  use flux_exchange_mod,       only: flux_down_from_atmos
  use flux_exchange_mod,       only: flux_up_to_atmos
  use flux_exchange_mod,       only: flux_land_to_ice
  use flux_exchange_mod,       only: flux_ice_to_ocean
  use flux_exchange_mod,       only: flux_ocean_to_ice
  use flux_exchange_mod,       only: flux_check_stocks, flux_init_stocks, flux_ice_to_ocean_stocks, flux_ocean_from_ice_stocks

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data

  use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_chksum
  use mpp_mod,                 only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod,                 only: stderr, stdlog, mpp_error, NOTE, FATAL, WARNING
  use mpp_mod,                 only: mpp_set_current_pelist, mpp_declare_pelist
  use mpp_mod,                 only: input_nml_file

  use mpp_io_mod,              only: mpp_open, mpp_close, mpp_io_clock_on
  use mpp_io_mod,              only: MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use mpp_domains_mod,         only: mpp_broadcast_domain

  use memutils_mod,            only: print_memuse_stats

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.F90,v 20.0 2013/12/13 23:27:07 fms Exp $'
  character(len=128) :: tag = '$Name: tikal $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land
  type   (ice_data_type) :: Ice
  ! allow members of ocean type to be aliased (ap)
  type (ocean_public_type), target :: Ocean
  type (ocean_state_type),  pointer :: Ocean_state => NULL()

  type(atmos_land_boundary_type)     :: Atmos_land_boundary
  type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
  type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
  type(land_ice_boundary_type)       :: Land_ice_boundary
  type(ice_ocean_boundary_type)      :: Ice_ocean_boundary
  type(ocean_ice_boundary_type)      :: Ocean_ice_boundary

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type (time_type) :: Time, Time_init, Time_end, &
                      Time_step_atmos, Time_step_cpld
  type(time_type) :: Time_atmos, Time_ocean
  integer :: num_atmos_calls, na
  integer :: num_cpld_calls, nc

!------ for intermediate restart
  type(restart_file_type), allocatable :: Ice_bc_restart(:), Ocn_bc_restart(:)
  character(len=64),       allocatable :: ice_bc_restart_file(:), ocn_bc_restart_file(:) 
  integer                              :: num_ice_bc_restart=0, num_ocn_bc_restart=0
  type(time_type)                      :: Time_restart, Time_restart_current, Time_start
  character(len=32)                    :: timestamp

! ----- coupled model initial date -----

  integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)
  integer :: calendar_type = INVALID_CALENDAR

!-----------------------------------------------------------------------
!------ namelist interface -------

! <NAMELIST NAME="coupler_nml">
!   <DATA NAME="current_date"  TYPE="integer, dimension(6)"  DEFAULT="0">
!     The date that the current integration starts with. 
!   </DATA>
!   <DATA NAME="force_date_from_namelist"  TYPE="logical"  DEFAULT=".false.">
!     Flag that determines whether the namelist variable current_date should 
!     override the date in the restart file INPUT/coupler.res. If the restart 
!     file does not exist then force_date_from_namelist has not effect, the value of current_date 
!     will be used.
!   </DATA>
!   <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
!     The calendar type used by the current integration. Valid values are consistent 
!     with the time_manager module: 'julian', 'noleap', or 'thirty_day'. The value 
!     'no_calendar' can not be used because the time_manager's date  function are used. 
!     All values must be lowercase.
!   </DATA>
!   <DATA NAME="months "  TYPE="integer"  DEFAULT="0">
!     The number of months that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="days "  TYPE="integer"  DEFAULT="0">
!     The number of days that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="hours"  TYPE="integer"  DEFAULT="0">
!     The number of hours that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="minutes "  TYPE="integer"  DEFAULT="0">
!     The number of minutes that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="seconds"  TYPE="integer"  DEFAULT="0">
!     The number of seconds that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="dt_atmos"  TYPE="integer"  DEFAULT="0">
!     Atmospheric model time step in seconds, including the fast coupling with 
!     land and sea ice. 
!   </DATA>
!   <DATA NAME="dt_cpld"  TYPE="integer"  DEFAULT="0">
!     Time step in seconds for coupling between ocean and atmospheric models: 
!     must be an integral multiple of dt_atmos and dt_ocean. This is the "slow" timestep.
!   </DATA>
!  <DATA NAME="do_atmos, do_ocean, do_ice, do_land, do_flux" TYPE="logical">
!  If true (default), that particular model component (atmos, etc.) is run.
!  If false, the execution of that component is skipped. This is used when
!  ALL the output fields sent by that component to the coupler have been
!  overridden using the data_override feature. For advanced users only:
!  if you're not sure, you should leave these values at TRUE.
!  </DATA> 
!  <DATA NAME="concurrent" TYPE="logical">
!  If true, the ocean executes concurrently with the atmosphere-land-ocean
!   on a separate set of PEs.
!  If false (default), the execution is serial: call atmos... followed by
!  call ocean...
!  If using concurrent execution, you must set one of
!   atmos_npes and ocean_npes, see below.
!  </DATA> 
!  <DATA NAME="atmos_npes, ocean_npes" TYPE="integer">
!  If concurrent is set to true, we use these to set the list of PEs on which
!   each component runs.
!  At least one of them must be set to a number between 0 and NPES.
!  If exactly one of these two is set non-zero, the other is set to the
!   remainder from NPES.
!  If both are set non-zero they must add up to NPES.
!  </DATA> 
!  <DATA NAME="atmos_nthreads, ocean_nthreads" TYPE="integer">
!  We set here the number of OpenMP threads to use
!  separately for each component (default 1)
!  </DATA> 
!  <DATA NAME="use_lag_fluxes" TYPE="logical">
!  If true, then mom4 is forced with SBCs from one coupling timestep ago
!  If false, then mom4 is forced with most recent SBCs.
!  For a leapfrog MOM coupling with dt_cpld=dt_ocean, lag fluxes
!  can be shown to be stable and current fluxes to be unconditionally unstable.
!  For dt_cpld>dt_ocean there is probably sufficient damping.
!  use_lag_fluxes is set to TRUE by default.
!  </DATA>
!  <DATA NAME="restart_interval" TYPE="integer, dimension(6)"  DEFAULT="0">
!     The time interval that write out intermediate restart file. The format is (yr,mo,day,hr,min,sec).
!     When restart_interval is all zero, no intermediate restart file will be written out.
!   </DATA>
!   <NOTE>
!     <PRE>
!     1.If no value is set for current_date, start_date, or calendar (or default value specified) then the value from restart
!       file "INPUT/coupler.res" will be used. If neither a namelist value or restart file value exist the program will fail. 
!     2.The actual run length will be the sum of months, days, hours, minutes, and seconds. A run length of zero is not a
!       valid option. 
!     3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>
! </NAMELIST>

  integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/)
  integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '
  logical :: force_date_from_namelist = .false.  ! override restart values for date
  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & ice/land
  integer :: dt_cpld  = 0  ! fluxes passed between ice & ocean


  integer :: atmos_npes=0, ocean_npes=0, ice_npes=0, land_npes=0
  integer :: atmos_nthreads=1, ocean_nthreads=1
  logical :: do_atmos =.true., do_land =.true., do_ice =.true., do_ocean=.true.
  logical :: do_flux =.true.
  logical :: concurrent=.FALSE.
  logical :: use_lag_fluxes=.TRUE.
  logical :: do_chksum=.FALSE.
  integer :: check_stocks = 0 ! -1: never 0: at end of run only n>0: every n coupled steps

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, months, days, hours,      &
                         minutes, seconds, dt_cpld, dt_atmos, do_atmos,              &
                         do_land, do_ice, do_ocean, do_flux, atmos_npes, ocean_npes, &
                         ice_npes, land_npes, atmos_nthreads, ocean_nthreads, &
                         concurrent, use_lag_fluxes, do_chksum, &
                         check_stocks, restart_interval

  integer :: initClock, mainClock, termClock

  integer :: newClock0, newClock1, newClock2, newClock3, newClock4, newClock5, newClock6, newClock7, &
             newClock8, newClock9, newClock10, newClock11, newClock12, newClock13, newClock14, newClocka, &
             newClockb, newClockc, newClockd, newClocke, newClockf, newClockg, newClockh

  integer :: id_atmos_model_init, id_land_model_init, id_ice_model_init
  integer :: id_ocean_model_init, id_flux_exchange_init

  character(len=80) :: text
  character(len=48), parameter                    :: mod_name = 'coupler_main_mod'
 
  integer :: ensemble_id = 1 , outunit
  integer, allocatable :: ensemble_pelist(:, :) 

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_main'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!#######################################################################

  call mpp_init()
!these clocks are on the global pelist
  initClock = mpp_clock_id( 'Initialization' )
  call mpp_clock_begin(initClock)
  
  call fms_init
  call constants_init

  call coupler_init
  if(do_chksum) call coupler_chksum('coupler_init+', 0)

  call mpp_set_current_pelist()

  call mpp_clock_end (initClock) !end initialization

  call mpp_clock_begin(mainClock) !begin main loop

!-----------------------------------------------------------------------
!------ ocean/slow-ice integration loop ------

     if(check_stocks >= 0) then
        call mpp_set_current_pelist()
        call flux_init_stocks(Time, Atm, Land, Ice, Ocean_state)
     endif

if( Atm%pe )then
 call mpp_set_current_pelist(Atm%pelist)
 newClock1 = mpp_clock_id( 'generate_sfc_xgrid' )
endif
call mpp_set_current_pelist()
newClock2 = mpp_clock_id( 'flux_ocean_to_ice' )
newClock3 = mpp_clock_id( 'flux_ice_to_ocean' )
newClock4 = mpp_clock_id( 'flux_check_stocks' ) 
if( Atm%pe )then
 call mpp_set_current_pelist(Atm%pelist)
 newClock5 = mpp_clock_id( 'ATM' )
 newClock6  = mpp_clock_id( '  ATM: update_ice_model_slow_up' )
 newClock7  = mpp_clock_id( '  ATM: atmos loop' )
 newClocka  = mpp_clock_id( '     A-L: atmos_tracer_driver_gather_data' )
 newClockb  = mpp_clock_id( '     A-L: sfc_boundary_layer' )
 newClockc  = mpp_clock_id( '     A-L: update_atmos_model_down' )
 newClockd  = mpp_clock_id( '     A-L: flux_down_from_atmos' )
 newClocke  = mpp_clock_id( '     A-L: update_land_model_fast' )
 newClockf  = mpp_clock_id( '     A-L: update_ice_model_fast' )
 newClockg  = mpp_clock_id( '     A-L: flux_up_to_atmos' )
 newClockh  = mpp_clock_id( '     A-L: update_atmos_model_up' )
 newClock8  = mpp_clock_id( '  ATM: update_land_model_slow' )
 newClock9  = mpp_clock_id( '  ATM: flux_land_to_ice' )
 newClock10 = mpp_clock_id( '  ATM: update_ice_model_slow_dn' )
 newClock11 = mpp_clock_id( '  ATM: flux_ice_to_ocean_stocks' )
endif
if( Ocean%is_ocean_pe )then
 call mpp_set_current_pelist(Ocean%pelist)
 newClock12 = mpp_clock_id( 'OCN' )
endif
call mpp_set_current_pelist()
newClock13 = mpp_clock_id( 'intermediate restart' )
newClock14 = mpp_clock_id( 'final flux_check_stocks' )

  do nc = 1, num_cpld_calls
     if(do_chksum) call coupler_chksum('top_of_coupled_loop+', nc)
     if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        call mpp_clock_begin(newClock1)
        call generate_sfc_xgrid( Land, Ice )
        call mpp_clock_end(newClock1)
     end if
     call mpp_set_current_pelist()
     if(do_chksum) then
       if (Atm%pe) then 
         call mpp_set_current_pelist(Atm%pelist)
         call atmos_ice_land_chksum('MAIN_LOOP-', nc)
       endif
       if (Ocean%is_ocean_pe) then 
         call mpp_set_current_pelist(Ocean%pelist)
         call ocean_chksum('MAIN_LOOP-', nc)
       endif
       call mpp_set_current_pelist()
     endif  

     ! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
     ! points when running concurrently. The calls are placed next to each other in
     ! concurrent mode to avoid multiple synchronizations within the main loop.
     ! This is only possible in the serial case when use_lag_fluxes.
     call mpp_clock_begin(newClock2)
     call flux_ocean_to_ice( Time, Ocean, Ice, Ocean_ice_boundary )
     call mpp_clock_end(newClock2)
     if(do_chksum) then
       call coupler_chksum('flux_ocn2ice+', nc)
       if (Atm%pe) then 
         call mpp_set_current_pelist(Atm%pelist)
         call atmos_ice_land_chksum('fluxocn2ice+', nc)
       endif
       if (Ocean%is_ocean_pe) then 
         call mpp_set_current_pelist(Ocean%pelist)
         call ocean_public_type_chksum('fluxocn2ice+', nc, Ocean)
       endif
       call mpp_set_current_pelist()
     endif

     ! Update Ice_ocean_boundary; first iteration is supplied by restart     
     if( use_lag_fluxes )then
        call mpp_clock_begin(newClock3)
        call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
        call mpp_clock_end(newClock3)
     end if

     ! Update Ice_ocean_boundary; first iteration is supplied by restart     
     ! To print the value of frazil heat flux at the right time the following block
     ! needs to sit here rather than at the end of the coupler loop.
     if(check_stocks > 0) then
        call mpp_clock_begin(newClock4)
        if(check_stocks*((nc-1)/check_stocks) == nc-1 .AND. nc > 1) then
           call mpp_set_current_pelist()
           call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
        endif
        call mpp_clock_end(newClock4)
     endif

     if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        call mpp_clock_begin(newClock5)
        call mpp_clock_begin(newClock6)
        if (do_ice .AND. Ice%pe) then
           if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Ice%pelist)
           call update_ice_model_slow_up( Ocean_ice_boundary, Ice )
        endif
        if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
        call mpp_clock_end(newClock6)
        if(do_chksum) call atmos_ice_land_chksum('update_ice_slow_up+', nc)

        !-----------------------------------------------------------------------
        !   ------ atmos/fast-land/fast-ice integration loop -------

        call mpp_clock_begin(newClock7)
        do na = 1, num_atmos_calls
           if(do_chksum) call atmos_ice_land_chksum('top_of_atmos_loop-', (nc-1)*num_atmos_calls+na)

           Time_atmos = Time_atmos + Time_step_atmos

           if (do_atmos) then
              call mpp_clock_begin(newClocka)
              call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
              call mpp_clock_end(newClocka)
           endif

           if (do_flux) then
              call mpp_clock_begin(newClockb)
              call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, &
                   Atm, Land, Ice, Land_ice_atmos_boundary )
              if(do_chksum)  call atmos_ice_land_chksum('sfc+', (nc-1)*num_atmos_calls+na)
              call mpp_clock_end(newClockb)
           end if

           !      ---- atmosphere down ----

           if (do_atmos) then
              call mpp_clock_begin(newClockc)
              call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
              call mpp_clock_end(newClockc)
           endif
           if(do_chksum) call atmos_ice_land_chksum('update_atmos_down+', (nc-1)*num_atmos_calls+na)

           call mpp_clock_begin(newClockd)
           call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
                Land_ice_atmos_boundary, &
                Atmos_land_boundary, &
                Atmos_ice_boundary )
           call mpp_clock_end(newClockd)
           if(do_chksum) call atmos_ice_land_chksum('flux_down_from_atmos+', (nc-1)*num_atmos_calls+na)

           !      --------------------------------------------------------------

           !      ---- land model ----

           call mpp_clock_begin(newClocke)
           if (do_land .AND. land%pe) then
              if(land_npes .NE. atmos_npes) call mpp_set_current_pelist(Land%pelist)  
              call update_land_model_fast( Atmos_land_boundary, Land )
           endif
           if(land_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
           call mpp_clock_end(newClocke)
           if(do_chksum) call atmos_ice_land_chksum('update_land_fast+', (nc-1)*num_atmos_calls+na)

           !      ---- ice model ----
           call mpp_clock_begin(newClockf)
           if (do_ice .AND. Ice%pe) then
              if(ice_npes .NE. atmos_npes)call mpp_set_current_pelist(Ice%pelist)
              call update_ice_model_fast( Atmos_ice_boundary, Ice )
           endif
           if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
           call mpp_clock_end(newClockf)
           if(do_chksum) call atmos_ice_land_chksum('update_ice_fast+', (nc-1)*num_atmos_calls+na)

           !      --------------------------------------------------------------
           !      ---- atmosphere up ----

           call mpp_clock_begin(newClockg)
           call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary, &
                & Atmos_land_boundary, Atmos_ice_boundary )
           call mpp_clock_end(newClockg)
           if(do_chksum) call atmos_ice_land_chksum('flux_up2atmos+', (nc-1)*num_atmos_calls+na)

           call mpp_clock_begin(newClockh)
           if (do_atmos) &
                call update_atmos_model_up( Land_ice_atmos_boundary, Atm )
           call mpp_clock_end(newClockh)
           if(do_chksum) call atmos_ice_land_chksum('update_atmos_up+', (nc-1)*num_atmos_calls+na)

           !--------------

        enddo
        call mpp_clock_end(newClock7)

        call mpp_clock_begin(newClock8)
        !   ------ end of atmospheric time step loop -----
        if (do_land .AND. Land%pe) then
           if(land_npes .NE. atmos_npes) call mpp_set_current_pelist(Land%pelist)
           call update_land_model_slow(Atmos_land_boundary,Land)
        endif
        if(land_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
        !-----------------------------------------------------------------------
        call mpp_clock_end(newClock8)
        if(do_chksum) call atmos_ice_land_chksum('update_land_slow+', nc)

        !
        !     need flux call to put runoff and p_surf on ice grid
        !
        call mpp_clock_begin(newClock9)
        call flux_land_to_ice( Time, Land, Ice, Land_ice_boundary )
        call mpp_clock_end(newClock9)
        if(do_chksum) call atmos_ice_land_chksum('fluxlnd2ice+', nc)

        Atmos_ice_boundary%p = 0.0 ! call flux_atmos_to_ice_slow ?

        !   ------ slow-ice model ------

        if (do_ice) then 
           call mpp_clock_begin(newClock10)
           if( Ice%pe ) then
              if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Ice%pelist)
              call update_ice_model_slow_dn( Atmos_ice_boundary, &
                   & Land_ice_boundary, Ice )
           endif
           if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
           call mpp_clock_end(newClock10)
           if(do_chksum) call atmos_ice_land_chksum('update_ice_slow_dn+', nc)

           call mpp_clock_begin(newClock11)
           if( Ice%pe ) then
              if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Ice%pelist)
              call flux_ice_to_ocean_stocks(Ice)
           endif
           if(ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
           call mpp_clock_end(newClock11)
           if(do_chksum) call atmos_ice_land_chksum('fluxice2ocn_stocks+', nc)
        endif
        Time = Time_atmos
        call mpp_clock_end(newClock5)
     end if                     !Atm%pe block

     if( .NOT.use_lag_fluxes )then !this will serialize
        call mpp_set_current_pelist()
        call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
     end if

     if( Ocean%is_ocean_pe )then
        call mpp_set_current_pelist(Ocean%pelist)
        call mpp_clock_begin(newClock12)

        if (do_chksum) call ocean_chksum('update_ocean_model-', nc)
        ! update_ocean_model since fluxes don't change here

        if (do_ocean) &
          call update_ocean_model( Ice_ocean_boundary, Ocean_state,  Ocean, &
                                   Time_ocean, Time_step_cpld )

        if (do_chksum) call ocean_chksum('update_ocean_model+', nc)
        ! Get stocks from "Ice_ocean_boundary" and add them to Ocean stocks.
        ! This call is just for record keeping of stocks transfer and
        ! does not modify either Ocean or Ice_ocean_boundary
        call flux_ocean_from_ice_stocks(Ocean_state, Ocean, Ice_ocean_boundary)

        Time_ocean = Time_ocean +  Time_step_cpld

        !-----------------------------------------------------------------------
        Time = Time_ocean

        call mpp_clock_end(newClock12)
     end if

!rabcall mpp_clock_begin(newClock13)
     !--- write out intermediate restart file when needed.
     if( Time >= Time_restart ) then
        Time_restart_current = Time
        Time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
             restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
        timestamp = date_to_string(time_restart_current)
        outunit= stdout()
        write(outunit,*) '=> NOTE from program coupler: intermediate restart file is written and ', &
             trim(timestamp),' is appended as prefix to each restart file name'
        if( Atm%pe )then        
           call atmos_model_restart(Atm, timestamp)
           call land_model_restart(timestamp)
           call ice_model_restart(Ice, timestamp)
        endif
        if( Ocean%is_ocean_pe) then
           call ocean_model_restart(Ocean_state, timestamp)
        endif
        call coupler_restart(Time, Time_restart_current, timestamp)
     end if

     !--------------
     if(do_chksum) call coupler_chksum('MAIN_LOOP+', nc)
     write( text,'(a,i4)' )'Main loop at coupling timestep=', nc
     call print_memuse_stats(text)
!rabcall mpp_clock_end(newClock13)


  enddo

     call mpp_set_current_pelist()
call mpp_clock_begin(newClock14)
  if(check_stocks >= 0) then
     call mpp_set_current_pelist()
     call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
  endif
call mpp_clock_end(newClock14)

! Need final update of Ice_ocean_boundary for concurrent restart
!  if( concurrent )then
!      call mpp_set_current_pelist()
!      call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
!  endif

  call mpp_set_current_pelist()
!-----------------------------------------------------------------------
  call mpp_clock_end(mainClock)
  call mpp_clock_begin(termClock)

  if(do_chksum) call coupler_chksum('coupler_end-', nc)
  call coupler_end

  call mpp_clock_end(termClock)

  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine coupler_init

    use ensemble_manager_mod, only : ensemble_manager_init, get_ensemble_id,ensemble_pelist_setup
    use ensemble_manager_mod, only : get_ensemble_size, get_ensemble_pelist

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

    character(len=64), parameter    :: sub_name = 'coupler_init'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: warn_header =                                &
         '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: unit,  ierr, io,    m, i, outunit, logunit, errunit
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes

    integer :: ens_siz(6), ensemble_size

    integer :: atmos_pe_start=0, atmos_pe_end=0, &
               ocean_pe_start=0, ocean_pe_end=0
    integer :: n
    integer :: diag_model_subset=DIAG_ALL
    logical :: other_fields_exist
    character(len=256) :: err_msg
    integer :: date_restart(6)
    character(len=64)  :: filename, fieldname
    integer :: id_restart, l
    integer :: omp_get_thread_num, omp_get_num_threads
    integer :: get_cpu_affinity, base_cpu
    character(len=8)  :: walldate
    character(len=10) :: walltime
    character(len=5)  :: wallzone
    integer           :: wallvalues(8)
!-----------------------------------------------------------------------

    outunit = stdout()
    errunit = stderr()
    logunit = stdlog()
    
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Entering coupler_init at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

!----- write version to logfile -------
    call write_version_number(version, tag)

!----- read namelist -------

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, coupler_nml, iostat=io)
    ierr = check_nml_error (io, 'coupler_nml')
#else
    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call mpp_close(unit)
#endif

    if(ocean_nthreads /= 1) then
       call error_mesg ('coupler_init', 'OpenMP threading is not currently available for the ocean, '// &
            'Resetting variable ocean_nthreads to 1.', NOTE )
        ocean_nthreads = 1
    endif


!---- when concurrent is set true and mpp_io_nml io_clock_on is set true, the model
!---- will crash with error message "MPP_CLOCK_BEGIN: cannot change pelist context of a clock",
!---- so need to make sure it will not happen
    if(concurrent) then
       if(mpp_io_clock_on()) then
          call error_mesg ('program coupler', 'when coupler_nml variable concurrent is set to true, '// &
              'mpp_io_nml variable io_clock_non can not be set to true.', FATAL )
       endif
    endif
!----- read date and calendar type from restart file -----

    if( file_exist('INPUT/coupler.res') )then
!Balaji: currently written in binary, needs form=MPP_NATIVE
        call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY )
        read( unit,*,err=999 )calendar_type
        read( unit,* )date_init
        read( unit,* )date
        goto 998 !back to fortran-4
!read old-style coupler.res
999     call mpp_close(unit)
        call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY, form=MPP_NATIVE )
        read(unit)calendar_type
        read(unit)date
998     call mpp_close(unit)
    else
        force_date_from_namelist = .true.
    endif

!----- use namelist value (either no restart or override flag on) ---

    if ( force_date_from_namelist ) then

        if ( sum(current_date) <= 0 ) then
            call error_mesg ('program coupler',  &
                 'no namelist value for base_date or current_date', FATAL)
        else
            date      = current_date
        endif

!----- override calendar type with namelist value -----

        select case( uppercase(trim(calendar)) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        end select

    endif

    call set_calendar_type (calendar_type, err_msg)
    if(err_msg /= '') then
      call mpp_error(FATAL, 'ERROR in coupler_init: '//trim(err_msg))
    endif

    if( concurrent .AND. .NOT.use_lag_fluxes )call mpp_error( WARNING, &
            'coupler_init: you have set concurrent=TRUE and use_lag_fluxes=FALSE &
            & in coupler_nml. When not using lag fluxes, components &
            & will synchronize at two points, and thus run serially.' )


    !Check with the ensemble_manager module for the size of ensemble
    !and PE counts for each member of the ensemble.
    !
    !NOTE: ensemble_manager_init renames all the output files (restart and diagnostics)
    !      to show which ensemble member they are coming from.
    !      There also need to be restart files for each member of the ensemble in INPUT.
    !
    !NOTE: if the ensemble_size=1 the input/output files will not be renamed.
    !
    
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting initializing ensemble_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call ensemble_manager_init() ! init pelists for ensembles
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing ensemble_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    ens_siz = get_ensemble_size()   
    ensemble_size = ens_siz(1)      
    npes = ens_siz(2)              

    !Check for the consistency of PE counts
    if( concurrent )then
!atmos_npes + ocean_npes must equal npes
        if( atmos_npes.EQ.0 )atmos_npes = npes - ocean_npes
        if( ocean_npes.EQ.0 )ocean_npes = npes - atmos_npes
!both must now be non-zero
        if( atmos_npes.EQ.0 .OR. ocean_npes.EQ.0 ) &
             call mpp_error( FATAL, 'coupler_init: atmos_npes or ocean_npes must be specified for concurrent coupling.' )
        if( atmos_npes+ocean_npes.NE.npes ) &
             call mpp_error( FATAL, 'coupler_init: atmos_npes+ocean_npes must equal npes for concurrent coupling.' )
    else                        !serial timestepping
        if( (atmos_npes.EQ.0) .and. (do_atmos .or. do_land .or. do_ice) ) atmos_npes = npes
        if( (ocean_npes.EQ.0) .and. (do_ocean) ) ocean_npes = npes
        if( max(atmos_npes,ocean_npes).EQ.npes )then !overlapping pelists
            ! do nothing
        else                    !disjoint pelists
            if( atmos_npes+ocean_npes.NE.npes ) call mpp_error( FATAL,  &
                 'coupler_init: atmos_npes+ocean_npes must equal npes for serial coupling on disjoint pelists.' )
        end if
    end if    

    if( land_npes == 0 ) land_npes = atmos_npes
    if( ice_npes  == 0 ) ice_npes  = atmos_npes    
    if(land_npes > atmos_npes) call mpp_error(FATAL, 'coupler_init: land_npes > atmos_npes')
    if(ice_npes  > atmos_npes) call mpp_error(FATAL, 'coupler_init: ice_npes > atmos_npes')

    allocate( Atm%pelist  (atmos_npes) )
    allocate( Ocean%pelist(ocean_npes) )
    allocate( Land%pelist (land_npes) )
    allocate( Ice%pelist  (ice_npes) )

    !Set up and declare all the needed pelists
    call ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
                               Atm%pelist, Ocean%pelist, Land%pelist, Ice%pelist)

!set up affinities based on threads

    ensemble_id = get_ensemble_id() 
 
    allocate(ensemble_pelist(1:ensemble_size,1:npes))   
    call get_ensemble_pelist(ensemble_pelist) 

    Atm%pe            = ANY(Atm%pelist   .EQ. mpp_pe()) 
    Ocean%is_ocean_pe = ANY(Ocean%pelist .EQ. mpp_pe())  
    Ice%pe            = ANY(Ice%pelist   .EQ. mpp_pe())  
    Land%pe           = ANY(Land%pelist  .EQ. mpp_pe()) 
 
    !Why is the following needed?
    if( Atm%pe )then
!$      call omp_set_num_threads(atmos_nthreads)
        call mpp_set_current_pelist( Atm%pelist )
!$      base_cpu = get_cpu_affinity()
!$OMP PARALLEL
!$        call set_cpu_affinity( base_cpu + omp_get_thread_num() )
!$OMP END PARALLEL
    end if

    if( Ocean%is_ocean_pe )then
       call mpp_set_current_pelist( Ocean%pelist )
!           call omp_set_num_threads(ocean_nthreads)
!           base_cpu = get_cpu_affinity()
!   !$OMP PARALLEL
!           call set_cpu_affinity( base_cpu + omp_get_thread_num() )
!   !$OMP END PARALLEL
       end if

   !--- initialization clock
    if( Atm%pe )then
       call mpp_set_current_pelist(Atm%pelist)
       id_atmos_model_init = mpp_clock_id( '  Init: atmos_model_init ' )
    endif
    if( Land%pe )then
       call mpp_set_current_pelist(Land%pelist)
       id_land_model_init  = mpp_clock_id( '  Init: land_model_init ' )
    endif
    if( Ice%pe )then
       call mpp_set_current_pelist(Ice%pelist)
       id_ice_model_init   = mpp_clock_id( '  Init: ice_model_init ' )
    endif
    if( Ocean%is_ocean_pe )then
       call mpp_set_current_pelist(Ocean%pelist)
       id_ocean_model_init = mpp_clock_id( '  Init: ocean_model_init ' )
    endif
    call mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    id_flux_exchange_init = mpp_clock_id( '  Init: flux_exchange_init' )

    call mpp_set_current_pelist()
    mainClock = mpp_clock_id( 'Main loop' )
    termClock = mpp_clock_id( 'Termination' )
    
    !Write out messages on root PEs
    if(mpp_pe().EQ.mpp_root_pe() )then
       write( text,'(a,2i6,a,i2.2)' )'Atmos PE range: ', Atm%pelist(1)  , Atm%pelist(atmos_npes)  ,&
            ' ens_', ensemble_id
       call mpp_error( NOTE, 'coupler_init: '//trim(text) )
       if (ocean_npes .gt. 0) then   ! only if ocean is active (cjg)
         write( text,'(a,2i6,a,i2.2)' )'Ocean PE range: ', Ocean%pelist(1), Ocean%pelist(ocean_npes), &
              ' ens_', ensemble_id
         call mpp_error( NOTE, 'coupler_init: '//trim(text) )
       else
         write( text,'(a,i2.2)' )'Ocean PE range is not set (do_ocean=.false. and concurrent=.false.) for ens_', &
               ensemble_id
         call mpp_error( NOTE, 'coupler_init: '//trim(text) )
       end if
       write( text,'(a,2i6,a,i2.2)' )'Land PE range: ', Land%pelist(1)  , Land%pelist(land_npes)  ,&
            ' ens_', ensemble_id
       call mpp_error( NOTE, 'coupler_init: '//trim(text) )
       write( text,'(a,2i6,a,i2.2)' )'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes), &
            ' ens_', ensemble_id
       call mpp_error( NOTE, 'coupler_init: '//trim(text) )

       if( concurrent )then
          call mpp_error( NOTE, 'coupler_init: Running with CONCURRENT coupling.' )

          write( logunit,'(a)' )'Using concurrent coupling...'
          write( logunit,'(a,4i4)' ) &
               'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
               Atm%pelist(1)  , Atm%pelist(atmos_npes), Ocean%pelist(1), Ocean%pelist(ocean_npes) 
       else
          call mpp_error( NOTE, 'coupler_init: Running with SERIAL coupling.' )
       end if
       if( use_lag_fluxes )then
          call mpp_error( NOTE, 'coupler_init: Sending LAG fluxes to ocean.' )
       else
          call mpp_error( NOTE, 'coupler_init: Sending most recent fluxes to ocean.' )
       end if
    endif

!----- write namelist to logfile -----
    if( mpp_pe() == mpp_root_pe() )write( logunit, nml=coupler_nml )

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe().EQ.mpp_root_pe() ) &
         write( logunit, 16 )date(1),trim(month_name(date(2))),date(3:6)
16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

!jwd Fork here is somewhat dangerous. It relies on "no side effects" from
!    diag_manager_init. diag_manager_init or this section should be 
!    re-architected to guarantee this or remove this assumption.
!    For instance, what follows assumes that get_base_date has the same
!    time for both Atm and Ocean pes. While this should be the case, the
!    possible error condition needs to be checked

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        if(atmos_npes /= npes)diag_model_subset = DIAG_OTHER  ! change diag_model_subset from DIAG_ALL
    elseif( Ocean%is_ocean_pe )then  ! Error check above for disjoint pelists should catch any problem
        call mpp_set_current_pelist(Ocean%pelist)
        if(ocean_npes /= npes)diag_model_subset = DIAG_OCEAN  ! change diag_model_subset from DIAG_ALL
    end if
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize diag_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call diag_manager_init(DIAG_MODEL_SUBSET=diag_model_subset)   ! initialize diag_manager for processor subset output
    call print_memuse_stats( 'diag_manager_init' )
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing diag_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
!-----------------------------------------------------------------------
!------ reset pelist to "full group" ------

    call mpp_set_current_pelist()
!----- always override initial/base date with diag_manager value -----

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6)  )

!----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------

    Time_init = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    Time      = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

    Time_start = Time

!----- compute the ending time -----

    Time_end = Time
    do m=1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do
    Time_end   = Time_end + set_time(hours*3600+minutes*60+seconds, days)
    !Need to pass Time_end into diag_manager for multiple thread case.
    call diag_manager_set_time_end(Time_end)

    Run_length = Time_end - Time

!--- get the time that last intermediate restart file was written out.
    if (file_exist('INPUT/coupler.intermediate.res')) then
       call mpp_open(unit,'INPUT/coupler.intermediate.res',action=MPP_RDONLY)
       read(unit,*) date_restart
       call mpp_close(unit)
    else
       date_restart = date
    endif

    Time_restart_current = Time
    if(ALL(restart_interval ==0)) then
       Time_restart = increment_date(Time_end, 0, 0, 10, 0, 0, 0)   ! no intermediate restart
    else
       Time_restart = set_date(date_restart(1), date_restart(2), date_restart(3),  &
                               date_restart(4), date_restart(5), date_restart(6) )
       Time_restart = increment_date(Time_restart, restart_interval(1), restart_interval(2), &
            restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
       if(Time_restart <= Time) call mpp_error(FATAL, &
            '==>Error from program coupler: The first intermediate restart time is no larger than the start time')
    end if

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

    call mpp_open( unit, 'time_stamp.out', nohdrs=.TRUE. )

    month = month_name(date(2))
    if ( mpp_pe().EQ.mpp_root_pe() ) write (unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
         date(4), date(5), date(6))
    month = month_name(date(2))
    if ( mpp_pe().EQ.mpp_root_pe() ) write (unit,20) date, month(1:3)

    call mpp_close(unit)

20  format (6i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step_cpld  = set_time (dt_cpld ,0)
    Time_step_atmos = set_time (dt_atmos,0)

!----- determine maximum number of iterations per loop ------

    num_cpld_calls  = Run_length      / Time_step_cpld
    num_atmos_calls = Time_step_cpld  / Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time ) call error_mesg ('program coupler',  &
         'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * Time_step_cpld  /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of coupled time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the atmos time step', FATAL)

!
!       Initialize the tracer manager. This needs to be done on all PEs,
!       before the individual models are initialized.
!

    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize tracer_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call tracer_manager_init
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing tracer_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
!
!       Initialize the coupler types
!

    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize coupler_types at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call coupler_types_init
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing coupler_types at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

!-----------------------------------------------------------------------
!------ initialize component models ------
!------ grid info now comes from grid_spec file

    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Beginning to initialize component models at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
!---- atmosphere ----
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize atmospheric model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif

        call mpp_clock_begin(id_atmos_model_init)
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos )
        call mpp_clock_end(id_atmos_model_init)
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing atmospheric model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call print_memuse_stats( 'atmos_model_init' )
        call data_override_init(Atm_domain_in = Atm%domain)
     endif
!---- land ----------
     if( Land%pe ) then
        call mpp_set_current_pelist(Land%pelist)
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize land model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call mpp_clock_begin(id_land_model_init)
        call land_model_init( Atmos_land_boundary, Land, Time_init, Time, &
             Time_step_atmos, Time_step_cpld )
        call mpp_clock_end(id_land_model_init)
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing land model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call print_memuse_stats( 'land_model_init' )
        call data_override_init(Land_domain_in = Land%domain)
     endif
!---- ice -----------
     if( Ice%pe ) then
        call mpp_set_current_pelist(Ice%pelist)
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize ice model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call mpp_clock_begin(id_ice_model_init)
        call ice_model_init( Ice, Time_init, Time, Time_step_atmos, Time_step_cpld )
        call mpp_clock_end(id_ice_model_init)
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing ice model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call print_memuse_stats( 'ice_model_init' )
        call data_override_init(Ice_domain_in = Ice%domain)
    end if
    if( Ocean%is_ocean_pe )then
        call mpp_set_current_pelist(Ocean%pelist)
!---- ocean ---------
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize ocean model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call mpp_clock_begin(id_ocean_model_init)
        call ocean_model_init( Ocean, Ocean_state, Time_init, Time )
        call mpp_clock_end(id_ocean_model_init)
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing ocean model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call print_memuse_stats( 'ocean_model_init' )
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize data_override at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call data_override_init(Ocean_domain_in = Ocean%domain )
        if( mpp_pe().EQ.mpp_root_pe() ) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing data_override at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
    end if
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing component models at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))

    call mpp_broadcast_domain(Ice%domain)
    call mpp_broadcast_domain(Ocean%domain)
!-----------------------------------------------------------------------
!---- initialize flux exchange module ----
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize flux_exchange at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call mpp_clock_begin(id_flux_exchange_init)
    call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, Ocean_state,&
         atmos_ice_boundary, land_ice_atmos_boundary, &
         land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
         dt_atmos=dt_atmos, dt_cpld=dt_cpld)
    call mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    call mpp_clock_end(id_flux_exchange_init)
    call mpp_set_current_pelist()
    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finsihed initializing flux_exchange at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

    Time_atmos = Time
    Time_ocean = Time

!
!       read in extra fields for the air-sea gas fluxes
!

    if ( Ice%pe ) then
      call mpp_set_current_pelist(Ice%pelist)
      allocate(Ice_bc_restart(Ice%ocean_fluxes%num_bcs))
      allocate(ice_bc_restart_file(Ice%ocean_fluxes%num_bcs))
      do n = 1, Ice%ocean_fluxes%num_bcs  !{
        if(Ice%ocean_fluxes%bc(n)%num_fields .LE. 0) cycle
        filename = trim(Ice%ocean_fluxes%bc(n)%ice_restart_file)
        do l = 1, num_ice_bc_restart
           if(trim(filename) == ice_bc_restart_file(l)) exit
        end do
        if(l>num_ice_bc_restart) then
           num_ice_bc_restart = num_ice_bc_restart + 1
           ice_bc_restart_file(l) = trim(filename)
        end if
        filename = 'INPUT/'//trim(filename)
        other_fields_exist = .false.
        do m = 1, Ice%ocean_fluxes%bc(n)%num_fields  !{
          fieldname = trim(Ice%ocean_fluxes%bc(n)%field(m)%name)
          id_restart = register_restart_field(Ice_bc_restart(l), ice_bc_restart_file(l), &
                       fieldname, Ice%ocean_fluxes%bc(n)%field(m)%values, Ice%domain    )
          if (field_exist(filename, fieldname, Ice%domain) ) then
            other_fields_exist = .true.
            write (outunit,*) trim(note_header), ' Reading restart info for ',         &
                 trim(fieldname), ' from ',  trim(filename)
            call read_data(filename, fieldname, Ice%ocean_fluxes%bc(n)%field(m)%values, Ice%domain)
          elseif (other_fields_exist) then
            call mpp_error(FATAL, trim(error_header) // ' Couldn''t find field ' //     &
                 trim(fieldname) // ' in file ' //trim(filename))
          endif
        enddo  !} m
      enddo  !} n
    endif
    if ( Ocean%is_ocean_pe ) then
      call mpp_set_current_pelist(Ocean%pelist)
      allocate(Ocn_bc_restart(Ocean%fields%num_bcs))
      allocate(ocn_bc_restart_file(Ocean%fields%num_bcs))
      do n = 1, Ocean%fields%num_bcs  !{
        if(Ocean%fields%bc(n)%num_fields .LE. 0) cycle
        filename = trim(Ocean%fields%bc(n)%ocean_restart_file)
        do l = 1, num_ocn_bc_restart
           if(trim(filename) == ocn_bc_restart_file(l)) exit
        end do
        if(l>num_ocn_bc_restart) then
           num_ocn_bc_restart = num_ocn_bc_restart + 1
           ocn_bc_restart_file(l) = trim(filename)
        end if
        filename = 'INPUT/'//trim(filename)
        other_fields_exist = .false.
        do m = 1, Ocean%fields%bc(n)%num_fields  !{
          fieldname = trim(Ocean%fields%bc(n)%field(m)%name)
          id_restart = register_restart_field(Ocn_bc_restart(l), Ocn_bc_restart_file(l), &
                       fieldname, Ocean%fields%bc(n)%field(m)%values, Ocean%domain    )
          if (field_exist(filename, fieldname, Ocean%domain) ) then
            other_fields_exist = .true.
            write (outunit,*) trim(note_header), ' Reading restart info for ',         &
                 trim(fieldname), ' from ', trim(filename)
            call read_data(filename, fieldname, Ocean%fields%bc(n)%field(m)%values, Ocean%domain)
          elseif (other_fields_exist) then
            call mpp_error(FATAL, trim(error_header) // ' Couldn''t find field ' //     &
                 trim(fieldname) // ' in file ' //trim(filename))
          endif
        enddo  !} m
      enddo  !} n
    endif

    call mpp_set_current_pelist()

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

    call mpp_open( unit, 'RESTART/file' )
    call mpp_close(unit, MPP_DELETE)

    ! Call to daig_grid_end to free up memory used during regional
    ! output setup
    CALL diag_grid_end()

!-----------------------------------------------------------------------
    if (Atm%pe) then 
      call mpp_set_current_pelist(Atm%pelist)
      call atmos_ice_land_chksum('coupler_init+', 0)
    endif
    if (Ocean%is_ocean_pe) then 
      call mpp_set_current_pelist(Ocean%pelist)
      call ocean_chksum('coupler_init+', nc)
    endif
    call mpp_set_current_pelist()
    call print_memuse_stats('coupler_init')

    if( mpp_pe().EQ.mpp_root_pe() ) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Exiting coupler_init at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
  end subroutine coupler_init

!#######################################################################

  subroutine coupler_end

!-----------------------------------------------------------------------

    if (Atm%pe) then 
      call mpp_set_current_pelist(Atm%pelist)
      call atmos_ice_land_chksum('coupler_end', 0)
    endif
    if (Ocean%is_ocean_pe) then 
      call mpp_set_current_pelist(Ocean%pelist)
      call ocean_chksum('coupler_end', 0)
    endif
    call mpp_set_current_pelist()

!----- check time versus expected ending time ----

    if (Time /= Time_end) call error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

!-----------------------------------------------------------------------
!the call to fms_io_exit has been moved here
!this will work for serial code or concurrent (disjoint pelists)
!but will fail on overlapping but unequal pelists
    if( Ocean%is_ocean_pe )then
        call mpp_set_current_pelist(Ocean%pelist)
        call ocean_model_end (Ocean, Ocean_state, Time)
    end if
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        call atmos_model_end (Atm)
    endif
    if( Land%pe ) then
        call mpp_set_current_pelist(Land%pelist)
        call  land_model_end (Atmos_land_boundary, Land)
    endif
    if( Ice%pe ) then
        call mpp_set_current_pelist(Ice%pelist)
        call   ice_model_end (Ice)
    end if

    !----- write restart file ------
    call coupler_restart(Time, Time_restart_current)

    call fms_io_exit
    call diag_manager_end (Time)
    call mpp_set_current_pelist()

!-----------------------------------------------------------------------

  end subroutine coupler_end

  !--- writing restart file that contains running time and restart file writing time.
  subroutine coupler_restart(Time_run, Time_res, time_stamp)
    type(time_type),   intent(in)           :: Time_run, Time_res
    character(len=*), intent(in),  optional :: time_stamp
    character(len=128)                      :: file_run, file_res
    integer :: yr, mon, day, hr, min, sec, date(6), unit, n

    call mpp_set_current_pelist()

    ! write restart file
    if(present(time_stamp)) then
       file_run = 'RESTART/'//trim(time_stamp)//'.coupler.res'
       file_res = 'RESTART/'//trim(time_stamp)//'.coupler.intermediate.res'
    else
       file_run = 'RESTART/coupler.res'
       file_res = 'RESTART/coupler.intermediate.res'
    endif

    !----- compute current date ------
    call get_date (Time_run, date(1), date(2), date(3),  &
                   date(4), date(5), date(6))
    call mpp_open( unit, file_run, nohdrs=.TRUE. )
    if ( mpp_pe().EQ.mpp_root_pe() )then
        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        write( unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write( unit, '(6i6,8x,a)' )date, &
             'Current model time: year, month, day, hour, minute, second'
    end if
    call mpp_close(unit)

    if(Time_res > Time_start) then
       call mpp_open( unit, file_res, nohdrs=.TRUE. )
       if ( mpp_pe().EQ.mpp_root_pe() )then
          call get_date(Time_res ,yr,mon,day,hr,min,sec)
          write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
               'Current intermediate restart time: year, month, day, hour, minute, second'
       end if
       call mpp_close(unit)     
    end if

    if( Ocean%is_ocean_pe )then
        call mpp_set_current_pelist(Ocean%pelist)
        do n = 1, num_ocn_bc_restart
           call save_restart(Ocn_bc_restart(n), time_stamp)
        enddo
    endif
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        do n = 1, num_ice_bc_restart
           call save_restart(Ice_bc_restart(n), time_stamp)
        enddo
    endif


  end subroutine coupler_restart

!--------------------------------------------------------------------------

  subroutine coupler_chksum(id, timestep)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep

    type :: tracer_ind_type
       integer :: atm, ice, lnd ! indices of the tracer in the respective models
    end type tracer_ind_type
    integer                            :: n_atm_tr, n_lnd_tr, n_exch_tr
    integer                            :: n_atm_tr_tot, n_lnd_tr_tot
    integer                            :: i, tr, n, m, outunit
    type(tracer_ind_type), allocatable :: tr_table(:)
    character(32) :: tr_name

    call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
                             num_prog=n_atm_tr)
    call get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, &
                             num_prog=n_lnd_tr)

    ! assemble the table of tracer number translation by matching names of
    ! prognostic tracers in the atmosphere and surface models; skip all atmos.
    ! tracers that have no corresponding surface tracers.
    allocate(tr_table(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
       call get_tracer_names( MODEL_ATMOS, i, tr_name )
       tr_table(n)%atm = i
       tr_table(n)%ice = get_tracer_index ( MODEL_ICE,  tr_name )
       tr_table(n)%lnd = get_tracer_index ( MODEL_LAND, tr_name )
       if(tr_table(n)%ice/=NO_TRACER.or.tr_table(n)%lnd/=NO_TRACER) &
            n = n+1
    enddo
    n_exch_tr = n-1

100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

    if( Atm%pe )then
       call mpp_set_current_pelist(Atm%pelist)

       outunit = stdout()
       write(outunit,*) 'BEGIN CHECKSUM(Atm):: ', id, timestep
       write(outunit,100) 'atm%t_bot', mpp_chksum(atm%t_bot)
       write(outunit,100) 'atm%z_bot', mpp_chksum(atm%z_bot)
       write(outunit,100) 'atm%p_bot', mpp_chksum(atm%p_bot)
       write(outunit,100) 'atm%u_bot', mpp_chksum(atm%u_bot)
       write(outunit,100) 'atm%v_bot', mpp_chksum(atm%v_bot)
       write(outunit,100) 'atm%p_surf', mpp_chksum(atm%p_surf)
       write(outunit,100) 'atm%gust', mpp_chksum(atm%gust)
       do tr = 1,n_exch_tr
          n = tr_table(tr)%atm
          if(n /= NO_TRACER ) then
             call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
             write(outunit,100) 'atm%'//trim(tr_name), mpp_chksum(Atm%tr_bot(:,:,n))
          endif
       enddo
   
       write(outunit,100) 'land%t_surf', mpp_chksum(land%t_surf)
       write(outunit,100) 'land%t_ca', mpp_chksum(land%t_ca)
       write(outunit,100) 'land%rough_mom', mpp_chksum(land%rough_mom)
       write(outunit,100) 'land%rough_heat', mpp_chksum(land%rough_heat)
       write(outunit,100) 'land%rough_scale', mpp_chksum(land%rough_scale)
       do tr = 1,n_exch_tr
          n = tr_table(tr)%lnd
          if(n /= NO_TRACER ) then
             call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
             write(outunit,100) 'land%'//trim(tr_name), mpp_chksum(Land%tr(:,:,:,n))
          endif
       enddo
   
       write(outunit,100) 'ice%t_surf', mpp_chksum(ice%t_surf)
       write(outunit,100) 'ice%rough_mom', mpp_chksum(ice%rough_mom)
       write(outunit,100) 'ice%rough_heat', mpp_chksum(ice%rough_heat)
       write(outunit,100) 'ice%rough_moist', mpp_chksum(ice%rough_moist)
       write(outunit,*) 'STOP CHECKSUM(Atm):: ', id, timestep
   
    !endif

    !if( Ocean%is_ocean_pe )then
        !call mpp_set_current_pelist(Ocean%pelist)

       write(outunit,*) 'BEGIN CHECKSUM(Ice):: ', id, timestep
       do n = 1, ice%ocean_fields%num_bcs  !{
          do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
             !write(outunit,101) 'ice%', m, n, mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
             write(outunit,101) 'ice%',trim(ice%ocean_fields%bc(n)%name), &
                  trim(ice%ocean_fields%bc(n)%field(m)%name), mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
          enddo  !} m
       enddo  !} n
       write(outunit,*) 'STOP CHECKSUM(Ice):: ', id, timestep

    endif

    deallocate(tr_table)

    call mpp_set_current_pelist()

  end subroutine coupler_chksum

  !#######################################################################

  subroutine atmos_ice_land_chksum(id, timestep)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep

! This subroutine calls subroutine that will print out checksums of the elements 
! of the appropriate type. 
! For coupled models typically these types are not defined on all processors.
! It is assumed that the appropriate pelist has been set before entering this routine.
! This can be achieved in the following way.
!       if (Atm%pe) then 
!         call mpp_set_current_pelist(Atm%pelist)
!         call atmos_ice_land_chksum('MAIN_LOOP-', nc)
!       endif
! If you are on the global pelist before you enter this routine using the above call, 
! you can return to the global pelist by invoking
!       call mpp_set_current_pelist()
! after you exit. This is only necessary if you need to return to the global pelist.

     call atmos_data_type_chksum(     id, timestep, Atm)
     call lnd_ice_atm_bnd_type_chksum(id, timestep, Land_ice_atmos_boundary)

     if(Ice%pe) then
        call mpp_set_current_pelist(Ice%pelist)
        call ice_data_type_chksum(       id, timestep, Ice)
        call atm_ice_bnd_type_chksum(    id, timestep, Atmos_ice_boundary)
        call ocn_ice_bnd_type_chksum(    id, timestep, Ocean_ice_boundary)
     endif
     if(Land%pe) then
        call mpp_set_current_pelist(Land%pelist)
        call land_data_type_chksum(      id, timestep, Land)
        call atm_lnd_bnd_type_chksum(    id, timestep, Atmos_land_boundary)
     endif
        
     call mpp_set_current_pelist(Atm%pelist)

  end subroutine atmos_ice_land_chksum

  subroutine ocean_chksum(id, timestep)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep

! This subroutine calls subroutine that will print out checksums of the elements 
! of the appropriate type. 
! For coupled models typically these types are not defined on all processors.
! It is assumed that the appropriate pelist has been set before entering this routine.
! This can be achieved in the following way.
!       if (Ocean%is_ocean_pe) then 
!         call mpp_set_current_pelist(Ocean%pelist)
!         call ocean_chksum('MAIN_LOOP-', nc)
!       endif
! If you are on the global pelist before you enter this routine using the above call, 
! you can return to the global pelist by invoking
!       call mpp_set_current_pelist()
! after you exit. This is only necessary if you need to return to the global pelist.

        call ocean_public_type_chksum(id, timestep, Ocean)
        call ice_ocn_bnd_type_chksum( id, timestep, Ice_ocean_boundary)
        
  end subroutine ocean_chksum


  end program coupler_main

