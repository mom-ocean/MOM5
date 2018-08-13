#ifndef USE_ESMF_LIB
program main
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison
!</CONTACT>
!
!<CONTACT EMAIL="Dave.Bi@csiro.au"> Dave Bi (for OASIS3 hooks)
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh (for OASIS3 hooks)
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> V. Balaji
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Driver for ocean-only simulations and prototype setup for OASIS3 driver.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Driver for the ocean-only simulations. Similar to the FMS coupler, but 
! allows one to run the ocean model without compiling  other models. 
! Much simpler than the full FMS coupler. 
!
! This driver also provides the prototype hooks for using MOM with OASIS3,
! with this code surrounded by the cpp-preprocessor option "ifdef OASIS3".  
! The couping of MOM to OASIS3 has not been tested at GFDL. Rather, 
! CSIRO in Australia uses MOM with OASIS3, with Dave.Bi@csiro.au the primary 
! contact for questions regarding MOM and OASIS3. 
! </DESCRIPTION>
!
! <NAMELIST NAME="ocean_solo_nml">
!
!   <DATA NAME="date_init"  TYPE="integer, dimension(6)"  DEFAULT="0">
!     The date that the current integration starts with. If the restart file
!      ocean_solo.res is present, date_init will be taken from there.
!   </DATA>
!   <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
!     The calendar type used by the current integration. Valid values are consistent 
!     with the time_manager module: 'julian', 'gregorian', 'noleap', or 'thirty_day'. 
!     The value 'no_calendar' can not be used because the time_manager's date 
!     function are used.
!     
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
!   <DATA NAME="dt_cpld"  TYPE="integer"  DEFAULT="0">
!     Time step in seconds for coupling between ocean and atmospheric models: 
!     must be an integral multiple of dt_ocean. This is the "slow" timestep.
!     Note that for an ocean_solo model, the coupling to an "atmosphere" is the coupling 
!     to some data files.  In this case, dt_cpld represents the time which data is updated.
!     For example, if data is "daily", then dt_cpld=86400 should be chosen.  
!     If data is fixed, then dt_cpld of any integer of dt_ocean is fine, with
!     dt_cpld=86400 the default. 
!   </DATA>
!  <DATA NAME="n_mask" TYPE="integer">
!    number of region to be masked out. Its value should be less than MAX_PES.
!  </DATA>
!  <DATA NAME="mask_list(2,MAXPES)" TYPE="integer, dimension(MAX_MASK_REGION,2)">
!    The position of the region to be masked out. mask_list(1,:) is the x-layout position
!    and mask_list(2,:) is y-layout position.  
!  </DATA>
!  <DATA NAME="layout_mask" TYPE="integer, dimension(2)">
!   Processor domain layout for all the component model. layout_mask need to be set when and only 
!   when n_mask is greater than 0 ( some domain region is masked out ). When this namelist is set,
!   it will overload the layout in each component model. The default value is (0,0).
!   Currently we require all the component model has the same layout and same grid size.
!  </DATA>
!  <DATA NAME="restart_interval" TYPE="integer, dimension(6)"  DEFAULT="0">
!     The time interval that write out intermediate restart file. The format is (yr,mo,day,hr,min,sec).
!     When restart_interval is all zero, no intermediate restart file will be written out.
!   </DATA>
!
! </NAMELIST>
!
!   <NOTE>
!     <PRE>
!     1.The actual run length will be the sum of months, 
!       days, hours, minutes, and seconds. A run length of zero
!       is not a valid option. 
!     2.The run length must be an integral multiple of the coupling 
!       timestep dt_cpld. 
!     </PRE>
!   </NOTE>
!
  use constants_mod,            only: constants_init
  use data_override_mod,        only: data_override_init, data_override
  use diag_manager_mod,         only: diag_manager_init, diag_manager_end
  use field_manager_mod,        only: field_manager_init
  use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
  use fms_mod,                  only: close_file, file_exist, uppercase
  use fms_io_mod,               only: fms_io_exit
  use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain
  use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
  use mpp_mod,                  only: input_nml_file, mpp_error, FATAL, NOTE, mpp_pe, mpp_npes, mpp_set_current_pelist
  use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, MAXPES
  use time_interp_external_mod, only: time_interp_external_init
  use time_manager_mod,         only: set_calendar_type, time_type, increment_date
  use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
  use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )
  use time_manager_mod,         only: operator( * ), operator( /= ), operator( > )
  use time_manager_mod,         only: date_to_string

  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end
  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_types_mod,          only: ice_ocean_boundary_type

  implicit none

  type (ocean_public_type)               :: Ocean_sfc          
  type (ocean_state_type),       pointer :: Ocean_state => NULL()
  type(ice_ocean_boundary_type), target  :: Ice_ocean_boundary 

  ! define some time types 
  type(time_type) :: Time_init    ! initial time for experiment
  type(time_type) :: Time_start   ! start time for experiment
  type(time_type) :: Time_end     ! end time for experiment (as determined by dtts)
  type(time_type) :: Run_len      ! length of experiment 
  type(time_type) :: Time        
  type(time_type) :: Time_step_coupled
  type(time_type) :: Time_restart_init
  type(time_type) :: Time_restart
  type(time_type) :: Time_restart_current

  character(len=17) :: calendar = 'julian'

  integer :: dt_cpld  = 86400
  integer :: num_cpld_calls  = 0
  integer :: nc
  integer :: calendar_type=-1

  integer :: date_init(6)=0, date(6)
  integer :: date_restart(6)
  integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0

  integer :: isc,iec,jsc,jec
  integer :: unit, io_status, ierr

  integer :: flags=0, override_clock
  integer :: nfields 
  
  character(len=256) :: version = ''
  character(len=256) :: tag = ''

  character(len=9) :: month
  character(len=64):: timestamp

  integer :: n, m 
  integer :: layout_mask(2) = (/0 , 0/)
  integer :: n_mask = 0
  integer :: mask_list(2,MAXPES)
  integer, parameter :: mp = 2*MAXPES
  data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/
  integer :: restart_interval(6) = (/0,0,0,0,0,0/)
  integer :: mpi_comm_mom
  integer ::  stdoutunit,stdlogunit
  logical :: external_initialization

  namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds, dt_cpld, &
                            n_mask, layout_mask, mask_list, restart_interval

  call external_coupler_mpi_init(mpi_comm_mom, external_initialization)

  if ( external_initialization ) then
     call fms_init(mpi_comm_mom)
  else
     call fms_init()
  endif

  call constants_init

  flags = MPP_CLOCK_SYNC

  stdoutunit=stdout();stdlogunit=stdlog()

  ! provide for namelist over-ride
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_solo_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_solo_nml')
#else
  unit = open_namelist_file('input.nml')
  read  (unit, ocean_solo_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_solo_nml')
  call close_file (unit)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,'(/47x,a/)') '======== MODEL BEING DRIVEN BY OCEAN_SOLO_MOD ========'
  write (stdoutunit, ocean_solo_nml)  
  write (stdlogunit, ocean_solo_nml)

  write (stdlogunit,'(/,80("="),/(a))') trim(version), trim(tag)

  ! set the calendar 

  select case( uppercase(trim(calendar)) )
  case( 'GREGORIAN' )
     calendar_type = GREGORIAN
  case( 'JULIAN' )
     calendar_type = JULIAN
  case( 'NOLEAP' )
     calendar_type = NOLEAP
  case( 'THIRTY_DAY' )
     calendar_type = THIRTY_DAY_MONTHS
  case( 'NO_CALENDAR' )
     calendar_type = NO_CALENDAR
  case default
     call mpp_error(FATAL, &
     'ocean_solo: ocean_solo_nml entry calendar must be one of GREGORIAN|JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
  end select 

  ! get ocean_solo restart : this can override settings from namelist
  if (file_exist('INPUT/ocean_solo.res')) then
      call mpp_open(unit,'INPUT/ocean_solo.res',form=MPP_ASCII,action=MPP_RDONLY)
      read(unit,*) calendar_type 
      read(unit,*) date_init
      read(unit,*) date
      call mpp_close(unit)
  endif

  if (file_exist('INPUT/ocean_solo.intermediate.res')) then
      call mpp_open(unit,'INPUT/ocean_solo.intermediate.res',form=MPP_ASCII,action=MPP_RDONLY)
      read(unit,*) date_restart
      call mpp_close(unit)
  else
      date_restart = date_init
  endif
      
  call set_calendar_type (calendar_type)

!!$ initialize pelists for ocean ensembles set current pelist to ensemble member
!!$ need to call prior to diagnostics_init
!!$ code presently not supported (mjh)
!!$  call ocean_ensemble_init() 
                             
  call field_manager_init(nfields)

  call diag_manager_init()

  call time_interp_external_init()

  if (sum(date_init) <= 0) then
      call mpp_error(FATAL,&
      '==>Error from ocean_solo_mod: date_init must be set either in ocean_solo.res or in ocean_solo_nml')
  else
      Time_init  = set_date(date_init(1),date_init(2), date_init(3), &
           date_init(4),date_init(5),date_init(6))
  endif

  if (file_exist('INPUT/ocean_solo.res')) then
      Time_start =  set_date(date(1),date(2),date(3),date(4),date(5),date(6))
  else
      Time_start = Time_init
      date = date_init
  endif

  Time_end          = increment_date(Time_start, years, months, days, hours, minutes, seconds)
  Run_len           = Time_end - Time_start


  Time_step_coupled = set_time(dt_cpld, 0)
  num_cpld_calls    = Run_len / Time_step_coupled
  Time = Time_start

  Time_restart_init = set_date(date_restart(1), date_restart(2), date_restart(3),  &
                               date_restart(4), date_restart(5), date_restart(6) )
  Time_restart_current = Time_start
  if(ALL(restart_interval ==0)) then
     Time_restart = increment_date(Time_end, 1, 0, 0, 0, 0, 0)   ! no intermediate restart
  else
     Time_restart = increment_date(Time_restart_init, restart_interval(1), restart_interval(2), &
        restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
     if(Time_restart <= Time_start) call mpp_error(FATAL, &
       '==>Error from program ocean_solo: The first intermediate restart time is no larger than the start time')
  end if

  !-----------------------------------------------------------------------
  !------------------- some error checks ---------------------------------
  if ( num_cpld_calls * Time_step_coupled /= Run_len ) call mpp_error(FATAL, &
     '==>Error from program ocean_solo: run length must be multiple of cpld time step', FATAL)

  call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_APPEND,threading=MPP_SINGLE)

  month = month_name(date(2))
  if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

  call get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
  month = month_name(date(2))
  if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

  call mpp_close (unit)  

  !----- check the value of layout and setup the maskmap for domain layout.
  if( n_mask > 0 ) then
      if( layout_mask(1)*layout_mask(2) - n_mask .NE. mpp_npes() ) call mpp_error(FATAL, &
      '==>Error from program ocean_solo: layout(1)*layout(2) - n_mask should equal to npes when n_mask>0')
      call mpp_error(NOTE, &
      '==>Error from program ocean_solo: layout_mask and mask_list is set in ocean_solo_nml, ' // &
      'the value of layout_mask will override the layout specified in ocean_model_mod')

      allocate(Ocean_sfc%maskmap(layout_mask(1), layout_mask(2)) )
      Ocean_sfc%maskmap = .TRUE.
      do n=1, n_mask
         if (mask_list(1,n) .gt. layout_mask(1) ) &
              call mpp_error( FATAL, &
              'program ocean_solo: mask_list elements outside layout defines.' )
         if (mask_list(2,n) .gt. layout_mask(2) ) &
              call mpp_error( FATAL, &
              'program ocean_solo: mask_list elements outside layout defines.' )
         Ocean_sfc%maskmap(mask_list(1,n),mask_list(2,n)) = .false.
      enddo
  else
      if( layout_mask(1)*layout_mask(2) .NE. 0 ) call mpp_error(NOTE, &
           'program ocean_solo: when no region is masked out, layout_mask need not be set' )
  end if



  call ocean_model_init(Ocean_sfc, Ocean_state, Time_init, Time)

  call data_override_init(Ocean_domain_in = Ocean_sfc%domain)

  override_clock = mpp_clock_id('Override', flags=flags,grain=CLOCK_COMPONENT)
  
  call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)
  
  allocate ( Ice_ocean_boundary% u_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% v_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% t_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% q_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% salt_flux (isc:iec,jsc:jec),       &
             Ice_ocean_boundary% lw_flux (isc:iec,jsc:jec),         &
             Ice_ocean_boundary% sw_flux_vis_dir (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux_vis_dif (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux_nir_dir (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux_nir_dif (isc:iec,jsc:jec), &
             Ice_ocean_boundary% lprec (isc:iec,jsc:jec),           &
             Ice_ocean_boundary% fprec (isc:iec,jsc:jec),           &
             Ice_ocean_boundary% runoff (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% calving (isc:iec,jsc:jec),         &
             Ice_ocean_boundary% p (isc:iec,jsc:jec))

  Ice_ocean_boundary%u_flux          = 0.0
  Ice_ocean_boundary%v_flux          = 0.0
  Ice_ocean_boundary%t_flux          = 0.0
  Ice_ocean_boundary%q_flux          = 0.0
  Ice_ocean_boundary%salt_flux       = 0.0
  Ice_ocean_boundary%lw_flux         = 0.0
  Ice_ocean_boundary%sw_flux_vis_dir = 0.0
  Ice_ocean_boundary%sw_flux_vis_dif = 0.0
  Ice_ocean_boundary%sw_flux_nir_dir = 0.0
  Ice_ocean_boundary%sw_flux_nir_dif = 0.0
  Ice_ocean_boundary%lprec           = 0.0
  Ice_ocean_boundary%fprec           = 0.0
  Ice_ocean_boundary%runoff          = 0.0
  Ice_ocean_boundary%calving         = 0.0
  Ice_ocean_boundary%p               = 0.0

  call external_coupler_sbc_init(Ocean_sfc%domain, dt_cpld, Run_len)

  ! loop over the coupled calls 
  do nc=1, num_cpld_calls
     
     call mpp_clock_begin(override_clock)

     call ice_ocn_bnd_from_data(Ice_ocean_boundary)

     call mpp_clock_end(override_clock)

     call external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

     
     call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, Time, Time_step_coupled)

     Time = Time + Time_step_coupled

     if( Time >= Time_restart ) then
       Time_restart_current = Time
       Time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
             restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
       timestamp = date_to_string(time_restart_current)
        write(stdoutunit,*) '=> NOTE from program ocean_solo: intermediate restart file is written and ', &
             trim(timestamp),' is appended as prefix to each restart file name'
        call ocean_model_restart(Ocean_state, timestamp)
        call ocean_solo_restart(Time, Time_restart_current, timestamp)
     end if


     call external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

  enddo

  ! close some of the main components 

  call ocean_model_end(Ocean_sfc, Ocean_state, Time)

  call diag_manager_end(Time)

  ! need to reset pelist before calling mpp_clock_end
  call mpp_set_current_pelist() 

  ! write restart file
  call ocean_solo_restart(Time_end, Time_restart_current)

  call external_coupler_restart( dt_cpld, num_cpld_calls )

  call fms_io_exit

  call external_coupler_exit

  call fms_end

  call external_coupler_mpi_exit(mpi_comm_mom, external_initialization)

  write(stdoutunit,*) 'MOM: --- completed ---'

  contains

  !--- writing restart file that contains running time and restart file writing time.
  subroutine ocean_solo_restart(Time_run, Time_res, time_stamp)
    type(time_type),   intent(in)           :: Time_run, Time_res
    character(len=*),  intent(in), optional :: time_stamp
    character(len=128)                      :: file_run, file_res
    integer :: yr, mon, day, hr, min, sec

    ! write restart file
    if(present(time_stamp)) then
       file_run = 'RESTART/'//trim(time_stamp)//'.ocean_solo.res'
       file_res = 'RESTART/'//trim(time_stamp)//'.ocean_solo.intermediate.res'
    else
       file_run = 'RESTART/ocean_solo.res'
       file_res = 'RESTART/ocean_solo.intermediate.res'
    endif

    call mpp_open( unit, file_run, nohdrs=.TRUE. )
    if ( mpp_pe().EQ.mpp_root_pe() )then
       write( unit, '(i6,8x,a)' )calendar_type, &
            '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

       call get_date(Time_init,yr,mon,day,hr,min,sec)
       write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
            'Model start time:   year, month, day, hour, minute, second'
       call get_date(Time_run ,yr,mon,day,hr,min,sec)
       write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
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

end subroutine ocean_solo_restart

!====================================================================
! get forcing data from data_overide 
  subroutine ice_ocn_bnd_from_data(x)

      type (ice_ocean_boundary_type) :: x
      type(time_type)                :: Time_next

      Time_next = Time + Time_step_coupled
      call data_override('OCN', 't_flux',          x%t_flux         , Time_next)
      call data_override('OCN', 'u_flux',          x%u_flux         , Time_next)
      call data_override('OCN', 'v_flux',          x%v_flux         , Time_next)
      call data_override('OCN', 'q_flux',          x%q_flux         , Time_next)
      call data_override('OCN', 'salt_flux',       x%salt_flux      , Time_next)
      call data_override('OCN', 'lw_flux',         x%lw_flux        , Time_next)
      call data_override('OCN', 'sw_flux_vis_dir', x%sw_flux_vis_dir, Time_next)
      call data_override('OCN', 'sw_flux_vis_dif', x%sw_flux_vis_dif, Time_next)
      call data_override('OCN', 'sw_flux_nir_dir', x%sw_flux_nir_dir, Time_next)
      call data_override('OCN', 'sw_flux_nir_dif', x%sw_flux_nir_dif, Time_next)
      call data_override('OCN', 'lprec',           x%lprec          , Time_next)
      call data_override('OCN', 'fprec',           x%fprec          , Time_next)
      call data_override('OCN', 'runoff',          x%runoff         , Time_next)
      call data_override('OCN', 'calving',         x%calving        , Time_next)
      call data_override('OCN', 'p',               x%p              , Time_next)
            
  end subroutine ice_ocn_bnd_from_data


!-----------------------------------------------------------------------------------------
! 
! Subroutines  for enabling coupling to external programs through a third party coupler
! such as OASIS/PRISM.
! If no external coupler then these will mostly be dummy routines.
! These routines can also serve as spots to call other user defined routines
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
#ifdef OASIS3

! Here we provide some hooks for calling an interface between the OASIS3 coupler and MOM.
! The mom_oasis3_interface module is NOT general and it is expected that the user will 
! heavily modify it depending on the coupling strategy.
! For clarity all variables should be passed as arguments rather than as globals.
! This may require changes to the argument lists.

  subroutine external_coupler_mpi_init(mom_local_communicator, external_initialization)
! OASIS3/PRISM acts as the master and initializes MPI. Get a local communicator.
! need to initialize prism and get local communicator MPI_COMM_MOM first! 
  use mom_oasis3_interface_mod, only : mom_prism_init
  implicit none
  integer, intent(out) :: mom_local_communicator
  logical, intent(out) :: external_initialization
  mom_local_communicator = -100         ! Is there mpp_undefined parameter corresponding to MPI_UNDEFINED?
                                        ! probably wouldn't need logical flag.
  call mom_prism_init(mom_local_communicator)
  external_initialization = .true.
  end subroutine external_coupler_mpi_init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_sbc_init(Dom, dt_cpld, Run_len)
! Call to routine initializing arrays etc for transferring via coupler
! Perform sanity checks and make sure all inputs are compatible
  use mom_oasis3_interface_mod, only : coupler_init
  implicit none
  type(domain2d) :: Dom
  integer :: dt_cpld
  type(time_type) :: Run_len
  call coupler_init(Dom, dt_cpld, Run_len=Run_len)
  end  subroutine external_coupler_sbc_init
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
  subroutine external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
! Perform transfers before ocean time stepping
! May need special tratment on first call.

  use mom_oasis3_interface_mod, only : from_coupler, into_coupler 

  implicit none
  type (ice_ocean_boundary_type), intent(INOUT) :: Ice_ocean_boundary
  type (ocean_public_type) , intent(INOUT)        :: Ocean_sfc
  integer , intent(IN)                       :: nsteps, dt_cpld
  
  integer                        :: rtimestep ! Receive timestep
  integer                        :: stimestep ! Send timestep

  rtimestep = (nsteps-1) * dt_cpld   ! runtime in this run segment!
  stimestep = rtimestep
  call from_coupler( rtimestep, Ice_ocean_boundary )
  call into_coupler( stimestep, Ocean_sfc, before_ocean_update = .true.)
  end subroutine external_coupler_sbc_before
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
  subroutine external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
! Perform transfers after ocean time stepping

  use mom_oasis3_interface_mod, only : into_coupler

  implicit none
  type (ice_ocean_boundary_type) :: Ice_ocean_boundary
  type (ocean_public_type)         :: Ocean_sfc
  integer                        :: nsteps, dt_cpld
  
  integer                        :: stimestep ! Send timestep

  stimestep = nsteps * dt_cpld   ! runtime in this run segment!
  call into_coupler(stimestep, Ocean_sfc, before_ocean_update = .false.)
  end subroutine external_coupler_sbc_after
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  subroutine external_coupler_restart( dt_cpld, num_cpld_calls )
! Clean up as appropriate and write a restart
  use mom_oasis3_interface_mod, only : write_coupler_restart
  implicit none
  integer, intent(in)               :: dt_cpld, num_cpld_calls
  integer                           :: timestep
  timestep = num_cpld_calls * dt_cpld
  call write_coupler_restart(timestep, write_restart=.true.)
  end subroutine external_coupler_restart


!-----------------------------------------------------------------------------------------
  subroutine external_coupler_exit
! Clean up as appropriate
! Final call to external program
  use mom_oasis3_interface_mod, only : mom_prism_terminate
  call mom_prism_terminate
  end subroutine external_coupler_exit
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_mpi_exit(mom_local_communicator, external_initialization)
! mpp_exit wont call MPI_FINALIZE if mom_local_communicator /= MPI_COMM_WORLD
  implicit none
  integer, intent(in) :: mom_local_communicator
  logical, intent(in) :: external_initialization
  integer :: ierr
  call MPI_FINALIZE(ierr)
  return
  end subroutine external_coupler_mpi_exit
!-----------------------------------------------------------------------------------------

#else

! Dummy subroutines.

  subroutine external_coupler_mpi_init(mom_local_communicator, external_initialization)
  implicit none
  integer, intent(out) :: mom_local_communicator
  logical, intent(out) :: external_initialization
  external_initialization = .false.
  mom_local_communicator = -100         ! Is there mpp_undefined parameter corresponding to MPI_UNDEFINED?
                                        ! probably wouldn't need logical flag.
  return
  end subroutine external_coupler_mpi_init

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_sbc_init(Dom, dt_cpld, Run_len)
  implicit none
  type(domain2d) :: Dom
  integer :: dt_cpld
  type(time_type) :: Run_len
  return
  end  subroutine external_coupler_sbc_init

  subroutine external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
  implicit none
  type (ice_ocean_boundary_type), intent(INOUT) :: Ice_ocean_boundary
  type (ocean_public_type) , intent(INOUT)        :: Ocean_sfc
  integer , intent(IN)                       :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_before


  subroutine external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
  type (ice_ocean_boundary_type) :: Ice_ocean_boundary
  type (ocean_public_type)         :: Ocean_sfc
  integer                        :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_after

  subroutine external_coupler_restart( dt_cpld, num_cpld_calls )
  implicit none
  integer, intent(in)               :: dt_cpld, num_cpld_calls
  return
  end subroutine external_coupler_restart

  subroutine external_coupler_exit
  return
  end subroutine external_coupler_exit

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_mpi_exit(mom_local_communicator, external_initialization)
  implicit none
  integer, intent(in) :: mom_local_communicator
  logical, intent(in) :: external_initialization
  return
  end subroutine external_coupler_mpi_exit
!-----------------------------------------------------------------------------------------
#endif

end program main

#else

#include "ocean_solo_nuopc.inc"

#endif
