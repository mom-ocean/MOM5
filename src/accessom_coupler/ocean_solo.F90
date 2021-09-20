!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
!  
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison
!</CONTACT>
!
!<CONTACT EMAIL="Dave.Bi@csiro.au"> Dave Bi (for OASIS3 hooks)
!</CONTACT>
!
!<REVIEWER EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh (for OASIS3 hooks)
!</REVIEWER>
!
!<REVIEWER EMAIL="V.Balaji@noaa.gov"> V. Balaji
!</REVIEWER>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen Griffies 
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
! This driver also provides the prototype hooks for using MOM4p1 with OASIS3,
! with this code surrounded by the cpp-preprocessor option "ifdef OASIS3".  
! The couping of MOM4p1 to OASIS3 has not been tested at GFDL. Rather, 
! CSIRO in Australia uses MOM4p1 with OASIS3, with Dave.Bi@csiro.au the primary 
! contact for questions regarding MOM4p1 and OASIS3. 
! </DESCRIPTION>
!
! <NAMELIST NAME="ocean_solo_nml">
!
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

  use constants_mod,            only: constants_init, SECONDS_PER_HOUR
  use data_override_mod,        only: data_override_init, data_override
  use diag_manager_mod,         only: diag_manager_init, diag_manager_end
  use field_manager_mod,        only: field_manager_init
  use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
  use fms_mod,                  only: close_file, file_exist, uppercase
  use fms_io_mod,               only: fms_io_exit
  use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain
  use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
  use mpp_mod,                  only: mpp_init
  use mpp_mod,                  only: mpp_error, FATAL, NOTE, mpp_pe, mpp_npes, mpp_set_current_pelist, mpp_sync
  use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, MAXPES
  use mpp_mod,                  only: mpp_broadcast
  use time_interp_external_mod, only: time_interp_external_init
  use time_manager_mod,         only: set_calendar_type, time_type, increment_date
  use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name, print_time
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
  use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )
  use time_manager_mod,         only: operator( * ), operator( /= ), operator( > )
  use time_manager_mod,         only: date_to_string

  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end
  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_types_mod,          only: ice_ocean_boundary_type
  use ocean_util_mod,           only: write_chksum_2d

  use auscom_ice_parameters_mod, only: redsea_gulfbay_sfix, do_sfix_now, sfix_hours, int_sec
  use accessom2_mod, only : accessom2_type => accessom2
  use coupler_mod, only : coupler_type => coupler

  implicit none

  type (ocean_public_type)               :: Ocean_sfc
  type (ocean_state_type),       pointer :: Ocean_state
  type(ice_ocean_boundary_type), target  :: Ice_ocean_boundary
  type(accessom2_type) :: accessom2
  type(coupler_type) :: coupler

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
  type(time_type) :: Time_last_sfix 
  type(time_type) :: Time_sfix 
  integer :: sfix_seconds

  character(len=17) :: calendar = 'julian'

  integer :: num_cpld_calls  = 0
  integer :: nc
  integer :: calendar_type=-1

  integer :: date_init(6)=0, date(6)
  integer :: date_restart(6)
  integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_cpld  = 86400
  integer :: yy, mm, dd, hh, mimi, ss

  integer :: isc,iec,jsc,jec
  integer :: unit, io_status, ierr

  integer :: flags=0
  integer :: init_clock, main_clock, term_clock
  integer :: override_clock, coupler_init_clock
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
  integer ::  stdoutunit, stdlogunit, tmp_unit
  logical :: debug_this_module
  character(len=1024) :: accessom2_config_dir = '../'
  integer, dimension(6) :: date_array

  namelist /ocean_solo_nml/ n_mask, layout_mask, mask_list, restart_interval, &
                            debug_this_module, accessom2_config_dir

    ! Initialise floating point exception error handler.
    !call fpe_err_handler_init()
  debug_this_module = .false.

  open(newunit=tmp_unit, file='input.nml')
  read(tmp_unit, nml=ocean_solo_nml)
  close(tmp_unit)

  call coupler%init_begin('mom5xx', config_dir=trim(accessom2_config_dir))

  call mpp_init(localcomm=coupler%localcomm)
  init_clock = mpp_clock_id('Initialization')
  main_clock = mpp_clock_id('Main Loop')
  term_clock = mpp_clock_id('Termination')
  call mpp_clock_begin(init_clock)

  call fms_init(coupler%localcomm)

  call constants_init()
  flags = MPP_CLOCK_SYNC

  stdoutunit=stdout();stdlogunit=stdlog()

  ! provide for namelist over-ride
  unit = open_namelist_file('input.nml')
  read  (unit, ocean_solo_nml,iostat=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,'(/47x,a/)') '======== MODEL BEING DRIVEN BY OCEAN_SOLO_MOD ========'
  write (stdoutunit, ocean_solo_nml)  
  write (stdlogunit, ocean_solo_nml)
  ierr = check_nml_error(io_status,'ocean_solo_nml')
  call close_file (unit)

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

  ! Initialise libaccessom2
  call accessom2%init('mom5xx', config_dir=trim(accessom2_config_dir))

  if (mpp_pe() == mpp_root_pe()) then
    call accessom2%print_version_info()
  endif

  ! Tell libaccessom2 about any global configs/state

  ! Synchronise accessom2 'state' (i.e. configuration) between all models.
  call accessom2%sync_config(coupler)

  ! Use accessom2 configuration to set calendar
  if (index(accessom2%get_calendar_type(), 'noleap') > 0) then
      calendar_type = NOLEAP
  elseif (index(accessom2%get_calendar_type(), 'gregorian') > 0) then
      calendar_type = GREGORIAN
  else
    call mpp_error(FATAL, 'ocean_solo: unsupported calendar type')
  endif

  ! Use accessom2 configuration to initial date and runtime
  date_init(:) = accessom2%get_cur_exp_date_array()
  dt_cpld = accessom2%get_ice_ocean_timestep()
  years = 0
  months = 0
  days = 0
  hours = 0
  minutes = 0
  seconds = accessom2%get_total_runtime_in_seconds()

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
      date_restart = date
  endif

  call set_calendar_type (calendar_type)

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
  if(ALL(restart_interval == 0)) then
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

  call ocean_model_init(Ocean_sfc, Ocean_state, Time_init, Time, &
                        accessom2%get_ice_ocean_timestep())

  if (redsea_gulfbay_sfix) then
    ! This must be called after ocean_model_init so sfix_hours is read in from namelist
    sfix_seconds = sfix_hours * SECONDS_PER_HOUR
    ! Get current model time from Time_init in seconds (must be done like this otherwise
    ! can get an overflow in seconds)
    call get_time(Time-Time_init,seconds=seconds,days=days)
    ! The last sfix time has to be determined from absolute model time, to ensure reproducibility 
    ! across restarts

    ! Current time to nearest hour
    hours = days*24 + int(seconds/SECONDS_PER_HOUR)

    ! Time of last sfix 
    hours = int(hours / sfix_hours) * sfix_hours

    ! Convert to days + hours
    days = int(hours / 24)
    hours = hours - days*24

    Time_last_sfix = set_time(days=int(days),seconds=int(hours*SECONDS_PER_HOUR)) + Time_init
    Time_sfix = set_time(seconds=int(sfix_seconds))

    call print_time(Time_last_sfix,'Time_last_sfix: ')
  end if

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
             Ice_ocean_boundary% p (isc:iec,jsc:jec),               &
             Ice_ocean_boundary% aice(isc:iec,jsc:jec),             &
             Ice_ocean_boundary% mh_flux(isc:iec,jsc:jec),          &
             Ice_ocean_boundary% wfimelt(isc:iec,jsc:jec),          &
             Ice_ocean_boundary% wfiform(isc:iec,jsc:jec),          &
             Ice_ocean_boundary% licefw(isc:iec,jsc:jec),           &
             Ice_ocean_boundary% liceht(isc:iec,jsc:jec),           &
             Ice_ocean_boundary%wnd(isc:iec,jsc:jec))
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  allocate ( Ice_ocean_boundary%iof_nit(isc:iec,jsc:jec),           &
             Ice_ocean_boundary%iof_alg(isc:iec,jsc:jec))
#endif
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
  Ice_ocean_boundary%aice            = 0.0
  Ice_ocean_boundary%mh_flux         = 0.0
  Ice_ocean_boundary% wfimelt        = 0.0
  Ice_ocean_boundary% wfiform        = 0.0
  Ice_ocean_boundary%licefw          = 0.0
  Ice_ocean_boundary%liceht          = 0.0
  Ice_ocean_boundary%wnd             = 0.0
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  Ice_ocean_boundary%iof_nit         = 0.0
  Ice_ocean_boundary%iof_alg         = 0.0
#endif
  coupler_init_clock = mpp_clock_id('OASIS init', grain=CLOCK_COMPONENT)
  call mpp_clock_begin(coupler_init_clock)
  call external_coupler_sbc_init(Ocean_sfc%domain, dt_cpld, Run_len, &
                                 accessom2%get_coupling_field_timesteps())
  call mpp_clock_end(coupler_init_clock)
  call mpp_clock_end(init_clock)

  ! loop over the coupled calls
  call mpp_clock_begin(main_clock)
  do nc=1, num_cpld_calls

     call external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld)

     call mpp_clock_begin(override_clock)
     call ice_ocn_bnd_from_data(Ice_ocean_boundary)
     call mpp_clock_end(override_clock)

     if (debug_this_module) then
        call write_boundary_chksums(Ice_ocean_boundary)
     endif

    if (redsea_gulfbay_sfix) then
        if ((Time - Time_last_sfix) >= Time_sfix) then
            do_sfix_now = .true.
            Time_last_sfix = Time
        else
            do_sfix_now = .false.
        end if
    end if

     call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, Time, Time_step_coupled)

     Time = Time + Time_step_coupled
     if ( mpp_pe() == mpp_root_pe() ) then
        call accessom2%progress_date(int(dt_cpld))
     endif

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
  call mpp_clock_end(main_clock)

  call mpp_clock_begin(term_clock)

  call external_coupler_restart( dt_cpld, num_cpld_calls, Ocean_sfc)

  ! close some of the main components 
  call ocean_model_end(Ocean_sfc, Ocean_state, Time)

  call diag_manager_end(Time)

  ! need to reset pelist before calling mpp_clock_end
  ! call mpp_set_current_pelist()

  ! write restart file
  call ocean_solo_restart(Time_end, Time_restart_current)

  call fms_io_exit

  call coupler%deinit()
  ! Allow libaccessom2 to check that all models are synchronised at the end of
  ! the run.
  call get_date(Time, date_array(1), date_array(2), date_array(3), &
                date_array(4), date_array(5), date_array(6))
  call accessom2%deinit(cur_date_array=date_array)

  call mpp_clock_end(term_clock)

  call fms_end

  call external_coupler_mpi_exit(coupler%localcomm, .true.)

  print *, 'MOM5: --- completed ---'

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
      call data_override('OCN', 'aice',            x%aice           , Time_next)
      call data_override('OCN', 'mh_flux',         x%mh_flux        , Time_next)
    call mpp_sync()
            
end subroutine ice_ocn_bnd_from_data



! Here we provide some hooks for calling an interface between the OASIS3 coupler and MOM.
! The mom_oasis3_interface module is NOT general and it is expected that the user will 
! heavily modify it depending on the coupling strategy.
! For clarity all variables should be passed as arguments rather than as globals.
! This may require changes to the argument lists.

! NOTE: libaccessom2 makes most of these functions redundant.

subroutine external_coupler_sbc_init(Dom, dt_cpld, Run_len, &
                                     coupling_field_timesteps)
    ! Call to routine initializing arrays etc for transferring via coupler
    ! Perform sanity checks and make sure all inputs are compatible
    use mom_oasis3_interface_mod, only : coupler_init
    implicit none
    type(domain2d) :: Dom
    integer :: dt_cpld
    type(time_type) :: Run_len
    integer, dimension(:), intent(in) :: coupling_field_timesteps

    call coupler_init(Dom, dt_cpld=dt_cpld, Run_len=Run_len, &
                      coupling_field_timesteps=coupling_field_timesteps)
end  subroutine external_coupler_sbc_init

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
    call from_coupler( rtimestep, Ocean_sfc, Ice_ocean_boundary )
    call into_coupler( stimestep, Ocean_sfc, before_ocean_update = .true.)
end subroutine external_coupler_sbc_before

subroutine external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
    !Perform transfers after ocean time stepping

    use mom_oasis3_interface_mod, only : into_coupler

    implicit none
    type (ice_ocean_boundary_type) :: Ice_ocean_boundary
    type (ocean_public_type)         :: Ocean_sfc
    integer                        :: nsteps, dt_cpld

    integer                        :: stimestep ! Send timestep

    stimestep = nsteps * dt_cpld   ! runtime in this run segment!
    if (stimestep < num_cpld_calls*dt_cpld) call into_coupler(stimestep, Ocean_sfc, before_ocean_update = .false.)
end subroutine external_coupler_sbc_after

subroutine external_coupler_restart( dt_cpld, num_cpld_calls, Ocean_sfc)
    !Clean up as appropriate and write a restart
    use mom_oasis3_interface_mod, only : write_coupler_restart
    implicit none
    integer, intent(in)               :: dt_cpld, num_cpld_calls
    integer                           :: timestep
    type (ocean_public_type)         :: Ocean_sfc

    timestep = num_cpld_calls * dt_cpld
    call write_coupler_restart(timestep, Ocean_sfc, write_restart=.true.)
end subroutine external_coupler_restart

subroutine external_coupler_mpi_exit(mom_local_communicator, external_initialization)
    ! mpp_exit wont call MPI_FINALIZE if mom_local_communicator /= MPI_COMM_WORLD
    implicit none
    integer, intent(in) :: mom_local_communicator
    logical, intent(in) :: external_initialization
    integer :: ierr
    call MPI_FINALIZE(ierr)
    return
end subroutine external_coupler_mpi_exit


subroutine write_boundary_chksums(Ice_ocean_boundary)
    type(ice_ocean_boundary_type), intent(in) :: Ice_ocean_boundary

    call write_chksum_2d('Ice_ocean_boundary%u_flux', Ice_ocean_boundary%u_flux)
    call write_chksum_2d('Ice_ocean_boundary%v_flux', Ice_ocean_boundary%v_flux)
    call write_chksum_2d('Ice_ocean_boundary%t_flux', Ice_ocean_boundary%t_flux)
    call write_chksum_2d('Ice_ocean_boundary%q_flux', Ice_ocean_boundary%q_flux)
    call write_chksum_2d('Ice_ocean_boundary%salt_flux', Ice_ocean_boundary%salt_flux)
    call write_chksum_2d('Ice_ocean_boundary%lw_flux', Ice_ocean_boundary%lw_flux)
    call write_chksum_2d('Ice_ocean_boundary%sw_flux_vis_dir', Ice_ocean_boundary%sw_flux_vis_dir)
    call write_chksum_2d('Ice_ocean_boundary%sw_flux_vis_dif', Ice_ocean_boundary%sw_flux_vis_dif)
    call write_chksum_2d('Ice_ocean_boundary%sw_flux_nir_dir', Ice_ocean_boundary%sw_flux_nir_dir)
    call write_chksum_2d('Ice_ocean_boundary%sw_flux_nir_dif', Ice_ocean_boundary%sw_flux_nir_dif)
    call write_chksum_2d('Ice_ocean_boundary%lprec', Ice_ocean_boundary%lprec)
    call write_chksum_2d('Ice_ocean_boundary%fprec', Ice_ocean_boundary%fprec)
    call write_chksum_2d('Ice_ocean_boundary%runoff', Ice_ocean_boundary%runoff)
    call write_chksum_2d('Ice_ocean_boundary%calving', Ice_ocean_boundary%calving)
    call write_chksum_2d('Ice_ocean_boundary%p', Ice_ocean_boundary%p)
    call write_chksum_2d('Ice_ocean_boundary%aice', Ice_ocean_boundary%aice)
    call write_chksum_2d('Ice_ocean_boundary%mh_flux', Ice_ocean_boundary%mh_flux)
    call write_chksum_2d('Ice_ocean_boundary%wfimelt', Ice_ocean_boundary%wfimelt)
    call write_chksum_2d('Ice_ocean_boundary%wfiform', Ice_ocean_boundary%wfiform)
    call write_chksum_2d('Ice_ocean_boundary%wnd', Ice_ocean_boundary%wnd)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
    call write_chksum_2d('Ice_ocean_boundary%iof_nit', Ice_ocean_boundary%iof_nit)
    call write_chksum_2d('Ice_ocean_boundary%iof_alg', Ice_ocean_boundary%iof_alg)
#endif
end subroutine write_boundary_chksums

end program main
