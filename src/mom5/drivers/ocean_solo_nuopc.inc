!---------------------------------------------------------------
!------------------ Utilities ----------------------------------
!
!This is a set of utilities that convert time management data structures
!between ESMF and FMS.
!
!Author: Fei Liu
!Wed Feb 20 09:15:23 MST 2013
!---------------------------------------------------------------

module time_utils_mod

  use fms_mod,                  only: uppercase
  use mpp_mod,                  only: mpp_error, FATAL
  use time_manager_mod,         only: time_type, set_time, set_date, get_date
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
  use ESMF
  
  implicit none
  private

  !-------------------- interface blocks ---------------------
  interface fms2esmf_cal
    module procedure fms2esmf_cal_c
    module procedure fms2esmf_cal_i
  end interface fms2esmf_cal
  interface esmf2fms_time
    module procedure esmf2fms_time_t
    module procedure esmf2fms_timestep
  end interface esmf2fms_time

  public fms2esmf_cal
  public esmf2fms_time
  public fms2esmf_time
  public string_to_date

  contains

  !-------------------- module code ---------------------

  function fms2esmf_cal_c(calendar)
!    ! Return Value:
    type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_c
!    ! Arguments:
    character(len=*), intent(in)       :: calendar

    select case( uppercase(trim(calendar)) )
    case( 'GREGORIAN' )
       fms2esmf_cal_c = ESMF_CALKIND_GREGORIAN
    case( 'JULIAN' )
       fms2esmf_cal_c = ESMF_CALKIND_JULIAN
    case( 'NOLEAP' )
       fms2esmf_cal_c = ESMF_CALKIND_NOLEAP
    case( 'THIRTY_DAY' )
       fms2esmf_cal_c = ESMF_CALKIND_360DAY
    case( 'NO_CALENDAR' )
       fms2esmf_cal_c = ESMF_CALKIND_NOCALENDAR
    case default
       call mpp_error(FATAL, &
       'ocean_solo: ocean_solo_nml entry calendar must be one of GREGORIAN|JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
    end select 
  end function fms2esmf_cal_c

  function fms2esmf_cal_i(calendar)
!    ! Return Value:
    type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_i
!    ! Arguments:
    integer, intent(in)                :: calendar

    select case(calendar)
      case(THIRTY_DAY_MONTHS)
        fms2esmf_cal_i = ESMF_CALKIND_360DAY
      case(GREGORIAN)
        fms2esmf_cal_i = ESMF_CALKIND_GREGORIAN
      case(JULIAN)
        fms2esmf_cal_i = ESMF_CALKIND_JULIAN
      case(NOLEAP)
        fms2esmf_cal_i = ESMF_CALKIND_NOLEAP
      case(NO_CALENDAR)
        fms2esmf_cal_i = ESMF_CALKIND_NOCALENDAR
    end select
  end function fms2esmf_cal_i

  function esmf2fms_time_t(time)
    ! Return Value
    type(Time_type)                    :: esmf2fms_time_t
    ! Input Arguments
    type(ESMF_Time), intent(in)        :: time
    ! Local Variables
    integer                            :: yy, mm, dd, h, m, s
    type(ESMF_CALKIND_FLAG)            :: calkind

    integer                            :: rc

    call ESMF_TimeGet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
        calkindflag=calkind, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    esmf2fms_time_t = Set_date(yy, mm, dd, h, m, s)

  end function esmf2fms_time_t

  function esmf2fms_timestep(timestep)
    ! Return Value
    type(Time_type)                    :: esmf2fms_timestep
    ! Input Arguments
    type(ESMF_TimeInterval), intent(in):: timestep
    ! Local Variables
    integer                            :: s
    type(ESMF_CALKIND_FLAG)            :: calkind

    integer                            :: rc

    call ESMF_TimeIntervalGet(timestep, s=s, calkindflag=calkind, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    esmf2fms_timestep = set_time(s, 0)

  end function esmf2fms_timestep

  function fms2esmf_time(time, calkind)
    ! Return Value
    type(ESMF_Time)                    :: fms2esmf_time
    ! Input Arguments
    type(Time_type), intent(in)        :: time
    type(ESMF_CALKIND_FLAG), intent(in), optional :: calkind
    ! Local Variables
    integer                            :: yy, mm, d, h, m, s
    type(ESMF_CALKIND_FLAG)            :: l_calkind

    integer                            :: rc

    if(present(calkind)) then
      l_calkind = calkind
    else
      l_calkind = fms2esmf_cal(fms_get_calendar_type())
    endif

    call get_date(time, yy, mm, d, h, m, s)

    call ESMF_TimeSet(fms2esmf_time, yy=yy, mm=mm, d=d, h=h, m=m, s=s, &
        calkindflag=l_calkind, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end function fms2esmf_time

  function string_to_date(string, rc)
    character(len=15), intent(in)           :: string
    integer, intent(out), optional          :: rc
    type(time_type)                         :: string_to_date

    integer                                 :: yr,mon,day,hr,min,sec
    
    if(present(rc)) rc = ESMF_SUCCESS

    read(string, '(I4.4,I2.2,I2.2,".",I2.2,I2.2,I2.2)') yr, mon, day, hr, min, sec
    string_to_date = set_date(yr, mon, day, hr, min, sec)

  end function string_to_date

end module time_utils_mod

!--------------- MOM5 Ocean solo model -----------------
! This is the MOM5 ocean solo model component that's NUOPC compiant.
! The public ocean_solo_SS method sets up all the model services such as
! initialize, run and finalize.
!
! Author:  Fei.Liu@gmail.com
!

module ocean_solo_mod
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
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end
  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_types_mod,          only: ice_ocean_boundary_type

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS      => routine_SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance

  use time_utils_mod

  implicit none
  private
  public ocean_solo_SS

  type ocean_internalstate_type
    type(ocean_public_type),       pointer :: ocean_public_type_ptr
    type(ocean_state_type),        pointer :: ocean_state_type_ptr
    type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
  end type

  type ocean_internalstate_wrapper
    type(ocean_internalstate_type), pointer :: ptr
  end type

  contains
  !-----------------------------------------------------------------------
  !------------------- Solo Ocean code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine ocean_solo_SS(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call model_routine_SS(gcomp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
      userRoutine=ocean_solo_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    ! No need to change clock settings
    !call ESMF_MethodAdd(gcomp, label=model_label_SetClock, &
    !  userRoutine=SetClock, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    
    call ESMF_MethodAdd(gcomp, label=model_label_Advance, &
      userRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_VM)      :: vm
    
    type (ocean_public_type),      pointer :: Ocean_sfc          
    type (ocean_state_type),       pointer :: Ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary 

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


    integer :: num_cpld_calls  = 0
    integer :: nc
    integer :: calendar_type=-1


    integer :: isc,iec,jsc,jec
    integer :: unit, io_status, ierr

    integer :: flags=0, override_clock
    integer :: nfields 
    
    character(len=256) :: version = ''
    character(len=256) :: tag = ''

    character(len=9) :: month
    character(ESMF_MAXSTR):: timestamp

    integer :: n, m 

    integer :: date_init(6)=0, date(6)
    integer :: date_restart(6)
    character(len=17) :: calendar = 'julian'
    integer :: dt_cpld  = 86400
    integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
    integer :: layout_mask(2) = (/0 , 0/)
    integer :: n_mask = 0
    integer :: mask_list(2,MAXPES)
    integer, parameter :: mp = 2*MAXPES
    data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/
    integer :: restart_interval(6) = (/0,0,0,0,0,0/)

    namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds, dt_cpld, &
                              n_mask, layout_mask, mask_list, restart_interval

    integer :: mpi_comm_mom
    integer ::  stdoutunit,stdlogunit
    logical :: external_initialization

    type(ocean_internalstate_wrapper)      :: ocean_internalstate

    rc = ESMF_SUCCESS

    allocate(Ice_ocean_boundary)
    !allocate(Ocean_state) ! ocean_model_init allocate this pointer
    allocate(Ocean_sfc)
    allocate(ocean_internalstate%ptr)
    ocean_internalstate%ptr%ice_ocean_boundary_type_ptr => Ice_ocean_boundary
    ocean_internalstate%ptr%ocean_public_type_ptr       => Ocean_sfc

    ! To start we have a dummy ocean that will not import or export anything.
    !
    !! importable field: air_pressure_at_sea_level
    !call NUOPC_StateAdvertiseField(importState, &
    !  StandardName="air_pressure_at_sea_level", rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !
    !! importable field: isotropic_shortwave_radiance_in_air
    !call NUOPC_StateAdvertiseField(importState, &
    !  StandardName="isotropic_shortwave_radiance_in_air", rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    ! exportable field: sea_surface_temperature
    !call NUOPC_StateAdvertiseField(exportState, &
    !  StandardName="sea_surface_temperature", rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    !call external_coupler_mpi_init(mpi_comm_mom, external_initialization)

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

    !override_clock = mpp_clock_id('Override', flags=flags,grain=CLOCK_COMPONENT)
    
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

    ocean_internalstate%ptr%ocean_state_type_ptr => Ocean_state
    call ESMF_GridCompSetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Set some MOM_solo specific attributes to be used during finalize
    call ESMF_AttributeSet(gcomp, name="dt_cpld", value=dt_cpld, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeSet(gcomp, name="num_cpld_calls", value=num_cpld_calls, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeSet(gcomp, name="calendar_type", value=calendar_type, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    timestamp = date_to_string(time_restart)
    call ESMF_AttributeSet(gcomp, name="Time_restart", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    timestamp = date_to_string(time_restart_current)
    call ESMF_AttributeSet(gcomp, name="Time_restart_current", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    timestamp = date_to_string(time_end)
    call ESMF_AttributeSet(gcomp, name="Time_end", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    timestamp = date_to_string(time_init)
    call ESMF_AttributeSet(gcomp, name="Time_init", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    timestamp = date_to_string(time_start)
    call ESMF_AttributeSet(gcomp, name="Time_start", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_TimeInterval) :: stabilityTimeStep
    type(ESMF_Field)        :: field
    type(ESMF_Grid)         :: gridIn
    type(ESMF_Grid)         :: gridOut
    
    rc = ESMF_SUCCESS

    ! To start we have a dummy ocean that will not import or export anything.
    !
    !! create a Grid object for Fields
    !gridIn = NUOPC_GridCreateSimpleXY(10._ESMF_KIND_R8, 20._ESMF_KIND_R8, &
    !  100._ESMF_KIND_R8, 200._ESMF_KIND_R8, 10, 100, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !gridOut = gridIn ! for now out same as in
    !! importable field: air_pressure_at_sea_level
    !field = ESMF_FieldCreate(name="pmsl", grid=gridIn, &
    !  typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call NUOPC_StateRealizeField(importState, field=field, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !
    !! importable field: isotropic_shortwave_radiance_in_air
    !field = ESMF_FieldCreate(name="risw", grid=gridIn, &
    !  typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call NUOPC_StateRealizeField(importState, field=field, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    !! exportable field: sea_surface_temperature
    !field = ESMF_FieldCreate(name="sst", grid=gridOut, &
    !  typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call NUOPC_StateRealizeField(exportState, field=field, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  ! Ocean solo model uses same clock as parent gridComp
  !subroutine SetClock(gcomp, rc)
  !  type(ESMF_GridComp)  :: gcomp
  !  integer, intent(out) :: rc
  !  
  !  ! local variables
  !  type(ESMF_Clock)              :: clock
  !  type(ESMF_TimeInterval)       :: stabilityTimeStep

  !  rc = ESMF_SUCCESS
  !  
  !  ! query the Component for its clock, importState and exportState
  !  call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
  !  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  !    line=__LINE__, &
  !    file=__FILE__)) &
  !    return  ! bail out
  !    
  !  ! initialize internal clock
  !  ! here: parent Clock and stability timeStep determine actual model timeStep
  !  !TODO: stabilityTimeStep should be read in from configuation
  !  !TODO: or computed from internal Grid information
  !  call ESMF_TimeIntervalSet(stabilityTimeStep, m=60, rc=rc) 
  !  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  !    line=__LINE__, &
  !    file=__FILE__)) &
  !    return  ! bail out
  !  call NUOPC_GridCompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
  !  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  !    line=__LINE__, &
  !    file=__FILE__)) &
  !    return  ! bail out

  !  
  !end subroutine

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep

    type (ocean_public_type),      pointer :: Ocean_sfc          
    type (ocean_state_type),       pointer :: Ocean_state
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary 

    ! define some time types 
    type(time_type) :: Time        
    type(time_type) :: Time_step_coupled
    type(time_type) :: Time_restart
    type(time_type) :: Time_restart_current

    character(len=17) :: calendar = 'julian'

    integer :: dt_cpld  = 86400
    integer :: nc
    integer :: calendar_type=-1

    integer :: date_init(6)=0, date(6)
    integer :: date_restart(6)
    integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0

    integer :: isc,iec,jsc,jec

    integer :: flags=0, override_clock
    integer :: nfields 
    
    character(len=9) :: month
    character(len=64):: timestamp
    integer :: restart_interval(6) = (/1,0,0,0,0,0/)
    type(ocean_internalstate_wrapper)      :: ocean_internalstate


    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call NUOPC_ClockPrintCurrTime(clock, &
      "------>Advancing OCN from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call NUOPC_TimePrint(currTime + timeStep, &
      "--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Time = esmf2fms_time(currTime)
    Time_step_coupled = esmf2fms_time(timeStep)

    ! Retrieve some MOM_solo specific attribute
    call ESMF_AttributeGet(gcomp, name="dt_cpld", value=dt_cpld, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Time_restart", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time_restart = string_to_date(timestamp)
    call ESMF_AttributeGet(gcomp, name="Time_restart_current", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time_restart_current = string_to_date(timestamp)

    override_clock = mpp_clock_id('Override', flags=flags,grain=CLOCK_COMPONENT)
    call mpp_clock_begin(override_clock)

    call ice_ocn_bnd_from_data(Ice_ocean_boundary, Time, Time_step_coupled)

    call mpp_clock_end(override_clock)

    call external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

    call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, Time, Time_step_coupled)

    Time = Time + Time_step_coupled

    if( Time >= Time_restart ) then
      Time_restart_current = Time
      Time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
            restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
      timestamp = date_to_string(time_restart_current)
       write(*,*) '=> NOTE from program ocean_solo: intermediate restart file is written and ', &
            trim(timestamp),' is appended as prefix to each restart file name'
       call ocean_model_restart(Ocean_state, timestamp)
       call ocean_solo_restart(gcomp, Time, Time_restart_current, timestamp)
    end if

    call external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

    ! Set Time_restart_current
    timestamp = date_to_string(time_restart_current)
    call ESMF_AttributeSet(gcomp, name="Time_restart_current", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine 

  subroutine ocean_solo_finalize(gcomp, importState, exportState, clock, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    type (ocean_public_type),      pointer :: Ocean_sfc          
    type (ocean_state_type),       pointer :: Ocean_state
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    type(time_type)                        :: Time, Time_restart_current, Time_end
    type(ESMF_Time)                        :: currTime
    type(ESMF_VM)                          :: vm
    integer                                :: mpi_comm_mom, dt_cpld, num_cpld_calls
    character(len=64)                      :: timestamp

    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time = esmf2fms_time(currTime)

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm_mom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Retrieve some MOM_solo specific attributes
    call ESMF_AttributeGet(gcomp, name="dt_cpld", value=dt_cpld, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(gcomp, name="num_cpld_calls", value=num_cpld_calls, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Time_end", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time_end = string_to_date(timestamp)
    call ESMF_AttributeGet(gcomp, name="Time_restart_current", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time_restart_current= string_to_date(timestamp)

    ! close some of the main components 

    call ocean_model_end(Ocean_sfc, Ocean_state, Time)

    call diag_manager_end(Time)

    !! need to reset pelist before calling mpp_clock_end
    call mpp_set_current_pelist() 

    !! write restart file
    call ocean_solo_restart(gcomp, Time_end, Time_restart_current)

    call external_coupler_restart( dt_cpld, num_cpld_calls )

    call fms_io_exit

    call external_coupler_exit

    call fms_end

    call external_coupler_mpi_exit(mpi_comm_mom, .true.)

    write(*,*) 'MOM: --- completed ---'

  end subroutine


  !--- writing restart file that contains running time and restart file writing time.
  subroutine ocean_solo_restart(gcomp, Time_run, Time_res, time_stamp)
    type(ESMF_GridComp), intent(in)         :: gcomp
    type(time_type),   intent(in)           :: Time_run, Time_res
    character(len=*),  intent(in), optional :: time_stamp

    integer                                 :: unit, calendar_type, rc
    character(len=64)                       :: timestamp
    type(Time_type)                         :: Time_init, Time_start
    character(len=128)                      :: file_run, file_res
    integer                                 :: yr, mon, day, hr, min, sec

    call ESMF_AttributeGet(gcomp, name="calendar_type", value=calendar_type, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Time_init", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time_init = string_to_date(timestamp)
    call ESMF_AttributeGet(gcomp, name="Time_start", value=timestamp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time_start = string_to_date(timestamp)

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
  subroutine ice_ocn_bnd_from_data(x, Time, Time_step_coupled)

      type (ice_ocean_boundary_type) :: x
      type(Time_type), intent(in)    :: Time, Time_step_coupled

      type(Time_type)                :: Time_next

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

!-----------------------------------------------------------------------------------------

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

end module ocean_solo_mod


!----------------Earth System Model driver -------------------
! Standard 1 component earth system model driver NUOPC compliant
!
module esm_mod
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
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end
  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_types_mod,          only: ice_ocean_boundary_type


  use ESMF
  use NUOPC
  use NUOPC_Driver, only: &
    driver_routine_SS             => routine_SetServices, &
    driver_type_IS                => type_InternalState, &
    driver_label_IS               => label_InternalState, &
    driver_setMC                  => label_SetModelCount, &
#ifdef WITHPETLISTS
    driver_label_SetModelPetLists => label_SetModelPetLists, &
#endif
    driver_label_SetModelServices => label_SetModelServices

  use time_utils_mod

  use ocean_solo_mod

  implicit none
  private
  public SetServices

  contains

  !-----------------------------------------------------------------------
  !------------------- Earth System Model code starts here ---------------
  !-----------------------------------------------------------------------
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! NUOPC_Driver registers the generic methods
    call driver_routine_SS(gcomp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! attach specializing method(s)
#ifdef WITHPETLISTS
    call ESMF_MethodAdd(gcomp, label=driver_label_SetModelPetLists, &
      userRoutine=SetModelPetLists, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    call ESMF_MethodAdd(gcomp, label=driver_setMC, &
        userRoutine=setModelCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_MethodAdd(gcomp, label=driver_label_SetModelServices, &
      userRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

#ifdef WITHPETLISTS

  subroutine SetModelPetLists(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                       :: localrc
    type(driver_type_IS)          :: is
    integer                       :: petCount, i

    rc = ESMF_SUCCESS

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, driver_label_IS, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! get the petCount
    call ESMF_GridCompGet(gcomp, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! set petList for OCN -> second half of PETs
    allocate(is%wrap%modelPetLists(petCount))
    do i=0, petCount-1
      is%wrap%modelPetLists(i) = i ! PET labeling goes from 0 to petCount-1
    enddo

  end subroutine
#endif

  subroutine setModelCount(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                       :: localrc
    type(driver_type_IS)          :: is
    integer                       :: petCount, i

    rc = ESMF_SUCCESS

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, driver_label_IS, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    !Only ocean model exists
    is%wrap%modelCount=1 

  end subroutine

  subroutine SetModelServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                       :: localrc
    type(driver_type_IS)          :: is
    type(ESMF_Grid)               :: grid
    type(ESMF_Field)              :: field
    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeInterval
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Clock)              :: internalClock

    integer :: n, m, unit, io_status, ierr
    ! Read clock information from namelist
    integer :: date_init(6)=0, date(6)
    integer :: date_restart(6)
    character(len=17) :: calendar = 'julian'
    integer :: dt_cpld  = 86400
    integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
    integer :: layout_mask(2) = (/0 , 0/)
    integer :: n_mask = 0
    integer :: mask_list(2,MAXPES)
    integer, parameter :: mp = 2*MAXPES
    data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/
    integer :: restart_interval(6) = (/0,0,0,0,0,0/)

    namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds, dt_cpld, &
                              n_mask, layout_mask, mask_list, restart_interval

    rc = ESMF_SUCCESS
    
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_solo_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_solo_nml')
#else
    unit = open_namelist_file('input.nml')
    read  (unit, ocean_solo_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_solo_nml')
    call close_file (unit)
#endif

    write(*,*) 'Time is', date_init(1), date_init(2), date_init(3), &
    date_init(4), date_init(5), date_init(6)
    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, driver_label_IS, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! SetServices for OCN
    call ESMF_GridCompSetServices(is%wrap%modelComp(1), ocean_solo_SS, userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) &
      return  ! bail out

    ! set the model clock
    call ESMF_TimeIntervalSet(timeStep, s=dt_cpld, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_TimeSet(startTime, yy=date_init(1), mm=date_init(2), &
      dd=date_init(3), h=date_init(4), m=date_init(5), s=date_init(6), &
      calkindflag=fms2esmf_cal(calendar), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_TimeIntervalSet(timeInterval, mm=months, d=days, h=hours, &
        m=minutes, s=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    stopTime = startTime + timeInterval

    internalClock = ESMF_ClockCreate(name="Application Clock", &
      timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
      
    call ESMF_GridCompSet(gcomp, clock=internalClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
  end subroutine

end module esm_mod

program ocean_solo
!
!This is an updated ocean solo model application using NUOPC framework. 
!
!Author: Fei Liu
!Thu Feb  7 10:09:06 MST 2013
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
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type

  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end
  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_types_mod,          only: ice_ocean_boundary_type


  use ESMF
  use time_utils_mod
  use esm_mod

  implicit none


  !-------------------- main code starts here ---------------------
  integer                 :: rc, urc, mpi_comm_mom
  type(ESMF_VM)           :: vm
  type(ESMF_GridComp)     :: esmComp
  
  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_VMGet(vm, mpiCommunicator=mpi_comm_mom, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call fms_init(mpi_comm_mom)
  
  call ESMF_LogWrite("MOM5 ocean solo STARTING", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Create the earth system Component
  esmComp = ESMF_GridCompCreate(name="mom_solo", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  ! SetServices for the earth system Component
  call ESMF_GridCompSetServices(esmComp, SetServices, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  ! Call Initialize for the earth system Component
  call ESMF_GridCompInitialize(esmComp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  ! Call Run  for earth the system Component
  call ESMF_GridCompRun(esmComp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  ! Call Finalize for the earth system Component
  call ESMF_GridCompFinalize(esmComp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  ! Destroy the earth system Component
  call ESMF_GridCompDestroy(esmComp, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  call ESMF_LogWrite("MOM5 ocean solo FINISHED", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Finalize ESMF
  call ESMF_Finalize()

end program ocean_solo
