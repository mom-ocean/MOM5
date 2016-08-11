                    module aerosol_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  smf
! </REVIEWER>
! <OVERVIEW>
!  Code to initialize/allocate aerosol climatology
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes prescribed aerosol climatology from input file,
!  allocates necessary memory space and interpolate the aerosol climatology
!  to the model specification. Afterwards the memory space is deallocated,
!  the aerosol climatology information is freed.
! </DESCRIPTION>
!  shared modules:

use time_manager_mod,  only: time_type, time_manager_init, operator(+),&
                             set_date, operator(-), print_date, &
                             assignment(=), &
                             set_time, days_in_month, get_date, &
                             operator(>), operator(/=)
use diag_manager_mod,  only: diag_manager_init, get_base_time, &
                             send_data, register_diag_field,  &
                             register_static_field
use field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod,only: get_tracer_index,   &
                             get_tracer_names,   &
                             get_tracer_indices, &
                             get_number_tracers, &
                             MAX_TRACER_FIELDS,  &
                             query_method
use mpp_mod,           only: input_nml_file
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use interpolator_mod,  only: interpolate_type, interpolator_init, &
                             interpolator, interpolator_end, &
                             obtain_interpolator_time_slices, &    
                             unset_interpolator_time_flag, &
                             CONSTANT, INTERP_WEIGHTED_P
use mpp_io_mod,        only: mpp_open, mpp_close, MPP_RDONLY,   &
                             MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI,  &
                             MPP_SINGLE, mpp_io_init
use constants_mod,     only: constants_init, RADIAN, GRAV
use data_override_mod, only: data_override

!  shared radiation package modules:

use   rad_utilities_mod, only  : aerosol_type, rad_utilities_init, &
                                 get_radiative_param,              &
                                 atmos_input_type     

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    aersosol_mod provides aerosol information that is needed by a
!    model physics package. the initial use of aerosol_mod was to
!    provide climatological aerosol fields to the model radiation 
!    package for use in calculating radiative fluxes and heating rates; 
!    with the introduction of predicted aerosols as members of the 
!    model's tracer array, aerosol_mod became the mechanism to collect
!    and bundle those tracers which were to be seen as aerosol by the
!    radiation code. the introduction of the treatment of aerosol 
!    impacts on cloud properties (aerosol indirect effect) required that
!    aerosol_mod be modified to provide needed aerosol information to 
!    routines involved with cloud calculation.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id: aerosol.F90,v 20.0 2013/12/13 23:18:53 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


!-----------------------------------------------------------------------
!------  interfaces -------

public          &
        aerosol_init, aerosol_driver, aerosol_end, &
        aerosol_time_vary, aerosol_endts, aerosol_dealloc

!private         &


!---------------------------------------------------------------------
!------  namelist ------

 character(len=32)      ::      &
         aerosol_data_source = 'climatology'
                                   ! source of aerosol data, either
                                   ! 'climatology' file (default) or
                                   ! single column 'input' file or
                                   ! calculate a column for location and
                                   ! time specified ('calculate_column')
                                   ! or 'predicted' (calculated online)
integer, parameter     ::      &
        MAX_DATA_FIELDS = 100     ! maximum number of aerosol species
integer, parameter     ::      &
        MAX_AEROSOL_FAMILIES = 12 ! maximum number of aerosol families
character(len=64)      ::      &
        data_names(MAX_DATA_FIELDS) = '  ' 
                                  ! names of active aerosol species 
character(len=64)      ::      &
        filename = '  '           ! name of netcdf file containing 
                                  ! aerosol species to be activated
character(len=64)      ::      &
      family_names(MAX_AEROSOL_FAMILIES) = '  ' 
                                  ! names of active aerosol families 
logical, dimension(MAX_DATA_FIELDS) :: in_family1 = .false.
                                  ! aerosol n is in family 1 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family2 = .false.
                                  ! aerosol n is in family 2 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family3 = .false.
                                  ! aerosol n is in family 3 ?
logical,dimension(MAX_DATA_FIELDS) :: in_family4 = .false.
                                  ! aerosol n is in family 4 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family5 = .false.
                                  ! aerosol n is in family 5 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family6 = .false.
                                  ! aerosol n is in family 6 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family7 = .false.
                                  ! aerosol n is in family 7 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family8 = .false.
                                  ! aerosol n is in family 8 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family9 = .false.
                                  ! aerosol n is in family 9 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family10 = .false.
                                  ! aerosol n is in family 10 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family11 = .false.
                                  ! aerosol n is in family 11 ?
logical, dimension(MAX_DATA_FIELDS) :: in_family12 = .false.
                                  ! aerosol n is in family 12 ?
logical         :: use_aerosol_timeseries = .false.
                                  ! use a timeseries providing inter-
                                  ! annual aerosol variation ?
logical, dimension(MAX_DATA_FIELDS) :: time_varying_species = .true.
                                  ! this aerosol species is 
                                  ! time-varying  ?
integer, dimension(6,MAX_DATA_FIELDS) :: aerosol_dataset_entry  =  1
                      ! time in aerosol data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
logical, dimension(MAX_AEROSOL_FAMILIES) ::   &
                                    volc_in_fam_col_opt_depth = .false.
                      ! is the volcanic contribution to column optical
                      ! depth to be included for this family in the
                      ! netcdf output fields ?
real,dimension(2)   :: lonb_col = (/-999., -999./)
                      ! longitudes defining the region to use for column
                      ! data calculation
real,dimension(2)   :: latb_col = (/-999., -999./)
                      ! latitudes defining the region to use for column
                      ! data calculation
integer, dimension(6)  :: time_col = (/0,0,0,0,0,0/)
                      ! time to use for column data calculation


namelist /aerosol_nml/                            &
                           aerosol_data_source,   &
                           lonb_col, latb_col, time_col, &
                           data_names, filename,  &
                           family_names,   &
                           use_aerosol_timeseries, &
                           time_varying_species,  &
                           aerosol_dataset_entry,  &               
                           in_family1, in_family2, in_family3, &
                           in_family4, in_family5, in_family6, &
                           in_family7, in_family8, in_family9, &
                           in_family10, in_family11, in_family12, &
                           volc_in_fam_col_opt_depth
                           

!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----


!-------------------------------------------------------------------
!    specified_aerosol contains the column input aerosol concentration
!    ratio (kg/m**2).  used when aerosol_data_source = 'input'.
!-------------------------------------------------------------------
real, dimension (:), allocatable     ::  specified_aerosol
 
!---------------------------------------------------------------------
!   the following is an interpolate_type variable containing the
!   information about the aerosol species.
!---------------------------------------------------------------------
type(interpolate_type), dimension(:), allocatable  :: Aerosol_interp

!--------------------------------------------------------------------
!    miscellaneous variables
!--------------------------------------------------------------------
integer  :: id                               ! number of grid points in 
                                             ! x direction (on processor)
integer  :: jd                               ! number of grid points in 
                                             ! y direction (on processor)
logical  :: make_separate_calls=.false.      ! aerosol interpolation
                                             ! to be done one at a 
                                             ! time
type(time_type), dimension(:), allocatable   ::    &
                   Aerosol_time              ! time for which data is
                                             ! obtained from aerosol
                                             ! timeseries
logical  :: do_column_aerosol = .false.      ! using single column aero-
                                             ! sol data ?
logical  :: do_predicted_aerosol = .false.   ! using predicted aerosol fields?
logical  :: do_specified_aerosol = .false.   ! using specified aerosol fields 
                                             ! from a timeseries file?
integer  :: nfields=0                        ! number of active aerosol 
                                             ! species
integer  :: nfamilies=0                      ! number of active aerosol 
                                             ! families
logical  :: module_is_initialized = .false.  ! module has been 
                                             ! initialized  ?

type(time_type) :: Model_init_time  ! initial calendar time for model  
                                    ! [ time_type ]
type(time_type), dimension(:), allocatable ::   &
                   Aerosol_offset   ! difference between model initial
                                    ! time and aerosol timeseries app-
                                    ! lied at model initial time
                                    ! [ time_type ]
type(time_type), dimension(:), allocatable ::   &
                   Aerosol_entry    ! time in aerosol timeseries which
                                    ! is mapped to model initial time
                                    ! [ time_type ]
type(time_type)  ::   &
              Aerosol_column_time   ! time for which aerosol data is
                                    ! extracted from aerosol timeseries
                                    ! in 'calculate_columns' case
                                    ! [ time_type ]
logical, dimension(:), allocatable    ::     &
                   negative_offset 
                                    ! the model initial time is later 
                                    ! than the aerosol_dataset_entry 
                                    ! time  ?
integer , dimension(:), allocatable :: data_out_of_bounds, vert_interp
logical, dimension(:), allocatable :: using_fixed_year_data 
                                    ! are we using a fixed year
                                    !  of data from a timeseries file ?
integer, dimension (MAX_DATA_FIELDS) :: aerosol_tracer_index
                                    ! tracer index for each of the 
                                    ! prognostic tracer to be seen as 
                                    ! aerosols by the radiation package
real, dimension (MAX_DATA_FIELDS)    :: aerosol_tracer_scale_factor
                                    ! scaling factor for each of the 
                                    ! prognostic tracer to be seen as 
                                    ! aerosols by the radiation package
character(len=32), dimension (:),  &
                     allocatable     ::  tracer_names
logical, dimension(:), allocatable   :: being_overridden
                                    ! is a given aerosol field to be over-
                                    ! ridden based on model data_table?
logical                              :: output_override_info = .true.
                                    ! should override info about each 
                                    ! aerosol field be output (will be set 
                                    ! to .false. after first time step)
integer                              :: override_counter = 0
                                    ! used to count calls to aerosol_endts
                                    ! so that output_override_info may be
                                    ! set .false. after physics_up.

#include <netcdf.inc>
 
!---------------------------------------------------------------------
!---------------------------------------------------------------------


                           contains


 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------
! <SUBROUTINE NAME="aerosol_init">
!  <OVERVIEW>
!   Subroutine to initialize/interpolate prescribed aerosol climatology
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize/interpolate prescribed aerosol climatology
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_init(lonb, latb, aerosol_names)
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
!   2d array of model longitudes on cell corners in [radians]
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes on cell corners in [radians]
!  </IN>
!  <IN NAME="aerosol_names" TYPE="character">
!   names of the activated aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_init (lonb, latb, aerosol_names,   &
                         aerosol_family_names)

!-----------------------------------------------------------------------
!    aerosol_init is the constructor for aerosol_mod.
!-----------------------------------------------------------------------

real, dimension(:,:),            intent(in)  :: lonb,latb
character(len=64), dimension(:), pointer     :: aerosol_names
character(len=64), dimension(:), pointer     :: aerosol_family_names

!----------------------------------------------------------------------
!  intent(in) variables:
!
!       lonb           2d array of model longitudes on cell corners 
!                      [ radians ]
!       latb           2d array of model latitudes at cell corners 
!                      [ radians ]
!
!   pointer variables:
!
!       aerosol_names  names of the activated aerosol species
!       aerosol_family_names  names of the activated aerosol families
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:
      
      character(len=64)  :: data_names_predicted (MAX_DATA_FIELDS) = '  '
                              ! predicted aerosol names to be 
                              ! seen by radiation code
      logical :: flag,rad_forc_online, single_year_file
      character(len=80)       ::tr_rad_name, tr_clim_name
      character(len=80)       :: name,control
      real                    ::tr_rad_scale_factor
      integer   ::   unit, ierr, io, logunit
      integer   ::   ntrace
      integer   ::   n

!---------------------------------------------------------------------
!    local variables:
!
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!         n          do-loop index
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
      call mpp_io_init
      call fms_init
      call diag_manager_init
      call rad_utilities_init
      call time_manager_init
      call constants_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=aerosol_nml, iostat=io)
      ierr = check_nml_error(io,'aerosol_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_nml, iostat=io,  &
               end=10)
        ierr = check_nml_error(io,'aerosol_nml')
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
                        write (logunit, nml=aerosol_nml)

!---------------------------------------------------------------------
!    define the dimensions of the local processors portion of the grid.
!---------------------------------------------------------------------
      id    = size(lonb,1) - 1
      jd    = size(latb,2) - 1

!---------------------------------------------------------------------
!    case of single input aerosol field. when running standalone code
!    on other than FMS level structure, aerosol_data_source must be 
!    'input'.
!---------------------------------------------------------------------
      if (trim(aerosol_data_source)== 'input') then
        do_column_aerosol = .true.
        call obtain_input_file_data
        nfields = 1
        allocate (aerosol_names(nfields))
        aerosol_names (1) = 'total_aerosol'

!---------------------------------------------------------------------
!    case of predicted aerosols.
!---------------------------------------------------------------------
      else if (trim(aerosol_data_source) == 'predicted') then
        do_predicted_aerosol = .true.

!-----------------------------------------------------------------------
!    count number of activated aerosol species, which will be carried
!    in the model as tracers with an attribute of 'radiative_param'. 
!    define the names associated with these aerosols.
!-----------------------------------------------------------------------
        call get_number_tracers(MODEL_ATMOS, num_tracers= ntrace)
        allocate (tracer_names(ntrace))
        do n = 1, ntrace
          call get_tracer_names(MODEL_ATMOS,n,tracer_names(n))
          flag = query_method ('radiative_param', MODEL_ATMOS, &
                               n, name, control)
          if (flag) then
            call get_radiative_param(name,control,rad_forc_online, &
                                     tr_rad_name, tr_clim_name, &
                                     tr_rad_scale_factor)
            if (rad_forc_online) then
              nfields = nfields +1
              aerosol_tracer_index(nfields) = n
              data_names_predicted(nfields) = trim(tr_rad_name)
!             data_names(nfields)           = trim(tr_clim_name)
              data_names(nfields)           = trim(tr_rad_name)
              aerosol_tracer_scale_factor(nfields) = tr_rad_scale_factor
            endif
          endif
        end do

!----------------------------------------------------------------------
!    allocate and fill pointer arrays to return the names of the activ-
!    ated species and any activated families to the calling routine.
!---------------------------------------------------------------------
        allocate (aerosol_names(nfields))
        aerosol_names(:)        = data_names_predicted(1:nfields)

!----------------------------------------------------------------------
!    allocate and fill an array to indicate whether or not each aerosol 
!    field is  to be overridden.
!---------------------------------------------------------------------
        allocate (being_overridden(nfields))
        being_overridden(:) = .false.

!---------------------------------------------------------------------
!    case of 'climatology' and 'calculate_column' aerosol data source.
!---------------------------------------------------------------------
      else   ! (trim(aerosol_data_source) == 'input')     
        do_specified_aerosol = .true.

!---------------------------------------------------------------------
!    determine how many aerosols in the file are activated.
!--------------------------------------------------------------
        do n=1,MAX_DATA_FIELDS
          if (data_names(n) /= ' '  ) then
            nfields = n
          else
            exit
          endif
        end do
  
!---------------------------------------------------------------------
!    check for case of inconsistent nml specification -- case of re-
!    questing time-varying aerosol for a given aerosol, but indicating
!    that the time series is not to be used. in this case, a note will
!    be written to stdout indicating that the aerosol species will not
!    be time-varying. this is needed because the default nml settings  
!    are defined so as to allow backward compatibility with existing 
!    code and script settings, and lead to this conflict.
!---------------------------------------------------------------------
        do n=1, nfields
          if (.not. use_aerosol_timeseries) then
            if (time_varying_species(n)) then
              call error_mesg ('aerosol_mod', &
                  'inconsistent nml settings -- not using aerosol  &
                  &timeseries but requesting interannual variation of  &
                  & aerosol amount for '  // trim (data_names(n))  // &
                 ' -- this aerosol will NOT exhibit interannual &
                  &variation', NOTE)
              time_varying_species(n) = .false.
            endif
          endif
        end do

!---------------------------------------------------------------------
!    allocate and fill pointer arrays to return the names of the activ-
!    ated species and any activated families to the calling routine.
!--------------------------------------------------------------------
        allocate (aerosol_names(nfields))
        aerosol_names (:) = data_names(1:nfields)

!----------------------------------------------------------------------
!    allocate and initialize module variables.
!----------------------------------------------------------------------
        allocate (Aerosol_offset(nfields), Aerosol_entry(nfields), &
              negative_offset(nfields), using_fixed_year_data(nfields)) 
        Aerosol_offset = set_time (0,0)
        Aerosol_entry = set_time (0,0)
        negative_offset = .false.
        using_fixed_year_data = .false.

!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table) 
!----------------------------------------------------------------------
        Model_init_time = get_base_time()

!----------------------------------------------------------------------
!    define the array using_fixed_year_data. it will be .true. for a 
!    given aerosol species if the nml variable use_aerosol_timeseries 
!    is .true., and the nml variable time_varying_species for that
!    aerosol is .false., or if use_aerosol_timeseries is .false but a 
!    non-default aerosol_dataset_entry has been specified; otherwise it 
!    will be .false..
!----------------------------------------------------------------------
        do n=1,nfields           
          if (use_aerosol_timeseries) then
            if (time_varying_species(n)) then
              using_fixed_year_data(n) = .false.
            else
              using_fixed_year_data(n) = .true.
            endif
  
!---------------------------------------------------------------------
!    if no dataset entry point is supplied when an aerosol timeseries
!    file is being used, define the entry point as the model base time.
!---------------------------------------------------------------------
            if (aerosol_dataset_entry(1,n) == 1 .and. &
                aerosol_dataset_entry(2,n) == 1 .and. &
                aerosol_dataset_entry(3,n) == 1 .and. &
                aerosol_dataset_entry(4,n) == 1 .and. &
                aerosol_dataset_entry(5,n) == 1 .and. &
                aerosol_dataset_entry(6,n) == 1 ) then
              Aerosol_entry(n) = Model_init_time

!----------------------------------------------------------------------
!    if a dataset entry time is defined, compute the offset from model 
!    base time to aerosol_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
            else
              Aerosol_entry(n) = set_date (aerosol_dataset_entry(1,n), &
                                           aerosol_dataset_entry(2,n), &
                                           aerosol_dataset_entry(3,n), &
                                           aerosol_dataset_entry(4,n), &
                                           aerosol_dataset_entry(5,n), &
                                           aerosol_dataset_entry(6,n))
            endif

!----------------------------------------------------------------------
!    indicate that aerosol species n will be defined from the timeseries
!    file, and the relationship of the timeseries to the model calendar.
!--------------------------------------------------------------------
            call error_mesg ( 'aerosol_mod', &
               'PROCESSING AEROSOL TIMESERIES FOR ' // &
               trim(aerosol_names(n)), NOTE)
            call print_date (Aerosol_entry(n) ,   &
                str= ' Data from aerosol timeseries at time: ')
            call print_date (Model_init_time , str=' This data is &
                               &mapped to model time:')

!---------------------------------------------------------------------
!    indicate whether a single year of the aerosol climatology will be
!    repeated throughout the model run, or if the aerosol time behavior
!    will show interannual changes.
!---------------------------------------------------------------------
            if (using_fixed_year_data(n)) then
              call error_mesg ('aerosol_mod', &
                 'This annual cycle will be used every model year &
               & -- no interannual variation for '  &
                                     // trim(aerosol_names(n)), NOTE)
            else
              call error_mesg ('aerosol_mod', &
                      trim(aerosol_names(n)) //   &
                     ' will exhibit interannual variation as defined &
                     & in the climatology file ', NOTE)
            endif
!---------------------------------------------------------------------
!    define the offset between the aerosol timeseries and the model
!    calendar, and whether this is a positive or negative offset.
!--------------------------------------------------------------------
            Aerosol_offset(n) = Aerosol_entry(n) - Model_init_time
            if (Model_init_time > Aerosol_entry(n)) then
              negative_offset(n) = .true.
            else
              negative_offset(n) = .false.
            endif

!---------------------------------------------------------------------
!    if use_aerosol_timeseries is .false., then either data from a 
!    single year defined in a timeseries file is to be used throughout 
!    the model integration, or a non-specific single-year aerosol 
!    climatology file is to be used.
!---------------------------------------------------------------------
          else 

!---------------------------------------------------------------------
!    if no dataset entry has been specified, then a non-specific single
!    year climatology file is being used. set the variable 
!    using_fixed_year_data to be .false.. output a descriptive message
!    to stdout.
!---------------------------------------------------------------------
            if (aerosol_dataset_entry(1,n) == 1 .and. &
                aerosol_dataset_entry(2,n) == 1 .and. &
                aerosol_dataset_entry(3,n) == 1 .and. &
                aerosol_dataset_entry(4,n) == 1 .and. &
                aerosol_dataset_entry(5,n) == 1 .and. &
                aerosol_dataset_entry(6,n) == 1 ) then
              using_fixed_year_data(n) = .false.
              if (mpp_pe() == mpp_root_pe() ) then
                print *, 'Aerosol data for ', trim(aerosol_names(n)),   &
                   ' obtained from single year climatology file '
              endif

!---------------------------------------------------------------------
!    if a year has been specified for the dataset entry, then the data
!    will be coming from an aerosol timeseries file, but the same annual
!    aerosol variation will be used for each model year. set the var-
!    iable using_fixed_year_data to be .true..  define Aerosol_entry as
!    feb 01 of the year given by the first element of nml variable 
!    aerosol_dataset_entry. output a descriptive message to stdout.
!--------------------------------------------------------------------
            else
              using_fixed_year_data(n) = .true.
              Aerosol_entry(n) = set_date (aerosol_dataset_entry(1,n), &
                                           2, 1, 0, 0, 0)
              call error_mesg ('aerosol_mod', &
                  'Aerosol data is defined from a single annual cycle &
                  &for ' // trim(aerosol_names(n)) //   &
                  &' - no interannual variation', NOTE)
              if (mpp_pe() == mpp_root_pe() ) then
                print *, 'Aerosol data for ', trim(aerosol_names(n)),  &
                    ' obtained from aerosol timeseries &
                    &for year:', aerosol_dataset_entry(1,n)
              endif
            endif
          endif
        end do

!-----------------------------------------------------------------------
!    count number of activated aerosol families. allocate a pointer 
!    array to return the names of the activated species to the calling 
!    routine.
!-----------------------------------------------------------------------
!       do n=1,MAX_AEROSOL_FAMILIES
!         if (family_names(n) /= ' '  ) then
!           nfamilies = n
!         else
!           exit
!         endif
!       end do

!-----------------------------------------------------------------------
!    allocate and initialize variables needed for interpolator_mod if 
!    any aerosol species have been activated.
!-----------------------------------------------------------------------
        allocate (data_out_of_bounds(nfields))
        allocate (vert_interp       (nfields))
        data_out_of_bounds = CONSTANT
        vert_interp = INTERP_WEIGHTED_P

!----------------------------------------------------------------------
!    determine if separate calls to interpolator  must be made for 
!    each aerosol species, or if all variables in the file may be
!    interpolated together.  reasons for separate calls include differ-
!    ent data times desired for different aerosols, different vertical
!    interpolation procedures and different treatment of undefined
!    data.
!----------------------------------------------------------------------
          do n=2,nfields
            if (time_varying_species(n) .and.   &
               (.not. time_varying_species(n-1) ) ) then
              make_separate_calls = .true.
              exit
            endif
            if (using_fixed_year_data(n) .and.   &
                (.not. using_fixed_year_data(n-1) ) ) then
              make_separate_calls = .true.
              exit
            endif
            if (Aerosol_entry(n) /= Aerosol_entry(n-1)) then
              make_separate_calls = .true.
              exit
            endif
            if (data_out_of_bounds(n) /= data_out_of_bounds(n-1)) then
              make_separate_calls = .true.
              exit
            endif
            if (vert_interp       (n) /= vert_interp       (n-1)) then
              make_separate_calls = .true.
              exit
            endif
          end do
          if (make_separate_calls) then
            allocate (Aerosol_interp(nfields))  
            allocate (Aerosol_time  (nfields))  
          else
            allocate (Aerosol_interp(1))  
            allocate (Aerosol_time  (1))  
          endif

!----------------------------------------------------------------------
!    determine if the aerosol_data_source is specified as 
!    'calculate_column'. 
!--------------------------------------------------------------------
        if (trim(aerosol_data_source) == 'calculate_column') then

!----------------------------------------------------------------------
!    if the aerosol_data_source is specified as 'calculate_column', then
!    the aerosol fields will be obtained by averaging the aerosol fields
!    in the climatology over a specified latitude-longitude section at
!    a specified calendar time, and this profile will be used in all 
!    model columns. make sure the specified lats / lons / time are 
!    valid.
!-------------------------------------------------------------------
          do n=1,2
            if (lonb_col(n) < 0. .or. lonb_col(n) > 360.) then
              call error_mesg ('aerosol_mod', &
                  ' invalid value for lonb_col', FATAL)
            endif
            if (latb_col(n) < -90. .or. latb_col(n) > 90.) then
              call error_mesg ('aerosol_mod', &
                  ' invalid value for latb_col', FATAL)
            endif
          end do
          if (time_col(1) == 0) then
            call error_mesg ('aerosol_mod', &
                'invalid time specified for time_col', FATAL)
          endif

          if (.not. use_aerosol_timeseries) then
              call error_mesg ('aerosol_mod', &
         'must use_aerosol_timeseries when calculate_column is .true.', FATAL)
          endif 

          if (any(time_varying_species(1:nfields))) then
              call error_mesg ('aerosol_mod', &
                   'aerosol values must be fixed in time when &
                                   &calculate_column is .true.', FATAL)
          endif 


!----------------------------------------------------------------------
!    call interpolator_init to begin processing the aerosol climat-
!    ology file. define the valid time as a time_type, and output
!    informative messages.
!---------------------------------------------------------------------
    if (make_separate_calls) then
              call error_mesg ('aerosol_mod', &
         'make_separate_calls not allowed  for calculate_column', FATAL)
     else       
          call interpolator_init (Aerosol_interp(1)    , filename,  &
                                    spread(lonb_col/RADIAN,2,2),  &
                                    spread(latb_col/RADIAN,1,2),&
                                    data_names(:nfields),   &
                                    data_out_of_bounds=   &
                                                  data_out_of_bounds, &
                                    vert_interp=vert_interp,  &
                                    single_year_file = single_year_file)
     endif
          Aerosol_column_time = set_date (time_col(1), time_col(2), &
                                          time_col(3), time_col(4), &
                                          time_col(5), time_col(6))
          call print_date (Aerosol_column_time,  str= &
              ' Aerosol data used is from aerosol timeseries at time: ')
          if (mpp_pe() == mpp_root_pe() ) then
            print *, 'Aerosol data is averaged over latitudes',  &
                latb_col(1), ' to', latb_col(2), ' and longitudes',&
                lonb_col(1), ' to', lonb_col(2)
          endif

!----------------------------------------------------------------------
!    if 'calculate_column' is .false., then the aerosol fields will have
!    the appropriate horizontal variation. call interpolator_init to
!    begin processing the aerosol climatology file.    
!-------------------------------------------------------------------
        else  ! (calculate_column)
    if (make_separate_calls) then
       do n=1,nfields
          call interpolator_init (Aerosol_interp(n), filename, lonb, &
                                  latb, data_names(n:n     ),   &
                                  data_out_of_bounds=    &
                                                  data_out_of_bounds(n:n), &
                                  vert_interp=vert_interp(n:n), &
                                  single_year_file=single_year_file)
       end do
     else
          call interpolator_init (Aerosol_interp(1), filename, lonb, &
                                  latb, data_names(:nfields),   &
                                  data_out_of_bounds=    &
                                                  data_out_of_bounds, &
                                  vert_interp=vert_interp, &
                                  single_year_file=single_year_file)
     endif
        endif ! (calculate_column)

!--------------------------------------------------------------------
!    check for compatibility of nml options requested and the aerosol
!    data file which was read.
!--------------------------------------------------------------------
        if (single_year_file .and.  use_aerosol_timeseries) then
          call error_mesg ('aerosol_mod', &
               'aerosol input file is single-year, yet interannual &
                &variation of aerosol is requested', FATAL  )
        endif
        do n=1, nfields
          if (.not. use_aerosol_timeseries .and.   &
              .not. using_fixed_year_data(n) .and. &
              .not. single_year_file)  then
            call error_mesg ('aerosol_mod', &
              'aerosol input file contains a time-series, yet nml  &
               &settings indicate that a non-specific single-year &
               &climatology is to be used', FATAL)
          endif
          if (.not. use_aerosol_timeseries .and.   &
              using_fixed_year_data(n) .and. &
              single_year_file)  then
            call error_mesg ('aerosol_mod', &
             'aerosol input file is non-specific single-year file,  &
               &yet nml settings specify that a particular single-year &
               &climatology is to be used', FATAL)
          endif
        end do
      endif  ! ('aerosol_data_source == 'input')

!-----------------------------------------------------------------------
!    count number of activated aerosol families. allocate a pointer 
!    array to return the names of the activated species to the calling 
!    routine.
!-----------------------------------------------------------------------
      do n=1,MAX_AEROSOL_FAMILIES
        if (family_names(n) /= ' '  ) then
          nfamilies = n
        else
          exit
        endif
      end do

!---------------------------------------------------------------------
!    allocate and fill pointer arrays to return the names of any activ-
!    ated families to the calling routine.
!-------------------------------------------------------------------
      allocate (aerosol_family_names(nfamilies))
      aerosol_family_names (:) = family_names(1:nfamilies)

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------



end subroutine aerosol_init

!############################################################################

subroutine aerosol_time_vary (model_time)

!--------------------------------------------------------------------------- 
!   subroutine aerosol_time_vary makes sure the aerosol interpolate_type 
!   variable has access to the proper time levels of data in the aerosol
!---------------------------------------------------------------------------

type(time_type), intent(in) :: model_time


      integer :: n

!----------------------------------------------------------------------
!    be sure the proper time levels are in memory for the aerosol timeseries.
!----------------------------------------------------------------------
      if ( do_specified_aerosol) then
 
        if (make_separate_calls) then
!--------------------------------------------------------------------
!    if separate calls are required for each aerosol species, loop over
!    the individual species.
!--------------------------------------------------------------------
          do n=1,nfields

!--------------------------------------------------------------------
!    if the data timeseries is to be used for species n, define the
!    time for which data is desired, and then call interpolator to
!    verify the time levels bracketing the desired time are available.
!--------------------------------------------------------------------
            if (use_aerosol_timeseries) then
              if (time_varying_species(n)) then

!----------------------------------------------------------------------
!    define the Aerosol_time for aerosol n and check for the  
!    appropriate time slices.
!----------------------------------------------------------------------
                if (negative_offset(n)) then
                  Aerosol_time(n) = model_time - Aerosol_offset(n)
                else
                  Aerosol_time(n) = model_time + Aerosol_offset(n)
                endif
                call obtain_interpolator_time_slices (Aerosol_interp(n), &
                                                          Aerosol_time(n))     
              else
                call set_aerosol_time (model_time, Aerosol_entry(n), &
                                       Aerosol_time(n))
                call obtain_interpolator_time_slices (Aerosol_interp(n), &
                                                       Aerosol_time(n))     
              endif

!--------------------------------------------------------------------
!    if the data timeseries is not to be used for species n, define the
!    time for which data is desired, and then call interpolator to
!    obtain the data.
!--------------------------------------------------------------------
            else  ! (use_aerosol_timeseries)

!---------------------------------------------------------------------
!    if a fixed year has not been specified, obtain data relevant for
!    the current model year.
!---------------------------------------------------------------------
              if ( .not. using_fixed_year_data(n)) then
                call obtain_interpolator_time_slices (Aerosol_interp(n), &
                                                              model_time)     
                Aerosol_time(n) = model_time

!----------------------------------------------------------------------
!    if a fixed year has been specified, call set_aerosol_time to define
!    the Aerosol_time to be used for aerosol n. call interpolator to 
!    obtain the aerosol values and store the aerosol amounts in 
!    Aerosol%aerosol.
!----------------------------------------------------------------------
              else 
                call set_aerosol_time (model_time, Aerosol_entry(n), &
                                         Aerosol_time(n))
                call obtain_interpolator_time_slices (Aerosol_interp(n), &
                                                          Aerosol_time(n))     
              endif  
            endif ! (use_aerosol_timeseries)
          end do  !(nfields)

        else  ! (make_separate_calls)
           
!--------------------------------------------------------------------
!    if the data timeseries is to be used for species n, define the
!    time for which data is desired.
!--------------------------------------------------------------------
          if (use_aerosol_timeseries) then
            if (negative_offset(1)) then
              Aerosol_time(1) = model_time - Aerosol_offset(1)
            else
              Aerosol_time(1) = model_time + Aerosol_offset(1)
            endif

!--------------------------------------------------------------------
!    if 'calculate_column' is being used,  be sure the needed time slices
!    are in memory.
!--------------------------------------------------------------------
            if (trim(aerosol_data_source) == 'calculate_column') then
              call obtain_interpolator_time_slices (Aerosol_interp(1), &
                                Aerosol_column_time)
            else

!-------------------------------------------------------------------
!    since separate calls are not required, all aerosol species are 
!    either time_varying  or not.  be sure the needed time slices are available.
!--------------------------------------------------------------------
              if (time_varying_species(1)) then
                call obtain_interpolator_time_slices (Aerosol_interp(1), &
                                Aerosol_time(1))
              else
                call set_aerosol_time (model_time, Aerosol_entry(1), &
                                       Aerosol_time(1))
                call obtain_interpolator_time_slices (Aerosol_interp(1), &
                                Aerosol_time(1))
              endif    
            endif  ! (calculate_column)

!--------------------------------------------------------------------
!    if the data timeseries is not to be used for species n, define the
!    time for which data is desired, and verify the presence of the needed
!    bracketing time slices.
!--------------------------------------------------------------------
          else ! (use_aerosol_timeseries)

!---------------------------------------------------------------------
!    if a fixed year has not been specified, obtain data relevant for
!    the current model year. this data comes from a non-specific single-yr
!    climatology file.
!---------------------------------------------------------------------
            if (.not. using_fixed_year_data(1)) then
              call obtain_interpolator_time_slices (Aerosol_interp(1), &
                                  model_time)
              Aerosol_time(1) = model_time

!----------------------------------------------------------------------
!    if a fixed year has been specified, call set_aerosol_time to define
!    the Aerosol_time, then verify the needed time slices are available. 
!----------------------------------------------------------------------
            else
              call set_aerosol_time (model_time, Aerosol_entry(1), &
                                                             Aerosol_time(1))
              call obtain_interpolator_time_slices (Aerosol_interp(1), &
                                                            Aerosol_time(1))
            endif ! (using_fixed_year)
          endif ! (use_aerosol_timeseries)
        endif  ! (make_separate_calls   )
      endif ! (do_column_aerosol)
               
!-------------------------------------------------------------------- 



end subroutine aerosol_time_vary 



!####################################################################

subroutine aerosol_endts

     integer :: n

     if (allocated(Aerosol_interp)) then
       do n=1, size(Aerosol_interp,1)
         call unset_interpolator_time_flag (Aerosol_interp(n))
       end do
     endif
  
     override_counter = override_counter + 1
     if (override_counter == 2) then
       output_override_info = .false.
     endif

end subroutine aerosol_endts



!######################################################################
! <SUBROUTINE NAME="aerosol_driver">
!  <OVERVIEW>
!   Interpolate aerosol verical profile based on prescribed aerosol
!   climatology input and model set up.
!  </OVERVIEW>
!  <TEMPLATE>
!   call aerosol_driver (is, js, model_time, p_half, Aerosol)
!  </TEMPLATE>
!  <INOUT NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatology input
!  </INOUT>
!  <IN NAME="model_time" TYPE="time_type">
!   The internal model simulation time, i.e. Jan. 1 1982
!  </IN>
!  <IN NAME="tracer" TYPE="real">
!   4 dimensional array of tracers, last index is the number of all tracers
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   The array of model layer pressure values
!  </IN>
!  <IN NAME="is" TYPE="integer">
!   The longitude index of model physics window domain
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   The latitude index of model physics window domain
!  </IN>
!  <IN NAME="override_aerosols" TYPE="logical, optional">
!   use offline aerosols via data_override?
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_driver (is, js, model_time, tracer, &
                           p_half, p_flux, Aerosol, override_aerosols)

!-----------------------------------------------------------------------
!    aerosol_driver returns the names and concentrations of activated 
!    aerosol species at model grid points at the model_time to the 
!    calling routine in aerosol_type variable Aerosol. 
!-----------------------------------------------------------------------

integer,                  intent(in)     :: is,js
type(time_type),          intent(in)     :: model_time
real, dimension(:,:,:,:), intent(in)     :: tracer
real, dimension(:,:,:),   intent(in)  :: p_half, p_flux
type(aerosol_type),       intent(inout)  :: Aerosol
logical, optional,        intent(in)     :: override_aerosols

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       is, js           starting subdomain i,j indices of data in 
!                        the physics_window being integrated
!       model_time       time for which aerosol data is desired
!                        [ time_type ]
!       p_half           model pressure at interface levels 
!                        [ Pa ]
!      
!   intent(inout) variables:
!
!       Aerosol    aerosol_type variable. the following components will
!                  be returned from this routine:
!                   aerosol      concentration of each active aerosol 
!                                species at each model grid point
!                                [ kg / m**2 ]
!                   aerosol_names 
!                                names assigned to each active species
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension(1,1, size(p_half,3)-1,    &
                                               nfields) :: aerosol_data
      real, dimension(1,1, size(p_half,3))   :: p_half_col
      real, dimension(id,jd,size(p_half,3)-1) :: aerosol_proc
      integer         :: n, k, j, i, na            ! do-loop index
      integer         :: nn
      logical         :: do_override, used
      integer         :: ie, je

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    allocate an array to hold the activated aerosol names. allocate an
!    array which defines the members of the requested aerosol families.
!    allocate an array to hold the aerosol amounts for each species at
!    each grid point. 
!---------------------------------------------------------------------
      allocate (Aerosol%aerosol_names (nfields))
      allocate (Aerosol%family_members(nfields+1, nfamilies))
      allocate (Aerosol%aerosol(size(p_half,1),  &
                                size(p_half,2), &
                                size(p_half,3) - 1, nfields)) 
      ie = is + size(p_half,1) - 1
      je = js + size(p_half,2) - 1

      if (do_column_aerosol) then
 
!---------------------------------------------------------------------
!    here all aerosol is consolidated into a single variable.
!----------------------------------------------------------------------
        Aerosol%aerosol_names(1) = 'total_aerosol'
        do k=1, size(Aerosol%aerosol,3)
          Aerosol%aerosol(:,:,k,1) = specified_aerosol(k)
        end do
      else 
      
!--------------------------------------------------------------------
!    define an array to hold the activated aerosol names.
!---------------------------------------------------------------------
        Aerosol%aerosol_names = data_names(:nfields) 

!--------------------------------------------------------------------
!    define an array which defines the members of the requested aerosol
!    families.
!---------------------------------------------------------------------
        if (nfamilies > 0) then
          do n=1,nfamilies
            do na = 1, nfields
              select case(n)
                case (1)
                  Aerosol%family_members(na,1) = in_family1(na)
                case (2)
                  Aerosol%family_members(na,2) = in_family2(na)
                case (3)
                  Aerosol%family_members(na,3) = in_family3(na)
                case (4)
                  Aerosol%family_members(na,4) = in_family4(na)
                case (5)
                  Aerosol%family_members(na,5) = in_family5(na)
                case (6)
                  Aerosol%family_members(na,6) = in_family6(na)
                case (7)
                  Aerosol%family_members(na,7) = in_family7(na)
                case (8)
                  Aerosol%family_members(na,8) = in_family8(na)
                case (9)
                  Aerosol%family_members(na,9) = in_family9(na)
                case (10)
                  Aerosol%family_members(na,10) = in_family10(na)
                case (11)
                  Aerosol%family_members(na,11) = in_family11(na)
                case (12)
                  Aerosol%family_members(na,12) = in_family12(na)
                case DEFAULT
              end select
            end do
            if (volc_in_fam_col_opt_depth(n)) then
              Aerosol%family_members(nfields+1,n) = .true.
            else
              Aerosol%family_members(nfields+1,n) = .false.
            endif
          end do
        endif

!----------------------------------------------------------------------
         if ( do_specified_aerosol) then

!--------------------------------------------------------------------
!    if 'calculate_column' is being used, obtain the aerosol values for
!    each column, one at a time, using the pressure profile for that
!    column. this allows each column to see the same aerosol fields,
!    but distributed appropriately for its pressure structure.
!--------------------------------------------------------------------
              if (trim(aerosol_data_source) == 'calculate_column') then
                do j=1, size(p_half,2)
                  do i=1, size(p_half,1)
                    p_half_col(1,1,:) = p_flux(i,j,:)
                    call interpolator (Aerosol_interp(1),&
                                       Aerosol_column_time,  &
                                       p_half_col, aerosol_data, &
                                       Aerosol%aerosol_names(1), 1, 1  )
                    Aerosol%aerosol(i,j,:,:) = aerosol_data(1,1,:,:)
                  end do
                end do
              else

!--------------------------------------------------------------------
!    if separate calls are required for each aerosol species, loop over
!    the individual species.
!--------------------------------------------------------------------
                if (make_separate_calls) then
                  do n=1,nfields                
                    call interpolator (Aerosol_interp(n), Aerosol_time(n),  &
                                     p_flux, &
                                     Aerosol%aerosol(:,:,:,n),    &
                                     Aerosol%aerosol_names(n), is, js)
                  end do  !(nfields)

!----------------------------------------------------------------------
!    if separate calls are not required, use the first aerosol char-
!    acteristics to define Aerosol_time and make a single call to 
!    interpolator. store the aerosol amounts in Aerosol%aerosol.
!----------------------------------------------------------------------
                else      
                  call interpolator (Aerosol_interp(1), Aerosol_time(1),  &
                                    p_flux, &
                                    Aerosol%aerosol,    &
                                    Aerosol%aerosol_names(1), is, js)
                endif      
              endif ! (calculate_column)

!-------------------------------------------------------------------- 
!    for predicted aerosols (obtained from tracer array), assign the 
!    tracer to "Aerosol" if that TRACER has "radiative_param" and 
!    "online", both defined in the field_table. 
! ***********************WARNINGS**************************************
!    the tracers are assumed to be expressed in Mass Mixing Ratio (MMR), 
!    and are converted into mass column for the radiative code.
!    Conversions (e.g. OC -> OM, or SO4 -> (NH4)2SO4 ) can be done
!    via radiative_param attribute scale_factor (in field_table).
!------------------------------------------------------------------ 
        else                  ! (do_specified_aerosol')
          if (present(override_aerosols)) then
            do_override = override_aerosols
          else
            do_override = .false.
          end if

          do nn=1,nfields
            n = aerosol_tracer_index(nn)

            Aerosol%aerosol(:,:,:,nn) = tracer(:,:,:,n)
            if (do_override) then
              call data_override('ATM', TRIM(tracer_names(n))//'_aerosol',&
                                 aerosol_proc, model_time, override=used)
              if (used) then
                if (output_override_info) then
                  call error_mesg ('aerosol_mod', &
                       TRIM(tracer_names(n))//'_aerosol => '// &
                     TRIM(tracer_names(n)) // ' is being overridden', NOTE)
                  being_overridden(nn) = .true.
                endif
                Aerosol%aerosol(:,:,:,nn) = aerosol_proc(is:ie,js:je,:)
              else
                if (output_override_info) then
                  call error_mesg ('aerosol_mod', &
                    TRIM(tracer_names(n))//'_aerosol => '//  &
                       TRIM(tracer_names(n)) // ' not overridden', NOTE)
                else
                  if (being_overridden(nn)) then
                    call error_mesg ('aerosol_mod', &
                       TRIM(tracer_names(n))//'_aerosol => '//  &
                         TRIM(tracer_names(n)) // ' not overridden &
                                      &when override was requested', FATAL)
                  endif
                endif
              endif
            endif

            do k=1,size(Aerosol%aerosol,3)
              do j=1,size(Aerosol%aerosol,2)
                do i=1,size(Aerosol%aerosol,1)
                  Aerosol%aerosol(i,j,k,nn) =    &
                          MAX (0.0, Aerosol%aerosol(i,j,k,nn)) * &
                          aerosol_tracer_scale_factor(nn) * &
                          ( p_half(i,j,k+1)-p_half(i,j,k) )/GRAV
                end do
              end do
            end do
          end do
        endif   ! (do_specified_aerosol')
      endif  ! (do_column_aerosol)

!-------------------------------------------------------------------- 

end subroutine aerosol_driver



!#####################################################################
! <SUBROUTINE NAME="aerosol_end">
!  <OVERVIEW>
!   aerosol_end is the destructor for aerosol_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   aerosol_end is the destructor for aerosol_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosol_end

!----------------------------------------------------------------------
!    aerosol_end is the destructor for aerosol_mod.
!----------------------------------------------------------------------

      integer  :: n


!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    call interpolator_end to release the interpolate_type variable 
!    used in this module.
!---------------------------------------------------------------------
      if (do_specified_aerosol) then
          if (nfields > 0) then
            do n=1, size(Aerosol_interp,1)
              call interpolator_end (Aerosol_interp(n))
            end do        
          endif
          deallocate (Aerosol_time)
      endif

      if (allocated (specified_aerosol)) deallocate (specified_aerosol)
      if (allocated (Aerosol_offset   )) deallocate (Aerosol_offset   )
      if (allocated (Aerosol_entry    )) deallocate (Aerosol_entry    )
      if (allocated (negative_offset  )) deallocate (negative_offset  )
      if (allocated (data_out_of_bounds))   &
                                      deallocate (data_out_of_bounds  )
      if (allocated (vert_interp      )) deallocate (vert_interp      ) 
      if (allocated (using_fixed_year_data))  &
                                     deallocate (using_fixed_year_data)

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------



end subroutine aerosol_end


!#####################################################################

subroutine set_aerosol_time (Model_time, Entry, Aerosol_time)

type(time_type), intent(in)   :: Model_time, Entry
type(time_type), intent(out)  :: Aerosol_time

      integer :: mo_yr, yr, mo, dy, hr, mn, sc, dum, dayspmn

      call get_date (Model_time, mo_yr, mo, dy, hr, mn, sc)
      call get_date (Entry, yr, dum,dum,dum,dum,dum)
      if (mo ==2 .and. dy == 29) then
        dayspmn = days_in_month(Entry)
        if (dayspmn /= 29) then
          Aerosol_time = set_date (yr, mo, dy-1, hr, mn, sc)
        else
          Aerosol_time = set_date (yr, mo, dy, hr, mn, sc)
        endif
      else
        Aerosol_time = set_date (yr, mo, dy, hr, mn, sc)
      endif

!--------------------------------------------------------------------


end subroutine set_aerosol_time


!#####################################################################
! <SUBROUTINE NAME="aerosol_dealloc">
!  <OVERVIEW>
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_dealloc
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosol_dealloc (Aerosol)

!---------------------------------------------------------------------
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!---------------------------------------------------------------------

type(aerosol_type), intent(inout) :: Aerosol

!---------------------------------------------------------------------
!  intent(inout) variables:
! 
!      Aerosol       aerosol_type variable containing information on
!                    the activated aerosol species
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate the components of the aerosol_type variable.
!---------------------------------------------------------------------
      deallocate (Aerosol%aerosol)
      deallocate (Aerosol%aerosol_names)
      deallocate (Aerosol%family_members)
 
!----------------------------------------------------------------------


end subroutine aerosol_dealloc


!#####################################################################

! <SUBROUTINE NAME="obtain_input_file_data">
!  <OVERVIEW>
!   obtain_input_file_data reads an input file containing a single
!    column aerosol profile. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_input_file_data reads an input file containing a single
!    column aerosol profile. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  obtain_input_file_data
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine obtain_input_file_data 

!---------------------------------------------------------------------
!    obtain_input_file_data reads an input file containing a single
!    column aerosol profile.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer   :: iounit    ! unit to read file on
      integer   :: kmax_file ! number of levels of data in file
      integer   :: k         ! do-loop index
      character(len=31), dimension(200) :: dimnam
      integer(kind=4), dimension(200) :: dimsiz
      integer(kind=4)                 :: ncid, rcode, nvars, ndims, &
                                         ngatts, recdim
      integer    :: i, j
      integer, PARAMETER :: MAXDIMS = 10
      integer(kind=4), dimension(MAXDIMS) :: start, count, vdims
      integer(kind=4)                     :: ivarid, ntp, nvdim, nvs, &
                                             ndsize
      character(len=31)   dummy
      



!-------------------------------------------------------------------
!    determine if a netcdf input data file exists. if so, read the 
!    number of data records in the file.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/id1aero.nc') ) then
        ncid = ncopn ('INPUT/id1aero.nc', 0, rcode)
        call ncinq (ncid, ndims, nvars, ngatts, recdim, rcode)
        do i=1,ndims
          call ncdinq (ncid, i, dimnam(i), dimsiz(i), rcode)
          if (dimnam(i) == 'lev') then
            kmax_file = dimsiz(i)
          endif
        end do
             
!-------------------------------------------------------------------
!    allocate space for the input data. read the data set. close the 
!    file upon completion.
!---------------------------------------------------------------------
        allocate (specified_aerosol(kmax_file) )
        ivarid = ncvid(ncid, 'aerosol', rcode)
        call ncvinq (ncid, ivarid, dummy, ntp, nvdim, vdims, nvs, rcode)
        do j=1,nvdim
          call ncdinq (ncid, vdims(j), dummy, ndsize, rcode)
          start(j) = 1
          count(j) = ndsize
        end do
       call ncvgt (ncid, ivarid, start, count, specified_aerosol, rcode)

         call ncclos (ncid, rcode)

!-------------------------------------------------------------------
!    determine if the input data input file exists in ascii format. if 
!    so, read the number of data records in the file.
!---------------------------------------------------------------------
      else if (file_exist ( 'INPUT/id1aero') ) then
        iounit = open_namelist_file ('INPUT/id1aero')
        read (iounit,FMT = '(i4)') kmax_file

!-------------------------------------------------------------------
!    allocate space for the input data. read the data set. close the 
!    file upon completion.
!---------------------------------------------------------------------
         allocate (specified_aerosol(kmax_file) )
         read (iounit,FMT = '(5e18.10)')   &
                          (specified_aerosol(k),k=1,kmax_file)
         call close_file (iounit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
       else
         call error_mesg ( 'aerosol_mod', &
              'desired aerosol input file is not present',FATAL)
       endif

!----------------------------------------------------------------------


end subroutine obtain_input_file_data 


!###################################################################### 



                  end module aerosol_mod 



!=======================================================================



#ifdef test_aerosol

program main

use aerosol_mod
use mpp_mod
use mpp_io_mod
use mpp_domains_mod
use time_manager_mod
use diag_manager_mod
use rad_utilities_mod




implicit none

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer, parameter :: NLON=20, NLAT=10,NLEV=8
integer, parameter :: MAX_AERSOL_NAMES = 100
real :: latb(NLON+1,NLAT+1),lonb(NLON+1,NLAT+1),pi,phalf(NLON,NLAT,NLEV+1)
integer :: i,nspecies
type(time_type) :: model_time
character(len=64), dimension(MAX_AEROSOL_NAMES) :: names
type(aerosol_type)  :: Aerosol

pi = 4.*atan(1.)

call mpp_init
call mpp_io_init
call mpp_domains_init
call diag_manager_init
call set_calendar_type(JULIAN)

do i = 1,NLAT+1
   latb(:,i) = -90. + 180.*REAL(i-1)/REAL(NLAT)
end do
do i = 1,NLON+1
   lonb(i,:) = -180. + 360.*REAL(i-1)/REAL(NLON)
end do

latb(:,:) = latb(:,:) * pi/180.
lonb(:,:) = lonb(:,:) * pi/180.

do i = 1,NLEV+1
   phalf(:,:,i) = 101325. * REAL(i-1) / REAL(NLEV)
end do

model_time = set_date(1980,1,1,0,0,0)

call aerosol_init (lonb, latb, names)

call aerosol_driver (1,1,model_time, phalf, Aerosol)

call aerosol_dealloc (Aerosol)

call aerosol_end

call mpp_exit

end program main

#endif
