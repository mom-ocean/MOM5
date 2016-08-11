                       module ozone_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  This code supplies mass mixing ratios of ozone (g/g) to the 
!  sea_esf_rad radiation_package (and the original_fms_rad package).
! </OVERVIEW>
! <DESCRIPTION>
!  This code supplies mass mixing ratios of ozone (g/g) to the 
!  sea_esf_rad radiation_package (and the original_fms_rad package).
!   Recent changes allow provision for predicted ozone to be considered. 
!   This is a field passed in the tracer array, r 
! </DESCRIPTION>
!   shared modules:

use mpp_mod,             only: input_nml_file
use fms_mod,             only: open_namelist_file, file_exist,    &
                               check_nml_error, error_mesg,  &
                               fms_init, stdlog, &
                               write_version_number, FATAL, NOTE, &
                               WARNING, mpp_pe, mpp_root_pe, close_file
use fms_io_mod,          only: read_data
use time_manager_mod,    only: time_type,  &
                               time_manager_init, operator(+), &
                               set_date, operator(-), print_date, &
                               set_time, operator(>), get_date, days_in_month
use diag_manager_mod,    only: diag_manager_init, get_base_time
use time_interp_mod,     only: fraction_of_year, &
                               time_interp_init  
use constants_mod,       only: constants_init, radian
use interpolator_mod,    only: interpolate_type, interpolator_init, &
                               obtain_interpolator_time_slices, &
                               unset_interpolator_time_flag, &
                               interpolator, interpolator_end, &
                               CONSTANT, INTERP_WEIGHTED_P

!-------------------------------

use tracer_manager_mod,  only: get_tracer_index, NO_TRACER
use field_manager_mod,   only: MODEL_ATMOS

!---------------------------------

!   shared radiation package modules:

use rad_utilities_mod,   only:  rad_utilities_init,   &
                                radiative_gases_type,  &
                                atmos_input_type

!---------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!    ozone_mod supplies mass mixing ratios of ozone [ g/g ] to the 
!    model radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: ozone.F90,v 19.0 2012/01/06 20:21:21 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public       &
        ozone_init, ozone_time_vary, ozone_driver, ozone_endts, &
        ozone_end
private         &

!   called from ozone_init:
        obtain_input_file_data, &
        obtain_gfdl_zonal_ozone_data, &
        obtain_clim_zonal_ozone_data, &

!   called from ozone_driver:
        gfdl_zonal_ozone, &
        get_clim_ozone

interface gfdl_zonal_ozone
    module procedure geto3_2d, geto3_3d
end interface


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=24)  ::   basic_ozone_type = 'clim_zonal'     
                      ! label for ozone type, currently unused      
                      ! 'clim_zonal' or 'time_varying' or 'fixed_year'
character(len=32)  :: filename = 'o3.trend.nc'
                      ! name of file which contains the ozone data 
integer, parameter :: MAX_DATA_FIELDS = 1
character(len=32)  :: data_name(MAX_DATA_FIELDS) = 'ozone_1990'
                      ! name of variable in the data file to be used
character(len=24)  ::   ozone_data_source = 'fortuin_kelder' 
                      ! source for the ozone data being used, either 
                      ! 'input', 'gfdl_zonal_ozone', 'calculate_column'
                      ! for the date and location specified,
                      ! or externally-derived datasets:
                      ! 'fortuin_kelder', 'mozart_moztop_fk',
                      ! 'mozart_trop_fk'
character(len=24)  ::   clim_base_year = '1990'
                      ! year with which the ozone data set is assoc-
                      ! iated, used with fortuin_kelder( either '1979',
                      ! '1990' or '1997'), or mozart datasets
                      ! (presently either '1850' or '1860' or '1990')
character(len=24)  ::   trans_data_type = 'linear'
                      ! time interpolation method to be used if trans-
                      ! ient ozone is activated, not yet available
character(len=24)  ::   gfdl_zonal_ozone_type = 'seasonally_varying'
                      ! if gfdl_zonal_ozone is active, the time behavior
                      ! of ozone to use; either 'winter', 'summer', 
                      ! 'spring', 'autumn', annual_mean' or
                      ! 'seasonally_varying'
logical            ::   do_mcm_o3_clim = .false.
                      ! treat ozone as in the manabe climate model ?
integer, dimension(6) ::       &
                        ozone_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
                      ! time in ozone data set corresponding to model
                      ! initial time  (yr, mo, dy, hr, mn, sc)
real,dimension(2)  :: lonb_col = (/-999., -999./)
                      ! longitudes defining the region to use for column
                      ! data calculation
real,dimension(2)  :: latb_col = (/-999., -999./)
                      ! latitudes defining the region to use for column
                      ! data calculation
integer, dimension(6) :: time_col = (/0,0,0,0,0,0/)
                      ! time to use for column calculation
logical      ::  do_coupled_stratozone = .false. ! include the coupled
                                                 ! stratospheric ozone effects?


namelist /ozone_nml/             &
                       basic_ozone_type, &
                       ozone_data_source, &
                       lonb_col, latb_col, time_col, &
                       clim_base_year, &
                       trans_data_type, &
                       data_name, &
                       filename, &
                       gfdl_zonal_ozone_type, &
                       ozone_dataset_entry, &
                       do_mcm_o3_clim, &
                       do_coupled_stratozone

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!-------------------------------------------------------------------
!    qqo3 contains the column input ozone mass mixing ratio (g/g).
!    used when ozone_data_source = 'input'.
!-------------------------------------------------------------------
real, dimension (:), allocatable     ::  qqo3

!-------------------------------------------------------------------
!    o3data contains the zonal ozone mass mixing ratio (10**6 g/g)
!    at 19 specified latitudes, 81 specified pressures and up to 4
!    specified times. rstd contains the zonal ozone mass mixing ratio 
!    (10**6 g/g) at 19 specified latitudes, 81 specified pressures and 
!    at the currently specified time. ph is the array defining the 
!    interface levels in the zonal ozone data set. used when
!    ozone_data_source = 'gfdl_zonal_ozone'.
!-------------------------------------------------------------------
real, dimension (82)        ::    ph
real, dimension (19,81,4)   ::    o3data
real, dimension (19,81)     ::    rstd

!----------------------------------------------------------------------
!    O3 is an interpolate_type variable containing the relevant 
!    information about the ozone data set. used when
!    ozone_data_source = 'fortuin_kelder', 'mozart_moztop_fk',
!    'mozart_trop_fk', and others in the future.
!---------------------------------------------------------------------- 
type(interpolate_type),save         ::  O3_interp 

!---------------------------------------------------------------------
!    miscellaneous variables:

integer     ::  iseason=-1 ! flag indicating type of gfdl_zonal_ozone 
                           ! time variation to use
real        ::  current_fyear=-1.0  
                           ! flag to force interpolation on initial call
                           ! (used with gfdl_zonal_ozone)
logical     ::  do_gfdl_zonal_ozone=.false.
                           ! using gfdl zonal ozone data set ?
logical     ::  do_clim_zonal_ozone=.false.
                           ! using clim zonal ozone data set ?
logical     ::  do_column_input_ozone=.false.
                           ! using ozone column input data set ?
logical     ::  do_predicted_ozone=.false.
                           ! using predicted ozone ?
logical     ::  module_is_initialized=.false.  ! module initialized ?

type(time_type) :: Model_init_time  ! initial calendar time for model  
                                    ! [ time_type ]
type(time_type) :: Ozone_offset     ! difference between model initial
                                    ! time and ozone timeseries app-
                                    ! lied at model initial time
                                    ! [ time_type ]
type(time_type) :: Ozone_entry      ! time in ozone timeseries which
                                    ! is mapped to model initial time
                                    ! [ time_type ]
type(time_type) :: Ozone_time       ! time for which ozone profile is 
                                    ! valid
type(time_type) :: Ozone_column_time
                                    ! time for which ozone data is extr
                                    ! acted from the ozone timeseries
                                    ! in the 'calculate_columns' case
                                    ! [ time_type ]
logical    :: negative_offset = .false.
                            !  the model initial time is later than
                            !  the ozone_dataset_entry time  ?


include 'netcdf.inc'


!-------------------------------------------------------------------
!-------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! <SUBROUTINE NAME="ozone_init">
!  <OVERVIEW>
!   ozone_init is the constructor for ozone_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   ozone_init is the constructor for ozone_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ozone_init (latb, lonb)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   2d array of model longitudes at cell corners [radians]
!  </IN>
! </SUBROUTINE>
!
subroutine ozone_init (latb, lonb)

!---------------------------------------------------------------------
!    ozone_init is the constructor for ozone_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:,:),   intent(in) :: latb, lonb
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      2d array of model latitudes at cell corners 
!                 [ radians ]
!       lonb      2d array of model longitudes at cell corners 
!                 [ radians ]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:

       integer           ::  unit, ierr, io, logunit
       integer           ::  n, no3

!---------------------------------------------------------------------
!  local variables:
!
!         unit       io unit number 
!         ierr       error code
!         io         error status returned from io operation
!
!---------------------------------------------------------------------
 

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!-------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call time_manager_init   
      call diag_manager_init   
      call time_interp_init   
      call constants_init
 
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=ozone_nml, iostat=io)
      ierr = check_nml_error(io,"ozone_nml")
#else
!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=ozone_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'ozone_nml')
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
                        write (logunit, nml=ozone_nml)

!---------------------------------------------------------------------
!    check for a valid value of clim_base_year.
!---------------------------------------------------------------------
      if (trim(basic_ozone_type) == 'clim_zonal' ) then 
        if (trim(clim_base_year)  == '1990' .or.        &
            trim(clim_base_year)  == '1979' .or.        &
            trim(clim_base_year)  == '1850' .or.        &
            trim(clim_base_year)  == '1860' .or.        &
            trim(clim_base_year)  == '1997')            then
        else
          call error_mesg ('ozone_mod', &
              ' clim_base_year must be 1990, 1979, 1850, 1860 &
                                        &or 1997 at present', FATAL)
        endif
        Ozone_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *,   &
           'Ozone data is obtained from a clim_zonal ozone file &
             &for year ', trim(clim_base_year)
        endif

!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the ozone timeseries 
!    files are to be used.                                
!---------------------------------------------------------------------
      else if (trim(basic_ozone_type) == 'time_varying') then
        Model_init_time = get_base_time()
 
!---------------------------------------------------------------------
!    if a dataset entry point is not supplied, use the model base time 
!    as the entry point.
!---------------------------------------------------------------------
        if (ozone_dataset_entry(1) == 1 .and. &
            ozone_dataset_entry(2) == 1 .and. &
            ozone_dataset_entry(3) == 1 .and. &
            ozone_dataset_entry(4) == 0 .and. &
            ozone_dataset_entry(5) == 0 .and. &
            ozone_dataset_entry(6) == 0 ) then
          Ozone_entry = Model_init_time
 
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ozone_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        else
          Ozone_entry  = set_date (ozone_dataset_entry(1), &
                                   ozone_dataset_entry(2), &
                                   ozone_dataset_entry(3), &
                                   ozone_dataset_entry(4), &
                                   ozone_dataset_entry(5), &
                                   ozone_dataset_entry(6))
        endif
        call print_date (Ozone_entry , str='Data from ozone timeseries &
                                           &at time:')
        call print_date (Model_init_time , str='This data is mapped to &
                                                &model time:')
        Ozone_offset = Ozone_entry - Model_init_time
 
        if (Model_init_time > Ozone_entry) then
          negative_offset = .true.
        else
          negative_offset = .false.
        endif
      else if (trim(basic_ozone_type) == 'fixed_year') then
        if (ozone_dataset_entry(1) == 1 .and. &
            ozone_dataset_entry(2) == 1 .and. &
            ozone_dataset_entry(3) == 1 .and. &
            ozone_dataset_entry(4) == 0 .and. &
            ozone_dataset_entry(5) == 0 .and. &
            ozone_dataset_entry(6) == 0 ) then
           call error_mesg ('ozone_mod', &
            'must set ozone_dataset_entry when using  &
                                  &fixed_year ozone', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ozone_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        Ozone_entry  = set_date (ozone_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('ozone_mod', &
           'Ozone data is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'Ozone data obtained from ozone timeseries &
                    &for year:', ozone_dataset_entry(1)
        endif

!------------------------------------------------------------------
!    use predicted ozone
!-----------------------------------------------------------------
      else if (trim(basic_ozone_type) == 'predicted_ozone' ) then
        do_predicted_ozone = .true.
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'Using predicted ozone '
        endif
        no3= get_tracer_index(MODEL_ATMOS,'O3')
        if (no3 == NO_TRACER) &
          call error_mesg ('ozone_mod', &
            'Using predicted ozone but O3 tracer not present in field_table', FATAL)

      else
        call error_mesg ('ozone_mod', &
          'invalid specification of basic_ozone_type', FATAL)
      endif

!------------------------------------------------------------------
!    obtain the desired ozone data set based on the value of namelist
!    variable 'ozone_data_source'. set a logical flag to true to 
!    indicate the desired ozone data source. 
!-----------------------------------------------------------------
      if (trim(ozone_data_source) == 'input' ) then
        do_column_input_ozone = .true.
        call obtain_input_file_data 

      else if (trim(ozone_data_source) == 'gfdl_zonal_ozone' ) then
        do_gfdl_zonal_ozone = .true.
        if (trim(gfdl_zonal_ozone_type)      == 'winter' ) then
          iseason = 1   
        else if (trim(gfdl_zonal_ozone_type) == 'spring' ) then
          iseason = 2   
        else if (trim(gfdl_zonal_ozone_type) == 'summer' ) then
          iseason = 3   
        else if (trim(gfdl_zonal_ozone_type) == 'autumn' ) then
          iseason = 4   
        else if (trim(gfdl_zonal_ozone_type) == 'annual_mean' ) then
          iseason = 0   
        else if (trim(gfdl_zonal_ozone_type) ==    &
                                            'seasonally_varying' ) then
          iseason = 5   
        else
          call error_mesg ( 'ozone_mod', &
                    'improper specification of nml variable  '//&
                                     'gfdl_zonal_ozone_type', FATAL)
        endif
        call obtain_gfdl_zonal_ozone_data (iseason)

      else if (trim(ozone_data_source) == 'fortuin_kelder' ) then
        do_clim_zonal_ozone = .true.
        call obtain_clim_zonal_ozone_data (lonb, latb)

      else if (trim(ozone_data_source) == 'mozart_moztop_fk' ) then
        do_clim_zonal_ozone = .true.
        call obtain_clim_zonal_ozone_data (lonb, latb)

      else if (trim(ozone_data_source) == 'mozart_trop_fk' ) then
        do_clim_zonal_ozone = .true.
        call obtain_clim_zonal_ozone_data (lonb, latb)

      else if (trim(ozone_data_source) == 'calculate_column' ) then
        do_clim_zonal_ozone = .true.
        do n=1,2
          if (lonb_col(n) < 0. .or. lonb_col(n) > 360.) then
            call error_mesg ('ozone_mod', &
                   ' invalid value for lonb_col', FATAL)
          endif
          if (latb_col(n) < -90. .or. latb_col(n) > 90.) then
            call error_mesg ('ozone_mod', &
                ' invalid value for latb_col', FATAL)
          endif
        end do
        if (time_col(1) == 0) then
          call error_mesg ('ozone_mod', &
                'invalid time specified for time_col', FATAL)
        endif
        call interpolator_init (O3_interp, filename,  &
                                spread(lonb_col/RADIAN,2,2),  &
                                spread(latb_col/RADIAN,1,2),&
                                data_out_of_bounds=  (/CONSTANT/), &
                                data_names = data_name, &
                                vert_interp=(/INTERP_WEIGHTED_P/) )
        Ozone_column_time = set_date (time_col(1), time_col(2), &
                                      time_col(3), time_col(4), &
                                      time_col(5), time_col(6))
        call print_date (Ozone_column_time,   &
            str= ' Ozone data used is from ozone timeseries at time: ')
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'Ozone data is averaged over latitudes',  &
                  latb_col(1), ' to', latb_col(2), ' and longitudes',&
                  lonb_col(1), ' to', lonb_col(2)
        endif

      else if (trim(basic_ozone_type) == 'predicted_ozone' ) then
        call error_mesg ('ozone_mod', &
          'Using predicted ozone: no input data set is necessary ', NOTE)

      else
        call error_mesg ( 'ozone_mod',    &
               ' ozone_data_source is not an acceptable value.', FATAL)
      endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.
    


end subroutine ozone_init


!#################################################################### 
 
subroutine ozone_time_vary (model_time)
 
!----------------------------------------------------------------------
!     subroutine ozone_time_vary calculates time-dependent, 
!     space-independent variables needed by this module.
!---------------------------------------------------------------------

type(time_type),    intent(in)   :: model_time
 
!----------------------------------------------------------------------
!
!   local variables
 
      integer         :: yr, mo, dy, hr, mn, sc, dum
      integer         :: dayspmn, mo_yr

      if (do_clim_zonal_ozone) then


       if(trim(basic_ozone_type) == 'time_varying') then
!--------------------------------------------------------------------
!    define the time in the ozone data set from which data is to be 
!    taken. if ozone is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
       if (negative_offset) then
         Ozone_time = model_time - Ozone_offset
       else
         Ozone_time = model_time + Ozone_offset
       endif
     else if(trim(basic_ozone_type) == 'fixed_year') then
       call get_date (Ozone_entry, yr, dum,dum,dum,dum,dum)
       call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
       if (mo ==2 .and. dy == 29) then
         dayspmn = days_in_month(Ozone_entry)
         if (dayspmn /= 29) then
           Ozone_time = set_date (yr, mo, dy-1, hr, mn, sc)
         else
           Ozone_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         Ozone_time = set_date (yr, mo, dy, hr, mn, sc)
       endif
     else if(trim(basic_ozone_type) == 'clim_zonal') then
       Ozone_time = model_time
    endif

!--------------------------------------------------------------------
!    if 'calculate_column' is being used, obtain the ozone values for
!    each column, one at a time, using the pressure profile for that
!    column. this allows each column to see the same ozone fields,
!    but distributed appropriately for its pressure structure.
!--------------------------------------------------------------------
   if (trim(ozone_data_source) == 'calculate_column') then
     call obtain_interpolator_time_slices (O3_interp, Ozone_column_time)
   else
     call obtain_interpolator_time_slices (O3_interp, Ozone_time)
   endif

!----------------------------------------------------------------------

     endif
 
end subroutine ozone_time_vary



!######################################################################

subroutine ozone_endts


     call unset_interpolator_time_flag (O3_interp)


end subroutine ozone_endts




!######################################################################
! <SUBROUTINE NAME="ozone_driver">
!  <OVERVIEW>
!    ozone_driver obtains the current ozone distributions and returns 
!   them in Rad_gases%qo3.
!  </OVERVIEW>
!  <DESCRIPTION>
!   ozone_driver obtains the current ozone distributions and returns 
!   them in Rad_gases%qo3.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ozone_driver (is, ie, js, je, lat, Rad_time, Atmos_input, &
!                      Rad_gases )
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   latitude of model points  [ radians ]
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which the climatologically-determined,
!                   time-varying ozone field should apply
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing the atmospheric
!   input fields needed by the radiation package
!  </IN> 
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable which will return
!                   the ozone mass mixing ratio (g/g) to the calling
!                   routine
!  </INOUT>
! </SUBROUTINE>
!

subroutine ozone_driver (is, ie, js, je, lat, Rad_time, Atmos_input, &
                         r, Rad_gases )
 
!--------------------------------------------------------------------
!   ozone_driver obtains the current ozone distributions and returns 
!   them in Rad_gases%qo3.
!--------------------------------------------------------------------

integer,                    intent(in)    :: is, ie, js, je
real, dimension(:,:),       intent(in)    :: lat
type(time_type),            intent(in)    :: Rad_time
type(atmos_input_type),     intent(in)    :: Atmos_input
real, dimension(:,:,:,:),   intent(in)    :: r
type(radiative_gases_type), intent(inout) :: Rad_gases

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      Rad_time     time at which the climatologically-determined,
!                   time-varying ozone field should apply
!                   [ time_type (days, seconds) ] 
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!
!   intent(out) variables:
!
!      Rad_gases    radiative_gases_type variable which will return
!                   the ozone mass mixing ratio [ g/g ] to the calling
!                   routine
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!      phaf         pressure at model interface levels, normalized
!                   by the mean sea level pressure (101325 N/m**2).
!                   if mcm ozone scheme is used, this variable is
!                   defined differently. [ Pa ]
!      kmax         number of model layers
!      k            do-loop index
!
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,1), &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3))    :: phaf

      integer   ::   kmax
      integer   ::   k, j, i   !  do-loop indices
      integer   ::   noy

!--------------------------------------------------------------------
      if ( .not. module_is_initialized)    & 
        call error_mesg ('ozone_mod',  &         
           'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kmax = size (Rad_gases%qo3,3) 

!-----------------------------------------------------------------------
!    if column input ozone is being used, be certain the input column
!    has the same number of levels as the model grid. broadcast the
!    column data over the entire horizontal domain.
!-----------------------------------------------------------------------
      if (do_column_input_ozone) then 
        if (size(qqo3(:)) /= kmax) then
           call error_mesg ('ozone_mod', &
             'size of ozone profile in input file does not match '//&
                                              'model grid.', FATAL)
         endif
         do k=1,kmax
           Rad_gases%qo3(:,:,k) = qqo3(k)
         end do
!---------------------------------------------------------------------
!    if a specified ozone column is not being used, define the 
!    normalized pressure to be used to assign the ozone values. this
!    formulation varies between mcm and the fms models.
!---------------------------------------------------------------------
      else 
        do k=1,kmax+1
          if (do_mcm_o3_clim) then
            phaf(:,:,k) = 100000.*Atmos_input%phalf(:,:,k)/  &
                          Atmos_input%phalf(:,:,kmax+1)
          else
            phaf(:,:,k) = (Atmos_input%pflux(:,:,k))
!           phaf(:,:,k) = (Atmos_input%pflux(:,:,k))*101325./   &
!                         (Atmos_input%pflux(:,:,kmax +1))
          endif
        end do

!---------------------------------------------------------------------
!    if gfdl_zonal_ozone has been activated, call gfdl_zonal_ozone
!    to obtain the values of ozone at model grid points at the desired
!    time Rad_time.
!---------------------------------------------------------------------
        if (do_gfdl_zonal_ozone) then
          call gfdl_zonal_ozone (Rad_time, lat, phaf, Rad_gases%qo3)

!---------------------------------------------------------------------
!    if clim_zonal_ozone has been activated, call clim_zonal_ozone
!    to obtain the values of ozone at model grid points at the desired
!    time Rad_time.
!---------------------------------------------------------------------
        else if (do_clim_zonal_ozone) then
          call get_clim_ozone (is, js, Rad_time, phaf, Rad_gases%qo3)



!---------------------------------------------------------------------
!    if do_predicted_ozone has been activated, set ozone to the appropriate
!    prognostic tracer field;  here defined by tracer name 'O3'
!---------------------------------------------------------------------
        else if ( do_predicted_ozone ) then 

          noy= get_tracer_index(MODEL_ATMOS,'O3')

          do k= 1,kmax
            do j= 1,je-js+1    
              do i= 1,ie-is+1    
!                Rad_gases%qo3(i,j,k) =  r(i,j,k,noy) * 48./29.0 
                 Rad_gases%qo3(i,j,k) =  MAX(r(i,j,k,noy),1.e-20) * 48./29.
             end do
            end do
          end do

        endif
      endif

!--------------------------------------------------------------------




end subroutine ozone_driver



!#####################################################################
! <SUBROUTINE NAME="ozone_end">
!  <OVERVIEW>
!   ozone_end is the destructor for ozone_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   ozone_end is the destructor for ozone_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ozone_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine ozone_end

!---------------------------------------------------------------------
!    ozone_end is the destructor for ozone_mod.
!---------------------------------------------------------------------
     
!---------------------------------------------------------------------
!    deallocate any module variables that have been allocated.
!---------------------------------------------------------------------
      if (do_column_input_ozone) then
        deallocate (qqo3)
      endif

!---------------------------------------------------------------------
!    call interpolator_end to close out the O3_interp interpolate_type
!    variable.
!---------------------------------------------------------------------
      if (do_clim_zonal_ozone) then
        call interpolator_end (O3_interp)
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.



end subroutine ozone_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!######################################################################
! <SUBROUTINE NAME="obtain_input_file_data">
!  <OVERVIEW>
!   obtain_input_file_data reads an input file containing a single
!    column ozone profile. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_input_file_data reads an input file containing a single
!    column ozone profile. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  obtain_input_file_data
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine obtain_input_file_data 

!---------------------------------------------------------------------
!    obtain_input_file_data reads an input file containing a single
!    column ozone profile.
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
      character(len=31) ::  dummy


!-------------------------------------------------------------------
!    determine if a netcdf input data file exists. if so, read the
!    number of data records in the file.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/id1o3.nc') ) then
        ncid = ncopn ('INPUT/id1o3.nc', 0, rcode)
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
        allocate (qqo3(kmax_file) )
        ivarid = ncvid(ncid, 'ozone', rcode)
        call ncvinq (ncid, ivarid, dummy, ntp, nvdim, vdims, nvs, rcode)
        do j=1,nvdim
          call ncdinq (ncid, vdims(j), dummy, ndsize, rcode)
          start(j) = 1
          count(j) = ndsize
        end do
        call ncvgt (ncid, ivarid, start, count, qqo3, rcode)
        call ncclos (ncid, rcode)

!-------------------------------------------------------------------
!    determine if the input data input file exists in ascii format.
!    if so, read the number of data records in the file.
!---------------------------------------------------------------------
      else if (file_exist ( 'INPUT/id1o3') ) then
        iounit = open_namelist_file ('INPUT/id1o3')
        read (iounit,FMT = '(i4)') kmax_file

!-------------------------------------------------------------------
!    allocate space for the input data. read the data set. close the 
!    file upon completion.
!---------------------------------------------------------------------
         allocate (qqo3(kmax_file) )
         read (iounit,FMT = '(5e18.10)') (qqo3(k),k=1,kmax_file)
         call close_file (iounit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
       else
         call error_mesg ( 'ozone_mod', &
              'desired ozone input file is not present',FATAL)
       endif

!----------------------------------------------------------------------


end subroutine obtain_input_file_data 




!######################################################################
! <SUBROUTINE NAME="obtain_gfdl_zonal_ozone_data">
!  <OVERVIEW>
!   obtain_gfdl_zonal_ozone_data generates data at the appropriate time
!    from the basic fms_zonal_ozone input data set, allowing the use of
!    annual mean, fixed seasonal, or seasonally-varying ozone distrib-
!    utions.
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_gfdl_zonal_ozone_data generates data at the appropriate time
!    from the basic fms_zonal_ozone input data set, allowing the use of
!    annual mean, fixed seasonal, or seasonally-varying ozone distrib-
!    utions.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  obtain_gfdl_zonal_ozone_data (season)
!  </TEMPLATE>
!  <IN NAME="season" TYPE="integer">
!   scalar integer between 0-5, where 1-4 uses fixed
!                      data (1=winter, 2=spring, etc.), season=0 is 
!                      annual mean ozone, and season=5 is seasonally
!                      varying ozone
!  </IN> 
! </SUBROUTINE>
!
subroutine obtain_gfdl_zonal_ozone_data (season)

!---------------------------------------------------------------------
!    obtain_gfdl_zonal_ozone_data retrieves ozone data as requested
!    from the gfdl_zonal_ozone input data set. data corresponding to
!    annual mean, fixed seasonal, or seasonally-varying ozone distrib-
!    utions may be obtained.
!---------------------------------------------------------------------

integer, intent(in) :: season

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      season          scalar integer between 0-5, where 1-4 uses fixed
!                      data (1=winter, 2=spring, etc.), season=0 is 
!                      annual mean ozone, and season=5 is seasonally
!                      varying ozone
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:

      real, dimension (10,41)             ::  ro31, ro32
      real, dimension (19,41)             ::  duo3n
      real, dimension (19,81,4)           ::  data4
      real, dimension (82)                ::  ph3
      real, dimension (10,25)             ::  o3hi
      real, dimension (10,16)             :: o3lo1, o3lo2, o3lo3, o3lo4

      real                      ::  pref = 101325. 
      character(len=48)         ::  err_string 
      integer                   ::  iounit   
      integer                   ::  j, k, kk, n  

!-----------------------------------------------------------------------
!  local variables:
! 
!     ro31, ro32       ozone values for the proper time, at 41 standard
!                      levels and 10 standard latitudes (index 1 = 
!                      equator, index 10= pole, 10 deg resolution)
!     duo3n            array of ozone at proper time, at 41 standard 
!                      levels, over 19 global latitudes (10 deg 
!                      resolution), index 1 = north pole, index 10=
!                      equator, index 19 = south pole)
!     data4            array of ozone at proper time, at 81 levels, over
!                      19 global latitudes (10 deg resolution). if 
!                      seasonally-varying ozone is desired, there are
!                      4 such arrays, to which harmonic interpolation
!                      will be applied.
!     ph3              sigma levels at which zonal ozone data set data
!                      exists
!     o3hi             array containing ozone data at 25 high layers in
!                      which mixing ratios are invariant with season
!     o3lo1            array containing ozone data at 16 lower layers
!                      valid for nh spring, latitudinal indices from
!                      equator to pole
!     o3lo2            array containing ozone data at 16 lower layers
!                      valid for nh fall, latitudinal indices from
!                      equator to pole
!     o3lo3            array containing ozone data at 16 lower layers
!                      valid for nh winter, latitudinal indices from
!                      equator to pole
!     o3lo4            array containing ozone data at 16 lower layers
!                      valid for nh summer, latitudinal indices from
!                      equator to pole
!     pref             assumed surface pressure value used to convert 
!                      ph3 from sigma to pressure
!     err_string       part of error message if generated
!     iounit           unit number used to read input data file
!     j,k,kk,n         do-loop indices
!     
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure that the input argument is valid.
!---------------------------------------------------------------------
      if (season < 0 .or. season > 5) then
        write (err_string,9001) season
 9001   format ('invalid value of season=',i10)
        call error_mesg ('ozone_mod', err_string, FATAL)
      endif

!---------------------------------------------------------------------
!    save input argument as a module variable.
!---------------------------------------------------------------------
      iseason = season

!---------------------------------------------------------------------
!    determine if the zonal ozone input data file exists. if so, 
!    open the file and read the data set.  close the file upon 
!    completion. if it is not present, write an error message.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/zonal_ozone_data.nc')) then
        if(mpp_pe() == mpp_root_pe()) &
             call error_mesg('ozone_mod','Reading netCDF input data: zonal_ozone_data.nc',NOTE)
        call read_data('INPUT/zonal_ozone_data.nc', 'ph3', ph3, no_domain=.true.)
        call read_data('INPUT/zonal_ozone_data.nc', 'o3hi', o3hi, no_domain=.true.)
        call read_data('INPUT/zonal_ozone_data.nc', 'o3lo1', o3lo1, no_domain=.true.)
        call read_data('INPUT/zonal_ozone_data.nc', 'o3lo2', o3lo2, no_domain=.true.)
        call read_data('INPUT/zonal_ozone_data.nc', 'o3lo3', o3lo3, no_domain=.true.)
        call read_data('INPUT/zonal_ozone_data.nc', 'o3lo4', o3lo4, no_domain=.true.)
      else if (file_exist ( 'INPUT/zonal_ozone_data') ) then
        call error_mesg('ozone_mod','Reading native input data zonal_ozone_data no longer supported',FATAL)
      else
        call error_mesg ( 'ozone_mod', &
                'zonal_ozone_data data file is not present',FATAL)
      endif

!---------------------------------------------------------------------
!    define standard pressure interfaces (ph).
!---------------------------------------------------------------------
      ph(:) = ph3(:)*pref

!-----------------------------------------------------------------------
!    define the seasonally-invariant elements of arrays ro31, ro32 from
!    the values in o3hi.
!---------------------------------------------------------------------
      do k=1,25
         ro31(:,k) = o3hi(:,k)
         ro32(:,k) = o3hi(:,k)
      end do

!---------------------------------------------------------------------
!    define the ozone values at the lower levels to be seasonally
!    varying. northern hemisphere seasons are used to define the
!    indices, with n = 1 being winter, going to n = 4 being fall.
!--------------------------------------------------------------------
      do n=1,4

!---------------------------------------------------------------------
!    for nh spring or fall, obtain the o3lo1 and o3lo2 data (ro31 and 
!    ro32). 
!----------------------------------------------------------------------
        if (n == 2 .or. n == 4) then
          do k=26,41
            ro31(:,k) = o3lo1(:,k-25)
            ro32(:,k) = o3lo2(:,k-25)
          end do
        endif

!---------------------------------------------------------------------
!    for nh winter or summer, obtain the o3lo3 and o3lo4 data (ro31 and
!    ro32). 
!----------------------------------------------------------------------
        if (n == 1 .or. n == 3) then
          do k=26,41
            ro31(:,k) = o3lo3(:,k-25)
            ro32(:,k) = o3lo4(:,k-25)
          end do
        endif

!---------------------------------------------------------------------
!    define ozone values for both hemispheres -- indices 1->9 = nh,
!    index 10 = equator, indices 11-19 = sh. for nh spring (n=2)
!    and nh winter (n=1), nh values are contained in ro31, in
!    reverse latitudinal order. sh values are contained in ro32, going
!    from equator to south pole.
!----------------------------------------------------------------------
        if (n == 2 .or. n == 1) then
          do k=1,41
            do j=1,10
              duo3n(j  ,k) = ro31(11-j,k)
              duo3n(j+9,k) = ro32(j   ,k)
            end do
            duo3n(10 ,k) = 0.50*(ro31(1,k) + ro32(1,k))
          end do

!---------------------------------------------------------------------
!    for nh summer (n=3) and nh fall (n=4), nh values are 
!    contained in ro32, in reverse latitudinal order. sh values are 
!    contained in ro31, going from equator to south pole.
!---------------------------------------------------------------------
        else if(n == 4 .or. n == 3) then
          do k=1,41
            do j=1,10
              duo3n(j  ,k) = ro32(11-j,k)
              duo3n(j+9,k) = ro31(j   ,k)
            end do
            duo3n(10 ,k) = 0.50*(ro31(1,k) + ro32(1,k))
          end do
        endif

!-----------------------------------------------------------------------
!    vertical interp between original data using bessels half-point 
!    interpolation formula.
!-----------------------------------------------------------------------
        do kk=4,78,2
          k = kk/2
          o3data(:,kk,n) = 0.50*(duo3n(:,k) + duo3n(:,k+1)) -  &
                           (duo3n(:,k+2) - duo3n(:,k+1) -   &
                            duo3n(:,k) + duo3n(:,k-1))/16.
        end do
        o3data(:, 2,n) = 0.50*(duo3n(:,2) + duo3n(:,1))
        o3data(:,80,n) = 0.50*(duo3n(:,41) + duo3n(:,40))

!---------------------------------------------------------------------
!    put intermediate (unchanged) data into new array.                
!---------------------------------------------------------------------
        do kk=1,81,2
          k = (kk + 1)/2
          o3data(:,kk,n) = duo3n(:,k)
        end do
      end do  ! n loop

!-----------------------------------------------------------------------
!    prepare seasonal interpolation.
!-----------------------------------------------------------------------
      if (iseason == 5) then
        data4(:,:,1) = 0.25*(o3data(:,:,1) + o3data(:,:,2)  &
                           + o3data(:,:,3) + o3data(:,:,4))
        data4(:,:,2) = 0.50*(o3data(:,:,2) - o3data(:,:,4))
        data4(:,:,3) = 0.50*(o3data(:,:,1) - o3data(:,:,3))
        data4(:,:,4) = 0.25*(o3data(:,:,1) - o3data(:,:,2)  &
                           + o3data(:,:,3) - o3data(:,:,4))
        o3data(:,:,1) = data4(:,:,1)
        o3data(:,:,2) = data4(:,:,2)
        o3data(:,:,3) = data4(:,:,3)
        o3data(:,:,4) = data4(:,:,4)
      endif

!-----------------------------------------------------------------------
!    prepare annual mean data. store it into o3data with index 1. reset
!    the o3data index to be 1 so that this overwritten field will be
!    retrieved.
!-----------------------------------------------------------------------
      if (iseason == 0) then
        data4(:,:,1) = 0.25*(o3data(:,:,1) + o3data(:,:,2)  &
                           + o3data(:,:,3) + o3data(:,:,4))
        o3data(:,:,1) = data4(:,:,1)
        iseason = 1
      endif

!---------------------------------------------------------------------



end subroutine obtain_gfdl_zonal_ozone_data




!######################################################################
! <SUBROUTINE NAME="obtain_clim_zonal_ozone_data">
!  <OVERVIEW>
!   obtain_clim_zonal_ozone_data provides the necessary information 
!    to interpolator_mod so that the appropriate clim_ozone data may
!    be obtained later on when needed.
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_clim_zonal_ozone_data provides the necessary information 
!    to interpolator_mod so that the appropriate clim_ozone data may
!    be obtained later on when needed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  obtain_clim_zonal_ozone_data (lonb, latb) 
!  </TEMPLATE>
!  <IN NAME="lonb, latb" TYPE="real">
!       lonb      2d array of model longitudes at cell corners [radians]
!       latb      2d array of model latitudes at cell corners [radians]
!  </IN> 
! </SUBROUTINE>
!
subroutine obtain_clim_zonal_ozone_data (lonb, latb)

!----------------------------------------------------------------------
!    obtain_clim_zonal_ozone_data calls interpolator_init to supply
!    the necessary information to interpolator_mod to allow the appro-
!    priate clim_zonal_ozone data to be obtained later when needed.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in) :: lonb, latb

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       lonb      2d array of model longitudes at cell corners [radians]
!       latb      2d array of model latitudes at cell corners [radians]
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------
!    call interpolator_init to initialize an interp_type variable
!    O3_interp which will be used to retrieve interpolated ozone
!    data when requested.
!---------------------------------------------------------------------
        call interpolator_init (O3_interp, filename, lonb, &
                                latb, data_out_of_bounds=(/CONSTANT/), &
                                data_names = data_name, &
                                vert_interp=(/INTERP_WEIGHTED_P/) )


!----------------------------------------------------------------------


end subroutine obtain_clim_zonal_ozone_data



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    INTERFACE GFDL_ZONAL_OZONE
!
!
! call gfdl_zonal_ozone (Rad_time, lat, phaf, ozone)
!
!  separate routines exist within this interface for 3d and 2d
!  array input and output:
!
!  real, dimension(:,:),    intent(in)  :: lat
!  real, dimension(:,:,:),  intent(in)  :: phalf
!  real, dimension(:,:,:),  intent(out) :: ozone
! OR
!  real, dimension(:),    intent(in)    :: lat
!  real, dimension(:,:),  intent(in)    :: phalf
!  real, dimension(:,:),  intent(out)   :: ozone
!
!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!      Time         current model time [ time_type (days, seconds) ] 
!      lat          latitude of model points  [ radians ]
!      phalf        pressure at model layer interfaces [ Pa ]
!
!  intent(out) variables:
!
!      ozone        ozone mass mixing ratio at model levels [ g / g ]
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#######################################################################
! <SUBROUTINE NAME="geto3_3d">
!  <OVERVIEW>
!   geto3_3d retrieves an (i,j,k) array of ozone mass mixing ratio 
!    (g / g ) valid at the specified time to be returned to the 
!    calling routine.
!  </OVERVIEW>
!  <DESCRIPTION>
!   geto3_3d retrieves an (i,j,k) array of ozone mass mixing ratio 
!    (g / g ) valid at the specified time to be returned to the 
!    calling routine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call geto3_3d (Time, lat, phalf, ozone)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current model time [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   latitude of model points  [ radians ]
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model layer interfaces [ kg / (m s^2) ]
!  </IN>
!  <OUT NAME="ozone" TYPE="real">
!   ozone mass mixing ratio at model levels [ g / g ]
!  </OUT>  
! </SUBROUTINE>
!
subroutine geto3_3d (Time, lat, phalf, ozone)

!---------------------------------------------------------------------
!    geto3_3d retrieves an (i,j,k) array of ozone mass mixing ratio 
!    [ g / g ] valid at the specified time to be returned to the 
!    calling routine.
!---------------------------------------------------------------------

type(time_type),         intent(in)  :: Time
real, dimension(:,:),    intent(in)  :: lat
real, dimension(:,:,:),  intent(in)  :: phalf
real, dimension(:,:,:),  intent(out) :: ozone

!---------------------------------------------------------------------
!  local variables:
!
!      rlag         time lag of  valid time of seasonal ozone data from 
!                   start of calendar year [ years ]
!      profile      ozone profile in model grid column
!      dp           pressure difference between data layer interfaces
!      dp1, dp2, dp3, dp4          
!                   pressure differences used to assign data set layer
!                   ozone to the proper model layer
!      o3col        total ozone in a model pressure layer
!      rang         angular position of specified time from start of
!                   ozone year
!      rsin1        argument for harmonic interpolation
!      rcos1        argument for harmonic interpolation
!      rcos2        argument for harmonic interpolation
!      phd          model latitude, guaranteed to be between -90 and +90
!                   degrees
!      th           relative latitudinal distance of model grid point
!                   from nearest lower data set index, normalized by
!                   data set resolution 
!      fyear        fraction of the year which has passed at the 
!                   specified time
!      j1           nearest lower latitude index in data set to the
!                   model latitude
!      j2           j1 + 1; data from index j1 and j2 will be inter-
!                   polated to the model grid point
!
!------------------------------------------------------------------
      real                :: rlag
      real, dimension(81) :: profile, dp, dp1, dp2, dp3, dp4, o3col
      real                :: rang, rsin1, rcos1, rcos2, phd, th, fyear
      integer             :: j1,j2

      integer             :: i,j,k,l  ! various indices

!--------------------------------------------------------------------
!    when seasonally-varying ozone is desired, perform a harmonic time 
!    interpolation to obtain values at the specified time.
!--------------------------------------------------------------------
      if (iseason == 5) then
        if (do_mcm_o3_clim) then
          rlag = 12.6875/365.
        else
          rlag = 1./24.
        endif
        fyear = fraction_of_year (time)
        if (fyear /= current_fyear) then
          rang = 4.0*acos(0.0)*(fyear-rlag)
          rsin1 = sin(rang)
          rcos1 = cos(rang)
          rcos2 = cos(2.0*rang)
          rstd(:,:) = o3data(:,:,1) + rsin1*o3data(:,:,2) +  &
                      rcos1*o3data(:,:,3) + rcos2*o3data(:,:,4)
          current_fyear = fyear
        endif
!---------------------------------------------------------------------
!    otherwise, no interpolation is needed, use the specified seasonal 
!    data. if annual mean value has been specified, that data will be
!    at iseason = 1.
!---------------------------------------------------------------------
      else
        rstd(:,:) = o3data(:,:,iseason)
      endif

!---------------------------------------------------------------------
!    define the pressure increments of the standard data levels.
!---------------------------------------------------------------------
      do l=1,81
         dp (l) = ph(l+1) - ph(l)
      end do

!---------------------------------------------------------------------
!    perform horizontal interpolation into the data set. profile is
!    the vertical ozone data at the grid point.
!---------------------------------------------------------------------
      do j=1,size(lat,2)
        do i=1,size(lat,1)
          phd = max(min(lat(i,j)*radian, 90.), -90.)
          j1 = 10.000 - phd*0.10
          j1 = max(min(j1, 18), 1)
          j2 = j1 + 1
          th = (10-j1) - phd*0.10
          profile(:) = rstd(j1,:) + th*(rstd(j2,:) - rstd(j1,:))

!---------------------------------------------------------------------
!    now interpolate in the vertical to produce values at model
!    levels. calculate layer-mean ozone mixing ratio from data set for 
!    each model layer.
!---------------------------------------------------------------------
          do k=1,size(ozone,3)

!---------------------------------------------------------------------
!    calculate sums over data set layers to get model layer 
!    ozone. the mcm ozone data is obtained in a somewhat
!    different manner.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    mcm algorithm.
!--------------------------------------------------------------------
            if (do_mcm_o3_clim) then
              ozone(i,j,k) = 0.0
              do l=1,81
                if ((ph(l+1) >= phalf(i,j,k)) .and.   &
                    (ph(l) <= phalf(i,j,k+1))) then
                  if ((ph(l+1) < phalf(i,j,k+1)) .and.    &
                      (ph(l) < phalf(i,j,k)))   then
                    ozone(i,j,k) = ozone(i,j,k) + profile(l)*   &
                                   (ph(l+1) - phalf(i,j,k))
                  endif
                  if ((ph(l+1) < phalf(i,j,k+1)) .and.    &
                      (ph(l) >= phalf(i,j,k)))   then  
                    ozone(i,j,k) = ozone(i,j,k) + profile(l)*    &
                                   (ph(l+1) - ph(l))
                  endif
                  if ((ph(l+1) > phalf(i,j,k+1)) .and.    &
                      (ph(l) > phalf(i,j,k)))  then 
                    ozone(i,j,k) = ozone(i,j,k) + profile(l)*   &
                                   (phalf(i,j,k+1) - ph(l))
                  endif
                endif
              end do

!---------------------------------------------------------------------
!    fms algorithm.
!---------------------------------------------------------------------
            else 
              do l=1,81
                dp1(l) = ph(l+1) - phalf(i,j,k)
                dp2(l) = phalf(i,j,k+1) - ph(l)
                dp3(l) = ph(l+1) - phalf(i,j,k+1)
                dp4(l) = ph(l) - phalf(i,j,k)
              end do
              where (dp1(:) < 0.0) dp1(:)=0.0
              where (dp2(:) < 0.0) dp2(:)=0.0
              do l=1,81
                o3col(l) = 0.0
                if ( dp3(l) < 0.0 ) then
                  if ( dp4(l) < 0.0 ) then
                    o3col(l) = profile(l)*dp1(l)
                  else
                    o3col(l) = profile(l)*dp(l)
                  endif
                else
                  if ( dp4(l) < 0.0 ) then
                    o3col(l) = profile(l)*    &
                               (phalf(i,j,k+1) - phalf(i,j,k))
                  else
                    o3col(l) = profile(l)*dp2(l)
                  endif
                endif
              end do
              ozone(i,j,k) = sum(o3col(:))
            endif ! do_mcm_o3_clim if block

!---------------------------------------------------------------------
!    normalize by the model pressure depth to produce a mass 
!    mixing ratio.
!---------------------------------------------------------------------
            ozone(i,j,k) = ozone(i,j,k)/(phalf(i,j,k+1)-phalf(i,j,k))

!---------------------------------------------------------------------
!    code to cover case when surface pressure is greater than pref.
!---------------------------------------------------------------------
            if (.not.do_mcm_o3_clim .and. ph(82) < phalf(i,j,k+1)) then
              ozone(i,j,k) = profile(81)
            endif

!----------------------------------------------------------------------
!    code to cover case when model resolution is so fine that no value
!    of ph(l) in the ozone data array falls betwen phalf(k+1) and
!    phalf(k).   procedure is to simply grab the nearest value from
!    rdata (or profile in the fms code).
!----------------------------------------------------------------------
            if (do_mcm_o3_clim) then
              if (ozone(i,j,k) <= 0.0) then
                do l=1,81
                  if (ph(l) < phalf(i,j,k) .and.     &
                      ph(l+1) > phalf(i,j,k+1) ) then
                    ozone(i,j,k) = profile(l)
                  endif
                end do
              endif
            endif 
          end do ! k loop
        end do   ! i loop
      end do     ! j loop

!---------------------------------------------------------------------
!    convert units from micrograms/gram to kg/kg.
!---------------------------------------------------------------------
      ozone(:,:,:) = ozone(:,:,:)*1.e-6

!-----------------------------------------------------------------------



end subroutine geto3_3d




!#######################################################################
! <SUBROUTINE NAME="geto3_2d">
!  <OVERVIEW>
!   geto3_2d retrieves an (i,k) array of ozone mass mixing ratio 
!    (g / g ) valid at the specified time to be returned to the 
!    calling routine.
!  </OVERVIEW>
!  <DESCRIPTION>
!   geto3_2d retrieves an (i,k) array of ozone mass mixing ratio 
!    (g / g ) valid at the specified time to be returned to the 
!    calling routine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  geto3_2d(Time, lat, phalf, ozone)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current model time [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   latitude of model points  [ radians ]
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model layer interfaces [ kg / (m s^2) ]
!  </IN>
!  <OUT NAME="ozone" TYPE="real">
!   ozone mass mixing ratio at model levels [ g / g ]
!  </OUT>  
! </SUBROUTINE>
!
subroutine geto3_2d (Time, lat, phalf, ozone)

!---------------------------------------------------------------------
!    geto3_2d retrieves an (i,k) array of ozone mass mixing ratio 
!    [ g / g ] valid at the specified time to be returned to the 
!    calling routine.
!---------------------------------------------------------------------

type(time_type),       intent(in)  :: Time
real, dimension(:),    intent(in)  :: lat
real, dimension(:,:),  intent(in)  :: phalf
real, dimension(:,:),  intent(out) :: ozone

!---------------------------------------------------------------------
!  local variables:
!
!      lat3         2d equivalent of lat
!      phalf3       3d equivalent of phalf
!      ozone3       3d equivalent of ozone
!
!--------------------------------------------------------------------

      real,dimension (size(lat,1),1)                  :: lat3
      real,dimension (size(phalf,1),1, size(phalf,2)) :: phalf3
      real,dimension (size(ozone,1),1, size(ozone,2)) :: ozone3

!---------------------------------------------------------------------
!    add an extra dummy dimension to the input variables.
!---------------------------------------------------------------------
      lat3(:,1)     = lat(:)
      phalf3(:,1,:) = phalf(:,:)

!---------------------------------------------------------------------
!    call the 3d interface of this procedure.
!---------------------------------------------------------------------
      call geto3_3d (time, lat3, phalf3, ozone3)

!---------------------------------------------------------------------
!    remove the extra dummy dimension from the output variables.
!---------------------------------------------------------------------
      ozone(:,:) = ozone3(:,1,:)

!--------------------------------------------------------------------


end subroutine geto3_2d


 
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                END INTERFACE GFDL_ZONAL_OZONE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!#######################################################################
! <SUBROUTINE NAME="get_clim_ozone">
!  <OVERVIEW>
!   get_clim_ozone retrieves the clim_ozone field at the desired place 
!    and time from the o3.climatology.nc file by accessing 
!    interpolator_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   get_clim_ozone retrieves the clim_ozone field at the desired place 
!    and time from the o3.climatology.nc file by accessing 
!    interpolator_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  get_clim_ozone (is, js, model_time, p_half, model_data)
!  </TEMPLATE>
!  <IN NAME="model_time" TYPE="time_type">
!   time at which the climatologically-determined,
!                   time-varying ozone field should apply
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model layer interfaces [ kg / (m s^2) ]
!  </IN>
!  <OUT NAME="model_data" TYPE="real">
!   output field containing ozone field at desired time
!  </OUT>  
!  <IN NAME="is, js" TYPE="integer">
!   OPTIONAL: starting subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
! </SUBROUTINE>
!
subroutine get_clim_ozone (is, js, model_time, p_half, model_data)

!--------------------------------------------------------------------
!    get_clim_ozone retrieves the clim_ozone field for the requested
!    points at the desired time from the clim ozone data file by
!    accessing interpolator_mod.
!---------------------------------------------------------------------

integer,                intent(in)      :: is,js
type(time_type),        intent(in)      :: model_time
real, dimension(:,:,:), intent(in)      :: p_half
real, dimension(:,:,:), intent(out)     :: model_data

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js        starting subdomain i,j indices of data in 
!                   the physics_window being integrated
!      model_time   time at which the climatologically-determined,
!                   time-varying ozone field should apply
!                   [ time_type (days, seconds) ] 
!      p_half       pressure at model interface levels
!                   [ Pa ]
!
!   intent(out) variables:
!
!      model_data   output field containing ozone field at desired time
!                   [ g / g ]
!
!
!----------------------------------------------------------------------
!
!   local variables
 
      real, dimension(1,1, size(p_half,3)-1) :: ozone_data
      real, dimension(1,1, size(p_half,3)) :: p_half_col
      integer         :: i, j

!--------------------------------------------------------------------
!    if 'calculate_column' is being used, obtain the ozone values for
!    each column, one at a time, using the pressure profile for that
!    column. this allows each column to see the same ozone fields,
!    but distributed appropriately for its pressure structure.
!--------------------------------------------------------------------
     if (trim(ozone_data_source) == 'calculate_column') then
       do j=1, size(p_half,2)
         do i=1, size(p_half,1)
           p_half_col(1,1,:) = p_half(i,j,:)
           call interpolator (O3_interp, Ozone_column_time,  &
                              p_half_col, ozone_data, &
!                             trim(data_name(1)), is, js)
                              trim(data_name(1)), 1, 1  )
           model_data(i,j,:) = ozone_data(1,1,:)
         end do
       end do
     else

!--------------------------------------------------------------------
!    call interpolator to obtain data at the specified grid points and 
!    specified time.
!--------------------------------------------------------------------
      
      call interpolator (O3_interp, Ozone_time, p_half, model_data,  &
                         trim(data_name(1)), is, js)
     endif 

!----------------------------------------------------------------------


end subroutine get_clim_ozone



!#####################################################################



                       end module ozone_mod



