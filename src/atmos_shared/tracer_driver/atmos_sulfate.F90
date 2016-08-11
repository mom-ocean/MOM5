module atmos_sulfate_mod
! <DESCRIPTION>
!   This module is an implementation of sulfate chemistry. It contains
!   tracer emissions and chemistry. The chemistry is partly based on MOZART.
!   The change of concentration of SO2, DMS, SO4, MSA and H2O2 are
!   calculated using monthly mean concentration of OH, HO2, jH2O2, NO3, O3,
!   pH. The emissions include:
!     - DMS from seawater
!     - SO2 by fossil fuel, biomass burning, non-eruptive volcanoes and aircraft
!     - SO4 by fossil fuel
! </DESCRIPTION>
! <WARNING>
!  To save space only the actual month of input files are kept in memory.
!  This implies that the "atmos_sulfate_init" should be executed at the begining
!  of each month. In other words, the script should not run more than 1 month
!  without a restart.
! </WARNING>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file 
use                    fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog,                  &
                                       check_nml_error, error_mesg, &
                                       open_namelist_file, FATAL, NOTE, WARNING

use           time_manager_mod, only : time_type, &
                                       days_in_month, days_in_year, &
                                       set_date, set_time, get_date_julian, &
                                       print_date, get_date, &
                                       operator(>), operator(+), operator(-)
use time_interp_mod,            only:  fraction_of_year, &
                                       time_interp_init
use           diag_manager_mod, only : send_data,               &
                                       register_diag_field,     &
                                       register_static_field,   &
                                       diag_manager_init, get_base_time
use         tracer_manager_mod, only : get_tracer_index,        &
                                       set_tracer_atts
use          field_manager_mod, only : MODEL_ATMOS
use           interpolator_mod, only:  interpolate_type, interpolator_init, &
                                      obtain_interpolator_time_slices, &
                                      unset_interpolator_time_flag, &
                                       interpolator, interpolator_end,     &
                                       CONSTANT, INTERP_WEIGHTED_P
use              constants_mod, only : PI, GRAV, RDGAS, WTMAIR

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  atmos_sulfate_init, atmos_sulfate_end, &
        atmos_sulfate_time_vary, atmos_sulfate_endts, &
        atmos_DMS_emission, atmos_SOx_emission, atmos_SOx_chem

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: nSO4 = 0  ! tracer number for Sulfate               = SO4=
integer :: nDMS = 0  ! tracer number for Dimethyl sulfide      = CH3SCH3
integer :: nSO2 = 0  ! tracer number for Sulfur dioxide        = SO2
integer :: nMSA = 0  ! tracer number for Methane sulfonic acid = CH3SO3H
integer :: nH2O2= 0  ! tracer number for Hydrogen peroxyde     = H2O2

real , parameter :: WTM_S     = 32.0
real , parameter :: WTM_O3    = 48.0
real , parameter :: WTM_SO2   = 64.0
real , parameter :: WTM_SO4   = 96.0
real , parameter :: WTM_NH4_2SO4   = 132.00
real , parameter :: WTM_DMS   = 62.0
real , parameter :: WTM_MSA   = 96.0

!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH                  = 0
integer ::   id_HO2                 = 0
integer ::   id_NO3                 = 0
integer ::   id_jH2O2               = 0
integer ::   id_O3                  = 0
integer ::   id_pH                  = 0

integer ::   id_DMSo                = 0
integer ::   id_DMS_emis            = 0
integer ::   id_DMS_emis_cmip       = 0
integer ::   id_SO2_emis            = 0
integer ::   id_SO4_emis            = 0
integer ::   id_DMS_chem            = 0
integer ::   id_SO2_chem            = 0
integer ::   id_SO4_chem            = 0
integer ::   id_SO4_oh_prod         = 0
integer ::   id_SO4_o3_prod         = 0
integer ::   id_SO4_h2o2_prod       = 0
integer ::   id_MSA_chem            = 0
integer ::   id_H2O2_chem           = 0
integer ::   id_so2_aircraft        = 0
integer ::   id_so2_cont_volc       = 0
integer ::   id_so2_expl_volc       = 0
integer ::   id_so2_biobur          = 0
integer ::   id_so2_ship            = 0
integer ::   id_so2_road            = 0
integer ::   id_so2_domestic        = 0
integer ::   id_so2_industry        = 0
integer ::   id_so2_power           = 0
integer ::   id_so2_off_road        = 0
integer ::   id_so2_ff              = 0

type(interpolate_type),save         ::  gas_conc_interp
type(interpolate_type),save         ::  aerocom_emission_interp
type(interpolate_type),save         ::  gocart_emission_interp
type(interpolate_type),save         ::  anthro_emission_interp
type(interpolate_type),save         ::  biobur_emission_interp
type(interpolate_type),save         ::  ship_emission_interp
type(interpolate_type),save         ::  aircraft_emission_interp
! type(interpolate_type),save         ::  cont_volc_emission_interp
! type(interpolate_type),save         ::  expl_volc_emission_interp
! Initial calendar time for model
type(time_type) :: model_init_time
type(time_type), save :: gas_conc_offset
type(time_type), save :: anthro_offset
type(time_type), save :: biobur_offset
type(time_type), save :: ship_offset
type(time_type), save :: aircraft_offset
! type(time_type), save :: cont_volc_offset
! type(time_type), save :: expl_volc_offset

type(time_type), save :: gas_conc_entry
type(time_type), save :: anthro_entry
type(time_type), save :: biobur_entry
type(time_type), save :: ship_entry
type(time_type), save :: aircraft_entry
! type(time_type), save :: cont_volc_entry
! type(time_type), save :: expl_volc_entry

logical, save    :: gas_conc_negative_offset
logical, save    :: anthro_negative_offset
logical, save    :: biobur_negative_offset
logical, save    :: ship_negative_offset
logical, save    :: aircraft_negative_offset
! logical, save    :: cont_volc_negative_offset
! logical, save    :: expl_volc_negative_offset

integer, save    :: gas_conc_time_serie_type
integer, save    :: anthro_time_serie_type
integer, save    :: biobur_time_serie_type
integer, save    :: ship_time_serie_type
integer, save    :: aircraft_time_serie_type
! integer, save    :: cont_volc_time_serie_type
! integer, save    :: expl_volc_time_serie_type

real             :: critical_sea_fraction = 0.5 ! DMS flux from sea occurs
                                ! in grid cells with ocn_flx_fraction .gt.
                                !  this value
                             
character(len=80)  :: runtype = 'default'

character(len=80)  :: gocart_emission_filename = 'gocart_emission.nc'
character(len=80), dimension(6) :: gocart_emission_name
data gocart_emission_name/'DMSo','SO2_GEIA1','SO2_GEIA2', &
                       'SO4_GEIA1','SO4_GEIA2','SO2_biobur'/

integer, parameter :: num_volc_levels = 12
character(len=80)  :: cont_volc_source = ' '
character(len=80)  :: expl_volc_source = ' '
real :: volc_altitude_edges(num_volc_levels+1) = 1.e3 * (/ &
  0.,0.1,0.2,0.5,1.,2.,3.,4.,5.,6.,7.,8.,20. /) ! m

character(len=80)  :: aerocom_emission_filename = 'aerocom_emission.nc'
integer, parameter :: std_aerocom_emission=18, &
                      max_aerocom_emission=std_aerocom_emission+2*num_volc_levels
character(len=80), dimension(max_aerocom_emission)  :: aerocom_emission_name = (/ &
         'SO2_RoadTransport         ', 'SO2_Off-road              ', &
         'SO2_Domestic              ', 'SO2_Industry              ', &
         'SO2_International_Shipping', 'SO2_Powerplants           ', &
         'SO2_cont_volc             ', 'alt_cont_volc_low         ', &
         'alt_cont_volc_high        ', 'SO2_expl_volc             ', &
         'alt_expl_volc_low         ', 'alt_expl_volc_high        ', &
         'GFED_SO2_l1               ', 'GFED_SO2_l2               ', &
         'GFED_SO2_l3               ', 'GFED_SO2_l4               ', &
         'GFED_SO2_l5               ', 'GFED_SO2_l6               ', &
         'SO2_cont_volc_l01         ', 'SO2_cont_volc_l02         ', &
         'SO2_cont_volc_l03         ', 'SO2_cont_volc_l04         ', &
         'SO2_cont_volc_l05         ', 'SO2_cont_volc_l06         ', &
         'SO2_cont_volc_l07         ', 'SO2_cont_volc_l08         ', &
         'SO2_cont_volc_l09         ', 'SO2_cont_volc_l10         ', &
         'SO2_cont_volc_l11         ', 'SO2_cont_volc_l12         ', &
         'SO2_expl_volc_l01         ', 'SO2_expl_volc_l02         ', &
         'SO2_expl_volc_l03         ', 'SO2_expl_volc_l04         ', &
         'SO2_expl_volc_l05         ', 'SO2_expl_volc_l06         ', &
         'SO2_expl_volc_l07         ', 'SO2_expl_volc_l08         ', &
         'SO2_expl_volc_l09         ', 'SO2_expl_volc_l10         ', &
         'SO2_expl_volc_l11         ', 'SO2_expl_volc_l12         '/)

character(len=80)  :: gas_conc_source   = ' '
character(len=80)  :: gas_conc_filename = 'gas_conc_3D.nc'
character(len=80), dimension(6) :: gas_conc_name
data gas_conc_name/'OH','HO2','NO3','O3','jH2O2','pH'/
character(len=80)     :: gas_conc_time_dependency_type
integer, dimension(6) :: gas_conc_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)

character(len=80)  :: anthro_source   = ' '
character(len=80)  :: anthro_filename = 'aero_anthro_emission_1979_2006.nc'
character(len=80), dimension(2) :: anthro_emission_name
data anthro_emission_name/'so2_anthro','so4_anthro'/
character(len=80)     :: anthro_time_dependency_type
integer, dimension(6) :: anthro_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)

character(len=80)  :: biobur_source   = ' '
character(len=80)  :: biobur_filename = 'aero_biobur_emission_1979_2006.nc'
character(len=80), dimension(2) :: biobur_emission_name
data biobur_emission_name/'so2_biobur','so4_biobur'/
character(len=80)     :: biobur_time_dependency_type
integer, dimension(6) :: biobur_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)

character(len=80)  :: ship_source   = ' '
character(len=80)  :: ship_filename = 'aero_ship_emission_1979_2006.nc'
character(len=80), dimension(2) :: ship_emission_name
data ship_emission_name/'so2_ship','so4_ship'/
character(len=80)     :: ship_time_dependency_type
integer, dimension(6) :: ship_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)

character(len=80)  :: aircraft_source = ' '
character(len=80)  :: aircraft_filename = 'aircraft_emission.nc'
character(len=80)  :: aircraft_emission_name(1)
data aircraft_emission_name/'fuel'/
character(len=80)     :: aircraft_time_dependency_type
integer, dimension(6) :: aircraft_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
real :: so2_aircraft_EI = 1.e-3  ! kg of SO2/kg of fuel

namelist /simple_sulfate_nml/  &
       critical_sea_fraction,     &
      runtype,                         &
      aerocom_emission_filename, aerocom_emission_name,  &
      gocart_emission_filename, gocart_emission_name,  &
      gas_conc_source, gas_conc_name, gas_conc_filename,        &
        gas_conc_time_dependency_type, gas_conc_dataset_entry, &
      anthro_source, anthro_emission_name, anthro_filename,        &
        anthro_time_dependency_type, anthro_dataset_entry, &
      biobur_source, biobur_emission_name, biobur_filename,        &
        biobur_time_dependency_type, biobur_dataset_entry, &
      ship_source, ship_emission_name, ship_filename,        &
        ship_time_dependency_type, ship_dataset_entry, &
      aircraft_source, aircraft_emission_name, aircraft_filename, &
        aircraft_time_dependency_type, aircraft_dataset_entry, so2_aircraft_EI,&
      cont_volc_source, expl_volc_source

type(time_type) :: anthro_time, biobur_time, ship_time, aircraft_time
type(time_type)        :: gas_conc_time

!trim(runtype) 
!biomass_only; fossil_fuels_only, natural_only, anthrop

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_sulfate.F90,v 20.0 2013/12/13 23:24:05 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="atmos_sulfate_init">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the sulfate module.
!</DESCRIPTION>
 subroutine atmos_sulfate_init ( lonb, latb, nlev, axes, Time, mask)

!-----------------------------------------------------------------------
real, intent(in),    dimension(:,:)                 :: lonb, latb
integer, intent(in)                                 :: nlev
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
character(len=7), parameter :: mod_name = 'tracers'
integer :: n, m, nsulfate
!
!----------------------------------------------------------------------
!  local variables:

      integer   ::   ierr

!---------------------------------------------------------------------
!    local variables:
!
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!         n          do-loop index
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!
      integer  unit,io, logunit
      character(len=12) :: SOx_tracer(5)
!
!     1. DMS       = Dimethyl sulfide            = CH3SCH3
!     2. SO2       = Sulfur dioxide              = SO2     
!     3. SO4       = Sulfate                     = SO4=   
!     4. MSA       = Methane sulfonic acid       = CH3SO3H
!     5. H2O2      = Hydrogen peroxyde           = H2O2
!                                                                      
      data SOx_tracer/'simpleDMS', &
                      'simpleSO2', &
                      'simpleSO4', &
                      'simpleMSA', &
                      'simpleH2O2' /
      
      if (module_is_initialized) return


!---- write namelist ------------------

      call write_version_number(version, tagname)

!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=simple_sulfate_nml, iostat=io)
        ierr = check_nml_error(io,'simple_sulfate_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=simple_sulfate_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'simple_sulfate_nml')
        end do
10      call close_file (unit)
#endif
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=simple_sulfate_nml)


!----- set initial value of sulfate ------------

     do m=1,size(SOx_tracer)

       n = get_tracer_index(MODEL_ATMOS,SOx_tracer(m))
       if (n>0) then
         nsulfate=n
        call set_tracer_atts(MODEL_ATMOS,SOx_tracer(m),SOx_tracer(m),'vmr')
         if (nsulfate > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (logunit,30) SOx_tracer(m),nsulfate
       endif
     enddo


  30   format (A,' was initialized as tracer number ',i2)

!----------------------------------------------------------------------
!    initialize namelist entries
!----------------------------------------------------------------------
        gas_conc_offset = set_time (0,0)
        anthro_offset   = set_time (0,0)
        biobur_offset   = set_time (0,0)
        ship_offset     = set_time (0,0)
        aircraft_offset = set_time (0,0)

        gas_conc_entry  = set_time (0,0)
        anthro_entry    = set_time (0,0)
        biobur_entry    = set_time (0,0)
        ship_entry      = set_time (0,0)
        aircraft_entry  = set_time (0,0)

        gas_conc_negative_offset = .false.
        anthro_negative_offset   = .false.
        biobur_negative_offset   = .false.
        ship_negative_offset     = .false.
        aircraft_negative_offset = .false.

        gas_conc_time_serie_type = 1
        anthro_time_serie_type   = 1
        biobur_time_serie_type   = 1
        ship_time_serie_type     = 1
        aircraft_time_serie_type = 1
!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table)
!----------------------------------------------------------------------
        model_init_time = get_base_time()
       call interpolator_init (aerocom_emission_interp, &
                             trim(aerocom_emission_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = aerocom_emission_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )
       call interpolator_init (gocart_emission_interp, &
                             trim(gocart_emission_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = gocart_emission_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(gas_conc_time_dependency_type) == 'constant' ) then
        gas_conc_time_serie_type = 1
        gas_conc_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'gas_conc are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for gas_conc is selected.
!---------------------------------------------------------------------
      else if (trim(gas_conc_time_dependency_type) == 'time_varying') then
        gas_conc_time_serie_type = 3
        if (gas_conc_dataset_entry(1) == 1 .and. &
            gas_conc_dataset_entry(2) == 1 .and. &
            gas_conc_dataset_entry(3) == 1 .and. &
            gas_conc_dataset_entry(4) == 0 .and. &
            gas_conc_dataset_entry(5) == 0 .and. &
            gas_conc_dataset_entry(6) == 0 ) then
          gas_conc_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to gas_conc_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          gas_conc_entry  = set_date (gas_conc_dataset_entry(1), &
                                  gas_conc_dataset_entry(2), &
                                  gas_conc_dataset_entry(3), &
                                  gas_conc_dataset_entry(4), &
                                  gas_conc_dataset_entry(5), &
                                  gas_conc_dataset_entry(6))
        endif
        call print_date (gas_conc_entry , str= &
          'Data from gas_conc timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        gas_conc_offset = gas_conc_entry - model_init_time
        if (model_init_time > gas_conc_entry) then
          gas_conc_negative_offset = .true.
        else
          gas_conc_negative_offset = .false.
        endif
      else if (trim(gas_conc_time_dependency_type) == 'fixed_year') then
        gas_conc_time_serie_type = 2
        if (gas_conc_dataset_entry(1) == 1 .and. &
            gas_conc_dataset_entry(2) == 1 .and. &
            gas_conc_dataset_entry(3) == 1 .and. &
            gas_conc_dataset_entry(4) == 0 .and. &
            gas_conc_dataset_entry(5) == 0 .and. &
            gas_conc_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_sulfate_mod', &
            'must set gas_conc_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to gas_conc_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        gas_conc_entry  = set_date (gas_conc_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_sulfate_mod', &
           'gas_conc is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'gas_conc correspond to year :', &
                    gas_conc_dataset_entry(1)
        endif
     endif
     call interpolator_init (gas_conc_interp, trim(gas_conc_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = gas_conc_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )
      if (trim(anthro_source) .eq. 'do_anthro') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(anthro_time_dependency_type) == 'constant' ) then
        anthro_time_serie_type = 1
        anthro_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'anthro are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for anthro is selected.
!---------------------------------------------------------------------
      else if (trim(anthro_time_dependency_type) == 'time_varying') then
        anthro_time_serie_type = 3
        if (anthro_dataset_entry(1) == 1 .and. &
            anthro_dataset_entry(2) == 1 .and. &
            anthro_dataset_entry(3) == 1 .and. &
            anthro_dataset_entry(4) == 0 .and. &
            anthro_dataset_entry(5) == 0 .and. &
            anthro_dataset_entry(6) == 0 ) then
          anthro_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to anthro_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          anthro_entry  = set_date (anthro_dataset_entry(1), &
                                  anthro_dataset_entry(2), &
                                  anthro_dataset_entry(3), &
                                  anthro_dataset_entry(4), &
                                  anthro_dataset_entry(5), &
                                  anthro_dataset_entry(6))
        endif
        call print_date (anthro_entry , str= &
          'Data from anthro timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        anthro_offset = anthro_entry - model_init_time
        if (model_init_time > anthro_entry) then
          anthro_negative_offset = .true.
        else
          anthro_negative_offset = .false.
        endif
      else if (trim(anthro_time_dependency_type) == 'fixed_year') then
        anthro_time_serie_type = 2
        if (anthro_dataset_entry(1) == 1 .and. &
            anthro_dataset_entry(2) == 1 .and. &
            anthro_dataset_entry(3) == 1 .and. &
            anthro_dataset_entry(4) == 0 .and. &
            anthro_dataset_entry(5) == 0 .and. &
            anthro_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_sulfate_mod', &
            'must set anthro_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to anthro_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        anthro_entry  = set_date (anthro_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_sulfate_mod', &
           'anthro is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'anthro correspond to year :', &
                    anthro_dataset_entry(1)
        endif
     endif
     call interpolator_init (anthro_emission_interp,             &
                             trim(anthro_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = anthro_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
     endif ! end do_anthro

    if (trim(biobur_source) .eq. 'do_biobur') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(biobur_time_dependency_type) == 'constant' ) then
        biobur_time_serie_type = 1
        biobur_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'biobur are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for biobur is selected.
!---------------------------------------------------------------------
      else if (trim(biobur_time_dependency_type) == 'time_varying') then
        biobur_time_serie_type = 3
        if (biobur_dataset_entry(1) == 1 .and. &
            biobur_dataset_entry(2) == 1 .and. &
            biobur_dataset_entry(3) == 1 .and. &
            biobur_dataset_entry(4) == 0 .and. &
            biobur_dataset_entry(5) == 0 .and. &
            biobur_dataset_entry(6) == 0 ) then
          biobur_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to biobur_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          biobur_entry  = set_date (biobur_dataset_entry(1), &
                                  biobur_dataset_entry(2), &
                                  biobur_dataset_entry(3), &
                                  biobur_dataset_entry(4), &
                                  biobur_dataset_entry(5), &
                                  biobur_dataset_entry(6))
        endif
        call print_date (biobur_entry , str= &
          'Data from biobur timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        biobur_offset = biobur_entry - model_init_time
        if (model_init_time > biobur_entry) then
          biobur_negative_offset = .true.
        else
          biobur_negative_offset = .false.
        endif
      else if (trim(biobur_time_dependency_type) == 'fixed_year') then
        biobur_time_serie_type = 2
        if (biobur_dataset_entry(1) == 1 .and. &
            biobur_dataset_entry(2) == 1 .and. &
            biobur_dataset_entry(3) == 1 .and. &
            biobur_dataset_entry(4) == 0 .and. &
            biobur_dataset_entry(5) == 0 .and. &
            biobur_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_sulfate_mod', &
            'must set biobur_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to biobur_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        biobur_entry  = set_date (biobur_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_sulfate_mod', &
           'biobur is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'biobur correspond to year :', &
                    biobur_dataset_entry(1)
        endif
     endif
     call interpolator_init (biobur_emission_interp,             &
                             trim(biobur_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = biobur_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
     endif

     if (trim(ship_source) .eq. 'do_ship') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(ship_time_dependency_type) == 'constant' ) then
        ship_time_serie_type = 1
        ship_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ship are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for ship is selected.
!---------------------------------------------------------------------
      else if (trim(ship_time_dependency_type) == 'time_varying') then
        ship_time_serie_type = 3
        if (ship_dataset_entry(1) == 1 .and. &
            ship_dataset_entry(2) == 1 .and. &
            ship_dataset_entry(3) == 1 .and. &
            ship_dataset_entry(4) == 0 .and. &
            ship_dataset_entry(5) == 0 .and. &
            ship_dataset_entry(6) == 0 ) then
          ship_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ship_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          ship_entry  = set_date (ship_dataset_entry(1), &
                                  ship_dataset_entry(2), &
                                  ship_dataset_entry(3), &
                                  ship_dataset_entry(4), &
                                  ship_dataset_entry(5), &
                                  ship_dataset_entry(6))
        endif
        call print_date (ship_entry , str= &
          'Data from ship timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        ship_offset = ship_entry - model_init_time
        if (model_init_time > ship_entry) then
          ship_negative_offset = .true.
        else
          ship_negative_offset = .false.
        endif
      else if (trim(ship_time_dependency_type) == 'fixed_year') then
        ship_time_serie_type = 2
        if (ship_dataset_entry(1) == 1 .and. &
            ship_dataset_entry(2) == 1 .and. &
            ship_dataset_entry(3) == 1 .and. &
            ship_dataset_entry(4) == 0 .and. &
            ship_dataset_entry(5) == 0 .and. &
            ship_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_sulfate_mod', &
            'must set ship_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ship_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        ship_entry  = set_date (ship_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_sulfate_mod', &
           'ship is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ship correspond to year :', &
                    ship_dataset_entry(1)
        endif
     endif
     call interpolator_init (ship_emission_interp,             &
                             trim(ship_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = ship_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
    endif

    if (trim(aircraft_source) .eq. 'do_aircraft') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(aircraft_time_dependency_type) == 'constant' ) then
        aircraft_time_serie_type = 1
        aircraft_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'aircraft are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for aircraft is selected.
!---------------------------------------------------------------------
      else if (trim(aircraft_time_dependency_type) == 'time_varying') then
        aircraft_time_serie_type = 3
        if (aircraft_dataset_entry(1) == 1 .and. &
            aircraft_dataset_entry(2) == 1 .and. &
            aircraft_dataset_entry(3) == 1 .and. &
            aircraft_dataset_entry(4) == 0 .and. &
            aircraft_dataset_entry(5) == 0 .and. &
            aircraft_dataset_entry(6) == 0 ) then
          aircraft_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to aircraft_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          aircraft_entry  = set_date (aircraft_dataset_entry(1), &
                                  aircraft_dataset_entry(2), &
                                  aircraft_dataset_entry(3), &
                                  aircraft_dataset_entry(4), &
                                  aircraft_dataset_entry(5), &
                                  aircraft_dataset_entry(6))
        endif
        call print_date (aircraft_entry , str= &
          'Data from aircraft timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        aircraft_offset = aircraft_entry - model_init_time
        if (model_init_time > aircraft_entry) then
          aircraft_negative_offset = .true.
        else
          aircraft_negative_offset = .false.
        endif
      else if (trim(aircraft_time_dependency_type) == 'fixed_year') then
        aircraft_time_serie_type = 2
        if (aircraft_dataset_entry(1) == 1 .and. &
            aircraft_dataset_entry(2) == 1 .and. &
            aircraft_dataset_entry(3) == 1 .and. &
            aircraft_dataset_entry(4) == 0 .and. &
            aircraft_dataset_entry(5) == 0 .and. &
            aircraft_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_sulfate_mod', &
            'must set aircraft_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to aircraft_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        aircraft_entry  = set_date (aircraft_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_sulfate_mod', &
           'aircraft is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'aircraft correspond to year :', &
                    aircraft_dataset_entry(1)
        endif
     endif

     call interpolator_init (aircraft_emission_interp, &
                             trim(aircraft_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = aircraft_emission_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )
   endif

! Register diagnostic fields
   id_DMS_emis   = register_diag_field ( mod_name,                           &
                   'simpleDMS_emis', axes(1:2),Time,                         &
                   'simpleDMS_emis', 'kgS/m2/s',                             &
                    missing_value=-999.  )
   id_DMS_emis_cmip   = register_diag_field ( mod_name,                                   &
                   'simpleDMS_emis_cmip', axes(1:2),Time,                                 &
                   'simpleDMS_emis_cmip', 'kgDMS/m2/s',                                     &
                    missing_value=-999.  )
   id_SO2_emis   = register_diag_field ( mod_name,                           &
                   'simpleSO2_emis', axes(1:3),Time,                         &
                   'simpleSO2_emis', 'kgS/m2/s',                             &
                    missing_value=-999.  )
   id_SO4_emis   = register_diag_field ( mod_name,                           &
                   'simpleSO4_emis', axes(1:3),Time,                         &
                   'simpleSO4_emis', 'kgS/m2/s',                             &
                    missing_value=-999.  )
   id_DMSo       = register_diag_field ( mod_name,                           &
                   'DMSo',axes(1:2),Time,                                    &
                   'Dimethylsulfide seawater concentration',                 &
                   'nM/L')
   id_ph          = register_diag_field ( mod_name,                          &
                   'pH_simple_sulfate',axes(1:3),Time,                       &
                   'pH in simple-sulfate',                                   &
                   'none')
   id_O3           = register_diag_field ( mod_name,                         &
                   'O3_simple_sulfate',axes(1:3),Time,                       &
                   'O3 in simple-sulfate',                                   &
                   'none')
   id_SO2_aircraft = register_diag_field ( mod_name,                         &
                   'simpleSO2_aircraft_emis',axes(1:3),Time,                 &
                   'simpleSO2 emission by aircraft',                         &
                   'kgS/m2/s')
   id_SO2_biobur  = register_diag_field ( mod_name,                          &
                   'simpleSO2_biobur_emis',axes(1:3),Time,                   &
                   'simpleSO2 emission from biomass burning',                &
                   'kgS/m2/s')
   id_SO2_cont_volc = register_diag_field ( mod_name,                        &
                   'simpleSO2_cont_volc_emis',axes(1:3),Time,                &
                   'simpleSO2 emission from non-eruptive volcanoes',         &
                   'kgS/m2/s')
   id_SO2_expl_volc = register_diag_field ( mod_name,                        &
                   'simpleSO2_expl_volc_emis',axes(1:3),Time,                &
                   'simpleSO2 emission from eruptive volcanoes',             &
                   'kgS/m2/s')
   id_SO2_ship      = register_diag_field ( mod_name,                        &
                   'simpleSO2_ship_emis',axes(1:3),Time,                     &
                   'simpleSO2 emission from international shipping',         &
                   'kgS/m2/s')
   id_SO2_road      = register_diag_field ( mod_name,                        &
                   'simpleSO2_road_emis',axes(1:3),Time,                     &
                   'simpleSO2 emission from road transport',                 &
                   'kgS/m2/s')
   id_SO2_domestic = register_diag_field ( mod_name,                         &
                   'simpleSO2_domestic_emis',axes(1:3),Time,                 &
                   'simpleSO2 emission from domestic fossil fuel burning',   &
                   'kgS/m2/s')
   id_SO2_industry = register_diag_field ( mod_name,                         &
                   'simpleSO2_industry_emis',axes(1:3),Time,                 &
                   'simpleSO2 emission from industrial fossil fuel burning', &
                   'kgS/m2/s')
   id_SO2_power   = register_diag_field ( mod_name,                          &
                   'simpleSO2_power_emis',axes(1:3),Time,                    &
                   'simpleSO2 emission from power plants',                   &
                   'kgS/m2/s')
   id_SO2_off_road = register_diag_field ( mod_name,                         &
                   'simpleSO2_off_road_emis',axes(1:3),Time,                 &
                   'simpleSO2 emission from off-road transport',             &
                   'kgS/m2/s')
   id_SO2_ff      = register_diag_field ( mod_name,                          &
                   'simpleSO2_ff_emis',axes(1:3),Time,                       &
                   'simpleSO2 emission from fossil fuel burning',            &
                   'kgS/m2/s')
   id_NO3        = register_diag_field ( mod_name,                           &
                   'simpleNO3_diurnal',axes(1:3),Time,                       &
                   'Time varying NO3 concentration',                         &
                   'molec.cm-3')
   id_OH         = register_diag_field ( mod_name,                           &
                   'OH_simple_sulfate',axes(1:3),Time,                       &
                   'Varying Hydroxyl radical concentration',                 &
                   'molec.cm-3')
   id_jH2O2         = register_diag_field ( mod_name,                        &
                   'jH2O2_simple_sulfate',axes(1:3),Time,                    &
                   'Varying H2O2 photodissociation',                         &
                   's-1')
   id_HO2         = register_diag_field ( mod_name,                          &
                   'HO2_simple_sulfate',axes(1:3),Time,                      &
                   'Varying Hydroperoxyl radical concentration',             &
                   'molec.cm-3')
   id_DMS_chem   = register_diag_field ( mod_name,                           &
                   'simpleDMS_chem',axes(1:3),Time,                          &
                   'simpleDMS chemical production',                          &
                   'kgS/m2/s')
   id_SO2_chem   = register_diag_field ( mod_name,                           &
                   'simpleSO2_chem',axes(1:3),Time,                          &
                   'simpleSO2 chemical production',                          &
                   'kgS/m2/s')
   id_SO4_chem   = register_diag_field ( mod_name,                           &
                   'simpleSO4_chem',axes(1:3),Time,                          &
                   'simpleSO4 chemical production',                          &
                   'kgS/m2/s')
   id_SO4_oh_prod= register_diag_field ( mod_name,                           &
                   'simpleSO4_oh_prod',axes(1:3),Time,                       &
                   'simpleSO4 gas phase production',                         &
                   'kgS/m2/s')
   id_SO4_o3_prod= register_diag_field ( mod_name,                           &
                   'simpleSO4_o3_prod',axes(1:3),Time,                       &
                   'simpleSO4 aqueous phase production SO2+O3',              &
                   'kgS/m2/s')
   id_SO4_h2o2_prod= register_diag_field ( mod_name,                         &
                   'simpleSO4_h2o2_prod',axes(1:3),Time,                     &
                   'simpleSO4 aqueous phase production SO2+H2O2',            &
                   'kgS/m2/s')
   id_MSA_chem   = register_diag_field ( mod_name,                           &
                   'simpleMSA_chem',axes(1:3),Time,                          &
                   'simpleMSA chemical production',                          &
                   'kgS/m2/s')
   id_H2O2_chem   = register_diag_field ( mod_name,                          &
                   'simpleH2O2_chem',axes(1:3),Time,                         &
                   'simpleH2O2 chemical production',                         &
                   'kgH2O2/m2/s')

   call write_version_number(version, tagname)

   module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
 end subroutine atmos_sulfate_init



!######################################################################

subroutine atmos_sulfate_time_vary (model_time)


type(time_type), intent(in) :: model_time

      integer :: yr,dum, mo_yr, mo, dy, hr, mn, sc, dayspmn


      call obtain_interpolator_time_slices (gocart_emission_interp, &
                                                           model_time)

      call obtain_interpolator_time_slices (aerocom_emission_interp, &
                                                           model_time)

      if (trim(anthro_source) .eq. 'do_anthro') then
!--------------------------------------------------------------------
!    define the time in the anthro data set from which data is to be 
!    taken. if anthro is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
          if(anthro_time_serie_type .eq. 3) then
            if (anthro_negative_offset) then
              anthro_time = model_time - anthro_offset
            else
              anthro_time = model_time + anthro_offset
            endif
          else
            if(anthro_time_serie_type .eq. 2 ) then
              call get_date (anthro_entry, yr, dum,dum,dum,dum,dum)
              call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
              if (mo ==2 .and. dy == 29) then
                dayspmn = days_in_month(anthro_entry)
                if (dayspmn /= 29) then
                  anthro_time = set_date (yr, mo, dy-1, hr, mn, sc)
                else
                  anthro_time = set_date (yr, mo, dy, hr, mn, sc)
                endif
              else
                anthro_time = set_date (yr, mo, dy, hr, mn, sc)
              endif
            else
              anthro_time = model_time
            endif
          endif
          call obtain_interpolator_time_slices &
                      (anthro_emission_interp, anthro_time)

      endif 

      if (trim(biobur_source) .eq. 'do_biobur') then
!--------------------------------------------------------------------
!    define the time in the biobur data set from which data is to be 
!    taken. if biobur is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
          if(biobur_time_serie_type .eq. 3) then
            if (biobur_negative_offset) then
              biobur_time = model_time - biobur_offset
            else
              biobur_time = model_time + biobur_offset
            endif
          else
            if(biobur_time_serie_type .eq. 2 ) then
              call get_date (biobur_entry, yr, dum,dum,dum,dum,dum)
              call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
              if (mo ==2 .and. dy == 29) then
                dayspmn = days_in_month(biobur_entry)
                if (dayspmn /= 29) then
                  biobur_time = set_date (yr, mo, dy-1, hr, mn, sc)
                else
                  biobur_time = set_date (yr, mo, dy, hr, mn, sc)
                endif
              else
                biobur_time = set_date (yr, mo, dy, hr, mn, sc)
              endif
            else
              biobur_time = model_time
            endif
          endif
          call obtain_interpolator_time_slices &
                            (biobur_emission_interp, biobur_time)

      endif 

      if (trim(ship_source) .eq. 'do_ship') then
!--------------------------------------------------------------------
!    define the time in the ship data set from which data is to be 
!    taken. if ship is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
          if(ship_time_serie_type .eq. 3) then
            if (ship_negative_offset) then
              ship_time = model_time - ship_offset
            else
              ship_time = model_time + ship_offset
            endif
          else
            if(ship_time_serie_type .eq. 2 ) then
              call get_date (ship_entry, yr, dum,dum,dum,dum,dum)
              call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
              if (mo ==2 .and. dy == 29) then
                dayspmn = days_in_month(ship_entry)
                if (dayspmn /= 29) then
                  ship_time = set_date (yr, mo, dy-1, hr, mn, sc)
                else
                  ship_time = set_date (yr, mo, dy, hr, mn, sc)
                endif
              else
                ship_time = set_date (yr, mo, dy, hr, mn, sc)
              endif
            else
              ship_time = model_time
            endif
          endif
          call obtain_interpolator_time_slices &
              (ship_emission_interp, ship_time)

        endif
!
! Aircraft emissions
      if (trim(aircraft_source) .eq. 'do_aircraft') then
        if(aircraft_time_serie_type .eq. 3) then
          if (aircraft_negative_offset) then
            aircraft_time = model_time - aircraft_offset
          else
            aircraft_time = model_time + aircraft_offset
          endif
        else
          if(aircraft_time_serie_type .eq. 2 ) then
            call get_date (aircraft_entry, yr, dum,dum,dum,dum,dum)
            call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
            if (mo ==2 .and. dy == 29) then
              dayspmn = days_in_month(aircraft_entry)
              if (dayspmn /= 29) then
                aircraft_time = set_date (yr, mo, dy-1, hr, mn, sc)
              else
                aircraft_time = set_date (yr, mo, dy, hr, mn, sc)
              endif
            else
              aircraft_time = set_date (yr, mo, dy, hr, mn, sc)
            endif
          else
            aircraft_time = model_time
          endif
        endif

!
          call obtain_interpolator_time_slices &
                         (aircraft_emission_interp, aircraft_time)
      endif


!--------------------------------------------------------------------
!    define the time in the gas_conc data set from which data is to be 
!    taken. if gas_conc is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(gas_conc_time_serie_type .eq. 3) then
       if (gas_conc_negative_offset) then
         gas_conc_time = model_time - gas_conc_offset
       else
         gas_conc_time = model_time + gas_conc_offset
       endif
     else
       if(gas_conc_time_serie_type .eq. 2 ) then
         call get_date (gas_conc_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(gas_conc_entry)
           if (dayspmn /= 29) then
             gas_conc_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             gas_conc_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           gas_conc_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         gas_conc_time = model_time
       endif
     endif

      call obtain_interpolator_time_slices &
                          (gas_conc_interp, gas_conc_time)

end subroutine atmos_sulfate_time_vary 




!######################################################################

subroutine atmos_sulfate_endts                    


      call unset_interpolator_time_flag (gocart_emission_interp)

      call unset_interpolator_time_flag (aerocom_emission_interp)

      if (trim(anthro_source) .eq. 'do_anthro') then
          call unset_interpolator_time_flag (anthro_emission_interp)
      endif 

      if (trim(biobur_source) .eq. 'do_biobur') then
          call unset_interpolator_time_flag (biobur_emission_interp)
      endif 

      if (trim(ship_source) .eq. 'do_ship') then
          call unset_interpolator_time_flag (ship_emission_interp)
      endif

      if (trim(aircraft_source) .eq. 'do_aircraft') then
          call unset_interpolator_time_flag (aircraft_emission_interp)
      endif

      call unset_interpolator_time_flag (gas_conc_interp)


end subroutine atmos_sulfate_endts         



!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_sulfate_end">
!<OVERVIEW>
!  The destructor routine for the sulfate module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_sulfate_end
!</TEMPLATE>
 subroutine atmos_sulfate_end

        call interpolator_end (aerocom_emission_interp) 
        call interpolator_end (gocart_emission_interp) 
        call interpolator_end (gas_conc_interp) 
        call interpolator_end (anthro_emission_interp) 
        call interpolator_end (biobur_emission_interp) 
        call interpolator_end (ship_emission_interp) 
        call interpolator_end (aircraft_emission_interp) 
        module_is_initialized = .FALSE.

 end subroutine atmos_sulfate_end
!</SUBROUTINE>
!#######################################################################
!</SUBROUTINE>
!<SUBROUTINE NAME="atmos_DMS_emission">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate dimethyl sulfide emission form the ocean
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_DMS_emission (r, mask, axes, Time)
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
subroutine atmos_DMS_emission (lon, lat, area, ocn_flx_fraction, t_surf_rad, w10m, &
       pwt, DMS_dt, Time, Time_next, is,ie,js,je,kbot)
!
      real, intent(in),    dimension(:,:)           :: lon, lat
      real, intent(in),    dimension(:,:)           :: ocn_flx_fraction
      real, intent(in),    dimension(:,:)           :: t_surf_rad
      real, intent(in),    dimension(:,:)           :: w10m
      real, intent(in),    dimension(:,:)           :: area
      real, intent(in),    dimension(:,:,:)         :: pwt
      real, intent(out),   dimension(:,:,:)         :: DMS_dt
      type(time_type), intent(in)                   :: Time, Time_next    
      integer, intent(in)                           :: is, ie, js, je
      integer, intent(in), dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
      real, dimension(size(DMS_dt,1),size(DMS_dt,2)) :: DMSo, DMS_emis
      integer                                        :: i, j, id, jd, kd
      real                                           :: sst, Sc, conc, w10 
      real                                           :: ScCO2, Akw
      real, parameter                                :: Sc_min=1.

      id=size(dms_dt,1); jd=size(dms_dt,2); kd=size(dms_dt,3)

      dms_dt(:,:,:) =0.0

      DMSo(:,:)=0.0
      call interpolator(gocart_emission_interp, Time, DMSo, &
                       trim(gocart_emission_name(1)), is, js)
! --- Send the DMS data to the diag_manager for output.
      if (id_DMSo > 0 ) &
          used = send_data ( id_DMSo, DMSo, Time_next, is_in=is, js_in=js )

! ****************************************************************************
! *  If ocn_flx_fraction > critical_sea_fraction: DMS_emis = seawaterDMS * transfer velocity * ocn_flx_fraction
! ****************************************************************************
!

      do j = 1, jd
      do i = 1, id
       SST = t_surf_rad(i,j)-273.15     ! Sea surface temperature [Celsius]
       if (ocn_flx_fraction(i,j) .gt. critical_sea_fraction) then

!  < Schmidt number for DMS (Saltzman et al., 1993) >
        Sc = 2674.0 - 147.12*SST + 3.726*(SST**2) - 0.038*(SST**3)
        Sc = max(Sc_min, Sc)

! ****************************************************************************
! *  Calculate transfer velocity in cm/hr  (AKw)                             *
! *                                                                          *
! *  Tans et al. transfer velocity (1990) for CO2 at 25oC (Erickson, 1993)   *
! *                                                                          *
! *  Tans et al. assumed AKW=0 when W10<=3. I modified it to let             *
! *  DMS emit at low windseeds too. Chose 3.6m/s as the threshold.           *
! *                                                                          *
! *  Schmidt number for CO2:       Sc = 600  (20oC, fresh water)             *
! *                                Sc = 660  (20oC, seawater)                *
! *                                Sc = 428  (25oC, Erickson 93)             *
! ****************************************************************************
!

        CONC = DMSo(i,j)

        W10  = W10M(i,j)

! ---  Tans et al. (1990) -----------------
!       ScCO2 = 428.
!       if (W10 .le. 3.6) then
!        AKw = 1.0667 * W10
!       else
!        AKw = 6.4 * (W10 - 3.)
!       end if

! ---  Wanninkhof (1992) ------------------
!       ScCO2 = 660.
!       AKw = 0.31 * W10**2

! ---  Liss and Merlivat (1986) -----------
        ScCO2 = 600.
        if (W10 .le. 3.6) then
         AKw = 0.17 * W10
        else if (W10 .le. 13.) then
         AKw = 2.85 * W10 - 9.65
        else
         AKw = 5.90 * W10 - 49.3
        end if
!------------------------------------------

        if (W10 .le. 3.6) then
         AKw = AKw * ((ScCO2/Sc) ** 0.667)
        else
         AKw = AKw * sqrt(ScCO2/Sc)
        end if

! ****************************************************************************
! *  Calculate emission flux in kg/m2/s                                  *
! *                                                                          *
! *   AKw is in cm/hr:             AKw/100/3600    -> m/sec.                 *
! *   CONC is in nM/L (nM/dm3):    CONC*1E-12*1000 -> kmole/m3.              *
! *   WTM_DMS          : kgDMS/kmol.                                         *
! *   DMS_EMIS         : kgDMS/m2/s.                                         *
! ****************************************************************************
!
        DMS_emis(i,j) = AKw/100./3600. * CONC*1.e-12*1000.* WTM_DMS &
            * (ocn_flx_fraction(i,j))
!
       else               
        DMS_emis(i,j) = 0.

       end if           

      end do
      end do
!--------------------------------------------------------------------
! Update DMS concentration in level kd (where emission occurs)
!--------------------------------------------------------------------
      dms_dt(:,:,kd)=DMS_emis(:,:)/pwt(:,:,kd)* WTMAIR/WTM_DMS
!------------------------------------------------------------------
! DIAGNOSTICS:      DMS surface emission in kg/m2/s     
!--------------------------------------------------------------------
      if (id_DMS_emis > 0) then
        used = send_data ( id_DMS_emis, dms_emis*WTM_S/WTM_DMS, Time_next, &
              is_in=is,js_in=js )
      endif
      if (id_DMS_emis_cmip > 0) then
        used = send_data ( id_DMS_emis_cmip, dms_emis, Time_next, &
              is_in=is,js_in=js )
      endif

end subroutine atmos_DMS_emission
!#######################################################################
!</SUBROUTINE>
!<SUBROUTINE NAME="atmos_SO2_emission">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate SO2 emission from volcanoes, biomass burning,
! anthropogenic sources, aircraft.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_SO2_emission ()
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
subroutine atmos_SOx_emission (lon, lat, area, frac_land, &
       z_pbl, zhalf, phalf, pwt, SO2_dt, SO4_dt, model_time, diag_time, is,ie,js,je,kbot)
!
! This subroutine calculates the tendencies of SO2 and SO4 due to
! their emissions.
! The inventories are based from AEROCOM (cf. Dentener, ACPD, 2006)
! except the aircraft emission.
! The emission of SO4 is assumed to be fe=2.5% of all sulfur emission
! (cf. Dentener, ACPD, 2006). NB. Some authors consider 5%
!
      real, intent(in),    dimension(:,:)           :: lon, lat
      real, intent(in),    dimension(:,:)           :: frac_land
      real, intent(in),    dimension(:,:)           :: area
      real, intent(in),    dimension(:,:)           :: z_pbl
      real, intent(in),    dimension(:,:,:)         :: zhalf, phalf
      real, intent(in),    dimension(:,:,:)         :: pwt
      real, intent(out),   dimension(:,:,:)         :: SO2_dt, SO4_dt
      type(time_type), intent(in)                   :: model_time,diag_time
      integer, intent(in)                           :: is, ie, js, je
      integer, intent(in), dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      integer, parameter :: nlevel_fire = 6
      real, dimension(size(SO4_dt,1),size(SO4_dt,2),size(SO4_dt,3)) :: SO4_emis
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),size(SO2_dt,3)) :: SO2_emis
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),size(SO2_dt,3)) :: &
        so2_aircraft,so2_emis_cont_volc, so2_emis_expl_volc, so2_emis_biobur, &
        so2_emis_ship, so2_emis_road, so2_emis_domestic, so2_emis_industry, &
        so2_emis_power, so2_emis_off_road, so2_emis_ff
! Input emission fields
      real, dimension(size(SO2_dt,1),size(SO2_dt,2)) :: &
             SO2_ff1, SO2_ff2, SO4_ff1, SO4_ff2,&
             SO2_RoadTransport,                         &
             SO2_Off_road,                              &
             SO2_Domestic,                              &
             SO2_Industry,                              &
             SO2_Ship, SO4_ship,                        &
             SO2_Powerplants
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),num_volc_levels) :: &
             SO2_cont_volc,                             &
             SO2_expl_volc
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),nlevel_fire) :: &
             SO2_biobur
! Factors of vertical distribution of emissions
      real, dimension(size(SO2_dt,3)) :: fbb, fa1, fa2, fcv, fev
      real, dimension(size(SO2_dt,3),nlevel_fire) :: ff
! Lower altitude of injection of SO2 from wild fires 
! These values correspond to the AEROCOM input data (cf. Dentener, ACPD, 2006)
      real, dimension(nlevel_fire) :: &
             alt_fire_min=(/0.,100.,500.,1000.,2000.,3000./)
! Upper altitude of injection of SO2 from wild fires 
! These values correspond to the AEROCOM input data (cf. Dentener, ACPD, 2006)
      real, dimension(nlevel_fire) :: &
             alt_fire_max=(/100.,500.,1000.,2000.,3000.,6000./)
! Altitude of injection of surafce anthropogenic emissions
      real :: ze1
! Altitude of injection of SO2 from industries and power plants.   
      real :: ze2
! Emission factor for SO4
      real, parameter :: fe = 0.025

      real :: z1, z2, bltop, fbt, del
      integer  :: i, j, l, id, jd, kd, il, lf
      integer :: ivolc_lev

      id=size(SO2_dt,1); jd=size(SO2_dt,2); kd=size(SO2_dt,3)
!
! Initialize
!
      SO2_dt(:,:,:) = 0.0
      SO4_dt(:,:,:) = 0.0
      SO2_emis(:,:,:) = 0.0
      SO4_emis(:,:,:) = 0.0
! GOCART emissions
      SO2_ff1(:,:)=0.0
      SO2_ff2(:,:)=0.0
      SO4_ff1(:,:)=0.0
      SO4_ff2(:,:)=0.0
! AEROCOM emissions
      SO2_RoadTransport(:,:)=0.0
      SO2_Off_road(:,:)=0.0
      SO2_Domestic(:,:)=0.0
      SO2_Industry(:,:)=0.0
      SO2_Ship(:,:)=0.0
      SO4_Ship(:,:)=0.0
      SO2_Powerplants(:,:)=0.0
      SO2_aircraft(:,:,:)=0.0
      SO2_biobur(:,:,:)=0.0

      SO2_cont_volc(:,:,:)=0.0
      SO2_expl_volc(:,:,:)=0.0
! Arrays for output diagnostics
      so2_aircraft(:,:,:)=0.0
      so2_emis_cont_volc(:,:,:)=0.0
      so2_emis_expl_volc(:,:,:)=0.0
      so2_emis_biobur(:,:,:)=0.0
      so2_emis_ship(:,:,:)=0.0
      so2_emis_road(:,:,:)=0.0
      so2_emis_domestic(:,:,:)=0.0
      so2_emis_industry(:,:,:)=0.0
      so2_emis_power(:,:,:)=0.0
      so2_emis_off_road(:,:,:)=0.0
      so2_emis_ff(:,:,:)=0.0
!
      select case ( trim(runtype))
        case ('gocart')
          if (trim(anthro_source) .eq. 'do_anthro') then
            call interpolator(gocart_emission_interp, model_time, SO2_ff1, &
                         trim(gocart_emission_name(2)), is, js)
            call interpolator(gocart_emission_interp, model_time, SO2_ff2, &
                         trim(gocart_emission_name(3)), is, js)
            call interpolator(gocart_emission_interp, model_time, SO4_ff1, &
                         trim(gocart_emission_name(4)), is, js)
            call interpolator(gocart_emission_interp, model_time, SO4_ff2, &
                         trim(gocart_emission_name(5)), is, js)
          endif
          if (trim(biobur_source) .eq. 'do_biobur') then 
            call interpolator(gocart_emission_interp, model_time, &
                         SO2_biobur(:,:,1),trim(gocart_emission_name(6)), is, js)
          endif
        case ('aerocom')
          if (trim(anthro_source) .eq. 'do_anthro') then
            call interpolator(aerocom_emission_interp, model_time, &
                         SO2_RoadTransport,trim(aerocom_emission_name(1)),is, js)
            call interpolator(aerocom_emission_interp, model_time, SO2_Off_road, &
                         trim(aerocom_emission_name(2)), is, js)
            call interpolator(aerocom_emission_interp, model_time, SO2_Domestic, &
                         trim(aerocom_emission_name(3)), is, js)
            call interpolator(aerocom_emission_interp, model_time, SO2_Industry, &
                         trim(aerocom_emission_name(4)), is, js)
            call interpolator(aerocom_emission_interp, model_time, SO2_ship, &
                         trim(aerocom_emission_name(5)), is, js)
            call interpolator(aerocom_emission_interp, model_time, &
                         SO2_Powerplants, trim(aerocom_emission_name(6)), is, js)
          endif
          if (trim(biobur_source) .eq. 'do_biobur') then
! Wildfire emissions at 6 levels from 0 to 6 km
! (cf. AEROCOM web site or Dentener et al., ACPD, 2006)
            do il=1,nlevel_fire
              call interpolator(aerocom_emission_interp, model_time, &
                         SO2_biobur(:,:,il), &
                         trim(aerocom_emission_name(12+il)), is, js)
            enddo
          endif
        case default
          if (trim(anthro_source) .eq. 'do_anthro') then
            call interpolator(anthro_emission_interp, anthro_time, SO2_ff1(:,:), &
                     trim(anthro_emission_name(1)), is, js)
            call interpolator(anthro_emission_interp, anthro_time, SO4_ff1(:,:), &
                     trim(anthro_emission_name(2)), is, js)
          endif
          if (trim(biobur_source) .eq. 'do_biobur') then
            call interpolator(biobur_emission_interp, biobur_time, SO2_biobur(:,:,1), &
                     trim(biobur_emission_name(1)), is, js)
          endif
          if (trim(ship_source) .eq. 'do_ship') then
            call interpolator(ship_emission_interp, ship_time, SO2_ship(:,:), &
                     trim(ship_emission_name(1)), is, js)
            call interpolator(ship_emission_interp, ship_time, SO4_ship(:,:), &
                     trim(ship_emission_name(2)), is, js)
          endif

      end select
!
! Aircraft emissions
      if (trim(aircraft_source) .eq. 'do_aircraft') then
        call interpolator(aircraft_emission_interp, aircraft_time, &
                     phalf, SO2_aircraft, &
                     trim(aircraft_emission_name(1)), is, js)
      endif
!
! Continuous volcanoes
!
      if (trim(cont_volc_source) .eq. 'do_cont_volc') then
        do ivolc_lev = 1,num_volc_levels
           call interpolator( aerocom_emission_interp, model_time, SO2_cont_volc(:,:,ivolc_lev), &
                              trim(aerocom_emission_name(std_aerocom_emission+ivolc_lev)), is, js )
        end do
      endif
!
! Explosive volcanoes
!
      if (trim(expl_volc_source) .eq. 'do_expl_volc') then
        do ivolc_lev = 1,num_volc_levels
           call interpolator( aerocom_emission_interp, model_time, SO2_expl_volc(:,:,ivolc_lev), &
!                             trim(expl_volc_emission_name(ivolc_lev)), is, js )
                              trim(aerocom_emission_name(std_aerocom_emission+num_volc_levels+ivolc_lev)), is, js )
        end do
      endif

      do j = 1, jd
      do i = 1, id

! --- Assuming biomass burning emission within the PBL -------
        fbb(:) = 0.
        ze1=100.
        ze2=500.
        fbt=0.
        bltop = z_pbl(i,j)
        do l = kd,1,-1
          z1=zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          z2=zhalf(i,j,l)-zhalf(i,j,kd+1)
          if (bltop.lt.z1) exit
          if (bltop.ge.z2) fbb(l)=(z2-z1)/bltop
          if (bltop.gt.z1.and.bltop.lt.z2) fbb(l) = (bltop-z1)/bltop
        enddo
! --- Assuming anthropogenic source L1 emitted below Ze1, and L2
!     emitted between Ze1 and Ze2.
        ff(:,:)=0.
        if (runtype.eq.'aerocom') then
          do l = kd,2,-1
            Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
            Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
            do lf=1,nlevel_fire
              del=alt_fire_max(lf)-alt_fire_min(lf)
              if (del.gt.0. .and. &
                  Z1.lt.alt_fire_max(lf).and.Z2.gt.alt_fire_min(lf) ) then
                if (Z1.ge.alt_fire_min(lf)) then
                  if (Z2 .lt. alt_fire_max(lf)) then
                    ff(l,lf)=(Z2-Z1)/del
                  else
                    ff(l,lf)=(alt_fire_max(lf)-z1)/del
                  endif
                else
                  if (Z2.le.alt_fire_max(lf)) then
                    ff(l,lf) = (Z2-alt_fire_min(lf))/del
                  else
                    ff(l,lf)=1.
                  endif
                endif
              endif
            enddo
          enddo
        endif
! --- Volcanic SO2 source ----
! --- For continuous and explosive volcanoes, calculate the fraction of emission
! --- for each vertical level
      if (trim(cont_volc_source) == 'do_cont_volc') then
        do ivolc_lev = 1,num_volc_levels
          fcv(:)=0.
          do l = kd,2,-1
            Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
            Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
            del=volc_altitude_edges(ivolc_lev+1)-volc_altitude_edges(ivolc_lev)
            if (del>0. .and. &
                Z1<volc_altitude_edges(ivolc_lev+1) .and. Z2>volc_altitude_edges(ivolc_lev) ) then
              if (Z1 >= volc_altitude_edges(ivolc_lev)) then
                if (Z2 < volc_altitude_edges(ivolc_lev+1)) then
                  fcv(l)=(Z2-Z1)/del
                else
                  fcv(l)=(volc_altitude_edges(ivolc_lev+1)-Z1)/del
                endif
              else
                if (Z2 <= volc_altitude_edges(ivolc_lev+1)) then
                  fcv(l)=(Z2-volc_altitude_edges(ivolc_lev))/del
                else
                  fcv(l)=1.
                endif
              endif
            endif
          enddo
          so2_emis_cont_volc(i,j,:) = so2_emis_cont_volc(i,j,:) + fcv(:) * SO2_cont_volc(i,j,ivolc_lev)
        end do
      endif
! --- For explosive volcanoes, calculate the fraction of emission for
! --- each vertical levels
      if (trim(expl_volc_source) == 'do_expl_volc') then 
        do ivolc_lev = 1,num_volc_levels
          fev(:)=0.
          do l = kd,2,-1
            Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
            Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
            del=volc_altitude_edges(ivolc_lev+1)-volc_altitude_edges(ivolc_lev)
            if (del>0. .and. &
                Z1<volc_altitude_edges(ivolc_lev+1).and.Z2>volc_altitude_edges(ivolc_lev) ) then
              if (Z1 >= volc_altitude_edges(ivolc_lev)) then
                if (Z2 < volc_altitude_edges(ivolc_lev+1)) then
                  fev(l)=(Z2-Z1)/del
                else
                  fev(l)=(volc_altitude_edges(ivolc_lev+1)-Z1)/del
                endif
              else
                if (Z2 <= volc_altitude_edges(ivolc_lev+1)) then
                  fev(l)=(Z2-volc_altitude_edges(ivolc_lev))/del
                else
                  fev(l)=1.
                endif
              endif
            endif
          enddo
          so2_emis_expl_volc(i,j,:) = so2_emis_expl_volc(i,j,:) + fev(:) * SO2_expl_volc(i,j,ivolc_lev)
        end do
      endif
! --- For fosil fuel emissions, calculate the fraction of emission for
! --- each vertical levels
        fa1(:) = 0.
        fa2(:) = 0.
        do l = kd,2,-1
          Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
          if (Z2.ge.0.and.Z1.lt.ze1) then
            if (Z1.gt.0) then
              if (Z2.lt.ze1) then
                fa1(l)=(Z2-Z1)/ze1
              else
                fa1(l)=(ze1-Z1)/ze1
              endif
            else
              if (Z2.le.ze1) then
                fa1(l)=Z2/ze1
              else
                fa1(l)=1.
              endif
            endif
          endif

          if (Z2.ge.ze1.and.z1.lt.ze2) then
            if (Z1.gt.Ze1) then
              if (Z2.lt.ze2) then
                fa2(l)=(z2-z1)/(ze2-ze1)
              else
                fa2(l)=(ze2-z1)/(ze2-ze1)
              endif
            else
              if (Z2.le.ze2) then
                fa2(l)=(z2-ze1)/(ze2-ze1)
              else
                fa2(l)=1.
              endif
            endif
          endif
          if (Z1.gt.Ze2) exit
        enddo
! SO2_emis: [kgSO2/m2/s]
!       Assuming that 1g of SO2 is emitted from 1kg of fuel: 1.e-3
        SO2_emis(i,j,:) = so2_aircraft_EI * SO2_aircraft(i,j,:) &
             + so2_emis_cont_volc(i,j,:) + so2_emis_expl_volc(i,j,:)
!
        select case (trim(runtype))
          case ('aerocom')
            do lf = 1, nlevel_fire
              so2_emis_biobur(i,j,:) = so2_emis_biobur(i,j,:) + &
                                       ff(:,lf)*SO2_biobur(i,j,lf)
            enddo
            so2_emis_road(i,j,:)     = fa1(:) * SO2_RoadTransport(i,j)
            so2_emis_off_road(i,j,:) = fa1(:) * SO2_off_Road(i,j)
            so2_emis_domestic(i,j,:) = fa1(:) * SO2_domestic(i,j)
            so2_emis_ship(i,j,:)     = fa1(:) * SO2_ship(i,j)
            so2_emis_industry(i,j,:) = fa2(:) * SO2_industry(i,j)
            so2_emis_power(i,j,:)    = fa2(:) * SO2_Powerplants(i,j)

            SO2_emis(i,j,:) = SO2_emis(i,j,:) + so2_emis_biobur(i,j,:)

            so2_emis_ff(i,j,:) =                                  &
                 + so2_emis_road(i,j,:)                           &
                 + so2_emis_off_road(i,j,:)                       &
                 + so2_emis_domestic(i,j,:)                       &
                 + so2_emis_ship(i,j,:)                           &
                 + so2_emis_industry(i,j,:)                       &
                 + so2_emis_power(i,j,:)
             SO2_emis(i,j,:) = SO2_emis(i,j,:) + so2_emis_ff(i,j,:)
          case ('gocart')
!
! GOCART assumes continent based emission index for sulfate:
!    Anthropogenic SOx emission from GEIA 1985.
!    Assuming:   Europe:      5.0% SOx emission is SO4;
!                US + Canada: 1.4% SOx emission is SO4;
!                The rest:    2.5% SOx emission is SO4.
            so2_emis_ff(i,j,:)=fa1(:) * SO2_ff1(i,j) + fa2(:) * SO2_ff2(i,j)
            so2_emis_biobur(i,j,:) = fbb(:) * SO2_biobur(i,j,1)
            SO2_emis(i,j,:) = SO2_emis(i,j,:) &
               + so2_emis_biobur(i,j,:)       &
               + so2_emis_ff(i,j,:)
            SO4_emis(i,j,:) = &
               fa1(:) * SO4_ff1(i,j) + fa2(:) * SO4_ff2(i,j)
          case default
            so2_emis_ff(i,j,:)=fa1(:) * SO2_ff1(i,j) + fa2(:) * SO2_ff2(i,j)
            so2_emis_biobur(i,j,:) = fbb(:) * SO2_biobur(i,j,1)
            so2_emis_ship(i,j,:)     = fa1(:) * SO2_ship(i,j)
            SO2_emis(i,j,:) = SO2_emis(i,j,:) &
               + so2_emis_biobur(i,j,:)       &
!++lwh
               + so2_emis_ship(i,j,:)         &
!--lwh
               + so2_emis_ff(i,j,:)
            SO4_emis(i,j,:) = fa1(:)*(SO4_ff1(i,j)+SO4_ship(i,j))
        end select

      end do   ! end i loop
      end do   ! end j loop
!
      SO2_dt(:,:,:)= SO2_emis(:,:,:)/pwt(:,:,:)*WTMAIR/WTM_SO2
      SO4_dt(:,:,:)= SO4_emis(:,:,:)/pwt(:,:,:)*WTMAIR/WTM_SO4

!------------------------------------------------------------------
! DIAGNOSTICS:      SO2 and SO4 emission in kg/timestep
!--------------------------------------------------------------------
      if (id_so2_emis > 0) then
        used = send_data ( id_so2_emis, so2_emis*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js,ks_in=1)
      endif
      if (id_so2_aircraft > 0) then
        used = send_data ( id_so2_aircraft, &
              so2_aircraft*so2_aircraft_EI*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_cont_volc > 0) then
        used = send_data ( id_so2_cont_volc, so2_emis_cont_volc*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_expl_volc > 0) then
        used = send_data ( id_so2_expl_volc, so2_emis_expl_volc*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_biobur > 0) then
        used = send_data ( id_so2_biobur, so2_emis_biobur*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js,ks_in=1)
      endif
      if (id_so2_ship > 0) then
        used = send_data ( id_so2_ship, so2_emis_ship*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_road > 0) then
        used = send_data ( id_so2_road, so2_emis_road*WTM_S/WTM_so2,  &
               diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_domestic > 0) then
        used = send_data ( id_so2_domestic, so2_emis_domestic*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_industry > 0) then
        used = send_data ( id_so2_industry, so2_emis_industry*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_power > 0) then
        used = send_data ( id_so2_power, so2_emis_power*WTM_S/WTM_so2, &
              diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_off_road > 0) then
        used = send_data ( id_so2_Off_road, so2_emis_off_road*WTM_S/WTM_so2, &
               diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so2_ff > 0) then
        used = send_data ( id_so2_ff, so2_emis_ff*WTM_S/WTM_so2, &
               diag_time, is_in=is,js_in=js, ks_in=1)
      endif
      if (id_so4_emis > 0) then
        used = send_data ( id_so4_emis, so4_emis*WTM_S/WTM_so4, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif

end subroutine atmos_SOx_emission
!</SUBROUTINE>
!-----------------------------------------------------------------------
!#######################################################################
      subroutine atmos_SOx_chem(pwt,temp,pfull, phalf, dt, lwc, &
        jday,hour,minute,second,lat,lon, &
        SO2, SO4, DMS, MSA, H2O2, oh_vmr, &
        SO2_dt, SO4_dt, DMS_dt, MSA_dt, H2O2_dt, &
        model_time,diag_time,is,ie,js,je,kbot)
!
      real, intent(in)                   :: dt
      integer, intent(in)                :: jday, hour,minute,second
      real, intent(in),  dimension(:,:)  :: lat, lon  ! [radi
      real, intent(in), dimension(:,:,:) :: pwt
      real, intent(in), dimension(:,:,:) :: lwc
      real, intent(in), dimension(:,:,:) :: temp, pfull, phalf
      real, intent(in), dimension(:,:,:) :: SO2, SO4, DMS, MSA, H2O2
      real, intent(inout), dimension(:,:,:) :: oh_vmr
      real, intent(out),dimension(:,:,:) :: SO2_dt,SO4_dt,DMS_dt,MSA_dt,H2O2_dt

      type(time_type), intent(in)                    :: model_time,diag_time
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      integer :: i,j,k,id,jd,kd
      integer                                    :: istep, nstep
!!! Input fields from interpolator
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: pH
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: O3_mmr
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: no3_conc
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: oh_conc
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: jh2o2
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: ho2_conc
!!! Time varying fields
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: no3_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: oh_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) ::jh2o2_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: ho2_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: &
               SO4_oh_prod, SO4_o3_prod, SO4_h2o2_prod

      real, dimension(size(pfull,1),size(pfull,2)) :: &
               xu, dayl, h, hl, hc, hred
      real, dimension(size(pfull,1),size(pfull,2)) :: fac_NO3, fact_NO3
      real, dimension(size(pfull,1),size(pfull,2)) :: fac_OH , fact_OH 
      real, dimension(size(pfull,1),size(pfull,2)) :: fac_HO2            
      real, parameter                            :: A0 = 0.006918
      real, parameter                            :: A1 = 0.399912
      real, parameter                            :: A2 = 0.006758
      real, parameter                            :: A3 = 0.002697
      real, parameter                            :: B1 = 0.070257
      real, parameter                            :: B2 = 0.000907
      real, parameter                            :: B3 = 0.000148
      real                                       :: decl, hd, x
      real :: f, f1, tk, rho_air
      real :: SO2_0,SO4_0,MSA_0,DMS_0,H2O2_0    ! initial concentrations
      real :: xSO2,xSO4,xMSA,xDMS,xH2O2,xno3,xo3,xoh,xho2,xjh2o2 ! update conc.
      real :: rk0, rk1, rk2, rk3  ! kinetic rates
      real :: work1, xk, xe, x2, xph
      real :: heh2o2, h2o2g, rah2o2, px, heso2, so2g, heo3, o3g, rao3
      real :: pso4a, pso4b
      real :: xlwc, xhnm, ccc1, ccc2
      real :: pmsa, pso2, ph2o2    ! chemical production terms
      real :: ldms, lso2, lh2o2          ! chemical loss terms
      real :: o2
      real, parameter        :: small_value=1.e-21
      real, parameter        :: t0 = 298.
      real, parameter        :: Ra = 8314./101325.
      real, parameter        :: xkw = 1.e-14 ! water acidity
      real, parameter        :: const0 = 1.e3/6.022e23


! Local grid sizes
      id=size(pfull,1) ; jd=size(pfull,2) ; kd=size(pfull,3)

      so2_dt(:,:,:) = 0.0
      so4_dt(:,:,:) = 0.0
      dms_dt(:,:,:) = 0.0
      msa_dt(:,:,:) = 0.0
      h2o2_dt(:,:,:) = 0.0
      SO4_h2o2_prod(:,:,:)=0.0
      SO4_o3_prod(:,:,:)=0.0
      SO4_oh_prod(:,:,:)=0.0


      OH_conc(:,:,:)=1.e5  ! molec/cm3
      call interpolator(gas_conc_interp, gas_conc_time, phalf, OH_conc, &
                       trim(gas_conc_name(1)), is, js)

      HO2_conc(:,:,:)=1.e6  ! molec/cm3
      call interpolator(gas_conc_interp, gas_conc_time, phalf, HO2_conc, &
                       trim(gas_conc_name(2)), is, js)

      NO3_conc(:,:,:)=0.0  ! molec/cm3
      call interpolator(gas_conc_interp, gas_conc_time, phalf, NO3_conc, &
                       trim(gas_conc_name(3)), is, js)

      O3_mmr(:,:,:)=0  ! Ozone mass mixing ratio
      call interpolator(gas_conc_interp, gas_conc_time, phalf, O3_mmr, &
                       trim(gas_conc_name(4)), is, js)
      O3_mmr(:,:,:)=O3_mmr(:,:,:)*WTM_O3/WTMAIR

      jH2O2(:,:,:)=1.e-6 ! s-1
      call interpolator(gas_conc_interp, gas_conc_time, phalf, jH2O2, &
                       trim(gas_conc_name(5)), is, js)

      pH(:,:,:)=1.e-5
      call interpolator(gas_conc_interp, gas_conc_time, phalf, pH, &
                       trim(gas_conc_name(6)), is, js)

      x = 2. *pi *float(jday-1)/365.
      decl = A0 - A1*cos(  X) + B1*sin(  X) - A2*cos(2.*X) + B2*sin(2.*X) &
           - A3*cos(3.*X) + B3*sin(3.*X)
      xu(:,:) = -tan(lat(:,:))*tan(decl)

      where ( xu > -1 .and. xu < 1 ) dayl=acos(xu)/pi
      where ( xu <= -1 ) dayl = 1.
      where ( xu >= 1 ) dayl = 0.
!   Calculate normalization factors for OH and NO3 such that
!   the diurnal average respect the monthly input values.
      hd=0.
      fact_OH(:,:)  = 0.
      fact_NO3(:,:) = 0.
      nstep = int(24.*3600./dt)
      do istep=1,nstep
        hd=hd+dt/3600./24.
        hl(:,:) = pi*(1.-dayl(:,:))
        hc(:,:) = pi*(1.+dayl(:,:))
        h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
        where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fact_OH  = fact_OH + amax1(0.,sin(pi*hred)/2.)/nstep
        elsewhere
! Nightime
          fact_NO3 = fact_NO3 + amax1(0.,amin1(1.,(1.-dayl)))/nstep
        endwhere
      enddo

      hd=amax1(0.,amin1(1.,(hour+minute/60.+second/3600.)/24.))
      hl(:,:) = pi*(1.-dayl(:,:))
      hc(:,:) = pi*(1.+dayl(:,:))
      h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
      fac_OH(:,:)  = 0.
      fac_NO3(:,:) = 0.
      fac_HO2(:,:) = 1.
      where ( h.ge.hl .and. h.lt.hc )
! Daytime
          fac_NO3 = 0.
          hred=(h-hl)/(hc-hl)
          fac_OH  = amax1(0.,sin(pi*hred)/2.)/fact_OH
      elsewhere
! Nightime
          fac_NO3 = amax1(0.,amin1(1.,(1.-dayl)))/fact_NO3
          fac_OH  = 0.
      endwhere


! < Factor to convert AIRDEN from kgair/m3 to molecules/cm3: >
      f  = 1000. / WTMAIR * 6.022e23 * 1.e-6

      do k = 1, kd
      do j = 1, jd
      do i = 1, id
       tk    = temp(i,j,k)
       rho_air = pfull(i,j,k)/tk/RDGAS             ! Air density [kg/m3]
       xhnm  = rho_air * f
       O2    = xhnm * 0.21
       xlwc  = lwc(i,j,k)*rho_air *1.e-3
       DMS_0 = max(0.,DMS(i,j,k))
       MSA_0 = max(0.,MSA(i,j,k))
       SO4_0 = max(0.,SO4(i,j,k))
       SO2_0 = max(0.,SO2(i,j,k))
       H2O2_0= max(0.,H2O2(i,j,k))
       xSO2  = SO2_0
       xSO4  = SO4_0
       xH2O2 = H2O2_0
       xDMS  = DMS_0
       xMSA  = MSA_0
       xph   = max(1.e-7,       pH(i,j,k))
       xoh   = max(0.         , OH_conc(i,j,k)  *fac_OH(i,j))
       xho2  = max(0.         , HO2_conc(i,j,k) *fac_HO2(i,j))
       xjh2o2= max(0.         , jH2O2(i,j,k)    *fac_OH(i,j))
       xno3  = max(0.         , NO3_conc(i,j,k) *fac_NO3(i,j))
       xo3   = max(small_value, O3_mmr(i,j,k))
       oh_diurnal(i,j,k)=xoh
       oh_vmr(i,j,k)=xoh/xhnm
       no3_diurnal(i,j,k)=xno3
       ho2_diurnal(i,j,k)=xho2
       jh2o2_diurnal(i,j,k)=xjh2o2
! ****************************************************************************
! *  H2O2 production by HO2 + HO2 reactions
! ****************************************************************************
       PH2O2=(2.2e-13*exp(619./tk)+xhnm*1.9e-33*exp(980./tk))* xHO2**2 /xhnm
! ****************************************************************************
! *  H2O2 loss by OH and photodissociation
! ****************************************************************************
       LH2O2= ( 2.9e-12*exp(-160./tk)* xOH + xjH2O2 )
       if (LH2O2 .gt. 0.) then
         xH2O2= H2O2_0 * exp(-LH2O2*dt) + PH2O2*(1.-exp(-LH2O2*dt))/LH2O2
       else
         xH2O2= H2O2_0 + PH2O2 * dt
       endif
! ****************************************************************************
! *  (1) DMS + OH:  RK1 - addition channel;  RK2 - abstraction channel.      *
! ****************************************************************************
       rk1 = (1.7e-42 * exp(7810./TK) * O2) /   &
              (1. + 5.5e-31 * exp(7460./TK) * O2 ) * xoh
       rk2 = 1.2e-11*exp(-260./TK) * xoh
! ****************************************************************************
! *  (2) DMS + NO3 (only happen at night):                                   *
! ****************************************************************************
!  < XNO3 fields are in molecules/cm3.        >
        rk3 = 1.9e-13 * exp(500./TK) * xno3
! ****************************************************************************
! *  Update DMS concentration after gas phase chemistry                      *
! ****************************************************************************
       LDMS = RK1 + RK2 + RK3
       if ( LDMS .gt. 0. ) then
         xDMS = DMS_0 * exp( - LDMS*dt)
       endif
! ****************************************************************************
! *  Update MSA concentration after gas phase chemistry                      *
! ****************************************************************************
       PMSA = RK1*0.25 * xDMS
       xMSA = MSA_0 + PMSA * dt
! ****************************************************************************
! *  SO2 oxydation by OH
! ****************************************************************************
       PSO2 = ( RK1*0.75 + RK2 + RK3 ) * xDMS
       rk0 = 3.0E-31 * (300./TK)**3.3
       rk1 = rk0 * xhnm / 1.5e-12
       f1 = ( 1.+ ( log10(rk1) )**2 )**(-1)
       LSO2 = ( rk0 * xhnm / (1.+ rk1) ) * 0.6**f1 * xoh
! ****************************************************************************
! *  Update SO2 concentration after gas phase chemistry                      *
! ****************************************************************************
       xSO2 = SO2_0 + dt * ( PSO2 - LSO2 * SO2_0)
       if (xSO2 .lt. 0.) then
        xSO2 = SO2_0 * exp(-LSO2*dt) + PSO2 * (1.-exp(-LSO2*dt)) / LSO2
       end if
! ****************************************************************************
! *  Update SO4 concentration after gas phase chemistry                      *
! ****************************************************************************
       xso4 = SO4_0 + LSO2*xso2 * dt
! ****************************************************************************
! < Cloud chemistry (above 258K): >
       work1 = (t0 - tk)/(tk*t0)
!-----------------------------------------------------------------------
!         ... h2o2
!-----------------------------------------------------------------------
       xk = 7.4e4   *exp( 6621.* work1 )
       xe = 2.2e-12 *exp(-3730.* work1 )
       heh2o2  = xk*(1. + xe/xph)
       px = heh2o2 * Ra * tk * xlwc
       h2o2g = xh2o2 /(1.+px) 
!-----------------------------------------------------------------------
!         ... so2
!-----------------------------------------------------------------------
       xk = 1.23   * exp(3120. * work1 )
       xe = 1.7e-2 * exp(2090. * work1 ) 
       x2 = 6.0e-8 * exp(1120. * work1 )
       heso2 = xk*(1. + xe/xph *(1. + x2/xph) ) 
!       heso2 = 1.e2 ! xk*(1. + xe/xph *(1. + x2/xph) )
       px = heso2 * Ra * tk * xlwc
       so2g = xso2/(1.+px)
!-----------------------------------------------------------------------
!         ... o3
!-----------------------------------------------------------------------
       xk = 1.15e-2 * exp( 2560. * work1 )
       heo3 = xk
       px = heo3 * Ra * tk *xlwc
       o3g = xo3 / (1.+px) 
!-----------------------------------------------
!       ... Aqueous phase reaction rates
!           SO2 + H2O2 -> SO4
!           SO2 + O3   -> SO4
!-----------------------------------------------

!------------------------------------------------------------------------
!       ... S(IV) (HSO3) + H2O2
!------------------------------------------------------------------------
            rah2o2 = 8.e4 * EXP( -3650.*work1 )  / (.1 + xph)

!------------------------------------------------------------------------
!        ... S(IV)+ O3
!------------------------------------------------------------------------
            rao3   = 4.39e11 * EXP(-4131./tk)  &
                  + 2.56e3  * EXP(-996. /tk) /xph

!-----------------------------------------------------------------
!       ... Prediction after aqueous phase
!       so4
!       When Cloud is present
!
!       S(IV) + H2O2 = S(VI)
!       S(IV) + O3   = S(VI)
!
!       reference:
!           (1) Seinfeld
!           (2) Benkovitz
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!       ... S(IV) + H2O2 = S(VI)
!-----------------------------------------------------------------
       ccc1=0.
       ccc2=0.
       if( xlwc >= 1.e-8 ) then                    ! when cloud is present
               pso4a = rah2o2 * heh2o2*h2o2g  &
                             * heso2 *so2g            ! [M/s]
               pso4a = pso4a       &                    ! [M/s] =  [mole/L(w)/s]
                    * xlwc       &                    ! [mole/L(a)/s]
                    / const0     &                    ! [/L(a)/s]
                    / xhnm                            ! [mixing ratio/s]

          ccc1 = pso4a*dt
          ccc1 = max(min(ccc1,xso2,xh2o2), 0.)
          xso4 = xso4 + ccc1
          xh2o2 = max(xh2o2 - ccc1, small_value)
          xso2 =  max(xso2  - ccc1, small_value)
!          ccc1 = max(ccc1, 0.)
!          if( xh2o2 > xso2 ) then
!              if( ccc1 > xso2 ) then
!                  xso4  = xso4 + xso2
!                  xso2  = small_value
!                  xh2o2 = xh2o2 - xso2
!              else
!                  xso4  = xso4  + ccc1
!                  xh2o2 = xh2o2 - ccc1
!                  xso2  = xso2  - ccc1
!              end if
!          else
!               if( ccc1 > xh2o2 ) then
!                   xso4  = xso4 + xh2o2
!                   xso2  = xso2 - xh2o2
!                   xh2o2 = small_value
!               else
!                   xso4  = xso4  + ccc1
!                   xh2o2 = xh2o2 - ccc1
!                   xso2  = xso2  - ccc1
!               end if
!          end if


!-----------------------------------------------
!       ... S(IV) + O3 = S(VI)
!-----------------------------------------------
           pso4b = rao3 * heo3*o3g * heso2*so2g       ! [M/s]
           pso4b = pso4b        &        ! [M/s] =  [mole/L(w)/s]
                * xlwc        &        ! [mole/L(a)/s]
                / const0      &        ! [/L(a)/s]
                / xhnm                 ! [mixing ratio/s]
 
           ccc2 = pso4b*dt
            ccc2 = max(min(ccc2, xso2), 0.)               ! mozart2
            xso4 = xso4 + ccc2                           ! mozart2
            xso2 = max(xso2 - ccc2, small_value)         ! mozart2
!          ccc2 = max(ccc2, 0.)
!          if( ccc2 > xso2 ) then
!             xso4 = xso4 + xso2
!             xso2 = small_value
!          else
!             xso4 = xso4  + ccc2
!             xso2 = xso2  - ccc2
!             xso2 = max(xso2, small_value)
!          end if
       end if
       MSA_dt(i,j,k) = (xMSA-MSA_0)/dt
       DMS_dt(i,j,k) = (xDMS-DMS_0)/dt
       SO2_dt(i,j,k) = (xso2-SO2_0)/dt
       SO4_dt(i,j,k) = (xso4-SO4_0)/dt
       H2O2_dt(i,j,k)= (xh2o2-H2O2_0)/dt
       SO4_oh_prod(i,j,k)=LSO2*xso2
       SO4_o3_prod(i,j,k)=ccc2/dt
       SO4_h2o2_prod(i,j,k)=ccc1/dt
      end do
      end do
      end do
      if ( id_NO3 > 0) then
        used = send_data ( id_NO3, NO3_diurnal, &
                           diag_time,is_in=is,js_in=js,ks_in=1)
      endif
      if ( id_OH > 0) then
        used = send_data ( id_OH, OH_diurnal, &
                           diag_time, is_in=is, js_in=js,ks_in=1 )
      endif
      if ( id_HO2 > 0) then
        used = send_data ( id_HO2, HO2_diurnal, &
                           diag_time, is_in=is, js_in=js,ks_in=1 )
      endif
      if ( id_jH2O2 > 0) then
        used = send_data ( id_jH2O2, jH2O2_diurnal, &
                           diag_time, is_in=is, js_in=js,ks_in=1 )
      endif
      if (id_ph > 0) then
        used = send_data ( id_ph, ph, &
                           diag_time,is_in=is,js_in=js,ks_in=1)
      endif
      if (id_o3 > 0) then
        used = send_data ( id_o3, o3_mmr, &
                           diag_time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_SO2_chem > 0) then
        used = send_data ( id_SO2_chem, &
              SO2_dt*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_SO4_chem > 0) then
        used = send_data ( id_SO4_chem, &
              SO4_dt*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_SO4_oh_prod > 0) then
        used = send_data ( id_SO4_oh_prod, &
              SO4_oh_prod*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif
      if (id_SO4_o3_prod > 0) then
        used = send_data ( id_SO4_o3_prod, &
              SO4_o3_prod*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif
      if (id_SO4_h2o2_prod > 0) then
        used = send_data ( id_SO4_h2o2_prod, &
              SO4_h2o2_prod*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_DMS_chem > 0) then
        used = send_data ( id_DMS_chem, &
              DMS_dt*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_MSA_chem > 0) then
        used = send_data ( id_MSA_chem, &
              MSA_dt*pwt*WTM_S/WTMAIR, diag_time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_H2O2_chem > 0) then
        used = send_data ( id_H2O2_chem, &
              H2O2_dt*pwt, diag_time,is_in=is,js_in=js,ks_in=1)
      endif
end subroutine atmos_SOx_chem

end module atmos_sulfate_mod
