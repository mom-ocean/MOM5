module atmos_co2_mod
! <CONTACT EMAIL="rdslater@splash.princeton.edu">
!   Richard D. Slater
! </CONTACT>

! <REVIEWER EMAIL="none@nowhere.dot">
!   none yet
! </REVIEWER>


! <OVERVIEW>
! </OVERVIEW>

! <DESCRIPTION>
! </DESCRIPTION>


use mpp_mod, only: input_nml_file 
use              fms_mod, only : file_exist, write_version_number,    &
                                 mpp_pe, mpp_root_pe,                 &
                                 close_file, stdlog, stdout,          &
                                 check_nml_error, error_mesg,         &
                                 open_namelist_file, FATAL, NOTE, WARNING

use   tracer_manager_mod, only : get_tracer_index, tracer_manager_init
use    field_manager_mod, only : MODEL_ATMOS
use     diag_manager_mod, only : register_diag_field, send_data
use     time_manager_mod, only : time_type
use        constants_mod, only : WTMCO2, WTMAIR, WTMC
use    data_override_mod, only : data_override
use              mpp_mod, only : mpp_pe, mpp_root_pe


implicit none

private

!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_co2_sourcesink
public  atmos_co2_emissions
public  atmos_co2_rad
public  atmos_co2_gather_data
public  atmos_co2_flux_init
public  atmos_co2_init
public  atmos_co2_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

character(len=48), parameter    :: mod_name = 'atmos_co2_mod'

integer, save   :: ind_co2_flux = 0
integer, save   :: ind_co2  = 0
integer, save   :: ind_sphum = 0
integer         :: id_co2restore, id_pwt, id_co2_mol_emiss, id_co2_emiss_orig
real            :: radiation_co2_dvmr = -1

!---------------------------------------------------------------------
!-------- namelist  ---------

real     :: restore_tscale   = -1
integer  :: restore_klimit   = -1
logical  :: do_co2_restore   = .false.
logical  :: co2_radiation_override = .false.
logical  :: do_co2_emissions = .false.

namelist /atmos_co2_nml/  &
          do_co2_restore, restore_tscale, restore_klimit,  &
          co2_radiation_override, do_co2_emissions

!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make the
!  following changes.
!
!  Add an integer variable below for each additional tracer. This should
!  be initialized to zero. 
!
!  Add id_tracername for each additional tracer. These are used in
!  initializing and outputting the tracer fields.
!
!-----------------------------------------------------------------------

!PUBLIC VARIABLES
public :: co2_radiation_override, do_co2_emissions

logical :: module_is_initialized = .FALSE.


!---- version number -----
character(len=128) :: version = '$$'
character(len=128) :: tagname = '$$'
!-----------------------------------------------------------------------

contains

!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_sourcesink">
!<OVERVIEW>
!  A subroutine to calculate the internal sources and sinks of carbon dioxide.
!
! do_co2_restore   = logical to turn co2_restore on/off: default = .false.
! restore_co2_dvmr = partial pressure of co2 to which to restore  (mol/mol)
! restore_klimit   = atmospheric level to which to restore starting from top
! restore_tscale   = timescale in seconds with which to restore
!
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate the sources and sinks of carbon dixoide.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_co2_sourcesink (Time, dt,  pwt, co2, sphum, co2_restore)
!
!</TEMPLATE>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model timestep.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav (kg/m2)
!   </IN>
!   <IN NAME="co2" TYPE="real" DIM="(:,:,:)">
!     The array of the carbon dioxide mixing ratio (kg co2/kg moist air)
!   </IN>
!   <IN NAME="sphum" TYPE="real" DIM="(:,:,:)">
!     The array of the specific humidity mixing ratio (kg/kg)
!   </IN>

!   <OUT NAME="co2_restore" TYPE="real" DIM="(:,:,:)">
!     The array of the restoring tendency of the carbon dioxide mixing ratio.
!   </OUT>
!

subroutine atmos_co2_sourcesink(is, ie, js, je, Time, Time_next, dt, pwt, co2, sphum, co2_restore)

   integer, intent(in)                 :: is, ie, js, je
   type (time_type),      intent(in)   :: Time, Time_next
   real, intent(in)                    :: dt
   real, intent(in),  dimension(:,:,:) :: pwt          ! kg/m2
   real, intent(in),  dimension(:,:,:) :: co2          ! moist mmr
   real, intent(in),  dimension(:,:,:) :: sphum        
   real, intent(out), dimension(:,:,:) :: co2_restore
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_sourcesink'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


integer   :: i,j,k,id,jd,kd, logunit
logical   :: sent
logical   :: used
real      :: restore_co2_dvmr = -1

!-----------------------------------------------------------------------

id=size(co2,1); jd=size(co2,2); kd=min(size(co2,3),restore_klimit)

co2_restore(:,:,:)=0.0

logunit=stdlog()
if (ind_co2 > 0 .and. do_co2_restore) then

! input is in dvmr (mol/mol)
  call data_override('ATM', 'co2_dvmr_restore', restore_co2_dvmr, Time, override=used)
  if (.not. used) then
    call error_mesg (trim(error_header), ' data override needed for co2_dvmr_restore ', FATAL)
  endif
!  if (mpp_pe() == mpp_root_pe() ) &
!      write (logunit,*)' atmos_co2_sourcesink: mean restore co2_dvmr   = ', restore_co2_dvmr


  if (restore_tscale .gt. 0 .and. restore_co2_dvmr .ge. 0.0) then
! co2mmr = (wco2/wair) * co2vmr;  wet_mmr is approximated as dry_mmr * (1-Q)
    do k=1,kd
      do j=1,jd
        do i=1,id
! convert restore_co2_dvmr to wet mmr and get tendency
          co2_restore(i,j,k) = (restore_co2_dvmr * (WTMCO2/WTMAIR) * (1.0 - &
                           sphum(i,j,k)) - co2(i,j,k))/restore_tscale
        enddo
      enddo
    enddo

! restoring diagnostic in moles co2/m2/sec 
! pwt is moist air, so no need to divide by 1-sphum here
    if (id_co2restore > 0) sent = send_data (id_co2restore, co2_restore  *  &
                                         pwt / (WTMCO2*1.e-3), Time_next, is_in=is,js_in=js)
  endif

!else
!  if (mpp_pe() == mpp_root_pe() ) &
!      write (logunit,*)' atmos_co2_sourcesink: CO2 restoring not active: ',do_co2_restore
endif

!! add pwt as a diagnostic
if (id_pwt > 0) sent = send_data (id_pwt, pwt, Time_next, is_in=is,js_in=js)


end subroutine atmos_co2_sourcesink
!</SUBROUTINE >


!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_rad">

!<OVERVIEW>
! Subroutine to get global avg co2 to be used in radiation.
! input co2 field is from data override 
!</OVERVIEW>

 subroutine atmos_co2_rad(Time, radiation_co2_dvmr)

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!
   type (time_type),      intent(in)    :: Time
   real,                  intent(inout) :: radiation_co2_dvmr
!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_rad'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
integer :: logunit
logical   :: used
!-----------------------------------------------------------------------

logunit=stdlog()

if (ind_co2 > 0 .and. co2_radiation_override) then

! input is in dvmr (mol/mol)

  call data_override('ATM', 'co2_dvmr_rad', radiation_co2_dvmr, Time, override=used)
  if (.not. used) then
    call error_mesg (trim(error_header), ' data override needed for co2_dvmr_rad ', FATAL)
  endif
!  if (mpp_pe() == mpp_root_pe() ) &
!      write (logunit,*)' atmos_co2_rad       : mean radiation co2_dvmr = ', radiation_co2_dvmr

!else
!  if (mpp_pe() == mpp_root_pe() ) &
!      write (logunit,*)' atmos_co2_rad: CO2 radiation override not active: ',co2_radiation_override
endif


!-----------------------------------------------------------------------

end subroutine atmos_co2_rad
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_emissions">
!<OVERVIEW>
!  A subroutine to calculate the internal sources and sinks of carbon dioxide
!  from input co2 emissions data.
!
! do_co2_emissions   = logical to activate using co2 emissions: default = .false.
!
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate the sources and sinks of carbon dixoide from co2 emissions.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_co2_emissions (Time, dt,  pwt, co2, sphum, co2_emiss_dt, kbot)
!
!</TEMPLATE>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model timestep.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav (kg/m2)
!   </IN>
!   <IN NAME="co2" TYPE="real" DIM="(:,:,:)">
!     The array of the carbon dioxide mixing ratio (kg co2/kg moist air)
!   </IN>
!   <IN NAME="sphum" TYPE="real" DIM="(:,:,:)">
!     The array of the specific humidity mixing ratio (kg/kg)
!   </IN>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

!   <OUT NAME="co2_emiss_dt" TYPE="real" DIM="(:,:,:)">
!     The array of the restoring tendency of the carbon dioxide emissions
!   </OUT>
!
!

subroutine atmos_co2_emissions(is, ie, js, je, Time, Time_next, dt, pwt, co2, sphum, co2_emiss_dt, kbot)

   integer, intent(in)                 :: is, ie, js, je
   type (time_type),      intent(in)   :: Time, Time_next
   real, intent(in)                    :: dt
   real, intent(in),  dimension(:,:,:) :: pwt          ! kg/m2
   real, intent(in),  dimension(:,:,:) :: co2          ! moist mmr
   real, intent(in),  dimension(:,:,:) :: sphum        
   real, intent(out), dimension(:,:,:) :: co2_emiss_dt 
   integer, intent(in),  dimension(:,:), optional :: kbot
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_emissions'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


integer   :: i,j,k,id,jd,kd,kb, logunit
logical   :: sent
real, dimension(size(co2,1),size(co2,2),size(co2,3)) ::  co2_emis_source
real, dimension(size(co2,1),size(co2,2)) :: co2_emis2d
logical   :: used


!---------------------------------------------------------------------------------------
! Original IPCC AR5 historical CO2 emissions converted to kg C/m2/sec 
! Assume scenario input will be in same format/units
!---------------------------------------------------------------------------------------

id=size(co2,1); jd=size(co2,2); kd=size(co2,3)

co2_emis2d(:,:)=0.0
co2_emis_source(:,:,:)=0.0
co2_emiss_dt(:,:,:)=0.0

logunit=stdlog()
if (ind_co2 > 0 .and. do_co2_emissions) then


  call data_override('ATM', 'co2_emiss', co2_emis2d, Time, override=used, &
                     is_in=is, ie_in=ie, js_in=js, je_in=je)
  if (.not. used) then
    call error_mesg (trim(error_header), ' data override needed for co2 emission ', FATAL)
  endif

  if (id_co2_emiss_orig > 0) sent = send_data (id_co2_emiss_orig, co2_emis2d, Time_next, is_in=is,js_in=js)

! lowest model layer
    do j=1,jd
      do i=1,id
        co2_emis_source(i,j,kd) = co2_emis2d(i,j) * (WTMCO2/WTMC) / pwt(i,j,kd)
      enddo
    enddo
  
  co2_emiss_dt = co2_emis_source

! co2 mol emission diagnostic in moles CO2/m2/sec 
  if (id_co2_mol_emiss > 0) sent = send_data (id_co2_mol_emiss,   &
                 co2_emiss_dt(:,:,kd)*pwt(:,:,kd)/(WTMCO2*1.e-3), Time_next, is_in=is,js_in=js)

endif


end subroutine atmos_co2_emissions
!</SUBROUTINE >


!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_gather_data">
!<OVERVIEW>
!  A subroutine to gather fields needed for calculating the CO2 gas flux
!</OVERVIEW>
!

subroutine atmos_co2_gather_data (gas_fields, tr_bot)

use coupler_types_mod, only: coupler_2d_bc_type, ind_pcair

implicit none

!-----------------------------------------------------------------------

type(coupler_2d_bc_type), intent(inout) :: gas_fields
real, dimension(:,:,:), intent(in)      :: tr_bot

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_gather_data'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

! 2008/06/17 JPD/jgj: OCMIP calculation expects pco2 in dry vmr (mol/mol) units 
! atm co2 is in moist mass mixing ratio (kg co2/kg moist air)
! tr_bot: co2 bottom layer moist mass mixing ratio
! convert to dry_mmr and then to dry_vmr for ocean model.
! dry_mmr = wet_mmr / (1-Q); co2vmr = (wair/wco2) * co2mmr

if (ind_co2_flux .gt. 0) then
  gas_fields%bc(ind_co2_flux)%field(ind_pcair)%values(:,:) = (tr_bot(:,:,ind_co2) / &
       (1.0 - tr_bot(:,:,ind_sphum))) * (WTMAIR/gas_fields%bc(ind_co2_flux)%mol_wt)
endif

end subroutine atmos_co2_gather_data
!</SUBROUTINE >

!#######################################################################


!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_flux_init">

!<OVERVIEW>
! Subroutine to initialize the carbon dioxide flux
!</OVERVIEW>

 subroutine atmos_co2_flux_init

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: n

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_flux_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: logunit

if ( .not. module_is_initialized) then

!----- set initial value of carbon ------------

  call tracer_manager_init      ! need to call here since the ocean pes never call it
  n = get_tracer_index(MODEL_ATMOS,'co2')
  if (n > 0) then
    ind_co2 = n
    if (ind_co2 > 0) then
      logunit=stdout()
      write (logunit,*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
      logunit=stdlog()
      write (logunit,*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
    endif
  endif
  module_is_initialized = .TRUE.
endif

!
!       initialize coupler flux
!

if (ind_co2 > 0) then
  ind_co2_flux = aof_set_coupler_flux('co2_flux',                       &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',       &
       atm_tr_index = ind_co2,                                          &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),              &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
endif

!-----------------------------------------------------------------------

end subroutine atmos_co2_flux_init
!</SUBROUTINE>


!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_init">

!<OVERVIEW>
! Subroutine to initialize the carbon dioxide module.
!</OVERVIEW>

 subroutine atmos_co2_init (Time, axes)

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

type(time_type),  intent(in)                        :: Time
integer, dimension(3), intent(in)                   :: axes

!
!-----------------------------------------------------------------------
!     local variables
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!-----------------------------------------------------------------------
!
integer :: ierr, unit, io, logunit
integer :: n
real    :: missing_value = -1.e10
character(len=64) :: desc
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


     if (module_is_initialized) return

     call write_version_number (version, tagname)

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=atmos_co2_nml, iostat=io)
        ierr = check_nml_error(io,'atmos_co2_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=atmos_co2_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'atmos_co2_nml')
        end do
10      call close_file (unit)
#endif
      endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=atmos_co2_nml)

!----- set initial value of carbon ------------

n = get_tracer_index(MODEL_ATMOS,'co2')
if (n > 0) then
  ind_co2 = n
    logunit=stdout()
    write (logunit,*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
    logunit=stdlog()
    write (logunit,*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2

 ! initialize diagnostics

   desc = ' restoring tendency'

   id_co2restore  = register_diag_field ('atmos_co2_restoring', 'co2_restore', axes, Time, &
                   'CO2'//trim(desc), 'moles co2/m2/s',missing_value=missing_value)

   desc = ' pressure weighting array = dP/grav'
   id_pwt    = register_diag_field ('atmos_co2', 'pwt', axes, Time, &
                   trim(desc), 'kg/m2',missing_value=missing_value)

   desc = ' mol emission'
   id_co2_mol_emiss = register_diag_field ('atmos_co2_emissions', 'co2_mol_emission', axes(1:2), Time, &
                      'CO2'//trim(desc), 'moles co2/m2/s',missing_value=missing_value)

   desc = ' emission_orig'
   id_co2_emiss_orig = register_diag_field ('atmos_co2_emissions', 'co2_emissions_orig', axes(1:2), Time, &
                   'CO2'//trim(desc), 'kg C/m2/s',missing_value=missing_value)

!
!       get the index for sphum
!

  ind_sphum = get_tracer_index(MODEL_ATMOS,'sphum')
  if (ind_sphum .le. 0) then
    call error_mesg (trim(error_header), ' Could not find index for sphum', FATAL)
  endif

endif

logunit=stdlog()
if (.not.(ind_co2 > 0 .and. do_co2_restore)) then
   if (mpp_pe() == mpp_root_pe() ) &
     write (logunit,*)' CO2 restoring not active:do_co2_restore= ',do_co2_restore
endif

if (.not.(ind_co2 > 0 .and. co2_radiation_override)) then
   if (mpp_pe() == mpp_root_pe() ) &
     write (logunit,*)' CO2 radiation override not active:co2_radiation_override= ',co2_radiation_override
endif

if (.not.(ind_co2 > 0 .and. do_co2_emissions)) then
   if (mpp_pe() == mpp_root_pe() ) &
     write (logunit,*)' not using CO2 emissions: do_co2_emissions= ',do_co2_emissions
endif

call write_version_number (version, tagname)
module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

end subroutine atmos_co2_init
!</SUBROUTINE>


!<SUBROUTINE NAME ="atmos_co2_end">
subroutine atmos_co2_end

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_end'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   module_is_initialized = .FALSE.

end subroutine atmos_co2_end
!</SUBROUTINE>


end module atmos_co2_mod
