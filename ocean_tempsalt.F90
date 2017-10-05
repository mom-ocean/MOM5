module ocean_tempsalt_mod
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> R.A.S. Fiedler
!</CONTACT>
!
!<CONTACT EMAIL="David.Jackett@csiro.au"> David Jackett 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Setup the temperature and salinity tracer fields.  
! Convert between potential temperature and conservative temperature.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This routine sets up the temperature and salinity tracer fields.  
! It also performs a converstion between potential temperature and 
! conservative temperature, with the conversion valid over most 
! large-scale oceanographically relevant ranges.  
!
! 0psu <= salinity <= 40 psu
!
! -3C <= theta <= 40C    (theta=conservative temperature or potential temperature) 
!
! 0dbar <= pressure <= 8000dbar
!
!  Input variables are the following:
!
!  salinity
!
!  potential temperature (theta) in deg C
!  OR
!  conservative temperature (theta) in deg C
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Feistel (2003)
! A new extended Gibbs thermodynamic potential of seawater
! Progress in Oceanography. vol 58, pages 43-114.
! </REFERENCE>
!
! <REFERENCE>
! Jackett, McDougall, Feistel, Wright, and Griffies (2005)
! Algorithms for density, potential temperature,
! conservative temperature, and freezing temperature of
! seawater. Journal of Atmospheric and Oceanic 
! Technology, 2005 submitted.
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_tempsalt_nml">
!
!  <DATA NAME="teos10"  TYPE="logical">
!  For choosing whether to use the TEOS-10 equation of state.
!  This usage requires conservative temperature as the temperature variable,
!  and two salinity variables: Preformed Salinity and Absolute Salinity
!  anomaly. Default teos10=.false. 
!  </DATA>
!
!  <DATA NAME="do_fafmip_heat"  TYPE="logical">
!  Set to true to initialize the CMIP6/FAFMIP heat fields.  
!  Default do_fafmip_heat=.false. 
!  </DATA>
!
!  <DATA NAME="temperature_variable" TYPE="character">
!  For choosing the temperature variable used in the model.
!  Choices are 'conservative_temp' and 'potential_temp'.
!  Since conservative temperature is more accurate, it is the default. 
!  </DATA>
!
!  <DATA NAME="pottemp_equal_contemp" TYPE="logical">
!  For certain idealized cases where the difference between potential
!  temperature and conservative temperature is irrelevant. Default=.false. 
!  </DATA>
!
!  <DATA NAME="pottemp_2nd_iteration" TYPE="logical">
!  For taking extra iteration in computation of potential temperature
!  from conservative temperature and salinity. Default is true.  
!  </DATA>
!
!  <DATA NAME="reinit_ts_with_ideal" TYPE="logical">
!  For setting up an ideal temperature and salinity profile
!  that is generated in the model. This profile can be 
!  generated after the model has already been running, hence
!  the name "reinit" for "reinitialize."
!  </DATA>
!  <DATA NAME="reinit_ts_with_ideal_efold" UNITS="metre" TYPE="real">
!  For setting efolding of reinitialized temp and salinity profile. 
!  Default reinit_ts_with_ideal_efold=1000.
!  </DATA> 
!  <DATA NAME="reinit_ts_with_ideal_tvalue" UNITS="C" TYPE="real">
!  For setting the reinitialized temperature value using the 
!  ideal profile. Default reinit_ts_with_ideal_tvalue = 10.0
!  </DATA> 
!  <DATA NAME="reinit_ts_with_ideal_svalue" UNITS="psu" TYPE="real">
!  For setting the reinitialized temperature value using the 
!  ideal profile. Default reinit_ts_with_ideal_svalue = 30.0
!  </DATA> 
!
!  <DATA NAME="t_min" UNITS="deg C" TYPE="real">
!  Minimum potential/conservative temperature below which we gracefully bring down the model. 
!  </DATA> 
!  <DATA NAME="t_max" UNITS="deg C" TYPE="real">
!  Maximum potential/conservative temperature above which we gracefully bring down the model.  
!  </DATA> 
!  <DATA NAME="s_min" UNITS="ppt" TYPE="real">
!  Minimum salinity below which we gracefully bring down the model. 
!  </DATA> 
!  <DATA NAME="s_max" UNITS="ppt" TYPE="real">
!  Maximum salinity below which we gracefully bring down the model. 
!  </DATA> 
!
!  <DATA NAME="t_min_limit" UNITS="deg C" TYPE="real">
!  Minimum potential/conservative temperature below which will employ upwind advection 
!  instead of quicker, and horizontal diffusion instead of neutral physics. 
!  </DATA> 
!  <DATA NAME="t_max_limit" UNITS="deg C" TYPE="real">
!  Maximum potential/conservative temperature above which will employ upwind advection 
!  instead of quicker, and horizontal diffusion instead of neutral physics. 
!  </DATA> 
!  <DATA NAME="s_min_limit" UNITS="psu" TYPE="real">
!  Minimum salinity below which will employ upwind advection instead
!  of quicker, and horizontal diffusion instead of neutral physics.  
!  </DATA> 
!  <DATA NAME="s_max_limit" UNITS="psu" TYPE="real">
!  Maximum salinity below which will employ upwind advection instead
!  of quicker, and horizontal diffusion instead of neutral physics.  
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging the module.
!  </DATA>
!
!</NAMELIST>

use constants_mod,        only: kelvin
use fms_mod,              only: write_version_number, mpp_error, FATAL, WARNING
use fms_mod,              only: open_namelist_file, check_nml_error, close_file
use mpp_mod,              only: input_nml_file, stdout, stdlog
use mpp_mod,              only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,              only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: cp_ocean, CP_OCEAN_PRETEOS10, CP_OCEAN_TEOS10 
use ocean_tpm_util_mod,   only: otpm_set_prog_tracer, otpm_set_diag_tracer
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_time_type, ocean_options_type, ocean_thickness_type
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3

implicit none

private

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

character(len=128) :: version=&
       '$Id: ocean_tempsalt.F90,v 20.0 2013/12/14 00:17:14 fms Exp $'
character (len=128) :: tagname = &
       '$Name: tikal $'

private dentropy_dtheta
interface dentropy_dtheta
  module procedure dentropy_dtheta_field
  module procedure dentropy_dtheta_level
  module procedure dentropy_dtheta_point
end interface 

public ocean_tempsalt_init
public ocean_tempsalt_ideal_reinit
public contemp_from_pottemp
public pottemp_from_contemp
public tempsalt_check_range

! interfaces for contemp_from_pottemp
interface contemp_from_pottemp
   module procedure contemp_from_pottemp_field
   module procedure contemp_from_pottemp_level
   module procedure contemp_from_pottemp_point
end interface

! interfaces for pottemp_from_contemp
interface pottemp_from_contemp
   module procedure pottemp_from_contemp_field
   module procedure pottemp_from_contemp_level
   module procedure pottemp_from_contemp_point
end interface

! useful constants 
real,parameter :: one_over_40=2.5e-2   
real,parameter :: one_over_1600=6.25e-4   
real,parameter :: sfact=35./35.16504

! for diagnostics 
integer :: id_contemp=-1
integer :: id_pottemp=-1
logical :: used 

! for timing clocks 
integer :: id_clock_c2p_1
integer :: id_clock_c2p_2
integer :: id_clock_c2p_3

! preTEOS10 coefficients for computing contemp to pottemp
real ::  a0,  a1,  a2,  a3,  a4,  a5    
real ::  b0,  b1,  b2,  b3              

! preTEOS10 coefficients for computing pottemp to contemp
real ::  c0,  c1,  c2,  c3,  c4,  c5    
real ::  c6,  c7,  c8,  c9, c10, c11    
real :: c12, c13, c14, c15, c16, c17    
real :: c18, c19, c20, c21, c22, c23    

! TEOS10 coefficients for computing contemp from pottemp 
real ::  sfac
real ::  v0,  v1,  v2,  v3,  v4,  v5    
real ::  v6,  v7,  v8,  v9, v10, v11    
real :: v12, v13, v14, v15, v16, v17    
real :: v18, v19, v20, v21, v22, v23    
real :: v24, v25, v26, v27, v28, v29    

! preTEOS10 coefficients for computing d(entropy)/d(pottemp)
real ::  f0,  f1,  f2,  f3,  f4,  f5, f6    
real ::  f7,  f8,  f9, f10, f11, f12, f13    

! TEOS10 polynomial coefficients for computing d(entropy)/d(pottemp) 
real ::  g1,  g2,  g3,  g4,  g5, g6    
real ::  h1,  h2,  h3,  h4,  h5, h6    
real ::  h7,  h8,  h9, h10, h11, h12

! for testing the conversion between pottemp and contemp 
real :: t_test= 20.0               ! deg C
real :: s_test= 20.0               ! psu
real :: ct_ref= 20.4527496128276   ! deg C
real :: pt_ref= 19.55627910604363  ! deg C 

! conversion factor for salinity  (preformed -> absolute)
real :: Fdelta

logical :: module_is_initialized = .FALSE.

logical :: reinit_ts_with_ideal        = .false. 
real    :: reinit_ts_with_ideal_efold  = 1000.0 
real    :: reinit_ts_with_ideal_tvalue = 10.0
real    :: reinit_ts_with_ideal_svalue = 30.0
 
! set min/max temp and salinity range, outside 
! of which the model will be brought down.
real :: t_min=-5.0
real :: t_max=40.0
real :: s_min=0.0
real :: s_max=70.0

! min/max temp and salinity beyond which employ upwind 
! advection if set limit_tracer_range_quick=.true.
! (set in ocean_tracer_advect_nml) and/or 
! horizontal diffusion if neutral_physics_limit=.true. 
! (set in ocean_nphysicsA_nml, ocean_nphysicsB_nml, or ocean_nphysicsC_nml)
! or remove submesoscale fluxes if submeso_flux_limit=.true.
! (set in ocean_submesoscale_nml)
real :: t_min_limit=-2.0
real :: t_max_limit=32.0
real :: s_min_limit=1.0
real :: s_max_limit=42.0

logical :: debug_this_module     = .false.
logical :: pottemp_2nd_iteration = .true.    
logical :: pottemp_equal_contemp = .false.
logical :: teos10                = .false.
logical :: do_fafmip_heat        = .false.

character(len=32) :: temperature_variable='conservative_temp'  ! "conservative_temp" or "potential_temp"          

namelist /ocean_tempsalt_nml/ debug_this_module, temperature_variable, pottemp_2nd_iteration,     &
                              pottemp_equal_contemp,                                              &
                              t_min,       t_max,       s_min,       s_max,                       &
                              t_min_limit, t_max_limit, s_min_limit, s_max_limit,                 &
                              reinit_ts_with_ideal, reinit_ts_with_ideal_efold,                   &
                              reinit_ts_with_ideal_tvalue, reinit_ts_with_ideal_svalue,           &
                              teos10, do_fafmip_heat 
                
contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tempsalt_init">
!
! <DESCRIPTION>
! Initialize the temperature/salinity module.
! </DESCRIPTION>
!
subroutine ocean_tempsalt_init(Domain, Grid, Ocean_options, itemp, isalt, i_redist_heat, debug)

  type(ocean_domain_type),    intent(in), target   :: Domain
  type(ocean_grid_type),      intent(in), target   :: Grid
  type(ocean_options_type),   intent(inout)        :: Ocean_options
  integer,                    intent(inout)        :: itemp
  integer,                    intent(inout)        :: isalt
  integer,                    intent(inout)        :: i_redist_heat
  logical,                    intent(in), optional :: debug

  integer           :: ioun, io_status, ierr
  integer           :: index_temp       =-1
  integer           :: index_salt       =-1  
  integer           :: index_diag_temp  =-1
  integer           :: index_Fdelta     =-1
  integer           :: index_added_heat =-1
  integer           :: index_redist_heat=-1
  character(len=32) :: longname_temperature 
  real              :: test_value, diff_test

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (ocean_tempsalt_init): module initialized.')
  endif

  module_is_initialized = .TRUE.

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  Dom => Domain
  Grd => Grid

  call write_version_number(version, tagname)

  ! provide for namelist override of defaults
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_tempsalt_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_tempsalt_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_tempsalt_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_tempsalt_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_tempsalt_nml)
  write (stdlogunit,ocean_tempsalt_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
  endif
  if(debug_this_module) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_tempsalt with debug_this_module=.true.'  
  endif

  ! check for consistent use of heat capacity according to teos10 or preteos10 
  if(teos10) then 
      if(cp_ocean==CP_OCEAN_PRETEOS10) then 
          call mpp_error(FATAL,  '==>Error in ocean_tempsalt_mod: Must have' //  &
               ' cp_ocean = CP_OCEAN_TEOS10 with TEOS10 equation of state.' //   &
               ' Make code change inside ocean_parameters.F90 and recompile.')
      endif
  else
      if(cp_ocean==CP_OCEAN_TEOS10) then 
          call mpp_error(WARNING,  '==>Warning in ocean_tempsalt_mod: Using ' //          &
               ' cp_ocean = cp_ocean_teos10, yet not using TEOS10 equation of state.' //  &
               ' You may use earlier CP_OCEAN_PRETEOS10 value by changing ocean_parameters.F90')
      endif
  endif


  if(temperature_variable/='conservative_temp' .and. teos10 ) then
     call mpp_error(FATAL,  '==>Error in ocean_tempsalt_mod: Must have' //  &
     ' prognostic temperature variable = conservative_temp when using TEOS10')
  endif
  if(temperature_variable=='conservative_temp') then
     write(stdoutunit,'(a)') '==>Note from ocean_tempsalt_mod: MOM prognostic temp = conservative temperature.' 
     write(stdoutunit,'(a/)')'                                 MOM diagnostic temp = potential temperature.' 
     longname_temperature = 'Conservative temperature'
     Ocean_options%temperature_variable = 'Prognostic temperature variable is conservative temperature.'
  elseif(temperature_variable=='potential_temp') then 
     write(stdoutunit,'(a)') '==>Note from ocean_tempsalt_mod: MOM prognostic temp = potential temperature.' 
     write(stdoutunit,'(a/)')'                                 MOM diagnostic temp = conservative temperature.' 
     longname_temperature = 'Potential temperature'
     Ocean_options%temperature_variable = 'Prognostic temperature variable is potential temperature.'
  else 
     write(stdoutunit,'(/a/)')'==>Error in ocean_tempsalt_mod: prognostic temperature variable remains unchosen.' 
     call mpp_error(FATAL,  '==>Error in ocean_tempsalt_mod: prognostic temperature variable remains unchosen.')
  endif 

  if(pottemp_equal_contemp) then 
     write(stdoutunit,'(a)') '==>Note from ocean_tempsalt_mod: setting potential temp = conservative temp.' 
     write(stdoutunit,'(a)') '   This setting is useful for simulations where the difference is irrelevant.'
  endif 


  ! initialize temperature and salinity

  index_temp = otpm_set_prog_tracer('temp', 'required',                                    &
       caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                    &
       longname=longname_temperature, units='deg_C',                                       &
       conversion=cp_ocean, offset=kelvin, min_tracer=t_min, max_tracer=t_max,             &
       min_range=-10.0, max_range=500.0, flux_units='Watts/m^2', min_flux_range=-1.0e+16,  &
       max_flux_range=1.0e+16, min_tracer_limit=t_min_limit, max_tracer_limit=t_max_limit, &
       psom_limit=.false.,ppm_hlimiter=2,ppm_vlimiter=2,mdt_scheme=1,                      &
       restart_file='ocean_temp_salt.res.nc' )
  itemp = index_temp

  if ( teos10 ) then
  index_salt = otpm_set_prog_tracer('salt', 'required',                                      &
       caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                      &
       longname='Preformed Salinity', units='g/kg',                                          &
       conversion=0.001, min_tracer=s_min, max_tracer=s_max,                                 &
       min_range=-10.0, max_range=100.0, flux_units='kg/(sec*m^2)', min_flux_range=-1.0e+05, &
       max_flux_range=1.0e+05, min_tracer_limit=s_min_limit, max_tracer_limit=s_max_limit,   &
       psom_limit=.false.,ppm_hlimiter=2,ppm_vlimiter=2,mdt_scheme=1,                        &
       restart_file='ocean_temp_salt.res.nc' )
  isalt = index_salt 
  index_Fdelta = otpm_set_prog_tracer('fdelta', 'required',                                &
       caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                    &
       longname='Salinity conversion factor', units='none',                                &
       conversion=1.0, min_tracer=s_min, max_tracer=s_max,                                 &
       min_range=-1.0, max_range=1.0, flux_units='kg/(sec*m^2)', min_flux_range=-1.0e+04,  &
       max_flux_range=1.0e+04, min_tracer_limit=s_min_limit, max_tracer_limit=s_max_limit, &
       psom_limit=.false.,ppm_hlimiter=2,ppm_vlimiter=2,mdt_scheme=1,                      &
       restart_file='ocean_Fdelta.res.nc' )
  else
  index_salt = otpm_set_prog_tracer('salt', 'required',                                      &
       caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                      &
       longname='Practical Salinity', units='psu',                                           &
       conversion=0.001, min_tracer=s_min, max_tracer=s_max,                                 &
       min_range=-10.0, max_range=100.0, flux_units='kg/(sec*m^2)', min_flux_range=-1.0e+05, &
       max_flux_range=1.0e+05, min_tracer_limit=s_min_limit, max_tracer_limit=s_max_limit,   &
       psom_limit=.false.,ppm_hlimiter=2,ppm_vlimiter=2,mdt_scheme=1,                        &
       restart_file='ocean_temp_salt.res.nc' )
  isalt = index_salt 
  endif

  if(temperature_variable=='conservative_temp') then
      index_diag_temp = otpm_set_diag_tracer('pot_temp',                                     &
           caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                  &
           longname='Potential temperature', units='deg_C',                                  &
           conversion=cp_ocean, offset=kelvin, min_tracer=t_min, max_tracer=t_max,           &
           min_range=-10.0, max_range=500.0, const_init_tracer=.true.,const_init_value=0.0,  &
           restart_file='ocean_pot_temp.res.nc' )
  elseif(temperature_variable=='potential_temp') then
      index_diag_temp = otpm_set_diag_tracer('con_temp',                                     &
           caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                  &
           longname='Conservative temperature', units='deg_C',                               &
           conversion=cp_ocean, offset=kelvin, min_tracer=t_min, max_tracer=t_max,           &
           min_range=-10.0, max_range=500.0, const_init_tracer=.true.,const_init_value=0.0,  &
           restart_file='ocean_con_temp.res.nc' )
  endif

  i_redist_heat = -1   
  if(do_fafmip_heat) then
    Ocean_options%fafmip_heat = 'Included the two CMIP6/FAFMIP prognostic temperature tracers.'
    write(stdoutunit,'(a)') '==>Note from ocean_tempsalt_mod: Initializing FAFMIP temperature fields.' 

    index_added_heat = otpm_set_prog_tracer('added_heat', 'required',                        &
         caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                    &
         longname='added_heat', units='deg_C',                                               &
         conversion=cp_ocean, offset=kelvin, min_tracer=t_min, max_tracer=t_max,             &
         min_range=-10.0, max_range=500.0, flux_units='Watts/m^2', min_flux_range=-1.0e+16,  &
         max_flux_range=1.0e+16, min_tracer_limit=t_min_limit, max_tracer_limit=t_max_limit, &
         psom_limit=.false.,ppm_hlimiter=2,ppm_vlimiter=2,mdt_scheme=1,                      &
         restart_file='ocean_added_heat.res.nc' )

    index_redist_heat = otpm_set_prog_tracer('redist_heat', 'required',                      &
         caller='ocean_tempsalt_mod/ocean_tempsalt_init',                                    &
         longname='redist_heat', units='deg_C',                                              &
         conversion=cp_ocean, offset=kelvin, min_tracer=t_min, max_tracer=t_max,             &
         min_range=-10.0, max_range=500.0, flux_units='Watts/m^2', min_flux_range=-1.0e+16,  &
         max_flux_range=1.0e+16, min_tracer_limit=t_min_limit, max_tracer_limit=t_max_limit, &
         psom_limit=.false.,ppm_hlimiter=2,ppm_vlimiter=2,mdt_scheme=1,                      &
         restart_file='ocean_redist_heat.res.nc' )
    i_redist_heat = index_redist_heat

  else
    Ocean_options%fafmip_heat = 'Did NOT include CMIP6/FAFMIP prognostic temperature tracers.'
    write(stdoutunit,'(a)') '==>Note from ocean_tempsalt_mod: NOT Initializing FAFMIP temperature fields.' 
  endif 


  if(reinit_ts_with_ideal) then 
      write(stdoutunit,'(a)') '==>Note: reinitializing ocean_tempsalt with ideal internally generated profile.'
  endif 

  ! polynomial coefficients for conservative temp and potential temp 
  
  a0 = -1.446013646344788e-2
  a1 = -3.305308995852924e-3
  a2 =  1.062415929128982e-4
  a3 =  9.477566673794488e-1
  a4 =  2.166591947736613e-3
  a5 =  3.828842955039902e-3

  b0 =  1.000000000000000e+0
  b1 =  6.506097115635800e-4
  b2 =  3.830289486850898e-3
  b3 =  1.247811760368034e-6


  if (teos10) then
     sfac = 9.95306702338459d-01
  else 
     sfac = 1.0
  endif 


  ! TEOS10 coefficients for CT from PT and salinity 
  v0  =  2.50509288068125d-04    ! v0=1/cp_ocean 
  v1 =  61.01362420681071d0 
  v2 = 168776.46138048015d0 
  v3 = -2735.2785605119625d0 
  v4 = 2574.2164453821433d0
  v5 = -1536.6644434977543d0 
  v6 = 545.7340497931629d0
  v7 = -50.91091728474331d0 
  v8 = -18.30489878927802d0

  v9 = 268.5520265845071d0 
  v10 = -12019.028203559312d0
  v11 = 3734.858026725145d0 
  v12 = -2046.7671145057618d0
  v13 = 465.28655623826234d0 
  v14 = -0.6370820302376359d0 
  v15 = -10.650848542359153d0
  v16 = 937.2099110620707d0 
  v17 = 588.1802812170108d0
  v18 = 248.39476522971285d0 
  v19 = -3.871557904936333d0
  v20 = -2.6268019854268356d0
  v21 = -1687.914374187449d0 
  v22 = 246.9598888781377d0
  v23 = 123.59576582457964d0 
  v24 = -48.5891069025409d0
  v25 = 936.3206544460336d0
  v26 = -942.7827304544439d0 
  v27 = 369.4389437509002d0
  v28 = -33.83664947895248d0 
  v29 = -9.987880382780322d0

  ! preTEOS10 coefficients for CT from PT and salinity
  c0  =  2.50494524832013e-4    ! c0=1/cp_ocean 
  c1  =  61.013624165232955
  c2  =  168776.46138048015
  c3  = -2735.2785605119643
  c4  =  2574.2164453821442
  c5  = -1536.6644434977545
  c6  =  545.734049793163
  c7  = -50.910917284743334                
  c8  = -18.30489878927802                
  c9  =  416.31512917743896
  c10 =  937.9793807560891                
  c11 = -3140.435779506947
  c12 =  2975.170149976973
  c13 = -1760.137081144729
  c14 =  414.5655751783703               
  c15 =  2167.72082596016               
  c16 = -1224.5772800562902
  c17 =  326.3074029273967
  c18 =  50.6703824689518                
  c19 = -12694.10018182362
  c20 =  4405.71847182968
  c21 = -2132.9690185026416
  c22 =  303.91071982808035               
  c23 =  69.74975368852  


  ! TEOS10 coefficients for computing d(entropy)/d(pottemp) 
  g1 =  -24715.571866078
  g2 = 4420.4472249096725
  g3 = -1778.231237203896
  g4 = 1160.5182516851419
  g5 = -569.531539542516 
  g6 = 128.13429152494615

  h1 = 1760.062705994408 
  h2 = -86.1329351956084
  h3 = -137.1145018408982 
  h4 = 296.20061691375236
  h5 = -205.67709290374563 
  h6 = 49.9394019139016
  h7 = -60.136422517125 
  h8 = 10.50720794170734
  h9 = -1351.605895580406 
  h10 = 1097.1125373015109
  h11 = -433.20648175062206 
  h12 = 63.905091254154904

  ! preTEOS10 coefficients for computing d(entropy)/d(pottemp) 
  f0  =  24715.571866078
  f1  = -1858.920033948178
  f2  =  317.440355256842
  f3  = -405.1392883572282
  f4  =  202.6815298758072
  f5  =  1562.563716288858
  f6  = -1165.8752731900836
  f7  =  348.7487684426
  f8  = -4420.4472249096725
  f9  =  1778.231237203896
  f10 = -1160.5182516851419
  f11 =  569.531539542516
  f12 =  128.13429152494615
  f13 =  6.2500000000000000e-4

 
  if(debug_this_module .and. .not. pottemp_equal_contemp) then 
      write (stdoutunit,'(/a/)') 'TEMPERATURE TEST VALUES'
      if(teos10) then
         ct_ref = 2.04566805754950e+01
         pt_ref = 1.95524819680081e+01
      endif
      test_value = contemp_from_pottemp(20.0,20.0)
      diff_test  = abs(test_value-ct_ref)
      write (stdoutunit,'(a,f6.2,a,f6.2,/a,e22.16,/a,e22.16,/a,e22.16)') &
           's_test(psu) = ',s_test,', pottemp_test(C) = ',t_test,      &
           'contemp_from_pottemp(C)       = ',test_value,              &
           'reference value               = ',ct_ref,                  &
           'abs difference between values = ',diff_test

      test_value = pottemp_from_contemp(20.0,20.0)
      diff_test  = abs(test_value-pt_ref)
      write (stdoutunit,'(/a,f6.2,a,f6.2,/a,e22.16,/a,e22.16,/a,e22.16)') &
           's_test(psu) = ',s_test,', contemp_test(C) = ',t_test,       &
           'pottemp_from_contemp(C)       = ',test_value,               &
           'reference value               = ',pt_ref,                   &
           'abs difference between values = ',diff_test
  endif

  ! clocks 
  id_clock_c2p_1  = mpp_clock_id('(Ocean CT2PT: 1) '  ,grain=CLOCK_MODULE)
  id_clock_c2p_2  = mpp_clock_id('(Ocean CT2PT: 2) '  ,grain=CLOCK_MODULE)
  id_clock_c2p_3  = mpp_clock_id('(Ocean CT2PT: 3) '  ,grain=CLOCK_MODULE)

  
end subroutine ocean_tempsalt_init
! </SUBROUTINE> NAME="ocean_tempsalt_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_tempsalt_ideal_reinit">
!
! <DESCRIPTION>
! Reinitialize the temperature/salinity fields with internally generated
! idealized profile.
!
! This routine can be called at any time within a model integration.  
! The main application is to allow for there to be a nontrivial 
! Thickness established, and then to reinitialize the temperature 
! and salinity to be functions of the depth_st array.  When 
! s=zstar, this approach will produce a flat initialization in 
! zstar space, but nontrivial initialization in depth space.  
!
! User can customize a profile.  
!
! One example given here assumes profile is a function of depth_st,
! with the assumption that s=geopotential OR s=zstar.  
! To use this approach with pressure coordinate, need to alter the 
! reinit_ts_with_ideal_efold value, since depth_st is in MKS, not dbar. 
!
! </DESCRIPTION>
!
subroutine ocean_tempsalt_ideal_reinit(Thickness, index_temp, index_salt, T_prog)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  integer,                      intent(in)    :: index_temp
  integer,                      intent(in)    :: index_salt
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: i,j,k,n

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. reinit_ts_with_ideal) return 

  write(stdoutunit,'(/a/)') &
  '==>From ocean_tempsalt_ideal_reinit: reinitializing temp/salinity with ideal profile.'

  do n=1,3
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(index_temp)%field(i,j,k,n) = Grd%tmask(i,j,k) &
              *reinit_ts_with_ideal_tvalue*exp(-Thickness%depth_st(i,j,k)/reinit_ts_with_ideal_efold)   
              T_prog(index_salt)%field(i,j,k,n) = Grd%tmask(i,j,k) & 
              *reinit_ts_with_ideal_svalue*exp(-Thickness%depth_st(i,j,k)/reinit_ts_with_ideal_efold)    
           enddo
        enddo
     enddo
  enddo


end subroutine ocean_tempsalt_ideal_reinit
! </SUBROUTINE> NAME="ocean_tempsalt_ideal_reinit"



!#######################################################################
! <FUNCTION NAME="contemp_from_pottemp_field">
!
! <DESCRIPTION>
!
! Compute conservative temperature for all grid points. 
! Input is potential temperature (C) and salinity(psu or g/kg).
!
! contemp = potential_enthalpy(s,theta)/cp_ocean
!
! </DESCRIPTION>
!
function contemp_from_pottemp_field (salinity, theta)

  real, dimension(isd:,jsd:,:), intent(in) :: salinity
  real, dimension(isd:,jsd:,:), intent(in) :: theta
  real, dimension(isd:ied,jsd:jed,nk)      :: contemp_from_pottemp_field

  integer :: i, j, k
  real    :: t1, s1, sp5

  if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (contemp_from_pottemp_field): module not initialized')
  endif

  if(pottemp_equal_contemp) then  
     contemp_from_pottemp_field = theta 
     return 
  endif 

  if (teos10) then

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied

               t1  = one_over_40*theta(i,j,k)
               s1  = one_over_40*sfac*salinity(i,j,k)
               sp5 = sqrt(s1)

               contemp_from_pottemp_field(i,j,k) =                                  &
                    v0*(v1+t1*(v2+t1*(v3+t1*(v4+t1*(v5+t1*(v6+t1*(v7+v8*t1))))))  + &
                    s1*(v9+t1*(v10+t1*(v11+t1*(v12+t1*(v13+t1*(v14+t1*v15)))))    + &
                    sp5*(v16+t1*(v17+t1*(v18+t1*(v19+t1*v20)))                    + &
                    sp5*(v21+sp5*(v22+sp5*(v23+sp5*v24))                          + &
                    t1*(v25+t1*(v26+t1*(v27+t1*(v28+t1*v29))))))))


            enddo
         enddo
      enddo

  else   ! preTEOS10 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied

               t1  = one_over_40*theta(i,j,k)
               s1  = one_over_40*salinity(i,j,k)
               sp5 = sqrt(s1)                                       

               contemp_from_pottemp_field(i,j,k) =                                  &
                    c0*(c1+t1*(c2+t1*(c3+t1*(c4+t1*(c5+t1*(c6+(c7+c8*t1)*t1)))))  + &  
                    s1*(c9+sp5*(c10+sp5*(c11+sp5*(c12+sp5*(c13+sp5*c14)))+t1*(c15 + &
                    t1*(c16+t1*(c17+c18*t1))))+t1*(c19+t1*(c20+t1*(c21+t1*(c22+c23*t1))))))

            enddo
         enddo
      enddo

  endif

end function contemp_from_pottemp_field
! </FUNCTION> NAME="contemp_from_pottemp_field"


!#######################################################################
! <FUNCTION NAME="contemp_from_pottemp_level">
!
! <DESCRIPTION>
!
! Compute conservative temperature for one k-level.
! Input is potential temperature (C) and salinity(psu or g/kg).
!
! contemp = potential_enthalpy(s,theta)/cp_ocean
!
! </DESCRIPTION>
!
function contemp_from_pottemp_level(salinity, theta)

  real, dimension(isd:,jsd:), intent(in) :: salinity
  real, dimension(isd:,jsd:), intent(in) :: theta
  real, dimension(isd:ied,jsd:jed)       :: contemp_from_pottemp_level

  integer :: i, j
  real    :: t1, s1, sp5

  if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (contemp_from_pottemp_level): module not initialized')
  endif

  if(pottemp_equal_contemp) then  
     contemp_from_pottemp_level = theta 
     return 
  endif 

  if (teos10) then

      do j=jsd,jed
         do i=isd,ied

            t1  = one_over_40*theta(i,j)
            s1  = one_over_40*sfac*salinity(i,j)
            sp5 = sqrt(s1)

            contemp_from_pottemp_level(i,j) =                                    &
                 v0*(v1+t1*(v2+t1*(v3+t1*(v4+t1*(v5+t1*(v6+t1*(v7+v8*t1))))))  + &
                 s1*(v9+t1*(v10+t1*(v11+t1*(v12+t1*(v13+t1*(v14+t1*v15)))))    + &
                 sp5*(v16+t1*(v17+t1*(v18+t1*(v19+t1*v20)))                    + &
                 sp5*(v21+sp5*(v22+sp5*(v23+sp5*v24))                          + &
                 t1*(v25+t1*(v26+t1*(v27+t1*(v28+t1*v29))))))))

         enddo
      enddo

  else

      do j=jsd,jed
         do i=isd,ied

            t1  = one_over_40*theta(i,j)
            s1  = one_over_40*salinity(i,j)
            sp5 = sqrt(s1)

            contemp_from_pottemp_level(i,j) =                                      &
                 c0*(c1+t1*(c2+t1*(c3+t1*(c4+t1*(c5+t1*(c6+(c7+c8*t1)*t1)))))    + &  
                 s1*(c9+sp5*(c10+sp5*(c11+sp5*(c12+sp5*(c13+sp5*c14)))+t1*(c15   + &
                 t1*(c16+t1*(c17+c18*t1))))+t1*(c19+t1*(c20+t1*(c21+t1*(c22+c23*t1))))))


         enddo
      enddo

  endif

end function contemp_from_pottemp_level
! </FUNCTION> NAME="contemp_from_pottemp_level"


!#######################################################################
! <FUNCTION NAME="contemp_from_pottemp_point">
!
! <DESCRIPTION>
!
! Compute conservative temperature for one grid point.
! Input is potential temperature (C) and salinity(psu or g/kg).
!
! contemp = potential_enthalpy(s,theta)/cp_ocean
!
! </DESCRIPTION>
!
function contemp_from_pottemp_point(salinity, theta)

  real, intent(in) :: salinity
  real, intent(in) :: theta
  real             :: contemp_from_pottemp_point

  real :: t1,s1,sp5

  if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (contemp_from_pottemp_point): module not initialized')
  endif

  if(pottemp_equal_contemp) then  
     contemp_from_pottemp_point = theta 
     return 
  endif 

  t1  = one_over_40*theta
  s1  = one_over_40*sfac*salinity
  sp5 = sqrt(s1)

  if (teos10) then

     contemp_from_pottemp_point =                                          &
          v0*(v1+t1*(v2+t1*(v3+t1*(v4+t1*(v5+t1*(v6+t1*(v7+v8*t1))))))  +  &
              s1*(v9+t1*(v10+t1*(v11+t1*(v12+t1*(v13+t1*(v14+t1*v15)))))+  &
                  sp5*(v16+t1*(v17+t1*(v18+t1*(v19+t1*v20)))            +  &
                       sp5*(v21+sp5*(v22+sp5*(v23+sp5*v24))             +  &
                       t1*(v25+t1*(v26+t1*(v27+t1*(v28+t1*v29))))))))
  else

     contemp_from_pottemp_point =                                          &
          c0*(c1+t1*(c2+t1*(c3+t1*(c4+t1*(c5+t1*(c6+(c7+c8*t1)*t1)))))   + &  
          s1*(c9+sp5*(c10+sp5*(c11+sp5*(c12+sp5*(c13+sp5*c14)))+t1*(c15  + &
          t1*(c16+t1*(c17+c18*t1))))+t1*(c19+t1*(c20+t1*(c21+t1*(c22+c23*t1))))))
  endif


end function contemp_from_pottemp_point
! </FUNCTION> NAME="contemp_from_pottemp_point"


!#######################################################################
! <FUNCTION NAME="pottemp_from_contemp_field">
!
! <DESCRIPTION>
!
! Compute potential temperature from conservative temperature 
! for all grid points. Perform one extra iteration to get 
! precision to near computer precision.  
!
! Input is salinity (psu) and conservative temperature (C).
!
! Use wrk1, wrk2, and wrk3 so to not take 3-d arrays from stack. 
! 
! </DESCRIPTION>
!
function pottemp_from_contemp_field (salinity, ct)

  real, dimension(isd:,jsd:,:), intent(in) :: salinity
  real, dimension(isd:,jsd:,:), intent(in) :: ct
  real, dimension(isd:ied,jsd:jed,nk)      :: pottemp_from_contemp_field

  integer :: i, j, k
  real    :: a5ct, b3ct, ct_factor
  real    :: th0_num, rec_th0_den

  if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (pottemp_from_contemp_field): module not initialized')
  endif

  if(pottemp_equal_contemp) then
     pottemp_from_contemp_field = ct
     return 
  endif 

  call mpp_clock_begin(id_clock_c2p_1)
  ! first estimate (let wrk1=th0 and wrk2=dth_dct) 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           a5ct        = a5*ct(i,j,k)
           b3ct        = b3*ct(i,j,k)
           ct_factor   = a3 + a4*salinity(i,j,k) + a5ct
           th0_num     = a0 + salinity(i,j,k)*(a1+a2*salinity(i,j,k)) + ct(i,j,k)*ct_factor
           rec_th0_den = 1.0 / ( b0 + b1*salinity(i,j,k) + ct(i,j,k)*(b2+b3ct) )
           wrk1(i,j,k) = th0_num*rec_th0_den
           wrk2(i,j,k) = ( ct_factor + a5ct - (b2+b3ct+b3ct)*wrk1(i,j,k) )*rec_th0_den        
        enddo
     enddo
  enddo
  call mpp_clock_end(id_clock_c2p_1)

  if(pottemp_2nd_iteration) then 

      ! TEOS-10 (via McDougall & Barker)uses a modified NR scheme. 
      ! The first iteration only adds half the increment 
      ! but then calculates a derivative here before performing a full update. 
      if (teos10) then

         wrk3                       = contemp_from_pottemp(salinity,wrk1) - ct  !CT_diff
         pottemp_from_contemp_field = wrk1 - 0.5* wrk3*wrk2                     ! ptm
         wrk2 = cp_ocean/((pottemp_from_contemp_field + kelvin )*dentropy_dtheta(salinity,pottemp_from_contemp_field)) !1/dCT_dpt

         pottemp_from_contemp_field = wrk1 - wrk3*wrk2                          ! pt_tmp end of first iteration.
         pottemp_from_contemp_field = pottemp_from_contemp_field -              &
             (contemp_from_pottemp(salinity,pottemp_from_contemp_field) - ct)*wrk2
      else

      ! second estimate to get precision to roundoff (let wrk3=ct0)
      wrk3                       = contemp_from_pottemp(salinity,wrk1)  
      pottemp_from_contemp_field = wrk1 - (wrk3-ct)*wrk2
      wrk2                       = cp_ocean &
           /( (pottemp_from_contemp_field + kelvin) * dentropy_dtheta(salinity,pottemp_from_contemp_field) )
      pottemp_from_contemp_field =  pottemp_from_contemp_field  &
           - (contemp_from_pottemp(salinity,pottemp_from_contemp_field) - ct)*wrk2
      endif

  else 

      pottemp_from_contemp_field = wrk1

  endif


end function pottemp_from_contemp_field
! </FUNCTION> NAME="pottemp_from_contemp_field"


!#######################################################################
! <FUNCTION NAME="pottemp_from_contemp_level">
!
! <DESCRIPTION>
!
! Compute potential temperature from conservative temperature 
! over a k-level. Perform one extra iteration to get 
! precision to near computer precision.
!
! Input is salinity (psu) and conservative temperature (C).
!
! </DESCRIPTION>
!
function pottemp_from_contemp_level(salinity, ct)

  real, dimension(isd:,jsd:), intent(in) :: salinity
  real, dimension(isd:,jsd:), intent(in) :: ct
  real, dimension(isd:ied,jsd:jed)       :: pottemp_from_contemp_level

  real,dimension(isd:ied,jsd:jed) :: th0
  real,dimension(isd:ied,jsd:jed) :: ct0
  real,dimension(isd:ied,jsd:jed) :: dth_dct
  real,dimension(isd:ied,jsd:jed) :: pt_tmp

  integer :: i, j
  real    :: a5ct, b3ct,ct_factor
  real    :: th0_num, rec_th0_den

  if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (ct_from_pt_level): module must be initialized')
  endif

  if(pottemp_equal_contemp) then
     pottemp_from_contemp_level = ct
     return 
  endif 

  do j=jsd,jed
     do i=isd,ied
        a5ct         = a5*ct(i,j)
        b3ct         = b3*ct(i,j)
        ct_factor    = a3 + a4*salinity(i,j) + a5ct
        th0_num      = a0 + salinity(i,j)*(a1+a2*salinity(i,j)) + ct(i,j)*ct_factor
        rec_th0_den  = 1.0 / ( b0+b1*salinity(i,j) + ct(i,j)*(b2+b3ct) )
        th0(i,j)     = th0_num*rec_th0_den
        dth_dct(i,j) = ( ct_factor + a5ct - (b2+b3ct+b3ct)*th0(i,j) )*rec_th0_den        
     enddo
  enddo

  ! second estimate to get precision to roundoff
  if(pottemp_2nd_iteration) then 

      !TEOS-10 (via McDougall & Barker) uses a modified NR scheme. 
      ! The first iteration only adds half the increment 
      ! but then calculates a derivative here before performing a full update. 
      if (teos10) then 
         ct0(:,:)                   = contemp_from_pottemp(salinity,th0)   
         pt_tmp = th0 - 0.5*(ct0-ct)*dth_dct                               ! 1/2 iteration
         dth_dct                    = cp_ocean &
           /( (pt_tmp + kelvin) * dentropy_dtheta(salinity,pt_tmp) )       ! Derivative at half iteration.
         pt_tmp = th0 - (ct0-ct)*dth_dct                                   ! full iteration
         pottemp_from_contemp_level = pt_tmp - (contemp_from_pottemp(salinity,pt_tmp) - ct)*dth_dct
      else
        ct0(:,:)                   = contemp_from_pottemp(salinity,th0)
        pottemp_from_contemp_level = th0 - (ct0-ct)*dth_dct
        dth_dct                    = cp_ocean &
           /( (pottemp_from_contemp_level + kelvin) * dentropy_dtheta(salinity,pottemp_from_contemp_level) )
        pottemp_from_contemp_level = pottemp_from_contemp_level &
           -(contemp_from_pottemp(salinity,pottemp_from_contemp_level) - ct)*dth_dct
      endif

  else 

      pottemp_from_contemp_level = th0

  endif

end function pottemp_from_contemp_level
! </FUNCTION> NAME="pottemp_from_contemp_level"


!#######################################################################
! <FUNCTION NAME="pottemp_from_contemp_point">
!
! <DESCRIPTION>
!
! Compute potential temperature from conservative temperature at a point.
! Perform one extra iteration to get precision to near computer precision.
!
! Input is salinity (psu) and conservative temperature (C).
!
! </DESCRIPTION>
!
function pottemp_from_contemp_point(salinity, ct)

  real, intent(in) :: salinity
  real, intent(in) :: ct
  real             :: pottemp_from_contemp_point

  real    :: a5ct, b3ct,ct_factor
  real    :: s1,th0_num,rec_th0_den
  real    :: th0, ct0, dth_dct, pt_tmp

  if ( .not. module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_tempsalt_mod (ct_from_pt_point): module must be initialized')
  endif

  if(pottemp_equal_contemp) then
     pottemp_from_contemp_point = ct
     return 
  endif 


  s1          = salinity * sfact
  a5ct        = a5*ct
  b3ct        = b3*ct
  ct_factor   = a3 + a4*s1 + a5ct
  th0_num     = a0 + s1*(a1+a2*s1) + ct*ct_factor
  rec_th0_den = 1.0 / ( b0 + b1*s1 + ct*(b2+b3ct) )
  th0         = th0_num*rec_th0_den
  dth_dct     = ( ct_factor + a5ct - (b2+b3ct+b3ct)*th0 )*rec_th0_den

  ! second estimate to get precision to roundoff
  if(pottemp_2nd_iteration) then 

      !TEOS-10 Modified NR scheme according to McDougall & Barker.
      if (teos10) then
         ct0                        = contemp_from_pottemp(salinity,th0)
         pt_tmp                     = th0 - 0.5*(ct0-ct)*dth_dct
         dth_dct                    = cp_ocean /( (pt_tmp + kelvin) * dentropy_dtheta(salinity,pt_tmp) )
         pt_tmp                     = th0 - (ct0-ct)*dth_dct
         pottemp_from_contemp_point = pt_tmp -(contemp_from_pottemp(salinity,pt_tmp) - ct)*dth_dct
      else
         ct0                        = contemp_from_pottemp(salinity,th0)
         pottemp_from_contemp_point = th0 - (ct0-ct)*dth_dct
         dth_dct                    = cp_ocean &
           /( (pottemp_from_contemp_point + kelvin) * dentropy_dtheta(salinity,pottemp_from_contemp_point) )
         pottemp_from_contemp_point = pottemp_from_contemp_point &
           -(contemp_from_pottemp(salinity,pottemp_from_contemp_point) - ct)*dth_dct
      endif

  else 

      pottemp_from_contemp_point = th0

  endif


end function pottemp_from_contemp_point
! </FUNCTION> NAME="pottemp_from_contemp_point"



!#######################################################################
! <FUNCTION NAME="dentropy_dtheta_field">
!
! <DESCRIPTION>
!
!   d(entropy)/d(pottemp) at each grid point from twice differentiating 
!   the Gibbs potential in Feistel (2003), Prog. Ocean. 58, 43-114.
!   (pressure=0 since use potential temperature) 
!
!   salinity        : salinity (psu)
!   theta           : potential temperature (deg C, ITS-90)
!   dentropy_dtheta : d(entropy)/d(pottemp) J/(kg degC^2)
!
!   check value: dentropy_dtheta(35,20) = 13.63256369213874
!
! </DESCRIPTION>
!
function dentropy_dtheta_field(salinity,theta)

  real, dimension(isd:,jsd:,:), intent(in) :: salinity
  real, dimension(isd:,jsd:,:), intent(in) :: theta 
  real, dimension(isd:ied,jsd:jed,nk)      :: dentropy_dtheta_field

  integer :: i,j,k 
  real    :: x2, x, y, de_dt 

  if (teos10) then

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied

              x2 = one_over_40*salinity(i,j,k)
              x  = sqrt(x2)
              y  = one_over_40*theta(i,j,k)
              de_dt= g1 + y*(g2 + y*(g3 + y*(g4+y*(g5+y*g6))))          + &
                   x2*(h1+x*(h2+x*(h3+y*(h4+y*(h5+y*h6))) + y*(h7+y*h8))+ &
                       y*(h9+y*(h10+y*(h11 +y*h12))))
              dentropy_dtheta_field(i,j,k) = (-one_over_1600)*de_dt
           enddo
        enddo
     enddo

  else

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied

              x2 = one_over_40*salinity(i,j,k)
              x  = sqrt(x2)
              y  = one_over_40*theta(i,j,k)

              de_dt = f0 + x2*(f1 + x*(f2 + y*(f3  + f4*y)) + y*(f5 + y*(f6 + f7*y))) &
                         +  y*(f8 + y*(f9 + y*(f10 + (f11 - f12*y)*y)))

              dentropy_dtheta_field(i,j,k) = f13*de_dt

           enddo
        enddo
     enddo

  endif

end function dentropy_dtheta_field
! </FUNCTION> NAME="dentropy_dtheta_field"


!#######################################################################
! <FUNCTION NAME="dentropy_dtheta_level">
!
! <DESCRIPTION>
!
!   d(entropy)/d(pottemp) at k-level from twice differentiating 
!   the Gibbs potential in Feistel (2003), Prog. Ocean. 58, 43-114.
!   (pressure=0 since use potential temperature) 
!
!   salinity        : salinity (psu)
!   theta           : potential temperature (deg C, ITS-90)
!   dentropy_dtheta : d(entropy)/d(pottemp) J/(kg degC^2)
!
!   check value: dentropy_dtheta(35,20) = 13.63256369213874
!
! </DESCRIPTION>
!
function dentropy_dtheta_level(salinity,theta)

  real, dimension(isd:,jsd:), intent(in) :: salinity
  real, dimension(isd:,jsd:), intent(in) :: theta 
  real, dimension(isd:ied,jsd:jed)       :: dentropy_dtheta_level

  integer :: i,j 
  real    :: x2, x, y, de_dt 

  if (teos10) then

     do j=jsd,jed
        do i=isd,ied

           x2 = one_over_40*salinity(i,j)
           x  = sqrt(x2)
           y  = one_over_40*theta(i,j)
           de_dt= g1 + y*(g2 + y*(g3 + y*(g4+y*(g5+y*g6))))          + &
                x2*(h1+x*(h2+x*(h3+y*(h4+y*(h5+y*h6))) + y*(h7+y*h8))+ &
                    y*(h9+y*(h10+y*(h11 +y*h12))))
           dentropy_dtheta_level = (-one_over_1600)*de_dt
        enddo
     enddo


  else


     do j=jsd,jed
        do i=isd,ied
   
           x2 = one_over_40*salinity(i,j)
           x  = sqrt(x2)
           y  = one_over_40*theta(i,j)

           de_dt = f0 + x2*(f1 + x*(f2 + y*(f3  + f4*y)) + y*(f5 + y*(f6 + f7*y))) &
                      +  y*(f8 + y*(f9 + y*(f10 + (f11 - f12*y)*y)))
   
           dentropy_dtheta_level(i,j) = f13*de_dt

        enddo
     enddo

   endif

end function dentropy_dtheta_level
! </FUNCTION> NAME="dentropy_dtheta_level"


!#######################################################################
! <FUNCTION NAME="dentropy_dtheta_point">
!
! <DESCRIPTION>
!
!   d(entropy)/d(pottemp) at an (i,j,k) point from twice differentiating 
!   the Gibbs potential in Feistel (2003), Prog. Ocean. 58, 43-114.
!   (pressure=0 since use potential temperature) 
!
!   salinity        : salinity (psu)
!   theta           : potential temperature (deg C, ITS-90)
!   dentropy_dtheta : d(entropy)/d(pottemp) J/(kg degC^2)
!
!   check value: dentropy_dtheta(35,20) = 13.63256369213874
!
! </DESCRIPTION>
!
function dentropy_dtheta_point(salinity,theta)

  real, intent(in) :: salinity
  real, intent(in) :: theta 

  real :: dentropy_dtheta_point
  real :: x2, x, y, de_dt 

  x2 = one_over_40*sfac*salinity
  x  = sqrt(x2)
  y  = one_over_40*theta

  if (teos10) then
    de_dt= g1 + y*(g2 + y*(g3 + y*(g4+y*(g5+y*g6))))              + &
            x2*(h1+x*(h2+x*(h3+y*(h4+y*(h5+y*h6))) + y*(h7+y*h8)) + &
                 y*(h9+y*(h10+y*(h11 +y*h12))))
    dentropy_dtheta_point = (-one_over_1600)*de_dt
  else
     de_dt = f0 + x2*(f1 + x*(f2 + y*(f3  + f4*y)) + y*(f5 + y*(f6 + f7*y))) &
             +  y*(f8 + y*(f9 + y*(f10 + (f11 - f12*y)*y)))
     dentropy_dtheta_point = f13*de_dt
  endif 

end function dentropy_dtheta_point
! </FUNCTION> NAME="dentropy_dtheta_point"


!#######################################################################
! <SUBROUTINE NAME="tempsalt_check_range">
!
! <DESCRIPTION>
!
!   Check to see that temperature and salinity are within preset
!   range. If outside of the range, then bring the model down.
!
!   This check is particularly useful to ensure that the equation
!   of state is evaluated with physically sensible values of temp
!   and salinity.  
!
! </DESCRIPTION>
!
subroutine tempsalt_check_range(salinity,theta,depth)

  real, dimension(isd:,jsd:,:), intent(in) :: salinity
  real, dimension(isd:,jsd:,:), intent(in) :: theta 
  real, dimension(isd:,jsd:,:), intent(in) :: depth
  integer :: i,j,k
  
  character(len=512) :: err_msg

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(salinity(i,j,k) < s_min .or. salinity(i,j,k) > s_max) then
               write(err_msg,'(a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)')      &
               ' Error: salinity out of range with value ', salinity(i,j,k),                &
               ' at (i,j,k) = (',i+Dom%ioff,',',j+Dom%joff,',',k,                           &
               '),  (lon,lat,dpt) = (',Grd%xt(i,j),',',                                     &
               Grd%yt(i,j),',', depth(i,j,k),' m)'
               call mpp_error(FATAL, trim(err_msg))
           endif 
           if(theta(i,j,k) < t_min .or. theta(i,j,k) > t_max) then
               write(err_msg,'(a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)')      &
               ' Error: temperature out of range with value ', theta(i,j,k),                &
               ' at (i,j,k) = (',i+Dom%ioff,',',j+Dom%joff,',',k,                           &
               '),  (lon,lat,dpt) = (',Grd%xt(i,j),',',                                     &
               Grd%yt(i,j),',', depth(i,j,k),' m)'
               call mpp_error(FATAL, trim(err_msg))
           endif
        enddo
     enddo
  enddo


end subroutine tempsalt_check_range
! </SUBROUTINE> NAME="tempsalt_check_range"



end module ocean_tempsalt_mod
