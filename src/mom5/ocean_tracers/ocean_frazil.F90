module ocean_frazil_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> R.A.S. Fiedler
!</CONTACT>
!
!<CONTACT EMAIL="David.Jackett@csiro.au"> David Jackett
!</CONTACT>
!
!<OVERVIEW>
! This module computes the heating of seawater due to 
! frazil ice formation.
!</OVERVIEW>
!
!<DESCRIPTION>
! Frazil can generally form at any vertical level, although 
! it is common for climate models to assume it is formed 
! only at k=1.  In this case, pressure is assumed to be 
! atmospheric for purposes of computing the freezing 
! temperature of seawater.  
!
! The freezing temperature of seawater is computed one of two 
! possible ways: 
!
! (1) simple way uses a linear function of salinity
! and assumes zero (i.e. atmospheric) pressure,
!
! tfreeze(deg C) = a1*salinity(psu)
! with a1 = -0.054 
!
! (2) accurate way uses a nonlinear function of 
! salinity(psu) and gauge pressure(dbar), where 
! gauge pressure=absolute pressure - 10.1325 dbar.
!
! tfreeze (deg C) = tf_num/tf_den
! tf_num  =  a0 + s*(a1 + sqrt(s)*(a2 + sqrt(s)*a3)) + p*(a4 + p*(a5 + s*a6)) 
! tf_dem  =  b0 + p*(b1 + p*b2) + s*s*sqrt(s)*b3
!
! check value : fp_theta(35,200,'air-sat')  = -2.076426227617581 deg C
!               fp_theta(35,200,'air-free') = -2.074408175943127 deg C
!
! This method results in a more accurate freezing 
! temperature than the simpler approach.  It is also 
! important for ice-shelf modelling to include a 
! pressure dependence when the shelf penetrates
! into the water column. 
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
! seawater - 2010: Calculation and use of thermodynamic properties.  
! Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
! UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
! See sections 3.33 and 3.34 of this TEOS-10 Manual.
! </REFERENCE>
!
! <REFERENCE>
! "Updated algorithms for density, potential temperature,
! conservative temperature and freezing temperature of 
! seawater", Jackett, McDougall, Feistel, Wright, and Griffies
! Journal of Atmospheric and Oceanic Technology, in press 2005. 
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_frazil_nml">
!
!<DATA NAME="use_this_module" TYPE="logical">
! If true, then compute frazil heating. 
!</DATA> 
!<DATA NAME="debug_this_module" TYPE="logical">
! For debugging this module
!</DATA> 
!
!<DATA NAME="frazil_factor" UNITS="dimensionless" TYPE="real">
! This factor accounts for possibly different time stepping used 
! in the sea ice model relative to the ocean model.  If sea-ice 
! and ocean use same time stepping schemes, then frazil_factor=1.0.
! If sea-ice uses a twolevel scheme and ocean a threelevel leap-frog,
! then frazil_factor=0.5. Default is 1.0 since the GFDL sea ice model 
! SIS uses  a two-level time stepping scheme and MOM defaults to 
! a staggered two-level scheme. 
!</DATA> 
!
!<DATA NAME="freezing_temp_simple" TYPE="logical">
! To use the simplified freezing point temperature of seawater,
! as used in MOM4.0. 
! </DATA> 
!<DATA NAME="freezing_temp_accurate_preteos10" TYPE="logical">
! To use the accurate freezing point temperature of seawater,
! which is a nonlinear function of salinity and pressure. 
!</DATA> 
!<DATA NAME="freezing_temp_accurate_teos10" TYPE="logical">
! To use the accurate freezing point temperature of seawater 
! according to TEOS-10, which is a nonlinear function of  
! salinity and pressure. 
! This equation is recommended for use with ice-shelf modelling.  
!</DATA> 
!
!<DATA NAME="frazil_only_in_surface" TYPE="logical">
! For typical case where compute frazil heating only in 
! the surface grid cell.  Will assume the gauge 
! pressure is zero in this case when computing freezing
! temperature.  
!</DATA> 
!
!</NAMELIST>

use constants_mod,     only: epsln 
use diag_manager_mod,  only: register_diag_field
use fms_mod,           only: open_namelist_file, check_nml_error, close_file
use fms_mod,           only: FATAL, NOTE, stdout, stdlog
use mpp_mod,           only: input_nml_file, mpp_error 
use mpp_domains_mod,   only: mpp_global_sum, NON_BITWISE_EXACT_SUM

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value, cp_ocean
use ocean_parameters_mod, only: TWO_LEVEL, THREE_LEVEL
use ocean_tpm_util_mod,   only: otpm_set_diag_tracer
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_diag_tracer_type, ocean_density_type
use ocean_workspace_mod,  only: wrk1_2d, wrk1
use ocean_util_mod,       only: diagnose_2d, diagnose_3d

#if defined(ACCESS)
use auscom_ice_parameters_mod, only: pop_icediag, do_ice
use auscom_ice_mod,            only: AQICE
use auscom_ice_mod,            only: auscom_ice_formation_new
use diag_manager_mod,          only: send_data
#endif

implicit none

private 

#include <ocean_memory.h>

logical :: module_initialized = .false.

character(len=256) :: version='CVS $$'
character(len=256) :: tagname='Tag $Name: tikal $'

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

integer :: index_frazil=-1
integer :: index_temp=-1
integer :: index_salt=-1

! for FAFMIP heat tracers 
integer :: index_frazil_redist=-1
integer :: index_redist_heat  =-1

integer :: tendency
real    :: dtimer
real    :: cp_ocean_r

! for diagnostics 
logical :: used
integer :: id_frazil_2d=-1
integer :: id_frazil_3d=-1
integer :: id_temp_freeze=-1
! for FAFMIP heat tracers 
integer :: id_frazil_redist_2d=-1
integer :: id_frazil_redist_3d=-1

! for ascii output
integer :: unit=6

! Fraction of air saturation 
! Now either 0 or 1. do we want to change this?
real :: saturation_fraction

! coefficients in freezing temperature with preTEOS10 
real :: a0, a1, a2, a3, a4, a5, a6
real :: b0, b1, b2, b3
real :: c1, c2

! for TEOS10 freezing temperature 
real,parameter  ::  v0 =  0.017947064327968736d0
real,parameter  ::  v1 =  -6.076099099929818d0
real,parameter  ::  v2 =   4.883198653547851d0
real,parameter  ::  v3 =  -11.88081601230542d0
real,parameter  ::  v4 =   13.34658511480257d0
real,parameter  ::  v5 =  -8.722761043208607d0
real,parameter  ::  v6 =   2.082038908808201d0
real,parameter  ::  v7 =  -7.389420998107497d0
real,parameter  ::  v8 =  -2.110913185058476d0
real,parameter  ::  v9 =   0.2295491578006229d0
real,parameter  ::  v10 = -0.9891538123307282d0
real,parameter  ::  v11 = -0.08987150128406496d0
real,parameter  ::  v12 =  0.3831132432071728d0
real,parameter  ::  v13 =  1.054318231187074d0
real,parameter  ::  v14 =  1.065556599652796d0
real,parameter  ::  v15 = -0.7997496801694032d0
real,parameter  ::  v16 =  0.3850133554097069d0
real,parameter  ::  v17 = -2.078616693017569d0
real,parameter  ::  v18 =  0.8756340772729538d0
real,parameter  ::  v19 = -2.079022768390933d0
real,parameter  ::  v20 =  1.596435439942262d0
real,parameter  ::  v21 =  0.1338002171109174d0
real,parameter  ::  v22 =  1.242891021876471d0


public compute_frazil_heating
public ocean_frazil_init
private compute_frazil_redist_heating 

! nml defaults 
logical :: use_this_module         = .false.
logical :: debug_this_module       = .false.
logical :: freezing_temp_simple    = .false.
logical :: freezing_temp_preteos10 = .false.
logical :: freezing_temp_teos10    = .false.
logical :: air_saturated_water     = .true. 
logical :: frazil_only_in_surface  = .true. 
real    :: frazil_factor           = 1.0

namelist /ocean_frazil_nml/ use_this_module, debug_this_module,                                  &
                            freezing_temp_simple, freezing_temp_preteos10, freezing_temp_teos10, &
                            frazil_factor, air_saturated_water, frazil_only_in_surface

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_frazil_init">
!
! <DESCRIPTION>
! Initialization code for the frazil diagnostic tracer.
! </DESCRIPTION>
!
subroutine ocean_frazil_init (Domain, Grid, Time, Time_steps, Ocean_options, &
                              itemp, isalt, iredist, debug)  
   
  type(ocean_domain_type),     intent(in), target   :: Domain
  type(ocean_grid_type),       intent(in), target   :: Grid
  type(ocean_time_type),       intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)           :: Time_steps
  type(ocean_options_type),    intent(inout)        :: Ocean_options
  integer,                     intent(in)           :: itemp
  integer,                     intent(in)           :: isalt
  integer,                     intent(in)           :: iredist 
  logical,                     intent(in), optional :: debug

  integer :: ierr, ioun, io_status
  integer :: num_freeze
  real    :: s, sqrts, press
  real    :: tf_num, tf_den
  real    :: tfreeze, tfreeze_check 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if (module_initialized) then
    call mpp_error(FATAL, '==>Error: ocean_frazil_mod: module already initialized')
  endif

  write( stdlogunit,'(/a/)') trim(version)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_frazil_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_frazil_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_frazil_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_frazil_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_frazil_nml)  
  write (stdlogunit, ocean_frazil_nml)

  if(.not. use_this_module) then 
      write(stdoutunit,'(a)') &
      '==>Note: NOT running with frazil. Be sure there is no frazil entry in field_table.'  
      call mpp_error(NOTE, '==>Note: ocean_frazil_mod: NOT using frazil heating.')
      Ocean_options%frazil_ice = 'Did NOT use frazil sea ice formation.'
      return 
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 

  call mpp_error(NOTE, '==>Note from ocean_frazil_mod: USING frazil heating.')
  if(debug_this_module) then 
      write(stdoutunit,'(a)') '==>Note: running with debug_this_module=.true.'  
  endif

  num_freeze = 0
  if(freezing_temp_simple) then 
      write(stdoutunit,'(a)') &
      '==>Note: Using simple equation for seawater freezing temperature.'
      num_freeze = num_freeze + 1 
  endif
  if(freezing_temp_preteos10) then 
      write(stdoutunit,'(a)') &
      '==>Note: Using preTEOS10 equation for seawater freezing temperature.'
      if(air_saturated_water) then 
         write(stdoutunit,'(a)') &
         '==>Note: Assuming seawater is saturated w/ air for freezing temperature calculation.'
      endif 
      num_freeze = num_freeze + 1 
  endif
  if(freezing_temp_teos10) then 
      write(stdoutunit,'(a)') &
      '==>Note: Using TEOS10 equation for seawater freezing temperature.'
      if(air_saturated_water) then 
         write(stdoutunit,'(a)') &
         '==>Note: Assuming seawater is saturated w/ air for freezing temperature calculation.'
      endif 
      num_freeze = num_freeze + 1 
  endif
  if(num_freeze > 1) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_frazil_mod: Choose only one freezing equation for frazil heating.')
  endif 
  if(num_freeze == 0) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_frazil_mod: Choose one freezing equation for frazil heating.')
  endif 

  if(frazil_only_in_surface) then
      write(stdoutunit,'(a)') &
      'Assuming that frazil forms only in the surface(k=1) ocean grid cell.'
      write(stdoutunit,'(a)') &
      'Setting gauge pressure to zero when computing seawater freezing temperature.'
      Ocean_options%frazil_ice = 'Frazil sea ice formed through heat exchange only in ocean surface.'
  endif
  if(.not. frazil_only_in_surface) then
      if(freezing_temp_simple) then 
          call mpp_error(FATAL, &
               '==>Error from ocean_frazil_mod: Cannot choose simple freezing temp when have frazil forming at depths k > 1.')
      else 
          write(stdoutunit,'(a)') &
          'Assuming that frazil can form in any ocean grid cell, including cells beneath surface.'
          Ocean_options%frazil_ice = 'Frazil sea ice anywhere in ocean column with an accurate freezing temp equation.'
      endif
  endif

  index_frazil = otpm_set_diag_tracer('frazil',                                       &
     caller='ocean_frazil_mod/ocean_frazil_init',                                     &
     longname='frazil heating', units='J/m^2',                                        &
     conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20,                    &
     min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0, &
     restart_file='ocean_frazil.res.nc' )

  ! The diagnostic tracer "frazil_redist" is only used for FAFMIP simulations.  
  if(iredist > 0) then
    write(stdoutunit,'(a)') &
    'Initializing frazil_redist heat diagnostic tracer for FAFMIP purposes.'
    index_frazil_redist = otpm_set_diag_tracer('frazil_redist',                         &
       caller='ocean_frazil_mod/ocean_frazil_init',                                     &
       longname='frazil redistributed heat field heating', units='J/m^2',               &
       conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20,                    &
       min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0, &
       restart_file='ocean_frazil_redist.res.nc' )
  else
    write(stdoutunit,'(a)') &
          'NOT Initializing frazil_redist heat diagnostic tracer for FAFMIP purposes.'
  endif

  tendency          = Time_steps%tendency
  dtimer            = 1.0/(Time_steps%dtime_t + epsln)
  cp_ocean_r        = 1.0/cp_ocean
  index_salt        = isalt
  index_temp        = itemp
  index_redist_heat = iredist 

  Grd => Grid
  Dom => Domain

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif 

  write(stdoutunit,'(/a,f5.2)') '==>Note from ocean_frazil_mod: using frazil_factor= ',frazil_factor
  if(tendency ==TWO_LEVEL .and. frazil_factor==0.5) then 
      call mpp_error(FATAL,&
      '==>ocean_sbc_mod: with 2-level tendencies for MOM and SIS, set frazil_factor==1.0 ')
      write(stdoutunit,'(/a)') &
      '==>Error in ocean_sbc_mod: SIS is two-time ice level model.  You are also running'
      write(stdoutunit,'(a) ') &
      '  MOM with two-time level scheme. In this case, frazil_factor==1.0 must be set.'
      write(stdoutunit,'(a/) ')&
       '  With another ice model w/ different time stepping, set frazil_factor accordingly.'
  endif
  if(tendency ==THREE_LEVEL .and. frazil_factor==1.0) then 
      call mpp_error(FATAL,&
      '==>ocean_frazil_mod: with 3-level tendency for MOM and 2-level for SIS, set frazil_factor==.50 ')
      write(stdoutunit,'(/a)') &
      '==>Error in ocean_sbc_mod: SIS is two-time level ice model.  You are running MOM'
      write(stdoutunit,'(a) ') &
      '   with three-time level leap-frog. In this case, frazil_factor==.50 should be set.'
      write(stdoutunit,'(a/) ') &
      '   With another ice model w/ different time stepping, set frazil_factor accordingly.'
  endif


  ! coefficients in the preTEOS10 freezing point of seawater
  ! assume prognostic temperature variable is potential temp
  ! note that TEOS10 coefficients set as parameters above. 
  a0 =  2.5180516744541290e-03 
  a1 = -5.8545863698926184e-02
  a2 =  2.2979985780124325e-03
  a3 = -3.0086338218235500e-04
  a4 = -7.0023530029351803e-04
  a5 =  8.4149607219833806e-09
  a6 =  1.1845857563107403e-11

  b0 =  1.0000000000000000e+00
  b1 = -3.8493266309172074e-05
  b2 =  9.1686537446749641e-10
  b3 =  1.3632481944285909e-06

  if(air_saturated_water) then 
      saturation_fraction=1.0
      if(freezing_temp_preteos10) then 
          c1 = -2.5180516744541290e-03 
          c2 =  1.428571428571429e-05
          tfreeze_check = -2.076426227617581
      elseif(freezing_temp_teos10) then 
          c1 = 0.014289763856964d0
          c2 = 0.05700064989972d0
          tfreeze_check = -2.062635500704721
      endif
  else 
      saturation_fraction = 0.0
      c1 = 0.0
      c2 = 0.0
      if (freezing_temp_preteos10) then
          tfreeze_check = -2.074408175943127
      elseif(freezing_temp_teos10) then 
          tfreeze_check = -2.060735134187700
      endif
  endif


  if( freezing_temp_simple) then 
     a1 = -0.054
  elseif(freezing_temp_preteos10) then
     s       = 35.0
     press   = 200.0
     sqrts   = sqrt(s)
     tf_num  = a0 + s*(a1 + sqrts*(a2 + sqrts*a3)) + press*(a4 + press*(a5 + s*a6)) 
     tf_den  = b0 + press*(b1 + press*b2) + s*s*sqrts*b3 
     tfreeze = tf_num/tf_den + (c1 + s*c2) 
     write(stdoutunit,'(a,e24.16)')'Check value for freezing temperature(C) at (35psu,200dbar) = ',tfreeze 
     write(stdoutunit,'(a,e24.16)')'This value differs from published check value by ', tfreeze-tfreeze_check    
  elseif(freezing_temp_teos10 ) then
     s       = 1.e-2*35.0
     press   = 1.e-4*200.0
     sqrts   = sqrt(s)
     tfreeze =  v0 &
                + s*(v1 + sqrts*(v2 + sqrts*(v3 + sqrts*(v4 + sqrts*(v5 + v6*sqrts)))))         &
                + press*(v7 + press*(v8 + v9*press))                                            &
                + s*press*(v10 + press*(v12 + press*(v15 + v21*s)) + s*(v13 + v17*press+ v19*s) &
                + sqrts*(v11 + press*(v14 + v18*press)  + s*(v16 + v20*press+ v22*s)))         &
                - saturation_fraction*(1e-3)*(2.4 - c1*s)*(1 + c2*(1d0 - s/35.16504d0))
     write(stdoutunit,'(a,e24.16)')'Check value for freezing temperature(C) at (35psu,200dbar) = ',tfreeze
     write(stdoutunit,'(a,e24.16)')'This value differs from published check value by ', tfreeze-tfreeze_check
  endif 


  ! when ice forms only in k=1 cell, only require frazil saved as 2d field 
  id_frazil_2d = register_diag_field ('ocean_model', 'frazil_2d', Grd%tracer_axes(1:2), &
         Time%model_time, 'ocn frazil heat flux over time step', 'W/m^2',               &
         missing_value=missing_value, range=(/-1.e10,1.e10/))  

  ! when ice forms at any depth, require frazil to be saved as 3d field
  id_frazil_3d = register_diag_field ('ocean_model', 'frazil_3d', Grd%tracer_axes(1:3), &
         Time%model_time, 'ocn frazil heat flux over time step', 'W/m^2',               &
         missing_value=missing_value, range=(/-1.e10,1.e10/))  

  ! for FAFMIP frazil diagnostic tracer 
  id_frazil_redist_2d = register_diag_field ('ocean_model','frazil_redist_2d',Grd%tracer_axes(1:2), &
        Time%model_time, 'ocn frazil redist heat flux over time step (for FAFMIP)', 'W/m^2',        &
        missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_frazil_redist_3d = register_diag_field ('ocean_model', 'frazil_redist_3d', Grd%tracer_axes(1:3), &
         Time%model_time, 'ocn frazil heat flux over time step (for FAFMIP)', 'W/m^2',                &
         missing_value=missing_value, range=(/-1.e10,1.e10/))  
   
  id_temp_freeze = register_diag_field ('ocean_model', 'temp_freeze', Grd%tracer_axes(1:3), &
         Time%model_time, 'freezing temperature of seawater', 'deg C',                      &
         missing_value=missing_value, range=(/-1.e10,1.e10/))  

end subroutine ocean_frazil_init  
! </SUBROUTINE> NAME="ocean_frazil_init">


!#######################################################################
! <SUBROUTINE NAME="compute_frazil_heating">
!
! <DESCRIPTION>
! Compute ocean heating due to formation of frazil-ice (Joules/m^2)
!
! Note that "frazil_factor" accounts for possibly different time 
! stepping used in ocean model and the sea ice model.  With MOM 
! using a leap-frog, and the GFDL ocean model SIS using forward,
! then frazil_factor=0.5. If use recommended tendency=twolevel 
! in MOM, then frazil_factor=1.0
!
! </DESCRIPTION>
!
subroutine compute_frazil_heating (Time, Thickness, Dens, T_prog, T_diag)

  type(ocean_time_type),        intent(in)    :: Time 
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type), intent(inout) :: T_diag(:)

  integer  :: i,j,k
  integer  :: taup1, tau
  real     :: tf_num, tf_den, tfreeze
  real     :: s, sqrts
  real     :: press 

  if(.not. use_this_module) return 

#if defined(ACCESS)
  if (pop_icediag) then
    if ( do_ice ) then
      T_diag(index_frazil)%field = 0.0
      call auscom_ice_formation_new(Time,T_prog,Thickness, T_diag(index_frazil) )
    endif
    T_diag(index_frazil)%field=T_diag(index_frazil)%field * Grd%tmask
    if (id_frazil_2d > 0) then 
      used = send_data(id_frazil_2d, T_diag(index_frazil)%field(:,:,1)*dtimer, &
           Time%model_time, rmask=Grd%tmask(:,:,1), &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif
    if (id_frazil_3d > 0) then 
      used = send_data(id_frazil_3d, T_diag(index_frazil)%field(:,:,:)*dtimer, &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
  
    return
  endif
#endif

  taup1 = Time%taup1
  tau   = Time%tau

  if(freezing_temp_simple) then 

     k=1
      do j=jsc,jec
         do i=isc,iec 
            T_diag(index_frazil)%field(i,j,k) = 0.0
            if(Grd%tmask(i,j,k) > 0.0) then 
                tfreeze = a1*Dens%rho_salinity(i,j,k,taup1)
                if(T_prog(index_temp)%field(i,j,k,taup1) < tfreeze) then  
                    T_diag(index_frazil)%field(i,j,k) = &
                         (tfreeze-T_prog(index_temp)%field(i,j,k,taup1)) &
                         *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor 
                    T_prog(index_temp)%field(i,j,k,taup1) = tfreeze
                endif
            endif
#if defined(ACCESS)
            AQICE(i,j)=AQICE(i,j) - T_diag(index_frazil)%field(i,j,k)
            !'-' is required in routine ice_heatflux for merge operation.
#endif
         enddo
      enddo
      
  elseif(freezing_temp_preteos10) then  

      if(frazil_only_in_surface) then 
          k=1     
          do j=jsc,jec
             do i=isc,iec 
                T_diag(index_frazil)%field(i,j,k) = 0.0
                if(Grd%tmask(i,j,k) > 0.0) then 
                    s       = T_prog(index_salt)%field(i,j,k,taup1)              
                    sqrts   = sqrt(s)
                    tf_num  = a0 + s*(a1 + sqrts*(a2 + sqrts*a3))
                    tf_den  = b0 + s*s*sqrts*b3 
                    tfreeze = tf_num/tf_den  + (c1 + s*c2) 
                    if(T_prog(index_temp)%field(i,j,k,taup1) < tfreeze) then 
                        T_diag(index_frazil)%field(i,j,k) = &
                             (tfreeze-T_prog(index_temp)%field(i,j,k,taup1)) &
                             *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor 
                        T_prog(index_temp)%field(i,j,k,taup1) = tfreeze
                    endif
                endif
#if defined(ACCESS)
                AQICE(i,j)=AQICE(i,j) - T_diag(index_frazil)%field(i,j,1)
                !'-' is required in routine ice_heatflux for merge operation.
#endif
             enddo
          enddo

      else 
#if defined(ACCESS)
       do j=jsc,jec
             do i=isc,iec
                do k=1,nk
#else
          do k=1,nk 
             do j=jsc,jec
                do i=isc,iec 
#endif
                   T_diag(index_frazil)%field(i,j,k) = 0.0
                   if(Grd%tmask(i,j,k) > 0.0) then 
                       s       = Dens%rho_salinity(i,j,k,taup1)
                       press   = Dens%pressure_at_depth(i,j,k)              
                       sqrts   = sqrt(s)
                       tf_num  = a0 + s*(a1 + sqrts*(a2 + sqrts*a3)) + press*(a4 + press*(a5 + s*a6)) 
                       tf_den  = b0 + press*(b1 + press*b2) + s*s*sqrts*b3 
                       tfreeze = tf_num/tf_den + (c1 + s*c2) 
                       if(T_prog(index_temp)%field(i,j,k,taup1) < tfreeze) then 
                           T_diag(index_frazil)%field(i,j,k) = &
                                (tfreeze-T_prog(index_temp)%field(i,j,k,taup1)) &
                                *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor 
                           T_prog(index_temp)%field(i,j,k,taup1) = tfreeze
                       endif
                   endif
#if defined(ACCESS)
                   AQICE(i,j)=AQICE(i,j) - T_diag(index_frazil)%field(i,j,k)
                   !'-' is required in routine ice_heatflux for merge operation.
#endif
                enddo
             enddo
          enddo

      endif

  else  ! teos10 

      if(frazil_only_in_surface) then 
          k=1
          do j=jsc,jec
             do i=isc,iec
                T_diag(index_frazil)%field(i,j,k) = 0.0
                if(Grd%tmask(i,j,k) > 0.0) then
                    s       = 1.d-2*Dens%rho_salinity(i,j,k,taup1)
                    sqrts   = sqrt(s)
                    tfreeze =  v0 &
                               + s*(v1 + sqrts*(v2 + sqrts*(v3 + sqrts*(v4 + sqrts*(v5 + v6*sqrts)))))  &
                               - saturation_fraction*(1e-3)*(2.4 - c1*s)*(1 + c2*(1.d0 - s/35.16504d0))
                    if(T_prog(index_temp)%field(i,j,k,taup1) < tfreeze) then
                        T_diag(index_frazil)%field(i,j,k) = &
                             (tfreeze-T_prog(index_temp)%field(i,j,k,taup1)) &
                             *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor
                        T_prog(index_temp)%field(i,j,k,taup1) = tfreeze
                    endif
                endif
             enddo
          enddo

      else

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  T_diag(index_frazil)%field(i,j,k) = 0.0
                  if(Grd%tmask(i,j,k) > 0.0) then
                     s       = 1.d-2*Dens%rho_salinity(i,j,k,taup1)
                     sqrts   = sqrt(s)
                     press   = Dens%pressure_at_depth(i,j,k)              
                     tfreeze =  v0 &
                                + s*(v1 + sqrts*(v2 + sqrts*(v3 + sqrts*(v4 + sqrts*(v5 + v6*sqrts)))))         &
                                + press*(v7 + press*(v8 + v9*press))                                            &
                                + s*press*(v10 + press*(v12 + press*(v15 + v21*s)) + s*(v13 + v17*press+ v19*s) &
                                + sqrts*(v11 + press*(v14 + v18*press)  + s*(v16 + v20*press+ v22*s)))          &
                                - saturation_fraction*(1e-3)*(2.4 - c1*s)*(1 + c2*(1d0 - s/35.16504d0))
                     if(T_prog(index_temp)%field(i,j,k,taup1) < tfreeze) then
                         T_diag(index_frazil)%field(i,j,k) = &
                              (tfreeze-T_prog(index_temp)%field(i,j,k,taup1)) &
                              *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor
                         T_prog(index_temp)%field(i,j,k,taup1) = tfreeze
                     endif
                  endif
               enddo
            enddo
         enddo
     endif

  endif  ! endif for freezing_temp_simple


  if (id_frazil_2d > 0) then
     call diagnose_2d(Time, Grd, id_frazil_2d, T_diag(index_frazil)%field(:,:,1)*dtimer)
  endif
  if (id_frazil_3d > 0) then
     call diagnose_3d(Time, Grd, id_frazil_3d, T_diag(index_frazil)%field(:,:,:)*dtimer)
  endif

  if(id_temp_freeze > 0) then 
     wrk1(:,:,:) =0.0
     if(freezing_temp_simple) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec 
                 wrk1(i,j,k) = a1*T_prog(index_salt)%field(i,j,k,taup1)*Grd%tmask(i,j,k)
              enddo
           enddo
        enddo
     else 
        do k=1,nk 
           do j=jsc,jec
              do i=isc,iec 
                 if(Grd%tmask(i,j,k) > 0.0) then 
                    s       = T_prog(index_salt)%field(i,j,k,taup1)
                    press   = Dens%pressure_at_depth(i,j,k)              
                    sqrts   = sqrt(s)
                    tf_num  = a0 + s*(a1 + sqrts*(a2 + sqrts*a3)) + press*(a4 + press*(a5 + s*a6)) 
                    tf_den  = b0 + press*(b1 + press*b2) + s*s*sqrts*b3 
                    wrk1(i,j,k) = tf_num/tf_den + (c1 + s*c2) 
                 endif
              enddo
           enddo
        enddo
     endif
     call diagnose_3d(Time, Grd, id_temp_freeze, wrk1(:,:,:))
  endif

  if(index_frazil_redist > 0) then
     call compute_frazil_redist_heating(Time, Thickness, Dens, T_prog, T_diag)
  endif 

end subroutine compute_frazil_heating
! </SUBROUTINE> NAME="compute_frazil_heating">


!#######################################################################
! <SUBROUTINE NAME="compute_frazil_redist_heating">
!
! <DESCRIPTION>
! Compute ocean heating due to formation of frazil-ice (Joules/m^2)
! for the FAFMIP redistributed heat tracer.  
!
! Note that "frazil_factor" accounts for possibly different time 
! stepping used in ocean model and the sea ice model.  With MOM 
! using a leap-frog, and the GFDL ocean model SIS using forward,
! then frazil_factor=0.5. If use recommended tendency=twolevel 
! in MOM, then frazil_factor=1.0
!
! Note: unsure if the ACCESS ifedf options have been properly
! handled. Need ACCESS developers to check.  
! </DESCRIPTION>
!
subroutine compute_frazil_redist_heating (Time, Thickness, Dens, T_prog, T_diag)

  type(ocean_time_type),        intent(in)    :: Time 
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type), intent(inout) :: T_diag(:)

  integer  :: i,j,k
  integer  :: taup1, tau
  real     :: tf_num, tf_den, tfreeze
  real     :: s, sqrts
  real     :: press 

  if(.not. use_this_module) return 

  taup1 = Time%taup1
  tau   = Time%tau

  if(freezing_temp_simple) then 

     k=1
      do j=jsc,jec
         do i=isc,iec 
            T_diag(index_frazil_redist)%field(i,j,k) = 0.0
            if(Grd%tmask(i,j,k) > 0.0) then 
                tfreeze = a1*Dens%rho_salinity(i,j,k,taup1)
                if(T_prog(index_redist_heat)%field(i,j,k,taup1) < tfreeze) then  
                    T_diag(index_frazil_redist)%field(i,j,k) = &
                         (tfreeze-T_prog(index_redist_heat)%field(i,j,k,taup1)) &
                         *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor 
                    T_prog(index_redist_heat)%field(i,j,k,taup1) = tfreeze
                endif
            endif
         enddo
      enddo

  elseif(freezing_temp_preteos10) then  

      if(frazil_only_in_surface) then 
          k=1     
          do j=jsc,jec
             do i=isc,iec 
                T_diag(index_frazil_redist)%field(i,j,k) = 0.0
                if(Grd%tmask(i,j,k) > 0.0) then 
                    s       = T_prog(index_salt)%field(i,j,k,taup1)              
                    sqrts   = sqrt(s)
                    tf_num  = a0 + s*(a1 + sqrts*(a2 + sqrts*a3))
                    tf_den  = b0 + s*s*sqrts*b3 
                    tfreeze = tf_num/tf_den  + (c1 + s*c2) 
                    if(T_prog(index_redist_heat)%field(i,j,k,taup1) < tfreeze) then 
                        T_diag(index_frazil_redist)%field(i,j,k) = &
                             (tfreeze-T_prog(index_redist_heat)%field(i,j,k,taup1)) &
                             *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor 
                        T_prog(index_redist_heat)%field(i,j,k,taup1) = tfreeze
                    endif
                endif
             enddo
          enddo

      else 
#if defined(ACCESS)
       do j=jsc,jec
             do i=isc,iec
                do k=1,nk
#else
          do k=1,nk 
             do j=jsc,jec
                do i=isc,iec 
#endif
                   T_diag(index_frazil_redist)%field(i,j,k) = 0.0
                   if(Grd%tmask(i,j,k) > 0.0) then 
                       s       = Dens%rho_salinity(i,j,k,taup1)
                       press   = Dens%pressure_at_depth(i,j,k)              
                       sqrts   = sqrt(s)
                       tf_num  = a0 + s*(a1 + sqrts*(a2 + sqrts*a3)) + press*(a4 + press*(a5 + s*a6)) 
                       tf_den  = b0 + press*(b1 + press*b2) + s*s*sqrts*b3 
                       tfreeze = tf_num/tf_den + (c1 + s*c2) 
                       if(T_prog(index_redist_heat)%field(i,j,k,taup1) < tfreeze) then 
                           T_diag(index_frazil_redist)%field(i,j,k) = &
                                (tfreeze-T_prog(index_redist_heat)%field(i,j,k,taup1)) &
                                *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor 
                           T_prog(index_redist_heat)%field(i,j,k,taup1) = tfreeze
                       endif
                   endif
                enddo
             enddo
          enddo

      endif

  else  ! teos10 

      if(frazil_only_in_surface) then 
          k=1
          do j=jsc,jec
             do i=isc,iec
                T_diag(index_frazil_redist)%field(i,j,k) = 0.0
                if(Grd%tmask(i,j,k) > 0.0) then
                    s       = 1.d-2*Dens%rho_salinity(i,j,k,taup1)
                    sqrts   = sqrt(s)
                    tfreeze =  v0 &
                               + s*(v1 + sqrts*(v2 + sqrts*(v3 + sqrts*(v4 + sqrts*(v5 + v6*sqrts)))))  &
                               - saturation_fraction*(1e-3)*(2.4 - c1*s)*(1 + c2*(1.d0 - s/35.16504d0))
                    if(T_prog(index_redist_heat)%field(i,j,k,taup1) < tfreeze) then
                        T_diag(index_frazil_redist)%field(i,j,k) = &
                             (tfreeze-T_prog(index_redist_heat)%field(i,j,k,taup1)) &
                             *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor
                        T_prog(index_redist_heat)%field(i,j,k,taup1) = tfreeze
                    endif
                endif
             enddo
          enddo

      else

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  T_diag(index_frazil_redist)%field(i,j,k) = 0.0
                  if(Grd%tmask(i,j,k) > 0.0) then
                     s       = 1.d-2*Dens%rho_salinity(i,j,k,taup1)
                     sqrts   = sqrt(s)
                     press   = Dens%pressure_at_depth(i,j,k)              
                     tfreeze =  v0 &
                                + s*(v1 + sqrts*(v2 + sqrts*(v3 + sqrts*(v4 + sqrts*(v5 + v6*sqrts)))))         &
                                + press*(v7 + press*(v8 + v9*press))                                            &
                                + s*press*(v10 + press*(v12 + press*(v15 + v21*s)) + s*(v13 + v17*press+ v19*s) &
                                + sqrts*(v11 + press*(v14 + v18*press)  + s*(v16 + v20*press+ v22*s)))          &
                                - saturation_fraction*(1e-3)*(2.4 - c1*s)*(1 + c2*(1d0 - s/35.16504d0))
                     if(T_prog(index_redist_heat)%field(i,j,k,taup1) < tfreeze) then
                         T_diag(index_frazil_redist)%field(i,j,k) = &
                              (tfreeze-T_prog(index_redist_heat)%field(i,j,k,taup1)) &
                              *Thickness%rho_dzt(i,j,k,taup1)*cp_ocean*frazil_factor
                         T_prog(index_redist_heat)%field(i,j,k,taup1) = tfreeze
                     endif
                  endif
               enddo
            enddo
         enddo
     endif

  endif  ! endif for freezing_temp_simple


  if (id_frazil_redist_2d > 0) then
     call diagnose_2d(Time, Grd, id_frazil_redist_2d, T_diag(index_frazil_redist)%field(:,:,1)*dtimer)
  endif
  if (id_frazil_redist_3d > 0) then
     call diagnose_3d(Time, Grd, id_frazil_redist_3d, T_diag(index_frazil_redist)%field(:,:,:)*dtimer)
  endif
 

end subroutine compute_frazil_redist_heating
! </SUBROUTINE> NAME="compute_frazil_redist_heating">


     
end module ocean_frazil_mod
