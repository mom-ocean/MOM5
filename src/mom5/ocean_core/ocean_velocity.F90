module ocean_velocity_mod
#define COMP isc:iec,jsc:jec
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! S.M. Griffies 
! </CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</CONTACT>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! M.J. Harrison 
! </REVIEWER>
!
!<OVERVIEW>
! Time step the velocity field. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module steps the velocity field forward in time.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Durran, Numerical Methods for Wave Equations in Geophysical
! Fluid Dynamics (1999).
! </REFERENCE>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies, The MOM3 Manual (1999).
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and 
! A. Rosati, A Technical Guide to MOM4 (2004).
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Fundamentals of Ocean Climate Models (2004).
! Princeton University Press.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies (2012), Elements of MOM
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_velocity_nml">
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.  
!  </DATA> 
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="update_velocity_via_uprime" TYPE="logical">
!  When updating the velocity, this method first computes uprime
!  as the updated velocity minus the barotropic pressure gradient.
!  This approach is motivated from the rigid lid approach, in which 
!  the surface pressure was never used to update the barotropic
!  fields.  
!  With the explicit free surface, we have the choice to 
!  update the full velocity field, with the barotropic contributions
!  to the pressure field resulting from a time average in the 
!  external mode algorithm.  This approach is for testing only,
!  and it has been found to be unstable for many cases.  
!
!  update_velocity_via_uprime=.true. uses the older aproach, 
!  in which the udrho,vdrho fields are taken from the external
!  mode module. 
! 
!  update_velocity_via_uprime=.false. only takes the time averaged 
!  pressure from the external mode, and thus updates the full 
!  velocity and so recomputes the udrho,vdrho fields.  
!
!  Default update_velocity_via_uprime=.true.  
!  The case of update_velocity_via_uprime=.false. is for testing only.
!  It is not supported for general use.  
!  </DATA> 
! 
!  <DATA NAME="use_constant_velocity" TYPE="logical">
!  For running with time independent constant velocity.
!  For use with idealized cases. Default=.false. 
!  </DATA> 
!  <DATA NAME="constant_u" TYPE="real" UNITS="meter/sec" >
!  For running with use_constant_velocity. 
!  Set the i-velocity component to this value.  
!  Default constant_u=0.0 
!  </DATA> 
!  <DATA NAME="constant_v" TYPE="real" UNITS="meter/sec" >
!  For running with use_constant_velocity. 
!  Set the j-velocity component to this value.  
!  Default constant_v=0.0 
!  </DATA> 
!
!  <DATA NAME="zero_tendency" TYPE="logical">
!  For debugging. Will freeze the baroclinic velocity  fields. 
!  </DATA> 
!  <DATA NAME="zero_tendency_explicit_a" TYPE="logical">
!  For debugging. Will not use explicit-a part of the tendency. 
!  </DATA> 
!  <DATA NAME="zero_tendency_explicit_b" TYPE="logical">
!  For debugging. Will not use explicit-b part of the tendency. 
!  </DATA> 
!  <DATA NAME="zero_tendency_implicit" TYPE="logical">
!  For debugging. Will not use implicit part of the tendency. 
!  </DATA> 
!
!  <DATA NAME="truncate_velocity" TYPE="logical">
!  Truncate the baroclinic velocity to a maximum value.  Useful for cases where
!  the initial spin-up initiates spuriously large model velocities that would
!  otherwise cause the model to blow-up. Also can be used as a very simple 
!  "polar filter" in cases where have spherical coordinates and wish to avoid 
!  using the traditional polar filters.  
!  </DATA> 
!  <DATA NAME="truncate_velocity_value" UNITS="meter/sec" TYPE="real">
!  Speed above which will truncate the baroclinic velocity
!  </DATA> 
!  <DATA NAME="truncate_velocity_lat" UNITS="dimensionless" TYPE="real">
!  Latitude poleward of which we truncate the velocity. Useful in cases
!  when wish to truncate the velocity only in polar regions. Default is 0.0   
!  </DATA> 
!  <DATA NAME="truncate_verbose" TYPE="logical">
!  For verbose printout 
!  </DATA> 
!
!  <DATA NAME="max_cgint" TYPE="real">
!  Maximum internal gravity wave speed--used for diagnosing conservative
!  estimate of stable time steps.  
!  </DATA>
!
!  <DATA NAME="adams_bashforth_epsilon" UNITS="dimensionless" TYPE="real">
!  Dimensionless parameter for 2nd order Adams-Bashforth implementation of 
!  velocity advection.  Values between 0.5 and 1.0 are recommended.  
!  Value of 0.5 leads to second order accurate, but it is formally 
!  weakly unstable (Durran, Section 2.3.4).
!  </DATA> 
!  <DATA NAME="adams_bashforth_third" TYPE="logical">
!  For a third order treatment of the velocity advection.  
!  This is stable and so needs no temporal dissipation 
!  (Section 2.3.6 of Durran).  This is the model default. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: write_version_number, check_nml_error, close_file, open_namelist_file
use fms_mod,          only: file_exist, FATAL, WARNING, NOTE
use fms_io_mod,       only: field_size, restart_file_type, register_restart_field
use fms_io_mod,       only: save_restart, restore_state, reset_field_pointer
use mpp_domains_mod,  only: mpp_update_domains, mpp_global_sum, BGRID_NE, CGRID_NE, BITWISE_EXACT_SUM
use mpp_mod,          only: input_nml_file, mpp_error, mpp_pe, mpp_min, mpp_max, mpp_sum
use mpp_mod,          only: mpp_broadcast, stdout, stdlog
use time_interp_external_mod, only: time_interp_external, init_external_field

use ocean_bih_friction_mod,    only: bih_friction
use ocean_coriolis_mod,        only: coriolis_force_bgrid, coriolis_force_bgrid_implicit, coriolis_force_cgrid
use ocean_domains_mod,         only: get_local_indices
use ocean_form_drag_mod,       only: form_drag_accel
use ocean_lap_friction_mod,    only: lap_friction
use ocean_momentum_source_mod, only: momentum_source
use ocean_obc_mod,             only: ocean_obc_update_boundary
use ocean_operators_mod,       only: DIV_UD
use ocean_parameters_mod,      only: TWO_LEVEL, THREE_LEVEL, DEPTH_BASED
use ocean_parameters_mod,      only: missing_value, rho0r, grav, onehalf 
use ocean_parameters_mod,      only: MOM_BGRID, MOM_CGRID 
use ocean_pressure_mod,        only: pressure_force
use ocean_types_mod,           only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,           only: ocean_density_type, ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,           only: ocean_adv_vel_type, ocean_external_mode_type, ocean_options_type
use ocean_types_mod,           only: ocean_lagrangian_type
use ocean_util_mod,            only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_2d_u, diagnose_3d_u
use ocean_util_mod,            only: write_chksum_2d, write_chksum_3d
use ocean_velocity_advect_mod, only: horz_advection_of_velocity, vert_advection_of_velocity
use ocean_velocity_diag_mod,   only: kinetic_energy, potential_energy 
use ocean_vert_mix_mod,        only: vert_friction_bgrid, vert_friction_implicit_bgrid
use ocean_vert_mix_mod,        only: vert_friction_cgrid, vert_friction_implicit_cgrid
use ocean_workspace_mod,       only: wrk1, wrk2, wrk3, wrk1_v  

implicit none

private

! for diagnostics 
integer :: id_u(2)              =-1
integer :: id_u_on_depth(2)     =-1
integer :: id_usurf(2)          =-1
integer :: id_ubott(2)          =-1
integer :: id_speed             =-1
integer :: id_accel(2)          =-1
integer :: id_converge_rho_ud_t =-1
logical :: used

! for data over ride files
logical :: use_velocity_override  = .false.
integer :: id_velocity_u_override = -1
integer :: id_velocity_v_override = -1

!for restart file
integer                       :: id_restart_u = 0
integer                       :: id_restart_v = 0
integer                       :: id_restart_advu = 0
integer                       :: id_restart_advv = 0
integer                       :: id_restart_coru = 0
integer                       :: id_restart_corv = 0
integer                       :: id_restart_px = 0
integer                       :: id_restart_py = 0
type(restart_file_type), save :: Vel_restart
type(restart_file_type), save :: Adv_restart
type(restart_file_type), save :: Cor_restart

integer :: unit=6

character(len=128) :: &
     version='$Id: ocean_velocity.F90,v 20.0 2013/12/14 00:12:41 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: have_obc = .false.

#include <ocean_memory.h>

real    :: dtime 
integer :: time_index 
integer :: tendency

! for Bgrid or Cgrid 
integer :: horz_grid 

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_velocity_init
public ocean_velocity_end
public ocean_explicit_accel_a
public ocean_explicit_accel_b
public ocean_implicit_accel
public update_ocean_velocity_bgrid
public update_ocean_velocity_cgrid
public ocean_velocity_restart

private ocean_velocity_chksum
private check_gravity_wave_cfl
private velocity_truncate
private remap_s_to_depth

! Adams-Bashforth treatment of velocity advection 
logical :: adams_bashforth_third   = .true.
real    :: adams_bashforth_epsilon = 0.6
real    :: abtau_m0 =  23.0/12.0
real    :: abtau_m1 = -16.0/12.0
real    :: abtau_m2 =  5.0/12.0

logical :: update_velocity_via_uprime = .true. 
logical :: module_is_initialized      = .false.
logical :: zero_tendency              = .false.
logical :: zero_tendency_explicit_a   = .false.
logical :: zero_tendency_explicit_b   = .false.
logical :: zero_tendency_implicit     = .false.
logical :: truncate_velocity          = .false.
logical :: truncate_verbose           = .false.
real    :: truncate_velocity_lat      = 0.0
real    :: truncate_velocity_value    = 1.0

! for idealized constant velocity situations 
logical :: use_constant_velocity = .false.
real    :: constant_u            = 0.0
real    :: constant_v            = 0.0

! for CFL check on gravity waves 
real    :: max_cgint = 2.0   ! speed (m/s) used to compute CFL check for internal gravity waves  
real    :: max_dt_for_cgint  ! maximum time step (s) allowed by CFL to resolve internal gravity waves 

! for debugging 
logical :: debug_this_module = .false.

! for writing restarts 
logical :: write_a_restart   = .true. 

namelist /ocean_velocity_nml/  debug_this_module, write_a_restart, max_cgint,                       &
         zero_tendency, zero_tendency_explicit_a, zero_tendency_explicit_b, zero_tendency_implicit, &
         truncate_velocity, truncate_verbose, truncate_velocity_lat, truncate_velocity_value,       &
         adams_bashforth_third, adams_bashforth_epsilon,                                            &
         use_constant_velocity, constant_u, constant_v,                                             &
         update_velocity_via_uprime 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_init">
!
! <DESCRIPTION>
! Initialize terms for the velocity equation. 
! </DESCRIPTION>
!
subroutine ocean_velocity_init (Grid, Domain, Time, Time_steps, Ocean_options, &
                                Velocity, hor_grid, obc, use_blobs, introduce_blobs, &
                                velocity_override, debug)

  type(ocean_grid_type),   target, intent(in)    :: Grid
  type(ocean_domain_type), target, intent(in)    :: Domain
  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_time_steps_type),     intent(in)    :: Time_steps 
  type(ocean_options_type),        intent(inout) :: Ocean_options
  type(ocean_velocity_type),       intent(inout) :: Velocity
  integer,                         intent(in)    :: hor_grid 
  logical,                         intent(in)    :: use_blobs
  logical,                         intent(in)    :: introduce_blobs
  logical,                         intent(in)    :: obc
  logical,                         intent(in)    :: velocity_override
  logical,               optional, intent(in)    :: debug

  integer               :: ioun, io_status, ierr
  integer               :: m, n
  integer               :: taum1, tau, taup1
  integer               :: tau_m0, tau_m1
  integer, dimension(4) :: siz 
  character(len=128)    :: filename
  character(len=128)    :: velfilename

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_velocity_mod (ocean_velocity_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  have_obc  = obc
  tendency  = Time_steps%tendency
  dtime     = Time_steps%dtime_u
  horz_grid = hor_grid 

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_velocity_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_velocity_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_velocity_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_velocity_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_velocity_nml)  
  write (stdlogunit, ocean_velocity_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_velocity with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif 

  Dom => Domain
  Grd => Grid
  use_velocity_override = velocity_override


#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grid%nk
  allocate (Velocity%u(isd:ied,jsd:jed,nk,2,3))
  allocate (Velocity%smf(isd:ied,jsd:jed,2))
  allocate (Velocity%smf_bgrid(isd:ied,jsd:jed,2))
  allocate (Velocity%smf_cgrid(isd:ied,jsd:jed,2))
  allocate (Velocity%bmf(isd:ied,jsd:jed,2))
  allocate (Velocity%gamma(isd:ied,jsd:jed))
  allocate (Velocity%langmuirfactor(isd:ied,jsd:jed))
  allocate (Velocity%ustoke(isd:ied,jsd:jed))
  allocate (Velocity%vstoke(isd:ied,jsd:jed))
  allocate (Velocity%wavlen(isd:ied,jsd:jed))
  allocate (Velocity%cdbot_array(isd:ied,jsd:jed))
  allocate (Velocity%rossby_radius(isd:ied,jsd:jed))
  allocate (Velocity%stokes_depth(isd:ied,jsd:jed))
  allocate (Velocity%stokes_drift(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%press_force(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%accel(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%source(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%wrkv(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%advection(isd:ied,jsd:jed,nk,2,3))
  allocate (Velocity%coriolis(isd:ied,jsd:jed,nk,2,3))
  allocate (Velocity%vfrict_impl(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%lap_friction_bt(isd:ied,jsd:jed,2))
  allocate (Velocity%bih_friction_bt(isd:ied,jsd:jed,2))
  allocate (Velocity%current_wave_stress(isd:ied,jsd:jed))
#endif

  Velocity%smf                 = 0.0
  Velocity%smf_bgrid           = 0.0
  Velocity%smf_cgrid           = 0.0
  Velocity%bmf                 = 0.0
  Velocity%gamma               = 0.0
  Velocity%langmuirfactor      = 1.0
  Velocity%ustoke              = 0.0
  Velocity%vstoke              = 0.0
  Velocity%wavlen              = 1.0 
  Velocity%cdbot_array         = 0.0
  Velocity%rossby_radius       = 1.e5
  Velocity%stokes_depth        = 0.0
  Velocity%stokes_drift        = 0.0
  Velocity%accel               = 0.0
  Velocity%source              = 0.0
  Velocity%wrkv                = 0.0
  Velocity%press_force         = 0.0
  Velocity%advection           = 0.0
  Velocity%coriolis            = 0.0
  Velocity%vfrict_impl         = 0.0
  Velocity%lap_friction_bt     = 0.0
  Velocity%bih_friction_bt     = 0.0
  Velocity%current_wave_stress = 0.0

  ! register fields for diagnostic output

  id_u(1)    = register_diag_field ('ocean_model', 'u', Grd%vel_axes_u(1:3), Time%model_time, &
     'i-current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/),                 &
     standard_name='sea_water_x_velocity')
  id_u(2)    = register_diag_field ('ocean_model', 'v', Grd%vel_axes_v(1:3), Time%model_time, &
     'j-current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/),                 &
     standard_name='sea_water_y_velocity')

  id_u_on_depth(1) = register_diag_field ('ocean_model', 'u_on_depth', Grd%vel_axes_u_depth(1:3), Time%model_time, &
     'i-current mapped to depth surface', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))
  id_u_on_depth(2) = register_diag_field ('ocean_model', 'v_on_depth', Grd%vel_axes_v_depth(1:3), Time%model_time, &
     'j-current mapped to depth surface', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))

  id_usurf(1)  = register_diag_field ('ocean_model', 'usurf', Grd%vel_axes_u(1:2), Time%model_time, &
     'i-surface current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))
  id_usurf(2)  = register_diag_field ('ocean_model', 'vsurf', Grd%vel_axes_v(1:2), Time%model_time, &
     'j-surface current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))

  id_ubott(1)  = register_diag_field ('ocean_model', 'ubott', Grd%vel_axes_u(1:2), Time%model_time, &
     'i-bottom current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))
  id_ubott(2)  = register_diag_field ('ocean_model', 'vbott', Grd%vel_axes_v(1:2), Time%model_time, &
     'j-bottom current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))

  id_accel(1) = register_diag_field ('ocean_model', 'accel_u', Grd%vel_axes_u(1:3), Time%model_time, &
     'baroclinic forcing of rho*u', '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e9,1e9/))

  id_accel(2) = register_diag_field ('ocean_model', 'accel_v', Grd%vel_axes_v(1:3), Time%model_time, &
     'baroclinic forcing of rho*v', '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e9,1e9/))

  id_converge_rho_ud_t = register_diag_field ('ocean_model', 'converge_rho_ud_t', Grd%tracer_axes(1:2), &
                     Time%model_time,'convergence rho*ud on T cells computed in velocity_mod',          &
                     '(kg/m^3)*(m/s)', missing_value=missing_value, range=(/-10.0,10.0/))

  ! smg: following diags need some work 

  id_speed = register_diag_field ('ocean_model', 'speed', Grd%vel_axes_uv(1:3), Time%model_time, &
     'speed of horizontal current', 'm/sec', missing_value=missing_value, range=(/-10.0,10.0/))


  if (truncate_velocity) then
      write (stdoutunit,'(/a)')&
      '==>Warning: truncate_velocity=.true. Model will truncate baroclinc velocity'
      write (stdoutunit,'(a,f5.3)')&
      '   so that max horz velocity magnitude is (m/s) ',truncate_velocity_value 
      write (stdoutunit,'(a,f8.3)') &
      '   Truncation occurs for regions poleward of the latitude',truncate_velocity_lat 
      write (stdoutunit,'(a)')&
      '   This option may be useful during the spin-up phase of an experiment.'
      write (stdoutunit,'(a/)')&
      '   It is also of use for polar filtering. '
  endif

  if(tendency==THREE_LEVEL) then 
     abtau_m0= 1.0
     abtau_m1= 0.0
     abtau_m2= 0.0
  endif 

  if(tendency==TWO_LEVEL) then 
     write(stdoutunit,'(/a)')'==>Note from ocean_velocity_mod: use of twolevel time_tendency' 
     write(stdoutunit,'(a)') '   necessitates an Adams-Bashforth treatment of velocity advection.'

     if(adams_bashforth_third) then 
       write(stdoutunit,'(/a)')'   Using 3rd order Adams-Bashforth for velocity advection.'
       write(stdoutunit,'(a/)')'   This is the MOM default. '
       abtau_m0= 23.0/12.0
       abtau_m1=-16.0/12.0
       abtau_m2= 5.0/12.0
     else 
       write(stdoutunit,'(/a,f8.3/)')' Using 2nd order Adams-Bashforth with dimensionless AB parameter = ', &
                                     adams_bashforth_epsilon 
       abtau_m0= 1.0 + adams_bashforth_epsilon
       abtau_m1= 0.0 - adams_bashforth_epsilon
       abtau_m2= 0.0
       if(adams_bashforth_epsilon < 0.5 .or. adams_bashforth_epsilon > 1.0) then 
        call mpp_error(FATAL,&
        '==>Error from ocean_velocity_mod: adams_bashforth parameter must be between 0.5 & 1.0')
       endif 
     endif  

  endif 

  if(zero_tendency) then
     call mpp_error(NOTE, &
     '==>zero_tendency=.true. so will not time step the baroclinic velocity field ')
     Ocean_options%baroclinic_tendency = 'Did NOT time step the baroclinic velocity.'
  else
     Ocean_options%baroclinic_tendency = 'Time stepped the baroclinic velocity.'
  endif 
  if(zero_tendency_explicit_a) then
     call mpp_error(NOTE, &
     '==>zero_tendency_explicit_a=.true. so will not include explicit_a velocity tendency ')
  endif 
  if(zero_tendency_explicit_b) then
     call mpp_error(NOTE, &
     '==>zero_tendency_explicit_b=.true. so will not include explicit_b velocity tendency ')
  endif 
  if(zero_tendency_implicit) then
     call mpp_error(NOTE, &
     '==>zero_tendency_implicit=.true. so will not include implicit part of the velocity tendency ')
  endif 

  if(update_velocity_via_uprime) then 
     call mpp_error(NOTE, &
     '==>update_velocity_via_uprime=.true., so keep udrho from external mode solver.')
  else
     call mpp_error(WARNING, &
     '==>update_velocity_via_uprime=.false., so recompute udrho after updating full velocity. For testing only!')
  endif 
  if(horz_grid == MOM_CGRID .and. .not. update_velocity_via_uprime) then 
     call mpp_error(WARNING, &
     '==>update_velocity_via_uprime=.false. is not supported for Cgrid.  Resetting update_velocity_via_uprime to true.')
  endif 

  if(use_velocity_override) then
      
      call mpp_error(NOTE, &
           '==>use_velocity_override=.true. so will over-ride prognosed (u,v)-velocity with values from a file.')
      write(stdoutunit,'(a)') '==>ocean_velocity_mod: override applied to (u,v)-velocity.'
      
      filename = 'INPUT/velocity_u_override.nc'
      if (file_exist(trim(filename))) then
          id_velocity_u_override = init_external_field(filename, "u", domain=Dom%domain2d)
      endif
      if (id_velocity_u_override == -1) then 
          call mpp_error(FATAL,'==>Error in ocean_velocity_mod: failed to find u field in INPUT/velocity_u_override.nc') 
      endif     
      
      filename = 'INPUT/velocity_v_override.nc'
      if (file_exist(trim(filename))) then
          id_velocity_v_override = init_external_field(filename, "v", domain=Dom%domain2d)
      endif 
      if (id_velocity_v_override == -1) then 
          call mpp_error(FATAL,'==>Error in ocean_velocity_mod: failed to find v field in INPUT/velocity_v_override.nc') 
      endif

      Ocean_options%override_velocity = 'Used velocity override for (u,v) and (udrho,vdrho) fields.' 

  else

      Ocean_options%override_velocity = 'Did NOT use velocity override'
      
  endif 

  if(use_constant_velocity) then 
     call mpp_error(NOTE, &
     '==>use_constant_velocity=.true. so will hold the velocity to a constant.')
     zero_tendency=.true.
     Ocean_options%baroclinic_tendency = 'Did NOT time step velocity; kept it = constant.'
     write(stdoutunit,'(a,f6.2,a,f6.2,a)') &
     '==>Note from ocean_velocity_init: velocity fixed at (u,v)(m/s) = (',constant_u,',',constant_v,')'
  endif 

  call check_gravity_wave_cfl()
  
  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0

  ! initialize to small but nonzero velocity 
  do m=1,3 
     do n=1,2  
        Velocity%u(:,:,:,n,m) = epsln*Grd%umask(:,:,:)
     enddo
  enddo

  if(use_constant_velocity) then 
      do m=1,3 
         Velocity%u(:,:,:,1,m) = constant_u*Grd%umask(:,:,:)
         Velocity%u(:,:,:,2,m) = constant_v*Grd%umask(:,:,:)
      enddo
  endif

  filename = 'ocean_velocity.res.nc'
  velfilename = filename
  if(tendency==THREE_LEVEL) then
     id_restart_u = register_restart_field(Vel_restart, filename,'u',Velocity%u(:,:,:,1,tau), &
          Velocity%u(:,:,:,1,taup1), Dom%domain2d)
     id_restart_v = register_restart_field(Vel_restart, filename,'v',Velocity%u(:,:,:,2,tau), &
          Velocity%u(:,:,:,2,taup1), Dom%domain2d )
  elseif (tendency==TWO_LEVEL) then
     id_restart_u = register_restart_field(Vel_restart, filename,'u',Velocity%u(:,:,:,1,taup1), &
          Dom%domain2d )
     id_restart_v = register_restart_field(Vel_restart, filename,'v',Velocity%u(:,:,:,2,taup1), &
          Dom%domain2d )
     if (use_blobs .and. .not. introduce_blobs) then
        id_restart_px = register_restart_field(Vel_restart, filename,'pressurex',Velocity%press_force(:,:,:,1),&
             Dom%domain2d )
        id_restart_py = register_restart_field(Vel_restart, filename,'pressurey',Velocity%press_force(:,:,:,2),&
             Dom%domain2d )        
     endif
     filename = 'ocean_velocity_advection.res.nc'
     id_restart_advu = register_restart_field(Adv_restart, filename,'advectionu',Velocity%advection(:,:,:,1,tau_m1),&
          Velocity%advection(:,:,:,1,tau_m0), Dom%domain2d )
     id_restart_advv = register_restart_field(Adv_restart, filename,'advectionv',Velocity%advection(:,:,:,2,tau_m1),&
          Velocity%advection(:,:,:,2,tau_m0), Dom%domain2d )
  end if

  if(horz_grid == MOM_CGRID) then 
     filename = 'ocean_velocity_coriolis.res.nc'
     id_restart_coru = register_restart_field(Cor_restart, filename,'coriolisu',Velocity%coriolis(:,:,:,1,tau_m1),&
          Velocity%coriolis(:,:,:,1,tau_m0), Dom%domain2d )
     id_restart_corv = register_restart_field(Cor_restart, filename,'coriolisv',Velocity%coriolis(:,:,:,2,tau_m1),&
          Velocity%coriolis(:,:,:,2,tau_m0), Dom%domain2d )
  endif 


  filename = 'INPUT/ocean_velocity.res.nc'
  if (file_exist(trim(filename))) then 

      !--- read restart file data by calling restore_state
      call restore_state(Vel_restart)

      if(tendency==THREE_LEVEL) then

          write (stdoutunit,'(a)') '  Reading THREE_LEVEL restart for velocity from INPUT/ocean_velocity.res.nc'
          write (stdoutunit,'(a)') '  Expecting two time records for each restart field.'
          call mpp_update_domains(Velocity%u(:,:,:,1,tau), Velocity%u(:,:,:,2,tau),&
               Dom%domain2d,gridtype=BGRID_NE)
          call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1),&
               Dom%domain2d,gridtype=BGRID_NE)

          if(have_obc) then
            call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','n')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','t')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','t')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','n')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,1,tau), 'M','n')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,2,tau), 'M','t')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,1,tau), 'Z','t')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,2,tau), 'Z','n')
          endif

      elseif (tendency==TWO_LEVEL) then

          write (stdoutunit,'(/a)') '  Reading TWO_LEVEL restart for velocity from INPUT/ocean_velocity.res.nc'
          write (stdoutunit,'(a)')  '  Expecting only one time record for each restart field.'

          call field_size(filename,'u', siz)
          if (siz(4) > 1) then
            write(stdoutunit,'(/a)') &
            '==>ERROR: Attempt to read ocean_velocity.res.nc from a 3-level time scheme (2 time records)'
            write(stdoutunit,'(a)')  &
            '          when running MOM with 2-level timestepping (only need 1 time record in restart).'
            write(stdoutunit,'(a)')  &
            '          Reduce restart file to a single time record in order to avoid confusion.'
            call mpp_error(FATAL,  &
            'Reading 3-time lev ocean_velocity.res.nc (w/ 2 time records) while using 2-lev (needs only 1 record)')
          endif

          if(horz_grid == MOM_BGRID) then 
             call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1),&
                  Dom%domain2d,gridtype=BGRID_NE)
             if(use_blobs .and. .not. introduce_blobs) then
                call mpp_update_domains(Velocity%press_force(:,:,:,1), Velocity%press_force(:,:,:,2),&
                  Dom%domain2d,gridtype=BGRID_NE)
             elseif(use_blobs .and. introduce_blobs) then
                ! If we are introducing blobs, pressurex and pressurey will not have been in the
                ! restart file, but, they will be written to the restart file at the end of this run.
                id_restart_px = register_restart_field(Vel_restart, velfilename,'pressurex',&
                     Velocity%press_force(:,:,:,1),Dom%domain2d )
                id_restart_py = register_restart_field(Vel_restart, velfilename,'pressurey',&
                     Velocity%press_force(:,:,:,2),Dom%domain2d )        
             endif
          else 
             call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1),&
                  Dom%domain2d,gridtype=CGRID_NE)
             if(use_blobs .and. .not. introduce_blobs) then
                call mpp_update_domains(Velocity%press_force(:,:,:,1), Velocity%press_force(:,:,:,2),&
                  Dom%domain2d,gridtype=CGRID_NE)
             elseif(use_blobs .and. introduce_blobs) then
                ! If we are introducing blobs, pressurex and pressurey will not have been in the
                ! restart file, but, they will be written to the restart file at the end of this run.
                id_restart_px = register_restart_field(Vel_restart, velfilename,'pressurex',&
                     Velocity%press_force(:,:,:,1),Dom%domain2d )
                id_restart_py = register_restart_field(Vel_restart, velfilename,'pressurey',&
                     Velocity%press_force(:,:,:,2),Dom%domain2d )        
             endif
          endif 


          write (stdoutunit,'(a)') '  Finished reading restart for velocity field'

          if(.not. Time%init) then  
              write (stdoutunit,'(/a)') '  Reading TWO_LEVEL restart from INPUT/ocean_velocity_advection.res.nc'
              write (stdoutunit,'(a)')  '  Expecting two time records for the advection velocity fields.'
              call restore_state(Adv_restart)
              write (stdoutunit,'(a)') '  Finished reading restart for velocity advection operator'

              call mpp_update_domains(Velocity%advection(:,:,:,1,tau_m1),Dom%domain2d)
              call mpp_update_domains(Velocity%advection(:,:,:,1,tau_m0),Dom%domain2d)
              call mpp_update_domains(Velocity%advection(:,:,:,2,tau_m1),Dom%domain2d)
              call mpp_update_domains(Velocity%advection(:,:,:,2,tau_m0),Dom%domain2d)

          endif

          if(have_obc) then
            call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','n')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','t')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','t')
            call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','n')
         endif
      endif

 else

     if (Time%init) then
          call mpp_error(NOTE,&
         '==> From ocean_velocity_mod: Initializing velocity to zero since Time%init=.true. &
          &and did not find INPUT/ocean_velocity.res.nc.')
     endif 

     if (.NOT. Time%init) then
          call mpp_error(FATAL,&
         '==> Error from ocean_velocity_mod: Expecting INPUT/ocean_velocity.res.nc to exist.&
          &This file was not found and Time%init=.false.')
     endif 

 endif

 ! Coriolis force saved for Cgrid 
 if(horz_grid == MOM_CGRID) then 

    filename = 'INPUT/ocean_velocity_coriolis.res.nc'
    if (file_exist(trim(filename))) then 

        !--- read restart file data by calling restore_state
        write (stdoutunit,'(/a)') '  Reading Coriolis restart from INPUT/ocean_velocity_coriolis.res.nc'
        write (stdoutunit,'(a)')  '  Expecting two time records for the Coriolis force.'
        call restore_state(Cor_restart)
        write (stdoutunit,'(a)') '  Finished reading restart for Coriolis force.'

        call mpp_update_domains(Velocity%coriolis(:,:,:,1,tau_m1),Dom%domain2d)
        call mpp_update_domains(Velocity%coriolis(:,:,:,1,tau_m0),Dom%domain2d)
        call mpp_update_domains(Velocity%coriolis(:,:,:,2,tau_m1),Dom%domain2d)
        call mpp_update_domains(Velocity%coriolis(:,:,:,2,tau_m0),Dom%domain2d)

    else 

       if (Time%init) then
          call mpp_error(NOTE,&
         '==> From ocean_velocity_mod: Initializing Coriolis to zero since Time%init=.true. &
          &and did not find INPUT/ocean_velocity_coriolis.res.nc.')
       endif 
       if (.NOT. Time%init) then
            call mpp_error(FATAL,&
           '==> Error from ocean_velocity_mod: Expecting INPUT/ocean_velocity_coriolis.res.nc to exist.&
            &This file was not found and Time%init=.false.')
       endif 

    endif 

 endif 

     write(stdoutunit,*)' ' 
     write(stdoutunit,*) '===Initial velocity checksums ==>'
     call write_timestamp(Time%model_time)
 if(tendency==THREE_LEVEL) then 
     write(stdoutunit,*) 'From ocean_velocity_mod: initial velocity chksum (tau)'
     call ocean_velocity_chksum(Velocity, tau, write_advection=.false.)
     write(stdoutunit,*) 'From ocean_velocity_mod: initial velocity chksum (taup1)'
     call ocean_velocity_chksum(Velocity, taup1, write_advection=.false.)
 else
     write(stdoutunit,*) 'From ocean_velocity_mod: initial velocity chksum (taup1)'
     call ocean_velocity_chksum(Velocity, taup1, write_advection=.true.)
 endif

end subroutine ocean_velocity_init
! </SUBROUTINE> NAME="ocean_velocity_init"



!#######################################################################
! <SUBROUTINE NAME="check_gravity_wave_cfl">
!
! <DESCRIPTION>
! Check CFL for internal gravity waves. 
! </DESCRIPTION>
!
 subroutine check_gravity_wave_cfl()

    logical :: cfl_error=.false.
    real    :: gridsp, dtcg, max_dt_for_cgint0
    real    :: cfl_grid_factor, dtuv
    integer :: i, j, icg, jcg
    integer :: stdlogunit
    stdlogunit=stdlog()
   
   ! estimate the maximum timestep allowable for resolving internal gravity waves
    if(tendency==TWO_LEVEL)   then 
       cfl_grid_factor=1.0
       dtuv = dtime 
    elseif(tendency==THREE_LEVEL) then
       cfl_grid_factor=0.5
       dtuv = 0.5*dtime
    endif 
    icg = isc; jcg = jsc; gridsp = 1.0e20; max_dt_for_cgint = 1.0e6
    do j=jsc,jec
       do i=isc,iec
          if (Grd%kmu(i,j) > 0) then
              gridsp = 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j))
              gridsp = sqrt(1.0/gridsp) 
              dtcg = cfl_grid_factor*gridsp/(epsln+max_cgint)                       
              if (dtcg < max_dt_for_cgint) then
                  max_dt_for_cgint = dtcg; icg  = i; jcg  = j
              endif
          endif
       enddo
    enddo
    max_dt_for_cgint = nint(max_dt_for_cgint)
    max_dt_for_cgint = max_dt_for_cgint + 0.001*mpp_pe() ! to separate redundancies
    max_dt_for_cgint0 = max_dt_for_cgint
    call mpp_min (max_dt_for_cgint)

    ! show the most unstable location for baroclinic gravity waves
    if (max_dt_for_cgint == max_dt_for_cgint0 .and. nint(dtuv) > 0) then
        if (dtuv <= max_dt_for_cgint) then
            write (unit,'(/a,i4,a,i4,a,f9.2,a,f9.2,a)') &
            ' Baroclinic time step stability most nearly violated at U-cell (i,j) = (',&
                 icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
            write(unit,'(a,i6)')    '         The number of kmu-levels  at this point is ',Grd%kmu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dxu grid distance (m) at this point is ',Grd%dxu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dyu grid distance (m) at this point is ',Grd%dyu(icg,jcg) 
            cfl_error=.false.
        else
            write (unit,'(/a,i4,a,i4,a,f9.2,a,f9.2,a)') &
            '==>Error: Baroclinic time step stability violated at U-cell (i,j) = (',&
                 icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
            write(unit,'(a,i6)')    '         The number of kmu-levels  at this point is ',Grd%kmu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dxu grid distance (m) at this point is ',Grd%dxu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dyu grid distance (m) at this point is ',Grd%dyu(icg,jcg) 
            cfl_error=.true.
        endif
        write (unit,'(a,f5.2,a/a,f6.0,a,f6.0,a)')&
             '         Due to a specified maximum baroclinic gravity wave speed of ',max_cgint,' m/s.',&
             '         "dtuv" must be less than ',max_dt_for_cgint,' sec. "dtuv" = ', dtuv,' sec.'
        if(cfl_error) then  
          call mpp_error(FATAL, &
          '==>Error: time step instability detected for baroclinic gravity waves in ocean_model_mod')
        endif 
        write (stdlogunit,'(/a/)') &
        ' Note: A more appropriate maximum baroclinic gravity wave speed can be specified via namelist.'
    endif
    max_dt_for_cgint = nint(max_dt_for_cgint)
    
  end subroutine check_gravity_wave_cfl
! </SUBROUTINE> NAME="check_gravity_wave_cfl"


 
!#######################################################################
! <SUBROUTINE NAME="ocean_explicit_accel_a">
!
! <DESCRIPTION>
!
! Time explicit contributions to thickness weighted and density 
! weighted acceleration.
! 
! Omit here the Coriolis force and vertical friction. 
! They are omitted in order to facilitate the construction of the 
! vertically integrated forcing of the barotropic dynamics.  
! They will be added later in ocean_explicit_accel_b. 
!
! </DESCRIPTION>
!
subroutine ocean_explicit_accel_a(Velocity, Time, Adv_vel, Thickness, Dens, &
                                  L_system, rho, pme, river, upme, uriver)

  type(ocean_velocity_type),      intent(inout) :: Velocity
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river
  real, dimension(isd:,jsd:,:),   intent(in)    :: upme
  real, dimension(isd:,jsd:,:),   intent(in)    :: uriver 

  integer :: i, j, k, n
  integer :: taum1, tau, taup1
  integer :: tau_m0, tau_m1, tau_m2

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  tau_m2 = Time%tau_m2
  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0

  ! initialize to zero 
  do n=1,2
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Velocity%accel(i,j,k,n)            = 0.0 
              Velocity%source(i,j,k,n)           = 0.0 
              Velocity%press_force(i,j,k,n)      = 0.0 
              Velocity%advection(i,j,k,n,tau_m0) = 0.0 
              Velocity%coriolis(i,j,k,n,tau_m0)  = 0.0 
           enddo
        enddo
     enddo
  enddo
  
  if (.not. zero_tendency .and. .not. zero_tendency_explicit_a) then 

      ! fill Velocity%advection(tau_m0) 
      call horz_advection_of_velocity(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.false.)       
      call vert_advection_of_velocity(Time, Adv_vel, Velocity, &
                                      pme, river, upme, uriver, energy_analysis_step=.false.)       

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = -abtau_m0*Velocity%advection(i,j,k,n,tau_m0) & 
                                            -abtau_m1*Velocity%advection(i,j,k,n,tau_m1) & 
                                            -abtau_m2*Velocity%advection(i,j,k,n,tau_m2)  
               enddo
            enddo
         enddo
      enddo

      ! include contribution from horizontal pressure gradient 
      call pressure_force(Time, Thickness, Dens, Velocity, L_system, rho(:,:,:)) 
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + Velocity%press_force(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

      ! compute horizontal friction and add to Velocity%accel 
      call lap_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.false.)
      call bih_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.false.)

      ! compute momentum source/sinks and add to Velocity%accel 
      call momentum_source(Time, Thickness, Velocity, energy_analysis_step=.false.)

      ! compute parameterized form drag and add to Velocity%accel 
      call form_drag_accel(Time, Thickness, Velocity, energy_analysis_step=.false.)

  endif

  if (debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_velocity_mod: acceleration chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('explicit accel_a(1)', Velocity%accel(COMP,:,1))
      call write_chksum_3d('explicit accel_a(2)', Velocity%accel(COMP,:,2))
  endif

end subroutine ocean_explicit_accel_a
! </SUBROUTINE> NAME="ocean_explicit_accel_a"

 
!#######################################################################
! <SUBROUTINE NAME="ocean_explicit_accel_b">
!
! <DESCRIPTION>
!
! Add the contributions from the Coriolis force, computed explicitly 
! in time, and those from explicit vertical friction.  Add these to
! the thickness weighted and density weighted acceleration. 
!
! Note: no visc_cbu_form_drag is included here, since it must
! be handled via implicit vertical friction to maintain stability
! for general cases. 
!
! </DESCRIPTION>
!
subroutine ocean_explicit_accel_b(visc_cbu, visc_cbt, Time, Thickness, Adv_vel, Velocity)

  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbt
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_velocity_type),      intent(inout) :: Velocity
  
  integer :: stdoutunit 
  integer :: tau_m0, tau_m1, tau_m2

  tau_m2 = Time%tau_m2
  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0
  stdoutunit=stdout() 

  if (.not. zero_tendency .and. .not. zero_tendency_explicit_b) then 

      if(horz_grid == MOM_BGRID) then    
         call vert_friction_bgrid(Time, Thickness, Velocity, visc_cbu, energy_analysis_step=.false.)
         call coriolis_force_bgrid(Time, Thickness, Velocity, energy_analysis_step=.false.)
      else 
         call vert_friction_cgrid(Time, Thickness, Velocity, visc_cbt, energy_analysis_step=.false.)
         call coriolis_force_cgrid(Time, Adv_vel, Velocity, abtau_m0, abtau_m1, abtau_m2, &
                                   energy_analysis_step=.false.)
      endif 

      if (debug_this_module) then
          write(stdoutunit,*) ' ' 
          write(stdoutunit,*) 'From ocean_velocity_mod: explicit acceleration chksums'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('accel_b(1)', Velocity%accel(COMP,:,1))
          call write_chksum_3d('accel_b(2)', Velocity%accel(COMP,:,2))
      endif

  endif


end subroutine ocean_explicit_accel_b
! </SUBROUTINE> NAME="ocean_explicit_accel_b"


!#######################################################################
! <SUBROUTINE NAME="ocean_implicit_accel">
!
! <DESCRIPTION>
!
! Add the time implicit contributions from the Coriolis force
! and vertical friction. Add these to the thickness weighted 
! and density weighted acceleration. 
!
! Note the contribution from visc_cbu_form_drag is used for 
! implicit vertical friction. 
!
! Note there is no time implicit Coriolis for Cgrid.
!
! </DESCRIPTION>
!
subroutine ocean_implicit_accel(visc_cbu, visc_cbt, visc_cbu_form_drag, Time, Thickness, Adv_vel, Velocity)

  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(in)    :: visc_cbu_form_drag
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_velocity_type),      intent(inout) :: Velocity
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. zero_tendency .and. .not. zero_tendency_implicit) then 

      if(horz_grid == MOM_BGRID) then    
         call vert_friction_implicit_bgrid(visc_cbu, visc_cbu_form_drag, Time, Thickness, Velocity)
         call coriolis_force_bgrid_implicit(Time, Velocity) 
      else 
         call vert_friction_implicit_cgrid(visc_cbt, visc_cbu_form_drag, Time, Thickness, Adv_vel, Velocity)
      endif 

      if (debug_this_module) then
          write(stdoutunit,*) ' ' 
          write(stdoutunit,*) 'From ocean_velocity_mod: acceleration chksums after implicit update'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('accel(1)', Velocity%accel(COMP,:,1))
          call write_chksum_3d('accel(2)', Velocity%accel(COMP,:,2))
      endif
  endif

  ! since implicit_accel is the last piece of accel computed, it is 
  ! time to now output the full acceleration to diagnostics manager. 
  if(horz_grid == MOM_BGRID) then 
     call diagnose_3d_u(Time, Grd, id_accel(1), Velocity%accel(:,:,:,1))
     call diagnose_3d_u(Time, Grd, id_accel(2), Velocity%accel(:,:,:,2))
  else
     call diagnose_3d(Time, Grd, id_accel(1), Velocity%accel(:,:,:,1))
     call diagnose_3d(Time, Grd, id_accel(2), Velocity%accel(:,:,:,2))
  endif


end subroutine ocean_implicit_accel
! </SUBROUTINE> NAME="ocean_implicit_accel"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_velocity_bgrid">
!
! <DESCRIPTION>
! Update velocity components.
! There are two general methods available.
! 
! 1/ update baroclinic velocity; then add (udrho,vdrho) from external 
! mode algorithm to get the full velocity field.  This method is 
! analogous to the older approach with the rigid lid. 
!
! 2/ update the full velocity, so there is no need to introduce the 
! intermediate step with the baroclinic velocity.  To remain stable,
! we must use the time filtered barotropic pressure gradient.  
! We then diagnose the vertically integrated horizontal momentum
! (udrho,vdrho) and its convergence, since these fields are needed
! elsewhere.  
! 
! </DESCRIPTION>
!
subroutine update_ocean_velocity_bgrid(Time, Thickness, barotropic_split, &
                                       vert_coordinate_class, Ext_mode, Velocity)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  integer,                        intent(in)    :: barotropic_split
  integer,                        intent(in)    :: vert_coordinate_class
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity

  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(isd:ied,jsd:jed) :: tmpu
  real, dimension(isd:ied,jsd:jed) :: tmpv
  integer :: i, j, k, kbot, n
  integer :: taum1, tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1
  tmp   = 0.0

  if (barotropic_split > 1 .and. update_velocity_via_uprime) then  

      do n=1,2

         ! compute uprime 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%u(i,j,k,n,taup1) = Thickness%rho_dzur(i,j,k)           &
                       *(Thickness%rho_dzu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1) &
                       +dtime*Velocity%accel(i,j,k,n)*Grd%umask(i,j,k) )                
               enddo
            enddo
         enddo

         ! sum uprime*rho_dzu
         tmp=0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec  
                  tmp(i,j) = tmp(i,j) &
                       +Velocity%u(i,j,k,n,taup1)*Thickness%rho_dzu(i,j,k,taup1)*Grd%umask(i,j,k)
               enddo
            enddo
         enddo

         ! (vertical mean) - (value from barotropic integration) 
         do j=jsc,jec
            do i=isc,iec  
               tmp(i,j) = Grd%umask(i,j,1)*(tmp(i,j)-Ext_mode%udrho(i,j,n,taup1)) &
                    /(Thickness%mass_u(i,j,taup1)+epsln) 
            enddo
         enddo

         ! replace vertical mean with value from barotropic integration to update velocity 
         do k=1,nk 
            do j=jsc,jec
               do i=isc,iec  
                  Velocity%u(i,j,k,n,taup1) = (Velocity%u(i,j,k,n,taup1)-tmp(i,j))*Grd%umask(i,j,k) 
               enddo
            enddo
         enddo

      enddo  ! finish the n=1,2 do-loop 


  ! update full velocity directly with the full pressure gradient.  
  ! barotropic contribution to pressure gradient is time averaged
  ! if barotropic_split>1.  Unless dtuv is very small, this approach
  ! can be rather unstable.  User beware! 
  else  

      do n=1,2

         if(vert_coordinate_class==DEPTH_BASED) then 

             ! include gradient of surface pressure 
             do k=1,nk 
                do j=jsc,jec
                   do i=isc,iec 
                      Velocity%u(i,j,k,n,taup1) =                                                    &
                      (Thickness%rho_dzu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1)                      &
                      +dtime*(Velocity%accel(i,j,k,n)                                                &
                      -rho0r*Thickness%rho_dzu(i,j,k,tau)*Ext_mode%grad_ps(i,j,n)*Grd%umask(i,j,k))) &
                      *Thickness%rho_dzur(i,j,k)
                   enddo
                enddo
             enddo

         else 
  
            ! include gradient of bottom pressure
             do k=1,nk 
                do j=jsc,jec
                   do i=isc,iec 
                      Velocity%u(i,j,k,n,taup1) =                                                        &
                      (Thickness%rho_dzu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1)                          &
                      +dtime*(Velocity%accel(i,j,k,n)                                                    &
                      -rho0r*Thickness%rho_dzu(i,j,k,tau)*Ext_mode%grad_anompb(i,j,n)*Grd%umask(i,j,k))) &
                      *Thickness%rho_dzur(i,j,k)
                   enddo
                enddo
             enddo

         endif

         ! diagnose vertically integrated horizontal momentum 
         Ext_mode%udrho(:,:,n,taup1) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Ext_mode%udrho(i,j,n,taup1) = Ext_mode%udrho(i,j,n,taup1) &
                       +Velocity%u(i,j,k,n,taup1)*Thickness%rho_dzu(i,j,k,taup1)
               enddo
            enddo
         enddo

      enddo  ! end for n=1,2

      ! fill halos for vertically integrated horizontal momentum
      if(horz_grid == MOM_BGRID) then 
         call mpp_update_domains (Ext_mode%udrho(:,:,1,taup1),Ext_mode%udrho(:,:,2,taup1), &
                                  Dom%domain2d, gridtype=BGRID_NE)
      else 
         call mpp_update_domains (Ext_mode%udrho(:,:,1,taup1),Ext_mode%udrho(:,:,2,taup1), &
                                  Dom%domain2d, gridtype=CGRID_NE)
      endif 

      ! construct convergence of udrho at (tau+1)
      Ext_mode%conv_rho_ud_t(:,:,taup1) = -DIV_UD(Ext_mode%udrho(:,:,:,taup1),1,1)
      call mpp_update_domains (Ext_mode%conv_rho_ud_t(:,:,taup1), Dom%domain2d)

  endif   ! endif for splitting or uprime method 


  ! override update of velocity in the case of zero tendency 
  if(zero_tendency) then 
     Velocity%u(isc:iec,jsc:jec,:,:,taup1) = Velocity%u(isc:iec,jsc:jec,:,:,tau)
  endif 

  ! override update of u-velocity with value from a file 
  if(use_velocity_override .and. id_velocity_u_override > 0) then 
     wrk1(:,:,:) = 0.0
     call time_interp_external(id_velocity_u_override, Time%model_time, wrk1)
     Velocity%u(isc:iec,jsc:jec,:,1,taup1) = wrk1(isc:iec,jsc:jec,:)*Grd%umask(isc:iec,jsc:jec,:)
  endif 

  ! override update of v-velocity with value from a file 
  if(use_velocity_override .and. id_velocity_v_override > 0) then 
     wrk1(:,:,:) = 0.0
     call time_interp_external(id_velocity_v_override, Time%model_time, wrk1)
     Velocity%u(isc:iec,jsc:jec,:,2,taup1) = wrk1(isc:iec,jsc:jec,:)*Grd%umask(isc:iec,jsc:jec,:)
  endif 


  ! for debugging 
  do n=1,2

     if(truncate_velocity) then 
         call velocity_truncate(n, taup1, Velocity)
     endif

     if (debug_this_module .and. n==2) then        
         write(stdoutunit,*) ' ' 
         write(stdoutunit,*) 'From update_ocean_velocity at the u-prime step'
         call write_timestamp(Time%model_time)
         call write_chksum_3d('accel(1)', Velocity%accel(COMP,:,1))
         call write_chksum_3d('accel(2)', Velocity%accel(COMP,:,2))
         call write_chksum_3d('rho_dzu(taup1)', Thickness%rho_dzu(COMP,:,taup1))
         call write_chksum_3d('rho_dzu(taum1)', Thickness%rho_dzu(COMP,:,taum1))
         call write_chksum_3d('rho_dzur', Thickness%rho_dzur(COMP,:))
         call write_chksum_2d('udrho(1)', Ext_mode%udrho(COMP,1,taup1))
         call write_chksum_2d('udrho(2)', Ext_mode%udrho(COMP,2,taup1))
         call write_chksum_2d('mass_u', Thickness%mass_u(COMP,taup1))
         call ocean_velocity_chksum(Velocity, taum1, write_advection=.false.)
         call ocean_velocity_chksum(Velocity, taup1, write_advection=.false.)
     endif

  enddo


  !--- update velocity at the global halo points to make the gradient 0 across boundary
  if(have_obc) then
     call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','n')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','t')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','t')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','n')
  endif

  if (debug_this_module) then        
     write(stdoutunit,*) ' ' 
     write(stdoutunit,*) 'From update_ocean_velocity: velocity(taup1) chksum'
     call write_timestamp(Time%model_time)
     call ocean_velocity_chksum(Velocity, taup1, write_advection=.false.)
  endif

  ! send diagnostics to diagnostics manager 
  call diagnose_3d_u(Time, Grd, id_u(1), Velocity%u(:,:,:,1,tau))
  call diagnose_3d_u(Time, Grd, id_u(2), Velocity%u(:,:,:,2,tau))

  call diagnose_2d_u(Time, Grd, id_usurf(1), Velocity%u(:,:,1,1,tau))
  call diagnose_2d_u(Time, Grd, id_usurf(2), Velocity%u(:,:,1,2,tau))

  call diagnose_2d(Time, Grd, id_converge_rho_ud_t, Ext_mode%conv_rho_ud_t(:,:,tau))

  if(id_speed > 0) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = sqrt(Velocity%u(i,j,k,1,tau)**2 + Velocity%u(i,j,k,2,tau)**2)
            enddo
         enddo
      enddo
      call diagnose_3d_u(Time, Grd, id_speed, wrk1(:,:,:))
  endif

  if(id_ubott(1) > 0 .or. id_ubott(2) > 0) then 
      tmpu = 0.0 ; tmpv = 0.0
      do j=jsd,jed
         do i=isd,ied
            kbot = Grd%kmu(i,j)
            if(kbot > 1) then 
                tmpu(i,j) = Velocity%u(i,j,kbot,1,tau)
                tmpv(i,j) = Velocity%u(i,j,kbot,2,tau)
            endif
         enddo
      enddo
      call diagnose_2d_u(Time, Grd, id_ubott(1), tmpu(:,:))
      call diagnose_2d_u(Time, Grd, id_ubott(2), tmpv(:,:))
  endif

  if(id_u_on_depth(1) > 0) then 
    call remap_s_to_depth(Thickness, Time, Velocity%u(:,:,:,1,tau), 1)
  endif 
  if(id_u_on_depth(2) > 0) then 
    call remap_s_to_depth(Thickness, Time, Velocity%u(:,:,:,2,tau), 2)
  endif 


end subroutine update_ocean_velocity_bgrid
! </SUBROUTINE> NAME="update_ocean_velocity_bgrid"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_velocity_cgrid">
!
! <DESCRIPTION>
! Update velocity components, assuming Cgrid layout.
!
! Also assume splitting, so update baroclinic velocity; 
! then add (udrho,vdrho) from external mode algorithm to get the
! full velocity field.  This method is 
! analogous to the older approach with the rigid lid. 
!
! As we have already called update_ucell_thickness, the array
! rho_dzten has been udpated to taup1.  We use Adv_vel to obtain 
! the tau value of (u,v)*rho_dzten, since Adv_vel has been computed
! at the start of the time stepping using (u,v)(tau).
!
! </DESCRIPTION>
!
subroutine update_ocean_velocity_cgrid(Time, Thickness, Adv_vel, Ext_mode, Velocity)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity

  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(isd:ied,jsd:jed) :: tmpu
  real, dimension(isd:ied,jsd:jed) :: tmpv
  real :: rho_dz_vel_tau
  integer :: i, j, k, kbot, n
  integer :: taum1, tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1
  tmp   = 0.0

  do n=1,2

     ! compute uprime 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              rho_dz_vel_tau = (2-n)*Adv_vel%uhrho_et(i,j,k)*Grd%dxte_dxtr(i,j)*Grd%dyte_dytr(i,j) &
                              +(n-1)*Adv_vel%vhrho_nt(i,j,k)*Grd%dxtn_dxtr(i,j)*Grd%dytn_dytr(i,j)
              Velocity%u(i,j,k,n,taup1) = Grd%tmasken(i,j,k,n)*(rho_dz_vel_tau + dtime*Velocity%accel(i,j,k,n)*Grd%tmasken(i,j,k,n)) &
                                          /(epsln + Thickness%rho_dzten(i,j,k,n))                
           enddo
        enddo
     enddo

     ! sum uprime*rho_dzten
     tmp=0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec  
              tmp(i,j) = tmp(i,j) + Velocity%u(i,j,k,n,taup1)*Thickness%rho_dzten(i,j,k,n)*Grd%tmasken(i,j,k,n)
           enddo
        enddo
     enddo

     ! (vertical mean) - (value from barotropic integration) 
     do j=jsc,jec
        do i=isc,iec  
           tmp(i,j) = Grd%tmasken(i,j,1,n)*(tmp(i,j)-Ext_mode%udrho(i,j,n,taup1))/(Thickness%mass_en(i,j,n)+epsln) 
        enddo
     enddo

     ! replace vertical mean with value from barotropic integration to update velocity 
     do k=1,nk 
        do j=jsc,jec
           do i=isc,iec  
              Velocity%u(i,j,k,n,taup1) = (Velocity%u(i,j,k,n,taup1)-tmp(i,j))*Grd%tmasken(i,j,k,n) 
           enddo
        enddo
     enddo

  enddo  ! finish the n=1,2 do-loop 


  ! override update of velocity in the case of zero tendency 
  if(zero_tendency) then 
     Velocity%u(isc:iec,jsc:jec,:,:,taup1) = Velocity%u(isc:iec,jsc:jec,:,:,tau)
  endif 

  ! override update of u-velocity with value from a file 
  if(use_velocity_override .and. id_velocity_u_override > 0) then 
     wrk1(:,:,:) = 0.0
     call time_interp_external(id_velocity_u_override, Time%model_time, wrk1)
     Velocity%u(isc:iec,jsc:jec,:,1,taup1) = wrk1(isc:iec,jsc:jec,:)*Grd%tmasken(isc:iec,jsc:jec,:,1)
  endif 

  ! override update of v-velocity with value from a file 
  if(use_velocity_override .and. id_velocity_v_override > 0) then 
     wrk1(:,:,:) = 0.0
     call time_interp_external(id_velocity_v_override, Time%model_time, wrk1)
     Velocity%u(isc:iec,jsc:jec,:,2,taup1) = wrk1(isc:iec,jsc:jec,:)*Grd%tmasken(isc:iec,jsc:jec,:,2)
  endif 


  ! for debugging 
  do n=1,2

     if(truncate_velocity) then 
         call velocity_truncate(n, taup1, Velocity)
     endif

     if (debug_this_module .and. n==2) then        
         write(stdoutunit,*) ' ' 
         write(stdoutunit,*) 'From update_ocean_velocity at the u-prime step'
         call write_timestamp(Time%model_time)
         call write_chksum_3d('accel(1)', Velocity%accel(COMP,:,1))
         call write_chksum_3d('accel(2)', Velocity%accel(COMP,:,2))
         call write_chksum_3d('rho_dzten(1)', Thickness%rho_dzten(COMP,:,1))
         call write_chksum_3d('rho_dzten(2)', Thickness%rho_dzten(COMP,:,2))
         call write_chksum_2d('udrho(1)', Ext_mode%udrho(COMP,1,taup1))
         call write_chksum_2d('udrho(2)', Ext_mode%udrho(COMP,2,taup1))
         call write_chksum_2d('mass_en(1)', Thickness%mass_en(COMP,1))
         call write_chksum_2d('mass_en(2)', Thickness%mass_en(COMP,2))
         call ocean_velocity_chksum(Velocity, taum1, write_advection=.false.)
         call ocean_velocity_chksum(Velocity, taup1, write_advection=.false.)
     endif

  enddo


  !--- update velocity at the global halo points to make the gradient 0 across boundary
  if(have_obc) then
     call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','n')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','t')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','t')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','n')
  endif

  if (debug_this_module) then        
     write(stdoutunit,*) ' ' 
     write(stdoutunit,*) 'From update_ocean_velocity: velocity(taup1) chksum'
     call write_timestamp(Time%model_time)
     call ocean_velocity_chksum(Velocity, taup1, write_advection=.false.)
  endif

  ! send diagnostics to diagnostics manager 
  if (id_u(1) > 0) used = send_data (id_u(1), Velocity%u(:,:,:,1,tau), &
                       Time%model_time, rmask=Grd%mask(:,:,:),         &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_u(2) > 0) used = send_data (id_u(2), Velocity%u(:,:,:,2,tau), &
                       Time%model_time, rmask=Grd%mask(:,:,:),         &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_usurf(1) > 0) used = send_data (id_usurf(1), Velocity%u(:,:,1,1,tau), &
                           Time%model_time, rmask=Grd%mask(:,:,1),             &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_usurf(2) > 0) used = send_data (id_usurf(2), Velocity%u(:,:,1,2,tau), &
                           Time%model_time, rmask=Grd%mask(:,:,1),             &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call diagnose_2d(Time, Grd, id_converge_rho_ud_t, Ext_mode%conv_rho_ud_t(:,:,tau))

  if(id_speed > 0) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = sqrt(Velocity%u(i,j,k,1,tau)**2 + Velocity%u(i,j,k,2,tau)**2)
            enddo
         enddo
      enddo
      used = send_data (id_speed, wrk1(:,:,:),       &
             Time%model_time, rmask=Grd%mask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  if(id_ubott(1) > 0 .or. id_ubott(2) > 0) then 
      tmpu = 0.0 ; tmpv = 0.0
      do j=jsd,jed
         do i=isd,ied
            kbot = Grd%kmu(i,j)
            if(kbot > 1) then 
                tmpu(i,j) = Velocity%u(i,j,kbot,1,tau)
                tmpv(i,j) = Velocity%u(i,j,kbot,2,tau)
            endif
         enddo
      enddo
      if (id_ubott(1) > 0) used = send_data (id_ubott(1), tmpu(:,:),      &
                                  Time%model_time, rmask=Grd%mask(:,:,1), &
                                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (id_ubott(2) > 0) used = send_data (id_ubott(2), tmpv(:,:),      &
                                  Time%model_time, rmask=Grd%mask(:,:,1), &
                                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  if(id_u_on_depth(1) > 0) then 
    call remap_s_to_depth(Thickness, Time, Velocity%u(:,:,:,1,tau), 1)
  endif 
  if(id_u_on_depth(2) > 0) then 
    call remap_s_to_depth(Thickness, Time, Velocity%u(:,:,:,2,tau), 2)
  endif 


end subroutine update_ocean_velocity_cgrid
! </SUBROUTINE> NAME="update_ocean_velocity_cgrid"


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_restart">
! <DESCRIPTION>
!
!  Write out restart files registered through register_restart_file
!
! </DESCRIPTION>
subroutine ocean_velocity_restart(Time, Velocity, use_blobs, time_stamp)

   type(ocean_time_type),     intent(in)           :: Time
   type(ocean_velocity_type), intent(in)           :: Velocity
   logical,                   intent(in)           :: use_blobs
   character(len=*),          intent(in), optional :: time_stamp

   integer :: tau, taup1, tau_m1, tau_m0

   tau    = Time%tau
   taup1  = Time%taup1

   tau_m1 = Time%tau_m1
   tau_m0 = Time%tau_m0

   ! resetting the pointer to data
   if(tendency == THREE_LEVEL) then
      call reset_field_pointer(Vel_restart, id_restart_u, Velocity%u(:,:,:,1,tau), Velocity%u(:,:,:,1,taup1) )
      call reset_field_pointer(Vel_restart, id_restart_v, Velocity%u(:,:,:,2,tau), Velocity%u(:,:,:,2,taup1) )
   elseif (tendency==TWO_LEVEL) then
      call reset_field_pointer(Vel_restart, id_restart_u, Velocity%u(:,:,:,1,taup1) )
      call reset_field_pointer(Vel_restart, id_restart_v, Velocity%u(:,:,:,2,taup1) )
      if (use_blobs) then
         call reset_field_pointer(Vel_restart, id_restart_px, Velocity%press_force(:,:,:,1) )
         call reset_field_pointer(Vel_restart, id_restart_py, Velocity%press_force(:,:,:,2) )
      endif
   end if
   call save_restart(Vel_restart, time_stamp)

   if(tendency==TWO_LEVEL) then
      call reset_field_pointer(Adv_restart, id_restart_advu, Velocity%advection(:,:,:,1,tau_m1), &
           Velocity%advection(:,:,:,1,tau_m0) )
      call reset_field_pointer(Adv_restart, id_restart_advv, Velocity%advection(:,:,:,2,tau_m1), &
           Velocity%advection(:,:,:,2,tau_m0) )
      call save_restart(Adv_restart, time_stamp)
   end if

   if(horz_grid == MOM_CGRID) then 
      call reset_field_pointer(Cor_restart, id_restart_coru, Velocity%coriolis(:,:,:,1,tau_m1), &
           Velocity%coriolis(:,:,:,1,tau_m0) )
      call reset_field_pointer(Cor_restart, id_restart_corv, Velocity%coriolis(:,:,:,2,tau_m1), &
           Velocity%coriolis(:,:,:,2,tau_m0) )
      call save_restart(Cor_restart, time_stamp)
   endif 


end subroutine ocean_velocity_restart
! </SUBROUTINE> NAME="ocean_velocity_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_end">
!
! <DESCRIPTION>
!
! Write velocity field to a restart 
!
! </DESCRIPTION>
subroutine ocean_velocity_end(Time, Velocity, use_blobs)

  type(ocean_time_type),     intent(in)  :: Time
  type(ocean_velocity_type), intent(in)  :: Velocity
  logical,                   intent(in)  :: use_blobs
  
  integer :: tau, taup1
  integer :: tau_m0, tau_m1
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau    = Time%tau
  taup1  = Time%taup1

  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0

  call ocean_velocity_restart(Time, Velocity, use_blobs)

  if(tendency==THREE_LEVEL) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_velocity_mod: ending velocity chksum (tau)'
      call write_timestamp(Time%model_time)
      call ocean_velocity_chksum(Velocity, tau, write_advection=.false.)
      write(stdoutunit,*)' ' 
      call write_timestamp(Time%model_time)
      write(stdoutunit,*) 'From ocean_velocity_mod: ending velocity chksum (taup1)'
      call ocean_velocity_chksum(Velocity, taup1, write_advection=.false.)
  elseif(tendency==TWO_LEVEL) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_velocity_mod: ending velocity chksum (taup1)'
      call write_timestamp(Time%model_time)
      call ocean_velocity_chksum(Velocity, taup1, write_advection=.true.)
  endif
  
  return

end subroutine ocean_velocity_end
! </SUBROUTINE> NAME="ocean_velocity_end"


!#######################################################################
! <SUBROUTINE NAME="velocity_truncate">
!
! <DESCRIPTION>
! Truncate velocity so that either component 
! has magnitude no larger than nml specified value. 
!
! </DESCRIPTION>
subroutine velocity_truncate(n, taup1, Velocity)

  integer, intent(in)                      :: n
  integer, intent(in)                      :: taup1
  type(ocean_velocity_type), intent(inout) :: Velocity

  logical :: velocity_truncated 
  integer :: i,j,k

  integer :: stdoutunit 
  stdoutunit=stdout() 

  velocity_truncated = .false. 
  wrk1 = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if (      abs(Grd%yu(i,j)) >= truncate_velocity_lat  &
                .and. abs(Velocity%u(i,j,k,n,taup1)) > truncate_velocity_value) then
               wrk1(i,j,k)               = Velocity%u(i,j,k,n,taup1)
               Velocity%u(i,j,k,n,taup1) = sign(truncate_velocity_value,Velocity%u(i,j,k,n,taup1))
               velocity_truncated        = .true.  
           endif
        enddo
     enddo
  enddo

  if(truncate_verbose .and. velocity_truncated) then 
      write(stdoutunit,'(/a)') ' Summary of baroclinic velocity truncation points'
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if (abs(wrk1(i,j,k)) > 0.0) then 
                   write(stdoutunit,'(a,i1,a,i3,a,i3,a,i3,a,e12.3,a,f12.3,a,f12.3,a,i4)') &
                        'Truncated b/c velocity('                                            &
                        ,n,': ',i+Dom%ioff,',',j+Dom%joff,',',k,') = ',wrk1(i,j,k)           &
                        ,' (m/s) at (lon,lat,klev) = ',Grd%xu(i,j), ',', Grd%yu(i,j), ',', Grd%kmu(i,j)
               endif
            enddo
         enddo
      enddo
  endif

end subroutine velocity_truncate
! </SUBROUTINE> NAME="velocity_truncate"


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_chksum">
!
! <DESCRIPTION>
! Compute checksum for velocity components 
! </DESCRIPTION>
subroutine ocean_velocity_chksum(Velocity, index, write_advection)
 
  type(ocean_velocity_type), intent(in)           :: Velocity
  integer,                   intent(in)           :: index
  logical,                   intent(in), optional :: write_advection

  logical :: writeadvection

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. write_a_restart) then
     write(stdoutunit,'(/a)') '==>Warning from ocean_velocity_mod (ocean_velocity_end): NO restart written.'
     call mpp_error(WARNING,'==>Warning from ocean_velocity_mod (ocean_velocity_end): NO restart written.')
     return
  endif 

  if (PRESENT(write_advection)) then 
    writeadvection = write_advection
  else 
    writeadvection = .true.
  endif 

  call write_chksum_3d('Zonal velocity', Velocity%u(COMP,:,1,index)*Grd%tmasken(COMP,:,1))
  call write_chksum_3d('Meridional velocity', Velocity%u(COMP,:,2,index)*Grd%tmasken(COMP,:,2))

  if(tendency==TWO_LEVEL .and. writeadvection) then 
     call write_chksum_3d('Advection of u', Velocity%advection(COMP,:,1,index)*Grd%tmasken(COMP,:,1))
     call write_chksum_3d('Advection of v', Velocity%advection(COMP,:,2,index)*Grd%tmasken(COMP,:,2))
  endif 

  if(horz_grid == MOM_CGRID) then 
     call write_chksum_3d('x-Coriolis', Velocity%coriolis(COMP,:,1,index)*Grd%tmasken(COMP,:,1))
     call write_chksum_3d('y-Coriolis', Velocity%coriolis(COMP,:,2,index)*Grd%tmasken(COMP,:,2))
  endif 

  return
  
end subroutine ocean_velocity_chksum
! </SUBROUTINE> NAME="ocean_velocity_chksum"



!#######################################################################
! <SUBROUTINE NAME="remap_s_to_depth">
! <DESCRIPTION>
!
! Remap in the vertical from s-coordinate to depth and then send to 
! diagnostic manager.  
!
! This routine is mostly of use for terrain following vertical 
! coordinates, which generally deviate a lot from depth or pressure
! coordinates.  The zstar and pstar coordinates are very similar 
! to z or pressure, so there is no need to do the remapping for 
! purposes of visualization. 
!
! The routine needs to be made more general and faster.  
! It also has been found to be problematic, so it is NOT 
! recommended. It remains here as a template for a better algorithm. 
! Remapping methods in Ferret are much better.  
!
! Use rho_dzu weighting to account for nonBoussinesq.  
!
! Author: Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine remap_s_to_depth(Thickness, Time, array_in, nvelocity)

  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_time_type), intent(in)           :: Time
  real, intent(in),    dimension(isd:,jsd:,:) :: array_in
  integer, intent(in)                         :: nvelocity

  integer :: i, j, k, kk, tau

  tau  = Time%tau
  wrk1 = 0.0
  wrk2 = 0.0
  wrk3 = 0.0

  k=1
  do kk=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Thickness%depth_zu(i,j,kk) < Grd%zw(k)) then 
              wrk1(i,j,k) = wrk1(i,j,k) + array_in(i,j,kk)*Thickness%rho_dzu(i,j,kk,tau)
              wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzu(i,j,kk,tau)
           endif 
        enddo
     enddo
  enddo

  do k=1,nk-1
     do kk=1,nk
        do j=jsd,jed
           do i=isd,ied
              if(Grd%zw(k) <= Thickness%depth_zu(i,j,kk) .and. Thickness%depth_zu(i,j,kk) < Grd%zw(k+1)) then 
                 wrk1(i,j,k+1) = wrk1(i,j,k+1) + array_in(i,j,kk)*Thickness%rho_dzu(i,j,kk,tau)
                 wrk2(i,j,k+1) = wrk2(i,j,k+1) + Thickness%rho_dzu(i,j,kk,tau)
              endif
           enddo
        enddo
     enddo
  enddo

  k=nk
  do kk=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Thickness%depth_zu(i,j,kk) > Grd%zw(k) ) then 
               wrk1(i,j,k) = wrk1(i,j,k) + array_in(i,j,kk)*Thickness%rho_dzu(i,j,kk,tau)
               wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzu(i,j,kk,tau)
           endif
        enddo
     enddo
  enddo

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = Grd%umask_depth(i,j,k)*wrk1(i,j,k)/(wrk2(i,j,k)+epsln)
        enddo
     enddo
  enddo

  call diagnose_3d_u(Time, Grd, id_u_on_depth(nvelocity), wrk3(:,:,:))


end subroutine remap_s_to_depth
! </SUBROUTINE>  NAME="remap_s_to_depth"


end module ocean_velocity_mod
