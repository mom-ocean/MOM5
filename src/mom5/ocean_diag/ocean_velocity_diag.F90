module ocean_velocity_diag_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Numerical diagnostics for velocity related quantities. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains some diagnostics for velocity related quantities.
! Account is taken for either Bgrid or Cgrid layout of the velocity 
! and related discrete fields.  
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_velocity_diag_nml">
!
!  <DATA NAME="diag_step" UNITS="dimensionless" TYPE="integer">
!  Number of time steps between which compute the diagnostics.
!  </DATA> 
!  <DATA NAME="energy_diag_step" TYPE="integer">
!  Perform energy analysis every n timesteps (1==every time step). 
!  This diagnostic is expensive, so should be used sparingly during
!  production runs.  
!  </DATA> 
!
!  <DATA NAME="land_cell_num_max" UNITS="dimensionless" TYPE="integer">
!  Maximum number of land cells where will printout nonzero velocity points.  
!  Default land_cell_num_max=100.
!  </DATA> 
!  <DATA NAME="max_cfl_value" UNITS="dimensionless" TYPE="real">
!  Critical value for Courant number, above which the model will be brought down.
!  </DATA> 
!  <DATA NAME="large_cfl_value" UNITS="dimensionless" TYPE="real">
!  Large value for Courant number, above which will write some diagnostics. 
!  </DATA> 
!  <DATA NAME="verbose_cfl" TYPE="logical">
!  For printing out lots of information about regions of large Courant numbers.
!  </DATA> 
!
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
! For some debugging purposes 
!  </DATA> 
!
!</NAMELIST>

#include <fms_platform.h>

use constants_mod,    only: epsln, c2dbars, omega 
use diag_manager_mod, only: register_diag_field, register_static_field, need_data, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,          only: FATAL, stdout, stdlog
use mpp_domains_mod,  only: mpp_global_sum, mpp_update_domains
use mpp_domains_mod,  only: BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM, BGRID_NE
use mpp_mod,          only: input_nml_file, mpp_error, mpp_max, mpp_sum, mpp_pe 
use mpp_mod,          only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use time_manager_mod, only: time_type, increment_time

use ocean_bih_friction_mod,    only: bih_friction, bih_viscosity_check, bih_reynolds_check
use ocean_coriolis_mod,        only: coriolis_force_bgrid, coriolis_force_cgrid
use ocean_domains_mod,         only: get_local_indices
use ocean_form_drag_mod,       only: form_drag_accel
use ocean_lap_friction_mod,    only: lap_friction, lap_viscosity_check, lap_reynolds_check
use ocean_momentum_source_mod, only: momentum_source
use ocean_operators_mod,       only: FAX, FAY, REMAP_BT_TO_BU, BDX_ET, BDY_NT
use ocean_parameters_mod,      only: DEPTH_BASED, PRESSURE_BASED, missing_value
use ocean_parameters_mod,      only: ENERGETIC, FINITEVOLUME
use ocean_parameters_mod,      only: MOM_BGRID, MOM_CGRID 
use ocean_parameters_mod,      only: rho0, rho0r, grav, onehalf
use ocean_pressure_mod,        only: hydrostatic_pressure, geopotential_anomaly
use ocean_pressure_mod,        only: hydrostatic_pressure_blob, geopotential_anomaly_blob
use ocean_types_mod,           only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_types_mod,           only: ocean_adv_vel_type, ocean_external_mode_type
use ocean_types_mod,           only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,           only: ocean_lagrangian_type
use ocean_util_mod,            only: write_timestamp, matrix, diagnose_2d, diagnose_3d, diagnose_3d_u, diagnose_3d_en
use ocean_util_mod,            only: write_chksum_3d, write_chksum_2d
use ocean_velocity_advect_mod, only: horz_advection_of_velocity, vert_advection_of_velocity
use ocean_vert_mix_mod,        only: vert_friction_bgrid, vert_friction_cgrid
use ocean_workspace_mod,       only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk6 
use ocean_workspace_mod,       only: wrk1_v2d, wrk2_v2d  
use ocean_workspace_mod,       only: wrk1_2d, wrk2_2d, wrk3_2d
use ocean_workspace_mod,       only: wrk1_v, wrk2_v, wrk3_v

implicit none

private

#include <ocean_memory.h>

! useful constants
real :: dbars2Pa
real :: p5grav
real :: p5grav_rho0r
real :: grav_rho0r

real :: dtime 
real :: dtimer 
real :: acor 

! placeholders for Adams-Bashforth Coriolis
real :: abtau_m0=0.0
real :: abtau_m1=0.0
real :: abtau_m2=0.0

! for Bgrid or Cgrid 
integer :: horz_grid 
real, dimension(:,:,:), allocatable :: grd_area
real, dimension(:,:,:), allocatable :: mass_column

! Coriolis parameter on T-point 
real, dimension(:,:), allocatable :: coriolis_t

logical :: used

integer :: id_clock_energy_analysis
integer :: id_clock_press_conversion
integer :: id_clock_press_energy
integer :: id_clock_friction_energy
integer :: id_clock_vert_dissipation

! for potential energy 
integer :: id_pe_tot    =-1
integer :: id_pe_tot_rel=-1
integer :: id_coriolis_t=-1
logical :: potential_energy_first=.true.
real :: pe_tot_0  ! gravitational potential energy (joules) at initial timestep

! for kinetic energy
integer :: id_ke_tot                    =-1
integer :: id_kinetic_energy            =-1
integer :: id_kinetic_energy_intz       =-1
integer :: id_kinetic_energy_clinic     =-1
integer :: id_kinetic_energy_clinic_intz=-1
integer :: id_kinetic_energy_tropic     =-1

! for pressure energy maps
integer :: id_u_dot_grad_pint(2)
integer :: id_u_dot_grad_pint_xy
integer :: id_ubar_dot_grad_pext(2)
integer :: id_ubar_dot_grad_pext_xy

! for friction energy maps
integer :: id_u_dot_lapfrict_horz=-1
integer :: id_v_dot_lapfrict_horz=-1
integer :: id_u_dot_bihfrict_horz=-1
integer :: id_v_dot_bihfrict_horz=-1
integer :: id_u_dot_frict_vert=-1
integer :: id_v_dot_frict_vert=-1
integer :: id_ubar_dot_frict_bt=-1
integer :: id_vbar_dot_frict_bt=-1
integer :: id_vert_lap_diss=-1

! for topostrophy
integer :: id_topostrophy=-1 

! for vorticity
integer :: id_vorticity_z =-1

! for potential vorticity
integer :: id_n2_pv     =-1
integer :: id_pg_pv     =-1
integer :: id_vert_pv   =-1
integer :: id_pg_pvf    =-1
integer :: id_vert_pvf  =-1
integer :: id_bc_pvf    =-1
integer :: id_pvf       =-1
integer :: id_ri_balance=-1
logical :: diagnose_pv = .false.

! pressure conversion diagnostics written in energy analysis subroutine
real :: press_int      ! integrated power from u_int dot grad(p) 
real :: press_ext      ! integrated power from u_ext dot grad(p)
real :: press_conv     ! integrated power converted from u dot grap(p) into other terms
real :: press_conv_err ! error in converting pressure work into other terms   
real :: ucell_mass

! for Stokes-Coriolis force
integer :: id_stokes_force_x =-1
integer :: id_stokes_force_y =-1

! for output
integer :: unit=6

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

logical :: module_is_initialized = .FALSE.

character(len=128) :: version=&
     '$Id: ocean_velocity_diag.F90,v 20.0 2013/12/14 00:12:59 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

public  ocean_velocity_diag_init
public  ocean_velocity_diagnostics 

public  velocity_change
public  kinetic_energy
public  potential_energy

public  pressure_conversion
public  energy_analysis

private pressure_energy
private friction_energy
private vert_dissipation 
private velocity_land_cell_check
private cfl_check1_bgrid
private cfl_check1_cgrid
private cfl_check2_bgrid
private cfl_check2_cgrid
private compute_topostrophy
private compute_vorticity
private stokes_coriolis_force
private diagnose_kinetic_energy_maps
private diagnose_potential_vorticity 

integer :: global_sum_flag
integer :: diag_step         = -1
integer :: energy_diag_step  = -1
integer :: land_cell_num_max = 100  

! for CFL checks 
real    :: max_cfl_value   = 100.0           
real    :: large_cfl_value = 10.0 
logical :: verbose_cfl     =.false.

! internally set for diagnosing kinetic energy maps
logical :: compute_kinetic_energy_maps = .false.

logical :: do_bitwise_exact_sum = .false.
logical :: debug_this_module    = .false.


namelist /ocean_velocity_diag_nml/ diag_step, energy_diag_step, do_bitwise_exact_sum, debug_this_module, &
                                   max_cfl_value, large_cfl_value, verbose_cfl, land_cell_num_max


contains

!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_velocity_diag module containing subroutines
! diagnosing velocity related properties of the simulation.  These are 
! not terms in the equations, but rather they are diagnosed from 
! terms. 
! </DESCRIPTION>
!
subroutine ocean_velocity_diag_init(Grid, Domain, Time, Time_steps, hor_grid)

type(ocean_grid_type),       target, intent(in) :: Grid
type(ocean_domain_type),     target, intent(in) :: Domain
type(ocean_time_type),               intent(in) :: Time
type(ocean_time_steps_type),         intent(in) :: Time_steps
integer,                             intent(in) :: hor_grid 

integer :: ioun, io_status, ierr
integer :: i,j,k

integer :: stdoutunit,stdlogunit 
stdoutunit=stdout();stdlogunit=stdlog() 

if (module_is_initialized) return

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_velocity_diag_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_velocity_diag_nml')
#else
ioun = open_namelist_file()
read(ioun, ocean_velocity_diag_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_velocity_diag_nml')
call close_file(ioun)
#endif
write (stdoutunit,'(/)')
write (stdoutunit, ocean_velocity_diag_nml)
write (stdlogunit, ocean_velocity_diag_nml)

if (diag_step == 0)        diag_step = 1
if (energy_diag_step == 0) energy_diag_step = 1

if(do_bitwise_exact_sum) then
   global_sum_flag = BITWISE_EXACT_SUM
else
   global_sum_flag = NON_BITWISE_EXACT_SUM
endif

 if(debug_this_module) then 
   write(stdoutunit,*) &
   '==>Note: running ocean_velocity_diag_mod with debug_this_module=.true.'
 endif 

#ifndef MOM_STATIC_ARRAYS
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

Dom => Domain
Grd => Grid

p5grav       = 0.5*grav
grav_rho0r   = grav*rho0r
p5grav_rho0r = 0.5*grav*rho0r
dbars2Pa     = 1.0/c2dbars 

dtime        = Time_steps%dtime_u
dtimer       = 1.0/dtime
acor         = Time_steps%acor
horz_grid    = hor_grid 

! pressure energetics 
id_u_dot_grad_pint(1) = register_diag_field ('ocean_model', 'u_dot_grad_pint_x', &
     Grd%vel_axes_u(1:3), Time%model_time,                                       &
     'i-current times 3d i-pressure force', 'Watt',                              &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_u_dot_grad_pint(2) = register_diag_field ('ocean_model', 'u_dot_grad_pint_y', &
     Grd%vel_axes_v(1:3), Time%model_time,                                       &
     'j-current times 3d j-pressure force', 'Watt',                              &
     missing_value=missing_value, range=(/-1e15, 1e15/))

id_ubar_dot_grad_pext(1) = register_diag_field ('ocean_model', 'ubar_dot_grad_pext_x', &
     Grd%vel_axes_u(1:2), Time%model_time,                                             &
     '2d i-current times 2d i-pressure force', 'Watt',                                 &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_ubar_dot_grad_pext(2) = register_diag_field ('ocean_model', 'ubar_dot_grad_pext_y', &
     Grd%vel_axes_v(1:2), Time%model_time,                                             &
     '2d j-current times 2d j-pressure force', 'Watt',                                 &
     missing_value=missing_value, range=(/-1e15, 1e15/))

! friction energetics
id_u_dot_lapfrict_horz = register_diag_field ('ocean_model', 'u_dot_lapfrict_horz', &
     Grd%vel_axes_u(1:3), Time%model_time,                                          &
     'i-current times horizontal laplacian friction', 'Watt',                       &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_v_dot_lapfrict_horz = register_diag_field ('ocean_model', 'v_dot_lapfrict_horz', &
     Grd%vel_axes_v(1:3), Time%model_time,                                          &
     'j-current times horizontal laplacian friction', 'Watt',                       &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_u_dot_bihfrict_horz = register_diag_field ('ocean_model', 'u_dot_bihfrict_horz', &
     Grd%vel_axes_u(1:3), Time%model_time,                                          &
     'i-current times horizontal biharmonic friction', 'Watt',                      &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_v_dot_bihfrict_horz = register_diag_field ('ocean_model', 'v_dot_bihfrict_horz', &
     Grd%vel_axes_v(1:3), Time%model_time,                                          &
     'j-current times horizontal biharmonic friction', 'Watt',                      &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_u_dot_frict_vert = register_diag_field ('ocean_model', 'u_dot_frict_vert', &
     Grd%vel_axes_u(1:3), Time%model_time,                                    &
     'i-current times vertical friction', 'Watt',                             &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_v_dot_frict_vert = register_diag_field ('ocean_model', 'v_dot_frict_vert', &
     Grd%vel_axes_v(1:3), Time%model_time,                                    &
     'j-current times vertical friction', 'Watt',                             &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_ubar_dot_frict_bt = register_diag_field ('ocean_model', 'ubar_dot_frict_bt', &
     Grd%vel_axes_u(1:2), Time%model_time,                                      &
     'i-ubar current times horizontal friction applied just to ubar', 'Watt',   &
     missing_value=missing_value, range=(/-1e15, 1e15/))
id_vbar_dot_frict_bt = register_diag_field ('ocean_model', 'vbar_dot_frict_bt', &
     Grd%vel_axes_v(1:2), Time%model_time,                                      &
     'j-vbar current times horizontal friction applied just to vbar', 'Watt',   &
     missing_value=missing_value, range=(/-1e15, 1e15/))

id_topostrophy = register_diag_field ('ocean_model', 'topostrophy', Grd%vel_axes_uv(1:3),&
       Time%model_time, 'topostrophy = f * (zhat cross u) dot grad H', 'm/s^2',          &
       missing_value=missing_value, range=(/-1.e20,1.e20/))

id_stokes_force_x = register_diag_field ('ocean_model', 'stokes_force_x', Grd%vel_axes_u(1:3), &
     Time%model_time, 'Zonal component to Stokes Coriolis force', '(kg/m^3)*(m^2/s^2)',        &
      missing_value=missing_value, range=(/-1e9,1e9/))

id_stokes_force_y = register_diag_field ('ocean_model', 'stokes_force_y', Grd%vel_axes_v(1:3), &
      Time%model_time, 'Meridional component to Stokes Coriolis force', '(kg/m^3)*(m^2/s^2)',  &
      missing_value=missing_value, range=(/-1e9,1e9/))

! kinetic energy scalar
id_ke_tot = register_diag_field ('ocean_model', 'ke_tot', Time%model_time, &
            'Globally integrated ocean kinetic energy', '10^15 Joules',    &
            missing_value=missing_value, range=(/0.0,1e20/))

! potential energy
id_pe_tot = register_diag_field ('ocean_model', 'pe_tot', Time%model_time, &
            'Globally integrated ocean potential energy', '10^15 Joules',  &
            missing_value=missing_value, range=(/0.0,1e20/))
id_pe_tot_rel = register_diag_field ('ocean_model', 'pe_tot_rel', Time%model_time,          &
                'global potential energy - initial global potential energy', '10^15 Joules',&
                missing_value=missing_value, range=(/0.0,1e20/))

! Coriolis on T-point 
id_coriolis_t = register_static_field ('ocean_model', 'coriolis_t', Grd%tracer_axes(1:2), &
                                       'Coriolis parameter on T-cell', '1/s',             &
                                       missing_value=missing_value, range=(/-10.0,10.0/))
allocate (coriolis_t(isd:ied,jsd:jed))
do j=jsc,jec
   do i=isc,iec
      coriolis_t(i,j) = 2.0*omega*sin(Grd%phit(i,j))*Grd%tmask(i,j,1) 
   enddo
enddo
call diagnose_2d(Time, Grd, id_coriolis_t, coriolis_t(:,:))

! Ertel PV and its pieces
id_n2_pv = register_diag_field ('ocean_model', 'n2_pv', Grd%tracer_axes(1:3), Time%model_time,&
  'squared BV frequency for planetary geostrophic PV: N^2', '1/sec^2',                        &
  missing_value=missing_value, range=(/-1e6,1e6/))
id_pg_pv = register_diag_field ('ocean_model', 'pg_pv', Grd%tracer_axes(1:3), Time%model_time,&
  'planetary geostrophic PV: f*N^2', '1/sec^3', missing_value=missing_value, range=(/-1e6,1e6/))
id_vert_pv = register_diag_field ('ocean_model', 'vert_pv', Grd%tracer_axes(1:3), Time%model_time,&
  'vertical piece of Ertel PV: (f+zeta)*N^2', '1/sec^3', missing_value=missing_value, range=(/-1e6,1e6/))
id_pg_pvf = register_diag_field ('ocean_model', 'pg_pvf', Grd%tracer_axes(1:3), Time%model_time,&
  'planetary geostrophic PV times f: (f*N)^2', '1/sec^4', missing_value=missing_value, range=(/-1e6,1e6/))
id_vert_pvf = register_diag_field ('ocean_model', 'vert_pvf', Grd%tracer_axes(1:3), Time%model_time,&
  'vertical piece of Ertel PV times f: f*(f+zeta)*N^2', '1/sec^4',                                  &
  missing_value=missing_value, range=(/-1e6,1e6/))
id_bc_pvf = register_diag_field ('ocean_model', 'bc_pvf', Grd%tracer_axes(1:3), Time%model_time,&
  'baroclinic component of balanced Ertel PV times f: -|\nabla buoy|^2', '1/sec^4',             &
  missing_value=missing_value, range=(/-1e6,1e6/))
id_pvf = register_diag_field ('ocean_model', 'pvf', Grd%tracer_axes(1:3), Time%model_time,&
  'Ertel PV times f', '1/sec^4', missing_value=missing_value, range=(/-1e6,1e6/))
id_ri_balance = register_diag_field ('ocean_model', 'ri_balance', Grd%tracer_axes(1:3), Time%model_time,&
  'Balanced Richardson number', 'dimensionless', missing_value=missing_value, range=(/-1e6,1e6/))
if (id_pg_pv > 0)      diagnose_pv = .true.
if (id_vert_pv > 0)    diagnose_pv = .true.
if (id_pg_pvf > 0)     diagnose_pv = .true.
if (id_vert_pvf > 0)   diagnose_pv = .true.
if (id_bc_pvf > 0)     diagnose_pv = .true.
if (id_pvf > 0)        diagnose_pv = .true.
if (id_ri_balance > 0) diagnose_pv = .true.


if(horz_grid == MOM_BGRID) then 
   id_u_dot_grad_pint_xy = register_diag_field ('ocean_model', 'u_dot_grad_pint_xy', &
        Grd%vel_axes_uv(1:3), Time%model_time,                                       &
        'current dot 3d pressure force', 'Watt',                                     &
        missing_value=missing_value, range=(/-1e15, 1e15/))
   id_ubar_dot_grad_pext_xy = register_diag_field ('ocean_model', 'ubar_dot_grad_pext_xy',&
        Grd%vel_axes_uv(1:2), Time%model_time,                                            &
        '2d current dot 2d pressure force', 'Watt',                                       &
        missing_value=missing_value, range=(/-1e15, 1e15/))
   id_vert_lap_diss = register_diag_field ('ocean_model', 'vert_lap_diss', Grd%vel_axes_uv(1:3),&
          Time%model_time, 'Energy dissipation from vertical Laplacian friction', 'W/m^2',      &
          missing_value=missing_value, range=(/-1.e20,1.e20/),                                  &
          standard_name='ocean_kinetic_energy_dissipation_per_unit_area_due_to_vertical_friction')
   id_vorticity_z = register_diag_field ('ocean_model', 'vorticity_z', Grd%tracer_axes(1:3), Time%model_time,&
     'vertical vorticity component: v_x - u_y', '1/sec', missing_value=missing_value, range=(/-1e6,1e6/))
   id_kinetic_energy = register_diag_field ('ocean_model', 'kinetic_energy', Grd%vel_axes_uv(1:3),&
          Time%model_time, 'kinetic energy in a grid cell', '10^15 Joules',                       &
          missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_intz = register_diag_field ('ocean_model', 'kinetic_energy_intz',  &
          Grd%vel_axes_uv(1:2), Time%model_time, 'Vertically integrated kinetic energy',&
          '10^15 Joules', missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_intz > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_clinic = register_diag_field ('ocean_model', 'kinetic_energy_clinic',   &
          Grd%vel_axes_uv(1:3), Time%model_time, 'Baroclinic kinetic energy', '10^15 Joules',&
          missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_clinic > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_clinic_intz = register_diag_field ('ocean_model', 'kinetic_energy_clinic_intz',&
          Grd%vel_axes_uv(1:2), Time%model_time, 'Vertically integrated baroclinic kinetic energy', &
          '10^15 Joules', missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_clinic_intz > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_tropic = register_diag_field ('ocean_model', 'kinetic_energy_tropic',   &
          Grd%vel_axes_uv(1:2), Time%model_time, 'Barotropic kinetic energy', '10^15 Joules',&
          missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_tropic > 0) compute_kinetic_energy_maps = .true.

else  ! Cgrid
   id_u_dot_grad_pint_xy = register_diag_field ('ocean_model', 'u_dot_grad_pint_xy',&
        Grd%tracer_axes(1:3), Time%model_time,                                      &
        'current dot 3d pressure force', 'Watt',                                    &
        missing_value=missing_value, range=(/-1e15, 1e15/))
   id_ubar_dot_grad_pext_xy = register_diag_field ('ocean_model', 'ubar_dot_grad_pext_xy',&
        Grd%tracer_axes(1:2), Time%model_time,                                            &
        '2d current dot 2d pressure force', 'Watt',                                       &
        missing_value=missing_value, range=(/-1e15, 1e15/))
   id_vert_lap_diss = register_diag_field ('ocean_model', 'vert_lap_diss', Grd%tracer_axes(1:3),&
          Time%model_time, 'Energy dissipation from vertical Laplacian friction', 'W/m^2',      &
          missing_value=missing_value, range=(/-1.e20,1.e20/),                                  &
          standard_name='ocean_kinetic_energy_dissipation_per_unit_area_due_to_vertical_friction')
   id_vorticity_z = register_diag_field ('ocean_model', 'vorticity_z', Grd%vel_axes_uv(1:3), Time%model_time, &
     'vertical vorticity component: v_x - u_y', '1/sec', missing_value=missing_value, range=(/-1e6,1e6/))


   id_kinetic_energy = register_diag_field ('ocean_model', 'kinetic_energy', Grd%tracer_axes(1:3),&
          Time%model_time, 'kinetic energy in a grid cell', '10^15 Joules',                       &
          missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_intz = register_diag_field ('ocean_model', 'kinetic_energy_intz',  &
          Grd%tracer_axes(1:2), Time%model_time, 'Vertically integrated kinetic energy',&
          '10^15 Joules', missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_intz > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_clinic = register_diag_field ('ocean_model', 'kinetic_energy_clinic',   &
          Grd%tracer_axes(1:3), Time%model_time, 'Baroclinic kinetic energy', '10^15 Joules',&
          missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_clinic > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_clinic_intz = register_diag_field ('ocean_model', 'kinetic_energy_clinic_intz',&
          Grd%tracer_axes(1:2), Time%model_time, 'Vertically integrated baroclinic kinetic energy', &
          '10^15 Joules', missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_clinic_intz > 0) compute_kinetic_energy_maps = .true.

   id_kinetic_energy_tropic = register_diag_field ('ocean_model', 'kinetic_energy_tropic',   &
          Grd%tracer_axes(1:2), Time%model_time, 'Barotropic kinetic energy', '10^15 Joules',&
          missing_value=missing_value, range=(/0.0,1.e20/))
   if(id_kinetic_energy_tropic > 0) compute_kinetic_energy_maps = .true.

endif 


allocate (grd_area(isd:ied,jsd:jed,nk))
allocate (mass_column(isd:ied,jsd:jed,2))
mass_column = 0.0
if(horz_grid == MOM_BGRID) then 
   do k=1,nk
      do j=jsd,jed
         do i=isd,ied
            grd_area(i,j,k) = Grd%umask(i,j,k)*Grd%dau(i,j)
         enddo
      enddo
   enddo
else 
   do k=1,nk
      do j=jsd,jed
         do i=isd,ied
            grd_area(i,j,k) = Grd%tmask(i,j,k)*Grd%dat(i,j)
         enddo
      enddo
   enddo
endif 

id_clock_energy_analysis   = mpp_clock_id('(Velocity diag: energy analysis)'  ,grain=CLOCK_MODULE)
id_clock_press_conversion  = mpp_clock_id('(Velocity diag: press convert)'    ,grain=CLOCK_MODULE)
id_clock_press_energy      = mpp_clock_id('(Velocity diag: press energy)'     ,grain=CLOCK_MODULE)
id_clock_friction_energy   = mpp_clock_id('(Velocity diag: friction energy)'  ,grain=CLOCK_MODULE)
id_clock_vert_dissipation  = mpp_clock_id('(Velocity diag: vert dissipation'  ,grain=CLOCK_MODULE)


end subroutine ocean_velocity_diag_init
! </SUBROUTINE>  NAME="ocean_velocity_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_diagnostics">
!
! <DESCRIPTION>
! Call diagnostics related to the velocity. 
! </DESCRIPTION>
!
subroutine ocean_velocity_diagnostics(Time, Thickness, Dens, Ext_mode, Velocity)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity

  type(time_type) :: next_time 

  real :: pe_tot     = 0.0
  real :: pe_tot_rel = 0.0
  real :: ke_tot     = 0.0

  ! diagnostics that do not send_data to diag_manager  
  if (diag_step > 0) then 
     if(mod(Time%itt,diag_step) == 0) then
        call lap_viscosity_check
        call lap_reynolds_check(Time, Velocity)
        call bih_viscosity_check
        call bih_reynolds_check(Time, Velocity)
        call velocity_land_cell_check(Time, Velocity)
        call write_timestamp(Time%model_time)
        call kinetic_energy(Time, Thickness, Velocity, ke_tot, .false., .true.)
        call potential_energy(Time, Thickness, Dens, pe_tot, pe_tot_rel, .false., .true.)
        call cfl_check1_bgrid(Time, Thickness, Velocity)
        call cfl_check1_cgrid(Time, Thickness, Velocity)
        call cfl_check2_bgrid(Time, Thickness, Velocity)
        call cfl_check2_cgrid(Time, Thickness, Velocity)
     endif  
  endif

  ! diagnostics that send_data to diag_manager 
  next_time = increment_time(Time%model_time, int(dtime), 0)
  if (need_data(id_ke_tot,next_time)) then 
    call kinetic_energy(Time, Thickness, Velocity, ke_tot, .true., .false.)
  endif 
  if (need_data(id_pe_tot,next_time)) then 
    call potential_energy(Time, Thickness, Dens, pe_tot, pe_tot_rel, .true., .false.)
  endif 

  ! diagnose topography scalar 
  call compute_topostrophy(Time, Velocity)

  ! diagnose vorticity 
  call compute_vorticity(Time, Velocity)

  ! diagnostic Stokes-Coriolis force 
  call stokes_coriolis_force(Time, Thickness, Velocity)

  ! diagnostic kinetic energy maps
  call diagnose_kinetic_energy_maps(Time, Thickness, Velocity, Ext_mode)

  ! diagnostic potential vorticity and terms 
  call diagnose_potential_vorticity(Time, Velocity, Dens)
  
end subroutine ocean_velocity_diagnostics
! </SUBROUTINE>  NAME="ocean_velocity_diagnostics"


!#######################################################################
! <SUBROUTINE NAME="potential_energy">
!
! <DESCRIPTION>
!
! Compute gravitational potential energy (Joules) relative to z=0 
! taken with respect to the value at the initial time step. 
!
! </DESCRIPTION>
!
subroutine potential_energy (Time, Thickness, Dens, pe_tot, pe_tot_rel, diag_flag, write_flag)

  type(ocean_time_type),      intent(in)           :: Time
  type(ocean_thickness_type), intent(in)           :: Thickness
  type(ocean_density_type),   intent(in)           :: Dens
  real,                       intent(inout)        :: pe_tot
  real,                       intent(inout)        :: pe_tot_rel
  logical,                    intent(in), optional :: diag_flag 
  logical,                    intent(in), optional :: write_flag 

  integer :: i, j, k
  integer :: taup1
  integer :: stdoutunit 
  logical :: send_diagnostics
  logical :: write_diagnostics

  stdoutunit=stdout() 
  taup1 = Time%taup1

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  ! assume write_diagnostics=.true., unless write_flag says it is false.
  if (PRESENT(write_flag)) then 
    write_diagnostics = write_flag
  else 
    write_diagnostics = .true.
  endif 

  ! Calculate total potential energy
  pe_tot = 0.
  do k= 1, nk
    do j = jsc, jec
      do i = isc, iec
        pe_tot = pe_tot + Grd%dat(i, j) * Thickness%dzt(i, j, k) &
                     * Dens%rho(i, j, k, taup1) * Grd%tmask(i, j, k) &
                     * Thickness%depth_zt(i, j, k) * grav
      end do
    end do
  end do
  call mpp_sum(pe_tot)

  pe_tot_rel = 0.
  if (potential_energy_first) then
    potential_energy_first = .false.
    pe_tot_0 = pe_tot
  else
    ! computing pe_tot relative to initial time step value reduces roundoff errors
    pe_tot_rel = pe_tot_rel + (pe_tot - pe_tot_0)
  end if

  if(write_diagnostics) then 
     write(stdoutunit,'(/,1x,a,e20.12)') &
     'Gravitational potential energy relative to first time (Joules) = ', pe_tot_rel
     write(stdoutunit,'(/,1x,a,e20.12)') &
     'Gravitational potential energy                                 = ', pe_tot
  endif 

  if(send_diagnostics) then
    if (id_pe_tot_rel > 0) used = send_data (id_pe_tot_rel, pe_tot_rel*1e-15, Time%model_time)
  endif 
  if(send_diagnostics) then
    if (id_pe_tot > 0) used = send_data (id_pe_tot, pe_tot*1e-15, Time%model_time)
  endif 

end subroutine potential_energy
! </SUBROUTINE>  NAME="potential_energy"


!#######################################################################
! <SUBROUTINE NAME="kinetic_energy">
!
! <DESCRIPTION>
! Compute global integrated horizontal kinetic energy.
! </DESCRIPTION>
!
subroutine kinetic_energy (Time, Thickness, Velocity, ke_tot, diag_flag, write_flag)

  type(ocean_time_type),      intent(in)           :: Time
  type(ocean_thickness_type), intent(in)           :: Thickness
  type(ocean_velocity_type),  intent(in)           :: Velocity
  real,                       intent(inout)        :: ke_tot
  logical,                    intent(in), optional :: diag_flag 
  logical,                    intent(in), optional :: write_flag 

  logical :: send_diagnostics
  logical :: write_diagnostics
  real    :: mass_cell
  integer :: i, j, k, tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  ! assume write_diagnostics=.true., unless write_flag says it is false.
  if (PRESENT(write_flag)) then 
    write_diagnostics = write_flag
  else 
    write_diagnostics = .true.
  endif 

  tau = Time%tau

  ! total horizontal kinetic energy(joules) at "tau"
  ke_tot = 0.0
  if(horz_grid == MOM_BGRID) then 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              mass_cell = grd_area(i,j,k)*Thickness%rho_dzu(i,j,k,tau)
              ke_tot    = ke_tot + 0.5*mass_cell*(Velocity%u(i,j,k,1,tau)**2 + Velocity%u(i,j,k,2,tau)**2)
           enddo
        enddo
     enddo

  else  ! cgrid 

     ! ignore offset in velocity components 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              mass_cell = Grd%dat(i,j)*Thickness%rho_dzt(i,j,k,tau)*Grd%tmask(i,j,k) 
              ke_tot    = ke_tot + 0.5*mass_cell*(Velocity%u(i,j,k,1,tau)**2 + Velocity%u(i,j,k,2,tau)**2)
           enddo
        enddo
     enddo

  endif 

  call mpp_sum (ke_tot)

  if(write_diagnostics) then 
     write(stdoutunit,'(1x,a,e20.12)') 'Kinetic energy (Joules)                          = ', ke_tot
  endif 

  if(send_diagnostics) then 
    if (id_ke_tot > 0) used = send_data (id_ke_tot, ke_tot*1e-15, Time%model_time)
  endif 

end subroutine kinetic_energy
! </SUBROUTINE> NAME="kinetic_energy"


!#######################################################################
! <SUBROUTINE NAME="diagnose_kinetic_energy_maps">
!
! <DESCRIPTION>
! Compute maps of horizontal kinetic energy, which is the kinetic
! energy appropriate for a hydrostatic fluid.
!
! Note that for Cgrid calculation we ignore the
! offset in the horizontal velocity components.
!
! </DESCRIPTION>
!
subroutine diagnose_kinetic_energy_maps (Time, Thickness, Velocity, Ext_mode)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_velocity_type),      intent(in) :: Velocity
  type(ocean_external_mode_type), intent(in) :: Ext_mode

  integer :: i, j, k, tau
  real    :: uprime, vprime
  real, dimension(isd:ied,jsd:jed,2)    :: ubar
  real, dimension(isd:ied,jsd:jed,2)    :: mass_per_area
  real, dimension(isd:ied,jsd:jed,nk,2) :: mass_per_cell

  if(.not. compute_kinetic_energy_maps) return

  tau                    = Time%tau
  wrk1(:,:,:)            = 0.0
  wrk2(:,:,:)            = 0.0
  wrk1_2d(:,:)           = 0.0
  wrk2_2d(:,:)           = 0.0
  wrk3_2d(:,:)           = 0.0
  mass_per_area(:,:,:)   = 0.0
  mass_per_cell(:,:,:,:) = 0.0

  ! compute masses appropriate for Bgrid or Cgrid
  if(horz_grid == MOM_BGRID) then

     do j=jsd,jed
        do i=isd,ied
           mass_per_area(i,j,1) = Thickness%mass_u(i,j,tau)
           mass_per_area(i,j,2) = Thickness%mass_u(i,j,tau)
        enddo
     enddo
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              mass_per_cell(i,j,k,1) = grd_area(i,j,k)*Thickness%rho_dzu(i,j,k,tau)
              mass_per_cell(i,j,k,2) = mass_per_cell(i,j,k,1)
           enddo
        enddo
     enddo

  else ! Cgrid

     do j=jsd,jed
        do i=isd,ied
           mass_per_area(i,j,1) = Thickness%mass_en(i,j,1)
           mass_per_area(i,j,2) = Thickness%mass_en(i,j,2)
        enddo
     enddo
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              mass_per_cell(i,j,k,1) = grd_area(i,j,k)*Thickness%rho_dzten(i,j,k,1)
              mass_per_cell(i,j,k,2) = grd_area(i,j,k)*Thickness%rho_dzten(i,j,k,2)
           enddo
        enddo
     enddo

  endif


  ! vertical average of horizontal velocity and kinetic energy in this vertical average
  do j=jsd,jed
     do i=isd,ied
        ubar(i,j,1)  = Grd%tmasken(i,j,1,1)*Ext_mode%udrho(i,j,1,tau)/(mass_per_area(i,j,1)+epsln)
        ubar(i,j,2)  = Grd%tmasken(i,j,1,2)*Ext_mode%udrho(i,j,2,tau)/(mass_per_area(i,j,2)+epsln)
        wrk3_2d(i,j) = onehalf*grd_area(i,j,1)*(mass_per_area(i,j,1)*ubar(i,j,1)**2 + mass_per_area(i,j,2)*ubar(i,j,2)**2)
     enddo
  enddo

  ! kinetic energy in full flow and in baroclinic flow
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           uprime      = Velocity%u(i,j,k,1,tau) - ubar(i,j,1)
           vprime      = Velocity%u(i,j,k,2,tau) - ubar(i,j,2)
           wrk1(i,j,k) = onehalf*(mass_per_cell(i,j,k,1)*Velocity%u(i,j,k,1,tau)**2 + mass_per_cell(i,j,k,2)*Velocity%u(i,j,k,2,tau)**2)
           wrk2(i,j,k) = onehalf*(mass_per_cell(i,j,k,1)*uprime**2 + mass_per_cell(i,j,k,2)*vprime**2)
        enddo
     enddo
  enddo

  ! vertical integrals
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1_2d(i,j) = wrk1_2d(i,j) + wrk1(i,j,k)
           wrk2_2d(i,j) = wrk2_2d(i,j) + wrk2(i,j,k)
        enddo
     enddo
  enddo

  if(id_kinetic_energy > 0) then
     used = send_data (id_kinetic_energy, wrk1(:,:,:)*1e-15, Time%model_time, rmask=Grd%mask(:,:,:), &
            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_kinetic_energy_intz > 0) then
     used = send_data (id_kinetic_energy_intz, wrk1_2d(:,:)*1e-15, Time%model_time, rmask=Grd%mask(:,:,1), &
            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if(id_kinetic_energy_clinic > 0) then
     used = send_data (id_kinetic_energy_clinic, wrk2(:,:,:)*1e-15, Time%model_time, rmask=Grd%mask(:,:,:), &
            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_kinetic_energy_clinic_intz > 0) then
     used = send_data (id_kinetic_energy_clinic_intz, wrk2_2d(:,:)*1e-15, Time%model_time, rmask=Grd%mask(:,:,1), &
            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if(id_kinetic_energy_tropic > 0) then
     used = send_data (id_kinetic_energy_tropic, wrk3_2d(:,:)*1e-15, Time%model_time, rmask=Grd%mask(:,:,1), &
            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


end subroutine diagnose_kinetic_energy_maps
! </SUBROUTINE> NAME="diagnose_kinetic_energy_maps"


!#######################################################################
! <SUBROUTINE NAME="velocity_land_cell_check">
!
! <DESCRIPTION>
! See if there are any points over land with nonzero ocean velocity 
! </DESCRIPTION>
!
subroutine velocity_land_cell_check(Time, Velocity)

  implicit none
  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  integer :: i, j, k, num, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taup1 = Time%taup1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_velocity_diag_mod(velocity_land_cell_check): module needs initialization ')
  endif 

  call write_timestamp(Time%model_time)
  write (stdoutunit,'(//60x,a/)') ' Land cell summary:'
  write (stdoutunit,'(1x,a/)')'Locations (if any) where land cell velocity is non-zero...'
  num = 0

  if(horz_grid == MOM_BGRID) then 
     kloop: do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if (Grd%umask(i,j,k) == 0 .and. (Velocity%u(i,j,k,1,taup1) /= 0.0 .or. Velocity%u(i,j,k,2,taup1) /= 0.0)) then
                 num = num + 1
                 write (unit,9000) i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), Velocity%u(i,j,k,1,taup1), Velocity%u(i,j,k,2,taup1)
              endif
              if(num > land_cell_num_max) then
                 write (stdoutunit,'(1x,a/)')'More than land_cell_num_max land cells with non-zero velocity. Stop the check.'
                 exit kloop          
              endif 
           enddo
        enddo
     enddo kloop

  else  ! cgrid 

     kloop2: do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if (Grd%tmask(i,j,k) == 0 .and. (Velocity%u(i,j,k,1,taup1) /= 0.0 .or. Velocity%u(i,j,k,2,taup1) /= 0.0)) then
                 num = num + 1
                 write (unit,9000) i, j, k, Grd%xt(i,j), Grd%yt(i,j), Grd%zt(k), Velocity%u(i,j,k,1,taup1), Velocity%u(i,j,k,2,taup1)
              endif
              if(num > land_cell_num_max) then
                 write (stdoutunit,'(1x,a/)')'More than land_cell_num_max land cells with non-zero velocity. Stop the check.'
                 exit kloop2          
              endif 
           enddo
        enddo
     enddo kloop2

  endif 

  if (num > 0)  then 
     call mpp_error(FATAL,'==>Error: found nonzero ocean velocity over land points. ')
  endif 

9000  format(/' " =>Error: Land cell at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m) has u=',e10.3, ' v=',e10.3)

end subroutine velocity_land_cell_check
! </SUBROUTINE> NAME="velocity_land_cell_check"



!#######################################################################
! <SUBROUTINE NAME="velocity_change">
!
! <DESCRIPTION>
! Determine the number of points that have large single-time step 
! changes in the absolute value of the velocity.  
! </DESCRIPTION>
!
subroutine velocity_change(Time, Velocity, velocity_change_max, velocity_change_max_num)

  implicit none
  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity
  real,                      intent(in) :: velocity_change_max
  integer,                   intent(in) :: velocity_change_max_num

  real     :: abschange(2)
  integer  :: i, j, k, n, num
  integer  :: taum1, tau, taup1 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  call write_timestamp(Time%model_time)
  write (stdoutunit,'(//60x,a/)')    &
   ' Velocity change summary (test for leap-frog noise):'
  write (stdoutunit,'(1x,a,e12.6/)') &
   'Locations (if any) where abs(0.5*(u(taup1)+u(taum1))-u(tau)) (m/s) > ',velocity_change_max
  num = 0
  if(horz_grid == MOM_BGRID) then 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if(Grd%umask(i,j,k) > 0) then 
                 do n=1,2
                    abschange(n) = abs(0.5*(Velocity%u(i,j,k,n,taup1)+Velocity%u(i,j,k,n,taum1))-Velocity%u(i,j,k,n,tau))
                 enddo
                 if (abschange(1) > velocity_change_max .or. abschange(2) > velocity_change_max) then
                    num = num + 1
                    write (unit,9000) i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), abschange(1), abschange(2)
                 endif 
              endif
           enddo
        enddo
     enddo

  else 

     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if(Grd%tmask(i,j,k) > 0) then 
                 do n=1,2
                    abschange(n) = abs(0.5*(Velocity%u(i,j,k,n,taup1)+Velocity%u(i,j,k,n,taum1))-Velocity%u(i,j,k,n,tau))
                 enddo
                 if (abschange(1) > velocity_change_max .or. abschange(2) > velocity_change_max) then
                    num = num + 1
                    write (unit,9000) i, j, k, Grd%xt(i,j), Grd%yt(i,j), Grd%zt(k), abschange(1), abschange(2)
                 endif 
              endif
           enddo
        enddo
     enddo

  endif 

  call mpp_max(num)
  if (num > velocity_change_max_num) then 
    call mpp_error(FATAL, &
    'Found too many points where velocity changed too much over single time step.')
  endif 

9000  format(/' " =>Error: Ocean at (i,j,k) = ','(',i4,',',i4,',',i4,'),',                 &
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m) has excessive change(u)=',e10.3, &
          ' or excessive change(v)=',e10.3)

end subroutine velocity_change
! </SUBROUTINE> NAME="velocity_change"


!#######################################################################
! <SUBROUTINE NAME="compute_topostrophy">
!
! <DESCRIPTION>
!
! Diagnose topostrophy as per Greg Holloway.  
! 
! Stephen.Griffies
! March 2012 
!
! </DESCRIPTION>
!
subroutine compute_topostrophy (Time, Velocity)

  type(ocean_time_type),     intent(in)  :: Time
  type(ocean_velocity_type), intent(in)  :: Velocity

  integer :: i,j,k
  integer :: tau
  real    :: ub, vb

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_velocity_diag_mod(compute_topostrophy): module needs initialization')
  endif 
  tau = Time%tau

  if(id_topostrophy > 0) then 

     wrk1(:,:,:) = 0.0
     if(horz_grid == MOM_BGRID) then 

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                wrk1(i,j,k) = Grd%umask(i,j,k)*Grd%f(i,j) &
                              *(-Velocity%u(i,j,k,2,tau)*Grd%dht_dx(i,j)+Velocity%u(i,j,k,1,tau)*Grd%dht_dy(i,j))
              enddo
           enddo
        enddo

     else ! average c-grid velocity components to vorticity corner         

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                ub = onehalf*(Velocity%u(i,j,k,1,tau)+Velocity%u(i,j+1,k,1,tau))
                vb = onehalf*(Velocity%u(i,j,k,2,tau)+Velocity%u(i+1,j,k,2,tau))
                wrk1(i,j,k) = Grd%umask(i,j,k)*Grd%f(i,j)*(-vb*Grd%dht_dx(i,j)+ub*Grd%dht_dy(i,j))
              enddo
           enddo
        enddo

     endif 

     call diagnose_3d_u(Time, Grd, id_topostrophy, wrk1(:,:,:))
  endif 


end subroutine compute_topostrophy
! </SUBROUTINE>  NAME="compute_topostrophy"



!#######################################################################
! <SUBROUTINE NAME="compute_vorticity">
!
! <DESCRIPTION>
! Compute z-component to vorticity.
! </DESCRIPTION>
!
subroutine compute_vorticity(Time, Velocity)

  type(ocean_time_type),      intent(in)  :: Time
  type(ocean_velocity_type),  intent(in)  :: Velocity

  integer :: i, j, k
  integer :: tau 
  tau   = Time%tau

  if(id_vorticity_z > 0) then 
     wrk1(:,:,:) = 0.0 

     if(horz_grid == MOM_BGRID) then 

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = onehalf*((Velocity%u(i,j,k,2,tau)-Velocity%u(i-1,j,k,2,tau)    )*Grd%dxtnr(i,j)   &
                                       +(Velocity%u(i,j-1,k,2,tau)-Velocity%u(i-1,j-1,k,2,tau))*Grd%dxtnr(i,j-1) &
                                       -(Velocity%u(i,j,k,1,tau)-Velocity%u(i,j-1,k,1,tau)    )*Grd%dyter(i,j)   &
                                       -(Velocity%u(i-1,j,k,1,tau)-Velocity%u(i-1,j-1,k,1,tau))*Grd%dyter(i-1,j) &
                                              )
              enddo
           enddo
        enddo
        call diagnose_3d(Time, Grd, id_vorticity_z, wrk1(:,:,:))

     else ! cgrid 

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = (Velocity%u(i,j,k,2,tau)-Velocity%u(i-1,j,k,2,tau))*Grd%dxur(i,j) & 
                              -(Velocity%u(i,j,k,1,tau)-Velocity%u(i,j-1,k,1,tau))*Grd%dyur(i,j)
              enddo
           enddo
        enddo
        call diagnose_3d_u(Time, Grd, id_vorticity_z, wrk1(:,:,:))

     endif 

  endif ! endif for id_vorticity_z 


end subroutine compute_vorticity
! </SUBROUTINE> NAME="compute_vorticity"


!#######################################################################
! <SUBROUTINE NAME="diagnose_potential_vorticity">
!
! <DESCRIPTION>
! Diagnose pieces to the Ertel potential vorticity, and f * PV.
!
! Motivated in part by the discussion of submeso instabilities in 
! Thomas, Taylor, Ferrari, and Joyce, DSR (2013), page 96-110.
!
! --Assume the Boussinesq form.
!
! --Compute relative vorticity from the full horizontal velocity. 
!   Differences betwen the full relative vorticity and the geostrophic
!   relative vorticity are small for much of the ocean. We prefer to use
!   the full velocity rather than geostrophic velocity, since the diagnosis
!   of the geostrophic relative vorticity requires the Laplacian of pressure,
!   which can be noisy and cumbersome. It is far simpler to just make use
!   of the prognostic u,v fields.
!
! --In contrast, we assume the baroclinic portion of the PV takes on
!   the geostrophically balanced form.  Doing so simplifies the calculation
!   and should be sufficient for many purposes. This assumption is made by 
!   Thomas et al.  But it is an assumption that may need to be revisited
!   if detailed studies are made in the boundary layers, where geostrophy
!   is not well maintained.  
!
! Note: only coded for the B-grid form of relative vorticity.
! 
!   Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine diagnose_potential_vorticity(Time, Velocity, Dens)

  type(ocean_time_type),      intent(in)  :: Time
  type(ocean_velocity_type),  intent(in)  :: Velocity
  type(ocean_density_type),   intent(in)  :: Dens

  real    :: grav_rho0r_sq 
  integer :: i, j, k
  integer :: tau
  
  if(.not. diagnose_pv) return

  tau = Time%tau
  grav_rho0r_sq = grav_rho0r*grav_rho0r

  wrk1(:,:,:)  = 0.0 
  wrk2(:,:,:)  = 0.0 
  wrk3(:,:,:)  = 0.0 
  wrk4(:,:,:)  = 0.0
  wrk5(:,:,:)  = 0.0
  wrk6(:,:,:)  = 0.0


  ! planetary geostrophic PV 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = -grav_rho0r*Dens%drhodz_zt(i,j,k) ! N^2
           wrk2(i,j,k) = wrk1(i,j,k)*coriolis_t(i,j)       ! f*N^2
           wrk3(i,j,k) = wrk2(i,j,k)*coriolis_t(i,j)       ! (f*N)^2 
        enddo
     enddo
  enddo
  if(id_n2_pv  > 0) call diagnose_3d(Time, Grd, id_n2_pv,  wrk1(:,:,:))
  if(id_pg_pv  > 0) call diagnose_3d(Time, Grd, id_pg_pv,  wrk2(:,:,:))
  if(id_pg_pvf > 0) call diagnose_3d(Time, Grd, id_pg_pvf, wrk3(:,:,:))
 

  ! contribution from absolute vorticity 
  if(id_vert_pv > 0 .or. id_vert_pvf > 0 .or. id_pvf > 0) then
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk4(i,j,k) = onehalf*((Velocity%u(i,j,k,2,tau)-Velocity%u(i-1,j,k,2,tau)    )*Grd%dxtnr(i,j)   &
                                    +(Velocity%u(i,j-1,k,2,tau)-Velocity%u(i-1,j-1,k,2,tau))*Grd%dxtnr(i,j-1) &
                                    -(Velocity%u(i,j,k,1,tau)-Velocity%u(i,j-1,k,1,tau)    )*Grd%dyter(i,j)   &
                                    -(Velocity%u(i-1,j,k,1,tau)-Velocity%u(i-1,j-1,k,1,tau))*Grd%dyter(i-1,j) &
                                           )
              wrk5(i,j,k) = (coriolis_t(i,j) + wrk4(i,j,k))*wrk1(i,j,k)  ! (f + zeta)*N^2
              wrk6(i,j,k) = coriolis_t(i,j)*wrk5(i,j,k)                  ! f*(f + zeta)*N^2
           enddo
        enddo
     enddo
     if(id_vert_pv  > 0) call diagnose_3d(Time, Grd, id_vert_pv,  wrk5(:,:,:))
     if(id_vert_pvf > 0) call diagnose_3d(Time, Grd, id_vert_pvf, wrk6(:,:,:))
  endif        

  
  ! contribution from baroclinicity and total contribution 
  if(id_bc_pvf > 0 .or. id_pvf > 0 .or. id_ri_balance > 0) then
     wrk1(:,:,:) = 0.0 
     wrk2(:,:,:) = 0.0 
     wrk4(:,:,:) = 0.0 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = grav_rho0r_sq *                                &
                            (Dens%drhodx_zt(i,j,k)*Dens%drhodx_zt(i,j,k) + &
                             Dens%drhody_zt(i,j,k)*Dens%drhody_zt(i,j,k)   &
                            )
              wrk2(i,j,k) =  wrk6(i,j,k) - wrk1(i,j,k)        ! PV*f = f*(f + zeta)*N^2 - |nabla b|^2 
              wrk4(i,j,k) =  wrk3(i,j,k)/(wrk1(i,j,k)+epsln)  ! Ri_b = (f*N)^2/|nabla b|^2 
           enddo
        enddo
     enddo
     if(id_bc_pvf     > 0) call diagnose_3d(Time, Grd, id_bc_pvf,    -wrk1(:,:,:))
     if(id_pvf        > 0) call diagnose_3d(Time, Grd, id_pvf,        wrk2(:,:,:))
     if(id_ri_balance > 0) call diagnose_3d(Time, Grd, id_ri_balance, wrk4(:,:,:))
  endif 


end subroutine diagnose_potential_vorticity
! </SUBROUTINE> NAME="diagnose_potential_vorticity"




!#######################################################################
! <SUBROUTINE NAME="pressure_conversion">
!
! <DESCRIPTION>
! Perform pressure conversion error analysis.  This analysis should be 
! computed prior to update_ucell_thickness since we need to use dzu(tau)
! and dzten(tau) here rather than dzu(taup1) or dzten(taup1).  
!
! Account taken for Bgrid and Cgrid.  However, blobs have no yet
! been updated for Cgrid.
!
! </DESCRIPTION>
!
subroutine pressure_conversion (Time, Thickness, rho, Dens, Ext_mode, &
                                Adv_vel, Velocity, L_system,          &
                                vert_coordinate_class, use_blobs)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity
  type(ocean_lagrangian_type),    intent(in)    :: L_system
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho
  integer,                        intent(in)    :: vert_coordinate_class 
  logical,                        intent(in)    :: use_blobs

  integer :: tau, taum1, taup1
  integer :: i, j, k, n, kb

  real, dimension(isd:ied,jsd:jed,2)  :: ubar
  real, dimension(isd:ied,jsd:jed,nk) :: rholo_anom
  real, dimension(isd:ied,jsd:jed,nk) :: rhoup_anom
  real :: rholo, rhoup
  real :: column_mass, boxarea
  real :: uint, uext, fx
  real :: term, term1, term2


  if(id_u_dot_grad_pint(1)    > 0 .or. id_u_dot_grad_pint(2)    > 0 .or. &
     id_ubar_dot_grad_pext(1) > 0 .or. id_ubar_dot_grad_pext(2) > 0 .or. &
     id_u_dot_grad_pint_xy    > 0 .or. id_ubar_dot_grad_pext_xy > 0) then 
     call pressure_energy(Time, Thickness, Ext_mode, Velocity, vert_coordinate_class)
  endif 
  
  if (.not. (energy_diag_step > 0 .and. mod(Time%itt, energy_diag_step) == 0)) return 

  call mpp_clock_begin(id_clock_press_conversion)

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  if (use_blobs) then
     ucell_mass = mpp_global_sum(Dom%domain2d, Grd%dau(:,:)*Thickness%mass_uT(:,:,tau), global_sum_flag)
  else
     ucell_mass = mpp_global_sum(Dom%domain2d, Grd%dau(:,:)*Thickness%mass_u(:,:,tau), global_sum_flag)
  endif

  if(horz_grid == MOM_BGRID) then 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_u(i,j,tau)
           mass_column(i,j,2) = Thickness%mass_u(i,j,tau)
        enddo
     enddo
  else 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_en(i,j,1)
           mass_column(i,j,2) = Thickness%mass_en(i,j,2)
        enddo
     enddo
  endif 


 ! truncate out the lower order bits if using non_bitwise_reproducible sums
  if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) ucell_mass = real(ucell_mass, kind=FLOAT_KIND)

  press_conv     = 0.0
  press_int      = 0.0
  press_ext      = 0.0
  press_conv_err = 0.0

  ! place density anomaly in wrk2 since it is needed frequently 
  if (use_blobs) then
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              wrk2(i,j,k) = Grd%tmask(i,j,k)*(rho(i,j,k) - rho0) 
              rholo = (Thickness%dztlo(i,j,k)*rho(i,j,k) + L_system%rho_dztlo(i,j,k))/Thickness%dztloT(i,j,k)
              rhoup = (Thickness%dztup(i,j,k)*rho(i,j,k) + L_system%rho_dztup(i,j,k))/Thickness%dztupT(i,j,k)

              rholo_anom(i,j,k) = Grd%tmask(i,j,k)*(rholo - rho0)
              rhoup_anom(i,j,k) = Grd%tmask(i,j,k)*(rhoup - rho0)
           enddo
        enddo
     enddo
  else
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              wrk2(i,j,k) = Grd%tmask(i,j,k)*(rho(i,j,k) - rho0) 
           enddo
        enddo
     enddo
  endif
  ! place mass tendency in wrk3 as it is needed frequently 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = Grd%tmask(i,j,k)*(Thickness%rho_dzt_tendency(i,j,k)-Thickness%mass_source(i,j,k))
        enddo
     enddo
  enddo

  ! vertically averaged horizontal velocity (meter/sec)
  do j=jsd,jed
     do i=isd,ied
        ubar(i,j,1) = Grd%tmasken(i,j,1,1)*Ext_mode%udrho(i,j,1,tau)/(mass_column(i,j,1)+epsln)
        ubar(i,j,2) = Grd%tmasken(i,j,1,2)*Ext_mode%udrho(i,j,2,tau)/(mass_column(i,j,2)+epsln)
     enddo
  enddo


  ! work done by baroclinic pressure 
  ! press_force [=] Pa [=] kg/(m*s^2)
  ! press_conv  [=] kg*m^2/s^3 = Watt
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea   = grd_area(i,j,k)
              uext      = ubar(i,j,n)
              uint      = Velocity%u(i,j,k,n,tau) - uext  
              term      = Velocity%press_force(i,j,k,n)*boxarea
              press_int = press_int + uint*term
              press_ext = press_ext + uext*term
           enddo
        enddo
     enddo
  enddo

  ! work done by depth independent pressure on depth independent flow 
  ! grad_ps [=] Pa/m =  kg/(m^2*s^2)
  ! grad_mb [=] m/s^2 
  ! press_conv [=] kg*m^2/s^3 = Watt.
  ! since this term also appears in the conversion 
  ! of u dot grad(p), we include it with press_conv.
  if(vert_coordinate_class==DEPTH_BASED) then 
      do j=jsc,jec
         do i=isc,iec
            column_mass = grd_area(i,j,1)*mass_column(i,j,1)
            term        = rho0r*Ext_mode%grad_ps(i,j,1)*ubar(i,j,1)
            press_ext   = press_ext  - term*column_mass
            press_conv  = press_conv - term*column_mass

            column_mass = grd_area(i,j,1)*mass_column(i,j,2)
            term        = rho0r*Ext_mode%grad_ps(i,j,2)*ubar(i,j,2)
            press_ext   = press_ext  - term*column_mass
            press_conv  = press_conv - term*column_mass
         enddo
      enddo
  elseif(vert_coordinate_class==PRESSURE_BASED) then 
      do j=jsc,jec
         do i=isc,iec
            column_mass = grd_area(i,j,1)*mass_column(i,j,1)
            term        = Ext_mode%grad_anompb(i,j,1)*ubar(i,j,1) 
            press_ext   = press_ext  - term*column_mass
            press_conv  = press_conv - term*column_mass

            column_mass = grd_area(i,j,1)*mass_column(i,j,2)
            term        = Ext_mode%grad_anompb(i,j,2)*ubar(i,j,2) 
            press_ext   = press_ext  - term*column_mass
            press_conv  = press_conv - term*column_mass
         enddo
      enddo
  endif


  ! compute the conversion terms 
  if(vert_coordinate_class==DEPTH_BASED) then 

      if (use_blobs) then
         wrk1(:,:,:) = hydrostatic_pressure_blob(Thickness, rhoup_anom(:,:,:), rholo_anom(:,:,:)) 
      else
         wrk1(:,:,:) = hydrostatic_pressure(Thickness, wrk2(:,:,:)) 
      endif

      if(Thickness%method==ENERGETIC) then 

          do j=jsc,jec
             do i=isc,iec
                kb = Grd%kmt(i,j)
                if (kb /= 0) then

                    k = 1
                    term = -rho0r*Grd%tmask(i,j,k)*Grd%dat(i,j)*wrk1(i,j,k)*(wrk3(i,j,k)+Adv_vel%wrho_bt(i,j,k-1))
                    press_conv  = press_conv + term

                    fx = Grd%dat(i,j)*p5grav   
                    do k=2,kb
                       term1 = -fx*Adv_vel%wrho_bt(i,j,k-1)*Thickness%dzwt(i,j,k-1)*(wrk2(i,j,k-1) + wrk2(i,j,k))
                       term2 = -Grd%dat(i,j)*wrk1(i,j,k)*wrk3(i,j,k)
                       press_conv  =  press_conv + rho0r*Grd%tmask(i,j,k)*(term1 + term2)
                    enddo

                endif
             enddo
          enddo

      else 

          do j=jsc,jec
             do i=isc,iec
                kb = Grd%kmt(i,j)
                if (kb /= 0) then

                    k = 1
                    term = -rho0r*Grd%tmask(i,j,k)*Grd%dat(i,j)*wrk1(i,j,k)*(wrk3(i,j,k)+Adv_vel%wrho_bt(i,j,k-1))
                    press_conv  = press_conv + term

                    fx = Grd%dat(i,j)*grav   
                    do k=2,kb
                       term1 = -fx*Adv_vel%wrho_bt(i,j,k-1) &
                       *(Thickness%dztlo(i,j,k-1)*wrk2(i,j,k-1) + Thickness%dztup(i,j,k)*wrk2(i,j,k)) 
                       term2 = -Grd%dat(i,j)*wrk1(i,j,k)*wrk3(i,j,k)
                       press_conv  =  press_conv + rho0r*Grd%tmask(i,j,k)*(term1 + term2)
                    enddo

                endif
             enddo
          enddo

      endif

     ! gradient of geopotential (area sum vanishes on k-level if geodepth_zt is independent of (i,j))
      do k=1,nk
         wrk1_v(:,:,k,1) = BDX_ET((Adv_vel%uhrho_et(:,:,k))*FAX(wrk2(:,:,k)))*Grd%dat(:,:)
         wrk1_v(:,:,k,2) = BDY_NT((Adv_vel%vhrho_nt(:,:,k))*FAY(wrk2(:,:,k)))*Grd%dat(:,:)
         fx = grav*rho0r
         do j=jsc,jec
            do i=isc,iec
               term = -fx*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))*Thickness%geodepth_zt(i,j,k)
               press_conv = press_conv + term
            enddo
         enddo
      enddo

  elseif(vert_coordinate_class==PRESSURE_BASED) then 

     
      if (use_blobs) then
         wrk1(:,:,:) = geopotential_anomaly_blob(Thickness, rhoup_anom(:,:,:), rholo_anom(:,:,:)) 
      else
         wrk1(:,:,:) = geopotential_anomaly(Thickness, wrk2(:,:,:)) 
      endif


      if(Thickness%method==ENERGETIC) then 

          do j=jsc,jec
             do i=isc,iec
                kb = Grd%kmt(i,j)
                if (kb /= 0) then

                    k = 1
                    term = -Grd%tmask(i,j,k)*Grd%dat(i,j)*wrk1(i,j,k)*(wrk3(i,j,k)+Adv_vel%wrho_bt(i,j,k-1))
                    press_conv = press_conv + term

                    fx = Grd%dat(i,j)*p5grav_rho0r    
                    do k=2,kb
                       term1 = -fx*Adv_vel%wrho_bt(i,j,k-1)*Thickness%dzwt(i,j,k-1)*(wrk2(i,j,k-1) + wrk2(i,j,k))
                       term2 = -Grd%dat(i,j)*wrk1(i,j,k)*wrk3(i,j,k)
                       press_conv  =  press_conv + Grd%tmask(i,j,k)*(term1 + term2)
                    enddo

                endif
             enddo
          enddo

      else 

          do j=jsc,jec
             do i=isc,iec
                kb = Grd%kmt(i,j)
                if (kb /= 0) then

                    k = 1
                    term = -Grd%tmask(i,j,k)*Grd%dat(i,j)*wrk1(i,j,k)*(wrk3(i,j,k)+Adv_vel%wrho_bt(i,j,k-1))
                    press_conv = press_conv + term

                    fx = Grd%dat(i,j)*grav_rho0r    
                    do k=2,kb
                       term1 = -fx*Adv_vel%wrho_bt(i,j,k-1) &
                                *(Thickness%dztlo(i,j,k-1)*wrk2(i,j,k-1) + &
                                  Thickness%dztup(i,j,k)  *wrk2(i,j,k)) 
                       term2 = -Grd%dat(i,j)*wrk1(i,j,k)*wrk3(i,j,k)
                       press_conv  =  press_conv + Grd%tmask(i,j,k)*(term1 + term2)
                    enddo

                endif
             enddo
          enddo

      endif

      ! now for the "unmanipulated" piece 
      if(horz_grid == MOM_BGRID) then 

         do k=1,nk
            do j=jsd,jed
               do i=isd,iec
                  wrk2_v(i,j,k,1) = (Dens%pressure_at_depth(i+1,j,k)-Dens%pressure_at_depth(i,j,k))*dbars2Pa
               enddo
            enddo
            do j=jsd,jec
               do i=isd,ied
                  wrk2_v(i,j,k,2) = (Dens%pressure_at_depth(i,j+1,k)-Dens%pressure_at_depth(i,j,k))*dbars2Pa
               enddo
            enddo
         enddo
         do k=1,nk
            wrk1_v(:,:,k,1) = FAY(FAX(wrk2(:,:,k))*wrk2_v(:,:,k,1))*Grd%dyu(:,:)
            wrk1_v(:,:,k,2) = FAX(FAY(wrk2(:,:,k))*wrk2_v(:,:,k,2))*Grd%dxu(:,:)
            do j=jsc,jec
               do i=isc,iec
                  term = rho0r*Thickness%dzu(i,j,k) &
                  *(Velocity%u(i,j,k,1,tau)*wrk1_v(i,j,k,1) + Velocity%u(i,j,k,2,tau)*wrk1_v(i,j,k,2))
                  press_conv = press_conv + term
               enddo
            enddo
         enddo

      else   ! Cgrid 

         do k=1,nk
            do j=jsd,jed
               do i=isd,iec
                  wrk2_v(i,j,k,1) = Grd%dxter(i,j)*(Dens%pressure_at_depth(i+1,j,k)-Dens%pressure_at_depth(i,j,k))*dbars2Pa
               enddo
            enddo
            do j=jsd,jec
               do i=isd,ied
                  wrk2_v(i,j,k,2) = Grd%dytnr(i,j)*(Dens%pressure_at_depth(i,j+1,k)-Dens%pressure_at_depth(i,j,k))*dbars2Pa
               enddo
            enddo
         enddo
         do k=1,nk
            wrk1_v(:,:,k,1) = FAX(wrk2(:,:,k))*wrk2_v(:,:,k,1)
            wrk1_v(:,:,k,2) = FAY(wrk2(:,:,k))*wrk2_v(:,:,k,2)
            do j=jsc,jec
               do i=isc,iec
                  term = rho0r*Thickness%dzten(i,j,k,1)*Grd%dat(i,j)*Velocity%u(i,j,k,1,tau)*wrk1_v(i,j,k,1) &
                        +rho0r*Thickness%dzten(i,j,k,2)*Grd%dat(i,j)*Velocity%u(i,j,k,2,tau)*wrk1_v(i,j,k,2)
                  press_conv = press_conv + term
               enddo
            enddo
         enddo

     endif ! endif for horz_grid 


  endif   ! endif for vertical coordinates  


  call mpp_sum(press_int)
  call mpp_sum(press_ext)
  call mpp_sum(press_conv)
  press_conv_err = press_int + press_ext - press_conv 


  call mpp_clock_end(id_clock_press_conversion)


end subroutine pressure_conversion
! </SUBROUTINE> NAME="pressure_conversion"


!#######################################################################
! <SUBROUTINE NAME="pressure_energy">
!
! <DESCRIPTION>
! Diagnose u dot grad(p) for diagnostic purposes. These maps 
! when summed over all grid points will result in an energy 
! that is equal to pint+pext as computed in pressure_conversion.  
!
! Account taken of either Bgrid or Cgrid. 
!
! </DESCRIPTION>
!
subroutine pressure_energy (Time, Thickness, Ext_mode, Velocity, vert_coordinate_class)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity
  integer,                        intent(in)    :: vert_coordinate_class 

  integer :: tau
  integer :: i, j, k, n

  real, dimension(isd:ied,jsd:jed,2) :: ubar
  real :: column_mass, boxarea
  real :: term

  call mpp_clock_begin(id_clock_press_energy)

  tau = Time%tau
  wrk1_v(:,:,:,:) = 0.0
  wrk1_v2d(:,:,:) = 0.0

  if(horz_grid == MOM_BGRID) then 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_u(i,j,tau)
           mass_column(i,j,2) = Thickness%mass_u(i,j,tau)
        enddo
     enddo
  else 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_en(i,j,1)
           mass_column(i,j,2) = Thickness%mass_en(i,j,2)
        enddo
     enddo
  endif 

  ! vertically averaged horizontal velocity 
  do j=jsd,jed
     do i=isd,ied
        ubar(i,j,1) = Grd%tmasken(i,j,1,1)*Ext_mode%udrho(i,j,1,tau)/(mass_column(i,j,1)+epsln)
        ubar(i,j,2) = Grd%tmasken(i,j,1,2)*Ext_mode%udrho(i,j,2,tau)/(mass_column(i,j,2)+epsln)
     enddo
  enddo


  ! work done by horizontal baroclinic gradients
  ! press_force [=] Pa [=] kg/(m*s^2)
  ! wrk1_v [=] wrk2_v [=] kg*m^2/s^3 = Watt
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea         =  grd_area(i,j,k)
              term            =  Velocity%press_force(i,j,k,n)*boxarea
              wrk1_v(i,j,k,n) =  Velocity%u(i,j,k,n,tau)*term
           enddo
        enddo
     enddo
  enddo

  ! work done by depth independent pressure on depth independent flow 
  ! grad_ps [=] Pa/m =  kg/(m^2*s^2)
  ! grad_mb [=] m/s^2 
  ! wrk2_v [=] kg*m^2/s^3 = Watt
  if(vert_coordinate_class==DEPTH_BASED) then 
      do j=jsc,jec
         do i=isc,iec
            column_mass     = grd_area(i,j,1)*mass_column(i,j,1)
            term            = rho0r*Ext_mode%grad_ps(i,j,1)*ubar(i,j,1) 
            wrk1_v2d(i,j,1) =  -term*column_mass
            column_mass     = grd_area(i,j,1)*mass_column(i,j,2)
            term            = rho0r*Ext_mode%grad_ps(i,j,2)*ubar(i,j,2)
            wrk1_v2d(i,j,2) =  -term*column_mass
         enddo
      enddo
  elseif(vert_coordinate_class==PRESSURE_BASED) then 
      do j=jsc,jec
         do i=isc,iec
            column_mass     = grd_area(i,j,1)*mass_column(i,j,1)
            term            = Ext_mode%grad_anompb(i,j,1)*ubar(i,j,1) 
            wrk1_v2d(i,j,1) = -term*column_mass
            column_mass     = grd_area(i,j,1)*mass_column(i,j,2)
            term            = Ext_mode%grad_anompb(i,j,2)*ubar(i,j,2)            
            wrk1_v2d(i,j,2) = -term*column_mass
         enddo
      enddo
  endif

  
  if (id_u_dot_grad_pint(1) > 0) used = send_data (id_u_dot_grad_pint(1), wrk1_v(:,:,:,1), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                             &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_u_dot_grad_pint(2) > 0) used = send_data (id_u_dot_grad_pint(2), wrk1_v(:,:,:,2), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                             &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_u_dot_grad_pint_xy > 0) used = send_data (id_u_dot_grad_pint_xy,             &
                       wrk1_v(:,:,:,1)+wrk1_v(:,:,:,2),                               &
                       Time%model_time, rmask=Grd%mask(:,:,:),                        &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  if (id_ubar_dot_grad_pext(1) > 0) used = send_data (id_ubar_dot_grad_pext(1), wrk1_v2d(:,:,1), &
                       Time%model_time, rmask=Grd%mask(:,:,1),                                   &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_ubar_dot_grad_pext(2) > 0) used = send_data (id_ubar_dot_grad_pext(2), wrk1_v2d(:,:,2), &
                       Time%model_time, rmask=Grd%mask(:,:,1),                                   &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_ubar_dot_grad_pext_xy > 0) used = send_data (id_ubar_dot_grad_pext_xy, &
                       wrk1_v2d(:,:,1)+wrk1_v2d(:,:,2),                         &
                       Time%model_time, rmask=Grd%mask(:,:,1),                  &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  call mpp_clock_end(id_clock_press_energy)

end subroutine pressure_energy
! </SUBROUTINE> NAME="pressure_energy"


!#######################################################################
! <SUBROUTINE NAME="friction_energy">
!
! <DESCRIPTION>
! Diagnose u dot Friction for diagnostic purposes. 
!
! Account taken for either Bgrid or Cgrid.  
!
! NOTE:
!
! A) DO NOT split into baroclinic and barotropic pieces. Just compute 
! u dot F using full velocity field u. Otherwise, the calculation emulates 
! that done in energy_analysis subroutine. 
!
! B) DO NOT remove the effects from bottom drag and from surface stress. 
! The reason is that bottom drag and surface stress are incorporated to 
! the vertical friction operator, even when doing vertical friction 
! implicitly in time.  So it is tough to remove these effects in an 
! explicit diagnostic manner.  So the u dot vertical friction piece
! includes BOTH surface (i.e., winds) and bottom stress.
!
! </DESCRIPTION>
!
subroutine friction_energy (Time, Thickness, Adv_vel, Ext_mode, Velocity, visc_cbu, visc_cbt)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbt

  integer :: tau
  integer :: i, j, k, n

  real, dimension(isd:ied,jsd:jed,2) :: ubar
  real :: boxarea, term

  call mpp_clock_begin(id_clock_friction_energy)

  tau = Time%tau
  wrk2_v(:,:,:,:) = 0.0   ! for holding u dot friction 
                          ! (do NOT use wrk1_v, as that is used in friction modules)
  wrk1_v2d(:,:,:) = 0.0   ! for holding ubar dot horizontal stress acting JUST on ubar 

  if(horz_grid == MOM_BGRID) then 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_u(i,j,tau)
           mass_column(i,j,2) = Thickness%mass_u(i,j,tau)
        enddo
     enddo
  else 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_en(i,j,1)
           mass_column(i,j,2) = Thickness%mass_en(i,j,2)
        enddo
     enddo
  endif 

  ! vertically averaged horizontal velocity 
  do j=jsd,jed
     do i=isd,ied
        ubar(i,j,1) = Grd%tmasken(i,j,1,1)*Ext_mode%udrho(i,j,1,tau)/(mass_column(i,j,1)+epsln)
        ubar(i,j,2) = Grd%tmasken(i,j,1,2)*Ext_mode%udrho(i,j,2,tau)/(mass_column(i,j,2)+epsln)
     enddo
  enddo

  ! work done by laplacian friction on horizontal velocity field.
  ! lap_friction returns thickness weighted and density weighted
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! wrk1_v [=]  kg*m^2/s^3 = Watt
  call lap_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.true.)
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea         = grd_area(i,j,k)
              term            = Velocity%wrkv(i,j,k,n)*boxarea    
              wrk2_v(i,j,k,n) = Velocity%u(i,j,k,n,tau)*term  
           enddo
        enddo
     enddo
  enddo

  if (id_u_dot_lapfrict_horz > 0) used = send_data (id_u_dot_lapfrict_horz, wrk2_v(:,:,:,1), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                               &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_v_dot_lapfrict_horz > 0) used = send_data (id_v_dot_lapfrict_horz, wrk2_v(:,:,:,2), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                               &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  ! work done by biharmonic friction
  ! bih_friction returns thickness weighted and density weighted
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa.
  ! include contributions from both biharmonic friction and 
  ! side drag friction that may optionally be applied.   
  ! eng(5) [=]  kg*m^2/s^3 = Watt
  call bih_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.true.)
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea         = grd_area(i,j,k)
              term            = Velocity%wrkv(i,j,k,n)*boxarea    
              wrk2_v(i,j,k,n) = Velocity%u(i,j,k,n,tau)*term  
           enddo
        enddo
     enddo
  enddo

  if (id_u_dot_bihfrict_horz > 0) used = send_data (id_u_dot_bihfrict_horz, wrk2_v(:,:,:,1), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                               &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_v_dot_bihfrict_horz > 0) used = send_data (id_v_dot_bihfrict_horz, wrk2_v(:,:,:,2), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                               &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  ! contribution from lap_friction acting just on barotropic. 
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! wrk1_v2d [=]  kg*m^2/s^3 = Watt
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           boxarea         = grd_area(i,j,1)
           term            = Velocity%lap_friction_bt(i,j,n)*boxarea    
           wrk1_v2d(i,j,n) = ubar(i,j,n)*term
        enddo
     enddo
  enddo

  ! contribution from bih_friction acting just on the barotropic. 
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! eng(5) [=]  kg*m^2/s^3 = Watt
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           boxarea         = grd_area(i,j,1)
           term            = Velocity%bih_friction_bt(i,j,n)*boxarea    
           wrk1_v2d(i,j,n) = wrk1_v2d(i,j,n) + ubar(i,j,n)*term
        enddo
     enddo
  enddo

  if (id_ubar_dot_frict_bt > 0) used = send_data (id_ubar_dot_frict_bt, wrk1_v2d(:,:,1), &
                       Time%model_time, rmask=Grd%mask(:,:,1),                           &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_vbar_dot_frict_bt > 0) used = send_data (id_vbar_dot_frict_bt, wrk1_v2d(:,:,2), &
                       Time%model_time, rmask=Grd%mask(:,:,1),                           &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  ! work done by vertical friction, including wind stress and bottom stress. 
  ! implicit acceleration has been saved from previous call to ocean_implicit_accel.
  ! the implicit piece includes contributions from visc_cbu_form_drag, as well as visc_cbu, visc_cbt. 
  ! vert_frict returns thickness weighted and density weighted vertical friction (kg/m^3)*(m^2/s^2)
  ! wrk2_v [=]  kg*m^2/s^3 = Watt
  ! Note: no need to call the implicit operator again, since saved result in vfrict_impl. 
  if(horz_grid == MOM_BGRID) then 
     call vert_friction_bgrid(Time, Thickness, Velocity, visc_cbu, energy_analysis_step=.true.)
  else 
     call vert_friction_cgrid(Time, Thickness, Velocity, visc_cbt, energy_analysis_step=.true.)
  endif 
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea         = grd_area(i,j,k)
              term            = (Velocity%wrkv(i,j,k,n)+ Velocity%vfrict_impl(i,j,k,n))*boxarea 
              wrk2_v(i,j,k,n) = Velocity%u(i,j,k,n,tau)*term     
           enddo
        enddo
     enddo
  enddo
 
  if (id_u_dot_frict_vert > 0) used = send_data (id_u_dot_frict_vert, wrk2_v(:,:,:,1), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                         &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_v_dot_frict_vert > 0) used = send_data (id_v_dot_frict_vert, wrk2_v(:,:,:,2), &
                       Time%model_time, rmask=Grd%mask(:,:,:),                         &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  call mpp_clock_end(id_clock_friction_energy)

end subroutine friction_energy
! </SUBROUTINE> NAME="friction_energy"


!#######################################################################
! <SUBROUTINE NAME="vert_dissipation">
!
! <DESCRIPTION>
! Diagnose dissipation from vertical friction due just to viscosity.
!
! Units W/m^2
!
! Assumptions:
!
! 1/ Ignore bottom drag here...just concerned with viscosity. 
!
! 2/ Assume vertical friction is handled implicitly in time. 
! 
! </DESCRIPTION>
!
subroutine vert_dissipation (Time, Thickness, Velocity, visc_cbu, visc_cbt, visc_cbu_form_drag)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_velocity_type),      intent(in) :: Velocity
  real, dimension(isd:,jsd:,:),   intent(in) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(in) :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(in) :: visc_cbu_form_drag

  integer :: itime
  integer :: i, j, k, n
  real    :: du_dz 
  type(time_type) :: next_time 

  call mpp_clock_begin(id_clock_vert_dissipation)

  itime = Time%taup1

  next_time = increment_time(Time%model_time, int(dtime), 0)
  if (need_data(id_vert_lap_diss,next_time)) then 

      if(horz_grid == MOM_BGRID) then 

         wrk1(:,:,:) = 0.0
         do n=1,2
            do k=1,nk-1
               do j=jsc,jec
                  do i=isc,iec
                     du_dz = Grd%umask(i,j,k+1) &
                     *(Velocity%u(i,j,k,n,itime)-Velocity%u(i,j,k+1,n,itime))/Thickness%dzwu(i,j,k)
                     wrk1(i,j,k) = wrk1(i,j,k) &   
                     -rho0*Thickness%dzu(i,j,k)*du_dz*du_dz &
                      *(visc_cbu(i,j,k)+visc_cbu_form_drag(i,j,k,n))  
                  enddo
               enddo
            enddo
         enddo
         call diagnose_3d_u(Time, Grd, id_vert_lap_diss, wrk1(:,:,:))
      
      else  ! Cgrid

         wrk1(:,:,:) = 0.0
         do n=1,2
            do k=1,nk-1
               do j=jsc,jec
                  do i=isc,iec
                     du_dz = Grd%tmasken(i,j,k+1,n) &
                     *(Velocity%u(i,j,k,n,itime)-Velocity%u(i,j,k+1,n,itime))/(epsln+Thickness%dzten(i,j,k,n))
                     wrk1(i,j,k) = wrk1(i,j,k) &   
                     -rho0*Thickness%dzten(i,j,k,n)*du_dz*du_dz &
                      *(visc_cbt(i,j,k)+visc_cbu_form_drag(i,j,k,n))  
                  enddo
               enddo
            enddo
         enddo
         call diagnose_3d(Time, Grd, id_vert_lap_diss, wrk1(:,:,:))

      endif 

  endif

  call mpp_clock_end(id_clock_vert_dissipation)

end subroutine vert_dissipation
! </SUBROUTINE> NAME="vert_dissipation"



!#######################################################################
! <SUBROUTINE NAME="energy_analysis">
!
! <DESCRIPTION>
! Perform energy analysis by taking scalar product of horizontal
! velocity with the velocity equations and integrating over the ocean volume.
!
! Pressure conversions have already been computed in pressure_conversion 
! subroutine.  It is necessary to perform that analysis earlier than 
! the call to update_ucell_thickness inside ocean_model_mod, whereas 
! the remaining elements in the energy analysis can be called at the 
! end of the update for velocity.  
! </DESCRIPTION>
!
subroutine energy_analysis (Time, Thickness, Ext_mode, Adv_vel, Dens,    & 
                            pme, river, upme, uriver, visc_cbu, visc_cbt,&
                            visc_cbu_form_drag, Velocity)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river 
  real, dimension(isd:,jsd:,:),   intent(in)    :: upme
  real, dimension(isd:,jsd:,:),   intent(in)    :: uriver
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(in)    :: visc_cbu_form_drag

  integer :: i, j, k, n
  integer :: tau, taum1, taup1, tau_m0

  real :: boxarea
  real :: uint, uext, term, vel
  real :: rho_dz_vel_tau
  real :: ke_tot, pe_tot, pe_tot_rel, ke_per_mass
  
  real, dimension(isd:ied,jsd:jed,2) :: ubar
  real, dimension(isd:ied,jsd:jed)   :: pme_u
  real, dimension(isd:ied,jsd:jed)   :: river_u

  real, dimension(13) :: engint         ! internal mode energy integral components (Watt)
  real, dimension(13) :: engext         ! external mode energy integral components (Watt)
  real                :: compression    ! power from compression of grid cell and/or mass sources  
  real                :: enleak         ! difference between integrated advection of momentum 
                                        ! and power from compression and water forcing 
  integer :: stdoutunit 
  stdoutunit=stdout() 


  ! get 3d maps for energetics from friction 
  if(id_u_dot_lapfrict_horz  > 0 .or. id_v_dot_lapfrict_horz  > 0 .or. &
     id_u_dot_bihfrict_horz  > 0 .or. id_v_dot_bihfrict_horz  > 0 .or. &
     id_u_dot_frict_vert     > 0 .or. id_v_dot_frict_vert     > 0 .or. &
     id_ubar_dot_frict_bt    > 0 .or. id_vbar_dot_frict_bt    > 0) then 
     call friction_energy(Time, Thickness, Adv_vel, Ext_mode, Velocity, visc_cbu, visc_cbt)
  endif 

  ! get energy dissipation from vertical friction 
  call vert_dissipation(Time, Thickness, Velocity, visc_cbu, visc_cbt, visc_cbu_form_drag)

  if (.not. (energy_diag_step > 0 .and. mod(Time%itt, energy_diag_step) == 0)) return 

  call mpp_clock_begin(id_clock_energy_analysis)

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_velocity_diag_mod (energy_analysis): module must be initialized')
  endif 

  tau    = Time%tau
  taum1  = Time%taum1
  taup1  = Time%taup1
  tau_m0 = Time%tau_m0

  compression = 0.0
  engint(:)   = 0.0
  engext(:)   = 0.0

  ! place mass tendency in wrk3 as it is needed frequently 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = Grd%tmask(i,j,k)*(Thickness%rho_dzt_tendency(i,j,k)-Thickness%mass_source(i,j,k))
        enddo
     enddo
  enddo

  ! mass per horizontal area of a fluid column
  if(horz_grid == MOM_BGRID) then 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_u(i,j,tau)
           mass_column(i,j,2) = Thickness%mass_u(i,j,tau)
        enddo
     enddo
  else 
     do j=jsd,jed
        do i=isd,ied
           mass_column(i,j,1) = Thickness%mass_en(i,j,1)
           mass_column(i,j,2) = Thickness%mass_en(i,j,2)
        enddo
     enddo
  endif 

  ! compute vertical average of horizontal velocity
  do j=jsd,jed
    do i=isd,ied
      ubar(i,j,1) = Grd%tmasken(i,j,1,1)*Ext_mode%udrho(i,j,1,tau)/(mass_column(i,j,1)+epsln)
      ubar(i,j,2) = Grd%tmasken(i,j,1,2)*Ext_mode%udrho(i,j,2,tau)/(mass_column(i,j,2)+epsln)
    enddo
  enddo

  if (debug_this_module) then
      write(stdoutunit,*) ' '
      write(stdoutunit,*) '===chksums in energy analysis===' 
      write(stdoutunit,*) 'dtime(seconds) = ',dtime, 'and dtimer = ',dtimer

      call write_chksum_3d('Thickness%rho_dzu(taum1)', Thickness%rho_dzu(COMP,:,taum1)*Grd%umask(COMP,:))
      call write_chksum_3d('Thickness%rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('Thickness%rho_dzu(taup1)', Thickness%rho_dzu(COMP,:,taup1)*Grd%umask(COMP,:))
      call write_chksum_3d('Thickness%rho_dzten(1)', Thickness%rho_dzten(COMP,:,1)*Grd%tmasken(COMP,:,1))
      call write_chksum_3d('Thickness%rho_dzten(2)', Thickness%rho_dzten(COMP,:,2)*Grd%tmasken(COMP,:,2))
      call write_chksum_2d('mass_column(1)', mass_column(COMP,1)*Grd%umask(COMP,1))
      call write_chksum_2d('mass_column(2)', mass_column(COMP,2)*Grd%umask(COMP,1))
      call write_chksum_3d('Velocity%u(n=1,taum1)', Velocity%u(COMP,:,1,taum1)*Grd%tmasken(COMP,:,1))
      call write_chksum_3d('Velocity%u(n=1,tau)', Velocity%u(COMP,:,1,tau)*Grd%tmasken(COMP,:,1))
      call write_chksum_3d('Velocity%u(n=1,taup1)', Velocity%u(COMP,:,1,taup1)*Grd%tmasken(COMP,:,1))
      call write_chksum_3d('Velocity%u(n=2,taum1)', Velocity%u(COMP,:,2,taum1)*Grd%tmasken(COMP,:,2))
      call write_chksum_3d('Velocity%u(n=2,tau)', Velocity%u(COMP,:,2,tau)*Grd%tmasken(COMP,:,2))
      call write_chksum_3d('Velocity%u(n=2,taup1)', Velocity%u(COMP,:,2,taup1)*Grd%tmasken(COMP,:,2))
      call write_chksum_2d('ubar(n=1)', ubar(COMP,1)*Grd%tmasken(COMP,1,1))
      call write_chksum_2d('ubar(n=2)', ubar(COMP,2)*Grd%tmasken(COMP,1,2))

      write(stdoutunit,*)'================================='
  endif


  ! scalar product of u with the time rate of change 
  ! eng(1) [=] kg*m^2/s^3 = Watt
  if(horz_grid == MOM_BGRID) then 

     do n=1,2
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 boxarea   = grd_area(i,j,k)
                 uext      = ubar(i,j,n)
                 uint      = Velocity%u(i,j,k,n,tau) - uext
                 term      = (Thickness%rho_dzu(i,j,k,taup1)*Velocity%u(i,j,k,n,taup1)  &
                             -Thickness%rho_dzu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1)) &
                             *dtimer*boxarea
                 engint(1) = engint(1) + uint*term
                 engext(1) = engext(1) + uext*term
              enddo
           enddo
        enddo
     enddo

  else   ! Cgrid 

     do n=1,2
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 boxarea        = grd_area(i,j,k)*Grd%tmasken(i,j,k,n)
                 uext           = ubar(i,j,n)
                 uint           = Velocity%u(i,j,k,n,tau) - uext
                 rho_dz_vel_tau = (2-n)*Adv_vel%uhrho_et(i,j,k)*Grd%dxte_dxtr(i,j)*Grd%dyte_dytr(i,j) &
                                 +(n-1)*Adv_vel%vhrho_nt(i,j,k)*Grd%dxtn_dxtr(i,j)*Grd%dytn_dytr(i,j)
                 term      = (Thickness%rho_dzten(i,j,k,n)*Velocity%u(i,j,k,n,taup1) - rho_dz_vel_tau)*dtimer*boxarea
                 engint(1) = engint(1) + uint*term
                 engext(1) = engext(1) + uext*term
              enddo
           enddo
        enddo
     enddo 

  endif   ! horz_grid endif 


  ! work done by horizontal advection 
  ! horz_advection_of_velocity [=] thickness and density weighted advection 
  ! which has dimensions (kg/m^3)*(m^2/s^2) = Pa
  ! eng(2) [=] kg*m^2/s^3 = Watt
  call horz_advection_of_velocity(Time, Thickness, Adv_vel, Velocity, &
                                  energy_analysis_step=.true.)       
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd_area(i,j,k)
          uext      = ubar(i,j,n)
          uint      = Velocity%u(i,j,k,n,tau) - uext 
          term      = Velocity%wrkv(i,j,k,n)*boxarea  
          engint(2) = engint(2) - uint*term
          engext(2) = engext(2) - uext*term
        enddo
      enddo
    enddo
  enddo

  ! work done by vertical advection
  ! vertical_advection_of_velocity [=] thickness and density weighted advection
  ! which as dimensions (kg/m^3)*(m^2/s^2) = Pa
  ! eng(3) [=]  kg*m^2/s^3 = Watt
  call vert_advection_of_velocity(Time, Adv_vel, Velocity,  &
                                  pme, river, upme, uriver, &
                                  energy_analysis_step=.true.) 
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea   = grd_area(i,j,k)
              uext      = ubar(i,j,n)
              uint      = Velocity%u(i,j,k,n,tau) - uext   
              term      = Velocity%wrkv(i,j,k,n)*boxarea
              engint(3) = engint(3) - uint*term
              engext(3) = engext(3) - uext*term
           enddo
        enddo
     enddo
  enddo

  ! work done by laplacian friction
  ! lap_friction returns thickness weighted and density weighted
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! eng(4) [=]  kg*m^2/s^3 = Watt
  call lap_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.true.)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = grd_area(i,j,k)
          uext      = ubar(i,j,n)
          uint      = Velocity%u(i,j,k,n,tau) - uext
          term      = Velocity%wrkv(i,j,k,n)*boxarea    
          engint(4) = engint(4) + uint*term
          engext(4) = engext(4) + uext*term
        enddo
      enddo
    enddo
  enddo

  ! add contribution from lap_friction acting just on the barotropic. 
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! eng(4) [=]  kg*m^2/s^3 = Watt
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           boxarea   = grd_area(i,j,1)
           uext      = ubar(i,j,n)
           term      = Velocity%lap_friction_bt(i,j,n)*boxarea    
           engext(4) = engext(4) + uext*term
        enddo
     enddo
  enddo

  ! work done by biharmonic friction
  ! lap_friction returns thickness weighted and density weighted
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! eng(5) [=]  kg*m^2/s^3 = Watt
  call bih_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step=.true.)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = grd_area(i,j,k)
          uext      = ubar(i,j,n)
          uint      = Velocity%u(i,j,k,n,tau) - uext
          term      = Velocity%wrkv(i,j,k,n)*boxarea    
          engint(5) = engint(5) + uint*term
          engext(5) = engext(5) + uext*term
        enddo
      enddo
    enddo
  enddo

  ! add contribution from bih_friction acting just on the barotropic. 
  ! friction with dimensions (kg/m^3)*(m^2/s^2) = N/m^2 = Pa
  ! eng(5) [=]  kg*m^2/s^3 = Watt
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           boxarea   = grd_area(i,j,1)
           uext      = ubar(i,j,n)
           term      = Velocity%bih_friction_bt(i,j,n)*boxarea    
           engext(5) = engext(5) + uext*term
        enddo
     enddo
  enddo

  ! work done by vertical friction (contributions from wind and bottom drag removed below)
  ! implicit acceleration has been saved from previous call to ocean_implicit_accel.
  ! Greatbatch form drag implemented through vertical viscosity is part of implict accel. 
  ! vert_frict returns thickness weighted and density weighted vertical friction (kg/m^3)*(m^2/s^2)
  ! eng(6) [=]  kg*m^2/s^3 = Watt
  if(horz_grid == MOM_BGRID) then 
     call vert_friction_bgrid(Time, Thickness, Velocity, visc_cbu, energy_analysis_step=.true.)
  else 
     call vert_friction_bgrid(Time, Thickness, Velocity, visc_cbt, energy_analysis_step=.true.)
  endif 
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea   = grd_area(i,j,k)
              uext      = ubar(i,j,n)
              uint      = Velocity%u(i,j,k,n,tau) - uext  
              term      = (Velocity%wrkv(i,j,k,n)+ Velocity%vfrict_impl(i,j,k,n))*boxarea 
              engint(6) = engint(6) + uint*term
              engext(6) = engext(6) + uext*term
           enddo
        enddo
     enddo
  enddo


  ! work done by wind stress and bottom drag 
  ! smf [=] N/m^2 = Pa 
  ! bmf [=] N/m^2 = Pa 
  ! eng(7) [=] kg*m^2/s^3 = Watt
  do n=1,2

    ! wind stress 
    k = 1
    if(horz_grid == MOM_BGRID) then 
       do j=jsc,jec
          do i=isc,iec
             uext = ubar(i,j,n)
             uint = Velocity%u(i,j,k,n,tau) - uext
             term = grd_area(i,j,k)*Velocity%smf_bgrid(i,j,n)  
             engint(7) = engint(7) + uint*term    
             engext(7) = engext(7) + uext*term
          enddo
       enddo
    else 
       do j=jsc,jec
          do i=isc,iec
             uext = ubar(i,j,n)
             uint = Velocity%u(i,j,k,n,tau) - uext
             term = grd_area(i,j,k)*Velocity%smf_cgrid(i,j,n)  
             engint(7) = engint(7) + uint*term    
             engext(7) = engext(7) + uext*term
          enddo
       enddo
    endif 

    ! bottom stress
    do j=jsc,jec
       do i=isc,iec
          k = Grd%kmu(i,j)
          if (k /= 0) then
             uext = ubar(i,j,n)
             uint = Velocity%u(i,j,k,n,tau) - uext    
             term = grd_area(i,j,k)*Velocity%bmf(i,j,n)
             engint(8) = engint(8) - uint*term    
             engext(8) = engext(8) - uext*term
          endif
       enddo
    enddo

  enddo  ! n=1,2

  ! subtract wind stress and bottom stress contributions
  ! from vertical friction to facilitate diagnostics 
   engint(6) = engint(6) - engint(7) - engint(8) 
   engext(6) = engext(6) - engext(7) - engext(8) 


  ! work due to Coriolis force (zero on B_grid only when acor=0.0)
  ! coriolis returns thickness weighted and density weighted Coriolis force (kg/m^3)*(m^2/s^2) 
  ! eng(9) [=]  kg*m^2/s^3 = Watt
  if(horz_grid == MOM_BGRID) then 
     call coriolis_force_bgrid(Time, Thickness, Velocity, energy_analysis_step=.true.)
  else 
     call coriolis_force_cgrid(Time, Adv_vel, Velocity, abtau_m0, abtau_m1, abtau_m2, &
                               energy_analysis_step=.true.)
  endif 
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea   = grd_area(i,j,k)
              uext      = ubar(i,j,n)
              uint      = Velocity%u(i,j,k,n,tau) - uext  
              term      = Velocity%wrkv(i,j,k,n)*boxarea 
              engint(9) = engint(9) + uint*term   
              engext(9) = engext(9) + uext*term   
           enddo
        enddo
     enddo
  enddo

  ! power (Watt) from surface water fluxes entering (+power) or leaving (-power) ocean 
  if(horz_grid == MOM_BGRID) then 
     pme_u(:,:)   = REMAP_BT_TO_BU(pme)
     river_u(:,:) = REMAP_BT_TO_BU(river)
  else 
     pme_u(:,:)   = pme(:,:)
     river_u(:,:) = river(:,:)
  endif 
  k = 1
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           boxarea    = grd_area(i,j,k)
           vel        = Velocity%u(i,j,k,n,tau)
           uext       = ubar(i,j,n)
           uint       = vel - uext   
           term       = pme_u(i,j)*upme(i,j,n) + river_u(i,j)*uriver(i,j,n)
           term       = term*boxarea 
           engint(10) = engint(10) + uint*term
           engext(10) = engext(10) + uext*term
        enddo
     enddo
  enddo

  ! contribution (Watt) to kinetic energy conversion from surface water fluxes 
  k = 1
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           boxarea    = grd_area(i,j,k)
           vel        = Velocity%u(i,j,k,n,tau)
           uext       = ubar(i,j,n)
           uint       = vel - uext   
           term       = pme_u(i,j)*upme(i,j,n)+river_u(i,j)*uriver(i,j,n)
           term       = 0.5*term*boxarea 
           engint(11) = engint(11) + uint*term
           engext(11) = engext(11) + uext*term
        enddo
     enddo
  enddo

  ! contribution (Watt) to kinetic energy from momentum sources 
  call momentum_source(Time, Thickness, Velocity, energy_analysis_step=.true.)
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea    = grd_area(i,j,k)
              uext       = ubar(i,j,n)
              uint       = Velocity%u(i,j,k,n,tau) - uext
              term       = Velocity%wrkv(i,j,k,n)*boxarea    
              engint(12) = engint(12) + uint*term
              engext(12) = engext(12) + uext*term
           enddo
        enddo
     enddo
  enddo

  ! contribution (Watt) to kinetic energy from parameterized form drag via Aiki scheme 
  call form_drag_accel(Time, Thickness, Velocity, energy_analysis_step=.true.)
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea    = grd_area(i,j,k)
              uext       = ubar(i,j,n)
              uint       = Velocity%u(i,j,k,n,tau) - uext
              term       = Velocity%wrkv(i,j,k,n)*boxarea    
              engint(13) = engint(13) + uint*term
              engext(13) = engext(13) + uext*term
           enddo
        enddo
     enddo
  enddo

  ! rho_dzu_tendency and mass source are converted from kinetic energy advection.
  ! these terms arise from the ability of the grid cell volume/mass to change due 
  ! to changes in its thickness or through the advent of volume/mass sources.  
  ! note that some of the surface water work is encompassed here, hence the need
  ! to compute engint(11) and engext(11).
  if(horz_grid == MOM_BGRID) then 
     do k=1,nk
        wrk1(:,:,k) = Grd%umask(:,:,k)*REMAP_BT_TO_BU(wrk3(:,:,k))
        do j=jsc,jec
           do i=isc,iec
              boxarea     = grd_area(i,j,k)
              ke_per_mass = 0.5*(Velocity%u(i,j,k,1,tau)*Velocity%u(i,j,k,1,tau)  &
                                +Velocity%u(i,j,k,2,tau)*Velocity%u(i,j,k,2,tau))
              term        = boxarea*ke_per_mass*wrk1(i,j,k)
              compression = compression + term
           enddo
        enddo
     enddo

  else ! Cgrid 

     do k=1,nk
        wrk1(:,:,k) = Grd%tmask(:,:,k)*wrk3(:,:,k)
        do j=jsc,jec
           do i=isc,iec
              boxarea     = grd_area(i,j,k)
              ke_per_mass = 0.5*(Velocity%u(i,j,k,1,tau)*Velocity%u(i,j,k,1,tau)  &
                                +Velocity%u(i,j,k,2,tau)*Velocity%u(i,j,k,2,tau))
              term        = boxarea*ke_per_mass*wrk1(i,j,k)
              compression = compression + term
           enddo
        enddo
     enddo

  endif 

  call mpp_sum(compression)
  call mpp_sum(engint, size(engint(:)))
  call mpp_sum(engext, size(engext(:)))

  enleak = engint(2) + engint(3) - engint(11) + engext(2) + engext(3) - engext(11)- compression

  call kinetic_energy(Time, Thickness, Velocity, ke_tot, .false., .false.)
  call potential_energy(Time, Thickness, Dens, pe_tot, pe_tot_rel, .false., .false.) 

  write(stdoutunit,*) ' '
  write (stdoutunit,'(//60x,a)') &
  ' Globally integrated energy analysis over model time step '
  call write_timestamp(Time%model_time)

  if(acor > 0) then 
    write (stdoutunit,'(1x,/a)') &
    ' ==> Note: acor > 0 means Coriolis force will not conserve energy.'
  endif 
  if(Grd%tripolar) then 
    write (stdoutunit,'(/1x,/a)') &
    ' ==> NOTE: Energy conversion errors are nontrivial when using the tripolar grid. '
  endif 
  write(stdoutunit,'(/,1x,a,e20.12)')  &
   'Potential energy relative to first time (Joules) = ', pe_tot_rel
  write(stdoutunit,'(1x,a,e20.12)')    &
   'Potential energy                        (Joules) = ', pe_tot
  write(stdoutunit,'(1x,a,e20.12)')    &
   'Kinetic energy  (Joules)                         = ', ke_tot
  write (stdoutunit,9100) &
   'Globally integrated U dot momentum eqns (Watt) ', ucell_mass, Grd%ucella(1)
  write (stdoutunit,9101) ' time rate of change    ', engint(1)+engext(1), engint(1), engext(1)
  write (stdoutunit,9101) ' horizontal advection   ', engint(2)+engext(2), engint(2), engext(2)
  write (stdoutunit,9101) ' vertical advection     ', engint(3)+engext(3), engint(3), engext(3)
  write (stdoutunit,9101) ' pressure force         ', press_int+press_ext, press_int, press_ext
  write (stdoutunit,9101) ' laplacian friction     ', engint(4)+engext(4), engint(4), engext(4)
  write (stdoutunit,9101) ' biharmonic friction    ', engint(5)+engext(5), engint(5), engext(5)
  write (stdoutunit,9101) ' vertical friction      ', engint(6)+engext(6), engint(6), engext(6)
  write (stdoutunit,9101) ' work by wind           ', engint(7)+engext(7), engint(7), engext(7)
  write (stdoutunit,9101) ' work by bottom drag    ', engint(8)+engext(8), engint(8), engext(8)
  write (stdoutunit,9101) ' coriolis forces        ', engint(9)+engext(9), engint(9), engext(9)
  write (stdoutunit,9101) ' work by water flux     ', engint(10)+engext(10), engint(10), engext(10)
  write (stdoutunit,9101) ' work by momentum src   ', engint(12)+engext(12), engint(12), engext(12)
  write (stdoutunit,9101) ' work by Aiki drag      ', engint(13)+engext(13), engint(13), engext(13)

  if(horz_grid==MOM_BGRID) then
     write (stdoutunit,9110) press_conv, press_int+press_ext, press_conv_err, compression, enleak
  else
     write (stdoutunit,9111) press_conv, press_int+press_ext, press_conv_err
  endif

9100 format(/1x,a,1x,/,' ocean mass         =',e16.9,' kg',/, ' ocean surface area =',e16.9,&
          ' m^2',//,35x,'total(Watt)             internal(Watt)             external(Watt)')
9101 format(a25,3(3x,es24.17))
9110 format(/ ' power contributed by pressure conversion and mass tendencies           = ',es24.17,&
            /,' power from horizontal pressure force -u dot force(p)                   = ',es24.17,&
            /,' power mismatch between -u dot grad(p) and its conversions              = ',es24.17,&
            /,' power from grid cell mass changes due to compresibility and/or sources = ',es24.17,&
            /,' power mismatch between -u dot grad (v u) and its conversions           = ',es24.17/)
9111 format(/ ' power contributed by pressure conversion and mass tendencies           = ',es24.17,&
            /,' power from horizontal pressure force -u dot force(p)                   = ',es24.17,&
            /,' power mismatch between -u dot grad(p) and its conversions              = ',es24.17)
 

  call mpp_clock_end(id_clock_energy_analysis)

end subroutine energy_analysis
! </SUBROUTINE> NAME="energy_analysis"


!#######################################################################
! <SUBROUTINE NAME="cfl_check1_bgrid">
!
! <DESCRIPTION>
! Perform the first of two CFL checks on horizontal velocity. 
!
! Assume Bgrid here.  Cgrid calculation not affected much
! for purposes of the CFL check.
!
! Vectorized version from Russell.Fiedler@csiro.au computes cfl
! values at a single latitude. The location of the maximum at this
! latitude is calculated via the maxloc() intrinsic. The maximum 
! value for this processor is then updated if necessary.
! 
! </DESCRIPTION>
subroutine cfl_check1_bgrid(Time, Thickness, Velocity)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type),  intent(in) :: Velocity

  real     :: cflup,cflvp        ! percent of cfl criteria reached by velocity component
  real     :: cflum,cflvm        ! velocity component which comes closest to its cfl criteria
  integer  :: icflu,jcflu,kcflu  ! "i","j","k" coordinate of "cflum"
  integer  :: icflv,jcflv,kcflv  ! "i","j","k" coordinate of "cflvm"
  real     :: dtmax
  real     :: cflup0, cflvp0
  real     :: fudge
  integer  :: i, j, k, tau

  ! work array containing cfl values 
  real , dimension(isc:iec) :: tempcfl

  ! array required for maxloc() intrinsic function
  integer, dimension(1) :: itemp

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(horz_grid == MOM_CGRID) return 

  tau = Time%tau   

  write (stdoutunit,'(/60x,a/)') ' Horizontal CFL summary I:  '
  write (stdoutunit,'(1x,a/)')'Location of largest horizontal CFL.'

  ! to distinguish processors 
  fudge = 1 + 1.e-12*mpp_pe() 

  icflu=isc; jcflu=jsc; kcflu=1; cflup=epsln; cflum=0.0
  icflv=isc; jcflv=jsc; kcflv=1; cflvp=epsln; cflvm=0.0
  dtmax = dtime  

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tempcfl(i)=abs(Velocity%u(i,j,k,1,tau)*Grd%dxur(i,j))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflup) then
            cflup=tempcfl(itemp(1))
            icflu=itemp(1)
            jcflu=j
            kcflu=k
        endif
        do i=isc,iec
           tempcfl(i)=abs(Velocity%u(i,j,k,2,tau)*Grd%dyur(i,j))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflvp) then
            cflvp=tempcfl(itemp(1))
            icflv=itemp(1)
            jcflv=j
            kcflv=k
        endif
     enddo
  enddo

  ! multiply by 100.0 to convert to percentages 
  ! multiply by dtmax to convert to dimensionless CFL number 
  ! multipby cflwtp by rho0r for dimensional reasons as well
  cflup  = 100.0*cflup*dtmax
  cflvp  = 100.0*cflvp*dtmax

  cflum  = Velocity%u(icflu,jcflu,kcflu,1,tau)
  cflvm  = Velocity%u(icflv,jcflv,kcflv,2,tau)

  cflup   = cflup*fudge
  cflvp   = cflvp*fudge
  cflup0  = cflup
  cflvp0  = cflvp
  call mpp_max(cflup)
  call mpp_max(cflvp)

  if (cflup == cflup0) then

    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f12.3,a,f12.3,a,f10.3,a)')&
   ' u (',cflum,' m/s) is ',cflup/fudge,' % of CFL (',abs(100.0*cflum/(cflup/fudge)), &
   ' m/s) at (i,j,k) = (',icflu+Dom%ioff,',',jcflu+Dom%joff,',',kcflu,'),',           &
   ' (lon,lat,thk) = (',Grd%xu(icflu,jcflu),',',Grd%yu(icflu,jcflu),',',              &
   Thickness%depth_zu(icflu,jcflu,kcflu),' m)'

    write (unit,'(a,e12.3,a,e12.3,a,f10.3)')&
    ' where grid is (m) (dxu,dyu,dzu) = (', &
    Grd%dxu(icflu,jcflu),',',Grd%dyu(icflu,jcflu),',',Thickness%dzu(icflu,jcflu,kcflu)

  endif

  if (cflvp == cflvp0) then

    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f12.3,a,f12.3,a,f10.3,a)') &
    ' v (',cflvm,' m/s) is ',cflvp/fudge,' % of CFL (',abs(100.0*cflvm/(cflvp/fudge)), &
    ' m/s) at (i,j,k) = (',icflv+Dom%ioff,',',jcflv+Dom%joff,',',kcflv,'),',           &
    ' (lon,lat,thk) = (',Grd%xu(icflv,jcflv),',',Grd%yu(icflv,jcflv),',',              &
    Thickness%depth_zu(icflv,jcflv,kcflv),' m)'

    write (unit,'(a,e12.3,a,e12.3,a,f10.3)')          &
   ' where grid is (m) (dxu,dyu,dzu) = (',            &
    Grd%dxu(icflv,jcflv),',',Grd%dyu(icflv,jcflv),',',&
    Thickness%dzu(icflv,jcflv,kcflv)

  endif

end subroutine cfl_check1_bgrid
! </SUBROUTINE> NAME="cfl_check1_bgrid"


!#######################################################################
! <SUBROUTINE NAME="cfl_check1_cgrid">
!
! <DESCRIPTION>
! Perform the first of two CFL checks on horizontal velocity. 
!
! Assume Cgrid here.
!
! Vectorized version from Russell.Fiedler@csiro.au computes cfl
! values at a single latitude. The location of the maximum at this
! latitude is calculated via the maxloc() intrinsic. The maximum 
! value for this processor is then updated if necessary.
! 
! </DESCRIPTION>
subroutine cfl_check1_cgrid(Time, Thickness, Velocity)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type),  intent(in) :: Velocity

  real     :: cflup,cflvp        ! percent of cfl criteria reached by velocity component
  real     :: cflum,cflvm        ! velocity component which comes closest to its cfl criteria
  integer  :: icflu,jcflu,kcflu  ! "i","j","k" coordinate of "cflum"
  integer  :: icflv,jcflv,kcflv  ! "i","j","k" coordinate of "cflvm"
  real     :: dtmax
  real     :: cflup0, cflvp0
  real     :: fudge
  integer  :: i, j, k, tau

  ! work array containing cfl values 
  real , dimension(isc:iec) :: tempcfl

  ! array required for maxloc() intrinsic function
  integer, dimension(1) :: itemp

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(horz_grid == MOM_BGRID) return 

  tau = Time%tau   

  write (stdoutunit,'(/60x,a/)') ' Horizontal CFL summary I:  '
  write (stdoutunit,'(1x,a/)')'Location of largest horizontal CFL.'

  ! to distinguish processors 
  fudge = 1 + 1.e-12*mpp_pe() 

  icflu=isc; jcflu=jsc; kcflu=1; cflup=epsln; cflum=0.0
  icflv=isc; jcflv=jsc; kcflv=1; cflvp=epsln; cflvm=0.0
  dtmax = dtime  

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tempcfl(i)=abs(Velocity%u(i,j,k,1,tau)*Grd%dxter(i,j))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflup) then
            cflup=tempcfl(itemp(1))
            icflu=itemp(1)
            jcflu=j
            kcflu=k
        endif
        do i=isc,iec
           tempcfl(i)=abs(Velocity%u(i,j,k,2,tau)*Grd%dytnr(i,j))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflvp) then
            cflvp=tempcfl(itemp(1))
            icflv=itemp(1)
            jcflv=j
            kcflv=k
        endif
     enddo
  enddo

  ! multiply by 100.0 to convert to percentages 
  ! multiply by dtmax to convert to dimensionless CFL number 
  ! multipby cflwtp by rho0r for dimensional reasons as well
  cflup  = 100.0*cflup*dtmax
  cflvp  = 100.0*cflvp*dtmax

  cflum  = Velocity%u(icflu,jcflu,kcflu,1,tau)
  cflvm  = Velocity%u(icflv,jcflv,kcflv,2,tau)

  cflup   = cflup*fudge
  cflvp   = cflvp*fudge
  cflup0  = cflup
  cflvp0  = cflvp
  call mpp_max(cflup)
  call mpp_max(cflvp)

  if (cflup == cflup0) then

    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f12.3,a,f12.3,a,f10.3,a)')&
   ' u (',cflum,' m/s) is ',cflup/fudge,' % of CFL (',abs(100.0*cflum/(cflup/fudge)), &
   ' m/s) at (i,j,k) = (',icflu+Dom%ioff,',',jcflu+Dom%joff,',',kcflu,'),',           &
   ' (lon,lat,thk) = (',Grd%xt(icflu,jcflu),',',Grd%yt(icflu,jcflu),',',              &
   Thickness%depth_zt(icflu,jcflu,kcflu),' m)'

    write (unit,'(a,e12.3,a,e12.3,a,f10.3)')&
    ' where grid is (m) (dxt,dyt,dzt) = (', &
    Grd%dxt(icflu,jcflu),',',Grd%dyt(icflu,jcflu),',',Thickness%dzt(icflu,jcflu,kcflu)

  endif

  if (cflvp == cflvp0) then

    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f12.3,a,f12.3,a,f10.3,a)') &
    ' v (',cflvm,' m/s) is ',cflvp/fudge,' % of CFL (',abs(100.0*cflvm/(cflvp/fudge)), &
    ' m/s) at (i,j,k) = (',icflv+Dom%ioff,',',jcflv+Dom%joff,',',kcflv,'),',           &
    ' (lon,lat,thk) = (',Grd%xt(icflv,jcflv),',',Grd%yt(icflv,jcflv),',',              &
    Thickness%depth_zt(icflv,jcflv,kcflv),' m)'

    write (unit,'(a,e12.3,a,e12.3,a,f10.3)')          &
   ' where grid is (m) (dxu,dyu,dzu) = (',            &
    Grd%dxt(icflv,jcflv),',',Grd%dyt(icflv,jcflv),',',&
    Thickness%dzt(icflv,jcflv,kcflv)

  endif

end subroutine cfl_check1_cgrid
! </SUBROUTINE> NAME="cfl_check1_cgrid"


!#######################################################################
! <SUBROUTINE NAME="cfl_check2_bgrid">
!
! <DESCRIPTION>
! Perform the second of two CFL checks on horizontal velocity.  
!
! Assume Bgrid here. 
!
! Bring the model down if too many large Courant numbers detected.  
! </DESCRIPTION>
subroutine cfl_check2_bgrid(Time, Thickness, Velocity)

 type(ocean_time_type),      intent(in) :: Time
 type(ocean_thickness_type), intent(in) :: Thickness
 type(ocean_velocity_type),  intent(in) :: Velocity

  real :: dtmax, cflu, cflv
  real :: umax, pcflu, vmax, pcflv, scl
  real, dimension(isd:ied,jsd:jed) :: u, v
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (horz_grid == MOM_CGRID) return 
  tau = Time%tau   

  write (stdoutunit,'(/60x,a/)') ' Horizontal CFL summary II:  '
  write (stdoutunit,'(1x,a,f5.2/)') &
  'Locations (if any) where horizontal Courant number exceeds ', large_cfl_value

  dtmax = dtime  
  do k=1,nk
    u(:,:) = Velocity%u(:,:,k,1,tau)
    v(:,:) = Velocity%u(:,:,k,2,tau)
    do j=jsc,jec
      do i=isc,iec
        cflu  = abs((dtmax/Grd%dxu(i,j))*u(i,j))    
        cflv  = abs((dtmax/Grd%dyu(i,j))*v(i,j))

        if (cflu >= large_cfl_value .or. cflv >= large_cfl_value) then

          write (unit,'(/,a,i4,a1,i3,a,i3,a,f6.3)')      &
           ' Note: CFL horizontal velocity limit exceeded at (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',large_cfl_value 

          umax  = Grd%dxu(i,j)/dtmax            
          pcflu = abs(100.0*u(i,j)/umax)
          vmax  = Grd%dyu(i,j)/dtmax
          pcflv = abs(100.0*v(i,j)/vmax)
          write (unit,'(a,f8.2,a,g15.8,a)') &
          ' u reached   ', pcflu,' % of the local CFL limit (',umax,' m/s)'
          write (unit,'(a,f8.2,a,g15.8,a)') & 
          ' v reached   ', pcflv,' % of the local CFL limit (',vmax,' m/s)'
          write (unit,'(a,e12.3,a,e12.3,a,f10.3)')                         &
          ' where the grid cell has dimensions (metres) (dxu,dyu,dzu) = (',&
          Grd%dxu(i,j),',',Grd%dyu(i,j),',',Thickness%dzu(i,j,k)
        endif 

        if (cflu >= max_cfl_value .or. cflv >= max_cfl_value) then

           write (unit,'(/,a,i4,a1,i3,a,i3,a,f6.3)')     &
           ' Note: CFL horizontal velocity limit exceeded at (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',max_cfl_value 

           write(*,'(/a)') &
           '==>Error in ocean_velocity_diag: max_cfl_value exceeded for horizontal velocity.'
           write(*,'(a)')  &
           '   The model is bringing itself down in order for the user to determine '
           write(*,'(a)')  &
           '   whether the simulation is going to produce Inf or NaN, or whether it '
           write(*,'(a)')  &
           '   will successfully run through the present phase with large Courant numbers.'
           write(*,'(a)')  &
           '   It is often the case that early spin-ups produce very large Courant'
           write(*,'(a)')  &
           '   numbers, with the model later adjusting to stable behaviour with modest'
           write(*,'(a)')  &
           '   Courant numbers. To test if the model will reach stable behaviour,'
           write(*,'(a)')  &
           '   try significally increasing "max_cfl_value" in ocean_velocity_diag_nml.'

          if(verbose_cfl) then 

              write(*,'(/a,f6.3/)')                                                      & 
              ' Information about region where horizontal Courant number is larger than',&
              max_cfl_value
              is = max(isc,i-3)
              ie = min(iec,i+3)
              js = max(jsc,j-3)
              je = min(jec,j+3)
              scl = 0.0
              write (*,9100)'u (m/s)', Time%itt, k, Thickness%depth_zu(is,js,k), &
                   Grd%xu(is,j), Grd%xu(ie,j), Grd%yu(i,js), Grd%yu(i,je), scl
              call matrix (u(is:ie,js:je), is, ie, js, je, scl)

              write (*,9100) 'v (m/s)', Time%itt, k, Thickness%depth_zu(is,js,k), &
                   Grd%xu(is,j), Grd%xu(ie,j), Grd%yu(i,js), Grd%yu(i,je), scl
              call matrix (v(is:ie,js:je), is, ie, js, je, scl)

          endif

          call mpp_error(FATAL,&
          '==>Error in ocean_velocity_diag_mod: max_cfl_value exceeded.')

        endif

      enddo
    enddo
  enddo

9100 format(1x,a12,1x,'ts=',i10,1x,',k=',i3,', z(m)=',f8.1,', lon(deg):',f6.2,' --> ',f6.2,&
            ', lat(deg):',f6.2,' --> ',f6.2,', scaling=',1pg10.3)

end subroutine cfl_check2_bgrid
! </SUBROUTINE> NAME="cfl_check2_bgrid"


!#######################################################################
! <SUBROUTINE NAME="cfl_check2_cgrid">
!
! <DESCRIPTION>
! Perform the second of two CFL checks on horizontal velocity.  
!
! Assume Cgrid here. 
!
! Bring the model down if too many large Courant numbers detected.  
! </DESCRIPTION>
subroutine cfl_check2_cgrid(Time, Thickness, Velocity)

 type(ocean_time_type),      intent(in) :: Time
 type(ocean_thickness_type), intent(in) :: Thickness
 type(ocean_velocity_type),  intent(in) :: Velocity

  real :: dtmax, cflu, cflv
  real :: umax, pcflu, vmax, pcflv, scl
  real, dimension(isd:ied,jsd:jed) :: u, v
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (horz_grid == MOM_BGRID) return 
  tau = Time%tau   

  write (stdoutunit,'(/60x,a/)') ' Horizontal CFL summary II:  '
  write (stdoutunit,'(1x,a,f5.2/)') &
  'Locations (if any) where horizontal Courant number exceeds ', large_cfl_value

  dtmax = dtime  
  do k=1,nk
    u(:,:) = Velocity%u(:,:,k,1,tau)
    v(:,:) = Velocity%u(:,:,k,2,tau)
    do j=jsc,jec
      do i=isc,iec
        cflu  = abs((dtmax*Grd%dxter(i,j))*u(i,j))    
        cflv  = abs((dtmax*Grd%dytnr(i,j))*v(i,j))

        if (cflu >= large_cfl_value .or. cflv >= large_cfl_value) then

          write (unit,'(/,a,i4,a1,i3,a,i3,a,f6.3)')      &
           ' Note: CFL horizontal velocity limit exceeded at (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',large_cfl_value 

          umax  = Grd%dxt(i,j)/dtmax            
          pcflu = abs(100.0*u(i,j)/umax)
          vmax  = Grd%dyt(i,j)/dtmax
          pcflv = abs(100.0*v(i,j)/vmax)
          write (unit,'(a,f8.2,a,g15.8,a)') &
          ' u reached   ', pcflu,' % of the local CFL limit (',umax,' m/s)'
          write (unit,'(a,f8.2,a,g15.8,a)') & 
          ' v reached   ', pcflv,' % of the local CFL limit (',vmax,' m/s)'
          write (unit,'(a,e12.3,a,e12.3,a,f10.3)')                         &
          ' where the grid cell has dimensions (metres) (dxt,dyt,dzt) = (',&
          Grd%dxt(i,j),',',Grd%dyt(i,j),',',Thickness%dzt(i,j,k)
        endif 

        if (cflu >= max_cfl_value .or. cflv >= max_cfl_value) then

           write (unit,'(/,a,i4,a1,i3,a,i3,a,f6.3)')     &
           ' Note: CFL horizontal velocity limit exceeded at (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',max_cfl_value 

           write(*,'(/a)') &
           '==>Error in ocean_velocity_diag: max_cfl_value exceeded for horizontal velocity.'
           write(*,'(a)')  &
           '   The model is bringing itself down in order for the user to determine '
           write(*,'(a)')  &
           '   whether the simulation is going to produce Inf or NaN, or whether it '
           write(*,'(a)')  &
           '   will successfully run through the present phase with large Courant numbers.'
           write(*,'(a)')  &
           '   It is often the case that early spin-ups produce very large Courant'
           write(*,'(a)')  &
           '   numbers, with the model later adjusting to stable behaviour with modest'
           write(*,'(a)')  &
           '   Courant numbers. To test if the model will reach stable behaviour,'
           write(*,'(a)')  &
           '   try significally increasing "max_cfl_value" in ocean_velocity_diag_nml.'

          if(verbose_cfl) then 

              write(*,'(/a,f6.3/)')                                                      & 
              ' Information about region where horizontal Courant number is larger than',&
              max_cfl_value
              is = max(isc,i-3)
              ie = min(iec,i+3)
              js = max(jsc,j-3)
              je = min(jec,j+3)
              scl = 0.0
              write (*,9100)'u (m/s)', Time%itt, k, Thickness%depth_zt(is,js,k), &
                   Grd%xt(is,j), Grd%xt(ie,j), Grd%yt(i,js), Grd%yt(i,je), scl
              call matrix (u(is:ie,js:je), is, ie, js, je, scl)

              write (*,9100) 'v (m/s)', Time%itt, k, Thickness%depth_zt(is,js,k), &
                   Grd%xt(is,j), Grd%xt(ie,j), Grd%yt(i,js), Grd%yt(i,je), scl
              call matrix (v(is:ie,js:je), is, ie, js, je, scl)

          endif

          call mpp_error(FATAL,&
          '==>Error in ocean_velocity_diag_mod: max_cfl_value exceeded.')

        endif

      enddo
    enddo
  enddo

9100 format(1x,a12,1x,'ts=',i10,1x,',k=',i3,', z(m)=',f8.1,', lon(deg):',f6.2,' --> ',f6.2,&
            ', lat(deg):',f6.2,' --> ',f6.2,', scaling=',1pg10.3)


end subroutine cfl_check2_cgrid
! </SUBROUTINE> NAME="cfl_check2_cgrid"


!#######################################################################
! <SUBROUTINE NAME="stokes_coriolis_force">
!
! <DESCRIPTION>
!
! Diagnostic to compute thickness weighted and density
! weighted acceleration from the Stokes coriolis force, where the
! Stokes drift arises from surface ocean waves.  To obtain the
! Stokes drift requires coupling MOM to a surface wave model.
!
! Note: for a hydrostatic model, we should NOT
! include the Stokes-Coriolis force as part of the
! model prognostic equations. It is computed here only
! for diagnostic purposes.
!
! Assume stokes_drift is on U-grid point
! (Stephen.Griffies: to be revisited when couple to wave model).
!
! </DESCRIPTION>
!
subroutine stokes_coriolis_force(Time, Thickness, Velocity)

  type(ocean_time_type),       intent(in)    :: Time
  type(ocean_thickness_type),  intent(in)    :: Thickness
  type(ocean_velocity_type),   intent(inout) :: Velocity

  integer :: i, j, k
  integer :: tau
  tau = Time%tau
  wrk1_v(:,:,:,:) = 0.0

  if(id_stokes_force_x > 0 .or. id_stokes_force_y > 0) then
     if(horz_grid == MOM_BGRID) then
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1_v(i,j,k,1) =  Grd%f(i,j)*Velocity%stokes_drift(i,j,k,2)*Thickness%rho_dzu(i,j,k,tau)
                 wrk1_v(i,j,k,2) = -Grd%f(i,j)*Velocity%stokes_drift(i,j,k,1)*Thickness%rho_dzu(i,j,k,tau)
              enddo
           enddo
        enddo
     else   ! Cgrid
        call mpp_update_domains(Velocity%stokes_drift(:,:,:,1),Velocity%stokes_drift(:,:,:,2),Dom%domain2d,gridtype=BGRID_NE)
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1_v(i,j,k,1) =  onehalf &
                       *( Grd%f(i,j)*Velocity%stokes_drift(i,j,k,2)*Thickness%rho_dzu(i,j,k,tau) &
                         +Grd%f(i,j-1)*Velocity%stokes_drift(i,j-1,k,2)*Thickness%rho_dzu(i,j-1,k,tau))
                 wrk1_v(i,j,k,2) =  -onehalf &
                       *( Grd%f(i,j)*Velocity%stokes_drift(i,j,k,1)*Thickness%rho_dzu(i,j,k,tau) &
                         +Grd%f(i-1,j)*Velocity%stokes_drift(i-1,j,k,1)*Thickness%rho_dzu(i-1,j,k,tau))
              enddo
           enddo
        enddo
     endif
  endif

  call diagnose_3d_en(Time, Grd, id_stokes_force_x, id_stokes_force_y, wrk1_v(:,:,:,:))

end subroutine stokes_coriolis_force
! </SUBROUTINE> NAME="stokes_coriolis_force"


end module ocean_velocity_diag_mod


