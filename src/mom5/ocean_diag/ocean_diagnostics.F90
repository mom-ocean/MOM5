module ocean_diagnostics_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Routine that calls the various numerical diagnostics.
!</OVERVIEW>
!
!<DESCRIPTION>
! Routine that calls the various numerical diagnostics.
! </DESCRIPTION>
!
use diag_manager_mod,        only: need_data
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,                 only: FATAL, stdout, stdlog
use mpp_mod,                 only: input_nml_file, mpp_error
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use time_manager_mod,        only: time_type, increment_time

use ocean_adv_vel_diag_mod,  only: ocean_adv_vel_diag_init, ocean_adv_vel_diagnostics
use ocean_domains_mod,       only: get_local_indices
use ocean_tracer_diag_mod,   only: ocean_tracer_diag_init, ocean_tracer_diagnostics
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,         only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,         only: ocean_adv_vel_type, ocean_velocity_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,         only: ocean_external_mode_type, ocean_density_type
use ocean_types_mod,         only: ocean_thickness_type, ocean_lagrangian_type
use ocean_velocity_diag_mod, only: ocean_velocity_diag_init, ocean_velocity_diagnostics


implicit none

private

#include <ocean_memory.h>

integer :: num_prog_tracers = 0
integer :: index_temp = -1
integer :: index_salt = -1

! for diagnostics clocks 
integer :: id_adv_vel_diag
integer :: id_tracer_diag
integer :: id_velocity_diag

logical :: module_is_initialized = .FALSE.

character(len=128) :: version=&
     '$Id: ocean_diagnostics.F90,v 20.0 2013/12/14 00:12:51 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

public :: ocean_diag_init, ocean_diagnostics

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_diag module.
! </DESCRIPTION>
!
subroutine ocean_diag_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, T_diag, Dens, &
                           vert_coordinate_class, horz_grid, have_obc, cmip_units, use_blobs)

type(ocean_grid_type),        intent(in)    :: Grid
type(ocean_domain_type),      intent(in)    :: Domain
type(ocean_time_type),        intent(in)    :: Time
type(ocean_time_steps_type),  intent(in)    :: Time_steps
type(ocean_thickness_type),   intent(in)    :: Thickness
type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
type(ocean_diag_tracer_type), intent(in)    :: T_diag(:)
type(ocean_density_type),     intent(inout) :: Dens
integer,                      intent(in)    :: vert_coordinate_class 
integer,                      intent(in)    :: horz_grid 
logical,                      intent(in)    :: have_obc
logical,                      intent(in)    :: cmip_units
logical,                      intent(in)    :: use_blobs

integer :: n

integer :: stdoutunit,stdlogunit 
stdoutunit=stdout();stdlogunit=stdlog() 

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

#ifndef MOM_STATIC_ARRAYS
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

num_prog_tracers = size(T_prog, 1)

do n = 1, num_prog_tracers
   if (trim(T_prog(n)%name) == 'temp') index_temp = n
   if (trim(T_prog(n)%name) == 'salt') index_salt = n
enddo

if (index_temp < 1 .or. index_salt < 1) then 
   call mpp_error(FATAL,'==>Error in ocean_diagnostics_mod (ocean_diag_init): temp and/or salt not in tracer array')
endif 

id_adv_vel_diag  = mpp_clock_id('(Ocean diagnostics: adv_vel)'  ,grain=CLOCK_MODULE)
id_tracer_diag   = mpp_clock_id('(Ocean diagnostics: tracer)'   ,grain=CLOCK_MODULE)
id_velocity_diag = mpp_clock_id('(Ocean diagnostics: velocity)' ,grain=CLOCK_MODULE)

call ocean_adv_vel_diag_init (Grid, Domain, Time, Time_steps, T_prog, Dens, horz_grid, cmip_units)
call ocean_tracer_diag_init  (Grid, Domain, Time, Time_steps, Thickness, T_prog, T_diag, Dens, &
                              vert_coordinate_class, use_blobs, have_obc)
call ocean_velocity_diag_init(Grid, Domain, Time, Time_steps, horz_grid)


end subroutine ocean_diag_init
! </SUBROUTINE>  NAME="ocean_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_diagnostics">
!
! <DESCRIPTION>
! Call some ocean numerical diagnostics 
! </DESCRIPTION>
!
subroutine ocean_diagnostics(Time, Thickness, T_prog, T_diag, Adv_vel,&
                             Ext_mode, Dens, Velocity, &  
                             pme, melt, runoff, calving, visc_cbt, diff_cbt)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness 
  type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(in)    :: T_diag(:)
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(inout) :: Dens
  type(ocean_velocity_type),      intent(inout) :: Velocity

  real, dimension(isd:,jsd:),    intent(in) :: pme
  real, dimension(isd:,jsd:),    intent(in) :: melt
  real, dimension(isd:,jsd:),    intent(in) :: runoff
  real, dimension(isd:,jsd:),    intent(in) :: calving
  real, dimension(isd:,jsd:,:),  intent(in) :: visc_cbt
  real, dimension(isd:,jsd:,:,:),intent(in) :: diff_cbt
  
  if (size(T_prog,1) /= num_prog_tracers) then 
     call mpp_error(FATAL, '==>Error from ocean_diagnostics_mod (ocean_diagnostics): wrong size for tracer array')
  endif 

  call mpp_clock_begin(id_adv_vel_diag)
  call ocean_adv_vel_diagnostics(Time, Thickness, Adv_vel, T_prog, Dens, visc_cbt)
  call mpp_clock_end(id_adv_vel_diag)

  call mpp_clock_begin(id_tracer_diag)
  call ocean_tracer_diagnostics(Time, Thickness, T_prog, T_diag, Dens, &
                                Ext_mode, Velocity, Adv_vel, &
                                diff_cbt, pme, melt, runoff, calving)
  call mpp_clock_end(id_tracer_diag)

  call mpp_clock_begin(id_velocity_diag)
  call ocean_velocity_diagnostics(Time, Thickness, Dens, Ext_mode, Velocity)
  call mpp_clock_end(id_velocity_diag)

end subroutine ocean_diagnostics
! </SUBROUTINE>  NAME="ocean_diagnostics"


end module ocean_diagnostics_mod
