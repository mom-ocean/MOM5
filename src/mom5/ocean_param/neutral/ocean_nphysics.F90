module ocean_nphysics_mod
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<REVIEWER EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</REVIEWER>
!  
!<OVERVIEW>
!  
! Driver for ocean neutral physics.
!</OVERVIEW>
!
!<DESCRIPTION>
! Driver for ocean neutral physics.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_nphysics_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. 
!  Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  Default debug_this_module=.false.
!  </DATA> 
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="use_nphysicsA" TYPE="logical">
!  For using the nphysicsA method of neutral physics, based on that 
!  developed in MOM4.0.  This scheme is more robust and recommended for 
!  general use.  Default use_nphysicsA=.true. 
!  </DATA> 
!  <DATA NAME="use_nphysicsB" TYPE="logical">
!  For using the nphysicsB method of neutral physics.  This method is 
!  experimental, and is not recommended for general use.
!  Default use_nphysicsB=.false. 
!  </DATA> 
!  <DATA NAME="use_nphysicsC" TYPE="logical">
!  For using the nphysicsC method of neutral physics.  This method is 
!  experimental, and is not recommended for general use.  
!  Default use_nphysicsC=.false. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln, pi
use diag_manager_mod, only: register_diag_field
use fms_mod,          only: FATAL, WARNING, NOTE
use fms_mod,          only: file_exist
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_mod,          only: input_nml_file, mpp_error, stdout, stdlog
use mpp_domains_mod,  only: mpp_update_domains

use ocean_domains_mod,       only: get_local_indices
use ocean_nphysics_util_mod, only: ocean_nphysics_util_init
use ocean_nphysicsA_mod,     only: ocean_nphysicsA_init, ocean_nphysicsA_end, nphysicsA
use ocean_nphysicsA_mod,     only: ocean_nphysicsA_restart
use ocean_nphysicsB_mod,     only: ocean_nphysicsB_init, ocean_nphysicsB_end, nphysicsB
use ocean_nphysicsB_mod,     only: ocean_nphysicsB_restart
use ocean_nphysicsC_mod,     only: ocean_nphysicsC_init, ocean_nphysicsC_end, nphysicsC
use ocean_nphysicsC_mod,     only: ocean_nphysicsC_restart
use ocean_parameters_mod,    only: TERRAIN_FOLLOWING, missing_value, onefourth, grav, rho0r
use ocean_types_mod,         only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_thickness_type, ocean_density_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_util_mod,          only: write_note, write_warning


implicit none

public ocean_nphysics_init
public ocean_nphysics_end
public neutral_physics 
public ocean_nphysics_restart

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

logical :: used

logical :: use_nphysicsA     = .true.   
logical :: use_nphysicsB     = .false.   
logical :: use_nphysicsC     = .false.   
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.
logical :: write_a_restart   = .true. 

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics.F90,v 20.0 2013/12/14 00:14:34 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__

logical :: module_is_initialized = .FALSE.

namelist /ocean_nphysics_nml/ use_this_module, debug_this_module, write_a_restart, &
                              use_nphysicsA, use_nphysicsB, use_nphysicsC

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module. 
! </DESCRIPTION>
!
subroutine ocean_nphysics_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,    &
                               Ocean_options, vert_coordinate_type, vert_coordinate_class, &
                               cmip_units, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  integer,                      intent(in)           :: vert_coordinate_type
  integer,                      intent(in)           :: vert_coordinate_class
  logical,                      intent(in)           :: cmip_units
  logical,                      intent(in), optional :: debug

  integer :: ioun, io_status, ierr  
  integer :: num_schemes 
  real    :: agm_closure_lower_depth
  real    :: agm_closure_upper_depth
  real    :: agm_closure_buoy_freq
  real    :: smax 
  real    :: swidth 
 
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  
  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysics_mod (ocean_nphysics_init):already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  Dom => Domain
  Grd => Grid

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_nphysics_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_nphysics_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysics_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_nphysics_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysics_nml)  
  write (stdlogunit,ocean_nphysics_nml)

  if(use_this_module) then 
     call write_note(FILENAME,&
     'USING ocean_nphysics_mod.')
    if(vert_coordinate_type==TERRAIN_FOLLOWING) then 
       call mpp_error(WARNING, &
            '==>Warning: ocean_nphysics is NOT supported with TERRRAIN_FOLLOWING vertical coordinates.')
    endif
  else 
     call write_note(FILENAME,&
     'NOT using ocean_nphysics_mod.')
     Ocean_options%neutral_physics = 'Did NOT use neutral physics option.'
     return
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
     call write_note(FILENAME,&
     'running ocean_nphysics_mod with debug_this_module=.true.')
  endif

  if(.not. write_a_restart) then
     call write_warning(FILENAME,&
     'NO restart written.')
  endif 

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  call ocean_nphysics_util_init(Grid, Domain, Time, Time_steps, Dens, T_prog,  &
       agm_closure_lower_depth, agm_closure_upper_depth, agm_closure_buoy_freq,&
       smax, swidth, cmip_units, debug)

  num_schemes=0
  if(use_nphysicsA) then 
    num_schemes = num_schemes+1
    call write_note(FILENAME,&
    'USING ocean_nphysicsA.')
    Ocean_options%neutral_physics = 'Used neutral physics with nphysicsA algorithm.'
    call ocean_nphysicsA_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,&
         vert_coordinate_class, agm_closure_lower_depth, agm_closure_upper_depth,     &
         smax, swidth, cmip_units, debug)

  elseif(use_nphysicsB) then 
    num_schemes = num_schemes+1
    call write_note(FILENAME,&
    'USING ocean_nphysicsB.')
    Ocean_options%neutral_physics = 'Used neutral physics with nphysicsB algorithm.'
    call ocean_nphysicsB_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,&
         vert_coordinate_class, agm_closure_lower_depth, agm_closure_upper_depth,     &
         smax, swidth, cmip_units, debug)

  elseif(use_nphysicsC) then 
    num_schemes = num_schemes+1
    call write_note(FILENAME,&
    'USING ocean_nphysicsC.')
    Ocean_options%neutral_physics = 'Used neutral physics with nphysicsC algorithm.'
    call ocean_nphysicsC_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,&
         vert_coordinate_class, agm_closure_lower_depth, agm_closure_upper_depth,     & 
         smax, swidth, cmip_units, debug)
  endif 

  if(num_schemes > 1) then 
     call mpp_error(FATAL,  &
    '==>ocean_nphysics_mod: Can only enable one of the nphysics schemes: A, B, or C.')
  endif
  if(num_schemes==0) then 
     call write_warning(FILENAME,&
     'no nphysics scheme enabled. Choose one of nphysicsA, nphysicsB, or nphysicsC.')
  endif


end subroutine ocean_nphysics_init
! </SUBROUTINE>  NAME="ocean_nphysics_init"


!#######################################################################
! <SUBROUTINE NAME="neutral_physics">
! <DESCRIPTION>
!
! Call the relevant neutral physics scheme.
!
! </DESCRIPTION>
!
!
subroutine neutral_physics (Time, Thickness, Dens, rho, T_prog, &
               gm_diffusivity, surf_blthick, rossby_radius_raw)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(in)    :: surf_blthick
  real, dimension(isd:,jsd:,:),   intent(inout) :: gm_diffusivity
  real, dimension(isd:,jsd:),     intent(inout) :: rossby_radius_raw

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysics (neutral_physics): needs initialization')
  endif 

  if(use_nphysicsA) then 
      call nphysicsA(Time, Thickness, Dens, rho, T_prog, &
           surf_blthick, gm_diffusivity, rossby_radius_raw)
  elseif(use_nphysicsB) then 
      call nphysicsB(Time, Thickness, Dens, rho, T_prog, &
           surf_blthick, gm_diffusivity, rossby_radius_raw)
  elseif(use_nphysicsC) then 
      call nphysicsC(Time, Thickness, Dens, rho, T_prog, &
           surf_blthick, gm_diffusivity, rossby_radius_raw)
  endif

end subroutine neutral_physics
! </SUBROUTINE> NAME="neutral_physics"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_restart">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysics_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  if(use_nphysicsA) then 
     call ocean_nphysicsA_restart(time_stamp)
  endif 
  if(use_nphysicsB) then 
     call ocean_nphysicsB_restart(time_stamp)
  endif 
  if(use_nphysicsC) then 
     call ocean_nphysicsC_restart(time_stamp)
  endif   

end subroutine ocean_nphysics_restart
! </SUBROUTINE> NAME="ocean_nphysics_restart"



!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysics_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysics (ocean_nphysics_end): needs initialization')
  endif 

  if(.not. write_a_restart) then
     call write_warning(FILENAME,&
     'NO restart written.')
     return 
  endif 

  call ocean_nphysics_restart

  if(use_nphysicsA) then 
     call ocean_nphysicsA_end(Time)
  endif 
  if(use_nphysicsB) then 
     call ocean_nphysicsB_end(Time)
  endif 
  if(use_nphysicsC) then 
     call ocean_nphysicsC_end(Time)
  endif 

end subroutine ocean_nphysics_end
! </SUBROUTINE> NAME="ocean_nphysics_end"


end module ocean_nphysics_mod
