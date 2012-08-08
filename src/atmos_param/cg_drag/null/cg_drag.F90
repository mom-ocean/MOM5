                     module cg_drag_mod

use fms_mod,                only:  fms_init, mpp_pe, mpp_root_pe,  &
                                   file_exist, check_nml_error,  &
                                   error_mesg,  FATAL, WARNING, NOTE, &
                                   close_file, open_namelist_file, &
                                   stdlog, write_version_number, &
                                   open_restart_file
use time_manager_mod,       only:  time_manager_init, time_type
use diag_manager_mod,       only:  diag_manager_init,   &
                                   register_diag_field, send_data
use constants_mod,          only:  constants_init, PI, RDGAS, GRAV, CP_AIR
use column_diagnostics_mod, only:  column_diagnostics_init, &
                                   initialize_diagnostic_columns, &
                                   column_diagnostics_header, &
                                   close_column_diagnostics_units

!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    cg_drag_mod computes the convective gravity wave forcing on 
!    the zonal flow. the parameterization is described in Alexander and 
!    Dunkerton [JAS, 15 December 1999]. 
!--------------------------------------------------------------------
  

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: cg_drag.F90,v 18.0 2010/03/02 23:28:42 fms Exp $'
character(len=128)  :: tagname =  '$Name: siena_201207 $'



!---------------------------------------------------------------------
!-------  interfaces --------

public    cg_drag_init, cg_drag_calc, cg_drag_end, cg_drag_restart, &
          cg_drag_time_vary, cg_drag_endts


logical          :: module_is_initialized=.false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------



                        contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                      PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################

subroutine cg_drag_init (lonb, latb, pref, Time, axes)

!-------------------------------------------------------------------
!   cg_drag_init is the constructor for cg_drag_mod.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
real,    dimension(:,:), intent(in)    :: lonb, latb
real,    dimension(:),   intent(in)    :: pref
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb      array of model longitudes on cell boundaries [radians]
!       latb      array of model latitudes at cell boundaries [radians]
!       pref      array of reference pressures at full levels (plus
!                 surface value at nlev+1), based on 1013.25hPa pstar
!                 [ Pa ]
!       Time      current time (time_type)
!       axes      data axes for diagnostics
!
!------------------------------------------------------------------


!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)


!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('cg_drag_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cg_drag_init



!####################################################################

subroutine cg_drag_calc (is, js, lat, pfull, zfull, temp, uuu, vvv,   &
                         Time, delt, gwfcng_u, gwfcng_v)

!--------------------------------------------------------------------  
!    cg_drag_calc defines the arrays needed to calculate the convective
!    gravity wave forcing, calls gwfc to calculate the forcing, returns 
!    the desired output fields, and saves the values for later retrieval
!    if they are not calculated on every timestep.
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
integer,                intent(in)      :: is, js
real, dimension(:,:),   intent(in)      :: lat
real, dimension(:,:,:), intent(in)      :: pfull, zfull, temp, uuu, vvv
type(time_type),        intent(in)      :: Time
real           ,        intent(in)      :: delt
real, dimension(:,:,:), intent(out)     :: gwfcng_u, gwfcng_v

!-------------------------------------------------------------------
!    intent(in) variables:
!
!       is,js    starting subdomain i,j indices of data in 
!                the physics_window being integrated
!       lat      array of model latitudes at cell boundaries [radians]
!       pfull    pressure at model full levels [ Pa ]
!       zfull    height at model full levels [ m ]
!       temp     temperature at model levels [ deg K ]
!       uuu      zonal wind  [ m/s ]
!       Time     current time, needed for diagnostics [ time_type ]
!       delt     physics time step [ s ]
!
!    intent(out) variables:
!
!       gwfcng   time tendency for u eqn due to gravity-wave forcing
!                [ m/s^2 ]
!
!-------------------------------------------------------------------

      call error_mesg('cg_drag_calc', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cg_drag_calc



!###################################################################

subroutine cg_drag_end

!--------------------------------------------------------------------
!    cg_drag_end is the destructor for cg_drag_mod.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('cg_drag_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine cg_drag_end


!####################################################################
!dummy interface
subroutine cg_drag_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

end subroutine cg_drag_restart


!####################################################################
!dummy interface
subroutine  cg_drag_time_vary (delt)
real           ,        intent(in)      :: delt

end subroutine cg_drag_time_vary

!####################################################################
!dummy interface
subroutine cg_drag_endts

end subroutine cg_drag_endts

end module cg_drag_mod


