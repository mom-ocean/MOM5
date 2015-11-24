module ocean_vert_const_mod
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! Compute constant vertical viscosity and diffusivity. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes a time independent vertical viscosity and diffusivity. 
!</DESCRIPTION>
!
! <INFO>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_const_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false. 
!  </DATA> 
!
!  <DATA NAME="kappa_h" UNITS="m^2/sec" TYPE="real">
!  The constant vertical diffusivity.  Used for cases when wanting a space-time
!  independent diffusivity.  The "h" is historical and stands for "heat".
!  </DATA> 
!  <DATA NAME="kappa_m" UNITS="m^2/sec" TYPE="real">
!  The constant vertical viscosity.  Used for cases when wanting a space-time
!  independent viscosity.  
!  </DATA> 
!  <DATA NAME="diff_cbt_limit" UNITS="m^2/sec" TYPE="real">
!  The largest allowable vertical diffusivity.  Of use for cases where vertically unstable
!  columns are stabilized with a large vertical diffusivity.  
!  </DATA> 
!</NAMELIST>

use constants_mod,      only: pi
use diag_manager_mod,   only: register_diag_field
use fms_mod,            only: write_version_number, FATAL, NOTE, stdout, stdlog
use fms_mod,            only: open_namelist_file, check_nml_error, close_file
use mpp_io_mod,         only: mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII
use mpp_mod,            only: input_nml_file, mpp_error

use ocean_density_mod,    only: density_delta_z
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_density_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type
use ocean_workspace_mod,  only: wrk1
use ocean_util_mod,       only: diagnose_3d, diagnose_3d_u

implicit none

public ocean_vert_const_init
public vert_mix_const

private

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
     '=>Using: const/ocean_vert_const.F90 ($Id: ocean_vert_const.F90,v 20.0 2013/12/14 00:16:38 fms Exp $)'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: index_temp
integer :: index_salt

integer :: id_diff_cbt_const  =-1
integer :: id_visc_cbt_const  =-1
integer :: id_visc_cbu_const  =-1
integer :: id_density_delta_z =-1

integer :: num_prog_tracers
logical :: module_is_initialized = .FALSE.

logical :: use_this_module = .false.
real :: kappa_h        = 0.1e-4  ! constant vertical diffusivity (m^2/sec)
real :: kappa_m        = 1.0e-4  ! constant vertical viscosity (m^2/sec)
real :: diff_cbt_limit = 1.0     ! diffusivity used to vertically adjust (m^2/sec)

namelist /ocean_vert_const_nml/ use_this_module, kappa_h, kappa_m, diff_cbt_limit

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_const_init">
!
! <DESCRIPTION>
! Initialize the constant vertical diffusivity module.
! </DESCRIPTION>
!
subroutine ocean_vert_const_init (Grid, Domain, Time, Time_steps, T_prog)
 
  type(ocean_grid_type),              target, intent(in) :: Grid
  type(ocean_domain_type),            target, intent(in) :: Domain
  type(ocean_time_type),                      intent(in) :: Time
  type(ocean_time_steps_type),                intent(in) :: Time_steps 
  type(ocean_prog_tracer_type), dimension(:), intent(in) :: T_prog

  integer :: k, n, ioun, ierr, io_status

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_vert_const_mod (ocean_vert_const_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_vert_const_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_const_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_vert_const_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_const_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_vert_const_nml)  
  write (stdlogunit,ocean_vert_const_nml)

  if(kappa_h==0.0) then 
    write(stdoutunit,'(/a)')'==>USING constant vertical diffusivity with kappa_h = 0.0.'
  else 
    write(stdoutunit,'(/a)')'==>USING constant vertical diffusivity with kappa_h > 0.0.'
  endif
  if(kappa_m==0.0) then 
    write(stdoutunit,'(a/)')'==>USING constant vertical viscosity with kappa_m = 0.0.'
  else 
    write(stdoutunit,'(a/)')'==>USING constant vertical viscosity with kappa_m > 0.0.'
  endif

  write(stdoutunit,'(/a,f10.2)')&
   '==>Note from ocean_vert_const_mod: using forward time step for vert-frict of (secs)', &
                                  Time_steps%dtime_u 
  write(stdoutunit,'(/a,f10.2)')&
   '==>Note from ocean_vert_const_mod: using forward time step for vert-diff  of (secs)', &
                                  Time_steps%dtime_t 

  Dom => Domain
  Grd => Grid

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(use_this_module) then 
    call mpp_error(NOTE, '==>Note: USING ocean_vert_const_mod')
  else
    call mpp_error(NOTE, '==>Note: NOT using ocean_vert_const_mod')
    return 
  endif

  if(Time_steps%aidif < 1.0) then 

      do k=1,nk
         if ((Time_steps%dtime_t*kappa_h)/Grid%dzt(k)**2 >= 0.5) then
             write (stdoutunit,'(a,a,i3)')&
             '==> Warning: vertical stability criteria exceeded for vertical diff.',&
                  ' Use aidif=1.0, or a smaller "dtts" and/or "kappa_h" at level k=',k
         endif
      enddo

      do k=1,nk
         if ((Time_steps%dtime_u*kappa_m)/Grid%dzt(k)**2 >= 0.5) then
             write (stdoutunit,'(a,a,i3)')&
             '==> Warning: vertical stability criteria exceeded for vertical visc.',&
                  ' Use aidif=1.0, smaller dtuv, or smaller kappa_m at level k=',k
         endif
      enddo

  endif

  num_prog_tracers = size(T_prog)

  index_temp = -1
  index_salt = -1
  
  do n = 1, num_prog_tracers
     if (trim(T_prog(n)%name) == 'temp') index_temp = n
     if (trim(T_prog(n)%name) == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_vert_const_mod: temp or salt not present in tracer array')
  endif 

  ! register vertical diffusivity
  id_diff_cbt_const = -1
  id_diff_cbt_const = register_diag_field ('ocean_model', 'diff_cbt_const', Grid%tracer_axes(1:3), &
                      Time%model_time, 'vert diff_cbt from constant scheme',                       &
                      'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_visc_cbt_const = -1
  id_visc_cbt_const = register_diag_field ('ocean_model', 'visc_cbt_const', Grid%tracer_axes(1:3), &
                      Time%model_time, 'vert visc_cbt from constant scheme',                       &
                      'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_visc_cbu_const = -1
  id_visc_cbu_const = register_diag_field ('ocean_model', 'visc_cbu_const', Grid%vel_axes_uv(1:3), &
                      Time%model_time, 'vert visc_cbu from constant scheme',                       &
                      'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_density_delta_z = -1
  id_density_delta_z = register_diag_field ('ocean_model', 'density_delta_z', Grid%tracer_axes(1:3), &
                       Time%model_time, 'rho(k)-rho(k+1) to check stability for const diff_cbt',     &
                       'kg/m^3',missing_value=missing_value, range=(/-10.0,1e6/))

end subroutine ocean_vert_const_init
! </SUBROUTINE>  NAME="ocean_vert_const_init"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_const">
!
! <DESCRIPTION>
! This function computes the vertical diffusivity and viscosity.  
! These mixing coefficients are time independent but generally 
! arbitrary functions of space. 
! </DESCRIPTION>
!
  subroutine vert_mix_const(aidif, Time, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)

  real,                           intent(in)    :: aidif
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt

  integer :: i, j, k
  integer :: tau

  if(.not. use_this_module) then 
    return 
  endif 

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_vert_const: module must be initialized')
  endif 

  tau         = Time%tau
  wrk1(:,:,:) = 0.0

  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           visc_cbu(i,j,k)   = kappa_m*Grd%umask(i,j,k+1)
           visc_cbt(i,j,k)   = kappa_m*Grd%tmask(i,j,k+1)
           diff_cbt(i,j,k,1) = kappa_h*Grd%tmask(i,j,k+1)
        enddo
     enddo
  enddo

  ! set vertical diffusivity to diff_cbt_limit where gravitationally unstable
  if(aidif==1.0) then

      wrk1(:,:,:) = density_delta_z (Dens%rho(:,:,:,tau),                 &
                                     Dens%rho_salinity(:,:,:,tau),        &
                                     T_prog(index_temp)%field(:,:,:,tau), &
                                     Dens%pressure_at_depth(:,:,:))   

      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if (wrk1(i,j,k)*Grd%tmask(i,j,k+1) > 0.0) then
                   diff_cbt(i,j,k,1) = diff_cbt_limit
               endif
            enddo
         enddo
      enddo

  endif

  ! no distinction between temperature and salinity for this mixing scheme 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           diff_cbt(i,j,k,2) = diff_cbt(i,j,k,1)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_diff_cbt_const, diff_cbt(:,:,:,1))
  call diagnose_3d(Time, Grd, id_visc_cbt_const, visc_cbt(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_const, visc_cbu(:,:,:))
  call diagnose_3d(Time, Grd, id_density_delta_z, wrk1(:,:,:))

end subroutine vert_mix_const
! </SUBROUTINE> NAME="vert_mix_const"

end module ocean_vert_const_mod
