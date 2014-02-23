module ocean_vert_gotm_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="mschmidt@io-warnemuende.de"> Martin Schmidt 
!</CONTACT>
!
!<CONTACT EMAIL="Mike.Herzfeld@csiro.au"> Mike Herzfeld
!</CONTACT>
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> Russell Fiedler 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Vertical viscosity and diffusivity according GOTM.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains interfaces to initialize and invoke the 
! Generalized Ocean Turbulence Model (GOTM) parameterizations.
! Full documentation of the schemes available with GOTM 
! can be found at www.gotm.net.
!
! MOM is distributed with the basic routines from the 4.0
! release of GOTM.  Questions about GOTM should be directed 
! to the GOTM users group at www.gotm.net.  
!
! This module assumes a twolevel time stepping scheme is used 
! to update the turbulence scalar fields tke and diss. 
!
! Presently it has only been implemented assuming Bgrid.
! So it needs to be updated for Cgrid layout.  
!
! Here is a brief outline of the GOTM scheme:
!
! The non-conservative part of the tke 
! equation is P+B-diss, where P=shear production, 
! B=buoyancy production.  The non-conservative part of 
! the diss equation is a linear combination of tke, 
! P, B and diss. 
!
! The conservatie part of both tke and diss equations is 
! 3D advection.  
!
! So vertical shear and buoyancy contribute to the 
! source and sinks (respectively) of tke and dissipation. 
!
! The mixing coefficients are the product of a stability
! function, sqrt(tke), and turbulence length scale (the latter 
! a non-linear function of tke and diss). 
!
! In the hydro model (i.e., MOM), buoyancy fluxes are the 
! surface boundary conditions for vertical diffusion of temp (heat
! fluxes) and salt (freshwater fluxes), which in turn determine the
! density, and thus enter GOTM via the calculation of buoyancy production.
!
! Hence, buoyancy fluxes are NOT directly required as surface boundary 
! conditions for GOTM. 
!
! Wind stress provides the surface boundary condition for vertical 
! momentum mixing in MOM.  This information then enters GOTM 
! via the shear production calculation.  However, wind stress
! can enter the GOTM tke equation as the surface boundary 
! condition for vertical diffusion of tke (i.e. the prescribed, 
! or Dirichlet, condition). This information is typically used in  
! Mellor-Yamada turbulence models. An alternate no-flux (Neumann)
! condition is used in the k-e models, and so do not require wind stress.
!
! Boundary conditions for vertical diffusion of diss involve roughness,
! tke, and constants (both Neumann and Dirichlet). 
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Burchard, H., K. Bolding and M. R. Villarreal
! GOTM, a general ocean turbulence model. Theory
! implementation and test cases.
! European Communities, EUR 18745 EN, 1999
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_gotm_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is .false. 
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.  Default is .false. 
!  </DATA> 
!  <DATA NAME="do_turbulence_gotm" TYPE="logical">
!  For debugging.  If do_turbulence_gotm=.false., then
!  will not invoke the GOTM scheme. Will only advect 
!  tke and diss using 3d advection scheme. 
!  Default is .true., so that will invoke GOTM scheme.  
!  </DATA> 
!  <DATA NAME="do_advection_gotm" TYPE="logical">
!  For debugging.  If do_advection_gotm=.false., then
!  will not invoke the advection of tke and diss. 
!  Default is .true., so that will 3d advect tke and diss. 
!  </DATA> 
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="advect_gotm_method" TYPE="character">
!  For choosing how to advect the GOTM scalar fields tke and diss.
!  Options are advect_gotm_method='upwind' (the default) 
!              advect_gotm_method='sweby'
!  </DATA> 
!
!  <DATA NAME="diff_cbt_min" UNITS="m^2/sec" TYPE="real">
!  Background diffusivity.  Default is 1.0e-5. 
!  </DATA> 
!  <DATA NAME="visc_cbu_min" UNITS="m^2/sec" TYPE="real">
!  Background viscosity.  Default is 1.0e-5. 
!  </DATA> 
!  <DATA NAME="z0s" UNITS="m" TYPE="real">
!  Surface roughness length.  Default is 1m.
!  </DATA> 
!  <DATA NAME="z0b" UNITS="m" TYPE="real">
!  Bottom roughness length.  Default is .002m.
!  </DATA> 
!  <DATA NAME="min_tke" UNITS="m^2/s^2" TYPE="real">
!  Minimum turbulent kinetic energy.  Default=1.0e-6.
!  </DATA> 
!  <DATA NAME="min_diss" UNITS="m^2/s^3" TYPE="real">
!  Minimum energy dissipation.  Default=1.0e-10.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: FATAL, WARNING, stdout, stdlog
use fms_mod,          only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_mod,          only: file_exist 
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: restart_file_type, reset_field_pointer
use mpp_domains_mod,  only: mpp_update_domains, NUPDATE, EUPDATE, XUPDATE, YUPDATE
use mpp_mod,          only: input_nml_file, mpp_error

use mtridiagonal,     only: init_tridiagonal
use turbulence,       only: init_turbulence, do_turbulence, tke, eps, num, nuh
use turbulence,       only: L1d => l, cde

use ocean_domains_mod,    only: set_ocean_domain, get_local_indices
use ocean_parameters_mod, only: missing_value, onesixth, rho0r, grav 
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_adv_vel_type
use ocean_types_mod,      only: ocean_density_type, ocean_velocity_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_gotm_type, advect_gotm_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_3d_u, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_v
use ocean_obc_mod,        only: ocean_obc_mixing, ocean_obc_zero_boundary

implicit none

private

public  ocean_vert_gotm_init
public  ocean_vert_gotm_end
public  vert_mix_gotm_bgrid 
public  advect_gotm_compute
public  ocean_vert_gotm_restart

private advect_gotm_sweby
private advect_gotm_upwind 


#include <ocean_memory.h>


#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,nk)         :: visc_cbt_gotm ! vertical viscosity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk)         :: diff_cbt_gotm ! vertical diffusivity (m^2/s)
real, dimension(isd:ied,jsd:jed,nk)         :: drhodT        ! drho/dtheta     (kg/(m^3*C))
real, dimension(isd:ied,jsd:jed,nk)         :: drhodS        ! drho/dsalinity  (kg/(m^3*psu))
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_advect  ! advection mask
real, dimension(isd:ied,jsd:jed,nk)         :: flux_x        ! for x-advection
real, dimension(isd:ied,jsd:jed,nk)         :: flux_y        ! for y-advection
real, dimension(isd:ied,jsd:jed,nk)         :: flux_z        ! for z-advection
real, dimension(isd:ied,jsd:jed,nk)         :: gotm_tendency ! for advection tendency 

#else

real, dimension(:,:,:), allocatable :: visc_cbt_gotm ! vertical viscosity (m^2/s)
real, dimension(:,:,:), allocatable :: diff_cbt_gotm ! vertical diffusivity (m^2/s)
real, dimension(:,:,:), allocatable :: drhodT        ! drho/dtheta     (kg/m^3/C)
real, dimension(:,:,:), allocatable :: drhodS        ! drho/dsalinity  (kg/m^3/psu)
real, dimension(:,:,:), allocatable :: tmask_advect  ! advection mask
real, dimension(:,:,:), allocatable :: flux_x        ! for x-advection
real, dimension(:,:,:), allocatable :: flux_y        ! for y-advection
real, dimension(:,:,:), allocatable :: flux_z        ! for z-advection
real, dimension(:,:,:), allocatable :: gotm_tendency ! for advection tendency 

#endif

! number of GOTM fields that have time dimensions 
integer, parameter :: num_gotm = 2

! time step (secs) for evolving prognostic GOTM fields 
real :: dtime

! time step labels for gotm updates using two-level scheme 
integer :: tau_gotm
integer :: taup1_gotm

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), save    :: Dom_advect_gotm
type(advect_gotm_type), save, dimension(num_gotm) :: advect_gotm
type(ocean_gotm_type),  save, dimension(num_gotm) :: Gotm

integer :: index_tke
integer :: index_diss
integer :: num_prog_tracers=0
integer :: num_diag_tracers=0
integer :: index_temp
integer :: index_salt


! identification for diagnostic output
integer :: id_diff_cbt_gotm
integer :: id_visc_cbt_gotm
integer :: id_visc_cbu_gotm
integer :: id_gotm_nn
integer :: id_gotm_ss
logical :: used
integer, dimension(:), allocatable :: id_gotm_diag
integer, dimension(:), allocatable :: id_adv_flux_x
integer, dimension(:), allocatable :: id_adv_flux_y
integer, dimension(:), allocatable :: id_adv_flux_z
integer, dimension(:), allocatable :: id_adv_flux_x_int_z
integer, dimension(:), allocatable :: id_adv_flux_y_int_z
integer, dimension(:), allocatable :: id_gotm_errors

! for restart
integer  :: id_restart_tke
integer  :: id_restart_diss
integer  :: id_restart_diff
integer  :: id_restart_visc
type(restart_file_type), save :: Got_restart

character(len=256) :: version=&
     '$Id: ocean_vert_gotm.F90,v 20.0 2013/12/14 00:16:40 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: module_is_initialized = .FALSE.
integer :: advection_gotm_method = 1  ! internally set: 1=upwind, 2=sweby 
logical :: have_obc=.false.           ! for running with OBC

! nml parameters 
character(len=10) :: advect_gotm_method = 'upwind' ! options are "upwind" and "sweby" 
logical :: use_this_module    = .false.  ! logical switch for turning on the gotm scheme 
logical :: debug_this_module  = .false.  ! for debugging
logical :: do_turbulence_gotm = .true.   ! set to .false. when only wish to advect tke and diss. 
logical :: do_advection_gotm  = .true.   ! set to .false. when do not wish to advect tke and diss 
logical :: write_a_restart    = .true.   ! to write a restart 
logical :: map_velocity_gotm  = .false.  ! map velocity to t-points before calculating shear
logical :: map_production_gotm= .true.   ! calculate shear production, map to t-points and 
                                         ! calculate vertical shear at t-points
logical :: correct_adv_errors = .false.  ! set negative tke and diss to minimum value after advection
                                         ! needed for sweby advection in some cases
real :: diff_cbt_min          = 1.0e-5   ! background diffusivity (m^2/s)
real :: visc_cbu_min          = 1.0e-5   ! background viscosity (m^2/s)
real :: z0s                   = 1.0      ! surface roughness length (m)
real :: z0b                   = 0.002    ! bottom roughness length (m)
real :: min_tke               = 1.0e-6   ! minimum turbulent kinetic energy (m^2/s^2)
real :: min_diss              = 1.0e-10  ! minimum energy dissipation (m^2/s^3) 

namelist /ocean_vert_gotm_nml/ use_this_module, debug_this_module, write_a_restart,       & 
                               do_turbulence_gotm, do_advection_gotm, advect_gotm_method, &
                               diff_cbt_min, visc_cbu_min,                                &
                               z0s, z0b, min_tke, min_diss,                               &
                               map_velocity_gotm, map_production_gotm, correct_adv_errors

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_gotm_init">
!
! <DESCRIPTION>
! Initialization for the MOM wrapper to the GOTM vertical mixing scheme.
!
! For restarts: 
!
! We use twolevel time stepping scheme to update tke and diss,
! so only need to read in the taup1_gotm value. 
! call mpp_update_domains since need tke and diss in halos for 
! the advection calculation. 
!
! We need viscosity and diffusivity saved to restarts in order 
! to update the tke and diss fields within GOTM. 
!
! </DESCRIPTION>
!
subroutine ocean_vert_gotm_init (Grid, Domain, Time, Time_steps, T_prog, obc, debug)
  
  type(ocean_grid_type),     target, intent(in) :: Grid
  type(ocean_domain_type),   target, intent(in) :: Domain
  type(ocean_time_type),             intent(in) :: Time
  type(ocean_time_steps_type),       intent(in) :: Time_steps
  type(ocean_prog_tracer_type),      intent(in) :: T_prog(:)
  logical,                           intent(in) :: obc
  logical,                 optional, intent(in) :: debug

  integer :: i, j, k, m, n, kp1
  integer :: ioun, io_status, ierr
  integer :: unit=104
  logical :: opened 
  character(len=64) :: file_name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_vert_gotm_mod (ocean_vert_gotm_init): module is already initialized')
  endif 

  have_obc = obc 
  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_gotm_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_gotm_nml')
#else
  ioun =  open_namelist_file ()
  read  (ioun, ocean_vert_gotm_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_gotm_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_gotm_nml)  
  write (stdlogunit, ocean_vert_gotm_nml)

  if(use_this_module) then 
    write(stdoutunit,'(/1x,a)') '==> NOTE: USING GOTM vertical mixing scheme.'
  else 
    write(stdoutunit,'(/1x,a)')'==> NOTE: NOT using GOTM vertical mixing scheme.'
    return
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_vert_gotm_mod with debug_this_module=.true.'  
  endif 

  ! initial time step labels for gotm fields 
  tau_gotm   = 1
  taup1_gotm = 2

  ! use a two-level time scheme for tke and diss updates, 
  ! which means we update tke and diss using time step dtts.
  dtime=Time_steps%dtts 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_vert_gotm with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif

  if(.not. do_turbulence_gotm) then 
  write(stdoutunit,'(/a)') &
   '==>Note from ocean_vert_gotm_mod: NOT invoking GOTM turbulence code. Are you sure you want this?'
  endif 
  
  if (map_velocity_gotm .and. map_production_gotm) then
    call mpp_error(FATAL, &
    '==>Error: Only one of both map_velocity_gotm or map_production_gotm can be .true.')
  endif 
  
  if (map_velocity_gotm) then
  write(stdoutunit,'(/a)') &
   '==>Note from ocean_vert_gotm_mod: mapping velocity to T-points before calculating vertical shear.'
  endif 

  if (map_production_gotm) then
  write(stdoutunit,'(/a)') &
   '==>Note from ocean_vert_gotm_mod: calculating shear production at U-points before mapping to T-points.'
  endif 
  
  if(.not. do_advection_gotm) then 
  write(stdoutunit,'(/a)') &
   '==>Note from ocean_vert_gotm_mod: NOT advecting GOTM tke and diss. Are you sure you want this?'
  endif 

  if(do_advection_gotm) then 
      if(advect_gotm_method=='upwind') then 
          advection_gotm_method = 1
          write(stdoutunit,'(/a)') &
          '==>Note from ocean_vert_gotm_mod: advecting tke and diss with 3D first order upwind scheme.'
      elseif(advect_gotm_method=='sweby') then 
          advection_gotm_method = 2
          write(stdoutunit,'(/a)') &
          '==>Note from ocean_vert_gotm_mod: advecting tke and diss with 3D sweby advection scheme.'
      else
          call mpp_error(FATAL, &
          '==>Error: ocean_vert_gotm_mod: NO advection scheme chosen for GOTM fields.')
      endif
  endif

  if (Time_steps%aidif /= 1.0) then
    call mpp_error(FATAL,&
    '==>ocean_vert_gotm_init: "gotm vmix" must use aidif=1.0 for implicit vertical mixing.')
  endif

  ! temperature and salinity indicies
  index_temp=-1; index_salt=-1
  num_prog_tracers = size(T_prog(:))
  do n= 1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL, &
    '==>Error: ocean_vert_gotm_mod: temp and/or salt not present in tracer array')
  endif 

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
  allocate (visc_cbt_gotm(isd:ied,jsd:jed,nk))
  allocate (diff_cbt_gotm(isd:ied,jsd:jed,nk))
  allocate (drhodT(isd:ied,jsd:jed,nk))
  allocate (drhodS(isd:ied,jsd:jed,nk))
  allocate (tmask_advect(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate (flux_x(isd:ied,jsd:jed,nk))
  allocate (flux_y(isd:ied,jsd:jed,nk))
  allocate (flux_z(isd:ied,jsd:jed,nk))
  allocate (gotm_tendency(isd:ied,jsd:jed,nk))
  do n=1,num_gotm
     allocate(Gotm(n)%field(isd:ied,jsd:jed,nk,2)) 
     allocate(advect_gotm(n)%field(isc-2:iec+2,jsc-2:jec+2,nk))
  enddo
#endif

  ! set defaults for T-cell viscosity and diffusity 
  visc_cbt_gotm = 0.
  diff_cbt_gotm = 0.
  do k=1,nk
     kp1=min(k+1,nk)
     do j=jsd,jed
        do i=isd,ied
           visc_cbt_gotm(i,j,k) = Grd%tmask(i,j,kp1)*visc_cbu_min
           diff_cbt_gotm(i,j,k) = Grd%tmask(i,j,kp1)*diff_cbt_min
        enddo
     enddo
  enddo

  drhodT(:,:,:)        = 0.0
  drhodS(:,:,:)        = 0.0
  tmask_advect(:,:,:)  = 0.0
  flux_x(:,:,:)        = 0.0
  flux_y(:,:,:)        = 0.0
  flux_z(:,:,:)        = 0.0
  gotm_tendency(:,:,:) = 0.0

  ! diagnostic flags
  allocate (id_gotm_diag(num_gotm))
  allocate (id_adv_flux_x(num_gotm))
  allocate (id_adv_flux_y(num_gotm))
  allocate (id_adv_flux_z(num_gotm))
  allocate (id_adv_flux_x_int_z(num_gotm))
  allocate (id_adv_flux_y_int_z(num_gotm))
  allocate (id_gotm_errors(num_gotm))

  ! turbulence kinetic energy 
  index_tke = 1
  Gotm(index_tke)%name = 'tke'
  Gotm(index_tke)%longname = 'turbulence kinetic energy'
  Gotm(index_tke)%units = 'm^2/s^2'
  Gotm(index_tke)%field(:,:,:,:) = 0.0
  Gotm(index_tke)%min_value = min_tke

  ! energy dissipation
  index_diss = 2
  Gotm(index_diss)%name = 'diss'
  Gotm(index_diss)%longname = 'energy dissipation'
  Gotm(index_diss)%units = 'm^2/s^3'
  Gotm(index_diss)%field(:,:,:,:) = 0.0
  Gotm(index_diss)%min_value = min_diss

  do n=1,num_gotm
     advect_gotm(n)%field(:,:,:) = 0.0
     do m=1,2
        Gotm(n)%field(:,:,:,m) = Grd%tmask(:,:,:)*Gotm(n)%min_value
     enddo 
  enddo

  file_name = 'gotm.res.nc'
  id_restart_tke  = register_restart_field(Got_restart, file_name, 'field_tke', Gotm(index_tke)%field(:,:,:,taup1_gotm), &
               domain=Dom%domain2d)
  id_restart_diss = register_restart_field(Got_restart, file_name, 'field_diss', Gotm(index_diss)%field(:,:,:,taup1_gotm), &
               domain=Dom%domain2d)
  id_restart_diff = register_restart_field(Got_restart, file_name, 'diff_cbt_gotm', diff_cbt_gotm(:,:,:), &
               domain=Dom%domain2d)
  id_restart_visc = register_restart_field(Got_restart, file_name, 'visc_cbt_gotm', visc_cbt_gotm(:,:,:), &
               domain=Dom%domain2d)

  if(file_exist('INPUT/gotm.res.nc')) then
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'Reading GOTM restart for tke and diss'
      call restore_state(Got_restart)
      call mpp_update_domains(Gotm(index_tke)%field(:,:,:,taup1_gotm), Dom%domain2d)
      call mpp_update_domains(Gotm(index_diss)%field(:,:,:,taup1_gotm), Dom%domain2d)
      call mpp_update_domains(diff_cbt_gotm(:,:,:), Dom%domain2d)
      call mpp_update_domains(visc_cbt_gotm(:,:,:), Dom%domain2d)

      write(stdoutunit,*) 'From ocean_vert_gotm_mod: initial chksum ==>'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('start field_tke', Gotm(index_tke)%field(COMP,:,taup1_gotm)*Grd%tmask(COMP,:))
      call write_chksum_3d('start field_diss', Gotm(index_diss)%field(COMP,:,taup1_gotm)*Grd%tmask(COMP,:))
      call write_chksum_3d('start diff_cbt_gotm', diff_cbt_gotm(COMP,:)*Grd%tmask(COMP,:))
      call write_chksum_3d('start visc_cbt_gotm', visc_cbt_gotm(COMP,:)*Grd%tmask(COMP,:))
  endif


  ! read the GOTM input parameters 
  do 
     inquire(unit=unit, opened=opened)
     if(.NOT.opened)exit
     unit=unit+1
     if(unit ==500) then 
        call mpp_error( FATAL,'==>ocean_vert_gotm_mod: Unable to locate unit number for gotmturb.inp.')
     endif 
  enddo 
  call init_turbulence(unit,'INPUT/gotmturb.inp',nk)
  call init_tridiagonal(nk)

  ! initialize advection of scalar turbulence fields: tke and diss
  call set_ocean_domain(Dom_advect_gotm,Grd,xhalo=2,yhalo=2,name='advect_gotm', maskmap=Dom%maskmap)
  tmask_advect  = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmask_advect(i,j,k) = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_advect,Dom_advect_gotm%domain2d)  


  ! register diagnostic output 
  id_gotm_diag        = -1
  id_adv_flux_x       = -1
  id_adv_flux_y       = -1
  id_adv_flux_z       = -1
  id_adv_flux_x_int_z = -1
  id_adv_flux_y_int_z = -1
  id_diff_cbt_gotm    = -1
  id_visc_cbt_gotm    = -1
  id_visc_cbu_gotm    = -1
  id_gotm_errors      = -1

  do n=1,num_gotm

      id_gotm_diag(n) = register_diag_field ('ocean_model', trim(Gotm(n)%name),   &
                   Grd%tracer_axes_wt(1:3), Time%model_time, trim(Gotm(n)%name), &
                   &trim(Gotm(n)%units), missing_value=missing_value, range=(/-1e18,1e18/))
      id_adv_flux_x(n) = register_diag_field ('ocean_model', trim(Gotm(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho*dzt*dyt*u*Gotm_field',   &
                   &' kg/sec', missing_value=missing_value, range=(/-1e18,1e18/))
      id_adv_flux_y(n)   = register_diag_field ('ocean_model', trim(Gotm(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'rho*dzt*dxt*v*Gotm_field',     &
                    &' kg/sec', missing_value=missing_value, range=(/-1e18,1e18/))
      id_adv_flux_z(n) = register_diag_field ('ocean_model', trim(Gotm(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'rho*dxt*dyt*wt*Gotm_field',      &
                   'kg/sec', missing_value=missing_value, range=(/-1e18,1e18/))            
      id_adv_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(Gotm(n)%name)//'_xflux_adv_int_z', &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time, 'rho*dzt*dyt*u*Gotm_field',               &
                   &' kg m/sec', missing_value=missing_value, range=(/-1e18,1e18/))
      id_adv_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(Gotm(n)%name)//'_yflux_adv_int_z', &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time, 'rho*dzt*dxt*v*Gotm_field',               &
                   &' kg m/sec', missing_value=missing_value, range=(/-1e18,1e18/))
      id_gotm_errors(n) = register_diag_field ('ocean_model', trim(Gotm(n)%name)//'_advect_error',   &
                   Grd%tracer_axes_wt(1:3), Time%model_time, trim(Gotm(n)%name), &
                   &trim(Gotm(n)%units), missing_value=missing_value, range=(/-1e18,1e18/))
  enddo 

  id_gotm_nn = register_diag_field ('ocean_model', 'gotm_nn',                       &
       Grd%tracer_axes_wt(1:3), Time%model_time, 'squared buoy freq given to GOTM', &
       '1/s^2', missing_value=missing_value, range=(/-1e18,1e18/))

  id_gotm_ss = register_diag_field ('ocean_model', 'gotm_ss',                     &
       Grd%tracer_axes(1:3), Time%model_time, 'squared vert shear given to GOTM', &
       '1/s^2', missing_value=missing_value, range=(/-1e18,1e18/))

  id_diff_cbt_gotm = register_diag_field('ocean_model','diff_cbt_gotm',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert diffusivity from GOTM', 'm^2/sec',                            &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbt_gotm = register_diag_field('ocean_model','visc_cbt_gotm',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert viscosity from GOTM', 'm^2/sec',                              &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbu_gotm = register_diag_field('ocean_model','visc_cbu_gotm',Grd%vel_axes_uv(1:3),&
       Time%model_time, 'vert viscosity from GOTM', 'm^2/sec',                              &
       missing_value = missing_value, range=(/-1.e5,1.e5/))


end subroutine ocean_vert_gotm_init
! </SUBROUTINE> NAME="ocean_vert_gotm_init">


!#######################################################################
! <SUBROUTINE NAME="vert_mix_gotm_bgrid">
!
! <DESCRIPTION>
! This subroutine computes the vertical diffusivity and viscosity 
! according to the GOTM mixing model.
!
! The tke, diss, NN, and SS arrays are computed on tracer cells.
!
! Use smf_bgrid since this array uses the primary smf array read in from 
! the coupler in ocean_core/ocean_sbc.F90 when using the FMS coupler.
!
! </DESCRIPTION>
!
subroutine vert_mix_gotm_bgrid (Time, Thickness, Velocity, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_velocity_type),    intent(in) :: Velocity
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  type(ocean_density_type),     intent(in) :: Dens

  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt

  real, dimension(0:nk)  :: dz
  real, dimension(0:nk)  :: NN1d
  real, dimension(0:nk)  :: SS1d
  real                   :: active_cells
  real                   :: u_taus, u_taub
  real                   :: momflux, depth
  real                   :: rho_N2, rho_inv 
  real                   :: smf1_active, smf2_active
  real                   :: bmf1_active, bmf2_active
  real                   :: dzwur
  integer                :: i, j, k, n, m, kb, kp1
  integer                :: tau, tau_gotm_old 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL,&
      '==>Error from ocean_vert_gotm_bgrid: module must be initialized')
  endif 

  if(.not. use_this_module) return 

  ! ocean model time step label 
  tau = Time%tau

  ! cycle the gotm time step labels 
  tau_gotm_old = tau_gotm
  tau_gotm     = taup1_gotm
  taup1_gotm   = tau_gotm_old

  ! partial derivatives of density wrt to temperature and salinity 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           drhodT(i,j,k) = Dens%drhodT(i,j,k)
           drhodS(i,j,k) = Dens%drhodS(i,j,k)
        enddo
     enddo
  enddo

  ! vertical derivative of temperature and salinity at bottom of tracer cells 
  wrk3(:,:,:) = 0.0 
  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           wrk3(i,j,k) = 1.0/Thickness%dzwt(i,j,k) 
           wrk1(i,j,k) = Grd%tmask(i,j,kp1)*wrk3(i,j,k) &
                         *(T_prog(index_temp)%field(i,j,k,tau)-T_prog(index_temp)%field(i,j,kp1,tau)) 
           wrk2(i,j,k) = Grd%tmask(i,j,kp1)*wrk3(i,j,k) &
                         *(Dens%rho_salinity(i,j,k,tau)-Dens%rho_salinity(i,j,kp1,tau)) 
        enddo
     enddo
  enddo

  ! squared buoyancy frequency at bottom of T-cells 
  wrk1_v(:,:,:,1) = 0.0
  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           rho_inv = 1.0/(epsln + Dens%rho(i,j,k,tau) + Dens%rho(i,j,kp1,tau)) 
           rho_N2  = -grav*( (drhodT(i,j,k)+drhodT(i,j,kp1))*wrk1(i,j,k) &
                            +(drhodS(i,j,k)+drhodS(i,j,kp1))*wrk2(i,j,k) ) 
           wrk1_v(i,j,k,1) = rho_inv*rho_N2
        enddo
     enddo
  enddo
      
  if (map_production_gotm) then

      ! recalculate visc_cbu to get a stable scheme (Hans Burchard priv. comm.) 
      ! visc_cbt_gotm should be defined in the data domain
      visc_cbu = 0.
      do k=1,nk
         do j=jsd,jec
            do i=isd,iec
               active_cells =  Grd%tmask(i+1,j+1,k) + Grd%tmask(i+1,j,k) + &
                               Grd%tmask(i,j+1,k)   + Grd%tmask(i,j,k)   + epsln
               visc_cbu(i,j,k) = ( Grd%tmask(i+1,j+1,k)*visc_cbt_gotm(i+1,j+1,k) + &
                                   Grd%tmask(i+1,j,k)  *visc_cbt_gotm(i+1,j,k)   + &
                                   Grd%tmask(i,j+1,k)  *visc_cbt_gotm(i,j+1,k)   + &
                                   Grd%tmask(i,j,k)    *visc_cbt_gotm(i,j,k) ) / active_cells   
            enddo
         enddo
      enddo

      ! squared vertical shear at bottom of U-cells
      ! note the visc_cbu factor will be divided out during next step 
      wrk4(:,:,:)=0.0
      do k=1,nk-1
         do j=jsd,jec
            do i=isd,iec
               dzwur = 1/Thickness%dzwu(i,j,k) 
               wrk4(i,j,k) = Grd%umask(i,j,k+1)*dzwur*dzwur*visc_cbu(i,j,k)      &
                    *(                                                           &
                     ( Velocity%u(i,j,k,1,tau) - Velocity%u(i,j,k+1,1,tau) )**2  &
                    +( Velocity%u(i,j,k,2,tau) - Velocity%u(i,j,k+1,2,tau) )**2  &
                     )
            enddo
         enddo
      enddo

      ! compute the squared vertical shear at the bottom of T cells as 
      ! an average of the four nearest shears on the bottom of U cells.
      ! (do not consider land cells in the average...only active cells.)
      ! division by visc_cbt_gotm factor is needed due to wrk4 definition.  
      wrk1_v(:,:,:,2) = 0.0
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               active_cells = (Grd%umask(i,j,k+1)   + Grd%umask(i-1,j,k+1)   + &
                               Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1)) * visc_cbt_gotm(i,j,k) +epsln
               wrk1_v(i,j,k,2) = (wrk4(i,j,k)   + wrk4(i-1,j,k) + &
                                  wrk4(i,j-1,k) + wrk4(i-1,j-1,k))/active_cells
            enddo
         enddo
      enddo

  endif  ! map_production_gotm


  if (map_velocity_gotm) then
  
    ! squared vertical shear at bottom of U-cells as
    ! an average of the four nearest shears on the bottom of U cells.
    ! (do not consider land cells in the average...only active cells.)
    wrk2(:,:,:)=0.0
    wrk4(:,:,:)=0.0
    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
           active_cells = Grd%umask(i,j,k+1)   + Grd%umask(i-1,j,k+1)   +          &
                          Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1) + epsln
           wrk2(i,j,k)  = (Velocity%u(i,j,k,1,tau)   + Velocity%u(i-1,j,k,1,tau) + &
                           Velocity%u(i,j-1,k,1,tau) + Velocity%u(i-1,j-1,k,1,tau))/active_cells
           wrk4(i,j,k)  = (Velocity%u(i,j,k,2,tau)   + Velocity%u(i-1,j,k,2,tau) + &
                           Velocity%u(i,j-1,k,2,tau) + Velocity%u(i-1,j-1,k,2,tau))/active_cells
        enddo
      enddo
    enddo

    ! compute squared vertical shear at the bottom of T cells 
    wrk1_v(:,:,:,2) = 0.0
    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,2) = Grd%tmask(i,j,k+1)*wrk3(i,j,k)*wrk3(i,j,k) &
                         * (                                            &
                             ( wrk2(i,j,k) - wrk2(i,j,k+1) )**2         &
                            +( wrk4(i,j,k) - wrk4(i,j,k+1) )**2         &
                           )
        enddo
      enddo
    enddo

  endif  !map_velocity_gotm


  ! 1D gotm scheme called for each vertical column 
  ! use smf_bgrid since this is directly from the FMS coupler. 
  if(do_turbulence_gotm) then 

      do j=jsc,jec
         do i=isc,iec

            if(Grd%tmask(i,j,1) > 0.0) then 

                ! surface and bottom friction velocities on T-cell
                active_cells= Grd%umask(i,j,1)   + Grd%umask(i-1,j,1)   + &
                              Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln

                smf1_active = (Velocity%smf_bgrid(i,j,1)   + Velocity%smf_bgrid(i-1,j,1) + &
                               Velocity%smf_bgrid(i,j-1,1) + Velocity%smf_bgrid(i-1,j-1,1))/active_cells
                smf2_active = (Velocity%smf_bgrid(i,j,2)   + Velocity%smf_bgrid(i-1,j,2) + &
                               Velocity%smf_bgrid(i,j-1,2) + Velocity%smf_bgrid(i-1,j-1,2))/active_cells
                momflux = sqrt(smf1_active**2 + smf2_active**2)
                u_taus  = sqrt(momflux*rho0r) + epsln

                bmf1_active = (Velocity%bmf(i,j,1)   + Velocity%bmf(i-1,j,1) + &
                               Velocity%bmf(i,j-1,1) + Velocity%bmf(i-1,j-1,1))/active_cells
                bmf2_active = (Velocity%bmf(i,j,2)   + Velocity%bmf(i-1,j,2) + &
                               Velocity%bmf(i,j-1,2) + Velocity%bmf(i-1,j-1,2))/active_cells
                momflux = sqrt(bmf1_active**2 + bmf2_active**2)
                u_taub  = sqrt(momflux*rho0r) + epsln

                ! roughness lengths z0s and z0b have been read from nml 

                ! water column depth
                kb    = Grd%kmt(i,j)
                depth = Thickness%depth_zwt(i,j,kb)

                ! copy 3D variables into 1D arrays
                ! note the inverse ordering of k-label 
                do k=0,kb-1
                   m       = kb-k
                   tke(k)  = Gotm(index_tke)%field(i,j,m,tau_gotm)
                   eps(k)  = Gotm(index_diss)%field(i,j,m,tau_gotm)
                   L1d(k)  = cde*tke(k)**1.5/eps(k)
                   num(k)  = visc_cbt_gotm(i,j,m)
                   nuh(k)  = diff_cbt_gotm(i,j,m)
                   NN1d(k) = wrk1_v(i,j,m,1)
                   SS1d(k) = wrk1_v(i,j,m,2)
                   dz(k)   = Thickness%dzt(i,j,m)
                enddo
                ! surface and bottom conditions
                NN1d(0)  = NN1d(1)
                SS1d(0)  = SS1d(1)
                tke(kb)  = tke(kb-1)
                eps(kb)  = eps(kb-1)
                num(kb)  = num(kb-1)
                nuh(kb)  = nuh(kb-1)
                NN1d(kb) = NN1d(kb-1)
                SS1d(kb) = 0.0
                dz(kb)   = dz(kb-1)

                ! invoke GOTM
                call do_turbulence(kb, dtime, depth, u_taus, u_taub, &
                                   z0s, z0b, dz(0:kb), NN1d(0:kb), SS1d(0:kb))

                ! copy 1D arrays into 3D arrays
                do k=0,kb-1
                   m = kb-k
                   Gotm(index_tke)%field(i,j,m,taup1_gotm)  = tke(k)
                   Gotm(index_diss)%field(i,j,m,taup1_gotm) = eps(k)
                   visc_cbt_gotm(i,j,m) = num(k)
                   diff_cbt_gotm(i,j,m) = nuh(k)
                enddo

            endif  ! endif for Grd%tmask > 0

         enddo     ! end i-loop
      enddo        ! end j-loop

      if(have_obc) then
         call ocean_obc_mixing(visc_cbt_gotm, diff_cbt_gotm, &
              Gotm(index_tke)%field(:,:,:,taup1_gotm), Gotm(index_diss)%field(:,:,:,taup1_gotm))
      endif

      ! update visc_cbt_gotm to halos in preparation for 
      ! 4-point average to get viscosity on U-cell 
      call mpp_update_domains (visc_cbt_gotm(:,:,:), Dom%domain2d)

  endif

  ! perform 4-point average to get viscosity onto bottom of U-cell. 
  ! also get viscosity for C-grid velocities.  
  wrk1(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           diff_cbt(i,j,k,1) = diff_cbt_gotm(i,j,k)
           diff_cbt(i,j,k,2) = diff_cbt(i,j,k,1)
           visc_cbt(i,j,k)   = visc_cbt_gotm(i,j,k)
           active_cells =  Grd%tmask(i+1,j+1,k) + Grd%tmask(i+1,j,k) + &
                           Grd%tmask(i,j+1,k)   + Grd%tmask(i,j,k)   + epsln
           wrk1(i,j,k) = ( Grd%tmask(i+1,j+1,k)*visc_cbt_gotm(i+1,j+1,k) + &
                           Grd%tmask(i+1,j,k)  *visc_cbt_gotm(i+1,j,k)   + &
                           Grd%tmask(i,j+1,k)  *visc_cbt_gotm(i,j+1,k)   + &
                           Grd%tmask(i,j,k)    *visc_cbt_gotm(i,j,k) ) / active_cells   
           visc_cbu(i,j,k) = wrk1(i,j,k)
        enddo
     enddo
  enddo

  if(debug_this_module) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_vert_gotm_mod: checksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('diff_cbt_gotm', diff_cbt_gotm(COMP,:)*Grd%tmask(COMP,:))
      call write_chksum_3d('visc_cbt_gotm', visc_cbt_gotm(COMP,:)*Grd%tmask(COMP,:))
      call write_chksum_3d('tke', Gotm(index_tke)%field(COMP,:,tau_gotm)*Grd%tmask(COMP,:))
      call write_chksum_3d('diss', Gotm(index_diss)%field(COMP,:,tau_gotm)*Grd%tmask(COMP,:))
      call write_chksum_3d('NN', wrk1_v(COMP,:,1)*Grd%tmask(COMP,:))
      call write_chksum_3d('SS', wrk1_v(COMP,:,2)*Grd%tmask(COMP,:))
  endif

  ! send diagnostics
  do n=1,num_gotm
     call diagnose_3d(Time, Grd, id_gotm_diag(n), Gotm(n)%field(:,:,:,tau_gotm))
  enddo

  call diagnose_3d(Time, Grd, id_diff_cbt_gotm, diff_cbt_gotm(:,:,:))
  call diagnose_3d(Time, Grd, id_visc_cbt_gotm, visc_cbt_gotm(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_gotm, wrk1(:,:,:))
  call diagnose_3d(Time, Grd, id_gotm_nn, wrk1_v(:,:,:,1))
  call diagnose_3d(Time, Grd, id_gotm_ss, wrk1_v(:,:,:,2))


 end subroutine vert_mix_gotm_bgrid
 ! </SUBROUTINE> NAME="vert_mix_gotm_bgrid"



!#######################################################################
! <SUBROUTINE NAME="advect_gotm_compute">
!
! <DESCRIPTION>
! Wrapper for advection of GOTM scalar fields tke and diss. 
!
! Horizontally tke and diss are on tracer cells, so advection just as if 
! they were tracers. 
! Vertically they are between tracer cells. We do not average but shift 
! tke and diss upward to tracer points for vertical advection. At the bottom 
! tke and diss are define by the Dirichlet boundary condition.
!
! Since use a two-level time stepping scheme for tke and diss,
! it is necessary to advect these scalars with an upwind biased
! advection scheme.
!
! </DESCRIPTION>
!
subroutine advect_gotm_compute(Time, Adv_vel, Thickness, pme, river)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_adv_vel_type),   intent(in) :: Adv_vel
  type(ocean_thickness_type), intent(in) :: Thickness
  real, dimension(isd:,jsd:), intent(in) :: pme
  real, dimension(isd:,jsd:), intent(in) :: river

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from advect_gotm_compute: ocean_vert_gotm_mod not yet initialized')
  endif 

  if(.not. do_advection_gotm) return 

  ! update tke and diss to halos in preparation for advection 
  call mpp_update_domains (Gotm(index_tke)%field(:,:,:,taup1_gotm),  Dom%domain2d, complete=.false.)
  call mpp_update_domains (Gotm(index_diss)%field(:,:,:,taup1_gotm), Dom%domain2d, complete=.true.)

  if(advection_gotm_method==1) then 
    call advect_gotm_upwind(Time, Adv_vel, Thickness, pme, river)
  elseif(advection_gotm_method==2) then 
    call advect_gotm_sweby(Time, Adv_vel, Thickness, pme, river)
  endif 
!  call mpp_update_domains (Gotm(index_tke)%field(:,:,:,taup1_gotm),  Dom%domain2d)
!  call mpp_update_domains (Gotm(index_diss)%field(:,:,:,taup1_gotm), Dom%domain2d)

end subroutine advect_gotm_compute
! </SUBROUTINE> NAME="advect_gotm_compute"


!#######################################################################
! <SUBROUTINE NAME="advect_gotm_upwind">
!
! <DESCRIPTION>
! First order upwind to advect GOTM scalar turbulence quantities tke and diss.
! </DESCRIPTION>
!
subroutine advect_gotm_upwind(Time, Adv_vel, Thickness, pme, river)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_adv_vel_type),   intent(in) :: Adv_vel
  type(ocean_thickness_type), intent(in) :: Thickness
  real, dimension(isd:,jsd:), intent(in) :: pme
  real, dimension(isd:,jsd:), intent(in) :: river

  real, dimension(isd:ied,jsd:jed) :: fe, fn, ft1, ft2
  real, dimension(isd:ied,jsd:jed) :: tmp_flux
  real                             :: velocity, upos, uneg, wpos, wneg
  integer                          :: i, j, k, kp1, n
  integer                          :: tau 

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from advect_gotm_upwind: ocean_vert_gotm_mod not yet initialized')
  endif 

  flux_x  = 0.0
  flux_y  = 0.0
  flux_z  = 0.0
  tau     = Time%tau 

  do n=1,num_gotm

     gotm_tendency(:,:,:) = 0.0
     wrk1(:,:,:)          = 0.0
     wrk2(:,:,:)          = 0.0

     do k=1,nk 

        ! i-flux
        do j=jsc,jec
           do i=isc-1,iec
              velocity = 0.5*Adv_vel%uhrho_et(i,j,k)
              upos     = velocity + abs(velocity)
              uneg     = velocity - abs(velocity)
              fe(i,j)  = Grd%dyte(i,j)*(upos*Gotm(n)%field(i,j,k,taup1_gotm) + uneg*Gotm(n)%field(i+1,j,k,taup1_gotm)) &
                         *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
              flux_x(i,j,k) = fe(i,j)
           enddo
        enddo

        ! j-flux
        do j=jsc-1,jec
           do i=isc,iec
              velocity = 0.5*Adv_vel%vhrho_nt(i,j,k)
              upos     = velocity + abs(velocity)
              uneg     = velocity - abs(velocity)
              fn(i,j)  = Grd%dxtn(i,j)*(upos*Gotm(n)%field(i,j,k,taup1_gotm) + uneg*Gotm(n)%field(i,j+1,k,taup1_gotm)) &
                         *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
              flux_y(i,j,k) = fn(i,j)
           enddo
        enddo

        ! divergence of horizontal advection fluxes  
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = Grd%tmask(i,j,k)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
           enddo
        enddo

     enddo  ! k-loop


     ! vertical fluxes 
     ft1=0.0
     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsc,jec
           do i=isc,iec
              velocity = 0.5*Adv_vel%wrho_bt(i,j,k)
              wpos     = velocity + abs(velocity) 
              wneg     = velocity - abs(velocity) 
              ft2(i,j) = (wneg*Gotm(n)%field(i,j,k,taup1_gotm) + wpos*Gotm(n)%field(i,j,kp1,taup1_gotm)) &
                         *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1) 
              flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
              wrk2(i,j,k)   = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
              ft1(i,j)      = ft2(i,j)
           enddo
        enddo
     enddo

     ! include effects from pme and river advecting tke and diss into 
     ! the top cell.  Assume pme and river have tke and diss equal 
     ! to that in the top model grid cell.
     do j=jsc,jec
        do i=isc,iec
           gotm_tendency(i,j,1) = Grd%tmask(i,j,1)*(pme(i,j)+river(i,j))*Gotm(n)%field(i,j,1,taup1_gotm)
        enddo
     enddo

     ! add Gotm(n)%field*Thickness%mass_source to maintain compatibility
     ! between evolution of the scalars tke & diss, and the evolution of mass.  
     ! assume that where there is creation of mass, this mass has tke and diss
     ! equal to that already in the grid cell. 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              gotm_tendency(i,j,k) =   gotm_tendency(i,j,k) - wrk1(i,j,k) - wrk2(i,j,k) &
                                     + Gotm(n)%field(i,j,k,taup1_gotm)*Thickness%mass_source(i,j,k)
           enddo
        enddo
     enddo

     if (id_gotm_errors(n) > 0) then
     ! update the scalar Gotm field and store erroneous negative values.  
       do k=1,nk
         do j=jsc,jec
           do i=isc,iec
              Gotm(n)%field(i,j,k,taup1_gotm) =                                &
                   (Thickness%rho_dzt(i,j,k,tau)*Gotm(n)%field(i,j,k,taup1_gotm) &
                   + dtime*gotm_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
              wrk1(i,j,k) = min(Gotm(n)%field(i,j,k,taup1_gotm),0.)
           enddo
         enddo
       enddo
       call diagnose_3d(Time, Grd, id_gotm_errors(n), wrk1(:,:,:))
     else
     ! update the scalar Gotm field.  
       do k=1,nk
         do j=jsc,jec
           do i=isc,iec
              Gotm(n)%field(i,j,k,taup1_gotm) =                                &
                   (Thickness%rho_dzt(i,j,k,tau)*Gotm(n)%field(i,j,k,taup1_gotm) &
                   + dtime*gotm_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
           enddo
         enddo
       enddo
     endif
! set to minimum value if negative
     if (correct_adv_errors) then
       do k=1,nk
         do j=jsc,jec
           do i=isc,iec
              Gotm(n)%field(i,j,k,taup1_gotm) = max(Gotm(n)%field(i,j,k,taup1_gotm), Gotm(n)%min_value) 
           enddo
          enddo
       enddo
     endif

     ! diagnostics 

     call diagnose_3d(Time, Grd, id_adv_flux_x(n), flux_x(:,:,:))

     if (id_adv_flux_x_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) +  flux_x(i,j,k) 
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_adv_flux_x_int_z(n), tmp_flux(:,:))
     endif

     call diagnose_3d(Time, Grd, id_adv_flux_y(n), flux_y(:,:,:))

     if (id_adv_flux_y_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) +  flux_y(i,j,k) 
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_adv_flux_y_int_z(n), tmp_flux(:,:))
     endif

     call diagnose_3d(Time, Grd, id_adv_flux_z(n), flux_z(:,:,:))

  enddo  ! n-loop

end subroutine advect_gotm_upwind
! </SUBROUTINE> NAME="advect_gotm_upwind"



!#######################################################################
! <SUBROUTINE NAME="advect_gotm_sweby">
!
! <DESCRIPTION>
! Sweby scheme to advect GOTM scalar turbulence quantities tke and diss.
! </DESCRIPTION>
!
subroutine advect_gotm_sweby(Time, Adv_vel, Thickness, pme, river)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_adv_vel_type),   intent(in) :: Adv_vel
  type(ocean_thickness_type), intent(in) :: Thickness
  real, dimension(isd:,jsd:), intent(in) :: pme
  real, dimension(isd:,jsd:), intent(in) :: river

  real,dimension(isc:iec,jsc:jec)        :: ftp
  real,dimension(isc:iec,jsc:jec)        :: fbt
  real,dimension(isc:iec,jsc:jec)        :: wkm1
  real,dimension(isd:ied,jsd:jed)        :: tmp_flux

  integer                                :: i, j, k, n
  integer                                :: kp1, kp2, km1
  integer                                :: tau
  real                                   :: Rjm, Rj, Rjp, cfl, massflux
  real                                   :: d0, d1, thetaP, psiP 
  real                                   :: thetaM, psiM

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from advect_gotm_sweby: ocean_vert_gotm_mod not yet initialized')
  endif 

  flux_x  = 0.0
  flux_y  = 0.0
  flux_z  = 0.0
  tau     = Time%tau


  ! calculate flux at bottom face of the T-cells
  do n=1,num_gotm

     ftp           = 0.0
     fbt           = 0.0
     wkm1          = 0.0
     gotm_tendency = 0.0

     do k=1,nk

        do j=jsc,jec
           do i=isc,iec
              advect_gotm(n)%field(i,j,k) = Gotm(n)%field(i,j,k,taup1_gotm)
           enddo
        enddo

        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        km1 = max(k-1,1)

        do j=jsc,jec
           do i=isc,iec

              Rjp = (Gotm(n)%field(i,j,km1,taup1_gotm) - Gotm(n)%field(i,j,k,taup1_gotm))    &
                   *tmask_advect(i,j,km1)*tmask_advect(i,j,k)
              Rj  = (Gotm(n)%field(i,j,k,taup1_gotm) - Gotm(n)%field(i,j,kp1,taup1_gotm))    &
                   *tmask_advect(i,j,k)*tmask_advect(i,j,kp1)
              Rjm = (Gotm(n)%field(i,j,kp1,taup1_gotm) - Gotm(n)%field(i,j,kp2,taup1_gotm))  &
                   *tmask_advect(i,j,kp1)*tmask_advect(i,j,kp2)

              massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k)
              cfl = abs(Adv_vel%wrho_bt(i,j,k) * dtime  / Thickness%rho_dzt(i,j,k,tau))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              fbt(i,j) =  0.5 * ( ( massflux + abs(massflux) )       &
                   * ( Gotm(n)%field(i,j,kp1,taup1_gotm) + psiP * Rj )    &
                   + ( massflux - abs(massflux) )                    &
                   * ( Gotm(n)%field(i,j, k ,taup1_gotm) - psiM * Rj ) )  &
                   * tmask_advect(i,j,kp1) * tmask_advect(i,j,k)

              advect_gotm(n)%field(i,j,k) = advect_gotm(n)%field(i,j,k) &
                   + dtime / Thickness%rho_dzt(i,j,k,tau) * (           &
                   Grd%datr(i,j) * ( fbt(i,j) - ftp(i,j))               &
                   + Gotm(n)%field(i,j,k,taup1_gotm) *                       &
                   (wkm1(i,j) - Adv_vel%wrho_bt(i,j,k)) )

              flux_z(i,j,k) = fbt(i,j)
              ftp(i,j)      = fbt(i,j)

           enddo
        enddo

        ! update vertical velocity for next k-level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (advect_gotm(n)%field, Dom_advect_gotm%domain2d, flags=XUPDATE)

     call diagnose_3d(Time, Grd, id_adv_flux_z(n), flux_z(:,:,:))

  enddo ! end of n-loop for Gotm


  ! calculate flux at the eastern face of the T-cells
  do n=1,num_gotm

     do k=1,nk

        do j=jsc,jec
           do i=isc-1,iec

              Rjp = (advect_gotm(n)%field(i+2,j,k) - advect_gotm(n)%field(i+1,j,k))    &
                   *tmask_advect(i+2,j,k)*tmask_advect(i+1,j,k)
              Rj  = (advect_gotm(n)%field(i+1,j,k) - advect_gotm(n)%field( i ,j,k) )   &
                   *tmask_advect(i+1,j,k)*tmask_advect(i,j,k)
              Rjm = (advect_gotm(n)%field(i,j,k) - advect_gotm(n)%field(i-1,j,k))      &
                   *tmask_advect(i,j,k)*tmask_advect(i-1,j,k)

              massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)
              cfl = abs(Adv_vel%uhrho_et(i,j,k) * dtime * 2.0    &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau)) * Grd%dxte(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_x(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )  &
                   * ( advect_gotm(n)%field( i ,j,k) + psiP * Rj )   &
                   + ( massflux - abs(massflux) )                    &
                   * ( advect_gotm(n)%field(i+1,j,k) - psiM * Rj ) ) &
                   * tmask_advect(i,j,k) * tmask_advect(i+1,j,k)
           enddo
        enddo

        ! update the scalar 
        do j=jsc,jec
           do i=isc,iec
              advect_gotm(n)%field(i,j,k) = advect_gotm(n)%field(i,j,k)                         &
                   + dtime * tmask_advect(i,j,k) * Grd%datr(i,j) / Thickness%rho_dzt(i,j,k,tau) &
                   * (flux_x(i-1,j,k) - flux_x(i,j,k)                                           &
                   + Gotm(n)%field(i,j,k,taup1_gotm)* (                                           &
                     Grd%dyte(i,j) *   Adv_vel%uhrho_et( i ,j,k)                                &
                   - Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) ) )
           enddo
        enddo

     enddo ! end of k-loop


     call mpp_update_domains (advect_gotm(n)%field, Dom_advect_gotm%domain2d, flags=YUPDATE)

     call diagnose_3d(Time, Grd, id_adv_flux_x(n), flux_x(:,:,:))

     if (id_adv_flux_x_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) +  flux_x(i,j,k) 
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_adv_flux_x_int_z(n), tmp_flux(:,:))
     endif

  enddo ! end of n-loop for Gotm


  ! calculate flux at the northern face of the T-cells 
  do n=1,num_gotm

     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = 0.0
        enddo
     enddo

     do k=1,nk

        do j=jsc-1,jec
           do i=isc,iec

              Rjp = (advect_gotm(n)%field(i,j+2,k) - advect_gotm(n)%field(i,j+1,k))   &
                   *tmask_advect(i,j+2,k)*tmask_advect(i,j+1,k)
              Rj  = (advect_gotm(n)%field(i,j+1,k) - advect_gotm(n)%field(i,j,k))     &
                   *tmask_advect(i,j+1,k)*tmask_advect(i,j,k)
              Rjm = (advect_gotm(n)%field(i,j,k)   - advect_gotm(n)%field(i,j-1,k))   &
                   *tmask_advect(i,j,k)*tmask_advect(i,j-1,k)

              massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)
              cfl = abs(Adv_vel%vhrho_nt(i,j,k) * dtime * 2.0  &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau)) * Grd%dytn(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_y(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )  &
                   * ( advect_gotm(n)%field(i,j,k) + psiP * Rj )     &
                   + ( massflux - abs(massflux) )                    &
                   * ( advect_gotm(n)%field(i,j+1,k) - psiM * Rj ) ) &
                   * tmask_advect(i,j,k) * tmask_advect(i,j+1,k)

           enddo
        enddo

        ! update GOTM scalar field 
        do j=jsc,jec
           do i=isc,iec
              advect_gotm(n)%field(i,j,k) = advect_gotm(n)%field(i,j,k)                          &
                   + dtime * tmask_advect(i,j,k) * Grd%datr(i,j) / Thickness%rho_dzt(i,j,k,tau)  &
                   * (flux_y(i,j-1,k) - flux_y(i,j,k))
           enddo
        enddo

        ! calculate the overall tendency
        do j=jsc,jec
           do i=isc,iec
 
              advect_gotm(n)%field(i,j,k) = advect_gotm(n)%field(i,j,k)                        &
                   + dtime * Gotm(n)%field(i,j,k,taup1_gotm) / Thickness%rho_dzt(i,j,k,tau) * (  &
                   Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)                                          &
                   + Grd%datr(i,j)*(                                                           &
                   ( Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k)                               &
                   - Grd%dyte( i ,j) * Adv_vel%uhrho_et( i ,j,k))))

              gotm_tendency(i,j,k) = Thickness%rho_dzt(i,j,k,tau)                                     &
                                   *(advect_gotm(n)%field(i,j,k)-Gotm(n)%field(i,j,k,taup1_gotm))/dtime &
                                   *tmask_advect(i,j,k) 

           enddo
        enddo

        ! update vertical velocity for next level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     ! include effects from pme and river advecting tke and diss into 
     ! the top cell.  Assume pme and river have tke and diss equal 
     ! to that in the top model grid cell.
     do j=jsd,jed
        do i=isd,ied
           gotm_tendency(i,j,1) = gotm_tendency(i,j,1) + &
                                  Grd%tmask(i,j,1)*(pme(i,j)+river(i,j))*Gotm(n)%field(i,j,1,taup1_gotm)
        enddo
     enddo

     ! Note that we need to add the array 
     ! Gotm(n)%field*Thickness%mass_source to maintain compatibility
     ! between evolution of the scalars tke & diss, and the evolution of mass.  
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              gotm_tendency(i,j,k) = gotm_tendency(i,j,k) + &
                                     Gotm(n)%field(i,j,k,taup1_gotm)*Thickness%mass_source(i,j,k)
           enddo
        enddo
     enddo

     if (have_obc) call ocean_obc_zero_boundary(gotm_tendency, "T")

     if (id_gotm_errors(n) > 0) then
     ! update the scalar Gotm field and store erroneous negative values  
       do k=1,nk
         do j=jsc,jec
           do i=isc,iec
              Gotm(n)%field(i,j,k,taup1_gotm) =                                &
                   (Thickness%rho_dzt(i,j,k,tau)*Gotm(n)%field(i,j,k,taup1_gotm) &
                   + dtime*gotm_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
              wrk1(i,j,k) = min(Gotm(n)%field(i,j,k,taup1_gotm),0.)
           enddo
         enddo
       enddo
       call diagnose_3d(Time, Grd, id_gotm_errors(n), wrk1(:,:,:))
     else
     ! update the scalar Gotm field.  
       do k=1,nk
         do j=jsc,jec
           do i=isc,iec
              Gotm(n)%field(i,j,k,taup1_gotm) =                                &
                   (Thickness%rho_dzt(i,j,k,tau)*Gotm(n)%field(i,j,k,taup1_gotm) &
                   + dtime*gotm_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
           enddo
         enddo
       enddo
     endif
! set to minimum value if negative
     if (correct_adv_errors) then
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 Gotm(n)%field(i,j,k,taup1_gotm) = max(Gotm(n)%field(i,j,k,taup1_gotm), Gotm(n)%min_value) 
              enddo
           enddo
        enddo
     endif
     
     call diagnose_3d(Time, Grd, id_adv_flux_y(n), flux_y(:,:,:))

     if (id_adv_flux_y_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) + flux_y(i,j,k) 
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd,id_adv_flux_y_int_z(n), tmp_flux(:,:))
     endif

  enddo ! end of n-loop


end subroutine advect_gotm_sweby
! </SUBROUTINE> NAME="advect_gotm_sweby"


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_gotm_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_vert_gotm_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

  call reset_field_pointer(Got_restart, id_restart_tke,  Gotm(index_tke)%field(:,:,:,taup1_gotm) )  
  call reset_field_pointer(Got_restart, id_restart_diss, Gotm(index_diss)%field(:,:,:,taup1_gotm) )  
  call reset_field_pointer(Got_restart, id_restart_diff, diff_cbt_gotm(:,:,:) )  
  call reset_field_pointer(Got_restart, id_restart_visc, visc_cbt_gotm(:,:,:) )  

  call save_restart(Got_restart, time_stamp)

end subroutine ocean_vert_gotm_restart
! </SUBROUTINE> NAME="ocean_vert_gotm_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_gotm_end">
!
! <DESCRIPTION>
! Save the advection term for restarting the next time step. 
! </DESCRIPTION>
!
subroutine ocean_vert_gotm_end(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_vert_gotm (ocean_vert_gotm_end): NO restart written.'
    call mpp_error(WARNING,&
    '==>Warning from ocean_vert_gotm_mod (ocean_vert_gotm_end): NO restart written.')
    return
  endif 

  write(stdoutunit,*)' ' 
  write(stdoutunit,*) 'From ocean_vert_gotm_mod: ending tke and diss chksum ==>'
  call write_timestamp(Time%model_time)
  call write_chksum_3d('ending field_tke', Gotm(index_tke)%field(COMP,:,taup1_gotm)*Grd%tmask(COMP,:))
  call write_chksum_3d('ending field_diss', Gotm(index_diss)%field(COMP,:,taup1_gotm)*Grd%tmask(COMP,:))
  call write_chksum_3d('ending diff_cbt_gotm', diff_cbt_gotm(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('ending visc_cbt_gotm', visc_cbt_gotm(COMP,:)*Grd%tmask(COMP,:))

  return

end subroutine ocean_vert_gotm_end 
! </SUBROUTINE> NAME="ocean_vert_gotm_end"



end module ocean_vert_gotm_mod

