module ocean_sigma_transport_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted and density weighted time tendency for 
! tracer from transport within a bottom "sigma" layer.
! The advective portion of this routine is experimental,
! and has many problems.  It is retained in MOM for 
! exploratory use only.  Also note that the advection 
! contributes a lot of instability when running realistic
! simulations with pressure vertical coordinates. The 
! instability mechanism is unknown. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted and density weighted
! time tendency for tracer arising from 
!
! 1. Laplacian diffusion within a bottom turbulent boundary layer.
!
! 2. Upwind advection within this layer.  Advection velocities 
!    determined by model resolved velocity and parameterized 
!    downslope velocity. We use first order upwind tracer advection
!    to ensure positive definite tracer transport in the sigma 
!    layer.  As the sigma layer is a proxy for a bottom turbulent
!    boundary layer, the added mixing from the first order upwind
!    should be physically acceptable.  
!
! CAUTION: The advective portion of this algorithm has problems
! and it retained in MOM only for research purposes.  It 
! is NOT supported for general use.  
!
! The diffusivity used to determine the strength of the diffusion  
! is generally set to be a function of the local horizontal grid 
! spacing.  Diffusivity is the sum of an a priori background plus
! a velocity dependent diffusivity.  It is large if there is a  
! a heavier parcel living adjacent within the "sigma layer" above
! a lighter parcel. It is small otherwise. 
! 
! The advection is set to zero if the density is not downslope
! favorable.  That is, rho_{,x} * H_{,x} < 0 for downslope
! flow in the x-direction, and likewise in the y-direction.  
!
! The thickness of the bottom layer can span more than a single 
! bottom grid cell.  This feature allows the sigma
! layer thickness to undulate in time according to the convergence
! or divergence of mass within the sigma layer. 
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! A. Beckmann and R. Doscher, 1997: A method for improved
! representation of dense water spreading over 
! topography in geopotential--coordinate models
! Journal of Physical Oceanography, vol 27, 
! pages 581--59.
! </REFERENCE>
!
! <REFERENCE>
! R. Doscher and A. Beckmann, 2000:
! Effects of a bottom boundary layer parameterization 
! in a coarse-resolution model of the North Atlantic Ocean
! Journal of Atmospheric and Oceanic Technology, 
! vol 17 pages 698--707
! </REFERENCE>
!
! <REFERENCE>
! Campin and Goosse 1999: Parameterization of density-driven downsloping 
! flow for a coarse-resolution model in z-coordinate", Tellus 51A, 
! pages 412-430.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM (2012)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_sigma_transport_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging. 
!  </DATA> 
!  <DATA NAME="sigma_diffusion_on" TYPE="logical">
!  For using sigma diffusion. Default is true.
!  </DATA> 
!  <DATA NAME="sigma_advection_on" TYPE="logical">
!  For using sigma advection. Default is false.
!  </DATA> 
!  <DATA NAME="sigma_advection_sgs_only" TYPE="logical">
!  In many cases, adding the resolved transport to the 
!  sigma-advective transport produces a tremendous level of 
!  noise at the bottom.  The problem is that there are 
!  grid-scale features that may cause large jumps in whether
!  the velocity should be added or not, depending on the logic
!  of the scheme.  For this reason, it may be prudent to remove
!  the resolved velocity from that contributing to the sigma
!  transport scheme. Note that its removal from sigma transport
!  does not remove the contributions of the resolved velocity 
!  from the resolved advective transport arising from
!  ocean_tracer_advect_mod. It simply removes it from the 
!  added transport arising in the sigma transport module. 
!  Default is sigma_advection_sgs_only=.true. 
!  </DATA>   
!  <DATA NAME="sigma_advection_check" TYPE="logical">
!  If true, then will only include the resolved advection 
!  velocity in the sigma-layer if the direction of 
!  transport is downslope favorable for enhancing deep density.
!  IF false, then will include the velocity regardless. 
!  This option aims to reduce the large divergences
!  that occur for the case when only include the velocity 
!  if it is favorable for deep water getting more dense. 
!  Default is sigma_advection_check=.true. 
!  </DATA>   
!
!  <DATA NAME="thickness_sigma_layer" UNITS="meter" TYPE="real">
!  Initial thickness of the bottom sigma layer.   
!  </DATA>   
!  <DATA NAME="thickness_sigma_min" UNITS="meter" TYPE="real">
!  Minimum thickness of the bottom sigma layer.   
!  </DATA>   
!  <DATA NAME="thickness_sigma_max" UNITS="meter" TYPE="real">
!  Maximum thickness of the bottom sigma layer.   
!  </DATA>   

!  <DATA NAME="sigma_just_in_bottom_cell" TYPE="logical">
!  For just having sigma layer in the bottom cell, as in mom4p0. 
!  This option must be .false. in order to use sigma_advection_on=.true. 
!  Default sigma_just_in_bottom_cell=.true.
!  </DATA>   
!  <DATA NAME="tmask_sigma_on" TYPE="logical">
!  IF .true. then masks out fluxes passing into the sigma layer, except those 
!  associated with sigma transport. Typically set to .false.  
!  </DATA> 
!
!  <DATA NAME="sigma_diffusivity" UNITS="m^2/sec" TYPE="real">
!  Sigma tracer diffusivity for use if not using micom diffusivity.   
!  </DATA> 
!  <DATA NAME="sigma_diffusivity_ratio" UNITS="dimensionless" TYPE="real">
!  When flow along sigma surface is stable (i.e., heavy parcels are below lighter parcels)
!  then sigma diffusivity is reduced by sigma_diffusivity_ratio from the case where 
!  heavy parcels are above lighter parcels.  
!  </DATA>
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the sigma diffusivity is set according to a velocity scale 
!  times the grid spacing. 
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!
!  <DATA NAME="campingoose_mu" TYPE="real" UNITS="inverse seconds">
!  Dissipation rate for the bottom friction.  Campin and Goosse 
!  suggest campingoose_mu=10^-4
!  </DATA> 
!  <DATA NAME="campingoose_delta" TYPE="real" UNITS="dimensionless">
!  Fraction of a grid cell participating in the overflow process. 
!  Campin and Goosse suggest campingoose_delta=1/3. 
!  </DATA> 
!  <DATA NAME="sigma_umax" TYPE="real" UNITS="m/s">
!  Maximum downslope speed allowed in sigma layer. 
!  In some cases, the model will be unstable if sigma_umax
!  is too large.  
!  </DATA> 
!
!  <DATA NAME="smooth_sigma_velocity" TYPE="logical">
!  To smooth the sigma advective transport velocity. 
!  Default is smooth_sigma_velocity=.true. 
!  </DATA>   
!  <DATA NAME="smooth_sigma_thickness" TYPE="logical">
!  To smooth the sigma thickness. This may be needed especially 
!  for case with sigma advection, in which case the thickness 
!  can become noisy. Default is smooth_sigma_thickness=.true. 
!  </DATA>   
!  <DATA NAME="sigma_velmicom" TYPE="real" UNITS="m/s">
!  For smoothing the sigma_thickness, use this as velocity scale to
!  determine the thickness diffusivity.  
!  Default is smooth_velmicom = 0.2
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: grav, c2dbars, epsln
use diag_manager_mod, only: register_static_field, register_diag_field, send_data
use fms_mod,          only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,          only: file_exist
use fms_mod,          only: FATAL, WARNING, NOTE, stdout, stdlog
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: restart_file_type
use mpp_domains_mod,  only: mpp_update_domains, CGRID_NE
use mpp_mod,          only: input_nml_file, mpp_error

use ocean_domains_mod,     only: get_local_indices, set_ocean_domain
use ocean_density_mod,     only: density
use ocean_operators_mod,   only: FMX, FMY, FDX_T, FDY_T, BDX_ET, BDY_NT, LAP_T
use ocean_parameters_mod,  only: GEOPOTENTIAL, TERRAIN_FOLLOWING
use ocean_parameters_mod,  only: missing_value, rho0, rho0r
use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_time_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_adv_vel_type
use ocean_types_mod,       only: ocean_options_type, ocean_thickness_type, ocean_density_type
use ocean_util_mod,        only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_sum, write_chksum_2d
use ocean_tracer_util_mod, only: diagnose_3d_rho
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk6

implicit none

public ocean_sigma_transport_init
public sigma_transport
public ocean_sigma_transport_end
public ocean_sigma_transport_restart
private advect_sigma_upwind
private watermass_diag_init
private watermass_diag


private

! internally set for computing watermass diagnstics
logical :: compute_watermass_diag = .false. 

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_tracer_sigma  ! tracer concentration within sigma layer 
integer, dimension(:), allocatable  :: id_sigma_diff    ! thickness and density weighted time tendency 
                                                        ! for tracer from sigma diffusion
integer, dimension(:), allocatable  :: id_sigma_adv     ! thickness and density weighted time tendency 
                                                        ! for tracer from sigma advection
integer, dimension(:), allocatable  :: id_sigma_diff_2d ! thickness and density weighted time tendency 
                                                        ! for tracer from sigma diffusion within sigma layer 
integer, dimension(:), allocatable  :: id_sigma_adv_2d  ! thickness and density weighted time tendency 
                                                        ! for tracer from sigma advection within sigma layer 
integer, dimension(:), allocatable  :: id_sigma_smooth  ! thickness and density weighted time tendency 
                                                        ! for tracer from smoothing of the sigma thickness
integer, dimension(:), allocatable  :: id_sigma_diff_xflux       ! for i-flux from sigma diffusion
integer, dimension(:), allocatable  :: id_sigma_diff_yflux       ! for j-flux from sigma diffusion
integer, dimension(:), allocatable  :: id_sigma_diff_xflux_int_z ! for i-flux integrated in z
integer, dimension(:), allocatable  :: id_sigma_diff_yflux_int_z ! for j-flux integrated in z

integer, dimension(:), allocatable  :: id_sigma_adv_xflux       ! for i-flux from sigma advection 
integer, dimension(:), allocatable  :: id_sigma_adv_yflux       ! for j-flux from sigma advection
integer, dimension(:), allocatable  :: id_sigma_adv_xflux_int_z ! for i-flux integrated in z
integer, dimension(:), allocatable  :: id_sigma_adv_yflux_int_z ! for j-flux integrated in z

integer :: id_sigma_uhrho     =-1 ! i-directed sigma dzt_rho*advection component 
integer :: id_sigma_vhrho     =-1 ! j-directed sigma dzt_rho*advection component 
integer :: id_sigma_uhrho_sgs =-1 ! i-directed sigma dzt_rho*advection component from Campin and Goose 
integer :: id_sigma_vhrho_sgs =-1 ! j-directed sigma dzt_rho*advection component from Campin and Goose 
integer :: id_sigma_uhrho_res =-1 ! i-directed sigma dzt_rho*advection component from resolved u
integer :: id_sigma_vhrho_res =-1 ! j-directed sigma dzt_rho*advection component from resolved v

integer :: id_dtopog_dx      =-1
integer :: id_dtopog_dy      =-1 
integer :: id_tmask_sigma    =-1
integer :: id_thickness_sigma=-1
integer :: id_weight_sigma   =-1
integer :: id_diff_cet       =-1
integer :: id_diff_cnt       =-1

integer :: id_neut_rho_sigma          =-1
integer :: id_wdian_rho_sigma         =-1
integer :: id_tform_rho_sigma         =-1
integer :: id_neut_rho_sigma_on_nrho  =-1
integer :: id_wdian_rho_sigma_on_nrho =-1
integer :: id_tform_rho_sigma_on_nrho =-1

integer :: id_eta_tend_sigma     =-1
integer :: id_eta_tend_sigma_glob=-1

! for vertical coordinate 
integer :: vert_coordinate

! for restart
type(restart_file_type), save :: Sigma_restart

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, public, dimension(isd:ied,jsd:jed) :: tmask_sigma ! mask defining active sigma cells 
real, dimension(isd:ied,jsd:jed)   :: thickness_sigma   ! space-time dependent sigma layer thickness (m)
real, dimension(isd:ied,jsd:jed)   :: dtopog_dx         ! i-derivative of the T-cell bottom topography (x-slope)
real, dimension(isd:ied,jsd:jed)   :: dtopog_dy         ! j-derivative of the T-cell bottom topography (y-slope) 
real, dimension(isd:ied,jsd:jed)   :: diff_cet          ! sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(isd:ied,jsd:jed)   :: diff_cnt          ! sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(isd:ied,jsd:jed)   :: diff_max          ! maximum a priori diffusivity (m2/sec)
real, dimension(isd:ied,jsd:jed)   :: smooth_diff       ! diffusivity (m^2/sec) for smoothing thickness_sigma

real, dimension(isd:ied,jsd:jed)   :: sigma_uhrho_res   ! rho_dzt*i-velocity in sigma (m^2/s*kg/m^3)
real, dimension(isd:ied,jsd:jed)   :: sigma_vhrho_res   ! rho_dzt*j-velocity in sigma (m^2/s*kg/m^3)
real, dimension(isd:ied,jsd:jed)   :: sigma_uhrho_sgs   ! rho_dzt*i-velocity in sigma from SGS (m^2/s*kg/m^3)
real, dimension(isd:ied,jsd:jed)   :: sigma_vhrho_sgs   ! rho_dzt*j-velocity in sigma from SGS (m^2/s*kg/m^3)
real, dimension(isd:ied,jsd:jed)   :: sigma_uhrho       ! sum of sigma_uhrho_res + sigma_uhrho_sgs
real, dimension(isd:ied,jsd:jed)   :: sigma_vhrho       ! sum of sigma_uhrho_res + sigma_uhrho_sgs

#else

real, public, dimension(:,:), allocatable :: tmask_sigma     ! mask defining active sigma cells 
real, dimension(:,:),         allocatable :: thickness_sigma ! space-time dependent sigma layer thickness (m)
real, dimension(:,:),         allocatable :: dtopog_dx       ! i-derivative of T-cell bottom topography (x-slope)
real, dimension(:,:),         allocatable :: dtopog_dy       ! j-derivative of T-cell bottom topography (y-slope)
real, dimension(:,:),         allocatable :: diff_cet        ! sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(:,:),         allocatable :: diff_cnt        ! sigma diffusivity for northern face of T cell (m^2/s)
real, dimension(:,:),         allocatable :: diff_max        ! maximum a priori diffusivity (m2/sec)
real, dimension(:,:),         allocatable :: smooth_diff     ! diffusivity (m^2/sec) for smoothing thickness_sigma

real, dimension(:,:),         allocatable :: sigma_uhrho_res ! rho_dzt*i-velocity in sigma (m^2/s*kg/m^3)
real, dimension(:,:),         allocatable :: sigma_vhrho_res ! rho_dzt*j-velocity in sigma (m^2/s*kg/m^3)
real, dimension(:,:),         allocatable :: sigma_uhrho_sgs ! rho_dzt*i-velocity in sigma from SGS (m^2/s*kg/m^3)
real, dimension(:,:),         allocatable :: sigma_vhrho_sgs ! rho_dzt*j-velocity in sigma from SGS (m^2/s*kg/m^3)
real, dimension(:,:),         allocatable :: sigma_uhrho     ! sum of sigma_uhrho_res + sigma_uhrho_sgs
real, dimension(:,:),         allocatable :: sigma_vhrho     ! sum of sigma_vhrho_res + sigma_vhrho_sgs

#endif

real, dimension(:,:,:), allocatable :: tracer_sigma       ! tracer concentration in the sigma layer 
real, dimension(:,:),   allocatable :: rho_salinity_sigma ! density salinity concentration in the sigma layer 
real, dimension(:,:,:), allocatable :: fx_diff            ! rho*dzt weighted sigma diffusive flux in i-direction
real, dimension(:,:,:), allocatable :: fy_diff            ! rho*dzt weighted sigma diffusive flux in j-direction
real, dimension(:,:,:), allocatable :: fx_adv             ! rho*dzt weighted sigma advection flux in i-direction
real, dimension(:,:,:), allocatable :: fy_adv             ! rho*dzt weighted sigma advection flux in j-direction
real, dimension(:,:,:), allocatable :: tchg_adv           ! rho*dzt weighted tendency from sigma advection 
real, dimension(:,:,:), allocatable :: tchg_diff          ! rho*dzt weighted tendency from sigma diffusion 
real, dimension(:,:,:), allocatable :: tchg_smooth        ! rho*dzt weighted tendency from sigma smoothing

real, dimension(:,:,:), allocatable :: uhrho_et     ! local version of Adv_vel%uhrho_et
real, dimension(:,:,:), allocatable :: vhrho_nt     ! local version of Adv_vel%vhrho_nt

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: num_prog_tracers=-1
integer :: index_temp=-1
integer :: index_salt=-1

character(len=128) :: version=&
     '$Id: ocean_sigma_transport.F90,v 20.0 2013/12/14 00:14:30 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: unit=6  
logical :: module_is_initialized = .FALSE.

real    :: press_constant  ! a useful constant for pressure calculation 
real    :: overflow_speed  ! grav*campingoose_delta/campingoose_mu  (m/sec)
real    :: dtime_sigma     ! time step for updating the sigma layer thickness (sec) 

! for global area normalization
real    :: cellarea_r

! nml settings 
logical, public  :: tmask_sigma_on    = .false. ! mask out non-sigma diffusive fluxes passing into sigma layer
logical  :: use_this_module           = .false. ! must be true to have any sigma diffusion occur 
logical  :: debug_this_module         = .false. ! for debugging 
logical  :: verbose_init              = .true.  ! for verbose initialization printout
logical  :: sigma_diffusion_on        = .true.  ! for using sigma diffusion
logical  :: sigma_advection_on        = .false. ! for using sigma advection
logical  :: sigma_advection_sgs_only  = .true.  ! for NOT adding resolved velocity to sigma velocity
logical  :: sigma_advection_check     = .true.  ! to only include resolved velocity when downslope favorable

real     :: sigma_diffusivity         = 1.0e3   ! sigma tracer diffusivity (m^2/sec)
real     :: sigma_diffusivity_ratio   = 1.0e-6  ! ratio of min to max sigma diffusivities
real     :: thickness_sigma_max       = 100.0   ! maximum thickness (m) of bottom sigma layer
real     :: thickness_sigma_min       = 10.0    ! minimum thickness (m) of bottom sigma layer
real     :: thickness_sigma_layer     = 50.0    ! initial thickness (m) of bottom sigma layer
real     :: vel_micom                 = 0.50    ! constant velocity scale (m/s) for setting micom diffusivity
real     :: campingoose_mu            = 1.e-4   ! (sec^-1) friction for deriving SGS downslope velocity 
real     :: campingoose_delta         = 0.3333  ! fraction of a grid cell participating in dowslope advection
real     :: sigma_umax                = .10     ! (m/sec) max speed in sigma layer for downslope advection
logical  :: tracer_mix_micom          = .false. ! if true, diffusivity made a function of horz grid spacing
logical  :: sigma_just_in_bottom_cell = .true.  ! for sigma just in the bottom cell, as in mom4p0. 
logical  :: write_a_restart           = .true.  ! to reproduce across restarts when undulating thickness_sigma  
logical  :: smooth_sigma_thickness    = .true.  ! to smooth the sigma_thickness
logical  :: smooth_sigma_velocity     = .true.  ! to smooth the advective transport in the sigma layer 
real     :: smooth_velmicom           = 0.2     ! velocity scale (m/s) to smooth sigma_thickness & sigma velocity 

namelist /ocean_sigma_transport_nml/ use_this_module, debug_this_module, tmask_sigma_on, sigma_diffusion_on, &
                                     sigma_advection_on, sigma_advection_sgs_only, sigma_advection_check,    &
                                     thickness_sigma_layer, thickness_sigma_max, thickness_sigma_min,        & 
                                     sigma_diffusivity, sigma_diffusivity_ratio,                             &
                                     tracer_mix_micom, vel_micom,                                            &
                                     verbose_init, sigma_just_in_bottom_cell,                                &
                                     campingoose_mu, campingoose_delta, sigma_umax, write_a_restart,         &
                                     smooth_sigma_thickness, smooth_sigma_velocity, smooth_velmicom

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_sigma_transport_init">
!
! <DESCRIPTION>
! Initialize the sigma transport module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_sigma_transport_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, &
                                      dtime, ver_coordinate, vert_coordinate_type, debug)

  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  real,                         intent(in)           :: dtime
  integer,                      intent(in)           :: ver_coordinate
  integer,                      intent(in)           :: vert_coordinate_type
  logical,                      intent(in), optional :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, n, num
  integer :: id_restart
  real    :: dxdymn, asigma_crit
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_sigma_transport_mod (ocean_sigma_transport_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.
  vert_coordinate       = ver_coordinate 

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_sigma_transport_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_sigma_transport_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_sigma_transport_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_sigma_transport_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_sigma_transport_nml)  
  write (stdlogunit,ocean_sigma_transport_nml)

  if(use_this_module) then 
      call mpp_error(NOTE, &
      '==>Note from ocean_sigma_transport_mod: USING ocean_sigma_transport_mod.')
      Ocean_options%sigma_transport = 'Used ocean sigma transport option.'
      if(vert_coordinate_type == TERRAIN_FOLLOWING) then 
          call mpp_error(WARNING, &
          '==>Warning: ocean_sigma_transport_mod is NOT supported for TERRAIN_FOLLOWING vert coodinates.')
      endif
  else 
      call mpp_error(NOTE, &
      '==>Note from ocean_sigma_transport_mod: NOT using ocean_sigma_transport_mod.')
      Ocean_options%sigma_transport = 'Did NOT use ocean sigma transport option.'
      return
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_sigma_transport_mod with debug_this_module=.true.'  
  endif 

  write(stdoutunit,'(a,f10.2)')   &
   '==>Note: ocean_sigma_transport_mod: using forward time step of (secs)', dtime 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_sigma_transport_mod with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif 

  if(sigma_diffusion_on) then 
      write(stdoutunit,'(a)') '==>Note: ocean_sigma_transport_mod: sigma_diffusion_on=.true.'
  else
      write(stdoutunit,'(a)') '==>Note: ocean_sigma_transport_mod: sigma_diffusion_on=.false.'
  endif
  if(sigma_advection_on) then 
      write(stdoutunit,'(a)') '==>Note: ocean_sigma_transport_mod: sigma_advection_on=.true.'
      if(sigma_advection_check) then 
          write(stdoutunit,'(a)') '==>Note: ocean_sigma_transport_mod: sigma_advection_check=.true.'
      else
          write(stdoutunit,'(a)') '==>Note: ocean_sigma_transport_mod: sigma_advection_check=.false.'
      endif
  else
      write(stdoutunit,'(a)') '==>Note: ocean_sigma_transport_mod: sigma_advection_on=.false.'
  endif
  if(.not. sigma_diffusion_on .and. .not. sigma_advection_on) then 
      write(stdoutunit,'(a)') &
      '==>Note: ocean_sigma_transport_mod: sigma_advection_on=.false. AND sigma_diffusion_on=.false.'
  endif   
  if(sigma_advection_on .and. sigma_just_in_bottom_cell) then 
    call mpp_error(FATAL, &
   '==>Error ocean_sigma_transport_mod: sigma_just_in_bottom=.true. incompatible w/ sigma_advection_on=.true.')
  endif 
  if(sigma_just_in_bottom_cell) then 
     write(stdoutunit,'(a)') &
     '==>Note from ocean_sigma_transport_mod: using sigma_just_in_bottom=.true., as in mom4p0.'
  endif 
 
#ifndef MOM_STATIC_ARRAYS  
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif
  
  Dom => Domain
  Grd => Grid

  dtime_sigma      = dtime   
  num_prog_tracers = size(T_prog(:))
  press_constant   = rho0*grav*c2dbars
  overflow_speed   = grav*campingoose_delta/(epsln+campingoose_mu)
  cellarea_r       = 1.0/(epsln + Grd%tcellsurf)

  call set_ocean_domain(Dom_flux, Grid, xhalo=Dom%xhalo, yhalo=Dom%yhalo, &
                        name='flux dom sigma', maskmap=Dom%maskmap)

  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL, &
     '==>Error in ocean_sigma_transport_mod: temp and/or salt not identified in tracer array')
  endif 

  if(thickness_sigma_layer==0.0) then 
    call mpp_error(FATAL, &
     '==>Error in ocean_sigma_transport_mod: can only run this module with thickness_sigma_layer > 0.0')
  endif 
  write (stdoutunit,'(a,f12.4)') &
   '==>Note: ocean_sigma_transport_mod: initial thickness of sigma layer (m)   = ', thickness_sigma_layer 

  if(sigma_advection_on) then 
    write (stdoutunit,'(a)') &
     '==>Note: ocean_sigma_transport_mod: allowing sigma layer thickness to undulate in time as per sigma advection. '
    write (stdoutunit,'(a,f12.4)') &
     '==>Note: ocean_sigma_transport_mod: maximum thickness of sigma layer (m) = ', thickness_sigma_max
    write (stdoutunit,'(a,f12.4)') &
     '==>Note: ocean_sigma_transport_mod: minimum thickness of sigma layer (m) = ', thickness_sigma_min

    if(sigma_advection_sgs_only) then 
       write (stdoutunit,'(a)') &
       '==>Note: ocean_sigma_transport_mod: sigma-advective transport is due JUST to parameterized transport.'
    else
       write (stdoutunit,'(a)') &
       '==>Note: ocean_sigma_transport_mod: sigma-advective transport arises from resolved and SGS transport.'
    endif 

  endif 

#ifndef MOM_STATIC_ARRAYS
  allocate (diff_cet(isd:ied,jsd:jed))
  allocate (diff_cnt(isd:ied,jsd:jed))
  allocate (diff_max(isd:ied,jsd:jed))
  allocate (smooth_diff(isd:ied,jsd:jed))
  allocate (dtopog_dx(isd:ied,jsd:jed))
  allocate (dtopog_dy(isd:ied,jsd:jed))
  allocate (tmask_sigma(isd:ied,jsd:jed))
  allocate (thickness_sigma(isd:ied,jsd:jed))
  allocate (sigma_uhrho_res(isd:ied,jsd:jed))
  allocate (sigma_vhrho_res(isd:ied,jsd:jed))
  allocate (sigma_uhrho_sgs(isd:ied,jsd:jed))
  allocate (sigma_vhrho_sgs(isd:ied,jsd:jed))
  allocate (sigma_uhrho(isd:ied,jsd:jed))
  allocate (sigma_vhrho(isd:ied,jsd:jed))
#endif
  allocate (tracer_sigma(isd:ied,jsd:jed,num_prog_tracers))
  allocate (rho_salinity_sigma(isd:ied,jsd:jed))
  allocate (fx_diff(isd:ied,jsd:jed,num_prog_tracers))
  allocate (fy_diff(isd:ied,jsd:jed,num_prog_tracers))
  allocate (fx_adv(isd:ied,jsd:jed,num_prog_tracers))
  allocate (fy_adv(isd:ied,jsd:jed,num_prog_tracers))
  allocate (tchg_diff(isd:ied,jsd:jed,num_prog_tracers))
  allocate (tchg_adv(isd:ied,jsd:jed,num_prog_tracers))
  allocate (tchg_smooth(isd:ied,jsd:jed,num_prog_tracers))

  allocate (uhrho_et(isd:ied,jsd:jed,nk))
  allocate (vhrho_nt(isd:ied,jsd:jed,nk))
  uhrho_et(:,:,:) = 0.0
  vhrho_nt(:,:,:) = 0.0

  ! if wish to further mask out regions, do so here;
  ! need at least kmt=3 to make sense of the algorithm.
  tmask_sigma(:,:) = Grd%tmask(:,:,1)
  do j=jsd,jed
    do i=isd,ied
      if(Grd%kmt(i,j) < 3) tmask_sigma(i,j) = 0.0
    enddo
  enddo
  
  ! initial sigma layer thickness
  thickness_sigma(:,:) = thickness_sigma_layer*tmask_sigma(:,:) 

  ! tracer concentration in the sigma layer 
  tracer_sigma(:,:,:)     = 0.0
  rho_salinity_sigma(:,:) = 0.0

  ! sigma tracer fluxes 
  fx_diff(:,:,:) = 0.0
  fy_diff(:,:,:) = 0.0
  fx_adv(:,:,:)  = 0.0
  fy_adv(:,:,:)  = 0.0

  ! sigma tendencies 
  tchg_diff(:,:,:)   = 0.0
  tchg_adv(:,:,:)    = 0.0
  tchg_smooth(:,:,:) = 0.0

  ! sigma diffusivities 
  diff_cet(:,:) = 0.0
  diff_cnt(:,:) = 0.0
  diff_max(:,:) = sigma_diffusivity*tmask_sigma(:,:) 

  ! diffusivity for smoothing thickness of sigma layer 
  if(smooth_sigma_thickness) then 
     write (stdoutunit,'(a)') &
     '==>Note: ocean_sigma_transport_mod: smooth_sigma_thickness=.true. => diffuse sigma_thickness.'
     smooth_diff(:,:) =  tmask_sigma(:,:)*smooth_velmicom &
                        *(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))
  else 
     smooth_diff(:,:) = 0.0
  endif 

 
  ! rho_dzt * advective velocity in sigma layer 
  sigma_uhrho_res(:,:) = 0.0
  sigma_vhrho_res(:,:) = 0.0
  sigma_uhrho_sgs(:,:) = 0.0
  sigma_vhrho_sgs(:,:) = 0.0
  sigma_uhrho(:,:)     = 0.0
  sigma_vhrho(:,:)     = 0.0

  ! Micom diffusivity
  if(tracer_mix_micom) then
    do j=jsd,jed
      do i=isd,ied
        diff_max(i,j) = tmask_sigma(i,j)*vel_micom &
                        *(2.0*Grd%dxt(i,j)*Grd%dyt(i,j)/(Grd%dxt(i,j)+Grd%dyt(i,j)))
      enddo
    enddo
    if(verbose_init) then 
      do j=jsc,jec
        write (stdoutunit,'(a,i4,a,e14.7,a)')  &
         ' Laplacian diffusivity in sigma layer at (isc,',j,') = ',diff_max(isc,j),' m^2/s'
      enddo 
    endif 
  endif

  ! Gradient of bottom topography.
  dtopog_dx=0.0 ;  dtopog_dy=0.0
  do j=jsd,jed-1
     do i=isd,ied-1 
        dtopog_dx(i,j) = Grd%tmask(i+1,j,1)*(Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)
        dtopog_dy(i,j) = Grd%tmask(i,j+1,1)*(Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)
     enddo
  enddo

  id_restart = register_restart_field(Sigma_restart, 'ocean_sigma_transport.res.nc', 'thickness_sigma', &
               thickness_sigma(:,:), domain=Dom%domain2d )

  ! read in sigma thickness from restart 
  if(.NOT.file_exist('INPUT/ocean_sigma_transport.res.nc')) then
      if (.NOT. Time%init) call mpp_error(FATAL,'Expecting file INPUT/ocean_sigma_transport.res.nc to exist.&
           &This file was not found and Time%init=.false.')
  else
      call restore_state(Sigma_restart)
      call mpp_update_domains(thickness_sigma(:,:),Dom%domain2d)
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_sigma_transport_mod: initial thickness_sigma chksum'
      call write_timestamp(Time%model_time)
      call write_chksum_2d('thickness_sigma', thickness_sigma(COMP)*Grd%tmask(COMP,1))
  endif

  ! checks
  if (dtime /= 0.0) then
    num = 0
    do j=jsc,jec
      do i=isc,iec
        dxdymn = 2.0/(1.0/(Grd%dxt(i,j))**2 + 1.0/Grd%dyt(i,j)**2)
        asigma_crit = 0.5*dxdymn/dtime
        if (diff_max(i,j)*Grd%tmask(i,j,1)  > asigma_crit .and. num <= 10) then
          num = num + 1
          if (num == 1) write (unit,'(/,(1x,a))')&
          '==> Warning: Diffusive criteria exceeded. use a smaller "dtts", or "asigma".  Show 1st 10 violations:'
          write (stdoutunit,'(a,i4,a,i4,a,f6.2,a,f6.2,a,2(e14.7,a))') &
          ' at (i,j)= (',i,',',j,'), (lon,lat)= (', Grd%xt(i,j),',',Grd%yt(i,j),&
          '),  "ah" = ', diff_max(i,j),' m^2/s. the critical value =',asigma_crit,' m^2/s'
        endif
      enddo
    enddo
  endif

  ! diagnostics 

  id_dtopog_dx   = register_static_field ('ocean_model', 'dtopog_dx', Grd%tracer_axes(1:2), &
                  'X-derivative of bottom depth', 'm/m', missing_value=missing_value, range=(/-1.e3,1.e3/))
  id_dtopog_dy   = register_static_field ('ocean_model', 'dtopog_dy', Grd%tracer_axes(1:2), &
                  'Y-derivative of bottom depth', 'm/m', missing_value=missing_value, range=(/-1.e3,1.e3/))
  id_tmask_sigma = register_static_field ('ocean_model', 'tmask_sigma', Grd%tracer_axes(1:2), &
                   'Mask for bottom sigma layer', 'm', missing_value=missing_value, range=(/-1.e3,1.e5/))

  if (id_dtopog_dx  > 0) then 
    used = send_data (id_dtopog_dx, dtopog_dx(isc:iec,jsc:jec), Time%model_time)
  endif
  if (id_dtopog_dy  > 0) then 
    used = send_data (id_dtopog_dy, dtopog_dy(isc:iec,jsc:jec), Time%model_time)
  endif 
  if (id_tmask_sigma  > 0) then 
    used = send_data (id_tmask_sigma, tmask_sigma(isc:iec,jsc:jec), Time%model_time)
  endif 

  id_thickness_sigma = register_diag_field ('ocean_model', 'thickness_sigma', Grd%tracer_axes(1:2), &
                     Time%model_time, 'Thickness of BBL', 'm', missing_value=missing_value,         &
                     range=(/-1.e3,1.e3/))
  id_weight_sigma    = register_diag_field ('ocean_model', 'weight_sigma', Grd%tracer_axes(1:3), &
                     Time%model_time, 'weight of cell w/i sigma layer', 'dimensionless',         &
                     missing_value=missing_value, range=(/-10.0,10.0/))
  id_diff_cet = register_diag_field ('ocean_model', 'diff_cet_sigma', Grd%tracer_axes(1:2),   &
                Time%model_time, 'i-sigma diffusivity', 'm^2/s', missing_value=missing_value, &
                range=(/-10.0,1.e8/))
  id_diff_cnt = register_diag_field ('ocean_model', 'diff_cnt_sigma', Grd%tracer_axes(1:2),   &
                Time%model_time, 'j-sigma diffusivity', 'm^2/s', missing_value=missing_value, &
                range=(/-10.0,1.e8/))

  id_sigma_uhrho = register_diag_field ('ocean_model', 'sigma_uhrho', Grd%tracer_axes(1:2), &
                   Time%model_time, 'i-sigma uhrho', 'm^2/s * kg/m^3',                      &
                   missing_value=missing_value, range=(/-10.0,1.e8/))
  id_sigma_vhrho = register_diag_field ('ocean_model', 'sigma_vhrho', Grd%tracer_axes(1:2),  &
                   Time%model_time, 'j-sigma vhrho', 'm^2/s * kg/m^3',                       &
                   missing_value=missing_value, range=(/-1.e4,1.e4/))
  id_sigma_uhrho_sgs = register_diag_field ('ocean_model', 'sigma_uhrho_sgs', Grd%tracer_axes(1:2), &
                       Time%model_time, 'i-sigma uhrho from SGS', 'm^2/s * kg/m^3',                 &
                       missing_value=missing_value, range=(/-1.e4,1.e4/))
  id_sigma_vhrho_sgs = register_diag_field ('ocean_model', 'sigma_vhrho_sgs', Grd%tracer_axes(1:2), &
                       Time%model_time, 'j-sigma vhrho from SGS', 'm^2/s * kg/m^3',                 &
                       missing_value=missing_value, range=(/-1.e4,1.e4/))
  id_sigma_uhrho_res = register_diag_field ('ocean_model', 'sigma_uhrho_res', Grd%tracer_axes(1:2), &
                       Time%model_time, 'i-sigma uhrho via resolved flow', 'm^2/s * kg/m^3',        &
                       missing_value=missing_value, range=(/-1.e4,1.e4/))
  id_sigma_vhrho_res = register_diag_field ('ocean_model', 'sigma_vhrho_res', Grd%tracer_axes(1:2), &
                       Time%model_time, 'j-sigma vhrho via resolved flow', 'm^2/s * kg/m^3',        &
                       missing_value=missing_value, range=(/-1.e4,1.e4/))


  allocate (id_tracer_sigma(num_prog_tracers))
  allocate (id_sigma_diff(num_prog_tracers))
  allocate (id_sigma_adv(num_prog_tracers))
  allocate (id_sigma_diff_2d(num_prog_tracers))
  allocate (id_sigma_adv_2d(num_prog_tracers))
  allocate (id_sigma_smooth(num_prog_tracers))
  id_tracer_sigma  = -1
  id_sigma_diff    = -1
  id_sigma_adv     = -1
  id_sigma_diff_2d = -1
  id_sigma_adv_2d  = -1
  id_sigma_smooth  = -1
 
  do n=1,num_prog_tracers

     id_tracer_sigma(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma', &
                          Grd%tracer_axes(1:2), Time%model_time,                              &
                          trim(T_prog(n)%longname)//' in sigma layer', trim(T_prog(n)%units), &
                          missing_value=missing_value, range=(/-10.0,1.e10/))

     if (T_prog(n)%name == 'temp') then
        id_sigma_diff(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_diff', &
                           Grd%tracer_axes(1:3), Time%model_time,                                   &
                           'thk wghtd sigma-diffusion heating ', 'Watts/m^2',                       &
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_adv(n)  = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_adv', &
                           Grd%tracer_axes(1:3), Time%model_time,                                  &
                           'thk wghtd sigma-advection heating ', 'Watts/m^2',                      &
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_diff_2d(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_diff_2d', &
                           Grd%tracer_axes(1:2), Time%model_time,                                          &
                           'thk wghtd sigma-diffusion heating in sigma layer', 'Watts/m^2',                &
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_adv_2d(n)  = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_adv_2d', &
                           Grd%tracer_axes(1:2), Time%model_time,                                        &
                           'thk wghtd sigma-advection heating in sigma layer', 'Watts/m^2',              &
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_smooth(n)  = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_smooth', &
                           Grd%tracer_axes(1:2), Time%model_time,                                        &
                           'thk wghtd sigma-smoothing heating in sigma layer', 'Watts/m^2',              &
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
     else 
        id_sigma_diff(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_diff',          &
                           Grd%tracer_axes(1:3), Time%model_time,                                            &
                           'thk wghtd sigma-diffusion on '//trim(T_prog(n)%name), trim(T_prog(n)%flux_units),&
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_adv',            &
                           Grd%tracer_axes(1:3), Time%model_time,                                            &
                           'thk wghtd sigma-advection on '//trim(T_prog(n)%name), trim(T_prog(n)%flux_units),&
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_diff_2d(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_diff_2d',    &
                           Grd%tracer_axes(1:2), Time%model_time,                                            &
                           'thk wghtd sigma-diffusion in sigma layer on '//trim(T_prog(n)%name),             &
                           trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_adv_2d(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_adv_2d',      &
                           Grd%tracer_axes(1:2), Time%model_time,                                            &
                           'thk wghtd sigma-advection in sigma layer on '//trim(T_prog(n)%name),             &
                           trim(T_prog(n)%flux_units), &
                           missing_value=missing_value, range=(/-1.e16,1.e16/))
        id_sigma_smooth(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sigma_smooth',&
                           Grd%tracer_axes(1:2), Time%model_time,                                      &
                           'thk wghtd sigma-smoothing in sigma layer on '//trim(T_prog(n)%name),       &
                           trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e16,1.e16/))
     endif 
  enddo


  ! diagnostics for sigma diffusion 

  allocate (id_sigma_diff_xflux(num_prog_tracers))
  allocate (id_sigma_diff_yflux(num_prog_tracers))
  allocate (id_sigma_diff_xflux_int_z(num_prog_tracers))
  allocate (id_sigma_diff_yflux_int_z(num_prog_tracers))
  id_sigma_diff_xflux       = -1
  id_sigma_diff_yflux       = -1
  id_sigma_diff_xflux_int_z = -1
  id_sigma_diff_yflux_int_z = -1

  do n=1,num_prog_tracers
    if(n == index_temp) then 

      id_sigma_diff_xflux(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_sigma_diff_xflux',           &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,        &
              'cp*sigma_diff_xflux*dyt*rho_dzt*temp',              &
              'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
      id_sigma_diff_yflux(n) = register_diag_field ('ocean_model',     &
              trim(T_prog(n)%name)//'_sigma_diff_yflux',               &
              Grd%tracer_axes_flux_y(1:3),  Time%model_time,           &
              'cp*sigma_diff_yflux*dxt*rho_dzt*temp',                  &
              'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
      id_sigma_diff_xflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_sigma_diff_xflux_int_z',           &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,              &
              'vertical sum of cp*sigma_diff_xflux*dyt*rho_dzt*temp',    &
              'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
      id_sigma_diff_yflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_sigma_diff_yflux_int_z',           &
              Grd%tracer_axes_flux_y(1:2),  Time%model_time,             &
              'vertical sum of cp*sigma_diff_yflux*dxt*rho_dzt*temp',    &
              'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
  else

      id_sigma_diff_xflux(n) = register_diag_field ('ocean_model',              &
              trim(T_prog(n)%name)//'_sigma_diff_xflux',                        &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                     &
              'sigma_diff_xflux*dyt*rho_dzt*tracer for '//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                            & 
              range=(/-1.e16,1.e16/))
      id_sigma_diff_yflux(n) = register_diag_field ('ocean_model',              &
              trim(T_prog(n)%name)//'_sigma_diff_yflux',                        &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                     &
              'sigma_diff_yflux*dxt*rho_dzt*tracer for '//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                            &
              range=(/-1.e16,1.e16/))
      id_sigma_diff_xflux_int_z(n) = register_diag_field ('ocean_model',                        &
              trim(T_prog(n)%name)//'_sigma_diff_xflux_int_z',                                  &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                                     &
              'vertical sum of sigma_diff_xflux*dyt*rho_dzt*tracer for '//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                                            & 
              range=(/-1.e16,1.e16/))
      id_sigma_diff_yflux_int_z(n) = register_diag_field ('ocean_model',                        &
              trim(T_prog(n)%name)//'_sigma_diff_yflux_int_z',                                  &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                                     &
              'vertical sum of sigma_diff_yflux*dxt*rho_dzt*tracer for '//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                                            &
              range=(/-1.e16,1.e16/))
    endif 
  enddo

  ! diagnostics for sigma advection 

  allocate (id_sigma_adv_xflux(num_prog_tracers))
  allocate (id_sigma_adv_yflux(num_prog_tracers))
  allocate (id_sigma_adv_xflux_int_z(num_prog_tracers))
  allocate (id_sigma_adv_yflux_int_z(num_prog_tracers))
  id_sigma_adv_xflux       = -1
  id_sigma_adv_yflux       = -1
  id_sigma_adv_xflux_int_z = -1
  id_sigma_adv_yflux_int_z = -1

  do n=1,num_prog_tracers
    if(n == index_temp) then 

      id_sigma_adv_xflux(n) = register_diag_field ('ocean_model',  &
               trim(T_prog(n)%name)//'_sigma_adv_xflux',           &
               Grd%tracer_axes_flux_x(1:3), Time%model_time,       &
               'cp*sigma_adv_xflux*dyt*rho_dzt*temp',              &
               'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
      id_sigma_adv_yflux(n) = register_diag_field ('ocean_model',     &
              trim(T_prog(n)%name)//'_sigma_adv_yflux',               &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,           &
              'cp*sigma_adv_yflux*dxt*rho_dzt*temp', &
              'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
      id_sigma_adv_xflux_int_z(n) = register_diag_field ('ocean_model',  &
               trim(T_prog(n)%name)//'_sigma_adv_xflux_int_z',           &
               Grd%tracer_axes_flux_x(1:2), Time%model_time,             &
               'vertical sum of cp*sigma_adv_xflux*dyt*rho_dzt*temp',    &
               'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))
      id_sigma_adv_yflux_int_z(n) = register_diag_field ('ocean_model',  &
              trim(T_prog(n)%name)//'_sigma_adv_yflux_int_z',            &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,              &
              'vertical sum of cp*sigma_adv_yflux*dxt*rho_dzt*temp',     &
              'Watt', missing_value=missing_value, range=(/-1.e16,1.e16/))

  else

      id_sigma_adv_xflux(n) = register_diag_field ('ocean_model',             &
              trim(T_prog(n)%name)//'_sigma_adv_xflux',                       &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                   &
              'sigma_adv_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                          & 
              range=(/-1.e16,1.e16/))
      id_sigma_adv_yflux(n) = register_diag_field ('ocean_model',             &
              trim(T_prog(n)%name)//'_sigma_adv_yflux',                       &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                   &
              'sigma_adv_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                          &
              range=(/-1.e16,1.e16/))
      id_sigma_adv_xflux_int_z(n) = register_diag_field ('ocean_model',                       &
              trim(T_prog(n)%name)//'_sigma_adv_xflux_int_z',                                 &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                                   &
              'vertical sum of sigma_adv_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                                          & 
              range=(/-1.e16,1.e16/))
      id_sigma_adv_yflux_int_z(n) = register_diag_field ('ocean_model',                       &
              trim(T_prog(n)%name)//'_sigma_adv_yflux_int_z',                                 &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                                   &
              'vertical sum of sigma_adv_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              'kg/sec', missing_value=missing_value,                                          &
              range=(/-1.e16,1.e16/))
    endif 
  enddo

  call watermass_diag_init(Time, Dens)

end subroutine ocean_sigma_transport_init
! </SUBROUTINE> NAME="ocean_sigma_transport_init"


!#######################################################################
! <SUBROUTINE NAME="sigma_transport">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted and density 
! weighted time tendency for tracer arising from transport in a 
! bottom turbulent boundary layer. The result is stored in 
! tracer th_tendency. 
!
! NOTE: In this algorithm, we ideally wish to have advection 
! velocity components on full data domain.  Unfortunately, 
! from ocean_advection_velocity_mod, they are only known 
! on the following domains:
!
! Adv_vel%uhrho_et: (isd,ied) x (jsc,jed)  
! Adv_vel%vhrho_nt: (isc,ied) x (jsd,jed).  
!
! So to proceed with the sigma_transport algorithm, we 
! transfer into local arrays and then update.  These 
! updates may be avoided (possibly), but at the price
! of much more logic in the algorithm.  We choose to 
! instead do the updates and have less logic.  This 
! decision may need to be revisited. 
!
! CAUTION: The advective portion of this algorithm 
! has fundamental problems.  It is retained in MOM
! only for process physics research purposes.  It is
! NOT recommended for use in general use. 
! 
! </DESCRIPTION>
subroutine sigma_transport (Time, Thickness, Dens, T_prog, Adv_vel, bott_blthick)
 
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  real, dimension(isd:,jsd:),   intent(inout) :: bott_blthick

  real, dimension(isd:ied,jsd:jed) :: thickness_sigma_eff
  real, dimension(isd:ied,jsd:jed) :: thickness_sigma_old
  real, dimension(isd:ied,jsd:jed) :: dzt_sigma_tendency
  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(2)               :: density_test

  real :: depth_sigma
  real :: check, press
  real :: temp0, salt0, tempip1, saltip1, tempjp1, saltjp1
  real :: ratio_depth, rho_dzt_sigma

  integer :: i, j, k, n
  integer :: kip1, kjp1, kmt, ksigma
  integer :: taum1, tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. use_this_module) return

  if (size(T_prog(:)) /= num_prog_tracers) then 
      call mpp_error(FATAL, &
           '==>Error from ocean_sigma_transport_mod (sigma_transport): size mismatch for tracer array')
  endif
  if (.not. module_is_initialized ) then 
      call mpp_error(FATAL, &
           '==>Error from ocean_sigma_transport_mod (sigma_transport): module must be initialized')
  endif

  taum1 = Time%taum1
  tau   = Time%tau

  ! initialize properties within sigma layer 
  thickness_sigma_eff(:,:) = 0.0
  wrk1(:,:,:)              = 0.0
  tracer_sigma(:,:,:)      = 0.0
  rho_salinity_sigma(:,:)  = 0.0
  sigma_uhrho_res(:,:)     = 0.0     
  sigma_vhrho_res(:,:)     = 0.0     
  dzt_sigma_tendency(:,:)  = 0.0
  rho_dzt_sigma            = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           uhrho_et(i,j,k)       = Adv_vel%uhrho_et(i,j,k)
           vhrho_nt(i,j,k)       = Adv_vel%vhrho_nt(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(uhrho_et(:,:,:),vhrho_nt(:,:,:),Dom_flux%domain2d, gridtype=CGRID_NE) 

  if(sigma_just_in_bottom_cell) then 
  ! as in mom4p0, sigma layer is just within the bottom 
  ! grid cell, no matter how thin the cell may be. 

      thickness_sigma(:,:) = 0.0 
      do j=jsd,jed
         do i=isd,ied
            if(tmask_sigma(i,j)==1.0) then 
                k                        = Grd%kmt(i,j)
                thickness_sigma(i,j)     = min(thickness_sigma_layer,Thickness%dzt(i,j,k))
                thickness_sigma_eff(i,j) = thickness_sigma(i,j) 
                wrk1(i,j,k)              = 1.0
                do n=1,num_prog_tracers 
                   tracer_sigma(i,j,n)   = T_prog(n)%field(i,j,k,taum1)
                enddo
                rho_salinity_sigma(i,j)  = Dens%rho_salinity(i,j,k,taum1)
            endif 
         enddo
      enddo

  else 
  ! allow sigma layer to be within more than just the bottom cell.
  ! doing so requires some extra accounting to see where the sigma
  ! layer lives.   

      do j=jsd,jed
         do i=isd,ied
            if(tmask_sigma(i,j)==1.0) then 

                kmt           = Grd%kmt(i,j)
                depth_sigma   = Thickness%depth_zwt(i,j,kmt) - thickness_sigma(i,j)
                rho_dzt_sigma = 0.0

                ! boundary layer fills the total water column to top of k=1 cell 
                if(depth_sigma <= 0.0) then 
                    ratio_depth = 1.0
                    ksigma      = 1 
                    do k=ksigma,kmt
                       thickness_sigma_eff(i,j) = thickness_sigma_eff(i,j) + Thickness%dzt(i,j,k)
                       rho_dzt_sigma            = rho_dzt_sigma            + Thickness%rho_dzt(i,j,k,tau)
                       sigma_uhrho_res(i,j)     = sigma_uhrho_res(i,j)     + uhrho_et(i,j,k) 
                       sigma_vhrho_res(i,j)     = sigma_vhrho_res(i,j)     + vhrho_nt(i,j,k) 
                       do n=1,num_prog_tracers 
                          tracer_sigma(i,j,n) = tracer_sigma(i,j,n) &
                               +Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%field(i,j,k,taum1)
                       enddo
                       rho_salinity_sigma(i,j) = rho_salinity_sigma(i,j) &
                               +Thickness%rho_dzt(i,j,k,tau)*Dens%rho_salinity(i,j,k,taum1)
                    enddo
                    do k=ksigma,kmt
                       wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)/rho_dzt_sigma
                    enddo

                ! boundary layer is within k=1 cell 
                elseif(depth_sigma > 0.0 .and. &
                     depth_sigma <= Thickness%depth_zwt(i,j,1)) then 
                    ksigma=1 
                    do k=ksigma+1,kmt
                       thickness_sigma_eff(i,j) = thickness_sigma_eff(i,j) + Thickness%dzt(i,j,k)
                       rho_dzt_sigma            = rho_dzt_sigma            + Thickness%rho_dzt(i,j,k,tau)
                       sigma_uhrho_res(i,j)     = sigma_uhrho_res(i,j)     + uhrho_et(i,j,k) 
                       sigma_vhrho_res(i,j)     = sigma_vhrho_res(i,j)     + vhrho_nt(i,j,k) 
                       do n=1,num_prog_tracers 
                          tracer_sigma(i,j,n) = tracer_sigma(i,j,n) &
                               +Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%field(i,j,k,taum1)
                       enddo
                       rho_salinity_sigma(i,j) = rho_salinity_sigma(i,j) &
                               +Thickness%rho_dzt(i,j,k,tau)*Dens%rho_salinity(i,j,k,taum1)
                    enddo
                    k=ksigma
                    ratio_depth = (Thickness%depth_zwt(i,j,k)-depth_sigma) &
                                 /(Thickness%depth_zwt(i,j,k) + epsln)
                    thickness_sigma_eff(i,j) = thickness_sigma_eff(i,j) + ratio_depth*Thickness%dzt(i,j,k)
                    rho_dzt_sigma            = rho_dzt_sigma            + ratio_depth*Thickness%rho_dzt(i,j,k,tau)
                    sigma_uhrho_res(i,j)     = sigma_uhrho_res(i,j)     + ratio_depth*uhrho_et(i,j,k) 
                    sigma_vhrho_res(i,j)     = sigma_vhrho_res(i,j)     + ratio_depth*vhrho_nt(i,j,k) 
                    do n=1,num_prog_tracers 
                       tracer_sigma(i,j,n) = tracer_sigma(i,j,n) &
                            +ratio_depth*Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%field(i,j,k,taum1)
                    enddo
                    rho_salinity_sigma(i,j) = rho_salinity_sigma(i,j) &
                            +ratio_depth*Thickness%rho_dzt(i,j,k,tau)*Dens%rho_salinity(i,j,k,taum1)

                    wrk1(i,j,k)    = ratio_depth*Thickness%rho_dzt(i,j,k,tau)/rho_dzt_sigma
                    do k=ksigma+1,kmt
                       wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)/rho_dzt_sigma
                    enddo

                ! boundary layer is fully contained within bottom cell 
                elseif(depth_sigma >  Thickness%depth_zwt(i,j,kmt-1) .and. &
                       depth_sigma <= Thickness%depth_zwt(i,j,kmt)) then  
                    ksigma=kmt
                    k=ksigma
                    ratio_depth = (Thickness%depth_zwt(i,j,k)-depth_sigma) &
                                 /(Thickness%depth_zwt(i,j,k)-Thickness%depth_zwt(i,j,k-1)) 
                    thickness_sigma_eff(i,j) = ratio_depth*Thickness%dzt(i,j,k)
                    rho_dzt_sigma            = ratio_depth*Thickness%rho_dzt(i,j,k,tau)
                    sigma_uhrho_res(i,j)     = ratio_depth*uhrho_et(i,j,k) 
                    sigma_vhrho_res(i,j)     = ratio_depth*vhrho_nt(i,j,k) 
                    do n=1,num_prog_tracers 
                       tracer_sigma(i,j,n)   = rho_dzt_sigma*T_prog(n)%field(i,j,k,taum1)
                    enddo
                    rho_salinity_sigma(i,j)  = rho_dzt_sigma*Dens%rho_salinity(i,j,k,taum1)
                    wrk1(i,j,k) = 1.0

                ! boundary layer is no shallower than k=2 and is not fully in k=kmt
                else

                    kloop1:    do k=2,kmt-1
                       if(depth_sigma >  Thickness%depth_zwt(i,j,k-1) .and. &
                          depth_sigma <= Thickness%depth_zwt(i,j,k)) then  
                          ksigma = k
                          exit kloop1
                       endif
                    enddo kloop1
                    do k=ksigma+1,kmt
                       thickness_sigma_eff(i,j) = thickness_sigma_eff(i,j) + Thickness%dzt(i,j,k)
                       rho_dzt_sigma            = rho_dzt_sigma            + Thickness%rho_dzt(i,j,k,tau)
                       sigma_uhrho_res(i,j)     = sigma_uhrho_res(i,j)     + uhrho_et(i,j,k) 
                       sigma_vhrho_res(i,j)     = sigma_vhrho_res(i,j)     + vhrho_nt(i,j,k) 
                       do n=1,num_prog_tracers 
                          tracer_sigma(i,j,n) = tracer_sigma(i,j,n) &
                               +Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%field(i,j,k,taum1)
                       enddo
                       rho_salinity_sigma(i,j) = rho_salinity_sigma(i,j) &
                               +Thickness%rho_dzt(i,j,k,tau)*Dens%rho_salinity(i,j,k,taum1)
                    enddo
                    k=ksigma 
                    ratio_depth = (Thickness%depth_zwt(i,j,k)-depth_sigma) &
                                 /(Thickness%depth_zwt(i,j,k)-Thickness%depth_zwt(i,j,k-1)) 
                    rho_dzt_sigma            = rho_dzt_sigma            + ratio_depth*Thickness%rho_dzt(i,j,k,tau)
                    thickness_sigma_eff(i,j) = thickness_sigma_eff(i,j) + ratio_depth*Thickness%dzt(i,j,k)
                    sigma_uhrho_res(i,j)     = sigma_uhrho_res(i,j)     + ratio_depth*uhrho_et(i,j,k) 
                    sigma_vhrho_res(i,j)     = sigma_vhrho_res(i,j)     + ratio_depth*vhrho_nt(i,j,k) 
                    do n=1,num_prog_tracers 
                       tracer_sigma(i,j,n) = tracer_sigma(i,j,n) &
                            +ratio_depth*Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%field(i,j,k,taum1)
                    enddo
                    rho_salinity_sigma(i,j) = rho_salinity_sigma(i,j) &
                               +ratio_depth*Thickness%rho_dzt(i,j,k,tau)*Dens%rho_salinity(i,j,k,taum1)
                    wrk1(i,j,k)    = ratio_depth*Thickness%rho_dzt(i,j,k,tau)/rho_dzt_sigma
                    do k=ksigma+1,kmt
                       wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)/rho_dzt_sigma
                    enddo

                endif  ! endif for depth_sigma 

                ! divide by rho_dzt_sigma to produce a tracer concentration
                do n=1,num_prog_tracers 
                   tracer_sigma(i,j,n) = tracer_sigma(i,j,n)/rho_dzt_sigma
                enddo
                rho_salinity_sigma(i,j) = rho_salinity_sigma(i,j)/rho_dzt_sigma

            endif      ! endif for tmask_sigma 
         enddo         ! i-loop  
      enddo            ! j-loop 

  endif                ! endif for sigma_just_in_bottom_cell


  ! set the resolved transport to zero 
  if(sigma_advection_sgs_only) then 
      sigma_uhrho_res(:,:) = 0.0
      sigma_vhrho_res(:,:) = 0.0
  endif

  ! determine diffusivities and SGS velocity within bottom sigma layer.
  ! allow enhanced diffusion even if topography is flat (check=0).
  ! use deepest depth to approximate hydrostatic pressure for computing density. 

  diff_cet(:,:)        = 0.0
  diff_cnt(:,:)        = 0.0
  sigma_uhrho_sgs(:,:) = 0.0     
  sigma_vhrho_sgs(:,:) = 0.0     
  sigma_uhrho(:,:)     = 0.0     
  sigma_vhrho(:,:)     = 0.0     

  do j=jsd,jed-1
     do i=isd,ied-1

        if(tmask_sigma(i,j)==1.0) then 

            k=Grd%kmt(i,j)

            temp0   = tracer_sigma(i,j,  index_temp)
            salt0   = rho_salinity_sigma(i,j)
            tempip1 = tracer_sigma(i+1,j,index_temp)
            saltip1 = rho_salinity_sigma(i+1,j) 
            tempjp1 = tracer_sigma(i,j+1,index_temp)
            saltjp1 = rho_salinity_sigma(i,j+1)

            kip1  = max(1,Grd%kmt(i+1,j))
            press = press_constant*max(Thickness%depth_zt(i,j,k),Thickness%depth_zt(i+1,j,kip1))
            density_test(1) = density(salt0,temp0,press)
            density_test(2) = density(saltip1,tempip1,press)
            check = (density_test(2)-density_test(1))*dtopog_dx(i,j)        

            ! density favorable for downslope flow 
            if(check <= 0.0) then 
                diff_cet(i,j)        = diff_max(i,j)              
                sigma_uhrho_sgs(i,j) = -overflow_speed*check*thickness_sigma_eff(i,j)*sign(1.0,dtopog_dx(i,j))
            else 
                diff_cet(i,j) = diff_max(i,j)*sigma_diffusivity_ratio              
            endif

            ! density favorable for downslope flow & resolved velocity in downslope direction 
            if(check <= 0.0 .and. sigma_uhrho_res(i,j)*dtopog_dx(i,j) >= 0.0) then 
                diff_cet(i,j) = & 
                min(diff_cet(i,j)+abs(sigma_uhrho_res(i,j))*Grd%dxt(i,j)/Thickness%rho_dzt(i,j,k,tau), &
                2.0*diff_max(i,j))
                sigma_uhrho(i,j) = sigma_uhrho_sgs(i,j) + sigma_uhrho_res(i,j)
            else 
                sigma_uhrho(i,j) = sigma_uhrho_sgs(i,j)
            endif

            ! always include resolved velocity  
            if(.not. sigma_advection_check) then 
                sigma_uhrho(i,j) = sigma_uhrho_sgs(i,j) + sigma_uhrho_res(i,j)
            endif 

            kjp1  = max(1,Grd%kmt(i,j+1))
            press = press_constant*max(Thickness%depth_zt(i,j,k),Thickness%depth_zt(i,j+1,kjp1))
            density_test(1) = density(salt0,temp0,press)     
            density_test(2) = density(saltjp1,tempjp1,press)
            check = (density_test(2)-density_test(1))*dtopog_dy(i,j)        

            ! density favorable for downslope flow 
            if(check <= 0.0) then 
                diff_cnt(i,j)        = diff_max(i,j)              
                sigma_vhrho_sgs(i,j) = -overflow_speed*check*thickness_sigma_eff(i,j)*sign(1.0,dtopog_dy(i,j))
            else 
                diff_cnt(i,j) = diff_max(i,j)*sigma_diffusivity_ratio               
            endif

            ! density favorable for downslope flow & resolved velocity in downslope direction 
            if(check <= 0.0 .and. sigma_vhrho_res(i,j)*dtopog_dy(i,j) >= 0.0) then 
                diff_cnt(i,j) = &
                min(diff_cnt(i,j)+abs(sigma_vhrho_res(i,j))*Grd%dyt(i,j)/Thickness%rho_dzt(i,j,k,tau), &
                2.0*diff_max(i,j))     
                sigma_vhrho(i,j) = sigma_vhrho_sgs(i,j) + sigma_vhrho_res(i,j)
            else 
                sigma_vhrho(i,j) = sigma_vhrho_sgs(i,j)
            endif

            ! always include resolved velocity  
            if(.not. sigma_advection_check) then 
                sigma_vhrho(i,j) = sigma_vhrho_sgs(i,j) + sigma_vhrho_res(i,j)
            endif 

            ! check for huge sigma_uhrho and huge sigma_vhrho
            check = rho0*thickness_sigma_eff(i,j)*sigma_umax
            if(sigma_uhrho(i,j) < -check .or. check < sigma_uhrho(i,j)) then 
                sigma_uhrho(i,j) = check*sign(1.0,sigma_uhrho(i,j))
            endif
            if(sigma_vhrho(i,j) < -check .or. check < sigma_vhrho(i,j)) then 
                sigma_vhrho(i,j) = check*sign(1.0,sigma_vhrho(i,j))
            endif


        endif  ! endif for tmask_sigma 

     enddo     ! i-loop
  enddo        ! j-loop 

  ! update sigma layer thickness and compute uhrho and vhrho
  if(.not. sigma_advection_on) then 

      do j=jsd,jed
         do i=isd,ied
            if(tmask_sigma(i,j)==1.0) then 

                thickness_sigma(i,j) = thickness_sigma_eff(i,j)
                if(thickness_sigma_eff(i,j) > thickness_sigma_max) then 
                    thickness_sigma(i,j) = thickness_sigma_max
                endif
                if(thickness_sigma_eff(i,j) < thickness_sigma_min) then 
                    thickness_sigma(i,j) = thickness_sigma_min
                endif

            endif
         enddo
      enddo

  else 

      if(smooth_sigma_velocity) then 
          call mpp_update_domains(sigma_uhrho(:,:),sigma_vhrho(:,:),Dom_flux%domain2d, gridtype=CGRID_NE) 
          sigma_uhrho(:,:) = sigma_uhrho(:,:) + dtime_sigma*LAP_T(sigma_uhrho(:,:),smooth_diff(:,:))
          sigma_vhrho(:,:) = sigma_vhrho(:,:) + dtime_sigma*LAP_T(sigma_vhrho(:,:),smooth_diff(:,:))
      endif

      call mpp_update_domains(sigma_uhrho(:,:),sigma_vhrho(:,:),Dom_flux%domain2d, gridtype=CGRID_NE) 
      dzt_sigma_tendency(:,:) = -rho0r*tmask_sigma(:,:)*(BDX_ET(sigma_uhrho(:,:))+BDY_NT(sigma_vhrho(:,:))) 

      if(smooth_sigma_thickness) then 
          dzt_sigma_tendency(:,:) = dzt_sigma_tendency(:,:) &
                                  + dtime_sigma*LAP_T(thickness_sigma(:,:),smooth_diff(:,:))

          ! this process changes thickness diffusively, so add
          ! a term to the tracer to maintain compatibility.  
          do n=1,num_prog_tracers
             do j=jsd,jed
                do i=isd,ied
                   tmp(i,j) = thickness_sigma(i,j)*tracer_sigma(i,j,n)
                enddo
             enddo
             tmp(:,:) = LAP_T(tmp(:,:),smooth_diff(:,:), update=.true.)
             do j=jsc,jec
                do i=isc,iec
                   tchg_smooth(i,j,n) = dtime_sigma*tmp(i,j)
                enddo
             enddo
          enddo
      endif

      ! update the thickness_sigma array  
      do j=jsc,jec
         do i=isc,iec
            thickness_sigma_old(i,j) = thickness_sigma(i,j)
            thickness_sigma(i,j)     = thickness_sigma_old(i,j) + dtime_sigma*dzt_sigma_tendency(i,j)
         enddo
      enddo

      ! ensure thickness is within bounds 
      do j=jsc,jec
         do i=isc,iec
            if(thickness_sigma(i,j) < thickness_sigma_min) then 
                thickness_sigma(i,j) = thickness_sigma_min
            endif
            if(thickness_sigma(i,j) > thickness_sigma_max) then 
                thickness_sigma(i,j) = thickness_sigma_max
            endif
         enddo
      enddo
      call mpp_update_domains(thickness_sigma(:,:), Dom%domain2d)


  endif  ! endif for sigma_advection_on

  if(debug_this_module) then 
     write(stdoutunit,*) 'From ocean_sigma_transport_mod: intermediate thickness_sigma chksum'
     call write_timestamp(Time%model_time)
     call write_chksum_2d('thickness_sigma', thickness_sigma(COMP)*Grd%tmask(COMP,1))
  endif 

  ! for use in neutral physics 
  do j=jsd,jed
     do i=isd,ied
        bott_blthick(i,j) = thickness_sigma(i,j)
     enddo
  enddo


  ! compute time tendency from sigma advection 
  if(sigma_advection_on) then 
      call advect_sigma_upwind(T_prog(1:num_prog_tracers))

      tchg_adv(:,:,:) = 0.0
      do n=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec
               tchg_adv(i,j,n) = -tmask_sigma(i,j) &
                    *(fx_adv(i,j,n)-fx_adv(i-1,j,n)+fy_adv(i,j,n)-fy_adv(i,j-1,n))*Grd%datr(i,j)
            enddo
         enddo
      enddo
  endif


  ! compute time tendency from sigma diffusion 
  if(sigma_diffusion_on) then 

      fx_diff(:,:,:)   = 0.0
      fy_diff(:,:,:)   = 0.0
      tchg_diff(:,:,:) = 0.0

      call mpp_update_domains(diff_cet(:,:), Dom%domain2d)
      call mpp_update_domains(diff_cnt(:,:), Dom%domain2d)

      do n=1,num_prog_tracers

         ! diffusive flux 
         fx_diff(:,:,n) = rho0*diff_cet(:,:)*FDX_T(tracer_sigma(:,:,n)) &
                          *FMX(thickness_sigma_eff(:,:)*tmask_sigma(:,:))
         fy_diff(:,:,n) = rho0*diff_cnt(:,:)*FDY_T(tracer_sigma(:,:,n)) &
                          *FMY(thickness_sigma_eff(:,:)*tmask_sigma(:,:))

         !for redundancies at Arctic fold 
         if (Grd%tripolar) then 
             call mpp_update_domains(fx_diff(:,:,n),fy_diff(:,:,n), Dom_flux%domain2d, &
                                     gridtype=CGRID_NE, complete=T_prog(n)%complete) 
         endif

      enddo

      ! diffusive tendency  
      do n=1,num_prog_tracers
         tchg_diff(:,:,n) = tmask_sigma(:,:)*(BDX_ET(fx_diff(:,:,n)) + BDY_NT(fy_diff(:,:,n)))
      enddo

  endif


  do n=1,num_prog_tracers

     ! th_tendency has dimensions = (tracer concentration)*(m/s)*(kg/m^3)
     wrk2(:,:,:) = 0.0
     wrk3(:,:,:) = 0.0
     wrk4(:,:,:) = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk2(i,j,k)                  = wrk1(i,j,k)*tchg_diff(i,j,n)
              wrk3(i,j,k)                  = wrk1(i,j,k)*tchg_adv(i,j,n)
              wrk4(i,j,k)                  = wrk1(i,j,k)*tchg_smooth(i,j,n)
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k)  &
                                             + wrk2(i,j,k) + wrk3(i,j,k) + wrk4(i,j,k)
           enddo
        enddo
     enddo

     ! send tracer concentration to diagnostic manager 
     call diagnose_2d(Time, Grd, id_tracer_sigma(n), tracer_sigma(:,:,n))

     ! send tendency to diagnostic manager 
     if (id_sigma_diff(n) > 0) then 
        call diagnose_3d(Time, Grd, id_sigma_diff(n), T_prog(n)%conversion*wrk2(:,:,:))
     endif
     if (id_sigma_adv(n) > 0) then 
        call diagnose_3d(Time, Grd, id_sigma_adv(n), T_prog(n)%conversion*wrk3(:,:,:))
     endif
     if (id_sigma_diff_2d(n) > 0) then 
        call diagnose_2d(Time, Grd, id_sigma_diff_2d(n), T_prog(n)%conversion*tchg_diff(:,:,n))
     endif
     if (id_sigma_adv_2d(n) > 0) then 
        call diagnose_2d(Time, Grd, id_sigma_adv_2d(n), T_prog(n)%conversion*tchg_adv(:,:,n))
     endif


     ! minus sign accounts for MOM sign convention for diffusive flux.  
     if(id_sigma_diff_xflux(n) > 0) then 
         tmp(:,:) = -T_prog(n)%conversion*Grd%dyte(:,:)*fx_diff(:,:,n)
         do k=1,nk
            wrk2(:,:,k) = wrk1(:,:,k)*tmp(:,:) 
         enddo
         call diagnose_3d(Time, Grd, id_sigma_diff_xflux(n), wrk2(:,:,:))
     endif
     if(id_sigma_diff_yflux(n) > 0) then 
         tmp(:,:) = -T_prog(n)%conversion*Grd%dxtn(:,:)*fy_diff(:,:,n)
         do k=1,nk
            wrk2(:,:,k) = wrk1(:,:,k)*tmp(:,:) 
         enddo
         call diagnose_3d(Time, Grd, id_sigma_diff_yflux(n), wrk2(:,:,:))
     endif

     ! advective fluxes have minus sign built in
     if(id_sigma_adv_xflux(n) > 0) then 
         tmp(:,:) = T_prog(n)%conversion*Grd%dyte(:,:)*fx_adv(:,:,n)
         do k=1,nk
            wrk3(:,:,k) = wrk1(:,:,k)*tmp(:,:) 
         enddo
         call diagnose_3d(Time, Grd, id_sigma_adv_xflux(n), wrk3(:,:,:))
     endif
     if(id_sigma_adv_yflux(n) > 0) then 
         tmp(:,:) = T_prog(n)%conversion*Grd%dxtn(:,:)*fy_adv(:,:,n)
         do k=1,nk
            wrk3(:,:,k) = wrk1(:,:,k)*tmp(:,:) 
         enddo
         call diagnose_3d(Time, Grd, id_sigma_adv_yflux(n), wrk3(:,:,:))
     endif

     ! vertically integrated flux from sigma transport = flux in sigma layer. 
     if(id_sigma_diff_xflux_int_z(n) > 0) then 
         tmp(:,:) = -T_prog(n)%conversion*Grd%dyte(:,:)*fx_diff(:,:,n)
         call diagnose_2d(Time, Grd, id_sigma_diff_xflux_int_z(n), tmp(:,:))
     endif
     if(id_sigma_diff_yflux_int_z(n) > 0) then 
         tmp(:,:) = -T_prog(n)%conversion*Grd%dxtn(:,:)*fy_diff(:,:,n)
         call diagnose_2d(Time, Grd, id_sigma_diff_yflux_int_z(n), tmp(:,:))
     endif

     ! advective fluxes have minus sign built in
     if(id_sigma_adv_xflux_int_z(n) > 0) then 
         tmp(:,:) = T_prog(n)%conversion*Grd%dyte(:,:)*fx_adv(:,:,n)
         call diagnose_2d(Time, Grd, id_sigma_adv_xflux_int_z(n), tmp(:,:))
     endif
     if(id_sigma_adv_yflux_int_z(n) > 0) then 
         tmp(:,:) = T_prog(n)%conversion*Grd%dxtn(:,:)*fy_adv(:,:,n)
         call diagnose_2d(Time, Grd, id_sigma_adv_yflux_int_z(n), tmp(:,:))
     endif


  enddo ! end of num_prog_tracer loop


  ! tracer-independent diagnostics 

  if(id_diff_cet > 0) then 
      used = send_data (id_diff_cet, diff_cet(:,:),     &
             Time%model_time, rmask=tmask_sigma(:,:),   &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if(id_diff_cnt > 0) then 
      used = send_data (id_diff_cnt, diff_cnt(:,:),      &
             Time%model_time, rmask=tmask_sigma(:,:),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_thickness_sigma  > 0) then 
      used = send_data (id_thickness_sigma, thickness_sigma_eff(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),                &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  call diagnose_3d(Time, Grd, id_weight_sigma, wrk1(:,:,:))

  if (id_sigma_uhrho  > 0) then 
      used = send_data (id_sigma_uhrho, sigma_uhrho(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_sigma_vhrho  > 0) then 
      used = send_data (id_sigma_vhrho, sigma_vhrho(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_sigma_uhrho_sgs  > 0) then 
      used = send_data (id_sigma_uhrho_sgs, sigma_uhrho_sgs(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),            &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_sigma_vhrho_sgs  > 0) then 
      used = send_data (id_sigma_vhrho_sgs, sigma_vhrho_sgs(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),            &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_sigma_uhrho_res  > 0) then 
      used = send_data (id_sigma_uhrho_res, sigma_uhrho_res(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),            &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_sigma_vhrho_res  > 0) then 
      used = send_data (id_sigma_vhrho_res, sigma_vhrho_res(:,:), &
             Time%model_time,  rmask=tmask_sigma(:,:),            &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  call watermass_diag(Time, Dens)  


end subroutine sigma_transport
! </SUBROUTINE> 


!#######################################################################
! <SUBROUTINE NAME="ocean_sigma_transport_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_sigma_transport_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_sigma_transport_end: module needs initialization')
  endif 
  call save_restart(Sigma_restart, time_stamp)

end subroutine ocean_sigma_transport_restart
! </SUBROUTINE> NAME="ocean_sigma_transport_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_sigma_transport_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_sigma_transport_end(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_sigma_transport_end: module needs initialization')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_sigma_transport_mod (ocean_sigma_transport_end): NO restart written.'
    call mpp_error(WARNING,&
    '==>Warning from ocean_sigma_transport_mod (ocean_sigma_transport_end): NO restart written.')
    return
  endif 

  call ocean_sigma_transport_restart

  write(stdoutunit,*) ' '
  write(stdoutunit,*) 'From ocean_sigma_transport_mod: ending chksums'
  call write_timestamp(Time%model_time)
  call write_chksum_2d('thickness_sigma', thickness_sigma(COMP)*Grd%tmask(COMP,1))


end subroutine ocean_sigma_transport_end
! </SUBROUTINE> NAME="ocean_sigma_transport_end"


!#######################################################################
! <SUBROUTINE NAME="advect_sigma_upwind">
!
! <DESCRIPTION>
! First order upwind to advect tracers in sigma layer. 
! </DESCRIPTION>
!
subroutine advect_sigma_upwind(T_prog)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  real    :: velocity, upos, uneg
  integer :: i, j, n

  fx_adv(:,:,:)   = 0.0
  fy_adv(:,:,:)   = 0.0
  tchg_adv(:,:,:) = 0.0

  do n=1,num_prog_tracers

     ! i-advective flux
     do j=jsc,jec
        do i=isc-1,iec
           velocity      = 0.5*sigma_uhrho(i,j)
           upos          = velocity + abs(velocity)
           uneg          = velocity - abs(velocity)
           fx_adv(i,j,n) = Grd%dyte(i,j)*(upos*tracer_sigma(i,j,n) + uneg*tracer_sigma(i+1,j,n)) &
                           *tmask_sigma(i,j)*tmask_sigma(i+1,j)  
        enddo
     enddo

     ! j-advective flux
     do j=jsc-1,jec
        do i=isc,iec
           velocity      = 0.5*sigma_vhrho(i,j)
           upos          = velocity + abs(velocity)
           uneg          = velocity - abs(velocity)
           fy_adv(i,j,n) = Grd%dxtn(i,j)*(upos*tracer_sigma(i,j,n) + uneg*tracer_sigma(i,j+1,n)) &
                           *tmask_sigma(i,j)*tmask_sigma(i,j+1)
        enddo
     enddo

     ! to handle redundancies across the bipolar fold 
     if (Grd%tripolar) then 
         call mpp_update_domains(fx_adv(:,:,n),fy_adv(:,:,n), Dom_flux%domain2d, &
              gridtype=CGRID_NE, complete=T_prog(n)%complete) 
     endif


  enddo   ! enddo for num_prog_tracers 

end subroutine advect_sigma_upwind
! </SUBROUTINE> NAME="advect_sigma_upwind"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag = .false. 

  ! diagnostics for neutral density and dianeutral velocity component 
  id_neut_rho_sigma = register_diag_field ('ocean_model', 'neut_rho_sigma', &
                      Grd%tracer_axes(1:3), Time%model_time,                &
                      'update of neutral density from sigma transport',     &
                      '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_sigma > 0) compute_watermass_diag = .true.

  id_wdian_rho_sigma = register_diag_field ('ocean_model', 'wdian_rho_sigma', &
                         Grd%tracer_axes(1:3), Time%model_time,               &
                         'dianeutral mass transport due to sigma transport',  &
                         'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_sigma > 0) compute_watermass_diag = .true.

  id_tform_rho_sigma = register_diag_field ('ocean_model', 'tform_rho_sigma',                       &
                         Grd%tracer_axes(1:3), Time%model_time,                                     &
                         'watermass transform due to sigma transport on levels (pre-layer binning)',&
                         'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_sigma > 0) compute_watermass_diag = .true.

  id_neut_rho_sigma_on_nrho = register_diag_field ('ocean_model',                         &
     'neut_rho_sigma_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                &
     'update of neutral density from sigma transport as binned to neutral density layers',&
                           '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_sigma_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_rho_sigma_on_nrho = register_diag_field ('ocean_model',                           &
     'wdian_rho_sigma_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
     'dianeutral mass transport due to sigma transport as binned to neutral density layers', &
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_sigma_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_rho_sigma_on_nrho = register_diag_field ('ocean_model',                     &
     'tform_rho_sigma_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
     'watermass transform due to sigma transport as binned to neutral density layers', &
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_sigma_on_nrho > 0) compute_watermass_diag = .true.

  id_eta_tend_sigma= -1          
  id_eta_tend_sigma= register_diag_field ('ocean_model','eta_tend_sigma',&
       Grd%tracer_axes(1:2), Time%model_time,                            &
       'non-Bouss steric sea level tendency from sigma tendency', 'm/s', &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_sigma > 0) compute_watermass_diag=.true.

  id_eta_tend_sigma_glob= -1          
  id_eta_tend_sigma_glob= register_diag_field ('ocean_model', 'eta_tend_sigma_glob',&
       Time%model_time,                                                             &
       'global mean non-bouss steric sea level tendency from sigma tendency',       &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_sigma_glob > 0) compute_watermass_diag=.true.



  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_sigma_transport_mod w/ compute_watermass_diag=.true. to compute some watermass diagnostics.'  
  endif 


end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from sigma transport on the watermass transformation.
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: i,j,k,tau
  real, dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_sigma_transport (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return 
 
  tau=Time%tau

  ! send diagnostics for update of neutral density from sigma transport 
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  wrk6(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2(i,j,k) = wrk1(i,j,k)*(                                                            &
                Dens%drhodT(i,j,k)                                                                &
                *(tchg_diff(i,j,index_temp)+tchg_adv(i,j,index_temp)+tchg_smooth(i,j,index_temp)) &
                +Dens%drhodS(i,j,k)                                                               &
                *(tchg_diff(i,j,index_salt)+tchg_adv(i,j,index_salt)+tchg_smooth(i,j,index_salt)) &
                )
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk4(i,j,k) = wrk3(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk5(i,j,k) = wrk2(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk6(i,j,k) =-wrk2(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend 
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_sigma, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_sigma, wrk4(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_sigma, wrk5(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_sigma_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_sigma_on_nrho, wrk4)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_sigma_on_nrho, wrk5)

  if(id_eta_tend_sigma > 0 .or. id_eta_tend_sigma_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk6(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_sigma, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_sigma_glob, eta_tend, cellarea_r)
  endif


end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_sigma_transport_mod
