module ocean_thickness_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Determine thickness of grid cells.  
!</OVERVIEW>
!
!<DESCRIPTION>
! This module determines the thickness of grid cells. 
! Thicknesses are generally time dependent and functions
! of the vertical coordinate.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM (2012)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_thickness_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.
!  </DATA> 
!  <DATA NAME="debug_this_module_detail" TYPE="logical">
!  For debugging pressure coordinate models. Lots of grid information
!  printed.  
!  </DATA> 
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="linear_free_surface" TYPE="logical">
!  For debugging, set the thickness of top cell in 
!  geopotential model to time independent values.
!  This option is needed if use the kappa_sort diagnostic. 
!  Default linear_free_surface=.false.
!  </DATA> 
!
!  <DATA NAME="thickness_method" TYPE="character">
!  To determine whether use energetic method or finite volume method 
!  to compute the thickness of a grid cell. Options are 
!  thickness_method=energetic or thickness_method=finitevolume.
!  There is little overall difference in results for pbot and eta.
!  However, it has been found that for realistic bottom topography
!  simulations, the vertical velocity component is very noisy with 
!  the finitevolume approach.  So this approach is considered
!  experimental.  The default is thickness_method='energetic'.  
!  </DATA> 
!  <DATA NAME="full_step_topography" TYPE="logical">
!  For case where with to only have the dzt be determined by the full step
!  bottom topography. This nml option is provided only for backwards
!  compatibility with older mom experiments using the full step topog.  
!  </DATA> 
!  <DATA NAME="enforce_positive_dzt" TYPE="logical">
!  For cases where wish to run model even with negative thickness.
!  Default enforce_positive_dzt=.false. 
!  </DATA> 
!  <DATA NAME="epsilon_init_thickness" TYPE="dimensionless">
!  For determining how strict we are to check for the thickness of 
!  a column when initializing pressure based vertical coordinate models.
!  </DATA> 
!
!  <DATA NAME="thickness_dzt_min" UNITS="m" TYPE="real">
!  Minimum dzt set when enforce_positive_dzt set true.
!  Default thickness_dzt_min=1.0.  
!  </DATA> 
!
!  <DATA NAME="initialize_zero_eta" TYPE="logical">
!  For pressure-based models, we can (with some work) initialize the 
!  model to have a zero surface height.  The recomended approach is to 
!  allow the surface height to be whatever it wants to be, and let adjustments
!  smooth it over time.  Default initialize_zero_eta=.false. 
!  </DATA> 
!  <DATA NAME="thickness_dzt_min_init" UNITS="m" TYPE="real">
!  For determining a modified bottom depth array 
!  that is required to ensure pressure model, based on initial
!  in-situ density, retains a nontrivial bottom cell thickness
!  in the case when initialize_zero_eta=0.0
!  Default thickness_dzt_min_init=5.0.  
!  </DATA> 
!  <DATA NAME="rescale_mass_to_get_ht_mod" TYPE="logical">
!  Expedient to allow for the computation of ht_mod. 
!  in the case when initialize_zero_eta=.true.  Here, 
!  we run the pressure based model with a rescaled mass that 
!  is sufficient to maintain non-negative dzt, at least for 
!  a short period.  This allows for one to run a day integration 
!  to produce ht_mod. rescale_mass_to_get_ht_mod=.true. will 
!  produce spurious results in general due to problems with the 
!  pressure gradient computation.  So it is not recommended for 
!  more than initial day or so. Default rescale_mass_to_get_ht_mod=.false.  
!  </DATA> 
!
!  <DATA NAME="read_rescale_rho0_mask" TYPE="logical">
!  For reading in a basin mask of use to re-define rho0 in 
!  isolated regions such as the Black Sea.  This is used for 
!  modifying the definition of the pressure or pstar levels 
!  during the initialization of the thicknesses dst.  This 
!  approach is appropriate in general, but has only been tested
!  when modify the pressure levels within a fully enclosed basin.
!  Default read_rescale_rho0_mask=.false.
!  </DATA> 
!  <DATA NAME="rescale_rho0_mask_gfdl" TYPE="logical">
!  For specifying the rescale_rho0_mask based on reading in 
!  the GFDL regional mask. Default rescale_rho0_mask_gfdl=.false.
!  </DATA> 
!  <DATA NAME="rescale_rho0_basin_label" TYPE="real">
!  For rescaling rho0 in a basin with a number rescale_rho0_basin_label.
!  For the Black Sea using GFDL basin masks in OM3,
!  rescale_rho0_basin_label=7.0. Default rescale_rho0_basin_label=-1.0
!  </DATA> 
!  <DATA NAME="rescale_rho0_value" TYPE="logical">
!  Fractional value for rescaling rho0 in the a region.  
!  Default rescale_rho0_value=1.0. 
!  </DATA> 
!
!  <DATA NAME="depth_min_for_sigma" UNITS="m"  TYPE="real">
!  For sigma coordinates, have minimum depth so that have layers defined 
!  globally.  Masks will zero out results over land, but for numerics
!  it is useful to compute everywhere. Default depth_min_for_sigma=0.01.
!  </DATA> 
!
!  <DATA NAME="read_rho0_profile" TYPE="logical">
!  To read in an initial rho0(z) profile to assist in defining the 
!  initial settings for the pressure increments dst, for use in 
!  setting the pressure-based vertical coordinate grids.  Ideally,
!  this profile is determined by the level averaged density in 
!  the initial conditions.  Note that it is essential to have 
!  rho0_profile have a sensible value at all depths even if there
!  is no water there, since there are places where we divide by 
!  rho0_profile in rock.  Also, be mindful that with denser water 
!  at depth, the pressure levels will be coarser at depth than if 
!  using the trivial density profile rho0(k)=rho0. 
!  This option is experimental, so it is recommended that user
!  maintain the default read_rho0_profile=.false.  
!  </DATA> 
!
!  <DATA NAME="pbot0_simple" TYPE="logical">
!  For testing purposes, have this option compute pbot0=g*rho0*ht
!  with rho0= constant.  Default pbot0_simple=.false.
!  </DATA> 
!
!  <DATA NAME="update_dzwu_k0" TYPE="logical">
!  A bug in certain versions of MOM4p1 was present, whereby  
!  Thickness%dzwu(i,j,k=0) was never updated, except for GEOPOTENTIAL
!  vertical coordinates.  This logical, whose default is 
!  update_dzwu_k0=.true., is provided for legacy purposes.
!  To test the older results, have update_dzwu_k0=.false. 
!  </DATA> 
!
!  <DATA NAME="max_num_bad_print" TYPE="integer">
!  Maximum bad grid cells printout for identifying problematic simulations.
!  Default max_num_bad_print=25.
!  </DATA> 
! </NAMELIST>
!
use constants_mod,     only: epsln, c2dbars
use diag_manager_mod,  only: register_diag_field, register_static_field
use fms_mod,           only: write_version_number, error_mesg, FATAL, WARNING
use fms_mod,           only: read_data
use fms_mod,           only: open_namelist_file, close_file, check_nml_error, file_exist
use fms_io_mod,        only: register_restart_field, save_restart, restore_state
use fms_io_mod,        only: restart_file_type, reset_field_pointer 
use mpp_domains_mod,   only: mpp_update_domains, mpp_global_min, mpp_global_max, domain2d
use mpp_mod,           only: input_nml_file, stdout, stdlog, mpp_error, mpp_max, mpp_min, mpp_pe

use ocean_domains_mod,    only: get_local_indices
use ocean_grids_mod,      only: update_boundaries
use ocean_parameters_mod, only: GEOPOTENTIAL, ZSTAR, ZSIGMA, PRESSURE, PSTAR, PSIGMA
use ocean_parameters_mod, only: PRESSURE_BASED, DEPTH_BASED, TERRAIN_FOLLOWING
use ocean_parameters_mod, only: TWO_LEVEL, THREE_LEVEL
use ocean_parameters_mod, only: ENERGETIC, FINITEVOLUME
use ocean_parameters_mod, only: missing_value, rho0, rho0r, grav
use ocean_tracer_util_mod,only: dzt_min_max
use ocean_types_mod,      only: ocean_time_type, ocean_domain_type, ocean_external_mode_type
use ocean_types_mod,      only: ocean_grid_type, ocean_thickness_type, ocean_density_type
use ocean_types_mod,      only: ocean_time_steps_type, ocean_lagrangian_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_2d_u, diagnose_3d_u, diagnose_2d_en
use ocean_util_mod,       only: write_chksum_2d, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_2d

implicit none

private

#include <ocean_memory.h>

! for vertical coordinate choice 
integer :: vert_coordinate
integer :: vert_coordinate_class
integer :: vert_coordinate_type

! useful constant 
real :: grav_r 
real :: grav_rho0

! for remap operator
integer :: halo

! to write some output to all processors 
integer :: unit=6
 
! for diagnostics 
integer :: id_dst              = -1
integer :: id_dzt              = -1
integer :: id_dzten(2)         = -1
integer :: id_dzt_pbc          = -1
integer :: id_dztlo            = -1
integer :: id_dztup            = -1
integer :: id_dstlo            = -1
integer :: id_dstup            = -1
integer :: id_dzwt             = -1
integer :: id_dzt_dst          = -1
integer :: id_dzu              = -1
integer :: id_dzwu             = -1
integer :: id_rho_dzt          = -1
integer :: id_rho_dzten(2)     = -1
integer :: id_rho_dzu          = -1
integer :: id_rho_dzt_tendency = -1
integer :: id_pbot0            = -1
integer :: id_depth_st         = -1
integer :: id_depth_swt        = -1
integer :: id_depth_zt         = -1
integer :: id_geodepth_zt      = -1
integer :: id_geodepth_zwt     = -1
integer :: id_depth_zwt        = -1 
integer :: id_depth_zu         = -1
integer :: id_depth_zwu        = -1 
integer :: id_mass_u           = -1
integer :: id_mass_en(2)       = -1
integer :: id_mass_uE          = -1
integer :: id_mass_t           = -1
integer :: id_mass_tendency    = -1
integer :: id_thicku           = -1
integer :: id_thicken(2)       = -1
integer :: id_rescale_mass     = -1
integer :: id_rescale_rho0_mask= -1
integer :: id_ht_mod           = -1

integer :: id_dztL             = -1
integer :: id_dztE             = -1
integer :: id_dztloL           = -1
integer :: id_dztloE           = -1
integer :: id_dztupL           = -1
integer :: id_dztupE           = -1
integer :: id_dzwtL            = -1
integer :: id_dzwtE            = -1
integer :: id_dzuL             = -1
integer :: id_dzuE             = -1
integer :: id_dzwuL            = -1
integer :: id_dzwuE            = -1
integer :: id_rho_dztL         = -1
integer :: id_rho_dztE         = -1
integer :: id_rho_dzuL         = -1
integer :: id_rho_dzuE         = -1

logical :: used

logical :: use_blobs

real :: dtime

!for restart file
integer                       :: id_restart_rho_dzt  = 0
integer                       :: id_restart_rho_dztL = 0
integer                       :: id_restart_rho_dztT = 0
type(restart_file_type), save :: Thk_restart

! for discretization of time tendency ("threelevel" or "twolevel")
integer :: tendency    

! conversion for prinout
real :: convert_factor=1.0

! for setting dst with pressure-based vertical coordinate models 
real, allocatable, dimension(:) :: rho0_profile 

! for the basin mask used to rescale rho0
real, allocatable, dimension(:,:) :: rescale_rho0_mask 
real, allocatable, dimension(:,:) :: data


character(len=128) :: version=&
     '$Id: ocean_thickness.F90,v 20.0 2013/12/14 00:12:33 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

type(ocean_domain_type), pointer :: Dom =>NULL()

public ocean_thickness_init
public ocean_thickness_init_adjust
public ocean_thickness_end
public update_tcell_thick_blob
public update_E_thickness
public update_tcell_thickness
public update_ucell_thickness
public dzt_dst_update
public rho_dzt_tendency
public ocean_thickness_restart

private thickness_restart
private thickness_initialize 
private dst_land_adjust
private thickness_chksum 
private thickness_chksum_blobs
private thickness_details 
private REMAP_ZT_TO_ZU

logical :: module_is_initialized      = .false.
logical :: rescale_mass_to_get_ht_mod = .false.

character(len=32) :: thickness_method='energetic' 
logical :: linear_free_surface               = .false.
logical :: full_step_topography              = .false. 
logical :: enforce_positive_dzt              = .false.
logical :: debug_this_module                 = .false.
logical :: debug_this_module_detail          = .false.
logical :: write_a_restart                   = .true. 
logical :: read_rho0_profile                 = .false. 
logical :: pbot0_simple                      = .false.
logical :: initialize_zero_eta               = .false. 
logical :: read_rescale_rho0_mask            = .false.
logical :: rescale_rho0_mask_gfdl            = .false.
logical :: update_dzwu_k0                    = .true. 

integer :: max_num_bad_print        = 25
real    :: rescale_rho0_value       = 1.0
real    :: rescale_rho0_basin_label = -1.0
real    :: thickness_dzt_min_init   = 5.0     ! metre
real    :: thickness_dzt_min        = 1.0     ! metre
real    :: depth_min_for_sigma      = 0.01    ! metre
real    :: epsilon_init_thickness   = 1.e-5

namelist /ocean_thickness_nml/ debug_this_module, debug_this_module_detail, write_a_restart,     &
                               full_step_topography, initialize_zero_eta, enforce_positive_dzt,  &
                               depth_min_for_sigma, thickness_method, read_rho0_profile,         &
                               thickness_dzt_min, thickness_dzt_min_init,                        &
                               rescale_mass_to_get_ht_mod, pbot0_simple, epsilon_init_thickness, &
                               read_rescale_rho0_mask, rescale_rho0_mask_gfdl,                   &
                               rescale_rho0_basin_label, rescale_rho0_value, linear_free_surface,&
                               max_num_bad_print, update_dzwu_k0
contains


!#######################################################################
! <SUBROUTINE NAME="ocean_thickness_init">
!
! <DESCRIPTION>
! Initialize the thickness type. 
!
! For pressure-based vertical coordinates, this initialization here
! is preliminary. 
! </DESCRIPTION>
!
subroutine ocean_thickness_init  (Time, Time_steps, Domain, Grid, Ext_mode, Thickness, &
                                 ver_coordinate, ver_coordinate_class, ver_coordinate_type, &
                                 blobs, introduce_blobs, dtimein, debug)

  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_time_steps_type),    intent(in)           :: Time_steps 
  type(ocean_domain_type),        intent(in), target   :: Domain
  type(ocean_grid_type),          intent(inout)        :: Grid  
  type(ocean_external_mode_type), intent(in)           :: Ext_mode
  type(ocean_thickness_type),     intent(inout)        :: Thickness
  integer,                        intent(in)           :: ver_coordinate 
  integer,                        intent(in)           :: ver_coordinate_class
  integer,                        intent(in)           :: ver_coordinate_type
  logical,                        intent(in)           :: blobs
  logical,                        intent(in)           :: introduce_blobs
  real,                           intent(in)           :: dtimein
  logical,                        intent(in), optional :: debug

  integer :: ioun, io_status, ierr
  integer :: i,j,k
  logical :: error_flag=.false.

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then
    call mpp_error(FATAL, &
    '==>Error ocean_thickness_init: module has been initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_thickness_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_thickness_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_thickness_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_thickness_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_thickness_nml)  
  write (stdlogunit,ocean_thickness_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_thickness_mod with debug_this_module=.true.'  
  endif 
  if(debug_this_module_detail) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_thickness_mod with debug_this_module_detail=.true.'  
  endif 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_thickness with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif 

  grav_r    = 1.0/grav 
  grav_rho0 = grav*rho0
  dtime     = dtimein

  tendency              = Time_steps%tendency
  vert_coordinate       = ver_coordinate
  vert_coordinate_class = ver_coordinate_class
  vert_coordinate_type  = ver_coordinate_type

  if(vert_coordinate_class==DEPTH_BASED) then 
    convert_factor=1.0
  elseif(vert_coordinate_class==PRESSURE_BASED) then 
    convert_factor=c2dbars
  endif 

  use_blobs = blobs

  Dom  => Domain
  halo =  Dom%xhalo   ! (xhalo=yhalo assumed in MOM)

  if(enforce_positive_dzt) then 
    write(stdoutunit,'(/a)') '==>Warning: running with enforce_positive_dzt=.true. '
    write(stdoutunit,'(a)')  '            This option artifically truncates cell size to enforce dzt > 0.'
    write(stdoutunit,'(a,f10.4/)') '      Minimum dzt is set to ', thickness_dzt_min
  endif 

  if(vert_coordinate == ZSIGMA .or. vert_coordinate == PSIGMA) then
     write(stdoutunit,'(a,f12.4,a)') &
     '==>To remove division by zero with sigma-coord, add fictitious water of depth(m) ' &
     ,depth_min_for_sigma, ' over land. This water does not contribute to budgets.'
    if(Grid%tripolar) then 
       call mpp_error(WARNING, &
       '==>Warning ocean_thickness_mod: zsigma and psigma may have problems with tripolar. Testing incomplete.')
    endif 
  endif 

  if(vert_coordinate == GEOPOTENTIAL .and. linear_free_surface) then
     write(stdoutunit,'(a)') &
     '==>Running GEOPOTENTIAL model with linear_free_surface, so cell thicknesses are constant in time.' 
     call mpp_error(WARNING, &
       '==>Running GEOPOTENTIAL model with linear_free_surface, so cell thicknesses are constant in time.')
  endif 

  if(thickness_method=='energetic') then 
      Thickness%method = ENERGETIC
      write(stdoutunit,'(a)') &
      '==>Note: running ocean_thickness with thickness_method=energetic.'
  elseif(thickness_method=='finitevolume') then 
      Thickness%method = FINITEVOLUME
      write(stdoutunit,'(a)') &
      '==>Warning: ocean_thickness w/ thickness_method=finitevolume is experimental. Testing incomplete.'
       call mpp_error(WARNING, &
       '==>Warning ocean_thickness w/ thickness_method=finitevolume is experimental. Testing incomplete.')
  endif

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  allocate (Thickness%rho_dzt(isd:ied,jsd:jed,nk,3) )
  allocate (Thickness%rho_dzten(isd:ied,jsd:jed,nk,2) )
  allocate (Thickness%rho_dztr(isd:ied,jsd:jed,nk) )
  allocate (Thickness%rho_dzu(isd:ied,jsd:jed,nk,3) )
  allocate (Thickness%rho_dzur(isd:ied,jsd:jed,nk) )

  allocate (Thickness%rho_dzt_tendency(isd:ied,jsd:jed,nk) )

  allocate (Thickness%sea_lev(isd:ied,jsd:jed) )
  allocate (Thickness%dzt(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dzten(isd:ied,jsd:jed,nk,2) )
  allocate (Thickness%dzu(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dzwt(isd:ied,jsd:jed,0:nk) )
  allocate (Thickness%dzwu(isd:ied,jsd:jed,0:nk) )

  allocate (Thickness%dztlo(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dztup(isd:ied,jsd:jed,nk) )

  if (use_blobs) then
     allocate (Thickness%rho_dztL(isd:ied,jsd:jed,nk,3) )
     allocate (Thickness%rho_dztT(isd:ied,jsd:jed,nk,3) )
     
     allocate (Thickness%rho_dzuL(isd:ied,jsd:jed,nk,3) )
     allocate (Thickness%rho_dzuT(isd:ied,jsd:jed,nk,3) )

     allocate (Thickness%dztL(isd:ied,jsd:jed,nk) )
     allocate (Thickness%dztT(isd:ied,jsd:jed,nk,3) )
     
     allocate (Thickness%dzuL(isd:ied,jsd:jed,nk) )
     allocate (Thickness%dzuT(isd:ied,jsd:jed,nk) )
     
     allocate (Thickness%dzwtL(isd:ied,jsd:jed,0:nk) )
     allocate (Thickness%dzwtT(isd:ied,jsd:jed,0:nk) )
     
     allocate (Thickness%dzwuL(isd:ied,jsd:jed,0:nk) )
     allocate (Thickness%dzwuT(isd:ied,jsd:jed,0:nk) )
     
     allocate (Thickness%dztloL(isd:ied,jsd:jed,nk) )
     allocate (Thickness%dztloT(isd:ied,jsd:jed,nk) )
     
     allocate (Thickness%dztupL(isd:ied,jsd:jed,nk) )
     allocate (Thickness%dztupT(isd:ied,jsd:jed,nk) )
     
     allocate (Thickness%mass_uT(isd:ied,jsd:jed,3) )

     allocate (Thickness%blob_source(isd:ied,jsd:jed) )

     allocate (Thickness%geodepth_zwt(isd:ied,jsd:jed,nk) )
  endif

  allocate (Thickness%geodepth_zt(isd:ied,jsd:jed,nk) )
  allocate (Thickness%depth_zt(isd:ied,jsd:jed,nk) )
  allocate (Thickness%depth_zwt(isd:ied,jsd:jed,nk) )
  allocate (Thickness%depth_zu(isd:ied,jsd:jed,nk) )
  allocate (Thickness%depth_zwu(isd:ied,jsd:jed,nk) )

  allocate (Thickness%depth_st(isd:ied,jsd:jed,nk) )
  allocate (Thickness%depth_swt(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dst(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dswt(isd:ied,jsd:jed,0:nk) )
  allocate (Thickness%dzt_dst(isd:ied,jsd:jed,nk) )

  allocate (Thickness%dstlo(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dstup(isd:ied,jsd:jed,nk) )

  allocate (Thickness%pbot0(isd:ied,jsd:jed))
  allocate (Thickness%pbot0r(isd:ied,jsd:jed))

  allocate (Thickness%mass_source(isd:ied,jsd:jed,nk) )
  allocate (Thickness%mass_u(isd:ied,jsd:jed,3) )
  allocate (Thickness%mass_en(isd:ied,jsd:jed,2) )
  allocate (Thickness%thicku(isd:ied,jsd:jed,3) )
  allocate (Thickness%thicken(isd:ied,jsd:jed,2) )

#endif

  ! set rescale_rho0_mask=rescale_rho0_value in those regions where we modify 
  ! the rho0 value by a fraction. Otherwise rho0 is unscaled. 
  allocate (rescale_rho0_mask(isd:ied,jsd:jed) )
  rescale_rho0_mask(:,:) = 1.0
  if(read_rescale_rho0_mask .and. vert_coordinate_class==PRESSURE_BASED) then 

      allocate (data(isd:ied,jsd:jed) )
      data(:,:) = 0.0
      call read_data('INPUT/basin_mask','basin_mask',data,Domain%domain2d)

      ! the following is specific to the mask used at GFDL
      if(rescale_rho0_mask_gfdl) then 
          write(stdoutunit,'(a,f12.6)') &
          '==>Note: running ocean_thickness with rho0 in selected basin modified to rho0(kg/m3) = ',&
          rho0*rescale_rho0_value     
          do j=jsc,jec
             do i=isc,iec
                if(data(i,j)==rescale_rho0_basin_label) then 
                    rescale_rho0_mask(i,j)=rescale_rho0_value
                else 
                    rescale_rho0_mask(i,j)=1.0
                endif
             enddo
          enddo
      endif
      call mpp_update_domains(rescale_rho0_mask(:,:), Dom%domain2d)
  endif

  id_rescale_rho0_mask= register_static_field ('ocean_model', 'rescale_rho0_mask',        &
   Grid%tracer_axes(1:2),'fraction of rho0 used for specifying pressure and pstar levels',&
   'dimensionles',missing_value=missing_value, range=(/-1.e1,1.e8/))
  call diagnose_2d(Time, Grid, id_rescale_rho0_mask, rescale_rho0_mask(:,:))


  ! rho0_profile for setting dst with pressure based vertical coordinate models 
  allocate(rho0_profile(nk))
  rho0_profile(:) = rho0
  if(read_rho0_profile) then 

      write(stdoutunit,'(a)') &
       '==>Warning: ocean_thickness_mod: read rho0 profile to set dst w/ pressure models.'
      write(stdoutunit,'(a)') &
      '             This option is experimental, and NOT generally supported.'
      if(vert_coordinate_class==DEPTH_BASED) then 
          write(stdoutunit,'(a)') &
               '   Since using DEPTH_BASED vertical coordinates, rho0_profile = rho0.'
      endif

      if(vert_coordinate_class==PRESSURE_BASED) then 
          call read_data('INPUT/rho0_profile.nc','rho0_profile', rho0_profile)
          write(stdoutunit,'(a)') 'rho0_profile is used to define pressure grid and pressure gradients.'
          do k=1,nk
             write(stdoutunit,'(a,i4,a,e22.12)') 'rho0_profile(',k,') = ',rho0_profile(k)
             if(rho0_profile(k) <= 0.0) then 
                 call mpp_error(FATAL, &
                 '==>ocean_thickness_mod: rho0_profile must be > 0.0, even in rock, since we divide by rho0_profile.')
             endif
          enddo
      endif

  endif

  if(pbot0_simple) then 
     write(stdoutunit,'(/a)') '==>Warning from ocean_thickness_mod: pbot0_simple=.true.'
  endif 

  ! initialize vertical grid arrays
  call thickness_initialize(Grid, Thickness) 

  ! modify time dependent vertical grid arrays based on restart values   
  call thickness_restart(Time, Grid, Ext_mode, Thickness, introduce_blobs)

  if(debug_this_module) then 
    write(stdoutunit,*) ' '
    write(stdoutunit,*) ' From ocean_thickness_mod: initial thickness checksum'
    call write_timestamp(Time%model_time)
    if (use_blobs) then
       call thickness_chksum_blobs(Time, Grid, Thickness)
    else
       call thickness_chksum(Time, Grid, Thickness)
    endif
  endif 

  ! check that ht is deep enough to maintain chosen minimum size cell. 
  error_flag=.false. 
  do j=jsc,jec
     do i=isc,iec
        if(Grid%ht(i,j) > 0.0 .and. Grid%ht(i,j) < thickness_dzt_min) then 
            write(unit,'(a,i4,a,i4,a,e22.12,a,e22.12)') &
            '==>Error: ocean_thickness_init: ht(',i+Dom%ioff,',',j+Dom%joff,') = ',Grid%ht(i,j), &
            'is less than the chosen setting for thickness_dzt_min = ',thickness_dzt_min
            error_flag=.true.
        endif
     enddo
  enddo
  if(error_flag) then
      call mpp_error(FATAL,&
      '==>ocean_thickness_init: ocean is not deep enough to support setting for thickness_dzt_min')
  endif

  ! check that ht is deep enough to maintain chosen minimum size initial bottom cell. 
  if(vert_coordinate_class==PRESSURE_BASED) then 
      error_flag=.false. 
      do j=jsc,jec
         do i=isc,iec
            if(Grid%ht(i,j) > 0.0 .and. Grid%ht(i,j) < thickness_dzt_min_init) then 
                write(unit,'(a,i4,a,i4,a,e22.12,a,e22.12)') &
                '==>Error: ocean_thickness_init: ht(',i+Dom%ioff,',',j+Dom%joff,') = ',Grid%ht(i,j), &
                'is less than the chosen setting for thickness_dzt_min_init = ',thickness_dzt_min_init
                error_flag=.true.
            endif
         enddo
      enddo
      if(error_flag) then
          call mpp_error(FATAL,&
          '==>ocean_thickness_init: ocean is not deep enough to support setting for thickness_dzt_min_init')
      endif
  endif

  ! ids for diagnostic manager 

  id_dzu = register_diag_field ('ocean_model', 'dzu', Grid%vel_axes_uv(1:3), Time%model_time, &
                                'u-cell thickness', 'm',                                      &
                                missing_value=missing_value, range=(/-1e1,1e5/))

  id_dzwu = register_diag_field ('ocean_model', 'dzwu', Grid%vel_axes_uv(1:3), Time%model_time,&
                                 'thickness between velocity points', 'm',                     &
                                 missing_value=missing_value, range=(/-1e1,1e5/))

  id_dzt = register_diag_field ('ocean_model', 'dzt', Grid%tracer_axes(1:3), Time%model_time, &
                                't-cell thickness', 'm',                                      &
                                missing_value=missing_value, range=(/-1e1,1e5/),              &
                                standard_name='cell_thickness')

  id_dzten(1) = register_diag_field ('ocean_model', 'dzteast', Grid%tracer_axes(1:3),    &
                                Time%model_time, 'thickness at east face of T-cell', 'm',&
                                missing_value=missing_value, range=(/-1e1,1e5/))

  id_dzten(2) = register_diag_field ('ocean_model', 'dztnorth', Grid%tracer_axes(1:3),   &
                                Time%model_time,'thickness at north face of T-cell', 'm',&
                                missing_value=missing_value, range=(/-1e1,1e5/))

  id_dzt_pbc = register_diag_field ('ocean_model', 'dzt_pbc', Grid%tracer_axes(1:2), Time%model_time,&
                                    'ocean bottom t-cell total thickness', 'm',                      &
                                    missing_value=missing_value, range=(/-1e1,1e5/))

  id_dztlo = register_diag_field ('ocean_model', 'dztlo', Grid%tracer_axes(1:3), Time%model_time, &
                                  'distance from t-cell center to t-cell bottom', 'm',            &
                                  missing_value=missing_value, range=(/-1e1,1e5/))

  id_dztup = register_diag_field ('ocean_model', 'dztup', Grid%tracer_axes(1:3), Time%model_time,&
                                  'distance from t-cell center to t-cell top', 'm',              &
                                  missing_value=missing_value, range=(/-1e1,1e5/))

  id_dstlo = register_diag_field ('ocean_model', 'dstlo', Grid%tracer_axes(1:3), Time%model_time, &
                                's-distance from t-cell center to t-cell bottom', 's-units',      &
                                missing_value=missing_value, range=(/-1e5,1e5/))

  id_dstup = register_diag_field ('ocean_model', 'dstup', Grid%tracer_axes(1:3), Time%model_time, &
                                's-distance from t-cell center to t-cell top', 's-units',         &
                                missing_value=missing_value, range=(/-1e5,1e5/))

  id_dzwt = register_diag_field ('ocean_model', 'dzwt', Grid%tracer_axes(1:3), Time%model_time,  &
                                'thickness between tracer points', 'm', &
                                missing_value=missing_value, range=(/-1e1,1e5/))

  id_dst = register_diag_field ('ocean_model', 'dst', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean s-cell thickness', 's-units',                          &
                                missing_value=missing_value, range=(/-1e5,1e5/))
 
  id_geodepth_zt = register_diag_field ('ocean_model', 'geodepth_zt', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell depth relative to z=0', 'm',                                    &
                                 missing_value=missing_value, range=(/-1e1,1.e8/))
  if (use_blobs) then
     id_geodepth_zwt = register_diag_field ('ocean_model', 'geodepth_zwt', Grid%tracer_axes(1:3), Time%model_time, &
                                            'ocean bottom of t-cell depth relative to z=0', 'm',                   &
                                            missing_value=missing_value, range=(/-1e1,1.e8/))                       

     id_dzuL = register_diag_field ('ocean_model', 'dzuL', Grid%vel_axes_uv(1:3), Time%model_time, &
                                    'L system contribution to u-cell thickness', 'm',              &
                                    missing_value=missing_value, range=(/-1e1,1e5/))                 
     id_dzuE = register_diag_field ('ocean_model', 'dzuE', Grid%vel_axes_uv(1:3), Time%model_time, &
                                    'E system contribution to  u-cell thickness', 'm',             &
                                    missing_value=missing_value, range=(/-1e1,1e5/))                 

     id_dzwuL = register_diag_field ('ocean_model', 'dzwuL', Grid%vel_axes_uv(1:3), Time%model_time,    &
                                     'L system contribution to thickness between velocity points', 'm', &
                                     missing_value=missing_value, range=(/-1e1,1e5/))                   
     id_dzwuE = register_diag_field ('ocean_model', 'dzwuE', Grid%vel_axes_uv(1:3), Time%model_time,    &
                                     'E system contribution to thickness between velocity points', 'm', &
                                     missing_value=missing_value, range=(/-1e1,1e5/))                   
     id_dztL = register_diag_field ('ocean_model', 'dztL', Grid%tracer_axes(1:3), Time%model_time, &
                                    'L system contribution to t-cell thickness', 'm',              &
                                    missing_value=missing_value, range=(/-1e1,1e5/))                
     id_dztE = register_diag_field ('ocean_model', 'dztE', Grid%tracer_axes(1:3), Time%model_time, &
                                    'E system contribution to t-cell thickness', 'm',              &
                                    missing_value=missing_value, range=(/-1e1,1e5/))                
     id_dztloL = register_diag_field ('ocean_model', 'dztloL', Grid%tracer_axes(1:3), Time%model_time,              &
                                      'L system contribution to distance from t-cell center to t-cell bottom', 'm', &
                                      missing_value=missing_value, range=(/-1e1,1e5/))               
     id_dztloE = register_diag_field ('ocean_model', 'dztloE', Grid%tracer_axes(1:3), Time%model_time,              &
                                      'E system contribution to distance from t-cell center to t-cell bottom', 'm', &
                                      missing_value=missing_value, range=(/-1e1,1e5/))                 
     id_dztupL = register_diag_field ('ocean_model', 'dztupL', Grid%tracer_axes(1:3), Time%model_time,           &
                                      'L system contribution to distance from t-cell center to t-cell top', 'm', &
                                      missing_value=missing_value, range=(/-1e1,1e5/))                  
     id_dztupE = register_diag_field ('ocean_model', 'dztupE', Grid%tracer_axes(1:3), Time%model_time,           &
                                      'E system contribution to distance from t-cell center to t-cell top', 'm', &
                                      missing_value=missing_value, range=(/-1e1,1e5/))                  
     id_dzwtL = register_diag_field ('ocean_model', 'dzwtL', Grid%tracer_axes(1:3), Time%model_time,  &
                                     'L system contribution to thickness between tracer points', 'm', &
                                     missing_value=missing_value, range=(/-1e1,1e5/))                 
     id_dzwtE = register_diag_field ('ocean_model', 'dzwtE', Grid%tracer_axes(1:3), Time%model_time,  &
                                     'E system contribution to thickness between tracer points', 'm', &
                                     missing_value=missing_value, range=(/-1e1,1e5/))                 
     id_rho_dzuL = register_diag_field ('ocean_model', 'rho_dzuL', Grid%vel_axes_uv(1:3), Time%model_time, &
                                        'L system contribution to u-cell rho*thickness', '(kg/m^3)*m',     &
                                        missing_value=missing_value, range=(/-1.e8,1.e8/))                  
     id_rho_dzuE = register_diag_field ('ocean_model', 'rho_dzuE', Grid%vel_axes_uv(1:3), Time%model_time, &
                                        'E system contribution to u-cell rho*thickness', '(kg/m^3)*m',     &
                                        missing_value=missing_value, range=(/-1.e8,1.e8/))                  
     id_rho_dztL = register_diag_field ('ocean_model', 'rho_dztL', Grid%tracer_axes(1:3), Time%model_time, &
                                        'L system contribution to t-cell rho*thickness', '(kg/m^3)*m',     &
                                        missing_value=missing_value, range=(/-1.e8,1.e8/))                  
     id_rho_dztE = register_diag_field ('ocean_model', 'rho_dztE', Grid%tracer_axes(1:3), Time%model_time, &
                                        'E system contribution to t-cell rho*thickness', '(kg/m^3)*m',     &
                                        missing_value=missing_value, range=(/-1.e8,1.e8/))                  
     id_mass_uE  = register_diag_field ('ocean_model', 'mass_uE', Grid%vel_axes_uv(1:2), Time%model_time,   &
                                       'E system contribution to mass per area on C-cell column', 'kg/m^2', &
                                       missing_value=missing_value, range=(/-1e1,1e8/))                    
  endif
  if(vert_coordinate == GEOPOTENTIAL) then 

    id_depth_zt = register_diag_field ('ocean_model', 'depth_zt', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell depth relative to z=0', 'm',                                &
                                 missing_value=missing_value, range=(/-1e1,1.e8/))

    id_depth_st = register_diag_field ('ocean_model', 'depth_st', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell depth relative to z=0', 'm',                                &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_swt = register_diag_field ('ocean_model', 'depth_swt', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell bottom relative to z=0', 'm',                                 &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_zwt = register_diag_field ('ocean_model', 'depth_zwt', Grid%tracer_axes(1:3), Time%model_time, &
                                'depth to t-cell bottom relative to z=0', 'm',                              &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_zu = register_diag_field ('ocean_model', 'depth_zu', Grid%vel_axes_uv(1:3), Time%model_time, &
                                'ocean u-cell depth relative to z=0', 'm',                                &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_zwu = register_diag_field ('ocean_model', 'depth_zwu', Grid%vel_axes_uv(1:3), Time%model_time, &
                                'depth to ocean u-cell bottom reative to z=0', 'm',                         &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

  else 

    id_depth_zt = register_diag_field ('ocean_model', 'depth_zt', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell depth relative to column top', 'm',                         &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_st = register_diag_field ('ocean_model', 'depth_st', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell s-depth', 's-units',                                        &
                                 missing_value=missing_value, range=(/-1e8,1.e8/))

    id_depth_swt = register_diag_field ('ocean_model', 'depth_swt', Grid%tracer_axes(1:3), Time%model_time, &
                                'ocean t-cell bottom s-depth', 's-units',                                   &
                                 missing_value=missing_value, range=(/-1e8,1.e8/))

    id_depth_zwt = register_diag_field ('ocean_model', 'depth_zwt', Grid%tracer_axes(1:3), Time%model_time, &
                                'depth to t-cell bottom relative to column top', 'm',                       &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_zu = register_diag_field ('ocean_model', 'depth_zu', Grid%vel_axes_uv(1:3), Time%model_time, &
                                'ocean u-cell depth relative to column top', 'm',                         &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

    id_depth_zwu = register_diag_field ('ocean_model', 'depth_zwu', Grid%vel_axes_uv(1:3), Time%model_time, &
                                'depth to ocean u-cell bottom reative to column top', 'm',                  &
                                 missing_value=missing_value, range=(/-1e1,1.e10/))

  endif 


  id_dzt_dst = register_diag_field ('ocean_model', 'dzt_dst', Grid%tracer_axes(1:3), Time%model_time, &
                                    'ocean t-cell specific thickness', 'm/s-coordinate',              &
                                    missing_value=missing_value, range=(/-1.e8,1.e8/))

  id_rho_dzu = register_diag_field ('ocean_model', 'rho_dzu', Grid%vel_axes_uv(1:3), Time%model_time, &
                                    'u-cell rho*thickness', '(kg/m^3)*m',                             &
                                    missing_value=missing_value, range=(/-1.e8,1.e8/))

  id_rho_dzt = register_diag_field ('ocean_model', 'rho_dzt', Grid%tracer_axes(1:3), Time%model_time, &
                                    't-cell rho*thickness', '(kg/m^3)*m',                             &
                                    missing_value=missing_value, range=(/-1.e8,1.e8/),                &
                                    standard_name='sea_water_mass_per_unit_area')

  id_rho_dzten(1) = register_diag_field ('ocean_model', 'rho_dzteast', Grid%tracer_axes(1:3),        &
                                    Time%model_time, 't-cell rho*thickness at east face of T-cell',  &
                                    '(kg/m^3)*m', missing_value=missing_value, range=(/-1.e8,1.e8/))

  id_rho_dzten(2) = register_diag_field ('ocean_model', 'rho_dztnorth', Grid%tracer_axes(1:3),      &
                                    Time%model_time, 't-cell rho*thickness at north face of T-cell',&
                                    '(kg/m^3)*m',missing_value=missing_value, range=(/-1.e8,1.e8/))

  id_mass_u  = register_diag_field ('ocean_model', 'mass_u', Grid%vel_axes_uv(1:2), Time%model_time,&
                                    'mass per area on C-cell column', 'kg/m^2',                     &
                                    missing_value=missing_value, range=(/-1e1,1e8/))  

  id_mass_en(1) = register_diag_field ('ocean_model', 'mass_en1', Grid%tracer_axes(1:2), Time%model_time,&
                                       'vertical sum rho_dzten(1)', 'kg/m^2',                            &
                                       missing_value=missing_value, range=(/-1e1,1e8/))  

  id_mass_en(2) = register_diag_field ('ocean_model', 'mass_en2', Grid%tracer_axes(1:2), Time%model_time,&
                                       'vertical sum rho_dzten(2)', 'kg/m^2',                            &
                                       missing_value=missing_value, range=(/-1e1,1e8/))  

  id_thicku  = register_diag_field ('ocean_model', 'thicku', Grid%vel_axes_uv(1:2), Time%model_time, &
                                    'thickness from z=eta_u to z=-H of U-cell column', 'm',          &
                                    missing_value=missing_value, range=(/-1e1,1e8/))  

  id_thicken(1) = register_diag_field ('ocean_model', 'thicken1', Grid%vel_axes_u(1:2), &
                                       Time%model_time, 'vertical sum of dzten(1)', 'm',&
                                       missing_value=missing_value, range=(/-1e1,1e8/))  

  id_thicken(2) = register_diag_field ('ocean_model', 'thicken2', Grid%vel_axes_v(1:2), &
                                       Time%model_time, 'vertical sum of dzten(2)', 'm',&
                                       missing_value=missing_value, range=(/-1e1,1e8/))  

  id_rho_dzt_tendency = register_diag_field ('ocean_model', 'rho_dzt_tendency', Grid%tracer_axes(1:3), &
                                             Time%model_time, 'rho_dzt_tendency', '(kg/m^3)*(m/s)',    &
                                             missing_value=missing_value, range=(/-1e8,1e8/))  

  id_pbot0  = register_static_field ('ocean_model', 'pbot0', Grid%tracer_axes(1:2), &
                                     'reference bottom pressure t-cells', 'dbar',   &
                                     missing_value=missing_value, range=(/-1.e8,1.e8/))

  if(vert_coordinate_class==DEPTH_BASED) then 
      id_mass_t = register_diag_field ('ocean_model', 'mass_t', Grid%tracer_axes(1:3),&
                  Time%model_time, 'ocean t-cell volume times rho0', 'kg',            &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
      id_mass_tendency = register_diag_field ('ocean_model', 'mass_tendency', Grid%tracer_axes(1:3),&
                  Time%model_time, 'tendency for ocean t-cell volume times rho0', 'kg/s',           &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  else 
      id_mass_t = register_diag_field ('ocean_model', 'mass_t', Grid%tracer_axes(1:3),&
                  Time%model_time, 'ocean t-cell mass', 'kg',                         &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
      id_mass_tendency = register_diag_field ('ocean_model', 'mass_tendency', Grid%tracer_axes(1:3),&
                  Time%model_time, 'tendency for ocean t-cell mass', 'kg/s',                        &
                  missing_value=missing_value, range=(/-1.e20,1.e20/))
  endif

end subroutine ocean_thickness_init
! </SUBROUTINE> NAME="ocean_thickness_init"


!#######################################################################
! <SUBROUTINE NAME="thickness_initialize">
!
! <DESCRIPTION>
! Initialize vertical thicknesses of grid cells. 
! For Boussinesq models, this code assumes the 
! surface heights eta_t and eta_u are zero.  
!
! The values here are relevant for the time=0 initialization
! of the model. Some time independent arrays are also set here,
! but they are over-written if there is a restart file.
!
! For pressure based vertical coordinates, the results here
! assume density = rho0_profile(k). The values are readjusted in 
! ocean_thickness_init_adjust after we have determined the initial 
! in situ density.  This readjustment may involve enforcing an 
! initial eta field of zero value, or the default which allows
! for eta to initially be nonzero. 
!
! </DESCRIPTION>
!
subroutine thickness_initialize(Grid, Thickness)

 type(ocean_grid_type),      intent(inout) :: Grid  
 type(ocean_thickness_type), intent(inout) :: Thickness

 real                              :: sigma_sign, density_tmp
 integer                           :: i,j,k,n,kb,kmin

 ! initialization to zero  
 Thickness%sea_lev(:,:) = 0.0
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%rho_dzt_tendency(i,j,k) = 0.0
          Thickness%mass_source(i,j,k)      = 0.0
       enddo
    enddo
 enddo

 ! vertical increments for non-terrain following coordinates 
 if(vert_coordinate_type /= TERRAIN_FOLLOWING) then 

     do j=jsd,jed
        do i=isd,ied 
           Thickness%dzwt(i,j,0) = Grid%dzw(0)
           Thickness%dzwu(i,j,0) = Grid%dzw(0)
        enddo
     enddo

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzt(i,j,k)       = Grid%dzt(k)
              Thickness%depth_zt(i,j,k)  = Grid%zt(k)
              Thickness%dzwt(i,j,k)      = Grid%dzw(k)
              Thickness%depth_zwt(i,j,k) = Grid%zw(k)
           enddo
        enddo
     enddo

     ! modifications for partial step representation of bottom topography 
     if(full_step_topography) then 
         do j=jsd,jed
            do i=isd,ied
               kb = Grid%kmt(i,j)
               if (kb > 1) then
                   Grid%ht(i,j) = Thickness%depth_zwt(i,j,kb)
               endif
            enddo
         enddo
     else 
         do j=jsd,jed
            do i=isd,ied
               kb = Grid%kmt(i,j)
               if (kb > 1) then
                   Thickness%depth_zwt(i,j,kb)= Grid%ht(i,j)
                   Thickness%dzt(i,j,kb)      = &
                        Thickness%depth_zwt(i,j,kb)   - Thickness%depth_zwt(i,j,kb-1)
                   Thickness%depth_zt(i,j,kb) = &
                        Thickness%depth_zwt(i,j,kb-1) + Grid%fracdz(kb,0)*Thickness%dzt(i,j,kb)
                   Thickness%dzwt(i,j,kb-1)   = &
                        Thickness%depth_zt(i,j,kb)    - Thickness%depth_zt(i,j,kb-1)
                   Thickness%dzwt(i,j,kb)     = &
                        Thickness%depth_zwt(i,j,kb)   - Thickness%depth_zt(i,j,kb) 
               endif
            enddo
         enddo
     endif

 ! vertical increments for terrain following coordinates
 else 

     ! -1 <= ZSIGMA <= 0, so dst and dsw are positive 
     !  0 <= PSIGMA <= 1, so dst and dsw are negative 
     ! Recall s-coordinate is dimensionless when use ZSIGMA or PSIGMA.
     if(vert_coordinate==ZSIGMA) then  
         sigma_sign =-1.0
     elseif(vert_coordinate==PSIGMA) then
         sigma_sign = 1.0
     endif

     ! set depth_st and depth_swt according to Grid arrays
     k=0
     do j=jsd,jed
        do i=isd,ied
           Thickness%dswt(i,j,k) = Grid%dsw(k)
        enddo
     enddo
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%depth_st(i,j,k)  = Grid%st(k)
              Thickness%depth_swt(i,j,k) = Grid%sw(k)
              Thickness%dst(i,j,k)       = Grid%dst(k)
              Thickness%dswt(i,j,k)      = Grid%dsw(k)
           enddo
        enddo
     enddo

     ! now get the z-increments 
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzt_dst(i,j,k) = -sigma_sign*max(Grid%ht(i,j),depth_min_for_sigma)
              Thickness%dzt(i,j,k)     = Thickness%dzt_dst(i,j,k)*Thickness%dst(i,j,k)
           enddo
        enddo
     enddo

     wrk1(:,:,:) = epsln
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              if(Grid%tmask(i,j,k) > 0.0) then 
                  wrk1(i,j,k) = 1.0/Thickness%dzt_dst(i,j,k)
               endif
            enddo
         enddo
      enddo
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
                Thickness%dzwt(i,j,k) = 2.0*Thickness%dswt(i,j,k)/(wrk1(i,j,k)+wrk1(i,j,k+1))
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
             Thickness%dzwt(i,j,0)  = Thickness%dswt(i,j,0) *Thickness%dzt_dst(i,j,1)
             Thickness%dzwt(i,j,nk) = Thickness%dswt(i,j,nk)*Thickness%dzt_dst(i,j,nk)
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwu(i,j,0) = Thickness%dzwt(i,j,0)
         enddo
      enddo

      ! vertically sum the z-increments to get z-depths 
      k=1
      do j=jsd,jed
         do i=isd,ied
            Thickness%depth_zt(i,j,k)  = Thickness%dzwt(i,j,k-1)
            Thickness%depth_zwt(i,j,k) = Thickness%dzt(i,j,k)
         enddo
      enddo
      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%depth_zt(i,j,k)  = Thickness%depth_zt(i,j,k-1)  + Thickness%dzwt(i,j,k-1)
               Thickness%depth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k-1) + Thickness%dzt(i,j,k)
            enddo
         enddo
      enddo

  endif  ! endif for ZSIGMA and PSIGMA 


  ! half-thicknesses of T-cells 
  k=1
  do j=jsd,jed
     do i=isd,ied
        Thickness%dztlo(i,j,k) = Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k)
        Thickness%dztup(i,j,k) = Thickness%dzwt(i,j,k-1)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%dztlo(i,j,k) = Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k)
           Thickness%dztup(i,j,k) = Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)
        enddo
     enddo
  enddo

  ! U-cell thicknesses as minimum of surrounding T-cell thicknesses.
  ! This approach follows that used in all versions of MOM using 
  ! thicknesses that are functions of (i,j,k).  
  ! Also define thicknesses at east and north face of T-cell as 
  ! minimum of the two adjacent cells.  
  Thickness%dzu  = Thickness%dzt
  Thickness%dzwu = Thickness%dzwt
  Thickness%dzten(:,:,:,1) = Thickness%dzt(:,:,:)
  Thickness%dzten(:,:,:,2) = Thickness%dzt(:,:,:)
  do k=1,nk
     do j=jsd,jec
        do i=isd,iec
           Thickness%dzu(i,j,k)     = min(Thickness%dzt(i,j,k),   &
                                          Thickness%dzt(i+1,j,k), &
                                          Thickness%dzt(i,j+1,k), &
                                          Thickness%dzt(i+1,j+1,k))
           Thickness%dzwu(i,j,k)    = min(Thickness%dzwt(i,j,k),   &
                                          Thickness%dzwt(i+1,j,k), &
                                          Thickness%dzwt(i,j+1,k), &
                                          Thickness%dzwt(i+1,j+1,k))
           Thickness%dzten(i,j,k,1) = min(Thickness%dzt(i,j,k),    &
                                          Thickness%dzt(i+1,j,k))
           Thickness%dzten(i,j,k,2) = min(Thickness%dzt(i,j,k),    &
                                          Thickness%dzt(i,j+1,k))
        enddo
     enddo
  enddo
  if(update_dzwu_k0) then 
     k=0
     do j=jsd,jec
        do i=isd,iec
           Thickness%dzwu(i,j,k) = min(Thickness%dzwt(i,j,k),   &
                                       Thickness%dzwt(i+1,j,k), &
                                       Thickness%dzwt(i,j+1,k), &
                                       Thickness%dzwt(i+1,j+1,k))
        enddo
     enddo
  endif


 ! update across domain boundaries (dzten(1) same update as dxte; dzten(2) same update as dxtn)
 do k=1,nk
    call update_boundaries (Dom, Grid, Thickness%dzu(:,:,k)    , -1,-1,0,0) 
    call update_boundaries (Dom, Grid, Thickness%dzten(:,:,k,1), -1, 0,0,0) 
    call update_boundaries (Dom, Grid, Thickness%dzten(:,:,k,2),  0,-1,0,0) 
 enddo
 do k=0,nk
    call update_boundaries (Dom, Grid, Thickness%dzwu(:,:,k), -1,-1,0,0)
 enddo


 ! U-cell depths based on U-cell thicknesses 
 Thickness%depth_zu  = 0.0
 Thickness%depth_zwu = 0.0
 do j=jsd,jed
    do i=isd,ied
       Thickness%depth_zu(i,j,1)  = Thickness%dzwu(i,j,0)
       Thickness%depth_zwu(i,j,1) = Thickness%dzu(i,j,1)
    enddo
 enddo
 do k=2,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%depth_zu(i,j,k)  = Thickness%depth_zu(i,j,k-1)  + Thickness%dzwu(i,j,k-1)
          Thickness%depth_zwu(i,j,k) = Thickness%depth_zwu(i,j,k-1) + Thickness%dzu(i,j,k) 
       enddo
    enddo
 enddo

 ! total thickness of U-cell column 
 do n=1,3
    do j=jsd,jed
       do i=isd,ied
          k=Grid%kmu(i,j)
          if(k > 1) then 
              Thickness%thicku(i,j,n) = Thickness%depth_zwu(i,j,k) 
          else
              Thickness%thicku(i,j,n) = 0.0
          endif
       enddo
    enddo
 enddo

 ! total thickness of dzten columns 
 Thickness%thicken(:,:,:) = 0.0
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%thicken(i,j,1) = Thickness%thicken(i,j,1) + Thickness%dzten(i,j,k,1)
          Thickness%thicken(i,j,2) = Thickness%thicken(i,j,2) + Thickness%dzten(i,j,k,2)
       enddo
    enddo 
 enddo

 ! compute rho*dz values with rho0_profile*rescale_rho0_mask for initialization
 Thickness%mass_u(:,:,:) = 0.0
 do n=1,3
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             Thickness%rho_dzt(i,j,k,n) = rho0_profile(k)*rescale_rho0_mask(i,j)*Thickness%dzt(i,j,k)
             Thickness%rho_dzu(i,j,k,n) = rho0_profile(k)*rescale_rho0_mask(i,j)*Thickness%dzu(i,j,k)
             Thickness%mass_u(i,j,n)    = Thickness%mass_u(i,j,n)  &
                                          + Grid%umask(i,j,k)*Thickness%rho_dzu(i,j,k,n)
          enddo
       enddo
    enddo
 enddo
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%rho_dztr(i,j,k)    = 1.0/(Thickness%rho_dzt(i,j,k,1)+epsln)
          Thickness%rho_dzur(i,j,k)    = 1.0/(Thickness%rho_dzu(i,j,k,1)+epsln)
          Thickness%rho_dzten(i,j,k,1) = rho0_profile(k)*rescale_rho0_mask(i,j)*Thickness%dzten(i,j,k,1)
          Thickness%rho_dzten(i,j,k,2) = rho0_profile(k)*rescale_rho0_mask(i,j)*Thickness%dzten(i,j,k,2)
       enddo
    enddo
 enddo

 ! total mass per area of rho_dzten columns 
 Thickness%mass_en(:,:,:) = 0.0
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%mass_en(i,j,1) = Thickness%mass_en(i,j,1) + Thickness%rho_dzten(i,j,k,1)
          Thickness%mass_en(i,j,2) = Thickness%mass_en(i,j,2) + Thickness%rho_dzten(i,j,k,2)
       enddo
    enddo 
 enddo


 ! arrays for vertical s-coordinates; 
 ! ZSIGMA and PSIGMA already computed these arrays

 if(vert_coordinate==GEOPOTENTIAL .or. vert_coordinate==ZSTAR) then 
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzt_dst(i,j,k)   = 1.0
              Thickness%dst(i,j,k)       = Thickness%dzt(i,j,k)
              Thickness%depth_st(i,j,k)  = Thickness%depth_zt(i,j,k)
              Thickness%depth_swt(i,j,k) = Thickness%depth_zwt(i,j,k)
           enddo
        enddo
     enddo
     do k=0,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%dswt(i,j,k) = Thickness%dzwt(i,j,k)
           enddo
        enddo
     enddo

 elseif(vert_coordinate==PRESSURE .or. vert_coordinate==PSTAR) then       
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              density_tmp              = rho0_profile(k)*rescale_rho0_mask(i,j)
              Thickness%dzt_dst(i,j,k) = -grav_r/(density_tmp+epsln)
              Thickness%dst(i,j,k)     = -grav*density_tmp*Thickness%dzt(i,j,k)
           enddo
        enddo
     enddo
     do k=0,nk
        kmin = max(1,k)
        do j=jsd,jed
           do i=isd,ied
              Thickness%dswt(i,j,k) = -grav*rho0_profile(kmin)*rescale_rho0_mask(i,j)*Thickness%dzwt(i,j,k)
           enddo
        enddo
     enddo

     ! vertically sum the s-increments to get s-depths
     k=1
     do j=jsd,jed
        do i=isd,ied
           Thickness%depth_st(i,j,k)  = -Thickness%dswt(i,j,k-1)
           Thickness%depth_swt(i,j,k) = -Thickness%dst(i,j,k)
        enddo
     enddo
     do k=2,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%depth_st(i,j,k)  = Thickness%depth_st(i,j,k-1)  - Thickness%dswt(i,j,k-1)
              Thickness%depth_swt(i,j,k) = Thickness%depth_swt(i,j,k-1) - Thickness%dst(i,j,k)
           enddo
        enddo
     enddo

 endif

 ! half-thicknesses of T-cells in s-space 
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%dstlo(i,j,k) = Thickness%dztlo(i,j,k)/Thickness%dzt_dst(i,j,k)
          Thickness%dstup(i,j,k) = Thickness%dztup(i,j,k)/Thickness%dzt_dst(i,j,k)
       enddo
    enddo
 enddo

 ! initialize depth from z=0 to t-point.  
 ! this depth is used for computing the geopotential
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          Thickness%geodepth_zt(i,j,k)  = Thickness%depth_zt(i,j,k)
       enddo
    enddo
 enddo

 if (use_blobs) then
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             Thickness%geodepth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k)
          enddo
       enddo
    enddo
 endif

 ! tentative setting for reference bottom pressure (Pa) 
 ! on T-cells, as well as its inverse. Will recompute
 ! this field in ocean_thickness_init_adjust
 Thickness%pbot0(:,:)  = 0.0
 Thickness%pbot0r(:,:) = 0.0
 if(pbot0_simple) then 
     do j=jsd,jed
        do i=isd,ied
           Thickness%pbot0(i,j) = Grid%tmask(i,j,1)*grav*rho0*Grid%ht(i,j)
        enddo
     enddo
 else 
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%pbot0(i,j) = Thickness%pbot0(i,j) + Grid%tmask(i,j,k)*grav*Thickness%rho_dzt(i,j,k,1)
           enddo
        enddo
     enddo
 endif
 do j=jsd,jed
    do i=isd,ied
       if(Grid%kmt(i,j) > 1) then 
           Thickness%pbot0r(i,j) = 1.0/Thickness%pbot0(i,j)
       endif
    enddo
 enddo
 

 if (use_blobs) then
    ! To begin with, we assume a zero thickness for the Lagrangian system.
    ! Therefore, the Eulerian system occupies the full thicknesses.  This
    ! may change later in the initialisation, depending on whether we start
    ! from a restart with a non-zero thickness L system.
    Thickness%dztT(:,:,:,1)     = Thickness%dzt(:,:,:)      
    Thickness%dztT(:,:,:,2)     = Thickness%dzt(:,:,:)      
    Thickness%dztT(:,:,:,3)     = Thickness%dzt(:,:,:)      
    Thickness%dzuT(:,:,:)       = Thickness%dzu(:,:,:)      
    Thickness%dzwtT(:,:,:)      = Thickness%dzwt(:,:,:)     
    Thickness%dzwuT(:,:,:)      = Thickness%dzwu(:,:,:)     
    Thickness%dztloT(:,:,:)     = Thickness%dztlo(:,:,:)    
    Thickness%dztupT(:,:,:)     = Thickness%dztup(:,:,:)    
    Thickness%rho_dztT(:,:,:,:) = Thickness%rho_dzt(:,:,:,:)
    Thickness%rho_dzuT(:,:,:,:) = Thickness%rho_dzu(:,:,:,:)
    Thickness%mass_uT(:,:,:)    = Thickness%mass_u(:,:,:)   
    
    Thickness%dztL(:,:,:)       = 0.0
    Thickness%dzuL(:,:,:)       = 0.0
    Thickness%dzwtL(:,:,:)      = 0.0
    Thickness%dzwuL(:,:,:)      = 0.0
    Thickness%dztloL(:,:,:)     = 0.0
    Thickness%dztupL(:,:,:)     = 0.0
    Thickness%rho_dzuL(:,:,:,:) = 0.0
    Thickness%rho_dztL(:,:,:,:) = 0.0
 endif
 
end subroutine thickness_initialize
! </SUBROUTINE> NAME="thickness_initialize"


!#######################################################################
! <SUBROUTINE NAME="ocean_thickness_init_adjust">
!
! <DESCRIPTION>
!
! When initializing the pressure model, we can choose to initialize
! with a zero surface height, which requires some adjustment of 
! bottom depths, or allow the surface height to be nonzero.
!
! The default is a nonzero surface height via 
! initialize_zero_eta=.false.
!
!===============================================================
! The following comments refer to the possible issues arising 
! from enforcing eta=0 on initialiation via 
! initialize_zero_eta=.true.
!
! If we wish to enforce eta=0 on initialization, then we determine
! the depth of ocean needed to ensure a minimum bottom cell 
! thickness as set by nml paramter thickness_dzt_min_init. 
! This bottom depth is a function of the initial in-situ density 
! as well as the nml thickness_dzt_min_init.
!
! If Grid%ht is sufficient according to the specifications of 
! thickness_dzt_min_init and the initial in-situ density, then 
! model is likely to remain stable, in that there is very small 
! chance of losing bottom cell. 
! (assuming thickness_dzt_min_init is > 5.0 or so).
!
! If Grid%ht is too shallow, again given the initial rho and 
! pressure increments dst, then it is recommended that the user
! modify either the initial density (not easy) or the bottom 
! topography (easy). The modifications are often quite trivial,
! depending on thickness_dzt_min_init and how light certain regions 
! of the model are. The array "ht_mod" can be saved in the 
! diag_table and compared to Grid%ht to see what modifications 
! are required. If modificatations are to be made to the original
! grid_spec.nc file, simply use the NCO commands as follows:
!
! 1/ remove the original depth_t array from the original grid spec
! file, and write out a new grid spec file absent the depth_t array. 
! ncea -x -v depth_t grid_spec_old.nc grid_spec_new.nc
!
! 2/ append to the new grid spec file the modified bottom topography
! array ht_mod, which was generated in an earlier run of the model.
! ncea -A -v ht_mod history_file_with_ht_mod.nc grid_spec_new.nc
!
! 3/ rename ht_mod to depth_t to conform to standard MOM grid_spec
! name convention.
! ncvarrename new_grid.nc ht_mod depth_t 
!
! If there are no modifications required to the depth_t array, 
! then it is unlikely that the pressure model will evolve to the 
! point of loosing the bottom cell and thus the model will remain 
! stable. This conclusion is certainly based on the degree to which
! the model is run with changes in water masses.
!
!===============================================================
!
! NOTE:  We generally compute values for the grid increments over 
! land.  The land values that are in the global halos are not 
! going to be available from the restart files.  Hence, they 
! will need to be recomputed when restarting the model.  This 
! recalculation is done in dst_land_adjust.
!
! The land values are needed for use in computing U-cell grid 
! factors using remapping operators.  If we naively place zeros
! in the land, then, for example, rho_dzu next to land will be 
! wrong, as it will be based on averaging a non-zero interior
! rho_dzt value with a rho_dzt incorretly set to zero, and so 
! lead to spuriously small rho_dzu values next to land.  Note that 
! the placement of nonzero values of the grid values over land 
! is something that is trivially done for the z-models. More care
! is needed with pressure models.  
! 
! </DESCRIPTION>
!
subroutine ocean_thickness_init_adjust(Grid, Time, Dens, Ext_mode, Thickness)

  type(ocean_grid_type),          intent(inout)  :: Grid  
  type(ocean_time_type),          intent(in)     :: Time
  type(ocean_density_type),       intent(in)     :: Dens 
  type(ocean_external_mode_type), intent(inout)  :: Ext_mode
  type(ocean_thickness_type),     intent(inout)  :: Thickness

  logical :: warning_flag=.false.
  integer :: i, j, k, kb, m, n
  integer :: tau, taup1
  real    :: thick_tmp, density_tmp, delta_tmp
  real    :: max_fraction_differ, rescale_mass_min
  real, dimension(isd:ied,jsd:jed) :: fraction_differ 
  real, dimension(isd:ied,jsd:jed) :: rescale_mass
  real, dimension(isd:ied,jsd:jed) :: ht_mod

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (file_exist('INPUT/ocean_thickness.res.nc'))  then 
      if (id_pbot0 > 0) then 
         call diagnose_2d(Time, Grid, id_pbot0, c2dbars*Thickness%pbot0(:,:))
      endif
      return 
  endif

  tau   = Time%tau 
  taup1 = Time%taup1 

  ! for diagnostic purposes, compute initial bottom pressure 
  if(vert_coordinate_class==DEPTH_BASED .or. vert_coordinate_type==TERRAIN_FOLLOWING) then 

      if(.not. pbot0_simple) then 
          Thickness%pbot0(:,:)  = 0.0
          Thickness%pbot0r(:,:) = 0.0
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   Thickness%pbot0(i,j) = Thickness%pbot0(i,j) &
                        + Grid%tmask(i,j,k)*grav*Dens%rho(i,j,k,tau)*Thickness%dzt(i,j,k)
                enddo
             enddo
          enddo
          do j=jsd,jed
             do i=isd,ied
                if(Grid%kmt(i,j) > 1) then 
                    Thickness%pbot0r(i,j) = 1.0/Thickness%pbot0(i,j)
                endif
             enddo
          enddo
      endif

      if (id_pbot0 > 0) then 
         call diagnose_2d(Time, Grid, id_pbot0, c2dbars*Thickness%pbot0(:,:))
      endif

      return
  endif

  write(stdoutunit,'(/a/)') &
  '==>Note: ocean_thickness_init_adjust adjusting time=0 vgrid using in situ density.'

  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  ! save original dzt computed from thickness_initialize 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = Thickness%dzt(i,j,k)
           wrk3(i,j,k) = Thickness%dztlo(i,j,k)
           wrk4(i,j,k) = Thickness%dztup(i,j,k)
        enddo
     enddo
  enddo

  ! fill arrays everywhere
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Grid%tmask(i,j,k) > 0.0) then 
               density_tmp = Dens%rho(i,j,k,tau) 
           else 
               density_tmp = rho0_profile(k)*rescale_rho0_mask(i,j)
           endif
           Thickness%dzt_dst(i,j,k)   = -grav_r/(density_tmp+epsln)
           Thickness%dzt(i,j,k)       =  Thickness%dst(i,j,k)*Thickness%dzt_dst(i,j,k)
           Thickness%rho_dzt(i,j,k,:) =  density_tmp*Thickness%dzt(i,j,k)
           Thickness%rho_dztr(i,j,k)  =  1.0/(Thickness%rho_dzt(i,j,k,tau)+epsln)
           wrk1(i,j,k)                =  1.0/Thickness%dzt_dst(i,j,k)  
        enddo
     enddo
  enddo

  rescale_mass(:,:)    = 1.0 
  fraction_differ(:,:) = 0.0
  ht_mod(:,:)          = Grid%ht(:,:)


  ! code for enforcing eta=0 at initialization 
  if(initialize_zero_eta) then 

      ! check to see if column is too thick to add
      ! a bottom cell of thickness thickness_dzt_min_init,
      ! assuming the surface height is initially eta_t=0.  
      do j=jsd,jed
         do i=isd,ied
            kb=Grid%kmt(i,j) 
            if(kb>1 .and. Grid%ht(i,j) > thickness_dzt_min_init) then 
                thick_tmp = thickness_dzt_min_init
                do k=1,kb-1
                   thick_tmp = thick_tmp + Thickness%dzt(i,j,k) 
                enddo
                if(thick_tmp >= (1.0+epsilon_init_thickness)*Grid%ht(i,j)) then 
                    write(*,*)' ' 
                    write(*,'(a,i4,a,i4,a,f16.8,a,f16.8)')       &
                    '==>ocean_thickness_init_adjust: thick_tmp(' &
                    ,i+Dom%ioff,',',j+Dom%joff,') = ',thick_tmp, ' is >= Grid%ht = ',Grid%ht(i,j)
                    ht_mod(i,j)          = thick_tmp
                    rescale_mass(i,j)    = Grid%ht(i,j)/thick_tmp
                    fraction_differ(i,j) = (thick_tmp-Grid%ht(i,j))/Grid%ht(i,j)
                endif
            endif
         enddo
      enddo

      max_fraction_differ = mpp_global_max(Dom%domain2d,fraction_differ)
      write(stdoutunit,'(/a)') &
      '==>ocean_thickness_init_adjust: WARNING ABOUT INITIAL TOPOGRAPHY AND DENSITY AND VERTICAL LEVELS'
      write(stdoutunit,'(a)')  &
      'During initialization, sum(dzt[k=1,kmt-1]) resulted in columns that are too thick to fit a bottom cell with'
      write(stdoutunit,'(a,e24.12)') &
      'thickness equal to or greater than the namelist thickness_dzt_min_init = ',thickness_dzt_min_init 
      write(stdoutunit,'(a)') &
      'This situation may be of no consequence to the stability of the pressure-based model simulation, if '
      write(stdoutunit,'(a,e24.12,a)') &
      'the maximum fractional difference ',max_fraction_differ,' is small, and thickness_dzt_min_init is greater than roughly 5m.'
      write(stdoutunit,'(a)') &
      'If these conditions are NOT met, there is the distinct possibility of losing bottom cells should density' 
      write(stdoutunit,'(a)') &
      'evolve away from the initial conditions. A means for ensuring added model stability is to use the modified '
      write(stdoutunit,'(a)') &
      'bottom depths found by adding "ht_mod" to diag_table. These depths should replace "depth_t" in the original'
      write(stdoutunit,'(a)') &
      'grid specification file.  Use the NCO command "ncea" to extract the original depth_t and replace with ht_mod.'

      if(Time%init .and. rescale_mass_to_get_ht_mod) then 
          rescale_mass_min = mpp_global_min(Dom%domain2d,rescale_mass)
          if(rescale_mass_min < 1.0) then 
              write(stdoutunit,'(a)') &
                   '==>Warning: rescale_mass_to_get_ht_mod=.true. Pressure gradients are wrong. Set false to run long.'
              write(stdoutunit,'(a,e22.12)') &
                   '==>Warning: rescale_mass_min = ',rescale_mass_min 
              call mpp_error(WARNING, &
                   '==>ocean_thickness_init_adjust: rescale_mass_to_get_ht_mod=.true. Pressure gradients are wrong.')
          endif
      else
          rescale_mass_min = 1.0
      endif
      do k=1,nk
         rho0_profile(k) = rho0_profile(k)*rescale_mass_min 
      enddo

  ! endif for initialize_zero_eta
  endif


  ! refill grid arrays over full column.  reconsider bottom cell k=kmt later.  
  ! NOTE: this step just reproduces what was done above the initialize_zero_eta
  ! if-block, in the case when rescale_mass_min=1.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           density_tmp                = rho0_profile(k)*rescale_rho0_mask(i,j)
           Thickness%dst(i,j,k)       = -grav*density_tmp*wrk2(i,j,k)
           Thickness%dstlo(i,j,k)     = -grav*density_tmp*wrk3(i,j,k)
           Thickness%dstup(i,j,k)     = -grav*density_tmp*wrk4(i,j,k)
           if(Grid%tmask(i,j,k) > 0.0) then 
               density_tmp = Dens%rho(i,j,k,tau) 
           endif
           Thickness%dzt_dst(i,j,k)   = -grav_r/(density_tmp+epsln)
           Thickness%dzt(i,j,k)       =  Thickness%dst(i,j,k)*Thickness%dzt_dst(i,j,k)
           Thickness%dztlo(i,j,k)     =  Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
           Thickness%dztup(i,j,k)     =  Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
           Thickness%rho_dzt(i,j,k,:) =  density_tmp*Thickness%dzt(i,j,k)
           Thickness%rho_dztr(i,j,k)  =  1.0/(Thickness%rho_dzt(i,j,k,tau)+epsln)
           wrk1(i,j,k)                =  1.0/Thickness%dzt_dst(i,j,k)  
        enddo
     enddo
  enddo

  ! Recompute depth_st and depth_swt by vertically summing s-increments. 
  ! Note that this readjustment still maintains constant depth_st and constant
  ! depth_swt when k=constant, except at the bottom cell, since wrk2, wrk3, and
  ! wrk4 are only functions of k-level except at the bottom.  We revisit the 
  ! bottom cell issue below. 
  k=1
  do j=jsd,jed
     do i=isd,ied
        Thickness%depth_st(i,j,k)  = -Thickness%dstup(i,j,k)
        Thickness%depth_swt(i,j,k) = -Thickness%dst(i,j,k)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%depth_st(i,j,k)  = Thickness%depth_swt(i,j,k-1) - Thickness%dstup(i,j,k)
           Thickness%depth_swt(i,j,k) = Thickness%depth_st(i,j,k)    - Thickness%dstlo(i,j,k)
        enddo
     enddo
  enddo

  ! modifications for partial step representation of bottom topography 
  do j=jsd,jed
     do i=isd,ied
        kb=Grid%kmt(i,j) 
        if(kb>1) then 
            thick_tmp = 0.0
            do k=1,kb-1
               thick_tmp = thick_tmp + Thickness%dzt(i,j,k) 
            enddo
            if(thick_tmp > Grid%ht(i,j) .and. initialize_zero_eta) then 
                warning_flag=.true.
                write(unit,'(/a,i4,a,i4,a,i4, a, e22.12,a,i4,a,i4,a,e22.12)')                           &
                '==>Warning: ocean_thickness_init_adjust. with kb = ', kb, ' thick_tmp(',i+Dom%ioff,',',&
                j+Dom%joff,')= ',thick_tmp, ' is greater than ht('  ,i+Dom%ioff,',',                    &
                j+Dom%joff,')= ',Grid%ht(i,j)
                if(debug_this_module_detail) then  
                    call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, i, j, kb, &
                         'From ocean_thickness_init_adjust', unit)
                endif
            endif
            k=kb
            delta_tmp = Grid%ht(i,j) - thick_tmp 
            if(delta_tmp >= thickness_dzt_min_init) then 
                Thickness%dzt(i,j,k) = delta_tmp  
            else 
                Thickness%dzt(i,j,k) = thickness_dzt_min_init
            endif

            Thickness%dst(i,j,k)       = -grav*Dens%rho(i,j,k,tau)*Thickness%dzt(i,j,k)
            Thickness%rho_dzt(i,j,k,:) = Dens%rho(i,j,k,tau)*Thickness%dzt(i,j,k)
            Thickness%rho_dztr(i,j,k)  =  1.0/(Thickness%rho_dzt(i,j,k,tau)+epsln)
        endif
     enddo
  enddo

  if(debug_this_module_detail) then  
     call dzt_min_max(Time, Thickness, &
     'From ocean_thickness_init_adjust, dzt_min_max information')
  endif  

  if(warning_flag) then 
     call mpp_error(WARNING, &
     '==>ocean_thickness_init_adjust: there potentially remains a problem with dzt.')
  endif 

  ! compute new reference bottom pressure 
  if(.not. pbot0_simple) then 
      Thickness%pbot0(:,:)  = 0.0
      Thickness%pbot0r(:,:) = 0.0
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%pbot0(i,j) = Thickness%pbot0(i,j) + Grid%tmask(i,j,k)*grav*Thickness%rho_dzt(i,j,k,1)
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            if(Grid%kmt(i,j) > 1) then 
                Thickness%pbot0r(i,j) = 1.0/Thickness%pbot0(i,j)
            endif
         enddo
      enddo
  endif

  ! adjust bottom cell s-dimensions (recall dst<0 and dswt<0 for pressure coordinates).
  ! do not alter the intermediate cell s-thicknesses, as these remain independent of i,j.
  do j=jsd,jed
     do i=isd,ied
        kb=Grid%kmt(i,j) 
        if(kb>1) then 

            Thickness%depth_swt(i,j,kb) = Thickness%pbot0(i,j)
            Thickness%depth_st(i,j,kb)  = Thickness%depth_swt(i,j,kb-1) - Grid%fracdz(kb,0)*Thickness%dst(i,j,kb)
            Thickness%dswt(i,j,kb-1)    = -(Thickness%depth_st(i,j,kb)  - Thickness%depth_st(i,j,kb-1)) 
            Thickness%dswt(i,j,kb)      = -(Thickness%depth_swt(i,j,kb) - Thickness%depth_st(i,j,kb)) 

            Thickness%dstup(i,j,kb) = -(Thickness%depth_st(i,j,kb) -Thickness%depth_swt(i,j,kb-1))
            Thickness%dstlo(i,j,kb) = -(Thickness%depth_swt(i,j,kb)-Thickness%depth_st(i,j,kb))

            Thickness%dztup(i,j,kb) = Thickness%dstup(i,j,kb)*Thickness%dzt_dst(i,j,kb)
            Thickness%dztlo(i,j,kb) = Thickness%dstlo(i,j,kb)*Thickness%dzt_dst(i,j,kb)

        endif
     enddo
  enddo


  ! distances between tracer points 
  if(Thickness%method==ENERGETIC) then 
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzwt(i,j,k) = 2.0*Thickness%dswt(i,j,k)/(wrk1(i,j,k)+wrk1(i,j,k+1)) 
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwt(i,j,0) = Thickness%dswt(i,j,0)*Thickness%dzt_dst(i,j,1)
            kb=Grid%kmt(i,j)
            if(kb > 0) then 
                Thickness%dzwt(i,j,kb) = Thickness%dswt(i,j,kb)*Thickness%dzt_dst(i,j,kb)
            endif
         enddo
      enddo
  elseif(Thickness%method==FINITEVOLUME) then 
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwt(i,j,0) = Thickness%dztup(i,j,1)
         enddo
      enddo
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzwt(i,j,k) = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k+1) 
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            kb=Grid%kmt(i,j)
            if(kb > 0) then 
                Thickness%dzwt(i,j,kb) = Thickness%dztlo(i,j,kb)
            endif
         enddo
      enddo
  endif

  ! adjust remaining z-dimensions throughout
  ! column through vertical sum of z-increments. 
  do j=jsd,jed
     do i=isd,ied
        Thickness%depth_zt(i,j,1)  = Thickness%dzwt(i,j,0)
        Thickness%depth_zwt(i,j,1) = Thickness%dzt(i,j,1)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%depth_zt(i,j,k)  = Thickness%depth_zt(i,j,k-1)  + Thickness%dzwt(i,j,k-1)
           Thickness%depth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k-1) + Thickness%dzt(i,j,k)
        enddo
     enddo
  enddo


  ! U-cell and east/north thicknesses computed as in update_ucell_thickness
  do k=1,nk
     do j=jsd,jec
        do i=isd,iec
           Thickness%dzu(i,j,k)         = min(Thickness%dzt(i,j,k),   &
                                              Thickness%dzt(i+1,j,k), &
                                              Thickness%dzt(i,j+1,k), &
                                              Thickness%dzt(i+1,j+1,k))
           Thickness%dzwu(i,j,k)        = min(Thickness%dzwt(i,j,k),   &
                                              Thickness%dzwt(i+1,j,k), &
                                              Thickness%dzwt(i,j+1,k), &
                                              Thickness%dzwt(i+1,j+1,k))
           Thickness%dzten(i,j,k,1)     = min(Thickness%dzt(i,j,k),   &
                                              Thickness%dzt(i+1,j,k))
           Thickness%dzten(i,j,k,2)     = min(Thickness%dzt(i,j,k),   &
                                              Thickness%dzt(i,j+1,k))
           Thickness%rho_dzten(i,j,k,1) = min(Thickness%rho_dzt(i,j,k,1),&
                                              Thickness%rho_dzt(i+1,j,k,1))
           Thickness%rho_dzten(i,j,k,2) = min(Thickness%rho_dzt(i,j,k,1),&
                                              Thickness%rho_dzt(i,j+1,k,1))
        enddo
     enddo
  enddo
  if(update_dzwu_k0) then 
     k=0
     do j=jsd,jec
        do i=isd,iec
           Thickness%dzwu(i,j,k) = min(Thickness%dzwt(i,j,k),   &
                                       Thickness%dzwt(i+1,j,k), &
                                       Thickness%dzwt(i,j+1,k), &
                                       Thickness%dzwt(i+1,j+1,k))
        enddo
     enddo
  endif
  do m=1,3
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              Thickness%rho_dzu(i,j,k,m) = min(Thickness%rho_dzt(i,j,k,m),   &
                                               Thickness%rho_dzt(i+1,j,k,m), &
                                               Thickness%rho_dzt(i,j+1,k,m), &
                                               Thickness%rho_dzt(i+1,j+1,k,m))
           enddo
        enddo
     enddo
  enddo
  call mpp_update_domains(Thickness%dzu(:,:,:),       Dom%domain2d)
  call mpp_update_domains(Thickness%dzwu(:,:,:),      Dom%domain2d)
  call mpp_update_domains(Thickness%rho_dzu(:,:,:,:), Dom%domain2d)

  call mpp_update_domains(Thickness%dzten(:,:,:,1),     Dom%domain2d)
  call mpp_update_domains(Thickness%dzten(:,:,:,2),     Dom%domain2d)
  call mpp_update_domains(Thickness%rho_dzten(:,:,:,1), Dom%domain2d)
  call mpp_update_domains(Thickness%rho_dzten(:,:,:,2), Dom%domain2d)


  do k=1,nk 
     do j=jsd,jed
        do i=isd,ied
           Thickness%rho_dzur(i,j,k) = Grid%umask(i,j,k)/(Thickness%rho_dzu(i,j,k,taup1)+epsln)
        enddo
     enddo
  enddo
  do j=jsd,jed
     do i=isd,ied
        Thickness%depth_zu(i,j,1)  = Thickness%dzwu(i,j,0)
        Thickness%depth_zwu(i,j,1) = Thickness%dzu(i,j,1)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%depth_zu(i,j,k)  = Thickness%depth_zu(i,j,k-1)  + Thickness%dzwu(i,j,k-1)
           Thickness%depth_zwu(i,j,k) = Thickness%depth_zwu(i,j,k-1) + Thickness%dzu(i,j,k) 
        enddo
     enddo
  enddo

  ! thickness of U-cell column 
  do n=1,3
     do j=jsd,jed
        do i=isd,ied
           k=Grid%kmu(i,j)
           if(k > 1) then 
               Thickness%thicku(i,j,n) = Thickness%depth_zwu(i,j,k) 
           else
               Thickness%thicku(i,j,n) = 0.0
           endif
        enddo
     enddo
  enddo

  ! total thickness of dzten columns 
  Thickness%thicken(:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%thicken(i,j,1) = Thickness%thicken(i,j,1) + Thickness%dzten(i,j,k,1)
           Thickness%thicken(i,j,2) = Thickness%thicken(i,j,2) + Thickness%dzten(i,j,k,2)
        enddo
     enddo 
  enddo


  ! mass per area of U-cell column
  Thickness%mass_u(:,:,:) = 0.0 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%mass_u(i,j,1) = Thickness%mass_u(i,j,1) + Grid%umask(i,j,k)*Thickness%rho_dzu(i,j,k,tau)
           Thickness%mass_u(i,j,2) = Thickness%mass_u(i,j,1)
           Thickness%mass_u(i,j,3) = Thickness%mass_u(i,j,1)
        enddo
     enddo
  enddo

  ! mass per area of rho_dzte column
  Thickness%mass_en(:,:,:) = 0.0 
  do n=1,2
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%mass_en(i,j,n) = Thickness%mass_en(i,j,n) + Grid%tmasken(i,j,k,n)*Thickness%rho_dzten(i,j,k,n)
           enddo
        enddo  
     enddo
  enddo


  ! reinitialize the bottom pressure used for external mode calculation 
  do m=1,3
     Ext_mode%pbot_t(:,:,m)     = Thickness%pbot0(:,:)
     Ext_mode%anompb(:,:,m)     = Thickness%pbot0(:,:) - rho0*grav*Grid%ht(:,:)
     Ext_mode%anompb_bar(:,:,m) = Thickness%pbot0(:,:) - rho0*grav*Grid%ht(:,:)
     Ext_mode%pbot_u(:,:,m)     = REMAP_ZT_TO_ZU(Thickness%pbot0(:,:),Grid)
  enddo

  if (use_blobs) then
     ! If we have gotten to this stage, then we are not starting from a restart.
     ! If we are not starting from a thickness restart, then the L system must 
     ! globally have zero thickness.  Thus, the total thickness is equal to the 
     ! E system thickness.

     ! Earlier in the initialisation, the L system properties were set to 0.
     ! The T system was also set equal to the E system, however, some of the
     ! E system properties may have changed during the initialisation -- see 
     ! section 7.3 of Griffies, (2009).  So, we set the Total thickness equal 
     ! to the E system thickness again.  We leave the L system as 0, as that 
     ! system will not have been altered.

     Thickness%dztT(:,:,:,1)     = Thickness%dzt(:,:,:)
     Thickness%dztT(:,:,:,2)     = Thickness%dzt(:,:,:)
     Thickness%dztT(:,:,:,3)     = Thickness%dzt(:,:,:)
     Thickness%dzwtT(:,:,:)      = Thickness%dzwt(:,:,:)
     Thickness%dztloT(:,:,:)     = Thickness%dztlo(:,:,:)
     Thickness%dztupT(:,:,:)     = Thickness%dztup(:,:,:)
     Thickness%rho_dztT(:,:,:,:) =  Thickness%rho_dzt(:,:,:,:)
     
     Thickness%dzuT(:,:,:)       = Thickness%dzu(:,:,:)
     THickness%rho_dzuT(:,:,:,:) = Thickness%rho_dzu(:,:,:,:)
     Thickness%dzwuT(:,:,:)      = Thickness%dzwu(:,:,:)
     Thickness%mass_uT(:,:,:)    = Thickness%mass_u(:,:,:)
  endif

  id_rescale_mass= register_static_field ('ocean_model', 'rescale_mass', Grid%tracer_axes(1:2),        &
                                          'rescale_mass factor as function of Grid%ht and initial rho',&
                                          'dimensionless', missing_value=missing_value, range=(/-1.e8,1.e8/))
  call diagnose_2d(Time, Grid, id_rescale_mass, rescale_mass(:,:))

  id_ht_mod = register_static_field ('ocean_model', 'ht_mod', Grid%tracer_axes(1:2),  &
                                     'modified bottom depth for pressure model', 'm', &
                                     missing_value=missing_value, range=(/-1.e1,1.e8/))
  call diagnose_2d(Time, Grid, id_ht_mod, ht_mod(:,:))

  if (id_pbot0 > 0) then 
     call diagnose_2d(Time, Grid, id_pbot0, c2dbars*Thickness%pbot0(:,:))
  endif

  if(debug_this_module) then 
    write(stdoutunit,*) ' '
    write(stdoutunit,*) ' From ocean_thickness_init_adjust: readjusted initial thickness checksum'
    call write_timestamp(Time%model_time)
    if (use_blobs) then
       call thickness_chksum_blobs(Time, Grid, Thickness)
    else
       call thickness_chksum(Time, Grid, Thickness)
    endif
  endif 

  if(debug_this_module_detail) then 
     call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, isc, jsc, nk, &
                            'From ocean_thickness_init_adjust', stdoutunit)
     call dzt_min_max(Time, Thickness,'From ocean_thickness_init_adjust, dzt_min_max information')
  endif

  if(debug_this_module_detail) then 
      write(unit,'(a)')''
      write(unit,'(a)')'---------Bottom values for grid thicknesses---------'
      do j=jsd,jed
         do i=isd,ied
            kb=Grid%kmt(i,j) 
            if(kb>1) then 
                if(Thickness%depth_swt(i,j,kb) > 0.0) then 
                    write(unit,'(a)')''
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'depth_swt(kb-1)  (',i+Dom%ioff,',',j+Dom%joff,',',kb-1,') = ',Thickness%depth_swt(i,j,kb-1)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'depth_swt(kb)    (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%depth_swt(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'depth_st(kb)     (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%depth_st(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dswt(kb)         (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dswt(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dst(kb)          (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dst(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dstup(kb)        (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dstup(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dstlo(kb)        (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dstlo(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dstlo+dstup(kb)  (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',&
                         Thickness%dstlo(i,j,kb)+Thickness%dstup(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dzten(kb,1)      (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dzten(i,j,kb,1)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dzten(kb,2)      (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dzten(i,j,kb,2)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dzt(kb)          (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dzt(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dztup(kb)        (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztup(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dztlo(kb)        (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztlo(i,j,kb)
                    write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &
                         'dztlo+dztup(kb)  (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',&
                         Thickness%dztlo(i,j,kb)+Thickness%dztup(i,j,kb)
                    if (use_blobs) then
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                   
                            'dztL(kb)         (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztL(i,j,kb)  
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                   
                            'dztupL(kb)       (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztupL(i,j,kb)
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                   
                            'dztloL(kb)       (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztloL(i,j,kb)
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                   
                            'dztloL+dztupL(kb)(',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',&                       
                            Thickness%dztloL(i,j,kb)+Thickness%dztupL(i,j,kb)                                      
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                    
                            'dztT(kb)         (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztT(i,j,kb,tau)
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                    
                            'dztupT(kb)       (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztupT(i,j,kb) 
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                    
                            'dztloT(kb)       (',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',Thickness%dztloT(i,j,kb) 
                       write(unit,'(a,i4,a,i4,a,i4,a,f24.8)')  &                                                    
                            'dztloT+dztupT(kb)(',i+Dom%ioff,',',j+Dom%joff,',',kb,')   = ',&                        
                            Thickness%dztloT(i,j,kb)+Thickness%dztupT(i,j,kb)
                    endif
                endif
            endif
         enddo
      enddo
  endif


end subroutine ocean_thickness_init_adjust
! </SUBROUTINE> NAME="ocean_thickness_init_adjust"


!#######################################################################
! <SUBROUTINE NAME="thickness_restart">
!
! <DESCRIPTION>
! Read basic elements of thickness derived type from restart, 
! then set remaining elements of the type.  
!
! Note that some array elements may be time independent.
! Their values have already been set in thickness_initialize
! or ocean_thickness_init_adjust, and they should not be 
! over-written here (in particular, they should not be 
! reset here to zero).
!
! To ensure reproducibility across restarts, this routine 
! uses the same logic as update_tcell_thickness and 
! update_ucell_thickness. This is particularly important
! for the terrain following calculations.
!
! Note that ocean_thickness_init is called prior to 
! ocean_operators_init.  Hence, we cannot use any 
! operators in this subroutine.  
!
! </DESCRIPTION>
!
subroutine thickness_restart (Time, Grid, Ext_mode, Thickness, introduce_blobs)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid  
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_thickness_type),     intent(inout) :: Thickness
  logical,                        intent(in)    :: introduce_blobs

  character*128 file_name
  integer :: tau, taum1, taup1, id_restart
  integer :: i,j,k,n,kb
  real    :: density_tmp

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  file_name = 'ocean_thickness.res.nc'
  call register_thickness_restart(Time, introduce_blobs, file_name, Thickness)
  call register_thickness_restart_blobs(Time, introduce_blobs, file_name, Thickness)
  
  if (.NOT. file_exist('INPUT/ocean_thickness.res.nc')) then
     if (.NOT.Time%init) then 
        call mpp_error(FATAL,&
             'Expecting file INPUT/ocean_thickness.res.nc to exist.&
             &This file was not found and Time%init=.false.')
     endif
     return
  endif
  
  call restore_state(Thk_restart)

  if(tendency==THREE_LEVEL) then 
      write (stdoutunit,'(a)') '  Reading THREE_LEVEL restart for rho_dzt from INPUT/ocean_thickness.res.nc'
      write (stdoutunit,'(a)') '  If not an initial condition, then expect two time records for rho_dzt.'
      call mpp_update_domains(Thickness%rho_dzt(:,:,:,tau), Dom%domain2d)
      call mpp_update_domains(Thickness%rho_dzt(:,:,:,taup1), Dom%domain2d)

  elseif(tendency==TWO_LEVEL) then 

      write (stdoutunit,'(a)') '  Reading TWO_LEVEL restart for rho_dzt from INPUT/ocean_thickness.res.nc'
      write (stdoutunit,'(a)') '  Expect only one time record for rho_dzt.'
 
      ! initialize to nonzero to eliminate NaNs when mask 
      ! out processors and then change layout upon a restart 
      if (use_blobs) then
         if (introduce_blobs) then
            ! Introducing blobs mid-run, we set the total thickness to that of the E system.
            ! The L system is assumed to be zero and has already been set in thickness_initialize
            Thickness%dztT(:,:,:,1)     = Thickness%dzt(:,:,:)      
            Thickness%dztT(:,:,:,2)     = Thickness%dzt(:,:,:)      
            Thickness%dztT(:,:,:,3)     = Thickness%dzt(:,:,:)      
            Thickness%dzuT(:,:,:)       = Thickness%dzu(:,:,:)      
            Thickness%dzwtT(:,:,:)      = Thickness%dzwt(:,:,:)     
            Thickness%dzwuT(:,:,:)      = Thickness%dzwu(:,:,:)     
            Thickness%dztloT(:,:,:)     = Thickness%dztlo(:,:,:)    
            Thickness%dztupT(:,:,:)     = Thickness%dztup(:,:,:)    
            Thickness%rho_dztT(:,:,:,:) = Thickness%rho_dzt(:,:,:,:)
            Thickness%rho_dzuT(:,:,:,:) = Thickness%rho_dzu(:,:,:,:)
            Thickness%mass_uT(:,:,:)    = Thickness%mass_u(:,:,:)   
         endif

         where (Thickness%rho_dztT(:,:,:,taup1)==0.) Thickness%rho_dztT(:,:,:,taup1) = rho0*thickness_dzt_min_init
         call mpp_update_domains(Thickness%rho_dztL(:,:,:,taup1), Dom%domain2d)
         call mpp_update_domains(Thickness%rho_dztT(:,:,:,taup1), Dom%domain2d)
         Thickness%rho_dztL(:,:,:,tau) = Thickness%rho_dztL(:,:,:,taup1)
         Thickness%rho_dztT(:,:,:,tau) = Thickness%rho_dztT(:,:,:,taup1)
         Thickness%rho_dzt(:,:,:,taup1) = Thickness%rho_dztT(:,:,:,taup1) - Thickness%rho_dztL(:,:,:,taup1)
         Thickness%rho_dzt(:,:,:,tau)   = Thickness%rho_dztT(:,:,:,tau)   - Thickness%rho_dztL(:,:,:,tau)
      else
         where (Thickness%rho_dzt(:,:,:,taup1)==0.) Thickness%rho_dzt(:,:,:,taup1) = rho0*thickness_dzt_min_init
      endif

      call mpp_update_domains(Thickness%rho_dzt(:,:,:,taup1),  Dom%domain2d)
      Thickness%rho_dzt(:,:,:,tau)  = Thickness%rho_dzt(:,:,:,taup1)
  endif 

  ! initialize to nonzero to eliminate NaNs when mask 
  ! out processors and then change layout upon a restart 
  where (Thickness%dzt_dst(:,:,:) ==0.) Thickness%dzt_dst(:,:,:) = thickness_dzt_min_init
  where (Thickness%dstlo(:,:,:)   ==0.) Thickness%dstlo(:,:,:)   = thickness_dzt_min_init
  where (Thickness%dstup(:,:,:)   ==0.) Thickness%dstup(:,:,:)   = thickness_dzt_min_init
  where (Thickness%dst(:,:,:)     ==0.) Thickness%dst(:,:,:)     = thickness_dzt_min_init
  where (Thickness%dswt(:,:,:)    ==0.) Thickness%dswt(:,:,:)    = thickness_dzt_min_init
  
  call mpp_update_domains(Thickness%dzt_dst(:,:,:), Dom%domain2d)
  call mpp_update_domains(Thickness%dstlo(:,:,:),   Dom%domain2d)
  call mpp_update_domains(Thickness%dstup(:,:,:),   Dom%domain2d)
  call mpp_update_domains(Thickness%dst(:,:,:),     Dom%domain2d)
  call mpp_update_domains(Thickness%dswt(:,:,:),    Dom%domain2d)
  call mpp_update_domains(Thickness%pbot0(:,:),     Dom%domain2d)
  if (use_blobs) then
     call mpp_update_domains(Thickness%dztloL(:,:,:),  Dom%domain2d)
     call mpp_update_domains(Thickness%dztupL(:,:,:),  Dom%domain2d)
  endif

  if(vert_coordinate_class==PRESSURE_BASED) then 
      call mpp_update_domains(Thickness%depth_st(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Thickness%depth_swt(:,:,:), Dom%domain2d)
  endif

  ! for adjusting values of pressure-coordinate intervals over land 
  call dst_land_adjust(Grid, Thickness)

  if (use_blobs) then
     do k=1,nk                                                                         
        do j=jsd,jed                                                                   
           do i=isd,ied                                                                
              Thickness%dztL(i,j,k) = Thickness%dztloL(i,j,k) + Thickness%dztupL(i,j,k)
           enddo
        enddo
     enddo
     
     do j=jsd,jed                                                                                  
        do i=isd,ied                                                                               
           Thickness%dzwtL(i,j,0)      = Thickness%dztupL(i,j,1)                                  
           Thickness%dzwtL(i,j,1:nk-1) = Thickness%dztupL(i,j,2:nk) + Thickness%dztloL(i,j,1:nk-1)
           Thickness%dzwtL(i,j,nk)     = Thickness%dztloL(i,j,nk)                                 
        enddo
     enddo
  endif

  ! inverse pbot0 values 
  Thickness%pbot0r = 0.0
  do j=jsd,jed
     do i=isd,ied
        if(Grid%kmt(i,j) > 1) then 
            Thickness%pbot0r(i,j) = 1.0/Thickness%pbot0(i,j)
        endif
     enddo
  enddo


  ! T-cell dz values 
  if(vert_coordinate==GEOPOTENTIAL) then 

     if (use_blobs) then
        do j=jsd,jed
           do i=isd,ied
              Thickness%dztloT(i,j,1)     = Thickness%dstlo(i,j,1)
              Thickness%dztupT(i,j,1)     = Thickness%dstup(i,j,1)
              Thickness%dztT(i,j,1,taup1) = Thickness%dztloT(i,j,1) + Thickness%dztupT(i,j,1)
           enddo
        enddo
        do k=1,nk                                                                                    
           do j=jsd,jed                                                                              
              do i=isd,ied                                                                           
                 Thickness%dztlo(i,j,k)    = Thickness%dztloT(i,j,k) - Thickness%dztloL(i,j,k)
                 Thickness%dztup(i,j,k)    = Thickness%dztupT(i,j,k) - Thickness%dztupL(i,j,k)
                 Thickness%dzt(i,j,k)      = Thickness%dztT(i,j,k,taup1)   - Thickness%dztL(i,j,k)
                 Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)        
                 wrk1(i,j,k)               = 1.0                                               
              enddo
           enddo
        enddo
    else
        do j=jsd,jed
           do i=isd,ied
              Thickness%dztlo(i,j,1) = Thickness%dstlo(i,j,1)
              Thickness%dztup(i,j,1) = Thickness%dstup(i,j,1)
              Thickness%dzt(i,j,1)   = Thickness%dztlo(i,j,1) + Thickness%dztup(i,j,1)
            Thickness%rho_dztr(i,j,1) = 1.0/(Thickness%rho_dzt(i,j,1,taup1)+epsln)
            wrk1(i,j,:)               = 1.0
           enddo
        enddo
    endif 

  elseif(vert_coordinate==ZSIGMA) then 
      ! note that it is important to do the "max" step here,
      ! as the restart will not get the halos correct, so we must do 
      ! this calculation here to get dzt_dst in the halos.
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzt_dst(i,j,k)  = max(depth_min_for_sigma,Grid%ht(i,j)+Ext_mode%eta_t(i,j,taup1))
               Thickness%dztlo(i,j,k)    = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztup(i,j,k)    = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dzt(i,j,k)      = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
               Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               wrk1(i,j,k)               = 1.0/Thickness%dzt_dst(i,j,k)
            enddo
         enddo
      enddo

  elseif(vert_coordinate==PRESSURE .or. vert_coordinate==PSTAR) then  
     if (use_blobs) then
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 if(Grid%tmask(i,j,k)==0.0) then  
                    density_tmp              = rho0_profile(k)*rescale_rho0_mask(i,j)
                    Thickness%dzt_dst(i,j,k) = -grav_r/(density_tmp+epsln)
                 endif
              enddo
           enddo
        enddo
        call calculate_restart_thickness(Thickness,taup1)
     else
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 if(Grid%tmask(i,j,k)==0.0) then  
                    density_tmp              = rho0_profile(k)*rescale_rho0_mask(i,j)
                    Thickness%dzt_dst(i,j,k) = -grav_r/(density_tmp+epsln)
                 endif
                 Thickness%dztlo(i,j,k)    = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
                 Thickness%dztup(i,j,k)    = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
                 Thickness%dzt(i,j,k)      = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
                 Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
                 wrk1(i,j,k)               = 1.0/Thickness%dzt_dst(i,j,k)
              enddo
           enddo
        enddo
     endif

  elseif(vert_coordinate==PSIGMA) then 
      ! note that it is important to do the Grid%tmask(i,j,k)==0.0 steps,
      ! as the restart will not get the halos correct, so we must do 
      ! this calculation here to get dzt_dst in the halos. 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               if(Grid%tmask(i,j,k) == 0.0) then  
                  Thickness%dzt_dst(i,j,k)  = -depth_min_for_sigma
               endif 
               Thickness%dztlo(i,j,k)    = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztup(i,j,k)    = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dzt(i,j,k)      = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
               if(Grid%tmask(i,j,k) == 0.0) then
                  Thickness%rho_dzt(i,j,k,taup1) = rho0_profile(k)*Thickness%dzt(i,j,k)
               endif 
               Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               wrk1(i,j,k)               = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo

  else !(vert_coordinate==ZSTAR)
     if (use_blobs) then
        call calculate_restart_thickness(Thickness,taup1)
     else
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dztlo(i,j,k)    = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
                 Thickness%dztup(i,j,k)    = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
                 Thickness%dzt(i,j,k)      = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
                 Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
                 wrk1(i,j,k)               = 1.0/Thickness%dzt_dst(i,j,k)
              enddo
           enddo
        enddo
     endif
  endif

  if (use_blobs) then
     ! vertical distance (m) between t-cell points 
     if(Thickness%method==ENERGETIC) then 
        do k=1,nk-1
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwtT(i,j,k) = 2.0*Thickness%dswt(i,j,k)/(wrk1(i,j,k)+wrk1(i,j,k+1)) 
                 Thickness%dzwt(i,j,k)  = Thickness%dzwtT(i,j,k) - Thickness%dzwtL(i,j,k)
              enddo
           enddo
        enddo
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwtT(i,j,0) = Thickness%dswt(i,j,0)*Thickness%dzt_dst(i,j,1)
              Thickness%dzwt(i,j,0)  = Thickness%dzwtT(i,j,0) - Thickness%dzwtL(i,j,0)
              kb=Grid%kmt(i,j)
              if(kb > 0) then 
                 Thickness%dzwtT(i,j,kb) = Thickness%dswt(i,j,kb)*Thickness%dzt_dst(i,j,kb)
                 Thickness%dzwt(i,j,kb)  = Thickness%dzwtT(i,j,kb) - Thickness%dzwtL(i,j,kb)
              endif
           enddo
        enddo
     elseif(Thickness%method==FINITEVOLUME) then 
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwtT(i,j,0) = Thickness%dztupT(i,j,1)
              Thickness%dzwt(i,j,0)  = Thickness%dzwtT(i,j,1) - Thickness%dzwtL(i,j,1)
           enddo
        enddo
        do k=1,nk-1
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwtT(i,j,k) = Thickness%dztloT(i,j,k) + Thickness%dztupT(i,j,k+1) 
                 Thickness%dzwt(i,j,k)  = Thickness%dzwtT(i,j,k)  - Thickness%dzwtL(i,j,k)
              enddo
           enddo
        enddo
        do j=jsd,jed
           do i=isd,ied
              kb=Grid%kmt(i,j)
              if(kb > 0) then 
                 Thickness%dzwtT(i,j,kb) = Thickness%dztloT(i,j,kb)
                 Thickness%dzwt(i,j,kb)  = Thickness%dzwtT(i,j,kb) - Thickness%dzwtL(i,j,kb)
              endif
           enddo
        enddo
     endif

  else !do not use blobs
     ! vertical distance (m) between t-cell points 
     if(Thickness%method==ENERGETIC) then 
        do k=1,nk-1
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwt(i,j,k) = 2.0*Thickness%dswt(i,j,k)/(wrk1(i,j,k)+wrk1(i,j,k+1)) 
              enddo
           enddo
        enddo
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwt(i,j,0) = Thickness%dswt(i,j,0)*Thickness%dzt_dst(i,j,1)
              kb=Grid%kmt(i,j)
              if(kb > 0) then 
                 Thickness%dzwt(i,j,kb) = Thickness%dswt(i,j,kb)*Thickness%dzt_dst(i,j,kb)
              endif
           enddo
        enddo
     elseif(Thickness%method==FINITEVOLUME) then 
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwt(i,j,0) = Thickness%dztup(i,j,1)
           enddo
        enddo
        do k=1,nk-1
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwt(i,j,k) = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k+1) 
              enddo
           enddo
        enddo
        do j=jsd,jed
           do i=isd,ied
              kb=Grid%kmt(i,j)
              if(kb > 0) then 
                 Thickness%dzwt(i,j,kb) = Thickness%dztlo(i,j,kb)
              endif
           enddo
        enddo
     endif
  endif !use_blobs


  ! vertical depth (m) of t and w points 
  if(vert_coordinate /= GEOPOTENTIAL) then 

     if (use_blobs) then
        ! depth from z=eta to t-point and t-cell bottom 
        do j=jsd,jed
           do i=isd,ied
              Thickness%depth_zt(i,j,1)  = Thickness%dzwtT(i,j,0)
              Thickness%depth_zwt(i,j,1) = Thickness%dztT(i,j,1,taup1)
           enddo
        enddo
        do k=2,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%depth_zt(i,j,k)  = Thickness%depth_zt(i,j,k-1)  + Thickness%dzwtT(i,j,k-1)
                 Thickness%depth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k-1) + Thickness%dztT(i,j,k,taup1)
              enddo
           enddo
        enddo
        
        ! depth from z=0 to t point
        do j=jsd,jed
           do i=isd,ied
              Thickness%geodepth_zt(i,j,1)  = Thickness%depth_zt(i,j,1)  - Ext_mode%eta_t(i,j,taup1)
              Thickness%geodepth_zwt(i,j,1) = Thickness%depth_zwt(i,j,1) - Ext_mode%eta_t(i,j,taup1)
           enddo
        enddo
        do k=2,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%geodepth_zt(i,j,k)  = Thickness%geodepth_zt(i,j,k-1)  + Thickness%dzwtT(i,j,k-1)
                 Thickness%geodepth_zwt(i,j,k) = Thickness%geodepth_zwt(i,j,k-1) + Thickness%dztT(i,j,k,taup1)
              enddo
           enddo
        enddo
     else !do not use_blobs
        
        ! depth from z=eta to t-point and t-cell bottom         
        do j=jsd,jed
           do i=isd,ied
              Thickness%depth_zt(i,j,1)  = Thickness%dzwt(i,j,0)
              Thickness%depth_zwt(i,j,1) = Thickness%dzt(i,j,1)
           enddo
        enddo
        do k=2,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%depth_zt(i,j,k)  = Thickness%depth_zt(i,j,k-1)  + Thickness%dzwt(i,j,k-1)
                 Thickness%depth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k-1) + Thickness%dzt(i,j,k)
              enddo
           enddo
        enddo
        
        ! depth from z=0 to t point
        do j=jsd,jed
           do i=isd,ied
              Thickness%geodepth_zt(i,j,1) = Thickness%depth_zt(i,j,1) - Ext_mode%eta_t(i,j,taup1)
           enddo
        enddo
        do k=2,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%geodepth_zt(i,j,k) = Thickness%geodepth_zt(i,j,k-1) + Thickness%dzwt(i,j,k-1)
              enddo
           enddo
        enddo
        
     endif!use blobs
  endif ! not geopotential coordinates

  if (use_blobs) then
     ! Calculate the L system thickness for the U grid, based on the L system thickness
     ! of the T grid.  These are required for all coordinate systems supported by the 
     ! Lagrangian blob framework.
     do j=jsd,jec
        do i=isd,iec
           Thickness%dzwuL(i,j,0) = min(Thickness%dzwtL(i,  j,  0),&
                                        Thickness%dzwtL(i+1,j,  0),&
                                        Thickness%dzwtL(i,  j+1,0),&
                                        Thickness%dzwtL(i+1,j+1,0)) 
        enddo
     enddo
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              Thickness%dzuL(i,j,k) = min(Thickness%dztL(i,  j,  k), &
                                          Thickness%dztL(i+1,j,  k), &
                                          Thickness%dztL(i,  j+1,k), &
                                          Thickness%dztL(i+1,j+1,k))

              Thickness%dzwuL(i,j,k) = min(Thickness%dzwtL(i,  j,  k), &
                                           Thickness%dzwtL(i+1,j,  k), &
                                           Thickness%dzwtL(i,  j+1,k), &
                                           Thickness%dzwtL(i+1,j+1,k))

              Thickness%rho_dzuL(i,j,k,taup1) = min(Thickness%rho_dztL(i,  j,  k,taup1), &
                                                    Thickness%rho_dztL(i+1,j,  k,taup1), &
                                                    Thickness%rho_dztL(i,  j+1,k,taup1), &
                                                    Thickness%rho_dztL(i+1,j+1,k,taup1))
           enddo
        enddo
     enddo
     call mpp_update_domains(Thickness%dzuL(:,:,:),           Dom%domain2d)
     call mpp_update_domains(Thickness%dzwuL(:,:,:),          Dom%domain2d)
     call mpp_update_domains(Thickness%rho_dzuL(:,:,:,taup1), Dom%domain2d)
  endif

  ! U-cell thicknesses computed as in update_ucell_thickness
  if(vert_coordinate==GEOPOTENTIAL) then 

     if (use_blobs) then
        do k=1,nk
           do j=jsd,jec
              do i=isd,iec
                 Thickness%rho_dzuL(i,j,k,tau) = min(Thickness%rho_dztL(i,  j,  k,tau), &
                                                     Thickness%rho_dztL(i+1,j,  k,tau), &
                                                     Thickness%rho_dztL(i,  j+1,k,tau), &
                                                     Thickness%rho_dztL(i+1,j+1,k,tau))
              enddo
           enddo
        enddo
        call mpp_update_domains(Thickness%rho_dzuL(:,:,:,tau), Dom%domain2d)
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwuT(i,j,0)          = Thickness%depth_zt(i,j,1)  + Ext_mode%eta_u(i,j,taup1)
              Thickness%dzwu(i,j,0)           = Thickness%dzwuT(i,j,0) - Thickness%dzwuL(i,j,0)
              Thickness%dzuT(i,j,1)           = Thickness%depth_zwt(i,j,1) + Ext_mode%eta_u(i,j,taup1)
              Thickness%rho_dzuT(i,j,1,taup1) = rho0*Thickness%dzuT(i,j,1)
              Thickness%rho_dzuT(i,j,1,tau)   = rho0*(Thickness%depth_zwt(i,j,1) + Ext_mode%eta_u(i,j,tau))
           enddo
        enddo
        if(linear_free_surface) then  
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwuT(i,j,0)          = Thickness%depth_zt(i,j,1)  
                 Thickness%dzwu(i,j,0)           = Thickness%dzwuT(i,j,0) - Thickness%dzwuL(i,j,0)
                 Thickness%dzuT(i,j,1)           = Thickness%depth_zwt(i,j,1) 
                 Thickness%rho_dzuT(i,j,1,taup1) = rho0*Thickness%dzuT(i,j,1)
                 Thickness%rho_dzuT(i,j,1,tau)   = rho0*Thickness%depth_zwt(i,j,1) 
              enddo
           enddo
        endif !linear_free_surface
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzu(i,j,k)           = Thickness%dzuT(i,j,k)           - Thickness%dzuL(i,j,k)
                 Thickness%rho_dzu(i,j,k,taup1) = Thickness%rho_dzuT(i,j,k,taup1) - Thickness%rho_dzuL(i,j,k,taup1)
                 Thickness%rho_dzu(i,j,k,tau)   = Thickness%rho_dzuT(i,j,k,tau)   - Thickness%rho_dzuL(i,j,k,tau)
                 Thickness%rho_dzur(i,j,k)      = 1.0/(Thickness%rho_dzu(i,j,k,taup1)+epsln)
              enddo
           enddo
        enddo

     else !do not use_blobs

        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwu(i,j,0)          = Thickness%depth_zt(i,j,1)  + Ext_mode%eta_u(i,j,taup1)
              Thickness%dzu(i,j,1)           = Thickness%depth_zwt(i,j,1) + Ext_mode%eta_u(i,j,taup1)
              Thickness%rho_dzu(i,j,1,taup1) = rho0*Thickness%dzu(i,j,1)
              Thickness%rho_dzu(i,j,1,tau)   = rho0*(Thickness%depth_zwt(i,j,1) + Ext_mode%eta_u(i,j,tau))
              Thickness%rho_dzur(i,j,1)      = 1.0/(Thickness%rho_dzu(i,j,1,taup1)+epsln)
           enddo
        enddo
        if(linear_free_surface) then  
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwu(i,j,0)          = Thickness%depth_zt(i,j,1)  
                 Thickness%dzu(i,j,1)           = Thickness%depth_zwt(i,j,1) 
                 Thickness%rho_dzu(i,j,1,taup1) = rho0*Thickness%dzu(i,j,1)
                 Thickness%rho_dzu(i,j,1,tau)   = rho0*Thickness%depth_zwt(i,j,1) 
                 Thickness%rho_dzur(i,j,1)      = 1.0/(Thickness%rho_dzu(i,j,1,taup1)+epsln)
              enddo
           enddo
        endif
     endif

  else ! not geopotential

      if(vert_coordinate==ZSIGMA .or. vert_coordinate==PSIGMA) then 
          ! note: zsigma and psigma not supported by this implementation of 
          ! the blobs framework
          do k=1,nk 
             do j=jsd,jed
                do i=isd,ied
                   Thickness%dzu(i,j,k)           = Thickness%dzt(i,j,k)
                   Thickness%dzten(i,j,k,1)       = Thickness%dzt(i,j,k)
                   Thickness%dzten(i,j,k,2)       = Thickness%dzt(i,j,k)
                   Thickness%dzwu(i,j,k)          = Thickness%dzwt(i,j,k)
                   Thickness%rho_dzu(i,j,k,taup1) = Thickness%rho_dzt(i,j,k,taup1)
                enddo
             enddo
          enddo
          k=0
          do j=jsd,jed
             do i=isd,ied
                Thickness%dzwu(i,j,k) = Thickness%dzwt(i,j,k)
             enddo
          enddo

      else !not zsima or psigma

         if (use_blobs) then

            call mpp_update_domains(Thickness%rho_dztT(:,:,:,taup1), Dom%domain2d)
            do j=jsd,jec
               do i=isd,iec
                  Thickness%dzwuT(i,j,0) = min(Thickness%dzwtT(i,  j,  0),&
                                               Thickness%dzwtT(i+1,j,  0),&
                                               Thickness%dzwtT(i,  j+1,0),&
                                               Thickness%dzwtT(i+1,j+1,0)) 
               enddo
            enddo
            do k=1,nk
               do j=jsd,jec
                  do i=isd,iec
                     Thickness%dzuT(i,j,k) = min(Thickness%dztT(i,  j,  k,taup1), &
                                                 Thickness%dztT(i+1,j,  k,taup1), &
                                                 Thickness%dztT(i,  j+1,k,taup1), &
                                                 Thickness%dztT(i+1,j+1,k,taup1))

                     Thickness%dzwuT(i,j,k) = min(Thickness%dzwtT(i,  j,  k), &
                                                  Thickness%dzwtT(i+1,j,  k), &
                                                  Thickness%dzwtT(i,  j+1,k), &
                                                  Thickness%dzwtT(i+1,j+1,k))
                     
                     Thickness%rho_dzuT(i,j,k,taup1) = min(Thickness%rho_dztT(i,  j,  k,taup1), &
                                                           Thickness%rho_dztT(i+1,j,  k,taup1), &
                                                           Thickness%rho_dztT(i,  j+1,k,taup1), &
                                                           Thickness%rho_dztT(i+1,j+1,k,taup1))
                  enddo
               enddo
            enddo
              
            call mpp_update_domains(Thickness%dzuT(:,:,:),           Dom%domain2d)
            call mpp_update_domains(Thickness%dzwuT(:,:,:),          Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzuT(:,:,:,taup1), Dom%domain2d)

            Thickness%rho_dzuT(:,:,:,tau) = Thickness%rho_dzuT(:,:,:,taup1)

            Thickness%dzu(:,:,:)  = Thickness%dzuT(:,:,:)  - Thickness%dzuL(:,:,:)
            Thickness%dzwu(:,:,:) = Thickness%dzwuT(:,:,:) - Thickness%dzwuL(:,:,:)
            Thickness%rho_dzu(:,:,:,taup1) = Thickness%rho_dzuT(:,:,:,taup1) - Thickness%rho_dzuL(:,:,:,taup1)

         else  !do not use_blobs

            do k=1,nk
               do j=jsd,jec
                  do i=isd,iec
                     Thickness%dzu(i,j,k)  = min(Thickness%dzt(i,j,k),   &
                                                 Thickness%dzt(i+1,j,k), &
                                                 Thickness%dzt(i,j+1,k), &
                                                 Thickness%dzt(i+1,j+1,k))
                     Thickness%dzwu(i,j,k) = min(Thickness%dzwt(i,j,k),   &
                                                 Thickness%dzwt(i+1,j,k), &
                                                 Thickness%dzwt(i,j+1,k), &
                                                 Thickness%dzwt(i+1,j+1,k))
                     Thickness%rho_dzu(i,j,k,taup1) = min(Thickness%rho_dzt(i,j,k,taup1),   &
                                                          Thickness%rho_dzt(i+1,j,k,taup1), &
                                                          Thickness%rho_dzt(i,j+1,k,taup1), &
                                                          Thickness%rho_dzt(i+1,j+1,k,taup1))
                     Thickness%dzten(i,j,k,1)     = min(Thickness%dzt(i,j,k), &
                                                        Thickness%dzt(i+1,j,k))
                     Thickness%dzten(i,j,k,2)     = min(Thickness%dzt(i,j,k), &
                                                        Thickness%dzt(i,j+1,k))
                     Thickness%rho_dzten(i,j,k,1) = min(Thickness%rho_dzt(i,j,k,taup1),&
                                                        Thickness%rho_dzt(i+1,j,k,taup1))
                     Thickness%rho_dzten(i,j,k,2) = min(Thickness%rho_dzt(i,j,k,taup1),&
                                                        Thickness%rho_dzt(i,j+1,k,taup1))

                  enddo
               enddo
            enddo
            if(update_dzwu_k0) then 
               k=0
               do j=jsd,jec
                  do i=isd,iec
                     Thickness%dzwu(i,j,k) = min(Thickness%dzwt(i,j,k),   &
                                                 Thickness%dzwt(i+1,j,k), &
                                                 Thickness%dzwt(i,j+1,k), &
                                                 Thickness%dzwt(i+1,j+1,k))
                  enddo
               enddo
            endif
            
            call mpp_update_domains(Thickness%dzu(:,:,:),           Dom%domain2d)
            call mpp_update_domains(Thickness%dzwu(:,:,:),          Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzu(:,:,:,taup1), Dom%domain2d)

            call mpp_update_domains(Thickness%dzten(:,:,:,1),      Dom%domain2d)
            call mpp_update_domains(Thickness%dzten(:,:,:,2),      Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzten(:,:,:,1),  Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzten(:,:,:,2),  Dom%domain2d)

            
         endif  !use_blobs
          
      endif! sigma coords


      ! this next bit done for all coords except geopotential

      do k=1,nk 
         do j=jsd,jed
            do i=isd,ied
               Thickness%rho_dzur(i,j,k) = 1.0/(Thickness%rho_dzu(i,j,k,taup1)+epsln)
            enddo
         enddo
      enddo

      if (use_blobs) then

         do j=jsd,jed
            do i=isd,ied
               Thickness%depth_zu(i,j,1)  = Thickness%dzwuT(i,j,0)
               Thickness%depth_zwu(i,j,1) = Thickness%dzuT(i,j,1)
            enddo
         enddo
         do k=2,nk
            do j=jsd,jed
               do i=isd,ied
                  Thickness%depth_zu(i,j,k)  = Thickness%depth_zu(i,j,k-1)  + Thickness%dzwuT(i,j,k-1)
                  Thickness%depth_zwu(i,j,k) = Thickness%depth_zwu(i,j,k-1) + Thickness%dzuT(i,j,k) 
               enddo
            enddo
         enddo

      else !do not use blobs

         do j=jsd,jed
            do i=isd,ied
               Thickness%depth_zu(i,j,1)  = Thickness%dzwu(i,j,0)
               Thickness%depth_zwu(i,j,1) = Thickness%dzu(i,j,1)
            enddo
         enddo
         do k=2,nk
            do j=jsd,jed
               do i=isd,ied
                  Thickness%depth_zu(i,j,k)  = Thickness%depth_zu(i,j,k-1)  + Thickness%dzwu(i,j,k-1)
                  Thickness%depth_zwu(i,j,k) = Thickness%depth_zwu(i,j,k-1) + Thickness%dzu(i,j,k) 
               enddo
            enddo
         enddo

      endif !use_blobs


   endif  !geopotential


  ! thickness of U-cell column 
  do n=1,3
     do j=jsd,jed
        do i=isd,ied
           k=Grid%kmu(i,j)
           if(k > 1) then 
               Thickness%thicku(i,j,n) = Thickness%depth_zwu(i,j,k) 
           else
               Thickness%thicku(i,j,n) = 0.0
           endif
        enddo
     enddo
  enddo

  ! total thickness of dzten columns 
  Thickness%thicken(:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%thicken(i,j,1) = Thickness%thicken(i,j,1) + Thickness%dzten(i,j,k,1)
           Thickness%thicken(i,j,2) = Thickness%thicken(i,j,2) + Thickness%dzten(i,j,k,2)
        enddo
     enddo 
  enddo


  ! mass per area of U-cell column
  Thickness%mass_u(:,:,:)  = 0.0 
  do n=1,3
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%mass_u(i,j,n) = Thickness%mass_u(i,j,n) + Grid%umask(i,j,k)*Thickness%rho_dzu(i,j,k,n)
           enddo
        enddo
     enddo
  enddo

  ! mass per area of dzten column
  Thickness%mass_en(:,:,:)  = 0.0 
  do n=1,2
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%mass_en(i,j,n) = Thickness%mass_en(i,j,n) + Grid%tmasken(i,j,k,n)*Thickness%rho_dzten(i,j,k,n)
           enddo
        enddo
     enddo
  enddo

  if (use_blobs) then
     Thickness%mass_uT(:,:,:) = 0.0
     do n=1,3
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 Thickness%mass_uT(i,j,n) = Thickness%mass_uT(i,j,n) &
                                            + Grid%umask(i,j,k)*Thickness%rho_dzuT(i,j,k,n)
              enddo
           enddo
        enddo
     enddo
  endif

end subroutine thickness_restart
! </SUBROUTINE> NAME="thickness_restart"

!#######################################################################
! <SUBROUTINE NAME="register_thickness_restart">
!
! <DESCRIPTION>
! Registers fields for the thickness restart file that are not directly
! associated with the partitioning that arises due to blobs.
! </DESCRIPTION>
!
subroutine register_thickness_restart(Time, ignore_blobs, file_name, Thickness)
                                      

  type(ocean_time_type),          intent(in)    :: Time
  logical,                        intent(in)    :: ignore_blobs
  character(len=128),             intent(in)    :: file_name
  type(ocean_thickness_type),     intent(inout) :: Thickness

  integer :: tau, taum1, taup1, id_restart

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  if(tendency==THREE_LEVEL) then 
     id_restart_rho_dzt = register_restart_field(Thk_restart, file_name, 'rho_dzt', Thickness%rho_dzt(:,:,:,tau), &
          Thickness%rho_dzt(:,:,:,taup1), Dom%domain2d)
  elseif(tendency==TWO_LEVEL) then
     if (ignore_blobs .or. .not. use_blobs) then
        id_restart_rho_dzt = register_restart_field(Thk_restart,  file_name, 'rho_dzt',  &
             Thickness%rho_dzt(:,:,:,taup1), Dom%domain2d)
     endif
  endif
  id_restart = register_restart_field(Thk_restart, file_name, 'dzt_dst', &
       Thickness%dzt_dst(:,:,:), Dom%domain2d)
  id_restart = register_restart_field(Thk_restart, file_name, 'dstlo', &
       Thickness%dstlo(:,:,:), Dom%domain2d)
  id_restart = register_restart_field(Thk_restart, file_name, 'dstup', &
       Thickness%dstup(:,:,:), Dom%domain2d)
  id_restart = register_restart_field(Thk_restart, file_name, 'dswt', &
       Thickness%dswt(:,:,:), Dom%domain2d)
  id_restart = register_restart_field(Thk_restart, file_name, 'dst', &
       Thickness%dst(:,:,:), Dom%domain2d)
  id_restart = register_restart_field(Thk_restart, file_name, 'pbot0', &
       Thickness%pbot0(:,:), Dom%domain2d)
  if(vert_coordinate_class==PRESSURE_BASED) then 
     id_restart = register_restart_field(Thk_restart, file_name, 'depth_st', &
          Thickness%depth_st(:,:,:), Dom%domain2d)
     id_restart = register_restart_field(Thk_restart, file_name, 'depth_swt', &
          Thickness%depth_swt(:,:,:), Dom%domain2d)
  endif




end subroutine register_thickness_restart
! </SUBROUTINE> NAME="register_thickness_restart"

!#######################################################################
! <SUBROUTINE NAME="register_thickness_restart_blobs">
!
! <DESCRIPTION>
! Registers restart fields that are directly associated with the 
! partitioning that arises from use of the embedded Lagrangian model.
! </DESCRIPTION>
!
subroutine register_thickness_restart_blobs(Time, ignore_blobs, file_name, Thickness)
                                            

  type(ocean_time_type),          intent(in)    :: Time
  logical,                        intent(in)    :: ignore_blobs
  character(len=128),             intent(in)    :: file_name
  type(ocean_thickness_type),     intent(inout) :: Thickness

  integer :: tau, taum1, taup1, id_restart

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  if (use_blobs .and. .not. ignore_blobs) then
     id_restart_rho_dztL = register_restart_field(Thk_restart, file_name, 'rho_dztL', &
          Thickness%rho_dztL(:,:,:,taup1), Dom%domain2d)
     id_restart_rho_dztT = register_restart_field(Thk_restart, file_name, 'rho_dztT', &
          Thickness%rho_dztT(:,:,:,taup1), Dom%domain2d)
     id_restart = register_restart_field(Thk_restart, file_name, 'dztloL', &
          Thickness%dztloL(:,:,:), Dom%domain2d)
     id_restart = register_restart_field(Thk_restart, file_name, 'dztupL', &
          Thickness%dztupL(:,:,:), Dom%domain2d)
  endif

end subroutine register_thickness_restart_blobs
! </SUBROUTINE> NAME="register_thickness_restart_blobs"


!#######################################################################
! <SUBROUTINE NAME="calculate_restart_thickness">
!
! <DESCRIPTION>
! This subroutine calculates the total thickness and the E contribution
! to thickness.  It is used for the press, pstar and zstar vertical
! coordinates.
! </DESCRIPTION>
!
subroutine calculate_restart_thickness(Thickness,taup1)
  type(ocean_thickness_type), intent(inout) :: Thickness
  integer                   , intent(in)    :: taup1
  
  integer :: i,j,k

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%dztloT(i,j,k)     = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
           Thickness%dztupT(i,j,k)     = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
           Thickness%dztT(i,j,k,taup1) = Thickness%dztloT(i,j,k) + Thickness%dztupT(i,j,k) 
           
           Thickness%dztlo(i,j,k)   = Thickness%dztloT(i,j,k) - Thickness%dztloL(i,j,k)
           Thickness%dztup(i,j,k)   = Thickness%dztupT(i,j,k) - Thickness%dztupL(i,j,k)
           Thickness%dzt(i,j,k)     = Thickness%dztT(i,j,k,taup1)   - Thickness%dztL(i,j,k)
           
           Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
           
           wrk1(i,j,k) = 1.0/Thickness%dzt_dst(i,j,k)
        enddo
     enddo
  enddo
end subroutine calculate_restart_thickness
! </SUBROUTINE> NAME="calculate_restart_thickness"


!#######################################################################
! <SUBROUTINE NAME="dst_land_adjust">
!
! <DESCRIPTION>
! For computing the grid values inside land in global halos.
! These values are not zero, and so they are not available from   
! restart file. They are needed for pressure-based vertical 
! coordinate models in order to get the U-cell values at 
! global land boundaries via the remap operator. 
!
! This step is necessitated by the modifications to dst 
! and other arrays made in ocean_thickness_init_adjust.
! </DESCRIPTION>
!
subroutine dst_land_adjust(Grid, Thickness)

  type(ocean_grid_type),      intent(in)    :: Grid  
  type(ocean_thickness_type), intent(inout) :: Thickness

  integer :: i, j, k
  real    :: density_tmp  

  if(vert_coordinate_class==DEPTH_BASED) return 

  if(vert_coordinate_type==TERRAIN_FOLLOWING) return 

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           density_tmp = rho0_profile(k)*rescale_rho0_mask(i,j)
           if(Grid%tmask(i,j,k) == 0.0) then 
               Thickness%dst(i,j,k)       = -grav*density_tmp*Thickness%dzt(i,j,k)
               Thickness%dzt_dst(i,j,k)   = -grav_r/(density_tmp+epsln)
               Thickness%dzt(i,j,k)       =  Thickness%dst(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%rho_dzt(i,j,k,:) =  density_tmp*Thickness%dzt(i,j,k)
               Thickness%rho_dztr(i,j,k)  =  1.0/(Thickness%rho_dzt(i,j,k,1)+epsln)
           endif
        enddo
     enddo
  enddo

end subroutine dst_land_adjust
! </SUBROUTINE> NAME="dst_land_adjust"


!#######################################################################
! <SUBROUTINE NAME="update_tcell_thickness">
!
! <DESCRIPTION>
! Update time dependent thickness of T cells. When not using the 
! Lagrangian blob scheme.  For routines relevant to the blobs, we use:
! 1/ update_tcell_thic_blob
! 2/ update_E_thickness
! 3/ dzt_dst_update
! which are further down the module.
!
! Notes
!
! 1. For pressure-based coordinates, must use rho(tau) since 
!    rho(taup1) has not yet been computed. This is a limitation of 
!    the "z-like" algorithm approach used in MOM.  It is minor,
!    however, since rho(tau) is very close to rho(taup1).  Also, this 
!    approach comes with the advantage of rendering vertical 
!    physical processes (e.g., diffusion) linear for the general 
!    coordinates in MOM, just as with geopotential models.    
!
! 2. For all coordinates, need to place something reasonable over 
!    land for dz increments to preclude division by zero. 
!    For GEOPOTENTIAL and ZSTAR, we just use the time=zero value for dzt.  
!    For ZSIGMA we set a minimum depth for a fictitious layer over land.
!    For PRESSURE and PSTAR, we set rho=rho0_profile and use the time=0 dzt values.
!    For PSIGMA, we rho=rho0_profile and use the fictitious layer over land.  
!    The treatment of these land-points only contributes to the treatment
!    of the U-cell grid quantities in the halos and on land, as they are 
!    computed via the REMAP_ZT_TO_ZU operator. 
!
! 3. We need to update s-coordinate increments for GEOPOTENTIAL and 
!    PRESSURE. It is only for these two coordinates that we modify the 
!    endpoint values for the s-grid increments.  
!
! </DESCRIPTION>
!
subroutine update_tcell_thickness (Time, Grid, Ext_mode, Dens, Thickness)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid  
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(in)    :: Dens 
  type(ocean_thickness_type),     intent(inout) :: Thickness

  integer :: i, j, k, kb
  integer :: tau, taup1
  integer :: num_bad=0
  logical :: huge_undulations
  logical :: error_flag=.false.
  real    :: density_tmp
  character(len=128) :: errorstring

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1
  huge_undulations=.false.

  ! write tau values of the vertical grid increments before update to taup1 
  call diagnose_3d(Time, Grid, id_dst, Thickness%dst(:,:,:))
  call diagnose_3d(Time, Grid, id_dzt, Thickness%dzt(:,:,:))

  if (id_dzt_pbc > 0) then 
      wrk1_2d(:,:) = 0.0
      do j=jsd,jed
         do i=isd,ied
            kb = Grid%kmt(i,j)
            if (kb > 1) then
                wrk1_2d(i,j) = Thickness%dzt(i,j,kb)
            endif
         enddo
      enddo
      call diagnose_2d(Time, Grid, id_dzt_pbc, wrk1_2d(:,:))
  endif

  call diagnose_3d(Time, Grid, id_dztlo, Thickness%dztlo(:,:,:))
  call diagnose_3d(Time, Grid, id_dztup, Thickness%dztup(:,:,:))
  call diagnose_3d(Time, Grid, id_dzwt, Thickness%dzwt(:,:,1:nk))
  call diagnose_3d(Time, Grid, id_rho_dzt, Thickness%rho_dzt(:,:,:,tau))
  call diagnose_3d(Time, Grid, id_dzt_dst, Thickness%dzt_dst(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_zt, Thickness%depth_zt(:,:,:))
  call diagnose_3d(Time, Grid, id_geodepth_zt, Thickness%geodepth_zt(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_st, Thickness%depth_st(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_swt, Thickness%depth_swt(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_zwt, Thickness%depth_zwt(:,:,:))

  if (id_mass_t > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grid%dat(i,j)*Thickness%rho_dzt(i,j,k,tau)
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grid, id_mass_t, wrk1(:,:,:))
  endif


  ! update coordinate increments for GEOPOTENTIAL and PRESSURE.
  ! only for these two coordinates do we need to modify the 
  ! endpoint values for the s-grid increments.  Other MOM 
  ! coordinates have dst and dswt constant in time.  
  if(vert_coordinate==GEOPOTENTIAL) then 
      do j=jsd,jed
         do i=isd,ied
            Thickness%dswt(i,j,0)  = Thickness%depth_zt(i,j,1)  + Ext_mode%eta_t(i,j,taup1)
            Thickness%dst(i,j,1)   = Thickness%depth_zwt(i,j,1) + Ext_mode%eta_t(i,j,taup1)
            Thickness%dstup(i,j,1) = Thickness%dswt(i,j,0)
            Thickness%dstlo(i,j,1) = Thickness%dst(i,j,1) - Thickness%dstup(i,j,1)
         enddo
      enddo
      if(linear_free_surface) then  
          do j=jsd,jed
             do i=isd,ied
                Thickness%dswt(i,j,0)  = Thickness%depth_zt(i,j,1)  
                Thickness%dst(i,j,1)   = Thickness%depth_zwt(i,j,1) 
                Thickness%dstup(i,j,1) = Thickness%dswt(i,j,0)
                Thickness%dstlo(i,j,1) = Thickness%dst(i,j,1) - Thickness%dstup(i,j,1)
             enddo
          enddo
      endif
  elseif(vert_coordinate==PRESSURE) then 
      do j=jsd,jed
         do i=isd,ied
            Thickness%dswt(i,j,0)     = -Thickness%depth_st(i,j,1)  + Ext_mode%patm_t(i,j,taup1)
            Thickness%dst(i,j,1)      = -Thickness%depth_swt(i,j,1) + Ext_mode%patm_t(i,j,taup1)
            Thickness%dstup(i,j,1)    =  Thickness%dswt(i,j,0)
            Thickness%dstlo(i,j,1)    =  Thickness%dst(i,j,1) - Thickness%dstup(i,j,1)
            kb                        =  Grid%kmt(i,j)
            if(kb > 1) then 
              Thickness%dswt(i,j,kb)  =  Thickness%depth_st(i,j,kb)    - Ext_mode%pbot_t(i,j,taup1)
              Thickness%dst(i,j,kb)   =  Thickness%depth_swt(i,j,kb-1) - Ext_mode%pbot_t(i,j,taup1)
              Thickness%dstlo(i,j,kb) =  Thickness%dswt(i,j,kb)
              Thickness%dstup(i,j,kb) =  Thickness%dst(i,j,kb) - Thickness%dstlo(i,j,kb)
            endif 
         enddo
      enddo
  endif


  ! update specific thicknesses and vertical increments 
  if(vert_coordinate==GEOPOTENTIAL) then 
      do j=jsd,jed
         do i=isd,ied 
            Thickness%dztlo(i,j,1)         = Thickness%dstlo(i,j,1)
            Thickness%dztup(i,j,1)         = Thickness%dstup(i,j,1)
            Thickness%dzt(i,j,1)           = Thickness%dztlo(i,j,1) + Thickness%dztup(i,j,1)           
            Thickness%rho_dzt(i,j,1,taup1) = rho0*Thickness%dzt(i,j,1)
            Thickness%rho_dztr(i,j,1)      = 1.0/(Thickness%rho_dzt(i,j,1,taup1)+epsln)
            wrk1(i,j,:)                    = 1.0  
         enddo
      enddo
  elseif(vert_coordinate==ZSTAR) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzt_dst(i,j,k)       = 1.0 + Ext_mode%eta_t(i,j,taup1)*Grid%htr(i,j)
               Thickness%dztlo(i,j,k)         = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztup(i,j,k)         = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dzt(i,j,k)           = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
               Thickness%rho_dzt(i,j,k,taup1) = rho0*Thickness%dzt(i,j,k)
               Thickness%rho_dztr(i,j,k)      = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               wrk1(i,j,k)                    = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo
  elseif(vert_coordinate==ZSIGMA) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzt_dst(i,j,k)       = max(depth_min_for_sigma,Grid%ht(i,j)+Ext_mode%eta_t(i,j,taup1))
               Thickness%dztlo(i,j,k)         = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztup(i,j,k)         = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dzt(i,j,k)           = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
               wrk1(i,j,k)                    = 1.0/Thickness%dzt_dst(i,j,k)  
               Thickness%rho_dzt(i,j,k,taup1) = rho0*Thickness%dzt(i,j,k)
               Thickness%rho_dztr(i,j,k)      = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PRESSURE) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               ! leave values over land untouched 
               if(Grid%tmask(i,j,k) > 0.0) then 
                  Thickness%dzt_dst(i,j,k)       = -grav_r/(Dens%rho(i,j,k,tau)+epsln)
                  Thickness%dztlo(i,j,k)         = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
                  Thickness%dztup(i,j,k)         = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
                  Thickness%dzt(i,j,k)           = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
                  Thickness%rho_dzt(i,j,k,taup1) = Dens%rho(i,j,k,tau)*Thickness%dzt(i,j,k)
                  Thickness%rho_dztr(i,j,k)      = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               endif     
               wrk1(i,j,k) = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSTAR) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               ! leave values over land untouched 
               if(Grid%tmask(i,j,k) > 0.0) then 
                   Thickness%dzt_dst(i,j,k)       = -grav_r/(Dens%rho(i,j,k,tau)+epsln)                      &
                                                    *(Ext_mode%pbot_t(i,j,taup1)-Ext_mode%patm_t(i,j,taup1)) &
                                                    *Thickness%pbot0r(i,j)
                   Thickness%dztlo(i,j,k)         = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
                   Thickness%dztup(i,j,k)         = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
                   Thickness%dzt(i,j,k)           = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
                   Thickness%rho_dzt(i,j,k,taup1) = Dens%rho(i,j,k,tau)*Thickness%dzt(i,j,k)
                   Thickness%rho_dztr(i,j,k)      = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               endif 
               wrk1(i,j,k) = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSIGMA) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               if(Grid%tmask(i,j,k) > 0.0) then 
                   Thickness%dzt_dst(i,j,k) = -grav_r/(Dens%rho(i,j,k,tau)+epsln) &
                                              *(Ext_mode%pbot_t(i,j,taup1)-Ext_mode%patm_t(i,j,taup1))
                   density_tmp              = Dens%rho(i,j,k,tau)
               else
                  Thickness%dzt_dst(i,j,k) = -depth_min_for_sigma
                  density_tmp              = rho0_profile(k)
               endif 
               Thickness%dztlo(i,j,k)         = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztup(i,j,k)         = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dzt(i,j,k)           = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k) 
               Thickness%rho_dzt(i,j,k,taup1) = density_tmp*Thickness%dzt(i,j,k)
               Thickness%rho_dztr(i,j,k)      = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               wrk1(i,j,k)                    = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo
  endif


  ! vertical distance (m) between t-cell points 
  if(Thickness%method==ENERGETIC) then 
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzwt(i,j,k) = 2.0*Thickness%dswt(i,j,k)/(wrk1(i,j,k)+wrk1(i,j,k+1)) 
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwt(i,j,0) = Thickness%dswt(i,j,0)*Thickness%dzt_dst(i,j,1)
            kb=Grid%kmt(i,j)
            if(kb > 0) then 
                Thickness%dzwt(i,j,kb) = Thickness%dswt(i,j,kb)*Thickness%dzt_dst(i,j,kb)
            endif
         enddo
      enddo
  elseif(Thickness%method==FINITEVOLUME) then 
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwt(i,j,0) = Thickness%dztup(i,j,1)
         enddo
      enddo
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzwt(i,j,k) = Thickness%dztlo(i,j,k) + Thickness%dztup(i,j,k+1) 
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            kb=Grid%kmt(i,j)
            if(kb > 0) then 
                Thickness%dzwt(i,j,kb) = Thickness%dztlo(i,j,kb)
            endif
         enddo
      enddo
  endif

  ! vertical depth (m) of t and w points as 
  ! computed via vertical sum of z-increments 
  if(vert_coordinate /= GEOPOTENTIAL) then 

      ! compute depth from z=eta to t-point and to bottom of t-cell 
      do j=jsd,jed
         do i=isd,ied
            Thickness%depth_zt(i,j,1)  = Thickness%dzwt(i,j,0)
            Thickness%depth_zwt(i,j,1) = Thickness%dzt(i,j,1)
         enddo
      enddo
      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%depth_zt(i,j,k)  = Thickness%depth_zt(i,j,k-1)  + Thickness%dzwt(i,j,k-1)
               Thickness%depth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k-1) + Thickness%dzt(i,j,k)
            enddo
         enddo
      enddo

  endif 

  ! check for negative specific thickness, which occurs when eta_t is too negative, 
  ! in which case the free surface is deeper than the ocean bottom topography. 
  do j=jsc,jec
     do i=isc,iec
        if(Ext_mode%eta_t(i,j,taup1) < -Grid%ht(i,j)) then 
            write(stdoutunit,*) &
                 '==>Error from ocean_thickness: Surface undulations too negative; model unstable'
            write(stdoutunit,*) &
                 'eta_t(',i+Dom%ioff,',',j+Dom%joff,') = ',Ext_mode%eta_t(i,j,taup1)
            write(stdoutunit,*) &
                 'depth(',i+Dom%ioff,',',j+Dom%joff,') = ',Grid%ht(i,j)        
            huge_undulations=.true.
        endif
     enddo
  enddo
  if(huge_undulations) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_thickness_mod: Free surface penetrating rock! Model unstable.')
  endif
 

  ! check for negative thickness at top and bottom 
  error_flag=.false.
  num_bad=0
  do j=jsc,jec
     do i=isc,iec

        kb = Grid%kmt(i,j)
        if(kb>0 .and. num_bad > max_num_bad_print) then

            if(Thickness%dzt(i,j,1) < 0.0) then 
                num_bad = num_bad+1

                write(stdoutunit,'(/a,i3,a,i3,a,i3,a,e10.4)')      &
                     '==>Error: negative thickness at top: dzt(',&
                     i+Dom%ioff,',',j+Dom%joff,',',k,') =', Thickness%dzt(i,j,1)

                if(enforce_positive_dzt) then
                    write(stdoutunit,'(a)') 'Resetting Thickness%dzt to ', thickness_dzt_min
                    Thickness%dzt(i,j,1) = thickness_dzt_min
                else
                    write(errorstring,'(a,i3,a,i3,a,e10.4)') &
                    '(',i+Dom%ioff,',',j+Dom%joff,',1) =', Thickness%dzt(i,j,1)
                    call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, i, j, kb, &
                    'From update_tcell_thickness with dzt<0', unit)
                    error_flag=.true.
                endif

            elseif(Thickness%dzt(i,j,kb) < 0.0) then 
                num_bad = num_bad+1

                write(stdoutunit,'(/a,i3,a,i3,a,i3,a,e10.4)')    &
                '==>Error: negative thickness at bottom: dzt(',&
                i+Dom%ioff,',',j+Dom%joff,',',kb,') =', Thickness%dzt(i,j,kb)

                if(enforce_positive_dzt) then
                    write(stdoutunit,'(a)') 'Resetting Thickness%dzt to ', thickness_dzt_min
                    Thickness%dzt(i,j,kb) = thickness_dzt_min
                else
                    write(errorstring,'(a,i3,a,i3,a,i3,a,e10.4)') &
                    '(',i+Dom%ioff,',',j+Dom%joff,',',kb,') =', Thickness%dzt(i,j,kb)
                    call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, i, j, kb, &
                    'From update_tcell_thickness with dzt<0', unit)
                    error_flag=.true.
                endif

            endif

        endif

     enddo
  enddo

  if(error_flag) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_thickness_mod: dzt'//trim(errorstring)//'<0. grep for "dzt" to find location')
  endif


  if(debug_this_module_detail) then 
     call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, isc, jsc, nk, &
                            'From update_tcell_thickness', stdoutunit)
     call dzt_min_max(Time, Thickness, &
     'From update_tcell_thickness with dzt<0, dzt_min_max information')
  endif


end subroutine update_tcell_thickness
! </SUBROUTINE> NAME="update_tcell_thickness"


!#######################################################################
! <SUBROUTINE NAME="update_ucell_thickness">
!
! <DESCRIPTION>
! Update time dependent thickness of U cells.
!
! Notes
!
! The computation of the depth arrays is a bit ad hoc.  Here are 
! the methods and their rationale.
!
! 1. For terrain following coordinates, the remap operator will entrain 
!    the very shallow layer T-points over land into the four-point averaging.
!    This will cause the resulting U-cell depths next to land to be far shallower
!    than what they should be.  To avoid this problem, we set the U-cell depths the
!    same as the T-cell depths. This specification causes no problems for budgets 
!    or energetic balance since Grid%umask array is already defined according to  
!    the usual B-grid specification.   
!
! 2. For non-terrain and non-geopotential coordinates, compute U-cell thicknesses 
!    as min of surrounding T-cell thicknesses.  This method follows that 
!    used with earlier MOM versions for the bottom partial step topography. 
!    Experiments with ZSTAR reveal that we must use the min function for dzu 
!    computation, rather than the alternative of a remap operator.  If using 
!    the remap operator and nontrivial topography, then the velocity field can 
!    develop nontrivial noise.  We have found that can compute dzwu using the 
!    remap operator without introducing noise. But we choose to use the 
!    min function to maintain compatibility with traditional approach in 
!    geopotential vertical coordinate versions of MOM. 
!
! If use_blobs=.true. in the ocean_model_nml, we adopt the following:
!  *dzuT = min of surrounting T-cell thickness,
!  *dzuL = min of surrounting T-cell thickness,
!  *dzu  = *dzuT - *dztL
! </DESCRIPTION>
!
subroutine update_ucell_thickness (Time, Grid, Ext_mode, Thickness)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid  
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_thickness_type),     intent(inout) :: Thickness

  integer :: i, j, k, n
  integer :: tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  ! write tau values of grid increments before update to taup1 
  
  if (use_blobs) then
     call diagnose_3d_u(Time, Grid, id_dzu, Thickness%dzuT(:,:,:))
     call diagnose_3d_u(Time, Grid, id_dzuL, Thickness%dzuL(:,:,:))
     call diagnose_3d_u(Time, Grid, id_dzuE, Thickness%dzu(:,:,:))
     call diagnose_3d_u(Time, Grid, id_dzwu, Thickness%dzwuT(:,:,1:nk))
     call diagnose_3d_u(Time, Grid, id_dzwuL, Thickness%dzwuL(:,:,1:nk))
     call diagnose_3d_u(Time, Grid, id_dzwuE, Thickness%dzwu(:,:,1:nk))
     call diagnose_3d_u(Time, Grid, id_rho_dzu, Thickness%rho_dzuT(:,:,:,tau))
     call diagnose_3d_u(Time, Grid, id_rho_dzuL, Thickness%rho_dzuL(:,:,:,tau))
     call diagnose_3d_u(Time, Grid, id_rho_dzuE, Thickness%rho_dzu(:,:,:,tau))
     call diagnose_2d_u(Time, Grid, id_mass_u, Thickness%mass_uT(:,:,tau))
     call diagnose_2d_u(Time, Grid, id_mass_uE, Thickness%mass_u(:,:,tau))
  else
     call diagnose_3d_u(Time, Grid, id_dzu, Thickness%dzu(:,:,:))
     call diagnose_3d_u(Time, Grid, id_dzten(1), Thickness%dzten(:,:,:,1))
     call diagnose_3d_u(Time, Grid, id_dzten(2), Thickness%dzten(:,:,:,2))
     call diagnose_3d_u(Time, Grid, id_dzwu, Thickness%dzwu(:,:,1:nk))
     call diagnose_3d_u(Time, Grid, id_rho_dzu, Thickness%rho_dzu(:,:,:,tau))
     call diagnose_2d_u(Time, Grid, id_mass_u, Thickness%mass_u(:,:,tau))
     call diagnose_2d_en(Time, Grid, id_mass_en(1), id_mass_en(2), Thickness%mass_en(:,:,:))
  endif

  call diagnose_3d(Time, Grid, id_dzten(1), Thickness%dzten(:,:,:,1))
  call diagnose_3d(Time, Grid, id_dzten(2), Thickness%dzten(:,:,:,2))
  call diagnose_3d(Time, Grid, id_rho_dzten(1), Thickness%rho_dzten(:,:,:,1))
  call diagnose_3d(Time, Grid, id_rho_dzten(2), Thickness%rho_dzten(:,:,:,2))

  call diagnose_3d_u(Time, Grid, id_depth_zu, Thickness%depth_zu(:,:,:))
  call diagnose_3d_u(Time, Grid, id_depth_zwu, Thickness%depth_zwu(:,:,:))

  call diagnose_2d_u(Time, Grid, id_thicku, Thickness%thicku(:,:,tau))
  call diagnose_2d_en(Time, Grid, id_thicken(1), id_thicken(2), Thickness%thicken(:,:,:))

  if (use_blobs) then
     do j=jsd,jec
        do i=isd,iec
           Thickness%dzwuL(i,j,0) = min(Thickness%dzwtL(i,  j,  0),&
                                        Thickness%dzwtL(i+1,j,  0),&
                                        Thickness%dzwtL(i,  j+1,0),&
                                        Thickness%dzwtL(i+1,j+1,0))
        enddo
     enddo
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              Thickness%dzuL(i,j,k)  = min(Thickness%dztL( i,  j,  k), &
                                           Thickness%dztL( i+1,j,  k), &
                                           Thickness%dztL( i,  j+1,k), &
                                           Thickness%dztL( i+1,j+1,k))

              Thickness%dzwuL(i,j,k) = min(Thickness%dzwtL(i,  j,  k), &
                                           Thickness%dzwtL(i+1,j,  k), &
                                           Thickness%dzwtL(i,  j+1,k), &
                                           Thickness%dzwtL(i+1,j+1,k))

              Thickness%rho_dzuL(i,j,k,taup1) = min(Thickness%rho_dztL(i,  j,  k,taup1),&
                                                    Thickness%rho_dztL(i+1,j,  k,taup1),&
                                                    Thickness%rho_dztL(i,  j+1,k,taup1),&
                                                    Thickness%rho_dztL(i+1,j+1,k,taup1))
           enddo
        enddo
     enddo
     call mpp_update_domains(Thickness%dzuL(:,:,:),           Dom%domain2d)
     call mpp_update_domains(Thickness%dzwuL(:,:,:),          Dom%domain2d)
     call mpp_update_domains(Thickness%rho_dzuL(:,:,:,taup1), Dom%domain2d)
  endif !use_blob


  if(vert_coordinate==GEOPOTENTIAL) then

     if (use_blobs) then
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwuT(i,j,0)          = Thickness%depth_zt(i,j,1)  + Ext_mode%eta_u(i,j,taup1)
              Thickness%dzwu(i,j,0)           = Thickness%dzwuT(i,j,0) - Thickness%dzwuL(i,j,0)
              Thickness%dzuT(i,j,1)           = Thickness%depth_zwt(i,j,1) + Ext_mode%eta_u(i,j,taup1)
              Thickness%rho_dzuT(i,j,1,taup1) = rho0*Thickness%dzuT(i,j,1)
           enddo
        enddo
        if(linear_free_surface) then  
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwuT(i,j,0)          = Thickness%depth_zt(i,j,1)
                 Thickness%dzwu(i,j,0)           = Thickness%dzwuT(i,j,0) - Thickness%dzwuL(i,j,0)
                 Thickness%dzuT(i,j,1)           = Thickness%depth_zwt(i,j,1) 
                 Thickness%rho_dzuT(i,j,1,taup1) = rho0*Thickness%dzu(i,j,1)
              enddo
           enddo
           do k=1,nk
              do j=jsd,jed
                 do i=isd,ied
                    Thickness%dzu(i,j,k)           = Thickness%dzuT(i,j,k)           - Thickness%dzuL(i,j,k)
                    Thickness%rho_dzu(i,j,k,taup1) = Thickness%rho_dzuT(i,j,k,taup1) - Thickness%rho_dzuL(i,j,k,taup1)
                    Thickness%rho_dzur(i,j,k)      = 1.0/(Thickness%rho_dzu(i,j,k,taup1)+epsln)
                 enddo
              enddo
           enddo
        endif

     else !do not use_blobs

        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwu(i,j,0)          = Thickness%depth_zt(i,j,1)  + Ext_mode%eta_u(i,j,taup1)
              Thickness%dzu(i,j,1)           = Thickness%depth_zwt(i,j,1) + Ext_mode%eta_u(i,j,taup1)
              Thickness%rho_dzu(i,j,1,taup1) = rho0*Thickness%dzu(i,j,1)
              Thickness%rho_dzur(i,j,1)      = 1.0/(Thickness%rho_dzu(i,j,1,taup1)+epsln)
           enddo
        enddo
        if(linear_free_surface) then  
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzwu(i,j,0)          = Thickness%depth_zt(i,j,1)
                 Thickness%dzu(i,j,1)           = Thickness%depth_zwt(i,j,1) 
                 Thickness%rho_dzu(i,j,1,taup1) = rho0*Thickness%dzu(i,j,1)
                 Thickness%rho_dzur(i,j,1)      = 1.0/(Thickness%rho_dzu(i,j,1,taup1)+epsln)
              enddo
           enddo
        endif
     endif !use_blobs

  else ! not GEOPOTENTIAL

     ! Note: zsigma and psigma not supported for blobs
     if(vert_coordinate==ZSIGMA .or. vert_coordinate==PSIGMA) then 
        do k=1,nk 
           do j=jsd,jed
              do i=isd,ied
                 Thickness%dzu(i,j,k)           = Thickness%dzt(i,j,k)
                 Thickness%dzten(i,j,k,1)       = Thickness%dzt(i,j,k)
                 Thickness%dzten(i,j,k,2)       = Thickness%dzt(i,j,k)
                 Thickness%dzwu(i,j,k)          = Thickness%dzwt(i,j,k)
                 Thickness%rho_dzu(i,j,k,taup1) = Thickness%rho_dzt(i,j,k,taup1)
              enddo
           enddo
        enddo
        k=0
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzwu(i,j,k) = Thickness%dzwt(i,j,k)
           enddo
        enddo
        
      else !not zsigma or psigma

         if (use_blobs) then

            ! U-cell thicknesses as minimum of surrounding T-cell thicknesses
            do j=jsd,jec
               do i=isd,iec
                  Thickness%dzwuT(i,j,0) = min(Thickness%dzwtT(i,  j,  0),&
                                               Thickness%dzwtT(i+1,j,  0),&
                                               Thickness%dzwtT(i,  j+1,0),&
                                               Thickness%dzwtT(i+1,j+1,0))
               enddo
            enddo

            do k=1,nk
               do j=jsd,jec
                  do i=isd,iec
                     Thickness%dzuT(i,j,k)  = min(Thickness%dztT(i,  j,  k,taup1), &
                                                  Thickness%dztT(i+1,j,  k,taup1), &
                                                  Thickness%dztT(i,  j+1,k,taup1), &
                                                  Thickness%dztT(i+1,j+1,k,taup1))
                     
                     Thickness%dzwuT(i,j,k) = min(Thickness%dzwtT(i,  j,  k), &
                                                  Thickness%dzwtT(i+1,j,  k), &
                                                  Thickness%dzwtT(i,  j+1,k), &
                                                  Thickness%dzwtT(i+1,j+1,k))
                     
                     Thickness%rho_dzuT(i,j,k,taup1) = min(Thickness%rho_dztT(i,  j,  k,taup1), &
                                                           Thickness%rho_dztT(i+1,j,  k,taup1), &
                                                           Thickness%rho_dztT(i,  j+1,k,taup1), &
                                                           Thickness%rho_dztT(i+1,j+1,k,taup1))
                  enddo
               enddo
            enddo

            call mpp_update_domains(Thickness%dzuT(:,:,:),           Dom%domain2d)
            call mpp_update_domains(Thickness%dzwuT(:,:,:),          Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzuT(:,:,:,taup1), Dom%domain2d)
            
            Thickness%dzu(:,:,:)  = Thickness%dzuT(:,:,:)  - Thickness%dzuL(:,:,:)
            Thickness%dzwu(:,:,:) = Thickness%dzwuT(:,:,:) - Thickness%dzwuL(:,:,:)
            Thickness%rho_dzu(:,:,:,taup1) = Thickness%rho_dzuT(:,:,:,taup1) - Thickness%rho_dzuL(:,:,:,taup1)

         else !not use_blobs
            
            ! U-cell thicknesses as minimum of surrounding T-cell thicknesses
            do k=1,nk
               do j=jsd,jec
                  do i=isd,iec

                     Thickness%dzu(i,j,k)  = min(Thickness%dzt(i,j,k),   &
                                                 Thickness%dzt(i+1,j,k), &
                                                 Thickness%dzt(i,j+1,k), &
                                                 Thickness%dzt(i+1,j+1,k))
                     Thickness%dzwu(i,j,k) = min(Thickness%dzwt(i,j,k),   &
                                                 Thickness%dzwt(i+1,j,k), &
                                                 Thickness%dzwt(i,j+1,k), &
                                                 Thickness%dzwt(i+1,j+1,k))
                     Thickness%rho_dzu(i,j,k,taup1) = min(Thickness%rho_dzt(i,j,k,taup1),    &
                                                          Thickness%rho_dzt(i+1,j,k,taup1),  &
                                                          Thickness%rho_dzt(i,j+1,k,taup1),  &
                                                          Thickness%rho_dzt(i+1,j+1,k,taup1))

                     Thickness%dzten(i,j,k,1)     = min(Thickness%dzt(i,j,k),&
                                                        Thickness%dzt(i+1,j,k))
                     Thickness%dzten(i,j,k,2)     = min(Thickness%dzt(i,j,k),&
                                                        Thickness%dzt(i,j+1,k))
                     Thickness%rho_dzten(i,j,k,1) = min(Thickness%rho_dzt(i,j,k,taup1),&
                                                        Thickness%rho_dzt(i+1,j,k,taup1))
                     Thickness%rho_dzten(i,j,k,2) = min(Thickness%rho_dzt(i,j,k,taup1),&
                                                        Thickness%rho_dzt(i,j+1,k,taup1))

                  enddo
               enddo
            enddo
            if(update_dzwu_k0) then 
               k=0
               do j=jsd,jec
                  do i=isd,iec
                     Thickness%dzwu(i,j,k) = min(Thickness%dzwt(i,j,k),   &
                                                 Thickness%dzwt(i+1,j,k), &
                                                 Thickness%dzwt(i,j+1,k), &
                                                 Thickness%dzwt(i+1,j+1,k))
                  enddo
               enddo
            endif
            
            call mpp_update_domains(Thickness%dzu(:,:,:),           Dom%domain2d)
            call mpp_update_domains(Thickness%dzwu(:,:,:),          Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzu(:,:,:,taup1), Dom%domain2d)

            call mpp_update_domains(Thickness%dzten(:,:,:,1),       Dom%domain2d)
            call mpp_update_domains(Thickness%dzten(:,:,:,2),       Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzten(:,:,:,1),   Dom%domain2d)
            call mpp_update_domains(Thickness%rho_dzten(:,:,:,2),   Dom%domain2d)

            
         endif !use_blobs


      endif !zsigma or psigma


      do k=1,nk 
         do j=jsd,jed
            do i=isd,ied
               Thickness%rho_dzur(i,j,k) = Grid%umask(i,j,k)/(Thickness%rho_dzu(i,j,k,taup1)+epsln)
            enddo
         enddo
      enddo

      if (use_blobs) then
         do j=jsd,jed
            do i=isd,ied
               Thickness%depth_zu(i,j,1)  = Thickness%dzwuT(i,j,0)
               Thickness%depth_zwu(i,j,1) = Thickness%dzuT(i,j,1)
            enddo
         enddo
         do k=2,nk
            do j=jsd,jed
               do i=isd,ied
                  Thickness%depth_zu(i,j,k)  = Thickness%depth_zu(i,j,k-1)  + Thickness%dzwuT(i,j,k-1)
                  Thickness%depth_zwu(i,j,k) = Thickness%depth_zwu(i,j,k-1) + Thickness%dzuT(i,j,k) 
               enddo
            enddo
         enddo
         
      else !do not use blobs

         do j=jsd,jed
            do i=isd,ied
               Thickness%depth_zu(i,j,1)  = Thickness%dzwu(i,j,0)
               Thickness%depth_zwu(i,j,1) = Thickness%dzu(i,j,1)
            enddo
         enddo
         do k=2,nk
            do j=jsd,jed
               do i=isd,ied
                  Thickness%depth_zu(i,j,k)  = Thickness%depth_zu(i,j,k-1)  + Thickness%dzwu(i,j,k-1)
                  Thickness%depth_zwu(i,j,k) = Thickness%depth_zwu(i,j,k-1) + Thickness%dzu(i,j,k) 
               enddo
            enddo
         enddo
      endif !use_blobs

   endif !geopotential


  ! thickness of U-cell column 
  do j=jsd,jed
     do i=isd,ied
        k=Grid%kmu(i,j)
        if(k > 1) then 
            Thickness%thicku(i,j,taup1) = Thickness%depth_zwu(i,j,k) 
        else
            Thickness%thicku(i,j,taup1) = 0.0
        endif
     enddo
  enddo

  ! total thickness of dzten columns 
  Thickness%thicken(:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%thicken(i,j,1) = Thickness%thicken(i,j,1) + Thickness%dzten(i,j,k,1)
           Thickness%thicken(i,j,2) = Thickness%thicken(i,j,2) + Thickness%dzten(i,j,k,2)
        enddo
     enddo 
  enddo


  ! mass per area of U-cell column
  Thickness%mass_u(:,:,taup1)  = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%mass_u(i,j,taup1)  = Thickness%mass_u(i,j,taup1) + Grid%umask(i,j,k)*Thickness%rho_dzu(i,j,k,taup1)
        enddo
     enddo
  enddo

 ! total mass per area of rho_dzten columns 
 Thickness%mass_en(:,:,:) = 0.0
 do n=1,2
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             Thickness%mass_en(i,j,n) = Thickness%mass_en(i,j,n) + Grid%tmasken(i,j,k,n)*Thickness%rho_dzten(i,j,k,n)
          enddo
       enddo 
    enddo
  enddo

  if (use_blobs) then
     Thickness%mass_uT(:,:,taup1) = 0.0
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%mass_uT(i,j,taup1) = Thickness%mass_uT(i,j,taup1) &
                                             + Grid%umask(i,j,k)*Thickness%rho_dzuT(i,j,k,taup1)
           enddo
        enddo
     enddo
  endif

  if(debug_this_module) then 
    write(stdoutunit,*) ' '
    write(stdoutunit,*) ' From ocean_thickness_mod: thickness checksum'
    call write_timestamp(Time%model_time)
    if (use_blobs) then
       call thickness_chksum_blobs(Time, Grid, Thickness)
    else
       call thickness_chksum(Time, Grid, Thickness)
    endif
  endif 


end subroutine update_ucell_thickness
! </SUBROUTINE> NAME="update_ucell_thickness"


!#######################################################################
! <SUBROUTINE NAME="rho_dzt_tendency">
!
! <DESCRIPTION>
!  Compute tendency for rho_dzt. This tendency is a function of the 
!  vertical coordinate. 
! </DESCRIPTION>
!
subroutine rho_dzt_tendency(Time, Grid, Ext_mode, Thickness)
  
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid  
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_thickness_type),     intent(inout) :: Thickness 

  integer :: i, j, k
  integer :: tau, taup1

  tau   = Time%tau
  taup1 = Time%taup1

  if(vert_coordinate==GEOPOTENTIAL) then 
      k=1
      do j=jsd,jed
         do i=isd,ied       
            Thickness%rho_dzt_tendency(i,j,k) = Grid%tmask(i,j,k)*rho0*Ext_mode%deta_dt(i,j)
         enddo
      enddo
  elseif(vert_coordinate==ZSTAR) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied       
               Thickness%rho_dzt_tendency(i,j,k) = Grid%tmask(i,j,k)*rho0*Ext_mode%deta_dt(i,j) &
                                                  *Grid%htr(i,j)*Thickness%dst(i,j,k)
            enddo
         enddo
      enddo
  elseif(vert_coordinate==ZSIGMA) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied       
               Thickness%rho_dzt_tendency(i,j,k) = Grid%tmask(i,j,k)*rho0*Ext_mode%deta_dt(i,j) &
                                                  *Thickness%dst(i,j,k)
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PRESSURE) then 
      k=1
      do j=jsd,jed
         do i=isd,ied       
            Thickness%rho_dzt_tendency(i,j,k) = -grav_r*Grid%tmask(i,j,k)*Ext_mode%dpatm_dt(i,j)
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            k=Grid%kmt(i,j)       
            if(k>0) then  
                Thickness%rho_dzt_tendency(i,j,k) = grav_r*Grid%tmask(i,j,k)*Ext_mode%dpbot_dt(i,j)
            endif
         enddo
      enddo
  elseif(vert_coordinate==PSTAR) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied       
               Thickness%rho_dzt_tendency(i,j,k) = -grav_r*Grid%tmask(i,j,k)*Thickness%dst(i,j,k) &
                    *(Ext_mode%dpbot_dt(i,j)-Ext_mode%dpatm_dt(i,j))*Thickness%pbot0r(i,j)
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSIGMA) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied       
                   Thickness%rho_dzt_tendency(i,j,k) = -grav_r*Grid%tmask(i,j,k)*Thickness%dst(i,j,k)   &
                                                       *(Ext_mode%dpbot_dt(i,j)-Ext_mode%dpatm_dt(i,j))
            enddo
         enddo
      enddo
  endif 

  call diagnose_3d(Time, Grid, id_rho_dzt_tendency, Thickness%rho_dzt_tendency(:,:,:))

  if (id_mass_tendency > 0) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied       
               wrk1(i,j,k) = Grid%dat(i,j)*Thickness%rho_dzt_tendency(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grid, id_mass_tendency, wrk1(:,:,:))
  endif



end subroutine rho_dzt_tendency
! </SUBROUTINE> NAME="rho_dzt_tendency"


!#######################################################################
! <SUBROUTINE NAME="thickness_chksum">
!
! <DESCRIPTION>
!
! Compute checksum for thickness components .
! Only print checksums for fields that should agree across restarts.  
!
! </DESCRIPTION>
subroutine thickness_chksum(Time, Grid, Thickness)
 
  type(ocean_time_type),      intent(in) :: Time
  type(ocean_grid_type),      intent(in) :: Grid 
  type(ocean_thickness_type), intent(in) :: Thickness

  integer :: tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  if(tendency == THREE_LEVEL) then 
     call write_chksum_3d('Thickness%rho_dzt(tau)', Thickness%rho_dzt(COMP,:,tau)*Grid%tmask(COMP,:))
     call write_chksum_3d('Thickness%rho_dzt(taup1)', Thickness%rho_dzt(COMP,:,taup1)*Grid%tmask(COMP,:))
     call write_chksum_3d('Thickness%rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grid%umask(COMP,:))
     call write_chksum_3d('Thickness%rho_dzu(taup1)', Thickness%rho_dzu(COMP,:,taup1)*Grid%umask(COMP,:))
     call write_chksum_2d('Thickness%mass_u(tau)', Thickness%mass_u(COMP,tau)*Grid%umask(COMP,1))
     call write_chksum_2d('Thickness%mass_u(taup1)', Thickness%mass_u(COMP,taup1)*Grid%umask(COMP,1))
  else 
     call write_chksum_3d('Thickness%rho_dzt(taup1)', Thickness%rho_dzt(COMP,:,taup1)*Grid%tmask(COMP,:))
     call write_chksum_3d('Thickness%rho_dzu(taup1)', Thickness%rho_dzu(COMP,:,taup1)*Grid%umask(COMP,:))
     call write_chksum_2d('Thickness%mass_u(taup1)', Thickness%mass_u(COMP,taup1)*Grid%umask(COMP,1))
  endif
  call write_chksum_3d('Thickness%rho_dzten(1)', Thickness%rho_dzten(COMP,:,1)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzten(2)', Thickness%rho_dzten(COMP,:,2)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dztr', Thickness%rho_dztr(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzur', Thickness%rho_dzur(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzt_tendency', Thickness%rho_dzt_tendency(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dzt', Thickness%dzt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dzten(1)', Thickness%dzten(COMP,:,1)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dzten(2)', Thickness%dzten(COMP,:,2)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztlo', Thickness%dztlo(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztup', Thickness%dztup(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dzt_dst', Thickness%dzt_dst(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_2d('Thickness%dzwt(k=0)', Thickness%dzwt(COMP,0)*Grid%tmask(COMP,1))
  call write_chksum_3d('Thickness%dzwt(k=1:nk)', Thickness%dzwt(COMP,1:nk)*Grid%tmask(COMP,1:nk))
  call write_chksum_3d('Thickness%dzu', Thickness%dzu(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_2d('Thickness%dzwu(k=0)', Thickness%dzwu(COMP,0)*Grid%umask(COMP,1))
  call write_chksum_3d('Thickness%dzwu(k=1:nk)', Thickness%dzwu(COMP,1:nk)*Grid%umask(COMP,1:nk))
  call write_chksum_3d('Thickness%depth_zt', Thickness%depth_zt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%geodepth_zt', Thickness%geodepth_zt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_zu', Thickness%depth_zu(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%depth_zwt', Thickness%depth_zwt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_zwu', Thickness%depth_zwu(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%depth_st', Thickness%depth_st(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_swt', Thickness%depth_swt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dst', Thickness%dst(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dstlo', Thickness%dstlo(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dstup', Thickness%dstup(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_2d('Thickness%dswt(k=0)', Thickness%dswt(COMP,0)*Grid%tmask(COMP,1))
  call write_chksum_3d('Thickness%dswt(k=1:nk)', Thickness%dswt(COMP,1:nk)*Grid%tmask(COMP,1:nk))
  call write_chksum_2d('Thickness%pbot0', Thickness%pbot0(COMP)*Grid%tmask(COMP,1))
  call write_chksum_2d('Thickness%mass_en(1)', Thickness%mass_en(COMP,1)*Grid%tmasken(COMP,1,1))
  call write_chksum_2d('Thickness%mass_en(2)', Thickness%mass_en(COMP,2)*Grid%tmasken(COMP,1,2))
  
  write(stdoutunit,*) ' '
  
  return
  
end subroutine thickness_chksum
! </SUBROUTINE> NAME="thickness_chksum"


!#######################################################################
! <SUBROUTINE NAME="thickness_details">
!
! <DESCRIPTION>
!
! For debugging, we print here some details of the grid at a particular
! (i,j) point. 
!
! </DESCRIPTION>
subroutine thickness_details(Grid, Time, Thickness, Ext_mode, Dens, &
                             ipoint, jpoint, kb, filecaller, outunit)

  type(ocean_grid_type),          intent(in) :: Grid  
  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  type(ocean_density_type),       intent(in) :: Dens
  integer,                        intent(in) :: ipoint
  integer,                        intent(in) :: jpoint
  integer,                        intent(in) :: kb
  character(len=*),               intent(in) :: filecaller
  integer,                        intent(in) :: outunit

  integer :: i,j,k
  integer :: taup1, tau

  taup1 = Time%taup1
  tau   = Time%tau

  write(outunit,*) ' ' 
  write(outunit,'(a)') trim(filecaller)
  call write_timestamp(Time%model_time)

  i=ipoint 
  j=jpoint

  write(outunit,*) ' ' 
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'ht(',i+Dom%ioff,',',j+Dom%joff,')(metres)    = ',&
       Grid%ht(i,j) 
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'eta_t(',i+Dom%ioff,',',j+Dom%joff,')(metres) = ',&
       Ext_mode%eta_t(i,j,taup1)
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'pbot_t(',i+Dom%ioff,',',j+Dom%joff,')(dbars) = ',&
       Ext_mode%pbot_t(i,j,taup1)
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'patm_t(',i+Dom%ioff,',',j+Dom%joff,')(dbars) = ',&
       Ext_mode%patm_t(i,j,taup1)
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'pbot0(',i+Dom%ioff,',',j+Dom%joff,')(dbars)  = ',&
       Thickness%pbot0(i,j)*c2dbars
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'pbot0r(',i+Dom%ioff,',',j+Dom%joff,')(1/dbar)= ',&
       Thickness%pbot0r(i,j)/c2dbars
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'xt(',i+Dom%ioff,',',j+Dom%joff,')(longitude) = ',&
       Grid%xt(i,j)       
  write(outunit,'(a,i4,a,i4,a,f24.8)')                   &
       'yt(',i+Dom%ioff,',',j+Dom%joff,')(latitude) = ', &
       Grid%yt(i,j)       

  write(outunit,*) ' ' 
  do k=1,kb
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                        &
          'depth_st(',i+Dom%ioff,',',j+Dom%joff,',',k,')(s-units)  = ',&
          Thickness%depth_st(i,j,k)*convert_factor
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                        &
          'depth_swt(',i+Dom%ioff,',',j+Dom%joff,',',k,')(s-units) = ',&
          Thickness%depth_swt(i,j,k)*convert_factor
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                        &
          'rho_t(',i+Dom%ioff,',',j+Dom%joff,',',k,')(kg/m^3)      = ',&
          Dens%rho(i,j,k,tau)
  enddo

  write(outunit,*) ' ' 
  do k=1,kb
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                    &
          'dst  (',i+Dom%ioff,',',j+Dom%joff,',',k,')(s-units)= ',&
          Thickness%dst(i,j,k)*convert_factor
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                    &
          'dswt (',i+Dom%ioff,',',j+Dom%joff,',',k,')(s-units)= ',&
          Thickness%dswt(i,j,k)*convert_factor
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                    &
          'dstlo(',i+Dom%ioff,',',j+Dom%joff,',',k,')(s-units)= ',&
          Thickness%dstlo(i,j,k)*convert_factor
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                    &
          'dstup(',i+Dom%ioff,',',j+Dom%joff,',',k,')(s-units)= ',&
          Thickness%dstup(i,j,k)*convert_factor
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                    &
          'dstlo + dstup(',i+Dom%ioff,',',j+Dom%joff,',',k,') = ',&
          (Thickness%dstlo(i,j,k)+Thickness%dstup(i,j,k))*convert_factor
  enddo

  write(outunit,*) ' ' 
  do k=1,kb
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                      &
          'depth_zt(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres) = ',&
          Thickness%depth_zt(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                      &
          'depth_zu(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres) = ',&
          Thickness%depth_zu(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                      &
          'depth_zwt(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)= ',&
          Thickness%depth_zwt(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)')                      &
          'depth_zwu(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)= ',&
          Thickness%depth_zwu(i,j,k)
  enddo

  write(outunit,*) ' ' 
  do k=1,kb
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dzt(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)     = ',Thickness%dzt(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dzten(',i+Dom%ioff,',',j+Dom%joff,',',k,',1)(metres) = ',Thickness%dzten(i,j,k,1)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dzten(',i+Dom%ioff,',',j+Dom%joff,',',k,',2)(metres) = ',Thickness%dzten(i,j,k,2)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dzu(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)     = ',Thickness%dzu(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dzwt(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)    = ',Thickness%dzwt(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dzwu(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)    = ',Thickness%dzwu(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dztlo(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)   = ',Thickness%dztlo(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dztup(',i+Dom%ioff,',',j+Dom%joff,',',k,')(metres)   = ',Thickness%dztup(i,j,k)
     write(outunit,'(a,i4,a,i4,a,i4,a,f24.8)') &
          'dztlo+dztup(',i+Dom%ioff,',',j+Dom%joff,',',k,')     = ',&
           Thickness%dztlo(i,j,k)+Thickness%dztup(i,j,k)
  enddo

  
end subroutine thickness_details
! </SUBROUTINE> NAME="thickness_details"


!#######################################################################
! <SUBROUTINE NAME="ocean_thickness_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_thickness_restart(Time, Thickness, introduce_blobs, time_stamp)
  type(ocean_time_type),      intent(in)           :: Time
  type(ocean_thickness_type), intent(inout)        :: Thickness
  logical,                    intent(in)           :: introduce_blobs
  character(len=*),           intent(in), optional :: time_stamp

  character(len=128) :: filename
  integer :: tau, taup1

  tau    = Time%tau
  taup1  = Time%taup1

  if (introduce_blobs) then
     ! If we have introduced blobs this run, we need to change the variables that are saved
     ! in the restarts.  So, we create a new restart file object with all the fields we
     ! need saved at the end of this run.
     filename = 'ocean_thickness.res.nc'
     call register_thickness_restart_blobs(Time, .false., filename, Thickness)
  endif

  if(tendency==THREE_LEVEL) then 
     call reset_field_pointer(Thk_restart, id_restart_rho_dzt, Thickness%rho_dzt(:,:,:,tau), &
          Thickness%rho_dzt(:,:,:,taup1) )
  elseif(tendency==TWO_LEVEL) then
     if (use_blobs) then
        call reset_field_pointer(Thk_restart, id_restart_rho_dztL, Thickness%rho_dztL(:,:,:,taup1) )
        call reset_field_pointer(Thk_restart, id_restart_rho_dztT, Thickness%rho_dztT(:,:,:,taup1) )
     else
        call reset_field_pointer(Thk_restart, id_restart_rho_dzt,  Thickness%rho_dzt(:,:,:,taup1) )
     endif
  endif

  call save_restart(Thk_restart, time_stamp)

end subroutine ocean_thickness_restart
! </SUBROUTINE> NAME="ocean_thickness_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_thickness_end">
!
! <DESCRIPTION>
! Write basic elements of thickness derived type to restart 
! </DESCRIPTION>
!
subroutine ocean_thickness_end (Time, Grid, introduce_blobs, Thickness)

  type(ocean_time_type),      intent(in)           :: Time
  type(ocean_grid_type),      intent(in)           :: Grid  
  logical,                    intent(in)           :: introduce_blobs
  type(ocean_thickness_type), intent(inout)        :: Thickness

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_thickness_mod (ocean_thickness_end): module must be initialized')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_thickness_mod (ocean_thickness_end): NO restart written.'
    call mpp_error(WARNING, &
    '==>Warning from ocean_thickness_mod (ocean_thickness_end): NO restart written.')
    return
  endif 

  call ocean_thickness_restart(Time, Thickness, introduce_blobs)

  write(stdoutunit,*) ' '
  write(stdoutunit,*) 'From ocean_thickness_mod: ending thickness checksums'
  call write_timestamp(Time%model_time)
  if (use_blobs) then
     call thickness_chksum_blobs(Time, Grid, Thickness)
  else
     call thickness_chksum(Time, Grid, Thickness)
  endif

  nullify(Dom)

end subroutine ocean_thickness_end
! </SUBROUTINE> NAME="ocean_thickness_end"


!#######################################################################
! <FUNCTION NAME="REMAP_ZT_TO_ZU">
!
! <DESCRIPTION>
! REMAP_ZT_TO_ZU remaps a T-cell thickness or vertical depth on a 
! T-cell to the corresponding place on U-cell.
!
! This is the same operator as REMAP_BT_TO_BU.
! It is needed for ocean_thickness_init since this routine is called
! prior to ocean_operators_init.  This is a bit awkward, but 
! ocean_operators needs thickness values, so it must be called 
! after thickness is initialized. 
! </DESCRIPTION>
!
function REMAP_ZT_TO_ZU(a,Grid) 

    real, dimension(isd:,jsd:), intent(in) :: a
    type(ocean_grid_type),     intent(in)  :: Grid
  
    real, dimension(isd:ied,jsd:jed)       :: REMAP_ZT_TO_ZU
    integer :: i, j

    do j=jsc-halo,jec+halo-1
       do i=isc-halo,iec+halo-1
          REMAP_ZT_TO_ZU(i,j) = (a(i,j)    *Grid%dte(i,j)    *Grid%dus(i,j) &
                               + a(i+1,j)  *Grid%dtw(i+1,j)  *Grid%dus(i,j) &
                               + a(i,j+1)  *Grid%dte(i,j+1)  *Grid%dun(i,j) &
                               + a(i+1,j+1)*Grid%dtw(i+1,j+1)*Grid%dun(i,j) &
                                 )*Grid%daur(i,j)
       enddo
    enddo
    REMAP_ZT_TO_ZU(iec+halo,:) = 0.0
    REMAP_ZT_TO_ZU(:,jec+halo) = 0.0

end function REMAP_ZT_TO_ZU
! </FUNCTION> NAME="REMAP_ZT_TO_ZU"


!#######################################################################
! <SUBROUTINE NAME="dzt_dst_update">
!
! <DESCRIPTION>
! Calculate the quantity dzt_dst.  dzt_dst is required by the Lagrangian
! blob framework in order to calculate the grid cell the a blob resides
! in at taup1.  This information is required prior to the the earliest
! time that we can update the total thickness.  As such, we are
! motivated to separate this calculation when use_blobs=.true. 
! When use_blobs=.false. the calculation is conducted when the tcell 
! thickness is updated.

! We also update coordinate increments for GEOPOTENTIAL and PRESSURE. 
! Only for these two coordinates do we need to modify the endpoint 
! values for the s-grid increments.  Other MOM coordinates have dst 
! and dswt constant in time.  
!
! Note that at present, GEOPOTENTIAL coordinates are not supported
! by this implementation of the blob framework.
! </DESCRIPTION>
!
subroutine dzt_dst_update(Time, Grid, Ext_mode, Dens, Thickness)
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_thickness_type),     intent(inout) :: Thickness

  integer :: i,j,k,kb
  integer :: tau, taup1

  tau   = Time%tau
  taup1 = Time%taup1

  call diagnose_3d(Time, Grid, id_dst, Thickness%dst(:,:,:))
  call diagnose_3d(Time, Grid, id_dzt_dst, Thickness%dzt_dst(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_st, Thickness%depth_st(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_swt, Thickness%depth_swt(:,:,:))

  ! Note: ZSIGMA and PSIGMA not supported by this implementation of the Lagrangian blobs

  if (vert_coordinate==GEOPOTENTIAL) then
     do j=jsc,jed
        do i=isd,ied
           Thickness%dswt(i,j,0)  = Thickness%depth_zt(i,j,1)  + Ext_mode%eta_t(i,j,taup1)
           Thickness%dst(i,j,1)   = Thickness%depth_zwt(i,j,1) + Ext_mode%eta_t(i,j,taup1)
           Thickness%dstup(i,j,1) = Thickness%dswt(i,j,0)
           Thickness%dstlo(i,j,1) = Thickness%dst(i,j,1) - Thickness%dstup(i,j,1)
        enddo
     enddo
     if(linear_free_surface) then  
        do j=jsd,jed
           do i=isd,ied
              Thickness%dswt(i,j,0)  = Thickness%depth_zt(i,j,1)  
              Thickness%dst(i,j,1)   = Thickness%depth_zwt(i,j,1) 
              Thickness%dstup(i,j,1) = Thickness%dswt(i,j,0)
              Thickness%dstlo(i,j,1) = Thickness%dst(i,j,1) - Thickness%dstup(i,j,1)
           enddo
        enddo
     endif
     ! dzt_dst does not change for GEOPOTENTIAL

  elseif (vert_coordinate==ZSTAR) then
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%dzt_dst(i,j,k) = 1.0 + Ext_mode%eta_t(i,j,taup1)*Grid%htr(i,j)
           enddo
        enddo
     enddo

  elseif (vert_coordinate==PRESSURE) then
     do j=jsd,jed
        do i=isd,ied
           Thickness%dswt(i,j,0)     = -Thickness%depth_st(i,j,1)  + Ext_mode%patm_t(i,j,taup1)
           Thickness%dst(i,j,1)      = -Thickness%depth_swt(i,j,1) + Ext_mode%patm_t(i,j,taup1)
           Thickness%dstup(i,j,1)    =  Thickness%dswt(i,j,0)
           Thickness%dstlo(i,j,1)    =  Thickness%dst(i,j,1) - Thickness%dstup(i,j,1)
           kb                        =  Grid%kmt(i,j)
           if(kb > 1) then 
              Thickness%dswt(i,j,kb)  =  Thickness%depth_st(i,j,kb)    - Ext_mode%pbot_t(i,j,taup1)
              Thickness%dst(i,j,kb)   =  Thickness%depth_swt(i,j,kb-1) - Ext_mode%pbot_t(i,j,taup1)
              Thickness%dstlo(i,j,kb) =  Thickness%dswt(i,j,kb)
              Thickness%dstup(i,j,kb) =  Thickness%dst(i,j,kb) - Thickness%dstlo(i,j,kb)
           endif
        enddo
     enddo
     
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              if(Grid%tmask(i,j,k) > 0.0) then 
                 Thickness%dzt_dst(i,j,k) = -grav_r/(Dens%rhoT(i,j,k)+epsln)
              endif
           enddo
        enddo
     enddo

  elseif (vert_coordinate==PSTAR) then
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
               if(Grid%tmask(i,j,k) > 0.0) then 
                   Thickness%dzt_dst(i,j,k) = -grav_r/(Dens%rhoT(i,j,k)+epsln)                          &
                                               *(Ext_mode%pbot_t(i,j,taup1)-Ext_mode%patm_t(i,j,taup1)) &
                                               *Thickness%pbot0r(i,j)
                endif
           enddo
        enddo
     enddo
  endif

end subroutine dzt_dst_update
! </SUBROUTINE> NAME="dzt_dst_update"


!#######################################################################
! <SUBROUTINE NAME="update_tcell_thick_blob">
!
! <DESCRIPTION>
! The same principle is used as for update_tcell_thickness, however,
! this routine is specific to the Lagrangian blob scheme.
!
! The routine calculates the total tracer grid cell thickness for the
! various vertical coordinate systems (excluding ZSGIMA and PSIGMA).
! Major things to note are:
! 1/ dzt_dst is calculated previously using dzt_dst_update
! 2/ this routine does not calculate the L or E system thicknesses,
!    those are calculated separately.
! 3/ This routine is not called unless use_blobs=.true. in
!    ocean_model_nml.
!
! Note that the present implementation of the Lagrangian blobs does
! not support ZSIGMA or ZPRESSURE coordinates.  GEOPOTENTIAL is also
! not supported.
!
! Also note that we calculate the total density weighted thickness,
! rho_dztT, is calculated differently to how the density weighted
! thickness is calculated without the blobs, rho_dzt.
! rho_dztT(taup1) = rho_dztT(tau) + dtime*rho_dzt_tendency
! </DESCRIPTION>
!
subroutine update_tcell_thick_blob(Time, Grid, Ext_mode, Dens, Thickness)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid  
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(in)    :: Dens 
  type(ocean_thickness_type),     intent(inout) :: Thickness

  integer :: i, j, k, kb
  integer :: tau, taup1
  integer :: num_bad=0
  logical :: huge_undulations
  logical :: error_flag=.false.
  character(len=128) :: errorstring

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1
  huge_undulations=.false.

  ! write tau values of the vertical grid increments before update to taup1 
  
  call diagnose_3d(Time, Grid, id_dzt, Thickness%dztT(:,:,:,tau))
  call diagnose_3d(Time, Grid, id_dztL, Thickness%dztL(:,:,:))
  call diagnose_3d(Time, Grid, id_dztE, Thickness%dzt(:,:,:))

  if (id_dzt_pbc > 0) then 
      wrk1_2d(:,:) = 0.0
      do j=jsd,jed
         do i=isd,ied
            kb = Grid%kmt(i,j)
            if (kb > 1) then
                wrk1_2d(i,j) = Thickness%dztT(i,j,kb,tau)
            endif
         enddo
      enddo
      call diagnose_2d(Time, Grid, id_dzt_pbc, wrk1_2d(:,:))
  endif

  call diagnose_3d(Time, Grid, id_dztlo, Thickness%dztloT(:,:,:))
  call diagnose_3d(Time, Grid, id_dztloL, Thickness%dztloL(:,:,:))
  call diagnose_3d(Time, Grid, id_dztloE, Thickness%dztlo(:,:,:))
  
  call diagnose_3d(Time, Grid, id_dztup, Thickness%dztupT(:,:,:))
  call diagnose_3d(Time, Grid, id_dztupL, Thickness%dztupL(:,:,:))
  call diagnose_3d(Time, Grid, id_dztupE, Thickness%dztup(:,:,:))

  call diagnose_3d(Time, Grid, id_dzwt, Thickness%dzwtT(:,:,1:nk))
  call diagnose_3d(Time, Grid, id_dzwtL, Thickness%dzwtL(:,:,1:nk))
  call diagnose_3d(Time, Grid, id_dzwtE, Thickness%dzwt(:,:,1:nk))

  call diagnose_3d(Time, Grid, id_rho_dzt, Thickness%rho_dztT(:,:,:,tau))
  call diagnose_3d(Time, Grid, id_rho_dztL, Thickness%rho_dztL(:,:,:,tau))
  call diagnose_3d(Time, Grid, id_rho_dztE, Thickness%rho_dzt(:,:,:,tau))

  call diagnose_3d(Time, Grid, id_depth_zt, Thickness%depth_zt(:,:,:))
  call diagnose_3d(Time, Grid, id_geodepth_zt, Thickness%geodepth_zt(:,:,:))
  call diagnose_3d(Time, Grid, id_geodepth_zwt, Thickness%geodepth_zwt(:,:,:))
  call diagnose_3d(Time, Grid, id_depth_zwt, Thickness%depth_zwt(:,:,:))

  ! Note, the coordinate increments for GEOPOTENTIAL and PRESSURE
  ! have been updated in dzt_dst_update

  ! update specific thicknesses and vertical increments 
  if(vert_coordinate==GEOPOTENTIAL) then 
      do j=jsd,jed
         do i=isd,ied 
            Thickness%dztloT(i,j,1)     = Thickness%dstlo(i,j,1)*Thickness%dzt_dst(i,j,1)
            Thickness%dztupT(i,j,1)     = Thickness%dstup(i,j,1)*Thickness%dzt_dst(i,j,1)
            Thickness%dztT(i,j,1,taup1) = Thickness%dztloT(i,j,1) + Thickness%dztupT(i,j,1) 
            Thickness%rho_dztT(i,j,1,taup1) = Thickness%rho_dztT(i,j,1,tau) + dtime*Thickness%rho_dzt_tendency(i,j,1)
            wrk1(i,j,:) = 1.0  
         enddo
      enddo
   elseif(vert_coordinate==ZSTAR) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%dztloT(i,j,k)     = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztupT(i,j,k)     = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
               Thickness%dztT(i,j,k,taup1) = Thickness%dztloT(i,j,k) + Thickness%dztupT(i,j,k) 
               Thickness%rho_dztT(i,j,k,taup1) = Thickness%rho_dztT(i,j,k,tau) + dtime*Thickness%rho_dzt_tendency(i,j,k)
!!$               Thickness%rho_dztT(i,j,k,taup1) = rho0*Thickness%dztT(i,j,k,taup1)
               wrk1(i,j,k) = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PRESSURE .or. vert_coordinate==PSTAR) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               ! leave values over land untouched 
               if(Grid%tmask(i,j,k) > 0.0) then 
                  Thickness%dztloT(i,j,k)     = Thickness%dstlo(i,j,k)*Thickness%dzt_dst(i,j,k)
                  Thickness%dztupT(i,j,k)     = Thickness%dstup(i,j,k)*Thickness%dzt_dst(i,j,k)
                  Thickness%dztT(i,j,k,taup1) = Thickness%dztloT(i,j,k) + Thickness%dztupT(i,j,k) 
                  Thickness%rho_dztT(i,j,k,taup1) = Thickness%rho_dztT(i,j,k,tau) + dtime*Thickness%rho_dzt_tendency(i,j,k)
               endif     
               wrk1(i,j,k) = 1.0/Thickness%dzt_dst(i,j,k)  
            enddo
         enddo
      enddo
  endif


  ! vertical distance (m) between t-cell points 
  if(Thickness%method==ENERGETIC) then 
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzwtT(i,j,k) = 2.0*Thickness%dswt(i,j,k)/(wrk1(i,j,k)+wrk1(i,j,k+1)) 
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwtT(i,j,0) = Thickness%dswt(i,j,0)*Thickness%dzt_dst(i,j,1)

            kb=Grid%kmt(i,j)
            if(kb > 0) then 
               Thickness%dzwtT(i,j,kb) = Thickness%dswt(i,j,kb)*Thickness%dzt_dst(i,j,kb)
            endif
         enddo
      enddo
  elseif(Thickness%method==FINITEVOLUME) then 
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               Thickness%dzwtT(i,j,k) = Thickness%dztloT(i,j,k) + Thickness%dztupT(i,j,k+1) 
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwtT(i,j,0) = Thickness%dztupT(i,j,1)

            kb=Grid%kmt(i,j)
            if(kb > 0) then 
                Thickness%dzwtT(i,j,kb) = Thickness%dztloT(i,j,kb)
            endif
         enddo
      enddo
  endif

  ! vertical depth (m) of t and w points as 
  ! computed via vertical sum of z-increments 

  ! compute depth from z=eta to t-point and to bottom of t-cell 
  do j=jsd,jed
     do i=isd,ied
        Thickness%depth_zt(i,j,1)  = Thickness%dzwtT(i,j,0)
        Thickness%depth_zwt(i,j,1) = Thickness%dztT(i,j,1,taup1)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%depth_zt(i,j,k)  = Thickness%depth_zt(i,j,k-1)  + Thickness%dzwtT(i,j,k-1)
           Thickness%depth_zwt(i,j,k) = Thickness%depth_zwt(i,j,k-1) + Thickness%dztT(i,j,k,taup1)
        enddo
     enddo
  enddo

  ! Note: geodepth_zt and geodepth_zwt are not updated here because we do not yet 
  ! know eta_t for pressure-like coordinates.  Thus, geodepth_zt and geodepth_zwt 
  ! are calculated in eta_and_pbot_diagnose in the barotropic module

  ! check for negative specific thickness, which occurs when eta_t is too negative, 
  ! in which case the free surface is deeper than the ocean bottom topography. 
  do j=jsc,jec
     do i=isc,iec
        if(Ext_mode%eta_t(i,j,taup1) < -Grid%ht(i,j)) then 
            write(stdoutunit,*) &
                 '==>Error from ocean_thickness: Surface undulations too negative; model unstable'
            write(stdoutunit,*) &
                 'eta_t(',i+Dom%ioff,',',j+Dom%joff,') = ',Ext_mode%eta_t(i,j,taup1)
            write(stdoutunit,*) &
                 'depth(',i+Dom%ioff,',',j+Dom%joff,') = ',Grid%ht(i,j)        
            huge_undulations=.true.
        endif
     enddo
  enddo
  if(huge_undulations) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_thickness_mod: Free surface penetrating rock! Model unstable.')
  endif
 

  ! check for negative thickness at top and bottom 
  error_flag=.false.
  num_bad=0
  do j=jsc,jec
     do i=isc,iec

        kb = Grid%kmt(i,j)
        if(kb>0 .and. num_bad > max_num_bad_print) then

            if(Thickness%dztT(i,j,1,taup1) < 0.0) then 
                num_bad = num_bad+1

                write(stdoutunit,'(/a,i3,a,i3,a,i3,a,e10.4)')      &
                     '==>Error: negative thickness at top: dztT(',&
                     i+Dom%ioff,',',j+Dom%joff,',',k,') =', Thickness%dztT(i,j,1,taup1)

                if(enforce_positive_dzt) then
                    write(stdoutunit,'(a)') 'Resetting Thickness%dztT to ', thickness_dzt_min
                    Thickness%dztT(i,j,1,taup1) = thickness_dzt_min
                else
                    write(errorstring,'(a,i3,a,i3,a,e10.4)') &
                    '(',i+Dom%ioff,',',j+Dom%joff,',1) =', Thickness%dztT(i,j,1,taup1)
                    call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, i, j, kb, &
                    'From update_tcell_thickness with dzt<0', unit)
                    error_flag=.true.
                endif

            elseif(Thickness%dztT(i,j,kb,taup1) < 0.0) then 
                num_bad = num_bad+1

                write(stdoutunit,'(/a,i3,a,i3,a,i3,a,e10.4)')    &
                '==>Error: negative thickness at bottom: dztT(',&
                i+Dom%ioff,',',j+Dom%joff,',',kb,') =', Thickness%dztT(i,j,kb,taup1)

                if(enforce_positive_dzt) then
                    write(stdoutunit,'(a)') 'Resetting Thickness%dztT to ', thickness_dzt_min
                    Thickness%dztT(i,j,kb,taup1) = thickness_dzt_min
                else
                    write(errorstring,'(a,i3,a,i3,a,i3,a,e10.4)') &
                    '(',i+Dom%ioff,',',j+Dom%joff,',',kb,') =', Thickness%dztT(i,j,kb,taup1)
                    call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, i, j, kb, &
                    'From update_tcell_thickness with dztT<0', unit)
                    error_flag=.true.
                endif

            endif

        endif

     enddo
  enddo

  if(error_flag) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_thickness_mod: dztT'//trim(errorstring)//'<0. grep for "dzt" to find location')
  endif

end subroutine update_tcell_thick_blob
! </SUBROUTINE> NAME="update_tcell_thick_blob"


!#######################################################################
! <SUBROUTINE NAME="update_E_thickness">
!
! <DESCRIPTION>
! Calculates the E system thickness using the basic formulation,
! E_thickness = Total_thickness - L_thickness
!
! Note that this routine is called twice each time step.  Once from 
! update_ocean_model, and once from the update_ocean_tracer.  The 
! second call is required after the implicit adjustment involving
! blobs.
! </DESCRIPTION>
!
subroutine update_E_thickness(Time, Grid, Thickness, Dens, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_grid_type),          intent(in)    :: Grid
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(in)    :: Ext_mode

  integer :: i,j,k
  integer :: tau, taup1
  integer :: stdoutunit

  if (.not. use_blobs) return
  
  tau   = Time%tau
  taup1 = Time%taup1
  stdoutunit = stdout()

  if (vert_coordinate_class==DEPTH_BASED) then
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              Thickness%dztlo(i,j,k) = Thickness%dztloT(i,j,k)     - Thickness%dztloL(i,j,k)
              Thickness%dztup(i,j,k) = Thickness%dztupT(i,j,k)     - Thickness%dztupL(i,j,k)
              Thickness%dzt(i,j,k)   = Thickness%dztT(i,j,k,taup1) - Thickness%dztL(i,j,k)
              
              Thickness%rho_dzt(i,j,k,taup1) = Thickness%rho_dztT(i,j,k,taup1) - Thickness%rho_dztL(i,j,k,taup1)
              Thickness%rho_dztr(i,j,k)      = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
           enddo
        enddo
     enddo
  else !PRESSURE_BASED
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               ! leave values over land untouched 
               if(Grid%tmask(i,j,k) > 0.0) then 
                  Thickness%dztlo(i,j,k) = Thickness%dztloT(i,j,k)     - Thickness%dztloL(i,j,k)
                  Thickness%dztup(i,j,k) = Thickness%dztupT(i,j,k)     - Thickness%dztupL(i,j,k)
                  Thickness%dzt(i,j,k)   = Thickness%dztT(i,j,k,taup1) - Thickness%dztL(i,j,k)
                  
                  Thickness%rho_dzt(i,j,k,taup1) = Thickness%rho_dztT(i,j,k,taup1) - Thickness%rho_dztL(i,j,k,taup1)
                  Thickness%rho_dztr(i,j,k)       = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
               endif 
            enddo
         enddo
      enddo
   endif !vert_coordinate_class==DEPTH_BASED

  ! vertical distance (m) between t-cell points 
   do k=0,nk
      do j=jsd,jed
         do i=isd,ied
            Thickness%dzwt(i,j,k)  = Thickness%dzwtT(i,j,k) - Thickness%dzwtL(i,j,k)
         enddo
      enddo
   enddo

   if (debug_this_module) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               if (Thickness%dztT(i,j,k,taup1) < Thickness%dztL(i,j,k)) then
                  write (*,'(/,a)') 'Lagrangian thickness greater than total thickness'
                  print ('(3(a1,i3),a1,x,2(a5,es14.7,x))'), '(',i,',',j,',',k,')', &
                       'dztT=', Thickness%dztT(i,j,k,taup1), 'dztL=', Thickness%dztL(i,j,k)
                  call mpp_error(FATAL, &
                       '==>Error from ocean_thickness_mod (update_E_thickness): Lagrangian thickness greater than total thickness')
               endif
            enddo
         enddo
      enddo
   endif

   if(debug_this_module_detail) then 
      call thickness_details(Grid, Time, Thickness, Ext_mode, Dens, isc, jsc, nk, &
                             'From update_tcell_thickness', stdoutunit)
      call dzt_min_max(Time, Thickness, &
           'From update_E_thickness with dzt<0, dzt_min_max information')
   endif

end subroutine update_E_thickness
! </SUBROUTINE> NAME="update_E_thickness"


!#######################################################################
! <SUBROUTINE NAME="thickness_chksum_blobs">
!
! <DESCRIPTION>
! Compute checksum for thickness components for the Lagrangian and
! total thicknesses.
!
! Only print checksums for fields that should agree across restarts.  
!
! </DESCRIPTION>
subroutine thickness_chksum_blobs(Time, Grid, Thickness)
 
  type(ocean_time_type),      intent(in) :: Time
  type(ocean_grid_type),      intent(in) :: Grid 
  type(ocean_thickness_type), intent(in) :: Thickness

  integer :: tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  call write_chksum_3d('Thickness%rho_dztT(taup1)', Thickness%rho_dztT(COMP,:,taup1)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzt(taup1)', Thickness%rho_dzt(COMP,:,taup1)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dztL(taup1)', Thickness%rho_dztL(COMP,:,taup1)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzuT(taup1)', Thickness%rho_dzuT(COMP,:,taup1)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzu(taup1)', Thickness%rho_dzu(COMP,:,taup1)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzuL(taup1)', Thickness%rho_dzuL(COMP,:,taup1)*Grid%umask(COMP,:))
  call write_chksum_2d('Thickness%mass_uT(taup1)', Thickness%mass_uT(COMP,taup1)*Grid%umask(COMP,1))
  call write_chksum_2d('Thickness%mass_u(taup1)', Thickness%mass_u(COMP,taup1)*Grid%umask(COMP,1))
  call write_chksum_3d('Thickness%rho_dztr', Thickness%rho_dztr(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzur', Thickness%rho_dzur(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%rho_dzt_tendency', Thickness%rho_dzt_tendency(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztT', Thickness%dztT(COMP,:,taup1)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dzt', Thickness%dzt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztL', Thickness%dztL(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztloT', Thickness%dztloT(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztlo', Thickness%dztlo(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztloL', Thickness%dztloL(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztupT', Thickness%dztupT(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztup', Thickness%dztup(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dztupL', Thickness%dztupL(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dzt_dst', Thickness%dzt_dst(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_2d('Thickness%dzwtT(k=0)', Thickness%dzwtT(COMP,0)*Grid%tmask(COMP,1))
  call write_chksum_2d('Thickness%dzwt(k=0)', Thickness%dzwt(COMP,0)*Grid%tmask(COMP,1))
  call write_chksum_2d('Thickness%dzwtL(k=0)', Thickness%dzwtL(COMP,0)*Grid%tmask(COMP,1))
  call write_chksum_3d('Thickness%dzwtT(k=1:nk)', Thickness%dzwtT(COMP,1:nk)*Grid%tmask(COMP,1:nk))
  call write_chksum_3d('Thickness%dzwt(k=1:nk)', Thickness%dzwt(COMP,1:nk)*Grid%tmask(COMP,1:nk))
  call write_chksum_3d('Thickness%dzwtL(k=1:nk)', Thickness%dzwtL(COMP,1:nk)*Grid%tmask(COMP,1:nk))
  call write_chksum_3d('Thickness%dzuT', Thickness%dzuT(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%dzu', Thickness%dzu(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%dzuL', Thickness%dzuL(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_2d('Thickness%dzwuT(k=0)', Thickness%dzwuT(COMP,0)*Grid%umask(COMP,1))
  call write_chksum_2d('Thickness%dzwu(k=0)', Thickness%dzwu(COMP,0)*Grid%umask(COMP,1))
  call write_chksum_2d('Thickness%dzwuL(k=0)', Thickness%dzwuL(COMP,0)*Grid%umask(COMP,1))
  call write_chksum_3d('Thickness%dzwuT(k=1:nk)', Thickness%dzwuT(COMP,1:nk)*Grid%umask(COMP,1:nk))
  call write_chksum_3d('Thickness%dzwu(k=1:nk)', Thickness%dzwu(COMP,1:nk)*Grid%umask(COMP,1:nk))
  call write_chksum_3d('Thickness%dzwuL(k=1:nk)', Thickness%dzwuL(COMP,1:nk)*Grid%umask(COMP,1:nk))
  call write_chksum_3d('Thickness%depth_zt', Thickness%depth_zt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_zwt', Thickness%depth_zwt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%geodepth_zt', Thickness%geodepth_zt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%geodepth_zwt', Thickness%geodepth_zwt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_zwt', Thickness%depth_zwt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_zu', Thickness%depth_zu(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%depth_zwu', Thickness%depth_zwu(COMP,:)*Grid%umask(COMP,:))
  call write_chksum_3d('Thickness%depth_st', Thickness%depth_st(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%depth_swt', Thickness%depth_swt(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dst', Thickness%dst(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dstlo', Thickness%dstlo(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_3d('Thickness%dstup', Thickness%dstup(COMP,:)*Grid%tmask(COMP,:))
  call write_chksum_2d('Thickness%dswt(k=0)', Thickness%dswt(COMP,0)*Grid%tmask(COMP,1))
  call write_chksum_3d('Thickness%dswt(k=1:nk)', Thickness%dswt(COMP,1:nk)*Grid%tmask(COMP,1:nk))
  call write_chksum_2d('Thickness%pbot0', Thickness%pbot0(COMP)*Grid%tmask(COMP,1))
  
  write(stdoutunit,*) ' '
  
end subroutine thickness_chksum_blobs
! </SUBROUTINE> NAME="thickness_chksum_blobs"


end module ocean_thickness_mod
