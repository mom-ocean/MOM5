module ocean_tracer_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater (initialization)
!</CONTACT>
!
!<OVERVIEW>
! This module time steps the tracer fields.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module time steps the tracer fields.
! Initialization for the tracer packages is done as well. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2004)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Fundamentals of ocean climate models (2004)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_tracer_nml">
!
!  <DATA NAME="zero_tendency" TYPE="logical">
!  If true, then will freeze the tracer fields.
!  </DATA> 
!
!  <DATA NAME="zero_tracer_source" TYPE="logical">
!  To remove the T_prog%source contribution to tracer
!  evolution. For debugging purposes.  Default 
!  zero_tracer_source=.false.
!  </DATA> 
!
!  <DATA NAME="limit_age_tracer" TYPE="logical">
!  Limit the values of age tracer to be less than 
!  total run time and greater than zero. 
!  Default limit_age_tracer=.false.
!  </DATA> 
!
!  <DATA NAME="age_tracer_max_init" TYPE="real" UNITS="years">
!  Initial maximum age tracer. This nml provides the ability to 
!  start an integration with an age tracer that is not initialized
!  to zero, say if we took an age tracer from another spin-up. 
!  Default age_tracer_max_init=0.0.
!  </DATA> 
!
!  <DATA NAME="remap_depth_to_s_init" TYPE="logical">
!  For remapping initial tracer distributions, generally determined 
!  according to depth vertical coordinates using the mom preprocessing
!  schemes, onto s-coordinates.  This method is of use for initializing
!  terrain following coordinate simulations with mom.  
!  </DATA> 
!
!  <DATA NAME="frazil_heating_before_vphysics" TYPE="logical">
!  For computing frazil heating before the implicit vertical physics
!  (which includes boundary fluxes), and before vertical convection. 
!  This is the order that CM2.0 and CM2.1 performed their calculations
!  of frazil.  It is arguable that one should NOT do frazil until the
!  end of a time step, after vertical physics and after surface 
!  boundary fluxes. 
!  Default frazil_heating_before_vphysics=.false.
!  </DATA> 
!
!  <DATA NAME="frazil_heating_after_vphysics" TYPE="logical">
!  For computing frazil heating after the implicit vertical physics
!  (which includes boundary fluxes), and after vertical convection. 
!  This is the recommended method. 
!  Default frazil_heating_after_vphysics=.false.
!  </DATA> 
!
!  <DATA NAME="tmask_limit_ts_same" TYPE="logical">
!  tmask_limit is derived separately for the tracers.  However,
!  it may be appropriate to have the mask be the same for temp 
!  and salinity, in which case the neutral physics fluxes are 
!  self-consistent.  But for some cases, such as when running with 
!  linear eos, may not wish to have the temp and salinity coupled
!  when computing the mask.  
!  </DATA> 
!
!  <DATA NAME="compute_tmask_limit_on" TYPE="logical">
!  For updating the tmaks_limit array. This calculation is 
!  recommended for the following physics and advection schemes:
!  1/ quicker advection
!  2/ neutral physics
!  3/ submesoscale closure.
!  The default is compute_tmask_limit_on=.true., but if none
!  of the above schemes is used, then some time savings can be
!  realized by setting compute_tmask_limit_on=.false.
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging the tracer module
!  </DATA> 
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="interpolate_tprog_to_pbott" TYPE="logical">
!  To linear interpolate the initial conditions for prognostic
!  tracers to the partial bottom cells.  Default 
!  interpolate_tprog_to_pbott=.true. 
!  </DATA> 
!
!  <DATA NAME="interpolate_tdiag_to_pbott" TYPE="logical">
!  To linear interpolate the initial conditions for diagnostic 
!  tracers to the partial bottom cells.  Default 
!  interpolate_tdiag_to_pbott=.false. 
!  </DATA> 
!
!  <DATA NAME="inflow_nboundary" TYPE="logical">
!  For adding an inflow transport from the northern boundary
!  which brings in temp and salinity according to inflow data
!  files. Default is inflow_nboundary=.false.
!  </DATA> 
!
!  <DATA NAME="ocean_tpm_debug" TYPE="logical">
!  For debugging ocean tracer package manager.  
!  </DATA>
!
!  <DATA NAME="use_tempsalt_check_range" TYPE="logical">
!  To call a check to see that temperature and salinity
!  are within their pre-selected range.
!  Default use_tempsalt_check_range=.false. since this
!  check may incur some cost that users should be aware of. 
!  </DATA>

  
!</NAMELIST>
!
use constants_mod,     only: epsln, kelvin
use diag_manager_mod,  only: register_diag_field, send_data
use field_manager_mod, only: fm_string_len, fm_type_name_len
use field_manager_mod, only: fm_get_length
use field_manager_mod, only: fm_dump_list, fm_change_list, fm_loop_over_list
use fms_mod,           only: read_data, file_exist, field_exist
use fms_mod,           only: open_namelist_file, check_nml_error, close_file
use fms_mod,           only: FATAL, WARNING, NOTE, stdout, stdlog
use fms_io_mod,        only: register_restart_field, save_restart, restore_state
use fms_io_mod,        only: restart_file_type, reset_field_pointer, reset_field_name 
use fms_io_mod,        only: field_size
use fm_util_mod,       only: fm_util_get_real, fm_util_get_integer
use fm_util_mod,       only: fm_util_get_logical, fm_util_get_string
use fm_util_mod,       only: fm_util_get_string_array
use fm_util_mod,       only: fm_util_check_for_bad_fields, fm_util_set_value
use mpp_domains_mod,   only: mpp_update_domains
use mpp_domains_mod,   only: mpp_global_sum, NON_BITWISE_EXACT_SUM
use mpp_io_mod,        only: mpp_open, fieldtype
use mpp_io_mod,        only: MPP_NETCDF, MPP_OVERWR, MPP_ASCII, MPP_RDONLY, MPP_SINGLE, MPP_MULTI 
use mpp_mod,           only: input_nml_file, mpp_error, mpp_pe, mpp_root_pe, mpp_broadcast, ALL_PES, mpp_max, mpp_min, mpp_chksum 
use mpp_mod,           only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,           only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE, CLOCK_ROUTINE
use platform_mod,      only: i8_kind
use time_manager_mod,  only: time_type, set_time, increment_time, operator( + )

use transport_matrix_mod, only: do_transport_matrix, transport_matrix_store_explicit

use ocean_blob_mod,             only: ocean_blob_implicit, adjust_L_thickness
use ocean_convect_mod,          only: convection
use ocean_domains_mod,          only: get_local_indices
use ocean_density_mod,          only: neutral_density, update_ocean_density_salinity
use ocean_frazil_mod,           only: compute_frazil_heating
use ocean_bih_tracer_mod,       only: bih_tracer 

#ifdef USE_OCEAN_BGC
use ocean_generic_mod,          only: ocean_generic_get_field, ocean_generic_get_field_pointer
use ocean_generic_mod,          only: ocean_generic_set_pointer
#endif

use ocean_lap_tracer_mod,       only: lap_tracer 
use ocean_obc_mod,              only: ocean_obc_tracer, ocean_obc_update_boundary, ocean_obc_tracer_init
use ocean_parameters_mod,       only: ADVECT_UPWIND, ADVECT_2ND_ORDER, ADVECT_4TH_ORDER, ADVECT_6TH_ORDER
use ocean_parameters_mod,       only: ADVECT_QUICKER, ADVECT_QUICKMOM3, ADVECT_MDFL_SUP_B, ADVECT_MDFL_SWEBY
use ocean_parameters_mod,       only: ADVECT_PSOM, ADVECT_MDPPM, ADVECT_DST_LINEAR, ADVECT_MDMDT_TEST
use ocean_parameters_mod,       only: ADVECT_MDPPM_TEST, ADVECT_MDFL_SWEBY_TEST, ADVECT_DST_LINEAR_TEST
use ocean_parameters_mod,       only: missing_value, sec_in_yr_r, rho0, rho0r, cp_ocean
use ocean_parameters_mod,       only: TWO_LEVEL, THREE_LEVEL
use ocean_parameters_mod,       only: CONSERVATIVE_TEMP, POTENTIAL_TEMP
use ocean_parameters_mod,       only: QUASI_HORIZONTAL, TERRAIN_FOLLOWING
use ocean_parameters_mod,       only: GEOPOTENTIAL
use ocean_parameters_mod,       only: DEPTH_BASED
use ocean_passive_mod,          only: passive_tracer_init, update_tracer_passive 
use ocean_shortwave_mod,        only: ocean_irradiance_init
use ocean_tempsalt_mod,         only: contemp_from_pottemp, pottemp_from_contemp, tempsalt_check_range
use ocean_tpm_mod,              only: ocean_tpm_init
use ocean_tpm_util_mod,         only: otpm_set_tracer_package
use ocean_tracer_advect_mod,    only: horz_advect_tracer, vert_advect_tracer
use ocean_tracer_diag_mod,      only: send_tracer_variance
use ocean_tracer_util_mod,      only: rebin_onto_rho, diagnose_mass_of_layer 
use ocean_tracer_util_mod,      only: tracer_prog_chksum, tracer_diag_chksum, tracer_min_max
use ocean_thickness_mod,        only: update_E_thickness
use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
use ocean_types_mod,            only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,            only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,            only: ocean_adv_vel_type, ocean_external_mode_type
use ocean_types_mod,            only: ocean_public_type, ocean_density_type, ocean_options_type
use ocean_types_mod,            only: ocean_lagrangian_type, ocean_velocity_type, blob_diag_type
use ocean_util_mod,             only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_sum, write_chksum_3d
use ocean_tracer_util_mod,      only: diagnose_3d_rho
use ocean_vert_mix_mod,         only: vert_diffuse, vert_diffuse_implicit
use ocean_workspace_mod,        only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, wrk1_2d

implicit none

private 

logical :: prog_module_initialized = .false.
logical :: diag_module_initialized = .false.

character(len=256) :: version='CVS $Id: ocean_tracer.F90,v 20.0 2013/12/14 00:17:20 fms Exp $'
character(len=256) :: tagname='Tag $Name: tikal $'
character(len=48), parameter          :: mod_name = 'ocean_tracer_mod'

integer :: num_tracers       =0
integer :: num_prog_tracers  =0
integer :: num_diag_tracers  =0
integer :: num_family_tracers=0

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

! for time steps and implicit vertical mixing 
integer :: tendency = 0
real    :: dtime    = 0.0
real    :: dtimer   = 0.0
real    :: dtts     = 0.0
real    :: dtuv     = 0.0
real    :: cp_oceanr

integer :: index_temp        =-1
integer :: index_salt        =-1  
integer :: index_frazil      =-1  
integer :: index_temp_sq     =-1
integer :: index_salt_sq     =-1  
integer :: index_diag_temp   =-1
integer :: index_added_heat  =-1
integer :: index_redist_heat =-1

! for obc 
logical :: have_obc=.false.

! for global normalization
real    :: cellarea_r

! for vertical coordinate class 
integer :: vert_coordinate_class 

! for setting the temperature variable
integer :: prog_temp_variable=0

! for CMIP units (deg K rather than degC)
real :: cmip_offset = 0.0
character(len=32) :: temp_units='degrees C' 


! for possible inflow tracer settings 
real, dimension(:,:,:), allocatable :: mask_inflow 
real, dimension(:,:,:), allocatable :: salt_inflow 
real, dimension(:,:,:), allocatable :: temp_inflow 
  
! identification numbers for mpp clocks
integer :: id_clock_bih_tracer
integer :: id_clock_lap_tracer
integer :: id_clock_vert_diffuse
integer :: id_clock_vert_diffuse_implicit
integer :: id_clock_tracer_advect
integer :: id_clock_frazil
integer :: id_clock_convection
integer :: id_clock_blob
integer :: id_clock_blob_implicit
integer :: id_clock_adjust_L_thickness
integer :: id_clock_adjust_E_thickness

! for restart
integer                              :: num_tracer_restart = 0
integer,                 allocatable :: id_tracer_restart(:)
integer,                 allocatable :: id_tracer_type(:)
character(len=64),       allocatable :: tracer_restart_file(:)
type(restart_file_type), allocatable :: Tracer_restart(:)

!work array on neutral density space
integer :: neutralrho_nk
real, dimension(:,:,:),   allocatable :: nrho_work 
real, dimension(:,:,:),   allocatable :: nrho_work2 
real, dimension(:,:,:),   allocatable :: nrho_mass

! internally set for computing watermass diagnstics
logical :: compute_watermass_diag = .false. 


! for diagnostics 
logical :: used
integer, allocatable, dimension(:) :: id_eta_smooth
integer, allocatable, dimension(:) :: id_eta_smooth_on_nrho
integer, allocatable, dimension(:) :: id_pbot_smooth
integer, allocatable, dimension(:) :: id_prog
integer, allocatable, dimension(:) :: id_progT
integer, allocatable, dimension(:) :: id_prog_explicit
integer, allocatable, dimension(:) :: id_prog_rhodzt
integer, allocatable, dimension(:) :: id_prog_rhodztT
integer, allocatable, dimension(:) :: id_prog_int_rhodz
integer, allocatable, dimension(:) :: id_prog_int_rhodzT
integer, allocatable, dimension(:) :: id_prog_on_depth
integer, allocatable, dimension(:) :: id_tendency_conc
integer, allocatable, dimension(:) :: id_tendency_concL
integer, allocatable, dimension(:) :: id_tendency_concT
integer, allocatable, dimension(:) :: id_tendency
integer, allocatable, dimension(:) :: id_tendency_on_nrho
integer, allocatable, dimension(:) :: id_tendencyL
integer, allocatable, dimension(:) :: id_tendencyT
integer, allocatable, dimension(:) :: id_tendency_expl
integer, allocatable, dimension(:) :: id_tendency_explL
integer, allocatable, dimension(:) :: id_tendency_explT
integer, allocatable, dimension(:) :: id_surf_tracer
integer, allocatable, dimension(:) :: id_surf_tracerT
integer, allocatable, dimension(:) :: id_surf_tracer_sq
integer, allocatable, dimension(:) :: id_diag_surf_tracer
integer, allocatable, dimension(:) :: id_diag_surf_tracer_sq
integer, allocatable, dimension(:) :: id_bott_tracer
integer, allocatable, dimension(:) :: id_bott_tracerT
 integer, allocatable, dimension(:) :: id_diag
integer, allocatable, dimension(:) :: id_diag_total
integer, allocatable, dimension(:) :: id_tmask_limit

integer  :: id_neut_rho_tendency 
integer  :: id_neut_rho_tendency_on_nrho
integer  :: id_wdian_rho_tendency 
integer  :: id_wdian_rho_tendency_on_nrho 
integer  :: id_tform_rho_tendency 
integer  :: id_tform_rho_tendency_on_nrho

integer  :: id_neut_temp_tendency 
integer  :: id_neut_temp_tendency_on_nrho
integer  :: id_wdian_temp_tendency 
integer  :: id_wdian_temp_tendency_on_nrho 
integer  :: id_tform_temp_tendency 
integer  :: id_tform_temp_tendency_on_nrho

integer  :: id_neut_salt_tendency 
integer  :: id_neut_salt_tendency_on_nrho
integer  :: id_wdian_salt_tendency 
integer  :: id_wdian_salt_tendency_on_nrho 
integer  :: id_tform_salt_tendency 
integer  :: id_tform_salt_tendency_on_nrho

integer  :: id_neut_rho_smooth
integer  :: id_neut_rho_smooth_on_nrho
integer  :: id_wdian_rho_smooth
integer  :: id_wdian_rho_smooth_on_nrho
integer  :: id_tform_rho_smooth
integer  :: id_tform_rho_smooth_on_nrho
integer  :: id_eta_tend_smooth
integer  :: id_eta_tend_smooth_glob

integer  :: id_neut_temp_smooth
integer  :: id_neut_temp_smooth_on_nrho
integer  :: id_wdian_temp_smooth
integer  :: id_wdian_temp_smooth_on_nrho
integer  :: id_tform_temp_smooth
integer  :: id_tform_temp_smooth_on_nrho

integer  :: id_neut_salt_smooth
integer  :: id_neut_salt_smooth_on_nrho
integer  :: id_wdian_salt_smooth
integer  :: id_wdian_salt_smooth_on_nrho
integer  :: id_tform_salt_smooth
integer  :: id_tform_salt_smooth_on_nrho


integer  :: id_neut_rho_pme
integer  :: id_wdian_rho_pme
integer  :: id_tform_rho_pme
integer  :: id_neut_rho_pme_on_nrho
integer  :: id_wdian_rho_pme_on_nrho
integer  :: id_tform_rho_pme_on_nrho

integer  :: id_neut_rho_pbl_pme_kn
integer  :: id_neut_rho_pbl_pme_kn_on_nrho
integer  :: id_wdian_rho_pbl_pme_kn
integer  :: id_wdian_rho_pbl_pme_kn_on_nrho
integer  :: id_tform_rho_pbl_pme_kn
integer  :: id_tform_rho_pbl_pme_kn_on_nrho

integer  :: id_neut_temp_pbl_pme_kn
integer  :: id_neut_temp_pbl_pme_kn_on_nrho
integer  :: idwdian_temp_pbl_pme_kn
integer  :: idwdian_temp_pbl_pme_kn_on_nrho
integer  :: idtform_temp_pbl_pme_kn
integer  :: idtform_temp_pbl_pme_kn_on_nrho

integer  :: id_neut_salt_pbl_pme_kn
integer  :: id_neut_salt_pbl_pme_kn_on_nrho
integer  :: idwdian_salt_pbl_pme_kn
integer  :: idwdian_salt_pbl_pme_kn_on_nrho
integer  :: idtform_salt_pbl_pme_kn
integer  :: idtform_salt_pbl_pme_kn_on_nrho


integer  :: id_neut_rho_pbl_pme_pr
integer  :: id_neut_rho_pbl_pme_pr_on_nrho
integer  :: id_wdian_rho_pbl_pme_pr
integer  :: id_wdian_rho_pbl_pme_pr_on_nrho
integer  :: id_tform_rho_pbl_pme_pr
integer  :: id_tform_rho_pbl_pme_pr_on_nrho

integer  :: id_neut_temp_pbl_pme_pr
integer  :: id_neut_temp_pbl_pme_pr_on_nrho
integer  :: idwdian_temp_pbl_pme_pr
integer  :: idwdian_temp_pbl_pme_pr_on_nrho
integer  :: idtform_temp_pbl_pme_pr
integer  :: idtform_temp_pbl_pme_pr_on_nrho

integer  :: id_neut_salt_pbl_pme_pr
integer  :: id_neut_salt_pbl_pme_pr_on_nrho
integer  :: idwdian_salt_pbl_pme_pr
integer  :: idwdian_salt_pbl_pme_pr_on_nrho
integer  :: idtform_salt_pbl_pme_pr
integer  :: idtform_salt_pbl_pme_pr_on_nrho

integer  :: id_frazil_on_nrho

integer  :: id_neut_rho_frazil
integer  :: id_wdian_rho_frazil
integer  :: id_tform_rho_frazil
integer  :: id_neut_rho_frazil_on_nrho
integer  :: id_wdian_rho_frazil_on_nrho
integer  :: id_tform_rho_frazil_on_nrho
integer  :: id_eta_tend_frazil
integer  :: id_eta_tend_frazil_glob

integer  :: id_mass_t_on_nrho
integer  :: id_mass_t_tendency_on_nrho
integer  :: id_mass_nrho_layer 
integer  :: id_mass_nrho_tendency_layer

! for ascii output
integer :: unit=6

public  update_ocean_tracer
public  ocean_prog_tracer_init
public  ocean_diag_tracer_init
public  ocean_tracer_diagnostics_init
public  ocean_tracer_end
public  compute_tmask_limit
public  ocean_tracer_restart

private update_advection_only
private remap_s_to_depth
private remap_depth_to_s
private inflow_nboundary_init
private watermass_diag 
private send_tracer_diagnostics


!---------------nml settings---------------

logical :: zero_tendency                   = .false.
logical :: zero_tracer_source              = .false. 
logical :: debug_this_module               = .false.
logical :: tmask_limit_ts_same             = .true.
logical :: write_a_restart                 = .true. 
logical :: remap_depth_to_s_init           = .false. 
logical :: inflow_nboundary                = .false. 
logical :: interpolate_tprog_to_pbott      = .true. 
logical :: interpolate_tdiag_to_pbott      = .false. 
logical :: limit_age_tracer                = .false.
logical :: frazil_heating_before_vphysics  = .false.
logical :: frazil_heating_after_vphysics   = .false.
logical :: compute_tmask_limit_on          = .true.
logical :: use_tempsalt_check_range        = .false. 
real    :: age_tracer_max_init             = 0.0

! for tracer package manager 
logical :: ocean_tpm_debug = .false.

namelist /ocean_tracer_nml/ debug_this_module, zero_tendency, zero_tracer_source, write_a_restart,   &
                            ocean_tpm_debug, tmask_limit_ts_same, remap_depth_to_s_init,             &
                            inflow_nboundary, interpolate_tprog_to_pbott, interpolate_tdiag_to_pbott,&
                            limit_age_tracer, age_tracer_max_init,                                   &
                            frazil_heating_before_vphysics, frazil_heating_after_vphysics,           &
                            compute_tmask_limit_on, use_tempsalt_check_range

contains


!#######################################################################
! <FUNCTION NAME="ocean_prog_tracer_init">
!
! <DESCRIPTION>
! Initialization code for prognostic tracers, returning a pointer to 
! the T_prog array.
! </DESCRIPTION>
!
function ocean_prog_tracer_init (Grid, Thickness, Ocean_options, Domain, Time, Time_steps, &
                                 num_prog, vert_coordinate_type, obc, cmip_units, use_blobs, debug)   &
                                 result (T_prog)  
  
  type(ocean_grid_type),       intent(in), target   :: Grid
  type(ocean_thickness_type),  intent(in)           :: Thickness
  type(ocean_options_type),    intent(inout)        :: Ocean_options
  type(ocean_domain_type),     intent(in), target   :: Domain
  type(ocean_time_type),       intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)           :: Time_steps
  integer,                     intent(out)          :: num_prog
  integer,                     intent(in)           :: vert_coordinate_type 
  logical,                     intent(in)           :: obc
  logical,                     intent(in)           :: cmip_units
  logical,                     intent(in)           :: use_blobs
  logical,                     intent(in), optional :: debug

  ! return value 
  type(ocean_prog_tracer_type), dimension(:), pointer :: T_prog

  integer               :: i, j, k, n, kb, l
  integer               :: ierr, num_diag
  integer               :: tau, taum1, taup1
  integer               :: ioun, io_status
  integer, dimension(4) :: siz 
  integer               :: frazil_heating_order=0
  real                  :: fact
  character(len=32)     :: name
  character(len=128)    :: filename
  logical               :: initialize_as_a_passive_tracer=.false.
  character(len=33)     :: prog_name
  character(len=138)    :: prog_longname

  character(len=48),  parameter :: sub_name = 'ocean_prog_tracer_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  ! variables for tracer package
  integer                               :: ind
  real, dimension(2)                    :: range_array
  character(len=64)                     :: caller_str
  character(len=fm_string_len)          :: string_fm
  character(len=fm_type_name_len)       :: typ
  character(len=fm_string_len), pointer, dimension(:) :: good_list
  integer :: stdoutunit, stdlogunit
  stdoutunit=stdout();stdlogunit=stdlog()

  if (prog_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' Prognostic tracers already initialized')
  endif
  
  nullify(T_prog)

  write( stdlogunit,'(/a/)') trim(version)

  have_obc      = obc
  dtts          = Time_steps%dtts
  dtuv          = Time_steps%dtuv
  cp_oceanr     = 1.0/cp_ocean 

  ! provide for namelist over-ride
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_tracer_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_tracer_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_tracer_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_tracer_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_tracer_nml)  
  write (stdlogunit, ocean_tracer_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  if(cmip_units) then 
      cmip_offset = kelvin
      temp_units='degrees K' 
  else
      cmip_offset = 0.0  
      temp_units='degrees C' 
  endif

  if(zero_tendency) then 
      call mpp_error(NOTE, trim(note_header) // ' zero_tendency=true so will not time step tracer fields.')
      Ocean_options%tracer_tendency = 'Did NOT time step prognostic tracer fields.'      
  else
      Ocean_options%tracer_tendency = 'Time stepped the prognostic tracer fields.'      
  endif 

  if(zero_tracer_source) then 
      call mpp_error(NOTE, &
      trim(note_header) // ' zero_tracer_source=true so remove T_prog%source from evolution.')
  endif 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_tracer with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif 

  if(frazil_heating_before_vphysics) then 
      frazil_heating_order=frazil_heating_order+1
      write(stdoutunit,'(a)') '==>Note: frazil heating called before vertical physics and before boundary fluxes.'
      write(stdoutunit,'(a)') '         This method is retained for legacy purposes: it is NOT recommended for new runs. '
      write(stdoutunit,'(a)') ' '
  endif
  if(frazil_heating_after_vphysics) then 
      frazil_heating_order=frazil_heating_order+1
      write(stdoutunit,'(a)') '==>Note: frazil heating called after vertical physics and after boundary fluxes.'
      write(stdoutunit,'(a)') '         This is the recommended method. '
      write(stdoutunit,'(a)') ' '
  endif
  if(frazil_heating_order>1) then 
      write(stdoutunit,'(a)') '==>Error from ocean_tracer_mod: choose just one temporal order for frazil heating.'
      write(stdoutunit,'(a)') ' '
      call mpp_error(FATAL, &
      trim(note_header) // ' can only choose one temporal order for frazil heating.')
  endif
  if(frazil_heating_order==0) then 
      write(stdoutunit,'(a)') '==>Error from ocean_tracer_mod: MUST specify order for frazil heating in ocean_tracer.'
      write(stdoutunit,'(a)') ' '
      call mpp_error(FATAL, &
      trim(note_header) // '  Need to specify an order for calling frazil heating from ocean_tracer_mod.')
  endif

  ! some time step information 
  if (dtts /= dtuv) then
     write (stdoutunit,'(/a)') trim(warn_header) // ' Asynchronous timesteps (dtts > dtuv) imply inaccurate transients.'
     write (stdoutunit,'(a)')  '          and total tracer (i.e. heat content) is not conserved.'
     write (stdoutunit,'(a,f5.2,a/)') '            dtts =',dtts/dtuv,' times larger than dtuv.'         
  else
     call mpp_error(NOTE, trim(note_header) // ' Synchronous timesteps have been specified (dtts = dtuv).')
  endif

  tendency   = Time_steps%tendency
  dtime      = Time_steps%dtime_t
  dtimer     = 1.0/(dtime+epsln)

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  ! initialize clock ids 
  id_clock_bih_tracer            = mpp_clock_id('(Ocean tracer: bih tracer)   '     ,grain=CLOCK_MODULE)
  id_clock_lap_tracer            = mpp_clock_id('(Ocean tracer: lap tracer)   '     ,grain=CLOCK_MODULE)
  id_clock_vert_diffuse          = mpp_clock_id('(Ocean tracer: vert diffuse) '     ,grain=CLOCK_MODULE)
  id_clock_vert_diffuse_implicit = mpp_clock_id('(Ocean tracer: vert diffuse impl) ',grain=CLOCK_MODULE)
  id_clock_tracer_advect         = mpp_clock_id('(Ocean tracer: advection)    '     ,grain=CLOCK_MODULE)
  id_clock_frazil                = mpp_clock_id('(Ocean tracer: frazil)       '     ,grain=CLOCK_MODULE)
  id_clock_convection            = mpp_clock_id('(Ocean tracer: convection)   '     ,grain=CLOCK_MODULE)
  id_clock_blob                  = mpp_clock_id('(Ocean tracer: Lagrangian blobs) ' ,grain=CLOCK_MODULE)
  id_clock_blob_implicit         = mpp_clock_id('(Ocean tracer: implicit blobs)'    ,grain=CLOCK_ROUTINE)
  id_clock_adjust_L_thickness    = mpp_clock_id('(Ocean tracer: adjust the L thickness)',grain=CLOCK_ROUTINE)
  id_clock_adjust_E_thickness    = mpp_clock_id('(Ocean tracer: adjust the E thickness)',grain=CLOCK_ROUTINE)

  ! perform requested debugging for ocean tracer package manager
  if (ocean_tpm_debug) then  
    write (stdoutunit,*) ' '
    write (stdoutunit,*) 'Dumping field tree at start of ocean_prog_tracer_init'
    call write_timestamp(Time%model_time)
    if (.not. fm_dump_list('/', recursive = .true.)) then  
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer tree')
    endif
  endif  

  ! call routine to initialize the required tracer package
  if (otpm_set_tracer_package('required', caller=trim(mod_name)//'('//trim(sub_name)//')') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set required packages list')
  endif  
  
  ! perform requested debugging for ocean tracer package manager
  if (ocean_tpm_debug) then  !{
    write (stdoutunit,*) ' '
    write (stdoutunit,*) 'Dumping /ocean_mod field tree before ocean_tpm_init'
    call write_timestamp(Time%model_time)
    if (.not. fm_dump_list('/ocean_mod', recursive = .true.)) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer tree')
    endif  !}
  endif  !}

  ! call the initialization routine for the ocean tracer packages
  call ocean_tpm_init(Domain, Grid, Time, Time_steps, &
                      Ocean_options, debug) 

  ! perform requested debugging for ocean tracer package manager
  if (ocean_tpm_debug) then  
    write (stdoutunit,*) ' ' 
    write (stdoutunit,*) 'Dumping /ocean_mod field tree after ocean_tpm_init'
    call write_timestamp(Time%model_time)
    if (.not. fm_dump_list('/ocean_mod', recursive = .true.)) then  
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer tree')
    endif  
  endif  


  ! check for any errors in the number of fields in the tracer_packages list
  good_list => fm_util_get_string_array('/ocean_mod/GOOD/good_tracer_packages',   &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call fm_util_check_for_bad_fields('/ocean_mod/tracer_packages', good_list,    &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  else  
    call mpp_error(FATAL,trim(error_header) // ' Empty "good_tracer_packages" list')
  endif  

  ! check for any errors in the number of fields in the prog_tracers list
  good_list => fm_util_get_string_array('/ocean_mod/GOOD/good_prog_tracers',         &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call fm_util_check_for_bad_fields('/ocean_mod/prog_tracers', good_list,          &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  else  
    call mpp_error(FATAL,trim(error_header) // ' Empty "good_prog_tracers" list')
    endif  

  ! check for any errors in the number of fields in the namelists list
  good_list => fm_util_get_string_array('/ocean_mod/GOOD/good_namelists',            &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call fm_util_check_for_bad_fields('/ocean_mod/namelists', good_list,             &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  !else  
    !call mpp_error(FATAL,trim(error_header) // ' Empty "good_namelists" list')
  endif  

  ! get the number of tracers
  ! for now, we will not try to dynamically allocate the tracer arrays here
  num_prog_tracers = fm_get_length('/ocean_mod/prog_tracers')
  num_diag         = fm_get_length('/ocean_mod/diag_tracers')
  num_prog=num_prog_tracers
 
  ! allocate arrays based on the number of prognostic tracers
  allocate( T_prog            (num_prog_tracers) )
  allocate( id_eta_smooth     (num_prog_tracers) )
  allocate( id_eta_smooth_on_nrho(num_prog_tracers) )
  allocate( id_pbot_smooth    (num_prog_tracers) )
  allocate( id_prog           (num_prog_tracers) )
  allocate( id_prog_explicit  (num_prog_tracers) )
  allocate( id_prog_rhodzt    (num_prog_tracers) )
  allocate( id_prog_int_rhodz (num_prog_tracers) )
  allocate( id_prog_on_depth  (num_prog_tracers) )
  allocate( id_tendency_conc  (num_prog_tracers) )
  allocate( id_tendency       (num_prog_tracers) )
  allocate( id_tendency_on_nrho(num_prog_tracers) )
  allocate( id_tendency_expl  (num_prog_tracers) )
  allocate( id_surf_tracer    (num_prog_tracers) )
  allocate( id_surf_tracer_sq (num_prog_tracers) )
  allocate( id_bott_tracer    (num_prog_tracers) )

  if (use_blobs) then
     allocate( id_progT          (num_prog_tracers) )
     allocate( id_prog_rhodztT   (num_prog_tracers) )
     allocate( id_prog_int_rhodzT(num_prog_tracers) )
     allocate( id_surf_tracerT   (num_prog_tracers) )
     allocate( id_bott_tracerT   (num_prog_tracers) )
     allocate( id_tendency_concL (num_prog_tracers) )
     allocate( id_tendency_concT (num_prog_tracers) )
     allocate( id_tendencyL      (num_prog_tracers) )
     allocate( id_tendencyT      (num_prog_tracers) )
     allocate( id_tendency_explL (num_prog_tracers) )
     allocate( id_tendency_explT (num_prog_tracers) )
  endif

  allocate( Tracer_restart      (num_prog+num_diag) )
  allocate( id_tracer_restart   (num_prog+num_diag) )
  allocate( tracer_restart_file (num_prog+num_diag) )
  allocate( id_tracer_type      (num_prog+num_diag) )

  id_eta_smooth(:)     = -1
  id_eta_smooth_on_nrho(:)= -1
  id_pbot_smooth(:)    = -1
  id_prog(:)           = -1
  id_prog_explicit(:)  = -1
  id_prog_rhodzt(:)    = -1
  id_prog_int_rhodz(:) = -1
  id_prog_on_depth(:)  = -1
  id_tendency_conc(:)  = -1
  id_tendency(:)       = -1
  id_tendency_on_nrho(:)= -1
  id_tendency_expl(:)  = -1
  id_surf_tracer(:)    = -1
  id_surf_tracer_sq(:) = -1
  id_bott_tracer(:)    = -1

  if (use_blobs) then
     id_progT(:)          = -1
     id_prog_rhodztT(:)   = -1
     id_prog_int_rhodzT(:)= -1
     id_surf_tracerT(:)   = -1
     id_bott_tracerT(:)   = -1
     id_tendency_concL(:) = -1
     id_tendency_concT(:) = -1
     id_tendencyL(:)      = -1
     id_tendencyT(:)      = -1
     id_tendency_explL(:) = -1
     id_tendency_explT(:) = -1
  endif

  ! set logical to determine when to mpp_update tracers
  do n=1,num_prog_tracers-1
    T_prog(n)%complete=.false.
  enddo
  T_prog(num_prog_tracers)%complete=.true.

  !  dump the lists for the tracer packages
  write (stdoutunit,*)
  write (stdoutunit,*) 'Dumping tracer_packages tracer tree'
  if (.not. fm_dump_list('/ocean_mod/tracer_packages', recursive = .true.)) then  
    call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer_packages tracer tree')
  endif  

  write (stdoutunit,*)
  write (stdoutunit,*) 'Dumping prog_tracers tracer tree'
  if (.not. fm_dump_list('/ocean_mod/prog_tracers', recursive = .true.)) then  
    call mpp_error(FATAL, trim(error_header) // ' Problem dumping prog_tracers tracer tree')
  endif  

  write (stdoutunit,*)
  write (stdoutunit,*) 'Dumping namelists tracer tree'
  if (.not. fm_dump_list('/ocean_mod/namelists', recursive = .true.)) then 
    call mpp_error(FATAL, trim(error_header) // ' Problem dumping namelists tracer tree')
  endif  

  !### finished with initializing t_prog arrays
  !### finished with call to tracer startup routine


  ! set local array indices

  Grd => Grid
  Dom => Domain

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif

  do n=1,num_prog_tracers
#ifndef MOM_STATIC_ARRAYS
    allocate( T_prog(n)%field(isd:ied,jsd:jed,nk,3)) 
    allocate( T_prog(n)%th_tendency(isd:ied,jsd:jed,nk))
    allocate( T_prog(n)%tendency(isd:ied,jsd:jed,nk))
    allocate( T_prog(n)%source(isd:ied,jsd:jed,nk))
    allocate( T_prog(n)%eta_smooth(isd:ied,jsd:jed))
    allocate( T_prog(n)%pbot_smooth(isd:ied,jsd:jed))
    allocate( T_prog(n)%wrk1(isd:ied,jsd:jed,nk))
    allocate( T_prog(n)%tmask_limit(isd:ied,jsd:jed,nk))
    allocate( T_prog(n)%K33_implicit(isd:ied,jsd:jed,nk))
    if (use_blobs) then
       allocate( T_prog(n)%fieldT(isd:ied,jsd:jed,nk))
       allocate( T_prog(n)%sum_blob(isd:ied,jsd:jed,nk,3))
       allocate( T_prog(n)%tend_blob(isd:ied,jsd:jed,nk))
    endif
#endif
    T_prog(n)%field(:,:,:,:)        = 0.0
    T_prog(n)%th_tendency(:,:,:)    = 0.0
    T_prog(n)%tendency(:,:,:)       = 0.0
    T_prog(n)%source(:,:,:)         = 0.0
    T_prog(n)%eta_smooth(:,:)       = 0.0
    T_prog(n)%pbot_smooth(:,:)      = 0.0
    T_prog(n)%wrk1(:,:,:)           = 0.0
    T_prog(n)%tmask_limit(:,:,:)    = 0.0
    T_prog(n)%K33_implicit(:,:,:)   = 0.0
    T_prog(n)%neutral_physics_limit = .false.
    if (use_blobs) then
       !note that fieldT is initialised in ocean_blob_init
       T_prog(n)%fieldT(:,:,:)      = 0.0
       T_prog(n)%sum_blob(:,:,:,:)  = 0.0
       T_prog(n)%tend_blob(:,:,:)   = 0.0
    endif
  enddo

  if(inflow_nboundary) then 
    call inflow_nboundary_init
  endif 

  ! for global surface area normalization 
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

  ! fill the field table entries for the prognostic tracers 
  n = 0
  do while (fm_loop_over_list('/ocean_mod/prog_tracers', name, typ, ind))  !{

     if (typ .ne. 'list') then  !{

         call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')

     else  !}{

         n = n + 1  ! increment the array index

         if (n .ne. ind) then  !{
             write (stdoutunit,*) trim(warn_header), ' Tracer index, ', ind,   &
                  ' does not match array index, ', n, ' for ', trim(name)
         endif  !}

         ! save the name
         T_prog(n)%name = name

         if (.not. fm_change_list('/ocean_mod/prog_tracers/' // trim(name))) then  !{
             call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
         endif  !}

         caller_str = 'ocean_tracer_mod(ocean_prog_tracer_init)'

         ! save the units
         T_prog(n)%units = fm_util_get_string('units', caller = caller_str, scalar = .true.)

         ! save the type
         T_prog(n)%type = fm_util_get_string('type', caller = caller_str, scalar = .true.)

         ! save the longname
         T_prog(n)%longname = fm_util_get_string('longname', caller = caller_str, scalar = .true.)

         ! save the conversion
         T_prog(n)%conversion = fm_util_get_real('conversion', caller = caller_str, scalar = .true.)

         ! save the offset
         T_prog(n)%offset = fm_util_get_real('offset', caller = caller_str, scalar = .true.)

         ! get the min and max of the tracer
         T_prog(n)%min_tracer = fm_util_get_real('min_tracer', caller = caller_str, scalar = .true.)
         T_prog(n)%max_tracer = fm_util_get_real('max_tracer', caller = caller_str, scalar = .true.)

         ! get the min and max of the range for analysis
         T_prog(n)%min_range = fm_util_get_real('min_range', caller = caller_str, scalar = .true.)
         T_prog(n)%max_range = fm_util_get_real('max_range', caller = caller_str, scalar = .true.)

         ! get the flux unit
         T_prog(n)%flux_units = fm_util_get_string('flux_units', caller = caller_str, scalar = .true.)

         ! get the min and max of the flux range for analysis
         T_prog(n)%min_flux_range = fm_util_get_real('min_flux_range', caller = caller_str, scalar = .true.)
         T_prog(n)%max_flux_range = fm_util_get_real('max_flux_range', caller = caller_str, scalar = .true.)

         ! save the restart file
         T_prog(n)%restart_file = fm_util_get_string('restart_file', caller = caller_str, scalar = .true.)

         ! save flag for whether the tracer must have a value in the restart file
         T_prog(n)%const_init_tracer = fm_util_get_logical('const_init_tracer', caller = caller_str, scalar = .true.)

         ! save value to globally initialize this tracer (optional)
         T_prog(n)%const_init_value = fm_util_get_real('const_init_value', caller = caller_str, scalar = .true.)

         ! get the horizontal-advection-scheme
         string_fm = fm_util_get_string('horizontal-advection-scheme', caller = caller_str, scalar = .true.)

         select case (trim(string_fm)) 
         case ('upwind')
             T_prog(n)%horz_advect_scheme = ADVECT_UPWIND
         case ('2nd_order')
             T_prog(n)%horz_advect_scheme = ADVECT_2ND_ORDER
         case ('4th_order')
             T_prog(n)%horz_advect_scheme = ADVECT_4TH_ORDER
         case ('quicker')
             T_prog(n)%horz_advect_scheme = ADVECT_QUICKER
         case ('quickMOM3')
             T_prog(n)%horz_advect_scheme = ADVECT_QUICKMOM3
         case ('6th_order')
             T_prog(n)%horz_advect_scheme = ADVECT_6TH_ORDER              
         case ('mdfl_sup_b')
             T_prog(n)%horz_advect_scheme = ADVECT_MDFL_SUP_B
         case ('mdfl_sweby')
             T_prog(n)%horz_advect_scheme = ADVECT_MDFL_SWEBY
         case ('mdfl_sweby_test')
             T_prog(n)%horz_advect_scheme = ADVECT_MDFL_SWEBY_TEST
         case ('dst_linear')
             T_prog(n)%horz_advect_scheme = ADVECT_DST_LINEAR
         case ('dst_linear_test')
             T_prog(n)%horz_advect_scheme = ADVECT_DST_LINEAR_TEST
         case ('mdppm')
             T_prog(n)%horz_advect_scheme = ADVECT_MDPPM
             T_prog(n)%ppm_hlimiter = fm_util_get_integer('ppm_hlimiter', caller = caller_str, scalar = .true.)
         case ('mdppm_test')
             T_prog(n)%horz_advect_scheme = ADVECT_MDPPM_TEST
             T_prog(n)%ppm_hlimiter = fm_util_get_integer('ppm_hlimiter', caller = caller_str, scalar = .true.)
         case ('mdmdt_test')
             T_prog(n)%horz_advect_scheme = ADVECT_MDMDT_TEST
             T_prog(n)%mdt_scheme = fm_util_get_integer('mdt_scheme', caller = caller_str, scalar = .true.)
         case ('psom')
             T_prog(n)%horz_advect_scheme = ADVECT_PSOM
             ! set psom_limit when using psom
             T_prog(n)%psom_limit = fm_util_get_logical('psom_limit', caller = caller_str, scalar = .true.)
         case default
          call mpp_error(FATAL, trim(error_header) // ' Invalid horz-advect-scheme '  // trim(string_fm) // ' for ' // trim(name))
         end select

         ! get the vertical-advection-scheme
         string_fm = fm_util_get_string('vertical-advection-scheme', caller = caller_str, scalar = .true.)

         select case (trim(string_fm))    
         case ('upwind')
             T_prog(n)%vert_advect_scheme = ADVECT_UPWIND
         case ('2nd_order')
             T_prog(n)%vert_advect_scheme = ADVECT_2ND_ORDER
         case ('4th_order')
             T_prog(n)%vert_advect_scheme = ADVECT_4TH_ORDER
         case ('6th_order')
             T_prog(n)%vert_advect_scheme = ADVECT_6TH_ORDER              
         case ('quicker')
             T_prog(n)%vert_advect_scheme = ADVECT_QUICKER
         case ('quickMOM3')
             T_prog(n)%vert_advect_scheme = ADVECT_QUICKMOM3
         case ('mdfl_sup_b')
             T_prog(n)%vert_advect_scheme = ADVECT_MDFL_SUP_B
         case ('mdfl_sweby')
             T_prog(n)%vert_advect_scheme = ADVECT_MDFL_SWEBY
         case ('mdfl_sweby_test')
             T_prog(n)%vert_advect_scheme = ADVECT_MDFL_SWEBY_TEST
         case ('dst_linear')
             T_prog(n)%vert_advect_scheme = ADVECT_DST_LINEAR
         case ('dst_linear_test')
             T_prog(n)%vert_advect_scheme = ADVECT_DST_LINEAR_TEST
         case ('psom')
             T_prog(n)%vert_advect_scheme = ADVECT_PSOM
         case ('mdppm')
             T_prog(n)%vert_advect_scheme = ADVECT_MDPPM
             T_prog(n)%ppm_vlimiter = fm_util_get_integer('ppm_vlimiter', caller = caller_str, scalar = .true.)
         case ('mdppm_test')
             T_prog(n)%vert_advect_scheme = ADVECT_MDPPM_TEST
             T_prog(n)%ppm_vlimiter = fm_util_get_integer('ppm_vlimiter', caller = caller_str, scalar = .true.)
         case ('mdmdt_test')
             T_prog(n)%vert_advect_scheme = ADVECT_MDMDT_TEST
             T_prog(n)%mdt_scheme = fm_util_get_integer('mdt_scheme', caller = caller_str, scalar = .true.)
         case default
           call mpp_error(FATAL, trim(error_header) // ' Invalid vert-advect scheme ' // trim(string_fm) // ' for ' // trim(name))
         end select

         ! save the max_tracer_limit
         T_prog(n)%max_tracer_limit = fm_util_get_real('max_tracer_limit', caller = caller_str, scalar = .true.)

         ! save the min_tracer_limit
         T_prog(n)%min_tracer_limit = fm_util_get_real('min_tracer_limit', caller = caller_str, scalar = .true.)

         ! save flag for whether to transport tracer using only advection 
         T_prog(n)%use_only_advection = fm_util_get_logical('use_only_advection', caller = caller_str, scalar = .true.)
         if(T_prog(n)%use_only_advection) then
             call mpp_error(NOTE, &
             trim(note_header) // ' Will evolve '// trim(T_prog(n)%name) // ' via advection alone.') 
         endif

         ! check consistency in selection of three-dimensional advection schemes 
         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SUP_B .and. &
            T_prog(n)%vert_advect_scheme /= ADVECT_MDFL_SUP_B) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdfl_sup_b advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SWEBY .and. &
            T_prog(n)%vert_advect_scheme /= ADVECT_MDFL_SWEBY) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdfl_sweby advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_DST_LINEAR .and. &
            T_prog(n)%vert_advect_scheme /= ADVECT_DST_LINEAR) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' dst_linear advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_DST_LINEAR_TEST .and. &
            T_prog(n)%vert_advect_scheme /= ADVECT_DST_LINEAR_TEST) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' dst_linear_legacy advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_MDFL_SUP_B .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_MDFL_SUP_B) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdfl_sup_b advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_MDFL_SWEBY .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_MDFL_SWEBY) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdfl_sweby advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_MDFL_SWEBY_TEST .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_MDFL_SWEBY_TEST) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdfl_sweby_test advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_DST_LINEAR .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_DST_LINEAR) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' dst_linear advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_DST_LINEAR_TEST .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_DST_LINEAR_TEST) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' dst_linear_test advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_PSOM .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_PSOM) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' psom advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_MDPPM .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_MDPPM) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdppm advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_MDPPM_TEST .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_MDPPM_TEST) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdppm_test advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_MDMDT_TEST .and. &
            T_prog(n)%horz_advect_scheme /= ADVECT_MDMDT_TEST) then 
            call mpp_error(FATAL,&
            trim(error_header) // ' mdmdt_test advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_2ND_ORDER .or. &
            T_prog(n)%horz_advect_scheme==ADVECT_2ND_ORDER) then 
            if(tendency==TWO_LEVEL) then 
              call mpp_error(FATAL,&
              trim(error_header) // ' 2nd_order advection for ' // trim(name) // 'is unstable with TWO_LEVEL time') 
            endif 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_4TH_ORDER .or. &
            T_prog(n)%horz_advect_scheme==ADVECT_4TH_ORDER) then 
            if(tendency==TWO_LEVEL) then 
              call mpp_error(FATAL,&
              trim(error_header) // ' 4th_order advection for ' // trim(name) // 'is unstable with TWO_LEVEL time') 
            endif 
         endif

         if(T_prog(n)%vert_advect_scheme==ADVECT_6TH_ORDER .or. &
            T_prog(n)%horz_advect_scheme==ADVECT_6TH_ORDER) then 
            if(tendency==TWO_LEVEL) then 
              call mpp_error(FATAL,&
              trim(error_header) // ' 6th_order advection for ' // trim(name) // 'is unstable with TWO_LEVEL time') 
            endif 
         endif

     endif  !}
  enddo  !}


  allocate(id_tmask_limit(num_prog_tracers))

  do n=1,num_prog_tracers

      id_tmask_limit(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_tmask_limit', &
           Grd%tracer_axes(1:3), Time%model_time, &
           'use upwind (not quicker) and horz diff (not neutral)', &
           'none', missing_value=missing_value, range=(/-2.0,2.0/))

   if (use_blobs) then
      id_prog_rhodzt(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rhodztE', &
           Grd%tracer_axes(1:3),                                                                &
           Time%model_time, 'E system '//trim(T_prog(n)%longname)//' * rho_dzt',                &
           trim(T_prog(n)%units)//'*(kg/m^3)*m',                                                &
           missing_value=missing_value, range=(/-1e6,1e12/))
      
      id_prog_int_rhodz(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_int_rhodzE', &
           Grd%tracer_axes(1:2),                                                                      &
           Time%model_time, 'E system vertical sum of '//trim(T_prog(n)%longname)//' * rho_dzt',      &
           trim(T_prog(n)%units)//'*(kg/m^3)*m',                                                      &
           missing_value=missing_value, range=(/-1e20,1e20/))
      
   else !not use_blobs
      id_prog_rhodzt(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rhodzt', &
           Grd%tracer_axes(1:3),                                                               &
           Time%model_time, trim(T_prog(n)%longname)//' * rho_dzt',                            &
           trim(T_prog(n)%units)//'*(kg/m^3)*m',                                               &
           missing_value=missing_value, range=(/-1e6,1e12/))
      
      id_prog_int_rhodz(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_int_rhodz', &
           Grd%tracer_axes(1:2),                                                                     &
           Time%model_time, 'vertical sum of '//trim(T_prog(n)%longname)//' * rho_dzt',              &
           trim(T_prog(n)%units)//'*(kg/m^3)*m',                                                     &
           missing_value=missing_value, range=(/-1e20,1e20/))
   endif

   range_array(1) = T_prog(n)%min_range
   range_array(2) = T_prog(n)%max_range

   if (use_blobs) then
      prog_name     = trim(T_prog(n)%name)//'E'
      prog_longname = 'E system '//trim(T_prog(n)%longname)
   else
      prog_name = T_prog(n)%name
      prog_longname = T_prog(n)%longname
   endif

   if(T_prog(n)%longname=='Potential temperature') then 
       id_prog(n) = register_diag_field ('ocean_model',      &
            trim(prog_name),                                 &
            Grd%tracer_axes(1:3),                            &
            Time%model_time, prog_longname,                  &
            trim(temp_units),                                &
            missing_value=missing_value, range=range_array,  &
            standard_name='sea_water_potential_temperature')

       id_surf_tracer(n) = register_diag_field ('ocean_model', &
            'surface_'//trim(prog_name),                       &
            Grd%tracer_axes(1:2),                              &
            Time%model_time, trim(prog_longname),              &
            trim(temp_units),                                  &
            missing_value=missing_value, range=range_array,    &
            standard_name='sea_surface_temperature')

       id_surf_tracer_sq(n) = register_diag_field ('ocean_model', &
            'squared_surface_'//trim(prog_name),                  &
            Grd%tracer_axes(1:2),                                 &
            Time%model_time, 'squared '//trim(T_prog(n)%name),    &
            'squared '//trim(temp_units),                         &
            missing_value=missing_value, range=(/0.0,1e10/),      &
            standard_name='square_of_sea_surface_temperature')

   elseif(T_prog(n)%name=='salt') then 
       id_prog(n) = register_diag_field ('ocean_model',     &
            trim(prog_name),                                &
            Grd%tracer_axes(1:3),                           &
            Time%model_time, trim(prog_longname),           &
            trim(T_prog(n)%units),                          &
            missing_value=missing_value, range=range_array, &
            standard_name='sea_water_salinity')

       id_surf_tracer(n) = register_diag_field ('ocean_model', &
            'surface_'//trim(prog_name),                       &
            Grd%tracer_axes(1:2),                              &
            Time%model_time, trim(prog_longname),              &
            trim(T_prog(n)%units),                             &
            missing_value=missing_value, range=range_array,    &
            standard_name='sea_surface_salinity')

       id_surf_tracer_sq(n) = register_diag_field ('ocean_model', &
            'squared_surface_'//trim(prog_name),                  &
            Grd%tracer_axes(1:2),                                 &
            Time%model_time, 'squared '//trim(prog_longname),     &
            'squared '//trim(T_prog(n)%units),                    &
            missing_value=missing_value, range=(/0.0,1e10/),      &
            standard_name='square_of_sea_surface_salinity')

   elseif(T_prog(n)%name=='age_global') then 
       id_prog(n) = register_diag_field ('ocean_model',          &
            trim(prog_name),                                     &
            Grd%tracer_axes(1:3),                                &
            Time%model_time, trim(prog_longname),                &
            trim(T_prog(n)%units),                               &
            missing_value=missing_value, range=range_array,      &
            standard_name='sea_water_age_since_surface_contact')

       id_surf_tracer(n) = register_diag_field ('ocean_model', &
            'surface_'//trim(prog_name),                       &
            Grd%tracer_axes(1:2),                              &
            Time%model_time, trim(prog_longname),              &
            trim(T_prog(n)%units),                             &
            missing_value=missing_value, range=range_array)

       id_surf_tracer_sq(n) = register_diag_field ('ocean_model', &
            'squared_surface_'//trim(prog_name),                  &
            Grd%tracer_axes(1:2),                                 &
            Time%model_time, 'squared '//trim(prog_longname),     &
            'squared '//trim(T_prog(n)%units),                    &
            missing_value=missing_value, range=(/0.0,1e20/))

   elseif(T_prog(n)%name=='cfc_11') then 
       id_prog(n) = register_diag_field ('ocean_model',                &
            trim(prog_name),                                           &
            Grd%tracer_axes(1:3),                                      &
            Time%model_time, trim(prog_longname),                      &
            trim(T_prog(n)%units),                                     &
            missing_value=missing_value, range=range_array,            &
            standard_name='moles_per_unit_mass_of_cfc11_in_sea_water')

       id_surf_tracer(n) = register_diag_field ('ocean_model', &
            'surface_'//trim(prog_name),                       &
            Grd%tracer_axes(1:2),                              &
            Time%model_time, trim(prog_longname),              &
            trim(T_prog(n)%units),                             &
            missing_value=missing_value, range=range_array)

       id_surf_tracer_sq(n) = register_diag_field ('ocean_model', &
            'squared_surface_'//trim(prog_name),                  &
            Grd%tracer_axes(1:2),                                 &
            Time%model_time, 'squared '//trim(prog_longname),     &
            'squared '//trim(T_prog(n)%units),                    &
            missing_value=missing_value, range=(/0.0,1e20/))

   else

       id_prog(n) = register_diag_field ('ocean_model', trim(prog_name), &
            Grd%tracer_axes(1:3),                                        &
            Time%model_time, trim(prog_longname), trim(T_prog(n)%units), &
            missing_value=missing_value, range=range_array)

       id_surf_tracer(n) = register_diag_field ('ocean_model', &
            'surface_'//trim(prog_name),                       &
            Grd%tracer_axes(1:2),                              &
            Time%model_time, trim(prog_longname), &
            trim(T_prog(n)%units),                             &
            missing_value=missing_value, range=range_array)

       id_surf_tracer_sq(n) = register_diag_field ('ocean_model', &
            'squared_surface_'//trim(prog_name),                  &
            Grd%tracer_axes(1:2),                                 &
            Time%model_time, 'squared '//trim(prog_longname),     &
            'squared '//trim(T_prog(n)%units),                    &
            missing_value=missing_value, range=(/0.0,1e20/))

   endif

   id_prog_explicit(n) = register_diag_field ('ocean_model',                        &
        trim(T_prog(n)%name)//'_explicit',                                          &
        Grd%tracer_axes(1:3), Time%model_time,                                      &
        trim(T_prog(n)%longname)//' at taup1 updated just with explicit tendencies',&
        trim(T_prog(n)%units), missing_value=missing_value, range=range_array)
   id_prog_on_depth(n) = register_diag_field ('ocean_model',              &
        trim(T_prog(n)%name)//'_on_depth', Grd%tracer_axes_depth(1:3),    &
        Time%model_time, trim(T_prog(n)%longname), trim(T_prog(n)%units), &
        missing_value=missing_value, range=range_array)

   id_bott_tracer(n) = register_diag_field ('ocean_model', &
        'bottom_'//trim(prog_name),                        &
        Grd%tracer_axes(1:2),                              &
        Time%model_time, trim(prog_longname),              &
        trim(T_prog(n)%units),                             &
        missing_value=missing_value, range=range_array)
 
    ! register diagnostics for some tracer tendency terms

    if (T_prog(n)%max_flux_range .gt. T_prog(n)%min_flux_range) then  
      range_array(1) = T_prog(n)%min_flux_range
      range_array(2) = T_prog(n)%max_flux_range

      id_tendency_conc(n) = register_diag_field ('ocean_model',                               & 
           trim(prog_name)//'_tendency_conc', Grd%tracer_axes(1:3),                           &
           Time%model_time, 'time tendency of tracer concentration for '//trim(prog_longname),&
           trim(T_prog(n)%units)//' per second', missing_value=missing_value, range=range_array)
      id_tendency(n) = register_diag_field ('ocean_model',                    & 
           trim(prog_name)//'_tendency', Grd%tracer_axes(1:3),                &
           Time%model_time, 'time tendency for tracer '//trim(prog_longname), &
           trim(T_prog(n)%flux_units), missing_value=missing_value, range=range_array)
      id_eta_smooth(n) = register_diag_field ('ocean_model',                            &
           trim(prog_name)//'_eta_smooth', Grd%tracer_axes(1:2),                        &
           Time%model_time, 'surface smoother for ' // trim(prog_name),                 &
           trim(T_prog(n)%flux_units),                                                  &
           missing_value=missing_value, range=range_array)    
      id_pbot_smooth(n) = register_diag_field ('ocean_model',                           &
           trim(prog_name)//'_pbot_smooth', Grd%tracer_axes(1  :2),                     &
           Time%model_time, 'bottom smoother for ' // trim(prog_name),                  &
           trim(T_prog(n)%flux_units),                                                  &
           missing_value=missing_value, range=range_array)    

      ! The explicit tendency is the total system, whether blobs are used or not. 
      ! So, we make no distinction between T and E systems.
      id_tendency_expl(n) = register_diag_field ('ocean_model',                                  &
           trim(T_prog(n)%name)//'_tendency_expl', Grd%tracer_axes(1:3),                         &
           Time%model_time, 'explicit in time tendency for tracer ' // trim(T_prog(n)%longname), &
           trim(T_prog(n)%flux_units), missing_value=missing_value, range=range_array)    
    else  

      id_tendency_conc(n) = register_diag_field ('ocean_model',                          & 
           trim(prog_name)//'_tendency_conc', Grd%tracer_axes(1:3),                      &
           Time%model_time, 'tendency of tracer concentration for '//trim(prog_longname),&
           trim(T_prog(n)%units)//' per second', missing_value=missing_value, range=(/-1.e10,1e10/))
      id_tendency(n) = register_diag_field ('ocean_model',                              &
           trim(prog_name)//'_tendency', Grd%tracer_axes(1:3),                          &
           Time%model_time, 'time tendency for tracer '//trim(prog_longname),           &
           trim(T_prog(n)%flux_units), missing_value=missing_value)
      id_eta_smooth(n) = register_diag_field ('ocean_model',                                 &
           trim(T_prog(n)%name)//'_eta_smooth', Grd%tracer_axes(1:2),                        &
           Time%model_time, 'surface smoother for ' // trim(T_prog(n)%name),                 &
           trim(T_prog(n)%flux_units),                                                       &
           missing_value=missing_value)    
      id_pbot_smooth(n) = register_diag_field ('ocean_model',                                &
           trim(T_prog(n)%name)//'_pbot_smooth', Grd%tracer_axes(1:2),                       &
           Time%model_time, 'bottom smoother for ' // trim(T_prog(n)%name),                  &
           trim(T_prog(n)%flux_units),                                                       &
           missing_value=missing_value)    

      ! The explicit tendency is the total system, whether blobs are used or not. 
      ! So, we make no distinction between T and E systems.
      id_tendency_expl(n) = register_diag_field ('ocean_model',                        &
           trim(T_prog(n)%name)//'_tendency_expl', Grd%tracer_axes(1:3),                  &
           Time%model_time, 'explicit in time tendency for ' // trim(T_prog(n)%longname), &
           trim(T_prog(n)%flux_units), missing_value=missing_value)    
    endif  

    !Diagnostics for the Lagrangian blob framework
    if (use_blobs) then
        id_prog_rhodztT(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rhodzt',&
             Grd%tracer_axes(1:3),                                                               &
             Time%model_time, trim(T_prog(n)%longname)//' * rho_dzt',                            &
             trim(T_prog(n)%units)//'*(kg/m^3)*m',                                               &
             missing_value=missing_value, range=(/-1e6,1e12/))

        id_prog_int_rhodzT(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_int_rhodz',&
             Grd%tracer_axes(1:2),                                                                     &
             Time%model_time, trim(T_prog(n)%longname)//' * rho_dzt',                                  &
             trim(T_prog(n)%units)//'*(kg/m^3)*m',                                                     &
             missing_value=missing_value, range=(/-1e20,1e20/))

        if(T_prog(n)%longname=='Potential temperature') then 
            id_progT(n) = register_diag_field ('ocean_model',     &
                 trim(T_prog(n)%name),                            &
                 Grd%tracer_axes(1:3),                            &
                 Time%model_time, trim(T_prog(n)%longname),       &
                 trim(temp_units),                                &
                 missing_value=missing_value, range=range_array,  &
                 standard_name='sea_water_potential_temperature')

            id_surf_tracerT(n) = register_diag_field ('ocean_model', &
                 'surface_'//trim(T_prog(n)%name),                   &
                 Grd%tracer_axes(1:2),                               &
                 Time%model_time, trim(T_prog(n)%longname),          &
                 trim(temp_units),                                   &
                 missing_value=missing_value, range=range_array,     &
                 standard_name='sea_surface_temperature')

        elseif(T_prog(n)%name=='salt') then 
            id_progT(n) = register_diag_field ('ocean_model',    &
                 trim(T_prog(n)%name),                           &
                 Grd%tracer_axes(1:3),                           &
                 Time%model_time, trim(T_prog(n)%longname),      &
                 trim(T_prog(n)%units),                          &
                 missing_value=missing_value, range=range_array, &
                 standard_name='sea_water_salinity')

            id_surf_tracerT(n) = register_diag_field ('ocean_model', &
                 'surface_'//trim(T_prog(n)%name),                   &
                 Grd%tracer_axes(1:2),                               &
                 Time%model_time, trim(T_prog(n)%longname),          &
                 trim(T_prog(n)%units),                              &
                 missing_value=missing_value, range=range_array,     &
                 standard_name='sea_surface_salinity')
        elseif(T_prog(n)%name=='age_global') then 
            id_progT(n) = register_diag_field ('ocean_model',         &
                 trim(T_prog(n)%name),                                &
                 Grd%tracer_axes(1:3),                                &
                 Time%model_time, trim(T_prog(n)%longname),           &
                 trim(T_prog(n)%units),                               &
                 missing_value=missing_value, range=range_array,      &
                 standard_name='sea_water_age_since_surface_contact')

            id_surf_tracerT(n) = register_diag_field ('ocean_model', &
                 'surface_'//trim(T_prog(n)%name),                   &
                 Grd%tracer_axes(1:2),                               &
                 Time%model_time, trim(T_prog(n)%longname),          &
                 trim(T_prog(n)%units),                              &
                 missing_value=missing_value, range=range_array)

        elseif(T_prog(n)%name=='cfc_11') then 
            id_progT(n) = register_diag_field ('ocean_model',               &
                 trim(T_prog(n)%name),                                      &
                 Grd%tracer_axes(1:3),                                      &
                 Time%model_time, trim(T_prog(n)%longname),                 &
                 trim(T_prog(n)%units),                                     &
                 missing_value=missing_value, range=range_array,            &
                 standard_name='moles_per_unit_mass_of_cfc11_in_sea_water')

            id_surf_tracerT(n) = register_diag_field ('ocean_model', &
                 'surface_'//trim(T_prog(n)%name),                   &
                 Grd%tracer_axes(1:2),                               &
                 Time%model_time, trim(T_prog(n)%longname),          &
                 trim(T_prog(n)%units),                              &
                 missing_value=missing_value, range=range_array)

        else
            id_progT(n) = register_diag_field ('ocean_model',    &
                 trim(T_prog(n)%name),                           &
                 Grd%tracer_axes(1:3),                           &
                 Time%model_time, trim(T_prog(n)%longname),      &
                 trim(T_prog(n)%units),                          &
                 missing_value=missing_value, range=range_array)

            id_surf_tracerT(n) = register_diag_field ('ocean_model', &
                 'surface_'//trim(T_prog(n)%name),                   &
                 Grd%tracer_axes(1:2),                               &
                 Time%model_time, trim(T_prog(n)%longname),          &
                 trim(T_prog(n)%units),                              &
                 missing_value=missing_value, range=range_array)
        endif

        id_bott_tracerT(n) = register_diag_field ('ocean_model', &
             'bottom_'//trim(T_prog(n)%name),                    &
             Grd%tracer_axes(1:2),                               &
             Time%model_time, trim(T_prog(n)%longname),          &
             trim(T_prog(n)%units),                              &
             missing_value=missing_value, range=range_array)

        id_tendency_concT(n) = register_diag_field ('ocean_model',                              & 
             trim(T_prog(n)%name)//'_tendency_conc', Grd%tracer_axes(1:3),                      &
             Time%model_time, 'tendency of tracer concentration for '//trim(T_prog(n)%longname),&
             trim(T_prog(n)%units)//' per second', missing_value=missing_value, range=(/-1.e10,1e10/))
        id_tendency_concL(n) = register_diag_field ('ocean_model',                              & 
             trim(T_prog(n)%name)//'_tendency_concL', Grd%tracer_axes(1:3),                     &
             Time%model_time, 'L contribution to tendency of tracer concentration for '//trim(T_prog(n)%longname),&
             trim(T_prog(n)%units)//' per second', missing_value=missing_value, range=(/-1.e10,1e10/))

        id_tendencyT(n) = register_diag_field ('ocean_model',                          &
             trim(T_prog(n)%name)//'_tendency', Grd%tracer_axes(1:3),                  &
             Time%model_time, 'time tendency for tracer '//trim(T_prog(n)%longname),   &
             trim(T_prog(n)%flux_units), missing_value=missing_value)
        id_tendencyL(n) = register_diag_field ('ocean_model',                             &
             trim(T_prog(n)%name)//'L_tendency', Grd%tracer_axes(1:3),                    &
             Time%model_time, 'time tendency for tracer '//trim(T_prog(n)%longname)//'L', &
             trim(T_prog(n)%flux_units), missing_value=missing_value)

    endif  ! endif for use_blobs 


    if (T_prog(n)%name == 'temp')        index_temp        = n
    if (T_prog(n)%name == 'salt')        index_salt        = n
    if (T_prog(n)%name == 'temp_sq')     index_temp_sq     = n
    if (T_prog(n)%name == 'salt_sq')     index_salt_sq     = n
    if (T_prog(n)%name == 'redist_heat') index_redist_heat = n
    if (T_prog(n)%name == 'added_heat')  index_added_heat  = n
    if (T_prog(n)%longname == 'Conservative temperature') prog_temp_variable=CONSERVATIVE_TEMP
    if (T_prog(n)%longname == 'Potential temperature')    prog_temp_variable=POTENTIAL_TEMP

 enddo  ! n-loop 

  if(prog_temp_variable==CONSERVATIVE_TEMP) then
      write (stdoutunit,'(/1x,a)') &
      ' ==> Note from ocean_tracer_mod: prognostic temperature = conservative temperature.'
      write (stdoutunit,'(1x,a/)') &
      '                                 diagnostic temperature = potential temperature.'
  elseif(prog_temp_variable==POTENTIAL_TEMP) then
      write (stdoutunit,'(/1x,a)') &
      ' ==> Note from ocean_tracer_mod: prognostic temperature = potential temperature.'
      write (stdoutunit,'(1x,a/)') &
      '                                 diagnostic temperature = conservative temperature.'
  else
      call mpp_error(FATAL, &
      '==>Error in ocean_tracer_mod: model temperature variable remains unspecified')
  endif

  if (index_temp == -1 .or. index_salt == -1 ) then 
    call mpp_error(FATAL, trim(error_header) // &
    ' temp or salt not included in table entries.  MOM needs both.')
  endif 

#ifdef USE_OCEAN_BGC 
  !Point the %field4d_ptr and %tendency for "generic" tracers to the corresponding T_prog(n)%field and T_prog(n)%K33_implicit 
  !to utilize the already allocated memory in MOM and to avoid copying the arrays back and forth between generic tracers and MOM
  do n=1,num_prog_tracers
    if(T_prog(n)%type .eq. 'generic') then
       call ocean_generic_set_pointer(T_prog(n)%name, 'field', T_prog(n)%field, isd, jsd)
       !T_prog(n)%K33_implicit is used in vertdiff method below for calculating vertical diffusivity
       call ocean_generic_set_pointer(T_prog(n)%name, 'tendency', T_prog(n)%K33_implicit, isd, jsd)
    endif
  enddo
#endif

  ! read prognostic tracer initial conditions or restarts
  write (stdoutunit,*) ' '  
  write (stdoutunit,*) trim(note_header), &
   ' Reading prognostic tracer initial conditions or restarts'

  do n=1,num_prog_tracers  !{

    write (stdoutunit,*)
    write (stdoutunit,*) &
    'Initializing tracer number', n,' at time level tau. This tracer is called ',trim(T_prog(n)%name)

    if(.not. Time%init) then 
        if(tendency==TWO_LEVEL) then 
            write (stdoutunit,'(/a)')'Expecting only one time record from the tracer restart.' 
        elseif(tendency==THREE_LEVEL) then  
            write (stdoutunit,'(/a)')'Expecting two time records from the tracer restart.' 
        endif
    endif

    T_prog(n)%field(:,:,:,taum1) = 0.0

    ! initialization logic
    ! 
    ! 1. Always call passive_tracer_init once. 
    !    If Time%init=.false. then return w/o doing anything
    !    If not using passive tracers, then return w/o doing anything
    !
    ! 2. If filename exists with field in it, then fill the tracer with field in filename.  
    !
    ! 3. If filename exists but there is no field in it, then fill the tracer with const_init_value
    !    if const_init_tracer is true, otherwise end with a fatal error.
    !
    ! 4. If filename does not exist, then fill the tracer with const_init_value
    !    if const_init_tracer is true, otherwise end with a fatal error.
    !

    filename = T_prog(n)%restart_file    

    call passive_tracer_init(Time%init, T_prog(n), initialize_as_a_passive_tracer)
    do l = 1, num_tracer_restart
       if(trim(filename) == tracer_restart_file(l)) exit
    end do    
    if(l>num_tracer_restart) then
       num_tracer_restart = num_tracer_restart + 1
       tracer_restart_file(l) = trim(filename)
    end if
    id_tracer_type(n) = l
    if(tendency==THREE_LEVEL) then
       id_tracer_restart(n) = register_restart_field(Tracer_restart(l), filename, T_prog(n)%name, &
             T_prog(n)%field(:,:,:,tau), T_prog(n)%field(:,:,:,taup1), Domain%domain2d)
    else
       id_tracer_restart(n) = register_restart_field(Tracer_restart(l), filename, T_prog(n)%name, &
             T_prog(n)%field(:,:,:,tau), Domain%domain2d)
    end if
    filename = 'INPUT/'//trim(T_prog(n)%restart_file)

    if( Time%init .and. T_prog(n)%const_init_tracer ) then 
       write (stdoutunit,*) &
       'Initializing the tracer ',trim(T_prog(n)%name),' to the constant ',T_prog(n)%const_init_value 
       T_prog(n)%field(:,:,:,1) = T_prog(n)%const_init_value*Grd%tmask(:,:,:) 
       T_prog(n)%field(:,:,:,2) = T_prog(n)%const_init_value*Grd%tmask(:,:,:) 
       T_prog(n)%field(:,:,:,3) = T_prog(n)%const_init_value*Grd%tmask(:,:,:)
    elseif (field_exist(trim(filename), trim(T_prog(n)%name))) then 

       call field_size(filename,T_prog(n)%name, siz)
       if (tendency==TWO_LEVEL .and. siz(4) > 1) then
          write(stdoutunit,'(/a)') &
          '==>WARMING: Attempt to read restart for tracer from a 3-level time scheme (2 time records)'
          write(stdoutunit,'(a)') &
          '          when running MOM with 2-level timestepping (only need 1 time record in restart).'
          write(stdoutunit,'(a)')  &
          '          Reduce restart file to a single time record in order to avoid confusion.'
          call mpp_error(WARNING, &
          'Reading 3-time level ocean tracer (w/ 2 time records) while using 2-level (needs only 1 record)')
       endif

       write (stdoutunit,*) 'Reading restart for prog tracer ', trim(T_prog(n)%name), ' from file ', &
            trim(T_prog(n)%restart_file)
       call restore_state(Tracer_restart(l), id_tracer_restart(n))

      ! modify initial conditions 
       if (Time%init) then  

           if(remap_depth_to_s_init .and. vert_coordinate_type == QUASI_HORIZONTAL) then 
              call mpp_error(WARNING, trim(note_header) // &
              'No need to remap init-conditions from depth to s-coord with quasi-horizontal vert-coord.')
           endif 
           if(.not. remap_depth_to_s_init .and. vert_coordinate_type == TERRAIN_FOLLOWING) then 
              call mpp_error(WARNING, trim(note_header) // &
              'May need to remap init-conditions from depth to s-coord with terrain following vert-coord.')
           endif 
           if(remap_depth_to_s_init .and. vert_coordinate_type == TERRAIN_FOLLOWING) then 
              call mpp_error(NOTE, trim(note_header) // &
              'Remap init-tracer from depth to s-coord for terrain following vert-coord.')
           endif 

           ! remap the depth-based initial conditions onto s-coordinates 
           if(remap_depth_to_s_init) then 
               write (stdoutunit,*) &
               'After reading tracer init, remapping depth to s-coord for tracer ',trim(T_prog(n)%name)
               call remap_depth_to_s(Thickness, Time, T_prog(n)%field(:,:,:,tau))
           endif 

           ! modifications due to partial bottom cell
           if(vert_coordinate_type == QUASI_HORIZONTAL .and. interpolate_tprog_to_pbott) then 
               write (stdoutunit,*) &
               'After reading ic, linearly interpolate ',trim(T_prog(n)%name),' to partial cell bottom.'
               do j=jsc,jec
                  do i=isc,iec
                     kb = Grd%kmt(i,j)
                     if (kb .gt. 1) then
                         if(use_blobs) then
                            fact = Thickness%dzwtT(i,j,kb-1)/Grd%dzw(kb-1)
                         else
                            fact = Thickness%dzwt(i,j,kb-1)/Grd%dzw(kb-1)
                         endif
                         T_prog(n)%field(i,j,kb,tau) =                  &
                              T_prog(n)%field(i,j,kb-1,tau) -           &
                              fact*(T_prog(n)%field(i,j,kb-1,tau) -     &
                              T_prog(n)%field(i,j,kb,tau))
                     endif
 
                     ! zero out tracers in land points
                     do k = kb+1,nk
                        T_prog(n)%field(i,j,k,tau) = 0.0
                     enddo
 
                  enddo
               enddo
           endif

       endif

    else 

       if(.not. initialize_as_a_passive_tracer) then
          call mpp_error(FATAL,                                                          &
          trim(T_prog(n)%name) // ' is not initialized as an ideal passive tracer, ' //  &
          ' it is not initialized as constant, and it does not exist in the file '       &
          // trim(filename) //                                                           &
          'All tracers must have initialization specified.  There is no default.')
       endif 

    endif 

    write (stdoutunit,*) &
    'Completed initialization of tracer ',trim(T_prog(n)%name),' at time level tau'

    ! fill halos for tau tracer value 
    call mpp_update_domains(T_prog(n)%field(:,:,:,tau), Dom%domain2d)

    ! read in second time level when running with threelevel scheme 
    if(tendency==THREE_LEVEL .and. .not.Time%init) then 
      write (stdoutunit,*) 'Since Time%init=.false. and using threelevel tendency, read taup1 for ',trim(T_prog(n)%name)
      write (stdoutunit,*) 'Completed read of restart for ', trim(T_prog(n)%name),' at time taup1'
      call mpp_update_domains(T_prog(n)%field(:,:,:,taup1), Dom%domain2d)
    else
      ! initialize tracer at taup1 to tracer at tau 
      T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,tau) 
    endif 

  enddo  !} n-loop 

  write (stdoutunit,*) ' '  
  write (stdoutunit,*) trim(note_header), ' finished reading prognostic tracer restarts.'

  if(have_obc) call ocean_obc_tracer_init(Time, T_prog, num_prog_tracers, debug_this_module)

  prog_module_initialized = .true.

  if (debug_this_module) then
    do n=1,num_prog_tracers
      write(stdoutunit,'(a)') ' '
      write(stdoutunit,*) 'From ocean_tracer_mod: initial ',trim(T_prog(n)%name), ' chksum (tau) ==>'
      call tracer_prog_chksum(Time, T_prog(n), tau)
      write(stdoutunit,*) 'From ocean_tracer_mod: initial ',trim(T_prog(n)%name), ' chksum (taup1) ==>'
      call tracer_prog_chksum(Time, T_prog(n), taup1)
      call tracer_min_max(Time, Thickness, T_prog(n))
    enddo
  endif

end function ocean_prog_tracer_init  
! </FUNCTION> NAME="ocean_prog_tracer_init">


!#######################################################################
! <FUNCTION NAME="ocean_diag_tracer_init">
!
! <DESCRIPTION>
! Initialization code for diagnostic tracers, returning a pointer to the T_diag array
! </DESCRIPTION>
!
function ocean_diag_tracer_init (Time, Thickness, vert_coordinate_type, num_diag, use_blobs)   &
     result (T_diag)  !{
  
  type(ocean_time_type), intent(in), target :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  integer,                     intent(in)   :: vert_coordinate_type 
  integer, intent(out)                      :: num_diag
  logical, intent(in)                       :: use_blobs

!
!       Return type
!

  type(ocean_diag_tracer_type), dimension(:), pointer :: T_diag

  integer :: kb, i, j, k, n, l
  character(len=32)     :: name
  character(len=128)    :: filename
  real :: fact
  character(len=48), parameter  :: sub_name = 'ocean_diag_tracer_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  ! variables for tracer package
  integer                               :: ind, id_restart
  real, dimension(2)                    :: range_array
  character(len=64)                     :: caller_str
  character(len=fm_type_name_len)       :: typ
  character(len=fm_string_len), pointer, dimension(:) :: good_list
  integer :: stdoutunit
  stdoutunit=stdout()


  
  if (diag_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' Diagnostic tracers already initialized')
  endif

  nullify(T_diag)

  
  ! calls to internally generated (i.e., no field table entries) diagnostic tracers
  call ocean_irradiance_init

  
!
!       Check for any errors in the number of fields in the diag_tracers list
!

  good_list => fm_util_get_string_array('/ocean_mod/GOOD/good_diag_tracers', &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  !{
    call fm_util_check_for_bad_fields('/ocean_mod/diag_tracers', good_list,  &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  endif  !}


! Get the number of tracers
! For now, we will not try to dynamically allocate the tracer arrays here

  num_diag_tracers = fm_get_length('/ocean_mod/diag_tracers')

  num_diag = num_diag_tracers

!
!       Only do the following if there are diag tracers
!

  if (num_diag_tracers .gt. 0) then  !{

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), ' ', num_diag_tracers, ' diagnostic tracers requested.'

!
!       Allocate the T_diag array
!

    allocate( T_diag(num_diag_tracers) )
    allocate( id_diag(num_diag_tracers) )
    allocate( id_diag_total(num_diag_tracers) )
    allocate( id_diag_surf_tracer   (num_diag_tracers) )
    allocate( id_diag_surf_tracer_sq(num_diag_tracers) )
    id_diag(:)                = -1
    id_diag_total(:)          = -1
    id_diag_surf_tracer(:)    = -1
    id_diag_surf_tracer_sq(:) = -1
!
!       Dump the diagnostic tracer tree
!

    write (stdoutunit,*)
    write (stdoutunit,*) 'Dumping ocean diag field tree after reading diag tracer tree'
    if (.not. fm_dump_list('/ocean_mod/diag_tracers', recursive = .true.)) then  
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping diag_tracers tracer tree')
    endif  

  !### finished with initializing t_diag arrays
  !### finished with call to tracer startup routine


! allocate tracer type arrays and fill in tracer table information  
    if (num_diag_tracers .gt. 0 ) then
      do n=1,num_diag_tracers
#ifndef MOM_STATIC_ARRAYS
        allocate(T_diag(n)%field(isd:ied,jsd:jed,nk))
#endif
        T_diag(n)%field(:,:,:) = 0.0
      enddo
    endif

!
!       process the T_diag array
!

    n = 0
    do while (fm_loop_over_list('/ocean_mod/diag_tracers', name, typ, ind))  !{
     
      if (typ .ne. 'list') then  !{
         
        call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')
         
      else  !}{

        n = n + 1

        if (n .ne. ind) then  
          write (stdoutunit,*) trim(warn_header), ' Tracer index, ', ind,   &
               ' does not match array index, ', n, ' for ', trim(name)
        endif  

        ! save the name (required)
        T_diag(n)%name = name

        if (.not. fm_change_list('/ocean_mod/diag_tracers/' // trim(name))) then  
          call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
        endif  

        caller_str = 'ocean_tracer_mod(ocean_diag_tracer_init)'

        ! save the units (optional)
        T_diag(n)%units = fm_util_get_string('units', caller = caller_str, scalar = .true.)

        ! save the type (optional)
        T_diag(n)%type = fm_util_get_string('type', caller = caller_str, scalar = .true.)

        ! save the longname (optional)
        T_diag(n)%longname = fm_util_get_string('longname', caller = caller_str, scalar = .true.)

        ! save the conversion (optional)
        T_diag(n)%conversion = fm_util_get_real('conversion', caller = caller_str, scalar = .true.)

        ! save the offset (optional)
        T_diag(n)%offset = fm_util_get_real('offset', caller = caller_str, scalar = .true.)

        ! get the min and max of the tracer for analysis (optional)
        T_diag(n)%min_tracer = fm_util_get_real('min_tracer', caller = caller_str, scalar = .true.)
        T_diag(n)%max_tracer = fm_util_get_real('max_tracer', caller = caller_str, scalar = .true.)

        ! get the min and max of the range for analysis (optional)
        T_diag(n)%min_range = fm_util_get_real('min_range', caller = caller_str, scalar = .true.)
        T_diag(n)%max_range = fm_util_get_real('max_range', caller = caller_str, scalar = .true.)

        ! save the input restart file (optional)
        T_diag(n)%restart_file = fm_util_get_string('restart_file', caller = caller_str, scalar = .true.)

        ! save flag for whether to initialize this tracer (optional)
        T_diag(n)%const_init_tracer = fm_util_get_logical('const_init_tracer', caller = caller_str, scalar = .true.)

        ! save value to globally initialize this tracer (optional)
        T_diag(n)%const_init_value = fm_util_get_real('const_init_value', caller = caller_str, scalar = .true.)

      endif  !}
    enddo  !}


    ! register diagnostics for the tracer
    do n=1,num_diag_tracers 

       range_array(1) = T_diag(n)%min_range
       range_array(2) = T_diag(n)%max_range

       if(T_diag(n)%longname=='Potential temperature') then 
           id_diag(n) = register_diag_field ('ocean_model', trim(T_diag(n)%name), &
                Grd%tracer_axes(1:3),                                             &
                Time%model_time, trim(T_diag(n)%longname), trim(temp_units),      &
                missing_value=missing_value, range=range_array,                   &
                standard_name='sea_water_potential_temperature')
           id_diag_surf_tracer(n) = register_diag_field ('ocean_model',           &
                'surface_'//trim(T_diag(n)%name),                                 &
                Grd%tracer_axes(1:2),                                             &
                Time%model_time, trim(T_diag(n)%longname), trim(temp_units),      &
                missing_value=missing_value, range=range_array,                   &
                standard_name='sea_surface_temperature')
           id_diag_surf_tracer_sq(n) = register_diag_field ('ocean_model',        &
                'squared_surface_'//trim(T_diag(n)%name),                         &
                Grd%tracer_axes(1:2),                                             &
                Time%model_time, 'squared '//trim(T_diag(n)%longname),            &
                'squared '//trim(temp_units),                                     &
                missing_value=missing_value, range=range_array,                   &
                standard_name='square_of_sea_surface_temperature')

       elseif(T_diag(n)%longname=='Conservative temperature') then 
           id_diag(n) = register_diag_field ('ocean_model', trim(T_diag(n)%name), &
                Grd%tracer_axes(1:3),                                             &
                Time%model_time, trim(T_diag(n)%longname), trim(temp_units),      &
                missing_value=missing_value, range=range_array)
           id_diag_surf_tracer(n) = register_diag_field ('ocean_model',           &
                'surface_'//trim(T_diag(n)%name),                                 &
                Grd%tracer_axes(1:2),                                             &
                Time%model_time, trim(T_diag(n)%longname), trim(temp_units),      &
                missing_value=missing_value, range=range_array)
           id_diag_surf_tracer_sq(n) = register_diag_field ('ocean_model',        &
                'squared_surface_'//trim(T_diag(n)%name),                         &
                Grd%tracer_axes(1:2),                                             &
                Time%model_time, 'squared '//trim(T_diag(n)%longname),            &
                'squared '//trim(temp_units),                                     &
                missing_value=missing_value, range=range_array)

       else 

           id_diag(n) = register_diag_field ('ocean_model', trim(T_diag(n)%name), &
                Grd%tracer_axes(1:3),                                             &
                Time%model_time, trim(T_diag(n)%longname), trim(T_diag(n)%units), &
                missing_value=missing_value)
           id_diag_surf_tracer(n) = register_diag_field ('ocean_model',           &
                'surface_'//trim(T_diag(n)%name),                                 &
                Grd%tracer_axes(1:2),                                             &
                Time%model_time, trim(T_diag(n)%longname), trim(T_diag(n)%units), &
                missing_value=missing_value, range=range_array)
           id_diag_surf_tracer_sq(n) = register_diag_field ('ocean_model',        &
                'squared_surface_'//trim(T_diag(n)%name),                         &
                Grd%tracer_axes(1:2),                                             &
                Time%model_time, 'squared '//trim(T_diag(n)%longname),            &
                'squared '//trim(T_diag(n)%units),                                &
                missing_value=missing_value, range=range_array)

       endif

       id_diag_total(n) = register_diag_field ('ocean_model','total_ocean'//trim(T_diag(n)%name), &
            Time%model_time, 'Total ocean mass for tracer '//trim(T_diag(n)%name), 'kg/1e18',     &
            missing_value=missing_value, range=(/0.0,1e20/))

    enddo


#ifdef USE_OCEAN_BGC 
   !Get the %filed for "generic" tracers as it might have already been set.
   if (num_diag_tracers .gt. 0 ) then 
      do n=1,num_diag_tracers
         if(T_diag(n)%type .eq. 'generic') then
            !call ocean_generic_get_field(T_diag(n)%name,T_diag(n)%field)
            call ocean_generic_set_pointer(T_diag(n)%name, 'field', T_diag(n)%field, isd, jsd)
            
            !nnz: find a way to use their already allocated field pointer directly.
            !#ifndef MOM_STATIC_ARRAYS
            !!        deallocate(T_diag(n)%field)
            !!#endif
            !!             call ocean_generic_get_field_pointer(T_diag(n)%name,T_diag(n)%field)
            
         endif
      enddo
   endif
#endif
   
    ! read diagnostic tracer restarts 
    write (stdoutunit,*) trim(note_header), ' Reading diagnostic tracer initial conditions and/or restarts'
    do n=1,num_diag_tracers  !{

       write (stdoutunit,*)
       write (stdoutunit,*) &
       'Initializing tracer number', n,' at time level tau. This tracer is called ',trim(T_diag(n)%name)

       if (T_diag(n)%restart_file .eq. ' ') then  !{
          write (stdoutunit,*) 'Skipping tracer ', trim(T_diag(n)%name)
       else  !}{

         filename = T_diag(n)%restart_file
 
         do l = 1, num_tracer_restart
           if(trim(filename) == tracer_restart_file(l)) exit
         end do
         if (l>num_tracer_restart) then
           num_tracer_restart = num_tracer_restart + 1
           tracer_restart_file(l) = trim(filename)
         end if
         id_restart = register_restart_field(Tracer_restart(l), filename, T_diag(n)%name,       &
              T_diag(n)%field(:,:,:), Dom%domain2d)

         filename = 'INPUT/'//trim(T_diag(n)%restart_file)

         if (Time%init .and. T_diag(n)%const_init_tracer ) then 
           write (stdoutunit,*) &
                'Initializing diagnostic tracer ', trim(T_diag(n)%name),' to constant ', T_diag(n)%const_init_value 
           T_diag(n)%field(:,:,:) = T_diag(n)%const_init_value 
         elseif (field_exist(filename, T_diag(n)%name)) then 
           write (stdoutunit,*) 'Reading restart for diag tracer ', trim(T_diag(n)%name), ' from file ', &
                trim(T_diag(n)%restart_file)
           call restore_state(Tracer_restart(l), id_restart)

           write (stdoutunit,*) 'Completed read of ic/restart for diagnostic tracer ', trim(T_diag(n)%name)

           if (Time%init) then  !{

             if(remap_depth_to_s_init .and. vert_coordinate_type == QUASI_HORIZONTAL) then 
                call mpp_error(WARNING, trim(note_header) // &
                'No need to remap init-conditions from depth to s-coord with quasi-horizontal vert-coord.')
             endif 
             if(.not. remap_depth_to_s_init .and. vert_coordinate_type == TERRAIN_FOLLOWING) then 
                call mpp_error(WARNING, trim(note_header) // &
                'May need to remap init-conditions from depth to s-coord with terrain following vert-coord.')
             endif 
             if(remap_depth_to_s_init .and. vert_coordinate_type == TERRAIN_FOLLOWING) then 
                call mpp_error(NOTE, trim(note_header) // &
                'Remap init-tracer from depth to s-coord for terrain following vert-coord.')
             endif 

             ! remap the depth-based initial conditions onto s-coordinates 
             if(remap_depth_to_s_init) then 
                 write (stdoutunit,*) &
                 'After reading tracer init, remapping depth to s-coord for tracer ',trim(T_diag(n)%name)
                 call remap_depth_to_s(Thickness, Time, T_diag(n)%field(:,:,:))
             endif 

             ! modifications due to partial bottom cell
             if (vert_coordinate_type == QUASI_HORIZONTAL .and. interpolate_tdiag_to_pbott) then 

               ! linearly interpolate tracers to partial bottom cells
               write (stdoutunit,*) 'Linearly interpolate tracer ',trim(T_diag(n)%name),' to partial cell bottom.'
               do j = jsc, jec  !{
                  do i = isc, iec  !{
                     kb = Grd%kmt(i,j)
                     if (kb .gt. 1) then  !{
                         if(use_blobs) then
                            fact = Thickness%dzwtT(i,j,kb-1)/Grd%dzw(kb-1) 
                         else
                            fact = Thickness%dzwt(i,j,kb-1)/Grd%dzw(kb-1) 
                         endif
                         T_diag(n)%field(i,j,kb) =                      &
                               (1.-fact)*T_diag(n)%field(i,j,kb-1) + fact*T_diag(n)%field(i,j,kb)
                     endif  !}

                     ! zero out tracers in land points
                     do k = kb+1,nk  !{
                        T_diag(n)%field(i,j,k) = 0.0
                     enddo  !} k

                  enddo  !} i
               enddo  !} j

             endif

           endif  !}

        else
          call mpp_error(FATAL, trim(T_diag(n)%name) &
               // ' is not initialized as constant, and it does not exist in the file '       &
               // trim(filename) //                                                           &
               '. All tracers must have initialization specified.  There is no default.')
        endif

        call mpp_update_domains(T_diag(n)%field(:,:,:), Dom%domain2d)

      endif  !}

    enddo  !} n
    write (stdoutunit,*) trim(note_header), ' Finished reading diagnostic tracer restarts.'

    diag_module_initialized = .true.

    if (debug_this_module) then
      write(stdoutunit,'(a)') ' '
      do n=1,num_diag_tracers
        if (T_diag(n)%restart_file .ne. ' ') then  !{
          write(stdoutunit,*) 'From ocean_tracer_mod: initial ',trim(T_diag(n)%name), ' chksum (diag) ==>'
          call tracer_diag_chksum(Time, T_diag(n))
        endif  !}
      enddo
    endif

    do n=1,num_diag_tracers
       if (T_diag(n)%longname == 'Conservative temperature') index_diag_temp = n
       if (T_diag(n)%longname == 'Potential temperature')    index_diag_temp = n
    enddo
    if (index_diag_temp == -1) then 
      call mpp_error(FATAL, trim(error_header) // ' diagnostic temperature not included.  MOM needs it.')
    endif 


  else  !}{

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'No diagnostic tracers requested.'

  endif  !}



end function ocean_diag_tracer_init  !}
! </FUNCTION> NAME="ocean_diag_tracer_init">


!#######################################################################
! <SUBROUTINE NAME="update_ocean_tracer">
!
! <DESCRIPTION>
! Update value of tracer concentration to time taup1.
!
! Note that T_prog(n)%source is added at the very end 
! of the time step, after the rho_dzt factor has been 
! divided. 
! </DESCRIPTION>
!
subroutine update_ocean_tracer (Time, Dens, Adv_vel, Thickness, pme, diff_cbt, &
                                T_prog, T_diag, L_system,   &
                                Velocity, Ext_mode, EL_diag, use_blobs)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_density_type),       intent(inout) :: Dens
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_thickness_type),     intent(inout) :: Thickness
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:,:,:), intent(in)    :: diff_cbt

  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(blob_diag_type),           intent(inout) :: EL_diag(0:)
  logical,                        intent(in)    :: use_blobs
 
  type(time_type) :: next_time
  type(time_type) :: time_step

  real      :: age_tracer_max
  integer   :: i, j, k, kbot, n
  integer   :: taum1, tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  if (size(T_prog(:)) < num_prog_tracers) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_mod (update_ocean_tracer): T_prog array size too small')
  endif 

  if (zero_tendency) then

      ! hold tracer fields constant in time
      do n=1,num_prog_tracers
         T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,tau)
      enddo

  else 

     time_step = set_time(int(dtts),0)
     next_time = Time%model_time + time_step

     ! add contributions to T_prog%th_tendency terms arising
     ! from time explicit advection and diffusion
     do n=1,num_prog_tracers

        if(n==index_temp_sq .or. n==index_salt_sq) cycle

        call mpp_clock_begin(id_clock_tracer_advect)
          call horz_advect_tracer(Time, Adv_vel, Thickness, Dens, &
            T_prog(1:num_prog_tracers), T_prog(n), n, dtime)
          call vert_advect_tracer(Time, Adv_vel, Dens, Thickness, &
            T_prog(1:num_prog_tracers), T_prog(n), n, dtime)
        call mpp_clock_end(id_clock_tracer_advect)
        
        call mpp_clock_begin(id_clock_lap_tracer)
          call lap_tracer(Time, Thickness, T_prog(n), n)
        call mpp_clock_end(id_clock_lap_tracer)

        call mpp_clock_begin(id_clock_bih_tracer)
          call bih_tracer(Time, Thickness, T_prog(n), n)
        call mpp_clock_end(id_clock_bih_tracer)

        call mpp_clock_begin(id_clock_vert_diffuse)
          call vert_diffuse(Time, Thickness, n, T_prog(n), diff_cbt)
        call mpp_clock_end(id_clock_vert_diffuse)
       
        ! --add pme contributions to the k=1 cell.
        ! --river contributions are in T_prog%th_tendency from ocean_rivermix_mod
        !   so long as the rivermix module is enabled. 
        ! --add eta smoother contributions due to surface height smoothing.
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%th_tendency(i,j,1) = T_prog(n)%th_tendency(i,j,1) + &
               Grd%tmask(i,j,1)*(pme(i,j)*T_prog(n)%tpme(i,j) + T_prog(n)%eta_smooth(i,j) ) 
           enddo
        enddo

        ! --add bottom smoother contributions due to bottom pressure smoothing
        do j=jsc,jec
           do i=isc,iec
              kbot = Grd%kmt(i,j) 
              if(kbot>0) then
                 T_prog(n)%th_tendency(i,j,kbot) = T_prog(n)%th_tendency(i,j,kbot) + &
                                      Grd%tmask(i,j,kbot)*T_prog(n)%pbot_smooth(i,j) 
              endif
           enddo
        enddo

        ! update tracer concentration to taup1 using explicit tendency update 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec 
                 T_prog(n)%field(i,j,k,taup1) =                               &
                 (Thickness%rho_dzt(i,j,k,taum1)*T_prog(n)%field(i,j,k,taum1) &
                  + dtime*T_prog(n)%th_tendency(i,j,k))                       &
                  *Thickness%rho_dztr(i,j,k)
              enddo
           enddo
        enddo
        
        ! Add the contribution from entraining/detraining/destroyed blobs
        if (use_blobs) then
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec 
                    T_prog(n)%field(i,j,k,taup1) = T_prog(n)%field(i,j,k,taup1)       &
                         + dtime*T_prog(n)%tend_blob(i,j,k)*Thickness%rho_dztr(i,j,k)
                 enddo
              enddo
           enddo
        endif

        ! this side boundary condition does not conserve tracer...
        ! it simply inserts tracer to the domain with a specified value. 
        if(inflow_nboundary) then 
            do k=1,nk
               do j=jsc,jec
                  do i=isc,iec 
                     if(mask_inflow(i,j,k)==1.0) then 
                         T_prog(index_temp)%field(i,j,k,taup1) = temp_inflow(i,j,k)*Grd%tmask(i,j,k)
                         T_prog(index_salt)%field(i,j,k,taup1) = salt_inflow(i,j,k)*Grd%tmask(i,j,k)
                     endif
                  enddo
               enddo
            enddo
        endif

        if(have_obc) then
           call ocean_obc_tracer(T_prog(n)%field, Adv_vel%uhrho_et, Adv_vel%vhrho_nt, &
           Thickness, pme, taum1, tau, taup1, Time, n)
        endif

        if (debug_this_module) then
           write(stdoutunit,'(a)') ' '
           write(stdoutunit,*) 'From ocean_tracer_mod: tracers after explicit tendency update (taup1)'
           call tracer_prog_chksum(Time, T_prog(n), taup1)
           call tracer_min_max(Time, Thickness, T_prog(n))
        endif

        if (id_tendency_expl(n)> 0) then
           call diagnose_3d(Time, Grd, id_tendency_expl(n), T_prog(n)%th_tendency(:,:,:)*T_prog(n)%conversion)
        endif
        if (id_eta_smooth(n)> 0) then
           call diagnose_2d(Time, Grd, id_eta_smooth(n), T_prog(n)%eta_smooth(:,:)*T_prog(n)%conversion)
        endif
        if (id_eta_smooth_on_nrho(n)> 0) then
           wrk1(:,:,:) = 0.0
           k=1
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = T_prog(n)%eta_smooth(i,j)*T_prog(n)%conversion
              enddo
           enddo
           call diagnose_3d_rho(Time, Dens, id_eta_smooth_on_nrho(n), wrk1)
        end if
        if (id_pbot_smooth(n)> 0) then
           call diagnose_2d(Time, Grd, id_pbot_smooth(n), T_prog(n)%pbot_smooth(:,:)*T_prog(n)%conversion)
        endif
        call diagnose_3d(Time, Grd, id_prog_explicit(n), T_prog(n)%field(:,:,:,taup1))

     enddo   ! enddo for n=1,num_prog_tracers

     
     if (do_transport_matrix) then
       ! accumulate time explicit pieces for transport matrix 
       call transport_matrix_store_explicit(Time, T_prog, isd, ied, jsd, jed, nk, isc, iec, jsc, jec, dtts, Grd%tmask)
     endif

     ! compute ocean heating due to frazil ice formation
     ! This is when frazil was computed for CM2.0 and CM2.1 simulations at GFDL
     if(frazil_heating_before_vphysics) then 
         call mpp_clock_begin(id_clock_frazil)
         call update_ocean_density_salinity(T_prog, taup1, Dens)
         call compute_frazil_heating(Time, Thickness, Dens, &
              T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers))
         call mpp_clock_end(id_clock_frazil)
     endif

     ! compute time implicit vertical diffusion of tracer concentration 
     call mpp_clock_begin(id_clock_vert_diffuse_implicit)
       call vert_diffuse_implicit(diff_cbt, index_salt, Time, Thickness, Dens, T_prog(1:num_prog_tracers)) 
     call mpp_clock_end(id_clock_vert_diffuse_implicit)

     ! call the time implicit part of the Lagrangian blobs
     call mpp_clock_begin(id_clock_blob)
       ! Call the routines that form/move/destroy blobs implicitly in time
       call mpp_clock_begin(id_clock_blob_implicit)
         call ocean_blob_implicit(Time, Thickness, T_prog(1:num_prog_tracers), Dens, Adv_vel, Velocity, EL_diag(:))
       call mpp_clock_end(id_clock_blob_implicit)
       ! We need to update the Lagrangian system and the partition of Thickness.
       call mpp_clock_begin(id_clock_adjust_L_thickness)
         call adjust_L_thickness(Time, Thickness, T_prog, L_system)
       call mpp_clock_end(id_clock_adjust_L_thickness)
       ! We need to adjust the E system thickness to reflect the new L system thickness
       call mpp_clock_begin(id_clock_adjust_E_thickness)
         call update_E_thickness(Time, Grd, Thickness, Dens, Ext_mode)
       call mpp_clock_end(id_clock_adjust_E_thickness)
     call mpp_clock_end(id_clock_blob)

     ! convectively adjust water columns
     call mpp_clock_begin(id_clock_convection)
       call convection (Time, Thickness, T_prog(1:num_prog_tracers), Dens) 
     call mpp_clock_end(id_clock_convection)

     ! compute ocean heating due to frazil ice formation
     ! Computing frazil here ensures that the ocean temperature 
     ! at the end of a time step will never be below freezing. 
     if(frazil_heating_after_vphysics) then 
         call mpp_clock_begin(id_clock_frazil)
         call update_ocean_density_salinity(T_prog, taup1, Dens)
         call compute_frazil_heating(Time, Thickness, Dens, &
              T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers))
         call mpp_clock_end(id_clock_frazil)
     endif

     if (debug_this_module) then
         write(stdoutunit,*) 'From ocean_tracer_mod: tracers at end of update_ocean_tracer (taup1)'
         call tracer_prog_chksum(Time, T_prog(index_temp), taup1)
         call tracer_min_max(Time, Thickness, T_prog(index_temp))
         call tracer_prog_chksum(Time, T_prog(index_salt), taup1)
         call tracer_min_max(Time, Thickness, T_prog(index_salt))
   	 call tracer_prog_chksum(Time, T_prog(index_added_heat), taup1)
         call tracer_min_max(Time, Thickness, T_prog(index_added_heat))
	 call tracer_prog_chksum(Time, T_prog(index_redist_heat), taup1)
         call tracer_min_max(Time, Thickness, T_prog(index_redist_heat))
     endif

     if(have_obc) then
         do n=1,num_prog_tracers
            call mpp_update_domains(T_prog(n)%field(:,:,:,taup1), Dom%domain2d, complete=T_prog(n)%complete)
         enddo
         do n=1,num_prog_tracers
            call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taup1), 'T')
         enddo
     endif

     ! for numerically testing tracer advection schemes 
     do n=1,num_prog_tracers
        if(T_prog(n)%use_only_advection) then 
          call update_advection_only(Time, Adv_vel, Dens, Thickness, T_prog, n)
        endif 
     enddo    

    ! add tracer sources, which have dimensions tracer concentration per time 
    if(.not. zero_tracer_source) then 
        do n=1,num_prog_tracers
           do k=1,nk
              do j=jsd,jed
                 do i=isd,ied
                    T_prog(n)%field(i,j,k,taup1) = T_prog(n)%field(i,j,k,taup1) &
                                                 + dtts*T_prog(n)%source(i,j,k)
                 enddo
              enddo
           enddo
        enddo
    endif

    ! For age tracer, limit the lower bound at 0
    ! and upper bound at age_tracer_max.
    age_tracer_max = age_tracer_max_init + Time%itt0*dtts*sec_in_yr_r
    if(limit_age_tracer) then 
        do n=1,num_prog_tracers
           if(T_prog(n)%name(1:3) =='age') then 
               do k=1,nk
                  do j=jsd,jed
                     do i=isd,ied
                        T_prog(n)%field(i,j,k,taup1) = max(0.0,T_prog(n)%field(i,j,k,taup1))
                        T_prog(n)%field(i,j,k,taup1) = min(age_tracer_max,T_prog(n)%field(i,j,k,taup1))
                     enddo
                  enddo
               enddo
           endif
        enddo
    endif

    ! update squared tracer 
    call update_tracer_passive(num_prog_tracers, T_prog, T_diag, dtts, taup1)

  endif  ! endif for zero_tendency



  ! update for Dens%rho_salinity, now that all of the 
  ! updates have been applied to the prognostic salinity.  
  call update_ocean_density_salinity(T_prog, taup1, Dens)

  ! update the diagnostic temperature variable to taup1
  if(prog_temp_variable==CONSERVATIVE_TEMP) then
     T_diag(index_diag_temp)%field(:,:,:) = pottemp_from_contemp(Dens%rho_salinity(:,:,:,taup1),        &
                                                                 T_prog(index_temp)%field(:,:,:,taup1)) &
                                                                 *Grd%tmask(:,:,:)   
  else
     T_diag(index_diag_temp)%field(:,:,:) = contemp_from_pottemp(Dens%rho_salinity(:,:,:,taup1),        &
                                                                 T_prog(index_temp)%field(:,:,:,taup1)) &
                                                                 *Grd%tmask(:,:,:)
  endif

  
  ! perform range check for temperature and salinity
  if(use_tempsalt_check_range .or. debug_this_module) then
      call tempsalt_check_range(T_prog(index_salt)%field(:,:,:,taup1),&
                                T_prog(index_temp)%field(:,:,:,taup1),&
                                Thickness%depth_zt(:,:,:))
  endif 


  ! send some tracer diagnostics at time tau 
  call send_tracer_diagnostics(Time, T_prog, T_diag, Thickness, Dens, use_blobs)

  ! compute watermass diagnostics 
  call watermass_diag(Time, T_prog, T_diag, Dens, Thickness, pme)

  ! Update the combined tracer concentration to be at time taup1.
  ! Do so here so that diagnostics for fieldT are at tau. 
  if (use_blobs) then

     do n=1,num_prog_tracers
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 T_prog(n)%fieldT(i,j,k) = T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                                           + T_prog(n)%sum_blob(i,j,k,taup1)*Grd%datr(i,j)
                 T_prog(n)%fieldT(i,j,k) = Grd%tmask(i,j,k)*T_prog(n)%fieldT(i,j,k)/(epsln+Thickness%rho_dztT(i,j,k,taup1))
              enddo
           enddo
        enddo
        call mpp_update_domains(T_prog(n)%fieldT(:,:,:), Dom%domain2d, complete=T_prog(n)%complete)
      enddo

     if (debug_this_module) then
        write(stdoutunit, '(/,a)') 'Checksums for tracer fields (E and L systems)'
        do n=1,num_prog_tracers
           call write_chksum_3d('Combined '//trim(T_prog(n)%name), Grd%tmask(COMP,:)*T_prog(n)%fieldT(COMP,:))
           call write_chksum_3d('Sum of L system '//trim(T_prog(n)%name), Grd%tmask(COMP,:)*T_prog(n)%sum_blob(COMP,:,taup1))
        enddo
     endif

  endif


end subroutine update_ocean_tracer
! </SUBROUTINE> NAME="update_ocean_tracer">


!#######################################################################
! <SUBROUTINE NAME="update_advection_only">
!
! <DESCRIPTION>
!
! Redo tracer updates for those that use only advection--nothing else.
! This method is useful for testing advection schemes.
!
! T_prog(n)%use_only_advection==.true. ignores all boundary forcing
! and sources, so if T_prog(n)%stf or pme, rivers, sources
! are nonzero, tracer diagnostics will spuriously indicate 
! non-conservation.  
!
! Assume for these tests that 
!  (1) vertical advection is done fully explictly in time 
!  (2) pme, rivers, stf, btf, and other sources are zero 
!  (3) do not use advect_sweby_all
!
! </DESCRIPTION>
!
subroutine update_advection_only(Time, Adv_vel, Dens, Thickness, T_prog, n)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_thickness_type),   intent(in)    :: Thickness 
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  integer,                      intent(in)    :: n
  integer                                     :: i, j, k
  integer                                     :: taup1, taum1

  taup1 = Time%taup1  
  taum1 = Time%taum1

  T_prog(n)%th_tendency = 0.0 

  call horz_advect_tracer(Time, Adv_vel, Thickness, Dens, &
       T_prog(1:num_prog_tracers), T_prog(n), n, dtime, store_flux=.FALSE.)
  call vert_advect_tracer(Time, Adv_vel, Dens, Thickness, &
       T_prog(1:num_prog_tracers), T_prog(n), n, dtime)  

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           T_prog(n)%field(i,j,k,taup1) =                                    &
                (Thickness%rho_dzt(i,j,k,taum1)*T_prog(n)%field(i,j,k,taum1) &
               + dtime*T_prog(n)%th_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
        enddo
     enddo
  enddo

end subroutine update_advection_only
! </SUBROUTINE> NAME="update_advection_only">


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_tracer_restart(Time, T_prog, time_stamp)
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
  character(len=*),             intent(in), optional :: time_stamp
  integer :: tau, taup1, n

  tau    = Time%tau
  taup1  = Time%taup1

  do n=1,num_prog_tracers  
     if(tendency==THREE_LEVEL) then 
        call reset_field_pointer(Tracer_restart(id_tracer_type(n)), id_tracer_restart(n), T_prog(n)%field(:,:,:,tau), &
             T_prog(n)%field(:,:,:,taup1) )
     else
        call reset_field_pointer(Tracer_restart(id_tracer_type(n)), id_tracer_restart(n), T_prog(n)%field(:,:,:,taup1) )
     endif
  enddo

  do n=1, num_tracer_restart
     call save_restart(Tracer_restart(n), time_stamp)
  enddo


end subroutine ocean_tracer_restart
! </SUBROUTINE> NAME="ocean_tracer_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_end">
!
! <DESCRIPTION>
! Write ocean tracer restarts
! </DESCRIPTION>
!
subroutine ocean_tracer_end(Time, T_prog, T_diag, use_blobs)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in) :: T_diag(:)
  logical,                      intent(in) :: use_blobs

  integer            :: n, tau, taup1, unit
  integer(i8_kind)   :: chksum
  character(len=128) :: filename
  real               :: fieldsum
  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  if(.not. write_a_restart) then
     write(stdoutunit,'(/a)') '==>Warning from ocean_tracer_mod (ocean_tracer_end): NO restart written.'
     call mpp_error(WARNING,'==>Warning from ocean_tracer_mod (ocean_tracer_end): NO restart written.')
     return
  endif 

  ! open ascii file containing checksums used to check for reproducibility
  filename = 'RESTART/ocean_tracer.res'
  call mpp_open(unit, trim(filename),form=MPP_ASCII,&
     action=MPP_OVERWR,threading=MPP_SINGLE,fileset=MPP_SINGLE,nohdrs=.true.)

  call ocean_tracer_restart(Time,  T_prog)

  ! write the restart fields to the file T_prog(n)%restart_file
  ! write the chksums to the ascii file RESTART/ocean_tracer.res
  do n=1,num_prog_tracers

     if(tendency==THREE_LEVEL) then 
         write(stdoutunit,*) ' '
         write(stdoutunit,*) 'Ending ',trim(T_prog(n)%name), ' chksum (tau) ==>'
         call tracer_prog_chksum(Time, T_prog(n), tau, chksum)
         if (mpp_pe() == mpp_root_pe()) write(unit,*) trim(T_prog(n)%name)//'   tau chksum =', chksum     
     endif

     write(stdoutunit,*) ' '
     write(stdoutunit,*) 'Ending ',trim(T_prog(n)%name), ' chksum (taup1) ==>'
     call tracer_prog_chksum(Time, T_prog(n), taup1, chksum)
     if (mpp_pe() == mpp_root_pe()) write(unit,*) trim(T_prog(n)%name)//' taup1 chksum =', chksum
     if (use_blobs) then
        fieldsum = mpp_chksum(Grd%tmask(isc:iec,jsc:jec,:)*T_prog(n)%fieldT(isc:iec,jsc:jec,:))
        if (mpp_pe() == mpp_root_pe()) write(stdoutunit,*) trim(T_prog(n)%name)//' Total field  =', fieldsum
     endif

  enddo

  ! write the restart fields to the file T_diag(n)%restart_file
  ! write the chksums to the ascii file RESTART/ocean_tracer.res where writing restart
  do n=1,num_diag_tracers

    if (T_diag(n)%restart_file .ne. ' ') then  

      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'Ending ',trim(T_diag(n)%name), ' chksum (diag) ==>'
      call tracer_diag_chksum(Time, T_diag(n), chksum)
      if (mpp_pe() == mpp_root_pe()) then  
        write(unit,*) trim(T_diag(n)%name)//'       chksum =', chksum     
      endif 

    endif  

  enddo

  return

end subroutine ocean_tracer_end
! </SUBROUTINE> NAME="ocean_tracer_end">


!#######################################################################
! <SUBROUTINE NAME="compute_tmask_limit">
!
! <DESCRIPTION>
! Provide for possibility that quicker advection reverts to  
! first order upwind when tracer is outside a specified range.
! Likewise, may wish to revert neutral physics to horizontal diffusion. 
!
! For this purpose, we define a mask which is set to unity where 
! fluxes revert to first order upwind advection 
! (if using quicker) and horizontal diffusion (if using neutral). 
!
! This method is very ad hoc.  What is preferred for advection is to use
! a monotonic scheme, such as mdfl_sweby or mdppm.  For neutral physics,  
! no analogous monotonic scheme has been implemented.  Such could be 
! useful, especially for passive tracers. In the meantime, tmask_limit 
! provides a very rough limiter for neutral physics to help keep tracers 
! within specified bounds.
!
! </DESCRIPTION>
!
subroutine compute_tmask_limit(Time, T_prog)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: tau
  integer :: i, j, k, n 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. compute_tmask_limit_on) return 

  tau = Time%tau

  do n=1,num_prog_tracers 

     do k=1,nk 
        do j=jsc-1,jec
           do i=isc-1,iec

              T_prog(n)%tmask_limit(i,j,k) = 0.0   

              ! is tracer outside range? 
              if(Grd%tmask(i,j,k) == 1.0) then 
                  if(T_prog(n)%field(i,j,k,tau) < T_prog(n)%min_tracer_limit .or. &
                     T_prog(n)%field(i,j,k,tau) > T_prog(n)%max_tracer_limit) then 
                      T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif

              ! is east tracer outside range? 
              if(Grd%tmask(i+1,j,k) == 1.0) then 
                  if(T_prog(n)%field(i+1,j,k,tau) < T_prog(n)%min_tracer_limit .or.  &
                     T_prog(n)%field(i+1,j,k,tau) > T_prog(n)%max_tracer_limit) then 
                      T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif

              ! is north tracer outside range? 
              if(Grd%tmask(i,j+1,k) == 1.0) then 
                  if(T_prog(n)%field(i,j+1,k,tau) < T_prog(n)%min_tracer_limit .or.  &
                     T_prog(n)%field(i,j+1,k,tau) > T_prog(n)%max_tracer_limit) then 
                      T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif

           enddo
        enddo
     enddo

     ! is deeper tracer outside range? 
     do k=1,nk-1 
        do j=jsc-1,jec
           do i=isc-1,iec
              if(Grd%tmask(i,j,k+1) == 1.0) then 
                  if(T_prog(n)%field(i,j,k+1,tau) < T_prog(n)%min_tracer_limit .or.  &
                     T_prog(n)%field(i,j,k+1,tau) > T_prog(n)%max_tracer_limit) then 
                      T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif
           enddo
        enddo
     enddo

  enddo ! end of loop n=1,num_prog_tracers  


  ! insist that if temperature or salinity tmask_limit is 1.0, then both must be 1.0
  if(tmask_limit_ts_same) then 
      do k=1,nk 
         do j=jsc-1,jec
            do i=isc-1,iec
               if(       T_prog(index_temp)%tmask_limit(i,j,k)==1.0 &
                    .or. T_prog(index_salt)%tmask_limit(i,j,k)==1.0) then 
                   T_prog(index_temp)%tmask_limit(i,j,k)=1.0 
                   T_prog(index_salt)%tmask_limit(i,j,k)=1.0   
               endif
            enddo
         enddo
      enddo
  endif


  ! debugging and diagnostics  
  do n=1,num_prog_tracers 

     if(debug_this_module) then 
         write(stdoutunit,*) ' ' 
         call write_timestamp(Time%model_time)
         call write_chksum_3d('tmask_limit', T_prog(n)%tmask_limit(COMP,:))
     endif

     call diagnose_3d(Time, Grd, id_tmask_limit(n), T_prog(n)%tmask_limit(:,:,:))

  enddo


end subroutine compute_tmask_limit
! </SUBROUTINE> NAME="compute_tmask_limit"



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
! Use rho_dzt weighting to account for nonBoussinesq.  
!
! Author: Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine remap_s_to_depth(Thickness, Time, array_in, ntracer)

  type(ocean_thickness_type),    intent(in) :: Thickness
  type(ocean_time_type),         intent(in) :: Time
  real,  dimension(isd:,jsd:,:), intent(in) :: array_in
  integer,                       intent(in) :: ntracer

  integer :: i, j, k, kk, tau

  tau  = Time%tau
  wrk1 = 0.0
  wrk2 = 0.0
  wrk3 = 0.0

  k=1
  do kk=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Thickness%depth_zt(i,j,kk) < Grd%zw(k)) then 
              wrk1(i,j,k) = wrk1(i,j,k) + array_in(i,j,kk)*Thickness%rho_dzt(i,j,kk,tau)
              wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,kk,tau)
           endif 
        enddo
     enddo
  enddo

  do k=1,nk-1
     do kk=1,nk
        do j=jsd,jed
           do i=isd,ied
              if(Grd%zw(k) <= Thickness%depth_zt(i,j,kk) .and. Thickness%depth_zt(i,j,kk) < Grd%zw(k+1)) then 
                 wrk1(i,j,k+1) = wrk1(i,j,k+1) + array_in(i,j,kk)*Thickness%rho_dzt(i,j,kk,tau)
                 wrk2(i,j,k+1) = wrk2(i,j,k+1) + Thickness%rho_dzt(i,j,kk,tau)
              endif
           enddo
        enddo
     enddo
  enddo

  k=nk
  do kk=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Thickness%depth_zt(i,j,kk) > Grd%zw(k) ) then 
               wrk1(i,j,k) = wrk1(i,j,k) + array_in(i,j,kk)*Thickness%rho_dzt(i,j,kk,tau)
               wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,kk,tau)
           endif
        enddo
     enddo
  enddo

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = Grd%tmask_depth(i,j,k)*wrk1(i,j,k)/(wrk2(i,j,k)+epsln)
        enddo
     enddo
  enddo

  used = send_data (id_prog_on_depth(ntracer), wrk3(:,:,:), &
  Time%model_time,rmask=Grd%tmask_depth(:,:,:),             &
  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


end subroutine remap_s_to_depth
! </SUBROUTINE>  NAME="remap_s_to_depth"



!#######################################################################
! <SUBROUTINE NAME="remap_depth_to_s">
! <DESCRIPTION>
!
! Remap in the vertical from depth to s-coordinate.  This routine is 
! used for initializing terrain following coordinate models given an 
! initial tracer field generated assuming depth-like vertical coordinate. 
!
! This routine is of use for terrain following vertical coordinates,
! which generally deviate a lot from z, zstar, pressure, or pstar 
! coordinates.  
!
! Algorithm is very rudimentary and can be made better.  
!
! Author: Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine remap_depth_to_s(Thickness, Time, array)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_time_type),        intent(in)    :: Time
  real, dimension(isd:,jsd:,:), intent(inout) :: array

  integer :: i, j, k, kk

  integer :: stdoutunit 
  stdoutunit=stdout() 

  wrk1 = 0.0

  k=1
  do kk=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Thickness%depth_zt(i,j,kk) < Grd%zw(k)) then 
              wrk1(i,j,kk) = array(i,j,k)
           endif  
        enddo
     enddo
  enddo

  do k=1,nk-1
     do kk=1,nk
        do j=jsd,jed
           do i=isd,ied
              if(Grd%zw(k) <= Thickness%depth_zt(i,j,kk) .and. Thickness%depth_zt(i,j,kk) < Grd%zw(k+1)) then 
                 wrk1(i,j,kk) = array(i,j,k+1)
              endif
           enddo
        enddo
     enddo
  enddo

  k=nk
  do kk=1,nk
     do j=jsd,jed
        do i=isd,ied
           if(Thickness%depth_zt(i,j,kk) > Grd%zw(k) ) then 
               wrk1(i,j,kk) = array(i,j,k)
           endif
        enddo
     enddo
  enddo

  if(debug_this_module) then 
      write(stdoutunit,'(/a)')'==>from remap_depth_to_s'
      call write_timestamp(Time%model_time)
      do k=1,nk
         do j=jsc,jsc
            do i=isc,isc
               write(stdoutunit,'(a,i3,a,i3,a,i3,a,f12.4)') &
               'array_in (',i+Dom%ioff,',',j+Dom%joff,',',k,') = ',array(i,j,k)
               write(stdoutunit,'(a,i3,a,i3,a,i3,a,f12.4)') &
               'array_out(',i+Dom%ioff,',',j+Dom%joff,',',k,') = ',wrk1(i,j,k)
            enddo
         enddo
      enddo
  endif

  array(:,:,:) = wrk1(:,:,:)


end subroutine remap_depth_to_s
! </SUBROUTINE>  NAME="remap_depth_to_s"


!#######################################################################
! <SUBROUTINE NAME="inflow_nboundary_init">
!
! <DESCRIPTION>
! Initialize the mask, temp, and salt values at the northern boundary
! where we are specifying values for an inflow boundary condition.
!
! </DESCRIPTION>
!
subroutine inflow_nboundary_init

  character(len=128) :: filename

  integer :: stdoutunit 
  stdoutunit=stdout() 

  allocate(mask_inflow(isd:ied,jsd:jed,nk))
  allocate(temp_inflow(isd:ied,jsd:jed,nk))
  allocate(salt_inflow(isd:ied,jsd:jed,nk))

  write(stdoutunit,'(a)') &
  '==>Note: running with inflow_nboundary. Specify tracer at north boundary. This is nonconservative.'  

  mask_inflow(:,:,:) = 0.0
  temp_inflow(:,:,:) = 0.0
  salt_inflow(:,:,:) = 0.0

  filename = 'INPUT/mask_inflow.nc'
  if (file_exist(trim(filename))) then 
      call read_data(filename,'mask_inflow',mask_inflow,Dom%domain2d,timelevel=1)
      call mpp_update_domains(mask_inflow(:,:,:), Dom%domain2d)
  else 
     call mpp_error(FATAL,&
    '==>Error from ocean_tracer_mod: inflow_nboundary_init cannot find mask_inflow file.')
  endif 

  filename = 'INPUT/temp_inflow.nc'
  if (file_exist(trim(filename))) then 
      call read_data(filename,'temp_inflow',temp_inflow,Dom%domain2d,timelevel=1)
      call mpp_update_domains(temp_inflow(:,:,:), Dom%domain2d)
  else 
     call mpp_error(FATAL,&
    '==>Error from ocean_tracer_mod: inflow_nboundary_init cannot find temp_inflow file.')
  endif 

  filename = 'INPUT/salt_inflow.nc'
  if (file_exist(trim(filename))) then 
      call read_data(filename,'salt_inflow',salt_inflow,Dom%domain2d,timelevel=1)
      call mpp_update_domains(salt_inflow(:,:,:), Dom%domain2d)
  else 
     call mpp_error(FATAL,&
    '==>Error from ocean_tracer_mod: inflow_nboundary_init cannot find salt_inflow file.')
  endif 


end subroutine inflow_nboundary_init
! </SUBROUTINE> NAME="inflow_nboundary_init"

    
!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_diagnostics_init">
!
! <DESCRIPTION>
! Do some initialization required for diagnostics. 
! Density must already be initialized for this routine to be used.  
! Hence, this routine is called from ocean_model.F90.
! </DESCRIPTION>
!
subroutine ocean_tracer_diagnostics_init(Time, Dens, T_diag, T_prog, ver_coordinate_class)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_diag_tracer_type), intent(in) :: T_diag(:)
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: ver_coordinate_class
  
  integer :: n
  integer :: stdoutunit
  stdoutunit=stdout()

  compute_watermass_diag = .false. 
  vert_coordinate_class = ver_coordinate_class 

  neutralrho_nk = size(Dens%neutralrho_ref(:))
  allocate( nrho_work(isd:ied,jsd:jed,neutralrho_nk) )
  allocate( nrho_work2(isd:ied,jsd:jed,neutralrho_nk) )
  allocate( nrho_mass(isd:ied,jsd:jed,neutralrho_nk) )
  nrho_work(:,:,:)  = 0.0  
  nrho_work2(:,:,:) = 0.0  
  nrho_mass(:,:,:)  = 0.0  

  do n=1,num_diag_tracers
     if (T_diag(n)%name == 'frazil') index_frazil = n
  enddo

  do n=1,num_prog_tracers

     ! Neutral density binning of heat budget diagnostics
     id_tendency_on_nrho(n) = register_diag_field ('ocean_model',            &
          trim(T_prog(n)%name)//'_tendency_on_nrho', Dens%neutralrho_axes(1:3),   &
          Time%model_time, 'time tendency for tracer '//trim(T_prog(n)%longname)//" binned to neutral density", &
          trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e20,1.e20/))
     id_eta_smooth_on_nrho(n) = register_diag_field ('ocean_model',                    &
          trim(T_prog(n)%name)//'_eta_smooth_on_nrho', Dens%neutralrho_axes(1:3),           &
          Time%model_time, 'surface smoother for ' // trim(T_prog(n)%longname)//" binned to neutral density", &
          trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e20,1.e20/))
  enddo

  ! time tendency 
  id_neut_rho_tendency = register_diag_field ('ocean_model','neut_rho_tendency',      &  
    Grd%tracer_axes(1:3), Time%model_time,'total time tendency of locally ref potrho',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_tendency > 0) compute_watermass_diag=.true.

  id_neut_rho_tendency_on_nrho = register_diag_field ('ocean_model',         &
    'neut_rho_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
    'total time tendency of locally ref potrho as binned to neutral density',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_tendency_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_tendency = register_diag_field ('ocean_model','wdian_rho_tendency',& 
    Grd%tracer_axes(1:3), Time%model_time,                                        &
    'dianeutral mass transport due to time tendency of locally ref potrho',       &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_tendency > 0) compute_watermass_diag=.true.

  id_wdian_rho_tendency_on_nrho = register_diag_field ('ocean_model',                   &
    'wdian_rho_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
    'dianeutral mass transport due to total time tendency as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_tendency_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_tendency = register_diag_field ('ocean_model',                     &
    'tform_rho_tendency', Grd%tracer_axes(1:3), Time%model_time,                  &
    'water mass transformation from total tendency on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_tendency > 0) compute_watermass_diag=.true.

  id_tform_rho_tendency_on_nrho = register_diag_field ('ocean_model',            &
    'tform_rho_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,    &
    'water mass transformation from total tendency as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_tendency_on_nrho > 0) compute_watermass_diag=.true.


  ! time tendency from temperature effects 
  id_neut_temp_tendency = register_diag_field ('ocean_model','neut_temp_tendency',           &  
    Grd%tracer_axes(1:3), Time%model_time,'temp related time tendency of locally ref potrho',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_tendency > 0) compute_watermass_diag=.true.

  id_neut_temp_tendency_on_nrho = register_diag_field ('ocean_model',             &
    'neut_temp_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'temp related time tendency of locally ref potrho  binned to neutral density',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_tendency_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_temp_tendency = register_diag_field ('ocean_model','wdian_temp_tendency',    &   
    Grd%tracer_axes(1:3), Time%model_time,                                              &
    'temp related dianeutral mass transport due to time tendency of locally ref potrho',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_tendency > 0) compute_watermass_diag=.true.

  id_wdian_temp_tendency_on_nrho = register_diag_field ('ocean_model',                         &
    'wdian_temp_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
    'temp related dianeutral mass transport due to time tendency as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_tendency_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_temp_tendency = register_diag_field ('ocean_model',       &
    'tform_temp_tendency', Grd%tracer_axes(1:3), Time%model_time,    &
    'temp related water mass transformation from tendency on levels',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_tendency > 0) compute_watermass_diag=.true.

  id_tform_temp_tendency_on_nrho = register_diag_field ('ocean_model',               &
    'tform_temp_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
    'temp related water mass transformation from tendency binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_tendency_on_nrho > 0) compute_watermass_diag=.true.


  ! time tendency from salinity effects 
  id_neut_salt_tendency = register_diag_field ('ocean_model','neut_salt_tendency',           &  
    Grd%tracer_axes(1:3), Time%model_time,'salt related time tendency of locally ref potrho',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_tendency > 0) compute_watermass_diag=.true.

  id_neut_salt_tendency_on_nrho = register_diag_field ('ocean_model',             &
    'neut_salt_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'salt related time tendency of locally ref potrho  binned to neutral density',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_tendency_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_salt_tendency = register_diag_field ('ocean_model','wdian_salt_tendency',    &   
    Grd%tracer_axes(1:3), Time%model_time,                                              &
    'salt related dianeutral mass transport due to time tendency of locally ref potrho',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_tendency > 0) compute_watermass_diag=.true.

  id_wdian_salt_tendency_on_nrho = register_diag_field ('ocean_model',                         &
    'wdian_salt_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
    'salt related dianeutral mass transport due to time tendency as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_tendency_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_salt_tendency = register_diag_field ('ocean_model',       &
    'tform_salt_tendency', Grd%tracer_axes(1:3), Time%model_time,    &
    'salt related water mass transformation from tendency on levels',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_tendency > 0) compute_watermass_diag=.true.

  id_tform_salt_tendency_on_nrho = register_diag_field ('ocean_model',               &
    'tform_salt_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
    'salt related water mass transformation from tendency binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_tendency_on_nrho > 0) compute_watermass_diag=.true.



  ! smoothing operator 
  id_neut_rho_smooth = register_diag_field ('ocean_model','neut_rho_smooth',&  
    Grd%tracer_axes(1:3), Time%model_time,                                  &
    'time tendency of local ref potrho from eta/pbot smoother',             &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_smooth > 0) compute_watermass_diag=.true.

  id_neut_rho_smooth_on_nrho = register_diag_field ('ocean_model',                      &
   'neut_rho_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'time tendency of local ref potrho from eta/pbot smoother binned to neutral density',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_smooth = register_diag_field ('ocean_model','wdian_rho_smooth',& 
    Grd%tracer_axes(1:3), Time%model_time,                                    &
    'dianeutral mass transport due to eta/pbot smoother',                     &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_smooth > 0) compute_watermass_diag=.true.

  id_wdian_rho_smooth_on_nrho = register_diag_field ('ocean_model',                &
    'wdian_rho_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,        &
    'dianeutral mass transport due to eta/pbot smoother binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_smooth = register_diag_field ('ocean_model',       &
   'tform_rho_smooth', Dens%neutralrho_axes(1:3), Time%model_time,&
   'water mass transformation from smoother on levels',           &
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_smooth_on_nrho = register_diag_field ('ocean_model',       &
   'tform_rho_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
   'water mass transformation from smoother binned to neutral density',   &
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_eta_tend_smooth = register_diag_field ('ocean_model','eta_tend_smooth',&  
    Grd%tracer_axes(1:2), Time%model_time,                                  &
    'non-Bouss steric sea level tendency from eta/pbot smoother',           &
    'm/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_eta_tend_smooth > 0) compute_watermass_diag=.true.

  id_eta_tend_smooth_glob = register_diag_field ('ocean_model',              &
    'eta_tend_smooth_glob', Time%model_time,                                 &
    'global mean non-Bouss steric sea level tendency from eta/pbot smoother',&
    'm/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_eta_tend_smooth_glob > 0) compute_watermass_diag=.true.


  ! smoothing operator from temperature effects 
  id_neut_temp_smooth = register_diag_field ('ocean_model','neut_temp_smooth',&  
    Grd%tracer_axes(1:3), Time%model_time,                                    &
    'temp related time tendency of local ref potrho from eta/pbot smoother',  &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_smooth > 0) compute_watermass_diag=.true.

  id_neut_temp_smooth_on_nrho = register_diag_field ('ocean_model',                                  &
   'neut_temp_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                           &
   'temp related time tendency of local ref potrho from eta/pbot smoother binned to neutral density',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_temp_smooth = register_diag_field ('ocean_model','wdian_temp_smooth',& 
    Grd%tracer_axes(1:3), Time%model_time,                                      &
    'temp related dianeutral mass transport due to eta/pbot smoother',          &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_smooth > 0) compute_watermass_diag=.true.

  id_wdian_temp_smooth_on_nrho = register_diag_field ('ocean_model',                            &
    'wdian_temp_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                    &
    'temp related dianeutral mass transport due to eta/pbot smoother binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_temp_smooth = register_diag_field ('ocean_model',        &
   'tform_temp_smooth', Dens%neutralrho_axes(1:3), Time%model_time, &
   'temp related water mass transformation from smoother on levels',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_temp_smooth_on_nrho = register_diag_field ('ocean_model',                &
   'tform_temp_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
   'temp related water mass transformation from smoother binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_smooth_on_nrho > 0) compute_watermass_diag=.true.


  ! smoothing operator from salinity effects 
  id_neut_salt_smooth = register_diag_field ('ocean_model','neut_salt_smooth',&  
    Grd%tracer_axes(1:3), Time%model_time,                                    &
    'salt related time tendency of local ref potrho from eta/pbot smoother',  &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_smooth > 0) compute_watermass_diag=.true.

  id_neut_salt_smooth_on_nrho = register_diag_field ('ocean_model',                                  &
   'neut_salt_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                           &
   'salt related time tendency of local ref potrho from eta/pbot smoother binned to neutral density',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_salt_smooth = register_diag_field ('ocean_model','wdian_salt_smooth',& 
    Grd%tracer_axes(1:3), Time%model_time,                                      &
    'salt related dianeutral mass transport due to eta/pbot smoother',          &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_smooth > 0) compute_watermass_diag=.true.

  id_wdian_salt_smooth_on_nrho = register_diag_field ('ocean_model',                            &
    'wdian_salt_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                    &
    'salt related dianeutral mass transport due to eta/pbot smoother binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_salt_smooth = register_diag_field ('ocean_model',        &
   'tform_salt_smooth', Dens%neutralrho_axes(1:3), Time%model_time, &
   'salt related water mass transformation from smoother on levels',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_smooth_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_salt_smooth_on_nrho = register_diag_field ('ocean_model',                &
   'tform_salt_smooth_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
   'salt related water mass transformation from smoother binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_smooth_on_nrho > 0) compute_watermass_diag=.true.

  ! contribution of frazil to heat budget binned on neutral density
  id_frazil_on_nrho = register_diag_field ('ocean_model', 'frazil_on_nrho', Dens%neutralrho_axes(1:3), &
       Time%model_time, 'ocn frazil heat flux over time step binned to neutral density', 'W/m^2',&
       missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_frazil_on_nrho > 0) compute_watermass_diag=.true.

  ! contributions from frazil formation 
  id_neut_rho_frazil = register_diag_field ('ocean_model', 'neut_rho_frazil',&
    Grd%tracer_axes(1:3), Time%model_time,                                   &
    'update of local ref potrho from frazil formation',                      &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_frazil > 0) compute_watermass_diag=.true.

  id_neut_rho_frazil_on_nrho = register_diag_field ('ocean_model',                  &
    'neut_rho_frazil_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
    'update of local ref potrho from frazil formation as binned to neutral density',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_frazil_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_frazil = register_diag_field ('ocean_model', 'wdian_rho_frazil',&
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'dianeutral mass transport due to frazil formation',                       &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_frazil > 0) compute_watermass_diag=.true.

  id_wdian_rho_frazil_on_nrho = register_diag_field ('ocean_model',                  &
    'wdian_rho_frazil_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
    'dianeutral mass transport due to frazil formation as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_frazil_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_frazil_on_nrho = register_diag_field ('ocean_model',        &
    'tform_rho_frazil_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
    'water mass transformation from frazil as binned to neutral density',  &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_frazil_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_frazil = register_diag_field ('ocean_model',               &
    'tform_rho_frazil', Dens%neutralrho_axes(1:3), Time%model_time,       &
    'water mass transformation from frazil on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_frazil > 0) compute_watermass_diag=.true.

  id_eta_tend_frazil = register_diag_field ('ocean_model','eta_tend_frazil',&  
    Grd%tracer_axes(1:2), Time%model_time,                                  &
    'non-Bouss steric sea level tendency from frazil',                      &
    'm/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_eta_tend_frazil > 0) compute_watermass_diag=.true.

  id_eta_tend_frazil_glob = register_diag_field ('ocean_model',   &
    'eta_tend_frazil_glob', Time%model_time,                      &
    'global mean non-Bouss steric sea level tendency from frazil',&
    'm/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_eta_tend_frazil_glob > 0) compute_watermass_diag=.true.

  ! contributions from pme in flux-form 
  id_neut_rho_pme = register_diag_field ('ocean_model','neut_rho_pme',                 & 
    Grd%tracer_axes(1:3), Time%model_time,'time tendency of local ref potrho from pme',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pme > 0) compute_watermass_diag=.true.

  id_neut_rho_pme_on_nrho = register_diag_field ('ocean_model',               &
    'neut_rho_pme_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
    'time tendency of local ref potrho from pme as binned to neutral density',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pme_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_pme = register_diag_field ('ocean_model','wdian_rho_pme',& 
    Grd%tracer_axes(1:3), Time%model_time,                              &
    'dianeutral mass transport due to pme',                             &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pme > 0) compute_watermass_diag=.true.

  id_wdian_rho_pme_on_nrho = register_diag_field ('ocean_model',         &
    'wdian_rho_pme_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
    'dianeutral mass transport due to pme as binned to neutral density', &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pme_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_pme = register_diag_field ('ocean_model',               &
    'tform_rho_pme', Dens%neutralrho_axes(1:3), Time%model_time,       &
    'water mass transformation from pme on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pme > 0) compute_watermass_diag=.true.

  id_tform_rho_pme_on_nrho = register_diag_field ('ocean_model',        &
    'tform_rho_pme_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
    'water mass transformation from pme as binned to neutral density',  &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pme_on_nrho > 0) compute_watermass_diag=.true.



  ! contributions from pme in advective form using kinematic based method 
  id_neut_rho_pbl_pme_kn = register_diag_field ('ocean_model',       &
       'neut_rho_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time, &
       'kinematic advective-form material time derivative from pme', &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_pme_kn > 0) compute_watermass_diag=.true.

  id_neut_rho_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                           &
       'neut_rho_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                &
       'kinematic advective-form material time derivative from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_pbl_pme_kn = register_diag_field ('ocean_model',       &
       'wdian_rho_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time, &
       'kinematic advective-form dianeutral transport from pme',      &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_pme_kn > 0) compute_watermass_diag=.true.

  id_wdian_rho_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                      &
       'wdian_rho_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
       'kinematic advective-form dianeutral transport from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_pbl_pme_kn = register_diag_field ('ocean_model',                     &
       'tform_rho_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,               &
       'kinematic advective-form transform from pme pre-binning to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_pme_kn > 0) compute_watermass_diag=.true.

  id_tform_rho_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',            &
       'tform_rho_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
       'kinematic advective-form  transform from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.


  ! contributions from pme temp pieces in advective form using kinematic based method 
  id_neut_temp_pbl_pme_kn = register_diag_field ('ocean_model',                  &
       'neut_temp_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,            &
       'temp related kinematic advective-form material time derivative from pme',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_pme_kn > 0) compute_watermass_diag=.true.

  id_neut_temp_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                                       &
       'neut_temp_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                            &
       'temp related kinematic advective-form material time derivative from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.

  idwdian_temp_pbl_pme_kn = register_diag_field ('ocean_model',              &
       'wdian_temp_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,       &
       'temp related kinematic advective-form dianeutral transport from pme',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_temp_pbl_pme_kn > 0) compute_watermass_diag=.true.

  idwdian_temp_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                                   &
       'wdian_temp_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
       'temp related kinematic advective-form dianeutral transport from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_temp_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.

  idtform_temp_pbl_pme_kn = register_diag_field ('ocean_model',                                  &
       'tform_temp_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,                           &
       'temp related kinematic advective-form transform from pme pre-binning to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_temp_pbl_pme_kn > 0) compute_watermass_diag=.true.

  idtform_temp_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                        &
       'tform_temp_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
       'temp related kinematic advective-form transform from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_temp_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.


  ! contributions from pme salt pieces in advective form using kinematic based method 
  id_neut_salt_pbl_pme_kn = register_diag_field ('ocean_model',                  &
       'neut_salt_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,            &
       'salt related kinematic advective-form material time derivative from pme',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_pme_kn > 0) compute_watermass_diag=.true.

  id_neut_salt_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                                       &
       'neut_salt_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                            &
       'salt related kinematic advective-form material time derivative from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.

  idwdian_salt_pbl_pme_kn = register_diag_field ('ocean_model',              &
       'wdian_salt_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,       &
       'salt related kinematic advective-form dianeutral transport from pme',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_salt_pbl_pme_kn > 0) compute_watermass_diag=.true.

  idwdian_salt_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                                   &
       'wdian_salt_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
       'salt related kinematic advective-form dianeutral transport from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_salt_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.

  idtform_salt_pbl_pme_kn = register_diag_field ('ocean_model',                                  &
       'tform_salt_pbl_pme_kn', Grd%tracer_axes(1:3), Time%model_time,                           &
       'salt related kinematic advective-form transform from pme pre-binning to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_salt_pbl_pme_kn > 0) compute_watermass_diag=.true.

  idtform_salt_pbl_pme_kn_on_nrho = register_diag_field ('ocean_model',                        &
       'tform_salt_pbl_pme_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
       'salt related kinematic advective-form transform from pme as binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_salt_pbl_pme_kn_on_nrho > 0) compute_watermass_diag=.true.




  ! contributions from pme in advective form using process based method 
  id_neut_rho_pbl_pme_pr = register_diag_field ('ocean_model',    &
    'neut_rho_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time, &
    'process advective-form material time derivative from pme',   &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_pme_pr > 0) compute_watermass_diag=.true.

  id_neut_rho_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                      &
    'neut_rho_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
    'process advective-form material time derivative from pme as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_pbl_pme_pr = register_diag_field ('ocean_model',    &
    'wdian_rho_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time, &
    'process advective-form dianeutral transport from pme',        &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_pme_pr > 0) compute_watermass_diag=.true.

  id_wdian_rho_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                 &
    'wdian_rho_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
    'process advective-form dianeutral transport from pme as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_pbl_pme_pr = register_diag_field ('ocean_model',                           &
    'tform_rho_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,                        &
    'process advective-form water mass transform from pme pre-binning to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_pme_pr > 0) compute_watermass_diag=.true.

  id_tform_rho_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                 &
       'tform_rho_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,      &
       'process advective-form water mass transform from pme binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.



  ! contributions from pme temp terms in advective form using process based method 
  id_neut_temp_pbl_pme_pr = register_diag_field ('ocean_model',             &
    'neut_temp_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,          &
    'temp related process advective-form material time derivative from pme',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_pme_pr > 0) compute_watermass_diag=.true.

  id_neut_temp_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                               &
    'neut_temp_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'temp related process advective-form material time derivative from pme binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.

  idwdian_temp_pbl_pme_pr = register_diag_field ('ocean_model',         &
    'wdian_temp_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,     &
    'temp related process advective-form dianeutral transport from pme',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_temp_pbl_pme_pr > 0) compute_watermass_diag=.true.

  idwdian_temp_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                           &
    'wdian_temp_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'temp related process advective-form dianeutral transport from pme binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_temp_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.

  idtform_temp_pbl_pme_pr = register_diag_field ('ocean_model',         &
    'tform_temp_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,     &
    'temp related process advective-form water mass transform from pme',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_temp_pbl_pme_pr > 0) compute_watermass_diag=.true.

  idtform_temp_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                           &
    'tform_temp_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'temp related process advective-form water mass transform from pme binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_temp_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.



  ! contributions from pme salinity terms in advective form using process based method 
  id_neut_salt_pbl_pme_pr = register_diag_field ('ocean_model',             &
    'neut_salt_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,          &
    'salt related process advective-form material time derivative from pme',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_pme_pr > 0) compute_watermass_diag=.true.

  id_neut_salt_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                               &
    'neut_salt_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'salt related process advective-form material time derivative from pme binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.

  idwdian_salt_pbl_pme_pr = register_diag_field ('ocean_model',         &
    'wdian_salt_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,     &
    'salt related process advective-form dianeutral transport from pme',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_salt_pbl_pme_pr > 0) compute_watermass_diag=.true.

  idwdian_salt_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                           &
    'wdian_salt_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'salt related process advective-form dianeutral transport from pme binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idwdian_salt_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.

  idtform_salt_pbl_pme_pr = register_diag_field ('ocean_model',         &
    'tform_salt_pbl_pme_pr', Grd%tracer_axes(1:3), Time%model_time,     &
    'salt related process advective-form water mass transform from pme',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_salt_pbl_pme_pr > 0) compute_watermass_diag=.true.

  idtform_salt_pbl_pme_pr_on_nrho = register_diag_field ('ocean_model',                           &
    'tform_salt_pbl_pme_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'salt related process advective-form water mass transform from pme binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(idtform_salt_pbl_pme_pr_on_nrho > 0) compute_watermass_diag=.true.



  if(vert_coordinate_class==DEPTH_BASED) then
      id_mass_t_on_nrho = register_diag_field ('ocean_model',                &
           'mass_t_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
           'ocean t-cell volume (*rho0) as binned to neutral density layers',&
           'kg', missing_value=missing_value, range=(/-1.e20,1.e20/))
      if(id_mass_t_on_nrho > 0) compute_watermass_diag=.true.

      id_mass_t_tendency_on_nrho = register_diag_field ('ocean_model',                   &
           'mass_t_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,        &
           'tendency of ocean t-cell volume (*rho0) as binned to neutral density layers',&
           'kg/s', missing_value=missing_value, range=(/-1.e20,1.e20/))
      if(id_mass_t_tendency_on_nrho > 0) compute_watermass_diag=.true.

      id_mass_nrho_layer = register_diag_field ('ocean_model', 'mass_nrho_layer',                         &
           Dens%neutralrho_axes(1:3), Time%model_time,                                                    &
          'volume of a neutral rho layer (*rho0) computed via interpolation to get layer thickness', 'kg',&
           missing_value=missing_value, range=(/0.0,1.e20/))   
      if(id_mass_nrho_layer > 0) compute_watermass_diag=.true.

      id_mass_nrho_tendency_layer = register_diag_field ('ocean_model',                                          &
           'mass_nrho_tendency_layer', Dens%neutralrho_axes(1:3), Time%model_time,                               &
           'tendency of volume of a neutral rho layer (*rho0) computed via interpolation to get layer thickness',&
           'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
      if(id_mass_nrho_tendency_layer > 0) compute_watermass_diag=.true.

  else 

      id_mass_t_on_nrho = register_diag_field ('ocean_model',           &
           'mass_t_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
           'ocean t-cell volume as binned to neutral density layers',   &
           'kg', missing_value=missing_value, range=(/-1.e20,1.e20/))
      if(id_mass_t_on_nrho > 0) compute_watermass_diag=.true.

      id_mass_t_tendency_on_nrho = register_diag_field ('ocean_model',            &
           'mass_t_tendency_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
           'tendency of ocean t-cell volume as binned to neutral density layers', &
           'kg/s', missing_value=missing_value, range=(/-1.e20,1.e20/))
      if(id_mass_t_tendency_on_nrho > 0) compute_watermass_diag=.true.

      id_mass_nrho_layer = register_diag_field ('ocean_model', 'mass_nrho_layer',                 &
           Dens%neutralrho_axes(1:3), Time%model_time,                                            &
          'volume of a neutral rho layer computed via interpolation to get layer thickness', 'kg',&
           missing_value=missing_value, range=(/0.0,1.e20/))   
      if(id_mass_nrho_layer > 0) compute_watermass_diag=.true.

      id_mass_nrho_tendency_layer = register_diag_field ('ocean_model',                                  &
           'mass_nrho_tendency_layer', Dens%neutralrho_axes(1:3), Time%model_time,                       &
           'tendency of volume of a neutral rho layer computed via interpolation to get layer thickness',&
           'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
      if(id_mass_nrho_tendency_layer > 0) compute_watermass_diag=.true.

  endif


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_tracer_mod w/ compute_watermass_diag=.true.'  
  endif 



end subroutine ocean_tracer_diagnostics_init
! </SUBROUTINE> NAME="ocean_tracer_diagnostics_init"

!#######################################################################
! <SUBROUTINE NAME="send_tracer_diagnostics">
!
! <DESCRIPTION>
!
! For sending some tracer diagnostics
! 
! </DESCRIPTION>
!
subroutine send_tracer_diagnostics(Time, T_prog, T_diag, Thickness, Dens, use_blobs)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(in)    :: T_diag(:)
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  logical,                        intent(in)    :: use_blobs 

  integer :: i,j,k,kbot,n
  integer :: taum1,tau,taup1
  real    :: total_tracer

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  do n=1,num_prog_tracers

     if(id_prog(n) > 0) then
        if( n==index_temp) then
           wrk1(:,:,:) = 0.0
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = T_prog(n)%field(i,j,k,tau) + Grd%tmask(i,j,k)*cmip_offset 
                 enddo
              enddo
           enddo
           call diagnose_3d(Time, Grd, id_prog(n), wrk1(:,:,:))
        else
           call diagnose_3d(Time, Grd, id_prog(n), T_prog(n)%field(:,:,:,tau))
         endif
     endif  

     if(id_prog_rhodzt(n)  > 0) then
        call diagnose_3d(Time, Grd, id_prog_rhodzt(n), Thickness%rho_dzt(:,:,:,tau)*T_prog(n)%field(:,:,:,tau))
     endif

     ! vertical sum of tracer*rho_dzt 
     if(id_prog_int_rhodz(n) > 0) then 
         wrk1_2d(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%field(i,j,k,tau) 
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_prog_int_rhodz(n), wrk1_2d(:,:))
     endif

     if (id_prog_on_depth(n) > 0)  then 
        call remap_s_to_depth(Thickness, Time, T_prog(n)%field(:,:,:,tau), n)
     endif 

     if (id_surf_tracer(n) > 0) then
        if(n==index_temp ) then
           wrk1_2d(:,:) = 0.0
           do j=jsc,jec
              do i=isc,iec
                 wrk1_2d(i,j) = T_prog(n)%field(i,j,1,tau) + Grd%tmask(i,j,1)*cmip_offset 
              enddo
           enddo
           call diagnose_2d(Time, Grd, id_surf_tracer(n), wrk1_2d(:,:))
        else
           call diagnose_2d(Time, Grd, id_surf_tracer(n), T_prog(n)%field(:,:,1,tau))
        endif
     endif

     if (id_surf_tracer_sq(n) > 0) then
         if(n==index_temp) then 
             wrk1_2d(:,:) = 0.0
             do j=jsc,jec
                do i=isc,iec
                   wrk1_2d(i,j) = T_prog(n)%field(i,j,1,tau) + Grd%tmask(i,j,1)*cmip_offset 
                enddo
             enddo
             call diagnose_2d(Time, Grd, id_surf_tracer_sq(n), wrk1_2d(:,:)**2)
         else
            call diagnose_2d(Time, Grd, id_surf_tracer_sq(n), T_prog(n)%field(:,:,1,tau)**2)
         endif
     endif

     if (id_bott_tracer(n) > 0) then
         wrk1_2d(:,:) = 0.0
         do j = jsd,jed
            do i = isd,ied
               kbot = Grd%kmt(i,j)    
               if(kbot > 1) then 
                   wrk1_2d(i,j) = T_prog(n)%field(i,j,kbot,tau)
               endif
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_bott_tracer(n), wrk1_2d(:,:))
     endif

     ! time tendency for tracer mass per horizontal area 
     if (id_tendency(n) > 0 .or. id_tendency_on_nrho(n) > 0) then
         wrk1(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1(i,j,k) = dtimer*T_prog(n)%conversion                        &
                     *(T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                      -T_prog(n)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1))
                     
               enddo
            enddo
         enddo
         if (id_tendency(n) > 0) then
            call diagnose_3d(Time, Grd, id_tendency(n),wrk1(:,:,:))
         endif
         if (id_tendency_on_nrho(n) > 0) then
            call diagnose_3d_rho(Time, Dens, id_tendency_on_nrho(n),wrk1)
         endif
     endif

     ! time tendency for tracer concentration 
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(n)%tendency(i,j,k) = dtimer*(T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,k,taum1))
           enddo
        enddo
     enddo
     call diagnose_3d(Time, Grd, id_tendency_conc(n), T_prog(n)%tendency(:,:,:))


     ! send some diagnostics for the blob systems 
     if (use_blobs) then

        if (id_progT(n) > 0) then 
           if(n==index_temp) then
              wrk1(:,:,:) = 0.0
              do k=1,nk
                 do j=jsc,jec
                    do i=isc,iec
                       wrk1(i,j,k) = T_prog(n)%fieldT(i,j,k) + Grd%tmask(i,j,k)*cmip_offset 
                    enddo
                 enddo
              enddo
              call diagnose_3d(Time, Grd, id_progT(n), wrk1(:,:,:))
           else
              call diagnose_3d(Time, Grd, id_progT(n), T_prog(n)%fieldT(:,:,:))
           endif
        endif

        if(id_prog_rhodztT(n)  > 0) then
           call diagnose_3d(Time, Grd, id_prog_rhodzt(n), Thickness%rho_dztT(:,:,:,tau)*T_prog(n)%fieldT(:,:,:))
        endif

        ! vertical sum of tracer*rho_dzt 
        if(id_prog_int_rhodzT(n) > 0) then 
           wrk1_2d(:,:) = 0.0
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%rho_dztT(i,j,k,tau)*T_prog(n)%fieldT(i,j,k) 
                 enddo
              enddo
           enddo
           call diagnose_2d(Time, Grd, id_prog_int_rhodzT(n), wrk1_2d(:,:))
        endif

        if(id_surf_tracerT(n) > 0) then
           if(n==index_temp ) then
              wrk1_2d(:,:) = 0.0
              do j=jsc,jec
                 do i=isc,iec
                    wrk1_2d(i,j) = T_prog(n)%fieldT(i,j,1) + Grd%tmask(i,j,1)*cmip_offset 
                 enddo
              enddo
              call diagnose_2d(Time, Grd, id_surf_tracer_sq(n), wrk1_2d(:,:))
           else
              call diagnose_2d(Time, Grd, id_surf_tracer_sq(n), T_prog(n)%fieldT(:,:,1))
           endif
        endif

        if (id_bott_tracerT(n) > 0) then
           wrk1_2d(:,:) = 0.0
           do j = jsd,jed
              do i = isd,ied
                 kbot = Grd%kmt(i,j)    
                 if(kbot > 1) then 
                    wrk1_2d(i,j) = T_prog(n)%fieldT(i,j,kbot)*Grd%tmask(i,j,1)
                 endif
              enddo
           enddo
           call diagnose_2d(Time, Grd, id_bott_tracerT(n),wrk1_2d(:,:))
        endif
        
        if (id_tendencyT(n) > 0) then
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = dtimer*T_prog(n)%conversion                          &
                         *(T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                         -T_prog(n)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1)  &
                         +T_prog(n)%sum_blob(i,j,k,taup1)*Grd%datr(i,j)                &
                         -T_prog(n)%sum_blob(i,j,k,taum1)*Grd%datr(i,j) )
                    
                 enddo
              enddo
           enddo
           call diagnose_3d(Time, Grd, id_tendencyT(n), wrk1(:,:,:))
        endif

        if (id_tendencyL(n) > 0) then
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = dtimer*T_prog(n)%conversion            &
                         *(T_prog(n)%sum_blob(i,j,k,taup1)*Grd%datr(i,j) &
                         -T_prog(n)%sum_blob(i,j,k,taum1)*Grd%datr(i,j) )
                     
                 enddo
              enddo
           enddo
           call diagnose_3d(Time, Grd, id_tendencyL(n),wrk1(:,:,:))
        endif
        
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 wrk1(i,j,k) = ( T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                      +Grd%datr(i,j)*T_prog(n)%sum_blob(i,j,k,taup1) )                       &
                      /(Thickness%rho_dztT(i,j,k,taup1)+epsln)
                 wrk2(i,j,k) = ( T_prog(n)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1) &
                      +Grd%datr(i,j)*T_prog(n)%sum_blob(i,j,k,taum1) )                       &
                      /(Thickness%rho_dztT(i,j,k,taum1)+epsln)
                 T_prog(n)%tendency(i,j,k) = dtimer*(wrk1(i,j,k) - wrk2(i,j,k))
              enddo
           enddo
        enddo
       
        call diagnose_3d(Time, Grd, id_tendency_concT(n), T_prog(n)%tendency(:,:,:))
        
        if (id_tendency_concL(n) > 0) then
           do k=1,nk
              wrk1(:,:,k) = Grd%datr(:,:)*T_prog(n)%sum_blob(:,:,k,taup1)/Thickness%rho_dztL(:,:,k,taup1)
              wrk2(:,:,k) = Grd%datr(:,:)*T_prog(n)%sum_blob(:,:,k,taum1)/Thickness%rho_dztL(:,:,k,taum1)
           enddo
           call diagnose_3d(Time, Grd, id_tendency_concL(n),        &
                dtimer*(wrk1(:,:,:) - wrk2(:,:,:)))
        endif

     endif  !use_blobs

     call send_tracer_variance(Time, T_prog(n), Thickness, n)

  enddo   ! enddo for n-tracers 

  ! send diagnostic tracers to diag manager 
  do n=1,num_diag_tracers

     if(id_diag(n) > 0) then 
         if(n==index_temp) then 
             wrk1(:,:,:) = 0.0
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      wrk1(i,j,k) = T_diag(n)%field(i,j,k) + Grd%tmask(i,j,k)*cmip_offset 
                   enddo
                enddo
             enddo
             call diagnose_3d(Time, Grd, id_diag(n), wrk1(:,:,:))
         else
            call diagnose_3d(Time, Grd, id_diag(n), T_diag(n)%field(:,:,:))
         endif
     endif

     if (id_diag_surf_tracer(n) > 0) then 
         if(n==index_temp) then 
             wrk1_2d(:,:) = 0.0
             do j=jsc,jec
                do i=isc,iec
                   wrk1_2d(i,j) = T_diag(n)%field(i,j,1) + Grd%tmask(i,j,1)*cmip_offset 
                enddo
             enddo
             call diagnose_2d(Time, Grd, id_diag_surf_tracer(n), wrk1_2d(:,:))
         else
            call diagnose_2d(Time, Grd, id_diag_surf_tracer(n), T_diag(n)%field(:,:,1))
         endif

     endif

     if (id_diag_surf_tracer_sq(n) > 0) then
         if(n==index_temp) then 
             wrk1_2d(:,:) = 0.0
             do j=jsc,jec
                do i=isc,iec
                   wrk1_2d(i,j) = T_diag(n)%field(i,j,1) + Grd%tmask(i,j,1)*cmip_offset 
                enddo
             enddo
             call diagnose_2d(Time, Grd, id_diag_surf_tracer_sq(n), wrk1_2d(:,:)**2)
         else
            call diagnose_2d(Time, Grd, id_diag_surf_tracer_sq(n), T_diag(n)%field(:,:,1)**2)
         endif
     endif

     if(id_diag_total(n) > 0) then 
         total_tracer=0.0
         do k=1,nk
            wrk1_2d(:,:) = T_diag(n)%conversion*Grd%tmask(:,:,k)*Grd%dat(:,:) &
                          *Thickness%rho_dzt(:,:,k,tau)*T_diag(n)%field(:,:,k)
            if(have_obc) wrk1_2d(:,:) = wrk1_2d(:,:)*Grd%obc_tmask(:,:)
            total_tracer = total_tracer + mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
         enddo
         used = send_data (id_diag_total(n), total_tracer*1e-18, Time%model_time)
     endif

  enddo ! enddo for n=1,num_diag_tracers


end subroutine send_tracer_diagnostics 
! </SUBROUTINE>  NAME="send_tracer_diagnostics"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from interior ocean vertical mixing of 
! temp and salt on the watermass transformation diagnostics.  
!
! This routine is called prior to implicit update of the tracer
! fields, so that taup1 value contains only explicit in-time 
! tendencies.  The implicit in-time tendencies are diagnosed
! in this routine by various calls to invtri using same methods 
! as for the prognostic calculation. 
! 
! </DESCRIPTION>
!
subroutine watermass_diag(Time, T_prog, T_diag, Dens, Thickness, pme)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(in) :: T_diag(:)
  type(ocean_density_type),       intent(in) :: Dens
  type(ocean_thickness_type),     intent(in) :: Thickness
  real,  dimension(isd:,jsd:),    intent(in) :: pme

  integer :: i,j,k,taum1,tau,taup1

  real,  dimension(isd:ied,jsd:jed) :: eta_tend

  if(.not. compute_watermass_diag) return 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1


  ! full time tendency 
  if (id_neut_rho_tendency  > 0 .or. id_neut_rho_tendency_on_nrho  > 0 .or.  &
      id_wdian_rho_tendency > 0 .or. id_wdian_rho_tendency_on_nrho > 0 .or.  &
      id_tform_rho_tendency > 0 .or. id_tform_rho_tendency_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0
      wrk6(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Dens%drhodT(i,j,k)*                                          &
                    (T_prog(index_temp)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                    -T_prog(index_temp)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1))
               wrk2(i,j,k) = Dens%drhodS(i,j,k)*                                          &
                    (T_prog(index_salt)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                    -T_prog(index_salt)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1))    
               wrk3(i,j,k) = Grd%tmask(i,j,k)*dtimer*(wrk1(i,j,k)+wrk2(i,j,k))
               wrk4(i,j,k) = wrk3(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk5(i,j,k) = wrk4(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk6(i,j,k) = wrk3(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_tendency, wrk4(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_tendency, wrk5(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_tendency, wrk6(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_tendency_on_nrho, wrk4)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_tendency_on_nrho, wrk5)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_tendency_on_nrho, wrk6)

  endif

  ! time tendency arising from temperature effects 
  if (id_neut_temp_tendency  > 0 .or. id_neut_temp_tendency_on_nrho  > 0 .or.  &
      id_wdian_temp_tendency > 0 .or. id_wdian_temp_tendency_on_nrho > 0 .or.  &
      id_tform_temp_tendency > 0 .or. id_tform_temp_tendency_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0
      wrk6(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Dens%drhodT(i,j,k)*                                          &
                    (T_prog(index_temp)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                    -T_prog(index_temp)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1))
               wrk2(i,j,k) = 0.0
               wrk3(i,j,k) = Grd%tmask(i,j,k)*dtimer*(wrk1(i,j,k)+wrk2(i,j,k))
               wrk4(i,j,k) = wrk3(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk5(i,j,k) = wrk4(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk6(i,j,k) = wrk3(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_temp_tendency,wrk4(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_temp_tendency, wrk5(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_temp_tendency, wrk6(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_temp_tendency_on_nrho, wrk4)
      call diagnose_3d_rho(Time, Dens, id_wdian_temp_tendency_on_nrho, wrk5)
      call diagnose_3d_rho(Time, Dens, id_tform_temp_tendency_on_nrho, wrk6)
  endif


  ! time tendency arising from salinity effects 
  if (id_neut_salt_tendency  > 0 .or. id_neut_salt_tendency_on_nrho  > 0 .or.  &
      id_wdian_salt_tendency > 0 .or. id_wdian_salt_tendency_on_nrho > 0 .or.  &
      id_tform_salt_tendency > 0 .or. id_tform_salt_tendency_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0
      wrk6(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = 0.0
               wrk2(i,j,k) = Dens%drhodS(i,j,k)*                                          &
                    (T_prog(index_salt)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                    -T_prog(index_salt)%field(i,j,k,taum1)*Thickness%rho_dzt(i,j,k,taum1))    
               wrk3(i,j,k) = Grd%tmask(i,j,k)*dtimer*(wrk1(i,j,k)+wrk2(i,j,k))
               wrk4(i,j,k) = wrk3(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk5(i,j,k) = wrk4(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk6(i,j,k) = wrk3(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_salt_tendency, wrk4(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_salt_tendency, wrk5(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_salt_tendency, wrk6(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_salt_tendency_on_nrho, wrk4)
      call diagnose_3d_rho(Time, Dens, id_wdian_salt_tendency_on_nrho, wrk5)
      call diagnose_3d_rho(Time, Dens, id_tform_salt_tendency_on_nrho, wrk6)
  endif



  ! effects from smooth eta or pbot 
  if (id_neut_rho_smooth  > 0 .or. id_neut_rho_smooth_on_nrho  > 0 .or.  &
      id_wdian_rho_smooth > 0 .or. id_wdian_rho_smooth_on_nrho > 0 .or.  &
      id_tform_rho_smooth > 0 .or. id_tform_rho_smooth_on_nrho > 0 .or.  &
      id_eta_tend_smooth  > 0 .or. id_eta_tend_smooth_glob     > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0
      wrk6(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)                                &
                 *( Dens%drhodT(i,j,k)*T_prog(index_temp)%eta_smooth(i,j) &
                   +Dens%drhodS(i,j,k)*T_prog(index_salt)%eta_smooth(i,j))  
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            k=Grd%kmt(i,j)
            if(k > 0) then
                wrk2(i,j,k) = Grd%tmask(i,j,k)                                 &
                     *( Dens%drhodT(i,j,k)*T_prog(index_temp)%pbot_smooth(i,j) &
                       +Dens%drhodS(i,j,k)*T_prog(index_salt)%pbot_smooth(i,j))  
            endif
         enddo
      enddo
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk3(i,j,k) = Grd%tmask(i,j,k)*(wrk1(i,j,k)+wrk2(i,j,k))
               wrk4(i,j,k) = wrk3(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk5(i,j,k) = wrk4(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk6(i,j,k) = wrk3(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
               wrk1(i,j,k) =-wrk3(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend 
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_smooth, wrk4(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_smooth, wrk5(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_smooth, wrk6(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_smooth_on_nrho, wrk4)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_smooth_on_nrho, wrk5)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_smooth_on_nrho, wrk6)
      if(id_eta_tend_smooth > 0 .or. id_eta_tend_smooth_glob > 0) then
          eta_tend(:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   eta_tend(i,j) = eta_tend(i,j) + wrk1(i,j,k)
                enddo
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_smooth, eta_tend(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_smooth_glob, eta_tend, cellarea_r)
      endif

  endif


  ! effects from smooth eta or pbot due to temperature contributions 
  if (id_neut_temp_smooth  > 0 .or. id_neut_temp_smooth_on_nrho  > 0 .or.  &
      id_wdian_temp_smooth > 0 .or. id_wdian_temp_smooth_on_nrho > 0 .or.  &
      id_tform_temp_smooth > 0 .or. id_tform_temp_smooth_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0
      wrk6(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)                                &
                 *( Dens%drhodT(i,j,k)*T_prog(index_temp)%eta_smooth(i,j) &
                   +0.0)  
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            k=Grd%kmt(i,j)
            if(k > 0) then
                wrk2(i,j,k) = Grd%tmask(i,j,k)                                 &
                     *( Dens%drhodT(i,j,k)*T_prog(index_temp)%pbot_smooth(i,j) &
                       +0.0)  
            endif
         enddo
      enddo
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk3(i,j,k) = Grd%tmask(i,j,k)*(wrk1(i,j,k)+wrk2(i,j,k))
               wrk4(i,j,k) = wrk3(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk5(i,j,k) = wrk4(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk6(i,j,k) = wrk3(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_temp_smooth, wrk4(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_temp_smooth, wrk5(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_temp_smooth, wrk6(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_temp_smooth_on_nrho, wrk4)
      call diagnose_3d_rho(Time, Dens, id_wdian_temp_smooth_on_nrho, wrk5)
      call diagnose_3d_rho(Time, Dens, id_tform_temp_smooth_on_nrho, wrk6)
  endif


  ! effects from smooth eta or pbot due to salinty contributions 
  if (id_neut_salt_smooth  > 0 .or. id_neut_salt_smooth_on_nrho  > 0 .or.  &
      id_wdian_salt_smooth > 0 .or. id_wdian_salt_smooth_on_nrho > 0 .or.  &
      id_tform_salt_smooth > 0 .or. id_tform_salt_smooth_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      wrk5(:,:,:) = 0.0
      wrk6(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)                                &
                 *( 0.0                                                   &
                   +Dens%drhodS(i,j,k)*T_prog(index_salt)%eta_smooth(i,j))  
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            k=Grd%kmt(i,j)
            if(k > 0) then
                wrk2(i,j,k) = Grd%tmask(i,j,k)                                 &
                     *( 0.0                                                    &
                       +Dens%drhodS(i,j,k)*T_prog(index_salt)%pbot_smooth(i,j))  
            endif
         enddo
      enddo
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk3(i,j,k) = Grd%tmask(i,j,k)*(wrk1(i,j,k)+wrk2(i,j,k))
               wrk4(i,j,k) = wrk3(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk5(i,j,k) = wrk4(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk6(i,j,k) = wrk3(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_salt_smooth,wrk4(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_salt_smooth, wrk5(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_salt_smooth, wrk6(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_salt_smooth_on_nrho, wrk4)
      call diagnose_3d_rho(Time, Dens, id_wdian_salt_smooth_on_nrho, wrk5)
      call diagnose_3d_rho(Time, Dens, id_tform_salt_smooth_on_nrho, wrk6)
  endif


  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme. formulated here according to flux form of material time derivative.
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  if (id_neut_rho_pme  > 0 .or. id_neut_rho_pme_on_nrho  > 0 .or.  &
      id_wdian_rho_pme > 0 .or. id_wdian_rho_pme_on_nrho > 0 .or.  &
      id_tform_rho_pme > 0 .or. id_tform_rho_pme_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                  &
                 *( Dens%drhodT(i,j,k)*T_prog(index_temp)%tpme(i,j)  &
                   +Dens%drhodS(i,j,k)*T_prog(index_salt)%tpme(i,j))  
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_pme, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_pme, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_pme, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_pme_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_pme_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pme_on_nrho, wrk4)
  endif


  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme induced buoyancy flux. written here according to advective form 
  ! of material time derivative using the kinematic formulation. 
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  if (id_neut_rho_pbl_pme_kn  > 0 .or. id_neut_rho_pbl_pme_kn_on_nrho  > 0 .or.  &
      id_wdian_rho_pbl_pme_kn > 0 .or. id_wdian_rho_pbl_pme_kn_on_nrho > 0 .or.  &
      id_tform_rho_pbl_pme_kn > 0 .or. id_tform_rho_pbl_pme_kn_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                             &
                 *(Dens%drhodT(i,j,k)*(0.0-T_prog(index_temp)%field(i,j,k,tau)) &
                  +Dens%drhodS(i,j,k)*(0.0-T_prog(index_salt)%field(i,j,k,tau)))
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_pbl_pme_kn, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_pbl_pme_kn, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_pbl_pme_kn, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_pbl_pme_kn_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_pbl_pme_kn_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_pme_kn_on_nrho, wrk4)
  endif 
  

  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme induced heat flux. written here according to advective form 
  ! of material time derivative using the kinematic formulation. 
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  if (id_neut_temp_pbl_pme_kn > 0 .or. id_neut_temp_pbl_pme_kn_on_nrho > 0 .or.  &
      idwdian_temp_pbl_pme_kn > 0 .or. idwdian_temp_pbl_pme_kn_on_nrho > 0 .or.  &
      idtform_temp_pbl_pme_kn > 0 .or. idtform_temp_pbl_pme_kn_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                             &
                 *(Dens%drhodT(i,j,k)*(0.0-T_prog(index_temp)%field(i,j,k,tau)) &
                  +0.0)
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_temp_pbl_pme_kn, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, idwdian_temp_pbl_pme_kn, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, idtform_temp_pbl_pme_kn, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_temp_pbl_pme_kn_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, idwdian_temp_pbl_pme_kn_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, idtform_temp_pbl_pme_kn_on_nrho, wrk4)
  endif 
  

  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme induced salt flux. written here according to advective form 
  ! of material time derivative using the kinematic formulation. 
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  if (id_neut_salt_pbl_pme_kn > 0 .or. id_neut_salt_pbl_pme_kn_on_nrho > 0 .or.  &
      idwdian_salt_pbl_pme_kn > 0 .or. idwdian_salt_pbl_pme_kn_on_nrho > 0 .or.  &
      idtform_salt_pbl_pme_kn > 0 .or. idtform_salt_pbl_pme_kn_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                             &
                 *(0.0                                                          &
                  +Dens%drhodS(i,j,k)*(0.0-T_prog(index_salt)%field(i,j,k,tau)))
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_salt_pbl_pme_kn, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, idwdian_salt_pbl_pme_kn, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, idtform_salt_pbl_pme_kn, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_salt_pbl_pme_kn_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, idwdian_salt_pbl_pme_kn_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, idtform_salt_pbl_pme_kn_on_nrho, wrk4)
  endif 
  


  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme induced buoyancy flux. written here according to advective form 
  ! of material time derivative using the process formulation. 
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  if (id_neut_rho_pbl_pme_pr  > 0 .or. id_neut_rho_pbl_pme_pr_on_nrho  > 0 .or.  &
      id_wdian_rho_pbl_pme_pr > 0 .or. id_wdian_rho_pbl_pme_pr_on_nrho > 0 .or.  &
      id_tform_rho_pbl_pme_pr > 0 .or. id_tform_rho_pbl_pme_pr_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                                                      &
                 *(Dens%drhodT(i,j,k)*(T_prog(index_temp)%tpme(i,j)-T_prog(index_temp)%field(i,j,k,tau)) &
                  +Dens%drhodS(i,j,k)*(T_prog(index_salt)%tpme(i,j)-T_prog(index_salt)%field(i,j,k,tau)))
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_pbl_pme_pr, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_pbl_pme_pr, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_pbl_pme_pr, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_pbl_pme_pr_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_pbl_pme_pr_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_pme_pr_on_nrho, wrk4)
  endif 
  

  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme induced buoyancy flux. written here according to advective form 
  ! of material time derivative using the process formulation. 
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  ! temperature contributions are diagnosed here. 
  if (id_neut_temp_pbl_pme_pr > 0 .or. id_neut_temp_pbl_pme_pr_on_nrho > 0 .or.  &
      idwdian_temp_pbl_pme_pr > 0 .or. idwdian_temp_pbl_pme_pr_on_nrho > 0 .or.  &
      idtform_temp_pbl_pme_pr > 0 .or. idtform_temp_pbl_pme_pr_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                                                      &
                 *(Dens%drhodT(i,j,k)*(T_prog(index_temp)%tpme(i,j)-T_prog(index_temp)%field(i,j,k,tau)) &
                  +0.0)
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_temp_pbl_pme_pr, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, idwdian_temp_pbl_pme_pr, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, idtform_temp_pbl_pme_pr, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_temp_pbl_pme_pr_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, idwdian_temp_pbl_pme_pr_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, idtform_temp_pbl_pme_pr_on_nrho, wrk4)
  endif 
  

  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to pme induced buoyancy flux. written here according to advective form 
  ! of material time derivative using the process formulation. 
  ! note that river contribution is diagnosed in ocean_rivermix.F90.
  ! salinity contributions are diagnosed here. 
  if (id_neut_salt_pbl_pme_pr > 0 .or. id_neut_salt_pbl_pme_pr_on_nrho > 0 .or.  &
      idwdian_salt_pbl_pme_pr > 0 .or. idwdian_salt_pbl_pme_pr_on_nrho > 0 .or.  &
      idtform_salt_pbl_pme_pr > 0 .or. idtform_salt_pbl_pme_pr_on_nrho > 0) then

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%tmask(i,j,k)*pme(i,j)                                                      &
                 *(0.0                                                                                   &
                  +Dens%drhodS(i,j,k)*(T_prog(index_salt)%tpme(i,j)-T_prog(index_salt)%field(i,j,k,tau)))
            wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
            wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
            wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_salt_pbl_pme_pr, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, idwdian_salt_pbl_pme_pr, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, idtform_salt_pbl_pme_pr, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_salt_pbl_pme_pr_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, idwdian_salt_pbl_pme_pr_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, idtform_salt_pbl_pme_pr_on_nrho, wrk4)
  endif 
  


  ! time tendency for locally ref potrho and dianeutral velocity component due
  ! to frazil formation, which is a heating of the liquid ocean. 
  if(index_frazil > 0) then 

      if (id_neut_rho_frazil  > 0 .or. id_neut_rho_frazil_on_nrho  > 0 .or.  &
          id_wdian_rho_frazil > 0 .or. id_wdian_rho_frazil_on_nrho > 0 .or.  &
          id_tform_rho_frazil > 0 .or. id_tform_rho_frazil_on_nrho > 0 .or.  &
          id_eta_tend_frazil  > 0 .or. id_eta_tend_frazil_glob     > 0 .or.  &
          id_frazil_on_nrho   > 0) then

          wrk1(:,:,:) = 0.0
          wrk2(:,:,:) = 0.0
          wrk3(:,:,:) = 0.0
          wrk4(:,:,:) = 0.0
          wrk5(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = Grd%tmask(i,j,k)*dtimer*T_diag(index_frazil)%field(i,j,k)*cp_oceanr*Dens%drhodT(i,j,k)
                   wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                   wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
                   wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
                   wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) ! for eta_tend 
                enddo
             enddo
          enddo

          call diagnose_3d(Time, Grd, id_neut_rho_frazil, wrk2(:,:,:))
          call diagnose_3d(Time, Grd, id_wdian_rho_frazil, wrk3(:,:,:))
          call diagnose_3d(Time, Grd, id_tform_rho_frazil, wrk4(:,:,:))
          call diagnose_3d_rho(Time, Dens, id_neut_rho_frazil_on_nrho, wrk2)
          call diagnose_3d_rho(Time, Dens, id_wdian_rho_frazil_on_nrho, wrk3)
          call diagnose_3d_rho(Time, Dens, id_tform_rho_frazil_on_nrho, wrk4)
          call diagnose_3d_rho(Time, Dens, id_frazil_on_nrho, T_diag(index_frazil)%field*dtimer)
          if(id_eta_tend_frazil > 0 .or. id_eta_tend_frazil_glob > 0) then
              eta_tend(:,:) = 0.0
              do k=1,nk
                 do j=jsc,jec
                    do i=isc,iec
                       eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
                    enddo
                 enddo
              enddo
              call diagnose_2d(Time, Grd, id_eta_tend_frazil, eta_tend(:,:))
              call diagnose_sum(Time, Grd, Dom, id_eta_tend_frazil_glob, eta_tend, cellarea_r)
          endif
      endif 

  endif ! endif for index_frazil > 0
  

  ! some further layer diagnostics for mass changes.
  ! these diagnostics do not lead to robust results as a function of 
  ! layer binning. User should thus beware. 
  if(id_mass_t_on_nrho > 0 .or. id_mass_t_tendency_on_nrho > 0) then

      ! get mass of layers at time taup1 via binning mass_t. 
      nrho_work(:,:,:) = 0.0
      wrk1(:,:,:)      = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%dat(i,j)*Thickness%rho_dzt(i,j,k,taup1) 
            enddo
         enddo
      enddo
      wrk3(:,:,:) = Grd%tmask(:,:,:)*neutral_density(T_prog(index_salt)%field(:,:,:,taup1), &
                                                     T_prog(index_temp)%field(:,:,:,taup1))
      call rebin_onto_rho (Dens%neutralrho_bounds, wrk3, wrk1, nrho_work) 

      ! we diagnose mass_t_on_nrho at taup1 in order to be 
      ! consistent with mass_nrho_layer diagnosed below. 
      ! The mass computed with mass_nrho_layer method can 
      ! only get the mass(taup1) since that is the time 
      ! where the level thicknesses are known.  
      if(id_mass_t_on_nrho > 0 ) then 
          used = send_data (id_mass_t_on_nrho, nrho_work(:,:,:), &
               Time%model_time,                                  &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif

      ! bin mass_t at two adjacent time steps to estimate the layer tendency 
      if(id_mass_t_tendency_on_nrho > 0) then
          nrho_work2(:,:,:) = 0.0
          wrk1(:,:,:)       = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = Grd%dat(i,j)*Thickness%rho_dzt(i,j,k,tau)
                enddo
             enddo
          enddo

          ! get mass of layers at time taup1 via binning mass_t. 
          ! we already have the mass at time taup1 contained in nrho_work.
          ! note that if mass_t(tau)=mass_t(taup1), then only means for getting a 
          ! layer tendency is if the layer interfaces evolve.  Indeed, such is the 
          ! dominant effect, since interior grid cells change their mass only a 
          ! tiny amount, since interior ocean grid cells are nearly 3d non-divergent. 
          call rebin_onto_rho (Dens%neutralrho_bounds, Dens%neutralrho, wrk1, nrho_work2) 

          ! compute time tendency for layer mass.  
          do k=1,neutralrho_nk
             do j=jsc,jec
                do i=isc,iec
                   nrho_work(i,j,k) = (nrho_work(i,j,k) - nrho_work2(i,j,k))*dtimer
                enddo
             enddo
          enddo
          used = send_data (id_mass_t_tendency_on_nrho, nrho_work(:,:,:),&
               Time%model_time,                                          &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif

  endif ! endif corresponding to  if(id_mass_t_on_nrho > 0 .or. id_mass_t_tendency_on_nrho > 0) 


  ! estimate layer mass by estimating mass per area of the layer interfaces.
  ! note that we do not save nrho_mass in a restart file.  So there will be 
  ! a slight blip at the start of a new model segment.  This blip is not 
  ! so important, since we are only working with diagnostics here.  
  if(id_mass_nrho_layer > 0 .or. id_mass_nrho_tendency_layer > 0) then

      nrho_work(:,:,:) = nrho_mass(:,:,:)   ! save previous time step layer mass
      nrho_mass(:,:,:) = 0.0                ! for the updated layer mass 
      wrk3(:,:,:)      = Grd%tmask(:,:,:)*neutral_density(T_prog(index_salt)%field(:,:,:,taup1), &
                                                          T_prog(index_temp)%field(:,:,:,taup1))

      ! note that all thickness arrays are at taup1, so we 
      ! compute mass(taup1) and mass_of_layer(taup1).  
      call diagnose_mass_of_layer (Grd%dat, &
                                   Thickness%dzt, Thickness%dztlo, Thickness%dztup,  &
                                   Thickness%rho_dzt(:,:,:,taup1),                   &
                                   wrk3, neutralrho_nk, Dens%neutralrho_bounds, nrho_mass)    

      if(id_mass_nrho_layer > 0) then 
         used = send_data (id_mass_nrho_layer, nrho_mass(:,:,:), &
              Time%model_time,                                   &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif 

      if(id_mass_nrho_tendency_layer > 0) then

          ! first time step need to compute mass_of_layer using the binning 
          ! approach, since we do not have access to the tau values of all 
          ! relevant Thickness arrays.  
          if(Time%itt <= 1) then 
              nrho_work(:,:,:) = 0.0
              wrk1(:,:,:)      = 0.0
              do k=1,nk
                 do j=jsc,jec
                    do i=isc,iec
                       wrk1(i,j,k) = Grd%dat(i,j)*Thickness%rho_dzt(i,j,k,tau)
                    enddo
                 enddo
              enddo
              call rebin_onto_rho (Dens%neutralrho_bounds, Dens%neutralrho, wrk1, nrho_work) 
          endif 
          do k=1,neutralrho_nk
             do j=jsc,jec
                do i=isc,iec
                   nrho_work(i,j,k) = (nrho_mass(i,j,k) - nrho_work(i,j,k))*dtimer
                enddo
             enddo
          enddo
          used = send_data (id_mass_nrho_tendency_layer, nrho_work(:,:,:),&
               Time%model_time,                                           &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif

  endif

end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


    
end module ocean_tracer_mod
