module ocean_blob_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This is the main "driver" module for the Lagrangian blob scheme.  This
! module calls other modules that contain the individual parameterisations.
!
! Please note that the preprocessor option MOM_STATIC_ARRAYS is NOT
! supported.  This is because to run a model where the memeory statically
! allocated without the blob framework will incur a large, unnecessary
! increase in memory requirements.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module is the top-level module of the Lagrangian blob framework.
! It calls routines to form blobs, to integrate their properties, to
! transfer them from one dynamic regime to another, and also calculates
! the L system contribution towards grid cell thickness.
!
! This module also handles some framework wide variables, mostly
! associated with system wide diagnostics and the system wide accounting
! required to ensure the Eulerian model and the Lagrangian model can coexist.
!
! It should be noted that many of the parameterisations are not mutually
! exclusive.  As such, care should be excercised when creating namelists
! for experiments.
!</DESCRIPTION>
!
!<INFO>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM4p1 (2009)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
!</INFO>
!<NAMELIST NAME="ocean_blob_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module.
!  Default is use_this_module=.false.
!  </DATA>
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  Writes additional diagnostic data to fms.out.  This
!  also controls debug output for the other related blob
!  modules.
!  Default is debug_this_module=.false.
!  </DATA>
!
!  <DATA NAME="really_debug" TYPE="logical">
!  Be careful what you wish for, this outputs A LOT of
!  diagnostics to standard out!
!  Default is debug_this_module=.false.
!  </DATA>
!
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  When global sum outputs are done there is additional
!  computational expense to ensure that they are bitwise
!  the same across an arbitrary number of processors.
!  However, for debugging purposes, it can be useful
!  for global sums to be the same.  Note, that this differs
!  from bitwise_reproduction in that it do_bitwise_exact_sum
!  only applies to the mpp_global_sum diagnostic.
!  Note that this flag controls the output for all associated
!  blob modules.
!  Default is do_bitwise_exact_sum=.false.
!  </DATA>
!
!  <DATA NAME="bitwise_reproduction" TYPE="logical">
!  There is additional cost involved in ensuring that
!  results are reproducable across an arbitrary number of
!  processors and across restarts.
!  Bitwise reproduction is a very memory intensive operation and should
!  only be used for debugging.  For bitwise_reproduction=.true. We need
!  to process blobs and their histories in the same relative order
!  regardless of domain decomposition and restarts. To do so, we save the
!  "history" of each blob subcycle is saved to a number of arrays (which
!  can be a very memory intensive process) and process them in order.
!  Note that this flag controls reproducability for all associated
!  blob modules.
!  Bitwise reproducibility is only possible with the appropriate
!  compiler flags AND when the simulation is run on hardware that is
!  capable of producing bitwise reproduction.
!  Default is bitwise_reproduction=.false.
!  </DATA>
!
!  <DATA NAME="blob_small_mass" UNITS="kg" TYPE="real">
!  Will delete blobs of mass less than blob_small_mass.
!  Note that this variable is for all associated blob
!  modules.  The deletion of blobs is a conservative
!  action, any mass/tracer fields that are nonzero
!  have the remaining properties transferred back to the
!  Eulerian system.  So, in principle, blob_small_mass
!  can actually be a relatively large number, and the
!  model will remain conservative.  It has been found in
!  certain test cases (with very low tracer values) that
!  setting blob_small_mass to be very small (i.e. <1e2)
!  that roundoff error can cause non-trivial errors.  So,
!  it is recommended that blob_small_mass be no smaller than
!  than 1e3 kg (which is approximately 1.0m**3 -- a very small
!  blob!)
!  Default is blob_small_mass=1.e3
!  </DATA>
!
!  <DATA NAME="mass_prop_thickness" UNITS="dimensionless" TYPE="real">
!  Sets the maximum proportion of a grid cell that the
!  Lagrangian system may occupy.  This is actually calculated
!  separately (and therefore must be satisfied separately) for
!  the uppper and lower portions of a grid cell.
!  Default is blob_small_mass=0.7
!  </DATA>
!
!</NAMELIST>
!
use constants_mod,   only: epsln
use diag_manager_mod,only: register_diag_field, send_data
use fms_mod,         only: open_namelist_file, check_nml_error, file_exist
use fms_mod,         only: write_version_number, close_file, stderr
use fms_mod,         only: mpp_error, FATAL, WARNING, stdout, stdlog
use fms_io_mod,      only: register_restart_field, save_restart
use fms_io_mod,      only: restore_state, restart_file_type
use mpp_domains_mod, only: mpp_update_domains, BGRID_NE
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use mpp_domains_mod, only: NORTH_EAST, NORTH_WEST, SOUTH_EAST, SOUTH_WEST
use mpp_domains_mod, only: BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,         only: mpp_pe, mpp_sum, CLOCK_ROUTINE, NULL_PE
use mpp_mod,         only: mpp_clock_id, mpp_clock_begin, mpp_clock_end

use ocean_blob_diag_mod,           only: blob_diag_init, blob_diag, blob_diag_end
use ocean_blob_dynamic_bottom_mod, only: blob_dynamic_bottom_init, blob_dynamic_bottom_update
use ocean_blob_dynamic_bottom_mod, only: transfer_free_to_bottom, dynamic_bottom_form_new
use ocean_blob_dynamic_bottom_mod, only: blob_dynamic_bottom_end
use ocean_blob_dynamic_free_mod,   only: blob_dynamic_free_init, blob_dynamic_free_implicit
use ocean_blob_dynamic_free_mod,   only: blob_dynamic_free_update, transfer_bottom_to_free
use ocean_blob_dynamic_free_mod,   only: blob_dynamic_free_end
use ocean_blob_static_free_mod,    only: blob_static_free_init, blob_static_free
use ocean_blob_static_free_mod,    only: blob_static_free_end
use ocean_blob_static_bottom_mod,  only: blob_static_bottom_init, blob_overflow_like
use ocean_blob_static_bottom_mod,  only: blob_static_bottom_end
use ocean_blob_util_mod,           only: inq_var, get_int, get_double, def_var
use ocean_blob_util_mod,           only: blob_util_init, E_and_L_totals, blob_util_end
use ocean_blob_util_mod,           only: blob_chksum, blob_delete, put_att, put_double
use ocean_blob_util_mod,           only: put_int, insert_blob, lagrangian_system_chksum
use ocean_blob_util_mod,           only: write_blobs, kill_blob, free_blob_memory
use ocean_blob_util_mod,           only: allocate_interaction_memory
use ocean_density_mod,             only: density, compute_buoyfreq
use ocean_domains_mod,             only: set_ocean_domain
use ocean_parameters_mod,          only: grav, rho0, rho0r, DEPTH_BASED, GEOPOTENTIAL, ZSIGMA
use ocean_parameters_mod,          only: PRESSURE, PSIGMA, ZSTAR, PSTAR
use ocean_types_mod,               only: ocean_time_type, ocean_grid_type, blob_diag_type
use ocean_types_mod,               only: ocean_domain_type, ocean_options_type
use ocean_types_mod,               only: ocean_density_type, ocean_lagrangian_type
use ocean_types_mod,               only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,               only: ocean_blob_type, blob_grid_type, ocean_adv_vel_type
use ocean_types_mod,               only: ocean_external_mode_type, ocean_velocity_type
use ocean_util_mod,                only: write_timestamp, diagnose_3d, write_chksum_2d, write_chksum_3d, write_chksum_3d_int
use ocean_workspace_mod,           only: wrk1

implicit none

include 'netcdf.inc'

private

! set up the various module specific types

! set module types
type(ocean_grid_type),   pointer :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_blob_type),   pointer :: head_dynamic_free => NULL()
type(ocean_blob_type),   pointer :: head_dynamic_bott => NULL()
type(ocean_blob_type),   pointer :: head_static_free  => NULL()
type(ocean_blob_type),   pointer :: head_static_bott  => NULL()
type(blob_grid_type),    pointer :: Info

! The special blob domain type for the extended halo
type(ocean_domain_type), save :: Blob_domain

! general module variables
integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1
integer :: vert_coordinate       ! tye type of vertical coordinate
integer :: vert_coordinate_class ! the class of vertical coordinate
integer :: global_sum_flag
integer :: isd,ied,jsd,jed
integer :: isc,iec,jsc,jec
integer :: isg,ieg,jsg,jeg
integer :: isbd, iebd, jsbd, jebd
integer :: nk
logical :: blob_diagnostics, use_dyn_fre, use_dyn_bot
real    :: free_minstep, bott_minstep
integer :: free_total_ns, bott_total_ns
integer, parameter :: file_format_major_version=0
integer, parameter :: file_format_minor_version=1
integer, allocatable, dimension(:,:,:) :: blob_counter
real,    allocatable, dimension(:,:)   :: datdtimer

!for restart
type(restart_file_type), save :: gridded_restart

! CVS stuff
character(len=128)  :: version = '$$'
character (len=128) :: tagname = '$Name: tikal $'

! initialisation
logical :: module_is_initialized=.false.
logical :: use_this_module

! identification numbers for mpp clocks
integer :: id_clock_static_bottom
integer :: id_clock_static_free
integer :: id_clock_dyn_free_implicit
integer :: id_clock_dyn_free_update
integer :: id_clock_dyn_bott_update
integer :: id_clock_tsfr_f2b
integer :: id_clock_dyn_bott_new

! gridded diagnostic fields
real, allocatable, dimension(:,:,:) :: mass_in
real, allocatable, dimension(:,:,:) :: mass_out

! identification for gridded diagnostic fields
integer, allocatable, dimension(:) :: id_tend_blob
integer, allocatable, dimension(:) :: id_entrainment
integer, allocatable, dimension(:) :: id_detrainment
integer, allocatable, dimension(:) :: id_new
integer, allocatable, dimension(:) :: id_dstry
integer :: id_prop_cell_mass
integer :: id_mass_in
integer :: id_mass_out
integer :: id_dstryd
integer :: id_grounded
integer :: id_surfaced
integer :: id_detraind
integer :: id_nblob_free
integer :: id_nblob_bott

logical :: used

! set the public subroutines
public ocean_blob_init
public ocean_blob_update
public update_L_thickness
public calculate_rhoT
public ocean_blob_implicit
public ocean_blob_diagnose_depth
public ocean_blob_end
public init_blob_thickness
public adjust_L_thickness
public write_all_blobs
public ocean_blob_cell_update

! namelist variables and default values
logical :: debug_this_module   =.false.  ! extra output for debugging
logical :: really_debug        =.false.  ! lots of extra output for debuugging
logical :: do_bitwise_exact_sum=.false.  ! bitwise exact sums for global sums?
logical :: bitwise_reproduction=.false.  ! bitwise reproduction of results
real    :: blob_small_mass     = 1.0e3   ! minimum mass for a blob to exist (kg)
real    :: max_prop_thickness  = 0.7     ! maximum proportion of a grid cell that the L system may occupy (non-dimensional)

namelist /ocean_blob_nml/ blob_small_mass, debug_this_module, really_debug, &
                          do_bitwise_exact_sum, bitwise_reproduction,       &
                          max_prop_thickness

contains


!######################################################################
! <SUBROUTINE NAME="ocean_blob_init">
!
! <DESCRIPTION>
! Initialises the Lagrangian blob module by setting up module wide
! variables and calling initialisation scripts for the related modules,
! as well as the Lagrangian system itself.
!
! Infrastructure for communicating between PE's, interpolation of E
! system variables to a blob, communicating model wide namelist values
! and variables to other modules and picking up restarts is all done.
! </DESCRIPTION>
!
subroutine ocean_blob_init (Grid, Domain, Time, T_prog, Dens, Thickness,   &
                            L_system, Ext_mode, EL_diag, Ocean_options,    &
                            dtimein, ver_coordinate_class, ver_coordinate, &
                            use_blobs, introduce_blobs)

  type(ocean_grid_type),          intent(in), target :: Grid
  type(ocean_domain_type),        intent(in), target :: Domain
  type(ocean_time_type),          intent(in)         :: Time
  type(ocean_prog_tracer_type),   intent(inout)      :: T_prog(:)
  type(ocean_density_type),       intent(inout)      :: Dens
  type(ocean_thickness_type),     intent(inout)      :: Thickness
  type(ocean_lagrangian_type),    intent(inout)      :: L_system
  type(ocean_external_mode_type), intent(in)         :: Ext_mode
  type(blob_diag_type),           intent(inout)      :: EL_diag(0:)
  type(ocean_options_type),       intent(inout)      :: Ocean_options
  real,                           intent(in)         :: dtimein
  integer,                        intent(in)         :: ver_coordinate_class
  integer,                        intent(in)         :: ver_coordinate
  logical,                        intent(in)         :: use_blobs
  logical,                        intent(in)         :: introduce_blobs

  integer :: nfstatus
  integer :: m, n, i, j, k
  integer :: ioun, ierr, io_status
  integer :: stdoutunit,stdlogunit,stderrunit
  integer :: taup1

  stdoutunit=stdout(); stdlogunit=stdlog() ; stderrunit=stderr()

  taup1 = Time%taup1

  if ( module_is_initialized ) then
    call mpp_error(FATAL,&
    '==>Error in ocean_blob_mod (ocean_blob_init): module already '&
    //'initialized')
  endif

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_blob_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_blob_nml)
  write (stdlogunit,ocean_blob_nml)
  ierr = check_nml_error(io_status,'ocean_blob_nml')
  call close_file (ioun)

  module_is_initialized = .true.

  use_this_module = use_blobs

  if (.not. use_this_module) then
     write(stdoutunit,'(a)') &
          '==>Note: NOT using the Lagrangian buoyancy blobs scheme.'
     Ocean_options%lagrangian_blobs = 'Did NOT use Lagrangian blobs.'
     return
  endif

#ifdef MOM_STATIC_ARRAYS
  write(stderrunit,'(2a)') '==>Error: MOM_STATIC_ARRAYS preprocessor option chosen with use_blobs=.true.', &
       'use_blobs=.true. is not supported with static array allocation'
  call mpp_error(FATAL, &
     '==>Error: MOM_STATIC_ARRAYS preprocessor option chosen with use_blobs=.true. ')

#endif

  if (introduce_blobs .and. Time%init) then
     call mpp_error(WARNING,&
          'introduce_blobs=.TRUE. and Time%init=.TRUE.: introduce_blobs is only for introducing blobs '//&
          'to an existing run. introduce_blobs should usually be .FALSE. if starting from initial '//&
          'conditions.')
  endif
  if (bitwise_reproduction) then
     write(stdoutunit,'(a,2(/,a))')                                                                 &
          '==>Note: bitwise_reproduction=.true. in ocean_blob_nml',                                 &
          '   This is a VERY memory intensive condition.  YOU WILL PROBABLY RUN OUT OF MEMORY!',    &
          '   Unless you require bitwise reproducability for debugging purposes, it is STRONGLY '// &
          'recommended that you set bitwise_reproduction=.false.'
     call mpp_error(WARNING, &
         'bitwise_reproduction=.true., adds heaps to model cost; it is meant only for debugging.')
  endif

  Grd => Grid
  Dom => Domain

  isd=Dom%isd; ied=Dom%ied; jsd=Dom%jsd; jed=Dom%jed
  isc=Dom%isc; iec=Dom%iec; jsc=Dom%jsc; jec=Dom%jec
  isg=Dom%isg; ieg=Dom%ieg; jsg=Dom%jsg; jeg=Dom%jeg
  nk = Grd%nk

  ! For the dynamic blobs, we need to have a halo of size 2
  ! to be able to do the interpolation is implemented here
  call set_ocean_domain(Blob_domain, Grid, xhalo=2, yhalo=2, &
                        name='blobs', maskmap=Dom%maskmap)
  isbd=Blob_domain%isd; iebd=Blob_domain%ied
  jsbd=Blob_domain%jsd; jebd=Blob_domain%jed

  ! Allocate memory to an Information type
  allocate(Info)
  call mpp_get_neighbor_pe(Dom%domain2d, NORTH,      Info%pe_N)
  call mpp_get_neighbor_pe(Dom%domain2d, SOUTH,      Info%pe_S)
  call mpp_get_neighbor_pe(Dom%domain2d, EAST,       Info%pe_E)
  call mpp_get_neighbor_pe(Dom%domain2d, WEST,       Info%pe_W)
  call mpp_get_neighbor_pe(Dom%domain2d, NORTH_EAST, Info%pe_NE)
  call mpp_get_neighbor_pe(Dom%domain2d, NORTH_WEST, Info%pe_NW)
  call mpp_get_neighbor_pe(Dom%domain2d, SOUTH_EAST, Info%pe_SE)
  call mpp_get_neighbor_pe(Dom%domain2d, SOUTH_WEST, Info%pe_SW)
  Info%pe_this = mpp_pe()
  Info%nk_nj = nk*Grd%nj

  ! Variables required for interpolation
  allocate(Info%it(9))
  allocate(Info%jt(9))
  allocate(Info%iu(4))
  allocate(Info%ju(4))
  !1: (i,j), 2: (i+1,j), 3: (i+1,j+1) 4: (i,j+1), 5: (i-1,j+1)
  !          6: (i-1,j), 7: (i-1,j-1), 8:(i,j-1), 9: (i+1,j-1)
  m=1; Info%it(m)= 0; Info%jt(m)= 0
  m=2; Info%it(m)= 1; Info%jt(m)= 0
  m=3; Info%it(m)= 1; Info%jt(m)= 1
  m=4; Info%it(m)= 0; Info%jt(m)= 1
  m=5; Info%it(m)=-1; Info%jt(m)= 1
  m=6; Info%it(m)=-1; Info%jt(m)= 0
  m=7; Info%it(m)=-1; Info%jt(m)=-1
  m=8; Info%it(m)= 0; Info%jt(m)=-1
  m=9; Info%it(m)= 1; Info%jt(m)=-1

  !1: (i,j), 2: (i-1,j), 3: (i-1,j-1), 4:(i,j-1)
  m=1; Info%iu(m)= 0; Info%ju(m)= 0
  m=2; Info%iu(m)=-1; Info%ju(m)= 0
  m=3; Info%iu(m)=-1; Info%ju(m)=-1
  m=4; Info%iu(m)= 0; Info%ju(m)=-1

  allocate(Info%ht(isbd:iebd,jsbd:jebd,2))
  allocate(Info%hu(isbd:iebd,jsbd:jebd,2))
  Info%ht(isc:iec,jsc:jec,1) = Grd%h1t(isc:iec,jsc:jec)
  Info%ht(isc:iec,jsc:jec,2) = Grd%h2t(isc:iec,jsc:jec)
  Info%hu(isc:iec,jsc:jec,1) = Grd%h1u(isc:iec,jsc:jec)
  Info%hu(isc:iec,jsc:jec,2) = Grd%h2u(isc:iec,jsc:jec)
  call mpp_update_domains(Info%ht(:,:,1), Blob_domain%domain2d, complete=.false.)
  call mpp_update_domains(Info%ht(:,:,2), Blob_domain%domain2d, complete=.false.)
  call mpp_update_domains(Info%hu(:,:,1), Blob_domain%domain2d, complete=.false.)
  call mpp_update_domains(Info%hu(:,:,2), Blob_domain%domain2d, complete=.true.)

  ! Arrays that are needed for the horizontal interpolation scheme
  ! to handle solid wall boundaries
  allocate(Info%tidx(0:9,isbd:iebd,jsbd:jebd))
  allocate(Info%uidx(0:4,isbd:iebd,jsbd:jebd))
  ! Index 0 of the first dimension says how many of the surrounding
  ! cells are used in the interpolation.  9 cells for T grid variables
  ! and 4 cells for U grid variables.
  Info%tidx(0,:,:) = 9
  Info%uidx(0,:,:) = 4

  ! The rest of the dimension lists the index of which cells participate
  ! as described below.
  ! Info%tidx:
  !1: (i,j), 2: (i+1,j), 3: (i+1,j+1) 4: (i,j+1), 5: (i-1,j+1)
  !          6: (i-1,j), 7: (i-1,j-1), 8:(i,j-1), 9: (i+1,j-1)
  ! Info%uidx:
  !1: (i,j), 2: (i-1,j), 3: (i-1,j-1), 4:(i,j-1)
  do m=1,4
     Info%tidx(m,:,:) = m
     Info%uidx(m,:,:) = m
  enddo
  do m=5,9
     Info%tidx(m,:,:) = m
  enddo

  ! Now for the exceptions at solid boundaries
  ! We pad the rest of the first dimension with zeros
  ! for values that are unused.

  ! Firstly, do the edges of the compute domain
  if (Info%pe_E==NULL_PE) then
     Info%uidx(0,iec,:) = 2
     Info%tidx(0,iec,:) = 6

     Info%uidx(1,  iec,:) = 2
     Info%uidx(2,  iec,:) = 3
     Info%uidx(3:4,iec,:) = 0

     Info%tidx(1,  iec,:) = 1
     Info%tidx(2,  iec,:) = 4
     Info%tidx(3,  iec,:) = 5
     Info%tidx(4,  iec,:) = 6
     Info%tidx(5,  iec,:) = 7
     Info%tidx(6,  iec,:) = 8
     Info%tidx(7:9,iec,:) = 0
  endif
  if (Info%pe_N==NULL_PE) then
     Info%uidx(0,:,jec) = 2
     Info%tidx(0,:,jec) = 6

     Info%uidx(1,  :,jec) = 3
     Info%uidx(2,  :,jec) = 4
     Info%uidx(3:4,:,jec) = 0

     Info%tidx(1,  :,jec) = 1
     Info%tidx(2,  :,jec) = 2
     Info%tidx(3,  :,jec) = 6
     Info%tidx(4,  :,jec) = 7
     Info%tidx(5,  :,jec) = 8
     Info%tidx(6,  :,jec) = 9
     Info%tidx(7:9,:,jec) = 0
  endif
  if (Info%pe_W==NULL_PE) then
     Info%uidx(0,isc,:) = 2
     Info%tidx(0,isc,:) = 6

     Info%uidx(1,  isc,:) = 1
     Info%uidx(2,  isc,:) = 4
     Info%uidx(3:4,isc,:) = 0

     Info%tidx(1,  isc,:) = 1
     Info%tidx(2,  isc,:) = 2
     Info%tidx(3,  isc,:) = 3
     Info%tidx(4,  isc,:) = 4
     Info%tidx(5,  isc,:) = 8
     Info%tidx(6,  isc,:) = 9
     Info%tidx(7:9,isc,:) = 0
  endif
  if (Info%pe_S==NULL_PE) then
     Info%uidx(0,:,jsc) = 2
     Info%tidx(0,:,jsc) = 6

     Info%uidx(1,  :,jsc) = 1
     Info%uidx(2,  :,jsc) = 2
     Info%uidx(3:4,:,jsc) = 0

     Info%tidx(1,  :,jsc) = 1
     Info%tidx(2,  :,jsc) = 2
     Info%tidx(3,  :,jsc) = 3
     Info%tidx(4,  :,jsc) = 4
     Info%tidx(5,  :,jsc) = 5
     Info%tidx(6,  :,jsc) = 6
     Info%tidx(7:9,:,jsc) = 0
  endif

  ! Then, do the corners
  if (Info%pe_N==NULL_PE .and. Info%pe_E==NULL_PE) then
     Info%uidx(0,iec,jec) = 1
     Info%tidx(0,iec,jec) = 4

     Info%uidx(1,  iec,jec) = 3
     Info%uidx(2:4,iec,jec) = 0

     Info%tidx(1,  iec,jec) = 1
     Info%tidx(2,  iec,jec) = 6
     Info%tidx(3,  iec,jec) = 7
     Info%tidx(4,  iec,jec) = 8
     Info%tidx(5:9,iec,jec) = 0
  endif
  if (Info%pe_N==NULL_PE .and. Info%pe_W==NULL_PE) then
     Info%uidx(0,isc,jec) = 1
     Info%tidx(0,isc,jec) = 4

     Info%uidx(1,  isc,jec) = 4
     Info%uidx(2:4,isc,jec) = 0

     Info%tidx(1,  isc,jec) = 1
     Info%tidx(2,  isc,jec) = 2
     Info%tidx(3,  isc,jec) = 8
     Info%tidx(4,  isc,jec) = 9
     Info%tidx(5:9,isc,jec) = 0
  endif
  if (Info%pe_S==NULL_PE .and. Info%pe_W==NULL_PE) then
     Info%uidx(0,isc,jsc) = 1
     Info%tidx(0,isc,jsc) = 4

     Info%uidx(1,  isc,jsc) = 1
     Info%uidx(2:4,isc,jsc) = 0

     Info%tidx(1,  isc,jsc) = 1
     Info%tidx(2,  isc,jsc) = 2
     Info%tidx(3,  isc,jsc) = 3
     Info%tidx(4,  isc,jsc) = 4
     Info%tidx(5:9,isc,jsc) = 0
  endif
  if (Info%pe_S==NULL_PE .and. Info%pe_E==NULL_PE) then
     Info%uidx(0,iec,jsc) = 1
     Info%tidx(0,iec,jsc) = 4

     Info%uidx(1,  iec,jsc) = 2
     Info%uidx(2:4,iec,jsc) = 0

     Info%tidx(1,  iec,jsc) = 1
     Info%tidx(2,  iec,jsc) = 4
     Info%tidx(3,  iec,jsc) = 5
     Info%tidx(4,  isc,jsc) = 6
     Info%tidx(5:9,iec,jsc) = 0
  endif

  ! Then, zero out the halos that are solid boundaries
  if (Info%pe_E==NULL_PE) then
     Info%uidx(:,ied:iebd,jsc:jec) = 0
     Info%tidx(:,ied:iebd,jsc:jec) = 0
  endif
  if (Info%pe_N==NULL_PE) then
     Info%uidx(:,isc:iec,jed:jebd) = 0
     Info%tidx(:,isc:iec,jed:jebd) = 0
  endif
  if (Info%pe_W==NULL_PE) then
     Info%uidx(:,isbd:isd,jsc:jec) = 0
     Info%tidx(:,isbd:isd,jsc:jec) = 0
  endif
  if (Info%pe_S==NULL_PE) then
     Info%uidx(:,isc:iec,jsbd:jsd) = 0
     Info%tidx(:,isc:iec,jsbd:jsd) = 0
  endif
  if (Info%pe_NE==NULL_PE) then
     Info%uidx(:,ied:iebd,jed:jebd) = 0
     Info%tidx(:,ied:iebd,jed:jebd) = 0
  endif
  if (Info%pe_NW==NULL_PE) then
     Info%uidx(:,isbd:isd,jed:jebd) = 0
     Info%tidx(:,isbd:isd,jed:jebd) = 0
  endif
  if (Info%pe_SW==NULL_PE) then
     Info%uidx(:,isbd:isd,jsbd:jsd) = 0
     Info%tidx(:,isbd:isd,jsbd:jsd) = 0
  endif
  if (Info%pe_SE==NULL_PE) then
     Info%uidx(:,ied:iebd,jsbd:jsd) = 0
     Info%tidx(:,ied:iebd,jsbd:jsd) = 0
  endif

  allocate(L_system%rho_dztup(isd:ied,jsd:jed,nk))
  allocate(L_system%rho_dztlo(isd:ied,jsd:jed,nk))
  allocate(L_system%conv_blob(isd:ied,jsd:jed,nk))

  L_system%rho_dztup(:,:,:) = 0.0
  L_system%rho_dztlo(:,:,:) = 0.0
  L_system%conv_blob(:,:,:) = 0.0

  allocate( datdtimer(isd:ied,jsd:jed) )
  datdtimer(:,:) = 1.0/(Grd%dat(:,:)*dtimein)

  num_prog_tracers  = size(T_prog(:))
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  if (index_temp == -1 .or. index_salt == -1) then
     call mpp_error(FATAL, &
     '==>Error: temp and/or salt not identified in call to '&
     //'ocean_blobs_init')
  endif

  vert_coordinate       = ver_coordinate
  vert_coordinate_class = ver_coordinate_class

  if (vert_coordinate == GEOPOTENTIAL) then
     call mpp_error(FATAL, &
          '==>Error: blobs module has not had geopotential coordinates implemented fully: ocean_blobs_init')
  endif
  if (vert_coordinate == ZSIGMA) then
     call mpp_error(FATAL, &
          '==>Error: blobs module has not had zsigma coordinates implemented: ocean_blobs_init')
  endif
  if (vert_coordinate == PSIGMA) then
     call mpp_error(FATAL, &
          '==>Error: blobs module has not had psigma coordinates implemented: ocean_blobs_init')
  endif
  id_clock_static_bottom     = mpp_clock_id('(Ocean blob: static bottom)           ', grain=CLOCK_ROUTINE)
  id_clock_static_free       = mpp_clock_id('(Ocean blob: static free)             ', grain=CLOCK_ROUTINE)
  id_clock_dyn_free_implicit = mpp_clock_id('(Ocean blob: implicit dynamic free )  ', grain=CLOCK_ROUTINE)
  id_clock_dyn_free_update   = mpp_clock_id('(Ocean blob: update dynamic free )    ', grain=CLOCK_ROUTINE)
  id_clock_dyn_bott_update   = mpp_clock_id('(Ocean blob: update dynamic bottom )  ', grain=CLOCK_ROUTINE)
  id_clock_tsfr_f2b          = mpp_clock_id('(Ocean blob: transfer free to bottom )', grain=CLOCK_ROUTINE)
  id_clock_dyn_bott_new      = mpp_clock_id('(Ocean blob: form new dynamic bottom )', grain=CLOCK_ROUTINE)

  if(debug_this_module) then
     write(stdoutunit,'(a)') &
          '  ==>Note: running ocean_blobs_mod with debug_this_module=.true.'
  endif
  if(really_debug) then
     write(stdoutunit,'(a,/,a)') &
          '  ==>Note: running ocean_blobs_mod with really_debug=.true.',&
          '           Expect LOTS of output!!'
  endif

  if (do_bitwise_exact_sum) then
     global_sum_flag = BITWISE_EXACT_SUM
  else
     global_sum_flag = NON_BITWISE_EXACT_SUM
     if (debug_this_module) then
        write(stdoutunit,'(a)') &
             '==>Note: running ocean_blobs_mod with bitwise_flag=.false. '&
             //'Global sums used for diagnostics will not agree bitwise.'
     endif
  endif

  write(stdoutunit,'(a)') &
       '==>Note: Using the Lagrangian buoyancy blobs scheme.'
  Ocean_options%lagrangian_blobs = 'Did use Lagrangian blobs.'

  call write_version_number( version, tagname )

  ! Register diagnostic fields
  allocate(id_tend_blob(num_prog_tracers))
  do n=1,num_prog_tracers
     if (n==index_temp) then
        id_tend_blob(n) = register_diag_field('ocean_model', 'heat_LtoE', Grd%tracer_axes(1:3), &
             Time%model_time, 'heat transferred from the L to E systems', 'J')
     else
        id_tend_blob(n) = register_diag_field('ocean_model', trim(T_prog(n)%name)//'_LtoE', Grd%tracer_axes(1:3), &
             Time%model_time, trim(T_prog(n)%name)//' transferred from the L to E systems', 'kg')
     endif
  enddo
  id_prop_cell_mass = register_diag_field('ocean_model', 'prop_cell_mass', Grd%tracer_axes(1:3), &
       Time%model_time, 'proportion of grid cell mass taken by L system', 'dimensionless')
  id_mass_in  = register_diag_field('ocean_model', 'mass_in', Grd%tracer_axes(1:3), &
       Time%model_time, 'blob mass entering a grid cell', 'kg')
  id_mass_out = register_diag_field('ocean_model', 'mass_out', Grd%tracer_axes(1:3), &
       Time%model_time, 'blob mass exiting a grid cell', 'kg')
  id_dstryd   = register_diag_field('ocean_model', 'dstryd_blobs', Time%model_time, &
       'blobs destroyed because dztL>dztT', 'number of blobs')
  id_grounded = register_diag_field('ocean_model', 'grounded_blobs', Time%model_time, &
       'blobs destroyed due to grounding', 'number of blobs')
  id_surfaced = register_diag_field('ocean_model', 'surfaced_blobs', Time%model_time, &
       'blobs destroyed due to surfacing', 'number of blobs')
  id_detraind = register_diag_field('ocean_model', 'detrained_blobs', Time%model_time, &
       'blobs detrained to a small mass', 'number of blobs')
  id_nblob_free = register_diag_field('ocean_model', 'nblob_free', Time%model_time, &
       'number of dynamic free blobs')
  id_nblob_bott = register_diag_field('ocean_model', 'nblob_bott', Time%model_time, &
       'number of dynamic bottom blobs')

  ! Diagnostic arrays for mass convergence and divergence
  allocate(mass_in( isd:ied,jsd:jed,1:nk))
  allocate(mass_out(isd:ied,jsd:jed,1:nk))

  ! Diagnostic structure for transfer of properties between E and L systems
  allocate(id_entrainment(0:num_prog_tracers))
  allocate(id_detrainment(0:num_prog_tracers))
  allocate(id_new(        0:num_prog_tracers))
  allocate(id_dstry(      0:num_prog_tracers))

  !The 0 is for mass
  do n=0,num_prog_tracers
     ! Allocate the arrays
     allocate(EL_diag(n)%entrainment(isd:ied,jsd:jed,1:nk))
     allocate(EL_diag(n)%detrainment(isd:ied,jsd:jed,1:nk))
     allocate(EL_diag(n)%new(        isd:ied,jsd:jed,1:nk))
     allocate(EL_diag(n)%dstry(      isd:ied,jsd:jed,1:nk))

     EL_diag(n)%entrainment(:,:,:) = 0.0
     EL_diag(n)%detrainment(:,:,:) = 0.0
     EL_diag(n)%new(:,:,:)         = 0.0
     EL_diag(n)%dstry(:,:,:)       = 0.0

     if (n==0) then
        id_entrainment(n) = register_diag_field('ocean_model', 'blob_mass_ent',   Grd%tracer_axes(1:3), &
             Time%model_time, 'Mass transferred from the E system to the L system due to entrainment by blobs', 'kg')
        id_detrainment(n) = register_diag_field('ocean_model', 'blob_mass_det',   Grd%tracer_axes(1:3), &
             Time%model_time, 'Mass transferred from the L system to the E system due to detrainment by blobs', 'kg')
        id_new(n)         = register_diag_field('ocean_model', 'blob_mass_new',   Grd%tracer_axes(1:3), &
             Time%model_time, 'Mass transferred from the E system to the L system due to blob formation',       'kg')
        id_dstry(n)       = register_diag_field('ocean_model', 'blob_mass_dstry', Grd%tracer_axes(1:3), &
             Time%model_time, 'Mass transferred from the L system to the E system due to blob destruction',     'kg')
     elseif (n==index_temp) then
        id_entrainment(n) = register_diag_field('ocean_model', 'blob_heat_ent',   Grd%tracer_axes(1:3), &
             Time%model_time, 'Heat transferred from the E system to the L system due to entrainment by blobs', 'J')
        id_detrainment(n) = register_diag_field('ocean_model', 'blob_heat_det',   Grd%tracer_axes(1:3), &
             Time%model_time, 'Heat transferred from the L system to the E system due to detrainment by blobs', 'J')
        id_new(n)         = register_diag_field('ocean_model', 'blob_heat_new',   Grd%tracer_axes(1:3), &
             Time%model_time, 'Heat transferred from the E system to the L system due to blob formation',       'J')
        id_dstry(n)       = register_diag_field('ocean_model', 'blob_heat_dstry', Grd%tracer_axes(1:3), &
             Time%model_time, 'Heat transferred from the L system to the E system due to blob destruction',     'J')
     else
        id_entrainment(n) = register_diag_field('ocean_model', 'blob_'//trim(T_prog(n)%name)//'_ent',   Grd%tracer_axes(1:3), &
             Time%model_time, trim(T_prog(n)%name)//' transferred from the E system to the L system due to entrainment by blobs', 'kg')
        id_detrainment(n) = register_diag_field('ocean_model', 'blob_'//trim(T_prog(n)%name)//'_det',   Grd%tracer_axes(1:3), &
             Time%model_time, trim(T_prog(n)%name)//' transferred from the L system to the E system due to detrainment by blobs', 'kg')
        id_new(n)         = register_diag_field('ocean_model', 'blob_'//trim(T_prog(n)%name)//'_new',   Grd%tracer_axes(1:3), &
             Time%model_time, trim(T_prog(n)%name)//' transferred from the E system to the L system due to blob formation',       'kg')
        id_dstry(n)       = register_diag_field('ocean_model', 'blob_'//trim(T_prog(n)%name)//'_dstry', Grd%tracer_axes(1:3), &
             Time%model_time, trim(T_prog(n)%name)//' transferred from the L system to the E system due to blob destruction',     'kg')
     endif
  enddo

  ! Initialise the related blob modules, which actually run the various
  ! available parameterisations
  call blob_util_init(Grd, Dom, Info, Blob_domain, global_sum_flag, num_prog_tracers,&
                      index_temp, index_salt, blob_small_mass, dtimein,              &
                      bitwise_reproduction, vert_coordinate_class, vert_coordinate,  &
                      debug_this_module, really_debug)

  ! Initialise the diagnostics
  call blob_diag_init(T_prog(:), Grd, Dom, num_prog_tracers, dtimein, index_temp, &
                      vert_coordinate_class, blob_diagnostics)

  call blob_static_free_init(Grd, Dom, num_prog_tracers, index_temp, index_salt)

  call blob_static_bottom_init(Grd, Dom, Info, num_prog_tracers, index_temp, index_salt, dtimein)

  call blob_dynamic_bottom_init(Time, Grd, Dom, Blob_domain, Info,                        &
                                debug_this_module,really_debug, bitwise_reproduction,     &
                                num_prog_tracers, index_temp, index_salt, dtimein,        &
                                vert_coordinate_class, vert_coordinate, blob_diagnostics, &
                                bott_minstep, bott_total_ns, blob_small_mass, use_dyn_bot)

  call blob_dynamic_free_init(Grd, Dom, Time, T_prog(:), Blob_domain, Info, debug_this_module,  &
                              really_debug, bitwise_reproduction, num_prog_tracers, index_temp, &
                              index_salt, dtimein,vert_coordinate_class, vert_coordinate,       &
                              blob_diagnostics, free_minstep, free_total_ns, blob_small_mass, use_dyn_fre)

  ! Read in the restart files
  call readrestart

  call init_blob_thickness(Time, Thickness, T_prog(:), L_system)

  ! Calculate the total (combined L and E) tracer
  do n=1,num_prog_tracers
     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(n)%fieldT(i,j,k) = T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1) &
                                        + T_prog(n)%sum_blob(i,j,k,taup1)*Grd%datr(i,j)
              T_prog(n)%fieldT(i,j,k) = T_prog(n)%fieldT(i,j,k)/Thickness%rho_dztT(i,j,k,taup1)
           enddo
        enddo
     enddo
     call mpp_update_domains(T_prog(n)%fieldT(:,:,:), Dom%domain2d, complete=T_prog(n)%complete)
  enddo

  ! Compute the buoyancy frequency based on the full system
  call compute_buoyfreq(Time, Thickness, T_prog(index_salt)%fieldT(:,:,:), &
                        T_prog(index_temp)%fieldT(:,:,:), Dens, use_blobs)

  if (debug_this_module) then
     write(stdoutunit,'(/a)') 'From ocean_blobs_init: debug blob chksums'
     call write_timestamp(Time%model_time)
     call write_chksum_3d('Dens%rhoT', Grd%tmask(COMP,:)*Dens%rhoT(COMP,:))
     do n=1,num_prog_tracers
        call write_chksum_3d(trim(T_prog(n)%name)//' fieldT', Grd%tmask(COMP,:)*T_prog(n)%fieldT(COMP,:))
     enddo
     call blob_chksum(T_prog(:), head_static_free, head_static_bott, &
          head_dynamic_free, head_dynamic_bott, blob_counter)
     call lagrangian_system_chksum(L_system)

     if (really_debug) call write_all_blobs('taup1')
  endif


contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine for checking netcdf read errors          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine nferror(description)
    character(len=*) :: description
    if (nfstatus /= NF_NOERR) then
       write (stdoutunit, *) ' '
       write (stdoutunit, *) 'ocean_blob_util_mod, blob_util_init: problem '//trim(description)
       write (stdoutunit, *) 'error code =', nfstatus
       if (nfstatus==NF_EEDGE)        write(stdoutunit, *) '==>Start+count exceeds dimension bound'
       if (nfstatus==NF_EINVALCOORDS) write(stdoutunit, *) '==>Index exceeds dimension bound'
       if (nfstatus==NF_ENOTVAR)      write(stdoutunit, *) '==>Variable not found'

       call mpp_error(FATAL, 'ocean_blob_mod, ocean_blob_init: problem '//trim(description))
    endif
  end subroutine nferror

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that reads the restart files for blobs   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine readrestart

    type(ocean_blob_type), pointer        :: blob
    integer, allocatable, dimension(:)    :: tracerid
    integer, dimension(1)                 :: id_restart
    character(len=128) :: filename
    integer :: ierr, ncid, dimid, nblobs_in_file
    integer :: latid, lonid, geodepthid, massid, uid, vid, wid, stid
    integer :: dragid, entid, ageid, heightid
    integer :: hashid, numberid, typeid, iid, jid, kid
    integer :: stepid, h1id, h2id, nstepsid, mstepsid, densityid
    integer :: n, m, i, j, k, type, nblobs, fileno
    logical :: found_restart, distributed, next_file_exists

    ! We have two restart files, and we need to read them in separately.
    ! The first is to read in the gridded data.  The second is to read
    ! in the blobs and distribute them amongst the processors.

    ! Allocate memory to the blob counter
    allocate(blob_counter(isc:iec,jsc:jec,nk))
    blob_counter(:,:,:)    = 0

    ! Read in the gridded data
    filename = 'ocean_blobs_gridded.res.nc'
    id_restart(1) = register_restart_field(gridded_restart, filename, 'blob_counter', &
                                           blob_counter(:,:,:), domain=Dom%domain2d)

    if (.not. introduce_blobs) then
       if(file_exist('INPUT/'//trim(filename))) then
          write(stdoutunit, '(/a)') 'Reading in some gridded blobs data from '//trim(filename)
          call restore_state( gridded_restart, id_restart(1) )

       else
          if (.NOT.Time%init) then
             call mpp_error(FATAL,&
                  'Expecting file '//trim(filename)//' to exist.&
                  &This file was not found and Time%init=.false.')
          endif
       endif

       ! Open the blob file and read in the blob attributes.
       filename = 'INPUT/ocean_blobs.res.nc'
       inquire(file=trim(filename),exist=found_restart)
       if (found_restart) then
          distributed=.false.
          next_file_exists=.true.
       else
          write(filename(1:29), '("INPUT/ocean_blobs.res.nc.", I4.4)') 0
          inquire(file=trim(filename),exist=found_restart)
          if (found_restart) then
             distributed=.true.
             next_file_exists=.true.
          else
             if (.NOT.Time%init) then
                call mpp_error(FATAL,&
                     'Expecting file '//trim(filename)//' or similar to exist.&
                     &This file was not found and Time%init=.false.')
             endif
          endif
       endif

       fileno=0
       nblobs=0

       do while (next_file_exists)
          write(stdoutunit, '(/,a,/)') 'Reading in blobs restart from '//trim(filename)

          ierr = nf_open(trim(filename), NF_NOWRITE, ncid)
          if (ierr .ne. NF_NOERR) write(stderrunit,'(a)') 'blobs, readrestart: nf_open failed'

          ierr = nf_inq_unlimdim(ncid, dimid)
          if (ierr .ne. NF_NOERR) write(stderrunit,'(a)') 'blobs, readrestart: nf_inq_unlimdim failed'

          ierr = nf_inq_dimlen(ncid, dimid, nblobs_in_file)
          if (ierr .ne. NF_NOERR) write(stderrunit,*) 'blobs, readrestart: nf_inq_dimlen failed'

          iid        = inq_var(ncid, 'i')
          jid        = inq_var(ncid, 'j')
          kid        = inq_var(ncid, 'k')
          latid      = inq_var(ncid, 'lat')
          lonid      = inq_var(ncid, 'lon')
          geodepthid = inq_var(ncid, 'geodepth')
          stid       = inq_var(ncid, 'st')
          massid     = inq_var(ncid, 'mass')
          densityid  = inq_var(ncid, 'density')
          uid        = inq_var(ncid, 'u')
          vid        = inq_var(ncid, 'v')
          wid        = inq_var(ncid, 'w')
          entid      = inq_var(ncid, 'ent')
          dragid     = inq_var(ncid, 'drag')
          ageid      = inq_var(ncid, 'age')
          heightid   = inq_var(ncid, 'height')
          h1id       = inq_var(ncid, 'h1')
          h2id       = inq_var(ncid, 'h2')
          stepid     = inq_var(ncid, 'step')
          hashid     = inq_var(ncid, 'hash')
          numberid   = inq_var(ncid, 'number')
          typeid     = inq_var(ncid, 'type')
          nstepsid   = inq_var(ncid, 'nsteps')
          mstepsid   = inq_var(ncid, 'model_steps')
          if (fileno==0) allocate(tracerid(num_prog_tracers))
          do n=1, num_prog_tracers
             if (n==index_temp) then
                tracerid(n) = inq_var(ncid, 'heat')
             else
                tracerid(n) = inq_var(ncid, trim(T_prog(n)%name))
             endif
          enddo

          ! Cycle through the blobs.  Ignore empty entries and blobs that are not
          ! in this compute domain.  If a blob is on this compute domain, then
          ! we read in the saved data, and derive other data (such as field,
          ! grid cell etc.).  Then, we put it into the linked list.
          do m=1,nblobs_in_file
             ! get horizontal coordinate information of the blobs and type information
             type = get_int(ncid, typeid, m)
             i = get_int(ncid, iid, m)
             j = get_int(ncid, jid, m)

             ! if type==0, it is an empty entry, so we just ignore it
             if(type /= 0) then
                ! Check if we are in the compute domain. If we are not,
                ! then ignore the blob.  If we are, then read in the
                ! blob's data.
                if (isc<=i .and. i<=iec .and. jsc<=j .and. j<=jec) then
                   nblobs = nblobs + 1
                   allocate(blob)
                   blob%lon         = get_double(ncid, lonid,      m)
                   blob%lat         = get_double(ncid, latid,      m)
                   blob%geodepth    = get_double(ncid, geodepthid, m)
                   blob%st          = get_double(ncid, stid,       m)
                   blob%mass        = get_double(ncid, massid,     m)
                   blob%density     = get_double(ncid, densityid,  m)
                   blob%v(1)        = get_double(ncid, uid,        m)
                   blob%v(2)        = get_double(ncid, vid,        m)
                   blob%v(3)        = get_double(ncid, wid,        m)
                   blob%ent         = get_double(ncid, entid,      m)
                   blob%drag        = get_double(ncid, dragid,     m)
                   blob%age         = get_double(ncid, ageid,      m)
                   blob%height      = get_double(ncid, heightid,   m)
                   blob%h1          = get_double(ncid, h1id,       m)
                   blob%h2          = get_double(ncid, h2id,       m)
                   blob%step        = get_double(ncid, stepid,     m)
                   blob%hash        = get_int(ncid,    hashid,     m)
                   blob%number      = get_int(ncid,    numberid,   m)
                   blob%nsteps      = get_int(ncid,    nstepsid,   m)
                   blob%model_steps = get_int(ncid,    mstepsid,   m)
                   blob%k           = get_int(ncid,    kid,        m)
                   blob%i           = i
                   blob%j           = j
                   allocate(blob%tracer(num_prog_tracers))
                   allocate(blob%field(num_prog_tracers))
                   do n=1,num_prog_tracers
                      blob%tracer(n) = get_double(ncid, tracerid(n), m)
                      blob%field(n)  = blob%tracer(n) / blob%mass
                   enddo

                   i = blob%i
                   j = blob%j
                   k = blob%k

                   ! new blobs were created implicitly in time, and thus, their
                   ! properties belong to the next time step and have not been
                   ! taken into account in diagnostics
                   if(blob%new) then
                      EL_diag(0)%new(i,j,k) = EL_diag(0)%new(i,j,k) + blob%mass
                      do n=1,num_prog_tracers
                         EL_diag(n)%new(i,j,k) = EL_diag(n)%new(i,j,k) + blob%tracer(n)
                      enddo
                   endif

                   blob%depth = blob%geodepth + Ext_mode%eta_t(i,j,taup1)
                   blob%new   = .false.

                   blob%densityr = 1./blob%density

                   if (vert_coordinate_class == DEPTH_BASED) then
                      blob%volume = blob%mass * rho0r
                   else
                      blob%volume = blob%mass * blob%densityr
                   endif

                   ! Insert the blob into its correct place in the linked
                   ! Give an error if the type is unrecognised. type==1 is
                   ! a static free blob and type==2 is a static bottom blob,
                   ! type==3 is a dynamic free blob and type==4 is a dynamic
                   ! bottom blob.  type==0 means that it is an empty entry.
                   if (type==1) then
                      call insert_blob(blob, head_static_free)
                   elseif(type==2) then
                      call insert_blob(blob, head_static_bott)
                   elseif(type==3) then
                      call allocate_interaction_memory(blob, free_total_ns)
                      call insert_blob(blob, head_dynamic_free)
                   elseif(type==4) then
                      call allocate_interaction_memory(blob, bott_total_ns)
                      call insert_blob(blob, head_dynamic_bott)
                   else
                      write(stdoutunit,'(a)') &
                           'blob, readrestart, blob is of unidentified type'
                      call mpp_error(FATAL, &
                           'blob, readrestart, blob is of unidentified type')
                   endif

                   nullify(blob)
                endif !on this pe
             endif !type == 0
          enddo !n=1,nblobs_in_file

          if (distributed) then
             fileno=fileno+1
             write(filename(1:29), '("INPUT/ocean_blobs.res.nc.", I4.4)') fileno
             inquire(file=trim(filename),exist=found_restart)
             if (.not. found_restart) next_file_exists=.false.
          else
             next_file_exists=.false.
          endif

       enddo

       print('(a22,i3,a3,i10)'),  'number of blobs on PE ',Info%pe_this,' =',nblobs
       call mpp_sum(nblobs)
       write(stdoutunit,'(a,i10)') 'global number of blobs     =', nblobs

    endif

  end subroutine readrestart

end subroutine ocean_blob_init
! </SUBROUTINE>  NAME="ocean_blob_init"


!######################################################################
! <SUBROUTINE NAME="init_blob_thickness">
!
! <DESCRIPTION>
! Initialises the L_system thickness, based on the existing blobs
! </DESCRIPTION>
!
subroutine init_blob_thickness(Time, Thickness, T_prog, L_system)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_lagrangian_type),  intent(inout) :: L_system

  type(ocean_blob_type), pointer :: this
  real    :: dzt, rho_dzt
  integer :: i,j,k,n,taup1
  integer :: stdoutunit

  if (.not. use_this_module) return

  stdoutunit = stdout()
  taup1 = Time%taup1

     this=>head_dynamic_free

     if(associated(this)) then
        freecycle: do
           ! for convenience
           i  = this%i - Dom%ioff
           j  = this%j - Dom%joff
           k  = this%k

           ! We use the EOS density regardless of the coordinates, since
           ! the value of rho_dzt for the L-system is used for the
           ! calculation of hydrostatic pressure.
           dzt     = Grd%datr(i,j)*this%volume
           rho_dzt = dzt*this%density
           if (vert_coordinate_class==DEPTH_BASED) then
              if (this%st>-Thickness%depth_st(i,j,k)) then
                 L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt
              else
                 L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt
              endif
           else !PRESSURE_BASED
              if (this%st<Thickness%depth_st(i,j,k)) then
                 L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt
              else
                 L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt
              endif
           endif
           ! Note, we do not need to initialise the L system thickness, as this is already
           ! done in the thickness module from restarts

           ! The gridded tracer data
           do n=1,num_prog_tracers
              T_prog(n)%sum_blob(i,j,k,taup1) = T_prog(n)%sum_blob(i,j,k,taup1) + this%tracer(n)
           enddo

           this=>this%next
           if(.not.associated(this)) exit freecycle
        enddo freecycle
     endif

     this=>head_dynamic_bott
     if(associated(this)) then
        bottcycle: do
           ! for convenience
           i  = this%i - Dom%ioff
           j  = this%j - Dom%joff
           k  = this%k

           ! We use the EOS density regardless of the coordinates, since
           ! the value of rho_dzt for the L-system is used for the
           ! calculation of hydrostatic pressure.  Bottom blobs only
           ! contribute to the lower half of the grid cell thickness.
           rho_dzt = Grd%datr(i,j)*this%volume*this%density
           L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt*0.5
           L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt*0.5

           ! Note, we do not need to initialise the L system thickness, as this is already
           ! done in the thickness module from restarts

           ! The gridded tracer data
           do n=1,num_prog_tracers
              T_prog(n)%sum_blob(i,j,k,taup1) = T_prog(n)%sum_blob(i,j,k,taup1) + this%tracer(n)
           enddo

           this=>this%next
           if(.not.associated(this)) exit bottcycle
        enddo bottcycle
     endif

     call mpp_update_domains(L_system%rho_dztlo(:,:,:), Dom%domain2d, complete=.false.)
     call mpp_update_domains(L_system%rho_dztup(:,:,:), Dom%domain2d, complete=.false.)
     do n=1,num_prog_tracers
        call mpp_update_domains(T_prog(n)%sum_blob(:,:,:,taup1), Dom%domain2d, complete=T_prog(n)%complete)
     enddo

end subroutine init_blob_thickness
! </SUBROUTINE>  NAME="init_blob_thickness"


!#######################################################################
! <SUBROUTINE NAME="ocean_blob_update">
!
! <DESCRIPTION>
! Updates the Lagrangian blobs by calling routines that run during the
! explicit in time part of the time step.
! </DESCRIPTION>
!
subroutine ocean_blob_update(Time, Thickness, T_prog, Dens, Ext_mode, Adv_vel, Velocity, L_system, EL_diag)
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_adv_vel_type),       intent(inout) :: Adv_vel
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  type(blob_diag_type),           intent(inout) :: EL_diag(0:)

  type(ocean_blob_type), pointer :: this=>NULL()
  type(ocean_blob_type), pointer :: head=>NULL()
  type(ocean_blob_type), pointer :: free2bott=>NULL()
  type(ocean_blob_type), pointer :: bott2free=>NULL()
  integer :: tau, taup1
  integer :: k, i, j, n, p, s
  integer :: ngrounded, nsurfaced, nblob, ndetraind
  real, dimension(num_prog_tracers,isd:ied,jsd:jed,nk) :: tend_blob
  real, dimension(isd:ied,jsd:jed)                     :: blob_source
  real, dimension(isbd:iebd,jsbd:jebd,nk,2)            :: press_grad, u
  real, dimension(isbd:iebd,jsbd:jebd,nk)              :: rho
  real, dimension(isbd:iebd,jsbd:jebd,0:nk)            :: w

  if (.not. use_this_module) return

  tau   = Time%tau
  taup1 = Time%taup1

  if (.not. module_is_initialized ) then
     call mpp_error(FATAL, &
          '==>Error in ocean_blob_mod (Lagrangian blob model): module '&
          //' needs to be initialized')
  endif

  ! Note, that in order to help reduce roundoff error, we add the blob mass to
  ! conv_blob in the blob modules.  Then, at the end of ocean_blob_cell_update, we
  ! divide by (dat*dtime)
  Ext_mode%conv_blob(:,:)   = 0.0
  L_system%conv_blob(:,:,:) = 0.0
  blob_source(:,:)          = 0.0
  tend_blob(:,:,:,:)        = 0.0

  ! For diagnostics
  mass_in(:,:,:)  = 0.0
  mass_out(:,:,:) = 0.0
  ngrounded = 0
  nsurfaced = 0
  ndetraind = 0

  ! We need to mpp_update these variables for the dynamic blobs if we wish to maintain
  ! bitwise reproduction across processors.
  press_grad(:,:,:,:) = 0.0
  u(:,:,:,:)          = 0.0
  rho(:,:,:)          = 0.0
  w(:,:,:)            = 0.0

  if (Grd%tripolar) then
     do k=1,nk
        press_grad(isc:iec,jsc:jec,k,1) = Grd%cos_rot(isc:iec,jsc:jec) * Velocity%press_force(isc:iec,jsc:jec,k,1) + &
                                          Grd%sin_rot(isc:iec,jsc:jec) * Velocity%press_force(isc:iec,jsc:jec,k,2)
        press_grad(isc:iec,jsc:jec,k,2) = Grd%sin_rot(isc:iec,jsc:jec) * Velocity%press_force(isc:iec,jsc:jec,k,1) + &
                                          Grd%cos_rot(isc:iec,jsc:jec) * Velocity%press_force(isc:iec,jsc:jec,k,2)

        u(isc:iec,jsc:jec,k,1) = Grd%cos_rot(isc:iec,jsc:jec) * Velocity%u(isc:iec,jsc:jec,k,1,tau) + &
                                 Grd%sin_rot(isc:iec,jsc:jec) * Velocity%u(isc:iec,jsc:jec,k,2,tau)
        u(isc:iec,jsc:jec,k,2) = Grd%sin_rot(isc:iec,jsc:jec) * Velocity%u(isc:iec,jsc:jec,k,1,tau) + &
                                 Grd%cos_rot(isc:iec,jsc:jec) * Velocity%u(isc:iec,jsc:jec,k,2,tau)
     enddo
  else
     press_grad(isc:iec,jsc:jec,1:nk,1:2) = Velocity%press_force(isc:iec,jsc:jec,1:nk,1:2)
     u(         isc:iec,jsc:jec,1:nk,1:2) = Velocity%u(          isc:iec,jsc:jec,1:nk,1:2,tau)
  endif

  rho(isc:iec,jsc:jec,1:nk) = Dens%rho(isc:iec,jsc:jec,1:nk,tau)

  w(isc:iec,jsc:jec,0)    = Ext_mode%deta_dt(isc:iec,jsc:jec)
  w(isc:iec,jsc:jec,1:nk) = rho0r*Adv_vel%wrho_bt(isc:iec,jsc:jec,1:nk)

  call mpp_update_domains(press_grad(:,:,:,1), press_grad(:,:,:,2),Blob_domain%domain2d,gridtype=BGRID_NE, complete=.false.)
  call mpp_update_domains(u(:,:,:,1), u(:,:,:,2),                  Blob_domain%domain2d,gridtype=BGRID_NE, complete=.true.)
  call mpp_update_domains(rho(:,:,:), Blob_domain%domain2d)
  call mpp_update_domains(w(:,:,:),   Blob_domain%domain2d)

  if (id_nblob_free>0) then
     nblob = 0
     this=>head_dynamic_free
     do while(associated(this))
        nblob = nblob + 1
        this=>this%next
     enddo
     call mpp_sum(nblob)
     used = send_data(id_nblob_free, real(nblob), Time%model_time)
  endif

  if (id_nblob_bott>0) then
     nblob = 0
     this=>head_dynamic_bott
     do while(associated(this))
        nblob = nblob + 1
        this=>this%next
     enddo
 call mpp_sum(nblob)
     used = send_data(id_nblob_bott, real(nblob), Time%model_time)
  endif

  ! call the explicit in time parts of the parameterisations
  ! Static bottom blobs
  call mpp_clock_begin(id_clock_static_bottom)
  call blob_overflow_like(Time, Thickness, T_prog(:), Dens, head_static_bott, blob_counter)
  call mpp_clock_end(id_clock_static_bottom)

  ! Dynamic free blobs
  call mpp_clock_begin(id_clock_dyn_free_update)
  call blob_dynamic_free_update(Time, Thickness, T_prog(:), Ext_mode, Dens, L_system, &
                                tend_blob, blob_source, press_grad, u, w, rho,        &
                                head_dynamic_free, free2bott, mass_in, mass_out,      &
                                EL_diag(:), ngrounded, nsurfaced, ndetraind)
  call mpp_clock_end(id_clock_dyn_free_update)

  ! Dynamic bottom blobs
  call mpp_clock_begin(id_clock_dyn_bott_update)
  call blob_dynamic_bottom_update(head_dynamic_bott, bott2free, Time, Dens, Thickness, &
                                  T_prog, Ext_mode, L_system, tend_blob, blob_source,  &
                                  u, rho, mass_in, mass_out, EL_diag(:), ngrounded, ndetraind)
  call mpp_clock_end(id_clock_dyn_bott_update)

  call mpp_sum(ngrounded)
  call mpp_sum(nsurfaced)
  call mpp_sum(ndetraind)
  if (id_grounded>0) used = send_data(id_grounded, real(ngrounded), Time%model_time)
  if (id_surfaced>0) used = send_data(id_surfaced, real(nsurfaced), Time%model_time)
  if (id_detraind>0) used = send_data(id_detraind, real(ndetraind), Time%model_time)

  if (debug_this_module .and. bitwise_reproduction) then
     call entrainment_checksum(Time, head_dynamic_free, head_dynamic_bott, free2bott, bott2free)
  endif

  ! Bitwise reproduction is a very memory intensive operation and should only be used
  ! for debugging.  For bitwise_reproduction=.true. we need to process blobs and their
  ! histories in the same relative order regardless of domain decomposition and restarts.
  ! To do so, we save the "history" of each blob subcycle to a number of arrays (which
  ! can be a very memory intensive process) and process them in order.  When we set
  ! bitwise_reproduction=.false. we process the blob histories "in situ" without regard
  ! to order and thus obviate the need for the history arrays and consequently greatly
  ! reduce the memory requirements of the Lagrangian model.
  if (bitwise_reproduction) then
     ! Accumulate the mass transferred between the E and L
     ! systems as a result of entrainment and detrainment.
     ! Also, calculate the divergence as a result of the blobs.
     do p=3,6
        if(p==3) head=>head_dynamic_free
        if(p==4) head=>head_dynamic_bott
        if(p==5) head=>free2bott
        if(p==6) head=>bott2free
        if (associated(head)) then
           this=>head

           blobcycle: do
              i=this%di(0)
              j=this%dj(0)
              k=this%dk(0)

              if (isc<= i .and. i<=iec .and. jsc<= j .and. j<=jec) then

                 Ext_mode%conv_blob(i,j)   = Ext_mode%conv_blob(i,j)   - this%mass_out(0)
                 L_system%conv_blob(i,j,k) = L_system%conv_blob(i,j,k) - this%mass_out(0)

                 mass_out(i,j,k) = mass_out(i,j,k) + this%mass_out(0)
              endif

              if (this%nfrac_steps>0) then
                 do s=1,this%nfrac_steps
                    i=this%di(s)
                    j=this%dj(s)
                    k=this%dk(s)

                    if (isc<= i .and. i<=iec .and. jsc<= j .and. j<=jec) then

                       Ext_mode%conv_blob(i,j)   = Ext_mode%conv_blob(i,j)   - this%mass_out(s)
                       Ext_mode%conv_blob(i,j)   = Ext_mode%conv_blob(i,j)   + this%mass_in(s)

                       L_system%conv_blob(i,j,k) = L_system%conv_blob(i,j,k) - this%mass_out(s)
                       L_system%conv_blob(i,j,k) = L_system%conv_blob(i,j,k) + this%mass_in(s)

                       mass_out(i,j,k) = mass_out(i,j,k) + this%mass_out(s)
                       mass_in(i,j,k)  = mass_in(i,j,k)  + this%mass_in(s)

                       blob_source(i,j) = blob_source(i,j) - this%entrainment(s,0) - this%detrainment(s,0)
                       EL_diag(0)%entrainment(i,j,k) = EL_diag(0)%entrainment(i,j,k) + this%entrainment(s,0)
                       EL_diag(0)%detrainment(i,j,k) = EL_diag(0)%detrainment(i,j,k) - this%detrainment(s,0)

                       do n=1,num_prog_tracers
                          tend_blob(n,i,j,k) = tend_blob(n,i,j,k) - this%entrainment(s,n) - this%detrainment(s,n)
                          EL_diag(n)%entrainment(i,j,k) = EL_diag(n)%entrainment(i,j,k) + this%entrainment(s,n)
                          EL_diag(n)%detrainment(i,j,k) = EL_diag(n)%detrainment(i,j,k) - this%detrainment(s,n)
                       enddo
                    endif
                 enddo
              endif

              ! dmass and dtracer are for when a blob is destroyed (e.g. grounded, penetrated free surface)
              ! So, they should take the last recorded value of (i,j,k)
              i=this%i
              j=this%j
              k=this%k
              if (isc<= i .and. i<=iec .and. jsc<= j .and. j<=jec) then
                 blob_source(i,j)   = blob_source(i,j) - this%dmass
                 EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) - this%dmass
                 do n=1,num_prog_tracers
                    tend_blob(n,i,j,k) = tend_blob(n,i,j,k) - this%dtracer(n)
                    EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) - this%dtracer(n)
                 enddo
              endif
              this=>this%next
              if (.not. associated(this)) exit blobcycle
           enddo blobcycle
        endif
     enddo!p
  endif

  ! Trasfer blobs between lists:
  ! Free blobs that have interacted with topography.
  call mpp_clock_begin(id_clock_tsfr_f2b)
  call transfer_free_to_bottom(Time, Thickness, T_prog, free2bott, head_dynamic_bott, &
                               use_dyn_fre, EL_diag(:))
  call mpp_clock_end(id_clock_tsfr_f2b)

  ! Bottom blobs that have separated.
  call mpp_clock_begin(id_clock_tsfr_f2b)
  call transfer_bottom_to_free(Time, Thickness, T_prog(:), Dens, bott2free, &
                               head_dynamic_free, use_dyn_bot, EL_diag(:))
  call mpp_clock_end(id_clock_tsfr_f2b)

  ! Remove all blobs from dynamic free and bottom lists that have a small mass
  call blob_delete(Time, Thickness, T_prog(:), head_dynamic_free)
  call blob_delete(Time, Thickness, T_prog(:), head_dynamic_bott)

  ! Form new blobs as a result of the on-shelf/off-shelf creation condition
  call mpp_clock_begin(id_clock_dyn_bott_new)
  call dynamic_bottom_form_new(Time, Dens, T_prog, Thickness, Ext_mode, L_system, tend_blob, &
                               Blob_counter, u, mass_out, mass_in, head_dynamic_bott, EL_diag(:))
  call mpp_clock_end(id_clock_dyn_bott_new)

  ! We only finalise the horizontal (depth integrated) blob convergence.
  ! As we have yet to find the position of the vertical coordinate system
  ! for taup1, we do not yet know what the final cell a blob lies in.
  do j=jsd,jed
     do i=isd,ied
        Ext_mode%conv_blob(i,j) = Ext_mode%conv_blob(i,j)*datdtimer(i,j)
     enddo
  enddo
  call mpp_update_domains(Ext_mode%conv_blob(:,:), Dom%domain2d)

  Thickness%blob_source(:,:) = Thickness%blob_source(:,:) + blob_source(:,:)*datdtimer(:,:)
  do n=1,num_prog_tracers
     do k=1,nk
        T_prog(n)%tend_blob(:,:,k) = T_prog(n)%tend_blob(:,:,k) + tend_blob(n,:,:,k)*datdtimer(:,:)
     enddo
     call mpp_update_domains(T_prog(n)%tend_blob(:,:,:), Dom%domain2d, complete=T_prog(n)%complete)
  enddo

end subroutine ocean_blob_update
! </SUBROUTINE>  NAME="ocean_blob_update"

!#######################################################################
! <SUBROUTINE NAME="ocean_blob_cell_update">
!
! <DESCRIPTION>
!
! For bottom blobs, it diagnoses the blobs vertical position in the
! models native vertical coordinate.
!
! For free blobs, it searches for the vertical grid cell that a free
! blob resides in.
!
! We require the information regarding the vertical grid cell that a
! blob resides in prior to the calculation of total grid cell thickness
! and prior to the calculation of the vertical advection velocity.
!
! In order to figure out which grid cell the blob is in, we need to
! employ different strategies for different coordinate system (see notes
! for details).
!
! For any blobs that penetrate the free surface (which is unlikely, but
! not impossible) we immediately kill them, returning their properties
! to the surface grid cell of the (i,j) column that they belong to.
!
! We recalculate a blobs density and volume for taup1, based on the pressure
! at the centre of the grid cell that it resides in.
!
! </DESCRIPTION>
!
subroutine ocean_blob_cell_update(Time, Thickness, Dens, T_prog, Ext_mode, L_system, EL_diag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  type(blob_diag_type),           intent(inout) :: EL_diag(0:)

  type(ocean_blob_type), pointer :: this=>NULL()
  integer :: i,j,k,n,old_k,km1,tau,taup1
  integer :: stdoutunit,thru_surface
  real :: zwt_k, zwt_km1, dz, p

  tau   = Time%tau
  taup1 = Time%taup1

  thru_surface = 0

  this=>head_dynamic_bott
  if(associated(this)) then
     bottcycle: do
        i = this%i
        j = this%j
        k = this%k
        if (vert_coordinate==ZSTAR) then
           ! From table 6.1 of Griffies (2009): zstar = H(z-eta)/(H+eta)
           ! Here z=-this%geodepth
           this%st = Grd%ht(i,j) * ( -this%geodepth-Ext_mode%eta_t(i,j,taup1) )&
                                  / ( Grd%ht(i,j)+Ext_mode%eta_t(i,j,taup1))

        elseif(vert_coordinate==PRESSURE) then
           ! From table 6.2 of Griffies (2009): p
           ! Here p = p(k) + dz/dzt_dst
           ! where dz is the difference in depth between the blob and the grid cell centre
           dz = this%geodepth - Thickness%depth_zt(i,j,k)
           this%st = Dens%pressure_at_depth(i,j,k) + dz/Thickness%dzt_dst(i,j,k)

        elseif(vert_coordinate==PSTAR) then
           ! From table 6.2 of Griffies (2009): pstar = pb0(p-p_a)/(p_b-p_a)
           ! Here p = p(k) + dz/dzt_dst
           ! where dz is the difference in depth between the blob and the grid cell centre
           dz = this%geodepth - Thickness%depth_zt(i,j,k)
           p  = Dens%pressure_at_depth(i,j,k) + dz/Thickness%dzt_dst(i,j,k)
           this%st = Thickness%pbot0(i,j) * ( p - Ext_mode%patm_t(i,j,taup1) )  &
                                           /( Ext_mode%pbot_t(i,j,taup1) - Ext_mode%patm_t(i,j,taup1) )

        endif

        this => this%next
        if (.not. associated(this)) exit bottcycle
     enddo bottcycle
  endif


  this=>head_dynamic_free
  if(associated(this)) then
     blobcycle: do

        i = this%i
        j = this%j
        old_k = this%k
        ! We calculate the depth using eta_t(taup1) in ocean_blob_diagnose_depth.
        ! It is only done here for debugging purposes (to make sure that intermediate
        ! checksums remain consistent).
        this%depth = this%geodepth + Ext_mode%eta_t(i,j,tau)

        if (vert_coordinate==ZSTAR) then
           ! From table 6.1 of Griffies (2009): zstar = H(z-eta)/(H+eta)
           ! Here z=-this%geodepth
           this%st = Grd%ht(i,j) * ( -this%geodepth-Ext_mode%eta_t(i,j,taup1) )&
                                  / ( Grd%ht(i,j)+Ext_mode%eta_t(i,j,taup1))

           ! Note: -H<=zstar<=0 (table 6.1 of Griffies, 2009) and st=zstar
           k=1
           if (this%st>0.) then
              ! the blob has penetrated the free surface, so kill the blob
              thru_surface=thru_surface+1
              EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + this%mass
              do n=1,num_prog_tracers
                 EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + this%tracer(n)
              enddo
              call kill_blob(Thickness, T_prog, this, i, j, k)
              this=>this%next
              if(.not. associated(this)) exit blobcycle
              cycle blobcycle
           else ! the blob is at or below the surface.
              call find_depth()
           endif

        else !vert_coordinate==PSTAR .or. PRESSURE
           ! Note: p_a<=p<=p_b; 0<=pstar<=p_b^0
           k=Grd%kmt(i,j)
           zwt_k = -Grd%ht(i,j)
           if (-this%geodepth<zwt_k) then
              ! the blob has penetrated the bottom surface.  Raise a warning.
              this%geodepth = -zwt_k
           else
              press_cycle: do
                 km1=k-1
                 zwt_km1 = zwt_k + Thickness%dzt_dst(i,j,k)*Thickness%dst(i,j,k)
                 if (zwt_k <= -this%geodepth .and. -this%geodepth < zwt_km1) then
                    exit press_cycle
                 endif
                 k=km1
                 zwt_k=zwt_km1
                 if (k==0) then
                    ! the blob has penetrated the free surface, so kill the blob
                    thru_surface=thru_surface+1
                    EL_diag(0)%dstry(i,j,1) = EL_diag(0)%dstry(i,j,1) + this%mass
                    do n=1,num_prog_tracers
                       EL_diag(n)%dstry(i,j,1) = EL_diag(n)%dstry(i,j,1) + this%tracer(n)
                    enddo
                    call kill_blob(Thickness, T_prog, this, i, j, 1)
                    this=>this%next
                    if(.not. associated(this)) exit blobcycle
                    cycle blobcycle
                 endif
              enddo press_cycle
           endif
           ! Calculate the st of the blob for later use (in update_L_thickness and adjust_L_thickness)
           this%st = Thickness%depth_swt(i,j,k) - (zwt_k + this%geodepth)/Thickness%dzt_dst(i,j,k)
        endif

        this%k = k

        ! If the blob is in fact in a different cell, we take this into account with the
        ! calculate of divergence.
        if (old_k/=k) then
           L_system%conv_blob(i,j,old_k) = L_system%conv_blob(i,j,old_k) - this%mass
           L_system%conv_blob(i,j,k)     = L_system%conv_blob(i,j,k)     + this%mass

           ! Diagnostics
           mass_in( i,j,k)     = mass_in( i,j,k)     + this%mass
           mass_out(i,j,old_k) = mass_out(i,j,old_k) + this%mass
        endif

        if (this%mass>blob_small_mass) then
           ! Update some of the blob variables if they wont give silly numbers
           ! Blobs with mass<blob_small_mass may be required only for their history, so
           ! updating their properties is not important as they will be erased from
           ! memory soon enough.
           do n=1,num_prog_tracers
              this%field(n) = this%tracer(n)/this%mass
           enddo
           this%density  = density(this%field(index_salt), &
                                   this%field(index_temp), &
                                   Dens%pressure_at_depth(i,j,k))
           this%densityr = 1./this%density
           if (vert_coordinate_class == DEPTH_BASED) then
              this%volume = this%mass * rho0r
           else
              this%volume = this%mass * this%densityr
           endif
        endif

        this => this%next
        if (.not. associated(this)) exit blobcycle
     enddo blobcycle
  endif !associated(this)

  ! We may now complete the calculation of blob divergence.
  do k=1,nk
     L_system%conv_blob(:,:,k) = L_system%conv_blob(:,:,k)*datdtimer(:,:)
  enddo
  call mpp_update_domains(Thickness%blob_source(:,:), Dom%domain2d)
  call mpp_update_domains(L_system%conv_blob(:,:,:),  Dom%domain2d, complete=.false.)
  do n=1,num_prog_tracers
     call mpp_update_domains(T_prog(n)%tend_blob(:,:,:), Dom%domain2d, complete=T_prog(n)%complete)
  enddo

  ! Diagnostics
  call diagnose_3d(Time, Grd, id_mass_in,  mass_in(:,:,:))
  call diagnose_3d(Time, Grd, id_mass_out, mass_out(:,:,:))

  if (debug_this_module) then
     stdoutunit = stdout()
     call mpp_sum(thru_surface)
     write(stdoutunit, *) 'Free blobs destroyed from reaching surface =', thru_surface
     write(stdoutunit, '(a)') ' '
     call write_timestamp(Time%model_time)
     write(stdoutunit, *) 'Intermediate Blob Checksums (taup1) after explicit update'
     call blob_chksum(T_prog(:), head_static_free, head_static_bott, &
                      head_dynamic_free, head_dynamic_bott, blob_counter)
     call write_chksum_2d('Ext_mode%conv_blob', Ext_mode%conv_blob(COMP))
     call write_chksum_3d('L_system%conv_blob', L_system%conv_blob(COMP,:))
     if (vert_coordinate_class==DEPTH_BASED) then
        call write_chksum_2d('Ext_mode%eta_t(taup1)', Ext_mode%eta_t(COMP,taup1))
     else
        call write_chksum_2d('Ext_mode%eta_t(tau)', Ext_mode%eta_t(COMP,tau))
        call write_chksum_3d('Thickness%dst', Thickness%dzt_dst(COMP,:))
        call write_chksum_3d('Thickness%dzt_dst', Thickness%dzt_dst(COMP,:))
     endif
     call write_chksum_2d('Thickness%blob_source', Thickness%blob_source(COMP))
     call write_chksum_3d('T_prog(temp)%tend_blob', T_prog(index_temp)%tend_blob(COMP,:))
     call write_chksum_3d('T_prog(salt)%tend_blob', T_prog(index_salt)%tend_blob(COMP,:))

     if (really_debug) call write_all_blobs('tau')
  endif

  call blob_delete(Time, Thickness, T_prog(:), head_dynamic_free)

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that finds the k cell that a blob resides!
! in.  If a blob penetrates the bottom boundary, we raise an error.    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine find_depth()
    if (this%st<-Thickness%depth_swt(i,j,1)) then
       k=2
       depth_cycle: do
          if ( -Thickness%depth_swt(i,j,k)<=this%st .and. this%st<-Thickness%depth_swt(i,j,k-1)) then
             exit depth_cycle
          endif
          k=k+1
       enddo depth_cycle
    endif

  end subroutine find_depth

end subroutine ocean_blob_cell_update
! </SUBROUTINE>  NAME="ocean_blob_cell_update"

!#######################################################################
! <SUBROUTINE NAME="update_L_thickness">
!
! <DESCRIPTION>
! Calculates the contribution to thickness of all the blobs for the
! L system arrays in the Thickness strcture.
!
! The mass per unit area from the blobs is also calculated for the
! upper and lower part of each grid cell and stored in the L_system
! structure.  The calculate of mass per unit area is required for
! pressure calculations.
!
! We note that in DEPTH_BASED models, the value for density used is
! rho0, while it is the actual density in PRESSURE_BASED models. For
! the L_system mass per unit area, the actual density is used for
! all supported coordinates.
! </DESCRIPTION>
!
subroutine update_L_thickness(Time, Thickness, T_prog, L_system, EL_diag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  type(blob_diag_type),           intent(inout) :: EL_diag(0:)

  type(ocean_blob_type), pointer :: this=>NULL()
  integer :: i,j,k,n
  integer :: tau,taup1
  integer :: total_blobs, dstryd_blobs, diag_dstryd_blobs
  integer :: stdoutunit
  logical :: up
  real :: dzt, rho_dzt, th_rho_dzt, max_dzt, max_dzt_lo, max_dzt_up 
  real :: new_dzt, new_dzt_lo, new_dzt_up

  if (.not. use_this_module) return

  tau   = Time%tau
  taup1 = Time%taup1

  ! We find the blob's position in coordinate space and then use that to
  ! find which vertical grid cell it resides in at taup1.

  do n=1,num_prog_tracers
     T_prog(n)%sum_blob(:,:,:,taup1) = 0.0
  enddo
  Thickness%rho_dztL(:,:,:,taup1) = 0.0
  Thickness%dztupL(:,:,:) = 0.0
  Thickness%dztloL(:,:,:) = 0.0
  L_system%rho_dztup(:,:,:) = 0.0
  L_system%rho_dztlo(:,:,:) = 0.0

  ! We have two kinds of rho_dzt that we save.  One is for the pressure calculation,
  ! the other is for the thickness.  For the pressure calculation, we need to use
  ! the actual density, regardless of whether it is DEPTH_BASED or PRESSURE_BASED.
  ! The second type is for thickness calculations.  For calculation of thickness, it
  ! does depend on the type of vertical coordinate as to whether we use the reference
  ! density, or the actual density.

  this=>head_dynamic_free
  total_blobs  = 0
  dstryd_blobs = 0

  if(associated(this)) then
     freecycle: do
        total_blobs = total_blobs + 1

        i = this%i - Dom%ioff
        j = this%j - Dom%joff
        k = this%k

        ! Calculate this blobs contribution to thickness and mass per unit area
        dzt     = Grd%datr(i,j)*this%volume
        rho_dzt = dzt*this%density

        if (vert_coordinate_class==DEPTH_BASED) then
           ! For the Thickness structure, we use a reference density
           th_rho_dzt = dzt*rho0
           if (this%st>-Thickness%depth_st(i,j,k)) then
              up=.true.
           else
              up=.false.
           endif

        else !vert_coordinate_class==PRESSURE_BASED
           th_rho_dzt = rho_dzt
           if (this%st<Thickness%depth_st(i,j,k)) then
              up=.true.
           else
              up=.false.
           endif

        endif

        ! If the contribution to thickness from all the blobs in a grid cell exceeds
        ! the total thickness of the cell, we kill the blob and return all its
        ! properties to the E system so that we do not have a negative E system
        ! thickness.
        if (up) then
           new_dzt = Thickness%dztupL(i,j,k)+dzt
           max_dzt = max_prop_thickness*Thickness%dztupT(i,j,k)
           if (new_dzt > max_dzt) then
              EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + this%mass
              do n=1,num_prog_tracers
                 EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + this%tracer(n)
              enddo
              call kill_blob(Thickness,T_prog(:),this,i,j,k)
              dstryd_blobs = dstryd_blobs + 1
              this=>this%next
              if(.not.associated(this)) exit freecycle
              cycle freecycle
           endif
           Thickness%dztupL(i,j,k)   = new_dzt
           L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt
        else
           new_dzt = Thickness%dztloL(i,j,k)+dzt
           max_dzt = max_prop_thickness*Thickness%dztloT(i,j,k)
           if (new_dzt > max_dzt) then
              EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + this%mass
              do n=1,num_prog_tracers
                 EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + this%tracer(n)
              enddo
              call kill_blob(Thickness,T_prog(:),this,i,j,k)
              dstryd_blobs = dstryd_blobs + 1
              this=>this%next
              if(.not.associated(this)) exit freecycle
              cycle freecycle
           endif
           Thickness%dztloL(i,j,k)   = new_dzt
           L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt
        endif

        Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + th_rho_dzt

        ! Accumulate the blob tracer content in a cell
        do n=1,num_prog_tracers
           T_prog(n)%sum_blob(i,j,k,taup1) = T_prog(n)%sum_blob(i,j,k,taup1) + this%tracer(n)
        enddo

        this=>this%next
        if(.not.associated(this)) exit freecycle
     enddo freecycle
  endif
  nullify(this)

  diag_dstryd_blobs = dstryd_blobs

  if (debug_this_module) then
     stdoutunit = stdout()
     call mpp_sum(total_blobs); call mpp_sum(dstryd_blobs)
     write(stdoutunit, *) ' '
     write(stdoutunit, *) 'Statistics from ocean_blob_mod: update_L_thickness'
     write(stdoutunit, *) 'Total dynamic free blobs          =', total_blobs
     write(stdoutunit, *) 'Blobs destroyed because dztL>dztT =', dstryd_blobs
  endif

  total_blobs  = 0
  dstryd_blobs = 0

  this=>head_dynamic_bott
  if (associated(this)) then
     bottcycle: do
        total_blobs = total_blobs + 1

        i = this%i - Dom%ioff
        j = this%j - Dom%joff
        k = this%k

        ! Calculate this blobs contribution to thickness and mass per unit area
        dzt     = Grd%datr(i,j)*this%volume
        rho_dzt = dzt*this%density

        if (vert_coordinate_class==DEPTH_BASED) then
           ! For the Thickness structure, we use a reference density
           th_rho_dzt = dzt*rho0
        else !vert_coordinate_class==PRESSURE_BASED
           th_rho_dzt = rho_dzt
        endif

        ! Bottom blobs now contribute to the lower and upper half of the cell
        new_dzt_lo = Thickness%dztloL(i,j,k)+dzt*0.5
        new_dzt_up = Thickness%dztupL(i,j,k)+dzt*0.5
        max_dzt_lo = max_prop_thickness*Thickness%dztloT(i,j,k)
        max_dzt_up = max_prop_thickness*Thickness%dztupT(i,j,k)

        if ((new_dzt_lo > max_dzt_lo) .or. (new_dzt_up > max_dzt_up)) then
           EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + this%mass
           do n=1,num_prog_tracers
              EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + this%tracer(n)
           enddo
           call kill_blob(Thickness,T_prog(:),this,i,j,k)
           dstryd_blobs = dstryd_blobs + 1
           this=>this%next
           if(.not.associated(this)) exit bottcycle
           cycle bottcycle
        endif
        Thickness%dztloL(i,j,k)         = new_dzt_lo
        Thickness%dztupL(i,j,k)         = new_dzt_up
        L_system%rho_dztlo(i,j,k)       = L_system%rho_dztlo(i,j,k) + rho_dzt*0.5
        L_system%rho_dztup(i,j,k)       = L_system%rho_dztup(i,j,k) + rho_dzt*0.5
        Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + th_rho_dzt

        do n=1,num_prog_tracers
           T_prog(n)%sum_blob(i,j,k,taup1) = T_prog(n)%sum_blob(i,j,k,taup1) + this%tracer(n)
        enddo

        this=>this%next
        if(.not.associated(this)) exit bottcycle
     enddo bottcycle
  endif

  diag_dstryd_blobs = diag_dstryd_blobs + dstryd_blobs
  if (debug_this_module) then
     stdoutunit = stdout()
     call mpp_sum(total_blobs); call mpp_sum(dstryd_blobs)
     write(stdoutunit, *) 'Total dynamic bottom blobs        =', total_blobs
     write(stdoutunit, *) 'Blobs destroyed because dztL>dztT =', dstryd_blobs
  endif

  call mpp_update_domains(L_system%rho_dztup(:,:,:),       Dom%domain2d, complete=.false.)
  call mpp_update_domains(L_system%rho_dztlo(:,:,:),       Dom%domain2d, complete=.false.)
  call mpp_update_domains(Thickness%rho_dztL(:,:,:,taup1), Dom%domain2d, complete=.false.)
  call mpp_update_domains(Thickness%dztupL(:,:,:),         Dom%domain2d, complete=.false.)
  call mpp_update_domains(Thickness%dztloL(:,:,:),         Dom%domain2d, complete=.false.)
  do n=1,num_prog_tracers
     call mpp_update_domains(T_prog(n)%sum_blob(:,:,:,taup1), Dom%domain2d, complete=T_prog(n)%complete)
  enddo

  Thickness%dztL(:,:,:) = Thickness%dztloL(:,:,:) + Thickness%dztupL(:,:,:)

  Thickness%dzwtL(:,:,0)      = Thickness%dztupL(:,:,1)
  Thickness%dzwtL(:,:,1:nk-1) = Thickness%dztupL(:,:,2:nk) + Thickness%dztloL(:,:,1:nk-1)
  Thickness%dzwtL(:,:,nk)     = Thickness%dztloL(:,:,nk)

  call blob_delete(Time, Thickness, T_prog(:),head_static_free)
  call blob_delete(Time, Thickness, T_prog(:),head_static_bott)
  call blob_delete(Time, Thickness, T_prog(:),head_dynamic_free)
  call blob_delete(Time, Thickness, T_prog(:),head_dynamic_bott)

  call mpp_sum(diag_dstryd_blobs)
  if (id_dstryd>0) used = send_data(id_dstryd, real(diag_dstryd_blobs), Time%model_time)

end subroutine update_L_thickness
! </SUBROUTINE>  NAME="update_L_thickness"


!#######################################################################
! <SUBROUTINE NAME="calculate_rhoT">
!
! <DESCRIPTION>
! Calculates the density of the combined E and L systems.
!
! It needs to be done here so that the blobs are taken into account.
!
! </DESCRIPTION>
!
subroutine calculate_rhoT(Time, Dens, Thickness)
  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_density_type),   intent(inout) :: Dens
  type(ocean_thickness_type), intent(in)    :: Thickness

  real, dimension(isd:ied,jsd:jed,nk) :: rho_dzt
  type(ocean_blob_type), pointer :: this
  integer :: taup1,i,j,k,p

  if (.not. use_this_module) return

  taup1 = Time%taup1

  rho_dzt(:,:,:) = Dens%rho(:,:,:,taup1)*Thickness%dzt(:,:,:)

  do p=3,4
     if (p==3) this=>head_dynamic_free
     if (p==4) this=>head_dynamic_bott

     do while (associated(this))
        i=this%i
        j=this%j
        k=this%k

        rho_dzt(i,j,k) = rho_dzt(i,j,k) + this%mass*Grd%datr(i,j)

        this=>this%next
     enddo
  enddo

  call mpp_update_domains(rho_dzt(:,:,:), Dom%domain2d)

  ! Update rhoT to time taup1
  Dens%rhoT(:,:,:) = rho_dzt(:,:,:)/Thickness%dztT(:,:,:,taup1)

end subroutine calculate_rhoT
! </SUBROUTINE>  NAME="calculate_rhoT"


!#######################################################################
! <SUBROUTINE NAME="ocean_blob_implicit">
!
! <DESCRIPTION>
! Updates the Lagrangian blobs by calling routines that run during the
! implicit in time part of the time step.
! </DESCRIPTION>
!
subroutine ocean_blob_implicit(Time, Thickness, T_prog, Dens, Adv_vel, &
                               Velocity, EL_diag)
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_velocity_type),    intent(in)    :: Velocity
  type(blob_diag_type),         intent(inout) :: EL_diag(0:)

  integer :: n,taup1

  taup1 = Time%taup1

  if (.not. use_this_module) return
  if (.not. module_is_initialized ) then
     call mpp_error(FATAL, &
          '==>Error in ocean_blob_mod (Lagrangian blob model): module '&
          //' needs to be initialized')
  endif

  ! run the static open ocean convection scheme
  call mpp_clock_begin(id_clock_static_free)
  call blob_static_free(Time, Thickness, T_prog, Dens, head_static_free)
  call mpp_clock_end(id_clock_static_free)

  ! create new blobs arising from vertical instabilities
  call mpp_clock_begin(id_clock_dyn_free_update)
  if (use_dyn_fre) then
     call blob_dynamic_free_implicit(Time, Thickness, T_prog(:), Dens, Adv_vel, &
                                     Velocity, head_dynamic_free, blob_counter, &
                                     EL_diag(:))
  endif
  call mpp_clock_end(id_clock_dyn_free_update)

  if (id_prop_cell_mass>0) then
     wrk1(:,:,:) = Thickness%rho_dztL(:,:,:,taup1)/Thickness%rho_dztT(:,:,:,taup1)
     call diagnose_3d(Time, Grd, id_prop_cell_mass, wrk1(:,:,:))
  endif

  do n=0,num_prog_tracers
     call diagnose_3d(Time, Grd, id_entrainment(n), EL_diag(n)%entrainment(:,:,:))
     call diagnose_3d(Time, Grd, id_detrainment(n), EL_diag(n)%detrainment(:,:,:))
     call diagnose_3d(Time, Grd, id_new(n), EL_diag(n)%new(:,:,:))
     call diagnose_3d(Time, Grd, id_dstry(n), EL_diag(n)%dstry(:,:,:))
  enddo

end subroutine ocean_blob_implicit
! </SUBROUTINE>  NAME="ocean_blob_implicit"

!######################################################################
! <SUBROUTINE NAME="adjust_L_thickness">
!
! <DESCRIPTION>
! blob_thickness is called after new blobs are formed implicitly in
! time.  blob_thickness provides the same function as
! update_L_thickness, with a couple of subtle differences.
!
! At the point that this routine is called, we:
! 1/ Know the vertical position of an "old" blob relative to the geoid
!    and in the Eulerian model's native vertical coordinate.  We do
!    not know a blobs depth relative to the sea surface.
! 2/ Know the vertical position relative to the sea surface for "new"
!    blobs, but we do not know its position relative to the geoid or in
!    the Euerlian model's native vertical coordinate.
!
! So, blob_thickness uses the logical new, which is part of the blob
! derived type, to distinguish between new blobs and old blobs. New
! blobs have their position (upper or lower part of a grid cell)
! calculated using the depth relative to the sea surface height, while
! blobs that are not new have their position calculated using the
! Eulerian models native coordinate.
!
! The accumulation of thickness and mass per unit area are done in the
! same was as is done in update_L_thickness.
! </DESCRIPTION>
!
subroutine adjust_L_thickness(Time, Thickness, T_prog, L_system)

  type(ocean_time_type),                  intent(in)    :: Time
  type(ocean_thickness_type),             intent(inout) :: Thickness
  type(ocean_prog_tracer_type),           intent(inout) :: T_prog(:)
  type(ocean_lagrangian_type),            intent(inout) :: L_system

  type(ocean_blob_type), pointer :: this=>NULL()
  character(len=128) :: tname
  integer :: i,j,k,n
  integer :: taup1
  integer :: stdoutunit
  real    :: rho_dzt

  if (.not. use_this_module) return

  taup1 = Time%taup1

  Thickness%dztupL(:,:,:) = 0.0
  Thickness%dztloL(:,:,:) = 0.0
  do n=1,num_prog_tracers
     T_prog(n)%sum_blob(:,:,:,taup1) = 0.0
  enddo
  Thickness%rho_dztL(:,:,:,taup1) = 0.0
  L_system%rho_dztup(:,:,:) = 0.0
  L_system%rho_dztlo(:,:,:) = 0.0

  this=>head_dynamic_free
  if(associated(this)) then
     freecycle: do
        ! for convenience
        i  = this%i - Dom%ioff
        j  = this%j - Dom%joff
        k  = this%k

        ! We have two kinds of rho_dzt that we save.  One is for the pressure calculation,
        ! the other is for the thickness.  For the pressure calculation, we need to use
        ! the actual density, regardless of whether it is DEPTH_BASED or PRESSURE_BASED.
        ! The second type is for thickness calculations.  For calculation of thickness, it
        ! does depend on the type of vertical coordinate as to whether we use the reference
        ! density, or the actual density.

        ! This is for pressure calculations
        rho_dzt = Grd%datr(i,j)*this%volume*this%density

        ! For new blobs, we do not know what their vertical position in
        ! coordiante space is, but, we do know their depth.  For old blobs
        ! we do not know their depth, but we know their position in coordinate
        ! space.
        if (this%new) then
           if (Thickness%depth_zt(i,j,k)<=this%depth) then
              Thickness%dztloL(i,j,k)   = Thickness%dztloL(i,j,k)   + Grd%datr(i,j)*this%volume
              L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt
           else
              Thickness%dztupL(i,j,k)   = Thickness%dztupL(i,j,k)   + Grd%datr(i,j)*this%volume
              L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt
           endif
           ! This is for thickness calculations
           if(vert_coordinate_class == DEPTH_BASED)then
              Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + Grd%datr(i,j)*this%volume*rho0
           else
              Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + rho_dzt
           endif
        else!not a new blob
           if (vert_coordinate_class==DEPTH_BASED) then
              if (this%st>-Thickness%depth_st(i,j,k)) then
                 Thickness%dztupL(i,j,k)   = Thickness%dztupL(i,j,k)   + Grd%datr(i,j)*this%volume
                 L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt
              else
                 Thickness%dztloL(i,j,k)   = Thickness%dztloL(i,j,k)   + Grd%datr(i,j)*this%volume
                 L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt
              endif
              Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + Grd%datr(i,j)*this%volume*rho0
           else !vert_coordinate_class==PRESSURE_BASED
              if (this%st<Thickness%depth_st(i,j,k)) then
                 Thickness%dztupL(i,j,k)   = Thickness%dztupL(i,j,k)   + Grd%datr(i,j)*this%volume
                 L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt
              else
                 Thickness%dztloL(i,j,k)   = Thickness%dztloL(i,j,k)   + Grd%datr(i,j)*this%volume
                 L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt
              endif
              Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + rho_dzt
           endif
        endif !new blob

        do n=1,num_prog_tracers
           T_prog(n)%sum_blob(i,j,k,taup1) = T_prog(n)%sum_blob(i,j,k,taup1) + this%tracer(n)
        enddo

        this=>this%next
        if(.not.associated(this)) exit freecycle
     enddo freecycle
  endif

  this=>head_dynamic_bott
  if (associated(this)) then
     bottcycle: do
        ! for convenience
        i  = this%i - Dom%ioff
        j  = this%j - Dom%joff
        k  = this%k

        ! This is for pressure calculations
        rho_dzt = Grd%datr(i,j)*this%volume*this%density

        ! Bottom blobs only contribute to the lower half of the bottom cell.
        L_system%rho_dztlo(i,j,k) = L_system%rho_dztlo(i,j,k) + rho_dzt*0.5
        L_system%rho_dztup(i,j,k) = L_system%rho_dztup(i,j,k) + rho_dzt*0.5

        ! This is for thickness calculations
        if(vert_coordinate_class == DEPTH_BASED)then
           Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + Grd%datr(i,j)*this%volume*rho0
        else
           Thickness%rho_dztL(i,j,k,taup1) = Thickness%rho_dztL(i,j,k,taup1) + rho_dzt
        endif
        Thickness%dztloL(i,j,k) = Thickness%dztloL(i,j,k) + Grd%datr(i,j)*this%volume*0.5
        Thickness%dztupL(i,j,k) = Thickness%dztupL(i,j,k) + Grd%datr(i,j)*this%volume*0.5

        do n=1,num_prog_tracers
           T_prog(n)%sum_blob(i,j,k,taup1) = T_prog(n)%sum_blob(i,j,k,taup1) + this%tracer(n)
        enddo

        this=>this%next
        if(.not.associated(this)) exit bottcycle
     enddo bottcycle
  endif

  call mpp_update_domains(L_system%rho_dztup(:,:,:),       Dom%domain2d, complete=.false.)
  call mpp_update_domains(L_system%rho_dztlo(:,:,:),       Dom%domain2d, complete=.false.)
  call mpp_update_domains(Thickness%rho_dztL(:,:,:,taup1), Dom%domain2d, complete=.false.)
  call mpp_update_domains(Thickness%dztupL(:,:,:),         Dom%domain2d, complete=.false.)
  call mpp_update_domains(Thickness%dztloL(:,:,:),         Dom%domain2d, complete=.false.)
  do n=1,num_prog_tracers
     call mpp_update_domains(T_prog(n)%sum_blob(:,:,:,taup1), Dom%domain2d, complete=T_prog(n)%complete)
  enddo

  Thickness%dztL(:,:,:)       = Thickness%dztloL(:,:,:) + Thickness%dztupL(:,:,:)
  Thickness%dzwtL(:,:,0)      = Thickness%dztupL(:,:,1)
  Thickness%dzwtL(:,:,1:nk-1) = Thickness%dztupL(:,:,2:nk) + Thickness%dztloL(:,:,1:nk-1)
  Thickness%dzwtL(:,:,nk)     = Thickness%dztloL(:,:,nk)

  if (debug_this_module) then
     stdoutunit = stdout()
     write(stdoutunit, '(/,a)')  &
        'Lagrangian Tracer and Mass Checksums (taup1)'
     call write_timestamp(Time%model_time)
     call write_chksum_3d('blob convergence', L_system%conv_blob(COMP,:))
     do n=1,num_prog_tracers
        tname = trim(T_prog(n)%name)//' tend'
        call write_chksum_3d(tname(1:17), T_prog(n)%tend_blob(COMP,:))
     enddo
     call write_chksum_2d('mass', Grd%dat(COMP)*sum(Thickness%rho_dztL(COMP,:,taup1),3))
     do n=1,num_prog_tracers
        tname = trim(T_prog(n)%name)
        call write_chksum_3d(tname(1:18), T_prog(n)%sum_blob(COMP,:,taup1))
     enddo
     write(stdoutunit, '(a,/)') 'end Lagrangian Tracer and Mass Checksums'
  endif

end subroutine adjust_L_thickness
! </SUBROUTINE>  NAME="adjust_L_thickness"


!#######################################################################
! <SUBROUTINE NAME="ocean_blob_diagnose_depth">
! <DESCRIPTION>
! When blobs are created implicitly in time, the sea surface height has
! not been calcualted in pressure based models (but depth has).  However,
! the prognostic variable for blobs is the geodepth and not the depth.  So,
! for a new blob, we set the depth when it is formed and then calculate
! the geodepth here.  For existing blobs, we diagnose the depth using
! eta_t and geodepth.
!
! We note that in GEOPOTENTIAL coordinates, depth and geodepth are
! equivalent, and so we set them to be equivalent for the blobs too.
! </DESCRIPTION>
!
subroutine ocean_blob_diagnose_depth(Time, T_prog, Ext_mode, L_system)
  type(ocean_prog_tracer_type),    intent(in) :: T_prog(:)
  type(ocean_time_type),           intent(in) :: Time
  type(ocean_external_mode_type),  intent(in) :: Ext_mode
  type(ocean_lagrangian_type),     intent(in) :: L_system

  type(ocean_blob_type), pointer :: this => NULL()
  integer :: i,j,k,p
  integer :: taup1
  integer :: stdoutunit

  if (.not. use_this_module) return
  stdoutunit = stdout()

  taup1 = Time%taup1
  do p=1,4
     if(p==1) this=>head_static_free
     if(p==2) this=>head_static_bott
     if(p==3) this=>head_dynamic_free
     if(p==4) this=>head_dynamic_bott
     if(associated(this)) then
        geoblobcycle: do
           i = this%i - Dom%ioff
           j = this%j - Dom%joff
           k = this%k
           if (this%new) then
              this%geodepth = this%depth - Ext_mode%eta_t(i,j,taup1)
              this%new = .false.
           else
              this%depth = this%geodepth + Ext_mode%eta_t(i,j,taup1)
           endif
           this=>this%next
           if(.not.associated(this)) exit geoblobcycle
        enddo geoblobcycle
     endif
     nullify(this)
  enddo

  call blob_diag(Time, head_static_free,  T_prog(:), 1)
  call blob_diag(Time, head_static_bott,  T_prog(:), 2)
  call blob_diag(Time, head_dynamic_free, T_prog(:), 3)
  call blob_diag(Time, head_dynamic_bott, T_prog(:), 4)

  if (debug_this_module) then
     call write_timestamp(Time%model_time)
     write(stdoutunit, '(/,a)') 'Blob Checksums (taup1)'
     call blob_chksum(T_prog(:), head_static_free, head_static_bott, &
                      head_dynamic_free, head_dynamic_bott, blob_counter)
     write(stdoutunit, '(/,a)') 'Lagrangian System Checksums (taup1)'
     call lagrangian_system_chksum(L_system)

     if (really_debug) call write_all_blobs('taup1')
  endif

end subroutine ocean_blob_diagnose_depth
! </SUBROUTINE>  NAME="ocean_blob_diagnose_depth"

!######################################################################
! <SUBROUTINE NAME="ocean_blob_end">
!
! <DESCRIPTION>
! Writes restarts and do checksums for the blobs at the end of a run.
! </DESCRIPTION>
!
subroutine ocean_blob_end(Time, T_prog, L_system, time_stamp)

  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_prog_tracer_type),   intent(in)           :: T_prog(:)
  type(ocean_lagrangian_type),    intent(inout)        :: L_system
  character(len=*),               intent(in), optional :: time_stamp

  integer :: taup1
  integer :: stdoutunit,stderrunit

  if (.not. use_this_module) return
  stdoutunit = stdout() ; stderrunit = stderr()
  taup1 = Time%taup1

  write(stdoutunit,'(/a)') 'From ocean_blobs_mod: ending blob chksums'
  call write_timestamp(Time%model_time)
  call blob_chksum(T_prog(:), head_static_free, head_static_bott, &
                   head_dynamic_free, head_dynamic_bott, blob_counter)

  call writerestart

  if (really_debug) call write_all_blobs('taup1')

  call deallocateblobs

  call blob_static_free_end()
  call blob_static_bottom_end()
  call blob_dynamic_free_end()
  call blob_dynamic_bottom_end()
  call blob_util_end()
  call blob_diag_end(T_prog(:))

  nullify(Info)
  nullify(Dom)
  nullify(Grd)

  call write_timestamp(Time%model_time)
  call lagrangian_system_chksum(L_system)

  write(stdoutunit,'(/a)') 'Completed end of Lagrangian Blobs'

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that write a blob restart at the end of
! a model run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   subroutine writerestart

     type(ocean_blob_type), pointer :: this=>NULL()
     integer, allocatable, dimension(:) :: tracerid
     integer           :: latid, lonid, geodepthid, massid, uid, vid, wid, stid
     integer           :: hashid, numberid, typeid, iid, jid, kid
     integer           :: stepid, h1id, h2id, nstepsid, mstepsid
     integer           :: densityid, entid, dragid, ageid, heightid
     integer           :: mret, ncid, m_dim, n, m, p
     character(len=31) :: filename

     ! There are two restart files for the blobs.  1/ for gridded data required
     ! by the blobs (e.g. blob_counter) 2/ for the blobs themselves

     ! first, the gridded file
     call save_restart(gridded_restart, time_stamp)

     ! second, the blobs themselves
     write(filename(1:31),'("RESTART/ocean_blobs.res.nc.", I4.4)') mpp_pe()
     if (debug_this_module) write(stdoutunit, '(2a)') 'creating blob restart file:', filename

     mret = nf_create(filename, NF_CLOBBER, ncid)
     if (mret .ne. NF_NOERR) write(stderrunit,'(a)') 'blobs, writerestart: nf_create failed'

     ! Register dimensions
     mret = nf_def_dim(ncid, 'blob', NF_UNLIMITED, m_dim)
     if (mret .ne. NF_NOERR) write(stderrunit,'(a)') 'blobs, writerestart: nf_def_dim failed'

     ! Register variables
     latid       = def_var(ncid, 'lat',         NF_DOUBLE, m_dim)
     lonid       = def_var(ncid, 'lon',         NF_DOUBLE, m_dim)
     geodepthid  = def_var(ncid, 'geodepth',    NF_DOUBLE, m_dim)
     stid        = def_var(ncid, 'st',          NF_DOUBLE, m_dim)
     massid      = def_var(ncid, 'mass',        NF_DOUBLE, m_dim)
     densityid   = def_var(ncid, 'density',     NF_DOUBLE, m_dim)
     uid         = def_var(ncid, 'u',           NF_DOUBLE, m_dim)
     vid         = def_var(ncid, 'v',           NF_DOUBLE, m_dim)
     wid         = def_var(ncid, 'w',           NF_DOUBLE, m_dim)
     entid       = def_var(ncid, 'ent',         NF_DOUBLE, m_dim)
     dragid      = def_var(ncid, 'drag',        NF_DOUBLE, m_dim)
     ageid       = def_var(ncid, 'age',         NF_DOUBLE, m_dim)
     heightid    = def_var(ncid, 'height',      NF_DOUBLE, m_dim)
     h1id        = def_var(ncid, 'h1',          NF_DOUBLE, m_dim)
     h2id        = def_var(ncid, 'h2',          NF_DOUBLE, m_dim)
     stepid      = def_var(ncid, 'step',        NF_DOUBLE, m_dim)
     hashid      = def_var(ncid, 'hash',        NF_INT,    m_dim)
     numberid    = def_var(ncid, 'number',      NF_INT,    m_dim)
     nstepsid    = def_var(ncid, 'nsteps',      NF_INT,    m_dim)
     mstepsid    = def_var(ncid, 'model_steps', NF_INT,    m_dim)
     iid         = def_var(ncid, 'i',           NF_INT,    m_dim)
     jid         = def_var(ncid, 'j',           NF_INT,    m_dim)
     kid         = def_var(ncid, 'k',           NF_INT,    m_dim)
     typeid      = def_var(ncid, 'type',        NF_INT,    m_dim)
     allocate(tracerid(num_prog_tracers))
     do n=1,num_prog_tracers
        if (n==index_temp) then
           tracerid(n) = def_var(ncid, 'heat', NF_DOUBLE, m_dim)
        else
           tracerid(n) = def_var(ncid, trim(T_prog(n)%name), NF_DOUBLE, m_dim)
        endif
     enddo

     call put_att(ncid, latid,      'long_name', 'latitude')
     call put_att(ncid, latid,      'units',     'degrees_N')
     call put_att(ncid, lonid,      'long_name', 'longitude')
     call put_att(ncid, lonid,      'units',     'degrees_E')
     call put_att(ncid, geodepthid, 'long_name', 'depth_from_z=0')
     call put_att(ncid, geodepthid, 'units',     'm')
     call put_att(ncid, stid,       'long_name', 'vertical position')
     if (vert_coordinate_class == DEPTH_BASED) then
        call put_att(ncid, stid,       'units',     'm')
     else
        call put_att(ncid, stid,       'units',     'N/m^2')
     endif
     call put_att(ncid, massid,     'long_name', 'mass')
     call put_att(ncid, massid,     'units',     'kg')
     call put_att(ncid, densityid,  'long_name', 'density')
     call put_att(ncid, densityid,  'units',     'kg/m^3')
     call put_att(ncid, uid,        'long_name', 'zonal_velocity')
     call put_att(ncid, uid,        'units',     'm/s')
     call put_att(ncid, vid,        'long_name', 'meridional_velocity')
     call put_att(ncid, vid,        'units',     'm/s')
     call put_att(ncid, wid,        'long_name', 'vertical_velocity')
     call put_att(ncid, wid,        'units',     'm/s')
     call put_att(ncid, entid,      'long_name', 'entrainment_velocity')
     call put_att(ncid, entid,      'units',     'm/s')
     call put_att(ncid, dragid,     'long_name', 'coefficient of drag (free blobs only)')
     call put_att(ncid, dragid,     'units',     '1/s')
     call put_att(ncid, ageid,      'long_name', 'age of blob')
     call put_att(ncid, ageid,      'units',     's')
     call put_att(ncid, heightid,   'long_name', 'blob height')
     call put_att(ncid, heightid,   'units',     'm')
     call put_att(ncid, h1id,       'long_name', 'stretching_function1')
     call put_att(ncid, h1id,       'units',     'none')
     call put_att(ncid, h2id,       'long_name', 'stretching_function2')
     call put_att(ncid, h2id,       'units',     'none')
     call put_att(ncid, stepid,     'long_name', 'blob_step_size')
     call put_att(ncid, stepid,     'units',     'seconds')
     call put_att(ncid, hashid,     'long_name', 'gridcell_identify')
     call put_att(ncid, hashid,     'units',     'none')
     call put_att(ncid, numberid,   'long_name', 'gridcell_counter')
     call put_att(ncid, numberid,   'units',     'none')
     call put_att(ncid, nstepsid,   'long_name', 'total_blob_steps')
     call put_att(ncid, nstepsid,   'units',     'none')
     call put_att(ncid, mstepsid,   'long_name', 'total_E_model_steps')
     call put_att(ncid, mstepsid,   'units',     'none')
     call put_att(ncid, iid,        'long_name', 'logical zonal cell index')
     call put_att(ncid, iid,        'units',     'none')
     call put_att(ncid, jid,        'long_name', 'logical meridional cell index')
     call put_att(ncid, jid,        'units',     'none')
     call put_att(ncid, kid,        'long_name', 'logical vertical cell index')
     call put_att(ncid, kid,        'units',     'none')
     call put_att(ncid, typeid,     'long_name', 'blob_type')
     call put_att(ncid, typeid,     'units',     'none')
     do n=1,num_prog_tracers
        if (n==index_temp) then
           call put_att(ncid, tracerid(n), 'long_name', 'heat_content')
           call put_att(ncid, tracerid(n), 'units',     'J')
        elseif (T_prog(n)%name(1:3)=='age') then
           call put_att(ncid, tracerid(n), 'long_name', trim(T_prog(n)%name)//'_content')
           call put_att(ncid, tracerid(n), 'units',     'kg/1e18')
        else
           call put_att(ncid, tracerid(n), 'long_name', trim(T_prog(n)%name)//'_content')
           call put_att(ncid, tracerid(n), 'units',     'kg')
        endif
     enddo

     mret = nf_put_att_int(ncid, NF_GLOBAL, 'file_format_major_version', NF_INT, 1, file_format_major_version)
     mret = nf_put_att_int(ncid, NF_GLOBAL, 'file_format_minor_version', NF_INT, 2, file_format_minor_version)

     ! End define mode
     mret = nf_enddef(ncid)

     m = 0
     do p=1,4
        ! We need to be able to distinguish between the four types of blobs
        ! when we restart.  So, we alot each a type.  type==1 is a static free
        ! blob and type==2 is a static bottom blob, type==3 is a dynamic free
        ! blob and type==4 is a dynamic bottom blob.  type==0 means that it is
        ! an empty entry.
        nullify(this)
        if(p==1 .and. associated(head_static_free))  this=>head_static_free
        if(p==2 .and. associated(head_static_bott))  this=>head_static_bott
        if(p==3 .and. associated(head_dynamic_free)) this=>head_dynamic_free
        if(p==4 .and. associated(head_dynamic_bott)) this=>head_dynamic_bott

        if (associated(this)) then
           fullcycle: do

              m = m+1
              call put_double(ncid, latid,      m, this%lat)
              call put_double(ncid, lonid,      m, this%lon)
              call put_double(ncid, geodepthid, m, this%geodepth)
              call put_double(ncid, stid,       m, this%st)
              call put_double(ncid, massid,     m, this%mass)
              call put_double(ncid, densityid,  m, this%density)
              call put_double(ncid, uid,        m, this%v(1))
              call put_double(ncid, vid,        m, this%v(2))
              call put_double(ncid, wid,        m, this%v(3))
              call put_double(ncid, entid,      m, this%ent)
              call put_double(ncid, dragid,     m, this%drag)
              call put_double(ncid, ageid,      m, this%age)
              call put_double(ncid, heightid,   m, this%height)
              call put_double(ncid, h1id,       m, this%h1)
              call put_double(ncid, h2id,       m, this%h2)
              call put_double(ncid, stepid,     m, this%step)
              call put_int(   ncid, hashid,     m, this%hash)
              call put_int(   ncid, numberid,   m, this%number)
              call put_int(   ncid, nstepsid,   m, this%nsteps)
              call put_int(   ncid, mstepsid,   m, this%model_steps)
              call put_int(   ncid, typeid,     m, p)
              call put_int(   ncid, iid,        m, this%i)
              call put_int(   ncid, jid,        m, this%j)
              call put_int(   ncid, kid,        m, this%k)
              do n=1,num_prog_tracers
                 call put_double(ncid, tracerid(n), m, this%tracer(n))
              enddo

              this=>this%next
              if(.not. associated(this)) exit fullcycle
           enddo fullcycle
        endif
     enddo

     if (m==0) then
        call put_double(ncid, latid,      1, 0.0)
        call put_double(ncid, lonid,      1, 0.0)
        call put_double(ncid, geodepthid, 1, 0.0)
        call put_double(ncid, massid,     1, 0.0)
        call put_double(ncid, densityid,  1, 0.0)
        call put_double(ncid, uid,        1, 0.0)
        call put_double(ncid, vid,        1, 0.0)
        call put_double(ncid, wid,        1, 0.0)
        call put_double(ncid, entid,      1, 0.0)
        call put_double(ncid, dragid,     1, 0.0)
        call put_double(ncid, ageid,      1, 0.0)
        call put_double(ncid, heightid,   1, 0.0)
        call put_double(ncid, h1id,       1, 0.0)
        call put_double(ncid, h2id,       1, 0.0)
        call put_double(ncid, stepid,     1, 0.0)
        call put_int(   ncid, hashid,     1, 0)
        call put_int(   ncid, numberid,   1, 0)
        call put_int(   ncid, nstepsid,   1, 0)
        call put_int(   ncid, mstepsid,   1, 0)
        call put_int(   ncid, typeid,     1, 0)
        call put_int(   ncid, iid,        1, 0)
        call put_int(   ncid, jid,        1, 0)
        call put_int(   ncid, kid,        1, 0)
        do n=1,num_prog_tracers
           call put_double(ncid, tracerid(n), 1, 0.0)
        enddo
     endif

     ! Finish up the blobs restart
     mret = nf_close(ncid)
     if (mret .ne. NF_NOERR) write(stderrunit,'(a)') 'blobs, writerestart: nf_close failed'
     deallocate(tracerid)

     print('(a22,i3,a3,i10)'),  'number of blobs on PE ',Info%pe_this,' =',m
     call mpp_sum(m)
     write(stdoutunit,'(a,i10)') 'global number of blobs     =', m

   end subroutine writerestart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that empties the blob lists at the end
! of a model run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine deallocateblobs

    type(ocean_blob_type), pointer :: this=>NULL()
    type(ocean_blob_type), pointer :: next=>NULL()
    type(ocean_blob_type), pointer :: head=>NULL()
    integer :: p

    do p=1,4
       if(p==1) head=>head_static_free
       if(p==2) head=>head_static_bott
       if(p==3) head=>head_dynamic_free
       if(p==4) head=>head_dynamic_bott

       this=>head
       if (associated(this)) then
          freecycle: do
             next=>this%next
             call free_blob_memory(this)
             this=>next
             if (.not. associated(this)) exit freecycle
          enddo freecycle
          nullify(head)
       endif
    enddo

   end subroutine deallocateblobs

end subroutine ocean_blob_end
! </SUBROUTINE>  NAME="ocean_blob_end"

!######################################################################
! <SUBROUTINE NAME="write_all_blobs">
!
! <DESCRIPTION>
! A convenient subroutine for debugging that dumps blob details from
! every list
! </DESCRIPTION>
!
subroutine write_all_blobs(time)
  character(len=*), intent(in) :: time

  call write_blobs(head_static_free,  'static free',    time)
  call write_blobs(head_static_bott,  'static bottom',  time)
  call write_blobs(head_dynamic_free, 'dynamic free',   time)
  call write_blobs(head_dynamic_bott, 'dynamic bottom', time)

end subroutine write_all_blobs
! </SUBROUTINE>  NAME="write_all_blobs"


!######################################################################
! <SUBROUTINE NAME="entrainment_checksum">
!
! <DESCRIPTION>
! Do the entrainment checksums
! </DESCRIPTION>
!
subroutine entrainment_checksum(Time, hfree, hbott, f2b, b2f)
  type(ocean_time_type), intent(in) :: Time
  type(ocean_blob_type), pointer    :: hfree
  type(ocean_blob_type), pointer    :: hbott
  type(ocean_blob_type), pointer    :: f2b
  type(ocean_blob_type), pointer    :: b2f

  type(ocean_blob_type), pointer :: this=>NULL()
  integer, dimension(isd:ied,jsd:jed,nk) :: di, dj, dk
  real, dimension(isd:ied,jsd:jed,nk) :: dmass, in, out
  real, dimension(isd:ied,jsd:jed,nk,num_prog_tracers) :: dtracer
  integer :: i,j,k,p,s,n
  integer :: stdoutunit

  stdoutunit = stdout()

  if (really_debug) then
     do p=3,6
        write(stdoutunit, *) ' '
        if(p==3) then
           this=>hfree
           write(stdoutunit, *) 'Entrainment, detrainment and converence details for the Dynamic Free Blobs'
        elseif(p==4) then
           this=>hbott
           write(stdoutunit, *) 'Entrainment, detrainment and converence details for the Dynamic Bottom Blobs'
        elseif(p==5) then
           this=>f2b
           write(stdoutunit, *) 'Entrainment, detrainment and converence details for the Dynamic Free Blobs transferring to Bottom'
        elseif(p==6) then
           this=>b2f
           write(stdoutunit, *) 'Entrainment, detrainment and converence details for the Dynamic Bottom Blobs transferring to Free'
        endif
        write(stdoutunit, '(2(a10,x),3(a1,a3),a1,x,a4,x,5(x,a21))') 'hash','number',&
             '(','i',',','j',',','k',')','s', &
             'dmass','dheat','dsalt',         &
             'convergence','divergence'
        if (associated(this)) then
           blobcycle0: do
              i=this%di(s)
              j=this%dj(s)
              k=this%dk(s)
              print ('(2(i10,x),3(a1,i3),a1,x,4x,x,4(x,21x),(x,es21.14))'), this%hash, this%number, &
                   '(',i,',',j,',',k,')', 0, &
                   this%mass_out(s)
              if (this%nfrac_steps>0) then
                 do s=1,this%nfrac_steps
                    i=this%di(s)
                    j=this%dj(s)
                    k=this%dk(s)
                    print ('(2(i10,x),3(a1,i3),a1,x,i4,x,5(x,es21.14))'), this%hash, this%number, &
                         '(',i,',',j,',',k,')', s, &
                         this%entrainment(s,0)+this%detrainment(s,0), &
                         this%entrainment(s,index_temp)+this%detrainment(s,index_temp), &
                         this%entrainment(s,index_salt)+this%detrainment(s,index_salt),&
                         this%mass_in(s), this%mass_out(s)
                 enddo
              endif
              this=>this%next
              if (.not. associated(this)) exit blobcycle0
           enddo blobcycle0
        endif
     enddo
     write(stdoutunit, *) 'end of list'
  endif

  write(stdoutunit, *) ' '
  write(stdoutunit, *) 'Entrainment and convergence checksums for ocean_blob_mod'
  call write_timestamp(Time%model_time)
  do p=3,6
     if(p==3) then
        this=>hfree
        write(stdoutunit, *) '==>Checksum for dynamic free blobs'
     elseif(p==4) then
        this=>hbott
        write(stdoutunit, *) '==>Checksum for dynamic bottom blobs'
     elseif(p==5) then
        this=>f2b
        write(stdoutunit, *) '==>Checksum for dynamic free blobs interacting with bottom'
     elseif(p==6) then
        this=>b2f
        write(stdoutunit, *) '==>Checksum for dynamic bottom blobs separating from bottom'
     endif

     di(:,:,:)    = 0
     dj(:,:,:)    = 0
     dk(:,:,:)    = 0
     dmass(:,:,:) = 0.0
     in(:,:,:)    = 0.0
     out(:,:,:)   = 0.0
     dtracer(:,:,:,:) = 0.0

     if (associated(this)) then
        blobcycle: do
           i=this%di(0)
           j=this%dj(0)
           k=this%dk(0)

           if (isc<= i .and. i<=iec .and.&
               jsc<= j .and. j<=jec) then

              di(i,j,k) = di(i,j,k) + i
              dj(i,j,k) = dj(i,j,k) + j
              dk(i,j,k) = dk(i,j,k) + k

              out(i,j,k) = out(i,j,k) + this%mass_out(0)
           endif

           if (this%nfrac_steps>0) then
              do s=1,this%nfrac_steps
                 i=this%di(s)
                 j=this%dj(s)
                 k=this%dk(s)

                 if (isc<= i .and. i<=iec .and.&
                     jsc<= j .and. j<=jec) then

                    di(i,j,k) = di(i,j,k) + i
                    dj(i,j,k) = dj(i,j,k) + j
                    dk(i,j,k) = dk(i,j,k) + k

                    dmass(i,j,k) = dmass(i,j,k) + this%entrainment(s,0) + this%detrainment(s,0)

                    in(i,j,k)  = in(i,j,k) + this%mass_in(s)
                    out(i,j,k) = out(i,j,k) + this%mass_out(s)
                    do n=1,num_prog_tracers
                       dtracer(i,j,k,n) = dtracer(i,j,k,n) + this%entrainment(s,n) + this%detrainment(s,n)
                    enddo
                 endif
              enddo
              i=this%di(this%nfrac_steps)
              j=this%dj(this%nfrac_steps)
              k=this%dk(this%nfrac_steps)
              if (isc<= i .and. i<=iec .and.&
                  jsc<= j .and. j<=jec) then
                 dmass(i,j,k) = dmass(i,j,k) + this%dmass
                 do n=1,num_prog_tracers
                    dtracer(i,j,k,n) = dtracer(i,j,k,n) + this%dtracer(n)
                 enddo
              endif
              if (isc<= i .and. i<=iec .and. jsc<= j .and. j<=jec) in(i,j,k)  = in(i,j,k) + this%mass_in(this%nfrac_steps)
           endif

           this=>this%next
           if (.not. associated(this)) exit blobcycle
        enddo blobcycle
     endif
     call write_chksum_3d_int('zonal index (i)', di(COMP,:))
     call write_chksum_3d_int('meridional index (j)', dj(COMP,:))
     call write_chksum_3d_int('vertical index (k)', dk(COMP,:))
     call write_chksum_3d('mass entrainment/detrainment', dmass(COMP,:))
     call write_chksum_3d('mass convergence', in(COMP,:))
     call write_chksum_3d('mass divergence', out(COMP,:))
     call write_chksum_3d('salt entrainment/detrainment', dtracer(COMP,:,index_salt))
     call write_chksum_3d('temp entrainment/detrainment', dtracer(COMP,:,index_temp))
  enddo
  write(stdoutunit, *) 'End checksums'

end subroutine entrainment_checksum
! </SUBROUTINE>  NAME="entrainment_checksum"

end module ocean_blob_mod
! Ninja comment
