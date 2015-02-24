module ocean_blob_dynamic_bottom_mod
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This module runs the dynamic bottom blob implementation of the embedded 
! Lagrangian buoyancy blob framework. The module forms new dynamic 
! bottom blobs, integrates the properties of existing blobs, and handles
! the transfer of free blobs to bottom blobs.
!</OVERVIEW>
!
!<DESCRIPTION>
! Bottom blobs are formed using the subroutine dynamic_bottom_form_new,
! which is called from the main blob driver module.  Bottom blobs are
! formed explicitly in time, directly after the integration of existing
! blobs.
!
! The properties of a bottom blob are also integrated in this module,
! that is, position, velocity, mass and tracer content.  Position and
! velocity are integrated using an adaptive step Runge-Kutta scheme.
! There are several schemes available of varying order.
!
! The module also recieves blobs that are transferring from the free
! blob dynamic regime to the bottom blob dynamic regime (i.e. free
! blobs that have interacted with topography).
!</DESCRIPTION>
!
!<INFO>
!
! <REFERENCE>
!  Price, J.F., Baringer, M.O'N. (1994) Outflows and deep water production 
!  by marginal seas. Progress in Oceanography 33(3), 161-200.
! </REFERENCE>
!
! <REFERENCE>
!  Campin, J-.M., Goosse, H. (1999) Parameterization of a density-driven
!  downsloping flow for a coarse-resolution ocean model in z-coordinate.
!  Tellus 51A(3), 412-430.
! </REFERENCE>
!
! <REFERENCE>
!  Bogacki, P., Shampine, L.F. (1989) A 3(2) pair of Runge-Kutta formulas.
!  Applied Mathematical Letters 2(4), 321-325.
! </REFERENCE>
!
! <REFERENCE>
!  Cash, J.R., Karp, A.H. (1990) A variable order Runge-Kutta method for 
!  initial value problems with rapidly varying right-hand sides.
!  ACM Transactions on Mathematical Software 16(3), 201-222.
! </REFERENCE>
!
!</INFO>
!<NAMELIST NAME="ocean_blob_dynamic_bottom_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module.
!  Default is use_this_module=.false.
!  </DATA>
!
!  <DATA NAME="blob_overflow_mu" TYPE="real">
!  Frictional dissipation rate used for calculating initial
!  properties of bottom blobs.  Corresponds to mu in Campin and 
!  Goosse (1999).  Units are 1/s.
!  Default is blob_overflow_mu=1.0e-4
!  </DATA>
!
!  <DATA NAME="blob_overflow_delta" TYPE="real">
!  Fraction of a grid cell participating in an overflow event.
!  Corresponds to delta in Campin and Goosse (1999).  Dimensionless.
!  Default is blob_overflow_mu=1.0e-4
!  </DATA>
!
!  <DATA NAME="drag" TYPE="real">
!  Coefficient of drag used for bottom stress drag.
!  Corresponds to Cd in Price and Baringer (1994). Dimensionless.
!  Default is drag=3.0e-3
!  </DATA>
!
!  <DATA NAME="det_param" TYPE="real">
!  The detrainment parameter (kg m^2/s). 
!  Corresponds to Gamma in the notes.
!  Default is det_param=5.0e-8
!
!  </DATA>
!  <DATA NAME="max_detrainment" TYPE="real">
!  The Maximum allowable detrainment velocity (m/s).
!  Default is max_detrainment=1.0e-3
!  </DATA>
!
!  <DATA NAME="rel_error" TYPE="real">
!  Relative error for the RK scheme (dimensionless).
!  A smaller number is more accurate, but, 
!  is more computationally expensive. Corresponds to
!  zeta* in the notes.
!  Must be 0<rel_error<=1.0
!  Default is rel_error=0.01
!  </DATA>
!
!  <DATA NAME="safety_factor" TYPE="real">
!  Safety factor for the RK scheme (dimensionless).
!  A smaller number should reduce the number
!  of rejected steps, but, decreases the locally
!  extrapolated step.  Corresponds to varrho in 
!  the notes.
!  Must be 0<safety_factor<=1.0
!  Default is safety_factor=0.8
!  </DATA>
!
!  <DATA NAME="minstep" TYPE="real">
!  Minimum step size (in seconds) for a blob.
!  Default is minstep=9.0
!  </DATA>
!
!  <DATA NAME="elastic" TYPE="real">
!  The elasticity of a blob's collision with
!  the topography.  Corresponds to epsilon in 
!  the notes. Should have values 0<=elastic<=1.0
!  Values greater than 1 would be super-elastic,
!  and values less than 0 would send the blob
!  in the opposite direction than it should be 
!  going in.
!  Default is elastic=1.0
!  </DATA>
!
!  <DATA NAME="min_do_levels" TYPE="integer">
!  Minimum number of deep ocean levels for
!  overflows to be considered.  That is, how many
!  k levels lower should the deep ocean water column
!  be than the shelf/shallow ocean column.  Value
!  must be greater than 0.
!  Default is min_do_levels=1
!  </DATA>
!
!  <DATA NAME="rho_threshold" TYPE="real">
!  The density difference required before a blob
!  is formed.  rho_threshold must be greater than
!  zero.
!  Default is rho_threshold=0.01
!  </DATA>
!
!  <DATA NAME="large_speed" TYPE="real">
!  A value for error checking.  If the speed of a
!  blob exceeds large_speed in any of x,y,z then
!  a warning flag is raised.
!  Default is large_speed=10.0
!  </DATA>
!
!  <DATA NAME="no_rotation" TYPE="real">
!  Sets the coriolis parameter to zero regardless
!  of latitude
!  Default is no_rotation=.false.
!  </DATA>
!
!  <DATA NAME="critical_richardson" TYPE="real">
!  The critical Richardson number for the entrainment
!  velocity.  Default is based on Price and Baringer
!  (1994).
!  Default is critical_richardson=0.8
!  </DATA>
!
!</NAMELIST>
!

use constants_mod,    only: rad_to_deg, deg_to_rad, epsln, pi
use diag_manager_mod, only: register_diag_field, send_data

use fms_mod,         only: stdout, stdlog, open_namelist_file, FATAL, WARNING
use fms_mod,         only: mpp_error, check_nml_error, close_file
use mpp_mod,         only: NULL_PE, mpp_sum, mpp_max, mpp_recv, mpp_send, mpp_sync_self
use mpp_mod,         only: CLOCK_LOOP, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,         only: mpp_get_current_pelist, mpp_root_pe, mpp_npes
use mpp_domains_mod, only: mpp_update_domains

use ocean_blob_util_mod,  only: insert_blob, unlink_blob, blob_delete, interp_tcoeff, interp_ucoeff
use ocean_blob_util_mod,  only: check_ijcell, free_blob_memory, kill_blob, count_blob, write_blobs
use ocean_blob_util_mod,  only: allocate_interaction_memory, reallocate_interaction_memory
use ocean_blob_util_mod,  only: E_and_L_totals, check_cyclic
use ocean_density_mod,    only: density
use ocean_parameters_mod, only: onehalf, onethird, twothirds, rho0, rho0r, grav, omega_earth
use ocean_parameters_mod, only: PRESSURE_BASED, DEPTH_BASED
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, blob_diag_type
use ocean_types_mod,      only: ocean_domain_type, blob_grid_type, ocean_external_mode_type
use ocean_types_mod,      only: ocean_lagrangian_type, ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,      only: ocean_blob_type, ocean_prog_tracer_type, ocean_density_type
use ocean_util_mod,       only: write_timestamp

implicit none

private

! set module types 
type(ocean_grid_type),   pointer :: Grd   => NULL()
type(ocean_domain_type), pointer :: Dom   => NULL()
type(ocean_domain_type), pointer :: Bdom  => NULL()
type(blob_grid_type),    pointer :: Info  => NULL()

real :: dtime
real :: dtime_yr
real :: det_factor
real :: ent_factor

real, parameter :: sixonpi   = 6.0/pi
real :: two_omega

logical :: module_is_initialized=.false.

!RK coefficients
real, dimension(:,:), allocatable :: beta
real, dimension(:),   allocatable :: lgamma
real, dimension(:),   allocatable :: hgamma
real :: order
real :: orderp1_r
integer :: ns, nsp1, total_ns, total_nsp1, total_nsp2

! For topography
real, dimension(:,:), allocatable :: ht, dht_dx, dht_dy

!initial stretching function for blobs
real, dimension(:,:,:), allocatable :: blobh1
real, dimension(:,:,:), allocatable :: blobh2

integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1
integer :: isd,ied,jsd,jed
integer :: isc,iec,jsc,jec
integer :: isg,ieg,jsg,jeg
integer :: nk
integer :: isbd, iebd, jsbd, jebd

integer :: method

!clock variables
integer :: id_clock_rk

!Diagnostics
logical :: used
integer :: id_free_to_bot
integer :: id_new_blobs

!variables for initialising blobs
real :: overflow_factor
real :: dtime_rad_to_deg
real, dimension(:,:,:), allocatable :: topog_step
real, dimension(:,:,:), allocatable :: topog_slope
real, dimension(:,:,:), allocatable :: vert_factor
integer, dimension(4) :: ip, jp, ish, jsh, ieh, jeh

real, dimension(:,:,:), allocatable :: hb
real, dimension(:,:,:), allocatable :: xb
real, dimension(:,:,:), allocatable :: yb

real, dimension(:,:,:), allocatable :: umask !umask, but, halo=2

integer, dimension(:,:), allocatable :: kmt
integer, dimension(:,:), allocatable :: kmu
integer, dimension(5) :: it, jt
integer, dimension(4) :: iu, ju 

type, private :: blob_buffer_type
   integer :: numblobs
   integer :: size
   integer :: pe
   integer, allocatable, dimension(:,:)   :: integer_buff
   real,    allocatable, dimension(:,:)   :: real_buff
   integer, allocatable, dimension(:,:,:) :: history_integer_buff
   real,    allocatable, dimension(:,:,:) :: history_real_buff
end type blob_buffer_type

integer :: rea_buff_size, int_buff_size
integer :: hist_rea_buff_size, hist_int_buff_size
type(blob_buffer_type), pointer :: Ebuffer_out,  Ebuffer_in
type(blob_buffer_type), pointer :: Wbuffer_out,  Wbuffer_in
type(blob_buffer_type), pointer :: Nbuffer_out,  Nbuffer_in
type(blob_buffer_type), pointer :: Sbuffer_out,  Sbuffer_in
type(blob_buffer_type), pointer :: NEbuffer_out, NEbuffer_in
type(blob_buffer_type), pointer :: NWbuffer_out, NWbuffer_in
type(blob_buffer_type), pointer :: SEbuffer_out, SEbuffer_in
type(blob_buffer_type), pointer :: SWbuffer_out, SWbuffer_in
type(blob_buffer_type), pointer :: diagbuffer_out, diagbuffer_in

! buffers and buffer variables to send blobs to adjacent PE's
integer, parameter :: delta_buffer = 25   ! Size by which to increment the buffer

public blob_dynamic_bottom_init
public blob_dynamic_bottom_update
public transfer_free_to_bottom
public dynamic_bottom_form_new
public blob_dynamic_bottom_end

! variables controlled by other namelists
logical :: debug_this_module
logical :: really_debug
integer :: vert_coordinate_class ! the class of vertical coordinate
integer :: vert_coordinate
logical :: blob_diag
real    :: small_mass
logical :: bitwise_reproduction

! namelist defaults
logical :: use_this_module   = .false.
character(len=10) :: update_method = 'BS_RK3(2)'
real    :: blob_overflow_mu    = 1.0e-4  ! frictional dissipation rate (sec^-1) at bottom  
real    :: blob_overflow_delta = 0.3333  ! fraction of a grid cell participating in overflow 
real    :: drag                = 3.0e-3
logical :: enforce_big_blobs   = .false.
real    :: det_param           = 5.0e-8
real    :: max_detrainment     = 1.0e-3
real    :: rel_error           = 0.01
real    :: safety_factor       = 0.8
real    :: minstep             =  9.0
real    :: first_step          = 50.0
real    :: elastic             = 1.0
integer :: min_do_levels       = 1
real    :: rho_threshold       = 0.01
logical :: accept_free_blobs   = .true.
real    :: large_speed         = 10.0 !m/s
logical :: no_rotation         = .false.
real    :: critical_richardson = 0.8
real    :: blob_height         = 100.0 !m
real    :: blobs_south_of      = -45.0
real    :: blobs_north_of      =  45.0

namelist /ocean_blob_dynamic_bottom_nml/ use_this_module, update_method,   &
     blob_overflow_mu, blob_overflow_delta, drag, enforce_big_blobs,       &
     det_param, max_detrainment, rel_error, safety_factor, elastic,        &
     minstep, first_step, min_do_levels, rho_threshold, accept_free_blobs, &
     large_speed, no_rotation, critical_richardson, blob_height,           &
     blobs_south_of, blobs_north_of
     
contains
!#######################################################################
! <SUBROUTINE NAME="blob_dynamic_bottom_init">
!
! <DESCRIPTION>
! Initialises the dynamic free blobs by checking the namelist and also
! inherited namelists (from ocean_blob_nml).  Also sets up some useful
! constants (including spatially varying constants) - particularly for
! the formation of bottom blobs.  It also allocates memory to special 
! halo=2 masks and sets up the blob buffers for sending blobs from one 
! PE to another.
! </DESCRIPTION>
!
subroutine blob_dynamic_bottom_init(Time, Grid, Domain, Blob_domain, PE_Info,          &
                                    debug, big_debug, bitwise, num_tracers,            &
                                    itemp, isalt, dtimein, ver_coord_class, ver_coord, &
                                    blob_diagnostics, bott_minstep, bott_total_ns,     &
                                     smallmass, use_dyn_bot)

  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_domain_type),      intent(in), target :: Blob_domain
  type(blob_grid_type),         intent(in), target :: PE_Info
  logical, intent(in)  :: debug
  logical, intent(in)  :: big_debug
  logical, intent(in)  :: bitwise
  integer, intent(in)  :: num_tracers
  integer, intent(in)  :: itemp
  integer, intent(in)  :: isalt
  real,    intent(in)  :: dtimein
  integer, intent(in)  :: ver_coord_class
  integer, intent(in)  :: ver_coord
  logical, intent(in)  :: blob_diagnostics
  real,    intent(out) :: bott_minstep
  integer, intent(out) :: bott_total_ns
  real,    intent(in)  :: smallmass
  logical, intent(out) :: use_dyn_bot

  real, dimension(:,:), allocatable :: slope_x
  real, dimension(:,:), allocatable :: slope_y
  real, dimension(:,:), allocatable :: coeff1
  real, dimension(:,:), allocatable :: coeff2
  real, parameter                   :: secs_in_year_r = 1.0 / (86400.0 * 365.25)
  integer :: ioun, ierr, io_status
  integer :: stdoutunit,stdlogunit
  integer :: i,j,m,iip,jjp

  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
         '==>Error in ocean_blob_dynamic_free_mod (ocean_blob_dynamic_free_init):' &
         //' module already initialized')
  endif 

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_blob_dynamic_bottom_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_blob_dynamic_bottom_nml)  
  write (stdlogunit,ocean_blob_dynamic_bottom_nml)
  ierr = check_nml_error(io_status,'ocean_blob_dynamic_bottom_nml')
  call close_file (ioun)

  module_is_initialized = .true.
  
  use_dyn_bot = use_this_module

  Grd   => Grid
  Dom   => Domain
  Bdom  => Blob_domain
  Info  => PE_Info

  if (no_rotation) then
     two_omega = 0.0
  else
     two_omega = 2*omega_earth
  endif

  debug_this_module = debug

  index_temp           = itemp
  index_salt           = isalt
  num_prog_tracers     = num_tracers
  bitwise_reproduction = bitwise
  dtime                = dtimein
  small_mass           = smallmass
  bott_minstep         = minstep
  blob_diag            = blob_diagnostics
  really_debug         = big_debug

  id_free_to_bot = register_diag_field('ocean_model','free_to_bot', Time%model_time, &
       'free blobs interacting with topography', 'number of blobs')

  if (.not. use_this_module) return

  id_new_blobs   = register_diag_field('ocean_model','new_bottom_blobs', Time%model_time, &
       'new bottom blobs', 'number of blobs')

  dtime_rad_to_deg = dtime*rad_to_deg
  dtime_yr         = dtime*secs_in_year_r

  vert_coordinate_class = ver_coord_class
  vert_coordinate       = ver_coord

  isd=Dom%isd; ied=Dom%ied; jsd=Dom%jsd; jed=Dom%jed
  isc=Dom%isc; iec=Dom%iec; jsc=Dom%jsc; jec=Dom%jec
  isg=Dom%isg; ieg=Dom%ieg; jsg=Dom%jsg; jeg=Dom%jeg
  nk = Grd%nk

  isbd=Bdom%isd; iebd=Bdom%ied; jsbd=Bdom%jsd; jebd=Bdom%jed

  if (min_do_levels<1) then
     write(stdoutunit,'(a)')&
          '==>Error in ocean_blob_dynamic_bottom_mod (ocean_blob_dynamic_bottom_init):' &
          //' invalid value for min_do_levels chosen.  Should be greater than or equal to one.'
     write(stdoutunit,'(a,i3)') 'min_do_levels =',min_do_levels
     call mpp_error(FATAL,&
          '==>Error in ocean_blob_dynamic_bottom_mod (ocean_blob_dynamic_bottom_init):' &
          //' invalid value for min_do_levels chosen.  Should be greater than or equal to one.')
  endif

  allocate(kmt(isbd:iebd, jsbd:jebd))
  allocate(kmu(isbd:iebd, jsbd:jebd))
  kmt(isc:iec,jsc:jec) = Grd%kmt(isc:iec, jsc:jec)
  kmu(isc:iec,jsc:jec) = Grd%kmu(isc:iec, jsc:jec)
  call mpp_update_domains(kmt(:,:), Bdom%domain2d, complete=.false.)
  call mpp_update_domains(kmu(:,:), Bdom%domain2d, complete=.true.)

  allocate(ht(isbd:iebd, jsbd:jebd))
  allocate(dht_dx(isbd:iebd, jsbd:jebd))
  allocate(dht_dy(isbd:iebd, jsbd:jebd))
  ht(isc:iec,jsc:jec) = Grd%ht(isc:iec,jsc:jec)
  dht_dx(isd:iec,jsd:jec) = Grd%dht_dx(isd:iec,jsd:jec)
  dht_dy(isd:iec,jsd:jec) = Grd%dht_dy(isd:iec,jsd:jec)
  call mpp_update_domains(ht(:,:),     Bdom%domain2d, complete=.false.)
  call mpp_update_domains(dht_dx(:,:), Bdom%domain2d, complete=.false.)
  call mpp_update_domains(dht_dy(:,:), Bdom%domain2d, complete=.true.)

  ! special treatment for boundaries if they are NULL_PEs.
  ! Make sure that kmt==0, kmu==0, ht==0.0, dht_dx==0 
  ! and dht_dy==0 in solid boundaries.
  if (Info%pe_S==NULL_PE) then
     kmt(:,jsbd:jsd) = 0
     kmu(:,jsbd:jsd) = 0
     ht (:,jsbd:jsd) = 0.0
     dht_dx(:,jsbd:jsd)  = 0.0
     dht_dy(:,jsbd:jsd)  = 0.0
  endif
  if (Info%pe_N==NULL_PE) then
     kmt(:,jed:jebd) = 0
     kmu(:,jed:jebd) = 0
     ht (:,jed:jebd) = 0.0
     dht_dx(:,jed:jebd) = 0.0
     dht_dy(:,jed:jebd) = 0.0
  endif
  if (Info%pe_E==NULL_PE) then
     kmt(ied:iebd,:) = 0
     kmu(ied:iebd,:) = 0
     ht( ied:iebd,:) = 0.0
     dht_dx(ied:iebd,:) = 0.0
     dht_dy(ied:iebd,:) = 0.0
  endif
  if (Info%pe_W==NULL_PE) then
     kmt(isbd:isd,:) = 0
     kmu(isbd:isd,:) = 0
     ht (isbd:isd,:) = 0.0
     dht_dx(isbd:isd,:) = 0.0
     dht_dy(isbd:isd,:) = 0.0
  endif

  allocate(umask(isbd:iebd,jsbd:jebd,1:nk))
  umask(isc:iec,jsc:jec,1:nk) = Grd%umask(isc:iec,jsc:jec, 1:nk)
  call mpp_update_domains(umask(:,:,1:nk), Bdom%domain2d)
  ! special treatment for boundaries if they are NULL_PEs.
  ! Make sure that umask==0
  if (Info%pe_S==NULL_PE) umask(:,jsbd:jsd,:) = 0.0
  if (Info%pe_N==NULL_PE) umask(:,jed:jebd,:) = 0.0
  if (Info%pe_E==NULL_PE) umask(ied:iebd,:,:) = 0.0
  if (Info%pe_W==NULL_PE) umask(isbd:isd,:,:) = 0.0

  allocate(coeff1(isc:iec,jsc:jec))
  allocate(coeff2(isc:iec,jsc:jec))

  ! Calculate the depth, lon and lat of initial blob positions
  allocate(hb(isc:iec,jsc:jec,4))
  allocate(xb(isc:iec,jsc:jec,4))
  allocate(yb(isc:iec,jsc:jec,4))

  ! Calculate the stretching function in the places we form blobs
  ! We do not need to worry about special treatment for solid wall 
  ! because we will never form a blob at a solid wall boundary
  allocate(blobh1(isc:iec,jsc:jec,4))
  allocate(blobh2(isc:iec,jsc:jec,4))

  m=1
  coeff1(:,:) = Grd%dtw(isc+1:iec+1,jsc:jec)/Grd%dxte(isc:iec,jsc:jec)
  coeff2(:,:) = Grd%dte(isc  :iec  ,jsc:jec)/Grd%dxte(isc:iec,jsc:jec)

  hb(:,:,m) = Grd%ht(isc:iec,jsc:jec)*coeff1(:,:) + Grd%ht(isc+1:iec+1,jsc:jec)*coeff2(:,:)
  xb(:,:,m) = Grd%xt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%xt(isc+1:iec+1,jsc:jec)*coeff2(:,:) 
  yb(:,:,m) = Grd%yt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%yt(isc+1:iec+1,jsc:jec)*coeff2(:,:) 
  blobh1(:,:,m) = Grd%h1t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h1t(isc+1:iec+1,jsc:jec)*coeff2(:,:)
  blobh2(:,:,m) = Grd%h2t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h2t(isc+1:iec+1,jsc:jec)*coeff2(:,:)

  m=2
  coeff1(:,:) = Grd%dts(isc:iec,jsc+1:jec+1)/Grd%dytn(isc:iec,jsc:jec)
  coeff2(:,:) = Grd%dtn(isc:iec,jsc  :jec  )/Grd%dytn(isc:iec,jsc:jec)

  hb(:,:,m) = Grd%ht(isc:iec,jsc:jec)*coeff1(:,:) + Grd%ht(isc:iec,jsc+1:jec+1)*coeff2(:,:) 
  xb(:,:,m) = Grd%xt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%xt(isc:iec,jsc+1:jec+1)*coeff2(:,:) 
  yb(:,:,m) = Grd%yt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%yt(isc:iec,jsc+1:jec+1)*coeff2(:,:) 
  blobh1(:,:,m) = Grd%h1t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h1t(isc:iec,jsc+1:jec+1)*coeff2(:,:)
  blobh2(:,:,m) = Grd%h2t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h2t(isc:iec,jsc+1:jec+1)*coeff2(:,:)

  m=3
  coeff1(:,:) = Grd%dte(isc-1:iec-1,jsc:jec)/Grd%dxte(isc-1:iec-1,jsc:jec)
  coeff2(:,:) = Grd%dtw(isc  :iec,  jsc:jec)/Grd%dxte(isc-1:iec-1,jsc:jec)

  hb(:,:,m) = Grd%ht(isc:iec,jsc:jec)*coeff1(:,:) + Grd%ht(isc-1:iec-1,jsc:jec)*coeff2(:,:) 
  xb(:,:,m) = Grd%xt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%xt(isc-1:iec-1,jsc:jec)*coeff2(:,:) 
  yb(:,:,m) = Grd%yt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%yt(isc-1:iec-1,jsc:jec)*coeff2(:,:) 
  blobh1(:,:,m) = Grd%h1t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h1t(isc-1:iec-1,jsc:jec)*coeff2(:,:)
  blobh2(:,:,m) = Grd%h2t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h2t(isc-1:iec-1,jsc:jec)*coeff2(:,:)

  m=4
  coeff1(:,:) = Grd%dtn(isc:iec,jsc-1:jec-1)/Grd%dytn(isc:iec,jsc-1:jec-1)
  coeff2(:,:) = Grd%dts(isc:iec,jsc  :jec  )/Grd%dytn(isc:iec,jsc-1:jec-1)

  hb(:,:,m) = Grd%ht(isc:iec,jsc:jec)*coeff1(:,:) + Grd%ht(isc:iec,jsc-1:jec-1)*coeff2(:,:)
  xb(:,:,m) = Grd%xt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%xt(isc:iec,jsc-1:jec-1)*coeff2(:,:)
  yb(:,:,m) = Grd%yt(isc:iec,jsc:jec)*coeff1(:,:) + Grd%yt(isc:iec,jsc-1:jec-1)*coeff2(:,:)
  blobh1(:,:,m) = Grd%h1t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h1t(isc:iec,jsc-1:jec-1)*coeff2(:,:)
  blobh2(:,:,m) = Grd%h2t(isc:iec,jsc:jec)*coeff1(:,:) + Grd%h2t(isc:iec,jsc-1:jec-1)*coeff2(:,:)

  solver_method: select case (trim(update_method))
  case('BS_RK3(2)')
     method = 1

     allocate(beta(1:3,0:2)); beta(:,:) = 0.
     allocate(lgamma(0:3)); lgamma(:)   = 0. !2nd order coefficients
     allocate(hgamma(0:3)); hgamma(:)   = 0. !3rd order coefficients
     
     beta(1,0)   =    1/2.
     beta(2,0:1) = (/ 0.,    3/4.       /)
     beta(3,0:2) = (/ 2/9.,  1/3., 4/9. /)
     
     hgamma(0:3) = (/ 2/9.,  1/3., 4/9., 0.   /)
     lgamma(0:3) = (/ 7/24., 1/4., 1/3., 1/8. /)
     
     order     = 2.
     orderp1_r = 1/(order + 1.)

  case('DP_RK5(4)')
     method = 2

     write(stdoutunit,'(a,/,a)')&
          '==>Warning in ocean_blob_dynamic_bottom_mod (ocean_blob_dynamic_bottom_init):' &
          //' The Dormand-Prince 5(4) method chosen.  This scheme is buggy.',             &
          '   It is suggested to use the Cash-Karp scheme by setting the namelist '      &
          //'variable update_method=CK_RK5(4)'
     
     call mpp_error(WARNING,&
          '==>Warning in ocean_blob_dynamic_bottom_mod (ocean_blob_dynamic_bottom_init):' &
          //' The Dormand-Prince 5(4) method chosen.  This scheme is buggy.')

     ! The scheme is buggy because it produces large error estimates compared to the CK or
     ! BS schemes.  Not sure why it is so.
     ! The large error estimates mean that it takes a lot of sub-cycles to for each blob.

     allocate(beta(1:6,0:5)); beta(:,:) = 0.
     allocate(lgamma(0:6));   lgamma(:) = 0. !4th order coefficients
     allocate(hgamma(0:6));   hgamma(:) = 0. !5th order coefficients
     
     beta(1,0)   =        1/5.
     beta(2,0:1) = (/     3/40.  ,       9/40.  /)
     beta(3,0:2) = (/    44/45.  ,     -56/15.  ,    32/9.   /)
     beta(4,0:3) = (/ 19372/6561.,  -25360/2187., 64448/6561., -212/729. /)
     beta(5,0:4) = (/  9017/3168.,    -355/33.  , 46732/5247.,   49/176., -5103/18656./)
     beta(6,0:5) = (/    35/384. ,       0.     ,   500/1113.,  125/192., -2187/6784. , 11/84./)
     
     hgamma(0:6) = (/   35/384.  , 0.,  500/1113. , 125/192.,  -2187/6784.  ,  11/84. , 0.    /)
     lgamma(0:6) = (/ 5179/57600., 0., 7571/16695., 393/640., -92097/339200., 187/2100., 1/40. /)
     
     order     = 4.
     orderp1_r = 1/(order + 1.)

  case('CK_RK5(4)')
     method = 3

     allocate(beta(1:5,0:4)); beta(:,:) = 0.
     allocate(lgamma(0:5));   lgamma(:) = 0. !4th order coefficients
     allocate(hgamma(0:5));   hgamma(:) = 0. !5th order coefficients
     
     beta(1,0)   =        1/5.
     beta(2,0:1) = (/     3/40.   ,   9/40.  /)
     beta(3,0:2) = (/     3/10.   ,  -9/10. ,   6/5.     /)
     beta(4,0:3) = (/   -11/54.   ,   5/2.  , -70/27.   ,    35/27.    /)
     beta(5,0:4) = (/  1631/55296., 175/512., 575/13824., 44275/110592.,   253/4096. /)
     
     hgamma(0:5) = (/   37/378.  , 0.,   250/621.  ,   125/594.  ,  0.       , 512/1771. /)
     lgamma(0:5) = (/ 2825/27648., 0., 18575/48384., 13525/55296., 277/14336.,   1/4.    /)
     
     order     = 4.
     orderp1_r = 1/(order + 1.)

  case default 
     write(stdoutunit,'(a)')&
          '==>Error in ocean_blob_dynamic_bottom_mod (ocean_blob_dynamic_bottom_init):' &
          //' invalid solver chosen ('//update_method//').  Check update_method in namelist.'
     call mpp_error(FATAL,&
          '==>Error in ocean_blob_dynamic_bottom_mod (ocean_blob_dynamic_bottom_init):' &
          //' invalid solver chosen ('//update_method//').  Check update_method in namelist.')
  endselect solver_method

  ns            = ubound(lgamma,1)   ! number of partial steps in a fractional step
  nsp1          = ns+1               ! number of partial steps in a fractional step+1
  total_ns      = ceiling(dtime/minstep) !maximum number of partial steps in a full step
  total_nsp1    = total_ns+1
  total_nsp2    = total_nsp1+1
  bott_total_ns = total_ns

  ! This collects all the constant terms together in the expression
  ! for the rate of change of mass.
  ! For detrainment:
  ! dm/dt = rhoL A D (A=blob surface area, D=detrainment rate)
  !       = -rhoL A Gamma/|rhoL - rhoE|
  !       = -m**2/3 Gamma rhoL (36pi)**1/3 / (rhoL**2/3 |rhoL-rhoE|)
  ! For entrainment:
  ! dm/dt = rhoE A E (A=blob surface area, E=entrainment rate)
  !       = -m**2/3 rhoE (36pi)**1/3 / (rhoL**2/3)
  !
  ! In the Boussinesq case rhoE=rhoL=rho0 (excepting |rhoL-rhoE|)
  ! we also need to take into account that this is done for each PARTIAL step.
  if (vert_coordinate_class==DEPTH_BASED) then
     det_factor = det_param*((rho0*36*pi)**onethird)
     ent_factor = ((rho0*36*pi)**onethird)
  else !PRESSURE_BASED
     det_factor = det_param*((36*pi)**onethird)
     ent_factor = ((36*pi)**onethird)
  endif

  ! Allocate the buffers
  ! Things that dictate the size of the real buffer are:
  ! tracer content (num_prog_tracer), change in mass (1),
  ! change in tracer (num_prog_tracer; only needed for bitwise_reproduction),
  ! velocity(3), position (3),
  ! step size (1), blob time(1),  mass (1), stretching function (2), 
  ! gprime (1), and age (1).
  rea_buff_size = 2*num_prog_tracers+18
  ! Integer buffer is: ijk (3), model_steps (1), hash (1), number (1),
  ! nfrac_steps (1), nsteps (1)
  int_buff_size = 8

  if (bitwise_reproduction) then
     ! History real buffer is: entrainment (num_prog_tracers+1), detrainment (num_prog_tracers+1) 
     ! mass in/out (2).
     hist_rea_buff_size = 2*num_prog_tracers+3 !it is not 2*(num_prog_tracers+1)+2 because the dimension goes from 0 
     ! History integer buffer is: ijk (3)
     hist_int_buff_size = 3
  else
     ! History real buffer is: entrainment (num_prog_tracers+1), detrainment (num_prog_tracers+1) 
     hist_rea_buff_size = 2*num_prog_tracers + 1 !it is not 2*(num_prog_tracers+1) because the dimension goes from 0 
  endif
  call allocate_buffer(Ebuffer_out,  Info%pe_E);  call allocate_buffer(Ebuffer_in,  Info%pe_E)
  call allocate_buffer(Wbuffer_out,  Info%pe_W);  call allocate_buffer(Wbuffer_in,  Info%pe_W)
  call allocate_buffer(Nbuffer_out,  Info%pe_N);  call allocate_buffer(Nbuffer_in,  Info%pe_N)
  call allocate_buffer(Sbuffer_out,  Info%pe_S);  call allocate_buffer(Sbuffer_in,  Info%pe_S)
  call allocate_buffer(NEbuffer_out, Info%pe_NE); call allocate_buffer(NEbuffer_in, Info%pe_NE)
  call allocate_buffer(NWbuffer_out, Info%pe_NW); call allocate_buffer(NWbuffer_in, Info%pe_NW)
  call allocate_buffer(SEbuffer_out, Info%pe_SE); call allocate_buffer(SEbuffer_in, Info%pe_SE)
  call allocate_buffer(SWbuffer_out, Info%pe_SW); call allocate_buffer(SWbuffer_in, Info%pe_SW)

  if (Info%pe_this == mpp_root_pe()) then
     call allocate_buffer(diagbuffer_in, NULL_PE)
  else
     call allocate_buffer(diagbuffer_out, mpp_root_pe())
  endif

  id_clock_rk = mpp_clock_id('(Ocean dyn. bottom blob: RK scheme)', grain=CLOCK_LOOP)

  ! The rest of the initialisation is for  variables that 
  ! are used for calculating formation of blobs.  They are NOT
  ! used for the evolution of blobs.

  ! gravity * fraction of cell participating in overflow / friction
  overflow_factor = grav*blob_overflow_delta/blob_overflow_mu

  ! compute topographic slope arrays for the i-slope and j-slope
  ! slopes are centered on the i-face and j-face of tracer cells 
  allocate(slope_x(isd:ied,jsd:jed))
  allocate(slope_y(isd:ied,jsd:jed))
  slope_x(:,:) = 0.0 
  slope_y(:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        slope_x(i,j) = (Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)
        slope_y(i,j) = (Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)
        slope_x(i,j) = abs(slope_x(i,j))
        slope_y(i,j) = abs(slope_y(i,j))
     enddo
  enddo

  call mpp_update_domains(slope_x(:,:),Dom%domain2d)
  call mpp_update_domains(slope_y(:,:),Dom%domain2d)
  
  ! labels for the four tracer cells surrounding a 
  ! central tracer cell (moving counter-clockwise)
  m=1 ; ip(m)=1  ; jp(m)=0
  m=2 ; ip(m)=0  ; jp(m)=1
  m=3 ; ip(m)=-1 ; jp(m)=0
  m=4 ; ip(m)=0  ; jp(m)=-1

  ! compute directions from an (i,j) point where topography deepens.
  ! these directions may potentially have downslope flow. 
  ! insist that downslope flow occurs only when there are more kmt 
  ! cells in the adjacent column. 
  ! also insist that downslope flow does not involve k=1 cells. 
  allocate (topog_step(isd:ied,jsd:jed,4))
  allocate (topog_slope(isd:ied,jsd:jed,4))
  topog_step(:,:,:)  = 0.0
  topog_slope(:,:,:) = 0.0

  ! figure out which vertical cell in the deep ocean column a blob
  ! will be inserted into when it is created

  do m=1,4
     do j=jsc,jec
        do i=isc,iec
           iip = i+ip(m)
           jjp = j+jp(m)
           ! Make sure that we are in an area we want to form blobs
           if ( blobs_south_of>Grd%yt(iip,jjp) .or. &
                blobs_north_of<Grd%yt(iip,jjp) ) then
              ! check whether downslope flow is possible
              if (kmt(i,j)>1 .and. kmt(iip,jjp)>kmt(i,j)+min_do_levels) then
                 topog_step(i,j,m)=1.0
              endif
           endif
        enddo
     enddo
  enddo

  ! Block out the bipolar fold in order to ensure tracer 
  ! conservation.
  ! NOTE TO SELF: do we need this for the blobs?
  if(jec+Dom%joff==Dom%jeg) topog_step(:,jec,:) = 0.0

  do j=jsc,jec
     do i=isc,iec
        topog_slope(i,j,1) = slope_x(i,j)   !m=1
        topog_slope(i,j,2) = slope_y(i,j)   !m=2
        topog_slope(i,j,3) = slope_x(i-1,j) !m=3
        topog_slope(i,j,4) = slope_y(i,j-1) !m=4
     enddo
  enddo

  call mpp_update_domains(topog_step(:,:,:), Dom%domain2d)
  call mpp_update_domains(topog_slope(:,:,:),Dom%domain2d)

  deallocate(slope_x)
  deallocate(slope_y)

  allocate(vert_factor(isc:iec,jsc:jec,4))
  do m=1,4
     do i=isc,iec
        do j=jsc,jec
           vert_factor(i,j,m) = dtime*grav*(-1 + 1/sqrt( 1+topog_slope(i,j,m)**2 ))
        enddo
     enddo
  enddo

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that allocates memory to a buffer        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine allocate_buffer(buffer, pe)
    type(blob_buffer_type), pointer :: buffer
    integer, intent(in) :: pe

    allocate(buffer)
    buffer%pe = pe
    buffer%numblobs = 0
    if (buffer%pe /= NULL_PE) then
       buffer%size = delta_buffer
       
       allocate(buffer%integer_buff(int_buff_size, delta_buffer))
       allocate(buffer%real_buff(   rea_buff_size, delta_buffer))

       if(bitwise_reproduction) then
          allocate(buffer%history_integer_buff(hist_int_buff_size, 0:total_nsp1, delta_buffer))
          allocate(buffer%history_real_buff( 0:hist_rea_buff_size, 0:total_nsp1, delta_buffer))
       else
          allocate(buffer%history_real_buff(0:hist_rea_buff_size, 1, delta_buffer))
       endif
    endif
  end subroutine allocate_buffer

end subroutine blob_dynamic_bottom_init
! </SUBROUTINE>  NAME="blob_dynamic_bottom_init"

!######################################################################
! <SUBROUTINE NAME="blob_dynamic_bottom_update">
!
! <DESCRIPTION>
! This routine calls the routine to update blob positions.  When
! bitwise_reproduction=.false., it also figures out when to continue
! the integration of blobs that have changed PE's.
! </DESCRIPTION>
!
subroutine blob_dynamic_bottom_update(bottom_head, free, Time, Dens, Thickness, &
                                      T_prog, Ext_mode, L_system, tend_blob,    &
                                      blob_source, u, model_rho, mass_in,       &
                                      mass_out, EL_diag, ngrnd, ndetrn)
  type(ocean_blob_type),                               pointer       :: bottom_head
  type(ocean_blob_type),                               pointer       :: free
  type(ocean_time_type),                               intent(in)    :: Time
  type(ocean_density_type),                            intent(in)    :: Dens
  type(ocean_thickness_type),                          intent(inout) :: Thickness
  type(ocean_prog_tracer_type),                        intent(inout) :: T_prog(:)
  type(ocean_external_mode_type),                      intent(inout) :: Ext_mode
  type(ocean_lagrangian_type),                         intent(inout) :: L_system
  real,dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  real,dimension(isd:ied,jsd:jed),                     intent(inout) :: blob_source
  real, dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),       intent(in)    :: u
  real, dimension(isbd:iebd,jsbd:jebd,1:nk),           intent(in)    :: model_rho
  real, dimension(isd:ied,jsd:jed,1:nk),               intent(inout) :: mass_in
  real, dimension(isd:ied,jsd:jed,1:nk),               intent(inout) :: mass_out
  type(blob_diag_type), dimension(0:num_prog_tracers), intent(inout) :: EL_diag
  integer,                                             intent(inout) :: ngrnd
  integer,                                             intent(inout) :: ndetrn

  type(ocean_blob_type), pointer :: prev, this, next, buffer_head
  integer :: check_buffers
  integer :: i,j,k,n,iit,jjt,mm,tau
  integer :: stdoutunit
  real :: topog, tbigd, tdsq_r(9)

  if (.not. use_this_module) return

  tau=Time%tau

  nullify(free)
  nullify(this)

  ! Clear buffers for sending and receiving bottom blobs
  call clear_buffer(Ebuffer_out);  call clear_buffer( Ebuffer_in)
  call clear_buffer(Wbuffer_out);  call clear_buffer( Wbuffer_in)
  call clear_buffer(Nbuffer_out);  call clear_buffer( Nbuffer_in)
  call clear_buffer(Sbuffer_out);  call clear_buffer( Sbuffer_in)
  call clear_buffer(NEbuffer_out); call clear_buffer(NEbuffer_in)
  call clear_buffer(NWbuffer_out); call clear_buffer(NWbuffer_in)
  call clear_buffer(SEbuffer_out); call clear_buffer(SEbuffer_in)
  call clear_buffer(SWbuffer_out); call clear_buffer(SWbuffer_in)
        
  if (debug_this_module) then
     stdoutunit = stdout()
     write(stdoutunit, '(a)') ' '
     call write_timestamp(Time%model_time)
     write(stdoutunit,'(a)') 'From ocean_blob_dynamic_free_mod'
     write(stdoutunit,'(a)') 'Totals before free dynamic blob update (tau)'
     call E_and_L_totals(L_system,Thickness,T_prog(:),Time%tau)
     write(stdoutunit,'(a)') ' '
  endif

  ! Reset the blob time to zero for the beginning of the Eulerian time step
  if (associated(bottom_head)) then
     this=>bottom_head
     timecycle: do
        this%blob_time = 0.0
        this=>this%next
        if(.not.associated(this)) exit timecycle
     enddo timecycle
  endif

  call dynamic_update(bottom_head, free, Time, Dens, Thickness, T_prog, &
                      Ext_mode, L_system, tend_blob, blob_source, u,    &
                      model_rho, mass_in, mass_out, EL_diag(:), ngrnd, ndetrn)

  if (bitwise_reproduction) then
     ! Send, receive and unpack buffers for free blobs
     call send_buffer( Ebuffer_out)
     call send_buffer( Wbuffer_out)
     call send_buffer( Sbuffer_out)
     call send_buffer( Nbuffer_out)
     call send_buffer(NEbuffer_out)
     call send_buffer(NWbuffer_out)
     call send_buffer(SEbuffer_out)
     call send_buffer(SWbuffer_out)
     
     call mpp_sync_self()

     call receive_buffer( Wbuffer_in); call unpackbuffer(Time,  Wbuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer( Ebuffer_in); call unpackbuffer(Time,  Ebuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer( Nbuffer_in); call unpackbuffer(Time,  Nbuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer( Sbuffer_in); call unpackbuffer(Time,  Sbuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer(SWbuffer_in); call unpackbuffer(Time, SWbuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer(SEbuffer_in); call unpackbuffer(Time, SEbuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer(NWbuffer_in); call unpackbuffer(Time, NWbuffer_in, bottom_head, Dens, Ext_mode, Thickness)
     call receive_buffer(NEbuffer_in); call unpackbuffer(Time, NEbuffer_in, bottom_head, Dens, Ext_mode, Thickness)

     call mpp_sync_self()

  else
     ! When bitwise reproduction is not necessary, blobs can be packed into 
     ! buffers and sent to other PE's before they have finished a full E system
     ! step, i.e. this%blob_time < dtime.
     ! So, here, we unpack the blobs and continue to evolve them on the new PE.
     ! We need to take into account that a blob may traverse more than one PE
     ! in a given E system step, e.g.
     ! ------------------
     ! |       |        |
     ! | PE=2  +a  PE=3 |
     ! |      /|        |
     ! |     b |        |
     ! ------+-----------
     ! |    /  |        |
     ! |   c   |        |
     ! |       |        |
     ! |  PE=0 |   PE=1 |
     ! ------------------
     ! So, we need to keep on checking buffers and evolving blobs until all received 
     ! buffers are empty

     ! Send, receive and unpack buffers for bottom blobs
     call send_buffer( Ebuffer_out) 
     call send_buffer( Wbuffer_out)
     call send_buffer( Sbuffer_out)
     call send_buffer( Nbuffer_out)
     call send_buffer(NEbuffer_out)
     call send_buffer(NWbuffer_out)
     call send_buffer(SEbuffer_out)
     call send_buffer(SWbuffer_out)
     
     call mpp_sync_self()
     
     buffer_head => NULL()
     call receive_buffer( Wbuffer_in) 
     call unpackbuffer(Time,  Wbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer( Ebuffer_in)
     call unpackbuffer(Time,  Ebuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer( Nbuffer_in)
     call unpackbuffer(Time,  Nbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer( Sbuffer_in)
     call unpackbuffer(Time,  Sbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer(SWbuffer_in)
     call unpackbuffer(Time, SWbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer(SEbuffer_in)
     call unpackbuffer(Time, SEbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer(NWbuffer_in)
     call unpackbuffer(Time, NWbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
     call receive_buffer(NEbuffer_in) 
     call unpackbuffer(Time, NEbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)

     call mpp_sync_self()

     ! Clear buffers for sending and receiving free blobs
     call clear_buffer(Ebuffer_out);  call clear_buffer( Ebuffer_in)
     call clear_buffer(Wbuffer_out);  call clear_buffer( Wbuffer_in)
     call clear_buffer(Nbuffer_out);  call clear_buffer( Nbuffer_in)
     call clear_buffer(Sbuffer_out);  call clear_buffer( Sbuffer_in)
     call clear_buffer(NEbuffer_out); call clear_buffer(NEbuffer_in)
     call clear_buffer(NWbuffer_out); call clear_buffer(NWbuffer_in)
     call clear_buffer(SEbuffer_out); call clear_buffer(SEbuffer_in)
     call clear_buffer(SWbuffer_out); call clear_buffer(SWbuffer_in)

     check_buffers = 0
     this=>buffer_head
     if (associated(this)) then
        buff_cycle0: do

           if(this%blob_time == dtime) then
              ! Adjust the geodepth so we are sitting on topography
              topog = 0.0
              i=this%i; j=this%j; k=this%k
              call interp_tcoeff(i,j,(/this%h1,this%h2/),this%lon,this%lat,tdsq_r(:))
              tbigd = 0.0
              do mm=1,Info%tidx(0,i,j)
                 iit=i+Info%it(Info%tidx(mm,i,j))
                 jjt=j+Info%jt(Info%tidx(mm,i,j))
                 topog = topog + ht(iit,jjt)*tdsq_r(mm)
                 tbigd = tbigd + tdsq_r(mm)
              enddo
              tbigd = 1.0/tbigd
              
              this%geodepth =  topog*tbigd
              this%depth    =  this%geodepth + Ext_mode%eta_t(i,j,tau)
              if (vert_coordinate_class==DEPTH_BASED) then
                 this%st = -Thickness%depth_swt(i,j,k)
              else !PRESSURE_BASED
                 this%st = Thickness%depth_swt(i,j,k)
              endif
              
              call unlink_blob(this, buffer_head, prev, next)
              call insert_blob(this, bottom_head)
              this=>next
           else
              check_buffers = check_buffers + 1
              this=>this%next
           endif
           if (.not. associated(this)) exit buff_cycle0
        enddo buff_cycle0
     endif

     call mpp_sum(check_buffers)
     
     do while (check_buffers>0) 
        call dynamic_update(buffer_head, free, Time, Dens, Thickness, &
                            T_prog(:), Ext_mode, L_system, tend_blob, &
                            blob_source, u, model_rho, mass_in,       &
                            mass_out, EL_diag(:), ngrnd, ndetrn)

        ! Send, receive and unpack buffers for bottom blobs
        call send_buffer( Ebuffer_out) 
        call send_buffer( Wbuffer_out)
        call send_buffer( Sbuffer_out)
        call send_buffer( Nbuffer_out)
        call send_buffer(NEbuffer_out)
        call send_buffer(NWbuffer_out)
        call send_buffer(SEbuffer_out)
        call send_buffer(SWbuffer_out)

        call mpp_sync_self()
  
        call receive_buffer( Wbuffer_in) 
        call unpackbuffer(Time,  Wbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer( Ebuffer_in)
        call unpackbuffer(Time,  Ebuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer( Nbuffer_in)
        call unpackbuffer(Time,  Nbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer( Sbuffer_in)
        call unpackbuffer(Time,  Sbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer(SWbuffer_in)
        call unpackbuffer(Time, SWbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer(SEbuffer_in)
        call unpackbuffer(Time, SEbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer(NWbuffer_in)
        call unpackbuffer(Time, NWbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)
        call receive_buffer(NEbuffer_in) 
        call unpackbuffer(Time, NEbuffer_in, buffer_head, Dens, Ext_mode, Thickness, L_system, EL_diag, blob_source, tend_blob, mass_in)

        call mpp_sync_self()

        ! Clear buffers for sending and receiving free blobs
        call clear_buffer(Ebuffer_out);  call clear_buffer( Ebuffer_in)
        call clear_buffer(Wbuffer_out);  call clear_buffer( Wbuffer_in)
        call clear_buffer(Nbuffer_out);  call clear_buffer( Nbuffer_in)
        call clear_buffer(Sbuffer_out);  call clear_buffer( Sbuffer_in)
        call clear_buffer(NEbuffer_out); call clear_buffer(NEbuffer_in)
        call clear_buffer(NWbuffer_out); call clear_buffer(NWbuffer_in)
        call clear_buffer(SEbuffer_out); call clear_buffer(SEbuffer_in)
        call clear_buffer(SWbuffer_out); call clear_buffer(SWbuffer_in)

        check_buffers = 0
        
        if (associated(buffer_head)) then
           this => buffer_head
           timecheck: do
              if(this%blob_time == dtime) then
                 ! Adjust the geodepth so we are sitting on topography
                 topog = 0.0
                 i=this%i; j=this%j; k=this%k
                 call interp_tcoeff(i,j,(/this%h1,this%h2/),this%lon,this%lat,tdsq_r(:))
                 tbigd = 0.0
                 do mm=1,Info%tidx(0,i,j)
                    iit=i+Info%it(Info%tidx(mm,i,j))
                    jjt=j+Info%jt(Info%tidx(mm,i,j))
                    topog = topog + ht(iit,jjt)*tdsq_r(mm)
                    tbigd = tbigd + tdsq_r(mm)
                 enddo
                 tbigd = 1.0/tbigd
                 
                 this%geodepth =  topog*tbigd
                 this%depth    =  this%geodepth + Ext_mode%eta_t(i,j,tau)
                 if (vert_coordinate_class==DEPTH_BASED) then
                    this%st = -Thickness%depth_swt(i,j,k)
                 else !PRESSURE_BASED
                    this%st = Thickness%depth_swt(i,j,k)
                 endif
 
                 call unlink_blob(this, buffer_head, prev, next)
                 call insert_blob(this, bottom_head)
                 this=>next
              else
                 check_buffers = check_buffers + 1
                 this=>this%next
              endif
              if(.not. associated(this)) exit timecheck
           enddo timecheck
        endif
        call mpp_sum(check_buffers)

     enddo
  endif!bitwise_reproduction


  ! Check if a separating blob has also changed processors.  If it has, pack it into
  ! a buffer, send the buffer, receive other buffers and unpack them.

  ! Clear buffers for sending and receiving free blobs
  call clear_buffer(Ebuffer_out);  call clear_buffer( Ebuffer_in)
  call clear_buffer(Wbuffer_out);  call clear_buffer( Wbuffer_in)
  call clear_buffer(Nbuffer_out);  call clear_buffer( Nbuffer_in)
  call clear_buffer(Sbuffer_out);  call clear_buffer( Sbuffer_in)
  call clear_buffer(NEbuffer_out); call clear_buffer(NEbuffer_in)
  call clear_buffer(NWbuffer_out); call clear_buffer(NWbuffer_in)
  call clear_buffer(SEbuffer_out); call clear_buffer(SEbuffer_in)
  call clear_buffer(SWbuffer_out); call clear_buffer(SWbuffer_in)
        
  if (associated(free)) then
     this=>free
     freecycle: do
        i=this%i
        j=this%j
        k=this%k
        if     (i>iec .and. j>jec) then
           call packfreebuffer(NEbuffer_out)
        elseif (i<isc .and. j>jec) then
           call packfreebuffer(NWbuffer_out)
        elseif (i>iec .and. j<jsc) then
           call packfreebuffer(SEbuffer_out)
        elseif (i<isc .and. j<jsc) then
           call packfreebuffer(SWbuffer_out)
        elseif (i>iec)             then
           call packfreebuffer(Ebuffer_out)
        elseif (i<isc)             then
           call packfreebuffer(Wbuffer_out)
        elseif (j>jec)             then
           call packfreebuffer(Nbuffer_out)
        elseif (j<jsc)             then
           call packfreebuffer(Sbuffer_out)
        endif
        this=>this%next
        if(.not.associated(this)) exit freecycle
     enddo freecycle
  endif

  ! Send buffers for bottom blobs that have become free blobs
  call send_buffer( Ebuffer_out)
  call send_buffer( Wbuffer_out)
  call send_buffer( Sbuffer_out)
  call send_buffer( Nbuffer_out)
  call send_buffer(NEbuffer_out)
  call send_buffer(NWbuffer_out)
  call send_buffer(SEbuffer_out)
  call send_buffer(SWbuffer_out)

  call mpp_sync_self()
  
  ! Receive and unpack bottom blobs that have become free blobs
  call receive_buffer( Ebuffer_in); call unpackbuffer(Time,  Ebuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer( Wbuffer_in); call unpackbuffer(Time,  Wbuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer( Nbuffer_in); call unpackbuffer(Time,  Nbuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer( Sbuffer_in); call unpackbuffer(Time,  Sbuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer(NEbuffer_in); call unpackbuffer(Time, NEbuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer(NWbuffer_in); call unpackbuffer(Time, NWbuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer(SEbuffer_in); call unpackbuffer(Time, SEbuffer_in, free, Dens, Ext_mode, Thickness)
  call receive_buffer(SWbuffer_in); call unpackbuffer(Time, SWbuffer_in, free, Dens, Ext_mode, Thickness)

  call mpp_sync_self()

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that packs the buffer of bottom blobs    !
! that are transferring to free blobs and have crossed a processor     !
! boundary                                                             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine packfreebuffer(buffer)
    type(blob_buffer_type), pointer :: buffer
    real, dimension(0:num_prog_tracers) :: entrainment, detrainment
    integer :: s

    if (buffer%pe == Info%pe_this) then
       ! We have a cyclic grid and the blob is just going to
       ! end up back on the same processor.  So, we just adjust
       ! the (i,j) values for the blob.
       call check_cyclic(this, this%i, this%j, .true.)

       if (bitwise_reproduction) then
          do s=0,this%nfrac_steps
             call check_cyclic(this, this%di(s), this%dj(s), .false.)
          enddo
       endif

    else !going to a different PE
       if (bitwise_reproduction) then
          call packbuffer(this,buffer)
       else
          entrainment(:) = 0.0 
          detrainment(:) = 0.0 
          call packbuffer(this,buffer, entrainment(:), detrainment(:))
       endif
       this%mass=0.
       do n=1,num_prog_tracers
          this%tracer(n) = 0.
       enddo

    endif
  end subroutine packfreebuffer
  

end subroutine blob_dynamic_bottom_update
! </SUBROUTINE>  NAME="blob_dynamic_bottom_update"

!######################################################################
! <SUBROUTINE NAME="dynamic_update">
!
! <DESCRIPTION>
! This routine contains the RK scheme used to integrate the position 
! and velocity of blobs.  It also does many checks for (and 
! subsequently handles) things like grounding of blobs, blobs going to
! different PEs, blobs that separate from topography, blobs that 
! detrain to less than small_mass and blobs going outside the compute 
! domain.  
!
! It also does the interpolation of E system variables to a blob.
! </DESCRIPTION>
!
subroutine dynamic_update(head, free, Time, Dens, Thickness,    &
                          T_prog, Ext_mode, L_system,           &
                          tend_blob, blob_source, u, model_rho, &
                          mass_in, mass_out, EL_diag, ngrnd, ndetrn)
  type(ocean_blob_type),                               pointer       :: head
  type(ocean_blob_type),                               pointer       :: free
  type(ocean_time_type),                               intent(in)    :: Time
  type(ocean_density_type),                            intent(in)    :: Dens
  type(ocean_thickness_type),                          intent(inout) :: Thickness
  type(ocean_prog_tracer_type),                        intent(inout) :: T_prog(:)
  type(ocean_external_mode_type),                      intent(inout) :: Ext_mode   
  type(ocean_lagrangian_type),                         intent(inout) :: L_system
  real,dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  real,dimension(isd:ied,jsd:jed),                     intent(inout) :: blob_source
  real, dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),       intent(in)    :: u
  real, dimension(isbd:iebd,jsbd:jebd,1:nk),           intent(in)    :: model_rho
  real, dimension(isd:ied,jsd:jed,1:nk),               intent(inout) :: mass_in
  real, dimension(isd:ied,jsd:jed,1:nk),               intent(inout) :: mass_out
  type(blob_diag_type), dimension(0:num_prog_tracers), intent(inout) :: EL_diag
  integer,                                             intent(inout) :: ngrnd
  integer,                                             intent(inout) :: ndetrn

  type(ocean_blob_type), pointer :: this => NULL()
  type(ocean_blob_type), pointer :: prev => NULL()
  type(ocean_blob_type), pointer :: next => NULL()
  real,    dimension(1:6,0:ns) :: V
  integer, dimension(3) :: old_ijk
  logical, dimension(2) :: off
  real,    dimension(6) :: Xn, update, Xnp1, Xnp1_hat
  real,    dimension(3) :: old_lld, vel
  real,    dimension(2) :: h, old_h
  real,    dimension(num_prog_tracers) :: tracer, field
  real,    dimension(0:num_prog_tracers) :: entrainment, detrainment
  real, dimension(9) :: tdsq_r
  real, dimension(4) :: udsq_r
  logical :: go_e, go_ne, go_n, go_nw, go_w, go_sw, go_s, go_se
  integer :: total_blobs, leaving_pe, tfer_free, detrn_zero, move_lateral, ngrounded
  integer :: i, j, k, m, mm, n, r, s, ii
  integer :: iiu,jju,iit,jjt,ku
  integer :: n_frac_steps, tau
  integer :: stdoutunit
  logical :: change_pe, reached_end, advance_blob, accept_step, grounded
  real :: blob_time, tstep, lon, lat, geodepth, rhoL, rhoLr
  real :: rho, rhor, rhoE, volume
  real :: old_tstep, trunc, hstar, tstar, mass
  real :: Hx, Hy, bstress_hr, absv2, absv, uE, vE, absrelv2, absrelv
  real :: height, heightr, f, richardson, dmdt_ent, dmdt_det, ent, det, tstress
  real :: ent_hr, red_grav, rgd
  real :: ubigd, tbigd, denom, topog

  nullify(this)

  tau = Time%tau

  ! zero some diagnostics
  total_blobs  = 0
  leaving_pe   = 0
  tfer_free    = 0
  detrn_zero   = 0
  move_lateral = 0
  ngrounded    = 0

  if (associated(head)) then
     this=>head

     blobcycle: do
        ! time relative to the beginning of the time step
        blob_time = this%blob_time
        tstep     = this%step

        this%dmass      = 0.0
        this%dtracer(:) = 0.0
        if (bitwise_reproduction) then
           ! initialise the values of the history arrays
           this%di(:) = -1
           this%dj(:) = -1
           this%dk(:) = -1
           this%nfrac_steps      = 0
           this%entrainment(:,:) = 0.0 
           this%detrainment(:,:) = 0.0 
           this%mass_in(:)       = 0.0
           this%mass_out(:)      = 0.0
           
           this%di(0) = this%i
           this%dj(0) = this%j
           this%dk(0) = this%k
        endif

        ! Make adjustments to the lon-lat and i-j (if necessary)
        ! on a cyclic or tripolar grid.
        call check_cyclic(this, this%i, this%j, .true.)

        ! We don't want to alter the blob properties until 
        ! we have completed the full fractional step.  So, we save
        ! the variables into local variables in case the step is
        ! rejected and we need to try again with a different value
        ! for tstep.
        
        ! These variables get updated from one RK sub-step to
        ! the next, so, we need to save copies of them in case
        ! the step is rejected.
        i        = this%i - Dom%ioff
        j        = this%j - Dom%joff
        k        = this%k
        old_ijk  = (/ i, j, k/)
        lon      = this%lon
        lat      = this%lat
        geodepth = this%geodepth
        old_lld  = (/lon, lat, geodepth/)
        h(:)     = (/ this%h1, this%h2 /)
        old_h    = h
        
        ! These variables are only updated at the end
        ! of all the RK sub-steps, so, we do not need
        ! to save copies
        ent    = this%ent
        det    = this%det
        richardson = this%richardson
        vel(:) = this%v(:)
        mass   = this%mass
        volume = this%volume
        rhoL   = this%density
        rhoLr  = this%densityr
        if (vert_coordinate_class==DEPTH_BASED) then
           rho  = rho0
           rhor = rho0r
        else !PRESSURE_BASED
           rho  = rhoL
           rhor = rhoLr
        endif
        do n=1,num_prog_tracers
           tracer(n) = this%tracer(n)
           field(n)  = this%field(n)
        enddo
        height  = this%height
        heightr = 1.0/this%height
        n_frac_steps = 0
        
        reached_end  = .false.
        change_pe    = .false.
        advance_blob = .true.
        off(:)       = .false.
        grounded     = .false.
        
        ! Calculate the horizontal interpolation coefficients
        call interp_tcoeff(i,j,h(:),lon,lat,tdsq_r(:))
        call interp_ucoeff(i,j,h(:),lon,lat,udsq_r(:))        

        partialstep: do while (.not. reached_end)
           accept_step = .false.
           
           ! Interpolate the E system variables to the blobs
           ! Note, we treat the stretching function separately

           ! U-grid variables
           ! We want to treat the 
           ubigd = 0.0
           uE=0.0; vE=0.0
           Hx=0.0; Hy=0.0
           do mm=1,Info%uidx(0,i,j)
              iiu=i+Info%iu(Info%uidx(mm,i,j))
              jju=j+Info%ju(Info%uidx(mm,i,j))
              ku =kmu(iiu,jju)
              Hx = Hx + dht_dx(iiu,jju)*udsq_r(mm)
              Hy = Hy + dht_dy(iiu,jju)*udsq_r(mm)
              ! Because we interpolate from surrounding bottom cells, we make
              ! use of kmu.  kmu can have values of k=0 (which means, a land
              ! point).  When interpolating velocities from the U grid, we wish
              ! for them to have a zero value over land points.
              if (ku>0) then
                 uE = uE + umask(iiu,jju,ku)*u(iiu,jju,ku,1)*udsq_r(mm)
                 vE = vE + umask(iiu,jju,ku)*u(iiu,jju,ku,2)*udsq_r(mm)
              endif
              ubigd = ubigd + udsq_r(mm)
           enddo
           ubigd = 1.0/ubigd

           Hx = Hx*ubigd
           Hy = Hy*ubigd
           uE = uE*ubigd
           vE = vE*ubigd

           ! Put rhoE to the model density (rather than interpolate) so 
           ! free blobs are not prematurely created.
           rhoE = model_rho(i,j,k)

           ! speed of the blob
           absv2  = vel(1)**2 + vel(2)**2 + vel(3)**2
           absv   = sqrt(absv2)

           ! relative speed of the blob (to the E system)
           absrelv2 = (vel(1)-uE)**2 + (vel(2)-vE)**2 + vel(3)**2
           absrelv  = sqrt(absrelv2)

           m=0
           trystep: do while (.not. accept_step)
              prev => NULL()
              next => NULL()

              ! Load the old variables into some local variables that we will change.
              ! We want to keep the old variables in case the step is rejected.
              i = old_ijk(1)
              j = old_ijk(2)
              k = old_ijk(3)
              h(:) = old_h(:)
              lon  = old_lld(1)
              lat  = old_lld(2)
              
              ! Horizontal distances are calculated relative to the inital 
              ! position for this partial step.  Vertical distance is calculated
              ! relative to the initial geodepth of the step
              Xn(1) =  0.         ! x-position
              Xn(2) =  vel(1)     ! u-velocity
              Xn(3) =  0.         ! y-position
              Xn(4) =  vel(2)     ! v-velocity
              Xn(5) = -old_lld(3) ! z-position (relative to z=0)
              Xn(6) =  vel(3)     ! w-velocity

              ! Temporary variables for the RK scheme
              Xnp1(:)     = Xn(:)
              Xnp1_hat(:) = Xn(:)
              update(:)   = Xn(:)

              ! Temporary variables for the entrainment/detrainment history
              entrainment(:) = 0.
              detrainment(:) = 0.

              ! Friction from entrainment
              ent_hr = ent*heightr
              
              ! The bottom stress.  We omit a rho, because it cancels below
              bstress_hr = drag*absv*heightr
              tstress    = bstress_hr + ent_hr
              
              ! rotation
              f  = two_omega*sin(deg_to_rad*lat)
              
              ! The pseudo-non-hydrostatic (reduced gravity) and normal force terms
              denom    = 1.0/sqrt(Hx**2 + Hy**2 + 1)
              red_grav = grav*(rhoL - rhoE)*rhor
              rgd      = red_grav * denom

              ! Cycle through the RK sub-steps
              call mpp_clock_begin(id_clock_rk)
              do s=0,ns
                 ! Solve the equations of motion
                 V(1,s) =   update(2)                                           !xdot
                 V(2,s) = - tstress*update(2) + update(4)*f + rgd*Hx + uE*ent_hr!xddot
                 
                 V(3,s) =   update(4)                                           !ydot
                 V(4,s) = - update(2)*f - tstress*update(4) + rgd*Hy + vE*ent_hr!yddot
                 
                 V(5,s) =   update(6)                                           !zdot
                 V(6,s) = - tstress*update(6)               + rgd    - red_grav !zddot

                 ! calculate new update
                 if(s<ns) then
                    update(:) = Xn(:)
                    do r=0,s
                       update(:) = update(:) + tstep*beta(s+1,r)*V(:,r)
                    enddo
                 endif

                 ! The high and low order estimates
                 Xnp1(:)     = Xnp1(:)     + tstep*hgamma(s)*V(:,s)
                 Xnp1_hat(:) = Xnp1_hat(:) + tstep*lgamma(s)*V(:,s)

              enddo!s
              call mpp_clock_end(id_clock_rk)

              ! Estimate the error using the difference in the high and low order schemes.
              ii    = maxloc(abs((Xnp1(:) - Xnp1_hat(:))/(Xnp1(:)+epsln)),1)
              trunc = abs(Xnp1(ii) - Xnp1_hat(ii) + epsln)
              tstar = rel_error*abs(Xnp1(ii))
              hstar = tstep * (safety_factor*tstar/trunc)**(orderp1_r)

              ! Save the old time step (if this step fails, we will need it again)
              old_tstep = tstep

              ! Reject the step if the accuracy condition is not met. Redo the 
              ! trystep loop with an adjusted step.
              ! There are some caveats, listed below.
              if(trunc > tstar) then 
                 accept_step = .false.
                 if (hstar > minstep) then 
                    !If the suggested step is not smaller than the minimum step, 
                    !retry with the smaller step.
                    tstep = min(hstar, dtime - blob_time)
                 else
                    if (tstep<=minstep) then 
                       !If the step is smaller than the minimum step and 
                       !it is still rejected, we accept the step anyway.
                       accept_step = .true.
                    else
                       ! Otherwise, set the step to the minimum step size.
                       accept_step = .false.
                       tstep = min(minstep, dtime - blob_time)
                    endif
                 endif
                 
              else !Accept the step
                 accept_step = .true.
                 ! We want to impose some restrictions on the local extrapolation to ensure
                 ! we maintain stability and minimise the number of rejected steps:
                 ! 1/ only let the next step increase in size if the previous step was 
                 ! never rejected,
                 ! 2/ we don't let the next step be any larger than twice this step.
                 if (m==0 .and. hstar>=tstep) then
                    tstep = min(2*tstep, hstar)
                 elseif (hstar<tstep) then
                    tstep = max(minstep, hstar)
                 endif
                 ! 3/ We do not want tstep to be greater than dtime (the E system time step)
                 tstep = min(tstep, dtime)

                 ! There may be more adjustments to tstep later on, depending on
                 ! how close we are to reaching the Eulerian time step.
              endif

              if (accept_step) then
                 ! Do some rudimentary stability checking
                 if (abs(Xnp1(2))>large_speed .or. abs(Xnp1(4))>large_speed .or. abs(Xnp1(6))>large_speed) then
                    if (tstep<=minstep) then
                       call mpp_error(WARNING,&
                            '==>Warning in ocean_blob_dynamic_bottom_mod (dynamic_update): It looks like a blob is becoming unstable'//&
                            ' and the timestep is already as small as is permitted.  Suggest reducing minstep in the namelist.')
                    endif
                 endif

                 ! Update the grid stretching values
                 tbigd = 0.0
                 h(:)  = 0.0
                 do mm=1,Info%tidx(0,i,j)
                    iit=i+Info%it(Info%tidx(mm,i,j))
                    jjt=j+Info%jt(Info%tidx(mm,i,j))
                    ! Stretching function is defined over land points, so, we do 
                    ! not mask out land points.
                    h(:) = h(:) + Info%ht(iit,jjt,:)*tdsq_r(mm)
                    tbigd = tbigd + tdsq_r(mm)
                 enddo

                 tbigd = 1.0/tbigd
                 h(:) = h(:)*tbigd

                 n_frac_steps = n_frac_steps+1

                 ! Entrainment and detrainment calculations are based on the properties of the blob
                 ! at the beginning of this blob step.
                 richardson = grav*height*abs(rhoL-rhoE)/(rho*absrelv2 + epsln)
                 if (richardson <= critical_richardson) then
                    ent = absrelv*(0.08 - 0.1*richardson)/(1.+5.*richardson)
                 else
                    ent = 0.
                 endif
                 
                 det = -det_factor/abs(rhoL - rhoE + epsln)
                 
                 if (vert_coordinate_class==DEPTH_BASED) then
                    !det_factor = Gamma * (rho0*36pi)**1/3
                    dmdt_det = det * (mass**twothirds)
                    !ent_factor = (rho0*36*pi)**onethird
                    dmdt_ent = ent * ent_factor * (mass**onethird)
                    
                 else !PRESSURE_BASED
                    !det_factor = Gamma * (36pi)**1/3
                    dmdt_det = det * (mass**twothirds) * (rhoL**onethird)
                    !ent_factor = (36*pi)**onethird
                    dmdt_ent = ent * rhoE * ent_factor*(volume**onethird)
                    
                 endif

                 ! Make sure we do not exceed the maximum detrainment
                 dmdt_det = sign(min(abs(dmdt_det),abs(max_detrainment)),dmdt_det)

                 detrainment(0) = dmdt_det*old_tstep
                 entrainment(0) = dmdt_ent*old_tstep
                 do n=1,num_prog_tracers
                    ! Detrainment
                    detrainment(n) = detrainment(0)*field(n)
                    ! Entrainment
                    entrainment(n) = entrainment(0)*T_prog(n)%field(i,j,k,tau)
                 enddo
                 
                 ! Check and see if the blob will detrain to zero mass during this partial time step.
                 ! If it does, just dump all its properties in the original cell.  We also check to 
                 ! see if a blob would fully detrain (without entrainment).  It has been obvserved
                 ! that instabilities can ensue with a realistic EOS if a blob "head bangs" 
                 ! around a zero tracer content.
                 if ((entrainment(0) + detrainment(0) + mass) < small_mass .or. mass < abs(detrainment(0))) then
                    entrainment(0) = 0.0
                    detrainment(0) = -mass
                    do n=1,num_prog_tracers
                       entrainment(n) = 0.0
                       detrainment(n) = -tracer(n)
                    enddo
                    detrn_zero = detrn_zero + 1

                    ! If a blob detrains to less than small mass, we no longer
                    ! need to calculate its trajectory.  So that it stops being 
                    ! processed, we set the time to the end of full step.
                    old_tstep = 0.0
                    blob_time = dtime

                 else
                    ! Check if we have changed horizontal cells
                    call check_ijcell(Xnp1(1),Xnp1(3),i,j,h,old_lld(1:2),lon,lat,off(1:2))

                 endif

                 if (Grd%kmt(i,j) == 0) then
                    ! If the blob has grounded, we put it back in its previous cell
                    ! It will have all its properties returned to the E system later
                    ! Here, grounded means a blob has moved laterally into a zero 
                    ! depth (i.e. land) column, NOT that it has interacted with topography.
                    off(1:2) = .false.
                    grounded = .true.
                    ngrounded = ngrounded + 1
                    i = old_ijk(1)
                    j = old_ijk(2)
                    lon = old_lld(1)
                    lat = old_lld(2)
                 else
                    k = Grd%kmt(i,j)
                 endif

                 if (any(off(1:2))) change_pe = .true.

                 if (change_pe .and. .not. bitwise_reproduction) then
                    ! Check to see if the blob has left this pe
                    ! If we are not enforcing bitwise reproduction, we pack 
                    ! the blob into a buffer immediately, as we do not need
                    ! to worry about saving the histories and processing the
                    ! blob histories in order.
                    Ext_mode%conv_blob(old_ijk(1),old_ijk(2)) = &
                         Ext_mode%conv_blob(old_ijk(1),old_ijk(2)) - mass
                    
                    L_system%conv_blob(old_ijk(1),old_ijk(2),old_ijk(3)) = &
                         L_system%conv_blob(old_ijk(1),old_ijk(2),old_ijk(3)) - mass

                    mass_out(old_ijk(1),old_ijk(2),old_ijk(3)) = &
                         mass_out(old_ijk(1),old_ijk(2),old_ijk(3)) + mass

                    ! Save some variables
                    do n=1,num_prog_tracers
                       this%tracer(n) = tracer(n)
                    enddo
                    this%v(1)      = Xnp1(2)
                    this%v(2)      = Xnp1(4)
                    this%v(3)      = Xnp1(6)
                    this%lon       = lon
                    this%lat       = lat
                    this%geodepth  = -Xnp1(5)
                    this%blob_time = blob_time + old_tstep
                    if (this%blob_time == dtime) then
                       this%ent    = ent
                       this%det    = det
                       this%richardson = richardson
                       this%age    = this%age + dtime_yr
                       this%gprime = grav*(rhoE-rhoL)/rhoE
                       this%step = max(tstep,minstep)
                    elseif ((this%blob_time + 1.1*tstep) > dtime) then
                       this%step = dtime - this%blob_time
                    else
                       this%step      = max(tstep,minstep)
                    endif
                    this%mass = mass

                    this%i = i
                    this%j = j
                    this%k = k

                    if     (i>iec .and. j>jec) then
                       call packbuffer(this, NEbuffer_out, entrainment(:), detrainment(:))
                    elseif (i<isc .and. j>jec) then
                       call packbuffer(this, NWbuffer_out, entrainment(:), detrainment(:))
                    elseif (i>iec .and. j<jsc) then
                       call packbuffer(this, SEbuffer_out, entrainment(:), detrainment(:))
                    elseif (i<isc .and. j<jsc) then
                       call packbuffer(this, SWbuffer_out, entrainment(:), detrainment(:))
                    elseif (i>iec)             then
                       call packbuffer(this,  Ebuffer_out, entrainment(:), detrainment(:))
                    elseif (i<isc)             then
                       call packbuffer(this,  Wbuffer_out, entrainment(:), detrainment(:))
                    elseif (j>jec)             then
                       call packbuffer(this,  Nbuffer_out, entrainment(:), detrainment(:))
                    elseif (j<jsc)             then
                       call packbuffer(this,  Sbuffer_out, entrainment(:), detrainment(:))
                    endif

                    move_lateral = move_lateral + 1
                    leaving_pe   = leaving_pe   + 1

                    call unlink_blob(this, head, prev, next)
                    call free_blob_memory(this)
                    advance_blob=.false.

                    exit partialstep
                 endif!change_pe .and. .not. bitwise_reproduction

                 if (grounded) then
                    
                    this%nfrac_steps = n_frac_steps

                    old_tstep = 0.0
                    blob_time = dtime

                    if (bitwise_reproduction) then
                       this%di(n_frac_steps) = i
                       this%dj(n_frac_steps) = j
                       this%dk(n_frac_steps) = k
                       this%dmass = -mass
                       this%mass  = 0.0
                       mass       = 0.0
                       do n=1,num_prog_tracers
                          this%dtracer(n) = -tracer(n)
                          this%tracer(n)  = 0.0
                          tracer(n)       = 0.0
                       enddo

                       exit trystep
                    else
                       blob_source(this%i,this%j) = blob_source(this%i,this%j) + mass
                       this%mass        = 0.0
                       mass             = 0.0
                       do n=1,num_prog_tracers
                          tend_blob(n,this%i,this%j,this%k) = tend_blob(n,this%i,this%j,this%k) + tracer(n)
                          this%tracer(n)     = 0.0
                          tracer(n)          = 0.0
                       enddo

                       ! Diagnostics 
                       EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + mass 
                       do n=1,num_prog_tracers 
                          EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + tracer(n) 
                       enddo
                       exit trystep
                    endif
                 endif!grounded

                 if (bitwise_reproduction) then

                    ! Update the intermediate values of variables
                    ! If we detect that a blob changes cells, save the mass in and out
                    ! Note, the convention is that both are positive and we take care of
                    ! the sign when calculating the convergence.
                    if (old_ijk(1)/=i .or. old_ijk(2)/=j .or. old_ijk(3)/=k) then
                       this%mass_out(n_frac_steps-1) = mass
                       this%mass_in( n_frac_steps  ) = mass
                    endif
                    
                    ! Save the indices to history too (for divergence and ent/detrainment)
                    this%di(n_frac_steps) = i + Dom%ioff
                    this%dj(n_frac_steps) = j + Dom%joff
                    this%dk(n_frac_steps) = k
                 
                    ! save mass
                    mass                     = mass + entrainment(0) + detrainment(0)
                    this%entrainment(n_frac_steps,0) = entrainment(0)
                    this%detrainment(n_frac_steps,0) = detrainment(0)
                    
                    ! save tracer
                    do n=1,num_prog_tracers
                       tracer(n)                        = tracer(n) + entrainment(n) + detrainment(n)
                       this%entrainment(n_frac_steps,n) = entrainment(n)
                       this%detrainment(n_frac_steps,n) = detrainment(n)
                    enddo

                 else !not bitwise reproduction

                    ! If we detect that a blob changes cells, save the mass into the new
                    ! water column and and out of the old water column.
                    ! save mass
                    if (old_ijk(1)/=i .or. old_ijk(2)/=j .or. old_ijk(3)/=k) then
                       Ext_mode%conv_blob(old_ijk(1),old_ijk(2)) = &
                            Ext_mode%conv_blob(old_ijk(1),old_ijk(2)) - mass
                       Ext_mode%conv_blob(i,j) = &
                            Ext_mode%conv_blob(i,j) + mass

                       L_system%conv_blob(old_ijk(1),old_ijk(2),old_ijk(3)) = &
                            L_system%conv_blob(old_ijk(1),old_ijk(2),old_ijk(3)) - mass
                       L_system%conv_blob(i,j,k) = &
                            L_system%conv_blob(i,j,k) + mass

                       mass_out(old_ijk(1),old_ijk(2),old_ijk(3)) = mass_out(old_ijk(1),old_ijk(2),old_ijk(3)) + mass
                       mass_in(          i,         j,         k) = mass_in(          i,         j,         k) + mass
                    endif

                    mass             = mass + entrainment(0) + detrainment(0)
                    blob_source(i,j) = blob_source(i,j) - entrainment(0) - detrainment(0)
                    
                    ! save tracer
                    do n=1,num_prog_tracers
                       tracer(n)          = tracer(n) + entrainment(n) + detrainment(n)
                       tend_blob(n,i,j,k) = tend_blob(n,i,j,k) - entrainment(n) - detrainment(n)
                    enddo

                    do n=0,num_prog_tracers 
                       EL_diag(n)%detrainment(i,j,k) = EL_diag(n)%detrainment(i,j,k) - detrainment(n) 
                       EL_diag(n)%entrainment(i,j,k) = EL_diag(n)%entrainment(i,j,k) + entrainment(n) 
                    enddo 

                 endif !bitwise_reproduction

                 ! save other variables
                 geodepth = -Xnp1(5)
                 old_lld  =  (/lon, lat, geodepth/)
                 old_ijk  =  (/i,   j,   k       /)
                 old_h    =  h
                 vel(1)   = Xnp1(2)
                 vel(2)   = Xnp1(4)
                 vel(3)   = Xnp1(6)

                 !Avoid crazy numbers if the blob has detrained to be less than small mass
                 if(mass>small_mass) then
                    do n=1,num_prog_tracers
                       field(n) = tracer(n)/mass
                    enddo
                    rhoL  = density(field(index_salt), field(index_temp), Dens%pressure_at_depth(i,j,k))
                    rhoLr = 1.0/rhoL
                    if (vert_coordinate_class==PRESSURE_BASED) then
                       rho  = rhoL
                       rhor = rhoLr
                       !else DEPTH_BASED: rho=rho0, rhor=rho0r
                    endif

                    volume = mass*rhor
                    ! Update rhoE before the free blob check
                    ! to account for blobs that moved to shallower water
                    rhoE = model_rho(i,j,k)

                    ! Check whether a blob that is not small has separated
                    if (rhoL < rhoE) then
                       ! If so, save values to the blob, unlink it and send it to an
                       ! intermediate linked list, before turning it into a free blob.
                       
                       ! Update a counter
                       if (this%i /= i .or. this%j /= j) move_lateral = move_lateral + 1
                       
                       this%blob_time = dtime
                       this%step      = max(tstep, minstep)

                       this%i         = i+Dom%ioff
                       this%j         = j+Dom%joff
                       this%k         = k
                       this%geodepth  = geodepth
                       this%lon       = lon
                       this%lat       = lat
                       this%ent       = ent
                       this%det       = det
                       this%richardson=richardson
                       this%v(:)      = vel(:)
                       this%h1        = h(1)
                       this%h2        = h(2)
                       this%model_steps = this%model_steps + 1
                       
                       this%mass      = mass
                       this%volume    = volume
                       this%density   = rhoL
                       this%densityr  = rhoLr
                       this%gprime    = grav*(rhoE-rhoL)/rhoE !diagnostic
                       this%age       = this%age + dtime_yr
                       do n=1,num_prog_tracers
                          if(T_prog(n)%name(1:3) =='age') then
                             ! If it is an age tracer advance the age of the tracer
                             this%field(n)  = field(n) + dtime_yr
                             this%tracer(n) = this%field(n)*mass
                          else
                             this%tracer(n) = tracer(n)
                             this%field(n)  = field(n)
                          endif
                       enddo
                       this%model_steps = this%model_steps + 1
                       call unlink_blob(this, head, prev, next)
                       call insert_blob(this, free)
                       advance_blob = .false.
                       this%nfrac_steps = n_frac_steps
                       this=>next
                       tfer_free = tfer_free + 1
                       

                       exit partialstep
                    endif !rhoL<rhoE
                    
                 endif !mass>small_mass

                 this%nfrac_steps = n_frac_steps
              endif !accept_step
              m=m+1
           enddo trystep

           ! The step has been accepted.  We now decide whether we need
           ! to conduct any more steps, and if we do, whether we need to 
           ! adjust them so that the final step coincides with the Eulerian
           ! model time step.
           blob_time = blob_time + old_tstep
           this%blob_time = blob_time
           
           if (blob_time == dtime) then
              ! no more steps required.  Save all the new variables to the blob
              reached_end = .true.

              ! update a counter
              if (this%i /= i .or. this%j /= j) move_lateral = move_lateral + 1

              ! Update the blob variables
              this%step = max(tstep,minstep)
              this%i    = i
              this%j    = j
              this%k    = k
              this%lon  = lon
              this%lat  = lat
              this%v(:) = vel(:)
              this%h1   = h(1)
              this%h2   = h(2)

              this%mass     = mass
              this%volume   = volume
              this%density  = rhoL
              this%densityr = rhoLr
              this%ent      = ent
              this%det      = det
              this%richardson=richardson
              this%gprime   = grav*(rhoE-rhoL)/rhoE !diagnostic
              this%age      = this%age + dtime_yr
              do n=1,num_prog_tracers
                 if(T_prog(n)%name(1:3) =='age') then
                    ! If it is an age tracer advance the age of the tracer
                    this%field(n)  = field(n) + dtime_yr
                    this%tracer(n) = this%field(n)*mass
                 else
                    this%tracer(n) = tracer(n)
                    this%field(n)  = field(n)
                 endif
              enddo
              this%model_steps = this%model_steps + 1

              ! Adjust the geodepth so we are sitting on topography
              topog = 0.0
              call interp_tcoeff(i,j,h(:),lon,lat,tdsq_r(:))
              tbigd = 0.0
              do mm=1,Info%tidx(0,i,j)
                 iit=i+Info%it(Info%tidx(mm,i,j))
                 jjt=j+Info%jt(Info%tidx(mm,i,j))
                 topog = topog + ht(iit,jjt)*tdsq_r(mm)
                 tbigd = tbigd + tdsq_r(mm)
              enddo
              tbigd = 1.0/tbigd

              this%geodepth =  topog*tbigd
              this%depth    =  this%geodepth + Ext_mode%eta_t(i,j,tau)
              if (vert_coordinate_class==DEPTH_BASED) then
                 this%st = -Thickness%depth_swt(i,j,k)
              else !PRESSURE_BASED
                 this%st = Thickness%depth_swt(i,j,k)
              endif
        
              if (bitwise_reproduction) then
                 ! check to see if the blob has left this pe
                 ! NOTE: we assume the blob does not go beyond the halo, so, we 
                 ! don't pack them into buffers until the end of the time step
                 ! We do this to maintain bitwise reproduction.  If we are not
                 ! maintaining bitwise reproduction, we pack a blob into a buffer
                 ! as soon as a change in pe is detected.
                 !
                 ! We also consider what happens when a blob traverses three
                 ! PE's in a time step.  e.g.
                 ! ------------------
                 ! |       |        |
                 ! | PE=2  +a  PE=3 |
                 ! |      /|        |
                 ! |     b |        |
                 ! ------+-----------
                 ! |    /  |        |
                 ! |   c   |        |
                 ! |       |        |
                 ! |  PE=0 |   PE=1 |
                 ! ------------------
                 ! Maintaining bitwise reproducability requires that each
                 ! PE knows about the history of the blob.  In the instance
                 ! depicted above, the blobs history must be processed on PE
                 ! 3,2 and 0.  We achieve this by keeping the history on
                 ! 3, and sending the history to 2 and 0.  When unpacking
                 ! the buffer, the unpacking subroutine checks if the blob 
                 ! is on that PE.  If it is not, it will set the blob%mass 
                 ! and blob%tracer(n) to zero.

                 if (change_pe) then
                    go_ne = .false.
                    go_nw = .false.
                    go_sw = .false.
                    go_se = .false.
                    go_e  = .false.
                    go_n  = .false.
                    go_w  = .false.
                    go_s  = .false.
                    do s=0,this%nfrac_steps

                       ! If we are on a cyclic grid, and, there is only one
                       ! processor in the cyclic direction, we want to avoid
                       ! having the blob sent and received to the same list.
                       ! (which means it can appear twice in the list). 
                       ! This can cause non-conservation.
                       call check_cyclic(this, this%di(s), this%dj(s), .false.)

                       if     (iec<this%di(s) .and. jec<this%dj(s)) then
                          go_ne = .true.
                       elseif (this%di(s)<isc .and. jec<this%dj(s)) then
                          go_nw = .true.
                       elseif (this%di(s)<isc .and. this%dj(s)<jsc) then
                          go_sw = .true.
                       elseif (iec<this%di(s) .and. this%dj(s)<jsc) then
                          go_se = .true.
                       elseif (iec<this%di(s)) then
                          go_e  = .true.
                       elseif (jec<this%dj(s)) then
                          go_n = .true.
                       elseif (this%di(s)<isc) then
                          go_w  = .true.
                       elseif (this%dj(s)<jsc) then
                          go_s = .true.
                       endif
                    enddo
                    if (go_ne) call packbuffer(this, NEbuffer_out)
                    if (go_nw) call packbuffer(this, NWbuffer_out)
                    if (go_se) call packbuffer(this, SEbuffer_out)
                    if (go_sw) call packbuffer(this, SWbuffer_out)
                    if (go_e ) call packbuffer(this,  Ebuffer_out)
                    if (go_w ) call packbuffer(this,  Wbuffer_out)
                    if (go_n ) call packbuffer(this,  Nbuffer_out)
                    if (go_s ) call packbuffer(this,  Sbuffer_out)

                    if (.not. any((/go_ne, go_nw, go_se, go_sw, go_e, go_w, go_n, go_s/))) then
                       call check_cyclic(this, this%i, this%j, .true.)
                       i=this%i; j=this%j
                    endif

                    ! We leave the blobs history in the linked list, but, we 
                    ! remove all its properties.  It will be deleted next call
                    ! to blob_delete, but, we need its history to be processed 
                    ! before it is deleted.  The history is processed only after
                    ! all blobs have been stepped.
                    ! A blob may go from one PE and then back to the same PE within an
                    ! E system time step, e.g.
                    ! +----------------+
                    ! | PE=2  +a  PE=3 |
                    ! |      /|        |
                    ! |     b |        |
                    ! |      \|        |
                    ! |       +c       |
                    ! +-------+--------+
                    ! |       |        |
                    ! |       |        |
                    ! |       |        |
                    ! |       |        |
                    ! |  PE=0 |   PE=1 |
                    ! +----------------+
                    ! So, we only set its properties to zero if it is no longer on this PE.
                    if (i<isc .or. iec<i .or. j<jsc .or. jec<j) then
                       this%mass=0.
                       do n=1,num_prog_tracers
                          this%tracer(n)=0.
                       enddo
                    endif
                    leaving_pe = leaving_pe + 1
                 endif
              endif !bitwise_reproduction
              
           elseif( (blob_time + 1.1*tstep) > dtime) then
              ! If we are within 10% of dtime with the suggested tstep, then, 
              ! we extend the step slightly to save having to take a really 
              ! small step next time. We ensure that the step size is chosen 
              ! so that blob_time coincides with the model step after the next blob step.
              reached_end = .false.
              tstep = dtime - blob_time

              ! Update the interpolation coefficients to reflect the new position
              call interp_tcoeff(i,j,h(:),lon,lat,tdsq_r(:))
              call interp_ucoeff(i,j,h(:),lon,lat,udsq_r(:))        
           else
              ! More step(s) required.  No need to adjust the step size.
              reached_end = .false.

              ! Update the interpolation coefficients to reflect the new position
              call interp_tcoeff(i,j,h(:),lon,lat,tdsq_r(:))
              call interp_ucoeff(i,j,h(:),lon,lat,udsq_r(:))        
           endif !blob_time == dtime

        enddo partialstep

        total_blobs = total_blobs + 1

        if (advance_blob) then
           this=>this%next 
        else
           this=>next
        endif
        if (.not. associated(this)) exit blobcycle
        
     enddo blobcycle
  endif !associated(head)

  ngrnd  = ngrnd  + ngrounded
  ndetrn = ndetrn + detrn_zero

  if (debug_this_module) then
     stdoutunit = stdout()
     write (stdoutunit, '(/,a)') 'Dynamic Bottom Blob Statistics'
     call write_timestamp(Time%model_time)
     call mpp_sum(total_blobs)
     write (stdoutunit, *) 'Total Bottom Dynamic Blobs             =', total_blobs
     call mpp_sum(leaving_pe)
     write (stdoutunit, *) 'Bottom Dynamic Blobs Changing PEs      =', leaving_pe
     call mpp_sum(tfer_free)
     write (stdoutunit, *) 'Bottom Dynamic Blob Transfer to free   =', tfer_free
     call mpp_sum(detrn_zero)
     write (stdoutunit, *) 'Bottom Blobs detrained to zero mass    =', detrn_zero
     call mpp_sum(move_lateral)
     write (stdoutunit, *) 'Bottom Blobs changing water columns    =', move_lateral
     call mpp_sum(ngrounded)
     write (stdoutunit, *) 'Grounded Bottom Blobs                  =', ngrounded
  endif
end subroutine dynamic_update
! </SUBROUTINE>  NAME="dynamic_update"


!######################################################################
! <SUBROUTINE NAME="transfer_free_to_bottom">
!
! <DESCRIPTION>
! Takes free blobs that have interacted with topography and turns them
! into bottom blobs.
! </DESCRIPTION>
!
subroutine transfer_free_to_bottom(Time, Thickness, T_prog, new_head, &
                                   head, use_free, EL_diag)
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_blob_type),        pointer       :: new_head
  type(ocean_blob_type),        pointer       :: head
  logical,                      intent(in)    :: use_free
  type(blob_diag_type), dimension(0:num_prog_tracers) :: EL_diag

  type(ocean_blob_type), pointer :: this => NULL()
  type(ocean_blob_type), pointer :: prev => NULL()
  type(ocean_blob_type), pointer :: next => NULL()
  real,    dimension(3) :: v_par, v_perp, v_new
  real, dimension(9) :: tdsq_r
  real, dimension(4) :: udsq_r
  real :: ubigd, tbigd
  integer :: i,j,k,iiu,jju,iit,jjt,mm,n
  integer :: nblobs
  real    :: Hx, Hy
  real    :: topog
  real    :: mag_v, mag_v_par, constant
  integer :: stdoutunit

  nblobs = 0
  if (associated(new_head)) then
     this=>new_head
     if (use_this_module .and. accept_free_blobs) then
        ! If we are using dynamic bottom blobs, transfer free blobs
        ! that have interacted with topography to the dynamic 
        ! bottom blob list
        blobcycle: do
           i = this%i
           j = this%j
           k = this%k

           call unlink_blob(this, new_head, prev, next)
           call insert_blob(this, head)
           call reallocate_interaction_memory(this, head, total_ns)

           call interp_ucoeff(i,j,(/this%h1,this%h2/),this%lon,this%lat,udsq_r(:))
           ubigd = 0.0
           Hx=0.0; Hy=0.0
           do mm=1,Info%uidx(0,i,j)
              iiu=i+Info%iu(Info%uidx(mm,i,j))
              jju=j+Info%ju(Info%uidx(mm,i,j))
              Hx = Hx + dht_dx(iiu,jju)*udsq_r(mm)
              Hy = Hy + dht_dy(iiu,jju)*udsq_r(mm)
              ubigd = ubigd + udsq_r(mm)
           enddo
           ubigd = 1.0/ubigd
           Hx = Hx*ubigd
           Hy = Hy*ubigd

           constant  = this%v(1)*Hx + this%v(2)*Hy + this%v(3)
           constant  = constant/(Hx**2 + Hy**2 + 1)
           v_perp(:) = (/Hx*constant, Hy*constant, constant/)
           v_par(:)  = this%v(:) - v_perp(:)
           mag_v     = sqrt(this%v(1)**2 + this%v(2)**2 + this%v(3)**2)
           mag_v_par = sqrt( v_par(1)**2 + v_par(2)**2 + v_par(3)**2 )
           if (mag_v_par /= 0.0) then
              v_new(:) = elastic*v_par(:)*mag_v/mag_v_par
           else
              v_new(:) = 0.0
           endif

           this%v(:) = v_new(:)

           this%height = blob_height

           call interp_tcoeff(i,j,(/this%h1,this%h2/),this%lon,this%lat,tdsq_r(:))
           
           tbigd = 0.0
           topog = 0.0
           do mm=1,Info%tidx(0,i,j)
              iit=i+Info%it(Info%tidx(mm,i,j))
              jjt=j+Info%jt(Info%tidx(mm,i,j))
              topog = topog + ht(iit,jjt)*tdsq_r(mm)
              tbigd = tbigd+tdsq_r(mm)
           enddo
           tbigd = 1.0/tbigd

           this%geodepth = topog*tbigd
           if (vert_coordinate_class==DEPTH_BASED) then
              this%st = -Thickness%depth_swt(i,j,k)
           else !PRESSURE_BASED
              this%st = Thickness%depth_swt(i,j,k)
           endif
           
           this=>next
           nblobs = nblobs+1
           if (.not. associated(this)) exit blobcycle
        enddo blobcycle

     else !not use_this_module
        ! If we are not using dynamic bottom blobs, return their
        ! properties to the E system.
        blobcycle2: do
           i = this%i
           j = this%j
           k = this%k

           ! Diagnostics
           EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + this%mass
           do n=1,num_prog_tracers
              EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + this%tracer(n)
           enddo

           call kill_blob(Thickness, T_prog(:), this, i, j, k)
           this=>this%next
           nblobs = nblobs+1
           if (.not. associated(this)) exit blobcycle2
        enddo blobcycle2
     endif !use_this_module
     call blob_delete(Time, Thickness, T_prog(:), new_head)
  endif!associated(new_head)

  call mpp_sum(nblobs)
  if (id_free_to_bot>0) used = send_data(id_free_to_bot, real(nblobs), Time%model_time)

  if(debug_this_module .and. use_free) then
     stdoutunit = stdout()
     write(stdoutunit, '(/,a)') 'Free blobs interacting with topography'
     if (use_this_module) then
        write(stdoutunit, *) 'Free blobs transferred to bottom blobs = ', nblobs
     else
        write(stdoutunit, *) 'Free blobs destroyed from reaching bottom = ', nblobs
     endif
  endif

end subroutine transfer_free_to_bottom
! </SUBROUTINE>  NAME="transfer_free_to_bottom"


!######################################################################
! <SUBROUTINE NAME="dynamic_bottom_form_new">
!
! <DESCRIPTION>
! Initialises blobs that are formed when an on-shelf/off-shelf
! instability occurs.  The method used for determining an instability
! and the initial conditions are based on that of Campin and Goosse
! (1999).  
!
! When the density difference between the shallow ocean cell and the 
! deep ocean cell (referenced to the deep ocean cell) exceeds the
! namelist variable rho_threshold, a blob is formed.  The deep ocean
! cell is chosen based on which deep ocean cell (in the k plane)
! the blob topography intersects.
!
! After formation, the new blobs are added to the bottom blob linked
! list, and, their properties are integrated, starting at time taup1.
! </DESCRIPTION>
!
subroutine dynamic_bottom_form_new(Time, Dens, T_prog, Thickness, Ext_mode, &
                                   L_system, tend_blob, blob_counter, u,    &
                                   mass_out, mass_in, head, EL_diag)
  type(ocean_time_type),                  intent(in)    :: Time
  type(ocean_density_type),               intent(in)    :: Dens
  type(ocean_prog_tracer_type),           intent(inout) :: T_prog(:)
  type(ocean_thickness_type),             intent(inout) :: Thickness
  type(ocean_external_mode_type),         intent(inout) :: Ext_mode
  type(ocean_lagrangian_type),            intent(inout) :: L_system
  real, dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  integer, dimension(isc:iec,jsc:jec,nk),               intent(inout) :: blob_counter
  real, dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),        intent(in)    :: u
  real, dimension(isd:ied,jsd:jed,1:nk),                intent(inout) :: mass_out
  real, dimension(isd:ied,jsd:jed,1:nk),                intent(inout) :: mass_in
  type(ocean_blob_type), pointer                                      :: head
  type(blob_diag_type),  dimension(0:num_prog_tracers), intent(inout) :: EL_diag
  
  type(ocean_blob_type), pointer :: this, prev, next, new_head
  real, dimension(4) :: udsq_r
  real :: ubigd
  integer :: tau, i, iip, iiu, j, jjp, jju, k, m, mm, n, nblobs
  integer :: stdoutunit
  real :: overflow_flux, overflow_speed
  real :: blob_mass, do_dens, so_dens, drho, drho_rhor
  real :: uE, vE

  if (.not. use_this_module) return

  nblobs = 0

  tau   = Time%tau

  ! For diagnostics
  do n=1,num_prog_tracers
     T_prog(n)%wrk1(:,:,:) = 0.0
  enddo

  call clear_buffer(Ebuffer_out)
  call clear_buffer(Wbuffer_out)
  call clear_buffer(Nbuffer_out)
  call clear_buffer(Ebuffer_out)

  nullify(new_head)

  if (debug_this_module) then
     stdoutunit = stdout()
     write(stdoutunit, '(a)') ' '
     call write_timestamp(Time%model_time)
     write(stdoutunit, '(a)') 'From ocean_blob_dynamics_bottom_mod'
     write(stdoutunit, '(a)') 'Totals before bottom dynamic blob creation (tau)'
     call E_and_L_totals(L_system,Thickness,T_prog(:),tau)
     call mpp_sum(nblobs)
     write(stdoutunit, '(a,i4,/)') 'New bottom dynamic blobs (tau) = ', nblobs
  endif

  ! search for horizontal instability, including in certain parts of the halo
  do m=1,4
     do j=jsc,jec
        do i=isc,iec

           ! check whether downslope flow is possible and it is within
           ! the latitude bands of interest
           if (topog_step(i,j,m) == 1.0) then

              ! some convenient variables
              iip = i+ip(m)
              jjp = j+jp(m)
              k   = Grd%kmt(i,j)

              ! Check whether density of E system shelf water is more than deep water.
              so_dens = density(T_prog(index_salt)%field(i,j,k,tau), &
                                T_prog(index_temp)%field(i,j,k,tau), &
                                Dens%pressure_at_depth(iip,jjp,kmt(iip,jjp)))
              do_dens = Dens%rho(iip,jjp,kmt(iip,jjp),tau)
              if ((so_dens - do_dens) > rho_threshold) then

                 nullify(this)
                 
                 drho      = abs(so_dens-do_dens)
                 drho_rhor = drho/so_dens
                 
                 ! calculate the overflow velocity
                 overflow_speed = overflow_factor*topog_slope(i,j,m)*drho_rhor
                 
                 ! Overflow flux has units kg/(m*s) 
                 overflow_flux = overflow_factor * topog_slope(i,j,m)      & 
                                 *drho * Thickness%dzt(i,j,k) 

                 if (enforce_big_blobs) then
                    ! 0.2 = 0.8*1/4
                    blob_mass = 0.2*Grd%dat(i,j)*Thickness%rho_dzt(i,j,k,tau)

                 else
                    if    (m==1) then 
                       blob_mass = overflow_flux*Grd%dyte(i,j)*dtime 
                    elseif(m==2) then 
                       blob_mass = overflow_flux*Grd%dxtn(i,j)*dtime 
                    elseif(m==3) then 
                       blob_mass = overflow_flux*Grd%dyte(i-1,j)*dtime 
                    else  !m==4 
                       blob_mass = overflow_flux*Grd%dxtn(i,j-1)*dtime 
                    endif

                    ! There can be, theoretically, up to four overflow events take 
                    ! place from a single grid cell.  Therefore, we do not want a single 
                    ! overflow event to take away more than 1/4 of the mass of a grid cell. 
                    if (0.25*Thickness%rho_dzt(i,j,k,tau)*Grd%dat(i,j) < blob_mass) then 
                       blob_mass = 0.25*Thickness%rho_dzt(i,j,k,tau)*Grd%dat(i,j) 
                    endif
                 endif


                 ! Form the onshelf blobs.  Make sure they are in the 
                 ! compute domain.
                 allocate(this)
                 allocate(this%tracer(num_prog_tracers))
                 allocate(this%field(num_prog_tracers))
                 
                 this%mass = blob_mass
                 do n=1,num_prog_tracers
                    this%tracer(n)     = this%mass*T_prog(n)%field(i,j,k,tau)
                    this%field(n)      = this%tracer(n)/this%mass
                    tend_blob(n,i,j,k) = tend_blob(n,i,j,k) - this%tracer(n)
                 enddo
                 this%density = density(this%field(index_salt), &
                                        this%field(index_temp), &
                                        Dens%pressure_at_depth(iip,jjp,k))
                 this%densityr = 1.0/this%density

                 ! Calculate the divergence (convergence comes later)
                 Ext_mode%conv_blob(i,j)   = Ext_mode%conv_blob(i,j)   - this%mass
                 L_system%conv_blob(i,j,k) = L_system%conv_blob(i,j,k) - this%mass
                 mass_out(i,j,k)           = mass_out(i,j,k)           + this%mass
                 
                 this%height = blob_height

                 if (vert_coordinate_class == DEPTH_BASED) then
                    this%volume = this%mass*rho0r
                 else
                    this%volume = this%mass*this%densityr
                 endif
                 
                 ! Firstly, give the blob the grid cell of origin, so that it is
                 ! counted correctly.
                 this%i = i
                 this%j = j
                 this%k = k

                 ! We need to keep the blob in a separate list for the moment
                 ! so that we can maintain reproducibility
                 call count_blob(this, blob_counter)
                 call insert_blob(this, new_head)
                 call allocate_interaction_memory(this, total_ns)
                 nblobs = nblobs + 1

                 this%age = 0.0

                 ! Diagnostic
                 EL_diag(0)%new(i,j,k) = EL_diag(0)%new(i,j,k) + this%mass
                 do n=1,num_prog_tracers
                    EL_diag(n)%new(i,j,k) = EL_diag(n)%new(i,j,k) + this%tracer(n)
                 enddo

                 this%i = iip
                 this%j = jjp
                 
                 ! Find the initial position and stretching function
                 ! We will adjust the initial position after the initial velocity
                 ! has been calculated
                 this%h1  = blobh1(i,j,m)
                 this%h2  = blobh2(i,j,m)
                 this%lon = xb(i,j,m)
                 this%lat = yb(i,j,m)

                 this%blob_time = 0.0
                 
                 this%geodepth = hb(i,j,m)
                 this%depth    = this%geodepth + Ext_mode%eta_t(i,j,tau)
                 if (vert_coordinate_class==DEPTH_BASED) then
                    this%st = -Thickness%depth_swt(i,j,k)
                 else !PRESSURE_BASED
                    this%st = Thickness%depth_swt(i,j,k)
                 endif

                 call interp_ucoeff(i,j,(/this%h1,this%h2/),this%lon,this%lat,udsq_r(:))        
                 ! Interpolate the E system variables to the blobs
                 ! Note, we treat the stretching function separately
                 ubigd = 0.0
                 uE=0.0; vE=0.0
                 do mm=1,Info%uidx(0,i,j)
                    iiu=i+Info%iu(Info%uidx(mm,i,j))
                    jju=j+Info%ju(Info%uidx(mm,i,j))
                    ! Land points are treated as a zero
                    if(kmu(iiu,jju)>0) then
                       uE = uE + u(iiu,jju,kmu(iiu,jju),1)*udsq_r(mm)
                       vE = vE + u(iiu,jju,kmu(iiu,jju),2)*udsq_r(mm)
                    endif
                    ubigd = ubigd + udsq_r(mm)
                 enddo
                 ubigd = 1.0/ubigd
                 
                 uE = uE*ubigd
                 vE = vE*ubigd

                 ! Initial Velocity
                 if(m==1) then
                    this%v(1) = uE + overflow_speed
                    this%v(2) = vE
                 elseif(m==2) then
                    this%v(1) = uE 
                    this%v(2) = vE + overflow_speed
                 elseif (m==3) then
                    this%v(1) = uE - overflow_speed
                    this%v(2) = vE
                 else!m==4
                    this%v(1) = uE 
                    this%v(2) = vE - overflow_speed
                 endif
                 this%v(3) = vert_factor(i,j,m)*drho_rhor

                 ! Adjust the initial position (this guarantees we are in
                 ! the correct cell, and not exactly on the edge of two cells)
                 this%lon = this%lon + this%v(1)*dtime_rad_to_deg/this%h1
                 this%lat = this%lat + this%v(2)*dtime_rad_to_deg/this%h2

                 this%new  = .false.
                 this%step = first_step
                 
                 ! test if this has gone out to the E,W,N or S of the 
                 ! compute domain.  If it has, then pack it into an outward
                 ! buffer and remove it from this PE's list of blobs
                 if (iip>iec) then
                    call packbuffer(this, Ebuffer_out)
                    call unlink_blob(this, new_head, prev, next)
                    call free_blob_memory(this)
                 elseif (iip<isc) then
                    call packbuffer(this, Wbuffer_out)
                    call unlink_blob(this, new_head, prev, next)
                    call free_blob_memory(this)
                 elseif (jjp>jec) then
                    call packbuffer(this, Nbuffer_out)
                    call unlink_blob(this, new_head, prev, next)
                    call free_blob_memory(this)
                 elseif (jjp<jsc) then
                    call packbuffer(this, Sbuffer_out)
                    call unlink_blob(this, new_head, prev, next)
                    call free_blob_memory(this)
                 endif
              endif !rho_so>rho_do?
           endif !topog_step==1.0?
        enddo !i
     enddo !j
  enddo !m
  
  ! Now, send blobs to neighbouring PE's.  Note: since blobs can only move
  ! N,S,E or W, we do not need to worry about diagonal PE's 
  ! (i.e. NE,SE,SW,NW)
  call send_buffer(Ebuffer_out)
  call send_buffer(Wbuffer_out)
  call send_buffer(Nbuffer_out)
  call send_buffer(Sbuffer_out)
  
  call mpp_sync_self()

  ! Clear the incoming buffers
  call clear_buffer(Wbuffer_in)
  call clear_buffer(Ebuffer_in)
  call clear_buffer(Sbuffer_in)
  call clear_buffer(Nbuffer_in)
  
  ! Now receive blobs from neighbouring PE's
  call receive_buffer(Wbuffer_in)
  call receive_buffer(Ebuffer_in)
  call receive_buffer(Sbuffer_in)
  call receive_buffer(Nbuffer_in)
  
  ! Now unpack the buffer and add the received blobs to this PE's blob list
  call unpackbuffer(Time, Wbuffer_in, new_head, Dens, Ext_mode, Thickness)
  call unpackbuffer(Time, Ebuffer_in, new_head, Dens, Ext_mode, Thickness)
  call unpackbuffer(Time, Sbuffer_in, new_head, Dens, Ext_mode, Thickness)
  call unpackbuffer(Time, Nbuffer_in, new_head, Dens, Ext_mode, Thickness)

  if(associated(new_head)) then
     this=>new_head
     blobcycle: do
        i = this%i
        j = this%j
        k = this%k

        ! Calculate the convergence
        Ext_mode%conv_blob(i,j)   = Ext_mode%conv_blob(i,j)   + this%mass
        L_system%conv_blob(i,j,k) = L_system%conv_blob(i,j,k) + this%mass
        mass_in(i,j,k)            = mass_in(i,j,k)            + this%mass

        call unlink_blob(this, new_head, prev, next)
        call insert_blob(this, head)
        
        this=>next
        if(.not.associated(this)) exit blobcycle
     enddo blobcycle
  endif
  
  call mpp_sum(nblobs)
  if (id_new_blobs>0) used = send_data(id_new_blobs, real(nblobs), Time%model_time)

  if (debug_this_module) then
     call write_timestamp(Time%model_time)
     write(stdoutunit, '(a)') 'From ocean_blob_dynamics_bottom_mod'
     write(stdoutunit, '(a)') 'Totals after bottom dynamic blob creation (tau)'
     call E_and_L_totals(L_system,Thickness,T_prog(:),tau)
     write(stdoutunit, '(a,i4)') 'New bottom dynamic blobs (taup1) = ', nblobs
     write(stdoutunit, '(a)') ' '
  endif

end subroutine dynamic_bottom_form_new
! </SUBROUTINE>  NAME="dynamic_bottom_form_new"


!######################################################################
! <SUBROUTINE NAME="blob_dynamic_bottom_end">
!
! <DESCRIPTION>
! Clears memory to give a nice clean ending to the run.
! </DESCRIPTION>
!
subroutine blob_dynamic_bottom_end()

  if (.not. use_this_module) return

  call deallocate_buffer(Ebuffer_out);  call deallocate_buffer(Ebuffer_in)
  call deallocate_buffer(Wbuffer_out);  call deallocate_buffer(Wbuffer_in)
  call deallocate_buffer(Nbuffer_out);  call deallocate_buffer(Nbuffer_in)
  call deallocate_buffer(Sbuffer_out);  call deallocate_buffer(Sbuffer_in)
  call deallocate_buffer(NEbuffer_out); call deallocate_buffer(NEbuffer_in)
  call deallocate_buffer(NWbuffer_out); call deallocate_buffer(NWbuffer_in)
  call deallocate_buffer(SEbuffer_out); call deallocate_buffer(SEbuffer_in)
  call deallocate_buffer(SWbuffer_out); call deallocate_buffer(SWbuffer_in)

  nullify(Dom)
  nullify(Grd)
  nullify(Info)
  nullify(Bdom)

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that deallocates memory from a buffer    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine deallocate_buffer(buffer)
    type(blob_buffer_type), pointer :: buffer
    if (buffer%pe /= NULL_PE) deallocate(buffer)
  end subroutine deallocate_buffer
end subroutine blob_dynamic_bottom_end
! </SUBROUTINE>  NAME="blob_dynamic_bottom_end"


!######################################################################
! <SUBROUTINE NAME="packbuffer">
!
! <DESCRIPTION>
! Packs a buffer with all the information needed to send a blob from
! one PE to another.
! </DESCRIPTION>
!
subroutine packbuffer(blob, buffer, entrainment, detrainment)
  type(ocean_blob_type),  pointer    :: blob
  type(blob_buffer_type), pointer    :: buffer
  real,         optional, intent(in) :: entrainment(0:)
  real,         optional, intent(in) :: detrainment(0:)
  integer :: n, nb, s, npt
  integer :: stdoutunit

  stdoutunit = stdout()

  if (buffer%pe == NULL_PE) then
     write (stdoutunit, '(a)') 'Error: Trying to send blob to a NULL_PE'
     call mpp_error(FATAL, &
          '==>Error in ocean_blob_static_bottom_mod (packbuffer): '&
          //'Trying to send blob to a NULL_PE')
  endif
  buffer%numblobs = buffer%numblobs+1
  if (buffer%numblobs>buffer%size) call increase_buffer(buffer,buffer%numblobs)

  npt = num_prog_tracers
  nb  = buffer%numblobs
  
  ! Fill the real buffer
  do n=1,num_prog_tracers
     buffer%real_buff(2*n-1,nb) = blob%tracer(n)
     buffer%real_buff(2*n  ,nb) = blob%dtracer(n)
  enddo
  buffer%real_buff(2*npt+1 ,nb) = blob%v(1)
  buffer%real_buff(2*npt+2 ,nb) = blob%v(2)
  buffer%real_buff(2*npt+3 ,nb) = blob%v(3)
  buffer%real_buff(2*npt+4 ,nb) = blob%lon
  buffer%real_buff(2*npt+5 ,nb) = blob%lat
  buffer%real_buff(2*npt+6 ,nb) = blob%geodepth
  buffer%real_buff(2*npt+7 ,nb) = blob%step
  buffer%real_buff(2*npt+8 ,nb) = blob%mass
  buffer%real_buff(2*npt+9 ,nb) = blob%blob_time
  buffer%real_buff(2*npt+10,nb) = blob%h1
  buffer%real_buff(2*npt+11,nb) = blob%h2
  buffer%real_buff(2*npt+12,nb) = blob%ent
  buffer%real_buff(2*npt+13,nb) = blob%det
  buffer%real_buff(2*npt+14,nb) = blob%richardson
  buffer%real_buff(2*npt+15,nb) = blob%dmass
  buffer%real_buff(2*npt+16,nb) = blob%gprime
  buffer%real_buff(2*npt+17,nb) = blob%age
  buffer%real_buff(2*npt+18,nb) = blob%height

  ! Fill the integer buffer
  buffer%integer_buff(1,nb) = blob%i
  buffer%integer_buff(2,nb) = blob%j
  buffer%integer_buff(3,nb) = blob%k
  buffer%integer_buff(4,nb) = blob%model_steps
  buffer%integer_buff(5,nb) = blob%hash
  buffer%integer_buff(6,nb) = blob%number
  buffer%integer_buff(7,nb) = blob%nfrac_steps
  buffer%integer_buff(8,nb) = blob%nsteps

  if (bitwise_reproduction) then
     ! History buffers
     
     ! only some variables use s=0
     s=0
     buffer%history_integer_buff(1,s,nb) = blob%di(s)
     buffer%history_integer_buff(2,s,nb) = blob%dj(s)
     buffer%history_integer_buff(3,s,nb) = blob%dk(s)
     buffer%history_real_buff(npt+3,s,nb) = blob%mass_out(s)
     
     do s=1,blob%nfrac_steps
        ! Integer buffer
        buffer%history_integer_buff(1,s,nb) = blob%di(s)
        buffer%history_integer_buff(2,s,nb) = blob%dj(s)
        buffer%history_integer_buff(3,s,nb) = blob%dk(s)
        
        ! Real buffer
        do n=0,num_prog_tracers
           buffer%history_real_buff(2*n  ,s,nb) = blob%entrainment(s,n)
           buffer%history_real_buff(2*n+1,s,nb) = blob%detrainment(s,n)
        enddo
        buffer%history_real_buff(2*npt+2,s,nb) = blob%mass_in(s)
        buffer%history_real_buff(2*npt+3,s,nb) = blob%mass_out(s)
     enddo

  else
     ! dtracer is not parsed for new blobs, so we do
     ! not need to worry about it
     if (present(entrainment)) then
        do n=0,num_prog_tracers
           buffer%history_real_buff(2*n  ,1,nb) = entrainment(n)
           buffer%history_real_buff(2*n+1,1,nb) = detrainment(n)
        enddo
     endif
  endif

end subroutine packbuffer
! </SUBROUTINE>  NAME="packbuffer"

!######################################################################
! <SUBROUTINE NAME="unpackbuffer">
!
! <DESCRIPTION>
! Unpacks a received buffer.
! </DESCRIPTION>
!
subroutine unpackbuffer(Time, buffer, head, Dens, Ext_mode, Thickness, &
                        L_system, EL_diag, blob_source, tend_blob, mass_in)
  type(ocean_time_type),                    intent(in)    :: Time
  type(blob_buffer_type),                   pointer       :: buffer
  type(ocean_blob_type),                    pointer       :: head
  type(ocean_density_type),                 intent(in)    :: Dens
  type(ocean_external_mode_type),           intent(inout) :: Ext_mode
  type(ocean_thickness_type),               intent(in)    :: Thickness
  type(ocean_lagrangian_type),    optional, intent(inout) :: L_system
  type(blob_diag_type),           optional, intent(inout) :: EL_diag(0:)
  real, optional, dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  real, optional, dimension(isd:ied,jsd:jed),                     intent(inout) :: blob_source
  real, optional, dimension(isd:ied,jsd:jed,1:nk),                intent(inout) :: mass_in

  type(ocean_blob_type), pointer       :: blob
  real, dimension(0:num_prog_tracers) :: entrainment, detrainment
  integer :: n, nb, s, npt
  integer :: i,j,k

  if (buffer%pe /= NULL_PE .and. buffer%numblobs>0) then
     npt = num_prog_tracers
     
     do nb=1,buffer%numblobs
        allocate(blob)
        allocate(blob%tracer(num_prog_tracers))
        allocate(blob%field(num_prog_tracers))
        call allocate_interaction_memory(blob, total_ns)
        
        ! unpack the real buffer
        do n=1,num_prog_tracers
           blob%tracer(n)  = buffer%real_buff(2*n-1,nb)
           blob%dtracer(n) = buffer%real_buff(2*n  ,nb)
        enddo
        blob%v(1)       = buffer%real_buff(2*npt+ 1,nb)
        blob%v(2)       = buffer%real_buff(2*npt+ 2,nb)
        blob%v(3)       = buffer%real_buff(2*npt+ 3,nb)
        blob%lon        = buffer%real_buff(2*npt+ 4,nb)
        blob%lat        = buffer%real_buff(2*npt+ 5,nb)
        blob%geodepth   = buffer%real_buff(2*npt+ 6,nb)
        blob%step       = buffer%real_buff(2*npt+ 7,nb)
        blob%mass       = buffer%real_buff(2*npt+ 8,nb)
        blob%blob_time  = buffer%real_buff(2*npt+ 9,nb)
        blob%h1         = buffer%real_buff(2*npt+10,nb)
        blob%h2         = buffer%real_buff(2*npt+11,nb)
        blob%ent        = buffer%real_buff(2*npt+12,nb)
        blob%det        = buffer%real_buff(2*npt+13,nb)
        blob%richardson = buffer%real_buff(2*npt+14,nb)
        blob%dmass      = buffer%real_buff(2*npt+15,nb)
        blob%gprime     = buffer%real_buff(2*npt+16,nb)
        blob%age        = buffer%real_buff(2*npt+17,nb)
        blob%height     = buffer%real_buff(2*npt+18,nb)

        ! unpack the integer buffer
        blob%i           = buffer%integer_buff(1,nb)
        blob%j           = buffer%integer_buff(2,nb)
        blob%k           = buffer%integer_buff(3,nb)
        blob%model_steps = buffer%integer_buff(4,nb)
        blob%hash        = buffer%integer_buff(5,nb)
        blob%number      = buffer%integer_buff(6,nb)
        blob%nfrac_steps = buffer%integer_buff(7,nb)
        blob%nsteps      = buffer%integer_buff(8,nb)
        call check_cyclic(blob, blob%i, blob%j, .true.)
        
        i=blob%i; j=blob%j; k=blob%k

        ! We check to make sure that the blob belongs on this PE, or,
        ! whether we only need to process its history.  If it does not
        ! belong on this PE, set its mass and tracer to zero so that 
        ! it will be destroyed.  If it does not belong, then we do not 
        ! need to bother with finding tracer concentration, density,
        ! or volume.

        if (isc <= i .and. i<=iec .and. jsc<=j .and. j<=jec .and. blob%mass>0.0) then
        ! Derived variables
           do n=1,num_prog_tracers
              blob%field(n)  = blob%tracer(n)/blob%mass
           enddo

           blob%density = density(blob%field(index_salt), & 
                                  blob%field(index_temp), &
                                  Dens%pressure_at_depth(i,j,k))
           blob%densityr = 1./blob%density
        
           if (vert_coordinate_class == DEPTH_BASED) then
              blob%volume = blob%mass * rho0r 
           else
              blob%volume = blob%mass * blob%densityr
           endif

           if (vert_coordinate_class==DEPTH_BASED) then
              blob%st = -Thickness%depth_swt(i,j,k)
           else !PRESSURE_BASED
              blob%st = Thickness%depth_swt(i,j,k)
           endif
           blob%depth =  blob%geodepth + Ext_mode%eta_t(i,j,Time%tau)
           
        else
           blob%mass      = 0.0
           blob%tracer(:) = 0.0
        endif

        if (bitwise_reproduction) then
           ! History buffers
           blob%di(:) = -1
           blob%dj(:) = -1
           blob%dk(:) = -1
           blob%entrainment(:,:) = 0.0
           blob%detrainment(:,:) = 0.0
           blob%mass_in(:)       = 0.0
           blob%mass_out(:)      = 0.0
           
           ! only some variables use s=0
           s=0
           blob%di(s) = buffer%history_integer_buff(1,s,nb)
           blob%dj(s) = buffer%history_integer_buff(2,s,nb)
           blob%dk(s) = buffer%history_integer_buff(3,s,nb)
           blob%mass_out(s) = buffer%history_real_buff(npt+3,s,nb)
           
           if (blob%nfrac_steps>0) then
              do s=1,blob%nfrac_steps
                 ! Integer buffer
                 blob%di(s) = buffer%history_integer_buff(1,s,nb)
                 blob%dj(s) = buffer%history_integer_buff(2,s,nb)
                 blob%dk(s) = buffer%history_integer_buff(3,s,nb)
                 ! Real buffer
                 do n=0,num_prog_tracers
                    blob%entrainment(s,n) = buffer%history_real_buff(2*n  ,s,nb)
                    blob%detrainment(s,n) = buffer%history_real_buff(2*n+1,s,nb)
                 enddo
                 blob%mass_in(s)  = buffer%history_real_buff(2*npt+2,s,nb)
                 blob%mass_out(s) = buffer%history_real_buff(2*npt+3,s,nb)
              enddo
           endif
              
           do s=0,blob%nfrac_steps
              call check_cyclic(blob, blob%di(s), blob%dj(s), .false.)
           enddo

           call insert_blob(blob, head)

        else!not bitwise_reproduction
           ! blob_source is not parsed for new blobs, so, we do 
           ! not need to worry about it
           if(present(blob_source)) then !if blob_source is present, so is tend_blob and L_system
              Ext_mode%conv_blob(i,j)   = Ext_mode%conv_blob(i,j)   + blob%mass
              L_system%conv_blob(i,j,k) = L_system%conv_blob(i,j,k) + blob%mass
              mass_in(i,j,k)            = mass_in(i,j,k)            + blob%mass

              do n=0,num_prog_tracers 
                 entrainment(n) = buffer%history_real_buff(2*n  ,1,nb) 
                 detrainment(n) = buffer%history_real_buff(2*n+1,1,nb) 
              enddo 
              blob_source(i,j) = blob_source(i,j) - detrainment(0) - entrainment(0) 
              blob%mass        = blob%mass + detrainment(0) + entrainment(0) 
              EL_diag(0)%detrainment(i,j,k) = EL_diag(0)%detrainment(i,j,k) - detrainment(0) 
              EL_diag(0)%entrainment(i,j,k) = EL_diag(0)%entrainment(i,j,k) + entrainment(0) 

              do n=1,num_prog_tracers
                 tend_blob(n,i,j,k) = tend_blob(n,i,j,k) - detrainment(n) - entrainment(n) 
                 blob%tracer(n)     = blob%tracer(n) + detrainment(n) + entrainment(n) 
                 EL_diag(n)%detrainment(i,j,k) = EL_diag(n)%detrainment(i,j,k) - detrainment(n) 
                 EL_diag(n)%entrainment(i,j,k) = EL_diag(n)%entrainment(i,j,k) + entrainment(n)
               enddo

              if (blob%mass>small_mass) then
                 do n=1,num_prog_tracers
                    blob%field(n) = blob%tracer(n)/blob%mass
                 enddo

                 blob%density = density(blob%field(index_salt), & 
                                        blob%field(index_temp), &
                                        Dens%pressure_at_depth(i,j,k))
                 blob%densityr = 1./blob%density
                 
                 if (vert_coordinate_class == DEPTH_BASED) then
                    blob%volume = blob%mass * rho0r 
                 else
                    blob%volume = blob%mass * blob%densityr
                 endif
              endif
           endif

           ! If a blob is of zero mass, we are only interested in its detrainment
           ! properties.  So, to stop any confusion, we kill it after its 
           ! tendencies have been added to the gridded variables.
           if (blob%mass>0.0) then
              call insert_blob(blob, head)
           else
              call free_blob_memory(blob)
           endif

        endif
        nullify(blob)
     enddo
  endif

end subroutine unpackbuffer
! </SUBROUTINE>  NAME="unpackbuffer"

!######################################################################
! <SUBROUTINE NAME="increase_buffer">
!
! <DESCRIPTION>
! Increases the buffer size for sending blobs from one PE to another.
! </DESCRIPTION>
!
subroutine increase_buffer(buffer, newnum)
  type(blob_buffer_type), pointer :: buffer
  integer, intent(in) :: newnum

  ! local variables
  type(blob_buffer_type), pointer :: new_buffer
  integer :: newbuffsize, m

  newbuffsize = buffer%size
  allocate(new_buffer)
  m=1
  do while (newnum>newbuffsize)
     newbuffsize = buffer%size + m*delta_buffer
     m=m+1
  enddo

  allocate(new_buffer%integer_buff(int_buff_size,newbuffsize))
  allocate(new_buffer%real_buff(   rea_buff_size,newbuffsize))

  if (bitwise_reproduction) then
     allocate(new_buffer%history_integer_buff(hist_int_buff_size, 0:total_nsp1, newbuffsize))
     allocate(new_buffer%history_real_buff( 0:hist_rea_buff_size, 0:total_nsp1, newbuffsize))
  else
     allocate(new_buffer%history_real_buff( 0:hist_rea_buff_size, 1, newbuffsize))
  endif

  call clear_buffer(new_buffer)
  
  new_buffer%integer_buff(:,1:buffer%size) = buffer%integer_buff(:,:)
  new_buffer%real_buff(   :,1:buffer%size) = buffer%real_buff(:,:)

  if (bitwise_reproduction) new_buffer%history_integer_buff(:,:,1:buffer%size) = buffer%history_integer_buff(:,:,:)
  new_buffer%history_real_buff(   :,:,1:buffer%size) = buffer%history_real_buff(:,:,:)

  new_buffer%pe = buffer%pe
  new_buffer%size = newbuffsize
  new_buffer%numblobs = buffer%numblobs

  deallocate(buffer)

  nullify(buffer)
  buffer=>new_buffer

end subroutine increase_buffer
! </SUBROUTINE>  NAME="increase_buffer"

!######################################################################
! <SUBROUTINE NAME="send_buffer">
!
! <DESCRIPTION>
! Sends a buffer to an adjoining PE
! </DESCRIPTION>
!
subroutine send_buffer(buffer)
  type(blob_buffer_type), pointer :: buffer
  integer :: nb, n_frac_steps
  if (buffer%pe /= NULL_PE) then
     call mpp_send(buffer%numblobs, plen=1, to_pe=buffer%pe)
     if (buffer%numblobs>0) then
        ! We cycle through the blobs because sending them en-mass can cause the
        ! program to hang.
        do nb=1,buffer%numblobs
           call mpp_send(buffer%integer_buff(:,nb), int_buff_size, buffer%pe)
           call mpp_send(buffer%real_buff(:,nb),    rea_buff_size, buffer%pe)

           n_frac_steps = buffer%integer_buff(7,nb)
           if (bitwise_reproduction) then
              call mpp_send(buffer%history_integer_buff(:,0:n_frac_steps,nb), hist_int_buff_size*(n_frac_steps+1), buffer%pe)
              call mpp_send(buffer%history_real_buff(:,0:n_frac_steps,nb),    (1+hist_rea_buff_size)*(n_frac_steps+1), buffer%pe)
           else
              call mpp_send(buffer%history_real_buff(:,1,nb), 1+hist_rea_buff_size, buffer%pe)
           endif
        enddo
     endif
  endif
end subroutine send_buffer
! </SUBROUTINE>  NAME="send_buffer"

!######################################################################
! <SUBROUTINE NAME="receive_buffer">
!
! <DESCRIPTION>
! Receives a buffer from an adjoining PE
! </DESCRIPTION>
!
subroutine receive_buffer(buffer)
  type(blob_buffer_type), pointer :: buffer
  integer :: nb, incoming, n_frac_steps
  if (buffer%pe /= NULL_PE) then
     call mpp_recv(incoming, glen=1, from_pe=buffer%pe)
     if (incoming>0) then
        if (incoming>buffer%size) then
           call increase_buffer(buffer,incoming)
        endif
        do nb=1,incoming
           call mpp_recv(buffer%integer_buff(:,nb), int_buff_size, buffer%pe)
           call mpp_recv(buffer%real_buff(:,nb),    rea_buff_size, buffer%pe)

           n_frac_steps = buffer%integer_buff(7,nb)
           if (bitwise_reproduction) then
              call mpp_recv(buffer%history_integer_buff(:,0:n_frac_steps,nb), hist_int_buff_size*(n_frac_steps+1), buffer%pe)
              call mpp_recv(buffer%history_real_buff(:,0:n_frac_steps,nb), (1+hist_rea_buff_size)*(n_frac_steps+1), buffer%pe)
           else
              call mpp_recv(buffer%history_real_buff(:,1,nb), 1+hist_rea_buff_size,buffer%pe)
           endif
        enddo
     endif
     buffer%numblobs = incoming
  endif
end subroutine receive_buffer
! </SUBROUTINE>  NAME="receive_buffer"

!######################################################################
! <SUBROUTINE NAME="clear_buffer">
!
! <DESCRIPTION>
! Clears the contents of a buffer
! </DESCRIPTION>
!
subroutine clear_buffer(buffer)
  type(blob_buffer_type), pointer :: buffer
  buffer%numblobs          = 0
  if (buffer%pe /= NULL_PE) then
     buffer%integer_buff(:,:) = 0
     buffer%real_buff(:,:)    = 0.0
     if(allocated(buffer%history_integer_buff)) buffer%history_integer_buff(:,:,:) = 0
     if(allocated(buffer%history_real_buff))    buffer%history_real_buff(:,:,:)    = 0.0
  endif
  !note, we do not clear buffer%size
end subroutine clear_buffer
! </SUBROUTINE>  NAME="clear_buffer"

end module ocean_blob_dynamic_bottom_mod
