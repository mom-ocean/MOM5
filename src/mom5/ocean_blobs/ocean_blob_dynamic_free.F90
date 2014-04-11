module ocean_blob_dynamic_free_mod
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This module runs the dynamic free blob implementation of the embedded
! Lagrangian blob framework.  The module forms new dynamic free blobs,
! integrates the properties of existing blobs, and handles the transfer
! of bottom blobs to free blobs.
!</OVERVIEW>
!
!<DESCRIPTION>
! Free blobs are formed using the subroutine blob_dynamic_free_implicit,
! which is called from the blob driver module.  Free blobs must be formed 
! implicitly in time so that the surface forcing has already been applied.
! 
! The properties of free blobs are also integrated in this module, that is,
! position, velocity, mass and tracer content.  Position and velocity are
! integrated using an adaptive step Runge-Kutta scheme.  There are several
! schemes available of varying order.
!
! The module also receives blobs that are transferring from the bottom
! blob dynamic regime to the free blob regime (i.e. they have separated
! from the bottom boundary).
!</DESCRIPTION>
!
!<INFO>
!
! <REFERENCE>
!  Bogacki, P., Shampine, L.F., (1989) A 3(2) pair of Runge-Kutta formulas.
!  Applied Mathematical Letters 2(4), 321-325.
! </REFERENCE>
!
! <REFERENCE>
!  Cash, J.R., Karp, A.H. (1990) A variable order Runge-Kutta method for 
!  initial value problems with rapidly varying right-hand sides.
!  ACM Transactions on Mathematical Software 16(3), 201-222.
! </REFERENCE>
!
! <REFERENCE>
!  Griffies, S.M., Harrison, M.J., Pacanowski, R.C., Rosati, A. (2004)
!  A Technical Guide to MOM4.  GFDL Ocean Group Technical Report No. 5.
!  NOAA/Geophysical Fluid Dynamics Laboratory.
! </REFERENCE>
!
! <REFERENCE>
!  Marshall, J., Schott, F. (1999) Open-ocean convection: Observations, theory,
!  and models.  Reviews of Geophysics 37(1), 1-64.
! </REFERENCE>
!
!<NAMELIST NAME="ocean_blob_dynamic_free_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module.
!  Default is use_this_module=.false.
!  </DATA>
!
!  <DATA NAME="rayleigh_drag_new" TYPE="real">
!  Rayleigh drag coefficient (1/s) for new blobs that
!  are formed due to the vertical instability
!  criterion.  Corresponds to alpha in the notes.  
!  Default is rayleigh_drag_new=1.0e-5
!  </DATA>
!
!  <DATA NAME="rayleigh_drag_bot" TYPE="real">
!  Rayleigh drag coefficient (1/s) for bottom blobs
!  that become free blobs.  Corresponds to alpha in 
!  the notes.
!  Default is rayleigh_drag_bot=1.0e-7
!  </DATA>
!
!  <DATA NAME="update_method" TYPE="character">
!  Decide which method to use to integrate the 
!  blobs.  Choices are 'BS_RK3(2)' or 'CK_RK5(4)'
!  for the Bogaki-Shampine or Cash-Karp methods
!  respectively.
!  Default is update_method='CK_RK5(4)
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
!  <DATA NAME="size_fact" TYPE="real">
!  An Adjustment for blob size, 0<size_fact<=1.0
!  Corresponds to Lambda in the notes.
!  Default is size_fact=1.0
!  </DATA>
!
!  <DATA NAME="det_param" TYPE="real">
!  The detrainment parameter (kg m^2/s). 
!  Corresponds to Gamma in the notes.
!  Default is det_param=5.0e-8
!  </DATA>
!
!  <DATA NAME="max_detrainment" TYPE="real">
!  The Maximum allowable detrainment velocity (m/s).
!  Default is max_detrainment=1.0e-3
!  </DATA>
!
!  <DATA NAME="bv_freq_threshold" TYPE="real">
!  The buoyancy frequency threshold at which 
!  the scheme will start to create blobs, i.e.
!  blobs will be formed when N^2<bv_freq_threshold
!  Default is bv_freq_threshold=-1.0e-15
!  </DATA>
!
!  <DATA NAME="full_N2" TYPE="logical">
!  Whether to use the buoyancy frequency calculated
!  from the combined E and L system (true) or, from
!  the E system only (false).
!  Default is full_N2=.true.
!  </DATA>
!
!  <DATA NAME="large_speed" TYPE="real">
!  A value for error checking.  If the speed of a
!  blob exceeds large_speed in any of x,y,z then
!  a warning flag is raised.
!  Default is large_speed=10.0
!  </DATA>
!
!</NAMELIST>
!</INFO>
!
use constants_mod,    only: deg_to_rad, epsln, pi
use diag_manager_mod, only: register_diag_field, send_data

use fms_mod,         only: stdout, stdlog, open_namelist_file, WARNING, FATAL
use fms_mod,         only: mpp_error, check_nml_error, close_file
use mpp_mod,         only: mpp_send, mpp_recv, NULL_PE, mpp_sum, mpp_sync_self
use mpp_mod,         only: mpp_set_current_pelist
use mpp_mod,         only: CLOCK_LOOP, CLOCK_ROUTINE, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_domains_mod, only: mpp_update_domains, mpp_global_sum

use ocean_blob_util_mod,  only: E_and_L_totals, count_blob, reallocate_interaction_memory
use ocean_blob_util_mod,  only: insert_blob, allocate_interaction_memory, interp_tcoeff, interp_ucoeff
use ocean_blob_util_mod,  only: check_ijcell, check_kcell, kill_blob, free_blob_memory
use ocean_blob_util_mod,  only: blob_delete, unlink_blob, check_cyclic
use ocean_density_mod,    only: density, buoyfreq2
use ocean_parameters_mod, only: rho0, rho0r, grav, omega_earth
use ocean_parameters_mod, only: PRESSURE_BASED, DEPTH_BASED, onehalf, onethird, twothirds, onefourth
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_lagrangian_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_blob_type, ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_density_type, ocean_velocity_type, blob_diag_type
use ocean_types_mod,      only: ocean_external_mode_type, blob_grid_type, ocean_adv_vel_type
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3

implicit none

private

type(ocean_grid_type),   pointer :: Grd  => NULL()
type(ocean_domain_type), pointer :: Dom  => NULL()
type(ocean_domain_type), pointer :: Bdom => NULL() !A domain variable with halo=2
type(blob_grid_type),    pointer :: Info => NULL()

! Module wide variables that are inherited/derived from the main model
real, allocatable, dimension(:,:,:) :: umask    !umask, but, halo=2
real, allocatable, dimension(:,:,:) :: tmask    !tmask, but, halo=2
real, allocatable, dimension(:,:)   :: datdtime !Grd%dat*dtime
integer :: vert_coordinate_class
integer :: vert_coordinate
real :: dtime
real :: dtime_yr

! Useful variables
real :: two_omega
real :: p5_dtime
real :: grav_dtime
real :: det_factor

!RK coefficients
real, dimension(:,:), allocatable :: beta
real, dimension(:),   allocatable :: lgamma
real, dimension(:),   allocatable :: hgamma
real :: order
real :: orderp1_r

!Useful indexes
integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1
integer :: isd,ied,jsd,jed
integer :: isc,iec,jsc,jec
integer :: isg,ieg,jsg,jeg
integer :: nk
integer :: isbd, iebd, jsbd, jebd
integer :: method
integer :: iq(4), jq(4)
integer :: ns, nsp1, total_ns, total_nsp1, total_nsp2

integer :: id_clock_dyn_update
integer :: id_clock_part_cycle
integer :: id_clock_rk
integer :: id_clock_updatevars
integer :: id_clock_findEvars

!Diagnostics
integer, allocatable, dimension(:) :: id_tracer_new
logical :: used
integer :: id_bot_to_free
integer :: id_new_blobs

!Buffers, for sending blobs between compute domains
type, private :: blob_buffer_type
   integer :: numblobs
   integer :: size
   integer :: pe
   integer, allocatable, dimension(:,:) :: integer_buff
   real,    allocatable, dimension(:,:) :: real_buff
   integer, allocatable, dimension(:,:,:) :: history_integer_buff
   real,    allocatable, dimension(:,:,:) :: history_real_buff
end type blob_buffer_type

integer :: rea_buff_size      !size of the real buffer
integer :: int_buff_size      !size of the integer buffer
integer :: hist_rea_buff_size !size of the real history buffer
integer :: hist_int_buff_size !size of the integer history buffer

!The buffers themselves
type(blob_buffer_type), pointer :: Ebuffer_out,  Ebuffer_in
type(blob_buffer_type), pointer :: Wbuffer_out,  Wbuffer_in
type(blob_buffer_type), pointer :: Nbuffer_out,  Nbuffer_in
type(blob_buffer_type), pointer :: Sbuffer_out,  Sbuffer_in
type(blob_buffer_type), pointer :: NEbuffer_out, NEbuffer_in
type(blob_buffer_type), pointer :: NWbuffer_out, NWbuffer_in
type(blob_buffer_type), pointer :: SEbuffer_out, SEbuffer_in
type(blob_buffer_type), pointer :: SWbuffer_out, SWbuffer_in

integer, parameter :: delta_buffer = 25   ! Size by which to increment the buffer

public blob_dynamic_free_init
public blob_dynamic_free_implicit
public blob_dynamic_free_update
public transfer_bottom_to_free
public blob_dynamic_free_end

! Module wide variables that are controlled by the ocean_blob_nml and ocean_blob_diag_nml
logical :: debug_this_module
logical :: really_debug
logical :: bitwise_reproduction
logical :: module_is_initialized=.false.
logical :: blob_diag
real    :: small_mass

! namelist defaults
logical           :: use_this_module   = .false.
real              :: rayleigh_drag_new = 1.0e-5  ! Rayleigh drag coefficient (1/s) for newly formed blobs
real              :: rayleigh_drag_bot = 1.0e-7  ! Rayleigh drag coefficient (1/s) for bottom blobs that become free blobs
character(len=10) :: update_method     = 'CK_RK5(4)'
real              :: rel_error         = 0.01    ! Relative error for the RK scheme (non-dimensional, 0<rel_error<=1.0)
real              :: safety_factor     = 0.8     ! Safety factor for the RK scheme (non-dimensional, 0<safety_factor<=1.0)
real              :: minstep           =  9.0    ! Minimum step for a blob (seconds)
real              :: first_step        = 50.0    ! Size of first step of a blob
real              :: size_fact         = 1.0     ! Adjustment for blob size (non-dimensional, 0<size_fact<=1.0)
real              :: det_param         = 5.e-8   ! Detrainment parameter (Gamma; kg m**2 / s)
real              :: max_detrainment   = 1.e-3   ! Maximum detrainment velocity (m/s)
real              :: bv_freq_threshold = -1.e-15 ! The stability threshold at which to form blobs (1/s**2)
logical           :: full_N2           =.true.
real              :: large_speed       = 10.0    ! For stability checking (m/s)

namelist /ocean_blob_dynamic_free_nml/ use_this_module, rayleigh_drag_new, &
     rayleigh_drag_bot, update_method, rel_error, safety_factor, minstep,  &
     first_step, size_fact, det_param, max_detrainment, bv_freq_threshold, &
     full_N2, large_speed

contains
!#######################################################################
! <SUBROUTINE NAME="blob_dynamic_free_init">
!
! <DESCRIPTION>
! Initialises the dynamic free blobs by checking the namelist and also
! inherited namelists (from ocean_blob_nml).  Also sets up some useful
! constants, allocates memory to special halo=2 masks and sets up 
! the blob buffers for sending blobs from one PE to another.
! </DESCRIPTION>
!
subroutine blob_dynamic_free_init(Grid, Domain, Time, T_prog, Blob_domain,         &
                                  PE_info, debug, big_debug, bitwise, num_tracers, &
                                  itemp, isalt, dtimein, ver_coord_class,          &
                                  ver_coord,  blob_diagnostics, free_minstep,      &
                                  free_total_ns, smallmass, use_dyn_fre)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  type(ocean_domain_type),      intent(in), target :: Blob_domain
  type(blob_grid_type),         intent(in), target :: PE_info
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
  real,    intent(out) :: free_minstep
  integer, intent(out) :: free_total_ns
  real,    intent(in)  :: smallmass
  logical, intent(out) :: use_dyn_fre

  real, parameter :: secs_in_year_r = 1.0 / (86400.0 * 365.25)

  integer :: n
  integer :: ioun, ierr, io_status
  integer :: stdoutunit,stdlogunit
  character(32) :: myname, myunit

  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
         '==>Error in ocean_blob_dynamic_free_mod (ocean_blob_dynamic_free_init):' &
         //' module already initialized')
  endif 

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_blob_dynamic_free_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_blob_dynamic_free_nml)  
  write (stdlogunit,ocean_blob_dynamic_free_nml)
  ierr = check_nml_error(io_status,'ocean_blob_dynamic_free_nml')
  call close_file (ioun)

  module_is_initialized = .true.
  
  ! The way things are formulated, we cannot run the bottom blobs
  ! without running free blobs.  But, we can stop free blobs being
  ! formed by the vertical instability condition.  So, if dynamic
  ! bottom blos are selected, we overwrite use_this_module for 
  ! the dynamic free blobs, however, we ensure that the formation
  ! of free blobs by vertical instabily is not invoked.
  use_dyn_fre = use_this_module

  id_bot_to_free = register_diag_field('ocean_model', 'bot_to_free', Time%model_time, &
       'blobs separating from the bottom', 'number of blobs')

  if (.not. use_this_module) return

  id_new_blobs   = register_diag_field('ocean_model', 'new_free_blobs', Time%model_time, &
       'new free blobs', 'number of blobs')

  blob_diag            = blob_diagnostics
  free_minstep         = minstep
  index_temp           = itemp
  index_salt           = isalt
  num_prog_tracers     = num_tracers
  debug_this_module    = debug
  really_debug         = big_debug
  bitwise_reproduction = bitwise
  small_mass           = smallmass

  vert_coordinate_class = ver_coord_class
  vert_coordinate       = ver_coord

  Grd  => Grid
  Dom  => Domain
  Bdom => Blob_domain
  Info => PE_info

  isd=Dom%isd; ied=Dom%ied; jsd=Dom%jsd; jed=Dom%jed
  isc=Dom%isc; iec=Dom%iec; jsc=Dom%jsc; jec=Dom%jec
  isg=Dom%isg; ieg=Dom%ieg; jsg=Dom%jsg; jeg=Dom%jeg
  nk = Grd%nk

  isbd=Bdom%isd; iebd=Bdom%ied; jsbd=Bdom%jsd; jebd=Bdom%jed 

  !1: (i-1,j-1), 2:(i-1,j), 3:(i,j), 4:(i,j-1)
  iq(1)=-1; jq(1)=-1
  iq(2)=-1; jq(2)= 0
  iq(3)= 0; jq(3)= 0
  iq(4)= 0; jq(4)=-1 

  allocate(umask(isbd:iebd,jsbd:jebd,1:nk))
  allocate(tmask(isbd:iebd,jsbd:jebd,0:nk)) !we need tmask(0:nk) for the W grid
  umask(isc:iec,jsc:jec,1:nk) = Grd%umask(isc:iec,jsc:jec, 1:nk)
  tmask(isc:iec,jsc:jec,1:nk) = Grd%tmask(isc:iec,jsc:jec, 1:nk)
  tmask(isc:iec,jsc:jec,0)    = Grd%tmask(isc:iec,jsc:jec, 1)
  call mpp_update_domains(umask(:,:,:), Bdom%domain2d)
  call mpp_update_domains(tmask(:,:,:), Bdom%domain2d)
  ! special treatment for boundaries if they are NULL_PEs.
  ! Make sure that umask==0 and tmask==0
  if (Info%pe_S==NULL_PE) then
     umask(:,jsbd:jsd,:) = 0.0
     tmask(:,jsbd:jsd,:) = 0.0
  endif
  if (Info%pe_N==NULL_PE) then
     umask(:,jed:jebd,:) = 0.0
     tmask(:,jed:jebd,:) = 0.0
  endif
  if (Info%pe_E==NULL_PE) then
     umask(ied:iebd,:,:) = 0.0
     tmask(ied:iebd,:,:) = 0.0
  endif
  if (Info%pe_W==NULL_PE) then
     umask(isbd:isd,:,:) = 0.0
     tmask(isbd:isd,:,:) = 0.0
  endif

  allocate( datdtime(isd:ied,jsd:jed) )
  datdtime(:,:) = Grd%dat(:,:)*dtimein

  dtime      = dtimein
  dtime_yr   = dtime*secs_in_year_r
  grav_dtime = grav*dtime
  p5_dtime   = onehalf*dtime !for convenience
  two_omega  = 2.0*omega_earth

  solver_method: select case(trim(update_method))
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
          '==>Error in ocean_blob_dynamic_free_mod (ocean_blob_dynamic_free_init):' &
          //' invalid solver chosen ('//update_method//').  Check update_method in namelist.'
     call mpp_error(FATAL,&
          '==>Error in ocean_blob_dynamic_free_mod (ocean_blob_dynamic_free_init):' &
          //' invalid solver chosen ('//update_method//').  Check update_method in namelist.')
  endselect solver_method

  if (method==0) minstep=dtime
  ns            = ubound(lgamma,1) !number of partial steps in a fractional step
  nsp1          = ns+1             !number of partial steps in a fractional step+1
  total_ns      = ceiling(dtime/minstep)
  total_nsp1    = total_ns+1
  total_nsp2    = total_nsp1+1
  free_total_ns = total_ns

  ! This collects all the constant terms together in the expression
  ! for the rate of change of mass.
  ! dm/dt = rhoL A D (A=blob surface area, D=detrainment rate)
  !       = -rhoL A Gamma/|rhoL - rhoE|
  !       = -m**2/3 Gamma rhoL (36pi)**1/3 / (rhoL**2/3 |rhoL-rhoE|)
  !
  ! In the Boussinesq case rhoL=rho0 (outside |rhoL-rhoE|)
  ! we also need to take into account that this is done for each PARTIAL step.
  if (vert_coordinate_class==DEPTH_BASED) then
     ! here: det_factor = Gamma (rho0*36pi)**1/3 / number of partial steps
     det_factor = det_param*( (rho0*36*pi)**onethird )
  else !PRESSURE_BASED
     ! here: det_factor = Gamma (36pi)**1/3
     det_factor = det_param*( (36*pi)**onethird )
  endif

  ! Allocate the buffers
  ! Things that dictate the size of the real buffer are:
  ! tracer content (num_prog_tracer), change in mass (1),
  ! change in tracer (num_prog_tracer; only needed for bitwise_reproduction),
  ! velocity(3), position (3),
  ! step size (1), blob time (1), mass (1),  stretching function (2) 
  ! drag (1), dmass(1), gprime (1) and age (1).  
  rea_buff_size = 2*num_prog_tracers+15
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

  id_clock_dyn_update = mpp_clock_id('(Ocean dyn. free blob: update)    ',grain=CLOCK_ROUTINE)
  id_clock_part_cycle = mpp_clock_id('(Ocean dyn. free blob: part step) ',grain=CLOCK_LOOP)
  id_clock_rk         = mpp_clock_id('(Ocean dyn. free blob: RK scheme) ',grain=CLOCK_LOOP)
  id_clock_updatevars = mpp_clock_id('(Ocean dyn. free blob: updatevars)',grain=CLOCK_LOOP)
  id_clock_findEvars  = mpp_clock_id('(Ocean dyn. free blob: findEvars) ',grain=CLOCK_LOOP)

  allocate(id_tracer_new(num_prog_tracers))
  do n=1,num_prog_tracers
     if (n==index_temp) then
        myname = 'heat'
        myunit = 'J'
     else
        myname = T_prog(n)%name
        myunit = 'kg'
     endif
     id_tracer_new(n) = register_diag_field('ocean_model',   'new_free_blob_'//trim(myname), &
          Grd%tracer_axes(1:3), Time%model_time, trim(myname)//' transferred from the E to the L system via new free blobs', myunit)
  enddo

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that deallocates memory from a buffer    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine allocate_buffer(buffer, pe)
    type(blob_buffer_type), pointer    :: buffer
    integer,                intent(in) :: pe

    allocate(buffer)
    buffer%pe = pe
    buffer%numblobs = 0
    if (buffer%pe /= NULL_PE) then
       buffer%size = delta_buffer
       
       allocate(buffer%integer_buff(int_buff_size, delta_buffer))
       allocate(buffer%real_buff(   rea_buff_size, delta_buffer))
       
       if (bitwise_reproduction) then
          allocate(buffer%history_integer_buff(hist_int_buff_size, 0:total_nsp1, delta_buffer))
          allocate(buffer%history_real_buff( 0:hist_rea_buff_size, 0:total_nsp1, delta_buffer))
       else
          allocate(buffer%history_real_buff(0:hist_rea_buff_size, 1, delta_buffer))
       endif
    endif
  end subroutine allocate_buffer

end subroutine blob_dynamic_free_init
! </SUBROUTINE>  NAME="blob_dynamic_free_init"


!#######################################################################
! <SUBROUTINE NAME="blob_dynamic_free_implicit">
!
! <DESCRIPTION>
! Initialises dynamic blobs in vertical statically unstable regions.
! Due to the instability condition, blobs should be formed after the
! surface forcing has been applied (which is a major source of 
! instability in the water column).  The surface forcing is applied
! implicitly in time in MOM, therefore, we must form blobs implicitly
! in time.
!
! If N^2<bv_freq_threshold, then, two blobs are formed.  One rising
! and one sinking.  The rising blobs is destroyed immediately (after
! it has been moved up one cell) and its properties returned to the E
! system.  The sinking blob is added to a linked list, and its 
! properties integrated at a later time step.
! </DESCRIPTION>
!
subroutine blob_dynamic_free_implicit(Time, Thickness, T_prog, Dens, Adv_vel, &
                                      Velocity, head, blob_counter, EL_diag)
                                          
  type(ocean_time_type),                                intent(in)    :: Time
  type(ocean_thickness_type),                           intent(inout) :: Thickness
  type(ocean_prog_tracer_type),                         intent(inout) :: T_prog(:)
  type(ocean_density_type),                             intent(in)    :: Dens
  type(ocean_adv_vel_type),                             intent(in)    :: Adv_vel
  type(ocean_velocity_type),                            intent(in)    :: Velocity
  type(ocean_blob_type), pointer                                      :: head
  integer,               dimension(isc:iec,jsc:jec,nk), intent(inout) :: blob_counter
  type(blob_diag_type),  dimension(0:num_prog_tracers), intent(inout) :: EL_diag

  real, dimension(isd:ied,jsd:jed,nk) :: bvfreq2
  integer :: i, j, k, kp1, tau, taup1, n, nblobs
  real    :: mass_blob, rhodzt_blob, rhodztk, rhodztkp1
  real    :: rho_dzt_old, rho_dzt_dat_r
  real    :: small
  real    :: wE
  type(ocean_blob_type), pointer :: sink
  type(ocean_blob_type), pointer :: rise

  small  = 1e-3
  nblobs = 0

  ! No need for use_this_module test, as we test for 
  ! use_dyn_fre in ocean_blob_implicit

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
         '==>Error in ocean_blob_static_free_mod (Lagrangian blob model): '&
         //'module needs to be initialized')
  endif 

  taup1 = Time%taup1
  tau   = Time%tau

  ! Begin by computing the square of the Brunt-Vaisalla Frequency
  if (full_N2) then
     ! Use the combied L and E system for calculating N**2
     do k=1,nk
        wrk1(:,:,k) = ( Thickness%rho_dzt(:,:,k,taup1)*T_prog(index_salt)%field(:,:,k,taup1) &
                        + Grd%datr(:,:)*T_prog(index_salt)%sum_blob(:,:,k,taup1)      )/Thickness%rho_dztT(:,:,k,taup1)
        wrk2(:,:,k) = ( Thickness%rho_dzt(:,:,k,taup1)*T_prog(index_temp)%field(:,:,k,taup1) &
                        + Grd%datr(:,:)*T_prog(index_temp)%sum_blob(:,:,k,taup1)      )/Thickness%rho_dztT(:,:,k,taup1)
     enddo
     wrk3(:,:,:) = density(wrk1(:,:,:), wrk2(:,:,:), Dens%pressure_at_depth(:,:,:))
     bvfreq2(:,:,:) = buoyfreq2(Time, Thickness, Dens, wrk1(:,:,:), wrk2(:,:,:), wrk3(:,:,:), use_this_module) 

     ! We need the E system density for calculations below
     wrk3(:,:,:) = density(T_prog(index_salt)%field(:,:,:,taup1), T_prog(index_temp)%field(:,:,:,taup1), Dens%pressure_at_depth(:,:,:))

  else
     ! Only use the E system for calculating N**2
     wrk3(:,:,:) = density(T_prog(index_salt)%field(:,:,:,taup1), T_prog(index_temp)%field(:,:,:,taup1), Dens%pressure_at_depth(:,:,:))
     bvfreq2(:,:,:) = buoyfreq2(Time, Thickness, Dens, T_prog(index_salt)%field(:,:,:,taup1), &
                                T_prog(index_temp)%field(:,:,:,taup1), wrk3(:,:,:), use_this_module) 
  endif

  do k=1,nk-1
     kp1 = k+1
     do j=jsc,jec
        do i=isc,iec
           if (bvfreq2(i,j,k) < bv_freq_threshold) then
              ! some shorthand variables for convenience
              rhodztk   = Thickness%rho_dzt(i,j,k,  taup1)
              rhodztkp1 = Thickness%rho_dzt(i,j,kp1,taup1)

              ! compute the mass per unit of the two blobs.  Note that
              ! dat is the same for both k and kp1;
              ! mass=(dat*rhodztk * dat*rhodztkp1)/(dat(rhodztk + rhodztkp1))
              !     =dat*rhodztk*rhodztkp1/(rhodztk + rhodztkp1)
              ! mass per unit area = rhodztk*rhodztkp1/(rhodztk + rhodztkp1)
              ! size_fact is a non-dimensional namelist parameter that is 
              ! scales the size of the blob.  0.0<size_fact<=1.0
              rhodzt_blob = size_fact*(rhodztk * rhodztkp1)/(rhodztk + rhodztkp1)

              mass_blob = rhodzt_blob * Grd%dat(i,j)

              ! form the sinking dynamic blob
              allocate(sink)
              allocate(sink%tracer(num_prog_tracers))
              allocate(sink%field(num_prog_tracers))
              
              sink%i         = i
              sink%j         = j
              sink%k         = k
              sink%blob_time = 0.0
              sink%mass      = mass_blob

              ! Assign tracer content and tracer concentration to the blob
              do n=1,num_prog_tracers
                 ! sink%tracer is blob_mass*[concentration]
                 sink%tracer(n) = sink%mass*T_prog(n)%field(i,j,k,taup1)
                 ! even though sink%field is initially the concentration of the grid cell
                 ! we need to calculate it like this to maintain bitwise agreement
                 sink%field(n)  = sink%tracer(n)/sink%mass
              enddo

              ! Calculate the density of the blob at kp1
              sink%density  = density(sink%field(index_salt), &
                                      sink%field(index_temp), &
                                      Dens%pressure_at_depth(i,j,kp1))
              sink%densityr = 1.0/sink%density

              ! the intial vertical velocity
              ! We use rho0r whether it is depth or pressure based coordinates to calculate
              ! the E system vertical velocity.  This represents an error of about 2% in 
              ! pressure based coordinates.
              wE = rho0r * Adv_vel%wrho_bt(i,j,k)
              sink%v(3) = wE - grav_dtime*(sink%density-wrk3(i,j,kp1))*sink%densityr

              if (sink%v(3) > 0.0) then
                 ! the vertical velocity is positive and therefore the blob will want to 
                 ! rise instead of sink.  this can be because of either
                 ! 1/ the pseudo-non-hydrostatic term is positive (i.e. sink%density<rhoE)
                 ! 2/ wE > abs(pseudo-non-hydrostatic term).
                 ! If a blob does not want to sink, then it wants to stay in the original
                 ! grid cell, and there is thus no point in forming it, so
                 ! we dont bother forming it and accept that the water column will remain
                 ! unstable.
                 call free_blob_memory(sink)
                 cycle
              endif

              ! count the blob
              call count_blob(sink, blob_counter)

              ! insert the blob into the linked list
              call insert_blob(sink, head)
              nblobs = nblobs+1

              sink%age = 0.0
              
              ! take the sink blob away from the E systems rho_dzt
              Thickness%rho_dzt(i,j,k,taup1) = Thickness%rho_dzt(i,j,k,taup1) &
                                              - sink%mass * Grd%datr(i,j)

              ! update rho_dztr
              Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)

              ! We do not need to adjust the E system field for the k grid cell at this
              ! stage.

              ! we take the average of the four surrounding velocity grid cells for the 
              ! initial horizontal velocity of sink
              sink%v(1) = ( Velocity%u(i  ,j  ,k,1,tau) + Velocity%u(i-1,j  ,k,1,tau) &
                          + Velocity%u(i  ,j-1,k,1,tau) + Velocity%u(i-1,j-1,k,1,tau) &
                          )*onefourth

              sink%v(2) = ( Velocity%u(i  ,j  ,k,2,tau) + Velocity%u(i-1,j  ,k,2,tau) &
                          + Velocity%u(i  ,j-1,k,2,tau) + Velocity%u(i-1,j-1,k,2,tau) &
                          )*onefourth
              
              ! find the intial longitude and latitude
              sink%lon = Grd%xt(i,j) + dtime*sink%v(1)/Grd%h1t(i,j)
              sink%lat = Grd%yt(i,j) + dtime*sink%v(2)/Grd%h2t(i,j)

              sink%new  = .true.
              sink%step = first_step
              
              sink%model_steps = 0
              sink%nsteps      = 0

              sink%drag = rayleigh_drag_new
              
              ! Now create the rising blob
              allocate(rise)
              allocate(rise%tracer(num_prog_tracers))

              rise%mass = mass_blob
              ! Assign tracer content and tracer concentration to the blob
              do n=1,num_prog_tracers
                 ! rise%tracer is blob_mass*[concentration]
                 rise%tracer(n) = rise%mass*T_prog(n)%field(i,j,kp1,taup1)
              enddo

              ! take the rise blob away from the E systems rho_dzt
              Thickness%rho_dzt(i,j,kp1,taup1) = Thickness%rho_dzt(i,j,kp1,taup1) &
                                                 - rise%mass * Grd%datr(i,j)

              ! update rho_dztr
              Thickness%rho_dztr(i,j,kp1) = 1.0/(Thickness%rho_dzt(i,j,kp1,taup1)+epsln)
              
              ! Now we swap the two blobs.  sink goes to the kp1 cell and rise goes to the
              ! k cell.
              rise%k = k
              sink%k = kp1

              ! Calculate the depth and position the native vertical coordinate of sink.
              ! We recall that z and depth have opposite signs.
              sink%depth = Thickness%depth_zwt(i,j,k) - p5_dtime*sink%v(3)

              ! If the initial position of the blob is calculated to be deeper than
              ! the bottom of the kp1 cell (in coordinate space), we adjust the initial
              ! position to ensure that the blob remains in the kp1 cell.  We also adjust
              ! the initial velocity to reflect the new imposed position.
              ! Also calculate the volume
              if (vert_coordinate_class == DEPTH_BASED) then
                 sink%st = -( Thickness%depth_swt(i,j,k) + (sink%depth-Thickness%depth_zwt(i,j,k))/Thickness%dzt_dst(i,j,kp1) )
                 if (sink%st <= -Thickness%depth_swt(i,j,kp1)) then
                    sink%st    = -Thickness%depth_swt(i,j,kp1) + small
                    sink%depth =  Thickness%depth_zwt(i,j,kp1) - small*Thickness%dzt_dst(i,j,kp1)
                    sink%v(3)  = -( sink%depth-Thickness%depth_zwt(i,j,k) )/dtime
                 endif
                 sink%volume = sink%mass * rho0r
              else !PRESSURE_BASED
                 sink%st = Thickness%depth_swt(i,j,k) - (sink%depth-Thickness%depth_zwt(i,j,k))/Thickness%dzt_dst(i,j,kp1)
                 if (sink%st >= Thickness%depth_swt(i,j,kp1)) then
                    sink%st    =  Thickness%depth_swt(i,j,kp1) - small
                    sink%depth =  Thickness%depth_zwt(i,j,kp1) + small*Thickness%dzt_dst(i,j,kp1)
                    sink%v(3)  = -( sink%depth-Thickness%depth_zwt(i,j,k) )/dtime
                 endif
                 sink%volume = sink%mass * sink%densityr
              endif
              
              ! Add the blob tracer content to the total blob tracer content
              do n=1,num_prog_tracers
                 T_prog(n)%sum_blob(i,j,kp1,taup1) = T_prog(n)%sum_blob(i,j,kp1,taup1) + sink%tracer(n)
              enddo
              
              ! Now we take care of rise, returning its properties to the E system.

              ! Transfer mass from the Lagrangian system to the Eulerian System
              rho_dzt_old = Thickness%rho_dzt(i,j,k,taup1)
              Thickness%rho_dzt(i,j,k,taup1) = Thickness%rho_dzt(i,j,k,taup1)&
                                               + rise%mass*Grd%datr(i,j)
              Thickness%rho_dztr(i,j,k) = 1.0/(Thickness%rho_dzt(i,j,k,taup1)+epsln)
              rise%mass = rise%mass - rise%mass
              
              ! transfer tracer from the Lagrangian system to the Eulerian System
              rho_dzt_dat_r = Grd%datr(i,j)*Thickness%rho_dztr(i,j,k)
              do n=1,num_prog_tracers
                 T_prog(n)%field(i,j,k,taup1) = rho_dzt_old*Thickness%rho_dztr(i,j,k)*T_prog(n)%field(i,j,k,taup1) &
                                                + rise%tracer(n)*rho_dzt_dat_r
                 rise%tracer(n) = rise%tracer(n) - rise%tracer(n)
              enddo

              call free_blob_memory(rise)
              
              ! Add in some more structure to the sinking blob.  The structure is required
              ! for the E-L system interaction.
              call allocate_interaction_memory(sink, total_ns)

              ! Diagnostics
              EL_diag(0)%new(i,j,k) = EL_diag(0)%new(i,j,k) + sink%mass
              do n=1,num_prog_tracers
                 EL_diag(n)%new(i,j,k) = EL_diag(n)%new(i,j,k) + sink%tracer(n)
              enddo

              ! Need to set this to zero for diagnostics
              sink%ent = 0.0

              nullify(sink)
              
            endif !bv frequency small?
        enddo !i
     enddo !j
  enddo!k

  call mpp_sum(nblobs)
  if (id_new_blobs>0) used = send_data(id_new_blobs, real(nblobs), Time%model_time)

end subroutine blob_dynamic_free_implicit
! </SUBROUTINE>  NAME="blob_dynamic_free_implicit"


!######################################################################
! <SUBROUTINE NAME="blob_dynamic_free_update">
!
! <DESCRIPTION>
! This routine calls the routine to update blob positions.  When
! bitwise_reproduction=.false., it also figures out when to continue
! the integration of blobs that have changed PE's.
! </DESCRIPTION>
!
subroutine blob_dynamic_free_update(Time, Thickness, T_prog, Ext_mode, Dens,      &
                                    L_system, tend_blob, blob_source, press_grad, &
                                    u, w, model_rho, free_head, bottom, &
                                    mass_in, mass_out, EL_diag, ngrnd, nsfc, ndetrn)
  type(ocean_time_type),                               intent(in)    :: Time
  type(ocean_thickness_type),                          intent(inout) :: Thickness
  type(ocean_prog_tracer_type),                        intent(inout) :: T_prog(:)
  type(ocean_external_mode_type),                      intent(inout) :: Ext_mode
  type(ocean_density_type),                            intent(in)    :: Dens
  type(ocean_lagrangian_type),                         intent(inout) :: L_system
  real,dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  real,dimension(isd:ied,jsd:jed),                     intent(inout) :: blob_source
  real, dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),       intent(in)    :: press_grad
  real, dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),       intent(in)    :: u
  real, dimension(isbd:iebd,jsbd:jebd,1:nk),           intent(in)    :: w
  real, dimension(isbd:iebd,jsbd:jebd,1:nk),           intent(in)    :: model_rho
  type(ocean_blob_type),                               pointer       :: free_head
  type(ocean_blob_type),                               pointer       :: bottom
  real, dimension(isd:ied,jsd:jed,1:nk),               intent(inout) :: mass_in
  real, dimension(isd:ied,jsd:jed,1:nk),               intent(inout) :: mass_out
  type(blob_diag_type), dimension(0:num_prog_tracers), intent(inout) :: EL_diag
  integer,                                             intent(inout) :: ngrnd
  integer,                                             intent(inout) :: nsfc
  integer,                                             intent(inout) :: ndetrn

  type(ocean_blob_type), pointer :: prev, this, next, buffer_head
  integer :: check_buffers
  integer :: i,j,k,n
  integer :: stdoutunit

  if (.not. use_this_module) return
  
  nullify(bottom)
  nullify(this)

  ! Clear buffers for sending and receiving free blobs
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
  if (associated(free_head)) then
     this=>free_head
     timecycle: do
        this%blob_time = 0.0
        this=>this%next
        if(.not.associated(this)) exit timecycle
     enddo timecycle
  endif

  call mpp_clock_begin(id_clock_dyn_update)
  call dynamic_update(Time, Thickness, Ext_mode, L_system, Dens,     &
                      T_prog(:), tend_blob, blob_source, press_grad, &
                      u, w, model_rho, free_head, bottom,            &
                      mass_in, mass_out, EL_diag(:), ngrnd, nsfc, ndetrn)
  call mpp_clock_end(id_clock_dyn_update)

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
  
     call receive_buffer( Wbuffer_in); call unpackbuffer(Time,  Wbuffer_in, free_head, Dens)
     call receive_buffer( Ebuffer_in); call unpackbuffer(Time,  Ebuffer_in, free_head, Dens)
     call receive_buffer( Nbuffer_in); call unpackbuffer(Time,  Nbuffer_in, free_head, Dens)
     call receive_buffer( Sbuffer_in); call unpackbuffer(Time,  Sbuffer_in, free_head, Dens)
     call receive_buffer(SWbuffer_in); call unpackbuffer(Time, SWbuffer_in, free_head, Dens)
     call receive_buffer(SEbuffer_in); call unpackbuffer(Time, SEbuffer_in, free_head, Dens)
     call receive_buffer(NWbuffer_in); call unpackbuffer(Time, NWbuffer_in, free_head, Dens)
     call receive_buffer(NEbuffer_in); call unpackbuffer(Time, NEbuffer_in, free_head, Dens)

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

     ! Send, receive and unpack buffers for free blobs
     call send_buffer( Ebuffer_out)
     call send_buffer( Wbuffer_out)
     call send_buffer( Sbuffer_out)
     call send_buffer( Nbuffer_out)
     call send_buffer(NEbuffer_out)
     call send_buffer(NWbuffer_out)
     call send_buffer(SEbuffer_out)
     call send_buffer(SWbuffer_out)

     buffer_head=> NULL()
     call receive_buffer( Wbuffer_in) 
     call unpackbuffer(Time,  Wbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer( Ebuffer_in)
     call unpackbuffer(Time,  Ebuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer( Nbuffer_in)
     call unpackbuffer(Time,  Nbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer( Sbuffer_in)
     call unpackbuffer(Time,  Sbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(SWbuffer_in)
     call unpackbuffer(Time, SWbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(SEbuffer_in)
     call unpackbuffer(Time, SEbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(NWbuffer_in)
     call unpackbuffer(Time, NWbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(NEbuffer_in) 
     call unpackbuffer(Time, NEbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     
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
              call unlink_blob(this, buffer_head, prev, next)
              call insert_blob(this, free_head)
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
        call dynamic_update(Time, Thickness, Ext_mode, L_system, Dens,     &
                            T_prog(:), tend_blob, blob_source, press_grad, &
                            u, w, model_rho, buffer_head, bottom,          &
                            mass_in, mass_out, EL_diag(:), ngrnd, nsfc, ndetrn)

        call mpp_set_current_pelist()

        ! Send, receive and unpack buffers for free blobs
        call send_buffer( Ebuffer_out) 
        call send_buffer( Wbuffer_out)
        call send_buffer( Sbuffer_out)
        call send_buffer( Nbuffer_out)
        call send_buffer(NEbuffer_out)
        call send_buffer(NWbuffer_out)
        call send_buffer(SEbuffer_out)
        call send_buffer(SWbuffer_out)

        call receive_buffer( Wbuffer_in) 
        call unpackbuffer(Time,  Wbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer( Ebuffer_in)
        call unpackbuffer(Time,  Ebuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer( Nbuffer_in)
        call unpackbuffer(Time,  Nbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer( Sbuffer_in)
        call unpackbuffer(Time,  Sbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer(SWbuffer_in)
        call unpackbuffer(Time, SWbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer(SEbuffer_in)
        call unpackbuffer(Time, SEbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer(NWbuffer_in)
        call unpackbuffer(Time, NWbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
        call receive_buffer(NEbuffer_in) 
        call unpackbuffer(Time, NEbuffer_in, buffer_head, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)

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
                 call unlink_blob(this, buffer_head, prev, next)
                 call insert_blob(this, free_head)
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

  ! Pack and send buffers for free blobs that have become bottom blobs

  ! Clear buffers for sending and receiving bottom blobs
  call clear_buffer(Ebuffer_out);  call clear_buffer( Ebuffer_in)
  call clear_buffer(Wbuffer_out);  call clear_buffer( Wbuffer_in)
  call clear_buffer(Nbuffer_out);  call clear_buffer( Nbuffer_in)
  call clear_buffer(Sbuffer_out);  call clear_buffer( Sbuffer_in)
  call clear_buffer(NEbuffer_out); call clear_buffer(NEbuffer_in)
  call clear_buffer(NWbuffer_out); call clear_buffer(NWbuffer_in)
  call clear_buffer(SEbuffer_out); call clear_buffer(SEbuffer_in)
  call clear_buffer(SWbuffer_out); call clear_buffer(SWbuffer_in)

  if (associated(bottom)) then
     this=>bottom
     blobcycle: do

        i = this%i
        j = this%j
        k = this%k

        if     (i>iec .and. j>jec) then
           call packbottombuffer(NEbuffer_out)
        elseif (i<isc .and. j>jec) then
           call packbottombuffer(NWbuffer_out)
        elseif (i>iec .and. j<jsc) then
           call packbottombuffer(SEbuffer_out)
        elseif (i<isc .and. j<jsc) then
           call packbottombuffer(SWbuffer_out)
        elseif (i>iec)             then
           call packbottombuffer(Ebuffer_out)
        elseif (i<isc)             then
           call packbottombuffer(Wbuffer_out)
        elseif (j>jec)             then
           call packbottombuffer(Nbuffer_out)
        elseif (j<jsc)             then
           call packbottombuffer(Sbuffer_out)
        endif
        this=>this%next
        if(.not.associated(this)) exit blobcycle
     enddo blobcycle
  endif

  ! Send buffers for free blobs that have become bottom blobs
  call send_buffer( Ebuffer_out)
  call send_buffer( Wbuffer_out)
  call send_buffer( Sbuffer_out)
  call send_buffer( Nbuffer_out)
  call send_buffer(NEbuffer_out)
  call send_buffer(NWbuffer_out)
  call send_buffer(SEbuffer_out)
  call send_buffer(SWbuffer_out)

  ! Receive and unpack free blobs that have become bottom blobs
  if (bitwise_reproduction) then
     call receive_buffer( Ebuffer_in); call unpackbuffer(Time,  Ebuffer_in, bottom, Dens)
     call receive_buffer( Wbuffer_in); call unpackbuffer(Time,  Wbuffer_in, bottom, Dens)
     call receive_buffer( Nbuffer_in); call unpackbuffer(Time,  Nbuffer_in, bottom, Dens)
     call receive_buffer( Sbuffer_in); call unpackbuffer(Time,  Sbuffer_in, bottom, Dens)
     call receive_buffer(NEbuffer_in); call unpackbuffer(Time, NEbuffer_in, bottom, Dens)
     call receive_buffer(NWbuffer_in); call unpackbuffer(Time, NWbuffer_in, bottom, Dens)
     call receive_buffer(SEbuffer_in); call unpackbuffer(Time, SEbuffer_in, bottom, Dens)
     call receive_buffer(SWbuffer_in); call unpackbuffer(Time, SWbuffer_in, bottom, Dens)
  else
     call receive_buffer( Ebuffer_in); call unpackbuffer(Time,  Ebuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer( Wbuffer_in); call unpackbuffer(Time,  Wbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer( Nbuffer_in); call unpackbuffer(Time,  Nbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer( Sbuffer_in); call unpackbuffer(Time,  Sbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(NEbuffer_in); call unpackbuffer(Time, NEbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(NWbuffer_in); call unpackbuffer(Time, NWbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(SEbuffer_in); call unpackbuffer(Time, SEbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
     call receive_buffer(SWbuffer_in); call unpackbuffer(Time, SWbuffer_in, bottom, Dens, Ext_mode, L_system, blob_source, tend_blob, mass_in, EL_diag)
  endif
  call mpp_sync_self()

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that packs the buffer for free blobs that!
! have interacted with topography to become bottom blobs and have      !
! crossed a processor boundary.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine packbottombuffer(buffer)
    type(blob_buffer_type), pointer :: buffer
    real, dimension(num_prog_tracers) :: dtracer
    real, dimension(0:num_prog_tracers) :: entrainment, detrainment
    integer :: s

    dtracer(:) = 0.0
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
  end subroutine packbottombuffer

end subroutine blob_dynamic_free_update
! </SUBROUTINE>  NAME="blob_dynamic_free_update"


!######################################################################
! <SUBROUTINE NAME="dynamic_update">
!
! <DESCRIPTION>
! This routine contains the RK scheme used to integrate the position 
! and velocity of blobs.  It also does many checks for (and 
! subsequently handles) things like grounding of blobs, blobs going to
! different PEs, blobs that interact with topography, blobs that 
! detrain to less than small_mass and blobs going outside the compute 
! domain.  
!
! It also does the interpolation of E system variables to a blob.
! </DESCRIPTION>
!
subroutine dynamic_update(Time, Thickness, Ext_mode, L_system, Dens, T_prog,   &
                          tend_blob, blob_source, press_grad, u, w, model_rho, & 
                          head, bottom_head, mass_in, mass_out, EL_diag,       &
                          ngrnd, nsfc, ndetrn)
  type(ocean_time_type),                               intent(in)    :: Time
  type(ocean_thickness_type),                          intent(inout) :: Thickness
  type(ocean_external_mode_type),                      intent(inout) :: Ext_mode
  type(ocean_lagrangian_type),                         intent(inout) :: L_system
  type(ocean_density_type),                            intent(in)    :: Dens
  type(ocean_prog_tracer_type),                        intent(in)    :: T_prog(:)
  real,dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  real,dimension(isd:ied,jsd:jed),                     intent(inout) :: blob_source
  real,dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),        intent(in)    :: press_grad
  real,dimension(isbd:iebd,jsbd:jebd,1:nk,1:2),        intent(in)    :: u
  real,dimension(isbd:iebd,jsbd:jebd,0:nk),            intent(in)    :: w
  real,dimension(isbd:iebd,jsbd:jebd,1:nk),            intent(in)    :: model_rho
  type(ocean_blob_type),                               pointer       :: head
  type(ocean_blob_type),                               pointer       :: bottom_head
  real,dimension(isd:ied,jsd:jed,1:nk),                intent(inout) :: mass_in
  real,dimension(isd:ied,jsd:jed,1:nk),                intent(inout) :: mass_out
  type(blob_diag_type),                                intent(inout) :: EL_diag(0:)
  integer,                                             intent(inout) :: ngrnd
  integer,                                             intent(inout) :: nsfc
  integer,                                             intent(inout) :: ndetrn
  
  type(ocean_blob_type), pointer :: prev, this, next
  real,    dimension(1:6,0:ns) :: V
  real,    dimension(num_prog_tracers) :: tracer, field
  real,    dimension(0:num_prog_tracers) :: entrainment, detrainment
  integer, dimension(3)   :: old_ijk
  real,    dimension(3)   :: old_lld, vel
  real,    dimension(2)   :: h, old_h
  real,    dimension(6)   :: Xn, update, Xnp1, Xnp1_hat
  real,    dimension(9)   :: tdsq_r
  real,    dimension(4)   :: udsq_r
  real,    dimension(0:2) :: px, py, uE, vE, wE, rhoE
  logical :: go_e, go_ne, go_n, go_nw, go_w, go_sw, go_s, go_se
  integer :: ii, iit, jjt, iiu, jju, dk, kdk
  integer :: i,j,k,tau,mm,m,n,r,s
  integer :: n_frac_steps
  integer :: total_blobs, leaving_pe, tfer_bottom, detrn_zero, move_lateral, ngrounded
  integer :: stdoutunit
  logical :: accept_step, reached_end, off(3), grounded
  logical :: change_pe, advance_blob, special
  real :: rhoL, rhoLr, rho, rhor, mass, volume
  real :: f, fs
  real :: lon, lat, geodepth
  real :: tstep, old_tstep, blob_time
  real :: trunc, hstar, tstar
  real :: dmdt_det, dmdt_ent
  real :: ubigd, tbigd, dzwtT, dz(2), vcoeff(2)

  tau = Time%tau

  ! Set debugging counters to zero
  total_blobs  = 0
  leaving_pe   = 0
  tfer_bottom  = 0
  detrn_zero   = 0
  move_lateral = 0
  ngrounded    = 0

  if(associated(head)) then
     this=>head

     blob_cycle: do
        ! time relative to the beginning of the time step
        blob_time = this%blob_time
        tstep     = this%step

        this%dmass      = 0.0
        this%dtracer(:) = 0.0
        if(bitwise_reproduction) then
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


        ! We dont want to alter the blob properties until 
        ! we have completed the full partial step.  So, we save
        ! the variables into local variables in case we need
        ! to redo something.

        ! These variables get updated from one RK sub-step to
        ! the next, so, we need to save copies of them in case
        ! the step is rejected.
        i        = this%i
        j        = this%j
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
        n_frac_steps = 0

        ! Set flags to default values
        reached_end  = .false.
        change_pe    = .false.
        advance_blob = .true.
        off(:)       = .false.
        grounded     = .false.

        ! Calculate the horizontal interpolation coefficients
        call interp_tcoeff(i,j,h(:),lon,lat,tdsq_r(:))
        call interp_ucoeff(i,j,h(:),lon,lat,udsq_r(:))

        call mpp_clock_begin(id_clock_part_cycle)
        partialstep: do while (.not. reached_end)
           accept_step = .false.
           
           ! Interpolate the E system variables to the blobs
           ! Note, we treat the stretching function separately
           ! Also: the stretching function does not require vertical 
           ! interpolation

           ! Figure out some stuff for the vertical interpolation
           ! for the T and U grids.  The W grid is handled differently.
           if (k==1 .and. geodepth<Thickness%geodepth_zt(i,j,k)) then
              special = .true.
              dzwtT   = Thickness%dzwtT(i,j,k)
              dk      = +1
           elseif(k==Grd%kmt(i,j) .and. geodepth>Thickness%geodepth_zt(i,j,k)) then
              special = .true.
              dzwtT   = Thickness%dzwtT(i,j,k-1)
              dk      = -1
           else
              special = .false.
              if(geodepth<Thickness%geodepth_zt(i,j,k)) then
                 dk = -1
                 dzwtT = Thickness%dzwtT(i,j,k-1)
              else
                 dk = +1
                 dzwtT = Thickness%dzwtT(i,j,k)
              endif
           endif
           kdk = k+dk

           ! Now do the horizontal interpolation at two levels, k and k+dk
           ubigd = 0.0
           px(:)=0.0; py(:)=0.0 
           uE(:)=0.0; vE(:)=0.0;
           do mm=1,Info%uidx(0,i,j)
              iiu=i+Info%iu(Info%uidx(mm,i,j))
              jju=j+Info%ju(Info%uidx(mm,i,j))
              ! Land points are treated as zero for pressure gradient
              ! and velocity
              px(1) = px(1) + umask(iiu,jju,k  )*press_grad(iiu,jju,k  ,1)*udsq_r(mm)
              px(2) = px(2) + umask(iiu,jju,kdk)*press_grad(iiu,jju,kdk,1)*udsq_r(mm)
              py(1) = py(1) + umask(iiu,jju,k  )*press_grad(iiu,jju,k  ,2)*udsq_r(mm)
              py(2) = py(2) + umask(iiu,jju,kdk)*press_grad(iiu,jju,kdk,2)*udsq_r(mm)
              uE(1) = uE(1) + umask(iiu,jju,k  )*u(iiu,jju,k  ,1)*udsq_r(mm)
              uE(2) = uE(2) + umask(iiu,jju,kdk)*u(iiu,jju,kdk,1)*udsq_r(mm)
              vE(1) = vE(1) + umask(iiu,jju,k  )*u(iiu,jju,k  ,2)*udsq_r(mm)
              vE(2) = vE(2) + umask(iiu,jju,kdk)*u(iiu,jju,kdk,2)*udsq_r(mm)
              ubigd = ubigd + udsq_r(mm)
           enddo
           ubigd = 1.0/ubigd
           
           px(1) = px(1)*ubigd 
           px(2) = px(2)*ubigd 
           py(1) = py(1)*ubigd 
           py(2) = py(2)*ubigd 
           uE(1) = uE(1)*ubigd 
           uE(2) = uE(2)*ubigd 
           vE(1) = vE(1)*ubigd 
           vE(2) = vE(2)*ubigd 

           rhoE(:) = 0.0
           tbigd   = 0.0
           do mm=1,Info%tidx(0,i,j)
              iit=i+Info%it(Info%tidx(mm,i,j))
              jjt=j+Info%jt(Info%tidx(mm,i,j))
              ! We ignore density for land points
              if (tmask(iit,jjt,k  ) /= 0.0) then
                 rhoE(1) = rhoE(1) + model_rho(iit,jjt,k  )*tdsq_r(mm)
                 tbigd = tbigd + tdsq_r(mm)
              endif
           enddo

           rhoE(1) = rhoE(1)/tbigd

           tbigd   = 0.0
           do mm=1,Info%tidx(0,i,j)
              iit=i+Info%it(Info%tidx(mm,i,j))
              jjt=j+Info%jt(Info%tidx(mm,i,j))
              ! We ignore density for land points
              if (tmask(iit,jjt,kdk) /= 0.0) then
                 rhoE(2) = rhoE(2) + model_rho(iit,jjt,kdk)*tdsq_r(mm)
                 tbigd = tbigd + tdsq_r(mm)
              endif
           enddo

           rhoE(2) = rhoE(2)/tbigd

           ! Vertically interpolate the horizontally interpolated values
           if (special) then
              ! The blob is either between the top tracer point and the surface
              ! or the lowest wet tracer point and the bottom
              dz(1)     = abs(geodepth - Thickness%geodepth_zt(i,j,k))
              vcoeff(1) = (dzwtT+dz(1))/dzwtT

              px(0)   = px(2)   + ( px(1)-px(2)     )*vcoeff(1)
              py(0)   = py(2)   + ( py(1)-py(2)     )*vcoeff(1)
              uE(0)   = uE(2)   + ( uE(1)-uE(2)     )*vcoeff(1)
              vE(0)   = vE(2)   + ( vE(1)-vE(2)     )*vcoeff(1)
              rhoE(0) = rhoE(2) + ( rhoE(1)-rhoE(2) )*vcoeff(1)
           else
              ! The blob lies vertically between two tracer points
              dz(1) = abs(geodepth - Thickness%geodepth_zt(i,j,k  ))
              dz(2) = abs(geodepth - Thickness%geodepth_zt(i,j,kdk))
              
              vcoeff(1) = dz(1)/dzwtT
              vcoeff(2) = dz(2)/dzwtT
              
              px(0)   = px(1)*vcoeff(2)   + px(2)*vcoeff(1)
              py(0)   = py(1)*vcoeff(2)   + py(2)*vcoeff(1)
              uE(0)   = uE(1)*vcoeff(2)   + uE(2)*vcoeff(1)
              vE(0)   = vE(1)*vcoeff(2)   + vE(2)*vcoeff(1)
              rhoE(0) = rhoE(1)*vcoeff(2) + rhoE(2)*vcoeff(1)
           endif

           ! Now handle the W grid
           kdk = k-1
           wE(:) = 0.0
           tbigd = 0.0
           do mm=1,Info%tidx(0,i,j)
              iit=i+Info%it(Info%tidx(mm,i,j))
              jjt=j+Info%jt(Info%tidx(mm,i,j))
              ! Vertical velocity is treated as zero for land points
              ! Note, this is actually the velocity of the fluid crossing
              ! coordinate surfaces, NOT the true velocity.  Thus, we are assuming
              ! that the vertical velocity of coordinate surfaces is small compared to the
              ! true vertical velocity.
              wE(1) = wE(1) + tmask(iit,jjt,k  )*w(iit,jjt,k  )*tdsq_r(mm)
              wE(2) = wE(2) + tmask(iit,jjt,kdk)*w(iit,jjt,kdk)*tdsq_r(mm)
              tbigd = tbigd + tdsq_r(mm)
           enddo
           tbigd = 1.0/tbigd
           wE(1) = wE(1)*tbigd
           wE(2) = wE(2)*tbigd

           dz(1) = Thickness%geodepth_zwt(i,j,k) - geodepth
           if (kdk==0) then
              dz(2) = geodepth - Ext_mode%eta_t(i,j,tau)
           else
              dz(2) = geodepth - Thickness%geodepth_zwt(i,j,kdk)
           endif
           wE(0) = (wE(1)*dz(2) + wE(2)*dz(1))/Thickness%dztT(i,j,k,tau)

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

              ! rotation
              f  = two_omega*sin(deg_to_rad*lat)
              fs = two_omega*cos(deg_to_rad*lat)
                 
              ! Cycle through the RK sub-steps
              call mpp_clock_begin(id_clock_rk)
              do s=0,ns
                 ! dx/dt and d2x/dt2
                 V(1,s) =   update(2)
                 V(2,s) = - this%drag*update(2) + f*update(4) - fs*update(6) - px(0)*rhor + this%drag*uE(0)

                 ! dy/dt and d2y/dt2
                 V(3,s) =   update(4)
                 V(4,s) = - f*update(2) - this%drag*update(4) - py(0)*rhor + this%drag*vE(0)
                 
                 ! dz/dt and d2z/dt2
                 V(5,s) =   update(6)
                 V(6,s) =   fs*update(2) - this%drag*update(6) - grav*rhor*(rhoL - rhoE(0)) + this%drag*wE(0)
                 
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

              enddo !s
              call mpp_clock_end(id_clock_rk)
              
              ! Estimate the error using the difference in the high and low order schemes.
              ! We only use the position (not velocity) in the error estimate.
              ii    = maxloc(abs((Xnp1(:) - Xnp1_hat(:))/(Xnp1(:)+epsln)),1)
              trunc = abs(Xnp1(ii) - Xnp1_hat(ii) + epsln)
              tstar = rel_error*abs(Xnp1(ii))
              hstar = tstep * (safety_factor*tstar/trunc)**(orderp1_r)

              ! Save the tstep so that we can update blob_time
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
                            '==>Warning in ocean_blob_dynamic_free_mod (dynamic_update): It looks like a blob is becoming unstable'//&
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
                 h(:) = h(:)/tbigd

                 n_frac_steps = n_frac_steps+1

                 ! Entrainment and detrainment calculations are based on the properties of the blob
                 ! at the beginning of this blob step.
                 ! Calculate the detrainment
                 if (vert_coordinate_class==DEPTH_BASED) then
                    !det_factor = Gamma * (rho0*36pi)**1/3
                    dmdt_det = -det_factor * (mass**twothirds) / abs(rhoL - rhoE(0) + epsln)
                 else !PRESSURE_BASED
                    !det_factor = Gamma * (36pi)**1/3
                    dmdt_det = -det_factor * (mass**twothirds) * (rho**onethird) / abs(rhoL - rhoE(0) + epsln)
                 endif

                 ! Make sure we do not exceed the maximum detrainment
                 dmdt_det = sign(min(abs(dmdt_det),abs(max_detrainment)),dmdt_det)

                 ! Entrainment
                 dmdt_ent = 0.0

                 detrainment(0) = dmdt_det*old_tstep
                 entrainment(0) = dmdt_ent*old_tstep

                 do n=1,num_prog_tracers
                    ! Detrainment
                    detrainment(n) = detrainment(0)*field(n)
                    ! Entrainment
                    entrainment(n) = entrainment(0)*T_prog(n)%field(i,j,k,tau)
                 enddo

                 ! Check and see if the blob will detrain to zero mass during this step.
                 ! If it does, just dump all its properties in the original cell.
                 if ((entrainment(0)+detrainment(0)+mass) < small_mass) then
                    entrainment(0) = 0.0
                    detrainment(0) = -mass
                    do n=1,num_prog_tracers
                       entrainment(n) = 0.0
                       detrainment(n) = -tracer(n)
                    enddo
                    ! Update some counters
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

                 ! Check to see if the blob has grounded on a land column
                 if (Grd%kmt(i,j)==0) then
                    ! If the blob has grounded, we put it back in its previous cell
                    ! It will have all its properties returned to the E system later
                    ! Here, grounded means a blob has moved laterally into a zero 
                    ! depth (i.e. land) column, NOT that it has interacted with topography.
                    off(1:2)  = .false.
                    grounded  = .true.
                    ngrounded = ngrounded + 1
                    i = old_ijk(1)
                    j = old_ijk(2)
                    lon = old_lld(1)
                    lat = old_lld(2)

                 endif

                 if (any(off(1:2))) change_pe = .true.
                 
                 ! Check if we have changed vertical cells and whether we have interacted with
                 ! the lower or upper boundary
                 if (.not. grounded) then
                    call check_kcell(Time, Ext_mode, Thickness,-Xnp1(5),Xnp1(6),i,j,k,off(3))
                 endif

                 if(change_pe .and. .not. bitwise_reproduction) then
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
                       this%step = max(tstep, minstep)
                    elseif ((this%blob_time + 1.1*tstep) > dtime) then
                       this%step = dtime - this%blob_time
                    else
                       this%step = max(tstep,minstep)
                    endif
                    this%mass      = mass
                    
                    this%i = i
                    this%j = j
                    this%k = k

                    ! If the blob interacts with topography in the new cell, we need to put it in the bottom
                    ! blob list.  Later, it will be packed into a buffer and sent to a neighbouring PE.
                    if (off(3) .and. k==Grd%kmt(i,j)) then
                       this%age = this%age+dtime_yr
                       call unlink_blob(this, head, prev, next)
                       call insert_blob(this, bottom_head)
                       tfer_bottom  = tfer_bottom + 1
                       this%nfrac_steps = n_frac_steps
                    else
                       if     (i>iec .and. j>jec) then
                          call packbuffer(this, NEbuffer_out, entrainment(:), detrainment(:))
                       elseif (i<isc .and. j>jec) then
                          call packbuffer(this, NWbuffer_out, entrainment(:), detrainment(:))
                       elseif (i>iec .and. j<jsc) then
                          call packbuffer(this, SEbuffer_out, entrainment(:), detrainment(:))
                       elseif (i<isc .and. j<jsc) then
                          call packbuffer(this, SWbuffer_out, entrainment(:), detrainment(:))
                       elseif (i>iec)             then
                          call packbuffer(this,  Ebuffer_out,  entrainment(:), detrainment(:))
                       elseif (i<isc)             then
                          call packbuffer(this,  Wbuffer_out,  entrainment(:), detrainment(:))
                       elseif (j>jec)             then
                          call packbuffer(this,  Nbuffer_out,  entrainment(:), detrainment(:))
                       elseif (j<jsc)             then
                          call packbuffer(this,  Sbuffer_out, entrainment(:), detrainment(:))
                       endif

                       call unlink_blob(this, head, prev, next)
                       call free_blob_memory(this)
                    endif

                    move_lateral = move_lateral + 1
                    leaving_pe   = leaving_pe   + 1

                    advance_blob=.false.

                    exit partialstep
                 endif!change_pe .and. .not. bitwise_reproduction

                
                 ! A little about the following order of operations:
                 ! Detrained mass and tracer is "given" to the cell that the blob is
                 ! entering (i.e. the cell that the blob resides in at the end of
                 ! a sub-cycle).  The mass used for calculating convergence/divergence 
                 ! must therefore be the mass of the blob before detrainment occurs.
                 if (bitwise_reproduction) then

                    ! If we detect that a blob changes cells, save the mass into the new
                    ! water column and and out of the old water column.
                    ! Note, the convention is that both are positive and we take care of
                    ! the sign when calculating the convergence.
                    if (old_ijk(1)/=i .or. old_ijk(2)/=j .or. old_ijk(3)/=k) then
                       this%mass_out(n_frac_steps-1) = mass
                       this%mass_in( n_frac_steps  ) = mass
                    endif

                    ! Save the indices to history too (for divergence and ent/detrainment)
                    this%di(n_frac_steps) = i
                    this%dj(n_frac_steps) = j
                    this%dk(n_frac_steps) = k
                    
                    ! save mass
                    mass = mass + entrainment(0) + detrainment(0)
                    this%entrainment(n_frac_steps,0) = entrainment(0)
                    this%detrainment(n_frac_steps,0) = detrainment(0)
                    
                    ! save tracer
                    do n=1,num_prog_tracers
                       tracer(n) = tracer(n) + entrainment(n) + detrainment(n)
                       this%entrainment(n_frac_steps,n) = entrainment(n)
                       this%detrainment(n_frac_steps,n) = detrainment(n)
                    enddo
                    
                 else!not bitwise_reproduction
                    ! If we detect that a blob changes cells, save the mass into the new
                    ! water column in to and out of the old water column.
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

                    mass = mass + entrainment(0) + detrainment(0)
                    blob_source(i,j) = blob_source(i,j) - entrainment(0) - detrainment(0)
                    
                    ! save tracer
                    do n=1,num_prog_tracers
                       tracer(n) = tracer(n) + entrainment(n) + detrainment(n)
                       tend_blob(n,i,j,k) = tend_blob(n,i,j,k) - entrainment(n) - detrainment(n)
                    enddo

                    do n=0,num_prog_tracers 
                       EL_diag(n)%detrainment(i,j,k) = EL_diag(n)%detrainment(i,j,k) - detrainment(n) 
                       EL_diag(n)%entrainment(i,j,k) = EL_diag(n)%entrainment(i,j,k) + entrainment(n) 
                    enddo 

                 endif!bitwise_reproduction

                 ! save other variables
                 geodepth = -Xnp1(5)
                 old_lld  =  (/lon, lat, geodepth/)
                 old_ijk  =  (/i,   j,   k       /)
                 old_h    =  h
                 vel(1) = Xnp1(2)
                 vel(2) = Xnp1(4)
                 vel(3) = Xnp1(6)

                 !Avoid crazy numbers if the blob has detrained to be less than small mass
                 if(mass>small_mass) then
                    do n=1,num_prog_tracers
                       field(n) = tracer(n)/mass
                    enddo

                    rhoL  = density(field(index_salt), field(index_temp), Dens%pressure_at_depth(i,j,k))
                    rhoLr = 1./rhoL
                    if (vert_coordinate_class == PRESSURE_BASED) then
                       rho  = rhoL
                       rhor = rhoLr
                       !else DEPTH_BASED; rho=rho0; rhor=rho0r
                    endif
                 endif
                 volume = mass*rhor

                 !Now we handle the blobs that have penetrated the surface or bottom boundaries
                 if (off(3) .or. grounded) then
                    ! Update a counter
                    if (this%i /= i .or. this%j /= j) move_lateral = move_lateral + 1
                    if (.not. grounded .and. k==1) nsfc = nsfc+1
                    
                    ! Update the blob variables
                    this%step = max(tstep, minstep)
                    this%i    = i
                    this%j    = j
                    this%k    = k
                    this%lon  = lon
                    this%lat  = lat
                    this%geodepth = geodepth
                    this%v(:) = vel(:)
                    this%h1   = h(1)
                    this%h2   = h(2)
                    
                    this%mass     = mass
                    this%volume   = volume
                    this%density  = rhoL
                    this%densityr = rhoLr
                    this%gprime   = grav*(rhoE(0)-rhoL)/rhoE(0) !diagnostic
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
                    
                    if (k>=Grd%kmt(i,j) .and. .not. grounded) then !the blob has penetrated the bottom boundary
                       ! The geodepth and velocity will be adjusted in transfer_free_to_bottom
                       call unlink_blob(this, head, prev, next)
                       call insert_blob(this, bottom_head)
                       advance_blob=.false.
                       tfer_bottom  = tfer_bottom + 1
                       this%nfrac_steps = n_frac_steps
                       exit partialstep
                    else !the blob has penetrated the free surface, or, the blob has grounded
                       n_frac_steps = n_frac_steps+1

                       ! Diagnostics
                       EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + mass
                       do n=1,num_prog_tracers
                          EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + tracer(n)
                       enddo

                       if (bitwise_reproduction) then

                          this%nfrac_steps      = n_frac_steps
                          this%di(n_frac_steps) = this%i
                          this%dj(n_frac_steps) = this%j
                          this%dk(n_frac_steps) = this%k
                          this%dmass = -mass
                          this%mass  = 0.0
                          mass       = 0.0
                          do n=1,num_prog_tracers
                             this%dtracer(n) = -tracer(n)
                             this%tracer(n)  = 0.0
                             tracer(n)       = 0.0
                          enddo

                          ! If a blob changes PEs and it interacts with the surface boundary
                          ! we still need it to be packed into a buffer and sent to
                          ! the neighbouring PE for the history.  So, we fool the 
                          ! algorithm into thinking that the blob has completed its
                          ! full time step so that it will get packed into a buffer if necessary.
                          ! We do not need to do this procedure for the blob interacting with 
                          ! the bottom boundary, because it is handled separately in its own 
                          ! linked list.
                          old_tstep = 0.0
                          blob_time = dtime
                          exit trystep

                       else!not bitwise_reproduction
                          this%nfrac_steps = n_frac_steps

                          blob_source(this%i,this%j) = blob_source(this%i,this%j) + mass
                          this%mass = 0.0
                          mass      = 0.0
                          do n=1,num_prog_tracers
                             tend_blob(n,this%i,this%j,this%k) = tend_blob(n,this%i,this%j,this%k) + tracer(n)
                             this%tracer(n) = 0.0
                             tracer(n)      = 0.0
                          enddo

                          ! Diagnostics 
                          EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + mass 
                          do n=1,num_prog_tracers 
                             EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + tracer(n) 
                          enddo

                          ! For bitwise_reproduction=.false. we do not need to worry about
                          ! histories being stored in buffers, so, we just exit partialstep
                          !exit partialstep
                          old_tstep = 0.0
                          blob_time = dtime
                          exit trystep

                       endif!bitwise_reproduction
                    endif!k==kmt and not grounded
                 endif!off(3) or grounded

                 this%nfrac_steps = n_frac_steps
              endif !accept_step

              m=m+1
           enddo trystep
           
           ! The step has been accepted.  We now decide whether we need
           ! to conduct any more steps, and if we do, whether we need to 
           ! adjust tstep to ensure the final step coincides with the Eulerian
           ! model time step.
           blob_time      = blob_time + old_tstep
           this%blob_time = blob_time

           if (blob_time == dtime) then
              ! no more steps required.  Save all the new variables to the blob
              reached_end = .true.

              ! Update a counter
              if (this%i /= i .or. this%j /= j) move_lateral = move_lateral + 1

              ! Update the blob variables
              this%step = max(tstep,minstep)
              this%i    = i
              this%j    = j
              this%k    = k
              this%lon  = lon
              this%lat  = lat
              this%geodepth = geodepth
              this%v(:) = vel(:)
              this%h1   = h(1)
              this%h2   = h(2)

              this%mass     = mass
              this%volume   = volume
              this%density  = rhoL
              this%densityr = rhoLr
              this%gprime   = grav*(rhoE(0)-rhoL)/rhoE(0) !diagnostic
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

                 if(change_pe) then
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

                    ! We leave the blobs history in the linked list on this processor
                    ! but, we remove its properties.  It will be deleted next call
                    ! to blob_delete, but, we need its history to be processed 
                    ! before it is deleted.  The history is processed only after
                    ! all blobs have been stepped.
                    ! A blob may go from one PE and then back to the same PE within an
                    ! E system time step, e.g.
                    ! ------------------
                    ! | PE=2  +a  PE=3 |
                    ! |      /|        |
                    ! |     b |        |
                    ! |      \|        |
                    ! |       +c       |
                    ! ------+-----------
                    ! |       |        |
                    ! |       |        |
                    ! |       |        |
                    ! |       |        |
                    ! |  PE=0 |   PE=1 |
                    ! ------------------
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
              
           elseif( (blob_time + 1.1*tstep) > dtime ) then
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
           endif !blobtime==dtime

        enddo partialstep
        call mpp_clock_end(id_clock_part_cycle)
        total_blobs = total_blobs + 1

        if (advance_blob) then
           this=>this%next
        else
           this=>next
        endif
        if (.not. associated(this)) exit blob_cycle
        
     end do blob_cycle
  endif !blob associated?

  ngrnd  = ngrnd  + ngrounded
  ndetrn = ndetrn + detrn_zero

  if (debug_this_module .and. bitwise_reproduction) then
     stdoutunit = stdout()
     write (stdoutunit, '(/,a)') 'Dynamic Free Blob Statistics'
     call write_timestamp(Time%model_time)
     call mpp_sum(total_blobs)
     write (stdoutunit, *) 'Total Free Dynamic Blobs             =', total_blobs
     call mpp_sum(leaving_pe)
     write (stdoutunit, *) 'Free Dynamic Blobs Changing PEs      =', leaving_pe
     call mpp_sum(tfer_bottom)
     write (stdoutunit, *) 'Free Dynamic Blob Transfer to bottom =', tfer_bottom
     call mpp_sum(detrn_zero)
     write (stdoutunit, *) 'Free Blobs detrained to zero mass    =', detrn_zero
     call mpp_sum(move_lateral)
     write (stdoutunit, *) 'Free Blobs changing water columns    =', move_lateral
     call mpp_sum(ngrounded)
     write (stdoutunit, *) 'Grounded Free Blobs                  =', ngrounded
  endif
     
end subroutine dynamic_update
! </SUBROUTINE>  NAME="dynamic_upate"

!######################################################################
! <SUBROUTINE NAME="transfer_bottom_to_free">
!
! <DESCRIPTION>
! Takes bottom blobs that have separated from the bottom boundary and
! turns it into a free blob.
! </DESCRIPTION>
!
subroutine transfer_bottom_to_free(Time, Thickness, T_prog, Dens, &
                                   new_free, head, use_bottom,    &
                                   EL_diag)
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_blob_type),        pointer       :: new_free
  type(ocean_blob_type),        pointer       :: head
  logical,                      intent(in)    :: use_bottom
  type(blob_diag_type), dimension(0:num_prog_tracers) :: EL_diag

  type(ocean_blob_type), pointer :: this, next, prev
  integer :: i,j,k,tau,nblobs,n
  real, parameter :: small=1e-3
  integer :: stdoutunit

  tau = Time%tau
  nblobs = 0

  if (associated(new_free)) then
     this => new_free
     if (use_this_module) then
        blobcycle: do
           i = this%i
           j = this%j
           k = this%k
           
           ! Ensure that the blob is above the Eulerian topography
           ! so that the free blob module does not think that we have
           ! penetrated rock.
           if (vert_coordinate_class==DEPTH_BASED) then
              this%st       = -Thickness%depth_swt(i,j,k)    + small
              this%geodepth =  Thickness%geodepth_zwt(i,j,k) - small*Thickness%dzt_dst(i,j,k)
           else!PRESSURE_BASED
              this%st       = Thickness%depth_swt(i,j,k)    - small
              this%geodepth = Thickness%geodepth_zwt(i,j,k) + small*Thickness%dzt_dst(i,j,k)
           endif
           
           ! The drag coefficient
           this%drag = rayleigh_drag_bot

           ! Enforce a positive or zero velocity on the blob, so that any downward inertia 
           ! wont cause it to penetrate rock again straight away.
           if (this%v(3)<0) this%v(3)=0.
           
           ! Update some of the blob variables
           this%density  = density(this%field(index_salt), &
                                   this%field(index_temp), &
                                   Dens%pressure_at_depth(i,j,k))
           this%densityr = 1./this%density
           if (vert_coordinate_class == DEPTH_BASED) then
              this%volume = this%mass * rho0r 
           else
              this%volume = this%mass * this%densityr
           endif
           
           call unlink_blob(this, new_free, prev, next)
           call insert_blob(this, head)
           call reallocate_interaction_memory(this,head,total_ns)

           this=>next
           nblobs = nblobs+1
           if(.not.associated(this)) exit blobcycle
        enddo blobcycle

     else !not use_this_module
        ! If we are not using dynamic free blobs, return the 
        ! transferred blbos properties to the E system.
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
           nblobs=nblobs+1
           if(.not.associated(this)) exit blobcycle2
        enddo blobcycle2
     endif !use_this_module
     call blob_delete(Time, Thickness, T_prog(:), new_free)
  endif!associated(new_free)

  call mpp_sum(nblobs)
  if (id_bot_to_free>0) used = send_data(id_bot_to_free, real(nblobs), Time%model_time)

  if(debug_this_module .and. use_bottom) then
     stdoutunit = stdout()
     write(stdoutunit, '(/,a)') 'Bottom blobs separating from topograhy'
     write(stdoutunit, *) 'Bottom blobs transferred to free blobs = ', nblobs
  endif

end subroutine transfer_bottom_to_free

! </SUBROUTINE>  NAME="transfer_bottom_to_free"

!######################################################################
! <SUBROUTINE NAME="blob_dynamic_free_end">
!
! <DESCRIPTION>
! Clears memory to give a nice clean ending to the run.
! </DESCRIPTION>
!
subroutine blob_dynamic_free_end()

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
end subroutine blob_dynamic_free_end
! </SUBROUTINE>  NAME="blob_dynamic_free_end"


!######################################################################
! <SUBROUTINE NAME="packbuffer">
!
! <DESCRIPTION>
! Packs a buffer with all the information needed to send a blob from
! one PE to another.
! </DESCRIPTION>
!
subroutine packbuffer(blob,buffer,entrainment,detrainment)
  type(ocean_blob_type),  pointer    :: blob
  type(blob_buffer_type), pointer    :: buffer
  real,         optional, intent(in) :: entrainment(0:)
  real,         optional, intent(in) :: detrainment(0:)
  integer :: n, nb, s, npt
  integer :: stdoutunit

  stdoutunit = stdout()

  if (buffer%pe == NULL_PE) then
     write (stdoutunit, '(a)'), 'Error: Trying to send blob to a NULL_PE'
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
  buffer%real_buff(2*npt+12,nb) = blob%drag
  buffer%real_buff(2*npt+13,nb) = blob%dmass
  buffer%real_buff(2*npt+14,nb) = blob%gprime
  buffer%real_buff(2*npt+15,nb) = blob%age

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
     buffer%history_integer_buff(1,0,nb) = blob%di(0)
     buffer%history_integer_buff(2,0,nb) = blob%dj(0)
     buffer%history_integer_buff(3,0,nb) = blob%dk(0)
     
     buffer%history_real_buff(npt+3,0,nb) = blob%mass_out(0)
     if (blob%nfrac_steps>0) then
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
           buffer%history_real_buff(npt+2,s,nb) = blob%mass_in(s)
           buffer%history_real_buff(npt+3,s,nb) = blob%mass_out(s)
        enddo
     endif
  else
     do n=0,num_prog_tracers
        buffer%history_real_buff(2*n  ,1,nb) = entrainment(n)
        buffer%history_real_buff(2*n+1,1,nb) = detrainment(n)
     enddo
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
subroutine unpackbuffer(Time, buffer, head, Dens, Ext_mode, L_system, &
                        blob_source, tend_blob, mass_in, EL_diag)
  type(ocean_time_type),                    intent(in)    :: Time
  type(blob_buffer_type),                   pointer       :: buffer
  type(ocean_blob_type),                    pointer       :: head
  type(ocean_density_type),                 intent(in)    :: Dens
  type(ocean_external_mode_type), optional, intent(inout) :: Ext_mode
  type(ocean_lagrangian_type),    optional, intent(inout) :: L_system
  type(blob_diag_type),           optional, intent(inout) :: EL_diag(0:)
  real, optional, dimension(num_prog_tracers,isd:ied,jsd:jed,nk), intent(inout) :: tend_blob
  real, optional, dimension(isd:ied,jsd:jed),                     intent(inout) :: blob_source
  real, optional, dimension(isd:ied,jsd:jed,1:nk),                intent(inout) :: mass_in
  
  type(ocean_blob_type), pointer :: blob
  real, dimension(0:num_prog_tracers) :: entrainment, detrainment
  integer :: n, nb, s, npt
  integer :: i,j,k,tau

  tau = Time%tau

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
        blob%v(1)     = buffer%real_buff(2*npt+ 1,nb)
        blob%v(2)     = buffer%real_buff(2*npt+ 2,nb)
        blob%v(3)     = buffer%real_buff(2*npt+ 3,nb)
        blob%lon      = buffer%real_buff(2*npt+ 4,nb)
        blob%lat      = buffer%real_buff(2*npt+ 5,nb)
        blob%geodepth = buffer%real_buff(2*npt+ 6,nb)
        blob%step     = buffer%real_buff(2*npt+ 7,nb)
        blob%mass     = buffer%real_buff(2*npt+ 8,nb)
        blob%blob_time= buffer%real_buff(2*npt+ 9,nb)
        blob%h1       = buffer%real_buff(2*npt+10,nb)
        blob%h2       = buffer%real_buff(2*npt+11,nb)
        blob%drag     = buffer%real_buff(2*npt+12,nb)
        blob%dmass    = buffer%real_buff(2*npt+13,nb)
        blob%gprime   = buffer%real_buff(2*npt+14,nb)
        blob%age      = buffer%real_buff(2*npt+15,nb)
        
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

        ! If a blob has zero mass, we are only interested in its history arrays
        ! So, we don't need to bother with finding tracer concentration, density,
        ! or volume.
        if (isc <= i .and. i<=iec .and. jsc<=j .and. j<=jec .and. blob%mass>0) then
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
                 do n=1,num_prog_tracers
                    blob%entrainment(s,n) = buffer%history_real_buff(2*n  ,s,nb)
                    blob%detrainment(s,n) = buffer%history_real_buff(2*n+1,s,nb)
                 enddo
                 blob%mass_in(s)  = buffer%history_real_buff(npt+2,s,nb)
                 blob%mass_out(s) = buffer%history_real_buff(npt+3,s,nb)
              enddo

           endif
              
           do s=0,blob%nfrac_steps
              call check_cyclic(blob, blob%di(s), blob%dj(s), .false.)
           enddo

           call insert_blob(blob, head)

        else !not bitwise_reproduction
           ! blob_source is not always parsed, so, we do 
           ! not need to worry about it
           if(present(blob_source)) then !if blob_source is present, so is tend_blob, L_system and EL_diag
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

              ! Check if the blob being received has penetrated the surface boundary.
              ! If so, return its properties to the E system and do not add it to
              ! the list.
              if (blob%geodepth < -Ext_mode%eta_t(i,j,tau)) then
                 ! Diagnostics
                 EL_diag(0)%dstry(i,j,k) = EL_diag(0)%dstry(i,j,k) + blob%mass
                 do n=1,num_prog_tracers
                    EL_diag(n)%dstry(i,j,k) = EL_diag(n)%dstry(i,j,k) + blob%tracer(n)
                 enddo
              
                 blob_source(i,j) = blob_source(i,j) + blob%mass
                 blob%mass        = 0.0
                 do n=1,num_prog_tracers
                    tend_blob(n,i,j,k) = tend_blob(n,i,j,k) + blob%tracer(n)
                    blob%tracer(n)     = 0.0
                 enddo

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
  new_buffer%history_real_buff(:,:,1:buffer%size) = buffer%history_real_buff(:,:,:)

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
              call mpp_send(buffer%history_real_buff(:,0:n_frac_steps,nb),  (1+hist_rea_buff_size)*(n_frac_steps+1), buffer%pe)
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
  integer :: incoming, nb, n_frac_steps
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
              call mpp_recv(buffer%history_real_buff(:,0:n_frac_steps,nb),   (1+hist_rea_buff_size)*(n_frac_steps+1), buffer%pe)
           else
              call mpp_recv(buffer%history_real_buff(:,1,nb), 1+hist_rea_buff_size, buffer%pe)
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
     if (bitwise_reproduction) then
        buffer%history_integer_buff(:,:,:) = 0
        buffer%history_real_buff(:,:,:)    = 0.0
     endif
  endif
  !note, we do not clear buffer%size
end subroutine clear_buffer
! </SUBROUTINE>  NAME="clear_buffer"

end module ocean_blob_dynamic_free_mod
