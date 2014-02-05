module ocean_obc_barotrop_mod
#define COMP isc:iec,jsc:jec
  !
  ! <CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt </CONTACT>
  ! <CONTACT EMAIL="Mike.Herzfeld@csiro.au"> Mike Herzfeld </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matthew Harrison </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies </REVIEWER>
  !
  ! 
  ! <OVERVIEW>
  !   Open Boundary condition for MOM.
  ! </OVERVIEW>
  !
  ! <DESCRIPTION>
  !   This module can extrapolate data on the open lateral
  !   boundaries for MOM. Tracer and surface height 
  !   are extrapolated on the boundary by using implicit radiation 
  !   boundary conditions, velocities are calculated on the boundary
  !   from a linear equation (omitted advection equation). The 
  !   gradient of each field is supposed to be zero  between boundary
  !   points and the first points accross the boundary.
  !
  !   This scheme has been tested only with the following vertical coordinates:
  !   vertical_coordinate=='geopotential' 
  !   vertical_coordinate=='zstar'
  !   Notably, there is no OBC prescription for pressure coordinates. 
  !
  ! </DESCRIPTION>
#include <fms_platform.h>

  use constants_mod,            only: grav
  use data_override_mod,        only: data_override
  use diag_manager_mod,         only: register_diag_field, send_data, diag_axis_init
  use fms_mod,                  only: write_version_number, open_namelist_file, close_file, file_exist
  use fms_mod,                  only: stdout, stdlog, check_nml_error, FATAL, NOTE 
  use fms_io_mod,               only: register_restart_field, save_restart, restore_state
  use fms_io_mod,               only: reset_field_pointer, restart_file_type
  use mpp_io_mod,               only: mpp_open, mpp_close, MPP_MULTI, MPP_OVERWR, MPP_SINGLE
  use mpp_domains_mod,          only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,          only: domain2d, mpp_update_domains, mpp_get_compute_domains
  use mpp_domains_mod,          only: mpp_set_compute_domain, mpp_set_data_domain, mpp_set_global_domain
  use mpp_domains_mod,          only: mpp_define_domains
  use mpp_domains_mod,          only: BGRID_NE, SCALAR_PAIR

  use mpp_mod,                  only: input_nml_file, mpp_error, mpp_pe, mpp_chksum, mpp_npes, mpp_send, mpp_recv
  use mpp_mod,                  only: CLOCK_MODULE
  use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod,                  only: mpp_get_current_pelist, mpp_declare_pelist, mpp_sync_self
  use time_interp_external_mod, only: time_interp_external, init_external_field, get_external_field_size
  use time_manager_mod,         only: time_type
  use tracer_manager_mod,       only: get_tracer_names, get_tracer_indices, get_number_tracers
  use ocean_util_mod,           only: write_timestamp, write_chksum_2d

  use ocean_domains_mod,        only: get_local_indices, get_domain_offsets
  use ocean_parameters_mod,     only: missing_value, rho0, GEOPOTENTIAL
  use ocean_types_mod,          only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
  use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_adv_vel_type, ocean_options_type
  use ocean_types_mod,          only: ocean_external_mode_type, obc_flux
  use ocean_types_mod,          only: ocean_time_type, ocean_time_steps_type

implicit none
  private

interface
   integer function tm_scale_to_secs(buf, sec)
     character(len=16) :: buf
     real              :: sec
   end function tm_scale_to_secs
end interface


  public :: ocean_obc_barotrop_init, ocean_obc_end
  public :: ocean_obc_barotropic
  public :: ocean_obc_update_boundary, ocean_obc_prepare
  public :: ocean_obc_zero_boundary
  public :: ocean_obc_check_for_update
  public :: ocean_obc_adjust_divud
  public :: ocean_obc_damp_newton, ocean_obc_ud
  public :: mpp_update_domains_obc
  interface mpp_update_domains_obc
    module procedure mpp_update_domains_obc
  end interface mpp_update_domains_obc
!  public :: ocean_obc_restart

  !--- some module variables -------------------------------------------
  integer, parameter :: max_obc = 4             ! maximum number of open boundaries (increase if want more)
#ifdef USE_OCEAN_BGC
  integer, parameter :: max_prog_tracers=40     ! maximum number of prognostics tracers (increase if want more)
                                                ! it can be increased if needed.
#else                                                
  integer, parameter :: max_prog_tracers=10     ! maximum number of prognostics tracers (increase if want more)
#endif
  ! Boundary location flags
  integer, parameter :: WEST  = 1               ! western boundary
  integer, parameter :: EAST  = 2               ! eastern boundary
  integer, parameter :: SOUTH = 3               ! southern boundary
  integer, parameter :: NORTH = 4               ! northern boundary

  ! Boundary condition flags
  integer, parameter :: NOTHIN = 1              ! No OBC
  integer, parameter :: CLAMPD = 2              ! Clamped to zero
  integer, parameter :: NOGRAD = 4              ! No gradient
  integer, parameter :: INGRAD = 8              ! Interior cell no gradient
  integer, parameter :: LINEAR = 16             ! Linear extrapolation
  integer, parameter :: FILEIN = 32             ! Input from file
  integer, parameter :: MEANIN = 64             ! Mean obc value
  integer, parameter :: ORLANS = 128            ! Orlanski radiation
  integer, parameter :: CAMOBR = 256            ! Camerlengo radiation
  integer, parameter :: MILLER = 512            ! Miller & Thorpe radiation
  integer, parameter :: IOW    = 1024           ! Schmidt radiation
  integer, parameter :: GRAVTY = 2048           ! Gravity wave radiation
  integer, parameter :: UPSTRM = 4096           ! Upstream advection
  integer, parameter :: RAYMND = 8192           ! Raymond and Kuo radiation
  integer, parameter :: FLATHR = 16384          ! Flather radiation
  
  ! parameter for enhanced viscosity or diffusion
  integer, parameter :: NONE     = 0
  integer, parameter :: MODERATE = 1            ! enhance against interiour
  integer, parameter :: LARGE    = 2            ! enhance to maximum value (CFL)

  real, parameter    :: small = 1.0e-8          ! some small number
  real               :: dtuv, dtts, dtbt, dteta ! time step.
  real               :: dtime_e                 ! timestep (secs) (dtime_e=2*dteta for threelevel and dtime_e=dteta for twolevel)
  real               :: dtime_t                 ! timestep (secs) (dtime_t=2*dtts for threelevel and dtime_t=dtts for twolevel)
  integer            :: isc, iec, jsc, jec      ! computational domain decompsition
  integer            :: isd, ied, jsd, jed      ! data domain decompsition
  integer            :: isg, ieg, jsg, jeg      ! global domain decompsition
  integer            :: nk                      ! number of vertical levels
  integer            :: tendency=0              ! integer corresponding to the choice of time tendency 

  ! for restart
  integer                       :: id_restart(4) = 0
  type(restart_file_type), save :: Obc_restart

  type(ocean_domain_type), pointer :: Dom => NULL()  ! ocean domain
  type(ocean_grid_type),pointer    :: Grd => NULL()  ! ocean grid
  integer, allocatable             :: obc_out_unit(:)

  !--- derived types to specify the domain and grid --------------------
  type obc_bound_type
     integer           :: is, ie, js, je         ! index specifying bound location on current pe
     integer           :: isd,ied,jsd,jed        ! index specifying bound location on current pe of data domain
     integer           :: isg,ieg,jsg,jeg        ! global index specifying bound location 
     integer           :: iers,iere,jers,jere    ! global indecees of eta files
     integer           :: itrs,itre,jtrs,jtre    ! global indecees of tracer files
     integer           :: np                     ! specifies total number of boundary points
     integer           :: direction              ! its value can be "w"(west), "e"(east), "s"(south), "n"(north)     
     logical           :: on_bound               ! true means there are some points on the boundary on current pe.
     logical           :: report                 ! true means this PE writes a report on the configuration
     integer           :: bcond_nor              ! Normal velocity OBC
     integer           :: bcond_tan              ! Normal velocity OBC
     integer           :: bcond_eta              ! Surface elevation OB
     integer           :: bcond_ud               ! Normal depth integrated velocity OBC 
     integer           :: bcond_eta_t            ! Surface height OBC
     integer, pointer  :: bcond_tra(:)           ! Tracer OBC
     integer           :: bcond_mix              ! Vertical mixing coeff. OBC

     !-----------------------------------
     ! Computational domain boundary maps
     integer           :: nloc                   ! Number of boundary points

     ! Boundary i and j vectors containing the (i,j) locations for 
     ! the computational domain OBC.
     ! Note that for east and west OBCs   iloc = Obc%bound%is is constant, 
     ! and for north and south boundaries jloc = Obc%bound%js is constant.
     integer, pointer  :: iloc(:)    => NULL()   ! Boundary i location
     integer, pointer  :: jloc(:)    => NULL()   ! Boundary j location

     ! Boundary i and j vectors containing the interior cell immediately
     ! adjacent to the boundary, i.e.
     ! For western boundaries oi1 = iloc + 1, oj1 => jloc
     ! For eastern boundaries oi1 = iloc - 1, oj1 => jloc
     ! For southern boundaries oj1 = jloc + 1, oi1 => iloc
     ! For northern boundaries oj1 = jloc - 1, oi1 => iloc
     integer, pointer  :: oi1(:)     => NULL()   ! Boundary interior i location
     integer, pointer  :: oj1(:)     => NULL()   ! Boundary interior j location

     ! Boundary i and j vectors containing interior cells 2 cells distant
     ! from the boundary, i.e.
     ! For western boundaries oi2 = iloc + 2, oj1 => jloc
     ! For eastern boundaries oi2 = iloc - 2, oj1 => jloc
     ! For southern boundaries oj2 = jloc + 2, oi1 => iloc
     ! For northern boundaries oj2 = jloc - 2, oi1 => iloc
     integer, pointer  :: oi2(:)     => NULL()   ! Interior 2i location
     integer, pointer  :: oj2(:)     => NULL()   ! Interior 2j location

     ! Boundary vectors containing the variable index only for use with
     ! Obc%eta%data.
     ! For western and eastern boundaries, d1i_bt = 1, d1j_bt => jloc
     ! For southern and northern boundaries, d1i_bt => iloc, d1j_bt = 1
     integer, pointer  :: d1i_bt(:)  => NULL()   ! 1D i index pointer
     integer, pointer  :: d1j_bt(:)  => NULL()   ! 1D j index pointer

     ! Boundary vectors containing the face centered boundary location.
     ! For eastern and northern boundaries this corresponds to the 
     ! interior cell. This corresponds to the NVIE location (Herzfeld, 2008).
     ! For western boundaries icf => iloc, jcf => jloc
     ! For eastern boundaries icf => oi1, jcf => jloc
     ! For southern boundaries icf => iloc, jcf => jloc
     ! For northern boundaries icf => iloc, jcf => oi1
     integer, pointer  :: icf(:)     => NULL()   ! Boundary face i location
     integer, pointer  :: jcf(:)     => NULL()   ! Boundary face j location

     ! Boundary vectors containing the face outside the elevation location.
     ! For western and sorthern boundaries this is iloc - imap, jloc - jmap.
     ! This corresponds to the NVOE location (Herzfeld, 2008).
     ! For western boundaries ico_bt => iloc-imap, jco_bt => jlod
     ! For eastern boundaries ico_bt => iloc, jco_bt => jlod
     ! For southern boundaries ico_bt => iloc, jco_bt => jloc-jmap
     ! For northern boundaries ico_bt => iloc, jco_bt => jloc

     integer, pointer  :: ico_bt(:)     => NULL()   ! Outside face i location
     integer, pointer  :: jco_bt(:)     => NULL()   ! Outside face j location

     ! Boundary vectors containing the normal (nic) and tangential (tic)
     ! indicies to the boundary.
     ! For west and east boundaries nic => iloc, tic => jloc
     ! For south and west boundaries nic => jloc, tic => iloc
     integer, pointer  :: nic(:)     => NULL()   ! Pointer to normal index
     integer, pointer  :: tic(:)     => NULL()   ! Pointer to tangential index

     integer, pointer  :: ones_bt(:) => NULL()   ! Array of ones

     ! Boundary interior maps :
     ! For western boundaries imap = 1, jmap = 0
     ! For eastern boundaries imap = -1, jmap = 0
     ! For southern boundaries imap = 0, jmap = 1
     ! For northern boundaries imap = 0, jmap = -1
     integer           :: imap                   ! Boundary i interior map
     integer           :: jmap                   ! Boundary j interior map

     ! Boundary tangential maps :
     ! For western and eastern boundaries imapt = 0, jmapt = 1
     ! For northern and southern boundaries imapt = 1, jmapt = 0
     integer           :: imapt                   ! Boundary i tangential map
     integer           :: jmapt                   ! Boundary j tangential map

     ! Start and end coordinates for the boundary.
     ! Western and eastern boundaries : nsc = jsc, nec = jec
     ! Southern and northern boundaries : nsc = isc, nec = iec
     integer           :: nsc                    ! Start coordinate
     integer           :: nec                    ! End coordinate
     real              :: sign                   ! -1 for east/north

     !-----------------------------------
     ! Data domain boundary maps
     integer           :: nlod                   ! Number of boundary points

     ! Boundary i and j vectors containing the (i,j) locations for 
     ! the data domain OBC. Analogous to (iloc, jloc).
     integer, pointer  :: ilod(:)    => NULL()   ! Boundary i location
     integer, pointer  :: jlod(:)    => NULL()   ! Boundary j location

     ! Boundary i and j vectors containing the interior cell immediately
     ! adjacent to the boundary for the data domain. Analogous to (oi1, oj1).
     integer, pointer  :: di1(:)     => NULL()   ! Boundary interior i location
     integer, pointer  :: dj1(:)     => NULL()   ! Boundary interior j location

     ! Boundary i and j vectors containing interior cells 2 cells distant
     ! from the boundary for the data domain. Analogous to (oi2, oj2).
     integer, pointer  :: di2(:)     => NULL()   ! Interior 2i location
     integer, pointer  :: dj2(:)     => NULL()   ! Interior 2j location

     ! Boundary vectors containing the face centered boundary location for
     ! the data domain. Analogous to (icf, jcf).
     integer, pointer  :: idf(:)     => NULL()   ! Boundary face i location
     integer, pointer  :: jdf(:)     => NULL()   ! Boundary face j location

     ! Boundary vectors containing the face outside the elevation location for
     ! the data domain. Analogous to (ico, jco).
     integer, pointer  :: ido(:)     => NULL()   ! Outside face i location
     integer, pointer  :: jdo(:)     => NULL()   ! Outside face jlocation

     ! Boundary vectors containing the normal (nic) and tangential (tic)
     ! indicies to the data domain boundary. Analogous to (nic, tic).
     integer, pointer  :: nid(:)     => NULL()   ! Pointer to normal index
     integer, pointer  :: tid(:)     => NULL()   ! Pointer to tangential index

     ! Boundary vectors containing start and end indicies for the halo
     ! cells associated with the data domain boundary.
     ! For western boundaries : istr = max(isd,isd-2), iend = isd-1
     !                          jstr = jend => jlod
     ! For eastern boundaries : istr = isd+1, iend = min(ied,i+2)
     !                          jstr = jend => jlod
     ! For southern boundaries : jstr = max(jsd,j-2), jend = jsd-1
     !                           istr = iend => ilod
     ! For northern boundaries : jstr = jsd+1, jend = min(jed,j+2)
     !                           istr = iend => ilod
     integer, pointer  :: istr(:)    => NULL()   ! Boundary i halo start map
     integer, pointer  :: iend(:)    => NULL()   ! Boundary i halo end map
     integer, pointer  :: jstr(:)    => NULL()   ! Boundary j halo start map
     integer, pointer  :: jend(:)    => NULL()   ! Boundary j halo end map

     ! Start and end coordinates for the data domain boundary.
     ! Western and eastern boundaries : nsd = jsd, ned = jed
     ! Southern and northern boundaries : nsd = isd, ned = ied
     integer           :: nsd                    ! Start coordinate
     integer           :: ned                    ! End coordinate

     real,    pointer  :: dum(:)     => NULL()   ! 1D dummy array
     real,    pointer  :: pgrad(:)   => NULL()   ! barotropic phase speed
     real,    pointer  :: hadv (:,:) => NULL()   ! cross boundary horizontal tracer advection, 
                                                 ! updated in the tracer loop, so no tracer number needed
     logical, pointer  :: mask(:)    => NULL()   ! logical varible to indicate the grid is open or not.
     character(len=16) :: name                   ! 
     real              :: ctrop_max              ! maximum barotropic phase speed (should be 1 .. 2) 
     real              :: ctrop_min              ! minimum barotropic phase speed (should be 0 .. 1) 
     real              :: ctrop_inc              ! barotropic phase speed for incoming waves (times cgrid) (should be 0 .. 1) 
     real              :: ctrop_smooth           ! smooth parameter for barotropic phase speed (should be 0 .. 1) 
     integer           :: enh_pnts               ! number of points besides the boundary with enhanced viscosity
     real              :: enh_fac_v              ! enhancement factor of viscosity
     real              :: enh_fac_d              ! enhancement factor of diffusivity
     integer           :: obc_enhance_visc_back  ! enhance viscosity near boundary
     integer           :: obc_enhance_diff_back  ! enhance diffusion near boundary
     logical           :: obc_consider_convu     ! true to consider convu for d eta/dt
     logical           :: obc_vert_advel_t       ! consider vertical advection for tracers 
     logical           :: obc_vert_advel_u       ! consider vertical advection for velocity 
     logical           :: obc_adjust_forcing_bt  ! remove baroclinic press gradients from forcing_bt
     logical           :: obc_relax_eta          ! true to relax sea level
     logical           :: obc_relax_eta_profile  ! true to relax sea level comparing mean sea level with reference value
     logical           :: obc_damp_newton        ! true to apply newtonian damping
     logical           :: obc_ud_active          ! Set ud to active forcing (FLATHR or FILEIN)
     real              :: damp_factor            ! newton damping factor
     logical, pointer  :: obc_relax_tracer(:)    => NULL() ! true to relax a tracer
     logical, pointer  :: obc_consider_sources(:)=> NULL() ! true for valid source and SGS terms
     integer, pointer  :: obc_flow_relax(:)      => NULL() ! true to invoke flow relaxation
     logical, pointer  :: obc_tracer_no_inflow(:)=> NULL() ! true to treat a tracer with Orlanski scheme
     real,    pointer  :: buffer(:)              => NULL() ! used for bitwise average on the boundary
     integer, pointer  :: start(:)               => NULL() ! list of starting index on the boundary
     integer, pointer  :: end(:)                 => NULL() ! list of ending index on the boundary
     integer, pointer  :: pelist(:)              => NULL() ! pelist used for this boundary
     integer           :: index
     real, pointer     :: work(:)  => NULL()     ! storage buffer
     real, pointer     :: work2(:,:,:) => NULL() ! 2d - storage buffer (1 index dummy)
     real, pointer     :: dsur(:) => NULL()       !kk Gravity wave speed
     integer           :: id_xt, id_yt, id_xu, id_yu
     integer           :: id_ctrop, id_eta_data, id_rel_coef, id_transport
     integer, allocatable :: id_tracer_data(:), id_tracer_flux(:), id_cclin(:)
  end type obc_bound_type

  type obc_data_type_2d
     real, pointer     :: data(:,:) => NULL() ! boundary values of eta for relaxation
     character(len=128):: file, field         ! name of the inputfile and inputfield
     real              :: rel_coef_in         ! relaxation coefficient to data during inflow
     real              :: rel_coef_out        ! relaxation coefficient to data during outflow
     integer           :: id                  ! the data id for time interpolation
     integer           :: rel_pnts            ! defines width of the band near the boundary for relaxation     
  end type obc_data_type_2d

  type obc_data_type_3d
     real, pointer     :: data(:,:,:) => NULL() ! boundary values of eta for relaxation
     character(len=128):: file, field         ! name of the inputfile and inputfield
     real              :: rel_coef_in         ! relaxation coefficient to data during inflow
     real              :: rel_coef_out        ! relaxation coefficient to data during outflow
     integer           :: id                  ! the data id for time interpolation
     integer           :: rel_pnts            ! defines width of the band near the boundary for relaxation     
  end type obc_data_type_3d

  type ocean_obc_type
     integer                       :: nobc               ! number of boundaries
     type(obc_bound_type), pointer :: bound(:) => NULL() ! contains boundary information.
     type(obc_data_type_2d),  pointer :: eta(:)      => NULL() ! contains boundary data information for sea level.
     type(obc_data_type_3d),  pointer :: tracer(:,:) => NULL() ! contains boundary data information for tracers.
     type(obc_data_type_2d),  pointer :: ud(:) => NULL() ! contains boundary data information for normal depth 
                                                         ! integrated velocity component
     type(obc_data_type_3d),  pointer :: uvel(:,:) => NULL() ! contains boundary data information for normal velocity component     
  end type ocean_obc_type

  type(ocean_obc_type), save :: Obc

  !<INTERFACE NAME="obc_update_boundary">
  !   <DESCRIPTION>
  !   update field on the halo points at the global boundaries.
  !   </DESCRIPTION>
  !   <INOUT NAME="field">
  !      field to be update on the boundary
  !   </INOUT>
  interface ocean_obc_update_boundary
     module procedure ocean_obc_update_boundary_2d
     module procedure ocean_obc_update_boundary_3d
  end interface
  !</INTERFACE>

  !<INTERFACE NAME="obc_zero_boundary">
  !   <DESCRIPTION>
  !   set field at open boundaries to zero.
  !   </DESCRIPTION>
  !   <INOUT NAME="field">
  !      field to be set to zero on the boundary
  !   </INOUT>
  interface ocean_obc_zero_boundary
     module procedure ocean_obc_zero_boundary_2d
  end interface
  !</INTERFACE>

  !--- namelist interface ---------------------------------------------
  !<NAMELIST NAME="ocean_obc_nml">
  !  <DATA NAME="nobc" TYPE="integer" >
  !    number of open boundary condition. Its value should be less than max_obc. Increase max_obc if needed.
  !  </DATA>
  !  <DATA NAME="direction" TYPE="character(len=10), dimension(max_obc)" >
  !    open boundary direction. Each element value should be west, east, south or north.
  !  </DATA>
  !  <DATA NAME="is, ie, js, je" TYPE="integer, dimension(max_obc)" >
  !    open boundary position. 
  !  </DATA>
  !  <DATA NAME="name" TYPE="character(len=32), dimension(max_obc)" >
  !    type of open bounday.
  !  </DATA>
  !  <DATA NAME="obc_nor" TYPE="character, dimension(max_obc)" >
  !    Normal velocity OBC
  !  </DATA>
  !  <DATA NAME="obc_tan" TYPE="character, dimension(max_obc)" >
  !    Tangential velocity OBC
  !  </DATA>
  !  <DATA NAME="obc_eta" TYPE="character, dimension(max_obc)" >
  !    Surface elevation OBC
  !  </DATA>
  !  <DATA NAME="obc_tra" TYPE="character, dimension(max_obc,max_prog_tracers)" >
  !    Tracers OBC
  !  </DATA>
  !  <DATA NAME="obc_mix" TYPE="character, dimension(max_obc)" >
  !    Vertical mixing coefficient OBC
  !  </DATA>
  !  <DATA NAME="obc_relax_eta" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether relax eta or not. Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_consider_convu" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether to account for one 
  !    component of convu within the boundary. The appropriate behavior
  !    depends on the model configuration.
  !    Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_vert_advel_t" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether to account for vertical
  !    advection of tracers at the boundary. The appropriate behavior
  !    depends on the model configuration.
  !    Default value is .false. (Currently inactive)
  !  </DATA>
  !  <DATA NAME="obc_vert_advel_u" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether to account for vertical
  !    advection of momentum at the boundary. The appropriate behavior
  !    depends on the model configuration.
  !    Default value is .false. 
  !  </DATA>
  !  <DATA NAME="obc_relax_eta_profile" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether relax eta to a prescribed profile or not. 
  !    Default value is .false. In this case  only the average sea level is relaxed,
  !    the profile, hence the geostrophic current, is unchanged.
  !  </DATA>
  !  <DATA NAME="obc_relax_tracer" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether relax tracer or not. Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_flow_relax" TYPE="integer, dimension(max_obc)" >
  !    Integer variable specifying the flow relaxation zone 
  !    (flow realxation of Martinsen and Engedahl (1987). Default value is 1.
  !  </DATA>
  !  <DATA NAME="obc_consider_sources" TYPE="logical, dimension(max_obc)" >
  !    Logical variable specifying if source and SGS terms of the normal tracer
  !    scheme are valid. Default value is .false..
  !  </DATA>
  !  <DATA NAME="obc_tracer_no_inflow" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether apply orlanski obc on tracer or not. Default value is .false.
  !  </DATA>
  !  <DATA NAME="ctrop_max" TYPE="real, dimension(max_obc)" >
  !    Maximum value to clip diagnosed barotropic phase speed in terms of sqrt(gH).
  !    Should be about  1. 
  !  </DATA>
  !  <DATA NAME="ctrop_min" TYPE="real, dimension(max_obc)" >
  !    Minimum value to diagnosed barotropic phase speed in terms of sqrt(gH).
  !    Should be about  0.  Default is 0.1.
  !  </DATA>
  !  <DATA NAME="ctrop_inc" TYPE="real, dimension(max_obc)" >
  !    value to be set for barotropic phase speed if incoming waves are diagnosed.
  !    (in terms of sqrt(gH)) Should be about 0. Default is 0. 
  !  </DATA>
  !  <DATA NAME="rel_coef_eta_in" TYPE="real, dimension(max_obc)" >
  !    Relaxation coefficient to be used for incoming wave situation. 
  !  </DATA>
  !  <DATA NAME="rel_coef_eta_out" TYPE="real, dimension(max_obc)" >
  !    Relaxation coefficient to be used for outgoing wave situation. 
  !    Should be smaller then or equal to rel_coef_eta_in. 
  !  </DATA>
  !  <DATA NAME="filename_eta" TYPE="character, dimension(max_obc)" >
  !    Filename to read sea level data.
  !  </DATA>
  !  <DATA NAME="fieldname_eta" TYPE="character, dimension(max_obc)" >
  !    Fieldname  to read sea level data.
  !  </DATA>
  !  <DATA NAME="rel_eta_pnts" TYPE="integer, dimension(max_obc)" >
  !    Relax sea level at a stripe of rel_eta_pnts. Default = 1.
  !  </DATA>
  !  <DATA NAME="rel_coef_tracer_in" TYPE="real, dimension(max_obc,max_prog_tracers)" >
  !    Relaxation coefficient to be used for inflow situation. 
  !  </DATA>
  !  <DATA NAME="rel_coef_tracer_out" TYPE="real, dimension(max_obc,max_prog_tracers)" >
  !    Relaxation coefficient to be used for outflow situation. 
  !    Should be smaller then or equal to rel_coef_tracer_in. 
  !  </DATA>
  !  <DATA NAME="rel_clin_pnts" TYPE="integer, dimension(max_obc,max_prog_tracers)" >
  !    Relax a tracer at a stripe of rel_clin_pnts. Default = 1.
  !  </DATA>
  !  <DATA NAME="filename_tracer" TYPE="character, dimension(max_obc,max_prog_tracers)" >
  !    Filename to read a tracer. It is allowed to put all data for a boundary in one file.
  !  </DATA>
  !  <DATA NAME="fieldname_tracer" TYPE="character, dimension(max_obc,max_prog_tracers)" >
  !    Fieldname of a tracer. 
  !  </DATA>
  !  <DATA NAME="obc_enhance_visc_back" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether to enhance viscosity at the boundary. Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_enhance_diff_back" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether to enhance mixing at the boundary. Default value is .false.
  !  </DATA>
  !  <DATA NAME="enh_fac_v" TYPE="real, dimension(max_obc)" >
  !    'Safety factor' applied to maximum stable viscosity at the boundary. Default = 0.9
  !  </DATA>
  !  <DATA NAME="enh_fac_d" TYPE="real, dimension(max_obc)" >
  !    Factor applied to enhance mixing at the boundary. Default = 1.
  !  </DATA>
  !  <DATA NAME="enh_pnts" TYPE="integer, dimension(max_obc)" >
  !    Enhance viscosity and mixing at a stripe of enh_pnts
  !    decreasing with the distance from the boundary. Default = 1.
  !  </DATA>
  !  <DATA NAME="debug_phase_speed" TYPE="logical">
  !  Includes the phase speed into the model output.
  !  </DATA> 
  !  <DATA NAME="debug_this_module" TYPE="logical">
  !  For debugging.
  !  </DATA> 
  !  <DATA NAME="nt" TYPE="integer">
  !    number of tracers to use open boundary condition. 
  !  </DATA>
  !</NAMELIST>

  integer                               :: nobc= 0
  character(len=10), dimension(max_obc) :: direction='' ! open directions; to be consistent with is, ie, js, je
  integer, dimension(max_obc)           :: is=-999, ie=-999, js=-999, je=-999     ! boundary position
  integer, dimension(max_obc)           :: iers=-999, iere=-999, jers=-999, jere=-999     ! boundary position
  integer, dimension(max_obc)           :: itrs=-999, itre=-999, jtrs=-999, jtre=-999     ! boundary position
  character(len=32), dimension(max_obc) :: name 
  character(len=128), dimension(max_obc):: obc_nor = 'NOGRAD'
  character(len=128), dimension(max_obc):: obc_tan = 'NOGRAD'
  character(len=128), dimension(max_obc):: obc_eta = 'NOTHIN'
  character(len=128), dimension(max_obc):: obc_height = 'NOTHIN'
  character(len=128), dimension(max_obc):: obc_ud = 'NOGRAD'  
  character(len=128), dimension(max_obc):: obc_mix = 'NOGRAD'
  character(len=128), dimension(max_obc,max_prog_tracers):: obc_tra = 'NOGRAD'
  character(len=128), dimension(max_obc):: obc_enhance_visc_back = 'NONE'
  character(len=16),  dimension(max_obc):: obc_enhance_diff_back = 'NONE'
  logical, dimension(max_obc)           :: obc_damp_newton = .FALSE.
  logical, dimension(max_obc)           :: obc_consider_convu = .FALSE.
  logical, dimension(max_obc)           :: obc_vert_advel_t= .FALSE.
  logical, dimension(max_obc)           :: obc_vert_advel_u= .FALSE.
  logical, dimension(max_obc)           :: obc_adjust_forcing_bt = .FALSE.
  real,    dimension(max_obc)           :: rel_coef_eta_in  = 0.0           ! relaxation coefficient, specify in namelist
  real,    dimension(max_obc)           :: rel_coef_eta_out = 0.0           ! relaxation coefficient, specify in namelist
  real,    dimension(max_obc)           :: ctrop_max = 1.5                  ! will be multiplied with sqrt(g*H)
  real,    dimension(max_obc)           :: ctrop_min = 0.1                  ! will be multiplied with sqrt(g*H)
  real,    dimension(max_obc)           :: ctrop_inc = .0                   ! will be multiplied with cgrid
  real,    dimension(max_obc)           :: ctrop_smooth = 0.7
  real,    dimension(max_obc)           :: damp_factor = 1.
  integer, dimension(max_obc)           :: enh_pnts  = 1                     
  integer, dimension(max_obc)           :: rel_eta_pnts  = 1                     
  real,    dimension(max_obc)           :: enh_fac_v   = 0.9
  real,    dimension(max_obc)           :: enh_fac_d   = 1.0
  character(len=256), dimension(max_obc):: filename_eta
  character(len=32),  dimension(max_obc):: fieldname_eta
  character(len=256), dimension(max_obc):: filename_ud
  character(len=32),  dimension(max_obc):: fieldname_ud  
  logical,            dimension(max_obc,max_prog_tracers):: obc_relax_tracer = .FALSE.
  integer,            dimension(max_obc,max_prog_tracers):: obc_flow_relax = 1
  logical,            dimension(max_obc,max_prog_tracers):: obc_consider_sources = .FALSE.
  logical,            dimension(max_obc,max_prog_tracers):: obc_tracer_no_inflow = .FALSE.
  integer,            dimension(max_obc,max_prog_tracers):: rel_clin_pnts  = 1            ! relax at rel_clin_pnts                 
  real,               dimension(max_obc,max_prog_tracers):: rel_coef_tracer_in  = 0.0     ! relaxation coefficient
  real,               dimension(max_obc,max_prog_tracers):: rel_coef_tracer_out = 0.0     ! relaxation coefficient
  character(len=256), dimension(max_obc,max_prog_tracers):: filename_tracer
  character(len=32),  dimension(max_obc,max_prog_tracers):: fieldname_tracer
  logical                               :: debug_this_module = .FALSE.
  logical                               :: debug_phase_speed = .FALSE.    
  logical                               :: update_eta_tm1    = .FALSE. 
  integer, parameter                    :: mobcm1 = max_obc-1
  integer, parameter                    :: ntm2   = max_prog_tracers-2
  integer, parameter                    :: tobc   = max_obc*max_prog_tracers
  integer                               :: ne, te

  character(len=128)                    :: errorstring
  
  integer                               :: id_obc
  
  type(domain2d),allocatable, save      :: obc_eta_domain(:)
  type(domain2d),allocatable, save      :: obc_diag_domain(:)
  
  data (name(ne), ne=1,max_obc)          /'test_obc', mobcm1*'none'/
  data (fieldname_eta(ne), ne=1,max_obc) /'eta_t', mobcm1*'none'/
  data (filename_eta(ne), ne=1,max_obc)  /'obc_eta_t.nc', mobcm1*'none'/
  data (filename_ud(ne), ne=1,max_obc)  /'obc_ud.nc', mobcm1*'none'/
  data (fieldname_ud(ne), ne=1,max_obc) /'ud', mobcm1*'none'/  
  data ((fieldname_tracer(ne,te),  ne=1, max_obc), te=1, max_prog_tracers) &
           /tobc*'none'/
  data ((filename_tracer(ne,te),   ne=1, max_obc), te=1, max_prog_tracers) &
           /tobc*'none'/

  namelist /ocean_obc_nml/ nobc, direction, name, &
                           is, ie, js, je, &
                           iers, iere, jers, jere, &     
                           itrs, itre, jtrs, jtre, &    
                           obc_nor, obc_tan, obc_eta, obc_ud, obc_tra, &
                           obc_mix, rel_coef_eta_in, rel_coef_eta_out, &
                           rel_eta_pnts, rel_clin_pnts,      &
                           ctrop_max, ctrop_min, ctrop_inc, ctrop_smooth, &
                           filename_eta, fieldname_eta, &
                           filename_ud, fieldname_ud, &                           
                           obc_consider_convu, obc_adjust_forcing_bt, &
                           obc_vert_advel_t, obc_vert_advel_u, &
                           obc_enhance_visc_back, obc_enhance_diff_back, &
                           enh_pnts, enh_fac_v, enh_fac_d, &
                           obc_relax_tracer, obc_flow_relax, &
                           obc_consider_sources, obc_tracer_no_inflow, &
                           rel_coef_tracer_in, rel_coef_tracer_out, filename_tracer, fieldname_tracer,   &
                           debug_phase_speed, debug_this_module, obc_damp_newton, damp_factor
  logical              ::  south_west_corner=.false.
  logical              ::  south_east_corner=.false.
  logical              ::  north_west_corner=.false.
  logical              ::  north_east_corner=.false.
  integer              ::  i_sw, j_sw, i_nw, j_nw, i_se, j_se, i_ne, j_ne

  character(len=128)   :: version = '$ID$'
  character (len=128)  :: tagname = '$Name: tikal $'
  logical              :: module_is_initialized = .FALSE.
  real, allocatable    :: wrk2(:)         ! needed for enhanced diffusion
  real, allocatable    :: wrk3(:)         ! needed for enhanced diffusion
  real, allocatable    :: eta_tend_obc(:,:), rel_coef(:,:)   
  real, allocatable, dimension (:,:)        :: t_mask, u_mask, ht_bt
  real, allocatable, dimension (:,:),target :: dxu_bt, dyu_bt
  real, allocatable, dimension (:,:),target :: eta_t_tm1, eta_tm1
  real, allocatable, dimension (:,:)        :: ctrop      ! needed for restart of phase speed

contains

  !#######################################################################
  !  <SUBROUTINE NAME="ocean_obc_barotrop_init" >
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable that 
  !      contains domain decompostion and grid information.
  !   </DESCRIPTION>
  !   <IN NAME="dtts, dtuv, dtbt, dteta" TYPE="real"></IN>
  !    time step.
  !   <IN NAME="Domain" TYPE="type(ocean_domain_type)" >
  !      A derived data type that contains domain information for MOM.
  !   </IN>
  !   <IN NAME="Grid" TYPE="type(ocean_grid_type)" >
  !      A derived data type that contains grid information for MOM.
  !   </IN>
  !   <INOUT NAME="have_obc" TYPE="logical" >
  !      logical variable to indicate if there is any open boundary condition. 
  !      if true, open boudanry exists. 
  !   </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_barotrop_init(have_obc, Time, Time_steps, Domain, Grid, Ocean_options,   &
                           use_legacy_barotropic_halos, debug)
  !</PUBLICROUTINE>

    logical, intent(inout)                       :: have_obc
    type(ocean_time_type), intent(in)            :: Time
    type(ocean_time_steps_type), intent(in)      :: Time_steps
    type (ocean_domain_type),intent(in),  target :: Domain
    type(ocean_grid_type), intent(in), target    :: Grid
    type(ocean_options_type), intent(inout)      :: Ocean_options
    logical, intent(in), optional                :: debug
    logical, intent(in)                          :: use_legacy_barotropic_halos
    !--- some local variables ------------------------------------------
    integer :: m, n, i, j, unit, io_status, ierr, ioff, joff
    integer :: il, iu, jl, ju
    integer :: ni, nj, nsize
    integer :: west_at_pe, south_at_pe, east_at_pe, north_at_pe
    integer :: irig_s, ilef_s , jbou_s
    integer :: irig_n, ilef_n , jbou_n
    integer :: jlow_w, jup_w, ibou_w
    integer :: jlow_e, jup_e, ibou_e
    integer :: npes, pe, nlist, list, ntotal, sindex, eindex
    integer :: fld_size(4)

    !--- some local variables ------------------------------------------
    integer, allocatable :: isl(:), iel(:), jsl(:), jel(:), istart(:), iend(:), pelist(:)
    logical, allocatable :: on_bound(:)
    character*128        :: file_name
    character(len=5)     :: pe_name
    
    integer :: stdoutunit,stdlogunit
    stdoutunit=stdout();stdlogunit=stdlog()    

    if ( module_is_initialized ) then 
       call mpp_error(FATAL, '==>Error in ocean_obc_barotrop_mod (ocean_obc_barotrop_init): module already initialized')
    endif

    module_is_initialized = .TRUE.

    Dom => Domain
    Grd => Grid
    call mpp_get_global_domain(Domain%Domain2d, isg, ieg, jsg, jeg)

    if (PRESENT(debug)) debug_this_module = debug

    id_obc  = mpp_clock_id('(Ocean open boundaries) ',grain=CLOCK_MODULE)
    

    !--- read namelist and write out namelist --------------------------
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_obc_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_obc_nml')
#else
    unit = open_namelist_file()
    read  (unit, ocean_obc_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_obc_nml')
    call close_file (unit)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_obc_nml)
    write (stdlogunit, ocean_obc_nml)

    !--- if there is no open boundary, just return
    Obc%nobc = nobc
    if(Obc%nobc == 0) then
       Ocean_options%OBC = 'Did NOT use any open/radiating boundary conditions.'
       return
    endif
    have_obc = .true.
    Ocean_options%OBC = 'Used open/radiating boundary conditions.'

    call mpp_clock_begin(id_obc)
! This is the boundary specific io_unit, which will be defined only for 
! tasks with points at boundary m. Initialize with stdoutunit. All tasks with
! no OBC write (if they do) to stdout.  
    allocate(obc_out_unit(nobc))
    obc_out_unit(:) = stdoutunit
    !--- number of boundaries should not be over the max_obc
    if(nobc .gt. max_obc) call mpp_error(FATAL,'==>Error in ocean_obc_barotrop_mod: nml nobc is greater than maximum'// &
         ' allowed bounds max_obc. Modify nobc or increase "max_obc" at top of ocean_obc_barotrop_mod' )

    !--- cyclic condition and open boundary condition along zonal direction 
    !--- these two boundary conditions cannot exist at the same place 
    do m = 1, nobc
       if((trim(direction(m)) == 'west' .or. trim(direction(m)) == 'east')  .and. Grid%cyclic_x) &
            call mpp_error(FATAL, "==>Error in ocean_obc_barotrop_mod: when west or east boundary is open, "//&
                        "cyclic condition in i-direction cannot exist")
       if((trim(direction(m)) == 'south' .or. trim(direction(m)) == 'north')  .and. Grid%cyclic_y) &
            call mpp_error(FATAL, "==>Error in ocean_obc_barotrop_mod: when south or north boundary is open, "//&
                        "cyclic condition in j-direction cannot exist")
    enddo

    !--- get the domain decomposition -----------------------------------
    npes = mpp_npes()
    allocate(isl(0:npes-1), iel(0:npes-1), jsl(0:npes-1), jel(0:npes-1) )
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)

    call mpp_get_compute_domains(Domain%domain2d, xbegin=isl, xend=iel, ybegin=jsl, yend=jel)
    call get_domain_offsets(Domain, ioff, joff)
    do m = 1, nobc
       is(m) = is(m) - ioff
       ie(m) = ie(m) - ioff
       js(m) = js(m) - joff
       je(m) = je(m) - joff
    enddo

    !--- get time step -------------------------------------------------
    dtts    = Time_steps%dtts
    dtuv    = Time_steps%dtuv
    dtbt    = Time_steps%dtbt
    dteta   = Time_steps%dteta
    dtime_t = Time_steps%dtime_t
    dtime_e = Time_steps%dtime_e
    tendency= Time_steps%tendency

    ! number of vertical grid points 
    nk      = size(Grid%zw)
    allocate(wrk2(1:nk))
    allocate(wrk3(1:nk))

!kk Get data on _bt domain    
    allocate (t_mask(isd:ied,jsd:jed))
    allocate (u_mask(isd:ied,jsd:jed))
    allocate (dxu_bt(isd:ied,jsd:jed))
    allocate (dyu_bt(isd:ied,jsd:jed))
    allocate (ht_bt(isd:ied,jsd:jed))
    t_mask     = 0.0
    u_mask     = 0.0
    dxu_bt     = 0.0
    dyu_bt     = 0.0
    ht_bt      = 0.0
    il = lbound(Grd%tmask,1)
    jl = lbound(Grd%tmask,2)
    iu = ubound(Grd%tmask,1)
    ju = ubound(Grd%tmask,2)
    t_mask(il:iu,jl:ju)   = Grd%tmask(:,:,1)
    u_mask(il:iu,jl:ju)   = Grd%umask(:,:,1)
    dxu_bt(il:iu,jl:ju)   = Grd%dxu(:,:)
    dyu_bt(il:iu,jl:ju)   = Grd%dyu(:,:)
    ht_bt (il:iu,jl:ju)   = Grd%ht(:,:)
    call mpp_update_domains(t_mask,   Dom%domain2d)
    call mpp_update_domains(u_mask,   Dom%domain2d)
    call mpp_update_domains(dxu_bt,   Dom%domain2d)
    call mpp_update_domains(dyu_bt,   Dom%domain2d)
    call mpp_update_domains(ht_bt,    Dom%domain2d)

    if(ied-isd+1 /= size(t_mask,1) .or. jed-jsd+1 /= size(t_mask,2)) then
       write (0,'(a,7i6)') 'illegal mask size ',ied,isd,size(t_mask,1), jed,jsd,size(t_mask,2),mpp_pe()
       call mpp_error(FATAL,trim(errorstring)//'illegal mask size')
    end if

    !--- Initialize Obc data type --------------------------------------
    allocate(Obc%bound (nobc))
    allocate(obc_eta_domain(nobc))
    allocate(obc_diag_domain(nobc))
    do m = 1, nobc
       call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Dom%layout &
                               , obc_eta_domain(m), maskmap=Dom%maskmap)
       call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Dom%layout &
                               , obc_diag_domain(m), maskmap=Dom%maskmap ) 
    enddo

    Obc%bound%is = -999; Obc%bound%ie = -1000    ! dummy number to indicate not on the boundary
    Obc%bound%js = -999; Obc%bound%je = -1000
    Obc%bound%on_bound = .FALSE.
    Obc%bound%report   = .FALSE.
    Obc%bound%obc_relax_eta = .false.
    Obc%bound%obc_relax_eta_profile = .false.

    allocate(istart(0:npes-1), iend(0:npes-1), on_bound(0:npes-1) )
    allocate(pelist(0:npes-1) )
    call mpp_get_current_pelist(pelist)

    write(pe_name,'(a,i4.4)' )'.', mpp_pe()   

    do m = 1, nobc
       Obc%bound(m)%name                  = trim(name(m))

       select case(trim(direction(m)))
       case('west')
          Obc%bound(m)%direction = WEST
       case('east')
          Obc%bound(m)%direction = EAST
       case('south')
          Obc%bound(m)%direction = SOUTH
       case('north')
          Obc%bound(m)%direction = NORTH
       case default
          call mpp_error(FATAL,'each element of nml direction should be west, east, south or north')
       end select

       ! Set the normal OBC
       select case(trim(obc_nor(m)))
       case('NOTHIN')
          Obc%bound(m)%bcond_nor = NOTHIN
       case('CLAMPD')
          Obc%bound(m)%bcond_nor = CLAMPD
       case('NOGRAD')
          Obc%bound(m)%bcond_nor = NOGRAD
       case('INGRAD')
          Obc%bound(m)%bcond_nor = INGRAD
       case('LINEAR')
          Obc%bound(m)%bcond_nor = LINEAR
       case default
          call mpp_error(FATAL,'each element of nml obc_nor should be NOTHIN, CLAMPD, NOGRAD, INGRAD or LINEAR')
       end select

       ! Set the normal barotropic velocity OBC
       Obc%bound(m)%obc_ud_active = .false.
       Obc%bound(m)%bcond_ud = 0
       if (index(trim(obc_ud(m)),'NOTHIN') /= 0) &
            Obc%bound(m)%bcond_ud = Obc%bound(m)%bcond_ud + NOTHIN
       if (index(trim(obc_ud(m)),'FLATHR') /= 0) then
          Obc%bound(m)%bcond_ud = Obc%bound(m)%bcond_ud + FLATHR + NOGRAD
          Obc%bound(m)%obc_ud_active = .true.
       endif
       if (index(trim(obc_ud(m)),'FILEIN') /= 0) then
          Obc%bound(m)%bcond_ud = Obc%bound(m)%bcond_ud + FILEIN
          Obc%bound(m)%obc_ud_active = .true.
       endif
       if (Obc%bound(m)%bcond_ud .eq. 0) then
          select case(trim(obc_ud(m)))
          case('CLAMPD')
             Obc%bound(m)%bcond_ud = CLAMPD
          case('NOGRAD')
             Obc%bound(m)%bcond_ud = NOGRAD
          case('INGRAD')
             Obc%bound(m)%bcond_ud = INGRAD
          case('LINEAR')
             Obc%bound(m)%bcond_ud = LINEAR
          case default
             Obc%bound(m)%bcond_ud = Obc%bound(m)%bcond_nor
          end select
       endif
       if(Obc%bound(m)%bcond_ud == 0) &
            call mpp_error(FATAL,' Invalid barotropic velocity OBC for boundary '//trim(Obc%bound(m)%name))

       ! Set the tangential OBC
       select case(trim(obc_tan(m)))
       case('NOTHIN')
          Obc%bound(m)%bcond_tan = NOTHIN
       case('CLAMPD')
          Obc%bound(m)%bcond_tan = CLAMPD
       case('NOGRAD')
          Obc%bound(m)%bcond_tan = NOGRAD
       case('INGRAD')
          Obc%bound(m)%bcond_tan = INGRAD
       case('LINEAR')
          Obc%bound(m)%bcond_tan = LINEAR
       case default
          call mpp_error(FATAL,'each element of nml obc_tan should be NOTHIN, CLAMPD, NOGRAD, INGRAD or LINEAR')
       end select

       ! Set the eta OBC
       Obc%bound(m)%bcond_eta = 0
       if (index(trim(obc_eta(m)),'NOTHIN') /= 0) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + NOTHIN
       if (index(trim(obc_eta(m)),'FILEIN') /= 0) then
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + FILEIN
            Obc%bound(m)%obc_relax_eta = .true.
            Obc%bound(m)%obc_relax_eta_profile = .true.
       endif
       if (index(trim(obc_eta(m)),'MEANIN') /= 0) then
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + MEANIN
            Obc%bound(m)%obc_relax_eta = .true.
       endif
       if (index(trim(obc_eta(m)),'FLATHR') /= 0) then
          Obc%bound(m)%obc_relax_eta = .true.
          Obc%bound(m)%obc_relax_eta_profile = .true.
       endif
       if (index(trim(obc_eta(m)),'NOGRAD') /= 0) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + NOGRAD
       if (index(trim(obc_eta(m)),'ORLANS') /= 0) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + ORLANS
       if (index(trim(obc_eta(m)),'CAMOBR') /= 0) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + CAMOBR
       if (index(trim(obc_eta(m)),'GRAVTY') /= 0 ) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + GRAVTY
       if (index(trim(obc_eta(m)),'MILLER') /= 0 ) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + MILLER
       if (index(trim(obc_eta(m)),'IOW') /= 0) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + IOW
       if (index(trim(obc_eta(m)),'RAYMND') /= 0) &
            Obc%bound(m)%bcond_eta = Obc%bound(m)%bcond_eta + RAYMND
       call check_eta_OBC(Obc%bound(m)%bcond_eta, obc_eta(m))
       if (index(trim(obc_eta(m)),'RAYMND') /= 0) then
         if (.not.use_legacy_barotropic_halos) call mpp_error(FATAL,' barotropic OBC method RAYMND ' &
              //'does not work without use_legacy_barotropic_halos')
       endif
       

       ! Set the surface height OBC : this is currently always set to NOTHIN
       Obc%bound(m)%bcond_eta_t = 0
       if (index(trim(obc_height(m)),'NOTHIN') /= 0) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + NOTHIN
       if (index(trim(obc_height(m)),'ORLANS') /= 0) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + ORLANS
       if (index(trim(obc_height(m)),'CAMOBR') /= 0) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + CAMOBR
       if (index(trim(obc_height(m)),'GRAVTY') /= 0 ) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + GRAVTY
       if (index(trim(obc_height(m)),'MILLER') /= 0 ) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + MILLER
       if (index(trim(obc_height(m)),'IOW') /= 0) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + IOW
       if (index(trim(obc_height(m)),'RAYMND') /= 0) &
            Obc%bound(m)%bcond_eta_t = Obc%bound(m)%bcond_eta_t + RAYMND
       if(Obc%bound(m)%bcond_eta_t == 0) &
            call mpp_error(FATAL,' Invalid height OBC for boundary '//trim(Obc%bound(m)%name))

       Obc%bound(m)%ctrop_max             = ctrop_max(m)
       Obc%bound(m)%ctrop_min             = ctrop_min(m)
       Obc%bound(m)%ctrop_inc             = ctrop_inc(m)
       Obc%bound(m)%ctrop_smooth          = ctrop_smooth(m)
       Obc%bound(m)%enh_fac_v             = enh_fac_v(m)
       Obc%bound(m)%enh_fac_d             = enh_fac_d(m)
       Obc%bound(m)%obc_consider_convu    = obc_consider_convu(m)
       Obc%bound(m)%obc_damp_newton       = obc_damp_newton(m)
       Obc%bound(m)%damp_factor           = damp_factor(m)
       
       if (Obc%bound(m)%ctrop_max .lt. Obc%bound(m)%ctrop_min) then
         call mpp_error(FATAL,' ctrop_max <  ctrop_min for open boundary '//trim(Obc%bound(m)%name))
       endif

       if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
          !--- define mask to indicate if it is open or not.
          allocate(Obc%bound(m)%mask(jsd:jed))
          Obc%bound(m)%mask = .FALSE.
          if (is(m) .ge. isc .and. is(m) .le. iec) then
             do j = jsd, jed             
                if(j .ge. js(m) .and. j .le. je(m)) then 
                   Obc%bound(m)%mask(j) = .TRUE.
                endif   
             enddo
             do j = jsc, jec                 !kk changed from data to compute domain
                if(j .ge. js(m) .and. j .le. je(m)) then 
                   Obc%bound(m)%on_bound = .true.
                endif   
             enddo
             if(Obc%bound(m)%on_bound) then
!--- only the southernmost PE writes a report on the configuration
!___ if debugging is switched on each PE reports
                if (js(m) .ge. jsc .and. js(m) .le. jec) Obc%bound(m)%report = .true.
!--- open obc.out file to be ready for writing.
                if (debug_this_module .or. Obc%bound(m)%report) &
                     call mpp_open( obc_out_unit(m), 'obc_b_'//trim(Obc%bound(m)%name)//'.out', &
                                  action=MPP_OVERWR, threading=MPP_MULTI,                     &
                                  fileset=MPP_MULTI, nohdrs=.TRUE. )  
                Obc%bound(m)%is       = is(m)
                Obc%bound(m)%ie       = ie(m)
                Obc%bound(m)%js       = max(js(m),jsc)
                Obc%bound(m)%je       = min(je(m),jec)
                Obc%bound(m)%np       = je(m) - js(m) + 1
                Obc%bound(m)%isd      = Obc%bound(m)%is
                Obc%bound(m)%ied      = Obc%bound(m)%ie
                Obc%bound(m)%jsd      = Obc%bound(m)%js
                Obc%bound(m)%jed      = Obc%bound(m)%je
                Obc%bound(m)%isg      = is(m)
                Obc%bound(m)%ieg      = ie(m)
                Obc%bound(m)%jsg      = js(m)
                Obc%bound(m)%jeg      = je(m)
! Extend the data limit of OBC if the compute limit is at compute domain limit
                if(Obc%bound(m)%js > js(m) ) Obc%bound(m)%jsd = max(js(m),jsd)
                if(Obc%bound(m)%je < je(m) ) Obc%bound(m)%jed = min(je(m),jed)

                   if (debug_this_module .or. Obc%bound(m)%report) &
                write(obc_out_unit(m),'(a,10i6)') 'After Open ',mpp_pe(), &
                   is(m),ie(m),Obc%bound(m)%js,Obc%bound(m)%je,Obc%bound(m)%np,Obc%bound(m)%jsd,Obc%bound(m)%jed,jsd,jed

                allocate(Obc%bound(m)%pgrad(Obc%bound(m)%js:Obc%bound(m)%je))
! I no global domain information for reading boundary data arespecified in the 
! namelist assume, that the data fit exactly at the boundary
                if (iers(m) == -999) then 
                  Obc%bound(m)%iers = is(m) 
                else 
                  Obc%bound(m)%iers = iers(m) 
                endif
                if (iere(m) == -999) then 
                  Obc%bound(m)%iere = ie(m) 
                else 
                  Obc%bound(m)%iere = iere(m) 
                endif 
                if (jers(m) == -999) then 
                  Obc%bound(m)%jers = js(m) 
                else 
                  Obc%bound(m)%jers = jers(m) 
                endif 
                if (jere(m) == -999) then 
                  Obc%bound(m)%jere = je(m) 
                else 
                  Obc%bound(m)%jere = jere(m) 
                endif 
                if (Obc%bound(m)%jers > js(m)) &
                  call mpp_error(FATAL, "ocean_obc_barotrop_mod: start of eta data > start of boundary")
                if (Obc%bound(m)%jere < je(m)) &
                  call mpp_error(FATAL, "ocean_obc_barotrop_mod: end of eta data < end of boundary")
             endif
          end if
          !--- loop over pelist to figure out the boundary index on each pe.       
          on_bound = .false.
          do pe = 0, npes-1
             if (is(m) .ge. isl(pe) .and. is(m) .le. iel(pe)) then
                istart(pe) = max(js(m),jsl(pe)); iend(pe) = min(je(m),jel(pe))
                if(iend(pe) .GE. istart(pe) ) then
                   on_bound(pe) = .true.
                   if (debug_this_module .or. Obc%bound(m)%report) then
                      write(obc_out_unit(m),'(a,11i5)') &
                       'PE_1 ',pe,is(m),isl(pe),is(m),iel(pe),istart(pe),iend(pe),js(m),jsl(pe),je(m),jel(pe)
                   end if
                end if
             end if
          end do
          nlist = count(on_bound)
          allocate( Obc%bound(m)%start(nlist), obc%bound(m)%end(nlist)    )
          allocate( Obc%bound(m)%pelist(nlist) )
          list = 0
          ntotal = 0
          Obc%bound(m)%index = -1
          sindex = 100000; eindex = -100000
          do n = 0, npes-1
             if(on_bound(n)) then
                list = list+1
                Obc%bound(m)%pelist(list) = pelist(n)
                if(Obc%bound(m)%on_bound) then
                   Obc%bound(m)%start(list)  = istart(n)
                   Obc%bound(m)%end(list)    = iend(n)
                   sindex                    = min(sindex, istart(n))
                   eindex                    = max(eindex, iend(n))
                   ntotal                    = ntotal + iend(n) - istart(n) + 1
                   if( mpp_pe() == pelist(n)) then
                     Obc%bound(m)%index = list
                        if (debug_this_module .or. Obc%bound(m)%report) &
                     write(obc_out_unit(m),'(a,9i6)') &
                       'PE_2 ',n,list,Obc%bound(m)%start(list),Obc%bound(m)%end(list),sindex,eindex,ntotal,pelist(n),mpp_pe()
                   end if
                end if
             end if
          end do
          if(Obc%bound(m)%on_bound) then
             if(Obc%bound(m)%index == -1) then
                  call mpp_error(FATAL, &
                  "ocean_obc_barotrop_mod: This PE must be on this boundary and index should not be -1")
             end if
             if(Obc%bound(m)%js .NE. Obc%bound(m)%start(Obc%bound(m)%index) .OR. &
                  Obc%bound(m)%je .NE. Obc%bound(m)%end(Obc%bound(m)%index) ) &
                  call mpp_error(FATAL, "ocean_obc_barotrop_mod: inconsistency for the index setup")
             if(eindex-sindex+1 .NE. ntotal .OR. Obc%bound(m)%np .NE. ntotal) &
                    call mpp_error(FATAL, "ocean_obc_barotrop_mod: size mismatch 1")
             allocate(Obc%bound(m)%buffer(sindex:eindex) )
          end if
       else    ! north or south direction
          allocate(Obc%bound(m)%mask(isd:ied))
          Obc%bound(m)%mask = .FALSE.
          if(js(m) .le. jec .and. js(m) .ge. jsc) then
             do i = isd, ied
                if(i .ge. is(m) .and. i .le. ie(m)) then 
                   Obc%bound(m)%mask(i) = .TRUE.
                endif
             enddo
             do i = isc, iec                    !kk changed from data to compute domain
                if(i .ge. is(m) .and. i .le. ie(m)) then 
                   Obc%bound(m)%on_bound = .true.
                endif
             enddo
             if(Obc%bound(m)%on_bound) then
!--- open obc.out file to be ready for writing.
!--- only the southernmost PE writes a report on the configuration
!___ if debugging is switched on each PE reports
                if (is(m) .ge. isc .and. is(m) .le. iec) Obc%bound(m)%report = .true.
!--- open obc.out file to be ready for writing.
                if (debug_this_module .or. Obc%bound(m)%report) &
                    call mpp_open( obc_out_unit(m), 'obc_b_'//trim(Obc%bound(m)%name)//'.out', &
                                   action=MPP_OVERWR, threading=MPP_MULTI,                     &
                                   fileset=MPP_MULTI, nohdrs=.TRUE. )  
                Obc%bound(m)%js       = js(m)
                Obc%bound(m)%je       = je(m)
                Obc%bound(m)%is       = max(is(m),isc)
                Obc%bound(m)%ie       = min(ie(m),iec)
                Obc%bound(m)%np       = ie(m) - is(m) + 1
                Obc%bound(m)%isd      = Obc%bound(m)%is
                Obc%bound(m)%ied      = Obc%bound(m)%ie
                Obc%bound(m)%jsd      = Obc%bound(m)%js
                Obc%bound(m)%jed      = Obc%bound(m)%je
                Obc%bound(m)%isg      = is(m)
                Obc%bound(m)%ieg      = ie(m)
                Obc%bound(m)%jsg      = js(m)
                Obc%bound(m)%jeg      = je(m)
                if(Obc%bound(m)%is > is(m) ) Obc%bound(m)%isd = max(is(m),isd)
                if(Obc%bound(m)%ie < ie(m) ) Obc%bound(m)%ied = min(ie(m),ied)

                write(obc_out_unit(m),'(a,10i6)') 'After Open_i ',mpp_pe(), &
                   Obc%bound(m)%is,Obc%bound(m)%ie,js(m),je(m),Obc%bound(m)%np,Obc%bound(m)%isd,Obc%bound(m)%ied,isd,ied

                allocate(Obc%bound(m)%pgrad(Obc%bound(m)%is:Obc%bound(m)%ie))

! I no global domain information for reading boundary data arespecified in the 
! namelist assume, that the data fit exactly at the boundary
                if (iers(m) == -999) then 
                  Obc%bound(m)%iers = is(m) 
                else 
                  Obc%bound(m)%iers = iers(m) 
                endif
                if (iere(m) == -999) then 
                  Obc%bound(m)%iere = ie(m) 
                else 
                  Obc%bound(m)%iere = iere(m) 
                endif 
                if (jers(m) == -999) then 
                  Obc%bound(m)%jers = js(m) 
                else 
                  Obc%bound(m)%jers = jers(m) 
                endif 
                if (jere(m) == -999) then 
                  Obc%bound(m)%jere = je(m) 
                else 
                  Obc%bound(m)%jere = jere(m) 
                endif 
                if (Obc%bound(m)%iers > is(m)) &
                  call mpp_error(FATAL, "ocean_obc_barotrop_mod: start of eta data > start of boundary")
                if (Obc%bound(m)%iere < ie(m)) &
                  call mpp_error(FATAL, "ocean_obc_barotrop_mod: end of eta data < end of boundary")
             endif
          endif

          !--- loop over pelist to figure out the boundary index on each pe.       
          on_bound = .false.
          do pe = 0, npes-1
             if (js(m) .ge. jsl(pe) .and. js(m) .le. jel(pe)) then
                istart(pe) = max(is(m),isl(pe)); iend(pe)  = min(ie(m),iel(pe))
                if(iend(pe) .GE. istart(pe) ) on_bound(pe)    = .true.
             end if
          end do
          nlist = count(on_bound)
          allocate( Obc%bound(m)%start(nlist), obc%bound(m)%end(nlist)    )
          allocate( Obc%bound(m)%pelist(nlist) )
          list = 0
          ntotal = 0
          Obc%bound(m)%index = -1
          sindex = 100000; eindex = -100000
          do n = 0, npes-1
             if(on_bound(n)) then
                list = list+1
                Obc%bound(m)%pelist(list)  = pelist(n)
                if(Obc%bound(m)%on_bound) then
                   Obc%bound(m)%start(list)  = istart(n)
                   Obc%bound(m)%end(list)    = iend(n)
                   sindex                     = min(sindex, istart(n))
                   eindex                     = max(eindex, iend(n))
                   ntotal                     = ntotal + iend(n) - istart(n) + 1
                   if( mpp_pe() == pelist(n)) Obc%bound(m)%index = list
                end if
             end if
          end do
          if(Obc%bound(m)%on_bound) then
             if(Obc%bound(m)%index > -1 .NEQV. obc%bound(m)%on_bound) call mpp_error(FATAL, &
                  "ocean_obc_barotrop_mod: inconsistency for the obc boundary distribution")
             if(Obc%bound(m)%on_bound ) then
                if(Obc%bound(m)%is .NE. Obc%bound(m)%start(Obc%bound(m)%index) .OR. &
                     Obc%bound(m)%ie .NE. Obc%bound(m)%end(Obc%bound(m)%index) ) &
                     call mpp_error(FATAL, "ocean_obc_barotrop_mod: inconsistency for the index setup")
             end if
             if(eindex-sindex+1 .NE. ntotal .OR. Obc%bound(m)%np .NE. ntotal)  &
                  call mpp_error(FATAL, "ocean_obc_barotrop_mod: size mismatch 2")
             allocate(Obc%bound(m)%buffer(sindex:eindex))
          end if
       endif

! Read data for relaxation into the compute area only        
    enddo
    update_eta_tm1 = .false.
    do m=1, Obc%nobc
      if(iand(Obc%bound(m)%bcond_eta, ORLANS) == ORLANS) update_eta_tm1 =.true.
      if(iand(Obc%bound(m)%bcond_eta, MILLER) == MILLER) update_eta_tm1 =.true.
      if(iand(Obc%bound(m)%bcond_eta, CAMOBR) == CAMOBR) update_eta_tm1 =.true.
      if(iand(Obc%bound(m)%bcond_eta, RAYMND) == RAYMND) update_eta_tm1 =.true.
    enddo
    allocate(ctrop(isd:ied,jsd:jed));  ctrop = 0.
    allocate(eta_tm1(isd:ied,jsd:jed));  eta_tm1 = 0.
    allocate(eta_tend_obc(isd:ied,jsd:jed));  eta_tend_obc = 0.
    allocate(rel_coef    (isd:ied,jsd:jed));  rel_coef     = 0.
    allocate(eta_t_tm1   (isd:ied,jsd:jed));  eta_t_tm1    = 0.0

! Now check for corner points, which are located at two boundaries 
! (south-west, north-west, north-east, south-east)    
! This code needs to be generalised later. Currently, only one
! boundary for each direction is possible. In complex coastal geometry
! this may be to restrictive, but should be sufficient for the beginning ... 
    west_at_pe = 0; south_at_pe= 0; east_at_pe = 0; north_at_pe= 0
    do m = 1, nobc
      if (.not. Obc%bound(m)%on_bound) cycle
      if (Obc%bound(m)%direction == WEST)  west_at_pe  = m
      if (Obc%bound(m)%direction == EAST)  east_at_pe  = m
      if (Obc%bound(m)%direction == NORTH) north_at_pe = m
      if (Obc%bound(m)%direction == SOUTH) south_at_pe = m
    enddo
    if (west_at_pe  * south_at_pe .ne. 0) then
      ilef_s = Obc%bound(south_at_pe)%is
      jbou_s = Obc%bound(south_at_pe)%js
      jlow_w  = Obc%bound(west_at_pe)%js
      ibou_w = Obc%bound(west_at_pe)%is
      if (ilef_s .le. ibou_w .and. jlow_w .le. jbou_s) then 
        ! a southern boundary west of a western is nonsense,
        ! however there may be an overlap 
        south_west_corner = .true.
        i_sw = ibou_w
        j_sw = jbou_s
      endif
    endif 
    if (west_at_pe  * north_at_pe .ne. 0) then
      ilef_n = Obc%bound(north_at_pe)%is
      jbou_n = Obc%bound(north_at_pe)%js
      jup_w  = Obc%bound(west_at_pe)%je
      ibou_w = Obc%bound(west_at_pe)%is
      if (ilef_n .le. ibou_w .and. jup_w .ge. jbou_n) then 
        ! a northern boundary west of a western is nonsense,
        ! however there may be an overlap 
        north_west_corner = .true.
        i_nw = ibou_w
        j_nw = jbou_n
      endif
    endif 
    if (east_at_pe  * north_at_pe .ne. 0) then
      irig_n = Obc%bound(north_at_pe)%ie
      jbou_n = Obc%bound(north_at_pe)%js
      jup_e  = Obc%bound(east_at_pe)%je
      ibou_e = Obc%bound(east_at_pe)%is
      if (irig_n .ge. ibou_e .and. jup_e .ge. jbou_n) then 
        ! a northern boundary east of an eastern is nonsense,
        ! however there may be an overlap 
        north_east_corner = .true.
        i_ne = ibou_e
        j_ne = jbou_n
      endif
    endif 
    if (east_at_pe  * south_at_pe .ne. 0) then
      irig_s = Obc%bound(south_at_pe)%ie
      jbou_s = Obc%bound(south_at_pe)%js
      jlow_e  = Obc%bound(east_at_pe)%js
      ibou_e = Obc%bound(east_at_pe)%is
      if (irig_s .ge. ibou_e .and. jlow_e .le. jbou_s) then 
        ! a southern boundary east of an eastern is nonsense,
        ! however there may be an overlap 
        south_east_corner = .true.
        i_se = ibou_e
        j_se = jbou_s
      endif
    endif 

    !------------------------------------------------------
    !------------------------------------------------------
    ! Get the computational domain boundary index maps
    do m = 1, nobc

      if(.not. Obc%bound(m)%on_bound) cycle

      ! Allocate map memory
      ni = Obc%bound(m)%ie - Obc%bound(m)%is
      nj = Obc%bound(m)%je - Obc%bound(m)%js
      nsize = max(ni, nj) + 1
      allocate(Obc%bound(m)%iloc(nsize))
      allocate(Obc%bound(m)%jloc(nsize))
      Obc%bound(m)%sign = 1.0

      !------------------------------------------------------
      ! Western and eastern boundary maps
      if(Obc%bound(m)%direction .eq. WEST .or. &
        Obc%bound(m)%direction .eq. EAST ) then
        ! Interior mapping indicies
        if(Obc%bound(m)%direction .eq. WEST) then
          Obc%bound(m)%imap = 1
          Obc%bound(m)%jmap = 0
        else
          Obc%bound(m)%imap = -1
          Obc%bound(m)%jmap = 0
          Obc%bound(m)%sign = -1.0
        endif
        ! Tangential mapping indicies
        Obc%bound(m)%imapt = 0
        Obc%bound(m)%jmapt = 1

        Obc%bound(m)%nsc = jsc
        Obc%bound(m)%nec = jec

        n = 1
        i =  Obc%bound(m)%is
        allocate(Obc%bound(m)%oi1(nsize))
        allocate(Obc%bound(m)%oi2(nsize))
        allocate(Obc%bound(m)%work(Obc%bound(m)%js:Obc%bound(m)%je))
        if(i .ne. Obc%bound(m)%ie) &
          write(obc_out_unit(m),*) &
           'WARNING : Start and end i locations differ for WEST or EAST obc. is=',i,' : ie=',Obc%bound(m)%ie
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          Obc%bound(m)%iloc(n) = i
          Obc%bound(m)%jloc(n) = j
          Obc%bound(m)%oi1(n) = i + Obc%bound(m)%imap
          Obc%bound(m)%oi2(n) = i + 2 * Obc%bound(m)%imap
          n = n + 1
        enddo
        Obc%bound(m)%nloc = n - 1
        ! Set pointers
        Obc%bound(m)%nic => Obc%bound(m)%iloc
        Obc%bound(m)%tic => Obc%bound(m)%jloc
        Obc%bound(m)%oj1 => Obc%bound(m)%jloc
        Obc%bound(m)%oj2 => Obc%bound(m)%jloc
        Obc%bound(m)%icf => Obc%bound(m)%iloc
        if(Obc%bound(m)%direction==EAST) Obc%bound(m)%icf => Obc%bound(m)%oi1
        Obc%bound(m)%jcf => Obc%bound(m)%jloc
      endif

      !------------------------------------------------------
      ! Northern and southern boundaries
      if(Obc%bound(m)%direction .eq. NORTH .or. &
        Obc%bound(m)%direction .eq. SOUTH ) then
        ! Interior mapping indicies
        if(Obc%bound(m)%direction .eq. SOUTH) then
          Obc%bound(m)%imap = 0
          Obc%bound(m)%jmap = 1
        else
          Obc%bound(m)%imap = 0
          Obc%bound(m)%jmap = -1
          Obc%bound(m)%sign = -1.0
        endif
        ! Tangential mapping indicies
        Obc%bound(m)%imapt = 1
        Obc%bound(m)%jmapt = 0

        Obc%bound(m)%nsc = isc
        Obc%bound(m)%nec = iec

        ! Boundary index maps
        n = 1
        j =  Obc%bound(m)%js
        allocate(Obc%bound(m)%oj1(nsize))
        allocate(Obc%bound(m)%oj2(nsize))
        allocate(Obc%bound(m)%work(Obc%bound(m)%is:Obc%bound(m)%ie))
        if(j .ne. Obc%bound(m)%je) &
          write(obc_out_unit(m),*)        &
             'WARNING : Boundary',n,' start and end j locations differ for NORTH or SOUTH obc. js=' &
             ,j,' : je=',Obc%bound(m)%je
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          Obc%bound(m)%iloc(n) = i
          Obc%bound(m)%jloc(n) = j
          Obc%bound(m)%oj1(n) = j + Obc%bound(m)%jmap
          Obc%bound(m)%oj2(n) = j + 2 * Obc%bound(m)%jmap
          n = n + 1
        enddo
        Obc%bound(m)%nloc = n - 1
        ! Set pointers
        Obc%bound(m)%nic => Obc%bound(m)%jloc
        Obc%bound(m)%tic => Obc%bound(m)%iloc
        Obc%bound(m)%oi1 => Obc%bound(m)%iloc
        Obc%bound(m)%oi2 => Obc%bound(m)%iloc
        Obc%bound(m)%icf => Obc%bound(m)%iloc
        Obc%bound(m)%jcf => Obc%bound(m)%jloc
        if(Obc%bound(m)%direction==NORTH) Obc%bound(m)%jcf => Obc%bound(m)%oj1
        n = 1
      endif

!Needed for advection of tracers
      allocate(Obc%bound(m)%work2(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,nk))

    enddo

    !------------------------------------------------------
    !------------------------------------------------------
    ! Get the data domain boundary index maps
    do m = 1, nobc

      if(.not. Obc%bound(m)%on_bound) cycle

      ! Allocate map memory
      ni = Obc%bound(m)%ied - Obc%bound(m)%isd
      nj = Obc%bound(m)%jed - Obc%bound(m)%jsd
      nsize = max(ni, nj) + 1
      allocate(Obc%bound(m)%ilod(nsize))
      allocate(Obc%bound(m)%jlod(nsize))
      allocate(Obc%bound(m)%dsur(nsize))
      allocate(Obc%bound(m)%dum (nsize))
      allocate(Obc%bound(m)%ones_bt(nsize))
      obc%bound(m)%ones_bt = 1

      !------------------------------------------------------
      ! Western and eastern boundary maps
      if(Obc%bound(m)%direction .eq. WEST .or. &
        Obc%bound(m)%direction .eq. EAST ) then

        Obc%bound(m)%nsd = jsd
        Obc%bound(m)%ned = jed
        n = 1
        i =  Obc%bound(m)%isd
        allocate(Obc%bound(m)%di1(nsize))
        allocate(Obc%bound(m)%di2(nsize))
        allocate(Obc%bound(m)%istr(nsize))
        allocate(Obc%bound(m)%iend(nsize))
        allocate(Obc%bound(m)%ico_bt(nsize))

        if(i .ne. Obc%bound(m)%ied) &
          write(obc_out_unit(m),*) 'WARNING : Start and end i data locations differ for WEST or EAST obc.'
        do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
          Obc%bound(m)%ilod(n) = i
          Obc%bound(m)%jlod(n) = j
          Obc%bound(m)%di1(n) = i + Obc%bound(m)%imap
          Obc%bound(m)%di2(n) = i + 2 * Obc%bound(m)%imap
          if(Obc%bound(m)%direction==WEST) then
             Obc%bound(m)%istr(n) = max(isd,i-2)
             Obc%bound(m)%iend(n) = i - Obc%bound(m)%imap
          else
             Obc%bound(m)%istr(n) = i - Obc%bound(m)%imap
             Obc%bound(m)%iend(n) = min(ied,i+2)
          endif
!kk
          if(Obc%bound(m)%direction==WEST) then
             Obc%bound(m)%ico_bt(n) = i - Obc%bound(m)%imap
          else
             Obc%bound(m)%ico_bt(n) = i
          endif

!kk       csurf on data domain
          if(ht_bt(i,j) <= 0) then
            write(6,*) 'ht_bt(i,j) <= 0 W-O ',mpp_pe(),i,j,ht_bt(i,j)
          else
            Obc%bound(m)%dsur(n)  = sqrt(grav*ht_bt(i,j))
          end if
          n = n + 1
        enddo
        Obc%bound(m)%nlod = n - 1
        ! Set pointers
        Obc%bound(m)%nid => Obc%bound(m)%ilod
        Obc%bound(m)%tid => Obc%bound(m)%jlod
        Obc%bound(m)%dj1 => Obc%bound(m)%jlod
        Obc%bound(m)%dj2 => Obc%bound(m)%jlod
        Obc%bound(m)%idf => Obc%bound(m)%ilod
        if(Obc%bound(m)%direction==EAST) Obc%bound(m)%idf => Obc%bound(m)%di1
        Obc%bound(m)%jdf => Obc%bound(m)%jlod
        Obc%bound(m)%ido => Obc%bound(m)%ilod
        if(Obc%bound(m)%direction==WEST) Obc%bound(m)%ido => Obc%bound(m)%iend
        Obc%bound(m)%jdo => Obc%bound(m)%jlod
        Obc%bound(m)%jstr => Obc%bound(m)%jlod
        Obc%bound(m)%jend => Obc%bound(m)%jlod
        Obc%bound(m)%d1i_bt => Obc%bound(m)%ones_bt
        Obc%bound(m)%d1j_bt => Obc%bound(m)%jlod
        Obc%bound(m)%jco_bt => Obc%bound(m)%jlod
      endif

      !------------------------------------------------------
      ! Northern and southern boundaries
      if(Obc%bound(m)%direction .eq. NORTH .or. &
        Obc%bound(m)%direction .eq. SOUTH ) then

        Obc%bound(m)%nsd = isd
        Obc%bound(m)%ned = ied
        n = 1
        j =  Obc%bound(m)%jsd
        allocate(Obc%bound(m)%dj1(nsize))
        allocate(Obc%bound(m)%dj2(nsize))
        allocate(Obc%bound(m)%jstr(nsize))
        allocate(Obc%bound(m)%jend(nsize))
        allocate(Obc%bound(m)%jco_bt(nsize))
        if(j .ne. Obc%bound(m)%jed) &
          write(obc_out_unit(m),*) 'WARNING : Start and end j data locations differ for NORTH or SOUTH obc.'
        do i = Obc%bound(m)%isd, Obc%bound(m)%ied
          Obc%bound(m)%ilod(n) = i
          Obc%bound(m)%jlod(n) = j
          Obc%bound(m)%dj1(n) = j + Obc%bound(m)%jmap
          Obc%bound(m)%dj2(n) = j + 2 * Obc%bound(m)%jmap
          if(Obc%bound(m)%direction==SOUTH) then
             Obc%bound(m)%jstr(n) = max(jsd,j-2)
             Obc%bound(m)%jend(n) = j - Obc%bound(m)%jmap
          else
             Obc%bound(m)%jstr(n) = j - Obc%bound(m)%jmap
             Obc%bound(m)%jend(n) = min(jed,j+2)
          endif
!kk
          if(Obc%bound(m)%direction==SOUTH) then
             Obc%bound(m)%jco_bt(n) = j - Obc%bound(m)%jmap
          else
             Obc%bound(m)%jco_bt(n) = j
          endif

!kk       csurv on data domain
          if(ht_bt(i,j) <= 0) then
            write(6,*) 'ht_bt(i,j) <= 0 N-S ',mpp_pe(),i,j,ht_bt(i,j)
          else
            Obc%bound(m)%dsur(n)  = sqrt(grav*ht_bt(i,j))
          end if
          n = n + 1
        enddo
        Obc%bound(m)%nlod = n - 1
        ! Set pointers
        Obc%bound(m)%nid => Obc%bound(m)%jlod
        Obc%bound(m)%tid => Obc%bound(m)%ilod
        Obc%bound(m)%di1 => Obc%bound(m)%ilod
        Obc%bound(m)%di2 => Obc%bound(m)%ilod
        Obc%bound(m)%idf => Obc%bound(m)%ilod
        Obc%bound(m)%jdf => Obc%bound(m)%jlod
        if(Obc%bound(m)%direction==NORTH) Obc%bound(m)%jdf => Obc%bound(m)%dj1
        Obc%bound(m)%ido => Obc%bound(m)%ilod
        Obc%bound(m)%jdo => Obc%bound(m)%jlod
        if(Obc%bound(m)%direction==SOUTH) Obc%bound(m)%jdo => Obc%bound(m)%jend
        Obc%bound(m)%istr => Obc%bound(m)%ilod
        Obc%bound(m)%iend => Obc%bound(m)%ilod
        Obc%bound(m)%d1i_bt => Obc%bound(m)%ilod
        Obc%bound(m)%d1j_bt => Obc%bound(m)%ones_bt
        Obc%bound(m)%ico_bt => Obc%bound(m)%ilod
        n = 1
      endif
    enddo

    ! Get restart value for phase speed and relaxation coefficient    
    allocate(Obc%eta   (nobc))  
    do m = 1, nobc
       if (Obc%bound(m)%on_bound) then
          !         data needed for eta
          Obc%eta(m)%rel_coef_in                = rel_coef_eta_in(m)
          Obc%eta(m)%rel_coef_out               = rel_coef_eta_out(m)
          Obc%eta(m)%rel_pnts                   = rel_eta_pnts(m)
          Obc%eta(m)%file                       = trim(filename_eta(m))
          Obc%eta(m)%field                      = trim(fieldname_eta(m))
       endif
    enddo

    allocate(Obc%ud   (nobc))  
    do m = 1, nobc
       if (Obc%bound(m)%on_bound) then
          if(iand(Obc%bound(m)%bcond_ud, FILEIN) == FILEIN) then
             !         data needed for Flather 
             Obc%ud(m)%file                       = trim(filename_ud(m))
             Obc%ud(m)%field                      = trim(fieldname_ud(m))
          endif
       endif
    enddo    

    file_name = 'ocean_obc.res.nc'
    id_restart(1) = register_restart_field(Obc_restart, file_name, 'ctrop', ctrop(:,:), &
            domain=Dom%domain2d)
    id_restart(2) = register_restart_field(Obc_restart, file_name, 'rtrop', rel_coef(:,:), &
            domain=Dom%domain2d)
    id_restart(3) = register_restart_field(Obc_restart, file_name, 'eta_tm1', eta_tm1(:,:), &
            domain=Dom%domain2d)
    id_restart(4) = register_restart_field(Obc_restart, file_name, 'eta_t_tm1', eta_t_tm1(:,:), &
            domain=Dom%domain2d)
    if (file_exist('INPUT/ocean_obc.res.nc')) then
      write (stdoutunit,'(/a)') &
      '  Reading restart for phase speed from INPUT/ocean_obc.res.nc'
      call restore_state(Obc_restart)
      call mpp_update_domains(ctrop(:,:), Dom%domain2d)
      call mpp_update_domains(rel_coef(:,:), Dom%domain2d)
      call mpp_update_domains(eta_tm1(:,:), Dom%domain2d)
      call mpp_update_domains(eta_t_tm1(:,:), Dom%domain2d)
    
    else
      do m = 1, nobc
         if(.not. Obc%bound(m)%on_bound) cycle
      ! initialise the relaxation coefficient and phase speed
      ! For GRAVTY these will be never changed
         call phase_speed_GRAVTY(Obc%bound(m), Obc%eta(m))
      enddo   
    endif
    write(stdoutunit,*) ' '
    write(stdoutunit,*) &
      'From ocean_obc_barotrop_mod: initial obc chksums'
    call write_timestamp(Time%model_time)
    call write_chksum_2d('phase speed', ctrop(COMP))
    call write_chksum_2d('relaxation coefficient', rel_coef(COMP))
    call write_chksum_2d('eta_tm1', eta_tm1(COMP))
    call write_chksum_2d('eta_t_tm1', eta_t_tm1(COMP))
    
    write(stdoutunit,*) '-----------------------------------------------------------------'
    write(stdoutunit,*) 'The following setup for OBC has been found:'
    write(stdoutunit,*) 'Total number of OBC: ', nobc
    do m = 1, nobc
       write(stdoutunit,*) ' Setup of OBC ',m,', ',Obc%bound(m)%name,':'
       write(stdoutunit,*) ' direction: ', trim(direction(m))
       write(stdoutunit,*) ' ==> See files obc_'//trim(Obc%bound(m)%name)//'.out.* for more information.'
       write(stdoutunit,*) '     Enable debug_this_module for details on domain decomposition.'
       
       if(.not. Obc%bound(m)%on_bound) cycle
       
       write(obc_out_unit(m),*) 'Setup of OBC ',m,', ',Obc%bound(m)%name,':'
       write(obc_out_unit(m),*) ' direction     : ', trim(direction(m))
       write(obc_out_unit(m),*) ' points        : ', Obc%bound(m)%np
! if debugging is on each PE has an output file open otherwise this info is useless
       if(debug_this_module) then
         write(obc_out_unit(m),*) ' pe - name       :', trim(pe_name)
         write(obc_out_unit(m),*) ' ocean-domain, comp :', isc, iec, jsc, jec
         write(obc_out_unit(m),*) ' ocean-domain, data :', isd, ied, jsd, jed
         write(obc_out_unit(m),*) ' obc - domain, comp :', Obc%bound(m)%is, Obc%bound(m)%ie, Obc%bound(m)%js, Obc%bound(m)%je
         write(obc_out_unit(m),*) ' obc - domain, data :', Obc%bound(m)%isd, Obc%bound(m)%ied, Obc%bound(m)%jsd, Obc%bound(m)%jed
         if (south_west_corner) write(obc_out_unit(m),*) ' south-west corner points :', i_sw, j_sw
         if (south_east_corner) write(obc_out_unit(m),*) ' south-east corner points :', i_se, j_se
         if (north_east_corner) write(obc_out_unit(m),*) ' north-east corner points :', i_ne, j_ne
         if (north_west_corner) write(obc_out_unit(m),*) ' north-west corner points :', i_nw, j_nw
       endif
       if(.not. (Obc%bound(m)%report .or. debug_this_module)) cycle
       write(obc_out_unit(m),*) ' Normal 3D velocity OBC : ', trim(obc_nor(m))
       write(obc_out_unit(m),*) ' Tangential 3D velocity OBC : ', trim(obc_tan(m))
       write(obc_out_unit(m),*) ' Elevation OBC : ', trim(obc_eta(m))
       write(obc_out_unit(m),*) ' Normal 2D velocity OBC : ', trim(obc_ud(m))
       write(obc_out_unit(m),*) ' Vertical mixing coefficients OBC : ', trim(obc_mix(m))
       write(obc_out_unit(m),*) ' min phase speed (times sqrt(gH)) : ', Obc%bound(m)%ctrop_min
       write(obc_out_unit(m),*) ' max phase speed (times sqrt(gH)) : ', Obc%bound(m)%ctrop_max
       write(obc_out_unit(m),*) ' inc phase speed (times sqrt(gH)) : ', Obc%bound(m)%ctrop_inc
       write(obc_out_unit(m),*) ' smoothing parameter              : ', Obc%bound(m)%ctrop_smooth
       if (Obc%bound(m)%obc_damp_newton) then
          write(obc_out_unit(m),*) ' apply Newtonian damping of ',Obc%bound(m)%damp_factor
       else
          write(obc_out_unit(m),*) ' no Newtonian damping '
       endif
       if (Obc%bound(m)%obc_consider_convu) then 
          write(obc_out_unit(m),*) ' consider convu for d eta/dt'
       else  
          write(obc_out_unit(m),*) ' do not consider convu for d eta/dt'
       endif
       write(obc_out_unit(m),*) ' relax eta     :', Obc%bound(m)%obc_relax_eta
       if (Obc%bound(m)%obc_relax_eta_profile) then
          write(obc_out_unit(m),*) '  method        :', ' profile'
       else  
          write(obc_out_unit(m),*) '  method        :', ' average'
       endif
       write(obc_out_unit(m),*) '-----------------------------------------------------------------'
    enddo


    do m = 1, nobc
       call mpp_declare_pelist(Obc%bound(m)%pelist)
    enddo

    do m = 1, nobc
       if (Obc%bound(m)%on_bound) then
! This part of the domain definition is the same for tracers and eta
          call mpp_set_compute_domain(obc_eta_domain(m),&
                                      Obc%bound(m)%isd, Obc%bound(m)%ied, &
                                      Obc%bound(m)%jsd, Obc%bound(m)%jed, &
                                      Obc%bound(m)%ied-Obc%bound(m)%isd+1, Obc%bound(m)%jed-Obc%bound(m)%jsd+1)
          call mpp_set_data_domain   (obc_eta_domain(m),&
                                      Obc%bound(m)%isd, Obc%bound(m)%ied, &
                                      Obc%bound(m)%jsd, Obc%bound(m)%jed, &
                                      Obc%bound(m)%ied-Obc%bound(m)%isd+1, Obc%bound(m)%jed-Obc%bound(m)%jsd+1,&
                                       .false., .false.)
          call mpp_set_compute_domain(obc_diag_domain(m),&
                                      Obc%bound(m)%is, Obc%bound(m)%ie, &
                                      Obc%bound(m)%js, Obc%bound(m)%je, &
                                      Obc%bound(m)%ie-Obc%bound(m)%is+1, Obc%bound(m)%je-Obc%bound(m)%js+1)
          call mpp_set_data_domain   (obc_diag_domain(m),&
                                      Obc%bound(m)%is, Obc%bound(m)%ie, &
                                      Obc%bound(m)%js, Obc%bound(m)%je, &
                                      Obc%bound(m)%ie-Obc%bound(m)%is+1, Obc%bound(m)%je-Obc%bound(m)%js+1,&
                                       .false., .false.)
          if(Obc%bound(m)%obc_relax_eta) then
             if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
!kk               Changed from compute to data domain
                allocate(Obc%eta(m)%data(1,Obc%bound(m)%jsd:Obc%bound(m)%jed))
             else
                allocate(Obc%eta(m)%data(Obc%bound(m)%isd:Obc%bound(m)%ied,1))
             endif
             Obc%eta(m)%data = 0
             write(obc_out_unit(m),*) trim(Obc%eta(m)%file)
             write(obc_out_unit(m),*) trim(Obc%eta(m)%field)

             ! Set the global domain part for the eta input file
             call mpp_set_global_domain( obc_eta_domain(m), &
                                      xbegin=Obc%bound(m)%iers, ybegin=Obc%bound(m)%jers, &
                                      xend=Obc%bound(m)%iere, yend=Obc%bound(m)%jere, &
                                      xsize=Obc%bound(m)%iere-Obc%bound(m)%iers+1, &
                                      ysize=Obc%bound(m)%jere-Obc%bound(m)%jers+1)

             Obc%eta(m)%id   = init_external_field( &
                  Obc%eta(m)%file,Obc%eta(m)%field, &
                  domain=obc_eta_domain(m),        &
!                  desired_units='m', &
                  verbose=debug_this_module)
             fld_size(:)  = get_external_field_size(Obc%eta(m)%id)
             if(debug_this_module) write(obc_out_unit(m),*) 'fld_size for this PE',fld_size(:)
             if(fld_size(2) .ne. size(Obc%eta(m)%data,2)) then
                write(errorstring,'(I3,I3,1x)') fld_size(2), size(Obc%eta(m)%data,2)
                call mpp_error(FATAL,trim(errorstring)//'invalid dimension 2 of input field in '//trim(Obc%eta(m)%file))
             endif
         endif

         if(iand(Obc%bound(m)%bcond_ud, FILEIN) == FILEIN) then
              if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
!kk               Changed from compute to data domain
                  allocate(Obc%ud(m)%data(1,Obc%bound(m)%jsd:Obc%bound(m)%jed))

                  write(*,*) ' e/w pe,js,je=',mpp_pe(),Obc%bound(m)%js,Obc%bound(m)%je
              else
                  allocate(Obc%ud(m)%data(Obc%bound(m)%isd:Obc%bound(m)%ied,1))

                  write(*,*) ' n/s pe,js,je=',mpp_pe(),Obc%bound(m)%is,Obc%bound(m)%ie
              endif
              Obc%ud(m)%data = 0
              write(stdoutunit,*) trim(Obc%ud(m)%file)
              write(stdoutunit,*) trim(Obc%ud(m)%field)
! ! Set the global domain part for the eta input file
              call mpp_set_global_domain( obc_eta_domain(m), &
                                       xbegin=Obc%bound(m)%iers, ybegin=Obc%bound(m)%jers, &
                                       xend=Obc%bound(m)%iere, yend=Obc%bound(m)%jere, &
                                       xsize=Obc%bound(m)%iere-Obc%bound(m)%iers+1, &
                                       ysize=Obc%bound(m)%jere-Obc%bound(m)%jers+1)

              Obc%ud(m)%id   = init_external_field( &
                   Obc%ud(m)%file,Obc%ud(m)%field, &
                   domain=obc_eta_domain(m))

              fld_size(:)  = get_external_field_size(Obc%ud(m)%id)
              if(debug_this_module) write(obc_out_unit(m),*) 'fld_size for this PE',fld_size(:)
              if(fld_size(2) .ne. size(Obc%ud(m)%data,2)) then
                 write(errorstring,'(I3,I3)') fld_size(2), size(Obc%ud(m)%data,2)
                 call mpp_error(FATAL,trim(errorstring)//'invalid dimension 2 of input field in '//trim(Obc%ud(m)%file))
              endif
           endif

          write(obc_out_unit(m),*) '-----------------------------------------------------------------'
          write(obc_out_unit(m),*) 'The following scheme for OBC data input has been found:'
          if (Obc%bound(m)%obc_relax_eta) then
             write(obc_out_unit(m),*) '  sea level relaxation'
             write(obc_out_unit(m),*) '  rel_coef_in : ', Obc%eta(m)%rel_coef_in
             write(obc_out_unit(m),*) '  rel_coef_out: ', Obc%eta(m)%rel_coef_out
             write(obc_out_unit(m),*) '  rel_pnts    : ', Obc%eta(m)%rel_pnts
             write(obc_out_unit(m),*) '  input-file  : ', trim(Obc%eta(m)%file)
             write(obc_out_unit(m),*) '  input-name  : ', trim(Obc%eta(m)%field)
          else
             write(obc_out_unit(m),*) '  no sea level relaxation'
          endif
       endif
    enddo


! prepare diagnostics
    do m = 1, nobc
       if (Obc%bound(m)%on_bound) then
         call mpp_set_global_domain( obc_diag_domain(m), &
                                   xbegin=Obc%bound(m)%isg, ybegin=Obc%bound(m)%jsg, &
                                   xend  =Obc%bound(m)%ieg, yend  =Obc%bound(m)%jeg, &
                                   xsize =Obc%bound(m)%ieg - Obc%bound(m)%isg+1, &
                                   ysize =Obc%bound(m)%jeg - Obc%bound(m)%jsg+1)
         Obc%bound(m)%id_xt = diag_axis_init ('xt_'//trim(Obc%bound(m)%name),Grd%grid_x_t(Obc%bound(m)%isg:Obc%bound(m)%ieg), &
              'degrees_E','x','tcell longitude',set_name='ocean', Domain2=obc_diag_domain(m), aux='geolon_t')
         Obc%bound(m)%id_yt = diag_axis_init ('yt_'//trim(Obc%bound(m)%name),Grd%grid_y_t(Obc%bound(m)%jsg:Obc%bound(m)%jeg), &
              'degrees_N','y','tcell latitude',set_name='ocean', Domain2=obc_diag_domain(m), aux='geolat_t')
         
         Obc%bound(m)%id_ctrop = register_diag_field ('ocean_model', 'ctrop_p_'//trim(Obc%bound(m)%name),          &
                (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yt/), Time%model_time, 'barotr phase speed on open bounds',  &
                                 'm/s', missing_value=missing_value, range=(/-1.e8,1.e8/))
         Obc%bound(m)%id_eta_data = register_diag_field ('ocean_model', 'eta_data_'//trim(Obc%bound(m)%name),                &
                       (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yt/), Time%model_time, 'external sea level data on open bounds',&
                                 'm', missing_value=missing_value, range=(/-1.e8,1.e8/))
         Obc%bound(m)%id_rel_coef = register_diag_field ('ocean_model', 'rel_coeff_'//trim(Obc%bound(m)%name),           &
                                 (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yt/), Time%model_time, 'relaxation coefficient',  &
                                 'none', missing_value=missing_value, range=(/-1.e8,1.e8/))

       endif
    enddo
    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_barotrop_init
  !  </SUBROUTINE> NAME="ocean_obc_barotrop_init"


  !<FUNCTION> NAME="ocean_obc_check_for_update"
  !<PUBLICROUTINE>
  function ocean_obc_check_for_update()
  !</PUBLICROUTINE>
    logical    :: ocean_obc_check_for_update
    logical    :: do_update
    integer    :: m
    do_update=.false.
    do m = 1, nobc
      if(iand(Obc%bound(m)%bcond_eta, RAYMND) == RAYMND) do_update=.true.
!      if(iand(Obc%bound(m)%bcond_eta, ORLANS) == ORLANS) do_update=.true.
    enddo
    ocean_obc_check_for_update=do_update
  end function ocean_obc_check_for_update
  !  </FUNCTION> NAME="ocean_obc_check_for_update"


  !#####################################################################
  !  <SUBROUTINE NAME="ocean_obc_prepare">
  !   <DESCRIPTION>
  !      Prepares OBC  
  !      
  !   </DESCRIPTION>
  !<PUBLICROUTINE>
  subroutine ocean_obc_prepare(Time)
  !</PUBLICROUTINE>
    type(ocean_time_type), intent(in)             :: Time 

    integer                                       :: m, n, i, j, taum1, tau
    integer                                       :: nn, id, jd
    logical                                       :: used
    integer                                       :: istrt, iend, jstrt, jend, ii, jj
    real                                          :: cdtdxr, cdtdyr, etasum, dummy_array(1)
    real, dimension(:), pointer :: data
    real, dimension(isd:ied,jsd:jed,nk)           :: wrk

!    call mpp_clock_begin(id_obc)

!   prepare the new data needed for time interpolation of external data.
    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

!      Prepare sea level data for relaxation         
       if (Obc%bound(m)%obc_relax_eta) then
         

         call time_interp_external(Obc%eta(m)%id,Time%model_time,Obc%eta(m)%data,verbose=debug_this_module)
       endif

!      Prepare Flather
       if(iand(Obc%bound(m)%bcond_ud, FILEIN) == FILEIN) then
         call time_interp_external(Obc%ud(m)%id,Time%model_time,Obc%ud(m)%data,verbose=debug_this_module)
       endif     

    enddo

    ! Relax eta_t(taum1) towards the previous external data
    ! This is brought in the scheme like a smoother
    ! Data are added this way to fields
    ! store eta_t before obc apply for consistent tracer treatment
    ! relaxation coefficxient is from previous barotropic sequence,
    ! data are old from previous time step. They are updated later below.
    ! This subroutine should be called after smoothing eta_t
    
    taum1=Time%taum1
    eta_tend_obc = 0.
    
      do m = 1, nobc
        !--- if on current pe there is no point on the bound, then just return
        if(.not. Obc%bound(m)%on_bound ) cycle
        if(.not. Obc%bound(m)%obc_relax_eta) cycle
        if(Obc%bound(m)%id_eta_data > 0) then
          wrk = 0.
          if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
            do j = Obc%bound(m)%js, Obc%bound(m)%je
              wrk(Obc%bound(m)%is,j,1) = Obc%eta(m)%data(1,j)
            enddo
          endif
          if(Obc%bound(m)%direction == SOUTH .or. Obc%bound(m)%direction == NORTH) then
            do i = Obc%bound(m)%is, Obc%bound(m)%ie
              wrk(i,Obc%bound(m)%js,1) = Obc%eta(m)%data(i,1)
            enddo
          endif
          used = send_data(Obc%bound(m)%id_eta_data, &
               wrk(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1),  &
               Time%model_time,                               &
               rmask=t_mask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je))
        endif
      enddo
!      call mpp_clock_end(id_obc)

    return

  end subroutine ocean_obc_prepare
  !  </SUBROUTINE> NAME="ocean_obc_prepare"


  !#####################################################################
  !set the divergency of the vertically integrated velocity to zero if needed
  !<SUBROUTINE NAME="ocean_obc_adjust_divud" >
  !   <INOUT NAME="divud" TYPE="real, dimension(isd:,jsd:)"></INOUT>

  !<PUBLICROUTINE INTERFACE="ocean_obc_adjust_divud" >
  subroutine ocean_obc_adjust_divud(divud)
  !</PUBLICROUTINE>
    real, dimension(isd:,jsd:),   intent(inout) :: divud
    integer :: m
    
    call mpp_clock_begin(id_obc)
    
    do m = 1, nobc
       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle
       if(Obc%bound(m)%obc_consider_convu) then
          call ocean_obc_update_boundary_2d(divud, 'T')
       else
          call ocean_obc_zero_boundary_2d(divud, 'T')
          call ocean_obc_update_boundary_2d(divud, 'T')
          if(debug_this_module) then
             write(obc_out_unit(m),*) 'Setting div_ud to zero for OBC ', trim(Obc%bound(m)%name)
          endif
       endif

    enddo
    call mpp_clock_end(id_obc)

  end subroutine ocean_obc_adjust_divud
  !  </SUBROUTINE> NAME="ocean_obc_adjust_divud"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_damp_newton" >

  !<PUBLICROUTINE INTERFACE="ocean_obc_damp_newton" >
  subroutine ocean_obc_damp_newton(udrho_bt,forcing)
  !</PUBLICROUTINE>
    real, dimension(isd:ied,jsd:jed,2),   intent(in) :: udrho_bt
    real, dimension(isd:ied,jsd:jed,2),   intent(inout) :: forcing
    character(len=1)    :: grid_type
    integer :: m, n, nn, i, j
    
    call mpp_clock_begin(id_obc)
    
    do m = 1, nobc
       !--- if on current pe there is no point on the bound, then just return
      if(.not. Obc%bound(m)%on_bound) cycle
      if(.not. Obc%bound(m)%obc_damp_newton) cycle
      
      
      if( Obc%bound(m)%direction == EAST.or. &
           Obc%bound(m)%direction == WEST) then
          
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%idf(n)   ! i boundary location
             j = Obc%bound(m)%jdf(n)   ! j boundary location
             nn = Obc%bound(m)%tid(n)  ! Index tangential to boundary
             
             if(Obc%bound(m)%mask(nn)) forcing(i,j,2) = -1*udrho_bt(i,j,2) * &
                  Obc%bound(m)%damp_factor

          enddo

      else if (Obc%bound(m)%direction == NORTH.or. &
           Obc%bound(m)%direction == SOUTH) then

          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%idf(n)   ! i boundary location
             j = Obc%bound(m)%jdf(n)   ! j boundary location
             nn = Obc%bound(m)%tid(n)  ! Index tangential to boundary

             if(Obc%bound(m)%mask(nn)) forcing(i,j,1) = -1*udrho_bt(i,j,1) * &
                  Obc%bound(m)%damp_factor


          enddo

      endif

      if(debug_this_module) then
         write(obc_out_unit(m),*) 'Applying newtonian damping for OBC ', trim(Obc%bound(m)%name)
      endif

   enddo
   call mpp_clock_end(id_obc)

  end subroutine ocean_obc_damp_newton
  !  </SUBROUTINE> NAME="ocean_obc_damp_newton"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_ud" >

    !<PUBLICROUTINE INTERFACE="ocean_obc_ud" >
  subroutine ocean_obc_ud(eta_t, udrho)
    !</PUBLICROUTINE>
    real, dimension(isd:ied,jsd:jed),   intent(in) :: eta_t
    real, dimension(isd:ied,jsd:jed,2),   intent(inout) :: udrho
!    real                                              :: dt
!    real, dimension(isd:ied,jsd:jed,2),   intent(inout) :: forcing
    real :: ud_bound, eta_diff, tendency
    character(len=1)    :: grid_type
    integer :: m, n, nn, i, j, idf, jdf, ico, jco, vn, dd
    logical :: flather_on=.false. 
    real, dimension(:), pointer :: hdata
    real, dimension(:), pointer :: udata
    integer, dimension(:), pointer :: dn

    call mpp_clock_begin(id_obc)

    do m = 1, nobc
       !--- if on current pe there is no point on the bound, then just return
      if(.not. Obc%bound(m)%on_bound) cycle
      if(.not. Obc%bound(m)%obc_ud_active) cycle

      ! Set pointers and initialize
      if (Obc%bound(m)%direction == NORTH.or.Obc%bound(m)%direction == SOUTH) then
         vn = 2
         dn => Obc%bound(m)%ilod
      else
         vn = 1
         dn => Obc%bound(m)%jlod
      endif

      ! Initialise the outside elevation location. This is overwritten with data (below) or a NOGRAD
      ! condition (in ocean_obc_update_boundary() below) later.
      do n = 1, Obc%bound(m)%nlod
         ico = Obc%bound(m)%ico_bt(n)      ! i boundary outside face location
         jco = Obc%bound(m)%jco_bt(n)      ! j boundary outside face location
         nn = Obc%bound(m)%tid(n)          ! Index tangential to boundary
         if(Obc%bound(m)%mask(nn)) then
            udrho(ico,jco,vn) = 0.0
         endif
      enddo

      ! Set depth integrated velocity to FILEIN value. This assignment occurs on the outside elevation (i,j) node.
      ! Set bcond_eta = NOTHIN if this condition is used in isolation.
      if(iand(Obc%bound(m)%bcond_ud, FILEIN) == FILEIN) then
         if (Obc%bound(m)%direction == NORTH.or.Obc%bound(m)%direction == SOUTH) then
            udata => Obc%ud(m)%data(:,1)
         else
            udata => Obc%ud(m)%data(1,:)
         endif
         do n = 1, Obc%bound(m)%nlod
            ico = Obc%bound(m)%ico_bt(n)   ! i boundary outside face location
            jco = Obc%bound(m)%jco_bt(n)   ! j boundary outside face location
            nn = Obc%bound(m)%tid(n)       ! Index tangential to boundary
            dd = dn(n)                     ! Data index
            if(Obc%bound(m)%mask(nn)) then
               udrho(ico,jco,vn) = udata(dd)
            endif
         enddo
      endif

      ! Set depth integrated velocity to Flather condition. This is implemented 
      ! on the inside elevation node. The elevation external data is used if bcond_eta |= FLATHR.
      if(iand(Obc%bound(m)%bcond_ud, FLATHR) == FLATHR) then

         Obc%bound(m)%work = 0.0
         hdata => Obc%bound(m)%work
         if (Obc%bound(m)%direction == NORTH.or.Obc%bound(m)%direction == SOUTH) then
            if(iand(Obc%bound(m)%bcond_eta, FLATHR) == FLATHR) hdata => Obc%eta(m)%data(:,1)
         else
            if(iand(Obc%bound(m)%bcond_eta, FLATHR) == FLATHR) hdata => Obc%eta(m)%data(1,:)
         endif
         flather_on = .true. 

         do n = 1, Obc%bound(m)%nlod
            i = Obc%bound(m)%ilod(n)       ! i boundary location
            j = Obc%bound(m)%jlod(n)       ! j boundary location
            idf = Obc%bound(m)%idf(n)      ! i boundary inside face location
            jdf = Obc%bound(m)%jdf(n)      ! j boundary inside face location
            ico = Obc%bound(m)%ico_bt(n)      ! i boundary outside face location
            jco = Obc%bound(m)%jco_bt(n)      ! i boundary outside face location
            nn = Obc%bound(m)%tid(n)       ! Index tangential to boundary
            dd = dn(n)                     ! Data index
               
            if(Obc%bound(m)%mask(nn)) then
               eta_diff = eta_t(i,j)  - hdata(dd)
               ud_bound = udrho(ico,jco,vn) - Obc%bound(m)%sign * Obc%bound(m)%dsur(n) * eta_diff
               udrho(idf,jdf,vn) = ud_bound * rho0
            endif
         enddo
      endif
   enddo

   if(flather_on) then 
!      call mpp_update_domains (udrho(:,:,1),udrho(:,:,2), Dom%domain2d, gridtype=BGRID_NE)     
! guess this can be removed when working at the data domain or the large barotropic domain respectively???
   endif
   
   ! Set the depth integrated velocity on the outside elevation node, and tangential node.
   ! The normal depth averaged velocity currently uses the bcond_nor OBC.
   call ocean_obc_update_boundary(udrho(:,:,1),'M','u')
   call ocean_obc_update_boundary(udrho(:,:,2),'M','t')
   call ocean_obc_update_boundary(udrho(:,:,1),'Z','t')
   call ocean_obc_update_boundary(udrho(:,:,2),'Z','u')

   call mpp_clock_end(id_obc)

 end subroutine ocean_obc_ud
!  </SUBROUTINE> NAME="ocean_obc_ud"


  !#####################################################################
  !update surface height and the vertically integrated horizontal velocity
  !<SUBROUTINE NAME="ocean_obc_barotropic">
  !   <IN NAME="taum1, tau, taup1" TYPE="integer"> </IN>
  !   <INOUT NAME="eta" TYPE="real, dimension(isd:,jsd:,:)"></INOUT>

  !<PUBLICROUTINE INTERFACE="ocean_obc_barotropic" >
  subroutine ocean_obc_barotropic(eta, taum1, tau, taup1, tstep)
  !</PUBLICROUTINE>

    real, dimension(isd:,jsd:,:), target, intent(inout) :: eta
    integer, intent(in)                           :: taum1, tau, taup1
    real, intent(in)                              :: tstep

    real    :: cdtdxr, dt, etasum, dummy_array(1)
    integer :: i, j, m, n, nn, i1, j1, id, jd
    integer :: istrt, iend, jstrt, jend, ii, jj
    real, dimension(:), pointer :: data
    real, dimension(:,:), pointer :: dg
    real :: sign

  integer :: stdoutunit
  stdoutunit=stdout()

    call mpp_clock_begin(id_obc)

    !   Update from taup1 to taup1. This implies, that divuv must be set
    ! to zero at boundaries (obc_consider_convu=.false.), before regular
    ! points are updated.
    !   With obc_consider_convu=.true. the divergency
    !   in the boundaries direction is taken into account. Also fresh water
    ! flux is already taken into account with the regular scheme.
    

    do m = 1, nobc

       ! if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

       ! Set pointers
       if(Obc%bound(m)%direction == WEST .or. &
          Obc%bound(m)%direction == EAST) then
          dg => dxu_bt
       else
          dg => dyu_bt
       endif
       sign = Obc%bound(m)%sign
       dt = tstep

       ! Calculate the phase speed
       if(iand(Obc%bound(m)%bcond_eta, IOW) == IOW) then
          ! Original momp1 Orlanski like formulation by Martin Schmidt.
          call phase_speed_IOW(Obc%bound(m), Obc%eta(m), eta, taum1, taup1, dt)
          ! Get the sea level due to radiation
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)   ! i boundary location
             j = Obc%bound(m)%jlod(n)   ! j boundary location
             nn = Obc%bound(m)%tid(n)   ! Index tangential to boundary
             i1 = i  + Obc%bound(m)%imap
             j1 = j  + Obc%bound(m)%jmap

             cdtdxr = ctrop(i,j) * dt / dg(i,j)
             eta(i,j,taup1) = (eta(i,j,taup1) + cdtdxr * eta(i1,j1,taup1)) / &
                              (1. + cdtdxr)
          enddo
       else if(iand(Obc%bound(m)%bcond_eta, ORLANS) == ORLANS) then
          ! Orlanski (1976) implicit.
          ! Based on Chapman (1985) Table 1. eqn. 8 (ORI)
          call phase_speed_ORLANS(Obc%bound(m), Obc%eta(m), eta, tau, taup1, dt)
          ! Get the sea level due to radiation
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)   ! i boundary location
             j = Obc%bound(m)%jlod(n)   ! j boundary location
             i1 = i  + Obc%bound(m)%imap   ! i interior cell location
             j1 = j  + Obc%bound(m)%jmap   ! j interior cell location

             cdtdxr = ctrop(i,j) * dt / dg(i,j)
             eta(i,j,taup1) = (eta_tm1(i,j) * (1.0 - cdtdxr) + &
                  2.0 * cdtdxr * eta(i1,j1,tau)) / (1.0 + cdtdxr)
             ! Save the value of eta at the backward time step
             eta_tm1(i,j)   = eta(i,j,tau)
             eta_tm1(i1,j1) = eta(i1,j1,tau)
          enddo
          dt = dt * 2.0
       else if(iand(Obc%bound(m)%bcond_eta, GRAVTY) == GRAVTY) then
          ! Gravity wave radiation implicit.
          ! Based on Chapman (1985) Table 1. eqn. 4 (GWI)
          ! Gravity wave speed is a constant and is calculated in
          ! ocean_obc_barotrop_init. It is stored in ctrop.
          ! Hence values in ctrop are valid.
          ! Get the sea level due to radiation.
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)   ! i boundary location
             j = Obc%bound(m)%jlod(n)   ! j boundary location
             i1 = i  + Obc%bound(m)%imap   ! i interior cell location
             j1 = j  + Obc%bound(m)%jmap   ! j interior cell location

             cdtdxr = ctrop(i,j) * dt / dg(i,j)
             eta(i,j,taup1) = (eta(i,j,tau) + cdtdxr * eta(i1,j1,taup1)) / &
                              (1.0 + cdtdxr)
          enddo
       else if(iand(Obc%bound(m)%bcond_eta, CAMOBR) == CAMOBR) then
          ! Camerlengo & O'Brien (1980) implicit.
          ! Based on Chapman (1985) Table 1. eqn. 10 (MOI)
          call phase_speed_ORLANS(Obc%bound(m), Obc%eta(m), eta, tau, taup1, dt)
          ! Get the sea level due to radiation
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)   ! i boundary location
             j = Obc%bound(m)%jlod(n)   ! j boundary location
             i1 = i  + Obc%bound(m)%imap   ! i interior cell location
             j1 = j  + Obc%bound(m)%jmap   ! j interior cell location

             if (ctrop(i,j) .gt. 0.0) then
                eta(i,j,taup1) = eta(i1,j1,tau)
             else
                eta(i,j,taup1) = eta_tm1(i,j)
             endif
             ! Save the value of eta at the backward time step
             eta_tm1(i,j)   = eta(i,j,tau)
             eta_tm1(i1,j1) = eta(i1,j1,tau)
          enddo
          dt = dt * 2.0
       else if(iand(Obc%bound(m)%bcond_eta, MILLER) == MILLER) then
          ! Miller and Thorpe (1980) explicit.
          ! Based on Miller and Thorpe (1981) eqn. 15
          call phase_speed_MILLER(Obc%bound(m), Obc%eta(m), eta, tau, taup1, dt)
          ! Get the sea level due to radiation
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)   ! i boundary location
             j = Obc%bound(m)%jlod(n)   ! j boundary location
             i1 = i  + Obc%bound(m)%imap   ! i interior cell location
             j1 = j  + Obc%bound(m)%jmap   ! j interior cell location
             id = i  + 2*Obc%bound(m)%imap   ! 2xi interior cell location
             jd = j  + 2*Obc%bound(m)%jmap   ! 2xj interior cell location

             cdtdxr = ctrop(i,j) * dt / dg(i,j)
             eta(i,j,taup1) = eta(i,j,tau) - cdtdxr * &
                  (eta(i,j,tau) - eta(i1,j1,tau))
             ! Save the value of eta at the backward time step
             eta_tm1(i,j)   = eta(i,j,tau)
             eta_tm1(i1,j1) = eta(i1,j1,tau)
             eta_tm1(id,jd) = eta(id,jd,tau)
          enddo
       else if(iand(Obc%bound(m)%bcond_eta, RAYMND) == RAYMND) then
          ! Raymond and Kuo (1984) implicit.
          call phase_speed_RAYMND(Obc%bound(m), Obc%eta(m), eta, tau, taup1, dt)
          ! Get the sea level due to radiation
          do n = 1, Obc%bound(m)%nloc
             i = Obc%bound(m)%iloc(n)   ! i boundary location
             j = Obc%bound(m)%jloc(n)   ! j boundary location
             i1 = i  + Obc%bound(m)%imap   ! i interior cell location
             j1 = j  + Obc%bound(m)%jmap   ! j interior cell location
             id = i  + 2*Obc%bound(m)%imap   ! 2xi interior cell location
             jd = j  + 2*Obc%bound(m)%jmap   ! 2xj interior cell location

!  [m s-1][s m-1] = [nd]
             cdtdxr = ctrop(i,j) * dt / dg(i,j)
             eta(i,j,taup1) = (eta(i,j,tau) +  &
                  cdtdxr * eta(i1,j1,taup1) - Obc%bound(m)%dum(n)) / &
                  (1 + cdtdxr)
             ! Save the value of eta at the backward time step
             eta_tm1(i,j)   = eta(i,j,tau)
             eta_tm1(i1,j1) = eta(i1,j1,tau)
          enddo
       else if(iand(Obc%bound(m)%bcond_eta, NOGRAD) == NOGRAD) then
          ! Gravity wave radiation implicit.
          ! Based on Chapman (1985) Table 1. eqn. 4 (GWI)
          ! Gravity wave speed is a constant and is calculated in
          ! ocean_obc_barotrop_init. It is stored in ctrop.
          ! Hence values in ctrop are valid.
          ! Get the sea level due to radiation.
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)   ! i boundary location
             j = Obc%bound(m)%jlod(n)   ! j boundary location
             i1 = i  + Obc%bound(m)%imap   ! i interior cell location
             j1 = j  + Obc%bound(m)%jmap   ! j interior cell location

             eta(i,j,taup1) = eta(i1,j1,taup1)
          enddo
       endif

       ! Get the sea level relaxation component
       if(iand(Obc%bound(m)%bcond_eta, FILEIN) == FILEIN) then
          if(Obc%bound(m)%bcond_eta == FILEIN) then
             ! Force with the prescribed profile only
             do n = 1, Obc%bound(m)%nlod
                i = Obc%bound(m)%ilod(n)
                j = Obc%bound(m)%jlod(n)
                id = Obc%bound(m)%d1i_bt(n)
                jd = Obc%bound(m)%d1j_bt(n)

                ii = i
                jj = j
!kk                write(0,'(a,8i6,E10.2)') "kk_FILEIN ",n,i,j,id,jd,ii,jj,Obc%eta(m)%rel_pnts,Obc%eta(m)%data(id,jd)
                do nn = 1, Obc%eta(m)%rel_pnts
                   eta(ii,jj,taup1) = Obc%eta(m)%data(id,jd)
                   ii = ii + Obc%bound(m)%imap
                   jj = jj + Obc%bound(m)%jmap
                enddo
             enddo
          else
          ! Relax implicitly towards a prescribed profile
             do n = 1, Obc%bound(m)%nlod
                i = Obc%bound(m)%ilod(n)
                j = Obc%bound(m)%jlod(n)
                id = Obc%bound(m)%d1i_bt(n)
                jd = Obc%bound(m)%d1j_bt(n)

                cdtdxr = rel_coef(i,j) * dt
                ii = i
                jj = j
                do nn = 1, Obc%eta(m)%rel_pnts
                   eta(ii,jj,taup1) = (eta(ii,jj,taup1) + &
                        cdtdxr * Obc%eta(m)%data(id,jd)) / (1. + cdtdxr)
                   ii = ii + Obc%bound(m)%imap
                   jj = jj + Obc%bound(m)%jmap
                enddo
             enddo
          endif
       endif

       if(iand(Obc%bound(m)%bcond_eta, MEANIN) == MEANIN) then
          ! Relax mean sea level leaving the gradient profile in the
          ! boundary unchanged.
          if(Obc%bound(m)%on_bound) then
             n = Obc%bound(m)%nloc
             if(Obc%bound(m)%direction == WEST .or. &
                Obc%bound(m)%direction == EAST) then
                i = Obc%bound(m)%oi1(1)
                data => eta(i, Obc%bound(m)%jloc(1):Obc%bound(m)%jloc(n), &
                            taup1)
             else
                j = Obc%bound(m)%oj1(1)
                data => eta(Obc%bound(m)%iloc(1):Obc%bound(m)%iloc(n), j, &
                            taup1)
             endif
             etasum = boundary_average(obc%bound(m),data)
          else
             etasum = boundary_average(obc%bound(m), dummy_array)
          end if
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%ilod(n)
             j = Obc%bound(m)%jlod(n)
             id = Obc%bound(m)%d1i_bt(n)
             jd = Obc%bound(m)%d1j_bt(n)
             cdtdxr = rel_coef(i,j) * dt
             ii = i
             jj = j
             do nn = 1, Obc%eta(m)%rel_pnts
                eta(ii,j,taup1) = eta(ii,j,taup1) + &
                                     cdtdxr*(Obc%eta(m)%data(id,jd) - etasum)
                ii = ii + Obc%bound(m)%imap
                jj = jj + Obc%bound(m)%jmap
             enddo
          enddo
       endif
       ! Phase speed debugging information
       if(debug_phase_speed) then
          n = Obc%bound(m)%nlod/2     ! Print for a location halfway along the boundary
          i = Obc%bound(m)%ilod(n)
          j = Obc%bound(m)%jlod(n)
          i1 = i  + Obc%bound(m)%imap
          j1 = j  + Obc%bound(m)%jmap
          if(Obc%bound(m)%direction == WEST .or. &
             Obc%bound(m)%direction == EAST) then
             cdtdxr = dxu_bt(i1,j1)/dt
          else
             cdtdxr = dyu_bt(i1,j1)/dt
          endif
          write(obc_out_unit(m),*) 'Elevation phase speed, OBC ',Obc%bound(m)%name
          write(obc_out_unit(m),*) 'i = ',i,' j = ',j
          write(obc_out_unit(m),*) 'cgrid   = ',cdtdxr
          write(obc_out_unit(m),*) 'cmax    = ',min(Obc%bound(m)%ctrop_max * Obc%bound(m)%dsur(n), cdtdxr)
          write(obc_out_unit(m),*) 'cmin    = ',Obc%bound(m)%ctrop_min * Obc%bound(m)%dsur(n)
          write(obc_out_unit(m),*) 'ctrop   = ',ctrop(i,j)
       endif
    enddo

    ! Set the corner cells
    if ( south_west_corner ) then
       eta(i_sw,j_sw,taup1) = ( eta(i_sw+1,j_sw+1,taup1) &
                              + eta(i_sw+1,j_sw,  taup1) &
                              + eta(i_sw,  j_sw+1,taup1) ) /3.
    endif
    if ( south_east_corner ) then
       eta(i_se,j_se,taup1) = ( eta(i_se-1,j_se+1,taup1) &
                              + eta(i_se-1,j_se,  taup1) &
                              + eta(i_se,  j_se+1,taup1) ) /3.
    endif
    if ( north_west_corner ) then
       eta(i_nw,j_nw,taup1) = ( eta(i_nw+1,j_nw-1,taup1) &
                              + eta(i_nw+1,j_nw,  taup1) &
                              + eta(i_nw,  j_nw-1,taup1) ) /3.
    endif
    if ( north_east_corner ) then
       eta(i_ne,j_ne,taup1) = ( eta(i_ne-1,j_ne-1,taup1) &
                              + eta(i_ne-1,j_ne,  taup1) &
                              + eta(i_ne,  j_ne-1,taup1) ) /3.
    endif

    ! Update eta at global halo point to make the gradient accross boundary is 0
    call ocean_obc_update_boundary(eta(:,:,taup1), 'T')

    if(debug_this_module) then
       call write_chksum_2d('After ocean_obc_barotropic_dtbt, eta', eta(COMP,taup1))
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_barotropic
  !</SUBROUTINE> NAME="ocean_obc_barotropic"


  !######################################################################
  !--- update field on the halo points at the boundaries
  !<SUBROUTINE NAME="ocean_obc_update_boundary_2d" INTERFACE="ocean_obc_update_boundary">
  !   <INOUT NAME="field" TYPE="real, dimension(:,:)"></INOUT> 
  !   <IN NAME="grid_type" TYPE=" character(len=1)"></IN> 
  !   <IN NAME="update_type" TYPE=" character(len=1)"></IN> 
  !<DESCRIPTION>
  ! This subroutine sets values at the halo points beyond open boundaries
  ! grid_type can be:    'T' for variables at tracer point or 
  !                      'C' for variables at velocity points 
  !                      'Z' velocity points at zonal (north, south) boundaries
  !                      'M' velocity points at meridional (east, west) boundaries
  ! 
  ! update_type can be: 's' for setting the field beyond the 
  !                         boundary to the value at the boundary 
  !                         - for a velocity 'Z' at northern and southern 
  !                           boundary
  !                         - for a velocity 'M' at western and eastern 
  !                           boundary
  !                         - for gridtype 'T' and 'C' at all boundaries
  !                         This is equivalent to the NOGRAD OBC.
  !                     'x' for linear extrapolation of the field beyond the 
  !                         boundary from the value at the boundary and the 
  !                         first internal point
  !                         This is equivalent to the LINEAR OBC.
  !                     'i' for setting the field beyond the boundary to   
  !                         the value at the first internal point:
  !                         - for a velocity 'Z' at northern and southern 
  !                           boundary
  !                         - for a velocity 'M' at western and eastern 
  !                           boundary
  !                         - for gridtype 'T' and 'C' at all boundaries
  !                         This is equivalent to the INGRAD OBC.
  !                     'z' for setting the field beyond the boundary to   
  !                         zero:
  !                         - for a velocity 'Z' at northern and southern 
  !                           boundary
  !                         - for a velocity 'M' at western and eastern 
  !                           boundary
  !                         - for gridtype 'T' and 'C' at all boundaries
  !                         This is called CLAMPD OBC here.

  !                        
  !</DESCRIPTION>
  !<PUBLICROUTINE INTERFACE="ocean_obc_update_boundary" >
  subroutine ocean_obc_update_boundary_2d(field,grid_type,update_type)
  !</PUBLICROUTINE>

    real, dimension(isd:,jsd:), intent(inout)  :: field
    character(len=1),               intent(in) :: grid_type
    character(len=1), optional,     intent(in) :: update_type

    integer :: i, j, m, is, ie, js, je, i1, j1, n, ii, jj
    real    :: u1, u2, x1, x2
    integer, dimension(nobc) :: uptype     ! OBC type
    integer, dimension(:), pointer :: ilod, jlod, oi1, oj1
    integer, dimension(:), pointer :: istr, iend, jstr, jend

    ! Set the boundary condition
    uptype = NOGRAD
    do m = 1, nobc
      if(PRESENT(update_type)) then
        if (update_type == 's') then
          uptype(m) = NOGRAD
        else if (update_type == 'x') then
          uptype(m) = LINEAR
        else if (update_type == 'i') then
          uptype(m) = INGRAD
        else if (update_type == 'z') then
          uptype(m) = CLAMPD
        else if (update_type == 'n') then   ! Normal velocity components 
          uptype(m) = Obc%bound(m)%bcond_nor
        else if (update_type == 't') then   ! Tangential velocity components
          uptype(m) = Obc%bound(m)%bcond_tan
        else if (update_type == 'u') then   ! Nornal depth integrated velocity components
          uptype(m) = Obc%bound(m)%bcond_ud
        else if (update_type == 'c') then   ! Cell centered components
          uptype(m) = NOGRAD                ! Unimplemented : default NOGRAD 
        else
          call mpp_error(FATAL, &
         'ocean_obc_barotrop_mod(ocean_obc_update_boundary) : update_type= '// update_type// &
         ' should be either "n", "t", "s", "x", "z", "u",  or "i" ' )
        endif
      endif
    enddo

    if(grid_type .ne. 'T' .and. grid_type .ne. 'C'.and. grid_type .ne. 'Z'.and. grid_type .ne. 'M') call mpp_error(FATAL, &
     'ocean_obc_barotrop_mod(ocean_obc_update_boundary) : grid_type= '// grid_type//' should be either "T", "C", "Z" or "M" ')

    ! Loop through the boundaries and update
    do m = 1, nobc

      if(.not.Obc%bound(m)%on_bound) cycle
       if((Obc%bound(m)%direction == WEST .or. &
           Obc%bound(m)%direction == EAST) .and. grid_type == 'Z') cycle  
       if((Obc%bound(m)%direction == NORTH .or. &
           Obc%bound(m)%direction == SOUTH) .and. grid_type == 'M') cycle  

       ! Set pointers
       ilod => Obc%bound(m)%ilod
       jlod => Obc%bound(m)%jlod
       istr => Obc%bound(m)%istr
       jstr => Obc%bound(m)%jstr
       iend => Obc%bound(m)%iend
       jend => Obc%bound(m)%jend

       ! Implement the OBC
!       if (uptype(m) .eq. LINEAR) then  ! Linear extrapolation
       if (iand(uptype(m), LINEAR) == LINEAR) then   ! Linear extrapolation
          oi1 => Obc%bound(m)%di1
          oj1 => Obc%bound(m)%dj1
          if (Obc%bound(m)%direction == EAST .or. &
              Obc%bound(m)%direction == NORTH) then
             if(grid_type /= 'T') then
                ilod => Obc%bound(m)%di1
                jlod => Obc%bound(m)%dj1
                oi1 => Obc%bound(m)%di2
                oj1 => Obc%bound(m)%dj2
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          do n = 1, Obc%bound(m)%nlod
             i = ilod(n)    ! i boundary location
             j = jlod(n)    ! j boundary location
             i1 = oi1(n)    ! i interior cell location
             j1 = oj1(n)    ! j interior cell location
             is = istr(n)   ! Start of i halo loop
             ie = iend(n)   ! End of i halo loop
             js = jstr(n)   ! Start of j halo loop
             je = jend(n)   ! End of j halo loop
                u1 = field(i,j)
                u2 = field(i1,j1)
                x1 = Grd%xt(i,j)
                x2 = Grd%xt(i1,j1)
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = u1 + (u2-u1)*(Grd%xt(ii,jj) - x1)/(x2-x1)
                   enddo
                enddo
!!             endif
          enddo
!       else if (uptype(m) .eq. CLAMPD) then  ! Clamped to zero
       else if (iand(uptype(m), CLAMPD) == CLAMPD) then ! Clamped to zero
          if (Obc%bound(m)%direction == EAST .or. &
              Obc%bound(m)%direction == NORTH) then
             if(grid_type /= 'T') then
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          do n = 1, Obc%bound(m)%nlod
             is = istr(n)   ! Start of i halo loop
             ie = iend(n)   ! End of i halo loop
             js = jstr(n)   ! Start of j halo loop
             je = jend(n)   ! End of j halo loop
             do jj = js, je
                do ii = is, ie
                   field(ii,jj) = 0.0
                enddo
             enddo
          enddo
!       else if (uptype(m) .eq. NOGRAD) then   ! No-gradient
       else if (iand(uptype(m), NOGRAD) == NOGRAD) then ! No-gradient
          if (Obc%bound(m)%direction == EAST .or. &
              Obc%bound(m)%direction == NORTH) then
             if(grid_type /= 'T') then
                ilod => Obc%bound(m)%di1
                jlod => Obc%bound(m)%dj1
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          if(grid_type /= 'T') then
             do n = 1, Obc%bound(m)%nlod
                i = ilod(n)    ! Internal i boundary location
                j = jlod(n)    ! Internal j boundary location
                is = istr(n)   ! Start of i halo loop
                ie = iend(n)   ! End of i halo loop
                js = jstr(n)   ! Start of j halo loop
                je = jend(n)   ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)*u_mask(ii,jj)
                   enddo
                enddo
             enddo
          else
             do n = 1, Obc%bound(m)%nlod
                i = ilod(n)    ! Internal i boundary location
                j = jlod(n)    ! Internal j boundary location
                is = istr(n)   ! Start of i halo loop
                ie = iend(n)   ! End of i halo loop
                js = jstr(n)   ! Start of j halo loop
                je = jend(n)   ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)
                   enddo
                enddo
             enddo
          endif
!      else if (uptype(m) .eq. INGRAD) then  ! No-gradient to internal point
       else if (iand(uptype(m), INGRAD) == INGRAD) then ! No-gradient to internal point
          oi1 => Obc%bound(m)%di1
          oj1 => Obc%bound(m)%dj1
          if (Obc%bound(m)%direction == EAST .or. &
              Obc%bound(m)%direction == NORTH) then
             if(grid_type /= 'T') then
                oi1 => Obc%bound(m)%di2
                oj1 => Obc%bound(m)%dj2
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          if(grid_type /= 'T') then
             do n = 1, Obc%bound(m)%nlod
                i = oi1(n)     ! Internal i boundary location
                j = oj1(n)     ! Internal j boundary location
                is = istr(n)   ! Start of i halo loop
                ie = iend(n)   ! End of i halo loop
                js = jstr(n)   ! Start of j halo loop
                je = jend(n)   ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)*u_mask(ii,jj)
                   enddo
                enddo
             enddo
          else
             do n = 1, Obc%bound(m)%nlod
                i = oi1(n)     ! Internal i boundary location
                j = oj1(n)     ! Internal j boundary location
                is = istr(n)   ! Start of i halo loop
                ie = iend(n)   ! End of i halo loop
                js = jstr(n)   ! Start of j halo loop
                je = jend(n)   ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)
                   enddo
                enddo
             enddo
          endif
       else if (uptype(m) .eq. FILEIN) then   ! File input for velocity
          if (Obc%bound(m)%direction == EAST .or. &
              Obc%bound(m)%direction == NORTH) then
             if(grid_type /= 'T') then
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          if(grid_type == 'M' .or. grid_type == 'Z') then
             do n = 1, Obc%bound(m)%nlod
                i = Obc%bound(m)%ico_bt(n) ! i outside elevation face location
                j = Obc%bound(m)%jco_bt(n) ! j outside elevation face location
                is = istr(n)            ! Start of i halo loop
                ie = iend(n)            ! End of i halo loop
                js = jstr(n)            ! Start of j halo loop
                je = jend(n)            ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)*u_mask(ii,jj)
                   enddo
                enddo
             enddo
          endif
       endif

! There may be open points in the shadow of corners. These t-cells can be
! open or closed. It should be open, but who knows what happens ... 

      if ( south_west_corner ) then
         if(grid_type == 'T') then
            field(i_sw-1,j_sw-1) = field(i_sw,j_sw) * t_mask(i_sw,j_sw)
            field(max(isd,i_sw-2),j_sw-1) = field(i_sw,j_sw) * t_mask(i_sw,j_sw)
            field(i_sw-1,max(jsd,j_sw-2)) = field(i_sw,j_sw) * t_mask(i_sw,j_sw)
            field(max(isd,i_sw-2),max(jsd,j_sw-2)) = field(i_sw,j_sw) * t_mask(i_sw,j_sw)
         else
            field(i_sw-1,j_sw-1) = field(i_sw,j_sw) * u_mask(i_sw,j_sw)
            field(max(isd,i_sw-2),j_sw-1) = field(i_sw,j_sw) * u_mask(i_sw,j_sw)
            field(i_sw-1,max(jsd,j_sw-2)) = field(i_sw,j_sw) * u_mask(i_sw,j_sw)
            field(max(isd,i_sw-2),max(jsd,j_sw-2)) = field(i_sw,j_sw) * u_mask(i_sw,j_sw)
         endif
      endif
      if ( south_east_corner ) then
         if(grid_type == 'T') then
            field(i_se+1,j_se-1) = field(i_se,j_se) * t_mask(i_se,j_se)
            field(min(ied,i_se+2),j_se-1) = field(i_se,j_se) * t_mask(i_se,j_se)
            field(i_se+1,max(jsd,j_se-2)) = field(i_se,j_se) * t_mask(i_se,j_se)
            field(min(ied,i_se+2),max(jsd,j_se-2)) = field(i_se,j_se) * t_mask(i_se,j_se)
         else
            field(i_se+1,j_se-1) = field(i_se,j_se) * u_mask(i_se,j_se)
            field(min(ied,i_se+2),j_se-1) = field(i_se,j_se) * u_mask(i_se,j_se)
            field(i_se+1,max(jsd,j_se-2)) = field(i_se,j_se) * u_mask(i_se,j_se)
            field(min(ied,i_se+2),max(jsd,j_se-2)) = field(i_se,j_se) *u_mask(i_se,j_se)
         endif
      endif
      if ( north_west_corner ) then
         if(grid_type == 'T') then
            field(i_nw-1,j_nw+1) = field(i_nw,j_nw) * t_mask(i_nw,j_nw)
            field(i_nw-1,min(jed,j_nw+2)) = field(i_nw,j_nw) * t_mask(i_nw,j_nw)
            field(max(isd,i_nw-2),j_nw+1) = field(i_nw,j_nw) * t_mask(i_nw,j_nw)
            field(max(isd,i_nw-2),min(jed,j_nw+2)) = field(i_nw,j_nw) * t_mask(i_nw,j_nw)
         else
            field(i_nw-1,j_nw+1) = field(i_nw,j_nw) * u_mask(i_nw,j_nw)
            field(i_nw-1,min(jed,j_nw+2)) = field(i_nw,j_nw) * u_mask(i_nw,j_nw)
            field(max(isd,i_nw-2),j_nw+1) = field(i_nw,j_nw) * u_mask(i_nw,j_nw)
            field(max(isd,i_nw-2),min(jed,j_nw+2)) = field(i_nw,j_nw) * u_mask(i_nw,j_nw)
         endif
      endif
      if ( north_east_corner ) then
         if(grid_type == 'T') then
            field(i_ne+1,j_ne+1) = field(i_ne,j_ne) * t_mask(i_ne,j_ne)
            field(i_ne+1,min(jed,j_ne+2)) = field(i_ne,j_ne) * t_mask(i_ne,j_ne)
            field(min(ied,i_ne+2),j_ne+1) = field(i_ne,j_ne) * t_mask(i_ne,j_ne)
            field(min(ied,i_ne+2),min(jed,j_ne+2)) = field(i_ne,j_ne) * t_mask(i_ne,j_ne)
         else
            field(i_ne+1,j_ne+1) = field(i_ne,j_ne) * u_mask(i_ne,j_ne)
            field(i_ne+1,min(jed,j_ne+2)) = field(i_ne,j_ne) * u_mask(i_ne,j_ne)
            field(min(ied,i_ne+2),j_ne+1) = field(i_ne,j_ne) * u_mask(i_ne,j_ne)
            field(min(ied,i_ne+2),min(jed,j_ne+2)) = field(i_ne,j_ne) * u_mask(i_ne,j_ne)
         endif
      endif

   enddo

    return

  end subroutine ocean_obc_update_boundary_2d
  !</SUBROUTINE> NAME="ocean_obc_update_boundary_2d"


  !#####################################################################
  !--- update field on the halo points at the boundaries 
  !<SUBROUTINE NAME="ocean_obc_update_boundary_3d" INTERFACE="ocean_obc_update_boundary">
  !   <INOUT NAME="field" TYPE="real, dimension(:,:,:)"></INOUT> 
  subroutine ocean_obc_update_boundary_3d(field,grid_type,update_type)
    real, dimension(isd:,jsd:,:), intent(inout) :: field
    character(len=1),                intent(in) :: grid_type
    character(len=1), optional,      intent(in) :: update_type
    integer :: k, n3

    n3 = size(field,3)
    if(PRESENT(update_type)) then
      do k = 1, n3
        call ocean_obc_update_boundary_2d(field(:,:,k),grid_type,update_type)
      enddo
    else
      do k = 1, n3
        call ocean_obc_update_boundary_2d(field(:,:,k),grid_type)
      enddo
    endif
    return
  end subroutine ocean_obc_update_boundary_3d
  !</SUBROUTINE> NAME="ocean_obc_update_boundary_3d"


  !#####################################################################
  !--- set field on the boundaries to zero, do not change halos
  !<SUBROUTINE NAME="ocean_obc_zero_boundary_2d" INTERFACE="ocean_obc_zero_boundary">
  subroutine ocean_obc_zero_boundary_2d(field, grid_type)

    real, dimension(isd:,jsd:), intent(inout)  :: field
    character(len=1),           intent(in) :: grid_type

    integer :: i, j, m
    
    if(grid_type .ne. 'T' .and. grid_type .ne. 'C'.and. grid_type .ne. 'Z'.and. grid_type .ne. 'M') call mpp_error(FATAL, &
     'ocean_obc_barotrop_mod(ocean_obc_zero_boundary) : grid_type= '// grid_type//' should be either "T", "C", "Z" or "M" ')
    
    do m = 1, nobc
      !--- if on current pe there is no point on the bound, then just return
      if(.not.Obc%bound(m)%on_bound) cycle
            
      select case( Obc%bound(m)%direction )
      
      case (WEST)  
        if( grid_type == 'Z' ) cycle  

        i = Obc%bound(m)%is

        do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
          if(Obc%bound(m)%mask(j)) field(i,j) = 0.
        enddo

      case (EAST)
        if( grid_type == 'Z' ) cycle  
        
        if( grid_type == 'T' ) then
          i = Obc%bound(m)%is 
        else
          i = Obc%bound(m)%is - 1
        endif  
        do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
          if(Obc%bound(m)%mask(j)) field(i,j) = 0.
        enddo

      case (SOUTH)
        if( grid_type == 'M' ) cycle  

        j = Obc%bound(m)%js 
        do i = Obc%bound(m)%isd, Obc%bound(m)%ied
          if(Obc%bound(m)%mask(i)) field(i,j) = 0.
        enddo
        
      case (NORTH)
        if( grid_type == 'M' ) cycle  

        if( grid_type == 'T' ) then
          j = Obc%bound(m)%js 
        else
          j = Obc%bound(m)%js - 1
        endif  
        do i = Obc%bound(m)%isd, Obc%bound(m)%ied
          if(Obc%bound(m)%mask(i)) field(i,j) = 0.
        enddo
        
      end select

    enddo

    return

  end subroutine ocean_obc_zero_boundary_2d
  !</SUBROUTINE> NAME="ocean_obc_zero_boundary_2d"


  !#####################################################################
! <SUBROUTINE NAME="ocean_obc_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_obc_restart(time_stamp, ctrop_chksum)
  character(len=*), intent(in), optional    :: time_stamp
  integer(LONG_KIND), intent(out), optional :: ctrop_chksum

   if(present(ctrop_chksum)) ctrop_chksum = mpp_chksum(ctrop(COMP))
   call reset_field_pointer(Obc_restart, id_restart(1), ctrop(:,:) )

   call save_restart(Obc_restart, time_stamp)

end subroutine ocean_obc_restart
! </SUBROUTINE> NAME="ocean_obc_restart"


  !#####################################################################
  !--- release memory --------------------------------------------------
  !<SUBROUTINE NAME="ocean_obc_end">
  ! <DESCRIPTION>
  !    Destructor routine. Release memory.
  !   </DESCRIPTION>
  !   <OUT NAME="have_obc" TYPE="logical">
  !      Contains open boundary information
  !   </OUT>

  subroutine ocean_obc_end(Time)

    type(ocean_time_type), intent(in) :: Time
    integer(LONG_KIND)                :: ctrop_chksum
    integer :: m

  integer :: stdoutunit 
  stdoutunit=stdout() 

    do m = 1, nobc
      if(obc_out_unit(m) .NE. stdoutunit ) then
         call mpp_close(obc_out_unit(m))
      endif
    enddo   
    call ocean_obc_restart(ctrop_chksum=ctrop_chksum)
    
    write(stdoutunit,*) ' '
    write(stdoutunit,*) &
        'From ocean_obc_barotrop_mod: ending obc chksums'
    call write_timestamp(Time%model_time)

    write(stdoutunit,*) &
        'chksum for phase speed            = ', ctrop_chksum

    call write_chksum_2d('relaxation coefficient', rel_coef(COMP))
    call write_chksum_2d('eta_tm1', eta_tm1(COMP))
    call write_chksum_2d('eta_t_tm1', eta_t_tm1(COMP))
    
    deallocate( Obc%bound )

    module_is_initialized = .FALSE.
    
    nullify(Dom)
    
    return

  end subroutine ocean_obc_end
  !</SUBROUTINE> NAME="ocean_obc_end"


  !#######################################################################
  !<SUBROUTINE NAME="phase_speed_IOW">
  subroutine phase_speed_IOW(Bound, Data, eta, taum1, taup1, tstep, init)

    type(obc_bound_type),         intent(inout) :: Bound
    type(obc_data_type_2d) ,      intent(in)    :: Data
    real, dimension(isd:,jsd:,:), intent(in)    :: eta
    integer,                      intent(in)    :: taum1, taup1
    real, intent(in)                            :: tstep
    logical, intent(in), optional               :: init

    real    :: cgrid, detai, cmax, cmin, cinc, c1tmp, rfak
    integer :: i, j
    integer :: n, i1, i2, j1, j2
    real, dimension(:,:), pointer :: dg
    real :: dt
    real :: scale, f1, f2
    logical :: init_f

    ! Set variables for calculating the barotropic phase speed in the 
    !first call of the barotropic OBC. Usually, eta is Ext_mode%eta_t 
    ! or Ext_mode%eta_t_bar with splitting on.
    init_f = .false.
    if (PRESENT(init)) init_f = init
    
    if (init_f) then
       dt = dtbt
       scale = dtbt / tstep
       f1 = 0.0
       f2 = 1.0
    else
       dt = tstep
       scale = 1.0
       f1 = Bound%ctrop_smooth
       f2 = 1. - Bound%ctrop_smooth
    endif

    ! If on current pe there is no point on the bound, then just return
    if(.not. Bound%on_bound) return
    
    ! Set pointers
    if(Bound%direction == WEST .or. Bound%direction == EAST) then
       dg => dxu_bt
    else
       dg => dyu_bt
    endif
    
    do n = 1, Bound%nlod

       i = Bound%ilod(n)
       j = Bound%jlod(n)
       i1 = i  + Bound%imap
       j1 = j  + Bound%jmap
       i2 = i1 + Bound%imap
       j2 = j1 + Bound%jmap
 
       cgrid = dg(i1,j1)/dt
       detai = (eta(i1,j1,taup1) - eta(i2,j2,taup1) )
       ! clip with shallow water phase speed and CFL-phase speed
       cmax  = min(Bound%ctrop_max * Bound%dsur(n), cgrid)
       cmin  = Bound%ctrop_min * Bound%dsur(n)
       cinc  = Bound%ctrop_inc * cgrid
       if (ABS(detai).lt. small) then
          if (init_f) then
             ! badly defined phases -> move this feature out with 
             ! barotropic phase speed
             c1tmp = cmax
          else
             ! badly defined phases -> do (about) nothing
             c1tmp = 0.99*ctrop(i,j)
          endif
       else
          ! implicit time scheme for internal points
          c1tmp = cgrid*(eta(i1,j1,taum1)-eta(i1,j1,taup1))/detai * scale
       endif
       ! If incoming waves are detected:
       ! -> replace undesired incoming waves by outgoing waves
        if (c1tmp .le. 0.0) then
          ctrop(i,j) = cinc
       else
          c1tmp = min(cmax,c1tmp)  ! cmax is the phase speed limit
          c1tmp = max(cmin,c1tmp)  ! ensure minimum phase speed
          ! mix with previous phase values
          ctrop(i,j) = f1*ctrop(i,j) + f2*c1tmp
       endif
       ! Find the relaxation coefficient
       ! Relax with rel_coef_in for incoming waves and rel_coef_out 
       ! otherwise a smooth transition is needed            
        rfak = ctrop(i,j)/cmax
       rel_coef(i,j) = rfak * Data%rel_coef_out + (1.0-rfak) * Data%rel_coef_in
    enddo

!kk  Changes for Halos > 1 have to be done
   ! Corner cells
    if ( south_west_corner ) then
       rel_coef(i_sw,j_sw) =  (rel_coef(i_sw+1,j_sw)+rel_coef(i_sw,j_sw+1))*0.5
    endif
    if ( south_east_corner ) then
       rel_coef(i_se,j_se) =  (rel_coef(i_se-1,j_se)+rel_coef(i_se,j_se+1))*0.5
    endif
    if ( north_east_corner ) then
       rel_coef(i_ne,j_ne) =  (rel_coef(i_ne-1,j_ne)+rel_coef(i_ne,j_ne-1))*0.5
    endif
    if ( north_west_corner ) then
       rel_coef(i_nw,j_nw) =  (rel_coef(i_nw+1,j_nw)+rel_coef(i_nw,j_nw-1))*0.5
    endif

    return
  end subroutine phase_speed_IOW
  !</SUBROUTINE> NAME="phase_speed_IOW"


  !#######################################################################
  !<SUBROUTINE NAME="phase_speed_ORLANS">
  subroutine phase_speed_ORLANS(Bound, Data, eta, tau, taup1, tstep, init)

    type(obc_bound_type),         intent(inout) :: Bound
    type(obc_data_type_2d) ,      intent(in)    :: Data
    real, dimension(isd:,jsd:,:), intent(in)    :: eta
    integer,                      intent(in)    :: tau, taup1
    real, intent(in)                            :: tstep
    logical, intent(in), optional               :: init

    real    :: cgrid, detai, cmax, cmin, cinc, c1tmp, rfak
    integer :: i, j
    integer :: n, i1, i2, j1, j2
    real, dimension(:,:), pointer :: dg, eta_taum1
    real :: dt
    real :: scale, f1, f2
    logical :: init_f

    ! Set variables for calculating the barotropic phase speed in the 
    !first call of the barotropic OBC. Usually, eta is Ext_mode%eta_t 
    ! or Ext_mode%eta_t_bar with splitting on.
    init_f = .false.
    if (PRESENT(init)) init_f = init
    
    if (init_f) then
       dt = dtbt
       scale = dtbt / tstep
       f1 = 0.0
       f2 = 1.0
       eta_taum1 => eta_t_tm1
    else
       dt = tstep
       scale = 1.0
       f1 = Bound%ctrop_smooth
       f2 = 1. - Bound%ctrop_smooth
       eta_taum1 => eta_tm1
    endif

    ! If on current pe there is no point on the bound, then just return
    if(.not. Bound%on_bound) return

    ! Set pointers
    if(Bound%direction == WEST .or. Bound%direction == EAST) then
       dg => dxu_bt
    else
       dg => dyu_bt
    endif

    do n = 1, Bound%nlod

       i = Bound%ilod(n)
       j = Bound%jlod(n)
       i1 = i + Bound%imap
       j1 = j + Bound%jmap
       i2 = i1 + Bound%imap
       j2 = j1 + Bound%jmap

       cgrid = dg(i1,j1)/dt
       detai = eta(i1,j1,taup1) + eta_taum1(i1,j1) - 2.0 * eta(i2,j2,tau)
       ! clip with shallow water phase speed and CFL-phase speed
       cmax  = min(Bound%ctrop_max * Bound%dsur(n), cgrid)
       cmin  = Bound%ctrop_min * Bound%dsur(n)
       cinc  = Bound%ctrop_inc * cgrid
       if (ABS(detai).lt. small) then
          if (init_f) then
             ! badly defined phases -> move this feature out with 
             ! barotropic phase speed
             c1tmp = cmax
          else
             ! badly defined phases -> do (about) nothing
             c1tmp = ctrop(i,j)
          endif
       else
          ! implicit time scheme for internal points
          c1tmp = cgrid*(eta_taum1(i1,j1)-eta(i1,j1,taup1))/detai * scale
       endif
       ! If incoming waves are detected:
       ! -> replace undesired incoming waves by outgoing waves
       if (c1tmp .le. 0.0) then
          ctrop(i,j) = cinc
       else
          c1tmp = min(cmax,c1tmp)  ! cmax is the phase speed limit
          c1tmp = max(cmin,c1tmp)  ! ensure minimum phase speed
          ! mix with previous phase values
          ctrop(i,j) = f1*ctrop(i,j) + f2*c1tmp
       endif
       ! Find the relaxation coefficient
       ! Relax with rel_coef_in for incoming waves and rel_coef_out 
       ! otherwise a smooth transition is needed            
       rfak = ctrop(i,j)/cmax
       rel_coef(i,j) = rfak * Data%rel_coef_out + (1.0-rfak) * Data%rel_coef_in
    enddo

    ! Corner cells
    if ( south_west_corner ) then
       rel_coef(i_sw,j_sw) =  (rel_coef(i_sw+1,j_sw)+rel_coef(i_sw,j_sw+1))*0.5
    endif
    if ( south_east_corner ) then
       rel_coef(i_se,j_se) =  (rel_coef(i_se-1,j_se)+rel_coef(i_se,j_se+1))*0.5
    endif
    if ( north_east_corner ) then
       rel_coef(i_ne,j_ne) =  (rel_coef(i_ne-1,j_ne)+rel_coef(i_ne,j_ne-1))*0.5
    endif
    if ( north_west_corner ) then
       rel_coef(i_nw,j_nw) =  (rel_coef(i_nw+1,j_nw)+rel_coef(i_nw,j_nw-1))*0.5
    endif

    return
  end subroutine phase_speed_ORLANS
  !</SUBROUTINE> NAME="phase_speed_ORLANS"


  !#######################################################################
  !<SUBROUTINE NAME="phase_speed_GRAVTY">
  subroutine phase_speed_GRAVTY(Bound, Data)

    type(obc_bound_type),         intent(inout) :: Bound
    type(obc_data_type_2d),       intent(in)    :: Data
    
    integer :: i, j, n

    ! If on current pe there is no point on the bound, then just return
    if(.not. Bound%on_bound) return

    do n = 1, Bound%nlod
       i  = Bound%ilod(n)
       j  = Bound%jlod(n)
       ctrop(i,j) = Bound%dsur(n)
       ! Find the relaxation coefficient
       ! Relax with rel_coef_out 
       rel_coef(i,j) = Data%rel_coef_out
    enddo

    ! Corner cells
    if ( south_west_corner ) then
       rel_coef(i_sw,j_sw) =  (rel_coef(i_sw+1,j_sw)+rel_coef(i_sw,j_sw+1))*0.5
    endif
    if ( south_east_corner ) then
       rel_coef(i_se,j_se) =  (rel_coef(i_se-1,j_se)+rel_coef(i_se,j_se+1))*0.5
    endif
    if ( north_east_corner ) then
       rel_coef(i_ne,j_ne) =  (rel_coef(i_ne-1,j_ne)+rel_coef(i_ne,j_ne-1))*0.5
    endif
    if ( north_west_corner ) then
       rel_coef(i_nw,j_nw) =  (rel_coef(i_nw+1,j_nw)+rel_coef(i_nw,j_nw-1))*0.5
    endif

    return
  end subroutine phase_speed_GRAVTY
  !</SUBROUTINE> NAME="phase_speed_GRAVTY"


  !#######################################################################
  !<SUBROUTINE NAME="phase_speed_MILLER">
  subroutine phase_speed_MILLER(Bound, Data, eta, tau, taup1, tstep, init)

    type(obc_bound_type),         intent(inout) :: Bound
    type(obc_data_type_2d) ,      intent(in)    :: Data
    real, dimension(isd:,jsd:,:), intent(in)    :: eta
    integer,                      intent(in)    :: tau, taup1
    real, intent(in)                            :: tstep
    logical, intent(in), optional               :: init

    real    :: cgrid, cmax, cmin, cinc, c1tmp, rfak
    integer :: i, j
    integer :: n, i1, i2, j1, j2
    real, dimension(:,:), pointer :: dg, eta_taum1
    real :: dt
    real :: scale, f1, f2, r1, r2, r3
    logical :: init_f

    ! Set variables for calculating the barotropic phase speed in the 
    !first call of the barotropic OBC. Usually, eta is Ext_mode%eta_t 
    ! or Ext_mode%eta_t_bar with splitting on.
    init_f = .false.
    if (PRESENT(init)) init_f = init
    
    if (init_f) then
       dt = dtbt
       scale = dtbt / tstep
       f1 = 0.0
       f2 = 1.0
       eta_taum1 => eta_t_tm1
    else
       dt = tstep
       scale = 1.0
       f1 = Bound%ctrop_smooth
       f2 = 1. - Bound%ctrop_smooth
       eta_taum1 => eta_tm1
    endif

    ! If on current pe there is no point on the bound, then just return
    if(.not. Bound%on_bound) return

    ! Set pointers
    if(Bound%direction == WEST .or. Bound%direction == EAST) then
       dg => dxu_bt
    else
       dg => dyu_bt
    endif

    do n = 1, Bound%nlod

       i = Bound%ilod(n)
       j = Bound%jlod(n)
       i1 = i + Bound%imap
       j1 = j + Bound%jmap
       i2 = i1 + Bound%imap
       j2 = j1 + Bound%jmap

       cgrid = dg(i1,j1)/dt

       r1 = eta(i2,j2,tau) - eta(i1,j1,tau)
       r2 = eta_taum1(i1,j1) - eta_taum1(i,j)
       r3 = eta_taum1(i2,j2) - eta_taum1(i1,j1)
       ! clip with shallow water phase speed and CFL-phase speed
       cmax  = min(Bound%ctrop_max * Bound%dsur(n), cgrid)
       cmin  = Bound%ctrop_min * Bound%dsur(n)
       cinc  = Bound%ctrop_inc * cgrid
       if (ABS(r1).lt. small) then
          r1 = cmax
       else
          r1 = (eta(i1,j1,taup1) - eta(i1,j1,tau)) / r1
       endif
       if (ABS(r2).lt. small) then
          r2 = cmax
       else
          r2 = (eta(i,j,tau) - eta_taum1(i,j)) / r2
       endif
       if (ABS(r3).lt. small) then
          r3 = cmax
       else
          r3 = (eta(i1,j1,tau) - eta_taum1(i1,j1)) / r3
       endif
       c1tmp = cgrid * (r1 + r2 - r3) * scale;
       
       ! If incoming waves are detected:
       ! -> replace undesired incoming waves by outgoing waves
       if (c1tmp .le. 0.0) then
          ctrop(i,j) = cinc
       else
          c1tmp = min(cmax,c1tmp)  ! cmax is the phase speed limit
          c1tmp = max(cmin,c1tmp)  ! ensure minimum phase speed
          ! mix with previous phase values
          ctrop(i,j) = f1*ctrop(i,j) + f2*c1tmp
       endif
       ! Find the relaxation coefficient
       ! Relax with rel_coef_in for incoming waves and rel_coef_out 
       ! otherwise a smooth transition is needed            
       rfak = ctrop(i,j)/cmax
       rel_coef(i,j) = rfak * Data%rel_coef_out + (1.0-rfak) * data%rel_coef_in
    enddo

    ! Corner cells
    if ( south_west_corner ) then
       rel_coef(i_sw,j_sw) =  (rel_coef(i_sw+1,j_sw)+rel_coef(i_sw,j_sw+1))*0.5
    endif
    if ( south_east_corner ) then
       rel_coef(i_se,j_se) =  (rel_coef(i_se-1,j_se)+rel_coef(i_se,j_se+1))*0.5
    endif
    if ( north_east_corner ) then
       rel_coef(i_ne,j_ne) =  (rel_coef(i_ne-1,j_ne)+rel_coef(i_ne,j_ne-1))*0.5
    endif
    if ( north_west_corner ) then
       rel_coef(i_nw,j_nw) =  (rel_coef(i_nw+1,j_nw)+rel_coef(i_nw,j_nw-1))*0.5
    endif

    return
  end subroutine phase_speed_MILLER
  !</SUBROUTINE> NAME="phase_speed_MILLER"


  !#######################################################################
  !<SUBROUTINE NAME="phase_speed_RAYMND">
  subroutine phase_speed_RAYMND(Bound, Data, eta, tau, taup1, tstep, init)

    type(obc_bound_type),         intent(inout) :: Bound
    type(obc_data_type_2d) ,      intent(in)    :: Data
    real, dimension(isd:,jsd:,:), intent(in)    :: eta
    integer,                      intent(in)    :: tau, taup1
    real, intent(in)                            :: tstep
    logical, intent(in), optional               :: init


    real    :: cgrid, detat, cmax, cmin, cinc, c1tmp, rfak
    real    :: rx, ry, dx, dy, d, dmp, dmm, dm
    integer :: i, j
    integer :: n, i1, i2, j1, j2
    integer :: ip, jp, im, jm, i1p, j1p, i1m, j1m
    real, dimension(:,:), pointer :: dg, eta_taum1
    real :: dt
    real :: scale, f1, f2
    logical :: init_f

    ! Set variables for calculating the barotropic phase speed in the 
    !first call of the barotropic OBC. Usually, eta is Ext_mode%eta_t 
    ! or Ext_mode%eta_t_bar with splitting on.
    init_f = .false.
    if (PRESENT(init)) init_f = init
    
    if (init_f) then
       dt = dtbt
       scale = dtbt / tstep
       f1 = 0.0
       f2 = 1.0
       eta_taum1 => eta_t_tm1
    else
       dt = tstep
       scale = 1.0
       f1 = Bound%ctrop_smooth
       f2 = 1. - Bound%ctrop_smooth
       eta_taum1 => eta_tm1
    endif

    ! If on current pe there is no point on the bound, then just return
    if(.not. Bound%on_bound) return

    if(Bound%direction == WEST .or. Bound%direction == EAST) then
       dg => dxu_bt
    else
       dg => dyu_bt
    endif

    do n = 1, Bound%nloc

       i = Bound%iloc(n)
       j = Bound%jloc(n)
       i1 = i + Bound%imap
       j1 = j + Bound%jmap
       i2 = i1 + Bound%imap
       j2 = j1 + Bound%jmap
       ip = i + Bound%imapt
       jp = j + Bound%jmapt
       im = i - Bound%imapt
       jm = j - Bound%jmapt
       i1p = i1 + Bound%imapt
       j1p = j1 + Bound%jmapt
       i1m = i1 - Bound%imapt
       j1m = j1 - Bound%jmapt
       dmp = t_mask(i1p,j1p)
       dmm = t_mask(i1m,j1m)
       dm  = t_mask(i1,j1)

       cgrid = dg(i1,j1)/dt !! scaling factor for phase speeds
                            !! The maximum non-dimensional CFL  phase speed, c*
                            !! is unity where c* = c[m s-1]/cgrid[m s-1]        
                              

! clip with shallow water phase speed and CFL-phase speed
! csur = sqrt(gH)
       
       cmax  = min(Bound%ctrop_max * Bound%dsur(n), cgrid)
       cmin  = Bound%ctrop_min * Bound%dsur(n)
       cinc  = Bound%ctrop_inc * cgrid

       ! Difference in time
       detat = eta(i1,j1,taup1) - eta_taum1(i1,j1)

       ! Difference in space
       dx = eta(i1,j1,taup1) - eta(i2,j2,taup1)
       dy = detat * (eta(i1p,j1p,tau) - eta(i1m,j1m,tau))
       if (dy > 0.0) then
          dy = eta(i1,j1,tau) - eta(i1m,j1m,tau) ! [m]
       else
          dy = eta(i1p,j1p,tau) - eta(i1,j1,tau) 
       endif
! set dy to zero, if one of the tangent points is land. Otherwise either
! the first check is wrong or dy is not well defined.
! This does not fit well in the B-grid scheme ... 
       dy = dy * dm * dmp * dmm 
       d = dx * dx + dy * dy 

       ! Normal phase speed
       if (ABS(d).lt. small) then
          rx = cmax
       else
          rx = max(0.0, -detat * dx / d)          ! [nd]
!!!!          rx = min(max(0.0, -detat * dx / d), cmax)
       endif

      ! [m s-1]
       c1tmp = cgrid * rx * scale

       ! Tangential phase speed
       ry = 0.0
       if (rx /= 0 .and. d /= 0) ry = -detat * dy / d
       if (ry > 0.0) then
          ry = min(ry,Bound%ctrop_max)
          d = eta(i,j,tau) - eta(im,jm,tau)
       else
          ry = max(ry,-1*Bound%ctrop_max)          
          d = eta(ip,jp,tau) - eta(i,j,tau)
       endif
! The tangent component must be zero, if one of the neigbours is land. 
       Bound%dum(n) = ry *d * dm * dmp * dmm
       ! Note : eta(i,j,taup1) = (eta(i,j,tau) +  &
       !                          rx * eta(i1,j1,taup1) - ry * d) / &
       !                         (1 + rx)

       ! If incoming waves are detected:
       ! -> replace undesired incoming waves by outgoing waves
       if (c1tmp .le. 0.0) then
          ctrop(i,j) = cinc
       else
          c1tmp = min(cmax,c1tmp)  ! cmax is the phase speed limit
          c1tmp = max(cmin,c1tmp)  ! ensure minimum phase speed
          ! mix with previous phase values
          ctrop(i,j) = f1*ctrop(i,j) + f2*c1tmp
       endif
       ! Find the relaxation coefficient
       ! Relax with rel_coef_in for incoming waves and rel_coef_out 
       ! otherwise a smooth transition is needed            
       rfak = ctrop(i,j)/cmax
       rel_coef(i,j) = rfak * Data%rel_coef_out + (1.0-rfak) * data%rel_coef_in

    enddo

    ! Corner cells
    if ( south_west_corner ) then
       rel_coef(i_sw,j_sw) =  (rel_coef(i_sw+1,j_sw)+rel_coef(i_sw,j_sw+1))*0.5
    endif
    if ( south_east_corner ) then
       rel_coef(i_se,j_se) =  (rel_coef(i_se-1,j_se)+rel_coef(i_se,j_se+1))*0.5
    endif
    if ( north_east_corner ) then
       rel_coef(i_ne,j_ne) =  (rel_coef(i_ne-1,j_ne)+rel_coef(i_ne,j_ne-1))*0.5
    endif
    if ( north_west_corner ) then
       rel_coef(i_nw,j_nw) =  (rel_coef(i_nw+1,j_nw)+rel_coef(i_nw,j_nw-1))*0.5
    endif

    return
  end subroutine phase_speed_RAYMND
  !</SUBROUTINE> NAME="phase_speed_RAYMND"


  !#######################################################################
  !<FUNCTION NAME="boundary_average">
  function boundary_average(Bound, data)
     real,    dimension(:),    intent(in) :: data
     type(obc_bound_type),  intent(inout) :: Bound             ! lateral boundary data
     real                                 :: boundary_average     
     integer                              :: nlist, n, sendsize, index

     if(Bound%on_bound) then
        nlist = size(Bound%pelist(:))
        index = Bound%index
        sendsize = Bound%end(Bound%index)-Bound%start(index) + 1
        if(size(data) .NE. sendsize ) call mpp_error(FATAL, &
             "ocean_obc_barotrop_mod: size mismatch between data and boundary index")
        Bound%buffer(Bound%start(index):Bound%end(index)) = data(:)
        do n = 1, nlist
           if(mpp_pe() .NE. Bound%pelist(n) ) then
              call mpp_send(Bound%buffer(Bound%start(index)), plen=sendsize, to_pe=Bound%pelist(n))
           end if
        end do

        !--- receive the data
        do n = 1, nlist
           if(mpp_pe() .NE. Bound%pelist(n) ) then
              call mpp_recv(Bound%buffer(Bound%start(n)), glen=Bound%end(n)-Bound%start(n)+1, from_pe=Bound%pelist(n))
           end if
        end do
     end if

     call mpp_sync_self(Bound%pelist)
     if(Bound%on_bound) then
        boundary_average = sum(Bound%buffer)/Bound%np
     else
        boundary_average = 0
     end if

  end function boundary_average
  !</FUNCTION> NAME="boundary_average"


  !#######################################################################
  !<SUBROUTINE NAME="check_eta_OBC">
  subroutine check_eta_OBC(obc, name)
     integer, intent(in)                :: obc
     character(len=*), intent(in)       :: name
     integer                            :: n, m, nl, code;

     if(obc == 0) &
          call mpp_error(FATAL,' Invalid eta OBC: '//name)
     nl = 1
     n = 0
     code = 0
     if (index(trim(name),'|') /= 0) then
        nl = 2
        n = -1
     endif
     do m = 1, nl
        if (iand(code, NOTHIN) /= NOTHIN .and. index(trim(name),'NOTHIN') /= 0) then
           n = n + 1
           code = code + NOTHIN
        endif
        if (iand(code, FILEIN) /= FILEIN .and. index(trim(name),'FILEIN') /= 0) then
           n = n + 1
           code = code + FILEIN
        endif
        if (iand(code, MEANIN) /= MEANIN .and. index(trim(name),'MEANIN') /= 0) then
           n = n + 1
           code = code + MEANIN
        endif
        if (iand(code, FLATHR) /= FLATHR .and. index(trim(name),'FLATHR') /= 0) then
           n = n + 1
           code = code + FLATHR
        endif
        if (iand(code, NOGRAD) /= NOGRAD .and. index(trim(name),'NOGRAD') /= 0) then
           n = n + 1
           code = code + NOGRAD
        endif
        if (iand(code, ORLANS) /= ORLANS .and. index(trim(name),'ORLANS') /= 0) then
           n = n + 1
           code = code + ORLANS
        endif
        if (iand(code, CAMOBR) /= CAMOBR .and. index(trim(name),'CAMOBR') /= 0) then
           n = n + 1
           code = code + CAMOBR
        endif
        if (iand(code, GRAVTY) /= GRAVTY .and. index(trim(name),'GRAVTY') /= 0) then
           n = n + 1
           code = code + GRAVTY
        endif
        if (iand(code, MILLER) /= MILLER .and. index(trim(name),'MILLER') /= 0) then
           n = n + 1
           code = code + MILLER
        endif
        if (iand(code, RAYMND) /= RAYMND .and. index(trim(name),'RAYMND') /= 0) then
           n = n + 1
           code = code + RAYMND
        endif
        if (iand(code, IOW) /= IOW .and. index(trim(name),'IOW') /= 0) then
           n = n + 1
           code = code + IOW
        endif
     enddo
     if (n /= 1) &
          call mpp_error(FATAL,' Invalid eta OBC: '//name)
     
  end subroutine check_eta_OBC
  !</SUBROUTINE> NAME="check_eta_OBC"


  !#######################################################################
  !<SUBROUTINE NAME="mpp_update_domains_obc">
  subroutine mpp_update_domains_obc(Domain)
    type (ocean_domain_type),intent(inout)           :: Domain

! make an incomlete communication, it will be completed in the calling subroutine
      call mpp_update_domains(ctrop   , Domain%domain2d, complete=.false.)
      if (update_eta_tm1)   &
      call mpp_update_domains(eta_tm1 , Domain%domain2d, complete=.false.)
    return
  end subroutine mpp_update_domains_obc

end module ocean_obc_barotrop_mod
 !</SUBROUTINE> NAME="mpp_update_domains_obc"
