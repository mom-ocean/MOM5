module ocean_obc_mod
#define COMP isc:iec,jsc:jec
  !
  ! <CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt </CONTACT>
  ! <CONTACT EMAIL="Mike.Herzfeld@csiro.au"> Mike Herzfeld </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matthew Harrison </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies </REVIEWER>

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
  use mpp_domains_mod,          only: BGRID_NE, mpp_define_domains

  use mpp_mod,                  only: input_nml_file, mpp_error, mpp_pe
  use mpp_mod,                  only: CLOCK_MODULE
  use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use time_interp_external_mod, only: time_interp_external, init_external_field, get_external_field_size
  use time_manager_mod,         only: time_type
  use tracer_manager_mod,       only: get_tracer_names, get_tracer_indices, get_number_tracers
  use ocean_util_mod,           only: write_timestamp, write_chksum_2d, write_chksum_3d

  use ocean_domains_mod,        only: get_local_indices, get_domain_offsets
  use ocean_parameters_mod,     only: missing_value, rho0, GEOPOTENTIAL
  use ocean_types_mod,          only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
  use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_adv_vel_type, ocean_options_type
  use ocean_types_mod,          only: ocean_external_mode_type, obc_flux
  use ocean_types_mod,          only: ocean_time_type, ocean_time_steps_type
  use ocean_obc_barotrop_mod,   only: ocean_obc_prepare_baro => ocean_obc_prepare
  use ocean_obc_barotrop_mod,   only: ocean_obc_end_baro => ocean_obc_end

implicit none
  private

interface
   integer function tm_scale_to_secs(buf, sec)
     character(len=16) :: buf
     real              :: sec
   end function tm_scale_to_secs
end interface


  public :: ocean_obc_init, ocean_obc_end, ocean_obc_tracer
  public :: ocean_obc_mixing, ocean_obc_surface_height
  public :: ocean_obc_update_boundary, ocean_obc_prepare
  public :: ocean_obc_zero_boundary
  public :: ocean_obc_tracer_flux, ocean_obc_mass_flux
  public :: store_ocean_obc_tracer_flux
  public :: ocean_obc_tracer_init, ocean_obc_adjust_advel, ocean_obc_adjust_divud
  public :: store_ocean_obc_pressure_grad, ocean_obc_adjust_forcing_bt
  public :: ocean_obc_enhance_visc_back
  public :: ocean_obc_enhance_diff_back
  public :: ocean_obc_restart

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
  integer            :: vert_coordinate         ! used to limit vertical loops in tracer for GEOPOTENTIAL

  ! for restart

  type(ocean_domain_type), pointer :: Dom => NULL()  ! ocean domain
  type(ocean_grid_type), pointer   :: Grd => NULL()  ! ocean grid
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
     ! For western and eastern boundaries, d1i = 1, d1j => jloc
     ! For southern and northern boundaries, d1i => iloc, d1j = 1
     integer, pointer  :: d1i(:)     => NULL()   ! 1D i index pointer
     integer, pointer  :: d1j(:)     => NULL()   ! 1D j index pointer

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
     ! For western boundaries ico => iloc-imap, jcf => jloc
     ! For eastern boundaries ico => iloc, jcf => jloc
     ! For southern boundaries ico => iloc, jcf => jloc-jmap
     ! For northern boundaries ico => iloc, jcf => jloc
     integer, pointer  :: ico(:)     => NULL()   ! Outside face i location
     integer, pointer  :: jco(:)     => NULL()   ! Outside face j location

     ! Boundary vectors containing the normal (nic) and tangential (tic)
     ! indicies to the boundary.
     ! For west and east boundaries nic => iloc, tic => jloc
     ! For south and west boundaries nic => jloc, tic => iloc
     integer, pointer  :: nic(:)     => NULL()   ! Pointer to normal index
     integer, pointer  :: tic(:)     => NULL()   ! Pointer to tangential index

     integer, pointer  :: ones(:)    => NULL()   ! Array of ones

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
     real,    pointer  :: dumv(:,:)  => NULL()   ! 2D dummy array
     real,    pointer  :: ctrop(:)   => NULL()   ! barotropic phase speed
     real,    pointer  :: pgrad(:)   => NULL()   ! barotropic phase speed
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
     real, pointer     :: csur(:) => NULL()       ! Gravity wave speed
     integer           :: id_xt, id_yt, id_xu, id_yu
     integer           :: id_ctrop, id_eta_data, id_transport
     integer, allocatable :: id_tracer_data(:), id_tracer_flux(:), id_cclin(:)
     integer, allocatable :: id_tracer_flux_dif(:), id_tracer_flux_adv(:)
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
     module procedure ocean_obc_update_boundary_4d
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
     module procedure ocean_obc_zero_boundary_3d
  end interface
  !</INTERFACE>

  !<INTERFACE NAME="ocean_obc_enhance_visc_back">
  !   <DESCRIPTION>
  !   enhance viscosity near open boundaries
  !   </DESCRIPTION>
  interface ocean_obc_enhance_visc_back
     module procedure ocean_obc_enhance_visc_back_3d
     module procedure ocean_obc_enhance_visc_back_2d
  end interface
  !</INTERFACE>

  !<INTERFACE NAME="ocean_obc_enhance_diff_back">
  !   <DESCRIPTION>
  !   enhance diffusion near open boundaries
  !   </DESCRIPTION>
  interface ocean_obc_enhance_diff_back
     module procedure ocean_obc_enhance_diff_back_3d
     module procedure ocean_obc_enhance_diff_back_2d
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
  character(len=32),  dimension(max_obc,max_prog_tracers):: fieldname_tracer
  character(len=256), dimension(max_obc,max_prog_tracers):: filename_tracer
  logical                               :: debug_this_module = .FALSE.
  logical                               :: debug_phase_speed = .FALSE.
  integer, parameter                    :: mobcm1 = max_obc-1
  integer, parameter                    :: ntm2   = max_prog_tracers-2
  integer, parameter                    :: tobcm2 = max_obc*(max_prog_tracers-2)
  integer                               :: ne, te

  character(len=128)                    :: errorstring
  
  integer                               :: id_obc
  
  type(domain2d),allocatable, save      :: obc_diag_domain(:)
  type(domain2d),allocatable, save      :: obc_tracer_domain(:)
  
  data (name(ne), ne=1,max_obc)          /'test_obc', mobcm1*'none'/
  data (fieldname_eta(ne), ne=1,max_obc) /'eta_t', mobcm1*'none'/
  data (filename_eta(ne), ne=1,max_obc)  /'obc_eta_t.nc', mobcm1*'none'/
  data (filename_ud(ne), ne=1,max_obc)  /'obc_ud.nc', mobcm1*'none'/
  data (fieldname_ud(ne), ne=1,max_obc) /'ud', mobcm1*'none'/  
  data ((fieldname_tracer(ne,te),  ne=1, max_obc), te=1, max_prog_tracers) &
           /max_obc*'temp_obc', max_obc*'salt_obc', tobcm2*'none'/
  data ((filename_tracer(ne,te),   ne=1, max_obc), te=1, max_prog_tracers) &
           /max_obc*'INPUT/obc_tr.nc', max_obc*'INPUT/obc_tr.nc' ,tobcm2*'none'/

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
  integer              :: nt  = 0                         ! number of tracers
  real, allocatable    :: wrk(:,:,:)      ! needed for output of phase speed and other quantities
  real, allocatable    :: wrk1(:,:)       ! needed for restart of phase speed
  real, allocatable    :: wrk2(:)         ! needed for enhanced diffusion
  real, allocatable    :: wrk3(:)         ! needed for enhanced diffusion

contains

  !#######################################################################
  !  <SUBROUTINE NAME="ocean_obc_init" >
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
  subroutine ocean_obc_init(have_obc, Time_steps, Domain, Grid, Ocean_options, ver_coordinate, debug)
  !</PUBLICROUTINE>

    logical, intent(inout)                       :: have_obc
    type(ocean_time_steps_type), intent(in)      :: Time_steps
    type (ocean_domain_type),intent(in),  target :: Domain
    type(ocean_grid_type), intent(inout), target :: Grid
    type(ocean_options_type), intent(inout)      :: Ocean_options
    logical, intent(in), optional                :: debug
    integer, intent(in)                          :: ver_coordinate
    !--- some local variables ------------------------------------------
    integer :: m, n, nn, i, j, unit, io_status, ierr, ioff, joff
    integer :: ni, nj, nsize
    integer :: west_at_pe, south_at_pe, east_at_pe, north_at_pe
    integer :: irig_s, ilef_s , jbou_s
    integer :: irig_n, ilef_n , jbou_n
    integer :: jlow_w, jup_w, ibou_w
    integer :: jlow_e, jup_e, ibou_e
    integer :: pe
    character(len=5)     :: pe_name
    
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
       call mpp_error(FATAL, '==>Error in ocean_obc_mod (ocean_obc_init): module already initialized')
    endif

    module_is_initialized = .TRUE.

    Dom => Domain
    Grd => Grid
    call mpp_get_global_domain(Domain%Domain2d, isg, ieg, jsg, jeg)

    vert_coordinate = ver_coordinate

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

    !--- write out version information ---------------------------------
    call write_version_number( version, tagname )

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
    if(nobc .gt. max_obc) call mpp_error(FATAL,'==>Error in ocean_obc_mod: nml nobc is greater than maximum'// &
         ' allowed bounds max_obc. Modify nobc or increase "max_obc" at top of ocean_obc_mod' )

    !--- cyclic condition and open boundary condition along zonal direction 
    !--- these two boundary conditions cannot exist at the same place 
    do m = 1, nobc
       if((trim(direction(m)) == 'west' .or. trim(direction(m)) == 'east')  .and. Grid%cyclic_x) &
            call mpp_error(FATAL, "==>Error in ocean_obc_mod: when west or east boundary is open, "//&
                        "cyclic condition in i-direction cannot exist")
       if((trim(direction(m)) == 'south' .or. trim(direction(m)) == 'north')  .and. Grid%cyclic_y) &
            call mpp_error(FATAL, "==>Error in ocean_obc_mod: when south or north boundary is open, "//&
                        "cyclic condition in j-direction cannot exist")
    enddo

    !--- get the domain decomposition -----------------------------------
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
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
    
    !--- Initialize Obc data type --------------------------------------
    allocate(Obc%bound (nobc))
    allocate(obc_tracer_domain(nobc))
    allocate(obc_diag_domain(nobc))
    do m = 1, nobc
       call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Dom%layout &
                               , obc_tracer_domain(m), maskmap=Dom%maskmap ) 
       call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Dom%layout &
                               , obc_diag_domain(m), maskmap=Dom%maskmap ) 
    enddo

    Obc%bound%is = -999; Obc%bound%ie = -1000    ! dummy number to indicate not on the boundary
    Obc%bound%js = -999; Obc%bound%je = -1000
    Obc%bound%on_bound = .FALSE.
    Obc%bound%report   = .FALSE.
    Obc%bound%obc_relax_eta = .false.
    Obc%bound%obc_relax_eta_profile = .false.

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

       ! Set the vertical mixing coefficients OBC
       select case(trim(obc_mix(m)))
       case('NOTHIN')
          Obc%bound(m)%bcond_mix = NOTHIN
       case('CLAMPD')
          Obc%bound(m)%bcond_mix = CLAMPD
       case('NOGRAD')
          Obc%bound(m)%bcond_mix = NOGRAD
       case('INGRAD')
          Obc%bound(m)%bcond_mix = INGRAD
       case default
          call mpp_error(FATAL,'each element of nml obc_mix should be NOTHIN, CLAMPD, NOGRAD or INGRAD')
       end select

       Obc%bound(m)%ctrop_max             = ctrop_max(m)
       Obc%bound(m)%ctrop_min             = ctrop_min(m)
       Obc%bound(m)%ctrop_inc             = ctrop_inc(m)
       Obc%bound(m)%ctrop_smooth          = ctrop_smooth(m)
       Obc%bound(m)%enh_pnts              = enh_pnts(m)
       Obc%bound(m)%enh_fac_v             = enh_fac_v(m)
       Obc%bound(m)%enh_fac_d             = enh_fac_d(m)
       Obc%bound(m)%obc_consider_convu    = obc_consider_convu(m)
       Obc%bound(m)%obc_adjust_forcing_bt = obc_adjust_forcing_bt(m)
       Obc%bound(m)%obc_vert_advel_t      = obc_vert_advel_t(m)
       Obc%bound(m)%obc_vert_advel_u      = obc_vert_advel_u(m)
       Obc%bound(m)%obc_damp_newton       = obc_damp_newton(m)
       Obc%bound(m)%damp_factor           = damp_factor(m)
       select case(trim(obc_enhance_visc_back(m)))
       case('NONE')
         Obc%bound(m)%obc_enhance_visc_back = NONE
       case('MODERATE')
         Obc%bound(m)%obc_enhance_visc_back = MODERATE
       case('LARGE')
         Obc%bound(m)%obc_enhance_visc_back = LARGE
       case default
          call mpp_error(FATAL,'each element of obc_enhance_visc_back should be NONE, MODERATE or LARGE')
       end select
       select case(trim(obc_enhance_diff_back(m)))
       case('NONE')
         Obc%bound(m)%obc_enhance_diff_back = NONE
       case('MODERATE')
         Obc%bound(m)%obc_enhance_diff_back = MODERATE
       case('LARGE')
         Obc%bound(m)%obc_enhance_diff_back = LARGE
       case default
          call mpp_error(FATAL,'each element of obc_enhance_diff_back should be NONE, MODERATE or LARGE')
       end select
       
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
             do j = jsc, jec                  !kk changed from data to compute domain
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
                     call mpp_open( obc_out_unit(m), 'obc_'//trim(Obc%bound(m)%name)//'.out', action=MPP_OVERWR, &
                                  threading=MPP_MULTI, fileset=MPP_MULTI, nohdrs=.TRUE. )  
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
                allocate(Obc%bound(m)%ctrop(Obc%bound(m)%js:Obc%bound(m)%je))
                allocate(Obc%bound(m)%pgrad(Obc%bound(m)%js:Obc%bound(m)%je))
                Obc%bound(m)%ctrop    = 0.
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
                  call mpp_error(FATAL, "ocean_obc_mod: start of eta data > start of boundary")
                if (Obc%bound(m)%jere < je(m)) &
                  call mpp_error(FATAL, "ocean_obc_mod: end of eta data < end of boundary")
                if (itrs(m) == -999) then 
                  Obc%bound(m)%itrs = is(m) 
                else 
                  Obc%bound(m)%itrs = itrs(m) 
                endif 
                if (itre(m) == -999) then 
                  Obc%bound(m)%itre = ie(m) 
                else 
                  Obc%bound(m)%itre = itre(m) 
                endif 
                if (jtrs(m) == -999) then 
                  Obc%bound(m)%jtrs = js(m) 
                else 
                  Obc%bound(m)%jtrs = jtrs(m) 
                endif 
                if (jtre(m) == -999) then 
                  Obc%bound(m)%jtre = je(m) 
                else 
                  Obc%bound(m)%jtre = jtre(m) 
                endif 
                if (Obc%bound(m)%jtrs > js(m)) &
                  call mpp_error(FATAL, "ocean_obc_mod: start of tracer data > start of boundary")
                if (Obc%bound(m)%jtre < je(m)) &
                  call mpp_error(FATAL, "ocean_obc_mod: end of tracerdata < end of boundary")
             endif
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
             do i = isc, iec                  !kk changed from data to compute domain
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
                    call mpp_open( obc_out_unit(m), 'obc_'//trim(Obc%bound(m)%name)//'.out', action=MPP_OVERWR, &
                                   threading=MPP_MULTI, fileset=MPP_MULTI, nohdrs=.TRUE. )  
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
                allocate(Obc%bound(m)%ctrop(Obc%bound(m)%is:Obc%bound(m)%ie))
                allocate(Obc%bound(m)%pgrad(Obc%bound(m)%is:Obc%bound(m)%ie))
                Obc%bound(m)%ctrop    = 0.
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
                  call mpp_error(FATAL, "ocean_obc_mod: start of eta data > start of boundary")
                if (Obc%bound(m)%iere < ie(m)) &
                  call mpp_error(FATAL, "ocean_obc_mod: end of eta data < end of boundary")
                if (itrs(m) == -999) then 
                  Obc%bound(m)%itrs = is(m) 
                else 
                  Obc%bound(m)%itrs = itrs(m) 
                endif 
                if (itre(m) == -999) then 
                  Obc%bound(m)%itre = ie(m) 
                else 
                  Obc%bound(m)%itre = itre(m) 
                endif 
                if (jtrs(m) == -999) then 
                  Obc%bound(m)%jtrs = js(m) 
                else 
                  Obc%bound(m)%jtrs = jtrs(m) 
                endif 
                if (jtre(m) == -999) then 
                  Obc%bound(m)%jtre = je(m) 
                else 
                  Obc%bound(m)%jtre = jtre(m) 
                endif 
                if (Obc%bound(m)%itrs > is(m)) &
                  call mpp_error(FATAL, "ocean_obc_mod: start of tracer data > start of boundary")
                if (Obc%bound(m)%itre < ie(m)) &
                  call mpp_error(FATAL, "ocean_obc_mod: end of tracerdata < end of boundary")
             endif
          endif
       endif
! Check if region of enhancement is at one task
       if ((Obc%bound(m)%obc_enhance_visc_back/=NONE  .or.   &
            Obc%bound(m)%obc_enhance_diff_back/=NONE) .and.  &
            Obc%bound(m)%on_bound ) then
         if(Obc%bound(m)%direction == WEST) then
           if (Obc%bound(m)%is+Obc%bound(m)%enh_pnts > iec) &
             call mpp_error(FATAL, "ocean_obc_mod: is+enh_pnts exceeds iec")
         endif
         if(Obc%bound(m)%direction == EAST) then
           if (Obc%bound(m)%ie-Obc%bound(m)%enh_pnts < isc) &
             call mpp_error(FATAL, "ocean_obc_mod: ie-enh_pnts lower than isc")
         endif
         if(Obc%bound(m)%direction == SOUTH) then
           if (Obc%bound(m)%js+Obc%bound(m)%enh_pnts > jec) &
             call mpp_error(FATAL, "ocean_obc_mod: js+enh_pnts exceeds jec")
         endif
         if(Obc%bound(m)%direction == NORTH) then
           if (Obc%bound(m)%je-Obc%bound(m)%enh_pnts < jsc) &
             call mpp_error(FATAL, "ocean_obc_mod: je-enh_pnts lower than jsc")
         endif
         if (Obc%bound(m)%obc_enhance_diff_back==MODERATE .and. Obc%bound(m)%enh_fac_d < 1. ) then
             call mpp_error(NOTE, &
         "ocean_obc_mod: obc_enhance_diff_back=MODERATE but enh_fac_d < 1 gives no enhancement of diffusivity. Use enh_fac_d > 1")
         endif
         if (Obc%bound(m)%obc_enhance_diff_back==LARGE .and. Obc%bound(m)%enh_fac_d > 1. ) then
             call mpp_error(NOTE, &
             "ocean_obc_mod: obc_enhance_diff_back = LARGE and enh_fac_d > 1 exceeds CFL-criterion. Use enh_fac_d <= 1")
         endif
         if (Obc%bound(m)%obc_enhance_visc_back==MODERATE .and. Obc%bound(m)%enh_fac_v < 1. ) then
             call mpp_error(NOTE, &
             "ocean_obc_mod: obc_enhance_visc_back = MODERATE but enh_fac_v < 1 gives no enhancement of viscosity.")
         endif
         if (Obc%bound(m)%obc_enhance_visc_back==LARGE .and. Obc%bound(m)%enh_fac_v > 1. ) then
             call mpp_error(NOTE, &
             "ocean_obc_mod: obc_enhance_visc_back = LARGE and enh_fac_v > 1 exceeds CFL-criterion. Use enh_fac_d <= 1")
         endif
       endif

! Read data for relaxation into the compute area only        
    enddo


! Now check for corner points, which are located at two boundaries 
! (south-west, north-west, north-east, south-east)    
! This code needs to be generalised later. Currently, only one
! boundary for each direction is possible. In complex coastal geometry
! this may be to restrictive, but should be sufficient for the beginning ... 
    west_at_pe = 0
    south_at_pe= 0
    east_at_pe = 0
    north_at_pe= 0
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
      allocate(Obc%bound(m)%ones(nsize))
      allocate(Obc%bound(m)%dum(nsize))
      allocate(Obc%bound(m)%dumv(nsize, nk))
      Obc%bound(m)%ones = 1
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
        allocate(Obc%bound(m)%ico(nsize))
        allocate(Obc%bound(m)%work(Obc%bound(m)%js:Obc%bound(m)%je))
        allocate(Obc%bound(m)%csur(Obc%bound(m)%js:Obc%bound(m)%je))
        if(i .ne. Obc%bound(m)%ie) &
         write(obc_out_unit(m),*) 'WARNING: Start and end i locations differ for WEST or EAST obc. is=',i,' : ie=',Obc%bound(m)%ie
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          Obc%bound(m)%iloc(n) = i
          Obc%bound(m)%jloc(n) = j
          Obc%bound(m)%oi1(n) = i + Obc%bound(m)%imap
          Obc%bound(m)%oi2(n) = i + 2 * Obc%bound(m)%imap
          if(Obc%bound(m)%direction==WEST) then
             Obc%bound(m)%ico(n) = i - Obc%bound(m)%imap
          else
             Obc%bound(m)%ico(n) = i
          endif
          n = n + 1
        enddo
        Obc%bound(m)%nloc = n - 1
        ! Set pointers
        Obc%bound(m)%nic => Obc%bound(m)%iloc
        Obc%bound(m)%tic => Obc%bound(m)%jloc
        Obc%bound(m)%oj1 => Obc%bound(m)%jloc
        Obc%bound(m)%oj2 => Obc%bound(m)%jloc
        Obc%bound(m)%d1i => Obc%bound(m)%ones
        Obc%bound(m)%d1j => Obc%bound(m)%jloc
        Obc%bound(m)%icf => Obc%bound(m)%iloc
        if(Obc%bound(m)%direction==EAST) Obc%bound(m)%icf => Obc%bound(m)%oi1
        Obc%bound(m)%jcf => Obc%bound(m)%jloc
        Obc%bound(m)%jco => Obc%bound(m)%jloc
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
        allocate(Obc%bound(m)%jco(nsize))
        allocate(Obc%bound(m)%work(Obc%bound(m)%is:Obc%bound(m)%ie))
        allocate(Obc%bound(m)%csur(Obc%bound(m)%is:Obc%bound(m)%ie))
        if(j .ne. Obc%bound(m)%je) &
          write(obc_out_unit(m),*)        &
             'WARNING : Boundary',n,' start and end j locations differ for NORTH or SOUTH obc. js=' &
             ,j,' : je=',Obc%bound(m)%je
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          Obc%bound(m)%iloc(n) = i
          Obc%bound(m)%jloc(n) = j
          Obc%bound(m)%oj1(n) = j + Obc%bound(m)%jmap
          Obc%bound(m)%oj2(n) = j + 2 * Obc%bound(m)%jmap
          if(Obc%bound(m)%direction==SOUTH) then
             Obc%bound(m)%jco(n) = j - Obc%bound(m)%jmap
          else
             Obc%bound(m)%jco(n) = j
          endif
          n = n + 1
        enddo
        Obc%bound(m)%nloc = n - 1
        ! Set pointers
        Obc%bound(m)%nic => Obc%bound(m)%jloc
        Obc%bound(m)%tic => Obc%bound(m)%iloc
        Obc%bound(m)%oi1 => Obc%bound(m)%iloc
        Obc%bound(m)%oi2 => Obc%bound(m)%iloc
        Obc%bound(m)%d1i => Obc%bound(m)%iloc
        Obc%bound(m)%d1j => Obc%bound(m)%ones
        Obc%bound(m)%icf => Obc%bound(m)%iloc
        Obc%bound(m)%jcf => Obc%bound(m)%jloc
        if(Obc%bound(m)%direction==NORTH) Obc%bound(m)%jcf => Obc%bound(m)%oj1
        Obc%bound(m)%ico => Obc%bound(m)%iloc
        n = 1
      endif
      ! Initialise the gravity wave speed
      do n = 1, Obc%bound(m)%nloc
         i = Obc%bound(m)%oi1(n)
         j = Obc%bound(m)%oj1(n)
         nn = Obc%bound(m)%tic(n)
         Obc%bound(m)%csur(nn)  = sqrt(grav*Grd%ht(i,j))
      enddo

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
    allocate(wrk1(isd:ied,jsd:jed))
    wrk1 = 0.
    
    write(stdoutunit,*) '-----------------------------------------------------------------'
    write(stdoutunit,*) 'The following setup for OBC has been found:'
    write(stdoutunit,*) 'Total number of OBC: ', nobc
    do m = 1, nobc
       write(stdoutunit,*) ' Setup of OBC ',m,', ',Obc%bound(m)%name,':'
       write(stdoutunit,*) ' direction : ', trim(direction(m))
       write(stdoutunit,*) ' ==> See files obc_'//trim(Obc%bound(m)%name)//'.out.* for more information.'
       write(stdoutunit,*) '     Enable debug_this_module for details on domain decomposition.'
       
       if(.not. Obc%bound(m)%on_bound) cycle
       
       write(obc_out_unit(m),*) 'Setup of OBC ',m,', ',Obc%bound(m)%name,':'
       write(obc_out_unit(m),*) ' direction     : ', trim(direction(m))
       write(obc_out_unit(m),*) ' points        : ', Obc%bound(m)%np
! if debugging is on each PE has an output file open otherwise this info is useless
       if(debug_this_module) then
         write(obc_out_unit(m),*) ' pe - name          :', trim(pe_name)
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
       if (Obc%bound(m)%obc_vert_advel_t) then 
          write(obc_out_unit(m),*) ' with vertical tracer advection'
       else  
          write(obc_out_unit(m),*) ' no vertical tracer advection'
       endif
       if (Obc%bound(m)%obc_vert_advel_u) then 
          write(obc_out_unit(m),*) ' with vertical momentum advection'
       else  
          write(obc_out_unit(m),*) ' no vertical momentum advection'
       endif
       select case(Obc%bound(m)%obc_enhance_visc_back)
       case( NONE )
          write(obc_out_unit(m),*) ' no enhanced background viscosity '
       case( MODERATE )
          write(obc_out_unit(m),*) ' enhance background viscosity by ',Obc%bound(m)%enh_fac_v
          write(obc_out_unit(m),*) ' over ',Obc%bound(m)%enh_pnts,' points.'
       case( LARGE )
          write(obc_out_unit(m),*) ' enhance viscosity to ',Obc%bound(m)%enh_fac_v ,&
                            ' of maximum CFL viscosity'
          write(obc_out_unit(m),*) ' over ',Obc%bound(m)%enh_pnts,' points.'
       end select
       select case(Obc%bound(m)%obc_enhance_diff_back)
       case( NONE )
          write(obc_out_unit(m),*) ' no enhanced background mixing '
       case( MODERATE )
          write(obc_out_unit(m),*) ' enhance background mixing by ',Obc%bound(m)%enh_fac_d
          write(obc_out_unit(m),*) ' over ',Obc%bound(m)%enh_pnts,' points.'
       case( LARGE )
          write(obc_out_unit(m),*) ' enhance mixing to ',Obc%bound(m)%enh_fac_d, &
                            ' of maximum CFL diffusivity'
          write(obc_out_unit(m),*) ' over ',Obc%bound(m)%enh_pnts,' points.'
       end select
       if (Obc%bound(m)%obc_damp_newton) then
          write(obc_out_unit(m),*) ' apply Newtonian damping of ',Obc%bound(m)%damp_factor
       else
          write(obc_out_unit(m),*) ' no Newtonian damping '
       endif
       if (Obc%bound(m)%obc_adjust_forcing_bt) then
          write(obc_out_unit(m),*) ' remove baroclinic cross boundary pressure gradients from'
          write(obc_out_unit(m),*) ' vertically integrated forcing'
       else  
          write(obc_out_unit(m),*) ' do not remove baroclinic cross boundary pressure gradients from'
          write(obc_out_unit(m),*) ' vertically integrated forcing'
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

    call ocean_obc_check_topog(Grid%kmt)
    
    call ocean_obc_set_mask

    call mpp_clock_end(id_obc)
    
    return

  end subroutine ocean_obc_init
  !  </SUBROUTINE> NAME="ocean_obc_init"


  !#######################################################################
  !  <SUBROUTINE NAME="ocean_obc_tracer_init" >
  !   <DESCRIPTION>
  !      Allocates space and initializes all stuff for tracers at OBC
  !   </DESCRIPTION>
  !   <IN NAME="debug" TYPE="logical"></IN>
  !<PUBLICROUTINE>
  subroutine ocean_obc_tracer_init(Time, T_prog, num_prog_tracers, debug)
    !</PUBLICROUTINE>
    type(ocean_time_type), intent(in)                   :: Time
    type(ocean_prog_tracer_type), intent(inout), target :: T_prog(:)
    integer, intent(in)                                 :: num_prog_tracers
    logical, intent(in),        optional                :: debug
    integer :: fld_size(4)

    !--- some local variables ------------------------------------------
    integer :: m, n, taum1, tau, taup1, unit, ierr, io_status
    integer :: is, ie, js, je
    integer :: is_u, ie_u, js_u, je_u
    
  integer :: stdoutunit 
  stdoutunit=stdout() 

    if (PRESENT(debug)) debug_this_module = (debug.or.debug_this_module)

    !   Allocate space to store tracer flux through open boundary       
    !--- get the number of prog tracers
    nt = num_prog_tracers

    !--- read namelist again, to get information on tracers, now nt is known

    unit = open_namelist_file()
    read  (unit, ocean_obc_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_obc_nml')
    call close_file (unit)
         
    if(debug_this_module) write(stdoutunit,*) 'Using ', nt,' tracers in OBC'
    
    if(nt .gt. max_prog_tracers) call mpp_error(FATAL,'==>Error in ocean_obc_mod: num_prog_tracers is greater than '//&
     'maximum number of prog tracers allocated in obc. For more tracers, increase "max_prog_tracers" at top of ocean_obc_mod')

    if(nt == 0) then 
       call mpp_error(NOTE,'==>NOTE in ocean_obc_mod: number of prognostics tracers is 0')
       return
    endif

    call mpp_clock_begin(id_obc)    

    do n=1, nt
       allocate(T_prog(n)%otf(nobc))
    enddo
    do m = 1, nobc
       if(.not. Obc%bound(m)%on_bound) cycle

       ! Set the tracer OBC
       allocate(Obc%bound(m)%bcond_tra(nt))
       do n=1, nt
          Obc%bound(m)%bcond_tra(n) = 0
          if (index(trim(obc_tra(m,n)),'NOTHIN') /= 0) &
               Obc%bound(m)%bcond_tra(n) = Obc%bound(m)%bcond_tra(n) + NOTHIN
          if (index(trim(obc_tra(m,n)),'FILEIN') /= 0) &
               Obc%bound(m)%bcond_tra(n) = Obc%bound(m)%bcond_tra(n) + FILEIN
          if (index(trim(obc_tra(m,n)),'IOW') /= 0)    &
               Obc%bound(m)%bcond_tra(n) = Obc%bound(m)%bcond_tra(n) + IOW
          if (index(trim(obc_tra(m,n)),'NOGRAD') /= 0) &
               Obc%bound(m)%bcond_tra(n) = Obc%bound(m)%bcond_tra(n) + NOGRAD
          if (index(trim(obc_tra(m,n)),'UPSTRM') /= 0) &
               Obc%bound(m)%bcond_tra(n) = Obc%bound(m)%bcond_tra(n) + UPSTRM
       enddo

       select case(Obc%bound(m)%direction )
       case(WEST)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%js:Obc%bound(m)%je,nk))
          enddo
       case(EAST)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%js:Obc%bound(m)%je,nk))
          enddo
       case(SOUTH)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%is:Obc%bound(m)%ie,nk))
          enddo
       case(NORTH)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%is:Obc%bound(m)%ie,nk))
          enddo
       end select
    enddo
    
    write(stdoutunit,*) '-----------------------------------------------------------------'
    write(stdoutunit,*) 'The following setup for time interpolation of OBC data has been found:'

    allocate(Obc%tracer(nobc,nt))
    do m = 1, nobc
       if (Obc%bound(m)%on_bound) then
          !         data needed for tracer
          allocate(Obc%bound(m)%obc_relax_tracer(nt)) 
          allocate(Obc%bound(m)%obc_flow_relax(nt)) 
          allocate(Obc%bound(m)%obc_consider_sources(nt)) 
          allocate(Obc%bound(m)%obc_tracer_no_inflow(nt)) 
          do n=1, nt
             Obc%bound(m)%obc_relax_tracer(n)      = obc_relax_tracer(m,n)
             Obc%bound(m)%obc_flow_relax(n)        = obc_flow_relax(m,n)
             Obc%bound(m)%obc_consider_sources(n)  = obc_consider_sources(m,n)
             Obc%bound(m)%obc_tracer_no_inflow(n)  = obc_tracer_no_inflow(m,n)
             Obc%tracer(m,n)%file               = trim(filename_tracer(m,n))
             Obc%tracer(m,n)%field              = trim(fieldname_tracer(m,n))
             Obc%tracer(m,n)%rel_coef_in        = rel_coef_tracer_in(m,n)
             Obc%tracer(m,n)%rel_coef_out       = rel_coef_tracer_out(m,n)
             Obc%tracer(m,n)%rel_pnts           = rel_clin_pnts(m,n)
             write(obc_out_unit(m),*) trim(Obc%tracer(m,n)%file)
             write(obc_out_unit(m),*) trim(Obc%tracer(m,n)%field)
          enddo

          is = min(Obc%bound(m)%is, Obc%bound(m)%ie + 2 * Obc%bound(m)%imap)
          ie = max(Obc%bound(m)%is, Obc%bound(m)%ie + 2 * Obc%bound(m)%imap)
          js = min(Obc%bound(m)%js, Obc%bound(m)%je + 2 * Obc%bound(m)%jmap)
          je = max(Obc%bound(m)%js, Obc%bound(m)%je + 2 * Obc%bound(m)%jmap)

! This part of the domain definition is the same for tracers and eta
!          call mpp_set_compute_domain(obc_eta_domain(m),&
!                                      Obc%bound(m)%is, Obc%bound(m)%ie, &
!                                      Obc%bound(m)%js, Obc%bound(m)%je, &
!                                      Obc%bound(m)%ie-Obc%bound(m)%is+1, Obc%bound(m)%je-Obc%bound(m)%js+1)
!          call mpp_set_data_domain   (obc_eta_domain(m),&
!                                      Obc%bound(m)%is, Obc%bound(m)%ie, &
!                                      Obc%bound(m)%js, Obc%bound(m)%je, &
!                                      Obc%bound(m)%ie-Obc%bound(m)%is+1, Obc%bound(m)%je-Obc%bound(m)%js+1,&
!                                       .false., .false.)
          call mpp_set_compute_domain(obc_tracer_domain(m),&
                                      Obc%bound(m)%is, Obc%bound(m)%ie, &
                                      Obc%bound(m)%js, Obc%bound(m)%je, &
                                      Obc%bound(m)%ie-Obc%bound(m)%is+1, Obc%bound(m)%je-Obc%bound(m)%js+1)
          call mpp_set_data_domain   (obc_tracer_domain(m),&
                                      Obc%bound(m)%is, Obc%bound(m)%ie, &
                                      Obc%bound(m)%js, Obc%bound(m)%je, &
                                      Obc%bound(m)%ie-Obc%bound(m)%is+1, Obc%bound(m)%je-Obc%bound(m)%js+1,&
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

     
             !           allocate variables for time interpolation of external data
             !           if relaxation or upstream advection is switched on
          do n=1, nt
             if(Obc%bound(m)%obc_relax_tracer(n).or. &
                  iand(Obc%bound(m)%bcond_tra(n), FILEIN) == FILEIN) then
                if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
                   allocate(Obc%tracer(m,n)%data(1,Obc%bound(m)%js:Obc%bound(m)%je,nk))
                else
                   allocate(Obc%tracer(m,n)%data(Obc%bound(m)%is:Obc%bound(m)%ie,1,nk))
                endif
! Set the global domain part for the tracer input file
                call mpp_set_global_domain( obc_tracer_domain(m), &
                                      xbegin=Obc%bound(m)%itrs, ybegin=Obc%bound(m)%jtrs, &
                                      xend=Obc%bound(m)%itre, yend=Obc%bound(m)%jtre, &
                                      xsize=Obc%bound(m)%itre-Obc%bound(m)%itrs+1, &
                                      ysize=Obc%bound(m)%jtre-Obc%bound(m)%jtrs+1)   
                Obc%tracer(m,n)%data = 0
                Obc%tracer(m,n)%id = init_external_field(  &
                     Obc%tracer(m,n)%file,Obc%tracer(m,n)%field, &
                     domain=obc_tracer_domain(m),        &
                     verbose=debug_this_module)
                fld_size(:) = get_external_field_size(Obc%tracer(m,n)%id)
                if(fld_size(2) .ne. size(Obc%tracer(m,n)%data,2)) then
                   write(errorstring,'(I3,I3)') fld_size(2), size(Obc%tracer(m,n)%data,2) 
                   call mpp_error(FATAL,trim(errorstring)//'invalid dimension 2 of input field in '//trim(Obc%tracer(m,n)%file))
                endif
                if(fld_size(3) .ne. size(Obc%tracer(m,n)%data,3)) then
                   write(errorstring,'(I3,I3)') fld_size(3), size(Obc%tracer(m,n)%data,3) 
                   call mpp_error(FATAL,trim(errorstring)//'invalid dimension 3 of input field in '//trim(Obc%tracer(m,n)%file))
                endif
             endif
          enddo

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
          do n=1, nt
             write(obc_out_unit(m),'(a,I2.2,1x, a)') ' tracer  : ', n, trim(T_prog(n)%name)
             write(obc_out_unit(m),*) 'Tracer OBC : ', trim(obc_tra(m,n))
             if (Obc%bound(m)%obc_tracer_no_inflow(n)) then
                write(obc_out_unit(m),*) '  Outflow only allowed for this boundary'
             else
                write(obc_out_unit(m),*) '  Inflow and outflow allowed for this boundary'
             endif
             if (Obc%bound(m)%obc_consider_sources(n)) then
                write(obc_out_unit(m),*) '  Source and SGS terms are valid at the boundary.'
             else
                write(obc_out_unit(m),*) '  Source and SGS terms are not valid at the boundary.'
             endif
             if (Obc%bound(m)%obc_flow_relax(n) > 1) then
                write(obc_out_unit(m),*) '  Flow relaxation zone = ', Obc%bound(m)%obc_flow_relax(n)
             endif
             if (Obc%bound(m)%obc_relax_tracer(n)) then
                write(obc_out_unit(m),*) '  relax tracer  :', Obc%bound(m)%obc_relax_tracer(n)
                write(obc_out_unit(m),*) '  relax tracer: rel_pnts rel_coef_in rel_coef_out input-file       input-name'
                write(obc_out_unit(m),'(16x,I4,2x,F8.5,2x,F8.5,2x,a,2x,a)') &
                     Obc%tracer(m,n)%rel_pnts,&
                     Obc%tracer(m,n)%rel_coef_in, Obc%tracer(m,n)%rel_coef_out,& 
                     trim(obc%tracer(m,n)%file), trim(obc%tracer(m,n)%field)
             else
                write(obc_out_unit(m),*) '  no tracer relaxation '
             endif
          enddo
       endif
    enddo

!  Set tracer poutside the boundaries to defined values to avoid artificial diffusion   
    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1
    
    do n = 1, nt
      call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taum1), 'T')
      call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,tau), 'T')
      call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taup1), 'T')
    enddo
        
! prepare diagnostics    
    allocate(wrk(isd:ied,jsd:jed,nk))
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
         is_u = Obc%bound(m)%isg + Obc%bound(m)%imap
         ie_u = Obc%bound(m)%ieg + Obc%bound(m)%imap
         js_u = Obc%bound(m)%jsg + Obc%bound(m)%jmap
         je_u = Obc%bound(m)%jeg + Obc%bound(m)%jmap
         Obc%bound(m)%id_xu = diag_axis_init ('xu_'//trim(Obc%bound(m)%name),Grd%grid_x_u(is_u:ie_u), &
              'degrees_E','x','ucell longitude',set_name='ocean', Domain2=obc_diag_domain(m), aux='geolon_u')
         Obc%bound(m)%id_yu = diag_axis_init ('yu_'//trim(Obc%bound(m)%name),Grd%grid_y_u(js_u:je_u), &
              'degrees_N','y','ucell latitude',set_name='ocean', Domain2=obc_diag_domain(m), aux='geolat_u')
         Obc%bound(m)%id_ctrop = register_diag_field ('ocean_model', 'ctrop_p_'//trim(Obc%bound(m)%name),                  &
                          (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yt/), Time%model_time, 'barotr phase speed on open bounds',&
                         'm/s', missing_value=missing_value, range=(/-1.e8,1.e8/))   
         if (Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
            Obc%bound(m)%id_transport = register_diag_field ('ocean_model', 'obc_transport_int_z_'//trim(Obc%bound(m)%name),&
                                 (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt/),Time%model_time,                                &
                                 'z-integral of horizontal transport on obc',                                               &
                                 'm^2/s', missing_value=missing_value, range=(/-1.e18,1.e18/))
         else
            Obc%bound(m)%id_transport = register_diag_field ('ocean_model', 'obc_transport_int_z_'//trim(Obc%bound(m)%name),&
                                 (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu/),Time%model_time,                                &
                                 'z-integral of horizontal transport on obc',                                               &
                                 'm^2/s', missing_value=missing_value, range=(/-1.e18,1.e18/))
         endif

         allocate(Obc%bound(m)%id_tracer_flux(nt))
         allocate(Obc%bound(m)%id_tracer_flux_dif(nt))
         allocate(Obc%bound(m)%id_tracer_flux_adv(nt))
         allocate(Obc%bound(m)%id_cclin(nt))
         allocate(Obc%bound(m)%id_tracer_data(nt))
         do n = 1, nt
           Obc%bound(m)%id_cclin(n) = register_diag_field & 
                ('ocean_model', 'cclin_'//trim(T_prog(n)%name)//'_'//trim(Obc%bound(m)%name), &
                (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/), Time%model_time, 'phase speed on obc',&
                'm/s', missing_value=missing_value, range=(/-1.e18,1.e18/))
           Obc%bound(m)%id_tracer_data(n) = register_diag_field &
                ('ocean_model', 'tracer_data_'//trim(T_prog(n)%name)//'_'//trim(Obc%bound(m)%name),                        &
                (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/), Time%model_time, 'tracer data on obc',&
                trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e18,1.e18/))
           if (Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
             if(trim(T_prog(n)%name) == 'temp') then
               Obc%bound(m)%id_tracer_flux(n) = register_diag_field &
                    ('ocean_model', 'obc_heat_flux_'//trim(Obc%bound(m)%name),          &
                    (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/),&
                    Time%model_time, 'total horizontal tracer flux on obc',         &
                    'Watt/m', missing_value=missing_value, range=(/-1.e18,1.e18/))
               
               Obc%bound(m)%id_tracer_flux_dif(n) = register_diag_field &
                    ('ocean_model', 'obc_heat_flux_dif_'//trim(Obc%bound(m)%name),          &
                    (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/),&
                    Time%model_time, 'horizontal tracer diffusion flux on obc',         &
                    'Watt/m', missing_value=missing_value, range=(/-1.e18,1.e18/))

               Obc%bound(m)%id_tracer_flux_adv(n) = register_diag_field &
                    ('ocean_model', 'obc_heat_flux_adv_'//trim(Obc%bound(m)%name),          &
                    (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/),&
                    Time%model_time, 'horizontal tracer advection flux on obc',         &
                    'Watt/m', missing_value=missing_value, range=(/-1.e18,1.e18/))
             else
               Obc%bound(m)%id_tracer_flux(n) = register_diag_field & 
                    ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_flux_'//trim(Obc%bound(m)%name), &
                    (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/),             &
                    Time%model_time, 'total horizontal tracer flux on obc',                     &
                    ' kg/m/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
               
               Obc%bound(m)%id_tracer_flux_dif(n) = register_diag_field & 
                    ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_flux_dif_'//trim(Obc%bound(m)%name), &
                    (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/),             &
                    Time%model_time, 'horizontal tracer diffusion flux on obc',                     &
                    ' kg/m/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
               
               Obc%bound(m)%id_tracer_flux_adv(n) = register_diag_field & 
                    ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_flux_adv_'//trim(Obc%bound(m)%name), &
                    (/Obc%bound(m)%id_xu, Obc%bound(m)%id_yt,Grd%tracer_axes_depth(3)/),             &
                    Time%model_time, 'horizontal tracer advection flux on obc',                     &
                    ' kg/m/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
             endif
           else
             if(trim(T_prog(n)%name) == 'temp') then
               Obc%bound(m)%id_tracer_flux(n) = register_diag_field &
                    ('ocean_model', 'obc_heat_flux_'//trim(Obc%bound(m)%name),           &
                    (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu,Grd%tracer_axes_depth(3)/), &
                    Time%model_time, 'totel horizontal tracer flux on obc',         &
                    'Watt/m', missing_value=missing_value, range=(/-1.e18,1.e18/))
               
               Obc%bound(m)%id_tracer_flux_dif(n) = register_diag_field &
                    ('ocean_model', 'obc_heat_flux_dif_'//trim(Obc%bound(m)%name),           &
                    (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu,Grd%tracer_axes_depth(3)/), &
                    Time%model_time, 'horizontal tracer diffusion flux on obc',         &
                    'Watt/m', missing_value=missing_value, range=(/-1.e18,1.e18/))
               
               Obc%bound(m)%id_tracer_flux_adv(n) = register_diag_field &
                    ('ocean_model', 'obc_heat_flux_adv_'//trim(Obc%bound(m)%name),           &
                    (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu,Grd%tracer_axes_depth(3)/), &
                    Time%model_time, 'horizontal tracer advection flux on obc',         &
                    'Watt/m', missing_value=missing_value, range=(/-1.e18,1.e18/))
             else
                Obc%bound(m)%id_tracer_flux(n) = register_diag_field & 
                     ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_flux_'//trim(Obc%bound(m)%name), &
                     (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu,Grd%tracer_axes_depth(3)/),             &
                     Time%model_time, 'total horizontal tracer flux on obc',                     &
                     ' kg/m/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

                Obc%bound(m)%id_tracer_flux_dif(n) = register_diag_field & 
                     ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_flux_dif_'//trim(Obc%bound(m)%name), &
                     (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu,Grd%tracer_axes_depth(3)/),             &
                     Time%model_time, 'horizontal tracer diffusion flux on obc',                     &
                     ' kg/m/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

                Obc%bound(m)%id_tracer_flux_adv(n) = register_diag_field & 
                     ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_flux_adv_'//trim(Obc%bound(m)%name), &
                     (/Obc%bound(m)%id_xt, Obc%bound(m)%id_yu,Grd%tracer_axes_depth(3)/),             &
                     Time%model_time, 'horizontal tracer advection flux on obc',                     &
                     ' kg/m/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
             endif
           endif
         enddo
       endif
    enddo   
    
    call mpp_clock_end(id_obc)
    return

  end subroutine ocean_obc_tracer_init
  !  </SUBROUTINE> NAME="ocean_obc_tracer_init"


  !#####################################################################
  !  <SUBROUTINE NAME="ocean_obc_prepare">
  !   <DESCRIPTION>
  !      Prepares OBC  
  !      
  !   </DESCRIPTION>
  !<PUBLICROUTINE>
  subroutine ocean_obc_prepare(Time, T_prog)
  !</PUBLICROUTINE>
    type(ocean_time_type), intent(in)             :: Time 
    type(ocean_prog_tracer_type), intent(inout)   :: T_prog(:)

    integer                                       :: m, n

    call mpp_clock_begin(id_obc)

!   prepare the new data needed for time interpolation of external data.

    call ocean_obc_prepare_baro(Time)
    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

       do n=1, nt
         if(Obc%bound(m)%obc_relax_tracer(n).or. &
              iand(Obc%bound(m)%bcond_tra(n), FILEIN) == FILEIN) then
           call time_interp_external(Obc%tracer(m,n)%id,Time%model_time,Obc%tracer(m,n)%data,verbose=debug_this_module)
         endif
       enddo

       !--- set the tracer_flux to 0
       do n =1, nt
         T_prog(n)%otf(m)%flux(:,:) = 0.0
       enddo
   
    enddo

    call mpp_clock_end(id_obc)

    return

  end subroutine ocean_obc_prepare
  !  </SUBROUTINE> NAME="ocean_obc_prepare"


  !<SUBROUTINE> NAME="ocean_obc_surface_height"
  !<PUBLICROUTINE>
  subroutine ocean_obc_surface_height(Time, Ext_mode)
  !</PUBLICROUTINE>
    type(ocean_time_type), intent(in)             :: Time
    type(ocean_external_mode_type), intent(inout) :: Ext_mode
    integer                                       :: taum1, tau, taup1
    logical                                       :: used

    integer :: i, j, m, n, nn

  integer :: stdoutunit
  stdoutunit=stdout()

    call mpp_clock_begin(id_obc)

    taum1=Time%taum1;tau=Time%tau;taup1=Time%taup1

! The best choice is to set eta_t to eta_t_bar. Updating with a radiation condition
! introduces much noise.

    do m = 1, nobc

       ! if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

       do n = 1, Obc%bound(m)%nloc
          i = Obc%bound(m)%iloc(n)   ! i boundary location
          j = Obc%bound(m)%jloc(n)   ! j boundary location
          Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t_bar(i,j,tau)
       enddo
    enddo

    ! Set the corner cells
    if ( south_west_corner ) then
       Ext_mode%eta_t(i_sw,j_sw,taup1) = ( Ext_mode%eta_t(i_sw+1,j_sw+1,taup1) &
                                         + Ext_mode%eta_t(i_sw+1,j_sw,  taup1) &
                                         + Ext_mode%eta_t(i_sw,  j_sw+1,taup1) ) /3.
    endif
    if ( south_east_corner ) then
       Ext_mode%eta_t(i_se,j_se,taup1) = ( Ext_mode%eta_t(i_se-1,j_se+1,taup1) &
                                         + Ext_mode%eta_t(i_se-1,j_se,  taup1) &
                                         + Ext_mode%eta_t(i_se,  j_se+1,taup1) ) /3.
    endif
    if ( north_west_corner ) then
       Ext_mode%eta_t(i_nw,j_nw,taup1) = ( Ext_mode%eta_t(i_nw+1,j_nw-1,taup1) &
                                         + Ext_mode%eta_t(i_nw+1,j_nw,  taup1) &
                                         + Ext_mode%eta_t(i_nw,  j_nw-1,taup1) ) /3.
    endif
    if ( north_east_corner ) then
       Ext_mode%eta_t(i_ne,j_ne,taup1) = ( Ext_mode%eta_t(i_ne-1,j_ne-1,taup1) &
                                         + Ext_mode%eta_t(i_ne-1,j_ne,  taup1) &
                                         + Ext_mode%eta_t(i_ne,  j_ne-1,taup1) ) /3.
    endif

    ! Update eta at global halo point to make the gradient accross boundary is 0
    call ocean_obc_update_boundary(Ext_mode%eta_t(:,:,taup1), 'T')

! Write some diagnostics for quantities calculated at barotropic time steps. Hence, here
! only the last value of a barotropic sequence of the average can be written.
    do m = 1, nobc
       if(.not. Obc%bound(m)%on_bound) cycle
       if(Obc%bound(m)%id_ctrop > 0) then
         wrk = 0.
         do n = 1, Obc%bound(m)%nloc
        i = Obc%bound(m)%iloc(n)   ! i boundary location
        j = Obc%bound(m)%jloc(n)   ! j boundary location
        nn = Obc%bound(m)%tic(n)   ! Index tangential to boundary
            wrk(i,j,1) = Obc%bound(m)%ctrop(nn)
     enddo
         used = send_data(Obc%bound(m)%id_ctrop, &
                      wrk(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1),  &
                          Time%model_time,                               &
              rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1))
       endif

    enddo

    if(debug_this_module) then
       call write_chksum_2d('After ocean_obc_barotropic_dteta, eta', Ext_mode%eta_t(COMP,taup1))
    endif

    call mpp_clock_end(id_obc)

    return

  end subroutine ocean_obc_surface_height
  !</SUBROUTINE> NAME="ocean_obc_surface_height"


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
  !update vertical viscosity and diffusivity OBC's
  !<SUBROUTINE NAME="ocean_obc_mixing" >
  !   <INOUT NAME="visc_cbu" TYPE="real, dimension(isd:,jsd:,:)"></INOUT>
  !   <INOUT NAME="diff_cbt" TYPE="real, dimension(isd:,jsd:,:,2)"></INOUT>

  !<PUBLICROUTINE INTERFACE="ocean_obc_mixing" >
  subroutine ocean_obc_mixing(visc_cbt, diff_cbt, field1, field2)
    !</PUBLICROUTINE>

    real, dimension(isd:ied,jsd:jed,nk),   intent(inout) :: visc_cbt
    real, dimension(isd:ied,jsd:jed,nk),   intent(inout) :: diff_cbt
    real, dimension(isd:ied,jsd:jed,nk), optional, intent(inout) :: field1, field2

    integer :: i, j, k, m, n, i1, j1
!    real :: tm

!    tm_scale_to_secs('10 days',tm)

    call mpp_clock_begin(id_obc)

    do m = 1, nobc
       if(iand(Obc%bound(m)%bcond_mix, NOTHIN) == NOTHIN) then
          ! Update eta at global halo point to make a no-gradient
          call ocean_obc_update_boundary(visc_cbt, 'T')
          call ocean_obc_update_boundary(diff_cbt, 'T')
          cycle
       endif
       if(iand(Obc%bound(m)%bcond_mix, NOGRAD) == NOGRAD) then
          ! Set a no-gradient
          do k= 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
                visc_cbt(i,j,k) = visc_cbt(i1,j1,k)
                diff_cbt(i,j,k) = diff_cbt(i1,j1,k)
             enddo
          enddo
          ! Update eta at global halo point to make a no-gradient
          call ocean_obc_update_boundary(visc_cbt, 'T')
          call ocean_obc_update_boundary(diff_cbt, 'T')
          if (PRESENT(field1)) then
            do k= 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
                field1(i,j,k) = field1(i1,j1,k)
              enddo
            enddo
            call ocean_obc_update_boundary(field1, 'T')
         endif
         if (PRESENT(field2)) then
            do k= 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
                field2(i,j,k) = field2(i1,j1,k)
              enddo
            enddo
            call ocean_obc_update_boundary(field2, 'T')
         endif
       else if(iand(Obc%bound(m)%bcond_mix, INGRAD) == INGRAD) then
          ! Set an interior gradient
          do k= 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi2(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj2(n)   ! j interior cell location
                visc_cbt(i,j,k) = visc_cbt(i1,j1,k)
                diff_cbt(i,j,k) = diff_cbt(i1,j1,k)
             enddo
          enddo
          ! Update eta at global halo point to make a no-gradient
          call ocean_obc_update_boundary(visc_cbt, 'T', 'i')
          call ocean_obc_update_boundary(diff_cbt, 'T', 'i')
          if (PRESENT(field1)) then
            do k= 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi2(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj2(n)   ! j interior cell location
                field1(i,j,k) = field1(i1,j1,k)
              enddo
            enddo
            call ocean_obc_update_boundary(field1, 'T')
         endif
         if (PRESENT(field2)) then
            do k= 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi2(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj2(n)   ! j interior cell location
                field2(i,j,k) = field2(i1,j1,k)
              enddo
            enddo
            call ocean_obc_update_boundary(field2, 'T')
         endif
       else if(iand(Obc%bound(m)%bcond_mix, CLAMPD) == CLAMPD) then
          ! Set to zero. No vertical diffusion is performed in this case.
          do k= 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                visc_cbt(i,j,k) = 0.0
                diff_cbt(i,j,k) = 0.0
             enddo
          enddo
          ! Update eta at global halo point to make a no-gradient
          call ocean_obc_update_boundary(visc_cbt, 'T', 'z')
          call ocean_obc_update_boundary(diff_cbt, 'T', 'z')
          if (PRESENT(field1)) then
            do k= 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                field1(i,j,k) = 0.0
              enddo
            enddo
            call ocean_obc_update_boundary(field1, 'T')
         endif
         if (PRESENT(field2)) then
            do k= 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                field2(i,j,k) = 0.0
              enddo
            enddo
            call ocean_obc_update_boundary(field2, 'T')
         endif
       endif
    enddo

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_mixing
  !</SUBROUTINE> NAME="ocean_obc_mixing"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_adjust_advel">
  !   <DESCRIPTION>
  !      Subtract wrong vertical bottom velocity 
  !   </DESCRIPTION>
  !   <INOUT NAME="Adv_vel" TYPE="ocean_adv_vel_type" >
  !      Advection velocities 
  !   </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_adjust_advel(Adv_vel)
    !</PUBLICROUTINE>

    type(ocean_adv_vel_type), intent(inout) :: Adv_vel

    integer :: i, j, k, m, n, kmt

    call mpp_clock_begin(id_obc)

    do m = 1, nobc
       ! If obc_vert_advel_t is true (default) calculate a reasonable 
       ! approximation to Adv_vel%wrho_bt. This will be used for 
       ! Ext_mode%deta_dt and tracer advection as well. Adv_vel%wrho_bu 
       ! will be calculated for baroclinic velocity.
       ! If obc_vert_advel_t is false (from namelist) calculate a 
       ! reasonale approximation for Adv_vel%wrho_bt for Adv_vel%wrho_bu. 
       ! Then set Adv_vel%wrho_bt to zero.

       ! If on current pe there is no point on the bound, then just return
       if(Obc%bound(m)%on_bound) then
          ! Calculate the phase speed at west or east boundary
         if( Obc%bound(m)%obc_vert_advel_t .or. Obc%bound(m)%obc_vert_advel_u) then
           do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
             do i = Obc%bound(m)%isd, Obc%bound(m)%ied
                kmt = Grd%kmt(i,j) 
                do k=0, kmt
                   Adv_vel%wrho_bt(i,j,k) = Adv_vel%wrho_bt(i,j,k) - Adv_vel%wrho_bt(i,j,kmt)   !*float(k)/float(kmt) 
                enddo
             enddo
           enddo
         else
           do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
             do i = Obc%bound(m)%isd, Obc%bound(m)%ied
                Adv_vel%wrho_bt(i,j,:) = 0. 
             enddo
           enddo
         endif
       endif

 !      call mpp_update_domains(Adv_vel%wrho_bt(:,:,:), Dom%domain2d)
 !      call ocean_obc_update_boundary(Adv_vel%wrho_bt(:,:,:), 'T')

       if(.not. Obc%bound(m)%on_bound) cycle

       if( Obc%bound(m)%obc_vert_advel_u ) then
          do n = 1, Obc%bound(m)%nloc
             i = Obc%bound(m)%icf(n)   ! i boundary face location
             j = Obc%bound(m)%jcf(n)   ! j boundary face location

             do k=0, Grd%kmu(i,j)
                Adv_vel%wrho_bu(i,j,k) = (Adv_vel%wrho_bt(i,j,k) * &
                     Grd%dte(i,j)*Grd%dus(i,j) + &
                     Adv_vel%wrho_bt(i+1,j,k)* &
                     Grd%dtw(i+1,j)*Grd%dus(i,j) + &
                     Adv_vel%wrho_bt(i,j+1,k) * &
                     Grd%dte(i,j+1)*Grd%dun(i,j) + &
                     Adv_vel%wrho_bt(i+1,j+1,k) * &
                     Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
             enddo
          enddo
       else
          do n = 1, Obc%bound(m)%nlod
             i = Obc%bound(m)%idf(n)   ! i boundary location
             j = Obc%bound(m)%jdf(n)   ! j boundary location
             Adv_vel%wrho_bu(i,j,:) = 0. 
          enddo
       endif
             
       if( .not. Obc%bound(m)%obc_vert_advel_t .and. Obc%bound(m)%obc_vert_advel_u) then
          do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
             do i = Obc%bound(m)%isd, Obc%bound(m)%ied
                Adv_vel%wrho_bt(i,j,:) = 0. 
             enddo
          enddo
       endif

       call ocean_obc_update_boundary(Adv_vel%wrho_bt(:,:,:), 'T')
       call ocean_obc_update_boundary(Adv_vel%wrho_bu(:,:,:), 'C')

       if(debug_this_module) then
          write(obc_out_unit(m),*) '==> ocean_obc_adjust_advel: '
          write(obc_out_unit(m),*) '    The vertical bottom velocity has been removed for boundary',m,'!'
       endif
    enddo

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_adjust_advel
  !</SUBROUTINE> NAME="ocean_obc_adjust_advel"

  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_adjust_forcing_bt">
  !   <DESCRIPTION>
  !      Add wrong pressure gradient
  !   </DESCRIPTION>
  !   <INOUT NAME="Ext_mode" TYPE="ocean_external_mode_type" > </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_adjust_forcing_bt(Ext_mode)
    !</PUBLICROUTINE>
    type(ocean_external_mode_type), intent(inout) :: Ext_mode

    integer :: i, j, m, n, nn, np
    integer, dimension(:), pointer :: ip
    integer, dimension(:), pointer :: jp

    call mpp_clock_begin(id_obc)

    do m = 1, nobc

       if(.not. Obc%bound(m)%on_bound) cycle
       if(.not. Obc%bound(m)%obc_adjust_forcing_bt) cycle

       ! Set pointers
       np = 1;
       ip => Obc%bound(m)%iloc
       jp => Obc%bound(m)%jloc
       if (Obc%bound(m)%direction == NORTH .or. &
            Obc%bound(m)%direction == SOUTH) np = 2
       if (Obc%bound(m)%direction == EAST .or. &
            Obc%bound(m)%direction == NORTH) then
          ip => Obc%bound(m)%oi1
          jp => Obc%bound(m)%oj1
       endif

       do n = 1, Obc%bound(m)%nloc
          i = ip(n)                  ! i boundary location
          j = jp(n)                  ! j boundary location
          nn = Obc%bound(m)%tic(n)   ! Index tangential to boundary

          Ext_mode%forcing_bt(i,j,np) = Ext_mode%forcing_bt(i,j,np) + &
               Obc%bound(m)%pgrad(nn)

       enddo

       if(debug_this_module) then
          write(obc_out_unit(m),*) '==> ocean_obc_adjust_forcing_bt: '
          write(obc_out_unit(m),*) '    The baroclinic pressure gradient across boundary',m
          write(obc_out_unit(m),*) '    has been removed from the free surface mode!'
       endif
    enddo

    call mpp_clock_end(id_obc)

    return

  end subroutine ocean_obc_adjust_forcing_bt
  !</SUBROUTINE> NAME="ocean_obc_adjust_forcing_bt"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_enhance_diff_back_3d">
  !   <DESCRIPTION>
  !      enhance diffusion near open boundary
  !   </DESCRIPTION>
  !   <INOUT NAME="diff_cet" TYPE="real array 3D" > </INOUT >
  !   <INOUT NAME="diff_cnt" TYPE="real array 3D" > </INOUT >
  !<PUBLICROUTINE>
  subroutine ocean_obc_enhance_diff_back_3d(diff_cet, diff_cnt, scheme)
    !</PUBLICROUTINE>
    real, dimension(isd:,jsd:,1:), intent(inout) :: diff_cet, diff_cnt
    character(len=3), optional, intent(in) :: scheme
    integer :: i, j, m, n, nn, ii, jj
    real :: fact, diff_max, sf
    character(len=3) :: turb_scheme

    call mpp_clock_begin(id_obc)

    turb_scheme = 'lap'
    if (PRESENT(scheme)) turb_scheme = trim(scheme)

    do m = 1, nobc

       if(.not. Obc%bound(m)%on_bound) cycle

       select case(Obc%bound(m)%obc_enhance_diff_back)
       case( NONE )
         cycle
       case( MODERATE )
         sf = Obc%bound(m)%enh_fac_d
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
           ii = i
           jj = j
           fact = 1.0
           do nn = 1, Obc%bound(m)%enh_pnts
             diff_cet(ii,jj,:) = diff_cet(ii,jj,:) * (1.0 + sf / fact)
             diff_cnt(ii,jj,:) = diff_cnt(ii,jj,:) * (1.0 + sf / fact)
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
             fact = fact + 1.0
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_diff_back: '
           write(obc_out_unit(m),*) '    Diffusion is moderately enhanced near boundary ',&
                             trim(Obc%bound(m)%name),'.'
         endif
       case( LARGE )
         sf = Obc%bound(m)%enh_fac_d
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
          ! Get the maximum allowable viscosity
           diff_max = 1.0 / (Grd%dxu(i,j) * Grd%dxu(i,j)) + &
                      1.0 / (Grd%dyu(i,j) * Grd%dyu(i,j))
           if ( turb_scheme == 'bih' ) diff_max = 2.*diff_max**2
           diff_max = sf * (1.0 / (4.0 * diff_max * dtime_t))
          ! Get the viscosity in the interior
           ii = i + Obc%bound(m)%imap * Obc%bound(m)%enh_pnts
           jj = j + Obc%bound(m)%jmap * Obc%bound(m)%enh_pnts
           wrk2(:) = diff_cet(ii,jj,:)
           wrk3(:) = diff_cnt(ii,jj,:)
          ! Linearly interpolate
           ii = i
           jj = j
           wrk2(:)   = (wrk2(:) - diff_max)  / Obc%bound(m)%enh_pnts
           wrk3(:)   = (wrk3(:) - diff_max)  / Obc%bound(m)%enh_pnts
           do nn = 1, Obc%bound(m)%enh_pnts
             diff_cet(ii,jj,:) = wrk2(:) * nn + diff_max
             diff_cnt(ii,jj,:) = wrk3(:) * nn + diff_max
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_diff_back: '
           write(obc_out_unit(m),*) '    Background diffusion is substantially enhanced near boundary',&
                             trim(Obc%bound(m)%name),'.'
         endif

       end select
    enddo
    call ocean_obc_update_boundary(diff_cet, 'C')
    call ocean_obc_update_boundary(diff_cnt, 'C')

!   do not care about parallel issues, update is done in calling subroutine.

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_enhance_diff_back_3d
  !</SUBROUTINE> NAME="ocean_obc_enhance_diff_back_3d"

  
  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_enhance_diff_back_2d">
  !   <DESCRIPTION>
  !      enhance diffusivity near open boundary
  !   </DESCRIPTION>
  !   <IN NAME="aiso_back" TYPE="real array 2D" > </IN>
  !<PUBLICROUTINE>  
  subroutine ocean_obc_enhance_diff_back_2d(aiso_back, scheme)
    !</PUBLICROUTINE>
    real, dimension(isd:,jsd:), intent(inout) :: aiso_back
    character(len=3), optional, intent(in) :: scheme
    integer :: i, j, m, n, nn, ii, jj
    real :: fact, diff_max, aiso_i, sf, vfaciso
    character(len=3) :: turb_scheme

 
    call mpp_clock_begin(id_obc)

    turb_scheme = 'lap'
    if (PRESENT(scheme)) turb_scheme = trim(scheme)

    do m = 1, nobc

       if(.not. Obc%bound(m)%on_bound) cycle

       select case(Obc%bound(m)%obc_enhance_diff_back)
       case( NONE )
         cycle
       case( MODERATE )
         sf = Obc%bound(m)%enh_fac_d
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
           ii = i
           jj = j
           fact = 1.0
           do nn = 1, Obc%bound(m)%enh_pnts
             aiso_back(ii,jj) = aiso_back(ii,jj) * (1.0 + sf / fact)
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
             fact = fact + 1.0
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_diff_back: '
           write(obc_out_unit(m),*) '    Background diffusion is moderately enhanced near boundary',&
                             trim(Obc%bound(m)%name),'.'
         endif
       case( LARGE )
         sf = Obc%bound(m)%enh_fac_d
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
          ! Get the maximum allowable viscosity
           diff_max = 1.0 / (Grd%dxu(i,j) * Grd%dxu(i,j)) + &
                      1.0 / (Grd%dyu(i,j) * Grd%dyu(i,j))
           if ( turb_scheme == 'bih' ) diff_max = 2.*diff_max**2
           diff_max = sf * (1.0 / (4.0 * diff_max * dtime_t))
          ! Get the viscosity in the interior
           ii = i + Obc%bound(m)%imap * Obc%bound(m)%enh_pnts
           jj = j + Obc%bound(m)%jmap * Obc%bound(m)%enh_pnts
           aiso_i = aiso_back(ii,jj)
          ! Linearly interpolate
           ii = i
           jj = j
           vfaciso   = (aiso_i - diff_max)  / Obc%bound(m)%enh_pnts
           do nn = 1, Obc%bound(m)%enh_pnts
             aiso_back(ii,jj) = vfaciso * nn + diff_max
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_diff_back: '
           write(obc_out_unit(m),*) '    Background diffusion is substantially enhanced near boundary',&
                             trim(Obc%bound(m)%name),'.'
         endif

       end select

       
    enddo      

    call ocean_obc_update_boundary(aiso_back, 'C')    
!   do not care about parallel issues, update is done in calling subroutine.

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_enhance_diff_back_2d
  !</SUBROUTINE> NAME="ocean_obc_enhance_diff_back_2d"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_enhance_visc_back_2d">
  !   <DESCRIPTION>
  !      enhance viscosity near open boundary
  !   </DESCRIPTION>
  !   <INOUT NAME="aiso_back" TYPE="real array 2D" > </INOUT>
  !<PUBLICROUTINE>  
  subroutine ocean_obc_enhance_visc_back_2d(aiso_back, scheme)
    !</PUBLICROUTINE>
    real, dimension(isd:,jsd:), intent(inout) :: aiso_back
    character(len=3),  optional, intent(in) :: scheme
    integer :: i, j, m, n, nn, ii, jj
    real :: fact, visc_max, aiso_i, sf, vfaciso
    character(len=3) :: turb_scheme
 
    call mpp_clock_begin(id_obc)

    turb_scheme = 'lap'
    if (PRESENT(scheme)) turb_scheme = trim(scheme)

    do m = 1, nobc

       if(.not. Obc%bound(m)%on_bound) cycle
        
       select case(Obc%bound(m)%obc_enhance_visc_back)
       case( NONE )
         cycle
       case( MODERATE )
         sf = Obc%bound(m)%enh_fac_v
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
           ii = i
           jj = j
           fact = 1.0
           do nn = 1, Obc%bound(m)%enh_pnts
             aiso_back(ii,jj) = aiso_back(ii,jj) * (1.0 + sf / fact)
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
             fact = fact + 1.0
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_visc_back: '
           write(obc_out_unit(m),*) '    Background diffusion is moderately enhanced near boundary',&
                             trim(Obc%bound(m)%name),'.'
         endif
       case( LARGE )
         sf = Obc%bound(m)%enh_fac_v
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
          ! Get the maximum allowable viscosity
           visc_max = 1.0 / (Grd%dxu(i,j) * Grd%dxu(i,j)) + &
                      1.0 / (Grd%dyu(i,j) * Grd%dyu(i,j))
           if ( turb_scheme == 'bih' ) visc_max = 2.*visc_max**2
           visc_max = sf * (1.0 / (4.0 * visc_max * dtime_t))
          ! Get the viscosity in the interior
           ii = i + Obc%bound(m)%imap * Obc%bound(m)%enh_pnts
           jj = j + Obc%bound(m)%jmap * Obc%bound(m)%enh_pnts
           aiso_i = aiso_back(ii,jj)
          ! Linearly interpolate
           ii = i
           jj = j
           vfaciso   = (aiso_i - visc_max)  / Obc%bound(m)%enh_pnts
           do nn = 1, Obc%bound(m)%enh_pnts
             aiso_back(ii,jj) = vfaciso * nn + visc_max
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_visc_back: '
           write(obc_out_unit(m),*) '    Background diffusion is substantially enhanced near boundary',&
                               trim(Obc%bound(m)%name),'.'
         endif
       end select
    enddo      

    call ocean_obc_update_boundary(aiso_back, 'C')    
!   do not care about parallel issues, update is done in calling subroutine.

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_enhance_visc_back_2d
  !</SUBROUTINE> NAME="ocean_obc_enhance_visc_back_2d"

  
  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_enhance_visc_back_3d">
  !   <DESCRIPTION>
  !      enhance viscosity near open boundary. Maximum viscosity for
  !      stability is set on the boundary, linearly decreasing to the
  !      interior value at enh_pnts into the interior.
  !   </DESCRIPTION>
  !   <INOUT NAME="aiso_back" TYPE="real array 3D" > </INOUT>
  !   <INOUT NAME="aaniso_back" TYPE="real array 3D" > </INOUT>
  !<PUBLICROUTINE>  
  subroutine ocean_obc_enhance_visc_back_3d(aiso_back, aaniso_back, scheme)
    !</PUBLICROUTINE>
    real, dimension(isd:,jsd:,1:), intent(inout) :: aiso_back, aaniso_back
    character(len=3),  optional, intent(in) :: scheme
    integer :: i, j, m, n, nn, ii, jj
    real :: fact, visc_max, sf
    character(len=3) :: turb_scheme

    call mpp_clock_begin(id_obc)

    turb_scheme = 'lap'
    if (PRESENT(scheme)) turb_scheme = trim(scheme)

    do m = 1, nobc

       if(.not. Obc%bound(m)%on_bound) cycle

       select case(Obc%bound(m)%obc_enhance_visc_back)
       case( NONE )
         cycle
       case( MODERATE )
         sf = Obc%bound(m)%enh_fac_v
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
           ii = i
           jj = j
           fact = 1.0
           do nn = 1, Obc%bound(m)%enh_pnts
             aiso_back(ii,jj,:)   = aiso_back  (ii,jj,:) * (1.0 + sf / fact)
             aaniso_back(ii,jj,:) = aaniso_back(ii,jj,:) * (1.0 + sf / fact)
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
             fact = fact + 1.0
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_visc_back: '
           write(obc_out_unit(m),*) '    Background diffusion is moderately enhanced near boundary',&
                             trim(Obc%bound(m)%name),'.'
         endif
       case( LARGE )
         sf = Obc%bound(m)%enh_fac_v
         do n = 1, Obc%bound(m)%nlod
           i = Obc%bound(m)%idf(n)   ! i boundary location
           j = Obc%bound(m)%jdf(n)   ! j boundary location
          ! Get the maximum allowable viscosity
           visc_max = 1.0 / (Grd%dxu(i,j) * Grd%dxu(i,j)) + &
                      1.0 / (Grd%dyu(i,j) * Grd%dyu(i,j))
           if ( turb_scheme == 'bih' ) visc_max = 2.*visc_max**2
           visc_max = sf * (1.0 / (4.0 * visc_max * dtime_t))
          ! Get the viscosity in the interior
           ii = i + Obc%bound(m)%imap * Obc%bound(m)%enh_pnts
           jj = j + Obc%bound(m)%jmap * Obc%bound(m)%enh_pnts
           wrk2(:) = aiso_back  (ii,jj,:)
           wrk3(:) = aaniso_back(ii,jj,:)
          ! Linearly interpolate
           ii = i
           jj = j
           wrk2(:) = (wrk2(:) - visc_max) / Obc%bound(m)%enh_pnts
           wrk3(:) = (wrk3(:) - visc_max) / Obc%bound(m)%enh_pnts
           do nn = 1, Obc%bound(m)%enh_pnts
             aiso_back  (ii,jj,:) = wrk2(:) * nn + visc_max
             aaniso_back(ii,jj,:) = wrk3(:) * nn + visc_max
             ii = ii + Obc%bound(m)%imap
             jj = jj + Obc%bound(m)%jmap
           enddo
         enddo
         if(debug_this_module) then
           write(obc_out_unit(m),*) '==> ocean_obc_enhance_visc_back: '
           write(obc_out_unit(m),*) '    Background diffusion is substantially enhanced near boundary',&
                               trim(Obc%bound(m)%name),'.'
         endif
       end select
    enddo      
    call ocean_obc_update_boundary(aiso_back, 'C')    
    call ocean_obc_update_boundary(aaniso_back, 'C')    

!   do not care about parallel issues, update is done in calling subroutine.

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_enhance_visc_back_3d
  !</SUBROUTINE>NAME="ocean_obc_enhance_visc_back_3d"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_tracer">
  !   <DESCRIPTION>
  !      Extrapolate tracer on the open boundaries for ocean model and regional atmosphere model.
  !   </DESCRIPTION>
  !   <IN NAME="rho_dztr" TYPE="real, dimension(isc:,jsc:,:)"> 
  !     contains Thickness%rho_dztr from update_tracer
  !   </IN>
  !   <IN NAME="taum1, tau, taup1" TYPE="integer" >
  !     time step index
  !   </IN>
  !   <IN NAME="time" TYPE="type(time_type)">
  !     model time
  !   </IN>
  !   <IN NAME="n" TYPE="integer">
  !     tracer number
  !   </IN>
  !   <INOUT NAME="tracer" TYPE="real, dimension(isd:,jsd:,:,:)"> 
  !      Tracer field
  !   </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_tracer(tracer, adv_vet, adv_vnt, Thickness, pme, taum1, tau, taup1, time, tn)
    !</PUBLICROUTINE>
    real, dimension(isd:,jsd:,:,:), intent(inout) :: tracer            ! tracer
    real, dimension(isd:,jsd:,:), target, intent(in) :: adv_vet           ! advection velocity on east face of t-cell
    real, dimension(isd:,jsd:,:), target, intent(in) :: adv_vnt           ! advection velocity on north face of t-cell
    type(ocean_thickness_type), intent(in)           :: Thickness
    real, dimension(isd:,jsd:),        intent(in) :: pme               ! pme
    integer,                           intent(in) :: taum1, tau, taup1 ! time step index
    type(ocean_time_type),                   intent(in) :: time              ! model time
    integer, intent(in)                           :: tn                ! only when n=1, the max phase speed will be calculated

    !--- local variables -----------------------------------------------
    integer :: it, iu, m
    integer :: k, i, j, tlevel, istrt, iend, ii, jj, dbg
    integer :: i1, i2, j1, j2, n, nn, im, jm, id, jd, ic, jc
    real    :: cgrid, var, cmax, uout, uin, rel_var, adv_obc, adv, alpha, sign
    real, dimension(:,:), pointer :: cclin
    real, dimension(:,:), pointer :: dg
    real, dimension(:,:), pointer :: ddt
    real, dimension(:,:), pointer :: dur
    real, dimension(:,:,:), pointer :: adv_vel
    logical                         :: used

  integer :: stdoutunit 
  stdoutunit=stdout() 

    tlevel = taum1

    call mpp_clock_begin(id_obc)

    ! Loop through the boundaries
    do m = 1, nobc

       ! If on current pe there is no point on the bound, then just return
       if(.not.Obc%bound(m)%on_bound) cycle

       if(debug_this_module) then
          write(obc_out_unit(m),*) 'Doing update of tracer ',tn,' at obc ', trim(Obc%bound(m)%name)
       endif

       !----------------------------------------------------------------
       ! Set pointers
       !----------------------------------------------------------------
       if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
          dg => Grd%dxu
          ddt => Grd%dyte
          dur => Grd%dxur
          adv_vel => adv_vet
       else
          dg => Grd%dyu
          ddt => Grd%dxtn
          dur => Grd%dyur
          adv_vel => adv_vnt
       endif
       sign = Obc%bound(m)%sign
       cclin => Obc%bound(m)%dumv
       cclin = 0.0;
       dbg = 0

       !----------------------------------------------------------------
       ! Get the phase speed for tracers if required
       !----------------------------------------------------------------
       if(iand(Obc%bound(m)%bcond_tra(tn), ORLANS) == ORLANS) then
          do n = 1, Obc%bound(m)%nloc

             i = Obc%bound(m)%iloc(n)  ! i boundary location
             j = Obc%bound(m)%jloc(n)  ! j boundary location
             i1 = Obc%bound(m)%oi1(n)  ! i interior cell location
             j1 = Obc%bound(m)%oj1(n)  ! j interior cell location
             i2 = Obc%bound(m)%oi2(n)  ! i 2x interior cell location
             j2 = Obc%bound(m)%oj2(n)  ! j 2x interior cell location

             cgrid = dg(i1,j1)/dtime_t
             cmax  = - sign*cgrid

             do k = 1, nk
                var = sign *(tracer(i1,j1,k,taup1) + tracer(i1,j1,k,taum1) - &
                     2.0 * tracer(i2,j2,k,tau))
                if(abs(var) .lt. small) then
                   cclin(n,k) = cmax
                else
                   cclin(n,k) = -(tracer(i1,j1,k,taum1) - &
                                  tracer(i1,j1,k,taup1)) * cgrid / var
                   if(cclin(n,k)*sign .gt. 0.0) cclin(n,k) = 0.
                endif
                if(Obc%bound(m)%direction == WEST .or. &
                   Obc%bound(m)%direction == SOUTH)  then
                   cclin(n,k) = max(cmax,cclin(n,k))* Grd%tmask(i1,j1,k)* Grd%tmask(i2,j2,k)
                else
                   cclin(n,k) = min(cmax,cclin(n,k))* Grd%tmask(i1,j1,k)* Grd%tmask(i2,j2,k)
                endif
             enddo
          enddo
          dbg = 1
       endif

       if(iand(Obc%bound(m)%bcond_tra(tn), IOW) == IOW) then
          do n = 1, Obc%bound(m)%nloc
             i = Obc%bound(m)%iloc(n)   ! i boundary location
             j = Obc%bound(m)%jloc(n)   ! j boundary location
             i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
             j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
             i2 = Obc%bound(m)%oi2(n)   ! i 2x interior cell location
             j2 = Obc%bound(m)%oj2(n)   ! j 2x interior cell location

             cgrid = dg(i1,j1)/dtime_t
             cmax  = - sign*cgrid

             do k = 1, nk
                var = sign * (tracer(i2,j2,k,tau)-tracer(i1,j1,k,tau) )
                if(abs(var) .lt. small) then
                   cclin(n,k) = cmax
                else
                   cclin(n,k) = -(tracer(i1,j1,k,taup1) - &
                                  tracer(i1,j1,k,taum1)) * cgrid / var
                   if(cclin(n,k)*sign .gt. 0.0) cclin(n,k) = 0.
                endif
                if(Obc%bound(m)%direction == WEST .or. &
                   Obc%bound(m)%direction == SOUTH)  then
                   cclin(n,k) = max(cmax,cclin(n,k))* Grd%tmask(i1,j1,k)* Grd%tmask(i2,j2,k)
                else
                   cclin(n,k) = min(cmax,cclin(n,k))* Grd%tmask(i1,j1,k)* Grd%tmask(i2,j2,k)
                endif
             enddo
          enddo
          dbg = 1
       endif

       ! Phase speed debugging information
       if(Obc%bound(m)%id_cclin(tn) > 0) then
         wrk = 0.
         do n = 1, Obc%bound(m)%nloc
            i = Obc%bound(m)%iloc(n)   ! i boundary location
            j = Obc%bound(m)%jloc(n)   ! j boundary location
            wrk(i,j,:) = cclin(n,:)
         enddo
         used = send_data(Obc%bound(m)%id_cclin(tn), &
              wrk(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1:nk),  &
              Time%model_time,           &
              rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1:nk))
       endif
       if(debug_this_module .and. dbg == 1) then
          k = 1
          if(Obc%bound(m)%direction == WEST .or. &
             Obc%bound(m)%direction == EAST) then
             if(Obc%bound(m)%direction == WEST) then
                j = Obc%bound(m)%je
             else
                j = Obc%bound(m)%js
             endif
             write(obc_out_unit(m),*) 'k= 1, j = ',j
             write(obc_out_unit(m),*) 'cgrid   = ',cgrid
             write(obc_out_unit(m),*) 'cmax    = ',cmax
!kk             write(obc_out_unit(m),*) 'ctrop   = ',Obc%bound(m)%ctrop(j)
             write(obc_out_unit(m),*) 'clin    = ',cclin(1,k)
          else
             if(Obc%bound(m)%direction == SOUTH) then
                i = Obc%bound(m)%ie
             else
                i = Obc%bound(m)%is
             endif
             write(obc_out_unit(m),*) 'k= 1, j = ',j
             write(obc_out_unit(m),*) 'cgrid   = ',cgrid
             write(obc_out_unit(m),*) 'cmax    = ',cmax
!kk             write(obc_out_unit(m),*) 'ctrop   = ',Obc%bound(m)%ctrop(i)
             write(obc_out_unit(m),*) 'clin    = ',cclin(1,k)
          endif
       endif

       !----------------------------------------------------------------
       ! Perform the open boundary condition
       !----------------------------------------------------------------
       ! No-gradient
       if(iand(Obc%bound(m)%bcond_tra(tn), NOGRAD) == NOGRAD) then
         if(iand(Obc%bound(m)%bcond_tra(tn), UPSTRM) == UPSTRM) then
           do k= 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
                
                Obc%bound(m)%work2(i,j,k) = tracer(i1,j1,k,taup1)
             enddo
           enddo
         else
           do k= 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
                
                tracer(i,j,k,taup1) = tracer(i1,j1,k,taup1)
             enddo
           enddo
         endif
          ! File input
       else if(iand(Obc%bound(m)%bcond_tra(tn), FILEIN) == FILEIN) then
         if(iand(Obc%bound(m)%bcond_tra(tn), UPSTRM) == UPSTRM) then
           do k= 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                id = Obc%bound(m)%d1i(n)
                jd = Obc%bound(m)%d1j(n)

                Obc%bound(m)%work2(i,j,k) = Obc%tracer(m,tn)%data(id,jd,k)
             enddo
           enddo
         endif
       endif
       ! Upstream advection
       ! orlanski   : upwind advection for outflow, no advection for inflow
       ! with relax : upwind advection in any case, 
       !              for inflow with prescribed tracer
       ! at a N/E boundary:  for north/eastward flow, uout > 0, uin = 0
       !                     for south/westward flow, uout = 0, uin < 0
       ! at a S/W boundary:  for north/eastward flow, uout = 0, uin > 0
       !                     for south/westward flow, uout < 0, uin = 0
       if(iand(Obc%bound(m)%bcond_tra(tn), UPSTRM) == UPSTRM) then
          if (Obc%bound(m)%obc_consider_sources(tn)) then
             k = 1
! start all OBC updates from taup1. 
! source terms and SGS terms from the normal scheme take effect.            
             tlevel = taup1   
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                
                tracer(i,j,k,taup1) = tracer(i,j,k,taup1) &
                     + tracer(i,j,k,taum1)*Thickness%rho_dztr(i,j,k)* Grd%tmask(i,j,k) & 
                       *(Thickness%rho_dzt(i,j,k,taup1) -  Thickness%rho_dzt(i,j,k,taum1) &
                       -(Thickness%mass_source(i,j,k) + pme(i,j)) * dtime_t) 
             enddo
             if (vert_coordinate .ne. GEOPOTENTIAL) then
               do k= 2, nk
                 do n = 1, Obc%bound(m)%nloc
                   i = Obc%bound(m)%iloc(n)   ! i boundary location
                   j = Obc%bound(m)%jloc(n)   ! j boundary location
                   
                   tracer(i,j,k,taup1) = tracer(i,j,k,taup1) &
                        + tracer(i,j,k,taum1)*Thickness%rho_dztr(i,j,k)* Grd%tmask(i,j,k) & 
                          *(Thickness%rho_dzt(i,j,k,taup1) -  Thickness%rho_dzt(i,j,k,taum1) &
                          - Thickness%mass_source(i,j,k) * dtime_t) 
                 enddo
               enddo
             endif 
          else
! start all OBC updates from taum1. 
! all changes from the normal scheme are discarded             
             tlevel = taum1
          endif  
          do k = 1, nk
             do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                i1 = Obc%bound(m)%oi1(n)   ! i interior cell location
                j1 = Obc%bound(m)%oj1(n)   ! j interior cell location
                ic = Obc%bound(m)%icf(n)   ! i boundary face location
                jc = Obc%bound(m)%jcf(n)   ! j boundary face location
                im = i - Obc%bound(m)%imap
                jm = j - Obc%bound(m)%jmap
                
                adv = adv_vel(ic,jc,k) * Thickness%rho_dztr(i,j,k) ! The velocity in m/s
                uout  = 0.5 * (adv - sign * abs(adv))
                uin   = 0.5 * (adv + sign * abs(adv))
                if(Obc%bound(m)%obc_tracer_no_inflow(tn)) uin = 0.0


!!                adv_obc = - sign * Grd%datr(i,j) * &
!!                     (uin * ( Obc%bound(m)%work2(i,j,k)* ddt(im,jm) - &
!!                     tracer(i,j,k,tau)* ddt(i,j)) + &
!!                     uout* (tracer(i,j,k,tau)* ddt(i,j) - &
!!                     tracer(i1,j1,k,tau) * ddt(i1,j1)))
                adv_obc = - sign * Grd%datr(i,j) * ddt(i,j) * &
                     (uin * ( Obc%bound(m)%work2(i,j,k) - tracer(i,j,k,tau  ))  &
                    + uout* ( tracer(i,j,k,tau)         - tracer(i1,j1,k,tau)))

                ! Update the tracer
!                tracer(i,j,k,taup1) = Obc%bound(m)%work2(i,j,k) + dtime_t * ( &
                tracer(i,j,k,taup1) = tracer(i,j,k,tlevel) + dtime_t * ( &
                     - adv_obc ) 

                ! Add the wave-like contribution implicitly
                var = sign * dtime_t * dur(i,j) * cclin(n,k)
                tracer(i,j,k,taup1) = (tracer(i,j,k,taup1) - &
                                    tracer(i1,j1,k,taup1) * var) / (1. - var) &
                                    * Grd%tmask(i,j,k)

                ! Relax the updated solution to prescribed data
                ! W/S boundary:   adv + c > 0 -> inward flow, strong relaxation
                !                 adv + c < 0 -> outward flow, weak relaxation
                ! E/N boundary:   adv + c < 0 -> inward flow, strong relaxation
                !                 adv + c > 0 -> outward flow, weak relaxation
                if (Obc%bound(m)%obc_relax_tracer(tn) ) then
                   id = Obc%bound(m)%d1i(n)
                   jd = Obc%bound(m)%d1j(n)
                   if (sign*(adv + cclin(n,k)) > 0.) then
                      rel_var = dtime_t*Obc%tracer(m,tn)%rel_coef_in * Grd%tmask(i,j,k)
                   else
                      rel_var = dtime_t*Obc%tracer(m,tn)%rel_coef_out* Grd%tmask(i,j,k)
                   endif
                   ! Add the relaxation implicitly
                   ii = i
                   jj = j
                   do nn = 1, Obc%tracer(m,tn)%rel_pnts
                      tracer(ii,jj,k,taup1) = (tracer(ii,jj,k,taup1) + &
                           rel_var * Obc%tracer(m,tn)%data(id,jd,k)) / (1.0 + rel_var)&
                           *Grd%tmask(ii,jj,k)
                      ii = ii + Obc%bound(m)%imap
                      jj = jj + Obc%bound(m)%jmap
                   enddo
                endif
                ! Invoke the flow relaxation scheme of Martinsen and 
                ! Engedahl (1987), eqn 3 if required.
                if (Obc%bound(m)%obc_flow_relax(tn) > 1) then
                   ii = i1
                   jj = j1
                   do nn = 2, Obc%bound(m)%obc_flow_relax(tn)
                      alpha = (1.0 - tanh(0.5 * (nn - 1.0)))* Grd%tmask(i,j,k)
!                     Alternate representation, M&E(1987), eqn 4
!                      alpha = ((Obc%bound(m)%obc_flow_relax(tn) - nn + 1) / &
!                           Obc%bound(m)%obc_flow_relax(tn)) ** 2.0
                      tracer(ii,jj,k,taup1) = (alpha * tracer(i,j,k,taup1) + &
                           (1.0 - alpha) * tracer(ii,jj,k,taup1))*Grd%tmask(ii,jj,k)
                      ii = ii + Obc%bound(m)%imap
                      jj = jj + Obc%bound(m)%jmap
                   enddo
                endif
             enddo
          enddo
       else ! no advection
         if (Obc%bound(m)%obc_relax_tracer(tn) ) then
           do k = 1, nk
              do n = 1, Obc%bound(m)%nloc
                i = Obc%bound(m)%iloc(n)   ! i boundary location
                j = Obc%bound(m)%jloc(n)   ! j boundary location
                id = Obc%bound(m)%d1i(n)
                jd = Obc%bound(m)%d1j(n)
                if (sign*(cclin(n,k)) > 0.) then
                   rel_var = dtime_t*Obc%tracer(m,tn)%rel_coef_in * Grd%tmask(i,j,k)
                else
                   rel_var = dtime_t*Obc%tracer(m,tn)%rel_coef_out* Grd%tmask(i,j,k)
                endif
                ! Add the relaxation implicitly
                ii = i
                jj = j
                do nn = 1, Obc%tracer(m,tn)%rel_pnts
                   tracer(ii,jj,k,taup1) = (tracer(ii,jj,k,taup1) + &
                        rel_var * Obc%tracer(m,tn)%data(id,jd,k)) / (1.0 + rel_var)&
                        *Grd%tmask(ii,jj,k)
                   ii = ii + Obc%bound(m)%imap
                   jj = jj + Obc%bound(m)%jmap
                enddo
             enddo
           enddo
         endif
       endif
       if(iand(Obc%bound(m)%bcond_tra(tn), FILEIN) == FILEIN) then
         if(Obc%bound(m)%id_tracer_data(tn) > 0) then
           wrk = 0.
           do n = 1, Obc%bound(m)%nloc
              i = Obc%bound(m)%iloc(n)   ! i boundary location
              j = Obc%bound(m)%jloc(n)   ! j boundary location
              id = Obc%bound(m)%d1i(n)
              jd = Obc%bound(m)%d1j(n)
              wrk(i,j,:) = Obc%tracer(m,tn)%data(id,jd,:)*Grd%tmask(i,j,:)
           enddo
           used = send_data(Obc%bound(m)%id_tracer_data(tn), &
                wrk(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1:nk),  &
                Time%model_time,                   &
                rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1:nk))
         endif
       endif
 
! manage OBC-points at two boundaries
    enddo
    if ( south_west_corner ) then
       tracer(i_sw,j_sw,:,taup1) = ( tracer(i_sw+1,j_sw+1,:,taup1) &
                                   + tracer(i_sw+1,j_sw,  :,taup1) &
                                   + tracer(i_sw,  j_sw+1,:,taup1) ) /3.
    endif
    if ( south_east_corner ) then
       tracer(i_se,j_se,:,taup1) = ( tracer(i_se-1,j_se+1,:,taup1) &
                                   + tracer(i_se-1,j_se,  :,taup1) &
                                   + tracer(i_se,  j_se+1,:,taup1) ) /3.
    endif
    if ( north_west_corner ) then
       tracer(i_nw,j_nw,:,taup1) = ( tracer(i_nw+1,j_nw-1,:,taup1) &
                                   + tracer(i_nw+1,j_nw,  :,taup1) &
                                   + tracer(i_nw,  j_nw-1,:,taup1) ) /3.
    endif
    if ( north_east_corner ) then
       tracer(i_ne,j_ne,:,taup1) = ( tracer(i_ne-1,j_ne-1,:,taup1) &
                                   + tracer(i_ne-1,j_ne,  :,taup1) &
                                   + tracer(i_ne,  j_ne-1,:,taup1) ) /3.
    endif
 
    !--- update the tracer at global halo point to make the gradient accross boundary is 0
    !    but not yet - update domains first in tracer
    if(debug_this_module) then
       call write_chksum_3d('After ocean_obc_tracer, tracer', tracer(COMP,:,taup1))
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_tracer
  !</SUBROUTINE> NAME="ocean_obc_tracer"


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

    integer :: i, j, m, is, ie, js, je, i1, j1, n, nn, ii, jj
    real    :: u1, u2, x1, x2, v1, v2, y1, y2
    integer, dimension(nobc) :: uptype     ! OBC type
    integer, dimension(:), pointer  :: iloc, jloc, oi1, oj1
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
         'ocean_obc_mod(ocean_obc_update_boundary) : update_type= '// update_type// &
         ' should be either "n", "t", "s", "x", "z", "u",  or "i" ' )
        endif
      endif
    enddo

    if(grid_type .ne. 'T' .and. grid_type .ne. 'C'.and. grid_type .ne. 'Z'.and. grid_type .ne. 'M') call mpp_error(FATAL, &
         'ocean_obc_mod(ocean_obc_update_boundary) : grid_type= '// grid_type//' should be either "T", "C", "Z" or "M" ' )

    ! Loop through the boundaries and update
    do m = 1, nobc

      if(.not.Obc%bound(m)%on_bound) cycle
       if((Obc%bound(m)%direction == WEST .or. &
           Obc%bound(m)%direction == EAST) .and. grid_type == 'Z') cycle  
       if((Obc%bound(m)%direction == NORTH .or. &
           Obc%bound(m)%direction == SOUTH) .and. grid_type == 'M') cycle  

       ! Set pointers
       iloc => Obc%bound(m)%ilod
       jloc => Obc%bound(m)%jlod
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
                iloc => Obc%bound(m)%di1
                jloc => Obc%bound(m)%dj1
                oi1 => Obc%bound(m)%di2
                oj1 => Obc%bound(m)%dj2
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          do n = 1, Obc%bound(m)%nlod
             i = iloc(n)    ! i boundary location
             j = jloc(n)    ! j boundary location
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
                iloc => Obc%bound(m)%di1
                jloc => Obc%bound(m)%dj1
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          if(grid_type /= 'T') then
             do n = 1, Obc%bound(m)%nlod
                i = iloc(n)    ! Internal i boundary location
                j = jloc(n)    ! Internal j boundary location
                is = istr(n)   ! Start of i halo loop
                ie = iend(n)   ! End of i halo loop
                js = jstr(n)   ! Start of j halo loop
                je = jend(n)   ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)*Grd%umask(ii,jj,1)
                   enddo
                enddo
             enddo
          else
             do n = 1, Obc%bound(m)%nlod
                i = iloc(n)    ! Internal i boundary location
                j = jloc(n)    ! Internal j boundary location
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
                      field(ii,jj) = field(i,j)*Grd%umask(ii,jj,1)
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
                iloc => Obc%bound(m)%di1
                jloc => Obc%bound(m)%dj1
                istr => Obc%bound(m)%ilod
                jstr => Obc%bound(m)%jlod
             endif
          endif

          if(grid_type == 'M' .or. grid_type == 'Z') then
             do n = 1, Obc%bound(m)%nlod
                i = Obc%bound(m)%ico(n) ! i outside elevation face location
                j = Obc%bound(m)%jco(n) ! j outside elevation face location
                is = istr(n)            ! Start of i halo loop
                ie = iend(n)            ! End of i halo loop
                js = jstr(n)            ! Start of j halo loop
                je = jend(n)            ! End of j halo loop
                do jj = js, je
                   do ii = is, ie
                      field(ii,jj) = field(i,j)*Grd%umask(ii,jj,1)
                   enddo
                enddo
             enddo
          endif
       endif

! There may be open points in the shadow of corners. These t-cells can be
! open or closed. It should be open, but who knows what happens ... 

      if ( south_west_corner ) then
         if(grid_type == 'T') then
            field(i_sw-1,j_sw-1) = field(i_sw,j_sw) * Grd%tmask(i_sw,j_sw,1)
            field(max(isd,i_sw-2),j_sw-1) = field(i_sw,j_sw) * Grd%tmask(i_sw,j_sw,1)
            field(i_sw-1,max(jsd,j_sw-2)) = field(i_sw,j_sw) * Grd%tmask(i_sw,j_sw,1)
            field(max(isd,i_sw-2),max(jsd,j_sw-2)) = field(i_sw,j_sw) * Grd%tmask(i_sw,j_sw,1)
         else
            field(i_sw-1,j_sw-1) = field(i_sw,j_sw) * Grd%umask(i_sw,j_sw,1)
            field(max(isd,i_sw-2),j_sw-1) = field(i_sw,j_sw) * Grd%umask(i_sw,j_sw,1)
            field(i_sw-1,max(jsd,j_sw-2)) = field(i_sw,j_sw) * Grd%umask(i_sw,j_sw,1)
            field(max(isd,i_sw-2),max(jsd,j_sw-2)) = field(i_sw,j_sw) * Grd%umask(i_sw,j_sw,1)
         endif
      endif
      if ( south_east_corner ) then
         if(grid_type == 'T') then
            field(i_se+1,j_se-1) = field(i_se,j_se) * Grd%tmask(i_se,j_se,1)
            field(min(ied,i_se+2),j_se-1) = field(i_se,j_se) * Grd%tmask(i_se,j_se,1)
            field(i_se+1,max(jsd,j_se-2)) = field(i_se,j_se) * Grd%tmask(i_se,j_se,1)
            field(min(ied,i_se+2),max(jsd,j_se-2)) = field(i_se,j_se) * Grd%tmask(i_se,j_se,1)
         else
            field(i_se+1,j_se-1) = field(i_se,j_se) * Grd%umask(i_se,j_se,1)
            field(min(ied,i_se+2),j_se-1) = field(i_se,j_se) * Grd%umask(i_se,j_se,1)
            field(i_se+1,max(jsd,j_se-2)) = field(i_se,j_se) * Grd%umask(i_se,j_se,1)
            field(min(ied,i_se+2),max(jsd,j_se-2)) = field(i_se,j_se) * Grd%umask(i_se,j_se,1)
         endif
      endif
      if ( north_west_corner ) then
         if(grid_type == 'T') then
            field(i_nw-1,j_nw+1) = field(i_nw,j_nw) * Grd%tmask(i_nw,j_nw,1)
            field(i_nw-1,min(jed,j_nw+2)) = field(i_nw,j_nw) * Grd%tmask(i_nw,j_nw,1)
            field(max(isd,i_nw-2),j_nw+1) = field(i_nw,j_nw) * Grd%tmask(i_nw,j_nw,1)
            field(max(isd,i_nw-2),min(jed,j_nw+2)) = field(i_nw,j_nw) * Grd%tmask(i_nw,j_nw,1)
         else
            field(i_nw-1,j_nw+1) = field(i_nw,j_nw) * Grd%umask(i_nw,j_nw,1)
            field(i_nw-1,min(jed,j_nw+2)) = field(i_nw,j_nw) * Grd%umask(i_nw,j_nw,1)
            field(max(isd,i_nw-2),j_nw+1) = field(i_nw,j_nw) * Grd%umask(i_nw,j_nw,1)
            field(max(isd,i_nw-2),min(jed,j_nw+2)) = field(i_nw,j_nw) * Grd%umask(i_nw,j_nw,1)
         endif
      endif
      if ( north_east_corner ) then
         if(grid_type == 'T') then
            field(i_ne+1,j_ne+1) = field(i_ne,j_ne) * Grd%tmask(i_ne,j_ne,1)
            field(i_ne+1,min(jed,j_ne+2)) = field(i_ne,j_ne) * Grd%tmask(i_ne,j_ne,1)
            field(min(ied,i_ne+2),j_ne+1) = field(i_ne,j_ne) * Grd%tmask(i_ne,j_ne,1)
            field(min(ied,i_ne+2),min(jed,j_ne+2)) = field(i_ne,j_ne) * Grd%tmask(i_ne,j_ne,1)
         else
            field(i_ne+1,j_ne+1) = field(i_ne,j_ne) * Grd%umask(i_ne,j_ne,1)
            field(i_ne+1,min(jed,j_ne+2)) = field(i_ne,j_ne) * Grd%umask(i_ne,j_ne,1)
            field(min(ied,i_ne+2),j_ne+1) = field(i_ne,j_ne) * Grd%umask(i_ne,j_ne,1)
            field(min(ied,i_ne+2),min(jed,j_ne+2)) = field(i_ne,j_ne) * Grd%umask(i_ne,j_ne,1)
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
  !--- update field on the halo points at the boundaries 
  !<SUBROUTINE NAME="ocean_obc_update_boundary_4d" INTERFACE="obc_update_boundary">
  !   <INOUT NAME="field" TYPE="real, dimension(:,:,:,:)"></INOUT> 
  subroutine ocean_obc_update_boundary_4d(field,grid_type,update_type)
    real, dimension(isd:,jsd:,:,:), intent(inout) :: field
    character(len=1),                  intent(in) :: grid_type
    character(len=1), optional,        intent(in) :: update_type
    integer :: n3, n4, k, n

    call mpp_clock_begin(id_obc)
    
    n3 = size(field,3)
    n4 = size(field,4)
    
    if(PRESENT(update_type)) then
      do k = 1, n3
        do n = 1, n4
          call ocean_obc_update_boundary_2d(field(:,:,k,n),grid_type,update_type)
        enddo
      enddo
    else
      do k = 1, n3
        do n = 1, n4
          call ocean_obc_update_boundary_2d(field(:,:,k,n),grid_type)
        enddo
      enddo
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_update_boundary_4d
  !</SUBROUTINE> NAME="ocean_obc_update_boundary_4d"


  !#####################################################################
  !--- set field on the boundaries to zero, do not change halos
  !<SUBROUTINE NAME="ocean_obc_zero_boundary_2d" INTERFACE="ocean_obc_zero_boundary">
  subroutine ocean_obc_zero_boundary_2d(field, grid_type)

    real, dimension(isd:,jsd:), intent(inout)  :: field
    character(len=1),           intent(in) :: grid_type

    integer :: i, j, m
    
    if(grid_type .ne. 'T' .and. grid_type .ne. 'C'.and. grid_type .ne. 'Z'.and. grid_type .ne. 'M') call mpp_error(FATAL,&
         'ocean_obc_mod(ocean_obc_zero_boundary) : grid_type= '// grid_type//' should be either "T", "C", "Z" or "M" ' )
    
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
  !--- set field on the boundaries to zero, do not change halos
  !<SUBROUTINE NAME="ocean_obc_zero_boundary_3d" INTERFACE="ocean_obc_zero_boundary">
  subroutine ocean_obc_zero_boundary_3d(field, grid_type)

    real, dimension(isd:,jsd:,:), intent(inout)  :: field
    character(len=1),           intent(in) :: grid_type

    integer :: i, j, m
    integer :: k, n3

    n3 = size(field,3)
    do k = 1, n3
       call ocean_obc_zero_boundary_2d(field(:,:,k),grid_type)
    enddo
             
    return

  end subroutine ocean_obc_zero_boundary_3d
  !</SUBROUTINE> NAME="ocean_obc_zero_boundary_3d"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_check_topog">
  subroutine ocean_obc_check_topog(kmt)
    integer, dimension(isd:,jsd:), intent(inout) :: kmt
    
    integer :: m, ib, jb, i, j

    do m = 1, nobc
    
       if(.not. Obc%bound(m)%on_bound) cycle

       if(debug_this_module) write(obc_out_unit(m),*) 'check topog ', trim(Obc%bound(m)%name)
       select case( Obc%bound(m)%direction )
       case (WEST) 
          ib   = Obc%bound(m)%is
! check, if the open boundary is at a western domain boundary        
          ! if this is the global domain boundary everything is fine otherwise stop
          if (isg /= isc) then
            if ((ib-isc) < 0) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
          endif
          ! now check for the eastern domain boundary
          if ((ied-ib) <= 2) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
! check if kmt is the same for the boundary point and the two first inner points
          do j = jsd, jed
             if(Obc%bound(m)%mask(j)) then
                if(kmt(ib,j) /= kmt(ib+1,j)) call mpp_error(NOTE,'ocean_obc_mod: topography near western boundary is not smoothed') 
                if(kmt(ib,j) /= kmt(ib+2,j)) call mpp_error(NOTE,'ocean_obc_mod: topography near western boundary is not smoothed') 
             endif
          enddo
       case (EAST) 
          ib   = Obc%bound(m)%is
! check, if the open boundary is at a eastern domain boundary        
          ! if this is the global domain boundary everything is fine otherwise stop
          if (ieg /= iec) then
            if ((iec-ib) < 0) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
          endif
          ! now check for the western domain boundary
          if ((ib-isd) <= 2) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
! check if kmt is the same for the boundary point and the two first inner points
          do j = jsd, jed
             if(Obc%bound(m)%mask(j)) then
                if(kmt(ib,j) /= kmt(ib-1,j)) call mpp_error(NOTE,'ocean_obc_mod: topography near western boundary is not smoothed') 
                if(kmt(ib,j) /= kmt(ib-2,j)) call mpp_error(NOTE,'ocean_obc_mod: topography near western boundary is not smoothed') 
             endif
          enddo
       case (SOUTH) 
          jb   = Obc%bound(m)%js
! check, if the open boundary is at a southern domain boundary        
          ! if this is the global domain boundary everything is fine otherwise stop
          if (jsg /= jsc) then
            if ((jb-jsc) < 0) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
          endif
          ! now check for the northern domain boundary
          if ((jec-jb) <= 2) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
! check if kmt is the same for the boundary point and the two first inner points
          do i = isd, ied
             if(Obc%bound(m)%mask(i)) then
                if (kmt(i,jb) /= kmt(i,jb+1)) call mpp_error(NOTE, &
                     'ocean_obc_mod: topography near southern boundary is not smoothed') 
                if (kmt(i,jb) /= kmt(i,jb+2)) call mpp_error(NOTE, &
                     'ocean_obc_mod: topography near southern boundary is not smoothed') 
             endif
          enddo
       case (NORTH)
          jb   = Obc%bound(m)%js
! check, if the open boundary is at a northern domain boundary        
          ! if this is the global domain boundary everything is fine otherwise stop
          if (jeg /= jec) then
            if ((jec-jb) < 0) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
          endif
          ! now check for the southern domain boundary
          if ((jb-jsc) <= 2) call mpp_error(FATAL,'ocean_obc_mod, OBC at domain boundary')
! check if kmt is the same for the boundary point and the two first inner points
          do i = isd, ied
             if(Obc%bound(m)%mask(i)) then
                if (kmt(i,jb) /= kmt(i,jb-1)) call mpp_error(NOTE, &
                    'ocean_obc_mod: topography near northern boundary is not smoothed') 
                if (kmt(i,jb) /= kmt(i,jb-2)) call mpp_error(NOTE, &
                    'ocean_obc_mod: topography near northern boundary is not smoothed') 
             endif
          enddo
          
       end select
    enddo

  end subroutine ocean_obc_check_topog
  !</SUBROUTINE> NAME="ocean_obc_check_topog"


  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_set_mask">
  subroutine ocean_obc_set_mask
    integer :: i, j, m

    do j=jsd,jed
      do i=isd,ied
         Grd%obc_tmask(i,j) = min(1.0, float(Grd%kmt(i,j)))
         Grd%obc_umask(i,j) = min(1.0, float(Grd%kmu(i,j)))
      enddo
    enddo

!   Zero out OBC points in Grd%obc_xmask to exclude OBC points from diagnostics    

    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle
      
      select case( Obc%bound(m)%direction )
      case(WEST)
        do i = isd, Obc%bound(m)%is
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Grd%obc_tmask(i,j) = 0.
            Grd%obc_umask(i,j) = 0.
          enddo
        enddo
      case(EAST)
        do i = Obc%bound(m)%is, ied
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Grd%obc_tmask(i,j)   = 0.
            Grd%obc_umask(i-1,j) = 0.
          enddo
        enddo
      case(SOUTH)
        do j = jsd, Obc%bound(m)%js
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            Grd%obc_tmask(i,j) = 0.
            Grd%obc_umask(i,j) = 0.
          enddo
        enddo
      case(NORTH)
        do j = Obc%bound(m)%js, jed
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            Grd%obc_tmask(i,j)   = 0.
            Grd%obc_umask(i,j-1) = 0.
          enddo
        enddo
      end select
    enddo

    return
    
  end subroutine ocean_obc_set_mask
  !</SUBROUTINE> NAME="ocean_obc_set_mask"


  !#####################################################################
! <SUBROUTINE NAME="ocean_obc_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_obc_restart()

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

  subroutine ocean_obc_end(Time, have_obc)

    type(ocean_time_type), intent(in) :: Time
    logical, intent(inout)            :: have_obc
    integer :: m

    integer :: stdoutunit 
    stdoutunit=stdout() 

    call ocean_obc_end_baro(Time)
    do m = 1, nobc
      if(obc_out_unit(m) .NE. stdoutunit ) then
         call mpp_close(obc_out_unit(m))
      endif
    enddo   
    
    write(stdoutunit,*) ' '
    write(stdoutunit,*) 'Ending tracer OBC'
    call write_timestamp(Time%model_time)
    
    deallocate( Obc%bound, Obc%eta, Obc%ud, Obc%tracer )
    deallocate( wrk, wrk1, wrk2, wrk3 )

    module_is_initialized = .FALSE.
    have_obc = .FALSE.
    
    nullify(Grd)
    nullify(Dom)
    
    return

  end subroutine ocean_obc_end
  !</SUBROUTINE> NAME="ocean_obc_end"


  !#######################################################################
  !<SUBROUTINE NAME="ocean_obc_mass_flux">
  subroutine ocean_obc_mass_flux(Time, Ext_mode, mass_flux)
    type(ocean_time_type), intent(in)          :: Time
    type(ocean_external_mode_type), intent(in) :: Ext_mode
    real, dimension(isd:,jsd:), intent(inout)  :: mass_flux
    integer                                    :: i, j, m, tau
    real                                       :: uh, uhjm, vh, vhim
    logical                                    :: used

    tau = Time%tau 

    mass_flux = 0.
    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle

      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        j = Obc%bound(m)%js-1
        uhjm = Ext_mode%udrho(i,j,1,tau)*Grd%dyu(i,j)
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          uh = Ext_mode%udrho(i,j,1,tau)*Grd%dyu(i,j)
          mass_flux(i,j) = 0.5*(uh+uhjm)
          uhjm = uh
        enddo
        if(Obc%bound(m)%id_transport > 0) used = send_data(Obc%bound(m)%id_transport, &
             mass_flux(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je),  &
             Time%model_time,                               &
             rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1))
      case(EAST)
        i = Obc%bound(m)%is-1
        j = Obc%bound(m)%js-1
        uhjm = Ext_mode%udrho(i,j,1,tau)*Grd%dyu(i,j)
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          uh = Ext_mode%udrho(i,j,1,tau)*Grd%dyu(i,j)
          mass_flux(i,j) = - 0.5*(uh+uhjm)
          uhjm = uh
        enddo
        if(Obc%bound(m)%id_transport > 0) used = send_data(Obc%bound(m)%id_transport, &
             mass_flux(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je),  &
             Time%model_time,                               &
             rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1))
      case(SOUTH)
        j = Obc%bound(m)%js
        i = Obc%bound(m)%is-1
        vhim = Ext_mode%udrho(i,j,2,tau)*Grd%dxu(i,j)
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          vh = Ext_mode%udrho(i,j,2,tau)*Grd%dxu(i,j)
          mass_flux(i,j) = 0.5*(vh+vhim)
          vhim = vh
        enddo
        if(Obc%bound(m)%id_transport > 0) used = send_data(Obc%bound(m)%id_transport, &
             mass_flux(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je),  &
             Time%model_time,                               &
             rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1))
      case(NORTH)
        j = Obc%bound(m)%js-1
        i = Obc%bound(m)%is-1
        vhim = Ext_mode%udrho(i,j,2,tau)*Grd%dxu(i,j)
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          vh = Ext_mode%udrho(i,j,2,tau)*Grd%dxu(i,j)
          mass_flux(i,j) = - 0.5*(vh+vhim)
          vhim = vh
        enddo
        if(Obc%bound(m)%id_transport > 0) used = send_data(Obc%bound(m)%id_transport, &
             mass_flux(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je),  &
             Time%model_time,                               &
             rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,Obc%bound(m)%js:Obc%bound(m)%je,1))
      end select
    enddo

    return   
  end subroutine ocean_obc_mass_flux
  !</SUBROUTINE> NAME="ocean_obc_mass_flux"


  !#######################################################################
  !<SUBROUTINE NAME="ocean_obc_tracer_flux">
  subroutine ocean_obc_tracer_flux(Tracer, tracer_flux)

    type(ocean_prog_tracer_type),  intent(in) :: Tracer
    real, dimension(isd:,jsd:), intent(inout) :: tracer_flux
    integer                                   :: i, j, k, m
    

    tracer_flux = 0.
    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle
      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        do k=1,nk
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(j,k)
          enddo
        enddo
      case(EAST)
        i = Obc%bound(m)%is-1
        do k=1,nk
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(j,k)
          enddo
        enddo
      case(SOUTH)
        j = Obc%bound(m)%js
        do k=1,nk
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(i,k)
          enddo
        enddo
      case(NORTH)
        j = Obc%bound(m)%js-1
        do k=1,nk
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(i,k)
          enddo
        enddo
      end select
    enddo

    return   
  end subroutine ocean_obc_tracer_flux
  !</SUBROUTINE> NAME="ocean_obc_tracer_flux"


  !#######################################################################
  !<SUBROUTINE NAME="store_ocean_obc_tracer_flux">
  subroutine store_ocean_obc_tracer_flux(Time, Tracer, flux, n, dir, caller)
    type(ocean_time_type), intent(in)           :: Time
    type(ocean_prog_tracer_type), intent(inout) :: Tracer
    real, dimension(isd:,jsd:,:), intent(in)    :: flux
    integer, intent(in)                         :: n
    character(len=1),             intent(in)    :: dir
    character(len=3),             intent(in)    :: caller
    integer                                     :: i, j, m
    logical                                     :: used

! The tracer flux is accumulated for tracer diagnostics
! The output is written as snapshot or time mean
! dir is the direction of the flux not of the boundary
    if (dir == 'z') then
       do m = 1, nobc
          if(.not. Obc%bound(m)%on_bound) cycle
          select case( Obc%bound(m)%direction )
          case(WEST)
             i = Obc%bound(m)%is
             do j = Obc%bound(m)%js, Obc%bound(m)%je
                Tracer%otf(m)%flux(j,:) = Tracer%otf(m)%flux(j,:) + flux(i,j,:)
             enddo
             if(Obc%bound(m)%id_tracer_flux_dif(n) > 0 .and. caller == 'dif')             &
                  used = send_data(Obc%bound(m)%id_tracer_flux_dif(n), &
                  Tracer%conversion*flux(i:i,Obc%bound(m)%js:Obc%bound(m)%je,1:nk), &
                  Time%model_time, rmask=Grd%tmask(i:i,Obc%bound(m)%js:Obc%bound(m)%je,1:nk))
             if(Obc%bound(m)%id_tracer_flux_adv(n) > 0 .and. caller == 'adv')             &
                  used = send_data(Obc%bound(m)%id_tracer_flux_adv(n), &
                  Tracer%conversion*flux(i:i,Obc%bound(m)%js:Obc%bound(m)%je,1:nk), &
                  Time%model_time, rmask=Grd%tmask(i:i,Obc%bound(m)%js:Obc%bound(m)%je,1:nk))
          case(EAST)
             i = Obc%bound(m)%is-1
             do j = Obc%bound(m)%js, Obc%bound(m)%je
                Tracer%otf(m)%flux(j,:) = Tracer%otf(m)%flux(j,:) - flux(i,j,:) ! negative for eastward flux
             enddo
             if(Obc%bound(m)%id_tracer_flux(n) > 0) used = send_data(Obc%bound(m)%id_tracer_flux(n), &
                  Tracer%conversion*flux(i:i,Obc%bound(m)%js:Obc%bound(m)%je,1:nk), &
                  Time%model_time, rmask=Grd%tmask(i:i,Obc%bound(m)%is:Obc%bound(m)%ie,1:nk))
          end select
       enddo
    elseif (dir == 'm') then
       do m = 1, nobc
          if(.not. Obc%bound(m)%on_bound) cycle
          select case( Obc%bound(m)%direction )
          case(SOUTH)
             j = Obc%bound(m)%js
             do i = Obc%bound(m)%is, Obc%bound(m)%ie
                Tracer%otf(m)%flux(i,:) = Tracer%otf(m)%flux(i,:) + flux(i,j,:)
             enddo
             if(Obc%bound(m)%id_tracer_flux(n) > 0) used = send_data(Obc%bound(m)%id_tracer_flux(n), &
                  Tracer%conversion*flux(Obc%bound(m)%is:Obc%bound(m)%ie,j:j,1:nk), &
                  Time%model_time, rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,j:j,1:nk))
          case(NORTH)
             j = Obc%bound(m)%js-1
             do i = Obc%bound(m)%is, Obc%bound(m)%ie
                Tracer%otf(m)%flux(i,:) = Tracer%otf(m)%flux(i,:) - flux(i,j,:) ! negative for northward flux
             enddo
             if(Obc%bound(m)%id_tracer_flux(n) > 0) used = send_data(Obc%bound(m)%id_tracer_flux(n), &
                  Tracer%conversion*flux(Obc%bound(m)%is:Obc%bound(m)%ie,j:j,1:nk), &
                  Time%model_time, rmask=Grd%tmask(Obc%bound(m)%is:Obc%bound(m)%ie,j:j,1:nk))
          end select
       enddo
    endif
        
    return   
  end subroutine store_ocean_obc_tracer_flux
  !</SUBROUTINE> NAME="store_ocean_obc_tracer_flux"


  !#######################################################################
  !<SUBROUTINE NAME="store_ocean_obc_pressure_grad">
  ! Store the pressure gradient across the boundary.
  subroutine store_ocean_obc_pressure_grad(Thickness, pressure_gradient)

   type(ocean_thickness_type), intent(in)            :: Thickness
   real, dimension(isc:iec,jsc:jec,nk,2), intent(in) :: pressure_gradient
   integer                                           :: i, j, k, m

    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle

      Obc%bound(m)%pgrad(:) = 0

      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(j) = Obc%bound(m)%pgrad(j) + &
                                      pressure_gradient(i,j,k,1) * &
                                      Thickness%dzu(i,j,k) * Grd%umask(i,j,k)
          enddo
        enddo
      case(EAST)
        i = Obc%bound(m)%is-1
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(j) = Obc%bound(m)%pgrad(j) + &
                                      pressure_gradient(i,j,k,1) * &
                                      Thickness%dzu(i,j,k) * Grd%umask(i,j,k)
          enddo
        enddo
      case(SOUTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(i) = Obc%bound(m)%pgrad(i) + &
                                      pressure_gradient(i,j,k,2) * &
                                      Thickness%dzu(i,j,k) * Grd%umask(i,j,k)
          enddo
        enddo
      case(NORTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(i) = Obc%bound(m)%pgrad(i) + &
                                      pressure_gradient(i,j,k,2) * &
                                      Thickness%dzu(i,j,k) * Grd%umask(i,j,k)
          enddo
        enddo
      end select
    enddo
    
    return   
  end subroutine store_ocean_obc_pressure_grad
  !</SUBROUTINE> NAME="store_ocean_obc_pressure_grad"


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

end module ocean_obc_mod
