! The include file fms_platform.h will handle the conversion of POINTER to ALLOCATABLE arrays
! for derived type members. The conversion affects performance only and should not change 
! any numeric result. It is limited to member arrays that are used within MOM only
! and to arrays that are never associated (=>) with another array.

! Fortran 90 requires array members of derived type to have the POINTER attribute. 
! However, most Fortran 95 compilers now also support ALLOCATABLE array components 
! (a Fortran 2003 feature). This avoids the aliasing problem afflicting pointers. 
! Some compilers may require an additional switch (e.g. -fno-alias) to fully exploit
! the performance benefit of the conversion.
!
! Macros used from fms_platform.h:
! _ALLOCATABLE maps to either POINTER  or ALLOCATABLE
! _NULL        maps to either =>NULL() or "nothing"
#include <fms_platform.h>

module ocean_types_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module contains type declarations and default values for ocean model.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains type declarations and default values for ocean model.
! Multiple model realizations need to be distinguished by
! an ensemble_id for use by the diag_manager.
!</DESCRIPTION>
!
  use coupler_types_mod,only: coupler_2d_bc_type
  use field_manager_mod,only: fm_field_name_len, fm_string_len
  use fms_mod,          only: write_version_number
  use mpp_domains_mod,  only: domain2d
  use mpp_mod,          only: FATAL, mpp_error
  use time_manager_mod, only: time_type

  implicit none

  private

  logical :: module_is_initialized=.false.
  character(len=128) :: version = &
     '$Id: ocean_types.F90,v 20.0 2013/12/14 00:12:37 fms Exp $'
  character (len=128) :: tagname = &
     '$Name: tikal $'


  type, public :: obc_flux
     real, _ALLOCATABLE, dimension(:,:) :: flux  _NULL   ! flux through a boundary 
  end type obc_flux

  type, public :: ocean_time_steps_type
     character(len=32)  :: time_tendency      ! either "twolevel" or "threelevel" 
     character(len=32)  :: barotropic_scheme  ! for diagnostic printout 
     integer :: tendency ! either "TWO_LEVEL" or "THREE_LEVEL"
     real    :: aidif    ! aidif=1.0 for fully implicit vertical mixing and 0.0 for fully explicit
     real    :: acor     ! acor=0.0 for explicit Coriolis, 0.5 <= acor <= 1.0 for semi-implicit
     real    :: dtts     ! tracer timestep (seconds)
     real    :: dtuv     ! internal mode timestep (seconds)
     real    :: dtbt     ! barotropic timestep (seconds)
     real    :: dteta    ! ocean volume time step (eta_t) (seconds)
     real    :: dtime_t  ! 2*dtts  if threelevel, 1*dtts  if twolevel
     real    :: dtime_u  ! 2*dtuv  if threelevel, 1*dtuv  if twolevel
     real    :: dtime_e  ! 2*dteta if threelevel, 1*dteta if twolevel
     real    :: robert_asselin_param ! parameter used for time filtering with leap-frog scheme 
  end type ocean_time_steps_type

  type, public :: ocean_options_type
     integer            :: vertmix
     character(len=128) :: vert_mix
     character(len=72)  :: tidal_wave_mix
     character(len=72)  :: tidal_drag_mix
     character(len=128) :: leewave_mix
     character(len=72)  :: bryan_lewis_mix
     character(len=72)  :: hwf_mix
     character(len=72)  :: tanh_diff_cbt
     character(len=72)  :: read_diff_cbt
     character(len=72)  :: j09_diff_cbt
     character(len=72)  :: horz_bih_tracer
     character(len=72)  :: horz_lap_tracer
     character(len=72)  :: horz_lap_friction
     character(len=72)  :: horz_bih_friction
     character(len=128) :: coriolis
     character(len=72)  :: momentum_source
     character(len=72)  :: form_drag
     character(len=128) :: bottom_roughness
     character(len=72)  :: bmf_implicit
     character(len=72)  :: geothermal_heating
     character(len=72)  :: OBC
     character(len=72)  :: override_velocity
     character(len=72)  :: baroclinic_tendency
     character(len=72)  :: barotropic_tendency
     character(len=72)  :: tracer_tendency
     character(len=72)  :: equation_of_state
     character(len=72)  :: temperature_variable
     character(len=128) :: frazil_ice
     character(len=72)  :: convective_adjustment
     character(len=72)  :: neutral_physics
     character(len=72)  :: neutral_physics_new
     character(len=72)  :: submesoscale
     character(len=72)  :: sigma_transport
     character(len=72)  :: lagrangian_blobs
     character(len=128) :: overexchange
     character(len=72)  :: mixdownslope
     character(len=72)  :: overflow
     character(len=72)  :: overflow_OFP
     character(len=72)  :: shortwave
     character(len=72)  :: xlandmix
     character(len=128) :: xlandinsert
     character(len=72)  :: rivermix 
     character(len=72)  :: riverspread
     character(len=72)  :: passive_tracers
     character(len=72)  :: ocean_sponges_eta
     character(len=72)  :: ocean_sponges_tracer
     character(len=72)  :: ocean_sponges_velocity
     character(len=72)  :: ocean_ideal_surf_wave
     character(len=72)  :: fafmip_heat
  end type ocean_options_type



#include <ocean_memory.h>


#ifdef MOM_STATIC_ARRAYS
!################################################################################################

  type, public :: ocean_thickness_type

     integer :: method       ! energetic or finite volume  

     real, dimension(isd:ied,jsd:jed,nk,3) :: rho_dzt   ! rho(kg/m^3)*thickness (m) of T cell at 3 times
     real, dimension(isd:ied,jsd:jed,nk,2) :: rho_dzten ! rho(kg/m^3)*thickness (m) at east/north face of T-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: rho_dztr  ! 1.0/(rho*dzt) at time taup1
     real, dimension(isd:ied,jsd:jed,nk,3) :: rho_dzu   ! rho (kg/m^3) * thickness (m) of U cell at 3 times
     real, dimension(isd:ied,jsd:jed,nk)   :: rho_dzur  ! 1.0/(rho*dzu) at time taup1

     real, dimension(isd:ied,jsd:jed,nk)   :: rho_dzt_tendency ! rho_dzt tendency (kg/m^3)*(m/s)

     real, dimension(isd:ied,jsd:jed)      :: sea_lev     ! eta_t + patm/(rho0*grav) - eta_geoid - eta_tide (m) at time taup1 for coupler 
     real, dimension(isd:ied,jsd:jed,nk)   :: dzt         ! thickness (m) of T cell at time tau/taup1
     real, dimension(isd:ied,jsd:jed,nk,2) :: dzten       ! thickness (m) of east/north face of T cell at time tau/taup1
     real, dimension(isd:ied,jsd:jed,nk)   :: dzu         ! thickness (m) of U cell at time tau/taup1
     real, dimension(isd:ied,jsd:jed,0:nk) :: dzwt        ! vertical distance (m) between T points at tau/taup1 
     real, dimension(isd:ied,jsd:jed,0:nk) :: dzwu        ! vertical distance (m) between U points at tau/taup1

     real, dimension(isd:ied,jsd:jed,nk)   :: dztup       ! distance (m) from T-cell point to top of T-cell 
     real, dimension(isd:ied,jsd:jed,nk)   :: dztlo       ! distance (m) from T-cell point to bottom of T-cell 

     real, dimension(isd:ied,jsd:jed,nk)   :: geodepth_zt ! vertical distance (m) from z=0 to T-point 
     real, dimension(isd:ied,jsd:jed,nk)   :: depth_zt    ! vertical distance (m) from column top to T-point 
     real, dimension(isd:ied,jsd:jed,nk)   :: depth_zwt   ! vertical distance (m) from column top to T-bottom 
     real, dimension(isd:ied,jsd:jed,nk)   :: depth_zu    ! vertical distance (m) from column top to U-point
     real, dimension(isd:ied,jsd:jed,nk)   :: depth_zwu   ! vertical distance (m) from column top to U-bottom

     real, dimension(isd:ied,jsd:jed,nk)   :: depth_st    ! s-distance to T-point
     real, dimension(isd:ied,jsd:jed,nk)   :: depth_swt   ! s-distance to T-bottom
     real, dimension(isd:ied,jsd:jed,nk)   :: dst         ! vertical increment of s-coordinate on T-cell 
     real, dimension(isd:ied,jsd:jed,0:nk) :: dswt        ! s-coordinate T-point increment at time tau 
     real, dimension(isd:ied,jsd:jed,nk)   :: dzt_dst     ! specific thickness (metre/s-coordinate) on T-cell 

     real, dimension(isd:ied,jsd:jed,nk)   :: dstup       ! s-distance (s-units) from T-cell point to top of T-cell 
     real, dimension(isd:ied,jsd:jed,nk)   :: dstlo       ! s-distance (s-units) from T-cell point to bottom of T-cell 

     real, dimension(isd:ied,jsd:jed)      :: pbot0       ! reference bottom pressure (Pa) 
     real, dimension(isd:ied,jsd:jed)      :: pbot0r      ! inverse of reference bottom pressure (1/Pa) 

     real, dimension(isd:ied,jsd:jed,nk)   :: mass_source ! mass source (kg/m^3)*(m/sec) 
     real, dimension(isd:ied,jsd:jed,3)    :: mass_u      ! mass per area (kg/m^2) in a velocity column 
     real, dimension(isd:ied,jsd:jed,2)    :: mass_en     ! vertical sum rho_dzten (kg/m^2) 
     real, dimension(isd:ied,jsd:jed,3)    :: thicku      ! thickness on U-cell (metre) from z=eta_u to z=-H.
     real, dimension(isd:ied,jsd:jed,2)    :: thicken     ! sum of dzten (metre)


     ! The following arrays are never allocated with MOM_STATIC_ARRAYS,
     ! but they need to be declared allocatable in order to compile. 
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztL     _NULL ! L system contribution to rho_dztT (3 time levels)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztT     _NULL ! rho(kg/m^3)*thickness (m) of T cell at 3 times (total)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuL     _NULL ! L system contribution to rho_dzuT (3 time levels)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuT     _NULL ! rho (kg/m^3) * thickness (m) of U cell at 3 times (total)
     real, dimension(:,:,:),   _ALLOCATABLE :: dztL         _NULL ! L system contribution to dztT
     real, dimension(:,:,:,:), _ALLOCATABLE :: dztT         _NULL ! thickness (m) of T cell at time tau/taup1
     real, dimension(:,:,:),   _ALLOCATABLE :: dzuL         _NULL ! L system contribution to dzuT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzuT         _NULL ! thickness (m) of U cell at time tau/taup1
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwtL        _NULL ! L system contribution to dzwtT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwtT        _NULL ! vertical distance (m) between T points at tau/taup1
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwuL        _NULL ! L system contribution to dzwuT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwuT        _NULL ! vertical distance (m) between U points at tau/taup1
     real, dimension(:,:,:),   _ALLOCATABLE :: dztupL       _NULL ! L system contribution to dztupT
     real, dimension(:,:,:),   _ALLOCATABLE :: dztupT       _NULL ! distance (m) from T-cell point to top of T-cell 
     real, dimension(:,:,:),   _ALLOCATABLE :: dztloL       _NULL ! L system contribution to dztloT
     real, dimension(:,:,:),   _ALLOCATABLE :: dztloT       _NULL ! distance (m) from T-cell point to bottom of T-cell
     real, dimension(:,:,:),   _ALLOCATABLE :: mass_uT      _NULL ! mass per area (kg/m^2) in a velocity column
     real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zwt _NULL ! vert distance (m) from z=0 to bottom of T-cell
     real, dimension(:,:),     _ALLOCATABLE :: blob_source _NULL  ! mass source (sink) for the L system (E system)

  end type ocean_thickness_type


  type, public :: ocean_grid_type
     character(len=32) :: name

     ! geometry and topology and rotation
     logical                           :: cyclic_x         ! true if domain is cyclic in the i direction
     logical                           :: cyclic_y         ! true if domain is cyclic in the j direction
     logical                           :: tripolar         ! folded connectivity at "top" row in bipolar Arctic 
     logical                           :: mosaic           ! true when using a mosaic grid 
     logical                           :: beta_plane       ! beta plane Cartesian 
     logical                           :: f_plane          ! f-plane Cartesian
     real                              :: f_plane_latitude ! latitude where f_plane is centered  
     real, dimension(isd:ied,jsd:jed)  :: f                ! coriolis parameter at u-cell points (sec^-1)
     real, dimension(isd:ied,jsd:jed)  :: fstar            ! horizontal coriolis parameter at u-cell points (sec^-1)
     real, dimension(isd:ied,jsd:jed)  :: beta             ! df/dy at u-cell points (1/(sec*m))     
     real, dimension(isd:ied,jsd:jed)  :: beta_eff         ! beta plus topographic beta at u-cell (1/(sec*m))

     ! vertical grid information (time independent) 
     integer                             :: nk     ! number of vertical grid points 
     integer, dimension(isd:ied,jsd:jed) :: kmt    ! number of t-levels
     integer, dimension(isd:ied,jsd:jed) :: kmu    ! number of u-levels

     real, dimension(nk)                 :: zt     ! full cell depth (m) from surface to level k T-cell 
     real, dimension(nk)                 :: zw     ! full cell depth (m) from surface to bottom of k T-cell  
     real, dimension(nk)                 :: dzt    ! initial vertical resolution of T or U grid cells (m) 
     real, dimension(nk)                 :: dztlo  ! distance (m) from T-cell point to bottom of T-cell
     real, dimension(nk)                 :: dztup  ! distance (m) from T-cell point to top of T-cell
     real, dimension(0:nk)               :: dzw    ! initial vertical resolution of W grid cells (m) 
     real, dimension(0:nk)               :: dzwr   ! reciprocal of dzw (W cell vertical resolution)

     real, dimension(nk)                 :: st     ! full cell s-depth from surface to level k T-cell 
     real, dimension(nk)                 :: sw     ! full cell s-depth from surface to bottom of k T-cell  
     real, dimension(nk)                 :: dst    ! initial vertical s-resolution of T or U grid cells (m) 
     real, dimension(nk)                 :: dstlo  ! s-distance (s-units) from T-cell point to bottom of T-cell
     real, dimension(nk)                 :: dstup  ! s-distance (s-units) from T-cell point to top of T-cell
     real, dimension(0:nk)               :: dsw    ! initial vertical s-resolution of W grid cells (m) 

     real, dimension(nk,0:1)             :: fracdz ! fractional distance between grid point and cell top/bot 
     real, dimension(isd:ied,jsd:jed)    :: ht     ! depth to bottom of ocean (m) on t-cells from z=0 
     real, dimension(isd:ied,jsd:jed)    :: htr    ! inverse depth to bottom of ocean (m^-1) on t-cells
     real, dimension(isd:ied,jsd:jed)    :: hu     ! depth to bottom of ocean (m) on u-cells from z=0
     real, dimension(isd:ied,jsd:jed)    :: dht_dx ! d(ht)/dx on u-cells (m/m)
     real, dimension(isd:ied,jsd:jed)    :: dht_dy ! d(ht)/dy on u-cells (m/m)
     real, dimension(isd:ied,jsd:jed)    :: gradH  ! sqrt(dht_dx**2+dht_dyx**2) on u-cells (m/m)

     ! horizontal grid information (time independent)
     integer                          :: ni, nj     ! number of global points in the two horizontal directions
     real, dimension(isd:ied,jsd:jed) :: xt         ! longitude of the T grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: xu         ! longitude of the U grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: yt         ! latitude of the T grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: yu         ! latitude of the U grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: phiu       ! latitude of U grid point in radians
     real, dimension(isd:ied,jsd:jed) :: phit       ! latitude of T grid point in radians
     real, dimension(isd:ied,jsd:jed) :: h1t        ! metric factors in i grid direction
     real, dimension(isd:ied,jsd:jed) :: h1u        ! metric factors in j grid jirection
     real, dimension(isd:ied,jsd:jed) :: h2t        ! metric factors in i grid direction
     real, dimension(isd:ied,jsd:jed) :: h2u        ! metric factors in j grid jirection
     real, dimension(isd:ied,jsd:jed) :: dh2dx      ! (1/delta_y)*d(delta_y)/dx (1/m)
     real, dimension(isd:ied,jsd:jed) :: dh1dy      ! (1/delta_x)*d(delta_x)/dy (1/m)
     real, dimension(isd:ied,jsd:jed) :: dxt        ! longitudinal width of T-cells at grid point (m)
     real, dimension(isd:ied,jsd:jed) :: dxu        ! longitudinal width of U-cells at grid point (m)
     real, dimension(isd:ied,jsd:jed) :: dyt        ! latitudinal width of T-cells at grid point (m) 
     real, dimension(isd:ied,jsd:jed) :: dyu        ! latitudinal width of U-cells at grid point (m)
     real, dimension(isd:ied,jsd:jed) :: dau        ! area of U-cells (m^2)
     real, dimension(isd:ied,jsd:jed) :: dat        ! area of T-cells (m^2)
     real, dimension(isd:ied,jsd:jed) :: dat_frac   ! fraction of total wet area occuped by a T-cell (dimensionless)
     real, dimension(isd:ied,jsd:jed) :: dxte       ! long-width between grid points at i+1 and i in T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxtn       ! long-width of north face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dyte       ! lat-width of east face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dytn       ! lat-width between grid points at j+1 and j in T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxue       ! long-width between grid points at i+1 and i in U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxun       ! long-width of north face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dyue       ! lat-width of east face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dyun       ! lat-width between grid points at j+1 and j in U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: datnr      ! reciprocal area at north face of T-cell
     real, dimension(isd:ied,jsd:jed) :: dater      ! reciprocal area at east face of T-cell
     real, dimension(isd:ied,jsd:jed) :: dun        ! width from grid point to north face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dus        ! width from grid point to south face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: duw        ! width from grid point to west  face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: due        ! width from grid point to east  face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dtn        ! width from grid point to north face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dts        ! width from grid point to south face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dtw        ! width from grid point to west  face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dte        ! width from grid point to east  face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxtr       ! 1/dxt
     real, dimension(isd:ied,jsd:jed) :: dxur       ! 1/dxu
     real, dimension(isd:ied,jsd:jed) :: dytr       ! 1/dyt
     real, dimension(isd:ied,jsd:jed) :: dyur       ! 1/dyu
     real, dimension(isd:ied,jsd:jed) :: daur       ! 1/[area of U-cells (m^2)]
     real, dimension(isd:ied,jsd:jed) :: datr       ! 1/[area of T-cells (m^2)]
     real, dimension(isd:ied,jsd:jed) :: dxter      ! 1/dxte
     real, dimension(isd:ied,jsd:jed) :: dyter      ! 1/dyte
     real, dimension(isd:ied,jsd:jed) :: dxtnr      ! 1/dxtn
     real, dimension(isd:ied,jsd:jed) :: dytnr      ! 1/dytn
     real, dimension(isd:ied,jsd:jed) :: dxuer      ! 1/dxue
     real, dimension(isd:ied,jsd:jed) :: dxunr      ! 1/dxun
     real, dimension(isd:ied,jsd:jed) :: dyuer      ! 1/dyue
     real, dimension(isd:ied,jsd:jed) :: dyunr      ! 1/dyun
     real, dimension(isd:ied,jsd:jed) :: dyue_dxuer ! dyue/dxue 
     real, dimension(isd:ied,jsd:jed) :: dxun_dyunr ! dxun/dyun
     real, dimension(isd:ied,jsd:jed) :: dxt_dxter  ! dxt/dxte
     real, dimension(isd:ied,jsd:jed) :: dxte_dxtr  ! dxte/dxt
     real, dimension(isd:ied,jsd:jed) :: dxt_dxtnr  ! dxt/dxtn
     real, dimension(isd:ied,jsd:jed) :: dxtn_dxtr  ! dxtn/dxt
     real, dimension(isd:ied,jsd:jed) :: dyt_dyter  ! dyt/dyte
     real, dimension(isd:ied,jsd:jed) :: dyte_dytr  ! dyte/dyt
     real, dimension(isd:ied,jsd:jed) :: dyt_dytnr  ! dyt/dytn
     real, dimension(isd:ied,jsd:jed) :: dytn_dytr  ! dytn/dyt

     ! land/sea masks 
     real, dimension(isd:ied,jsd:jed)      :: obc_tmask   ! land/sea mask for T cell diagnostics with OBC 
     real, dimension(isd:ied,jsd:jed)      :: obc_umask   ! land/sea mask for B-grid U cell diagnostics with OBC
     real, dimension(isd:ied,jsd:jed,nk)   :: mask        ! land/sea tmask or umask depending on C/B grids 
     real, dimension(isd:ied,jsd:jed,nk)   :: tmask       ! land/sea mask for T cells based on s-coordinate
     real, dimension(isd:ied,jsd:jed,nk)   :: umask       ! land/sea mask for B-grid U cells based on s-coordinate
     real, dimension(isd:ied,jsd:jed,nk,2) :: tmasken     ! land/sea mask for east/north face of t-cell based on s-coordinate
     real, dimension(isd:ied,jsd:jed,nk)   :: tmask_depth ! based on depth-based vert_coordinate
     real, dimension(isd:ied,jsd:jed,nk)   :: umask_depth ! based on depth-based vert_coordinate

     ! grid areas and volumes 
     real                 :: tcellv    ! initial T cell volume m^3 (entire resting ocean)
     real                 :: ucellv    ! initial U cell volume m^3 (entire resting ocean)
     real                 :: tcellsurf ! T cell surface area (k=1) in bitwise reproducible form
     real                 :: ucellsurf ! U cell surface area (k=1) in bitwise reproducible form
     real, dimension(nk)  :: tcella    ! T cell surface area m^2 (entire ocean)
     real, dimension(nk)  :: ucella    ! U cell surface area m^2 (entire ocean)

     ! model grid information
     integer :: wet_t_points   ! total number of wet tracer points 
     integer :: wet_u_points   ! total number of wet B-grid velocity points 
     integer :: total_t_points ! total number of wet or dry tracer points  

     ! sine and cosine of rotation angles (clockwise) of velocity for tripolar 
     real, dimension(isd:ied,jsd:jed) :: sin_rot  
     real, dimension(isd:ied,jsd:jed) :: cos_rot  

     ! 1-d grid coordinates for COARDS NetCDF files
     real, dimension(ni) :: grid_x_t  
     real, dimension(nj) :: grid_y_t 
     real, dimension(ni) :: grid_x_u  
     real, dimension(nj) :: grid_y_u  

     ! axes id for diagnostic manager 
     integer, dimension(3)  :: tracer_axes        
     integer, dimension(3)  :: vel_axes_uv         
     integer, dimension(3)  :: vel_axes_u
     integer, dimension(3)  :: vel_axes_v
     integer, dimension(3)  :: vel_axes_wu    
     integer, dimension(3)  :: vel_axes_wt    
     integer, dimension(3)  :: tracer_axes_wt 
     integer, dimension(3)  :: tracer_axes_flux_x  
     integer, dimension(3)  :: tracer_axes_flux_y  
     integer, dimension(3)  :: vel_axes_flux_x    
     integer, dimension(3)  :: vel_axes_flux_y    

     ! axes id for diagnostic manager, appropriate
     ! when with to remap native vertical fields 
     ! to depth or pressure levels.  Appropropriate
     ! when vert_coordinate == ZSIGMA or PSIGMA, or other 
     ! vertical coordinats whose iso-surfaces are 
     ! not quasi-horizontal 
     integer, dimension(3)  :: tracer_axes_depth        
     integer, dimension(3)  :: vel_axes_uv_depth
     integer, dimension(3)  :: vel_axes_u_depth
     integer, dimension(3)  :: vel_axes_v_depth
     integer, dimension(3)  :: vel_axes_wu_depth    
     integer, dimension(3)  :: vel_axes_wt_depth
     integer, dimension(3)  :: tracer_axes_wt_depth
     integer, dimension(3)  :: tracer_axes_flux_x_depth
     integer, dimension(3)  :: tracer_axes_flux_y_depth
     integer, dimension(3)  :: vel_axes_flux_x_depth
     integer, dimension(3)  :: vel_axes_flux_y_depth

  end type ocean_grid_type


  type, public :: ocean_domain_type
     type(domain2d) :: domain2d           ! fms variable, used by mpp routines
     integer :: isc, iec, jsc, jec        ! computational domain indices 
     integer :: isd, ied, jsd, jed        ! local indices including halo, consistent with domain2d
     integer :: isg, ieg, jsg, jeg        ! global indices
     integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)
     integer :: xhalo, yhalo              ! halo sizes 
     integer :: xflags, yflags, layout(2) ! options to mpp_define_domains
     integer :: ioff , joff               ! index offset for absolute indices if MOM_STATIC_ARRAYS (0 otherwise)
     integer :: io_layout(2)              ! options to define io_domain.
     integer :: x_cyclic_offset           ! offset applied to x-direction cyclic boundary condition
     integer :: y_cyclic_offset           ! offset applied to y-direction cyclic boundary condition
     logical, pointer :: maskmap(:,:) =>NULL() ! option to mpp_define_domains
  end type ocean_domain_type

  type, public :: ocean_time_type
     type(time_type) :: model_time     ! fms variable
     type(time_type) :: Time_step      ! time step for tracers (and for ocean model)
     type(time_type) :: Time_init      ! time of initial conditions 
     integer :: calendar               ! calendar type defined by time_manager_mod
     logical :: init                   ! true at beginning of run (initial condition time)
     integer :: itt                    ! timestep counter measured relative to time at restart 
     integer :: itt0                   ! timestep counter measured relative to initial condition time
     integer :: taum1, tau, taup1      ! time level indices
     integer :: tau_m2, tau_m1, tau_m0 ! time level indices for Adams-Bashforth velocity advection and coriolis 
  end type ocean_time_type

  type, public :: ocean_adv_vel_type
     real, dimension(isd:ied,jsd:jed,nk)   :: uhrho_et   ! rho_dzu * advect vel (kg/(m*s)) on i-face of T-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: vhrho_nt   ! rho*dzu * advect vel (kg/(m*s)) on j-face of T-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: uhrho_eu   ! remapped uhrho_et (kg/(m*s)) onto i-face of U-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: vhrho_nu   ! remapped vhrho_nt (kg/(m*s)) onto j-face of U-cell
     real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bt    ! rho * vertical advect vel (kg/(m^2*s)) on T-bottom
     real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bu    ! remapped wrho_bt onto U-bottom
     real, dimension(isd:ied,jsd:jed,nk)   :: diverge_t  ! divergence on T-cell of horiz momentum per mass (kg/m3)*(m/s)
     real, dimension(isd:ied,jsd:jed,nk)   :: diverge_u  ! divergence on U-cell of horiz momentum per mass (kg/m3)*(m/s)
  end type ocean_adv_vel_type

  type , public :: ocean_density_type
     logical                               :: use_teos10           ! for using the TEOS-2010 equation of state recommendations 
     real, dimension(isd:ied,jsd:jed,nk,3) :: rho                  ! in situ density (kg/m^3) at time levels
     real, dimension(isd:ied,jsd:jed,nk,3) :: rho_salinity         ! salinity used in rho calculations (psu or g/kg) at time levels
     real, dimension(isd:ied,jsd:jed,nk)   :: rho_dztr_tau         ! rho_dztr at time tau 1/(kg/m^3) (for diagnostic uses) 
     real, dimension(isd:ied,jsd:jed,nk)   :: rho_fresh            ! in situ fresh water density (kg/m^3) at tau   
     real, dimension(isd:ied,jsd:jed,nk)   :: pressure_at_depth    ! hydrostatic pressure (including patm)
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodT               ! partial rho wrt theta (kg/(m^3 C)
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodS               ! partial rho wrt salinity (kg/(m^3 psu)
     real, dimension(isd:ied,jsd:jed,nk)   :: dpotrhodT            ! partial potrho wrt theta (kg/(m^3 C)
     real, dimension(isd:ied,jsd:jed,nk)   :: dpotrhodS            ! partial potrho wrt salinity (kg/(m^3 psu)
     real, dimension(isd:ied,jsd:jed,nk)   :: dpotrhodP            ! partial potrho wrt pressure (kg/(m^3 psu)
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodP               ! partial rho wrt pressure (kg/(m^3 Pa)
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodz_wt            ! d(neutral density)/dz (kg/m^4) at W-point
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodz_zt            ! d(neutral density)/dz (kg/m^4) at T-point
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodx_zt            ! d(neutral density)/dx (kg/m^4) at T-point
     real, dimension(isd:ied,jsd:jed,nk)   :: drhody_zt            ! d(neutral density)/dy (kg/m^4) at T-point
     real, dimension(isd:ied,jsd:jed,nk)   :: drhodz_diag          ! regularized drhodz_zt for diagnostics 
     real, dimension(isd:ied,jsd:jed,nk)   :: dTdz_zt              ! partial theta wrt z  (C/m) at T-point
     real, dimension(isd:ied,jsd:jed,nk)   :: dSdz_zt              ! partial salinity wrt z  (psu/m) at T-point
     real, dimension(isd:ied,jsd:jed,nk)   :: potrho               ! potential density (kg/m^3)
     real, dimension(isd:ied,jsd:jed,nk)   :: neutralrho           ! neutral density (kg/m^3)
     real, dimension(isd:ied,jsd:jed,nk)   :: watermass_factor        ! ratio (|grad nrho|/|grad local ref potrho|) / delta(gamma)
     real, dimension(isd:ied,jsd:jed,nk)   :: stratification_factor   ! rho*Area(h) / gamma_{,h}, w/ h=direction w/ max gamma strat
     real, dimension(isd:ied,jsd:jed)      :: mld_subduction          ! depth mixed layer base (m) for subduction diagnostics 
     integer, dimension(3)                 :: potrho_axes             ! axis ids for diagnosing potential density 
     integer, dimension(3)                 :: potrho_axes_flux_x      ! axis ids for diagnosing x-flux
     integer, dimension(3)                 :: potrho_axes_flux_y      ! axis ids for diagnosing y-flux
     integer, dimension(3)                 :: neutralrho_axes         ! axis ids for diagnosing neutral density 
     integer, dimension(3)                 :: neutralrho_axes_flux_x  ! axis ids for diagnosing x-flux 
     integer, dimension(3)                 :: neutralrho_axes_flux_y  ! axis ids for diagnosing y-flux 
     integer, dimension(3)                 :: theta_axes              ! axis ids for potential temperature 
     integer, dimension(3)                 :: theta_axes_flux_x       ! axis ids for diagnosing x-flux 
     integer, dimension(3)                 :: theta_axes_flux_y       ! axis ids for diagnosing y-flux 

     real, _ALLOCATABLE, dimension(:)      :: theta_ref     _NULL  ! partition vertical into theta classes 
     real, _ALLOCATABLE, dimension(:)      :: theta_bounds  _NULL  ! bounds for theta classes 
     real, _ALLOCATABLE, dimension(:)      :: potrho_ref _NULL     ! partition vertical into potrho classes 
     real, _ALLOCATABLE, dimension(:)      :: potrho_bounds _NULL  ! bounds for potrho classes 
     real, _ALLOCATABLE, dimension(:)      :: neutralrho_ref _NULL ! partition vertical into neutral density classes 
     real, _ALLOCATABLE, dimension(:)      :: neutralrho_bounds _NULL ! bounds for neutral density classes 

     ! The following array is never allocated with MOM_STATIC_ARRAYS,
     ! but it must be declared allocatable in order to compile. 
     real, _ALLOCATABLE, dimension(:,:,:)  :: rhoT              _NULL ! combined L and E in situ density (kg/m^3)

  end type ocean_density_type

  
  type, public :: ocean_prog_tracer_type
     character(len=32)  :: name
     character(len=128) :: units
     character(len=128) :: type
     character(len=128) :: longname

     logical :: use_only_advection     ! for testing purposes, evolve using ONLY advection
     logical :: neutral_physics_limit  ! neutral physics reduce to horz diffusion if tracer out of bounds
     logical :: complete               ! to determine if ready to do mpp updates

     integer :: sfc_flux_id=-1         ! index for time_interp_external
     integer :: horz_advect_scheme=-1  ! id for horizontal advection scheme
     integer :: vert_advect_scheme=-1  ! id for vertical advection scheme
     integer :: ppm_hlimiter=1         ! Limiter for use with PPM in horizontal
     integer :: ppm_vlimiter=1         ! Limiter for use with PPM in vertical
     integer :: mdt_scheme=4           ! Version of Multi-Dim. Modified Daru & Tenaud (MDMDT)

     type(obc_flux), _ALLOCATABLE, dimension(:) :: otf   _NULL ! flux through open boundaries, allocate nobc

     real, dimension(isd:ied,jsd:jed,nk,3) :: field            ! tracer concentration at 3 time levels
     real, dimension(isd:ied,jsd:jed,nk)   :: th_tendency      ! thickness weighted tracer tendency
     real, dimension(isd:ied,jsd:jed,nk)   :: tendency         ! for diagnostics: tendency concentration [concentration/sec]
     real, dimension(isd:ied,jsd:jed,nk)   :: source           ! tracer source [=tracer concentration per time]
     real, dimension(isd:ied,jsd:jed,nk)   :: wrk1             ! work array
     real, dimension(isd:ied,jsd:jed,nk)   :: tmask_limit      ! to limit advection &/or neutral physics fluxes 
     real, dimension(isd:ied,jsd:jed,nk)   :: K33_implicit     ! m^2/sec vert-diffusivity from neutral diffusion 
     real, dimension(isd:ied,jsd:jed,nk)   :: radiation        ! radiation absorbed within a cell [W/m^2]
     real, dimension(isd:ied,jsd:jed)      :: stf              ! surface tracer flux [rho*m/sec*tracer concen]
     real, dimension(isd:ied,jsd:jed)      :: btf              ! bottom tracer flux [rho*m/sec*tracer concen]
     real, dimension(isd:ied,jsd:jed)      :: tpme             ! tracer concentration in precip-evap
     real, dimension(isd:ied,jsd:jed)      :: triver           ! tracer concentration in river(=runoff+calving) water  
     real, dimension(isd:ied,jsd:jed)      :: trunoff          ! tracer concentration in liquid runoff from land  
     real, dimension(isd:ied,jsd:jed)      :: tcalving       ! tracer concentration in frozen runoff from land (e.g., calving ice)
     real, dimension(isd:ied,jsd:jed) :: runoff_tracer_flux  ! tracer flux in liquid runoff (e.g., kg*degC/(m^2 s) for temp)  
     real, dimension(isd:ied,jsd:jed) :: calving_tracer_flux ! tracer flux in solid  runoff (e.g., kg*psu/(m^2 s)  for salinity)
     real, dimension(isd:ied,jsd:jed) :: riverdiffuse        ! sets where to enhance diff_cbt according to rivers
     real, dimension(isd:ied,jsd:jed) :: eta_smooth          ! tendency [tracer*(kg/m^3)*(m/s)] from eta_t smoother  
     real, dimension(isd:ied,jsd:jed) :: pbot_smooth         ! tendency [tracer*(kg/m^3)*(m/s)] from pbot_t smoother 

     ! variables for prather second order moment advection
     logical :: psom_limit                             ! controls whether a limiter is placed on the prather flux
     real, dimension(:,:,:), _ALLOCATABLE :: s0  _NULL ! zeroth moment (mean)
     real, dimension(:,:,:), _ALLOCATABLE :: sx  _NULL ! 1st moment in i direction
     real, dimension(:,:,:), _ALLOCATABLE :: sxx _NULL ! 2nd moment in i direction
     real, dimension(:,:,:), _ALLOCATABLE :: sy  _NULL ! 1st moment in j direction
     real, dimension(:,:,:), _ALLOCATABLE :: syy _NULL ! 2nd moment in j direction
     real, dimension(:,:,:), _ALLOCATABLE :: sz  _NULL ! 1st moment in k direction
     real, dimension(:,:,:), _ALLOCATABLE :: szz _NULL ! 2nd moment in k direction
     real, dimension(:,:,:), _ALLOCATABLE :: sxy _NULL ! 2nd moment for coupling the i and j directions
     real, dimension(:,:,:), _ALLOCATABLE :: sxz _NULL ! 2nd moment for coupling the i and k directions
     real, dimension(:,:,:), _ALLOCATABLE :: syz _NULL ! 2nd moment for coupling the j and k directions

     ! The following arrays are never allocated with MOM_STATIC_ARRAYS,
     ! but they need to be declared allocatable in order to compile. 
     real, _ALLOCATABLE, dimension(:,:,:)   :: fieldT     _NULL ! Total tracer concentration
     real, _ALLOCATABLE, dimension(:,:,:,:) :: sum_blob   _NULL ! tracer content [concentration*kg] from the L system
     real, _ALLOCATABLE, dimension(:,:,:)   :: tend_blob  _NULL ! blob contribution to dat*th_tendency

     real                                  :: conversion         ! conversion of dimensions  
     real                                  :: offset             ! offset in dimensions (e.g., Celsius to Kelvin)
     real                                  :: min_tracer         ! min acceptable value--model stopped if less
     real                                  :: max_tracer         ! max acceptable value--model stopped if greater
     real                                  :: min_range          ! min value used for calls to diagnostic manager 
     real                                  :: max_range          ! max value used for calls to diagnostic manager 
     real                                  :: min_tracer_limit   ! min value used to limit quicker and neutral fluxes
     real                                  :: max_tracer_limit   ! max value used to limit quicker and neutral fluxes 
     real                                  :: min_flux_range     ! min and max values used for flux diagnostics
     real                                  :: max_flux_range     ! min and max values used for flux diagnostics
     real                                  :: const_init_value   ! value used for constant tracer init
     logical                               :: const_init_tracer  ! false (default) if the tracer must exist in the restart file
                                                                 !   otherwise will initialize with const_init_value
     character(len=128)                    :: flux_units         ! units for the tracer flux
     character(len=128)                    :: restart_file       ! name for restart file
  end type ocean_prog_tracer_type
  
  type, public :: ocean_diag_tracer_type
     character(len=32)  :: name
     character(len=128) :: units, type
     character(len=128) :: longname
     real, dimension(isd:ied,jsd:jed,nk) :: field              ! tracer concentration at single time level 
     real :: conversion                                        ! conversion from model dimensions to others 
     real :: offset                                            ! offset in dimensions (e.g., Celsius to Kelvin)
     real :: min_tracer, max_tracer                            ! min and max acceptable values used to error check 
     real :: min_range, max_range                              ! min and max values used for diagnostics
     logical :: const_init_tracer                              ! false (default) if the tracer must exist in the restart file
                                                               !   otherwise will initialize with const_init_value
     character(len=128) :: restart_file                        ! name for restart file
     real               :: const_init_value                    ! to initialize tracer when constant_init_tracer
  end type ocean_diag_tracer_type


  type, public :: ocean_velocity_type
     logical                                 :: bmf_implicit    ! is true when time stepping bmf implicitly
     real, dimension(isd:ied,jsd:jed,nk,2,3) :: u               ! horz velocity (m/s) in i,j directions at 3 time levels
     real, dimension(isd:ied,jsd:jed,2)      :: smf             ! momentum flux per mass into ocean surface at uv point (N/m^2)
     real, dimension(isd:ied,jsd:jed,2)      :: smf_bgrid       ! momentum flux per mass into ocean surface at Bgrid uv point (N/m^2)
     real, dimension(isd:ied,jsd:jed,2)      :: smf_cgrid       ! momentum flux per mass into ocean surface at Cgrid u/v points (N/m^2)
     real, dimension(isd:ied,jsd:jed,2)      :: bmf             ! momentum flux per mass into ocean bottom  (N/m^2)
     real, dimension(isd:ied,jsd:jed)        :: ustar           ! surface friction velocity (m/s)
     real, dimension(isd:ied,jsd:jed)        :: gamma           ! dimensionful bottom drag coefficient (kg/(m^2 sec))
     real, dimension(isd:ied,jsd:jed)        :: langmuirfactor  ! dimensionless langmuir turbulence enhancement factor (non dimensional)
     real, dimension(isd:ied,jsd:jed)        :: u10             ! 10m wind speed (m/s)
     real, dimension(isd:ied,jsd:jed)        :: ustoke          ! x-dir surface stokes drift (m/s)
     real, dimension(isd:ied,jsd:jed)        :: vstoke          ! y-dir surface stokes drift (m/s)
     real, dimension(isd:ied,jsd:jed)        :: wavlen          ! wave length (m)
     real, dimension(isd:ied,jsd:jed)        :: cdbot_array     ! dimensionless static bottom drag coefficient
     real, dimension(isd:ied,jsd:jed)        :: current_wave_stress !wave-current bottom stress for sediment dynamics (N/m^2)
     real, dimension(isd:ied,jsd:jed)        :: rossby_radius   ! first baroclinic rossby radius (m)
     real, dimension(isd:ied,jsd:jed)        :: stokes_depth    ! depth scale (m) used for exponential decay of surface wave Stokes velocity
     real, dimension(isd:ied,jsd:jed,nk,2)   :: stokes_drift    ! Stokes drift velocity (m/s) from surface wave model 
     real, dimension(isd:ied,jsd:jed,nk,2)   :: stokes_force    ! Coriolis force from Stokes drift (N/m^2)
     real, dimension(isd:ied,jsd:jed,nk,2)   :: press_force     ! thickness*density weighted (i,j)-directed press force (N/m^2)
     real, dimension(isd:ied,jsd:jed,nk,2)   :: accel           ! thickness*density weighted (i,j)-directed acceleration (N/m^2)
     real, dimension(isd:ied,jsd:jed,nk,2)   :: vfrict_impl     ! thickness*density weighted vertical friction (N/m^2)
     real, dimension(isd:ied,jsd:jed,nk,2)   :: source          ! thickness*density weighted velocity source/sink (N/m^2)
     real, dimension(isd:ied,jsd:jed,nk,2)   :: wrkv            ! work array
     real, dimension(isd:ied,jsd:jed,nk,2,3) :: advection       ! thickness weighted tendency by velocity advect (N/m^2)
     real, dimension(isd:ied,jsd:jed,nk,2,3) :: coriolis        ! thickness weighted tendency by Cgrid coriolis force (N/m^2)
     real, dimension(isd:ied,jsd:jed,2)      :: lap_friction_bt ! friction just on barotropic velocity (N/m^2)
     real, dimension(isd:ied,jsd:jed,2)      :: bih_friction_bt ! friction just on barotropic velocity (N/m^2)
  end type ocean_velocity_type


  type, public :: ocean_external_mode_type  
     real, dimension(isd:ied,jsd:jed,3)   :: eta_t        ! surface height on tracer cell center (m) 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_u        ! surface height on velocity cell center  (m) 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_t_bar    ! surface height on tracer cell time avg over ext-mode time steps (m) 
     real, dimension(isd:ied,jsd:jed)     :: deta_dt      ! surface height time tendency on t-cell (m/s)

     real, dimension(isd:ied,jsd:jed,3)   :: pbot_t         ! bottom pressure on tracer cell center (Pa) 
     real, dimension(isd:ied,jsd:jed,3)   :: pbot_u         ! bottom pressure on velocity cell center (Pa) 
     real, dimension(isd:ied,jsd:jed,3)   :: anompb         ! pbot_t-rho0*grav*ht = anomalous bottom pressure (Pa)
     real, dimension(isd:ied,jsd:jed,3)   :: anompb_bar     ! pbot_t-rho0*grav*ht (Pa) time avg over ext-mode time steps 
     real, dimension(isd:ied,jsd:jed)     :: dpbot_dt       ! bottom pressure time tendency on tracer cell (Pa/s)

     real, dimension(isd:ied,jsd:jed)     :: patm_u           ! atmospheric pressure on velocity cell (Pa)
     real, dimension(isd:ied,jsd:jed,3)   :: patm_t           ! atmospheric pressure on tracer cell (Pa)
     real, dimension(isd:ied,jsd:jed)     :: patm_for_sea_lev ! atmospheric pressure on tracer cell (Pa) passed to coupler
     real, dimension(isd:ied,jsd:jed)     :: dpatm_dt         ! atmospheric pressure time tendency on tracer cell (Pa/s)

     real, dimension(isd:ied,jsd:jed,3)   :: conv_rho_ud_t  ! convergence of sum_k(rho*dzt*u) on T-cell (kg/m^3)*(m/s)
     real, dimension(isd:ied,jsd:jed,2,3) :: udrho          ! vertically integrated and rho weighted horz velocity (kg/m^3)*(m^2/s)
     real, dimension(isd:ied,jsd:jed,2)   :: forcing_bt     ! depth integrated time change of velocity (without coriolis)
     real, dimension(isd:ied,jsd:jed)     :: ps             ! surface pressure (pressure at z=0 due to eta_t) at time tau
     real, dimension(isd:ied,jsd:jed,2)   :: grad_ps        ! rho0r * surface pressure gradient at time tau
     real, dimension(isd:ied,jsd:jed,2)   :: grad_anompb    ! gradient of anompb (Pa/m) on U-point

     real, dimension(isd:ied,jsd:jed,2)   :: press_force    ! pressure force/area (N/m^2) updated each barotropic time step 

     real, dimension(isd:ied,jsd:jed)     :: source         ! density weighted vertical integral of mass source (kg/m^2/sec)
     real, dimension(isd:ied,jsd:jed)     :: eta_smooth     ! (kg/m^3)*(m/s) from eta_t smoother  
     real, dimension(isd:ied,jsd:jed)     :: pbot_smooth    ! (kg/m^3)*(m/s) from pbot_t smoother 

     real, dimension(isd:ied,jsd:jed,3)   :: eta_nonbouss   ! diagnosed eta (m) including steric effects
     real, dimension(isd:ied,jsd:jed,3)   :: eta_nonsteric  ! diagnosed piece of eta_nonbouss (m) from non-steric effects 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_steric     ! diagnosed eta_nonbouss (m) from steric effects 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_dynamic    ! diagnosed eta_nonbouss (m) from "dynamic" effects 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_water      ! diagnosed eta_nonbouss (m) from water forcing at boundaries
     real, dimension(isd:ied,jsd:jed,3)   :: eta_source     ! diagnosed eta_nonbouss (m) from source term 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_surf_temp  ! diagnosed eta_nonbouss (m) arising from surface temp flux
     real, dimension(isd:ied,jsd:jed,3)   :: eta_surf_salt  ! diagnosed eta_nonbouss (m) arising from surface salt flux
     real, dimension(isd:ied,jsd:jed,3)   :: eta_surf_water ! diagnosed eta_nonbouss (m) arising from surface water flux
     real, dimension(isd:ied,jsd:jed,3)   :: eta_bott_temp  ! diagnosed eta_nonbouss (m) arising from bottom temp flux

     ! The following arrays are never allocated with MOM_STATIC_ARRAYS,
     ! but they need to be declared allocatable in order to compile. 
     real, _ALLOCATABLE, dimension(:,:)   :: conv_blob  _NULL ! water column divergence of the L system


  end type ocean_external_mode_type


  ! types used by neutral physics 
  type, public :: tracer_2d_type
     real, dimension(isd:ied,jsd:jed)     :: field
  end type tracer_2d_type

  type, public :: tracer_3d_0_nk_type
    real, dimension(isd:ied,jsd:jed,0:nk) :: field
  end type tracer_3d_0_nk_type

  type, public :: tracer_3d_1_nk_type
    real, dimension(isd:ied,jsd:jed,nk)   :: field
  end type tracer_3d_1_nk_type


  ! for gotm vertical mixing scheme 
  type, public :: ocean_gotm_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname
     real               :: min_value   ! minimum value
     real               :: max_value   ! maximum value
     real, dimension(isd:ied,jsd:jed,nk,2) :: field  ! scalar at 2 time levels
  end type ocean_gotm_type

  ! advection tendency for GOTM scalar fields 
  type, public :: advect_gotm_type              
    real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: field
  end type advect_gotm_type


#else
!############################################################################################
! not MOM_STATIC_ARRAYS

  type, public :: ocean_thickness_type
     integer :: method       ! energetic or finite volume  

     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzt   _NULL ! E system contribution to rho_dztT (3 time levels)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzten _NULL ! rho*dz at east/north face of T-cell
     real, dimension(:,:,:),   _ALLOCATABLE :: rho_dztr  _NULL ! 1.0/(rho*dzt) at time taup1 (E system)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztL  _NULL ! L system contribution to rho_dztT (3 time levels)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztT  _NULL ! rho(kg/m^3)*thickness (m) of T cell at 3 times (total)

     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzu   _NULL ! E system contribution to rho_dzuT (3 time levels)
     real, dimension(:,:,:),   _ALLOCATABLE :: rho_dzur  _NULL ! 1.0/(rho*dzu) at time taup1 (E system)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuL  _NULL ! L system contribution to rho_dzuT (3 time levels)
     real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuT  _NULL ! rho (kg/m^3) * thickness (m) of U cell at 3 times (total)

     real, dimension(:,:,:),   _ALLOCATABLE :: rho_dzt_tendency _NULL ! rho_dzt tendency (kg/m^3)*(m/s)

     real, dimension(:,:),     _ALLOCATABLE :: sea_lev _NULL ! eta_t + patm/(rho0*grav) - eta_geoid - eta_tide (m) at time taup1 for coupler 

     real, dimension(:,:,:),   _ALLOCATABLE :: dzt    _NULL ! E system contribution to dztT
     real, dimension(:,:,:,:), _ALLOCATABLE :: dzten  _NULL ! E system contribution to dzt at east/north face of T-cell
     real, dimension(:,:,:),   _ALLOCATABLE :: dztL   _NULL ! L system contribution to dztT
     real, dimension(:,:,:,:), _ALLOCATABLE :: dztT   _NULL ! thickness (m) of T cell at time tau/taup1

     real, dimension(:,:,:),   _ALLOCATABLE :: dzu    _NULL ! E system contribution to dzuT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzuL   _NULL ! L system contribution to dzuT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzuT   _NULL ! thickness (m) of U cell at time tau/taup1

     real, dimension(:,:,:),   _ALLOCATABLE :: dzwt   _NULL ! E system contribution to dzwtT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwtL  _NULL ! L system contribution to dzwtT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwtT  _NULL ! vertical distance (m) between T points at tau/taup1

     real, dimension(:,:,:),   _ALLOCATABLE :: dzwu   _NULL ! E system contribution to dzwuT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwuL  _NULL ! L system contribution to dzwuT
     real, dimension(:,:,:),   _ALLOCATABLE :: dzwuT  _NULL ! vertical distance (m) between U points at tau/taup1

     real, dimension(:,:,:),   _ALLOCATABLE :: dztup  _NULL ! E system contribution to dztupT
     real, dimension(:,:,:),   _ALLOCATABLE :: dztupL _NULL ! L system contribution to dztupT
     real, dimension(:,:,:),   _ALLOCATABLE :: dztupT _NULL ! distance (m) from T-cell point to top of T-cell 

     real, dimension(:,:,:),   _ALLOCATABLE :: dztlo  _NULL ! E system contribution to dztloT
     real, dimension(:,:,:),   _ALLOCATABLE :: dztloL _NULL ! L system contribution to dztloT
     real, dimension(:,:,:),   _ALLOCATABLE :: dztloT _NULL ! distance (m) from T-cell point to bottom of T-cell

     real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zt  _NULL ! vert distance (m) from z=0 to T-point
     real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zwt _NULL ! vert distance (m) from z=0 to bottom of T-cell
     real, dimension(:,:,:),   _ALLOCATABLE :: depth_zt     _NULL ! vert distance (m) from column top to T-point
     real, dimension(:,:,:),   _ALLOCATABLE :: depth_zwt    _NULL ! vert distance (m) from column top to T-bottom
     real, dimension(:,:,:),   _ALLOCATABLE :: depth_zu     _NULL ! vert distance (m) from column top to U-point
     real, dimension(:,:,:),   _ALLOCATABLE :: depth_zwu    _NULL ! vert distance (m) from column top to U-bottom

     real, dimension(:,:,:),   _ALLOCATABLE :: depth_st    _NULL ! s-distance to T cell grid point
     real, dimension(:,:,:),   _ALLOCATABLE :: depth_swt   _NULL ! s-distance to T cell grid bottom
     real, dimension(:,:,:),   _ALLOCATABLE :: dst         _NULL ! s-increment of T-cell at time tau 
     real, dimension(:,:,:),   _ALLOCATABLE :: dswt        _NULL ! s-increment of T-point at time tau 
     real, dimension(:,:,:),   _ALLOCATABLE :: dzt_dst     _NULL ! T-cell specific thickness (m/s-coordinate)

     real, dimension(:,:,:),   _ALLOCATABLE :: dstlo       _NULL ! s-distance (s-units) from T-cell point to top of T-cell 
     real, dimension(:,:,:),   _ALLOCATABLE :: dstup       _NULL ! s-distance (s-units) from T-cell point to bottom of T-cell 

     real, dimension(:,:),     _ALLOCATABLE :: pbot0       _NULL ! reference bottom pressure (Pa) 
     real, dimension(:,:),     _ALLOCATABLE :: pbot0r      _NULL ! inverse of reference bottom pressure (1/Pa) 

     real, dimension(:,:,:),   _ALLOCATABLE :: mass_source _NULL ! mass source (kg/m^3)*(m/sec)  
     real, dimension(:,:),     _ALLOCATABLE :: blob_source _NULL ! mass source (sink) for the L system (E system)
     real, dimension(:,:,:),   _ALLOCATABLE :: mass_u      _NULL ! mass per area (kg/m^2) in a velocity column
     real, dimension(:,:,:),   _ALLOCATABLE :: mass_en     _NULL ! vertical sum rho_dzten (kg/m^2) 
     real, dimension(:,:,:),   _ALLOCATABLE :: mass_uT     _NULL ! mass per area (kg/m^2) in a velocity column
     real, dimension(:,:,:),   _ALLOCATABLE :: thicku      _NULL ! thickness on U-cell (metre) from z=eta_u to z=-H.
     real, dimension(:,:,:),   _ALLOCATABLE :: thicken     _NULL ! sum of dzten (metre)

  end type ocean_thickness_type


  type, public :: ocean_grid_type
     character(len=32) :: name

     ! geometry and topology and rotation
     logical                            :: cyclic_x          ! true if domain is cyclic in the i direction
     logical                            :: cyclic_y          ! true if domain is cyclic in the j direction
     logical                            :: tripolar          ! folded connectivity at "top" row w/i bipolar Arctic 
     logical                            :: mosaic            ! true when using a mosaic grid 
     logical                            :: beta_plane        ! beta plane Cartesian 
     logical                            :: f_plane           ! f-plane Cartesian
     real                               :: f_plane_latitude  ! latitude where f_plane is centered  
     real, dimension(:,:), _ALLOCATABLE :: f        _NULL ! coriolis parameter at u-cell points (sec^-1)
     real, dimension(:,:), _ALLOCATABLE :: fstar    _NULL ! horizontal coriolis parameter at u-cell points (sec^-1)
     real, dimension(:,:), _ALLOCATABLE :: beta     _NULL ! df/dy at u-cell points (1/(sec*m))     
     real, dimension(:,:), _ALLOCATABLE :: beta_eff _NULL ! df/dy plus topographic beta at u-cell (1/(sec*m))

     ! vertical grid information (time independent) 
     integer                               :: nk           ! number of vertical grid points 
     integer, dimension(:,:), _ALLOCATABLE :: kmt    _NULL ! number of t-levels
     integer, dimension(:,:), _ALLOCATABLE :: kmu    _NULL ! number of u-levels

     real,    dimension(:),   _ALLOCATABLE :: zt     _NULL ! distance from surface to grid point in level k (m) 
     real,    dimension(:),   _ALLOCATABLE :: zw     _NULL ! distance from surface down to bottom of level k (m) 
     real,    dimension(:),   _ALLOCATABLE :: dzt    _NULL ! initial vertical resolution of T or U grid cells (m) 
     real,    dimension(:),   _ALLOCATABLE :: dztlo  _NULL ! z-distance (m) from T-cell point to bottom of T-cell
     real,    dimension(:),   _ALLOCATABLE :: dztup  _NULL ! z-distance (m) from T-cell point to top of T-cell
     real,    dimension(:),   _ALLOCATABLE :: dzw    _NULL ! initial vertical resolution of W grid cells (m) 
     real,    dimension(:),   _ALLOCATABLE :: dzwr   _NULL ! reciprocal of dzw (W cell vertical resolution)

     real,    dimension(:),   _ALLOCATABLE :: st     _NULL ! s-distance from surface to grid point in level k  
     real,    dimension(:),   _ALLOCATABLE :: sw     _NULL ! s-distance from surface down to bottom of level k  
     real,    dimension(:),   _ALLOCATABLE :: dst    _NULL ! initial s-vertical resolution of T or U grid cells 
     real,    dimension(:),   _ALLOCATABLE :: dstlo  _NULL ! s-distance (s-units) from T-cell point to bottom of T-cell
     real,    dimension(:),   _ALLOCATABLE :: dstup  _NULL ! s-distance (s-units) from T-cell point to top of T-cell
     real,    dimension(:),   _ALLOCATABLE :: dsw    _NULL ! initial s-vertical resolution of W grid cells

     real,    dimension(:,:), _ALLOCATABLE :: fracdz _NULL ! fractional distance between grid point & cell top/bot 
     real,    dimension(:,:), _ALLOCATABLE :: ht     _NULL ! depth to bottom of ocean (m) on t-cells
     real,    dimension(:,:), _ALLOCATABLE :: htr    _NULL ! inverse depth to bottom of ocean (m^-1) on t-cells  
     real,    dimension(:,:), _ALLOCATABLE :: hu     _NULL ! depth to bottom of ocean (m) on u-cells  
     real,    dimension(:,:), _ALLOCATABLE :: dht_dx _NULL ! d(ht)/dx on u-cells (m/m) 
     real,    dimension(:,:), _ALLOCATABLE :: dht_dy _NULL ! d(ht)/dy on u-cells (m/m) 
     real,    dimension(:,:), _ALLOCATABLE :: gradH  _NULL ! sqrt(dht_dx**2+dht_dyx**2) on u-cells (m/m)

     ! horizontal grid information (time independent)
     integer                            :: ni, nj           ! global points in the two horizontal directions
     real, dimension(:,:), _ALLOCATABLE :: xt         _NULL ! longitude of the T grid points in degrees
     real, dimension(:,:), _ALLOCATABLE :: xu         _NULL ! longitude of the U grid points in degrees
     real, dimension(:,:), _ALLOCATABLE :: yt         _NULL ! latitude of the T grid points in degrees
     real, dimension(:,:), _ALLOCATABLE :: yu         _NULL ! latitude of the U grid points in degrees
     real, dimension(:,:), _ALLOCATABLE :: phiu       _NULL ! latitude of U grid point in radians
     real, dimension(:,:), _ALLOCATABLE :: phit       _NULL ! latitude of T grid point in radians
     real, dimension(:,:), _ALLOCATABLE :: h1t        _NULL ! metric factors in i grid direction
     real, dimension(:,:), _ALLOCATABLE :: h1u        _NULL ! metric factors in j grid jirection
     real, dimension(:,:), _ALLOCATABLE :: h2t        _NULL ! metric factors in i grid direction
     real, dimension(:,:), _ALLOCATABLE :: h2u        _NULL ! metric factors in j grid jirection
     real, dimension(:,:), _ALLOCATABLE :: dh2dx      _NULL ! (1/delta_y)*d(delta_y)/dx (1/m)
     real, dimension(:,:), _ALLOCATABLE :: dh1dy      _NULL ! (1/delta_x)*d(delta_x)/dy (1/m)
     real, dimension(:,:), _ALLOCATABLE :: dxt        _NULL ! longitudinal width of T-cells at grid point (m)
     real, dimension(:,:), _ALLOCATABLE :: dxu        _NULL ! longitudinal width of U-cells at grid point (m)
     real, dimension(:,:), _ALLOCATABLE :: dyt        _NULL ! latitudinal width of T-cells at grid point (m) 
     real, dimension(:,:), _ALLOCATABLE :: dyu        _NULL ! latitudinal width of U-cells at grid point (m)
     real, dimension(:,:), _ALLOCATABLE :: dau        _NULL ! area of U-cells (m^2)
     real, dimension(:,:), _ALLOCATABLE :: dat        _NULL ! area of T-cells (m^2)
     real, dimension(:,:), _ALLOCATABLE :: dat_frac   _NULL ! fraction of total wet area occuped by a T-cell (dimensionless)
     real, dimension(:,:), _ALLOCATABLE :: dxte       _NULL ! i-width between i+1 and i points in T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxtn       _NULL ! i-width of north face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dyte       _NULL ! j-width of east  face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dytn       _NULL ! j-width between j+1 and j points in T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxue       _NULL ! i-width between i+1 and i points in U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxun       _NULL ! i-width of north face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dyue       _NULL ! j-width of east  face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dyun       _NULL ! j-width between j+1 and j points in U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: datnr      _NULL ! reciprocal area at north face of T-cell
     real, dimension(:,:), _ALLOCATABLE :: dater      _NULL ! reciprocal area at east face of T-cell
     real, dimension(:,:), _ALLOCATABLE :: dun        _NULL ! width from grid point to north face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dus        _NULL ! width from grid point to south face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: duw        _NULL ! width from grid point to west  face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: due        _NULL ! width from grid point to east  face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dtn        _NULL ! width from grid point to north face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dts        _NULL ! width from grid point to south face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dtw        _NULL ! width from grid point to west  face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dte        _NULL ! width from grid point to east  face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxtr       _NULL ! 1/dxt
     real, dimension(:,:), _ALLOCATABLE :: dxur       _NULL ! 1/dxu
     real, dimension(:,:), _ALLOCATABLE :: dytr       _NULL ! 1/dyt
     real, dimension(:,:), _ALLOCATABLE :: dyur       _NULL ! 1/dyu
     real, dimension(:,:), _ALLOCATABLE :: daur       _NULL ! 1/[area of U-cells (m^2)]
     real, dimension(:,:), _ALLOCATABLE :: datr       _NULL ! 1/[area of T-cells (m^2)]
     real, dimension(:,:), _ALLOCATABLE :: dxter      _NULL ! 1/dxte
     real, dimension(:,:), _ALLOCATABLE :: dyter      _NULL ! 1/dyte
     real, dimension(:,:), _ALLOCATABLE :: dxtnr      _NULL ! 1/dxtn
     real, dimension(:,:), _ALLOCATABLE :: dytnr      _NULL ! 1/dytn
     real, dimension(:,:), _ALLOCATABLE :: dxuer      _NULL ! 1/dxue
     real, dimension(:,:), _ALLOCATABLE :: dxunr      _NULL ! 1/dxun
     real, dimension(:,:), _ALLOCATABLE :: dyuer      _NULL ! 1/dyue
     real, dimension(:,:), _ALLOCATABLE :: dyunr      _NULL ! 1/dyun
     real, dimension(:,:), _ALLOCATABLE :: dyue_dxuer _NULL ! dyue/dxue 
     real, dimension(:,:), _ALLOCATABLE :: dxun_dyunr _NULL ! dxun/dyun
     real, dimension(:,:), _ALLOCATABLE :: dxt_dxter  _NULL ! dxt/dxte
     real, dimension(:,:), _ALLOCATABLE :: dxte_dxtr  _NULL ! dxte/dxt
     real, dimension(:,:), _ALLOCATABLE :: dxt_dxtnr  _NULL ! dxt/dxtn
     real, dimension(:,:), _ALLOCATABLE :: dxtn_dxtr  _NULL ! dxtn/dxt
     real, dimension(:,:), _ALLOCATABLE :: dyt_dyter  _NULL ! dyt/dyte
     real, dimension(:,:), _ALLOCATABLE :: dyte_dytr  _NULL ! dyte/dyt
     real, dimension(:,:), _ALLOCATABLE :: dyt_dytnr  _NULL ! dyt/dytn
     real, dimension(:,:), _ALLOCATABLE :: dytn_dytr  _NULL ! dytn/dyt

     ! land/sea masks 
     real, dimension(:,:),     _ALLOCATABLE :: obc_tmask   _NULL ! land/sea mask for T cell diagnostics with OBC 
     real, dimension(:,:),     _ALLOCATABLE :: obc_umask   _NULL ! land/sea mask for U cell diagnostics with OBC
     real, dimension(:,:,:),   _ALLOCATABLE :: mask        _NULL ! land/sea tmask or umask depending on C/B grids 
     real, dimension(:,:,:),   _ALLOCATABLE :: tmask       _NULL ! land/sea mask for T cells based on s-coordinate
     real, dimension(:,:,:),   _ALLOCATABLE :: umask       _NULL ! land/sea mask for U cells based on s-coordinate
     real, dimension(:,:,:,:), _ALLOCATABLE :: tmasken     _NULL ! land/sea mask for east/north of t-cell based on s-coordinate
     real, dimension(:,:,:),   _ALLOCATABLE :: tmask_depth _NULL ! based on depth-based vert_coordinate
     real, dimension(:,:,:),   _ALLOCATABLE :: umask_depth _NULL ! based on depth-based vert_coordinate

     ! grid areas and volumes 
     real                              :: tcellv          ! initial T cell volume m^3 (entire ocean) 
     real                              :: ucellv          ! initial U cell volume m^3 (entire ocean)
     real                              :: tcellsurf       ! T cell surface area (k=1) (bitwise reproducible)
     real                              :: ucellsurf       ! U cell surface area (k=1) (bitwise reproducible)
     real, dimension(:), _ALLOCATABLE  :: tcella   _NULL  ! T cell surface area m^2 (entire ocean)
     real, dimension(:), _ALLOCATABLE  :: ucella   _NULL  ! U cell surface area m^2 (entire ocean)

     ! model grid information
     integer :: wet_t_points   ! total number of wet tracer points 
     integer :: wet_u_points   ! total number of wet velocity points 
     integer :: total_t_points ! total number of wet or dry tracer points  

     ! sine and cosine of rotation angles (clockwise) of velocity for tripolar 
     real, dimension(:,:), _ALLOCATABLE :: sin_rot  _NULL  
     real, dimension(:,:), _ALLOCATABLE :: cos_rot  _NULL  

     ! 1-d grid coordinates for COARDS NetCDF files
     real, dimension(:), _ALLOCATABLE  :: grid_x_t _NULL 
     real, dimension(:), _ALLOCATABLE  :: grid_y_t _NULL 
     real, dimension(:), _ALLOCATABLE  :: grid_x_u _NULL 
     real, dimension(:), _ALLOCATABLE  :: grid_y_u _NULL 

     ! axes id for diagnostic manager 
     integer, dimension(3)  :: tracer_axes        
     integer, dimension(3)  :: vel_axes_uv         
     integer, dimension(3)  :: vel_axes_u
     integer, dimension(3)  :: vel_axes_v         
     integer, dimension(3)  :: vel_axes_wu    
     integer, dimension(3)  :: vel_axes_wt    
     integer, dimension(3)  :: tracer_axes_wt 
     integer, dimension(3)  :: tracer_axes_flux_x  
     integer, dimension(3)  :: tracer_axes_flux_y  
     integer, dimension(3)  :: vel_axes_flux_x    
     integer, dimension(3)  :: vel_axes_flux_y    

     ! axes id for diagnostic manager, appropriate
     ! when with to remap native vertical fields 
     ! to depth or pressure levels.  Appropropriate
     ! when vert_coordinate == ZSIGMA or PSIGMA, or other 
     ! vertical coordinats whose iso-surfaces are 
     ! not quasi-horizontal 
     integer, dimension(3)  :: tracer_axes_depth        
     integer, dimension(3)  :: vel_axes_uv_depth
     integer, dimension(3)  :: vel_axes_u_depth
     integer, dimension(3)  :: vel_axes_v_depth
     integer, dimension(3)  :: vel_axes_wu_depth    
     integer, dimension(3)  :: vel_axes_wt_depth
     integer, dimension(3)  :: tracer_axes_wt_depth
     integer, dimension(3)  :: tracer_axes_flux_x_depth
     integer, dimension(3)  :: tracer_axes_flux_y_depth
     integer, dimension(3)  :: vel_axes_flux_x_depth
     integer, dimension(3)  :: vel_axes_flux_y_depth


  end type ocean_grid_type
  
  type, public :: ocean_domain_type
     type(domain2d) :: domain2d           ! fms variable, used by mpp routines
     integer :: isc, iec, jsc, jec        ! computational domain indices 
     integer :: isd, ied, jsd, jed        ! local indices, consistent with domain2d
     integer :: isg, ieg, jsg, jeg        ! global indices
     integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)
     integer :: xhalo, yhalo              ! halo sizes 
     integer :: xflags, yflags, layout(2) ! options to mpp_define_domains
     integer :: ioff, joff                ! offset to get absolute indices when MOM_STATIC_ARRAYS (0 otherwise)
     integer :: io_layout(2)              ! options to define io_domain.
     integer :: x_cyclic_offset           ! offset applied to x-direction cyclic boundary condition
     integer :: y_cyclic_offset           ! offset applied to y-direction cyclic boundary condition
     logical, pointer :: maskmap(:,:) =>NULL() ! option to mpp_define_domains
  end type ocean_domain_type

  type, public :: ocean_time_type
     type(time_type) :: model_time     ! fms variable
     type(time_type) :: time_step      ! ocean tracer timestep
     type(time_type) :: Time_init      ! time of initial conditions 
     integer :: calendar               ! calendar type defined by time_manager_mod
     logical :: init                   ! true at beginning of run (initial condition time)
     integer :: itt                    ! timestep counter measured relative to time at restart 
     integer :: itt0                   ! timestep counter measured relative to initial condition time
     integer :: taum1, tau, taup1      ! time level indices
     integer :: tau_m2, tau_m1, tau_m0 ! time level indices for Adams-Bashforth velocity advection and coriolis 
  end type ocean_time_type

  type, public :: ocean_adv_vel_type
     real, _ALLOCATABLE, dimension(:,:,:)  :: uhrho_et   _NULL ! rho_dzu weight advect vel (kg/(m*s)) on T-cell i-face
     real, _ALLOCATABLE, dimension(:,:,:)  :: vhrho_nt   _NULL ! rho_dzu weight advect vel (kg/(m*s)) on T-cell j-face
     real, _ALLOCATABLE, dimension(:,:,:)  :: uhrho_eu   _NULL ! remapped uhrho_et (kg/(m*s)) onto i-face of U-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: vhrho_nu   _NULL ! remapped vhrho_nt (kg/(m*s)) onto j-face of U-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: wrho_bt    _NULL ! rho weight (kg/(m^2*s)) vert advect vel on T-bottom
     real, _ALLOCATABLE, dimension(:,:,:)  :: wrho_bu    _NULL ! remapped wrho_bt onto U-bottom
     real, _ALLOCATABLE, dimension(:,:,:)  :: diverge_t  _NULL ! divergence on T-cell of horiz momentum per mass (kg/m3)*(m/s)
     real, _ALLOCATABLE, dimension(:,:,:)  :: diverge_u  _NULL ! divergence on U-cell of horiz momentum per mass (kg/m3)*(m/s)
  end type ocean_adv_vel_type


  type, public ::  ocean_density_type
     logical                                :: use_teos10              ! for using the TEOS-2010 equation of state recommendations 
     real, _ALLOCATABLE, dimension(:,:,:,:) :: rho               _NULL ! in situ density (kg/m^3) at time levels
     real, _ALLOCATABLE, dimension(:,:,:,:) :: rho_salinity      _NULL ! salinity used in density calculations (psu or g/kg) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: rho_dztr_tau      _NULL ! rho_dztr at time tau 1/(kg/m^3) (for diagnostic uses) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: rhoT              _NULL ! combined L and E in situ density (kg/m^3)
     real, _ALLOCATABLE, dimension(:,:,:)   :: rho_fresh         _NULL ! in situ fresh water density (kg/m^3) at tau
     real, _ALLOCATABLE, dimension(:,:,:)   :: pressure_at_depth _NULL ! hydrostatic pressure (including patm)
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodT            _NULL ! partial rho wrt theta (kg/(m^3 C)
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodS            _NULL ! partial rho wrt salinity (kg/(m3 psu)
     real, _ALLOCATABLE, dimension(:,:,:)   :: dpotrhodT         _NULL ! partial potrho wrt theta (kg/(m^3 C)
     real, _ALLOCATABLE, dimension(:,:,:)   :: dpotrhodS         _NULL ! partial potrho wrt salinity (kg/(m3 psu)
     real, _ALLOCATABLE, dimension(:,:,:)   :: dpotrhodP         _NULL ! partial potrho wrt pressure (kg/(m3 psu)
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodP            _NULL ! partial rho wrt pressure (kg/(m3 Pa)
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodz_wt         _NULL ! d(neutral rho)/dz (kg/m^4) at W-point
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodz_zt         _NULL ! d(neutral rho)/dz (kg/m^4) at T-point
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodx_zt         _NULL ! d(neutral rho)/dx (kg/m^4) at T-point
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhody_zt         _NULL ! d(neutral rho)/dy (kg/m^4) at T-point
     real, _ALLOCATABLE, dimension(:,:,:)   :: drhodz_diag       _NULL ! regularized drhodz_zt for diagnostics 
     real, _ALLOCATABLE, dimension(:,:,:)   :: dTdz_zt           _NULL ! partial theta wrt z  (C/m) at T-point
     real, _ALLOCATABLE, dimension(:,:,:)   :: dSdz_zt           _NULL ! partial salinity wrt z  (psu/m) at T-point
     real, _ALLOCATABLE, dimension(:,:,:)   :: potrho            _NULL ! potential density (kg/m^3)
     real, _ALLOCATABLE, dimension(:,:,:)   :: neutralrho        _NULL ! neutral density (kg/m^3)
     real, _ALLOCATABLE, dimension(:,:,:)   :: watermass_factor      _NULL ! ratio |grad nrho|/|grad local ref potrho|*/delta(gamma)
     real, _ALLOCATABLE, dimension(:,:,:)   :: stratification_factor _NULL ! ratio |grad nrho|/|grad local ref potrho|*/delta(gamma)
     real, _ALLOCATABLE, dimension(:,:)     :: mld_subduction        _NULL ! depth mixed layer base (m) for subduction diagnostics 
     real, _ALLOCATABLE, dimension(:)       :: theta_ref             _NULL ! partition vertical into theta classes 
     real, _ALLOCATABLE, dimension(:)       :: theta_bounds          _NULL ! bounds for theta classes 
     real, _ALLOCATABLE, dimension(:)       :: potrho_ref            _NULL ! partition vertical into potrho classes 
     real, _ALLOCATABLE, dimension(:)       :: potrho_bounds         _NULL ! bounds for potrho classes 
     real, _ALLOCATABLE, dimension(:)       :: neutralrho_ref        _NULL ! partition vertical into neutral density classes 
     real, _ALLOCATABLE, dimension(:)       :: neutralrho_bounds     _NULL ! bounds for neutral density classes 
     integer, dimension(3)                  :: potrho_axes             ! axis ids for diagnosing potential density 
     integer, dimension(3)                  :: potrho_axes_flux_x          ! axis ids for diagnosing x-flux
     integer, dimension(3)                  :: potrho_axes_flux_y          ! axis ids for diagnosing y-flux
     integer, dimension(3)                  :: neutralrho_axes         ! axis ids for diagnosing neutral density 
     integer, dimension(3)                  :: neutralrho_axes_flux_x      ! axis ids for diagnosing x-flux 
     integer, dimension(3)                  :: neutralrho_axes_flux_y      ! axis ids for diagnosing y-flux 
     integer, dimension(3)                  :: theta_axes              ! axis ids for potential temperature 
     integer, dimension(3)                  :: theta_axes_flux_x           ! axis ids for diagnosing x-flux 
     integer, dimension(3)                  :: theta_axes_flux_y           ! axis ids for diagnosing y-flux 
  end type ocean_density_type
  
  type, public :: ocean_prog_tracer_type
     character(len=32)  :: name
     character(len=128) :: units
     character(len=128) :: type
     character(len=128) :: longname

     logical :: use_only_advection      ! for testing purposes, evolve using ONLY advection
     logical :: neutral_physics_limit   ! revert neutral physics to horz diffusion where tracer out of bounds
     logical :: complete                ! to determine if ready to do mpp updates

     integer :: sfc_flux_id=-1          ! index for time_interp_external
     integer :: horz_advect_scheme=-1   ! id for horizontal advection scheme
     integer :: vert_advect_scheme=-1   ! id for vertical advection scheme
     integer :: id_obc                  ! id to identify tracer in OBC-subroutines
     integer :: ppm_hlimiter=1          ! Limiter for use with PPM in horizontal
     integer :: ppm_vlimiter=1          ! Limiter for use with PPM in vertical
     integer :: mdt_scheme=4            ! Version of Multi-Dim. Modified Daru & Tenaud (MDMDT)

     real, _ALLOCATABLE, dimension(:,:,:,:) :: field          _NULL ! tracer concentration at 3 time levels (E system)
     real, _ALLOCATABLE, dimension(:,:,:)   :: fieldT         _NULL ! Total tracer concentration
     real, _ALLOCATABLE, dimension(:,:,:)   :: th_tendency    _NULL ! thickness weighted tracer tendency
     real, _ALLOCATABLE, dimension(:,:,:)   :: tendency       _NULL ! for diagnostics: tendency concentration [concentration/sec]
     real, _ALLOCATABLE, dimension(:,:,:)   :: source         _NULL ! tracer source [=tracer concentration per time]
     real, _ALLOCATABLE, dimension(:,:)     :: eta_smooth     _NULL ! tendency from eta_t smoother  
     real, _ALLOCATABLE, dimension(:,:)     :: pbot_smooth    _NULL ! tendency from pbot_t smoother  

     real, _ALLOCATABLE, dimension(:,:,:)   :: wrk1                _NULL ! work array
     real, _ALLOCATABLE, dimension(:,:,:)   :: tmask_limit         _NULL ! to limit advective and/or neutral physics flux
     real, _ALLOCATABLE, dimension(:,:,:)   :: K33_implicit        _NULL ! m^2/sec vert-diffusivity from neutral diffusion 
     real, _ALLOCATABLE, dimension(:,:,:)   :: radiation           _NULL ! radiation absorbed within a cell [W/m^2]
     real, _ALLOCATABLE, dimension(:,:)     :: stf                 _NULL ! surface tracer flux [rho*m/sec*tracer concen]
     real, _ALLOCATABLE, dimension(:,:)     :: btf                 _NULL ! bottom tracer flux [rho*m/sec*tracer concen]
     real, _ALLOCATABLE, dimension(:,:)     :: tpme                _NULL ! tracer concentration in precip-evap
     real, _ALLOCATABLE, dimension(:,:)     :: triver              _NULL ! tracer concentration in river(=runoff+calving) 
     real, _ALLOCATABLE, dimension(:,:)     :: trunoff             _NULL ! tracer concentration in river runoff 
     real, _ALLOCATABLE, dimension(:,:)     :: tcalving            _NULL ! tracer concentration in calving lang ice
     real, _ALLOCATABLE, dimension(:,:)     :: runoff_tracer_flux  _NULL ! flux in liquid runoff (e.g., kg*degC/(m^2 s) for temp)
     real, _ALLOCATABLE, dimension(:,:)     :: calving_tracer_flux _NULL ! flux in solid runoff (e.g., kg*psu/(m^2 s) for salt)
     real, _ALLOCATABLE, dimension(:,:)     :: riverdiffuse        _NULL ! where to enhance diff_cbt according to rivers
     real, _ALLOCATABLE, dimension(:,:)     :: flux_int            _NULL ! integrated sfc tracer flux for diagnostics
     real, _ALLOCATABLE, dimension(:,:,:,:) :: sum_blob            _NULL ! tracer content [concentration*kg] from the L system
     real, _ALLOCATABLE, dimension(:,:,:)   :: tend_blob           _NULL ! blob contribution to dat*th_tendency

     ! variables for prather second order moment advection
     logical :: psom_limit                             ! controls whether a limiter is placed on the prather flux
     real, dimension(:,:,:), _ALLOCATABLE :: s0  _NULL ! zeroth moment (mean)
     real, dimension(:,:,:), _ALLOCATABLE :: sx  _NULL ! 1st moment in i direction
     real, dimension(:,:,:), _ALLOCATABLE :: sxx _NULL ! 2nd moment in i direction
     real, dimension(:,:,:), _ALLOCATABLE :: sy  _NULL ! 1st moment in j direction
     real, dimension(:,:,:), _ALLOCATABLE :: syy _NULL ! 2nd moment in j direction
     real, dimension(:,:,:), _ALLOCATABLE :: sz  _NULL ! 1st moment in k direction
     real, dimension(:,:,:), _ALLOCATABLE :: szz _NULL ! 2nd moment in k direction
     real, dimension(:,:,:), _ALLOCATABLE :: sxy _NULL ! 2nd moment for coupling the i and j directions
     real, dimension(:,:,:), _ALLOCATABLE :: sxz _NULL ! 2nd moment for coupling the i and k directions
     real, dimension(:,:,:), _ALLOCATABLE :: syz _NULL ! 2nd moment for coupling the j and k directions

     type(obc_flux), _ALLOCATABLE, dimension(:) :: otf      _NULL ! flux through open boundaries, allocate nobc

     real                              :: conversion            ! conversion between  dimensions  
     real                              :: offset                ! offset in dimensions (e.g., Celsius to Kelvin)
     real                              :: min_tracer            ! min acceptable value used for error checking 
     real                              :: max_tracer            ! max acceptable value used for error checking 
     real                              :: min_range             ! min value used for calls to diagnostic manager
     real                              :: max_range             ! max value used for calls to diagnostic manager
     real                              :: min_tracer_limit      ! min value used to limit quicker & neutral fluxes
     real                              :: max_tracer_limit      ! max value used to limit quicker & neutral fluxes 
     real                              :: min_flux_range        ! min and max values used for flux diagnostics
     real                              :: max_flux_range        ! min and max values used for flux diagnostics
     real                              :: const_init_value      ! value used to initialize constant tracer
     logical                           :: const_init_tracer     ! false (default) if the tracer must exist in the restart file
                                                                !   otherwise will initialize with const_init_value
     character(len=128)                :: flux_units            ! units for the tracer flux
     character(len=128)                :: restart_file          ! name for restart file
  end type ocean_prog_tracer_type

  type, public :: ocean_diag_tracer_type
     character(len=32)  :: name
     character(len=128) :: units, type
     character(len=128) :: longname
     real, dimension(:,:,:), pointer  :: field       ! tracer concentration at single time level 
     real :: conversion                              ! conversion between dimensions  
     real :: offset                                  ! offset in dimensions (e.g., Celsius to Kelvin)
     real :: min_tracer, max_tracer                  ! min and max acceptable values used for error checking 
     real :: min_range, max_range                    ! min and max values used for diagnostics
     logical :: const_init_tracer                    ! false (default) if the tracer must exist in the restart file
                                                     !   otherwise will initialize with const_init_value
     character(len=128) :: restart_file              ! name for restart file
     real               :: const_init_value          ! value used to initialize constant tracer
  end type ocean_diag_tracer_type

  type, public :: ocean_velocity_type
     logical                                  :: bmf_implicit      ! is true when time stepping bmf implicitly
     real, _ALLOCATABLE, dimension(:,:,:,:,:) :: u               _NULL ! horz velocity (m/s) in i,j directions at 3 time levels
     real, _ALLOCATABLE, dimension(:,:,:)     :: smf             _NULL ! momentum flux into ocn surface (N/m^2) at uv point
     real, _ALLOCATABLE, dimension(:,:,:)     :: smf_bgrid       _NULL ! momentum flux into ocn surface (N/m^2) at Bgrid uv point
     real, _ALLOCATABLE, dimension(:,:,:)     :: smf_cgrid       _NULL ! momentum flux into ocn surface (N/m^2) at Cgrid u/v points
     real, _ALLOCATABLE, dimension(:,:,:)     :: bmf             _NULL ! momentum flux per mass into ocean bottom  (N/m^2)
     real, _ALLOCATABLE, dimension(:,:)       :: ustar           _NULL ! surface friction velocity (m/s)
     real, _ALLOCATABLE, dimension(:,:)       :: gamma           _NULL ! dimensionful bottom drag coefficient (kg/(m^2 sec))
     real, _ALLOCATABLE, dimension(:,:)       :: langmuirfactor  _NULL ! Langmuir turbulence enhancement factor
     real, _ALLOCATABLE, dimension(:,:)       :: u10             _NULL ! 10m wind speed (m/s)
     real, _ALLOCATABLE, dimension(:,:)       :: ustoke          _NULL ! x-dir surface stokes drift
     real, _ALLOCATABLE, dimension(:,:)       :: vstoke          _NULL ! y-dir surface stokes drift
     real, _ALLOCATABLE, dimension(:,:)       :: wavlen          _NULL ! wave length
     real, _ALLOCATABLE, dimension(:,:)       :: cdbot_array     _NULL ! dimensionless static bottom drag coefficient
     real, _ALLOCATABLE, dimension(:,:)       :: current_wave_stress _NULL !wave-current bottom stress for sediment dynamics (N/m^2)
     real, _ALLOCATABLE, dimension(:,:)       :: rossby_radius   _NULL ! first baroclinic rossby radius (m)
     real, _ALLOCATABLE, dimension(:,:)       :: stokes_depth    _NULL ! depth scale (m) for exp decay of surface wave Stokes vel
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: stokes_drift    _NULL ! Stokes drift velocity (m/s) from surface wave model 
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: stokes_force    _NULL ! Coriolis force from Stokes drift velocity (N/m2)
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: press_force     _NULL ! rho*dz*horz (i,j)-directed press force (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: accel           _NULL ! rho*dz*velocity (i,j)-directed acceleration (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: vfrict_impl     _NULL ! rho*dz*vertical friction (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: source          _NULL ! thickness*density weighted velocity source/sink (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:,:)   :: wrkv            _NULL ! work array 
     real, _ALLOCATABLE, dimension(:,:,:,:,:) :: advection       _NULL ! rho*dz*velocity advection tendency (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:,:,:) :: coriolis        _NULL ! rho*dz*velocity Cgrid coriolis tendency (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:)     :: lap_friction_bt _NULL ! friction just on barotropic (N/m^2)
     real, _ALLOCATABLE, dimension(:,:,:)     :: bih_friction_bt _NULL ! friction just on barotropic (N/m^2)
  end type ocean_velocity_type


  type, public :: ocean_external_mode_type
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_t        _NULL ! surface height on tracer cell center (m) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_u        _NULL ! surface height on tracer cell center (m)
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_t_bar    _NULL ! eta_t time averaged over ext-mode time steps (m)
     real, _ALLOCATABLE, dimension(:,:)     :: deta_dt      _NULL ! surface height time tendency on t-cell (m/s)

     real, _ALLOCATABLE, dimension(:,:,:)   :: pbot_t     _NULL ! bottom pressure on tracer cell center (Pa) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: pbot_u     _NULL ! bottom pressure on velocity cell center (Pa) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: anompb     _NULL ! pbot_t_rho0*grav*ht (Pa) = anomalous bottom pressure (m^2/s^2)
     real, _ALLOCATABLE, dimension(:,:,:)   :: anompb_bar _NULL ! pbot_t-rho0*grav*ht (Pa) time avg over ext-mode time steps 
     real, _ALLOCATABLE, dimension(:,:)     :: dpbot_dt   _NULL ! bottom pressure time tendency on tracer cell (Pa/s)

     real, _ALLOCATABLE, dimension(:,:,:)   :: patm_t           _NULL ! ice plus atmospheric press on tracer cell (Pa)
     real, _ALLOCATABLE, dimension(:,:)     :: patm_for_sea_lev _NULL ! ice plus atmospheric press on tracer cell, to coupler (Pa)
     real, _ALLOCATABLE, dimension(:,:)     :: patm_u           _NULL ! ice plus atmospheric press on velocity cell (Pa)
     real, _ALLOCATABLE, dimension(:,:)     :: dpatm_dt         _NULL ! time tendency of patm on T-cell(Pa/s)

     real, _ALLOCATABLE, dimension(:,:,:)   :: forcing_bt _NULL ! depth integrated forcing of barotropic (w/o coriolis)
     real, _ALLOCATABLE, dimension(:,:,:,:) :: udrho      _NULL ! vertically integrated & rho wghted horz velocity (kg/m^3)*(m^2/s)
     real, _ALLOCATABLE, dimension(:,:,:)   :: conv_rho_ud_t  _NULL ! convergence of sum_k(rho*dzt*u) on T-cell (kg/m^3)*(m/s)
     real, _ALLOCATABLE, dimension(:,:)     :: ps             _NULL ! surface pressure (pressure at z=0 due to eta_t) at time tau
     real, _ALLOCATABLE, dimension(:,:,:)   :: grad_ps        _NULL ! rho0r * horizontal surface pressure gradient at time tau
     real, _ALLOCATABLE, dimension(:,:,:)   :: grad_anompb    _NULL ! gradient of anompb (Pa/m)
     real, _ALLOCATABLE, dimension(:,:,:)   :: press_force    _NULL ! pressure force/area (N/m^2) updated each barotropic time step  

     real, _ALLOCATABLE, dimension(:,:)     :: conv_blob      _NULL ! water column divergence of the L system

     real, _ALLOCATABLE, dimension(:,:)     :: source         _NULL ! density weighted vertical integral of mass source (kg/m^2/sec)
     real, _ALLOCATABLE, dimension(:,:)     :: eta_smooth     _NULL ! (kg/m^3)*(m/s) from eta_t smoother 
     real, _ALLOCATABLE, dimension(:,:)     :: pbot_smooth    _NULL ! (kg/m^3)*(m/s) from pbot_t smoother 

     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_nonbouss   _NULL ! diagnosed eta (m) including steric effects
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_nonsteric  _NULL ! diagnosed piece of eta_nonbouss (m) from non-steric effects 
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_steric     _NULL ! diagnosed eta_nonbouss (m) from steric effects 
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_dynamic    _NULL ! diagnosed eta_nonbouss (m) from "dynamic" effects 
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_water      _NULL ! diagnosed eta_nonbouss (m) from water forcing at boundaries 
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_source     _NULL ! diagnosed eta_nonbouss (m) from source term  
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_surf_temp  _NULL ! diagnosed eta_nonbouss (m) arising from surface temp flux
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_surf_salt  _NULL ! diagnosed eta_nonbouss (m) arising from surface salt flux
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_surf_water _NULL ! diagnosed eta_nonbouss (m) arising from surface water flux
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_bott_temp  _NULL ! diagnosed eta_nonbouss (m) arising from bottom temp flux
  end type ocean_external_mode_type


  ! types used by neutral physics 
  type, public :: tracer_2d_type
    real, _ALLOCATABLE, dimension(:,:)   :: field  _NULL
  end type tracer_2d_type

  type, public :: tracer_3d_0_nk_type
    real, _ALLOCATABLE, dimension(:,:,:) :: field  _NULL
  end type tracer_3d_0_nk_type

  type, public :: tracer_3d_1_nk_type
    real, _ALLOCATABLE, dimension(:,:,:) :: field  _NULL
  end type tracer_3d_1_nk_type


  ! for gotm vertical mixing scheme 
  type, public :: ocean_gotm_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname
     real               :: min_value ! minimum value
     real               :: max_value ! maximum value
     real, _ALLOCATABLE, dimension(:,:,:,:) :: field  _NULL ! scalar at 2 time levels
  end type ocean_gotm_type

  ! advection tendency for GOTM scalar fields 
  type, public  :: advect_gotm_type  
    real, dimension(:,:,:), pointer :: field => NULL() 
  end type advect_gotm_type


#endif
!############################################################################################
! end of STATIC_MEMORY

  type, public :: ocean_lagrangian_type
     real, _ALLOCATABLE, dimension(:,:,:) :: rho_dztlo _NULL ! L system denisty*volume/dat
     real, _ALLOCATABLE, dimension(:,:,:) :: rho_dztup _NULL ! L system denisty*volume/dat
     real, _ALLOCATABLE, dimension(:,:,:) :: conv_blob _NULL ! the convergence of blobs in a water column; 
                                                             ! used for surface height and bottom pressure
  end type ocean_lagrangian_type


  type, public :: ice_ocean_boundary_type
     real, pointer, dimension(:,:) :: u_flux          =>NULL() ! i-directed wind stress into ocean (Pa) 
     real, pointer, dimension(:,:) :: v_flux          =>NULL() ! j-directed wind stress into ocean (Pa) 
     real, pointer, dimension(:,:) :: t_flux          =>NULL() ! sensible heat flux into ocean (W/m2) 
     real, pointer, dimension(:,:) :: q_flux          =>NULL() ! specific humidity flux (kg/m2/s)
     real, pointer, dimension(:,:) :: salt_flux       =>NULL() ! salt flux (kg/m2/s)
     real, pointer, dimension(:,:) :: lw_flux         =>NULL() ! long wave radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_vis_dir =>NULL() ! direct visible sw radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_vis_dif =>NULL() ! diffuse visible sw radiation (W/m2) 
     real, pointer, dimension(:,:) :: sw_flux_nir_dir =>NULL() ! direct near IR sw radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_nir_dif =>NULL() ! diffuse near IR sw radiation (W/m2) 
     real, pointer, dimension(:,:) :: lprec           =>NULL() ! mass flux of liquid precip (kg/m2/s)
     real, pointer, dimension(:,:) :: fprec           =>NULL() ! mass flux of frozen precip (kg/m2/s)
     real, pointer, dimension(:,:) :: runoff          =>NULL() ! mass flux of liquid runoff (kg/m2/s) 
     real, pointer, dimension(:,:) :: calving         =>NULL() ! mass flux of frozen runoff (kg/m2/s) 
     real, pointer, dimension(:,:) :: runoff_hflx     =>NULL() ! heat flux, relative to 0C, of liquid land water into ocean (W/m2) 
     real, pointer, dimension(:,:) :: calving_hflx    =>NULL() ! heat flux, relative to 0C, of frozen land water into ocean (W/m2) 
     real, pointer, dimension(:,:) :: p               =>NULL() ! pressure of overlying sea ice and atmosphere (Pa)
     real, pointer, dimension(:,:) :: mi              =>NULL() ! mass of overlying sea ice 
     real, pointer, dimension(:,:) :: langmuirfactor  =>NULL() ! langmuir turbulence boost factor (non-dimensional)
     real, pointer, dimension(:,:) :: ustoke          =>NULL() ! x-dir surface stokes drift
     real, pointer, dimension(:,:) :: vstoke          =>NULL() ! y-dir surface stokes drift
     real, pointer, dimension(:,:) :: wavlen          =>NULL() ! wave length
#if defined(ACCESS_CM) || defined(ACCESS_OM)
     real, pointer, dimension(:,:) :: aice             =>NULL() !  ice fraction
     real, pointer, dimension(:,:) :: mh_flux          =>NULL() ! heat flux from melting ice (W/m^2)
     real, pointer, dimension(:,:) :: wfimelt          =>NULL() ! water flux from melting ice (kg/m^2/s)
     real, pointer, dimension(:,:) :: wfiform          =>NULL() ! water flux from forming ice (kg/m^2/s)
     real, pointer, dimension(:,:) :: licefw           =>null() ! waterflux into ocean (kg/m2/s) off Antarctica and Greenland
     real, pointer, dimension(:,:) :: liceht           =>null() ! heatflux due to land ice melt (W/m2)
#endif
#if defined(ACCESS_CM)
     real, pointer, dimension(:,:) :: co2              =>NULL() ! co2
#endif
     real, pointer, dimension(:,:) :: wnd              =>NULL() ! wind speed
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
     real, pointer, dimension(:,:) :: iof_nit              =>NULL() ! ice-ocean flux of nitrate
     real, pointer, dimension(:,:) :: iof_alg              =>NULL() ! ice-ocean flux of algae
#endif
     integer :: xtype                                          ! REGRID, REDIST or DIRECT

     type(coupler_2d_bc_type)      :: fluxes                   ! array of fields used for additional tracers
  end type ice_ocean_boundary_type


  ! for communication with FMS coupler and ESMF coupler 
  type, public ::  ocean_public_type 
     type(domain2d) :: Domain

     ! ESMF requires these arrays to be pointers 
     real, pointer, dimension(:,:)    :: t_surf  =>NULL() ! SST on t-cell (degrees Kelvin)
     real, pointer, dimension(:,:)    :: s_surf  =>NULL() ! SSS on t-cell (psu)
     real, pointer, dimension(:,:)    :: u_surf  =>NULL() ! i-directed surface ocean velocity on u-cell (m/s)
     real, pointer, dimension(:,:)    :: v_surf  =>NULL() ! j-directed surface ocean velocity on u-cell (m/s)
     real, pointer, dimension(:,:)    :: sea_lev =>NULL() ! eta_t + patm/(rho0*grav) - eta_geoid - eta_tide (m) 
     real, pointer, dimension(:,:)    :: frazil  =>NULL() ! accumulated heating (J/m^2) from 
                                                          ! frazil formation in the ocean 
     real, pointer, dimension(:,:)    :: area    =>NULL() ! T-cell area.
#if defined(ACCESS_CM) || defined(ACCESS_OM)
     real, pointer, dimension(:,:,:)  :: gradient =>NULL() ! x/y slopes of sea surface.
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
     real, pointer, dimension(:,:)    :: n_surf =>NULL() ! sea surface nitrate (mmol m-3)
     real, pointer, dimension(:,:)    :: alg_surf =>NULL() ! sea surface algae (mmol m-3)
#endif
#if defined(ACCESS_CM)
     real, pointer, dimension(:,:)    :: co2     =>NULL() ! co2 ( )
     real, pointer, dimension(:,:)    :: co2flux =>NULL() ! co2 flux ()
#endif
     logical, pointer, dimension(:,:) :: maskmap =>NULL()! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used.
     integer, pointer, dimension(:) :: pelist  =>NULL()  ! Used for flux-exchange.
     integer                        :: avg_kount         ! Used for accumulating averages - can be omitted?
     logical                        :: is_ocean_pe       ! .true. on processors that run the ocean model.
     integer, dimension(3)          :: axes              ! for diagnostics     
     type(coupler_2d_bc_type)       :: fields            ! array of fields used for additional tracers
  end type ocean_public_type

type, public :: ocean_blob_type
   ! Note that if an allocatable array is added, free_blob_memory in the ocean_blob_util module 
   ! must be edited to ensure complete deallocation of blobs occurs
   integer :: i, j, k         ! blobs current tracer grid box
   integer :: m, kdw          ! indices for the overflow scheme
   integer :: kup             ! indices for the overflow scheme
   integer :: hash, number    ! for sorting blobs
   integer :: model_steps     ! number of Eulerian model steps
   integer :: nsteps          ! total number of steps
   logical :: sink            ! whether the blob sinks or rises
   logical :: new             ! if the blob is newly created
   real    :: h1, h2          ! stretching coefficients
   real    :: lat, lon        ! current horizontal location of blob
   real    :: depth, geodepth ! current depth of blob
   real    :: st
   real    :: mass            ! current mass of blob
   real    :: density         ! current density of a blob
   real    :: densityr        ! 1./density
   real    :: volume          ! current volume of a blob
   real    :: height          ! height of the blob (m)
   real    :: step            ! time step size
   real    :: blob_time       ! blob time relative to last E system time step
   real    :: ent             ! entrainment velocity (m/s)
   real    :: det             ! detrainment velocity (m/s)
   real    :: richardson      ! Richardson number
   real    :: drag            ! Rayleigh drag coefficient (1/s)
   real    :: gprime          ! The reduced gravity (for diagnostics; m/s^2)
   real    :: age             ! Age of the blob (for diagnostics; seconds)
   integer :: nfrac_steps     ! number of fractional steps taken in one E step
   type(ocean_blob_type), pointer :: next ! next blob in the list
   type(ocean_blob_type), pointer :: prev ! previous blob in the list
   real,    _ALLOCATABLE, dimension(:)   :: tracer      ! current tracer content of blob
   real,    _ALLOCATABLE, dimension(:)   :: field       ! current tracer concentration of blob
   real,                  dimension(3)   :: v           ! blob's velocity (u,v,w)
   integer, _ALLOCATABLE, dimension(:)   :: di, dj, dk  ! grid cells where properties need to be communicated to E system
   real                                  :: dmass       ! change in mass
   real,    _ALLOCATABLE, dimension(:)   :: dtracer     ! change in tracer
   real,    _ALLOCATABLE, dimension(:,:) :: entrainment ! changes due to entrainment
   real,    _ALLOCATABLE, dimension(:,:) :: detrainment ! changes due to detrainment
   real,    _ALLOCATABLE, dimension(:)   :: mass_in(:)  ! mass into a grid cell
   real,    _ALLOCATABLE, dimension(:)   :: mass_out(:) ! mass out of a grid cell
end type ocean_blob_type

type, public :: blob_grid_type
   integer :: nk_nj                           ! nk*nj
   integer :: pe_this, pe_E, pe_W, pe_N, pe_S ! MPI PE identifiers
   integer :: pe_NE, pe_NW, pe_SE, pe_SW      ! MPI PE identifiers
   real,    _ALLOCATABLE, dimension(:,:,:) :: ht
   real,    _ALLOCATABLE, dimension(:,:,:) :: hu
   integer, _ALLOCATABLE, dimension(:,:,:) :: uidx
   integer, _ALLOCATABLE, dimension(:,:,:) :: tidx
   integer, _ALLOCATABLE, dimension(:)     :: it
   integer, _ALLOCATABLE, dimension(:)     :: jt
   integer, _ALLOCATABLE, dimension(:)     :: iu
   integer, _ALLOCATABLE, dimension(:)     :: ju
   real,    _ALLOCATABLE, dimension(:)     :: minlon
   real,    _ALLOCATABLE, dimension(:)     :: maxlon
   real,    _ALLOCATABLE, dimension(:)     :: minlat
   real,    _ALLOCATABLE, dimension(:)     :: maxlat
end type blob_grid_type

type, public :: blob_diag_type
   real, _ALLOCATABLE, dimension(:,:,:) :: entrainment
   real, _ALLOCATABLE, dimension(:,:,:) :: detrainment
   real, _ALLOCATABLE, dimension(:,:,:) :: new
   real, _ALLOCATABLE, dimension(:,:,:) :: dstry
end type blob_diag_type



public ocean_types_init
 
contains 


!#######################################################################
! <SUBROUTINE NAME="ocean_types_init">
!
! <DESCRIPTION>
! Initialize the ocean types. 
! </DESCRIPTION>
!
  subroutine ocean_types_init()

    if (module_is_initialized) then 
       call mpp_error( FATAL, '==>Error: ocean_types_init: module already initialized')
    endif
    module_is_initialized = .true.

    call write_version_number(version, tagname)

    return

  end subroutine ocean_types_init
! </SUBROUTINE> NAME="ocean_types_init"


end module ocean_types_mod


