module land_properties_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Christopher Milly
! </CONTACT> 

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Elena Shevliakova
! </REVIEWER>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Sergey Malyshev
! </REVIEWER>

! <OVERVIEW>
!   Contains land properties namelist variables and procedures relating
!   to the properties of the land.
! </OVERVIEW>

! <DESCRIPTION>
!     Initialization and calculation of the land property data. The input
!     cover field and implied glacier field are obtained and the glacier mask
!     and cover type are redefined, if necessary. The properties that depend
!     on cover and ground types are assigned. The land properties diagnostics
!     are initialized and the static fields are sent to the diagnostic manager
!     for output.
!
!     Updates the rapidly changing parameters. Computes the albedo of
!     the hypothetical no-snow and deep-snow surfaces and uses snow mass to
!     blend snow-free and deep-snow albedo values. Regrids integer index data
!     to any output grid from a uniformly spaced grid using a maximum-area
!     voting scheme. Includes calculation of total area of each land cover
!     type within the territory and the determination of the type occupying
!     the biggest area within the territory.
! </DESCRIPTION>

use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain, mpp_get_global_domain
use time_manager_mod,   only: time_type, set_date, print_date, get_date,     &
                              operator(-), operator(>)
use mpp_io_mod,         only: mpp_open, mpp_RDONLY
use fms_mod,            only: error_mesg, file_exist, open_namelist_file,    &
                              check_nml_error, stdlog, write_version_number, &
                              mpp_pe, mpp_root_pe, close_file, read_data,    &
                              FATAL, NOTE, open_restart_file, set_domain,    &
                              field_size, stdout
use constants_mod,      only: tfreeze, pi
use horiz_interp_mod,   only: horiz_interp_type, horiz_interp_new, horiz_interp, &
                              horiz_interp_del
use numerics_mod,       only: is_latlon, expand_cell
use climap_albedo_mod,  only: get_climap_albedo,     get_climap_glacier,     &
                              get_climap_albedo_mcm, get_climap_glacier_mcm
use diag_manager_mod,  only : register_static_field, send_data, get_base_time
use topography_mod,     only: get_topog_stdev
implicit none
private

! ==== public interface =====================================================
public :: land_properties_init,land_properties_end
public :: update_land_properties_slow
public :: update_land_properties_fast

public :: regrid_discrete_field
! ==== end public interface =================================================

! <NAMELIST NAME="land_properties_nml">
! n_dim_ground_types        size of ground (soil) parameter lookup tables
! n_dim_cover_types         size of cover (vegetation) parameter lookup tables
! n_map_ground_types        number of ground types in input map file
! n_map_cover_types         number of cover types in input map file

!   <DATA NAME="n_dim_ground_types" TYPE="integer, parameter" DEFAULT="11">
!     Size of ground (soil) parameter lookup tables
!   </DATA>

!   <DATA NAME="n_dim_cover_types" TYPE="integer, parameter" DEFAULT="14">
!     Size of cover (vegetation) parameter lookup tables
!   </DATA>

!   <DATA NAME="n_map_ground_types" TYPE="integer, parameter" DEFAULT="10">
!     Number of ground types in input map file
!   </DATA>

!   <DATA NAME="n_map_cover_types" TYPE="integer, parameter" DEFAULT="11">
!     Number of cover types in input map file
!   </DATA>

integer, parameter   :: n_dim_ground_types = 11  ! = 10+ssl
integer, parameter   :: n_dim_cover_types  = 14  ! = 11+tun+ssl_veg+ssl_ice
integer, parameter   :: n_map_ground_types = 10
integer, parameter   :: n_map_cover_types  = 11

!---- namelist - with default values ----------------------------------
!
! do_all_mcm       run Manabe Climate Model land surface       [logical]
!          Setting this global control option to TRUE causes
!          specification all the following (regardless of default or input
!          settings, which are then ignored):
!              veg_to_use      = 'cons_ssl', 
!              soil_to_use     = 'cons_ssl', 
!              use_climap      = .true.  (this forced by veg_to_use='cons_ssl'), 
!              use_climap_mcm  = .true., 
!              do_mcm_masking  = .true., 
!              use_single_geo  = .true., 
!              geo_res_time    = res_time_ssl, 
!              factor_root     = 1., 
!              factor_rough    = 1., 
!              z_root_min      = 0., 
!              max_snow        = max_snow_ssl
! veg_to_use          choice of method for defining vegetation   [character]
!          'variable' - use input map to index vegetation parameter vectors
!          'constant' - use veg_index_constant to index veg parameter vectors
!          'cons_tun' - use tuned global constant vegetation
!          'cons_ssl' - use Manabe Climate Model-like vegetation (forces use_climap to be true)
! soil_to_use         choice of method for defining soil         [character]
!          'variable' - use input map to index soil parameter vectors
!          'constant' - use soil_index_constant to index soil parameter vectors
!          'cons_ssl' - use Manabe Climate Model-like soil
! use_glaciers        false to remove glaciers from land          [logical]
! use_climap          true to override default albedo by CLIMAP   [logical]
!                     and to invoke use of CLIMAP to define glacier
!                     locations (if use_glaciers)
! use_desert_albedo_map   true to override default snow-free albedo of desert only by
!                         albedo map (SRB)
! do_mcm_masking      true to use Manabe Climate Model snow-albedo fn. [logical]
! use_single_basin    true to avoid using basin maps              [logical]
! use_single_geo      true for global constant gw residence time  [logical]
! i_dest0             if use_single_geo is set to .true., all the river
!                     discharge is put in a single grid cell. i_dest0 is the
!                     longitude index of this grid cell.  
! j_dest0             if use_single_geo is set to .true., all the river
!                     discharge is put in a single grid cell. j_dest0 is the
!                     latitude index of this grid cell.            
! veg_index_constant  veg index used when veg_to_use='constant'          -
! soil_index_constant soil index used when soil_to_use='constant'        -
! soil_index_ice_substitute   ground type to be substituted for
!                     ice ground type when such type is over-ruled due to
!                     absence of glacier cover type. if this is set to 
!                     zero, then the ground type is temporarily marked as
!                     ocean, and, if land is present, ground type will then
!                     be assigned based on neighbor cells through regrid_discrete_field.
! geo_res_time        time constant when use_single_geo is true          s
! t_range             T range over which snow/glacier albedo varies      K
! factor_root         global factor for critical_root_density            -
! factor_stomata      global factor for veg_rs_min_vec                   -
! factor_rough        global factor for rough_momentum_vec               -
! z_root_min          lower bound for root-zone depth                    m
! max_snow            value of snow above which 'snow runoff' occurs  kg/m**3
! sfc_heat_factor     "fudge" factor for heat capacity and thermal
!                     conductivity in surface soil layers                -
! num_sfc_layers      number of surface layers for sfc_heat_factor       -
! dynamic_cover_type    set to true if cover type forcing data varies with
!                       time                               [logical]
! read_old_ascii_cover  read the original ASCII file of cover type
!                                                          [logical]
! cover_dataset_init_year  initial year in the cover_type dataset [integer]
! cover_dataset_entry      beginning time for reading the cover_type
!                          dataset (yr, mo, dy, hr, mn, sc)       [integer]


!   <DATA NAME="do_all_mcm" TYPE="logical" DEFAULT=".false.">
!     Run Manabe Climate Model land surface. Setting this global control
!     option to TRUE causes specification all the following (regardless of
!     default or input settings, which are then ignored):
!              veg_to_use      = 'cons_ssl'
!              soil_to_use     = 'cons_ssl'
!              use_climap      = .true.  (this forced by veg_to_use='cons_ssl')
!              use_climap_mcm  = .true.
!              do_mcm_masking  = .true.
!              use_single_geo  = .true.
!              geo_res_time    = res_time_ssl
!              factor_root     = 1.
!              factor_rough    = 1.
!              z_root_min      = 0.
!              max_snow        = max_snow_ssl
!   </DATA>

!   <DATA NAME="veg_to_use" TYPE="character*8" DEFAULT="'variable'">
!     Choice of method for defining vegetation.
!          'variable' - use input map to index vegetation parameter vectors
!          'constant' - use veg_index_constant to index veg parameter vectors
!          'cons_tun' - use tuned global constant vegetation
!          'cons_ssl' - use Manabe Climate Model-like vegetation (forces use_climap to be true)
!   </DATA>

!   <DATA NAME="soil_to_use" TYPE="character*8" DEFAULT="'variable'">
!     Choice of method for defining soil.
!          'variable' - use input map to index soil parameter vectors
!          'constant' - use soil_index_constant to index soil parameter vectors
!          'cons_ssl' - use Manabe Climate Model-like soil
!   </DATA>

!   <DATA NAME="use_glaciers" TYPE="logical" DEFAULT=".true.">
!     False to remove glaciers from land
!   </DATA>

!   <DATA NAME="use_climap" TYPE="logical" DEFAULT=".false.">
!     True to override default albedo by CLIMAP and to invoke use of CLIMAP
!     to define glacier locations (if use_glaciers)
!   </DATA>

!   <DATA NAME="use_desert_albedo_map" TYPE="logical" DEFAULT=".false.">
!     true to override default snow-free albedo of desert only by
!     albedo map (SRB)
!   </DATA>

!   <DATA NAME="use_climap_mcm" TYPE="logical" DEFAULT=".false.">
!     Run Manabe Climate Model land surface
!   </DATA>

!   <DATA NAME="do_mcm_masking" TYPE="logical" DEFAULT=".false.">
!     True to use Manabe Climate Model snow-albedo function
!   </DATA>

!   <DATA NAME="use_single_basin" TYPE="logical" DEFAULT=".false.">
!     True to avoid using basin maps
!   </DATA>

!   <DATA NAME="use_single_geo" TYPE="logical" DEFAULT=".false.">
!     True for global constant groundwater residence time
!   </DATA>

!   <DATA NAME="i_dest0" TYPE="integer" DEFAULT="1">
!     If use_single_geo is set to .true., all the river discharge is put in a
!     single grid cell. i_dest0 is the longitude index of this grid cell. 
!   </DATA>

!   <DATA NAME="j_dest0" TYPE="integer" DEFAULT="1">
!     If use_single_geo is set to .true., all the river discharge is put in a
!     single grid cell. j_dest0 is the latitude index of this grid cell. 
!   </DATA>

!   <DATA NAME="veg_index_constant" TYPE="integer" DEFAULT="3">
!     Veg index used when veg_to_use='constant'
!   </DATA>

!   <DATA NAME="soil_index_constant" TYPE="integer" DEFAULT="2">
!     Soil index used when soil_to_use='constant'
!   </DATA>

!   <DATA NAME="soil_index_ice_substitute" TYPE="integer" DEFAULT="1">
!     Ground type to be substituted for ice ground type when such type is
!     over-ruled due to absence of glacier cover type. if this is set to zero,
!     then the ground type is temporarily marked as ocean, and, if land is
!     present, ground type will then be assigned based on neighbor cells
!     through regrid_discrete_field.
!   </DATA>

!   <DATA NAME="geo_res_time" UNITS="s" TYPE="real" DEFAULT="60.*86400.">
!     Time constant when use_single_geo is true
!   </DATA>

!   <DATA NAME="t_range" UNITS="K" TYPE="real" DEFAULT="10.0">
!     Temperature range over which snow/glacier albedo varies
!   </DATA>

!   <DATA NAME="factor_root" TYPE="real" DEFAULT="1.0">
!     Global factor for critical_root_density
!   </DATA>

!   <DATA NAME="factor_stomata" TYPE="real" DEFAULT="1.0">
!     Global factor for veg_rs_min_vec
!   </DATA>

!   <DATA NAME="factor_rough" TYPE="real" DEFAULT="1.0">
!     Global factor for rough_momentum_vec
!   </DATA>

!   <DATA NAME="z_root_min" UNITS="m" TYPE="real" DEFAULT="0.01">
!     lower bound for root-zone depth
!   </DATA>

!   <DATA NAME="max_snow" UNITS="kg/m3" TYPE="real" DEFAULT="1000.">
!     Value of snow above which 'snow runoff' occurs
!   </DATA>

!   <DATA NAME="sfc_heat_factor" TYPE="real" DEFAULT="1.">
!     "fudge" factor for heat capacity and thermal
!      conductivity in surface soil layers
!   </DATA>

!   <DATA NAME="num_sfc_layers" TYPE="integer" DEFAULT="0">
!     number of surface layers for sfc_heat_factor
!   </DATA>

!   <DATA NAME="dynamic_cover_type" TYPE="logical" DEFAULT="true">
!     Set to true if cover type forcing data varies with time.
!   </DATA>

!   <DATA NAME="read_old_ascii_cover" TYPE="logical" DEFAULT="false">
!     Set to true if reading the original ASCII static cover type forcing data
!   </DATA>

!   <DATA NAME="cover_dataset_init_year" TYPE="integer" DEFAULT="1860">
!     The initial year in the cover_type dataset.
!   </DATA>

!   <DATA NAME="cover_dataset_entry" TYPE="integer" DEFAULT="(/1,1,1,0,0,0/)">
!     Beginning time for reading the cover_type dataset (yr,mo,dy,hr,mn,sc).
!   </DATA>

logical                       :: do_all_mcm            = .false.
character*8                   :: veg_to_use            = 'variable'
character*8                   :: soil_to_use           = 'variable'
logical                       :: use_glaciers          = .true.
logical                       :: use_climap            = .false.
logical                       :: use_desert_albedo_map = .false.
logical                       :: use_climap_mcm        = .false.
logical                       :: do_mcm_masking        = .false.
logical                       :: use_single_basin      = .false.
logical                       :: use_single_geo        = .false.
integer                       :: i_dest0               = 1
integer                       :: j_dest0               = 1
integer                       :: veg_index_constant    = 3    ! mixed broadleaf/needleleaf
integer                       :: soil_index_constant   = 2    ! medium texture mineral soil
integer                   :: soil_index_ice_substitute = 1
logical :: reconcile_lakes   = .false. ! if true, lake distribution 
           ! in cover and ground type fields is reconciled based on the cover type
integer :: soil_index_lake_substitute = 1 ! substitute for the ground lake type,
           ! used only if reconcile_lakes is true
real                          :: geo_res_time          = 60.*86400.   ! two months
real                          :: t_range               = 10.0
real                          :: factor_root           =  1.0
real                          :: factor_stomata        =  1.0
real                          :: factor_rough          =  1.0
real                          :: z_root_min            =  0.01
real                          :: max_snow              = 1000.  ! (1 m water equiv)
real                          :: sfc_heat_factor       = 1.
integer                       :: num_sfc_layers        = 0
logical :: dynamic_cover_type = .true.  ! true if cover type forcing data
                                        ! varies with time. Must also set
                                        ! cover_dataset_entry.
logical :: read_old_ascii_cover = .false.       ! true if reading original
                                               ! static ACII cover type forcing.
integer, dimension(6) ::     &
       cover_dataset_entry  = (/ 1860, 1, 1, 0, 0, 0 /) 
          ! beginning time for reading the cover_type data. Must be set when
          ! using read_old_ascii_cover = false.

integer :: cover_dataset_init_year = 1860     ! initial year for the cover_type
                                              ! dataset.

!   <DATA NAME="veg_rs_min_vec" UNITS="s/m" TYPE="real" DIM="n_dim_cover_types">
!     Minimum bulk stomatal resistance. For default values, refer to the table
!     in the Public Code section below.
!   </DATA>
!   <DATA NAME="veg_zeta_vec" UNITS="m" TYPE="real" DIM="n_dim_cover_types">
!     Depth scale of root distribution. For default values, refer to the table
!     in the Public Code section below.
!   </DATA>
!   <DATA NAME="veg_root_mass_vec" UNITS="kg/m2" TYPE="real" DIM="n_dim_cover_types">
!     Root biomass areal density. For default values, refer to the table
!     in the Public Code section below.
!   </DATA>
!   <DATA NAME="rough_momentum_vec" UNITS="m" TYPE="real" DIM="n_dim_cover_types">
!     Roughness length for momentum. For default values, refer to the table
!     in the Public Code section below.
!   </DATA>
!   <DATA NAME="k_over_B_vec" TYPE="real" DIM="n_dim_cover_types">
!     ln (z_0_momentum / z_0_scalar). For default values, refer to the table
!     in the Public Code section below.
!   </DATA>
!   <DATA NAME="crit_snowmass_vec" UNITS="kg/m3" TYPE="real" DIM="n_dim_cover_types">
!     Snow amount that half hides surface. For default values, refer to the
!     table in the Public Code section below.
!   </DATA>
!   <DATA NAME="min_nosnow_alb_vec" TYPE="real" DIM="n_dim_cover_types">
!     Snow-free albedo at freezing point. For default values, refer to
!     the table in the Public Code section below. 
!   </DATA>
!   <DATA NAME="max_nosnow_alb_vec" TYPE="real" DIM="n_dim_cover_types">
!     Snow-free albedo at t_range below freezing. For default values, refer to
!     the table in the Public Code section below. 
!   </DATA>
!   <DATA NAME="min_snow_alb_vec" TYPE="real" DIM="n_dim_cover_types">
!     Snow albedo at freezing point.  For default values, refer to
!     the table in the Public Code section below. 
!   </DATA>
!   <DATA NAME="max_snow_alb_vec" TYPE="real" DIM="n_dim_cover_types">
!     Snow albedo at t_range below freezing.  For default values, refer to
!     the table in the Public Code section below. 
!   </DATA>
!   <DATA NAME="soil_awc_vec" UNITS="kg/m3" TYPE="real" DIM="n_dim_ground_types">
!     Available water capacity. For default values, refer to the table in
!     the Public Code section below.
!   </DATA>
!   <DATA NAME="soil_therm_cap_vec" UNITS="J/(K m3)" TYPE="real" DIM="n_dim_ground_types">
!     Volumetric heat capacity. For default values, refer to the table in the
!     Public Code section below.
!   </DATA>
!   <DATA NAME="soil_therm_dif_vec" UNITS="m2/s" TYPE="real" DIM="n_dim_ground_types">
!     Thermal diffusivity. For default values, refer to the table in the
!     Public Code section below.
!   </DATA>
!   <DATA NAME="use_topo_rough" TYPE="logical" DEFAULT="false">
!     If true, the topographic momentum drag scaling scheme is used
!   </DATA>
!   <DATA NAME="max_topo_rough" TYPE="real" DEFAULT="100" UNITS="m">
!     Maximum of topographic "roughness length" used for momentum drag scaling
!   </DATA>
!   <DATA NAME="topo_rough_factor" TYPE="real" DEFAULT="1.0">
!     Scaling factor to convert topography variance to topographic 
!     "roughness length"
!   </DATA>
!   <DATA NAME="topo_rough_source" TYPE="caharacter(len=16)" DEFAULT="'computed'">
!     Source of the sub-grid topography variance data for topographic momentum drag scaling. 
!     'computed' means that the variance is calculated based on high-resolution 
!     topography data. 'input' means that the data will be provided in specified file
!     (NetCDF of IEEE binary)
!   </DATA>
!   <DATA NAME="topo_rough_file" TYPE="character(len=256)" DEFAULT="INPUT/mg_drag.data.nc">
!     Name of the file to be used as an input for sub-grid topography variance data. 
!     The file can be either NetCDF (in this case variable name can also be specified), or
!     IEEE.
!   </DATA>
!   <DATA NAME="topo_rough_var" TYPE="character(len=128)" DEFAULT="ghprime">
!     Name of the NetCDF variable to be used as a topography variance field. Ignored if
!     the file specified in topo_rough_file is not NetCDF file.
!   </DATA>
! </NAMELIST>


!   <PUBLICCOMMENT>
! the following namelist vectors contain properties indexed by cover (veg) type:
!  0      ocean
!  1 (BE) broadleaf evergreen trees
!  2 (BD) broadleaf deciduous trees
!  3 (BN) broadleaf/needleleaf trees
!  4 (NE) needleleaf evergreen trees
!  5 (ND) needleleaf deciduous trees
!  6 (G)  grassland
!  7 (D)  desert
!  8 (T)  tundra
!  9 (A)  agriculture
! 10 (I)  ice
! 11 (L)  lake
! 12 (TV) tuned global vegetation
! 13 (SV) Manabe Climate Model-like vegetation
! 14 (SI) Manabe Climate Model-like ice/glacier

! veg_rs_min_vec       minimum bulk stomatal resistance           s/m
! veg_zeta_vec         depth scale of root distribution            m
! veg_root_mass_vec    root biomass areal density               kg/m**2
! rough_momentum_vec   roughness length for momentum               m
! k_over_B_vec         ln (z_0_momentum / z_0_scalar)              -
! crit_snowmass_vec    snow amount that half hides surface      kg/m**3
! min_nosnow_alb_vec   snow-free albedo at freezing point          -
! max_nosnow_alb_vec   snow-free albedo at t_range below freezing       -
! min_snow_alb_vec     snow albedo at freezing point               -
! max_snow_alb_vec     snow albedo at t_range below freezing       -

real, dimension(n_dim_cover_types) :: &
!                      BE     BD     BN     NE     ND      G      D      T      A      I      L     TV     SV     SI
veg_rs_min_vec    =(/  43.6,  131., 87.1,  69.7,  218.,  56.6,  .01,   170.,  56.6, .01,   .01,     67.,  .01,   .01/),& 
veg_zeta_vec      =(/  .26,   .29,  .35,   .17,   .17,   .26,   .35,   .11,   .25,   0.0,   1.0,   .35,   .35,   0.0/),&
veg_root_mass_vec =(/  4.9,   4.2,  4.3,   2.9,   2.9,   1.4,  .762,   1.2,   .15,   0.0,   1.0,  .362,  .762,   0.0/),&
rough_momentum_vec=(/ 2.65,  .90,   1.2,   .90,   .80,   .07,   .01,   .07,   .40,   .01, 1.4e-4,  1.0,  .045,  .045/),&
k_over_B_vec      =(/   2.,    2.,    2.,    2.,    2.,   2.,    2.,    2.,    2.,    2.,  0.25,    2.,    0.,    0. /),&
crit_snowmass_vec =(/  60.,   10.,   25.,   40.,   40.,   5.,    5.,    5.,    5.,    5.,    5.,  100.,    5.,    5. /),&
min_nosnow_alb_vec=(/0.149, 0.130, 0.132, 0.126, 0.143, 0.182, 0.333, 0.139, 0.160, 0.650, 0.06,  0.12,   999.,  0.55/),&
max_nosnow_alb_vec=(/0.149, 0.130, 0.132, 0.126, 0.143, 0.182, 0.333, 0.139, 0.160, 0.800, 0.06,  0.12,   999.,  0.65/),&
min_snow_alb_vec  =(/0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.06,  0.450, 0.450, 0.650/),&
max_snow_alb_vec  =(/0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.06,  0.600, 0.600, 0.800/)
!                      BE     BD     BN     NE     ND      G      D      T      A      I      L     TV     SV     SI


! the following vectors contain properties indexed by ground (soil) type:
!  0       ocean
!  1 (C)   coarse soil
!  2 (M)   medium soil
!  3 (F)   fine soil
!  4 (CM)  coarse/medium mix
!  5 (CF)  coarse/fine mix
!  6 (MF)  medium/fine mix
!  7 (CMF) coarse/medium/fine mix
!  8 (P)   organic soil (peat)
!  9 (I)   ice
! 10 (L)   lake
! 11 (MCM)  Manabe Climate Model

! soil_awc_vec         available water capacity                 kg/m**3
! soil_therm_cap_vec   volumetric heat capacity                J/(K m**3)
! soil_therm_dif_vec   thermal diffusivity                       m**2/s

real, dimension(n_dim_ground_types) :: &
!   C      M      F      CM     CF     MF    CMF     P      I      L      MCM
soil_awc_vec      = &
(/   63.,  132.,  109.,   98.,   86.,  120.,  101.,  445., 1000., 1000.,  150./),&
soil_therm_cap_vec= &
(/ 1.8e6, 2.0e6, 2.6e6, 1.9e6, 2.2e6, 2.3e6, 2.1e6, 3.0e6, 1.6e6, 8.4e7, 1.0/),&
soil_therm_dif_vec=&
(/8.3e-7,4.0e-7,5.2e-7,6.2e-7,6.8e-7,4.6e-7,5.8e-7,1.3e-7,1.1e-6, 1.0,  2.0e-7/)
!   C      M      F      CM     CF     MF    CMF     P      I      L      MCM
! </PUBLICCOMMENT>


logical     :: use_topo_rough    = .false.
real        :: max_topo_rough    = 100 ! m
real        :: topo_rough_factor = 1.0
character(len=16) :: topo_rough_source = 'computed'
character(len=256):: topo_rough_file   = 'INPUT/mg_drag.data.nc'
character(len=128):: topo_rough_var    = 'ghprime'

namelist /land_properties_nml/ do_all_mcm,         veg_to_use,         &
                               soil_to_use,        use_glaciers,       &
                               use_climap,         use_desert_albedo_map,  &
                               do_mcm_masking,     &
                               use_single_basin,   use_single_geo,     &
                               i_dest0,            j_dest0,            &
                               veg_index_constant, soil_index_constant,&
                               soil_index_ice_substitute,              &
                               geo_res_time,       t_range,            &
                               factor_root,        factor_stomata,     &
                               factor_rough,       z_root_min,         &
                               max_snow,                               &
                               veg_rs_min_vec,     veg_zeta_vec,       &
                               veg_root_mass_vec,                      &
                               rough_momentum_vec, k_over_B_vec,       &
                               crit_snowmass_vec,                      &
                               min_nosnow_alb_vec, max_nosnow_alb_vec, &
                               min_snow_alb_vec,   max_snow_alb_vec,   &
                               soil_awc_vec,                           &
                               soil_therm_cap_vec, soil_therm_dif_vec, &
                               use_topo_rough, max_topo_rough, topo_rough_factor, &
                               sfc_heat_factor, num_sfc_layers,        &
                               dynamic_cover_type, read_old_ascii_cover, &
                               cover_dataset_init_year, cover_dataset_entry, &
                               reconcile_lakes, soil_index_lake_substitute

! ---- private data ----------------------------------------------------------
logical :: module_is_initialized =.FALSE.

character(len=128) :: version = '$Id: land_properties.F90,v 19.0 2012/01/06 20:39:32 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!---- other module variables and named constants

! nlon, nlat             numbers of longitudes, latitudes in land grid
! npart                  number of partitions of each land grid cell
! ground_type            lon/lat/tile array of ground type index
! albedo_min_no_snow     lon/lat/tile array of min snow-free land albedo
! albedo_max_no_snow     lon/lat/tile array of max snow-free land albedo
! critical_root_density  root density at the bottom of the root zone (kg/m**3)
! max_snow_ssl           Manabe Climate Model value of max_snow              (kg/m**3)
! res_time_ssl           Manabe Climate Model value of groundwater res. time    (s)
! veg_index_desert       location of desert parameters in veg vecs
! veg_index_tundra       location of tundra parameters in veg vecs
! veg_index_ice          location of standard ice parameters in veg vecs
! veg_index_tuned        location of tuned parameters in veg vecs
! veg_index_ssl_veg      location of ssl veg parameters in veg vecs
! veg_index_ssl_ice      location of ssl ice parameters in veg vecs
! soil_index_coarse      location of coarse soil parameters in soil vecs
! soil_index_ice         location of ice parameters in soil vecs
! soil_index_ssl         location of ssl soil parameters in soil vecs
! topo_stdev             standard deviation of topography
! model_init_time        model base time - used to calculate cover type offset
! cover_entry            cover type dataset entry time 

integer                                :: nlon, nlat, npart
integer, pointer, dimension(:,:,:)     :: ground_type =>NULL()
real, pointer, dimension(:,:,:)        :: albedo_min_no_snow =>NULL()
real, pointer, dimension(:,:,:)        :: albedo_max_no_snow =>NULL()
real                          :: critical_root_density   =  0.125
real                          :: critical_glacier_albedo =  0.79
real                          :: max_snow_ssl            =  200.
real                          :: res_time_ssl            = 60.*86400.
integer                       :: veg_index_desert      =  7
integer                       :: veg_index_tundra      =  8
integer                       :: veg_index_ice         = 10
integer                       :: veg_index_lake        = 11
integer                       :: veg_index_tuned       = 12
integer                       :: veg_index_ssl_veg     = 13
integer                       :: veg_index_ssl_ice     = 14
integer                       :: soil_index_coarse     =  1
integer                       :: soil_index_ice        =  9
integer                       :: soil_index_lake       = 10
integer                       :: soil_index_ssl        = 11
real, allocatable :: topo_stdev(:,:)

type(time_type) :: model_init_time, cover_entry
type(horiz_interp_type), save :: Interp_gw

! ---- name of the diagnostic "module" ---------------------------------------
character(len=*), parameter :: diag_mod_name = "soil"
! ---- diagnostic fields ids -------------------------------------------------
integer :: id_alb_max, id_alb_min
integer :: id_ground_type
integer :: id_topo_stdev

include 'netcdf.inc'
contains

!#######################################################################

! <SUBROUTINE NAME="land_properties_init">

!   <OVERVIEW>
!     Initialize land property data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Reads and re-grids each land input data field and allocates arrays.
!     The input cover field and implied glacier field are obtained, if either
!     will be needed. Then, the ice cover type is reset to tundra. It
!     may be changed back to ice according to actual glacier specification
!     later. The glacier mask and cover type are redefined, if necessary.
!
!     The appropriate cover type to glacier cells is assigned, if any. The
!     ground (soil) type is defined and the groundwater residence time is
!     calculated. The properties that depend on cover and ground types are
!     assigned. If requested, the snow-free albedo is changed to the CLIMAP
!     array.
!
!     The land properties diagnostics are initialized and the static fields 
!     are sent to the diagnostic manager for output.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call land_properties_init (lonb,   latb,  land,  time,  domain, &
!     glacier, lake, rough_momentum, rough_heat, soil_therm_con,  &
!     soil_therm_cap, veg_rs_min, max_water, tau_groundwater, max_snow_out, &
!     id_lon, id_lat   )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine land_properties_init ( &
     lonb,   latb,  land,  time,  domain, &
     glacier,         &
     lake,            &
     rough_momentum,  &
     rough_heat,      &
     rough_scale,     &  ! topographic roughness for drag scaling
     soil_therm_con,  &
     soil_therm_cap,  &
     veg_rs_min,      &
     max_water,       &
     tau_groundwater, &
     max_snow_out,    &
     cover_type,      &
     id_lon, id_lat )

  real,            intent(in) :: lonb(:,:), latb(:,:) ! corners of the cells,
                                                  ! radian
  logical,         intent(in) :: land(:,:,:)      ! land mask
  type(time_type), intent(in) :: time             ! current time
  type(domain2d),  intent(in) :: domain           ! our domain 

  logical, intent(out), dimension(:,:,:) :: glacier, lake ! glacier, lake masks
  real   , intent(out), dimension(:,:,:) :: &
       rough_momentum, &                  ! roughness length for momentum
       rough_heat,     &                
       rough_scale,    &
       veg_rs_min,     & 
       max_water,      &
       tau_groundwater                    ! groundwater residence time
  real   , intent(out), dimension(:,:,:,:) :: &
       soil_therm_con, &
       soil_therm_cap
  real   , intent(out)  :: max_snow_out
  integer, intent(out), dimension(:,:,:) :: cover_type
  integer, intent(in) :: id_lon, id_lat   ! ids of diagnostic axes
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer, allocatable :: integer_input_1(:,:)
  real,    allocatable :: real_input     (:,:)
  real,    allocatable :: tmp(:,:,:)  ! temp array for read data for cover_type
  integer, allocatable :: cover_type_input (:,:,:)  ! cover_type input data
  integer, allocatable :: nlon_input (:)     ! input data lon
  integer, allocatable :: nlat_input (:)     ! input data lat
  real,    allocatable :: nblon_input (:)    ! input data western boundary

  real, dimension(size(land,1),size(land,2)):: stdev ! standard deviation of
                                                     ! orography 
  logical got_stdev

  character*80            input_format
  integer :: unit, io, ierr, i_cover_type, i_ground_type
  integer :: ii, jj, num_lon_input, num_lat_input, k
  real    :: wb_degrees ! western boundary of the input grid cells
  real    :: sb, wb, dx, dy, z_root
  integer :: is,ie, js,je ! boundaries of our domain
  integer :: gnlon, gnlat
  integer :: siz(4)       ! size of field in field_size call
  integer :: model_yr, model_mo, model_dy, model_hr, model_mn, model_sc
  integer :: init_yr, init_mo, init_dy, init_hr, init_mn, init_sc
  integer :: kk           ! time loop integer for cover type data
  integer :: outunit

  module_is_initialized = .TRUE.
  outunit = stdout()
  
  ! get the size of our domain
  call mpp_get_compute_domain ( domain, is,ie,js,je )

  call mpp_get_global_domain ( domain, xsize=gnlon, ysize=gnlat)

!---- read and write namelist and version -----------------
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_properties_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_properties_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  
  call write_version_number(version, tagname)
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_properties_nml)
     call close_file (unit)
  endif

  if (do_all_mcm) then
     veg_to_use        = 'cons_ssl'
     soil_to_use       = 'cons_ssl'
     do_mcm_masking    = .true.
     use_single_geo    = .true.
     geo_res_time      = res_time_ssl
     factor_root       = 1.
     factor_rough      = 1.
     z_root_min        = 0.
     max_snow          = max_snow_ssl
     use_climap_mcm    = .true.
  endif
  if (veg_to_use .eq. 'cons_ssl') then
     use_climap = .true.
  endif

   ! allocate private arrays
   nlon = size(lonb,1) - 1
   nlat = size(latb,2) - 1
   npart = size(land,3)
   allocate (ground_type   (nlon, nlat, npart))
   allocate (albedo_min_no_snow(nlon, nlat, npart))
   allocate (albedo_max_no_snow(nlon, nlat, npart))
   allocate (topo_stdev(nlon,nlat))

   ! Read and re-grid each land input data field

sb = -pi/2.

!  <NOTE>
!---- Cover (vegetation) type. There are 12 possible cases: veg_to_use can
! take 4 values, and glaciers can be (1) absent, (2) based on input cover
! field, or (3) based on climap albedo field. (use of climap albedo forces
! use of climap to locate glaciers, if glaciers are used.) the code
! ensures consistency between cover_type field and glacier mask, with
! the latter taking precedence: where glacier based on input cover field
! must be removed in deference to climap or if no glaciers are used,
! such points are assigned tundra (or global constant) type.
!  </NOTE>

!--- first get input cover field and implied glacier field, if either of
! these will be needed. then re-set ice cover type to tundra; it may be
! changed back to ice according to actual glacier specification later.
if ( veg_to_use.eq.'variable' &
     .or. (use_glaciers.and..not.use_climap) ) then

! check that dynamic_cover_type and read_old_ascii_cover both are not true
   if (dynamic_cover_type .and. read_old_ascii_cover) &
      call error_mesg('land_properties_init',                  &
       'The namelist variables dynamic_cover_type and read_old_ascii_cover cannot both have .true. values',FATAL)

! If dynamic_cover_type_netcdf is chosen, then the cover_type will be obtained
! from a netCDF file of annual data. The initial year of the data is 
! specified in the namelist option cover_dataset_init_year.
   if (dynamic_cover_type) then
     if(.not.file_exist('INPUT/cover_type_field.nc')) &
       call error_mesg('land_properties_init','Cannot find file INPUT/cover_type_field.nc',FATAL)
! <ERROR MSG="Cannot find file INPUT/cover_type_field.nc" STATUS="FATAL">
!   The cover type field file cannot be found. Provide this file or set up
!   namelist parameters so it is not necessary. To do the latter, set
!   veg_to_use to something other than 'variable' in the namelist
!   land_properties_nml.
! </ERROR>

     call error_mesg('land_properties_init', 'Using dynamic cover type forcing data',NOTE)

! obtain the size on lon in the input data
     call field_size('INPUT/cover_type_field.nc', 'lon', siz)
     allocate (nlon_input(siz(1)))
     num_lon_input = size(nlon_input)

! obtain the size of lat in the input data
     call field_size('INPUT/cover_type_field.nc', 'lat', siz)
     allocate (nlat_input(siz(1)))
     num_lat_input = size(nlat_input)

! allocate arrays, read data as real and convert back to integer
     allocate (tmp(num_lon_input, num_lat_input, npart)) ! real temp array
     allocate (cover_type_input (num_lon_input, num_lat_input, npart))  ! input data

! check to make sure cover_dataset_entry has been entered
     if (cover_dataset_entry(1) == 1 .and. &
         cover_dataset_entry(2) == 1 .and. &
         cover_dataset_entry(3) == 1 .and. &
         cover_dataset_entry(4) == 0 .and. &
         cover_dataset_entry(5) == 0 .and. &
         cover_dataset_entry(6) == 0 )     &
       call error_mesg ('land_properties_mod', &
             'must set cover_dataset_entry when using dynamic cover type inputs', FATAL)
! <ERROR MSG="must set cover_dataset_entry when using dynamic cover type inputs" STATUS="FATAL">
! When specifying dynamic cover type forcing data inputs, the
! cover_dataset_entry time must be set.
! </ERROR>

! get the model initial time, current model time and cover type dataset time
     model_init_time = get_base_time()
     call get_date(model_init_time, init_yr, init_mo, init_dy, init_hr, init_mn, init_sc)
     call get_date(time, model_yr, model_mo, model_dy, model_hr, model_mn, model_sc)
     cover_entry  = set_date (cover_dataset_entry(1), &
                              cover_dataset_entry(2), &
                              cover_dataset_entry(3), &
                              cover_dataset_entry(4), &
                              cover_dataset_entry(5), &
                              cover_dataset_entry(6))


! calculate the initial timelevel to be read from the netcdf file
     kk = ( model_yr - init_yr ) +     &
          ( cover_dataset_entry(1) - cover_dataset_init_year )
     if (mpp_pe() == mpp_root_pe()) then
       write (outunit,'(a45,i4)')    &
          'The cover type dataset begins with year:', cover_dataset_init_year
       write (outunit,'(a40,i4)')    &
          'The cover type data year being read is:', cover_dataset_entry(1)
       write (outunit,'(a60,i4)')    &
          'This cover type data is mapped to current model year:', model_yr
       write (outunit,'(a42,i4)')    &
          'The model init year:', init_yr    
     endif

! read the initial cover type data from the netcdf file
     call read_data('INPUT/cover_type_field.nc','cover',tmp(:,:,1),  &
        timelevel=kk+1, no_domain=.true.)
     cover_type_input(:,:,1) = int(tmp(:,:,1))
   
! obtain the western boundary in the input data
     call field_size('INPUT/cover_type_field.nc', 'blon', siz)
     allocate (nblon_input(siz(1)))
     call read_data('INPUT/cover_type_field.nc','blon',nblon_input(:),timelevel=1,no_domain=.true.)

! regrid the data to the land model grid
     wb_degrees = (nblon_input(1))
     wb = pi*wb_degrees/180.
     dx = 2.*pi/float(num_lon_input)
     dy = pi/float(num_lat_input)     
     call regrid_discrete_field(cover_type_input(:,:,1), wb, sb, dx, dy, lonb,&
        latb, n_map_cover_types, land(:,:,1), cover_type(:,:,1) )
     deallocate (cover_type_input)


! if using the read_old_ascii_cover option, read the ASCII data file for cover
! type
   elseif (read_old_ascii_cover) then
     call error_mesg('land_properties_init', 'Using static ASCII cover type dataset',NOTE)
! <ERROR MSG="Using the static ASCII cover type dataset." STATUS="NOTE">
!   Using the read_old_ascii_cover option. Reading the ASCII cover type dataset.
! </ERROR>

     call mpp_open (unit, 'INPUT/cover_type_field', action = MPP_RDONLY)
     read (unit,*) num_lon_input, num_lat_input, wb_degrees
     read (unit,*) input_format
     allocate ( integer_input_1 (num_lon_input, num_lat_input) )
     do jj = 1, num_lat_input
        read(unit,input_format) (integer_input_1(ii,jj),ii=1,num_lon_input)
     end do
     close(unit)

     wb = pi*wb_degrees/180.
     dx = 2.*pi/float(num_lon_input)
     dy = pi/float(num_lat_input)
     call regrid_discrete_field(integer_input_1, wb, sb, dx, dy, lonb, latb, &
        n_map_cover_types, land(:,:,1), cover_type(:,:,1) )
     deallocate (integer_input_1)
     call close_file(unit)
  
! if not using dynamic_cover_type nor the read_old_ascii_cover option, use the
! static netcdf cover_type dataset if it exists
   elseif(file_exist('INPUT/cover_type_field.nc')) then

     call error_mesg('land_properties_init', 'Using the static netcdf cover type dataset',NOTE)
! <ERROR MSG="Using the static NetCDF cover type dataset." STATUS="NOTE">
!   Using the static NetCDF cover type dataset.
! </ERROR>

! obtain the size on lon in the input data
     call field_size('INPUT/cover_type_field.nc', 'lon', siz)
     allocate (nlon_input(siz(1)))
     num_lon_input = size(nlon_input)

! obtain the size of lat in the input data
     call field_size('INPUT/cover_type_field.nc', 'lat', siz)
     allocate (nlat_input(siz(1)))
     num_lat_input = size(nlat_input)

! allocate arrays, read data as real and convert back to integer
     allocate (tmp(num_lon_input, num_lat_input, npart)) ! real temp array
     allocate (cover_type_input (num_lon_input, num_lat_input, npart))  ! input data

     if (cover_dataset_entry(1) == 1 .and. &
         cover_dataset_entry(2) == 1 .and. &
         cover_dataset_entry(3) == 1 .and. &
         cover_dataset_entry(4) == 0 .and. &
         cover_dataset_entry(5) == 0 .and. &
         cover_dataset_entry(6) == 0 )     &
      call error_mesg ('land_properties_mod', &
     'must set cover_dataset_entry when using static_cover_type_netcdf', FATAL)
! <ERROR MSG="must set cover_dataset_entry when using time-varying cover type inputs" STATUS="FATAL">
! When specifying static_cover_type_netcdf forcing data inputs, the
! cover_dataset_entry time must be set.
! </ERROR>

! get the model initial time and cover type dataset time
     model_init_time = get_base_time()
     cover_entry  = set_date (cover_dataset_entry(1), &
                              cover_dataset_entry(2), &
                              cover_dataset_entry(3), &
                              cover_dataset_entry(4), &
                              cover_dataset_entry(5), &
                              cover_dataset_entry(6))

     if (mpp_pe() == mpp_root_pe()) then
       write (outunit,'(a45,i4)')      &
              'The cover type dataset begins with year:', cover_dataset_init_year
       write (outunit,'(a40,i4)')      &
              'The cover type data year being read is:', cover_dataset_entry(1)
     endif

! calculate the static timelevel to read from the NetCDF file 
     kk = cover_dataset_entry(1) - cover_dataset_init_year
     call read_data('INPUT/cover_type_field.nc','cover',tmp(:,:,1),timelevel=kk+1, no_domain=.true.)
     cover_type_input(:,:,1) = int(tmp(:,:,1))
   
! obtain the western boundary in the input data
     call field_size('INPUT/cover_type_field.nc', 'blon', siz)
     allocate (nblon_input(siz(1)))
     call read_data('INPUT/cover_type_field.nc','blon',nblon_input(:),timelevel=1, no_domain=.true.)

! regrid the data onto the land grid
     wb_degrees = (nblon_input(1))
     wb = pi*wb_degrees/180.
     dx = 2.*pi/float(num_lon_input)
     dy = pi/float(num_lat_input)     
     call regrid_discrete_field(cover_type_input(:,:,1), wb, sb, dx, dy, lonb,&
        latb, n_map_cover_types, land(:,:,1), cover_type(:,:,1) )
     deallocate (cover_type_input)

! if the static cover type forcing data option is chosen and the NetCDF file
! does not exist, read the GSWP2 cover type field in ASCII format
   else 
     call error_mesg('land_properties_init','Cannot find the netcdf file INPUT/cover_type_field. Using the ASCII file.',NOTE)
! <ERROR MSG="Cannot find the netCDF file INPUT/cover_type_field. Using the ASCII file." STATUS="NOTE">
!   The netCDF cover type field file cannot be found. Using the ASCII file.
! </ERROR>

     call mpp_open (unit, 'INPUT/cover_type_field', action = MPP_RDONLY)
     read (unit,*) num_lon_input, num_lat_input, wb_degrees
     read (unit,*) input_format
     allocate ( integer_input_1 (num_lon_input, num_lat_input) )
     do jj = 1, num_lat_input
        read(unit,input_format) (integer_input_1(ii,jj),ii=1,num_lon_input)
     end do
     close(unit)

     wb = pi*wb_degrees/180.
     dx = 2.*pi/float(num_lon_input)
     dy = pi/float(num_lat_input)
     call regrid_discrete_field(integer_input_1, wb, sb, dx, dy, lonb, latb, &
        n_map_cover_types, land(:,:,1), cover_type(:,:,1) )
     deallocate (integer_input_1)
     call close_file(unit)
  endif
endif

!--- update cover type
where (.not.land(:,:,1)) cover_type(:,:,1) = 0
  glacier(:,:,1) = .false.
where (cover_type(:,:,1).eq.veg_index_ice)
  glacier(:,:,1) = .true.
  cover_type(:,:,1) = veg_index_tundra
endwhere

!--- now re-define glacier mask, if necessary
if (.not.use_glaciers) then
   glacier(:,:,1) = .false.
else if (use_climap) then
  if(use_climap_mcm) then
     call get_climap_glacier_mcm(lonb, latb, gnlon, gnlat, is, is+size(land,1)-1, &
                js, js+size(land,2)-1, critical_glacier_albedo, glacier(:,:,1))
  else
     call get_climap_glacier(lonb, latb, critical_glacier_albedo, glacier(:,:,1))
  endif
  where (.not.land(:,:,1)) glacier(:,:,1) = .false.
endif

!--- next re-define cover type, if necessary, not yet worrying
!--- about glacier type
if (veg_to_use.eq.'constant') then
   cover_type(:,:,1) = 0
   where (land(:,:,1)) cover_type(:,:,1) = veg_index_constant
else if (veg_to_use.eq.'cons_tun') then
   cover_type(:,:,1) = 0
   where (land(:,:,1)) cover_type(:,:,1) = veg_index_tuned
else if (veg_to_use.eq.'cons_ssl') then
   cover_type(:,:,1) = 0
   where (land(:,:,1)) cover_type(:,:,1) = veg_index_ssl_veg
endif
!--- finally, assign appropriate cover type to glacier cells (if any)
if (veg_to_use.ne.'cons_ssl') then
   where (glacier(:,:,1)) cover_type(:,:,1) = veg_index_ice
else
   where (glacier(:,:,1)) cover_type(:,:,1) = veg_index_ssl_ice
endif

! call error_mesg ('land_properties_init', 'cover type defined', NOTE)
!  <NOTE>
!---- Ground (soil) type. To force consistency with glacier, already
! defined, any ice soil types are first re-set to coarse soil.
! ice types are then located on the basis of glacier. note that a distinct
! ice type is used only if soil_to_use='variable'; this is in contrast
! to the analogous treatment of vegetation types.
!  </NOTE>
if (soil_to_use .eq. 'variable') then
   if (.not.file_exist('INPUT/ground_type_field')) &
       call error_mesg('land_properties_init','Cannot find file INPUT/ground_type_field',FATAL)
! <ERROR MSG="Cannot find file INPUT/ground_type_field" STATUS="FATAL">
!   The ground (soil) type field file cannot be found. Provide this file or
!   set up namelist parameters so it is not necessary. To do the latter, set
!   soil_to_use to something other than 'variable' in the namelist
!   land_properties_nml.
! </ERROR>
   call mpp_open (unit, 'INPUT/ground_type_field', action = MPP_RDONLY)
   read (unit,*) num_lon_input, num_lat_input, wb_degrees
   read (unit,*) input_format
   allocate ( integer_input_1 (num_lon_input, num_lat_input) )
   do jj = 1, num_lat_input
      read(unit,input_format) (integer_input_1(ii,jj),ii=1,num_lon_input)
   end do
   close(unit)
   wb = pi*wb_degrees/180.
   dx = 2.*pi/float(num_lon_input)
   dy = pi/float(num_lat_input)
   if (soil_index_ice_substitute.eq.0) then
      where (integer_input_1.eq.soil_index_ice) integer_input_1 = 0
   endif
   call regrid_discrete_field(integer_input_1, wb, sb, dx, dy, lonb, latb, &
        n_map_ground_types, land(:,:,1), ground_type(:,:,1) )
   deallocate (integer_input_1)
   where (.not.land(:,:,1)) ground_type(:,:,1) = 0
   if (soil_index_ice_substitute.ne.0) then
      where (ground_type(:,:,1).eq.soil_index_ice) &
           ground_type(:,:,1) = soil_index_ice_substitute
   endif
   where (glacier(:,:,1)) ground_type(:,:,1) = soil_index_ice

   ! reconcile lake distribution in cover type and ground type fields
   if (reconcile_lakes) then
      where (ground_type(:,:,1).eq.soil_index_lake) &
           ground_type(:,:,1) = soil_index_lake_substitute
      where (cover_type(:,:,1).eq.veg_index_lake) &
           ground_type(:,:,1) = soil_index_lake
   endif
   
   call close_file(unit)
else if (soil_to_use .eq. 'constant') then
   ground_type(:,:,1) = 0
   where (land(:,:,1)) ground_type(:,:,1) = soil_index_constant
else if (soil_to_use .eq. 'cons_ssl') then
   ground_type(:,:,1) = 0
   where (land(:,:,1)) ground_type(:,:,1) = soil_index_ssl
endif

! call error_mesg ('land_properties_init', 'ground type defined', NOTE)

!---- groundwater residence time

if (use_single_geo) then
   where (land(:,:,1)) tau_groundwater(:,:,1) = geo_res_time
else
   if(.not.file_exist('INPUT/groundwater_residence_time_field')) &
      call error_mesg('land_properties_init','Cannot find file INPUT/groundwater_residence_time_field',FATAL)
! <ERROR MSG="Cannot find file INPUT/groundwater_residence_time_field" STATUS="FATAL">
! The groundwater residence time field file cannot be found. Provide this file
! or set up namelist parameters so it is not necessary. To do the latter, set
! use_single_geo to .true. in the namelist land_properties_nml.
! </ERROR>
   call mpp_open (unit,'INPUT/groundwater_residence_time_field', action = MPP_RDONLY)
   read (unit,*) num_lon_input, num_lat_input, wb_degrees
   read (unit,*) input_format
   allocate ( real_input (num_lon_input, num_lat_input) )
   do jj = 1, num_lat_input
      read(unit,input_format) (real_input(ii,jj),ii=1,num_lon_input)
   end do
   call close_file(unit)
   wb = pi*wb_degrees/180.
   dx = 2.*pi/float(num_lon_input)
   dy = pi/float(num_lat_input)
   call create_horiz_interp_new (Interp_gw, num_lon_input, num_lat_input, &
                                 wb, sb, dx, dy, lonb, latb)
   call horiz_interp ( Interp_gw, real_input, tau_groundwater(:,:,1) )
   deallocate (real_input)
endif

! call error_mesg ('land_properties_init', 'groundwater residence time defined', NOTE)

!---- river destination pointers


lake(:,:,1) = .false.
where (ground_type(:,:,1).eq.soil_index_lake) lake(:,:,1) = .true.

do ii = 2, npart
   glacier        (:,:,ii)  = glacier        (:,:,1)
   lake           (:,:,ii)  = lake           (:,:,1)
   cover_type     (:,:,ii)  = cover_type     (:,:,1)
   ground_type    (:,:,ii)  = ground_type    (:,:,1)
   tau_groundwater(:,:,ii)  = tau_groundwater(:,:,1)
end do
   
!---- assign properties that depend only on cover type -----------
veg_rs_min = 0
rough_momentum = 0
do i_cover_type = 1, n_dim_cover_types
   where (cover_type.eq.i_cover_type)
      albedo_min_no_snow  = min_nosnow_alb_vec(i_cover_type)
      albedo_max_no_snow  = max_nosnow_alb_vec(i_cover_type)
      veg_rs_min          = factor_stomata * veg_rs_min_vec    (i_cover_type)
      rough_momentum      = factor_rough   * rough_momentum_vec(i_cover_type)
      rough_heat          = factor_rough   * rough_momentum_vec(i_cover_type) &
           * exp(-k_over_B_vec(i_cover_type))
   endwhere
enddo

!---- compute or read the topographic roughness ----
if (use_topo_rough) then
   call set_topo_rough (domain, lonb, latb, topo_stdev)
   stdev = min(topo_stdev*topo_rough_factor,max_topo_rough)
   do k = 1,size(rough_momentum,3)
      rough_scale(:,:,k) = max(stdev,rough_momentum(:,:,k))
   end do
else
   rough_scale = rough_momentum
endif

!---- assign properties that depend only on ground type ----------
do i_ground_type = 1, n_dim_ground_types
   where (ground_type.eq.i_ground_type)
      soil_therm_con(:,:,:,1)  = soil_therm_dif_vec(i_ground_type) &
                       *soil_therm_cap_vec(i_ground_type)
      soil_therm_cap(:,:,:,1)  = soil_therm_cap_vec(i_ground_type)
   endwhere
enddo

do ii=2,size(soil_therm_con,4)
  where (land(:,:,:))
     soil_therm_con(:,:,:,ii) = soil_therm_con(:,:,:,1)
     soil_therm_cap(:,:,:,ii) = soil_therm_cap(:,:,:,1)
  endwhere
enddo

do ii=1,num_sfc_layers
  where (land(:,:,:))
     soil_therm_con(:,:,:,ii) = sfc_heat_factor * soil_therm_con(:,:,:,ii)
     soil_therm_cap(:,:,:,ii) = sfc_heat_factor * soil_therm_cap(:,:,:,ii)
  endwhere
enddo

!---- assign properties that depend on cover and ground types ----
max_water = 0
do i_cover_type = 1, n_dim_cover_types
   if ((i_cover_type .ne. veg_index_ice) .and. &
        (i_cover_type .ne. veg_index_ssl_ice)) then
      z_root = veg_zeta_vec(i_cover_type) &
           *log( veg_root_mass_vec(i_cover_type)/ &
           (veg_zeta_vec(i_cover_type)*factor_root*critical_root_density) )
      z_root = max(z_root_min, z_root)
      do i_ground_type = 1, n_dim_ground_types
         where (cover_type.eq.i_cover_type.and.ground_type.eq.i_ground_type)
            max_water = z_root * soil_awc_vec(i_ground_type)
         endwhere
      enddo
   endif
enddo
   
!---- change snow-free albedo to CLIMAP array if requested -------
!
if (use_climap) then
   allocate (real_input(nlon, nlat))
   if(use_climap_mcm) then
     call get_climap_albedo_mcm(lonb, latb, gnlon, gnlat, is, is+size(land,1)-1, &
                                js, js+size(land,2)-1, 1, real_input)
   else
     call get_climap_albedo(lonb, latb, 1, real_input)
   endif
   do ii = 1, npart
      albedo_min_no_snow(:,:,ii) = real_input
      albedo_max_no_snow(:,:,ii) = real_input
      if(use_climap_mcm) then
        where(glacier(:,:,ii))
          albedo_min_no_snow(:,:,ii) = min_nosnow_alb_vec(veg_index_ssl_ice)
          albedo_max_no_snow(:,:,ii) = max_nosnow_alb_vec(veg_index_ssl_ice)
        endwhere
      endif
   enddo
   deallocate (real_input)
elseif (use_desert_albedo_map) then
   allocate (real_input(nlon, nlat))
   call get_climap_albedo(lonb, latb, 5, real_input)
   do ii = 1, npart
      where (cover_type(:,:,ii).eq.veg_index_desert)
          albedo_min_no_snow(:,:,ii) = real_input
          albedo_max_no_snow(:,:,ii) = real_input
        endwhere
   enddo
   deallocate (real_input)

endif

max_snow_out = max_snow

! initialize land properties diagnostics
call init_land_properties_diag(id_lon, id_lat, time)

! send static fields for diagnostic manager for output
call diag_static(time, land)


end subroutine land_properties_init
! </SUBROUTINE>


!#######################################################################

! <SUBROUTINE NAME="update_land_properties_slow">

!   <OVERVIEW>
!     Updates slowly changing parameters, such as cover type.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates slowly changing parameters, such as cover type.
!   </DESCRIPTION>

!   <PUBLICROUTINE>
subroutine update_land_properties_slow ( &
     lonb,   latb,  land,  time,  domain, &
     glacier,         &
     lake,            &
     rough_momentum,  &
     rough_heat,      &
     rough_scale,     &  ! topographic roughness for drag scaling
     soil_therm_con,  &
     soil_therm_cap,  &
     veg_rs_min,      &
     max_water,       &
     tau_groundwater, &
     max_snow_out,    &
     cover_type )

  real,            intent(in) :: lonb(:,:), latb(:,:) ! corners of the cells,
                                                  ! radian
  logical,         intent(in) :: land(:,:,:)      ! land mask
  type(time_type), intent(in) :: time             ! current time
  type(domain2d),  intent(in) :: domain           ! our domain 

  logical, intent(out), dimension(:,:,:) :: glacier  ! glacier mask
  logical, intent(out), dimension(:,:,:) :: lake     ! lake mask
  real   , intent(out), dimension(:,:,:) :: &
       rough_momentum, &                  ! roughness length for momentum
       rough_heat,     &                
       rough_scale,    &
       veg_rs_min,     & 
       max_water,      &
       tau_groundwater                    ! groundwater residence time
  real   , intent(out), dimension(:,:,:,:) :: &
       soil_therm_con, &
       soil_therm_cap
  real   , intent(out)  :: max_snow_out
  integer, intent(inout), dimension(:,:,:)  :: cover_type
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer, allocatable :: integer_input_1(:,:)
  real,    allocatable :: real_input     (:,:)
  real,    allocatable :: tmp(:,:,:)  ! temp array for read data for cover_type
  integer, allocatable :: cover_type_input (:,:,:)  ! cover_type input data
  integer, allocatable :: nlon_input (:)     ! input data lon
  integer, allocatable :: nlat_input (:)     ! input data lat
  real,    allocatable :: nblon_input (:)    ! input data western boundary

  real, dimension(size(land,1),size(land,2)):: stdev ! standard deviation of
                                                     ! orography 
  logical got_stdev

  character*80            input_format
  integer :: unit, ierr, i_cover_type, i_ground_type, outunit
  integer :: ii, jj, num_lon_input, num_lat_input, k
  real    :: wb_degrees ! western boundary of the input grid cells
  real    :: sb, wb, dx, dy, z_root
  integer :: is, js     ! boundaries of our domain
  integer :: gnlon, gnlat
  integer :: siz(4)       ! size of field in field_size call
  integer :: model_yr, model_mo, model_dy, model_hr, model_mn, model_sc
  integer :: init_yr, init_mo, init_dy, init_hr, init_mn, init_sc
  integer :: kk           ! time loop integer for cover type data

! call only if using time-varying cover type
if (dynamic_cover_type) then

! Read and re-grid each land input data field

sb = -pi/2.

!  <NOTE>
!---- Cover (vegetation) type. There are 12 possible cases: veg_to_use can
! take 4 values, and glaciers can be (1) absent, (2) based on input cover
! field, or (3) based on climap albedo field. (use of climap albedo forces
! use of climap to locate glaciers, if glaciers are used.) the code
! ensures consistency between cover_type field and glacier mask, with
! the latter taking precedence: where glacier based on input cover field
! must be removed in deference to climap or if no glaciers are used,
! such points are assigned tundra (or global constant) type.
!  </NOTE>

!--- first get input cover field and implied glacier field, if either of
! these will be needed. then re-set ice cover type to tundra; it may be
! changed back to ice according to actual glacier specification later.
! If dynamic_cover_type_netcdf is chosen, then the cover_type will be obtained
! from a netCDF file of annual data. The initial year of the dataset is
! specified by the namelist option cover_dataset_init_year.
if ( veg_to_use.eq.'variable' &
     .or. (use_glaciers.and..not.use_climap) ) then
     if(.not.file_exist('INPUT/cover_type_field.nc')) &
       call error_mesg('land_properties_slow','Cannot find file INPUT/cover_type_field.nc',FATAL)

! <ERROR MSG="Cannot find file INPUT/cover_type_field.nc" STATUS="FATAL">
!   The cover type field file cannot be found. Provide this file or set up
!   namelist parameters so it is not necessary. To do the latter, set
!   veg_to_use to something other than 'variable' in the namelist
!   land_properties_nml.
! </ERROR>

! obtain the size on lon in the input data
     call field_size('INPUT/cover_type_field.nc', 'lon', siz)
     allocate (nlon_input(siz(1)))
     num_lon_input = size(nlon_input)

! obtain the size of lat in the input data
     call field_size('INPUT/cover_type_field.nc', 'lat', siz)
     allocate (nlat_input(siz(1)))
     num_lat_input = size(nlat_input)

! allocate arrays, read data as real and convert back to integer
     allocate (tmp(num_lon_input, num_lat_input, npart)) ! real temp array
     allocate (cover_type_input (num_lon_input, num_lat_input, npart))  ! input data
  
! get the current year of the model
     call get_date(model_init_time, init_yr, init_mo, init_dy, init_hr, init_mn, init_sc)
     call get_date(time, model_yr, model_mo, model_dy, model_hr, model_mn, model_sc)

! calculate the timelevel read from the netcdf file
     kk = ( model_yr - init_yr ) +     &
          ( cover_dataset_entry(1) - cover_dataset_init_year )
     if (mpp_pe() == mpp_root_pe()) then
       outunit = stdout()
       write (outunit,'(a45,i4)')      &
          'The cover type dataset begins with year:', cover_dataset_init_year
       write (outunit,'(a40,i4)')      &
          'The cover type data year being read is:', cover_dataset_entry(1)
       write (outunit,'(a60,i4)')      &
          'This cover type data is mapped to current model year:', model_yr
       write (outunit,'(a42,i4)')      &
          'The model init year:', init_yr     
     endif

! read the cover_type data from the netcdf file and covert it back to integer
     call read_data('INPUT/cover_type_field.nc','cover',tmp(:,:,1),timelevel=kk+1,no_domain=.true.)
     cover_type_input(:,:,1) = int(tmp(:,:,1))

! obtain the western boundary in the input data
     call field_size('INPUT/cover_type_field.nc', 'blon', siz)
     allocate (nblon_input(siz(1)))
     call read_data('INPUT/cover_type_field.nc','blon',nblon_input(:),timelevel=1,no_domain=.true.)

! regrid the data to the land grid and update the cover type
     wb_degrees = (nblon_input(1))
     wb = pi*wb_degrees/180.
     dx = 2.*pi/float(num_lon_input)
     dy = pi/float(num_lat_input)     
     call regrid_discrete_field(cover_type_input(:,:,1), wb, sb, dx, dy, lonb,&
        latb, n_map_cover_types, land(:,:,1), cover_type(:,:,1) )
     deallocate (cover_type_input)
     where (.not.land(:,:,1)) cover_type(:,:,1) = 0
     glacier(:,:,1) = .false.
     where (cover_type(:,:,1).eq.veg_index_ice)
       glacier(:,:,1) = .true.
       cover_type(:,:,1) = veg_index_tundra
     endwhere
endif

!--- next re-define cover type, if necessary, not yet worrying
!--- about glacier type
if (veg_to_use.eq.'constant') then
  cover_type(:,:,1) = 0
  where (land(:,:,1)) cover_type(:,:,1) = veg_index_constant
  else if (veg_to_use.eq.'cons_tun') then
   cover_type(:,:,1) = 0
   where (land(:,:,1)) cover_type(:,:,1) = veg_index_tuned
  else if (veg_to_use.eq.'cons_ssl') then
   cover_type(:,:,1) = 0
   where (land(:,:,1)) cover_type(:,:,1) = veg_index_ssl_veg
endif

!--- finally, assign appropriate cover type to glacier cells (if any)
if (veg_to_use.ne.'cons_ssl') then
  where (glacier(:,:,1)) cover_type(:,:,1) = veg_index_ice
else
  where (glacier(:,:,1)) cover_type(:,:,1) = veg_index_ssl_ice
endif

! call error_mesg ('land_properties_slow', 'cover type defined', NOTE)
!  <NOTE>
!---- Ground (soil) type. Any ice soil types are first re-set to coarse soil.
! ice types are then located on the basis of glacier. note that a distinct
! ice type is used only if soil_to_use='variable'; this is in contrast
! to the analogous treatment of vegetation types.
!  </NOTE>
if (soil_to_use .eq. 'variable') then
  if (.not.file_exist('INPUT/ground_type_field')) &
    call error_mesg('land_properties_slow','Cannot find file INPUT/ground_type_field',FATAL)
! <ERROR MSG="Cannot find file INPUT/ground_type_field" STATUS="FATAL">
!   The ground (soil) type field file cannot be found. Provide this file or
!   set up namelist parameters so it is not necessary. To do the latter, set
!   soil_to_use to something other than 'variable' in the namelist
!   land_properties_nml.
! </ERROR>
  call mpp_open (unit, 'INPUT/ground_type_field', action = MPP_RDONLY)
  read (unit,*) num_lon_input, num_lat_input, wb_degrees
  read (unit,*) input_format
  allocate ( integer_input_1 (num_lon_input, num_lat_input) )
  do jj = 1, num_lat_input
    read(unit,input_format) (integer_input_1(ii,jj),ii=1,num_lon_input)
  end do
  close(unit)

  wb = pi*wb_degrees/180.
  dx = 2.*pi/float(num_lon_input)
  dy = pi/float(num_lat_input)

  if (soil_index_ice_substitute.eq.0) then
    where (integer_input_1.eq.soil_index_ice) integer_input_1 = 0
  endif

  call regrid_discrete_field(integer_input_1, wb, sb, dx, dy, lonb, latb,&
     n_map_ground_types, land(:,:,1), ground_type(:,:,1) )
  deallocate (integer_input_1)
  where (.not.land(:,:,1)) ground_type(:,:,1) = 0
  if (soil_index_ice_substitute.ne.0) then
    where (ground_type(:,:,1).eq.soil_index_ice) &
       ground_type(:,:,1) = soil_index_ice_substitute
  endif
  where (glacier(:,:,1)) ground_type(:,:,1) = soil_index_ice
  call close_file(unit)
  else if (soil_to_use .eq. 'constant') then
    ground_type(:,:,1) = 0
    where (land(:,:,1)) ground_type(:,:,1) = soil_index_constant 
  else if (soil_to_use .eq. 'cons_ssl') then
    ground_type(:,:,1) = 0
    where (land(:,:,1)) ground_type(:,:,1) = soil_index_ssl
  endif

  ! reconcile lake distribution in cover type and ground type fields
  if (reconcile_lakes) then
     where (ground_type(:,:,1).eq.soil_index_lake) &
          ground_type(:,:,1) = soil_index_lake_substitute
     where (cover_type(:,:,1).eq.veg_index_lake) &
          ground_type(:,:,1) = soil_index_lake
  endif

  lake(:,:,1) = .false.
  where (ground_type(:,:,1).eq.soil_index_lake) lake(:,:,1) = .true.

! call error_mesg ('land_properties_slow', 'ground type defined', NOTE)

!---- groundwater residence time

  if (use_single_geo) then
    where (land(:,:,1)) tau_groundwater(:,:,1) = geo_res_time
  else
    if(.not.file_exist('INPUT/groundwater_residence_time_field')) &
     call error_mesg('land_properties_slow','Cannot find file INPUT/groundwater_residence_time_field',FATAL)
! <ERROR MSG="Cannot find file INPUT/groundwater_residence_time_field"
! STATUS="FATAL">
!   The groundwater residence time field file cannot be found. Provide this
!   file or set up namelist parameters so it is not necessary. To do the
!   latter, set use_single_geo to .true. in the namelist land_properties_nml.
! </ERROR>
     call mpp_open (unit,'INPUT/groundwater_residence_time_field', action = MPP_RDONLY)
     read (unit,*) num_lon_input, num_lat_input, wb_degrees
     read (unit,*) input_format
     allocate ( real_input (num_lon_input, num_lat_input) )
     do jj = 1, num_lat_input
       read(unit,input_format) (real_input(ii,jj),ii=1,num_lon_input)
     end do
     call close_file(unit)
     wb = pi*wb_degrees/180.
     dx = 2.*pi/float(num_lon_input)
     dy = pi/float(num_lat_input)
     ! initialized already
     !call create_horiz_interp_new (Interp_gw, num_lon_input, num_lat_input, &
     !                              wb, sb, dx, dy, lonb, latb )
     call horiz_interp ( Interp_gw, real_input, tau_groundwater(:,:,1) )
     deallocate (real_input)
endif

! call error_mesg ('land_properties_slow', 'groundwater residence time defined', NOTE)

!---- river destination pointers

   do ii = 2, npart
    lake           (:,:,ii)  = lake           (:,:,1)
    cover_type     (:,:,ii)  = cover_type     (:,:,1)
    ground_type    (:,:,ii)  = ground_type    (:,:,1)
    tau_groundwater(:,:,ii)  = tau_groundwater(:,:,1)
   end do
   
!---- assign properties that depend only on cover type -----------
   veg_rs_min = 0
   rough_momentum = 0
   do i_cover_type = 1, n_dim_cover_types
   where (cover_type.eq.i_cover_type)
      albedo_min_no_snow  = min_nosnow_alb_vec(i_cover_type)
      albedo_max_no_snow  = max_nosnow_alb_vec(i_cover_type)
      veg_rs_min          = factor_stomata * veg_rs_min_vec    (i_cover_type)
      rough_momentum      = factor_rough   * rough_momentum_vec(i_cover_type)
      rough_heat          = factor_rough   * rough_momentum_vec(i_cover_type) &
           * exp(-k_over_B_vec(i_cover_type))
   endwhere
   enddo

!---- compute or read the topographic roughness ----
if (use_topo_rough) then
   call set_topo_rough (domain, lonb, latb, topo_stdev)
   stdev = min(topo_stdev*topo_rough_factor,max_topo_rough)
   do k = 1,size(rough_momentum,3)
      rough_scale(:,:,k) = max(stdev,rough_momentum(:,:,k))
   end do
else
   rough_scale = rough_momentum
endif

!---- assign properties that depend only on ground type ----------
do i_ground_type = 1, n_dim_ground_types
   where (ground_type.eq.i_ground_type)
      soil_therm_con(:,:,:,1)  = soil_therm_dif_vec(i_ground_type) &
                       *soil_therm_cap_vec(i_ground_type)
      soil_therm_cap(:,:,:,1)  = soil_therm_cap_vec(i_ground_type)
   endwhere
enddo

do ii=2,size(soil_therm_con,4)
  where (land(:,:,:))
     soil_therm_con(:,:,:,ii) = soil_therm_con(:,:,:,1)
     soil_therm_cap(:,:,:,ii) = soil_therm_cap(:,:,:,1)
  endwhere
enddo

do ii=1,num_sfc_layers
  where (land(:,:,:))
     soil_therm_con(:,:,:,ii) = sfc_heat_factor * soil_therm_con(:,:,:,ii)
     soil_therm_cap(:,:,:,ii) = sfc_heat_factor * soil_therm_cap(:,:,:,ii)
  endwhere
enddo

!---- assign properties that depend on cover and ground types ----
max_water = 0
do i_cover_type = 1, n_dim_cover_types
   if ((i_cover_type .ne. veg_index_ice) .and. &
        (i_cover_type .ne. veg_index_ssl_ice)) then
      z_root = veg_zeta_vec(i_cover_type) &
           *log( veg_root_mass_vec(i_cover_type)/ &
           (veg_zeta_vec(i_cover_type)*factor_root*critical_root_density) )
      z_root = max(z_root_min, z_root)
      do i_ground_type = 1, n_dim_ground_types
         where (cover_type.eq.i_cover_type.and.ground_type.eq.i_ground_type)
            max_water = z_root * soil_awc_vec(i_ground_type)
         endwhere
      enddo
   endif
enddo
   
!---- change snow-free albedo to CLIMAP array if requested -------
!
if (use_climap) then
   allocate (real_input(nlon, nlat))
   if(use_climap_mcm) then
     call get_climap_albedo_mcm(lonb, latb, gnlon, gnlat, is, is+size(land,1)-1, &
                                js, js+size(land,2)-1, 1, real_input)
   else
     call get_climap_albedo(lonb, latb, 1, real_input)
   endif
   do ii = 1, npart
      albedo_min_no_snow(:,:,ii) = real_input
      albedo_max_no_snow(:,:,ii) = real_input
      if(use_climap_mcm) then
        where(glacier(:,:,ii))
          albedo_min_no_snow(:,:,ii) = min_nosnow_alb_vec(veg_index_ssl_ice)
          albedo_max_no_snow(:,:,ii) = max_nosnow_alb_vec(veg_index_ssl_ice)
        endwhere
      endif
   enddo
   deallocate (real_input)
elseif (use_desert_albedo_map) then
   allocate (real_input(nlon, nlat))
   call get_climap_albedo(lonb, latb, 5, real_input)
   do ii = 1, npart
      where (cover_type(:,:,ii).eq.veg_index_desert)
          albedo_min_no_snow(:,:,ii) = real_input
          albedo_max_no_snow(:,:,ii) = real_input
        endwhere
   enddo
   deallocate (real_input)

endif

max_snow_out = max_snow
endif

end subroutine update_land_properties_slow
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="update_land_properties_fast">

!   <OVERVIEW>
!     Updates the rapidly changing parameters, such as albedo.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates the rapidly changing parameters. Computes the albedo of
!     hypothetical no-snow and deep-snow surfaces and uses snow mass to blend
!     snow-free and deep-snow albedo values.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_land_properties_fast (snowmass, t_sfc, land, albedo)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_land_properties_fast (snowmass, t_sfc, land, albedo, cover_type)

!-----------------------------------------------------------------------
!
! INPUT
!    snowmass     = mass of snow on the ground (in kg/(m**2))
!    t_sfc        = surface temperature (in degrees kelvin)
!    land         = logical land mask
! 
!  OUTPUT
!    albedo       = surface albedo

real,    intent(in),  dimension(:,:,:) :: snowmass, t_sfc
logical, intent(in) , dimension(:,:,:) :: land
real,    intent(out), dimension(:,:,:) :: albedo
integer, intent(in),  dimension(:,:,:) :: cover_type
!-----------------------------------------------------------------------
!   </PUBLICROUTINE>

integer :: i_cover_type
real    :: t_crit
real   , dimension(size(t_sfc,1),size(t_sfc,2),size(t_sfc,3)) ::  &
                   albedo_no_snow, albedo_hi_snow, blend

!---- albedo update ----------------------------------------------------

!--- compute albedo of hypothetical no-snow and deep-snow surfaces

t_crit = tfreeze - t_range
blend = (t_sfc - t_crit) / t_range
where (blend .lt. 0.) blend = 0.
where (blend .gt. 1.) blend = 1.

where (land(:,:,:))
   albedo_no_snow = albedo_max_no_snow                            &
    + blend * (albedo_min_no_snow - albedo_max_no_snow)
endwhere

do i_cover_type = 1, n_dim_cover_types
  where (cover_type.eq.i_cover_type)
    albedo_hi_snow = max_snow_alb_vec(i_cover_type) + blend    &
      * (min_snow_alb_vec(i_cover_type)                        &
         - max_snow_alb_vec(i_cover_type) )
  endwhere
enddo

!--- use snow mass to blend snow-free and deep-snow albedo values

if (do_mcm_masking) then
    blend = 0.
    do i_cover_type = 1, n_dim_cover_types
      where (cover_type.eq.i_cover_type) &
           blend = 0.5 * sqrt(snowmass/crit_snowmass_vec(i_cover_type))
      enddo
    where (blend .gt. 1.) blend = 1.
  else
    blend = 0.
    do i_cover_type = 1, n_dim_cover_types
      where (cover_type.eq.i_cover_type) &
           blend = snowmass/(snowmass+crit_snowmass_vec(i_cover_type))
      enddo
  endif

where (cover_type.gt. 0) albedo = albedo_no_snow  &
             + (albedo_hi_snow - albedo_no_snow) * blend

end subroutine update_land_properties_fast
! </SUBROUTINE>

! <SUBROUTINE NAME="regrid_discrete_field">

!   <OVERVIEW>
!     Regrids integer index data to any output grid from a uniformly spaced
!     grid using a maximum-area voting scheme.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Regrids integer index data to any output grid from a uniformly spaced
!     grid using a maximum-area voting scheme.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call regrid_discrete_field (data_in, wb_in, sb_in, dlon_in, dlat_in, &
!     lon_out, lat_out, ntype, mask_out, data_out, mask_in)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine regrid_discrete_field (data_in, wb_in, sb_in, dlon_in, dlat_in, &
     lon_out, lat_out, ntype, &   
     mask_out, data_out, mask_in)

!-----------------------------------------------------------------------
!  input:
!  -----
!     data_in     input data; dimensioned by mdim x ndim
!                      stored from south to north
!     wb_in      longitude corresponding to western boundary of box i=1
!     sb_in      latitude corresponding to southern boundary of box j=1
!     dlon_in    x axis grid spacing in degrees of longitude
!     dlat_in    y axis grid spacing in degrees of latitude
!
!     lon_out   longitudes of output data at grid box corners
!                  dimensioned by size(data_out,1)+1 by size(data_out,2)+1
!     lat_out   latitudes of output data at grid box corners
!                  dimensioned by size(data_out,1)+1 by size(data_out,2)+1
!     ntype     number of land cover types specified
!     mask_out  output mask that specifies where the data are defined
!
!  output:
!  ------
!     data_out     output number of land cover type
!
!  optional
!  --------
!     mask_in   input mask;  =0.0 for data points that should not
!               be used or have no data; has the same size as data_in
!
!-----------------------------------------------------------------------
  integer, intent(in)           :: ntype
  integer, intent(in)           :: data_in(:,:)
  real,    intent(in)           :: sb_in,wb_in, dlat_in, dlon_in
  real,    intent(in)           :: lon_out(:,:), lat_out(:,:)
  logical, intent(in)           :: mask_out(:,:)
  integer, intent(out)          :: data_out(:,:)
  real,    intent(in), optional :: mask_in(:,:)
!   </PUBLICROUTINE>

  if (is_latlon(lon_out,lat_out)) then
     call regrid_discrete_field_latlon (data_in, wb_in, sb_in, dlon_in, dlat_in, &
                   lon_out(:,1), lat_out(1,:), ntype, mask_out, data_out, mask_in)
  else
     call regrid_discrete_field_cube   (data_in, wb_in, sb_in, dlon_in, dlat_in, &
                   lon_out(:,:), lat_out(:,:), ntype, mask_out, data_out, mask_in)
  endif

end subroutine regrid_discrete_field
! </SUBROUTINE>


! <SUBROUTINE NAME="regrid_discrete_field_latlon">
subroutine regrid_discrete_field_latlon (data_in, wb_in, sb_in, dlon_in, dlat_in, &
                                         lon_out, lat_out, ntype, &   
                                         mask_out, data_out, mask_in)
  integer, intent(in)           :: ntype
  integer, intent(in)           :: data_in(:,:)
  real,    intent(in)           :: sb_in,wb_in, dlat_in, dlon_in
  real,    intent(in)           :: lon_out(:), lat_out(:)
  logical, intent(in)           :: mask_out(:,:)
  integer, intent(out)          :: data_out(:,:)
  real,    intent(in), optional :: mask_in(:,:)

  ! --- local vars -----------------------------------------------------
  real,    dimension(size(data_in,2)) :: area_in
  real,    dimension(size(lat_out(:))) :: ph
!  real,    dimension(ntype) :: sumtype
!  real,    dimension(size(data_out,1),size(data_out,2)) :: totarea
  integer, dimension(size(data_out,1)) :: is,ie
  integer, dimension(size(data_out,2)) :: js,je,jsave
  real,    dimension(size(data_out,1)) :: fis,fie,facis,facie
  real,    dimension(size(data_out,2)) :: fjs,fje,facjs,facje,fsave

  integer  i,j,mdim,ndim,m360,nlon,nlat
  integer  ismin,iemax, iis,iie,jjs,jje, kk
  real     ratlat,ratlon,hpie,phs,phn, dsph
  real     ffacis,ffacie,ffacjs,ffacje
  logical  flip, enlarge

  integer, allocatable :: data(:,:)
  real,    allocatable :: area(:,:)

  hpie= pi/2.0
  mdim=size(data_in,1); ndim=size(data_in,2)
  m360=int(4.*hpie/dlon_in + 0.001)
  if (mdim < m360) call error_mesg ('regrid_discrete_field',  &
       'inner dimension for input data is too small.', FATAL)
  
!   <ERROR MSG="inner dimension for input data is too small." STATUS="FATAL">
!   </ERROR>  

  !   --- input (hires) grid resolution ---
  ratlat=1.0/dlat_in
  ratlon=1.0/dlon_in

  !   --- area of input (hires) grid boxes ---
  do j=1,ndim
     phs=sb_in + float(j-1)*dlat_in
     phn=sb_in + float(j)  *dlat_in
     phs=min(phs, hpie); phs=max(phs,-hpie)
     phn=min(phn, hpie); phn=max(phn,-hpie)
     dsph=sin(phn)-sin(phs)
     area_in(j) = dsph * dlon_in
  enddo

  !-----------------------------------------------------------------------

  nlon=size(data_out,1); nlat=size(data_out,2)

  !***********************************************************************

  !------ set up latitudinal indexing ------
  !------ make sure output grid goes south to north ------
  
  if (lat_out(1) < lat_out(nlat+1)) then
     ph(1:nlat+1) = lat_out(1:nlat+1)
     flip = .false.
  else
     ph(1:nlat+1) = lat_out(nlat+1:1:-1)
     flip = .true.
  endif

  fjs(1:nlat)=(ph(1:nlat  )-sb_in)*ratlat+1.0
  fje(1:nlat)=(ph(2:nlat+1)-sb_in)*ratlat+1.0

  js=fjs
  where (fjs < 0.) js=js-1
  where (js  < 1 ) js=1

  je=fje
  where (fje < 0.  ) je=je-1
  where (je  > ndim) je=ndim

  facjs=float(js+1)-fjs
  facje=fje-float(je)

  if (flip) then
     jsave=js; js(1:nlat)=jsave(nlat:1:-1)
     jsave=je; je(1:nlat)=jsave(nlat:1:-1)
     fsave=facjs; facjs(1:nlat)=fsave(nlat:1:-1)
     fsave=facje; facje(1:nlat)=fsave(nlat:1:-1)
  endif

  !------ set up longitudinal indexing ------

  fis(1:nlon)=(lon_out(1:nlon  )-wb_in)*ratlon+1.0
  fie(1:nlon)=(lon_out(2:nlon+1)-wb_in)*ratlon+1.0

  is=fis
  where (fis < 0.) is=is-1
  ie=fie
  where (fie < 0.) ie=ie-1

  facis=float(is+1)-fis
  facie=fie-float(ie)
  
  !----- allocate expanded data arrays ------
  
  ismin=min(minval(is),1)
  iemax=max(maxval(ie),mdim)

  allocate (data(ismin:iemax,1:ndim), area(ismin:iemax,1:ndim))

  data(1:mdim,1:ndim) = data_in(1:mdim,1:ndim)

  if (ismin < 1) then
     data(ismin:0, 1:ndim) = data_in(ismin+mdim:mdim, 1:ndim)
  endif
  if (iemax > mdim) then
     data(mdim+1:iemax, 1:ndim) =  data_in(1:iemax-mdim, 1:ndim)
  endif
  
  do j = 1,ndim
     area(:,j) = area_in(j)
  enddo
  
  if (present(mask_in)) then
     area(1:mdim,1:ndim) = area(1:mdim,1:ndim)*mask_in(1:mdim,1:ndim)
     if (ismin < 1)  &
          area(ismin:0, 1:ndim) = area(ismin+mdim:mdim, 1:ndim)
     if (iemax > mdim)  &
          area(mdim+1:iemax, 1:ndim) = area(1:iemax-mdim, 1:ndim)
  endif
  !-----------------------------------------------------------------------

  !     totarea=0.
  do j=1,nlat
     do i=1,nlon
        iis = is(i)
        iie = ie(i)
        jjs = js(j)
        jje = je(j)
        ffacis = facis(i)
        ffacie = facie(i)
        ffacjs = facjs(j)
        ffacje = facje(j)
        kk = 1
522     enlarge = .false.
        if (mask_out(i,j)) call typemax (data(iis:iie,jjs:jje), &
             area(iis:iie,jjs:jje), &
             ffacis,ffacie,ffacjs,ffacje,&
             ntype,data_out(i,j),enlarge)
        !         totarea(i,j)=sum(area(iis:iie,jjs:jje))
        
        if(is(i)-kk.le.1.and.ie(i)+kk.ge.mdim.and. &
             js(j)-kk.le.1.and.je(j)+kk.ge.ndim) goto 523
        
        if(enlarge) then
           if(is(i)-kk.ge.1) iis = is(i)-kk        
           if(ie(i)+kk.le.mdim) iie = ie(i)+kk
           if(js(j)-kk.ge.1) jjs = js(j)-kk
           if(je(j)+kk.le.ndim) jje = je(j)+kk
           ffacis = 1.
           ffacie = 1.
           ffacjs = 1.
           ffacje = 1.          
           kk=kk+1
           go to 522
        end if
523     continue
     enddo
  enddo

!      sumtype = 0.
!      do j=1,nlat
!      do i=1,nlon
!        do kk = 1, ntype
!          if(data_out(i,j).eq.kk) sumtype(kk)=sumtype(kk)+totarea(i,j)
!        end do  
!      enddo
!      enddo
!write(*,*)'ntype=',ntype
!do i = 1, ntype
!  write(*,*)'type,areatype=',i,sumtype(i)
!end do


  deallocate(data, area)

end subroutine regrid_discrete_field_latlon
! </SUBROUTINE>

! <SUBROUTINE NAME="create_horiz_interp_new">

  subroutine create_horiz_interp_new (Interp, nx, ny, wb, sb, dx, dy, lonb_out, latb_out, is_latlon_out)
  type(horiz_interp_type), intent(inout) :: Interp
  integer, intent(in) :: nx, ny
  real,    intent(in) :: wb, sb, dx, dy, lonb_out(:,:), latb_out(:,:)
  logical, intent(in), optional :: is_latlon_out
  real    :: lonb_in(nx+1), latb_in(ny+1), tpi
  integer :: i, j

     tpi = 2.*PI

     ! longitude boundaries
     do i = 1, nx+1
       lonb_in(i) = wb + float(i-1)*dx
     enddo
       if (abs(lonb_in(nx+1)-lonb_in(1)-tpi) < epsilon(lonb_in)) &
               lonb_in(nx+1)=lonb_in(1)+tpi

     ! latitude boundaries
     do j = 2, ny
       latb_in(j) = sb + float(j-1)*dy
     enddo
       latb_in(1)    = -0.5*PI
       latb_in(ny+1) =  0.5*PI

     call horiz_interp_new (Interp, lonb_in, latb_in, &
                                    lonb_out, latb_out, is_latlon_out=is_latlon_out)

  end subroutine create_horiz_interp_new

! </SUBROUTINE>


! <SUBROUTINE NAME="typemax">

!   <OVERVIEW>
!     Calculation of total area of each land cover type within the territory
!     and the determination of the type occupying the biggest area within the
!     territory.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Calculation of total area of each land cover type within the territory
!     and the determination of the type occupying the biggest area within the
!     territory.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call typemax(data,area,facis,facie,facjs,facje,ntype,ntypemax,enlarge)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine typemax(data,area,facis,facie,facjs,facje,ntype,ntypemax,enlarge)

  integer, intent(in)    :: data(:,:)
  real,    intent(in)    :: facis,facie,facjs,facje
  real,    intent(in)    :: area(:,:)
  integer, intent(in)    :: ntype   
  integer, intent(out)   :: ntypemax ! number of land cover type to be assigned
                                     ! to the entire grid cell; 
  logical, intent(inout) :: enlarge  ! says if it's necessary to repeat the
                                     ! process for a larger area
!   </PUBLICROUTINE>

  ! --- local vars -----------------------------------------------------
  integer itype,id,jd,i,j  
  real, dimension(size(area,1),size(area,2)) :: wt
  real, dimension(ntype) :: sumarea     


  id=size(area,1); jd=size(area,2)
  
  wt = area
  wt( 1,:)=wt( 1,:)*facis
  wt(id,:)=wt(id,:)*facie
  wt(:, 1)=wt(:, 1)*facjs
  wt(:,jd)=wt(:,jd)*facje

! Calculation of total area of each land cover type within the territory:
  sumarea = 0.
  do i = 1, id
  do j = 1, jd
    do itype = 1, ntype
      if(int(data(i,j)).eq.itype) sumarea(itype)=sumarea(itype)+wt(i,j)
    end do
  end do
  end do

! Determination of the type occupying the biggest area within the territory:
  ntypemax = 0
    do itype = 1, ntype
      if(maxval(sumarea).eq.sumarea(itype).and.sumarea(itype).gt.0.) &
                ntypemax = itype
    end do

! Check if the land cover type is defined or not (ntypemax = 0  means  
! ocean or undefined type)
  if(ntypemax.eq.0) enlarge = .true.

end subroutine typemax
! </SUBROUTINE>

! <SUBROUTINE NAME="regrid_discrete_field_cube">

!   <OVERVIEW>
!     Regrids integer index data to any output grid from a uniformly spaced
!     grid using a maximum-area voting scheme.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Regrids integer index data to any output grid from a uniformly spaced
!     grid using a maximum-area voting scheme.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call regrid_discrete_field_cube (data_in, wb_in, sb_in, dlon_in, dlat_in, &
!     lon_out, lat_out, ntype, mask_out, data_out, mask_in)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine regrid_discrete_field_cube (data_in, wb_in, sb_in, dlon_in, dlat_in, &
     lon_out, lat_out, ntype, mask_out, data_out, mask_in)

!-----------------------------------------------------------------------
!  input:
!  -----
!     data_in     input data; dimensioned by mdim x ndim
!                      stored from south to north
!     wb_in      longitude corresponding to western boundary of box i=1
!     sb_in      latitude corresponding to southern boundary of box j=1
!     dlon_in    x axis grid spacing in degrees of longitude
!     dlat_in    y axis grid spacing in degrees of latitude
!
!     lon_out   longitudes of output data at grid box boundaries
!                  dimensioned by size(data_out,1)+1
!     lat_out   latitudes of output data at grid box boundaries
!                  dimensioned by size(data_out,2)+1
!     ntype     number of land cover types specified
!     mask_out  output mask that specifies where the data are defined
!
!  output:
!  ------
!     data_out     output number of land cover type
!
!  optional
!  --------
!     mask_in   input mask;  =0.0 for data points that should not
!               be used or have no data; has the same size as data_in
!
!-----------------------------------------------------------------------
  integer, intent(in)           :: ntype
  integer, intent(in)           :: data_in(:,:)
  real,    intent(in)           :: sb_in, wb_in, dlat_in, dlon_in
  real,    intent(in)           :: lon_out(:,:), lat_out(:,:)
  logical, intent(in)           :: mask_out(:,:)
  integer, intent(out)          :: data_out(:,:)
  real,    intent(in), optional :: mask_in(:,:)
!   </PUBLICROUTINE>

  ! --- local vars -----------------------------------------------------
  real :: lonb(2,2), latb(2,2)
  integer :: data1(1,1)
  logical :: mask1(1,1)
  integer :: i, j, np

  integer :: n, mdim, ndim, m360
  real    :: hpie

  type(horiz_interp_type) :: Interp


  ! setup interpolation
  mdim = size(data_in,1)
  ndim = size(data_in,2)

  hpie = pi/2.0
  m360 = int(4.*hpie/dlon_in + 0.001)
  if (mdim < m360) call error_mesg ('regrid_discrete_field_cube',  &
       'inner dimension for input data is too small.', FATAL)
!   <ERROR MSG="inner dimension for input data is too small." STATUS="FATAL">
!   </ERROR>  

  call create_horiz_interp_new (Interp, mdim, ndim, wb_in, sb_in, dlon_in, dlat_in, lon_out, lat_out)

  !-----------------------------------------------------------------------
  ! loop thru number of types
  ! at each grid box find the type with the largest weight (i.e., area)

  call regrid_discrete_field_base (Interp, data_in, ntype, mask_out, data_out, mask_in)
  call horiz_interp_del (Interp)

  !-----------------------------------------------------------------------
  !if (count(mask_out.and.data_out==0) > 0) then
  !   print *, 'Bad discrete points, pe, nbad = ',mpp_pe(),count(mask_out.and.data_out==0)
  !endif
  ! search for missing/bad points
  do j = 1, size(data_out,2)
  do i = 1, size(data_out,1)
     if (mask_out(i,j) .and. data_out(i,j) == 0) then
        ! Initial grid box
        lonb = lon_out(i:i+1,j:j+1)
        latb = lat_out(i:i+1,j:j+1)
        mask1 = mask_out(i,j)
        np = 0
        do
           ! grow the grid box
           call expand_cell (lonb,latb,1.5)
           ! create a new interp type
           call create_horiz_interp_new (Interp, mdim, ndim, wb_in, sb_in, dlon_in, dlat_in, lonb, latb, is_latlon_out=.false.)
           call regrid_discrete_field_base (Interp, data_in, ntype, mask1, data1, mask_in)
           call horiz_interp_del (Interp)
           np = np+1
           if (data1(1,1) /= 0) then
              data_out(i,j)=data1(1,1)
             !print *, 'Fixed discrete point: pe, np = ',mpp_pe(),np
              exit
           endif
           if (np > 10) then
              print *, 'orig lon = ', lon_out(i:i+1,j:j+1)
              print *, 'orig lat = ', lat_out(i:i+1,j:j+1)
              print *, 'curr lon = ', lonb
              print *, 'curr lat = ', latb
              call error_mesg ('regrid_discrete_field_cube', 'failed to fix discrete point', FATAL)
           endif
        enddo
     endif
  enddo
  enddo

  !-----------------------------------------------------------------------

end subroutine regrid_discrete_field_cube
! </SUBROUTINE>

! <SUBROUTINE NAME="regrid_discrete_field_base">

subroutine regrid_discrete_field_base (Intrp, data_in, ntype, mask_out, data_out, mask_in)
  type(horiz_interp_type), intent(inout) :: Intrp
  integer, intent(in)           :: ntype
  integer, intent(in)           :: data_in(:,:)
  logical, intent(in)           :: mask_out(:,:)
  integer, intent(out)          :: data_out(:,:)
  real,    intent(in), optional :: mask_in(:,:)

  ! --- local vars ---
  real,    dimension(size(data_in,1), size(data_in,2))  :: wt_in
  real,    dimension(size(data_out,1),size(data_out,2)) :: wt_max, wt_out
  integer :: n

  ! loop thru number of types
  ! at each grid box find the type with the largest weight (i.e., area)

  data_out = 0
  wt_max = 0.
  do n = 1, ntype
     where (data_in == n)
        wt_in = 1.
     elsewhere
        wt_in = 0.
     endwhere
     if (present(mask_in)) then
        wt_in = wt_in * mask_in
     endif

     call horiz_interp (Intrp, wt_in, wt_out)

     where (mask_out .and. wt_out > wt_max)
        wt_max = wt_out
        data_out = n
     endwhere
  enddo

end subroutine regrid_discrete_field_base

! </SUBROUTINE>

! <SUBROUTINE NAME="set_topo_rough">

 subroutine set_topo_rough (domain, lonb, latb, topo_stdev)
 type(domain2d), intent(in) :: domain
 real,           intent(in),  dimension(:,:) :: lonb, latb
 real,           intent(out), dimension(:,:) :: topo_stdev
 logical :: got_stdev
 integer :: unit, nc
 character(len=128) :: topo_file_base

!   <ERROR MSG="could not read topography data" STATUS="FATAL">
!     get_topog_stdev failed to provide topography variance data.
!   </ERROR>  
!   <ERROR MSG="input file for topography standard deviation ... does not exist" STATUS="FATAL">
!     topo_rough_source is set to 'input', but input file name either
!     not specified or specified incorrectly, so the program cannot 
!     find it.
!   </ERROR>
!   <ERROR MSG="... is not a valid value for topo_rough_source" STATUS="FATAL">
!     specified value of namelist parameter topo_rough_source is invalid; 
!     valid values are 'computed', 'input' or 'input/computed'.
!   </ERROR>

   if (trim(topo_rough_source) == 'computed') then
      got_stdev = get_topog_stdev(lonb,latb,topo_stdev)
      if (.not.got_stdev) &
           call error_mesg ('land_properties_mod', &
           'could not read topography data', FATAL)
   else if (trim(topo_rough_source)=='input' .or. trim(topo_rough_source)=='input/computed') then
      nc = len_trim(topo_rough_file)
      topo_file_base = topo_rough_file(1:nc-3)
      if (file_exist(topo_rough_file) .or. file_exist(trim(topo_file_base)//'.tile1.nc')) then
         call set_domain(domain)
         call read_data (topo_rough_file,topo_rough_var,topo_stdev)
      else if (file_exist(topo_file_base)) then
         unit = open_restart_file(topo_file_base,'read')
         call read_data(unit,topo_stdev)
         call close_file(unit)
      else
         if (trim(topo_rough_source)=='input') then
             call error_mesg('land_properties_mod',                    &
                     'input file for topography standard deviation "'// &
                      trim(topo_rough_file)//'" does not exist', FATAL)
         else
             got_stdev = get_topog_stdev(lonb,latb,topo_stdev)
             if (.not.got_stdev) &
                       call error_mesg ('land_properties_mod', &
                       'could not read topography data', FATAL)
         endif
      endif
   else
      call error_mesg('land_properties_mod','"'//trim(topo_rough_source)//&
           '" is not a valid value for topo_rough_source', FATAL)
   endif

 end subroutine set_topo_rough

! </SUBROUTINE>

! <SUBROUTINE NAME="init_land_properties_diag">

!   <OVERVIEW>
!      Initializes land properties diagnostics.
!   </OVERVIEW>

!   <DESCRIPTION>
!      Initializes land properties diagnostics.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call init_land_properties_diag (id_lon, id_lat, Time)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine init_land_properties_diag (id_lon, id_lat, Time)

  integer,         intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer,         intent(in) :: id_lat  ! ID of land longitude (X) axis
  type(time_type), intent(in) :: Time    ! current time
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: axes(2)

  axes = (/id_lon, id_lat/)
  id_alb_max = register_static_field(diag_mod_name, "albedo_max_no_snow", axes, &
       "max. snow-free land albedo","dimensionless", missing_value=-1.0)
  id_alb_min = register_static_field(diag_mod_name, "albedo_min_no_snow", axes, &
       "min. snow-free land albedo","dimensionless", missing_value=-1.0)
  id_ground_type = register_static_field(diag_mod_name,"ground_type", axes, &
       "land surface ground type", "dimensionless", missing_value = -1.0)

end subroutine init_land_properties_diag
! </SUBROUTINE>

! <DIAGFIELDS>
!   <NETCDF NAME="albedo_max_no_snow" UNITS="dimensionless">
!     Maximum snow-free land albedo
!   </NETCDF>
!   <NETCDF NAME="albedo_min_no_snow" UNITS="dimensionless">
!     Minimum snow-free land albedo
!   </NETCDF>
!   <NETCDF NAME="cover_type" UNITS="dimensionless">
!     Land surface cover type
!   </NETCDF>
!   <NETCDF NAME="ground_type" UNITS="dimensionless">
!     Land surface ground type
!   </NETCDF>
! </DIAGFIELDS>


! <SUBROUTINE NAME="diag_static">

!   <OVERVIEW>
!     Sends static fields to diagnostic output.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Sends static fields to diagnostic output.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call diag_static ( Time,mask )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine diag_static ( Time,mask )

  type(time_type), intent(in) :: Time
  logical,         intent(in) :: mask(:,:,:)
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  logical :: used
  real    :: dummy (size(mask,1), size(mask,2))

  if (id_alb_max > 0) &
       used = send_data ( id_alb_max, albedo_max_no_snow(:,:,1), Time, mask=mask(:,:,1) )
  if (id_alb_min > 0) &
       used = send_data ( id_alb_min, albedo_min_no_snow(:,:,1), Time, mask=mask(:,:,1) )
  if (id_ground_type > 0) then
     dummy = ground_type(:,:,1)
     used = send_data (id_ground_type, dummy, Time, mask=mask(:,:,1))
  endif

end subroutine diag_static
! </SUBROUTINE>


!=============================================================================

! <SUBROUTINE NAME="land_properties_end">

!   <OVERVIEW>
!     Deallocates the land property data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Deallocates the land property data.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call land_properties_end()
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine land_properties_end()
!   </PUBLICROUTINE>

  module_is_initialized =.FALSE.

  deallocate(ground_type, albedo_min_no_snow, albedo_max_no_snow)

end subroutine land_properties_end
! </SUBROUTINE>

end module land_properties_mod

