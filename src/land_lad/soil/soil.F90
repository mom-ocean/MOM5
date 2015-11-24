module soil_mod

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
!   Module containing processes relating to the soil.
! </OVERVIEW>

! <DESCRIPTION>
!   Soil data type is defined and describes the characteristics of the soil.
!   The soil module and the state of soil is initialized. The soil data is
!   updated on the fast and slow time-scale. Contains updates to the "fast"
!   boundary data that the atmosphere sees and the "slow" part of boundary
!   data for the atmosphere. Sends tile-averaged data to the diagnostics
!   manager.
! </DESCRIPTION>

use time_manager_mod,    only: time_type, get_time, get_date, increment_time
use mpp_domains_mod,     only: domain2d, mpp_get_compute_domain
use fms_mod,             only: file_exist,  close_file, read_data,   &
                               error_mesg, FATAL, NOTE, set_domain, mpp_pe,    &
                               check_nml_error, write_version_number,          &
                               open_namelist_file, open_restart_file,          &
                               write_data, stdlog, mpp_root_pe

use fms_io_mod,          only: get_restart_io_mode

use rivers_mod,          only: rivers_type, rivers_init, update_rivers_fast,   &
                               rivers_end, update_rivers_slow,                 &
                               update_rivers_bnd_slow
use land_types_mod,      only: land_data_type

use land_properties_mod, only: land_properties_init,                          &
                               update_land_properties_fast,                   &
                               update_land_properties_slow
use constants_mod,       only: hlf, hlv, tfreeze

use diag_manager_mod,    only: diag_axis_init, register_diag_field, send_data, &
                               register_static_field

use data_override_mod,   only: data_override

#ifdef SCM
use scm_forc_mod,        only: do_specified_tskin, do_specified_land, TSKIN
#endif

use stock_constants_mod, only: ISTOCK_WATER, ISTOCK_HEAT


implicit none
private


! ==== public interfaces =====================================================
public soil_type            ! soil data type
public soil_init            ! initialize the soil module, and the state of soil
public soil_end             ! finish soil calculations

public update_soil_fast     ! fast time-scale soil update
public update_soil_slow     ! slow time-scale soil update
public update_soil_bnd_fast ! updates "fast" boundary data that atmosphere sees
public update_soil_bnd_slow ! updates "slow" part of boundary data for atmosphere

public send_averaged_data   ! sends tile-averaged data to diagnostics manager
public soil_stock_pe        ! calculate and return total amount of requested quantitiy per PE

! ==== end of public interface ===============================================

! <TYPE NAME="soil_type">

!   <DESCRIPTION>
!     Public data type describing soil characteristics.
!   </DESCRIPTION>

type soil_type
!   private                ! opaque type slm, Mar 19 2002: made it public

!   <DATA NAME="domain" TYPE="domain2d" DIM="2">
!     Domain of computation
!   </DATA>
   type(domain2d)::domain ! our domain of computation

   ! geometry of the grid

!   <DATA NAME="is" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="ie" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="js" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="je" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="n_levels" TYPE="integer">
!     Number of levels in the soil
!   </DATA>
!   <DATA NAME="n_tiles" TYPE="integer">
!     Maximum number of tiles in the soil
!   </DATA>
   integer is,ie,js,je;   ! computational domain bounds
   integer n_levels       ! number of levels in the soil
   integer n_tiles        ! max number of tiles in the soil

!   <DATA NAME="area" UNITS="m2" TYPE="real, pointer" DIM="2">
!     Land area for each cell
!   </DATA>
!   <DATA NAME="frac" TYPE="real, pointer" DIM="3">
!     Tile area fraction for each tile
!   </DATA>
!   <DATA NAME="mask" TYPE="logical, pointer" DIM="3">
!     Land mask: true where data are present 
!   </DATA>
!   <DATA NAME="dz" UNITS="m" TYPE="real, pointer" DIM="1">
!     Thickness of layers
!   </DATA>
!   <DATA NAME="z" UNITS="m" TYPE="real, pointer" DIM="1">
!     Full level
!   </DATA>
!   <DATA NAME="zz" UNITS="m" TYPE="real, pointer" DIM="1">
!     Half level
!   </DATA>
!   <DATA NAME="max_fusion" UNITS="J/m3" TYPE="real, pointer" DIM="1">
!     Max amount of energy stored in frozen water, J/m3
!   </DATA>
   real,    pointer :: area(:,:) =>NULL()    ! land area for each cell
   real,    pointer :: frac(:,:,:) =>NULL()  ! tile area fraction for each tile
   logical, pointer :: mask(:,:,:) =>NULL()  ! land mask: true where data are present
   real,    pointer :: dz(:) =>NULL()        ! thickness of layers, m
   real,    pointer :: z(:) =>NULL(),  zz(:) =>NULL() ! full (z) and half (zz) levels, m

   real,    pointer :: max_fusion(:) =>NULL() ! max amount of fusion energy, J/m3
!   <NOTE>
!     NOTE: original intent was to use frac as a mask, meaning that the land 
!     does not exist where fractional area is zero. However, original version
!     of exchange grid requires sum of fractional areas to be == 1 for every
!     grid point. In addition, it is probably useful to have separate mask,
!     just in case somebody ever needs a tile with zero fractional area.
!   </NOTE>

!   <DATA NAME="time" TYPE="time_type">
!     Current time
!   </DATA>
!   <DATA NAME="dt" UNITS="s" TYPE="real">
!     Fast time step
!   </DATA>
!   <DATA NAME="dt_slow" UNITS="s" TYPE="real">
!     Slow time step
!   </DATA>
   type(time_type) :: time         ! current time
   real            :: dt           ! fast time step, s
   real            :: dt_slow      ! slow time step, s

!   <DATA NAME="temp" UNITS="K" TYPE="real, pointer" DIM="4">
!     Temperature (i,j,tile,level)
!   </DATA>
!   <DATA NAME="fusion" UNITS="J/m3" TYPE="real, pointer" DIM="4">
!      Energy needed to melt frozen freezable water (i,j,tile,level)
!   </DATA>
   ! soil state variables
   real, pointer :: temp(:,:,:,:) =>NULL() !  temperature (i,j,tile,level), degK
   real, pointer :: fusion(:,:,:,:) =>NULL()  ! energy needed to melt frozen freezable water, J/m3
 
!   <DATA NAME="snow" UNITS="kg/m2" TYPE="real, pointer" DIM="3">
!     Snow water mass (i,j,tile)
!   </DATA>
!   <DATA NAME="water" UNITS="kg/m2" TYPE="real, pointer" DIM="3">
!     Root-zone water (i,j,tile)
!   </DATA>
!   <DATA NAME="groundwater" UNITS="kg/m2" TYPE="real, pointer" DIM="3">
!     Groundwater (i,j,tile)
!   </DATA>
!   <DATA NAME="drainage" UNITS="kg/(m2 s)" TYPE="real, pointer" DIM="3">
!     Flux from water to groundwater (i,j,tile)
!   </DATA>
!   <DATA NAME="drainage_accum" UNITS="kg/m2" TYPE="real, pointer" DIM="3">
!     Drainage, accumulated for slow time steps (i,j,tile)
!   </DATA>
!   <DATA NAME="calving" UNITS="kg/(m2 s)" TYPE="real, pointer" DIM="3">
!     Horizontal snow flux divergence (i,j,tile)
!   </DATA>
!   <DATA NAME="calving_accum" UNITS="kg/m2" TYPE="real, pointer" DIM="3">
!     Calving, accumulated for slow time steps (i,j,tile)
!   </DATA>
!   <DATA NAME="albedo" TYPE="real, pointer" DIM="3">
!     Snow-adjusted land (soil?) albedo (i,j,tile)
!   </DATA>
!   <DATA NAME="albedo_vis_dir" TYPE="real, pointer" DIM="3">
!     Direct visible surface albedo (i,j,tile)
!   </DATA>
!   <DATA NAME="albedo_nir_dir" TYPE="real, pointer" DIM="3">
!     Direct near-infrared surface albedo (i,j,tile)
!   </DATA>
!   <DATA NAME="albedo_vis_dif" TYPE="real, pointer" DIM="3">
!     Diffuse visible surface albedo (i,j,tile)
!   </DATA>
!   <DATA NAME="albedo_nir_dif" TYPE="real, pointer" DIM="3">
!     Diffuse near-infrared albedo (i,j,tile)
!   </DATA>
!   <DATA NAME="rough_mom" UNITS="m" TYPE="real, pointer" DIM="3">
!     Momentum roughness length (i,j,tile)
!   </DATA>
!   <DATA NAME="rough_heat" UNITS="m" TYPE="real, pointer" DIM="3">
!     Heat roughness length (i,j,tile)
!   </DATA>
!   <DATA NAME="rough_scale" UNITS="m" TYPE="real, pointer" DIM="3">
!     Scale momentum drag coefficient (i,j,tile)
!   </DATA>
!   <DATA NAME="stomatal" UNITS="s/m" TYPE="real, pointer" DIM="3">
!     Non-water-stressed bulk stomatal resistance (i,j,tile)
!   </DATA>
!   <DATA NAME="max_water" UNITS="kg/m2" TYPE="real, pointer" DIM="3">
!     Water capacity of root zone (i,j,tile)
!   </DATA>
!   <DATA NAME="tcon" UNITS="W/(m K s)" TYPE="real, pointer" DIM="4">
!     Thermal conductivity of the ground (i,j,tile,level)
!   </DATA>
!   <DATA NAME="rho_cap" UNITS="J/(m3 K)" TYPE="real, pointer" DIM="4">
!     Volumetric heat capacity of the ground (i,j,tile,level)
!   </DATA>
!   <DATA NAME="tau_groundwater" UNITS="s" TYPE="real, pointer" DIM="3">
!     Groundwater residence time (i,j,tile)
!   </DATA>
   real, pointer, dimension(:,:,:) :: &  ! (lon, lat, tile)
        snow =>NULL(),           & ! snow water mass, kg/m2
        water =>NULL(),          & ! root-zone water, kg/m2
        groundwater =>NULL(),    & ! groundawter, kg/m2
        drainage =>NULL(),       & ! flux from water to groundwater
        drainage_accum =>NULL(), & ! drainage, accumulated for slow time steps
        calving =>NULL(),        & ! horizontal snow flux divergence
        calving_accum =>NULL(),  & ! calving, accumulated for slow time steps
        albedo =>NULL(),         & ! snow-adjusted land (soil?) albedo -- should be revised
        albedo_vis_dir =>NULL(),     &
        albedo_nir_dir =>NULL(),     &
        albedo_vis_dif =>NULL(),     &
        albedo_nir_dif =>NULL(),     &
        rough_mom =>NULL(),      & ! momentum roughness length
        rough_heat =>NULL(),     & ! heat roughness length
        rough_scale =>NULL(),    & ! scale for momentum drag coefficient
        stomatal =>NULL(),       & ! non-water-stressed bulk stomatal resistance
        max_water =>NULL(),      & ! water capacity of root zone, kg/m2
        tau_groundwater =>NULL()   ! groundwater
   real, pointer, dimension (:,:,:,:) :: &  ! (lon, lat, tile, level)
        tcon =>NULL(),           & ! thermal conductivity of the ground
        rho_cap =>NULL()           ! volumetric heat capacity of the ground, J/(m3 K)

!   <DATA NAME="glacier" TYPE="logical, pointer" DIM="3">
!     Glacier mask, true if land is a glacier (lon, lat, tile)
!   </DATA>
!   <DATA NAME="lake" TYPE="logical, pointer" DIM="3">
!     Lake mask, true where lake on land (lon, lat, tile)
!   </DATA>
   logical, pointer, dimension (:,:,:) :: & ! (lon, lat, tile)
        glacier =>NULL(),        & ! true if land is a glacier
        lake =>NULL()              ! lake mask (true where lake on land)

!   <DATA NAME="cover_type" TYPE="integer, pointer" DIM="3">
!     Land surface cover type (lon, lat, tile)
!   </DATA>
!   <DATA NAME="ground" TYPE="integer, pointer" DIM="3">
!     Ground type (lon, lat, tile)
!   </DATA>
   integer, pointer, dimension (:,:,:) :: & ! (lon, lat, tile)
        cover_type =>NULL(),          & ! land cover type
        ground =>NULL()            ! ground type

!   <DATA NAME="max_snow" TYPE="real">
!     Maximum snow capacity
!   </DATA>
!   <DATA NAME="conserve_glacier_mass" TYPE="logical">
!     When true, prevent glacier sublimation and do not allow glacier melt
!     to run off
!   </DATA>
!   <DATA NAME="aa" TYPE="real, pointer" DIM="1">
!     Coefficients for diffusion
!   </DATA>
!   <DATA NAME="cc" TYPE="real, pointer" DIM="1">
!     Coefficients for diffusion
!   </DATA>
   real :: max_snow       ! max snow capacity
   logical :: conserve_glacier_mass ! when true, prevent glacier sublimation and don't
                                    ! allow glacier melt to run off
   logical :: conserve_glacier_snow_mass ! when true, don't reset negative snow
                                    ! masses to zero on top of glacier
   real, pointer :: aa(:) =>NULL(), cc(:) =>NULL() ! coefficients for diffusion

!   <DATA NAME="Rivers" TYPE="rivers_type">
!     State of rivers and runoff
!   </DATA>
   type(rivers_type) :: Rivers ! state of rivers and runoff

end type soil_type
! </TYPE>


! some names, for information only
character(len=*),    parameter :: module_name = 'soil'
character(len=128),  parameter :: version     = '$Id: soil.F90,v 16.0 2008/07/30 22:30:32 fms Exp $'
character(len=128),  parameter :: tagname         = '$Name: tikal $'

! ---- module constants ------------------------------------------------------
integer, parameter :: max_lev = 50

! <NAMELIST NAME="soil_nml">
!   <DATA NAME="n_levels" TYPE="integer" DEFAULT="5">
!     Default number of soil levels
!   </DATA>
!   <DATA NAME="dz(max_lev)" UNITS="m" TYPE="real" DIM="1" DEFAULT="(/ 0.005, 0.045, 0.10, 0.35, 1.0, 45*0. /)">
!     Default thickness of layers, from top down
!   </DATA>
!   <DATA NAME="init_temp_no_ice" UNITS="K" TYPE="real" DEFAULT="288">
!     Initial no-ice T if no restart
!   </DATA>
!   <DATA NAME="init_temp_ice" UNITS="K" TYPE="real" DEFAULT="260">
!     Initial ice T if no restart, used where snow or glacier
!   </DATA>
!   <DATA NAME="init_water" TYPE="real" DEFAULT="0">
!     Initial soil water is set to min(init_water, max_water) if no restart
!   </DATA>
!   <DATA NAME="init_snow" TYPE="real" DEFAULT="0">
!     Initial snow pack is set to min(init_snow, max_snow) if no restart
!   </DATA>
!   <DATA NAME="init_groundwater" TYPE="real" DEFAULT="0">
!     Initial groundwater storage if no restart
!   </DATA>
!   <DATA NAME="conserve_glacier_mass" TYPE="logical" DEFAULT=".false.">
!     When true, prevent glacier sublimation and don't allow glacier melt to
!     run off
!   </DATA>
!   <DATA NAME="freezable_water(max_lev)" UNITS="kg/m3" TYPE="real" DIM="1" DEFAULT="(/0,index=1,max_lev/)">
!     Amount of freezable water in soil, from top down
!   </DATA>
!   <DATA NAME="init_frozen_water" UNITS="kg/m3" TYPE="real" DEFAULT="0">
!     Initial amount of frozen water in soil if no restart
!   </DATA>
integer :: index
! ---- namelist variables and their default values ---------------------------
logical :: do_netcdf_restart = .true.
integer :: n_levels    = 5 ! default number of soil levels
real    :: dz(max_lev) = & ! default thickness of layers, from top down, m
     (/ 0.005, 0.045, 0.10, 0.35, 1.0, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., &
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
real    :: init_temp_no_ice     = 288.        ! initial no-ice T if no restart
real    :: init_temp_ice        = 260.        ! initial ice T if no restart,
                                              ! used where snow or glacier
real    :: init_water           =   0.        ! initial soil water is set to
                                              ! min(init_water, max_water)
                                              ! if no restart
real    :: init_snow            =   0.        ! initial snow pack is set to
                                              ! min(init_snow, max_snow)
                                              ! if no restart
real    :: init_groundwater     =   0.        ! initial groundwater storage
                                              ! if no restart
logical :: conserve_glacier_mass= .false.     ! when true, prevent glacier sublimation
                                              ! and don't allow glacier melt to run off
logical :: conserve_glacier_snow_mass= .false.   ! when true, don't reset negative snow
                                              ! masses to zero on top of glacier
real    :: freezable_water(max_lev) = (/(0.0,index=1,max_lev)/)      
         ! amount of water that can be frozen (but cannot participate in other hydrologycal processes)
real    :: init_frozen_water   =  0 

namelist /soil_nml/ n_levels, dz, &
     init_temp_no_ice, init_temp_ice, init_water, init_snow, &
     init_groundwater, conserve_glacier_mass, conserve_glacier_snow_mass, &
     freezable_water, init_frozen_water
! </NAMELIST>


logical :: module_is_initialized = .FALSE. ! indicates whether the module was initialized


! ---- diagnostic field ids --------------------------------------------------
! diag fields
integer :: id_zhalf, id_zfull
integer :: id_snow, id_water, id_temp
integer :: id_albedo, id_albedo_vis_dir, id_albedo_nir_dir
integer ::            id_albedo_vis_dif, id_albedo_nir_dif
integer :: id_drainage, id_calving, id_precip, id_snowfall, id_evapor, id_sublim
integer :: id_smelt, id_gmelt, id_watsno
integer :: id_sens, id_lhf
integer :: id_lw, id_sw
integer :: id_surface_runoff, id_surface_runoff_snow
integer :: id_groundwater
! static fields
integer :: id_area
integer :: id_lfrac
integer :: id_glacier
integer :: id_cover_type
integer :: id_tcon, id_rho_cap
integer :: id_stomatal
integer :: id_max_water
integer :: id_rough_heat, id_rough_mom
integer :: id_tau_gw
integer :: id_frozen
integer :: id_hlf
integer :: id_hlv
! <DIAGFIELDS>
!   <NETCDF NAME="zhalf" UNITS="m">
!     Half level
!   </NETCDF>
!   <NETCDF NAME="zfull" UNITS="m">
!     Full level
!   </NETCDF>
!   <NETCDF NAME="water" UNITS="kg/m2">
!     Mass of water in the bucket
!   </NETCDF>
!   <NETCDF NAME="snow" UNITS="kg/m2">
!     Mass of snow on the ground
!   </NETCDF>
!   <NETCDF NAME="temp" UNITS="degK">
!     Temperature
!   </NETCDF>
!   <NETCDF NAME="albedo" UNITS="dimensionless">
!     Albedo
!   </NETCDF>
!   <NETCDF NAME="albedo_vis_dir" UNITS="dimensionless">
!     Direct visible surface albedo
!   </NETCDF>
!   <NETCDF NAME="albedo_nir_dir" UNITS="dimensionless">
!     Direct near-infrared surface albedo
!   </NETCDF>
!   <NETCDF NAME="albedo_vis_dif" UNITS="dimensionless">
!     Diffuse visible surface albedo
!   </NETCDF>
!   <NETCDF NAME="albedo_nir_dif" UNITS="dimensionless">
!     Diffuse near-infrared surface albedo
!   </NETCDF>
!   <NETCDF NAME="drainage" UNITS="kg/(m2 s)">
!     Drainage rate
!   </NETCDF>
!   <NETCDF NAME="calving" UNITS="kg/(m2 s)">
!     Snow sweep rate
!   </NETCDF>
!   <NETCDF NAME="precip" UNITS="kg/(m2 s)">
!     Total precipitation rate
!   </NETCDF>
!   <NETCDF NAME="snowfall" UNITS="kg/(m2 s)">
!     Snowfall rate
!   </NETCDF>
!   <NETCDF NAME="evap" UNITS="kg/(m2 s)">
!     Evaporation rate
!   </NETCDF>
!   <NETCDF NAME="sublim" UNITS="kg/(m2 s)">
!     Sublimation rate
!   </NETCDF>
!   <NETCDF NAME="smelt" UNITS="kg/(m2 s)">
!     Snow melt rate
!   </NETCDF>
!   <NETCDF NAME="gmelt" UNITS="kg/(m2 s)">
!     Glacier melt rate
!   </NETCDF>
!   <NETCDF NAME="sens" UNITS="W/m2">
!     Sensible heat flux
!   </NETCDF>
!   <NETCDF NAME="latent" UNITS="W/m2">
!     Latent heat flux
!   </NETCDF>
!   <NETCDF NAME="flw" UNITS="W/m2">
!     Net longwave radiative flux
!   </NETCDF>
!   <NETCDF NAME="fsw" UNITS="W/m2">
!     Net shortwave radiative flux
!   </NETCDF>
!   <NETCDF NAME="wroff" UNITS="kg/(m2 s)">
!     Surface runoff of water
!   </NETCDF>
!   <NETCDF NAME="sroff" UNITS="kg/(m2 s)">
!     Surface runoff of snow
!   </NETCDF>
!   <NETCDF NAME="groundwater" UNITS="kg/m2">
!     Mass of water below bucket
!   </NETCDF>
!   <NETCDF NAME="area" UNITS="m2">
!     Land area
!   </NETCDF>
!   <NETCDF NAME="lfrac" UNITS="dimensionless">
!     Fraction of land in cell
!   </NETCDF>
!   <NETCDF NAME="glacier" UNITS="dimensionless">
!     Glacier logical mask
!   </NETCDF>
!   <NETCDF NAME="cover_type" UNITS="dimensionless">
!     Land surface cover type
!   </NETCDF>
!   <NETCDF NAME="tcon" UNITS="W/(m K s)">
!     Ground thermal conductivity
!   </NETCDF>
!   <NETCDF NAME="rho_cap" UNITS="J/(m3 K)">
!     Ground volumetric heat capacity
!   </NETCDF>
!   <NETCDF NAME="stomatal" UNITS="s/m">
!     Non-water-stressed stomatal resistance
!   </NETCDF>
!   <NETCDF NAME="max_water" UNITS="kg/m2">
!     Root-zone water capacity
!   </NETCDF>
!   <NETCDF NAME="rough_mom" UNITS="m">
!     Momentum roughness length
!   </NETCDF>
!   <NETCDF NAME="rough_heat" UNITS="m">
!     Scalar roughness length
!   </NETCDF>
!   <NETCDF NAME="tau_gw" UNITS="s">
!     Groundwater residence time
!   </NETCDF>
!   <NETCDF NAME="frozen" UNITS="kg/m3">
!     Amount of frozen water in soil
!   </NETCDF>
!   <NETCDF NAME="hlf" UNITS="J/kg">
!     Latent heat of fusion
!   </NETCDF>
!   <NETCDF NAME="hlv" UNITS="J/kg">
!     Latent heat of evaporation
!   </NETCDF>
! </DIAGFIELDS>
! ---- end of diagnostic field ids -------------------------------------------

! <INTERFACE NAME="send_averaged_data">

!   <OVERVIEW>
!     Interface to tile-averaged diagnostic routines
!   </OVERVIEW>

!   <DESCRIPTION>
!     Interface to tile-averaged diagnostic routines
!   </DESCRIPTION>
! ---- interface to tile-averaged diagnostic routines ------------------------
interface send_averaged_data
   module procedure send_averaged_data2d
   module procedure send_averaged_data3d
end interface
! </INTERFACE>

! for diagnostics only
integer :: i 
integer, parameter :: iwatch = 133, jwatch=84

contains


! <SUBROUTINE NAME="soil_init">

!   <OVERVIEW>
!     Initializes the state of the soil.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Reads the namelist, which is assumed to be named soil_nml and located
!     in the file input.nml. Sets the number of tiles in the soil and the
!     level number for the soil data. Allocate soil data.  Sets up
!     time-related values and vertical layering parameters. Sets up initial
!     land area, fractional area, and mask. Reconciles the fractional area with
!     mask: the sum of fractional areas has to be 1 everywhere, even where
!     there is no land. Initializes accumulated values, soil diagnostics and
!     rivers submodule. Assigns initial values to dynamic land properties
!     (albedo). Send static fields to diagnostic manager for output and copies
!     glacier mass conservation flag to the soil data stricture.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call soil_init ( soil, gblon, gblat, garea, gfrac,  &
!     time, dt_fast, dt_slow, domain, frac, mask, &
!     id_lon, id_lat, dz1 )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine soil_init ( soil, gblon, gblat, garea, gfrac,  &
     time, dt_fast, dt_slow, domain, frac, mask, &
     id_lon, id_lat, dz1 )

  type(soil_type), intent(inout) :: soil        ! soil data to initialize
  real,            intent(in)  :: gblon(:,:)  ! lon corners of the grid cells
  real,            intent(in)  :: gblat(:,:)  ! lat corners of the grid cells
  real,            intent(in)  :: garea(:,:)  ! full area of the land grid cells
  real,            intent(in)  :: gfrac(:,:)  ! fraction of grid cells covered by land
  type(time_type), intent(in)  :: time        ! time origin
  type(time_type), intent(in)  :: dt_fast     ! fast time step
  type(time_type), intent(in)  :: dt_slow     ! slow time step
  type(domain2d),  intent(in)  :: domain      ! domain assigned for us
  ! may be the values below should be replaced with something more
  ! general, like "grid_type" ?
  real,            intent(in)  :: frac(:,:,:) ! fractional area of tiles
  logical,         intent(in)  :: mask(:,:,:) ! true if land
  integer,         intent(in)  :: id_lon      ! ID of longitude land diag axis
  integer,         intent(in)  :: id_lat      ! ID of latitude land diag axis

  real, optional,  intent(in)  :: dz1(:)      ! thickness of soil layers
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: unit         ! unit for various i/o
  integer :: ierr         ! error code, returned by i/o routines
  integer :: io           ! i/o status for the namelist
  integer :: day, sec     ! for computation of real time step
  integer :: k            ! level iterator
  real    :: init_fusion  ! initial value stored in frozen water, J/m3
  integer :: model_yr, model_month, model_day        ! used for get_date
  integer :: model_hour, model_minute, model_second  ! used for get_date
  integer :: this_year    ! the current model year
  character(len=64) :: fname = 'INPUT/soil.res.nc'
  character(len=64) :: lvltag
  logical :: is_data_read = .false.

  module_is_initialized =.TRUE.   
   
  ! read namelist, which is assumed to be named soil_nml and located in the 
  ! file input.nml
  if ( file_exist( 'input.nml' ) ) then
     unit = open_namelist_file ( )
     ierr = 1
     do while ( ierr /= 0 )
        read ( unit,  nml = soil_nml, iostat = io, end = 10 ) 
        ierr = check_nml_error ( io, 'soil_nml' )
     enddo
10   continue
     call close_file( unit )
  endif
  call get_restart_io_mode(do_netcdf_restart)

  ! write version information to a log file
  call write_version_number(version, tagname)

  !  write the namelist to a log file
  if( mpp_pe()==0 ) then
     unit = stdlog( )
     write (unit, nml=soil_nml)
     call close_file (unit)
  endif
  
  ! copy specified domain to our data
  soil % domain = domain

  ! get the size of our domain
  call mpp_get_compute_domain ( soil%domain, &
       soil%is, soil%ie, soil%js, soil%je )

  ! set the number of tiles in the soil
  soil%n_tiles = size(frac,3)

  ! set up level number for the soil data
  if (present(dz1)) then
     soil%n_levels = size(dz1(:)) ! values from the argument list
  else
     soil%n_levels = n_levels ! default value
  endif

  ! allocate soil data 
  allocate ( &
       soil% dz             (soil%n_levels),   &
       soil% z              (soil%n_levels),   &
       soil% zz             (soil%n_levels+1), &
       soil% aa             (soil%n_levels+1), &
       soil% cc             (soil%n_levels+1), &
       soil% max_fusion     (soil%n_levels),   &
       soil% area           (soil%is:soil%ie, soil%js:soil%je), & 
       soil% frac           (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), & 
       soil% mask           (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), & 
       soil% temp           (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles,  &
                             soil%n_levels), &
       soil% fusion         (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles,  &
                             soil%n_levels), &
       soil% snow           (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% water          (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% groundwater    (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% drainage       (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% drainage_accum (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% calving        (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% calving_accum  (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% albedo         (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% albedo_vis_dir (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% albedo_nir_dir (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% albedo_vis_dif (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% albedo_nir_dif (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% rough_mom      (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% rough_heat     (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% rough_scale    (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% stomatal       (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% max_water      (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), & 
       stat = ierr )
  if (ierr /= 0) &
       call error_mesg (module_name, 'Cannot allocate memory for soil', FATAL)
  allocate ( &
       soil% tcon           (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles,  &
                             soil%n_levels), & 
       soil% rho_cap        (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles,  &
                             soil%n_levels), & 
       soil% tau_groundwater(soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% glacier        (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% lake           (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% cover_type     (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), &
       soil% ground         (soil%is:soil%ie, soil%js:soil%je, soil%n_tiles), stat = ierr )

  if (ierr /= 0) &
       call error_mesg (module_name, 'Cannot allocate memory for soil', FATAL)
!   <ERROR MSG="Cannot allocate memory for soil" STATUS="FATAL">
!   </ERROR>

  ! set up time-related values
  soil % time = time
  call get_time(dt_fast, sec, day); soil%dt      = day*86400.0+sec
  call get_time(dt_slow, sec, day); soil%dt_slow = day*86400.0+sec

  ! get current model year from the soil%time
  call get_date(soil%Time, model_yr, model_month, model_day, model_hour, model_minute, model_second)
  this_year = model_yr

  ! set up level number for the soil data
  if (present(dz1)) then
     soil%dz       = dz1(1:size(soil%dz(:)))
  else
     soil%dz       = dz (1:size(soil%dz(:)))
  endif

  ! set up vertical layering parameters
  soil%zz(1) = 0
  do k = 1, soil%n_levels
     soil%zz(k+1) = soil%zz(k) + soil%dz(k)
  end do
  do k = 1, soil%n_levels
     soil%z(k) = 0.5*(soil%zz(k+1) + soil%zz(k))
  end do

  do k = 1, soil%n_levels-1
     soil%aa(k) = - 2.0*soil%dt/(soil%dz(k)*(soil%dz(k+1)+soil%dz(k)))
  end do

  do k = 2, soil%n_levels
     soil%cc(k) = - 2.0*soil%dt/(soil%dz(k)*(soil%dz(k)+soil%dz(k-1)))
  end do

  soil%aa(soil%n_levels) = 0.0
  soil%cc(1)             = 0.0

  ! set up initial land area, fractional area, and mask
  soil%area = garea(soil%is:soil%ie, soil%js:soil%je)       &
       * gfrac(soil%is:soil%ie, soil%js:soil%je)
  soil%frac = frac
  soil%mask = mask

  ! initialize max amount of possible fusion 
  soil%max_fusion = freezable_water(1:size(soil%dz(:)))*hlf

  ! read restart file, if it exists
  if( file_exist('INPUT/soil.res.nc') .or. file_exist('INPUT/soil.res.tile1.nc') ) then
     if (mpp_pe() == mpp_root_pe()) call error_mesg ('soil_mod', &
            'Reading NetCDF formatted restart file: INPUT/soil.res.nc', NOTE)
     call set_domain( soil%domain )
     call read_data ( fname, 'frac', soil%frac)
     do k = 1, size(soil%temp, 4)
        if(k < 10)   write(lvltag, '(i1)') k
        if(k >= 10)  write(lvltag, '(i2)') k
        if(k >= 100) write(lvltag, '(i3)') k
        call read_data ( fname, 'soil_temp_'//trim(lvltag), soil%temp(:,:,:,k) )
     enddo
     call read_data ( fname, 'water', soil%water )
     call read_data ( fname, 'snow', soil%snow )
     call read_data ( fname, 'groundwater', soil%groundwater )
     do k = 1, size(soil%fusion, 4)
        if(k < 10)   write(lvltag, '(i1)') k
        if(k >= 10)  write(lvltag, '(i2)') k
        if(k >= 100) write(lvltag, '(i3)') k
        call read_data ( fname, 'fusion_'//trim(lvltag), soil%fusion(:,:,:,k) )
     enddo
     soil%mask = (soil%frac >= 0)   
     is_data_read = .true.
  else if (file_exist('INPUT/soil.res')) then
     if (mpp_pe() == mpp_root_pe()) call error_mesg ('soil_mod', &
            'Reading native formatted restart file.', NOTE)
     unit = open_restart_file('INPUT/soil.res', 'read')
     call set_domain( soil%domain )
     call read_data ( unit, soil%frac )
     call read_data ( unit, soil%temp )
     call read_data ( unit, soil%water )
     call read_data ( unit, soil%snow )
     call read_data ( unit, soil%groundwater )
     call read_data ( unit, soil%fusion )
     call close_file( unit )
     is_data_read = .true. 
     ! set up mask and fractional area
     soil%mask = (soil%frac >= 0)
  else
     ! set up initial land area, fractional area, and mask
     soil%area  = garea(soil%is:soil%ie, soil%js:soil%je)   &
          * gfrac(soil%is:soil%ie, soil%js:soil%je)
     soil%frac  = frac
     soil%mask  = mask
  endif

  ! reconcile fractional area with mask: the sum of fractional areas
  ! has to be 1 everywhere, even where there is no land. 
  where ( .not.soil%mask )             soil%frac        = 0.0
  where ( .not.any(soil%mask, DIM=3) ) soil%frac(:,:,1) = 1.0

  ! initialize accumulated values
  soil%drainage_accum = 0
  soil%calving_accum  = 0

  ! initialize soil diagnostics -- must be before calls to initialization of 
  ! sub-modules because they use id_lon and id_lat, defined by this call
  call init_soil_diag ( id_lon, id_lat, soil%z, soil%zz, Time )

  call land_properties_init ( &
       gblon(soil%is:(soil%ie+1),soil%js:(soil%je+1)), &
       gblat(soil%is:(soil%ie+1),soil%js:(soil%je+1)), &
       soil%mask, time, soil%domain, &
       soil%glacier,         &
       soil%lake,            &
       soil%rough_mom,       &
       soil%rough_heat,      &
       soil%rough_scale,     &
       soil%tcon,            &
       soil%rho_cap,         &
       soil%stomatal,        &
       soil%max_water,       &
       soil%tau_groundwater, &
       soil%max_snow,        &
       soil%cover_type,      &
       id_lon, id_lat )

  if (.not.is_data_read) then
     ! if we did not read the restart file, initialize data with reasonable values
     ! ("cold start"). Note that we could not do that before initialization of land
     ! properties; on the other hand the initialization of land properties cannot be
     ! done before reading restart file since mask is being read from restart.
     Soil%water = init_water

     Soil%snow = init_snow

     Soil%groundwater = init_groundwater
     init_fusion = init_frozen_water*hlf
     do k = 1, n_levels
        Soil%temp(:,:,:,k) = init_temp_no_ice
        Soil%fusion(:,:,:,k) = init_fusion
        where (Soil%glacier .or. Soil%snow>0)
           Soil%temp(:,:,:,k) = init_temp_ice
        endwhere
     enddo
  endif

#ifdef SCM
  !--- for single column model -------------------------------------!
  !--- initialize surface temperature to observed value ------------!  
  if (do_specified_tskin .and. do_specified_land) then     
       Soil%temp(:,:,:,1) = TSKIN
  end if   
  !-----------------------------------------------------------------!
#endif

  ! initialize rivers submodule
  call rivers_init &
     ( Soil%rivers, gblon, gblat, garea, gfrac, time, dt_fast, dt_slow, domain, &
       id_lon, id_lat )

  !  assign initial values to dynamic land properties (albedo)
!! either include all the albedo components in this call, or define 
!!  the albedo components upon return. 
  call update_land_properties_fast &
       (soil%snow(:,:,:), soil%temp(:,:,:,1), soil%mask, soil%albedo(:,:,:), &
        soil%cover_type(:,:,:))

!      (soil%snow(:,:,:), soil%temp(:,:,:,1), soil%mask,   &
!        soil%albedo(:,:,:),  &
!        soil%albedo_vis_dir(:,:,:),  &
!        soil%albedo_nir_dir(:,:,:),  &
!        soil%albedo_vis_dif(:,:,:),  &
!        soil%albedo_nir_dif(:,:,:) )
!! FOR NOW, set all values to be the same:
   soil%albedo_vis_dir = soil%albedo
   soil%albedo_nir_dir = soil%albedo
   soil%albedo_vis_dif = soil%albedo
   soil%albedo_nir_dif = soil%albedo

  ! send static fields to diagnostic manager for output
  call diag_static(soil, gfrac(soil%is:soil%ie, soil%js:soil%je))

  ! copy glacier mass conservation flag to the soil data stricture
  soil%conserve_glacier_mass = conserve_glacier_mass
  soil%conserve_glacier_snow_mass = conserve_glacier_snow_mass

end subroutine soil_init
! </SUBROUTINE>


! <SUBROUTINE NAME="init_soil_diag">

!   <OVERVIEW>
!     Initializes diagnostic for soil.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Initializes diagnostics for soil. Defines vertical axis, array of
!     axis indices and diagnostic and static fields.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call init_soil_diag ( id_lon, id_lat, zfull, zhalf, Time )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine init_soil_diag ( id_lon, id_lat, zfull, zhalf, Time )

  integer,         intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer,         intent(in) :: id_lat  ! ID of land longitude (X) axis
  real,            intent(in) :: zfull(:)! Full levels, m
  real,            intent(in) :: zhalf(:)! Half levels, m
  type(time_type), intent(in) :: Time    ! Current time
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: axes(3)

  ! define vertical axis
  id_zhalf = diag_axis_init ( &
       'zhalf', zhalf, 'meters', 'z', 'half level',  -1, set_name='land' )
  id_zfull = diag_axis_init ( &
       'zfull', zfull, 'meters', 'z', 'full level',  -1, set_name='land', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define diagnostic fields
  id_water = register_diag_field ( module_name, 'water', axes(1:2),  &
       Time, 'mass of water in bucket', 'kg/m2', missing_value=-100.0 )
  id_snow  = register_diag_field ( module_name, 'snow',  axes(1:2),  &
       Time, 'mass of snow on ground', 'kg/m2',  missing_value=-100.0 )
  id_temp  = register_diag_field ( module_name, 'temp',  axes,       &
       Time, 'temperature',            'degK',  missing_value=-100.0 )
  id_frozen = register_diag_field ( module_name, 'frozen',  axes,       &
       Time, 'amount of frozen water', 'kg/m3',  missing_value=-100.0 )
  id_albedo = register_diag_field ( module_name, 'albedo',  axes(1:2), &
       Time, 'albedo', 'dimensionless',  missing_value=-100.0 )
  id_albedo_vis_dir = register_diag_field ( module_name, 'albedo_vis_dir',  axes(1:2), &
       Time, 'albedo_vis_dir', 'dimensionless',  missing_value=-100.0 )
  id_albedo_vis_dif = register_diag_field ( module_name, 'albedo_vis_dif',  axes(1:2), &
       Time, 'albedo_vis_dif', 'dimensionless',  missing_value=-100.0 )
  id_albedo_nir_dir = register_diag_field ( module_name, 'albedo_nir_dir',  axes(1:2), &
       Time, 'albedo_nir_dir', 'dimensionless',  missing_value=-100.0 )
  id_albedo_nir_dif = register_diag_field ( module_name, 'albedo_nir_dif',  axes(1:2), &
       Time, 'albedo_nir_dif', 'dimensionless',  missing_value=-100.0 )
  id_drainage = register_diag_field ( module_name, 'drainage',  axes(1:2), &
       Time, 'drainage rate',          'kg/(m2 s)',missing_value=-100.0 )
  id_calving = register_diag_field ( module_name, 'calving',  axes(1:2), &
       Time, 'snow sweep rate',     'kg/(m2 s)',missing_value=-100.0 )
  id_precip = register_diag_field ( module_name, 'precip',  axes(1:2), &
       Time, 'total precipitation rate',  'kg/(m2 s)', missing_value=-100.0 )
  id_snowfall = register_diag_field ( module_name, 'snowfall',  axes(1:2), &
       Time, 'snowfall rate',          'kg/(m2 s)', missing_value=-100.0 )
  id_evapor = register_diag_field ( module_name, 'evap',  axes(1:2), &
       Time, 'evaporation rate',       'kg/(m2 s)', missing_value=-100.0 )
  id_sublim = register_diag_field ( module_name, 'sublim',axes(1:2), &
       Time, 'sublimation rate',       'kg/(m2 s)', missing_value=-100.0 )
  id_smelt = register_diag_field ( module_name, 'smelt',  axes(1:2),  &
       Time, 'snow melt rate',         'kg/(m2 s)',     missing_value=-100.0 )
  id_gmelt = register_diag_field ( module_name, 'gmelt',  axes(1:2),  &
       Time, 'glacier melt rate',      'kg/(m2 s)', missing_value=-100.0 )
  id_watsno = register_diag_field ( module_name, 'watsno',  axes(1:2),  &
       Time, 'rate of water theft by snow',    'kg/(m2 s)', missing_value=-100.0 )
  id_sens = register_diag_field ( module_name, 'sens',  axes(1:2),  &
       Time, 'sensible heat flux',      'W/m2',  missing_value=-999.0 )
  id_lhf = register_diag_field ( module_name, 'latent',  axes(1:2),  &
       Time, 'latent heat flux',        'W/m2',  missing_value=-999.0 )
  id_lw   = register_diag_field ( module_name, 'flw',  axes(1:2),      &
       Time, 'net longwave radiative flux', 'W/m2', missing_value=-999.0 )
  id_sw   = register_diag_field ( module_name, 'fsw',  axes(1:2),      &
       Time, 'net shortwave radiative flux', 'W/m2', missing_value=-999.0 )
  id_surface_runoff= register_diag_field ( module_name, 'wroff',  axes(1:2), &
       Time, 'surface runoff of water',    'kg/(m2 s)', missing_value=-999.0 )
  id_surface_runoff_snow= register_diag_field ( module_name, 'sroff', axes(1:2), &
       Time, 'surface runoff of snow',     'kg/(m2 s)', missing_value=-999.0 )
  id_groundwater   = register_diag_field ( module_name, 'groundwater',  axes(1:2),&
       Time, 'mass of water below bucket', 'kg/m2', missing_value=-999.0 )
  id_cover_type = register_diag_field( module_name, 'cover_type', axes(1:2),    &
       Time, 'Land surface cover type', 'dimensionless', missing_value=-1.0)
  id_hlf  = register_diag_field ( module_name, 'hlf', long_name='Latent heat of fusion', units='J/kg' )
  id_hlv  = register_diag_field ( module_name, 'hlv', long_name='Latent heat of evaporation', units='J/kg' )


  ! define static diagnostic fields
  id_area  = register_static_field ( module_name, 'area',  axes(1:2), &
       'land area', 'm2', missing_value=-999.0 )
  id_lfrac = register_static_field ( module_name, 'lfrac', axes(1:2), &
       'fraction of land in cell','none', missing_value=-999.0 ) 
  id_glacier  = register_static_field ( module_name, 'glacier',  axes(1:2), &
       'glacier logical mask', 'dimensionless', missing_value=-999.0 )
  id_tcon     = register_static_field ( module_name, 'tcon',  axes, &
       'ground thermal conductivity', 'W/(m K s)', missing_value=-999.0 )
  id_rho_cap  = register_static_field ( module_name, 'rho_cap',  axes, &
       'ground volumetric heat capacity', 'J/(m3 K)', missing_value=-999.0 )
  id_stomatal = register_static_field ( module_name,'stomatal', axes(1:2), &
       'non-water-stressed stomatal resistance', 's/m', missing_value=-999.0 )
  id_max_water = register_static_field ( module_name, 'max_water', axes(1:2),&
       'root-zone water capacity', 'kg/m2', missing_value=-999.0 )
  id_rough_mom = register_static_field ( module_name, 'rough_mom',axes(1:2),&
       'momentum roughness length', 'm', missing_value=-999.0 )
  id_rough_heat = register_static_field ( module_name, 'rough_heat',axes(1:2),&
       'scalar roughness length', 'm', missing_value=-999.0 )
  id_tau_gw = register_static_field ( module_name, 'tau_gw',axes(1:2),&
       'groundwater residence time', 's', missing_value=-999.0 )
  
end subroutine init_soil_diag
! </SUBROUTINE>


! <SUBROUTINE NAME="soil_end">

!   <OVERVIEW>
!     Deallocates soil data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Deallocates soil data.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call soil_end ( soil )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine soil_end ( soil )

  type(soil_type), intent(inout)      :: soil    ! data to finish using
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer unit ! file unit number to write
  real    frac(soil%is:soil%ie, soil%js:soil%je, soil%n_tiles)
  character(len=64) :: fname = 'RESTART/soil.res.nc'
  character(len=64) :: lvltag
  integer           :: k

  module_is_initialized =.FALSE.
  ! finish using rivers data
  call rivers_end ( soil%rivers )

  ! save restart file
  frac=soil%frac
  where ( .not.soil%mask ) frac = -1

  call set_domain(soil%domain)

  if( do_netcdf_restart ) then
     if (mpp_pe() == mpp_root_pe()) call error_mesg ('soil_mod', &
          'Writing NetCDF formatted restart file: RESTART/soil.res.nc', NOTE)
     call write_data ( fname, 'frac', frac )
     do k = 1, size(soil%temp, 4)
        if(k < 10)   write(lvltag, '(i1)') k
        if(k >= 10)  write(lvltag, '(i2)') k
        if(k >= 100) write(lvltag, '(i3)') k
        call write_data ( fname, 'soil_temp_'//trim(lvltag), soil%temp(:,:,:,k) )
     enddo
     call write_data ( fname, 'water', soil%water )
     call write_data ( fname, 'snow', soil%snow )
     call write_data ( fname, 'groundwater', soil%groundwater )
     do k = 1, size(soil%fusion, 4)
        if(k < 10)   write(lvltag, '(i1)') k
        if(k >= 10)  write(lvltag, '(i2)') k
        if(k >= 100) write(lvltag, '(i3)') k
        call write_data ( fname, 'fusion_'//trim(lvltag), soil%fusion(:,:,:,k) )
     enddo
  else
     if (mpp_pe() == mpp_root_pe()) call error_mesg ('soil_mod', &
          'Writing native formatted restart file.', NOTE)
  unit = open_restart_file ('RESTART/soil.res', 'write')

  call write_data ( unit,frac )
  call write_data ( unit,soil%temp )
  call write_data ( unit,soil%water )
  call write_data ( unit,soil%snow )
  call write_data ( unit,soil%groundwater )
  call write_data ( unit,soil%fusion )

  call close_file ( unit )
  endif

  ! deallocate data
  deallocate ( &
       soil % dz,   &
       soil % z,    &
       soil % zz,   &
       soil % aa,   &
       soil % cc,   &
       soil % max_fusion,   &

       soil % area, & 
       soil % frac, & 
       soil % mask, & 

       soil % temp, &

       soil % snow, &
       soil % water, &
       soil % groundwater, &
       soil % drainage, &
       soil % drainage_accum, &
       soil % calving, &
       soil % calving_accum, &
       soil % albedo, &
       soil % albedo_vis_dir, &
       soil % albedo_nir_dir, &
       soil % albedo_vis_dif, &
       soil % albedo_nir_dif, &
       soil % rough_mom, &
       soil % rough_heat, &
       soil % rough_scale, &
       soil % stomatal, &
       soil % max_water, & 
       soil % tcon, & 
       soil % rho_cap, & 
       soil % tau_groundwater, &

       soil % glacier, &
       soil % lake, &

       soil % cover_type, &
       soil % ground, &
       soil % fusion  )
  
end subroutine soil_end
! </SUBROUTINE>


! <SUBROUTINE NAME="update_soil_fast">

!   <OVERVIEW>
!     Fast time-scale soil updater.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Calculate surface temperature with implicit correction via upward part
!     of tridiagonal elimination for implicit correction. The surface
!     temperature is corrected for snow melt before proceeding to complete
!     downward part of tridiagonal elimination. If the temperature is above
!     freezing and snow or glacial ice is present, compute potential melt.
!
!     In non-glaciated regions, if there is not enough snow, melt only what's
!     there and increase surface temperature above freezing. In glaciated
!     regions, if there is not enough snow, melt only what's there and then
!     melt the glacier. Otherwise go ahead and use all of the potential snow
!     melt, and leave the surface temperature at freezing.
!
!     Accumulate snow melt water into bucket in non-glacier regions. Using
!     this snow melt-modified surface temperature, correct all fluxes and
!     compute temperature increment and new temperature for output. Complete
!     the tridiagonal solver.
!     
!     Bucket hydrology involves treating unglaciated cells first. Tentative
!     changes are based on time-step evaporation and precipitation. If snow
!     cover went negative during this step, take the necessary mass from the
!     soil (water) and re-zero snow. Then put excess soil water into soil
!     drainage. Then do the glaciated cells. soil water is not present, and
!     glacier mass is not tracked. For all land cells, put excess snow into
!     calving. Accumulate drainage for transfer to groundwater on the slow
!     time step. Similarly, accumulate calving. Compute the new albedos and 
!     update rivers.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_soil_fast ( soil, &
!     fsw, flw, sens, evap, dhdt, dedt, drdt, lprec, fprec )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_soil_fast ( soil, &
     fsw, flw, sens, evap, dhdt, dedt, drdt, lprec, fprec )

  type(soil_type), intent(inout) :: soil ! soil state to update
  real, intent(in) :: &
       fsw  (soil%is:soil%ie,soil%js:soil%je,soil%n_tiles), & ! shortwave flux
       flw  (soil%is:soil%ie,soil%js:soil%je,soil%n_tiles), & ! longwave flux
       sens (soil%is:soil%ie,soil%js:soil%je,soil%n_tiles), & ! sensible heat
                                                              ! flux
       evap (soil%is:soil%ie,soil%js:soil%je,soil%n_tiles), & ! evaporation
       dhdt (:,:,:), & ! derivative of sens over T
       dedt (:,:,:), & ! derivative of evap over T
       drdt (:,:,:), & ! derivative of LW radiation over T
       lprec(:,:,:), & ! liquid prec
       fprec(:,:,:)    ! solid prec (snow)
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  real, dimension(size(soil%frac,1), size(soil%frac,2), size(soil%frac,3)) :: &
       t_surf,        &
       snow_melt,     &
       glacier_melt,  &
       watsno,        &
       rho_cap_dz,    &
       latent,        & ! spec heat of sublimation/evaporation
       sublim,        &
       t_to_snow,     & 
       phase_change,  &
       sens_new,      &
       evap_new,      &
       flw_new,       &
       aaa, ccc, bbb, &
       num, denom,    &
       flux, deriv
  logical, dimension(size(fsw,1), size(fsw,2), size(soil%frac,3)) :: &
       not_glacier,   &
       ice_or_snow,   &
       snow_on_ground
  real, dimension(&
       size(fsw,1), size(fsw,2), size(soil%temp,3), size(soil%temp,4))::&
       e,             &
       f,             &
       del_t

  integer :: k, k2
  integer :: num_lev        ! number of levels in the soil
  

  ! check input argument consistency -- should be disabled after debugging
  if ( size(fsw,1) /= size(soil%area,1) .or. size(fsw,2) /= size(soil%area,2) )  &
       call error_mesg ( module_name, &
       'update_land_fast:: input arguments have incorrect horizontal dimensions', &
       FATAL)

!   <ERROR MSG="update_land_fast:: input arguments have incorrect horizontal dimensions" STATUS="FATAL">
!   </ERROR>

  not_glacier  = soil%mask .and. .not.soil%glacier
  ice_or_snow  = soil%mask .and. (soil%snow > 0 .or. soil%glacier)

  ! when snow or ice is present, include latent heat of fusion when computing 
  ! energy equivalent of evaporation (sublimation) rate 

  latent = hlv
  where (ice_or_snow) latent = hlf + hlv

  ! some initial values
  t_surf       = soil%temp(:,:,:,1)
  snow_melt    = 0.0
  glacier_melt = 0.0
  evap_new     = 0.0
  sens_new     = 0.0
  flw_new      = 0.0
  sublim       = 0.0
  watsno       = 0.0

  num_lev = size(soil%temp,4)

  ! some useful combinations
  where(soil%mask) 
     ! factor for converting flux into temperature change in top layer
     rho_cap_dz = soil%dt/(soil%rho_cap(:,:,:,1)*soil%dz(1))
     ! factor for converting temperature change into energy equivalent 
     ! of snow melt
     t_to_snow  = soil%rho_cap(:,:,:,1)*soil%dz(1)/hlf
     ! temperature change in top layer due to flux into soil
     flux  = (flw + fsw - (sens + latent*evap))*rho_cap_dz  
     ! derivative  of 'flux' with respect to surface T -- a non-dimensional quantity
     deriv = -(dhdt + latent*dedt + drdt)*rho_cap_dz 
  endwhere

  if(num_lev > 1) then
     ! explicit tendency due to diffusion
     

! pcm NOTE: I REALIZE THERE ARE MORE EFFICIENT WAYS TO PROGRAM THIS !!!
     
     do k = 2, num_lev - 1
        where(soil%mask) &
             del_t(:,:,:,k) =  (    &
!            - soil%diff*(  soil%aa(k)*(soil%temp(:,:,:,k+1) - soil%temp(:,:,:,k))  &
!                         - soil%cc(k)*(soil%temp(:,:,:,k)   - soil%temp(:,:,:,k-1)))
             - (  (dz(k)+dz(k+1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k+1)/soil%tcon(:,:,:,k+1))) &
             *(  soil%aa(k)*(soil%temp(:,:,:,k+1) - soil%temp(:,:,:,k)))  &
             + (  (dz(k)+dz(k-1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k-1)/soil%tcon(:,:,:,k-1))) &
             *(  soil%cc(k)*(soil%temp(:,:,:,k)   - soil%temp(:,:,:,k-1))) &
                               ) / soil%rho_cap(:,:,:,k)

     enddo

     k = 1
     k2 = num_lev
     where ( soil%mask ) 
!       aaa = soil%diff*soil%aa(1)
!       ccc = soil%diff*soil%cc(num_lev)
        
        aaa = (  (dz(k)+dz(k+1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k+1)/soil%tcon(:,:,:,k+1))) &
                       *soil%aa(k) / soil%rho_cap(:,:,:,k)
        ccc = (  (dz(k2)+dz(k2-1))/(dz(k2)/soil%tcon(:,:,:,k2)+dz(k2-1)/soil%tcon(:,:,:,k2-1))) &
                       *soil%cc(k2) / soil%rho_cap(:,:,:,k2)
        
        del_t(:,:,:,1)= &
             -aaa*(soil%temp(:,:,:,2) - soil%temp(:,:,:,1))
        del_t(:,:,:,num_lev) =  &
             ccc*(soil%temp(:,:,:,num_lev) - soil%temp(:,:,:,num_lev-1))

        e(:,:,:,num_lev) = 0.0   ! boundary condition for tridiagonal elimination
        f(:,:,:,num_lev) = 0.0   ! boundary condition for tridiagonal elimination
     endwhere

     ! upward part of tridiagonal elimination for implicit correction
     !   (one could save a few multiplications by defining some additional arrays)

     k = num_lev
        where(soil%mask) 
!          aaa = soil%diff*soil%aa(k)
!          ccc = soil%diff*soil%cc(k)
           aaa =  0
           ccc =  ( (dz(k)+dz(k-1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k-1)/soil%tcon(:,:,:,k-1)))&
                   * soil%cc(k)/soil%rho_cap(:,:,:,k)
           bbb = 1.0 - aaa - ccc
           denom = aaa*e(:,:,:,k) + bbb
           e(:,:,:,k-1) = -ccc/denom
           f(:,:,:,k-1) = (del_t(:,:,:,k) - aaa*f(:,:,:,k))/denom
        endwhere

     do k = num_lev-1, 2, -1
        where(soil%mask) 
!          aaa = soil%diff*soil%aa(k)
!          ccc = soil%diff*soil%cc(k)
           aaa =  ( (dz(k)+dz(k+1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k+1)/soil%tcon(:,:,:,k+1)))&
                   * soil%aa(k)/soil%rho_cap(:,:,:,k)
           ccc =  ( (dz(k)+dz(k-1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k-1)/soil%tcon(:,:,:,k-1)))&
                   * soil%cc(k)/soil%rho_cap(:,:,:,k)
           bbb = 1.0 - aaa - ccc
           denom = aaa*e(:,:,:,k) + bbb
           e(:,:,:,k-1) = -ccc/denom
           f(:,:,:,k-1) = (del_t(:,:,:,k) - aaa*f(:,:,:,k))/denom
        endwhere
     end do

     ! surface temperature with implicit correction

     k = 1
     where(soil%mask)
        aaa   =  ( (dz(k)+dz(k+1))/(dz(k)/soil%tcon(:,:,:,k)+dz(k+1)/soil%tcon(:,:,:,k+1)))&
                  *soil%aa(1)/soil%rho_cap(:,:,:,1)
        num   = del_t(:,:,:,1) - aaa*f(:,:,:,1)  + flux
        denom = 1.0 + aaa*(e(:,:,:,1) -1.0) - deriv
        
        t_surf = soil%temp(:,:,:,1) + num/denom
     endwhere

  else  ! one level case
     where(soil%mask)
        num   = flux
        denom = 1.0 - deriv
        t_surf = soil%temp(:,:,:,1) + num/denom
     endwhere
  end if

  ! correct this surface temperature for snow melt before proceeding to 
  !    complete "downward" part of tridiagonal elimination

  ! if temperature is above freezing and snow or glacial ice is present, 
  !     compute potential melt

  where (ice_or_snow .and. t_surf > tfreeze) &
       snow_melt = (t_surf - tfreeze)*denom*t_to_snow


  ! in non-glaciated region, 
  !    if there is not enough snow, melt only what's there and increase 
  !    surface temperature above freezing

  where (not_glacier .and. snow_melt > soil%snow ) 
     snow_melt  = soil%snow
     t_surf     = soil%temp(:,:,:,1) + (num - snow_melt/t_to_snow)/denom
     soil%snow  = 0.0
     soil%water = soil%water + snow_melt
  endwhere 

  ! in glaciated region, 
  !    if there is not enough snow, melt only what's there and then
  !     melt the glacier

  where (soil%glacier .and. snow_melt > soil%snow) 
     t_surf       = tfreeze
     glacier_melt = snow_melt - soil%snow
     snow_melt    = soil%snow
     soil%snow    = 0.0
  endwhere

  ! otherwise go ahead and use all of the potential snow melt, and leave the 
  !   surface temperature at freezing

  ! accumulate snow melt water into bucket in non-glacier regions

  where (not_glacier .and. snow_melt > 0.0 .and. snow_melt <= soil%snow)
     soil%water = soil%water + snow_melt
  endwhere

  where (soil%mask .and. snow_melt > 0 .and. snow_melt <= soil%snow)
     t_surf    = tfreeze
     soil%snow = soil%snow - snow_melt
  endwhere

  !  using this snow melt-modified surface temperature, correct all fluxes
  !   and compute temperature increment and new temperature for output 
  where (soil%mask)
     del_t(:,:,:,1)     = t_surf - soil%temp(:,:,:,1)
     soil%temp(:,:,:,1) = t_surf
     sens_new = sens + dhdt*del_t(:,:,:,1)
     evap_new = evap + dedt*del_t(:,:,:,1)
     flw_new  = flw  - drdt*del_t(:,:,:,1)
  endwhere

  where (ice_or_snow) sublim = evap_new

  ! completing the tridiagonal solver
  if( num_lev > 1) then
     do k = 1, num_lev-1
        where (soil%mask)
           del_t(:,:,:,k+1)     = e(:,:,:,k)*del_t(:,:,:,k) + f(:,:,:,k)
           soil%temp(:,:,:,k+1) = soil%temp(:,:,:,k+1) + del_t(:,:,:,k+1)
        endwhere
     end do
  end if

  do k=1, num_lev
     where (soil%mask .and. soil%fusion(:,:,:,k)>0.and.soil%temp(:,:,:,k)>tfreeze)
        phase_change = min (soil%fusion(:,:,:,k),(soil%temp(:,:,:,k)-tfreeze)*soil%rho_cap(:,:,:,k))
        soil%fusion(:,:,:,k) = soil%fusion(:,:,:,k) - phase_change
        soil%temp(:,:,:,k) = soil%temp(:,:,:,k) - phase_change/soil%rho_cap(:,:,:,k)
     elsewhere (soil%mask .and.soil%fusion(:,:,:,k)<soil%max_fusion(k) .and. soil%temp(:,:,:,k)<tfreeze) 
        phase_change = min (soil%max_fusion(k)-soil%fusion(:,:,:,k),(tfreeze-soil%temp(:,:,:,k))*soil%rho_cap(:,:,:,k))
        soil%fusion(:,:,:,k) = soil%fusion(:,:,:,k) + phase_change
        soil%temp(:,:,:,k) = soil%temp(:,:,:,k) + phase_change/soil%rho_cap(:,:,:,k)
     endwhere
  enddo

  !--- bucket hydrology  (snow_melt has already been taken care of) ----
  snow_on_ground = soil%snow > 0.0

  !--- first treat unglaciated cells...
  !
  ! tentative changes based on time-step evaporation and precipitation.
  ! if snow cover went negative during step, take necessary mass from
  ! soil (water) and re-zero snow. (this does not conserve latent heat?)
  ! then put excess soil water into soil drainage.
  !
  where(not_glacier .and. .not. snow_on_ground) 
     soil%water = soil%water + (lprec - evap_new)*soil%dt
     soil%snow  = fprec*soil%dt
     where(ice_or_snow) watsno = evap_new
  end where
  where (not_glacier .and. snow_on_ground)  
     soil%water = soil%water  + lprec*soil%dt 
     soil%snow  = soil%snow   + (fprec - evap_new)*soil%dt
     soil%water = soil%water  + min(soil%snow, 0.0) 
     watsno     = -min(soil%snow, 0.0)/soil%dt
     soil%snow  = max(soil%snow, 0.0)
  end where

  where (not_glacier) 
     soil%drainage = max((soil%water-soil%max_water)/soil%dt, 0.0)
     soil%water    = min( soil%water,soil%max_water)
  end where
  where (soil%lake)
     soil%drainage = soil%drainage + (soil%water - soil%max_water)/soil%dt
     soil%water= soil%max_water
  end where

  !--- then do the glaciated cells. soil water is not present, and
  !    glacier mass is not tracked.
  if (soil%conserve_glacier_snow_mass) then
      where(soil%glacier) soil%snow = soil%snow  + (fprec - evap_new)*soil%dt
    else
      where(soil%glacier .and. .not. snow_on_ground) &
           soil%snow  = fprec*soil%dt  ! seems to asssume evap_new is not<0 !!
      where (soil%glacier .and. snow_on_ground)  
         soil%snow = soil%snow  + (fprec - evap_new)*soil%dt
         soil%snow = max(soil%snow, 0.0)
      end where
    endif
  if (soil%conserve_glacier_mass) then
      where (soil%glacier) 
        soil%drainage = snow_melt/soil%dt + lprec
        soil%water    = soil%max_water ! to avoid spurious "bone dry" on glaciers
      endwhere
    else
      where (soil%glacier) 
        soil%drainage = (snow_melt + glacier_melt)/soil%dt + lprec
        soil%water    = soil%max_water ! to avoid spurious "bone dry" on glaciers
      endwhere
    endif

  !--- for all land cells, put excess snow into calving.
  where (soil%mask)
     soil%calving = max((soil%snow - soil%max_snow)/soil%dt, 0.0)
     soil%snow    = min(soil%snow, soil%max_snow)
  endwhere

  !--- accumulate drainage for transfer to groundwater on slow time step.
  !    similarly, accumulate calving
  where (soil%mask)
     soil%drainage_accum = soil%drainage_accum &
          + soil%drainage*soil%dt
     soil%calving_accum  = soil%calving_accum &
          + soil%calving*soil%dt
  endwhere

  !--- end of bucket hydrology -----------------------------------------

  ! increment time
  soil%Time = increment_time(soil%Time, int(soil%dt), 0)

  !--- override snow and water, if necessary ---------------------------
  call data_override('LND', 'water', soil%water, soil%Time)
  call data_override('LND', 'snow',  soil%snow,  soil%Time)

#ifdef SCM
  !--- for single column model keep top level --------------------------
  !--- temperature at observed value -----------------------------------
  if (do_specified_tskin .and. do_specified_land) then     
     soil%temp(:,:,:,1) = TSKIN
  end if
#endif

  !  compute new albedos
!! either include all the albedo components in this call, or define
!! the albedo components upon return.
  call update_land_properties_fast &
       (soil%snow, soil%temp(:,:,:,1), soil%mask, soil%albedo, soil%cover_type)
!      (soil%snow(:,:,:), soil%temp(:,:,:,1), soil%mask,   &
!        soil%albedo(:,:,:),  &
!        soil%albedo_vis_dir(:,:,:),  &
!        soil%albedo_nir_dir(:,:,:),  &
!        soil%albedo_vis_dif(:,:,:),  &
!        soil%albedo_nir_dif(:,:,:) ) 

!! FOR NOW, set all values to be the same:
     soil%albedo_vis_dir = soil%albedo
     soil%albedo_nir_dir = soil%albedo
     soil%albedo_vis_dif = soil%albedo
     soil%albedo_nir_dif = soil%albedo        

  ! update rivers
  call update_rivers_fast ( soil%rivers )

  call diag_fast(soil, lprec+fprec, fprec, evap_new, sublim, &
       snow_melt/soil%dt, &
       glacier_melt/soil%dt, watsno, sens_new, evap_new*latent, flw_new, fsw )

end subroutine update_soil_fast
! </SUBROUTINE>


! <SUBROUTINE NAME="update_soil_slow">

!   <OVERVIEW>
!     Slow time-scale soil state updater: groundwater dynamics and routing of 
!     river flows (and calving). 
!   </OVERVIEW>

!   <DESCRIPTION>
!     Slow time-scale soil state updater: groundwater dynamics and routing of 
!     river flows (and calving). Calculations include snow runoff, liquid water
!     runoff and surface runoff. Ground water is modified and the accumulated
!     values for the next interval are cleaned up. Calls to average_tiles and
!     update_rivers_slow.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_soil_slow(soil)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_soil_slow(soil, blon, blat)

  type(soil_type), intent(inout) :: soil  ! state to update
  real, intent(inout) :: blon(:,:)
  real, intent(inout) :: blat(:,:)
!   </PUBLICROUTINE>

  ! ATTENTION: check and recheck to make sure area measures are consistent!

  ! ---- local vars ----------------------------------------------------------
  real, dimension(soil%is:soil%ie,soil%js:soil%je,size(soil%frac,3)):: &
       c0, c1, c2, & ! runoff calculation coefficients
       x             ! temporary value for calculations
  real :: runoff_w (soil%is:soil%ie,soil%js:soil%je) ! liquid runoff
  real :: runoff_s (soil%is:soil%ie,soil%js:soil%je) ! snow runoff
  integer :: model_yr, model_month, model_day
  integer :: model_hour, model_minute, model_second
  integer, save :: this_year = -1

  ! calculate "snow runoff"
  call average_tiles ( soil%calving_accum/soil%dt_slow, soil%frac, soil%mask, runoff_s )

  ! calculate liquid water runoff
  where (soil%mask)
     ! calculate surface runoff
     c0 = soil%dt_slow/soil%tau_groundwater
     c1 = exp(-c0)
     c2 = (1-c1)/c0
     x  = ((1-c1)*soil%groundwater + (1-c2)*soil%drainage_accum) / soil%dt_slow
     
     ! modify ground water
     soil%groundwater    = c1 * soil%groundwater + c2 * soil%drainage_accum

     ! clean up accumulated values for the next interval
     soil%drainage_accum = 0.
     soil%calving_accum  = 0.
  endwhere
  
  call average_tiles ( x, soil%frac, soil%mask, runoff_w)

! get the current year of the model
  call get_date(soil%Time, model_yr, model_month, model_day, model_hour, model_minute, model_second)

! initialize this_year
  if (this_year == -1) this_year = model_yr

! read the cover type forcing data every year
  if (model_yr /= this_year) then
     this_year = model_yr 
     call update_land_properties_slow ( &
         blon, blat, &
         soil%mask, soil%Time, soil%domain, &
         soil%glacier,         &
         soil%lake,            &
         soil%rough_mom,       &
         soil%rough_heat,      &
         soil%rough_scale,     &
         soil%tcon,            & 
         soil%rho_cap,         &
         soil%stomatal,        &
         soil%max_water,       &
         soil%tau_groundwater, &
         soil%max_snow,        &
         soil%cover_type )
  endif 

  call update_rivers_slow ( soil % Rivers, runoff_w, runoff_s )

  ! diagnostic section
  call diag_slow ( soil, runoff_w, runoff_s )
  
end subroutine update_soil_slow
! </SUBROUTINE>


! <SUBROUTINE NAME="update_soil_bnd_fast">

!   <OVERVIEW>
!     Updates boundary data for the atmosphere on the fast time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates boundary data for the atmosphere on the fast time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_soil_bnd_fast ( soil, bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_soil_bnd_fast ( soil, bnd )

  type(soil_type),      intent(in)  :: soil ! current soil state
  type(land_data_type), intent(inout) :: bnd  ! output boundary data
!   </PUBLICROUTINE>

  bnd % t_surf     = soil % temp(:,:,:,1)
  bnd % t_ca       = soil % temp(:,:,:,1) ! vegetation can override this,
                                          ! if necessary
  bnd % albedo     = soil % albedo
  bnd % albedo_vis_dir = soil%albedo_vis_dir
  bnd % albedo_nir_dir = soil%albedo_nir_dir
  bnd % albedo_vis_dif = soil%albedo_vis_dif
  bnd % albedo_nir_dif = soil%albedo_nir_dif
  

end subroutine update_soil_bnd_fast
! </SUBROUTINE>


! <SUBROUTINE NAME="update_soil_bnd_slow">

!   <OVERVIEW>
!     Updates boundary data for the atmosphere on the slow time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates boundary data for the atmosphere on the slow time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_soil_bnd_slow ( soil, bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_soil_bnd_slow ( soil, bnd )

  type(soil_type),      intent(in)  :: soil ! current soil state
  type(land_data_type), intent(inout) :: bnd  ! output boundary data
!   </PUBLICROUTINE>

  bnd % tile_size  = soil % frac
  bnd % mask       = soil % mask

  bnd % rough_mom  = soil % rough_mom
  bnd % rough_heat = soil % rough_heat
  bnd % rough_scale= soil % rough_scale

  call update_rivers_bnd_slow (soil%Rivers, bnd)

end subroutine update_soil_bnd_slow
! </SUBROUTINE>

! <SUBROUTINE NAME="diag_fast">

!   <OVERVIEW>
!     Calculate and return total amount of requested quantitiy per PE.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Compute and return stock values of mass, heat, etc.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call soil_stock_pe &
!     ( soil, index, value )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine soil_stock_pe(soil,index,value)
  type(soil_type), intent(in) :: soil      ! soil state
  integer        , intent(in)  :: index ! ID of the stock to calculate
  real           , intent(out) :: value ! calculated value of the stock

  ! ---- local vars
  integer :: i,j,k,l,num_lev

  select case(index)
  case(ISTOCK_WATER)
     value = 0
     do k = 1, soil%n_tiles
     do j = soil%js,soil%je
     do i = soil%is,soil%ie
        if(.not.soil%mask(i,j,k)) cycle
        value = value + soil%area(i,j)*soil%frac(i,j,k)*&
             (soil%water(i,j,k)+soil%groundwater(i,j,k)+soil%snow(i,j,k))
     enddo
     enddo
     enddo
  case(ISTOCK_HEAT)
    value = 0
    num_lev = size(soil%temp,4)
    do k = 1, soil%n_tiles
    do j = soil%js,soil%je
    do i = soil%is,soil%ie
       if(.not.soil%mask(i,j,k)) cycle
       do l = 1, num_lev
          value = value + soil%area(i,j)*soil%frac(i,j,k)*soil%dz(l)* &
            (-soil%fusion(i,j,k,l)&
             +soil%rho_cap(i,j,k,l)*(soil%temp(i,j,k,l)-tfreeze) )
       enddo
       value = value - hlf*soil%area(i,j)*soil%frac(i,j,k)*soil%snow(i,j,k)
    enddo
    enddo
    enddo
  case default 
    value = 0
  end select

end subroutine soil_stock_pe
! </SUBROUTINE>


! <SUBROUTINE NAME="diag_fast">

!   <OVERVIEW>
!     Sends data for diagnostics on fast time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Sends data for diagnostics on fast time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call diag_fast &
!     ( soil, prec, fprec, evap, sublim, smelt, gmelt, sens, lhf, flw, fsw )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine diag_fast &
     ( soil, prec, fprec, evap, sublim, smelt, gmelt, watsno, sens, lhf, flw, fsw )

  type(soil_type), intent(in) :: soil      ! soil state
  real,            intent(in), dimension(:,:,:) :: &
       prec,  &   ! total precipitation, 
       fprec, &   ! frozen precipitation (snow)
       evap,  &   ! evaporation
       sublim,&   ! sublimation
       smelt, &   ! snow melt
       gmelt, &   ! glacier melt
       watsno,&   ! water stolen by snow to satisfy sublimation
       sens,  &   ! sensible heat flux, W/m2
       lhf,   &   ! latent heat flux, W/m2
       flw,   &   ! net longwave flux, W/m2
       fsw        ! shortwave flux, W/m2
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  logical, pointer, save :: mask2 (:,:) =>NULL()
  logical :: used


  ! set up masks
  mask2 => soil%mask(:,:,1)


  ! send data to diagnostics manager ---

  ! tile-averaged data
  if ( id_water > 0 )    used = send_averaged_data &
       ( id_water,    soil%water,      soil%frac, soil%Time, mask=soil%mask )
  if ( id_snow  > 0 )    used = send_averaged_data &
       ( id_snow ,    soil%snow,       soil%frac, soil%Time, mask=soil%mask )
  if ( id_temp  > 0 )    used = send_averaged_data &
       ( id_temp ,    soil%temp,       soil%frac, soil%Time, mask=soil%mask )
  if ( id_frozen > 0 )    used = send_averaged_data &
       ( id_frozen ,  soil%fusion/hlf, soil%frac, soil%Time, mask=soil%mask )
  if ( id_albedo > 0 )   used = send_averaged_data &
       ( id_albedo,   soil%albedo,     soil%frac, soil%Time, mask=soil%mask )
  if ( id_albedo_vis_dir > 0 )   used = send_averaged_data &
       ( id_albedo_vis_dir,   soil%albedo_vis_dir,     soil%frac, soil%Time, mask=soil%mask )
  if ( id_albedo_nir_dir > 0 )   used = send_averaged_data &
       ( id_albedo_nir_dir,   soil%albedo_nir_dir,     soil%frac, soil%Time, mask=soil%mask )
  if ( id_albedo_vis_dif > 0 )   used = send_averaged_data &
       ( id_albedo_vis_dif,   soil%albedo_vis_dif,     soil%frac, soil%Time, mask=soil%mask )
  if ( id_albedo_nir_dif > 0 )   used = send_averaged_data &
       ( id_albedo_nir_dif,   soil%albedo_nir_dif,     soil%frac, soil%Time, mask=soil%mask )
  if ( id_drainage > 0 ) used = send_averaged_data &
       ( id_drainage, soil%drainage,   soil%frac, soil%Time, mask=soil%mask )
  if ( id_calving > 0 )  used = send_averaged_data &
       ( id_calving,  soil%calving,    soil%frac, soil%Time, mask=soil%mask )
  if ( id_precip > 0 )   used = send_averaged_data &
       ( id_precip,  prec,             soil%frac, soil%Time, mask=soil%mask )
  if ( id_snowfall > 0 ) used = send_averaged_data &
       ( id_snowfall, fprec,           soil%frac, soil%Time, mask=soil%mask )
  if ( id_evapor > 0 )   used = send_averaged_data &
       ( id_evapor,   evap,            soil%frac, soil%Time, mask=soil%mask )
  if ( id_sublim > 0 )   used = send_averaged_data &
       ( id_sublim,   sublim,          soil%frac, soil%Time, mask=soil%mask )
  if ( id_smelt > 0 )    used = send_averaged_data &
       ( id_smelt, smelt,              soil%frac, soil%Time, mask=soil%mask )
  if ( id_gmelt > 0 )    used = send_averaged_data &
       ( id_gmelt,    gmelt,           soil%frac, soil%Time, mask=soil%mask )
  if ( id_watsno > 0 )    used = send_averaged_data &
       ( id_watsno,    watsno,         soil%frac, soil%Time, mask=soil%mask )
  if ( id_sens  > 0 )    used = send_averaged_data &
       ( id_sens ,    sens,            soil%frac, soil%Time, mask=soil%mask )
  if ( id_lhf  > 0 )     used = send_averaged_data &
       ( id_lhf ,     lhf,             soil%frac, soil%Time, mask=soil%mask )
  if ( id_lw    > 0 )    used = send_averaged_data &
       ( id_lw   ,    flw,             soil%frac, soil%Time, mask=soil%mask )
  if ( id_sw    > 0 )    used = send_averaged_data &
       ( id_sw   ,    fsw,             soil%frac, soil%Time, mask=soil%mask )


end subroutine diag_fast
! </SUBROUTINE>


! <SUBROUTINE NAME="diag_slow">

!   <OVERVIEW>
!      Diagnostics on slow time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!      Diagnostics on slow time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call diag_slow (soil, runoff_w, runoff_s )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine diag_slow (soil, runoff_w, runoff_s )

  type(soil_type), intent(in)  :: soil           ! current soil state
  real,            intent(in)  :: runoff_w(:,:)  ! liquid runoff
  real,            intent(in)  :: runoff_s(:,:)  ! snow runoff
!   </PUBLICROUTINE>

  ! --- local data ----------------------------------------------------------
  logical :: mask2(size(soil%frac,1),size(soil%frac,2))
  logical :: used
  real    :: dummy(size(soil%mask,1), size(soil%mask,2))

  mask2 = ANY(soil%mask(:,:,:),3)

  if (id_surface_runoff>0) &
     used = send_data (id_surface_runoff, runoff_w, soil%Time, mask=mask2)
  if (id_surface_runoff_snow>0) &
     used = send_data (id_surface_runoff_snow, runoff_s, soil%Time, mask=mask2)
  if (id_groundwater>0) &
     used = send_averaged_data &
     (id_groundwater, soil%groundwater, soil%frac, soil%Time, mask=soil%mask)
  if (id_cover_type>0) then
     dummy = soil%cover_type(:,:,1)
     used = send_data (id_cover_type, dummy, soil%Time, mask=mask2)
  endif

end subroutine
! </SUBROUTINE>


! <SUBROUTINE NAME="diag_static">

!   <OVERVIEW>
!      Diagnostics of the static variables.
!   </OVERVIEW>

!   <DESCRIPTION>
!      Diagnostics of the static variables.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call diag_static(soil, frac)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine diag_static(soil, frac)

  type(soil_type), intent(in) :: soil            ! soil state to diagnose
  real,            intent(in) :: frac(:,:)       ! fraction of land in each cell
!   </PUBLICROUTINE>

  ! --- local data ----------------------------------------------------------
  real     :: dummy(size(soil%frac,1),size(soil%frac,2))
  logical  :: mask2(size(soil%frac,1),size(soil%frac,2))
  logical  :: mask3(size(soil%frac,1),size(soil%frac,2),size(soil%tcon,4))
  logical  :: used

  ! --- update to averaged data later

  mask2 = ANY(soil%mask(:,:,:),3)
  mask3 = SPREAD(mask2,3,size(mask3,3))

  dummy = 0.
  if ( id_area > 0 ) &
       used = send_data ( id_area, soil%area(:,:), soil%Time, mask=mask2 )
  if ( id_lfrac > 0 ) &
       used = send_data ( id_lfrac, frac(:,:), soil%Time, mask=mask2 )
  where (soil%glacier(:,:,1)) dummy = 1.
  if ( id_glacier > 0 ) &
       used = send_data ( id_glacier, dummy(:,:), soil%Time, mask=mask2 )
  if ( id_tcon > 0 )&
       used = send_data ( id_tcon, soil%tcon(:,:,1,:), soil%Time, mask=mask3 )
  if ( id_rho_cap > 0 )&
       used = send_data ( id_rho_cap, soil%rho_cap(:,:,1,:), soil%Time, &
       mask=mask3 )
  if ( id_stomatal > 0 )&
       used = send_data ( id_stomatal, soil%stomatal(:,:,1), soil%Time, &
       mask=mask2 )
  if ( id_max_water > 0 )&
       used = send_data ( id_max_water, soil%max_water(:,:,1), soil%Time, &
       mask=mask2 )
  if ( id_rough_mom > 0 )&
       used = send_data ( id_rough_mom, soil%rough_mom(:,:,1), soil%Time, &
       mask=mask2 )
  if ( id_rough_heat > 0 )&
       used = send_data ( id_rough_heat, soil%rough_heat(:,:,1), soil%Time, &
       mask=mask2 )
  if ( id_tau_gw > 0 )&
       used = send_data ( id_tau_gw, soil%tau_groundwater(:,:,1), soil%Time, &
       mask=mask2 )
  if ( id_hlf > 0 )&
       used = send_data ( id_hlf, hlf, soil%Time )
  if ( id_hlv > 0 )&
       used = send_data ( id_hlv, hlv, soil%Time )

end subroutine diag_static
!   </SUBROUTINE>



! <FUNCTION NAME="send_averaged_data2d" INTERFACE="send_averaged_data">

!   <OVERVIEW>
!     Averages the data over tiles and sends them to diagnostics. 
!   </OVERVIEW>

!   <DESCRIPTION>
!     Averages the data over tiles and sends them to diagnostics. 
!   </DESCRIPTION>

!   <TEMPLATE>
!     value=send_averaged_data2d ( id, field, area, time, mask )
!   </TEMPLATE>

!   <PUBLICROUTINE INTERFACE="send_averaged_data">
function send_averaged_data2d ( id, field, area, time, mask )

  integer, intent(in)          :: id             ! id of the diagnostic field 
  real,    intent(in)          :: field(:,:,:)   ! field to average and send
  real,    intent(in)          :: area (:,:,:)   ! area of tiles (== averaging 
                                                 ! weights), arbitrary units
  type(time_type), intent(in)  :: time           ! current time
  logical, intent(in),optional :: mask (:,:,:)   ! land mask
!   </PUBLICROUTINE>

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data2d

  ! --- local vars -----------------------------------------------------------
  real  :: out(size(field,1), size(field,2))

  call average_tiles( field, area, mask, out )
  send_averaged_data2d = send_data( id, out, time, mask=ANY(mask,DIM=3) )
end function send_averaged_data2d
!   </FUNCTION>


! <FUNCTION NAME="send_averaged_data3d" INTERFACE="send_averaged_data">

!   <OVERVIEW>
!     Averages the data over tiles and sends them to diagnostics. 
!   </OVERVIEW>

!   <DESCRIPTION>
!     Averages the data over tiles and sends them to diagnostics. 
!   </DESCRIPTION>

!   <TEMPLATE>
!     value=send_averaged_data3d ( id, field, area, time, mask )
!   </TEMPLATE>

!   <PUBLICROUTINE INTERFACE="send_averaged_data">
function send_averaged_data3d( id, field, area, time, mask )

  integer, intent(in)          :: id              ! id of the diagnostic field
  real,    intent(in)          :: field(:,:,:,:)  ! (lon, lat, tile, lev) field 
                                                  ! to average and send
  real,    intent(in)          :: area (:,:,:)    ! (lon, lat, tile) tile areas 
                                                  ! ( == averaging weights), 
                                                  ! arbitrary units
  type(time_type), intent(in)  :: time            ! current time
  logical, intent(in),optional :: mask (:,:,:)    ! (lon, lat, tile) land mask
!   </PUBLICROUTINE>

  ! --- return value ---------------------------------------------------------
  logical  :: send_averaged_data3d     

  ! --- local vars -----------------------------------------------------------
  real    :: out(size(field,1), size(field,2), size(field,4))
  logical :: mask3(size(field,1), size(field,2), size(field,4))
  integer :: it

  do it=1,size(field,4)
     call average_tiles( field(:,:,:,it), area, mask, out(:,:,it) )
  enddo

  mask3(:,:,1) = ANY(mask,DIM=3)
  do it = 2, size(field,4)
     mask3(:,:,it) = mask3(:,:,1)
  enddo

  send_averaged_data3d = send_data( id, out, time, mask=mask3 )
end function send_averaged_data3d
!   </FUNCTION>


! <SUBROUTINE NAME="average_tiles">

!   <OVERVIEW>
!     Average 2-dimensional field over tiles.
!   </OVERVIEW>

!   <DESCRIPTION>
!    Average 2-dimensional field over tiles.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call average_tiles ( x, area, mask, out )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine average_tiles ( x, area, mask, out )

  real,    intent(in)  :: x   (:,:,:) ! (lon, lat, tile) field to average
  real,    intent(in)  :: area(:,:,:) ! (lon, lat, tile) fractional area
  logical, intent(in)  :: mask(:,:,:) ! (lon, lat, tile) land mask
  real,    intent(out) :: out (:,:)   ! (lon, lat)       result of averaging
!   </PUBLICROUTINE>

! --- local vars --------------------------------------------------------------
  integer  :: it                      ! iterator over tile number
  real     :: s(size(x,1),size(x,2))  ! area accumulator

  s(:,:)   = 0.0
  out(:,:) = 0.0

  do it = 1,size(area,3)
     where (mask(:,:,it)) 
        out(:,:) = out(:,:) + x(:,:,it)*area(:,:,it)
        s(:,:)   = s(:,:) + area(:,:,it)
     endwhere
  enddo

  where( s(:,:) > 0 ) &
       out(:,:) = out(:,:)/s(:,:)

end subroutine average_tiles
! </SUBROUTINE>

end module soil_mod
