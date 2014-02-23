#include <fms_platform.h>

module soil_tile_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : &
     write_version_number, file_exist, check_nml_error, &
     close_file, stdlog, read_data, error_mesg, FATAL
use constants_mod, only : &
     pi, tfreeze, rvgas, grav, dens_h2o, hlf
use land_constants_mod, only : &
     NBANDS
use land_io_mod, only : &
     init_cover_field, print_netcdf_error
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_SOIL, register_tile_selector

implicit none
private

! ==== netcdf declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== public interfaces =====================================================
public :: soil_pars_type
public :: soil_prog_type
public :: soil_tile_type

public :: new_soil_tile, delete_soil_tile
public :: soil_tiles_can_be_merged, merge_soil_tiles
public :: soil_is_selected
public :: get_soil_tile_tag
public :: soil_tile_stock_pe
public :: soil_tile_heat, soil_tile_carbon

public :: read_soil_data_namelist
public :: soil_cover_cold_start

public :: soil_data_radiation
public :: soil_data_diffusion
public :: soil_data_thermodynamics
public :: soil_data_hydraulic_properties
public :: soil_data_psi_for_rh
public :: soil_data_gw_hydraulics
public :: soil_data_gw_hydraulics_ar5
public :: soil_data_vwc_for_init_only
public :: soil_data_init_derive_subsurf_pars
public :: soil_data_init_derive_subsurf_pars_ar5
public :: soil_ave_temp  ! calculate average soil temeperature
public :: soil_ave_theta0! calculate average soil moisture, pcm based on available water, zeta input
public :: soil_ave_theta1! calculate average soil moisture, ens based on all water
public :: soil_theta     ! returns array of soil moisture, for all layers
public :: soil_psi_stress ! return soil-water-stress index
public :: g_iso, g_vol, g_geo, g_RT
public :: num_storage_pts
public :: gw_zeta_s, gw_flux_table, gw_area_table
public :: gw_scale_length, gw_scale_relief, gw_scale_soil_depth, slope_exp
public :: num_zeta_pts, num_tau_pts
public :: log_tau, log_zeta_s, log_rho_table, gw_scale_perm, aspect
public :: use_alpha, z_ref, k_macro_constant, use_tau_fix

public :: max_lev, psi_wilt
! =====end of public interfaces ==============================================
interface new_soil_tile
   module procedure soil_tile_ctor
   module procedure soil_tile_copy_ctor
end interface

! ==== module constants ======================================================
character(len=*), parameter   :: &
     version     = '$Id: soil_tile.F90,v 20.0 2013/12/13 23:30:50 fms Exp $', &
     tagname     = '$Name: tikal $', &
     module_name = 'soil_tile_mod'

integer, parameter :: max_lev          = 100 
integer, parameter :: n_dim_soil_types = 9       ! size of lookup table
real,    parameter :: small            = 1.e-4
real,    parameter :: t_ref            = 293
real,    parameter :: g_RT             = grav / (rvgas*t_ref)
real,    parameter :: sigma_max        = 2.2
real,    parameter :: K_rel_min        = 1.e-12

! from the modis brdf/albedo product user's guide:
real            :: g_iso  = 1.
real            :: g_vol  = 0.189184
real            :: g_geo  = -1.377622
real            :: g0_iso = 1.0
real            :: g1_iso = 0.0
real            :: g2_iso = 0.0
real            :: g0_vol = -0.007574
real            :: g1_vol = -0.070987
real            :: g2_vol =  0.307588
real            :: g0_geo = -1.284909
real            :: g1_geo = -0.166314
real            :: g2_geo =  0.041840

! geohydrology option selector (not used by lm2)
integer, parameter, public ::   &
     GW_LM2            = 0, &
     GW_LINEAR         = 1, &
     GW_HILL_AR5       = 2, &
     GW_HILL           = 3, &
     GW_TILED          = 4

! ==== types =================================================================
type :: soil_pars_type
  real vwc_wilt
  real vwc_fc
  real vwc_sat
  real vlc_min
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  real heat_capacity_dry
  real thermal_cond_dry
  real thermal_cond_sat
  real thermal_cond_exp
  real thermal_cond_scale
  real thermal_cond_weight
  real refl_dry_dir(NBANDS)
  real refl_dry_dif(NBANDS)
  real refl_sat_dir(NBANDS)
  real refl_sat_dif(NBANDS)
  real f_iso_dry(NBANDS)
  real f_vol_dry(NBANDS)
  real f_geo_dry(NBANDS)
  real f_iso_sat(NBANDS)
  real f_vol_sat(NBANDS)
  real f_geo_sat(NBANDS)
  real emis_dry
  real emis_sat
  real z0_momentum
  real tau_groundwater
  real rsa_exp         ! riparian source-area exponent
  real hillslope_length
  real hillslope_relief
  real hillslope_zeta_bar
  real soil_e_depth
  real zeta
  real tau
  real k_sat_gw
  real k_sat_sfc
  integer storage_index
  real tfreeze
end type soil_pars_type


type :: soil_prog_type
  real wl
  real ws
  real T
  real groundwater
  real groundwater_T
end type soil_prog_type


type :: soil_tile_type
   integer :: tag ! kind of the soil
   type(soil_pars_type)               :: pars
   type(soil_prog_type), pointer :: prog(:)
   real,                 pointer :: w_fc(:)
   real,                 pointer :: w_wilt(:)
   real,                 pointer :: d_trans(:)
   real,                 pointer :: alpha(:)
   real,                 pointer :: k_macro(:)
   real :: Eg_part_ref
   real :: z0_scalar
   real :: geothermal_heat_flux
   real :: psi_x0 = -1000.
   ! data that were local to soil.f90
   real,                 pointer :: uptake_frac(:)
   real,                 pointer :: heat_capacity_dry(:)
   real,                 pointer :: e(:),f(:)
   ! added to avoid recalculation of soil hydraulics in case of Darcy uptake
   real          :: uptake_T
   real, pointer :: psi(:) ! soil water potential
   real,                 pointer ::  gw_flux_norm(:)
   real,                 pointer ::  gw_area_norm(:)

   ! soil carbon
   real, pointer :: fast_soil_C(:) => NULL() ! fast soil carbon pool, (kg C/m2), per layer
   real, pointer :: slow_soil_C(:) => NULL() ! slow soil carbon pool, (kg C/m2), per layer
   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real, pointer :: asoil_in(:)    => NULL()
   real, pointer :: fsc_in(:)      => NULL()
   real, pointer :: ssc_in(:)      => NULL()
end type soil_tile_type

! ==== module data ===========================================================
integer :: gw_option

real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

!---- namelist ---------------------------------------------------------------
real    :: psi_wilt              = -150.0  ! matric head at wilting
real    :: comp                  = 0.001  ! m^-1, dThdPsi at saturation
real    :: k_over_B              = 2         ! reset to 0 for MCM
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
real    :: z_sfc_layer           = 0.0
real    :: sub_layer_tc_fac      = 1.0
real    :: z_sub_layer_min       = 0.0
real    :: z_sub_layer_max       = 0.0
real    :: freeze_factor         = 1.0
real    :: aspect                = 1.0
real    :: zeta_mult             = 1.0  ! multiplier for root depth scale in stress index
integer :: num_l                 = 18        ! number of soil levels
real    :: dz(max_lev)           = (/ &
    0.02, 0.04, 0.04, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, &
    0.2,   0.4,  0.4,  0.4,  0.4, 0.4,  1.,  1.,  1., &
    0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
                                              ! thickness (m) of model layers,
                                              ! from top down
logical :: use_lm2_awc           = .false.
logical :: lm2                   = .false.
logical :: use_alt_psi_for_rh    = .false.
logical :: use_tau_fix           = .false.
logical :: use_sat_fix           = .false.
logical :: use_comp_for_ic       = .false.
! ---- remainder are used only for cold start ---------
character(32):: soil_to_use     = 'single-tile'
       ! 'multi-tile' for multiple soil types per grid cell, a tile per type
       ! 'single-tile' for geographically varying soil with single type per
       !     model grid cell [default]
       ! 'uniform' for global constant soil, e.g., to reproduce MCM
character(32) :: geohydrology_to_use = 'hill_ar5'
       ! 'lm2' for lm2 hydrology -- must be consistent with lm2 flag
       ! 'linear'
       ! 'hill' -- for hillslope geohydrology
       ! 'hill_ar5' -- hillslope geohydrology reproducing AR5 simulations

logical :: use_mcm_albedo        = .false.   ! .true. for CLIMAP albedo inputs
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
logical :: use_alpha             = .false.   ! for vertical change in soil properties
integer :: soil_index_constant   = 9         ! index of global constant soil,
                                             ! used when use_single_soil
real    :: gw_res_time           = 60.*86400 ! mean groundwater residence time,
                                             ! used when use_single_geo
real    :: rsa_exp_global        = 1.5
real    :: gw_scale_length       = 1.0
real    :: gw_scale_relief       = 1.0
real    :: gw_scale_soil_depth   = 1.0
real    :: slope_exp             = 0.0
real    :: gw_scale_perm         = 1.0
real    :: k_macro_constant      = 0.0
real    :: log_rho_max           = 2.0
real    :: z_ref                 = 0.0       ! depth where [psi/k]_sat = [psi/k]_sat_ref
real    :: geothermal_heat_flux_constant = 0.0  ! true continental average is ~0.065 W/m2
real    :: Dpsi_min_const        = -1.e16

real, dimension(n_dim_soil_types) :: &
  dat_w_sat=&
  (/ 0.380, 0.445, 0.448, 0.412, 0.414, 0.446, 0.424, 0.445, 0.445   /),&
  dat_awc_lm2=&
  (/ 0.063, 0.132, 0.109, 0.098, 0.086, 0.120, 0.101, 0.445, 0.150   /),&
  dat_k_sat_ref=&
  (/ 0.021, .0036, .0018, .0087, .0061, .0026, .0051, .0036, .0036   /),&
  dat_psi_sat_ref=&
  (/ -.059, -0.28, -0.27, -0.13, -0.13, -0.27, -0.16, -0.28, -0.28   /),&
  dat_chb=&
  (/   3.5,   6.4,  11.0,   4.8,   6.3,   8.4,   6.3,   6.4,   6.4   /),&
!  dat_heat_capacity_ref =&
!  (/ 1.8e6, 2.0e6, 2.6e6, 1.9e6, 2.2e6, 2.3e6, 2.1e6, 3.0e6,   1.0   /),&
! previous (ref) values were based on typical water contents
! following dry values are based on w_min=(1-w_sat) w_org=0
! except for peat, where            w_org-(1-w_sat) w_min=0
! microscopic rho*c for w_min is 2650*733 and for w_org is 1300*1926
! (brutsaert 1982 evaporation into the atmosphere p 146)
! ignored air
  dat_heat_capacity_dry =&
  (/ 1.2e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.4e6,   1.0   /),&
!  dat_thermal_cond_ref =&
!  (/   1.5,   0.8,  1.35,  1.15, 1.475, 1.075, 1.217,  0.39, 2.e-7   /),&
! previous (ref) values were based on typical water contents
! following dry and sat values and interpolating exponents are based on
! computations after deVries. i computed C M and F functions for
! unfrozen soil in
! spreadsheet Research\LaD2\soil thermal conductivity. the dry and
! sat values come right out of those computations. the exponent was
! done by eye. curves look like typical literature curves.
! TEMP: still need to treat freezing, maybe import deVries model into code.
  dat_thermal_cond_dry =&
  (/  0.14,  0.21,  0.20,  .175, 0.170, 0.205, 0.183,  0.05, 2.e-7   /),&
  dat_thermal_cond_sat =&
  (/  2.30,  1.50,  1.50,  1.90, 1.900, 1.500, 1.767,  0.50, 2.e-7   /),&
  dat_thermal_cond_scale =&
  (/  15.0,  0.50,   10.,  2.74,  12.2,  2.24,  4.22,   1.0,   1.0   /),&
  dat_thermal_cond_exp =&
  (/   3.0,   5.0,   6.0,   4.0,   4.5,   5.5, 4.667,   1.0,   1.0   /),&
  dat_thermal_cond_weight =&
  (/  0.20,  0.70,   0.7,  0.45, 0.450, 0.700, 0.533,   1.0,   1.0   /),&
  dat_emis_dry=&
  (/ 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,   1.0   /),&
  dat_emis_sat=&
  (/ 0.980, 0.975, 0.970, .9775, 0.975, .9725, 0.975, 0.975,   1.0   /),&
  dat_z0_momentum=&
  (/  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01, 0.045   /),&
  dat_tf_depr=&
  (/  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,   0.0   /)
  !Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
real :: dat_refl_dry_dir(n_dim_soil_types,NBANDS); data dat_refl_dry_dir &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /     ! NIR
real :: dat_refl_dry_dif(n_dim_soil_types,NBANDS); data dat_refl_dry_dif &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /     ! NIR
real :: dat_refl_sat_dir(n_dim_soil_types,NBANDS); data dat_refl_sat_dir &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /     ! NIR
real :: dat_refl_sat_dif(n_dim_soil_types,NBANDS); data dat_refl_sat_dif &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /     ! NIR
  !Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
integer, dimension(n_dim_soil_types) :: &
  input_cover_types=&
  (/ 1,     2,     3,     4,     5,     6,     7,     8,     100   /)
character(len=4), dimension(n_dim_soil_types) :: &
  tile_names=&
  (/'c   ','m   ','f   ','cm  ','cf  ','mf  ','cmf ','peat','mcm ' /)

namelist /soil_data_nml/ psi_wilt, &
     soil_to_use, tile_names, input_cover_types, &
     comp, k_over_B,             &
     rate_fc, sfc_heat_factor, z_sfc_layer, &
     sub_layer_tc_fac, z_sub_layer_min, z_sub_layer_max, freeze_factor, &
     aspect, zeta_mult, num_l,  dz,                      &
     use_lm2_awc,    lm2, use_alt_psi_for_rh, &
     use_tau_fix, use_sat_fix, use_comp_for_ic, use_mcm_albedo,            &
     use_single_geo, geohydrology_to_use, &
     use_alpha, &
     soil_index_constant,         &
     gw_res_time,            rsa_exp_global,      &
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth, &
     slope_exp, &
     gw_scale_perm, z_ref, geothermal_heat_flux_constant, &
     Dpsi_min_const, &
     k_macro_constant, log_rho_max, &
     dat_w_sat,               dat_awc_lm2,     &
     dat_k_sat_ref,            &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_dry,       dat_thermal_cond_dry,   &
     dat_thermal_cond_sat,        dat_thermal_cond_exp,   &
     dat_thermal_cond_scale,        dat_thermal_cond_weight,   &
     dat_refl_dry_dir,            dat_refl_sat_dir,              &
     dat_refl_dry_dif,            dat_refl_sat_dif,              &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,           dat_tf_depr
!---- end of namelist --------------------------------------------------------

real    :: gw_hillslope_length   = 1000.
real    :: gw_hillslope_relief   =  100.
real    :: gw_hillslope_zeta_bar =    0.5
real    :: gw_soil_e_depth       =    4.
real    :: gw_perm               = 3.e-14   ! nominal permeability, m^2
real :: gw_flux_table(26,31), gw_area_table(26,31)
real, allocatable :: log_rho_table(:,:,:)
real, allocatable :: log_deficit_list(:)
real, allocatable :: log_zeta_s(:)
real, allocatable :: log_tau(:)

real, dimension(31 ) :: gw_zeta_s       = &
  (/ 1.0000000e-5, 1.5848932e-5, 2.5118864e-5, 3.9810717e-5, 6.3095737e-5, &
     1.0000000e-4, 1.5848932e-4, 2.5118864e-4, 3.9810717e-4, 6.3095737e-4, &
     1.0000000e-3, 1.5848932e-3, 2.5118864e-3, 3.9810717e-3, 6.3095737e-3, &
     1.0000000e-2, 1.5848932e-2, 2.5118864e-2, 3.9810717e-2, 6.3095737e-2, &
     1.0000000e-1, 1.5848932e-1, 2.5118864e-1, 3.9810717e-1, 6.3095737e-1, &
     1.0000000e+0, 1.5848932e+0, 2.5118864e+0, 3.9810717e+0, 6.3095737e+0, &
     1.0000000e+1 /)

real, dimension(26) :: gw_storage_norm = &
  (/ 0.,      0.04000, 0.08000, 0.12000, 0.16000, 0.20000, &
     0.24000, 0.28000, 0.32000, 0.36000, 0.40000, 0.44000, &
     0.48000, 0.52000, 0.56000, 0.60000, 0.64000, 0.68000, &
     0.72000, 0.76000, 0.80000, 0.84000, 0.88000, 0.92000, &
     0.96000, 1.00000   /)
real, dimension(26) :: gw_flux_norm_zeta_s_04 = &
  (/ 0.0e000, 7.04e-6, 1.14e-5, 1.85e-5, 3.01e-5, 4.89e-5, &
     6.95e-5, 8.10e-5, 9.26e-5, 1.42e-4, 2.93e-4, 6.14e-4, &
     1.25e-3, 2.47e-3, 4.76e-3, 8.98e-3, 1.66e-2, 3.02e-2, &
     5.41e-2, 9.56e-2, 1.67e-1, 2.88e-1, 4.92e-1, 8.36e-1, &
     1.53e+0, 1.00e+1   /)
real, dimension(26) :: gw_area_norm_zeta_s_04 = &
  (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     3.48e-1, 1.00000  /)

integer :: num_sfc_layers, sub_layer_min, sub_layer_max
integer :: num_storage_pts, num_zeta_pts, num_tau_pts

real    :: zfull (max_lev)
real    :: zhalf (max_lev+1)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_soil_data_namelist(soil_num_l, soil_dz, soil_single_geo, &
                                   soil_gw_option )
  integer, intent(out) :: soil_num_l
  real,    intent(out) :: soil_dz(:)
  logical, intent(out) :: soil_single_geo
  integer, intent(out) :: soil_gw_option
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i, rcode, ncid, varid, dimids(3)
  real    :: z

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=soil_data_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=soil_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=soil_data_nml)
  
  ! register selector for all soil tiles
  call register_tile_selector('soil', long_name='soil',&
       tag = SEL_SOIL, idata1 = 0 )
  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_soil_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_SOIL, idata1 = i )
  enddo
  z = 0
  num_sfc_layers = 0
  sub_layer_min = 0
  sub_layer_max = 0

  zhalf(1) = 0
  do i = 1, num_l
    zhalf(i+1) = zhalf(i) + dz(i)
    zfull(i)   = 0.5*(zhalf(i+1) + zhalf(i))
    
    if (z < z_sub_layer_min+1.e-4) sub_layer_min = i
    z = z + dz(i)
    if (z < z_sfc_layer+1.e-4) num_sfc_layers = i
    if (z < z_sub_layer_max+1.e-4) sub_layer_max = i
  enddo

!!$  write (*,*) 'min/max index of layers whose thermal cond is scaled:',sub_layer_min,sub_layer_max

  if (trim(geohydrology_to_use)=='lm2') then
     gw_option = GW_LM2
  else if (trim(geohydrology_to_use)=='linear') then
     gw_option = GW_LINEAR
  else if (trim(geohydrology_to_use)=='hill_ar5') then
     gw_option = GW_HILL_AR5
  else if (trim(geohydrology_to_use)=='hill') then
     gw_option = GW_HILL
  else
     call error_mesg('read_soil_data_namelist',&
        'option geohydrology_to_use="'//trim(geohydrology_to_use)//'" is invalid, use '// &
        ' "lm2", "linear", "hill", or "hill_ar5"', &
        FATAL)
  endif
  
  if ((gw_option==GW_LM2).neqv.lm2) then
     call error_mesg('read_soil_data_namelist',&
        'geohydrlogy/LM2 option conflict: geohydrology_to_use must be consistent with the LM2 flag', &
        FATAL)
  endif

  if (gw_option==GW_HILL_AR5.and..not.use_single_geo) then
     num_storage_pts = 26
     call read_data('INPUT/geohydrology_table.nc', 'gw_flux_norm', &
                 gw_flux_table, no_domain=.true.)
     call read_data('INPUT/geohydrology_table.nc', 'gw_area_norm', &
                 gw_area_table, no_domain=.true.)
  else if (gw_option==GW_HILL.or.gw_option==GW_HILL_AR5) then
     __NF_ASRT__(nf_open('INPUT/geohydrology_table.nc',NF_NOWRITE,ncid))
     __NF_ASRT__(nf_inq_varid(ncid,'log_rho',varid))
     __NF_ASRT__(nf_inq_vardimid(ncid,varid,dimids))
     __NF_ASRT__(nf_inq_dimlen(ncid,dimids(1),num_storage_pts))
     __NF_ASRT__(nf_inq_dimlen(ncid,dimids(2),num_tau_pts))
     __NF_ASRT__(nf_inq_dimlen(ncid,dimids(3),num_zeta_pts))
     __NF_ASRT__(nf_close(ncid))
     
     allocate (log_rho_table(num_storage_pts,  &
                                 num_tau_pts, num_zeta_pts))
     allocate (log_deficit_list(num_storage_pts))
     allocate (log_tau(num_tau_pts))
     allocate (log_zeta_s(num_zeta_pts))
  
     call read_data('INPUT/geohydrology_table.nc', 'log_rho', &
                 log_rho_table, no_domain=.true.)
     call read_data('INPUT/geohydrology_table.nc', 'log_deficit', &
                 log_deficit_list, no_domain=.true.)
     call read_data('INPUT/geohydrology_table.nc', 'log_tau', &
                 log_tau, no_domain=.true.)
     call read_data('INPUT/geohydrology_table.nc', 'log_zeta_s', &
                 log_zeta_s, no_domain=.true.)
  endif


  ! set up output arguments
  soil_num_l      = num_l
  soil_dz         = dz
  soil_single_geo = use_single_geo
  soil_gw_option  = gw_option

end subroutine 


! ============================================================================
function soil_tile_ctor(tag) result(ptr)
  type(soil_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  allocate( ptr%prog(num_l))
  allocate( ptr%w_fc              (num_l),  &
            ptr%w_wilt            (num_l),  &
            ptr%d_trans           (num_l),  &
            ptr%alpha             (num_l),  &
            ptr%k_macro           (num_l),  &
            ptr%uptake_frac       (num_l),  &
            ptr%heat_capacity_dry (num_l),  &
            ptr%e                 (num_l),  &
            ptr%f                 (num_l),  &
            ptr%psi               (num_l),  &
            ptr%gw_flux_norm      (num_storage_pts),  &
            ptr%gw_area_norm      (num_storage_pts),  &
            ptr%fast_soil_C       (num_l),  &
            ptr%slow_soil_C       (num_l),  &
            ptr%fsc_in            (num_l),  & 
            ptr%ssc_in            (num_l),  & 
            ptr%asoil_in          (num_l)   )
  call soil_data_init_0d(ptr)
end function soil_tile_ctor


! ============================================================================
function soil_tile_copy_ctor(soil) result(ptr)
  type(soil_tile_type), pointer :: ptr ! return value
  type(soil_tile_type), intent(in) :: soil ! tile to copy

  allocate(ptr)
  ptr = soil ! copy all non-pointer members
  ! allocate storage for tile data
  allocate( ptr%prog(num_l))
  allocate( ptr%w_fc              (num_l),  &
            ptr%w_wilt            (num_l),  &
            ptr%d_trans           (num_l),  &
            ptr%alpha             (num_l),  &
            ptr%k_macro           (num_l),  &
            ptr%uptake_frac       (num_l),  &
            ptr%heat_capacity_dry (num_l),  &
            ptr%e                 (num_l),  &
            ptr%f                 (num_l),  &
            ptr%psi               (num_l),  &
            ptr%gw_flux_norm      (num_storage_pts),  &
            ptr%gw_area_norm      (num_storage_pts),  &
            ptr%fast_soil_C       (num_l),  &
            ptr%slow_soil_C       (num_l),  &
            ptr%fsc_in            (num_l),  & 
            ptr%ssc_in            (num_l),  & 
            ptr%asoil_in          (num_l)   )
  ! copy all pointer members
  ptr%prog(:) = soil%prog(:)
  ptr%w_fc(:) = soil%w_fc(:)
  ptr%w_wilt(:) = soil%w_wilt(:)
  ptr%d_trans(:) = soil%d_trans(:)
  ptr%alpha(:) = soil%alpha(:)
  ptr%k_macro(:) = soil%k_macro(:)
  ptr%uptake_frac(:) = soil%uptake_frac(:)
  ptr%uptake_T = soil%uptake_T
  ptr%heat_capacity_dry(:) = soil%heat_capacity_dry(:)
  ptr%e(:) = soil%e(:)
  ptr%f(:) = soil%f(:)
  ptr%psi(:) = soil%psi(:)
  ptr%gw_flux_norm(:) = soil%gw_flux_norm(:)
  ptr%gw_area_norm(:) = soil%gw_area_norm(:)
  ptr%fast_soil_C(:)  = soil%fast_soil_C(:)
  ptr%slow_soil_C(:)  = soil%slow_soil_C(:)
  ptr%fsc_in(:)       = soil%fsc_in(:)
  ptr%ssc_in(:)       = soil%ssc_in(:)
  ptr%asoil_in(:)     = soil%asoil_in(:)
end function soil_tile_copy_ctor


! ============================================================================
subroutine delete_soil_tile(ptr)
  type(soil_tile_type), pointer :: ptr

  deallocate(ptr%prog)
  deallocate(ptr%w_fc, ptr%w_wilt, ptr%d_trans, ptr%alpha, ptr%k_macro, &
             ptr%uptake_frac,&
             ptr%heat_capacity_dry, ptr%e, ptr%f, ptr%psi, &
             ptr%gw_flux_norm, ptr%gw_area_norm, &
             ptr%fast_soil_C, ptr%slow_soil_C, &
             ptr%fsc_in, ptr%ssc_in, ptr%asoil_in )
  deallocate(ptr)
end subroutine delete_soil_tile


! ============================================================================
subroutine soil_data_init_0d(soil)
  type(soil_tile_type), intent(inout) :: soil
  
!  real tau_groundwater
!  real rsa_exp         ! riparian source-area exponent
  integer :: k, i, l, code, m_zeta, m_tau
  real alpha_inf_sq, alpha_sfc_sq
  real :: single_log_zeta_s, single_log_tau, frac_zeta, frac_tau
  
  k = soil%tag

  soil%pars%vwc_sat           = dat_w_sat            (k)
  soil%pars%awc_lm2           = dat_awc_lm2          (k)
  soil%pars%k_sat_ref         = dat_k_sat_ref        (k)
  soil%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  soil%pars%chb               = dat_chb              (k)
  soil%pars%heat_capacity_dry = dat_heat_capacity_dry(k)
  soil%pars%thermal_cond_dry  = dat_thermal_cond_dry (k)
  soil%pars%thermal_cond_sat  = dat_thermal_cond_sat (k)
  soil%pars%thermal_cond_exp  = dat_thermal_cond_exp (k)
  soil%pars%thermal_cond_scale  = dat_thermal_cond_scale (k)
  soil%pars%thermal_cond_weight  = dat_thermal_cond_weight (k)
  soil%pars%refl_dry_dir      = dat_refl_dry_dir     (k,:)
  soil%pars%refl_dry_dif      = dat_refl_dry_dif     (k,:)
  soil%pars%refl_sat_dir      = dat_refl_sat_dir     (k,:)
  soil%pars%refl_sat_dif      = dat_refl_sat_dif     (k,:)
  soil%pars%emis_dry          = dat_emis_dry         (k)
  soil%pars%emis_sat          = dat_emis_sat         (k)
  soil%pars%z0_momentum       = dat_z0_momentum      (k)
  soil%pars%tfreeze           = tfreeze - dat_tf_depr(k)
  soil%pars%rsa_exp           = rsa_exp_global
  soil%pars%tau_groundwater   = gw_res_time
  soil%pars%hillslope_length  = gw_hillslope_length*gw_scale_length
  soil%pars%hillslope_zeta_bar= gw_hillslope_zeta_bar
  soil%pars%hillslope_relief  = gw_hillslope_relief*gw_scale_relief
  soil%pars%soil_e_depth      = gw_soil_e_depth*gw_scale_soil_depth
  soil%pars%k_sat_gw          = gw_perm*gw_scale_perm*9.8e9  ! m^2 to kg/(m2 s)
  soil%pars%storage_index     = 1
  soil%alpha                  = 1.0
  soil%fast_soil_C(:)         = 0.0
  soil%slow_soil_C(:)         = 0.0
  soil%asoil_in(:)            = 0.0
  soil%fsc_in(:)              = 0.0
  soil%ssc_in(:)              = 0.0
  do l = 1, num_l
    soil%k_macro(l) = k_macro_constant &
                         * exp(zfull(l)/soil%pars%soil_e_depth)
  enddo

  select case(gw_option)
  case(GW_HILL_AR5)
     if (use_single_geo) then
        soil%gw_flux_norm = gw_flux_norm_zeta_s_04
        soil%gw_area_norm = gw_area_norm_zeta_s_04
     endif
     soil%pars%zeta = soil%pars%soil_e_depth / soil%pars%hillslope_relief
  case(GW_HILL)
     if (use_single_geo) &
         call soil_data_init_derive_subsurf_pars ( soil )
  case default
     ! do nothing
  end select

  ! ---- derived constant soil parameters
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     soil%w_wilt(:) = 0.15
     soil%w_fc  (:) = 0.15 + soil%pars%awc_lm2
  else
     soil%w_wilt(:) = soil%pars%vwc_sat &
          *(soil%pars%psi_sat_ref/(psi_wilt*soil%alpha(:)))**(1/soil%pars%chb)
     soil%w_fc  (:) = soil%pars%vwc_sat &
          *(rate_fc/(soil%pars%k_sat_ref*soil%alpha(:)**2))**(1/(3+2*soil%pars%chb))
  endif

  soil%pars%vwc_wilt = soil%w_wilt(1)
  soil%pars%vwc_fc   = soil%w_fc  (1)

  soil%pars%vlc_min = soil%pars%vwc_sat*K_rel_min**(1/(3+2*soil%pars%chb))

  soil%z0_scalar = soil%pars%z0_momentum * exp(-k_over_B)
  
  soil%geothermal_heat_flux = geothermal_heat_flux_constant

end subroutine soil_data_init_0d

! ============================================================================
subroutine soil_data_init_derive_subsurf_pars ( soil )
  type(soil_tile_type), intent(inout) :: soil
  integer i, l, code, m_zeta, m_tau
  real :: alpha_inf_sq, alpha_sfc_sq
  real :: single_log_zeta_s, single_log_tau, frac_zeta, frac_tau

  if (use_alpha) then
      soil%pars%k_sat_sfc = (soil%pars%k_sat_ref-soil%pars%k_sat_gw) &
                * exp(z_ref/soil%pars%soil_e_depth) + soil%pars%k_sat_gw
      alpha_inf_sq = soil%pars%k_sat_gw / soil%pars%k_sat_ref
      alpha_sfc_sq = soil%pars%k_sat_sfc / soil%pars%k_sat_ref
      do l = 1, num_l
        soil%alpha(l) = sqrt(alpha_inf_sq+(alpha_sfc_sq-alpha_inf_sq) &
                    *exp(-zfull(l)/soil%pars%soil_e_depth))
        enddo
    else
      soil%pars%k_sat_sfc = soil%pars%k_sat_ref
      soil%alpha = 1.0
    endif

  soil%pars%tau =    &
    (soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length) &
     / ((soil%pars%k_sat_sfc+k_macro_constant)*soil%pars%soil_e_depth)

  soil%d_trans(1) = 1.
  do l = 1, num_l-1
    soil%d_trans(l+1) = exp(-zhalf(l+1)/soil%pars%soil_e_depth)
    soil%d_trans(l)   = soil%d_trans(l) - soil%d_trans(l+1)
    enddo
  soil%d_trans = soil%d_trans &
       * (soil%pars%k_sat_sfc + k_macro_constant ) &
       * soil%pars%soil_e_depth
  soil%d_trans(num_l) = soil%d_trans(num_l) &
       + soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length

  if (soil%pars%hillslope_relief.le.0.) soil%pars%hillslope_relief = 1.e-10
  soil%pars%zeta = soil%pars%soil_e_depth &
                                        / soil%pars%hillslope_relief
  single_log_zeta_s = log10(soil%pars%zeta)
  single_log_zeta_s = max(single_log_zeta_s, log_zeta_s(1))
  single_log_zeta_s = min(single_log_zeta_s, log_zeta_s(num_zeta_pts))
  m_zeta = num_zeta_pts / 2
  code = 0
  do while (code.eq.0)
    if (single_log_zeta_s .lt. log_zeta_s(m_zeta)) then
        m_zeta = m_zeta - 1
      else if (single_log_zeta_s .gt. log_zeta_s(m_zeta+1)) then
        m_zeta = m_zeta + 1
      else
        code = 1
      endif
    enddo
  frac_zeta = (single_log_zeta_s - log_zeta_s(m_zeta)) &
     / (log_zeta_s(m_zeta+1) - log_zeta_s(m_zeta))
  single_log_tau = log10( soil%pars%tau )
  single_log_tau = max(single_log_tau, log_tau(1))
  single_log_tau = min(single_log_tau, log_tau(num_tau_pts))
  m_tau = num_tau_pts / 2
  code = 0
  do while (code.eq.0)
    if (single_log_tau .lt. log_tau(m_tau)) then
        m_tau = m_tau - 1
      else if (single_log_tau .gt. log_tau(m_tau+1)) then
        m_tau = m_tau + 1
      else
        code = 1
      endif
    enddo
  frac_tau = (single_log_tau - log_tau(m_tau)) &
     / (log_tau(m_tau+1) - log_tau(m_tau))
  do i = 1, num_storage_pts
    soil%gw_flux_norm(i) = min ( log_rho_max,                      &
       (1.-frac_zeta)*(1.-frac_tau)*log_rho_table(i , m_tau  , m_zeta  ) &
     + (   frac_zeta)*(1.-frac_tau)*log_rho_table(i , m_tau  , m_zeta+1) &
     + (1.-frac_zeta)*(   frac_tau)*log_rho_table(i , m_tau+1, m_zeta  ) &
     + (   frac_zeta)*(   frac_tau)*log_rho_table(i , m_tau+1, m_zeta+1) )
    enddo

end subroutine soil_data_init_derive_subsurf_pars

! ============================================================================
subroutine soil_data_init_derive_subsurf_pars_ar5 ( soil )
  type(soil_tile_type), intent(inout) :: soil
  integer i, code, m_zeta
  real :: frac_zeta
  if (soil%pars%hillslope_relief.le.0.) &
     soil%pars%hillslope_relief =      &
        soil%pars%soil_e_depth / 10.**log_zeta_s(num_zeta_pts)
  soil%pars%zeta = soil%pars%soil_e_depth &
     / soil%pars%hillslope_relief
  soil%pars%zeta = max(soil%pars%zeta, gw_zeta_s(1))
  soil%pars%zeta = min(soil%pars%zeta, gw_zeta_s(31))
  m_zeta = 31 / 2
  code = 0
  do while (code.eq.0)
    if (soil%pars%zeta .lt. gw_zeta_s(m_zeta)) then
        m_zeta = m_zeta - 1
      else if (soil%pars%zeta .gt. gw_zeta_s(m_zeta+1)) then
        m_zeta = m_zeta + 1
      else
        code = 1
      endif
    enddo
  frac_zeta = (soil%pars%zeta - gw_zeta_s(m_zeta)) &
                / (gw_zeta_s(m_zeta+1) - gw_zeta_s(m_zeta))
  do i = 1, num_storage_pts
    soil%gw_flux_norm(i) = gw_flux_table(i,m_zeta) &
         + frac_zeta*(gw_flux_table(i,m_zeta+1)-gw_flux_table(i,m_zeta))
    soil%gw_area_norm(i) = gw_area_table(i,m_zeta) &
         + frac_zeta*(gw_area_table(i,m_zeta+1)-gw_area_table(i,m_zeta))
    enddo
end subroutine soil_data_init_derive_subsurf_pars_ar5

! ============================================================================
function soil_cover_cold_start(land_mask, lonb, latb) result (soil_frac)
! creates and initializes a field of fractional soil coverage
  logical, intent(in) :: land_mask(:,:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  real,    pointer    :: soil_frac (:,:,:) ! output: map of soil fractional coverage

  allocate( soil_frac(size(land_mask,1),size(land_mask,2),n_dim_soil_types))

  call init_cover_field(soil_to_use, 'INPUT/ground_type.nc', 'cover','frac', &
       lonb, latb, soil_index_constant, input_cover_types, soil_frac)
  
end function 


! =============================================================================
function soil_tiles_can_be_merged(soil1,soil2) result(response)
  logical :: response
  type(soil_tile_type), intent(in) :: soil1,soil2

  response = (soil1%tag==soil2%tag)
end function


! =============================================================================
subroutine merge_soil_tiles(s1,w1,s2,w2)
  type(soil_tile_type), intent(in) :: s1
  type(soil_tile_type), intent(inout) :: s2
  real                , intent(in) :: w1,w2

  ! ---- local vars
  real    :: x1, x2 ! normalized relative weights
  real    :: gw, HEAT1, HEAT2 ! temporaries for groundwater and heat
  integer :: i
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! combine state variables
  do i = 1,num_l
     ! calculate heat content at this level for both source tiles
     HEAT1 = &
          (s1%heat_capacity_dry(i)*dz(i)+clw*s1%prog(i)%Wl+csw*s1%prog(i)%Ws)* &
          (s1%prog(i)%T-tfreeze)
     HEAT2 = &
          (s2%heat_capacity_dry(i)*dz(i)+clw*s2%prog(i)%Wl+csw*s2%prog(i)%Ws)* &
          (s2%prog(i)%T-tfreeze)
     ! merge the amounts of water
     s2%prog(i)%Wl = x1*s1%prog(i)%Wl + x2*s2%prog(i)%Wl
     s2%prog(i)%Ws = x1*s1%prog(i)%Ws + x2*s2%prog(i)%Ws
     ! if the dry heat capacity of merged soil is to be changed, do it here
     ! ...
     ! calculate the merged temperature based on heat content
     s2%prog(i)%T = tfreeze + (x1*HEAT1+x2*HEAT2)/ &
          (s2%heat_capacity_dry(i)*dz(i)+clw*s2%prog(i)%Wl+csw*s2%prog(i)%Ws)

     ! calculate combined groundwater content
     gw = s1%prog(i)%groundwater*x1 + s2%prog(i)%groundwater*x2
     ! calculate combined groundwater temperature
     if (gw/=0) then
        s2%prog(i)%groundwater_T = ( &
             s1%prog(i)%groundwater*x1*(s1%prog(i)%groundwater_T-tfreeze) + &
             s2%prog(i)%groundwater*x2*(s2%prog(i)%groundwater_T-tfreeze)   &
             ) / gw + tfreeze
     else
        s2%prog(i)%groundwater_T = &
             s1%prog(i)%groundwater_T*x1 + s2%prog(i)%groundwater_T*x2
     endif
     s2%prog(i)%groundwater = gw

  enddo
  s2%uptake_T    = s1%uptake_T*x1 + s2%uptake_T*x2
  ! merge soil carbon
  s2%fast_soil_C(:) = s1%fast_soil_C(:)*x1 + s2%fast_soil_C(:)*x2
  s2%slow_soil_C(:) = s1%slow_soil_C(:)*x1 + s2%slow_soil_C(:)*x2
  s2%asoil_in(:)    = s1%asoil_in(:)*x1 + s2%asoil_in(:)*x2
  s2%fsc_in(:)      = s1%fsc_in(:)*x1 + s2%fsc_in(:)*x2
  s2%ssc_in(:)      = s1%ssc_in(:)*x1 + s2%ssc_in(:)*x2
end subroutine

! =============================================================================
! returns true if tile fits the specified selector
function soil_is_selected(soil, sel)
  logical soil_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(soil_tile_type),      intent(in) :: soil

  soil_is_selected = (sel%idata1==0).or.(sel%idata1==soil%tag)
end function


! ============================================================================
! returns tag of the tile
function get_soil_tile_tag(soil) result(tag)
  integer :: tag
  type(soil_tile_type), intent(in) :: soil
  
  tag = soil%tag
end function



! ============================================================================
! compute average soil temperature with a given depth scale
function soil_ave_temp(soil, depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0 
  do k = 1, num_l
     w = dz(k) * exp(-zfull(k)/depth)
     A = A + soil%prog(k)%T * w
     N = N + w
     if (zhalf(k+1).gt.depth) exit
  enddo
  A = A/N
end function soil_ave_temp


! ============================================================================
! compute average soil moisture with a given depth scale
function soil_ave_theta0(soil, zeta) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: zeta ! root depth scale
  real :: zeta2  !  adjusted length scale for soil-water-stress
  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0 
  zeta2 = zeta*zeta_mult
  do k = 1, num_l
     w = exp(-zhalf(k)/zeta2)-exp(-zhalf(k+1)/zeta2)
     A = A + max(soil%prog(k)%wl/(dens_h2o*dz(k))-soil%w_wilt(k),0.0)/&
          (soil%w_fc(k)-soil%w_wilt(k)) * w
     N = N + w
  enddo
  A = A/N
end function soil_ave_theta0


! ============================================================================
function soil_ave_theta1(soil, depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0 
  do k = 1, num_l
     w = dz(k) * exp(-zfull(k)/depth)
     A = A +min(max(soil%prog(k)%wl/(dens_h2o*dz(k)),0.0)/&
          (soil%pars%vwc_sat),1.0) * w
     N = N + w
     if (zhalf(k+1).gt.depth) exit
  enddo
  A = A/N
end function soil_ave_theta1


! ============================================================================
! returns array of soil moisture
function soil_theta(soil) result (theta1)
  type(soil_tile_type), intent(in) :: soil
  real :: theta1(num_l)

  theta1(:) = min(max(soil%prog(:)%wl/(dens_h2o*dz(:)),0.0)/(soil%pars%vwc_sat),1.0)
end function soil_theta


! ============================================================================
function soil_psi_stress(soil, zeta) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real,                 intent(in) :: zeta  ! root-mass depth scale
  real :: zeta2  !  adjusted length scale for soil-water-stress
  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0
  zeta2 = zeta*zeta_mult
  do k = 1, num_l
     w = exp(-zhalf(k)/zeta2)-exp(-zhalf(k+1)/zeta2)
     A = A + w * soil%psi(k)/psi_wilt
     N = N + w
     enddo
  A = A/N

end function soil_psi_stress

! ============================================================================
! compute bare-soil albedo, bare-soil emissivity, bare-soil roughness
! for scalar transport, and beta function
subroutine soil_data_radiation ( soil, cosz, use_brdf, soil_alb_dir, soil_alb_dif, soil_emis )
  type(soil_tile_type), intent(in)  :: soil
  real,                 intent(in)  :: cosz
  logical,              intent(in)  :: use_brdf
  real,                 intent(out) :: soil_alb_dir(NBANDS), soil_alb_dif(NBANDS), soil_emis
  ! ---- local vars
  real :: soil_sfc_vlc, blend, dry_value(NBANDS), sat_value(NBANDS)
  real :: zenith_angle, zsq, zcu

  soil_sfc_vlc  = soil%prog(1)%wl/(dens_h2o*dz(1))
  blend         = max(0., min(1., soil_sfc_vlc/soil%pars%vwc_sat))
  if (use_brdf) then
      zenith_angle = acos(cosz)
      zsq = zenith_angle*zenith_angle
      zcu = zenith_angle*zsq
      dry_value =  soil%pars%f_iso_dry*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                 + soil%pars%f_vol_dry*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                 + soil%pars%f_geo_dry*(g0_geo+g1_geo*zsq+g2_geo*zcu)
      sat_value =  soil%pars%f_iso_sat*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                 + soil%pars%f_vol_sat*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                 + soil%pars%f_geo_sat*(g0_geo+g1_geo*zsq+g2_geo*zcu)
    else
      dry_value = soil%pars%refl_dry_dir
      sat_value = soil%pars%refl_sat_dir
    endif
  soil_alb_dir  = dry_value              + blend*(sat_value             -dry_value)
  soil_alb_dif  = soil%pars%refl_dry_dif + blend*(soil%pars%refl_sat_dif-soil%pars%refl_dry_dif)
  soil_emis     = soil%pars%emis_dry     + blend*(soil%pars%emis_sat    -soil%pars%emis_dry    )
end subroutine soil_data_radiation


! ============================================================================
! compute bare-soil albedo, bare-soil emissivity, bare-soil roughness
! for scalar transport, and beta function
subroutine soil_data_diffusion ( soil, soil_z0s, soil_z0m )
  type(soil_tile_type), intent(in)  :: soil
  real,                 intent(out) :: soil_z0s, soil_z0m

  soil_z0s = soil%z0_scalar
  soil_z0m = soil%pars%z0_momentum
end subroutine soil_data_diffusion

! ============================================================================
! compute soil thermodynamic properties.
subroutine soil_data_thermodynamics ( soil, vlc, vsc, &
                                      soil_E_max, thermal_cond)
  type(soil_tile_type), intent(inout) :: soil
  real,                 intent(in)  :: vlc(:)
  real,                 intent(in)  :: vsc(:)
  real,                 intent(out) :: soil_E_max
  real,                 intent(out) :: thermal_cond(:)
  real s, w, a, n, f

  integer l

  ! assign some index of water availability for snow-free soil

  soil_E_max = (soil%pars%k_sat_ref*soil%alpha(1)**2) &
               * (-soil%pars%psi_sat_ref/soil%alpha(1)) &
               * ((4.+soil%pars%chb)*vlc(1)/ &
                ((3.+soil%pars%chb)*soil%pars%vwc_sat))**(3.+soil%pars%chb) &
                / ((1.+3./soil%pars%chb)*dz(1))

     w = soil%pars%thermal_cond_weight
     a = soil%pars%thermal_cond_scale
     n = soil%pars%thermal_cond_exp
  do l = 1, num_sfc_layers
     soil%heat_capacity_dry(l) = sfc_heat_factor*soil%pars%heat_capacity_dry
     s = (vlc(l)+vsc(l))/soil%pars%vwc_sat
     thermal_cond(l)      = sfc_heat_factor * &
          ( soil%pars%thermal_cond_dry+ &
            (soil%pars%thermal_cond_sat-soil%pars%thermal_cond_dry) &
            *(w*s +(1-w)*(1+a**n)*(s**n)/(1+(a*s)**n))    )
     f = 1.
     if (vlc(l)+vsc(l).gt.0.) f = 1.+(freeze_factor-1.)*vsc(l)/(vlc(l)+vsc(l))
     thermal_cond(l) = f * thermal_cond(l)
  enddo
  do l = num_sfc_layers+1, num_l
     soil%heat_capacity_dry(l) = soil%pars%heat_capacity_dry
     s = (vlc(l)+vsc(l))/soil%pars%vwc_sat
     thermal_cond(l)  = &
          ( soil%pars%thermal_cond_dry+ &
            (soil%pars%thermal_cond_sat-soil%pars%thermal_cond_dry) &
            *(w*s +(1-w)*(1+a**n)*(s**n)/(1+(a*s)**n))    )
     f = 1.
     if (vlc(l)+vsc(l).gt.0.) f = 1.+(freeze_factor-1.)*vsc(l)/(vlc(l)+vsc(l))
     thermal_cond(l) = f * thermal_cond(l)
  enddo
  
  ! this is an additional factor intended for tuning annual T range in
  ! high latitudes. presumably other locations are insensitive to this
  ! global parameter, since they don't have freeze/thaw. this really is just a fudge.
  do l = sub_layer_min, sub_layer_max
    thermal_cond(l) = sub_layer_tc_fac * thermal_cond(l)
    enddo
  
end subroutine soil_data_thermodynamics


! ============================================================================
! compute soil hydraulic properties: wrapper to get all parameters for
! richards equation for full column. note that psi_for_rh is not passed.
subroutine soil_data_hydraulic_properties (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, hyd_cond, DKDP
  real,                        intent(out) :: &
      DPsi_min, DPsi_max
  ! ---- local vars ----------------------------------------------------------
  real ::  psi_for_rh

  call soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )

end subroutine soil_data_hydraulic_properties


! ============================================================================
! wrapper to compute psi and the potentially different psi used for surface
! relative humidity computation.
! use of this wrapper allows us to avoid passing all the other properties
! back to soil module. (they're computed anyway.)
subroutine soil_data_psi_for_rh (soil, vlc, vsc, psi, psi_for_rh  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: psi
  real,                        intent(out) :: psi_for_rh
  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) ::  DThDP, hyd_cond, DKDP
  real :: DPsi_min, DPsi_max

  call soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )

end subroutine soil_data_psi_for_rh


! ============================================================================
! compute soil hydraulic properties.
subroutine soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, hyd_cond, DKDP
  real,                        intent(out) :: &
      DPsi_min, DPsi_max, psi_for_rh
  ! ---- local vars ----------------------------------------------------------
  integer l
  real :: vlc_loc, vlc_k, psi_k, sigma, B, por, psi_s, k_sat, alt_psi_for_rh
  real :: alpha_sq, f_psi
  logical flag
  
  ! ---- T-dependence of hydraulic properties --------------------------------
  ! k_sat   = soil%pars%k_sat0   !  * mu(t0)/mu(t), where mu is dynamic viscosity
  ! psi_sat = soil%pars%psi_sat0 !  * exp(c*(psi-psi0)), where c~+/-(?)0.0068
                     ! better approach would be to adopt air entrapment model
                     ! or at least to scale against surface tension model
  ! ---- water and ice dependence of hydraulic properties --------------------
  ! ---- (T-dependence can be added later)
  hyd_cond=1;DThDP=1;psi=1
  DKDP=0
  flag = .false.
  do l = 1, size(vlc)
    alpha_sq = soil%alpha(l)**2
    hyd_cond(l) = (soil%pars%k_sat_ref*alpha_sq)*  &
               (vlc(l)/soil%pars%vwc_sat)**(3+2*soil%pars%chb)
    if (hyd_cond(l).lt.1.e-12*soil%pars%k_sat_ref*alpha_sq) then
        vlc_loc     = soil%pars%vwc_sat*(1.e-12)**(1./(3+2*soil%pars%chb))
        hyd_cond(l) = 1.e-12*soil%pars%k_sat_ref*alpha_sq
        if (l.eq.1) flag = .true.
        if (vsc(l).eq.0.) then
            DThDP   (l) = -vlc_loc     &
                      *(vlc_loc/soil%pars%vwc_sat)**soil%pars%chb &
                      /(soil%pars%psi_sat_ref*soil%pars%chb/soil%alpha(l))
            psi     (l) = (soil%pars%psi_sat_ref/soil%alpha(l)) &
                      *(soil%pars%vwc_sat/vlc_loc)**soil%pars%chb &
                      + (vlc(l)-vlc_loc)/DThDP   (l)
            if (l.eq.1.and.vlc(1).gt.0.) then
                alt_psi_for_rh = &
                     (soil%pars%psi_sat_ref/soil%alpha(l)) &
                     *(soil%pars%vwc_sat/vlc(1)   )**soil%pars%chb
              else if (l.eq.1.and.vlc(1).le.0.) then
                alt_psi_for_rh = -1.e10
              endif
          else
            psi     (l) = ((soil%pars%psi_sat_ref/soil%alpha(l)) / 2.2) &
                    *(soil%pars%vwc_sat/vlc_loc   )**soil%pars%chb
            DThDP   (l) = 0.
            if (l.eq.1) alt_psi_for_rh = -1.e10
          endif
      else
          if (vsc(l).eq.0.) then
              if (vlc(l).le.soil%pars%vwc_sat) then
                  psi     (l) = (soil%pars%psi_sat_ref/soil%alpha(l)) &
                          *(soil%pars%vwc_sat/vlc(l))**soil%pars%chb
                  DKDP    (l) = -(2+3/soil%pars%chb)*hyd_cond(l) &
                                                 /psi(l)
                  DThDP   (l) = -vlc(l)/(psi(l)*soil%pars%chb)
                else
                  psi(l) = soil%pars%psi_sat_ref/soil%alpha(l) &
                         + (vlc(l)-soil%pars%vwc_sat)/comp
	          f_psi = min(max(1.-psi(l)/soil%pars%psi_sat_ref/soil%alpha(l),0.),1.)
                  DThDP(l) = comp
                  hyd_cond(l) = soil%pars%k_sat_ref*alpha_sq &
	                           + f_psi * soil%k_macro(l)
                endif
            else
              psi     (l) = ((soil%pars%psi_sat_ref/soil%alpha(l)) / 2.2) &
                        *(soil%pars%vwc_sat/vlc(l))**soil%pars%chb
              DThDP   (l) = 0.
            endif
      endif
    enddo

  if (use_alt_psi_for_rh .and. flag) then
      psi_for_rh = alt_psi_for_rh
    else
      psi_for_rh = psi(1)
    endif

  if (DThDP(1).ne.0.) then
      DPsi_min =            -vlc(1) /DThDP(1)
      DPsi_max = (soil%pars%vwc_sat-vlc(1))/DThDP(1)
    else
      Dpsi_min = Dpsi_min_const
      DPsi_max = -psi(1)
    endif

end subroutine soil_data_hydraulics

! ============================================================================
! compute soil hydraulic properties (undeveloped code)
subroutine soil_data_hydraulics_EXPERIMENTAL (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, hyd_cond, DKDP
  real,                        intent(out) :: &
      DPsi_min, DPsi_max, psi_for_rh
  ! ---- local vars ----------------------------------------------------------
  integer l
  real :: vlc_loc, vlc_k, psi_k, sigma, B, por, psi_s, k_sat, alt_psi_for_rh
  real :: alpha_sq, f_psi
  logical flag
  
  ! ---- T-dependence of hydraulic properties --------------------------------
  ! k_sat   = soil%pars%k_sat0   !  * mu(t0)/mu(t), where mu is dynamic viscosity
  ! psi_sat = soil%pars%psi_sat0 !  * exp(c*(psi-psi0)), where c~+/-(?)0.0068
                     ! better approach would be to adopt air entrapment model
                     ! or at least to scale against surface tension model


  ! ---- water and ice dependence of hydraulic properties --------------------
  ! ---- (T-dependence can be added later)
  hyd_cond=1;DThDP=1;psi=1
    B     = soil%pars%chb
    por   = soil%pars%vwc_sat
    vlc_k = soil%pars%vlc_min
    psi_s = soil%pars%psi_sat_ref
    psi_k = psi_s*(por/vlc_k)**B
    k_sat = soil%pars%k_sat_ref
    do l = 1, num_l
      vlc_loc = max(vlc(l), vlc_k)
      ! sigma is an adjustment to surface tension for presence of ice
      sigma = 1. + (sigma_max-1.)*min(vsc(l)/vlc(l),1.)
      if (vlc(l).lt.vlc_k) then  ! very dry, no ice, sigma=1 in this case
          DThDP(l) = -vlc_k/(B*psi_k)
          psi  (l) = psi_s*(por/vlc_k)**B + (vlc(l)-vlc_k)/DThDP(l)
        else if (vlc(l).lt.por-vsc(l)) then  ! unsaturated, maybe with ice
          psi  (l) = (psi_s/sigma)*(por/vlc(l))**B
          DThDP(l) = -vlc(l)/(B*psi(l))
        else  ! no air is present in this case
          psi  (l) = (psi_s/sigma)*(por/(por-vsc(l)))**B &
                             + (vlc(l)+vsc(l)-por)/comp
          DThDP(l) = comp
        endif
      hyd_cond(l) = k_sat*(vlc_loc/por)**(3+2*B)
      DKDP(l) = 0
      enddo
    DPsi_min =            -vlc(1) /DThDP(1)
    if (vsc(1).gt.0.) DPsi_min = (vlc_k-vlc(1)) /DThDP(1)
    DPsi_max = (por-vsc(1)-vlc(1))/DThDP(1)
    psi_for_rh = psi(1)

end subroutine soil_data_hydraulics_EXPERIMENTAL

! ============================================================================
subroutine soil_data_gw_hydraulics_ar5(soil, storage_normalized, &
                                                     gw_flux, sat_frac )
  type(soil_tile_type), intent(inout)  :: soil
  real,                 intent(in)  :: storage_normalized
  real,                 intent(out) :: gw_flux
  real,                 intent(out) :: sat_frac

  integer :: code, m
  real :: recharge_normalized, frac

  ! storage_normalized is the fraction of soil above drainage base elevation
  ! that is below the water table
  code = 0
  m = soil%pars%storage_index
  do while (code.eq.0)
    if (storage_normalized .lt. gw_storage_norm(m)) then
        m = m - 1
      else if (storage_normalized .gt. gw_storage_norm(m+1)) then
        m = m + 1
      else
        code = 1
      endif
    enddo
  if (m.lt.1.or.m.gt.num_storage_pts-1) then
      write(*,*) '!!! *** m=',m, ' is outside the table in soil_data_gw_hydraulics_ar5 *** !!!'
      write(*,*) 'num_storage_pts=',num_storage_pts
      write(*,*) 'storage_normalized=',storage_normalized
      write(*,*) 'interval bounds:',gw_storage_norm(m),gw_storage_norm(m+1)
    endif
  frac = (storage_normalized-gw_storage_norm(m)) &
           /(gw_storage_norm(m+1)-gw_storage_norm(m))
  sat_frac = soil%gw_area_norm(m) &
               + frac*(soil%gw_area_norm(m+1)-soil%gw_area_norm(m))
  recharge_normalized = soil%gw_flux_norm(m) &
               + frac*(soil%gw_flux_norm(m+1)-soil%gw_flux_norm(m))
  gw_flux = recharge_normalized * soil%pars%k_sat_ref * soil%pars%soil_e_depth &
                * soil%pars%hillslope_relief &
                   / (soil%pars%hillslope_length * soil%pars%hillslope_length)
  soil%pars%storage_index = m
  
end subroutine soil_data_gw_hydraulics_ar5

! ============================================================================
subroutine soil_data_gw_hydraulics(soil, deficit_normalized, &
                                                     gw_flux, sat_frac )
  type(soil_tile_type), intent(inout)  :: soil
  real,                 intent(in)  :: deficit_normalized
  real,                 intent(out) :: gw_flux
  real,                 intent(out) :: sat_frac

  integer :: code, m
  real :: recharge_normalized, frac, log_deficit_normalized

  ! deficit_normalized is the fraction of soil above drainage base elevation
  ! that is NOT saturated/contributing to horizontal flow
  
  IF (deficit_normalized .LE. 0.) THEN
    m = 0
  ELSE IF (deficit_normalized .GE. 1.) THEN
    m = num_storage_pts
  ELSE 
    log_deficit_normalized = log10(deficit_normalized)
    code = 0
    m = soil%pars%storage_index
    do while (code.eq.0.and.m.gt.0.and.m.lt.num_storage_pts)
      if (log_deficit_normalized .lt. log_deficit_list(m)) then
          m = m - 1
        else if (log_deficit_normalized .gt. log_deficit_list(m+1)) then
          m = m + 1
        else
          code = 1
        endif
      enddo
  ENDIF

  IF (m.eq.0) THEN
      recharge_normalized = 10.**soil%gw_flux_norm(1)
    ELSE IF (m.lt.num_storage_pts) then
      frac = (log_deficit_normalized-log_deficit_list(m)) &
               /(log_deficit_list(m+1)-log_deficit_list(m))
      recharge_normalized = 10.**(soil%gw_flux_norm(m) &
                   + frac*(soil%gw_flux_norm(m+1)-soil%gw_flux_norm(m)))
    ELSE
      recharge_normalized = 0.
    ENDIF

  if (use_tau_fix) then
      gw_flux = recharge_normalized * &
         ((soil%pars%k_sat_sfc+k_macro_constant) * soil%pars%soil_e_depth + &
	           soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length) &
         * soil%pars%hillslope_relief &
         / (soil%pars%hillslope_length * soil%pars%hillslope_length)
    else
      gw_flux = recharge_normalized * &
         ((soil%pars%k_sat_sfc+k_macro_constant) * soil%pars%soil_e_depth) &
         * soil%pars%hillslope_relief &
         / (soil%pars%hillslope_length * soil%pars%hillslope_length)
    endif

  if (recharge_normalized.gt.1.) then
      sat_frac = 1. - 1./recharge_normalized
    else
      sat_frac = 0.
    endif
  if (use_sat_fix) then
      if (recharge_normalized.gt.1.+soil%pars%tau) then
          sat_frac = 1. - (1.+soil%pars%tau)/recharge_normalized
        else
          sat_frac = 0.
        endif
    endif

  soil%pars%storage_index = min(max(1,m),num_storage_pts-1)
  
!          if(is_watch_point()) then
!	         write(*,*) 'm            ', m
!		 write(*,*) 'code         ',code
!		 write(*,*) 'num_storage_pts', num_storage_pts
!		 write(*,*) 'frac         ', frac
!		 write(*,*) 'recharge_normalized', recharge_normalized
!		 write(*,*) 'deficit_normalized ', deficit_normalized
!		 write(*,*) 'log_deficit_normalized',log_deficit_normalized
!		 write(*,*) 'gw_flux      ',gw_flux
!		 write(*,*) 'm,log_deficit_list,gw_flux_norm'
!		 do m=1,num_storage_pts
!		   write(*,*) m, log_deficit_list(m),soil%gw_flux_norm(m)
!		   enddo
!               endif

  
end subroutine soil_data_gw_hydraulics

! ============================================================================
subroutine soil_data_vwc_for_init_only (soil, psi, vwc)
  type(soil_tile_type), intent(in) :: soil
  real,                 intent(in) :: psi(:)
  real,                 intent(out):: vwc(:)
  vwc = soil%pars%vwc_sat* &
    ((soil%pars%psi_sat_ref/soil%alpha)/min(psi,(soil%pars%psi_sat_ref/soil%alpha))) &
     ** (1./ soil%pars%chb)
  if (use_comp_for_ic) then
      vwc = vwc + comp*max(psi-soil%pars%psi_sat_ref/soil%alpha,0.)
    endif
end subroutine soil_data_vwc_for_init_only

! ============================================================================
subroutine soil_tile_stock_pe (soil, twd_liq, twd_sol  )
  type(soil_tile_type),  intent(in)    :: soil
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n
  
  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(soil%prog)
    twd_liq = twd_liq + soil%prog(n)%wl + soil%prog(n)%groundwater
    twd_sol = twd_sol + soil%prog(n)%ws
    enddo

end subroutine soil_tile_stock_pe


! ============================================================================
! returns soil tile heat content, J/m2
function soil_tile_heat (soil) result(heat) ; real heat
  type(soil_tile_type),  intent(in)  :: soil

  integer :: i

  heat = 0
  do i = 1, num_l
     heat = heat + &
          (soil%heat_capacity_dry(i)*dz(i)+clw*soil%prog(i)%Wl+csw*soil%prog(i)%Ws)&
                           *(soil%prog(i)%T-tfreeze) + &
          clw*soil%prog(i)%groundwater*(soil%prog(i)%groundwater_T-tfreeze) - &
          hlf*soil%prog(i)%ws
  enddo
end function

! ============================================================================
! returns soil tile carbon content, kg C/m2
function soil_tile_carbon (soil); real soil_tile_carbon
  type(soil_tile_type),  intent(in)  :: soil

  soil_tile_carbon = sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:))
end function

end module soil_tile_mod
