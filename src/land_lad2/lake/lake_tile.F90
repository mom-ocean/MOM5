#include <fms_platform.h>

module lake_tile_mod

use mpp_domains_mod, only : &
     domain2d, mpp_get_compute_domain, mpp_global_field

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : &
     write_version_number, file_exist, check_nml_error, &
     read_data, close_file, stdlog
use constants_mod, only : &
     pi, tfreeze, hlf
use land_constants_mod, only : &
     NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_LAKE, register_tile_selector

implicit none
private

! ==== public interfaces =====================================================
public :: lake_pars_type
public :: lake_prog_type
public :: lake_tile_type

public :: new_lake_tile, delete_lake_tile
public :: lake_tiles_can_be_merged, merge_lake_tiles
public :: lake_is_selected
public :: get_lake_tile_tag
public :: lake_tile_stock_pe
public :: lake_tile_heat

public :: read_lake_data_namelist
public :: lake_cover_cold_start

public :: lake_data_radiation
public :: lake_data_diffusion
public :: lake_data_thermodynamics

public :: max_lev
public :: lake_width_inside_lake
public :: large_lake_sill_width
! =====end of public interfaces ==============================================
interface new_lake_tile
   module procedure lake_tile_ctor
   module procedure lake_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), private, parameter   :: &
     version     = '$Id: lake_tile.F90,v 19.0 2012/01/06 20:40:53 fms Exp $', &
     tagname     = '$Name: siena_201207 $', &
     module_name = 'lake_tile_mod'

integer, parameter :: max_lev          = 80
integer, parameter :: n_dim_lake_types = 1  ! size of lookup table
real,    parameter :: psi_wilt         = -150.  ! matric head at wilting
real,    parameter :: comp             = 0.001  ! m^-1

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

! ==== types =================================================================
type :: lake_pars_type
  real w_sat
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  real alpha              ! *** REPLACE LATER BY alpha(layer)
  real heat_capacity_ref
  real thermal_cond_ref
  real refl_dry_dir(NBANDS)
  real refl_dry_dif(NBANDS)
  real refl_sat_dir(NBANDS)
  real refl_sat_dif(NBANDS)
  real emis_dry
  real emis_sat
  real z0_momentum
  real depth_sill
  real width_sill
  real whole_area
  real connected_to_next
  real tau_groundwater
  real rsa_exp         ! riparian source-area exponent
end type lake_pars_type

type :: lake_prog_type
  real dz
  real wl
  real ws
  real T
  real K_z
  real groundwater
  real groundwater_T
end type lake_prog_type

type :: lake_tile_type
   integer :: tag ! kind of the lake
   type(lake_prog_type), pointer :: prog(:)
   type(lake_pars_type)          :: pars
   real,                 pointer :: w_fc(:)
   real,                 pointer :: w_wilt(:)
   real :: Eg_part_ref
   real :: z0_scalar
   real, pointer :: e(:),f(:)
   real, pointer :: heat_capacity_dry(:)
end type lake_tile_type

! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

!---- namelist ---------------------------------------------------------------
real    :: lake_width_inside_lake = 1.e5
real    :: large_lake_sill_width = 200.
real    :: min_lake_frac         = 0.
real    :: max_lake_rh           = 1.
real    :: k_over_B              = 0.25      ! reset to 0 for MCM
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
integer, public :: num_l                 = 18           ! number of lake levels
real    :: dz(max_lev)           = (/ &
    0.02, 0.04, 0.04, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, &
    0.2,   0.4,  0.4,  0.4,  0.4, 0.4,  1.,  1.,  1., &
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., &
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., &
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., &
    0.,0.,0.,0.,0. /)
                                              ! thickness (m) of model layers,
                                              ! from top down
logical :: use_lm2_awc           = .false.
  integer :: n_map_1st_lake_type = 10

! from analysis of modis data (ignoring temperature dependence):
  real :: f_iso_ice(NBANDS) = (/ 0.056, 0.131 /)
  real :: f_vol_ice(NBANDS) = (/ 0.017, 0.053 /)
  real :: f_geo_ice(NBANDS) = (/ 0.004, 0.010 /)
  real :: f_iso_liq(NBANDS) = (/ 0.056, 0.131 /)
  real :: f_vol_liq(NBANDS) = (/ 0.017, 0.053 /)
  real :: f_geo_liq(NBANDS) = (/ 0.004, 0.010 /)
  real :: refl_ice_dif(NBANDS), refl_liq_dif(NBANDS)

! ---- remainder are used only for cold start ---------
logical :: round_frac_down       = .false.  ! when false, any lake_frac < min_lake_frac
                                            ! is set to min_lake_frac.
                                            ! when true, any lake_frac < min_lake_frac
                                            ! is set to 0.
                                            
character(len=16):: lake_to_use     = 'single-tile'
       ! 'multi-tile' for tiled soil [default]
       ! 'single-tile' for geographically varying soil with single type per
       !     model grid cell
       ! 'uniform' for global constant soil, e.g., to reproduce MCM
       ! 'from-rivers' to get the fraction of lakes from river module
logical :: use_single_lake       = .false.   ! true for single global lake,
                                             ! e.g., to recover MCM
logical :: use_mcm_albedo        = .false.   ! .true. for CLIMAP albedo inputs
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
integer :: lake_index_constant   = 1         ! index of global constant lake,
                                             ! used when use_single_lake
real    :: gw_res_time           = 60.*86400 ! mean groundwater residence time,
                                             ! used when use_single_geo
real    :: rsa_exp_global        = 1.5
real, dimension(n_dim_lake_types) :: &
  dat_w_sat             =(/ 1.000   /),&
  dat_awc_lm2           =(/ 1.000   /),&
  dat_k_sat_ref         =(/ 0.021   /),&
  dat_psi_sat_ref       =(/ -.059   /),&
  dat_chb               =(/   3.5   /),&
  dat_heat_capacity_ref =(/ 8.4e7   /),&
  dat_thermal_cond_ref  =(/ 8.4e7   /),&
  dat_emis_dry          =(/ 0.950   /),&
  dat_emis_sat          =(/ 0.980   /),&
  dat_z0_momentum       =(/ 1.4e-4  /),&
  dat_tf_depr           =(/  0.00   /)
real, dimension(n_dim_lake_types, NBANDS) :: &
     dat_refl_dry_dif, dat_refl_dry_dir, &
     dat_refl_sat_dif, dat_refl_sat_dir
data &
                   !  VIS    NIR
  dat_refl_dry_dif / 0.060, 0.060 /, &
  dat_refl_dry_dir / 0.060, 0.060 /, &
  dat_refl_sat_dir / 0.060, 0.060 /, &
  dat_refl_sat_dif / 0.060, 0.060 /
integer, dimension(n_dim_lake_types) :: &
  input_cover_types     =(/ 10 /)
character(len=4), dimension(n_dim_lake_types) :: &
  tile_names            =(/ 'lake' /)

namelist /lake_data_nml/ lake_width_inside_lake, &
     large_lake_sill_width, &
     min_lake_frac, round_frac_down, max_lake_rh, &
     lake_to_use,input_cover_types, tile_names, &
     k_over_B,         &
     rate_fc, sfc_heat_factor,        &
     num_l,                   dz,                      &
     use_lm2_awc,    n_map_1st_lake_type, &
     use_single_lake,           use_mcm_albedo,            &
     use_single_geo,            lake_index_constant,         &
     gw_res_time,            rsa_exp_global,      &
     dat_w_sat,               dat_awc_lm2,     &
     dat_k_sat_ref,            &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_ref,         dat_thermal_cond_ref,   &
     dat_refl_dry_dir,            dat_refl_sat_dir,              &
     dat_refl_dry_dif,            dat_refl_sat_dif,              &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,           dat_tf_depr, &
     f_iso_ice, f_vol_ice, f_geo_ice, f_iso_liq, f_vol_liq, f_geo_liq 


!---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_lake_data_namelist(lake_n_lev)
  integer, intent(out) :: lake_n_lev
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  real    :: z

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=lake_data_nml, iostat=io)
     ierr = check_nml_error(io, 'lake_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=lake_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'lake_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=lake_data_nml)

  ! initialize global module data here

  refl_ice_dif = g_iso*f_iso_ice + g_vol*f_vol_ice + g_geo*f_geo_ice
  refl_liq_dif = g_iso*f_iso_liq + g_vol*f_vol_liq + g_geo*f_geo_liq

  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_lake_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_LAKE, idata1 = i )
  enddo

  ! set up output arguments
  lake_n_lev = num_l
end subroutine 


! ============================================================================
function lake_tile_ctor(tag) result(ptr)
  type(lake_tile_type), pointer :: ptr ! return value
  integer, intent(in)           :: tag ! kind of lake

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  allocate(ptr%prog   (num_l),  &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry(num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )
  call init_lake_data_0d(ptr)

end function lake_tile_ctor


! ============================================================================
function lake_tile_copy_ctor(lake) result(ptr)
  type(lake_tile_type), pointer    :: ptr  ! return value
  type(lake_tile_type), intent(in) :: lake ! tile to copy
  
  allocate(ptr)
  ptr=lake ! copy all non-pointer data
  ! allocate storage for tile data
  allocate(ptr%prog   (num_l),  &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry(num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )
  ! copy all allocatable data
  ptr%prog(:)   = lake%prog(:)
  ptr%w_fc(:)   = lake%w_fc(:)
  ptr%w_wilt(:) = lake%w_wilt(:)
  ptr%e(:)      = lake%e(:)
  ptr%f(:)      = lake%f(:)
  ptr%heat_capacity_dry(:) = lake%heat_capacity_dry(:)
end function lake_tile_copy_ctor


! ============================================================================
subroutine delete_lake_tile(ptr)
  type(lake_tile_type), pointer :: ptr

  deallocate(ptr%prog, ptr%w_fc, ptr%w_wilt,ptr%heat_capacity_dry, ptr%e, ptr%f)

  deallocate(ptr)
end subroutine delete_lake_tile


subroutine init_lake_data_0d(lake)
  type(lake_tile_type), intent(inout) :: lake

!  real tau_groundwater
!  real rsa_exp         ! riparian source-area exponent

  integer :: k
  k = lake%tag

  lake%pars%w_sat             = dat_w_sat            (k)
  lake%pars%awc_lm2           = dat_awc_lm2          (k)
  lake%pars%k_sat_ref         = dat_k_sat_ref        (k)
  lake%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  lake%pars%chb               = dat_chb              (k)
  lake%pars%alpha             = 1
  lake%pars%heat_capacity_ref = dat_heat_capacity_ref(k)
  lake%pars%thermal_cond_ref  = dat_thermal_cond_ref (k)
  lake%pars%refl_dry_dir      = dat_refl_dry_dir     (k,:)
  lake%pars%refl_dry_dif      = dat_refl_dry_dif     (k,:)
  lake%pars%refl_sat_dir      = dat_refl_sat_dir     (k,:)
  lake%pars%refl_sat_dif      = dat_refl_sat_dif     (k,:)
  lake%pars%emis_dry          = dat_emis_dry         (k)
  lake%pars%emis_sat          = dat_emis_sat         (k)
  lake%pars%z0_momentum       = dat_z0_momentum      (k)

  lake%pars%tau_groundwater   = 86400.*30.
  lake%pars%rsa_exp           = rsa_exp_global

  ! -------- derived constant lake parameters --------
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     lake%w_wilt(:) = 0.15
     lake%w_fc  (:) = 0.15 + lake%pars%awc_lm2
  else
     lake%w_wilt(:) = lake%pars%w_sat &
          *(lake%pars%psi_sat_ref/(psi_wilt*lake%pars%alpha))**(1/lake%pars%chb)
     lake%w_fc  (:) = lake%pars%w_sat &
          *(rate_fc/(lake%pars%k_sat_ref*lake%pars%alpha**2))**(1/(3+2*lake%pars%chb))
  endif

  ! below made use of phi_e from parlange via entekhabi
  lake%Eg_part_ref = (-4*lake%w_fc(1)**2*lake%pars%k_sat_ref*lake%pars%psi_sat_ref*lake%pars%chb &
       /(pi*lake%pars%w_sat)) * (lake%w_fc(1)/lake%pars%w_sat)**(2+lake%pars%chb)   &
       *(2*pi/(3*lake%pars%chb**2*(1+3/lake%pars%chb)*(1+4/lake%pars%chb)))/2

  lake%z0_scalar = lake%pars%z0_momentum * exp(-k_over_B)

end subroutine 


! ============================================================================
function lake_cover_cold_start(land_mask, lonb, latb, domain) result (lake_frac)
! creates and initializes a field of fractional lake coverage
  logical, intent(in) :: land_mask(:,:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:)! boundaries of the grid cells
  real,    pointer    :: lake_frac (:,:,:) ! output: map of lake fractional coverage
  type(domain2d), intent(in) :: domain

  allocate( lake_frac(size(land_mask,1),size(land_mask,2),n_dim_lake_types))

  if (trim(lake_to_use)=='from-rivers') then
     lake_frac = 0.0
     call read_data('INPUT/river_data.nc', 'lake_frac', lake_frac(:,:,1), &
          domain=domain)
     ! make sure 'missing values' don't get into the result
     where (lake_frac < 0) lake_frac = 0
     where (lake_frac > 1) lake_frac = 1
  else
     call init_cover_field(lake_to_use, 'INPUT/ground_type.nc', 'cover','frac', &
          lonb, latb, lake_index_constant, input_cover_types, lake_frac)
  endif
  
  if (round_frac_down) then
      where (lake_frac.gt.0. .and. lake_frac.lt.min_lake_frac) lake_frac = 0.
    else
      where (lake_frac.gt.0. .and. lake_frac.lt.min_lake_frac) lake_frac = min_lake_frac
    endif
  
end function 

! =============================================================================
function lake_tiles_can_be_merged(lake1,lake2) result(response)
  logical :: response
  type(lake_tile_type), intent(in) :: lake1,lake2

  response = (lake1%tag==lake2%tag)
end function

! =============================================================================
! combine two lake tiles with specified weights; the results goes into the 
! second one
! THIS NEEDS TO BE REVISED FOR TILE-DEPENDENT DZ
subroutine merge_lake_tiles(t1,w1,t2,w2)
  type(lake_tile_type), intent(in)    :: t1
  type(lake_tile_type), intent(inout) :: t2
  real, intent(in) :: w1, w2 ! relative weights

  ! ---- local vars
  real    :: x1, x2 ! normalized relative weights
  real    :: HEAT1, HEAT2 ! temporaries for heat
  real    :: C1, C2 ! heat capacities
  real    :: gw
  integer :: i

WRITE (*,*) 'SORRY, BUT merge_lake_tiles NEEDS TO BE REVISED TO ALLOW FOR ', &
            'HORIZONTALLY VARYING VERTICAL DISCRETIZATION'

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! combine state variables
  do i = 1,num_l
    ! calculate "dry" heat capacities:
    C1 = sfc_heat_factor*t1%pars%heat_capacity_ref
    C2 = sfc_heat_factor*t2%pars%heat_capacity_ref
    ! calculate heat content at this level for both source tiles
    HEAT1 = &
    (C1*dz(i)+clw*t1%prog(i)%wl+csw*t1%prog(i)%ws) * (t1%prog(i)%T-tfreeze)
    HEAT2 = &
    (C2*dz(i)+clw*t2%prog(i)%wl+csw*t2%prog(i)%ws) * (t2%prog(i)%T-tfreeze)
    ! calculate (and assign) combined water mass
    t2%prog(i)%wl = t1%prog(i)%wl*x1 + t2%prog(i)%wl*x2
    t2%prog(i)%ws = t1%prog(i)%ws*x1 + t2%prog(i)%ws*x2
    ! if dry heat capacity of combined lake is to be changed, update it here
    ! ...
    ! calculate combined temperature, based on total heat content and combined
    ! heat capacity
    t2%prog(i)%T = (HEAT1*x1+HEAT2*x2) / &
      (C2*dz(i)+clw*t2%prog(i)%wl+csw*t2%prog(i)%ws) + tfreeze

    ! calculate combined groundwater content
    gw = t1%prog(i)%groundwater*x1 + t2%prog(i)%groundwater*x2
    ! calculate combined groundwater temperature
    if(gw/=0) then
       t2%prog(i)%groundwater_T = ( &
            t1%prog(i)%groundwater*x1*(t1%prog(i)%groundwater_T-tfreeze) + &
            t2%prog(i)%groundwater*x2*(t2%prog(i)%groundwater_T-tfreeze)   &
            ) / gw + Tfreeze
    else
       t2%prog(i)%groundwater_T = &
            x1*t1%prog(i)%groundwater_T + x2*t2%prog(i)%groundwater_T
    endif
    t2%prog(i)%groundwater = gw
  enddo

end subroutine

! =============================================================================
! returns true if tile fits the specified selector
function lake_is_selected(lake, sel)
  logical lake_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(lake_tile_type),      intent(in) :: lake

  lake_is_selected = (sel%idata1 == lake%tag)
end function


! ============================================================================
! returns tag of the tile
function get_lake_tile_tag(lake) result(tag)
  integer :: tag
  type(lake_tile_type), intent(in) :: lake
  
  tag = lake%tag
end function


! ============================================================================
! compute bare-lake albedos and bare-lake emissivity
subroutine lake_data_radiation ( lake, cosz, use_brdf, &
                                  lake_alb_dir, lake_alb_dif, lake_emis )
  type(lake_tile_type), intent(in)  :: lake
  real,                 intent(in)  :: cosz
  logical,              intent(in)  :: use_brdf
  real,                 intent(out) :: lake_alb_dir(:), lake_alb_dif(:), lake_emis

  ! ---- local vars
  real :: lake_sfc_vlc, blend
  real :: liq_value_dir(NBANDS), ice_value_dir(NBANDS)
  real :: liq_value_dif(NBANDS), ice_value_dif(NBANDS)
  real :: zenith_angle, zsq, zcu

  ! ---- radiation properties
  lake_sfc_vlc = lake%prog(1)%wl/lake%prog(1)%dz
  blend        = lake_sfc_vlc/lake%pars%w_sat
  if (use_brdf) then
      zenith_angle = acos(cosz)
      zsq = zenith_angle*zenith_angle
      zcu = zenith_angle*zsq
      liq_value_dir =  f_iso_liq*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                     + f_vol_liq*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                     + f_geo_liq*(g0_geo+g1_geo*zsq+g2_geo*zcu)
      ice_value_dir =  f_iso_ice*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                     + f_vol_ice*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                     + f_geo_ice*(g0_geo+g1_geo*zsq+g2_geo*zcu)
      liq_value_dif = refl_liq_dif
      ice_value_dif = refl_ice_dif
    else
      liq_value_dir = lake%pars%refl_sat_dir
      ice_value_dir = lake%pars%refl_dry_dir
      liq_value_dif = lake%pars%refl_sat_dif
      ice_value_dif = lake%pars%refl_dry_dif
    endif
  lake_alb_dir = ice_value_dir + blend*(liq_value_dir-ice_value_dir)
  lake_alb_dif = ice_value_dif + blend*(liq_value_dif-ice_value_dif)
  lake_emis = lake%pars%emis_dry   + blend*(lake%pars%emis_sat-lake%pars%emis_dry  )
end subroutine

! ============================================================================
! compute bare-lake roughness
subroutine lake_data_diffusion ( lake,lake_z0s, lake_z0m )
  type(lake_tile_type), intent(in)  :: lake
  real,                 intent(out) :: lake_z0s, lake_z0m

  ! ---- surface roughness
  lake_z0s = lake%z0_scalar
  lake_z0m = lake%pars%z0_momentum
end subroutine

! ============================================================================
! compute lake thermodynamic properties.
subroutine lake_data_thermodynamics ( lake_pars, lake_depth, &
     lake_rh, heat_capacity_dry, thermal_cond)
  type(lake_pars_type), intent(in)  :: lake_pars
  real,                 intent(in)  :: lake_depth
  real,                 intent(out) :: lake_rh
  real,                 intent(out) :: heat_capacity_dry(:)
  real,                 intent(out) :: thermal_cond(:)

  ! ---- local vars
  integer l

! ----------------------------------------------------------------------------

  lake_rh = min(max_lake_rh, max(lake_depth/lake_pars%depth_sill,0.))

  do l = 1, num_l
     heat_capacity_dry(l) = lake_pars%heat_capacity_ref
     thermal_cond(l)  = lake_pars%thermal_cond_ref
  enddo

end subroutine

! ============================================================================
subroutine lake_tile_stock_pe (lake, twd_liq, twd_sol  )
  type(lake_tile_type),  intent(in)    :: lake
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n
  
  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(lake%prog)
    twd_liq = twd_liq + lake%prog(n)%wl + lake%prog(n)%groundwater
    twd_sol = twd_sol + lake%prog(n)%ws
    enddo

end subroutine lake_tile_stock_pe


! ============================================================================
! returns lake tile heat content, J/m2
function lake_tile_heat (lake) result(heat) ; real heat
  type(lake_tile_type),  intent(in)  :: lake

  integer :: i

  heat = 0
  do i = 1, num_l
     heat = heat + &
          (lake%heat_capacity_dry(i)*dz(i) + clw*lake%prog(i)%wl + csw*lake%prog(i)%ws)&
                           *(lake%prog(i)%T-tfreeze) + &
          clw*lake%prog(i)%groundwater*(lake%prog(i)%groundwater_T-tfreeze) - &
          hlf*lake%prog(i)%ws
  enddo
end function

end module lake_tile_mod
