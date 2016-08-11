#include <fms_platform.h>

module glac_tile_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : &
     write_version_number, file_exist, check_nml_error, &
     close_file, stdlog
use constants_mod, only : &
     pi, tfreeze, hlf
use land_constants_mod, only : NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_GLAC, register_tile_selector

implicit none
private

! ==== public interfaces =====================================================
public :: glac_prog_type
public :: glac_pars_type
public :: glac_tile_type

public :: new_glac_tile, delete_glac_tile
public :: glac_tiles_can_be_merged, merge_glac_tiles
public :: glac_is_selected
public :: get_glac_tile_tag
public :: glac_tile_stock_pe
public :: glac_tile_heat

public :: read_glac_data_namelist
public :: glac_cover_cold_start

public :: glac_data_radiation
public :: glac_data_diffusion
public :: glac_data_thermodynamics
public :: glac_data_hydraulics

public :: max_lev
! =====end of public interfaces ==============================================
interface new_glac_tile
   module procedure glac_tile_ctor
   module procedure glac_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), parameter   :: &
     version     = '$Id: glac_tile.F90,v 20.0 2013/12/13 23:29:33 fms Exp $', &
     tagname     = '$Name: tikal $', &
     module_name = 'glac_tile_mod'

integer, parameter :: max_lev          = 30 ! max number of levels in glacier
integer, parameter :: n_dim_glac_types = 1  ! size of lookup table
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
type :: glac_pars_type
  real w_sat
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  real alpha              ! *** REPLACE LATER BY alpha(layer)
  real heat_capacity_ref
  real thermal_cond_ref
  real :: refl_max_dir(NBANDS), refl_max_dif(NBANDS) ! max reflectance for the direct and diffuse light
  real :: refl_min_dir(NBANDS), refl_min_dif(NBANDS) ! min reflectance for the direct and diffuse light
  real emis_dry
  real emis_sat
  real z0_momentum
  real tau_groundwater
  real rsa_exp         ! riparian source-area exponent
  real tfreeze
end type glac_pars_type

type :: glac_prog_type
  real wl
  real ws
  real T
  real groundwater
  real groundwater_T
end type glac_prog_type

type :: glac_tile_type
   integer :: tag ! kind of the glacier
   type(glac_pars_type)               :: pars
   type(glac_prog_type), pointer :: prog(:)
   real,                 pointer :: w_fc(:)
   real,                 pointer :: w_wilt(:)
   real :: Eg_part_ref
   real :: z0_scalar
   real :: geothermal_heat_flux

   real, pointer :: heat_capacity_dry(:)
   real, pointer :: e(:), f(:)
end type glac_tile_type
! NOTE: When adding or modifying fields of this types, don't forget to change
! the operations on tile (e.g. copy) accordingly
! ==== module data ===========================================================

!---- namelist ---------------------------------------------------------------
real    :: k_over_B              = 2         ! reset to 0 for MCM
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
integer :: z_sfc_layer           = 0
integer :: num_l                 = 18        ! number of glacier levels
real    :: dz(max_lev)           = (/ &
    0.02, 0.04, 0.04, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, &
    0.2,   0.4,  0.4,  0.4,  0.4, 0.4,  1.,  1.,  1., &
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
                                             ! thickness (m) of model layers,
                                             ! from top down
logical :: use_lm2_awc           = .false.
logical :: use_lad1_glac         = .false.
real    :: t_range               = 10.0

! from analysis of modis data (ignoring temperature dependence):
  real :: f_iso_cold(NBANDS) = (/ 0.177, 0.265 /)
  real :: f_vol_cold(NBANDS) = (/ 0.100, 0.126 /)
  real :: f_geo_cold(NBANDS) = (/ 0.027, 0.032 /)
  real :: f_iso_warm(NBANDS) = (/ 0.177, 0.265 /)
  real :: f_vol_warm(NBANDS) = (/ 0.100, 0.126 /)
  real :: f_geo_warm(NBANDS) = (/ 0.027, 0.032 /)
  real :: refl_cold_dif(NBANDS), refl_warm_dif(NBANDS)

! ---- remainder are used only for cold start ---------
character(len=16):: glac_to_use  = 'single-tile'
       ! 'multi-tile' for tiled soil [default]
       ! 'single-tile' for geographically varying glacier with single type per
       !     model grid cell
       ! 'uniform' for global constant soil, e.g., to reproduce MCM
logical :: use_single_glac       = .false.   ! true for single global glacier,
                                             ! e.g., to recover MCM
logical :: use_mcm_albedo        = .false.   ! .true. for CLIMAP albedo inputs
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
integer :: glac_index_constant   = 1         ! index of global constant glacier,
                                             ! used when use_single_glacier
real    :: gw_res_time           = 60.*86400 ! mean groundwater residence time,
                                             ! used when use_single_geo
real    :: rsa_exp_global        = 1.5
real    :: geothermal_heat_flux_constant = 0.0  ! true continental average is ~0.065 W/m2
real, dimension(n_dim_glac_types) :: &
     dat_w_sat             =(/ 1.000   /),&
     dat_awc_lm2           =(/ 1.000   /),&
     dat_k_sat_ref         =(/ 0.021   /),&
     dat_psi_sat_ref       =(/ -.059   /),&
     dat_chb               =(/   3.5   /),&
     dat_heat_capacity_ref =(/ 1.6e6   /),&
     dat_thermal_cond_ref  =(/   1.8   /),&
     dat_emis_dry          =(/ 0.950   /),&
     dat_emis_sat          =(/ 0.980   /),&
     dat_z0_momentum       =(/  0.01   /),&
     dat_tf_depr           =(/  0.00   /)
real, dimension(n_dim_glac_types, NBANDS) :: &
     dat_refl_max_dir, dat_refl_max_dif, &
     dat_refl_min_dir, dat_refl_min_dif
                      !  VIS    NIR
data dat_refl_max_dir / 0.800, 0.800 /, & 
     dat_refl_max_dif / 0.800, 0.800 /, & 
     dat_refl_min_dir / 0.650, 0.650 /, &
     dat_refl_min_dif / 0.650, 0.650 /
integer, dimension(n_dim_glac_types) :: &
     input_cover_types     =(/ 9 /)
character(len=4), dimension(n_dim_glac_types) :: &
     tile_names            =(/'glac'/)
real, public :: &
     cpw = 1952.0, &  ! specific heat of water vapor at constant pressure
     clw = 4218.0, &  ! specific heat of water (liquid)
     csw = 2106.0     ! specific heat of water (ice)

namelist /glac_data_nml/ &
     glac_to_use, tile_names, input_cover_types, &
     k_over_B,  &
     rate_fc, sfc_heat_factor, z_sfc_layer,          &
     num_l, dz,             &
     use_lm2_awc,   use_lad1_glac, &
     use_single_glac,      use_mcm_albedo,            &
     use_single_geo,        glac_index_constant,         &
     gw_res_time,      rsa_exp_global,  geothermal_heat_flux_constant, &
     dat_w_sat,             dat_awc_lm2,     &
     dat_k_sat_ref,         &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_ref,         dat_thermal_cond_ref,   &
     dat_refl_max_dir,  dat_refl_max_dif,  &
     dat_refl_min_dir,  dat_refl_min_dif,  &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,           dat_tf_depr, &
     tile_names, input_cover_types, &
     f_iso_cold, f_vol_cold, f_geo_cold, f_iso_warm, f_vol_warm, f_geo_warm 

! ---- end of namelist

integer :: num_sfc_layers

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine read_glac_data_namelist(glac_n_lev, glac_dz)
  integer, intent(out) :: glac_n_lev
  real,    intent(out) :: glac_dz(:)
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  real    :: z

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=glac_data_nml, iostat=io)
     ierr = check_nml_error(io, 'glac_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=glac_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'glac_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write (unit, nml=glac_data_nml)

  refl_cold_dif = g_iso*f_iso_cold + g_vol*f_vol_cold + g_geo*f_geo_cold
  refl_warm_dif = g_iso*f_iso_warm + g_vol*f_vol_warm + g_geo*f_geo_warm

  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_glac_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_GLAC, idata1 = i )
  enddo

  ! initialize num_sfc_layers
  z = 0
  num_sfc_layers = 0
  do i = 1, num_l
     z = z + dz(i)
     if (z < z_sfc_layer+1.e-4) num_sfc_layers = i
  enddo

  ! set up output arguments
  glac_n_lev = num_l
  glac_dz = dz
end subroutine 


! ============================================================================
function glac_tile_ctor(tag) result(ptr)
  type(glac_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  allocate(ptr%prog   (num_l),  &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry (num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )

  ! set initial values of the tile data
  call glacier_data_init_0d(ptr)
end function glac_tile_ctor
  

! ============================================================================
function glac_tile_copy_ctor(glac) result(ptr)
  type(glac_tile_type), intent(in) :: glac ! return value
  type(glac_tile_type), pointer    :: ptr  ! return value

  allocate(ptr)
  ptr = glac ! copy all non-allocatable data
  ! allocate storage for tile data
  allocate(ptr%prog   (num_l),  &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry (num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )
  ! copy allocatable data
   ptr%prog(:)   = glac%prog(:)
   ptr%w_fc(:)   = glac%w_fc(:)
   ptr%w_wilt(:) = glac%w_wilt(:)
   ptr%heat_capacity_dry(:) = glac%heat_capacity_dry(:) 
   ptr%e(:)      = glac%e(:)
   ptr%f(:)      = glac%f(:)
 end function glac_tile_copy_ctor


! ============================================================================
subroutine delete_glac_tile(ptr)
  type(glac_tile_type), pointer :: ptr

  deallocate(ptr%prog, ptr%w_fc, ptr%w_wilt, ptr%heat_capacity_dry, &
        ptr%e,  ptr%f)
  deallocate(ptr)
end subroutine delete_glac_tile


! ============================================================================
subroutine glacier_data_init_0d(glac)
  type(glac_tile_type), intent(inout) :: glac

  integer :: k
  k = glac%tag

  glac%pars%w_sat             = dat_w_sat            (k)
  glac%pars%awc_lm2           = dat_awc_lm2          (k)
  glac%pars%k_sat_ref         = dat_k_sat_ref        (k)
  glac%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  glac%pars%chb               = dat_chb              (k)
  glac%pars%alpha             = 1
  glac%pars%heat_capacity_ref = dat_heat_capacity_ref(k)
  glac%pars%thermal_cond_ref  = dat_thermal_cond_ref (k)
  glac%pars%refl_max_dir      = dat_refl_max_dir     (k,:)
  glac%pars%refl_max_dif      = dat_refl_max_dif     (k,:)
  glac%pars%refl_min_dir      = dat_refl_min_dir     (k,:)
  glac%pars%refl_min_dif      = dat_refl_min_dif     (k,:)
  glac%pars%emis_dry          = dat_emis_dry         (k)
  glac%pars%emis_sat          = dat_emis_sat         (k)
  glac%pars%z0_momentum       = dat_z0_momentum      (k)
  glac%pars%tfreeze           = tfreeze - dat_tf_depr(k)

  glac%pars%tau_groundwater   = 86400*30 ! 30 days
  glac%pars%rsa_exp           = rsa_exp_global

  ! initialize derived data
  if (use_lm2_awc) then
     glac%w_wilt(:) = 0.15
     glac%w_fc  (:) = 0.15 + glac%pars%awc_lm2
  else
     glac%w_wilt(:) = glac%pars%w_sat &
          *(glac%pars%psi_sat_ref/(psi_wilt*glac%pars%alpha))**(1/glac%pars%chb)
     glac%w_fc  (:) = glac%pars%w_sat &
          *(rate_fc/(glac%pars%k_sat_ref*glac%pars%alpha**2))**(1/(3+2*glac%pars%chb))
  endif

  ! below made use of phi_e from parlange via entekhabi
  glac%Eg_part_ref  = (-4*glac%w_fc(1)**2*glac%pars%k_sat_ref*glac%pars%psi_sat_ref*glac%pars%chb &
       /(pi*glac%pars%w_sat)) * (glac%w_fc(1)/glac%pars%w_sat)**(2+glac%pars%chb)   &
       *(2*pi/(3*glac%pars%chb**2*(1+3/glac%pars%chb)*(1+4/glac%pars%chb)))/2

  glac%z0_scalar = glac%pars%z0_momentum * exp(-k_over_B)
  glac%geothermal_heat_flux = geothermal_heat_flux_constant
  
end subroutine glacier_data_init_0d


! ============================================================================
function glac_cover_cold_start(land_mask, lonb, latb) result (glac_frac)
  logical, intent(in) :: land_mask(:,:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  real,    pointer    :: glac_frac (:,:,:) ! output: map of fractional coverage

  allocate( glac_frac(size(land_mask,1),size(land_mask,2),n_dim_glac_types))

  call init_cover_field(glac_to_use, 'INPUT/ground_type.nc', 'cover','frac', &
       lonb, latb, glac_index_constant, input_cover_types, glac_frac)
  
end function 

! =============================================================================
function glac_tiles_can_be_merged(glac1,glac2) result(response)
  logical :: response
  type(glac_tile_type), intent(in) :: glac1,glac2

  response = (glac1%tag==glac2%tag)
end function

! =============================================================================
! combine two glacier states with specified weights; the result goes into
! the second one.
! t1 does not change; the properties of t2 (that is, heat capacity, etc) do not 
! change either -- that might be not good in the long run
subroutine merge_glac_tiles(g1,w1,g2,w2)
  type(glac_tile_type), intent(in)    :: g1
  type(glac_tile_type), intent(inout) :: g2     
  real                , intent(in)    :: w1, w2 ! relative weights
  
  ! ---- local vars
  real    :: x1, x2 ! normalized relative weights
  real    :: gw, HEAT1, HEAT2 ! temporaries for groundwater and heat
  integer :: i
  real    :: C1(num_l), C2(num_l) ! heat capacities of dry glacier
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1
  
  ! calculate dry heat capacities of the glaciers
  call glac_dry_heat_cap(g1,C1)
  call glac_dry_heat_cap(g2,C2)

  ! combine state variables
  do i = 1,num_l
    ! calculate heat content at this level for both source tiles
    HEAT1 = &
    (C1(i)*dz(i)+clw*g1%prog(i)%wl+csw*g1%prog(i)%ws) * (g1%prog(i)%T-tfreeze)
    HEAT2 = &
    (C2(i)*dz(i)+clw*g2%prog(i)%wl+csw*g2%prog(i)%ws) * (g2%prog(i)%T-tfreeze)
    ! calculate (and assign) combined water mass
    g2%prog(i)%wl = g1%prog(i)%wl*x1 + g2%prog(i)%wl*x2
    g2%prog(i)%ws = g1%prog(i)%ws*x1 + g2%prog(i)%ws*x2
    ! if dry heat capacity of combined glacier is to be changed, update it here
    ! ...
    ! calculate combined temperature, based on total heat content and combined
    ! heat capacity
    g2%prog(i)%T = (HEAT1*x1+HEAT2*x2) / &
      (C2(i)*dz(i)+clw*g2%prog(i)%wl+csw*g2%prog(i)%ws) + tfreeze

    ! calculate combined groundwater content
    gw = g1%prog(i)%groundwater*x1 + g2%prog(i)%groundwater*x2
    ! calculate combined groundwater temperature
    if(gw/=0) then
       g2%prog(i)%groundwater_T = ( &
            g1%prog(i)%groundwater*x1*(g1%prog(i)%groundwater_T-tfreeze) + &
            g2%prog(i)%groundwater*x2*(g2%prog(i)%groundwater_T-tfreeze)   &
            ) / gw + Tfreeze
    else
       g2%prog(i)%groundwater_T = &
            x1*g1%prog(i)%groundwater_T + x2*g2%prog(i)%groundwater_T
    endif
    g2%prog(i)%groundwater = gw
  enddo
  
end subroutine

! =============================================================================
! returns true if tile fits the specified selector
function glac_is_selected(glac, sel)
  logical glac_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(glac_tile_type),      intent(in) :: glac

  glac_is_selected = (sel%idata1 == glac%tag)
end function 


! ============================================================================
! returns tag of the tile
function get_glac_tile_tag(glac) result(tag)
  integer :: tag
  type(glac_tile_type), intent(in) :: glac
  
  tag = glac%tag
end function

! ============================================================================
! compute glacier thermodynamic properties.
subroutine glac_dry_heat_cap ( glac, heat_capacity_dry)
  type(glac_tile_type), intent(in)  :: glac
  real,                 intent(out) :: heat_capacity_dry(:)

  ! ---- local vars
  integer l

  ! these will eventually be functions of water content (or psi) and T.
  do l = 1, num_sfc_layers
     heat_capacity_dry(l) = sfc_heat_factor*glac%pars%heat_capacity_ref
  enddo
  do l = num_sfc_layers+1, num_l
     heat_capacity_dry(l) = glac%pars%heat_capacity_ref
  enddo

end subroutine 


! ============================================================================
! compute bare-glacier albedos and bare-glacier emissivity
subroutine glac_data_radiation ( glac, cosz, use_brdf, &
                                 glac_refl_dir, glac_refl_dif, glac_emis )
  type(glac_tile_type), intent(in) :: glac
  real,                 intent(in) :: cosz
  logical,              intent(in) :: use_brdf
  real,                 intent(out):: glac_refl_dir(:), glac_refl_dif(:), glac_emis

  ! ---- local vars 
  real  :: blend, t_crit
  real :: warm_value_dir(NBANDS), cold_value_dir(NBANDS)
  real :: warm_value_dif(NBANDS), cold_value_dif(NBANDS)
  real :: zenith_angle, zsq, zcu

  t_crit = tfreeze - t_range
  blend = (glac%prog(1)%T - t_crit) / t_range
  if (blend < 0.0) blend = 0.0
  if (blend > 1.0) blend = 1.0
  if (use_brdf) then
      zenith_angle = acos(cosz)
      zsq = zenith_angle*zenith_angle
      zcu = zenith_angle*zsq
      warm_value_dir = f_iso_warm*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                     + f_vol_warm*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                     + f_geo_warm*(g0_geo+g1_geo*zsq+g2_geo*zcu)
      cold_value_dir = f_iso_cold*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                     + f_vol_cold*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                     + f_geo_cold*(g0_geo+g1_geo*zsq+g2_geo*zcu)
      warm_value_dif = refl_warm_dif
      cold_value_dif = refl_cold_dif
    else
      warm_value_dir = glac%pars%refl_min_dir
      cold_value_dir = glac%pars%refl_max_dir
      warm_value_dif = glac%pars%refl_min_dif
      cold_value_dif = glac%pars%refl_max_dif
    endif
  glac_refl_dir = cold_value_dir + blend*(warm_value_dir-cold_value_dir)
  glac_refl_dif = cold_value_dif + blend*(warm_value_dif-cold_value_dif)
  glac_emis     = glac%pars%emis_dry + blend*(glac%pars%emis_sat-glac%pars%emis_dry)

end subroutine glac_data_radiation


! ============================================================================
! compute bare-glacier roughness
subroutine glac_data_diffusion ( glac, glac_z0s, glac_z0m )
  type(glac_tile_type), intent(in) :: glac
  real,                 intent(out):: glac_z0s, glac_z0m

  ! ---- surface roughness ---------------------------------------------------
  glac_z0s = glac%z0_scalar
  glac_z0m = glac%pars%z0_momentum

end subroutine glac_data_diffusion


! ============================================================================
! compute glacier thermodynamic properties.
subroutine glac_data_thermodynamics ( glac_pars, vlc_sfc, vsc_sfc, &  
     glac_rh, heat_capacity_dry, thermal_cond)
  type(glac_pars_type), intent(in)  :: glac_pars
  real,                 intent(in)  :: vlc_sfc
  real,                 intent(in)  :: vsc_sfc
  real,                 intent(out) :: glac_rh
  real,                 intent(out) :: heat_capacity_dry(:)
  real,                 intent(out) :: thermal_cond(:)

  ! ---- local vars
  integer l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  glac_rh = 1

  ! these will eventually be functions of water content (or psi) and T.
  do l = 1, num_sfc_layers
     heat_capacity_dry(l) = sfc_heat_factor*glac_pars%heat_capacity_ref
     thermal_cond(l)  = sfc_heat_factor*glac_pars%thermal_cond_ref
  enddo
  do l = num_sfc_layers+1, num_l
     heat_capacity_dry(l) = glac_pars%heat_capacity_ref
     thermal_cond(l)  = glac_pars%thermal_cond_ref
  enddo

end subroutine 


! ============================================================================
! compute glacier hydraulic properties.
subroutine glac_data_hydraulics (glac, vlc, vsc, &
     psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, tau_gw, &
     glac_w_fc  )
  type(glac_tile_type),        intent(in)  :: glac
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
       psi, DThDP, hyd_cond, DKDP, glac_w_fc
  real,                        intent(out) :: &
       DPsi_min, DPsi_max, tau_gw

  ! ---- local vars
  integer l
  real,  dimension(num_l) :: vlc_loc
  real :: real1
  logical :: logical1

  ! ---- T-dependence of hydraulic properties --------------------------------
  ! k_sat   = glac%pars%k_sat0   !  * mu(t0)/mu(t), where mu is dynamic viscosity
  ! psi_sat = glac%pars%psi_sat0 !  * exp(c*(psi-psi0)), where c~+/-(?)0.0068
  ! better approach would be to adopt air entrapment model
  ! or at least to scale against surface tension model


  ! ---- water and ice dependence of hydraulic properties --------------------
  ! ---- (T-dependence can be added later)
  hyd_cond = 1; DThDP = 1
  do l = 1, num_l
    hyd_cond(l) = (glac%pars%k_sat_ref*glac%pars%alpha**2)*  &
         ! * mu(T)/mu(t_ref), where mu is dynamic viscosity
          (vlc(l)/glac%pars%w_sat)**(3+2*glac%pars%chb)
    if (hyd_cond(l).lt.1.e-12*glac%pars%k_sat_ref) then
      vlc_loc (l) = glac%pars%w_sat*(1.e-12)**(1./(3+2*glac%pars%chb))
      hyd_cond(l) = 1.e-12*glac%pars%k_sat_ref
      if (vsc(l).eq.0.) then
        DThDP   (l) = -vlc_loc(l)  &
             *(vlc_loc(l)/glac%pars%w_sat)**glac%pars%chb &
             /(glac%pars%psi_sat_ref*glac%pars%chb)
        psi     (l) = (glac%pars%psi_sat_ref/glac%pars%alpha) &
             *(glac%pars%w_sat/vlc_loc(l))**glac%pars%chb &
             + (vlc(l)-vlc_loc (l))/DThDP   (l)
        DKDP    (l) = 0.
      else
        psi     (l) = ((glac%pars%psi_sat_ref/glac%pars%alpha) / 2.2) &
             *(glac%pars%w_sat/vlc_loc(l))**glac%pars%chb
        DKDP    (l) = 0.
        DThDP   (l) = 0.
      endif
    else
      if (vsc(l).eq.0.) then
        real1 = glac%pars%w_sat - vlc(l)
        logical1 = real1.ge.0.
        if (logical1) then
           !              where (vlc(l).le.glac%pars%w_sat)
           psi     (l) = (glac%pars%psi_sat_ref/glac%pars%alpha) &
                *(glac%pars%w_sat/vlc(l))**glac%pars%chb
           DKDP    (l) = -(2+3/glac%pars%chb)*hyd_cond(l) &
                /psi(l)
           DThDP   (l) = -vlc(l)/(psi(l)*glac%pars%chb)
        else
           psi(l) = glac%pars%psi_sat_ref &
                + (vlc(l)-glac%pars%w_sat)/comp
           DThDP(l) = comp
           hyd_cond(l) = glac%pars%k_sat_ref
           DKDP(l) = 0.
        endif
      else
        psi     (l) = ((glac%pars%psi_sat_ref/glac%pars%alpha) / 2.2) &
             *(glac%pars%w_sat/vlc(l))**glac%pars%chb
        DKDP    (l) = 0.
        DThDP   (l) = 0.
      endif
    endif
  enddo

  if (DThDP(1).ne.0.) then
    DPsi_min =            -vlc(1) /DThDP(1)
    DPsi_max = (glac%pars%w_sat-vlc(1))/DThDP(1)
  else
    Dpsi_min = -1.e16
    DPsi_max = -psi(1)
  endif

  glac_w_fc = glac%w_fc

  ! ---- groundwater ---------------------------------------------------------
  tau_gw = glac%pars%tau_groundwater

end subroutine glac_data_hydraulics

! ============================================================================
subroutine glac_tile_stock_pe (glac, twd_liq, twd_sol  )
  type(glac_tile_type),  intent(in)    :: glac
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n
  
  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(glac%prog)
    twd_liq = twd_liq + glac%prog(n)%wl + glac%prog(n)%groundwater
    twd_sol = twd_sol + glac%prog(n)%ws
    enddo

end subroutine glac_tile_stock_pe

! ============================================================================
! returns glacier heat content, J/m2
function glac_tile_heat (glac) result(heat) ; real heat
  type(glac_tile_type),  intent(in)  :: glac

  integer :: i

  heat = 0
  do i = 1, num_l
     heat = heat + &
          (glac%heat_capacity_dry(i)*dz(i) + clw*glac%prog(i)%wl + csw*glac%prog(i)%ws)&
                           *(glac%prog(i)%T-tfreeze) + &
          clw*glac%prog(i)%groundwater*(glac%prog(i)%groundwater_T-tfreeze) - &
          hlf*glac%prog(i)%ws
  enddo
end function

end module glac_tile_mod
