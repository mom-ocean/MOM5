#include <fms_platform.h>

module snow_tile_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : &
     write_version_number, file_exist, check_nml_error, &
     close_file, stdlog
use constants_mod,only: tfreeze, hlf
use land_constants_mod, only : &
     NBANDS
use land_tile_selectors_mod, only : &
     tile_selector_type

implicit none
private

! ==== public interfaces =====================================================
public :: snow_prog_type
public :: snow_tile_type

public :: new_snow_tile, delete_snow_tile
public :: snow_tiles_can_be_merged, merge_snow_tiles
public :: snow_is_selected
public :: get_snow_tile_tag
public :: snow_tile_stock_pe
public :: snow_tile_heat

public :: read_snow_data_namelist

public :: snow_data_thermodynamics
public :: snow_data_hydraulics
public :: snow_data_area
public :: snow_data_radiation
public :: snow_data_diffusion

public :: max_lev
public :: cpw, clw, csw
! ==== end of public interfaces ==============================================
interface new_snow_tile
   module procedure snow_tile_ctor
   module procedure snow_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'snow_tile_mod' ,&
     version     = '$Id: snow_tile.F90,v 19.0 2012/01/06 20:42:44 fms Exp $' ,&
     tagname     = '$Name: siena_201207 $'
integer, parameter :: max_lev = 10
real   , parameter :: t_range = 10.0 ! degK

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
type :: snow_prog_type
  real wl
  real ws
  real T
end type snow_prog_type


type :: snow_tile_type
   integer :: tag ! kind of the tile
   type(snow_prog_type), pointer :: prog(:)
   real,                 pointer :: e(:), f(:)
end type snow_tile_type

! ==== module data ===========================================================

!---- namelist ---------------------------------------------------------------
logical :: use_mcm_masking       = .false.   ! MCM snow mask fn
real    :: w_sat                 = 670.
real    :: psi_sat               = -0.06
real    :: k_sat                 = 0.02
real    :: chb                   = 3.5
real    :: thermal_cond_ref      = 0.3
real    :: depth_crit            = 0.0167
real    :: z0_momentum           = 0.001
real    :: refl_snow_max_dir(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real    :: refl_snow_max_dif(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real    :: refl_snow_min_dir(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real    :: refl_snow_min_dif(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real    :: emis_snow_max         = 0.95      ! reset to 1 for MCM
real    :: emis_snow_min         = 0.90      ! reset to 1 for MCM
real    :: k_over_B              = 2         ! reset to 0 for MCM
integer :: num_l                 = 3         ! number of snow levels
real    :: dz(max_lev)           = (/0.1,0.8,0.1,0.,0.,0.,0.,0.,0.,0./)
                                              ! rel. thickness of model layers,
                                              ! from top down
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
real    :: mc_fict = 10. * 4218 ! additional (fictitious) soil heat capacity (for numerical stability?).
! from analysis of modis data (ignoring temperature dependence):
  real :: f_iso_cold(NBANDS) = (/ 0.354, 0.530 /)
  real :: f_vol_cold(NBANDS) = (/ 0.200, 0.252 /)
  real :: f_geo_cold(NBANDS) = (/ 0.054, 0.064 /)
  real :: f_iso_warm(NBANDS) = (/ 0.354, 0.530 /)
  real :: f_vol_warm(NBANDS) = (/ 0.200, 0.252 /)
  real :: f_geo_warm(NBANDS) = (/ 0.054, 0.064 /)
  real :: refl_cold_dif(NBANDS), refl_warm_dif(NBANDS)

namelist /snow_data_nml/use_mcm_masking,    w_sat,                 &
                    psi_sat,                k_sat,                 &
                    chb,                                           &
                    thermal_cond_ref,       depth_crit,            &
                    z0_momentum,                                   &
                    refl_snow_max_dir,    refl_snow_min_dir,   &
                    refl_snow_max_dif,    refl_snow_min_dif,   &
                    emis_snow_max,          emis_snow_min,         &
                    k_over_B,             &
                    num_l,                   dz, cpw, clw, csw, mc_fict, &
     f_iso_cold, f_vol_cold, f_geo_cold, f_iso_warm, f_vol_warm, f_geo_warm 
     
!---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_snow_data_namelist(snow_num_l, snow_dz, snow_mc_fict)
  integer, intent(out) :: snow_num_l
  real,    intent(out) :: snow_dz(:)
  real,    intent(out) :: snow_mc_fict

  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=snow_data_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=snow_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=snow_data_nml)

  ! initialize global module data here

  ! set up output arguments
  snow_num_l = num_l
  snow_dz    = dz
  snow_mc_fict = mc_fict

  refl_cold_dif = g_iso*f_iso_cold + g_vol*f_vol_cold + g_geo*f_geo_cold
  refl_warm_dif = g_iso*f_iso_warm + g_vol*f_vol_warm + g_geo*f_geo_warm

end subroutine 

! ============================================================================
function snow_tile_ctor(tag) result(ptr)
  type(snow_tile_type), pointer :: ptr ! return value
  integer, optional, intent(in) :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = 0 ; if(present(tag)) ptr%tag = tag
  ! allocate storage for tile data
  allocate(ptr%prog(num_l))
  allocate(ptr%e(num_l))
  allocate(ptr%f(num_l))

end function snow_tile_ctor

! ============================================================================
function snow_tile_copy_ctor(snow) result(ptr)
  type(snow_tile_type), pointer :: ptr ! return value
  type(snow_tile_type), intent(in) :: snow ! tile to copy

  allocate(ptr)
  ! copy all non-pointer members
  ptr = snow
  ! allocate storage for tile data
  allocate(ptr%prog(num_l))
  allocate(ptr%e(num_l))
  allocate(ptr%f(num_l))
  ! copy all pointer members
  ptr%prog(:) = snow%prog(:)
  ptr%e(:) = snow%e(:)
  ptr%f(:) = snow%f(:)
end function snow_tile_copy_ctor

! ============================================================================
subroutine delete_snow_tile(snow)
  type(snow_tile_type), pointer :: snow

  deallocate(snow%prog)
  deallocate(snow%e)
  deallocate(snow%f)
  deallocate(snow)
end subroutine delete_snow_tile

! =============================================================================
function snow_tiles_can_be_merged(snow1,snow2) result(response)
  logical :: response
  type(snow_tile_type), intent(in) :: snow1,snow2

  response = .TRUE.
end function

! =============================================================================
subroutine merge_snow_tiles(snow1, w1, snow2, w2)
  type(snow_tile_type), intent(in)    :: snow1
  type(snow_tile_type), intent(inout) :: snow2
  real                , intent(in)    :: w1, w2 ! relative weights
  
  ! ---- local vars
  real    :: x1, x2 ! normalized weights
  real    :: HEAT1, HEAT2
  integer :: i
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1-x1
  
  do i = 1, num_l
    HEAT1 = (mc_fict*dz(i)+clw*snow1%prog(i)%wl+csw*snow1%prog(i)%ws)*(snow1%prog(i)%T-tfreeze)
    HEAT2 = (mc_fict*dz(i)+clw*snow2%prog(i)%wl+csw*snow2%prog(i)%ws)*(snow2%prog(i)%T-tfreeze)
    snow2%prog(i)%wl = snow1%prog(i)%wl*x1 + snow2%prog(i)%wl*x2
    snow2%prog(i)%ws = snow1%prog(i)%ws*x1 + snow2%prog(i)%ws*x2
    if (snow2%prog(i)%wl/=0.or.snow2%prog(i)%ws/=0) then
       snow2%prog(i)%T  = (HEAT1*x1+HEAT2*x2)/&
            (mc_fict*dz(i)+clw*snow2%prog(i)%wl+csw*snow2%prog(i)%ws)+tfreeze
    else
       snow2%prog(i)%T  = snow1%prog(i)%T*x1 + snow2%prog(i)%T*x2
    endif
  enddo
end subroutine

! =============================================================================
! returns true if tile fits the specified selector
function snow_is_selected(snow, sel)
  logical snow_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(snow_tile_type),      intent(in) :: snow

  snow_is_selected = .TRUE.
end function

! ============================================================================
! retruns tag of the tile
function get_snow_tile_tag(snow) result(tag)
  integer :: tag
  type(snow_tile_type), intent(in) :: snow
  
  tag = snow%tag
end function

! ============================================================================
! compute snow thermodynmamic properties.
subroutine snow_data_thermodynamics ( snow_rh, thermal_cond)
  real, intent(out) :: snow_rh
  real, intent(out) :: thermal_cond(:)

  ! snow surface assumed to have air at saturation
  snow_rh = 1

  ! these will eventually be functions of water contents and T.
  thermal_cond  = thermal_cond_ref

end subroutine 


! ============================================================================
! compute snow hydraulic properties (assumed dependent only on wl)
subroutine snow_data_hydraulics (wl, ws, psi, hyd_cond )
  real, intent(in),  dimension(:) :: wl, ws
  real, intent(out), dimension(:) :: psi, hyd_cond

  ! ---- local vars 
  integer :: l
  
  do l = 1, num_l
    psi     (l) = psi_sat *(w_sat/(wl(l)+ws(l)))**chb
    hyd_cond(l) = k_sat*(wl(l)/w_sat)**(3+2*chb)
  enddo

end subroutine snow_data_hydraulics


! ============================================================================
! compute snow area
subroutine snow_data_area ( snow_depth, snow_area )
    real, intent(in)  :: snow_depth
    real, intent(out) :: snow_area

  snow_area = 0.
  if (use_mcm_masking) then
     snow_area = min(1., 0.5*sqrt(max(0.,snow_depth)/depth_crit))
  else
     snow_area = max(0.,snow_depth) / (max(0.,snow_depth) + depth_crit)
  endif

end subroutine


! ============================================================================
subroutine snow_data_radiation(snow_T, snow_refl_dir, snow_refl_dif, snow_emis,&
                                  cosz, use_brdf)
  real, intent(in) :: snow_T
  real, intent(out):: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_emis
  real, optional :: cosz
  logical :: use_brdf

  ! ---- local vars
  real :: blend
  real :: warm_value_dir(NBANDS), cold_value_dir(NBANDS)
  real :: warm_value_dif(NBANDS), cold_value_dif(NBANDS)
  real :: zenith_angle, zsq, zcu

  blend = max(0.,min(1.,1.-(tfreeze-snow_T)/t_range))
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
     warm_value_dir = refl_snow_min_dir
     cold_value_dir = refl_snow_max_dir
     warm_value_dif = refl_snow_min_dif
     cold_value_dif = refl_snow_max_dif
  endif
  snow_refl_dir = cold_value_dir + blend*(warm_value_dir-cold_value_dir)
  snow_refl_dif = cold_value_dif + blend*(warm_value_dif-cold_value_dif)
  snow_emis     = emis_snow_max + blend*(emis_snow_min-emis_snow_max  )
end subroutine


! ============================================================================
subroutine snow_data_diffusion(snow_z0s, snow_z0m)
  real, intent(out):: snow_z0s, snow_z0m

  snow_z0m =  z0_momentum
  snow_z0s =  z0_momentum * exp(-k_over_B)
end subroutine

! ============================================================================
subroutine snow_tile_stock_pe (snow, twd_liq, twd_sol  )
  type(snow_tile_type),  intent(in)    :: snow
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n
  
  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(snow%prog)
    twd_liq = twd_liq + snow%prog(n)%wl
    twd_sol = twd_sol + snow%prog(n)%ws
    enddo

end subroutine snow_tile_stock_pe

! ============================================================================
! returns snow heat content, J/m2
function snow_tile_heat (snow) result(heat) ; real heat
  type(snow_tile_type), intent(in)  :: snow

  integer :: i

  heat = 0
  do i = 1,num_l
     heat = heat - snow%prog(i)%ws*hlf &
        + (mc_fict*dz(i) + clw*snow%prog(i)%wl + csw*snow%prog(i)%ws)  &
                                      * (snow%prog(i)%T-tfreeze) 
  enddo
end function

end module snow_tile_mod
