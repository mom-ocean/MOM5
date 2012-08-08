module cana_tile_mod

use land_tile_selectors_mod, only : &
     tile_selector_type
use constants_mod, only : &
     cp_air, tfreeze

implicit none
private

! ==== public interfaces =====================================================
public :: cana_prog_type
public :: cana_tile_type

public :: new_cana_tile, delete_cana_tile
public :: cana_tiles_can_be_merged, merge_cana_tiles
public :: get_cana_tile_tag
public :: cana_is_selected

public :: cana_tile_stock_pe
public :: cana_tile_carbon
public :: cana_tile_heat

! public data:
real, public :: canopy_air_mass = 0.0    ! mass of wet air in the canopy air 
                                         ! space for heat and water vapor, kg/m2
real, public :: canopy_air_mass_for_tracers = 0.0 ! mass of wet air in the canopy air 
                                         ! space for tracers other than water vapor, kg/m2
! Water vapor is bundled with heat and not with other tracers because it is
! tightly coupled with the heat capacity of the canopy air and therefore with
! the equations for heat. We assume that other tracers do not contribute to
! the heat capacity of the canopy air.
real, public :: cpw             = 1952.0 ! specific heat of water vapor at constant pressure, J/(kg K)
! ==== end of public interfaces ==============================================
interface new_cana_tile
   module procedure cana_tile_ctor
   module procedure cana_tile_copy_ctor
end interface

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: cana_tile.F90,v 18.0 2010/03/02 23:36:42 fms Exp $', &
     tagname = '$Name: siena_201207 $'

! ==== data types ======================================================
type :: cana_prog_type
  real T
  real q
  real :: co2 ! co2 concentration in canopy air, kg CO2/kg of wet air
end type cana_prog_type

type :: cana_tile_type
   type(cana_prog_type) :: prog
end type cana_tile_type

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! =============================================================================
function cana_tile_ctor() result(ptr)
  type(cana_tile_type), pointer :: ptr ! return value

  allocate(ptr)
end function cana_tile_ctor

! =============================================================================
function cana_tile_copy_ctor(cana) result(ptr)
  type(cana_tile_type), pointer :: ptr ! return value
  type(cana_tile_type), intent(in) :: cana ! return value

  allocate(ptr)
  ptr = cana
end function cana_tile_copy_ctor

! =============================================================================
subroutine delete_cana_tile(cana)
  type(cana_tile_type), pointer :: cana

  deallocate(cana)
end subroutine delete_cana_tile

! =============================================================================
function cana_tiles_can_be_merged(cana1,cana2) result(response)
  logical :: response
  type(cana_tile_type), intent(in) :: cana1,cana2

  response = .TRUE.
end function

! =============================================================================
subroutine merge_cana_tiles(cana1,w1,cana2,w2)
  type(cana_tile_type), intent(in)    :: cana1
  type(cana_tile_type), intent(inout) :: cana2
  real                , intent(in)    :: w1, w2
  
  ! ---- local vars
  real :: x1,x2 ! normalized weights
  real :: HEAT1, HEAT2 ! heat content of the tiles
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1-x1
  HEAT1 = canopy_air_mass*(cp_air+(cpw - cp_air)*cana1%prog%q)*cana1%prog%T
  HEAT2 = canopy_air_mass*(cp_air+(cpw - cp_air)*cana2%prog%q)*cana2%prog%T

  cana2%prog%q = cana1%prog%q*x1+cana2%prog%q*x2
  if (canopy_air_mass > 0) then
     cana2%prog%T = (HEAT1*x1+HEAT2*x2)/&
          (canopy_air_mass*(cp_air+(cpw - cp_air)*cana2%prog%q))
  else
     cana2%prog%T = cana1%prog%T*x1+cana2%prog%T*x2
  endif

  cana2%prog%co2 = cana1%prog%co2*x1+cana2%prog%co2*x2
end subroutine

! =============================================================================
! returns tag of the tile
function get_cana_tile_tag(cana) result(tag)
  integer :: tag
  type(cana_tile_type), intent(in) :: cana
  
  tag = 1
end function

! =============================================================================
! returns true if tile fits the specified selector
function cana_is_selected (cana, sel)
  logical cana_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(cana_tile_type),      intent(in) :: cana

  cana_is_selected = .TRUE.
end function

! =============================================================================
subroutine cana_tile_stock_pe (cana, twd_liq, twd_sol)
  type(cana_tile_type), intent(in) :: cana
  real, intent(out) :: twd_liq, twd_sol

  twd_liq = canopy_air_mass*cana%prog%q; twd_sol = 0
end subroutine

! =============================================================================
function cana_tile_heat (cana) result(heat) ; real heat
  type(cana_tile_type), intent(in) :: cana
  
  heat = canopy_air_mass*(cp_air+(cpw - cp_air)*cana%prog%q)*(cana%prog%T-tfreeze)
end function

! =============================================================================
function cana_tile_carbon (cana) result(c) ; real c
  type(cana_tile_type), intent(in) :: cana

  c = canopy_air_mass_for_tracers * cana%prog%co2
end function 

end module cana_tile_mod
