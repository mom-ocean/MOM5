module land_tile_mod

use fms_mod, only : error_mesg, FATAL

use land_constants_mod, only : NBANDS
use glac_tile_mod, only : &
     glac_tile_type, new_glac_tile, delete_glac_tile, glac_is_selected, &
     glac_tiles_can_be_merged, merge_glac_tiles, get_glac_tile_tag, &
     glac_tile_stock_pe, glac_tile_heat
use lake_tile_mod, only : &
     lake_tile_type, new_lake_tile, delete_lake_tile, lake_is_selected, &
     lake_tiles_can_be_merged, merge_lake_tiles, get_lake_tile_tag, &
     lake_tile_stock_pe, lake_tile_heat
use soil_tile_mod, only : &
     soil_tile_type, new_soil_tile, delete_soil_tile, soil_is_selected, &
     soil_tiles_can_be_merged, merge_soil_tiles, get_soil_tile_tag, &
     soil_tile_stock_pe, soil_tile_carbon, soil_tile_heat
use cana_tile_mod, only : &
     cana_tile_type, new_cana_tile, delete_cana_tile, cana_is_selected, &
     cana_tiles_can_be_merged, merge_cana_tiles, get_cana_tile_tag, &
     cana_tile_stock_pe, cana_tile_carbon, cana_tile_heat
use vegn_tile_mod, only : &
     vegn_tile_type, new_vegn_tile, delete_vegn_tile, vegn_is_selected, &
     vegn_tiles_can_be_merged, merge_vegn_tiles, get_vegn_tile_tag, &
     vegn_tile_stock_pe, vegn_tile_carbon, vegn_tile_heat
use snow_tile_mod, only : &
     snow_tile_type, new_snow_tile, delete_snow_tile, snow_is_selected, &
     snow_tiles_can_be_merged, merge_snow_tiles, get_snow_tile_tag, &
     snow_tile_stock_pe, snow_tile_heat
use land_tile_selectors_mod, only : tile_selector_type, &
     SEL_SOIL, SEL_VEGN, SEL_LAKE, SEL_GLAC, SEL_SNOW, SEL_CANA
use tile_diag_buff_mod, only : &
     diag_buff_type, new_diag_buff, delete_diag_buff

implicit none
private
! ==== public interfaces =====================================================
public :: land_tile_type
public :: land_tile_list_type
public :: land_tile_enum_type
public :: diag_buff_type

! operations with tile
public :: new_land_tile, delete_land_tile
public :: land_tiles_can_be_merged, merge_land_tiles

public :: get_tile_tags ! returns the tags of the sub-model tiles
public :: get_tile_water ! returns liquid and frozen water masses
public :: land_tile_carbon ! returns total carbon in the tile
public :: land_tile_heat ! returns tile heat content

! operations with tile lists and tile list enumerators
public :: land_tile_list_init, land_tile_list_end
public :: first_elmt, tail_elmt
public :: operator(==), operator(/=) ! comparison of two enumerators
public :: next_elmt, prev_elmt ! enumerator advance operations
public :: current_tile ! returns pointer to the tile at a position
public :: insert  ! inserts a tile at a given position, or appends it to a list
public :: erase   ! erases tile at current position
public :: remove  ! removes tile at current position, but does not delete it
public :: get_elmt_indices ! returns i,j,k of current element

public :: empty   ! returns true if the list of tiles is empty
public :: nitems  ! count of items in list

public :: tile_is_selected

public :: print_land_tile_info
public :: print_land_tile_statistics
! ==== end of public interfaces ==============================================
interface new_land_tile
   module procedure land_tile_ctor
   module procedure land_tile_copy_ctor
end interface

interface first_elmt
   module procedure land_tile_list_begin_0d
   module procedure land_tile_list_begin_2d
end interface
interface tail_elmt
   module procedure land_tile_list_end_0d
   module procedure land_tile_list_end_2d
end interface

interface operator(==)
   module procedure enums_are_equal
end interface
interface operator(/=)
   module procedure enums_are_not_equal
end interface

interface insert
   module procedure insert_at_position, insert_in_list
end interface
interface remove
   module procedure remove_at_position, remove_all_from_list
end interface
interface erase
   module procedure erase_at_position, erase_all_from_list
end interface
interface nitems
   module procedure n_items_in_list
end interface


! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: land_tile.F90,v 20.0 2013/12/13 23:29:28 fms Exp $', &
     tagname = '$Name: tikal $'

! ==== data types ============================================================
! land_tile_type describes the structure of the land model tile; basically
! it is a container for tile-specific data, plus some information common to 
! all of them: fraction of tile area, etc.
type :: land_tile_type
   integer :: tag = 0   ! defines type of the tile 

   real    :: frac      ! fractional tile area, dimensionless
   type(glac_tile_type), pointer :: glac => NULL() ! glacier model data
   type(lake_tile_type), pointer :: lake => NULL() ! lake model data
   type(soil_tile_type), pointer :: soil => NULL() ! soil model data
   type(snow_tile_type), pointer :: snow => NULL() ! snow data
   type(cana_tile_type), pointer :: cana => NULL() ! canopy air data
   type(vegn_tile_type), pointer :: vegn => NULL() ! vegetation model data

   type(diag_buff_type), pointer :: diag => NULL() ! diagnostic data storage
   
   ! data that are carried over from the previous time step
   real :: Sg_dir(NBANDS), Sg_dif(NBANDS) ! fractions of downward direct and 
       ! diffuse short-wave radiation absorbed by ground and snow
   real :: Sv_dir(NBANDS), Sv_dif(NBANDS) ! fractions of downward direct and 
       ! diffuse radiation absorbed by the vegetation.
   real :: land_refl_dir(NBANDS), land_refl_dif(NBANDS)
   
   real :: land_d, land_z0m, land_z0s
   real :: surf_refl_lw ! long-wave reflectivity of the ground surface (possibly snow-covered)
   real :: vegn_refl_lw ! black background long-wave reflectivity of the vegetation canopy
   real :: vegn_tran_lw ! black background long-wave transmissivity of the vegetation canopy

   real :: lwup     = 200.0  ! upward long-wave flux from the entire land (W/m2), the result of
           ! the implicit time step -- used in update_bc_fast to return to the flux exchange.
   real :: e_res_1  = 0.0 ! energy residual in canopy air EB equation
   real :: e_res_2  = 0.0 ! energy residual in canopy EB equation
   real :: runon_l  = 0.0 ! water discharged by rivers into the tile, kg/(m2 s)
   real :: runon_s  = 0.0 ! snow discharged by rivers into the tile, kg/(m2 s)
   real :: runon_H  = 0.0 ! heat carried by water discharged by rivers into the tile, W/m2
   real :: runon_Hl  = 0.0 ! heat carried by water discharged by rivers into the tile, W/m2
   real :: runon_Hs  = 0.0 ! heat carried by water discharged by rivers into the tile, W/m2
end type land_tile_type

! tile_list_type provides a container for the tiles
type :: land_tile_list_type
   private
   type(land_tile_list_node_type), pointer :: head => NULL()
end type land_tile_list_type

! land_tile_enum_type provides a enumerator of tiles -- a data structure
! that allows to walk through all tiles in a container (or a 2D array of 
! containers) without bothering with details of container implementation
type :: land_tile_enum_type
   private
   type(land_tile_list_type), pointer :: &
        tiles(:,:) => NULL()  ! pointer to array of tiles to walk -- may be disassociated
   integer :: i=0,j=0 ! indices in the above array
   integer :: k=0     ! number of the current tile in its container
   integer :: io,jo   ! offsets of indices (to keep track of non-1 ubounds of tiles array)
   type(land_tile_list_node_type), pointer :: node => NULL() ! pointer to the current container node
end type land_tile_enum_type

! private type -- used internally to implement tile lists
type :: land_tile_list_node_type
   type(land_tile_list_node_type), pointer :: prev => NULL()
   type(land_tile_list_node_type), pointer :: next => NULL()
   type(land_tile_type), pointer :: data => NULL()
end type land_tile_list_node_type

! ==== module data ===========================================================
integer :: n_created_land_tiles = 0 ! total number of created tiles
integer :: n_deleted_land_tiles = 0 ! total number of deleted tiles


contains 

! #### land_tile_type and operations #########################################

! ============================================================================
! tile constructor: given a list of sub-model tile tags, creates a land tile
! calls sub-tile constructors from individual component models
function land_tile_ctor(frac,glac,lake,soil,vegn,tag) result(tile)
  real   , optional, intent(in) :: frac ! fractional area of tile
  integer, optional, intent(in) :: &
               glac,lake,soil,vegn ! kinds of respective tiles
  integer, optional, intent(in) :: tag  ! general tile tag
  type(land_tile_type), pointer :: tile ! return value

  ! ---- local vars
  integer :: glac_, lake_, soil_, vegn_
  
  ! initialize internal variables
  glac_ = -1 ; if(present(glac)) glac_ = glac
  lake_ = -1 ; if(present(lake)) lake_ = lake
  soil_ = -1 ; if(present(soil)) soil_ = soil
  vegn_ = -1 ; if(present(vegn)) vegn_ = vegn
  
  allocate(tile)
  ! fill common fields
  tile%frac = 0.0 ; if(present(frac)) tile%frac = frac
  tile%tag  = 0   ; if(present(tag))  tile%tag  = tag

  ! create sub-model tiles
  tile%cana => new_cana_tile()
  if(glac_>=0) tile%glac => new_glac_tile(glac_)
  if(lake_>=0) tile%lake => new_lake_tile(lake_)
  tile%snow => new_snow_tile()
  if(soil_>=0) tile%soil => new_soil_tile(soil_)
  if(vegn_>=0) tile%vegn => new_vegn_tile(vegn_)

  ! create a buffer for diagnostic output
  tile%diag=>new_diag_buff()

  ! increment total number of created files for tile statistics
  n_created_land_tiles = n_created_land_tiles + 1

end function land_tile_ctor


! ============================================================================
function land_tile_copy_ctor(t) result(tile)
  type(land_tile_type), intent(in) :: t    ! tile to copy
  type(land_tile_type), pointer :: tile ! return value

  allocate(tile)
  tile = t ! copy all non-pointer members
  if (associated(t%glac)) tile%glac=>new_glac_tile(t%glac)
  if (associated(t%lake)) tile%lake=>new_lake_tile(t%lake)
  if (associated(t%soil)) tile%soil=>new_soil_tile(t%soil)
  if (associated(t%snow)) tile%snow=>new_snow_tile(t%snow)
  if (associated(t%cana)) tile%cana=>new_cana_tile(t%cana)
  if (associated(t%vegn)) tile%vegn=>new_vegn_tile(t%vegn)

  if (associated(t%diag)) tile%diag=>new_diag_buff(t%diag)
end function land_tile_copy_ctor


! ============================================================================
! tile destructor -- releases memory occupied by the tile;
! calls sub-model tile destructors to free the memory of the components
subroutine delete_land_tile(tile)
  type(land_tile_type), pointer :: tile ! tile to delete

  if (.not.associated(tile)) return

  if (associated(tile%glac)) call delete_glac_tile(tile%glac)
  if (associated(tile%lake)) call delete_lake_tile(tile%lake)
  if (associated(tile%soil)) call delete_soil_tile(tile%soil)
  if (associated(tile%snow)) call delete_snow_tile(tile%snow)
  if (associated(tile%cana)) call delete_cana_tile(tile%cana)
  if (associated(tile%vegn)) call delete_vegn_tile(tile%vegn)
  
  ! deallocate diagnostic storage
  call delete_diag_buff(tile%diag)

  ! release the tile memory
  deallocate(tile)
  
  ! increment the number of deleted files for tile statistics
  n_deleted_land_tiles = n_deleted_land_tiles + 1

end subroutine delete_land_tile


! ============================================================================
! returns tags of the component model tiles
subroutine get_tile_tags(tile,land,glac,lake,soil,snow,cana,vegn)
   type(land_tile_type), intent(in)  :: tile
   integer, optional,    intent(out) :: land,glac,lake,soil,snow,cana,vegn

   if(present(land)) land=tile%tag
   if(present(glac)) then
      glac=-HUGE(glac)
      if (associated(tile%glac)) glac=get_glac_tile_tag(tile%glac)
   endif
   if(present(lake)) then
      lake=-HUGE(lake)
      if (associated(tile%lake)) lake=get_lake_tile_tag(tile%lake)
   endif
   if(present(soil)) then
      soil=-HUGE(soil)
      if (associated(tile%soil)) soil=get_soil_tile_tag(tile%soil)
   endif
   if(present(snow)) then
      snow=-HUGE(snow)
      if (associated(tile%snow)) snow=get_snow_tile_tag(tile%snow)
   endif
   if(present(cana)) then
      cana=-HUGE(cana)
      if (associated(tile%cana)) cana=get_cana_tile_tag(tile%cana)
   endif
   if(present(vegn)) then
      vegn=-HUGE(vegn)
      if (associated(tile%vegn)) vegn=get_vegn_tile_tag(tile%vegn)
   endif
end subroutine


! ============================================================================
! returns totals water and ice masses associated with tile
subroutine get_tile_water(tile, lmass, fmass)
  type(land_tile_type), intent(in) :: tile
  real, intent(out) :: lmass, fmass ! liquid and solid water masses, kg/m2

  ! ---- local vars
  real :: lm, fm

  lmass = 0; fmass = 0
  if (associated(tile%cana)) then
     call cana_tile_stock_pe(tile%cana, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%glac)) then
     call glac_tile_stock_pe(tile%glac, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%lake)) then
     call lake_tile_stock_pe(tile%lake, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%soil)) then
     call soil_tile_stock_pe(tile%soil, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%snow)) then
     call snow_tile_stock_pe(tile%snow, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif
  if (associated(tile%vegn)) then
     call vegn_tile_stock_pe(tile%vegn, lm, fm)
     lmass = lmass+lm ; fmass = fmass + fm
  endif

end subroutine


! ============================================================================
! returns total tile carbon, kg C/m2
function land_tile_carbon(tile) result(carbon) ; real carbon
  type(land_tile_type), intent(in) :: tile

  carbon = 0
  if (associated(tile%cana)) &
     carbon = carbon + cana_tile_carbon(tile%cana)
  if (associated(tile%vegn)) &
     carbon = carbon + vegn_tile_carbon(tile%vegn)
  if (associated(tile%soil)) &
     carbon = carbon + soil_tile_carbon(tile%soil)
end function 


! ============================================================================
! returns total heat content of the tile
function land_tile_heat(tile) result(heat) ; real heat
  type(land_tile_type), intent(in) :: tile

  heat = 0
  if (associated(tile%cana)) &
       heat = heat+cana_tile_heat(tile%cana)
  if (associated(tile%glac)) &
       heat = heat+glac_tile_heat(tile%glac)
  if (associated(tile%lake)) &
       heat = heat+lake_tile_heat(tile%lake)
  if (associated(tile%soil)) &
       heat = heat+soil_tile_heat(tile%soil)
  if (associated(tile%snow)) &
       heat = heat+snow_tile_heat(tile%snow)
  if (associated(tile%vegn)) &
       heat = heat+vegn_tile_heat(tile%vegn)
end function


! ============================================================================
! returns true if two land tiles can be merged 
function land_tiles_can_be_merged(tile1,tile2) result (answer)
   logical :: answer ! returned value
   type(land_tile_type), intent(in) :: tile1, tile2
   
   ! make sure that the two tiles have the same components. For 
   ! uniformity every component is checked, even though snow and
   ! cana are always present in current design
   answer = (associated(tile1%glac).eqv.associated(tile2%glac)).and. &
            (associated(tile1%lake).eqv.associated(tile2%lake)).and. &
            (associated(tile1%soil).eqv.associated(tile2%soil)).and. &
            (associated(tile1%snow).eqv.associated(tile2%snow)).and. &
            (associated(tile1%cana).eqv.associated(tile2%cana)).and. &
            (associated(tile1%vegn).eqv.associated(tile2%vegn))
     
   if (answer.and.associated(tile1%glac)) &
      answer = answer.and.glac_tiles_can_be_merged(tile1%glac,tile2%glac)
   if (answer.and.associated(tile1%lake)) &
      answer = answer.and.lake_tiles_can_be_merged(tile1%lake,tile2%lake)
   if (answer.and.associated(tile1%soil)) &
      answer = answer.and.soil_tiles_can_be_merged(tile1%soil,tile2%soil)
   if (answer.and.associated(tile1%cana)) &
      answer = answer.and.cana_tiles_can_be_merged(tile1%cana,tile2%cana)
   if (answer.and.associated(tile1%snow)) &
      answer = answer.and.snow_tiles_can_be_merged(tile1%snow,tile2%snow)
   if (answer.and.associated(tile1%vegn)) &
      answer = answer.and.vegn_tiles_can_be_merged(tile1%vegn,tile2%vegn)
   
end function

! ============================================================================
! merges the two tiles, putting merged state into the second tile. The first
! tile is unchanged
subroutine merge_land_tiles(tile1,tile2)
  type(land_tile_type), intent(in)    :: tile1
  type(land_tile_type), intent(inout) :: tile2

  ! ---- local vars
  real :: x1,x2

  if(associated(tile1%glac)) &
       call merge_glac_tiles(tile1%glac, tile1%frac, tile2%glac, tile2%frac)
  if(associated(tile1%lake)) &
       call merge_lake_tiles(tile1%lake, tile1%frac, tile2%lake, tile2%frac)
  if(associated(tile1%soil)) &
       call merge_soil_tiles(tile1%soil, tile1%frac, tile2%soil, tile2%frac)
  
  if(associated(tile1%cana)) &
       call merge_cana_tiles(tile1%cana, tile1%frac, tile2%cana, tile2%frac)
  if(associated(tile1%snow)) &
       call merge_snow_tiles(tile1%snow, tile1%frac, tile2%snow, tile2%frac)

  if(associated(tile1%vegn)) &
       call merge_vegn_tiles(tile1%vegn, tile1%frac, tile2%vegn, tile2%frac)

  ! calculate normalized weights
  x1 = tile1%frac/(tile1%frac+tile2%frac)
  x2 = 1.0 - x1

#define __MERGE__(field) tile2%field = x1*tile1%field + x2*tile2%field
  __MERGE__(lwup)
  __MERGE__(e_res_1)
  __MERGE__(e_res_2)
  __MERGE__(runon_l)
  __MERGE__(runon_s)
  __MERGE__(runon_H)
  __MERGE__(runon_Hl)
  __MERGE__(runon_Hs)
#undef __MERGE__

  tile2%frac = tile1%frac + tile2%frac
end subroutine

! #### tile container ########################################################

! ============================================================================
! tile list constructor: initializes essential innards of tile collection
! for future use. In current implementation, it is safe to call this function
! on a tile list more then once
subroutine land_tile_list_init(list)
  type(land_tile_list_type), intent(inout) :: list

  if (.not.associated(list%head)) then
     allocate(list%head)
     list%head%prev=>list%head
     list%head%next=>list%head
  endif
end subroutine land_tile_list_init

! ============================================================================
! tile list destructor: destroys the list of tiles. NOTE that it also destroys
! all the tiles that are still in the list.
subroutine land_tile_list_end(list)
  type(land_tile_list_type), intent(inout) :: list

  if(associated(list%head)) then
     call erase(list)
     deallocate(list%head)
  endif
end subroutine land_tile_list_end

! ============================================================================
subroutine check_tile_list_inited(list)
  type(land_tile_list_type), intent(in) :: list

  if (.not.associated(list%head)) &
     call error_mesg('land_tile_mod','tile container was not initialized before use', FATAL)
     
end subroutine


! ============================================================================
! returns true is the list is empty
function empty(list)
  logical empty
  type(land_tile_list_type), intent(in) :: list

  empty = .not.associated(list%head)
  if (.not.empty) &
       empty = associated(list%head%next,list%head)

end function empty

! ============================================================================
! returns the number of items currently stored in the list
function n_items_in_list(list) result (n)
  type(land_tile_list_type), intent(in) :: list
  integer :: n

  type(land_tile_list_node_type), pointer :: node

  n=0; 
  if(.not.associated(list%head)) return

  node => list%head%next
  do while ( .not.(associated(node,list%head)) )
     n = n+1
     node => node%next
  enddo
end function n_items_in_list

! ============================================================================
subroutine insert_in_list(tile,list)
  type(land_tile_type),           pointer :: tile
  type(land_tile_list_type), intent(inout) :: list

  call insert_at_position(tile,tail_elmt(list))

end subroutine insert_in_list


! ============================================================================
subroutine remove_all_from_list(list)
  type(land_tile_list_type), intent(inout) :: list
  
  type(land_tile_enum_type) :: ce
  ce=first_elmt(list)
  do while(ce/=tail_elmt(list))
     call remove_at_position(ce)
  enddo
end subroutine remove_all_from_list


! ============================================================================
subroutine erase_all_from_list(list)
  type(land_tile_list_type), intent(inout) :: list
  
  type(land_tile_enum_type) :: ce
  ce=first_elmt(list)
  do while(ce/=tail_elmt(list))
     call erase_at_position(ce)
  enddo
end subroutine erase_all_from_list



! #### tile container enumerator #############################################

! ============================================================================
! returns enumerator pointing to the first element of the container
function land_tile_list_begin_0d(list) result(ce)
  type(land_tile_enum_type) :: ce  ! return value
  type(land_tile_list_type), intent(in) :: list

  call check_tile_list_inited(list)
  ce%node=>list%head%next
  ce%i = 1 ; ce%j = 1 ; ce%k = 1 
end function


! ============================================================================
! returns enumerator pointing to the first element of the 2D array of 
! containers 
function land_tile_list_begin_2d(tiles, is, js) result(ce)
  type(land_tile_enum_type) :: ce  ! return value
  type(land_tile_list_type), intent(in), target :: tiles(:,:)
  integer, intent(in), optional :: is,js ! origin of the array
  
  integer :: i,j

  ! list up pointer to the array of containers
  ce%tiles=>tiles

  ! initialize offsets of indices
  ce%io = 0; ce%jo = 0;
  if(present(is)) ce%io = is-lbound(tiles,1)
  if(present(js)) ce%jo = js-lbound(tiles,2)

  ! initialize current position in the array of containers -- find
  ! first non-empty container and list the pointer to the current
  ! container node
  ce%k = 1
  do j = lbound(tiles,2),ubound(tiles,2)
  do i = lbound(tiles,1),ubound(tiles,1)
     call check_tile_list_inited(tiles(i,j))
     ce%node => tiles(i,j)%head%next
     ce%i = i ; ce%j = j
     if(associated(ce%node%data)) return
  enddo
  enddo
end function


! ============================================================================
! returns enumerator pointing to the end of container: actually the next element 
! behind the last element of the container
function land_tile_list_end_0d(list) result (ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_list_type), intent(in) :: list
  
  call check_tile_list_inited(list)
  ce%node=>list%head
  ce%i = 1 ; ce%j = 1 ; ce%k = nitems(list)+1
end function


! ============================================================================
! returns enumerator pointing to the end of 2D array of containers: actually 
! the next element behind the last element of the last container
function land_tile_list_end_2d(tiles,is,js) result (ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_list_type), intent(in), target :: tiles(:,:)
  integer, intent(in), optional :: is,js ! lower boundaries of the array

  ! list up pointer to the array of containers
  ce%tiles=>tiles

  ! initialize offsets of indices
  ce%io = 0; ce%jo = 0;
  if(present(is)) ce%io = is-lbound(tiles,1)
  if(present(js)) ce%jo = js-lbound(tiles,2)

  ! initialize current position in the array of containers 
  ce%i = ubound(tiles,1)
  ce%j = ubound(tiles,2)
  ce%k = nitems(tiles(ce%i,ce%j))+1

  ! list the pointer to the current tile
  call check_tile_list_inited(tiles(ce%i,ce%j))
  ce%node=>tiles(ce%i,ce%j)%head

end function


! ============================================================================
! returns enumerator pointing to the next element of the container.
function next_elmt(pos0) result(ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_enum_type), intent(in) :: pos0

  integer :: is,ie,js,je

  ce = pos0
  ce%node => ce%node%next ; ce%k = ce%k+1
  if(associated(ce%tiles)) then
     is = lbound(ce%tiles,1); ie = ubound(ce%tiles,1)
     js = lbound(ce%tiles,2); je = ubound(ce%tiles,2)
     do while(.not.associated(ce%node%data))
        ce%k = 1; ! reset tile index
        if(ce%i<ie)then
           ce%i = ce%i+1
        else if(ce%j<je) then
           ce%i = is
           ce%j = ce%j + 1
        else
           return
        endif
        call check_tile_list_inited(ce%tiles(ce%i,ce%j))
        ce%node => ce%tiles(ce%i,ce%j)%head%next
     enddo
  endif
end function


! ============================================================================
! returns enumerator pointing to the previous element of the container.
function prev_elmt(pos0) result(ce)
  type(land_tile_enum_type) :: ce ! return value
  type(land_tile_enum_type), intent(in) :: pos0

  integer :: is,ie,js,je

  ce = pos0
  ce%node => ce%node%prev ; ce%k = ce%k - 1
  if(associated(ce%tiles)) then
     is = lbound(ce%tiles,1); ie = ubound(ce%tiles,1)
     js = lbound(ce%tiles,2); je = ubound(ce%tiles,2)
     do while(.not.associated(ce%node%data))
        ce%k = 1; ! reset tile index
        if(ce%i>is)then
           ce%i = ce%i - 1
        else if(ce%j>js) then
           ce%i = ie
           ce%j = ce%j - 1
        else
           return
        endif
        call check_tile_list_inited(ce%tiles(ce%i,ce%j))
        ce%node => ce%tiles(ce%i,ce%j)%head%prev
        ce%k    =  nitems(ce%tiles(ce%i,ce%j))
     enddo
  endif

end function

! ============================================================================
! returns TRUE if both enums refer to the same list node (and, therefore, tile)
! or if both do not refer to anything. 
function enums_are_equal(pos1,pos2) result(ret)
  logical :: ret ! return value
  type(land_tile_enum_type), intent(in) :: pos1,pos2
  
  if(associated(pos1%node)) then
     ret = associated(pos1%node,pos2%node)
  else
     ret = .not.associated(pos2%node)
  endif
end function enums_are_equal

! ============================================================================
! returns TRUE if two enumerators are not equal
function enums_are_not_equal(pos1,pos2) result(ret)
  logical :: ret ! return value
  type(land_tile_enum_type), intent(in) :: pos1,pos2

  ret=.not.enums_are_equal(pos1,pos2)
end function enums_are_not_equal

! ============================================================================
! returns pointer to the tile currently addressed by the enumerator
function current_tile(ce) result(ptr)
  type(land_tile_type), pointer :: ptr ! return value
  type(land_tile_enum_type), intent(in) :: ce 

  ptr => ce%node%data
end function 

! ============================================================================
! returns indices corresponding to the enumerator; for enumerator associated
! with a single tile list (not with 2D array of lists) returned i and j are 
! equal to 1
subroutine get_elmt_indices(ce,i,j,k)
  type(land_tile_enum_type), intent(in) :: ce 
  integer, intent(out), optional :: i,j,k

  if (present(i)) i = ce%i+ce%io
  if (present(j)) j = ce%j+ce%jo
  if (present(k)) k = ce%k

end subroutine

! ============================================================================
! inserts tile at the position indicated by enumerator: in fact right in front 
! of it. 
subroutine insert_at_position(tile,ce)
  type(land_tile_type),         pointer :: tile
  type(land_tile_enum_type), intent(in) :: ce

  ! local vars
  type(land_tile_list_node_type), pointer :: node,n,p

  allocate(node)
  node%data=>tile

  n=>ce%node  ; p=>n%prev

  node%next=>n ; node%prev=>p
  n%prev=>node ; p%next=>node
  
end subroutine insert_at_position

! ============================================================================
subroutine remove_at_position(enum)
  type(land_tile_enum_type), intent(inout) :: enum

  type(land_tile_list_node_type),pointer :: n,p
  type(land_tile_enum_type) :: next
  
  if(.not.associated(enum%node)) &
     call error_mesg('remove_at_position','attempt to remove tail element of a list', FATAL)

  next = next_elmt(enum)
  
  n => enum%node%next
  p => enum%node%prev

  n%prev=>p ; p%next=>n
  deallocate(enum%node)
  
  enum=next
  if(enum%k>1) enum%k = enum%k-1

end subroutine remove_at_position

! ============================================================================
subroutine erase_at_position(ce)
  type(land_tile_enum_type), intent(inout) :: ce

  type(land_tile_type), pointer :: tile

  tile=>current_tile(ce)
  call remove_at_position(ce)
  call delete_land_tile(tile)

end subroutine erase_at_position


! ============================================================================
function tile_is_selected(tile, sel)
! returns true if the tile fits specified selector
  logical :: tile_is_selected
  type(land_tile_type)    , intent(in) :: tile
  type(tile_selector_type), intent(in) :: sel

  tile_is_selected = .FALSE.
  select case(sel%tag)
  case(SEL_SOIL)
     if(associated(tile%soil)) & 
          tile_is_selected = soil_is_selected(tile%soil,sel)
  case(SEL_VEGN)
     if(associated(tile%vegn)) &
          tile_is_selected = vegn_is_selected(tile%vegn,sel)
  case(SEL_LAKE)
     if(associated(tile%lake)) &
          tile_is_selected = lake_is_selected(tile%lake,sel)
  case(SEL_GLAC)
     if(associated(tile%glac)) &
          tile_is_selected = glac_is_selected(tile%glac,sel)
  case(SEL_SNOW)
     if(associated(tile%snow)) &
          tile_is_selected = snow_is_selected(tile%snow,sel)
  case(SEL_CANA)
     if(associated(tile%cana)) &
          tile_is_selected = cana_is_selected(tile%cana,sel)
  case default
     tile_is_selected=.true.
  end select

end function tile_is_selected


! ============================================================================
subroutine print_land_tile_info(tile)
  type(land_tile_type), intent(in) :: tile
  
  write(*,'("(tag =",i3,", frac =",f7.4)',advance='no') tile%tag, tile%frac
  if(associated(tile%lake)) write(*,'(a,i3)',advance='no')', lake =',tile%lake%tag
  if(associated(tile%soil)) write(*,'(a,i3)',advance='no')', soil =',tile%soil%tag
  if(associated(tile%glac)) write(*,'(a,i3)',advance='no')', glac =',tile%glac%tag
  if(associated(tile%snow)) write(*,'(a,i3)',advance='no')', snow =',tile%snow%tag
  if(associated(tile%cana)) write(*,'(a)',advance='no')', cana'
  if(associated(tile%vegn)) write(*,'(a)',advance='no')', vegn'
  write(*,'(")")',advance='no')
  
end subroutine


! ============================================================================
subroutine print_land_tile_statistics()
  write(*,*)'Total number of created land_tiles =',n_created_land_tiles
  write(*,*)'Total number of deleted land_tiles =',n_deleted_land_tiles
end subroutine

end module land_tile_mod
