module land_utils_mod

use land_tile_mod, only : land_tile_type, land_tile_enum_type, land_tile_list_type, &
     first_elmt, tail_elmt, next_elmt, operator(/=), get_elmt_indices, current_tile

implicit none
private
! ==== public interfaces =====================================================
public :: put_to_tiles_r0d_fptr
public :: put_to_tiles_r1d_fptr
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: land_utils.F90,v 17.0 2009/07/21 03:02:46 fms Exp $', &
     tagname = '$Name: tikal $'

contains

! ============================================================================
subroutine put_to_tiles_r0d_fptr(x2d, tile_map, fptr)
  real, intent(in)                         :: x2d     (:,:)
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       real                , pointer :: ptr  ! returned pointer to the data
     end subroutine fptr 
  end interface

  integer :: i,j
  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr     ! pointer to the data element within a tile

  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     call get_elmt_indices(ce,i,j)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if (associated(ptr)) ptr=x2d(i,j)
     ce=next_elmt(ce)
  enddo
end subroutine


! ============================================================================
subroutine put_to_tiles_r1d_fptr(x2d, tile_map, fptr)
  real, intent(in)                         :: x2d     (:,:,:)
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       real                , pointer :: ptr(:) ! returned pointer to the data
     end subroutine fptr 
  end interface

  integer :: i,j
  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr(:)  ! pointer to the data element within a tile

  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     call get_elmt_indices(ce,i,j)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if (associated(ptr)) ptr(:)=x2d(i,j,:)
     ce=next_elmt(ce)
  enddo
end subroutine

end module land_utils_mod
