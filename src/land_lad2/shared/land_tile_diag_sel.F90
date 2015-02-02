module land_tile_selectors_mod

use fms_mod, only : error_mesg, WARNING

implicit none
private

! ==== public interface ======================================================
public :: tile_selector_type      !

! selector tags
public :: SEL_SOIL, SEL_VEGN, SEL_LAKE, SEL_GLAC, SEL_SNOW, SEL_CANA

public :: tile_selectors_init     ! initialize module
public :: tile_selectors_end      ! clean up ufter ourselves

public :: register_tile_selector  ! register selector for diag field
public :: selector_suffix         ! return suffix for the field name

public :: get_selector            ! array of selectors
public :: get_n_selectors         ! number of available selectors
! ==== end of public interface ===============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_tile_selectors_mod', &
     version     = '$Id: land_tile_diag_sel.F90,v 18.0 2010/03/02 23:37:10 fms Exp $', &
     tagname     = '$Name: tikal $'

integer, parameter :: SEL_LEN           = 16  ! max length of the selector name
integer, parameter :: SEL_LONG_NAME_LEN = 128 ! max name of the selector long name
integer, parameter :: INIT_SEL_SIZE     = 1   ! initial size of the array of selectors

! tags for tile-specific diagnostic selectors
integer, parameter :: SEL_SOIL = 1
integer, parameter :: SEL_VEGN = 2
integer, parameter :: SEL_LAKE = 3
integer, parameter :: SEL_GLAC = 4
integer, parameter :: SEL_SNOW = 5
integer, parameter :: SEL_CANA = 6

! ==== derived types =========================================================
type :: tile_selector_type
   character(len=SEL_LEN)           :: name =''          ! name of the selector
   character(len=SEL_LONG_NAME_LEN) :: long_name = ''    ! long name of the selector
   logical, pointer                 :: mask(:) => NULL() ! mask (selector cache)
   integer :: tag = 0 ! tag of the model component
   integer :: idata1=0, idata2=0 ! integer data
   integer :: rdata1=0, rdata2=0 ! real data
end type tile_selector_type

! ==== module private data ===================================================
logical :: module_is_initialized = .false.

! ==== module public data ====================================================
! array of registered selectors
type(tile_selector_type), pointer :: selectors(:) => NULL()
integer :: n_selectors = 0 ! number of registered selectors

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine tile_selectors_init()

  if (module_is_initialized) return

  module_is_initialized = .true.

  allocate (selectors(INIT_SEL_SIZE))
  n_selectors = 0 ! initialize number of regitered selectors
  ! register couple of default selectors (for all tiles and for each tile)
end subroutine tile_selectors_init


! ============================================================================
subroutine tile_selectors_end()

  integer :: i

  module_is_initialized = .false.
  do i = 1,n_selectors
     if(associated(selectors(i)%mask)) deallocate(selectors(i)%mask)
  enddo
  deallocate(selectors)
  n_selectors = 0
end subroutine tile_selectors_end


! ============================================================================
! registers a selector to be used for diagnostic output 
subroutine register_tile_selector( name, long_name, tag, idata1, idata2, rdata1, rdata2 )
  character(len=*), intent(in) :: name
  character(len=*), intent(in), optional :: long_name
  integer, intent(in), optional :: tag
  integer, intent(in), optional :: idata1, idata2
  real,    intent(in), optional :: rdata1, rdata2

  ! ---- local vars
  type(tile_selector_type), pointer :: new_selectors(:)
  character(len=SEL_LEN) :: name_
  integer :: i

  ! check for conflict of names -- presumably, if the selector was already
  ! registered, then it is an error to register it again
  name_=name
  do i = 1, n_selectors
     if (trim(name_)==trim(selectors(i)%name)) then
        call error_mesg(module_name,'attempt to register selector "'&
             //trim(name)//'" which has already been registered',WARNING)
        return ! just skip it 
     endif
  enddo

  ! allocate additional space for selectors if necessary 
  if(n_selectors >= size(selectors)) then
     allocate(new_selectors(max(n_selectors*2,1)))
     new_selectors(1:n_selectors) = selectors(1:n_selectors)
     deallocate(selectors)
     selectors => new_selectors
  endif

  ! set up the selector values
  n_selectors = n_selectors + 1
  selectors(n_selectors)%name = name
  if (present(long_name)) &
       selectors(n_selectors)%long_name = long_name
  if (present(tag)) selectors(n_selectors)%tag = tag
  if (present(idata1)) selectors(n_selectors)%idata1 = idata1
  if (present(idata2)) selectors(n_selectors)%idata2 = idata2
  if (present(rdata1)) selectors(n_selectors)%rdata1 = rdata1
  if (present(rdata2)) selectors(n_selectors)%rdata2 = rdata2
end subroutine register_tile_selector


! ============================================================================
! returns total number of selectors
function get_n_selectors()
  integer :: get_n_selectors
  get_n_selectors = n_selectors
end function 

! ============================================================================
! retrns n-th selector
function get_selector(n) result (sel)
  type(tile_selector_type) :: sel
  integer, intent(in) :: n

  if (n<1.or.n>get_n_selectors()) return
  sel = selectors(n)

end function 


! ============================================================================
! returns variable suffix for given selector
function selector_suffix(selector)
  character(len=SEL_LEN+1) :: selector_suffix
  type(tile_selector_type), intent(in) :: selector
  
  if(trim(selector%name)=='') then
     selector_suffix = ''
  else
     selector_suffix = '_'//trim(selector%name)
  endif
end function selector_suffix


! ============================================================================
subroutine update_tile_selectors()

  integer :: i
  do i = 1,n_selectors
     if(associated(selectors(i)%mask)) deallocate(selectors(i)%mask)
  enddo
end subroutine update_tile_selectors



end module land_tile_selectors_mod
