module land_tile_diag_mod

use mpp_mod,            only : mpp_sum
use time_manager_mod,   only : time_type
use diag_axis_mod,      only : get_axis_length
use diag_manager_mod,   only : register_diag_field, register_static_field, &
     send_data
use diag_util_mod,      only : log_diag_field_info
use fms_mod,            only : write_version_number, error_mesg, string, FATAL

use land_tile_selectors_mod, only : tile_selectors_init, tile_selectors_end, &
     tile_selector_type, register_tile_selector, selector_suffix, &
     get_n_selectors, get_selector
use land_tile_mod,      only : land_tile_type, diag_buff_type, &
     land_tile_list_type, first_elmt, tail_elmt, next_elmt, get_elmt_indices, &
     land_tile_enum_type, operator(/=), current_tile, &
     tile_is_selected
use land_data_mod,      only : lnd
use tile_diag_buff_mod, only : diag_buff_type, realloc_diag_buff

implicit none
private


! ==== public interface ======================================================
public :: tile_diag_init
public :: tile_diag_end

public :: diag_buff_type

public :: register_tiled_diag_field
public :: register_tiled_static_field
public :: add_tiled_diag_field_alias
public :: add_tiled_static_field_alias
public :: send_tile_data
public :: send_tile_data_r0d_fptr, send_tile_data_r1d_fptr
public :: send_tile_data_i0d_fptr

public :: dump_tile_diag_fields

! codes of tile aggregaton operations
public :: OP_AVERAGE, OP_SUM

interface send_tile_data
   module procedure send_tile_data_0d
   module procedure send_tile_data_1d
end interface
! ==== end of public interface ===============================================


! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'lan_tile_diag_mod', &
     version     = '$Id: land_tile_diag.F90,v 20.0 2013/12/13 23:29:55 fms Exp $', &
     tagname     = '$Name: tikal $'

integer, parameter :: INIT_FIELDS_SIZE     = 1     ! initial size of the fields array
integer, parameter :: BASE_TILED_FIELD_ID  = 65536 ! base value for tiled field 
! ids, to distinguish them from regular diagnostic fields. All IDs of tiled
! (that is, registered by register_*tiled_field functions are larger than 
! BASE_TILED_FIELD_ID)
integer, parameter :: MIN_DIAG_BUFFER_SIZE = 1     ! min size of the per-tile diagnostic buffer
! operations used for tile data aggregation
integer, parameter :: OP_AVERAGE = 0 ! weighted average of tile values
integer, parameter :: OP_SUM     = 1 ! sum of all tile values


! ==== derived types =========================================================
type :: tiled_diag_field_type
   integer, pointer :: ids(:) => NULL()
   integer :: offset ! offset of the field data in the buffer
   integer :: size   ! size of the field data in the per-tile buffers
   integer :: op     ! aggregation operation
   logical :: static ! if true, the diag field is static
   integer :: n_sends! number of data points sent to the field since last dump
   integer :: alias = 0 ! ID of the first alias in the chain
   character(32) :: module,name ! for debugging purposes only
end type tiled_diag_field_type


! ==== module data ===========================================================
logical :: module_is_initialized = .false.

! list of registered fields
type(tiled_diag_field_type), pointer :: fields(:) => NULL()
integer :: n_fields       = 0 ! current number of diag fields
integer :: current_offset = 1 ! current total size of the diag fields per tile



contains



! ============================================================================
subroutine tile_diag_init()

  if (module_is_initialized) return

  module_is_initialized = .true.
  call write_version_number(version, tagname)

  ! initialize diag selectors
  call tile_selectors_init()
  call register_tile_selector('')

  ! initailze global data
  allocate(fields(INIT_FIELDS_SIZE))
  n_fields       = 0
  current_offset = 1

end subroutine tile_diag_init



! ============================================================================
subroutine tile_diag_end()

  integer :: i

  ! deallocate global data
  do i = 1, n_fields
     deallocate(fields(i)%ids)
  end do
  deallocate(fields)

  ! destroy selectors
  call tile_selectors_end()

  module_is_initialized = .false.

end subroutine tile_diag_end


! ============================================================================
function register_tiled_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code
  
  id = reg_field(.false., module_name, field_name, init_time, axes, long_name, &
         units, missing_value, range, op=op)

end function

! ============================================================================
function register_tiled_static_field(module_name, field_name, axes, &
     long_name, units, missing_value, range, require, op) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  integer,          intent(in), optional :: op ! aggregation operation code
  
  ! --- local vars
  type(time_type) :: init_time

  id = reg_field(.true., module_name, field_name, init_time, axes, long_name, &
         units, missing_value, range, require, op)

end function


! ============================================================================
subroutine add_tiled_static_field_alias(id0, module_name, field_name, axes, &
     long_name, units, missing_value, range, op)
  integer,          intent(inout) :: id0 ! id of the original diag field on input;
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code

  ! --- local vars
  type(time_type) :: init_time

  call reg_field_alias(id0, .true., module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)
end subroutine


! ============================================================================
subroutine add_tiled_diag_field_alias(id0, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)
  integer,          intent(inout) :: id0 ! id of the original diag field on input;
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code

  call reg_field_alias(id0, .false., module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)
end subroutine

! ============================================================================
subroutine reg_field_alias(id0, static, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)


  integer,          intent(inout) :: id0 ! id of the original diag field on input;
  logical,          intent(in) :: static
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code
  
  ! local vars
  integer :: id1
  integer :: ifld0, ifld1
  
  if (id0>0) then
    ifld0 = id0-BASE_TILED_FIELD_ID
    if (ifld0<1.or.ifld0>n_fields) &
       call error_mesg(module_name, 'incorrect index ifld0 '//string(ifld0)//&
                    ' in definition of tiled diag field alias "'//&
                    trim(module_name)//'/'//trim(field_name)//'"', FATAL)
    id1 = reg_field(static, module_name, field_name, init_time, axes, long_name, &
          units, missing_value, range, op=op, offset=fields(ifld0)%offset)
    if (id1>0) then
      ifld1 = id1-BASE_TILED_FIELD_ID
      ! check that sizes of the fields are identical
      if (fields(ifld0)%size/=fields(ifld1)%size) &
         call error_mesg(module_name, 'sizes of diag field "'//           &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//    &
           '" and its alias "'//trim(module_name)//'/'//trim(field_name)//&
           '" are not the same', FATAL)
      ! check that "static" status of the fields is the same
      if(fields(ifld0)%static.and..not.fields(ifld1)%static) &
         call error_mesg(module_name,                                     &
           'attempt to register non-static alias"'//                      &
           trim(module_name)//'/'//trim(field_name)//                     &
           '" of static field "'//                                        &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//'"',&
           FATAL)
      if(.not.fields(ifld0)%static.and.fields(ifld1)%static) &
         call error_mesg(module_name,                                     &
           'attempt to register static alias"'//                          &
           trim(module_name)//'/'//trim(field_name)//                     &
           '" of non-static field "'//                                    &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//'"',&
           FATAL)

      ! copy alias field from the original into the alias, to preserve the chain
      fields(ifld1)%alias = fields(ifld0)%alias
      ! update alias field in the head of alias chain
      fields(ifld0)%alias = ifld1
    endif
  else
    ! the "main" field has not been registered, so simply redister the alias
    ! as a diag field
    id0 = reg_field(static, module_name, field_name, init_time, axes, long_name, &
          units, missing_value, range, op=op)
  endif
end subroutine

! ============================================================================
! provides unified interface for registering a diagnostic field with full set
! of selectors
function reg_field(static, module_name, field_name, init_time, axes, &
     long_name, units, missing_value, range, require, op, offset) result(id)
 
  integer :: id

  logical,          intent(in) :: static
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  integer,          intent(in), optional :: op
  integer,          intent(in), optional :: offset

  ! ---- local vars
  integer, pointer :: diag_ids(:) ! ids returned by FMS diag manager for each selector
  integer :: i
  integer :: isel    ! selector iterator
  integer :: n_selectors ! number of registered diagnostic tile selectors
  type(tiled_diag_field_type), pointer :: new_fields(:)
  type(tile_selector_type) :: sel
  ! ---- global vars: n_fields, fields, current_offset -- all used and updated

  ! log diagnostic field information
  call log_diag_field_info ( module_name, trim(field_name), axes, long_name, units,&
                             missing_value, range, dynamic=.not.static )
  ! go through all possible selectors and try to register a diagnostic field 
  ! with the name derived from field name and selector; if any of the 
  ! registrations succeeds, return a tiled field id, otherwise return 0.
  ! Note that by design one of the selectors have empty name and selects all
  ! the tiles.
  id = 0
  n_selectors = get_n_selectors()
  allocate(diag_ids(n_selectors))
  
  do isel = 1, n_selectors
     ! try to register field+selector pair with FMS diagnostics manager
     sel = get_selector(isel)
     diag_ids(isel) = reg_field_set(static, sel, module_name, field_name, axes, &
          init_time, long_name, units, missing_value, range, require)

  enddo
  
  if(any(diag_ids>0)) then
     ! if any of the field+selector pairs was found for this field, an entry
     ! must be added to the table of tile diagnostic fields

     ! if there is not enough slots in the field table to add another one,
     ! allocate more space
     if(n_fields>=size(fields)) then
        allocate(new_fields(max(2*n_fields,1)))
        new_fields(1:n_fields) = fields(1:n_fields)
        deallocate(fields)
        fields => new_fields
     endif
     ! add the current field to the field table
     n_fields = n_fields+1
     id       = n_fields
     ! set the array of FMS diagnostic field IDs for each selector
     fields(id)%ids => diag_ids
     ! set the field offset in the diagnostic buffers
     if (present(offset)) then
        fields(id)%offset = offset
     else
        fields(id)%offset = current_offset
     endif  
     ! calculate field size per tile and increment current offset to
     ! reserve space in per-tile buffers. We assume that the first two axes 
     ! are horizontal coordinates, so their size is not taken into account
     fields(id)%size = 1
     do i = 3, size(axes(:))
        fields(id)%size = fields(id)%size * get_axis_length(axes(i))
     enddo
     ! if offset is present in the list of the arguments, it means that we don't
     ! want to increase the current_offset -- this is an alias field
     if (.not.present(offset)) &
        current_offset = current_offset + fields(id)%size
     ! store the code of the requested tile aggregation operation
     if(present(op)) then
        fields(id)%op = op
     else
        fields(id)%op = OP_AVERAGE
     endif
     ! store the static field flag
     fields(id)%static = static
     ! zero out the number of data points ent to the field
     fields(id)%n_sends = 0
     ! store the name of the field -- for now, only to be able to see what it is 
     ! in the debugger
     fields(id)%module=module_name 
     fields(id)%name=field_name
     ! increment the field id by some (large) number to distinguish it from the 
     ! IDs of regular FMS diagnostic fields
     id = id + BASE_TILED_FIELD_ID
  else
     deallocate(diag_ids)
  endif

end function


! ============================================================================
! provides unified interface for registering a diagnostic field with a given
! selector, whether static or time-dependent
function reg_field_set(static, sel, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, require) result (id)

  integer :: id 

  logical,          intent(in) :: static
  type(tile_selector_type), intent(in) :: sel
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require

  character(len=128) :: fname
  character(len=128) :: lname

  ! form field name as concatenation of name of the field and selector suffix
  fname = trim(field_name)//trim(selector_suffix(sel))
  ! form long name as concatenation of specified long name (if present) and
  ! selector long name
  lname = ''
  if(present(long_name)) lname=long_name
  if(trim(sel%long_name)/='') &
     lname = trim(lname)//' ('//trim(sel%long_name)//')'

  ! try registering diagnostic field with FMS diagnostic manager.
  if (static) then
     id = register_static_field ( module_name, fname,   &
          axes, lname, units, missing_value, range, require, do_not_log=.TRUE. )
  else
     id = register_diag_field ( module_name,  fname,   &
          axes, init_time, lname, units, missing_value, range, &
          mask_variant=.true., do_not_log=.TRUE. )
  endif

end function


! ============================================================================
subroutine send_tile_data_0d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x
  type(diag_buff_type), intent(inout) :: buffer
  
  integer :: idx, i

  if (id <= 0) return

  ! reallocate diagnostic buffer according to the current number and size of 
  ! tiled diag fields
  call realloc_diag_buff(buffer,current_offset)

  ! calculate offset for the current diagnostic field in the buffer
  i = id - BASE_TILED_FIELD_ID 
  idx = fields(i)%offset
  
  ! store the diagnostic data
  buffer%data(idx) = x
  buffer%mask(idx) = .TRUE.

  ! increment sent data counter
  fields(i)%n_sends = fields(i)%n_sends + 1
  ! increment sent data counter in all aliases
  do while(fields(i)%alias>0)
    i=fields(i)%alias
    fields(i)%n_sends = fields(i)%n_sends + 1
  enddo
end subroutine

! ============================================================================
subroutine send_tile_data_1d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x(:)
  type(diag_buff_type), intent(inout) :: buffer

  integer :: is, ie, i
  if (id <= 0) return

  ! reallocate diagnostic buffer according to the current number and size of 
  ! tiled diag fields
  call realloc_diag_buff(buffer, current_offset)
  
  ! calculate offset for the current diagnostic field in the buffer
  i = id - BASE_TILED_FIELD_ID ! index in the array of fields
  is = fields(i)%offset ; ie = is+fields(i)%size-1

  ! store the data
  buffer%data(is:ie) = x(:)
  buffer%mask(is:ie) = .TRUE.

  ! increment sent data counter
  fields(i)%n_sends = fields(i)%n_sends + 1
  ! increment sent data counter in all aliases
  do while(fields(i)%alias>0)
    i=fields(i)%alias
    fields(i)%n_sends = fields(i)%n_sends + 1
  enddo
end subroutine

! NOTE: 2-d fields can be handled similarly to 1-d with reshape

! ============================================================================
subroutine send_tile_data_r0d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       real                , pointer :: ptr  ! returned pointer to the data
     end subroutine fptr 
  end interface

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine


! ============================================================================
subroutine send_tile_data_r1d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       real                , pointer :: ptr(:)  ! returned pointer to the data
     end subroutine fptr 
  end interface

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr(:)     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine


! ============================================================================
subroutine send_tile_data_i0d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  ! subroutine returning the pointer to the tile data
  interface
     subroutine fptr(tile, ptr)
       use land_tile_mod, only : land_tile_type
       type(land_tile_type), pointer :: tile ! input
       integer             , pointer :: ptr  ! returned pointer to the data
     end subroutine fptr 
  end interface

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  integer             , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,real(ptr),tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine


! ============================================================================
subroutine dump_tile_diag_fields(tiles, time)
  type(land_tile_list_type), intent(in) :: tiles(:,:) ! 
  type(time_type)          , intent(in) :: time       ! current time

  ! ---- local vars
  integer :: ifld ! field number
  integer :: isel ! selector number
  type(land_tile_enum_type)     :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: total_n_sends(n_fields)
  ! ---- local static variables -- saved between calls
  logical :: first_dump = .TRUE.

  total_n_sends(:) = fields(1:n_fields)%n_sends
  call mpp_sum(total_n_sends, n_fields, pelist=lnd%pelist)

!$OMP parallel do schedule(dynamic) default(shared) private(ifld,isel)
  do ifld = 1, n_fields
     if (total_n_sends(ifld) == 0) cycle ! no data to send 
     do isel = 1, get_n_selectors()
        if (fields(ifld)%ids(isel) <= 0) cycle
        call dump_diag_field_with_sel ( fields(ifld)%ids(isel), tiles, &
             fields(ifld), get_selector(isel), time )
     enddo
  enddo
  ! zero out the number of data points sent to the field 
  fields(1:n_fields)%n_sends=0

  ! all the data are sent to the output, so set the data presence tag to FALSE 
  ! in all diag buffers in preparation for the next time step
  ce = first_elmt(tiles)
  te = tail_elmt (tiles)
  do while(ce /= te)
    tile => current_tile(ce)       ! get the pointer to the current tile
    tile%diag%mask(:) = .FALSE.
    ce = next_elmt(ce)            ! move to the next position
  enddo
  ! reset the first_dump flag
  first_dump = .FALSE.

end subroutine

! ============================================================================
subroutine dump_diag_field_with_sel(id, tiles, field, sel, time)
  integer :: id
  type(land_tile_list_type),   intent(in) :: tiles(:,:)
  type(tiled_diag_field_type), intent(in) :: field
  type(tile_selector_type)   , intent(in) :: sel
  type(time_type)            , intent(in) :: time ! current time
   
  ! ---- local vars
  integer :: i,j ! iterators
  integer :: is,ie,js,je,ks,ke ! array boundaries
  logical :: used ! value returned from send_data (ignored)
  real, allocatable :: buffer(:,:,:), weight(:,:,:)
  type(land_tile_enum_type)     :: ce, te
  type(land_tile_type), pointer :: tile
  
  ! calculate array boundaries
  is = lbound(tiles,1); ie = ubound(tiles,1)
  js = lbound(tiles,2); je = ubound(tiles,2)
  ks = field%offset   ; ke = field%offset + field%size - 1
  
  ! allocate and initialize temporary buffers
  allocate(buffer(is:ie,js:je,ks:ke), weight(is:ie,js:je,ks:ke))
  buffer(:,:,:) = 0.0
  weight(:,:,:) = 0.0
  
  ! accumulate data
  ce = first_elmt(tiles, is=is, js=js)
  te = tail_elmt (tiles)
  do while(ce /= te)
    tile => current_tile(ce)      ! get the pointer to current tile
    call get_elmt_indices(ce,i,j) ! get the indices of current tile
    ce = next_elmt(ce)           ! move to the next position
    
    if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
    if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
    select case (field%op)
    case (OP_AVERAGE)
       where(tile%diag%mask(ks:ke)) 
          buffer(i,j,:) = buffer(i,j,:) + tile%diag%data(ks:ke)*tile%frac
          weight(i,j,:) = weight(i,j,:) + tile%frac
       end where
    case (OP_SUM)
       where(tile%diag%mask(ks:ke)) 
          buffer(i,j,:) = buffer(i,j,:) + tile%diag%data(ks:ke)
          weight(i,j,:) = 1
       end where
    end select
  enddo

  ! normalize accumulated data
  where (weight>0) buffer=buffer/weight
  
  ! send diag field
  used = send_data ( id, buffer, time, mask=weight>0 )   

  ! clean up temporary data
  deallocate(buffer,weight)

end subroutine


end module land_tile_diag_mod
