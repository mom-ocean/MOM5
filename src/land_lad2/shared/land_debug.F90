module land_debug_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: &
     error_mesg, file_exist, check_nml_error, stdlog, &
     write_version_number, close_file, mpp_pe, mpp_npes, mpp_root_pe, FATAL, NOTE
use time_manager_mod, only : &
     get_date
use grid_mod, only: &
     get_grid_ntiles
use land_data_mod, only : &
     lnd

! NOTE TO SELF: the "!$" sentinels are not comments: they are compiled if OpenMP 
! support is turned on
!$ use omp_lib, only: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM

implicit none
private

! ==== public interfaces =====================================================
public :: land_debug_init
public :: land_debug_end

public :: set_current_point
public :: get_current_point
public :: current_i, current_j, current_k, current_face
public :: is_watch_point
public :: is_watch_cell
public :: get_watch_point

public :: check_temp_range
! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'land_debug',&
    version     = '$Id: land_debug.F90,v 19.0 2012/01/06 20:41:29 fms Exp $',&
    tagname     = '$Name: siena_201207 $'

! ==== module variables ======================================================
integer, allocatable :: current_debug_level(:)
integer :: mosaic_tile = 0
integer, allocatable :: curr_i(:), curr_j(:), curr_k(:)

!---- namelist ---------------------------------------------------------------
integer :: watch_point(4)=(/0,0,0,1/) ! coordinates of the point of interest, i,j,tile,mosaic_tile
real    :: temp_lo = 120.0 ! lower limit of "reasonable" temperature range, deg K
real    :: temp_hi = 373.0 ! upper limit of "reasonable" temperature range, deg K
namelist/land_debug_nml/ watch_point, temp_lo, temp_hi


contains

! ============================================================================
subroutine land_debug_init()
  ! ---- local vars
  integer :: unit, ierr, io, ntiles
  integer :: max_threads

  call write_version_number(version, tagname)
  
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=land_debug_nml, iostat=io)
  ierr = check_nml_error(io, 'land_debug_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_debug_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_debug_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=land_debug_nml)
  endif
  ! set number of our mosaic tile 
  call get_grid_ntiles('LND',ntiles)
  mosaic_tile = ntiles*mpp_pe()/mpp_npes() + 1  ! assumption

  ! set number of threads and allocate by-thread arrays
    max_threads = 1  
!$  max_threads = OMP_GET_MAX_THREADS()
  allocate(curr_i(max_threads),curr_j(max_threads),curr_k(max_threads))
  allocate(current_debug_level(max_threads))
  current_debug_level(:) = 0
end subroutine land_debug_init

! ============================================================================
subroutine land_debug_end()
  deallocate(curr_i,curr_j,curr_k)
  deallocate(current_debug_level)
end subroutine

! ============================================================================
subroutine set_current_point(i,j,k)
  integer, intent(in) :: i,j,k

  integer :: thread
    thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1

  curr_i(thread) = i ; curr_j(thread) = j ; curr_k(thread) = k

  current_debug_level(thread) = 0
  if ( watch_point(1)==i.and. &
       watch_point(2)==j.and. &
       watch_point(3)==k.and. &
       watch_point(4)==mosaic_tile) then
     current_debug_level(thread) = 1
  endif
end subroutine set_current_point

! ============================================================================
subroutine get_current_point(i,j,k,face)
  integer, intent(out), optional :: i,j,k,face

  integer :: thread
    thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1

  if (present(i)) i = curr_i(thread)
  if (present(j)) j = curr_j(thread)
  if (present(k)) k = curr_k(thread)
  if (present(face)) face = mosaic_tile
end subroutine get_current_point

! ============================================================================
integer function current_i()
  integer :: thread
    thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  current_i = curr_i(thread)
end function

integer function current_j()
  integer :: thread
    thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  current_j = curr_j(thread)
end function

integer function current_k()
  integer :: thread
    thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  current_k = curr_k(thread)
end function

integer function current_face() ; current_face = mosaic_tile ; end function

! ============================================================================
function is_watch_point()
  logical :: is_watch_point

  integer :: thread
    thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  is_watch_point = (current_debug_level(thread) > 0)
end function is_watch_point

! ============================================================================
! returns true, if the watch point is within the grid cell, regardless of
! the tile number
function is_watch_cell()
  logical :: is_watch_cell
  is_watch_cell = ( current_i() == watch_point(1) &
              .and. current_j() == watch_point(2) &
              .and. mosaic_tile == watch_point(4) )
end function is_watch_cell


! ============================================================================
subroutine get_watch_point(i,j,k,face)
  integer, intent(out), optional :: i,j,k,face
  if (present(i)) i = watch_point(1)
  if (present(j)) j = watch_point(2)
  if (present(k)) k = watch_point(3)
  if (present(face)) face = watch_point(4)
end subroutine get_watch_point

! ============================================================================
! checks if the temperature within reasonable range, and prints a message
! if it isn't
subroutine check_temp_range(temp, tag, varname)
  real, intent(in) :: temp ! temperature to check
  character(*), intent(in) :: tag ! tag to print
  character(*), intent(in) :: varname ! name of the variable for printout

  ! ---- local vars
  integer :: y,mo,d,h,m,s ! components of date
  integer :: thread


  if(temp_lo<temp.and.temp<temp_hi) then
     return
  else
     thread = 1
!$   thread = OMP_GET_THREAD_NUM()+1
     call get_date(lnd%time,y,mo,d,h,m,s)
     write(*,'(a," : ",a,g,4(x,a,i4),x,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
          trim(tag), trim(varname)//' out of range: value=', &
         temp,'at i=',curr_i(thread),'j=',curr_j(thread),'tile=',curr_k(thread),'face=',mosaic_tile, &
         'time=',y,mo,d,h,m,s
  endif
end subroutine 


end module land_debug_mod
