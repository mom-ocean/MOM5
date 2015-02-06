#include <fms_platform.h>

module ocean_residency_range_mod  !{
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Ocean age tracer module
!</OVERVIEW>
!
!<DESCRIPTION>
!       This module handles the case where the specified region is determined
!       from a range of another tracer. Any tracer, either prognostic or
!       diagnostic, may be used. The tracer field, mutiplied by an optional
!       scaling factor, will be multiplied by the current value in residency
!       mask array.
!
!       To use this module, in the field table namelist, you must do the following:
!
!       1) set "module_name" to "ocean_residency_range"
!       2) "strings" should have two values, the first is "tracer_range", and the
!          second to the name of a prognostic or diagnostic variable
!       3) set "params" to pairs of points defining the range of the variable
!          to select. For each pair of values, the first value must be less than
!          the second value. The grid cells will satisfy the union of all of the
!          grid cells which satisfy
!               "low value" < "cell value" <= "high value"
!          for at least one pair of points
!       4) no "flags" should be set
!
!       Any prognostic or diagnostic variable may be used. If you wish to allow
!       another variable to be used, you can just create a diagnostic variable
!       for it.
!
!       For an overview of the ocean residency modules, see ocean_residency.F90.
!
!</DESCRIPTION>
!
! $Id: ocean_residency_range.F90,v 20.0 2013/12/14 00:17:12 fms Exp $
!

!
!       modules
!

use field_manager_mod,   only: fm_field_name_len
use field_manager_mod,   only: fm_get_index
use mpp_mod,             only: stdout, mpp_error, FATAL

use ocean_residency_meta_mod, only: ocean_residency_get_instances
use ocean_residency_meta_mod, only: ocean_residency_instance_type
use ocean_residency_meta_mod, only: ocean_residency_set_region_3d

use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_diag_tracer_type

!
!       force all variables to be "typed"
!

implicit none

!
!       Set all variables to be private by default

private

!
!       Public routines
!

public  :: ocean_residency_range_source
public  :: ocean_residency_range_start

!
!       Private routines
!

!
!       Public parameters
!


!
!       Private parameters
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_residency_range'
character(len=48), parameter                    :: mod_name = 'ocean_residency_range_mod'

!
!       Public variables
!

logical, public :: do_ocean_residency_range

!
!       Private types
!

type, private                                                   :: instance_extra_region_type
  integer                                                       :: tracer_index = 0
  logical                                                       :: is_prog_tracer = .true.
end type instance_extra_region_type

type, private                                                   :: instance_extra_type
  type(instance_extra_region_type), dimension(:), _ALLOCATABLE  :: region  _NULL
end type instance_extra_type

!
!       Private variables
!

type(instance_extra_type), dimension(:), allocatable            :: instance_extra
type(ocean_residency_instance_type), dimension(:), pointer      :: instance
integer                                                         :: num_instances = 0

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_range_source">
!
! <DESCRIPTION>
!       Calculate the mask for the tracer range(s)
!
!       Arguments:
!
!         Input:
!
!              isd:  low dimension of first index
!              ied:  high dimension of first index
!              jsd:  low dimension of second index
!              jed:  high dimension of second index
!               nk:  dimension of third index
!            taum1:  tau-1 time level index
!           T_prog:  array of ocean prognostic types
!           T_diag:  array of ocean diagnostic types
!         grid_kmt:  indices of the bottom grid cell
! </DESCRIPTION>
!

subroutine ocean_residency_range_source(isd, ied, jsd, jed, nk, taum1, T_prog, T_diag, grid_kmt)  !{

!
!       modules
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: taum1
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_diag_tracer_type), dimension(:), intent(inout)       :: T_diag
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt

!
!       local parameters
!

!
!       local variables
!

integer :: n
integer :: nn

!
!       set the source values for the residency tracers
!

do n = 1, num_instances  !{

  do nn = 1, instance(n)%num_regions  !{

    if (instance_extra(n)%region(nn)%tracer_index .ne. 0) then  !{

      if (instance_extra(n)%region(nn)%is_prog_tracer) then  !{
        call ocean_residency_set_region_3d(isd, ied, jsd, jed, nk, instance(n)%region(nn)%mask, &
             t_prog(instance_extra(n)%region(nn)%tracer_index)%field(:,:,:,taum1),              &
             instance(n)%region(nn)%params, grid_kmt,                                           &
             restore_region_value = 1.0, integrate_region_value = 0.0,                          &
             swap = instance(n)%region(nn)%swap_module)
      else  !}{
        call ocean_residency_set_region_3d(isd, ied, jsd, jed, nk, instance(n)%region(nn)%mask, &
             t_diag(instance_extra(n)%region(nn)%tracer_index)%field,                           &
             instance(n)%region(nn)%params, grid_kmt,                                           &
             restore_region_value = 1.0, integrate_region_value = 0.0,                          &
             swap = instance(n)%region(nn)%swap_module)
      endif  !}

    endif  !}

  enddo  !} nn

enddo  !} n

return

end subroutine ocean_residency_range_source   !}
! </SUBROUTINE> NAME="ocean_residency_range_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_range_start">
!
! <DESCRIPTION>
!       Start the ocean residency tracer range package
!
! </DESCRIPTION>
!
subroutine ocean_residency_range_start  !{

!
!       modules
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!


!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_range_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       local variables
!

integer :: i
integer :: n
integer :: nn

  integer :: stdoutunit 
  stdoutunit=stdout() 

instance => ocean_residency_get_instances(package_name)

write (stdoutunit,*)
if (.not. associated(instance)) then  !{
  write (stdoutunit,*) trim(note_header), ' No instances'
  num_instances = 0
  do_ocean_residency_range = .false.
else  !}{
  num_instances = size(instance)
  if (num_instances .eq. 1) then  !{
    write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instance'
  else  !}{
    write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instances'
  endif  !}
  do_ocean_residency_range = .true.
endif  !}

!
!       return if no tracers
!

if (.not. do_ocean_residency_range) then
  return
endif

!
!       Allocate array for extra instance information
!

allocate( instance_extra(num_instances) )
do n = 1, num_instances  !{
  allocate( instance_extra(n)%region(instance(n)%num_regions) )
enddo  !} n

do n = 1, num_instances  !{

  do nn = 1, instance(n)%num_regions  !{

!
!       Only process items for this module
!

    if (instance(n)%region(nn)%module_name .eq. package_name) then  !{

!
!       Check that a type was given
!

      if (.not. associated(instance(n)%region(nn)%strings)) then  !{
        call mpp_error(FATAL, trim(error_header) // ' No strings set for ' // trim(instance(n)%name))
      endif  !}


!
!       Check the types
!

      if (instance(n)%region(nn)%strings(1) .eq. 'tracer_range') then  !{

!
!       tracer_range: works with a region defined as
!               params(1) < T_prog(strings(2)) <= params(2)
!       if swap_module is true, then the region is inverted
!

        if (size(instance(n)%region(nn)%strings) .ne. 2) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Wrong number of strings set for ' //            &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (.not. associated(instance(n)%region(nn)%params)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' No params array set for ' //                    &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        elseif (mod(size(instance(n)%region(nn)%params),2) .ne. 0) then  !}{
          call mpp_error(FATAL, trim(error_header) // ' Odd number of params set for ' //               &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        else  !}{
          do i = 1, size(instance(n)%region(nn)%params), 2  !{
            if (instance(n)%region(nn)%params(i) .ge. instance(n)%region(nn)%params(i+1)) then  !{
              write (stdoutunit, *) trim(error_header), instance(n)%region(nn)%params(i), ' >= ',         &
                   instance(n)%region(nn)%params(i+1)
              call mpp_error(FATAL, trim(error_header) // ' Params not increasing for ' //              &
                   trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
            endif  !}
          enddo  !} i
        endif  !}

        if (associated(instance(n)%region(nn)%flags)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Flags set for ' // trim(instance(n)%name) //    &
               ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        instance_extra(n)%region(nn)%tracer_index = fm_get_index('/ocean_mod/prog_tracers/' //          &
             trim(instance(n)%region(nn)%strings(2)))
        if (instance_extra(n)%region(nn)%tracer_index .gt. 0) then  !{
          instance_extra(n)%region(nn)%is_prog_tracer = .true.
        else  !}{
          instance_extra(n)%region(nn)%tracer_index = fm_get_index('/ocean_mod/diag_tracers/' //        &
               trim(instance(n)%region(nn)%strings(2)))
          if (instance_extra(n)%region(nn)%tracer_index .gt. 0) then  !{
            instance_extra(n)%region(nn)%is_prog_tracer = .false.
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Could not get the index for tracer "' //      &
                 trim(instance(n)%region(nn)%strings(2)))
          endif  !}
        endif  !}

      else  !}{

!
!       Unknown type
!

        call mpp_error(FATAL, trim(error_header) // ' Unrecognized type (' //                           &
             trim(instance(n)%region(nn)%strings(1)) // ') for ' // trim(instance(n)%name))

      endif  !}

    endif  !}

  enddo  !} nn

enddo  !} n

return

end subroutine ocean_residency_range_start  !}
! </SUBROUTINE> NAME="ocean_residency_range_start"

end module ocean_residency_range_mod  !}
