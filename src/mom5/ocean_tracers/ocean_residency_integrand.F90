module ocean_residency_integrand_mod  !{
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
!       This module contains the subroutines to set up fields other
!       than reciprocal time over which to integrate.
!
!       To use this module, in the field table namelist, you must do the following:
!
!       1) set "module_name_int" to "ocean_residency_integrand"
!       2) "int_strings" should have one or two values, the first is
!          the name of a prognostic or diagnostic variable to use as the integrand,
!          and if the second value is given and is "average", then the array is assumed
!          to be defined at the top of the grid cells, and the average of the top
!          and bottom will be used for each grid cell (top of bottommost grid cell),
!          otherwise the value at each level will be used
!       3) "int_params" may be set to a single value to scale the integrand,
!          at most one value may be set
!       4) no "int_flags" should be set
!
!       For an overview of the ocean residency modules, see ocean_residency.F90.
!
!</DESCRIPTION>
!
! $Id: ocean_residency_integrand.F90,v 20.0 2013/12/14 00:17:06 fms Exp $
!

!
!       modules
!

use field_manager_mod,   only: fm_field_name_len
use field_manager_mod,   only: fm_get_index
use mpp_mod,             only: stdout, mpp_error, FATAL
use fms_mod,             only: lowercase

use ocean_residency_meta_mod, only: ocean_residency_get_instances
use ocean_residency_meta_mod, only: ocean_residency_instance_type

use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_diag_tracer_type

!
!       force all variables to be "typed"
!

implicit none

!
!       Set all variables to be private by default
!

private

!
!       Public routines
!

public  :: ocean_residency_integrand_source
public  :: ocean_residency_integrand_start

!
!       Private routines
!

!
!       Public parameters
!


!
!       Private parameters
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_residency_integrand'
character(len=48), parameter                    :: mod_name = 'ocean_residency_integrand_mod'

!
!       Public variables
!

logical, public :: do_ocean_residency_integrand

!
!       Private types
!

type, private                           :: instance_extra_type
  integer                               :: tracer_index = 0
  logical                               :: is_prog_tracer = .true.
  logical                               :: average = .true.
  real                                  :: scale = 1.0
end type instance_extra_type

!
!       Private variables
!

type(instance_extra_type), dimension(:), allocatable            :: instance_extra
type(ocean_residency_instance_type), dimension(:), pointer      :: instance
integer                                                         :: num_instances = 0

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_integrand_source">
!
! <DESCRIPTION>
!       Set the selected integrand field to be used in the "integrate_region"
!       as defined in ocean_residency.F90. The selected prognostic
!       or diagnostic variable will be multiplied by the residency mask
!       (which is usually either 1 or 0, but is not so required) and a
!       user-selected scale factor (default of 1). The variable may be either
!       the average of the surrounding points ( C(k) = (C(k) + C(k+1))/2 )
!       or just the same indexical value.
! </DESCRIPTION>
!

subroutine ocean_residency_integrand_source(isd, ied, jsd, jed, nk, taum1, T_prog, T_diag, grid_tmask)  !{

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
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask

!
!       local parameters
!

!
!       local variables
!

integer                                 :: tracer_index
integer                                 :: i
integer                                 :: j
integer                                 :: k
integer                                 :: n

!
!       set the source values for the residency tracers
!

do n = 1, num_instances  !{

  if (instance(n)%integrate_region_value .ne. 0.0) then  !{

    tracer_index = instance_extra(n)%tracer_index

    if (instance_extra(n)%is_prog_tracer) then  !{

!
!       The variable is a prognostic variable
!

      if (instance_extra(n)%average) then  !{

        do k = 1, nk - 1  !{
          do j = jsd, jed  !{
            do i = isd, ied  !{
              if (instance(n)%index(i,j,k) .eq. 0) then  !{
                instance(n)%mask(i,j,k) = instance(n)%mask(i,j,k) * instance_extra(n)%scale *   &
                     (T_prog(tracer_index)%field(i,j,k,taum1) * grid_tmask(i,j,k) +             &
                      T_prog(tracer_index)%field(i,j,k+1,taum1) * grid_tmask(i,j,k+1)) * 0.5
              endif  !}
            enddo  !} i
          enddo  !} j
        enddo  !} k
        do j = jsd, jed  !{
          do i = isd, ied  !{
            if (instance(n)%index(i,j,k) .eq. 0) then  !{
              instance(n)%mask(i,j,nk) = instance(n)%mask(i,j,nk) * instance_extra(n)%scale *   &
                   T_prog(tracer_index)%field(i,j,nk,taum1) * grid_tmask(i,j,nk) * 0.5
            endif  !}
          enddo  !} i
        enddo  !} j

      else  !}{

        do k = 1, nk  !{
          do j = jsd, jed  !{
            do i = isd, ied  !{
              if (instance(n)%index(i,j,k) .eq. 0) then  !{
                instance(n)%mask(i,j,k) = instance(n)%mask(i,j,k) * instance_extra(n)%scale *   &
                     T_prog(tracer_index)%field(i,j,k,taum1) * grid_tmask(i,j,k)
              endif  !}
            enddo  !} i
          enddo  !} j
        enddo  !} k

      endif  !}

    else  !}{

!
!       The variable is a diagnostic variable
!

      if (instance_extra(n)%average) then  !{

        do k = 1, nk - 1  !{
          do j = jsd, jed  !{
            do i = isd, ied  !{
              if (instance(n)%index(i,j,k) .eq. 0) then  !{
                instance(n)%mask(i,j,k) = instance(n)%mask(i,j,k) * instance_extra(n)%scale *   &
                     (T_diag(tracer_index)%field(i,j,k) * grid_tmask(i,j,k) +                   &
                      T_diag(tracer_index)%field(i,j,k+1) * grid_tmask(i,j,k+1)) * 0.5
              endif  !}
            enddo  !} i
          enddo  !} j
        enddo  !} k
        do j = jsd, jed  !{
          do i = isd, ied  !{
            if (instance(n)%index(i,j,k) .eq. 0) then  !{
              instance(n)%mask(i,j,nk) = instance(n)%mask(i,j,nk) * instance_extra(n)%scale *   &
                   T_diag(tracer_index)%field(i,j,nk) * grid_tmask(i,j,nk) * 0.5
            endif  !}
          enddo  !} i
        enddo  !} j

      else  !}{

        do k = 1, nk  !{
          do j = jsd, jed  !{
            do i = isd, ied  !{
              if (instance(n)%index(i,j,k) .eq. 0) then  !{
                instance(n)%mask(i,j,k) = instance(n)%mask(i,j,k) * instance_extra(n)%scale *   &
                     T_diag(tracer_index)%field(i,j,k) * grid_tmask(i,j,k)
              endif  !}
            enddo  !} i
          enddo  !} j
        enddo  !} k

      endif  !}

    endif  !}

  endif  !}

enddo  !} n

return

end subroutine ocean_residency_integrand_source  !}
! </SUBROUTINE> NAME="ocean_residency_integrand_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_integrand_start">
!
! <DESCRIPTION>
!       Start the ocean residency integrand package
!
! </DESCRIPTION>
!
subroutine ocean_residency_integrand_start  !{

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

character(len=64), parameter    :: sub_name = 'ocean_residency_integrand_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       local variables
!

integer :: n

  integer :: stdoutunit 
  stdoutunit=stdout() 

instance => ocean_residency_get_instances(package_name, integrand = .true.)

write (stdoutunit,*)
if (.not. associated(instance)) then  !{
  write (stdoutunit,*) trim(note_header), ' No instances'
  num_instances = 0
  do_ocean_residency_integrand = .false.
else  !}{
  num_instances = size(instance)
  if (num_instances .eq. 1) then  !{
    write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instance'
  else  !}{
    write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instances'
  endif  !}
  do_ocean_residency_integrand = .true.
endif  !}

!
!       return if no tracers
!

if (.not. do_ocean_residency_integrand) then
  return
endif

!
!       Allocate array for extra instance information
!

allocate( instance_extra(num_instances) )

do n = 1, num_instances  !{

!
!       Check that a type was given
!

  if (.not. associated(instance(n)%int_strings)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' No int_strings set for ' // trim(instance(n)%name))
  endif  !}

!
!       Check the types
!

  if (associated(instance(n)%int_strings)) then  !{

!
!       Define the tracer field to use as the integrand.
!       The first int_strings has the name of the field to use, and the
!       second int_strings, if given, sets whether to do a vertical average of
!       the field or just use the field. By default, averaging is not done.
!
!       An optional int_param may be given to scale to field. If not given,
!       then a value of one is used.
!

    if (size(instance(n)%int_strings) .gt. 2) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Wrong number of int_strings set for ' //    &
           trim(instance(n)%name))
    elseif (size(instance(n)%int_strings) .eq. 1) then  !}{
      instance_extra(n)%average = .false.
    else  !}{
      instance_extra(n)%average = lowercase(instance(n)%int_strings(2)) .eq. 'average'
    endif  !}

    if (associated(instance(n)%int_params)) then  !{
      if (size(instance(n)%int_params) .ne. 1) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Too many int_params set for ' //          &
             trim(instance(n)%name) // ' (' // trim(instance(n)%int_strings(1)) // ')')
      else  !}{
        instance_extra(n)%scale = instance(n)%int_params(1)
      endif  !}
    else  !}{
      instance_extra(n)%scale = 1.0
    endif  !}

    if (associated(instance(n)%int_flags)) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Flags set for ' //                          &
           trim(instance(n)%name) // ' (' // trim(instance(n)%int_strings(1)) // ')')
    endif  !}

    instance_extra(n)%tracer_index = fm_get_index('/ocean_mod/prog_tracers/' //                 &
         trim(instance(n)%int_strings(1)))
    if (instance_extra(n)%tracer_index .gt. 0) then  !{
      instance_extra(n)%is_prog_tracer = .true.
    else  !}{
      instance_extra(n)%tracer_index = fm_get_index('/ocean_mod/diag_tracers/' //               &
           trim(instance(n)%int_strings(1)))
      if (instance_extra(n)%tracer_index .gt. 0) then  !{
        instance_extra(n)%is_prog_tracer = .false.
      else  !}{
        call mpp_error(FATAL, trim(error_header) // ' Could not get the index for tracer "' //  &
             trim(instance(n)%int_strings(1)))
      endif  !}
    endif  !}

  endif  !}

enddo  !} n

return

end subroutine ocean_residency_integrand_start  !}
! </SUBROUTINE> NAME="ocean_residency_integrand_start"

end module ocean_residency_integrand_mod  !}
