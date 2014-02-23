module ocean_residency_ml_mod  !{
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
!       by a 2-d array of depths, typically mixed layer depths.
!
!       To use this module, in the field table namelist, you must do the following:
!
!       1) set "module_name" to "ocean_residency_ml"
!       2) "strings" should have one value: one of "kpp_bl", "mld_buoyancy", or
!          "mld_potrho"
!       3) no "params" should be set
!       4) no "flags" should be set
!
!       Currently, the following mixed layers are supported:
!
!          kpp_bl:  KPP mixed layer
!    mld_buoyancy:  mixed layer defined by a change in buoyancy
!                   (mld in diagnostic output)
!      mld_potrho:  mixed layer defined by a change in potential density
!                   (depth_of_potrho in diagnostic output)
!
!       There should be no "params" or "flags" set, and only one element of
!       "strings".
!
!       For an overview of the ocean residency modules, see ocean_residency.F90.
!
!</DESCRIPTION>
!
! $Id: ocean_residency_ml.F90,v 20.0 2013/12/14 00:17:10 fms Exp $
!

!
!       modules
!

use field_manager_mod,        only: fm_field_name_len, fm_get_index
use mpp_mod,                  only: stdout, mpp_error, FATAL

use ocean_residency_meta_mod, only: ocean_residency_get_instances
use ocean_residency_meta_mod, only: ocean_residency_instance_type
use ocean_residency_meta_mod, only: ocean_residency_set_region_2d
use ocean_types_mod,          only: ocean_thickness_type, ocean_prog_tracer_type
use ocean_types_mod,          only: ocean_time_type, ocean_density_type
use ocean_tracer_diag_mod,    only: calc_mixed_layer_depth, calc_potrho_mixed_layer

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

public  :: ocean_residency_ml_source
public  :: ocean_residency_ml_start

!
!       Private routines
!

!
!       Public parameters
!


!
!       Private parameters
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_residency_ml'
character(len=48), parameter                    :: mod_name = 'ocean_residency_ml_mod'

!
!       Public variables
!

logical, public :: do_ocean_residency_ml

!
!       Private types
!


!
!       Private variables
!

type(ocean_residency_instance_type), dimension(:), pointer      :: instance
integer                                                         :: num_instances = 0
integer                                                         :: indtemp, indsal

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_ml_source">
!
! <DESCRIPTION>
!       Calculate the mask for the mixed layer.
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
!           T_prog:  array of ocean prognostic types
!           T_diag:  array of ocean diagnostic types
!             Time:  ocean time type
!        Thickness:  ocean thickness type
!             Dens:  ocean density type
!        depth_zwt:  depth of bottom of grid cell
!       hblt_depth:  array of depths of KPP boundary layer
! </DESCRIPTION>
!

subroutine ocean_residency_ml_source(isd, ied, jsd, jed, nk,    &
     T_prog, Time, Thickness, Dens, depth_zwt, hblt_depth)  !{

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
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_time_type), intent(in)                               :: Time
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, dimension(isd:,jsd:,:), intent(in)                        :: depth_zwt
real, dimension(isd:,jsd:), intent(in)                          :: hblt_depth

!
!       local parameters
!

!
!       local variables
!

integer                                 :: n
integer                                 :: nn
real, dimension(isd:ied,jsd:jed)        :: mld

!
!       set the source values for the age tracers
!

do n = 1, num_instances  !{

  do nn = 1, instance(n)%num_regions  !{

    if (associated(instance(n)%region(nn)%strings)) then  !{

      if (instance(n)%region(nn)%strings(1) .eq. 'kpp_bl') then  !{
      
        mld(:,:) = hblt_depth(:,:)

      elseif (instance(n)%region(nn)%strings(1) .eq. 'mld_buoyancy') then  !}{

        call calc_mixed_layer_depth(Thickness,                                  &
             t_prog(indsal)%field(isd:ied,jsd:jed,:,time%tau),                  &
             t_prog(indtemp)%field(isd:ied,jsd:jed,:,time%tau),                 &
             dens%rho(isd:ied,jsd:jed,:,Time%tau),                              &
             dens%pressure_at_depth(isd:ied,jsd:jed,:),         &
             mld)

      elseif (instance(n)%region(nn)%strings(1) .eq. 'mld_potrho') then  !}{

        call calc_potrho_mixed_layer(thickness, dens,                     &
             potrho_mix_depth = mld)

      endif  !}

      call ocean_residency_set_region_2d(isd, ied, jsd, jed, nk, mld,           &
           instance(n)%region(nn)%mask, depth_zwt,                              &
           restore_region_value = 1.0, integrate_region_value = 0.0,            &
           swap = instance(n)%region(nn)%swap_module)

    endif  !}

  enddo  !} nn

enddo  !} n

return

end subroutine ocean_residency_ml_source   !}
! </SUBROUTINE> NAME="ocean_residency_ml_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_ml_start">
!
! <DESCRIPTION>
!       Start the ocean residency mixed layer package
!
! </DESCRIPTION>
!
subroutine ocean_residency_ml_start  !{

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

character(len=64), parameter    :: sub_name = 'ocean_residency_ml_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       local variables
!

integer :: n
integer :: nn

  integer :: stdoutunit 
  stdoutunit=stdout() 

instance => ocean_residency_get_instances(package_name)

write (stdoutunit,*)
if (.not. associated(instance)) then  !{
  write (stdoutunit,*) trim(note_header), ' No instances'
  num_instances = 0
  do_ocean_residency_ml = .false.
else  !}{
  num_instances = size(instance)
  if (num_instances .eq. 1) then  !{
    write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instance'
  else  !}{
    write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instances'
  endif  !}
  do_ocean_residency_ml = .true.
endif  !}

!
!       return if no tracers
!

if (.not. do_ocean_residency_ml) then
  return
endif

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
!       Whole grid boxes which contain at least part
!       of the mixed layer are used.
!       If swap_module is true, then the region is inverted
!

      if (instance(n)%region(nn)%strings(1) .eq. 'kpp_bl') then  !{

!
!       kpp_bl: KPP boundary layer
!

        if (size(instance(n)%region(nn)%strings) .ne. 1) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Wrong number of strings set for ' //            &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (associated(instance(n)%region(nn)%params)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Params set for ' //                             &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (associated(instance(n)%region(nn)%flags)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Flags set for ' //                              &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

      elseif (instance(n)%region(nn)%strings(1) .eq. 'mld_buoyancy') then  !}{

!
!       mld_buoyancy: mixed layer depth defined by a change in buoyancy
!               from the surface value.
!

        if (size(instance(n)%region(nn)%strings) .ne. 1) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Wrong number of region(nn)%strings set for ' // &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (associated(instance(n)%region(nn)%params)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Params set for ' //                             &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (associated(instance(n)%region(nn)%flags)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Flags set for ' //                              &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}


      elseif (instance(n)%region(nn)%strings(1) .eq. 'mld_potrho') then  !}{

!
!       mld_potrho: mixed layer depth defined by a change in potential
!               density from the surface value.
!

        if (size(instance(n)%region(nn)%strings) .ne. 1) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Wrong number of region(nn)%strings set for ' // &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (associated(instance(n)%region(nn)%params)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Params set for ' //                             &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
        endif  !}

        if (associated(instance(n)%region(nn)%flags)) then  !{
          call mpp_error(FATAL, trim(error_header) // ' Flags set for ' //                              &
               trim(instance(n)%name) // ' (' // trim(instance(n)%region(nn)%strings(1)) // ')')
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

indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then  !{
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif  !}
                                                                                                                          
indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then  !{
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif  !}

return

end subroutine ocean_residency_ml_start  !}
! </SUBROUTINE> NAME="ocean_residency_ml_start"

end module ocean_residency_ml_mod  !}
