module ocean_residency_mod  !{
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Ocean residency module
!</OVERVIEW>
!
!<DESCRIPTION>
!
!       This module is a superset of the ocean age tracer module. It may
!       be used to reproduce the age tracer, but it can also do much
!       more. Unlike the ocean age tracer module, here you may specify
!       a 3-d field by specifying a series of "rectangular prisms".
!       The grid cells which occupy this field may vary with time due
!       to the variations in the thickness of the surface layer.
!
!       You may also specify a 3-d field by choosing one of three mixed-layers
!       (KPP, density-derived or buoyancy-derived). You may also specify a
!       region based on a a range of a prognostic or diagnostic ocean
!       tracer (such as all points with a temperature of 10-20 degrees).
!
!       You may specify either the mixed-layer method or the tracer range along
!       with the geographic specification. If multiple methods are used, then 
!       the resultant field is the intersection of the two methods. In fact, the
!       geographic method is always in use, and defaults to the whole ocean.
!
!       There is an option of using the inverse of any method. This is sometimes
!       easier to use than explicitly specifying the inverse. For instance,
!       to get the temperatures outside of the range 10-20 degrees, one could
!       specify the 10-20 degree range and set <swap> to true, or, specify two
!       different regions, one less than 10 degrees and one greater than 20 degrees.
!
!       By default, the values inside the specified region are set to a
!       specified value each time-step (default is 0), and outside of this region
!       the field is integrated over time in units of years (integrand is
!       1/(365.25*86400) by default, but this value can be changed. The inner region
!       can be forced to be 0 (or the user-specified inner value) at each time-step
!       (the default), or it can be restored to this value at a user-specified
!       rate (given in days), or, it can be left alone (not integrated or set to
!       any value).
!       
!       Finally, since we are just integrating a specified field by a constant value
!       (by default), it is simple to make that field a simple function of another
!       tracer. This is the final option. You may specify a prognostic or diagnostic
!       variable to be the integrand, and have it scaled by a constant. One case
!       which has been done is to integrate irradiance in the mixed-layer.
!       
!       This module is split into several different modules. This particular module
!       is the only one called from outside, and it makes appropriate calls to 
!       the other modules to implement its features. This was done so that it would
!       be easier to expand the features of the ocean residency, without cluttering
!       the code too much. Ideally, this could all be done with classes, but there
!       is currently no support for classes in Fortran (pre-F2003). The current
!       modules are:
!
!       ocean_residency.F90             This module, the control center, and also does the
!                                       integrations and resetting in the inner regions.
!       ocean_residency_meta.F90        Has specifications of the ocean_residency type
!                                       and other utility routines. This is separate
!                                       to stop circular references
!       ocean_residency_ml.F90          Handles setting mixed-layer masks
!       ocean_residency_range.F90       Handles setting the masks for the tracer ranges
!       ocean_residency_integrand.F90   Handles setting the integrand to a non-constant
!                                       value
!
!       field_table namelist inputs:
!
!       The following 6 arrays specify the bounds of boxes, the intersection of which
!       will define the region where the integrand is set to restore_region_value. ALl arrays
!       must have the same number of elements. If none are specified, then no region will
!       be used, unless swap is true, which is a quick way of selecting the entire ocean.
!       This may be useful when one really wants to use one of the other methods of selecting,
!       such as mixed layer depth, or tracer range.
!
!       Longitudes will be shifted to lie in the range 0-360 degrees. If the eastern side is
!       greater than the western, then the selected region will consist of those grid cells from
!       the eastern value to 360, and 0 to the western value. If the northern value is less than
!       the southern value, or if the top depth is greater than the bottom depth, then
!       an error occurs, and the model stops.
!
!       Three special cases exist for the depth (which currently must be in meters).
!       If the bottom value is less than or equal to zero, then the top box is selected. 
!       If the top value is greater than the maximum depth, then the bottom box is selected.
!       If the top value is negative, then grid cells within "the absolute value of the top value"  
!       from the bottom are selected.
!
!       Note that the geographic specification is used for all residency tracers, so either
!       values for these 6 arrays must be given so as to select a region, or else swap must
!       be set to true and no arrays given.
!
!               east_bnd: array of boundary points of the eastern side of the box,
!                         in degrees longitude (default: NULL)
!              north_bnd: array of boundary points of the northern side of the box,
!                         in degrees latitude (default: NULL)
!              south_bnd: array of boundary points of the southern side of the box,
!                         in degrees latitude (default: NULL)
!               west_bnd: array of boundary points of the western side of the box,
!                         in degrees longitude (default: NULL)
!                top_bnd: array of boundary points of the top side of the box,
!                         in meters (default: NULL)
!             bottom_bnd: array of boundary points of the bottom side of the box,
!                         in meters (default: NULL)
!
!                   swap: if true, then select the inverse of the specified geographic region,
!                         otherwise, just use the specified region (default: false)
!                restore: restoring value for values in the defined regions
!                         negative => do nothing, 0 => force to integrate_region_value,
!                         positive => time scale in days to force to integrate_region_value
!                         (default: 0.0)
!   restore_region_value: value to set the mask to for grid cells within
!                         the specified region (default: 0.0)
! integrate_region_value: value to set the mask to for grid cells outside
!                         the specified region (default: secs_in_year_r)
!            swap_module: if true, then select the inverse of the region from the specified module,
!                         otherwise, just use the specified region (default: false)
!            module_name: if set, then it will be used to select the alternate
!                         method of selecting the region (default: ' ' -- only geographic
!                         selection is used)
!
!       For the following arrays, see the different extra modules for required and
!       possible values for the module selected.
!
!                 params: an array of real parameters which may be used for the method
!                         being used for selection (default: NULL)
!                  flags: an array of real parameters which may be used for the method
!                         being used for selection (default: NULL)
!                strings: an array of real parameters which may be used for the method
!                         being used for selection (default: NULL)
!
!       For the following arrays, see the different extra modules for required and
!       possible values for the module selected to specify the integrand.

!        int_module_name: if set, then the name of the module used to set
!                         the integrand (default: ' ' -- integrate time, in years)
!             int_params: an array of real parameters which may be used for the integrand
!                         (default: NULL)
!              int_flags: an array of real parameters which may be used for the integrand
!                         (default: NULL)
!            int_strings: an array of real parameters which may be used for the integrand
!                         (default: NULL)
!
!----------------------------------------------------------------------------------------
!
!       Sample field table entries:
!       ---------------------------
!
! "tracer_packages","ocean_mod","ocean_residency"
! names = age_surface, age_bottom_inv, kppbl_nil, kppbl_14d, kppbl_frc, kppbl_irr_14d, temp_15_20
! horizontal-advection-scheme = mdppm
! vertical-advection-scheme = mdppm
! units = yr
! min_tracer_limit=0.0
! /
!
!       This is the same as the old age tracer with all surface
!       values forced to zero
! 
! "namelists","ocean_mod","ocean_residency/age_surface"
! south_bnd = -90.0
! north_bnd =  90.0
! west_bnd =   0.0
! east_bnd = 360.0
! top_bnd = 0.0
! bottom_bnd = 0.0
! /
!
!       This integrates the age in the bottom box and forces to
!       zero everywhere else (note that swap is true)
! 
! "namelists","ocean_mod","ocean_residency/age_bottom_inv"
! south_bnd = -90.0
! north_bnd =  90.0
! west_bnd =   0.0
! east_bnd = 360.0
! top_bnd = 10000.0
! bottom_bnd = 10000.0
! swap = t
! /
!
!       This integrates the age of the water in the KPP
!       boundary layer, but lets the age outside of this region
!       keep its value until it again mixes with the boundary layer.
!       The module_name is set, and strings is set to pick the
!       type of mixed layer desired
!       (note that the global geographic area is explicitly specified,
!       and that the restoring time scale is negative)
! 
! "namelists","ocean_mod","ocean_residency/kppbl_nil"
! south_bnd = -90.0
! north_bnd =  90.0
! west_bnd =   0.0
! east_bnd = 360.0
! top_bnd = 0.0
! bottom_bnd = 10000.0
! restore = -1.0
! module_name = ocean_residency_ml
! strings = kpp_bl
! swap_module = t
! /
!
!       This is the same as above, but forces the age outside
!       the boundary layer to 0 with a 14 day time scale
!       (note that the global geographic region is specified by
!       setting swap to true, also note the value of restore)
! 
! "namelists","ocean_mod","ocean_residency/kppbl_14d"
! swap = t
! restore = 14.0
! module_name = ocean_residency_ml
! strings = kpp_bl
! swap_module = t
! /
!
!       This is the same as above, but forces age to zero outside
!       the boundary layer (note that restore did not need to
!       be explicitly specified, as the default is zero)
! 
! "namelists","ocean_mod","ocean_residency/kppbl_frc"
! swap = t
! restore = 0.0
! module_name = ocean_residency_ml
! strings = kpp_bl
! swap_module = t
! /
!
!       The following integrates irradiance in the boundary
!       layer (note that the units needed to be changed for
!       netCDF output purposes)
! 
! "prog_tracers","ocean_mod","residency_kppbl_irr_14d"
! units = W-yr/m^2
! /
! 
! "namelists","ocean_mod","ocean_residency/kppbl_irr_14d"
! swap = t
! restore = 14.0
! module_name = ocean_residency_ml
! strings = kpp_bl
! swap_module = t
! int_module_name = ocean_residency_integrand
! int_strings = irr
! /
!
!       This specifies the region as the area with
!       a temperature range of between 15 and 20 degrees
!       (note that the params holds the variable
!       range)
! 
! "namelists","ocean_mod","ocean_residency/temp_15_20"
! swap = t
! module_name = ocean_residency_range
! strings = tracer_range, temp
! params = 15.0, 20.0
! /
!
!</DESCRIPTION>
!
! $Id: ocean_residency.F90,v 20.0 2013/12/14 00:17:04 fms Exp $
!

!
!       modules
!

use time_manager_mod,         only: time_type
use field_manager_mod,        only: fm_string_len, fm_path_name_len, fm_field_name_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use mpp_mod,                  only: stdout, mpp_error, FATAL, mpp_pe, mpp_root_pe
use diag_manager_mod,         only: register_diag_field, send_data

use ocean_tpm_util_mod,       only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer
use fm_util_mod,              only: fm_util_set_value
use fm_util_mod,              only: fm_util_start_namelist, fm_util_end_namelist
use fm_util_mod,              only: fm_util_get_string, fm_util_get_string_array
use fm_util_mod,              only: fm_util_get_logical, fm_util_get_logical_array
use fm_util_mod,              only: fm_util_get_integer, fm_util_get_real, fm_util_get_real_array
use fm_util_mod,              only: fm_util_check_for_bad_fields

use ocean_types_mod,               only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,               only: ocean_thickness_type, ocean_time_type, ocean_density_type
use ocean_residency_meta_mod,      only: ocean_residency_set_region_geog
use ocean_residency_meta_mod,      only: instance, num_instances, secs_in_year_r, sec_per_day
use ocean_residency_ml_mod,        only: do_ocean_residency_ml
use ocean_residency_ml_mod,        only: ocean_residency_ml_start, ocean_residency_ml_source
use ocean_residency_range_mod,     only: do_ocean_residency_range
use ocean_residency_range_mod,     only: ocean_residency_range_start, ocean_residency_range_source
use ocean_residency_integrand_mod, only: do_ocean_residency_integrand
use ocean_residency_integrand_mod, only: ocean_residency_integrand_start, ocean_residency_integrand_source

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

public  :: ocean_residency_init
public  :: ocean_residency_source
public  :: ocean_residency_start
public  :: ocean_residency_tracer

!
!       Private routines
!

!
!       Public parameters
!


!
!       Private parameters
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_residency'
character(len=48), parameter                    :: mod_name = 'ocean_residency_mod'
character(len=48), parameter                    :: diag_name = 'ocean_residency'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocean_residency.res.nc'

!
!       Public variables
!

logical, public :: do_ocean_residency

!
!       Public types
!


!
!       Private variables
!

integer :: package_index

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
! </DESCRIPTION>
!

subroutine ocean_residency_init  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer                                                 :: n
integer                                                 :: nn
integer                                                 :: num_regions
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_field_name_len+3)                      :: long_suffix
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

  integer :: stdoutunit 
  stdoutunit=stdout() 

!
!       Initialize the ocean residency package
!

package_index = otpm_set_tracer_package(package_name,                           &
     units = 'yr', flux_units = 'm',                                            &
     min_tracer_limit = 0.0, max_tracer_limit = 1.0e+20,                        &
     restart_file = default_restart_file,                                       &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!
!       Check whether to use this package
!

path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
num_instances = fm_get_length(path_to_names)
if (num_instances .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif  !}

!
!       Check some things
!

write (stdoutunit,*)
if (num_instances .eq. 0) then  !{
  write (stdoutunit,*)                                            &
       trim(note_header), ' No instances'
  do_ocean_residency = .false.
else  !}{
  write (stdoutunit,*) trim(note_header), ' ', num_instances, ' instances'
  do_ocean_residency = .true.
endif  !}

!
!       Return if we don't want to use this package
!

if (.not. do_ocean_residency) then  !{
  return
endif  !}

!
!       allocate the instance array
!

allocate ( instance(num_instances) )

!
!       loop over the names, saving them into the instance array
!

do n = 1, num_instances  !{

  if (fm_get_value(path_to_names, name, index = n)) then  !{
    instance(n)%name = name
  else  !}{
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         ' Bad field name for index ' // trim(name))
  endif  !}

enddo  !} n

do n = 1, num_instances  !{

!
!       determine the tracer name for this instance
!

  name = instance(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

  instance(n)%tracer_index = otpm_set_prog_tracer('residency' // trim(suffix), package_name,    &
       longname = 'Residency' // trim(long_suffix), units = 'yr',                               &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')

enddo  !} n

!
!       set up the density diagnostic tracers
!

n = otpm_set_diag_tracer('rho_in_situ', longname = 'In situ density', units = 'kg/m^3',         &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
n = otpm_set_diag_tracer('rho_neutral', longname = 'Neutral density', units = 'kg/m^3',         &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
n = otpm_set_diag_tracer('rho_pot_0', longname = 'Potential density-0', units = 'kg/m^3',       &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
n = otpm_set_diag_tracer('rho_pot_1', longname = 'Potential density-1', units = 'kg/m^3',       &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
n = otpm_set_diag_tracer('rho_pot_2', longname = 'Potential density-2', units = 'kg/m^3',       &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
n = otpm_set_diag_tracer('rho_pot_3', longname = 'Potential density-3', units = 'kg/m^3',       &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
n = otpm_set_diag_tracer('rho_pot_4', longname = 'Potential density-4', units = 'kg/m^3',       &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!
!-----------------------------------------------------------------------
!       Set up the instance residency namelists
!-----------------------------------------------------------------------
!

!
!       Add the package name to the list of good namelists, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif  !}

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, num_instances  !{

!
!       create the instance namelist
!

  call fm_util_start_namelist(package_name, instance(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('num_regions', 1)
  call fm_util_set_value('integrate_region_value', secs_in_year_r)
  call fm_util_set_value('union', .true.)
  call fm_util_set_value('int_module_name', ' ')
  call fm_util_set_value('int_params', 0.0, index = 0)
  call fm_util_set_value('int_flags', .false., index = 0)
  call fm_util_set_value('int_strings', ' ', index = 0)

!
!       get the number of regions so that we may set up the extra namelists
!

  num_regions =  fm_util_get_integer('num_regions', scalar = .true.)

!
!       create namelists for each region
!

  if (num_regions .le. 0) then
    call mpp_error(FATAL,trim(error_header) // ' num_regions is non-positive for instance ' //   &
         trim(instance(n)%name) // ' in package ' // trim(package_name))
  endif

  do nn = 1, num_regions  !{

    if (num_regions .eq. 1) then
      suffix = ' '
    else
      write (long_suffix,*) nn
      suffix = '_' // long_suffix(scan(long_suffix,'0123456789'):)
    endif

    call fm_util_set_value('restore' // suffix, 0.0)
    call fm_util_set_value('restore_region_value' // suffix, 0.0)
    call fm_util_set_value('swap' // suffix, .false.)
    call fm_util_set_value('swap_module' // suffix, .false.)
    call fm_util_set_value('west_bnd' // suffix, 0.0, index = 0)
    call fm_util_set_value('east_bnd' // suffix, 0.0, index = 0)
    call fm_util_set_value('south_bnd' // suffix, 0.0, index = 0)
    call fm_util_set_value('north_bnd' // suffix, 0.0, index = 0)
    call fm_util_set_value('top_bnd' // suffix, 0.0, index = 0)
    call fm_util_set_value('bottom_bnd' // suffix, 0.0, index = 0)
    call fm_util_set_value('module_name' // suffix, ' ')
    call fm_util_set_value('params' // suffix, 0.0, index = 0)
    call fm_util_set_value('flags' // suffix, .false., index = 0)
    call fm_util_set_value('strings' // suffix, ' ', index = 0)

  enddo  !} nn

  call fm_util_end_namelist(package_name, instance(n)%name, check = .true., caller = caller_str)

enddo  !} n

!
!       Check for any errors in the number of fields in the namelists for this package
!

good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',     &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then  !{
  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,                   &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif  !}

return

end subroutine ocean_residency_init  !}
! </SUBROUTINE> NAME="ocean_residency_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_source">
!
! <DESCRIPTION>
!       Calculate the source arrays for the tracer packages
! </DESCRIPTION>
!

subroutine ocean_residency_source(isc, iec, jsc, jec, isd, ied, jsd, jed, nk,   &
     T_prog, T_diag, Time, Thickness, Dens, grid_xt, grid_yt, grid_zw,          &
     grid_tmask, grid_kmt, hblt_depth)

!
!       modules
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
integer, intent(in)                                             :: nk
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_diag_tracer_type), dimension(:), intent(inout)       :: T_diag
type(ocean_time_type), intent(in)                               :: Time
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, dimension(isd:,jsd:), intent(in)                          :: grid_xt
real, dimension(isd:,jsd:), intent(in)                          :: grid_yt
real, dimension(:), intent(in)                                  :: grid_zw
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt
real, intent(in), dimension(isd:,jsd:)                          :: hblt_depth

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_source'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       local variables
!

integer :: i
integer :: ind
integer :: j
integer :: k
integer :: n
integer :: nn
logical :: used
logical :: good

  integer :: stdoutunit 
  stdoutunit=stdout() 

!
!       set the source values for the residency tracers
!

do n = 1, num_instances  !{

!
!       set the values via the input values
!

  do nn = 1, instance(n)%num_regions  !{

    call ocean_residency_set_region_geog(isd, ied, jsd, jed, nk, instance(n)%region(nn)%mask,   &
         grid_xt, grid_yt, grid_zw(nk), Thickness%depth_zt, Thickness%depth_zwt,                &
         instance(n)%region(nn)%num_geog_regions,                                               &
         instance(n)%region(nn)%west_bnd, instance(n)%region(nn)%east_bnd,                      &
         instance(n)%region(nn)%south_bnd, instance(n)%region(nn)%north_bnd,                    &
         instance(n)%region(nn)%top_bnd, instance(n)%region(nn)%bottom_bnd,                     &
         grid_kmt, t_prog(instance(n)%tracer_index)%name,                                       &
         restore_region_value = 1.0, integrate_region_value = 0.0,                              &
         swap = instance(n)%region(nn)%swap,                                                    &
         initialize = .true.)

  enddo  !} nn

enddo  !} n

!
!       set the source values for the residency tracers which need to be set in other
!       modules
!

if (do_ocean_residency_ml) then  !{
  call ocean_residency_ml_source(isd, ied, jsd, jed, nk,                &
       t_prog, time, thickness, dens, Thickness%depth_zwt, hblt_depth)
endif  !}

if (do_ocean_residency_range) then  !{
  call ocean_residency_range_source(isd, ied, jsd, jed, nk, Time%taum1, t_prog, t_diag, grid_kmt)
endif  !}

!
!       save out the mask fields for each region
!

do n = 1, num_instances  !{

  do nn = 1, instance(n)%num_regions  !{

    do k = 1, nk
      do j = jsd, jed
        do i = isd, ied
          instance(n)%region(nn)%mask(i,j,k) = instance(n)%region(nn)%mask(i,j,k) * grid_tmask(i,j,k)
        enddo
      enddo
    enddo

    if (instance(n)%region(nn)%id_restore_region .gt. 0) then  !{
      used = send_data(instance(n)%region(nn)%id_restore_region,        &
           instance(n)%region(nn)%mask(:,:,:),                          &
           Time%model_time, rmask = grid_tmask(:,:,:),                  &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif  !}

  enddo  !} nn

enddo  !} n

!
!       merge the mask fields into one field
!

do n = 1, num_instances  !{

  if (instance(n)%num_regions .eq. 1) then  !{

    do k = 1, nk  !{
      do j = jsd, jed  !{
        do i = isd, ied  !{
          if (grid_tmask(i,j,k) .lt. 0.5) then  !{
            instance(n)%mask(i,j,k) = 0.0
            instance(n)%index(i,j,k) = -1
          elseif (instance(n)%region(1)%mask(i,j,k) .eq. 1.0) then  !}{
            instance(n)%mask(i,j,k) = instance(n)%region(1)%restore_region_value
            instance(n)%index(i,j,k) = 1
          else  !}{
            instance(n)%mask(i,j,k) = instance(n)%integrate_region_value
            instance(n)%index(i,j,k) = 0
          endif  !}
        enddo  !} i
      enddo  !} j
    enddo  !} k

  else  !}{

    do k = 1, nk  !{
      do j = jsd, jed  !{
        do i = isd, ied  !{
          if (grid_tmask(i,j,k) .gt. 0.5) then  !{

!
!       we're looking to set the restore_region_value for any point that has at least
!       one region set
!

            if (instance(n)%union) then  !{

!
!       find the first region with a value set
!

              ind = 0
              do nn = 1, instance(n)%num_regions  !{
                if (instance(n)%region(nn)%mask(i,j,k) .eq. 1.0) then  !{
                  ind = nn
                  exit
                endif  !}
              enddo  !} nn

!
!       if ind == 0, then this point is out of all regions
!

              if (ind .ne. 0) then  !{

!
!       Quit with an error if a
!       region with a different restore_region_value also has a value set
!

                do nn = ind + 1, instance(n)%num_regions  !{
                  if (instance(n)%region(nn)%mask(i,j,k) .eq. 1.0) then  !{
                    if (instance(n)%region(nn)%restore_region_value .ne.                                &
                        instance(n)%region(ind)%restore_region_value) then  !{
                      write (stdoutunit,*) trim(error_header), ' Grid point ', i, ', ', j, ', ',k,        &
                           ' in two regions for "', trim(instance(n)%name), '"'
                      call mpp_error(FATAL,trim(error_header) // ' Grid point in two regions for "' //  &
                           trim(instance(n)%name) // '"')
                    elseif (instance(n)%region(nn)%restore .ne. instance(n)%region(ind)%restore) then  !}{
                      write (stdoutunit,*) trim(error_header), ' Grid point ', i, ', ', j, ', ',k,        &
                           ' has different restore values for "', trim(instance(n)%name), '"'
                      call mpp_error(FATAL,trim(error_header) // ' Grid point has different restore values for "' //    &
                           trim(instance(n)%name) // '"')
                    endif  !}
                  endif  !}
                enddo  !} nn
              endif  !}

            else  !}{     ! not a union, therefore an intersection

!
!       Find the first region with a value set
!

              ind = 0
              do nn = 1, instance(n)%num_regions  !{
                if (instance(n)%region(nn)%mask(i,j,k) .eq. 1.0) then  !{
                  ind = nn
                  exit
                endif  !}
              enddo  !} nn

!
!       If ind == 0, then this point is out of all regions
!

              if (ind .gt. 0) then  !{

!
!       Otherwise, check whether this value is set for all regions with
!       the same restore_region_value. If so, use that value, otherwise
!       use the integrate_region_value. Also, quit with an error if a
!       region with a different restore_region_value also has a value set
!

                good = .true.
                do nn = ind + 1, instance(n)%num_regions  !{
                  if (instance(n)%region(nn)%mask(i,j,k) .eq. 1.0) then  !{
                    if (instance(n)%region(nn)%restore_region_value .ne.                                &
                        instance(n)%region(ind)%restore_region_value) then  !{
                      write (stdoutunit,*) trim(error_header), ' Grid point ', i, ', ', j, ', ',k,        &
                           ' in two regions for "', trim(instance(n)%name), '"'
                      call mpp_error(FATAL,trim(error_header) // ' Grid point in two regions for "' //  &
                           trim(instance(n)%name) // '"')
                    elseif (instance(n)%region(nn)%restore_region_value .ne.                            &
                            instance(n)%region(ind)%restore_region_value) then  !}{
                      write (stdoutunit,*) trim(error_header), ' Grid point ', i, ', ', j, ', ',k,        &
                           ' has different restore values for "', trim(instance(n)%name), '"'
                      call mpp_error(FATAL,trim(error_header) // ' Grid point has different restore values for "' //    &
                           trim(instance(n)%name) // '"')
                    else  !}{
                      good = good .and. instance(n)%region(nn)%mask(i,j,k) .eq. 1.0
                    endif  !}
                  endif  !}
                enddo  !} nn
                if (.not. good) then  !{
                  ind = 0
                endif  !}
              endif  !}

            endif  !}

            instance(n)%index(i,j,k) = ind
            if (ind .gt. 0) then
              instance(n)%mask(i,j,k) = instance(n)%region(ind)%restore_region_value
            else
              instance(n)%mask(i,j,k) = instance(n)%integrate_region_value
            endif

          else  !}{

            instance(n)%index(i,j,k) = -1
            instance(n)%mask(i,j,k) = 0.0

          endif  !}
        enddo  !} i
      enddo  !} j
    enddo  !} k

  endif  !}

enddo  !} n

!
!       set the integrate_region value if something other than time is to be integrated
!

if (do_ocean_residency_integrand) then  !{
  call ocean_residency_integrand_source(isd, ied, jsd, jed, nk, Time%taum1, t_prog, t_diag, grid_tmask)
endif  !}

do n = 1, num_instances  !{

  do k = 1, nk  !{
    do j = jsd, jed  !{
      do i = isd, ied  !{
        t_prog(instance(n)%tracer_index)%source(i,j,k) = t_prog(instance(n)%tracer_index)%source(i,j,k) +   &
             instance(n)%mask(i,j,k)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!       Save the diagnostic for the merged mask array
!

  if (instance(n)%id_restore_region .gt. 0) then  !{
    used = send_data(instance(n)%id_restore_region,           &
         instance(n)%mask(:,:,:),                             &
         Time%model_time, rmask = grid_tmask(:,:,:),          &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif  !}

enddo  !} n

!
!       Restore to the restore_region values
!

do n = 1, num_instances  !{

  if (instance(n)%some_restore) then  !{

    do k = 1, nk
      do j = jsd, jed  !{
        do i = isd, ied  !{
          ind = instance(n)%index(i,j,k)
          !if (instance(n)%region(ind)%mask(i,j,k) .eq. 1.0) then  !{}
          if (ind .gt. 0) then  !{
            if (instance(n)%region(ind)%restore .gt. 0.0) then  !{
              instance(n)%mask(i,j,k) =                                                 &
                   (instance(n)%region(ind)%restore_region_value -                      &
                    t_prog(instance(n)%tracer_index)%field(i,j,k,Time%tau)) /           &
                   (instance(n)%region(ind)%restore * sec_per_day) * grid_tmask(i,j,k)
              t_prog(instance(n)%tracer_index)%source(i,j,k) =                          &
                   t_prog(instance(n)%tracer_index)%source(i,j,k) + instance(n)%mask(i,j,k)
            endif  !}
          endif  !}
        enddo  !}i
      enddo  !} j
    enddo  !} k

  endif  !}

enddo  !} n

return

end subroutine ocean_residency_source 
! </SUBROUTINE> NAME="ocean_residency_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_start">
!
! <DESCRIPTION>
!       Start the ocean residency package
!
!       Residency surface area specification
!
!  west_bnd     : western longitude of residency region
!  east_bnd     : eastern longitude of residency region
!  south_bnd    : southern latitude of residency region
!  north_bnd    : northern latitude of residency region
!  top_bnd      : top depth of residency region
!  bottom_bnd   : bottom depth of residency region
!
!       To set the volumes, a number of namelists are read,
!       each containing the above values. You may specify up to
!       num_geog_region rectangular cubes bounded by
!       (west_bnd, east_bnd, north_bnd, south_bnd, top_bnd, bottom_bnd).
!       Any grid box whose center is in one of these volumes will
!       be considered to be part of the volume where the
!       residency is reset to zero every time-step.
!
!       north_bnd may not equal south_bnd, and west_bnd may not equal east_bnd
!
!       top_depth may equal bottom_depth. In that case, then whatever vertical
!       box contains that depth will define the vertical range for the box
!
!       If south_bnd > north_bnd, then nothing will be done for that rectangle
!
!       The initial surface area is empty, with the default rectangle
!       setting the surface area to be empty
!
!       More than num_geog_regions rectanglar volumes may be used to specify
!       the volume by using more than one namelist
! </DESCRIPTION>
!
subroutine ocean_residency_start(isd, ied, jsd, jed, nk, model_time, grid_tracer_axes)  !{

!
!       modules
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: nk
type(time_type), intent(in)                             :: model_time
integer, dimension(:), intent(in)                       :: grid_tracer_axes

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       local variables
!

integer                                 :: i
integer                                 :: j
integer                                 :: k
integer                                 :: n
integer                                 :: nn
character(len=256)                      :: caller_str
integer                                 :: len_w
integer                                 :: len_e
integer                                 :: len_s
integer                                 :: len_n
integer                                 :: len_t
integer                                 :: len_b
character(len=fm_field_name_len+3)      :: long_suffix
character(len=fm_field_name_len+1)      :: suffix
character(len=24)                       :: num_suffix
character(len=24)                       :: num_long_suffix
character(len=24)                       :: number_str
logical                                 :: found_error

  integer :: stdoutunit 
  stdoutunit=stdout() 

!
!-----------------------------------------------------------------------
!       save the instance namelist values
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, num_instances !{

  call fm_util_start_namelist(package_name, instance(n)%name, caller = caller_str)

  instance(n)%num_regions            =  fm_util_get_integer      ('num_regions', scalar = .true.)
  instance(n)%integrate_region_value =  fm_util_get_real         ('integrate_region_value', scalar = .true.)
  instance(n)%union                  =  fm_util_get_logical      ('union', scalar = .true.)
  instance(n)%int_module_name        =  fm_util_get_string       ('int_module_name', scalar = .true.)
  instance(n)%int_params             => fm_util_get_real_array   ('int_params')
  instance(n)%int_flags              => fm_util_get_logical_array('int_flags')
  instance(n)%int_strings            => fm_util_get_string_array ('int_strings')

!
!       create namelists for each region
!

  if (instance(n)%num_regions .le. 0) then
    call mpp_error(FATAL,trim(error_header) // ' num_regions is non-positive for instance ' //   &
         trim(instance(n)%name) // ' in package ' // trim(package_name))
  endif

!
!       allocate storage for this instance and
!       set all of the values to the default
!

  allocate( instance(n)%region(instance(n)%num_regions) )
  do nn = 1, instance(n)%num_regions  !{
    allocate( instance(n)%region(nn)%mask(isd:ied,jsd:jed,nk) )
    do k = 1, nk
      do j = jsd, jed
        do i = isd, ied
          instance(n)%region(nn)%mask(i,j,k) = 0.0
        enddo
      enddo
    enddo
  enddo  !} nn
  allocate( instance(n)%mask(isd:ied,jsd:jed,nk) )
  do k = 1, nk
    do j = jsd, jed
      do i = isd, ied
        instance(n)%mask(i,j,k) = 0.0
      enddo
    enddo
  enddo
  allocate( instance(n)%index(isd:ied,jsd:jed,nk) )
  do k = 1, nk
    do j = jsd, jed
      do i = isd, ied
        instance(n)%index(i,j,k) = 0
      enddo
    enddo
  enddo
  !instance(n)%mask => instance(n)%region(1)%mask

  instance(n)%some_restore = .false.
  instance(n)%some_fix = .false.
  do nn = 1, instance(n)%num_regions  !{

    if (instance(n)%num_regions .eq. 1) then
      suffix = ' '
    else
      write (long_suffix,*) nn
      suffix = '_' // long_suffix(scan(long_suffix,'0123456789'):)
    endif

    instance(n)%region(nn)%restore              =  fm_util_get_real         ('restore' // suffix, scalar = .true.)
    instance(n)%region(nn)%restore_region_value =  fm_util_get_real         ('restore_region_value' // suffix, scalar = .true.)
    instance(n)%region(nn)%swap                 =  fm_util_get_logical      ('swap' // suffix, scalar = .true.)
    instance(n)%region(nn)%swap_module          =  fm_util_get_logical      ('swap_module' // suffix, scalar = .true.)
    instance(n)%region(nn)%west_bnd             => fm_util_get_real_array   ('west_bnd' // suffix)
    instance(n)%region(nn)%east_bnd             => fm_util_get_real_array   ('east_bnd' // suffix)
    instance(n)%region(nn)%south_bnd            => fm_util_get_real_array   ('south_bnd' // suffix)
    instance(n)%region(nn)%north_bnd            => fm_util_get_real_array   ('north_bnd' // suffix)
    instance(n)%region(nn)%top_bnd              => fm_util_get_real_array   ('top_bnd' // suffix)
    instance(n)%region(nn)%bottom_bnd           => fm_util_get_real_array   ('bottom_bnd' // suffix)
    instance(n)%region(nn)%module_name          =  fm_util_get_string       ('module_name' // suffix, scalar = .true.)
    instance(n)%region(nn)%params               => fm_util_get_real_array   ('params' // suffix)
    instance(n)%region(nn)%flags                => fm_util_get_logical_array('flags' // suffix)
    instance(n)%region(nn)%strings              => fm_util_get_string_array ('strings' // suffix)

    instance(n)%some_restore = instance(n)%some_restore .or. instance(n)%region(nn)%restore .gt. 0.0
    instance(n)%some_fix = instance(n)%some_fix .or. instance(n)%region(nn)%restore .eq. 0.0

  enddo  !} nn

  call fm_util_end_namelist(package_name, instance(n)%name, caller = caller_str)

!
!       Check some things
!

  do nn = 1, instance(n)%num_regions  !{

    if (.not. associated(instance(n)%region(nn)%west_bnd) .and.      &
        .not. associated(instance(n)%region(nn)%east_bnd) .and.      &
        .not. associated(instance(n)%region(nn)%south_bnd) .and.     &
        .not. associated(instance(n)%region(nn)%north_bnd) .and.     &
        .not. associated(instance(n)%region(nn)%top_bnd) .and.       &
        .not. associated(instance(n)%region(nn)%bottom_bnd)) then  !{

      write (stdoutunit,*) trim(note_header),               &
           ' No region specified, assuming it will be specified externally: ', trim(instance(n)%name), ' region: ', nn
      instance(n)%region(nn)%num_geog_regions = 0

    elseif (.not. associated(instance(n)%region(nn)%west_bnd) .or.   &
        .not. associated(instance(n)%region(nn)%east_bnd) .or.       &
        .not. associated(instance(n)%region(nn)%south_bnd) .or.      &
        .not. associated(instance(n)%region(nn)%north_bnd) .or.      &
        .not. associated(instance(n)%region(nn)%top_bnd) .or.        &
        .not. associated(instance(n)%region(nn)%bottom_bnd)) then  !}{

      call mpp_error(FATAL, trim(error_header) // ' Some regions not specified: ' // trim(instance(n)%name))

    else  !}{

      len_w = size(instance(n)%region(nn)%west_bnd)
      len_e = size(instance(n)%region(nn)%east_bnd)
      len_s = size(instance(n)%region(nn)%south_bnd)
      len_n = size(instance(n)%region(nn)%north_bnd)
      len_t = size(instance(n)%region(nn)%top_bnd)
      len_b = size(instance(n)%region(nn)%bottom_bnd)

      if (len_e .ne. len_w .or. len_e .ne. len_s .or. len_e .ne. len_n .or.         &
          len_e .ne. len_t .or. len_e .ne. len_b) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal: ' // trim(instance(n)%name))
      endif  !}
      instance(n)%region(nn)%num_geog_regions = len_w

    endif  !}

  enddo  !} nn

!
!       register the fields
!

  if (instance(n)%name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // instance(n)%name
    long_suffix = ' (' // trim(instance(n)%name) // ')'
  endif  !}

  instance(n)%id_change = register_diag_field(trim(diag_name),                  &
       'residency_change' // trim(suffix), grid_tracer_axes(1:3),               &
       model_time, 'Residency change' // trim(long_suffix), ' ',                &
       missing_value = -1.0e+10)

  instance(n)%id_restore_region = register_diag_field(trim(diag_name),          &
       'residency_in_merged_regions' // trim(suffix), grid_tracer_axes(1:3),    &
       model_time, 'Residency in merged regions' // trim(long_suffix), ' ',     &
       missing_value = -1.0e+10)

  do nn = 1, instance(n)%num_regions  !{

    if (instance(n)%num_regions .eq. 1) then
      num_suffix = ' '
    else
      write (number_str,*) nn
      num_suffix = '_' // number_str(scan(number_str,'0123456789'):)
      num_long_suffix = ' ' // number_str(scan(number_str,'0123456789'):)
    endif

    instance(n)%region(nn)%id_restore_region = register_diag_field(trim(diag_name),                     &
         'residency_restore_region' // trim(num_suffix) // trim(suffix), grid_tracer_axes(1:3),         &
         model_time, 'Residency in region' // trim(num_long_suffix) // trim(long_suffix), ' ',          &
         missing_value = -1.0e+10)

  enddo  !} nn

enddo  !} n

!
!       call the start routines for instances controlled by external modules
!

call ocean_residency_ml_start

call ocean_residency_range_start

call ocean_residency_integrand_start

!
!       check that the external tracers are all accounted for
!

found_error = .false.
do n = 1, num_instances  !{
  if (instance(n)%int_module_name .ne. ' ' .and. .not. instance(n)%int_found) then  !}{
    found_error = .true.
    write (stdoutunit,*) trim(error_header), ' Instance "', trim(instance(n)%name),      &
         '" not found with integrand module name "', trim(instance(n)%int_module_name), '"'
  endif  !}
  do nn = 1, instance(n)%num_regions  !{
    if (instance(n)%region(nn)%module_name .ne. ' ' .and. .not. instance(n)%region(nn)%found) then  !{
      found_error = .true.
      write (stdoutunit,*) trim(error_header), ' Instance "', trim(instance(n)%name),      &
           '" not found with module name "', trim(instance(n)%region(nn)%module_name), '"'
    endif  !}
  enddo  !} nn
enddo  !} n
if (found_error) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Some external instances not found')
endif  !}

return

end subroutine ocean_residency_start  !}
! </SUBROUTINE> NAME="ocean_residency_start"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_tracer">
!
! <DESCRIPTION>
!       Subroutine to do calculations needed every time-step after
!       the continuity equation has been integrated
! </DESCRIPTION>
!

subroutine ocean_residency_tracer(isc, iec, jsc, jec,           &
     isd, ied, jsd, jed, nk, T_prog, grid_tmask, taup1, model_time, dtts) !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: taup1
type(time_type), intent(in)                                     :: model_time
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
real, intent(in)                                                :: dtts

!
!       local parameters
!

!
!       local variables
!

integer :: i
integer :: j
integer :: k
integer :: n
integer :: ind
logical :: used

!
!       fix the restore_region values
!

do n = 1, num_instances  !{

  if (instance(n)%some_fix) then  !{

    do k = 1, nk
      do j = jsd, jed  !{
        do i = isd, ied  !{
          ind = instance(n)%index(i,j,k)
          !if (instance(n)%region(ind)%mask(i,j,k) .eq. 1.0) then  !{}
          if (ind .gt. 0) then  !{
            if (instance(n)%region(ind)%restore .eq. 0.0) then  !{
              instance(n)%mask(i,j,k) =                                                 &
                   (instance(n)%region(ind)%restore_region_value -                      &
                    t_prog(instance(n)%tracer_index)%field(i,j,k,taup1)) /              &
                   dtts * grid_tmask(i,j,k)
              t_prog(instance(n)%tracer_index)%field(i,j,k,taup1) = instance(n)%region(ind)%restore_region_value
            endif  !}
          endif  !}
        enddo  !}i
      enddo  !} j
    enddo  !} k

  endif  !}

  if (instance(n)%id_change .gt. 0) then
    used = send_data(instance(n)%id_change,                   &
         instance(n)%mask(:,:,:),                             &
         model_time, rmask = grid_tmask(:,:,:),               &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

enddo  !} n

return

end subroutine ocean_residency_tracer  !}
! </SUBROUTINE> NAME="ocean_residency_tracer"

end module ocean_residency_mod  !}
