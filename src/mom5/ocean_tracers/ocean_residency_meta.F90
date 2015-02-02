module ocean_residency_meta_mod  !{
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
!       This module contains the meta definitions, subroutines and functions
!       to use in the ocean residency package. These routines are used by all
!       of the ocean residency modules and need to be defined in a separate
!       module so that there are no circular references between the modules.
!
!       For an overview of the ocean residency modules, see ocean_residency.F90.
!
!</DESCRIPTION>
!
! $Id: ocean_residency_meta.F90,v 20.0 2013/12/14 00:17:08 fms Exp $
!

!
!       modules
!

use field_manager_mod, only: fm_string_len, fm_field_name_len
use fm_util_mod,       only: fm_util_default_caller
use mpp_mod,           only: stdout, mpp_error, FATAL

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

public  :: ocean_residency_get_instances
public  :: ocean_residency_set_region_geog
public  :: ocean_residency_set_region_2d
public  :: ocean_residency_set_region_3d

!
!       Private routines
!

!
!       Public parameters
!

real, public, parameter                         :: secs_in_year_r = 1.0 / (86400.0 * 365.25)
real, public, parameter                         :: sec_per_day = 86400.0

!
!       Private parameters
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_residency_meta'
character(len=48), parameter                    :: mod_name = 'ocean_residency_meta_mod'
character(len=48), parameter                    :: diag_name = 'ocean_residency_meta'

!
!       Public types
!

type, public                                                    :: ocean_residency_instance_type
  character(len=fm_field_name_len)                              :: name = ' '
  integer                                                       :: tracer_index = 0
  logical                                                       :: some_restore = .false.
  logical                                                       :: some_fix = .false.
  integer                                                       :: num_regions = 1
  type(ocean_residency_region_type), dimension(:), pointer      :: region => NULL()
  real                                                          :: integrate_region_value = secs_in_year_r
  integer, dimension(:,:,:), pointer                            :: index => NULL()
  real, dimension(:,:,:), pointer                               :: mask => NULL()
  logical                                                       :: union = .true.
  integer                                                       :: id_restore_region = -1
  integer                                                       :: id_change = -1
  logical                                                       :: int_found = .false.
  character(len=128)                                            :: int_module_name = ' '
  real, dimension(:), pointer                                   :: int_params => NULL()
  logical, dimension(:), pointer                                :: int_flags => NULL()
  character(len=fm_string_len), dimension(:), pointer           :: int_strings => NULL()
end type ocean_residency_instance_type

type, public                                                    :: ocean_residency_region_type
  character(len=fm_field_name_len)                              :: name = ' '
  real, dimension(:,:,:), pointer                               :: mask => NULL()
  real                                                          :: restore = 0.0
  real                                                          :: restore_region_value = 0.0
  logical                                                       :: swap = .false.
  logical                                                       :: swap_module = .false.
  logical                                                       :: found = .false.
  real, dimension(:), pointer                                   :: east_bnd => NULL()
  real, dimension(:), pointer                                   :: north_bnd => NULL()
  real, dimension(:), pointer                                   :: south_bnd => NULL()
  real, dimension(:), pointer                                   :: west_bnd => NULL()
  real, dimension(:), pointer                                   :: top_bnd => NULL()
  real, dimension(:), pointer                                   :: bottom_bnd => NULL()
  integer                                                       :: num_geog_regions = 0
  integer                                                       :: id_restore_region = -1
  character(len=128)                                            :: module_name = ' '
  real, dimension(:), pointer                                   :: params => NULL()
  logical, dimension(:), pointer                                :: flags => NULL()
  character(len=fm_string_len), dimension(:), pointer           :: strings => NULL()
end type ocean_residency_region_type

!
!       Public variables
!

type(ocean_residency_instance_type), dimension(:), allocatable, save, public    :: instance
integer, save, public                                                           :: num_instances = 0

!
!       Private variables
!


!
!        Interface definitions for overloaded routines
!

contains


!#######################################################################
! <FUNCTION NAME="ocean_residency_get_instances">
!
! <DESCRIPTION>
!       Return an array of instances which have the given module_name
!       This is used by modules, such as the mixed layer module, to obtain
!       a list of instances which use that module. Then, the module
!       needs only to loop over those instance to perform its required tasks.
! </DESCRIPTION>
!

function ocean_residency_get_instances(module_name, integrand)            &
     result(array)  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

character(len=*), intent(in)    :: module_name
logical, intent(in), optional   :: integrand

!
!-----------------------------------------------------------------------
!       Result
!-----------------------------------------------------------------------
!

type(ocean_residency_instance_type), dimension(:), pointer      :: array

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_get_instances'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

logical                                 :: found
integer                                 :: n
integer                                 :: ind
integer                                 :: nn
integer                                 :: num_elements
logical                                 :: do_integrand

if (module_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' No module name given')
endif  !}

!
!       check whether we want to look at modules for defining the region
!       (integrand = f or not given) or whether we are looking to set
!       the field we're integrating over to something other than time
!       (integrand = t)
!

if (present(integrand)) then  !{
  do_integrand = integrand
else  !}{
  do_integrand = .false.
endif  !}

!
!       return the mask array pointer for tracer "tracer_index"
!

num_elements = 0
do n = 1, num_instances  !{
  if (do_integrand) then  !{
    if (module_name .eq. instance(n)%int_module_name) then
      num_elements = num_elements + 1
    endif
  else  !}{
    do nn = 1, instance(n)%num_regions  !{
      if (module_name .eq. instance(n)%region(nn)%module_name) then
        num_elements = num_elements + 1
        exit
      endif
    enddo  !} nn
  endif  !}
enddo  !} n

if (num_elements .eq. 0) then  !{

  nullify(array)

elseif (num_elements .lt. 0) then  !}{

  call mpp_error(FATAL, trim(error_header) // ' Should never get here for ' // trim(module_name))

else  !}{

  allocate( array(num_elements) )
  ind = 0
  do n = 1, num_instances  !{
    found = .false.
    if (do_integrand) then  !{
      if (module_name .eq. instance(n)%int_module_name) then
        found = .true.
        instance(n)%int_found = .true.
      endif
    else  !}{
      do nn = 1, instance(n)%num_regions  !{
        if (module_name .eq. instance(n)%region(nn)%module_name) then
          found = .true.
          instance(n)%region(nn)%found = .true.
        endif
      enddo  !} nn
    endif  !}
    if (found) then  !{
      ind = ind + 1
      allocate ( array(ind)%region(instance(n)%num_regions) )
      array(ind)%num_regions            =  instance(n)%num_regions
      array(ind)%region                 => instance(n)%region
      array(ind)%name                   =  instance(n)%name
      array(ind)%tracer_index           =  instance(n)%tracer_index
      array(ind)%integrate_region_value =  instance(n)%integrate_region_value
      array(ind)%index                  => instance(n)%index
      array(ind)%mask                   => instance(n)%mask
      array(ind)%union                  =  instance(n)%union
      array(ind)%id_change              =  instance(n)%id_change
      array(ind)%int_module_name        =  instance(n)%int_module_name
      array(ind)%int_params             => instance(n)%int_params
      array(ind)%int_flags              => instance(n)%int_flags
      array(ind)%int_strings            => instance(n)%int_strings
    endif  !}
  enddo  !} n

  if (num_elements .ne. ind) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Number of elements do not match for ' // trim(module_name))
  endif  !}

endif  !}

return

end function ocean_residency_get_instances  !}
! </FUNCTION> NAME="ocean_residency_get_instances"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_set_region_2d">
!
! <DESCRIPTION>
!       Given a 2-d field of depths, determine the grid cells which fall inside (above) and
!       outside of that range, then set the mask appropriately, as to whether
!       we are interested in those points inside (swap=false) or outside (swap=true) of the region.
!
!       Note: if the output array is not initialized, then it is assumed that it is already
!             filled with a region, and the resulting region will be the intersection of the
!             existing region and the newly specified region.
!
!       Arguments:
!
!         Input:
!
!                    isd:  low dimension of first index
!                    ied:  high dimension of first index
!                    jsd:  low dimension of second index
!                    jed:  high dimension of second index
!                     nk:  dimension of third index
!                control:  array of depths that specifies the region
!              depth_zwt:  depth of bottom of grid cell
!   restore_region_value:  if supplied, the value to assign to array in the user-specified
!                          regions (default: 0.0)
! integrate_region_value:  if supplied, the value to assign to array outside of the
!                          user-specified regions (default: secs_in_year_r)
!                   swap:  if supplied and true then change the invert the defined region (default: true)
!             initialize:  if supplied and true, then initialize the region to the integrate_region_value (default: false)
!                 caller:  if supplied, use for traceback of error messages (default: fm_default_caller)
!
!         Input/Output:
!                  array:  3-d array which will contain the restore_region- and integrate_region-
!                          values.
!
! </DESCRIPTION>
!

subroutine ocean_residency_set_region_2d(isd, ied, jsd, jed, nk, control, array, depth_zwt,     &
     restore_region_value, integrate_region_value, swap, initialize, caller)  !{

implicit none

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'ocean_residency_set_region_2d'

!
!       arguments
!

integer, intent(in)                             :: isd
integer, intent(in)                             :: ied
integer, intent(in)                             :: jsd
integer, intent(in)                             :: jed
integer, intent(in)                             :: nk
real, intent(in), dimension(isd:,jsd:)          :: control
real, dimension(isd:,jsd:,:), intent(inout)     :: array
real, dimension(isd:,jsd:,:), intent(in)        :: depth_zwt
real, intent(in), optional                      :: restore_region_value
real, intent(in), optional                      :: integrate_region_value
logical, intent(in), optional                   :: swap
logical, intent(in), optional                   :: initialize
character(len=*), intent(in), optional          :: caller

!
!       Local variables
!

character(len=256)      :: error_header
character(len=128)      :: caller_str
integer                 :: i
integer                 :: j
integer                 :: k
real                    :: restore_region_value_use
real                    :: integrate_region_value_use
logical                 :: swap_use
logical                 :: assume_restore_region

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
  
!
! check that the array is as long as depth dimension
!

if (size(array,3) .lt. nk) then  
  call mpp_error(FATAL, trim(error_header) // ' Dimension 3 of output array too small')
endif

!
!       set some default values in case they have not been set in the argument list
!

if (present(restore_region_value)) then  !{
  restore_region_value_use = restore_region_value
else  !}{
  restore_region_value_use = 0.0
endif  !}

if (present(integrate_region_value)) then  !{
  integrate_region_value_use = integrate_region_value
else  !}{
  integrate_region_value_use = secs_in_year_r
endif  !}

if (present(swap)) then  !{
  swap_use = swap
else  !}{
  swap_use = .false.
endif  !}

!
!       If initialize is true, then set all values to integrate_region_value,
!       otherwise, array values are assumed to be initialized with the appropriate
!       distribution of in- and out-region-values and the final value will be
!       the intersection between this value and that determined below (i.e., both values
!       must be "restore_region_value" to stay "restore_region_value", otherwise it becomes "integrate_region_value")
!       If we initialize, we need to set assume_restore_region to true so that when
!       calculating the region below, the intersection will have values in it.
!       However, we don't want to initialize to restore_region_value, or else the behavior
!       when no regions are found will be to have all locations restore_region.
!

if (present(initialize)) then  !{
  assume_restore_region = initialize
  if (initialize) then  !{
    do k = 1, nk  !{
      do j = jsd, jed  !{
        do i = isd, ied  !{
          array(i,j,k) = integrate_region_value_use
        enddo  !} i
      enddo  !} j
    enddo  !} k
  endif  !}
else  !}{
  assume_restore_region = .false.
endif  !}

do j = jsd, jed  !{
  do i = isd, ied  !{
    if (((control(i,j) .gt. 0.0) .eqv. (.not. swap_use)) .and.   &
        (assume_restore_region .or. array(i,j,1) .eq. restore_region_value_use)) then  !{
      array(i,j,1) = restore_region_value_use
    else  !}{
      array(i,j,1) = integrate_region_value_use
    endif  !}
  enddo  !} i
enddo  !} j
do k = 2, nk  !{
  do j = jsd, jed  !{
    do i = isd, ied  !{
      if (((control(i,j) .gt. depth_zwt(i,j,k-1)) .eqv. (.not. swap_use)) .and.   &
          (assume_restore_region .or. array(i,j,k) .eq. restore_region_value_use)) then  !{
        array(i,j,k) = restore_region_value_use
      else  !}{
        array(i,j,k) = integrate_region_value_use
      endif  !}
    enddo  !} i
  enddo  !} j
enddo  !} k

return

end subroutine  ocean_residency_set_region_2d  !}
! </SUBROUTINE> NAME="ocean_residency_set_region_2d"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_set_region_geog">
!
! <DESCRIPTION>
!       Set up an array covering the model domain with a user-specified
!       value, in user-specified regions. There are a given number of
!       3-d regions specified by the values south_bnd, north_bnd, west_bnd,
!       east_bnd, top_bnd and bottom_bnd.
!       The longitudes are for a cyclic domain, and if west_bnd and east_bnd
!       are on opposite sides of the cut, the correct thing will
!       be done. east_bnd is considered to be east of west_bnd, so if east_bnd is
!       less than west_bnd, then the region east of east_bnd to the cut will be
!       filled, and the region from the cut to west_bnd will be filled.
!
!       If the bottom bound is less than or equal to zero, then the top model box
!       will be chosen. If the top bound is greater than or equal to the maximum
!       model depth, then the bottom box will be chosed. Otherwise, if the grid
!       cell center depth falls between top bound and bottom bound, then those cells
!       shall be chosen.
!
!       For longitude and latitude, if the grid cell center lies within the
!       rectabgle defined by (west_bnd,south_bnd) and (east_bnd,north_bnd), then
!       the whole grid cell is inside the region.
!
!       Arrays of coordinates may be specified for irregular regions.
!       The final region is the union of the multiple sets
!       of coordinates. If swap is true, then the inverse of the defined
!       region will be set.
!
!       Note: if the output array is not initialized, then it is assumed that it is already
!             filled with a region, and the resulting region will be the intersection of the
!             existing region and the newly specified region.
!
!       Arguments:
!
!         Input:
!
!                    isd:  low dimension of first index
!                    ied:  high dimension of first index
!                    jsd:  low dimension of second index
!                    jed:  high dimension of second index
!                     nk:  dimension of third index
!                grid_xt:  array of coordinates in the x-direction (typically longitude)
!                grid_yt:  array of coordinates in the y-direction (typically latitude)
!              max_depth:  maximum depth of the model
!               depth_zt:  depth of center of grid cell
!              depth_zwt:  depth of bottom of grid cell
!       num_geog_regions:  number of user-specified regions which will be filled
!            west_bnd_in:  1-d array of western (starting) longitudes for the regions
!            east_bnd_in:  1-d array of eastern (ending) longitudes for the regions
!              south_bnd:  1-d array of southern (starting) latitudes for the regions
!              north_bnd:  1-d array of northern (ending) latitudes for the regions
!                top_bnd:  1-d array of southern (starting) depths for the regions
!             bottom_bnd:  1-d array of northern (ending) depths for the regions
!                             Note: if south_bnd >= north_bnd, then nothing is done
!                                   for that region
!                    kmt:  array of indices for bottom grid cells
!                   name:  character variable used in informative messages
!   restore_region_value:  if supplied, the value to assign to array in the user-specified
!                          regions (default: 0.0)
! integrate_region_value:  if supplied, the value to assign to array outside of the
!                          user-specified regions (default: secs_in_year_r)
!                   swap:  if supplied and true then change the invert the defined region (default: true)
!
!         Input/Output:
!                  array:  3-d array which will contain the restore_region- and integrate_region-
!                          values.
! </DESCRIPTION>
!
!

subroutine ocean_residency_set_region_geog(isd, ied, jsd, jed, nk, array,       &
                     grid_xt, grid_yt, max_depth, depth_zt, depth_zwt,          &
                     num_geog_regions,                                          &
                     west_bnd_in, east_bnd_in, south_bnd, north_bnd,            &
                     top_bnd, bottom_bnd, kmt, name,                            &
                     restore_region_value, integrate_region_value, swap,        &
                     initialize, caller)  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                             :: isd
integer, intent(in)                             :: ied
integer, intent(in)                             :: jsd
integer, intent(in)                             :: jed
integer, intent(in)                             :: nk
real, dimension(isd:,jsd:,:), intent(inout)     :: array
real, dimension(isd:,jsd:), intent(in)          :: grid_xt
real, dimension(isd:,jsd:), intent(in)          :: grid_yt
real, dimension(isd:,jsd:,:), intent(in)        :: depth_zt
real, dimension(isd:,jsd:,:), intent(in)        :: depth_zwt
real, intent(in)                                :: max_depth
integer, intent(in)                             :: num_geog_regions
real, dimension(:), intent(in)                  :: west_bnd_in
real, dimension(:), intent(in)                  :: east_bnd_in
real, dimension(:), intent(in)                  :: south_bnd
real, dimension(:), intent(in)                  :: north_bnd
real, dimension(:), intent(in)                  :: top_bnd
real, dimension(:), intent(in)                  :: bottom_bnd
integer, dimension(isd:,jsd:), intent(in)       :: kmt
real, intent(in), optional                      :: restore_region_value
real, intent(in), optional                      :: integrate_region_value
logical, intent(in), optional                   :: swap
character*(*), intent(in)                       :: name
logical, intent(in), optional                   :: initialize
character(len=*), intent(in), optional          :: caller

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_set_region_geog'

!
!       local variables
!

integer                         :: i
integer                         :: j
integer                         :: k
integer                         :: n
real, dimension(:), allocatable :: west_bnd
real, dimension(:), allocatable :: east_bnd
character(len=256)              :: error_header
character(len=128)              :: caller_str
real                            :: restore_region_value_use
real                            :: integrate_region_value_use
logical                         :: swap_use
logical                         :: restore_region
logical                         :: assume_restore_region

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       set some default values in case they have not been set in the argument list
!

if (present(restore_region_value)) then  !{
  restore_region_value_use = restore_region_value
else  !}{
  restore_region_value_use = 0.0
endif  !}

if (present(integrate_region_value)) then  !{
  integrate_region_value_use = integrate_region_value
else  !}{
  integrate_region_value_use = secs_in_year_r
endif  !}

if (present(swap)) then  !{
  swap_use = swap
else  !}{
  swap_use = .false.
endif  !}

!
!       check whether any regions have been specified
!

if (num_geog_regions .le. 0) then  !{

!
!       no regions specified, return all points out of region
!       unless swap is true
!

  if (swap_use) then  !{

    do k = 1, nk  !{
      do j = jsd, jed  !{
        do i = isd, ied  !{
          array(i,j,k) = restore_region_value_use
        enddo  !} i
      enddo  !} j
    enddo  !} k
  
  else  !}{

    do k = 1, nk  !{
      do j = jsd, jed  !{
        do i = isd, ied  !{
          array(i,j,k) = integrate_region_value_use
        enddo  !} i
      enddo  !} j
    enddo  !} k
  
  endif  !}

else  !}{

!
!       regions are specified, check that all array large enough
!

  if (size(west_bnd_in) .lt. num_geog_regions) then  !{
    call mpp_error(FATAL, trim(error_header) // ' west_bnd_in array too small for ' // trim(name))
  elseif (size(east_bnd_in) .lt. num_geog_regions) then  !}{
    call mpp_error(FATAL, trim(error_header) // ' east_bnd_in array too small for ' // trim(name))
  elseif (size(north_bnd) .lt. num_geog_regions) then  !}{
    call mpp_error(FATAL, trim(error_header) // ' north_bnd array too small for ' // trim(name))
  elseif (size(south_bnd) .lt. num_geog_regions) then  !}{
    call mpp_error(FATAL, trim(error_header) // ' south_bnd array too small for ' // trim(name))
  elseif (size(top_bnd) .lt. num_geog_regions) then  !}{
    call mpp_error(FATAL, trim(error_header) // ' top_bnd array too small for ' // trim(name))
  elseif (size(bottom_bnd) .lt. num_geog_regions) then  !}{
    call mpp_error(FATAL, trim(error_header) // ' bottom_bnd array too small for ' // trim(name))
  endif  !}

!
! loop over the regions, checking for errors
!

  do n = 1, num_geog_regions  !{

    if (north_bnd(n) .lt. south_bnd(n)) then  !{

      call mpp_error(FATAL, trim(error_header) //       &
           ' North bound less than south bound for ' // trim(name))

    elseif (top_bnd(n) .gt. bottom_bnd(n)) then  !}{

      call mpp_error(FATAL, trim(error_header) //       &
           ' Top depth deeper than bottom depth for ' // trim(name))

    endif  !}

  enddo  !} n

!
!       Find all boxes where the center lies in the
!       3-d rectangular region
!

!
!       If initialize is true, then set all values to integrate_region_value,
!       otherwise, array values are assumed to be initialized with the appropriate
!       distribution of in- and out-region-values and the final value will be
!       the intersection between this value and that determined below (i.e., both values
!       must be "restore_region_value" to stay "restore_region_value", otherwise it becomes "integrate_region_value")
!       If we initialize, we need to set assume_restore_region to true so that when
!       calculating the region below, the intersection will have values in it.
!       However, we don't want to initialize to restore_region_value, or else the behavior
!       when no regions are found will be to have all locations restore_region.
!

  if (present(initialize)) then  !{
    assume_restore_region = initialize
    if (initialize) then  !{
      do k = 1, nk  !{
        do j = jsd, jed  !{
          do i = isd, ied  !{
            array(i,j,k) = integrate_region_value_use
          enddo  !} i
        enddo  !} j
      enddo  !} k
    endif  !}
  else  !}{
    assume_restore_region = .false.
  endif  !}

!
!       save the longitudes in case they need to be modified
!

  allocate( west_bnd(num_geog_regions) )
  allocate( east_bnd(num_geog_regions) )

  do n = 1, num_geog_regions  !{

    west_bnd(n) = west_bnd_in(n)
    east_bnd(n) = east_bnd_in(n)

!
!       make sure that the longitudes are in the range [0,360]
!

    do while (west_bnd(n) .gt. 360.0)  !{
      west_bnd(n) = west_bnd(n) - 360.0
    enddo  !}
    do while (west_bnd(n) .lt. 0.0)  !{
      west_bnd(n) = west_bnd(n) + 360.0
    enddo  !}
    do while (east_bnd(n) .gt. 360.0)  !{
      east_bnd(n) = east_bnd(n) - 360.0
    enddo  !}
    do while (east_bnd(n) .lt. 0.0)  !{
      east_bnd(n) = east_bnd(n) + 360.0
    enddo  !}

  enddo  !} n

  do k = 1, nk  !{
    do j = jsd, jed  !{
      do i = isd, ied  !{
        if (k .le. kmt(i,j)) then  !{
          restore_region = .false.
          do n = 1, num_geog_regions  !{
            if ((((top_bnd(n) .lt. depth_zt(i,j,k) .and.                        &
                   depth_zt(i,j,k) .le. bottom_bnd(n)) .or.                     &
                  (bottom_bnd(n) .le. 0.0 .and. k .eq. 1) .or.                  &
                  (top_bnd(n) .ge. max_depth .and. k .eq. kmt(i,j)) .or.        &
                  (top_bnd(n) .lt. 0.0 .and.                                    &
                   depth_zt(i,j,k) .ge.                                         &
                   depth_zwt(i,j,kmt(i,j)) + top_bnd(n))) .and.                 &
                 north_bnd(n) .gt. grid_yt(i,j) .and.                           &
                 grid_yt(i,j) .ge. south_bnd(n) .and.                           &
                 lon_between(grid_xt(i,j), west_bnd(n), east_bnd(n)))) then  !{
              restore_region = .true.
            endif  !}
          enddo  !} n
          restore_region = restore_region .eqv. (.not. swap_use)
          if ((assume_restore_region .or. array(i,j,k) .eq. restore_region_value_use) .and. restore_region) then  !{
            array(i,j,k) = restore_region_value_use
          else  !}{
            array(i,j,k) = integrate_region_value_use
          endif  !}
        endif  !}
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!       Clean up
!

  deallocate( west_bnd )
  deallocate( east_bnd )

endif  !}

return

contains

!
!       Return true if w <= x_in < e, taking into account the
!       periodicity of longitude.
!
!       x_in    value to test
!       w       west longitude of boundary
!       e       east longitude of boundary
!

function lon_between(x_in, w, e)  !{

implicit none

!
!       function definition
!

logical :: lon_between

!
!       arguments
!

real, intent(in)                :: x_in
real, intent(in)                :: w
real, intent(in)                :: e

!
!       local variables
!

real                    :: x

!
!       Save input values so we may modify them safely
!

x = x_in

!
!       make sure that all longitudes are in the range [0,360]
!

do while (x .gt. 360.0)  !{
  x = x - 360.0
enddo  !}
do while (x .lt. 0.0)  !{
  x = x + 360.0
enddo  !}
 
if (w .gt. e) then  !{
  lon_between = w .le. x .or. x .lt. e
else  !}{
  lon_between = w .le. x .and. x .lt. e
endif  !}

return

end function  lon_between  !}

end subroutine  ocean_residency_set_region_geog  !}
! </SUBROUTINE> NAME="ocean_residency_set_region_geog"


!#######################################################################
! <SUBROUTINE NAME="ocean_residency_set_region_3d">
!
! <DESCRIPTION>
!       Set up an array where the a grid box is in the region if the value
!       of the specified property (temperature, say) is within the given bounds.
!       Multiple values for the range may be given, and the resulting mask
!       will be the union of the multiple regions. If swap is true, then the
!       inverse of the selected region will be set.
!
!       Note: if the output array is not initialized, then it is assumed that it is already
!             filled with a region, and the resulting region will be the intersection of the
!             existing region and the newly specified region.
!
!       Arguments:
!
!         Input:
!
!                    isd:  low dimension of first index
!                    ied:  high dimension of first index
!                    jsd:  low dimension of second index
!                    jed:  high dimension of second index
!                     nk:  dimension of third index
!       num_geog_regions:  number of user-specified regions which will be filled
!                 bounds:  1-d array of pairs of bounding values. The first value in
!                          the pair must be less than the second value
!                    kmt:  array of indices for bottom grid cells
!                   name:  character variable used in informative messages
!   restore_region_value:  if supplied, the value to assign to array in the user-specified
!                          regions (default: 0.0)
! integrate_region_value:  if supplied, the value to assign to array outside of the
!                          user-specified regions (default: secs_in_year_r)
!                   swap:  if supplied and true then change the invert the defined region (default: true)
!             initialize:  if supplied and true, then initialize the region to the integrate_region_value (default: false)
!                 caller:  if supplied, use for traceback of error messages (default: fm_default_caller)
!
!         Input/Output:
!                  array:  3-d array which will contain the restore_region- and integrate_region-
!                          values.
! </DESCRIPTION>
!
!

subroutine ocean_residency_set_region_3d(isd, ied, jsd, jed, nk, array,         &
                     variable, bounds, kmt,                                &
                     restore_region_value, integrate_region_value, swap,        &
                     initialize, caller)  !{

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
real, dimension(isd:,jsd:,:), intent(inout)             :: array
real, dimension(isd:,jsd:,:), intent(in)                :: variable
real, dimension(:), intent(in)                          :: bounds
integer, dimension(isd:,jsd:), intent(in)               :: kmt
real, intent(in), optional                              :: restore_region_value
real, intent(in), optional                              :: integrate_region_value
logical, intent(in), optional                           :: swap
logical, intent(in), optional                           :: initialize
character(len=*), intent(in), optional                  :: caller

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_residency_set_region_3d'

!
!       local variables
!

character(len=132)      :: error_string
integer                 :: i
integer                 :: j
integer                 :: k
integer                 :: n
logical                 :: restore_region
character(len=256)      :: error_header
character(len=128)      :: caller_str
real                    :: restore_region_value_use
real                    :: integrate_region_value_use
logical                 :: swap_use
logical                 :: assume_restore_region

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = fm_util_default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that an even number of bounds have been given
!

if (mod(size(bounds),2) .ne. 0) then
  call mpp_error(FATAL, trim(error_header) // ' Odd number of bounds given')
endif

!
! check that the bounds are correctly set (i.e., lesser value first)
!

do n = 1, size(bounds), 2  !{
  if (bounds(n+1) .le. bounds(n)) then  !{
    write (error_string,*) ' Bounds out of order: ', bounds(n+1), ' <= ', bounds(n)
    call mpp_error(FATAL, trim(error_header) // trim(error_string))
  endif  !}
enddo  !} n

!
!       set some default values in case they have not been set in the argument list
!

if (present(restore_region_value)) then  !{
  restore_region_value_use = restore_region_value
else  !}{
  restore_region_value_use = 0.0
endif  !}

if (present(integrate_region_value)) then  !{
  integrate_region_value_use = integrate_region_value
else  !}{
  integrate_region_value_use = secs_in_year_r
endif  !}

if (present(swap)) then  !{
  swap_use = swap
else  !}{
  swap_use = .false.
endif  !}

!
!       If initialize is true, then set all values to integrate_region_value,
!       otherwise, array values are assumed to be initialized with the appropriate
!       distribution of in- and out-region-values and the final value will be
!       the intersection between this value and that determined below (i.e., both values
!       must be "restore_region_value" to stay "restore_region_value", otherwise it becomes "integrate_region_value")
!       If we initialize, we need to set assume_restore_region to true so that when
!       calculating the region below, the intersection will have values in it.
!       However, we don't want to initialize to restore_region_value, or else the behavior
!       when no regions are found will be to have all locations restore_region.
!

if (present(initialize)) then  !{
  assume_restore_region = initialize
  if (initialize) then  !{
    do k = 1, nk  !{
      do j = jsd, jed  !{
        do i = isd, ied  !{
          array(i,j,k) = integrate_region_value_use
        enddo  !} i
      enddo  !} j
    enddo  !} k
  endif  !}
else  !}{
  assume_restore_region = .false.
endif  !}

!
!       find all boxes where the center lies in the region
!

do k = 1, nk  !{
  do j = jsd, jed  !{
    do i = isd, ied  !{
      if (k .le. kmt(i,j)) then  !{
        restore_region = .false.
        do n = 1, size(bounds), 2  !{
          if ((bounds(n+1) .ge. variable(i,j,k) .and. variable(i,j,k) .gt. bounds(n))) then  !{
            restore_region = .true.
            exit
          endif  !}
        enddo  !} n
        restore_region = restore_region .eqv. (.not. swap_use)
        if ((assume_restore_region .or. array(i,j,k) .eq. restore_region_value_use) .and. restore_region) then  !{
          array(i,j,k) = restore_region_value_use
        else  !}{
          array(i,j,k) = integrate_region_value_use
        endif  !}
      endif  !}
    enddo  !} i
  enddo  !} j
enddo  !} k

return

end subroutine  ocean_residency_set_region_3d  !}
! </SUBROUTINE> NAME="ocean_residency_set_region_3d"

end module ocean_residency_meta_mod  !}
