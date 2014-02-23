#include <fms_platform.h>

! ----------------------------------------------------------------
!                   GNU General Public License                        
! This file is a part of MOM.                                                                 
!                                                                      
! MOM is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Jennifer Simeon
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Eric Galbraith
!</REVIEWER>
!
!<OVERVIEW>
! Add-on ocean biogeochemistry module
!</OVERVIEW>
!
!<DESCRIPTION>
! This module has simple implementation of preformed phosphate.
! Where,
!    po4_pre=po4, if z <= mld
!
! It is an optional package that requires TOPAZ, ocmip2_biotic, or ocean_bgc_restore to be running.
!
! Various mixed layer depth options are available and can be
! set via namelist
!
!         1 = kpp mixed layer (default)
!         2 = buoyancy criteria defined mixed layer, where this buoyancy
!              references in situ density
!         3 = buoyancy criteria defined mixed layer, where this buoyancy 
!              references potential density
!
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! No reference yet.
! </REFERENCE>
!
! </INFO>
!

module  ocean_po4_pre_mod

use diag_manager_mod, only: register_diag_field, diag_axis_init
use field_manager_mod, only: fm_get_index
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use fms_mod,                  only: field_exist
use time_manager_mod,         only: get_date, time_type
use mpp_domains_mod,          only: domain2d

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,    only: ocean_thickness_type, ocean_time_type, ocean_density_type
use ocean_tracer_diag_mod, only: calc_mixed_layer_depth, calc_potrho_mixed_layer

implicit none

private

public  :: ocean_po4_pre_bbc
public  :: ocean_po4_pre_end
public  :: ocean_po4_pre_init
public  :: ocean_po4_pre_flux_init
public  :: ocean_po4_pre_sbc
public  :: ocean_po4_pre_source
public  :: ocean_po4_pre_start
public  :: ocean_po4_pre_tracer
public  :: ocean_po4_pre_init_sfc
public  :: ocean_po4_pre_avg_sfc
public  :: ocean_po4_pre_sum_sfc
public  :: ocean_po4_pre_zero_sfc
public  :: ocean_po4_pre_sfc_end

private :: allocate_arrays
private :: set_array

character(len=32), parameter              :: package_name = 'ocean_po4_pre'
character(len=48), parameter              :: mod_name = 'ocean_po4_pre_mod'
character(len=fm_string_len), parameter   :: default_restart_file =   'ocean_po4_pre.res.nc'
character(len=fm_string_len), parameter   :: default_phosphate_name = 'po4'

type mask_region_type
  real, dimension(:,:,:), pointer       :: mask => NULL()
  real, dimension(:), pointer           :: elon => NULL()
  real, dimension(:), pointer           :: nlat => NULL()
  real, dimension(:), pointer           :: slat => NULL()
  real, dimension(:), pointer           :: wlon => NULL()
  logical                               :: coastal_only
  real                                  :: factor
  logical, dimension(:), pointer        :: t_mask => NULL()
end type mask_region_type


type po4_pre_type

  integer                               :: ind_po4_pre = -1
  character(len=fm_string_len)          :: restart_file
  character(len=fm_field_name_len)      :: name
  integer                               :: mld_option 
  type(mask_region_type)                :: po4_pre_mask
  real, _ALLOCATABLE, dimension(:,:)    :: ml_depth  _NULL
  character(len=fm_string_len)          :: phosphate_name

end type po4_pre_type

logical, public :: do_ocean_po4_pre

integer                                 :: package_index
logical                                 :: module_initialized = .false.

character(len=128) :: version = '$Id: ocean_po4_pre.F90,v 20.0 2013/12/14 00:09:36 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

type(po4_pre_type), allocatable, dimension(:)    :: po4_pre
integer                                          :: instances
integer                                          :: indtemp, indsal, indpo4

contains

!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays
! </DESCRIPTION>
!

subroutine allocate_arrays(isd, ied, jsd, jed)

integer, intent(in)     :: isd
integer, intent(in)     :: ied
integer, intent(in)     :: jsd
integer, intent(in)     :: jed

integer :: n

do n = 1, instances
       allocate( po4_pre(n)%ml_depth(isd:ied,jsd:jed) )
       po4_pre(n)%ml_depth(:,:) = 0.0
enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_po4_pre_bbc

end subroutine  ocean_po4_pre_bbc
! </SUBROUTINE> NAME="ocean_po4_pre_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_end">
!
! <DESCRIPTION>
!     Clean up various PO4_PRE quantities for this run.
! </DESCRIPTION>

subroutine ocean_po4_pre_end(isc, iec, jsc, jec, nk, isd, jsd,  &
                             t_prog, grid_dat, grid_tmask, rho_dzt, taup1)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: nk
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: t_prog
integer, intent(in)                                     :: taup1
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

integer :: i, j, k, n
real    :: total_po4_pre

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       integrate the total concentrations of some tracers
!       for the end of the run

!       Use taup1 time index for the end of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
do n = 1, instances
  total_po4_pre = 0.0
  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_po4_pre = total_po4_pre +                                   &
             t_prog(po4_pre(n)%ind_po4_pre)%field(i,j,k,taup1) *      &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_po4_pre)

  write (stdoutunit,*) '  Instance ', trim(po4_pre(n)%name)
  write (stdoutunit,                                              &
         '(/'' Total preformed phosphate  = '',es19.12,'' Gmol-PO4'')')        &
              total_po4_pre * 1.0e-09
enddo

end subroutine  ocean_po4_pre_end
! </SUBROUTINE> NAME="ocean_po4_pre_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>

subroutine ocean_po4_pre_sbc

end subroutine  ocean_po4_pre_sbc
! </SUBROUTINE> NAME="ocean_po4_pre_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocean_po4_pre_flux_init

end subroutine  ocean_po4_pre_flux_init
! </SUBROUTINE> NAME="ocean_po4_pre_flux_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>
subroutine ocean_po4_pre_init

character(len=64), parameter    :: sub_name = 'ocean_po4_pre_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_field_name_len+3)                      :: long_suffix
logical, dimension(12)                                  :: t_mask
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! Initialize the package
package_index = otpm_set_tracer_package(package_name,                 &
     restart_file = default_restart_file,   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

! Check whether to use this package
path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then
  call mpp_error(FATAL, trim(error_header) // 'Could not get number of instances')
endif

! Check some things
write (stdoutunit,*)
if (instances .eq. 0) then
  write (stdoutunit,*) trim(note_header), ' No instances'
  do_ocean_po4_pre = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif

    do_ocean_po4_pre = .true.
endif

module_initialized = .true.

!       Return if we don't want to use this package,
!       after changing the list back
if (.not. do_ocean_po4_pre) then
  return
endif

! allocate storage for po4_pre array
allocate ( po4_pre(instances) )

! loop over the names, saving them into the po4_pre array
do n = 1, instances

  if (fm_get_value(path_to_names, name, index = n)) then
    po4_pre(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = po4_pre(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

  ! PO4_pre
  po4_pre(n)%ind_po4_pre = otpm_set_prog_tracer('po4_pre' // suffix, package_name,   &
       longname = 'Preformed Phosphate' // trim(long_suffix),                      &
       units = 'mol/kg', flux_units = 'mol/m-2/s',             &
       caller = caller_str)
enddo

!       Process the namelists

! Add the package name to the list of good namelists, to be used
! later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

!       Set up the instance namelists
t_mask(:) = .true.

do n = 1, instances
   ! create the instance namelist
  call fm_util_start_namelist(package_name, po4_pre(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('restart_file', default_restart_file)
  call fm_util_set_value('mld_option', 1)
  call fm_util_set_value('phosphate_name', default_phosphate_name)

  call fm_util_end_namelist(package_name, po4_pre(n)%name, check = .true., caller = caller_str)

  ! create some sub-namelists
  call fm_util_start_namelist(trim(package_name), trim(po4_pre(n)%name) // '+po4_pre_mask',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('factor', 0.0)
  call fm_util_set_value('coastal_only', .false.)
  call fm_util_set_value('t_mask', t_mask, size(t_mask))
  call fm_util_set_value('wlon', 0.0, index = 0)
  call fm_util_set_value('elon', 0.0, index = 0)
  call fm_util_set_value('slat', 0.0, index = 0)
  call fm_util_set_value('nlat', 0.0, index = 0)

  call fm_util_end_namelist(trim(package_name), trim(po4_pre(n)%name) // '+po4_pre_mask', caller = caller_str)
enddo

! Check for any errors in the number of fields in the namelists for this package
good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then
  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif

return

end subroutine ocean_po4_pre_init
! </SUBROUTINE> NAME="ocean_po4_pre_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocean_po4_pre_start
! </DESCRIPTION>

subroutine ocean_po4_pre_init_sfc

end subroutine ocean_po4_pre_init_sfc
! </SUBROUTINE> NAME="ocean_po4_pre_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_po4_pre_sum_sfc

end subroutine ocean_po4_pre_sum_sfc
! </SUBROUTINE> NAME="ocean_po4_pre_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_po4_pre_zero_sfc

end subroutine ocean_po4_pre_zero_sfc
! </SUBROUTINE> NAME="ocean_po4_pre_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_po4_pre_avg_sfc

end subroutine ocean_po4_pre_avg_sfc
! </SUBROUTINE> NAME="ocean_po4_pre_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_sfc_end">
!
! <DESCRIPTION>
!       Finish up stuff for surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_po4_pre_sfc_end

end subroutine ocean_po4_pre_sfc_end
! </SUBROUTINE> NAME="ocean_po4_pre_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_source">
!
! <DESCRIPTION>
!     compute the source terms for the PO4_PREs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>

subroutine ocean_po4_pre_source

end subroutine  ocean_po4_pre_source
! </SUBROUTINE> NAME="ocean_po4_pre_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants for a given run
! and allocate diagnostic arrays
! </DESCRIPTION>

subroutine ocean_po4_pre_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,      &
     T_prog, taup1, grid_dat, grid_tmask, grid_kmt,                 &
     grid_xt, grid_yt, rho_dzt)

character(len=64), parameter    :: sub_name = 'ocean_po4_pre_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
integer, intent(in)                                     :: taup1
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
integer, dimension(isd:,jsd:), intent(in)               :: grid_kmt
real, dimension(isd:,jsd:), intent(in)                  :: grid_xt
real, dimension(isd:,jsd:), intent(in)                  :: grid_yt
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

!       Global values to apply the following inhibitions
!       and depletions
!
!  coastal_only : if true, then only apply the changes in
!                 coastal boxes
!  t_mask_len   : parameter giving the number of elements in
!                 the time mask per year (eg., 12 would
!                 imply monthly)
!  t_mask_array : logical array controlling whether to apply
!                 the following inhibitions and depletions to
!                 each time-period (true means set the masks,
!                 false means use the defaults everywhere)
!  num_reg      : number of regions
!  factor       : factor by which to scale the field
!               : in the selected regions
!  wlon : western longitude of region
!  elon : eastern longitude of region
!  slat : southern latitude of region
!  nlat : northern latitude of region
!  mask(imt,jmt)  : mask array (0.0 - alternate, 1.0 - normal)
!
!       Set up a mask array using wlon,elon,nlat,slat
!       (any box with its lon,lat inside the box bounded by
!       wlon,elon,nlat,slat value in mask set to factor).

integer                                 :: i, j, k, l, n
integer                                 :: done
character(len=256)                      :: caller_str
integer                                 :: len_w
integer                                 :: len_e
integer                                 :: len_s
integer                                 :: len_n
real                                    :: total_po4_pre
character(len=25)                       :: po4_tree_path

  integer :: stdoutunit 
  stdoutunit=stdout() 

write(stdoutunit,*) 
write(stdoutunit,*) trim(note_header),                     &
                  'Starting ', trim(package_name), ' module'

!     dynamically allocate the global PO4_PRE arrays
call allocate_arrays(isd, ied, jsd, jed)

!       read in the namelists for each instance
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances
  call fm_util_start_namelist(package_name, po4_pre(n)%name, caller = caller_str)

  po4_pre(n)%restart_file           = fm_util_get_string ('restart_file', scalar = .true.)
  po4_pre(n)%mld_option             = fm_util_get_integer ('mld_option', scalar = .true.)
  po4_pre(n)%phosphate_name         = fm_util_get_string ('phosphate_name', scalar = .true.)


  call fm_util_end_namelist(package_name, po4_pre(n)%name, caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(po4_pre(n)%name) // '+po4_pre_mask', caller = caller_str)

  po4_pre(n)%po4_pre_mask%factor        =  fm_util_get_real          ('factor', scalar = .true.)
  po4_pre(n)%po4_pre_mask%coastal_only  =  fm_util_get_logical       ('coastal_only', scalar = .true.)
  po4_pre(n)%po4_pre_mask%wlon          => fm_util_get_real_array    ('wlon')
  po4_pre(n)%po4_pre_mask%elon          => fm_util_get_real_array    ('elon')
  po4_pre(n)%po4_pre_mask%slat          => fm_util_get_real_array    ('slat')
  po4_pre(n)%po4_pre_mask%nlat          => fm_util_get_real_array    ('nlat')
  po4_pre(n)%po4_pre_mask%t_mask        => fm_util_get_logical_array ('t_mask')

  call fm_util_end_namelist(trim(package_name), trim(po4_pre(n)%name) // '+po4_pre_mask', caller = caller_str)
enddo

!       read in the po4_pre_mask namelist data
do n = 1, instances
  if (associated(po4_pre(n)%po4_pre_mask%wlon)) then
    len_w = size(po4_pre(n)%po4_pre_mask%wlon)
  else
    len_w = 0
  endif
  if (associated(po4_pre(n)%po4_pre_mask%elon)) then
    len_e = size(po4_pre(n)%po4_pre_mask%elon)
  else
    len_e = 0
  endif
  if (associated(po4_pre(n)%po4_pre_mask%slat)) then
    len_s = size(po4_pre(n)%po4_pre_mask%slat)
  else
    len_s = 0
  endif
  if (associated(po4_pre(n)%po4_pre_mask%nlat)) then
    len_n = size(po4_pre(n)%po4_pre_mask%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal')
  endif

  if (size(po4_pre(n)%po4_pre_mask%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12')
  endif

  ! set all of the values to the default
  po4_pre(n)%po4_pre_mask%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process po4_pre_mask array for ', trim(po4_pre(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
  done = 0
  do l = 1, 12
    if (po4_pre(n)%po4_pre_mask%t_mask(l)) then
      if (done .eq. 0) then
         ! set the values via the input values, saving this time index
         ! afterwards
        write (stdoutunit,*) trim(note_header), ' Assigning month ', l
        call set_array(po4_pre(n)%po4_pre_mask%mask(:,:,l), isd, ied, jsd, jed,                                   &
                grid_xt, grid_yt, grid_kmt,                                        &
                len_w, po4_pre(n)%po4_pre_mask%wlon, po4_pre(n)%po4_pre_mask%elon, &
                po4_pre(n)%po4_pre_mask%slat, po4_pre(n)%po4_pre_mask%nlat,        &
                po4_pre(n)%po4_pre_mask%factor, 1.0,                                                          &
                T_prog(po4_pre(n)%ind_po4_pre)%name, po4_pre(n)%po4_pre_mask%coastal_only)
        done = l
      else
         ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          po4_pre(n)%po4_pre_mask%mask(:,:,l) = po4_pre(n)%po4_pre_mask%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

!       integrate the total concentrations of some tracers
!       for the start of the run
indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif
                                                                                                                          
indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif
      
do n = 1,instances

  po4_tree_path = '/ocean_mod/prog_tracers/' // trim(po4_pre(n)%phosphate_name) 
  write(stdoutunit,*)  po4_tree_path

  indpo4 = fm_get_index('/ocean_mod/prog_tracers/' //                 &
  trim(po4_pre(n)%phosphate_name))

  if (indpo4 .le. 0) then
   call mpp_error(FATAL,trim(error_header) // ' Could not get the phosphate index check if running with an ocean biology model')

 endif

enddo

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at start of run'

do n = 1, instances
  total_po4_pre = 0.0
  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_po4_pre = total_po4_pre +                         &
             T_prog(po4_pre(n)%ind_po4_pre)%field(i,j,k,taup1) *  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)

      enddo
    enddo
  enddo

  call mpp_sum(total_po4_pre)
  write (stdoutunit,*) '  Instance ', trim(po4_pre(n)%name)
  write (stdoutunit,                                              &
         '(/'' Total Preformed Phosphate  = '',es19.12,'' (mol/kg)*Gmol-PO4_pre'')')       &
              total_po4_pre * 1.0e-09
enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), 'Tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocean_po4_pre_start
! </SUBROUTINE> NAME="ocean_po4_pre_start"


!#######################################################################
! <SUBROUTINE NAME="ocean_po4_pre_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!
subroutine ocean_po4_pre_tracer(isc, iec, jsc, jec, isd, ied, jsd, jed, nk,     &
                                t_prog, Time, Thickness, Dens, grid_zt, hblt_depth)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
integer, intent(in)                                             :: nk
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: t_prog
type(ocean_time_type), intent(in)                               :: Time
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, dimension(:), intent(in)                                  :: grid_zt
real, dimension(isd:,jsd:), intent(in)                          :: hblt_depth

integer :: i, j, k, n

! set Preformed phosphate
do n = 1, instances

   if (po4_pre(n)%mld_option == 1) then
      po4_pre(n)%ml_depth = hblt_depth
   elseif (po4_pre(n)%mld_option == 2) then
      call calc_mixed_layer_depth(Thickness,                            &
           t_prog(indsal)%field(isd:ied,jsd:jed,:,Time%tau),            &
           t_prog(indtemp)%field(isd:ied,jsd:jed,:,Time%tau),           &
           Dens%rho(isd:ied,jsd:jed,:,Time%tau),                        &
           Dens%pressure_at_depth(isd:ied,jsd:jed,:), &
           po4_pre(n)%ml_depth)
   else
      call calc_potrho_mixed_layer(Thickness, Dens,               &
            potrho_mix_depth= po4_pre(n)%ml_depth)
   endif

  do j = jsc, jec
    do i = isc, iec
      do k = 1,nk
         if (grid_zt(k) <= po4_pre(n)%ml_depth(i,j)) then
           t_prog(po4_pre(n)%ind_po4_pre)%field(i,j,k,Time%taup1) = t_prog(indpo4)%field(i,j,k,Time%taup1) 
         endif
      enddo
    enddo
  enddo

enddo

return

end subroutine  ocean_po4_pre_tracer
! </SUBROUTINE> NAME="ocean_po4_pre_tracer"

!#######################################################################
! <SUBROUTINE NAME="set_array">
!
! <DESCRIPTION>
!       Set up an array covering the model domain with a user-specified
!       value, in user-specified regions. There are a given number of
!       2-d regions specified by the values slat, nlat, wlon and elon.
!       The longitudes are for a cyclic domain, and if wlon and elon
!       are on opposite sides of the cut, the correct thing will
!       be done. Elon is considered to be east of wlon, so if elon is
!       less than wlon, then the region east of elon to the cut will be
!       filled, and the region from the cut to wlon will be filled.
!
!       After setting up the array in this routine, it may prove useful
!       to allow fine-tuning the settings via an array in a namelist.
!
!       Arguments:
!         Input:
!      num_regions = number of user-specified regions which will be
!                    filled
!
!             wlon = 1-d array of western (starting) longitudes for the
!                    rectangular regions
!
!             elon = 1-d array of eastern (ending) longitudes for the
!                    rectangular regions
!
!             slat = 1-d array of southern (starting) latitudes for the
!                    rectangular regions
!
!             nlat = 1-d array of northern (ending) latitudes for the
!                    rectangular regions
!
!                       Note: if slat >= nlat, then nothing is done
!                             for that region
!
!        set_value = the value to assign to array in the user-specified
!                    regions
!
!      unset_value = the value to assign to array outside of the
!                    user-specified regions
!
!             name = character variable used in informative messages
!
!     coastal_only = true to limit changes only to coastal points
!                    (i.e., at least one bordering point is land)
!
!         Output:
!
!            array = 2-d array which will contain the set- and unset-
!                    values. The array is assumed to have a border
!                    one unit wide on all edges, ala MOM. A cyclic
!                    boundary condition will be set if requested.
! </DESCRIPTION>
!
subroutine set_array(array, isd, ied, jsd, jed,                 &
                     xt, yt, kmt,                               &
                     num_regions, wlon_in, elon_in, slat, nlat, &
                     set_value, unset_value, name,              &
                     coastal_only)

integer, intent(in)                             :: isd
integer, intent(in)                             :: ied
integer, intent(in)                             :: jsd
integer, intent(in)                             :: jed
integer, intent(in)                             :: num_regions
real, dimension(isd:ied,jsd:jed), intent(out)   :: array
logical, intent(in)                             :: coastal_only
real, dimension(num_regions), intent(in)        :: elon_in
integer, dimension(isd:ied,jsd:jed), intent(in) :: kmt
character(len=*), intent(in)                    :: name
real, dimension(num_regions), intent(in)        :: nlat
real, intent(in)                                :: set_value
real, dimension(num_regions), intent(in)        :: slat
real, intent(in)                                :: unset_value
real, dimension(num_regions), intent(in)        :: wlon_in
real, dimension(isd:ied,jsd:jed), intent(in)    :: xt
real, dimension(isd:ied,jsd:jed), intent(in)    :: yt

character(len=64), parameter    :: sub_name = 'set_array'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: i, j, n
real, dimension(:), allocatable :: wlon
real, dimension(:), allocatable :: elon

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! save the longitudes in case they need to be modified
allocate(wlon(num_regions))
allocate(elon(num_regions))

wlon(:) = wlon_in(:)
elon(:) = elon_in(:)

! loop over the regions, applying changes as necessary
do n = 1, num_regions

  if (nlat(n) .ge. slat(n)) then
    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header),                          &
                       trim(name), ' region: ', n

    ! make sure that all longitudes are in the range [0,360]
    do while (wlon(n) .gt. 360.0)
      wlon(n) = wlon(n) - 360.0
    enddo
    do while (wlon(n) .lt. 0.0)
      wlon(n) = wlon(n) + 360.0
    enddo
    do while (elon(n) .gt. 360.0)
      elon(n) = elon(n) - 360.0
    enddo
    do while (elon(n) .lt. 0.0)
      elon(n) = elon(n) + 360.0
    enddo
    ! if the southern and northern latitudes are the same, then
    ! find the grid box which encompasses them ...
    if (slat(n) .eq. nlat(n)) then

     call mpp_error(FATAL, trim(error_header) //                &
                    'Equal latitudes not supported')

    elseif (wlon(n) .eq. elon(n)) then

     call mpp_error(FATAL, trim(error_header) //                &
                    'Equal longitudes not supported')
    else
       ! ... else find all boxes where the center lies in the
       ! rectangular region
      do j = jsd, jed
        do i = isd, ied
          if (nlat(n) .ge. yt(i,j) .and.                        &
              slat(n) .le. yt(i,j) .and.                        &
              lon_between(xt(i,j), wlon(n), elon(n))) then
            array(i,j) = set_value
          endif
        enddo
      enddo

    endif

  endif

enddo

!       if desired only apply mask to coastal regions
if (coastal_only) then
  do j = jsd, jed
    do i = isd, ied
      if (kmt(i,j) .ne. 0 .and.                         &
          array(i,j) .eq. set_value) then
!       if all the surrounding points are ocean, then this is not
!       a coastal point, therefore reset the mask
        if (kmt(i-1,j) .ne. 0 .and.                     &
            kmt(i+1,j) .ne. 0 .and.                     &
            kmt(i,j-1) .ne. 0 .and.                     &
            kmt(i,j+1) .ne. 0) then
          array(i,j) = unset_value
        endif
      endif
    enddo
  enddo
endif

!       clean up
deallocate(wlon)
deallocate(elon)

return

contains

!       Return true if w <= x_in <= e, taking into account the
!       periodicity of longitude.
!
!       x_in    = value to test
!
!       w       = west longitude of boundary
!
!       e       = east longitude of boundary
function lon_between(x_in, w, e)

implicit none

logical :: lon_between

real, intent(in)                :: x_in
real, intent(in)                :: w
real, intent(in)                :: e

real                    :: x

! Save input values so we may modify them safely
x = x_in

! make sure that all longitudes are in the range [0,360]
do while (x .gt. 360.0)
  x = x - 360.0
enddo
do while (x .lt. 0.0)
  x = x + 360.0
enddo
 
if (w .gt. e) then
  lon_between = w .le. x .or. x .le. e
else
  lon_between = w .le. x .and. x .le. e
endif

return

end function  lon_between

end subroutine  set_array
! </SUBROUTINE> NAME="set_array"


end module  ocean_po4_pre_mod
