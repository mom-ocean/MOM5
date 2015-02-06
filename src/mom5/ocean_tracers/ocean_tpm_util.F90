module ocean_tpm_util_mod  !{
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Ocean tracer package module pointers module
!</OVERVIEW>
!
!<DESCRIPTION>
! This module allocates a suite of variables used in ocean_tpm
!</DESCRIPTION>
!
! <INFO>
! </INFO>
!

use field_manager_mod, only: fm_string_len, fm_path_name_len
use field_manager_mod, only: fm_get_type, fm_get_index, fm_get_length
use field_manager_mod, only: fm_get_current_list, fm_new_list, fm_change_list
use field_manager_mod, only: fm_new_value, fm_get_value
use field_manager_mod, only: fm_exists
use fms_mod,           only: FATAL, stdout
use mpp_mod,           only: mpp_error
use fm_util_mod,       only: fm_util_default_caller
use fm_util_mod,       only: fm_util_check_for_bad_fields
use fm_util_mod,       only: fm_util_set_caller, fm_util_reset_caller, fm_util_set_no_overwrite
use fm_util_mod,       only: fm_util_reset_no_overwrite, fm_util_set_good_name_list, fm_util_reset_good_name_list
use fm_util_mod,       only: fm_util_get_string_array, fm_util_set_value, fm_util_get_index_string

implicit none

private

!
!       Public routines
!

public  otpm_set_tracer_package
public  otpm_set_prog_tracer
public  otpm_set_diag_tracer

!
!       Public variables
!

!
!       Private routines
!

private  set_prog_value
private  set_prog_value_integer
private  set_prog_value_logical
private  set_prog_value_real
private  set_prog_value_string
private  check_ocean_mod

!
!       private parameters
!

character(len=48), parameter            :: mod_name = 'ocean_tpm_util_mod'

character(len=fm_string_len), parameter :: default_units = ' '
character(len=fm_string_len), parameter :: default_type  = ' '
real, parameter                         :: default_conversion = 1.0
real, parameter                         :: default_offset = 0.0
real, parameter                         :: default_min_tracer = -1.0e+20
real, parameter                         :: default_max_tracer = +1.0e+20
logical, parameter                      :: default_use_only_advection = .false.
real, parameter                         :: default_min_range = 1.0
real, parameter                         :: default_max_range = 0.0
character(len=fm_string_len), parameter :: default_restart_file = 'ocean_tracer.res.nc'
logical, parameter                      :: default_const_init_tracer = .false.
real, parameter                         :: default_const_init_value = 0.0
character(len=fm_string_len), parameter :: default_flux_units = ' '
real, parameter                         :: default_min_flux_range = 1.0
real, parameter                         :: default_max_flux_range = 0.0
real, parameter                         :: default_min_tracer_limit = -1.0e+20
real, parameter                         :: default_max_tracer_limit = +1.0e+20
character(len=fm_string_len), parameter :: default_vert_adv_scheme = 'mdfl_sweby'
character(len=fm_string_len), parameter :: default_horiz_adv_scheme = 'mdfl_sweby'
logical, parameter                      :: default_psom_limit = .false.
integer, parameter                      :: default_ppm_hlimiter = 2
integer, parameter                      :: default_ppm_vlimiter = 2
integer, parameter                      :: default_mdt_scheme = 1

!
!       Private variables
!

character(len=128) :: version = '$Id: ocean_tpm_util.F90,v 20.0 2013/12/14 00:17:18 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!
!        Interface definitions for overloaded routines
!

interface  set_prog_value  !{
  module procedure  set_prog_value_integer
  module procedure  set_prog_value_logical
  module procedure  set_prog_value_real
  module procedure  set_prog_value_string
end interface  !}

contains


!#######################################################################
! <SUBROUTINE NAME="check_ocean_mod">
!
! <DESCRIPTION>
! Be sure that the /ocean_mod hierarchy has been initialized.
! </DESCRIPTION>
!

subroutine check_ocean_mod(caller)  !{

implicit none

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'check_ocean_mod'

!
!       arguments
!

character(len=*), intent(in), optional  :: caller

!
!       Local variables
!

character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
logical, save                                           :: ocean_mod_initialized = .false.

!
!       return if ocean_mod has been initialized
!

if (ocean_mod_initialized) then
  return
endif

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
  
!
! make sure that /ocean_mod exists, just in case there were no inputs in the field table
!

if (fm_new_list('/ocean_mod') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "ocean_mod" list')
endif  

if (fm_new_list('/ocean_mod/GOOD') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "GOOD" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'GOOD', append = .true.)

if (fm_new_list('/ocean_mod/tracer_packages') .le. 0) then 
  call mpp_error(FATAL, trim(error_header) // ' Could not set "tracer packages" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'tracer_packages', append = .true.)

if (fm_new_list('/ocean_mod/prog_tracers') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "prog_tracers" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'prog_tracers', append = .true.)

if (fm_new_list('/ocean_mod/diag_tracers') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "diag_tracers" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'diag_tracers', append = .true.)

if (fm_new_list('/ocean_mod/namelists') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "namelists" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'namelists', append = .true.)
  
if (fm_new_list('/ocean_mod/xland_mix') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "xland_mix" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'xland_mix', append = .true.)

if (fm_new_list('/ocean_mod/xland_insert') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "xland_insert" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'xland_insert', append = .true.)

if (fm_new_list('/ocean_mod/diff_cbt_enhance') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "diff_cbt_enhance" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'diff_cbt_enhance', append = .true.)

if (fm_new_list('/ocean_mod/riverspread') .le. 0) then  
   call mpp_error(FATAL, trim(error_header) // ' Could not set the "riverspread" list')
endif
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'riverspread', append = .true.)

if (fm_new_list('/ocean_mod/rayleigh_damp_table') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "rayleigh_damp_table" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'rayleigh_damp_table', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_info') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_info" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_info', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_src') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_src" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_src', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_int') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_int" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_int', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_ent') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_ent" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_ent', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_01') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_01" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_01', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_02') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_02" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_02', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_03') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_03" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_03', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_04') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_04" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_04', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_05') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_05" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_05', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_06') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_06" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_06', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_07') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_07" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_07', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_08') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_08" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_08', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_09') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_09" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_09', append = .true.)

if (fm_new_list('/ocean_mod/overflow_ofp_prd_line_10') .le. 0) then  
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "overflow_ofp_prd_line_10" list')
endif  
call fm_util_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'overflow_ofp_prd_line_10', append = .true.)


!
! Initialize the good_namelists variable as it may not be otherwise set
!

call fm_util_set_value('/ocean_mod/GOOD/good_namelists', ' ', index = 0)

!
! Check for any errors in the number of fields in the ocean_mod list
!

good_list => fm_util_get_string_array('/ocean_mod/GOOD/good_ocean_mod_list',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then  
  call fm_util_check_for_bad_fields('/ocean_mod', good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else  
  call mpp_error(FATAL,trim(error_header) // ' Empty "good_ocean_mod_list" list')
endif 

ocean_mod_initialized = .true.

return

end subroutine  check_ocean_mod  !}
! </SUBROUTINE> NAME="check_ocean_mod"


!#######################################################################
! <FUNCTION NAME="otpm_set_tracer_package">
!
! <DESCRIPTION>
! Set the values for a tracer package and return its index (0 on error)
!
! </DESCRIPTION>
!

function otpm_set_tracer_package(name, caller, units, conversion, offset,       &
     min_tracer, max_tracer, use_only_advection,                                &
     min_range, max_range, restart_file,                                        &
     const_init_tracer,                                                         &
     const_init_value, flux_units,                                              &
     min_flux_range, max_flux_range,                                            &
     min_tracer_limit, max_tracer_limit,                                        &
     psom_limit, ppm_hlimiter, ppm_vlimiter, mdt_scheme,                        &
     vert_adv_scheme, horiz_adv_scheme)                                         &
         result (pack_index)  !{

implicit none

!
!       Return type
!

integer :: pack_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
character(len=*), intent(in), optional  :: units
character(len=*), intent(in), optional  :: flux_units
character(len=*), intent(in), optional  :: vert_adv_scheme
character(len=*), intent(in), optional  :: horiz_adv_scheme
real, intent(in), optional              :: conversion
character(len=*), intent(in), optional  :: restart_file
logical, intent(in), optional           :: const_init_tracer
logical, intent(in), optional           :: use_only_advection
real, intent(in), optional              :: max_range
real, intent(in), optional              :: min_range
real, intent(in), optional              :: max_flux_range
real, intent(in), optional              :: min_flux_range
real, intent(in), optional              :: max_tracer
real, intent(in), optional              :: min_tracer
real, intent(in), optional              :: max_tracer_limit
real, intent(in), optional              :: min_tracer_limit
real, intent(in), optional              :: offset
real, intent(in), optional              :: const_init_value  
logical, intent(in), optional           :: psom_limit
integer, intent(in), optional           :: ppm_hlimiter
integer, intent(in), optional           :: ppm_vlimiter
integer, intent(in), optional           :: mdt_scheme

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_tracer_package'

!
!       Local variables
!

character(len=fm_path_name_len)                         :: current_list
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_path_name_len)                         :: tracer_package_name
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
logical                                                 :: add_name

  integer :: stdoutunit 
  stdoutunit=stdout() 

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

call check_ocean_mod(caller = trim(sub_name) // caller_str)

write (stdoutunit,*)
write (stdoutunit,*) trim(note_header), ' Processing tracer package ', trim(name)

!
!       Check whether the package already exists. If so, then use that package
!

tracer_package_name = '/ocean_mod/tracer_packages/' // trim(name) // '/'
pack_index = fm_get_index(tracer_package_name)

if (pack_index .le. 0) then  !{

!
!       Set a new tracer package and get its index
!

  pack_index = fm_new_list(tracer_package_name)
  if (pack_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set tracer package')
  endif  !}

endif  !}

!
!       Change to the new tracer package, first saving the current list
!

current_list = fm_get_current_list()
if (current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list(tracer_package_name)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
endif  !}

!
!       Set the array in which to save the valid names for this list,
!       used later for a consistency check. This is used in the fm_util_set_value
!       routines to make the list of valid values
!

call fm_util_set_good_name_list('/ocean_mod/GOOD/tracer_packages/' // trim(name) // '/good_list')

!
!       Set other defaults for the fm_util_set_value routines
!

call fm_util_set_no_overwrite(.true.)
call fm_util_set_caller(caller_str)

!
!       Set the default number of instances (always zero)
!

call fm_util_set_value('names', ' ', index = 0)

!
!       Set various values to given values, or to defaults if not given
!

if (present(units)) then  !{
  call fm_util_set_value('units', units)
else  !}{
  call fm_util_set_value('units', default_units, no_create = .true.)
endif  !}

if (present(conversion)) then  !{
  call fm_util_set_value('conversion', conversion)
else  !}{
  call fm_util_set_value('conversion', default_conversion, no_create = .true.)
endif  !}

if (present(offset)) then  !{
  call fm_util_set_value('offset', offset)
else  !}{
  call fm_util_set_value('offset', default_offset, no_create = .true.)
endif  !}

if (present(min_tracer)) then  !{
  call fm_util_set_value('min_tracer', min_tracer)
else  !}{
  call fm_util_set_value('min_tracer', default_min_tracer, no_create = .true.)
endif  !}

if (present(max_tracer)) then  !{
  call fm_util_set_value('max_tracer', max_tracer)
else  !}{
  call fm_util_set_value('max_tracer', default_max_tracer, no_create = .true.)
endif  !}

if (present(min_range)) then  !{
  call fm_util_set_value('min_range', min_range)
else  !}{
  call fm_util_set_value('min_range', default_min_range, no_create = .true.)
endif  !}

if (present(max_range)) then  !{
  call fm_util_set_value('max_range', max_range)
else  !}{
  call fm_util_set_value('max_range', default_max_range, no_create = .true.)
endif  !}

if (present(use_only_advection)) then  !{
  call fm_util_set_value('use_only_advection', use_only_advection)
else  !}{
  call fm_util_set_value('use_only_advection', default_use_only_advection, no_create = .true.)
endif  !}

if (present(restart_file)) then  !{
  call fm_util_set_value('restart_file', restart_file)
else  !}{
  call fm_util_set_value('restart_file', default_restart_file, no_create = .true.)
endif  !}

if (present(const_init_tracer)) then  !{
  call fm_util_set_value('const_init_tracer', const_init_tracer)
else  !}{
  call fm_util_set_value('const_init_tracer', default_const_init_tracer, no_create = .true.)
endif  !}

if (present(const_init_value)) then  !{
  call fm_util_set_value('const_init_value', const_init_value)
else  !}{
  call fm_util_set_value('const_init_value', default_const_init_value, no_create = .true.)
endif  !}

if (present(psom_limit)) then  !{
  call fm_util_set_value('psom_limit', psom_limit)
else  !}{
  call fm_util_set_value('psom_limit', default_psom_limit, no_create = .true.)
endif  !}

if (present(ppm_hlimiter)) then  !{
  call fm_util_set_value('ppm_hlimiter', ppm_hlimiter)
else  !}{
  call fm_util_set_value('ppm_hlimiter', default_ppm_hlimiter, no_create = .true.)
endif  !}

if (present(ppm_vlimiter)) then  !{
  call fm_util_set_value('ppm_vlimiter', ppm_vlimiter)
else  !}{
  call fm_util_set_value('ppm_vlimiter', default_ppm_vlimiter, no_create = .true.)
endif  !}

if (present(mdt_scheme)) then  !{
  call fm_util_set_value('mdt_scheme', mdt_scheme)
else  !}{
  call fm_util_set_value('mdt_scheme', default_mdt_scheme, no_create = .true.)
endif  !}

if (present(flux_units)) then  !{
  call fm_util_set_value('flux_units', flux_units)
else  !}{
  call fm_util_set_value('flux_units', default_flux_units, no_create = .true.)
endif  !}

if (present(min_flux_range)) then  !{
  call fm_util_set_value('min_flux_range', min_flux_range)
else  !}{
  call fm_util_set_value('min_flux_range', default_min_flux_range, no_create = .true.)
endif  !}

if (present(max_flux_range)) then  !{
  call fm_util_set_value('max_flux_range', max_flux_range)
else  !}{
  call fm_util_set_value('max_flux_range', default_max_flux_range, no_create = .true.)
endif  !}

if (present(min_tracer_limit)) then  !{
  call fm_util_set_value('min_tracer_limit', min_tracer_limit)
else  !}{
  call fm_util_set_value('min_tracer_limit', default_min_tracer_limit, no_create = .true.)
endif  !}

if (present(max_tracer_limit)) then  !{
  call fm_util_set_value('max_tracer_limit', max_tracer_limit)
else  !}{
  call fm_util_set_value('max_tracer_limit', default_max_tracer_limit, no_create = .true.)
endif  !}

if (present(vert_adv_scheme)) then  !{
  call fm_util_set_value('vertical-advection-scheme', vert_adv_scheme)
else  !}{
  call fm_util_set_value('vertical-advection-scheme', default_vert_adv_scheme, no_create = .true.)
endif  !}

if (present(horiz_adv_scheme)) then  !{
  call fm_util_set_value('horizontal-advection-scheme', horiz_adv_scheme)
else  !}{
  call fm_util_set_value('horizontal-advection-scheme', default_horiz_adv_scheme, no_create = .true.)
endif  !}

!
!       Reset the defaults for the fm_util_set_value calls
!

call fm_util_reset_good_name_list
call fm_util_reset_no_overwrite
call fm_util_reset_caller

!
!       Change back to the saved current list
!

if (.not. fm_change_list(current_list)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
endif  !}

!
!       Check for any errors in the number of fields in this list
!

if (caller_str .eq. ' ') then  !{
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
endif  !}
good_list => fm_util_get_string_array('/ocean_mod/GOOD/tracer_packages/' // trim(name) // '/good_list', &
     caller = caller_str)
if (associated(good_list)) then  !{
  call fm_util_check_for_bad_fields('/ocean_mod/tracer_packages/' // trim(name), good_list, caller = caller_str)
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
endif  !}

!
!       Add the package name to the list of good packages (if not already there), to be used
!       later for a consistency check
!

if (fm_exists('/ocean_mod/GOOD/good_tracer_packages')) then  !{
  add_name = fm_util_get_index_string('/ocean_mod/GOOD/good_tracer_packages', name,           &
     caller = caller_str) .le. 0                ! true if name does not exist in string array
else  !}{
  add_name = .true.                             ! always add to new list
endif  !}
if (add_name) then  !{
  if (fm_new_value('/ocean_mod/GOOD/good_tracer_packages', name, append = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                                         &
         ' Could not add ' // trim(name) // ' to "good_tracer_packages" list')
  endif  !}
endif  !}

return

end function otpm_set_tracer_package  !}
! </FUNCTION> NAME="otpm_set_tracer_package"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_integer">
!
! <DESCRIPTION>
! Set an integer value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_integer(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
integer, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_integer'

!
!       Local variables
!

integer                         :: integer_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine fm_util_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) value given in argument list

  call fm_util_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, integer_value)) then  !{

    call fm_util_set_value(name, integer_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_integer  !}
! </SUBROUTINE> NAME="set_prog_value_integer"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_logical">
!
! <DESCRIPTION>
! Set a logical value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_logical(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
logical, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_logical'

!
!       Local variables
!

logical                         :: logical_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine fm_util_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) value given in argument list

  call fm_util_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, logical_value)) then  !{

    call fm_util_set_value(name, logical_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_logical  !}
! </SUBROUTINE> NAME="set_prog_value_logical"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_real">
!
! <DESCRIPTION>
! Set a real value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_real(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
real, intent(in)                        :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_real'

!
!       Local variables
!

real                            :: real_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine fm_util_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) value given in argument list

  call fm_util_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, real_value)) then  !{

    call fm_util_set_value(name, real_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_real  !}
! </SUBROUTINE> NAME="set_prog_value_real"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_string">
!
! <DESCRIPTION>
! Set a string value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_string(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
character(len=*), intent(in)            :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_string'

!
!       Local variables
!

character(len=fm_string_len)    :: string_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine fm_util_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) value given in argument list

  call fm_util_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the fm_util_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, string_value)) then  !{

    call fm_util_set_value(name, string_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_string  !}
! </SUBROUTINE> NAME="set_prog_value_string"


!#######################################################################
! <FUNCTION NAME="otpm_set_prog_tracer">
!
! <DESCRIPTION>
! Set the values for a prog tracer and return its index (0 on error)
! </DESCRIPTION>
!
function otpm_set_prog_tracer(name, package, caller, longname, units,   &
     type, conversion, offset,                                          &
     min_tracer, max_tracer, use_only_advection,                        &
     min_range, max_range, restart_file,                                &
     const_init_tracer,                                                 &
     const_init_value, flux_units,                                      &
     min_flux_range, max_flux_range,                                    &
     min_tracer_limit, max_tracer_limit,                                &
     psom_limit, ppm_hlimiter, ppm_vlimiter, mdt_scheme,                &
     vert_adv_scheme, horiz_adv_scheme)                                 &
         result (prog_index)  !{

implicit none

!
!       Return type
!

integer :: prog_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package
character(len=*), intent(in), optional  :: caller
character(len=*), intent(in), optional  :: units
character(len=*), intent(in), optional  :: type
character(len=*), intent(in), optional  :: flux_units
character(len=*), intent(in), optional  :: longname
character(len=*), intent(in), optional  :: vert_adv_scheme
character(len=*), intent(in), optional  :: horiz_adv_scheme
real, intent(in), optional              :: conversion
character(len=*), intent(in), optional  :: restart_file
logical, intent(in), optional           :: const_init_tracer
logical, intent(in), optional           :: use_only_advection
real, intent(in), optional              :: max_range
real, intent(in), optional              :: min_range
real, intent(in), optional              :: max_flux_range
real, intent(in), optional              :: min_flux_range
real, intent(in), optional              :: max_tracer
real, intent(in), optional              :: min_tracer
real, intent(in), optional              :: max_tracer_limit
real, intent(in), optional              :: min_tracer_limit
real, intent(in), optional              :: offset
real, intent(in), optional              :: const_init_value  
logical, intent(in), optional           :: psom_limit  
integer, intent(in), optional           :: ppm_hlimiter
integer, intent(in), optional           :: ppm_vlimiter
integer, intent(in), optional           :: mdt_scheme

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_prog_tracer'

!
!       Local variables
!

character(len=fm_path_name_len)                         :: current_list
character(len=fm_path_name_len)                         :: package_name
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_path_name_len)                         :: tracer_name
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
logical                                                 :: add_name

  integer :: stdoutunit 
  stdoutunit=stdout() 

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

call check_ocean_mod(caller = trim(sub_name) // caller_str)

!
!       check the package name
!

if (package .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package given')
endif  !}
package_name = '/ocean_mod/tracer_packages/' // trim(package) // '/'
if (fm_get_type(package_name) .ne. 'list') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Package does not exist or is not a list: ' // trim(package))
endif  !}

!
!       Begin processing
!

write (stdoutunit,*)
write (stdoutunit,*) trim(note_header), ' Processing prog tracer ', trim(name)

!
!       Check whether the tracer already exists. If so, then use that tracer
!

tracer_name = '/ocean_mod/prog_tracers/' // trim(name) // '/'
prog_index = fm_get_index(tracer_name)

if (prog_index .le. 0) then  !{

!
!       Set a new prog tracer and get its index
!

  prog_index = fm_new_list(tracer_name)
  if (prog_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set prog tracer ' // trim(name))
  endif  !}

endif  !}

!
!       Change to the new tracer, first saving the current list
!

current_list = fm_get_current_list()
if (current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list(tracer_name)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
endif  !}

!
!       Set the array in which to save the valid names for this list,
!       used later for a consistency check. This is used in the fm_util_set_value
!       routines to make the list of valid values
!

call fm_util_set_good_name_list('/ocean_mod/GOOD/prog_tracers/' // trim(name) // '/good_list')

!
!       Set other defaults for the fm_util_set_value routines
!

call fm_util_set_caller(caller_str)

!
!       When the following is set to true, fm_util_set_value will not overwrite
!       any values already set in the tracer tree If there is no
!       value present, then a new entry will be created in the tracer tree and
!       the value supplied will be set.
!

call fm_util_set_no_overwrite(.true.)

!
!       Set various values to given values, or to defaults if not given
!

!
!       The longname is distinct here in that there is no option for a package
!       default. Hence, the precedence of values is:
!               1) a value set in the field table
!               2) an optional argument given to this subroutine
!               3) the tracer name
!

if (present(longname)) then  !{
  call fm_util_set_value('longname', longname)
else  !}{
  call fm_util_set_value('longname', name)
endif  !}

!
!       The precedence of values to use in set_prog_value is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value passed to it in the argument list
!       This subroutine will preferentially supply the given optional
!       argument over the module-wide default value
!

if (present(units)) then  !{
  call set_prog_value('units', package_name, units)
else  !}{
  call set_prog_value('units', package_name, default_units)
endif  !}

if (present(type)) then  !{
   call set_prog_value('type', package_name, type)
else  !}{
   call set_prog_value('type', package_name, default_type)
endif  !}

if (present(conversion)) then  !}{
  call set_prog_value('conversion', package_name, conversion)
else  !}{
  call set_prog_value('conversion', package_name, default_conversion)
endif  !}

if (present(offset)) then  !}{
  call set_prog_value('offset', package_name, offset)
else  !}{
  call set_prog_value('offset', package_name, default_offset)
endif  !}

if (present(min_tracer)) then  !}{
  call set_prog_value('min_tracer', package_name, min_tracer)
else  !}{
  call set_prog_value('min_tracer', package_name, default_min_tracer)
endif  !}

if (present(max_tracer)) then  !}{
  call set_prog_value('max_tracer', package_name, max_tracer)
else  !}{
  call set_prog_value('max_tracer', package_name, default_max_tracer)
endif  !}

if (present(min_range)) then  !}{
  call set_prog_value('min_range', package_name, min_range)
else  !}{
  call set_prog_value('min_range', package_name, default_min_range)
endif  !}

if (present(max_range)) then  !}{
  call set_prog_value('max_range', package_name, max_range)
else  !}{
  call set_prog_value('max_range', package_name, default_max_range)
endif  !}

if (present(use_only_advection)) then  !}{
  call set_prog_value('use_only_advection', package_name, use_only_advection)
else  !}{
  call set_prog_value('use_only_advection', package_name, default_use_only_advection)
endif  !}

if (present(restart_file)) then  !}{
  call set_prog_value('restart_file', package_name, restart_file)
else  !}{
  call set_prog_value('restart_file', package_name, default_restart_file)
endif  !}

if (present(const_init_tracer)) then  !}{
  call set_prog_value('const_init_tracer', package_name, const_init_tracer)
else  !}{
  call set_prog_value('const_init_tracer', package_name, default_const_init_tracer)
endif  !}

if (present(const_init_value)) then  !}{
  call set_prog_value('const_init_value', package_name, const_init_value)
else  !}{
  call set_prog_value('const_init_value', package_name, default_const_init_value)
endif  !}

if (present(psom_limit)) then  !}{
  call set_prog_value('psom_limit', package_name, psom_limit)
else  !}{
  call set_prog_value('psom_limit', package_name, default_psom_limit)
endif  !}

if (present(ppm_hlimiter)) then  !}{
  call set_prog_value('ppm_hlimiter', package_name, ppm_hlimiter)
else  !}{
  call set_prog_value('ppm_hlimiter', package_name, default_ppm_hlimiter)
endif  !}

if (present(ppm_vlimiter)) then  !}{
  call set_prog_value('ppm_vlimiter', package_name, ppm_vlimiter)
else  !}{
  call set_prog_value('ppm_vlimiter', package_name, default_ppm_vlimiter)
endif  !}

if (present(mdt_scheme)) then  !}{
  call set_prog_value('mdt_scheme', package_name, mdt_scheme)
else  !}{
  call set_prog_value('mdt_scheme', package_name, default_mdt_scheme)
endif  !}

if (present(flux_units)) then  !}{
  call set_prog_value('flux_units', package_name, flux_units)
else  !}{
  call set_prog_value('flux_units', package_name, default_flux_units)
endif  !}

if (present(min_flux_range)) then  !}{
  call set_prog_value('min_flux_range', package_name, min_flux_range)
else  !}{
  call set_prog_value('min_flux_range', package_name, default_min_flux_range)
endif  !}

if (present(max_flux_range)) then  !}{
  call set_prog_value('max_flux_range', package_name, max_flux_range)
else  !}{
  call set_prog_value('max_flux_range', package_name, default_max_flux_range)
endif  !}

if (present(min_tracer_limit)) then  !}{
  call set_prog_value('min_tracer_limit', package_name, min_tracer_limit)
else  !}{
  call set_prog_value('min_tracer_limit', package_name, default_min_tracer_limit)
endif  !}

if (present(max_tracer_limit)) then  !}{
  call set_prog_value('max_tracer_limit', package_name, max_tracer_limit)
else  !}{
  call set_prog_value('max_tracer_limit', package_name, default_max_tracer_limit)
endif  !}

if (present(vert_adv_scheme)) then  !}{
  call set_prog_value('vertical-advection-scheme', package_name, vert_adv_scheme)
else  !}{
  call set_prog_value('vertical-advection-scheme', package_name, default_vert_adv_scheme)
endif  !}

if (present(horiz_adv_scheme)) then  !}{
  call set_prog_value('horizontal-advection-scheme', package_name, horiz_adv_scheme)
else  !}{
  call set_prog_value('horizontal-advection-scheme', package_name, default_horiz_adv_scheme)
endif  !}

!
!       Reset the defaults for the fm_util_set_value calls
!

call fm_util_reset_good_name_list
call fm_util_reset_no_overwrite
call fm_util_reset_caller

!
!       Change back to the saved current list
!

if (.not. fm_change_list(current_list)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
endif  !}

!
!       Check for any errors in the number of fields in this list
!

if (caller_str .eq. ' ') then  !{
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
endif  !}
good_list => fm_util_get_string_array('/ocean_mod/GOOD/prog_tracers/' // trim(name) // '/good_list', &
     caller = caller_str)
if (associated(good_list)) then  !{
  call fm_util_check_for_bad_fields('/ocean_mod/prog_tracers/' // trim(name), good_list, caller = caller_str)
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
endif  !}

!
!       Add the tracer name to the list of good tracers (if not already there), to be used
!       later for a consistency check
!

if (fm_exists('/ocean_mod/GOOD/good_prog_tracers')) then  !{
  add_name = fm_util_get_index_string('/ocean_mod/GOOD/good_prog_tracers', name,           &
     caller = caller_str) .le. 0                ! true if name does not exist in string array
else  !}{
  add_name = .true.                             ! always add to new list
endif  !}
if (add_name) then  !{
  if (fm_new_value('/ocean_mod/GOOD/good_prog_tracers', name, append = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                                         &
         ' Could not add ' // trim(name) // ' to "good_prog_tracers" list')
  endif  !}
endif  !}

return

end function otpm_set_prog_tracer  !}
! </FUNCTION> NAME="otpm_set_prog_tracer"


!#######################################################################
! <FUNCTION NAME="otpm_set_diag_tracer">
!
! <DESCRIPTION>
! Set the values for a diag tracer and return its index (0 on error)
! </DESCRIPTION>
!
function otpm_set_diag_tracer(name, caller, longname, units,    &
     type, conversion, offset, min_tracer, max_tracer,          &
     min_range, max_range, restart_file,                        &
     const_init_tracer,                                         &
     const_init_value)                                          &
         result (diag_index)  !{

implicit none

!
!       Return type
!

integer :: diag_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
character(len=*), intent(in), optional  :: units
character(len=*), intent(in), optional  :: type
character(len=*), intent(in), optional  :: longname
real, intent(in), optional              :: conversion
character(len=*), intent(in), optional  :: restart_file
logical, intent(in), optional           :: const_init_tracer
real, intent(in), optional              :: max_range
real, intent(in), optional              :: min_range
real, intent(in), optional              :: max_tracer
real, intent(in), optional              :: min_tracer
real, intent(in), optional              :: offset
real, intent(in), optional              :: const_init_value  

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_diag_tracer'

!
!       Local variables
!

character(len=fm_path_name_len)                         :: current_list
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
character(len=fm_path_name_len)                         :: tracer_name
logical                                                 :: add_name

  integer :: stdoutunit 
  stdoutunit=stdout() 

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
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

call check_ocean_mod(caller = trim(sub_name) // caller_str)

!
!       Check whether the tracer already exists. If so, then use that tracer
!

tracer_name = '/ocean_mod/diag_tracers/' // trim(name) // '/'
diag_index = fm_get_index(tracer_name)

if (diag_index .le. 0) then  !{

!
!       Set a new diag tracer and get its index
!

  diag_index = fm_new_list(tracer_name)
  if (diag_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set diag tracer ' // trim(name))
  endif  !}

endif  !}

write (stdoutunit,*)
write (stdoutunit,*) trim(note_header), ' Processing diag tracer ', trim(name)

!
!       Change to the new tracer, first saving the current list
!

current_list = fm_get_current_list()
if (current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list(tracer_name)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
endif  !}

!
!       Set the array in which to save the valid names for this list,
!       used later for a consistency check. This is used in the fm_util_set_value
!       routines to make the list of valid values
!

call fm_util_set_good_name_list('/ocean_mod/GOOD/diag_tracers/' // trim(name) // '/good_list')

!
!       Set other defaults for the fm_util_set_value routines
!

call fm_util_set_no_overwrite(.true.)
call fm_util_set_caller(caller_str)

!
!       Set various values to given values, or to defaults if not given
!

if (present(longname)) then  !{
  call fm_util_set_value('longname', longname)
else  !}{
  call fm_util_set_value('longname', name)
endif  !}

if (present(units)) then  !{
  call fm_util_set_value('units', units)
else  !}{
  call fm_util_set_value('units', default_units)
endif  !}

if (present(type)) then  !{
   call fm_util_set_value('type', type)
else  !}{
   call fm_util_set_value('type', default_type)
endif  !}

if (present(conversion)) then  !{
  call fm_util_set_value('conversion', conversion)
else  !}{
  call fm_util_set_value('conversion', default_conversion)
endif  !}

if (present(offset)) then  !{
  call fm_util_set_value('offset', offset)
else  !}{
  call fm_util_set_value('offset', default_offset)
endif  !}

if (present(min_tracer)) then  !{
  call fm_util_set_value('min_tracer', min_tracer)
else  !}{
  call fm_util_set_value('min_tracer', default_min_tracer)
endif  !}

if (present(max_tracer)) then  !{
  call fm_util_set_value('max_tracer', max_tracer)
else  !}{
  call fm_util_set_value('max_tracer', default_max_tracer)
endif  !}

if (present(min_range)) then  !{
  call fm_util_set_value('min_range', min_range)
else  !}{
  call fm_util_set_value('min_range', default_min_range)
endif  !}

if (present(max_range)) then  !{
  call fm_util_set_value('max_range', max_range)
else  !}{
  call fm_util_set_value('max_range', default_max_range)
endif  !}

if (present(restart_file)) then  !{
  call fm_util_set_value('restart_file', restart_file)
else  !}{
  call fm_util_set_value('restart_file', ' ')
endif  !}

if (present(const_init_tracer)) then  !{
  call fm_util_set_value('const_init_tracer', const_init_tracer)
else  !}{
  call fm_util_set_value('const_init_tracer', default_const_init_tracer)
endif  !}

if (present(const_init_value)) then  !{
  call fm_util_set_value('const_init_value', const_init_value)
else  !}{
  call fm_util_set_value('const_init_value', default_const_init_value)
endif  !}

!
!       Reset the defaults for the fm_util_set_value calls
!

call fm_util_reset_good_name_list
call fm_util_reset_no_overwrite
call fm_util_reset_caller

!
!       Change back to the saved current list
!

if (.not. fm_change_list(current_list)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
endif  !}

!
!       Check for any errors in the number of fields in this list
!

if (caller_str .eq. ' ') then  !{
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
endif  !}
good_list => fm_util_get_string_array('/ocean_mod/GOOD/diag_tracers/' // trim(name) // '/good_list', &
     caller = caller_str)
if (associated(good_list)) then  !{
  call fm_util_check_for_bad_fields(tracer_name, good_list, caller = caller_str)
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
endif  !}

!
!       Add the tracer name to the list of good tracers (if not already there), to be used
!       later for a consistency check
!

if (fm_exists('/ocean_mod/GOOD/good_diag_tracers')) then  !{
  add_name = fm_util_get_index_string('/ocean_mod/GOOD/good_diag_tracers', name,           &
     caller = caller_str) .le. 0                ! true if name does not exist in string array
else  !}{
  add_name = .true.                             ! always add to new list
endif  !}
if (add_name) then  !{
  if (fm_new_value('/ocean_mod/GOOD/good_diag_tracers', name, append = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                                         &
         ' Could not add ' // trim(name) // ' to "good_diag_tracers" list')
  endif  !}
endif  !}

return

end function otpm_set_diag_tracer  !}
! </FUNCTION> NAME="otpm_set_diag_tracer"


end module ocean_tpm_util_mod  !}
