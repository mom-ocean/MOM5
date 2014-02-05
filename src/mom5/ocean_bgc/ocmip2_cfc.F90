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
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: CFC module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 CFC
!       simulations as outlined in the CFC-HOWTO documentation,
!       revision 1.6, 1999/04/29.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/CFC/HOWTO-CFC.html
! </REFERENCE>
! </INFO>
!

module ocmip2_cfc_mod

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux
use time_manager_mod,   only: time_type
use field_manager_mod,  only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,  only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,            only: field_exist
use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use constants_mod,      only: WTMCFC11, WTMCFC12

use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use mpp_mod,            only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use coupler_types_mod,  only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_grid_type, ocean_time_type
use diag_manager_mod,   only: register_diag_field
use field_manager_mod,  only: fm_get_index
use ocean_util_mod,     only: diagnose_2d, diagnose_2d_comp

implicit none

private

public  :: ocmip2_cfc_bbc
public  :: ocmip2_cfc_end
public  :: ocmip2_cfc_init
public  :: ocmip2_cfc_flux_init
public  :: ocmip2_cfc_sbc
public  :: ocmip2_cfc_source
public  :: ocmip2_cfc_start
public  :: ocmip2_cfc_init_sfc
public  :: ocmip2_cfc_avg_sfc
public  :: ocmip2_cfc_sum_sfc
public  :: ocmip2_cfc_zero_sfc
public  :: ocmip2_cfc_sfc_end

private :: allocate_arrays

character(len=fm_field_name_len), parameter     :: package_name = 'ocmip2_cfc'
character(len=48), parameter                    :: mod_name = 'ocmip2_cfc_mod'
character(len=48), parameter                    :: diag_name = 'ocean_ocmip2_cfc'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocmip2_cfc.res.nc'
character(len=fm_string_len), parameter         :: default_ice_restart_file = 'ice_ocmip2_cfc.res.nc'
character(len=fm_string_len), parameter         :: default_ocean_restart_file = 'ocmip2_cfc_airsea_flux.res.nc'

integer, parameter :: max_cfc_rec = 1200

type cfc_type

  real                                  :: sc_11_0
  real                                  :: sc_11_1
  real                                  :: sc_11_2
  real                                  :: sc_11_3
  real                                  :: d1_11
  real                                  :: d2_11
  real                                  :: d3_11
  real                                  :: d4_11
  real                                  :: e1_11
  real                                  :: e2_11
  real                                  :: e3_11
  real                                  :: sc_12_0
  real                                  :: sc_12_1
  real                                  :: sc_12_2
  real                                  :: sc_12_3
  real                                  :: d1_12
  real                                  :: d2_12
  real                                  :: d3_12
  real                                  :: d4_12
  real                                  :: e1_12
  real                                  :: e2_12
  real                                  :: e3_12
  integer                               :: id_sc_11 = -1
  integer                               :: id_alpha_11 = -1
  integer                               :: id_sc_12 = -1
  integer                               :: id_alpha_12 = -1
  integer                               :: ind_cfc_11
  integer                               :: ind_cfc_12
  integer                               :: ind_cfc_11_flux
  integer                               :: ind_cfc_12_flux
  character(len=fm_field_name_len)      :: name
  real, _ALLOCATABLE, dimension(:,:)    :: sc_11  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: sc_12  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha_11  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha_12  _NULL
  integer                               :: id_sfc_flux_cfc_11 = -1
  integer                               :: id_sfc_flux_cfc_12 = -1

end type cfc_type

logical, public :: do_ocmip2_cfc

type(cfc_type), allocatable, dimension(:)       :: cfc
integer                                         :: instances
integer                                         :: package_index
logical                                         :: module_initialized = .false.
real, allocatable, dimension(:,:)               :: sc_no_term
integer                                         :: indsal
integer                                         :: indtemp

character(len=128) :: version = '$Id: ocmip2_cfc.F90,v 20.0 2013/12/14 00:09:42 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

contains

!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays
! </DESCRIPTION>
!

subroutine allocate_arrays(isc, iec, jsc, jec)

integer, intent(in)     :: isc
integer, intent(in)     :: iec
integer, intent(in)     :: jsc
integer, intent(in)     :: jec

integer :: n

allocate( sc_no_term(isc:iec,jsc:jec) )

!       allocate cfc array elements
do n = 1, instances
  allocate( cfc(n)%sc_11(isc:iec,jsc:jec) )
  allocate( cfc(n)%alpha_11(isc:iec,jsc:jec) )
  allocate( cfc(n)%sc_12(isc:iec,jsc:jec) )
  allocate( cfc(n)%alpha_12(isc:iec,jsc:jec) )
enddo

sc_no_term(:,:) = 0.0

do n = 1, instances
  cfc(n)%sc_11(:,:) = 0.0
  cfc(n)%alpha_11(:,:) = 0.0
  cfc(n)%sc_12(:,:) = 0.0
  cfc(n)%alpha_12(:,:) = 0.0
enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!
subroutine ocmip2_cfc_bbc

end subroutine  ocmip2_cfc_bbc
! </SUBROUTINE> NAME="ocmip2_cfc_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_end">
!
! <DESCRIPTION>
!     Clean up various CFC quantities for this run.
! </DESCRIPTION>
!
subroutine ocmip2_cfc_end(isc, iec, jsc, jec, nk, isd, jsd,  &
     T_prog, grid_dat, grid_tmask, rho_dzt, taup1)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
integer, intent(in)                                     :: taup1
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_end'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: i, j, k, n
real    :: total_cfc_11
real    :: total_cfc_12

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       integrate the total concentrations of some tracers
!       for the end of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at end of run'

do n = 1, instances

  total_cfc_11 = 0.0
  total_cfc_12 = 0.0

  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_cfc_11 = total_cfc_11 +                           &
             t_prog(cfc(n)%ind_cfc_11)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_cfc_12 = total_cfc_12 +                           &
             t_prog(cfc(n)%ind_cfc_12)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_cfc_11)
  call mpp_sum(total_cfc_12)

  write (stdoutunit,*) '  Instance ', trim(cfc(n)%name)
  write (stdoutunit,                                      &
       '(/'' Total CFC-11  = '',es19.12,'' Gmol'')')    &
       total_cfc_11 * 1.0e-09
  write (stdoutunit,                                      &
       '(/'' Total CFC-12  = '',es19.12,'' Gmol'')')    &
       total_cfc_12 * 1.0e-09

enddo

return
end subroutine  ocmip2_cfc_end
! </SUBROUTINE> NAME="ocmip2_cfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!
subroutine ocmip2_cfc_sbc(isc, iec, jsc, jec,     &
     isc_bnd, jsc_bnd,                                 &
     T_prog, Grid, Time, ice_ocean_boundary_fluxes)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: isc_bnd
integer, intent(in)                                             :: jsc_bnd
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
type(coupler_2d_bc_type), intent(in)                            :: ice_ocean_boundary_fluxes

integer :: i_bnd_off
integer :: j_bnd_off
integer :: i, j, n

!     use the surface fluxes from the coupler
!       stf is in mol/m^2/s, flux from coupler is positive upwards
i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
      t_prog(cfc(n)%ind_cfc_11)%stf(i,j) =                              &
            -ice_ocean_boundary_fluxes%bc(cfc(n)%ind_cfc_11_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
      t_prog(cfc(n)%ind_cfc_12)%stf(i,j) =                              &
            -ice_ocean_boundary_fluxes%bc(cfc(n)%ind_cfc_12_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d(Time, Grid, cfc(n)%id_sfc_flux_cfc_11, t_prog(cfc(n)%ind_cfc_11)%stf(:,:))
   call diagnose_2d(Time, Grid, cfc(n)%id_sfc_flux_cfc_12, t_prog(cfc(n)%ind_cfc_12)%stf(:,:))
enddo

return

end subroutine  ocmip2_cfc_sbc
! </SUBROUTINE> NAME="ocmip2_cfc_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>
subroutine ocmip2_cfc_flux_init

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_flux_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=256)                                      :: caller_str

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       First, perform some initialization if this module has not been
!       initialized because the normal initialization routine will
!       not have been called as part of the normal ocean model
!       initialization if this is an Atmosphere pe of a coupled
!       model running in concurrent mode
if (.not. module_initialized) then

   ! Initialize the package
  package_index = otpm_set_tracer_package(package_name,            &
       restart_file = default_restart_file,                        &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')

  ! Check whether to use this package
  path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
  instances = fm_get_length(path_to_names)
  if (instances .lt. 0) then
    call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
  endif

  write (stdoutunit,*)
  if (instances .eq. 0) then
    write (stdoutunit,*) trim(note_header), ' No instances'
    do_ocmip2_cfc = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocmip2_cfc = .true.
  endif

  module_initialized = .true.

endif

!       Return if we don't want to use this package
if (.not. do_ocmip2_cfc) then
  return
endif

if (.not. allocated(cfc)) then

   ! allocate storage for cfc array
  allocate ( cfc(instances) )

  ! loop over the names, saving them into the cfc array
  do n = 1, instances

    if (fm_get_value(path_to_names, name, index = n)) then
      cfc(n)%name = name
    else
      write (name,*) n
      call mpp_error(FATAL, trim(error_header) //        &
           'Bad field name for index ' // trim(name))
    endif

  enddo

endif

! Set up the ocean-atmosphere gas flux fields
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = cfc(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // name
  endif

  ! Coupler fluxes
  cfc(n)%ind_cfc_11_flux = aof_set_coupler_flux('cfc_11_flux' // suffix,                        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       mol_wt = WTMCFC11, param = (/ 9.36e-07, 9.7561e-06 /),                                   &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)

  cfc(n)%ind_cfc_12_flux = aof_set_coupler_flux('cfc_12_flux' // suffix,                        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       mol_wt = WTMCFC12, param = (/ 9.36e-07, 9.7561e-06 /),                                   &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)
enddo

return

end subroutine  ocmip2_cfc_flux_init
!</SUBROUTINE> NAME="ocmip2_cfc_flux_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>
subroutine ocmip2_cfc_init

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!     Schmidt number coefficients 
!      Use coefficients given by Zheng et al (1998), JGR vol 103, C1
!         for CFC11 and CFC12
real, parameter :: sc_11_0_def = 3501.8
real, parameter :: sc_11_1_def = -210.31
real, parameter :: sc_11_2_def =    6.1851
real, parameter :: sc_11_3_def =   -0.07513

real, parameter :: sc_12_0_def = 3845.4
real, parameter :: sc_12_1_def = -228.95
real, parameter :: sc_12_2_def =    6.1908
real, parameter :: sc_12_3_def =   -0.067430

!     Solubility coefficients for alpha in mol/l/atm
!      (1) for CFC11, (2) for CFC12
!     after Warner and Weiss (1985) DSR, vol 32 for CFC11 and CFC12
real, parameter :: d1_11_def = -229.9261
real, parameter :: d2_11_def =  319.6552
real, parameter :: d3_11_def =  119.4471
real, parameter :: d4_11_def =   -1.39165
real, parameter :: e1_11_def =   -0.142382
real, parameter :: e2_11_def =    0.091459
real, parameter :: e3_11_def =   -0.0157274
 
real, parameter :: d1_12_def = -218.0971
real, parameter :: d2_12_def =  298.9702
real, parameter :: d3_12_def =  113.8049
real, parameter :: d4_12_def =   -1.39165
real, parameter :: e1_12_def =   -0.143566
real, parameter :: e2_12_def =    0.091015
real, parameter :: e3_12_def =   -0.0153924

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_field_name_len+3)                      :: long_suffix
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       Check which tracer packages have been turned on

!       Initialize the ocmip2 cfc package
package_index = otpm_set_tracer_package(package_name,           &
     caller=trim(mod_name) // '(' // trim(sub_name) // ')',     &
     restart_file=default_restart_file )

!       Check whether to use this package
path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif

if (instances .eq. 0) then
  write (stdoutunit,*) trim(note_header), ' No instances'
  do_ocmip2_cfc = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocmip2_cfc = .true.
endif

module_initialized = .true.

!       Return if we don't want to use this package,
!       after changing the list back
if (.not. do_ocmip2_cfc) then
  return
endif

! after reading tracer tree
!       allocate storage for cfc array
allocate ( cfc(instances) )

! loop over the names, saving them into the cfc array
do n = 1, instances
  if (fm_get_value(path_to_names, name, index = n)) then
    cfc(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //                 &
         ' Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
do n = 1, instances

  name = cfc(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

  ! CFC-11
  cfc(n)%ind_cfc_11 = otpm_set_prog_tracer('cfc_11' // suffix, package_name,    &
       longname = 'CFC-11' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                              &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

  ! CFC-12
  cfc(n)%ind_cfc_12 = otpm_set_prog_tracer('cfc_12' // suffix, package_name,    &
       longname = 'CFC-12' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                              &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

enddo

!       Add the package name to the list of good namelists, to be used
!       later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

! Set up the *global* CFC namelist
caller_str=trim(mod_name) // '(' // trim(sub_name) // ')'

! Set up the instance CFC namelists
do n = 1, instances
  call fm_util_start_namelist(package_name, cfc(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('sc_11_0', sc_11_0_def)
  call fm_util_set_value('sc_11_1', sc_11_1_def)
  call fm_util_set_value('sc_11_2', sc_11_2_def)
  call fm_util_set_value('sc_11_3', sc_11_3_def)

  call fm_util_set_value('sc_12_0', sc_12_0_def)
  call fm_util_set_value('sc_12_1', sc_12_1_def)
  call fm_util_set_value('sc_12_2', sc_12_2_def)
  call fm_util_set_value('sc_12_3', sc_12_3_def)

  call fm_util_set_value('d1_11', d1_11_def)
  call fm_util_set_value('d2_11', d2_11_def)
  call fm_util_set_value('d3_11', d3_11_def)
  call fm_util_set_value('d4_11', d4_11_def)

  call fm_util_set_value('d1_12', d1_12_def)
  call fm_util_set_value('d2_12', d2_12_def)
  call fm_util_set_value('d3_12', d3_12_def)
  call fm_util_set_value('d4_12', d4_12_def)

  call fm_util_set_value('e1_11', e1_11_def)
  call fm_util_set_value('e2_11', e2_11_def)
  call fm_util_set_value('e3_11', e3_11_def)

  call fm_util_set_value('e1_12', e1_12_def)
  call fm_util_set_value('e2_12', e2_12_def)
  call fm_util_set_value('e3_12', e3_12_def)

  call fm_util_end_namelist(package_name, cfc(n)%name, check = .true., caller = caller_str)
enddo

!       Check for any errors in the number of fields in the namelists for this package
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

end subroutine ocmip2_cfc_init
! </SUBROUTINE> NAME="ocmip2_cfc_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocmip2_cfc_start
! </DESCRIPTION>
subroutine ocmip2_cfc_init_sfc(isc, iec, jsc, jec, isd, jsd,      &
     isc_bnd, jsc_bnd, Ocean_fields, T_prog, rho, taum1, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer         :: i, j, n
integer :: i_bnd_off
integer :: j_bnd_off
integer         :: ind
real            :: sal
real            :: ta
real            :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
   ! CFC-11 flux
  ind = cfc(n)%ind_cfc_11_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

!     Calculate solubilities
!       Use Warner and Weiss (1985) DSR, vol 32, final result
!       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
!
!       use Bullister and Wisegavger for CCl4
!
!       the factor 1.0e+03 is for the conversion from mol/(l * atm) 
!       to mol/(m3 * atm) 
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      cfc(n)%alpha_11(i,j) =                                                            &
           exp(cfc(n)%d1_11 + cfc(n)%d2_11 / ta + cfc(n)%d3_11 * log(ta) +              &
               cfc(n)%d4_11* ta * ta +                                                  &
               sal * ((cfc(n)%e3_11 * ta + cfc(n)%e2_11) * ta + cfc(n)%e1_11)) *        &
           1.0e+03 * grid_tmask(i,j,1)
    enddo
  enddo

  !     Calculate Schmidt numbers
  !      use coefficients given by Zheng et al (1998), JGR vol 103, C1
    do j = jsc, jec
      do i = isc, iec
        cfc(n)%sc_11(i,j) = cfc(n)%sc_11_0 + t_prog(indtemp)%field(i,j,1,taum1) *       &
             (cfc(n)%sc_11_1 + t_prog(indtemp)%field(i,j,1,taum1) *                     &
              (cfc(n)%sc_11_2 + t_prog(indtemp)%field(i,j,1,taum1) * cfc(n)%sc_11_3)) * &
             grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (cfc(n)%sc_11(i,j) + epsln)) * grid_tmask(i,j,1)
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             cfc(n)%alpha_11(i,j) * sc_no_term(i,j)
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             t_prog(cfc(n)%ind_cfc_11)%field(i,j,1,taum1) * rho(i,j,1,taum1) * sc_no_term(i,j) 
      enddo
    enddo

  endif

  ! CFC-12 flux
  ind = cfc(n)%ind_cfc_12_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then
!     Calculate solubilities
!       Use Warner and Weiss (1985) DSR, vol 32, final result
!       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
!
!       use Bullister and Wisegavger for CCl4
!
!       the factor 1.0e+03 is for the conversion from mol/(l * atm) 
!       to mol/(m3 * atm) 
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      cfc(n)%alpha_12(i,j) =                                                            &
           exp(cfc(n)%d1_12 + cfc(n)%d2_12 / ta + cfc(n)%d3_12 * log(ta) +              &
               cfc(n)%d4_12* ta * ta +                                                  &
               sal * ((cfc(n)%e3_12 * ta + cfc(n)%e2_12) * ta + cfc(n)%e1_12)) *        &
           1.0e+03 * grid_tmask(i,j,1)
    enddo
  enddo

  ! Calculate Schmidt numbers
  ! use coefficients given by Zheng et al (1998), JGR vol 103, C1
    do j = jsc, jec
      do i = isc, iec
        cfc(n)%sc_12(i,j) = cfc(n)%sc_12_0 + t_prog(indtemp)%field(i,j,1,taum1) *       &
             (cfc(n)%sc_12_1 + t_prog(indtemp)%field(i,j,1,taum1) *                     &
              (cfc(n)%sc_12_2 + t_prog(indtemp)%field(i,j,1,taum1) * cfc(n)%sc_12_3)) * &
             grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (cfc(n)%sc_12(i,j) + epsln)) * grid_tmask(i,j,1)
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             cfc(n)%alpha_12(i,j) * sc_no_term(i,j)
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             t_prog(cfc(n)%ind_cfc_12)%field(i,j,1,taum1) * rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo

  endif

enddo

return

end subroutine ocmip2_cfc_init_sfc
! </SUBROUTINE> NAME="ocmip2_cfc_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_cfc_sum_sfc(isc, iec, jsc, jec, isd, jsd,        &
     isc_bnd, jsc_bnd,                                         &
     Ocean_fields, T_prog, rho, taum1, grid_tmask, Grid, Time)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
type(ocean_grid_type), intent(in)                       :: Grid
type(ocean_time_type), intent(in)                       :: Time


integer         :: i, j, n
integer :: i_bnd_off
integer :: j_bnd_off
integer         :: ind
real            :: sal
real            :: ta
real            :: epsln=1.0e-30
logical, save   :: done = .false.
logical, save   :: need = .false.

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances

   ! CFC-11 flux
  ind = cfc(n)%ind_cfc_11_flux

!     Calculate solubilities
!       Use Warner and Weiss (1985) DSR, vol 32, final result
!       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
!
!       use Bullister and Wisegavger for CCl4
!
!       the factor 1.0e+03 is for the conversion from mol/(l * atm) 
!       to mol/(m3 * atm) 
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      cfc(n)%alpha_11(i,j) =                                                            &
           exp(cfc(n)%d1_11 + cfc(n)%d2_11 / ta + cfc(n)%d3_11 * log(ta) +              &
               cfc(n)%d4_11* ta * ta +                                                  &
               sal * ((cfc(n)%e3_11 * ta + cfc(n)%e2_11) * ta + cfc(n)%e1_11)) *        &
           1.0e+03 * grid_tmask(i,j,1)
    enddo
  enddo

  !     Calculate Schmidt numbers
  !      use coefficients given by Zheng et al (1998), JGR vol 103, C1
  do j = jsc, jec
    do i = isc, iec
      cfc(n)%sc_11(i,j) = cfc(n)%sc_11_0 + t_prog(indtemp)%field(i,j,1,taum1) *         &
           (cfc(n)%sc_11_1 + t_prog(indtemp)%field(i,j,1,taum1) *                       &
            (cfc(n)%sc_11_2 + t_prog(indtemp)%field(i,j,1,taum1) * cfc(n)%sc_11_3)) *   &
           grid_tmask(i,j,1)
      sc_no_term(i,j) = sqrt(660.0 / (cfc(n)%sc_11(i,j) + epsln)) * grid_tmask(i,j,1)
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           cfc(n)%alpha_11(i,j) * sc_no_term(i,j)
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           t_prog(cfc(n)%ind_cfc_11)%field(i,j,1,taum1) * rho(i,j,1,taum1) * sc_no_term(i,j)
    enddo
  enddo

  ! CFC-12 flux
  ind = cfc(n)%ind_cfc_12_flux

!     Calculate solubilities
!       Use Warner and Weiss (1985) DSR, vol 32, final result
!       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
!
!       use Bullister and Wisegavger for CCl4
!
!       the factor 1.0e+03 is for the conversion from mol/(l * atm) 
!       to mol/(m3 * atm) 
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      cfc(n)%alpha_12(i,j) =                                                            &
           exp(cfc(n)%d1_12 + cfc(n)%d2_12 / ta + cfc(n)%d3_12 * log(ta) +              &
               cfc(n)%d4_12* ta * ta +                                                  &
               sal * ((cfc(n)%e3_12 * ta + cfc(n)%e2_12) * ta + cfc(n)%e1_12)) *        &
           1.0e+03 * grid_tmask(i,j,1)
    enddo
  enddo

  ! Calculate Schmidt numbers
  ! use coefficients given by Zheng et al (1998), JGR vol 103, C1
  do j = jsc, jec
    do i = isc, iec
      cfc(n)%sc_12(i,j) = cfc(n)%sc_12_0 + t_prog(indtemp)%field(i,j,1,taum1) *         &
           (cfc(n)%sc_12_1 + t_prog(indtemp)%field(i,j,1,taum1) *                       &
            (cfc(n)%sc_12_2 + t_prog(indtemp)%field(i,j,1,taum1) * cfc(n)%sc_12_3)) *   &
           grid_tmask(i,j,1)
      sc_no_term(i,j) = sqrt(660.0 / (cfc(n)%sc_12(i,j) + epsln)) * grid_tmask(i,j,1)
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           cfc(n)%alpha_12(i,j) * sc_no_term(i,j)
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           t_prog(cfc(n)%ind_cfc_12)%field(i,j,1,taum1) * rho(i,j,1,taum1) * sc_no_term(i,j)
    enddo
  enddo

enddo

! Save variables for diagnostics

!       set up the grid mask on the computational grid so that we
!       will not need to implicitly copy arrays in the following
!       subroutine calls
if (.not. done) then
  need = .false.
  do n = 1, instances
    need = need .or.                    &
         cfc(n)%id_alpha_11 .gt. 0 .or. &
         cfc(n)%id_sc_11 .gt. 0 .or.    &
         cfc(n)%id_alpha_12 .gt. 0 .or. &
         cfc(n)%id_sc_12 .gt. 0
  enddo
  done = .true.
endif

if (need) then
  do n = 1, instances
     call diagnose_2d_comp(Time, Grid, cfc(n)%id_alpha_11, cfc(n)%alpha_11(:,:))
     call diagnose_2d_comp(Time, Grid, cfc(n)%id_sc_11, cfc(n)%sc_11(:,:))
     call diagnose_2d_comp(Time, Grid, cfc(n)%id_alpha_12, cfc(n)%alpha_12(:,:))
     call diagnose_2d_comp(Time, Grid, cfc(n)%id_sc_12, cfc(n)%sc_12(:,:))
  enddo
endif

return

end subroutine ocmip2_cfc_sum_sfc
! </SUBROUTINE> NAME="ocmip2_cfc_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_cfc_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

do n = 1, instances
  ind = cfc(n)%ind_cfc_11_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = cfc(n)%ind_cfc_12_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0
enddo

return

end subroutine ocmip2_cfc_zero_sfc
! </SUBROUTINE> NAME="ocmip2_cfc_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_cfc_avg_sfc(isc, iec, jsc, jec, isd, jsd,        &
     isc_bnd, jsc_bnd, Ocean_fields, Ocean_avg_kount, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
integer                                                 :: Ocean_avg_kount
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer :: i_bnd_off
integer :: j_bnd_off
integer :: i, j
integer         :: n
integer         :: ind
real            :: divid

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

divid = 1./float(Ocean_avg_kount)

do n = 1, instances

  ind = cfc(n)%ind_cfc_11_flux

  do j = jsc, jec
    do i = isc, iec
      if (Grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo

  ind = cfc(n)%ind_cfc_12_flux

  do j = jsc, jec
    do i = isc, iec
      if (Grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo

enddo

return

end subroutine ocmip2_cfc_avg_sfc
! </SUBROUTINE> NAME="ocmip2_cfc_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_cfc_sfc_end

end subroutine ocmip2_cfc_sfc_end
! </SUBROUTINE> NAME="ocmip2_cfc_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_source">
!
! <DESCRIPTION>
!     compute the source terms for the CFCs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!
subroutine ocmip2_cfc_source

end subroutine  ocmip2_cfc_source
! </SUBROUTINE> NAME="ocmip2_cfc_source"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!
subroutine ocmip2_cfc_start(isc, iec, jsc, jec, nk, isd, jsd,         &
     T_prog, taup1, model_time, grid_dat, grid_tmask, grid_tracer_axes, rho_dzt)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: jsd
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
integer, intent(in)                                     :: taup1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
integer, dimension(3), intent(in)                       :: grid_tracer_axes
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

character(len=fm_field_name_len+3)      :: long_suffix
integer                                 :: i, j, k, n
character(len=fm_field_name_len+1)      :: suffix
character(len=256)                      :: caller_str
real                                    :: total_cfc_11
real                                    :: total_cfc_12

  integer :: stdoutunit 
  stdoutunit=stdout() 

write(stdoutunit,*) 
write(stdoutunit,*) trim(note_header),                     &
                  ' Starting ', trim(package_name), ' module'

! Determine indices for temperature and salinity
indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif

indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif

! dynamically allocate the global CFC arrays
call allocate_arrays(isc, iec, jsc, jec)

! save the *global* namelist values
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances
  call fm_util_start_namelist(package_name, cfc(n)%name, caller = caller_str)

  cfc(n)%sc_11_0 =    fm_util_get_real   ('sc_11_0', scalar = .true.)
  cfc(n)%sc_11_1 =    fm_util_get_real   ('sc_11_1', scalar = .true.)
  cfc(n)%sc_11_2 =    fm_util_get_real   ('sc_11_2', scalar = .true.)
  cfc(n)%sc_11_3 =    fm_util_get_real   ('sc_11_3', scalar = .true.)
  cfc(n)%sc_12_0 =    fm_util_get_real   ('sc_12_0', scalar = .true.)
  cfc(n)%sc_12_1 =    fm_util_get_real   ('sc_12_1', scalar = .true.)
  cfc(n)%sc_12_2 =    fm_util_get_real   ('sc_12_2', scalar = .true.)
  cfc(n)%sc_12_3 =    fm_util_get_real   ('sc_12_3', scalar = .true.)

  cfc(n)%d1_11 =    fm_util_get_real   ('d1_11', scalar = .true.)
  cfc(n)%d2_11 =    fm_util_get_real   ('d2_11', scalar = .true.)
  cfc(n)%d3_11 =    fm_util_get_real   ('d3_11', scalar = .true.)
  cfc(n)%d4_11 =    fm_util_get_real   ('d4_11', scalar = .true.)
  cfc(n)%d1_12 =    fm_util_get_real   ('d1_12', scalar = .true.)
  cfc(n)%d2_12 =    fm_util_get_real   ('d2_12', scalar = .true.)
  cfc(n)%d3_12 =    fm_util_get_real   ('d3_12', scalar = .true.)
  cfc(n)%d4_12 =    fm_util_get_real   ('d4_12', scalar = .true.)

  cfc(n)%e1_11 =    fm_util_get_real   ('e1_11', scalar = .true.)
  cfc(n)%e2_11 =    fm_util_get_real   ('e2_11', scalar = .true.)
  cfc(n)%e3_11 =    fm_util_get_real   ('e3_11', scalar = .true.)
  cfc(n)%e1_12 =    fm_util_get_real   ('e1_12', scalar = .true.)
  cfc(n)%e2_12 =    fm_util_get_real   ('e2_12', scalar = .true.)
  cfc(n)%e3_12 =    fm_util_get_real   ('e3_12', scalar = .true.)

  call fm_util_end_namelist(package_name, cfc(n)%name, caller = caller_str)
enddo

! Set up analyses

! register the fields
do n = 1, instances

  if (cfc(n)%name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // cfc(n)%name
    long_suffix = ' (' // trim(cfc(n)%name) // ')'
  endif

  cfc(n)%id_sfc_flux_cfc_11 = register_diag_field(trim(diag_name),                      &
       'sfc_flux_cfc_11' // trim(suffix), grid_tracer_axes(1:2),                        &
       model_time, 'Surface Flux - CFC-11' // trim(long_suffix), 'mol m^-2 s^-1',       &
       missing_value = -1.0e+10)

  cfc(n)%id_sfc_flux_cfc_12 = register_diag_field(trim(diag_name),                      &
       'sfc_flux_cfc_12' // trim(suffix), grid_tracer_axes(1:2),                        &
       model_time, 'Surface Flux - CFC-12' // trim(long_suffix), 'mol m^-2 s^-1',       &
       missing_value = -1.0e+10)

  cfc(n)%id_sc_11 = register_diag_field(trim(diag_name),                &
       'sc_11'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Schmidt number - CFC-11'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  cfc(n)%id_alpha_11 = register_diag_field(trim(diag_name),             &
       'alpha_11'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'Solubility CFC-11' // trim(long_suffix), 'mol m^-3 atm^-1',     &
       missing_value = -1.0e+10)

  cfc(n)%id_sc_12 = register_diag_field(trim(diag_name),                &
       'sc_12'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Schmidt number - CFC-12'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  cfc(n)%id_alpha_12 = register_diag_field(trim(diag_name),             &
       'alpha_12'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'Solubility CFC-12' // trim(long_suffix), 'mol m^-3 atm^-1',     &
       missing_value = -1.0e+10)
enddo

!       integrate the total concentrations of some tracers
!       for the start of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
!
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at start of run'

do n = 1, instances

  total_cfc_11 = 0.0
  total_cfc_12 = 0.0

  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_cfc_11 = total_cfc_11 +                           &
             t_prog(cfc(n)%ind_cfc_11)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_cfc_12 = total_cfc_12 +                           &
             t_prog(cfc(n)%ind_cfc_12)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_cfc_11)
  call mpp_sum(total_cfc_12)

  write (stdoutunit,*) '  Instance ', trim(cfc(n)%name)
  write (stdoutunit,                                      &
       '(/'' Total CFC-11  = '',es19.12,'' Gmol'')')    &
       total_cfc_11 * 1.0e-09
  write (stdoutunit,                                      &
       '(/'' Total CFC-12  = '',es19.12,'' Gmol'')')    &
       total_cfc_12 * 1.0e-09

enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), ' Tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocmip2_cfc_start
! </SUBROUTINE> NAME="ocmip2_cfc_start"


end module  ocmip2_cfc_mod
