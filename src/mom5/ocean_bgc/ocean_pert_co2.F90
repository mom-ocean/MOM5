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
! Ocean perturbation CO2 module, based on Sarmiento, Orr and Siegenthaler, 1992
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the Ocean perturbation CO2
!       simulations as outlined by "A Perturbation Simulation of
!       CO2 Uptake in an Ocean General Circulation Model", Jorge L. Sarmiento,
!       James C. Orr and Ulrich Siegenthaler, 1992, JGR, 97,
!       pp 3621-3645.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
!       A Perturbation Simulation of
!       CO2 Uptake in an Ocean General Circulation Model, Jorge L. Sarmiento,
!       James C. Orr and Ulrich Siegenthaler, 1992, JGR, 97,
!       pp 3621-3645.
! </REFERENCE>
!
! </INFO>
!
module  ocean_pert_co2_mod

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux
use diag_manager_mod,   only: register_diag_field, diag_axis_init
use field_manager_mod,  only: fm_get_index
use time_manager_mod,   only: time_type
use field_manager_mod,  only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,  only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,            only: field_exist
use mpp_mod,            only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use time_manager_mod,   only: get_date
use mpp_domains_mod,    only: domain2d
use constants_mod,      only: WTMCO2

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use coupler_types_mod,  only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_grid_type, ocean_time_type
use ocmip2_co2calc_mod, only: ocmip2_co2_alpha
use ocean_util_mod,     only: diagnose_2d, diagnose_2d_comp

implicit none

private

public  :: ocean_pert_co2_end
public  :: ocean_pert_co2_init
public  :: ocean_pert_co2_flux_init
public  :: ocean_pert_co2_sbc
public  :: ocean_pert_co2_source
public  :: ocean_pert_co2_start
public  :: ocean_pert_co2_init_sfc
public  :: ocean_pert_co2_avg_sfc
public  :: ocean_pert_co2_sum_sfc
public  :: ocean_pert_co2_zero_sfc

private :: allocate_arrays

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_pert_co2'
character(len=48), parameter                    :: mod_name = 'ocean_pert_co2_mod'
character(len=48), parameter                    :: diag_name = 'ocean_pert_co2'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocean_pert_co2.res.nc'
character(len=fm_string_len), parameter         :: default_ice_restart_file = 'ice_perturbation_co2.res.nc'
character(len=fm_string_len), parameter         :: default_ocean_restart_file = 'ocean_pert_co2_airsea_flux.res.nc'

!  pert_tco2_global           = global annual surface mean perturbation TCO2 concentration
!  pert_tco2_global_wrk       = work variable used in calculation of
!                         pert_tco2_global
!  sal_global           = surface global annual mean salinity
!                         concentration (PSU)
!  sal_global_wrk       = work variable used in calculation of
!                         sal_global
!  do_pert_co2_virtual_flux  = true to compute virtual flux for perturbation TCO2
type instance_type

  real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL
  real                                  :: pert_tco2_global = 0.0
  real                                  :: pert_tco2_global_wrk = 0.0
  logical                               :: do_pert_co2_virtual_flux
  character(len=fm_string_len)          :: implementation
  real                                  :: global_wrk_duration = 0.0
  integer                               :: id_alpha = -1
  integer                               :: id_csurf = -1
  integer                               :: id_pco2surf = -1
  integer                               :: id_z0 = -1
  integer                               :: id_z1 = -1
  integer                               :: id_sfc_flux_pert_co2 = -1
  integer                               :: id_vstf_pert_tco2 = -1
  integer                               :: ind_pert_tco2
  integer                               :: ind_co2_flux
  character(len=fm_field_name_len)      :: name
  real, _ALLOCATABLE, dimension(:,:)    :: pco2surf  _NULL
  real                                  :: sal_global = 35.0
  real                                  :: sal_global_wrk = 0.0
  integer                               :: id_sc_co2 = -1
  real, _ALLOCATABLE, dimension(:,:)    :: sc_co2  _NULL
  real                                  :: sc_co2_0
  real                                  :: sc_co2_1
  real                                  :: sc_co2_2
  real                                  :: sc_co2_3
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_pert_tco2  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: z0  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: z1  _NULL

end type instance_type

logical, public :: do_ocean_pert_co2

integer                                 :: indsal
integer                                 :: indtemp
integer                                 :: package_index
logical                                 :: module_initialized = .false.

character(len=128) :: version = '$Id: ocean_pert_co2.F90,v 20.0 2013/12/14 00:09:34 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!       Calculated parameters (with possible initial input values):
!
!  global_wrk_duration  = total time during calculation of global
!                         variables

real, allocatable, dimension(:,:)               :: sc_no_term
type(instance_type), allocatable, dimension(:)  :: instance
integer                                         :: instances

contains

!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays
! </DESCRIPTION>
!
subroutine allocate_arrays(isc, iec, jsc, jec, nk, isd, ied, jsd, jed)

integer, intent(in)     :: isc
integer, intent(in)     :: iec
integer, intent(in)     :: jsc
integer, intent(in)     :: jec
integer, intent(in)     :: isd
integer, intent(in)     :: ied
integer, intent(in)     :: jsd
integer, intent(in)     :: jed
integer, intent(in)     :: nk

integer :: i, j, n

! global variables
allocate( sc_no_term(isc:iec,jsc:jec) )

! initialize some arrays
sc_no_term(isc:iec,jsc:jec) = 1.0

!       allocate instance array elements
do n = 1, instances
  allocate( instance(n)%csurf(isc:iec,jsc:jec) )
  allocate( instance(n)%alpha(isc:iec,jsc:jec) )
  allocate( instance(n)%pco2surf(isc:iec,jsc:jec) )
  allocate( instance(n)%z0(isc:iec,jsc:jec) )
  allocate( instance(n)%z1(isc:iec,jsc:jec) )
  if (instance(n)%do_pert_co2_virtual_flux) then
    allocate( instance(n)%vstf_pert_tco2(isc:iec,jsc:jec) )
  endif
  allocate( instance(n)%sc_co2(isc:iec,jsc:jec) )
enddo

!       initialize instance array elements
do n = 1, instances

  do j = jsc, jec
    do i = isc, iec
      instance(n)%csurf(i,j) = 0.0
      instance(n)%alpha(i,j) = 0.0
      instance(n)%pco2surf(i,j) = 0.0
      instance(n)%sc_co2(i,j) = 0.0
      instance(n)%z0(i,j) = 0.0
      instance(n)%z1(i,j) = 0.0
    enddo
  enddo
  if (instance(n)%do_pert_co2_virtual_flux) then
    do j = jsc, jec
      do i = isc, iec
        instance(n)%vstf_pert_tco2(i,j) = 0.0
      enddo
    enddo
  endif

enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!
subroutine ocean_pert_co2_bbc

end subroutine  ocean_pert_co2_bbc
! </SUBROUTINE> NAME="ocean_pert_co2_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_end">
!
! <DESCRIPTION>
!     Clean up various quantities for this run.
! </DESCRIPTION>
!
subroutine ocean_pert_co2_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
     T_prog, grid_dat, grid_tmask, mpp_domain2d, rho_dzt, taup1)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
integer, intent(in)                                     :: taup1
real, dimension(isd:,jsd:), intent(in)            :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)         :: grid_tmask
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)       :: rho_dzt

integer                                 :: i
integer                                 :: j
integer                                 :: k
integer                                 :: n
real                                    :: total_pert_tco2

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       integrate the total concentrations of some tracers
!       for the end of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
do n = 1, instances

  total_pert_tco2 = 0.0

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        total_pert_tco2 = total_pert_tco2 +                             &
             t_prog(instance(n)%ind_pert_tco2)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_pert_tco2)

  write (stdoutunit,*) '  Instance ', trim(instance(n)%name)
  write (stdoutunit,                                                      &
       '(/'' Total pert TCO2  = '',es19.12,'' Gmol-C'')')               &
       total_pert_tco2 * 1.0e-09
enddo

write(stdoutunit,*)

return
end subroutine  ocean_pert_co2_end
! </SUBROUTINE> NAME="ocean_pert_co2_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_pert_co2_sbc(isc, iec, jsc, jec,       &
     isc_bnd, jsc_bnd,                      &
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

integer :: i, j, k, n, m
integer :: i_bnd_off
integer :: j_bnd_off
integer :: kz
logical :: used

!     use the surface fluxes from the coupler
!       stf is in mol/m^2/s, flux from coupler is positive upwards
i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
      t_prog(instance(n)%ind_pert_tco2)%stf(i,j) =                      &
            -ice_ocean_boundary_fluxes%bc(instance(n)%ind_co2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo

!     add in the virtual fluxes as defined by equations (2) and (3)
!     in the OCMIP2 ABIOTIC HOWTO.
!       Note: the factor of 1000 is to convert the delta salinity from
!             model units to PSU

!       Save variables for diagnostics
do n = 1, instances
   call diagnose_2d_comp(Time, Grid, instance(n)%id_sc_co2, instance(n)%sc_co2(:,:))
   call diagnose_2d(Time, Grid, instance(n)%id_sfc_flux_pert_co2, t_prog(instance(n)%ind_pert_tco2)%stf(:,:))
enddo

return

end subroutine  ocean_pert_co2_sbc
! </SUBROUTINE> NAME="ocean_pert_co2_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocean_pert_co2_flux_init

character(len=64), parameter    :: sub_name = 'ocean_pert_co2_flux_init'
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
       restart_file = default_restart_file,     &
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
    do_ocean_pert_co2 = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocean_pert_co2 = .true.
  endif

  module_initialized = .true.

endif

! Return if we don't want to use this package
if (.not. do_ocean_pert_co2) then
  return
endif

if (.not. allocated(instance)) then

   ! allocate storage for instance array
  allocate ( instance(instances) )

  ! loop over the names, saving them into the instance array
  do n = 1, instances

    if (fm_get_value(path_to_names, name, index = n)) then
      instance(n)%name = name
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

  name = instance(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // name
  endif

  ! Coupler fluxes
  instance(n)%ind_co2_flux = aof_set_coupler_flux('pert_co2_flux' // suffix,                    &
       flux_type = 'air_sea_gas_flux', implementation = 'linear',                               &
       mol_wt = WTMCO2, param = (/ 4.033e-10, 2.0, 1.0 /),                                                       &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)

  instance(n)%implementation = fm_util_get_string ('/coupler_mod/fluxes/pert_co2_flux' // trim(suffix) //       &
       '/implementation', scalar = .true.)
  if (.not. (instance(n)%implementation .eq. 'linear' .or.              &
             instance(n)%implementation .eq. 'ocmip2' .or.              &
             instance(n)%implementation .eq. 'ocmip2_data')) then
    call mpp_error(FATAL, trim(error_header) //                         &
         'Unsupported flux implementation: "' // trim(instance(n)%implementation) // '"')
  endif
enddo

return

end subroutine  ocean_pert_co2_flux_init
! </SUBROUTINE> NAME="ocean_pert_co2_flux_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocean_pert_co2_init

character(len=64), parameter    :: sub_name = 'ocean_pert_co2_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_field_name_len+3)                      :: long_suffix
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! Initialize the package
package_index = otpm_set_tracer_package(package_name,            &
     restart_file = default_restart_file,   &
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
  do_ocean_pert_co2 = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocean_pert_co2 = .true.
endif

module_initialized = .true.

! Return if we don't want to use this package
if (.not. do_ocean_pert_co2) then
  return
endif

! allocate storage for instance array
allocate ( instance(instances) )

! loop over the names, saving them into the instance array
do n = 1, instances

  if (fm_get_value(path_to_names, name, index = n)) then
    instance(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = instance(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

  ! perturbation TCO2
  instance(n)%ind_pert_tco2 = otpm_set_prog_tracer('pert_tco2' // suffix,       &
       package_name,                                                            &
       longname = 'perturbation TCO2' // trim(long_suffix),                     &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                             &
       caller = caller_str)

enddo

!       Process the namelists

!       Add the package name to the list of good namelists, to be used
!       later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

! Set up the instance namelists
do n = 1, instances
   ! create the instance namelist
  call fm_util_start_namelist(package_name, instance(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('sal_global', 35.0)                            ! PSU
  call fm_util_set_value('do_pert_co2_virtual_flux', .false.)
  call fm_util_set_value('pert_tco2_global', 2.0)                       ! mol/m^3

  ! New Wanninkhof numbers
  call fm_util_set_value('sc_co2_0', 2068.9)
  call fm_util_set_value('sc_co2_1', -118.63)
  call fm_util_set_value('sc_co2_2', 2.9311)
  call fm_util_set_value('sc_co2_3', -0.027)

  call fm_util_end_namelist(package_name, instance(n)%name, check = .true., caller = caller_str)

enddo

! Check for any errors in the number of fields in the namelists for this package
good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',     &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then
  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,                   &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif

return

end subroutine ocean_pert_co2_init
! </SUBROUTINE> NAME="ocean_pert_co2_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocean_pert_co2_start
! </DESCRIPTION>
subroutine ocean_pert_co2_init_sfc(isc, iec, jsc, jec, isd, ied, jsd, jed,  &
     isc_bnd, jsc_bnd,                                         &
     Ocean_fields, T_prog, rho, taum1, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer :: i, j, m, n
integer :: i_bnd_off
integer :: j_bnd_off
integer :: nn
integer :: ind

real    :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances

   ! CO2 flux
  ind = instance(n)%ind_co2_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then
     !  Compute the moist air compensation term (placed in alpha) and the
     !  surface delta ocean perturbation CO2

     ! z0 has units of ppm/(umol/kg) and z1 has units of 1/(umol/kg),
     ! and we need to convert then to (kg/kg)/(mol/m^3) and 1/(mol/m^3), respectively
     
     ! Changed 0.31618 to 0.031618, as it appears to have been an error in the paper
    if (instance(n)%implementation .eq. 'linear') then
      do j = jsc, jec
        do i = isc, iec
          instance(n)%alpha(i,j) = (1.0 - exp(20.1050 -                                           &
               0.0097982 * (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) -                        &
               6163.10 / (t_prog(indtemp)%field(i,j,1,taum1) + 273.15))) * 1.0e+06 * grid_tmask(i,j,1)
        enddo
      enddo
    else
      call ocmip2_co2_alpha(                                                    &
           isd, jsd, isc, iec, jsc, jec,                              &
           t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                      &
           t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1), grid_tmask(isd:ied,jsd:jed,1), instance(n)%alpha)
      do j = jsc, jec
        do i = isc, iec
          instance(n)%sc_co2(i,j) =                                             &
             instance(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *        &
             (instance(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *       &
              (instance(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *      &
               instance(n)%sc_co2_3)) * grid_tmask(i,j,1)
          sc_no_term(i,j) = sqrt(660.0 / (instance(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)
        enddo
      enddo
    endif
    do j = jsc, jec
      do i = isc, iec
        instance(n)%z0(i,j) = (1.7561 -                                                         &
             0.031618 * t_prog(indtemp)%field(i,j,1,taum1) +                                    &
             0.0004444 * t_prog(indtemp)%field(i,j,1,taum1)**2) * grid_tmask(i,j,1)
        instance(n)%z1(i,j) = (0.004096 -                                                       &
             7.7086e-05 * t_prog(indtemp)%field(i,j,1,taum1) +                                  &
             6.10e-07 * t_prog(indtemp)%field(i,j,1,taum1)**2) * grid_tmask(i,j,1)
        instance(n)%pco2surf(i,j) = (instance(n)%z0(i,j) *                                      &
             t_prog(instance(n)%ind_pert_tco2)%field(i,j,1,taum1) * 1.0e+06 / 1024.5 /          &
             (1.0 - instance(n)%z1(i,j) *                                                       &
              t_prog(instance(n)%ind_pert_tco2)%field(i,j,1,taum1) * 1.0e+06 / 1024.5)) * grid_tmask(i,j,1)
        instance(n)%csurf(i,j) = instance(n)%alpha(i,j) * instance(n)%pco2surf(i,j) * 1.0e-06
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             instance(n)%alpha(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             instance(n)%csurf(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo

  endif

enddo

return

end subroutine ocean_pert_co2_init_sfc
! </SUBROUTINE> NAME="ocean_pert_co2_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_pert_co2_sum_sfc(isc, iec, jsc, jec, isd, ied, jsd, jed,   &
     isc_bnd, jsc_bnd,                            &
     Ocean_fields, T_prog, rho, taum1, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer :: i, j, n, nn
integer :: i_bnd_off
integer :: j_bnd_off
integer :: ind

real    :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances

    ind = instance(n)%ind_co2_flux

    !  Compute the moist air compensation term (placed in alpha) and the
    !  surface delta ocean perturbation CO2

    !       z0 has units of ppm/(umol/kg) and z1 has units of 1/(umol/kg),
    !       and we need to convert then to (kg/kg)/(mol/m^3) and 1/(mol/m^3), respectively

    !       Changed 0.31618 to 0.031618, as it appears to have been an error in the paper
    if (instance(n)%implementation .eq. 'linear') then
      do j = jsc, jec
        do i = isc, iec
          instance(n)%alpha(i,j) = (1.0 - exp(20.1050 -                                           &
               0.0097982 * (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) -                        &
               6163.10 / (t_prog(indtemp)%field(i,j,1,taum1) + 273.15))) * 1.0e+06 * grid_tmask(i,j,1)
        enddo
      enddo
    else
      call ocmip2_co2_alpha(                                                    &
           isd, jsd, isc, iec, jsc, jec,                              &
           t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                      &
           t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1), grid_tmask(isd:ied,jsd:jed,1), instance(n)%alpha)
      do j = jsc, jec
        do i = isc, iec
          instance(n)%sc_co2(i,j) =                                             &
             instance(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *        &
             (instance(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *       &
              (instance(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *      &
               instance(n)%sc_co2_3)) * grid_tmask(i,j,1)
          sc_no_term(i,j) = sqrt(660.0 / (instance(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)
        enddo
      enddo
    endif
    do j = jsc, jec
      do i = isc, iec
        instance(n)%z0(i,j) = (1.7561 -                                                         &
             0.031618 * t_prog(indtemp)%field(i,j,1,taum1) +                                    &
             0.0004444 * t_prog(indtemp)%field(i,j,1,taum1)**2) * grid_tmask(i,j,1)
        instance(n)%z1(i,j) = (0.004096 -                                                       &
             7.7086e-05 * t_prog(indtemp)%field(i,j,1,taum1) +                                  &
             6.10e-07 * t_prog(indtemp)%field(i,j,1,taum1)**2) * grid_tmask(i,j,1)
        instance(n)%pco2surf(i,j) = (instance(n)%z0(i,j) *                                      &
             t_prog(instance(n)%ind_pert_tco2)%field(i,j,1,taum1) * 1.0e+06 / 1024.5 /          &
             (1.0 - instance(n)%z1(i,j) *                                                       &
              t_prog(instance(n)%ind_pert_tco2)%field(i,j,1,taum1) * 1.0e+06 / 1024.5)) * grid_tmask(i,j,1)
        instance(n)%csurf(i,j) = instance(n)%alpha(i,j) * instance(n)%pco2surf(i,j) * 1.0e-06
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +            &
             instance(n)%alpha(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +            &
             instance(n)%csurf(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo
enddo

return

end subroutine ocean_pert_co2_sum_sfc
! </SUBROUTINE> NAME="ocean_pert_co2_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocean_pert_co2_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

do n = 1, instances
  ind = instance(n)%ind_co2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0
enddo

return

end subroutine ocean_pert_co2_zero_sfc
! </SUBROUTINE> NAME="ocean_pert_co2_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_pert_co2_avg_sfc(isc, iec, jsc, jec, isd, jsd,     &
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
integer :: i, j, n
integer :: ind
real    :: divid

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

divid = 1./float(Ocean_avg_kount)

do n = 1, instances

  ind = instance(n)%ind_co2_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo
enddo

return

end subroutine ocean_pert_co2_avg_sfc
! </SUBROUTINE> NAME="ocean_pert_co2_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>
subroutine ocean_pert_co2_sfc_end 

end subroutine ocean_pert_co2_sfc_end 
! </SUBROUTINE> NAME="ocean_pert_co2_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_source">
!
! <DESCRIPTION>
!     compute the source terms, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!
subroutine ocean_pert_co2_source(Grid, Time)

type(ocean_grid_type), intent(in) :: Grid
type(ocean_time_type), intent(in) :: Time

integer :: n

do n = 1, instances
   call diagnose_2d_comp(Time, Grid, instance(n)%id_alpha, instance(n)%alpha(:,:))
   call diagnose_2d_comp(Time, Grid, instance(n)%id_csurf, instance(n)%csurf(:,:))
   call diagnose_2d_comp(Time, Grid, instance(n)%id_pco2surf, instance(n)%pco2surf(:,:))
   call diagnose_2d_comp(Time, Grid, instance(n)%id_z0, instance(n)%z0(:,:))
   call diagnose_2d_comp(Time, Grid, instance(n)%id_z1, instance(n)%z1(:,:))
enddo

return

end subroutine  ocean_pert_co2_source
! </SUBROUTINE> NAME="ocean_pert_co2_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_pert_co2_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!
subroutine ocean_pert_co2_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, taup1, model_time, grid_dat, grid_tmask,                           &
     grid_tracer_axes, rho_dzt)

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
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
integer, dimension(:), intent(in)                       :: grid_tracer_axes
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocean_pert_co2_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                 :: i, j, k, n
character(len=fm_field_name_len+1)      :: suffix
character(len=fm_field_name_len+3)      :: long_suffix
character(len=256)                      :: caller_str
real                                    :: total_pert_tco2

  integer :: stdoutunit 
  stdoutunit=stdout() 

write(stdoutunit,*) 
write(stdoutunit,*) trim(note_header),                     &
                  'Starting ', trim(package_name), ' module'

!       Determine indices for temperature and salinity
indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif

indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif

! dynamically allocate the global arrays
call allocate_arrays(isc, iec, jsc, jec, nk, isd, ied, jsd, jed)

! save the *global* namelist values
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

! read in the namelists for each instance
do n = 1, instances
  call fm_util_start_namelist(package_name, instance(n)%name, caller = caller_str)

  instance(n)%sal_global                = fm_util_get_real   ('sal_global', scalar = .true.)
  instance(n)%do_pert_co2_virtual_flux  = fm_util_get_logical('do_pert_co2_virtual_flux', scalar = .true.)
  instance(n)%pert_tco2_global          = fm_util_get_real   ('pert_tco2_global', scalar = .true.)
  instance(n)%sc_co2_0                  = fm_util_get_real   ('sc_co2_0', scalar = .true.)
  instance(n)%sc_co2_1                  = fm_util_get_real   ('sc_co2_1', scalar = .true.)
  instance(n)%sc_co2_2                  = fm_util_get_real   ('sc_co2_2', scalar = .true.)
  instance(n)%sc_co2_3                  = fm_util_get_real   ('sc_co2_3', scalar = .true.)

  call fm_util_end_namelist(package_name, instance(n)%name, caller = caller_str)
enddo

!     Set up analyses

!       register the instance fields
do n = 1, instances

  if (instance(n)%name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // instance(n)%name
    long_suffix = ' (' // trim(instance(n)%name) // ')'
  endif

  instance(n)%id_sc_co2 = register_diag_field(trim(diag_name),          &
       'sc_co2_' // trim(suffix), grid_tracer_axes(1:2),                &
       model_time, 'Schmidt number - CO2' // trim(long_suffix), ' ',    &
       missing_value = -1.0e+10)
  instance(n)%id_alpha = register_diag_field(trim(diag_name),           &
       'alpha' // trim(suffix), grid_tracer_axes(1:2),                  &
       model_time, 'Alpha CO2' // trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  instance(n)%id_csurf = register_diag_field(trim(diag_name),           &
       'csurf' // trim(suffix), grid_tracer_axes(1:2),                  &
       model_time, 'CO2* water' // trim(long_suffix), ' ',              &
       missing_value = -1.0e+10)

  instance(n)%id_pco2surf = register_diag_field(trim(diag_name),        &
       'pco2surf' // trim(suffix), grid_tracer_axes(1:2),               &
       model_time, 'Oceanic pCO2' // trim(long_suffix), 'ppm',          &
       missing_value = -1.0e+10)

  instance(n)%id_z0 = register_diag_field(trim(diag_name),              &
       'z0' // trim(suffix), grid_tracer_axes(1:2),                     &
       model_time, 'z0' // trim(long_suffix), ' ',                      &
       missing_value = -1.0e+10)

  instance(n)%id_z1 = register_diag_field(trim(diag_name),              &
       'z1' // trim(suffix), grid_tracer_axes(1:2),                     &
       model_time, 'z1' // trim(long_suffix), ' ',                      &
       missing_value = -1.0e+10)

  instance(n)%id_sfc_flux_pert_co2 = register_diag_field(trim(diag_name),       &
       'sfc_flux_pert_co2' // trim(suffix), grid_tracer_axes(1:2),              &
       model_time, 'CO2 surface flux' // trim(long_suffix), 'mol m^-1 s^-1',    &
       missing_value = -1.0e+10)
enddo

!       integrate the total concentrations of some tracers
!       for the start of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
do n = 1, instances

  total_pert_tco2 = 0.0

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        total_pert_tco2 = total_pert_tco2 +                             &
             t_prog(instance(n)%ind_pert_tco2)%field(i,j,k,taup1) *       &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_pert_tco2)

  write (stdoutunit,*) '  Instance ', trim(instance(n)%name)
  write (stdoutunit,                                                      &
       '(/'' Total pert TCO2  = '',es19.12,'' Gmol-C'')')               &
       total_pert_tco2 * 1.0e-09
enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), 'ocean_pert CO2 tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocean_pert_co2_start
! </SUBROUTINE> NAME="ocean_pert_co2_start"

end module  ocean_pert_co2_mod
