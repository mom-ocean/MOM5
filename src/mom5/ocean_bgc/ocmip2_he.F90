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
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Jennifer Simeon
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Eric Galbraith
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Daniele Bianchi
!</REVIEWER>
!
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: HE module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 HE
!       simulations as outlined in the Helium-HOWTO documentation,
!       revision 1.6, 1999/04/29.
!       Modified Jan 2008 b1d. Separated atmospheric and mantle component
!       and added a factor to chance the source strength
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Helium/HOWTO-Helium.html
! </REFERENCE>
! </INFO>
!

module  ocmip2_he_mod

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux
use time_manager_mod,         only: get_date, time_type
use time_interp_external_mod, only: time_interp_external, init_external_field
use diag_manager_mod,   only: register_diag_field
use field_manager_mod,  only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,  only: fm_get_length, fm_get_value, fm_new_value, fm_get_index
use fms_mod,            only: field_exist
use fms_io_mod,         only: register_restart_field, save_restart, restore_state
use fms_io_mod,         only: restart_file_type
use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use mpp_mod,            only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use coupler_types_mod,  only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_grid_type, ocean_time_type
use mpp_domains_mod,    only: domain2d
use ocean_util_mod,     only: diagnose_2d, diagnose_2d_comp, diagnose_3d_comp

implicit none

private

public  :: ocmip2_he_bbc
public  :: ocmip2_he_end
public  :: ocmip2_he_init
public  :: ocmip2_he_flux_init
public  :: ocmip2_he_sbc
public  :: ocmip2_he_source
public  :: ocmip2_he_start
public  :: ocmip2_he_tracer
public  :: ocmip2_he_init_sfc
public  :: ocmip2_he_avg_sfc
public  :: ocmip2_he_sum_sfc
public  :: ocmip2_he_zero_sfc
public  :: ocmip2_he_sfc_end
public  :: ocmip2_he_restart

private :: allocate_arrays

character(len=fm_field_name_len), parameter     :: package_name = 'ocmip2_he'
character(len=48), parameter                    :: mod_name = 'ocmip2_he_mod'
character(len=48), parameter              :: diag_name = 'ocmip2_he'
character(len=fm_string_len), parameter   :: default_restart_file =  'ocmip2_he.res.nc'
character(len=fm_string_len), parameter   :: default_ice_restart_file   =    'ice_ocmip2_he.res.nc'
character(len=fm_string_len), parameter   :: default_ocean_restart_file =  'ocmip2_he_airsea_flux.res.nc'

integer, parameter :: max_he_rec = 1200

type he_type

  real                                  :: a1_4
  real                                  :: a2_4
  real                                  :: a3_4
  real                                  :: a4_4
  real                                  :: d1_4
  real                                  :: d2_4
  real                                  :: d3_4
  real                                  :: e1_4
  real                                  :: e2_4
  real                                  :: e3_4
  real                                  :: he4_sourcefac
  real                                  :: he3_sourcefac
  integer                               :: id_sc_3 = -1
  integer                               :: id_alpha_3_atm = -1
  integer                               :: id_alpha_3_man = -1
  integer                               :: id_sc_4 = -1
  integer                               :: id_alpha_4_atm = -1
  integer                               :: id_alpha_4_man = -1
  integer                               :: id_jhe3_man = -1
  integer                               :: id_jhe4_man = -1
  integer                               :: id_jhe_depth = -1
  integer                               :: ind_he_3_atm
  integer                               :: ind_he_3_man
  integer                               :: ind_he_4_atm
  integer                               :: ind_he_4_man
  integer                               :: ind_jhe_3_man
  integer                               :: ind_jhe_4_man
  integer                               :: ind_he_3_atm_flux
  integer                               :: ind_he_3_man_flux
  integer                               :: ind_he_4_atm_flux
  integer                               :: ind_he_4_man_flux
  character(len=fm_field_name_len)      :: name
  character(len=fm_string_len)          :: restart_file
  real, _ALLOCATABLE, dimension(:,:)    :: sc_3  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: sc_4  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha_3_atm  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha_3_man  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha_4_atm  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha_4_man  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jhe3_man _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jhe4_man _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: jhe_depth _NULL
  integer                               :: id_sfc_flux_he_3_atm = -1
  integer                               :: id_sfc_flux_he_4_atm = -1
  integer                               :: id_sfc_flux_he_3_man = -1
  integer                               :: id_sfc_flux_he_4_man = -1

end type he_type

logical, public :: do_ocmip2_he

type(he_type), allocatable, dimension(:)       :: he
integer                                         :: instances
integer                                         :: package_index
logical                                         :: module_initialized = .false.
real, allocatable, dimension(:,:)               :: sc_no_term
integer                                         :: indsal
integer                                         :: indtemp

character(len=128) :: version = '$Id: ocmip2_he.F90,v 20.0 2013/12/14 00:09:46 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

integer                               :: src_he3_id
character*128                         :: src_he3_file    
character*128                         :: src_he3_name
real, allocatable, dimension(:,:)     :: src_he3_t

integer                               :: src_he4_id
character*128                         :: src_he4_file    
character*128                         :: src_he4_name
real, allocatable, dimension(:,:)     :: src_he4_t

integer                               :: src_he_depth_id
character*128                         :: src_he_depth_file    
character*128                         :: src_he_depth_name
real, allocatable, dimension(:,:)     :: src_he_depth_t

! for restart
type(restart_file_type), allocatable :: Top_restart(:)

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
integer, intent(in)     :: nk
integer, intent(in)     :: isd
integer, intent(in)     :: ied
integer, intent(in)     :: jsd
integer, intent(in)     :: jed

integer :: n

allocate( sc_no_term(isc:iec,jsc:jec) )
allocate( src_he3_t(isd:ied,jsd:jed) )
allocate( src_he4_t(isd:ied,jsd:jed) )
allocate( src_he_depth_t(isd:ied,jsd:jed) )

! allocate he array elements
do n = 1, instances
  allocate( he(n)%sc_3(isc:iec,jsc:jec) )
  allocate( he(n)%alpha_3_atm(isc:iec,jsc:jec) )
  allocate( he(n)%alpha_3_man(isc:iec,jsc:jec) )
  allocate( he(n)%jhe3_man(isc:iec,jsc:jec,nk) )
  allocate( he(n)%sc_4(isc:iec,jsc:jec) )
  allocate( he(n)%alpha_4_atm(isc:iec,jsc:jec) )
  allocate( he(n)%alpha_4_man(isc:iec,jsc:jec) )
  allocate( he(n)%jhe4_man(isc:iec,jsc:jec,nk) )
  allocate( he(n)%jhe_depth(isc:iec,jsc:jec) )
enddo

! initialize some arrays
sc_no_term(:,:) = 0.0
src_he3_t(:,:) = 0.0
src_he4_t(:,:) = 0.0
src_he_depth_t(:,:) = 0.0


do n = 1, instances
  he(n)%sc_3(:,:) = 0.0
  he(n)%alpha_3_atm(:,:) = 0.0
  he(n)%alpha_3_man(:,:) = 0.0
  he(n)%jhe3_man(:,:,:) = 0.0
  he(n)%sc_4(:,:) = 0.0
  he(n)%alpha_4_atm(:,:) = 0.0
  he(n)%alpha_4_man(:,:) = 0.0
  he(n)%jhe4_man(:,:,:) = 0.0
  he(n)%jhe_depth(:,:) = 0.0
enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_bbc">
!
! <DESCRIPTION>
!       Called each time-step
!     calculate the bottom boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_he_bbc

end subroutine  ocmip2_he_bbc
! </SUBROUTINE> NAME="ocmip2_he_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_end">
!
! <DESCRIPTION>
!       Called once at the end of the run
!     Clean up various HE quantities for this run.
! </DESCRIPTION>
!
subroutine ocmip2_he_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,   &
     T_prog, grid_dat, grid_tmask, rho_dzt, taup1)

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
real, dimension(isd:ied,jsd:jed), intent(in)            :: grid_dat
real, dimension(isd:ied,jsd:jed,nk), intent(in)         :: grid_tmask
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_he_end'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: i, j, k, n
real    :: total_he_3_atm
real    :: total_he_4_atm
real    :: total_he_3_man
real    :: total_he_4_man

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       integrate the total concentrations of some tracers
!       for the end of the run
total_he_3_atm = 0.0
total_he_3_man = 0.0
total_he_4_atm = 0.0
total_he_4_man = 0.0

!       Use tau time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at end of run'

do n = 1, instances

  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_he_3_atm = total_he_3_atm +                           &
             t_prog(he(n)%ind_he_3_atm)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_he_4_atm = total_he_4_atm +                           &
             t_prog(he(n)%ind_he_4_atm)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_he_3_man = total_he_3_man +                           &
             t_prog(he(n)%ind_he_3_man)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_he_4_man = total_he_4_man +                           &
             t_prog(he(n)%ind_he_4_man)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_he_3_atm)
  call mpp_sum(total_he_4_atm)
  call mpp_sum(total_he_3_man)
  call mpp_sum(total_he_4_man)

  write (stdoutunit,*) '  Instance ', trim(he(n)%name)
  write (stdoutunit,                                              &
       '(/'' Total Atmospheric HE-3  = '',es19.12,'' mol '')')        &
       total_he_3_atm 
  write (stdoutunit,                                              &
       '(/'' Total Atmospheric HE-4  = '',es19.12,'' mol '')')        &
       total_he_4_atm 
  write (stdoutunit,                                              &
       '(/'' Total Mantle HE-3  = '',es19.12,'' mol '')')        &
       total_he_3_man 
  write (stdoutunit,                                              &
       '(/'' Total Mantle HE-4  = '',es19.12,'' mol '')')        &
       total_he_4_man 

enddo

return
end subroutine  ocmip2_he_end
! </SUBROUTINE> NAME="ocmip2_he_end"

!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocmip2_he_restart(time_stamp)
  character(len=*),             intent(in), optional :: time_stamp
  integer :: n

  do n=1, instances
     call save_restart(Top_restart(n), time_stamp)
  end do

end subroutine ocmip2_he_restart
! </SUBROUTINE> NAME="ocmip2_he_restart"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_sbc">
!
! <DESCRIPTION>
!     Called each time-step
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!
subroutine ocmip2_he_sbc(isc, iec, jsc, jec, &
     isc_bnd, jsc_bnd,      &
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
logical :: used

!     use the surface fluxes from the coupler
!       stf is in mol/m^2/s, flux from coupler is positive upwards
i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
       t_prog(he(n)%ind_he_3_atm)%stf(i,j) =                           &
             -ice_ocean_boundary_fluxes%bc(he(n)%ind_he_3_atm_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
       t_prog(he(n)%ind_he_4_atm)%stf(i,j) =                           &
             -ice_ocean_boundary_fluxes%bc(he(n)%ind_he_4_atm_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
       t_prog(he(n)%ind_he_3_man)%stf(i,j) =                           &
             -ice_ocean_boundary_fluxes%bc(he(n)%ind_he_3_man_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
       t_prog(he(n)%ind_he_4_man)%stf(i,j) =                           &
             -ice_ocean_boundary_fluxes%bc(he(n)%ind_he_4_man_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d(Time, Grid, he(n)%id_sfc_flux_he_3_atm, t_prog(he(n)%ind_he_3_atm)%stf(:,:))
   call diagnose_2d(Time, Grid, he(n)%id_sfc_flux_he_4_atm, t_prog(he(n)%ind_he_4_atm)%stf(:,:))
   call diagnose_2d(Time, Grid, he(n)%id_sfc_flux_he_3_man, t_prog(he(n)%ind_he_3_man)%stf(:,:))
   call diagnose_2d(Time, Grid, he(n)%id_sfc_flux_he_4_man, t_prog(he(n)%ind_he_4_man)%stf(:,:))
enddo

return

end subroutine  ocmip2_he_sbc
! </SUBROUTINE> NAME="ocmip2_he_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_flux_init">
!
! <DESCRIPTION>
!     Called once at the beginning of the run
!     Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocmip2_he_flux_init

character(len=64), parameter    :: sub_name = 'ocmip2_he_flux_init'
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
    do_ocmip2_he = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocmip2_he = .true.
  endif

  module_initialized = .true.

endif

! Return if we don't want to use this package
if (.not. do_ocmip2_he) then
  return
endif

if (.not. allocated(he)) then

   ! allocate storage for he array
  allocate ( he(instances) )

  ! loop over the names, saving them into the he array
  do n = 1, instances

    if (fm_get_value(path_to_names, name, index = n)) then
      he(n)%name = name
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

  name = he(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // name
  endif

  ! Coupler fluxes

  ! NOTE: set param = (/ 1, 1 /) if reading in ocmip2 input files
  ! can reset in flux field table
  he(n)%ind_he_3_atm_flux = aof_set_coupler_flux('he_3_atm_flux' // suffix,                        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                    &
       ocean_restart_file = default_ocean_restart_file,                                &
       caller = caller_str)

  he(n)%ind_he_4_atm_flux = aof_set_coupler_flux('he_4_atm_flux' // suffix,                        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                    &
       ocean_restart_file = default_ocean_restart_file,                                &
       caller = caller_str)

  he(n)%ind_he_3_man_flux = aof_set_coupler_flux('he_3_man_flux' // suffix,                        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                    &
       ocean_restart_file = default_ocean_restart_file,                                &
       caller = caller_str)

  he(n)%ind_he_4_man_flux = aof_set_coupler_flux('he_4_man_flux' // suffix,                        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                    &
       ocean_restart_file = default_ocean_restart_file,                                &
       caller = caller_str)
enddo

return

end subroutine  ocmip2_he_flux_init
! </SUBROUTINE> NAME="ocmip2_he_flux_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_init">
!
! <DESCRIPTION>
!       Called once at the beginning of the run
!       Set up any extra fields needed by the tracer packages
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocmip2_he_init

character(len=64), parameter    :: sub_name = 'ocmip2_he_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!     Schmidt number coefficients 
!       Use Wanninkhof (1992) JGR, vol 97, pp 7373-7381
!         for He-4 (and He-3)

real, parameter :: a1_4_def = 410.14
real, parameter :: a2_4_def =  20.503 
real, parameter :: a3_4_def =   0.53175 
real, parameter :: a4_4_def =   0.0060111 

! Solubility coefficients for alpha in mol/(m3 atm)
! (1) for He-4 
! after Wanninkhof (1992) JGR, vol 97, pp 7373-7381
real, parameter :: d1_4_def =  -34.6261
real, parameter :: d2_4_def =   43.0285 
real, parameter :: d3_4_def =   14.1391

real, parameter :: e1_4_def = -0.042340 
real, parameter :: e2_4_def =  0.022624 
real, parameter :: e3_4_def = -0.0033120

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_string_len)                            :: string
character(len=fm_field_name_len+3)                      :: long_suffix
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

  integer :: stdoutunit 
  stdoutunit=stdout() 

! Check which tracer packages have been turned on

! Initialize the ocmip2 he package
package_index = otpm_set_tracer_package(package_name,                 &
     restart_file = default_restart_file,                             &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!       Check whether to use this package
path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then 
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif

if (instances .eq. 0) then
  write (stdoutunit,*) trim(note_header), ' No instances'
  do_ocmip2_he = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocmip2_he = .true.
endif

module_initialized = .true.

! Return if we don't want to use this package,
! after changing the list back
if (.not. do_ocmip2_he) then 
  return
endif

! after reading tracer tree
!       allocate storage for he array
allocate ( he(instances) )

! loop over the names, saving them into the he array
do n = 1, instances

  if (fm_get_value(path_to_names, name, index = n)) then
    he(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //                 &
         ' Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
do n = 1, instances

  name = he(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

  ! NOTE: Coupler wants fluxes in mol/m2/s for MOM4.1. jes.4jun08
  ! HE-3
  he(n)%ind_he_3_atm = otpm_set_prog_tracer('he_3_atm' // suffix, package_name,    &
       longname = 'HE-3 Atmospheric' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                             &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

  he(n)%ind_he_3_man = otpm_set_prog_tracer('he_3_man' // suffix, package_name,    &
       longname = 'HE-3 Mantle' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                             &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

  ! HE-4
  he(n)%ind_he_4_atm = otpm_set_prog_tracer('he_4_atm' // suffix, package_name,    &
       longname = 'HE-4 Atmospheric' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                             &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

  he(n)%ind_he_4_man = otpm_set_prog_tracer('he_4_man' // suffix, package_name,    &
       longname = 'HE-4 Mantle' // trim(long_suffix),                                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                             &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

enddo

! Add the package name to the list of good namelists, to be used
! later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

! Set up the *global* HE namelist
caller_str=trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call fm_util_set_value('src_he_depth_file', 'INPUT/src_he_depth.nc')
call fm_util_set_value('src_he_depth_name', 'V1')

call fm_util_set_value('src_he3_file', 'INPUT/src_he3.nc')
call fm_util_set_value('src_he3_name', 'V2')

call fm_util_set_value('src_he4_file', 'INPUT/src_he4.nc')
call fm_util_set_value('src_he4_name', 'V3')


call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

! Set up the instance HE namelists
do n = 1, instances

  call fm_util_start_namelist(package_name, he(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('a1_4', a1_4_def)
  call fm_util_set_value('a2_4', a2_4_def)
  call fm_util_set_value('a3_4', a3_4_def)
  call fm_util_set_value('a4_4', a4_4_def)

  call fm_util_set_value('d1_4', d1_4_def)
  call fm_util_set_value('d2_4', d2_4_def)
  call fm_util_set_value('d3_4', d3_4_def)

  call fm_util_set_value('e1_4', e1_4_def)
  call fm_util_set_value('e2_4', e2_4_def)
  call fm_util_set_value('e3_4', e3_4_def)

  ! Override the default value, 1, in the field table namelist entries
  call fm_util_set_value('he3_sourcefac', 1.0)
  call fm_util_set_value('he4_sourcefac', 1.0)

  call fm_util_set_value('restart_file', default_restart_file)

  call fm_util_end_namelist(package_name, he(n)%name, check = .true., caller = caller_str)

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

end subroutine ocmip2_he_init
! </SUBROUTINE> NAME="ocmip2_he_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_init_sfc">
!
! <DESCRIPTION>
!       Called once at the beginning of the run
!       Initialize surface fields for flux calculations
!
! </DESCRIPTION>
subroutine ocmip2_he_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,      &
     isc_bnd, jsc_bnd, Ocean_fields, T_prog, rho,            &
     taum1, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
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
real, dimension(isd:ied,jsd:jed,nk), intent(in)         :: grid_tmask

integer         :: i, j, n
integer         :: i_bnd_off
integer         :: j_bnd_off
integer         :: ind
real            :: sal
real            :: ta
real            :: term1, term2, term3

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd


do n = 1, instances
   ! HE flux : ATMOSPHERIC
  ind = he(n)%ind_he_3_atm_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),  &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

!     Calculate solubilities
!       He-4 Sol: Use Wanninkhof (1992) JGR, vol 97, pp 7373-7381
!       He-3 Sol: He-4 Sol * 0.984, See Weiss, 1977; Top et al., 1987;
!                                       Fuchs et al., 1987.
!       Equation for alpha is given in volumetric units (mol/(l atm)).
!       Solubilities output in mol/(m3 * atm) 
!        molar volume = 22414 cm3/mol = 22.4e-3 m3/mol = 22.4e-3 l/mol 
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      he(n)%alpha_4_atm(i,j) =                                                            &
           (exp( he(n)%d1_4 + &
                (he(n)%d2_4/ta) + & 
                (he(n)%d3_4 * (log(ta))) + &
               (sal * ((he(n)%e3_4 * ta * ta) + (he(n)%e2_4 * ta) + he(n)%e1_4))) &
             /(22.4e-3))* grid_tmask(i,j,1)

      he(n)%alpha_3_atm(i,j) = (he(n)%alpha_4_atm(i,j) * 0.984)* grid_tmask(i,j,1)

    enddo
  enddo

  ! Calculate Schmidt numbers
  ! Use Wanninkhof (1992) JGR, vol 97, pp 7373-7381
    do j = jsc, jec
      do i = isc, iec
      term1 = he(n)%a2_4 * t_prog(indtemp)%field(i,j,1,taum1)
      term2 = he(n)%a3_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)
      term3 = he(n)%a4_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)

      he(n)%sc_4(i,j) = he(n)%a1_4 -                                       &
                        term1 + &
                        term2 - &
                        term3 
                       
      he(n)%sc_3(i,j) = he(n)%sc_4(i,j) / 1.15

      sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_3(i,j) + 1.0e-40) * grid_tmask(i,j,1)

    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         he(n)%alpha_3_atm(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         t_prog(he(n)%ind_he_3_atm)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo

  endif

  ! HE-3 flux : MANTLE(Equilibrates with zero atmospheric concentration)
  ind = he(n)%ind_he_3_man_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),  &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then
  ! Calculate solubilities
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      he(n)%alpha_4_man(i,j) =                                                            &
           (exp( he(n)%d1_4 + &
                (he(n)%d2_4/ta) + & 
                (he(n)%d3_4 * (log(ta))) + &
               (sal * ((he(n)%e3_4 * ta * ta) + (he(n)%e2_4 * ta) + he(n)%e1_4))) &
             /(22.4e-3))* grid_tmask(i,j,1)

      he(n)%alpha_3_man(i,j) = (he(n)%alpha_4_man(i,j) * 0.984)* grid_tmask(i,j,1)

    enddo
  enddo

  ! Calculate Schmidt numbers
    do j = jsc, jec
      do i = isc, iec
      term1 = he(n)%a2_4 * t_prog(indtemp)%field(i,j,1,taum1)
      term2 = he(n)%a3_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)
      term3 = he(n)%a4_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)

      he(n)%sc_4(i,j) = he(n)%a1_4 -                                       &
                        term1 + &
                        term2 - &
                        term3 
                       
      he(n)%sc_3(i,j) = he(n)%sc_4(i,j) / 1.15

      sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_3(i,j) + 1.0e-40) * grid_tmask(i,j,1)

    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         he(n)%alpha_3_man(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         t_prog(he(n)%ind_he_3_man)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)

      enddo
    enddo

  endif

  ! HE-4 flux ATMOSPHERIC
  ind = he(n)%ind_he_4_atm_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),  &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then
     ! Solubilities and Schmidt numbers done with He-3 flux (atmos)

     ! Set the He-4 bc values
    do j = jsc, jec
      do i = isc, iec
        sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_4(i,j) + 1.0e-40)* grid_tmask(i,j,1)

    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         he(n)%alpha_4_atm(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         t_prog(he(n)%ind_he_4_atm)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo
  endif

  ! HE-4 flux MANTLE
  ind = he(n)%ind_he_4_man_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),  &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then
     ! Solubilities and Schmidt numbers done with He-3 flux (mantle)

     ! Set the He-4 bc values
    do j = jsc, jec
      do i = isc, iec
        sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_4(i,j) + 1.0e-40)* grid_tmask(i,j,1)
    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         he(n)%alpha_4_man(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         t_prog(he(n)%ind_he_4_man)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)

      enddo
    enddo

  endif

enddo

return

end subroutine ocmip2_he_init_sfc
! </SUBROUTINE> NAME="ocmip2_he_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_sum_sfc">
!
! <DESCRIPTION>
!       Called for FMS coupler
!       ocean_tpm_sum_sfc: Accumulate data for the coupler
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocmip2_he_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,       &
     isc_bnd, jsc_bnd,                                         &
     Ocean_fields, T_prog, rho, taum1, grid_tmask, Grid, Time)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
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
real, dimension(isd:ied,jsd:jed,nk), intent(in)         :: grid_tmask
type(ocean_grid_type), intent(in)                       :: Grid
type(ocean_time_type), intent(in)                       :: Time

integer         :: i, j, n
integer         :: i_bnd_off
integer         :: j_bnd_off
integer         :: ind
real            :: sal
real            :: ta
logical         :: used
real            :: term1, term2, term3

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
   ! HE flux : ATMOSPHERIC
  ind = he(n)%ind_he_3_atm_flux

!     Calculate solubilities
!       He-4 Sol: Use Wanninkhof (1992) JGR, vol 97, pp 7373-7381
!       He-3 Sol: He-4 Sol * 0.984, See Weiss, 1977; Top et al., 1987;
!                                       Fuchs et al., 1987.
!       Equation for alpha is given in volumetric units (mol/(l atm)).
!       Solubilities output in mol/(m3 * atm) 
!        molar volume = 22414 cm3/mol = 22.4e-3 m3/mol = 22.4e-3 l/mol 
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      he(n)%alpha_4_atm(i,j) =                                                            &
           (exp( he(n)%d1_4 + &
                (he(n)%d2_4/ta) + & 
                (he(n)%d3_4 * (log(ta))) + &
               (sal * ((he(n)%e3_4 * ta * ta) + (he(n)%e2_4 * ta) + he(n)%e1_4))) &
             /(22.4e-3))* grid_tmask(i,j,1)

      he(n)%alpha_3_atm(i,j) = (he(n)%alpha_4_atm(i,j) * 0.984)* grid_tmask(i,j,1)
    enddo
  enddo

  ! Calculate Schmidt numbers
  ! Use Wanninkhof (1992) JGR, vol 97, pp 7373-7381
    do j = jsc, jec
      do i = isc, iec
      term1 = he(n)%a2_4 * t_prog(indtemp)%field(i,j,1,taum1)
      term2 = he(n)%a3_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)
      term3 = he(n)%a4_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)

      he(n)%sc_4(i,j) = he(n)%a1_4 -                                       &
                        term1 + &
                        term2 - &
                        term3 
                       
      he(n)%sc_3(i,j) = he(n)%sc_4(i,j) / 1.15

      sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_3(i,j) + 1.0e-40) * grid_tmask(i,j,1)

    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +         &
        he(n)%alpha_3_atm(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +         &
      t_prog(he(n)%ind_he_3_atm)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)

      enddo
    enddo

  ! HE-3 flux : MANTLE(Equilibrates with zero atmospheric concentration)
  ind = he(n)%ind_he_3_man_flux

  ! Calculate solubilities
  do j = jsc, jec
    do i = isc, iec
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      he(n)%alpha_4_man(i,j) =                                                            &
           (exp( he(n)%d1_4 + &
                (he(n)%d2_4/ta) + & 
                (he(n)%d3_4 * (log(ta))) + &
               (sal * ((he(n)%e3_4 * ta * ta) + (he(n)%e2_4 * ta) + he(n)%e1_4))) &
             /(22.4e-3))* grid_tmask(i,j,1)

      he(n)%alpha_3_man(i,j) = (he(n)%alpha_4_man(i,j) * 0.984)* grid_tmask(i,j,1)
    enddo
  enddo

  ! Calculate Schmidt numbers
    do j = jsc, jec
      do i = isc, iec
      term1 = he(n)%a2_4 * t_prog(indtemp)%field(i,j,1,taum1)
      term2 = he(n)%a3_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)
      term3 = he(n)%a4_4 * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)    &
                         * t_prog(indtemp)%field(i,j,1,taum1)

      he(n)%sc_4(i,j) = he(n)%a1_4 -                                       &
                        term1 + &
                        term2 - &
                        term3 
                       
      he(n)%sc_3(i,j) = he(n)%sc_4(i,j) / 1.15

      sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_3(i,j) + 1.0e-40) * grid_tmask(i,j,1)

    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
       Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +        &
         he(n)%alpha_3_man(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
       Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +        &
         t_prog(he(n)%ind_he_3_man)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)

      enddo
    enddo

  ! HE-4 flux ATMOSPHERIC
  ind = he(n)%ind_he_4_atm_flux

  ! Solubilities and Schmidt numbers done with He-3 flux (atmos)

  ! Set the He-4 bc values
    do j = jsc, jec
      do i = isc, iec
        sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_4(i,j) + 1.0e-40)* grid_tmask(i,j,1)

    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
         he(n)%alpha_4_atm(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
         t_prog(he(n)%ind_he_4_atm)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)

      enddo
    enddo

  ! HE-4 flux MANTLE
  ind = he(n)%ind_he_4_man_flux

  ! Solubilities and Schmidt numbers done with He-3 flux (mantle)

  ! Set the He-4 bc values
    do j = jsc, jec
      do i = isc, iec
        sc_no_term(i,j) = sqrt(660.0 / he(n)%sc_4(i,j) + 1.0e-40)* grid_tmask(i,j,1)
    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
         he(n)%alpha_4_man(i,j) * sc_no_term(i,j)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
         t_prog(he(n)%ind_he_4_man)%field(i,j,1,taum1)* rho(i,j,1,taum1) * sc_no_term(i,j)

      enddo
    enddo

enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d_comp(Time, Grid, he(n)%id_alpha_3_atm, he(n)%alpha_3_atm(:,:))
   call diagnose_2d_comp(Time, Grid, he(n)%id_alpha_3_man, he(n)%alpha_3_man(:,:))
   call diagnose_2d_comp(Time, Grid, he(n)%id_sc_3, he(n)%sc_3(:,:))
   call diagnose_2d_comp(Time, Grid, he(n)%id_alpha_4_atm, he(n)%alpha_4_atm(:,:))
   call diagnose_2d_comp(Time, Grid, he(n)%id_alpha_4_man, he(n)%alpha_4_man(:,:))
   call diagnose_2d_comp(Time, Grid, he(n)%id_sc_4, he(n)%sc_4(:,:))
enddo

return

end subroutine ocmip2_he_sum_sfc
! </SUBROUTINE> NAME="ocmip2_he_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_zero_sfc">
!
! <DESCRIPTION>
!               Zero out the fields for the coupler to allow
!               for accumulation for the next time period
! </DESCRIPTION>
subroutine ocmip2_he_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

do n = 1, instances

  ind = he(n)%ind_he_3_atm_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = he(n)%ind_he_4_atm_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = he(n)%ind_he_3_man_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = he(n)%ind_he_4_man_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

enddo

return

end subroutine ocmip2_he_zero_sfc
! </SUBROUTINE> NAME="ocmip2_he_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_avg_sfc">
!
! <DESCRIPTION>
!       Called for FMS coupler
!       ocean_tpm_avg_sfc: Take the time-mean of the fields for the coupler
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocmip2_he_avg_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,       &
     isc_bnd, jsc_bnd, Ocean_fields, Ocean_avg_kount, grid_tmask)

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: jsc_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
integer                                                 :: Ocean_avg_kount
real, dimension(isd:ied,jsd:jed,nk), intent(in)         :: grid_tmask

integer :: i, j, n
integer :: i_bnd_off
integer :: j_bnd_off
integer :: ind
real    :: divid


i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

divid = 1./float(Ocean_avg_kount)

do n = 1, instances

  ind = he(n)%ind_he_3_atm_flux

  do j = jsc, jec
    do i = isc, iec
     if (Grid_tmask(i,j,1) == 1.0) then
       Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =        &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
     endif
    enddo
  enddo

  ind = he(n)%ind_he_4_atm_flux

  do j = jsc, jec
    do i = isc, iec
     if (Grid_tmask(i,j,1) == 1.0) then
    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
     endif
    enddo
  enddo

  ind = he(n)%ind_he_3_man_flux

  do j = jsc, jec
    do i = isc, iec
     if  (Grid_tmask(i,j,1) == 1.0) then
    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
     endif
    enddo
  enddo

  ind = he(n)%ind_he_4_man_flux

  do j = jsc, jec
    do i = isc, iec
     if (Grid_tmask(i,j,1) == 1.0) then
    Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
    Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
     endif
    enddo
  enddo

enddo

return

end subroutine ocmip2_he_avg_sfc
! </SUBROUTINE> NAME="ocmip2_he_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_sfc_end">
!
! <DESCRIPTION>
!       Called for FMS coupler
!       ocean_tpm_sfc_end: Save out fields for the restart.
! </DESCRIPTION>
subroutine ocmip2_he_sfc_end

end subroutine ocmip2_he_sfc_end
! </SUBROUTINE> NAME="ocmip2_he_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_source">
!
! <DESCRIPTION>
!     compute the source terms for the HEs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! As described by J-C Dutay et al.'s  Helium HOWTO.
! Mantle Helium has a source due to emission of helium rich waters
! along ocean ridges on the seafloor.  Globally integrated, the source
! term amounts to 1000 moles of He-3 per year. Regionally, sources are
! partitioned as a function of ridge position, length and spreading rate.
! A loss term exists in the air-sea flux of mantle helium.
!
! The loss term is calculated in the subroutine ocmip2_he_sbc.
! </DESCRIPTION>
!
subroutine ocmip2_he_source(isc, iec, jsc, jec, nk, isd, ied, jsd, jed, t_prog, &
     depth_zt, dzt, model_time, grid_tmask, Grid, Time, grid_kmt)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: t_prog
real, dimension(isd:ied,jsd:jed,nk), intent(in)                 :: depth_zt
real, dimension(isd:ied,jsd:jed,nk), intent(in)                 :: dzt
type(time_type), intent(in)                                     :: model_time
real, dimension(isd:ied,jsd:jed,nk), intent(in)                 :: grid_tmask
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
integer, dimension(isd:ied,jsd:jed), intent(in)                 :: grid_kmt

character(len=64), parameter    :: sub_name = 'ocmip2_he_source'

character(len=256)                      :: caller_str
integer :: index1, index2
integer :: k,i,j,n
logical :: used

call time_interp_external(src_he3_id, model_time, src_he3_t)
call time_interp_external(src_he4_id, model_time, src_he4_t)
call time_interp_external(src_he_depth_id, model_time, src_he_depth_t)

! calculate the source terms for HEs
caller_str=trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

   call fm_util_start_namelist(package_name, he(n)%name, caller = caller_str)
   call fm_util_end_namelist(package_name, he(n)%name, caller = caller_str)

   index1 = he(n)%ind_he_3_man
   index2 = he(n)%ind_he_4_man

   ! these source factors are used to vary the global integral of helium injection
   ! as default use 1.0 * 1.064070463 which with the current source function
   ! gives a global integral flux of 500 mol/year (half the canonical for OCMIP2)

   ! Now in the field table namelist entry
   !   he(n)%he3_sourcefac = 1.0 * 1.064070463 * 0.5
   !   he(n)%he4_sourcefac = 1.0 * 1.064070463 * 0.5
    do j = jsc, jec
      do i = isc, iec
        if (src_he3_t(i,j) .gt. 0 .and. grid_tmask(i,j,1) == 1 ) then
         do k = 1, grid_kmt(i,j)

             if (k .lt. grid_kmt(i,j)) then
               if (depth_zt(i,j,k)-(dzt(i,j,k)/2) .lt. src_he_depth_t(i,j) .and. &
                   depth_zt(i,j,k)+(dzt(i,j,k)/2) .gt. src_he_depth_t(i,j)) then
                 t_prog(index1)%th_tendency(i,j,k) =  &
                                (t_prog(index1)%th_tendency(i,j,k) &
                              + src_he3_t(i,j) * he(n)%he3_sourcefac ) * grid_tmask(i,j,k)
                 t_prog(index2)%th_tendency(i,j,k) =  &
                                (t_prog(index2)%th_tendency(i,j,k) &
                              + src_he4_t(i,j) * he(n)%he4_sourcefac ) * grid_tmask(i,j,k)

                 ! Source term units for history file is mol/m2/s 
                 he(n)%jhe3_man(i,j,k)=(src_he3_t(i,j) * he(n)%he3_sourcefac)/dzt(i,j,k)
                 he(n)%jhe4_man(i,j,k)=(src_he4_t(i,j) * he(n)%he4_sourcefac)/dzt(i,j,k)
                 he(n)%jhe_depth(i,j)=depth_zt(i,j,k)

              else
                 t_prog(index1)%th_tendency(i,j,k) =  &
                                t_prog(index1)%th_tendency(i,j,k)
                 t_prog(index2)%th_tendency(i,j,k) =  &
                                t_prog(index2)%th_tendency(i,j,k)

                 he(n)%jhe3_man(i,j,k)=0
                 he(n)%jhe4_man(i,j,k)=0

              endif

            elseif (k .eq. grid_kmt(i,j)) then
               if (depth_zt(i,j,k)+(dzt(i,j,k)/2) .lt. src_he_depth_t(i,j)) then
                 t_prog(index1)%th_tendency(i,j,k) =  &
                                (t_prog(index1)%th_tendency(i,j,k) &
                              + src_he3_t(i,j) * he(n)%he3_sourcefac) &
                              * grid_tmask(i,j,k)

                 t_prog(index2)%th_tendency(i,j,k) =  &
                                (t_prog(index2)%th_tendency(i,j,k) &
                              + src_he4_t(i,j) * he(n)%he4_sourcefac) &
                              * grid_tmask(i,j,k)
                 he(n)%jhe3_man(i,j,k)=(src_he3_t(i,j) * he(n)%he3_sourcefac)/dzt(i,j,k)
                 he(n)%jhe4_man(i,j,k)=(src_he4_t(i,j) * he(n)%he4_sourcefac)/dzt(i,j,k)
                 he(n)%jhe_depth(i,j)=depth_zt(i,j,k)

               else
                 t_prog(index1)%th_tendency(i,j,k) =  &
                                t_prog(index1)%th_tendency(i,j,k)
                 t_prog(index2)%th_tendency(i,j,k) =  &
                                t_prog(index2)%th_tendency(i,j,k)

                 he(n)%jhe3_man(i,j,k)=0
                 he(n)%jhe4_man(i,j,k)=0

               endif
            endif
         enddo

       else
             t_prog(index1)%th_tendency(i,j,:) = t_prog(index1)%th_tendency(i,j,:) * grid_tmask(i,j,:)
             t_prog(index2)%th_tendency(i,j,:) = t_prog(index2)%th_tendency(i,j,:) * grid_tmask(i,j,:)
             he(n)%jhe3_man(i,j,:)=0
             he(n)%jhe4_man(i,j,:)=0
             he(n)%jhe_depth(i,j)=0

       endif

      enddo
    enddo

  ! Save variables for diagnostics
  call diagnose_3d_comp(Time, Grid, he(n)%id_jhe3_man, he(n)%jhe3_man(:,:,:))
  call diagnose_3d_comp(Time, Grid, he(n)%id_jhe4_man, he(n)%jhe4_man(:,:,:))
  call diagnose_2d_comp(Time, Grid, he(n)%id_jhe_depth, he(n)%jhe_depth(:,:))

enddo

return

end subroutine  ocmip2_he_source
! </SUBROUTINE> NAME="ocmip2_he_source"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!
subroutine ocmip2_he_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,         &
     T_prog, taup1, model_time, grid_dat, grid_tmask, grid_tracer_axes, mpp_domain2d, rho_dzt)

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
real, dimension(isd:ied,jsd:jed), intent(in)            :: grid_dat
real, dimension(isd:ied,jsd:jed,nk), intent(in)         :: grid_tmask
integer, dimension(3), intent(in)                       :: grid_tracer_axes
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_he_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

character(len=fm_field_name_len+3)      :: long_suffix
integer                                 :: i, j, k, n
character(len=fm_field_name_len+1)      :: suffix
character(len=256)                      :: caller_str
real                                    :: total_he_3_atm
real                                    :: total_he_4_atm
real                                    :: total_he_3_man
real                                    :: total_he_4_man

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

! dynamically allocate the global HE arrays
call allocate_arrays(isc, iec, jsc, jec, nk, isd, ied, jsd, jed)

! save the *global* namelist values
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str)

  src_he3_file  =  fm_util_get_string('src_he3_file', scalar = .true.)
  src_he3_name  =  fm_util_get_string('src_he3_name', scalar = .true.)
  src_he4_file  =  fm_util_get_string('src_he4_file', scalar = .true.)
  src_he4_name  =  fm_util_get_string('src_he4_name', scalar = .true.)
  src_he_depth_file  =  fm_util_get_string('src_he_depth_file', scalar = .true.)
  src_he_depth_name  =  fm_util_get_string('src_he_depth_name', scalar = .true.)


call fm_util_end_namelist(package_name, '*global*', caller = caller_str)

do n = 1, instances

  call fm_util_start_namelist(package_name, he(n)%name, caller = caller_str)

  he(n)%a1_4 =    fm_util_get_real('a1_4', scalar = .true.)
  he(n)%a2_4 =    fm_util_get_real('a2_4', scalar = .true.)
  he(n)%a3_4 =    fm_util_get_real('a3_4', scalar = .true.)
  he(n)%a4_4 =    fm_util_get_real('a4_4', scalar = .true.)

  he(n)%d1_4 =    fm_util_get_real('d1_4', scalar = .true.)
  he(n)%d2_4 =    fm_util_get_real('d2_4', scalar = .true.)
  he(n)%d3_4 =    fm_util_get_real('d3_4', scalar = .true.)

  he(n)%e1_4 =    fm_util_get_real('e1_4', scalar = .true.)
  he(n)%e2_4 =    fm_util_get_real('e2_4', scalar = .true.)
  he(n)%e3_4 =    fm_util_get_real('e3_4', scalar = .true.)

  he(n)%he4_sourcefac = fm_util_get_real('he4_sourcefac', scalar = .true.)
  he(n)%he3_sourcefac = fm_util_get_real('he3_sourcefac', scalar = .true.)

  he(n)%restart_file  = fm_util_get_string('restart_file', scalar = .true.)

  call fm_util_end_namelist(package_name, he(n)%name, caller = caller_str)

enddo

! Open up the files for boundary conditions
src_he3_id = init_external_field(src_he3_file,          &
                                     src_he3_name,          &
                                     domain = mpp_domain2d)
if (src_he3_id .eq. 0) then
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open He-3 source file: ' //          &
       trim(src_he3_file))
endif

src_he4_id = init_external_field(src_he4_file,          &
                                     src_he4_name,          &
                                     domain = mpp_domain2d)
if (src_he4_id .eq. 0) then
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open He-3 source file: ' //          &
       trim(src_he4_file))
endif

src_he_depth_id = init_external_field(src_he_depth_file,          &
                                     src_he_depth_name,          &
                                     domain = mpp_domain2d)
if (src_he_depth_id .eq. 0) then
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open He-3 source file: ' //          &
       trim(src_he_depth_file))
endif

! Set up analyses

! register the fields
do n = 1, instances

  if (he(n)%name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // he(n)%name
    long_suffix = ' (' // trim(he(n)%name) // ')'
  endif

  he(n)%id_sfc_flux_he_3_man = register_diag_field('ocean_model',                  &
       'sfc_flux_he_3_man'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Surface Flux - HE-3 man.'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  he(n)%id_sfc_flux_he_3_atm = register_diag_field('ocean_model',                  &
       'sfc_flux_he_3_atm'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Surface Flux - HE-3 atm.'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  he(n)%id_sfc_flux_he_4_man = register_diag_field('ocean_model',                  &
       'sfc_flux_he_4_man'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Surface Flux - HE-4 man.'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  he(n)%id_sfc_flux_he_4_atm = register_diag_field('ocean_model',                  &
       'sfc_flux_he_4_atm'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Surface Flux - HE-4 atm.'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  he(n)%id_sc_3 = register_diag_field('ocean_model',                  &
       'sc_3'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Schmidt number - HE-3'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  he(n)%id_alpha_3_atm = register_diag_field('ocean_model',               &
       'alpha_3_atm'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'Solubility HE-3 Atmospheric'//trim(long_suffix), 'mol m^-3 pptv^-1',      &
       missing_value = -1.0e+10)

  he(n)%id_alpha_3_man = register_diag_field('ocean_model',               &
       'alpha_3_man'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'Solubility HE-3 Mantle'//trim(long_suffix), 'mol m^-3 pptv^-1',      &
       missing_value = -1.0e+10)

  he(n)%id_sc_4 = register_diag_field('ocean_model',                  &
       'sc_4'//trim(suffix), grid_tracer_axes(1:2),                    &
       model_time,                                                      &
       'Schmidt number - HE-4'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  he(n)%id_alpha_4_atm = register_diag_field('ocean_model',               &
       'alpha_4_atm'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'Solubility HE-4 Atmospheric'//trim(long_suffix), 'mol m^-3 pptv^-1',      &
       missing_value = -1.0e+10)

  he(n)%id_alpha_4_man = register_diag_field('ocean_model',               &
       'alpha_4_man'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'Solubility HE-4 Mantle'//trim(long_suffix), 'mol m^-3 pptv^-1',      &
       missing_value = -1.0e+10)

  he(n)%id_jhe3_man = register_diag_field('ocean_model',               &
       'jhe3_man'//trim(suffix), grid_tracer_axes(1:3),                 &
       model_time,                                                      &
       'HE-3 Mantle Source term'//trim(long_suffix), 'mol m^-2 s^-1',      &
       missing_value = -1.0e+10)

  he(n)%id_jhe4_man = register_diag_field('ocean_model',               &
       'jhe4_man'//trim(suffix), grid_tracer_axes(1:3),                 &
       model_time,                                                      &
       'HE-4 Mantle Source term'//trim(long_suffix), 'mol m^-2 s^-1',      &
       missing_value = -1.0e+10)

  he(n)%id_jhe_depth = register_diag_field('ocean_model',               &
       'jhe_depth'//trim(suffix), grid_tracer_axes(1:2),                 &
       model_time,                                                      &
       'HE Injection depth'//trim(long_suffix), 'm',      &
       missing_value = -1.0e+10)

enddo

!       integrate the total concentrations of some tracers
!       for the start of the run
total_he_3_atm = 0.0
total_he_4_atm = 0.0
total_he_3_man = 0.0
total_he_4_man = 0.0

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at start of run'

do n = 1, instances

  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_he_3_atm = total_he_3_atm +                           &
             t_prog(he(n)%ind_he_3_atm)%field(i,j,k,taup1) * &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_he_4_atm = total_he_4_atm +                           &
             t_prog(he(n)%ind_he_4_atm)%field(i,j,k,taup1) *  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_he_3_man = total_he_3_man +                           &
             t_prog(he(n)%ind_he_3_man)%field(i,j,k,taup1) *  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_he_4_man = total_he_4_man +                           &
             t_prog(he(n)%ind_he_4_man)%field(i,j,k,taup1) *  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_he_3_atm)
  call mpp_sum(total_he_4_atm)
  call mpp_sum(total_he_3_man)
  call mpp_sum(total_he_4_man)

  write (stdoutunit,*) '  Instance ', trim(he(n)%name)
  write (stdoutunit,                                              &
       '(/'' Total HE-3 Atmospheric  = '',es19.12,'' mol'')')        &
       total_he_3_atm 
  write (stdoutunit,                                              &
       '(/'' Total HE-4  Atmospheric= '',es19.12,'' mol'')')        &
       total_he_4_atm 
  write (stdoutunit,                                              &
       '(/'' Total HE-3 Mantle = '',es19.12,'' mol'')')        &
       total_he_3_man 
  write (stdoutunit,                                              &
       '(/'' Total HE-4 Mantle = '',es19.12,'' mol'')')        &
       total_he_4_man 

enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), ' Tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocmip2_he_start
! </SUBROUTINE> NAME="ocmip2_he_start"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_he_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!

subroutine ocmip2_he_tracer

end subroutine  ocmip2_he_tracer
! </SUBROUTINE> NAME="ocmip2_he_tracer"

end module  ocmip2_he_mod
