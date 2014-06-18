#include <fms_platform.h>
module  ocmip2_abiotic_mod
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: Abiotic module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 Abiotic
!       simulations as outlined in the Abiotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Abiotic/HOWTO-Abiotic.html
! </REFERENCE>
!
! <REFERENCE>
! Press, W. H., S. A. Teukosky, W. T. Vetterling, B. P. Flannery, 1992. 
! Numerical Recipes in FORTRAN, Second Edition, Cambridge University Press. 
! </REFERENCE>
!
! <REFERENCE>
! Enting, I.G., T. M. L. Wigley, M. Heimann, 1994. Future Emissions 
! and concentrations of carbon dioxide: key ocean / atmosphere / 
! land analyses, CSIRO Aust. Div. Atmos. Res. Tech. Pap. No. 31, 
! 118 pp.
! </REFERENCE>
! </INFO>
!
!------------------------------------------------------------------
!
!       Module ocmip2_abiotic_mod
!
!       Implementation of routines to solve the OCMIP-2 Abiotic
!       simulations as outlined in the Abiotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!
!------------------------------------------------------------------
!
use atmos_ocean_fluxes_mod,   only: aof_set_coupler_flux
use diag_manager_mod,         only: register_diag_field, diag_axis_init
use field_manager_mod,        only: fm_get_index
use time_manager_mod,         only: time_type
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,                  only: field_exist, file_exist
use fms_io_mod,               only: register_restart_field, save_restart, restore_state
use fms_io_mod,               only: restart_file_type
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL, NOTE
use time_manager_mod,         only: get_date
use time_interp_external_mod, only: time_interp_external
use time_interp_external_mod, only: init_external_field
use mpp_domains_mod,          only: domain2d
use mpp_domains_mod,          only: mpp_global_sum, BITWISE_EXACT_SUM
use constants_mod,            only: WTMCO2

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use coupler_types_mod,  only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use ocean_types_mod,    only: ocean_prog_tracer_type
use ocmip2_co2calc_mod, only: ocmip2_co2calc
use ocean_types_mod,    only: ocean_grid_type, ocean_time_type
use ocean_util_mod,     only: diagnose_2d, diagnose_2d_comp, diagnose_3d_comp

implicit none

private

public  :: ocmip2_abiotic_end
public  :: ocmip2_abiotic_init
public  :: ocmip2_abiotic_flux_init
public  :: ocmip2_abiotic_sbc
public  :: ocmip2_abiotic_source
public  :: ocmip2_abiotic_start
public  :: ocmip2_abiotic_init_sfc
public  :: ocmip2_abiotic_avg_sfc
public  :: ocmip2_abiotic_sum_sfc
public  :: ocmip2_abiotic_zero_sfc
public  :: ocmip2_abiotic_sfc_end
public  :: ocmip2_abiotic_restart
public  :: ocmip2_abiotic_tracer

private :: allocate_arrays

character(len=fm_field_name_len), parameter     :: package_name = 'ocmip2_abiotic'
character(len=48), parameter                    :: mod_name = 'ocmip2_abiotic_mod'
character(len=48), parameter                    :: diag_name = 'ocean_ocmip2_abiotic'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocmip2_abiotic.res.nc'
character(len=fm_string_len), parameter         :: default_local_restart_file = 'ocmip2_abiotic_local.res.nc'
character(len=fm_string_len), parameter         :: default_ice_restart_file = 'ice_ocmip2_abiotic.res.nc'
character(len=fm_string_len), parameter         :: default_ocean_restart_file = 'ocmip2_abiotic_airsea_flux.res.nc'

!  sio4_const           = SiO4 concentration (mol/kg)
!  po4_const            = PO4 concentration (mol/kg)
!  dic_global           = global annual surface mean DIC concentration
!  dic_global_wrk       = work variable used in calculation of
!                         dic_global
!  di14c_global         = global annual surface mean DI14C concentration
!  di14c_global_wrk     = work variable used in calculation of
!                         DI14C_global
!  sal_global           = surface global annual mean salinity
!                         concentration (PSU)
!  sal_global_wrk       = work variable used in calculation of
!                         sal_global
type abiotic_type

  real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: c14surf  _NULL
  real                                  :: di14c_global = 0.0
  real, _ALLOCATABLE, dimension(:,:)    :: di14c_global_wrk  _NULL
  real                                  :: dic_global = 0.0
  real, _ALLOCATABLE, dimension(:,:)    :: dic_global_wrk  _NULL
  character(len=128)                    :: frac_14catm_file = ' '
  character(len=128)                    :: frac_14catm_name = ' '
  integer                               :: frac_14catm_id = 0
  real, _ALLOCATABLE, dimension(:,:)    :: frac_14catm  _NULL
  real                                  :: frac_14catm_const = 0.0
  character(len=fm_string_len)          :: local_restart_file
  real                                  :: global_wrk_duration = 0.0
  real, _ALLOCATABLE, dimension(:,:)    :: htotal  _NULL
  integer                               :: id_alpha = -1
  integer                               :: id_csurf = -1
  integer                               :: id_c14surf = -1
  integer                               :: id_htotal = -1
  integer                               :: id_alk = -1
  integer                               :: id_po4 = -1
  integer                               :: id_sio4 = -1
  integer                               :: id_frac_14catm = -1
  integer                               :: id_jdi14c = -1
  integer                               :: id_pco2surf = -1
  integer                               :: id_sfc_flux_co2 = -1
  integer                               :: id_sfc_flux_14co2 = -1
  integer                               :: id_vstf_di14c = -1
  integer                               :: id_vstf_dic = -1
  integer                               :: ind_di14c
  integer                               :: ind_dic
  integer                               :: ind_co2_flux
  integer                               :: ind_14co2_flux
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdi14c  _NULL
  real                                  :: half_life
  real                                  :: lambda
  real                                  :: alkbar
  integer                               :: id_sc_co2 = -1
  real, _ALLOCATABLE, dimension(:,:)    :: sc_co2  _NULL
  real                                  :: sc_co2_0
  real                                  :: sc_co2_1
  real                                  :: sc_co2_2
  real                                  :: sc_co2_3
  character(len=fm_field_name_len)      :: name
  real, _ALLOCATABLE, dimension(:,:)    :: pco2surf  _NULL
  real                                  :: sal_global = 35.0
  real, _ALLOCATABLE, dimension(:,:)    :: sal_global_wrk  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: po4  _NULL
  real                                  :: po4_const
  real, _ALLOCATABLE, dimension(:,:)    :: sio4  _NULL
  real                                  :: sio4_const
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_di14c  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_dic  _NULL

end type abiotic_type

logical, public :: do_ocmip2_abiotic

integer                 :: indsal
integer                 :: indtemp
integer                 :: package_index
logical                 :: module_initialized = .false.

character(len=128)      :: version = '$Id: ocmip2_abiotic.F90,v 20.0 2013/12/14 00:09:38 fms Exp $'
character(len=128)      :: tagname = '$Name: tikal $'

!       Input parameters:
!
!  htotal_in            = default value for htotal for an initial run
!  htotal_scale_lo      = scaling parameter to chose htotallo
!  htotal_scale_hi      = scaling parameter to chose htotalhi
real                                    :: htotal_in
real, allocatable, dimension(:,:)       :: htotal_scale_hi
real                                    :: htotal_scale_hi_in
real, allocatable, dimension(:,:)       :: htotal_scale_lo
real                                    :: htotal_scale_lo_in

!       Calculated parameters (with possible initial input values):
!
!  global_wrk_duration  = total time during calculation of global
!                         variables
real, allocatable, dimension(:,:)               :: sc_no_term
type(abiotic_type), allocatable, dimension(:)   :: abiotic
integer                                         :: instances

! for restart
integer                              :: num_restart = 0
type(restart_file_type), allocatable :: restart(:)

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

integer :: i, j, k, l, m, n

! global variables
allocate( sc_no_term(isc:iec,jsc:jec) )
allocate( htotal_scale_lo(isc:iec,jsc:jec) )
allocate( htotal_scale_hi(isc:iec,jsc:jec) )

! initialize some arrays
sc_no_term(:,:) = 0.0
htotal_scale_lo(:,:) = 0.0
htotal_scale_hi(:,:) = 0.0

! allocate abiotic array elements
do n = 1, instances
  allocate( abiotic(n)%htotal(isc:iec,jsc:jec) )
  allocate( abiotic(n)%csurf(isc:iec,jsc:jec) )
  allocate( abiotic(n)%c14surf(isc:iec,jsc:jec) )
  allocate( abiotic(n)%frac_14catm(isc:iec,jsc:jec) )
  allocate( abiotic(n)%alpha(isc:iec,jsc:jec) )
  allocate( abiotic(n)%pco2surf(isc:iec,jsc:jec) )
  allocate( abiotic(n)%po4(isd:ied,jsd:jed) )
  allocate( abiotic(n)%sio4(isd:ied,jsd:jed) )
  allocate( abiotic(n)%sc_co2(isc:iec,jsc:jec) )
  allocate( abiotic(n)%sal_global_wrk(isc:iec,jsc:jec) )
  allocate( abiotic(n)%jdi14c(isc:iec,jsc:jec,nk) )
enddo

! initialize abiotic array elements
do n = 1, instances

  do j = jsd, jed
    do i = isd, ied
      abiotic(n)%po4(i,j) = 0.0
      abiotic(n)%sio4(i,j) = 0.0
    enddo
  enddo
  do j = jsc, jec
    do i = isc, iec
      abiotic(n)%sc_co2(i,j) = 0.0
      abiotic(n)%htotal(i,j) = 0.0
      abiotic(n)%csurf(i,j) = 0.0
      abiotic(n)%c14surf(i,j) = 0.0
      abiotic(n)%alpha(i,j) = 0.0
      abiotic(n)%pco2surf(i,j) = 0.0
      abiotic(n)%sal_global_wrk(i,j) = 0.0
    enddo
  enddo
  do j = jsc, jec
    do i = isc, iec
      abiotic(n)%frac_14catm(i,j) = abiotic(n)%frac_14catm_const
    enddo
  enddo
  do j = jsc, jec
    do i = isc, iec
      do k = 1, nk
        abiotic(n)%jdi14c(i,j,k) = 0.0
      enddo
    enddo
  enddo

enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_abiotic_bbc

end subroutine  ocmip2_abiotic_bbc
! </SUBROUTINE> NAME="ocmip2_abiotic_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_end">
!
! <DESCRIPTION>
!     Clean up various ABIOTIC quantities for this run.
! </DESCRIPTION>
!

subroutine ocmip2_abiotic_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
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
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_abiotic_end'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                 :: i, j, k, n
integer                                 :: lun
character(len=fm_field_name_len+1)      :: suffix
real                                    :: total_di14c
real                                    :: total_dic
real                                    :: total_di14c_bitwise
real                                    :: total_dic_bitwise
real, dimension(isd:ied,jsd:jed,nk)     :: wrk

  integer :: stdoutunit 
  stdoutunit=stdout() 
!       integrate the total concentrations of some tracers
!       for the end of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
do n = 1, instances

  total_dic = 0.0
  total_di14c = 0.0

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        total_dic = total_dic +                                         &
             t_prog(abiotic(n)%ind_dic)%field(i,j,k,taup1) *            &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_di14c = total_di14c +                                     &
             t_prog(abiotic(n)%ind_di14c)%field(i,j,k,taup1) *          &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_dic)
  call mpp_sum(total_di14c)

  write (stdoutunit,*) '  Instance ', trim(abiotic(n)%name)
  write (stdoutunit,                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol-C'')')             &
       total_dic * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DI14C  = '',es19.12,'' Gmol-C'')')           &
       total_di14c * 1.0e-09

  do k = 1, nk
    do j = jsd, jed
      do i = isd, ied
        wrk(i,j,k) =                                                            &
             t_prog(abiotic(n)%ind_dic)%field(i,j,k,taup1) *                    &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo
  total_dic_bitwise = mpp_global_sum(mpp_domain2d, wrk, BITWISE_EXACT_SUM)

  do k = 1, nk
    do j = jsd, jed
      do i = isd, ied
        wrk(i,j,k) =                                                            &
             t_prog(abiotic(n)%ind_di14c)%field(i,j,k,taup1) *                  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo
  total_di14c_bitwise = mpp_global_sum(mpp_domain2d, wrk, BITWISE_EXACT_SUM)

  write (stdoutunit,*) '  Instance ', trim(abiotic(n)%name), ' bitwise exact sum'
  write (stdoutunit,                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol-C'')')             &
       total_dic_bitwise * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DI14C  = '',es19.12,'' Gmol-C'')')           &
       total_di14c_bitwise * 1.0e-09

enddo

! save out additional information for a restart
write(stdoutunit,*)

write(stdoutunit,*) trim(note_header),                            &
     'Writing additional restart information for instances'

call ocmip2_abiotic_restart

write (stdoutunit,*) trim(note_header),                           &
     'Done writing additional restart information for instances'

do n = 1, instances

  write (stdoutunit,'(/1x,a,es16.9,a,a)')                         &
        'Annual, global, surface mean salinity = ',             &
        abiotic(n)%sal_global, ' (PSU) for instance ',          &
        trim(abiotic(n)%name)

enddo

return
end subroutine  ocmip2_abiotic_end
! </SUBROUTINE> NAME="ocmip2_abiotic_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocmip2_abiotic_restart(time_stamp)
  character(len=*),             intent(in), optional :: time_stamp
  integer :: n

  do n=1, num_restart
     call save_restart(restart(n), time_stamp)
  end do

end subroutine ocmip2_abiotic_restart
! </SUBROUTINE> NAME="ocmip2_abiotic_restart"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_abiotic_sbc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,       &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     T_prog, taum1, Grid, Time, ice_ocean_boundary_fluxes)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
integer, intent(in)                                             :: isc_bnd
integer, intent(in)                                             :: iec_bnd
integer, intent(in)                                             :: jsc_bnd
integer, intent(in)                                             :: jec_bnd
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
integer, intent(in)                                             :: taum1
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
      t_prog(abiotic(n)%ind_dic)%stf(i,j) =                                     &
            -ice_ocean_boundary_fluxes%bc(abiotic(n)%ind_co2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
      t_prog(abiotic(n)%ind_di14c)%stf(i,j) =                                   &
            -ice_ocean_boundary_fluxes%bc(abiotic(n)%ind_14co2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo
enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d(Time, Grid, abiotic(n)%id_sfc_flux_co2, t_prog(abiotic(n)%ind_dic)%stf(:,:))
   call diagnose_2d(Time, Grid, abiotic(n)%id_sfc_flux_14co2, t_prog(abiotic(n)%ind_di14c)%stf(:,:))
enddo

return

end subroutine  ocmip2_abiotic_sbc
! </SUBROUTINE> NAME="ocmip2_abiotic_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>
subroutine ocmip2_abiotic_flux_init

character(len=64), parameter    :: sub_name = 'ocmip2_abiotic_flux_init'
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

  ! Check some things
  write (stdoutunit,*)
  if (instances .eq. 0) then
    write (stdoutunit,*) trim(note_header), ' No instances'
    do_ocmip2_abiotic = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocmip2_abiotic = .true.
  endif

  module_initialized = .true.

endif


! Return if we don't want to use this package
if (.not. do_ocmip2_abiotic) then
  return
endif

if (.not. allocated(abiotic)) then

   ! allocate storage for abiotic array
  allocate ( abiotic(instances) )

  ! loop over the names, saving them into the abiotic array
  do n = 1, instances

    if (fm_get_value(path_to_names, name, index = n)) then
      abiotic(n)%name = name
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

  name = abiotic(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // name
  endif

  ! Coupler fluxes
  abiotic(n)%ind_co2_flux = aof_set_coupler_flux('co2_flux' // suffix,                          &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)

  abiotic(n)%ind_14co2_flux = aof_set_coupler_flux('c14o2_flux' // suffix,                      &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)
enddo

return

end subroutine  ocmip2_abiotic_flux_init
! </SUBROUTINE> NAME="ocmip2_abiotic_flux_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocmip2_abiotic_init

character(len=64), parameter    :: sub_name = 'ocmip2_abiotic_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

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
package_index = otpm_set_tracer_package(package_name,            &
     restart_file = default_restart_file,                        &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

! Check whether to use this package
path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif

! Check some things
write (stdoutunit,*)
if (instances .eq. 0) then
  write (stdoutunit,*) trim(note_header), ' No instances'
  do_ocmip2_abiotic = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocmip2_abiotic = .true.
endif

module_initialized = .true.

! Return if we don't want to use this package
if (.not. do_ocmip2_abiotic) then
  return
endif

! allocate storage for abiotic array
allocate ( abiotic(instances) )

! loop over the names, saving them into the abiotic array
do n = 1, instances

  if (fm_get_value(path_to_names, name, index = n)) then
    abiotic(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = abiotic(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

  ! DIC
  abiotic(n)%ind_dic = otpm_set_prog_tracer('dic' // suffix,            &
       package_name,                                                    &
       longname = 'DIC' // trim(long_suffix),                           &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                      &
       caller = caller_str)

  ! DI14C
  abiotic(n)%ind_di14c = otpm_set_prog_tracer('di14c' // suffix,        &
       package_name,                                                    &
       longname = 'DI14C' // trim(long_suffix),                         &
       units = 'mol/kg', flux_units = 'mol/m^2/s',                      &
       caller = caller_str)

enddo

! Process the namelists

!       Add the package name to the list of good namelists, to be used
!       later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

! Set up the *global* namelist
call fm_util_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call fm_util_set_value('htotal_scale_lo_in', 0.01 )
call fm_util_set_value('htotal_scale_hi_in', 100.0)
call fm_util_set_value('htotal_in', 1.0e-08)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

! Set up the instance namelists
do n = 1, instances
   ! create the instance namelist
  call fm_util_start_namelist(package_name, abiotic(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('sal_global', 35.0)                            ! PSU
  call fm_util_set_value('half_life', 5730.0)                           ! a
  call fm_util_set_value('alkbar', 2.310e-03)                           ! eq/kg
  call fm_util_set_value('frac_14catm_file', ' ')
  call fm_util_set_value('frac_14catm_name', ' ')
  call fm_util_set_value('frac_14catm_const', abiotic(n)%frac_14catm_const)
  call fm_util_set_value('po4_const', 5.0e-07)                          ! mol/kg
  call fm_util_set_value('sio4_const', 7.5e-06)                         ! mol/kg
  call fm_util_set_value('dic_global', 2.0e-03)                         ! mol/kg
  call fm_util_set_value('di14c_global', 2.0e-03)                       ! mol/kg
  call fm_util_set_value('local_restart_file', default_local_restart_file)

  ! New Wanninkhof numbers
  call fm_util_set_value('sc_co2_0', 2068.9)
  call fm_util_set_value('sc_co2_1', -118.63)
  call fm_util_set_value('sc_co2_2', 2.9311)
  call fm_util_set_value('sc_co2_3', -0.027)

  call fm_util_end_namelist(package_name, abiotic(n)%name, check = .true., caller = caller_str)

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

end subroutine ocmip2_abiotic_init
! </SUBROUTINE> NAME="ocmip2_abiotic_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocmip2_abiotic_start
! </DESCRIPTION>

subroutine ocmip2_abiotic_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,  &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)

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
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer :: i, j, m, n, nn
integer :: i_bnd_off
integer :: j_bnd_off
integer :: ind

real    :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
   ! CO2 flux
  ind = abiotic(n)%ind_co2_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,         &
         grid_tmask(isd:ied,jsd:jed,1),                                 &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
         t_prog(abiotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),     &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1) *                &
              abiotic(n)%alkbar / abiotic(n)%sal_global,                &
         abiotic(n)%po4,                                                &
         abiotic(n)%sio4,                                               &
         htotal_scale_lo, htotal_scale_hi, abiotic(n)%htotal,           &
         co2star = abiotic(n)%csurf, alpha = abiotic(n)%alpha,          &
         pco2surf = abiotic(n)%pco2surf)

!  Compute the Schmidt number of CO2 in seawater using the 
!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!  7373-7382).
    do j = jsc, jec
      do i = isc, iec
        abiotic(n)%sc_co2(i,j) =                                                &
             abiotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *         &
             (abiotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *        &
              (abiotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *       &
               abiotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (abiotic(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             abiotic(n)%alpha(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             abiotic(n)%csurf(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo

    ind = abiotic(n)%ind_14co2_flux

    ! calculate interpolated frac_14catm (fractionation of atmospheric 14CO2)
    if (abiotic(n)%frac_14catm_file .ne. ' ') then
      call time_interp_external(abiotic(n)%frac_14catm_id, model_time, abiotic(n)%frac_14catm)
    endif

    do j = jsc, jec
      do i = isc, iec
        abiotic(n)%c14surf(i,j) = abiotic(n)%csurf(i,j) *                       &
             t_prog(abiotic(n)%ind_di14c)%field(i,j,1,taum1) /                  &
             (t_prog(abiotic(n)%ind_dic)%field(i,j,1,taum1) + 1.0e-40)
        Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) = &
             abiotic(n)%alpha(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j) *       &
             (1.0 + abiotic(n)%frac_14catm(i,j) * 1.0e-03)
        Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) = &
             abiotic(n)%c14surf(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
      enddo
    enddo

  endif

enddo

return

end subroutine ocmip2_abiotic_init_sfc
! </SUBROUTINE> NAME="ocmip2_abiotic_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocmip2_abiotic_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)

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
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

integer :: i, j, n, nn
integer :: i_bnd_off
integer :: j_bnd_off
integer :: ind

real    :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances

  ind = abiotic(n)%ind_co2_flux

  call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,           &
       grid_tmask(isd:ied,jsd:jed,1),                                   &
       t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                  &
       t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                   &
       t_prog(abiotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),       &
       t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1) *                  &
            abiotic(n)%alkbar / abiotic(n)%sal_global,                  &
       abiotic(n)%po4,                                                  &
       abiotic(n)%sio4,                                                 &
       htotal_scale_lo, htotal_scale_hi, abiotic(n)%htotal,             &
       co2star = abiotic(n)%csurf, alpha = abiotic(n)%alpha,            &
       pco2surf = abiotic(n)%pco2surf)

!  Compute the Schmidt number of CO2 in seawater using the 
!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!  7373-7382).
  do j = jsc, jec
    do i = isc, iec
      abiotic(n)%sc_co2(i,j) =                                                  &
             abiotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *         &
             (abiotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *        &
              (abiotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *       &
               abiotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
      sc_no_term(i,j) = sqrt(660.0 / (abiotic(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           abiotic(n)%alpha(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           abiotic(n)%csurf(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
    enddo
  enddo

  ind = abiotic(n)%ind_14co2_flux

  ! calculate interpolated frac_14catm (fractionation of atmospheric 14CO2)
  if (abiotic(n)%frac_14catm_file .ne. ' ') then
    call time_interp_external(abiotic(n)%frac_14catm_id, model_time, abiotic(n)%frac_14catm)
  endif

  do j = jsc, jec
    do i = isc, iec
      abiotic(n)%c14surf(i,j) = abiotic(n)%csurf(i,j) *                                 &
           t_prog(abiotic(n)%ind_di14c)%field(i,j,1,taum1) /                            &
           (t_prog(abiotic(n)%ind_dic)%field(i,j,1,taum1) + 1.0e-40)
      Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +      &
           abiotic(n)%alpha(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j) *                                    &
           (1.0 + abiotic(n)%frac_14catm(i,j) * 1.0e-03)
      Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =           &
           Ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +      &
           abiotic(n)%c14surf(i,j) * rho(i,j,1,taum1) * sc_no_term(i,j)
    enddo
  enddo

enddo

return

end subroutine ocmip2_abiotic_sum_sfc
! </SUBROUTINE> NAME="ocmip2_abiotic_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_abiotic_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

do n = 1, instances
  ind = abiotic(n)%ind_co2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = abiotic(n)%ind_14co2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0
enddo

return

end subroutine ocmip2_abiotic_zero_sfc
! </SUBROUTINE> NAME="ocmip2_abiotic_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_abiotic_avg_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd, Ocean_fields, Ocean_avg_kount, grid_tmask)

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
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
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

  ind = abiotic(n)%ind_co2_flux

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

  ind = abiotic(n)%ind_14co2_flux

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

end subroutine ocmip2_abiotic_avg_sfc
! </SUBROUTINE> NAME="ocmip2_abiotic_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_abiotic_sfc_end

end subroutine ocmip2_abiotic_sfc_end
! </SUBROUTINE> NAME="ocmip2_abiotic_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_source">
!
! <DESCRIPTION>
!     compute the source terms for the ABIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!
subroutine ocmip2_abiotic_source(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, taum1, Grid, Time, rho_dzt)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
integer, intent(in)                                             :: taum1
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho_dzt

integer :: i, j, k, n
logical :: used

! calculate the source terms for ABIOTICs

! Loop over multiple instances
do n = 1, instances
  ! DI14C

  ! compute DI14C decay
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        abiotic(n)%jdi14c(i,j,k) = t_prog(abiotic(n)%ind_di14c)%field(i,j,k,taum1) *    &
             abiotic(n)%lambda
        t_prog(abiotic(n)%ind_di14c)%th_tendency(i,j,k) =                               &
             t_prog(abiotic(n)%ind_di14c)%th_tendency(i,j,k) -                          &
             abiotic(n)%jdi14c(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_sc_co2, abiotic(n)%sc_co2(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_alpha, abiotic(n)%alpha(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_csurf, abiotic(n)%csurf(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_c14surf, abiotic(n)%c14surf(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_pco2surf, abiotic(n)%pco2surf(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_htotal, abiotic(n)%htotal(:,:))
   if (abiotic(n)%id_alk .gt. 0) then
      call diagnose_2d(Time, Grid, abiotic(n)%id_alk, t_prog(indsal)%field(:,:,1,taum1) * abiotic(n)%alkbar / abiotic(n)%sal_global)
   endif
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_po4, abiotic(n)%po4(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_sio4, abiotic(n)%sio4(:,:))
   call diagnose_2d_comp(Time, Grid, abiotic(n)%id_frac_14catm, abiotic(n)%frac_14catm(:,:))
   call diagnose_3d_comp(Time, Grid, abiotic(n)%id_jdi14c, abiotic(n)%jdi14c(:,:,:))
enddo

return

end subroutine  ocmip2_abiotic_source
! </SUBROUTINE> NAME="ocmip2_abiotic_source"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!
subroutine ocmip2_abiotic_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, taup1, model_time, grid_dat, grid_tmask, grid_kmt,                 &
     grid_xt, grid_yt, grid_zt, grid_zw, grid_dzt, grid_tracer_axes,            &
     mpp_domain2d, rho_dzt)

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
integer, dimension(isd:,jsd:), intent(in)               :: grid_kmt
real, dimension(isd:,jsd:), intent(in)                  :: grid_xt
real, dimension(isd:,jsd:), intent(in)                  :: grid_yt
real, dimension(nk), intent(in)                         :: grid_zt
real, dimension(nk), intent(in)                         :: grid_zw
real, dimension(nk), intent(in)                         :: grid_dzt
integer, dimension(3), intent(in)                       :: grid_tracer_axes
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_abiotic_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

real, parameter :: sperd = 24.0 * 3600.0
real, parameter :: spery = 365.25 * sperd

integer                                         :: i, j, k, n
character(len=fm_field_name_len+1)              :: suffix
character(len=fm_field_name_len+3)              :: long_suffix
character(len=256)                              :: caller_str
real                                            :: total_di14c
real                                            :: total_dic
real                                            :: total_di14c_bitwise
real                                            :: total_dic_bitwise
real, dimension(isd:ied,jsd:jed,nk)             :: wrk
character(len=fm_string_len), allocatable       :: local_restart_file(:)
logical                                         :: fld_exist
integer                                         :: l
integer                                         :: ind
integer                                         :: id_restart

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

! dynamically allocate the global ABIOTIC arrays
call allocate_arrays(isc, iec, jsc, jec, nk, isd, ied, jsd, jed)

!       save the *global* namelist values
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str)

htotal_scale_lo_in =  fm_util_get_real   ('htotal_scale_lo_in', scalar = .true.)
htotal_scale_hi_in =  fm_util_get_real   ('htotal_scale_hi_in', scalar = .true.)
htotal_in          =  fm_util_get_real   ('htotal_in', scalar = .true.)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str)
      
! set default values for htotal_scale bounds
htotal_scale_lo(:,:) = htotal_scale_lo_in
htotal_scale_hi(:,:) = htotal_scale_hi_in

!       read in the namelists for each instance
do n = 1, instances
  call fm_util_start_namelist(package_name, abiotic(n)%name, caller = caller_str)

  abiotic(n)%frac_14catm_file       = fm_util_get_string ('frac_14catm_file', scalar = .true.)
  abiotic(n)%frac_14catm_name       = fm_util_get_string ('frac_14catm_name', scalar = .true.)
  abiotic(n)%frac_14catm_const      = fm_util_get_real   ('frac_14catm_const', scalar = .true.)
  abiotic(n)%alkbar                 = fm_util_get_real   ('alkbar', scalar = .true.)
  abiotic(n)%po4_const              = fm_util_get_real   ('po4_const', scalar = .true.)
  abiotic(n)%sio4_const             = fm_util_get_real   ('sio4_const', scalar = .true.)
  abiotic(n)%local_restart_file     = fm_util_get_string ('local_restart_file', scalar = .true.)
  abiotic(n)%sal_global             = fm_util_get_real   ('sal_global', scalar = .true.)
  abiotic(n)%dic_global             = fm_util_get_real   ('dic_global', scalar = .true.)
  abiotic(n)%di14c_global           = fm_util_get_real   ('di14c_global', scalar = .true.)
  abiotic(n)%half_life              = fm_util_get_real   ('half_life', scalar = .true.)
  abiotic(n)%sc_co2_0               = fm_util_get_real   ('sc_co2_0', scalar = .true.)
  abiotic(n)%sc_co2_1               = fm_util_get_real   ('sc_co2_1', scalar = .true.)
  abiotic(n)%sc_co2_2               = fm_util_get_real   ('sc_co2_2', scalar = .true.)
  abiotic(n)%sc_co2_3               = fm_util_get_real   ('sc_co2_3', scalar = .true.)

  call fm_util_end_namelist(package_name, abiotic(n)%name, caller = caller_str)
enddo
      
do n = 1, instances

!     Open the frac_14catm (fractionation of atmospheric 14CO2) file
!
!       If the file name is blank, then the 14C fractionation is assumed to
!       be added to the atmospheric concentration
  if (abiotic(n)%frac_14catm_file .ne. ' ') then
    abiotic(n)%frac_14catm_id = init_external_field(abiotic(n)%frac_14catm_file,        &
         abiotic(n)%frac_14catm_name, domain = mpp_domain2d, use_comp_domain = .true.)
    if (abiotic(n)%frac_14catm_id .eq. 0) then
      call mpp_error(FATAL, trim(error_header) //                                       &
           ' Could not open frac_14catm_file file: ' //                                 &
           trim(abiotic(n)%frac_14catm_file) // ' for ' // trim(abiotic(n)%frac_14catm_name))
    endif
  else
    call mpp_error(NOTE, trim(error_header) //                                          &
         ' Using constant field for atmospheric 14C for instance ' // trim(abiotic(n)%name))
    do j = jsc, jec
      do i = isc, iec
        abiotic(n)%frac_14catm(i,j) = abiotic(n)%frac_14catm_const
      enddo
    enddo
  endif

enddo

!       Read in additional information for a restart.
!
!       We must process all of the instances before restoring any files
!       as all fields must be registered before the fields are
!       restored, and fields from different instances may be in the
!       same file.
!
!       Note that the restart file names here must be different from
!       those for the tracer values.
allocate(restart(instances))
allocate(local_restart_file(instances))

write(stdoutunit,*)

do n = 1, instances

   ! Set the suffix for this instance (if instance name is "_",
   ! then use a blank suffix).
  if (abiotic(n)%name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // abiotic(n)%name
  endif

  ! Check whether we are already using this restart file, if so,
  ! we do not want to duplicate it in the list of restart files
  ! since we only read each restart file once.
  ind = 0
  do l = 1, num_restart
    if (abiotic(n)%local_restart_file == local_restart_file(l)) then
      ind = l
      exit
    endif
  end do

  if (ind .eq. 0) then
    num_restart = num_restart + 1
    ind = num_restart
    local_restart_file(ind) = trim(abiotic(n)%local_restart_file)
  end if

  ! Check whether the field already exists in the restart file.
  ! If not, then set a default value.
  fld_exist = field_exist('INPUT/' // trim(abiotic(n)%local_restart_file), 'htotal' // trim(suffix) )

  if ( fld_exist ) then
    write (stdoutunit,*) trim(note_header),                       &
         'Reading additional information for instance ',        &
         ': Initializing instance ', trim(abiotic(n)%name)
  else
    write (stdoutunit,*) trim(note_header),                       &
         'Initializing instance ', trim(abiotic(n)%name)
    abiotic(n)%htotal(:,:) = htotal_in
    ! abiotic(n)%sal_global is set via the namelist
    abiotic(n)%sal_global_wrk(:,:) = 0.0
    abiotic(n)%global_wrk_duration = 0.0
  endif

  ! Register the field for restart
  id_restart = register_restart_field(restart(ind), abiotic(n)%local_restart_file,              &
                    'htotal' // trim(suffix), abiotic(n)%htotal,                                &
                    domain=mpp_domain2d, mandatory=fld_exist )
  id_restart = register_restart_field(restart(ind), abiotic(n)%local_restart_file,              &
                    'sal_global' // trim(suffix), abiotic(n)%sal_global,                        &
                    domain = mpp_domain2d, mandatory=fld_exist )
  id_restart = register_restart_field(restart(ind), abiotic(n)%local_restart_file,              &
                    'sal_global_wrk' // trim(suffix), abiotic(n)%sal_global_wrk,                &
                    domain=mpp_domain2d, mandatory=fld_exist )
  id_restart = register_restart_field(restart(ind), abiotic(n)%local_restart_file,              &
                    'global_wrk_duration' // trim(suffix), abiotic(n)%global_wrk_duration,      &
                    domain = mpp_domain2d, mandatory=fld_exist )

enddo

! Restore the restart fields if the file exists
do l = 1, num_restart
  if (file_exist('INPUT/' // trim(local_restart_file(l)))) then
    call restore_state(restart(l))
  end if
end do

! Print the surface salinities
do n = 1, instances
  write (stdoutunit,'(/1x,a,es16.9,a,a)')                         &
        'Annual, global, surface mean salinity = ',             &
        abiotic(n)%sal_global, ' (PSU) for instance ',          &
        trim(abiotic(n)%name)
enddo

deallocate(local_restart_file)

! initialize some arrays which are held constant for this
! simulation
do n = 1, instances
  abiotic(n)%po4(:,:) = abiotic(n)%po4_const
  abiotic(n)%sio4(:,:) = abiotic(n)%sio4_const
enddo

do n = 1, instances
  if (abiotic(n)%half_life .gt. 0.0) then
    abiotic(n)%lambda = log(2.0) / (abiotic(n)%half_life * spery)
  else
    call mpp_error(FATAL,trim(error_header) // ' Half-life <= 0')
  endif
enddo

! Set up analyses

! register the instance fields
do n = 1, instances

  if (abiotic(n)%name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // abiotic(n)%name
    long_suffix = ' (' // trim(abiotic(n)%name) // ')'
  endif

  abiotic(n)%id_sc_co2 = register_diag_field(trim(diag_name),           &
       'sc_co2' // trim(suffix), grid_tracer_axes(1:2),                 &
       model_time, 'Schmidt number - CO2' // trim(long_suffix), ' ',    &
       missing_value = -1.0e+10)
  abiotic(n)%id_alpha = register_diag_field(trim(diag_name),            &
       'alpha' // trim(suffix), grid_tracer_axes(1:2),                  &
       model_time, 'Alpha CO2' // trim(long_suffix), 'mol/kg/atm',      &
       missing_value = -1.0e+10)

  abiotic(n)%id_csurf = register_diag_field(trim(diag_name),            &
       'csurf' // trim(suffix), grid_tracer_axes(1:2),                  &
       model_time, 'CO2* water' // trim(long_suffix), 'mol/kg',         &
       missing_value = -1.0e+10)

  abiotic(n)%id_c14surf = register_diag_field(trim(diag_name),          &
       'c14surf' // trim(suffix), grid_tracer_axes(1:2),                &
       model_time, 'CO2* water' // trim(long_suffix), 'mol/kg',         &
       missing_value = -1.0e+10)

  abiotic(n)%id_pco2surf = register_diag_field(trim(diag_name),         &
       'pco2surf' // trim(suffix), grid_tracer_axes(1:2),               &
       model_time, 'Oceanic pCO2' // trim(long_suffix), 'ppm',          &
       missing_value = -1.0e+10)

  abiotic(n)%id_sfc_flux_co2 = register_diag_field(trim(diag_name),                     &
       'sfc_flux_co2' // trim(suffix), grid_tracer_axes(1:2),                           &
       model_time, 'CO2 surface flux' // trim(long_suffix), 'mol m^-2 s^-1',            &
       missing_value = -1.0e+10)

  abiotic(n)%id_sfc_flux_14co2 = register_diag_field(trim(diag_name),                   &
       'sfc_flux_14co2' // trim(suffix), grid_tracer_axes(1:2),                         &
       model_time, '14CO2 surface flux' // trim(long_suffix), 'mol m^-2 s^-1',          &
       missing_value = -1.0e+10)

  abiotic(n)%id_htotal = register_diag_field(trim(diag_name),           &
       'htotal' // trim(suffix), grid_tracer_axes(1:2),                 &
       model_time, 'H+ ion concentration' // trim(long_suffix), ' ',    &
       missing_value = -1.0e+10)

  abiotic(n)%id_alk = register_diag_field(trim(diag_name),              &
       'alk' // trim(suffix), grid_tracer_axes(1:2),                    &
       model_time, 'ALK' // trim(long_suffix), 'eq kg^-1',              &
       missing_value = -1.0e+10)

  abiotic(n)%id_po4 = register_diag_field(trim(diag_name),              &
       'po4' // trim(suffix), grid_tracer_axes(1:2),                    &
       model_time, 'PO4' // trim(long_suffix), 'mol kg^-1',             &
       missing_value = -1.0e+10)

  abiotic(n)%id_sio4 = register_diag_field(trim(diag_name),             &
       'sio4' // trim(suffix), grid_tracer_axes(1:2),                   &
       model_time, 'SiO4' // trim(long_suffix), 'mol kg^-1',            &
       missing_value = -1.0e+10)

  abiotic(n)%id_frac_14catm = register_diag_field(trim(diag_name),              &
       'frac_14catm' // trim(suffix), grid_tracer_axes(1:2),                    &
       model_time, '14C fractionation' // trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)
  abiotic(n)%id_jdi14c = register_diag_field(trim(diag_name),                   &
       'jdi14c' // trim(suffix), grid_tracer_axes(1:3),                         &
       model_time, 'Restoring production' // trim(long_suffix), 'mol/kg/s',     &
       missing_value = -1.0e+10)
enddo

!       integrate the total concentrations of some tracers
!       for the start of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
do n = 1, instances

  total_dic = 0.0
  total_di14c = 0.0

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        total_dic = total_dic +                                 &
             t_prog(abiotic(n)%ind_dic)%field(i,j,k,taup1) *    &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_di14c = total_di14c +                             &
             t_prog(abiotic(n)%ind_di14c)%field(i,j,k,taup1) *  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_dic)
  call mpp_sum(total_di14c)

  write (stdoutunit,*) '  Instance ', trim(abiotic(n)%name)
  write (stdoutunit,                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol-C'')')             &
       total_dic * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DI14C  = '',es19.12,'' Gmol-C'')')           &
       total_di14c * 1.0e-09

  do k = 1, nk
    do j = jsd, jed
      do i = isd, ied
        wrk(i,j,k) =                                                            &
             t_prog(abiotic(n)%ind_dic)%field(i,j,k,taup1) *                    &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo
  total_dic_bitwise = mpp_global_sum(mpp_domain2d, wrk, BITWISE_EXACT_SUM)

  do k = 1, nk
    do j = jsd, jed
      do i = isd, ied
        wrk(i,j,k) =                                                            &
             t_prog(abiotic(n)%ind_di14c)%field(i,j,k,taup1) *                  &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo
  total_di14c_bitwise = mpp_global_sum(mpp_domain2d, wrk, BITWISE_EXACT_SUM)

  write (stdoutunit,*) '  Instance ', trim(abiotic(n)%name), ' bitwise exact sum'
  write (stdoutunit,                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol-C'')')             &
       total_dic_bitwise * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DI14C  = '',es19.12,'' Gmol-C'')')           &
       total_di14c_bitwise * 1.0e-09

enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), 'Abiotic tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocmip2_abiotic_start
! </SUBROUTINE> NAME="ocmip2_abiotic_start"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_abiotic_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!
subroutine ocmip2_abiotic_tracer(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, grid_dat, grid_tmask, grid_tcella, taum1, dtts, end_of_year)

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
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
real, dimension(nk), intent(in)                         :: grid_tcella
integer, intent(in)                                     :: taum1
real, intent(in)                                        :: dtts
logical, intent(in)                                     :: end_of_year
integer :: i, j, n
real    :: temp

! accumulate global annual means
do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
      abiotic(n)%sal_global_wrk(i,j) = abiotic(n)%sal_global_wrk(i,j) +         &
           t_prog(indsal)%field(i,j,1,taum1) *                                  &
           grid_tmask(i,j,1) * grid_dat(i,j) * dtts
    enddo
  enddo
  abiotic(n)%global_wrk_duration = abiotic(n)%global_wrk_duration + dtts
enddo

! calculate global means of at the end of the year
if (end_of_year) then

  do n = 1, instances
    temp = 0.0
    do j = jsc, jec
      do i = isc, iec
        temp = temp + abiotic(n)%sal_global_wrk(i,j)
      enddo
    enddo
    call mpp_sum(temp)
    abiotic(n)%sal_global = temp / abiotic(n)%global_wrk_duration / grid_tcella(1)
  enddo

  ! reset work variables to zero
  do n = 1, instances
    do j = jsc, jec
      do i = isc, iec
        abiotic(n)%sal_global_wrk(i,j) = 0.0
      enddo
    enddo
    abiotic(n)%global_wrk_duration = 0.0
  enddo

endif

return

end subroutine  ocmip2_abiotic_tracer
! </SUBROUTINE> NAME="ocmip2_abiotic_tracer"

end module  ocmip2_abiotic_mod
