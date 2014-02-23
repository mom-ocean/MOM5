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
! Ocean Carbon Model Intercomparison Study II: Biotic module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 Biotic
!       simulations as outlined in the Biotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Biotic/HOWTO-Biotic.html
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

!
!------------------------------------------------------------------
!
!       Module ocmip2_biotic_mod
!
!       Implementation of routines to solve the OCMIP-2 Biotic
!       simulations as outlined in the Biotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!
!------------------------------------------------------------------
!
module  ocmip2_biotic_mod

use time_manager_mod,         only: time_type
use diag_manager_mod,         only: send_data, register_diag_field, diag_axis_init
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value, fm_get_index
use fms_mod,                  only: field_exist, file_exist
use fms_io_mod,               only: register_restart_field, save_restart, restore_state
use fms_io_mod,               only: restart_file_type
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use time_manager_mod,         only: get_date
use time_interp_external_mod, only: time_interp_external, init_external_field
use mpp_domains_mod,          only: domain2d
use constants_mod,            only: WTMCO2, WTMO2

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use coupler_types_mod,  only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_grid_type, ocean_time_type
use ocmip2_co2calc_mod, only: ocmip2_co2calc
use ocean_util_mod,     only: diagnose_2d, diagnose_2d_comp, diagnose_3d_comp

implicit none

private

public  :: ocmip2_biotic_bbc
public  :: ocmip2_biotic_end
public  :: ocmip2_biotic_init
public  :: ocmip2_biotic_flux_init
public  :: ocmip2_biotic_sbc
public  :: ocmip2_biotic_source
public  :: ocmip2_biotic_start
public  :: ocmip2_biotic_init_sfc
public  :: ocmip2_biotic_avg_sfc
public  :: ocmip2_biotic_sum_sfc
public  :: ocmip2_biotic_zero_sfc
public  :: ocmip2_biotic_sfc_end
public  :: ocmip2_biotic_restart

private :: allocate_arrays
private :: locate
private :: set_array

character(len=fm_field_name_len), parameter     :: package_name = 'ocmip2_biotic'
character(len=48), parameter                    :: mod_name = 'ocmip2_biotic_mod'
character(len=48), parameter                    :: diag_name = 'ocean_ocmip2_biotic'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocmip2_biotic.res.nc'
character(len=fm_string_len), parameter         :: default_local_restart_file = 'ocmip2_biotic_local.res.nc'
character(len=fm_string_len), parameter         :: default_ice_restart_file = 'ice_ocmip2_biotic.res.nc'
character(len=fm_string_len), parameter         :: default_ocean_restart_file = 'ocmip2_biotic_airsea_flux.res.nc'

! coefficients for O2 saturation
real, parameter :: a_0 = 2.00907
real, parameter :: a_1 = 3.22014
real, parameter :: a_2 = 4.05010
real, parameter :: a_3 = 4.94457
real, parameter :: a_4 = -2.56847e-01
real, parameter :: a_5 = 3.88767
real, parameter :: b_0 = -6.24523e-03
real, parameter :: b_1 = -7.37614e-03
real, parameter :: b_2 = -1.03410e-02
real, parameter :: b_3 = -8.17083e-03
real, parameter :: c_0 = -4.88682e-07

!  sio2_const           = SiO2 concentration (mol/m^3)
!  add_phosphate        : if true, then add sufficient PO4 to keep
!                         the predicted PO4 the same as if no depletion
!                         or changed uptake rate were in effect

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

type biotic_type

  logical                               :: soft_tissue_pump
  logical                               :: add_phosphate
  real                                  :: bio_tau
  real                                  :: c_2_p
  real                                  :: ca_remin_depth
  real                                  :: caco3_2_c
  real                                  :: compensation_depth
  real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fc  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fca  _NULL
  character(len=fm_string_len)          :: local_restart_file
  real, _ALLOCATABLE, dimension(:,:)    :: flux_caco3  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: flux_poc  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: htotal  _NULL
  integer                               :: id_alpha = -1
  integer                               :: id_csurf = -1
  integer                               :: id_fc = -1
  integer                               :: id_fca = -1
  integer                               :: id_flux_caco3 = -1
  integer                               :: id_flux_poc = -1
  integer                               :: id_htotal = -1
  integer                               :: id_jca = -1
  integer                               :: id_jdop = -1
  integer                               :: id_jo2 = -1
  integer                               :: id_jpo4 = -1
  integer                               :: id_jpo4_add = -1
  integer                               :: id_jprod = -1
  integer                               :: id_pco2surf = -1
  integer                               :: id_sfc_flux_co2 = -1
  integer                               :: id_sfc_flux_o2 = -1
  integer                               :: id_vstf_alk = -1
  integer                               :: id_vstf_dic = -1
  integer                               :: id_vstf_dop = -1
  integer                               :: id_vstf_o2 = -1
  integer                               :: id_vstf_po4 = -1
  integer                               :: ind_alk
  integer                               :: ind_dic
  integer                               :: ind_dop
  integer                               :: ind_o2
  integer                               :: ind_po4
  integer                               :: ind_co2_flux
  integer                               :: ind_o2_flux
  real, _ALLOCATABLE, dimension(:,:,:)  :: jca  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jo2  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpo4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpo4_add  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_norm  _NULL
  real                                  :: kappa
  integer                               :: km_c
  real                                  :: martin_coeff
  real                                  :: martin_coeff_alt
  real                                  :: n_2_p
  character(len=fm_field_name_len)      :: name
  real                                  :: o2_min
  real                                  :: o_2_p
  real, _ALLOCATABLE, dimension(:,:)    :: pco2surf  _NULL
  real                                  :: r_bio_tau
  type(mask_region_type)                :: r_bio_tau_a
  type(mask_region_type)                :: nut_depl
  type(mask_region_type)                :: norm_remin
  type(mask_region_type)                :: no_caco3
  integer                               :: id_sc_co2 = -1
  integer                               :: id_sc_o2 = -1
  real, _ALLOCATABLE, dimension(:,:)    :: sc_co2  _NULL
  real                                  :: sc_co2_0
  real                                  :: sc_co2_1
  real                                  :: sc_co2_2
  real                                  :: sc_co2_3
  real, _ALLOCATABLE, dimension(:,:)    :: sc_o2  _NULL
  real                                  :: sc_o2_0
  real                                  :: sc_o2_1
  real                                  :: sc_o2_2
  real                                  :: sc_o2_3
  real                                  :: sigma
  real, _ALLOCATABLE, dimension(:,:)    :: sio2  _NULL
  real                                  :: sio2_const
  real                                  :: stp_alkalinity
  real                                  :: stp_salinity
  real                                  :: stp_temperature
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_alk  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_dic  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_dop  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_o2  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: vstf_po4  _NULL
  real, _ALLOCATABLE, dimension(:)      :: zfca  _NULL
  real, _ALLOCATABLE, dimension(:)      :: zforg  _NULL
  real, _ALLOCATABLE, dimension(:)      :: zforg_alt  _NULL

end type biotic_type

logical, public :: do_ocmip2_biotic

integer                                 :: indsal
integer                                 :: indtemp
integer                                 :: package_index
logical                                 :: module_initialized = .false.

character(len=128) :: version = '$Id: ocmip2_biotic.F90,v 20.0 2013/12/14 00:09:40 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

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

! Calculated parameters (with possible initial input values):
integer                                         :: id_o2_sat
integer                                         :: km_c_max
integer                                         :: po4_star_id
real, allocatable, dimension(:,:,:)             :: po4_star_t
character*128                                   :: po4_star_file    
character*128                                   :: po4_star_name
real, allocatable, dimension(:,:)               :: sc_no_term
type(biotic_type), allocatable, dimension(:)    :: biotic
integer                                         :: instances
real, allocatable, dimension(:,:)               :: o2_saturation
real, allocatable, dimension(:)                 :: tk
real, allocatable, dimension(:)                 :: ts
real, allocatable, dimension(:)                 :: ts2
real, allocatable, dimension(:)                 :: ts3
real, allocatable, dimension(:)                 :: ts4
real, allocatable, dimension(:)                 :: ts5
real, allocatable, dimension(:)                 :: tt

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

integer :: i, j, k, l, n

allocate( sc_no_term(isc:iec,jsc:jec) )
allocate( htotal_scale_lo(isc:iec,jsc:jec) )
allocate( htotal_scale_hi(isc:iec,jsc:jec) )
allocate( o2_saturation(isc:iec,jsc:jec) )
allocate( tt(isc:iec) )
allocate( tk(isc:iec) )
allocate( ts(isc:iec) )
allocate( ts2(isc:iec) )
allocate( ts3(isc:iec) )
allocate( ts4(isc:iec) )
allocate( ts5(isc:iec) )
!allocate( po4_star_t(isd:ied,jsd:jed,km_c_max) )
!       this should be dimensioned as above, but the time_interp routine
!       requires that the array dimensions match the datasets dimensions
allocate( po4_star_t(isd:ied,jsd:jed,nk) )

! initialize some arrays
sc_no_term(:,:) = 0.0
htotal_scale_lo(:,:) = 0.0
htotal_scale_hi(:,:) = 0.0
o2_saturation(:,:) = 0.0
tt(:) = 0.0
tk(:) = 0.0
ts(:) = 0.0
ts2(:) = 0.0
ts3(:) = 0.0
ts4(:) = 0.0
ts5(:) = 0.0
po4_star_t(:,:,:) = 0.0

! allocate biotic array elements
do n = 1, instances

  allocate( biotic(n)%sc_co2(isc:iec,jsc:jec) )
  allocate( biotic(n)%sc_o2(isc:iec,jsc:jec) )
  allocate( biotic(n)%htotal(isc:iec,jsc:jec) )
  allocate( biotic(n)%csurf(isc:iec,jsc:jec) )
  allocate( biotic(n)%alpha(isc:iec,jsc:jec) )
  allocate( biotic(n)%pco2surf(isc:iec,jsc:jec) )
  allocate( biotic(n)%sio2(isd:ied,jsd:jed) )
  allocate( biotic(n)%zforg(nk) )
  allocate( biotic(n)%zforg_alt(nk) )
  allocate( biotic(n)%zfca(nk) )
  allocate( biotic(n)%jprod(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpo4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jdop(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jo2(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jca(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fc(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fca(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%flux_poc(isc:iec,jsc:jec) )
  allocate( biotic(n)%flux_caco3(isc:iec,jsc:jec) )
  allocate( biotic(n)%nut_depl%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%no_caco3%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%norm_remin%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%r_bio_tau_a%mask(isc:iec,jsc:jec,12) )

  if (biotic(n)%add_phosphate) then
    allocate( biotic(n)%jprod_norm(isc:iec,jsc:jec,nk) )
    allocate( biotic(n)%jpo4_add(isc:iec,jsc:jec,nk) )
  endif

enddo

! initialize biotic array elements
do n = 1, instances

  do j = jsd, jed
    do i = isd, ied
      biotic(n)%sio2(i,j) = 0.0
      biotic(n)%sc_co2(i,j) = 0.0
      biotic(n)%sc_o2(i,j) = 0.0
    enddo
  enddo
  do j = jsc, jec
    do i = isc, iec
      biotic(n)%htotal(i,j) = 0.0
      biotic(n)%csurf(i,j) = 0.0
      biotic(n)%alpha(i,j) = 0.0
      biotic(n)%pco2surf(i,j) = 0.0
    enddo
  enddo
  do k = 1, nk
    biotic(n)%zforg(k) = 0.0
    biotic(n)%zforg_alt(k) = 0.0
    biotic(n)%zfca(k) = 0.0
  enddo
  do j = jsc, jec
    do i = isc, iec
      do k = 1, nk
        biotic(n)%jprod(i,j,k) = 0.0
        biotic(n)%jpo4(i,j,k) = 0.0
        biotic(n)%jdop(i,j,k) = 0.0
        biotic(n)%jo2(i,j,k) = 0.0
        biotic(n)%jca(i,j,k) = 0.0
        biotic(n)%fc(i,j,k) = 0.0
        biotic(n)%fca(i,j,k) = 0.0
      enddo
    enddo
  enddo
  do j = jsc, jec
    do i = isc, iec
      biotic(n)%flux_poc(i,j) = 0.0
      biotic(n)%flux_caco3(i,j) = 0.0
    enddo
  enddo
  do j = jsc, jec
    do i = isc, iec
      do l = 1, 12
        biotic(n)%nut_depl%mask(i,j,l) = 0.0
        biotic(n)%no_caco3%mask(i,j,l) = 0.0
        biotic(n)%norm_remin%mask(i,j,l) = 0.0
        biotic(n)%r_bio_tau_a%mask(i,j,l) = 0.0
      enddo
    enddo
  enddo

  if (biotic(n)%add_phosphate) then
    do j = jsc, jec
      do i = isc, iec
        do k = 1, nk
          biotic(n)%jprod_norm(i,j,k) = 0.0
          biotic(n)%jpo4_add(i,j,k) = 0.0
        enddo
      enddo
    enddo
  endif

enddo

return
end subroutine  allocate_arrays
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="locate">
!
! <DESCRIPTION>
!     After Numerical recipes:
!
!     Given an array XX of length N, and a given value of X, returns a
!     value of J such that X is between XX(J) and XX(J+1).  XX must be
!     monotonic, either increasing or decreasing. J=0 or J=N is
!     returned to indicate that X is out of range.      

!       New features:
!
!       If "period" is specified, then the array, xx, is considered
!       to be periodic with a period of "period". If "x_in" is out
!       of range, then add or subtract "period" once to attempt to 
!       make "x_in" be in range.
!
!       If "nearest" is specified, and true, then return "j" such
!       that it is the element of "xx" which is nearest to the value
!       of "x_in" (where "x_in" may have been modified by the value
!       "period", above). With this option, "j" will be in
!       the range 1 <= j <= n.
! </DESCRIPTION>
!
subroutine locate(xx , n, x_in, j, period, nearest)

integer, intent(in)             :: n
real, intent(in)                :: x_in
real, dimension(n), intent(in)  :: xx
integer, intent(out)            :: j
real, optional, intent(in)      :: period
logical, optional, intent(in)   :: nearest

integer :: jl, ju, jm
real    :: x, xt
logical :: increasing

increasing = xx(1) .lt. xx(n)

if (present(period)) then
  if (increasing) then
    ! increasing array
    if (x_in .lt. xx(1)) then
      ! original value less than start, therefore add period
      xt = x_in + period
      if (xt .gt. xx(n)) then
        ! new value greater than end
        if (abs(x_in - xx(1)) .gt. abs(xt - xx(n))) then
          ! new value closer to end than original value to start
          ! use new value
          x = xt
        else
          ! original value closer to start than new value to end
          ! use original value
          x = x_in
        endif
      else
        ! new value in range
        ! use new value
        x = xt
      endif
    elseif (x_in .gt. xx(n)) then
      ! original value greater than end, therefore subtract period
      xt = x_in - period
      if (xt .lt. xx(1)) then
        ! new value less than start
        if (abs(xt - xx(1)) .lt. abs(x_in - xx(n))) then
          ! new value closer to start than original value to end
          ! use new value
          x = xt
        else
          ! original value closer to end than new value to start
          ! use original value
          x = x_in
        endif
      else
        ! new value in range
        ! use new value
        x = xt
      endif
    else
      ! original value in range
      ! use original value
      x = x_in
    endif
  else
    ! decreasing array
    if (x_in .gt. xx(1)) then
      ! original value greater than start, therefore subtract period
      xt = x_in - period
      if (xt .lt. xx(n)) then
        ! new value less than end
        if (abs(x_in - xx(1)) .gt. abs(xt - xx(n))) then
          ! new value closer to end than original value to start
          ! use new value
          x = xt
        else
          ! original value closer to start than new value to end
          ! use original value
          x = x_in
        endif
      else
        ! new value in range
        ! use new value
        x = xt
      endif
    elseif (x_in .lt. xx(n)) then
      ! original value less than end, therefore add period
      xt = x_in + period
      if (xt .gt. xx(1)) then
        ! new value greater than start
        if (abs(xt - xx(1)) .lt. abs(x_in - xx(n))) then
          ! new value closer to start than original value to end
          ! use new value
          x = xt
        else
          ! original value closer to end than new value to start
          ! use original value
          x = x_in
        endif
      else
        ! new value in range
        ! use new value
        x = xt
      endif
    else
      ! original value in range
      ! use original value
      x = x_in
    endif
  endif
else
  ! no period specified
  ! use original value
  x = x_in
endif

jl = 0
ju = n+1
10 continue
if (ju - jl .gt. 1) then
  jm = (ju + jl) / 2
  if (increasing .eqv. (x .gt. xx(jm))) then
    jl = jm
  else
    ju = jm
  endif
  go to 10
endif
j = jl

if (present(nearest)) then
  if (nearest) then
    if (j .eq. 0) then
      j = 1
    elseif (j .lt. n) then
      if (abs(x - xx(j)) .gt. abs(x - xx(j+1))) then
        j = j + 1
      endif
    endif
  endif
endif

return
end subroutine  locate
! </SUBROUTINE> NAME="locate"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!
subroutine ocmip2_biotic_bbc(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, grid_kmt)

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt

integer  :: i, j, n, kz

!   set the bottom flux of the column for phosphate to reflect a
!   regenerative flux from the sediments where the compensation
!   depth is greater than the bottom depth
do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
      kz = grid_kmt(i,j)
      if (kz .lt. biotic(n)%km_c .and. kz .gt. 0) then
        t_prog(biotic(n)%ind_po4)%btf(i,j) =                    &
             -biotic(n)%flux_poc(i,j)
        t_prog(biotic(n)%ind_o2)%btf(i,j)  =                    &
             biotic(n)%o_2_p * biotic(n)%flux_poc(i,j)
        t_prog(biotic(n)%ind_dic)%btf(i,j) =                    &
             -(biotic(n)%c_2_p * (1.0 + biotic(n)%caco3_2_c) *  &
             biotic(n)%flux_poc(i,j))
        t_prog(biotic(n)%ind_alk)%btf(i,j) =                    &
             (biotic(n)%n_2_p - 2.0 * biotic(n)%c_2_p *         &
             biotic(n)%caco3_2_c) * biotic(n)%flux_poc(i,j)
      endif
    enddo
  enddo
enddo

return

end subroutine  ocmip2_biotic_bbc
! </SUBROUTINE> NAME="ocmip2_biotic_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>
!
subroutine ocmip2_biotic_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
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

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_end'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                 :: i, j, k, n
character(len=fm_field_name_len+1)      :: suffix
real                                    :: total_alkalinity
real                                    :: total_dic
real                                    :: total_dop
real                                    :: total_o2
real                                    :: total_phosphate

  integer :: stdoutunit 
  stdoutunit=stdout() 

! integrate the total concentrations of some tracers
! for the end of the run
total_phosphate = 0.0
total_dop = 0.0
total_o2 = 0.0
total_dic = 0.0
total_alkalinity = 0.0

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at end of run'

do n = 1, instances
  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_phosphate = total_phosphate +                     &
             t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_dop = total_dop +                                 &
             t_prog(biotic(n)%ind_dop)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_o2 = total_o2 +                                   &
             t_prog(biotic(n)%ind_o2)%field(i,j,k,taup1) *      &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_dic = total_dic +                                 &
             t_prog(biotic(n)%ind_dic)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_alkalinity = total_alkalinity +                   &
             t_prog(biotic(n)%ind_alk)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_phosphate)
  call mpp_sum(total_dop)
  call mpp_sum(total_o2)
  call mpp_sum(total_dic)
  call mpp_sum(total_alkalinity)

  write (stdoutunit,*) '  Instance ', trim(biotic(n)%name)
  write (stdoutunit,                                              &
       '(/'' Total phosphate  = '',es19.12,'' Gmol'')')         &
       total_phosphate * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DOP  = '',es19.12,'' Gmol'')')               &
       total_DOP * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total O2  = '',es19.12,'' Gmol'')')                &
       total_o2 * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol'')')               &
       total_dic * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total alkalinity  = '',es19.12,'' Geq'')')         &
       total_alkalinity * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total phosphorous  = '',es19.12,'' Gmol'')')       &
       (total_phosphate + total_dop) * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total real O2  = '',es19.12,'' Gmol'')')           &
       (total_o2 + biotic(n)%o_2_p * total_phosphate) * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total Carbon  = '',es19.12,'' Gmol'')')            &
       (total_dic + biotic(n)%c_2_p * total_dop) * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total real alkalinity  = '',es19.12,'' Geq'')')    &
       (total_alkalinity + biotic(n)%n_2_p * total_phosphate) * 1.0e-09
enddo

! save out additional information for a restart

write(stdoutunit,*)
call ocmip2_biotic_restart

do n = 1, instances

  write(stdoutunit,*) trim(note_header),                          &
       'Writing additional restart information for instance ',  &
       trim(biotic(n)%name)

  write (stdoutunit,*) trim(note_header),                         &
       'Done writing additional restart information for instance ',&
       trim(biotic(n)%name)

enddo

return
end subroutine  ocmip2_biotic_end
! </SUBROUTINE> NAME="ocmip2_biotic_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocmip2_biotic_restart(time_stamp)
  character(len=*),             intent(in), optional :: time_stamp
  integer :: n

  do n=1, num_restart
     call save_restart(restart(n), time_stamp)
  end do

end subroutine ocmip2_biotic_restart
! </SUBROUTINE> NAME="ocmip2_biotic_restart"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
subroutine ocmip2_biotic_sbc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
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

integer :: i, j, n
logical :: used

!     use the surface fluxes from the coupler
!       stf is in mol/m^2/s, flux from coupler is positive upwards
do n = 1, instances
  t_prog(biotic(n)%ind_dic)%stf(isc:iec,jsc:jec) =                      &
        -ice_ocean_boundary_fluxes%bc(biotic(n)%ind_co2_flux)%field(ind_flux)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)
  t_prog(biotic(n)%ind_o2)%stf(isc:iec,jsc:jec) =                       &
        -ice_ocean_boundary_fluxes%bc(biotic(n)%ind_o2_flux)%field(ind_flux)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)
enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d(Time, Grid, biotic(n)%id_sfc_flux_co2, t_prog(biotic(n)%ind_dic)%stf(:,:))
   call diagnose_2d(Time, Grid, biotic(n)%id_sfc_flux_o2, t_prog(biotic(n)%ind_o2)%stf(:,:))
enddo

return

end subroutine  ocmip2_biotic_sbc
! </SUBROUTINE> NAME="ocmip2_biotic_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocmip2_biotic_flux_init

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_flux_init'
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
    do_ocmip2_biotic = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocmip2_biotic = .true.
  endif

  module_initialized = .true.

endif

! Return if we don't want to use this package
if (.not. do_ocmip2_biotic) then
  return
endif

if (.not. allocated(biotic)) then

   ! allocate storage for biotic array
  allocate ( biotic(instances) )

  ! loop over the names, saving them into the biotic array
  do n = 1, instances

    if (fm_get_value(path_to_names, name, index = n)) then
      biotic(n)%name = name
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

  name = biotic(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // name
  endif

  ! Coupler fluxes
  biotic(n)%ind_co2_flux = aof_set_coupler_flux('co2_flux' // suffix,                           &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)

  biotic(n)%ind_o2_flux = aof_set_coupler_flux('o2_flux' // suffix,                             &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
       mol_wt = WTMO2, param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
       ice_restart_file = default_ice_restart_file,                                             &
       ocean_restart_file = default_ocean_restart_file,                                         &
       caller = caller_str)
enddo

return

end subroutine  ocmip2_biotic_flux_init
!</SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocmip2_biotic_init

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

real, parameter :: rho_avg = 1024.5
real, parameter :: sperd = 24.0 * 3600.0
real, parameter :: spery = 365.25 * sperd

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
  do_ocmip2_biotic = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocmip2_biotic = .true.
endif

module_initialized = .true.

! Return if we don't want to use this package
if (.not. do_ocmip2_biotic) then
  return
endif

! allocate storage for biotic array
allocate ( biotic(instances) )

! loop over the names, saving them into the biotic array
do n = 1, instances

  if (fm_get_value(path_to_names, name, index = n)) then
    biotic(n)%name = name
  else
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif

enddo

! Set up the field input
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances

  name = biotic(n)%name
  if (name(1:1) .eq. '_') then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif

  ! PO4
  biotic(n)%ind_po4 = otpm_set_prog_tracer('po4' // suffix,     &
       package_name,                                            &
       longname = 'Phosphate' // trim(long_suffix),             &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! DOP
  biotic(n)%ind_dop = otpm_set_prog_tracer('dop' // suffix,     &
       package_name,                                            &
       longname = 'DOP' // trim(long_suffix),                   &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! DIC
  biotic(n)%ind_dic = otpm_set_prog_tracer('dic' // suffix,     &
       package_name,                                            &
       longname = 'DIC' // trim(long_suffix),                   &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! O2
  biotic(n)%ind_o2 = otpm_set_prog_tracer('o2' // suffix,       &
       package_name,                                            &
       longname = 'Oxygen' // trim(long_suffix),                &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! ALK
  biotic(n)%ind_alk = otpm_set_prog_tracer('alk' // suffix,     &
       package_name,                                            &
       longname = 'Alkalinity' // trim(long_suffix),            &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

enddo

! Process the namelists

! Add the package name to the list of good namelists, to be used
! later for a consistency check
if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif

! Set up the *global* namelist
call fm_util_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call fm_util_set_value('po4_star_file', 'INPUT/po4_star_ocmip2.nc')
call fm_util_set_value('po4_star_name', 'po4_star')
call fm_util_set_value('htotal_scale_lo_in', 0.01 )
call fm_util_set_value('htotal_scale_hi_in', 100.0)
call fm_util_set_value('htotal_in', 1.0e-08)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

! Set up the instance namelists

t_mask(:) = .true.

do n = 1, instances
   ! create the instance namelist
  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('compensation_depth', 75.0)                  ! m
  call fm_util_set_value('add_phosphate', .false.)
  call fm_util_set_value('martin_coeff', 0.9)
  call fm_util_set_value('martin_coeff_alt', 0.0)
  call fm_util_set_value('ca_remin_depth', 3500.0)
  call fm_util_set_value('soft_tissue_pump', .false.)
  call fm_util_set_value('stp_temperature', 10.0)
  call fm_util_set_value('stp_salinity', 34.7)
  call fm_util_set_value('stp_alkalinity', 2370.0 * 1024.5 * 1.0e-06)    ! alkalinity is in ueq/kg, converted to eq/m^3
  call fm_util_set_value('sio2_const', 7.7e-03)                          ! mol/m^3
  call fm_util_set_value('local_restart_file', default_local_restart_file)
  call fm_util_set_value('n_2_p', 16.0)
  call fm_util_set_value('c_2_p', 117.0)
  call fm_util_set_value('o_2_p', 170.0)
  call fm_util_set_value('o2_min', 4.0 * rho_avg * 1.0e-06)
  call fm_util_set_value('bio_tau', 30.0 * sperd)
  call fm_util_set_value('sigma', 0.67)
  call fm_util_set_value('kappa', 1.0 / 0.5 / spery)
  call fm_util_set_value('caco3_2_c', 0.07)

  ! New Wanninkhof numbers
  call fm_util_set_value('sc_co2_0', 2068.9)
  call fm_util_set_value('sc_co2_1', -118.63)
  call fm_util_set_value('sc_co2_2', 2.9311)
  call fm_util_set_value('sc_co2_3', -0.027)
  call fm_util_set_value('sc_o2_0', 1929.7)
  call fm_util_set_value('sc_o2_1', -117.46)
  call fm_util_set_value('sc_o2_2', 3.116)
  call fm_util_set_value('sc_o2_3', -0.0306)

  call fm_util_end_namelist(package_name, biotic(n)%name, check = .true., caller = caller_str)

  ! create some sub-namelists
  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('factor', 0.0)
  call fm_util_set_value('coastal_only', .false.)
  call fm_util_set_value('t_mask', t_mask, size(t_mask))
  call fm_util_set_value('wlon', 0.0, index = 0)
  call fm_util_set_value('elon', 0.0, index = 0)
  call fm_util_set_value('slat', 0.0, index = 0)
  call fm_util_set_value('nlat', 0.0, index = 0)

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin', caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('factor', 0.0)
  call fm_util_set_value('coastal_only', .false.)
  call fm_util_set_value('t_mask', t_mask, size(t_mask))
  call fm_util_set_value('wlon', 0.0, index = 0)
  call fm_util_set_value('elon', 0.0, index = 0)
  call fm_util_set_value('slat', 0.0, index = 0)
  call fm_util_set_value('nlat', 0.0, index = 0)

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3', caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('factor', 0.0)
  call fm_util_set_value('coastal_only', .false.)
  call fm_util_set_value('t_mask', t_mask, size(t_mask))
  call fm_util_set_value('wlon', 0.0, index = 0)
  call fm_util_set_value('elon', 0.0, index = 0)
  call fm_util_set_value('slat', 0.0, index = 0)
  call fm_util_set_value('nlat', 0.0, index = 0)

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl', caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('factor', 0.0)
  call fm_util_set_value('coastal_only', .false.)
  call fm_util_set_value('t_mask', t_mask, size(t_mask))
  call fm_util_set_value('wlon', 0.0, index = 0)
  call fm_util_set_value('elon', 0.0, index = 0)
  call fm_util_set_value('slat', 0.0, index = 0)
  call fm_util_set_value('nlat', 0.0, index = 0)

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a', caller = caller_str)

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

end subroutine ocmip2_biotic_init
! </SUBROUTINE> NAME="ocmip2_biotic_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocmip2_biotic_start
! </DESCRIPTION>
subroutine ocmip2_biotic_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,   &
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

integer         :: i, j, n
integer         :: ind
real            :: epsln=1.0e-30

do n = 1, instances

  ! CO2 flux
  ind = biotic(n)%ind_co2_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,         &
         grid_tmask(isd:ied,jsd:jed,1),                                 &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
         t_prog(biotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_alk)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_po4)%field(isd:ied,jsd:jed,1,taum1),      &
         biotic(n)%sio2,                                                &
         htotal_scale_lo, htotal_scale_hi, biotic(n)%htotal,            &
         co2star = biotic(n)%csurf, alpha = biotic(n)%alpha,            &
         pco2surf = biotic(n)%pco2surf)

    !  Compute the Schmidt number of CO2 in seawater using the 
    !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    !  7373-7382).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_co2(i,j) =                                                 &
             biotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *          &
             (biotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *         &
              (biotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *        &
               biotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_co2(i,j) + epsln))
      enddo
    enddo
    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =   &
         biotic(n)%alpha(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =   &
         biotic(n)%csurf(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)

  endif

  ! O2 flux
  ind = biotic(n)%ind_o2_flux
  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then

!  Compute the oxygen saturation concentration at 1 atm total
!  pressure in mol/m^3 given the temperature (t, in deg C) and
!  the salinity (s, in permil)
!
!  From Garcia and Gordon (1992), Limnology and Oceonography.
!  The formula used is from page 1310, eq (8).
!
!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!  *** It shouldn't be there.                                ***
!
!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol/m^3
    do j = jsc, jec
      do i = isc, iec
        tt(i) = 298.15 - t_prog(indtemp)%field(i,j,1,taum1)
        tk(i) = 273.15 + t_prog(indtemp)%field(i,j,1,taum1)
        ts(i) = log(tt(i) / tk(i))
        ts2(i) = ts(i) * ts(i)
        ts3(i) = ts2(i) * ts(i)
        ts4(i) = ts3(i) * ts(i)
        ts5(i) = ts4(i) * ts(i)
        o2_saturation(i,j) =                                            &
             exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                         &
                 a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +                 &
                 t_prog(indsal)%field(i,j,1,taum1) *                    &
                 (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +           &
                  c_0*t_prog(indsal)%field(i,j,1,taum1)))
      enddo
    enddo

    ! convert from ml/l to mol/m^3
    do j = jsc, jec
      do i = isc, iec
        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
      enddo
    enddo

    !  Compute the Schmidt number of O2 in seawater using the 
    !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
    !  Cycles, 12, 141-163).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_o2(i,j) =                                                  &
             biotic(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *           &
             (biotic(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *          &
              (biotic(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *         &
               biotic(n)%sc_o2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_o2(i,j) + epsln))
      enddo
    enddo
    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         o2_saturation(isc:iec,jsc:jec) * sc_no_term(isc:iec,jsc:jec)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         t_prog(biotic(n)%ind_o2)%field(isc:iec,jsc:jec,1,taum1) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)

  endif
enddo

return

end subroutine ocmip2_biotic_init_sfc
! </SUBROUTINE> NAME="ocmip2_biotic_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocmip2_biotic_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
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

integer         :: i, j, n
integer         :: ind
real            :: epsln=1.0e-30

do n = 1, instances

    ind = biotic(n)%ind_co2_flux

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,         &
         grid_tmask(isd:ied,jsd:jed,1),                                 &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
         t_prog(biotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_alk)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_po4)%field(isd:ied,jsd:jed,1,taum1),      &
         biotic(n)%sio2,                                                &
         htotal_scale_lo, htotal_scale_hi, biotic(n)%htotal,            &
         co2star = biotic(n)%csurf, alpha = biotic(n)%alpha,            &
         pco2surf = biotic(n)%pco2surf)

    !  Compute the Schmidt number of CO2 in seawater using the 
    !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    !  7373-7382).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_co2(i,j) =                                                 &
             biotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *          &
             (biotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *         &
              (biotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *        &
               biotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_co2(i,j) + epsln))
      enddo
    enddo
    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
         biotic(n)%alpha(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
         biotic(n)%csurf(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)

    ind = biotic(n)%ind_o2_flux

!  Compute the oxygen saturation concentration at 1 atm total
!  pressure in mol/m^3 given the temperature (t, in deg C) and
!  the salinity (s, in permil)
!
!  From Garcia and Gordon (1992), Limnology and Oceonography.
!  The formula used is from page 1310, eq (8).
!
!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!  *** It shouldn't be there.                                ***
!
!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol/m^3
    do j = jsc, jec
      do i = isc, iec
        tt(i) = 298.15 - t_prog(indtemp)%field(i,j,1,taum1)
        tk(i) = 273.15 + t_prog(indtemp)%field(i,j,1,taum1)
        ts(i) = log(tt(i) / tk(i))
        ts2(i) = ts(i) * ts(i)
        ts3(i) = ts2(i) * ts(i)
        ts4(i) = ts3(i) * ts(i)
        ts5(i) = ts4(i) * ts(i)
        o2_saturation(i,j) =                                        &
             exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                     &
                 a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +             &
                 t_prog(indsal)%field(i,j,1,taum1) *                &
                 (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +       &
                  c_0*t_prog(indsal)%field(i,j,1,taum1)))
      enddo
    enddo

    ! convert from ml/l to mol/m^3
    do j = jsc, jec
      do i = isc, iec
        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
      enddo
    enddo

    !  Compute the Schmidt number of O2 in seawater using the 
    !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
    !  Cycles, 12, 141-163).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_o2(i,j) =                                                    &
             biotic(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *           &
             (biotic(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *          &
              (biotic(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *         &
               biotic(n)%sc_o2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_o2(i,j) + epsln))
      enddo
    enddo
    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
         o2_saturation(isc:iec,jsc:jec) * sc_no_term(isc:iec,jsc:jec)
    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
         t_prog(biotic(n)%ind_o2)%field(isc:iec,jsc:jec,1,taum1) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)

enddo

return

end subroutine ocmip2_biotic_sum_sfc
! </SUBROUTINE> NAME="ocmip2_biotic_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_biotic_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

do n = 1, instances

  ind = biotic(n)%ind_co2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = biotic(n)%ind_o2_flux

  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

enddo

return

end subroutine ocmip2_biotic_zero_sfc
! </SUBROUTINE> NAME="ocmip2_biotic_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocmip2_biotic_avg_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
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

integer         :: n
integer         :: ind
real            :: divid

divid = 1./float(Ocean_avg_kount)

do n = 1, instances

  ind = biotic(n)%ind_co2_flux

  where (grid_tmask(isc:iec,jsc:jec,1) == 1.0)
    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
  endwhere

  ind = biotic(n)%ind_o2_flux

  where (grid_tmask(isc:iec,jsc:jec,1) == 1.0)
    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
  endwhere

enddo

return

end subroutine ocmip2_biotic_avg_sfc
! </SUBROUTINE> NAME="ocmip2_biotic_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>

subroutine ocmip2_biotic_sfc_end

end subroutine ocmip2_biotic_sfc_end
! </SUBROUTINE> NAME="ocmip2_biotic_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!
subroutine ocmip2_biotic_source(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, taum1, grid_zw, grid_ht, grid_tmask, Grid, Time, rho_dzt)
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
real, dimension(nk), intent(in)                                 :: grid_zw
real, dimension(isd:,jsd:), intent(in)                          :: grid_ht
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho_dzt

integer :: i, j, k, n
integer :: ind_po4
integer :: ind_dop
integer :: ind_dic
integer :: ind_alk
integer :: ind_o2
integer :: km_c
logical :: used
integer :: day
integer :: month
integer :: year
integer :: hour
integer :: minute
integer :: second

  integer :: stdoutunit 
  stdoutunit=stdout() 

! get the model month
call get_date(Time%model_time, year, month, day,                &
              hour, minute, second)

! calculate the source terms for BIOTICs

! calculate interpolated PO4_star
call time_interp_external(po4_star_id, Time%model_time, po4_star_t)

! Loop over multiple instances
do n = 1, instances

  ind_po4 = biotic(n)%ind_po4
  ind_dic = biotic(n)%ind_dic
  ind_dop = biotic(n)%ind_dop
  ind_o2 = biotic(n)%ind_o2
  ind_alk = biotic(n)%ind_alk
  km_c = biotic(n)%km_c

  ! Production

  ! compute PO4 restoring term and correct for partial
  ! production in the bottom box
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        if (t_prog(ind_po4)%field(i,j,k,taum1) .gt.             &
            po4_star_t(i,j,k) *                                 &
            biotic(n)%nut_depl%mask(i,j,month)) then
          biotic(n)%jprod(i,j,k) =                              &
               (t_prog(ind_po4)%field(i,j,k,taum1) -            &
                po4_star_t(i,j,k) *                             &
                biotic(n)%nut_depl%mask(i,j,month)) *           &
               biotic(n)%r_bio_tau_a%mask(i,j,month) * grid_tmask(i,j,k)
        else
          biotic(n)%jprod(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jprod(i,j,km_c) = biotic(n)%jprod(i,j,km_c) *   &
           (min(grid_ht(i,j),biotic(n)%compensation_depth) - grid_zw(km_c-1)) /   &
           rho_dzt(i,j,km_c,taum1)
    enddo
  enddo

  ! Normal production (used to maintain PO4 concentrations assuming no
  ! other changes)
  if (biotic(n)%add_phosphate) then
    ! compute PO4 restoring term and correct for partial
    ! production in the bottom box
    do k = 1, km_c
      do j = jsc, jec
        do i = isc, iec
          if (t_prog(ind_po4)%field(i,j,k,taum1) .gt.           &
              po4_star_t(i,j,k)) then
            biotic(n)%jprod_norm(i,j,k) =                       &
                 (t_prog(ind_po4)%field(i,j,k,taum1) -          &
                  po4_star_t(i,j,k)) *                          &
                 biotic(n)%r_bio_tau * grid_tmask(i,j,k)
          else
            biotic(n)%jprod_norm(i,j,k) = 0.0
          endif
        enddo
      enddo
    enddo

    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jprod_norm(i,j,km_c) =                        &
             biotic(n)%jprod_norm(i,j,km_c) *                   &
             (min(grid_ht(i,j),biotic(n)%compensation_depth) - grid_zw(km_c-1)) / &
             rho_dzt(i,j,km_c,taum1)
      enddo
    enddo
  endif

  ! Particle flux
  do j = jsc, jec
    do i = isc, iec
      biotic(n)%flux_poc(i,j) = (1.0 - biotic(n)%sigma) *       &
           biotic(n)%jprod(i,j,1) * rho_dzt(i,j,1,taum1)
    enddo
  enddo

  do k = 2, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%flux_poc(i,j) = biotic(n)%flux_poc(i,j) +     &
             (1.0 - biotic(n)%sigma) * biotic(n)%jprod(i,j,k) * &
             rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! calculate the flux at the base of each layer below the
  ! compensation depth
  do k = km_c, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%fc(i,j,k) = biotic(n)%flux_poc(i,j) *         &
             grid_tmask(i,j,k) *                                &
             (biotic(n)%norm_remin%mask(i,j,month) *            &
              biotic(n)%zforg(k) +                              &
              (1.0 - biotic(n)%norm_remin%mask(i,j,month)) *    &
              biotic(n)%zforg_alt(k))
      enddo
    enddo
  enddo

  ! Calcium Carbonate

  ! calculate the formation of CaCO3 above the compensation
  ! depth
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jca(i,j,k) =                                  &
             -biotic(n)%caco3_2_c * biotic(n)%c_2_p *           &
             (1.0 - biotic(n)%sigma) *                          &
             biotic(n)%jprod(i,j,k) * grid_tmask(i,j,k) *       &
             biotic(n)%no_caco3%mask(i,j,month)
      enddo
    enddo
  enddo

  ! calculate the flux of CaCO3 at compensation depth and at
  ! bottom of each level below
  do j = jsc, jec
    do i = isc, iec
      biotic(n)%flux_caco3(i,j) =                               &
           biotic(n)%caco3_2_c * biotic(n)%c_2_p *              &
           biotic(n)%flux_poc(i,j) * biotic(n)%no_caco3%mask(i,j,month)
    enddo
  enddo

  do k=km_c,nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%fca(i,j,k) = biotic(n)%flux_caco3(i,j) *      &
             biotic(n)%zfca(k) * grid_tmask(i,j,k)
      enddo
    enddo
  enddo

  ! calculate the dissolution of CaCO3 below the compensation
  ! depth
  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jca(i,j,km_c) = biotic(n)%jca(i,j,km_c) +       &
           (biotic(n)%flux_caco3(i,j) * grid_tmask(i,j,km_c) -  &
            biotic(n)%fca(i,j,km_c) * grid_tmask(i,j,km_c+1)) / &
           rho_dzt(i,j,km_c,taum1)
    enddo
  enddo

  do k = km_c + 1, nk - 1
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jca(i,j,k) =                                  &
             (biotic(n)%fca(i,j,k-1) * grid_tmask(i,j,k) -      &
              biotic(n)%fca(i,j,k) * grid_tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jca(i,j,nk) = biotic(n)%fca(i,j,nk-1) *         &
           grid_tmask(i,j,nk) / rho_dzt(i,j,nk,taum1)
    enddo
  enddo

  ! PO4
  if (biotic(n)%add_phosphate) then
    do k = 1, km_c
      do j = jsc, jec
        do i = isc, iec
          biotic(n)%jpo4_add(i,j,k) =                           &
               biotic(n)%jprod(i,j,k) - biotic(n)%jprod_norm(i,j,k)
        enddo
      enddo
    enddo
  endif

  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jpo4(i,j,k) = -biotic(n)%jprod(i,j,k) +       &
             biotic(n)%kappa * t_prog(ind_dop)%field(i,j,k,taum1)
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jpo4(i,j,km_c) = biotic(n)%jpo4(i,j,km_c) +     &
           (biotic(n)%flux_poc(i,j) * grid_tmask(i,j,km_c) -    &
            biotic(n)%fc(i,j,km_c) * grid_tmask(i,j,km_c+1)) /  &
           rho_dzt(i,j,km_c,taum1)
    enddo
  enddo

  do k = km_c + 1, nk - 1
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jpo4(i,j,k) =                                 &
             (biotic(n)%fc(i,j,k-1) * grid_tmask(i,j,k) -       &
              biotic(n)%fc(i,j,k) * grid_tmask(i,j,k+1)) /      &
             rho_dzt(i,j,k,taum1) +                              &
             biotic(n)%kappa * t_prog(ind_dop)%field(i,j,k,taum1)
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jpo4(i,j,nk) =                                  &
           biotic(n)%fc(i,j,nk-1) * grid_tmask(i,j,nk) /        &
           rho_dzt(i,j,nk,taum1) +                               &
           biotic(n)%kappa * t_prog(ind_dop)%field(i,j,nk,taum1)
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_po4)%th_tendency(i,j,k) =                         &
             t_prog(ind_po4)%th_tendency(i,j,k) +                    &
             biotic(n)%jpo4(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  if (biotic(n)%add_phosphate) then
    do k = 1, km_c
      do j = jsc, jec
        do i = isc, iec
          t_prog(ind_po4)%th_tendency(i,j,k) =                       &
               t_prog(ind_po4)%th_tendency(i,j,k) +                  &
               biotic(n)%jpo4_add(i,j,k) * rho_dzt(i,j,k,taum1)
        enddo
      enddo
    enddo
  endif

  ! DOP
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdop(i,j,k) =                                 &
             biotic(n)%sigma * biotic(n)%jprod(i,j,k) -         &
             biotic(n)%kappa * t_prog(ind_dop)%field(i,j,k,taum1)
      enddo
    enddo
  enddo

  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdop(i,j,k) = -biotic(n)%kappa *              &
                                t_prog(ind_dop)%field(i,j,k,taum1)
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_dop)%th_tendency(i,j,k) =                         &
             t_prog(ind_dop)%th_tendency(i,j,k) +                    &
             biotic(n)%jdop(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! O2
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        if (t_prog(ind_o2)%field(i,j,k,taum1) .gt.              &
            biotic(n)%o2_min) then
          biotic(n)%jo2(i,j,k) = -biotic(n)%o_2_p *             &
                                 biotic(n)%jpo4(i,j,k)
        else
          biotic(n)%jo2(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_o2)%th_tendency(i,j,k) =                     &
             t_prog(ind_o2)%th_tendency(i,j,k) +                &
             biotic(n)%jo2(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! DIC
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_dic)%th_tendency(i,j,k) =                    &
             t_prog(ind_dic)%th_tendency(i,j,k) +               &
             (biotic(n)%c_2_p * biotic(n)%jpo4(i,j,k) +         &
              biotic(n)%jca(i,j,k)) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! ALK
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_alk)%th_tendency(i,j,k) =                    &
             t_prog(ind_alk)%th_tendency(i,j,k) +               &
             (-biotic(n)%n_2_p * biotic(n)%jpo4(i,j,k) +        &
              2.0 * biotic(n)%jca(i,j,k)) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

enddo

! Save variables for diagnostics
call diagnose_2d_comp(Time, Grid, id_o2_sat, o2_saturation(:,:))
do n = 1, instances
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_sc_co2, biotic(n)%sc_co2(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_sc_o2, biotic(n)%sc_o2(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_alpha, biotic(n)%alpha(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_csurf, biotic(n)%csurf(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_pco2surf, biotic(n)%pco2surf(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_flux_poc, biotic(n)%flux_poc(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_flux_caco3, biotic(n)%flux_caco3(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_htotal, biotic(n)%htotal(:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod, biotic(n)%jprod(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jca, biotic(n)%jca(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jpo4, biotic(n)%jpo4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jdop, biotic(n)%jdop(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jo2, biotic(n)%jo2(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fc, biotic(n)%fc(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fca, biotic(n)%fca(:,:,:))
   if (biotic(n)%add_phosphate) then
      if (biotic(n)%id_jpo4_add .gt. 0) then
         used = send_data(biotic(n)%id_jpo4_add,           &
              biotic(n)%jpo4_add(:,:,1:km_c),              &
              Time%model_time, rmask = Grid%tmask(isc:iec,jsc:jec,1:km_c))
      endif
  endif
enddo

return

end subroutine  ocmip2_biotic_source
! </SUBROUTINE> NAME="ocmip2_biotic_source"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
subroutine ocmip2_biotic_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,      &
     T_prog, taup1, model_time, grid_dat, grid_tmask, grid_kmt,                 &
     grid_xt, grid_yt, grid_zt, grid_zw, grid_dzt, grid_name, grid_tracer_axes, &
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
real, dimension(:), intent(in)                          :: grid_zt
real, dimension(:), intent(in)                          :: grid_zw
real, dimension(:), intent(in)                          :: grid_dzt
character(len=*), intent(in)                            :: grid_name
integer, dimension(3), intent(in)                       :: grid_tracer_axes
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

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
!  
integer                                         :: done
integer                                         :: i, j, k, l, n
character(len=fm_field_name_len+1)              :: suffix
character(len=fm_field_name_len+3)              :: long_suffix
character(len=256)                              :: caller_str
integer                                         :: len_w
integer                                         :: len_e
integer                                         :: len_s
integer                                         :: len_n
real                                            :: total_alkalinity
real                                            :: total_dic
real                                            :: total_dop
real                                            :: total_o2
real                                            :: total_phosphate
character(len=fm_string_len), allocatable       :: local_restart_file(:)
integer                                         :: ind
logical                                         :: fld_exist
integer                                         :: id_restart

  integer :: stdoutunit 
  stdoutunit=stdout() 

write(stdoutunit,*) 
write(stdoutunit,*) trim(note_header),                     &
                  'Starting ', trim(package_name), ' module'

! Determine indices for temperature and salinity
indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif

indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif

! dynamically allocate the global BIOTIC arrays
call allocate_arrays(isc, iec, jsc, jec, nk, isd, ied, jsd, jed)

! save the *global* namelist values
caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str)

po4_star_file      =  fm_util_get_string ('po4_star_file', scalar = .true.)
po4_star_name      =  fm_util_get_string ('po4_star_name', scalar = .true.)
htotal_scale_lo_in =  fm_util_get_real   ('htotal_scale_lo_in', scalar = .true.)
htotal_scale_hi_in =  fm_util_get_real   ('htotal_scale_hi_in', scalar = .true.)
htotal_in          =  fm_util_get_real   ('htotal_in', scalar = .true.)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str)
      
! Open up the PO4 file for restoring
po4_star_id = init_external_field(po4_star_file,                &
                                  po4_star_name,                &
                                  domain = mpp_domain2d)
if (po4_star_id .eq. 0) then
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open po4_star file: ' //                      &
       trim(po4_star_file))
endif

! set default values for htotal_scale bounds
htotal_scale_lo(:,:) = htotal_scale_lo_in
htotal_scale_hi(:,:) = htotal_scale_hi_in

! read in the namelists for each instance
do n = 1, instances

  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str)

  biotic(n)%compensation_depth     = fm_util_get_real   ('compensation_depth', scalar = .true.)
  biotic(n)%add_phosphate          = fm_util_get_logical('add_phosphate', scalar = .true.)
  biotic(n)%martin_coeff           = fm_util_get_real   ('martin_coeff', scalar = .true.)
  biotic(n)%martin_coeff_alt       = fm_util_get_real   ('martin_coeff_alt', scalar = .true.)
  biotic(n)%ca_remin_depth         = fm_util_get_real   ('ca_remin_depth', scalar = .true.)
  biotic(n)%soft_tissue_pump       = fm_util_get_logical('soft_tissue_pump', scalar = .true.)
  biotic(n)%stp_temperature        = fm_util_get_real   ('stp_temperature', scalar = .true.)
  biotic(n)%stp_salinity           = fm_util_get_real   ('stp_salinity', scalar = .true.)
  biotic(n)%stp_alkalinity         = fm_util_get_real   ('stp_alkalinity', scalar = .true.)
  biotic(n)%sio2_const             = fm_util_get_real   ('sio2_const', scalar = .true.)
  biotic(n)%local_restart_file     = fm_util_get_string ('local_restart_file', scalar = .true.)
  biotic(n)%n_2_p                  = fm_util_get_real   ('n_2_p', scalar = .true.)
  biotic(n)%c_2_p                  = fm_util_get_real   ('c_2_p', scalar = .true.)
  biotic(n)%o_2_p                  = fm_util_get_real   ('o_2_p', scalar = .true.)
  biotic(n)%o2_min                 = fm_util_get_real   ('o2_min', scalar = .true.)
  biotic(n)%bio_tau                = fm_util_get_real   ('bio_tau', scalar = .true.)
  biotic(n)%sigma                  = fm_util_get_real   ('sigma', scalar = .true.)
  biotic(n)%kappa                  = fm_util_get_real   ('kappa', scalar = .true.)
  biotic(n)%caco3_2_c              = fm_util_get_real   ('caco3_2_c', scalar = .true.)
  biotic(n)%sc_co2_0               = fm_util_get_real   ('sc_co2_0', scalar = .true.)
  biotic(n)%sc_co2_1               = fm_util_get_real   ('sc_co2_1', scalar = .true.)
  biotic(n)%sc_co2_2               = fm_util_get_real   ('sc_co2_2', scalar = .true.)
  biotic(n)%sc_co2_3               = fm_util_get_real   ('sc_co2_3', scalar = .true.)
  biotic(n)%sc_o2_0                = fm_util_get_real   ('sc_o2_0', scalar = .true.)
  biotic(n)%sc_o2_1                = fm_util_get_real   ('sc_o2_1', scalar = .true.)
  biotic(n)%sc_o2_2                = fm_util_get_real   ('sc_o2_2', scalar = .true.)
  biotic(n)%sc_o2_3                = fm_util_get_real   ('sc_o2_3', scalar = .true.)

  call fm_util_end_namelist(package_name, biotic(n)%name, caller = caller_str)

  biotic(n)%r_bio_tau = 1.0 / biotic(n)%bio_tau

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin', caller = caller_str)

  biotic(n)%norm_remin%factor        =  fm_util_get_real          ('factor', scalar = .true.)
  biotic(n)%norm_remin%coastal_only  =  fm_util_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%norm_remin%wlon          => fm_util_get_real_array    ('wlon')
  biotic(n)%norm_remin%elon          => fm_util_get_real_array    ('elon')
  biotic(n)%norm_remin%slat          => fm_util_get_real_array    ('slat')
  biotic(n)%norm_remin%nlat          => fm_util_get_real_array    ('nlat')
  biotic(n)%norm_remin%t_mask        => fm_util_get_logical_array ('t_mask')

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin', caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3', caller = caller_str)

  biotic(n)%no_caco3%factor          =  fm_util_get_real          ('factor', scalar = .true.)
  biotic(n)%no_caco3%coastal_only    =  fm_util_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%no_caco3%wlon            => fm_util_get_real_array    ('wlon')
  biotic(n)%no_caco3%elon            => fm_util_get_real_array    ('elon')
  biotic(n)%no_caco3%slat            => fm_util_get_real_array    ('slat')
  biotic(n)%no_caco3%nlat            => fm_util_get_real_array    ('nlat')
  biotic(n)%no_caco3%t_mask          => fm_util_get_logical_array ('t_mask')

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3', caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl', caller = caller_str)

  biotic(n)%nut_depl%factor          =  fm_util_get_real          ('factor', scalar = .true.)
  biotic(n)%nut_depl%coastal_only    =  fm_util_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%nut_depl%wlon            => fm_util_get_real_array    ('wlon')
  biotic(n)%nut_depl%elon            => fm_util_get_real_array    ('elon')
  biotic(n)%nut_depl%slat            => fm_util_get_real_array    ('slat')
  biotic(n)%nut_depl%nlat            => fm_util_get_real_array    ('nlat')
  biotic(n)%nut_depl%t_mask          => fm_util_get_logical_array ('t_mask')

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl', caller = caller_str)

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a', caller = caller_str)

  biotic(n)%r_bio_tau_a%factor       =  fm_util_get_real          ('factor', scalar = .true.)
  biotic(n)%r_bio_tau_a%coastal_only =  fm_util_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%r_bio_tau_a%wlon         => fm_util_get_real_array    ('wlon')
  biotic(n)%r_bio_tau_a%elon         => fm_util_get_real_array    ('elon')
  biotic(n)%r_bio_tau_a%slat         => fm_util_get_real_array    ('slat')
  biotic(n)%r_bio_tau_a%nlat         => fm_util_get_real_array    ('nlat')
  biotic(n)%r_bio_tau_a%t_mask       => fm_util_get_logical_array ('t_mask')

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a', caller = caller_str)

enddo

! calculate the index for the box containing the compensation depth
km_c_max = 0
do n = 1, instances
  call locate(grid_zw, nk, biotic(n)%compensation_depth,        &
              biotic(n)%km_c, nearest = .true.)
  if (grid_zw(biotic(n)%km_c) .lt.                              &
      biotic(n)%compensation_depth) then
    biotic(n)%km_c = biotic(n)%km_c + 1
  endif

  write (stdoutunit,*) trim(note_header),                         &
                     'The compensation depth for instance ',    &
                     n, ', ', biotic(n)%compensation_depth,     &
                     ' m, occurs in box ', biotic(n)%km_c ,     &
                     ' between depths ',                        &
                     grid_zw(biotic(n)%km_c-1),                 &
                     ' m and ', grid_zw(biotic(n)%km_c), ' m'
  km_c_max = max(km_c_max, biotic(n)%km_c)
enddo

! read in the norm_remin namelist data
do n = 1, instances
  if (associated(biotic(n)%norm_remin%wlon)) then
    len_w = size(biotic(n)%norm_remin%wlon)
  else
    len_w = 0
  endif
  if (associated(biotic(n)%norm_remin%elon)) then
    len_e = size(biotic(n)%norm_remin%elon)
  else
    len_e = 0
  endif
  if (associated(biotic(n)%norm_remin%slat)) then
    len_s = size(biotic(n)%norm_remin%slat)
  else
    len_s = 0
  endif
  if (associated(biotic(n)%norm_remin%nlat)) then
    len_n = size(biotic(n)%norm_remin%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif

  if (size(biotic(n)%norm_remin%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif

  ! set all of the values to the default
  biotic(n)%norm_remin%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process norm_remin array for ', trim(biotic(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
    done = 0
    do l = 1, 12
      if (biotic(n)%norm_remin%t_mask(l)) then
        if (done .eq. 0) then

           ! set the values via the input values, saving this time index
           ! afterwards
          write (stdoutunit,*) 'Assigning month ', l
          call set_array(biotic(n)%norm_remin%mask(:,:,l),        &
                  isd, ied, jsd, jed, grid_xt, grid_yt, grid_kmt, &
                  len_w, biotic(n)%norm_remin%wlon,               &
                  biotic(n)%norm_remin%elon,                      &
                  biotic(n)%norm_remin%slat,                      &
                  biotic(n)%norm_remin%nlat,                      &
                  biotic(n)%norm_remin%factor, 1.0,               &
                  'Normal remineralization', biotic(n)%norm_remin%coastal_only)
          done = l
        else
           ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          biotic(n)%norm_remin%mask(:,:,l) = biotic(n)%norm_remin%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

! read in the no_caco3 namelist data
do n = 1, instances
  if (associated(biotic(n)%no_caco3%wlon)) then
    len_w = size(biotic(n)%no_caco3%wlon)
  else
    len_w = 0
  endif
  if (associated(biotic(n)%no_caco3%elon)) then
    len_e = size(biotic(n)%no_caco3%elon)
  else
    len_e = 0
  endif
  if (associated(biotic(n)%no_caco3%slat)) then
    len_s = size(biotic(n)%no_caco3%slat)
  else
    len_s = 0
  endif
  if (associated(biotic(n)%no_caco3%nlat)) then
    len_n = size(biotic(n)%no_caco3%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif

  if (size(biotic(n)%no_caco3%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif

  ! set all of the values to the default
  biotic(n)%no_caco3%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process no_caco3 array for ', trim(biotic(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
    done = 0
    do l = 1, 12
      if (biotic(n)%no_caco3%t_mask(l)) then
        if (done .eq. 0) then
           ! set the values via the input values, saving this time index
           ! afterwards
          write (stdoutunit,*) 'Assigning month ', l
          call set_array(biotic(n)%no_caco3%mask(:,:,l),                &
                  isd, ied, jsd, jed, grid_xt, grid_yt, grid_kmt,       &
                  len_w, biotic(n)%no_caco3%wlon,                       &
                  biotic(n)%no_caco3%elon,                              &
                  biotic(n)%no_caco3%slat,                              &
                  biotic(n)%no_caco3%nlat,                              &
                  biotic(n)%no_caco3%factor, 1.0,                       &
                  'Carbonate inhibition', biotic(n)%no_caco3%coastal_only)
          done = l
        else

           ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          biotic(n)%no_caco3%mask(:,:,l) = biotic(n)%no_caco3%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

! read in the nut_depl namelist data
do n = 1, instances
  if (associated(biotic(n)%nut_depl%wlon)) then
    len_w = size(biotic(n)%nut_depl%wlon)
  else
    len_w = 0
  endif
  if (associated(biotic(n)%nut_depl%elon)) then
    len_e = size(biotic(n)%nut_depl%elon)
  else
    len_e = 0
  endif
  if (associated(biotic(n)%nut_depl%slat)) then
    len_s = size(biotic(n)%nut_depl%slat)
  else
    len_s = 0
  endif
  if (associated(biotic(n)%nut_depl%nlat)) then
    len_n = size(biotic(n)%nut_depl%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif

  if (size(biotic(n)%nut_depl%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif

  ! set all of the values to the default
  biotic(n)%nut_depl%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process nut_depl array for ', trim(biotic(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
    done = 0
    do l = 1, 12
      if (biotic(n)%nut_depl%t_mask(l)) then
        if (done .eq. 0) then
           ! set the values via the input values, saving this time index
           ! afterwards
          write (stdoutunit,*) 'Assigning month ', l
          call set_array(biotic(n)%nut_depl%mask(:,:,l),                &
                  isd, ied, jsd, jed, grid_xt, grid_yt, grid_kmt,       &
                  len_w, biotic(n)%nut_depl%wlon,                       &
                  biotic(n)%nut_depl%elon,                              &
                  biotic(n)%nut_depl%slat,                              &
                  biotic(n)%nut_depl%nlat,                              &
                  biotic(n)%nut_depl%factor, 1.0,                       &
                  'Nutrient depletion', biotic(n)%nut_depl%coastal_only)
          done = l
        else
           ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          biotic(n)%nut_depl%mask(:,:,l) = biotic(n)%nut_depl%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

! read in the r_bio_tau_a namelist data
do n = 1, instances

  if (associated(biotic(n)%r_bio_tau_a%wlon)) then
    len_w = size(biotic(n)%r_bio_tau_a%wlon)
  else
    len_w = 0
  endif
  if (associated(biotic(n)%r_bio_tau_a%elon)) then
    len_e = size(biotic(n)%r_bio_tau_a%elon)
  else
    len_e = 0
  endif
  if (associated(biotic(n)%r_bio_tau_a%slat)) then
    len_s = size(biotic(n)%r_bio_tau_a%slat)
  else
    len_s = 0
  endif
  if (associated(biotic(n)%r_bio_tau_a%nlat)) then
    len_n = size(biotic(n)%r_bio_tau_a%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif

  if (size(biotic(n)%r_bio_tau_a%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif

  ! set all of the values to the default
  biotic(n)%r_bio_tau_a%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process r_bio_tau_a array for ', trim(biotic(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
    done = 0
    do l = 1, 12
      if (biotic(n)%r_bio_tau_a%t_mask(l)) then
        if (done .eq. 0) then
           ! set the values via the input values, saving this time index
           ! afterwards
          write (stdoutunit,*) 'Assigning month ', l
          call set_array(biotic(n)%r_bio_tau_a%mask(:,:,l),             &
                  isd, ied, jsd, jed, grid_xt, grid_yt, grid_kmt,       &
                  len_w, biotic(n)%r_bio_tau_a%wlon,                    &
                  biotic(n)%r_bio_tau_a%elon,                           &
                  biotic(n)%r_bio_tau_a%slat,                           &
                  biotic(n)%r_bio_tau_a%nlat,                           &
                  biotic(n)%r_bio_tau_a%factor, 1.0,                    &
                  'Primary production limitation', biotic(n)%r_bio_tau_a%coastal_only)
          done = l
        else
           ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          biotic(n)%r_bio_tau_a%mask(:,:,l) = biotic(n)%r_bio_tau_a%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

! multiply by the restoring factor
do n = 1, instances
  biotic(n)%r_bio_tau_a%mask(:,:,:) =                                &
       biotic(n)%r_bio_tau * biotic(n)%r_bio_tau_a%mask(:,:,:)
enddo

! initialize special arrays for remineralization
do n = 1, instances
  do k = 1, nk
    biotic(n)%zforg(k) =                                        &
         (grid_zw(k) /                                          &
          biotic(n)%compensation_depth) ** (-biotic(n)%martin_coeff)
    biotic(n)%zfca(k) =                                         &
         exp(-(grid_zw(k) -                                     &
               biotic(n)%compensation_depth) / biotic(n)%ca_remin_depth)
  enddo
enddo

! initialize special arrays for alternate remineralization
do n = 1, instances
  do k = 1, nk
    biotic(n)%zforg_alt(k) = (grid_zw(k) /                      &
         biotic(n)%compensation_depth) ** (-biotic(n)%martin_coeff_alt)
  enddo
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
  if (biotic(n)%name(1:1) .eq. '_') then
    suffix = ' '
  else
    suffix = '_' // biotic(n)%name
  endif

  ! Check whether we are already using this restart file, if so,
  ! we do not want to duplicate it in the list of restart files
  ! since we only read each restart file once.
  ind = 0
  do l = 1, num_restart
    if (biotic(n)%local_restart_file == local_restart_file(l)) then
      ind = l
      exit
    endif
  end do

  if (ind .eq. 0) then
    num_restart = num_restart + 1
    ind = num_restart
    local_restart_file(ind) = trim(biotic(n)%local_restart_file)
  end if

  ! Check whether the field already exists in the restart file.
  ! If not, then set a default value.
  fld_exist = field_exist('INPUT/' // trim(biotic(n)%local_restart_file), 'htotal' // trim(suffix) )

  if ( fld_exist ) then
    write (stdoutunit,*) trim(note_header),                       &
         'Reading additional information for instance ',        &
         ': Initializing instance ', trim(biotic(n)%name)
  else
    write (stdoutunit,*) trim(note_header),                       &
         'Initializing instance ', trim(biotic(n)%name)
    biotic(n)%htotal(:,:) = htotal_in
  endif

  ! Register the field for restart
  id_restart = register_restart_field(restart(ind), biotic(n)%local_restart_file,       &
                    'htotal' // trim(suffix), biotic(n)%htotal,                         &
                    domain=mpp_domain2d, mandatory=fld_exist )

enddo

! Restore the restart fields if the file exists
do l = 1, num_restart
  if (file_exist('INPUT/' // trim(local_restart_file(l)))) then
    call restore_state(restart(l))
  end if
end do

deallocate(local_restart_file)

! initialize some arrays which are held constant for this
! simulation
do n = 1, instances
  biotic(n)%sio2(:,:) = biotic(n)%sio2_const
enddo

! Set up analyses

! register the global fields
suffix = '_' // package_name
long_suffix = ' (' // trim(package_name) // ')'

id_o2_sat = register_diag_field(trim(diag_name),                        &
     'o2_saturation' // trim(suffix), grid_tracer_axes(1:2),            &
     model_time, 'O2 saturation' // trim(long_suffix), 'mol/m^3/atm',   &
     missing_value = -1.0e+10)

! register the instance fields
do n = 1, instances

  if (instances .eq. 1) then
    suffix = ' '
    long_suffix = ' '
  else
    suffix = '_' // biotic(n)%name
    long_suffix = ' (' // trim(biotic(n)%name) // ')'
  endif

  biotic(n)%id_sc_co2 = register_diag_field(trim(diag_name),                    &
       'sc_co2' // trim(suffix), grid_tracer_axes(1:2),                         &
       model_time, 'Schmidt number - CO2' // trim(long_suffix), ' ',            &
       missing_value = -1.0e+10)

  biotic(n)%id_sc_o2 = register_diag_field(trim(diag_name),                     &
       'sc_o2' // trim(suffix), grid_tracer_axes(1:2),                          &
       model_time, 'Schmidt number - O2' // trim(long_suffix), ' ',             &
       missing_value = -1.0e+10)

  biotic(n)%id_alpha = register_diag_field(trim(diag_name),                     &
       'alpha' // trim(suffix), grid_tracer_axes(1:2),                          &
       model_time, 'Alpha CO2' // trim(long_suffix), 'mol/kg/atm',              &
       missing_value = -1.0e+10)

  biotic(n)%id_csurf = register_diag_field(trim(diag_name),                     &
       'csurf' // trim(suffix), grid_tracer_axes(1:2),                          &
       model_time, 'CO2* water' // trim(long_suffix), 'mol/kg',                 &
       missing_value = -1.0e+10)

  biotic(n)%id_pco2surf = register_diag_field(trim(diag_name),                  &
       'pco2surf' // trim(suffix), grid_tracer_axes(1:2),                       &
       model_time, 'Oceanic pCO2' // trim(long_suffix), 'ppm',                  &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_poc = register_diag_field(trim(diag_name),                  &
       'flux_poc' // trim(suffix), grid_tracer_axes(1:2),                       &
       model_time, 'POC flux' // trim(long_suffix), 'mol/m^2/s',                &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_caco3 = register_diag_field(trim(diag_name),                &
       'flux_caco3' // trim(suffix), grid_tracer_axes(1:2),                     &
       model_time, 'CaCO3 flux' // trim(long_suffix), 'mol/m^2/s',              &
       missing_value = -1.0e+10)

  biotic(n)%id_htotal = register_diag_field(trim(diag_name),                    &
       'htotal' // trim(suffix), grid_tracer_axes(1:2),                         &
       model_time, 'H+ ion concentration' // trim(long_suffix), ' ',            &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod = register_diag_field(trim(diag_name),                     &
       'jprod' // trim(suffix), grid_tracer_axes(1:3),                          &
       model_time, 'Restoring production' // trim(long_suffix), 'mol/kg/s',     &
       missing_value = -1.0e+10)

  biotic(n)%id_jca = register_diag_field(trim(diag_name),                       &
       'jca' // trim(suffix), grid_tracer_axes(1:3),                            &
       model_time, 'CaCO3 change' // trim(long_suffix), 'mol/kg/s',             &
       missing_value = -1.0e+10)

  biotic(n)%id_jpo4 = register_diag_field(trim(diag_name),                      &
       'jpo4' // trim(suffix), grid_tracer_axes(1:3),                           &
       model_time, 'PO4 source' // trim(long_suffix), 'mol/kg/s',               &
       missing_value = -1.0e+10)

  biotic(n)%id_jdop = register_diag_field(trim(diag_name),                      &
       'jdop' // trim(suffix), grid_tracer_axes(1:3),                           &
       model_time, 'DOP source' // trim(long_suffix), 'mol/kg/s',               &
       missing_value = -1.0e+10)

  biotic(n)%id_jo2 = register_diag_field(trim(diag_name),                       &
       'jo2' // trim(suffix), grid_tracer_axes(1:3),                            &
       model_time, 'O2 source' // trim(long_suffix), 'mol/kg/s',                &
       missing_value = -1.0e+10)

  biotic(n)%id_fc = register_diag_field(trim(diag_name),                        &
       'fc' // trim(suffix), grid_tracer_axes(1:3),                             &
       model_time, 'POP change' // trim(long_suffix), 'mol/m^2/s',              &
       missing_value = -1.0e+10)

  biotic(n)%id_fca = register_diag_field(trim(diag_name),                       &
       'fca' // trim(suffix), grid_tracer_axes(1:3),                            &
       model_time, 'CaCO3 change' // trim(long_suffix), 'mol/m^2/s',            &
       missing_value = -1.0e+10)

  if (biotic(n)%add_phosphate) then
    biotic(n)%id_jpo4_add = register_diag_field(trim(diag_name),                &
         'jpo4_add' // trim(suffix), grid_tracer_axes(1:3),                     &
         model_time, 'Additional phosphate' // trim(long_suffix), 'mol/kg/s',   &
         missing_value = -1.0e+10)
  endif

enddo

! integrate the total concentrations of some tracers
! for the start of the run
total_phosphate = 0.0
total_dop = 0.0
total_o2 = 0.0
total_dic = 0.0
total_alkalinity = 0.0

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at start of run'

do n = 1, instances
  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_phosphate = total_phosphate +                     &
             t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_dop = total_dop +                                 &
             t_prog(biotic(n)%ind_dop)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_o2 = total_o2 +                                   &
             t_prog(biotic(n)%ind_o2)%field(i,j,k,taup1) *      &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_dic = total_dic +                                 &
             t_prog(biotic(n)%ind_dic)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_alkalinity = total_alkalinity +                   &
             t_prog(biotic(n)%ind_alk)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo
    enddo
  enddo

  call mpp_sum(total_phosphate)
  call mpp_sum(total_dop)
  call mpp_sum(total_o2)
  call mpp_sum(total_dic)
  call mpp_sum(total_alkalinity)

  write (stdoutunit,*) '  Instance ', trim(biotic(n)%name)
  write (stdoutunit,                                              &
       '(/'' Total phosphate  = '',es19.12,'' Gmol'')')         &
       total_phosphate * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DOP  = '',es19.12,'' Gmol'')')               &
       total_DOP * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total O2  = '',es19.12,'' Gmol'')')                &
       total_o2 * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol'')')               &
       total_dic * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total alkalinity  = '',es19.12,'' Geq'')')         &
       total_alkalinity * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total phosphorous  = '',es19.12,'' Gmol'')')       &
       (total_phosphate + total_dop) * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total real O2  = '',es19.12,'' Gmol'')')           &
       (total_o2 + biotic(n)%o_2_p * total_phosphate) * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total Carbon  = '',es19.12,'' Gmol'')')            &
       (total_dic + biotic(n)%c_2_p * total_dop) * 1.0e-09
  write (stdoutunit,                                              &
       '(/'' Total real alkalinity  = '',es19.12,'' Geq'')')    &
       (total_alkalinity + biotic(n)%n_2_p * total_phosphate) * 1.0e-09
enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), 'Tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocmip2_biotic_start
! </SUBROUTINE> NAME="ocmip2_biotic_start"

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
integer, dimension(isd:,jsd:), intent(in)       :: kmt
character(len=*), intent(in)                    :: name
real, dimension(num_regions), intent(in)        :: nlat
real, intent(in)                                :: set_value
real, dimension(num_regions), intent(in)        :: slat
real, intent(in)                                :: unset_value
real, dimension(num_regions), intent(in)        :: wlon_in
real, dimension(isd:,jsd:), intent(in)          :: xt
real, dimension(isd:,jsd:), intent(in)          :: yt

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

! if desired only apply mask to coastal regions
if (coastal_only) then
  do j = jsd, jed
    do i = isd, ied
      if (kmt(i,j) .ne. 0 .and.                         &
          array(i,j) .eq. set_value) then
         ! if all the surrounding points are ocean, then this is not
         ! a coastal point, therefore reset the mask
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

! clean up
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

real, intent(in)                :: x_in
real, intent(in)                :: w
real, intent(in)                :: e

logical :: lon_between
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

end module  ocmip2_biotic_mod
