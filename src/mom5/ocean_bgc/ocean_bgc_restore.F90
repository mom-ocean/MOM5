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
!------------------------------------------------------------------
!
!       Module ocean_bgc_restore_mod
!
!       Implementation of routines to solve the OCMIP-2 Biotic
!       simulations as outlined in the Biotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!
!------------------------------------------------------------------
!

module  ocean_bgc_restore_mod

use time_manager_mod,         only: time_type
use diag_manager_mod,         only: register_diag_field, diag_axis_init
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

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer
use fm_util_mod,        only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod,        only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod,        only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod,        only: fm_util_start_namelist, fm_util_end_namelist
use coupler_types_mod,  only: ind_alpha, ind_csurf, coupler_2d_bc_type, ind_flux
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type, ocean_grid_type, ocean_time_type
use ocmip2_co2calc_mod, only: ocmip2_co2calc
use ocean_util_mod,     only: diagnose_2d, diagnose_2d_comp, diagnose_3d_comp

implicit none

private

public  :: ocean_bgc_restore_bbc
public  :: ocean_bgc_restore_end
public  :: ocean_bgc_restore_init
public  :: ocean_bgc_restore_flux_init
public  :: ocean_bgc_restore_sbc
public  :: ocean_bgc_restore_source
public  :: ocean_bgc_restore_start
public  :: ocean_bgc_restore_init_sfc
public  :: ocean_bgc_restore_avg_sfc
public  :: ocean_bgc_restore_sum_sfc
public  :: ocean_bgc_restore_zero_sfc
public  :: ocean_bgc_restore_sfc_end
public  :: ocean_bgc_restore_restart


private :: allocate_arrays
private :: locate
private :: set_array


character(len=32), parameter              :: package_name = 'ocean_bgc_restore'
character(len=48), parameter              :: mod_name = 'ocean_bgc_restore_mod'
character(len=48), parameter              :: diag_name = 'ocean_bgc_restore'
character(len=fm_string_len), parameter   :: default_restart_file = 'ocean_bgc_restore.res.nc'
character(len=fm_string_len), parameter   :: default_local_restart_file = 'ocean_bgc_restore_local.res.nc'
character(len=fm_string_len), parameter   :: default_ice_restart_file = 'ice_bgc_restore.res.nc'
character(len=fm_string_len), parameter   :: default_ocean_restart_file = 'ocean_bgc_restore_airsea_flux.res.nc'

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

  real                                  :: bio_tau
  real                                  :: bio_tau_don
  real                                  :: bio_tau_dop
  real                                  :: bio_tau_fix
  real                                  :: bio_tau_ldoc
  real                                  :: bio_tau_nh4
  real                                  :: bio_tau_nitrif_d
  real                                  :: bio_tau_nitrif_s
  real                                  :: c_2_n
  real                                  :: ca_remin_depth
  real                                  :: si_remin_depth
  real                                  :: compensation_depth
  real, _ALLOCATABLE, dimension(:,:)    :: comp_depth_frac  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL
  logical                               :: fe_ballast_assoc
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpon  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fpop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fcaco3  _NULL
  real                                  :: fdetL0
  real                                  :: fdetS0
  character(len=fm_string_len)          :: local_restart_file
  real, _ALLOCATABLE, dimension(:,:)    :: flux_caco3  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: flux_pon  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: flux_pop  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: flux_sio2  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fracl  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: fsio2  _NULL
  real                                  :: gamma_det
  real                                  :: global_wrk_duration = 0.0
  real, _ALLOCATABLE, dimension(:,:)    :: htotal  _NULL
  integer                               :: id_alpha = -1
  integer                               :: id_csurf = -1
  integer                               :: id_fpon = -1
  integer                               :: id_fpop = -1
  integer                               :: id_fcaco3 = -1
  integer                               :: id_flux_caco3 = -1
  integer                               :: id_flux_pon = -1
  integer                               :: id_flux_pop = -1
  integer                               :: id_flux_sio2 = -1
  integer                               :: id_fracl = -1
  integer                               :: id_fsio2 = -1
  integer                               :: id_htotal = -1
  integer                               :: id_jcaco3 = -1
  integer                               :: id_jdenit = -1
  integer                               :: id_jdon = -1
  integer                               :: id_jdop = -1
  integer                               :: id_jfe_ads = -1
  integer                               :: id_jfe_des = -1
  integer                               :: id_jfe_graz = -1
  integer                               :: id_jfe_sink = -1
  integer                               :: id_jldoc = -1
  integer                               :: id_jnh4 = -1
  integer                               :: id_jnh4_graz = -1
  integer                               :: id_jno3 = -1
  integer                               :: id_jo2 = -1
  integer                               :: id_jpo4 = -1
  integer                               :: id_jpo4_graz = -1
  integer                               :: id_jpofe = -1
  integer                               :: id_jpon = -1
  integer                               :: id_jpop = -1
  integer                               :: id_jprod_alk = -1
  integer                               :: id_jprod_fed = -1
  integer                               :: id_jprod_n_fix = -1
  integer                               :: id_jprod_nh4 = -1
  integer                               :: id_jprod_no3 = -1
  integer                               :: id_jprod_n_norm = -1
  integer                               :: id_jprod_p_norm = -1
  integer                               :: id_jprod_p_fix = -1
  integer                               :: id_jprod_po4 = -1
  integer                               :: id_jprod_pofe = -1
  integer                               :: id_jprod_pon = -1
  integer                               :: id_jprod_pop = -1
  integer                               :: id_jprod_sio4 = -1
  integer                               :: id_jsio4 = -1
  integer                               :: id_pco2surf = -1
  integer                               :: id_sfc_flux_co2 = -1
  integer                               :: id_sfc_flux_o2 = -1
  integer                               :: id_sfc_flux_fed = -1
  integer                               :: ind_alk
  integer                               :: ind_dic
  integer                               :: ind_don
  integer                               :: ind_dop
  integer                               :: ind_fed
  integer                               :: ind_fep
  integer                               :: ind_ldoc
  integer                               :: ind_nh4
  integer                               :: ind_no3
  integer                               :: ind_o2
  integer                               :: ind_po4
  integer                               :: ind_sio4
  integer                               :: ind_co2_flux
  integer                               :: ind_o2_flux
  real, _ALLOCATABLE, dimension(:,:,:)  :: jcaco3  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdenit  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdon  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_ads  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_des  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_graz  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_sink  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jldoc  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jnh4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jnh4_graz  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jno3  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jo2  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpo4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpo4_graz  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpofe  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpon  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_alk  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_fed  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_n_fix  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_nh4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_no3  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_n_norm  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_p_norm  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_p_fix  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_po4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pofe  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pon  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_sio4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jsio4  _NULL
  real                                  :: kappa_eppley
  real                                  :: kappa_remin
  integer                               :: km_c
  real                                  :: kfe_bal
  real                                  :: kfe_des
  real                                  :: kfe_max_prime
  real                                  :: kfe_org
  real                                  :: martin_coeff
  real                                  :: mass_2_n
  real                                  :: n_2_p
  real                                  :: n_2_p_fix
  character(len=fm_field_name_len)      :: name
  real                                  :: o2_min
  real                                  :: o_2_c
  real                                  :: o_2_no3
  real                                  :: o_2_nh4
  real                                  :: o_2_nitrif
  real, _ALLOCATABLE, dimension(:,:)    :: pco2surf  _NULL
  real                                  :: phi_wet
  real                                  :: phi_dry
  real                                  :: phi_don
  real                                  :: phi_dop
  real                                  :: phi_ldoc
  real                                  :: Prodstar
  real, _ALLOCATABLE, dimension(:)      :: r_1plusintzscale_si  _NULL
  real, _ALLOCATABLE, dimension(:)      :: r_1plusintzscale_ca  _NULL
  real, _ALLOCATABLE, dimension(:)      :: r_1plusintzscale_n  _NULL
  real                                  :: r_bio_tau
  real                                  :: r_bio_tau_don
  real                                  :: r_bio_tau_dop
  real                                  :: r_bio_tau_fix
  real                                  :: r_bio_tau_ldoc
  real                                  :: r_bio_tau_nh4
  real                                  :: r_bio_tau_nitrif_d
  real                                  :: r_bio_tau_nitrif_s
  type(mask_region_type)                :: r_bio_tau_prod
  type(mask_region_type)                :: nut_depl
  type(mask_region_type)                :: norm_remin
  type(mask_region_type)                :: no_caco3
  real, _ALLOCATABLE, dimension(:)      :: r_intzscale_n  _NULL
  real                                  :: r_wsink
  logical                               :: remin_density
  logical                               :: remin_lability
  logical                               :: remin_ocmip2
  logical                               :: remin_protection
  logical                               :: remin_simple
  logical                               :: remin_temp
  logical                               :: remin_viscosity
  logical                               :: remin_zoop_resp
  real                                  :: rpcaco3
  real                                  :: rpsio2
  logical                               :: soft_tissue_pump
  real                                  :: stp_alkalinity
  real                                  :: stp_salinity
  real                                  :: stp_temperature
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
  real                                  :: wsink
  real, _ALLOCATABLE, dimension(:)      :: zforg  _NULL

end type biotic_type

logical, public :: do_ocean_bgc_restore

integer                                 :: indsal
integer                                 :: indtemp
integer                                 :: package_index
logical                                 :: module_initialized = .false.

character(len=128) :: version = '$Id: ocean_bgc_restore.F90,v 20.0 2013/12/14 00:09:28 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

! Input parameters:
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
!
!  global_wrk_duration  = total time during calculation of global
!                         variables
character*128                                   :: alk_star_file    
integer                                         :: alk_star_id
character*128                                   :: alk_star_name
real, allocatable, dimension(:,:,:)             :: alk_star_t
integer                                         :: dep_wet_id
character*128                                   :: dep_wet_file    
character*128                                   :: dep_wet_name
real, allocatable, dimension(:,:)               :: dep_wet_t
integer                                         :: dep_dry_id
character*128                                   :: dep_dry_file    
character*128                                   :: dep_dry_name
real, allocatable, dimension(:,:)               :: dep_dry_t
character*128                                   :: fed_star_file    
integer                                         :: fed_star_id
character*128                                   :: fed_star_name
real, allocatable, dimension(:,:,:)             :: fed_star_t
integer                                         :: id_o2_sat
integer                                         :: km_c_max
character*128                                   :: no3_star_file    
integer                                         :: no3_star_id
character*128                                   :: no3_star_name
real, allocatable, dimension(:,:,:)             :: no3_star_t
character*128                                   :: po4_star_file    
integer                                         :: po4_star_id
character*128                                   :: po4_star_name
real, allocatable, dimension(:,:,:)             :: po4_star_t
real, allocatable, dimension(:,:)               :: sc_no_term
character*128                                   :: sio4_star_file    
integer                                         :: sio4_star_id
character*128                                   :: sio4_star_name
real, allocatable, dimension(:,:,:)             :: sio4_star_t
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
integer                                  :: num_restart = 0
type(restart_file_type), allocatable     :: restart(:)

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

allocate( dep_wet_t(isd:ied,jsd:jed) )
allocate( dep_dry_t(isd:ied,jsd:jed) )
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
allocate( alk_star_t(isd:ied,jsd:jed,nk) )
allocate( fed_star_t(isd:ied,jsd:jed,nk) )
allocate( no3_star_t(isd:ied,jsd:jed,nk) )
allocate( po4_star_t(isd:ied,jsd:jed,nk) )
allocate( sio4_star_t(isd:ied,jsd:jed,nk) )

! initialize some arrays
dep_wet_t(:,:) = 0.0
dep_dry_t(:,:) = 0.0
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
alk_star_t(:,:,:) = 0.0
fed_star_t(:,:,:) = 0.0
no3_star_t(:,:,:) = 0.0
po4_star_t(:,:,:) = 0.0
sio4_star_t(:,:,:) = 0.0

! allocate biotic array elements
do n = 1, instances

  allocate( biotic(n)%sc_co2(isc:iec,jsc:jec) )
  allocate( biotic(n)%sc_o2(isc:iec,jsc:jec) )
  allocate( biotic(n)%htotal(isc:iec,jsc:jec) )
  allocate( biotic(n)%csurf(isc:iec,jsc:jec) )
  allocate( biotic(n)%alpha(isc:iec,jsc:jec) )
  allocate( biotic(n)%comp_depth_frac(isc:iec,jsc:jec) )
  allocate( biotic(n)%pco2surf(isc:iec,jsc:jec) )
  allocate( biotic(n)%r_1plusintzscale_si(nk) )
  allocate( biotic(n)%r_1plusintzscale_ca(nk) )
  allocate( biotic(n)%r_1plusintzscale_n(nk) )
  allocate( biotic(n)%r_intzscale_n(nk) )
  allocate( biotic(n)%jfe_ads(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jfe_des(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jfe_graz(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jfe_sink(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jnh4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jnh4_graz(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jno3(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpo4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpo4_graz(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpofe(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpon(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpop(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_alk(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_fed(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_n_fix(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_nh4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_no3(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_n_norm(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_p_norm(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_p_fix(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_po4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_pofe(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_pon(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_pop(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jprod_sio4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jsio4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jdenit(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jdon(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jdop(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jldoc(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jo2(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jcaco3(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fpon(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fpop(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fracl(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fsio2(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fcaco3(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%flux_pon(isc:iec,jsc:jec) )
  allocate( biotic(n)%flux_pop(isc:iec,jsc:jec) )
  allocate( biotic(n)%flux_sio2(isc:iec,jsc:jec) )
  allocate( biotic(n)%flux_caco3(isc:iec,jsc:jec) )
  allocate( biotic(n)%nut_depl%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%no_caco3%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%norm_remin%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%r_bio_tau_prod%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%zforg(nk) )

enddo

! initialize bgc_restore array elements
do n = 1, instances

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%sc_co2(i,j)           = 0.0
      biotic(n)%sc_o2(i,j)            = 0.0
      biotic(n)%alpha(i,j)            = 0.0
      biotic(n)%comp_depth_frac(i,j)  = 0.0
      biotic(n)%pco2surf(i,j)         = 0.0
      biotic(n)%csurf(i,j)            = 0.0
      biotic(n)%htotal(i,j)           = 0.0
      do k = 1, nk
        biotic(n)%jprod_alk(i,j,k)    = 0.0
        biotic(n)%jprod_fed(i,j,k)    = 0.0
        biotic(n)%jprod_n_fix(i,j,k)  = 0.0
        biotic(n)%jprod_nh4(i,j,k)    = 0.0
        biotic(n)%jprod_no3(i,j,k)    = 0.0
        biotic(n)%jprod_n_norm(i,j,k) = 0.0
        biotic(n)%jprod_p_norm(i,j,k) = 0.0
        biotic(n)%jprod_p_fix(i,j,k)  = 0.0
        biotic(n)%jprod_po4(i,j,k)    = 0.0
        biotic(n)%jprod_pofe(i,j,k)   = 0.0
        biotic(n)%jprod_pon(i,j,k)    = 0.0
        biotic(n)%jprod_pop(i,j,k)    = 0.0
        biotic(n)%jprod_sio4(i,j,k)   = 0.0
        biotic(n)%jfe_ads(i,j,k)      = 0.0
        biotic(n)%jfe_des(i,j,k)      = 0.0
        biotic(n)%jfe_graz(i,j,k)     = 0.0
        biotic(n)%jfe_sink(i,j,k)     = 0.0
        biotic(n)%jpo4(i,j,k)         = 0.0
        biotic(n)%jpofe(i,j,k)        = 0.0
        biotic(n)%jpon(i,j,k)         = 0.0
        biotic(n)%jpop(i,j,k)         = 0.0
        biotic(n)%jsio4(i,j,k)        = 0.0
        biotic(n)%jdenit(i,j,k)       = 0.0
        biotic(n)%jdon(i,j,k)         = 0.0
        biotic(n)%jdop(i,j,k)         = 0.0
        biotic(n)%jldoc(i,j,k)        = 0.0
        biotic(n)%jo2(i,j,k)          = 0.0
        biotic(n)%jcaco3(i,j,k)       = 0.0
        biotic(n)%fpon(i,j,k)         = 0.0
        biotic(n)%fpop(i,j,k)         = 0.0
        biotic(n)%fracl(i,j,k)        = 0.0
        biotic(n)%fsio2(i,j,k)        = 0.0
        biotic(n)%fcaco3(i,j,k)       = 0.0
      enddo
    enddo
  enddo

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
! <SUBROUTINE NAME="ocean_bgc_restore_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>

subroutine ocean_bgc_restore_bbc(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, grid_kmt)

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
  if(biotic(n)%remin_ocmip2) then
    do j = jsc, jec
      do i = isc, iec
        kz = grid_kmt(i,j)
        if (kz .le. biotic(n)%km_c .and. kz .gt. 0) then
          t_prog(biotic(n)%ind_po4)%btf(i,j) =                    &
             t_prog(biotic(n)%ind_po4)%btf(i,j) -               &
             biotic(n)%flux_pop(i,j)
          t_prog(biotic(n)%ind_nh4)%btf(i,j) =                    &
             t_prog(biotic(n)%ind_nh4)%btf(i,j) -               &
             biotic(n)%flux_pon(i,j)
          t_prog(biotic(n)%ind_sio4)%btf(i,j) =                   &
             t_prog(biotic(n)%ind_sio4)%btf(i,j) -              &
             biotic(n)%flux_sio2(i,j)
          t_prog(biotic(n)%ind_o2)%btf(i,j)  =                    &
             t_prog(biotic(n)%ind_o2)%btf(i,j) +                &
             biotic(n)%o_2_nh4 * biotic(n)%flux_pon(i,j)
          t_prog(biotic(n)%ind_dic)%btf(i,j) =                    &
             t_prog(biotic(n)%ind_dic)%btf(i,j) +               &
             biotic(n)%flux_caco3(i,j) - biotic(n)%c_2_n *      &
             biotic(n)%flux_pon(i,j)
          t_prog(biotic(n)%ind_alk)%btf(i,j) =                    &
             t_prog(biotic(n)%ind_alk)%btf(i,j) +               &
             biotic(n)%flux_pon(i,j) - 2.0 * biotic(n)%flux_caco3(i,j)
        endif
      enddo
    enddo
  endif
enddo

return

end subroutine  ocean_bgc_restore_bbc
! </SUBROUTINE> NAME="ocean_bgc_restore_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>

subroutine ocean_bgc_restore_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     T_prog, T_diag, grid_dat, grid_tmask, mpp_domain2d, rho_dzt, taup1)

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
type(ocean_diag_tracer_type), intent(in), dimension(:)  :: T_diag
integer, intent(in)                                     :: taup1
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocean_bgc_restore_end'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer                                 :: i, j, k, n
integer                                 :: lun
character(len=fm_field_name_len+1)      :: suffix
real                                    :: total_alkalinity
real                                    :: total_ammonia
real                                    :: total_dic
real                                    :: total_don
real                                    :: total_dop
real                                    :: total_fediss
real                                    :: total_fepart
real                                    :: total_ldoc
real                                    :: total_o2
real                                    :: total_nitrate
real                                    :: total_phosphate
real                                    :: total_silicate

  integer :: stdoutunit 
  stdoutunit=stdout() 

!       integrate the total concentrations of some tracers
!       for the end of the run

!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
do n = 1, instances
  total_alkalinity = 0.0
  total_ammonia = 0.0
  total_dic = 0.0
  total_don = 0.0
  total_dop = 0.0
  total_fediss = 0.0
  total_fepart = 0.0
  total_ldoc = 0.0
  total_nitrate = 0.0
  total_o2 = 0.0
  total_phosphate = 0.0
  total_silicate = 0.0
  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_nitrate = total_nitrate +                         &
             t_prog(biotic(n)%ind_no3)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_ammonia = total_ammonia +                         &
             t_prog(biotic(n)%ind_nh4)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_phosphate = total_phosphate +                     &
             t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_fediss = total_fediss +                           &
             t_prog(biotic(n)%ind_fed)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_fepart = total_fepart +                           &
             t_diag(biotic(n)%ind_fep)%field(i,j,k) *           &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_silicate = total_silicate +                       &
             t_prog(biotic(n)%ind_sio4)%field(i,j,k,taup1) *    &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_don = total_don +                                 &
             t_prog(biotic(n)%ind_don)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_dop = total_dop +                                 &
             t_prog(biotic(n)%ind_dop)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_ldoc = total_ldoc +                               &
             t_prog(biotic(n)%ind_ldoc)%field(i,j,k,taup1) *    &
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

  call mpp_sum(total_nitrate)
  call mpp_sum(total_ammonia)
  call mpp_sum(total_phosphate)
  call mpp_sum(total_fediss)
  call mpp_sum(total_fepart)
  call mpp_sum(total_silicate)
  call mpp_sum(total_don)
  call mpp_sum(total_dop)
  call mpp_sum(total_ldoc)
  call mpp_sum(total_o2)
  call mpp_sum(total_dic)
  call mpp_sum(total_alkalinity)

  write (stdoutunit,*) '  Instance ', trim(biotic(n)%name)
  write (stdoutunit,                                              &
         '(/'' Total nitrate  = '',es19.12,'' Gmol-N'')')       &
              total_nitrate * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total ammonia  = '',es19.12,'' Gmol-N'')')       &
              total_ammonia * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total phosphate  = '',es19.12,'' Gmol-P'')')     &
              total_phosphate * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total fediss  = '',es19.12,'' Gmol-Fe'')')       &
              total_fediss * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total fepart  = '',es19.12,'' Gmol-Fe'')')       &
              total_fepart * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total silicate  = '',es19.12,'' Gmol-Si'')')     &
              total_silicate * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total DON  = '',es19.12,'' Gmol-C'')')           &
              total_DON * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total DOP  = '',es19.12,'' Gmol-P'')')           &
              total_DOP * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total LDOC  = '',es19.12,'' Gmol-C'')')          &
              total_LDOC * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total O2  = '',es19.12,'' Gmol-O'')')            &
              total_o2 * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total DIC  = '',es19.12,'' Gmol-C'')')           &
              total_dic * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total alkalinity  = '',es19.12,'' Geq'')')       &
              total_alkalinity * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total nitrogen  = '',es19.12,'' Gmol-N'')')      &
              (total_nitrate + total_don) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total phosphorus  = '',es19.12,'' Gmol-P'')')    &
              (total_phosphate + total_dop) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total real O2  = '',es19.12,'' Gmol-O'')')       &
              (total_o2 + biotic(n)%o_2_no3 * total_nitrate) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total Carbon  = '',es19.12,'' Gmol-C'')')        &
              (total_dic + biotic(n)%c_2_n * total_don + total_ldoc) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total real alkalinity  = '',es19.12,'' Geq'')')  &
              (total_alkalinity + total_nitrate) * 1.0e-09
enddo

! save out additional information for a restart
write(stdoutunit,*)

call ocean_bgc_restore_restart
do n = 1, instances

  write(stdoutunit,*) trim(note_header),                          &
       'Writing additional restart information for instance ',  &
       trim(biotic(n)%name)

  write (stdoutunit,*) trim(note_header),                         &
       'Done writing additional restart information for instance ',&
       trim(biotic(n)%name)

enddo

return
end subroutine  ocean_bgc_restore_end
! </SUBROUTINE> NAME="ocean_bgc_restore_end"

!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_bgc_restore_restart(time_stamp)
  character(len=*),  intent(in), optional :: time_stamp
  integer :: n

  do n=1, num_restart
     call save_restart(restart(n), time_stamp)
  end do

end subroutine ocean_bgc_restore_restart
! </SUBROUTINE> NAME="ocean_bgc_restore_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
subroutine ocean_bgc_restore_sbc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     T_prog, tau, Time, Grid, ice_ocean_boundary_fluxes)

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
integer, intent(in)                                             :: tau
type(ocean_time_type), intent(in)                               :: Time
type(ocean_grid_type), intent(in)                               :: Grid
type(coupler_2d_bc_type), intent(in)                            :: ice_ocean_boundary_fluxes

integer :: i, j, n
integer :: i_bnd_off
integer :: j_bnd_off
logical :: used

!     calculate interpolated iron deposition
call time_interp_external(dep_dry_id, Time%model_time, dep_dry_t)

! use the surface fluxes from the coupler
! stf is in mol/m^2/s, flux from coupler is positive upwards
i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances
  do j = jsc, jec
    do i = isc, iec
      t_prog(biotic(n)%ind_dic)%stf(i,j) =                      &
            -ice_ocean_boundary_fluxes%bc(biotic(n)%ind_co2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
      t_prog(biotic(n)%ind_o2)%stf(i,j) =                       &
            -ice_ocean_boundary_fluxes%bc(biotic(n)%ind_o2_flux)%field(ind_flux)%values(i-i_bnd_off,j-j_bnd_off)
    enddo
  enddo

!     surface iron deposition fluxes from Ginoux et al (JGR, Sources 
!     and distributions of dust aerosols simulated with the GOCART 
!     model, 2001) and converted by Gregg et al (GRL, Ocean primary
!     producion and climate: Global decadal changes, in press) assuming a 2%
!     solubility of Fe in dust.
!          For additional source from extraterrestrial dust after Johnson
!     (GBC, Iron supply and demand in the upper ocean: Is extraterrestrial
!     dust a significant source of bioavailable iron?), add a constant flux
!     of 9.51e-15.
  do j = jsc, jec
    do i = isc, iec
      t_prog(biotic(n)%ind_fed)%stf(i,j) =  dep_dry_t(i,j)
    enddo
  enddo

enddo

! Save variables for diagnostics
do n = 1, instances
   call diagnose_2d(Time, Grid, biotic(n)%id_sfc_flux_co2, t_prog(biotic(n)%ind_dic)%stf(:,:))
   call diagnose_2d(Time, Grid, biotic(n)%id_sfc_flux_o2, t_prog(biotic(n)%ind_o2)%stf(:,:))
   call diagnose_2d(Time, Grid, biotic(n)%id_sfc_flux_fed, t_prog(biotic(n)%ind_fed)%stf(:,:))
enddo

return

end subroutine  ocean_bgc_restore_sbc
! </SUBROUTINE> NAME="ocean_bgc_restore_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocean_bgc_restore_flux_init

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux

character(len=64), parameter    :: sub_name = 'ocean_bgc_restore_flux_init'
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

  ! First, perform some initialization if this module has not been
  ! initialized because the normal initialization routine will
  ! not have been called as part of the normal ocean model
  ! initialization if this is an Atmosphere pe of a coupled
  ! model running in concurrent mode

if (.not. module_initialized) then 

   ! Initialize the package
  package_index = otpm_set_tracer_package(package_name,            &
       restart_file = default_restart_file,                             &
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
    do_ocean_bgc_restore = .false.
  else
    if (instances .eq. 1) then
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif
    do_ocean_bgc_restore = .true.
  endif

  module_initialized = .true.

endif

! Return if we don't want to use this package
if (.not. do_ocean_bgc_restore) then
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

end subroutine  ocean_bgc_restore_flux_init
!</SUBROUTINE> NAME="ocean_bgc_restore_flux_init"

!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocean_bgc_restore_init

character(len=64), parameter    :: sub_name = 'ocean_bgc_restore_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

real, parameter :: rho_avg = 1024.5
real, parameter :: sperd = 24.0 * 3600.0

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

! Initialize the restoring package
package_index = otpm_set_tracer_package(package_name,            &
     restart_file = default_restart_file,                             &
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
  do_ocean_bgc_restore = .false.
else
  if (instances .eq. 1) then
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
  else
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  endif
  do_ocean_bgc_restore = .true.
endif

module_initialized = .true.

! Return if we don't want to use this package,
! after changing the list back
if (.not. do_ocean_bgc_restore) then
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

  ! NO3
  biotic(n)%ind_no3 = otpm_set_prog_tracer('no3' // suffix, package_name,   &
       longname = 'Nitrate' // trim(long_suffix),                      &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! NH4
  biotic(n)%ind_nh4 = otpm_set_prog_tracer('nh4' // suffix, package_name,   &
       longname = 'Ammonia' // '(' // trim(name) // ')',               &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! PO4
  biotic(n)%ind_po4 = otpm_set_prog_tracer('po4' // suffix, package_name,   &
       longname = 'Phosphate' // trim(long_suffix),             &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! Dissolved Fe (assumed to be all available to phytoplankton)
  biotic(n)%ind_fed = otpm_set_prog_tracer('fed' // suffix, package_name,   &
       longname = 'Dissolved Iron' // '(' // trim(name) // ')',        &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! Fep (Sinking detrital/particulate iron)
  biotic(n)%ind_fep = otpm_set_diag_tracer('fep' // suffix, package_name,   &
       longname = 'Particulate Iron' // '(' // trim(name) // ')',      &
       restart_file = default_restart_file, units = 'mol/kg')

  ! SiO4
  biotic(n)%ind_sio4 = otpm_set_prog_tracer('sio4' // suffix, package_name, &
       longname = 'Silicate' // '(' // trim(name) // ')',              &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! DON
  biotic(n)%ind_don = otpm_set_prog_tracer('don' // suffix, package_name,   &
       longname = 'DON' // '(' // trim(name) // ')',                   &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! DOP
  biotic(n)%ind_dop = otpm_set_prog_tracer('dop' // suffix, package_name,   &
       longname = 'DOP' // '(' // trim(name) // ')',                   &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! LDOC
  biotic(n)%ind_ldoc = otpm_set_prog_tracer('ldoc' // suffix, package_name, &
       longname = 'labile DOC' // '(' // trim(name) // ')',            &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! DIC
  biotic(n)%ind_dic = otpm_set_prog_tracer('dic' // suffix, package_name,   &
       longname = 'DIC' // '(' // trim(name) // ')',                   &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! O2
  biotic(n)%ind_o2 = otpm_set_prog_tracer('o2' // suffix, package_name,     &
       longname = 'Oxygen' // '(' // trim(name) // ')',                &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
       caller = caller_str)

  ! ALK (Total carbonate alkalinity)
  biotic(n)%ind_alk = otpm_set_prog_tracer('alk' // suffix, package_name,   &
       longname = 'Alkalinity' // '(' // trim(name) // ')',            &
       units = 'mol/kg', flux_units = 'mol/m^2/s',             &
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

call fm_util_set_value('dep_dry_file', 'INPUT/fe_dep_bc.nc')
call fm_util_set_value('dep_dry_name', 'fe_dep')
call fm_util_set_value('dep_wet_file', 'INPUT/fe_dep_bc.nc')
call fm_util_set_value('dep_wet_name', 'fe_dep')
call fm_util_set_value('no3_star_file', 'INPUT/no3_star.nc')
call fm_util_set_value('no3_star_name', 'no3_star')
call fm_util_set_value('po4_star_file', 'INPUT/po4_star.nc')
call fm_util_set_value('po4_star_name', 'po4_star')
call fm_util_set_value('sio4_star_file', 'INPUT/sio4_star.nc')
call fm_util_set_value('sio4_star_name', 'sio4_star')
call fm_util_set_value('alk_star_file', 'INPUT/alk_star.nc')
call fm_util_set_value('alk_star_name', 'alk_star')
call fm_util_set_value('fed_star_file', 'INPUT/fed_star.nc')
call fm_util_set_value('fed_star_name', 'fed_star')
call fm_util_set_value('htotal_scale_lo_in', 0.01 )    ! scale
call fm_util_set_value('htotal_scale_hi_in', 100.0)    ! scale
call fm_util_set_value('htotal_in', 1.0e-08)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

! Set up the instance namelists

t_mask(:) = .true.

do n = 1, instances

  ! create the instance namelist
  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('compensation_depth', 75.0)                  ! m
  call fm_util_set_value('martin_coeff', 0.9)
  call fm_util_set_value('ca_remin_depth', 3500.0)
  call fm_util_set_value('si_remin_depth', 1000.0)
  call fm_util_set_value('soft_tissue_pump', .false.)
  call fm_util_set_value('stp_temperature', 10.0)
  call fm_util_set_value('stp_salinity', 34.7)
  ! alkalinity is in ueq/kg, converted to eq/m^3
  call fm_util_set_value('stp_alkalinity', 2370.0 * 1024.5 * 1.0e-06)
  call fm_util_set_value('local_restart_file', default_local_restart_file)
  call fm_util_set_value('fe_ballast_assoc', .false.)
  call fm_util_set_value('phi_dry', 0.03)
  call fm_util_set_value('phi_wet', 0.03)
  call fm_util_set_value('kfe_org', 0.005/sperd)
  call fm_util_set_value('kfe_bal', 0.005/sperd)
  call fm_util_set_value('kfe_des', 0.0/sperd)
  call fm_util_set_value('kfe_max_prime', 1.0/sperd)
  call fm_util_set_value('mass_2_n', 117.0/16.0*12.0*1.87)
  call fm_util_set_value('n_2_p', 16.0)
  call fm_util_set_value('n_2_p_fix', 40.0)
  call fm_util_set_value('c_2_n', 117.0/16.0)
  call fm_util_set_value('o_2_c', 170.0/16.0)
  call fm_util_set_value('o_2_no3', 255.0/16.0)
  call fm_util_set_value('o_2_nh4', 207.0/16.0)
  call fm_util_set_value('o_2_nitrif', 3.0)
  call fm_util_set_value('o2_min', 4.0 * rho_avg * 1.0e-06)
  call fm_util_set_value('bio_tau', 30.0 * sperd)
  call fm_util_set_value('bio_tau_don', 30.0*365.0*sperd)
  call fm_util_set_value('bio_tau_dop', 30.0*365.0*sperd)
  call fm_util_set_value('bio_tau_fix', 90.0*sperd)
  call fm_util_set_value('bio_tau_ldoc', 0.5*365.0*sperd)
  call fm_util_set_value('bio_tau_nh4', 3.0*sperd)
  call fm_util_set_value('bio_tau_nitrif_d', 10.0*sperd)
  call fm_util_set_value('bio_tau_nitrif_s', 1e20*sperd)
  call fm_util_set_value('kappa_eppley', 0.063)
  call fm_util_set_value('kappa_remin', -0.032)
  call fm_util_set_value('Prodstar', 0.37/1000*16/117/sperd)
  call fm_util_set_value('fdets0', 0.14)
  call fm_util_set_value('fdetl0', 0.74)
  call fm_util_set_value('phi_don', 0.04)
  call fm_util_set_value('phi_dop', 0.04)
  call fm_util_set_value('phi_ldoc', 0.04)
  call fm_util_set_value('gamma_det', 0.0116/sperd)
  call fm_util_set_value('remin_density', .false.)
  call fm_util_set_value('remin_lability', .false.)
  call fm_util_set_value('remin_ocmip2', .false.)
  call fm_util_set_value('remin_protection', .true.)
  call fm_util_set_value('remin_simple', .false.)
  call fm_util_set_value('remin_temp', .false.)
  call fm_util_set_value('remin_viscosity', .false.)
  call fm_util_set_value('remin_zoop_resp', .false.)
  call fm_util_set_value('rpcaco3', 0.070/12*16/117*100)
  call fm_util_set_value('rpsio2', 0.026/12*16/117*60)
  call fm_util_set_value('wsink', 3.0/sperd)
  ! Wanninkhof numbers
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

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_prod',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call fm_util_set_value('factor', 0.0)
  call fm_util_set_value('coastal_only', .false.)
  call fm_util_set_value('t_mask', t_mask, size(t_mask))
  call fm_util_set_value('wlon', 0.0, index = 0)
  call fm_util_set_value('elon', 0.0, index = 0)
  call fm_util_set_value('slat', 0.0, index = 0)
  call fm_util_set_value('nlat', 0.0, index = 0)

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_prod', caller = caller_str)

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

end subroutine ocean_bgc_restore_init
! </SUBROUTINE> NAME="ocean_bgc_restore_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocean_bgc_restore_start
! </DESCRIPTION>

subroutine ocean_bgc_restore_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,       &
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
  ind = biotic(n)%ind_co2_flux
  if (.not. field_exist('INPUT/'//trim(ocean_fields%bc(ind)%ocean_restart_file),   &
                        ocean_fields%bc(ind)%field(ind_alpha)%name)) then

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,         &
         grid_tmask(isd:ied,jsd:jed,1),                                 &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
         t_prog(biotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_alk)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_po4)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_sio4)%field(isd:ied,jsd:jed,1,taum1),     &
         htotal_scale_lo, htotal_scale_hi, biotic(n)%htotal,            &
         co2star = biotic(n)%csurf, alpha = biotic(n)%alpha,            &
         pco2surf = biotic(n)%pco2surf)

    ! Compute the Schmidt number of CO2 in seawater using the 
    ! formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    ! 7373-7382).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_co2(i,j) =                                                 &
             biotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *          &
             (biotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *         &
              (biotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *        &
               biotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1)
                                              
        ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) = &
             biotic(n)%alpha(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
        ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) = &
             biotic(n)%csurf(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo

  endif

  ! O2 flux
  ind = biotic(n)%ind_o2_flux
  if (.not. field_exist('INPUT/'//trim(ocean_fields%bc(ind)%ocean_restart_file),    &
                        ocean_fields%bc(ind)%field(ind_alpha)%name)) then

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
                  c_0*t_prog(indsal)%field(i,j,1,taum1))) * grid_tmask(i,j,1)
      enddo
    enddo

    ! convert from ml/l to mol/m^3
    do j = jsc, jec
      do i = isc, iec
        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
      enddo
    enddo

    ! Compute the Schmidt number of O2 in seawater using the 
    ! formulation proposed by Keeling et al. (1998, Global Biogeochem.
    ! Cycles, 12, 141-163).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_o2(i,j) =                                                  &
             biotic(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *           &
             (biotic(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *          &
              (biotic(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *         &
               biotic(n)%sc_o2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_o2(i,j) + epsln)) * grid_tmask(i,j,1)
        ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             o2_saturation(i,j) * sc_no_term(i,j)
        ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             t_prog(biotic(n)%ind_o2)%field(i,j,1,taum1) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo

  endif
enddo

return

end subroutine ocean_bgc_restore_init_sfc
! </SUBROUTINE> NAME="ocean_bgc_restore_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_bgc_restore_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
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

integer :: i
integer :: i_bnd_off
integer :: j_bnd_off
integer :: j
integer :: n
integer :: nn
integer :: ind

real    :: epsln=1.0e-30

i_bnd_off = isc - isc_bnd
j_bnd_off = jsc - jsc_bnd

do n = 1, instances

    ind = biotic(n)%ind_co2_flux

    call ocmip2_co2calc(isd, jsd, isc, iec, jsc, jec,         &
         grid_tmask(isd:ied,jsd:jed,1),                                 &
         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
         t_prog(biotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_alk)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_po4)%field(isd:ied,jsd:jed,1,taum1),      &
         t_prog(biotic(n)%ind_sio4)%field(isd:ied,jsd:jed,1,taum1),     &
         htotal_scale_lo, htotal_scale_hi, biotic(n)%htotal,            &
         co2star = biotic(n)%csurf, alpha = biotic(n)%alpha,            &
         pco2surf = biotic(n)%pco2surf)

    ! Compute the Schmidt number of CO2 in seawater using the 
    ! formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    ! 7373-7382).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_co2(i,j) =                                                 &
             biotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *          &
             (biotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *         &
              (biotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *        &
               biotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_co2(i,j) + epsln)) * grid_tmask(i,j,1) 
        ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +    &
             biotic(n)%alpha(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
        ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +    &
             biotic(n)%csurf(i,j) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo

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
                  c_0*t_prog(indsal)%field(i,j,1,taum1))) * grid_tmask(i,j,1) 
      enddo
    enddo

    ! convert from ml/l to mol/m^3
    do j = jsc, jec
      do i = isc, iec
        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
      enddo
    enddo

    ! Compute the Schmidt number of O2 in seawater using the 
    ! formulation proposed by Keeling et al. (1998, Global Biogeochem.
    ! Cycles, 12, 141-163).
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%sc_o2(i,j) =                                                  &
             biotic(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *           &
             (biotic(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *          &
              (biotic(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *         &
               biotic(n)%sc_o2_3)) * grid_tmask(i,j,1)
        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_o2(i,j) + epsln)) * grid_tmask(i,j,1) 
        ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =         &
             ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) +    &
             o2_saturation(i,j) * sc_no_term(i,j)
        ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =         &
             ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) +    &
             t_prog(biotic(n)%ind_o2)%field(i,j,1,taum1) * sc_no_term(i,j) * rho(i,j,1,taum1)
      enddo
    enddo

enddo

return

end subroutine ocean_bgc_restore_sum_sfc
! </SUBROUTINE> NAME="ocean_bgc_restore_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_bgc_restore_zero_sfc(Ocean_fields)

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

integer         :: n
integer         :: ind

do n = 1, instances

  ind = biotic(n)%ind_co2_flux

  ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

  ind = biotic(n)%ind_o2_flux

  ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
  ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

enddo

return

end subroutine ocean_bgc_restore_zero_sfc
! </SUBROUTINE> NAME="ocean_bgc_restore_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine ocean_bgc_restore_avg_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
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

divid = 1./float(ocean_avg_kount)

do n = 1, instances

  ind = biotic(n)%ind_co2_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo

  ind = biotic(n)%ind_o2_flux

  do j = jsc, jec
    do i = isc, iec
      if (grid_tmask(i,j,1) == 1.0) then
        ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             ocean_fields%bc(ind)%field(ind_alpha)%values(i-i_bnd_off,j-j_bnd_off) * divid
        ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) =                 &
             ocean_fields%bc(ind)%field(ind_csurf)%values(i-i_bnd_off,j-j_bnd_off) * divid
      endif
    enddo
  enddo

enddo

return

end subroutine ocean_bgc_restore_avg_sfc
! </SUBROUTINE> NAME="ocean_bgc_restore_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_bgc_restore_sfc_end

end subroutine ocean_bgc_restore_sfc_end
! </SUBROUTINE> NAME="ocean_bgc_restore_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>

subroutine ocean_bgc_restore_source(isc, iec, jsc, jec, nk, isd, ied, jsd, jed, &
     T_prog, T_diag, taum1, model_time, Grid, Time, grid_kmt, rho_dzt, dtts)

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
type(ocean_diag_tracer_type), intent(inout), dimension(:)       :: T_diag
integer, intent(in)                                             :: taum1
type(time_type), intent(in)                                     :: model_time
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt
real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho_dzt
real, intent(in)                                                :: dtts

integer :: i, j, k, kmax, n
integer :: ind_no3
integer :: ind_nh4
integer :: ind_po4
integer :: ind_fed
integer :: ind_fep
integer :: ind_sio4
integer :: ind_alk
integer :: ind_ldoc
integer :: ind_don
integer :: ind_dop
integer :: ind_dic
integer :: ind_o2
integer :: km_c
logical :: used
integer :: day
integer :: month
integer :: year
integer :: hour
integer :: minute
integer :: second

real :: jtot
real :: SoverPstar2
real :: expkT
real :: fpon_protected

  integer :: stdoutunit 
  stdoutunit=stdout() 


! get the model month
call get_date(model_time, year, month, day,                &
              hour, minute, second)

! calculate the source terms for BIOTICs

! calculate interpolated NO3_star, PO4_star, SiO4_star and Alk_star
call time_interp_external(no3_star_id, model_time, no3_star_t)
call time_interp_external(po4_star_id, model_time, po4_star_t)
call time_interp_external(fed_star_id, model_time, fed_star_t)
call time_interp_external(sio4_star_id, model_time, sio4_star_t)
call time_interp_external(alk_star_id, model_time, alk_star_t)

! Loop over multiple instances
do n = 1, instances

  ind_no3 = biotic(n)%ind_no3
  ind_nh4 = biotic(n)%ind_nh4
  ind_po4 = biotic(n)%ind_po4
  ind_fed = biotic(n)%ind_fed
  ind_fep = biotic(n)%ind_fep
  ind_sio4 = biotic(n)%ind_sio4
  ind_alk = biotic(n)%ind_alk
  ind_ldoc = biotic(n)%ind_ldoc
  ind_don = biotic(n)%ind_don
  ind_dop = biotic(n)%ind_dop
  ind_dic = biotic(n)%ind_dic
  ind_o2 = biotic(n)%ind_o2
  km_c = biotic(n)%km_c

  ! Production

  ! compute NO3 restoring term and correct for partial
  ! production in the bottom box
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        if (t_prog(ind_no3)%field(i,j,k,taum1) .gt.             &
            no3_star_t(i,j,k) *                                 &
            biotic(n)%nut_depl%mask(i,j,month)) then
          biotic(n)%jprod_no3(i,j,k) =                              &
               (t_prog(ind_no3)%field(i,j,k,taum1) -            &
                no3_star_t(i,j,k) *                             &
                biotic(n)%nut_depl%mask(i,j,month)) *           &
               biotic(n)%r_bio_tau_prod%mask(i,j,month) * Grid%tmask(i,j,k)
        else
          biotic(n)%jprod_no3(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jprod_no3(i,j,km_c) = biotic(n)%jprod_no3(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
    enddo
  enddo

  ! compute NH4 restoring term and correct for partial
  ! production in the bottom box
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        if (t_prog(ind_nh4)%field(i,j,k,taum1) .gt. 0.0) then
          biotic(n)%jprod_nh4(i,j,k) =                              &
               t_prog(ind_nh4)%field(i,j,k,taum1) *           &
               biotic(n)%r_bio_tau_nh4 * Grid%tmask(i,j,k)
        else
          biotic(n)%jprod_nh4(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jprod_nh4(i,j,km_c) = biotic(n)%jprod_nh4(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
    enddo
  enddo

  ! compute PO4 restoring term and correct for partial
  ! production in the bottom box
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        if (t_prog(ind_po4)%field(i,j,k,taum1) .gt.             &
            po4_star_t(i,j,k) *                                 &
            biotic(n)%nut_depl%mask(i,j,month)) then
          biotic(n)%jprod_po4(i,j,k) =                              &
               (t_prog(ind_po4)%field(i,j,k,taum1) -            &
                po4_star_t(i,j,k) *                             &
                biotic(n)%nut_depl%mask(i,j,month)) *           &
               biotic(n)%r_bio_tau_prod%mask(i,j,month) * Grid%tmask(i,j,k)
        else
          biotic(n)%jprod_po4(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jprod_po4(i,j,km_c) = biotic(n)%jprod_po4(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
    enddo
  enddo

  ! compute Fed restoring term and correct for partial
  ! production in the bottom box
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        if (t_prog(ind_fed)%field(i,j,k,taum1) .gt.             &
            fed_star_t(i,j,k) *                                 &
            biotic(n)%nut_depl%mask(i,j,month)) then
          biotic(n)%jprod_fed(i,j,k) =                              &
               (t_prog(ind_fed)%field(i,j,k,taum1) -            &
                fed_star_t(i,j,k) *                             &
                biotic(n)%nut_depl%mask(i,j,month)) *           &
               biotic(n)%r_bio_tau_prod%mask(i,j,month) * Grid%tmask(i,j,k)
        else
          biotic(n)%jprod_fed(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec 
      biotic(n)%jprod_fed(i,j,km_c) = biotic(n)%jprod_fed(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
    enddo
  enddo


  ! compute nitrogen fixation
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
       biotic(n)%jprod_n_norm(i,j,k) = biotic(n)%jprod_nh4(i,j,k) +     &
         biotic(n)%jprod_no3(i,j,k)
       biotic(n)%jprod_p_norm(i,j,k) = biotic(n)%jprod_po4(i,j,k)
       biotic(n)%jprod_p_fix(i,j,k) = 0.0
       biotic(n)%jprod_n_fix(i,j,k) = 0.0
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%jprod_p_norm(i,j,km_c) = biotic(n)%jprod_p_norm(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
      biotic(n)%jprod_p_fix(i,j,km_c) = biotic(n)%jprod_p_fix(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
      biotic(n)%jprod_n_fix(i,j,km_c) = biotic(n)%jprod_n_fix(i,j,km_c) *   &
           biotic(n)%comp_depth_frac(i,j)
    enddo 
  enddo

  ! compute SiO4 restoring term and correct for partial
  ! production in the bottom box
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
          if (t_prog(ind_sio4)%field(i,j,k,taum1) .gt.          &
              sio4_star_t(i,j,k)) then
            biotic(n)%jprod_sio4(i,j,k) =                       &
                 (t_prog(ind_sio4)%field(i,j,k,taum1) -         &
                  sio4_star_t(i,j,k)) *                         &
                 biotic(n)%r_bio_tau * Grid%tmask(i,j,k)
          else
            biotic(n)%jprod_sio4(i,j,k) = 0.0
        endif
      enddo
    enddo
  enddo

    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jprod_sio4(i,j,km_c) =                        &
             biotic(n)%jprod_sio4(i,j,km_c) *                   &
             biotic(n)%comp_depth_frac(i,j)
      enddo
    enddo

    ! compute Alk restoring term and correct for partial
    ! production in the bottom box
    do k = 1, km_c
      do j = jsc, jec
        do i = isc, iec
          if ((t_prog(ind_alk)%field(i,j,k,taum1) +            &
            t_prog(ind_no3)%field(i,j,k,taum1)) * 35.0 /       &
            (t_prog(indsal)%field(i,j,k,taum1) + 1e-40) .gt.             &
              alk_star_t(i,j,k)) then
            biotic(n)%jprod_alk(i,j,k) =                       &
                 ((t_prog(ind_alk)%field(i,j,k,taum1) +        &
            t_prog(ind_no3)%field(i,j,k,taum1)) * 35.0 /       &
            (t_prog(indsal)%field(i,j,k,taum1) + 1e-40) -         &
                  alk_star_t(i,j,k)) *                         &
                 biotic(n)%r_bio_tau * Grid%tmask(i,j,k)
          else
            biotic(n)%jprod_alk(i,j,k) = 0.0
          endif
        enddo
      enddo
    enddo

    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jprod_alk(i,j,km_c) =                        &
             biotic(n)%jprod_alk(i,j,km_c) *                   &
             biotic(n)%comp_depth_frac(i,j)
      enddo
    enddo

    ! Do not allow production in bottom level
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jprod_no3(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_nh4(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_po4(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_n_norm(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_p_norm(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_n_fix(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_p_fix(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_fed(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_sio4(i,j,grid_kmt(i,j)) = 0.0
        biotic(n)%jprod_alk(i,j,grid_kmt(i,j)) = 0.0
      enddo
    enddo

  ! Food Web Processing
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        jtot=max(biotic(n)%jprod_n_norm(i,j,k) +                                &
            biotic(n)%jprod_n_fix(i,j,k),1.0e-30)
          expkT = exp(biotic(n)%kappa_eppley * t_prog(indtemp)%field(i,j,k,taum1))
        SoverPstar2=-0.5 + 0.5 * sqrt(1.0 + 4.0 * jtot /                        &
            (expkT * biotic(n)%Prodstar))
        biotic(n)%fracl(i,j,k) = SoverPstar2 / (1.0 + SoverPstar2) *            &
            Grid%tmask(i,j,k)
        biotic(n)%jprod_pon(i,j,k) = max(0.0,jtot * (exp(biotic(n)%kappa_remin *&
             t_prog(indtemp)%field(i,j,k,taum1)) *                              &
            (biotic(n)%fdetl0 * biotic(n)%fracl(i,j,k) + biotic(n)%fdets0 *     &
            (1.0 - biotic(n)%fracl(i,j,k))) - biotic(n)%phi_don)) *             &
            Grid%tmask(i,j,k)
        biotic(n)%jprod_pop(i,j,k) = biotic(n)%jprod_pon(i,j,k) / jtot *        &
            (biotic(n)%jprod_p_norm(i,j,k) + biotic(n)%jprod_p_fix(i,j,k))
        biotic(n)%jprod_pofe(i,j,k) = biotic(n)%jprod_pon(i,j,k) / jtot *       &
            biotic(n)%jprod_fed(i,j,k)
        biotic(n)%jldoc(i,j,k) = jtot * biotic(n)%phi_ldoc * Grid%tmask(i,j,k)
        biotic(n)%jdon(i,j,k) = jtot * biotic(n)%phi_don * Grid%tmask(i,j,k)
        biotic(n)%jdop(i,j,k) = (biotic(n)%jprod_p_norm(i,j,k) +                &
            biotic(n)%jprod_p_fix(i,j,k)) * biotic(n)%phi_dop *                 &
            Grid%tmask(i,j,k)
        biotic(n)%jnh4_graz(i,j,k) = (jtot - biotic(n)%jprod_pon(i,j,k) -       &
            biotic(n)%jdon(i,j,k)) * Grid%tmask(i,j,k)
        biotic(n)%jpo4_graz(i,j,k) = (biotic(n)%jprod_p_norm(i,j,k) +           &
            biotic(n)%jprod_p_fix(i,j,k) - biotic(n)%jprod_pop(i,j,k) -         &
            biotic(n)%jdop(i,j,k)) * Grid%tmask(i,j,k)
        biotic(n)%jfe_graz(i,j,k) = (biotic(n)%jprod_fed(i,j,k)                &
            - biotic(n)%jprod_pofe(i,j,k)) * Grid%tmask(i,j,k)
      enddo
    enddo
  enddo

  ! Accumulate first level of flux
  do j=jsc,jec
    do i=isc,iec
      biotic(n)%fsio2(i,j,1) = biotic(n)%jprod_sio4(i,j,1) *        &
        rho_dzt(i,j,1,taum1) * Grid%tmask(i,j,1)
      biotic(n)%fcaco3(i,j,1) = 0.5 * biotic(n)%jprod_alk(i,j,1) *  &
        rho_dzt(i,j,1,taum1) * Grid%tmask(i,j,1)
      biotic(n)%fpon(i,j,1) = biotic(n)%jprod_pon(i,j,1) *          &
         rho_dzt(i,j,1,taum1) * Grid%tmask(i,j,1)
      biotic(n)%fpop(i,j,1) = biotic(n)%jprod_pop(i,j,1) *          &
         rho_dzt(i,j,1,taum1) * Grid%tmask(i,j,1)
      biotic(n)%jfe_ads(i,j,1)=min(biotic(n)%kfe_max_prime,         &
        (biotic(n)%kfe_org / 2.0 * biotic(n)%fpon(i,j,1) *          &
        biotic(n)%mass_2_n + biotic(n)%kfe_bal / 2.0 *              &
        (biotic(n)%fsio2(i,j,1) * 60.0 + biotic(n)%fcaco3(i,j,1) *  &
        100.0)) * biotic(n)%r_wsink) *                              &
        max(0.0,t_prog(ind_fed)%field(i,j,1,taum1))
      biotic(n)%jfe_des(i,j,1) = biotic(n)%kfe_des *                &
        max(0.0,t_diag(ind_fep)%field(i,j,1))
      biotic(n)%jpon(i,j,1) = 0.0
      biotic(n)%jpop(i,j,1) = 0.0
      biotic(n)%jpofe(i,j,1) = 0.0
      biotic(n)%jsio4(i,j,1) = 0.0
      biotic(n)%jcaco3(i,j,1) = 0.0
    enddo
  enddo

  do k=2,nk-1
    do j=jsc,jec
      do i=isc,iec
        biotic(n)%fsio2(i,j,k) = biotic(n)%fsio2(i,j,k-1) *                    &
          biotic(n)%r_1plusintzscale_si(k) * Grid%tmask(i,j,k)

        ! Calculate regeneration term
        biotic(n)%jsio4(i,j,k) = (biotic(n)%fsio2(i,j,k-1) * Grid%tmask(i,j,k) &
          - biotic(n)%fsio2(i,j,k) * Grid%tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)

        ! Add production within box to flux
        biotic(n)%fsio2(i,j,k) = biotic(n)%fsio2(i,j,k) +                      &
          biotic(n)%jprod_sio4(i,j,k) *                                        &
          rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
      enddo
    enddo
  enddo

  do k=2,nk-1
    do j=jsc,jec
      do i=isc,iec
        biotic(n)%fcaco3(i,j,k) = biotic(n)%fcaco3(i,j,k-1) *                  &
          biotic(n)%r_1plusintzscale_ca(k) * Grid%tmask(i,j,k)

        ! Calculate regeneration term
        biotic(n)%jcaco3(i,j,k) = (biotic(n)%fcaco3(i,j,k-1) * Grid%tmask(i,j,k) &
          - biotic(n)%fcaco3(i,j,k) * Grid%tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)

        ! Add production within box to flux
        biotic(n)%fcaco3(i,j,k) = biotic(n)%fcaco3(i,j,k) + 0.5 *              &
          biotic(n)%jprod_alk(i,j,k) *                                         &
          rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
        enddo
      enddo
    enddo

if(biotic(n)%remin_ocmip2) then
  ! OCMIP2 Interior Remineralization Scheme

  ! F(z) = F75 * (z / 75)^-0.9
  do k = 2, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%fpon(i,j,k) = biotic(n)%fpon(i,j,k-1) +       &
          biotic(n)%jprod_pon(i,j,k) *                          &
          rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
        biotic(n)%fpop(i,j,k) = biotic(n)%fpop(i,j,k-1) +       &
          biotic(n)%jprod_pop(i,j,k) *                          &
          rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
      enddo
    enddo
  enddo

  do j = jsc, jec
    do i = isc, iec
      biotic(n)%flux_pon(i,j) = biotic(n)%fpon(i,j,km_c)
      biotic(n)%flux_pop(i,j) = biotic(n)%fpop(i,j,km_c)
    enddo
  enddo

  do k = km_c, nk
    do j = jsc, jec
      do i = isc, iec
        ! Calculate the flux at the base of each layer below the
        ! compensation depth
        biotic(n)%fpon(i,j,k) = biotic(n)%flux_pon(i,j) *         &
             Grid%tmask(i,j,k) * biotic(n)%zforg(k)
        biotic(n)%fpop(i,j,k) = biotic(n)%flux_pop(i,j) *         &
             Grid%tmask(i,j,k) * biotic(n)%zforg(k)

        ! Calculate regeneration term
        biotic(n)%jpon(i,j,k) = (biotic(n)%fpon(i,j,k-1) * Grid%tmask(i,j,k)   &
          - biotic(n)%fpon(i,j,k) * Grid%tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)
        biotic(n)%jpop(i,j,k) = (biotic(n)%fpop(i,j,k-1) * Grid%tmask(i,j,k)   &
          - biotic(n)%fpop(i,j,k) * Grid%tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo
elseif(biotic(n)%remin_protection) then
  ! Ballast Protection Interior Remineralization Scheme
  ! remin = g * max(0.0 , Forg - rpcaco3 * Fcaco3 - rpsio2 * Fsio2)
  do k=2,nk-1
    do j=jsc,jec
      do i=isc,iec
         ! Remineralization of unprotected organic material and
         ! previously protected particulate organic material
        fpon_protected = biotic(n)%rpsio2 * biotic(n)%fsio2(i,j,k-1) +         &
          biotic(n)%rpcaco3 * biotic(n)%fcaco3(i,j,k-1) + 1.0e-30
        biotic(n)%fpon(i,j,k) = min(biotic(n)%fpon(i,j,k-1),                  &
          ((biotic(n)%fpon(i,j,k-1) + biotic(n)%r_intzscale_n(k) *              &
          fpon_protected) * biotic(n)%r_1plusintzscale_n(k)) -                 &
          (biotic(n)%rpsio2 * (biotic(n)%fsio2(i,j,k-1) -                      &
          biotic(n)%fsio2(i,j,k)) + biotic(n)%rpcaco3 *                        &
          (biotic(n)%fcaco3(i,j,k-1) - biotic(n)%fcaco3(i,j,k))) *             &
          min(biotic(n)%fpon(i,j,k-1),fpon_protected) /                        &
          fpon_protected * biotic(n)%r_1plusintzscale_n(k)) * Grid%tmask(i,j,k)

        ! Apply N change to P assuming equal partitioning between protected,
        ! previously protected and unprotected particulate organic material
        biotic(n)%fpop(i,j,k) = biotic(n)%fpon(i,j,k) /                        &
          max(biotic(n)%fpon(i,j,k-1),1e-30) * biotic(n)%fpop(i,j,k-1)

        ! Calculate regeneration term
        biotic(n)%jpon(i,j,k) = (biotic(n)%fpon(i,j,k-1) * Grid%tmask(i,j,k)   &
          - biotic(n)%fpon(i,j,k) * Grid%tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)
        biotic(n)%jpop(i,j,k) = (biotic(n)%fpop(i,j,k-1)  * Grid%tmask(i,j,k)  &
          - biotic(n)%fpop(i,j,k) * Grid%tmask(i,j,k+1)) / rho_dzt(i,j,k,taum1)

        ! Add production within box to flux
        biotic(n)%fpon(i,j,k) = biotic(n)%fpon(i,j,k) +                    &
          biotic(n)%jprod_pon(i,j,k) *                                     &
          rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
        biotic(n)%fpop(i,j,k) = biotic(n)%fpop(i,j,k) +                    &
          biotic(n)%jprod_pop(i,j,k) *                                     &
          rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
    enddo
  enddo
  enddo
elseif (biotic(n)%remin_density) then
elseif (biotic(n)%remin_lability) then
elseif (biotic(n)%remin_temp) then
elseif (biotic(n)%remin_viscosity) then
elseif (biotic(n)%remin_zoop_resp) then
endif


! Choose between associating particulate Fe with ballast and organic matter, 
! or just with organic matter
if (biotic(n)%fe_ballast_assoc) then
  do k=2,nk-1
    do j = jsc, jec
      do i = isc, iec
         ! Apply N change to Fe incorporating adsorption and desorption
        biotic(n)%jfe_ads(i,j,k)=min(biotic(n)%kfe_max_prime , (biotic(n)%kfe_org  &
          / 2.0 * (biotic(n)%fpon(i,j,k-1)+biotic(n)%fpon(i,j,k)) *            &
          biotic(n)%mass_2_n + biotic(n)%kfe_bal / 2.0 *                         &
          ((biotic(n)%fsio2(i,j,k-1) + biotic(n)%fsio2(i,j,k)) *               &
          60.0 + (biotic(n)%fcaco3(i,j,k-1) + biotic(n)%fcaco3(i,j,k)) * 100.0))&
          * biotic(n)%r_wsink ) * max(0.0,t_prog(ind_fed)%field(i,j,k,taum1))
        biotic(n)%jfe_des(i,j,k)=biotic(n)%kfe_des *                             &
          max(0.0,t_diag(ind_fep)%field(i,j,k))
        biotic(n)%jpofe(i,j,k) = (biotic(n)%jpon(i,j,k) *                &
          biotic(n)%mass_2_n + biotic(n)%jsio4(i,j,k) * 60.0 +           &
          biotic(n)%jcaco3(i,j,k) * 100.0) /                             &
          max(1.0e-30,biotic(n)%fpon(i,j,k-1) * biotic(n)%mass_2_n +     &
          biotic(n)%fsio2(i,j,k-1) * 60.0 + biotic(n)%fcaco3(i,j,k-1)    &
          * 100.0) * max(0.0,t_diag(ind_fep)%field(i,j,k))               &
          * biotic(n)%wsink *                                            &
          (1.0 - Grid%tmask(i,j,k) + Grid%tmask(i,j,k+1))
      enddo
    enddo
  enddo
else
  do k=2,nk-1
    do j=jsc,jec
      do i=isc,iec
         ! Apply N change to Fe incorporating adsorption and desorption
        biotic(n)%jfe_ads(i,j,k)=min(biotic(n)%kfe_max_prime , (biotic(n)%kfe_org  &
          / 2.0 * (biotic(n)%fpon(i,j,k-1)+biotic(n)%fpon(i,j,k)) *            &
          biotic(n)%mass_2_n + biotic(n)%kfe_bal / 2.0 *                         &
          ((biotic(n)%fsio2(i,j,k-1) + biotic(n)%fsio2(i,j,k)) *               &
          60.0 + (biotic(n)%fcaco3(i,j,k-1) + biotic(n)%fcaco3(i,j,k)) * 100.0))&
          * biotic(n)%r_wsink ) * max(0.0,t_prog(ind_fed)%field(i,j,k,taum1))
        biotic(n)%jfe_des(i,j,k)=biotic(n)%kfe_des *                             &
          max(0.0,t_diag(ind_fep)%field(i,j,k))
        biotic(n)%jpofe(i,j,k) = biotic(n)%jpon(i,j,k) /                 &
          max(1.0e-30,biotic(n)%fpon(i,j,k-1)) *                         &
          max(0.0,t_diag(ind_fep)%field(i,j,k)) * biotic(n)%wsink *       &
          (1.0 - Grid%tmask(i,j,k) + Grid%tmask(i,j,k+1))
        enddo
      enddo
    enddo
  endif


  do j=jsc,jec
    do i=isc,iec
      biotic(n)%fsio2(i,j,nk) = 0.0
      biotic(n)%fcaco3(i,j,nk) = 0.0
      biotic(n)%fpon(i,j,nk) = 0.0
      biotic(n)%fpop(i,j,nk) = 0.0
      biotic(n)%jfe_ads(i,j,nk) = min(biotic(n)%kfe_max_prime ,  &
        (biotic(n)%kfe_org / 2.0 * (biotic(n)%fpon(i,j,nk-1) +   &
        biotic(n)%fpon(i,j,nk)) * biotic(n)%mass_2_n +         &
        biotic(n)%kfe_bal / 2.0 * ((biotic(n)%fsio2(i,j,nk-1) +  &
        biotic(n)%fsio2(i,j,nk)) * 60.0 +                      &
        (biotic(n)%fcaco3(i,j,nk-1) + biotic(n)%fcaco3(i,j,nk))&
        * 100.0)) * biotic(n)%r_wsink ) *                      &
        max(0.0,t_prog(ind_fed)%field(i,j,nk,taum1))
      biotic(n)%jfe_des(i,j,nk)=biotic(n)%kfe_des *              &
        max(0.0,t_diag(ind_fep)%field(i,j,nk))

      ! Calculate regeneration term
      biotic(n)%jsio4(i,j,nk) = (biotic(n)%fsio2(i,j,nk-1) -   &
        biotic(n)%fsio2(i,j,nk)) / rho_dzt(i,j,nk,taum1)
      biotic(n)%jcaco3(i,j,nk) = (biotic(n)%fcaco3(i,j,nk-1) - &
        biotic(n)%fcaco3(i,j,nk)) / rho_dzt(i,j,nk,taum1)
      biotic(n)%jpon(i,j,nk) = (biotic(n)%fpon(i,j,nk-1) -     &
        biotic(n)%fpon(i,j,nk)) / rho_dzt(i,j,nk,taum1)
      biotic(n)%jpop(i,j,nk) = (biotic(n)%fpop(i,j,nk-1) -     &
        biotic(n)%fpop(i,j,nk)) / rho_dzt(i,j,nk,taum1)
      biotic(n)%jpofe(i,j,nk) = 0.0
    enddo
  enddo

    do j = jsc, jec
      do i = isc, iec
      kmax=min(km_c,grid_kmt(i,j))
      biotic(n)%flux_sio2(i,j) = biotic(n)%fsio2(i,j,kmax)
      biotic(n)%flux_caco3(i,j) = biotic(n)%fcaco3(i,j,kmax)
      biotic(n)%flux_pon(i,j) = biotic(n)%fpon(i,j,kmax)
      biotic(n)%flux_pop(i,j) = biotic(n)%fpop(i,j,kmax)
      enddo
    enddo

    do j = jsc, jec
      do i = isc, iec
      biotic(n)%jfe_sink(i,j,1)=                               &
        - max(0.0,t_diag(ind_fep)%field(i,j,1)) *        &
        biotic(n)%wsink / rho_dzt(i,j,1,taum1) * Grid%tmask(i,j,1)
      enddo
    enddo

  do k = 2, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jfe_sink(i,j,k)=                             &
          (max(0.0,t_diag(ind_fep)%field(i,j,k-1)) -     &
          max(0.0,t_diag(ind_fep)%field(i,j,k))) *       &
          biotic(n)%wsink / rho_dzt(i,j,k,taum1) * Grid%tmask(i,j,k)
      enddo
    enddo
  enddo

  ! CALCULATE SOURCE/SINK TERMS FOR EACH TRACER

  ! NO3
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
         biotic(n)%jdenit(i,j,k) = 0.0
        biotic(n)%jno3(i,j,k) =  biotic(n)%r_bio_tau_nitrif_s *             &
        max(0.0,t_prog(ind_nh4)%field(i,j,k,taum1)) -                       &
        biotic(n)%jprod_no3(i,j,k) - biotic(n)%jdenit(i,j,k)
      enddo
    enddo
  enddo

  do k=km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdenit(i,j,k) = 0.0
        biotic(n)%jno3(i,j,k) = biotic(n)%r_bio_tau_nitrif_d *                &
        max(0.0,t_prog(ind_nh4)%field(i,j,k,taum1)) -                       &
          biotic(n)%jdenit(i,j,k)
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_no3)%th_tendency(i,j,k) = t_prog(ind_no3)%th_tendency(i,j,k) +  &
          biotic(n)%jno3(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! NH4
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jnh4(i,j,k) = biotic(n)%jnh4_graz(i,j,k)               &
          - biotic(n)%jprod_nh4(i,j,k)                                   &
          + biotic(n)%r_bio_tau_don * max(0.0,t_prog(ind_don)%field(i,j,k,taum1)) &
          - biotic(n)%r_bio_tau_nitrif_s *                               &
          max(0.0,t_prog(ind_nh4)%field(i,j,k,taum1))                    &
          +(1.0 - biotic(n)%phi_don)*biotic(n)%jpon(i,j,k)  
      enddo
    enddo
  enddo

  do k=km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jnh4(i,j,k) =(1.0 - biotic(n)%phi_don) *               &
          biotic(n)%jpon(i,j,k)                                          &
          + biotic(n)%r_bio_tau_don * max(0.0,t_prog(ind_don)%field(i,j,k,taum1)) &
          - biotic(n)%r_bio_tau_nitrif_d *                               &
          max(0.0,t_prog(ind_nh4)%field(i,j,k,taum1)) 
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_nh4)%th_tendency(i,j,k) = t_prog(ind_nh4)%th_tendency(i,j,k) +  &
          biotic(n)%jnh4(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! PO4
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jpo4(i,j,k) = - biotic(n)%jprod_p_norm(i,j,k) -        &
          biotic(n)%jprod_p_fix(i,j,k) + biotic(n)%jpo4_graz(i,j,k) +    &
          biotic(n)%r_bio_tau_dop * max(0.0,t_prog(ind_dop)%field(i,j,k,taum1)) + &
         (1.0 - biotic(n)%phi_dop) * biotic(n)%jpop(i,j,k)
      enddo
    enddo
  enddo

  do k=km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jpo4(i,j,k) =  (1.0 - biotic(n)%phi_dop) *             &
          biotic(n)%jpop(i,j,k) + biotic(n)%r_bio_tau_dop *              &
          max(0.0,t_prog(ind_dop)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_po4)%th_tendency(i,j,k) = t_prog(ind_po4)%th_tendency(i,j,k) +  &
             biotic(n)%jpo4(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! Fe dissolved
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_fed)%th_tendency(i,j,k) = t_prog(ind_fed)%th_tendency(i,j,k) +  &
             (biotic(n)%jfe_graz(i,j,k) - biotic(n)%jprod_fed(i,j,k)    &
             - biotic(n)%jfe_ads(i,j,k) + biotic(n)%jfe_des(i,j,k) +    &
             biotic(n)%jpofe(i,j,k)) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  do k=km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_fed)%th_tendency(i,j,k) = t_prog(ind_fed)%th_tendency(i,j,k) + &
          (biotic(n)%jfe_des(i,j,k) - biotic(n)%jfe_ads(i,j,k) +        &
          biotic(n)%jpofe(i,j,k)) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! Fe particulate
  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_diag(ind_fep)%field(i,j,k) = t_diag(ind_fep)%field(i,j,k) + &
             (biotic(n)%jprod_pofe(i,j,k) + biotic(n)%jfe_ads(i,j,k) -  &
             biotic(n)%jfe_des(i,j,k) - biotic(n)%jpofe(i,j,k) +        &
             biotic(n)%jfe_sink(i,j,k)) * dtts
      enddo
    enddo
  enddo

  ! LDOC
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jldoc(i,j,k) = biotic(n)%jldoc(i,j,k)                    &
          - biotic(n)%r_bio_tau_ldoc * max(0.0,t_prog(ind_ldoc)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jldoc(i,j,k) =                                           &
          - biotic(n)%r_bio_tau_ldoc * max(0.0,t_prog(ind_ldoc)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_ldoc)%th_tendency(i,j,k) = t_prog(ind_ldoc)%th_tendency(i,j,k) + &
          biotic(n)%jldoc(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! DON
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdon(i,j,k) = biotic(n)%jdon(i,j,k)                      &
        +  biotic(n)%phi_don * biotic(n)%jpon(i,j,k)                       &
          - biotic(n)%r_bio_tau_don * max(0.0,t_prog(ind_don)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdon(i,j,k) = biotic(n)%phi_don * biotic(n)%jpon(i,j,k)  &
          - biotic(n)%r_bio_tau_don * max(0.0,t_prog(ind_don)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_don)%th_tendency(i,j,k) = t_prog(ind_don)%th_tendency(i,j,k) +      &
          biotic(n)%jdon(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! DOP
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdop(i,j,k) = biotic(n)%jdop(i,j,k)                      &
          +  biotic(n)%phi_dop * biotic(n)%jpop(i,j,k)                     &
          - biotic(n)%r_bio_tau_dop * max(0.0,t_prog(ind_dop)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jdop(i,j,k) = biotic(n)%phi_dop * biotic(n)%jpop(i,j,k)  &
          - biotic(n)%r_bio_tau_dop * max(0.0,t_prog(ind_dop)%field(i,j,k,taum1))
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_dop)%th_tendency(i,j,k) = t_prog(ind_dop)%th_tendency(i,j,k) +      &
          biotic(n)%jdop(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! SiO4
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        biotic(n)%jsio4(i,j,k) = biotic(n)%jsio4(i,j,k) -                  &
          biotic(n)%jprod_sio4(i,j,k)
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_sio4)%th_tendency(i,j,k) = t_prog(ind_sio4)%th_tendency(i,j,k) +  &
          biotic(n)%jsio4(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! ALK
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_alk)%th_tendency(i,j,k) = t_prog(ind_alk)%th_tendency(i,j,k) + &
          (2.0 * biotic(n)%jcaco3(i,j,k) - biotic(n)%jprod_alk(i,j,k)   &
          - biotic(n)%jno3(i,j,k) + biotic(n)%jnh4(i,j,k)) *            &
          rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_alk)%th_tendency(i,j,k) = t_prog(ind_alk)%th_tendency(i,j,k) + &
          (2.0 * biotic(n)%jcaco3(i,j,k) - biotic(n)%jno3(i,j,k) +      &
          biotic(n)%jnh4(i,j,k)) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! O2

  ! O2 production from nitrate, ammonia and nitrogen fixation and
  ! O2 consumption from production of NH4 from non-sinking particles,
  ! and DOM
  do k = 1, km_c
    do j =jsc, jec
      do i = isc, iec
        biotic(n)%jo2(i,j,k) = (biotic(n)%o_2_no3 *                        &
          biotic(n)%jprod_no3(i,j,k)                                       &
          + biotic(n)%o_2_nh4 * (biotic(n)%jprod_nh4(i,j,k) +              &
          biotic(n)%jprod_n_fix(i,j,k))) * Grid%tmask(i,j,k)

        ! If O2 is present
        if (t_prog(ind_o2)%field(i,j,k,taum1) .gt. biotic(n)%o2_min)       &
             then
          biotic(n)%jo2(i,j,k) = biotic(n)%jo2(i,j,k) - biotic(n)%o_2_nh4 *&
            (biotic(n)%jnh4_graz(i,j,k) + biotic(n)%r_bio_tau_don *        &
            max(0.0,t_prog(ind_don)%field(i,j,k,taum1)) + (1.0 - biotic(n)%phi_don) &
            * biotic(n)%jpon(i,j,k)) +                                     &
            biotic(n)%o_2_c * biotic(n)%jldoc(i,j,k) -                     &
            biotic(n)%o_2_nitrif * biotic(n)%r_bio_tau_nitrif_s *          &
            max(0.0,t_prog(ind_nh4)%field(i,j,k,taum1))
        endif
      enddo
    enddo
  enddo

  ! O2 consumption from production of NH4 from sinking particles and DOM
  ! O2 consumption from nitrification
  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
      ! If O2 is present
      if (t_prog(ind_o2)%field(i,j,k,taum1) .gt. biotic(n)%o2_min) then
        biotic(n)%jo2(i,j,k) = - biotic(n)%o_2_nh4 *                       &
          ((1.0 - biotic(n)%phi_don) * biotic(n)%jpon(i,j,k) +             &
          biotic(n)%r_bio_tau_don * max(0.0,t_prog(ind_don)%field(i,j,k,taum1))) +  &
            biotic(n)%o_2_c * biotic(n)%jldoc(i,j,k) -                     &
            biotic(n)%o_2_nitrif * biotic(n)%r_bio_tau_nitrif_d *          &
            max(0.0,t_prog(ind_nh4)%field(i,j,k,taum1))
        else
          biotic(n)%jo2(i,j,k) = 0.0
          endif
      enddo
    enddo
  enddo

  do k = 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_o2)%th_tendency(i,j,k) = t_prog(ind_o2)%th_tendency(i,j,k) + &
          biotic(n)%jo2(i,j,k) * rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  ! DIC
  do k = 1, km_c
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_dic)%th_tendency(i,j,k) = t_prog(ind_dic)%th_tendency(i,j,k) + &
          (biotic(n)%c_2_n * (biotic(n)%jno3(i,j,k) +                             &
          biotic(n)%jnh4(i,j,k)) - biotic(n)%jldoc(i,j,k) +                       &
          biotic(n)%jcaco3(i,j,k) - 0.5 * biotic(n)%jprod_alk(i,j,k)) *           &
          rho_dzt(i,j,k,taum1)
      enddo
    enddo
  enddo

  do k = km_c + 1, nk
    do j = jsc, jec
      do i = isc, iec
        t_prog(ind_dic)%th_tendency(i,j,k) = t_prog(ind_dic)%th_tendency(i,j,k) + &
          (biotic(n)%c_2_n * (biotic(n)%jno3(i,j,k)                               &
          + biotic(n)%jnh4(i,j,k)) - biotic(n)%jldoc(i,j,k) +                     &
          biotic(n)%jcaco3(i,j,k)) * rho_dzt(i,j,k,taum1)
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
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_flux_pon, biotic(n)%flux_pon(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_flux_pop, biotic(n)%flux_pop(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_flux_sio2, biotic(n)%flux_sio2(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_flux_caco3, biotic(n)%flux_caco3(:,:))
   call diagnose_2d_comp(Time, Grid, biotic(n)%id_htotal, biotic(n)%htotal(:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_alk, biotic(n)%jprod_alk(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_fed, biotic(n)%jprod_fed(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_n_fix, biotic(n)%jprod_n_fix(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_no3, biotic(n)%jprod_no3(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_nh4, biotic(n)%jprod_nh4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_p_fix, biotic(n)%jprod_p_fix(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_po4, biotic(n)%jprod_po4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_pofe, biotic(n)%jprod_pofe(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_pon, biotic(n)%jprod_pon(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_pop, biotic(n)%jprod_pop(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jprod_sio4,biotic(n)%jprod_sio4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jcaco3,biotic(n)%jcaco3(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jfe_ads, biotic(n)%jfe_ads(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jfe_des, biotic(n)%jfe_des(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jfe_graz, biotic(n)%jfe_graz(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jfe_sink, biotic(n)%jfe_sink(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jno3, biotic(n)%jno3(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jnh4, biotic(n)%jnh4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jnh4_graz, biotic(n)%jnh4_graz(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jpo4, biotic(n)%jpo4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jpo4_graz, biotic(n)%jpo4_graz(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jpofe, biotic(n)%jpofe(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jpon, biotic(n)%jpon(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jpop, biotic(n)%jpop(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jsio4, biotic(n)%jsio4(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jdenit, biotic(n)%jdenit(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jdon, biotic(n)%jdon(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jdop, biotic(n)%jdop(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jldoc, biotic(n)%jldoc(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_jo2, biotic(n)%jo2(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fpon, biotic(n)%fpon(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fpop, biotic(n)%fpop(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fracl, biotic(n)%fracl(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fsio2, biotic(n)%fsio2(:,:,:))
   call diagnose_3d_comp(Time, Grid, biotic(n)%id_fcaco3, biotic(n)%fcaco3(:,:,:))
enddo

return

end subroutine  ocean_bgc_restore_source
! </SUBROUTINE> NAME="ocean_bgc_restore_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_bgc_restore_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants for a given run
! and allocate diagnostic arrays
! </DESCRIPTION>
subroutine ocean_bgc_restore_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,  &
     T_prog, T_diag, taup1, model_time, grid_dat, grid_tmask, grid_kmt,         &
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
type(ocean_diag_tracer_type), dimension(:), intent(in)  :: T_diag
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
character(len=*), intent(in)                            :: grid_name
integer, dimension(3), intent(in)                       :: grid_tracer_axes
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

character(len=64), parameter    :: sub_name = 'ocean_bgc_restore_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

!----------------------------------------------------------------------
!
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
!----------------------------------------------------------------------

integer                                         :: done
integer                                         :: i, j, k, l, m, n
character(len=fm_field_name_len)                :: name
character(len=fm_field_name_len+1)              :: suffix
character(len=fm_field_name_len+3)              :: long_suffix
character(len=256)                              :: caller_str
integer                                         :: len_w
integer                                         :: len_e
integer                                         :: len_s
integer                                         :: len_n
real                                            :: total_alkalinity
real                                            :: total_ammonia
real                                            :: total_dic
real                                            :: total_don
real                                            :: total_dop
real                                            :: total_fediss
real                                            :: total_fepart
real                                            :: total_ldoc
real                                            :: total_o2
real                                            :: total_nitrate
real                                            :: total_phosphate
real                                            :: total_silicate
character(len=fm_string_len), allocatable       :: local_restart_file(:)
logical                                         :: fld_exist
integer                                         :: ind
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

dep_dry_file       =  fm_util_get_string ('dep_dry_file', scalar = .true.)
dep_dry_name       =  fm_util_get_string ('dep_dry_name', scalar = .true.)
dep_wet_file       =  fm_util_get_string ('dep_wet_file', scalar = .true.)
dep_wet_name       =  fm_util_get_string ('dep_wet_name', scalar = .true.)
no3_star_file      =  fm_util_get_string ('no3_star_file', scalar = .true.)
no3_star_name      =  fm_util_get_string ('no3_star_name', scalar = .true.)
po4_star_file      =  fm_util_get_string ('po4_star_file', scalar = .true.)
po4_star_name      =  fm_util_get_string ('po4_star_name', scalar = .true.)
sio4_star_file     =  fm_util_get_string ('sio4_star_file', scalar = .true.)
sio4_star_name     =  fm_util_get_string ('sio4_star_name', scalar = .true.)
alk_star_file      =  fm_util_get_string ('alk_star_file', scalar = .true.)
alk_star_name      =  fm_util_get_string ('alk_star_name', scalar = .true.)
fed_star_file      =  fm_util_get_string ('fed_star_file', scalar = .true.)
fed_star_name      =  fm_util_get_string ('fed_star_name', scalar = .true.)
htotal_scale_lo_in =  fm_util_get_real   ('htotal_scale_lo_in', scalar = .true.)
htotal_scale_hi_in =  fm_util_get_real   ('htotal_scale_hi_in', scalar = .true.)
htotal_in          =  fm_util_get_real   ('htotal_in', scalar = .true.)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str)

! Open up the Alk file for restoring
alk_star_id = init_external_field(alk_star_file,                &
                                  alk_star_name,                &
                                  domain = mpp_domain2d)
if (alk_star_id .eq. 0) then
  call mpp_error(FATAL,                                         &
       trim(sub_name) //                                        &
       ': Error: could not open alk_star file: ' //             &
       trim(alk_star_file))
endif

! Open up the Fed file for restoring
fed_star_id = init_external_field(fed_star_file,                &
                                  fed_star_name,                &
                                  domain = mpp_domain2d)
if (alk_star_id .eq. 0) then
  call mpp_error(FATAL,                                         &
       trim(sub_name) //                                        &
       ': Error: could not open fed_star file: ' //             &
       trim(fed_star_file))
endif

! Open up the NO3 file for restoring
no3_star_id = init_external_field(no3_star_file,                &
                                  no3_star_name,                &
                                  domain = mpp_domain2d)
if (no3_star_id .eq. 0) then
  call mpp_error(FATAL,                                         &
       trim(sub_name) //                                        &
       ': Error: could not open no3_star file: ' //             &
       trim(no3_star_file))
endif

! Open up the PO4 file for restoring
po4_star_id = init_external_field(po4_star_file,                &
                                  po4_star_name,                &
                                  domain = mpp_domain2d)
if (po4_star_id .eq. 0) then
  call mpp_error(FATAL,                                         &
       trim(sub_name) //                                        &
       ': Error: could not open po4_star file: ' //             &
       trim(po4_star_file))
endif

! Open up the SiO4 file for restoring
sio4_star_id = init_external_field(sio4_star_file,              &
                                  sio4_star_name,               &
                                  domain = mpp_domain2d)
if (sio4_star_id .eq. 0) then
  call mpp_error(FATAL,                                         &
       trim(sub_name) //                                        &
       ': Error: could not open sio4_star file: ' //            &
       trim(sio4_star_file))
endif

! Open up the files for boundary conditions
dep_wet_id = init_external_field(dep_wet_file,          &
                                     dep_wet_name,          &
                                     domain = mpp_domain2d)
if (dep_wet_id .eq. 0) then
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open wet deposition file: ' //          &
       trim(dep_wet_file))
endif

dep_dry_id = init_external_field(dep_dry_file,          &
                                     dep_dry_name,          &
                                     domain = mpp_domain2d)
if (dep_dry_id .eq. 0) then
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open dry deposition file: ' //          &
       trim(dep_dry_file))
endif

! set default values for htotal_scale bounds
htotal_scale_lo(:,:) = htotal_scale_lo_in
htotal_scale_hi(:,:) = htotal_scale_hi_in


! read in the namelists for each instance
do n = 1, instances

  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str)

  biotic(n)%compensation_depth     = fm_util_get_real   ('compensation_depth', scalar = .true.)
  biotic(n)%martin_coeff           = fm_util_get_real   ('martin_coeff', scalar = .true.)
  biotic(n)%ca_remin_depth         = fm_util_get_real   ('ca_remin_depth', scalar = .true.)
  biotic(n)%si_remin_depth         = fm_util_get_real   ('si_remin_depth', scalar = .true.)
  biotic(n)%soft_tissue_pump       = fm_util_get_logical('soft_tissue_pump', scalar = .true.)
  biotic(n)%stp_temperature        = fm_util_get_real   ('stp_temperature', scalar = .true.)
  biotic(n)%stp_salinity           = fm_util_get_real   ('stp_salinity', scalar = .true.)
  biotic(n)%stp_alkalinity         = fm_util_get_real   ('stp_alkalinity', scalar = .true.)
  biotic(n)%local_restart_file     = fm_util_get_string ('local_restart_file', scalar = .true.)
  biotic(n)%fe_ballast_assoc       = fm_util_get_logical('fe_ballast_assoc', scalar = .true.)
  biotic(n)%phi_dry                = fm_util_get_real   ('phi_dry', scalar = .true.)
  biotic(n)%phi_wet                = fm_util_get_real   ('phi_wet', scalar = .true.)
  biotic(n)%kfe_org                = fm_util_get_real   ('kfe_org', scalar = .true.)
  biotic(n)%kfe_bal                = fm_util_get_real   ('kfe_bal', scalar = .true.)
  biotic(n)%kfe_des                = fm_util_get_real   ('kfe_des', scalar = .true.)
  biotic(n)%kfe_max_prime          = fm_util_get_real   ('kfe_max_prime', scalar = .true.)
  biotic(n)%mass_2_n               = fm_util_get_real   ('mass_2_n', scalar = .true.)
  biotic(n)%n_2_p                  = fm_util_get_real   ('n_2_p', scalar = .true.)
  biotic(n)%n_2_p_fix              = fm_util_get_real   ('n_2_p_fix', scalar = .true.)
  biotic(n)%c_2_n                  = fm_util_get_real   ('c_2_n', scalar = .true.)
  biotic(n)%o_2_c                  = fm_util_get_real   ('o_2_c', scalar = .true.)
  biotic(n)%o_2_no3                = fm_util_get_real   ('o_2_no3', scalar = .true.)
  biotic(n)%o_2_nh4                = fm_util_get_real   ('o_2_nh4', scalar = .true.)
  biotic(n)%o_2_nitrif             = fm_util_get_real   ('o_2_nitrif', scalar = .true.)
  biotic(n)%o2_min                 = fm_util_get_real   ('o2_min', scalar = .true.)
  biotic(n)%bio_tau                = fm_util_get_real   ('bio_tau', scalar = .true.)
  biotic(n)%bio_tau_don            = fm_util_get_real   ('bio_tau_don', scalar = .true.)
  biotic(n)%bio_tau_dop            = fm_util_get_real   ('bio_tau_dop', scalar = .true.)
  biotic(n)%bio_tau_fix            = fm_util_get_real   ('bio_tau_fix', scalar = .true.)
  biotic(n)%bio_tau_ldoc           = fm_util_get_real   ('bio_tau_ldoc', scalar = .true.)
  biotic(n)%bio_tau_nh4            = fm_util_get_real   ('bio_tau_nh4', scalar = .true.)
  biotic(n)%bio_tau_nitrif_d       = fm_util_get_real   ('bio_tau_nitrif_d', scalar = .true.)
  biotic(n)%bio_tau_nitrif_s       = fm_util_get_real   ('bio_tau_nitrif_s', scalar = .true.)
  biotic(n)%kappa_eppley           = fm_util_get_real   ('kappa_eppley', scalar = .true.)
  biotic(n)%kappa_remin            = fm_util_get_real   ('kappa_remin', scalar = .true.)
  biotic(n)%Prodstar               = fm_util_get_real   ('Prodstar', scalar = .true.)
  biotic(n)%fdets0                 = fm_util_get_real   ('fdets0', scalar = .true.)
  biotic(n)%fdetl0                 = fm_util_get_real   ('fdetl0', scalar = .true.)
  biotic(n)%phi_don                = fm_util_get_real   ('phi_don', scalar = .true.)
  biotic(n)%phi_dop                = fm_util_get_real   ('phi_dop', scalar = .true.)
  biotic(n)%phi_ldoc               = fm_util_get_real   ('phi_ldoc', scalar = .true.)
  biotic(n)%gamma_det              = fm_util_get_real   ('gamma_det', scalar = .true.)
  biotic(n)%remin_density          = fm_util_get_logical('remin_density', scalar = .true.)
  biotic(n)%remin_lability         = fm_util_get_logical('remin_lability', scalar = .true.)
  biotic(n)%remin_ocmip2           = fm_util_get_logical('remin_ocmip2', scalar = .true.)
  biotic(n)%remin_protection       = fm_util_get_logical('remin_protection', scalar = .true.)
  biotic(n)%remin_simple           = fm_util_get_logical('remin_simple', scalar = .true.)
  biotic(n)%remin_temp             = fm_util_get_logical('remin_temp', scalar = .true.)
  biotic(n)%remin_viscosity        = fm_util_get_logical('remin_viscosity', scalar = .true.)
  biotic(n)%remin_zoop_resp        = fm_util_get_logical('remin_zoop_resp', scalar = .true.)
  biotic(n)%rpcaco3                = fm_util_get_real   ('rpcaco3', scalar = .true.)
  biotic(n)%rpsio2                 = fm_util_get_real   ('rpsio2', scalar = .true.)
  biotic(n)%wsink                  = fm_util_get_real   ('wsink', scalar = .true.)
  biotic(n)%sc_co2_0               = fm_util_get_real   ('sc_co2_0', scalar = .true.)
  biotic(n)%sc_co2_1               = fm_util_get_real   ('sc_co2_1', scalar = .true.)
  biotic(n)%sc_co2_2               = fm_util_get_real   ('sc_co2_2', scalar = .true.)
  biotic(n)%sc_co2_3               = fm_util_get_real   ('sc_co2_3', scalar = .true.)
  biotic(n)%sc_o2_0                = fm_util_get_real   ('sc_o2_0', scalar = .true.)
  biotic(n)%sc_o2_1                = fm_util_get_real   ('sc_o2_1', scalar = .true.)
  biotic(n)%sc_o2_2                = fm_util_get_real   ('sc_o2_2', scalar = .true.)
  biotic(n)%sc_o2_3                = fm_util_get_real   ('sc_o2_3', scalar = .true.)

  call fm_util_end_namelist(package_name, biotic(n)%name, caller = caller_str)

  biotic(n)%r_bio_tau          = 1.0 / biotic(n)%bio_tau
  biotic(n)%r_bio_tau_don      = 1.0 / biotic(n)%bio_tau_don
  biotic(n)%r_bio_tau_dop      = 1.0 / biotic(n)%bio_tau_dop
  biotic(n)%r_bio_tau_ldoc     = 1.0 / biotic(n)%bio_tau_ldoc
  biotic(n)%r_bio_tau_nh4      = 1.0 / biotic(n)%bio_tau_nh4
  biotic(n)%r_bio_tau_fix      = 1.0 / biotic(n)%bio_tau_fix
  biotic(n)%r_wsink            = 1.0 / biotic(n)%wsink
  biotic(n)%r_bio_tau_nitrif_s = 1.0 / biotic(n)%bio_tau_nitrif_s
  biotic(n)%r_bio_tau_nitrif_d = 1.0 / biotic(n)%bio_tau_nitrif_d

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

  call fm_util_start_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_prod', caller = caller_str)

  biotic(n)%r_bio_tau_prod%factor       =  fm_util_get_real          ('factor', scalar = .true.)
  biotic(n)%r_bio_tau_prod%coastal_only =  fm_util_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%r_bio_tau_prod%wlon         => fm_util_get_real_array    ('wlon')
  biotic(n)%r_bio_tau_prod%elon         => fm_util_get_real_array    ('elon')
  biotic(n)%r_bio_tau_prod%slat         => fm_util_get_real_array    ('slat')
  biotic(n)%r_bio_tau_prod%nlat         => fm_util_get_real_array    ('nlat')
  biotic(n)%r_bio_tau_prod%t_mask       => fm_util_get_logical_array ('t_mask')

  call fm_util_end_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_prod', caller = caller_str)

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

! read in the r_bio_tau_prod namelist data
do n = 1, instances

  if (associated(biotic(n)%r_bio_tau_prod%wlon)) then
    len_w = size(biotic(n)%r_bio_tau_prod%wlon)
  else
    len_w = 0
  endif
  if (associated(biotic(n)%r_bio_tau_prod%elon)) then
    len_e = size(biotic(n)%r_bio_tau_prod%elon)
  else
    len_e = 0
  endif
  if (associated(biotic(n)%r_bio_tau_prod%slat)) then
    len_s = size(biotic(n)%r_bio_tau_prod%slat)
  else
    len_s = 0
  endif
  if (associated(biotic(n)%r_bio_tau_prod%nlat)) then
    len_n = size(biotic(n)%r_bio_tau_prod%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif

  if (size(biotic(n)%r_bio_tau_prod%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif

  ! set all of the values to the default
  biotic(n)%r_bio_tau_prod%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process r_bio_tau_prod array for ', trim(biotic(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
    done = 0
    do l = 1, 12
      if (biotic(n)%r_bio_tau_prod%t_mask(l)) then
        if (done .eq. 0) then

          ! set the values via the input values, saving this time index
          ! afterwards
          write (stdoutunit,*) 'Assigning month ', l
          call set_array(biotic(n)%r_bio_tau_prod%mask(:,:,l),             &
                  isd, ied, jsd, jed, grid_xt, grid_yt, grid_kmt,       &
                  len_w, biotic(n)%r_bio_tau_prod%wlon,                    &
                  biotic(n)%r_bio_tau_prod%elon,                           &
                  biotic(n)%r_bio_tau_prod%slat,                           &
                  biotic(n)%r_bio_tau_prod%nlat,                           &
                  biotic(n)%r_bio_tau_prod%factor, 1.0,                    &
                  'Primary production limitation', biotic(n)%r_bio_tau_prod%coastal_only)
          done = l
        else

          ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          biotic(n)%r_bio_tau_prod%mask(:,:,l) = biotic(n)%r_bio_tau_prod%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

! multiply by the restoring factor
do n = 1, instances
  biotic(n)%r_bio_tau_prod%mask(:,:,:) =                                &
       biotic(n)%r_bio_tau * biotic(n)%r_bio_tau_prod%mask(:,:,:)
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

! read in the r_bio_tau_prod namelist data
do n = 1, instances

  if (associated(biotic(n)%r_bio_tau_prod%wlon)) then
    len_w = size(biotic(n)%r_bio_tau_prod%wlon)
  else
    len_w = 0
  endif
  if (associated(biotic(n)%r_bio_tau_prod%elon)) then
    len_e = size(biotic(n)%r_bio_tau_prod%elon)
  else
    len_e = 0
  endif
  if (associated(biotic(n)%r_bio_tau_prod%slat)) then
    len_s = size(biotic(n)%r_bio_tau_prod%slat)
  else
    len_s = 0
  endif
  if (associated(biotic(n)%r_bio_tau_prod%nlat)) then
    len_n = size(biotic(n)%r_bio_tau_prod%nlat)
  else
    len_n = 0
  endif

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif

  if (size(biotic(n)%r_bio_tau_prod%t_mask) .ne. 12) then
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif

  ! set all of the values to the default
  biotic(n)%r_bio_tau_prod%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then

    write (stdoutunit,*)
    write (stdoutunit,*) trim(note_header), 'Process r_bio_tau_prod array for ', trim(biotic(n)%name)
    write (stdoutunit,*)

    ! set values for this time-level
    done = 0
    do l = 1, 12
      if (biotic(n)%r_bio_tau_prod%t_mask(l)) then
        if (done .eq. 0) then
          ! set the values via the input values, saving this time index
          ! afterwards
          write (stdoutunit,*) 'Assigning month ', l
          call set_array(biotic(n)%r_bio_tau_prod%mask(:,:,l),             &
                  isd, ied, jsd, jed, grid_xt, grid_yt, grid_kmt,       &
                  len_w, biotic(n)%r_bio_tau_prod%wlon,                    &
                  biotic(n)%r_bio_tau_prod%elon,                           &
                  biotic(n)%r_bio_tau_prod%slat,                           &
                  biotic(n)%r_bio_tau_prod%nlat,                           &
                  biotic(n)%r_bio_tau_prod%factor, 1.0,                    &
                  'Primary production limitation', biotic(n)%r_bio_tau_prod%coastal_only)
          done = l
        else
          ! Duplicate the values for a previous time-level
          write (stdoutunit,*) 'Duplicating month ', done, ' as ', l
          biotic(n)%r_bio_tau_prod%mask(:,:,l) = biotic(n)%r_bio_tau_prod%mask(:,:,done)
        endif
      endif
    enddo
  endif

enddo

! multiply by the restoring factor
do n = 1, instances
  biotic(n)%r_bio_tau_prod%mask(:,:,:) =                                &
       biotic(n)%r_bio_tau * biotic(n)%r_bio_tau_prod%mask(:,:,:)
enddo

! initialize special arrays for remineralization
do n = 1, instances
  do k = 1, nk
    biotic(n)%r_intzscale_n(k) = biotic(n)%gamma_det /        &
        biotic(n)%wsink * grid_dzt(k)
    biotic(n)%r_1plusintzscale_n(k) = 1.0 / (1.0 +            &
        biotic(n)%r_intzscale_n(k))
    biotic(n)%r_1plusintzscale_si(k) = 1.0 / (1.0 +           &
        grid_dzt(k) / biotic(n)%si_remin_depth)
    biotic(n)%r_1plusintzscale_ca(k) = 1.0 / (1.0 +           &
        grid_dzt(k) / biotic(n)%ca_remin_depth)
    biotic(n)%zforg(k) = (grid_zw(k) /                        &
          biotic(n)%compensation_depth) ** (-biotic(n)%martin_coeff)
  enddo
enddo

! Read in additional information for a restart.

! We must process all of the instances before restoring any files
! as all fields must be registered before the fields are
! restored, and fields from different instances may be in the
! same file.

! Note that the restart file names here must be different from
! those for the tracer values.

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

! Set up analyses

! register the fields

suffix = '_' // package_name
long_suffix = ' (' // trim(package_name) // ')'

id_o2_sat = register_diag_field(trim(diag_name),                        &
     'o2_saturation' // trim(suffix), grid_tracer_axes(1:2),            &
     model_time, 'O2 saturation' // trim(long_suffix), ' ',             &
     missing_value = -1.0e+10)

do n = 1, instances

  if (biotic(n)%name(1:1) .eq. '_') then
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

  biotic(n)%id_sfc_flux_co2 = register_diag_field(trim(diag_name),              &
       'sfc_flux_co2' // trim(suffix), grid_tracer_axes(1:2),                   &
       model_time, 'Surface Flux - CO2' // trim(long_suffix), 'mol m^-2 s^-1',  &
       missing_value = -1.0e+10)

  biotic(n)%id_sfc_flux_o2 = register_diag_field(trim(diag_name),               &
       'sfc_flux_o2' // trim(suffix), grid_tracer_axes(1:2),                    &
       model_time, 'Surface Flux - O2' // trim(long_suffix), 'mol m^-2 s^-1',   &
       missing_value = -1.0e+10)

  biotic(n)%id_sfc_flux_fed = register_diag_field(trim(diag_name),              &
       'sfc_flux_fed' // trim(suffix), grid_tracer_axes(1:2),                   &
       model_time, 'Surface Flux - Fed' // trim(long_suffix), 'mol m^-2 s^-1',  &
       missing_value = -1.0e+10)

  biotic(n)%id_alpha = register_diag_field(trim(diag_name),                     &
       'alpha'// trim(suffix), grid_tracer_axes(1:2),                           &
       model_time, 'Alpha CO2' // trim(long_suffix), ' ',                       &
       missing_value = -1.0e+10)

  biotic(n)%id_csurf = register_diag_field(trim(diag_name),                     &
       'csurf'// trim(suffix), grid_tracer_axes(1:2),                           &
       model_time, 'CO2* water' // trim(long_suffix), ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_pco2surf = register_diag_field(trim(diag_name),                  &
       'pco2surf'// trim(suffix), grid_tracer_axes(1:2),                        &
       model_time, 'Oceanic pCO2' // trim(long_suffix), ' ',                    &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_pon = register_diag_field(trim(diag_name),                  &
       'flux_pon'// trim(suffix), grid_tracer_axes(1:2),                        &
       model_time, 'PON flux' // trim(long_suffix), ' ',                        &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_pop = register_diag_field(trim(diag_name),                  &
       'flux_pop'// trim(suffix), grid_tracer_axes(1:2),                        &
       model_time, 'POP flux' // trim(long_suffix), ' ',                        &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_sio2 = register_diag_field(trim(diag_name),                 &
       'flux_sio2'// trim(suffix), grid_tracer_axes(1:2),                       &
       model_time, 'Si flux' // trim(long_suffix), ' ',                         &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_caco3 = register_diag_field(trim(diag_name),                &
       'flux_caco3'// trim(suffix), grid_tracer_axes(1:2),                      &
       model_time, 'CaCO3 flux' // trim(long_suffix), ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_htotal = register_diag_field(trim(diag_name),                    &
       'htotal'// trim(suffix), grid_tracer_axes(1:2),                          &
       model_time, 'H+ ion concentration' // trim(long_suffix), ' ',            &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_alk = register_diag_field(trim(diag_name),                         &
       'jprod_alk'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'Restoring alkalinity-based production' // trim(long_suffix), ' ',   &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_fed = register_diag_field(trim(diag_name),                         &
       'jprod_fed'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'Restoring iron-based production' // trim(long_suffix), ' ',         &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_n_fix = register_diag_field(trim(diag_name),                       &
       'jprod_n_fix'// trim(suffix), grid_tracer_axes(1:3),                             &
       model_time, 'Nitrogen fixation' // trim(long_suffix), ' ',                       &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_no3 = register_diag_field(trim(diag_name),                         &
       'jprod_no3'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'Restoring NO3-based production' // trim(long_suffix), ' ',          &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_nh4 = register_diag_field(trim(diag_name),                         &
       'jprod_nh4'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'NH4-based production' // trim(long_suffix), ' ',                    &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_p_fix = register_diag_field(trim(diag_name),                       &
       'jprod_p_fix'// trim(suffix), grid_tracer_axes(1:3),                             &
       model_time, 'PO4 in nitrogen fixation' // trim(long_suffix), ' ',                &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_po4 = register_diag_field(trim(diag_name),                         &
       'jprod_po4'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'Restoring PO4-based production' // trim(long_suffix), ' ',          &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_pofe = register_diag_field(trim(diag_name),                        &
       'jprod_pofe'// trim(suffix), grid_tracer_axes(1:3),                              &
       model_time, 'Detrital iron production' // trim(long_suffix), ' ',                &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_pon = register_diag_field(trim(diag_name),                         &
       'jprod_pon'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'Detrital nitrogen production' // trim(long_suffix), ' ',            &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_pop = register_diag_field(trim(diag_name),                         &
       'jprod_pop'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'Detrital phosphorus production' // trim(long_suffix), ' ',          &
       missing_value = -1.0e+10)

  biotic(n)%id_jprod_sio4 = register_diag_field(trim(diag_name),                        &
       'jprod_sio4'// trim(suffix), grid_tracer_axes(1:3),                              &
       model_time, 'Restoring Si-based production' // trim(long_suffix), ' ',           &
       missing_value = -1.0e+10)

  biotic(n)%id_jcaco3 = register_diag_field(trim(diag_name),                            &
       'jcaco3'// trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'CaCO3 change' // trim(long_suffix), ' ',                            &
       missing_value = -1.0e+10)

  biotic(n)%id_jfe_ads = register_diag_field(trim(diag_name),                           &
       'jfe_ads'// trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'Iron adsorption' // trim(long_suffix), ' ',                         &
       missing_value = -1.0e+10)

  biotic(n)%id_jfe_des = register_diag_field(trim(diag_name),                           &
       'jfe_des'// trim(suffix), grid_tracer_axes(1:3),                                 &
       model_time, 'Iron desorption' // trim(long_suffix), ' ',                         &
       missing_value = -1.0e+10)

  biotic(n)%id_jfe_graz = register_diag_field(trim(diag_name),                          &
       'jfe_graz'// trim(suffix), grid_tracer_axes(1:3),                                &
       model_time, 'Dissolved iron source from grazing' // trim(long_suffix), ' ',      &
       missing_value = -1.0e+10)

  biotic(n)%id_jfe_sink = register_diag_field(trim(diag_name),                          &
       'jfe_sink'// trim(suffix), grid_tracer_axes(1:3),                                &
       model_time, 'Particulate iron sinking' // trim(long_suffix), ' ',                &
       missing_value = -1.0e+10)

  biotic(n)%id_jno3 = register_diag_field(trim(diag_name),                              &
       'jno3'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'NO3 source' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_jnh4 = register_diag_field(trim(diag_name),                              &
       'jnh4'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'NH4 source' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_jnh4_graz = register_diag_field(trim(diag_name),                         &
       'jnh4_graz'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'NH4 source from grazing' // trim(long_suffix), ' ',                 &
       missing_value = -1.0e+10)

  biotic(n)%id_jpo4_graz = register_diag_field(trim(diag_name),                         &
       'jpo4_graz'// trim(suffix), grid_tracer_axes(1:3),                               &
       model_time, 'PO4 source from grazing' // trim(long_suffix), ' ',                 &
       missing_value = -1.0e+10)

  biotic(n)%id_jpo4 = register_diag_field(trim(diag_name),                              &
       'jpo4'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'PO4 source' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_jpofe = register_diag_field(trim(diag_name),                             &
       'jpofe'// trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Loss of sinking iron' // trim(long_suffix), ' ',                    &
       missing_value = -1.0e+10)

  biotic(n)%id_jpon = register_diag_field(trim(diag_name),                              &
       'jpon'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Loss of sinking nitrogen' // trim(long_suffix), ' ',                &
       missing_value = -1.0e+10)

  biotic(n)%id_jpop = register_diag_field(trim(diag_name),                              &
       'jpop'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'Loss of sinking phosphorus' // trim(long_suffix), ' ',              &
       missing_value = -1.0e+10)

  biotic(n)%id_jsio4 = register_diag_field(trim(diag_name),                             &
       'jsio4'// trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'SiO4 source' // trim(long_suffix), ' ',                             &
       missing_value = -1.0e+10)

  biotic(n)%id_jdenit = register_diag_field(trim(diag_name),                            &
       'jdenit'// trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'Denitrification' // trim(long_suffix), ' ',                         &
       missing_value = -1.0e+10)

  biotic(n)%id_jdon = register_diag_field(trim(diag_name),                              &
       'jdon'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'DON source' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_jdop = register_diag_field(trim(diag_name),                              &
       'jdop'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'DOP source' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_jldoc = register_diag_field(trim(diag_name),                             &
       'jldoc'// trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Labile DOC source' // trim(long_suffix), ' ',                       &
       missing_value = -1.0e+10)

  biotic(n)%id_jo2 = register_diag_field(trim(diag_name),                               &
       'jo2'// trim(suffix), grid_tracer_axes(1:3),                                     &
       model_time, 'O2 source' // trim(long_suffix), ' ',                               &
       missing_value = -1.0e+10)

  biotic(n)%id_fpon = register_diag_field(trim(diag_name),                              &
       'fpon'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'PON change' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_fpop = register_diag_field(trim(diag_name),                              &
       'fpop'// trim(suffix), grid_tracer_axes(1:3),                                    &
       model_time, 'POP change' // trim(long_suffix), ' ',                              &
       missing_value = -1.0e+10)

  biotic(n)%id_fracl = register_diag_field(trim(diag_name),                             &
       'fracl'// trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Fraction large phytoplankton' // trim(long_suffix), ' ',            &
       missing_value = -1.0e+10)

  biotic(n)%id_fsio2 = register_diag_field(trim(diag_name),                             &
       'fsio2'// trim(suffix), grid_tracer_axes(1:3),                                   &
       model_time, 'Si change' // trim(long_suffix), ' ',                               &
       missing_value = -1.0e+10)

  biotic(n)%id_fcaco3 = register_diag_field(trim(diag_name),                            &
       'fcaco3'// trim(suffix), grid_tracer_axes(1:3),                                  &
       model_time, 'CaCO3 change' // trim(long_suffix), ' ',                            &
       missing_value = -1.0e+10)

enddo

! integrate the total concentrations of some tracers
! for the start of the run

! Use taup1 time index for the start of a run, and taup1 time
! index for the end of a run so that we are integrating the
! same time level and should therefore get identical results

do n = 1, instances
  total_alkalinity = 0.0
  total_ammonia = 0.0
  total_dic = 0.0
  total_don = 0.0
  total_dop = 0.0
  total_fediss = 0.0
  total_fepart = 0.0
  total_ldoc = 0.0
  total_nitrate = 0.0
  total_o2 = 0.0
  total_phosphate = 0.0
  total_silicate = 0.0
  do k = 1,nk
    do j = jsc, jec
      do i = isc, iec
        total_nitrate = total_nitrate +                         &
             t_prog(biotic(n)%ind_no3)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_ammonia = total_ammonia +                         &
             t_prog(biotic(n)%ind_nh4)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_phosphate = total_phosphate +                     &
             t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_fediss = total_fediss +                           &
             t_prog(biotic(n)%ind_fed)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_fepart = total_fepart +                           &
             t_diag(biotic(n)%ind_fep)%field(i,j,k) *           &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_silicate = total_silicate +                       &
             t_prog(biotic(n)%ind_sio4)%field(i,j,k,taup1) *    &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_don = total_don +                                 &
             t_prog(biotic(n)%ind_don)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_dop = total_dop +                                 &
             t_prog(biotic(n)%ind_dop)%field(i,j,k,taup1) *     &
             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
        total_ldoc = total_ldoc +                               &
             t_prog(biotic(n)%ind_ldoc)%field(i,j,k,taup1) *    &
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

  call mpp_sum(total_nitrate)
  call mpp_sum(total_ammonia)
  call mpp_sum(total_phosphate)
  call mpp_sum(total_fediss)
  call mpp_sum(total_fepart)
  call mpp_sum(total_silicate)
  call mpp_sum(total_don)
  call mpp_sum(total_dop)
  call mpp_sum(total_ldoc)
  call mpp_sum(total_o2)
  call mpp_sum(total_dic)
  call mpp_sum(total_alkalinity)

  write (stdoutunit,*) '  Instance ', trim(biotic(n)%name)
  write (stdoutunit,                                              &
         '(/'' Total nitrate  = '',es19.12,'' Gmol-N'')')       &
              total_nitrate * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total ammonia  = '',es19.12,'' Gmol-N'')')       &
              total_ammonia * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total phosphate  = '',es19.12,'' Gmol-P'')')     &
              total_phosphate * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total fediss  = '',es19.12,'' Gmol-Fe'')')       &
              total_fediss * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total fepart  = '',es19.12,'' Gmol-Fe'')')       &
              total_fepart * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total silicate  = '',es19.12,'' Gmol-Si'')')     &
              total_silicate * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total DON  = '',es19.12,'' Gmol-C'')')           &
              total_DON * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total DOP  = '',es19.12,'' Gmol-P'')')           &
              total_DOP * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total LDOC  = '',es19.12,'' Gmol-C'')')          &
              total_LDOC * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total O2  = '',es19.12,'' Gmol-O'')')            &
              total_o2 * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total DIC  = '',es19.12,'' Gmol-C'')')           &
              total_dic * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total alkalinity  = '',es19.12,'' Geq'')')       &
              total_alkalinity * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total nitrogen  = '',es19.12,'' Gmol-N'')')      &
              (total_nitrate + total_don) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total phosphorus  = '',es19.12,'' Gmol-P'')')    &
              (total_phosphate + total_dop) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total real O2  = '',es19.12,'' Gmol-O'')')       &
              (total_o2 + biotic(n)%o_2_no3 * total_nitrate) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total Carbon  = '',es19.12,'' Gmol-C'')')        &
              (total_dic + biotic(n)%c_2_n * total_don + total_ldoc) * 1.0e-09
  write (stdoutunit,                                              &
         '(/'' Total real alkalinity  = '',es19.12,'' Geq'')')  &
              (total_alkalinity + total_nitrate) * 1.0e-09
enddo

write(stdoutunit,*)
write(stdoutunit,*) trim(note_header), 'Tracer runs initialized'
write(stdoutunit,*)

return

end subroutine  ocean_bgc_restore_start
! </SUBROUTINE> NAME="ocean_bgc_restore_start"

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
real, dimension(isd:,jsd:), intent(out)         :: array
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

!
!       Return true if w <= x_in <= e, taking into account the
!       periodicity of longitude.
!
!       x_in    = value to test
!
!       w       = west longitude of boundary
!
!       e       = east longitude of boundary
!

function lon_between(x_in, w, e)

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

end module  ocean_bgc_restore_mod
