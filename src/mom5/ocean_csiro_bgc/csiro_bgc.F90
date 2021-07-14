!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! 
!<CONTACT EMAIL="Richard.Matear@csiro.au"> Richard Matear
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Chamberlain@csiro.au"> Matt Chamberlain
!</CONTACT>
!
!<REVIEWER EMAIL=""> 
!</REVIEWER>
!
!<OVERVIEW>
! CSIRO v3 bgc model
!</OVERVIEW>
!
!<DESCRIPTION>
!       Generic setup of the CSIRO BGC model; NO3-based NPZD cycle, 
!         with coupled carbon cycle
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! </REFERENCE>
!
! <REFERENCE>
! </REFERENCE>
!
! </INFO>
!
!------------------------------------------------------------------
!
!       Module csiro_bgc_mod
!
!	csiro bgc model version 3       
!
!------------------------------------------------------------------

module  csiro_bgc_mod  !{

!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------

!----------------------------------------------------------------------
!
!       Modules
!
!----------------------------------------------------------------------

use diag_manager_mod,         only: send_data
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use field_manager_mod,        only: fm_get_index
use fms_mod,                  only: write_data, write_version_number
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use mpp_mod,                  only: mpp_pe, mpp_root_pe, mpp_sync
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use fms_io_mod,               only: register_restart_field, save_restart, restore_state
use fms_io_mod,               only: restart_file_type, reset_field_pointer

use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type
use mpp_domains_mod,       only: mpp_global_field, mpp_global_sum, BITWISE_EXACT_SUM

use time_manager_mod,         only: get_date
use time_interp_external_mod, only: time_interp_external

use ocean_parameters_mod,     only: rho0

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use fm_util_mod, only: fm_util_check_for_bad_fields, fm_util_set_value
use fm_util_mod, only: fm_util_get_string, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod, only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod, only: fm_util_start_namelist, fm_util_end_namelist
! use fm_util_mod, only: domain, grid, time, dtts
! use fm_util_mod, only: isc, iec, jsc, jec, nk, isd, ied, jsd, jed 
! use fm_util_mod, only: taum1, tau, taup1 
! use fm_util_mod, only: t_prog, t_diag
! use fm_util_mod, only: indsal, indtemp
! use fm_util_mod, only: end_of_year, end_of_month

use ocean_types_mod,    only: ocean_thickness_type
use ocean_types_mod,    only: ocean_density_type
use ocean_types_mod,    only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,    only: ocean_time_type

use ocean_types_mod,    only: ocean_public_type

use ocean_tracer_diag_mod,   only: calc_mixed_layer_depth

!----------------------------------------------------------------------
!       force all variables to be "typed"
!----------------------------------------------------------------------

implicit none

!----------------------------------------------------------------------
!
!       Make all routines and variables private by default
!
!----------------------------------------------------------------------
private

!----------------------------------------------------------------------
!
!       Public routines
!
!----------------------------------------------------------------------
! order of call  init, start, loop(sbc, source, bbc, tracer), end
public  :: csiro_bgc_bbc
public  :: csiro_bgc_end
public  :: csiro_bgc_init
public  :: csiro_bgc_sbc
public  :: csiro_bgc_source
public  :: csiro_bgc_start
public  :: csiro_bgc_tracer
public  :: csiro_bgc_virtual_fluxes

!----------------------------------------------------------------------
!
!       Private routines
!
!----------------------------------------------------------------------

private :: allocate_arrays

!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------

character(len=fm_field_name_len), parameter     :: package_name = 'csiro_bgc'
character(len=48), parameter                    :: mod_name = 'csiro_bgc_mod'
character(len=fm_string_len), parameter         :: default_file_in = 'INPUT/csiro_bgc.res.nc'
character(len=fm_string_len), parameter         :: default_file_out = 'RESTART/csiro_bgc.res.nc'

integer, parameter                              :: ntr_bmax = 10

!-------------------------------------------------------
! private types
!--------------------------------------------------------

type biotic_type  !{
! set-up tracer indexes
  integer                               :: ntr_bgc = ntr_bmax
  integer,  dimension(ntr_bmax)         :: ind_bgc
  integer,  dimension(ntr_bmax)         :: id_bgc_stf
  integer,  dimension(ntr_bmax)         :: id_bgc_btf
  integer,  dimension(ntr_bmax)         :: id_bgc_vstf
  integer,  dimension(ntr_bmax)         :: id_bgc_src

  logical                               :: init
  character(len=fm_field_name_len)      :: name
! arrays for air-sea fluxes
  real                                      :: sal_global
  real, allocatable, dimension(:)           :: bgc_global
  real, allocatable, dimension(:,:)         :: htotal
  real, allocatable, dimension(:,:)         :: alpha
  real, allocatable, dimension(:,:)         :: csat
  real, allocatable, dimension(:,:)         :: csat_csurf
  real, allocatable, dimension(:,:)         :: csat_acsurf
  real, allocatable, dimension(:,:)         :: csurf
  real, allocatable, dimension(:,:)         :: acsurf
  real, allocatable, dimension(:,:)         :: dpco2
  real, allocatable, dimension(:,:)         :: pco2surf
  real, allocatable, dimension(:,:)         :: paco2surf
  real, allocatable, dimension(:,:)         :: pco2atm
  real, allocatable, dimension(:,:)         :: paco2atm
  real, allocatable, dimension(:,:)         :: det_sediment       ! mmol(NO3) m-2 in DET sitting at base of column as sediment. 
  real, allocatable, dimension(:,:)         :: caco3_sediment     ! mmol(CaCO3) m-2 sitting at base of column as sediment. 
  real, allocatable, dimension(:,:)         :: det_sed_remin      ! mmol(NO3) m-2 s-1, rate of remineralisation of DET in sediment.  
  real, allocatable, dimension(:,:)         :: caco3_sed_remin    ! mmol m-2 s-1, rate of remineralisation of CaCO3 in sediment.  
  real, allocatable, dimension(:,:)         :: det_sed_depst      ! mmol(NO3) m-2 s-1, rate of deposition of DET in sediment.  
  real, allocatable, dimension(:,:)         :: caco3_sed_depst    ! mmol m-2 s-1, rate of deposition of CaCO3 in sediment.  
  real, allocatable, dimension(:,:)         :: sio2
  real, allocatable, dimension(:,:)         :: po4
  real, allocatable, dimension(:,:,:)       :: vstf
end type biotic_type  !}

!----------------------------------------------------------------------
!       Public variables
!----------------------------------------------------------------------

logical, public :: do_csiro_bgc

!----------------------------------------------------------------------
!       Private variables
!----------------------------------------------------------------------
integer                                 :: package_index

! set the tracer index for the various tracers
integer :: id_po4, id_dic, id_alk, id_o2, id_no3, id_phy, id_det, id_zoo &
      , id_caco3, id_adic, id_fe, id_caco3_sediment, id_det_sediment
! internal pointer to make reading the code easier
integer,public :: ind_po4 = -1
integer,public :: ind_dic = -1 
integer,public :: ind_alk = -1
integer,public :: ind_o2 = -1
integer,public :: ind_no3 = -1
integer,public :: ind_phy = -1
integer,public :: ind_det = -1
integer,public :: ind_zoo = -1
integer,public :: ind_caco3 = -1
integer,public :: ind_adic = -1
integer,public :: ind_fe = -1
character*6  :: qbio_model
integer      :: bio_version    ! version of the bgc module to use
logical      :: zero_floor     ! apply hard floor to bgc tracers 
logical      :: sw_thru_ice    ! make this true in a coupled model, so bgc knows swflx is already modified for the presense of ice.  
logical      :: gasx_from_file ! use gasx exchange coefficients from provided files. mac, may13.  
logical      :: ice_file4gasx  ! make this true in a ocean-only model, option to control whether to use ice file or the model ice when determining masking effect of ice on co2 gas exchange. 
logical      :: use_access_co2=.false. ! determine whether to use ACCESS atmospheric CO2 or a CO2 file to drive the ocean's anthropogenic CO2 tracer. mac, may13.  

integer  :: id_clock_csiro_obgc

integer                                 :: atmpress_id
character*128                           :: atmpress_file    
character*32                            :: atmpress_name    
real, allocatable, dimension(:,:)       :: fice_t
integer                                 :: id_light_limit = -1
integer                                 :: id_adic_intmld = -1
integer                                 :: id_dic_intmld = -1
integer                                 :: id_o2_intmld = -1
integer                                 :: id_no3_intmld = -1
integer                                 :: id_fe_intmld = -1
integer                                 :: id_phy_intmld = -1
integer                                 :: id_det_intmld = -1
integer                                 :: id_pprod_gross_intmld = -1
integer                                 :: id_npp_intmld = -1
integer                                 :: id_radbio_intmld = -1
integer                                 :: id_wdet100avg = -1
integer                                 :: id_radbio1 = -1
integer                                 :: id_radbio3d = -1
integer                                 :: id_wdet100 = -1
integer                                 :: id_npp1 = -1
integer                                 :: id_npp2d = -1
integer                                 :: id_npp3d = -1
integer                                 :: id_pprod_gross = -1
integer                                 :: id_pprod_gross_2d = -1
integer                                 :: id_zprod_gross = -1
integer                                 :: id_kw_o2  = -1
integer                                 :: id_o2_sat = -1
integer                                 :: id_sc_o2  = -1
integer                                 :: id_pco2 = -1, id_paco2 = -1
integer                                 :: id_co2_sat = -1, id_aco2_sat = -1
integer                                 :: id_caco3_sed_remin, id_det_sed_remin
integer                                 :: id_caco3_sed_depst, id_det_sed_depst
integer                                 :: id_total_aco2_flux, id_total_co2_flux
real, allocatable, dimension(:,:)       :: kw_co2 
real, allocatable, dimension(:,:)       :: kw_o2
real, allocatable, dimension(:,:)       :: patm_t
integer                                 :: pistonveloc_id
integer                                 :: aco2_id
integer                                 :: seaicefract_id
character*128                           :: pistonveloc_file
character*32                            :: pistonveloc_name
character*128                           :: aco2_file
character*32                            :: aco2_name
real, allocatable, dimension(:,:)       :: sc_o2
real, allocatable, dimension(:,:)       :: sc_co2
real, allocatable, dimension(:,:)       :: o2_saturation
character*128                           :: seaicefract_file
character*32                            :: seaicefract_name
character*128                           :: dust_file
character*32                            :: dust_name
integer                                 :: dust_id
real, allocatable, dimension(:,:)       :: dust_t

real, allocatable, dimension(:,:)       :: xkw_t
real, allocatable, dimension(:,:)       :: aco2
real, allocatable, dimension(:,:)             :: htotalhi
real, allocatable, dimension(:,:)             :: htotallo

type(biotic_type), allocatable, dimension(:)    :: biotic
integer                                         :: instances
real, allocatable, dimension(:)                 :: tk
real, allocatable, dimension(:)                 :: ts
real, allocatable, dimension(:)                 :: ts2
real, allocatable, dimension(:)                 :: ts3
real, allocatable, dimension(:)                 :: ts4
real, allocatable, dimension(:)                 :: ts5
real, allocatable, dimension(:)                 :: tt


! rjm biotic parameters 
integer :: n_eudepth=4
real :: p_k = 0.1
real :: rain_ratio = 8.48
real :: s_npp = -1
real :: ratio_poc(ntr_bmax)
real :: ratio_pic(ntr_bmax)

real, allocatable, dimension(:,:,:) :: poc
real, allocatable, dimension(:,:,:) :: pic
real, allocatable, dimension(:,:) :: poc_tot
real, allocatable, dimension(:,:) :: pic_tot
real, allocatable, dimension(:,:) :: pmax_growth
real, allocatable, dimension(:,:) :: pp
real, allocatable, dimension(:) :: fmin_poc
real, allocatable, dimension(:) :: fmin_pic
real, allocatable, dimension(:,:,:) :: biotr
real, allocatable, dimension(:,:) :: light_limit
real, allocatable, dimension(:,:) :: adic_intmld,dic_intmld,o2_intmld,no3_intmld,fe_intmld,phy_intmld,det_intmld
real, allocatable, dimension(:,:) :: pprod_gross_intmld,npp_intmld,radbio_intmld,wdet100avg
real, allocatable, dimension(:,:,:) :: radbio3d
real, allocatable, dimension(:,:) :: wdet100
real, allocatable, dimension(:,:) :: npp2d
real, allocatable, dimension(:,:,:) :: npp3d
real, allocatable, dimension(:,:,:) :: pprod_gross
real, allocatable, dimension(:,:) :: pprod_gross_2d
real, allocatable, dimension(:,:,:) :: zprod_gross
real, allocatable, dimension(:) :: ray
real, allocatable, dimension(:) :: dummy
real :: dummy1
real, allocatable, dimension(:) :: tracer_sources
real, allocatable, dimension(:,:) :: area_k
real, allocatable, dimension(:,:) :: tmp

integer :: global_sum_flag      ! flag for mpp_global_sum


integer                                 :: alphabio_id
real, allocatable, dimension(:,:)       :: alphabio
integer                                 :: parbio_id
real, allocatable, dimension(:,:)       :: parbio
integer                                 :: kwbio_id
real, allocatable, dimension(:,:)       :: kwbio
integer                                 :: kcbio_id
real, allocatable, dimension(:,:)       :: kcbio
integer                                 :: abio_id
real, allocatable, dimension(:,:)       :: abio
integer                                 :: bbio_id
real, allocatable, dimension(:,:)       :: bbio
integer                                 :: cbio_id
real, allocatable, dimension(:,:)       :: cbio
integer                                 :: k1bio_id
real, allocatable, dimension(:,:)       :: k1bio
integer                                 :: muepbio_id
real, allocatable, dimension(:,:)       :: muepbio
integer                                 :: muepsbio_id
real, allocatable, dimension(:,:)       :: muepsbio
integer                                 :: gam1bio_id
real, allocatable, dimension(:,:)       :: gam1bio
integer                                 :: gbio_id
real, allocatable, dimension(:,:)       :: gbio
integer                                 :: epsbio_id
real, allocatable, dimension(:,:)       :: epsbio
integer                                 :: muezbio_id
real, allocatable, dimension(:,:)       :: muezbio
integer                                 :: gam2bio_id
real, allocatable, dimension(:,:)       :: gam2bio
integer                                 :: muedbio_id
real, allocatable, dimension(:,:)       :: muedbio
integer                                 :: muecaco3_id
real, allocatable, dimension(:,:)       :: muecaco3
integer                                 :: muedbio_sed_id
real, allocatable, dimension(:,:)       :: muedbio_sed
integer                                 :: muecaco3_sed_id
real, allocatable, dimension(:,:)       :: muecaco3_sed
integer                                 :: wdetbio_id
real, allocatable, dimension(:,:)       :: wdetbio
integer                                 :: wcaco3_id
real, allocatable, dimension(:,:)       :: wcaco3
integer                                 :: nat_co2_id
real, allocatable, dimension(:,:)       :: nat_co2
integer                                 :: tscav_fe_id
real, allocatable, dimension(:,:)       :: tscav_fe
integer                                 :: fe_bkgnd_id
real, allocatable, dimension(:,:)       :: fe_bkgnd
integer                                 :: f_inorg_id
real, allocatable, dimension(:,:)       :: f_inorg


! for extra restart file(s)
integer                          :: id_restart(2)=0
type(restart_file_type), save    :: sed_restart


! Tracer names

character(5), dimension(10) :: tracer_name

!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------

contains

!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays
! </DESCRIPTION>
!

subroutine allocate_arrays (isc, iec, jsc, jec, isd, ied, jsd, jed, nk)  !{

!       arguments

!       local variables

integer, intent(in)                                        :: isd, isc
integer, intent(in)                                        :: ied, iec
integer, intent(in)                                        :: jsd, jsc
integer, intent(in)                                        :: jed, jec
integer, intent(in)                                        :: nk
integer :: m
integer :: k
integer :: j
integer :: i
integer :: n
integer :: ntr_bgc

!-----------------------------------------------------------------------
!     start executable code
!-----------------------------------------------------------------------

allocate( xkw_t(isd:ied,jsd:jed) )
allocate( patm_t(isd:ied,jsd:jed) )
allocate( fice_t(isd:ied,jsd:jed) )
allocate( aco2(isd:ied,jsd:jed) )
allocate( dust_t(isd:ied,jsd:jed) )

allocate( sc_o2(isc:iec,jsc:jec) )
allocate( sc_co2(isc:iec,jsc:jec) )
allocate( kw_o2(isc:iec,jsc:jec) )
allocate( kw_co2(isc:iec,jsc:jec) )
allocate( htotalhi(isd:ied,jsd:jed) )
allocate( htotallo(isd:ied,jsd:jed) )


allocate( o2_saturation(isc:iec,jsc:jec) )

allocate( tt(isc:iec) )
allocate( tk(isc:iec) )
allocate( ts(isc:iec) )
allocate( ts2(isc:iec) )
allocate( ts3(isc:iec) )
allocate( ts4(isc:iec) )
allocate( ts5(isc:iec) )

allocate( poc(isc:iec,jsc:jec,nk) )
allocate( pic(isc:iec,jsc:jec,nk) )
allocate( poc_tot(isc:iec,jsc:jec) )
allocate( pic_tot(isc:iec,jsc:jec) )
allocate( pmax_growth(isc:iec,nk) )
allocate( pp(isc:iec,nk) )
allocate( fmin_poc(nk) )
allocate( fmin_pic(nk) )
  ntr_bgc = biotic(1)%ntr_bgc
allocate( ray(nk) )
allocate( biotr(isc:iec,nk,ntr_bgc) )
allocate( light_limit(isc:iec,jsc:jec) )
allocate( adic_intmld(isc:iec,jsc:jec) )
allocate( dic_intmld(isc:iec,jsc:jec) )
allocate( o2_intmld(isc:iec,jsc:jec) )
allocate( no3_intmld(isc:iec,jsc:jec) )
allocate( phy_intmld(isc:iec,jsc:jec) )
allocate( det_intmld(isc:iec,jsc:jec) )
allocate( pprod_gross_intmld(isc:iec,jsc:jec) )
allocate( npp_intmld(isc:iec,jsc:jec) )
allocate( radbio_intmld(isc:iec,jsc:jec) )
allocate( wdet100avg(isc:iec,jsc:jec) )
allocate( radbio3d(isc:iec,jsc:jec,nk) )
allocate( wdet100(isc:iec,jsc:jec) )
allocate( npp2d(isc:iec,jsc:jec) )
allocate( npp3d(isc:iec,jsc:jec,nk) )
allocate( pprod_gross(isc:iec,jsc:jec,nk) )
allocate( pprod_gross_2d(isc:iec,jsc:jec) )
allocate( zprod_gross(isc:iec,jsc:jec,nk) )

allocate (tmp(isd:ied,jsd:jed) )
allocate ( tracer_sources(0:nk) )
allocate(area_k(isd:ied,jsd:jed) )    
allocate( dummy(isc:iec) )

!       initialize some arrays

xkw_t(:,:)         = 0.0
patm_t(:,:)        = 0.0
fice_t(:,:)        = 0.0
dust_t(:,:)        = 0.0
aco2(:,:)          = 0.0 
sc_co2(:,:)        = 0.0
kw_co2(:,:)        = 0.0
sc_o2(:,:)         = 0.0
kw_o2(:,:)         = 0.0
o2_saturation(:,:) = 0.0

tt(:)              = 0.0
tk(:)              = 0.0
ts(:)              = 0.0
ts2(:)             = 0.0
ts3(:)             = 0.0
ts4(:)             = 0.0
ts5(:)             = 0.0

!       allocate biotic array elements

do n = 1, instances  !{
! rjm problem -  the allocation of the arrays is done
! after the initialization routine which cause problems
! solution - move the call to array allocation to the init subroutine
!  ntr_bgc = biotic(n)%ntr_bgc
!  allocate( biotic(n)%ind_bgc(1:ntr_bgc) )

  allocate( biotic(n)%htotal(isd:ied,jsd:jed) )
  allocate( biotic(n)%alpha(isd:ied,jsd:jed) )
  allocate( biotic(n)%csat(isd:ied,jsd:jed) )
  allocate( biotic(n)%csat_csurf(isd:ied,jsd:jed) )
  allocate( biotic(n)%csat_acsurf(isd:ied,jsd:jed) )
  allocate( biotic(n)%csurf(isd:ied,jsd:jed) )
  allocate( biotic(n)%acsurf(isd:ied,jsd:jed) )
  allocate( biotic(n)%dpco2(isd:ied,jsd:jed) )
  allocate( biotic(n)%pco2atm(isd:ied,jsd:jed) )
  allocate( biotic(n)%paco2atm(isd:ied,jsd:jed) )
  allocate( biotic(n)%pco2surf(isd:ied,jsd:jed) )
  allocate( biotic(n)%paco2surf(isd:ied,jsd:jed) )
  allocate( biotic(n)%caco3_sediment(isd:ied,jsd:jed) )
  allocate( biotic(n)%det_sediment(isd:ied,jsd:jed) )
  allocate( biotic(n)%caco3_sed_remin(isd:ied,jsd:jed) )
  allocate( biotic(n)%det_sed_remin(isd:ied,jsd:jed) )
  allocate( biotic(n)%caco3_sed_depst(isd:ied,jsd:jed) )
  allocate( biotic(n)%det_sed_depst(isd:ied,jsd:jed) )
  allocate( biotic(n)%sio2(isd:ied,jsd:jed) )
  allocate( biotic(n)%po4(isd:ied,jsd:jed) )

  ntr_bgc = biotic(n)%ntr_bgc
  allocate(biotic(n)%bgc_global(1:ntr_bgc) )
  allocate(biotic(n)%vstf(isc:iec,jsc:jec,1:ntr_bgc) )

enddo  !}

! allocate biotic parameters

allocate( alphabio(isd:ied,jsd:jed) )
allocate( parbio(isd:ied,jsd:jed) )
allocate( kwbio(isd:ied,jsd:jed) )
allocate( kcbio(isd:ied,jsd:jed) )
allocate( abio(isd:ied,jsd:jed) )
allocate( bbio(isd:ied,jsd:jed) )
allocate( cbio(isd:ied,jsd:jed) )
allocate( k1bio(isd:ied,jsd:jed) )
allocate( muepbio(isd:ied,jsd:jed) )
allocate( muepsbio(isd:ied,jsd:jed) )
allocate( gam1bio(isd:ied,jsd:jed) )
allocate( gbio(isd:ied,jsd:jed) )
allocate( epsbio(isd:ied,jsd:jed) )
allocate( muezbio(isd:ied,jsd:jed) )
allocate( gam2bio(isd:ied,jsd:jed) )
allocate( muedbio(isd:ied,jsd:jed) )
allocate( muecaco3(isd:ied,jsd:jed) )
allocate( muedbio_sed(isd:ied,jsd:jed) )
allocate( muecaco3_sed(isd:ied,jsd:jed) )
allocate( wdetbio(isd:ied,jsd:jed) )
allocate( wcaco3(isd:ied,jsd:jed) )
allocate( nat_co2(isd:ied,jsd:jed) )
allocate( tscav_fe(isd:ied,jsd:jed) )
allocate( fe_bkgnd(isd:ied,jsd:jed) )
allocate( f_inorg(isd:ied,jsd:jed) )


!       initialize some arrays

do n = 1, instances  !{
 biotic(n)%htotal(:,:)  = 1.e-8
  biotic(n)%sio2(:,:) = 35. *1e-3
  biotic(n)%bgc_global(:) = 0.  ! this will make vstf zero
  biotic(n)%sal_global = 35.
  biotic(n)%vstf(:,:,:) = 0.
  biotic(n)%caco3_sediment(:,:) = 0.0
  biotic(n)%det_sediment(:,:) = 0.0
enddo  !} n

return
end subroutine  allocate_arrays  !}
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_bbc">
!
! <DESCRIPTION>
!     compute bottom boundary conditions e.g. fluxes due to
!     interaction with the sediment 
! </DESCRIPTION>
!

subroutine csiro_bgc_bbc (isc, iec, jsc, jec, T_prog, grid, Time)  !{

integer, intent(in)                             :: isc, iec
integer, intent(in)                             :: jsc, jec
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time


!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_bbc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer  :: i, j, n, k, nn
integer  :: ind_temp
real     :: fbc 
logical  :: used

! =====================================================================
!     begin executable code
! =====================================================================

!
!! rjm -- relax Fe to 1 at ocean bottom with 1 day time-constant
!! Huh? This is setting it to zero. RASF
!! Should we delete commented out code?       RASF
!! Yes, fe%btf() values are left as zero...
!!  bottom values for fe are assigned in csiro_bgc_tracer subroutine...
!!  so commenting out this block.  mac, nov10.
!
!do n = 1, instances  !{
! if (id_fe.ne.0) then 
!    ind_fe  = biotic(n)%ind_bgc(id_fe)
!    do j = jsc, jec  !{
!        do i = isc, iec  !{
!           k = grid%kmt(i,j)
!           t_prog(ind_fe)%btf(i,j) = rho0 * -1./86400. *0 ! negative into the ocean
!!grid%tmask(i,j,k) &
!!            / (-1.)  *            &
!!         (1  - max(0.0,t_prog(ind_fe)%field(i,j,k,taum1)))
!     enddo  !} i
!    enddo  !} j
!  endif
!enddo  !} n
!

!! Trial fe source within csiro_bgc_bbc. mac, nov12.
!do n = 1, instances  !{
! if (id_fe.ne.0) then 
!    do j = jsc, jec  !{
!        do i = isc, iec  !{
!           k = grid%kmt(i,j)
!           if (grid%zw(k) .le. 200) &
!t_prog(ind_fe)%btf(i,j) = -1.0 * rho0 * 1000. / 86400.
!!             t_prog(ind_fe)%btf(i,j) = -1.0 * rho0 * (0.999 - t_prog(ind_fe)%field(i,j,k,time%taum1)) * &
!!               grid%dzt(k) / (0.1 * 86400)  ! flux ~ density * diff(concentration) * dz / time_scale
!           ! setting time scale initial to 0.1 of a day, since at present, the concentration of the bottom is fixed in csiro_bgc_tracer;
!           !  and want the time scale to be shorter (stronger) than the scavenging to 0.6 within csiro_bgc_source.  
!           ! mac, nov12.  
!     enddo  !} i
!    enddo  !} j
!  endif
!enddo  !} n





 ! find remineralisation rate of sediment tracers.  mac, nov12.  
 call time_interp_external(muedbio_sed_id, time%model_time, muedbio_sed)
 call time_interp_external(muecaco3_sed_id, time%model_time, muecaco3_sed)
 ind_temp = fm_get_index('/ocean_mod/prog_tracers/temp')

 do n = 1, instances  !{
   do j = jsc, jec  !{
     do i = isc, iec  !{
       k = grid%kmt(i,j)
       if (k .gt. 0) then
         fbc= bbio(i,j)**(cbio(i,j)*T_prog(ind_temp)%field(i,j,k,time%taum1))
         biotic(n)%det_sed_remin(i,j) = muedbio_sed(i,j)*fbc*biotic(n)%det_sediment(i,j)
         biotic(n)%caco3_sed_remin(i,j) = muecaco3_sed(i,j)*fbc*biotic(n)%caco3_sediment(i,j)
         
! remineralisation of sediments to supply nutrient fields.  
! NB, btf values are positive from the water column into the sediment.  mac, nov12.  
         T_prog(ind_no3)%btf(i,j) = -1.0 * rho0 * biotic(n)%det_sed_remin(i,j)
         if (id_o2 .ne. 0) &
           T_prog(ind_o2)%btf(i,j)  = -172./16. * T_prog(ind_no3)%btf(i,j)
         if (id_dic .ne. 0) &
           T_prog(ind_dic)%btf(i,j)  = 106./16. * T_prog(ind_no3)%btf(i,j) - &
             rho0 * biotic(n)%caco3_sed_remin(i,j)
         if (id_adic .ne. 0) &
           T_prog(ind_adic)%btf(i,j)  = T_prog(ind_dic)%btf(i,j)
         if (id_fe .ne. 0) &
           T_prog(ind_fe)%btf(i,j)  = 2.0e-2 * T_prog(ind_no3)%btf(i,j)
         if (id_alk .ne. 0) &
           T_prog(ind_alk)%btf(i,j)  = -2.0 * rho0 * biotic(n)%caco3_sed_remin(i,j) - &
             T_prog(ind_no3)%btf(i,j)

       endif !} k > 0
     enddo  !} i
   enddo  !} j
 enddo  !} n
 
 

 ! send rate of remineralisation of sediment tracers to output. mac, nov12.
 do n = 1, instances  !{
  if (id_caco3_sed_remin .gt. 0) then
     used = send_data(id_caco3_sed_remin, biotic(n)%caco3_sed_remin(isc:iec,jsc:jec),          &
        time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (id_det_sed_remin .gt. 0) then
     used = send_data(id_det_sed_remin, biotic(n)%det_sed_remin(isc:iec,jsc:jec),          &
        time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
 enddo

do n = 1, instances  !{
 do nn=1,biotic(n)%ntr_bgc
! save the tracer bottom fluxes
  if (biotic(n)%id_bgc_btf(nn) .gt. 0) then
    used = send_data(biotic(n)%id_bgc_btf(nn),                        &
         t_prog(biotic(n)%ind_bgc(nn))%btf(isc:iec,jsc:jec)/rho0,                    &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
 enddo  !} nn
enddo  !} n

 
return
end subroutine  csiro_bgc_bbc  !}
! </SUBROUTINE> NAME="csiro_bgc_bbc"

!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>
!

subroutine csiro_bgc_end(isc, iec, jsc, jec, taup1, Thickness, T_prog, grid)  !{ 

type(ocean_thickness_type), intent(in) :: Thickness
integer, intent(in)                             :: isc
integer, intent(in)                             :: iec
integer, intent(in)                             :: jsc
integer, intent(in)                             :: jec
integer, intent(in)                             :: taup1
type(ocean_prog_tracer_type), dimension(:), intent(in)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid

!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_end'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: i
integer :: j
integer :: k
integer :: lun
integer :: n
real    :: total_det
real    :: total_zoo
real    :: total_phy
real    :: total_o2
real    :: total_no3

! =====================================================================
!     begin executable code
! =====================================================================

!       integrate the total concentrations of some tracers
!       for the end of the run

total_no3      = 0.0
total_phy      = 0.0
total_o2       = 0.0
total_zoo      = 0.0
total_det      = 0.0

do n = 1, instances  !{
  do k = 1,grid%nk !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        total_no3 = total_no3 +                                 &
             t_prog(biotic(n)%ind_bgc(1))%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,taup1)
        total_phy = total_phy +                                 &
             t_prog(biotic(n)%ind_bgc(2))%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,taup1)
        total_o2 = total_o2 +                                   &
             t_prog(biotic(n)%ind_bgc(3))%field(i,j,k,taup1) *      &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,taup1)
        total_zoo = total_zoo +                                 &
             t_prog(biotic(n)%ind_bgc(4))%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,taup1)
        total_det = total_det +                   &
             t_prog(biotic(n)%ind_bgc(5))%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,taup1)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  call mpp_sum(total_no3)
  call mpp_sum(total_phy)
  call mpp_sum(total_o2)
  call mpp_sum(total_zoo)
  call mpp_sum(total_det)

  write (stdout(),*) '  Instance ', trim(biotic(n)%name)
  write (stdout(),                                              &
       '(/'' Total NO3  = '',es19.12,'' Gmol N'')')     &
       total_no3 * 1.0e-12
  write (stdout(),                                              &
       '(/'' Total PHY  = '',es19.12,'' Gmol N'')')         &
       total_phy * 1.0e-12
  write (stdout(),                                              &
       '(/'' Total O2  = '',es19.12,'' Gmol O2'')')          &
       total_o2 * 1.0e-12
  write (stdout(),                                              &
       '(/'' Total ZOO  = '',es19.12,'' Gmol N'')')         &
       total_zoo * 1.0e-12
  write (stdout(),                                              &
       '(/'' Total DET  = '',es19.12,'' Gmol N'')')     &
       total_det * 1.0e-12
  write (stdout(),                                              &
       '(/'' Total N  = '',es19.12,'' Gmol N'')') &
       (total_no3+total_phy+total_zoo+total_det) * 1.0e-12

!
!  Write out extra restart file(s)
!

  call reset_field_pointer(sed_restart, id_restart(1), biotic(n)%caco3_sediment(:,:))
  call reset_field_pointer(sed_restart, id_restart(2), biotic(n)%det_sediment(:,:))

  call save_restart(sed_restart)

enddo  !} n

return
end subroutine  csiro_bgc_end  !}
! </SUBROUTINE> NAME="csiro_bgc_end"


!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine csiro_bgc_sbc(isc, iec, jsc, jec, isd, ied, jsd, jed, &
    T_prog, aice, wnd, grid, time, use_waterflux, salt_restore_as_salt_flux, atm_co2, co2flux, sfc_co2, iof_nit, iof_alg)

use ocmip2_co2calc_mod
use mpp_mod, only : mpp_sum
use time_interp_external_mod, only: time_interp_external
use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date
use field_manager_mod,    only: fm_get_index

integer, intent(in)                             :: isc, iec
integer, intent(in)                             :: jsc, jec
integer, intent(in)                             :: isd, ied
integer, intent(in)                             :: jsd, jed
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
real, intent(in), dimension(isd:ied,jsd:jed)                    :: aice, wnd
real, intent(in), dimension(isd:ied,jsd:jed), optional          :: iof_nit, iof_alg
logical, intent(in) :: use_waterflux, salt_restore_as_salt_flux

real, intent(in), dimension(isd:ied,jsd:jed), optional          :: atm_co2
real, intent(out), dimension(isd:ied,jsd:jed), optional         :: co2flux
real, intent(out), dimension(isd:ied,jsd:jed), optional         :: sfc_co2

real :: total_co2_flux, total_aco2_flux
logical :: used

!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_sbc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!       coefficients for O2 saturation
!

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

integer :: i
integer :: j
integer :: k
integer :: n
integer :: nn
integer :: ntr_bgc
integer :: ind_trc
integer :: m
integer :: kz
integer :: hour
integer :: day
real    :: days_in_this_year
integer :: minute
integer :: month
integer :: num_days
integer :: second
integer :: year
integer :: indtemp, indsal

! =====================================================================
!     begin executable code
! =====================================================================

  call get_date(time%model_time,                                &
                year, month, day, hour, minute, second)
  num_days = days_in_year(time%model_time)
  days_in_this_year = 0.0
  do m = 1, month - 1
    days_in_this_year = days_in_this_year +                     &
                       days_in_month(set_date(year, m, 1))
  enddo
  days_in_this_year = days_in_this_year + day - 1 + hour/24.0 + &
                      minute/1440.0 + second/86400.0

!---------------------------------------------------------------------
!     calculate interpolated xkw, seaice fraction and atmospheric
!       pressure & solar radiation
!---------------------------------------------------------------------
if (gasx_from_file) then
 call time_interp_external(pistonveloc_id, time%model_time, xkw_t)
else
 do j=jsc, jec
  do i=isc, iec
! the relations 0.31 * u^2 comes from Wanninkhof 1992, and is the coefficient for steady wind speed; the equation is 0.39 * u^2 for instaneous wind speeds. I don't know what exactly the timescale is between steady and instaneous...
! the 3.6e5 is a conversion of units; cm/hr to m/s.
! mac, may13.
   xkw_t(i,j)=(0.31 * wnd(i,j)**2.0 ) /3.6e5   
  enddo ! i
 enddo ! j
endif ! if (gasx_from_file)

call time_interp_external(nat_co2_id, time%model_time, nat_co2)
if (gasx_from_file) then
        call time_interp_external(atmpress_id, time%model_time, patm_t)
else !use the sea level pressure from the forcing (convert Pa to atm)
        !THIS HAS NOT BEEN IMPLEMENTED YET. READ INPUT FILE FOR NOW...
        call time_interp_external(atmpress_id, time%model_time, patm_t)
        !patm_t(isc:iec,jsc:jec) = patm(isc:iec,jsc:jec)/101325
endif
call time_interp_external(dust_id, time%model_time, dust_t)
if (id_adic .ne. 0) then
! The atmospheric co2 value for the anthropogenic+natural carbon tracer
! is either read from a file or a value from the access atmospheric model, 
! as determined by the flag use_access_co2.  mac, may13.  
 if (use_access_co2) then 
  aco2(isc:iec,jsc:jec) = atm_co2(isc:iec,jsc:jec)
 else
  call time_interp_external(aco2_id, time%model_time, aco2)
 endif
endif
if (ice_file4gasx) then
 call time_interp_external(seaicefract_id, time%model_time, fice_t)
else
 fice_t(isc:iec,jsc:jec) = aice(isc:iec,jsc:jec)
endif

indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
indsal = fm_get_index('/ocean_mod/prog_tracers/salt')

! -------------------------------------------------------------------
! rjm: start by zeroing the fluxes
if ((.not. use_waterflux) .or. (salt_restore_as_salt_flux)) then  !{
 do n = 1, instances  !{
  ntr_bgc = biotic(n)%ntr_bgc
    do nn=1,ntr_bgc
     ind_trc =  biotic(n)%ind_bgc(nn)
     t_prog(ind_trc)%stf(:,:) = 0     
    enddo
 enddo
endif

!---------------------------------------------------------------------
!  Compute the Schmidt number of CO2 in seawater using the 
!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!  7373-7382).
!---------------------------------------------------------------------
!
! no point doing the following calculations if no dic and o2 tracer
if (id_dic .eq. 0 .and. id_o2.eq. 0) return

do j = jsc, jec  !{
  do i = isc, iec  !{
    sc_co2(i,j) = 2073.1 + t_prog(indtemp)%field(i,j,1,time%taum1) * &
         (-125.62 + t_prog(indtemp)%field(i,j,1,time%taum1) *        &
          (3.6276 + t_prog(indtemp)%field(i,j,1,time%taum1) *        &
           (-0.043219))) * grid%tmask(i,j,1)
  enddo  !} i
enddo  !} j 

!---------------------------------------------------------------------
!  Compute the Schmidt number of O2 in seawater using the 
!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!  Cycles, 12, 141-163).
!---------------------------------------------------------------------

do j = jsc, jec  !{
  do i = isc, iec  !{
    sc_o2(i,j) = 1638.0 + t_prog(indtemp)%field(i,j,1,time%taum1) *  &
         (-81.83 + t_prog(indtemp)%field(i,j,1,time%taum1) *         &
          (1.483 + t_prog(indtemp)%field(i,j,1,time%taum1) *         &
           (-0.008004))) * grid%tmask(i,j,1)
  enddo  !} i
enddo  !} j 


!---------------------------------------------------------------------
!     calculate csurf, csat and csat - csurf via the routine co2calc
!        input and output units are in mol/m^3
!---------------------------------------------------------------------

!   determine the atmospheric pCO2 concetration for this time-step
do n = 1, instances  !{
  do j = jsc, jec
   do i = isc, iec
    biotic(n)%pco2atm(i,j) = nat_co2(i,j)
   enddo
  enddo
enddo  !} n

if (id_adic .ne. 0) then
 do n = 1, instances  !{
  do j = jsc, jec
   do i = isc, iec
    biotic(n)%paco2atm(i,j) = aco2(i,j)
   enddo
  enddo
 enddo  !} n
endif

!call init_ocmip2_co2calc(                                       &
!     time%model_time, isc, iec, jsc, jec, 1,                    &
!     t_prog(indtemp)%field(isc:iec,jsc:jec,1,time%taum1),            &
!     t_prog(indsal)%field(isc:iec,jsc:jec,1,time%taum1))

do n = 1, instances  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      htotallo(i,j) =  .1
      htotalhi(i,j) =  100.0
    enddo  !} i
  enddo  !} j 

! ocmip2+co2cal expects tracers in mol/m3 !!!!
  if (id_dic.ne.0) then
    biotic(n)%po4(:,:) = t_prog(ind_dic)%field(isd:ied,jsd:jed,1,time%taum1)*0
    if (id_no3.ne.0) biotic(n)%po4(:,:) = t_prog(ind_no3)%field(isd:ied,jsd:jed,1,time%taum1)/16.*1e-3
    if (id_po4.ne.0) biotic(n)%po4(:,:) = t_prog(ind_po4)%field(isd:ied,jsd:jed,1,time%taum1)*1e-3

    call ocmip2_co2calc(isd, jsd, &
       isc, iec, jsc, jec,                     &
       grid%tmask(isd:ied,jsd:jed,1),                           &
       t_prog(indtemp)%field(isd:ied,jsd:jed,1,time%taum1),     &
       t_prog(indsal)%field(isd:ied,jsd:jed,1,time%taum1),     &
       t_prog(ind_dic)%field(isd:ied,jsd:jed,1,time%taum1)*1e-3,&
       t_prog(ind_alk)%field(isd:ied,jsd:jed,1,time%taum1)*1e-3,&
       biotic(n)%po4(isd:ied,jsd:jed),&
       biotic(n)%sio2(isd:ied,jsd:jed),                       &
       htotallo(isc:iec,jsc:jec), htotalhi(isc:iec,jsc:jec), &
       biotic(n)%htotal(isc:iec,jsc:jec),                      &
       biotic(n)%csurf(isc:iec,jsc:jec),                    &
       alpha=biotic(n)%alpha(isc:iec,jsc:jec) ,                   &            
       pco2surf = biotic(n)%pco2surf(isc:iec,jsc:jec),            &
       scale= 1.0/1024.5 )

    if (id_adic .ne. 0) then ! calculate CO2 flux including anthropogenic CO2
      call ocmip2_co2calc(isd, jsd,                     &
       isc, iec, jsc, jec,                    &
       grid%tmask(isd:ied,jsd:jed,1),                           &
       t_prog(indtemp)%field(isd:ied,jsd:jed,1,time%taum1),     &
       t_prog(indsal)%field(isd:ied,jsd:jed,1,time%taum1),     &
       t_prog(ind_adic)%field(isd:ied,jsd:jed,1,time%taum1)*1e-3,&
       t_prog(ind_alk)%field(isd:ied,jsd:jed,1,time%taum1)*1e-3,&
       biotic(n)%po4(isd:ied,jsd:jed),&
       biotic(n)%sio2(isd:ied,jsd:jed),                       &
       htotallo(isc:iec,jsc:jec), htotalhi(isc:iec,jsc:jec), &
       biotic(n)%htotal(isc:iec,jsc:jec),                      &
       biotic(n)%acsurf(isc:iec,jsc:jec),                    &
       alpha=biotic(n)%alpha(isc:iec,jsc:jec) ,                   &            
       pco2surf = biotic(n)%paco2surf(isc:iec,jsc:jec),           &
       scale = 1.0/1024.5 )
    endif  ! no adic
  endif  ! no dic
enddo  !} n 


!---------------------------------------------------------------------
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
!  o2_saturation is defiend between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol/m^3
!---------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    tt(i) = 298.15 - t_prog(indtemp)%field(i,j,1,time%taum1)
    tk(i) = 273.15 + t_prog(indtemp)%field(i,j,1,time%taum1)
    ts(i) = log(tt(i) / tk(i))
    ts2(i) = ts(i) * ts(i)
    ts3(i) = ts2(i) * ts(i)
    ts4(i) = ts3(i) * ts(i)
    ts5(i) = ts4(i) * ts(i)
    o2_saturation(i,j) =                                        &
         exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                     &
             a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +             &
             t_prog(indsal)%field(i,j,1,time%taum1) *                &
             (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +       &
              c_0*t_prog(indsal)%field(i,j,1,time%taum1)))
  enddo  !} i
enddo  !} j 

!chd
!chd       convert from ml/l to mmol/m^3
!chd

do j = jsc, jec  !{
  do i = isc, iec  !{
    o2_saturation(i,j) = o2_saturation(i,j) *                   &
                        (1000.*1000.0/22391.6)
  enddo  !} i
enddo  !} j 

!
!---------------------------------------------------------------------
!     calculate piston-velocities
!      including  effect of sea-ice
!      xkw is given in cm/s (converted in read_biotic_bc), therefore
!      kw is also in cm/s
!---------------------------------------------------------------------

do j = jsc, jec  !{
  do i = isc, iec  !{
    kw_co2(i,j) = (1.0 - fice_t(i,j)) * xkw_t(i,j) *            &
                   sqrt(660.0/sc_co2(i,j)) * grid%tmask(i,j,1)

    kw_o2(i,j) = (1.0 - fice_t(i,j)) * xkw_t(i,j) *             &
                  sqrt(660.0/sc_o2(i,j)) * grid%tmask(i,j,1)
  enddo  !} i
enddo  !} j 

!
!chd---------------------------------------------------------------------
!chd     calculate surface fluxes for BIOTICs
!chd       kw is in m/s, o2_sat is in mmol/m^3,  
!chd       hence stf should be in  mmol/m^2/s
!chd---------------------------------------------------------------------
!

total_co2_flux = 0.0

if (id_dic.ne.0) then
  do n = 1, instances  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{ 
       ! These calculations used to be done in ocmip2_co2calc in mom4p1_2007, 
       !  but not with mom4p1_2009.  
       ! Some extra calculations now required in csiro_bgc_sbc for co2 air-sea 
       ! alpha is in units of mol/m3/atm (atm -> partial pressure of CO2), and 1e6 is to convert micro-atm to atm.
       !  flux.   mac, apr11.  
       biotic(n)%csat(i,j) = biotic(n)%pco2atm(i,j) / 1e6 * biotic(n)%alpha(i,j) * patm_t(i,j)
       biotic(n)%csat_csurf(i,j) = biotic(n)%csat(i,j) - biotic(n)%csurf(i,j)
       
       t_prog(ind_dic)%stf(i,j) = rho0 * kw_co2(i,j) *        &
        biotic(n)%csat_csurf(i,j)*1e3 !convert from  mol/m^3 to mmol/m^3

       total_co2_flux = total_co2_flux + kw_co2(i,j) *        &
        biotic(n)%csat_csurf(i,j) * 3.7843e-7 * grid%dat(i,j) * grid%tmask(i,j,1) ! convert from  mol/s to Pg/year (12.0*1e-15*86400*365=3.78e-7)
      enddo  !} i
    enddo  !} j 
  enddo  !} n 
endif

if (id_total_co2_flux .gt. 0) then
 call mpp_sum(total_co2_flux)
 used = send_data(id_total_co2_flux,total_co2_flux,Time%model_time)
endif


total_aco2_flux = 0.0 

if (id_adic.ne.0) then
  do n = 1, instances  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
       ! These calculations used to be done in ocmip2_co2calc in mom4p1_2007, 
       !  but not with mom4p1_2009.  
       ! Some extra calculations now required in csiro_bgc_sbc for co2 air-sea 
       !  flux.   mac, apr11.  
       ! Note, csat has been used as temporary variable, alpha is only f(T,S),
       !  csurf had been used as a temporary variable, but now needs to be saved
       !  for preindustrial vs anthropogenic.
       biotic(n)%csat(i,j) = biotic(n)%paco2atm(i,j) / 1e6 * biotic(n)%alpha(i,j) * patm_t(i,j)
       biotic(n)%csat_acsurf(i,j) = biotic(n)%csat(i,j) - biotic(n)%acsurf(i,j)
       
       t_prog(ind_adic)%stf(i,j) = rho0 * kw_co2(i,j) *        &
        biotic(n)%csat_acsurf(i,j)*1e3 !convert from  mol/m^3 to mmol/m^3

       total_aco2_flux = total_aco2_flux + kw_co2(i,j) *        &
        biotic(n)%csat_acsurf(i,j) * 3.7843e-7 * grid%dat(i,j) * grid%tmask(i,j,1) ! convert from  mol/s to Pg/year (12.0*1e-15*86400*365=3.78e-7)
! send the anthropogenic Pco2 and co2 flux into the ocean back to the atmospheric model.  mac, may13.  
!RASF avoid using Ocean_sfc. 
!       if(present(Ocean_sfc)) then
!          Ocean_sfc%co2flux(i,j) = kw_co2(i,j) *        &
!                                   biotic(n)%csat_acsurf(i,j)*0.04401 !convert from  mol/m^2/s to kg(CO2)/m^2/s, 0.04401 kg(CO2)/mole
!          Ocean_sfc%co2(i,j) = biotic(n)%paco2surf(i,j) 
!       endif
      enddo  !} i
    enddo  !} j 
    if(present(co2flux) .and. present(sfc_co2)) then
      do j = jsc, jec  !{
        do i = isc, iec  !{
          co2flux(i,j) = kw_co2(i,j) *        &
                                   biotic(n)%csat_acsurf(i,j)*0.04401 !convert from  mol/m^2/s to kg(CO2)/m^2/s, 0.04401 kg(CO2)/mole
          sfc_co2(i,j) = biotic(n)%paco2surf(i,j) 
        enddo  !} i
      enddo  !} j 
    endif
    
  enddo  !} n 
endif

if (id_total_aco2_flux .gt. 0) then
 call mpp_sum(total_aco2_flux)
 used = send_data(id_total_aco2_flux,total_aco2_flux,Time%model_time)
endif


if (id_o2.ne.0) then
  do n = 1, instances  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_o2)%stf(i,j) =   rho0 *     kw_o2(i,j) *        &
            (o2_saturation(i,j) * patm_t(i,j) -                     &
             t_prog(ind_o2)%field(i,j,1,time%taum1))
      enddo  !} i
    enddo  !} j 
  enddo  !} n 
endif

!rjm: Do atmospheric iron input
! for Fe units of nmol/m3, Fe input must be in nmol/m2/s
if (id_fe.ne.0) then
  do n = 1, instances  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_fe)%stf(i,j) =  rho0 * dust_t(i,j)
      enddo  !} i
    enddo  !} j
  enddo  !} n
endif

if (.not. use_waterflux)  then  !{
! rjm - One only needs to compute virtual fluxes if waterflux is not used
! NB, when this routine is called, t_prog(ind_sal)%stf() only has applied fluxes
!  and not any restoring fluxes.  
! csiro_bgc_virtual_fluxes() now calculated BGC virtual fluxes if any virtual flux 
!  is use to restore salinity.   mac, dec12.  
  do n = 1, instances  !{
    ntr_bgc = biotic(n)%ntr_bgc

    do nn=1,ntr_bgc
      ind_trc =  biotic(n)%ind_bgc(nn) 

      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf(i,j,nn) =                    &
               t_prog(indsal)%stf(i,j) *                &
               biotic(n)%bgc_global(nn) / biotic(n)%sal_global
          t_prog(ind_trc)%stf(i,j) =                  &
               t_prog(ind_trc)%stf(i,j) +             &
               biotic(n)%vstf(i,j,nn)
        enddo  !} i
      enddo  !} j
    enddo !} nn
  enddo  !} n 
endif  !}
   
return

end subroutine  csiro_bgc_sbc  !}
! </SUBROUTINE> NAME="csiro_bgc_sbc"

!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_virtual_fluxes">
!
! <DESCRIPTION>
!      Apply virtual fluxes to BGC tracers when sea surface salinity
!      is restored as a salt flux within ocean_sbc//flux_adjust()
!      This subroutine only adds the virtual flux associated with
!      any virtual salt fluxes from salinity restoring; any virtual fluxes
!      from applied forcing (i.e. when use_waterflux = .false.) are applied to 
!      bgc virtual fluxes with the csiro_bgc_sbc() subroutine.  
!      mac, dec12.
! </DESCRIPTION>

subroutine csiro_bgc_virtual_fluxes(isc, iec, jsc, jec, isd, ied, jsd, jed, salt_vstf, T_prog)  !{

!-----------------------------------------------------------------------
!       arguments
!-----------------------------------------------------------------------

integer, intent(in)                          :: isc, iec, jsc, jec, isd, ied, jsd, jed
real, intent(in), dimension(isd:ied,jsd:jed) :: salt_vstf    ! virtual salt flux.  
type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog


!-----------------------------------------------------------------------
!       local definitions
!-----------------------------------------------------------------------

integer  :: i, j, n, nn, ntr_bgc, ind_trc

!-----------------------------------------------------------------------
!       start executable code
!-----------------------------------------------------------------------

  do n = 1, instances  !{
    ntr_bgc = biotic(n)%ntr_bgc

    do nn=1,ntr_bgc
      ind_trc =  biotic(n)%ind_bgc(nn) 

      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf(i,j,nn) =             &
               salt_vstf(i,j) *                &
               biotic(n)%bgc_global(nn) / biotic(n)%sal_global
          t_prog(ind_trc)%stf(i,j) =                  &
               t_prog(ind_trc)%stf(i,j) +             &
               biotic(n)%vstf(i,j,nn)
        enddo  !} i
      enddo  !} j
    enddo !} nn
  enddo  !} n 


end subroutine  csiro_bgc_virtual_fluxes  !}
! </SUBROUTINE> NAME="csiro_bgc_virtual_fluxes"


!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine csiro_bgc_init  !{

!-----------------------------------------------------------------------
!       local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=128)              :: version = "git:/mom_mac/src/mom5csiro_bgc.F90"
character(len=128)              :: tagname = "mac, sep16."

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_string_len)                            :: string
character(len=fm_field_name_len+3)                      :: long_suffix
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

integer                                                 :: nn, ntr_bgc
real                                                    :: min_range=0.0, max_range=1.e4
real                                                    :: sum_ntr = 0.0
character(len=8)                                        :: bgc_trc='tracer00'

!       Initialize the csiro_bgc package

!package_index = otpm_set_tracer_package(package_name,            &
!     file_in = default_file_in, file_out = default_file_out,     &
!     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
package_index = otpm_set_tracer_package(package_name,            &
       restart_file = default_file_in,      &
      caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!       Check whether to use this package

path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif  !}

!       Check some things

if (instances .eq. 0) then  !{
  write (stdout(),*) trim(note_header), ' No instances'
  do_csiro_bgc = .false.
else  !}{
  if (instances .eq. 1) then  !{
    write (stdout(),*) trim(note_header), ' ', instances, ' instance'
  else  !}{
    write (stdout(),*) trim(note_header), ' ', instances, ' instances'
    call mpp_error(FATAL, trim(error_header) // ' CSIRO_bgc code not yet ready for instances > 1, mac, nov12')
  endif  !}
  do_csiro_bgc = .true.
endif  !}

!       Return if we don't want to use this package,
!       after changing the list back

if (.not. do_csiro_bgc) then  !{
  return
endif  !}

!       allocate storage for biotic array
allocate ( biotic(instances) )

!       loop over the names, saving them into the biotic array

do n = 1, instances  !{
  if (fm_get_value(path_to_names, name, index = n)) then  !{
    biotic(n)%name = name
  else  !}{
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif  !}
enddo  !}

 
!-----------------------------------------------------------------------
!       Process the namelists
!-----------------------------------------------------------------------

!       Add the package name to the list of good namelists, to be used
!       later for a consistency check

if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif  !}

!-----------------------------------------------------------------------
!       Set up the *global* namelist
!-----------------------------------------------------------------------

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

  call fm_util_set_value('atmpress_file', 'INPUT/atmpress_ocmip2.nc')
  call fm_util_set_value('atmpress_name', 'atmpress')
  call fm_util_set_value('pistonveloc_file', 'INPUT/pistonveloc_ocmip2.nc')
  call fm_util_set_value('pistonveloc_name', 'pistonveloc')
  call fm_util_set_value('aco2_file', 'INPUT/co2_b35_2D.nc')
  call fm_util_set_value('aco2_name', 'co2')
  call fm_util_set_value('seaicefract_file', 'INPUT/f_ice_ocmip2.nc')
  call fm_util_set_value('seaicefract_name', 'f_ice')
  call fm_util_set_value('dust_file', 'INPUT/dust.nc')
  call fm_util_set_value('dust_name', 'DUST')

! additional information
  call fm_util_set_value('s_npp', -.1)           ! scale factor for NP
  call fm_util_set_value('qbio_model', 'bio_vx') ! version name 
  call fm_util_set_value('bio_version',1)        ! version of the bgc module
  call fm_util_set_value('zero_floor', .false.)  ! apply hard floor to bgc tracers
  call fm_util_set_value('sw_thru_ice', .true.)  ! is shortwave flux modified by an ice model?
  call fm_util_set_value('gasx_from_file', .true.)! use file with gas exchange coefficients?
  call fm_util_set_value('ice_file4gasx', .true.)! use file with ice cover for gas exchange?
! Set defaults to use the ACCESS atmospheric CO2 for the anthropogenic CO2 in the ocean if this is compiled for ACCESS-CM/ESM, otherwise to read CO2 from files provided.  mac, may13.  
#if defined(ACCESS_CM) 
  call fm_util_set_value('use_access_co2', .true.)! use access model co2 values.  mac, may13.  
#else
  call fm_util_set_value('use_access_co2', .false.)! use file for atmospheric co2 values.  
#endif
  call fm_util_set_value('id_po4',0)      
  call fm_util_set_value('id_dic',0)      
  call fm_util_set_value('id_adic',0)      
  call fm_util_set_value('id_alk',0)      
  call fm_util_set_value('id_o2',0)      
  call fm_util_set_value('id_no3',0)      
  call fm_util_set_value('id_phy',0)      
  call fm_util_set_value('id_zoo',0)      
  call fm_util_set_value('id_det',0)      
  call fm_util_set_value('id_caco3',0)      
  call fm_util_set_value('id_fe',0)      

  atmpress_file      =  fm_util_get_string ('atmpress_file', scalar = .true.)
  atmpress_name      =  fm_util_get_string ('atmpress_name', scalar = .true.)
  pistonveloc_file   =  fm_util_get_string ('pistonveloc_file', scalar = .true.)
  pistonveloc_name   =  fm_util_get_string ('pistonveloc_name', scalar = .true.)
  aco2_file          =  fm_util_get_string ('aco2_file', scalar = .true.)
  aco2_name          =  fm_util_get_string ('aco2_name', scalar = .true.)
  seaicefract_file   =  fm_util_get_string ('seaicefract_file', scalar = .true.)
  seaicefract_name   =  fm_util_get_string ('seaicefract_name', scalar = .true.)
  dust_file   =  fm_util_get_string ('dust_file', scalar = .true.)
  dust_name   =  fm_util_get_string ('dust_name', scalar = .true.)

  qbio_model   =  fm_util_get_string ('qbio_model', scalar = .true.)
  s_npp   =  fm_util_get_real ('s_npp', scalar = .true.)
  bio_version   =  fm_util_get_integer ('bio_version', scalar = .true.)
  ! rjm: Tracer to use
  zero_floor = fm_util_get_logical ('zero_floor', scalar = .true.)
  sw_thru_ice = fm_util_get_logical ('sw_thru_ice', scalar = .true.)
  gasx_from_file = fm_util_get_logical ('gasx_from_file', scalar = .true.)
  ice_file4gasx = fm_util_get_logical ('ice_file4gasx', scalar = .true.)
  use_access_co2 = fm_util_get_logical ('use_access_co2', scalar = .true.)

  id_dic   =   fm_util_get_integer ('id_dic', scalar = .true.)
  id_adic  =   fm_util_get_integer ('id_adic', scalar = .true.)
  id_alk   =   fm_util_get_integer ('id_alk', scalar = .true.)
  id_po4   =   fm_util_get_integer ('id_po4', scalar = .true.)
  id_o2    =   fm_util_get_integer ('id_o2', scalar = .true.)
  id_no3   =   fm_util_get_integer ('id_no3', scalar = .true.)
  id_zoo   =   fm_util_get_integer ('id_zoo', scalar = .true.)
  id_phy   =   fm_util_get_integer ('id_phy', scalar = .true.)
  id_det   =   fm_util_get_integer ('id_det', scalar = .true.)
  id_caco3 =   fm_util_get_integer ('id_caco3', scalar = .true.)
  id_fe    =   fm_util_get_integer ('id_fe', scalar = .true.)


call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)


sum_ntr = min(1,id_dic)+min(1,id_adic)+min(1,id_po4)+min(1,id_alk) &
   +min(1,id_o2)+min(1,id_no3)+min(1,id_phy)+min(1,id_zoo)+min(1,id_det) &
   +min(1,id_caco3)+min(1,id_fe)
if (mpp_pe() == mpp_root_pe() ) print*,'csiro_bgc_init: Number bgc tracers = ',sum_ntr


!-----------------------------------------------------------------------
!       Set up the prognostic tracers.
!-----------------------------------------------------------------------

! RASF: Use sensible names for tracers
! For default case names will be simply adic, o2 etc
! If we have multiple instances the names will take the form
! adic_instancename, o2_instancename
! Different input files etc can be set in the field_table like for temp, salt

do n = 1, instances  !{

  biotic(n)%ntr_bgc = min(1,id_dic)+min(1,id_adic)+min(1,id_po4)+min(1,id_alk) &
     +min(1,id_o2)+min(1,id_no3)+min(1,id_phy)+min(1,id_zoo)+min(1,id_det) &
     +min(1,id_caco3)+min(1,id_fe)
  if (mpp_pe() == mpp_root_pe() ) print*,'Number bgc tracers = ',biotic(n)%ntr_bgc
      

  name = biotic(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ''
    long_suffix = ''
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

  ntr_bgc=biotic(n)%ntr_bgc

 do nn=1,ntr_bgc
    if ( nn == id_no3 ) then
          bgc_trc='no3'
          min_range=-1e-5
          max_range=100.0
    else if ( nn ==  id_phy ) then
          bgc_trc='phy'
          min_range=-1e-7
          max_range=10.0
    else if (nn == id_o2 ) then
          bgc_trc='o2'
          min_range=-1e-7
          max_range=10.0
    else if (nn == id_zoo ) then
          bgc_trc='zoo'
          min_range=-1e-7
          max_range=10.0
    else if (nn == id_det ) then
          bgc_trc='det'
          min_range=0.0
          max_range=10.0
    else if (nn == id_caco3 ) then
          bgc_trc='caco3'
          min_range=-1e-5
          max_range=10.0
    else if (nn == id_dic ) then
          bgc_trc='dic'
          min_range=0.0
          max_range=3000.0
    else if (nn == id_alk ) then
          bgc_trc='alk'
          min_range=0.0
          max_range=3000.0
    else if (nn == id_adic ) then
          bgc_trc='adic'
          min_range=-1e-5
          max_range=3000.0
    else if (nn == id_fe ) then
          bgc_trc='fe'
          min_range=-0.0001
          max_range=10.0
    else if (nn == id_po4 ) then
          bgc_trc='po4'
          min_range=0.0
          max_range=100.0
    else 
          bgc_trc(1:6)='dummy'
          write(bgc_trc(7:8),'(i2.2)') nn
          min_range=0.0
          max_range=100.0
    endif 
          
    
    biotic(n)%ind_bgc(nn) = otpm_set_prog_tracer(trim(bgc_trc) // trim(suffix),     &
       package_name,                                                  &
       longname = trim(bgc_trc) // trim(long_suffix),                      &
       units = 'mmol/m^3', flux_units = 'mmol/m^2/s',                 &
       caller = trim(mod_name)//'('//trim(sub_name)//')',              &
       min_range=min_range,max_range=max_range)
 enddo  !} nn
enddo  !} n

!-----------------------------------------------------------------------
!       Set up the instance namelists
!-----------------------------------------------------------------------


do n = 1, instances  !{

!       create the instance namelist

  call fm_util_start_namelist(package_name, trim(biotic(n)%name) , caller = caller_str, no_overwrite = .true., &
       check = .true.)


  call fm_util_set_value('file_in', default_file_in)
  call fm_util_set_value('file_out', default_file_out)
  call fm_util_set_value('init', .false.)

  call fm_util_end_namelist(package_name, biotic(n)%name, check = .true., caller = caller_str)

enddo  !} n

!       Check for any errors in the number of fields in the namelists for this package

good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then  !{
  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif  !}


! possible tracers in model
do n = 1, instances  !{

   ind_dic = biotic(n)%ind_bgc(id_dic)
   ind_adic = biotic(n)%ind_bgc(id_adic)
   ind_alk = biotic(n)%ind_bgc(id_alk)
   if(id_po4 .ne. 0 ) ind_po4 = biotic(n)%ind_bgc(id_po4)
   ind_o2  = biotic(n)%ind_bgc(id_o2)

   ind_no3  = biotic(n)%ind_bgc(id_no3)

   ind_zoo  = biotic(n)%ind_bgc(id_zoo)
   ind_phy  = biotic(n)%ind_bgc(id_phy)
   ind_det  = biotic(n)%ind_bgc(id_det)

   ind_caco3= biotic(n)%ind_bgc(id_caco3)
   if(id_fe .ne. 0 ) ind_fe= biotic(n)%ind_bgc(id_fe)

enddo  !} n

! Write SVN version numbers to both the logfile and the standard output
! mac, mar12.  
call write_version_number(version, tagname)

if (mpp_pe() == mpp_root_pe() ) then
 print*
 print*,"Repo. details of csiro_bgc compiled into executable"
 print*,"git:/mom_mac/src/mom5csiro_bgc.F90"
 print*,"mac, sep16."
 print*
endif

id_clock_csiro_obgc = mpp_clock_id('(CSIRO ocean biogeochemistry)',grain=CLOCK_ROUTINE)

return

end subroutine csiro_bgc_init  !}
! </SUBROUTINE> NAME="csiro_bgc_init"


!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!

subroutine csiro_bgc_source(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, grid, time, dtts, Thickness, Dens, swflx, sw_frac_zt)  !{

use field_manager_mod,        only: fm_get_index

!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc, iec
integer, intent(in)                                             :: jsc, jec
integer, intent(in)                                             :: isd, ied
integer, intent(in)                                             :: jsd, jed
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
real, intent(in)                                                :: dtts
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, intent(in), dimension(isd:ied,jsd:jed)                    :: swflx        ! short wave radiation flux (W/m^2)
real, intent(in), dimension(isd:,jsd:,:)                        :: sw_frac_zt        ! fraction of short wave radiation flux (none)


!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_source'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: i
integer :: j
integer :: k
integer :: n
logical :: used
integer :: day
integer :: month
integer :: year
integer :: hour
integer :: minute
integer :: second
integer :: indtemp

integer :: nn, ntr_bgc

! include "bio_v1.par"

real    :: total_trc
!
! =====================================================================
!     begin executable code
! =====================================================================

!       get the model month

call mpp_clock_begin(id_clock_csiro_obgc)

call get_date(time%model_time, year, month, day,                &
              hour, minute, second)

!-----------------------------------------------------------------------
!     calculate the source terms for BIOTICs
!-----------------------------------------------------------------------


 select case(bio_version)
case(0)
! call bio_v0a(Thickness)
!include "bio_v0.f90"

case(1)
! call bio_v1a(Thickness)
!include "bio_v1.f90"

case(2)
!include "bio_v2.f90"

case(3)
call bio_v3(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, grid, time, dtts, Thickness, Dens, swflx, sw_frac_zt)

case default
if (mpp_pe() == mpp_root_pe() )print*,'rjm do no bio version'

end select

!-----------------------------------------------------------------------
!       Save variables for diagnostics
!-----------------------------------------------------------------------
if (id_pco2 .gt. 0) then
  used = send_data(id_pco2, biotic(1)%pco2surf(isc:iec,jsc:jec),            &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_paco2 .gt. 0) then
  used = send_data(id_paco2, biotic(1)%paco2surf(isc:iec,jsc:jec),            &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_co2_sat .gt. 0) then
  used = send_data(id_co2_sat, biotic(1)%csat_csurf(isc:iec,jsc:jec),   &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_aco2_sat .gt. 0) then
  used = send_data(id_aco2_sat, biotic(1)%csat_acsurf(isc:iec,jsc:jec),   &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

if (id_sc_o2 .gt. 0) then
  used = send_data(id_sc_o2, sc_o2(isc:iec,jsc:jec),            &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_o2_sat .gt. 0) then
  used = send_data(id_o2_sat, o2_saturation(isc:iec,jsc:jec),   &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_kw_o2 .gt. 0) then
  used = send_data(id_kw_o2, kw_o2(isc:iec,jsc:jec),            &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

! Gross production of phytoplankton

if (id_pprod_gross .gt. 0) then
  used = send_data(id_pprod_gross, pprod_gross(isc:iec,jsc:jec,:),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
endif

if (id_pprod_gross_2d .gt. 0) then
  pprod_gross_2d(:,:)=0.0
  do k=1,grid%nk
     do j=jsc,jec
        do i=isc,iec
           pprod_gross_2d(i,j)=pprod_gross_2d(i,j) + pprod_gross(i,j,k)*Thickness%dzt(i,j,k)
        enddo
     enddo
  enddo
  used = send_data(id_pprod_gross_2d, pprod_gross_2d(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

!det export at 100 m
if (id_wdet100 .gt. 0) then
  wdet100(:,:) = wdetbio(isc:iec,jsc:jec)*t_prog(ind_det)%field(isc:iec,jsc:jec,minloc(grid%zt(:)-100,dim=1),time%taum1)
  used = send_data(id_wdet100, wdet100(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

! Net primary productivity

! at each depth
if (id_npp3d .gt. 0) then
  used = send_data(id_npp3d, npp3d(isc:iec,jsc:jec,:),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
endif
! depth integrated
if (id_npp2d .gt. 0) then
  npp2d(:,:)=0.0
  do k=1,grid%nk
        npp2d(isc:iec,jsc:jec) = npp2d(isc:iec,jsc:jec) + npp3d(isc:iec,jsc:jec,k)*Thickness%dzt(isc:iec,jsc:jec,k)
  enddo
  used = send_data(id_npp2d, npp2d(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
! at surface
if (id_npp1 .gt. 0) then
  used = send_data(id_npp1, npp3d(isc:iec,jsc:jec,1),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

! Gross production of zooplankton

if (id_zprod_gross .gt. 0) then
  used = send_data(id_zprod_gross, zprod_gross(isc:iec,jsc:jec,:),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
endif

! growth limit of phytoplankton for light.

if (id_light_limit .gt. 0) then
  used = send_data(id_light_limit, light_limit(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

! mixed-layer-integrated quantities

if (id_adic_intmld .gt. 0) then
  used = send_data(id_adic_intmld, adic_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_dic_intmld .gt. 0) then
  used = send_data(id_dic_intmld, dic_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_o2_intmld .gt. 0) then
  used = send_data(id_o2_intmld, o2_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_no3_intmld .gt. 0) then
  used = send_data(id_no3_intmld, no3_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_fe_intmld .gt. 0) then
  used = send_data(id_fe_intmld, fe_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_phy_intmld .gt. 0) then
  used = send_data(id_phy_intmld, phy_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_det_intmld .gt. 0) then
  used = send_data(id_det_intmld, det_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_pprod_gross_intmld .gt. 0) then
  used = send_data(id_pprod_gross_intmld, pprod_gross_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_npp_intmld .gt. 0) then
  used = send_data(id_npp_intmld, npp_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_radbio_intmld .gt. 0) then
  used = send_data(id_radbio_intmld, radbio_intmld(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

!detritus sinking at 100 m
if (id_wdet100avg .gt. 0) then
  used = send_data(id_wdet100avg, wdet100avg(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

! PAR for phytoplankton at surface.

if (id_radbio1 .gt. 0) then
  used = send_data(id_radbio1, radbio3d(isc:iec,jsc:jec,1),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

! PAR for phytoplankton at all depths.

if (id_radbio3d .gt. 0) then
  used = send_data(id_radbio3d, radbio3d(isc:iec,jsc:jec,:),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
endif

do n = 1, instances  !{

 ntr_bgc=biotic(n)%ntr_bgc
 do nn=1,ntr_bgc

! save the tracer surface fluxes
  if (biotic(n)%id_bgc_stf(nn) .gt. 0) then
    used = send_data(biotic(n)%id_bgc_stf(nn),                        &
         t_prog(biotic(n)%ind_bgc(nn))%stf(isc:iec,jsc:jec)/rho0,                    &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
! save the tracer sources
  if (biotic(n)%id_bgc_src(nn) .gt. 0) then
    used = send_data(biotic(n)%id_bgc_src(nn),                        &
         t_prog(biotic(n)%ind_bgc(nn))%th_tendency(isc:iec,jsc:jec,:)/     &
         Thickness%rho_dzt(isc:iec,jsc:jec,:,time%tau),               &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
 enddo  !} nn
 
 ! rate of deposition of sinking tracers to sediment
 if (id_caco3_sed_depst .gt. 0) then
    used = send_data(id_caco3_sed_depst, biotic(n)%caco3_sed_depst(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
 endif
 if (id_det_sed_depst .gt. 0) then
    used = send_data(id_det_sed_depst, biotic(n)%det_sed_depst(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
 endif
 
 
enddo  !} n


call mpp_clock_end(id_clock_csiro_obgc)

return
end subroutine  csiro_bgc_source  !}
! </SUBROUTINE> NAME="csiro_bgc_source"


!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!

subroutine csiro_bgc_start (time, domain, grid) !{

use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date
use time_interp_external_mod, only: init_external_field
use diag_manager_mod, only: register_diag_field, diag_axis_init
use fms_mod, only : read_data

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

type(ocean_time_type), intent(in)                               :: Time
type(ocean_domain_type), intent(in)                             :: Domain
type(ocean_grid_type), intent(in)                               :: Grid

!
!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

character(len=fm_string_len)            :: file_in
character(len=fm_string_len)            :: file_out
integer                                 :: index
logical                                 :: init
integer                                 :: check
integer                                 :: m, n
integer                                 :: second
integer                                 :: minute
integer                                 :: hour
integer                                 :: day
integer                                 :: month
integer                                 :: year
character(len=fm_field_name_len)        :: name
real                                    :: days_in_this_year
integer                                 :: num_days
character(len=fm_field_name_len+1)      :: str
integer, dimension(:)                   :: tracer_axes_1(3)
integer                                 :: id_zt_1
character(len=256)                      :: caller_str

integer :: nn, ntr_bgc
character(len=5)        :: bgc_stf='stf00'
character(len=5)        :: bgc_btf='btf00'
character(len=5)        :: bgc_src='src00'
character(len=64)       :: name1, name2, name3, name4 


! =====================================================================
!       begin of executable code
! =====================================================================

write(stdout(),*) 
write(stdout(),*) trim(note_header),                     &
                  'Starting ', trim(package_name), ' module'

!-----------------------------------------------------------------------
!     dynamically allocate the global BIOTIC arrays
!-----------------------------------------------------------------------

call allocate_arrays(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Domain%isd, Domain%ied, Domain%jsd, Domain%jed, grid%nk) 

      
!-----------------------------------------------------------------------
!       Open up the files for boundary conditions
!-----------------------------------------------------------------------

atmpress_id = init_external_field(atmpress_file,                &
                                  atmpress_name,                &
                                  domain = Domain%domain2d)
if (atmpress_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open atmpress file: ' //                      &
       trim(atmpress_file))
endif  !}

if (gasx_from_file) then
   pistonveloc_id = init_external_field(pistonveloc_file,          &
                                     pistonveloc_name,          &
                                     domain = Domain%domain2d)
   if (pistonveloc_id .eq. 0) then  !{
     call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open pistonveloc file: ' //                   &
       trim(pistonveloc_file))
   endif  !}
endif

!RASF I think the ifdafs are redundant
if (id_adic .ne. 0 .and. .not. use_access_co2) then
 aco2_id = init_external_field(aco2_file,                    &
                                      aco2_name,              &
                                      domain = Domain%domain2d)
 if (aco2_id .eq. 0) then  !{
   call mpp_error(FATAL, trim(error_header) //                   &
        'Could not open aco2 file: ' //                        &
        trim(aco2_file))
 endif  !}
endif

if (ice_file4gasx) then
   seaicefract_id = init_external_field(seaicefract_file,          &
                                     seaicefract_name,          &
                                     domain = Domain%domain2d)
   if (seaicefract_id .eq. 0) then  !{
     call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open seaicefract file: ' //                   &
       trim(seaicefract_file))
   endif  !}
endif

dust_id = init_external_field(dust_file,          &
                                     dust_name,          &
                                     domain = Domain%domain2d)
if (dust_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open dust file: ' //                   &
       trim(dust_file))
endif  !}


alphabio_id = init_external_field("INPUT/bgc_param.nc",          &
        "alphabio", domain = Domain%domain2d)
       
if (alphabio_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open bgc_param.nc file' )
endif  !}

parbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "parbio", domain = Domain%domain2d)
kwbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "kwbio", domain = Domain%domain2d)
kcbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "kcbio", domain = Domain%domain2d)
abio_id = init_external_field("INPUT/bgc_param.nc",          &
        "abio", domain = Domain%domain2d)
bbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "bbio", domain = Domain%domain2d)
cbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "cbio", domain = Domain%domain2d)
k1bio_id = init_external_field("INPUT/bgc_param.nc",          &
        "k1bio", domain = Domain%domain2d)
muepbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "muepbio", domain = Domain%domain2d)
muepsbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "muepsbio", domain = Domain%domain2d)
gam1bio_id = init_external_field("INPUT/bgc_param.nc",          &
        "gam1bio", domain = Domain%domain2d)
gbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "gbio", domain = Domain%domain2d)
epsbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "epsbio", domain = Domain%domain2d)
muezbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "muezbio", domain = Domain%domain2d)
gam2bio_id = init_external_field("INPUT/bgc_param.nc",          &
        "gam2bio", domain = Domain%domain2d)
muedbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "muedbio", domain = Domain%domain2d)
muecaco3_id = init_external_field("INPUT/bgc_param.nc",          &
        "muecaco3", domain = Domain%domain2d)
muedbio_sed_id = init_external_field("INPUT/bgc_param.nc",          &
        "muedbio_sed", domain = Domain%domain2d)
muecaco3_sed_id = init_external_field("INPUT/bgc_param.nc",          &
        "muecaco3_sed", domain = Domain%domain2d)
wdetbio_id = init_external_field("INPUT/bgc_param.nc",          &
        "wdetbio", domain = Domain%domain2d)
wcaco3_id = init_external_field("INPUT/bgc_param.nc",          &
        "wcaco3", domain = Domain%domain2d)
nat_co2_id = init_external_field("INPUT/bgc_param.nc",          &
        "nat_co2", domain = Domain%domain2d)
tscav_fe_id = init_external_field("INPUT/bgc_param.nc",          &
        "tscav_fe", domain = Domain%domain2d)
fe_bkgnd_id = init_external_field("INPUT/bgc_param.nc",          &
        "fe_bkgnd", domain = Domain%domain2d)
f_inorg_id = init_external_field("INPUT/bgc_param.nc",          &
        "f_inorg", domain = Domain%domain2d)

! ---------------------------------
!
! Read extra restart field(s)
!
! ---------------------------------

do n = 1, instances !{
 id_restart(1) = register_restart_field(sed_restart, "csiro_bgc_sediment.res.nc", "caco3_sediment", biotic(n)%caco3_sediment(:,:), domain=Domain%domain2d)
 id_restart(2) = register_restart_field(sed_restart, "csiro_bgc_sediment.res.nc", "det_sediment", biotic(n)%det_sediment(:,:), domain=Domain%domain2d)
enddo !} n

call restore_state(sed_restart)


!-----------------------------------------------------------------------
!     do some calendar calculation for the next section,
!      calculate the number of days in this year as well as
!      those that have already passed
!-----------------------------------------------------------------------

call get_date(time%model_time,                                  &
              year, month, day, hour, minute, second)
num_days = days_in_year(time%model_time)
days_in_this_year = 0.0
do m = 1, month - 1
  days_in_this_year = days_in_this_year +                       &
                     days_in_month(set_date(year, m, 1))
enddo
days_in_this_year = days_in_this_year + day - 1 + hour/24.0 +   &
                   minute/1440.0 + second/86400.0

!-----------------------------------------------------------------------
!       read in additional information for a restart
!-----------------------------------------------------------------------

write(stdout(),*)

!-----------------------------------------------------------------------
!
!       initialize some arrays which are held constant for this
!       simulation
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       register the fields
!-----------------------------------------------------------------------
id_pco2 = register_diag_field('ocean_model',                   &
     'pco2', grid%tracer_axes(1:2),                            &
     Time%model_time, 'pCO2', ' ',               &
     missing_value = -1.0e+10)
id_paco2 = register_diag_field('ocean_model',                   &
     'paco2', grid%tracer_axes(1:2),                            &
     Time%model_time, 'pCO2 inc. anthropogenic', ' ',               &
     missing_value = -1.0e+10)
id_co2_sat = register_diag_field('ocean_model',                  &
     'co2_saturation', grid%tracer_axes(1:2),                    &
     Time%model_time, 'CO2 saturation', 'mmol/m^3',            &
     missing_value = -1.0e+10)
id_aco2_sat = register_diag_field('ocean_model',                  &
     'aco2_saturation', grid%tracer_axes(1:2),                    &
     Time%model_time, 'CO2 saturation inc. anthropogenic', 'mmol/m^3',            &
     missing_value = -1.0e+10)

id_sc_o2 = register_diag_field('ocean_model',                   &
     'sc_o2', grid%tracer_axes(1:2),                            &
     Time%model_time, 'Schmidt number - O2', ' ',               &
     missing_value = -1.0e+10)

id_o2_sat = register_diag_field('ocean_model',                  &
     'o2_saturation', grid%tracer_axes(1:2),                    &
     Time%model_time, 'O2 saturation', 'mmolO2/m^3',            &
     missing_value = -1.0e+10)

id_kw_o2 = register_diag_field('ocean_model',                   &
     'kw_o2', grid%tracer_axes(1:2),                            &
     Time%model_time, 'Piston velocity - O2', ' ',              &
     missing_value = -1.0e+10)

id_light_limit = register_diag_field('ocean_model','light_limit', &
     grid%tracer_axes(1:2),Time%model_time, 'Integrated light limitation of phytoplankton growth', &
     ' ',missing_value = -1.0e+10)

id_adic_intmld = register_diag_field('ocean_model','adic_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated adic', &
     ' ',missing_value = -1.0e+10)
id_dic_intmld = register_diag_field('ocean_model','dic_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated dic', &
     ' ',missing_value = -1.0e+10)     
id_o2_intmld = register_diag_field('ocean_model','o2_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated o2', &
     ' ',missing_value = -1.0e+10)     
id_no3_intmld = register_diag_field('ocean_model','no3_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated no3', &
     ' ',missing_value = -1.0e+10)     
id_fe_intmld = register_diag_field('ocean_model','fe_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated fe', &
     ' ',missing_value = -1.0e+10)     
id_phy_intmld = register_diag_field('ocean_model','phy_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated phy', &
     ' ',missing_value = -1.0e+10)     
id_det_intmld = register_diag_field('ocean_model','det_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated det', &
     ' ',missing_value = -1.0e+10)     
id_pprod_gross_intmld = register_diag_field('ocean_model','pprod_gross_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated pprod_gross', &
     ' ',missing_value = -1.0e+10)     
id_npp_intmld = register_diag_field('ocean_model','npp_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated npp', &
     ' ',missing_value = -1.0e+10)     
id_radbio_intmld = register_diag_field('ocean_model','radbio_intmld', &
     grid%tracer_axes(1:2),Time%model_time, 'MLD-integrated radbio', &
     ' ',missing_value = -1.0e+10)     

id_wdet100avg = register_diag_field('ocean_model','wdet100avg', &
     grid%tracer_axes(1:2),Time%model_time, 'detritus sinking at 100 m', &
     ' ',missing_value = -1.0e+10)     

id_radbio1 = register_diag_field('ocean_model','radbio1', &
     grid%tracer_axes(1:2),Time%model_time, 'Photosynthetically active radiation for phytoplankton growth at surface', &
     'W m-2',missing_value = -1.0e+10)

id_radbio3d = register_diag_field('ocean_model','radbio3d', &
     grid%tracer_axes(1:3),Time%model_time, 'Photosynthetically active radiation for phytoplankton growth', &
     'W m-2',missing_value = -1.0e+10)

id_wdet100 = register_diag_field('ocean_model','wdet100', &
     grid%tracer_axes(1:2),Time%model_time, 'detritus export at 100 m (det*sinking rate)', &
     'mmolN/m^2/s',missing_value = -1.0e+10)

id_npp3d = register_diag_field('ocean_model','npp3d', &
     grid%tracer_axes(1:3),Time%model_time, 'Net primary productivity', &
     'mmolN/m^3/s',missing_value = -1.0e+10)

id_npp2d = register_diag_field('ocean_model','npp2d', &
     grid%tracer_axes(1:2),Time%model_time, 'Vertically integrated net primary productivity', &
     'mmolN/m^2/s',missing_value = -1.0e+10)

id_npp1 = register_diag_field('ocean_model','npp1', &
     grid%tracer_axes(1:2),Time%model_time, 'Net primary productivity in the first ocean layer', &
     'mmolN/m^2/s',missing_value = -1.0e+10)

id_pprod_gross = register_diag_field('ocean_model','pprod_gross', &
     grid%tracer_axes(1:3),Time%model_time, 'Gross PHY production', &
     'mmolN/m^3/s',missing_value = -1.0e+10)

id_pprod_gross_2d = register_diag_field('ocean_model','pprod_gross_2d', &
     grid%tracer_axes(1:2),Time%model_time, 'Vertically integrated Gross PHY production', &
     'mmolN/m^2/s',missing_value = -1.0e+10)

id_zprod_gross = register_diag_field('ocean_model','zprod_gross', &
     grid%tracer_axes(1:3),Time%model_time, 'Gross ZOO production', &
     'mmolN/m^3/s',missing_value = -1.0e+10)

id_caco3_sediment = register_diag_field('ocean_model','caco3_sediment', &
     grid%tracer_axes(1:2),Time%model_time, 'Accumulated CaCO3 in sediment at base of water column', &
     'mmolN/m^2',missing_value = -1.0e+10)

id_det_sediment = register_diag_field('ocean_model','det_sediment', &
     grid%tracer_axes(1:2),Time%model_time, 'Accumulated DET in sediment at base of water column', &
     'mmolN/m^2',missing_value = -1.0e+10)

id_caco3_sed_remin = register_diag_field('ocean_model','caco3_sed_remin', &
     grid%tracer_axes(1:2),Time%model_time, 'Rate of remineralisation of CaCO3 in accumulated sediment', &
     'mmolN/m^2',missing_value = -1.0e+10)

id_det_sed_remin = register_diag_field('ocean_model','det_sed_remin', &
     grid%tracer_axes(1:2),Time%model_time, 'Rate of remineralisation of DET in accumulated sediment', &
     'mmolN/m^2',missing_value = -1.0e+10)

id_caco3_sed_depst = register_diag_field('ocean_model','caco3_sed_depst', &
     grid%tracer_axes(1:2),Time%model_time, 'Rate of deposition of CaCO3 to sediment at base of water column', &
     'mmolN/m^2',missing_value = -1.0e+10)

id_det_sed_depst = register_diag_field('ocean_model','det_sed_depst', &
     grid%tracer_axes(1:2),Time%model_time, 'Rate of deposition of DET to sediment at base of water column', &
     'mmolN/m^2',missing_value = -1.0e+10)

id_total_co2_flux = register_diag_field('ocean_model','total_co2_flux', &
     Time%model_time, 'Total surface flux of inorganic C (natural) into ocean', &
     'Pg/yr',missing_value = -1.0e+30)

id_total_aco2_flux = register_diag_field('ocean_model','total_aco2_flux', &
     Time%model_time, 'Total surface flux of inorganic C (natural + anthropogenic) into ocean', &
     'Pg/yr',missing_value = -1.0e+30)

do n = 1, instances  !{

  if (instances .eq. 1) then  !{
    str = ' '
  else  !}{
    str = '_' // biotic(n)%name
  endif  !}

!  use generic definition for 
! surface tracer fluxes -stf1, stf2
! virtual tracer flux - vstf1
! sources - src1
 ntr_bgc=biotic(n)%ntr_bgc

 do nn=1,ntr_bgc
  if (nn .ge. 10) then
   write(bgc_stf(4:5),'(i2)') nn
   write(bgc_src(4:5),'(i2)') nn
   write(bgc_btf(4:5),'(i2)') nn
  else
    write(bgc_stf(4:4),'(i1)') 0
    write(bgc_src(4:4),'(i1)') 0
    write(bgc_btf(4:4),'(i1)') 0
    write(bgc_stf(5:5),'(i1)') nn
    write(bgc_src(5:5),'(i1)') nn
    write(bgc_btf(5:5),'(i1)') nn
  endif

! default field names...
  name1 = bgc_stf // ' flux into ocean'
  name2 = bgc_stf // ' virtual flux into ocean'
  name3 = bgc_src // ' source term'
  name4 = bgc_btf // ' flux into sediment'
! specified field names
  if (nn .eq. id_no3) then 
   name1 = 'Flux into ocean - nitrate'
   name2 = 'Virtual flux into ocean - nitrate'
   name3 = 'Source term - nitrate'
   name4 = 'Flux into sediment - nitrate'
  endif
  if (nn .eq. id_phy) then 
   name1 = 'Flux into ocean - phytoplankton'
   name2 = 'Virtual flux into ocean - phytoplankton'
   name3 = 'Source term - phytoplankton'
   name4 = 'Flux into sediment - phytoplankton'
  endif
  if (nn .eq. id_zoo) then 
   name1 = 'Flux into ocean - zooplankton'
   name2 = 'Virtual flux into ocean - zooplankton'
   name3 = 'Source term - zooplankton'
   name4 = 'Flux into sediment - zooplankton'
  endif
  if (nn .eq. id_det) then 
   name1 = 'Flux into ocean - detritus'
   name2 = 'Virtual flux into ocean - detritus'
   name3 = 'Source term - detritus'
   name4 = 'Flux into sediment - detritus'
  endif
  if (nn .eq. id_o2) then 
   name1 = 'Flux into ocean - oxygen'
   name2 = 'Virtual flux into ocean - oxygen'
   name3 = 'Source term - oxygen'
   name4 = 'Flux into sediment - oxygen'
  endif
  if (nn .eq. id_po4) then 
   name1 = 'Flux into ocean - phosphate'
   name2 = 'Virtual flux into ocean - phosphate'
   name3 = 'Source term - phosphate'
   name4 = 'Flux into sediment - phosphate'
  endif
  if (nn .eq. id_dic) then 
   name1 = 'Flux into ocean - DIC, PI'
   name2 = 'Virtual flux into ocean - DIC, PI'
   name3 = 'Source term - DIC, PI'
   name4 = 'Flux into sediment - DIC, PI'
  endif
  if (nn .eq. id_adic) then 
   name1 = 'Flux into ocean - DIC, inc. anth.'
   name2 = 'Virtual flux into ocean - DIC, inc. anth.'
   name3 = 'Source term - DIC, inc. anth.'
   name4 = 'Flux into sediment - DIC, inc. anth.'
  endif
  if (nn .eq. id_alk) then 
   name1 = 'Flux into ocean - alkalinity'
   name2 = 'Virtual flux into ocean - alkalinity'
   name3 = 'Source term - alkalinity'
   name4 = 'Flux into sediment - alkalinity'
  endif
  if (nn .eq. id_caco3) then 
   name1 = 'Flux into ocean - calcium carbonate'
   name2 = 'Virtual flux into ocean - calcium carbonate'
   name3 = 'Source term - calcium carbonate'
   name4 = 'Flux into sediment - calcium carbonate'
  endif
  if (nn .eq. id_fe) then
   name1 = 'Flux into ocean - iron'
   name2 = 'Virtual flux into ocean - iron'
   name3 = 'Source term - iron'
   name4 = 'Flux into sediment - iron'
  endif
  if (mpp_pe() == mpp_root_pe() )print*,'rjm bio',bgc_stf,'v'//bgc_stf

  biotic(n)%id_bgc_stf(nn) = register_diag_field('ocean_model',       &
       bgc_stf//str, grid%tracer_axes(1:2),                     &
       Time%model_time, name1, 'mmol/m^2/s',    &
!       Time%model_time, bgc_stf//'flux into ocean', 'mmol/m^2/s',    &
       missing_value = -1.0e+10)

  biotic(n)%id_bgc_vstf(nn) = register_diag_field('ocean_model',       &
       'v'//bgc_stf//str, grid%tracer_axes(1:2),                     &
       Time%model_time, name2, 'mmol/m^3/s',    &
!       Time%model_time, bgc_stf//'virtual flux into ocean', 'mmol/m^3/s',    &
       missing_value = -1.0e+10)

  biotic(n)%id_bgc_src(nn) = register_diag_field('ocean_model',       &
       bgc_src//str, grid%tracer_axes(1:3),                     &
       Time%model_time, name3, 'mmolN/m^3/s',    &
!       Time%model_time, bgc_src, 'mmolN/m^3/s',    &
       missing_value = -1.0e+10)

  biotic(n)%id_bgc_btf(nn) = register_diag_field('ocean_model',       &
       bgc_btf//str, grid%tracer_axes(1:2),                     &
       Time%model_time, name4, 'mmol/m^2/s',    &
!       Time%model_time, bgc_btf//'flux into sediment', 'mmol/m^2/s',    &
       missing_value = -1.0e+10)

 enddo  !} nn

enddo  !} n

!-----------------------------------------------------------------------
!     give info
!-----------------------------------------------------------------------

do n = 1, instances  !{
 biotic(n)%sal_global=34.6
 if (id_dic.ne.0) biotic(n)%bgc_global(id_dic)=1966.
 if (id_adic.ne.0) biotic(n)%bgc_global(id_adic)=1966.
 if (id_alk.ne.0) biotic(n)%bgc_global(id_alk)=2303.
 if (id_no3.ne.0) biotic(n)%bgc_global(id_no3)=4.67
enddo
  
write(stdout(),*)
write(stdout(),*) trim(note_header), 'Tracer runs initialized'
write(stdout(),*)


return
end subroutine  csiro_bgc_start  !}
! </SUBROUTINE> NAME="csiro_bgc_start"


!#######################################################################
! <SUBROUTINE NAME="csiro_bgc_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!

subroutine csiro_bgc_tracer (isc, iec, jsc, jec, t_prog, grid, time, dtts) !{

use mpp_mod, only : mpp_sum


type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
integer, intent(in)                             :: isc, iec
integer, intent(in)                             :: jsc, jec
real, intent(in)                                :: dtts

!-----------------------------------------------------------------------
!     local definitions
!-----------------------------------------------------------------------

character(len=64), parameter    :: sub_name = 'csiro_bgc_tracer'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

integer :: i
integer :: j
integer :: n
integer :: k
integer :: nn, ntr_bgc, ind_trc
real    :: total_trc
integer :: indsal
logical :: used


indsal = fm_get_index('/ocean_mod/prog_tracers/salt')

!-----------------------------------------------------------------------
!     accumulate global annual means
!-----------------------------------------------------------------------


do n = 1, instances  !{
!  use generic definition for bgc tracers -tracer1, tracer2
 ntr_bgc=biotic(n)%ntr_bgc

 ! Use zero_floor = .false. to skip this floor to BGC tracers and allow the time_tendency calculation to non-negative values bring any small negatives up.
 ! This then improves the 'mismatch' in diagnostic output.
 ! mac, aug11.
 if (zero_floor) then
  do nn=1,ntr_bgc


   do k = 1,grid%nk  !{
     do j = jsc, jec  !{
       do i = isc, iec  !{
        t_prog(biotic(n)%ind_bgc(nn))%field(i,j,k,time%taup1)= & 
        max(0.,t_prog(biotic(n)%ind_bgc(nn))%field(i,j,k,time%taup1) )
       enddo  !} i
     enddo  !} j
   enddo  !} k

  enddo !} nn 
 endif ! zero_floor

! comment out this sediment source of fe while testing equivalent code in csiro_bgc_bbc.  mac, nov12.  

! rjm bottom Fe fix
! mac aug10, only apply this fix when the water is <= 200 m deep.  
 if (id_fe.ne.0) then
    do j = jsc, jec  !{
      do i = isc, iec  !{
         if (grid%kmt(i,j) .gt. 0) then
            k = grid%kmt(i,j)
            if (grid%zw(k) .le. 200) &
               t_prog(biotic(n)%ind_bgc(id_fe))%field(i,j,k,time%taup1)= 0.999
         endif
      enddo  !} i
    enddo  !} j
 endif

! apply deposition and remineralisation rates to sediment tracers.  mac, nov12. 
 do j = jsc, jec  !{
   do i = isc, iec  !{
     if (grid%kmt(i,j) .gt. 0) then
       biotic(n)%det_sediment(i,j) = biotic(n)%det_sediment(i,j) + dtts* &
         ( biotic(n)%det_sed_depst(i,j) -   &
         biotic(n)%det_sed_remin(i,j) )
       biotic(n)%caco3_sediment(i,j) = biotic(n)%caco3_sediment(i,j) + dtts* &
         ( biotic(n)%caco3_sed_depst(i,j) -   &
         biotic(n)%caco3_sed_remin(i,j) )
     endif
   enddo  !} i
 enddo  !} j

 ! send sediment tracers to output. mac, nov12. 
 if (id_caco3_sediment .gt. 0) then
    used = send_data(id_caco3_sediment, biotic(n)%caco3_sediment(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
 endif
 if (id_det_sediment .gt. 0) then
    used = send_data(id_det_sediment, biotic(n)%det_sediment(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
 endif

enddo  !} n

return
end subroutine  csiro_bgc_tracer  !}
! </SUBROUTINE> NAME="csiro_bgc_tracer"

include "bio_v3.inc" 

end module  csiro_bgc_mod  !}
