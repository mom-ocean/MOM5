module vegetation_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: write_version_number, error_mesg, NOTE,FATAL, file_exist, close_file, &
                   check_nml_error, stdlog 
use mpp_mod, only: mpp_sum, mpp_max, mpp_pe, mpp_root_pe
use time_manager_mod, only: time_type, time_type_to_real, get_date, operator(-)
use constants_mod,    only: tfreeze, rdgas, rvgas, hlv, hlf, cp_air, PI
use sphum_mod, only: qscomp
use nf_utils_mod, only: nfu_def_var, nfu_get_var, nfu_put_var, nfu_inq_var

use vegn_tile_mod, only: vegn_tile_type, &
     vegn_seed_demand, vegn_seed_supply, vegn_add_bliving, &
     cpw, clw, csw
use soil_tile_mod, only: soil_tile_type, soil_ave_temp, &
                         soil_ave_theta0, soil_ave_theta1, soil_psi_stress
use land_constants_mod, only : NBANDS, BAND_VIS, d608, mol_C, mol_CO2, mol_air, &
     seconds_per_year
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=), &
     get_elmt_indices
use land_tile_diag_mod, only : &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, diag_buff_type
use land_data_mod,      only : land_state_type, lnd
use land_io_mod, only : read_field
use land_tile_io_mod, only : &
     create_tile_out_file, &
     write_tile_data_r0d_fptr, write_tile_data_i0d_fptr, write_tile_data_r1d_fptr, &
     read_tile_data_r0d_fptr,  read_tile_data_i0d_fptr,  read_tile_data_r1d_fptr, &
     print_netcdf_error, get_input_restart_name
use vegn_data_mod, only : SP_C4GRASS, LEAF_ON, LU_NTRL, read_vegn_data_namelist, &
     tau_drip_l, tau_drip_s, T_transp_min, cold_month_threshold, soil_carbon_depth_scale, &
     fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
     N_HARV_POOLS, HARV_POOL_NAMES
use vegn_cohort_mod, only : vegn_cohort_type, vegn_phys_prog_type, &
     update_species,&
     vegn_data_heat_capacity, vegn_data_intrcptn_cap, &
     get_vegn_wet_frac, vegn_data_cover
use canopy_air_mod, only : cana_turbulence
     
use cohort_io_mod, only :  read_create_cohorts, create_cohort_dimension, &
     read_cohort_data_r0d_fptr,  read_cohort_data_i0d_fptr,&
     write_cohort_data_r0d_fptr, write_cohort_data_i0d_fptr
use land_debug_mod, only : is_watch_point, set_current_point, check_temp_range
use vegn_radiation_mod, only : vegn_radiation_init, vegn_radiation
use vegn_photosynthesis_mod, only : vegn_photosynthesis_init, vegn_photosynthesis
use static_vegn_mod, only : read_static_vegn_namelist, static_vegn_init, static_vegn_end, &
     read_static_vegn
use vegn_dynamics_mod, only : vegn_dynamics_init, vegn_carbon_int, vegn_growth, &
     vegn_daily_npp, vegn_phenology, vegn_biogeography
use vegn_disturbance_mod, only : vegn_nat_mortality, vegn_disturbance, update_fuel
use vegn_harvesting_mod, only : &
     vegn_harvesting_init, vegn_harvesting_end, vegn_harvesting

implicit none
private

! ==== public interfaces =====================================================
public :: read_vegn_namelist
public :: vegn_init
public :: vegn_end
public :: save_vegn_restart

public :: vegn_get_cover
public :: vegn_radiation
public :: vegn_diffusion

public :: vegn_step_1
public :: vegn_step_2
public :: vegn_step_3

public :: update_vegn_slow
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegetation.F90,v 20.0 2013/12/13 23:30:58 fms Exp $', &
   tagname = '$Name: tikal $', &
   module_name = 'vegn'
! values for internal selector of CO2 option used for photosynthesis
integer, parameter :: VEGN_PHOT_CO2_PRESCRIBED  = 1
integer, parameter :: VEGN_PHOT_CO2_INTERACTIVE = 2


! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2               = .false.
real    :: init_Wl           = 0
real    :: init_Ws           = 0
real    :: init_Tv           = 288.
real    :: init_cohort_bl    = 0.05 ! initial biomass of leaves, kg C/m2
real    :: init_cohort_blv   = 0.0  ! initial biomass of labile store, kg C/m2
real    :: init_cohort_br    = 0.05 ! initial biomass of fine roots, kg C/m2
real    :: init_cohort_bsw   = 0.05 ! initial biomass of sapwood, kg C/m2
real    :: init_cohort_bwood = 0.05 ! initial biomass of heartwood, kg C/m2
real    :: init_cohort_cmc   = 0.0  ! initial intercepted water
character(32) :: rad_to_use = 'big-leaf' ! or 'two-stream'
character(32) :: snow_rad_to_use = 'ignore' ! or 'paint-leaves'
character(32) :: photosynthesis_to_use = 'simple' ! or 'leuning'
character(32) :: co2_to_use_for_photosynthesis = 'prescribed' ! or 'interactive'
   ! specifies what co2 concentration to use for photosynthesis calculations: 
   ! 'prescribed'  : a prescribed value is used, equal to co2_for_photosynthesis
   !      specified below.
   ! 'interactive' : concentration of co2 in canopy air is used
real    :: co2_for_photosynthesis = 350.0e-6 ! concentration of co2 for photosynthesis 
   ! calculations, mol/mol. Ignored if co2_to_use_for_photosynthesis is not 'prescribed'
character(32) :: soil_decomp_option = 'use_ave_t_and_theta' ! or 'use_layer_t_and_theta'
logical :: do_cohort_dynamics   = .TRUE. ! if true, do vegetation growth
logical :: do_patch_disturbance = .TRUE. ! 
logical :: do_phenology         = .TRUE. 
logical :: xwilt_available      = .TRUE.
logical :: do_biogeography      = .TRUE.
logical :: do_seed_transport    = .TRUE.
real    :: min_Wl=-1.0, min_Ws=-1.0 ! threshold values for condensation numerics, kg/m2:
   ! if water or snow on canopy fall below these values, the derivatives of
   ! condensation are set to zero, thereby prohibiting switching from condensation to
   ! evaporation in one time step.
real    :: tau_smooth_ncm = 0.0 ! Time scale for ncm smoothing
   ! (low-pass filtering), years. 0.0 retrieves previous behavior (no smoothing)
real :: rav_lit_0         = 0.0 ! constant litter resistance to vapor
real :: rav_lit_vi        = 0.0 ! litter resistance to vapor per LAI+SAI
real :: rav_lit_fsc       = 0.0 ! litter resistance to vapor per fsc
real :: rav_lit_ssc       = 0.0 ! litter resistance to vapor per ssc
real :: rav_lit_bwood     = 0.0 ! litter resistance to vapor per bwood
namelist /vegn_nml/ &
    lm2, init_Wl, init_Ws, init_Tv, cpw, clw, csw, &
    init_cohort_bl, init_cohort_blv, init_cohort_br, init_cohort_bsw, &
    init_cohort_bwood, init_cohort_cmc, &
    rad_to_use, snow_rad_to_use, photosynthesis_to_use, &
    co2_to_use_for_photosynthesis, co2_for_photosynthesis, &
    soil_decomp_option, &
    do_cohort_dynamics, do_patch_disturbance, do_phenology, &
    xwilt_available, &
    do_biogeography, do_seed_transport, &
    min_Wl, min_Ws, tau_smooth_ncm, &
    rav_lit_0, rav_lit_vi, rav_lit_fsc, rav_lit_ssc, rav_lit_bwood
    
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time ! *** NOT YET USED
real            :: delta_time      ! fast time step
real            :: dt_fast_yr      ! fast time step in years
integer         :: vegn_phot_co2_option = -1 ! internal selector of co2 option 
                                   ! used for photosynthesis
! diagnostic field ids
integer :: id_vegn_type, id_temp, id_wl, id_ws, id_height, id_lai, id_sai, id_leaf_size, &
   id_root_density, id_root_zeta, id_rs_min, id_leaf_refl, id_leaf_tran,&
   id_leaf_emis, id_snow_crit, id_stomatal, id_an_op, id_an_cl, &
   id_bl, id_blv, id_br, id_bsw, id_bwood, id_btot, id_species, id_status, &
   id_con_v_h, id_con_v_v, id_fuel, id_harv_pool(N_HARV_POOLS), &
   id_harv_rate(N_HARV_POOLS), id_t_harv_pool, id_t_harv_rate, &
   id_csmoke_pool, id_csmoke_rate, id_fsc_in, id_fsc_out, id_ssc_in, &
   id_ssc_out, id_veg_in, id_veg_out, id_fsc_pool, id_fsc_rate, &
   id_ssc_pool, id_ssc_rate, id_t_ann, id_t_cold, id_p_ann, id_ncm, &
   id_lambda, id_afire, id_atfall, id_closs, id_cgain, id_wdgain, id_leaf_age, &
   id_phot_co2, id_theph, id_psiph, id_evap_demand
! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_vegn_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  logical :: use_static_veg ! if true, switch off vegetation dynamics

  call read_vegn_data_namelist()
  call read_static_vegn_namelist(use_static_veg)

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=vegn_nml, iostat=io)
    ierr = check_nml_error(io, 'vegn_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=vegn_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'vegn_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()

  ! switch off vegetation dynamics if static vegetation is set
  if (use_static_veg) then
     call error_mesg('vegn_init', &
          'use_static_veg=.TRUE., switching off vegetation dynamics', NOTE)
     write(unit,*)'use_static_veg=.TRUE., switching off vegetation dynamics'
     do_cohort_dynamics   = .FALSE.
     do_patch_disturbance = .FALSE.
     do_phenology         = .FALSE. 
     do_biogeography      = .FALSE.
     do_seed_transport    = .FALSE.
  endif

  if (mpp_pe() == mpp_root_pe()) then
     write(unit, nml=vegn_nml)
  endif

  ! convert symbolic names of photosynthesis CO2 options into numeric IDs to
  ! speed up selection during run-time
  if (trim(co2_to_use_for_photosynthesis)=='prescribed') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_PRESCRIBED
  else if (trim(co2_to_use_for_photosynthesis)=='interactive') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_INTERACTIVE
  else
     call error_mesg('vegn_init',&
          'vegetation photosynthesis option co2_to_use_for_photosynthesis="'//&
          trim(co2_to_use_for_photosynthesis)//'" is invalid, use "prescribed" or "interactive"',&
          FATAL)
  endif

  ! ---- initialize vegetation radiation options
  call vegn_radiation_init(rad_to_use, snow_rad_to_use)

  ! ---- initialize vegetation photosynthesis options
  call vegn_photosynthesis_init(photosynthesis_to_use)

end subroutine read_vegn_namelist


! ============================================================================
! initialize vegetation
subroutine vegn_init ( id_lon, id_lat, id_band )
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in) :: id_band ! ID of spectral band axis

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! current and tail tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  type(vegn_cohort_type), pointer :: cohort! pointer to initial cohort for cold-start
  integer :: n_accum
  integer :: nmn_acm
  character(len=256) :: restart_file_name_1, restart_file_name_2
  logical :: restart_1_exists, restart_2_exists
  real, allocatable :: t_ann(:,:),t_cold(:,:),p_ann(:,:),ncm(:,:) ! buffers for biodata reading 
  logical :: did_read_biodata = .FALSE.
  integer :: i,j ! indices of current tile

  module_is_initialized = .TRUE.

  ! ---- make module copy of time and calculate time step ------------------
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year


  ! ---- initialize vegn state ---------------------------------------------
  n_accum = 0
  nmn_acm = 0
  call get_input_restart_name('INPUT/vegn1.res.nc',restart_1_exists, restart_file_name_1)
  call get_input_restart_name('INPUT/vegn2.res.nc',restart_2_exists, restart_file_name_2)
  if (restart_1_exists) then
     call error_mesg('vegn_init',&
          'reading NetCDF restarts "'//trim(restart_file_name_1)//&
                            '" and "'//trim(restart_file_name_2)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name_1,NF_NOWRITE,unit))
     ! read the cohort index and generate appropriate number of cohorts
     ! for each vegetation tile
     call read_create_cohorts(unit)
     
     ! read cohort data
     call read_cohort_data_r0d_fptr(unit, 'tv', cohort_tv_ptr )
     call read_cohort_data_r0d_fptr(unit, 'wl', cohort_wl_ptr )
     call read_cohort_data_r0d_fptr(unit, 'ws', cohort_ws_ptr )
     __NF_ASRT__(nf_close(unit))     

     __NF_ASRT__(nf_open(restart_file_name_2,NF_NOWRITE,unit))
     ! read global variables
     __NF_ASRT__(nfu_get_var(unit,'n_accum',n_accum))
     __NF_ASRT__(nfu_get_var(unit,'nmn_acm',nmn_acm))

     call read_cohort_data_i0d_fptr(unit, 'species', cohort_species_ptr )
     call read_cohort_data_r0d_fptr(unit, 'hite', cohort_height_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bl', cohort_bl_ptr )
     call read_cohort_data_r0d_fptr(unit, 'blv', cohort_blv_ptr )
     call read_cohort_data_r0d_fptr(unit, 'br', cohort_br_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bsw', cohort_bsw_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bwood', cohort_bwood_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bliving', cohort_bliving_ptr )
     call read_cohort_data_i0d_fptr(unit, 'status', cohort_status_ptr )
     if(nfu_inq_var(unit,'leaf_age')==NF_NOERR) &
          call read_cohort_data_r0d_fptr(unit,'leaf_age',cohort_leaf_age_ptr)
     call read_cohort_data_r0d_fptr(unit, 'npp_prev_day', cohort_npp_previous_day_ptr )

     if(nfu_inq_var(unit,'landuse')==NF_NOERR) &
          call read_tile_data_i0d_fptr(unit,'landuse',vegn_landuse_ptr)
     call read_tile_data_r0d_fptr(unit,'age',vegn_age_ptr)
     call read_tile_data_r0d_fptr(unit,'fsc_pool',vegn_fsc_pool_ptr)
     call read_tile_data_r0d_fptr(unit,'fsc_rate',vegn_fsc_rate_ptr)
     call read_tile_data_r0d_fptr(unit,'ssc_pool',vegn_ssc_pool_ptr)
     call read_tile_data_r0d_fptr(unit,'ssc_rate',vegn_ssc_rate_ptr)
     ! monthly-mean values
     call read_tile_data_r0d_fptr(unit,'tc_av', vegn_tc_av_ptr)
     if(nfu_inq_var(unit,'theta_av_phen')==NF_NOERR) then
        call read_tile_data_r0d_fptr(unit,'theta_av_phen', vegn_theta_av_phen_ptr)
        call read_tile_data_r0d_fptr(unit,'theta_av_fire', vegn_theta_av_fire_ptr)
        call read_tile_data_r0d_fptr(unit,'psist_av', vegn_psist_av_ptr)
     else
        call read_tile_data_r0d_fptr(unit,'theta_av', vegn_theta_av_phen_ptr)
        call read_tile_data_r0d_fptr(unit,'theta_av', vegn_theta_av_fire_ptr)
	! psist_av remains at initial value (equal to 0)
     endif
     call read_tile_data_r0d_fptr(unit,'tsoil_av', vegn_tsoil_av_ptr)
     call read_tile_data_r0d_fptr(unit,'precip_av', vegn_precip_av_ptr)
     call read_tile_data_r0d_fptr(unit,'lambda', vegn_lambda_ptr)
     call read_tile_data_r0d_fptr(unit,'fuel', vegn_fuel_ptr)
     ! annual-mean values
     call read_tile_data_r0d_fptr(unit,'t_ann', vegn_t_ann_ptr)
     call read_tile_data_r0d_fptr(unit,'t_cold', vegn_t_cold_ptr)
     call read_tile_data_r0d_fptr(unit,'p_ann', vegn_p_ann_ptr)
     call read_tile_data_r0d_fptr(unit,'ncm', vegn_ncm_ptr)
     ! accumulated values for annual averaging
     call read_tile_data_r0d_fptr(unit,'t_ann_acm', vegn_t_ann_acm_ptr)
     call read_tile_data_r0d_fptr(unit,'t_cold_acm', vegn_t_cold_acm_ptr)
     call read_tile_data_r0d_fptr(unit,'p_ann_acm', vegn_p_ann_acm_ptr)
     call read_tile_data_r0d_fptr(unit,'ncm_acm', vegn_ncm_acm_ptr)
     ! burned carbon pool and rate
     if(nfu_inq_var(unit,'csmoke_pool')==NF_NOERR) &
          call read_tile_data_r0d_fptr(unit,'csmoke_pool',vegn_csmoke_pool_ptr)
     if(nfu_inq_var(unit,'csmoke_rate')==NF_NOERR) &
          call read_tile_data_r0d_fptr(unit,'csmoke_rate',vegn_csmoke_rate_ptr)
     ! harvesting pools and rates
     do i = 1, N_HARV_POOLS
        if (nfu_inq_var(unit,trim(HARV_POOL_NAMES(i))//'_harv_pool')==NF_NOERR) &
             call read_tile_data_r1d_fptr(unit,trim(HARV_POOL_NAMES(i))//'_harv_pool',vegn_harv_pool_ptr,i)
        if (nfu_inq_var(unit,trim(HARV_POOL_NAMES(i))//'_harv_rate')==NF_NOERR) &
             call read_tile_data_r1d_fptr(unit,trim(HARV_POOL_NAMES(i))//'_harv_rate',vegn_harv_rate_ptr,i)
     enddo

     __NF_ASRT__(nf_close(unit))
  else
     call error_mesg('vegn_init',&
          'cold-starting vegetation',&
          NOTE)
  endif
  ! read climatological fields for initialization of species distribution
  if (file_exist('INPUT/biodata.nc'))then
     allocate(&
          t_ann (lnd%is:lnd%ie,lnd%js:lnd%je),&
          t_cold(lnd%is:lnd%ie,lnd%js:lnd%je),&
          p_ann (lnd%is:lnd%ie,lnd%js:lnd%je),&
          ncm   (lnd%is:lnd%ie,lnd%js:lnd%je) )
     call read_field( 'INPUT/biodata.nc','T_ANN', &
          lnd%lon, lnd%lat, t_ann, interp='nearest')
     call read_field( 'INPUT/biodata.nc','T_COLD', &
          lnd%lon, lnd%lat, t_cold, interp='nearest')
     call read_field( 'INPUT/biodata.nc','P_ANN', &
          lnd%lon, lnd%lat, p_ann, interp='nearest')
     call read_field( 'INPUT/biodata.nc','NCM', &
          lnd%lon, lnd%lat, ncm, interp='nearest')
     did_read_biodata = .TRUE.
  endif
  ! Go through all tiles and initialize the cohorts that have not been initialized yet --
  ! this allows to read partial restarts. Also initialize accumulation counters to zero
  ! or the values from the restarts.
  te = tail_elmt(lnd%tile_map)
  ce = first_elmt(lnd%tile_map, is=lnd%is, js=lnd%js)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     call get_elmt_indices(ce,i,j)
     ce=next_elmt(ce)       ! advance position to the next tile
     if (.not.associated(tile%vegn)) cycle

     tile%vegn%n_accum = n_accum
     tile%vegn%nmn_acm = nmn_acm

     if (tile%vegn%n_cohorts>0) cycle ! skip initialized tiles
     
     ! create and initialize cohorts for this vegetation tile
     ! for now, just create a new cohort with default values of biomasses
     tile%vegn%n_cohorts = 1
     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
     cohort => tile%vegn%cohorts(1)
     cohort%prog%Wl = init_Wl
     cohort%prog%Ws = init_Ws
     cohort%prog%Tv = init_Tv
     
     cohort%bl      = init_cohort_bl
     cohort%blv     = init_cohort_blv
     cohort%br      = init_cohort_br
     cohort%bsw     = init_cohort_bsw
     cohort%bwood   = init_cohort_bwood
     cohort%bliving = cohort%bl+cohort%br+cohort%blv+cohort%bsw
     cohort%npp_previous_day = 0.0
     cohort%status  = LEAF_ON
     cohort%leaf_age = 0.0
     if(did_read_biodata.and.do_biogeography) then
        call update_species(cohort,t_ann(i,j),t_cold(i,j),p_ann(i,j),ncm(i,j),LU_NTRL)
     else
        cohort%species = tile%vegn%tag
     endif
  enddo
    
  ! initialize carbon integrator
  call vegn_dynamics_init ( id_lon, id_lat, lnd%time, delta_time, soil_decomp_option )

  ! initialize static vegetation
  call static_vegn_init ()
  call read_static_vegn ( lnd%time )

  ! initialize harvesting options
  call vegn_harvesting_init()

  ! initialize vegetation diagnostic fields
  call vegn_diag_init ( id_lon, id_lat, id_band, lnd%time )

  ! ---- diagnostic section
  ce = first_elmt(lnd%tile_map, is=lnd%is, js=lnd%js)
  te  = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     tile => current_tile(ce)
     ce=next_elmt(ce)     
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetation tiles
     ! send the data
     call send_tile_data(id_vegn_type,  real(tile%vegn%tag), tile%diag)
  enddo

  if (allocated(t_ann))  deallocate(t_ann)
  if (allocated(t_cold)) deallocate(t_cold)
  if (allocated(p_ann))  deallocate(p_ann)
  if (allocated(ncm))    deallocate(ncm)

end subroutine vegn_init

! ============================================================================
subroutine vegn_diag_init ( id_lon, id_lat, id_band, time )
  integer        , intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer        , intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer        , intent(in) :: id_band ! ID of spectral band axis
  type(time_type), intent(in) :: time    ! initial time for diagnostic fields
  
  ! ---- local vars
  integer :: i

  id_vegn_type = register_tiled_static_field ( module_name, 'vegn_type',  &
       (/id_lon,id_lat/), 'vegetation type', missing_value=-1.0 )

  id_temp = register_tiled_diag_field ( module_name, 'temp',  &
       (/id_lon,id_lat/), time, 'canopy temperature', 'degK', missing_value=-1.0 )
  id_wl = register_tiled_diag_field ( module_name, 'wl',  &
       (/id_lon,id_lat/), time, 'canopy liquid water content', 'kg/m2', missing_value=-1.0 )
  id_ws = register_tiled_diag_field ( module_name, 'ws',  &
       (/id_lon,id_lat/), time, 'canopy solid water content', 'kg/m2', missing_value=-1.0 )

  id_height = register_tiled_diag_field ( module_name, 'height',  &
       (/id_lon,id_lat/), time, 'vegetation height', 'm', missing_value=-1.0 )
  id_lai    = register_tiled_diag_field ( module_name, 'lai',  &
       (/id_lon,id_lat/), time, 'leaf area index', 'm2/m2', missing_value=-1.0 )
  id_sai    = register_tiled_diag_field ( module_name, 'sai',  &
       (/id_lon,id_lat/), time, 'stem area index', 'm2/m2', missing_value=-1.0 )
  id_leaf_size = register_tiled_diag_field ( module_name, 'leaf_size',  &
       (/id_lon,id_lat/), time, missing_value=-1.0 )
  id_root_density = register_tiled_diag_field ( module_name, 'root_density',  &
       (/id_lon,id_lat/), time, 'total biomass below ground', 'kg/m2', missing_value=-1.0 )
  id_root_zeta = register_tiled_diag_field ( module_name, 'root_zeta',  &
       (/id_lon,id_lat/), time, 'e-folding depth of root biomass', 'm',missing_value=-1.0 )
  id_rs_min = register_tiled_diag_field ( module_name, 'rs_min',  &
       (/id_lon,id_lat/), time, missing_value=-1.0 )
  id_leaf_refl = register_tiled_diag_field ( module_name, 'leaf_refl',  &
       (/id_lon,id_lat,id_band/), time, 'reflectivity of leaf', missing_value=-1.0 )
  id_leaf_tran = register_tiled_diag_field ( module_name, 'leaf_tran',  &
       (/id_lon,id_lat,id_band/), time, 'transmittance of leaf', missing_value=-1.0 )
  id_leaf_emis = register_tiled_diag_field ( module_name, 'leaf_emis',  &
       (/id_lon,id_lat/), time, 'leaf emissivity', missing_value=-1.0 )
  id_snow_crit = register_tiled_diag_field ( module_name, 'snow_crit',  &
       (/id_lon,id_lat/), time, missing_value=-1.0 )
  id_stomatal = register_tiled_diag_field ( module_name, 'stomatal_cond',  &
       (/id_lon,id_lat/), time, 'vegetation stomatal conductance', missing_value=-1.0 )
  id_evap_demand = register_tiled_diag_field ( module_name, 'evap_demand',  &
       (/id_lon,id_lat/), time, 'plant evaporative water demand',&
       'kg/(m2 s)', missing_value=-1e20 )
  id_an_op = register_tiled_diag_field ( module_name, 'an_op',  &
       (/id_lon,id_lat/), time, 'net photosynthesis with open stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20 )
  id_an_cl = register_tiled_diag_field ( module_name, 'an_cl',  &
       (/id_lon,id_lat/), time, 'net photosynthesis with closed stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20 )

  id_bl = register_tiled_diag_field ( module_name, 'bl',  &
       (/id_lon,id_lat/), time, 'biomass of leaves', 'kg C/m2', missing_value=-1.0 )
  id_blv = register_tiled_diag_field ( module_name, 'blv',  &
       (/id_lon,id_lat/), time, 'biomass in labile store', 'kg C/m2', missing_value=-1.0 )
  id_br = register_tiled_diag_field ( module_name, 'br',  &
       (/id_lon,id_lat/), time, 'biomass of fine roots', 'kg C/m2', missing_value=-1.0 )
  id_bsw = register_tiled_diag_field ( module_name, 'bsw',  &
       (/id_lon,id_lat/), time, 'biomass of sapwood', 'kg C/m2', missing_value=-1.0 )
  id_bwood = register_tiled_diag_field ( module_name, 'bwood',  &
       (/id_lon,id_lat/), time, 'biomass of heartwood', 'kg C/m2', missing_value=-1.0 )
  id_btot = register_tiled_diag_field ( module_name, 'btot',  &
       (/id_lon,id_lat/), time, 'total biomass', 'kg C/m2', missing_value=-1.0 )
  id_fuel = register_tiled_diag_field ( module_name, 'fuel',  &
       (/id_lon,id_lat/), time, 'mass of fuel', 'kg C/m2', missing_value=-1.0 )
  id_lambda = register_tiled_diag_field (module_name, 'lambda',(/id_lon,id_lat/), &
       time, 'drought', 'months', missing_value=-100.0)

  id_species = register_tiled_diag_field ( module_name, 'species',  &
       (/id_lon,id_lat/), time, 'vegetation species number', missing_value=-1.0 )
  id_status = register_tiled_diag_field ( module_name, 'status',  &
       (/id_lon,id_lat/), time, 'status of leaves', missing_value=-1.0 )
  id_theph = register_tiled_diag_field ( module_name, 'theph',  &
       (/id_lon,id_lat/), time, 'theta for phenology', missing_value=-1.0 )
  id_psiph = register_tiled_diag_field ( module_name, 'psiph',  &
       (/id_lon,id_lat/), time, 'psi stress for phenology', missing_value=-1.0 )
  id_leaf_age = register_tiled_diag_field ( module_name, 'leaf_age',  &
       (/id_lon,id_lat/), time, 'age of leaves since bud burst', 'days', missing_value=-1.0 )!ens

  id_con_v_h = register_tiled_diag_field ( module_name, 'con_v_h', (/id_lon,id_lat/), &
       time, 'conductance for sensible heat between canopy and canopy air', &
       'm/s', missing_value=-1.0 )
  id_con_v_v = register_tiled_diag_field ( module_name, 'con_v_v', (/id_lon,id_lat/), &
       time, 'conductance for water vapor between canopy and canopy air', &
       'm/s', missing_value=-1.0 )

  id_cgain = register_tiled_diag_field ( module_name, 'cgain', (/id_lon,id_lat/), &
       time, 'carbon gain', 'kg C', missing_value=-100.0 )
  id_closs = register_tiled_diag_field ( module_name, 'closs', (/id_lon,id_lat/), &
       time, 'carbon loss', 'kg C', missing_value=-100.0 )
  id_wdgain = register_tiled_diag_field ( module_name, 'wdgain', (/id_lon,id_lat/), &
       time, 'wood biomass gain', 'kg C', missing_value=-100.0 )

  id_t_ann  = register_tiled_diag_field ( module_name, 't_ann', (/id_lon,id_lat/), &
       time, 'annual mean temperature', 'degK', missing_value=-999.0 )
  id_t_cold  = register_tiled_diag_field ( module_name, 't_cold', (/id_lon,id_lat/), &
       time, 'average temperature of the coldest month', 'degK', missing_value=-999.0 )
  id_p_ann  = register_tiled_diag_field ( module_name, 'p_ann', (/id_lon,id_lat/), &
       time, 'annual mean precipitation', 'kg/(m2 s)', missing_value=-999.0 )
  id_ncm = register_tiled_diag_field ( module_name, 'ncm', (/id_lon,id_lat/), &
       time, 'number of cold months', 'dimensionless', missing_value=-999.0 )

  id_t_harv_pool = register_tiled_diag_field( module_name, 'harv_pool', (/id_lon,id_lat/), &
       time, 'total harvested carbon', 'kg C/m2', missing_value=-999.0)
  id_t_harv_rate = register_tiled_diag_field( module_name, 'harv_rate', (/id_lon,id_lat/), &
       time, 'total rate of release of harvested carbon to the atmosphere', &
       'kg C/(m2 year)', missing_value=-999.0)
  do i = 1,N_HARV_POOLS
     id_harv_pool(i) = register_tiled_diag_field( module_name, &
          trim(HARV_POOL_NAMES(i))//'_harv_pool', (/id_lon,id_lat/), time, &
          'harvested carbon', 'kg C/m2', missing_value=-999.0)
     id_harv_rate(i) = register_tiled_diag_field( module_name, &
          trim(HARV_POOL_NAMES(i))//'_harv_rate', (/id_lon,id_lat/), time, &
          'rate of release of harvested carbon to the atmosphere', 'kg C/(m2 year)', &
          missing_value=-999.0)
  enddo

  id_fsc_pool = register_tiled_diag_field (module_name, 'fsc_pool', (/id_lon, id_lat/), &
       time, 'intermediate pool of fast soil carbon', 'kg C/m2', missing_value=-999.0)
  id_fsc_rate = register_tiled_diag_field (module_name, 'fsc_rate', (/id_lon, id_lat/), &
       time, 'rate of conversion of fsc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_ssc_pool = register_tiled_diag_field (module_name, 'ssc_pool', (/id_lon, id_lat/), &
       time, 'intermediate pool of slow soil carbon', 'kg C/m2', missing_value=-999.0)
  id_ssc_rate = register_tiled_diag_field (module_name, 'ssc_rate', (/id_lon, id_lat/), &
       time, 'rate of conversion of ssc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)

  id_csmoke_pool = register_tiled_diag_field ( module_name, 'csmoke', (/id_lon, id_lat/), &
       time, 'carbon lost through fire', 'kg C/m2', missing_value=-999.0)
  id_csmoke_rate = register_tiled_diag_field ( module_name, 'csmoke_rate', (/id_lon, id_lat/), &
       time, 'rate of release of carbon lost through fire to the atmosphere', &
       'kg C/(m2 yr)', missing_value=-999.0)

  id_ssc_in = register_tiled_diag_field ( module_name, 'ssc_in',  (/id_lon, id_lat/), &
     time,  'soil slow carbon in', 'kg C/m2', missing_value=-999.0 )
  id_ssc_out = register_tiled_diag_field ( module_name, 'ssc_out',  (/id_lon, id_lat/), &
     time,  'soil slow carbon out', 'kg C/m2', missing_value=-999.0 )
  id_fsc_in = register_tiled_diag_field ( module_name, 'fsc_in',  (/id_lon, id_lat/), &
     time,  'soil fast carbon in', 'kg C/m2', missing_value=-999.0 )
  id_fsc_out = register_tiled_diag_field ( module_name, 'fsc_out',  (/id_lon, id_lat/), &
     time,  'soil fast carbon out', 'kg C/m2', missing_value=-999.0 )
  id_veg_in = register_tiled_diag_field ( module_name, 'veg_in',  (/id_lon, id_lat/), &
     time,  'vegetation carbon in', 'kg C/m2', missing_value=-999.0 )
  id_veg_out = register_tiled_diag_field ( module_name, 'veg_out',  (/id_lon, id_lat/), &
     time,  'vegetation carbon out', 'kg C/m2', missing_value=-999.0 )

  id_afire = register_tiled_diag_field (module_name, 'afire', (/id_lon,id_lat/), &
       time, 'area been fired', missing_value=-100.0)
  id_atfall = register_tiled_diag_field (module_name, 'atfall',(/id_lon,id_lat/), &
       time, 'area been disturbed', missing_value=-100.0)

  id_phot_co2 = register_tiled_diag_field (module_name, 'qco2_phot',(/id_lon,id_lat/), &
       time, 'CO2 mixing ratio for photosynthesis calculations', 'mol CO2/mol dry air', &
       missing_value=-1.0)
end subroutine


! ============================================================================
! write restart file and release memory
subroutine vegn_end ()

  module_is_initialized =.FALSE.

  ! finalize harvesting
  call vegn_harvesting_end ()

  ! finalize static vegetation, if necessary
  call static_vegn_end ()
end subroutine vegn_end


! ============================================================================
subroutine save_vegn_restart(tile_dim_length,timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit ! restart file unit 
  integer :: ierr, i
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: n_accum, nmn_acm

  call error_mesg('vegn_end','writing NetCDF restart',NOTE)
  ! create output file, including internal structure necessary for tile output
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'vegn1.res.nc', &
          lnd%coord_glon, lnd%coord_glat, vegn_tile_exists, tile_dim_length)
  ! create compressed dimension for vegetation cohorts -- must be called even
  ! if restart has not been created, because it calls mpp_max and that should 
  ! be called on all PEs to work
  call create_cohort_dimension(unit)

  call write_cohort_data_r0d_fptr(unit,'tv',cohort_tv_ptr,'vegetation temperature','degrees_K')
  call write_cohort_data_r0d_fptr(unit,'wl',cohort_wl_ptr,'vegetation liquid water content','kg/m2')
  call write_cohort_data_r0d_fptr(unit,'ws',cohort_ws_ptr,'vegetation solid water content','kg/m2')
  ! close output file
  __NF_ASRT__(nf_close(unit))


  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'vegn2.res.nc', &
          lnd%coord_glon, lnd%coord_glat, vegn_tile_exists, tile_dim_length )
  ! create compressed dimension for vegetation cohorts -- see note above
  call create_cohort_dimension(unit)

  ! store global variables
  ! find first tile and get n_accum and nmn_acm from it
  n_accum = 0; nmn_acm = 0
  ce = first_elmt(lnd%tile_map) ; te = tail_elmt(lnd%tile_map)
  do while ( ce /= te )
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(associated(tile%vegn)) then
        n_accum = tile%vegn%n_accum
        nmn_acm = tile%vegn%nmn_acm
     endif
  enddo
  ! n_accum and nmn_acm are currently the same for all tiles; we only call mpp_max
  ! to handle the situation when there are no tiles in the current domain
  call mpp_max(n_accum); call mpp_max(nmn_acm)
  
  if(mpp_pe()==lnd%io_pelist(1)) then
     ierr = nf_redef(unit)
     __NF_ASRT__(nfu_def_var(unit,'n_accum',NF_INT,long_name='number of accumulated steps'))
     __NF_ASRT__(nfu_def_var(unit,'nmn_acm',NF_INT,long_name='number of accumulated months'))
     ierr = nf_enddef(unit)
     __NF_ASRT__(nfu_put_var(unit,'n_accum',n_accum))
     __NF_ASRT__(nfu_put_var(unit,'nmn_acm',nmn_acm))
  end if
  
  call write_cohort_data_i0d_fptr(unit,'species', cohort_species_ptr, 'vegetation species')
  call write_cohort_data_r0d_fptr(unit,'hite', cohort_height_ptr, 'vegetation height','m')
  call write_cohort_data_r0d_fptr(unit,'bl', cohort_bl_ptr, 'biomass of leaves per individual','kg C/m2')
  call write_cohort_data_r0d_fptr(unit,'blv', cohort_blv_ptr, 'biomass of virtual leaves (labile store) per individual','kg C/m2')
  call write_cohort_data_r0d_fptr(unit,'br', cohort_br_ptr, 'biomass of fine roots per individual','kg C/m2')
  call write_cohort_data_r0d_fptr(unit,'bsw', cohort_bsw_ptr, 'biomass of sapwood per individual','kg C/m2')
  call write_cohort_data_r0d_fptr(unit,'bwood', cohort_bwood_ptr, 'biomass of heartwood per individual','kg C/m2')
  call write_cohort_data_r0d_fptr(unit,'bliving', cohort_bliving_ptr, 'total living biomass per individual','kg C/m2')
!     call write_cohort_data_r0d_fptr(unit,'tleaf', cohort_tleaf_ptr, 'leaf temperature','degK')
  call write_cohort_data_i0d_fptr(unit,'status', cohort_status_ptr, 'leaf status')
  call write_cohort_data_r0d_fptr(unit,'leaf_age',cohort_leaf_age_ptr, 'age of leaves since bud burst', 'days')

!     call write_cohort_data_r0d_fptr(unit,'intercept_l', cohort_cmc_ptr, 'intercepted water per cohort','kg/m2')
  call write_cohort_data_r0d_fptr(unit,'npp_prev_day', cohort_npp_previous_day_ptr, 'previous day NPP','kg C/(m2 year)')

  call write_tile_data_i0d_fptr(unit,'landuse',vegn_landuse_ptr,'vegetation land use type')
  call write_tile_data_r0d_fptr(unit,'age',vegn_age_ptr,'vegetation age', 'yr')
  call write_tile_data_r0d_fptr(unit,'fsc_pool',vegn_fsc_pool_ptr,'intermediate pool for fast soil carbon input', 'kg C/m2')
  call write_tile_data_r0d_fptr(unit,'fsc_rate',vegn_fsc_rate_ptr,'conversion rate of fsc_pool to fast soil carbon', 'kg C/(m2 yr)')
  call write_tile_data_r0d_fptr(unit,'ssc_pool',vegn_ssc_pool_ptr,'intermediate pool for slow soil carbon input', 'kg C/m2')
  call write_tile_data_r0d_fptr(unit,'ssc_rate',vegn_ssc_rate_ptr,'conversion rate of ssc_pool to slow soil carbon', 'kg C/(m2 yr)')

  ! monthly-mean values
  call write_tile_data_r0d_fptr(unit,'tc_av', vegn_tc_av_ptr,'average canopy air temperature','degK')
  call write_tile_data_r0d_fptr(unit,'theta_av_phen', vegn_theta_av_phen_ptr,'average soil moisture for phenology')
  call write_tile_data_r0d_fptr(unit,'theta_av_fire', vegn_theta_av_fire_ptr,'average soil moisture for fire')
  call write_tile_data_r0d_fptr(unit,'psist_av', vegn_psist_av_ptr,'average soil-water-stress index')
  call write_tile_data_r0d_fptr(unit,'tsoil_av', vegn_tsoil_av_ptr,'average bulk soil temperature for soil carbon','degK')
  call write_tile_data_r0d_fptr(unit,'precip_av', vegn_precip_av_ptr,'average total precipitation','kg/(m2 s)')
  call write_tile_data_r0d_fptr(unit,'lambda', vegn_lambda_ptr,'dryness parameter')
  call write_tile_data_r0d_fptr(unit,'fuel', vegn_fuel_ptr,'fuel density','kg C/m2')
  ! annual-mean values
  call write_tile_data_r0d_fptr(unit,'t_ann', vegn_t_ann_ptr,'average annual canopy air temperature','degK')
  call write_tile_data_r0d_fptr(unit,'t_cold', vegn_t_cold_ptr,'average canopy air temperature of coldest month','degK')
  call write_tile_data_r0d_fptr(unit,'p_ann', vegn_p_ann_ptr,'average annual precipitation','kg/(m2 s)')
  call write_tile_data_r0d_fptr(unit,'ncm', vegn_ncm_ptr,'number of cold months')
  ! accumulated values for annual averaging
  call write_tile_data_r0d_fptr(unit,'t_ann_acm', vegn_t_ann_acm_ptr,'accumulated annual canopy air temperature','degK')
  call write_tile_data_r0d_fptr(unit,'t_cold_acm', vegn_t_cold_acm_ptr,'accumulated temperature of coldest month','degK')
  call write_tile_data_r0d_fptr(unit,'p_ann_acm', vegn_p_ann_acm_ptr,'accumulated precipitation','kg/(m2 s)')
  call write_tile_data_r0d_fptr(unit,'ncm_acm', vegn_ncm_acm_ptr,'accumulated number of cold months')

  ! burned carbon pool and rate
  call write_tile_data_r0d_fptr(unit,'csmoke_pool',vegn_csmoke_pool_ptr,'carbon lost through fires', 'kg C/m2')
  call write_tile_data_r0d_fptr(unit,'csmoke_rate',vegn_csmoke_rate_ptr,'rate of release of carbon lost through fires to the atmosphere', 'kg C/(m2 yr)')

  ! harvesting pools and rates
  do i = 1, N_HARV_POOLS
     call write_tile_data_r1d_fptr(unit, trim(HARV_POOL_NAMES(i))//'_harv_pool', &
          vegn_harv_pool_ptr, i, 'harvested carbon','kg C/m2')
     call write_tile_data_r1d_fptr(unit, trim(HARV_POOL_NAMES(i))//'_harv_rate', &
          vegn_harv_rate_ptr, i, 'rate of release of harvested carbon to the atmosphere','kg C/(m2 yr)')
  enddo
     

  __NF_ASRT__(nf_close(unit))

end subroutine save_vegn_restart


! ============================================================================
subroutine vegn_get_cover(vegn, snow_depth, vegn_cover)
  type(vegn_tile_type), intent(inout)  :: vegn ! it is only inout because vegn%data%cover
                                    ! changes cohort; can it be avoided?
  real,                 intent(in)  :: snow_depth
  real,                 intent(out) :: vegn_cover

  real :: vegn_cover_snow_factor

  call vegn_data_cover(vegn%cohorts(1), snow_depth, vegn_cover, vegn_cover_snow_factor)
  
end subroutine vegn_get_cover


! ============================================================================
subroutine vegn_diffusion ( vegn, vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf)
  type(vegn_tile_type), intent(in) :: vegn
  real,                intent(out) :: &
       vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf
  
  vegn_cover  = vegn%cohorts(1)%cover
  vegn_lai    = vegn%cohorts(1)%lai
  vegn_sai    = vegn%cohorts(1)%sai
  vegn_height = vegn%cohorts(1)%height
  vegn_d_leaf = vegn%cohorts(1)%leaf_size

end subroutine vegn_diffusion


! ============================================================================
subroutine vegn_step_1 ( vegn, soil, diag, &
        p_surf, ustar, drag_q, &
        SWdn, RSv, precip_l, precip_s, &
        land_d, land_z0s, land_z0m, grnd_z0s, &
        soil_beta, soil_water_supply, &
        cana_T, cana_q, cana_co2_mol, &
        ! output
        con_g_h, con_g_v, & ! aerodynamic conductance between canopy air and canopy, for heat and vapor flux
        vegn_T,vegn_Wl,  vegn_Ws,           & ! temperature, water and snow mass of the canopy
        vegn_ifrac,                         & ! intercepted fraction of liquid and frozen precipitation
        vegn_lai,                           & ! leaf area index
        drip_l, drip_s,                     & ! water and snow drip rate from precipitation, kg/(m2 s)
        vegn_hcap,                          & ! vegetation heat capacity
        Hv0,   DHvDTv,   DHvDTc,            & ! sens heat flux
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
        Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf  ) ! sublimation of intercepted snow
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(soil_tile_type), intent(inout) :: soil ! soil data
  ! TODO: possibly move calculation of soil-related stuff from calling subroutine to here
  !       now since we have soil tiled passed to us
  type(diag_buff_type), intent(inout) :: diag ! diagnostic buffer
  real, intent(in) :: &
       p_surf,    & ! surface pressure, N/m2
       ustar,     & ! friction velocity, m/s
       drag_q,    & ! bulk drag coefficient for specific humidity
       SWdn(NBANDS), & ! downward SW radiation at the top of the canopy, W/m2
       RSv (NBANDS), & ! net SW radiation balance of the canopy, W/m2
       precip_l, precip_s, & ! liquid and solid precipitation rates, kg/(m2 s)
       land_d, land_z0s, land_z0m, & ! land displacement height and roughness, m
       grnd_z0s, & ! roughness of ground surface (including snow effect)
       soil_beta, & ! relative water availability
       soil_water_supply, & ! max rate of water supply to the roots, kg/(m2 s)
       cana_T,    & ! temperature of canopy air, deg K
       cana_q,    & ! specific humidity of canopy air, kg/kg
       cana_co2_mol ! co2 mixing ratio in the canopy air, mol CO2/mol dry air
  ! output -- coefficients of linearized expressions for fluxes
  real, intent(out) ::   &
       vegn_T,vegn_Wl,  vegn_Ws,& ! temperature, water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid and frozen precipitation
       vegn_lai, & ! vegetation leaf area index
       drip_l, drip_s, & ! water and snow drip rate from precipitation, kg/(m2 s)
       vegn_hcap, & ! total vegetation heat capacity, including intercepted water and snow
       con_g_h, con_g_v, & ! aerodynamic conductance between ground and canopy air
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf    ! sublimation of intercepted snow
  
  ! ---- local vars 
  real :: &
       ft,DftDwl,DftDwf, & ! fraction of canopy not covered by intercepted water/snow, and its' 
                    ! derivatives w.r.t. intercepted water masses 
       fw,DfwDwl,DfwDwf, & ! fraction of canopy covered by intercepted water, and its' 
                    ! derivatives w.r.t. intercepted water masses 
       fs,DfsDwl,DfsDwf, & ! fraction of canopy covered by intercepted snow, and its' 
                    ! derivatives w.r.t. intercepted water masses
       stomatal_cond, & ! integral stomatal conductance of canopy
       con_v_h, con_v_v, & ! aerodyn. conductance between canopy and CAS, for heat and vapor
       rav_lit,   & ! additional resistance of litter to vapor transport
       total_cond, &! overall conductance from inside stomata to canopy air 
       qvsat,     & ! sat. specific humidity at the leaf T
       DqvsatDTv, & ! derivative of qvsat w.r.t. leaf T
       rho,       & ! density of canopy air
       phot_co2,  & ! co2 mixing ratio for photosynthesis, mol CO2/mol dry air
       evap_demand, & ! evaporative water demand, kg/(m2 s)
       photosynt, & ! photosynthesis
       photoresp    ! photo-respiration
  type(vegn_cohort_type), pointer :: cohort
  
  ! get the pointer to the first (and, currently, the only) cohort
  cohort => vegn%cohorts(1)

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_1 input ####'
     __DEBUG3__(p_surf, ustar, drag_q)
     __DEBUG1__(SWdn)
     __DEBUG1__(RSv) 
     __DEBUG2__(precip_l, precip_s)
     __DEBUG4__(land_d, land_z0s, land_z0m, grnd_z0s)
     __DEBUG2__(soil_beta, soil_water_supply)
     __DEBUG3__(cana_T, cana_q, cana_co2_mol)
     write(*,*)'#### end of vegn_step_1 input ####'
     __DEBUG3__(cohort%height, cohort%lai, cohort%sai)
     __DEBUG2__(cohort%cover,cohort%leaf_size)
     __DEBUG1__(cohort%prog%Tv)
  endif

  ! check the range of input temperature
  call check_temp_range(cohort%prog%Tv,'vegn_step_1','cohort%prog%Tv', lnd%time) 

  ! calculate the fractions of intercepted precipitation
  vegn_ifrac = cohort%cover

  ! get the lai
  vegn_lai = cohort%lai

  ! calculate the aerodynamic conductance coefficients
  call cana_turbulence(ustar, &
     cohort%cover, cohort%height, cohort%lai, cohort%sai, cohort%leaf_size, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v)

  ! take into account additional resistance of litter to the water vapor flux.
  ! not a good parameterization, but just using for sensitivity analyses now.
  ! ignores differing biomass and litter turnover rates.
  rav_lit = rav_lit_0 + rav_lit_vi * (cohort%lai+cohort%sai) &
                      + rav_lit_fsc * soil%fast_soil_C(1) &
                      + rav_lit_ssc * soil%slow_soil_C(1) &
                      + rav_lit_bwood * cohort%bwood
  con_g_v = con_g_v/(1.0+rav_lit*con_g_v)

  ! calculate the vegetation photosynthesis and associated stomatal conductance
  if (vegn_phot_co2_option == VEGN_PHOT_CO2_INTERACTIVE) then
     phot_co2 = cana_co2_mol
  else 
     phot_co2 = co2_for_photosynthesis
  endif
  call vegn_photosynthesis ( vegn, &
     SWdn(BAND_VIS), RSv(BAND_VIS), cana_q, phot_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply, &
     evap_demand, stomatal_cond, photosynt, photoresp )

  call get_vegn_wet_frac ( cohort, fw, DfwDwl, DfwDwf, fs, DfsDwl, DfsDwf )
  ! transpiring fraction and its derivatives
  ft     = 1 - fw - fs
  DftDwl = - DfwDwl - DfsDwl
  DftDwf = - DfwDwf - DfsDwf
  call qscomp(cohort%prog%Tv, p_surf, qvsat, DqvsatDTv)

  rho = p_surf/(rdgas*cana_T *(1+d608*cana_q))
  
  ! get the vegetation temperature
  vegn_T  =  cohort%prog%Tv
  ! get the amount of intercepted water and snow
  vegn_Wl =  cohort%prog%Wl
  vegn_Ws =  cohort%prog%Ws
  ! calculate the drip rates
  drip_l  = max(vegn_Wl,0.0)/tau_drip_l
  drip_s  = max(vegn_Ws,0.0)/tau_drip_s
  ! correct the drip rates so that the amount of water and snow accumulated over time step 
  ! is no larger then the canopy water-holding capacity
  drip_l = max((vegn_Wl+precip_l*delta_time*vegn_ifrac-cohort%Wl_max)/delta_time,drip_l)
  drip_s = max((vegn_Ws+precip_s*delta_time*vegn_ifrac-cohort%Ws_max)/delta_time,drip_s)

  ! calculate the total heat capacity
  call vegn_data_heat_capacity (cohort, vegn_hcap)
  vegn_hcap = vegn_hcap + clw*cohort%prog%Wl + csw*cohort%prog%Ws
  ! calculate the coefficient of sensible heat flux linearization
  Hv0     =  2*rho*cp_air*con_v_h*(cohort%prog%Tv - cana_T)
  DHvDTv  =  2*rho*cp_air*con_v_h
  DHvDTc  = -2*rho*cp_air*con_v_h
  ! calculate the coefficients of the transpiration linearization
  if(con_v_v==0.and.stomatal_cond==0) then
     total_cond = 0.0
  else
     total_cond = stomatal_cond*con_v_v/(stomatal_cond+con_v_v)
  endif

  if(qvsat>cana_q)then
     ! flux is directed from the surface: transpiration is possible, and the
     ! evaporation of intercepted water depends on the fraction of wet/snow
     ! covered canopy.

     ! prohibit transpiration if leaf temperature below some predefined minimum
     ! typically (268K, but check namelist)
     if(cohort%prog%Tv < T_transp_min) total_cond = 0 
     ! calculate the transpiration linearization coefficients
     Et0     =  rho*total_cond*ft*(qvsat - cana_q)
     DEtDTv  =  rho*total_cond*ft*DqvsatDTv
     DEtDqc  = -rho*total_cond*ft
     DEtDwl  =  rho*total_cond*DftDwl*(qvsat - cana_q)
     DEtDwf  =  rho*total_cond*DftDwf*(qvsat - cana_q)
     ! calculate the coefficients of the intercepted liquid evaporation linearization
     Eli0    =  rho*con_v_v*fw*(qvsat - cana_q)
     DEliDTv =  rho*con_v_v*fw*DqvsatDTv
     DEliDqc = -rho*con_v_v*fw
     DEliDwl =  rho*con_v_v*DfwDwl*(qvsat-cana_q)
     DEliDwf =  rho*con_v_v*DfwDwf*(qvsat-cana_q)
     ! calculate the coefficients of the intercepted snow evaporation linearization
     Efi0    =  rho*con_v_v*fs*(qvsat - cana_q)
     DEfiDTv =  rho*con_v_v*fs*DqvsatDTv
     DEfiDqc = -rho*con_v_v*fs
     DEfiDwl =  rho*con_v_v*DfsDwl*(qvsat-cana_q)
     DEfiDwf =  rho*con_v_v*DfsDwf*(qvsat-cana_q)
  else
     ! Flux is directed TOWARD the surface: no transpiration (assuming plants do not
     ! take water through stomata), and condensation does not depend on the fraction
     ! of wet canopy -- dew formation occurs on the entire surface

     ! prohibit transpiration:
     Et0     = 0
     DEtDTv  = 0; DEtDwl = 0; DEtDwf = 0;
     DEtDqc  = 0
     ! calculate dew or frost formation rates, depending on the temperature
     Eli0    = 0; Efi0    = 0
     DEliDTv = 0; DEfiDTv = 0
     DEliDqc = 0; DEfiDqc = 0
     DEliDwl = 0; DEfiDwl = 0
     DEliDwf = 0; DEfiDwf = 0
     ! calculate the coefficients of the intercepted liquid condensation linearization
     if(vegn_T >= tfreeze) then
        Eli0    =  rho*con_v_v*(qvsat - cana_q)
        DEliDTv =  rho*con_v_v*DqvsatDTv
        DEliDqc = -rho*con_v_v
     else
        ! calculate the coefficients of the intercepted snow condensation linearization
        Efi0    =  rho*con_v_v*(qvsat - cana_q)
        DEfiDTv =  rho*con_v_v*DqvsatDTv
        DEfiDqc = -rho*con_v_v
     endif
     ! prohibit switching from condensation to evaporation if the water content
     ! is below certain threshold
     if (vegn_Wl < min_Wl) then
        Eli0 = 0 ; DEliDTv = 0 ; DEliDqc = 0 ; DEliDwl = 0 ; DEliDwf = 0
     endif
     if (vegn_Ws < min_Ws) then
        Efi0 = 0 ; DEfiDTv = 0 ; DEfiDqc = 0 ; DEfiDwl = 0 ; DEfiDwf = 0
     endif
        
  endif
  ! ---- diagnostic section
  call send_tile_data(id_evap_demand, evap_demand, diag)
  call send_tile_data(id_stomatal, stomatal_cond, diag)
  call send_tile_data(id_an_op, cohort%An_op, diag)
  call send_tile_data(id_an_cl, cohort%An_cl, diag)
  call send_tile_data(id_con_v_h, con_v_h, diag)
  call send_tile_data(id_con_v_v, con_v_v, diag)
  call send_tile_data(id_phot_co2, phot_co2, diag)

end subroutine vegn_step_1


! ============================================================================
! Given the surface solution, substitute it back into the vegetation equations 
! to determine new vegetation state.
subroutine vegn_step_2 ( vegn, diag, &
     delta_Tv, delta_wl, delta_wf, &
     vegn_melt, &
     vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy, kg/(m2 s)
     vegn_ovfl_Hl, vegn_ovfl_Hs  ) ! heat flux carried from canopy by overflow, W/(m2 s)

  ! ---- arguments 
  type(vegn_tile_type) , intent(inout) :: vegn
  type(diag_buff_type) , intent(inout) :: diag
  real, intent(in) :: &
       delta_Tv, & ! change in vegetation temperature, degK
       delta_wl, & ! change in intercepted liquid water mass, kg/m2
       delta_wf    ! change in intercepted frozen water mass, kg/m2 
  real, intent(out) :: &
       vegn_melt, &
       vegn_ovfl_l,   vegn_ovfl_s,   & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs      ! heat flux from canopy due to overflow

  ! ---- local variables
  real :: &
     vegn_Wl_max, &  ! max. possible amount of liquid water in the canopy
     vegn_Ws_max, &  ! max. possible amount of solid water in the canopy
     mcv, &
     cap0, melt_per_deg, &
     Wl, Ws  ! positively defined amounts of water and snow on canopy
  type(vegn_cohort_type), pointer :: cohort
  
  ! get the pointer to the first (and, currently, the only) cohort
  cohort => vegn%cohorts(1)

  if (is_watch_point()) then
     write(*,*)'#### vegn_step_2 input ####'
     __DEBUG3__(delta_Tv, delta_wl, delta_wf)
     __DEBUG1__(cohort%prog%Tv)
  endif

  ! update vegetation state
  cohort%prog%Tv = cohort%prog%Tv + delta_Tv
  cohort%prog%Wl = cohort%prog%Wl + delta_wl
  cohort%prog%Ws = cohort%prog%Ws + delta_wf 

  call vegn_data_intrcptn_cap(cohort, vegn_Wl_max, vegn_Ws_max)
  call vegn_data_heat_capacity(cohort, mcv)


  ! ---- update for evaporation and interception -----------------------------
  cap0 = mcv + clw*cohort%prog%Wl + csw*cohort%prog%Ws

  if(is_watch_point()) then
     write (*,*)'#### vegn_step_2 #### 1'
     __DEBUG1__(cap0)
     __DEBUG1__(cohort%prog%Tv)
     __DEBUG2__(cohort%prog%Wl, cohort%prog%Ws)
  endif
  ! melt on the vegetation should probably be prohibited altogether, since
  ! the amount of melt or freeze calculated this way is severely underestimated 
  ! (depending on the overall vegetation heat capacity) which leads to extended 
  ! periods when the canopy temperature is fixed at freezing point.
  if (lm2) then 
     vegn_melt = 0
  else
     ! ---- freeze/melt of intercepted water
     ! heat capacity of leaf + intercepted water/snow _can_ go below zero if the 
     ! total water content goes below zero as a result of implicit time step.
     ! If it does, we just prohibit melt, setting it to zero.
     if(cap0 > 0)then
        melt_per_deg = cap0 / hlf
        if (cohort%prog%Ws>0 .and. cohort%prog%Tv>tfreeze) then
           vegn_melt =  min(cohort%prog%Ws, (cohort%prog%Tv-tfreeze)*melt_per_deg)
        else if (cohort%prog%Wl>0 .and. cohort%prog%Tv<tfreeze) then
           vegn_melt = -min(cohort%prog%Wl, (tfreeze-cohort%prog%Tv)*melt_per_deg)
        else
           vegn_melt = 0
        endif
        cohort%prog%Ws = cohort%prog%Ws - vegn_melt
        cohort%prog%Wl = cohort%prog%Wl + vegn_melt
        if (vegn_melt/=0) &
             cohort%prog%Tv = tfreeze + (cap0*(cohort%prog%Tv-tfreeze) - hlf*vegn_melt) &
             / ( cap0 + (clw-csw)*vegn_melt )
        vegn_melt = vegn_melt / delta_time
     else
        vegn_melt = 0
     endif
  endif

  if(is_watch_point()) then
     write (*,*)'#### vegn_step_2 #### 1'
     __DEBUG1__(cap0)
     __DEBUG1__(cohort%prog%Tv)
     __DEBUG3__(vegn_melt, cohort%prog%Wl, cohort%prog%Ws)
  endif

  ! ---- update for overflow -------------------------------------------------
  Wl = max(cohort%prog%Wl,0.0); Ws = max(cohort%prog%Ws,0.0)
  vegn_ovfl_l = max (0.,Wl-vegn_Wl_max)/delta_time
  vegn_ovfl_s = max (0.,Ws-vegn_Ws_max)/delta_time
  vegn_ovfl_Hl = clw*vegn_ovfl_l*(cohort%prog%Tv-tfreeze)
  vegn_ovfl_Hs = csw*vegn_ovfl_s*(cohort%prog%Tv-tfreeze)

  cohort%prog%Wl = cohort%prog%Wl - vegn_ovfl_l*delta_time
  cohort%prog%Ws = cohort%prog%Ws - vegn_ovfl_s*delta_time

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_2 output #####'
     __DEBUG3__(vegn_melt, vegn_ovfl_l, vegn_ovfl_s)
     __DEBUG2__(vegn_ovfl_Hl,vegn_ovfl_Hs)
  endif

  ! ---- diagnostic section
  call send_tile_data(id_temp,   cohort%prog%Tv, diag)
  call send_tile_data(id_wl,     cohort%prog%Wl, diag)
  call send_tile_data(id_ws,     cohort%prog%Ws, diag)

  call send_tile_data(id_height, cohort%height, diag)
  call send_tile_data(id_lai, cohort%lai, diag)
  call send_tile_data(id_sai, cohort%sai, diag)
  call send_tile_data(id_leaf_size, cohort%leaf_size, diag)
  call send_tile_data(id_root_density, cohort%root_density, diag)
  call send_tile_data(id_root_zeta, cohort%root_zeta, diag)
  call send_tile_data(id_rs_min, cohort%rs_min, diag)
  call send_tile_data(id_leaf_refl, cohort%leaf_refl, diag)
  call send_tile_data(id_leaf_tran, cohort%leaf_tran, diag)
  call send_tile_data(id_leaf_emis, cohort%leaf_emis, diag)
  call send_tile_data(id_snow_crit, cohort%snow_crit, diag)
  
end subroutine vegn_step_2


! ============================================================================
! do the vegetation calculations that require updated (end-of-timestep) values 
! of prognostic land variables
subroutine vegn_step_3(vegn, soil, cana_T, precip, vegn_fco2, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: cana_T ! canopy temperature, deg K
  real, intent(in) :: precip ! total (rain+snow) precipitation, kg/(m2 s)
  real, intent(out) :: vegn_fco2 ! co2 flux from vegetation, kg CO2/(m2 s)
  type(diag_buff_type), intent(inout) :: diag
  
  ! ---- local vars
  real :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real :: theta ! average soil wetness, unitless
  real :: psist ! psi stress index
  real :: depth_ave! depth for averaging soil moisture based on Jackson function for root distribution  
  real :: percentile = 0.95

  tsoil = soil_ave_temp (soil,soil_carbon_depth_scale)
  ! depth for 95% of root according to Jackson distribution 
  depth_ave = -log(1.-percentile)*vegn%cohorts(1)%root_zeta

  theta = soil_ave_theta1(soil, depth_ave)

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_3 drought input ####'
     __DEBUG3__(depth_ave, tsoil, theta)
  endif

  call vegn_carbon_int(vegn, soil, tsoil, theta, diag)
  ! decrease, if necessary, csmoke spending rate so that csmoke pool
  ! is never depleted below zero
  vegn%csmoke_rate = max( 0.0, &
       min( vegn%csmoke_rate, &
            vegn%csmoke_pool/dt_fast_yr)&
       )
  ! update smoke pool -- stored amount of carbon lost to fire
  vegn%csmoke_pool = vegn%csmoke_pool - &
       vegn%csmoke_rate*dt_fast_yr
  ! decrease harvested rates so that pools are not depleted below zero  
  vegn%harv_rate(:) = max( 0.0, &
                           min(vegn%harv_rate(:), vegn%harv_pool(:)/dt_fast_yr) &
                         )
  ! update harvested pools -- amounts of stored harvested carbon by category
  vegn%harv_pool(:) = vegn%harv_pool(:) - &
       vegn%harv_rate(:)*dt_fast_yr
  ! --- calculate total co2 flux from vegetation
  vegn_fco2 = -vegn%nep + vegn%csmoke_rate + sum(vegn%harv_rate(:))
  ! --- convert it to kg CO2/(m2 s)
  vegn_fco2 = vegn_fco2*mol_CO2/(mol_C*seconds_per_year)

  ! --- accumulate values for climatological averages
  vegn%tc_av     = vegn%tc_av + cana_T
  vegn%tsoil_av  = vegn%tsoil_av + tsoil
  vegn%precip_av = vegn%precip_av + precip
  if (xwilt_available) then
     theta = soil_ave_theta1(soil,depth_ave)
  else
     theta = soil_ave_theta0(soil,vegn%cohorts(1)%root_zeta)
  endif
  vegn%theta_av_phen = vegn%theta_av_phen + theta
  vegn%theta_av_fire = vegn%theta_av_fire + soil_ave_theta1(soil,depth_ave)
  psist = soil_psi_stress(soil,vegn%cohorts(1)%root_zeta)
  vegn%psist_av  = vegn%psist_av + psist

  vegn%n_accum   = vegn%n_accum+1
  
  call send_tile_data(id_theph, theta, diag)
  call send_tile_data(id_psiph, psist, diag)

end subroutine vegn_step_3


! ============================================================================
! update slow components of the vegetation model
subroutine update_vegn_slow( )

  ! ---- local vars ----------------------------------------------------------
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: i,j,k ! current point indices
  integer :: ii ! pool iterator
  integer :: n ! number of cohorts
  real    :: weight_ncm ! low-pass filter value for the number of cold months

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)

  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js) ; te = tail_elmt(lnd%tile_map)
  do while ( ce /= te )
     call get_elmt_indices(ce,i,j,k) ; call set_current_point(i,j,k) ! this is for debug output only
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     if (day1 /= day0) then
        call vegn_daily_npp(tile%vegn)
     endif

     ! monthly averaging
     if (month1 /= month0) then
        ! compute averages from accumulated monthly values 
        tile%vegn%tc_av     = tile%vegn%tc_av     / tile%vegn%n_accum
        tile%vegn%tsoil_av  = tile%vegn%tsoil_av  / tile%vegn%n_accum
        tile%vegn%theta_av_phen  = tile%vegn%theta_av_phen  / tile%vegn%n_accum
        tile%vegn%theta_av_fire  = tile%vegn%theta_av_fire  / tile%vegn%n_accum
	tile%vegn%psist_av  = tile%vegn%psist_av  / tile%vegn%n_accum
        tile%vegn%precip_av = tile%vegn%precip_av / tile%vegn%n_accum
        ! accumulate annual values
        tile%vegn%p_ann_acm = tile%vegn%p_ann_acm+tile%vegn%precip_av
        tile%vegn%t_ann_acm = tile%vegn%t_ann_acm+tile%vegn%tc_av
        if ( tile%vegn%tc_av < cold_month_threshold ) & 
             tile%vegn%ncm_acm = tile%vegn%ncm_acm+1
        tile%vegn%t_cold_acm = min(tile%vegn%t_cold_acm, tile%vegn%tc_av)

        tile%vegn%nmn_acm = tile%vegn%nmn_acm+1 ! increase the number of accumulated months
     endif

     ! annual averaging
     if (year1 /= year0) then
        ! The ncm smoothing is coded as a low-pass exponential filter. See, for example
        ! http://en.wikipedia.org/wiki/Low-pass_filter
        weight_ncm = 1/(1+tau_smooth_ncm)
        if(tile%vegn%nmn_acm /= 0) then
           ! calculate annual averages from accumulated values
           tile%vegn%p_ann  = tile%vegn%p_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_ann  = tile%vegn%t_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_cold = tile%vegn%t_cold_acm
           tile%vegn%ncm    = weight_ncm*tile%vegn%ncm_acm + (1-weight_ncm)*tile%vegn%ncm
           ! reset accumulated values
           tile%vegn%ncm_acm    = 0
           tile%vegn%p_ann_acm  = 0
           tile%vegn%t_ann_acm  = 0
           tile%vegn%t_cold_acm = HUGE(tile%vegn%t_cold_acm)
        endif
!!$        call calc_miami_npp(tile%vegn)
        tile%vegn%nmn_acm = 0
     endif

     if (year1 /= year0 .and. do_biogeography) then
        call vegn_biogeography(tile%vegn)
     endif

     if (month1 /= month0.and.do_patch_disturbance) then
        call update_fuel(tile%vegn,tile%soil%w_wilt(1)/tile%soil%pars%vwc_sat)
        ! assume that all layers are the same soil type and wilting is vertically homogeneous
     endif

     if (day1 /= day0 .and. do_cohort_dynamics) then
        n = tile%vegn%n_cohorts
        call send_tile_data(id_cgain,sum(tile%vegn%cohorts(1:n)%carbon_gain),tile%diag)
        call send_tile_data(id_closs,sum(tile%vegn%cohorts(1:n)%carbon_loss),tile%diag)
        call send_tile_data(id_wdgain,sum(tile%vegn%cohorts(1:n)%bwood_gain),tile%diag)
        call vegn_growth(tile%vegn)
        call vegn_nat_mortality(tile%vegn,tile%soil,86400.0)
     endif

     if  (month1 /= month0 .and. do_phenology) then
        call vegn_phenology (tile%vegn, tile%soil)
        ! assume that all layers are the same soil type and wilting is vertically homogeneous
     endif

     if (year1 /= year0 .and. do_patch_disturbance) then
        call vegn_disturbance(tile%vegn, tile%soil, seconds_per_year)
     endif

     if (year1 /= year0) then
        call vegn_harvesting(tile%vegn)
        tile%vegn%fsc_rate = tile%vegn%fsc_pool/fsc_pool_spending_time
        tile%vegn%ssc_rate = tile%vegn%ssc_pool/ssc_pool_spending_time
        where(harvest_spending_time(:)>0)
           tile%vegn%harv_rate(:) = &
                tile%vegn%harv_pool(:)/harvest_spending_time(:)
        elsewhere
           tile%vegn%harv_rate(:) = 0.0
        end where
     endif

     ! ---- diagnostic section
     call send_tile_data(id_t_ann,   tile%vegn%t_ann,   tile%diag)
     call send_tile_data(id_t_cold,  tile%vegn%t_cold,  tile%diag)
     call send_tile_data(id_lambda,  tile%vegn%lambda,  tile%diag)
     call send_tile_data(id_p_ann,   tile%vegn%p_ann,   tile%diag)
     call send_tile_data(id_ncm,     real(tile%vegn%ncm), tile%diag)
     call send_tile_data(id_afire,   tile%vegn%disturbance_rate(1), tile%diag)
     call send_tile_data(id_atfall,  tile%vegn%disturbance_rate(0), tile%diag)

     do ii = 1,N_HARV_POOLS
        call send_tile_data(id_harv_pool(ii),tile%vegn%harv_pool(ii),tile%diag)
        call send_tile_data(id_harv_rate(ii),tile%vegn%harv_rate(ii),tile%diag)
     enddo
     call send_tile_data(id_t_harv_pool,sum(tile%vegn%harv_pool(:)),tile%diag)
     call send_tile_data(id_t_harv_rate,sum(tile%vegn%harv_rate(:)),tile%diag)
     call send_tile_data(id_csmoke_pool,tile%vegn%csmoke_pool,tile%diag)
     call send_tile_data(id_csmoke_rate,tile%vegn%csmoke_rate,tile%diag)
     call send_tile_data(id_fsc_pool,tile%vegn%fsc_pool,tile%diag)
     call send_tile_data(id_fsc_rate,tile%vegn%fsc_rate,tile%diag)
     call send_tile_data(id_ssc_pool,tile%vegn%ssc_pool,tile%diag)
     call send_tile_data(id_ssc_rate,tile%vegn%ssc_rate,tile%diag)

     n=tile%vegn%n_cohorts
     call send_tile_data(id_bl,      sum(tile%vegn%cohorts(1:n)%bl),     tile%diag)
     call send_tile_data(id_blv,     sum(tile%vegn%cohorts(1:n)%blv),    tile%diag)
     call send_tile_data(id_br,      sum(tile%vegn%cohorts(1:n)%br),     tile%diag)
     call send_tile_data(id_bsw,     sum(tile%vegn%cohorts(1:n)%bsw),    tile%diag)
     call send_tile_data(id_bwood,   sum(tile%vegn%cohorts(1:n)%bwood),  tile%diag)
     call send_tile_data(id_btot,    sum(tile%vegn%cohorts(1:n)%bl    &
                                        +tile%vegn%cohorts(1:n)%blv   &
                                        +tile%vegn%cohorts(1:n)%br    &
                                        +tile%vegn%cohorts(1:n)%bsw   &
                                        +tile%vegn%cohorts(1:n)%bwood ), tile%diag)
     call send_tile_data(id_fuel,    tile%vegn%fuel, tile%diag)
     call send_tile_data(id_species, real(tile%vegn%cohorts(1)%species), tile%diag)
     call send_tile_data(id_status,  real(tile%vegn%cohorts(1)%status),  tile%diag)
     call send_tile_data(id_leaf_age,real(tile%vegn%cohorts(1)%leaf_age),  tile%diag)!ens

     ! carbon budget tracking
     call send_tile_data(id_fsc_in,  sum(tile%soil%fsc_in(:)),  tile%diag)
     call send_tile_data(id_fsc_out, tile%vegn%fsc_out, tile%diag)
     call send_tile_data(id_ssc_in,  sum(tile%soil%ssc_in(:)),  tile%diag)
     call send_tile_data(id_ssc_out, tile%vegn%ssc_out, tile%diag)
     call send_tile_data(id_veg_in,  tile%vegn%veg_in,  tile%diag)
     call send_tile_data(id_veg_out, tile%vegn%veg_out, tile%diag)
     ! ---- end of diagnostic section

     ! reset averages and number of steps to 0 before the start of new month
     if (month1 /= month0) then
        tile%vegn%n_accum  = 0
        tile%vegn%tc_av    = 0.
        tile%vegn%tsoil_av = 0.
        tile%vegn%theta_av_phen = 0.
        tile%vegn%theta_av_fire = 0.
        tile%vegn%psist_av = 0.
        tile%vegn%precip_av= 0.
     endif

     !reset fuel and drought months before the start of new year
     if (year1 /= year0) then
        tile%vegn%lambda     = 0
        tile%vegn%fuel       = 0
     endif

  enddo

  ! seed transport
  if (year1 /= year0 .and. do_seed_transport) then
     call vegn_seed_transport()
  endif

  ! override with static vegetation
  if(day1/=day0) &
       call  read_static_vegn(lnd%time)
end subroutine update_vegn_slow


! ============================================================================
subroutine vegn_seed_transport()

  ! local vars
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: i,j ! current point indices
  real :: total_seed_supply
  real :: total_seed_demand
  real :: f_supply ! fraction of the supply that gets spent
  real :: f_demand ! fraction of the demand that gets satisfied

  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js) ; te = tail_elmt(lnd%tile_map)
  total_seed_supply = 0.0; total_seed_demand = 0.0
  do while ( ce /= te )
     call get_elmt_indices(ce,i,j)
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     total_seed_supply = total_seed_supply + vegn_seed_supply(tile%vegn)*tile%frac*lnd%area(i,j)
     total_seed_demand = total_seed_demand + vegn_seed_demand(tile%vegn)*tile%frac*lnd%area(i,j)
  enddo
  ! sum totals globally
  call mpp_sum(total_seed_demand, pelist=lnd%pelist)
  call mpp_sum(total_seed_supply, pelist=lnd%pelist)
  ! if either demand or supply are zeros we don't need (or can't) transport anything
  if (total_seed_demand==0.or.total_seed_supply==0)then
     return
  end if

  ! calculate the fraction of the supply that's going to be used
  f_supply = MIN(total_seed_demand/total_seed_supply, 1.0)
  ! calculate the fraction of the demand that's going to be satisfied
  f_demand = MIN(total_seed_supply/total_seed_demand, 1.0)
  ! note that either f_supply or f_demand is 1; the mass conservation law in the
  ! following calculations is satisfied since 
  ! f_demand*total_seed_demand - f_supply*total_seed_supply == 0

  ! redistribute part (or possibly all) of the supply to satisfy part (or possibly all) 
  ! of the demand
  ce = first_elmt(lnd%tile_map) ; te = tail_elmt(lnd%tile_map)
  do while ( ce /= te )
     call get_elmt_indices(ce,i,j)
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body
     
     call vegn_add_bliving(tile%vegn, &
          f_demand*vegn_seed_demand(tile%vegn)-f_supply*vegn_seed_supply(tile%vegn))
  enddo
end subroutine vegn_seed_transport


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_VEGN_ACCESSOR_0D(xtype,x) subroutine vegn_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x;endif;end subroutine

#define DEFINE_VEGN_ACCESSOR_1D(xtype,x) subroutine vegn_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x;endif;end subroutine

#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;end subroutine

#define DEFINE_COHORT_COMPONENT_ACCESSOR(xtype,component,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%component%x;end subroutine

DEFINE_VEGN_ACCESSOR_0D(integer,landuse)
DEFINE_VEGN_ACCESSOR_0D(real,age)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate)
DEFINE_VEGN_ACCESSOR_0D(real,tc_av)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_phen)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_fire)
DEFINE_VEGN_ACCESSOR_0D(real,psist_av)
DEFINE_VEGN_ACCESSOR_0D(real,tsoil_av)
DEFINE_VEGN_ACCESSOR_0D(real,precip_av)
DEFINE_VEGN_ACCESSOR_0D(real,fuel)
DEFINE_VEGN_ACCESSOR_0D(real,lambda)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold)
DEFINE_VEGN_ACCESSOR_0D(real,ncm)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold_acm)
DEFINE_VEGN_ACCESSOR_0D(real,ncm_acm)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_pool)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_rate)

DEFINE_VEGN_ACCESSOR_1D(real,harv_pool)
DEFINE_VEGN_ACCESSOR_1D(real,harv_rate)

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(integer,status)
DEFINE_COHORT_ACCESSOR(real,leaf_age)
DEFINE_COHORT_ACCESSOR(real,npp_previous_day)

DEFINE_COHORT_COMPONENT_ACCESSOR(real,prog,tv)
DEFINE_COHORT_COMPONENT_ACCESSOR(real,prog,wl)
DEFINE_COHORT_COMPONENT_ACCESSOR(real,prog,ws)

DEFINE_COHORT_ACCESSOR(real,height)

end module vegetation_mod
