! ============================================================================
! top-level core of the Land Dynamics (LaD) model code
! ============================================================================
module land_model_mod

#include "shared/debug.inc"
#include "shared/concat.inc"

use time_manager_mod, only : time_type, get_time, increment_time, time_type_to_real, &
     operator(+)
use mpp_domains_mod, only : domain2d, mpp_get_ntile_count

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use mpp_mod, only : mpp_max, mpp_sum
use fms_mod, only : write_version_number, error_mesg, FATAL, WARNING, NOTE, mpp_pe, &
     mpp_root_pe, file_exist, check_nml_error, close_file, &
     stdlog, stderr, mpp_clock_id, mpp_clock_begin, mpp_clock_end, string, &
     stdout, CLOCK_FLAG_DEFAULT, CLOCK_COMPONENT, CLOCK_ROUTINE
use field_manager_mod, only : MODEL_LAND
use data_override_mod, only : data_override
use diag_manager_mod, only : diag_axis_init, register_static_field, &
     register_diag_field, send_data
use constants_mod, only : radius, hlf, hlv, hls, tfreeze, pi, rdgas, rvgas, cp_air, &
     stefan
use astronomy_mod, only : diurnal_solar
use sphum_mod, only : qscomp
use tracer_manager_mod, only : NO_TRACER

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, mol_air, mol_C, mol_co2
use glacier_mod, only : read_glac_namelist, glac_init, glac_end, glac_get_sfc_temp, &
     glac_radiation, glac_diffusion, glac_step_1, glac_step_2, save_glac_restart
use lake_mod, only : read_lake_namelist, lake_init, lake_end, lake_get_sfc_temp, &
     lake_radiation, lake_diffusion, lake_step_1, lake_step_2, save_lake_restart
use soil_mod, only : read_soil_namelist, soil_init, soil_end, soil_get_sfc_temp, &
     soil_radiation, soil_diffusion, soil_step_1, soil_step_2, soil_step_3, &
     save_soil_restart
use snow_mod, only : read_snow_namelist, snow_init, snow_end, snow_get_sfc_temp, &
     snow_radiation, snow_diffusion, snow_get_depth_area, snow_step_1, snow_step_2, &
     save_snow_restart
use vegetation_mod, only : read_vegn_namelist, vegn_init, vegn_end, vegn_get_cover, &
     vegn_radiation, vegn_diffusion, vegn_step_1, vegn_step_2, vegn_step_3, &
     update_vegn_slow, save_vegn_restart
use cana_tile_mod, only : canopy_air_mass, canopy_air_mass_for_tracers, cana_tile_heat
use canopy_air_mod, only : read_cana_namelist, cana_init, cana_end, cana_state,&
     cana_step_1, cana_step_2, cana_radiation, cana_roughness, &
     save_cana_restart
use river_mod, only : river_init, river_end, update_river, river_stock_pe, &
     save_river_restart
use topo_rough_mod, only : topo_rough_init, topo_rough_end, update_topo_rough
use soil_tile_mod, only : soil_cover_cold_start, soil_tile_stock_pe, &
                          soil_tile_heat
use vegn_tile_mod, only : vegn_cover_cold_start, vegn_data_rs_min, &
                          update_derived_vegn_data, vegn_tile_stock_pe, &
                          vegn_tile_heat
use lake_tile_mod, only : lake_cover_cold_start, lake_tile_stock_pe, &
                          lake_tile_heat
use glac_tile_mod, only : glac_pars_type, glac_cover_cold_start, &
                          glac_tile_stock_pe, glac_tile_heat
use snow_tile_mod, only : snow_tile_stock_pe, snow_tile_heat
use land_numerics_mod, only : ludcmp, lubksb, nearest, &
     horiz_remap_type, horiz_remap_new, horiz_remap, horiz_remap_del, &
     horiz_remap_print
use land_io_mod, only : read_land_io_namelist, input_buf_size
use land_tile_mod, only : land_tile_type, land_tile_list_type, &
     land_tile_enum_type, new_land_tile, insert, nitems, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=), &
     get_elmt_indices, get_tile_tags, land_tile_carbon
use land_data_mod, only : land_data_type, atmos_land_boundary_type, &
     land_state_type, land_data_init, land_data_end, lnd, &
     dealloc_land2cplr, realloc_land2cplr, &
     dealloc_cplr2land, realloc_cplr2land, &
     land_data_type_chksum, atm_lnd_bnd_type_chksum
use nf_utils_mod,  only : nfu_inq_var, nfu_inq_dim, nfu_get_var
use land_utils_mod, only : put_to_tiles_r0d_fptr
use land_tile_io_mod, only : print_netcdf_error, create_tile_out_file, &
    read_tile_data_r0d_fptr, write_tile_data_r0d_fptr, &
    write_tile_data_i0d_fptr, get_input_restart_name
use land_tile_diag_mod, only : tile_diag_init, tile_diag_end, &
    register_tiled_diag_field, send_tile_data, dump_tile_diag_fields, &
    add_tiled_diag_field_alias, &
    OP_AVERAGE, OP_SUM
use land_debug_mod, only : land_debug_init, land_debug_end, set_current_point, &
     is_watch_point, get_watch_point, check_temp_range, current_face, &
     get_current_point
use static_vegn_mod, only : write_static_vegn
use land_transitions_mod, only : &
     land_transitions_init, land_transitions_end, land_transitions, &
     save_land_transitions_restart
use stock_constants_mod, only: ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT

implicit none
private

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! saves the land model restart(s)
public update_land_model_fast   ! time-step integration
public update_land_model_slow   ! time-step integration
public atmos_land_boundary_type ! data from coupler to land
public land_data_type           ! data from land to coupler
public land_data_type_chksum    ! routine to print checksums for land_data_type
public atm_lnd_bnd_type_chksum  ! routine to print checksums for atmos_land_boundary_type

public :: Lnd_stock_pe          ! return stocks of conservative quantities
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land', &
     version     = '$Id: land_model.F90,v 20.0 2013/12/13 23:29:26 fms Exp $', &
     tagname     = '$Name: tikal $'

! ==== module variables ======================================================

! ---- namelist --------------------------------------------------------------
logical :: use_old_conservation_equations  = .false.
logical :: lm2                             = .false.
logical :: do_age                          = .false.
logical :: give_stock_details              = .false.
logical :: use_tfreeze_in_grnd_latent      = .false.
logical :: use_atmos_T_for_precip_T        = .false.
logical :: use_atmos_T_for_evap_T          = .false.
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
real    :: min_sum_lake_frac = 1.e-8
real    :: min_frac = 0.0 ! minimum fraction of soil, lake, and glacier that is not discarded on cold start
real    :: gfrac_tol         = 1.e-6
real    :: discharge_tol = -1.e20
real    :: con_fac_large = 1.e6
real    :: con_fac_small = 1.e-6
integer :: num_c = 0
real    :: tau_snow_T_adj = -1.0 ! time scale of snow temperature adjustment
              ! for the snow-free surface (s); negative means no adjustment
logical :: prohibit_negative_canopy_water = .FALSE. ! if true, then in case of negative canopy
              ! water the evaporation is fixed and the equations are re-solved.
              ! Default retrievs old behavior.
character(16) :: nearest_point_search = 'global' ! specifies where to look for
              ! nearest points for missing data, "global" or "face"
logical :: print_remapping = .FALSE. ! if true, full land cover remapping
              ! information is printed on the cold start
integer :: layout(2) = (/0,0/)
integer :: io_layout(2) = (/0,0/)
namelist /land_model_nml/ use_old_conservation_equations, &
                          lm2, do_age, give_stock_details, &
                          use_tfreeze_in_grnd_latent, &
                          use_atmos_T_for_precip_T, &
                          use_atmos_T_for_evap_T, &
                          cpw, clw, csw, min_sum_lake_frac, min_frac, &
                          gfrac_tol, discharge_tol, &
                          con_fac_large, con_fac_small, num_c, &
                          tau_snow_T_adj, prohibit_negative_canopy_water, &
                          nearest_point_search, print_remapping, &
                          layout, io_layout
! ---- end of namelist -------------------------------------------------------

logical  :: module_is_initialized = .FALSE.
logical  :: stock_warning_issued  = .FALSE.
logical  :: update_cana_co2 ! if false, cana_co2 is not updated during the model run.
character(len=256) :: grid_spec_file="INPUT/grid_spec.nc" 
real     :: delta_time ! duration of main land time step (s)
integer  :: num_species
integer  :: num_phys = 2
real,    allocatable :: frac           (:,:)    ! fraction of land in cells
logical, allocatable :: river_land_mask(:,:), missing_rivers(:,:)
real,    allocatable :: no_riv(:,:)

! ---- diag field IDs --------------------------------------------------------
integer :: &
 ! COLUMN        VEGN        SNOW      GLAC/LAKE/SOIL  CANOPY-AIR  RIVER
  id_VWS,                                               id_VWSc,           &
  id_LWS,      id_LWSv,     id_LWSs,     id_LWSg,                          &
  id_FWS,      id_FWSv,     id_FWSs,     id_FWSg,                          &
  id_HS,       id_HSv,      id_HSs,      id_HSg,        id_HSc,            &
  id_precip,                                                               &
  id_hprec,                                                                &
  id_lprec,    id_lprecv,   id_lprecs,   id_lprecg,                        &
  id_hlprec,   id_hlprecv,  id_hlprecs,  id_hlprecg,                       &
  id_fprec,    id_fprecv,   id_fprecs,                                     &
  id_hfprec,   id_hfprecv,  id_hfprecs,                                    &
  id_evap,                                                                 &
  id_hevap,                                                                &
  id_levap,    id_levapv,   id_levaps,   id_levapg,                        &
  id_hlevap,   id_hlevapv,  id_hlevaps,  id_hlevapg,                       &
  id_fevap,    id_fevapv,   id_fevaps,   id_fevapg,                        &
  id_hfevap,   id_hfevapv,  id_hfevaps,  id_hfevapg,                       &
  id_runf,                                                                 &
  id_hrunf,                                                                &
  id_lrunf,                 id_lrunfs,   id_lrunfg,                        &
  id_hlrunf,                id_hlrunfs,  id_hlrunfg,                       &
  id_frunf,                 id_frunfs,                                     &
  id_hfrunf,                id_hfrunfs,                                    &
  id_melt,     id_meltv,    id_melts,    id_meltg,                         &
  id_fsw,      id_fswv,     id_fsws,     id_fswg,                          &
  id_flw,      id_flwv,     id_flws,     id_flwg,                          &
  id_sens,     id_sensv,    id_senss,    id_sensg,                         &
!
  id_e_res_1,  id_e_res_2,  id_cd_m,     id_cd_t,                          &
  id_cellarea, id_landarea, id_landfrac, id_no_riv,                        &
  id_geolon_t, id_geolat_t,                                                &
  id_frac,     id_area,     id_ntiles,                                     &
  id_dis_liq,  id_dis_ice,  id_dis_heat, id_dis_sink,                      &
  id_z0m,      id_z0s,      id_con_g_h,                                    &
  id_transp,                id_wroff,    id_sroff,                         &
  id_htransp,  id_huptake,  id_hroff,    id_gsnow,    id_gequil,           &
  id_grnd_flux,                                                            &
  id_soil_water_supply,     id_levapg_max,                                 &
  id_water,    id_snow,                                                    &
  id_Trad,     id_Tca,      id_qca,      id_qco2,     id_qco2_dvmr,        &
  id_swdn_dir, id_swdn_dif, id_swup_dir, id_swup_dif, id_lwdn,             &
  id_fco2,                                                                 &
  id_vegn_cover,    id_cosz,                                               &
  id_albedo_dir,    id_albedo_dif,                                         &
  id_vegn_refl_dir, id_vegn_refl_dif, id_vegn_refl_lw,                     &
  id_vegn_tran_dir, id_vegn_tran_dif, id_vegn_tran_lw,                     &
  id_vegn_sctr_dir,                                                        &
  id_subs_refl_dir, id_subs_refl_dif, id_subs_emis, id_grnd_T, id_total_C

! ---- global clock IDs
integer :: landClock, landFastClock, landSlowClock


! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains


! ============================================================================
subroutine land_model_init &
     (cplr2land, land2cplr, time_init, time, dt_fast, dt_slow)
! initialize land model using grid description file as an input. This routine
! reads land grid boundaries and area of land from a grid description file

! NOTES: theoretically, the grid description file can specify any regular
! rectangular grid for land, not just lon/lat grid. Therefore the variables
! "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
! boundaries of the grid.
!   However, at this time the module land_properties assumes that grid _is_
! lon/lat and therefore the entire module also have to assume that the land 
! grid is lon/lat.
!   lon/lat grid is also assumed for the diagnostics, but this is probably not
! so critical. 
  type(atmos_land_boundary_type), intent(inout) :: cplr2land ! boundary data
  type(land_data_type)          , intent(inout) :: land2cplr ! boundary data
  type(time_type), intent(in) :: time_init ! initial time of simulation (?)
  type(time_type), intent(in) :: time      ! current time
  type(time_type), intent(in) :: dt_fast   ! fast time step
  type(time_type), intent(in) :: dt_slow   ! slow time step

  ! ---- local vars ----------------------------------------------------------
  integer :: ncid, varid
  integer :: unit, ierr, io
  integer :: id_lon, id_lat, id_band     ! IDs of land diagnostic axes
  logical :: used                        ! return value of send_data diagnostics routine
  integer :: i,j,k
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce, te
  character(len=256) :: restart_file_name
  logical :: restart_exists
  ! IDs of local clocks
  integer :: landInitClock

  module_is_initialized = .TRUE.

  ! [1] print out version number
  call write_version_number (version, tagname)

  ! initialize land model clocks
  landClock      = mpp_clock_id('Land'               ,CLOCK_FLAG_DEFAULT,CLOCK_COMPONENT)
  landFastClock  = mpp_clock_id('Update-Land-Fast'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landSlowClock  = mpp_clock_id('Update-Land-Slow'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landInitClock  = mpp_clock_id('Land init'          ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)

  call mpp_clock_begin(landInitClock)

  ! [ ] initialize land debug output
  call land_debug_init()

  ! [ ] initialize tile-specific diagnostics internals
  call tile_diag_init()

  ! [2] read namelists
  ! [2.1] read land model namelist
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_model_nml, iostat=io)
     ierr = check_nml_error(io, 'land_model_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_model_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_model_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_model_nml)
     call close_file (unit)
  endif
  ! [2.2] read sub-model namelists: then need to be read before initialization
  ! because they can affect the way cover and tiling is initialized on cold start.
  call read_land_io_namelist()
  call read_soil_namelist()
  call read_vegn_namelist()
  call read_lake_namelist()
  call read_glac_namelist()
  call read_snow_namelist()
  call read_cana_namelist()

  ! [ ] initialize land state data, including grid geometry and processor decomposition
  call land_data_init(layout, io_layout, time, dt_fast, dt_slow)
  delta_time  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience

  ! calculate land fraction
  allocate(frac(lnd%is:lnd%ie,lnd%js:lnd%je))
  frac = lnd%area/lnd%cellarea

  ! [5] initialize tiling
  call get_input_restart_name('INPUT/land.res.nc',restart_exists,restart_file_name)
  if(restart_exists) then
     call error_mesg('land_model_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     ! read map of tiles -- retrieve information from 
     call land_cover_warm_start(restart_file_name,lnd)
     ! initialize land model data
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,ncid))
     if (nf_inq_varid(ncid,'lwup',varid)==NF_NOERR) &
          call read_tile_data_r0d_fptr(ncid,'lwup',land_lwup_ptr)
     if (nf_inq_varid(ncid,'e_res_1',varid)==NF_NOERR) &
          call read_tile_data_r0d_fptr(ncid,'e_res_1',land_e_res_1_ptr)
     if (nf_inq_varid(ncid,'e_res_2',varid)==NF_NOERR) &
          call read_tile_data_r0d_fptr(ncid,'e_res_2',land_e_res_2_ptr)
     __NF_ASRT__(nf_close(ncid))
  else
     ! initialize map of tiles -- construct it by combining tiles
     ! from component models
     call error_mesg('land_model_init',&
          'cold-starting land cover map',&
          NOTE)
     call land_cover_cold_start(lnd)
  endif

  ! [6] initialize land model diagnostics -- must be before *_data_init so that
  ! *_data_init can write static fields if necessary
  call land_diag_init( lnd%coord_glonb, lnd%coord_glatb, lnd%coord_glon, lnd%coord_glat, time, lnd%domain, &
       id_lon, id_lat, id_band )
  ! set the land diagnostic axes ids for the flux exchange
  land2cplr%axes = (/id_lon,id_lat/)
  ! send some static diagnostic fields to output
  if ( id_cellarea > 0 ) used = send_data ( id_cellarea, lnd%cellarea, lnd%time )
  if ( id_landarea > 0 ) used = send_data ( id_landarea, lnd%area, lnd%time )
  if ( id_landfrac > 0 ) used = send_data ( id_landfrac, frac,     lnd%time )
  if ( id_geolon_t > 0 ) used = send_data ( id_geolon_t, lnd%lon*180.0/PI, lnd%time )
  if ( id_geolat_t > 0 ) used = send_data ( id_geolat_t, lnd%lat*180.0/PI, lnd%time )

  ! [7] initialize individual sub-models
  num_species = num_phys + num_c
  if (do_age) num_species = num_species + 1

  call soil_init ( id_lon, id_lat, id_band )
  call vegn_init ( id_lon, id_lat, id_band )
  call lake_init ( id_lon, id_lat )
  call glac_init ( id_lon, id_lat )
  call snow_init ( id_lon, id_lat )
  call cana_init ( id_lon, id_lat )
  call topo_rough_init( lnd%time, lnd%lonb, lnd%latb, &
       lnd%domain, id_lon, id_lat)
  allocate (river_land_mask(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate ( missing_rivers(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate ( no_riv        (lnd%is:lnd%ie,lnd%js:lnd%je))
  call river_init( lnd%lon, lnd%lat, &
                   lnd%time, lnd%dt_fast, lnd%domain,     &
                   frac, &
                   id_lon, id_lat,                        &
                   river_land_mask                        )
  missing_rivers = frac.gt.0. .and. .not.river_land_mask
  no_riv = 0.
  where (missing_rivers) no_riv = 1.
  if ( id_no_riv > 0 ) used = send_data( id_no_riv, no_riv, lnd%time )
  call land_transitions_init (id_lon, id_lat) 
  ! [8] initialize boundary data
  ! [8.1] allocate storage for the boundary data 
  call realloc_land2cplr ( land2cplr )
  call realloc_cplr2land ( cplr2land )
  ! [8.2] set the land mask to FALSE everywhere -- update_land_bc_fast
  ! will set it to true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0
  ! [8.3] get the current state of the land boundary for the coupler
  ce = first_elmt(lnd%tile_map,                  &
              is=lbound(cplr2land%t_flux,1), &
              js=lbound(cplr2land%t_flux,2)  )
  te = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     ! calculate indices of the current tile in the input arrays;
     ! assume all the cplr2land components have the same lbounds
     call get_elmt_indices(ce,i,j,k)
     ! set this point coordinates as current for debug output
     call set_current_point(i,j,k)
     ! get pointer to current tile
     tile => current_tile(ce)
     ! advance enumerator to the next tile
     ce=next_elmt(ce)

     call update_land_bc_fast (tile, i,j,k, land2cplr, is_init=.true.)
  enddo

  ! [8.4] update topographic roughness scaling
  call update_land_bc_slow( land2cplr )

  ! mask error checking
  do j=lnd%js,lnd%je
  do i=lnd%is,lnd%ie
     if(frac(i,j)>0.neqv.ANY(land2cplr%mask(i,j,:))) then
        call error_mesg('land_model_init','masks are not equal',FATAL)
     endif
  enddo
  enddo

  ! [9] check the properties of co2 exchange with the atmosphere and set appropriate
  ! flags
  if (canopy_air_mass_for_tracers==0.and.lnd%ico2==NO_TRACER) then
     call error_mesg('land_model_init', &
          'canopy_air_mass_for_tracers is set to zero, and CO2 exchange with the atmosphere is not set up: '// &
          'canopy air CO2 concentration will not be updated',NOTE)
     update_cana_co2 = .FALSE.
  else
     update_cana_co2 = .TRUE.
  end if

  call mpp_clock_end(landInitClock)

end subroutine land_model_init


! ============================================================================
subroutine land_model_end (cplr2land, land2cplr)
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  
  module_is_initialized = .FALSE.

  call error_mesg('land_model_end','writing NetCDF restart',NOTE)
  call land_model_restart()

  ! we still want to call the *_end procedures for component models, even
  ! if the number of tiles in this domain is zero, in case they are doing 
  ! something else besides saving the restart, of if they want to save
  ! restart anyway
  call land_transitions_end()
  call glac_end ()
  call lake_end ()
  call soil_end ()
  call snow_end ()
  call vegn_end ()
  call cana_end ()
  call topo_rough_end()
  call river_end()

  deallocate(frac)
  deallocate(river_land_mask, missing_rivers, no_riv)
  call dealloc_land2cplr(land2cplr, dealloc_discharges=.TRUE.)
  call dealloc_cplr2land(cplr2land)

  call tile_diag_end()

  ! deallocate tiles
  call land_data_end()

  ! finish up the land debugging diagnostics
  call land_debug_end
  
end subroutine land_model_end


! ============================================================================
! write land model restarts 
subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp ! timestamp to add to the file name
  
  ! ---- local vars
  integer :: tile_dim_length ! length of tile dimension in output files 
                             ! global max of number of tiles per gridcell 
  integer :: i,j,k
  integer :: unit ! netcdf id of the restart file
  character(256) :: timestamp_

  ! [1] count all land tiles and determine the length of tile dimension
  ! sufficient for the current domain
  tile_dim_length = 0
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     k = nitems(lnd%tile_map(i,j))
     tile_dim_length = max(tile_dim_length,k)
  enddo
  enddo

  ! [2] calculate the tile dimension length by taking the max across all domains
  call mpp_max(tile_dim_length)
  if (tile_dim_length==0) then
     call error_mesg('land_model_restart',&
       'No land points exist (tile_dim_length=0), therefore no land restarts will be saved',&
       WARNING)
     return
  endif
   
  ! [3] create tile output file
  timestamp_=''
  if (present(timestamp)) then
     if(trim(timestamp)/='') timestamp_=trim(timestamp)//'.'
  endif
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp_)//'land.res.nc', &
       lnd%coord_glon, lnd%coord_glat, land_tile_exists, tile_dim_length)
     
  ! [4] write data fields
  ! write fractions and tile tags
  call write_tile_data_r0d_fptr(unit,'frac',land_frac_ptr,'fractional area of tile')
  call write_tile_data_i0d_fptr(unit,'glac',glac_tag_ptr,'tag of glacier tiles')
  call write_tile_data_i0d_fptr(unit,'lake',lake_tag_ptr,'tag of lake tiles')
  call write_tile_data_i0d_fptr(unit,'soil',soil_tag_ptr,'tag of soil tiles')
  call write_tile_data_i0d_fptr(unit,'vegn',vegn_tag_ptr,'tag of vegetation tiles')
  ! write the upward long-wave flux 
  call write_tile_data_r0d_fptr(unit,'lwup',land_lwup_ptr,'upward long-wave flux')
  ! write energy residuals
  call write_tile_data_r0d_fptr(unit,'e_res_1',land_e_res_1_ptr,&
       'energy residual in canopy air energy balance equation', 'W/m2')
  call write_tile_data_r0d_fptr(unit,'e_res_2',land_e_res_2_ptr,&
       'energy residual in canopy energy balance equation', 'W/m2')
  
  ! [5] close file
  __NF_ASRT__(nf_close(unit))

  ! [6] save component models' restarts
  call save_land_transitions_restart(timestamp_)
  call save_glac_restart(tile_dim_length,timestamp_)
  call save_lake_restart(tile_dim_length,timestamp_)
  call save_soil_restart(tile_dim_length,timestamp_)
  call save_snow_restart(tile_dim_length,timestamp_)
  call save_vegn_restart(tile_dim_length,timestamp_)
  call save_cana_restart(tile_dim_length,timestamp_)
  call save_river_restart(timestamp_)

end subroutine land_model_restart

! ============================================================================
subroutine land_cover_cold_start(lnd)
  type(land_state_type), intent(inout) :: lnd

  ! ---- local vars
  real, dimension(:,:,:), pointer :: &
       glac, soil, lake, vegn ! arrays of fractions for respective sub-models
  logical, dimension(lnd%ie-lnd%is+1,lnd%je-lnd%js+1) :: &
       land_mask, valid_data, invalid_data
  integer :: iwatch,jwatch,kwatch,face
  integer :: i,j
  integer :: ps,pe ! boundaries of PE list for remapping
  type(horiz_remap_type) :: map

  ! calculate the global land mask
  land_mask = lnd%area > 0

  ! get the global maps of fractional covers for each of the sub-models
  glac=>glac_cover_cold_start(land_mask,lnd%lonb,lnd%latb)
  lake=>lake_cover_cold_start(land_mask,lnd%lonb,lnd%latb,lnd%domain)
  soil=>soil_cover_cold_start(land_mask,lnd%lonb,lnd%latb)
  vegn=>vegn_cover_cold_start(land_mask,lnd%lonb,lnd%latb)

  ! remove any input lake fraction in coastal cells
  where (frac.lt. 1.-gfrac_tol) lake(:,:,1) = 0.
  ! NOTE that the lake area in the coastal cells can be set to non-zero
  ! again by the "ground fraction reconciliation code" below. Strictly
  ! speaking the above line of code should be replaced with the section
  ! commented out with "!-zero" below, but we preserve the old way to avoid
  ! backward incompatibility with older runs. This needs updating in the
  ! future when the decision about what to do with lakes in coastal cells is
  ! made.

  ! reconcile ground fractions with the land mask within compute domain
  valid_data = land_mask.and.(sum(glac,3)+sum(lake,3)+sum(soil,3)>0)
  invalid_data = land_mask.and..not.valid_data

  call get_watch_point(iwatch,jwatch,kwatch,face)
  if (face==lnd%face.and.(lnd%is<=iwatch.and.iwatch<=lnd%ie).and.(lnd%js<=jwatch.and.jwatch<=lnd%je)) then
     write(*,*)'###### land_cover_cold_start: input data #####'
     write(*,'(99(a,i4.2,x))')'iwatch=',iwatch,'jwatch=',jwatch,'face=',lnd%face
     write(*,'(99(a,g23.16,x))')'lon=',lnd%lon(iwatch,jwatch)*180/PI,'lat=',lnd%lat(iwatch,jwatch)*180/PI
     ! calculate local compute domain indices; we assume glac,lake,soil,vegn all
     ! have the same lbounds
     i = iwatch-lnd%is+lbound(glac,1); j = jwatch-lnd%js+lbound(glac,2)
     __DEBUG2__(lnd%is,lnd%js)
     write(*,'(a,99(a,i4.2,x))')'local indices:','i=',i,'j=',j
     __DEBUG3__(frac(iwatch,jwatch),land_mask(i,j),valid_data(i,j))
     __DEBUG1__(glac(i,j,:))
     __DEBUG1__(lake(i,j,:))
     __DEBUG1__(soil(i,j,:))
     __DEBUG1__(vegn(i,j,:))
  endif

  if (trim(nearest_point_search)=='global') then
     ps=0 ; pe=size(lnd%pelist)-1
  else if (trim(nearest_point_search)=='face') then
     ! this assumes that the number of PEs is divisible by the number of
     ! mosaic faces. lnd%pelist starts with 0
     ps = size(lnd%pelist)/lnd%nfaces*(lnd%face-1)
     pe = size(lnd%pelist)/lnd%nfaces*lnd%face - 1
  else
     call error_mesg('land_cover_cold_start',&
          'option nearest_point_search="'//trim(nearest_point_search)//&
          '" is illegal, use "global" or "face"',&
          FATAL)
  endif
  call horiz_remap_new(invalid_data,valid_data,lnd%lon,lnd%lat,lnd%domain,&
          lnd%pelist(ps:pe),map)
  if (print_remapping) call horiz_remap_print(map,'land cover remap:')
  call horiz_remap(map,lnd%domain,glac)
  call horiz_remap(map,lnd%domain,lake)
  call horiz_remap(map,lnd%domain,soil)
  call horiz_remap_del(map)

!-zero  ! remove any input lake fraction in coastal cells
!-zero  do j = lnd%js,lnd%je
!-zero  do i = lnd%is,lnd%ie
!-zero     call set_current_point(i,j,1)
!-zero     if (frac(i,j) < 1-gfrac_tol) then
!-zero        lake(i,j,:) = 0.0
!-zero        if(is_watch_point())then
!-zero           write(*,*)'###### land_cover_cold_start: lake fraction is set to zero #####'
!-zero        endif
!-zero     endif
!-zero  enddo
!-zero  enddo
  
  ! reconcile vegetation fractions with the land mask within compute domain
  valid_data = sum(vegn,3) > 0
  invalid_data = .FALSE.
  do j = 1,size(land_mask,2)
  do i = 1,size(land_mask,1)
     if(.not.land_mask(i,j)) cycle ! skip ocean points
     if(valid_data(i,j)) cycle ! don't need to do anything with valid points
     if(sum(glac(i,j,:))+sum(lake(i,j,:))>=1) &
          cycle                ! skip points fully covered by glaciers or lakes
     invalid_data(i,j)=.TRUE.
  enddo
  enddo
  call horiz_remap_new(invalid_data,valid_data,lnd%lon,lnd%lat,lnd%domain,&
       lnd%pelist(ps:pe),map)
  if (print_remapping) call horiz_remap_print(map,'vegetation cover remap:')
  call horiz_remap(map,lnd%domain,vegn)
  call horiz_remap_del(map)
  
  ! create tiles
  do j = 1,size(land_mask,2)
  do i = 1,size(land_mask,1)
     if(.not.land_mask(i,j)) cycle ! skip ocean points
     call set_current_point(i+lnd%is-1,j+lnd%js-1,1)
     call land_cover_cold_start_0d &
          (lnd%tile_map(i+lnd%is-1,j+lnd%js-1),glac(i,j,:),lake(i,j,:),soil(i,j,:),vegn(i,j,:))
     if(nitems(lnd%tile_map(i+lnd%is-1,j+lnd%js-1))==0) then
        call error_mesg('land_cover_cold_start',&
             'No tiles were created for a valid land point at i='&
             //trim(string(lnd%is+i-1))//' j='//trim(string(lnd%js+j-1))//' face='//trim(string(lnd%face)), FATAL)
     endif
  enddo
  enddo

  deallocate(glac,lake,soil,vegn)
  
end subroutine land_cover_cold_start

! ============================================================================
subroutine land_cover_cold_start_0d (set,glac0,lake0,soil0,vegn0)
  type(land_tile_list_type), intent(inout) :: set 
  real, dimension(:)       , intent(in) :: &
       glac0,lake0,soil0,vegn0 ! fractions of area

  ! ---- local vars
  real :: glac(size(glac0(:))), lake(size(lake0(:))), &
          soil(size(soil0(:))), vegn(size(vegn0(:)))
  type(land_tile_type), pointer :: tile
  integer :: i,j,k
  real :: factor ! normalizing factor for the tile areas
  real :: frac
  type(land_tile_enum_type) :: first_non_vegn ! position of first non-vegetated tile in the list

  glac = glac0; lake = lake0; soil = soil0; vegn = vegn0
  if (sum(glac)>1) &
       glac=glac/sum(glac)
  if (sum(lake)+sum(glac)>1)&
       lake = lake*(1-sum(glac))/sum(lake)
  if (sum(lake)<min_sum_lake_frac) lake=0
  if (sum(soil)+sum(glac)+sum(lake)>1)&
       soil = soil*(1-sum(lake)-sum(glac))/sum(soil)
  ! make sure that the sum of the fractions of the soil, lake, and glaciers are 
  ! either one or zero
  factor = sum(soil)+sum(glac)+sum(lake)
  if(factor>0)then
     glac = glac/factor
     lake = lake/factor
     soil = soil/factor
  endif

  ! remove soil/glac/lake fractions that are too small
  if (min_frac>0) then
     where (glac<min_frac) glac = 0
     where (lake<min_frac) lake = 0
     where (soil<min_frac) soil = 0
     ! do the renormalization again
     factor = sum(soil)+sum(glac)+sum(lake)
     if(factor>0)then
	glac = glac/factor
	lake = lake/factor
	soil = soil/factor
     endif
  endif

  if(is_watch_point()) then
     write(*,*)'#### land_cover_cold_start_0d input data ####'
     __DEBUG1__(glac0)
     __DEBUG1__(lake0)
     __DEBUG1__(soil0)
     __DEBUG1__(vegn0)
     __DEBUG1__(factor)
     write(*,*)'#### land_cover_cold_start_0d renormlaized fractions ####'
     __DEBUG1__(glac)
     __DEBUG1__(lake)
     __DEBUG1__(soil)
     __DEBUG1__(vegn)
  endif

  do i = 1,size(glac)
     if (glac(i)>0) then
        tile => new_land_tile(frac=glac(i),glac=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created glac tile: frac=',glac(i),' tag=',i
        endif
     endif
  enddo
  do i = 1,size(lake)
     if (lake(i)>0) then
        tile => new_land_tile(frac=lake(i),lake=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created lake tile: frac=',lake(i),' tag=',i
        endif
     endif
  enddo

  factor = sum(soil)*sum(vegn)
  if (factor/=0) factor = 1/factor
  factor = factor*(1-sum(glac)-sum(lake))
  ! vegetation tiles, if any, are inserted in front of non-vegetated tiles;
  ! this really doesn't matter except for the static vegetation override
  ! case with the data saved by lm3v -- there the vegetation tiles are
  ! in front, so it works more consistently where lad2 has more than 
  ! one tile (e.g. glac/soil or lake/soil), if lad2 vegetation tiles are 
  ! also in front of the list. 
  first_non_vegn=first_elmt(set)
  do i = 1,size(soil)
  do j = 1,size(vegn)
     frac = soil(i)*vegn(j)*factor
     if(frac>0) then
        tile  => new_land_tile(frac=frac,soil=i,vegn=j)
        call insert(tile,first_non_vegn)
        if(is_watch_point()) then
           write(*,*)'created soil tile: frac=', frac, ' soil tag=',i, ' veg tag=',j
        endif
     endif
  enddo
  enddo

end subroutine land_cover_cold_start_0d

! ============================================================================
! reads the land restart file and restores the tiling structure from this file
subroutine land_cover_warm_start ( restart_file_name, lnd )
  character(len=*), intent(in) :: restart_file_name
  type(land_state_type), intent(inout) :: lnd
    
  ! ---- local vars
  integer, allocatable :: idx(:) ! compressed tile index
  integer, allocatable :: glac(:), lake(:), soil(:), snow(:), cana(:), vegn(:) ! tile tags
  real,    allocatable :: frac(:) ! fraction of land covered by tile
  integer :: ncid ! unit number of the input file
  integer :: ntiles    ! total number of land tiles in the input file
  integer :: bufsize   ! size of the input buffer
  integer :: dimids(1) ! id of tile dimension
  character(NF_MAX_NAME) :: tile_dim_name ! name of the tile dimension and respective variable
  integer :: i,j,k,it
  type(land_tile_type), pointer :: tile;
  integer :: start, count ! slab for reading
  ! netcdf variable IDs
  integer :: id_idx, id_frac, id_glac, id_lake, id_soil, id_vegn
  
  __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,ncid))
  ! allocate the input data
  __NF_ASRT__(nfu_inq_var(ncid,'frac',id=id_frac,varsize=ntiles,dimids=dimids))
   ! allocate input buffers for compression index and the variable
  bufsize=min(input_buf_size,ntiles)
  allocate(idx (bufsize), glac(bufsize), lake(bufsize), soil(bufsize), &
           snow(bufsize), cana(bufsize), vegn(bufsize), frac(bufsize)  )
  ! get the name of the fist (and only) dimension of the variable 'frac' -- this
  ! is supposed to be the compressed dimension, and associated variable will
  ! hold the compressed indices
  __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),name=tile_dim_name))
  __NF_ASRT__(nfu_inq_var(ncid,tile_dim_name,id=id_idx))
  ! get the IDs of the varables to read
  __NF_ASRT__(nfu_inq_var(ncid,'glac',id=id_glac))
  __NF_ASRT__(nfu_inq_var(ncid,'lake',id=id_lake))
  __NF_ASRT__(nfu_inq_var(ncid,'soil',id=id_soil))
  __NF_ASRT__(nfu_inq_var(ncid,'vegn',id=id_vegn))
  
  do start = 1,ntiles,bufsize
    count = min(bufsize,ntiles-start+1)
    ! read the compressed tile indices
    __NF_ASRT__(nf_get_vara_int(ncid,id_idx,(/start/),(/count/),idx))
    ! read input data -- fractions and tags
    __NF_ASRT__(nf_get_vara_double(ncid,id_frac,(/start/),(/count/),frac))
    __NF_ASRT__(nf_get_vara_int(ncid,id_glac,(/start/),(/count/),glac))
    __NF_ASRT__(nf_get_vara_int(ncid,id_lake,(/start/),(/count/),lake))
    __NF_ASRT__(nf_get_vara_int(ncid,id_soil,(/start/),(/count/),soil))
    __NF_ASRT__(nf_get_vara_int(ncid,id_vegn,(/start/),(/count/),vegn))
  
    ! create tiles
    do it = 1,count
       k = idx(it)
       if (k<0) cycle ! skip negative indices
       i = modulo(k,lnd%nlon)+1; k = k/lnd%nlon
       j = modulo(k,lnd%nlat)+1; k = k/lnd%nlat
       k = k + 1
       if (i<lnd%is.or.i>lnd%ie) cycle
       if (j<lnd%js.or.j>lnd%je) cycle
       ! the size of the tile set at the point (i,j) must be equal to k
       tile=>new_land_tile(frac=frac(it),&
                glac=glac(it),lake=lake(it),soil=soil(it),vegn=vegn(it))
       call insert(tile,lnd%tile_map(i,j))
    enddo
  enddo
  __NF_ASRT__(nf_close(ncid))
  deallocate(idx, glac, lake, soil, snow, cana, vegn, frac)
end subroutine


! ============================================================================
subroutine update_land_model_fast ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars 
  real :: fco2_0,Dfco2Dq , & ! co2 flux from canopy air to the atmosphere
       ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
       ISa_dn_dif(NBANDS), & ! downward diffuse sw radiation at the top of the canopy
       cana_q

  ! variables for stock calculations
  real :: &
     cana_VMASS, cana_HEAT,             &
     vegn_LMASS, vegn_FMASS, vegn_HEAT, &
     snow_LMASS, snow_FMASS, snow_HEAT, &
     subs_LMASS, subs_FMASS, subs_HEAT, &
     glac_LMASS, glac_FMASS, glac_HEAT, &
     lake_LMASS, lake_FMASS, lake_HEAT, &
     soil_LMASS, soil_FMASS, soil_HEAT 

  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je) :: &
       runoff,           & ! total (liquid+snow) runoff accumulated over tiles in cell
       runoff_snow,      & ! runoff snow accumulated over tiles in cell
       runoff_heat,      & ! runoff heat accumulated over tiles in cell
       heat_frac_liq,    & ! fraction of runoff heat in liquid
       discharge_l,      & ! discharge of liquid water to ocean
       discharge_sink      ! container to collect small/negative values for later accounting
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je,num_species) :: &
       runoff_c,         & ! runoff of tracers accumulated over tiles in cell
       discharge_c         ! discharge of tracers to ocean
  logical :: used          ! return value of send_data diagnostics routine
  real, allocatable :: runoff_1d(:),runoff_snow_1d(:),runoff_heat_1d(:)
  integer :: i,j,k     ! lon, lat, and tile indices
  integer :: i_species ! river tracer iterator
  integer :: i1        ! index used to iterate over grid cells efficiently
  integer :: is,ie,js,je ! horizontal bounds of the override buffer
  type(land_tile_enum_type) :: ce, te ! tile enumarator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  
  ! variables for data override
  real, allocatable :: phot_co2_data(:,:)  ! buffer for data
  logical           :: phot_co2_overridden ! flag indicating successful override
  

  ! start clocks
  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landFastClock)

  ! to avoid output of static vegetation after the transitions worked and
  ! changed the tiling structure, static vegetation output is done here.
  call write_static_vegn()

  ! override data at the beginning of the time step
  is=lbound(cplr2land%t_flux,1) ; ie = is+size(cplr2land%t_flux,1)-1
  js=lbound(cplr2land%t_flux,2) ; je = js+size(cplr2land%t_flux,2)-1
  allocate(phot_co2_data(is:ie,js:je))
  call data_override('LND','phot_co2',phot_co2_data,lnd%time, &
       override=phot_co2_overridden)
  
  ! clear the runoff values, for accumulation over the tiles
  runoff = 0 ; runoff_snow = 0 ; runoff_heat = 0  ; runoff_c = 0

  ! main tile loop
!$OMP parallel do schedule(dynamic) default(shared) private(i1,i,j,k,ce,te,tile,fco2_0,Dfco2Dq,ISa_dn_dir,ISa_dn_dif)
  do i1 = 0,(ie-is+1)*(je-js+1)-1
     i = mod(i1,ie-is+1)+is
     j = i1/(ie-is+1)+js
!     __DEBUG4__(is,js,i-is+lnd%is,j-js+lnd%js)
     ce = first_elmt(lnd%tile_map(i-is+lnd%is,j-js+lnd%js))
     te = tail_elmt (lnd%tile_map(i-is+lnd%is,j-js+lnd%js))
     k = 0 
     do while (ce/=te)
        k = k+1 ; tile=>current_tile(ce) ; ce = next_elmt(ce)
   
        ! set this point coordinates as current for debug output
        call set_current_point(i-is+lnd%is,j-js+lnd%js,k)
   
        if (lnd%ico2/=NO_TRACER) then
           fco2_0  = cplr2land%tr_flux(i,j,k, lnd%ico2)
           Dfco2Dq = cplr2land%dfdtr  (i,j,k, lnd%ico2)
        else
           fco2_0  = 0
           Dfco2Dq = 0
        endif
        ISa_dn_dir(BAND_VIS) = cplr2land%sw_flux_down_vis_dir(i,j,k)
        ISa_dn_dir(BAND_NIR) = cplr2land%sw_flux_down_total_dir(i,j,k)&
                              -cplr2land%sw_flux_down_vis_dir(i,j,k)
        ISa_dn_dif(BAND_VIS) = cplr2land%sw_flux_down_vis_dif(i,j,k)
        ISa_dn_dif(BAND_NIR) = cplr2land%sw_flux_down_total_dif(i,j,k)&
                              -cplr2land%sw_flux_down_vis_dif(i,j,k)
   
        call update_land_model_fast_0d(tile, i,j,k, land2cplr, &
           cplr2land%lprec(i,j,k),  cplr2land%fprec(i,j,k), cplr2land%tprec(i,j,k), &
           cplr2land%t_flux(i,j,k), cplr2land%dhdt(i,j,k), &
           cplr2land%tr_flux(i,j,k, lnd%isphum), cplr2land%dfdtr(i,j,k, lnd%isphum), &
           fco2_0, Dfco2Dq, &
           ISa_dn_dir, ISa_dn_dif, cplr2land%lwdn_flux(i,j,k), &
           cplr2land%ustar(i,j,k), cplr2land%p_surf(i,j,k), cplr2land%drag_q(i,j,k), &
           phot_co2_overridden, phot_co2_data(i,j),&
           runoff(i,j), runoff_heat(i,j), runoff_snow(i,j) &
        )
        ! some of the diagnostic variables are sent from here, purely for coding 
        ! convenince: the compute domain-level 2d and 3d vars are generaly not 
        ! available inside update_land_model_fast_0d, so the diagnostics for those 
        ! was left here.
        call send_tile_data(id_area, tile%frac*lnd%area(i,j),        tile%diag)
        call send_tile_data(id_z0m,  land2cplr%rough_mom(i,j,k),     tile%diag)
        call send_tile_data(id_z0s,  land2cplr%rough_heat(i,j,k),    tile%diag)
        call send_tile_data(id_Trad, land2cplr%t_surf(i,j,k),        tile%diag)
        call send_tile_data(id_Tca,  land2cplr%t_ca(i,j,k),          tile%diag)
        call send_tile_data(id_qca,  land2cplr%tr(i,j,k,lnd%isphum), tile%diag)
	call send_tile_data(id_cd_m, cplr2land%cd_m(i,j,k),          tile%diag)
	call send_tile_data(id_cd_t, cplr2land%cd_t(i,j,k),          tile%diag)
     enddo
  enddo
  
  ! set values of tracer fluxes
  runoff_c(:,:,1) = runoff_snow
  runoff_c(:,:,2) = runoff_heat
  do i_species = num_phys+1, num_species
    runoff_c(:,:,i_species) = 0        ! age, species
    enddo

!=================================================================================
  ! update river state
  call update_river(runoff, runoff_c, discharge_l, discharge_c)
!=================================================================================

  discharge_l = discharge_l/lnd%cellarea
  do i_species = 1, num_species
    discharge_c(:,:,i_species) =  discharge_c(:,:,i_species)/lnd%cellarea
    enddo

  ! pass through to ocean the runoff that was not seen by river module because of land_frac diffs.
  ! need to multiply by gfrac to spread over whole cell
  where (missing_rivers) discharge_l = (runoff-runoff_c(:,:,1))*frac
  do i_species = 1, num_species
    where (missing_rivers) &
     discharge_c(:,:,i_species) = runoff_c(:,:,i_species)*frac
  enddo

  ! don't send negatives or insignificant values to ocean. put them in the sink instead.
  ! this code does not seem necessary, and default discharge_tol value should be used.
  discharge_sink = 0.
  where (discharge_l.le.discharge_tol)
      discharge_sink = discharge_sink + discharge_l
      discharge_l    = 0.
    endwhere
  where (discharge_c(:,:,1).le.discharge_tol)
      discharge_sink     = discharge_sink + discharge_c(:,:,1)
      discharge_c(:,:,1) = 0.
    endwhere

  ! find phase partitioning ratio for discharge sensible heat flux
  where (discharge_l.gt.0. .or. discharge_c(:,:,1).gt.0.)
      heat_frac_liq = clw*discharge_l / (clw*discharge_l+csw*discharge_c(:,:,1))
    elsewhere
      heat_frac_liq = 1.
    endwhere

  ! scale up fluxes sent to ocean to compensate for non-ocean fraction of discharge cell.
  ! split heat into liquid and solid streams
  where (frac.lt.1.) 
      land2cplr%discharge           = discharge_l        / (1-frac)
      land2cplr%discharge_snow      = discharge_c(:,:,1) / (1-frac)
      land2cplr%discharge_heat      = heat_frac_liq*discharge_c(:,:,2) / (1-frac)
      land2cplr%discharge_snow_heat =               discharge_c(:,:,2) / (1-frac) &
                                     - land2cplr%discharge_heat
    endwhere

  ce = first_elmt(lnd%tile_map,                  &
              is=lbound(cplr2land%t_flux,1), &
              js=lbound(cplr2land%t_flux,2)  )
  te = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     call get_elmt_indices(ce,i,j,k)
     tile => current_tile(ce)
     ce=next_elmt(ce)
     cana_VMASS = 0. ;                   cana_HEAT = 0.
     vegn_LMASS = 0. ; vegn_FMASS = 0. ; vegn_HEAT = 0.
     snow_LMASS = 0. ; snow_FMASS = 0. ; snow_HEAT = 0.
     subs_LMASS = 0. ; subs_FMASS = 0. ; subs_HEAT = 0.
     glac_LMASS = 0. ; glac_FMASS = 0. ; glac_HEAT = 0.
     lake_LMASS = 0. ; lake_FMASS = 0. ; lake_HEAT = 0.
     soil_LMASS = 0. ; soil_FMASS = 0. ; soil_HEAT = 0.
     if (associated(tile%cana)) then
         call cana_state( tile%cana, cana_q=cana_q )
         cana_VMASS = canopy_air_mass*cana_q
         cana_HEAT  = cana_tile_heat(tile%cana)
       endif
     if (associated(tile%vegn)) then
         call vegn_tile_stock_pe(tile%vegn, vegn_LMASS, vegn_FMASS)
         vegn_HEAT = vegn_tile_heat(tile%vegn)
       endif
     if(associated(tile%snow)) then
         call snow_tile_stock_pe(tile%snow, snow_LMASS, snow_FMASS)
         snow_HEAT = snow_tile_heat(tile%snow)
     endif
     if (associated(tile%glac)) then
         call glac_tile_stock_pe(tile%glac, subs_LMASS, subs_FMASS)
         subs_HEAT  = glac_tile_heat(tile%glac)
         glac_LMASS = subs_LMASS
         glac_FMASS = subs_FMASS
         glac_HEAT  = subs_HEAT
       else if (associated(tile%lake)) then
         call lake_tile_stock_pe(tile%lake, subs_LMASS, subs_FMASS)
         subs_HEAT  = lake_tile_heat(tile%lake)
         lake_LMASS = subs_LMASS
         lake_FMASS = subs_FMASS
         lake_HEAT  = subs_HEAT
       else if (associated(tile%soil)) then
         call soil_tile_stock_pe(tile%soil, subs_LMASS, subs_FMASS)
         subs_HEAT  = soil_tile_heat(tile%soil)
         soil_LMASS = subs_LMASS
         soil_FMASS = subs_FMASS
         soil_HEAT  = subs_HEAT
       endif

     call send_tile_data(id_VWS,  cana_VMASS, tile%diag)
     call send_tile_data(id_VWSc, cana_VMASS, tile%diag)
     call send_tile_data(id_LWS,  vegn_LMASS+snow_LMASS+subs_LMASS, tile%diag)
     call send_tile_data(id_LWSv, vegn_LMASS, tile%diag)
     call send_tile_data(id_LWSs, snow_LMASS, tile%diag)
     call send_tile_data(id_LWSg, subs_LMASS, tile%diag)
     call send_tile_data(id_FWS,  vegn_FMASS+snow_FMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_FWSv, vegn_FMASS, tile%diag)
     call send_tile_data(id_FWSs, snow_FMASS, tile%diag)
     call send_tile_data(id_FWSg, subs_FMASS, tile%diag)
     call send_tile_data(id_HS,  vegn_HEAT+snow_HEAT+subs_HEAT+cana_HEAT, tile%diag)
     call send_tile_data(id_HSv, vegn_HEAT, tile%diag)
     call send_tile_data(id_HSs, snow_HEAT, tile%diag)
     call send_tile_data(id_HSg, subs_HEAT, tile%diag)
     call send_tile_data(id_HSc, cana_HEAT, tile%diag)
     call send_tile_data(id_water, subs_LMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_snow,  snow_LMASS+snow_FMASS, tile%diag)
     enddo

  ! advance land model time
  lnd%time = lnd%time + lnd%dt_fast

  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(lnd%tile_map, lnd%time)

  if (id_dis_liq > 0)  used = send_data (id_dis_liq,  discharge_l,        lnd%time) 
  if (id_dis_ice > 0)  used = send_data (id_dis_ice,  discharge_c(:,:,1), lnd%time) 
  if (id_dis_heat > 0) used = send_data (id_dis_heat, discharge_c(:,:,2), lnd%time) 
  if (id_dis_sink > 0) used = send_data (id_dis_sink, discharge_sink,     lnd%time) 

  ! deallocate override buffer
  deallocate(phot_co2_data)

  call mpp_clock_end(landFastClock)
  call mpp_clock_end(landClock)
end subroutine update_land_model_fast


! ============================================================================
subroutine update_land_model_fast_0d(tile, i,j,k, land2cplr, &
   precip_l, precip_s, atmos_T, &
   Ha0, DHaDTc, Ea0, DEaDqc, fco2_0, Dfco2Dq,&
   ISa_dn_dir, ISa_dn_dif, ILa_dn, &
   ustar, p_surf, drag_q, &
   phot_co2_overridden, phot_co2_data, &
   runoff, runoff_heat, runoff_snow &
   )
  type (land_tile_type), pointer :: tile
  type(land_data_type), intent(inout) :: land2cplr
  integer, intent(in) :: i,j,k ! coordinates
  real, intent(in) :: &
       precip_l, precip_s, & ! liquid and solid precipitation, kg/(m2 s)
       atmos_T, &        ! incoming precipitation temperature (despite its name), deg K
       Ha0,   DHaDTc, &  ! sensible heat flux from the canopy air to the atmosphere 
       Ea0,   DEaDqc, &  ! water vapor flux from canopy air to the atmosphere
       fco2_0,Dfco2Dq,&  ! co2 flux from canopy air to the atmosphere
       ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
       ISa_dn_dif(NBANDS), & ! downward diffuse sw radiation at the top of the canopy
       ILa_dn,             & ! downward lw radiation at the top of the canopy
       ustar,              & ! friction velocity, m/s
       p_surf,             & ! surface pressure, Pa
       drag_q,             & !
       phot_co2_data         ! data input for the CO2 for photosynthesis

  logical, intent(in):: phot_co2_overridden
  real, intent(inout) :: &
        runoff, runoff_heat, runoff_snow

  ! ---- local constants
  ! indices of variables and equations for implicit time stepping solution :
  integer, parameter :: iqc=1, iTc=2, iTv=3, iwl=4, iwf=5

  ! ---- local vars 
  real :: A(5,5),B0(5),B1(5),B2(5) ! implicit equation matrix and right-hand side vectors
  real :: A00(5,5),B10(5),B00(5) ! copy of the above, only for debugging
  integer :: indx(5) ! permutation vector
  ! linearization coefficients of various fluxes between components of land
  ! surface scheme
  real :: &
       G0,    DGDTg,  &  ! ground heat flux 
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux from vegetation
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf, & ! sublimation of intercepted snow
       Hg0,   DHgDTg,   DHgDTc, & ! linearization of the sensible heat flux from ground
       Eg0,   DEgDTg,   DEgDqc, DEgDpsig, & ! linearization of evaporation from ground
       flwv0,  DflwvDTg,  DflwvDTv,& ! linearization of net LW radiation to the canopy
       flwg0,  DflwgDTg,  DflwgDTv,& ! linearization of net LW radiation to the canopy
       vegn_drip_l, vegn_drip_s, & ! drip rate of water and snow, respectively, kg/(m2 s)
       vegn_lai
  
  ! increments of respective variables over time step, results of the implicit
  ! time step:
  real :: delta_qc, delta_Tc, delta_Tv, delta_wl, delta_ws, delta_Tg, delta_psig
  real :: flwg ! updated value of long-wave ground energy balance
  real :: vegn_emis_lw, surf_emis_lw ! emissivities of ground and surface
  real :: vegn_emsn,    surf_emsn    ! emission by vegetation and surface, respectively
  real :: denom ! denominator in the LW radiative balance calculations
  real :: sum0, sum1

  real :: &
       grnd_T, gT, & ! ground temperature and its value used for sensible heat advection
       vegn_T, vT, & ! vegetation (canopy) temperature
       cana_T, cT, & ! canopy air temperature
       evap_T, eT, & ! temperature assigned to vapor going between land and atmosphere
       soil_uptake_T, & ! average temperature of water taken up by the vegetation
       vegn_Wl,  vegn_Ws, & ! water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid or frozen precipitation
       vegn_hcap,      & ! vegetation heat capacity, including intercepted water and snow
       vegn_fco2, & ! co2 flux from the vegetation, kg CO2/(m2 s)
       hlv_Tv, hlv_Tu, & ! latent heat of vaporization at vegn and uptake temperatures, respectively 
       hls_Tv, &         ! latent heat of sublimation at vegn temperature
       grnd_rh,        & ! explicit relative humidity at ground surface
       grnd_rh_psi,    & ! psi derivative of relative humidity at ground surface
       grnd_liq, grnd_ice, grnd_subl, &
       grnd_tf, &  ! temperature of freezing on the ground
       grnd_latent, &
       grnd_flux, &
       grnd_E_min, &
       grnd_E_max, &
       soil_E_min, &
       soil_E_max, &
       soil_beta, &
       RSv(NBANDS), & ! net short-wave radiation balance of the canopy, W/m2
       con_g_h, con_g_v, & ! turbulent cond. between ground and canopy air, for heat and vapor respectively
       snow_area, &
       cana_q, & ! specific humidity of canopy air
       cana_co2, & ! co2 moist mixing ratio in canopy air, kg CO2/kg wet air
       cana_co2_mol, & ! co2 dry mixing ratio in canopy air, mol CO2/mol dry air
       fswg, evapg, sensg, &
       subs_G, subs_G2, Mg_imp, snow_G_Z, snow_G_TZ, &
       snow_avrg_T, delta_T_snow,  & ! vertically-average snow temperature and it's change due to s
       vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs, & ! heat flux from canopy due to overflow
       delta_fprec, & ! correction of below-canopy solid precip in case it's average T > tfreeze 

       hprec,              & ! sensible heat flux carried by precipitation
       hevap,              & ! sensible heat flux carried by total evapotranspiration
       land_evap,          & ! total vapor flux from land to atmosphere
       land_sens,          & ! turbulent sensible heat flux from land to atmosphere
       vegn_flw,vegn_sens,snow_sens,snow_levap,snow_fevap,snow_melt,&
       snow_lprec, snow_hlprec,snow_lrunf,vegn_levap,vegn_fevap,vegn_uptk,&
       vegn_fsw, vegn_melt,vegn_lprec,vegn_fprec,vegn_hlprec,vegn_hfprec,&
       precip_T,pT,snow_fsw,snow_flw,snow_frunf,snow_hlrunf,&
       snow_hfrunf,subs_fsw,subs_flw,subs_sens,&
       subs_DT, subs_M_imp, subs_evap, snow_Tbot, snow_Cbot, snow_C, subs_levap,&
       subs_fevap,subs_melt,subs_lrunf,subs_hlrunf,&
       subs_Ttop,subs_Ctop, subs_subl, new_T
  real :: soil_water_supply ! supply of water to roots, per unit active root biomass, kg/m2
  real :: snow_T, snow_rh, snow_liq, snow_ice, snow_subl
  integer :: ii, jj ! indices for debug output
  integer :: ierr
  logical :: conserve_glacier_mass, snow_active, redo_leaf_water
  integer :: canopy_water_step
  real :: subs_z0m, subs_z0s, snow_z0m, snow_z0s, grnd_z0s

  soil_uptake_T = tfreeze ! just to avoid using un-initialized values
  soil_water_supply = 0.0
  if (associated(tile%glac)) then
     call glac_step_1 ( tile%glac, &
          grnd_T, grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
          snow_G_Z, snow_G_TZ, conserve_glacier_mass  )
     grnd_E_min = -HUGE(grnd_E_min)
     grnd_E_max =  HUGE(grnd_E_max)
     grnd_rh_psi = 0
  else if (associated(tile%lake)) then
     call lake_step_1 ( ustar, p_surf, &
          lnd%lat(i,j), tile%lake, &
          grnd_T, grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
          snow_G_Z, snow_G_TZ)
     grnd_E_min = -HUGE(grnd_E_min)
     grnd_E_max =  HUGE(grnd_E_max)
     grnd_rh_psi = 0
  else if (associated(tile%soil)) then
     call soil_step_1 ( tile%soil, tile%vegn, tile%diag, &
          grnd_T, soil_uptake_T, soil_beta, soil_water_supply, soil_E_min, soil_E_max, &
          grnd_rh, grnd_rh_psi, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
          snow_G_Z, snow_G_TZ)
     grnd_E_min = soil_E_min
     grnd_E_max = soil_E_max
     grnd_liq = 0 ! sorry, but solver cannot handle implicit melt anymore
     grnd_ice = 0 ! sorry, but solver cannot handle implicit melt anymore
                  ! no big loss, it's just the surface layer anyway
  else
     call get_current_point(face=ii)
     call error_mesg('update_land_model_fast','none of the surface tiles exist at ('//&
          trim(string(i))//','//trim(string(j))//','//trim(string(k))//&
          ', face='//trim(string(ii))//')',FATAL)
  endif

  subs_subl = grnd_subl

  call snow_step_1 ( tile%snow, snow_G_Z, snow_G_TZ, &
       snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, G0, DGDTg )
  if (snow_active) then
     grnd_T    = snow_T;   grnd_rh   = snow_rh;   grnd_liq  = snow_liq
     grnd_rh_psi = 0
     grnd_ice  = snow_ice; grnd_subl = snow_subl; grnd_tf   = tfreeze
     grnd_E_min = -HUGE(grnd_E_min)
     grnd_E_max =  HUGE(grnd_E_max)
  endif

  call cana_state(tile%cana, cana_T, cana_q, cana_co2)

  if (associated(tile%vegn)) then
  ! Calculate net short-wave radiation input to the vegetation
     RSv    = tile%Sv_dir*ISa_dn_dir + tile%Sv_dif*ISa_dn_dif
     call soil_diffusion(tile%soil, subs_z0s, subs_z0m)
     call snow_diffusion(tile%snow, snow_z0s, snow_z0m)
     grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))
     
     ! cana_co2 is moist mass mixing ratio [kg CO2/kg wet air], convert it to dry
     ! volumetric mixing ratio [mol CO2/mol dry air] 
     cana_co2_mol = cana_co2*mol_air/mol_CO2/(1-cana_q)
     if (phot_co2_overridden) cana_co2_mol = phot_co2_data
     call vegn_step_1 ( tile%vegn, tile%soil, tile%diag, &
        p_surf, &
        ustar, &
        drag_q, &
        ISa_dn_dir+ISa_dn_dif, RSv, precip_l, precip_s, &
        tile%land_d, tile%land_z0s, tile%land_z0m, grnd_z0s, & 
        soil_beta, soil_water_supply,&
        cana_T, cana_q, cana_co2_mol, &
        ! output
        con_g_h, con_g_v, &
        vegn_T, vegn_Wl, vegn_Ws, & ! temperature, water and snow mass on the canopy
        vegn_ifrac, vegn_lai, &
        vegn_drip_l, vegn_drip_s,& 
        vegn_hcap, & ! total vegetation heat capacity (including intercepted water/snow)
        Hv0,   DHvDTv,   DHvDTc,            & 
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & 
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & 
        Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf  ) 
	if (LM2) then
	   con_g_h = con_g_h * con_fac_large
	   if (snow_active) then
	      con_g_v = con_g_v * con_fac_large
	   else
              con_g_v = con_g_v * con_fac_small
	   endif
	endif
  else
     RSv    = 0
     con_g_h = con_fac_large ; con_g_v = con_fac_large
     if(associated(tile%glac).and.conserve_glacier_mass.and..not.snow_active) &
          con_g_v = con_fac_small
     vegn_T  = cana_T ; vegn_Wl = 0 ; vegn_Ws = 0
     vegn_ifrac  = 0 ; vegn_lai    = 0
     vegn_drip_l = 0 ; vegn_drip_s = 0
     vegn_hcap = 1.0
     Hv0 =0;  DHvDTv =0;  DHvDTc=0;
     Et0 =0;  DEtDTv =0;  DEtDqc=0;   DEtDwl=0;   DEtDwf=0
     Eli0=0;  DEliDTv=0;  DEliDqc=0;  DEliDwl=0;  DEliDwf=0 
     Esi0=0;  DEsiDTv=0;  DEsiDqc=0;  DEsiDwl=0;  DEsiDwf=0
  endif
  ! calculate net shortwave for ground and canopy
  fswg     = SUM(tile%Sg_dir*ISa_dn_dir + tile%Sg_dif*ISa_dn_dif)
  vegn_fsw = SUM(RSv)
  
  call cana_step_1 (tile%cana, p_surf, con_g_h, con_g_v,   &
       grnd_t, grnd_rh, grnd_rh_psi, &
       Hg0,  DHgDTg, DHgDTc, Eg0, DEgDTg, DEgDqc, DEgDpsig)

! [X.X] using long-wave optical properties, calculate the explicit long-wave 
!       radiative balances and their derivatives w.r.t. temperatures
  vegn_emis_lw = 1 - tile%vegn_refl_lw - tile%vegn_tran_lw
  surf_emis_lw = 1 - tile%surf_refl_lw

  denom = 1-tile%vegn_refl_lw*tile%surf_refl_lw

  vegn_emsn = vegn_emis_lw * stefan * vegn_T**4
  surf_emsn = surf_emis_lw * stefan * grnd_T**4

  flwv0 = ILa_dn * vegn_emis_lw*(1+tile%vegn_refl_lw*tile%surf_refl_lw/denom) &
       + vegn_emsn * (tile%surf_refl_lw*vegn_emis_lw/denom-2) &
       + surf_emsn * vegn_emis_lw/denom
  DflwvDTg = vegn_emis_lw/denom                      * surf_emis_lw * stefan * 4 * grnd_T**3
  DflwvDTv = (tile%surf_refl_lw*vegn_emis_lw/denom-2)* vegn_emis_lw * stefan * 4 * vegn_T**3

  flwg0 = (ILa_dn*tile%vegn_tran_lw + vegn_emsn)*(1-tile%surf_refl_lw)/denom &
       - surf_emsn*(1-tile%vegn_refl_lw)/denom
  DflwgDTg = -(1-tile%vegn_refl_lw)/denom * surf_emis_lw * stefan * 4 * grnd_T**3
  DflwgDTv =  (1-tile%surf_refl_lw)/denom * vegn_emis_lw * stefan * 4 * vegn_T**3

! [X.0] calculate the latent heats of vaporization at appropriate temperatures
  if (use_tfreeze_in_grnd_latent) then
    grnd_latent = hlv + hlf*grnd_subl
  else
    grnd_latent = hlv + (cpw-clw)*(grnd_T-tfreeze) &
               + (hlf + (clw-csw)*(grnd_T-tfreeze)) * grnd_subl
  endif
  if (use_atmos_T_for_precip_T) then
    precip_T = atmos_T
  else
    precip_T = cana_T
  endif
  if (use_atmos_T_for_evap_T) then
    evap_T = atmos_T
  else
    evap_T = cana_T
  endif
  if (use_old_conservation_equations) then
    hlv_Tv = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T
    hls_Tv = hlv + hlf - (cpw-csw)*tfreeze + cpw*vegn_T
    hlv_Tu = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T - clw*soil_uptake_T
    pT = precip_T
    cT = cana_T
    eT = evap_T
    gT = grnd_T
    vT = vegn_T
  else
    hlv_Tv = hlv    + cpw*(vegn_T-tfreeze)
    hls_Tv = hlf    + hlv_Tv
    hlv_Tu = hlv_Tv - clw*(soil_uptake_T-tfreeze)
    pT = precip_T-tfreeze
    cT = cana_T-tfreeze
    eT = evap_T-tfreeze
    gT = grnd_T-tfreeze
    vT = vegn_T-tfreeze
  endif

  do canopy_water_step = 1,2
     if(is_watch_point()) then
        write(*,*)'#### input data for the matrix ####'
        __DEBUG1__(delta_time)
        __DEBUG4__(vegn_T,vT,vegn_Wl,vegn_Ws)
        __DEBUG3__(grnd_T,gT,grnd_rh)
        __DEBUG3__(cana_T,cT,cana_q)
        __DEBUG2__(evap_T,eT)
        __DEBUG2__(vegn_emis_lw,surf_emis_lw)
        __DEBUG2__(vegn_emsn,surf_emsn)
        __DEBUG4__(precip_l, vegn_drip_l, pT, precip_T)
        __DEBUG2__(precip_s, vegn_drip_s)
        __DEBUG2__(vegn_ifrac, vegn_lai)
        __DEBUG1__(ILa_dn)
        __DEBUG2__(ISa_dn_dir(1),ISa_dn_dir(2))
        __DEBUG2__(ISa_dn_dif(1),ISa_dn_dif(2))
        __DEBUG2__(fswg, vegn_fsw)
        __DEBUG1__(vegn_hcap)
        __DEBUG3__(hlv_Tv, hlv_Tu, hls_Tv)
        __DEBUG2__(G0, DGDTg)
        __DEBUG2__(Ha0, DHaDTc)
        __DEBUG2__(Ea0, DEaDqc)
        __DEBUG3__(Hv0, DHvDTv, DHvDTc)
        __DEBUG5__(Et0,  DEtDTv,  DEtDqc,  DEtDwl,  DEtDwf)
        __DEBUG5__(Eli0, DEliDTv, DEliDqc, DEliDwl, DEliDwf)
        __DEBUG5__(Esi0, DEsiDTv, DEsiDqc, DEsiDwl, DEsiDwf)
        __DEBUG3__(Hg0, DHgDTg, DHgDTc)
        __DEBUG3__(Eg0, DEgDTg, DEgDqc)
        __DEBUG3__(flwv0, DflwvDTg, DflwvDTv)
        __DEBUG3__(flwg0, DflwgDTg, DflwgDTv)
        __DEBUG2__(tile%e_res_1,tile%e_res_2)
     endif

! [X.1] form the system of equations for implicit scheme, such that A*X = B1*delta_Tg+B2*delta_psig+B0
! [X.1.1] equation of canopy air mass balance
     A(iqc,iqc) = canopy_air_mass/delta_time-DEtDqc-DEliDqc-DEsiDqc-DEgDqc+DEaDqc
     A(iqc,iTc) = 0
     A(iqc,iTv) = -DEtDTv-DEliDTv-DEsiDTv
     A(iqc,iwl) = -DEtDwl-DEliDwl-DEsiDwl
     A(iqc,iwf) = -DEtDwf-DEliDwf-DEsiDwf
     B0(iqc)  = Esi0+Eli0+Et0+Eg0-Ea0
     B1(iqc)  = DEgDTg
     B2(iqc)  = DEgDpsig
! [X.1.2] equation of canopy air energy balance
#ifdef USE_DRY_CANA_MASS
     A(iTc,iqc) = canopy_air_mass*cpw*cT/delta_time &
#else
     A(iTc,iqc) = canopy_air_mass*(cpw-cp_air)*cT/delta_time &
#endif
       - cpw*vT*(DEtDqc+DEliDqc+DEsiDqc) - cpw*gT*DEgDqc + cpw*eT*DEaDqc
#ifdef USE_DRY_CANA_MASS
     A(iTc,iTc) = canopy_air_mass*cp_air/delta_time-DHvDTc-DHgDTc+DHaDTc
#else
     A(iTc,iTc) = canopy_air_mass*(cp_air+cana_q*(cpw-cp_air))/delta_time-DHvDTc-DHgDTc+DHaDTc
#endif
     A(iTc,iTv) = -DHvDTv-cpw*vT*(DEtDTv+DEliDTv+DEsiDTv)
     A(iTc,iwl) =        -cpw*vT*(DEtDwl+DEliDwl+DEsiDwl)
     A(iTc,iwf) =        -cpw*vT*(DEtDwf+DEliDwf+DEsiDwf)
     B0(iTc)  = Hv0 + Hg0 - Ha0 + cpw*(vT*(Et0+Eli0+Esi0)+gT*Eg0-eT*Ea0) - tile%e_res_1
     B1(iTc)  = DHgDTg + cpw*gT*DEgDTg
     B2(iTc)  =          cpw*gT*DEgDpsig
! [X.1.3] equation of canopy energy balance
     A(iTv,iqc) = hlv_Tu*DEtDqc + hlv_Tv*DEliDqc + hls_Tv*DEsiDqc
     A(iTv,iTc) = DHvDTc
     A(iTv,iTv) = vegn_hcap/delta_time-DflwvDTv + DHvDTv + &
          hlv_Tu*DEtDTv + hlv_Tv*DEliDTv + hls_Tv*DEsiDTv + clw*vegn_drip_l + csw*vegn_drip_s
     A(iTv,iwl) = clw*vT/delta_time + hlv_Tu*DEtDwl + hlv_Tv*DEliDwl + hls_Tv*DEsiDwl
     A(iTv,iwf) = csw*vT/delta_time + hlv_Tu*DEtDwf + hlv_Tv*DEliDwf + hls_Tv*DEsiDwf
     B0(iTv)  = vegn_fsw + flwv0 - Hv0 - hlv_Tu*Et0 - Hlv_Tv*Eli0 - hls_Tv*Esi0 &
          + clw*precip_l*vegn_ifrac*pT + csw*precip_s*vegn_ifrac*pT &
          - clw*vegn_drip_l*vT - csw*vegn_drip_s*vT - tile%e_res_2
     B1(iTv)  = DflwvDTg
     B2(iTv)  = 0
! [X.1.4] equation of intercepted liquid water mass balance
     A(iwl,iqc) = DEliDqc
     A(iwl,iTc) = 0
     A(iwl,iTv) = DEliDTv
     A(iwl,iwl) = 1.0/delta_time + DEliDwl
     A(iwl,iwf) = DEliDwf
     B0(iwl)  = -Eli0 + precip_l*vegn_ifrac - vegn_drip_l
     B1(iwl)  = 0
     B2(iwl)  = 0
! [X.1.5] equation of intercepted frozen water mass balance
     A(iwf,iqc) = DEsiDqc
     A(iwf,iTc) = 0
     A(iwf,iTv) = DEsiDTv
     A(iwf,iwl) = DEsiDwl
     A(iwf,iwf) = 1.0/delta_time + DEsiDwf
     B0(iwf)  = -Esi0 + precip_s*vegn_ifrac - vegn_drip_s
     B1(iwf)  = 0
     B2(iwf)  = 0
! [X.1.6] if LAI becomes zero (and, therefore, all fluxes from vegetation and their 
! derivatives must be zero too) we get a degenerate case. Still, the drip may be non-zero
! because some water may remain from before leaf drop, and non-zero energy residual can be
! carried over from the previous time step.
! To prevent temperature from going haywire in those cases, we simply replace the equations 
! of canopy energy and mass balance with the following:
! vegn_T + delta_Tv = cana_T + delta_Tc
! delta_Wl = -vegn_drip_l*delta_time
! delta_Ws = -vegn_drip_s*delta_time
! the residual vegn_Wl and vegn_Ws, if any, are taken care of by the overflow calculations 
     if(vegn_hcap==0) then
        ! vegn_T + delta_Tv = cana_T + delta_Tc
        A(iTv,:)   = 0
        A(iTv,iTc) = -1
        A(iTv,iTv) = +1
        B0(iTv) = cana_T - vegn_T
        B1(iTv) = 0
        ! delta_Wl = -vegn_drip_l*delta_time
        A(iwl,:)   = 0
        A(iwl,iwl) = 1
        B0(iwl) = -vegn_drip_l*delta_time
        B1(iwl) = 0
        ! delta_Ws = -vegn_drip_s*delta_time
        A(iwf,:)   = 0
        A(iwf,iwf) = 1
        B0(iwf) = -vegn_drip_s*delta_time
        B1(iwf) = 0
     endif

     if(is_watch_point()) then
        write(*,*)'#### A, B0, B1, B2 ####'
        do ii = 1, size(A,1)
           write(*,'(99g23.16)')(A(ii,jj),jj=1,size(A,2)),B0(ii),B1(ii),B2(ii)
        enddo
     endif

     A00 = A
     B00 = B0
     B10 = B1

! [X.2] solve the system for free terms and delta_Tg and delta_psig terms, getting
!       linear equation for delta_Tg and delta_psig
     call ludcmp(A,indx, ierr)
     if (ierr/=0)&
          write(*,*) 'Matrix is singular',i,j,k
     call lubksb(A,indx,B0)
     call lubksb(A,indx,B1)
     call lubksb(A,indx,B2)

     if(is_watch_point()) then
        write(*,*)'#### solution: B0, B1, B2 ####'
        do ii = 1, size(A,1)
           __DEBUG3__(B0(ii),B1(ii),B2(ii))
        enddo
!!$        write(*,*)'#### solution check ####'
!!$        do ii = 1, size(A,1)
!!$           sum0 = 0; sum1 = 0;
!!$           do jj = 1, size(A,2)
!!$              sum0 = sum0 + A00(ii,jj)*B0(jj)
!!$              sum1 = sum1 + A00(ii,jj)*B1(jj)
!!$           enddo
!!$           write(*,'(99g)')sum0-B00(ii),sum1-B10(ii)
!!$        enddo
     endif
! the result of this solution is a set of expressions for delta_xx in terms
! of delta_Tg and delta_psig: 
! delta_xx(i) = B0(i) + B1(i)*delta_Tg + B2(i)*delta_psig. Note that A, B0, B1 and B2
! are destroyed in the process: A is replaced with LU-decomposition, and
! B0, B1, B2 are replaced with solutions

     ! solve the non-linear equation for energy balance at the surface.

     call land_surface_energy_balance( &
          grnd_T, grnd_liq, grnd_ice, grnd_latent, grnd_Tf, grnd_E_min, &
          grnd_E_max, fswg, &
          flwg0 + b0(iTv)*DflwgDTv, DflwgDTg + b1(iTv)*DflwgDTv, b2(iTv)*DflwgDTv, &
          Hg0   + b0(iTc)*DHgDTc,   DHgDTg   + b1(iTc)*DHgDTc,   b2(iTc)*DHgDTc,   &
          Eg0   + b0(iqc)*DEgDqc,   DEgDTg   + b1(iqc)*DEgDqc,   DEgDpsig + b2(iqc)*DEgDqc,   &
          G0,                       DGDTg, &
          ! output
          delta_Tg, delta_psig, Mg_imp )

! [X.5] calculate final value of other tendencies
     delta_qc = B0(iqc) + B1(iqc)*delta_Tg + B2(iqc)*delta_psig
     delta_Tc = B0(iTc) + B1(iTc)*delta_Tg + B2(iTc)*delta_psig
     delta_Tv = B0(iTv) + B1(iTv)*delta_Tg + B2(iTv)*delta_psig
     delta_wl = B0(iwl) + B1(iwl)*delta_Tg + B2(iwl)*delta_psig
     delta_ws = B0(iwf) + B1(iwf)*delta_Tg + B2(iwf)*delta_psig

! [X.6] calculate updated values of energy balance components used in further 
!       calculations
     flwg       = flwg0 + DflwgDTg*delta_Tg + DflwgDTv*delta_Tv
     evapg      = Eg0   + DEgDTg*delta_Tg   + DEgDpsig*delta_psig + DEgDqc*delta_qc
     sensg      = Hg0   + DHgDTg*delta_Tg   + DHgDTc*delta_Tc
     grnd_flux  = G0    + DGDTg*delta_Tg
     vegn_sens  = Hv0   + DHvDTv*delta_Tv   + DHvDTc*delta_Tc
     vegn_levap = Eli0  + DEliDTv*delta_Tv  + DEliDqc*delta_qc + DEliDwl*delta_wl + DEliDwf*delta_ws
     vegn_fevap = Esi0  + DEsiDTv*delta_Tv  + DEsiDqc*delta_qc + DEsiDwl*delta_wl + DEsiDwf*delta_ws
     vegn_uptk  = Et0   + DEtDTv*delta_Tv   + DEtDqc*delta_qc  + DEtDwl*delta_wl  + DEtDwf*delta_ws
     vegn_flw   = flwv0 + DflwvDTv*delta_Tv + DflwvDTg*delta_Tg
     land_evap  = Ea0   + DEaDqc*delta_qc
     land_sens  = Ha0   + DHaDTc*delta_Tc
! [X.7] calculate energy residuals due to cross-product of time tendencies
#ifdef USE_DRY_CANA_MASS
     tile%e_res_1 = canopy_air_mass*cpw*delta_qc*delta_Tc/delta_time
#else
     tile%e_res_1 = canopy_air_mass*(cpw-cp_air)*delta_qc*delta_Tc/delta_time
#endif
     tile%e_res_2 = delta_Tv*(clw*delta_Wl+csw*delta_Ws)/delta_time
! calculate the final value upward long-wave radiation flux from the land, to be 
! returned to the flux exchange.
     tile%lwup = ILa_dn - vegn_flw - flwg 

     if(is_watch_point())then
        write(*,*)'#### ground balance'
        __DEBUG2__(fswg,flwg)
        __DEBUG2__(sensg,evapg*grnd_latent)
        __DEBUG1__(grnd_flux)
        __DEBUG1__(Mg_imp)
        write(*,*)'#### implicit time steps'
        __DEBUG3__(delta_Tg, grnd_T,  grnd_T+delta_Tg )
        __DEBUG1__(delta_psig                         )
        __DEBUG3__(delta_qc, cana_q,  cana_q+delta_qc )
        __DEBUG3__(delta_Tc, cana_T,  cana_T+delta_Tc )
        __DEBUG3__(delta_Tv, vegn_T,  vegn_T+delta_Tv )
        __DEBUG3__(delta_wl, vegn_Wl, vegn_Wl+delta_wl)
        __DEBUG3__(delta_ws, vegn_Ws, vegn_Ws+delta_ws)
        __DEBUG2__(tile%e_res_1, tile%e_res_2)
        write(*,*)'#### resulting fluxes'
        __DEBUG4__(flwg, evapg, sensg, grnd_flux)
        __DEBUG3__(vegn_levap,vegn_fevap,vegn_uptk)
        __DEBUG2__(vegn_sens,vegn_flw)
        __DEBUG1__(Ea0+DEaDqc*delta_qc)
        __DEBUG2__(tile%cana%prog%q,cana_q)
     endif

     if (.not.prohibit_negative_canopy_water) exit ! do no corrections
     redo_leaf_water = .FALSE.
     if (vegn_Wl+delta_wl<0) then
        redo_leaf_water = .TRUE.
        Eli0 = vegn_Wl/delta_time + precip_l*vegn_ifrac - vegn_drip_l 
        DEliDTv = 0.0;  DEliDqc = 0.0
        DEliDwl = 0.0;  DEliDwf = 0.0
     endif
     if (vegn_Ws+delta_ws<0) then
        redo_leaf_water = .TRUE.
        Esi0 = vegn_Ws/delta_time + precip_s*vegn_ifrac - vegn_drip_s
        DEsiDTv = 0.0;  DEsiDqc = 0.0
        DEsiDwl = 0.0;  DEsiDwf = 0.0
     endif
     if (.not.redo_leaf_water) exit ! from loop
  enddo ! canopy_water_step
  
  call cana_step_2 ( tile%cana, delta_Tc, delta_qc )

  if(associated(tile%vegn)) then
     call vegn_step_2 ( tile%vegn, tile%diag, &
          delta_Tv, delta_wl, delta_ws, &
          vegn_melt,  &
          vegn_ovfl_l,   vegn_ovfl_s, &
          vegn_ovfl_Hl, vegn_ovfl_Hs )
     ! calculate total amount of liquid and solid precipitation below the canopy
     vegn_lprec  = (1-vegn_ifrac)*precip_l + vegn_drip_l + vegn_ovfl_l
     vegn_fprec  = (1-vegn_ifrac)*precip_s + vegn_drip_s + vegn_ovfl_s
     ! calculate heat carried by liquid and solid precipitation below the canopy
     vegn_hlprec = clw*((1-vegn_ifrac)*precip_l*(precip_T-tfreeze) &
                      + vegn_drip_l*(vegn_T+delta_Tv-tfreeze)) &
                      + vegn_ovfl_Hl
     vegn_hfprec = csw*((1-vegn_ifrac)*precip_s*(precip_T-tfreeze) &
                      + vegn_drip_s*(vegn_T+delta_Tv-tfreeze)) &
                      + vegn_ovfl_Hs
     ! make sure the temperature of the snow falling below canopy is below freezing
     ! this correction was introduced in an attempt to fix the problem with fictitious 
     ! heat accumulating in near-zero-mass snow; however it does not seem to make a 
     ! difference.
     if(vegn_hfprec>0)then
        ! solid precipitation from vegetation carries positive energy -- we can't have
        ! that, because that would bring snow T above tfreeze, so convert excess to 
        ! liquid
        delta_fprec = min(vegn_fprec,vegn_hfprec/hlf)
        vegn_fprec = vegn_fprec - delta_fprec
        vegn_lprec = vegn_lprec + delta_fprec
        vegn_hfprec = vegn_hfprec - hlf*delta_fprec
        ! we don't need to correct the vegn_hlprec since the temperature of additional
        ! liquid precip is tfreeze, and therefore its contribution to vegn_hlprec is
        ! exactly zero
     endif
     ! possibly we need to correct for the opposite situation: negative energy carried
     ! by liquid precipitation.
  else
     vegn_lprec  = precip_l
     vegn_fprec  = precip_s
     vegn_hlprec = precip_l*clw*(precip_T-tfreeze)
     vegn_hfprec = precip_s*csw*(precip_T-tfreeze)
     ! the fields below are only used in diagnostics
     vegn_melt   = 0
     vegn_fsw    = 0
  endif

  call snow_step_2 ( tile%snow, &
       snow_subl, vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
       delta_Tg, Mg_imp, evapg, fswg, flwg, sensg, &
       use_tfreeze_in_grnd_latent, &
       ! output:
       subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
       snow_fsw, snow_flw, snow_sens, &
       snow_levap, snow_fevap, snow_melt, &
       snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
       snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T )
  if(is_watch_point()) then
     write(*,*) 'subs_M_imp', subs_M_imp
  endif

  if (snow_active) then
     subs_G = snow_G_Z+snow_G_TZ*subs_DT
  else
     subs_G = 0
  endif
     
  if (associated(tile%glac)) then
     call glac_step_2 &
          ( tile%glac, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
          subs_DT, subs_M_imp, subs_evap, &
          subs_levap, subs_fevap, &
          subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop )
  else if (associated(tile%lake)) then
     call lake_step_2 &
          ( tile%lake, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
          subs_DT, subs_M_imp, subs_evap, &
          use_tfreeze_in_grnd_latent, subs_levap, subs_fevap, &
          subs_melt, subs_Ttop, subs_Ctop )
     subs_lrunf = 0.
     subs_hlrunf = 0.
  else if (associated(tile%soil)) then
     call soil_step_2 &
          ( tile%soil, tile%vegn, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
          vegn_uptk, subs_DT, subs_M_imp, subs_evap, &
          use_tfreeze_in_grnd_latent, &
          ! output:
          subs_levap, subs_fevap, &
          subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop )
  endif
     
! TEMP FIX: MAIN PROG SHOULD NOT TOUCH CONTENTS OF PROG VARS. ******
! ALSO, DIAGNOSTICS IN COMPONENT MODULES SHOULD _FOLLOW_ THIS ADJUSTMENT******
  if (LM2) then
     tile%snow%T = subs_Ttop
     subs_G2 = 0.
  else
     if (sum(tile%snow%ws(:))>0)then
        new_T = (subs_Ctop*subs_Ttop +snow_Cbot*snow_Tbot) &
                        / (subs_Ctop+snow_Cbot)
        tile%snow%T(size(tile%snow%T)) = new_T
        if(associated(tile%glac)) tile%glac%prog(1)%T = new_T
        if(associated(tile%lake)) tile%lake%prog(1)%T = new_T
        if(associated(tile%soil)) tile%soil%prog(1)%T = new_T
        subs_G2 = subs_Ctop*(new_T-subs_Ttop)/delta_time
     else
        if(tau_snow_T_adj>=0) then
           delta_T_snow = subs_Ctop*(subs_Ttop-snow_avrg_T)/&
                (subs_Ctop*tau_snow_T_adj/delta_time+subs_Ctop+snow_C)
           tile%snow%T(:) = snow_avrg_T + delta_T_snow

           new_T = subs_Ttop-snow_C/subs_Ctop*delta_T_snow
           if(associated(tile%glac)) tile%glac%prog(1)%T = new_T
           if(associated(tile%lake)) tile%lake%prog(1)%T = new_T
           if(associated(tile%soil)) tile%soil%prog(1)%T = new_T
           subs_G2 = subs_Ctop*(new_T-subs_Ttop)/delta_time
        else
           subs_G2 = 0.
        endif
     endif
  endif

  vegn_fco2 = 0
  if (associated(tile%vegn)) then
     ! do the calculations that require updated land surface prognostic variables
     call vegn_step_3 (tile%vegn, tile%soil, tile%cana%prog%T, precip_l+precip_s, &
          vegn_fco2, tile%diag)
     ! if vegn is present, then soil must be too
     call soil_step_3(tile%soil, tile%diag)
  endif
  ! update co2 concentration in the canopy air. It would be more consistent to do that
  ! in the same place and fashion as the rest of prognostic variables: that is, have the
  ! vegn_step_1 (and perhaps other *_step_1 procedures) calculate fluxes and their
  ! derivatives, then solve the linear equation(s), and finally have cana_step_2 update
  ! the concentration.
  if(update_cana_co2) then
     tile%cana%prog%co2 = tile%cana%prog%co2 + &
          (vegn_fco2 - fco2_0)/(canopy_air_mass_for_tracers/delta_time+Dfco2Dq)
  endif
  if(is_watch_point())then
     __DEBUG1__(tile%cana%prog%co2)
     __DEBUG3__(fco2_0,Dfco2Dq,vegn_fco2)
  endif
  
  call update_land_bc_fast (tile, i,j,k, land2cplr)

  runoff      = runoff      + (snow_frunf  + subs_lrunf  + snow_lrunf )*tile%frac
  runoff_heat = runoff_heat + (snow_hfrunf + subs_hlrunf + snow_hlrunf)*tile%frac
  runoff_snow = runoff_snow + snow_frunf*tile%frac
  hprec = (clw*precip_l+csw*precip_s)*(precip_T-tfreeze)
  hevap = cpw*land_evap*(evap_T-tfreeze)

  ! ---- diagnostic section ----------------------------------------------
  call send_tile_data(id_frac,    tile%frac,                          tile%diag)
  call send_tile_data(id_ntiles,  1.0,                                tile%diag)     
  call send_tile_data(id_precip,  precip_l+precip_s,                  tile%diag)
  call send_tile_data(id_hprec,   hprec,                              tile%diag)
  call send_tile_data(id_lprec,   precip_l,                           tile%diag)
  call send_tile_data(id_lprecv,  precip_l-vegn_lprec,                tile%diag)
  call send_tile_data(id_lprecs,  vegn_lprec-snow_lprec,              tile%diag)
  call send_tile_data(id_lprecg,  snow_lprec,                         tile%diag)
  call send_tile_data(id_hlprec,  clw*precip_l*(precip_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlprecv, clw*precip_l*(precip_T-tfreeze)-vegn_hlprec, &
                                                                      tile%diag)
  call send_tile_data(id_hlprecs, vegn_hlprec-snow_hlprec,            tile%diag)
  call send_tile_data(id_hlprecg, snow_hlprec,                        tile%diag)
  call send_tile_data(id_fprec,   precip_s,                           tile%diag)
  call send_tile_data(id_fprecv,  precip_s-vegn_fprec,                tile%diag)
  call send_tile_data(id_fprecs,  vegn_fprec,                         tile%diag)
  call send_tile_data(id_hfprec,  csw*precip_s*(precip_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfprecv, csw*precip_s*(precip_T-tfreeze)-vegn_hfprec, &
                                                                      tile%diag)
  call send_tile_data(id_hfprecs, vegn_hfprec,                        tile%diag)
  call send_tile_data(id_evap,    land_evap,                          tile%diag)
  call send_tile_data(id_hevap,   hevap,                              tile%diag)
  call send_tile_data(id_levap,   vegn_levap+snow_levap+subs_levap+vegn_uptk, &
                                                                      tile%diag)
  call send_tile_data(id_levapv,  vegn_levap,                         tile%diag)
  call send_tile_data(id_levaps,  snow_levap,                         tile%diag)
  call send_tile_data(id_levapg,  subs_levap,                         tile%diag)
  call send_tile_data(id_hlevap,  cpw*vegn_levap*(vegn_T-tfreeze) &
                                    +cpw*snow_levap*(snow_T-tfreeze) &
                                    +cpw*subs_levap*(grnd_T-tfreeze), tile%diag)
  call send_tile_data(id_hlevapv, cpw*vegn_levap*(vegn_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlevaps, cpw*snow_levap*(snow_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlevapg, cpw*subs_levap*(grnd_T-tfreeze),    tile%diag)
  call send_tile_data(id_fevap,   vegn_fevap+snow_fevap+subs_fevap,   tile%diag)
  call send_tile_data(id_fevapv,  vegn_fevap,                         tile%diag)
  call send_tile_data(id_fevaps,  snow_fevap,                         tile%diag)
  call send_tile_data(id_fevapg,  subs_fevap,                         tile%diag)
  call send_tile_data(id_hfevap,  cpw*vegn_fevap*(vegn_T-tfreeze) &
                                    +cpw*snow_fevap*(snow_T-tfreeze) &
                                    +cpw*subs_fevap*(grnd_T-tfreeze), tile%diag)
  call send_tile_data(id_hfevapv, cpw*vegn_fevap*(vegn_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfevaps, cpw*snow_fevap*(snow_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfevapg, cpw*subs_fevap*(grnd_T-tfreeze),    tile%diag)
  call send_tile_data(id_runf,    snow_lrunf+snow_frunf+subs_lrunf,   tile%diag)
  call send_tile_data(id_hrunf,   snow_hlrunf+snow_hfrunf+subs_hlrunf,tile%diag)
  call send_tile_data(id_lrunf,   snow_lrunf+subs_lrunf,              tile%diag)
  call send_tile_data(id_lrunfs,  snow_lrunf,                         tile%diag)
  call send_tile_data(id_lrunfg,  subs_lrunf,                         tile%diag)
  call send_tile_data(id_hlrunf,  snow_hlrunf+subs_hlrunf,            tile%diag)
  call send_tile_data(id_hlrunfs, snow_hlrunf,                        tile%diag)
  call send_tile_data(id_hlrunfg, subs_hlrunf,                        tile%diag)
  call send_tile_data(id_frunf,   snow_frunf,                         tile%diag)
  call send_tile_data(id_frunfs,  snow_frunf,                         tile%diag)
  call send_tile_data(id_hfrunf,  snow_hfrunf,                        tile%diag)
  call send_tile_data(id_hfrunfs, snow_hfrunf,                        tile%diag)
  call send_tile_data(id_melt,    vegn_melt+snow_melt+subs_melt,      tile%diag)
  call send_tile_data(id_meltv,   vegn_melt,                          tile%diag)
  call send_tile_data(id_melts,   snow_melt,                          tile%diag)
  call send_tile_data(id_meltg,   subs_melt,                          tile%diag)
  call send_tile_data(id_fsw,     vegn_fsw+snow_fsw+subs_fsw,         tile%diag)
  call send_tile_data(id_fswv,    vegn_fsw,                           tile%diag)
  call send_tile_data(id_fsws,    snow_fsw,                           tile%diag)
  call send_tile_data(id_fswg,    subs_fsw,                           tile%diag)
  call send_tile_data(id_flw,     vegn_flw+snow_flw+subs_flw,         tile%diag)
  call send_tile_data(id_flwv,    vegn_flw,                           tile%diag)
  call send_tile_data(id_flws,    snow_flw,                           tile%diag)
  call send_tile_data(id_flwg,    subs_flw,                           tile%diag)
  call send_tile_data(id_sens,    land_sens,                          tile%diag)
  call send_tile_data(id_sensv,   vegn_sens,                          tile%diag)
  call send_tile_data(id_senss,   snow_sens,                          tile%diag)
  call send_tile_data(id_sensg,   subs_sens,                          tile%diag)
  call send_tile_data(id_e_res_1, tile%e_res_1,                       tile%diag)
  call send_tile_data(id_e_res_2, tile%e_res_2,                       tile%diag)
  call send_tile_data(id_con_g_h, con_g_h,                            tile%diag)
  call send_tile_data(id_transp,  vegn_uptk,                          tile%diag)
  call send_tile_data(id_wroff,   snow_lrunf+subs_lrunf,              tile%diag)
  call send_tile_data(id_sroff,   snow_frunf,                         tile%diag)
  call send_tile_data(id_htransp, cpw*vegn_uptk*(vegn_T-tfreeze),     tile%diag)
  call send_tile_data(id_huptake, clw*vegn_uptk*(soil_uptake_T-tfreeze), &
                                                                      tile%diag)
  call send_tile_data(id_hroff,   snow_hlrunf+subs_hlrunf+snow_hfrunf, &
                                                                      tile%diag)
  call send_tile_data(id_gsnow,   subs_G,                             tile%diag)
  call send_tile_data(id_gequil,  subs_G2,                            tile%diag)
  call send_tile_data(id_grnd_flux, grnd_flux,                        tile%diag)
  call send_tile_data(id_soil_water_supply, soil_water_supply,        tile%diag)
  if(grnd_E_max.lt.0.5*HUGE(grnd_E_Max)) &
      call send_tile_data(id_levapg_max, grnd_E_max,                  tile%diag)
  call send_tile_data(id_qco2,    tile%cana%prog%co2,                 tile%diag)
  call send_tile_data(id_qco2_dvmr,&
       tile%cana%prog%co2*mol_air/mol_co2/(1-tile%cana%prog%q),       tile%diag)
  call send_tile_data(id_fco2,    vegn_fco2*mol_C/mol_co2,            tile%diag)
  call send_tile_data(id_swdn_dir, ISa_dn_dir,                        tile%diag)
  call send_tile_data(id_swdn_dif, ISa_dn_dif,                        tile%diag)
  call send_tile_data(id_swup_dir, ISa_dn_dir*tile%land_refl_dir,     tile%diag)
  call send_tile_data(id_swup_dif, ISa_dn_dif*tile%land_refl_dif,     tile%diag)
  call send_tile_data(id_lwdn,     ILa_dn,                            tile%diag)
  call send_tile_data(id_subs_emis,surf_emis_lw,                      tile%diag)

  if (id_total_C > 0) &
      call send_tile_data(id_total_C, land_tile_carbon(tile),         tile%diag)
end subroutine update_land_model_fast_0d


! ============================================================================
subroutine update_land_model_slow ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  integer :: i,j,k
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce, te

  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landSlowClock)

  call land_transitions( lnd%time )
  call update_vegn_slow( )
  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(lnd%tile_map, lnd%time)

  ! land_transitions may have changed the number of tiles per grid cell: reallocate 
  ! boundary conditions, if necessary
  call realloc_cplr2land( cplr2land )
  call realloc_land2cplr( land2cplr )
  ! set the land mask to FALSE everywhere -- update_land_bc_fast will set it to 
  ! true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0

  ! get the current state of the land boundary for the coupler
  ce = first_elmt(lnd%tile_map,                  &
              is=lbound(cplr2land%t_flux,1), &
              js=lbound(cplr2land%t_flux,2)  )
  te = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     ! calculate indices of the current tile in the input arrays;
     ! assume all the cplr2land components have the same lbounds
     call get_elmt_indices(ce,i,j,k)
     ! set this point coordinates as current for debug output
     call set_current_point(i,j,k)
     ! get pointer to current tile
     tile => current_tile(ce)
     ! advance enumerator to the next tile
     ce=next_elmt(ce)

     call update_land_bc_fast (tile, i,j,k, land2cplr, is_init=.true.)
  enddo

  call update_land_bc_slow( land2cplr )

  call mpp_clock_end(landClock)
  call mpp_clock_end(landSlowClock)
end subroutine update_land_model_slow


! ============================================================================
! solve for surface temperature. ensure that melt does not exceed available
! snow or soil ice (which would create energy-balance problems). also ensure
! that melt and temperature are consistent and that evaporation from liquid
! soil water does not exceed exfiltration rate limit, if one is supplied.
! because the possible combinations of active constraints has multiplied
! greatly, we do not allow phase change for any surface (i.e., soil) at which
! we might apply a constraint on Eg.
subroutine land_surface_energy_balance ( &
     ! surface parameters
     grnd_T,          & ! ground temperature
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface, kg/m2
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_E_min,      & ! Eg floor of 0 if condensation is prohibited
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     ! components of the ground energy balance linearization. Note that those
     ! are full derivatives, which include the response of the other surface
     ! scheme parameters to the change in ground temperature.
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, DflwgDpsig, & ! net long-wave
     Hg0,   DHgDTg,   DHgDpsig,   & ! sensible heat
     Eg0,   DEgDTg,   DEgDpsig,   & ! latent heat
     G0,    DGDTg,    & ! sub-surface heat 
     delta_Tg,        & ! surface temperature change for the time step
     delta_psig,      & 
     Mg_imp          )  ! implicit melt, kg/m2

  real, intent(in) :: & 
     grnd_T,          & ! ground temperature
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_E_min,      & ! Eg floor of 0 if condensation is prohibited
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, DflwgDpsig, & ! net long-wave
     Hg0,   DHgDTg,   DHgDpsig,   & ! sensible heat
     Eg0,   DEgDTg,   DEgDpsig,   & ! latent heat
     G0,    DGDTg                   ! sub-surface heat 
  real, intent(out) :: &
     delta_Tg,        & ! change in surface temperature
     delta_psig,      & ! change in surface soil-water matric head
     Mg_imp             ! mass of surface ice melted (or water frozen) during the 
                        ! time step, kg/m2

  real :: grnd_B     ! surface energy balance
  real :: grnd_DBDTg ! full derivative of grnd_B w.r.t. surface temperature
  real :: grnd_DBDpsig ! full derivative of grnd_B w.r.t. surface soil-water matric head
  real :: grnd_E_force, Eg_trial, Eg_check, determinant

  grnd_B      = fswg + flwg0      - Hg0      - grnd_latent*Eg0      - G0
  grnd_DBDTg  =        DflwgDTg   - DHgDTg   - grnd_latent*DEgDTg   - DGDTg
  grnd_DBDpsig =       DflwgDpsig - DHgDpsig - grnd_latent*DEgDpsig

  if(is_watch_point())then
     write(*,*)'#### ground balance input'
     __DEBUG1__(grnd_T)
     __DEBUG2__(grnd_liq, grnd_ice)
     __DEBUG1__(grnd_latent)
     __DEBUG1__(grnd_Tf)
     __DEBUG1__(grnd_E_min)
     __DEBUG1__(grnd_E_max)
     __DEBUG1__(fswg)
     __DEBUG3__(flwg0, DflwgDTg,DflwgDpsig)
     __DEBUG3__(Hg0,   DHgDTg,  DHgDpsig  )
     __DEBUG3__(Eg0,   DEgDTg,  DEgDpsig  )
     __DEBUG2__(G0,    DGDTg)
     write(*,*)'#### end of ground balance input'
     __DEBUG3__(grnd_B, grnd_DBDTg, grnd_DBDpsig)
  endif

  ! determine the ground temperature change under the assumptions that
  ! (1) no phase change occurs at surface (always true for soil now), and
  ! (2) delta_psig is zero (always true for snow, lake, glacier now)
  delta_Tg = - grnd_B/grnd_DBDTg
  delta_psig = 0
  ! calculate phase change on the ground, if necessary
  if     (grnd_ice>0.and.grnd_T+delta_Tg>grnd_Tf) then ! melt > 0
     Mg_imp =  min(grnd_ice,  grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  elseif (grnd_liq>0.and.grnd_T+delta_Tg<grnd_Tf) then ! melt < 0
     Mg_imp = -min(grnd_liq, -grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  else
     Mg_imp = 0
  endif
  ! adjust temperature change for the phase change
  delta_Tg = -(grnd_B - Mg_imp*hlf/delta_time)/grnd_DBDTg
  Eg_trial = Eg0 + DEgDTg*delta_Tg

  if(is_watch_point())then
     write(*,*)'#### ground balance solution with psig constant:'
     __DEBUG2__(grnd_B, grnd_DBDTg)
     __DEBUG3__(Mg_imp, delta_Tg, delta_psig)
     __DEBUG3__(grnd_E_min, Eg_trial, grnd_E_max)
  endif

  ! Solution above (assuming no change in psig) is acceptable if 
  ! it does not imply unacceptable value of Eg. this is always
  ! true for lake, glacier, snow. If implied Eg is outside prescribed
  ! bounds (a possibility only for soil), then Eg is set at bound (grnd_E_force)
  ! and the required Tg and psig are found. To accomplish this,
  ! we solve the system
  ! grnd_B + grnd_DBDTg*delta_Tg + grnd_DBDpsig*delta_psig = 0
  ! Eg0    +     DEgDTg*delta_Tg +     DEgDpsig*delta_psig = grnd_E_force
  ! There is no need to revisit the solution for phase change, which
  ! is only done explicitly for soil.

  if (Eg_trial.lt.grnd_E_min .or. Eg_trial.gt.grnd_E_max) then
      grnd_E_force = max(grnd_E_min, Eg_trial)
      grnd_E_force = min(grnd_E_max, grnd_E_force)
      determinant = grnd_DBDTg*DEgDpsig-grnd_DBDpsig*DEgDTg
      delta_Tg   = - (DEgDpsig*grnd_B + grnd_DBDpsig*(grnd_E_force-Eg0))/determinant
      delta_psig =   (DEgDTg  *grnd_B + grnd_DBDTg  *(grnd_E_force-Eg0))/determinant
      Eg_check = Eg0 + DEgDTg*delta_Tg + DEgDpsig*delta_psig
      Mg_imp = 0
         if(is_watch_point())then
            write(*,*)'#### trial solution violated Eg limit, new soln:'
            __DEBUG2__(grnd_B,grnd_DBDTg)
            __DEBUG3__(Mg_imp,delta_Tg,delta_psig)
            __DEBUG3__(grnd_E_min, Eg_check, grnd_E_max)
         endif
    endif
end subroutine land_surface_energy_balance


! ============================================================================
subroutine update_land_bc_fast (tile, i,j,k, land2cplr, is_init)
  type(land_tile_type), intent(inout) :: tile
  integer             , intent(in) :: i,j,k
  type(land_data_type), intent(inout) :: land2cplr
  logical, optional :: is_init

  ! ---- local vars
  real :: &
         grnd_T, subs_z0m, subs_z0s, &
                 snow_z0s, snow_z0m, &
         snow_area, snow_depth

  real :: subs_refl_dir(NBANDS), subs_refl_dif(NBANDS) ! direct and diffuse albedos
  real :: subs_refl_lw ! reflectance for thermal radiation
  real :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS) ! direct and diffuse albedos of snow
  real :: snow_refl_lw ! snow reflectance for thermal radiation
  real :: snow_emis ! snow emissivity
  real :: grnd_emis ! ground emissivity
  ! NOTE :  grnd_emis is used only to satisfy xxxx_radiation interfaces; its value is ignored, but
  ! 1-refl is used instead. snow_emis is used in the the vegn_radiation, but it shouldn't be since
  ! properties of intercepted snowpack are, in general, different from the snow on the ground 
  real :: snow_area_rad ! "snow area for radiation calculations" -- introduced
                        ! to reproduce lm2 behavior
  real :: vegn_refl_lw, vegn_tran_lw ! reflectance and transmittance of vegetation for thermal radiation
  real :: vegn_refl_dif(NBANDS), vegn_tran_dif(NBANDS) ! reflectance and transmittance of vegetation for diffuse light
  real :: vegn_refl_dir(NBANDS), vegn_tran_dir(NBANDS) ! reflectance and transmittance of vegetation for direct light
  real :: vegn_tran_dir_dir(NBANDS) ! (?)
  real :: &
         vegn_Tv,     &
         vegn_cover,  &
         vegn_height, vegn_lai, vegn_sai, vegn_d_leaf, cana_co2
  logical :: do_update

  real :: cosz    ! cosine of solar zenith angle
  real :: fracday ! daytime fraction of time interval
  real :: rrsun   ! earth-sun distance (r) relative to semi-major axis
                  ! of orbital ellipse (a) : (a/r)**2
  integer :: face ! for debugging
  vegn_Tv = 0

  do_update = .not.present(is_init)

  ! on initialization the albedos are calculated for the current time step ( that is, interval
  ! lnd%time, lnd%time+lnd%dt_fast); in the course of the run this subroutine is called
  ! at the end of time step (but before time is advanced) to calculate the radiative properties 
  ! for the _next_ time step
  if (do_update) then
     call diurnal_solar(lnd%lat(i,j), lnd%lon(i,j), lnd%time+lnd%dt_fast, &
          cosz, fracday, rrsun, lnd%dt_fast)
  else
     call diurnal_solar(lnd%lat(i,j), lnd%lon(i,j), lnd%time, &
          cosz, fracday, rrsun, lnd%dt_fast)
  endif
  
  if (associated(tile%glac)) then
     call glac_radiation(tile%glac, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call glac_diffusion(tile%glac, subs_z0s, subs_z0m )
  else if (associated(tile%lake)) then
     call lake_radiation(tile%lake, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call lake_diffusion(tile%lake, subs_z0s, subs_z0m )
  else if (associated(tile%soil)) then
     call soil_radiation(tile%soil, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call soil_diffusion(tile%soil, subs_z0s, subs_z0m )
  else
     call get_current_point(face=face)
     call error_mesg('update_land_model_fast','none of the surface tiles exist at ('//&
             trim(string(i))//','//trim(string(j))//','//trim(string(k))//&
             ', face='//trim(string(face))//')',FATAL)
  endif

  call snow_radiation ( tile%snow%T(1), cosz, snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
  call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
  call snow_diffusion ( tile%snow, snow_z0s, snow_z0m )

  if (associated(tile%vegn)) then
     call update_derived_vegn_data(tile%vegn)
     ! USE OF SNOWPACK RAD PROPERTIES FOR INTERCEPTED SNOW IS ERRONEOUS,
     ! NEEDS TO BE CHANGED. TEMPORARY.
     call vegn_radiation ( tile%vegn, cosz, snow_depth, snow_refl_dif, snow_emis, &
                   vegn_refl_dif, vegn_tran_dif, &
                   vegn_refl_dir, vegn_tran_dir, vegn_tran_dir_dir, &
                   vegn_refl_lw, vegn_tran_lw)
     ! (later see if we can remove vegn_cover from c-a-radiation...) TEMPORARY
     call vegn_get_cover ( tile%vegn, snow_depth, vegn_cover)
     call vegn_diffusion ( tile%vegn, vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf)
  else
     ! set radiative properties for null vegetation
     vegn_refl_dif     = 0
     vegn_tran_dif     = 1
     vegn_refl_dir     = 0
     vegn_tran_dir     = 0
     vegn_tran_dir_dir = 1
     vegn_refl_lw      = 0 
     vegn_tran_lw      = 1
     ! set cover for null vegetation
     vegn_cover        = 0
     ! set other parameters for null vegetation
     vegn_height       = 0
     vegn_lai          = 0
     vegn_sai          = 0
     vegn_d_leaf       = 0
  endif

  ! store the values of long-wave optical properties to be used in the update_land_model_fast
  tile%surf_refl_lw = subs_refl_lw  + (snow_refl_lw  - subs_refl_lw ) * snow_area
  tile%vegn_refl_lw = vegn_refl_lw
  tile%vegn_tran_lw = vegn_tran_lw

  
  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 1 ####'
     __DEBUG3__(cosz, fracday, rrsun)
     __DEBUG2__(vegn_lai,vegn_sai)
     __DEBUG1__(subs_refl_dif)
     __DEBUG1__(subs_refl_dir)
     __DEBUG1__(vegn_refl_dif)
     __DEBUG1__(vegn_tran_dif)
     __DEBUG1__(vegn_refl_dir)
     __DEBUG1__(vegn_tran_dir)
     __DEBUG1__(vegn_tran_dir_dir)
     __DEBUG2__(vegn_refl_lw,vegn_tran_lw)
     write(*,*) '#### update_land_bc_fast ### end of checkpoint 1 ####'
  endif

  snow_area_rad = snow_area
  if (lm2) then
     if(associated(tile%glac)                    ) snow_area_rad = 0
     if(associated(tile%soil).and.vegn_cover>0.01) snow_area_rad = 1
  endif

  call cana_radiation( lm2, &
       subs_refl_dir, subs_refl_dif, subs_refl_lw, &
       snow_refl_dir, snow_refl_dif, snow_refl_lw, &
       snow_area_rad,  &
       vegn_refl_dir, vegn_refl_dif, vegn_tran_dir, vegn_tran_dif, &
       vegn_tran_dir_dir, vegn_refl_lw, vegn_tran_lw, &
       vegn_cover, &
       tile%Sg_dir, tile%Sg_dif, tile%Sv_dir, tile%Sv_dif, &
       tile%land_refl_dir, tile%land_refl_dif )

  call cana_roughness( lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover,  vegn_height, vegn_lai, vegn_sai, &
     tile%land_d, tile%land_z0m, tile%land_z0s, &
     associated(tile%lake).or.associated(tile%glac))

  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 2 ####'
     write(*,*) 'Sg_dir', tile%Sg_dir
     write(*,*) 'Sg_dif', tile%Sg_dif
     write(*,*) 'Sv_dir', tile%Sv_dir
     write(*,*) 'Sv_dif', tile%Sv_dif
     write(*,*) 'land_albedo_dir', tile%land_refl_dir
     write(*,*) 'land_albedo_dif', tile%land_refl_dif
     write(*,*) 'land_z0m', tile%land_z0m
     write(*,*) '#### update_land_bc_fast ### end of checkpoint 2 ####'
  endif

  land2cplr%t_surf         (i,j,k) = tfreeze
  land2cplr%t_ca           (i,j,k) = tfreeze
  land2cplr%tr             (i,j,k, lnd%isphum) = 0.0
  land2cplr%albedo         (i,j,k) = 0.0
  land2cplr%albedo_vis_dir (i,j,k) = 0.0
  land2cplr%albedo_nir_dir (i,j,k) = 0.0
  land2cplr%albedo_vis_dif (i,j,k) = 0.0
  land2cplr%albedo_nir_dif (i,j,k) = 0.0
  land2cplr%rough_mom      (i,j,k) = 0.1
  land2cplr%rough_heat     (i,j,k) = 0.1

  ! Calculate radiative surface temperature. lwup can't be calculated here
  ! based on the available temperatures because it's a result of the implicit 
  ! time step: lwup = lwup0 + DlwupDTg*delta_Tg + ..., so we have to carry it
  ! from the update_land_fast 
  ! Consequence: since update_landbc_fast is called once at the very beginning of 
  ! every run (before update_land_fast is called) lwup from the previous step 
  ! must be stored in the in the restart for reproducibility
  land2cplr%t_surf(i,j,k) = ( tile%lwup/stefan ) ** 0.25

  if (associated(tile%glac)) call glac_get_sfc_temp(tile%glac, grnd_T)
  if (associated(tile%lake)) call lake_get_sfc_temp(tile%lake, grnd_T)
  if (associated(tile%soil)) call soil_get_sfc_temp(tile%soil, grnd_T)
  if (snow_area > 0)         call snow_get_sfc_temp(tile%snow, grnd_T)

  ! set the boundary conditions for the flux exchange
  land2cplr%mask           (i,j,k) = .TRUE.
  land2cplr%tile_size      (i,j,k) = tile%frac

  call cana_state ( tile%cana, land2cplr%t_ca(i,j,k), &
                               land2cplr%tr(i,j,k,lnd%isphum), cana_co2)
!  land2cplr%t_ca           (i,j,k) = tile%cana%prog%T
!  land2cplr%tr             (i,j,k,lnd%isphum) = tile%cana%prog%q
  if(lnd%ico2/=NO_TRACER) then
     land2cplr%tr(i,j,k,lnd%ico2) = cana_co2
  endif
  land2cplr%albedo_vis_dir (i,j,k) = tile%land_refl_dir(BAND_VIS)
  land2cplr%albedo_nir_dir (i,j,k) = tile%land_refl_dir(BAND_NIR)
  land2cplr%albedo_vis_dif (i,j,k) = tile%land_refl_dif(BAND_VIS)
  land2cplr%albedo_nir_dif (i,j,k) = tile%land_refl_dif(BAND_NIR)
  land2cplr%albedo         (i,j,k) = SUM(tile%land_refl_dir + tile%land_refl_dif)/4 ! incorrect, replace with proper weighting later
  land2cplr%rough_mom      (i,j,k) = tile%land_z0m
  land2cplr%rough_heat     (i,j,k) = tile%land_z0s

  if(is_watch_point()) then
     write(*,*)'#### update_land_bc_fast ### output ####'
     write(*,*)'land2cplr%mask',land2cplr%mask(i,j,k)
     write(*,*)'land2cplr%tile_size',land2cplr%tile_size(i,j,k)
     write(*,*)'land2cplr%t_surf',land2cplr%t_surf(i,j,k)
     write(*,*)'land2cplr%t_ca',land2cplr%t_ca(i,j,k)
     write(*,*)'land2cplr%albedo',land2cplr%albedo(i,j,k)
     write(*,*)'land2cplr%rough_mom',land2cplr%rough_mom(i,j,k)
     write(*,*)'land2cplr%rough_heat',land2cplr%rough_heat(i,j,k)
     write(*,*)'land2cplr%tr',land2cplr%tr(i,j,k,:)
     write(*,*)'#### update_land_bc_fast ### end of output ####'
  endif

  ! ---- diagnostic section
  call send_tile_data(id_vegn_cover, vegn_cover, tile%diag)
  call send_tile_data(id_cosz, cosz, tile%diag)
  call send_tile_data(id_albedo_dir, tile%land_refl_dir, tile%diag)
  call send_tile_data(id_albedo_dif, tile%land_refl_dif, tile%diag)
  call send_tile_data(id_vegn_refl_dir, vegn_refl_dir,     tile%diag)
  call send_tile_data(id_vegn_refl_dif, vegn_refl_dif, tile%diag)
  call send_tile_data(id_vegn_refl_lw,  vegn_refl_lw, tile%diag)
  call send_tile_data(id_vegn_tran_dir, vegn_tran_dir_dir, tile%diag)
  call send_tile_data(id_vegn_tran_dif, vegn_tran_dif, tile%diag)
  call send_tile_data(id_vegn_tran_lw,  vegn_tran_lw, tile%diag)
  call send_tile_data(id_vegn_sctr_dir, vegn_tran_dir,     tile%diag)
  call send_tile_data(id_subs_refl_dir, subs_refl_dir, tile%diag)
  call send_tile_data(id_subs_refl_dif, subs_refl_dif, tile%diag)
  call send_tile_data(id_grnd_T,     grnd_T,     tile%diag)

  ! --- debug section
  call check_temp_range(land2cplr%t_ca(i,j,k),'update_land_bc_fast','T_ca',lnd%time)

end subroutine update_land_bc_fast


! ============================================================================
subroutine update_land_bc_slow (land2cplr)
  type(land_data_type), intent(inout) :: land2cplr

  ! ---- local vars 
  integer :: i,j,k,face ! coordinates of the watch point, for debug printout

  call update_topo_rough(land2cplr%rough_scale)
  where (land2cplr%mask) &
       land2cplr%rough_scale = max(land2cplr%rough_mom,land2cplr%rough_scale)
  call get_watch_point(i,j,k,face)
  if ( lnd%face==face.and.             &
       lnd%is<=i.and.i<=lnd%ie.and.    &
       lnd%js<=j.and.j<=lnd%je.and.    &
       k<=size(land2cplr%rough_scale,3)) then
     write(*,*)'#### update_land_bc_slow ### output ####'
     write(*,*)'land2cplr%rough_scale',land2cplr%rough_scale(i,j,k)
     write(*,*)'#### update_land_bc_slow ### end of output ####'
  endif
    
end subroutine update_land_bc_slow


! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)

type(land_data_type), intent(in)  :: bnd 
integer             , intent(in)  :: index
real                , intent(out) :: value ! Domain water (Kg) or heat (Joules)

integer :: i,j,n
type(land_tile_enum_type)     :: ce,te
type(land_tile_type), pointer :: tile
character(len=128) :: message
integer :: is,ie,js,je
real :: area_factor, river_value
! *_twd are tile water densities (kg water per m2 of tile)
real twd_gas_cana,               twd_liq_glac, twd_sol_glac, &
     twd_liq_lake, twd_sol_lake, twd_liq_soil, twd_sol_soil, &
     twd_liq_snow, twd_sol_snow, twd_liq_vegn, twd_sol_vegn
! *_gcwd are grid-cell water densities (kg water per m2 of land in grid cell)
real gcwd_cana, gcwd_glac, gcwd_lake, gcwd_soil, gcwd_snow, gcwd_vegn
! v_* are global masses of water
real v_cana, v_glac, v_lake, v_soil, v_snow, v_vegn
real cana_q, a_globe

value = 0.0
v_cana = 0.
v_glac = 0.
v_lake = 0.
v_soil = 0.
v_snow = 0.
v_vegn = 0.
if(.not.bnd%pe) return
is = lnd%is
ie = lnd%ie
js = lnd%js
je = lnd%je

! The following is a dirty getaround
if(lnd%cellarea(is,js) < 1.0) then
  area_factor = 4*pi*radius**2 ! lnd%area is fraction of globe
else
  area_factor = 1.0 ! lnd%area is actual area (m**2)
endif

select case(index)
case(ISTOCK_WATER)
  do j = js, je
  do i = is, ie
    ce = first_elmt(lnd%tile_map(i,j))
    te = tail_elmt (lnd%tile_map(i,j))
    gcwd_cana = 0.0; gcwd_glac = 0.0; gcwd_lake = 0.0
    gcwd_soil = 0.0; gcwd_snow = 0.0; gcwd_vegn = 0.0
    do while(ce /= te)
      tile => current_tile(ce)
      twd_gas_cana = 0.0
      twd_liq_glac = 0.0 ; twd_sol_glac = 0.0
      twd_liq_lake = 0.0 ; twd_sol_lake = 0.0
      twd_liq_soil = 0.0 ; twd_sol_soil = 0.0
      twd_liq_snow = 0.0 ; twd_sol_snow = 0.0
      twd_liq_vegn = 0.0 ; twd_sol_vegn = 0.0
      if(associated(tile%cana)) then
        call cana_state ( tile%cana, cana_q=cana_q )
        twd_gas_cana = canopy_air_mass*cana_q
        endif
      if(associated(tile%glac)) &
        call glac_tile_stock_pe(tile%glac, twd_liq_glac, twd_sol_glac)
      if(associated(tile%lake)) &
        call lake_tile_stock_pe(tile%lake, twd_liq_lake, twd_sol_lake)
      if(associated(tile%soil)) &
        call soil_tile_stock_pe(tile%soil, twd_liq_soil, twd_sol_soil)
      if(associated(tile%snow)) &
        call snow_tile_stock_pe(tile%snow, twd_liq_snow, twd_sol_snow)
      if(associated(tile%vegn)) &
        call vegn_tile_stock_pe(tile%vegn, twd_liq_vegn, twd_sol_vegn)
      gcwd_cana = gcwd_cana +  twd_gas_cana                 * tile%frac
      gcwd_glac = gcwd_glac + (twd_liq_glac + twd_sol_glac) * tile%frac
      gcwd_lake = gcwd_lake + (twd_liq_lake + twd_sol_lake) * tile%frac
      gcwd_soil = gcwd_soil + (twd_liq_soil + twd_sol_soil) * tile%frac
      gcwd_snow = gcwd_snow + (twd_liq_snow + twd_sol_snow) * tile%frac
      gcwd_vegn = gcwd_vegn + (twd_liq_vegn + twd_sol_vegn) * tile%frac
      ce=next_elmt(ce)
    enddo
    v_cana = v_cana + gcwd_cana * lnd%area(i,j)*area_factor
    v_glac = v_glac + gcwd_glac * lnd%area(i,j)*area_factor
    v_lake = v_lake + gcwd_lake * lnd%area(i,j)*area_factor
    v_soil = v_soil + gcwd_soil * lnd%area(i,j)*area_factor
    v_snow = v_snow + gcwd_snow * lnd%area(i,j)*area_factor
    v_vegn = v_vegn + gcwd_vegn * lnd%area(i,j)*area_factor
  enddo
  enddo
  value  = v_cana + v_glac + v_lake + v_soil + v_snow + v_vegn
a_globe = 4. * pi * radius**2
case(ISTOCK_HEAT)
  if(.not.stock_warning_issued) then
    call error_mesg('Lnd_stock_pe','Heat stock not yet implemented',NOTE)
    stock_warning_issued = .true.
  endif
! do j = js, je
! do i = is, ie
!   ce = first_elmt(lnd%tile_map(i,j))
!   te = tail_elmt (lnd%tile_map(i,j))
!   grid_cell_heat_density = 0.0
!   do while(ce /= te)
!     tile => current_tile(ce)
!     tile_heat_density = 0.0
!     if(associated(tile%soil)) then
!       do n=1, size(tile%soil%prog)
!       tile_heat_density = tile_heat_density + (tile%soil%prog(n)%T-tfreeze)* &
!                  (tile%soil%heat_capacity_dry(n)*dz(n) + &
!                   clw*tile%soil%prog(n)%wl             + &
!                   csw*tile%soil%prog(n)%ws)
!       enddo
!       tile_heat_density = tile_heat_density + clw*soil%prog(1)%groundwater*(soil%prog(1)%groundwater_T-tfreeze) ! Why is this outside n loop?
!     endif
!     grid_cell_heat_density = grid_cell_heat_density + tile_heat_density * tile%frac
!     ce=next_elmt(ce)
!   enddo
!   grid_cell_heat = grid_cell_heat_density * lnd%area(i,j)*area_factor
!   value = value + grid_cell_heat
! enddo
! enddo
case(ISTOCK_SALT)
  if(.not.stock_warning_issued) then
    call error_mesg('Lnd_stock_pe','Salt stock not yet implemented',NOTE)
    stock_warning_issued = .true.
  endif
case default
  write(message,'(i2,a,i2,a,i2,a,i2,a)') &
  index,' is an invalid stock index. Must be ISTOCK_WATER or ISTOCK_HEAT or ISTOCK_SALT (',ISTOCK_WATER,' or ',ISTOCK_HEAT,' or ', ISTOCK_SALT,')'
  call error_mesg('Lnd_stock_pe',message,FATAL)
end select

call river_stock_pe(index, river_value)
value = value + river_value

if (index.eq.ISTOCK_WATER.and.give_stock_details) then
    call mpp_sum(river_value, pelist=lnd%pelist)
    call mpp_sum(v_cana, pelist=lnd%pelist)
    call mpp_sum(v_glac, pelist=lnd%pelist)
    call mpp_sum(v_lake, pelist=lnd%pelist)
    call mpp_sum(v_soil, pelist=lnd%pelist)
    call mpp_sum(v_snow, pelist=lnd%pelist)
    call mpp_sum(v_vegn, pelist=lnd%pelist)
    write (message,'(a,f10.5)') 'total land storage:',v_cana/a_globe+v_glac/a_globe+ &
        v_lake/a_globe+v_soil/a_globe+v_snow/a_globe+v_vegn/a_globe+river_value/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '...canopy air:',v_cana/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '......glacier:',v_glac/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........lake:',v_lake/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........soil:',v_soil/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........snow:',v_snow/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........vegn:',v_vegn/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.......rivers:',river_value/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
endif

end subroutine Lnd_stock_pe


! ============================================================================
! initialize horizontal axes for land grid so that all sub-modules can use them,
! instead of creating their own
subroutine land_diag_init(clonb, clatb, clon, clat, time, domain, &
     id_lon, id_lat, id_band)
  real, intent(in) :: &
       clonb(:), clatb(:), & ! longitudes and latitudes of grid cells vertices,
                             ! specified for the global grid
       clon(:),  clat(:)     ! lon and lat of respective grid cell centers
  type(time_type), intent(in) :: time ! initial time for diagnostic fields
  type(domain2d), intent(in)  :: domain
  integer, intent(out) :: &
       id_lon, id_lat, id_band   ! IDs of respective diag. manager axes

  ! ---- local vars ----------------------------------------------------------
  integer :: id_lonb, id_latb ! IDs for cell boundaries
  integer :: nlon, nlat       ! sizes of respective axes
  integer :: axes(2)          ! array of axes for 2-D fields
  integer :: i

  nlon = size(clon)
  nlat = size(clat)

  if(mpp_get_ntile_count(lnd%domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define longitude axes and its edges
     id_lonb = diag_axis_init ( &
          'lonb', clonb, 'degrees_E', 'X', 'longitude edges', &
          set_name='land', domain2=domain )
     id_lon  = diag_axis_init (                                                &
          'lon',  clon, 'degrees_E', 'X',  &
          'longitude', set_name='land',  edges=id_lonb, domain2=domain )

     ! define latitude axes and its edges
     id_latb = diag_axis_init ( &
          'latb', clatb, 'degrees_N', 'Y', 'latitude edges',  &
          set_name='land',  domain2=domain   )
     id_lat = diag_axis_init (                                                &
          'lat',  clat, 'degrees_N', 'Y', &
          'latitude', set_name='land', edges=id_latb, domain2=domain   )
  else
     id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
          'T-cell longitude', set_name='land',  domain2=domain, aux='geolon_t' )
     id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
          'T-cell latitude', set_name='land',  domain2=domain, aux='geolat_t' )
  endif
  id_band = diag_axis_init (                                                &
       'band',  (/1.0,2.0/), 'unitless', 'Z', &
       'spectral band', set_name='land' )

  ! set up an array of axes, for convenience
  axes = (/id_lon, id_lat/)

  ! register auxilary coordinate variables

  id_geolon_t = register_static_field ( module_name, 'geolon_t', axes, &
       'longitude of grid cell centers', 'degrees_E', missing_value = -1.0e+20 )
  id_geolat_t = register_static_field ( module_name, 'geolat_t', axes, &
       'latitude of grid cell centers', 'degrees_N', missing_value = -1.0e+20 )

  ! register static diagnostic fields

  id_cellarea = register_static_field ( module_name, 'cell_area', axes, &
       'total area in grid cell', 'm2', missing_value=-1.0 )
  id_landarea = register_static_field ( module_name, 'land_area', axes, &
       'land area in grid cell', 'm2', missing_value=-1.0 )
  id_landfrac = register_static_field ( module_name, 'land_frac', axes, &
       'fraction of land in grid cell','unitless', missing_value=-1.0 ) 
  id_no_riv = register_static_field ( module_name, 'no_riv', axes, &
       'indicator of land without rivers','unitless', missing_value=-1.0 ) 

  ! register regular (dynamic) diagnostic fields

  id_VWS = register_tiled_diag_field ( module_name, 'VWS', axes, time, &
             'vapor storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_VWSc    = register_tiled_diag_field ( module_name, 'VWSc', axes, time, &
             'vapor mass in canopy air', 'kg/m2', missing_value=-1.0e+20 )
  id_LWS = register_tiled_diag_field ( module_name, 'LWS', axes, time, &
             'liquid storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSv    = register_tiled_diag_field ( module_name, 'LWSv', axes, time, &
             'liquid interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSs    = register_tiled_diag_field ( module_name, 'LWSs', axes, time, &
             'liquid storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSg    = register_tiled_diag_field ( module_name, 'LWSg', axes, time, &
             'liquid ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWS = register_tiled_diag_field ( module_name, 'FWS', axes, time, &
             'frozen storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSv    = register_tiled_diag_field ( module_name, 'FWSv', axes, time, &
             'frozen interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSs    = register_tiled_diag_field ( module_name, 'FWSs', axes, time, &
             'frozen storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSg    = register_tiled_diag_field ( module_name, 'FWSg', axes, time, &
             'frozen ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_HS = register_tiled_diag_field ( module_name, 'HS', axes, time, &
             'land heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSv     = register_tiled_diag_field ( module_name, 'HSv', axes, time, &
             'interception heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSs     = register_tiled_diag_field ( module_name, 'HSs', axes, time, &
             'heat storage in snowpack', 'J/m2', missing_value=-1.0e+20 )
  id_HSg     = register_tiled_diag_field ( module_name, 'HSg', axes, time, &
             'ground heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSc     = register_tiled_diag_field ( module_name, 'HSc', axes, time, &
             'canopy-air heat storage', 'J/m2', missing_value=-1.0e+20 )

  id_dis_liq   = register_diag_field ( module_name, 'dis_liq', axes, &
       time, 'liquid discharge to ocean', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_dis_ice   = register_diag_field ( module_name, 'dis_ice', axes, &
       time, 'ice discharge to ocean', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_dis_heat   = register_diag_field ( module_name, 'dis_heat', axes, &
       time, 'heat of mass discharge to ocean', 'W/m2', missing_value=-1.0e+20 )
  id_dis_sink   = register_diag_field ( module_name, 'dis_sink', axes, &
       time, 'burial rate of small/negative discharge', 'kg/(m2 s)', missing_value=-1.0e+20 )

  id_precip = register_tiled_diag_field ( module_name, 'precip', axes, time, &
             'precipitation rate', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hprec = register_tiled_diag_field ( module_name, 'hprec', axes, time, &
             'sensible heat of precipitation', 'W/m2', missing_value=-1.0e+20 )
  id_lprec = register_tiled_diag_field ( module_name, 'lprec_l', axes, time, &
             'rainfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lprecv = register_tiled_diag_field ( module_name, 'lprecv', axes, time, &
             'net rainfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecs = register_tiled_diag_field ( module_name, 'lprecs', axes, time, &
             'rainfall to snow, minus drainage', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecg = register_tiled_diag_field ( module_name, 'lprecg', axes, time, &
             'effective rainfall to ground sfc', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hlprec  = register_tiled_diag_field ( module_name, 'hlprec', axes, time, &
             'total liq precipitation heat', 'W/m2', missing_value=-1.0e+20)
  id_hlprecv = register_tiled_diag_field ( module_name, 'hlprecv', axes, time, &
             'net liq heat to vegetation', 'W/m2', missing_value=-1.0e+20)
  id_hlprecs = register_tiled_diag_field ( module_name, 'hlprecs', axes, time, &
             'net liq heat to snow', 'W/m2', missing_value=-1.0e+20)
  id_hlprecg = register_tiled_diag_field ( module_name, 'hlprecg', axes, time, &
             'net liq heat to ground sfc', 'W/m2', missing_value=-1.0e+20)
  id_fprec = register_tiled_diag_field ( module_name, 'fprec_l', axes, time, &
             'snowfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fprecv = register_tiled_diag_field ( module_name, 'fprecv', axes, time, &
             'net snowfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fprecs = register_tiled_diag_field ( module_name, 'fprecs', axes, time, &
             'effective snowfall to snowpack', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hfprec = register_tiled_diag_field ( module_name, 'hfprec', axes, time, &
             'sens heat of snowfall', 'W/m2', missing_value=-1.0e+20)
  id_hfprecv = register_tiled_diag_field ( module_name, 'hfprecv', axes, time, &
             'net sol heat to vegetation', 'W/m2', missing_value=-1.0e+20)
  id_hfprecs = register_tiled_diag_field ( module_name, 'hfprecs', axes, time, &
             'net sol heat to snow', 'W/m2', missing_value=-1.0e+20)
  id_evap = register_tiled_diag_field ( module_name, 'evap', axes, time, &
             'vapor flux up from land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hevap = register_tiled_diag_field ( module_name, 'hevap', axes, time, &
             'sensible heat of evap', 'W/m2', missing_value=-1.0e+20 )
  id_levap   = register_tiled_diag_field ( module_name, 'levap', axes, time, &
             'vapor flux from all liq (inc Tr)', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_levapv = register_tiled_diag_field ( module_name, 'levapv', axes, time, &
             'vapor flux leaving intercepted liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levaps = register_tiled_diag_field ( module_name, 'levaps', axes, time, &
             'vapor flux leaving snow liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levapg = register_tiled_diag_field ( module_name, 'levapg', axes, time, &
             'vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hlevap = register_tiled_diag_field ( module_name, 'hlevap', axes, time, &
             'vapor flux heat from liq source', 'W/m2', missing_value=-1.0e+20)
  id_hlevapv = register_tiled_diag_field ( module_name, 'hlevapv', axes, time, &
             'vapor heat from liq interc', 'W/m2', missing_value=-1.0e+20)
  id_hlevaps = register_tiled_diag_field ( module_name, 'hlevaps', axes, time, &
             'vapor heat from snow liq', 'W/m2', missing_value=-1.0e+20)
  id_hlevapg = register_tiled_diag_field ( module_name, 'hlevapg', axes, time, &
             'vapor heat from ground liq', 'W/m2', missing_value=-1.0e+20)
  id_fevap   = register_tiled_diag_field ( module_name, 'fevap', axes, time, &
             'vapor flux from all ice', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fevapv = register_tiled_diag_field ( module_name, 'fevapv', axes, time, &
             'vapor flux leaving vegn ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevaps = register_tiled_diag_field ( module_name, 'fevaps', axes, time, &
             'vapor flux leaving snow ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevapg = register_tiled_diag_field ( module_name, 'fevapg', axes, time, &
             'vapor flux leaving ground ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hfevap = register_tiled_diag_field ( module_name, 'hfevap', axes, time, &
             'vapor flux heat from solid source', 'W/m2', missing_value=-1.0e+20)
  id_hfevapv = register_tiled_diag_field ( module_name, 'hfevapv', axes, time, &
             'vapor heat from sol interc', 'W/m2', missing_value=-1.0e+20)
  id_hfevaps = register_tiled_diag_field ( module_name, 'hfevaps', axes, time, &
             'vapor heat from snow sol', 'W/m2', missing_value=-1.0e+20)
  id_hfevapg = register_tiled_diag_field ( module_name, 'hfevapg', axes, time, &
             'vapor heat from ground sol', 'W/m2', missing_value=-1.0e+20)
  id_runf   = register_tiled_diag_field ( module_name, 'runf', axes, time, &
             'total runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hrunf   = register_tiled_diag_field ( module_name, 'hrunf', axes, time, &
             'sensible heat of total runoff', 'W/m2', missing_value=-1.0e+20 )
  id_lrunf   = register_tiled_diag_field ( module_name, 'lrunf', axes, time, &
             'total rate of liq runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfs  = register_tiled_diag_field ( module_name, 'lrunfs', axes, time, &
             'rate of liq runoff via calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfg  = register_tiled_diag_field ( module_name, 'lrunfg', axes, time, &
             'rate of liq runoff, ground', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hlrunf  = register_tiled_diag_field ( module_name, 'hlrunf', axes, time, &
             'heat of liq runoff', 'W/m2', missing_value=-1.0e+20 )
  id_hlrunfs  = register_tiled_diag_field ( module_name, 'hlrunfs', axes, time, &
             'heat of liq runoff from snow pack', 'W/m2', missing_value=-1.0e+20 )
  id_hlrunfg  = register_tiled_diag_field ( module_name, 'hlrunfg', axes, time, &
             'heat of liq surface runoff', 'W/m2', missing_value=-1.0e+20 )
  id_frunf   = register_tiled_diag_field ( module_name, 'frunf', axes, time, &
             'total rate of solid runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_frunfs  = register_tiled_diag_field ( module_name, 'frunfs', axes, time, &
             'rate of solid calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hfrunf  = register_tiled_diag_field ( module_name, 'hfrunf', axes, time, &
             'heat of total ice runoff', 'W/m2', missing_value=-1.0e+20 )
  id_hfrunfs  = register_tiled_diag_field ( module_name, 'hfrunfs', axes, time, &
             'heat of sol snow runoff', 'W/m2', missing_value=-1.0e+20 )
  id_melt    = register_tiled_diag_field ( module_name, 'melt', axes, time, &
             'total rate of melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltv   = register_tiled_diag_field ( module_name, 'meltv', axes, time, &
             'rate of melt, interception', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_melts   = register_tiled_diag_field ( module_name, 'melts', axes, time, &
             'rate of snow melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltg   = register_tiled_diag_field ( module_name, 'meltg', axes, time, &
             'rate of substrate thaw', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fsw     = register_tiled_diag_field ( module_name, 'fsw', axes, time, &
             'net sw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_fswv    = register_tiled_diag_field ( module_name, 'fswv', axes, time, &
             'net sw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_fsws    = register_tiled_diag_field ( module_name, 'fsws', axes, time, &
             'net sw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_fswg    = register_tiled_diag_field ( module_name, 'fswg', axes, time, &
             'net sw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_flw     = register_tiled_diag_field ( module_name, 'flw', axes, time, &
             'net lw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_flwv    = register_tiled_diag_field ( module_name, 'flwv', axes, time, &
             'net lw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_flws    = register_tiled_diag_field ( module_name, 'flws', axes, time, &
             'net lw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_flwg    = register_tiled_diag_field ( module_name, 'flwg', axes, time, &
             'net lw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_cd_m    = register_tiled_diag_field ( module_name, 'cd_m', axes, time, &
       'drag coefficient for momentum', missing_value=-1e20)
  id_cd_t    = register_tiled_diag_field ( module_name, 'cd_t', axes, time, &
       'drag coefficient for heat and tracers', missing_value=-1e20)
  id_sens    = register_tiled_diag_field ( module_name, 'sens', axes, time, &
             'sens heat flux from land', 'W/m2', missing_value=-1.0e+20 )
  id_sensv   = register_tiled_diag_field ( module_name, 'sensv', axes, time, &
             'sens heat flux from vegn', 'W/m2', missing_value=-1.0e+20 )
  id_senss   = register_tiled_diag_field ( module_name, 'senss', axes, time, &
             'sens heat flux from snow', 'W/m2', missing_value=-1.0e+20 )
  id_sensg   = register_tiled_diag_field ( module_name, 'sensg', axes, time, &
             'sens heat flux from ground', 'W/m2', missing_value=-1.0e+20 )
  id_e_res_1 = register_tiled_diag_field ( module_name, 'e_res_1', axes, time, &
       'canopy air energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_e_res_2 = register_tiled_diag_field ( module_name, 'e_res_2', axes, time, &
       'canopy energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_frac = register_tiled_diag_field(module_name,'frac', axes,&
       time, 'fraction of land area', 'unitless', missing_value=-1.0, op=OP_SUM )
  id_area = register_tiled_diag_field(module_name,'area', axes,&
       time, 'area in the grid cell', 'm2', missing_value=-1.0, op=OP_SUM )
  id_ntiles = register_tiled_diag_field(module_name,'ntiles',axes,  &
       time, 'number of tiles', 'unitless', missing_value=-1.0, op=OP_SUM)
  id_z0m     = register_tiled_diag_field ( module_name, 'z0m', axes, time, &
             'momentum roughness of land', 'm', missing_value=-1.0e+20 )
  id_z0s     = register_tiled_diag_field ( module_name, 'z0s', axes, time, &
             'scalar roughness of land', 'm', missing_value=-1.0e+20 )
  id_con_g_h = register_tiled_diag_field ( module_name, 'con_g_h', axes, time, &
       'conductance for sensible heat between ground surface and canopy air', &
       'm/s', missing_value=-1.0 )
  id_transp  = register_tiled_diag_field ( module_name, 'transp', axes, time, &
             'transpiration; = uptake by roots', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_wroff = register_tiled_diag_field ( module_name, 'wroff', axes, time, &
             'rate of liquid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_sroff = register_tiled_diag_field ( module_name, 'sroff', axes, time, &
             'rate of solid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_htransp = register_tiled_diag_field ( module_name, 'htransp', axes, time, &
             'heat of transpired vapior', 'W/m2', missing_value=-1.0e+20 )
  id_huptake = register_tiled_diag_field ( module_name, 'huptk', axes, time, &
             'heat of soil water uptake', 'W/m2', missing_value=-1.0e+20 )
  id_hroff = register_tiled_diag_field ( module_name, 'hroff', axes, time, &
             'sensible heat of runoff', 'W/m2', missing_value=-1.0e+20 )
  id_gsnow   = register_tiled_diag_field ( module_name, 'gsnow', axes, time, &
             'sens heat into ground from snow', 'W/m2', missing_value=-1.0e+20 )
  call add_tiled_diag_field_alias ( id_gsnow, module_name, 'gflux', axes, time, &
             'obsolete, please use "gsnow" instead', 'W/m2', missing_value=-1.0e+20 )
  id_gequil   = register_tiled_diag_field ( module_name, 'gequil', axes, time, &
             'snow-subs equilibration flux', 'W/m2', missing_value=-1.0e+20 )
  id_grnd_flux = register_tiled_diag_field ( module_name, 'grnd_flux', axes, time, &
             'sensible heat into ground from surface', 'W/m2', missing_value=-1.0e+20 )
  id_soil_water_supply = register_tiled_diag_field ( module_name, 'soil_water_supply', axes, time, &
       'maximum rate of soil water supply to vegetation', 'kg/(m2 s)', missing_value=-1e20)
  id_levapg_max = register_tiled_diag_field ( module_name, 'Eg_max', axes, time, &
             'soil_water limit on vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_water = register_tiled_diag_field ( module_name, 'water', axes, time, &
             'column-integrated soil water', 'kg/m2', missing_value=-1.0e+20 )
  id_snow = register_tiled_diag_field ( module_name, 'snow', axes, time, &
             'column-integrated snow water', 'kg/m2', missing_value=-1.0e+20 )
  id_Trad    = register_tiled_diag_field ( module_name, 'Trad', axes, time, &
             'radiative sfc temperature', 'degK', missing_value=-1.0e+20 )
  id_Tca     = register_tiled_diag_field ( module_name, 'Tca', axes, time, &
             'canopy-air temperature', 'degK', missing_value=-1.0e+20 )
  id_qca     = register_tiled_diag_field ( module_name, 'qca', axes, time, &
             'canopy-air specific humidity', 'kg/kg', missing_value=-1.0 )
  id_qco2    = register_tiled_diag_field ( module_name, 'qco2', axes, time, &
             'canopy-air CO2 moist mass mixing ratio', 'kg/kg', missing_value=-1.0 )
  id_qco2_dvmr = register_tiled_diag_field ( module_name, 'qco2_dvmr', axes, time, &
             'canopy-air CO2 dry volumetric mixing ratio', 'mol CO2/mol air', missing_value=-1.0 )
  id_fco2    = register_tiled_diag_field ( module_name, 'fco2', axes, time, &
             'flux of CO2 to canopy air', 'kg C/(m2 s)', missing_value=-1.0 )
  id_swdn_dir = register_tiled_diag_field ( module_name, 'swdn_dir', (/id_lon,id_lat,id_band/), time, &
       'downward direct short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_swdn_dif = register_tiled_diag_field ( module_name, 'swdn_dif', (/id_lon,id_lat,id_band/), time, &
       'downward diffuse short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_swup_dir = register_tiled_diag_field ( module_name, 'swup_dir', (/id_lon,id_lat,id_band/), time, &
       'direct short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_swup_dif = register_tiled_diag_field ( module_name, 'swup_dif', (/id_lon,id_lat,id_band/), time, &
       'diffuse short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_lwdn = register_tiled_diag_field ( module_name, 'lwdn', axes, time, &
       'downward long-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_vegn_cover = register_tiled_diag_field ( module_name, 'vegn_cover', axes, time, &
             'fraction covered by vegetation', missing_value=-1.0 )
  id_cosz = register_tiled_diag_field ( module_name, 'coszen', axes, time, &
       'cosine of zenith angle', missing_value=-2.0 )
  id_albedo_dir = register_tiled_diag_field ( module_name, 'albedo_dir', &
       (/id_lon,id_lat,id_band/), time, &
       'land surface albedo for direct light', missing_value=-1.0 )
  id_albedo_dif = register_tiled_diag_field ( module_name, 'albedo_dif', &
       (/id_lon,id_lat,id_band/), time, &
       'land surface albedo for diffuse light', missing_value=-1.0 )
  id_vegn_refl_dir = register_tiled_diag_field(module_name, 'vegn_refl_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'black-background canopy reflectivity for direct light',missing_value=-1.0)
  id_vegn_refl_dif = register_tiled_diag_field(module_name, 'vegn_refl_dif', &
       (/id_lon, id_lat, id_band/), time, &
       'black-background canopy reflectivity for diffuse light',missing_value=-1.0)
  id_vegn_refl_lw = register_tiled_diag_field ( module_name, 'vegn_refl_lw', axes, time, &
       'canopy reflectivity for thermal radiation', missing_value=-1.0)
  id_vegn_tran_dir = register_tiled_diag_field(module_name, 'vegn_tran_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'part of direct light that passes through canopy unscattered',missing_value=-1.0)
  id_vegn_tran_dif = register_tiled_diag_field(module_name, 'vegn_tran_dif', &
       (/id_lon, id_lat, id_band/), time, &
       'black-background canopy transmittance for diffuse light',missing_value=-1.0)
  id_vegn_tran_lw = register_tiled_diag_field ( module_name, 'vegn_tran_lw', axes, time, &
       'canopy transmittance for thermal radiation', missing_value=-1.0)
  id_vegn_sctr_dir = register_tiled_diag_field(module_name, 'vegn_sctr_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'part of direct light scattered downward by canopy',missing_value=-1.0)
  id_subs_refl_dir = register_tiled_diag_field(module_name, 'subs_refl_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'substrate reflectivity for direct light',missing_value=-1.0)
  id_subs_refl_dif = register_tiled_diag_field(module_name, 'subs_refl_dif', &
       (/id_lon, id_lat, id_band/), time, &
       'substrate reflectivity for diffuse light',missing_value=-1.0)
  id_subs_emis = register_tiled_diag_field(module_name, 'subs_emis', &
       (/id_lon, id_lat/), time, &
       'substrate emissivity for long-wave radiation',missing_value=-1.0)
  id_grnd_T = register_tiled_diag_field ( module_name, 'Tgrnd', axes, time, &
       'ground surface temperature', 'degK', missing_value=-1.0 )
  id_total_C = register_tiled_diag_field ( module_name, 'Ctot', axes, time, &
       'total land carbon', 'kg C/m2', missing_value=-1.0 )
end subroutine land_diag_init

! the code below defines the accessor routines that are used to access fields of the 
! tile data structure in collective operations, like restart i/o. Fore example, a statement
! DEFINE_LAND_ACCESSOR_0D(real,lwup)
! defines a function "land_lwup_ptr" that, given a tile, returns a pointer to the field
! called "lwup" in this tile. The procedure implementing a collective operation would
! enumerate all the tiles within the domain, call accessor routine for each of them, and
! get or set the value pointed to by the accessor routine.
#define DEFINE_LAND_ACCESSOR_0D(xtype,x) subroutine CONCAT3(land_,x,_ptr(t,p));\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))p=>t%x;end subroutine

DEFINE_LAND_ACCESSOR_0D(real,frac)
DEFINE_LAND_ACCESSOR_0D(real,lwup)
DEFINE_LAND_ACCESSOR_0D(real,e_res_1)
DEFINE_LAND_ACCESSOR_0D(real,e_res_2)

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function land_tile_exists(tile)
  type(land_tile_type), pointer :: tile
  land_tile_exists = associated(tile)
end function land_tile_exists

#define DEFINE_TAG_ACCESSOR(x) subroutine  CONCAT2(x,_tag_ptr(t,p));\
type(land_tile_type),pointer::t;integer,pointer::p;p=>NULL();if(associated(t))\
then;if (associated(t%x)) p=>t%x%tag;endif;end subroutine

DEFINE_TAG_ACCESSOR(glac)
DEFINE_TAG_ACCESSOR(lake)
DEFINE_TAG_ACCESSOR(soil)
DEFINE_TAG_ACCESSOR(vegn)

end module land_model_mod
   
