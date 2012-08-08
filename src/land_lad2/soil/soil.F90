! ============================================================================
! soil model module
! ============================================================================
module soil_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: error_mesg, file_exist, check_nml_error, &
     stdlog, write_version_number, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use mpp_io_mod,         only: mpp_open, MPP_RDONLY
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o, PI
use horiz_interp_mod,   only: horiz_interp

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR
use soil_tile_mod, only : &
     soil_tile_type, soil_pars_type, soil_prog_type, read_soil_data_namelist, &
     soil_data_radiation, soil_data_diffusion, soil_data_thermodynamics, &
     soil_data_hydraulics, soil_data_gw_hydraulics, & ! soil_data_gw_tables, &
     soil_data_vwc_sat, &
     max_lev, psi_wilt, cpw, clw, csw, g_iso, g_vol, g_geo, g_RT, &
     num_storage_pts, num_zeta_s_pts, gw_zeta_s, gw_flux_table, gw_area_table, &
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth

use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, get_elmt_indices, &
     operator(/=)
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr
use land_data_mod,      only : land_state_type, lnd
use land_io_mod, only : read_field
use land_tile_io_mod, only : create_tile_out_file, write_tile_data_r0d_fptr,& 
     write_tile_data_r1d_fptr,read_tile_data_r0d_fptr, read_tile_data_r1d_fptr,&
     print_netcdf_error, get_input_restart_name, sync_nc_files
use nf_utils_mod, only : nfu_def_dim, nfu_put_att, nfu_inq_var
use vegn_tile_mod, only : vegn_tile_type, vegn_uptake_profile, vegn_root_properties
use land_debug_mod, only : is_watch_point, get_current_point
use uptake_mod, only : UPTAKE_LINEAR, UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN, &
     uptake_init, &
     darcy2d_uptake, darcy2d_uptake_solver, &
     darcy2d_uptake_lin, darcy2d_uptake_solver_lin
implicit none
private

! ==== public interfaces =====================================================
public :: read_soil_namelist
public :: soil_init
public :: soil_end
public :: save_soil_restart

public :: soil_get_sfc_temp
public :: soil_radiation
public :: soil_diffusion
public :: soil_step_1
public :: soil_step_2
! =====end of public interfaces ==============================================



! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'soil',&
    version     = '$Id: soil.F90,v 17.0.2.2.2.2 2011/12/16 19:01:57 pjp Exp $',&
    tagname     = '$Name: siena_201207 $'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_E_min            = .false.     ! prohibit condensation
logical :: use_E_max            = .true.      ! theoretical effiltration capacity flag
real    :: init_temp            = 288.        ! cold-start soil T
real    :: init_w               = 150.        ! cold-start w(l)/dz(l)
real    :: init_groundwater     =   0.        ! cold-start gw storage
real    :: lrunf_ie_min         = -1.0e-4     ! trigger for clip and runoff
real    :: lrunf_ie_tol         =  1.e-12
character(len=16) :: albedo_to_use = ''       ! or 'albedo-map' or 'brdf-maps'
character(len=24) :: uptake_to_use = 'linear' ! or 'darcy2d', or 'darcy2d-linearized'
logical :: uptake_oneway        = .false.     ! if true, roots can't loose water to soil
logical :: uptake_from_sat      = .true.      ! if false, the uptake from saturated soil is prohibited
logical :: unconditional_sweep  = .false.
logical :: allow_negative_rie   = .false.
logical :: baseflow_where_frozen = .false.
logical :: write_when_flagged   = .false.
logical :: bypass_richards_when_stiff = .true.
logical :: corrected_lm2_gw     = .true.
real    :: active_layer_drainage_acceleration = 0.

namelist /soil_nml/ lm2, use_E_min, use_E_max,           &
                    init_temp,      &
                    init_w,       &
                    init_groundwater, lrunf_ie_min, lrunf_ie_tol, &
                    cpw, clw, csw, &
                    albedo_to_use, &
                    uptake_to_use, uptake_oneway, uptake_from_sat, &
                    unconditional_sweep, allow_negative_rie, &
                    baseflow_where_frozen, &
                    write_when_flagged, &
                    bypass_richards_when_stiff, corrected_lm2_gw, &
                    active_layer_drainage_acceleration
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf = .false.
type(time_type) :: time
real            :: delta_time
logical         :: use_single_geo, use_geohydrology
integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)
real            :: Eg_min

integer         :: uptake_option = -1 

! ---- diagnostic field IDs
integer :: id_lwc, id_swc, id_psi, id_temp, &
    id_ie, id_sn, id_bf, id_nu, id_hie, id_hsn, id_hbf, id_hnu, &
    id_heat_cap, id_thermal_cond, id_type, id_tau_gw, id_slope_l, &
    id_slope_Z, id_zeta_bar, id_e_depth, id_vwc_sat, id_vwc_fc, &
    id_vwc_wilt, id_K_sat, id_w_fc, &
    id_refl_dry_dif, id_refl_dry_dir, id_refl_sat_dif, id_refl_sat_dir, &
    id_f_iso_dry, id_f_vol_dry, id_f_geo_dry, &
    id_f_iso_sat, id_f_vol_sat, id_f_geo_sat, &
    id_evap, id_uptk_n_iter, id_uptk, id_uptk_residual, id_excess, &
    id_psi_bot, id_sat_frac, id_stor_frac, id_sat_depth, &
    id_uptk_old, id_psi_bot_old, id_sat_frac_old, id_stor_frac_old, &
    id_sat_depth_old, id_slope_Z_old, id_e_depth_old, &
    id_vwc_wilt_old, id_vwc_fc_old, id_vwc_sat_old, id_K_sat_old

! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_soil_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_soil_data_namelist(num_l,dz,use_single_geo,use_geohydrology)

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=soil_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=soil_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;   
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

  ! ---- convert symbolic names of the uptake options into numeric IDs to speed up
  ! the selection during run-time
  if (trim(uptake_to_use)=='linear') then
     uptake_option = UPTAKE_LINEAR
  else if (trim(uptake_to_use)=='darcy2d') then
     uptake_option = UPTAKE_DARCY2D
  else if (trim(uptake_to_use)=='darcy2d-linearized') then
     uptake_option = UPTAKE_DARCY2D_LIN
  else 
     call error_mesg('soil_init',&
          'soil uptake option uptake_to_use="'//&
          trim(uptake_to_use)//'" is invalid, use "linear", "darcy2d" or "darcy2d-linearized"',&
          FATAL)
  endif

  if (use_E_min) then
      Eg_min = 0.
    else
      Eg_min = -HUGE(Eg_min)
    endif

end subroutine read_soil_namelist


! ============================================================================
! initialize soil model
subroutine soil_init ( id_lon, id_lat, id_band )
  integer, intent(in)  :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in)  :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in)  :: id_band ! ID of spectral band axis

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce  ! tail and current tile list elements
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  real, allocatable :: gw_param(:,:), gw_param2(:,:), albedo(:,:,:) ! input data buffers for respective variables
  real, allocatable :: f_iso(:,:,:), f_vol(:,:,:), f_geo(:,:,:), refl_dif(:,:,:)

  integer :: i, code, m
  real :: zeta_s, frac
 character(len=256) :: restart_file_name
  logical :: restart_exists

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  call uptake_init(num_l,dz,zfull)

  ! -------- initialize soil state --------
  te = tail_elmt (lnd%tile_map)
  ce = first_elmt(lnd%tile_map)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)        ! advance position to the next tile
     
     if (.not.associated(tile%soil)) cycle
     
     if (init_temp.ge.tile%soil%pars%tfreeze) then
        tile%soil%prog(1:num_l)%wl = init_w*dz(1:num_l)
        tile%soil%prog(1:num_l)%ws = 0
     else
        tile%soil%prog(1:num_l)%wl = 0
        tile%soil%prog(1:num_l)%ws = init_w*dz(1:num_l)
     endif
     tile%soil%prog%T             = init_temp
     tile%soil%prog%groundwater   = init_groundwater
     tile%soil%prog%groundwater_T = init_temp
     
     tile%soil%uptake_T           = init_temp
  enddo

  call get_input_restart_name('INPUT/soil.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('soil_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call read_tile_data_r1d_fptr(unit, 'temp'         , soil_T_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'           , soil_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'           , soil_ws_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater'  , soil_groundwater_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater_T', soil_groundwater_T_ptr)
     if(nfu_inq_var(unit, 'uptake_T')==NF_NOERR) &
          call read_tile_data_r0d_fptr(unit, 'uptake_T', soil_uptake_T_ptr)
          
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('soil_init',&
          'cold-starting soil',&
          NOTE)
  endif
  
  ! initialize soil model diagnostic fields
  call soil_diag_init ( id_lon, id_lat, id_band )
  
  ! read groundwater parameters, if requested
  if (.not.use_single_geo) then
      if (.not.use_geohydrology) then
          allocate(gw_param(lnd%is:lnd%ie,lnd%js:lnd%je))
          call read_field( 'INPUT/groundwater_residence.nc','tau', &
               lnd%lon, lnd%lat, gw_param, interp='bilinear' )
          call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_tau_groundwater_ptr )
          deallocate(gw_param)
        else
          allocate(gw_param (lnd%is:lnd%ie,lnd%js:lnd%je))
          allocate(gw_param2(lnd%is:lnd%ie,lnd%js:lnd%je))
          call read_field( 'INPUT/geohydrology.nc','hillslope_length', &
               lnd%lon, lnd%lat, gw_param, interp='bilinear' )
          call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, lnd%tile_map, soil_hillslope_length_ptr )
          call read_field( 'INPUT/geohydrology.nc','slope', &
               lnd%lon, lnd%lat, gw_param2, interp='bilinear' )
          gw_param = gw_param*gw_param2
          call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, lnd%tile_map, soil_hillslope_relief_ptr )
          call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
               lnd%lon, lnd%lat, gw_param, interp='bilinear' )
          call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_zeta_bar_ptr )
          call read_field( 'INPUT/geohydrology.nc','soil_e_depth', &
               lnd%lon, lnd%lat, gw_param, interp='bilinear' )
          call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, lnd%tile_map, soil_soil_e_depth_ptr )
          deallocate(gw_param, gw_param2)
          te = tail_elmt (lnd%tile_map)
          ce = first_elmt(lnd%tile_map)
          do while(ce /= te)
            tile=>current_tile(ce)  ! get pointer to current tile
            ce=next_elmt(ce)        ! advance position to the next tile
            if (.not.associated(tile%soil)) cycle
            if (tile%soil%pars%hillslope_relief.le.0.) &
                tile%soil%pars%hillslope_relief =      &
                   tile%soil%pars%soil_e_depth / gw_zeta_s(num_zeta_s_pts)
            zeta_s = tile%soil%pars%soil_e_depth / tile%soil%pars%hillslope_relief
            zeta_s = max(zeta_s, gw_zeta_s(1))
            zeta_s = min(zeta_s, gw_zeta_s(num_zeta_s_pts))
            m = num_zeta_s_pts / 2
            code = 0
            do while (code.eq.0)
              if (zeta_s .lt. gw_zeta_s(m)) then
                  m = m - 1
                else if (zeta_s .gt. gw_zeta_s(m+1)) then
                  m = m + 1
                else
                 code = 1
                endif
              enddo
            frac = (zeta_s - gw_zeta_s(m)) / (gw_zeta_s(m+1) - gw_zeta_s(m))
            do i = 1, num_storage_pts
              tile%soil%pars%gw_flux_norm(i) = gw_flux_table(i,m) &
                   + frac*(gw_flux_table(i,m+1)-gw_flux_table(i,m))
              tile%soil%pars%gw_area_norm(i) = gw_area_table(i,m) &
                   + frac*(gw_area_table(i,m+1)-gw_area_table(i,m))
               enddo
            enddo
        endif
    endif

  ! set dry soil albedo values, if requested
  if (trim(albedo_to_use)=='albedo-map') then
     allocate(albedo(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_VIS',&
          lnd%lon, lnd%lat, albedo(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_NIR',&
          lnd%lon, lnd%lat, albedo(:,:,BAND_NIR),'bilinear')
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_dry_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo doesn't depend on soil wetness
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_sat_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_sat_dif_ptr )
     deallocate(albedo)
  else if (trim(albedo_to_use)=='brdf-maps') then
     use_brdf = .true.
     allocate(   f_iso(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     allocate(   f_vol(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     allocate(   f_geo(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     allocate(refl_dif(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
! *********************** ????????? *************
! sergey-- these are hig-res maps. is 'bilinear' the best option to use? i simply
! copied it from the albedo-map option.
! *********************** ????????? *************
     call read_field( 'INPUT/soil_brdf.nc','f_iso_vis',&
          lnd%lon, lnd%lat, f_iso(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_vis',&
          lnd%lon, lnd%lat, f_vol(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_vis',&
          lnd%lon, lnd%lat, f_geo(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_iso_nir',&
          lnd%lon, lnd%lat, f_iso(:,:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_nir',&
          lnd%lon, lnd%lat, f_vol(:,:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_nir',&
          lnd%lon, lnd%lat, f_geo(:,:,BAND_NIR),'bilinear')
     refl_dif = g_iso*f_iso + g_vol*f_vol + g_geo*f_geo
     call put_to_tiles_r1d_fptr( f_iso,    lnd%tile_map, soil_f_iso_dry_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    lnd%tile_map, soil_f_vol_dry_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    lnd%tile_map, soil_f_geo_dry_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, lnd%tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo doesn't depend on soil wetness
     call put_to_tiles_r1d_fptr( f_iso,    lnd%tile_map, soil_f_iso_sat_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    lnd%tile_map, soil_f_vol_sat_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    lnd%tile_map, soil_f_geo_sat_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, lnd%tile_map, soil_refl_sat_dif_ptr )
     deallocate(f_iso, f_vol, f_geo, refl_dif)
  else if (trim(albedo_to_use)=='') then
     ! do nothing, that is leave soil albedo parameters as defined based on the data table
  else
     call error_mesg('soil_init',&
          'option albedo_to_use="'// trim(albedo_to_use)//&
          '" is invalid, use "albedo-map", "brdf-maps", or nothing ("")',&
          FATAL)
  endif
  
  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_tau_gw,       lnd%tile_map, soil_tau_groundwater_ptr)
  call send_tile_data_r0d_fptr(id_slope_l,      lnd%tile_map, soil_hillslope_length_ptr)
  call send_tile_data_r0d_fptr(id_slope_Z,      lnd%tile_map, soil_hillslope_relief_ptr)
  call send_tile_data_r0d_fptr(id_zeta_bar,     lnd%tile_map, soil_hillslope_zeta_bar_ptr)
  call send_tile_data_r0d_fptr(id_e_depth,      lnd%tile_map, soil_soil_e_depth_ptr)
  call send_tile_data_r0d_fptr(id_vwc_wilt,     lnd%tile_map, soil_vwc_wilt_ptr)
  call send_tile_data_r0d_fptr(id_vwc_fc,       lnd%tile_map, soil_vwc_fc_ptr)
  call send_tile_data_r0d_fptr(id_vwc_sat,      lnd%tile_map, soil_vwc_sat_ptr)
  call send_tile_data_r0d_fptr(id_K_sat,        lnd%tile_map, soil_k_sat_ref_ptr)
  call send_tile_data_r1d_fptr(id_w_fc,         lnd%tile_map, soil_w_fc_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dir, lnd%tile_map, soil_refl_dry_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dif, lnd%tile_map, soil_refl_dry_dif_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dir, lnd%tile_map, soil_refl_sat_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dif, lnd%tile_map, soil_refl_sat_dif_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_dry, lnd%tile_map, soil_f_iso_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_dry, lnd%tile_map, soil_f_vol_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_dry, lnd%tile_map, soil_f_geo_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_sat, lnd%tile_map, soil_f_iso_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_sat, lnd%tile_map, soil_f_vol_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_sat, lnd%tile_map, soil_f_geo_sat_ptr)
  call send_tile_data_i0d_fptr(id_type,         lnd%tile_map, soil_tag_ptr)

  call send_tile_data_r0d_fptr(id_slope_Z_old,      lnd%tile_map, soil_hillslope_relief_ptr)
  call send_tile_data_r0d_fptr(id_e_depth_old,      lnd%tile_map, soil_soil_e_depth_ptr)
  call send_tile_data_r0d_fptr(id_vwc_wilt_old,     lnd%tile_map, soil_vwc_wilt_ptr)
  call send_tile_data_r0d_fptr(id_vwc_fc_old,       lnd%tile_map, soil_vwc_fc_ptr)
  call send_tile_data_r0d_fptr(id_vwc_sat_old,      lnd%tile_map, soil_vwc_sat_ptr)
  call send_tile_data_r0d_fptr(id_K_sat_old,        lnd%tile_map, soil_k_sat_ref_ptr)

end subroutine soil_init


! ============================================================================
subroutine soil_diag_init ( id_lon, id_lat, id_band )
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in) :: id_band ! ID of spectral band axis  

  ! ---- local vars
  integer :: axes(3)
  integer :: id_zhalf, id_zfull

  ! define vertical axis and its' edges
  id_zhalf = diag_axis_init ( &
       'zhalf_soil', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='soil' )
  id_zfull = diag_axis_init ( &
       'zfull_soil', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='soil', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define diagnostic fields
  id_lwc = register_tiled_diag_field ( module_name, 'soil_liq', axes,  &
       Time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'soil_ice',  axes,  &
       Time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_psi = register_tiled_diag_field ( module_name, 'soil_psi', axes,  &
       Time, 'soil-water matric head', 'm', missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'soil_T',  axes,       &
       Time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'soil_rie',  axes(1:2),  &
       Time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'soil_rsn',  axes(1:2),  &
       Time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'soil_rbf',  axes(1:2),  &
       Time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_nu  = register_tiled_diag_field ( module_name, 'soil_rnu',  axes(1:2),  &
       Time, 'numerical runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'soil_hie',  axes(1:2), &
       Time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'soil_hsn',  axes(1:2), &
       Time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'soil_hbf',  axes(1:2), &
       Time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )
  id_hnu  = register_tiled_diag_field ( module_name, 'soil_hnu',  axes(1:2), &
       Time, 'heat nu runoff',          'W/m2',  missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'soil_evap',  axes(1:2), &
       Time, 'soil evap',            'kg/(m2 s)',  missing_value=-100.0 )
  id_excess  = register_tiled_diag_field ( module_name, 'sfc_excess',  axes(1:2),  &
       Time, 'sfc excess pushed down',    'kg/(m2 s)',  missing_value=-100.0 )

  id_uptk_n_iter  = register_tiled_diag_field ( module_name, 'uptake_n_iter',  axes(1:2), &
       Time, 'number of iterations for soil uptake',  missing_value=-100.0 )
  id_uptk = register_tiled_diag_field ( module_name, 'soil_uptk', axes, &
       Time, 'uptake of water by roots', 'kg/(m2 s)',  missing_value=-100.0 )
  id_psi_bot = register_tiled_diag_field ( module_name, 'soil_psi_n', axes(1:2), &
       Time, 'psi at bottom of soil column', 'm',  missing_value=-100.0 )
  id_sat_frac = register_tiled_diag_field ( module_name, 'soil_fsat', axes(1:2), &
       Time, 'fraction of soil area saturated at surface', '-',  missing_value=-100.0 )
  id_stor_frac = register_tiled_diag_field ( module_name, 'soil_fgw', axes(1:2), &
       Time, 'groundwater storage frac above base elev', '-',  missing_value=-100.0 )
  id_sat_depth = register_tiled_diag_field ( module_name, 'soil_wtdep', axes(1:2), &
       Time, 'depth below sfc to saturated soil', 'm',  missing_value=-100.0 )

  ! ---- the following fields are for compatibility with older diag tables
  id_uptk_old = register_tiled_diag_field ( module_name, 'uptake', axes, &
       Time, 'uptake of water by roots (obsolete, use "soil_uptk" instead)', &
       'kg/(m2 s)',  missing_value=-100.0 )
  id_psi_bot_old = register_tiled_diag_field ( module_name, 'psi_bot', axes(1:2), &
       Time, 'psi at bottom of soil column (obsolete, use "soil_psi_n" instead)', &
       'm',  missing_value=-100.0 )
  id_sat_frac_old = register_tiled_diag_field ( module_name, 'sat_frac', axes(1:2), &
       Time, 'fraction of soil area saturated at surface (obsolete, use "soil_fsat" instead)',&
       '-',  missing_value=-100.0 )
  id_stor_frac_old = register_tiled_diag_field ( module_name, 'stor_frac', axes(1:2), &
       Time, 'groundwater storage frac above base elev (obsolete, use "soil_fgw" instead)',&
       '-',  missing_value=-100.0 )
  id_sat_depth_old = register_tiled_diag_field ( module_name, 'sat_depth', axes(1:2), &
       Time, 'depth below sfc to saturated soil (obsolete, use "soil_wtdep" instead)', &
       'm',  missing_value=-100.0 )
  ! ---- end of compatibility section

  id_heat_cap = register_tiled_diag_field ( module_name, 'soil_heat_cap',  &
       axes, Time, 'heat capacity of dry soil','J/(m3 K)', missing_value=-100.0 )
  id_thermal_cond =  register_tiled_diag_field ( module_name, 'soil_tcon', &
       axes, Time, 'soil thermal conductivity', 'W/(m K)',  missing_value=-100.0 )
  
  id_type = register_tiled_static_field ( module_name, 'soil_type',  &
       axes(1:2), 'soil type', missing_value=-1.0 )
  id_tau_gw = register_tiled_static_field ( module_name, 'tau_gw',  &
       axes(1:2), 'groundwater residence time', 's', missing_value=-100.0 )
  id_slope_l = register_tiled_static_field ( module_name, 'slope_l',  &
       axes(1:2), 'hillslope length', 'm', missing_value=-100.0 )
  id_slope_Z = register_tiled_static_field ( module_name, 'soil_rlief',  &
       axes(1:2), 'hillslope relief', 'm', missing_value=-100.0 )
  id_zeta_bar = register_tiled_static_field ( module_name, 'zeta_bar',  &
       axes(1:2), 'hillslope zeta bar', '-', missing_value=-100.0 )
  id_e_depth = register_tiled_static_field ( module_name, 'soil_depth',  &
       axes(1:2), 'soil e-folding depth', 'm', missing_value=-100.0 )
  id_vwc_wilt = register_tiled_static_field ( module_name, 'soil_wilt',  &
       axes(1:2), 'wilting water content', '-', missing_value=-100.0 )
  id_vwc_fc = register_tiled_static_field ( module_name, 'soil_fc',  &
       axes(1:2), 'field capacity', '-', missing_value=-100.0 )
  id_vwc_sat = register_tiled_static_field ( module_name, 'soil_sat',  &
       axes(1:2), 'soil porosity', '-', missing_value=-100.0 )
  id_K_sat = register_tiled_static_field ( module_name, 'soil_Ksat',  &
       axes(1:2), 'soil sat. hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_w_fc = register_tiled_static_field ( module_name, 'w_fc',  &
       axes, 'soil field capacity', missing_value=-1.0 )
  id_refl_dry_dir = register_tiled_static_field ( module_name, 'refl_dry_dir',  &
       (/id_lon, id_lat, id_band/), 'reflectance of dry soil for direct light', &
       missing_value=-1.0 )
  id_refl_dry_dif = register_tiled_static_field ( module_name, 'refl_dry_dif',  &
       (/id_lon, id_lat, id_band/), 'reflectance of dry soil for diffuse light', &
       missing_value=-1.0 )
  id_refl_sat_dir = register_tiled_static_field ( module_name, 'refl_sat_dir',  &
       (/id_lon, id_lat, id_band/), 'reflectance of saturated soil for direct light', &
       missing_value=-1.0 )
  id_refl_sat_dif = register_tiled_static_field ( module_name, 'refl_sat_dif',  &
       (/id_lon, id_lat, id_band/), 'reflectance of saturated soil for diffuse light', &
       missing_value=-1.0 )
  id_f_iso_dry = register_tiled_static_field ( module_name, 'f_iso_dry',  &
       (/id_lon, id_lat, id_band/), 'isotropic brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_vol_dry = register_tiled_static_field ( module_name, 'f_vol_dry',  &
       (/id_lon, id_lat, id_band/), 'volumetric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_geo_dry = register_tiled_static_field ( module_name, 'f_geo_dry',  &
       (/id_lon, id_lat, id_band/), 'geometric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_iso_sat = register_tiled_static_field ( module_name, 'f_iso_sat',  &
       (/id_lon, id_lat, id_band/), 'isotropic brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_vol_sat = register_tiled_static_field ( module_name, 'f_vol_sat',  &
       (/id_lon, id_lat, id_band/), 'volumetric brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_geo_sat = register_tiled_static_field ( module_name, 'f_geo_sat',  &
       (/id_lon, id_lat, id_band/), 'geometric brdf weight, saturated soil', &
       missing_value=-1.0 )

  ! the following fields are for compatibility with older diag tables only
  id_slope_Z_old = register_tiled_static_field ( module_name, 'slope_Z',  &
       axes(1:2), 'hillslope relief (obsolete, use "soil_rlief" instead)',&
       'm', missing_value=-100.0 )
  id_e_depth_old = register_tiled_static_field ( module_name, 'e_depth',  &
       axes(1:2), 'soil e-folding depth (obsolete, use "soil_depth" instead)', &
       'm', missing_value=-100.0 )
  id_vwc_wilt_old = register_tiled_static_field ( module_name, 'vwc_wilt',  &
       axes(1:2), 'wilting water content (obsolete, use "soil_wilt" instead)', &
       '-', missing_value=-100.0 )
  id_vwc_fc_old = register_tiled_static_field ( module_name, 'vwc_fc',  &
       axes(1:2), 'field capacity (obsolete, use "soil_fc" instead)', &
       '-', missing_value=-100.0 )
  id_vwc_sat_old = register_tiled_static_field ( module_name, 'vwc_sat',  &
       axes(1:2), 'soil porosity (obsolete, use "soil_sat")', &
       '-', missing_value=-100.0 )
  id_K_sat_old = register_tiled_static_field ( module_name, 'K_sat',  &
       axes(1:2), 'soil sat. hydraulic conductivity (obsolte, use "soil_Ksat" instead)', &
       'kg /(m2 s)', missing_value=-100.0 )

end subroutine soil_diag_init


! ============================================================================
subroutine soil_end ()

  module_is_initialized =.FALSE.

end subroutine soil_end


! ============================================================================
subroutine save_soil_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit

  call error_mesg('soil_end','writing NetCDF restart',NOTE)
  ! create output file, including internal structure necessary for output
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'soil.res.nc', &
          lnd%coord_glon, lnd%coord_glat, soil_tile_exists, tile_dim_length )
  ! in addition, define vertical coordinate
  if (mpp_pe()==lnd%io_pelist(1)) then
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
     __NF_ASRT__(nfu_put_att(unit,'zfull','positive','down'))
  endif
  call sync_nc_files(unit)
        
  ! write out fields
  call write_tile_data_r1d_fptr(unit,'temp'         ,soil_T_ptr   ,'zfull','soil temperature','degrees_K')
  call write_tile_data_r1d_fptr(unit,'wl'           ,soil_wl_ptr  ,'zfull','liquid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'ws'           ,soil_ws_ptr  ,'zfull','solid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'groundwater'  ,soil_groundwater_ptr  ,'zfull')
  call write_tile_data_r1d_fptr(unit,'groundwater_T',soil_groundwater_T_ptr ,'zfull')
  call write_tile_data_r0d_fptr(unit,'uptake_T',     soil_uptake_T_ptr, 'temperature of transpiring water', 'degrees_K')
  
  ! close file
  __NF_ASRT__(nf_close(unit))

end subroutine save_soil_restart


! ============================================================================
subroutine soil_get_sfc_temp ( soil, soil_T )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_T

  soil_T= soil%prog(1)%T
end subroutine soil_get_sfc_temp


! ============================================================================
! compute soil radiative properties
subroutine soil_radiation ( soil, cosz, &
     soil_refl_dir, soil_refl_dif, soil_refl_lw, soil_emis )
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)  :: cosz
  real, intent(out) :: soil_refl_dir(NBANDS), soil_refl_dif(NBANDS), soil_refl_lw, soil_emis

  call soil_data_radiation ( soil, cosz, use_brdf, soil_refl_dir, soil_refl_dif, soil_emis )
  soil_refl_lw = 1 - soil_emis
  if(any(soil_refl_dif<0).or.any(soil_refl_dif>1).or.&
     any(soil_refl_dir<0).or.any(soil_refl_dir>1)) then
    write(*,*)'soil_refl is out of range'
    write(*,*)'soil_refl_dif=',soil_refl_dif
    write(*,*)'soil_refl_dir=',soil_refl_dir
  endif
end subroutine soil_radiation


! ============================================================================
! compute soil roughness
subroutine soil_diffusion ( soil, soil_z0s, soil_z0m )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_z0s, soil_z0m

  call soil_data_diffusion ( soil, soil_z0s, soil_z0m )
end subroutine soil_diffusion


! ============================================================================
! compute beta function
! after Manabe (1969), but distributed vertically.
subroutine soil_data_beta ( soil, vegn, soil_beta, soil_water_supply, &
                            soil_uptake_T, soil_rh, soil_rh_psi )
  type(soil_tile_type), intent(inout)  :: soil
  type(vegn_tile_type), intent(in)     :: vegn
  real, intent(out) :: soil_beta
  real, intent(out) :: soil_water_supply ! max rate of water supply to roots, kg/(m2 s)
  real, intent(out) :: soil_uptake_T ! an estimate of temperature of the water 
             ! taken up by transpiration. In case of 'linear' uptake it is an exact
             ! value; in case of 'darcy*' treatments the actual uptake profile
             ! is calculated only in step 2, so the value returned is an estimate  
  real, intent(out) :: soil_rh
  real, intent(out) :: soil_rh_psi

  ! ---- local vars
  integer :: l
  real, dimension(num_l) :: &
       uptake_frac_max, & ! root distribution
       vegn_uptake_term, &
       vlc, vsc, & ! volumetric fractions of water and ice in the layer
       DThDP, hyd_cond, DKDP, soil_w_fc, & ! soil hydraulic parameters (not used)
       VRL, & ! vertical distribution of volumetric root length, m/m3
       u, du ! uptake and its derivative (the latter is not used)
  real :: DPsi_min, DPsi_max, tau_gw, psi_for_rh
  real :: gw_length, gw_relief, gw_zeta_bar, gw_e_depth, K_sat ! soil hydraulic parameters (not used)
  real :: K_r, r_r ! root properties
  real :: z  !  soil depth

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

  vlc=0;vsc=0
  do l = 1, num_l
    vlc(l) = max(0., soil%prog(l)%wl / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%prog(l)%ws / (dens_h2o*dz(l)))
    enddo
  
  soil%uptake_frac = 0
  do l = 1, num_l
     soil%uptake_frac(l) = uptake_frac_max(l) &
          * max(0.0, min(1.0,(vlc(l)-soil%w_wilt(l))/&
               (0.75*(soil%w_fc(l)-soil%w_wilt(l)))))
  enddo
  soil_beta = sum(soil%uptake_frac)
  do l = 1, num_l
     if (soil_beta /= 0) then
          soil%uptake_frac(l) = soil%uptake_frac(l) / soil_beta
     else
          soil%uptake_frac(l) = uptake_frac_max(l)
     endif
  enddo
  if (lm2) soil%uptake_frac = uptake_frac_max

  ! calculate soil hydraulic properties, in particular psi_for_rh -- we don't use 
  ! anything else in this subroutine. this moved out of 'case' because
  ! we need psi unconditionally now for soil_rh
  
  call soil_data_hydraulics ( soil, vlc, vsc, &
       soil%psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, tau_gw, &
       psi_for_rh, soil_w_fc )
  soil_rh = exp(psi_for_rh*g_RT)
  soil_rh_psi = g_RT*soil_rh

  ! calculate total water supply
  select case (uptake_option)
  case(UPTAKE_LINEAR)
     soil_water_supply = 0
     z = 0
     do l = 1, num_l
        soil_water_supply = soil_water_supply + &
          vegn_uptake_term(l)*max(0.0,soil%prog(l)%wl/dz(l)-soil%w_wilt(l)*dens_h2o)
        z = z + dz(l)
     enddo
     soil_water_supply = z * soil_water_supply
     soil_water_supply = soil_water_supply/delta_time
     soil_uptake_T = sum(soil%uptake_frac*soil%prog%T)
  case(UPTAKE_DARCY2D)
     call vegn_root_properties (vegn, dz(1:num_l), VRL, K_r, r_r)
     call darcy2d_uptake ( soil, psi_wilt, VRL, K_r, r_r, uptake_oneway,&
          uptake_from_sat, u, du )
     soil_water_supply = max(0.0,sum(u))
     soil_uptake_T = soil%uptake_T
  case(UPTAKE_DARCY2D_LIN)
     call vegn_root_properties (vegn, dz(1:num_l), VRL, K_r, r_r)
     call darcy2d_uptake_lin ( soil, psi_wilt, VRL, K_r, r_r, uptake_oneway, &
          uptake_from_sat, u, du)
     soil_water_supply = max(0.0,sum(u))
     soil_uptake_T = soil%uptake_T
  end select
end subroutine soil_data_beta


! ============================================================================
! update soil properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate soil-heat conduction equation upward from bottom of soil
! to surface, delivering linearization of surface ground heat flux.
subroutine soil_step_1 ( soil, vegn, diag, &
                         soil_T, soil_uptake_T, soil_beta, soil_water_supply, &
                         soil_E_min, soil_E_max, &
                         soil_rh, soil_rh_psi, soil_liq, soil_ice, soil_subl, soil_tf, &
                         soil_G0, soil_DGDT )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: &
       soil_T, &    ! temperature of the upper layer of the soil, degK
       soil_uptake_T, & ! estimate of the temperature of the water taken up by transpiration
       soil_beta, &
       soil_water_supply, & ! supply of water to vegetation per unit total active root biomass, kg/m2 
       soil_E_min, &
       soil_E_max, &
       soil_rh,   & ! soil surface relative humidity
       soil_rh_psi,& ! derivative of soil_rh w.r.t. soil surface matric head
       soil_liq,  & ! amount of liquid water available for implicit freeze (=0)
       soil_ice,  & ! amount of ice available for implicit melt (=0)
       soil_subl, & ! part of sublimation in water vapor flux, dimensionless [0,1]
       soil_tf,   & ! soil freezing temperature, degK
       soil_G0, soil_DGDT ! linearization of ground heat flux
  ! ---- local vars
  real :: bbb, denom, dt_e
  real, dimension(num_l) :: aaa, ccc, thermal_cond, heat_capacity, vlc, vsc
  integer :: l

  if(is_watch_point()) then
     write(*,*) 'soil%tag', soil%tag
     write(*,*) 'soil%pars%k_sat_ref', soil%pars%k_sat_ref 
     write(*,*) 'soil%pars%psi_sat_ref', soil%pars%psi_sat_ref
     write(*,*) 'soil%pars%chb', soil%pars%chb
     write(*,*) 'soil%pars%w_sa', soil%pars%vwc_sat
  endif
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  soil_T = soil%prog(1)%T
  call soil_data_beta ( soil, vegn, soil_beta, soil_water_supply, soil_uptake_T, &
                        soil_rh, soil_rh_psi )

  do l = 1, num_l
    vlc(l) = max(0.0, soil%prog(l)%wl / (dens_h2o * dz(l)))
    vsc(l) = max(0.0, soil%prog(l)%ws / (dens_h2o * dz(l)))
    enddo
  call soil_data_thermodynamics ( soil, vlc, vsc,  &  
                                  soil_E_max, thermal_cond )
  if (.not.use_E_max) soil_E_max =  HUGE(soil_E_max)
  soil_E_min = Eg_min

  do l = 1, num_l
     heat_capacity(l) = soil%heat_capacity_dry(l) *dz(l) &
          + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws
  enddo

  soil_liq  = max(soil%prog(1)%wl, 0.)
  soil_ice  = max(soil%prog(1)%ws, 0.)
  if (soil_liq + soil_ice > 0) then
     soil_subl = soil_ice / (soil_liq + soil_ice)
  else
     soil_subl = 0
  endif
  soil_liq = 0
  soil_ice = 0

  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz(l+1)/thermal_cond(l+1) &
                     + dz(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(soil%prog(num_l)%T - soil%prog(num_l-1)%T)
     soil%e(num_l-1) = -aaa(num_l)/denom
     soil%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*soil%e(l)
        dt_e = - ( ccc(l)*(soil%prog(l+1)%T - soil%prog(l)%T  ) &
                  -aaa(l)*(soil%prog(l)%T   - soil%prog(l-1)%T) )
        soil%e(l-1) = -aaa(l)/denom
        soil%f(l-1) = (dt_e - ccc(l)*soil%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     soil_G0   = ccc(1)*(soil%prog(2)%T- soil%prog(1)%T + soil%f(1)) / denom
     soil_DGDT = (1 - ccc(1)*(1-soil%e(1))) / denom   
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     soil_G0    = 0.
     soil_DGDT  = 1. / denom
  end if
  
  ! set soil freezing temperature
  soil_tf = soil%pars%tfreeze

  if(is_watch_point()) then
     write(*,*) '#### soil_step_1 checkpoint 1 ####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'T       ', soil_T
     write(*,*) 'uptake_T', soil_uptake_T
     write(*,*) 'beta    ', soil_beta
     write(*,*) 'E_max   ', soil_E_max
     write(*,*) 'rh      ', soil_rh
     write(*,*) 'liq     ', soil_liq
     write(*,*) 'ice     ', soil_ice
     write(*,*) 'subl    ', soil_subl
     write(*,*) 'G0      ', soil_G0
     write(*,*) 'DGDT    ', soil_DGDT
     __DEBUG1__(soil_water_supply)
  endif

  call send_tile_data(id_thermal_cond, thermal_cond, diag)

end subroutine soil_step_1


! ============================================================================
! apply boundary flows to soil water and move soil water vertically.
  subroutine soil_step_2 ( soil, vegn, diag, soil_subl, snow_lprec, snow_hlprec,  &
                           vegn_uptk, &
                           subs_DT, subs_M_imp, subs_evap, &
                           use_tfreeze_in_grnd_latent, &
                           soil_levap, soil_fevap, soil_melt, &
                           soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
     soil_subl     !
  real, intent(in) :: &
     snow_lprec, &
     snow_hlprec, &
     vegn_uptk, &
     subs_DT,       &!
     subs_M_imp,       &! rate of phase change of non-evaporated soil water
     subs_evap
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
     soil_levap, soil_fevap, soil_melt, &
     soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) :: del_t, eee, fff, &
             psi, DThDP, hyd_cond, DKDP, K, DKDPm, DKDPp, grad, &
             vlc, vsc, dW_l, u_minus, u_plus, DPsi, soil_w_fc, soil_vwc_sat
  real, dimension(num_l+1) :: flow, infilt
  real, dimension(num_l  ) :: div, dq, div_active
  real      :: &
     lprec_eff, hlprec_eff, tflow, hcap,cap_flow, dheat, &
     melt_per_deg, melt, adj, &
     liq_to_extract, ice_to_extract, heat_of_extract, &
     liq_to_extract_here, ice_to_extract_here, &
     lrunf_sn,lrunf_ie,lrunf_bf,lrunf_nu, hlrunf_sn,hlrunf_ie,hlrunf_bf,hlrunf_nu, &
     Qout, DQoutDP, tau_gw, gw_length, gw_relief, gw_zeta_bar, gw_e_depth, K_sat, &
     c0, c1, c2, x, aaa, bbb, ccc, ddd, xxx, sat_frac, z_sat, &
     gw_flux, storage_frac, depth_to_saturation, &
     Dpsi_min, Dpsi_max, psi_for_rh, &
     liq_frac, excess_wat, excess_liq, excess_ice, h1, h2, summax, &
     space_avail, liq_placed, ice_placed, excess_t, dW_l_internal, w_to_move_up
  logical :: stiff, flag
  real, dimension(num_l-1) :: del_z
  integer :: n_iter, l, ipt, jpt, kpt, fpt, l_internal, l_max_active_layer
  real :: jj,dpsi_alt
  real :: &
       VRL(num_l), & ! volumetric root length, m/m3
       K_r, & ! root membrame permeability, kg/(m3 s)
       r_r, & ! root radius, m
       uptake(num_l),   & ! uptake by roots per layer, kg/(m2 s)
       uptake_tot,      & ! total uptake, kg/(m2 s)
       uptake_pos,      & ! sum of the positive uptake, kg/(m2 s) 
       uptake_T_new, & ! updated average temperature of uptaken water, deg K
       uptake_T_corr,& ! correction for uptake temperature, deg K
       Tu ! temperature of water taken up from (or added to) a layer, deg K
  
  jj = 1.
  flag = .false.
  
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 1 #####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'subs_evap    ', subs_evap
     write(*,*) 'snow_lprec   ', snow_lprec
     write(*,*) 'uptake  ', vegn_uptk
     write(*,*) 'subs_M_imp   ', subs_M_imp
     write(*,*) 'theta_s ', soil%pars%vwc_sat
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') 'level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws+soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws,&
             ' gw=', soil%prog(l)%groundwater
     enddo
  endif

  ! ---- record fluxes ---------
  soil_levap  = subs_evap*(1-soil_subl)
  soil_fevap  = subs_evap*   soil_subl
  soil_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  soil%prog(1)%T = soil%prog(1)%T + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = soil%e(l) * del_t(l) + soil%f(l)
      soil%prog(l+1)%T = soil%prog(l+1)%T + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 2 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') 'level=',l, 'T=', soil%prog(l)%T, &
             'del_t=', del_t(l), 'e=', soil%e(l), 'f=', soil%f(l)
     enddo
  endif

  ! ---- extract evap from soil and do implicit melt --------------------
  IF(LM2) THEN
    do l = 1, num_l
      soil%prog(l)%wl = soil%prog(l)%wl &
                      - soil%uptake_frac(l)*soil_levap*delta_time
    enddo
  ELSE
    soil%prog(1)%wl = soil%prog(1)%wl - soil_levap*delta_time
    soil%prog(1)%ws = soil%prog(1)%ws - soil_fevap*delta_time
  ENDIF
  hcap = soil%heat_capacity_dry(1)*dz(1) &
                       + clw*soil%prog(1)%wl + csw*soil%prog(1)%ws
  ! T adjustment for nonlinear terms (del_T)*(del_W)
  dheat = delta_time*(clw*soil_levap+csw*soil_fevap)*del_T(1)
  ! take out extra heat not claimed in advance for evaporation
  if (use_tfreeze_in_grnd_latent) dheat = dheat &
          - delta_time*((cpw-clw)*soil_levap+(cpw-csw)*soil_fevap) &
                                 *(soil%prog(1)%T-del_T(1)-tfreeze)
  soil%prog(1)%T  = soil%prog(1)%T  + dheat/hcap
  soil%prog(1)%wl = soil%prog(1)%wl + subs_M_imp
  soil%prog(1)%ws = soil%prog(1)%ws - subs_M_imp
  soil%prog(1)%T  = tfreeze + (hcap*(soil%prog(1)%T-tfreeze) ) &
                              / ( hcap + (clw-csw)*subs_M_imp )

  ! calculate actual vertical distribution of uptake
  select case(uptake_option)
  case ( UPTAKE_LINEAR )
     uptake_T_corr = 0
     n_iter = 0
     uptake = soil%uptake_frac*vegn_uptk
  case ( UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN )     
     ! for Darcy-flow uptake, find the root water potential to satify actual
     ! transpiration by the vegetation
     call vegn_root_properties (vegn, dz(1:num_l), VRL, K_r, r_r)
     
     if ( uptake_option==UPTAKE_DARCY2D ) then
        call darcy2d_uptake_solver ( soil, vegn_uptk, VRL, K_r, r_r, &
             uptake_oneway, uptake_from_sat, uptake, n_iter)
     else
        call darcy2d_uptake_solver_lin ( soil, vegn_uptk, VRL, K_r, r_r, &
             uptake_oneway, uptake_from_sat, uptake, n_iter )
     endif

     uptake_pos = sum(uptake(:),mask=uptake(:)>0)
     if (uptake_pos > 0) then
        ! calculate actual temperature of uptake
        uptake_T_new  = sum(uptake*soil%prog%T,mask=uptake>0)/uptake_pos
        ! and temperature correction
        uptake_T_corr = soil%uptake_T - uptake_T_new
        if(is_watch_point()) then
           __DEBUG3__(soil%uptake_T, uptake_T_new, uptake_T_corr)
        endif
        ! save new uptake for the next time step to serve as an estimate of uptake 
        ! temperature
        soil%uptake_T    = uptake_T_new
     else
        uptake_T_corr = 0.0
        ! and don't change the soil%uptake_T
     endif
  case default
     call error_mesg('soil_step_2', 'invalid soil uptake option', FATAL)
  end select

  if (is_watch_point())then
     write(*,*) ' ##### soil_step_2 checkpoint 2.1 #####'
     __DEBUG2__(vegn_uptk,sum(uptake))
     do l = 1,num_l
        write(*,'(a,i2.2,100(2x,a,g))')'level=',l, &
             'uptake=',uptake(l),'dwl=',-uptake(l)*delta_time,&
             'wl=',soil%prog(l)%wl,'new wl=',soil%prog(l)%wl - uptake(l)*delta_time
     enddo
  endif

  call send_tile_data(id_uptk_n_iter, real(n_iter), diag)
  call send_tile_data(id_uptk, uptake, diag)
  call send_tile_data(id_uptk_old, uptake, diag)

  ! update temperature and water content of soil due to root uptake processes 
  do l = 1, num_l
     ! calculate the temperature of water that is taken from the layer (or added 
     ! to the layer), including energy balance correction 
     if (uptake(l) > 0) then
        Tu = soil%prog(l)%T + uptake_T_corr
     else
        Tu = soil%uptake_T + uptake_T_corr
     endif
     ! heat capacity of the layer
     hcap = soil%heat_capacity_dry(l)*dz(l) &
          + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws

     soil%prog(l)%T = soil%prog(l)%T - &
          uptake(l)*delta_time*clw*( Tu-soil%prog(l)%T ) / &
          ( hcap - uptake(l)*delta_time*clw )
     soil%prog(l)%wl = soil%prog(l)%wl - uptake(l)*delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws+soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws
     enddo
  endif
   lrunf_ie=0; lrunf_sn=0; lrunf_bf=0; lrunf_nu=0
  hlrunf_ie=0;hlrunf_sn=0;hlrunf_bf=0;hlrunf_nu=0
  ! ---- push down any excess surface water, with heat ---------------------
  call soil_data_vwc_sat(soil, soil_vwc_sat)
  liq_frac=0;excess_wat=0;excess_liq=0;excess_ice=0;h1=0;h2=0;liq_frac=0
  l = 1
  summax = max(0.,soil%prog(l)%wl)+max(0.,soil%prog(l)%ws)
  if (summax > 0) then
     liq_frac = max(0.,soil%prog(l)%wl) / summax
  else
     liq_frac = 1
  endif
  excess_wat = max(0., soil%prog(l)%wl + soil%prog(l)%ws &
       - dens_h2o*dz(l)*soil_vwc_sat(l) )
  excess_liq = excess_wat*liq_frac
  excess_ice = excess_wat-excess_liq
  excess_t   = soil%prog(l)%T
  soil%prog(l)%wl = soil%prog(l)%wl - excess_liq
  soil%prog(l)%ws = soil%prog(l)%ws - excess_ice
  call send_tile_data(id_excess, excess_wat/delta_time, diag)

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.001 #####'
     write(*,*) ' level=', l,&
          ' summax =', summax,&
          ' liq_frac =', liq_frac,&
          ' soil_vwc_sat =', soil_vwc_sat(l),&
          ' excess_liq =', excess_liq,&
          ' excess_ice =', excess_ice, &
          ' dens_h2o=', dens_h2o, &
          ' dz(l)=',dz(l),&
          'friday am'
  endif

  do l = 2, num_l
     if (excess_liq+excess_ice>0) then
        space_avail = dens_h2o*dz(l)*soil_vwc_sat(l) &
             - (soil%prog(l)%wl + soil%prog(l)%ws)
        liq_placed = max(min(space_avail, excess_liq), 0.)
        ice_placed = max(min(space_avail-liq_placed, excess_ice), 0.)
        h1 = (soil%heat_capacity_dry(l)*dz(l) &
             + csw*soil%prog(l)%ws + clw*soil%prog(l)%wl)
        h2 = liq_placed*clw+ice_placed*csw
        soil%prog(l)%T = (h1 * soil%prog(l)%T &
             + h2 * excess_T )  / (h1+h2)
        soil%prog(l)%wl = soil%prog(l)%wl + liq_placed
        soil%prog(l)%ws = soil%prog(l)%ws + ice_placed
        excess_liq = excess_liq - liq_placed
        excess_ice = excess_ice - ice_placed
     endif
  enddo

! to avoid adding frozen runoff to soil interface, melt all remaining
! excess ice, even if it results in supercooled liquid runoff
   lrunf_nu = (excess_liq+excess_ice) / delta_time
  hlrunf_nu = (  excess_liq*clw*(excess_T-tfreeze)  &
               + excess_ice*csw*(excess_T-tfreeze)  &
               - hlf*excess_ice                   ) / delta_time

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.01 #####'
     write(*,*) ' lrunf_nu',lrunf_nu
     write(*,*) 'hlrunf_nu',hlrunf_nu
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,100(x,a,g))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws
     enddo
  endif

  ! ---- fetch soil hydraulic properties -------------------------------------
  vlc=0;vsc=0
  do l = 1, num_l
    vlc(l) = max(0., soil%prog(l)%wl / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%prog(l)%ws / (dens_h2o*dz(l)))
  enddo
  call soil_data_hydraulics (soil, vlc, vsc, &
                   psi, DThDP, hyd_cond, DKDP, Dpsi_min, Dpsi_max, tau_gw, &
                   psi_for_rh, soil_w_fc )

  IF (lm2) THEN ! ********************************

     if(is_watch_point()) then
        write(*,*) ' ##### soil_step_2 checkpoint 3.1 #####'
        do l = 1, num_l
           write(*,'(x,a,x,i2.2,100(x,a,g))')'level=', l, 'vlc', vlc(l), 'K  ', hyd_cond(l)
        enddo
     endif
  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
    flow=1
    flow(1)  = 0
    do l = 1, num_l
      infilt(l) = soil%uptake_frac(l)*snow_lprec *delta_time
      flow(l+1) = max(0., soil%prog(l)%wl + flow(l) &
            + infilt(l) - soil_w_fc(l)*dz(l)*dens_h2o)
      dW_l(l) = flow(l) - flow(l+1) + infilt(l)
      soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
    enddo
    do l = 1, num_l
      flow(l) = flow(l) + infilt(l)
    enddo
    dW_l=0
    dpsi=0
    lrunf_bf = lrunf_bf + flow(num_l)/delta_time
  ELSE   ! ********************************
    IF (.NOT.USE_GEOHYDROLOGY) THEN
        div = 0.
        IF (CORRECTED_LM2_GW) THEN
            do l = 1, num_l
              if (vlc(l) .ge. soil_vwc_sat(l) .and. vsc(l).le.0.) &
                  div(l) = 0.15*dens_h2o*dz(l)/tau_gw
              enddo
          ELSE
            do l = 1, num_l
              if ((vsc(l)+vlc(l)) .ge. soil_vwc_sat(l)) &
                  div(l) = 0.15*dens_h2o*dz(l)*(vlc(l)/(vsc(l)+vlc(l)))/tau_gw
              enddo
          ENDIF
        z_sat = 0.
        do l=num_l,1,-1
           if(vsc(l)+vlc(l).le.soil_vwc_sat(l)) exit
           z_sat = z_sat + dz(l)
        enddo
        sat_frac = min((z_sat/zhalf(num_l+1))**soil%pars%rsa_exp,1.)
      ELSE
        call soil_data_gw_hydraulics(soil, zfull(num_l), psi(num_l), &
                gw_flux, sat_frac, storage_frac, depth_to_saturation)
        dq = 0.
        z_sat = 0.
        l = num_l
        div_active = 0.
        IF (BASEFLOW_WHERE_FROZEN) THEN
           do l=num_l,1,-1
              if(psi(l).le.0.) exit
              dq(l) = dz(l)*vlc(l)/(vlc(l)+vsc(l))
              z_sat = z_sat + dz(l)
           enddo
          ELSE
             do l=num_l,1,-1
                if(psi(l).le.0.) exit
                if (vsc(l).le.0.) dq(l) = dz(l)
                z_sat = z_sat + dz(l)
             enddo
         !   IF (DRAIN_ACTIVE_LAYER) THEN
                l_max_active_layer = 0
                do l=1,num_l
                   if(vsc(l).gt.0.) exit
                   l_max_active_layer = l
                enddo
                if (l_max_active_layer.lt.num_l .and. l_max_active_layer.gt.0) then
                    do l = 1, l_max_active_layer
                       if(vlc(l).gt.0) &
                           div_active(l) = hyd_cond(l) * soil%pars%hillslope_relief*dz(l) &
                                / (soil%pars%hillslope_length*soil%pars%hillslope_length)
                      enddo
                  endif
         !     ENDIF
          ENDIF
        div = 0.
        if (z_sat.gt.0.) div = (dq / z_sat) * gw_flux
        where (div.eq.0.) div = div_active*active_layer_drainage_acceleration
        call send_tile_data(id_psi_bot, psi(num_l), diag)
        call send_tile_data(id_psi_bot_old, psi(num_l), diag)
        call send_tile_data(id_sat_frac, sat_frac, diag)
        call send_tile_data(id_sat_frac_old, sat_frac, diag)
        call send_tile_data(id_stor_frac, storage_frac, diag)
        call send_tile_data(id_stor_frac_old, storage_frac, diag)
        if (depth_to_saturation .ge. -0.5) call send_tile_data(id_sat_depth, depth_to_saturation, diag)
        if (depth_to_saturation .ge. -0.5) call send_tile_data(id_sat_depth_old, depth_to_saturation, diag)
      ENDIF
    lrunf_bf = lrunf_bf + sum(div)
    if (snow_lprec.ne.0.) then
      lrunf_sn = sat_frac * snow_lprec
      hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
    else
      lrunf_sn = 0.
      hlrunf_sn = 0.
    endif


    if(is_watch_point()) then
       do l = 1, num_l
          write(*,'(a,1x,i2.2,100(2x,g))')'div,vsc,psi,dz',l,div(l),vsc(l),psi(l),dz(l)
       enddo
       write(*,*)'lrunf_bf',lrunf_bf
       write(*,*)'tau_gw',tau_gw
       write(*,*)'dens_h2o',dens_h2o
    endif
  ! ---- soil-water flow ----------------------------------------------------
    flow = 0
    stiff = all(DThDP.eq.0)
    lprec_eff = snow_lprec - lrunf_sn
    hlprec_eff = snow_hlprec - hlrunf_sn
IF(stiff .AND. BYPASS_RICHARDS_WHEN_STIFF) THEN   ! BYPASS_RICHARDS_WHEN_STIFF
flow = 0.
div  = 0.
dW_l = 0.
lrunf_ie = lprec_eff
hlrunf_ie = hlprec_eff
lrunf_bf = 0.
hlrunf_bf =0.
psi=-zfull
dpsi=0.
ELSE                                              ! BYPASS_RICHARDS_WHEN_STIFF
    flow(1) = delta_time*lprec_eff
    do l = 1, num_l-1
      del_z(l) = zfull(l+1)-zfull(l)
      K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
      DKDPm(l) = 0. !0.5*DKDP(l)
      DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = jj*(psi(l+1)-psi(l))/del_z(l) - 1
    enddo

    if(is_watch_point()) then
       write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
       do l = 1, num_l
          write(*,'(x,a,x,i2.2,x,a,100(x,g))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
               DThDP(l),&
               hyd_cond(l),&
               psi(l),&
               DKDP(l)
       enddo
       do l = 1, num_l-1
          write(*,'(a,i2.2,1x,a,100(2x,g))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
               K(l),&
               DKDPm(l),&
               DKDPp(l),&
               grad(l)
       enddo
    endif


    l = num_l
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    aaa =     - ( jj* K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
!      where (stiff)
    bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
    ddd = - K(l-1) *grad(l-1) - div(l)
!        elsewhere
!          Qout = hyd_cond(l) ! gravity drainage
!          DQoutDP = DKDP(l)  ! gravity drainage
!          Qout = 0.                ! no drainage
!          DQoutDP = 0.             ! no drainage
!          where (psi(l).gt.0.) ! linear baseflow from gw
!              Qout = 0.15*psi(l)/tau_gw
!              DQoutDP = 0.15/tau_gw
!            elsewhere
!              Qout = 0.
!              DQoutDP = 0.
!            endwhere
!          bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
!                      -DQoutDP )
!          ddd = -Qout - K(l-1) *grad(l-1)
!        endwhere
    eee(l-1) = -aaa/bbb
    fff(l-1) =  ddd/bbb
  
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g))') 'l,a,b, ,d', l,aaa, bbb,ddd
    endif

    do l = num_l-1, 2, -1
      xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
      aaa = - ( jj*K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
      bbb = xxx-( -jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                  -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
      ccc =   - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
      ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                            - div(l)
      eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
      fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      if(is_watch_point()) then
         write(*,'(a,i2.2,100(2x,g))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
      endif
    enddo
  
    l = 1
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    bbb = xxx - ( -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =     - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =          flow(1)/delta_time +    K(l)     *grad(l) &
                            - div(l)

IF (bbb+ccc*eee(l) .NE. 0.) THEN
    dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
    if (.not.stiff .and. dPsi(l).gt.Dpsi_min .and. dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
      else
      	if (stiff) then
            dPsi(l) = - psi(l)
          else
            dPsi(l) = min (dPsi(l), Dpsi_max)
            dPsi(l) = max (dPsi(l), Dpsi_min)
          endif
        flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                      - K(l)*grad(l))*delta_time
        lrunf_ie = lprec_eff - flow(l)/delta_time
        if (.not.allow_negative_rie.and.lrunf_ie.lt.-lrunf_ie_tol) then
            flag = .true.
            call get_current_point(ipt,jpt,kpt,fpt)
            dpsi_alt = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
            write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'rie= ',lrunf_ie,' reset to 0'
            write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi=',dPsi(l), ' reset to ',dpsi_alt
            write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi_min/max=',dPsi_min,dPsi_max
            dPsi(l) = dpsi_alt
            lrunf_ie = 0.
            flow(l) = lprec_eff*delta_time
          endif
      endif
  ELSE
      	if (stiff) then
            dPsi(l) = - psi(l)
          else
            dPsi(l) = Dpsi_max
          endif
        flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                      - K(l)*grad(l))*delta_time
        lrunf_ie = lprec_eff - flow(l)/delta_time
        if (.not.allow_negative_rie.and.lrunf_ie.lt.-lrunf_ie_tol) then
            flag = .true.
            call get_current_point(ipt,jpt,kpt,fpt)
       ! next change will not change answers in previous runs, since old version would crash
       ! the only time this point was reached was when DThDP was zero everywhere.
       !     dpsi_alt = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
            dpsi_alt = 0.
            write(*,*) 'note 2: at point ',ipt,jpt,kpt,fpt,'rie=',lrunf_ie,' reset to 0'
            write(*,*) 'note 2: at point ',ipt,jpt,kpt,fpt,'dPsi=',dPsi(l),' reset to',dpsi_alt
            write(*,*) 'note 2: at point ',ipt,jpt,kpt,fpt,'dPsi_min/max=',dPsi_min,dPsi_max
            dPsi(l) = dpsi_alt
            lrunf_ie = 0.
            flow(l) = lprec_eff*delta_time
          endif
  ENDIF
      
    if(is_watch_point().or.(flag.and.write_when_flagged)) then
       write(*,'(a,i2.2,100(2x,g))') 'l,  b,c,d', l, bbb,ccc,ddd
       write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
       write(*,*) 'ie,sn,bf:', lrunf_ie,lrunf_sn,lrunf_bf
       do l = 1, num_l-1
          write(*,'(a,i2.2,100(2x,g))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
       enddo
       write(*,*) 'DThDP(1)', DThDP(1)
       write(*,*) 'K(1)', K(1)
       write(*,*) 'grad(1)', grad(1)
       write(*,*) 'ddd(1)', ddd
       write(*,*) 'ccc(1)', ccc
       write(*,*) 'bbb(1)', bbb
       write(*,*) 'dPsi(1)', dPsi(1)
       write(*,*) 'Psi(1)', Psi(1)
       write(*,*) 'div(1)', div(1)
    endif

    do l = 2, num_l
      dPsi(l) = eee(l-1)*dPsi(l-1) + fff(l-1)
    enddo
  
    l_internal = 1
    dW_l_internal = -1.e20
    do l = 1, num_l-1
      flow(l+1) = delta_time*( &
           -K(l)*(grad(l)&
           +jj*(DPsi(l+1)-DPsi(l))/ del_z(l)) &
           -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                           DKDPm(l)*Dpsi(l) )  )
      dW_l(l) = flow(l) - flow(l+1) - div(l)*delta_time
      if (flag .and. l.gt.1. .and. dW_l(l).gt.dW_l_internal) then
          l_internal = l
          dW_l_internal = dW_l(l)
        endif
      enddo
    flow(num_l+1) = 0.
    dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                            - div(num_l)*delta_time
    if (flag .and. dW_l(num_l).gt.dW_l_internal) then
        l_internal = num_l
        dW_l_internal = dW_l(num_l)
      endif

    if(is_watch_point().or.(flag.and.write_when_flagged)) then
       write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
       do l = 1, num_l
          write(*,'(i2.2,100(2x,a,g))') l,&
               ' dW_l=', dW_l(l),&
               ' flow=', flow(l),&
               ' div=', div(l)
       enddo
    endif

    if (flag) then
        w_to_move_up = min(dW_l_internal, -(soil%prog(1)%wl+dW_l(1)))
        w_to_move_up = max(w_to_move_up, 0.)
        write(*,*) 'l_internal=',l_internal
        write(*,*) 'dW_l(l_internal)=',dW_l(l_internal)
        write(*,*) 'soil%prog(1)%wl+dW_l(1)',soil%prog(1)%wl+dW_l(1)
        write(*,*) 'w_to_move_up=',w_to_move_up
        if (l_internal.gt.1) then
            dW_l(1) = dW_l(1) + w_to_move_up
            dW_l(l_internal) = dW_l(l_internal) - w_to_move_up
            do l = 2, l_internal
              flow(l) = flow(l) - w_to_move_up
              enddo
          endif
      endif

    if(is_watch_point().or.(flag.and.write_when_flagged)) then
       write(*,*) ' ##### soil_step_2 checkpoint 3.22 #####'
       do l = 1, num_l
          write(*,'(i2.2,100(2x,a,g))') l,&
               ' dW_l=', dW_l(l),&
               ' flow=', flow(l),&
               ' div=', div(l)
       enddo
    endif

! In rare situations where lrunf_ie is large and negative, clip any liquid supersaturation
! layer by layer and recompute lrunf_ie (this is not good, since it ignores 'comp'):
  IF (lrunf_ie < lrunf_ie_min) THEN
       call get_current_point(ipt,jpt,kpt,fpt)
       write(*,*) 'note: at point ',ipt,jpt,kpt,fpt,' clip triggered by lrunf_ie=',lrunf_ie
       do l = num_l, 1, -1
          adj = max(dW_l(l)+soil%prog(l)%ws+soil%prog(l)%wl &
               - soil_vwc_sat(l)*dz(l)*dens_h2o, 0. )

          if(is_watch_point()) then
             write(*,*) '3.22 l=', l,&
                  ' soil_prog%wl=',soil%prog(l)%wl,  &
                  ' soil_prog%ws=',soil%prog(l)%ws , &
                  ' soil_vwc_sat=', soil_vwc_sat(l), &
                  ' dz=', dz(l), &
                  ' adj=', adj
          endif

          adj = min(adj, max(0.,soil%prog(l)%wl))

          if(is_watch_point()) then
             write(*,*) '3.23 l=', l, ' adj=', adj
          endif

          dW_l(l) = dW_l(l) - adj
          flow(l) = flow(l+1) + dW_l(l) + div(l)*delta_time
       enddo
       lrunf_ie = lprec_eff - flow(1)/delta_time

  ELSE IF (UNCONDITIONAL_SWEEP) THEN
  ! Sweep and fill upward any liquid supersaturation
  ! USE OF THIS EXPERIMENTAL CODE IS NOT RECOMMENDED. 
  ! CONFLICTS WITH COMP.NE.0 !!!
       excess_liq = 0.
       do l = num_l, 1, -1
         adj = dW_l(l)+soil%prog(l)%ws+soil%prog(l)%wl &
                        - soil_vwc_sat(l)*dz(l)*dens_h2o
         if (adj.gt.0.) then  ! collect excess liquid
             adj = min(adj, max(0.,soil%prog(l)%wl))
             dW_l(l) = dW_l(l) - adj
             excess_liq = excess_liq + adj
           else if (adj.lt.0.) then  ! deposit collected liquid
             if (excess_liq.gt.0.) then
                 adj = min(-adj, excess_liq)
                 dW_l(l) = dW_l(l) + adj
                 excess_liq = excess_liq - adj
               endif
           endif
         flow(l) = flow(l+1) + dW_l(l) + div(l)*delta_time
         enddo
       lrunf_ie = lprec_eff - flow(1)/delta_time
  ENDIF
       
       do l = 1, num_l
         soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
         enddo

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g))') l, &
             'Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)), &
             'wl=', soil%prog(l)%wl, &
             'ws=', soil%prog(l)%ws, &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

ENDIF                                              ! BYPASS_RICHARDS_WHEN_STIFF
  ENDIF ! ************************************

  if  (snow_lprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4 ***** '
     write(*,*) ' tfreeze', tfreeze
     write(*,*) '  tflow ', tflow
     write(*,*) ' snow_hlprec', snow_hlprec
  endif


! Upstream weighting of advection. Preserving u_plus here for now.
  u_minus = 1.
  where (flow.lt.0.) u_minus = 0.
  do l = 1, num_l-1
    u_plus(l) = 1. - u_minus(l+1)
    enddo
  hcap = (soil%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*soil%prog(num_l)%ws)/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + soil%prog(num_l)%wl - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(soil%prog(num_l)%T-soil%prog(num_l-1)%T) / bbb

  do l = num_l-1, 2, -1
    hcap = (soil%heat_capacity_dry(l)*dz(l) &
                              + csw*soil%prog(l)%ws)/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + soil%prog(l)%wl - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(soil%prog(l)%T-soil%prog(l-1)%T)    &
                       + ccc*(soil%prog(l)%T-soil%prog(l+1)%T)    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo
    
  hcap = (soil%heat_capacity_dry(1)*dz(1) + csw*soil%prog(1)%ws)/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + soil%prog(1)%wl - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(soil%prog(1)%T-tflow          ) &
                     + ccc*(soil%prog(1)%T-soil%prog(2)%T) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  soil%prog(1)%T = soil%prog(1)%T + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', soil%prog(1)%T
  endif

  do l = 1, num_l-1
    del_t(l+1) = eee(l)*del_t(l) + fff(l)
    soil%prog(l+1)%T = soil%prog(l+1)%T + del_t(l+1)
  enddo

  tflow = soil%prog(num_l)%T

!  do l = 1, num_l
!    where (mask)
!        hcap = soil%heat_capacity_dry(l)*dz(l) &
!                 + clw*(soil%prog(l)%wl-dW_l(l)) + csw*soil%prog(l)%ws
!        cap_flow = clw*flow(l)
!        soil%prog(l)%T = (hcap*soil%prog(l)%T + cap_flow*tflow) &
!                         /(hcap                 + cap_flow      )
!        tflow  = soil%prog(l)%T
!      endwhere
!    enddo

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.5 ***** '
     write(*,*) 'hcap', hcap
!    write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,*) 'level=', l, ' T', soil%prog(l)%T
     enddo
  endif

  ! ---- groundwater ---------------------------------------------------------
  ! THIS T AVERAGING IS WRONG, BECAUSE IT NEGLECTS THE MEDIUM  ***
  ! ALSO, FREEZE-THAW IS NEEDED!
  ! PROBABLY THIS SECTION WILL BE DELETED ANYWAY, WITH GW TREATED ABOVE.
  IF (lm2) THEN
    do l = 1, 1      !TEMPORARY LAYER THING !!!!***
      if (soil%prog(l)%groundwater + flow(num_l+1) .ne. 0.) then ! TEMP FIX
          soil%prog(l)%groundwater_T =    &
           + (soil%prog(l)%groundwater*soil%prog(l)%groundwater_T &
              + flow(num_l+1)*tflow) &
            /(soil%prog(l)%groundwater + flow(num_l+1))
      endif
      c0 = delta_time/tau_gw
      c1 = exp(-c0)
      c2 = (1-c1)/c0
      x  = (1-c1)*soil%prog(l)%groundwater/delta_time &
                          + (1-c2)*flow(num_l+1)/delta_time
      soil%prog(l)%groundwater = c1 * soil%prog(l)%groundwater &
                                + c2 * flow(num_l+1)
      soil_lrunf  = x
      soil_hlrunf = x*clw*(soil%prog(l)%groundwater_T-tfreeze)
    enddo
  ELSE
    if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
      hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
    else if (flow(1).lt.0. ) then
      hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(soil%prog(1)%T-tfreeze)
    else
      hlrunf_ie = 0.
    endif
    hlrunf_bf = hlrunf_bf + clw*sum(div*(soil%prog%T-tfreeze))


    soil_lrunf  =  lrunf_sn +  lrunf_ie +  lrunf_bf +  lrunf_nu
    soil_hlrunf = hlrunf_sn + hlrunf_ie + hlrunf_bf + hlrunf_nu
  ENDIF

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = soil%heat_capacity_dry(l)*dz(l) &
             + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws
    melt_per_deg = hcap/hlf
    if       (soil%prog(l)%ws>0 .and. soil%prog(l)%T>soil%pars%tfreeze) then
      melt =  min(soil%prog(l)%ws, (soil%prog(l)%T-soil%pars%tfreeze)*melt_per_deg)
    else if (soil%prog(l)%wl>0 .and. soil%prog(l)%T<soil%pars%tfreeze) then
      melt = -min(soil%prog(l)%wl, (soil%pars%tfreeze-soil%prog(l)%T)*melt_per_deg)
    else
      melt = 0
    endif
    soil%prog(l)%wl = soil%prog(l)%wl + melt
    soil%prog(l)%ws = soil%prog(l)%ws - melt
    soil%prog(l)%T = tfreeze &
       + (hcap*(soil%prog(l)%T-tfreeze) - hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
    soil_melt = soil_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 5 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws,&
             ' gw=', soil%prog(l)%groundwater
     enddo
  endif


  soil_Ttop = soil%prog(1)%T
  soil_Ctop = soil%heat_capacity_dry(1)*dz(1) &
    + clw*soil%prog(1)%wl + csw*soil%prog(1)%ws


! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!  

  ! ---- increment time and do diagnostics -----------------------------------
  time = increment_time(time, int(delta_time), 0)
  
  ! ---- diagnostic section
   call send_tile_data(id_temp, soil%prog%T, diag)
   if (id_lwc > 0) call send_tile_data(id_lwc,  soil%prog%wl/dz(1:num_l), diag)
   if (id_swc > 0) call send_tile_data(id_swc,  soil%prog%ws/dz(1:num_l), diag)
   if (id_psi > 0) call send_tile_data(id_psi,  psi+dPsi, diag)
   call send_tile_data(id_ie,   lrunf_ie, diag)
   call send_tile_data(id_sn,   lrunf_sn, diag)
   call send_tile_data(id_bf,   lrunf_bf, diag)
   call send_tile_data(id_nu,   lrunf_nu, diag)
   call send_tile_data(id_hie,  hlrunf_ie, diag)
   call send_tile_data(id_hsn,  hlrunf_sn, diag)
   call send_tile_data(id_hbf,  hlrunf_bf, diag)
   call send_tile_data(id_hnu,  hlrunf_nu, diag)
   if (id_evap > 0) call send_tile_data(id_evap,  soil_levap+soil_fevap, diag)

   call send_tile_data(id_heat_cap, soil%heat_capacity_dry, diag)

end subroutine soil_step_2


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function soil_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   soil_tile_exists = associated(tile%soil)
end function soil_tile_exists


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_SOIL_ACCESSOR_0D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;end subroutine
#define DEFINE_SOIL_ACCESSOR_1D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_0D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_1D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;end subroutine

DEFINE_SOIL_ACCESSOR_1D(real,w_fc)
DEFINE_SOIL_ACCESSOR_0D(real,uptake_T)
DEFINE_SOIL_ACCESSOR_0D(integer,tag)

DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau_groundwater)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_length)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_relief)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_zeta_bar)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,soil_e_depth)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_wilt)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_fc)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_ref)

DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_sat)

DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,T)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,wl)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,ws)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,groundwater)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,groundwater_T)

end module soil_mod



