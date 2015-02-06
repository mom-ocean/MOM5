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
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR
use soil_tile_mod, only : GW_LM2, GW_LINEAR, GW_HILL_AR5, GW_HILL, GW_TILED, &
     soil_tile_type, soil_pars_type, soil_prog_type, read_soil_data_namelist, &
     soil_data_radiation, soil_data_diffusion, soil_data_thermodynamics, &
     soil_data_hydraulic_properties, soil_data_psi_for_rh, &
     soil_data_gw_hydraulics, soil_data_gw_hydraulics_ar5, &
     soil_data_vwc_for_init_only, &
     soil_data_init_derive_subsurf_pars, &
     soil_data_init_derive_subsurf_pars_ar5,&
     max_lev, psi_wilt, cpw, clw, csw, g_iso, g_vol, g_geo, g_RT, aspect,&
     num_storage_pts, gw_zeta_s, gw_flux_table, gw_area_table, &
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth, &
     slope_exp, &
     num_zeta_pts, num_tau_pts, &
     log_rho_table, log_zeta_s, log_tau, gw_scale_perm, &
     z_ref, k_macro_constant, use_tau_fix

use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, get_elmt_indices, &
     operator(/=)
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr, &
     add_tiled_diag_field_alias, add_tiled_static_field_alias
use land_data_mod, only : land_state_type, lnd
use land_io_mod, only : read_field
use land_tile_io_mod, only : create_tile_out_file, write_tile_data_r0d_fptr,& 
     write_tile_data_r1d_fptr, read_tile_data_r0d_fptr, read_tile_data_r1d_fptr,&
     print_netcdf_error, get_input_restart_name, sync_nc_files
use nf_utils_mod, only : nfu_def_dim, nfu_put_att, nfu_inq_var
use vegn_tile_mod, only : vegn_tile_type, vegn_uptake_profile, vegn_root_properties
use land_debug_mod, only : is_watch_point, get_current_point, check_var_range
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
public :: soil_step_3
! =====end of public interfaces ==============================================



! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'soil',&
    version     = '$Id: soil.F90,v 20.0 2013/12/13 23:30:48 fms Exp $',&
    tagname     = '$Name: tikal $'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_E_min            = .false.     ! prohibit condensation
logical :: use_E_max            = .true.      ! theoretical effiltration capacity flag
real    :: init_temp            = 288.        ! cold-start soil T
real    :: init_w               = 150.        ! cold-start w(l)/dz(l)
real    :: init_wtdep           = -1.         ! positive value activates hydrostatic IC,
                                              ! overriding init_w
real    :: init_groundwater     =   0.        ! cold-start gw storage
real    :: lrunf_ie_min         = -1.0e-4     ! trigger for clip and runoff
real    :: lrunf_ie_tol         =  1.e-12
character(len=16) :: albedo_to_use = ''       ! or 'albedo-map' or 'brdf-maps'
character(len=24) :: uptake_to_use = 'linear' ! or 'darcy2d', or 'darcy2d-linearized'
logical :: uptake_oneway        = .false.     ! if true, roots can't loose water to soil
logical :: uptake_from_sat      = .true.      ! if false, the uptake from saturated soil is prohibited
logical :: allow_negative_rie   = .false.
logical :: baseflow_where_frozen = .false.
logical :: write_when_flagged   = .false.
logical :: bypass_richards_when_stiff = .true.
logical :: corrected_lm2_gw     = .true.
logical :: use_stiff_bug        = .false.
logical :: fix_z_bot            = .false.
logical :: update_psi           = .false.
logical :: consistent_d_trans   = .false.
logical :: fix_interp           = .false.
logical :: use_new_dq           = .false.
logical :: use_fringe           = .false.
logical :: push_down_sfc_excess = .true.
logical :: lrunf_from_div       = .true.
real    :: active_layer_drainage_acceleration = 0.
real    :: hlf_factor           = 1.
real    :: gw_flux_max          = 1.e10
real    :: log_rho_max          = 100.
real    :: aquifer_heat_cap     = 0.         ! in equivalent liquid water amount, kg/m2
logical :: write_soil_carbon_restart = .FALSE. ! indicates whether to write
                        ! information for soil carbon acceleration

namelist /soil_nml/ lm2, use_E_min, use_E_max,           &
                    init_temp,      &
                    init_w,   init_wtdep,    &
                    init_groundwater, lrunf_ie_min, lrunf_ie_tol, &
                    cpw, clw, csw, &
                    albedo_to_use, &
                    uptake_to_use, uptake_oneway, uptake_from_sat, &
                    allow_negative_rie, &
                    baseflow_where_frozen, &
                    write_when_flagged, &
                    bypass_richards_when_stiff, corrected_lm2_gw, &
                    use_stiff_bug, fix_z_bot, update_psi, &
                    consistent_d_trans, fix_interp, use_new_dq, use_fringe, &
                    push_down_sfc_excess, lrunf_from_div, &
                    active_layer_drainage_acceleration, hlf_factor, &
                    gw_flux_max, log_rho_max, aquifer_heat_cap, &
                    write_soil_carbon_restart
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf = .false.
type(time_type) :: time
real            :: delta_time
logical         :: use_single_geo
integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)
real            :: Eg_min

integer         :: uptake_option = -1 
integer         :: gw_option = -1 

! ---- diagnostic field IDs
integer :: id_fast_soil_C, id_slow_soil_C, id_fsc, id_ssc, &
    id_lwc, id_swc, id_psi, id_temp, &
    id_ie, id_sn, id_bf, id_if, id_al, id_nu, id_sc, &
    id_hie, id_hsn, id_hbf, id_hif, id_hal, id_hnu, id_hsc, &
    id_heat_cap, id_thermal_cond, id_type, id_tau_gw, id_slope_l, &
    id_slope_Z, id_zeta_bar, id_e_depth, id_vwc_sat, id_vwc_fc, &
    id_vwc_wilt, id_K_sat, id_K_gw, id_w_fc, &
    id_refl_dry_dif, id_refl_dry_dir, id_refl_sat_dif, id_refl_sat_dir, &
    id_f_iso_dry, id_f_vol_dry, id_f_geo_dry, &
    id_f_iso_sat, id_f_vol_sat, id_f_geo_sat, &
    id_evap, id_uptk_n_iter, id_uptk, id_psi_x0, id_uptk_residual, &
    id_excess, id_deficit, id_deficit_2, id_deficit_3, id_zeta, id_tau, &
    id_psi_bot, id_sat_frac, id_stor_frac, id_sat_depth, id_sat_dept2, &
    id_cf_1, id_cf_3, id_wt_1, id_wt_2, id_wt_2a, id_wt_3, id_wt2_3, &
    id_div_bf, id_div_if, id_div_al, &
    id_z_cap, id_active_layer

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

  call read_soil_data_namelist(num_l,dz,use_single_geo,gw_option)

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

  ! ---- convert symbolic names of options into numeric IDs to speed up
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
  integer :: unit, unit1  ! unit numbers for various i/o
  type(land_tile_enum_type)     :: te,ce  ! tail and current tile list elements
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  real, allocatable :: gw_param(:,:), gw_param2(:,:), albedo(:,:,:) ! input data buffers for respective variables
  real, allocatable :: f_iso(:,:,:), f_vol(:,:,:), f_geo(:,:,:), refl_dif(:,:,:)

  integer :: i
  real :: psi(num_l), mwc(num_l)
  character(len=256) :: restart_file_name
  logical :: restart_exists

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  call uptake_init(num_l,dz,zfull)

  ! -------- initialize soil model diagnostic fields
  call soil_diag_init ( id_lon, id_lat, id_band )
  
  ! -------- read spatially distributed fields for groundwater parameters, if requested
  if (.not.use_single_geo) then
     select case (gw_option)
     case (GW_LINEAR,GW_LM2)
        allocate(gw_param(lnd%is:lnd%ie,lnd%js:lnd%je))
        call read_field( 'INPUT/groundwater_residence.nc','tau', lnd%lon, lnd%lat, &
             gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_tau_groundwater_ptr )
        deallocate(gw_param)
     case (GW_HILL, GW_HILL_AR5)
        allocate(gw_param (lnd%is:lnd%ie,lnd%js:lnd%je))
        allocate(gw_param2(lnd%is:lnd%ie,lnd%js:lnd%je))
        call read_field( 'INPUT/geohydrology.nc','hillslope_length',  lnd%lon, lnd%lat, &
          gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, lnd%tile_map, soil_hillslope_length_ptr )
        call read_field( 'INPUT/geohydrology.nc','slope', lnd%lon, lnd%lat, &
          gw_param2, interp='bilinear' )
        gw_param = gw_param*gw_param2
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, lnd%tile_map, soil_hillslope_relief_ptr )
        call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
          lnd%lon, lnd%lat, gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param, lnd%tile_map, soil_hillslope_zeta_bar_ptr )
        call read_field( 'INPUT/geohydrology.nc','soil_e_depth', &
          lnd%lon, lnd%lat, gw_param, interp='bilinear' )
        if (slope_exp.gt.0.01) then
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                  lnd%tile_map, soil_soil_e_depth_ptr )
        else
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, lnd%tile_map, soil_soil_e_depth_ptr )
        endif
        if (gw_option /= GW_HILL_AR5) then
            call read_field( 'INPUT/geohydrology.nc','perm', lnd%lon, lnd%lat, &
                 gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, lnd%tile_map, &
                                            soil_k_sat_gw_ptr )
        endif
        deallocate(gw_param, gw_param2)
        te = tail_elmt (lnd%tile_map)
        ce = first_elmt(lnd%tile_map)
        do while(ce /= te)
            tile=>current_tile(ce)  ! get pointer to current tile
            ce=next_elmt(ce)        ! advance position to the next tile
            if (.not.associated(tile%soil)) cycle
            select case (gw_option)
            case (GW_HILL)
                call soil_data_init_derive_subsurf_pars(tile%soil)
            case (GW_HILL_AR5)
                call soil_data_init_derive_subsurf_pars_ar5(tile%soil)
            end select
        enddo
     end select ! gw_option
  endif ! single geo

  ! -------- set dry soil albedo values, if requested
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
          '" is invalid, use "albedo-map", "brdf-maps", or empty line ("")',&
          FATAL)
  endif
  
  ! -------- initialize soil state --------
  te = tail_elmt (lnd%tile_map)
  ce = first_elmt(lnd%tile_map)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)        ! advance position to the next tile
     if (.not.associated(tile%soil)) cycle
     if (init_wtdep .gt. 0.) then
         psi = zfull - init_wtdep
	 call soil_data_vwc_for_init_only(tile%soil, psi, mwc)
	 mwc = mwc * dens_h2o
       else if (init_w .ge. 0.) then
         mwc = init_w
       else ! negative init_w is to be intrepreted as prescribed saturation
         mwc = -init_w*tile%soil%pars%vwc_sat*dens_h2o
       endif
     if (init_temp.ge.tile%soil%pars%tfreeze) then
         tile%soil%prog%wl = mwc*dz(1:num_l)
         tile%soil%prog%ws = 0
       else
         tile%soil%prog%wl = 0
         tile%soil%prog%ws = mwc*dz(1:num_l)
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
     if(nfu_inq_var(unit, 'fsc')==NF_NOERR) then 
        call read_tile_data_r1d_fptr(unit,'fsc',soil_fast_soil_C_ptr)
        call read_tile_data_r1d_fptr(unit,'ssc',soil_slow_soil_C_ptr)
     else
        ! try to read fsc and ssc from vegetation restart
        call get_input_restart_name('INPUT/vegn2.res.nc',restart_exists,restart_file_name)
        if (restart_exists) then
           __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit1))
           ! read old (scalar) fsc and ssc into the first element of the fast_soil_C
           ! and slow_soil_C arrays
           call read_tile_data_r1d_fptr(unit1,'fsc',soil_fast_soil_C_ptr,1)
           call read_tile_data_r1d_fptr(unit1,'ssc',soil_slow_soil_C_ptr,1)
        endif
     endif
          
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('soil_init', 'cold-starting soil', NOTE)
  endif

  ! read soil carbon restart, if present
  call get_input_restart_name('INPUT/soil_carbon.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
     call read_tile_data_r1d_fptr(unit,'asoil_in',soil_asoil_in_ptr)
     call read_tile_data_r1d_fptr(unit,'fsc_in',soil_fsc_in_ptr)
     call read_tile_data_r1d_fptr(unit,'ssc_in',soil_ssc_in_ptr)
     __NF_ASRT__(nf_close(unit))     
  endif
  
  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_tau_gw,       lnd%tile_map, soil_tau_groundwater_ptr)
  call send_tile_data_r0d_fptr(id_slope_l,      lnd%tile_map, soil_hillslope_length_ptr)
  call send_tile_data_r0d_fptr(id_slope_Z,      lnd%tile_map, soil_hillslope_relief_ptr)
  call send_tile_data_r0d_fptr(id_zeta_bar,     lnd%tile_map, soil_hillslope_zeta_bar_ptr)
  call send_tile_data_r0d_fptr(id_e_depth,      lnd%tile_map, soil_soil_e_depth_ptr)
  call send_tile_data_r0d_fptr(id_zeta,         lnd%tile_map, soil_zeta_ptr)
  call send_tile_data_r0d_fptr(id_tau,          lnd%tile_map, soil_tau_ptr)
  call send_tile_data_r0d_fptr(id_vwc_wilt,     lnd%tile_map, soil_vwc_wilt_ptr)
  call send_tile_data_r0d_fptr(id_vwc_fc,       lnd%tile_map, soil_vwc_fc_ptr)
  call send_tile_data_r0d_fptr(id_vwc_sat,      lnd%tile_map, soil_vwc_sat_ptr)
  call send_tile_data_r0d_fptr(id_K_sat,        lnd%tile_map, soil_k_sat_ref_ptr)
  call send_tile_data_r0d_fptr(id_K_gw,         lnd%tile_map, soil_k_sat_gw_ptr)
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
  id_fast_soil_C = register_tiled_diag_field ( module_name, 'fast_soil_C', axes,  &
       Time, 'fast soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_slow_soil_C = register_tiled_diag_field ( module_name, 'slow_soil_C', axes,  &
       Time, 'slow soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_fsc = register_tiled_diag_field ( module_name, 'fsc', axes(1:2),  &
       Time, 'total fast soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_ssc = register_tiled_diag_field ( module_name, 'ssc', axes(1:2),  &
       Time, 'total slow soil carbon', 'kg C/m2', missing_value=-100.0 )
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
  id_if  = register_tiled_diag_field ( module_name, 'soil_rif',  axes(1:2),  &
       Time, 'interflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_al  = register_tiled_diag_field ( module_name, 'soil_ral',  axes(1:2),  &
       Time, 'active layer flow',    'kg/(m2 s)',  missing_value=-100.0 )
  id_nu  = register_tiled_diag_field ( module_name, 'soil_rnu',  axes(1:2),  &
       Time, 'numerical runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_sc  = register_tiled_diag_field ( module_name, 'soil_rsc',  axes(1:2),  &
       Time, 'lm2 groundwater runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'soil_hie',  axes(1:2), &
       Time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'soil_hsn',  axes(1:2), &
       Time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'soil_hbf',  axes(1:2), &
       Time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )
  id_hif  = register_tiled_diag_field ( module_name, 'soil_hif',  axes(1:2), &
       Time, 'heat if runf',            'W/m2',  missing_value=-100.0 )
  id_hal  = register_tiled_diag_field ( module_name, 'soil_hal',  axes(1:2), &
       Time, 'heat al runf',            'W/m2',  missing_value=-100.0 )
  id_hnu  = register_tiled_diag_field ( module_name, 'soil_hnu',  axes(1:2), &
       Time, 'heat nu runoff',          'W/m2',  missing_value=-100.0 )
  id_hsc  = register_tiled_diag_field ( module_name, 'soil_hsc',  axes(1:2), &
       Time, 'heat sc runoff',          'W/m2',  missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'soil_evap',  axes(1:2), &
       Time, 'soil evap',            'kg/(m2 s)',  missing_value=-100.0 )
  id_excess  = register_tiled_diag_field ( module_name, 'sfc_excess',  axes(1:2),  &
       Time, 'sfc excess pushed down',    'kg/(m2 s)',  missing_value=-100.0 )

  id_uptk_n_iter  = register_tiled_diag_field ( module_name, 'uptake_n_iter',  axes(1:2), &
       Time, 'number of iterations for soil uptake',  missing_value=-100.0 )
  id_uptk = register_tiled_diag_field ( module_name, 'soil_uptk', axes, &
       Time, 'uptake of water by roots', 'kg/(m2 s)',  missing_value=-100.0 )
  id_psi_x0 = register_tiled_diag_field ( module_name, 'soil_psix0', axes(1:2), &
       Time, 'xylem potential at z=0', 'm',  missing_value=-100.0 )
  id_deficit = register_tiled_diag_field ( module_name, 'soil_def', axes(1:2), &
       Time, 'groundwater storage deficit', '-',  missing_value=-100.0 )
  id_deficit_2 = register_tiled_diag_field ( module_name, 'soil_def2', axes(1:2), &
       Time, 'groundwater storage deficit2', '-',  missing_value=-100.0 )
  id_deficit_3 = register_tiled_diag_field ( module_name, 'soil_def3', axes(1:2), &
       Time, 'groundwater storage deficit3', '-',  missing_value=-100.0 )
  id_psi_bot = register_tiled_diag_field ( module_name, 'soil_psi_n', axes(1:2), &
       Time, 'psi at bottom of soil column', 'm',  missing_value=-100.0 )
  id_sat_frac = register_tiled_diag_field ( module_name, 'soil_fsat', axes(1:2), &
       Time, 'fraction of soil area saturated at surface', '-',  missing_value=-100.0 )
  id_stor_frac = register_tiled_diag_field ( module_name, 'soil_fgw', axes(1:2), &
       Time, 'groundwater storage frac above base elev', '-',  missing_value=-100.0 )
  id_sat_depth = register_tiled_diag_field ( module_name, 'soil_wtdep', axes(1:2), &
       Time, 'depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_sat_dept2 = register_tiled_diag_field ( module_name, 'soil_wtdp2', axes(1:2), &
       Time, 'alt depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_z_cap = register_tiled_diag_field ( module_name, 'soil_zcap', axes(1:2), &
       Time, 'depth below sfc to capillary fringe', 'm',  missing_value=-100.0 )

  id_div_bf = register_tiled_diag_field ( module_name, 'soil_dvbf', axes, &
       Time, 'baseflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_if = register_tiled_diag_field ( module_name, 'soil_dvif', axes, &
       Time, 'interflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_al = register_tiled_diag_field ( module_name, 'soil_dval', axes, &
       Time, 'active-layer flow by layer', 'kg/(m2 s)',  missing_value=-100.0 )

  id_cf_1 = register_tiled_diag_field ( module_name, 'soil_cf_1', axes(1:2), &
       Time, 'soil_cf_1', 'm',  missing_value=-100.0 )
  id_cf_3 = register_tiled_diag_field ( module_name, 'soil_cf_3', axes(1:2), &
       Time, 'soil_cf_1', 'm',  missing_value=-100.0 )
  id_wt_1 = register_tiled_diag_field ( module_name, 'soil_wt_1', axes(1:2), &
       Time, 'soil_wt_1', 'm',  missing_value=-100.0 )
  id_wt_2 = register_tiled_diag_field ( module_name, 'soil_wt_2', axes(1:2), &
       Time, 'soil_wt_2', 'm',  missing_value=-100.0 )
  id_wt_2a = register_tiled_diag_field ( module_name, 'soil_wt_2a', axes(1:2), &
       Time, 'soil_wt_2a', 'm',  missing_value=-100.0 )
  id_wt_3 = register_tiled_diag_field ( module_name, 'soil_wt_3', axes(1:2), &
       Time, 'soil_wt_3', 'm',  missing_value=-100.0 )
  id_wt2_3 = register_tiled_diag_field ( module_name, 'soil_wt2_3', axes(1:2), &
       Time, 'soil_wt2_3', 'm',  missing_value=-100.0 )

  id_active_layer = register_tiled_diag_field ( module_name, 'soil_alt', axes(1:2), &
       Time, 'active-layer thickness', 'm',  missing_value=-100.0 )
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
  id_zeta = register_tiled_static_field ( module_name, 'soil_zeta',      &
       axes(1:2), 'soil depth/topo relief', '-',  missing_value=-100.0 )
  id_tau = register_tiled_static_field ( module_name, 'soil_tau',        &
       axes(1:2), 'gw transmissivity/soil transmissivity', '-',  missing_value=-100.0 )
  id_vwc_wilt = register_tiled_static_field ( module_name, 'soil_wilt',  &
       axes(1:2), 'wilting water content', '-', missing_value=-100.0 )
  id_vwc_fc = register_tiled_static_field ( module_name, 'soil_fc',  &
       axes(1:2), 'field capacity', '-', missing_value=-100.0 )
  id_vwc_sat = register_tiled_static_field ( module_name, 'soil_sat',  &
       axes(1:2), 'soil porosity', '-', missing_value=-100.0 )
  id_K_sat = register_tiled_static_field ( module_name, 'soil_Ksat',  &
       axes(1:2), 'soil sat. hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_K_gw  = register_tiled_static_field ( module_name, 'soil_K_gw',  &
       axes(1:2), 'deep hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
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
  call add_tiled_static_field_alias ( id_slope_Z, module_name, 'slope_Z',  &
       axes(1:2), 'hillslope relief (obsolete, use "soil_rlief" instead)',&
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_e_depth, module_name, 'e_depth',  &
       axes(1:2), 'soil e-folding depth (obsolete, use "soil_depth" instead)', &
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_wilt, module_name, 'vwc_wilt',  &
       axes(1:2), 'wilting water content (obsolete, use "soil_wilt" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_fc, module_name, 'vwc_fc',  &
       axes(1:2), 'field capacity (obsolete, use "soil_fc" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_sat, module_name, 'vwc_sat',  &
       axes(1:2), 'soil porosity (obsolete, use "soil_sat")', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_K_sat, module_name, 'K_sat',  &
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
  call write_tile_data_r1d_fptr(unit,'fsc',          soil_fast_soil_C_ptr,'zfull','fast soil carbon', 'kg C/m2')
  call write_tile_data_r1d_fptr(unit,'ssc',          soil_slow_soil_C_ptr,'zfull','slow soil carbon', 'kg C/m2')
  
  ! close file
  __NF_ASRT__(nf_close(unit))

  if (write_soil_carbon_restart) then
     call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'soil_carbon.res.nc', &
          lnd%coord_glon, lnd%coord_glat, soil_tile_exists, tile_dim_length )
     ! in addition, define vertical coordinate
     if (mpp_pe()==lnd%io_pelist(1)) then
        __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
        __NF_ASRT__(nfu_put_att(unit,'zfull','positive','down'))
     endif
     call sync_nc_files(unit)

     call write_tile_data_r1d_fptr(unit,'asoil_in',soil_asoil_in_ptr,'zfull','aerobic activity modifier', 'unitless')
     call write_tile_data_r1d_fptr(unit,'fsc_in',soil_fsc_in_ptr,'zfull','fast soil carbon input', 'kg C/m2')
     call write_tile_data_r1d_fptr(unit,'ssc_in',soil_ssc_in_ptr,'zfull','slow soil carbon input', 'kg C/m2')
     __NF_ASRT__(nf_close(unit))
  endif


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
  call check_var_range(soil_refl_dir(BAND_VIS), 0.0, 1.0, 'soil_radiation', 'soil_refl_dir(VIS)', lnd%time, FATAL)
  call check_var_range(soil_refl_dir(BAND_NIR), 0.0, 1.0, 'soil_radiation', 'soil_refl_dir(NIR)', lnd%time, FATAL)
  call check_var_range(soil_refl_dif(BAND_VIS), 0.0, 1.0, 'soil_radiation', 'soil_refl_dif(VIS)', lnd%time, FATAL)
  call check_var_range(soil_refl_dif(BAND_NIR), 0.0, 1.0, 'soil_radiation', 'soil_refl_dif(NIR)', lnd%time, FATAL)
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
       VRL, & ! vertical distribution of volumetric root length, m/m3
       u, du ! uptake and its derivative (the latter is not used)
  real :: psi_for_rh
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

  ! calculate relative humidity at soil surface
  call soil_data_psi_for_rh ( soil, vlc, vsc, soil%psi, psi_for_rh )
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
     dt_e = aaa(num_l)*(soil%prog(num_l)%T - soil%prog(num_l-1)%T) &
               + soil%geothermal_heat_flux * delta_time / heat_capacity(num_l)
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
  real, dimension(num_l)   :: del_t, psi, DThDP, hyd_cond, DKDP, &
                                 vlc, vsc, dW_l, DPsi
  real, dimension(num_l+1) :: flow, infilt
  real, dimension(num_l  ) :: div, div_bf, div_if, div_al, dq, div_active
    real      :: &
     lprec_eff, hlprec_eff, tflow, hcap,cap_flow, dheat, &
     melt_per_deg, melt, &
     lrunf_sn,lrunf_ie,lrunf_bf,lrunf_if,lrunf_al,lrunf_nu,lrunf_sc, d_GW, &
     hlrunf_sn,hlrunf_ie,hlrunf_bf,hlrunf_if,hlrunf_al,hlrunf_nu,hlrunf_sc, &
     c0, c1, c2, Dpsi_min, Dpsi_max, &
     sat_area_frac, sat_thick, sum_trans, &
     gw_flux, depth_to_wt, depth_to_wt2_3, depth_to_wt_apparent, &
     depth_to_gw_flow_3, deficit, z_bot, &
     active_layer_thickness, depth_to_cf, d_psi, d_psi_s, psi_star, &
     depth_to_cf_1, depth_to_cf_3, &
     depth_to_wt_1, depth_to_wt_2, depth_to_wt_2a, depth_to_wt_3, &
     storage_2, deficit_2, deficit_3
  logical :: stiff
  real :: zimh, ziph, dTr_g(num_l), dTr_s(num_l)
  integer :: n_iter, l, l_max_active_layer
  real :: &
       VRL(num_l), & ! volumetric root length, m/m3
       K_r, & ! root membrame permeability, kg/(m3 s)
       r_r, & ! root radius, m
       uptake(num_l),   & ! uptake by roots per layer, kg/(m2 s)
       uptake_tot,      & ! total uptake, kg/(m2 s)
       uptake_pos,      & ! sum of the positive uptake, kg/(m2 s) 
       uptake_T_new, & ! updated average temperature of uptaken water, deg K
       uptake_T_corr,& ! correction for uptake temperature, deg K
       Tu,           & ! temperature of water taken up from (or added to) a layer, deg K
       psi_x0          ! water potential inside roots (in xylem) at zero depth, m

  ! --------------------------------------------------------------------------

  
  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 1 #####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'subs_evap    ', subs_evap
     write(*,*) 'snow_lprec   ', snow_lprec
     write(*,*) 'uptake  ', vegn_uptk
     write(*,*) 'subs_M_imp   ', subs_M_imp
     write(*,*) 'theta_s ', soil%pars%vwc_sat
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') 'level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws+soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws,&
             ' gw=', soil%prog(l)%groundwater
     enddo
    endif
  !.........................................................................

  ! ---- record fluxes -----------------------------------------------------
  soil_levap  = subs_evap*(1-soil_subl)
  soil_fevap  = subs_evap*   soil_subl
  soil_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution ------------
  del_t(1) = subs_DT
  soil%prog(1)%T = soil%prog(1)%T + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = soil%e(l) * del_t(l) + soil%f(l)
      soil%prog(l+1)%T = soil%prog(l+1)%T + del_t(l+1)
      enddo
    endif

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 2 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') 'level=',l, 'T=', soil%prog(l)%T, &
             'del_t=', del_t(l), 'e=', soil%e(l), 'f=', soil%f(l)
       enddo
    endif
  !.........................................................................

  ! ---- extract evap from soil, adjusting T, and do implicit melt ---------
  IF (LM2) THEN  ! (extract surface E--is there any?--uniformly from bucket)
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

  ! ---- calculate actual uptake and update its T --------------------------
  select case(uptake_option)
  case ( UPTAKE_LINEAR )
     uptake_T_corr = 0
     n_iter = 0
     uptake = soil%uptake_frac*vegn_uptk
     soil%psi_x0 = 0.
  case ( UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN )     
     ! for Darcy-flow uptake, find the root water potential to satify actual
     ! transpiration by the vegetation
     call vegn_root_properties (vegn, dz(1:num_l), VRL, K_r, r_r)
     
     if ( uptake_option==UPTAKE_DARCY2D ) then
        call darcy2d_uptake_solver ( soil, vegn_uptk, VRL, K_r, r_r, &
             uptake_oneway, uptake_from_sat, uptake, psi_x0, n_iter)
     else
        call darcy2d_uptake_solver_lin ( soil, vegn_uptk, VRL, K_r, r_r, &
             uptake_oneway, uptake_from_sat, uptake, psi_x0, n_iter )
     endif
     soil%psi_x0 = psi_x0

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

  !.........................................................................
  if (is_watch_point())then
      write(*,*) ' ##### soil_step_2 checkpoint 2.1 #####'
      __DEBUG2__(vegn_uptk,sum(uptake))
      do l = 1,num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))')'level=',l, &
             'uptake=',uptake(l),'dwl=',-uptake(l)*delta_time,&
             'wl=',soil%prog(l)%wl,'new wl=',soil%prog(l)%wl - uptake(l)*delta_time
        enddo
    endif
  !.........................................................................

  call send_tile_data(id_uptk_n_iter, real(n_iter), diag)
  call send_tile_data(id_uptk, uptake, diag)
  call send_tile_data(id_psi_x0, psi_x0, diag)

  ! ---- perform the uptake ------------------------------------------------
  do l = 1, num_l
    ! calculate the temperature of water that is taken from the layer (or added 
    ! to the layer), including energy balance correction 
    if (uptake(l) > 0) then
        Tu = soil%prog(l)%T + uptake_T_corr
      else
        Tu = soil%uptake_T + uptake_T_corr
      endif
    hcap = soil%heat_capacity_dry(l)*dz(l) &
          + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws
    soil%prog(l)%T = soil%prog(l)%T - &
          uptake(l)*delta_time*clw*( Tu-soil%prog(l)%T ) / &
          ( hcap - uptake(l)*delta_time*clw )
    soil%prog(l)%wl = soil%prog(l)%wl - uptake(l)*delta_time
    enddo

  !.........................................................................
  if(is_watch_point()) then
      write(*,*) ' ##### soil_step_2 checkpoint 3 #####'
      do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws+soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws
        enddo
    endif
  !.........................................................................

  ! ---- push down any excess surface water, with heat ---------------------
  IF (PUSH_DOWN_SFC_EXCESS) THEN
      CALL SOIL_PUSH_DOWN_EXCESS ( soil, diag, lrunf_nu, hlrunf_nu )
    ELSE
      lrunf_nu=0; hlrunf_nu=0
    ENDIF

  ! ---- fetch soil hydraulic properties -----------------------------------
  do l = 1, num_l
    vlc(l) = max(0., soil%prog(l)%wl / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%prog(l)%ws / (dens_h2o*dz(l)))
  enddo
  call soil_data_hydraulic_properties (soil, vlc, vsc, &
                   psi, DThDP, hyd_cond, DKDP, Dpsi_min, Dpsi_max )

  ! ---- compute various measures of water table depth ---------------------
    sat_thick = 0.
    do l=num_l,1,-1
       if(vsc(l)+vlc(l).le.soil%pars%vwc_sat) exit
       sat_thick = sat_thick + dz(l)
    enddo
    depth_to_cf_1 = zhalf(num_l+1) - sat_thick
    depth_to_wt_1 = depth_to_cf_1 - soil%pars%psi_sat_ref/soil%alpha(max(l,1))

    depth_to_wt_2 = zfull(num_l)-psi(num_l)

    depth_to_wt_2a = 0.
    do l=1,num_l
      if (soil%prog(l)%wl+soil%prog(l)%ws .lt. &
                       soil%pars%vwc_sat*dens_h2o*dz(l)) then
          depth_to_wt_2a = depth_to_wt_2a + dz(l)
          if (l.eq.num_l) depth_to_wt_2a = -1.
      else
          exit
      endif
    enddo
    storage_2 = 1 - depth_to_wt_2  &
             /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)
    storage_2 = min( max( 0., storage_2 ) , 1.)
    deficit_2 = 1 - storage_2

        if (vsc(num_l).gt.0.) then   ! permafrost
	    depth_to_wt2_3 = 0.
	    depth_to_cf_3 = 0.
            depth_to_wt_3 = 0.
	  else                       ! liquid water at depth
	    depth_to_cf_3 = 0.
            if (use_fringe) then
	        do l = num_l, 1, -1
		  if ( l.eq.num_l .and. psi(l).le.soil%pars%psi_sat_ref/soil%alpha(l) ) then
		      depth_to_cf_3 = zfull(l) + soil%pars%psi_sat_ref/soil%alpha(l) - psi(l)
		      exit
		    else if (psi(l).le.soil%pars%psi_sat_ref/soil%alpha(l) ) then
		      d_psi = psi(l+1) - psi(l)
		      d_psi_s = (soil%pars%psi_sat_ref/soil%alpha(l+1)) &
		               -(soil%pars%psi_sat_ref/soil%alpha(l))
	              psi_star = (psi(l)*d_psi_s - &
		                  d_psi*(soil%pars%psi_sat_ref/soil%alpha(l)))&
				  / (d_psi_s - d_psi)
		      depth_to_cf_3 = zfull(l) + (zfull(l+1)-zfull(l)) &
		                            * (psi_star-psi(l)) / d_psi
		      exit
		    else if (l.eq.1) then
		      d_psi = psi(l+1) - psi(l)
		      d_psi_s = (soil%pars%psi_sat_ref/soil%alpha(l+1)) &
		               -(soil%pars%psi_sat_ref/soil%alpha(l))
	              psi_star = (psi(l)*d_psi_s - &
		                  d_psi*(soil%pars%psi_sat_ref/soil%alpha(l)))&
				  / (d_psi_s - d_psi)
		      depth_to_cf_3 = zfull(l) + (zfull(l+1)-zfull(l)) &
		                            * (psi_star-psi(l)) / d_psi
		      depth_to_cf_3 = max(0.,depth_to_cf_3)
		    endif
		  enddo
	      endif
            depth_to_wt_3 = max(0., zhalf(num_l+1)-(psi(num_l)+dz(num_l)/2.))
	    depth_to_wt2_3 = depth_to_wt_3
	  endif

	    if (use_fringe) then
	        depth_to_gw_flow_3 = depth_to_cf_3
              else
	        depth_to_gw_flow_3 = depth_to_wt_3
	      endif
        deficit_3 = depth_to_gw_flow_3  &
             /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)

  ! ---- get saturated area and column flow divergences --------------------
    SELECT CASE(gw_option)
    
    CASE(GW_LM2)

        div_bf=0; div_if=0; div_al=0; sat_area_frac = 0

    CASE(GW_LINEAR)

        IF (CORRECTED_LM2_GW) THEN
            do l = 1, num_l
              if (vlc(l) .ge. soil%pars%vwc_sat .and. vsc(l).le.0.) &
                  div_bf(l) = 0.15*dens_h2o*dz(l)/soil%pars%tau_groundwater
              enddo
          ELSE
            do l = 1, num_l
              if ((vsc(l)+vlc(l)) .ge. soil%pars%vwc_sat) &
                  div_bf(l) = 0.15*dens_h2o*dz(l)*(vlc(l)/(vsc(l)+vlc(l)))  &
		                      /soil%pars%tau_groundwater
              enddo
          ENDIF
        div_if = 0
	div_al = 0
	sat_thick = zhalf(num_l+1) - depth_to_cf_1
        sat_area_frac = min((sat_thick/zhalf(num_l+1))**soil%pars%rsa_exp,1.)

    CASE(GW_HILL_AR5)

	call soil_data_gw_hydraulics_ar5(soil, storage_2, &
                                                 gw_flux, sat_area_frac)
        dq = 0.
        sat_thick = 0.
        do l=num_l,1,-1
          if(psi(l).le.0.) exit
          if (vsc(l).le.0.) dq(l) = dz(l)
          sat_thick = sat_thick + dz(l)
          enddo
        div_bf = 0.
        if (sat_thick.gt.0.) div_bf = (dq/sat_thick)*gw_flux

        div_active = 0.
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

        div_al = 0
        where (div_bf.eq.0.) div_al = div_active*active_layer_drainage_acceleration
        div_if = 0

    CASE(GW_HILL)
      
        if (vsc(num_l).gt.0.) then   ! permafrost
	    sat_area_frac = 0.
            div_bf = 0.
            div_if = 0.
	  else                       ! liquid water at depth
            call soil_data_gw_hydraulics(soil, deficit_3, &
                                               gw_flux, sat_area_frac)
            gw_flux = min(gw_flux, gw_flux_max)
             dTr_g = 0.
             dTr_s = 0.
             dTr_g(num_l) = 1.
             l = num_l
             ziph = sum(dz)
             zimh = ziph - dz(num_l)
             if (depth_to_gw_flow_3 .lt. zimh) then
                 dTR_g(l) = dz(l)
                 dTr_s(l) = (exp(-zimh/soil%pars%soil_e_depth))
                 do l = num_l-1, 1, -1
                   if (vsc(l).gt.0.) exit
                       ziph = zimh
                       zimh = ziph - dz(l)
                       if (depth_to_gw_flow_3 .lt. zimh) then
                           dTR_g(l) = dz(l)
                           dTr_s(l) = exp(-zimh/soil%pars%soil_e_depth)-exp(-ziph/soil%pars%soil_e_depth)
                         else if (depth_to_gw_flow_3 .lt. ziph) then
                           dTR_g(l) =(ziph-depth_to_gw_flow_3)
                           dTr_s(l) = exp(-depth_to_gw_flow_3/soil%pars%soil_e_depth)-exp(-ziph/soil%pars%soil_e_depth)
                         else
                           exit
                         endif
                   enddo
               endif
             sum_trans = sum(dTr_g)
             if (sum_trans.ne.0.) then
                 dTR_g = dTR_g / sum_trans
                 dTR_g = dTR_g * soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length
               endif
             dTR_s = dTR_s * (soil%pars%k_sat_sfc+k_macro_constant)*soil%pars%soil_e_depth
             sum_trans = sum(dTR_g) + sum(dTr_s)
             if (sum_trans.ne.0.) then
                 div_bf = gw_flux * dTR_g /sum_trans
                 div_if = gw_flux * dTR_s /sum_trans
               else
                 div_bf = 0.
                 div_if = 0.
               endif
	  endif

	div_al = 0
        l_max_active_layer = 0   ! "active layer" either over permafrost or perched
        do l=1,num_l
          if(vsc(l).gt.0.) exit
          l_max_active_layer = l
          enddo
        if (l_max_active_layer.lt.num_l .and. l_max_active_layer.gt.0) then
            do l = 1, l_max_active_layer
              div_al(l) = hyd_cond(l) * soil%pars%hillslope_relief*dz(l) &
                   / (soil%pars%hillslope_length*soil%pars%hillslope_length)
              enddo
          endif

    END SELECT

    div = div_bf + div_if + div_al
    lrunf_bf = sum(div_bf)
    lrunf_if = sum(div_if)
    lrunf_al = sum(div_al)

    if (snow_lprec.ne.0.) then
      lrunf_sn = sat_area_frac * snow_lprec
      hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
    else
      lrunf_sn = 0.
      hlrunf_sn = 0.
    endif
    hlrunf_ie=0
    lprec_eff = snow_lprec - lrunf_sn
    hlprec_eff = snow_hlprec - hlrunf_sn

    if(is_watch_point()) then
       do l = 1, num_l
          write(*,'(a,1x,i2.2,100(2x,g23.16))')'div_ac,div_bf,div_if,div_al,div', &
	                 l,div_active(l),div_bf(l),div_if(l),div_al(l),div(l)
       enddo
       do l = 1, num_l
          write(*,'(a,1x,i2.2,100(2x,g23.16))')'vsc,psi,dz',l,vsc(l),psi(l),dz(l)
       enddo
       write(*,*)'lrunf_bf',lrunf_bf
       write(*,*)'tau_gw',soil%pars%tau_groundwater
       write(*,*)'dens_h2o',dens_h2o
    endif

  ! ---- soil-water flow ----------------------------------------------------
  IF (LM2) THEN
      flow(1) = 0
      do l = 1, num_l
        infilt(l) = soil%uptake_frac(l)*lprec_eff *delta_time
        flow(l+1) = max(0., soil%prog(l)%wl + flow(l) &
              + infilt(l) - soil%w_fc(l)*dz(l)*dens_h2o)
        dW_l(l) = flow(l) - flow(l+1) + infilt(l)
        soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
        enddo
      do l = 1, num_l
        flow(l) = flow(l) + infilt(l)
        enddo
      dW_l=0
      dpsi=0
      c0 = delta_time/soil%pars%tau_groundwater
      c1 = exp(-c0)
      c2 = (1-c1)/c0
      l = 1
      d_GW = c1 * soil%prog(l)%groundwater + c2 * flow(num_l+1) &
	                    - soil%prog(l)%groundwater
      soil%prog(l)%groundwater = soil%prog(l)%groundwater + d_GW
      lrunf_sc  = (1-c1)*soil%prog(l)%groundwater/delta_time &
                          + (1-c2)*flow(num_l+1)/delta_time
      lrunf_ie=0
    ELSE
      lrunf_sc = 0
      d_GW = 0
      stiff = all(DThDP.eq.0)
      IF(stiff .AND. BYPASS_RICHARDS_WHEN_STIFF) THEN
          flow = 0.
          dW_l = 0.
          div  = 0; div_bf=0; div_if=0; div_al=0
          lrunf_bf = 0; lrunf_if = 0; lrunf_al = 0
          lrunf_ie = lprec_eff
          hlrunf_ie = hlprec_eff
          if (use_stiff_bug) then
              psi=-zfull
            else
              psi=zfull
            endif
          dpsi=0.
        ELSE
          CALL RICHARDS(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, stiff, &
                 dPsi, dW_l, flow, lrunf_ie)
        ENDIF
    ENDIF

  ! ---- heat advection by water flow ---------------------------------------
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

  call advection(soil, flow, dW_l, tflow, d_GW, snow_lprec, snow_hlprec)

  if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
      hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
    else if (flow(1).lt.0. ) then
      hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(soil%prog(1)%T-tfreeze)
    else
      hlrunf_ie = 0.
    endif

  hlrunf_bf = clw*sum(div_bf*(soil%prog%T-tfreeze))
  hlrunf_if = clw*sum(div_if*(soil%prog%T-tfreeze))
  hlrunf_al = clw*sum(div_al*(soil%prog%T-tfreeze))
  hlrunf_sc = clw*lrunf_sc  *(soil%prog(1)%groundwater_T-tfreeze)
  if (lrunf_from_div) then
      soil_lrunf  =  lrunf_sn +  lrunf_ie +  sum(div) +  lrunf_nu +  lrunf_sc
      soil_hlrunf = hlrunf_sn + hlrunf_ie +  clw*sum(div*(soil%prog%T-tfreeze)) &
                                                      + hlrunf_nu + hlrunf_sc
    else
      soil_lrunf  =  lrunf_sn +  lrunf_ie +  lrunf_bf +  lrunf_if &
                              +  lrunf_al +  lrunf_nu +  lrunf_sc
      soil_hlrunf = hlrunf_sn + hlrunf_ie + hlrunf_bf + hlrunf_if &
                              + hlrunf_al + hlrunf_nu + hlrunf_sc
    endif

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = soil%heat_capacity_dry(l)*dz(l) &
             + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws
    melt_per_deg = hcap/(hlf_factor*hlf)
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
       + (hcap*(soil%prog(l)%T-tfreeze) - hlf_factor*hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
    soil_melt = soil_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 5 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws,&
             ' gw=', soil%prog(l)%groundwater
     enddo
  endif

  active_layer_thickness = 0.
  do l = 1, num_l
    if (soil%prog(l)%ws.gt.0.) then
        active_layer_thickness = active_layer_thickness &
	  + dz(l)*soil%prog(l)%wl/(soil%prog(l)%wl+soil%prog(l)%ws)
        exit
      endif
    active_layer_thickness = active_layer_thickness + dz(l)
    enddo

  soil_Ttop = soil%prog(1)%T
  soil_Ctop = soil%heat_capacity_dry(1)*dz(1) &
    + clw*soil%prog(1)%wl + csw*soil%prog(1)%ws

if (update_psi) soil%psi=psi+dPsi

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
!    call send_tile_data(id_deficit, deficit, diag)
!    call send_tile_data(id_sat_depth, depth_to_wt_3, diag)
!    call send_tile_data(id_sat_dept2, depth_to_wt2_3, diag)
!    call send_tile_data(id_z_cap, depth_to_cf_3, diag)
!    if (depth_to_wt_2a .ge. -0.5) &
!	                call send_tile_data(id_sat_depth, depth_to_wt_2a, diag)
    call send_tile_data(id_cf_1, depth_to_cf_1, diag)
    call send_tile_data(id_cf_3, depth_to_cf_3, diag)
    call send_tile_data(id_wt_1, depth_to_wt_1, diag)
    call send_tile_data(id_wt_2, depth_to_wt_2, diag)
    call send_tile_data(id_wt_2a, depth_to_wt_2a, diag)
    call send_tile_data(id_wt_3, depth_to_wt_3, diag)
    call send_tile_data(id_wt2_3, depth_to_wt2_3, diag)
    call send_tile_data(id_deficit_2, deficit_2, diag)
    call send_tile_data(id_deficit_3, deficit_3, diag)
    call send_tile_data(id_sat_frac, sat_area_frac, diag)
    call send_tile_data(id_div_bf, div_bf, diag)
    call send_tile_data(id_div_if, div_if, diag)
    call send_tile_data(id_div_al, div_al, diag)

   call send_tile_data(id_ie,   lrunf_ie, diag)
   call send_tile_data(id_sn,   lrunf_sn, diag)
   call send_tile_data(id_bf,   lrunf_bf, diag)
   call send_tile_data(id_if,   lrunf_if, diag)
   call send_tile_data(id_al,   lrunf_al, diag)
   call send_tile_data(id_nu,   lrunf_nu, diag)
   call send_tile_data(id_sc,   lrunf_sc, diag)
   call send_tile_data(id_hie,  hlrunf_ie, diag)
   call send_tile_data(id_hsn,  hlrunf_sn, diag)
   call send_tile_data(id_hbf,  hlrunf_bf, diag)
   call send_tile_data(id_hif,  hlrunf_if, diag)
   call send_tile_data(id_hal,  hlrunf_al, diag)
   call send_tile_data(id_hnu,  hlrunf_nu, diag)
   call send_tile_data(id_hsc,  hlrunf_sc, diag)
   if (id_evap > 0) call send_tile_data(id_evap,  soil_levap+soil_fevap, diag)

   call send_tile_data(id_heat_cap, soil%heat_capacity_dry, diag)
   call send_tile_data(id_active_layer, active_layer_thickness, diag)

end subroutine soil_step_2

! ============================================================================

subroutine soil_step_3(soil, diag)
  type(soil_tile_type), intent(in) :: soil
  type(diag_buff_type), intent(inout) :: diag

  if(id_fast_soil_C>0) call send_tile_data(id_fast_soil_C, soil%fast_soil_C(:)/dz(:), diag)
  if(id_slow_soil_C>0) call send_tile_data(id_slow_soil_C, soil%slow_soil_C(:)/dz(:), diag)
  if (id_fsc > 0)      call send_tile_data(id_fsc, sum(soil%fast_soil_C(:)), diag)
  if (id_ssc > 0)      call send_tile_data(id_ssc, sum(soil%slow_soil_C(:)), diag)
end subroutine soil_step_3

! ============================================================================
  
  subroutine soil_push_down_excess ( soil, diag, lrunf_nu, hlrunf_nu )
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: lrunf_nu, hlrunf_nu

  ! ---- local vars ----------------------------------------------------------
  real      :: &
     liq_frac, excess_wat, excess_liq, excess_ice, excess_t, &
     h1, h2, summax, space_avail, liq_placed, ice_placed
  integer :: l

  liq_frac=0;excess_wat=0;excess_liq=0;excess_ice=0;h1=0;h2=0
  l = 1
  summax = max(0.,soil%prog(l)%wl)+max(0.,soil%prog(l)%ws)
  if (summax > 0) then
     liq_frac = max(0.,soil%prog(l)%wl) / summax
  else
     liq_frac = 1
  endif
  excess_wat = max(0., soil%prog(l)%wl + soil%prog(l)%ws &
       - dens_h2o*dz(l)*soil%pars%vwc_sat )
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
          ' soil%pars%vwc_sat =', soil%pars%vwc_sat,&
          ' excess_liq =', excess_liq,&
          ' excess_ice =', excess_ice, &
          ' dens_h2o=', dens_h2o, &
          ' dz(l)=',dz(l)
    endif

  do l = 2, num_l
     if (excess_liq+excess_ice>0) then
        space_avail = dens_h2o*dz(l)*soil%pars%vwc_sat &
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
               - hlf_factor*hlf*excess_ice                   ) / delta_time

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.01 #####'
     write(*,*) ' lrunf_nu',lrunf_nu
     write(*,*) 'hlrunf_nu',hlrunf_nu
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,100(x,a,g23.16))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws
     enddo
  endif
end subroutine soil_push_down_excess

! ============================================================================

  subroutine RICHARDS(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, stiff, &
                 dPsi, dW_l, flow, lrunf_ie)
  type(soil_tile_type), intent(inout)   :: soil
  real, intent(in),  dimension(num_l)   :: psi, DThDP, hyd_cond, DKDP, div
  real, intent(in)                      :: lprec_eff, Dpsi_min, Dpsi_max
  logical, intent(in)                   :: stiff
  real, intent(out), dimension(num_l)   :: dPsi, dW_l
  real, intent(out), dimension(num_l+1) :: flow
  real, intent(out)                     :: lrunf_ie
  ! ---- local vars ----------------------------------------------------------
  integer l, ipt, jpt, kpt, fpt, l_internal
  real, dimension(num_l-1) :: del_z, K, DKDPm, DKDPp, grad, eee, fff
  real aaa, bbb, ccc, ddd, xxx, dpsi_alt, dW_l_internal, w_to_move_up, adj
  logical flag

    flag = .false.
    flow(1) = delta_time*lprec_eff
    do l = 1, num_l-1
      del_z(l) = zfull(l+1)-zfull(l)
      K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
      DKDPm(l) = 0. !0.5*DKDP(l)
      DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = (psi(l+1)-psi(l))/del_z(l) - 1
    enddo

    if(is_watch_point()) then
       write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
       do l = 1, num_l
          write(*,'(x,a,x,i2.2,x,a,100(x,g23.16))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
               DThDP(l),&
               hyd_cond(l),&
               psi(l),&
               DKDP(l)
       enddo
       do l = 1, num_l-1
          write(*,'(a,i2.2,1x,a,100(2x,g23.16))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
               K(l),&
               DKDPm(l),&
               DKDPp(l),&
               grad(l)
       enddo
    endif


    l = num_l
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    aaa =     - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
!      where (stiff)
    bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
    ddd = - K(l-1) *grad(l-1) - div(l)
!        elsewhere
!          Qout = hyd_cond(l) ! gravity drainage
!          DQoutDP = DKDP(l)  ! gravity drainage
!          Qout = 0.                ! no drainage
!          DQoutDP = 0.             ! no drainage
!          where (psi(l).gt.0.) ! linear baseflow from gw
!              Qout = 0.15*psi(l)/soil%pars%tau_groundwater
!              DQoutDP = 0.15/soil%pars%tau_groundwater
!            elsewhere
!              Qout = 0.
!              DQoutDP = 0.
!            endwhere
!          bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
!                      -DQoutDP )
!          ddd = -Qout - K(l-1) *grad(l-1)
!        endwhere
    eee(l-1) = -aaa/bbb
    fff(l-1) =  ddd/bbb
  
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b, ,d', l,aaa, bbb,ddd
    endif

    do l = num_l-1, 2, -1
      xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
      aaa = - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
      bbb = xxx-( -K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                  -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
      ccc =   - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
      ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                            - div(l)
      eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
      fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      if(is_watch_point()) then
         write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
      endif
    enddo
  
    l = 1
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    bbb = xxx - ( -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =     - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =          flow(1)/delta_time +    K(l)     *grad(l) &
                            - div(l)

IF (bbb+ccc*eee(l) .NE. 0.) THEN
    dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
    if(is_watch_point()) then
       write(*,*) 'bbb+ccc*eee(l) .NE. 0.'
       write(*,*) 'bbb', bbb
       write(*,*) 'ccc', ccc
       write(*,*) 'ddd', ddd
       write(*,*) 'eee(l)', eee(l)
       write(*,*) 'fff(l)', fff(l)
       write(*,*) 'dPsi(l)', dPsi(l)
       write(*,*) 'stiff', stiff
       write(*,*) 'dPsi(l)', dPsi(l)
       write(*,*) 'Dpsi_min', Dpsi_min
       write(*,*) 'Dpsi_max', Dpsi_max
       endif
    if (.not.stiff .and. dPsi(l).gt.Dpsi_min .and. dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
      else
      	if (stiff) then
            dPsi(l) = - psi(l)
          else
	    if (dPsi(l).lt.Dpsi_min) then
                flag = .true.
                call get_current_point(ipt,jpt,kpt,fpt)
                write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'computed dPsi(1) too negative'
                write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi(1)=',dPsi(l)
                write(*,*) 'note 1: at point ',ipt,jpt,kpt,fpt,'dPsi_min =',dPsi_min
              endif
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
    if(is_watch_point()) then
       write(*,*) 'bbb+ccc*eee(l) .EQ. 0.'
       write(*,*) 'bbb', bbb
       write(*,*) 'ccc', ccc
       write(*,*) 'ddd', ddd
       write(*,*) 'eee(l)', eee(l)
       write(*,*) 'fff(l)', fff(l)
       write(*,*) 'dPsi(l)', dPsi(l)
       write(*,*) 'stiff', stiff
       write(*,*) 'dPsi(l)', dPsi(l)
       write(*,*) 'Dpsi_min', Dpsi_min
       write(*,*) 'Dpsi_max', Dpsi_max
       endif
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
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,  b,c,d', l, bbb,ccc,ddd
       write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
       write(*,*) 'ie:', lrunf_ie
       do l = 1, num_l-1
          write(*,'(a,i2.2,100(2x,g23.16))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
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
           +(DPsi(l+1)-DPsi(l))/ del_z(l)) &
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
          write(*,'(i2.2,100(2x,a,g23.16))') l,&
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
          write(*,'(i2.2,100(2x,a,g23.16))') l,&
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
               - soil%pars%vwc_sat*dz(l)*dens_h2o, 0. )

          if(is_watch_point()) then
             write(*,*) '3.22 l=', l,&
                  ' soil_prog%wl=',soil%prog(l)%wl,  &
                  ' soil_prog%ws=',soil%prog(l)%ws , &
                  ' soil%pars%vwc_sat=', soil%pars%vwc_sat, &
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

  ENDIF
       
       do l = 1, num_l
         soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
         enddo

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l, &
             'Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)), &
             'wl=', soil%prog(l)%wl, &
             'ws=', soil%prog(l)%ws, &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

end subroutine richards

! ============================================================================
  subroutine advection(soil, flow, dW_l, tflow, d_GW, snow_lprec, snow_hlprec)
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in), dimension(:) :: flow
  real, intent(in), dimension(:) :: dW_l
  real, intent(in) :: &
     tflow, &
     d_GW, &
     snow_lprec, &
     snow_hlprec
  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l)   :: u_minus, u_plus, del_t
  real, dimension(num_l-1) :: eee, fff
  real hcap, aaa, bbb, ccc
  integer l
  
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

  ! (lumped=lm2 groundwater stored in l=1 prog variable, liquid only)
  if (soil%prog(1)%groundwater.ne. 0.) soil%prog(1)%groundwater_T =    &
       + ((aquifer_heat_cap+soil%prog(1)%groundwater-d_GW)  &
	                         *soil%prog(1)%groundwater_T &
        + flow(num_l+1)*soil%prog(num_l)%T) &
         /((aquifer_heat_cap+soil%prog(1)%groundwater-d_GW) + flow(num_l+1))

end subroutine advection

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
DEFINE_SOIL_ACCESSOR_1D(real,fast_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,slow_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,asoil_in)
DEFINE_SOIL_ACCESSOR_1D(real,fsc_in)
DEFINE_SOIL_ACCESSOR_1D(real,ssc_in)

DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau_groundwater)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_length)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_relief)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_zeta_bar)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,soil_e_depth)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,zeta)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_gw)
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

!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,T)
!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,wl)
!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,ws)
!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,groundwater)
!DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,groundwater_T)

! ============================================================================
  subroutine soil_T_ptr(t,p)
  type(land_tile_type),pointer::t
  real,pointer::p(:)
  integer :: n

  p=>NULL()
  if(associated(t))then
    if(associated(t%soil))then
      n=size(t%soil%prog(:))
      p(1:n)=>t%soil%prog(1:n)%T
    endif
  endif
  end subroutine
 
! ============================================================================
  subroutine soil_wl_ptr(t,p)
  type(land_tile_type),pointer::t
  real,pointer::p(:)
  integer :: n

  p=>NULL()
  if(associated(t))then
  if(associated(t%soil))then
  n=size(t%soil%prog(:))
  p(1:n)=>t%soil%prog(1:n)%wl
  endif
  endif
  end subroutine
 
! ============================================================================
  subroutine soil_ws_ptr(t,p)
  type(land_tile_type),pointer::t
  real,pointer::p(:)
  integer :: n

  p=>NULL()
  if(associated(t))then
  if(associated(t%soil))then
  n=size(t%soil%prog(:))
  p(1:n)=>t%soil%prog(1:n)%ws
  endif
  endif
  end subroutine

! ============================================================================
  subroutine soil_groundwater_ptr(t,p)
  type(land_tile_type),pointer::t
  real,pointer::p(:)
  integer :: n

  p=>NULL()
  if(associated(t))then
  if(associated(t%soil))then
  n=size(t%soil%prog(:))
  p(1:n)=>t%soil%prog(1:n)%groundwater
  endif
  endif
  end subroutine

! ============================================================================
  subroutine soil_groundwater_T_ptr(t,p)
  type(land_tile_type),pointer::t
  real,pointer::p(:)
  integer :: n

  p=>NULL()
  if(associated(t))then
  if(associated(t%soil))then
  n=size(t%soil%prog(:))
  p(1:n)=>t%soil%prog(1:n)%groundwater_T
  endif
  endif
  end subroutine

end module soil_mod
