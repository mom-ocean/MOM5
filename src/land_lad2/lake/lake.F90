! ============================================================================
! lake model module
! ============================================================================
module lake_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, read_data, check_nml_error, &
     stdlog, write_version_number, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init, register_diag_field,           &
                              register_static_field, send_data
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o, PI, grav, vonkarm, &
                              rdgas

use land_constants_mod, only : &
     NBANDS
use land_io_mod, only : read_field
use lake_tile_mod, only : &
     lake_tile_type, lake_pars_type, lake_prog_type, read_lake_data_namelist, &
     lake_data_radiation, lake_data_diffusion, &
     lake_data_thermodynamics, &
     max_lev, cpw,clw,csw, lake_width_inside_lake, large_lake_sill_width, &
     lake_specific_width, n_outlet, outlet_face, outlet_i, outlet_j, outlet_width
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=)
use land_tile_diag_mod, only : register_tiled_static_field, &
     register_tiled_diag_field, send_tile_data, diag_buff_type, &
     send_tile_data_r0d_fptr, add_tiled_static_field_alias
use land_data_mod,      only : land_state_type, lnd
use land_tile_io_mod, only : print_netcdf_error, create_tile_out_file, &
     read_tile_data_r1d_fptr, write_tile_data_r1d_fptr, sync_nc_files, &
     get_input_restart_name
use nf_utils_mod, only : nfu_inq_var, nfu_def_dim, nfu_put_att
use land_debug_mod, only: is_watch_point
use land_utils_mod, only : put_to_tiles_r0d_fptr

implicit none
private

! ==== public interfaces =====================================================
public :: read_lake_namelist
public :: lake_init
public :: lake_end
public :: save_lake_restart

public :: lake_get_sfc_temp
public :: lake_radiation
public :: lake_diffusion
public :: lake_step_1
public :: lake_step_2

public :: large_dyn_small_stat
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'lake',&
    version     = '$Id: lake.F90,v 20.0 2013/12/13 23:29:37 fms Exp $',&
    tagname     = '$Name: tikal $'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
real    :: init_temp            = 288.        ! cold-start lake T
real    :: init_w               = 1000.      ! cold-start w(l)/dz(l)
real    :: init_groundwater     =   0.        ! cold-start gw storage
logical :: use_rh_feedback      = .true.
logical :: make_all_lakes_wide  = .false.
logical :: large_dyn_small_stat = .true.
logical :: relayer_in_step_one  = .false.
logical :: float_ice_to_top     = .false.
logical :: wind_penetrates_ice  = .false.
real    :: min_rat              = 0.4
logical :: do_stratify          = .true.
character(len=16):: albedo_to_use = ''  ! or 'brdf-params'
real    :: K_z_large            = 1.
real    :: K_z_background       = 0.
real    :: K_z_min              = 0.
real    :: K_z_factor           = 1.
real    :: c_drag               = 1.2e-3
real    :: lake_depth_max       = 1.e10
real    :: lake_depth_min       = 1.99
real    :: max_plain_slope      = -1.e10

namelist /lake_nml/ init_temp, init_w,       &
                    init_groundwater, use_rh_feedback, cpw, clw, csw, &
                    make_all_lakes_wide, large_dyn_small_stat, &
                    relayer_in_step_one, float_ice_to_top, &
                    min_rat, do_stratify, albedo_to_use, K_z_large, &
		    K_z_background, K_z_min, K_z_factor, &
		    lake_depth_max, lake_depth_min, max_plain_slope
!---- end of namelist --------------------------------------------------------
real    :: K_z_molec            = 1.4e-7
real    :: tc_molec             = 0.59052 ! dens_h2o*clw*K_z_molec
real    :: tc_molec_ice         = 2.5

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf
type(time_type) :: time
real            :: delta_time

integer         :: num_l              ! # of water layers
real, allocatable:: zfull (:)    ! diag axis, dimensionless layer number
real, allocatable:: zhalf (:)
real            :: max_rat

! ---- diagnostic field IDs
integer :: id_lwc, id_swc, id_temp, id_ie, id_sn, id_bf, id_hie, id_hsn, id_hbf
integer :: id_evap, id_dz, id_wl, id_ws, id_K_z, id_silld, id_sillw, id_backw
integer :: id_back1
! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_lake_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_lake_data_namelist(num_l)

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=lake_nml, iostat=io)
     ierr = check_nml_error(io, 'lake_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=lake_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'lake_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=lake_nml)
  endif

  ! ---- set up vertical discretization
  allocate (zhalf(num_l+1), zfull(num_l))
  zhalf(1) = 0
  do l = 1, num_l;   
     zhalf(l+1) = zhalf(l) + 1.
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo
  
  max_rat = 1. / min_rat
  if (trim(albedo_to_use)=='brdf-params') then
     use_brdf = .true.
  else if (trim(albedo_to_use)=='') then
     use_brdf = .false.
  else
     call error_mesg('lake_init',&
          'option albedo_to_use="'// trim(albedo_to_use)//&
          '" is invalid, use "brdf-params", or nothing ("")',&
          FATAL)
  endif

end subroutine read_lake_namelist


! ============================================================================
! initialize lake model
subroutine lake_init ( id_lon, id_lat )
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis

  ! ---- local vars 
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! last and current tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  character(len=256) :: restart_file_name
  logical :: restart_exists
  real, allocatable :: buffer(:,:),bufferc(:,:),buffert(:,:)
  integer i

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

allocate(buffer (lnd%is:lnd%ie,lnd%js:lnd%je))
allocate(bufferc(lnd%is:lnd%ie,lnd%js:lnd%je))
allocate(buffert(lnd%is:lnd%ie,lnd%js:lnd%je))

IF (LARGE_DYN_SMALL_STAT) THEN

call read_data('INPUT/river_data.nc', 'connected_to_next', bufferc(:,:), lnd%domain)
call put_to_tiles_r0d_fptr(bufferc, lnd%tile_map, lake_connected_to_next_ptr)

call read_data('INPUT/river_data.nc', 'whole_lake_area', buffer(:,:), lnd%domain)
call put_to_tiles_r0d_fptr(buffer, lnd%tile_map, lake_whole_area_ptr)

call read_data('INPUT/river_data.nc', 'lake_depth_sill', buffer(:,:),  lnd%domain)
buffer = min(buffer, lake_depth_max)
buffer = max(buffer, lake_depth_min)
call put_to_tiles_r0d_fptr(buffer,  lnd%tile_map, lake_depth_sill_ptr)

! lake_tau is just used here as a flag for 'large lakes'
! sill width of -1 is a flag saying not to allow transient storage
call read_data('INPUT/river_data.nc', 'lake_tau', buffert(:,:),  lnd%domain)
buffer = -1.
!where (bufferc.gt.0.5) buffer = lake_width_inside_lake
where (bufferc.lt.0.5 .and. buffert.gt.1.) buffer = large_lake_sill_width
if (lake_specific_width) then
    do i = 1, n_outlet
      if(lnd%face.eq.outlet_face(i).and.lnd%is.le.outlet_i(i).and.lnd%ie.ge.outlet_i(i) &
                                 .and.lnd%js.le.outlet_j(i).and.lnd%je.ge.outlet_j(i)) &
        buffer(outlet_i(i),outlet_j(i)) = outlet_width(i)
      enddo
endif
call put_to_tiles_r0d_fptr(buffer, lnd%tile_map, lake_width_sill_ptr)

buffer = 1.e8
if (max_plain_slope.gt.0.) &
   call read_data('INPUT/river_data.nc', 'max_slope_to_next', buffer(:,:), lnd%domain)
call read_data('INPUT/river_data.nc', 'travel', buffert(:,:), lnd%domain)
bufferc = 0.
where (buffer.lt.max_plain_slope .and. buffert.gt.1.5) bufferc = 1.
call put_to_tiles_r0d_fptr(bufferc, lnd%tile_map, lake_backwater_ptr)
bufferc = 0
where (buffer.lt.max_plain_slope .and. buffert.lt.1.5) bufferc = 1.
call put_to_tiles_r0d_fptr(bufferc, lnd%tile_map, lake_backwater_1_ptr)

ELSE
call read_data('INPUT/river_data.nc', 'whole_lake_area', bufferc(:,:), lnd%domain)

call read_data('INPUT/river_data.nc', 'lake_depth_sill', buffer(:,:), lnd%domain)
where (bufferc.eq.0.)                      buffer = 0.
where (bufferc.gt.0..and.bufferc.lt.2.e10) buffer = max(2., 2.5e-4*sqrt(bufferc))
call put_to_tiles_r0d_fptr(buffer,  lnd%tile_map, lake_depth_sill_ptr)
call put_to_tiles_r0d_fptr(bufferc, lnd%tile_map, lake_whole_area_ptr)

buffer = 4. * buffer
where (bufferc.gt.2.e10) buffer = min(buffer, 60.)
call read_data('INPUT/river_data.nc', 'connected_to_next', bufferc(:,:), lnd%domain)
call put_to_tiles_r0d_fptr(bufferc, lnd%tile_map, lake_connected_to_next_ptr)

where (bufferc.gt.0.5) buffer=lake_width_inside_lake
if (make_all_lakes_wide) buffer = lake_width_inside_lake
call put_to_tiles_r0d_fptr(buffer, lnd%tile_map, lake_width_sill_ptr)
ENDIF

deallocate (buffer, bufferc, buffert)

  ! -------- initialize lake state --------
  te = tail_elmt (lnd%tile_map)
  ce = first_elmt(lnd%tile_map)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)        ! advance position to the next tile
     
     if (.not.associated(tile%lake)) cycle
     
     tile%lake%prog%dz = tile%lake%pars%depth_sill/num_l
     if (init_temp.ge.tfreeze) then
        tile%lake%prog%wl = init_w*tile%lake%prog%dz
        tile%lake%prog%ws = 0
     else
        tile%lake%prog%wl = 0
        tile%lake%prog%ws = init_w*tile%lake%prog%dz
     endif
     tile%lake%prog%T             = init_temp
     tile%lake%prog%groundwater   = init_groundwater
     tile%lake%prog%groundwater_T = init_temp
  enddo

  call get_input_restart_name('INPUT/lake.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('lake_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     if(nfu_inq_var(unit, 'dz')==NF_NOERR) &
          call read_tile_data_r1d_fptr(unit, 'dz', lake_dz_ptr  )
     call read_tile_data_r1d_fptr(unit, 'temp'         , lake_temp_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'           , lake_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'           , lake_ws_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater'  , lake_gw_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater_T', lake_gwT_ptr)
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('lake_init',&
          'cold-starting lake',&
          NOTE)
  endif

  call lake_diag_init ( id_lon, id_lat )
  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_sillw, lnd%tile_map, lake_width_sill_ptr)
  call send_tile_data_r0d_fptr(id_silld, lnd%tile_map, lake_depth_sill_ptr)
  call send_tile_data_r0d_fptr(id_backw, lnd%tile_map, lake_backwater_ptr)
  call send_tile_data_r0d_fptr(id_back1, lnd%tile_map, lake_backwater_1_ptr)
end subroutine lake_init


! ============================================================================
subroutine lake_end ()

  deallocate (zfull, zhalf)
  module_is_initialized =.FALSE.

end subroutine lake_end


! ============================================================================
subroutine save_lake_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit

  call error_mesg('lake_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'lake.res.nc', &
          lnd%coord_glon, lnd%coord_glat, lake_tile_exists, tile_dim_length)

  ! in addition, define vertical coordinate
  if (mpp_pe()==lnd%io_pelist(1)) then
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
     __NF_ASRT__(nfu_put_att(unit,'zfull','positive','down'))
  endif
  call sync_nc_files(unit)
  
  ! write out fields
  call write_tile_data_r1d_fptr(unit,'dz'           ,lake_dz_ptr,  'zfull','layer thickness','m')
  call write_tile_data_r1d_fptr(unit,'temp'         ,lake_temp_ptr,'zfull','lake temperature','degrees_K')
  call write_tile_data_r1d_fptr(unit,'wl'           ,lake_wl_ptr  ,'zfull','liquid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'ws'           ,lake_ws_ptr  ,'zfull','solid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'groundwater'  ,lake_gw_ptr  ,'zfull')
  call write_tile_data_r1d_fptr(unit,'groundwater_T',lake_gwT_ptr ,'zfull')
  
  ! close file
  __NF_ASRT__(nf_close(unit))

end subroutine save_lake_restart


! ============================================================================
subroutine lake_get_sfc_temp(lake, lake_T)
  type(lake_tile_type), intent(in) :: lake
  real, intent(out) :: lake_T

  lake_T = lake%prog(1)%T
end subroutine lake_get_sfc_temp


! ============================================================================
! compute lake-only radiation properties
subroutine lake_radiation ( lake, cosz, &
     lake_refl_dir, lake_refl_dif, lake_refl_lw, lake_emis )
  type(lake_tile_type), intent(in) :: lake
  real, intent(in) :: cosz
  real, intent(out) :: lake_refl_dir(NBANDS), lake_refl_dif(NBANDS), lake_refl_lw, lake_emis

  call lake_data_radiation ( lake, cosz, use_brdf, lake_refl_dir, lake_refl_dif, lake_emis )
  lake_refl_lw = 1 - lake_emis
end subroutine lake_radiation


! ============================================================================
! compute lake-only roughness parameters
subroutine lake_diffusion ( lake, lake_z0s, lake_z0m )
  type(lake_tile_type), intent(in) :: lake
  real, intent(out) :: lake_z0s, lake_z0m

  call lake_data_diffusion ( lake, lake_z0s, lake_z0m )
end subroutine lake_diffusion

! ============================================================================
! update lake properties explicitly for time step.
! integrate lake-heat conduction equation upward from bottom of lake
! to surface, delivering linearization of surface ground heat flux.
subroutine lake_step_1 ( u_star_a, p_surf, latitude, lake, &
                         lake_T, &
                         lake_rh, lake_liq, lake_ice, lake_subl, lake_tf, lake_G0, &
                         lake_DGDT )

  real, intent(in)   :: u_star_a, p_surf, latitude
  type(lake_tile_type), intent(inout) :: lake
  real, intent(out)  :: &
       lake_T, &
       lake_rh, lake_liq, lake_ice, lake_subl, &
       lake_tf, & ! freezing temperature of lake, degK
       lake_G0, &
       lake_DGDT

  ! ---- local vars
  real                  :: bbb, denom, dt_e, tc_dz_eff
  real                  :: z_cum, z_mid, dz_mid, rho_t_mid, k_neutral
  real, dimension(num_l):: aaa, ccc, thermal_cond, heat_capacity, dz_alt, &
                            z_alt, rho_t
  integer               :: l
  real                  :: k_star, N_sq, Ri, u_star, z_liq, z_ice, rho_a
  real                  :: lake_depth, lshc1, lshc2
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  if(is_watch_point()) then
     write(*,*) 'lake_step_1 checkpoint 1'
     write(*,*) 'mask    ', .true.   
     write(*,*) 'T       ', lake_T
     write(*,*) 'rh      ', lake_rh
     write(*,*) 'liq     ', lake_liq
     write(*,*) 'ice     ', lake_ice
     write(*,*) 'subl    ', lake_subl
     write(*,*) 'G0      ', lake_G0
     write(*,*) 'DGDT    ', lake_DGDT
    do l = 1, num_l
      write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
                 ' dz=', lake%prog(l)%dz,&
                 ' T =', lake%prog(l)%T,&
                 ' wl=', lake%prog(l)%wl,&
                 ' ws=', lake%prog(l)%ws, &
                 'K_z=', lake%prog(l)%K_z
      enddo
  endif


  if (relayer_in_step_one) call lake_relayer ( lake )

  lake%prog%K_z = 0.
  lake_T = lake%prog(1)%T
  if (use_rh_feedback) then
      lake_depth = (sum(lake%prog(:)%wl)+sum(lake%prog(:)%ws)) / DENS_H2O
    else
      lake_depth = lake%pars%depth_sill
    endif
  call lake_data_thermodynamics ( lake%pars, lake_depth, lake_rh, &
                                  lake%heat_capacity_dry, thermal_cond )
! Ignore air humidity in converting atmospheric friction velocity to lake value
  rho_a = p_surf/(rdgas*lake_T)
! No momentum transfer through ice cover
  if (lake%prog(1)%ws.le.0. .or. wind_penetrates_ice) then
      u_star = u_star_a*sqrt(rho_a/dens_h2o)
      k_star = 2.79e-5*sqrt(sin(abs(latitude)))*u_star**(-1.84)
      k_star = k_star*(c_drag/1.2e-3)**1.84
    else
      u_star = 0.
      k_star = 1.
    endif
! k_star from B. Henderson-Sellers (1985, Appl. Math. Mod., 9)
!  k_star = 2.79e-5*sqrt(sin(abs(latitude)))*u_star**(-1.84)
  z_cum = 0.
  do l = 1, num_l
    heat_capacity(l) = lake%heat_capacity_dry(l) * lake%prog(l)%dz &
            + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws
    dz_alt(l) = (lake%prog(l)%wl + lake%prog(l)%ws)/dens_h2o
    z_alt(l) = z_cum + 0.5*dz_alt(l)
    z_cum = z_cum + dz_alt(l)
! rho_t from hostetler and bartlein (1990), citing Heggen (1983)
! is a call available in fms?
    rho_t(l) = 1. - 1.9549e-5*abs(lake%prog(l)%T-277.)**1.68
    enddo

  lake_liq  = max(lake%prog(1)%wl, 0.)
  lake_ice  = max(lake%prog(1)%ws, 0.)
  if (lake_ice > 0) then
     lake_subl = 1
  else
     lake_subl = 0
  endif

  if(num_l > 1) then
    if (do_stratify) then
        do l = 1, num_l-1
          if (lake%prog(l)%ws.le.0..and.lake%prog(l+1)%ws.le.0.) then
              dz_mid = z_alt(l+1)-z_alt(l)
              z_mid = 0.5 * (z_alt(l)+z_alt(l+1))
              rho_t_mid = 0.5*(rho_t(l)+rho_t(l+1))
              if (k_star*z_mid .lt. 10.) then
                  k_neutral = vonkarm * u_star * z_mid * exp (-k_star * z_mid)
                else
                  k_neutral = 0.
                endif
              N_sq = (grav / rho_t_mid) * (rho_t(l+1)-rho_t(l)) /dz_mid
              if (N_sq .gt. 0. .and. k_neutral.ne.0.) then
                  ! stability function from B. Henderson-Sellers (1985)
                  Ri = 0.05*(-1. + sqrt(1.+40.*N_sq*(vonkarm*z_mid/u_star)**2 &
                                                   *exp(2.*k_star*z_mid)))
                  lake%prog(l)%K_z = k_neutral / (1. + 37.*Ri*Ri) + K_z_molec
                else if (k_neutral.eq.0.) then
                  lake%prog(l)%K_z = K_z_molec
                else  ! arbitrary constant for unstable mixing
                  lake%prog(l)%K_z = K_z_large
                endif
	      if (lake%pars%depth_sill.gt.2.01) &
	          lake%prog(l)%K_z = K_z_factor &
		   * max(lake%prog(l)%K_z + K_z_background, K_z_min)
              aaa(l+1) = - lake%prog(l)%K_z * delta_time / (dz_alt(l+1)*dz_mid)
              ccc(l)   = - lake%prog(l)%K_z * delta_time / (dz_alt(l  )*dz_mid)
            else
              z_liq = 0.5*(lake%prog(l)%wl+lake%prog(l+1)%wl)/dens_h2o
              z_ice = 0.5*(lake%prog(l)%ws+lake%prog(l+1)%ws)/dens_h2o
              tc_dz_eff = 1. / (z_liq/tc_molec + z_ice/tc_molec_ice)
              aaa(l+1) = - tc_dz_eff * delta_time / heat_capacity(l+1)
              ccc(l)   = - tc_dz_eff * delta_time / heat_capacity(l)
            endif
          enddo
      else
        do l = 1, num_l-1
          tc_dz_eff = 2 / ( dz_alt(l+1)/thermal_cond(l+1) &
             + dz_alt(l)/thermal_cond(l)   )
          aaa(l+1) = - tc_dz_eff * delta_time / heat_capacity(l+1)
          ccc(l)   = - tc_dz_eff * delta_time / heat_capacity(l)
          enddo
      endif

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(lake%prog(num_l)%T - lake%prog(num_l-1)%T) &
               + lake%geothermal_heat_flux * delta_time / heat_capacity(num_l)
     lake%e(num_l-1) = -aaa(num_l)/denom
     lake%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*lake%e(l)
        dt_e = - ( ccc(l)*(lake%prog(l+1)%T - lake%prog(l)%T  ) &
                  -aaa(l)*(lake%prog(l)%T   - lake%prog(l-1)%T) )
        lake%e(l-1) = -aaa(l)/denom
        lake%f(l-1) = (dt_e - ccc(l)*lake%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     lake_G0    = ccc(1)*(lake%prog(2)%T- lake%prog(1)%T &
          + lake%f(1)) / denom
     lake_DGDT  = (1 - ccc(1)*(1-lake%e(1))) / denom   
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     lake_G0    = 0.
     lake_DGDT  = 1. / denom
  end if
  
  ! set the freezing temperature of the lake
  lake_tf = tfreeze
  
  if(is_watch_point()) then
     write(*,*) 'lake_step_1 checkpoint 2'
     write(*,*) 'mask    ', .true.   
     write(*,*) 'T       ', lake_T
     write(*,*) 'rh      ', lake_rh
     write(*,*) 'liq     ', lake_liq
     write(*,*) 'ice     ', lake_ice
     write(*,*) 'subl    ', lake_subl
     write(*,*) 'G0      ', lake_G0
     write(*,*) 'DGDT    ', lake_DGDT
    do l = 1, num_l
      write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
                 ' dz=', lake%prog(l)%dz,&
                 ' T =', lake%prog(l)%T,&
                 ' wl=', lake%prog(l)%wl,&
                 ' ws=', lake%prog(l)%ws, &
                 'K_z=', lake%prog(l)%K_z
      enddo
  endif

end subroutine lake_step_1


! ============================================================================
! apply boundary flows to lake water and move lake water vertically.
  subroutine lake_step_2 ( lake, diag, lake_subl, snow_lprec, snow_hlprec,  &
                           subs_DT, subs_M_imp, subs_evap, &
                           use_tfreeze_in_grnd_latent, &
                           lake_levap, lake_fevap, lake_melt, &
                           lake_Ttop, lake_Ctop )
  type(lake_tile_type), intent(inout) :: lake
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
     lake_subl     !
  real, intent(in) :: &
     snow_lprec, &
     snow_hlprec, &
     subs_DT,       &!
     subs_M_imp,       &! rate of phase change of non-evaporated lake water
     subs_evap
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
     lake_levap, lake_fevap, lake_melt, &
     lake_Ttop, lake_Ctop

  ! ---- local vars
  real, dimension(num_l) :: del_t, eee, fff, &
             psi, DThDP, hyd_cond, DKDP, K, DKDPm, DKDPp, grad, &
             dW_l, u_minus, u_plus, DPsi, lake_w_fc
  real, dimension(num_l+1) :: flow
  real, dimension(num_l  ) :: div
  real :: ice_to_move, h_upper, h_lower, h_to_move_up, &
     lprec_eff, hlprec_eff, hcap, dheat, &
     melt_per_deg, melt, lshc1, lshc2
  real, dimension(num_l-1) :: del_z
  real :: jj
  integer :: l

  jj = 1.
  
  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 1 ***** '
    write(*,*) 'mask    ', .true.   
    write(*,*) 'subs_evap    ', subs_evap   
    write(*,*) 'snow_lprec   ', snow_lprec  
    write(*,*) 'subs_M_imp   ', subs_M_imp   
    write(*,*) 'theta_s ', lake%pars%w_sat
    do l = 1, num_l
      write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
                 ' dz=', lake%prog(l)%dz,&
                 ' T =', lake%prog(l)%T,&
                 ' Th=', (lake%prog(l)%ws &
                         +lake%prog(l)%wl)/(dens_h2o*lake%prog(l)%dz),&
                 ' wl=', lake%prog(l)%wl,&
                 ' ws=', lake%prog(l)%ws,&
                 ' gw=', lake%prog(l)%groundwater
      enddo
  endif

  ! ---- record fluxes ---------
  lake_levap  = subs_evap*(1-lake_subl)
  lake_fevap  = subs_evap*   lake_subl
  lake_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  lake%prog(1)%T = lake%prog(1)%T + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = lake%e(l) * del_t(l) + lake%f(l)
      lake%prog(l+1)%T = lake%prog(l+1)%T + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 2 ***** '
    do l = 1, num_l
       write(*,*) 'level=', l, 'T', lake%prog(l)%T
    enddo
  endif

  ! ---- extract evap from lake and do implicit melt --------------------
  lake%prog(1)%wl = lake%prog(1)%wl - lake_levap*delta_time
  lake%prog(1)%ws = lake%prog(1)%ws - lake_fevap*delta_time
  hcap = lake%heat_capacity_dry(1)*lake%prog(1)%dz &
                     + clw*lake%prog(1)%wl + csw*lake%prog(1)%ws
  ! T adjustment for nonlinear terms (del_T)*(del_W)
  dheat = delta_time*(clw*lake_levap+csw*lake_fevap)*del_T(1)
  ! take out extra heat not claimed in advance for evaporation
  if (use_tfreeze_in_grnd_latent) dheat = dheat &
          - delta_time*((cpw-clw)*lake_levap+(cpw-csw)*lake_fevap) &
                                 *(lake%prog(1)%T-del_T(1)-tfreeze)
  lake%prog(1)%T  = lake%prog(1)%T  + dheat/hcap
  lake%prog(1)%wl = lake%prog(1)%wl + subs_M_imp
  lake%prog(1)%ws = lake%prog(1)%ws - subs_M_imp
  lake%prog(1)%T  = tfreeze + (hcap*(lake%prog(1)%T-tfreeze) ) &
                            / ( hcap + (clw-csw)*subs_M_imp )

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 2.1 ***** '
     do l = 1, num_l
        write(*,*) 'level=', l, 'T', lake%prog(l)%T
     enddo
  endif

  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
  ! note that only liquid inputs are received by  lake from snow pack. any
  ! snow fall just creates a snow pack on top  of lake, even if lake is not
  ! frozen. but snow pack on top of unfrozen lake will interact thermally,
  ! so that either lake freezes or snow melts and falls in.
    flow=1
    flow(1)  = snow_lprec *delta_time
    do l = 1, num_l
      flow(l+1) = 0
      dW_l(l) = flow(l) - flow(l+1)
      lake%prog(l)%wl = lake%prog(l)%wl + dW_l(l)
    enddo

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 3.3 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' wl=', lake%prog(l)%wl,&
             'flow=', flow(l)
     enddo
  endif
  
  hcap = lake%heat_capacity_dry(1)*lake%prog(1)%dz &
                     + clw*(lake%prog(1)%wl-dW_l(1)) + csw*lake%prog(1)%ws
  lake%prog(1)%T = tfreeze + (hcap*(lake%prog(1)%T-tfreeze) +  &
                                 snow_hlprec*delta_time) &
                            / ( hcap + clw*dW_l(1) )

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 3.4 ***** '
    write(*,*) ' tfreeze', tfreeze
    write(*,*) ' snow_hlprec', snow_hlprec
  endif

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = lake%heat_capacity_dry(l)*lake%prog(l)%dz &
             + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws
    melt_per_deg = hcap/hlf
    if (lake%prog(l)%ws>0 .and. lake%prog(l)%T>tfreeze) then
      melt =  min(lake%prog(l)%ws, (lake%prog(l)%T-tfreeze)*melt_per_deg)
    else if (lake%prog(l)%wl>0 .and. lake%prog(l)%T<tfreeze) then
      melt = -min(lake%prog(l)%wl, (tfreeze-lake%prog(l)%T)*melt_per_deg)
    else
      melt = 0
    endif
    lake%prog(l)%wl = lake%prog(l)%wl + melt
    lake%prog(l)%ws = lake%prog(l)%ws - melt
    lake%prog(l)%T = tfreeze &
       + (hcap*(lake%prog(l)%T-tfreeze) - hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
    lake_melt = lake_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 5 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' dz=', lake%prog(l)%dz,&
             ' T =', lake%prog(l)%T,&
             ' Th=', (lake%prog(l)%ws +lake%prog(l)%wl)/(dens_h2o*lake%prog(l)%dz),&
             ' wl=', lake%prog(l)%wl,&
             ' ws=', lake%prog(l)%ws,&
             ' gw=', lake%prog(l)%groundwater
     enddo
  endif

  if (.not.relayer_in_step_one) call lake_relayer ( lake )

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 6 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' dz=', lake%prog(l)%dz,&
             ' T =', lake%prog(l)%T,&
             ' Th=', (lake%prog(l)%ws +lake%prog(l)%wl)/(dens_h2o*lake%prog(l)%dz),&
             ' wl=', lake%prog(l)%wl,&
             ' ws=', lake%prog(l)%ws,&
             ' gw=', lake%prog(l)%groundwater
     enddo
  endif

  if (float_ice_to_top) then
      do l = num_l, 2, -1
        if (lake%prog(l)%ws .gt. 0. .and. lake%prog(l-1)%wl .gt. 0.) then
            ice_to_move = min(lake%prog(l)%ws, lake%prog(l-1)%wl)
            h_upper = (clw*lake%prog(l-1)%wl+csw*lake%prog(l-1)%ws)*lake%prog(l-1)%T
            h_lower = (clw*lake%prog(l  )%wl+csw*lake%prog(l  )%ws)*lake%prog(l  )%T
            lake%prog(l-1)%wl = lake%prog(l-1)%wl - ice_to_move
            lake%prog(l-1)%ws = lake%prog(l-1)%ws + ice_to_move
            lake%prog(l  )%wl = lake%prog(l  )%wl + ice_to_move
            lake%prog(l  )%ws = lake%prog(l  )%ws - ice_to_move
            h_to_move_up = ice_to_move*(csw*lake%prog(l)%T-clw*lake%prog(l-1)%T)
            h_upper  = h_upper + h_to_move_up
            h_lower  = h_lower - h_to_move_up
            lake%prog(l-1)%T = h_upper / (clw*lake%prog(l-1)%wl+csw*lake%prog(l-1)%ws)
            lake%prog(l  )%T = h_lower / (clw*lake%prog(l  )%wl+csw*lake%prog(l  )%ws)
          endif
        enddo
    endif

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 7 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' dz=', lake%prog(l)%dz,&
             ' T =', lake%prog(l)%T,&
             ' Th=', (lake%prog(l)%ws +lake%prog(l)%wl)/(dens_h2o*lake%prog(l)%dz),&
             ' wl=', lake%prog(l)%wl,&
             ' ws=', lake%prog(l)%ws,&
             ' gw=', lake%prog(l)%groundwater
     enddo
  endif


  lake_Ttop = lake%prog(1)%T
  lake_Ctop = lake%heat_capacity_dry(1)*lake%prog(1)%dz &
       + clw*lake%prog(1)%wl + csw*lake%prog(1)%ws

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!  

  ! ---- increment time
  time = increment_time(time, int(delta_time), 0)

  ! ---- diagnostic section
  call send_tile_data (id_dz,   lake%prog%dz,     diag )
  call send_tile_data (id_temp, lake%prog%T,     diag )
  call send_tile_data (id_wl,  lake%prog(1:num_l)%wl, diag )
  call send_tile_data (id_ws,  lake%prog(1:num_l)%ws, diag )
  call send_tile_data (id_lwc,  lake%prog(1:num_l)%wl/lake%prog(1:num_l)%dz, diag )
  call send_tile_data (id_swc,  lake%prog(1:num_l)%ws/lake%prog(1:num_l)%dz, diag )
  call send_tile_data (id_K_z,  lake%prog(1:num_l)%K_z,        diag )
  call send_tile_data (id_evap, lake_levap+lake_fevap, diag )

end subroutine lake_step_2


! ============================================================================
!
  subroutine lake_relayer ( lake )
  type(lake_tile_type), intent(inout) :: lake

  ! ---- local vars
  integer :: l, l_lowest_thin_layer, l_highest_thick_layer
  real :: new_dz, new_ws, new_wl, new_h, new_T, liq_frac

! now check whether we need to re-layer the lake.
  if ( (lake%prog(1)%wl+lake%prog(1)%ws) &
      /(lake%prog(2)%wl+lake%prog(2)%ws) .gt. max_rat) then
      ! top layer has grown too thick. join two lower layers, and split
      ! top layer into two layers. in special case, just join and
      ! re-split top two layers.
      l_lowest_thin_layer = num_l
      do l = 2, num_l-1
        if (lake%prog(l)%dz.lt.0.99*lake%prog(num_l)%dz) l_lowest_thin_layer = l
        enddo
      if (l_lowest_thin_layer.gt.2) then
          new_dz = lake%prog(l_lowest_thin_layer)%dz &
                  + lake%prog(l_lowest_thin_layer-1)%dz
          new_wl = lake%prog(l_lowest_thin_layer)%wl &
                  + lake%prog(l_lowest_thin_layer-1)%wl
          new_ws = lake%prog(l_lowest_thin_layer)%ws &
                  + lake%prog(l_lowest_thin_layer-1)%ws
          new_h  = ( clw*lake%prog(l_lowest_thin_layer)%wl &
                   + csw*lake%prog(l_lowest_thin_layer)%ws) &
                   *     lake%prog(l_lowest_thin_layer)%T &
                 + ( clw*lake%prog(l_lowest_thin_layer-1)%wl &
                   + csw*lake%prog(l_lowest_thin_layer-1)%ws) &
                   *     lake%prog(l_lowest_thin_layer-1)%T
          new_T = new_h / (clw*new_wl+csw*new_ws)
          lake%prog(l_lowest_thin_layer)%dz = new_dz
          lake%prog(l_lowest_thin_layer)%wl = new_wl
          lake%prog(l_lowest_thin_layer)%ws = new_ws
          lake%prog(l_lowest_thin_layer)%T  = new_T
          do l = l_lowest_thin_layer-1, 3, -1
            lake%prog(l)%dz = lake%prog(l-1)%dz
            lake%prog(l)%wl = lake%prog(l-1)%wl
            lake%prog(l)%ws = lake%prog(l-1)%ws
            lake%prog(l)%T  = lake%prog(l-1)%T
            enddo
          liq_frac = lake%prog(1)%wl / (lake%prog(1)%wl+lake%prog(1)%ws)
          lake%prog(2)%wl =     liq_frac *DENS_H2O*lake%prog(2)%dz
          lake%prog(2)%ws = (1.-liq_frac)*DENS_H2O*lake%prog(2)%dz
          lake%prog(2)%T  = lake%prog(1)%T
          lake%prog(1)%dz = lake%prog(2)%dz
          lake%prog(1)%wl = lake%prog(1)%wl - lake%prog(2)%wl
          lake%prog(1)%ws = lake%prog(1)%ws - lake%prog(2)%ws
        else
          new_wl = lake%prog(1)%wl + lake%prog(2)%wl
          new_ws = lake%prog(1)%ws + lake%prog(2)%ws
          new_h  = ( clw*lake%prog(1)%wl + csw*lake%prog(1)%ws)   &
                                                 * lake%prog(1)%T &
                 + ( clw*lake%prog(2)%wl + csw*lake%prog(2)%ws)    &
                                                 * lake%prog(2)%T
          new_T  = new_h / (clw*new_wl+csw*new_ws)
          liq_frac = new_wl / (new_wl+new_ws)
          lake%prog(2)%dz = lake%prog(3)%dz
          lake%prog(2)%wl =     liq_frac *DENS_H2O*lake%prog(2)%dz
          lake%prog(2)%ws = (1.-liq_frac)*DENS_H2O*lake%prog(2)%dz
          lake%prog(2)%T  = new_T
          lake%prog(1)%dz = lake%prog(2)%dz
          lake%prog(1)%wl = new_wl - lake%prog(2)%wl
          lake%prog(1)%ws = new_ws - lake%prog(2)%ws
          lake%prog(1)%T  = new_T
        endif
    else if(  (lake%prog(1)%wl+lake%prog(1)%ws) &
             /(lake%prog(2)%wl+lake%prog(2)%ws) .lt. min_rat) then
      ! top layer has grown too thin. join with next layer down, and split
      ! a lower layer to maintain number of layers.  in special case, just
      ! join and re-split top two layers.
      l_highest_thick_layer = 2
      do l = num_l, 3, -1
        if (lake%prog(l)%dz.gt.1.01*lake%prog(2)%dz) l_highest_thick_layer = l
        enddo
      new_wl = lake%prog(1)%wl + lake%prog(2)%wl
      new_ws = lake%prog(1)%ws + lake%prog(2)%ws
      new_h  = ( clw*lake%prog(1)%wl + csw*lake%prog(1)%ws)   &
                                             * lake%prog(1)%T &
             + ( clw*lake%prog(2)%wl + csw*lake%prog(2)%ws)    &
                                             * lake%prog(2)%T
      new_T  = new_h / (clw*new_wl+csw*new_ws)
      if (l_highest_thick_layer.gt.2) then
          lake%prog(1)%dz = lake%prog(2)%dz
          lake%prog(1)%wl = new_wl
          lake%prog(1)%ws = new_ws
          lake%prog(1)%T  = new_T
          do l = 2, l_highest_thick_layer-2
            lake%prog(l)%dz = lake%prog(l+1)%dz
            lake%prog(l)%wl = lake%prog(l+1)%wl
            lake%prog(l)%ws = lake%prog(l+1)%ws
            lake%prog(l)%T  = lake%prog(l+1)%T
            enddo
          new_dz = lake%prog(l_highest_thick_layer)%dz / 2.
          new_wl = lake%prog(l_highest_thick_layer)%wl / 2.
          new_ws = lake%prog(l_highest_thick_layer)%ws / 2.
          new_T  = lake%prog(l_highest_thick_layer)%T
          lake%prog(l_highest_thick_layer-1)%dz = new_dz
          lake%prog(l_highest_thick_layer-1)%wl = new_wl
          lake%prog(l_highest_thick_layer-1)%ws = new_ws
          lake%prog(l_highest_thick_layer-1)%T  = new_T
          lake%prog(l_highest_thick_layer)%dz = new_dz
          lake%prog(l_highest_thick_layer)%wl = new_wl
          lake%prog(l_highest_thick_layer)%ws = new_ws
          lake%prog(l_highest_thick_layer)%T  = new_T
        else
          liq_frac = new_wl / (new_wl+new_ws)
          lake%prog(2)%dz = lake%prog(3)%dz / 2.
          lake%prog(2)%wl =     liq_frac *DENS_H2O*lake%prog(2)%dz
          lake%prog(2)%ws = (1.-liq_frac)*DENS_H2O*lake%prog(2)%dz
          lake%prog(2)%T  = new_T
          lake%prog(1)%dz = lake%prog(2)%dz
          lake%prog(1)%wl = new_wl - lake%prog(2)%wl
          lake%prog(1)%ws = new_ws - lake%prog(2)%ws
          lake%prog(1)%T  = new_T          
        endif
    endif
  end subroutine lake_relayer

! ============================================================================
subroutine lake_diag_init ( id_lon, id_lat )
  integer,         intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer,         intent(in) :: id_lat  ! ID of land longitude (X) axis

  ! ---- local vars
  integer :: axes(3)
  integer :: id_zhalf, id_zfull

  ! define vertical axis
  id_zhalf = diag_axis_init ( &
       'zhalf_lake', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='lake' )
  id_zfull = diag_axis_init ( &
       'zfull_lake', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='lake', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define static diagnostic fields
  id_sillw = register_tiled_static_field ( module_name, 'lake_width', &
       axes(1:2), 'lake width at outflow', 'm', missing_value=-100.0 )
  id_silld = register_tiled_static_field ( module_name, 'lake_depth', &
       axes(1:2), 'lake depth below sill', 'm', missing_value=-100.0 )
  id_backw = register_tiled_static_field ( module_name, 'backwater', &
       axes(1:2), 'backwater flag', '-', missing_value=-100.0 )
  id_back1 = register_tiled_static_field ( module_name, 'backwater_1', &
       axes(1:2), 'backwater1 flag', '-', missing_value=-100.0 )
  ! define dynamic diagnostic fields
  id_dz  = register_tiled_diag_field ( module_name, 'lake_dz', axes,         &
       Time, 'nominal layer thickness', 'm', missing_value=-100.0 )
  id_wl  = register_tiled_diag_field ( module_name, 'lake_wl', axes,         &
       Time, 'liquid water mass', 'kg/m2', missing_value=-100.0 )
  id_ws  = register_tiled_diag_field ( module_name, 'lake_ws', axes,         &
       Time, 'solid water mass', 'kg/m2', missing_value=-100.0 )
  id_lwc  = register_tiled_diag_field ( module_name, 'lake_liq',  axes,       &
       Time, 'bulk density of liquid water', 'kg/m3',  missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'lake_ice',  axes,       &
       Time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'lake_T',  axes,        &
       Time, 'temperature',            'degK',  missing_value=-100.0 )
  id_K_z  = register_tiled_diag_field ( module_name, 'lake_K_z', axes,         &
       Time, 'vertical diffusivity', 'm2/s', missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'lake_evap',  axes(1:2),  &
       Time, 'lake evap',            'kg/(m2 s)',  missing_value=-100.0 )
       
  call add_tiled_static_field_alias (id_silld, module_name, 'sill_depth', &
       axes(1:2), 'obsolete, pls use lake_depth (static)','m', &
       missing_value=-100.0 )
  call add_tiled_static_field_alias (id_sillw, module_name, 'sill_width', &
       axes(1:2), 'obsolete, pls use lake_width (static)','m', &
       missing_value=-100.0 )

end subroutine lake_diag_init

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function lake_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   lake_tile_exists = associated(tile%lake)
end function lake_tile_exists


! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine lake_dz_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   integer :: n
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) then
        n = size(tile%lake%prog)
        ptr(1:n) => tile%lake%prog(1:n)%dz
      endif
   endif
end subroutine lake_dz_ptr

subroutine lake_temp_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   integer :: n
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) then
        n = size(tile%lake%prog)
        ptr(1:n) => tile%lake%prog(1:n)%T
      endif
   endif
end subroutine lake_temp_ptr

subroutine lake_wl_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   integer :: n
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) then
        n = size(tile%lake%prog)
        ptr(1:n) => tile%lake%prog(1:n)%wl
      endif
   endif
end subroutine lake_wl_ptr

subroutine lake_ws_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   integer :: n
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) then
        n = size(tile%lake%prog)
        ptr(1:n) => tile%lake%prog(1:n)%ws
      endif
   endif
end subroutine lake_ws_ptr

subroutine lake_gw_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   integer :: n
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) then
        n = size(tile%lake%prog)
        ptr(1:n) => tile%lake%prog(1:n)%groundwater
      endif
   endif
end subroutine lake_gw_ptr

subroutine lake_gwT_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   integer :: n
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) then
        n = size(tile%lake%prog)
        ptr(1:n) => tile%lake%prog(1:n)%groundwater_T
      endif
   endif
end subroutine lake_gwT_ptr

subroutine lake_connected_to_next_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%connected_to_next
   endif
end subroutine lake_connected_to_next_ptr

subroutine lake_depth_sill_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%depth_sill
   endif
end subroutine lake_depth_sill_ptr

subroutine lake_whole_area_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%whole_area
   endif
end subroutine lake_whole_area_ptr

subroutine lake_width_sill_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%width_sill
   endif
end subroutine lake_width_sill_ptr

subroutine lake_backwater_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%backwater
   endif
end subroutine lake_backwater_ptr

subroutine lake_backwater_1_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%backwater_1
   endif
end subroutine lake_backwater_1_ptr

end module lake_mod



