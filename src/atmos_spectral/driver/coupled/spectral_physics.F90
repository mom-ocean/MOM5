module spectral_physics_mod

use fms_mod,               only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number, set_domain, file_exist, &
                                 field_size, read_data, write_data, fms_init

use constants_mod,         only: grav, pi

use time_manager_mod,      only: time_type, set_time, get_time, operator(-), operator(/=), time_manager_init

use press_and_geopot_mod,  only: pressure_variables, compute_pressures_and_heights

use transforms_mod,        only: get_grid_boundaries, get_deg_lon, get_deg_lat, get_wts_lat, &
                                 get_grid_domain, get_lon_max, get_lat_max, grid_domain, area_weighted_global_mean

use global_integral_mod,   only: mass_weighted_global_integral

use spectral_dynamics_mod, only: get_reference_sea_level_press, get_num_levels

use physics_driver_mod,    only: physics_driver_init, physics_driver_down, physics_driver_up, physics_driver_end, &
                                 surf_diff_type, do_moist_in_phys_up, get_diff_t, get_radturbten, zero_radturbten

use moist_processes_mod,   only: moist_processes_init, moist_processes, moist_processes_end

use mcm_moist_processes_mod, only: mcm_moist_processes, mcm_moist_processes_init, mcm_moist_processes_end

use tracer_type_mod,       only: tracer_type

use  field_manager_mod,    only: MODEL_ATMOS

use tracer_manager_mod,    only: get_number_tracers, get_tracer_index, NO_TRACER

implicit none
private

public :: spectral_physics_init, spectral_physics_down, spectral_physics_up, &
          spectral_physics_end, surf_diff_type, spectral_physics_moist

character(len=128), parameter :: version = &
'$Id: spectral_physics.F90,v 13.0 2006/03/28 21:17:25 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: tikal $'

integer, parameter :: num_time_levels=2

real, allocatable, dimension(:,:      ) :: rad_lon_2d, rad_lat_2d, area_2d
real, allocatable, dimension(:,:,:    ) :: diff_cu_mo
logical, allocatable, dimension(:,:   ) :: convect
real, allocatable, dimension(:,:,:,:  ) :: diag_tracers
integer :: num_levels, num_tracers, nhum
integer :: is, ie, js, je
logical :: module_is_initialized = .false.
logical :: do_mcm_moist_processes

contains

!------------------------------------------------------------------------------------------------

subroutine spectral_physics_init(Time, axes, Surf_diff, nhum_in, p_half, do_mcm_moist_processes_in)

type(time_type), intent(in) :: Time
integer, intent(in),    dimension(:) :: axes
type(surf_diff_type), intent(inout) :: Surf_diff
integer, intent(in) :: nhum_in
real, intent(in), dimension(:,:,:) :: p_half
logical, intent(in) :: do_mcm_moist_processes_in
real, allocatable, dimension(:,:,:,:) :: grid_tracers

real, allocatable, dimension(:) :: rad_lon, rad_lat, wts_lat, lon_boundaries, lat_boundaries
real, dimension(2) :: radiation_ref_press_surf = (/ 101325., 81060. /)

real, allocatable, dimension(:,:) :: radiation_ref_press, rconvect
real, allocatable, dimension(:)   :: p_half_1d, ln_p_half_1d, shalf
real, allocatable, dimension(:)   :: p_full_1d, ln_p_full_1d, sfull

real :: reference_sea_level_press
integer :: i, j, num_diag, lon_max, lat_max
character(len=4) :: ch1, ch2, ch3, ch4, ch5, ch6
character(len=64) :: file
integer, dimension(4) :: siz
logical :: do_donner

if(module_is_initialized) return

call write_version_number(version, tagname)

call fms_init
call time_manager_init

nhum = nhum_in
do_mcm_moist_processes = do_mcm_moist_processes_in

call get_grid_domain(is, ie, js, je)
allocate(rad_lon(is:ie), rad_lat(js:je), wts_lat(js:je))
allocate(lon_boundaries(ie-is+2), lat_boundaries(je-js+2))

call get_num_levels(num_levels)
allocate(radiation_ref_press(num_levels+1,2))
allocate(p_half_1d(num_levels+1), ln_p_half_1d(num_levels+1), shalf(num_levels+1))
allocate(p_full_1d(num_levels  ), ln_p_full_1d(num_levels  ), sfull(num_levels  ))

allocate (rad_lon_2d(is:ie,js:je))
allocate (rad_lat_2d(is:ie,js:je))
allocate (   area_2d(is:ie,js:je))

call pressure_variables(p_half_1d,ln_p_half_1d,radiation_ref_press(1:num_levels,1),ln_p_full_1d,radiation_ref_press_surf(1))
call pressure_variables(p_half_1d,ln_p_half_1d,radiation_ref_press(1:num_levels,2),ln_p_full_1d,radiation_ref_press_surf(2))
radiation_ref_press(num_levels+1,:) = radiation_ref_press_surf

call get_reference_sea_level_press(reference_sea_level_press)
call pressure_variables(p_half_1d, ln_p_half_1d, p_full_1d, ln_p_full_1d, reference_sea_level_press)
shalf = p_half_1d/reference_sea_level_press
sfull = p_full_1d/reference_sea_level_press

call get_lon_max(lon_max)
call get_lat_max(lat_max)
call get_deg_lon(rad_lon)
call get_deg_lat(rad_lat)
call get_wts_lat(wts_lat)
rad_lon = pi*rad_lon/180.
rad_lat = pi*rad_lat/180.
do j=js,je
  rad_lat_2d(:,j) = rad_lat(j)
  area_2d(:,j)    = wts_lat(j)/(2.*lon_max)
enddo
do i=is,ie
  rad_lon_2d(i,:) = rad_lon(i)
enddo

call get_grid_boundaries(lon_boundaries, lat_boundaries)

call get_number_tracers(MODEL_ATMOS, num_diag=num_diag)
if(num_diag > 0) then
  call error_mesg('spectral_physics_init', &
   'This version of the spectral atmospheric model not coded to handle diagnostic tracers',FATAL)
endif
allocate (diag_tracers(is:ie,js:je,num_levels,num_diag))
diag_tracers = 0.

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)

call set_domain(grid_domain)

if(do_mcm_moist_processes) then
  call mcm_moist_processes_init( ie-is+1, je-js+1, num_levels, axes, Time)
else
  call moist_processes_init(ie-is+1, je-js+1, num_levels, lon_boundaries, lat_boundaries, &
                            radiation_ref_press(:,1), axes, Time, do_donner)
endif

allocate(grid_tracers(is:ie, js:je, num_levels, num_tracers))
grid_tracers = 0.

call physics_driver_init(Time, lon_boundaries, lat_boundaries, axes, radiation_ref_press, grid_tracers, Surf_diff, p_half)

if(sum(grid_tracers) /= 0.) then
  call error_mesg('spectral_physics_init','This version of the spectral atmospheric model not coded to handle'// &
                  ' initialization of tracer fields by physics_driver_init',FATAL)
endif

deallocate(rad_lon, rad_lat, wts_lat, lon_boundaries, lat_boundaries, grid_tracers)
deallocate(radiation_ref_press, p_half_1d, ln_p_half_1d, shalf, p_full_1d, ln_p_full_1d, sfull)

allocate (diff_cu_mo(is:ie,js:je,num_levels))
allocate ( convect(is:ie,js:je))
file = 'INPUT/spectral_physics.res.nc'
if(file_exist(trim(file))) then
  call field_size(trim(file), 'diff_cu_mo', siz)
  if(siz(1) /= lon_max .or. siz(2) /= lat_max .or. siz(3) /= num_levels) then
    write(ch1,'(i4)') siz(1)
    write(ch2,'(i4)') siz(2)
    write(ch3,'(i4)') siz(3)
    write(ch4,'(i4)') lon_max
    write(ch5,'(i4)') lat_max
    write(ch6,'(i4)') num_levels
    call error_mesg('spectral_physics_init','Resolution of restart data is incorrect.'// &
    ' Restart data: lon_max='//ch1//', lat_max='//ch2//', num_levels='//ch3// &
    '    Should be: lon_max='//ch4//', lat_max='//ch5//', num_levels='//ch6, FATAL)
  endif
  call read_data(trim(file), 'diff_cu_mo', diff_cu_mo, grid_domain)
  allocate (rconvect(is:ie,js:je))
  call read_data(trim(file), 'convect',    rconvect,   grid_domain) ! No interface for reading/writing netcdf logicals
  where(rconvect == 1.)
    convect = .true.
  elsewhere
    convect = .false.
  endwhere
  deallocate (rconvect)
else
  diff_cu_mo = 0.
  convect = .false.
endif

module_is_initialized = .true.

return
end subroutine spectral_physics_init
!------------------------------------------------------------------------------------------------

subroutine spectral_physics_down(Time_prev, Time, Time_next, previous, current,                                   &
                        p_half, p_full, z_half, z_full, psg, ug, vg, tg, grid_tracers,                            &
                        frac_land, rough_mom, albedo, t_surf, u_star, b_star, q_star, dtau_du, dtau_dv, tau_x, tau_y, &
                        albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif,                           &
                        dt_ug, dt_vg, dt_tg, dt_tracers, flux_sw, flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir, &
                        flux_sw_down_vis_dif, flux_sw_down_total_dir, flux_sw_down_total_dif, flux_sw_vis,        &
                        flux_sw_vis_dir, flux_sw_vis_dif, flux_lw, coszen, gust, Surf_diff)

type(time_type), intent(in) :: Time_prev, Time, Time_next
integer, intent(in)         :: previous, current
real,    intent(in),    dimension(:,:,:    ) :: p_full, z_full
real,    intent(in),    dimension(:,:,:    ) :: p_half, z_half
real,    intent(in),    dimension(:,:,:    ) :: psg
real,    intent(in),    dimension(:,:,:,:  ) :: ug, vg, tg
real,    intent(inout), dimension(:,:,:,:,:) :: grid_tracers
real,    intent(in),    dimension(:,:      ) :: frac_land, rough_mom, albedo, t_surf, u_star, b_star, q_star, dtau_du, dtau_dv
real,    intent(inout), dimension(:,:      ) :: tau_x, tau_y
real,    intent(in),    dimension(:,:      ) :: albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif
real,    intent(inout), dimension(:,:,:    ) :: dt_ug, dt_vg, dt_tg
real,    intent(inout), dimension(:,:,:,:  ) :: dt_tracers
real,    intent(out),   dimension(:,:      ) :: flux_sw, flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir
real,    intent(out),   dimension(:,:      ) :: flux_sw_down_vis_dif, flux_sw_down_total_dir, flux_sw_down_total_dif
real,    intent(out),   dimension(:,:      ) :: flux_sw_vis, flux_sw_vis_dir, flux_sw_vis_dif, flux_lw, coszen, gust
type(surf_diff_type), intent(inout) :: Surf_diff

real, dimension(size(dt_tracers,4)) :: gavg_rrv
integer :: nco2

!**************************************************************************************

if(.not.module_is_initialized) then
  call error_mesg('spectral_physics_down','spectral_physics module is not initialized', FATAL)
end if

nco2 = get_tracer_index(MODEL_ATMOS, 'co2')
gavg_rrv = 0.0
if(nco2 /= NO_TRACER) then
  gavg_rrv(nco2) = mass_weighted_global_integral(grid_tracers(:,:,:,current,nco2), psg(:,:,current))/ &
                    (area_weighted_global_mean(psg(:,:,current))/grav)
endif

if(do_moist_in_phys_up()) then
  call physics_driver_down(1, ie-is+1, 1, je-js+1, Time_prev, Time, Time_next,   &
            rad_lat_2d,    rad_lon_2d,  area_2d,                                 &
                p_half,      p_full, z_half, z_full,                             &
                    ug(:,:,:,current), vg(:,:,:,current),                        &
                    tg(:,:,:,current), grid_tracers(:,:,:,current,nhum),         &
          grid_tracers(:,:,:,current,:), ug(:,:,:,previous), vg(:,:,:,previous), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,nhum),       &
          grid_tracers(:,:,:,previous,:),                                        &
             frac_land, rough_mom,      albedo,                                  &
               albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif,   &
                t_surf,    u_star,      b_star, q_star,                          &
               dtau_du, dtau_dv,     tau_x,       tau_y,                         &
               dt_ug,       dt_vg,         dt_tg,                                &
               dt_tracers(:,:,:,nhum), dt_tracers, flux_sw(:,:),                 &
               flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir,                   &
               flux_sw_down_vis_dif, flux_sw_down_total_dir,                     &
               flux_sw_down_total_dif, flux_sw_vis,                              &
               flux_sw_vis_dir, flux_sw_vis_dif,                                 &
               flux_lw, coszen, gust, Surf_diff, gavg_rrv)
else
  call physics_driver_down(1, ie-is+1, 1, je-js+1, Time_prev, Time, Time_next,   &
            rad_lat_2d,    rad_lon_2d,  area_2d,                                 &
                p_half,      p_full, z_half, z_full,                             &
                    ug(:,:,:,current), vg(:,:,:,current),                        &
                    tg(:,:,:,current), grid_tracers(:,:,:,current,nhum),         &
          grid_tracers(:,:,:,current,:), ug(:,:,:,previous), vg(:,:,:,previous), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,nhum),       &
          grid_tracers(:,:,:,previous,:),                                        &
             frac_land, rough_mom,      albedo,                                  &
               albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif,   &
                t_surf,    u_star,      b_star, q_star,                          &
               dtau_du, dtau_dv,     tau_x,       tau_y,                         &
                 dt_ug,       dt_vg,         dt_tg,                              &
                 dt_tracers(:,:,:,nhum), dt_tracers, flux_sw,                    &
               flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir,                   &
               flux_sw_down_vis_dif, flux_sw_down_total_dir,                     &
               flux_sw_down_total_dif, flux_sw_vis,                              &
               flux_sw_vis_dir, flux_sw_vis_dif,                                 &
               flux_lw, coszen, gust, Surf_diff, gavg_rrv, diff_cum_mom=diff_cu_mo, &
               moist_convect=convect)
endif

return
end subroutine spectral_physics_down
!------------------------------------------------------------------------------------------------
subroutine spectral_physics_up(Time_prev, Time, Time_next, previous, current, p_half, p_full,     &
                               z_half, z_full, wg_full, ug, vg, tg, grid_tracers,                 &
                               frac_land, dt_ug, dt_vg, dt_tg, dt_tracers, Surf_diff, lprec, fprec, gust)

type(time_type), intent(in) :: Time_prev, Time, Time_next
integer,         intent(in) :: previous, current
real, intent(in),    dimension(:,:,:    ) :: p_full, z_full, wg_full
real, intent(in),    dimension(:,:,:    ) :: p_half, z_half
real, intent(in),    dimension(:,:,:,:  ) :: ug, vg, tg
real, intent(in),    dimension(:,:,:,:,:) :: grid_tracers
real, intent(in),    dimension(:,:      ) :: frac_land
real, intent(inout), dimension(:,:,:    ) :: dt_ug, dt_vg, dt_tg
real, intent(inout), dimension(:,:,:,:  ) :: dt_tracers
type(surf_diff_type), intent(inout) :: Surf_diff
real, intent(out),   dimension(:,:) :: lprec, fprec
real, intent(inout), dimension(:,:) :: gust

if(.not.module_is_initialized) then
  call error_mesg('spectral_physics_up','spectral_physics module is not initialized', FATAL)
end if

call physics_driver_up(1, ie-is+1, 1, je-js+1, Time_prev, Time, Time_next,    &
                       rad_lat_2d, rad_lon_2d, area_2d,                       &
                       p_half, p_full, z_half, z_full,                        &
                       wg_full, ug(:,:,:,current), vg(:,:,:,current),         &
                       tg(:,:,:,current), grid_tracers(:,:,:,current,nhum),   &
                       grid_tracers(:,:,:,current,:),                         &
                       ug(:,:,:,previous), vg(:,:,:,previous),                &
                       tg(:,:,:,previous), grid_tracers(:,:,:,previous,nhum), &
                       grid_tracers(:,:,:,previous,:),                        &
                       frac_land, dt_ug, dt_vg, dt_tg,                        &
                       dt_tracers(:,:,:,nhum), dt_tracers, Surf_diff,         &
                       lprec, fprec, gust)

return
end subroutine spectral_physics_up
!------------------------------------------------------------------------------------------------
subroutine spectral_physics_moist(Time_next, delta_t, frac_land, p_half, p_full, z_half, z_full, wg_full, &
                                  tg, qg, ug, vg, tracers, lprec, fprec, gust)

type(time_type), intent(in) :: Time_next
real, intent(in)                        :: delta_t
real, intent(in),    dimension(:,:)     :: frac_land
real, intent(in),    dimension(:,:,:)   :: p_half, p_full, z_half, z_full, wg_full
real, intent(inout), dimension(:,:,:)   :: tg, qg, ug, vg
real, intent(inout), dimension(:,:,:,:) :: tracers
real, intent(inout), dimension(:,:)     :: lprec, fprec, gust
real, dimension(size(tg,1), size(tg,2), size(tg,3)) :: dt_tg, dt_qg, dt_ug, dt_vg
real, dimension(size(tracers,1), size(tracers,2), size(tracers,3), size(tracers,4)) :: dt_tracers
real, dimension(size(gust,1), size(gust,2)) :: gust_cv

dt_tg=0.; dt_qg=0.; dt_ug=0.; dt_vg=0.; dt_tracers=0.

if(do_moist_in_phys_up()) then
  if ( do_mcm_moist_processes ) then
    call error_mesg('spectral_physics_moist','do_mcm_moist_processes cannot be .true. when moist_processes'// &
                    ' is called by physics_driver_up', FATAL)
  endif
else
  if ( do_mcm_moist_processes ) then
    call mcm_moist_processes(1, ie-is+1, 1, je-js+1, Time_next, delta_t, p_half, p_full, tg, tracers(:,:,:,nhum), lprec, fprec)
  else
    call moist_processes(1, ie-is+1, 1, je-js+1, Time_next, delta_t, frac_land, p_half, p_full, z_half, z_full, wg_full, &
                         get_diff_t(), get_radturbten(),                                                                 &
                         tg, qg, tracers, ug, vg, tg, qg, tracers, ug, vg, dt_tg, dt_qg, dt_tracers, dt_ug, dt_vg,       &
                         diff_cu_mo, convect, lprec, fprec, gust_cv, area_2d, rad_lat_2d)
    gust = sqrt( gust*gust + gust_cv*gust_cv)
    call zero_radturbten()
    tg = tg + dt_tg*delta_t
    qg = qg + dt_qg*delta_t
    ug = ug + dt_ug*delta_t
    vg = vg + dt_vg*delta_t
    tracers = tracers + dt_tracers*delta_t
  endif
endif

return
end subroutine spectral_physics_moist
!------------------------------------------------------------------------------------------------

subroutine spectral_physics_end(Time)
character(len=64) :: file

type(time_type), intent(in) :: Time
real, dimension(is:ie,js:je)  :: rconvect

if(.not.module_is_initialized) return

file = 'RESTART/spectral_physics.res.nc'
call write_data(trim(file), 'diff_cu_mo', diff_cu_mo, grid_domain)
where(convect)
  rconvect = 1.
elsewhere
  rconvect = 0.
endwhere
call write_data(trim(file), 'convect', rconvect, grid_domain) ! pjp: No interface for reading/writing netcdf logicals

deallocate(rad_lon_2d, rad_lat_2d, area_2d, diff_cu_mo, convect)
call physics_driver_end(Time)
if(do_mcm_moist_processes) then
  call mcm_moist_processes_end
else
  call moist_processes_end
endif
module_is_initialized = .false.

return
end subroutine spectral_physics_end
!------------------------------------------------------------------------------------------------

end module spectral_physics_mod
