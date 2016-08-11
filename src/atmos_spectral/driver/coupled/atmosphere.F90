module atmosphere_mod

use               mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC

use               fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number, set_domain, &
                                 stdlog, close_file, open_namelist_file, check_nml_error, file_exist, field_size, &
                                 read_data, write_data

use    physics_driver_mod, only: do_moist_in_phys_up

use     field_manager_mod, only: MODEL_ATMOS

use    tracer_manager_mod, only: get_number_tracers, get_tracer_index

use  spectral_physics_mod, only: spectral_physics_init, spectral_physics_down, spectral_physics_up, &
                                 spectral_physics_end, surf_diff_type, spectral_physics_moist

use         constants_mod, only: grav

use        transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, get_deg_lon, get_deg_lat, &
                                 get_wts_lat, get_grid_boundaries, compute_ucos_vcos, divide_by_cos, get_lon_max, &
                                 get_lat_max, get_grid_domain, grid_domain

use  press_and_geopot_mod, only: pressure_variables, compute_pressures_and_heights, compute_z_bot

use      time_manager_mod, only: time_type, set_time, get_time, operator(+), operator(<), operator(-)

use spectral_dynamics_mod, only: spectral_dynamics_init, spectral_dynamics, spectral_dynamics_end, get_num_levels, &
                                 complete_robert_filter, complete_update_of_future, &
                                 get_axis_id, spectral_diagnostics, get_initial_fields

use       mpp_domains_mod, only: domain2d

use       tracer_type_mod, only: tracer_type

implicit none
private

character(len=128), parameter :: version = &
'$Id: atmosphere.F90,v 17.0 2009/07/21 03:00:28 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: tikal $'

public :: atmosphere_init, atmosphere_down, atmosphere_up, atmosphere_end, atmosphere_domain
public :: atmosphere_resolution, atmosphere_boundary, get_bottom_mass, get_bottom_wind, get_atmosphere_axes
public :: surf_diff_type, atmosphere_restart
integer :: seconds, days, num_tracers, num_levels, nhum
logical :: dry_model

integer, parameter :: num_time_levels=2
integer :: phyclock, dynclock

real, allocatable, dimension(:,:,:) :: p_half, p_full, z_half, z_full, wg_full
type(tracer_type), allocatable, dimension(:) :: tracer_attributes
real,    allocatable, dimension(:,:,:,:,:) :: grid_tracers
real,    allocatable, dimension(:,:,:    ) :: psg
real,    allocatable, dimension(:,:,:,:  ) :: ug, vg, tg

real, allocatable, dimension(:,:    ) :: dt_psg, z_bot
real, allocatable, dimension(:,:,:  ) :: dt_ug, dt_vg, dt_tg
real, allocatable, dimension(:,:,:,:) :: dt_tracers

! dt_real is the atmospheric time step, converted to a real number.
! delta_t is passed to physics. It is twice dt_real, except for
! the first time step of a cold start, when it equals dt_real.
real :: delta_t, dt_real

integer :: is, ie, js, je
integer :: previous, current, future
logical :: module_is_initialized=.false., atmos_domain_is_computed=.false.

type(time_type) :: Time_step, Time_prev, Time_next

logical :: do_mcm_moist_processes = .false.

namelist / atmosphere_nml / do_mcm_moist_processes

!------------------------------------------------------------------------------------------------

contains

!####################################################################################################################

subroutine atmosphere_init(Time_init, Time, Time_step_in, Surf_diff)

type(time_type),      intent(in)    :: Time_init, Time, Time_step_in
type(surf_diff_type), intent(inout) :: Surf_diff

integer :: ierr, io, unit, lon_max, lat_max, ntr, nt
integer, dimension(4) :: siz
character(len=64) :: file, tr_name
character(len=4) :: ch1,ch2,ch3,ch4

if(module_is_initialized) return

dynclock = mpp_clock_id('Dynamics', flags=MPP_CLOCK_SYNC)
phyclock = mpp_clock_id('Physics',  flags=MPP_CLOCK_SYNC)

unit = open_namelist_file()
ierr=1
do while (ierr /= 0)
  read(unit, nml=atmosphere_nml, iostat=io, end=20)
  ierr = check_nml_error (io, 'atmosphere_nml')
  enddo
20 call close_file (unit)

call write_version_number(version, tagname)
if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=atmosphere_nml)
!-----------------------------------------------------------------------------------------

!  because the time step is used in different ways,
!  it must exist as a real and time_type variable

call get_time(Time_step_in, seconds, days)
dt_real   = float(86400*days + seconds)
Time_step = Time_step_in

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
allocate (tracer_attributes(num_tracers))

call spectral_dynamics_init(Time, Time_step, tracer_attributes, dry_model, nhum)
atmos_domain_is_computed = .true.
call get_grid_domain(is, ie, js, je)
call get_num_levels(num_levels)

allocate (p_half       (is:ie, js:je, num_levels+1))
allocate (z_half       (is:ie, js:je, num_levels+1))
allocate (p_full       (is:ie, js:je, num_levels))
allocate (z_full       (is:ie, js:je, num_levels))
allocate (wg_full      (is:ie, js:je, num_levels))
allocate (psg          (is:ie, js:je, num_time_levels))
allocate (ug           (is:ie, js:je, num_levels, num_time_levels))
allocate (vg           (is:ie, js:je, num_levels, num_time_levels))
allocate (tg           (is:ie, js:je, num_levels, num_time_levels))
allocate (grid_tracers (is:ie, js:je, num_levels, num_time_levels, num_tracers ))
allocate (dt_psg       (is:ie, js:je))
allocate (dt_ug        (is:ie, js:je, num_levels))
allocate (dt_vg        (is:ie, js:je, num_levels))
allocate (dt_tg        (is:ie, js:je, num_levels))
allocate (dt_tracers   (is:ie, js:je, num_levels, num_tracers ))
allocate (z_bot        (is:ie, js:je))

p_half=0.; z_half=0.; p_full=0.; z_full=0.; wg_full=0.
psg=0.; ug=0.; vg=0.; tg=0.; grid_tracers=0.
dt_psg=0.; dt_ug=0.; dt_vg=0.; dt_tg=0.; dt_tracers=0.

file = 'INPUT/atmosphere.res.nc'
if(file_exist(trim(file))) then
  call get_lon_max(lon_max)
  call get_lat_max(lat_max)
  call field_size(trim(file), 'ug', siz)
  if(lon_max /= siz(1) .or. lat_max /= siz(2)) then
    write(ch1,'(i4)') siz(1)
    write(ch2,'(i4)') siz(2)
    write(ch3,'(i4)') lon_max
    write(ch4,'(i4)') lat_max
    call error_mesg('atmosphere_init','Resolution of restart data does not match resolution specified on namelist.'// &
    ' Restart data: lon_max='//ch1//', lat_max='//ch2//'  Namelist: lon_max='//ch3//', lat_max='//ch4, FATAL)
  endif
  call read_data(trim(file), 'previous', previous, no_domain=.true.)
  call read_data(trim(file), 'current',  current,  no_domain=.true.)
  do nt=1,num_time_levels
    call read_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'psg', psg(:,:,  nt), grid_domain, timelevel=nt)
    do ntr = 1,num_tracers
      tr_name = trim(tracer_attributes(ntr)%name)
      call read_data(trim(file), trim(tr_name), grid_tracers(:,:,:,nt,ntr), grid_domain, timelevel=nt)      
    enddo ! end loop over tracers
  enddo ! end loop over time levels
  call read_data(trim(file), 'wg_full', wg_full, grid_domain)
else if(file_exist('INPUT/atmosphere.res')) then
  call error_mesg('atmosphere_init', &
                  'Binary restart file, INPUT/atmosphere.res, is not supported by this version of atmosphere.f90',FATAL)
else
  previous = 1; current = 1
  call get_initial_fields(ug(:,:,:,1), vg(:,:,:,1), tg(:,:,:,1), psg(:,:,1), grid_tracers(:,:,:,1,:))
endif

call spectral_physics_init(Time, get_axis_id(), Surf_diff, nhum, p_half, do_mcm_moist_processes)

if(dry_model) then
  call compute_pressures_and_heights(tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half)
else
  call compute_pressures_and_heights( &
       tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half, grid_tracers(:,:,:,current,nhum))
endif

call compute_z_bot(psg(:,:,current), tg(:,:,num_levels,current), z_bot, grid_tracers(:,:,num_levels,current,nhum))

module_is_initialized = .true.

return
end subroutine atmosphere_init
!#################################################################################################################################

subroutine atmosphere_down(Time, frac_land, t_surf, albedo,                                              &
                      albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif,                    &
                      rough_mom, u_star, b_star, q_star, dtau_du, dtau_dv, tau_x, tau_y,                 &
                      gust, coszen, flux_sw, flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir,             &
                      flux_sw_down_vis_dif, flux_sw_down_total_dir, flux_sw_down_total_dif, flux_sw_vis, &
                      flux_sw_vis_dir, flux_sw_vis_dif, flux_lw, Surf_diff )

type(time_type),      intent(in) :: Time
real,                 intent(in),    dimension(:,:) :: frac_land, t_surf, albedo
real,                 intent(in),    dimension(:,:) :: albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif 
real,                 intent(in),    dimension(:,:) :: rough_mom, u_star, b_star, q_star, dtau_du, dtau_dv
real,                 intent(inout), dimension(:,:) :: tau_x, tau_y
real,                 intent(out),   dimension(:,:) :: flux_sw, flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir
real,                 intent(out),   dimension(:,:) :: flux_sw_down_vis_dif, flux_sw_down_total_dir, flux_sw_down_total_dif
real,                 intent(out),   dimension(:,:) :: flux_sw_vis, flux_sw_vis_dir, flux_sw_vis_dif, flux_lw, coszen, gust
type(surf_diff_type), intent(inout)                 :: Surf_diff

if(.not.module_is_initialized) then
  call error_mesg('atmosphere_down','atmosphere module has not been initialized.', FATAL)
endif

dt_psg = 0.
dt_ug  = 0.
dt_vg  = 0.
dt_tg  = 0.
dt_tracers = 0.

if(current == previous) then
  delta_t = dt_real
  Time_prev = Time
else
  delta_t = 2*dt_real
  Time_prev = Time - Time_step
endif
Time_next = Time + Time_step

call mpp_clock_begin(phyclock)
call spectral_physics_down(Time_prev, Time, Time_next, previous, current, p_half, p_full, z_half, z_full, psg,      &
                        ug, vg, tg, grid_tracers, frac_land, rough_mom, albedo, t_surf, u_star, b_star, q_star,     &
                        dtau_du, dtau_dv, tau_x, tau_y, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif, &
                        dt_ug, dt_vg, dt_tg, dt_tracers, flux_sw, flux_sw_dir, flux_sw_dif,                         &
                        flux_sw_down_vis_dir, flux_sw_down_vis_dif, flux_sw_down_total_dir, flux_sw_down_total_dif, &
                        flux_sw_vis, flux_sw_vis_dir, flux_sw_vis_dif, flux_lw, coszen,gust, Surf_diff)
call mpp_clock_end(phyclock)

return
end subroutine atmosphere_down
!#################################################################################################################################

subroutine atmosphere_up(Time, frac_land, Surf_diff, lprec, fprec, gust)

type(time_type),      intent(in)                          :: Time
real,                 intent(in),  dimension(is:ie,js:je) :: frac_land
type(surf_diff_type), intent(inout)                       :: Surf_diff
real,                 intent(out), dimension(is:ie,js:je) :: lprec, fprec, gust

real, dimension(is:ie,js:je,1:num_levels+1) :: ln_p_half
real, dimension(is:ie,js:je,1:num_levels  ) :: ln_p_full

if(.not.module_is_initialized) then
  call error_mesg('atmosphere_up','atmosphere module has not been initialized.', FATAL)
endif

call mpp_clock_begin(phyclock)
call spectral_physics_up(Time_prev, Time, Time_next, previous, current, p_half, p_full, z_half, z_full, wg_full, ug, vg, tg, &
                         grid_tracers, frac_land, dt_ug, dt_vg, dt_tg, dt_tracers, Surf_diff, lprec, fprec, gust)
call mpp_clock_end(phyclock)

if(previous == current) then
  future = num_time_levels + 1 - current
else
  future = previous
endif

call mpp_clock_begin(dynclock)
call spectral_dynamics(Time, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                       tg(:,:,:,future), tracer_attributes, grid_tracers(:,:,:,:,:), future, &
                       dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, p_full, p_half, z_full)
call mpp_clock_end(dynclock)

if(do_moist_in_phys_up()) then
  if ( do_mcm_moist_processes ) then
    call error_mesg('atmosphere_up','do_mcm_moist_processes cannot be .true. when moist_processes'// &
                    ' is called by physics_driver_up', FATAL)
  endif
else
  call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, psg(:,:,future))
  call spectral_physics_moist(Time_next, delta_t, frac_land, p_half, p_full, z_half, z_full, wg_full, &
                              tg(:,:,:,future), grid_tracers(:,:,:,future,nhum), ug(:,:,:,future), vg(:,:,:,future), &
                              grid_tracers(:,:,:,future,:), lprec, fprec, gust)
  call complete_update_of_future(psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), tg(:,:,:,future), &
                              tracer_attributes, grid_tracers(:,:,:,future,:))
endif
call complete_robert_filter(tracer_attributes)

call spectral_diagnostics(Time_next, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                          tg(:,:,:,future), wg_full, grid_tracers(:,:,:,:,:), future)

previous = current
current  = future

call compute_pressures_and_heights( &
            tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half, grid_tracers(:,:,:,current,nhum))
call compute_z_bot(psg(:,:,current), tg(:,:,num_levels,current), z_bot, grid_tracers(:,:,num_levels,current,nhum))

return
end subroutine atmosphere_up
!####################################################################################################################

subroutine get_bottom_mass (t_bot, tr_bot, p_bot, z_bot_out, p_surf)

real, intent(out), dimension(:,:) :: t_bot, p_bot, z_bot_out, p_surf
real, intent(out), dimension(:,:,:) :: tr_bot

if(.not.module_is_initialized) then
  call error_mesg('get_bottom_mass','atmosphere module has not been initialized.', FATAL)
endif

t_bot     = tg(:,:,num_levels, previous)
tr_bot    = grid_tracers(:,:,num_levels, previous,:)
p_bot     = p_full(:,:,num_levels)
p_surf    = psg(:,:,previous)
z_bot_out = z_bot

return
end subroutine get_bottom_mass
!####################################################################################################################

subroutine get_bottom_wind (u_bot, v_bot)

real, intent(out), dimension(:,:) :: u_bot, v_bot

if(.not.module_is_initialized) then
  call error_mesg('get_bottom_wind','atmosphere module has not been initialized.', FATAL)
endif

u_bot = ug(:,:,num_levels, previous)
v_bot = vg(:,:,num_levels, previous)

return
end subroutine get_bottom_wind
!####################################################################################################################

subroutine atmosphere_resolution(num_lon_out, num_lat_out, global)

integer, intent(out)          :: num_lon_out, num_lat_out
logical, intent(in), optional :: global
logical :: global_tmp

if(.not.module_is_initialized) then
  call error_mesg('atmosphere_resolution','atmosphere module has not been initialized.', FATAL)
endif

if (present(global)) then
  global_tmp = global
else
  global_tmp = .false.
endif

if(global_tmp) then
  call get_lon_max(num_lon_out)
  call get_lat_max(num_lat_out)
else
  num_lon_out = ie-is+1
  num_lat_out = je-js+1
endif

return
end subroutine atmosphere_resolution
!####################################################################################################################

subroutine get_atmosphere_axes(axes_out)
integer, intent(out), dimension(:) :: axes_out

if(.not.module_is_initialized) then
  call error_mesg('get_atmosphere_axes','atmosphere module has not been initialized.', FATAL)
endif

axes_out = get_axis_id()

return
end subroutine get_atmosphere_axes
!####################################################################################################################

subroutine atmosphere_boundary(lon_boundaries, lat_boundaries, global)

real,    intent(out), dimension(:) :: lon_boundaries, lat_boundaries
logical, intent(in),  optional     :: global

logical :: global_tmp

if(.not.module_is_initialized) then
  call error_mesg('atmosphere_boundary','atmosphere module has not been initialized.', FATAL)
endif

if(present(global)) then
  global_tmp = global
else
  global_tmp = .false.
endif
call get_grid_boundaries(lon_boundaries, lat_boundaries, global_tmp)

return
end subroutine atmosphere_boundary
!####################################################################################################################
subroutine atmosphere_domain(domain)
type (domain2d), intent(out) :: domain

if(.not.atmos_domain_is_computed) then
  call error_mesg('atmosphere_domain','spec_mpp has not been initialized.', FATAL)
endif

domain = grid_domain

end subroutine atmosphere_domain
!####################################################################################################################

subroutine atmosphere_end(Time)
type(time_type), intent(in) :: Time
integer :: ntr, nt
character(len=64) :: file, tr_name

if(.not.module_is_initialized) return

file='RESTART/atmosphere.res'
call write_data(trim(file), 'previous', previous, no_domain=.true.)
call write_data(trim(file), 'current',  current,  no_domain=.true.)
do nt=1,num_time_levels
  call write_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'psg', psg(:,:,  nt), grid_domain)
  do ntr = 1,num_tracers
    tr_name = trim(tracer_attributes(ntr)%name)
    call write_data(trim(file), tr_name, grid_tracers(:,:,:,nt,ntr), grid_domain)
  enddo
enddo
call write_data(trim(file), 'wg_full', wg_full, grid_domain)

deallocate(dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers)

call set_domain(grid_domain)
call spectral_physics_end(Time)
call spectral_dynamics_end(tracer_attributes, Time)

module_is_initialized = .false.

return
end subroutine atmosphere_end

  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  dummy routine.
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call error_mesg ('atmosphere_restart in atmosphere_mod', &
                     'intermediate restart capability is not implemented for this model', FATAL)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>


!####################################################################################################################

end module atmosphere_mod
