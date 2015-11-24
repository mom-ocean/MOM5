module fv_physics_mod

! Note: the option "use_tendency = .true." will not be reproducible
!       this is because winds on the physics grid need to be saved
!       in the restart file. This is a to-do item since there is no
!       clear benefit at this point for doing so.

!-----------------------------------------------------------------------
!
!        interface for FV dynamics with atmospheric physics
!
!-----------------------------------------------------------------------
use      time_manager_mod, only: time_type, get_time, operator(-)
use               fms_mod, only: error_mesg, FATAL, write_version_number, &
                                 clock_flag_default, NOTE, WARNING
use    physics_driver_mod, only: physics_driver_init,   &
                                 physics_driver_end,    &
                                 physics_driver_moist_init, &
                                 physics_driver_moist_end, &
                                 physics_driver_down_time_vary, &
                                 physics_driver_up_time_vary,  &
                                 physics_driver_down_endts,  &
                                 physics_driver_up_endts, &
                                 physics_driver_down,   &
                                 physics_driver_up,     &
                                 surf_diff_type
use     field_manager_mod, only: MODEL_ATMOS
use    tracer_manager_mod, only: get_tracer_index, NO_TRACER

!use    fv_pack,       only: nlon, mlat, nlev, beglat, endlat,          &
!                            nt_phys, nt_prog, ncnst, beglon, endlon,   &
!                            u, v, pt, q, ua, va, delp, phis, pe, peln, &
!                            omga, rlat, rlon, area, rlonb, rlatb,      &
!                            u_phys, v_phys, t_phys, q_phys,            &
!                            use_tendency, get_eta_level
use    fv_pack,       only: nlon, mlat, beglat, endlat,          &
                            nt_phys, nt_prog, beglon, endlon,   &
                            rlat, rlon, area, rlonb, rlatb,      &
                            use_tendency, get_eta_level, fv_domain
use     timingModule, only: timing_on, timing_off

use    atmos_co2_mod, only: atmos_co2_rad, co2_radiation_override
use  atmos_nudge_mod, only: atmos_nudge_init, atmos_nudge_end
use    constants_mod, only: rdgas, grav, rvgas, WTMAIR, WTMCO2

use mpp_mod, only: stdout, stderr, mpp_error, mpp_clock_id, &
     mpp_clock_begin, mpp_clock_end, CLOCK_MODULE_DRIVER
use update_fv_phys_mod, only: update_fv_phys
use fv_arrays_mod, only: fv_array_limits, fv_array_sync
use    mpp_mod,       only: mpp_pe, mpp_root_pe
use   mpp_domains_mod,only: mpp_global_sum, BITWISE_EXACT_SUM

#ifdef _OPENMP
use omp_lib
#endif
!-----------------------------------------------------------------------

implicit none
private

public   fv_physics_down, fv_physics_up,  &
         fv_physics_init, fv_physics_end
public   surf_diff_type

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id: fv_physics.F90,v 20.0 2013/12/13 23:08:17 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

   real zvir, rrg, ginv
   logical :: do_atmos_nudge
!RSH make module variable:
   integer :: sphum

   integer :: isw, iew, jsw, jew !window start/end in global index space
   integer :: nx_win, ny_win !iew-isw+1, jew-jsw+1 (window sizes)
   integer :: nx_dom, ny_dom !ie-is+1, je-js+1 (compute domain sizes)
!MPP clocks
   integer :: id_fv_physics_down, id_fv_physics_up, id_update_fv_phys

! Auto arrays:
   real, allocatable, dimension(:,:,:) :: p_full, z_full, p_half, z_half
   integer, allocatable, dimension(:)  :: physics_window_x, physics_window_y
   integer :: numthreads, ny_per_thread, num_phys_windows
!#include "omp.h"

contains


  subroutine fv_physics_init (axes, Time, window, Surf_diff)
#include "fv_arrays.h"
!-----------------------------------------------------------------------
!
!   axes      = array of axis indices for diagnostics (x,y,pf,ph)
!   Time      = current time (time_type)
!
!-----------------------------------------------------------------------
    integer,               intent(in)    :: axes(4)
    type (time_type),      intent(in)    :: Time
    integer,               intent(in)    :: window(2)
    type (surf_diff_type), intent(inout) :: Surf_diff
!-----------------------------------------------------------------------
    real, dimension(nlev+1,2) :: pref
!-----------------------------------------------------------------------
!RSHinteger :: i, j, k, sphum
    integer :: i, j, k
    integer :: ios
    character(len=80) evalue
!-----------------------------------------------------------------------
    real p_edge(nlon,beglat:endlat,nlev+1)  ! used by atmos_tracer_driver_init
    real phalf(nlev+1)
    character(len=132) :: text
    real, allocatable :: rlonb2d(:,:), rlatb2d(:,:)
    integer :: unit
#include "fv_point.inc"

    zvir = rvgas/rdgas - 1.
    ginv = 1./ grav
    rrg = rdgas / grav        

!----- write version to logfile --------

    call write_version_number(version, tagname)

!---------- reference profile -----------

    pref(nlev+1,1) = 101325.
    pref(nlev+1,2) = 81060.

    call get_eta_level ( nlev, pref(nlev+1,1), pref(1,1), phalf )
    call get_eta_level ( nlev, pref(nlev+1,2), pref(1,2), phalf )

!------- pressure at model layer interfaces -----

    do k=1,nlev+1
       do j=beglat,endlat
          do i=beglon,endlon
             p_edge(i,j,k) = pe(i,k,j)
          enddo
       enddo
    enddo

!---------- initialize physics -------
    allocate( rlonb2d( size(rlonb), size(rlatb) ) )
    allocate( rlatb2d( size(rlonb), size(rlatb) ) )
    do j = 1, size(rlatb)
       rlonb2d(:,j) = rlonb(:)
    end do
    do i = 1, size(rlonb)
       rlatb2d(i,:) = rlatb(:)
    end do

    call physics_driver_init(Time, &
         rlonb2d(beglon:endlon+1,beglat:endlat+1), &
         rlatb2d(beglon:endlon+1,beglat:endlat+1), &
         rlonb2d(beglon:endlon+1,beglat:endlat+1), & !dummy argument needed for interface change
         rlatb2d(beglon:endlon+1,beglat:endlat+1), & !dummy argument needed for interface change
         axes, pref, q(beglon:endlon,beglat:endlat,:,1:ncnst), &
         Surf_diff, p_edge(beglon:endlon,:,:) )
!    call physics_driver_init(Time, rlonb, rlatb(beglat:endlat+1), &
!                             axes, pref, q(:,beglat:endlat,:,1:ncnst), &
!                             Surf_diff, p_edge )

! Specific humidity is assumed to be q(:,:,:,1)
    sphum = get_tracer_index (MODEL_ATMOS, 'sphum' )
    if(sphum /= 1) call error_mesg('fv_physics_init:','sphum /= 1', FATAL)

!--- initialize nudging module ---
    call atmos_nudge_init ( Time, axes(1:3), flag=do_atmos_nudge )
!physics window
    nx_win = window(1)
    ny_win = window(2)
    nx_dom = endlon - beglon + 1
    ny_dom = endlat - beglat + 1
    if( nx_win.LE.0 )nx_win = nx_dom
    if( ny_win.LE.0 )ny_win = ny_dom
    if( mod(nx_dom,nx_win).NE.0 )then
        write( text,'(a,2i4)' )'FV_PHYSICS_INIT: atmosphere_nml problem,'// &
             ' physics_window must divide domain size: ',  nx_win, nx_dom
        call mpp_error( FATAL, text )
    end if
    if( mod(ny_dom,ny_win).NE.0 )then
        write( text,'(a,2i4)' )'FV_PHYSICS_INIT: atmosphere_nml problem,'// &
             ' physics_window must divide domain size: ',  ny_win, ny_dom
        call mpp_error( FATAL, text )
    end if
    allocate( p_full(beglon:endlon,beglat:endlat,nlev) )
    allocate( z_full(beglon:endlon,beglat:endlat,nlev) )
    allocate( p_half(beglon:endlon,beglat:endlat,nlev+1) )
    allocate( z_half(beglon:endlon,beglat:endlat,nlev+1) )
!MPP clocks
    id_fv_physics_down = mpp_clock_id( 'FV_PHYSICS_DOWN', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )
    id_fv_physics_up = mpp_clock_id( 'FV_PHYSICS_UP', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )
    id_update_fv_phys = mpp_clock_id( 'UPDATE_FV_PHYS', &
         flags=clock_flag_default, grain=CLOCK_MODULE_DRIVER )

    numthreads = 1
!$OMP PARALLEL
!$OMP MASTER
!$       numthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL

    if ( mpp_pe()==mpp_root_pe() ) then
      unit = stdout()
      write(unit,*) 'Starting Threads : ', numthreads
    endif
    ny_per_thread = max (1,ny_win/numthreads)


    num_phys_windows = (nx_dom/nx_win)*(ny_dom/ny_win) 
    write(text,'(a,2i4)') 'num_phys_windows, numthreads',num_phys_windows,numthreads
    call error_mesg ('fv_physics_init', trim(text), NOTE)
    allocate(physics_window_x(num_phys_windows))
    allocate(physics_window_y(num_phys_windows))
    i = 1 
    do jsw = beglat,endlat,ny_win
       do isw = beglon,endlon,nx_win
          physics_window_x(i) =isw
          physics_window_y(i) =jsw
          i = i + 1
       enddo
    enddo
    
    
      end subroutine fv_physics_init


  subroutine fv_physics_down(dt_phys,Time_prev, Time, Time_next, &
       frac_land,   albedo,        &
       albedo_vis_dir, albedo_nir_dir, &
       albedo_vis_dif, albedo_nir_dif, &
       rough_vel,   t_surf,        &
       u_star, b_star, q_star,     &
       dtau_du, dtau_dv, tau_x, tau_y, &
       flux_sw,                    &
       flux_sw_dir, flux_sw_dif,   &
       flux_sw_down_vis_dir,       &
       flux_sw_down_vis_dif,       &
       flux_sw_down_total_dir,     &
       flux_sw_down_total_dif,     &
       flux_sw_vis, flux_sw_vis_dir, &
       flux_sw_vis_dif,            &
       flux_lw, coszen,            &
       gust, Surf_diff, frac_open_sea )
#include "fv_arrays.h"

!-----------------------------------------------------------------------
!
!   Time_prev =  time at the previous time level, tau-1 (time_type)
!   Time      =  time at the current time level,  tau   (time_type)
!   Time_next =  time at the next time level,     tau+1 (time_type)
!
!   NOTE: for a two time level scheme (e.g., forward-backward scheme)
!         Time_prev = Time.
!
!-----------------------------------------------------------------------
    real,    intent(in)           :: dt_phys
    type(time_type),intent(in)    :: Time_prev, Time, Time_next

    real, intent(in),    dimension(beglon:endlon,beglat:endlat) :: &
         frac_land,  albedo,   &
         albedo_vis_dir, albedo_nir_dir, &
         albedo_vis_dif, albedo_nir_dif, &
         rough_vel,  t_surf,   &
         u_star,     b_star, &
         q_star,     dtau_du, dtau_dv, frac_open_sea

    real, intent(inout), dimension(beglon:endlon,beglat:endlat) :: &
         tau_x, tau_y
    real, intent(out),   dimension(beglon:endlon,beglat:endlat) :: &
         flux_sw,  flux_sw_dir, flux_sw_dif,&
         flux_sw_down_vis_dir,       &
         flux_sw_down_vis_dif,       &
         flux_sw_down_total_dir,     &
         flux_sw_down_total_dif,     &
         flux_sw_vis,                &
         flux_sw_vis_dir,            &
         flux_sw_vis_dif,            &
         flux_lw, coszen, gust

    type(surf_diff_type), intent(inout):: Surf_diff

!-----------------------------------------------------------------------
    integer :: i, j, k, m, idx, phys_loop
    real :: rdt, gavg_rrv(nt_prog)
    integer :: is_w, ie_w, js_w, je_w
    real    :: dt 
    integer :: sec, day

#include "fv_point.inc"

!---------------------------- do physics -------------------------------

    rdt = 1. / dt_phys

!----------------------------------------------------------------------
! obtain pressure-weighted global mean co2 dry volume mixing ratio for
! use by radiation package.
!----------------------------------------------------------------------
    gavg_rrv = 0.
! check if to override predicted global pressure-weighted rad co2
    idx = get_tracer_index(MODEL_ATMOS, 'co2')
    if(idx /= NO_TRACER .and. co2_radiation_override) then
      call atmos_co2_rad(Time, gavg_rrv(idx))
    elseif (idx /= NO_TRACER) then
      call compute_g_avg(gavg_rrv, 'co2')
    endif  
!---------------------------------------------------------------------
! compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time (Time_next-Time_prev, sec, day)
    dt = real(sec+day*86400)
 
!---------------------------------------------------------------------
! call physics_driver_down_time_vary to do the time-dependent, spatially
! independent calculations before entering windows / threads loop. 
!--------------------------------------------------------------------- 
    call physics_driver_down_time_vary (Time, Time_next, gavg_rrv, dt)


    call mpp_clock_begin(id_fv_physics_down)

!isw,iew,jsw,jew are in global index space
    call compute_p_z( beglon, beglat, nx_dom, ny_dom )
!$OMP parallel do schedule(dynamic) default(shared) private(phys_loop, isw, iew, jsw, jew, m, k, j, i)
    do phys_loop = 1, size(physics_window_y)!num_phys_windows
       jsw = physics_window_y(phys_loop)
       jew = jsw + ny_win - 1
       isw = physics_window_x(phys_loop)
       iew = isw + nx_win - 1
          if ( use_tendency ) then
              do k=1,nlev
                 do j = jsw,jew
                    do i=isw,iew
                       u_dt(i,j,k) = (ua(i,j,k) - u_phys(i,j,k)) * rdt
                       v_dt(i,j,k) = (va(i,j,k) - v_phys(i,j,k)) * rdt
                       t_dt(i,j,k) = (pt(i,j,k) - t_phys(i,j,k)) * rdt
                    end do
                 enddo
              enddo

              do m=1,nt_prog
                 do k=1,nlev
                    do j = jsw,jew
                       do i=isw,iew
                          q_dt(i,j,k,m) = (q(i,j,k,m)-q_phys(i,j,k,m)) * rdt
                       enddo
                    end do
                 enddo
              enddo
!isw-beglon+1 are wrt domain index space
              call physics_driver_down                                       &
                   ( isw-beglon+1, iew-beglon+1, jsw-beglat+1, jew-beglat+1, &
                   Time_prev, Time, Time_next                              , &
                   rlat(isw:iew,jsw:jew) , rlon(isw:iew,jsw:jew)           , &
                   area(isw:iew,jsw:jew) , &
                   p_half(isw-is_w+1:iew-is_w+1,jsw-js_w+1:jew-js_w+1,:), &
                   p_full(isw-is_w+1:iew-is_w+1,jsw-js_w+1:jew-js_w+1,:), &
                   z_half(isw-is_w+1:iew-is_w+1,jsw-js_w+1:jew-js_w+1,:), &
                   z_full(isw-is_w+1:iew-is_w+1,jsw-js_w+1:jew-js_w+1,:), &
                   p_half(isw-is_w+1:iew-is_w+1,jsw-js_w+1:jew-js_w+1,:), &
                   u_phys(isw:iew,jsw:jew,:), v_phys(isw:iew,jsw:jew,:)    , &
                   t_phys(isw:iew,jsw:jew,:), q_phys(isw:iew,jsw:jew,:,1)  , &
                   q_phys(isw:iew,jsw:jew,:,:)                             , &
                   u_phys(isw:iew,jsw:jew,:), v_phys(isw:iew,jsw:jew,:)    , &
                   t_phys(isw:iew,jsw:jew,:), q_phys(isw:iew,jsw:jew,:,1)  , &
                   q_phys(isw:iew,jsw:jew,:,:)                             , &
                   frac_land(isw:iew,jsw:jew), rough_vel(isw:iew,jsw:jew)  , &
                   frac_open_sea(isw:iew,jsw:jew)                          , &
                   albedo   (isw:iew,jsw:jew)                              , &
                   albedo_vis_dir(isw:iew,jsw:jew)                         , &
                   albedo_nir_dir(isw:iew,jsw:jew)                         , &
                   albedo_vis_dif(isw:iew,jsw:jew)                         , &
                   albedo_nir_dif(isw:iew,jsw:jew)                         , &
                   t_surf   (isw:iew,jsw:jew), u_star   (isw:iew,jsw:jew)  , &
                   b_star   (isw:iew,jsw:jew), q_star   (isw:iew,jsw:jew)  , &
                   dtau_du  (isw:iew,jsw:jew), dtau_dv  (isw:iew,jsw:jew)  , &
                   tau_x    (isw:iew,jsw:jew), tau_y    (isw:iew,jsw:jew)  , &
                   u_dt     (isw:iew,jsw:jew,:), v_dt(isw:iew,jsw:jew,:)   , &
                   t_dt     (isw:iew,jsw:jew,:), q_dt(isw:iew,jsw:jew,:,1) , &
                   q_dt     (isw:iew,jsw:jew,:,1:nt_prog)                  , &
                   flux_sw               (isw:iew,jsw:jew)                 , &
                   flux_sw_dir           (isw:iew,jsw:jew)                 , &
                   flux_sw_dif           (isw:iew,jsw:jew)                 , &
                   flux_sw_down_vis_dir  (isw:iew,jsw:jew)                 , &
                   flux_sw_down_vis_dif  (isw:iew,jsw:jew)                 , &
                   flux_sw_down_total_dir(isw:iew,jsw:jew)                 , &
                   flux_sw_down_total_dif(isw:iew,jsw:jew)                 , &
                   flux_sw_vis           (isw:iew,jsw:jew)                 , &
                   flux_sw_vis_dir       (isw:iew,jsw:jew)                 , &
                   flux_sw_vis_dif       (isw:iew,jsw:jew)                 , &
                   flux_lw               (isw:iew,jsw:jew)                 , &
                   coszen                (isw:iew,jsw:jew)                 , &
                   gust                  (isw:iew,jsw:jew)                 , &
                   Surf_diff, gavg_rrv                                     )
          else
              do k=1,nlev
                 do j = jsw,jew
                    do i=isw,iew
                       u_dt(i,j,k) = 0.
                       v_dt(i,j,k) = 0.
                       t_dt(i,j,k) = 0.
                    enddo
                 end do
              enddo

              do m=1,nt_prog
                 do k=1,nlev
                    do j = jsw,jew
                       do i=isw,iew
                          q_dt(i,j,k,m) = 0.
                       enddo
                    end do
                 enddo
              enddo
!isw-beglon+1 are wrt domain index space
              call physics_driver_down                                       &
                   ( isw-beglon+1, iew-beglon+1, jsw-beglat+1, jew-beglat+1, &
                   Time_prev, Time, Time_next                              , &
                   rlat(isw:iew,jsw:jew), rlon(isw:iew,jsw:jew), &
                   area(isw:iew,jsw:jew) , &
                   p_half(isw:iew,jsw:jew,:),&
                   p_full(isw:iew,jsw:jew,:),&
                   z_half(isw:iew,jsw:jew,:),&
                   z_full(isw:iew,jsw:jew,:),&
                   p_half(isw:iew,jsw:jew,:),&
                   ua (isw:iew,jsw:jew,:), va (isw:iew,jsw:jew,:)          , &
                   pt (isw:iew,jsw:jew,:), q  (isw:iew,jsw:jew,:,1)        , &
                   q  (isw:iew,jsw:jew,:,:)                                , &
                   ua (isw:iew,jsw:jew,:), va (isw:iew,jsw:jew,:)          , &
                   pt (isw:iew,jsw:jew,:), q  (isw:iew,jsw:jew,:,1)        , &
                   q  (isw:iew,jsw:jew,:,:)                                , &
                   frac_land(isw:iew,jsw:jew), rough_vel(isw:iew,jsw:jew), &
                   frac_open_sea(isw:iew,jsw:jew)                          , &
                   albedo   (isw:iew,jsw:jew)     ,&
                   albedo_vis_dir(isw:iew,jsw:jew),&
                   albedo_nir_dir(isw:iew,jsw:jew),&
                   albedo_vis_dif(isw:iew,jsw:jew),&
                   albedo_nir_dif(isw:iew,jsw:jew),&
                   t_surf   (isw:iew,jsw:jew), u_star   (isw:iew,jsw:jew),&
                   b_star   (isw:iew,jsw:jew), q_star   (isw:iew,jsw:jew),&
                   dtau_du  (isw:iew,jsw:jew), dtau_dv  (isw:iew,jsw:jew),&
                   tau_x    (isw:iew,jsw:jew), tau_y    (isw:iew,jsw:jew),&
                   u_dt     (isw:iew,jsw:jew,:), v_dt(isw:iew,jsw:jew,:)   , &
                   t_dt     (isw:iew,jsw:jew,:), q_dt(isw:iew,jsw:jew,:,1) , &
                   q_dt     (isw:iew,jsw:jew,:,1:nt_prog)                  , &
                   flux_sw               (isw:iew,jsw:jew),&
                   flux_sw_dir           (isw:iew,jsw:jew),&
                   flux_sw_dif           (isw:iew,jsw:jew),&
                   flux_sw_down_vis_dir  (isw:iew,jsw:jew),&
                   flux_sw_down_vis_dif  (isw:iew,jsw:jew),&
                   flux_sw_down_total_dir(isw:iew,jsw:jew),&
                   flux_sw_down_total_dif(isw:iew,jsw:jew),&
                   flux_sw_vis           (isw:iew,jsw:jew),&
                   flux_sw_vis_dir       (isw:iew,jsw:jew),&
                   flux_sw_vis_dif       (isw:iew,jsw:jew),&
                   flux_lw               (isw:iew,jsw:jew),&
                   coszen                (isw:iew,jsw:jew),&
                   gust                  (isw:iew,jsw:jew),&
                   Surf_diff, gavg_rrv                                     )
endif
       enddo

    call physics_driver_down_endts (1, 1) !Note that these arguments are not used yet.

    call fv_array_sync
    call mpp_clock_end(id_fv_physics_down)


  end subroutine fv_physics_down


  subroutine fv_physics_up(dt_phys, Time_prev, Time, Time_next, &
       frac_land, Surf_diff, lprec, fprec, gust, u_star, b_star, q_star )
#include "fv_arrays.h"
!-----------------------------------------------------------------------
!
!   Time_prev =  time at the previous time level, tau-1 (time_type)
!   Time      =  time at the current time level,  tau   (time_type)
!   Time_next =  time at the next time level,     tau+1 (time_type)
!
!   NOTE: for a two time level scheme (e.g., forward-backward scheme)
!         Time_prev = Time.
!
!-----------------------------------------------------------------------
    real,           intent(in)    :: dt_phys
    type(time_type),intent(in)    :: Time_prev, Time, Time_next
    real, intent(in),  dimension(beglon:endlon,beglat:endlat)  :: frac_land
    real, intent(out), dimension(beglon:endlon,beglat:endlat)  :: lprec, fprec
    real, intent(inout), dimension(beglon:endlon,beglat:endlat) :: gust
    type(surf_diff_type), intent(inout)                         :: Surf_diff
    real,    intent(in),dimension(beglon:endlon,beglat:endlat)  :: u_star, b_star, q_star

    integer :: is_w, ie_w, js_w, je_w, phys_loop
    integer :: sec, day, npz
    real    :: dt

#include "fv_point.inc"

    call mpp_clock_begin(id_fv_physics_up)
    npz = size(p_full,3)
!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time (Time_next-Time_prev, sec, day)
    dt = real(sec+day*86400)
 
    call physics_driver_up_time_vary (Time, Time_next, dt)

    call compute_p_z( beglon, beglat, nx_dom, ny_dom )
    call physics_driver_moist_init (nx_dom, ny_dom,  npz, nt_prog, ncnst) 
!$OMP parallel  do default(shared) private(phys_loop, isw, iew, jsw, jew)
    do phys_loop = 1, size(physics_window_y)!num_phys_windows
       jsw = physics_window_y(phys_loop)
       jew = jsw + ny_win - 1  
       isw = physics_window_x(phys_loop)
       iew = isw + nx_win - 1
          if ( use_tendency ) then

              call physics_driver_up (isw-beglon+1, iew-beglon+1, jsw-beglat+1, jew-beglat+1, &
                   Time_prev, Time, Time_next             , &
                   rlat     (isw:iew,jsw:jew)             , &
                   rlon     (isw:iew,jsw:jew)             , &
                   area     (isw:iew,jsw:jew)             , &
              p_half(isw:iew,jsw:jew,:), p_full(isw:iew,jsw:jew,:), &
              z_half(isw:iew,jsw:jew,:), z_full(isw:iew,jsw:jew,:), &
                   omga     (isw:iew,jsw:jew,:)           , &
                   u_phys   (isw:iew,jsw:jew,:)           , &
                   v_phys   (isw:iew,jsw:jew,:)           , &
                   t_phys   (isw:iew,jsw:jew,:)           , &
                   q_phys   (isw:iew,jsw:jew,:,1)         , &
                   q_phys   (isw:iew,jsw:jew,:,:)         , & ! cjg: pass all tracers
                   u_phys   (isw:iew,jsw:jew,:)           , &
                   v_phys   (isw:iew,jsw:jew,:)           , &
                   t_phys   (isw:iew,jsw:jew,:)           , &
                   q_phys   (isw:iew,jsw:jew,:,1)         , &
                   q_phys   (isw:iew,jsw:jew,:,:)         , & ! cjg: pass all tracers
                   frac_land(isw:iew,jsw:jew)             , &
                   u_star   (isw:iew,jsw:jew)             , &
                   b_star   (isw:iew,jsw:jew)             , &
                   q_star   (isw:iew,jsw:jew)             , &
                   u_dt     (isw:iew,jsw:jew,:)           , &
                   v_dt     (isw:iew,jsw:jew,:)           , &
                   t_dt     (isw:iew,jsw:jew,:)           , &
                   q_dt     (isw:iew,jsw:jew,:,1)         , &
                   q_dt     (isw:iew,jsw:jew,:,1:nt_prog) , &
                   Surf_diff                              , &
                   lprec    (isw:iew,jsw:jew)             , &
                   fprec    (isw:iew,jsw:jew)             , &
                   gust     (isw:iew,jsw:jew)             )
          else
              call physics_driver_up (isw-beglon+1, iew-beglon+1, jsw-beglat+1, jew-beglat+1, &
                   Time_prev, Time, Time_next             , &
                   rlat     (isw:iew,jsw:jew)             , &
                   rlon     (isw:iew,jsw:jew)             , &
                   area     (isw:iew,jsw:jew)             , &
              p_half(isw:iew,jsw:jew,:), p_full(isw:iew,jsw:jew,:), &
              z_half(isw:iew,jsw:jew,:), z_full(isw:iew,jsw:jew,:), &
                   omga     (isw:iew,jsw:jew,:)           , &
                   ua       (isw:iew,jsw:jew,:)           , &
                   va       (isw:iew,jsw:jew,:)           , &
                   pt       (isw:iew,jsw:jew,:)           , &
                   q        (isw:iew,jsw:jew,:,1)         , &
                   q        (isw:iew,jsw:jew,:,:)         , & ! cjg: pass all tracers
                   ua       (isw:iew,jsw:jew,:)           , &
                   va       (isw:iew,jsw:jew,:)           , &
                   pt       (isw:iew,jsw:jew,:)           , &
                   q        (isw:iew,jsw:jew,:,1)         , &
                   q        (isw:iew,jsw:jew,:,:)         , & ! cjg: pass all tracers
                   frac_land(isw:iew,jsw:jew)             , &
                   u_star   (isw:iew,jsw:jew)             , &
                   b_star   (isw:iew,jsw:jew)             , &
                   q_star   (isw:iew,jsw:jew)             , &
                   u_dt     (isw:iew,jsw:jew,:)           , &
                   v_dt     (isw:iew,jsw:jew,:)           , &
                   t_dt     (isw:iew,jsw:jew,:)           , &
                   q_dt     (isw:iew,jsw:jew,:,1)         , &
                   q_dt     (isw:iew,jsw:jew,:,1:nt_prog) , &
                   Surf_diff                              , &
                   lprec    (isw:iew,jsw:jew)             , &
                   fprec    (isw:iew,jsw:jew)             , &
                   gust     (isw:iew,jsw:jew)              )
          endif

    enddo
    call physics_driver_moist_end

    call physics_driver_up_endts(1, 1) !Note that these arguments are not used yet.

    call fv_array_sync()
    call mpp_clock_end(id_fv_physics_up)
!    call timing_on('update_fv')
    call mpp_clock_begin(id_update_fv_phys)
    call update_fv_phys( dt_phys, nt_prog, &
         .true., do_atmos_nudge, Time_next)
    call fv_array_sync()
    call mpp_clock_end(id_update_fv_phys)
!    call timing_off('update_fv')

  end subroutine fv_physics_up



subroutine fv_physics_end (Time)

!-----------------------------------------------------------------------
   type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
!  NOTE: this is not the dynamics time
!-----------------------------------------------------------------------


    call physics_driver_end (Time)
    call atmos_nudge_end

!-----------------------------------------------------------------------

end subroutine fv_physics_end



subroutine compute_p_z (istart, jstart, isiz, jsiz )
#include "fv_arrays.h"


  integer, intent(in):: istart, jstart, isiz, jsiz
!  real,    intent(out), dimension(isiz,jsiz,nlev)   :: p_full, z_full
!  real,    intent(out), dimension(isiz,jsiz,nlev+1) :: p_half, z_half

! local
  integer i,j,k,id,jd
  real    tvm
#include "fv_point.inc"

!----------------------------------------------------
! Compute pressure and height at full and half levels
!----------------------------------------------------

  do j = 1,jsiz
     jd = j + jstart - 1
     do i=1,isiz
        id = i + istart - 1
!        z_half(i,j,nlev+1) = phis(id,jd) * ginv
        z_half(id,jd,nlev+1) = phis(id,jd) * ginv
     enddo
  end do

  do k=1,nlev+1
     do j = 1,jsiz
        jd = j + jstart - 1
        do i=1,isiz
           id = i + istart - 1
!           p_half(i,j,k) = pe(id,k,jd)
           p_half(id,jd,k) = pe(id,k,jd)
        enddo
     enddo
  end do

  if ( use_tendency ) then
      do k=nlev,1,-1
         do j = 1,jsiz
            jd = j + jstart - 1
            do i=1,isiz
               id = i + istart - 1
               tvm = rrg*t_phys(id,jd,k)*(1.+zvir*q_phys(id,jd,k,1))
! Testing for M30_v1
!                    tvm = rrg*pt(id,jd,k)*(1.+zvir*q(id,jd,k,1))
!               p_full(i,j,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
!               z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
!               z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(id,k+1,jd)-peln(id,k,jd))
               p_full(id,jd,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
               z_full(id,jd,k) = z_half(id,jd,k+1) + tvm*(1.-p_half(id,jd,k)/p_full(id,jd,k))
               z_half(id,jd,k) = z_half(id,jd,k+1) + tvm*(peln(id,k+1,jd)-peln(id,k,jd))
            enddo
         end do
      enddo
  else
      do k=nlev,1,-1
         do j = 1,jsiz
            jd = j + jstart - 1
            do i=1,isiz
               id = i + istart - 1
               tvm = rrg*pt(id,jd,k)*(1.+zvir*q(id,jd,k,1))
!               p_full(i,j,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
!               z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
!               z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(id,k+1,jd)-peln(id,k,jd))
               p_full(id,jd,k) = delp(id,jd,k)/(peln(id,k+1,jd)-peln(id,k,jd))
               z_full(id,jd,k) = z_half(id,jd,k+1) + tvm*(1.-p_half(id,jd,k)/p_full(id,jd,k))
               z_half(id,jd,k) = z_half(id,jd,k+1) + tvm*(peln(id,k+1,jd)-peln(id,k,jd))
            enddo
         end do
       enddo
     endif

end subroutine compute_p_z
subroutine compute_g_avg(rrv, tracer_name)
#include "fv_arrays.h"

  real,          intent(inout) :: rrv(nt_prog)
  character(len=*), intent(in) :: tracer_name

  real psfc_sum(beglon:endlon,beglat:endlat,1),  &
       qp_sum(beglon:endlon,beglat:endlat,1), &
       qp, s1, s2
  integer j, i, k, idx
#include "fv_point.inc"
  psfc_sum = 0.
  qp_sum = 0.
  idx = get_tracer_index(MODEL_ATMOS, trim(tracer_name))
  if(idx /= NO_TRACER) then
     
     ! Subroutine argument check      
!     if(.not. allocated(pe)) then
!        call error_mesg("compute_g_avg", "pressure array is not available", FATAL)
!     else
!        if(use_tendency) then
!           if(.not. allocated(q_phys)) call error_mesg("compute_g_avg", &
!                "q_phys is not available when use_tendency = .true.", &
!                FATAL)
!        else
!           if(.not. allocated(q)) call error_mesg("compute_g_avg", &
!                "q is not available when use_tendency = .false.", &
!                FATAL)
!        endif
!     endif
     
!---------------------------------------------------------------------
!  define pressure-weighted column mean value of dry mass mixing 
!  ratio  for tracer idx. assumption is that the tracer field q_phys
!  is a moist mass mixing ratio. convert to dry mass mixing ratio by 
!  dividing by (1 - qh2o).
!---------------------------------------------------------------------
     do j=beglat,endlat
        do i = beglon, endlon
           psfc_sum(i,j,1) = pe(i,nlev+1,j)*area(i,j)
           qp = 0.0
           do k = 2, nlev+1
              if(use_tendency) then
                 qp = qp + (q_phys(i,j,k-1,idx)/  &
                                (1.0 - q_phys(i,j,k-1,sphum))) * &
                                             (pe(i,k,j) - pe(i,k-1,j))
              else
                 qp = qp + (q(i,j,k-1, idx)/ &
                                (1.0 - q_phys(i,j,k-1,sphum))) * &
                                            (pe(i,k,j) - pe(i,k-1,j))
              endif
           enddo
           qp_sum(i,j,1) = qp * area(i, j)
        enddo
     enddo
     
!---------------------------------------------------------------------
!    compute global sum of pressure-weighted tracer dry mass mixing 
!    ratio and global pressure.
!---------------------------------------------------------------------
     s1 = mpp_global_sum(fv_domain, psfc_sum, flags=BITWISE_EXACT_SUM)
     s2 = mpp_global_sum(fv_domain, qp_sum, flags=BITWISE_EXACT_SUM)

!---------------------------------------------------------------------
!    produce global mean tracer dry mass mixing ratio.
!---------------------------------------------------------------------
     rrv(idx) = s2 / s1
     !       if(mpp_pe() == mpp_root_pe()) print *, 'Note from PE: ', mpp_pe(), tracer_name, ' dry MASS mixing ratio = ', rrv

!---------------------------------------------------------------------
!    convert the tracer dry mass mixing ratio to the dry volume 
!    mixing ratio.
!---------------------------------------------------------------------

     if (trim(tracer_name).eq.'co2') then
       rrv(idx) = rrv(idx)*WTMAIR/WTMCO2
     !       if(mpp_pe() == mpp_root_pe()) print *, 'Note from PE: ', mpp_pe(), tracer_name, ' dry VOLUME mixing ratio = ', rrv
     endif

  endif
  
end subroutine compute_g_avg

end module fv_physics_mod
