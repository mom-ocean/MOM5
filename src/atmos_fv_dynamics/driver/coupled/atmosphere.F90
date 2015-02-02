module atmosphere_mod

!-----------------------------------------------------------------------
!
!         interface for fv dynamical core and physics
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------


use time_manager_mod, only: time_type, get_time, set_time, operator(+)

use fms_mod, only: file_exist, open_namelist_file, &
                               error_mesg, FATAL,              &
                               check_nml_error, stdlog,        &
                               write_version_number,           &
                               mpp_pe, mpp_root_pe,            &
                               close_file, set_domain,         &
                               mpp_clock_id, mpp_clock_begin,  &
                               mpp_clock_end, CLOCK_SUBCOMPONENT, &
                               clock_flag_default

use     fv_physics_mod, only: fv_physics_down, fv_physics_up,  &
                              fv_physics_init, fv_physics_end, &
                              surf_diff_type

use    mpp_domains_mod, only: domain2d
use  field_manager_mod, only: MODEL_ATMOS
use   diag_manager_mod, only: diag_send_complete
use      constants_mod, only: omega, cp_air, rdgas, grav, rvgas, kappa, radius, pstd_mks

use            fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                              endlon, rlonb, rlatb,  cold_start, ncnst, &
                              pnats, consv_te, ptop, fv_init, fv_domain, &
                              fv_end, change_time, restart_format, area, &
                              ak, bk
use     fv_diagnostics, only: fv_diag_init, fv_diag, fv_time
use       timingModule, only: timing_on, timing_off
use fv_restart_mod, only: fv_restart, write_fv_rst
use fv_dynamics_mod, only: fv_dynamics
use fv_arrays_mod, only: fv_print_chksums
use mpp_mod, only: mpp_error, input_nml_file
use tracer_manager_mod, only: get_tracer_index, NO_TRACER
use xgrid_mod, only: grid_box_type

!-----------------------------------------------------------------------

implicit none
private

public  atmosphere_down,       atmosphere_up,       &
        atmosphere_init,       atmosphere_end,      &
        atmosphere_resolution, atmosphere_boundary, &
        atmosphere_cell_area,  atmosphere_restart,  &
        get_atmosphere_axes,   atmosphere_domain

public  get_bottom_mass,  get_bottom_wind,  get_stock_pe

public  surf_diff_type

integer sec
integer seconds, days
integer liq_wat, ice_wat
!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 20.0 2013/12/13 23:08:15 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,0]

   integer, dimension(2) :: physics_window = (/0,0/)

! Note: the default size of window(1) is chosen to work for N30, N45, ...

   namelist /atmosphere_nml/ physics_window

!-----------------------------------------------------------------------
!---- private data ----

  type    (time_type) :: Time_step_atmos
  real                :: dt_atmos
  integer, dimension(4)              :: atmos_axes
  integer :: id_dynam, id_phys_down, id_phys_up, id_fv_diag

!-----------------------------------------------------------------------

contains

 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff, Grid_box)

 type (time_type),     intent(in)    :: Time_init, Time, Time_step
 type(surf_diff_type), intent(inout) :: Surf_diff
 type(grid_box_type),  intent(inout) :: Grid_box

  integer :: unit, ierr, io
  integer :: ss, ds, log_unit

!----- read namelist -----

#ifdef INTERNAL_FILE_NML
      read (input_nml_file,atmosphere_nml,iostat=io)
      ierr = check_nml_error (io, 'atmosphere_nml')
#else
    if (file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read (unit, nml=atmosphere_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'atmosphere_nml')
        enddo
 10     call close_file (unit)
    endif
#endif

!----- write version and namelist to log file -----

    call write_version_number ( version, tagname )
    if ( mpp_pe() == mpp_root_pe() ) then
       log_unit = stdlog()
       write (log_unit, nml=atmosphere_nml)
    endif

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----

   call fv_init( sec )
   call fv_restart( days, seconds )

    if ( cold_start .or. change_time ) then
        fv_time = time
    else
        fv_time = set_time (seconds, days)
        call get_time (Time, ss,  ds)

        if( seconds /= ss .or. days /= ds )   call  error_mesg         &
            ('FV_init:','Time inconsistent between fv_rst and INPUT/atmos_model.res', FATAL)
    endif

!----- initialize atmos_axes and fv_dynamics diagnostics

    call fv_diag_init( atmos_axes, Time )
!    call fv_print_chksums( 'after fv_diag_init' )
   

!----- initialize physics interface -----
!----- initialize domains for reading global physics data -----

    call set_domain ( fv_domain )

    call fv_physics_init (atmos_axes, Time, physics_window, Surf_diff)
!    call fv_print_chksums( 'after fv_physics_init' )


!  --- initialize clocks for dynamics, physics_down and physics_up

    id_dynam     = mpp_clock_id ('FV dynamical core',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    id_phys_down = mpp_clock_id ('Physics_down',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    id_phys_up   = mpp_clock_id ('Physics_up',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    id_fv_diag   = mpp_clock_id ('FV Diag',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

    call fv_print_chksums( 'Exiting  atmosphere_init' )

    liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
    ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
 end subroutine atmosphere_init

 subroutine atmosphere_down (Time,    frac_land,             &
      t_surf,  albedo,                &
      albedo_vis_dir, albedo_nir_dir, &
      albedo_vis_dif, albedo_nir_dif, &
      rough_mom,                      &
      u_star,  b_star, q_star,        &
      dtau_du, dtau_dv, tau_x, tau_y, &
      frac_open_sea,                  &
      gust, coszen, flux_sw,          &
      flux_sw_dir, flux_sw_dif,       &
      flux_sw_down_vis_dir,           &
      flux_sw_down_vis_dif,           &
      flux_sw_down_total_dir,         &
      flux_sw_down_total_dif,         &
      flux_sw_vis,                    &
      flux_sw_vis_dir,                &
      flux_sw_vis_dif,                &
      flux_lw,                        &
      Surf_diff                       )
#include "fv_arrays.h"
!
!        Time = time at the current time level
!

   type(time_type),intent(in)    :: Time

   real, intent(in),    dimension(:,:) :: frac_land,                 &
        t_surf, albedo,            &
        albedo_vis_dir, albedo_nir_dir, &
        albedo_vis_dif, albedo_nir_dif, &
        rough_mom, u_star, b_star, &
        q_star, dtau_du, dtau_dv
   real, intent(in),    dimension(:,:) :: frac_open_sea
   real, intent(inout), dimension(:,:) :: tau_x,  tau_y
   real, intent(out),   dimension(:,:) :: gust, coszen, flux_sw, &
        flux_sw_dir, flux_sw_dif, &
        flux_sw_down_vis_dir,     &
        flux_sw_down_total_dir,   &
        flux_sw_down_vis_dif,     &
        flux_sw_down_total_dif,   &
        flux_sw_vis,              &
        flux_sw_vis_dir,          &
        flux_sw_vis_dif, flux_lw
   type(surf_diff_type), intent(inout) :: Surf_diff

   type(time_type) :: Time_prev, Time_next
   real zvir
#include "fv_point.inc"

   zvir = rvgas/rdgas - 1.
   call fv_print_chksums( 'Entering  atmosphere_down' )

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- dynamics -----

   call timing_on('fv_dynamics')
   call mpp_clock_begin (id_dynam)
#ifndef USE_LIMA
   call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
        ncnst,   pnats,  .false., consv_te,            &
        u,       v,      delp,    pt,       q,         &
        ps,      pe,     pk,      pkz,      phis,      &
        omga,    peln,   ptop,    omega,    sec,       &
        zvir,    cp_air, rdgas,   kappa,  radius, ua, va, Time_next )
#else
   call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
        ncnst,   pnats,  .false., consv_te,            &
        u,       v,      delp,    pt,       q,         &
        ps,      pe,     pk,      pkz,      phis,      &
        omga,    peln,   ptop,    omega,    sec,       &
        zvir,    cp_air, rdgas,   kappa,  radius, ua, va )
#endif
   call mpp_clock_end (id_dynam)
   call timing_off('fv_dynamics')

!---- call physics -----

   call timing_on('fv_physics_down')
   call mpp_clock_begin (id_phys_down)
   call fv_physics_down (dt_atmos, Time_prev, Time, Time_next,     &
        frac_land, albedo,              &
        albedo_vis_dir, albedo_nir_dir, &
        albedo_vis_dif, albedo_nir_dif, &
        rough_mom,  t_surf,             &
        u_star,  b_star, q_star,        &
        dtau_du, dtau_dv, tau_x, tau_y, &
        flux_sw, flux_sw_dir,           &
        flux_sw_dif,                    &
        flux_sw_down_vis_dir,           &
        flux_sw_down_vis_dif,           &
        flux_sw_down_total_dir,         &
        flux_sw_down_total_dif,         &
        flux_sw_vis, flux_sw_vis_dir,   &
        flux_sw_vis_dif, flux_lw,       &
        coszen, gust, Surf_diff, frac_open_sea )
   call mpp_clock_end (id_phys_down)
   call timing_off('fv_physics_down')

   call fv_print_chksums( 'Exiting  atmosphere_down' )
 end subroutine atmosphere_down


 subroutine atmosphere_up (Time,  frac_land, Surf_diff, lprec, fprec, gust, u_star, b_star, q_star )

   type(time_type),intent(in)        :: Time
   real, intent(in),  dimension(:,:) :: frac_land
   type(surf_diff_type), intent(inout) :: Surf_diff
   real, intent(out), dimension(:,:) :: lprec,   fprec
   real, intent(inout), dimension(:,:) :: gust
   real,dimension(:,:),    intent(in) :: u_star, b_star, q_star

   type(time_type) :: Time_prev, Time_next
   real zvir

   zvir = rvgas/rdgas - 1.
   call fv_print_chksums( 'Entering  atmosphere_up' )
   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

   call timing_on('fv_physics_up')
   call mpp_clock_begin (id_phys_up)
   call fv_physics_up( dt_atmos, Time_prev, Time, Time_next,      &
        frac_land,  Surf_diff,           &
        lprec,      fprec,      gust, u_star, b_star, q_star     )
   call mpp_clock_end (id_phys_up)
   call timing_off('fv_physics_up')

!---- diagnostics for FV dynamics -----
   call timing_on('FV_DIAG')
   call mpp_clock_begin(id_fv_diag)
   fv_time = Time + Time_step_atmos
   call get_time (fv_time, seconds,  days)

   call fv_diag(fv_time, nlon, mlat, nlev, beglat, endlat, &       
        ncnst, zvir, dt_atmos, .false.)                          
   call timing_off('FV_DIAG')     
   call mpp_clock_end(id_fv_diag)

   ! Indicate to diag_manager to write diagnostics to file (if needed)
   ! This is needed for a threaded run.
   call diag_send_complete(Time_step_atmos)

   call fv_print_chksums( 'Exiting  atmosphere_up' )
#ifdef PSET_DEBUG
   call mpp_error( FATAL, 'ATMOSPHERE_UP is set to abort if -DPSET_DEBUG!' )
#endif
 end subroutine atmosphere_up


 subroutine atmosphere_end (Time, Grid_box)

 type (time_type), intent(in) :: Time
 type(grid_box_type),  intent(inout) :: Grid_box

!----- initialize domains for writing global physics data -----

    call set_domain ( fv_domain )
    call get_time (Time, seconds,  days)
    call write_fv_rst( 'RESTART/fv_rst.res', days, seconds, grav, &
         restart_format )

    call fv_end(days, seconds)      

    call fv_physics_end(Time) 


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

!---------------------------------------------------------------
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------

 subroutine atmosphere_resolution (i_size, j_size, global)

   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .TRUE.
   if( PRESENT(global) )local = .NOT.global

   if( local )then
       i_size = endlon - beglon + 1
       j_size = endlat - beglat + 1
   else
       i_size = nlon
       j_size = mlat
   end if
 end subroutine atmosphere_resolution


 subroutine atmosphere_cell_area  (area_out)
    real, dimension(:,:),  intent(out) :: area_out

    area_out(1:size(area_out,1), 1:size(area_out,2)) =  &
                                   area (beglon:endlon, beglat:endlat)

 end subroutine atmosphere_cell_area


!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------

 subroutine atmosphere_boundary (blon, blat, global)

    real,    intent(out)          :: blon(:,:), blat(:,:)   ! radian
    logical, intent(in), optional :: global
! Local:
    integer i,j
    logical :: local

    local = .TRUE.
    if( PRESENT(global) )local = .NOT.global

    if( local )then
        do i = beglon,endlon+1
           blon(i-beglon+1,:) = rlonb(i)
        end do
        do j = beglat,endlat+1
           blat(:,j-beglat+1) = rlatb(j)
        end do
    else
        do i=1,nlon+1
           blon(i,:) = rlonb(i)
        end do
        do j=1,mlat+1
           blat(:,j) = rlatb(j)
        end do
    end if

 end subroutine atmosphere_boundary




 subroutine atmosphere_domain (Domain)
 type(domain2d), intent(inout) :: Domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   Domain = fv_domain

 end subroutine atmosphere_domain


 subroutine get_atmosphere_axes ( axes )
! returns the axis indices associated with the coupling grid

   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----

     if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                           'get_atmosphere_axes in atmosphere_mod', &
                           'size of argument is incorrect', FATAL   )

     axes (1:size(axes(:))) = atmos_axes (1:size(axes(:)))

 
 end subroutine get_atmosphere_axes



 subroutine get_bottom_mass (t_bot, tr_bot, p_bot, z_bot, p_surf, slp)

#include "fv_arrays.h"
! returns temp, sphum, pres, height at the lowest model level
!         and surface pressure and sea level pressure

   real, intent(out), dimension(beglon:endlon,beglat:endlat)  &
        :: t_bot, p_bot, z_bot, p_surf, slp
   real, intent(out), dimension(beglon:endlon,beglat:endlat,ncnst-pnats):: tr_bot
   integer :: i, j, k, kr
   real zvir, rrg, sigtop, sigbot
   real, dimension(beglon:endlon,beglat:endlat) :: tref
   real, parameter :: tlaps = 6.5e-3
#include "fv_point.inc"

   rrg  = rdgas / grav
   zvir = rvgas/rdgas - 1.


   ! determine 0.8 sigma reference level
   sigtop = ak(1)/pstd_mks+bk(1)
   do k = 1, nlev 
      sigbot = ak(k+1)/pstd_mks+bk(k+1)
      if (sigbot+sigtop > 1.6) then
         kr = k
         exit
      endif   
      sigtop = sigbot
   enddo

!$omp parallel do default(shared)    &
!$omp private (i, j)
     do j = beglat, endlat
        do i = beglon, endlon
           p_surf(i,j) =  ps(i,j)
           t_bot(i,j) =  pt(i,j,nlev)

           p_bot(i,j) = delp(i,j,nlev)/(peln(i,nlev+1,j)-peln(i,nlev,j))
           z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*q(i,j,nlev,1))*  &
                  (1. - pe(i,nlev,j)/p_bot(i,j))
           ! sea level pressure
           tref(i,j) = pt(i,j,kr)*(delp(i,j,kr)/((peln(i,kr+1,j)-peln(i,kr,j))*ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = ps(i,j)*(1.+tlaps*phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
! Copy tracers
!$omp parallel do default(shared)    &
!$omp private (i, j, k)
     do k = 1,ncnst-pnats
        do j = beglat,endlat
           do i = beglon,endlon
              tr_bot(i,j,k) = q(i,j,nlev,k)
           enddo
        enddo
     enddo

 end subroutine get_bottom_mass



 subroutine get_bottom_wind (u_bot, v_bot)
#include "fv_arrays.h"
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------

   real, intent(out), dimension(beglon:,beglat:) :: u_bot, v_bot

   integer i, j
#include "fv_point.inc"
!Balaji: this cannot work unless beglon=1: corrected declaration Lbounds
   do j=beglat,endlat
      do i=beglon,endlon
         u_bot(i,j) = u_srf(i,j)
         v_bot(i,j) = v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind

 subroutine get_stock_pe(index, value)
#include "fv_arrays.h"
    integer, intent(in) :: index
    real, intent(out)   :: value
#ifdef USE_STOCK
    include 'stock.inc' 
#endif
    real wm(beglon:endlon, beglat:endlat)
    integer i,j,k
#include "fv_point.inc"
   
    select case (index)
#ifdef USE_STOCK
    case (ISTOCK_WATER)
#else
    case (1)
#endif
     
!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j = beglat, endlat
        do k=1,nlev
           do i = beglon, endlon
! Note: There is a check in fv_pack.F90 that ensures that tracer number one
!       is sphum and that the cloud water and ice tracers exist.
              wm(i,j) = wm(i,j) + delp(i,j,k)*(q(i,j,k,1)+q(i,j,k,liq_wat)+q(i,j,k,ice_wat))
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j = beglat, endlat
        do i = beglon, endlon
           value = value + wm(i,j)*area(i,j)
        enddo
     enddo
     value = value/grav

    case default
     value = 0.0
    end select

 end subroutine get_stock_pe 

end module atmosphere_mod
