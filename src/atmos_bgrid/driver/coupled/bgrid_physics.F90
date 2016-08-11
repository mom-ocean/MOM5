
module bgrid_physics_mod

!-----------------------------------------------------------------------
!
!        interface for b-grid dynamics with atmospheric physics
!
!-----------------------------------------------------------------------

use bgrid_core_driver_mod, only: bgrid_dynam_type
use bgrid_halo_mod,        only: update_halo, NORTH, EAST, WEST, SOUTH, &
                                 TEMP, UWND, VWND
use bgrid_prog_var_mod,    only: prog_var_type
use bgrid_vert_mod,        only: vert_grid_type, compute_height, &
                                 compute_pres_full, compute_pres_half, &
				 compute_pres_depth
use bgrid_horiz_mod,       only: horiz_grid_type, get_horiz_grid_bound, TGRID

use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID
use bgrid_integrals_mod,   only: global_integral
use       bgrid_masks_mod, only: grid_mask_type

use      time_manager_mod, only: time_type, get_time, operator(-)
use               fms_mod, only: error_mesg, FATAL, write_version_number , mpp_pe
use    physics_driver_mod, only: physics_driver_init,  &
                                 physics_driver_end,   &
                                 physics_driver_moist_init, &
                                 physics_driver_moist_end, &
                                 physics_driver_down_time_vary, &
                                 physics_driver_up_time_vary,  &
                                 physics_driver_down_endts,  &
                                 physics_driver_up_endts, &
                                 physics_driver_down,  &
                                 physics_driver_up,    &
                                 surf_diff_type
use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index, NO_TRACER

!-----------------------------------------------------------------------

implicit none
private

public   bgrid_physics_down, bgrid_physics_up,  &
         bgrid_physics_init, bgrid_physics_end
public   surf_diff_type

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id: bgrid_physics.F90,v 19.0.2.1 2013/12/18 23:34:26 Niki.Zadeh Exp $'
character(len=128) :: tagname = '$Name:  $'
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:) :: u_dt, v_dt

contains

!#######################################################################

subroutine bgrid_physics_down (window, dt_phys,            &
                               Time_prev, Time, Time_next, &
                               Hgrid, Vgrid, Dynam,        &
                               Var, Var_dt,                &
                               frac_land,   albedo,        &
                               albedo_vis_dir, albedo_nir_dir,     &
                               albedo_vis_dif, albedo_nir_dif,     &
                               rough_vel,   t_surf,        &
                               u_star, b_star, q_star,     &
                               dtau_du, dtau_dv, tau_x, tau_y,     &
                               flux_sw,    &
                               flux_sw_dir,    &
                               flux_sw_dif,    &
                               flux_sw_down_vis_dir,  &
                               flux_sw_down_vis_dif,  &
                               flux_sw_down_total_dir,         &
                               flux_sw_down_total_dif,         &
                               flux_sw_vis, &
                               flux_sw_vis_dir, &
                               flux_sw_vis_dif, &
                               flux_lw, coszen,   &
                               gust, Surf_diff, frac_open_sea )

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
  integer, intent(in)                :: window(2)
  real,    intent(in)                :: dt_phys
       type(time_type),intent(in)    :: Time_prev, Time, Time_next
type (horiz_grid_type),intent(inout) :: Hgrid
type  (vert_grid_type),intent(in)    :: Vgrid
type(bgrid_dynam_type),intent(in)    :: Dynam
type   (prog_var_type),intent(inout) :: Var
type   (prog_var_type),intent(inout) :: Var_dt

   ! Note: Var is intent inout because of diagnostic tracers

  real, intent(in),    dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) :: &
                                           frac_land,  albedo,   &
                                     albedo_vis_dir, albedo_nir_dir, &
                                     albedo_vis_dif, albedo_nir_dif, &
                                           rough_vel,  t_surf,   &
                                           u_star,     b_star,   &
                                           q_star, dtau_du, dtau_dv, frac_open_sea
  real, intent(inout), dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) :: tau_x,   tau_y
  real, intent(out),   dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) ::  &
                                    flux_sw,        &
                                    flux_sw_dir,   &
                                    flux_sw_dif,   &
                                    flux_lw, coszen,    &
                                    flux_sw_down_vis_dir, &
                                    flux_sw_down_vis_dif, &
                                    flux_sw_down_total_dir, &
                                    flux_sw_down_total_dif, &
                                    flux_sw_vis, &
                                    flux_sw_vis_dir, &
                                    flux_sw_vis_dif, &
                                    gust

  type(surf_diff_type), intent(inout) :: Surf_diff
!-----------------------------------------------------------------------
  integer :: j, k, n, is, ie, js, je, i1, i2, j1, j2, sphum, nt, ntp
  integer :: ix, jx, idim, jdim
  real    :: dt 
  integer :: sec, day
!-----------------------------------------------------------------------

   real, dimension(window(1),window(2),Vgrid%nlev) :: p_full, z_full

   real, dimension(window(1),window(2),Vgrid%nlev+1) :: p_half, z_half
   real, dimension(1,1,Vgrid%nlev+1) :: phalfgrey

   real, dimension(Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: uh, vh

   real, dimension(window(1),window(2)) :: pssl_new, area
   real, dimension(size(Var%r, 4))      :: gavg_rrv
   real, dimension(1,1)                 :: psurf

!-----------------------------------------------------------------------
!---------------------------- do physics -------------------------------

    idim = window(1)
    jdim = window(2)

!   --- momentum and momentum tendency on mass grid ---
!      note: tendency is saved for "physics_up" call 

    call update_halo (Hgrid, UWND, Var_dt%u, halos=SOUTH+WEST)
    call update_halo (Hgrid, VWND, Var_dt%v, halos=SOUTH+WEST)

    if (Dynam%Masks%sigma) then
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID,    Var%u,    Var%v, uh  , vh  )
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, Var_dt%u, Var_dt%v, u_dt, v_dt)
    else
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, &
                          Var%u,    Var%v, uh  , vh  , mask_inp=Dynam%Masks%Vel%mask)
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, &
                       Var_dt%u, Var_dt%v, u_dt, v_dt, mask_inp=Dynam%Masks%Vel%mask)
    endif

!   --- loop through physics windows ---

    sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    if (sphum <= 0) call error_mesg ('bgrid_physics_mod', &
               'specific humidity tracer not found', FATAL)
    ntp = Var_dt%ntrace
    nt  = Var   %ntrace
    js = Hgrid%Tmp%js

    gavg_rrv = 0.
    call compute_g_avg(Hgrid, Vgrid, Var, Dynam%Masks, gavg_rrv, 'co2')

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


    do while ( js <= Hgrid%Tmp%je )

       je = min ( js+jdim-1, Hgrid%Tmp%je )
       jx = je-js+1
       is = Hgrid%Tmp%is

    do while ( is <= Hgrid%Tmp%ie )

       ie = min ( is+idim-1, Hgrid%Tmp%ie )
       ix = ie-is+1

!      ---- pass updated surface pressure ----
       do j = 1, jx
          pssl_new(1:ix,j) = Var   %pssl(is:ie,js+j-1) + &
                             Var_dt%pssl(is:ie,js+j-1) * dt_phys
          area    (1:ix,j) = Hgrid%Tmp%area(js+j-1)
       enddo

       call compute_pres_full (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_full(1:ix,1:jx,:))
       call compute_pres_half (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_half(1:ix,1:jx,:))

!      compute a reference profile of pressure based on psurf = 1000 hPa
       psurf = reshape ( (/ 100000. /), (/ 1, 1 /) )
       call compute_pres_half (Vgrid, psurf, phalfgrey)


!      ----- compute height (in meters) -------
!      ----- spec humidity assumed to be tracer #1 --------

       call compute_height (Vgrid, Dynam%fisl(is:ie,js:je),     &
                                   Var%t (is:ie,js:je,:),       &
                                   Var%r (is:ie,js:je,:,sphum), &
                                   p_full (1:ix,1:jx,:),        &
                                   p_half (1:ix,1:jx,:),        &
                                   z_full (1:ix,1:jx,:),        &
                                   z_half (1:ix,1:jx,:),        &
                                   Dynam%Masks%Tmp%mask(is:ie,js:je,:) )

!      ---- j-axis indices in the global physics grid ----

       j1 = js-Hgrid%Tmp%js+1; j2 = j1+(je-js)
       i1 = is-Hgrid%Tmp%is+1; i2 = i1+(ie-is)

!-----------------------------------------------------------------------
!-------------------------- call physics -------------------------------
!------------ (need to add leap-frog option for uh,vh) -----------------
!-----------------------------------------------------------------------
  if ( .not. Dynam%Masks%sigma ) then
!------------ eta coordinate -------------------------------------------
      call physics_driver_down (i1, i2, j1, j2                 ,&
                           Time_prev, Time, Time_next          ,&
         Hgrid%Tmp%aph (is:ie,js:je), Hgrid%Tmp%alm (is:ie,js:je),&
         area     ( 1:ix, 1:jx),    p_half   ( 1:ix, 1:jx,:)     ,&
         p_full   ( 1:ix, 1:jx,:),  z_half   ( 1:ix, 1:jx,:)     ,&
         z_full   ( 1:ix, 1:jx,:)                                ,&
         phalfgrey                                               ,&
         uh       (is:ie,js:je,:)                                ,&
         vh       (is:ie,js:je,:)     ,  Var%t(is:ie,js:je,:)         ,&
         Var%r(is:ie,js:je,:,sphum)   ,  Var%r(is:ie,js:je,:, :   )   ,&
         uh       (is:ie,js:je,:)     ,  vh       (is:ie,js:je,:)     ,&
         Var%t (is:ie,js:je,:)        ,  Var%r (is:ie,js:je,:,sphum)  ,&
         Var%r (is:ie,js:je,:, :   )  ,                                &
         frac_land(is:ie,js:je)       , rough_vel(is:ie,js:je)        ,&
         frac_open_sea(is:ie,js:je)   ,                                &
         albedo   (is:ie,js:je)       , albedo_vis_dir (is:ie,js:je)  ,&
        albedo_nir_dir (is:ie,js:je) , albedo_vis_dif (is:ie,js:je)   ,&
        albedo_nir_dif (is:ie,js:je) , t_surf   (is:ie,js:je)         ,&
         u_star   (is:ie,js:je)   , b_star   (is:ie,js:je)   ,&
        q_star   (is:ie,js:je)   , dtau_du  (is:ie,js:je),   dtau_dv (is:ie,js:je),&
        tau_x    (is:ie,js:je)   , tau_y    (is:ie,js:je)   ,&
         u_dt     (is:ie,js:je,:) , v_dt     (is:ie,js:je,:) ,&
         Var_dt%t (is:ie,js:je,:) , Var_dt%r (is:ie,js:je,:,sphum),&
         Var_dt%r (is:ie,js:je,:,1:ntp) , flux_sw  (is:ie,js:je)   ,&
                           flux_sw_dir  (is:ie,js:je)              ,&
                           flux_sw_dif  (is:ie,js:je)              ,&
                           flux_sw_down_vis_dir  (is:ie,js:je)     ,&
                           flux_sw_down_vis_dif  (is:ie,js:je)     ,&
                           flux_sw_down_total_dir  (is:ie,js:je)   ,&
                           flux_sw_down_total_dif  (is:ie,js:je)   ,&
                           flux_sw_vis  (is:ie,js:je)          ,&
                           flux_sw_vis_dir  (is:ie,js:je)          ,&
                           flux_sw_vis_dif  (is:ie,js:je)          ,&
                           flux_lw  (is:ie,js:je)              ,&
                           coszen   (is:ie,js:je)              ,&
                           gust     (is:ie,js:je)              ,&
                           Surf_diff                           ,&
                           gavg_rrv                            ,&
                      mask=Dynam%Masks%Tmp%mask(is:ie,js:je,:) ,&
                      kbot=Dynam%Masks%Tmp%kbot(is:ie,js:je)    )
  else
!------------- sigma coordinate ----------------------------------------
      call physics_driver_down (i1, i2, j1, j2                 ,&
                           Time_prev, Time, Time_next          ,&
         Hgrid%Tmp%aph (is:ie,js:je), Hgrid%Tmp%alm (is:ie,js:je),&
         area     ( 1:ix, 1:jx),    p_half   ( 1:ix, 1:jx,:)     ,&
         p_full   ( 1:ix, 1:jx,:),  z_half   ( 1:ix, 1:jx,:)     ,&
         z_full   ( 1:ix, 1:jx,:)                                ,&
         phalfgrey                                               ,&
         uh       (is:ie,js:je,:)                                ,&
         vh       (is:ie,js:je,:)   ,  Var%t(is:ie,js:je,:)         ,&
         Var%r(is:ie,js:je,:,sphum)  ,  Var%r(is:ie,js:je,:, :   )   ,&
         uh       (is:ie,js:je,:)   ,  vh       (is:ie,js:je,:)     ,&
         Var%t (is:ie,js:je,:) ,  Var%r (is:ie,js:je,:,sphum)  ,&
         Var%r (is:ie,js:je,:, :   ) ,  &
         frac_land(is:ie,js:je)   , rough_vel(is:ie,js:je)  ,&
         frac_open_sea(is:ie,js:je)   ,                                &
         albedo   (is:ie,js:je)       , albedo_vis_dir (is:ie,js:je) ,&
        albedo_nir_dir (is:ie,js:je) , albedo_vis_dif (is:ie,js:je) ,&
        albedo_nir_dif (is:ie,js:je) , t_surf   (is:ie,js:je)       ,&
         u_star   (is:ie,js:je)   , b_star   (is:ie,js:je)   ,&
        q_star   (is:ie,js:je)   , dtau_du  (is:ie,js:je)   ,dtau_dv (is:ie,js:je),&
        tau_x    (is:ie,js:je)   , tau_y    (is:ie,js:je)   ,&
         u_dt     (is:ie,js:je,:) , v_dt     (is:ie,js:je,:) ,&
         Var_dt%t (is:ie,js:je,:) , Var_dt%r (is:ie,js:je,:,sphum),&
         Var_dt%r (is:ie,js:je,:,1:ntp) , flux_sw  (is:ie,js:je)   ,&
!                          Time_prev, Time, Time_next          ,&
!                          Hgrid%Tmp%aph (is:ie,js:je)         ,&
!                          Hgrid%Tmp%alm (is:ie,js:je)         ,&
!                          area     ( 1:ix, 1:jx)              ,&
!                          p_half   ( 1:ix, 1:jx,:)            ,&
!                          p_full   ( 1:ix, 1:jx,:)            ,&
!                          z_half   ( 1:ix, 1:jx,:)            ,&
!                          z_full   ( 1:ix, 1:jx,:)            ,&
!                          uh       (is:ie,js:je,:)            ,&
!                          vh       (is:ie,js:je,:)            ,&
!                          Var%t(is:ie,js:je,:)                ,&
!                          Var%r(is:ie,js:je,:,sphum)          ,&
!                          Var%r(is:ie,js:je,:,1:ntp)          ,&
!                          uh       (is:ie,js:je,:)            ,&
!                          vh       (is:ie,js:je,:)            ,&
!                          Var%t (is:ie,js:je,:)               ,&
!                          Var%r (is:ie,js:je,:,sphum)         ,&
!                          Var%r (is:ie,js:je,:,1:ntp)         ,&
!                          Var%r (is:ie,js:je,:,ntp+1:nt)      ,&
!                          frac_land(is:ie,js:je)              ,&
!                          rough_vel(is:ie,js:je)              ,&
!                          albedo   (is:ie,js:je)              ,&
!                          albedo_vis_dir (is:ie,js:je)              ,&
!                          albedo_nir_dir (is:ie,js:je)              ,&
!                          albedo_vis_dif (is:ie,js:je)              ,&
!                          albedo_nir_dif (is:ie,js:je)              ,&
!                          t_surf   (is:ie,js:je)              ,&
!                          u_star   (is:ie,js:je)              ,&
!                          b_star   (is:ie,js:je)              ,&
!                          q_star   (is:ie,js:je)              ,&
!                          dtau_dv  (is:ie,js:je)              ,&
!                          tau_x    (is:ie,js:je)              ,&
!                          tau_y    (is:ie,js:je)              ,&
!                          u_dt     (is:ie,js:je,:)            ,&
!                          v_dt     (is:ie,js:je,:)            ,&
!                          Var_dt%t (is:ie,js:je,:)            ,&
!                          Var_dt%r (is:ie,js:je,:,sphum)      ,&
!                          Var_dt%r (is:ie,js:je,:,1:ntp)      ,&
!                          flux_sw  (is:ie,js:je)              ,&
                           flux_sw_dir  (is:ie,js:je)              ,&
                           flux_sw_dif  (is:ie,js:je)              ,&
                           flux_sw_down_vis_dir  (is:ie,js:je)     ,&
                           flux_sw_down_vis_dif  (is:ie,js:je)     ,&
                           flux_sw_down_total_dir  (is:ie,js:je)   ,&
                           flux_sw_down_total_dif  (is:ie,js:je)   ,&
                           flux_sw_vis  (is:ie,js:je)          ,&
                           flux_sw_vis_dir  (is:ie,js:je)          ,&
                           flux_sw_vis_dif  (is:ie,js:je)          ,&
                           flux_lw  (is:ie,js:je)              ,&
                           coszen   (is:ie,js:je)              ,&
                           gust     (is:ie,js:je)              ,&
                           Surf_diff                           ,&
                           gavg_rrv                            )
  endif

        is = is + idim

     enddo

        js = js + jdim

     enddo


  call physics_driver_down_endts (is-Hgrid%Tmp%is+1, js-Hgrid%Tmp%js+1)

! halo rows for tendencies do not need updating until after physics_up
! update halos for diagnostic tracers only

  if (ntp < nt) call update_halo (Hgrid, TEMP, Var%r(:,:,:,ntp+1:nt))

!-----------------------------------------------------------------------

end subroutine bgrid_physics_down

!#######################################################################

subroutine bgrid_physics_up (window, dt_phys,            &
                             Time_prev, Time, Time_next, &
                             Hgrid, Vgrid, Dynam,        &
                             Var, Var_dt, omega,         &
                             frac_land, Surf_diff, lprec, fprec, gust, &
                             u_star, b_star, q_star )

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
  integer, intent(in)                :: window(2)
  real,    intent(in)                :: dt_phys
       type(time_type),intent(in)    :: Time_prev, Time, Time_next
type (horiz_grid_type),intent(inout) :: Hgrid
type  (vert_grid_type),intent(in)    :: Vgrid
type(bgrid_dynam_type),intent(in)    :: Dynam
type   (prog_var_type),intent(inout)    :: Var
type   (prog_var_type),intent(inout) :: Var_dt

real, intent(in), dimension(Hgrid%ilb:Hgrid%iub, &
                            Hgrid%jlb:Hgrid%jub, &
                            Vgrid%nlev) :: omega

  real, intent(in),  dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) :: frac_land
type(surf_diff_type), intent(inout) :: Surf_diff
  real, intent(out), dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) :: lprec, fprec
  real, intent(inout), dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) :: gust
  real,    intent(in), dimension(Hgrid%Tmp%is:,Hgrid%Tmp%js:) :: u_star, b_star, q_star

!-----------------------------------------------------------------------
  integer :: j, k, n, is, ie, js, je, i1, i2, j1, j2, sphum, ntp, npz,nt
  integer :: ix, jx, idim, jdim
  integer :: sec, day
  real    :: dt
!-----------------------------------------------------------------------

   real, dimension(window(1),window(2),Vgrid%nlev) :: p_full, z_full

   real, dimension(window(1),window(2),Vgrid%nlev+1) :: p_half, z_half

   real, dimension(Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: uh, vh, uh_dt, vh_dt

   real, dimension(window(1),window(2)) :: pssl_new, area

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
    call get_time (Time_next-Time_prev, sec, day)
    dt = real(sec+day*86400)
 
    call physics_driver_up_time_vary (Time, Time_next, dt)
 

!-----------------------------------------------------------------------
!---------------------------- do physics -------------------------------

    idim = window(1)
    jdim = window(2)

!   --- momentum and previous momentum tendency on mass grid ---

    if (Dynam%Masks%sigma) then
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID,    Var%u,    Var%v, uh   , vh   )
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, Var_dt%u, Var_dt%v, uh_dt, vh_dt)
    else
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, &
                          Var%u,    Var%v, uh   , vh   , mask_inp=Dynam%Masks%Vel%mask)
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, &
                       Var_dt%u, Var_dt%v, uh_dt, vh_dt, mask_inp=Dynam%Masks%Vel%mask)
    endif

!   --- loop through physics windows ---

    sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    if (sphum <= 0) call error_mesg ('bgrid_physics_mod', &
               'specific humidity tracer not found', FATAL)
    ntp = Var_dt%ntrace
    js = Hgrid%Tmp%js
    npz = size(p_full,3)
    nt = Var%ntrace
    call physics_driver_moist_init (window(1), window(2),  npz, ntp, nt) 
    do while ( js <= Hgrid%Tmp%je )

       je = min ( js+jdim-1, Hgrid%Tmp%je )
       jx = je-js+1
       is = Hgrid%Tmp%is

    do while ( is <= Hgrid%Tmp%ie )

       ie = min ( is+idim-1, Hgrid%Tmp%ie )
       ix = ie-is+1

!      ---- pass updated surface pressure ----
       do j = 1, jx
          pssl_new(1:ix,j) = Var   %pssl(is:ie,js+j-1) + &
                             Var_dt%pssl(is:ie,js+j-1) * dt_phys
          area    (1:ix,j) = Hgrid%Tmp%area(js+j-1)
       enddo

       call compute_pres_full (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_full(1:ix,1:jx,:))
       call compute_pres_half (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_half(1:ix,1:jx,:))

       call compute_height (Vgrid, Dynam%fisl(is:ie,js:je),     &
                                   Var%t (is:ie,js:je,:),       &
                                   Var%r (is:ie,js:je,:,sphum), &
                                   p_full (1:ix,1:jx,:),        &
                                   p_half (1:ix,1:jx,:),        &
                                   z_full (1:ix,1:jx,:),        &
                                   z_half (1:ix,1:jx,:),        &
                                   Dynam%Masks%Tmp%mask(is:ie,js:je,:) )

!      ---- j-axis indices in the global physics grid ----

       j1 = js-Hgrid%Tmp%js+1; j2 = j1+(je-js)
       i1 = is-Hgrid%Tmp%is+1; i2 = i1+(ie-is)

!-----------------------------------------------------------------------
!-------------------------- call physics -------------------------------
!------------ (need to add leap-frog option for uh,vh) -----------------
!-----------------------------------------------------------------------
  if ( .not. Dynam%Masks%sigma ) then
!------------ eta coordinate -------------------------------------------
      call physics_driver_up (i1, i2, j1, j2                      ,&
                              Time_prev, Time, Time_next          ,&
                              Hgrid%Tmp%aph (is:ie,js:je)         ,&
                              Hgrid%Tmp%alm (is:ie,js:je)         ,&
                              area     ( 1:ix, 1:jx)              ,&
                              p_half   ( 1:ix, 1:jx,:)            ,&
                              p_full   ( 1:ix, 1:jx,:)            ,&
                              z_half   ( 1:ix, 1:jx,:)            ,&
                              z_full   ( 1:ix, 1:jx,:)            ,&
                              omega    (is:ie,js:je,:)            ,&
                              uh       (is:ie,js:je,:)            ,&
                              vh       (is:ie,js:je,:)            ,&
                              Var%t(is:ie,js:je,:)                ,&
                              Var%r(is:ie,js:je,:,sphum)          ,&
                              Var%r(is:ie,js:je,:,:)          ,&
                              uh       (is:ie,js:je,:)            ,&
                              vh       (is:ie,js:je,:)            ,&
                              Var%t (is:ie,js:je,:)               ,&
                              Var%r (is:ie,js:je,:,sphum)         ,&
                              Var%r (is:ie,js:je,:,:)         ,&
                              frac_land(is:ie,js:je)              ,&
                              u_star(is:ie,js:je)                 ,&
                              b_star(is:ie,js:je)                 ,&
                              q_star(is:ie,js:je)                 ,&
                              u_dt     (is:ie,js:je,:)            ,&
                              v_dt     (is:ie,js:je,:)            ,&
                              Var_dt%t (is:ie,js:je,:)            ,&
                              Var_dt%r (is:ie,js:je,:,sphum)      ,&
                              Var_dt%r (is:ie,js:je,:,:)      ,&
                              Surf_diff                           ,&
                              lprec    (is:ie,js:je)              ,&
                              fprec    (is:ie,js:je)              ,&
                              gust     (is:ie,js:je)              ,&
                         mask=Dynam%Masks%Tmp%mask(is:ie,js:je,:) ,&
                         kbot=Dynam%Masks%Tmp%kbot(is:ie,js:je)    )
  else
!------------- sigma coordinate ----------------------------------------
      call physics_driver_up (i1, i2, j1, j2                      ,&
                              Time_prev, Time, Time_next          ,&
                              Hgrid%Tmp%aph (is:ie,js:je)         ,&
                              Hgrid%Tmp%alm (is:ie,js:je)         ,&
                              area     ( 1:ix, 1:jx)              ,&
                              p_half   ( 1:ix, 1:jx,:)            ,&
                              p_full   ( 1:ix, 1:jx,:)            ,&
                              z_half   ( 1:ix, 1:jx,:)            ,&
                              z_full   ( 1:ix, 1:jx,:)            ,&
                              omega    (is:ie,js:je,:)            ,&
                              uh       (is:ie,js:je,:)            ,&
                              vh       (is:ie,js:je,:)            ,&
                              Var%t(is:ie,js:je,:)                ,&
                              Var%r(is:ie,js:je,:,sphum)          ,&
                              Var%r(is:ie,js:je,:,1:ntp)          ,&
                              uh       (is:ie,js:je,:)            ,&
                              vh       (is:ie,js:je,:)            ,&
                              Var%t (is:ie,js:je,:)               ,&
                              Var%r (is:ie,js:je,:,sphum)         ,&
                              Var%r (is:ie,js:je,:,1:ntp)         ,&
                              frac_land(is:ie,js:je)              ,&
                              u_star(is:ie,js:je)                 ,&
                              b_star(is:ie,js:je)                 ,&
                              q_star(is:ie,js:je)                 ,&
                              u_dt     (is:ie,js:je,:)            ,&
                              v_dt     (is:ie,js:je,:)            ,&
                              Var_dt%t (is:ie,js:je,:)            ,&
                              Var_dt%r (is:ie,js:je,:,sphum)      ,&
                              Var_dt%r (is:ie,js:je,:,1:ntp)      ,&
                              Surf_diff                           ,&
                              lprec    (is:ie,js:je)              ,&
                              fprec    (is:ie,js:je)              ,&
                              gust     (is:ie,js:je) )
  endif

        is = is + idim

     enddo

        js = js + jdim

     enddo

     call physics_driver_moist_end
     call physics_driver_up_endts (is-Hgrid%Tmp%is+1, js-Hgrid%Tmp%js+1)

!-----------------------------------------------------------------------
! compute momentum tendencies on mass grid for physics only
!    udt(phys) = udt(current) - udt(before physics)

     uh(:,:,:) = u_dt(:,:,:) - uh_dt(:,:,:)
     vh(:,:,:) = v_dt(:,:,:) - vh_dt(:,:,:)

!  update halos of momentum tendencies on mass grid
!  then move momentum tendencies to momentum grid

     call update_halo (Hgrid, TEMP, uh, halos=NORTH+EAST)
     call update_halo (Hgrid, TEMP, vh, halos=NORTH+EAST)

     call change_grid (Hgrid, TEMP_GRID, WIND_GRID, uh, vh, uh, vh)

     uh(:,Hgrid%jub,:) = 0.0   ! zero out unused polar halo row
     vh(:,Hgrid%jub,:) = 0.0   ! no harm done when not polar row

!---- update momentum tendencies ----

     Var_dt%u = Var_dt%u + uh * Dynam%Masks%Vel%mask
     Var_dt%v = Var_dt%v + vh * Dynam%Masks%Vel%mask

!---- update all halo rows ----

     call update_halo (Hgrid, TEMP, Var_dt%t)
     call update_halo (Hgrid, TEMP, Var_dt%r)
     call update_halo (Hgrid, UWND, Var_dt%u)
     call update_halo (Hgrid, VWND, Var_dt%v)

!-----------------------------------------------------------------------

end subroutine bgrid_physics_up

!#######################################################################

subroutine bgrid_physics_init (axes, Time, Hgrid, Vgrid, Dynam, &
                               Var, Surf_diff)

!-----------------------------------------------------------------------
!
!   axes      = array of axis indices for diagnostics (x,y,pf,ph)
!   Time      = current time (time_type)
!   Hgrid     = horizontal grid constants
!   Vgrid     = vertical grid constants
!   Dynam     = current state of the dynamical core
!   Var       = prognostic variables
!
!-----------------------------------------------------------------------
integer,               intent(in)    :: axes(4)
type (time_type),      intent(in)    :: Time
type (horiz_grid_type),intent(inout) :: Hgrid
type  (vert_grid_type),intent(in)    :: Vgrid
type(bgrid_dynam_type),intent(in)    :: Dynam
type   (prog_var_type),intent(inout) :: Var
type (surf_diff_type), intent(inout) :: Surf_diff
!-----------------------------------------------------------------------
   real, dimension(Hgrid%Tmp%is:Hgrid%Tmp%ie+1)   :: lonb
   real, dimension(Hgrid%Tmp%js:Hgrid%Tmp%je+1)   :: latb
   real, dimension(Hgrid%Tmp%is:Hgrid%Tmp%ie, &
                   Hgrid%Tmp%js:Hgrid%Tmp%je,Vgrid%nlev+1) :: phalf
   real, dimension(Vgrid%nlev+1,2) :: pref
   real, dimension(2,1,Vgrid%nlev) :: pref_full
   real, dimension(2,1)            :: pref_sl
   real, dimension(Hgrid%Tmp%is:Hgrid%Tmp%ie+1,Hgrid%Tmp%js:Hgrid%Tmp%je+1) :: lonb2d, latb2d
!-----------------------------------------------------------------------
   integer :: i, j, k, n, nt, is, ie, js, je, unit
!-----------------------------------------------------------------------

    nt = Var%ntrace

!----- write version to logfile --------

   call write_version_number(version, tagname)

!---------- get local grid box edges ---------

    call get_horiz_grid_bound (Hgrid, TGRID, lonb, latb)
    do j = Hgrid%Tmp%js,Hgrid%Tmp%je+1
       lonb2d(:,j) = lonb(:)
    end do
    do i = Hgrid%Tmp%is,Hgrid%Tmp%ie+1
       latb2d(i,:) = latb(:)
    end do

!---------- reference profile -----------

    pref_sl = reshape ( (/ 101325., 81060. /), (/ 2, 1 /) )
    call compute_pres_full (Vgrid, pref_sl, pref_full)
    pref(1:Vgrid%nlev,1) = pref_full(1,1,:)
    pref(1:Vgrid%nlev,2) = pref_full(2,1,:)
    pref(Vgrid%nlev+1,1) = pref_sl(1,1)
    pref(Vgrid%nlev+1,2) = pref_sl(2,1)

!------- pressure at model layer interfaces -----

    is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
    js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je

    call compute_pres_half (Vgrid, Var%pssl(is:ie,js:je), phalf)

!---------- initialize physics -------

    if (Dynam%Masks%sigma) then
        call physics_driver_init (Time, lonb2d, latb2d, lonb2d, latb2d, axes, pref,  &
                                  Var%r(is:ie,js:je,:,1:nt),     &
                                  Surf_diff, phalf               )
    else
        call physics_driver_init (Time, lonb2d, latb2d, lonb2d, latb2d, axes, pref,   &
                                  Var%r(is:ie,js:je,:,1:nt),      &
                                  Surf_diff, phalf,               &
                        mask=Dynam%Masks%Tmp%mask(is:ie,js:je,:), &
                        kbot=Dynam%Masks%Vel%kbot(is:ie,js:je)    )
    endif

!       ---- boundaries for tracers ----
    call update_halo (Hgrid, TEMP, Var%r)

!       ---- storage for global data ----

    allocate (u_dt(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,Vgrid%nlev))
    allocate (v_dt(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,Vgrid%nlev))

!-----------------------------------------------------------------------

end subroutine bgrid_physics_init

!#######################################################################

subroutine bgrid_physics_end (Time)

!-----------------------------------------------------------------------
   type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
!  NOTE: this is not the dynamics time
!-----------------------------------------------------------------------

    call physics_driver_end (Time)

!-----------------------------------------------------------------------

end subroutine bgrid_physics_end

!#######################################################################

subroutine compute_g_avg(Hgrid, Vgrid, Var, Masks, rrv, tracer_name)

  type (horiz_grid_type),intent(inout) :: Hgrid
  type  (vert_grid_type),intent(in)    :: Vgrid
  type   (prog_var_type),intent(in)    :: Var
  type (grid_mask_type), intent(in)    :: Masks

  real, dimension(:),    intent(inout) :: rrv
  character(len=*), intent(in) :: tracer_name

  real, dimension (Hgrid%ilb:Hgrid%iub,          &
                   Hgrid%jlb:Hgrid%jub, Var%nlev) :: dpde, avg

  real avgps
  real psfc_sum, qp_sum, qp
  integer j, i, k, idx

    psfc_sum = 0.
    qp_sum = 0.
    idx = get_tracer_index(MODEL_ATMOS, trim(tracer_name))
    if(idx /= NO_TRACER) then
      call compute_pres_depth (Vgrid, Var%pssl, dpde)
      avg = Var%r(:,:,:,idx) * dpde
      avgps  = global_integral (Hgrid, 1, Var%ps,  do_exact=.true.)
      rrv(idx) = global_integral(Hgrid, 1, avg, Masks, .true.) / avgps
    endif

end subroutine compute_g_avg

end module bgrid_physics_mod

