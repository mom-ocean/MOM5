!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! A SEA ICE MODEL for coupling through the GFDL exchange grid;  this module    !
! manages fluxes, diagnostics, and ice timesteps; sea ice dynamics and         !
! thermodynamics are performed in ice_[dyn|thm].f90 - Mike Winton (Michael.Winton)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_model_mod

  use mpp_mod,          only: mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod,  only: mpp_update_domains, BGRID_NE, CGRID_NE
  use fms_mod,          only: error_mesg
  use diag_manager_mod, only: send_data
  use time_manager_mod, only: time_type, operator(+), get_date, get_time
  use time_manager_mod, only: operator(-), set_date
  use astronomy_mod,    only: universal_time, orbital_time, diurnal_solar, daily_mean_solar
  use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
  use constants_mod,    only: hlv, hlf, Tfreeze, grav, STEFAN
  use ocean_albedo_mod, only: compute_ocean_albedo            ! ice sets ocean surface
  use ocean_rough_mod,  only: compute_ocean_roughness         ! properties over water
  use ice_type_mod,     only: ice_data_type, ice_model_init, ice_model_end, hlim,&
                              mom_rough_ice, heat_rough_ice, atmos_winds, kmelt, &
                              slab_ice, spec_ice, ice_bulk_salin, id_cn, id_hi,  &
                              id_hs, id_t1, id_t2, id_ts, id_sh, id_lh, id_sw,   &
                              id_swdn, id_lw, id_lwdn, id_snofl, id_rain,        &
                              id_runoff_hflx, id_calving_hflx,                   &
                              id_runoff, id_calving, id_evap, id_saltf, id_tmelt,&
                              id_mi, id_bmelt, id_bheat, id_frazil, id_alb,      &
                              id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna,      &
                              id_sigi, id_sigii, id_stren, id_ui, id_vi, id_fax, &
                              id_fay, id_fix, id_fiy, id_fcx, id_fcy, id_fwx,    &
                              id_fwy, id_sn2ic,  id_ext, id_slp, id_sst, id_sss, &
                              id_ssh, id_uo, id_vo, id_e2m, id_qflim, id_qfres,  &
                              id_ix_trans, id_iy_trans,                          &
                              do_ice_restore, do_ice_limit, max_ice_limit,       &
                              ice_restore_timescale, do_init, h2o, heat,         &
                              conservation_check, slp2ocean, iceClock, verbose,  &
                              iceClock1, iceClock2, iceClock3, &
                              iceClock4, iceClock5, iceClock6, &
                              iceClock7, iceClock8, iceClock9, &
                              iceClocka, iceClockb, iceClockc, &
                              id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir,    &
                              id_sw_vis_dif, id_sw_nir_dir, id_sw_nir_dif,       &
                              cm2_bugs, ice_stock_pe, do_icebergs, ice_model_restart, &
                              add_diurnal_sw, id_mib, ice_data_type_chksum,      &
                              id_ustar, id_vstar, channel_viscosity, smag_ocn,   &
                              ssh_gravity, chan_cfl_limit, id_vocean, id_uocean, &
                              id_vchan, id_uchan, id_wnd
  use ice_type_mod,     only: do_sun_angle_for_alb,              &
			      id_alb_vis_dir, id_alb_vis_dif,    &
	                      id_alb_nir_dir, id_alb_nir_dif

  use ice_grid_mod,     only: uv_to_t, t_to_uv, vel_t_to_uv, cut_check, tripolar_grid
  use ice_grid_mod,     only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  use ice_grid_mod,     only: ice_advect, ice_avg, all_avg, ice_line, slab_ice_advect
  use ice_grid_mod,     only: geo_lon, geo_lat, cell_area, sin_rot, cos_rot, latitude
  use ice_spec_mod,     only: get_sea_surface
  use ice_grid_mod,     only: dte, dtn, dxv, dyv, dxt, dyt, dt_adv, wett
  use ice_grid_mod,     only: reproduce_siena_201303
  !
  ! the following two modules are the work horses of the sea ice model
  !
  use ice_thm_mod,      only: ice_optics, ice_thm_param, ice3lay_temp, ice3lay_resize
  use ice_thm_mod,      only: thm_pack, thm_unpack, DI, DS, MU_TS, TFI, e_to_melt
  use ice_dyn_mod,      only: ice_dynamics, ice_dyn_param, strain_angle, ice_strength, sigI, sigII
  use ice_bergs,        only: icebergs_run, icebergs_incr_mass

  implicit none
  private

  public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
  public :: ice_model_init, ice_model_end, ice_bottom_to_ice_top, update_ice_model_fast, ice_stock_pe, cell_area
  public :: update_ice_model_slow, update_ice_model_slow_up, update_ice_model_slow_dn ! for new coupler
  public :: ice_model_restart  ! for intermediate restart
  public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum, &
            lnd_ice_bnd_type_chksum, ice_data_type_chksum

  interface update_ice_model_fast ! overload to support old interface
     module procedure update_ice_model_fast_new
     !  module procedure update_ice_model_fast_old !Balaji: no need to expose old interface
  end interface
  !
  ! the following three types are for data exchange with the new coupler
  ! they are defined here but declared in coupler_main and allocated in flux_init
  !
  type :: ocean_ice_boundary_type
     real, dimension(:,:),   pointer :: u         =>NULL()
     real, dimension(:,:),   pointer :: v         =>NULL()
     real, dimension(:,:),   pointer :: t         =>NULL()
     real, dimension(:,:),   pointer :: s         =>NULL()
     real, dimension(:,:),   pointer :: frazil    =>NULL()
     real, dimension(:,:),   pointer :: sea_level =>NULL()
     real, dimension(:,:,:), pointer :: data      =>NULL() ! collective field for "named" fields above
     integer                         :: xtype              ! REGRID, REDIST or DIRECT used by coupler
     type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
  end type 

  type :: atmos_ice_boundary_type 
     real, dimension(:,:,:), pointer :: u_flux  =>NULL()
     real, dimension(:,:,:), pointer :: v_flux  =>NULL()
     real, dimension(:,:,:), pointer :: u_star  =>NULL()
     real, dimension(:,:,:), pointer :: t_flux  =>NULL()
     real, dimension(:,:,:), pointer :: q_flux  =>NULL()
     real, dimension(:,:,:), pointer :: lw_flux =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_vis_dir =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_vis_dif =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_nir_dir =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_nir_dif =>NULL()
     real, dimension(:,:,:), pointer :: lprec   =>NULL()
     real, dimension(:,:,:), pointer :: fprec   =>NULL()
     real, dimension(:,:,:), pointer :: dhdt    =>NULL()
     real, dimension(:,:,:), pointer :: dedt    =>NULL()
     real, dimension(:,:,:), pointer :: drdt    =>NULL()
     real, dimension(:,:,:), pointer :: coszen  =>NULL()
     real, dimension(:,:,:), pointer :: p       =>NULL()
     real, dimension(:,:,:), pointer :: data    =>NULL()
     integer                         :: xtype
     type(coupler_3d_bc_type)        :: fluxes     ! array of fluxes used for additional tracers
  end type

  type :: land_ice_boundary_type
     real, dimension(:,:),   pointer :: runoff  =>NULL()
     real, dimension(:,:),   pointer :: calving =>NULL()
     real, dimension(:,:),   pointer :: runoff_hflx  =>NULL()
     real, dimension(:,:),   pointer :: calving_hflx =>NULL()
     real, dimension(:,:,:), pointer :: data    =>NULL() ! collective field for "named" fields above
     integer                         :: xtype            ! REGRID, REDIST or DIRECT used by coupler
  end type

  real, parameter :: T_sw_freeze = Tfreeze-1.8 ! seawater freezing temperature (K)
  logical :: first_time = .true.               ! first time ice_bottom_to_ice_top

contains

  !#######################################################################
  !
  ! new coupler interface to provide ocean surface data to atmosphere
  !
  subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
    type (ice_data_type),          intent(inout) :: Ice

    call mpp_clock_begin(iceClock)
    call mpp_clock_begin(iceClock1)
    call ice_bottom_to_ice_top (Ice, Ocean_boundary%t, Ocean_boundary%u, Ocean_boundary%v,        &
                                Ocean_boundary%frazil, Ocean_boundary, Ocean_boundary%s, Ocean_boundary%sea_level  )
    call mpp_clock_end(iceClock1)
    call mpp_clock_end(iceClock)

  end subroutine update_ice_model_slow_up
  !
  ! new coupler interface to do slow ice processes:  dynamics, transport, mass
  !
  subroutine update_ice_model_slow_dn ( Atmos_boundary, Land_boundary, Ice )
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
    type(land_ice_boundary_type),  intent(inout) :: Land_boundary
    type (ice_data_type),          intent(inout) :: Ice

    call mpp_clock_begin(iceClock)
    call mpp_clock_begin(iceClock2)
    call update_ice_model_slow (Ice, Land_boundary%runoff, Land_boundary%calving, Land_boundary%runoff_hflx, Land_boundary%calving_hflx, Atmos_boundary%p )
    call mpp_clock_end(iceClock2)
    call mpp_clock_end(iceClock)

  end subroutine update_ice_model_slow_dn

  subroutine update_ice_model_fast_new ( Atmos_boundary, Ice )
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
    type (ice_data_type),          intent(inout) :: Ice

    call mpp_clock_begin(iceClock)
    call mpp_clock_begin(iceClock3)
    call update_ice_model_fast_old (Ice, Atmos_boundary%fluxes,  &
                                         Atmos_boundary%u_flux,  &
                                         Atmos_boundary%v_flux,  &
                                         Atmos_boundary%u_star,  &
                                         Atmos_boundary%sw_flux_nir_dir, &
                                         Atmos_boundary%sw_flux_nir_dif, &
                                         Atmos_boundary%sw_flux_vis_dir, &
                                         Atmos_boundary%sw_flux_vis_dif, &
                                         Atmos_boundary%lw_flux, &
                                         Atmos_boundary%t_flux,  &
                                         Atmos_boundary%q_flux,  &
                                         Atmos_boundary%dhdt,    &
                                         Atmos_boundary%dedt,    &
                                         Atmos_boundary%drdt,    &
                                         Atmos_boundary%lprec,   &
                                         Atmos_boundary%fprec,   &
                                         Atmos_boundary%coszen,  &
                                         Atmos_boundary%p        )
    call mpp_clock_end(iceClock3)
    call mpp_clock_end(iceClock)

  end subroutine update_ice_model_fast_new

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! zero_top_quantities - zero fluxes to begin summing in ice fast physics       !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine zero_top_quantities ( Ice )
    type (ice_data_type), intent(inout)  :: Ice

    integer :: i, j, k, n, m

    Ice % avg_kount = 0

    do k = 1, km
       do j = jsd, jed
          do i = isd, ied
             Ice % flux_u_top(i,j,k)  = 0.0
             Ice % flux_v_top(i,j,k)  = 0.0
          enddo
       enddo
    enddo

    do k = 1, km
       do j = jsc, jec
          do i = isc, iec
             Ice % flux_t_top(i,j,k)          = 0.0
             Ice % flux_q_top(i,j,k)          = 0.0
             Ice % flux_lw_top(i,j,k)         = 0.0
             Ice % flux_lh_top(i,j,k)         = 0.0
             Ice % flux_sw_nir_dir_top(i,j,k) = 0.0
             Ice % flux_sw_nir_dif_top(i,j,k) = 0.0
             Ice % flux_sw_vis_dir_top(i,j,k) = 0.0
             Ice % flux_sw_vis_dif_top(i,j,k) = 0.0
             Ice % lprec_top(i,j,k)           = 0.0
             Ice % fprec_top(i,j,k)           = 0.0
             do n = 1, Ice%ocean_fluxes_top%num_bcs  !{
               do m = 1, Ice%ocean_fluxes_top%bc(n)%num_fields  !{
                 Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) = 0.0
               enddo  !} m
             enddo  !} n
          enddo
       enddo
    enddo

    do j = jsc, jec
       do i = isc, iec
          Ice % lwdn(i,j)        = 0.0
          Ice % swdn(i,j)        = 0.0
       enddo
    enddo

    return

  end subroutine zero_top_quantities

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! sum_top_quantities - sum fluxes for later use by ice/ocean slow physics      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine sum_top_quantities ( Ice, Atmos_boundary_fluxes, flux_u,  flux_v, flux_t,  flux_q, &
         flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif,&
         flux_lw, lprec,   fprec, flux_lh )
    type (ice_data_type),            intent(inout)  :: Ice
    type(coupler_3d_bc_type),        intent(inout)  :: Atmos_boundary_fluxes
    real, intent(in), dimension(isc:iec,jsc:jec,km) :: flux_u,  flux_v, flux_t, flux_q
    real, intent(in), dimension(isc:iec,jsc:jec,km) :: flux_lw, lprec, fprec, flux_lh
    real, intent(in), dimension(isc:iec,jsc:jec,km) :: flux_sw_nir_dir
    real, intent(in), dimension(isc:iec,jsc:jec,km) :: flux_sw_nir_dif
    real, intent(in), dimension(isc:iec,jsc:jec,km) :: flux_sw_vis_dir
    real, intent(in), dimension(isc:iec,jsc:jec,km) :: flux_sw_vis_dif
    real,dimension(isc:iec,jsc:jec)                 :: tmp
    integer                                         :: i, j, k, m, n

    if (Ice % avg_kount == 0) call zero_top_quantities (Ice)

    do k = 1, km
       do j = jsc, jec
          do i = isc, iec
             Ice % flux_u_top(i,j,k)  = Ice % flux_u_top(i,j,k)  + flux_u(i,j,k)
             Ice % flux_v_top(i,j,k)  = Ice % flux_v_top(i,j,k)  + flux_v(i,j,k)
             Ice % flux_t_top(i,j,k)  = Ice % flux_t_top(i,j,k)  + flux_t(i,j,k)
             Ice % flux_q_top(i,j,k)  = Ice % flux_q_top(i,j,k)  + flux_q(i,j,k)
             Ice % flux_sw_nir_dir_top(i,j,k) = Ice % flux_sw_nir_dir_top(i,j,k) + flux_sw_nir_dir(i,j,k)
             Ice % flux_sw_nir_dif_top(i,j,k) = Ice % flux_sw_nir_dif_top(i,j,k) + flux_sw_nir_dif(i,j,k)
             Ice % flux_sw_vis_dir_top(i,j,k) = Ice % flux_sw_vis_dir_top(i,j,k) + flux_sw_vis_dir(i,j,k)
             Ice % flux_sw_vis_dif_top(i,j,k) = Ice % flux_sw_vis_dif_top(i,j,k) + flux_sw_vis_dif(i,j,k)
             Ice % flux_lw_top(i,j,k) = Ice % flux_lw_top(i,j,k) + flux_lw(i,j,k)
             Ice % lprec_top(i,j,k)   = Ice % lprec_top(i,j,k)   + lprec(i,j,k)
             Ice % fprec_top(i,j,k)   = Ice % fprec_top(i,j,k)   + fprec(i,j,k)
             Ice % flux_lh_top(i,j,k) = Ice % flux_lh_top(i,j,k) + flux_lh(i,j,k)
          enddo
       enddo
    enddo

    do n = 1, Ice%ocean_fluxes_top%num_bcs  !{
       do m = 1, Ice%ocean_fluxes_top%bc(n)%num_fields  !{
          do k = 1, km
             do j = jsc, jec
                do i = isc, iec
                   Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) = Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) +     &
                        Atmos_boundary_fluxes%bc(n)%field(m)%values(i,j,k)
                enddo  !} m
             enddo  !} n
          enddo
       enddo
    enddo

    if (id_lwdn > 0) then
       tmp = all_avg(flux_lw+STEFAN*Ice%t_surf**4, Ice%part_size(isc:iec,jsc:jec,:))
       do j = jsc, jec
          do i = isc, iec
             if(Ice%mask(i,j) ) Ice%lwdn(i,j) = Ice%lwdn(i,j) + tmp(i,j)
          enddo
       enddo
    endif

    if (id_swdn > 0) then
       tmp = all_avg(flux_sw_vis_dir/(1-Ice%albedo_vis_dir)+flux_sw_vis_dif/(1-Ice%albedo_vis_dif)+ &
             flux_sw_nir_dir/(1-Ice%albedo_nir_dir)+flux_sw_nir_dif/(1-Ice%albedo_nir_dif), &
             Ice%part_size(isc:iec,jsc:jec,:) )
       do j = jsc, jec
          do i = isc, iec
             if(Ice%mask(i,j) ) Ice%swdn(i,j) = Ice%swdn(i,j) + tmp(i,j)
          enddo
       enddo
    endif

    Ice % avg_kount = Ice % avg_kount + 1

  end subroutine sum_top_quantities

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! avg_top_quantities - time average fluxes for ice and ocean slow physics      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine avg_top_quantities ( Ice )
    type (ice_data_type), intent(inout)  :: Ice

    real    :: u, v, divid
    integer :: i, j, k, m, n
    logical :: sent
    !
    ! compute average fluxes
    !
    if (Ice % avg_kount == 0) call error_mesg ('avg_top_quantities', &
         'no ocean model fluxes have been averaged', 3)

    divid = 1./float(Ice%avg_kount)

    do k=1,km
       do j = jsc, jec
          do i = isc, iec
             u                       = Ice % flux_u_top(i,j,k)  * divid
             v                       = Ice % flux_v_top(i,j,k)  * divid
             Ice % flux_u_top(i,j,k) = u*cos_rot(i,j)-v*sin_rot(i,j) ! rotate stress from lat/lon
             Ice % flux_v_top(i,j,k) = v*cos_rot(i,j)+u*sin_rot(i,j) ! to ocean coordinates
          end do
       end do
    end do

    if (atmos_winds) then ! put wind stress on u,v points and change sign to +down
       call mpp_update_domains(Ice % flux_u_top, Ice % flux_v_top, Domain  )
       do k=1,km
          call vel_t_to_uv( -Ice%flux_u_top(:,:,k),-Ice%flux_v_top(:,:,k), &
               Ice%flux_u_top_bgrid(isc:iec,jsc:jec,k), Ice%flux_v_top_bgrid(isc:iec,jsc:jec,k) )
       end do
    else
       if(reproduce_siena_201303) then
          Ice%flux_u_top_bgrid(isc:iec,jsc:jec,:) = Ice%flux_u_top(isc:iec,jsc:jec,:)
          Ice%flux_v_top_bgrid(isc:iec,jsc:jec,:) = Ice%flux_v_top(isc:iec,jsc:jec,:)
       else
          call mpp_update_domains(Ice % flux_u_top, Ice % flux_v_top, Domain  )
          do k=1,km
             call vel_t_to_uv( Ice%flux_u_top(:,:,k),Ice%flux_v_top(:,:,k), &
                  Ice%flux_u_top_bgrid(isc:iec,jsc:jec,k), Ice%flux_v_top_bgrid(isc:iec,jsc:jec,k) )
          end do
       endif
    endif

    do k = 1, km
       do j = jsc, jec
          do i = isc, iec
             Ice % flux_t_top(i,j,k)  = Ice % flux_t_top(i,j,k)  * divid
             Ice % flux_q_top(i,j,k)  = Ice % flux_q_top(i,j,k)  * divid
             Ice % flux_sw_nir_dir_top(i,j,k) = Ice % flux_sw_nir_dir_top(i,j,k) * divid
             Ice % flux_sw_nir_dif_top(i,j,k) = Ice % flux_sw_nir_dif_top(i,j,k) * divid
             Ice % flux_sw_vis_dir_top(i,j,k) = Ice % flux_sw_vis_dir_top(i,j,k) * divid
             Ice % flux_sw_vis_dif_top(i,j,k) = Ice % flux_sw_vis_dif_top(i,j,k) * divid
             Ice % flux_lw_top(i,j,k) = Ice % flux_lw_top(i,j,k) * divid
             Ice % fprec_top(i,j,k)   = Ice % fprec_top(i,j,k)   * divid
             Ice % lprec_top(i,j,k)   = Ice % lprec_top(i,j,k)   * divid
             Ice % flux_lh_top(i,j,k) = Ice % flux_lh_top(i,j,k) * divid
             do n = 1, Ice%ocean_fluxes_top%num_bcs  !{
               do m = 1, Ice%ocean_fluxes_top%bc(n)%num_fields  !{
                 Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) = Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) * divid
               enddo  !} m
             enddo  !} n
          enddo
       enddo
    enddo

    do j = jsc, jec
       do i = isc, iec
          Ice % lwdn(i,j) = Ice % lwdn(i,j)* divid
          Ice % swdn(i,j) = Ice % swdn(i,j)* divid
       enddo
    enddo
    !
    ! Flux diagnostics
    !
    if (id_sh>0) sent = send_data(id_sh, all_avg(Ice%flux_t_top,Ice%part_size(isc:iec,jsc:jec,:)),     &
         Ice%Time, mask=Ice%mask)
    if (id_lh>0) sent = send_data(id_lh, all_avg(Ice%flux_lh_top,Ice%part_size(isc:iec,jsc:jec,:)),    &
         Ice%Time, mask=Ice%mask)
    if (id_evap>0) sent = send_data(id_evap, all_avg(Ice%flux_q_top,Ice%part_size(isc:iec,jsc:jec,:)), &
         Ice%Time, mask=Ice%mask)
    if (id_sw>0) sent = send_data(id_sw, all_avg(Ice%flux_sw_vis_dir_top+Ice%flux_sw_vis_dif_top+ &
                        Ice%flux_sw_nir_dir_top+Ice%flux_sw_nir_dif_top,Ice%part_size(isc:iec,jsc:jec,:)), &
                        Ice%Time, mask=Ice%mask)
    if (id_lw>0) sent = send_data(id_lw, all_avg(Ice%flux_lw_top,Ice%part_size(isc:iec,jsc:jec,:)),    &
         Ice%Time, mask=Ice%mask)
    if (id_snofl>0) sent = send_data(id_snofl, all_avg(Ice%fprec_top,Ice%part_size(isc:iec,jsc:jec,:)),&
         Ice%Time, mask=Ice%mask)
    if (id_rain>0) sent = send_data(id_rain, all_avg(Ice%lprec_top,Ice%part_size(isc:iec,jsc:jec,:)),  &
         Ice%Time, mask=Ice%mask)
    if (id_lwdn>0) sent = send_data(id_lwdn, Ice%lwdn, Ice%Time, mask=Ice%mask)
    if (id_swdn>0) sent = send_data(id_swdn, Ice%swdn, Ice%Time, mask=Ice%mask)
    if (id_sw_vis>0) sent = send_data(id_sw_vis, all_avg(Ice%flux_sw_vis_dif_top+Ice%flux_sw_vis_dir_top, &
                            Ice%part_size(isc:iec,jsc:jec,:)), Ice%Time, mask=Ice%mask)
    if (id_sw_nir_dir>0) sent = send_data(id_sw_nir_dir, all_avg(Ice%flux_sw_nir_dir_top,Ice%part_size(isc:iec,jsc:jec,:)),    &
         Ice%Time, mask=Ice%mask)
    if (id_sw_nir_dif>0) sent = send_data(id_sw_nir_dif, all_avg(Ice%flux_sw_nir_dif_top,Ice%part_size(isc:iec,jsc:jec,:)),    &
         Ice%Time, mask=Ice%mask)
    if (id_sw_vis_dir>0) sent = send_data(id_sw_vis_dir, all_avg(Ice%flux_sw_vis_dir_top,Ice%part_size(isc:iec,jsc:jec,:)),    &
         Ice%Time, mask=Ice%mask)
    if (id_sw_vis_dif>0) sent = send_data(id_sw_vis_dif, all_avg(Ice%flux_sw_vis_dif_top,Ice%part_size(isc:iec,jsc:jec,:)),    &
         Ice%Time, mask=Ice%mask)

    !
    ! set count to zero and fluxes will be zeroed before the next sum
    !
    Ice % avg_kount = 0

  end subroutine avg_top_quantities

  subroutine ice_top_to_ice_bottom (Ice, part_size, part_size_uv)
    type (ice_data_type), intent(inout) :: Ice
    real, dimension (:,:,:), intent(in) :: part_size, part_size_uv
    integer                             :: m, n

    Ice % flux_u  = all_avg( Ice % flux_u_top_bgrid(isc:iec,jsc:jec,:) , part_size_uv )
    Ice % flux_v  = all_avg( Ice % flux_v_top_bgrid(isc:iec,jsc:jec,:) , part_size_uv )
    Ice % flux_t  = all_avg( Ice % flux_t_top , part_size )
    Ice % flux_q  = all_avg( Ice % flux_q_top , part_size )
    Ice % flux_sw_nir_dir = all_avg( Ice % flux_sw_nir_dir_top, part_size )
    Ice % flux_sw_nir_dif = all_avg( Ice % flux_sw_nir_dif_top, part_size )
    Ice % flux_sw_vis_dir = all_avg( Ice % flux_sw_vis_dir_top, part_size )
    Ice % flux_sw_vis_dif = all_avg( Ice % flux_sw_vis_dif_top, part_size )
    Ice % flux_lw = all_avg( Ice % flux_lw_top, part_size )
    Ice % fprec   = all_avg( Ice % fprec_top  , part_size )
    Ice % lprec   = all_avg( Ice % lprec_top  , part_size )
    Ice % flux_lh = all_avg( Ice % flux_lh_top, part_size )
    do n = 1, Ice%ocean_fluxes%num_bcs  !{
      do m = 1, Ice%ocean_fluxes%bc(n)%num_fields  !{
        Ice%ocean_fluxes%bc(n)%field(m)%values =                &
             all_avg(Ice%ocean_fluxes_top%bc(n)%field(m)%values, part_size)
      enddo  !} m
    enddo  !} n

  end subroutine ice_top_to_ice_bottom

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_bottom_to_ice_top - prepare surface state for atmosphere fast physics    !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_bottom_to_ice_top (Ice, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                    frazil_ice_bot, Ocean_ice_boundary, s_surf_ice_bot, sea_lev_ice_bot )
    type (ice_data_type),                    intent(inout) :: Ice
    real, dimension(isc:iec,jsc:jec),           intent(in) :: t_surf_ice_bot, u_surf_ice_bot
    real, dimension(isc:iec,jsc:jec),           intent(in) :: v_surf_ice_bot, frazil_ice_bot
    type(ocean_ice_boundary_type),           intent(inout) :: Ocean_ice_boundary
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: s_surf_ice_bot, sea_lev_ice_bot

    real, dimension(isc:iec,jsc:jec) :: sst, tmp
    real                             :: u, v
    integer                          :: i, j, k, m, n
    logical                          :: sent
    real, parameter                  :: LI = hlf

    !
    ! pass ocean state through ice on first partition
    !
    if (.not. spec_ice) & ! otherwise, already set by update_ice_model_slow
         Ice%t_surf(:,:,1) = t_surf_ice_bot

    if (do_init) then
       call get_sea_surface(Ice%Time, Ice%t_surf(:,:,1), Ice%part_size(isc:iec,jsc:jec,1:2), &
            Ice%h_ice    (isc:iec,jsc:jec,2  ) )
       call mpp_update_domains(Ice%part_size(:,:,1:2), Domain ) ! these two updates cannot be combined
       call mpp_update_domains(Ice%h_ice(:,:,2), Domain )       ! these two updates cannot be combined
       Ice%part_size_uv(:,:,1) = 1.0
       do k=2,km
          call t_to_uv(Ice%part_size(:,:,k),Ice%part_size_uv(:,:,k))
       enddo
       do k=2,km
          do j=jsc,jec
             do i=isc,iec
                Ice%part_size_uv (i,j,1) = Ice%part_size_uv(i,j,1)-Ice%part_size_uv (i,j,k)
             enddo
          end do
       end do
       do_init = .false.
    end if

    if (first_time .and. conservation_check) then
       h2o    = 0.0
       h2o(1) = sum(cell_area*all_avg(DI*Ice%h_ice(isc:iec,jsc:jec,:)+   &
                DS*Ice%h_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)))
       heat   = 0.0

       do k=2,km
          do j=jsc,jec
             do i=isc,iec
                if ((Ice%part_size(i,j,k)>0.0.and.Ice%h_ice(i,j,k)>0.0)) then
                   if (slab_ice) then
                      heat(1) = heat(1) - cell_area(i,j) * Ice%part_size(i,j,k) * Ice%h_ice(i,j,2)*DI*LI
                   else
                      heat(1) = heat(1) - cell_area(i,j) * Ice%part_size(i,j,k) *        &
                                e_to_melt(Ice%h_snow(i,j,k), Ice%h_ice(i,j,k)/2,         &
                                Ice%t_ice1(i,j,k), Ice%h_ice(i,j,k)/2, Ice%t_ice2(i,j,k))
                   end if
                end if
             end do
          end do
       end do
       first_time = .false.;
    end if
    sst = Ice%t_surf(:,:,1)

    if (present(s_surf_ice_bot)) then
       do j = jsc, jec
          do i = isc, iec     
             Ice%s_surf(i,j) = s_surf_ice_bot(i,j)
          enddo
       enddo
    else
       do j = jsc, jec
          do i = isc, iec  
             Ice%s_surf(i,j) = 0.0
          enddo
       enddo
    end if

    do j = jsc,jec
       do i = isc,iec
          Ice%frazil(i,j) = frazil_ice_bot(i,j)
       enddo
    enddo

!       transfer the ocean state for extra tracer fluxes
!

    do n = 1, Ocean_ice_boundary%fields%num_bcs  !{
      do m = 1, Ocean_ice_boundary%fields%bc(n)%num_fields  !{
        Ice%ocean_fields%bc(n)%field(m)%values(:,:,1) = Ocean_ice_boundary%fields%bc(n)%field(m)%values
      enddo  !} m
    enddo  !} n

    if (present(sea_lev_ice_bot)) then
       do j = jsc, jec
          do i = isc, iec  
             Ice%sea_lev(i,j) = sea_lev_ice_bot(i,j)
          enddo
       enddo
    else
       do j = jsc, jec
          do i = isc, iec
             Ice%sea_lev(i,j) = 0.0
          enddo
       enddo
    end if

    call mpp_update_domains(Ice%sea_lev, Domain)

    do k = 2, km
       do j = jsc,jec
          do i = isc,iec  
             Ice % tmelt(i,j,k)   = 0.0
             Ice % bmelt(i,j,k)   = 0.0
          enddo
       enddo
    enddo

    tmp = ice_avg(Ice%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:) )

    do j = jsc, jec
       do i = isc, iec
          if( tmp(i,j) > 0.0) then
             Ice%bheat(i,j) = kmelt*(sst(i,j)-Tfreeze+MU_TS*Ice%s_surf(i,j))
          else
             Ice%bheat(i,j) = 0.0
          endif
       enddo
    enddo

    do j = jsc, jec
       do i = isc, iec
          Ice%ice_mask(i,j,1) = .false.
       enddo
    enddo

    do k = 2, km
       do j = jsc, jec
          do i = isc, iec
             if(Ice%h_ice(i,j,k) > 0.0) then
                Ice%ice_mask(i,j,k) = .true.
             else
                Ice%ice_mask(i,j,k) = .false.
             endif
          enddo
       enddo
    enddo

    do k=2,km
       do j=jsc, jec
          do i=isc, iec
             if (Ice%ice_mask(i,j,k)) then
                call ice_optics(Ice%albedo(i,j,k), Ice%pen(i,j,k), Ice%trn(i,j,k), &
                     Ice%h_snow(i,j,k), Ice%h_ice(i,j,k), &
                     Ice%t_surf(i,j,k)-Tfreeze, -MU_TS*Ice%s_surf(i,j) )
             end if
          end do
       end do
    end do

!!! THIS NEEDS TO BE DEFINED PROPERLY -- ARE THEY NEEDED ON RESTART ??
!!!  FOR NOW, simply set all values to that of Ice%albedo, In general,
!! the ability to read from restart file and initialize when none
!! present should be incorporated above.
    if (cm2_bugs) then
       do k = 1, km 
          do j = jsc, jec
             do i = isc, iec
                Ice%albedo_vis_dir(i,j,k) = Ice%albedo(i,j,k)
                Ice%albedo_nir_dir(i,j,k) = Ice%albedo(i,j,k)
                Ice%albedo_vis_dif(i,j,k) = Ice%albedo(i,j,k)
                Ice%albedo_nir_dif(i,j,k) = Ice%albedo(i,j,k)
             enddo
          enddo
       enddo
    else
       do k = 2, km 
          do j = jsc, jec
             do i = isc, iec
                Ice%albedo_vis_dir(i,j,k) = Ice%albedo(i,j,k)
                Ice%albedo_nir_dir(i,j,k) = Ice%albedo(i,j,k)
                Ice%albedo_vis_dif(i,j,k) = Ice%albedo(i,j,k)
                Ice%albedo_nir_dif(i,j,k) = Ice%albedo(i,j,k)
             enddo
          enddo
       enddo
    end if

    do j = jsc, jec
       do i = isc, iec
          Ice%u_ocn(i,j) = u_surf_ice_bot(i,j) ! need under-ice current
          Ice%v_ocn(i,j) = v_surf_ice_bot(i,j) ! for water drag term
       enddo
    enddo

    if(reproduce_siena_201303) then
       call mpp_update_domains(Ice%u_ocn, Ice%v_ocn, Domain)
    else
       call mpp_update_domains(Ice%u_ocn, Ice%v_ocn, Domain, gridtype=BGRID_NE)
    endif

    ! put ocean and ice velocities into Ice%u_surf/v_surf on t-cells
    call uv_to_t(Ice%u_ocn, Ice%u_surf(:,:,1))
    call uv_to_t(Ice%v_ocn, Ice%v_surf(:,:,1))

    call uv_to_t(Ice%u_ice, Ice%u_surf(:,:,2))
    call uv_to_t(Ice%v_ice, Ice%v_surf(:,:,2))

    do k=1,2
       do j = jsc, jec
          do i = isc, iec
             u = Ice % u_surf(i,j,k)
             v = Ice % v_surf(i,j,k)
             Ice % u_surf(i,j,k) =  u*cos_rot(i,j)+v*sin_rot(i,j) ! rotate velocity from ocean
             if (cm2_bugs) then
                Ice % v_surf(i,j,k) =  -v*cos_rot(i,j)+u*sin_rot(i,j) ! coord. to lat/lon coord.
             else
                Ice % v_surf(i,j,k) =  v*cos_rot(i,j)-u*sin_rot(i,j) ! coord. to lat/lon coord.
             endif
          end do
       end do
    end do

    do k=3,km
       do j = jsc, jec
          do i = isc, iec
             Ice%u_surf(i,j,k) = Ice%u_surf(i,j,2)  ! same ice flow on all ice partitions
             Ice%v_surf(i,j,k) = Ice%v_surf(i,j,2)  !
          end do
       end do
    end do
    !
    ! Pre-timestep diagnostics
    !
    if (id_sst>0) sent = send_data(id_sst, sst-Tfreeze, Ice%Time, mask=Ice%mask)
    if (id_sss>0) sent = send_data(id_sss, Ice%s_surf , Ice%Time, mask=Ice%mask)
    if (id_ssh>0) sent = send_data(id_ssh, Ice%sea_lev(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
    if (id_uo >0) sent = send_data(id_uo , Ice%u_ocn(isc:iec,jsc:jec)  , Ice%Time, mask=Ice%mask)
    if (id_vo >0) sent = send_data(id_vo , Ice%v_ocn(isc:iec,jsc:jec)  , Ice%Time, mask=Ice%mask)
    if (id_bheat>0) sent = send_data(id_bheat, Ice%bheat, Ice%Time, mask=Ice%mask)

  end subroutine ice_bottom_to_ice_top

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! update_ice_model_fast - records fluxes (in Ice) and calculates ice temp. on  !
  !                         (fast) atmospheric timestep (see coupler_main.f90)   !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine update_ice_model_fast_old (Ice, Atmos_boundary_fluxes, flux_u,  flux_v, u_star, &
       flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif,&
       flux_lw, flux_t, flux_q, dhdt, dedt, drdt, lprec, fprec, coszen, p_surf )

    type (ice_data_type),             intent(inout) :: Ice
    type(coupler_3d_bc_type),         intent(inout) :: Atmos_boundary_fluxes
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_u,  flux_v  ! surface stress (+ up)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: u_star           ! friction velocity (for ocean roughness)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_lw          ! net longwave radiation (+ down)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_t           ! sensible heat flux (+ up)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_q           ! specific humidity flux (+up)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_nir_dir ! net near IR direct shortwave radiation (+ down)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_nir_dif ! net near IR diffuse shortwave radiation (+ down)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_vis_dir ! net visible direct shortwave radiation (+ down)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_vis_dif ! net visible diffuse shortwave radiation (+ down)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: dhdt, dedt, drdt ! d(flux)/d(surf_temp) (+ up)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: lprec, fprec     ! liquid & frozen precip. rate [kg/(m^2*sec)](+ down)
    real, dimension(isc:iec,jsc:jec,km), intent(in) :: coszen           ! cosine of the zenith angle
    real, dimension(isc:iec,jsc:jec,km), intent(in), optional :: p_surf

    real, dimension(isc:iec,jsc:jec,km) :: flux_t_new, flux_q_new, flux_lh_new, flux_lw_new
    real, dimension(isc:iec,jsc:jec,km) :: flux_sw_nir_dir_new
    real, dimension(isc:iec,jsc:jec,km) :: flux_sw_nir_dif_new
    real, dimension(isc:iec,jsc:jec,km) :: flux_sw_vis_dir_new
    real, dimension(isc:iec,jsc:jec,km) :: flux_sw_vis_dif_new
    real, dimension(isc:iec,jsc:jec,km) :: flux_u_new, flux_v_new, lprec_new, fprec_new
    integer                             :: dy, sc, i, j, k
    real                                :: dt_fast, ts_new, dts, hf, hfd, latent
    logical                             :: sent
    real                                :: gmt, time_since_ae, cosz, rrsun, fracday, fracday_dt_ice, fracday_day
    real, dimension(isc:iec,jsc:jec)    :: diurnal_factor
    real                                :: rad, cosz_day, cosz_dt_ice, rrsun_day, rrsun_dt_ice
    type (time_type)                    :: Dt_ice
    real, dimension(isc:iec,jsc:jec)    :: cosz_alb

    if (id_alb>0) sent = send_data(id_alb, all_avg(Ice%albedo,Ice%part_size(isc:iec,jsc:jec,:)), &
         Ice%Time, mask=Ice%mask)
    !
    ! set up local copies of fluxes for modification
    !
    do k = 1, km
       do j = jsc, jec
          do i = isc, iec
             flux_u_new(i,j,k)  = flux_u(i,j,k)
             flux_v_new(i,j,k)  = flux_v(i,j,k)
             flux_t_new(i,j,k)  = flux_t(i,j,k)
             flux_q_new(i,j,k)  = flux_q(i,j,k)
             flux_lh_new(i,j,k) = hlv*flux_q(i,j,k)
             flux_lw_new(i,j,k) = flux_lw(i,j,k)
             flux_sw_nir_dir_new(i,j,k) = flux_sw_nir_dir(i,j,k)
             flux_sw_nir_dif_new(i,j,k) = flux_sw_nir_dif(i,j,k)
             flux_sw_vis_dir_new(i,j,k) = flux_sw_vis_dir(i,j,k)
             flux_sw_vis_dif_new(i,j,k) = flux_sw_vis_dif(i,j,k)
             lprec_new(i,j,k)   = lprec(i,j,k)
             fprec_new(i,j,k)   = fprec(i,j,k)
          enddo
       enddo
    enddo
    if (add_diurnal_sw .or. do_sun_angle_for_alb) then
!---------------------------------------------------------------------
!    extract time of day (gmt) from time_type variable time with
!    function universal_time.
!---------------------------------------------------------------------
      gmt = universal_time(Ice%Time)
!---------------------------------------------------------------------
!    extract the time of year relative to the northern hemisphere
!    autumnal equinox (time_since_ae) from time_type variable 
!    time using the function orbital_time.
!---------------------------------------------------------------------
      time_since_ae = orbital_time(Ice%Time)
!--------------------------------------------------------------------
!    call diurnal_solar_2d to calculate astronomy fields
!    convert geo_lon and geo_lat to radians
!    Per Rick Hemler:
!      call daily_mean_solar to get cosz (over a day)
!      call diurnal_solar with dtime=Dt_ice to get cosz over Dt_ice
!      diurnal_factor = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice/cosz_day*fracday_day*rrsun_day
!--------------------------------------------------------------------
      rad = acos(-1.)/180.
      Dt_ice = Ice%Time_step_fast
    endif
    if (add_diurnal_sw) then
      do j = jsc, jec
        do i = isc, iec
          call diurnal_solar(geo_lat(i,j)*rad, geo_lon(i,j)*rad, Ice%time, cosz=cosz_dt_ice,     &
                             fracday=fracday_dt_ice, rrsun=rrsun_dt_ice, dt_time=Dt_ice)
          call daily_mean_solar (geo_lat(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
          diurnal_factor(i,j) = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice /   &
                                       max(1e-30, cosz_day*fracday_day*rrsun_day)
        enddo
      enddo

      do k = 1, km
        do j = jsc, jec
          do i = isc, iec
            flux_sw_nir_dir_new(i,j,k) = flux_sw_nir_dir(i,j,k) * diurnal_factor(i,j)
            flux_sw_nir_dif_new(i,j,k) = flux_sw_nir_dif(i,j,k) * diurnal_factor(i,j)
            flux_sw_vis_dir_new(i,j,k) = flux_sw_vis_dir(i,j,k) * diurnal_factor(i,j)
            flux_sw_vis_dif_new(i,j,k) = flux_sw_vis_dif(i,j,k) * diurnal_factor(i,j)
          enddo
        enddo
      enddo
    endif

    Ice%p_surf = 0.0
    if (present(p_surf)) &
      Ice%p_surf = all_avg(p_surf,Ice%part_size(isc:iec,jsc:jec,:))

    !
    ! implicit update of ice surface temperature
    !
    call get_time(Ice%Time_step_fast,sc,dy); dt_fast = 864e2*dy+sc

    do k = 2, km
       do j = jsc, jec
          do i=isc, iec
             if (Ice%ice_mask(i,j,k)) then
                if (Ice%h_snow(i,j,k)>0.0 .or. slab_ice) then
                   flux_lh_new(i,j,k) = flux_lh_new(i,j,k) + hlf*flux_q(i,j,k)
                   latent             = hlv+hlf
                else
                   flux_lh_new(i,j,k) = flux_lh_new(i,j,k)+hlf*flux_q(i,j,k)*(1-TFI/Ice%t_ice1(i,j,k))
                   latent             = hlv+hlf*(1-TFI/Ice%t_ice1(i,j,k))
                end if
                hfd = dhdt(i,j,k) + dedt(i,j,k)*latent + drdt(i,j,k)
                hf  = flux_t(i,j,k) + flux_q(i,j,k)*latent - flux_lw(i,j,k)                  &
                      - (1-Ice%pen(i,j,k))* (flux_sw_vis_dir(i,j,k)+flux_sw_vis_dif(i,j,k)+  &
                        flux_sw_nir_dir(i,j,k)+flux_sw_nir_dif(i,j,k))                       &
                      - hfd*(Ice%t_surf(i,j,k)-Tfreeze)
                call ice3lay_temp(Ice%h_snow(i,j,k), Ice%h_ice(i,j,k), Ice%t_ice1(i,j,k),    &
                                  Ice%t_ice2(i,j,k), ts_new, hf, hfd, Ice%pen(i,j,k)*        &
                                  (1-Ice%trn(i,j,k))*(flux_sw_vis_dir(i,j,k)+                &
                                  flux_sw_vis_dif(i,j,k)+flux_sw_nir_dir(i,j,k)+             &
                                  flux_sw_nir_dif(i,j,k)), -MU_TS*Ice%s_surf(i,j),           &
                                  Ice%bheat(i,j), dt_fast, Ice%tmelt(i,j,k), Ice%bmelt(i,j,k))
                dts                = ts_new-(Ice%t_surf(i,j,k)-Tfreeze)
                Ice%t_surf(i,j,k)  = Ice%t_surf(i,j,k)  + dts
                flux_t_new(i,j,k)  = flux_t_new(i,j,k)  + dts * dhdt(i,j,k)
                flux_q_new(i,j,k)  = flux_q_new(i,j,k)  + dts * dedt(i,j,k)
                flux_lh_new(i,j,k) = flux_lh_new(i,j,k) + dts * dedt(i,j,k) * latent
                flux_lw_new(i,j,k) = flux_lw_new(i,j,k) - dts * drdt(i,j,k)
             end if
          end do
       end do
    end do

    call compute_ocean_roughness (Ice%mask, u_star(:,:,1), Ice%rough_mom(:,:,1), &
                                  Ice%rough_heat(:,:,1), Ice%rough_moist(:,:,1)  )

    if(cm2_bugs) then
       call compute_ocean_albedo (Ice%mask, coszen(:,:,1), Ice%albedo(:,:,1), latitude )
       Ice%albedo_vis_dir(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_nir_dir(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_vis_dif(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_nir_dif(:,:,1) = Ice%albedo(:,:,1)
    elseif (do_sun_angle_for_alb) then
      call diurnal_solar(geo_lat*rad, geo_lon*rad, Ice%time, cosz=cosz_alb,	&
			     fracday=diurnal_factor, rrsun=rrsun_dt_ice, dt_time=Dt_ice)  !diurnal_factor as dummy
      call compute_ocean_albedo (Ice%mask, cosz_alb(:,:), Ice%albedo_vis_dir(:,:,1),&
                                   Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                                   Ice%albedo_nir_dif(:,:,1), latitude )
    else
      call compute_ocean_albedo (Ice%mask, coszen(:,:,1), Ice%albedo_vis_dir(:,:,1),&
                                 Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                                 Ice%albedo_nir_dif(:,:,1), latitude )
    endif

    if (id_alb_vis_dir>0) sent = send_data(id_alb_vis_dir, all_avg(Ice%albedo_vis_dir,Ice%part_size(isc:iec,jsc:jec,:)), &
         Ice%Time, mask=Ice%mask)
    if (id_alb_vis_dif>0) sent = send_data(id_alb_vis_dif, all_avg(Ice%albedo_vis_dif,Ice%part_size(isc:iec,jsc:jec,:)), &
         Ice%Time, mask=Ice%mask)
    if (id_alb_nir_dir>0) sent = send_data(id_alb_nir_dir, all_avg(Ice%albedo_nir_dir,Ice%part_size(isc:iec,jsc:jec,:)), &
         Ice%Time, mask=Ice%mask)
    if (id_alb_nir_dif>0) sent = send_data(id_alb_nir_dif, all_avg(Ice%albedo_nir_dif,Ice%part_size(isc:iec,jsc:jec,:)), &
         Ice%Time, mask=Ice%mask)

    call sum_top_quantities ( Ice, Atmos_boundary_fluxes, flux_u_new,  flux_v_new, flux_t_new, &
      flux_q_new, flux_sw_nir_dir_new, flux_sw_nir_dif_new, flux_sw_vis_dir_new,               &
      flux_sw_vis_dif_new, flux_lw_new, lprec_new,   fprec_new,  flux_lh_new )

    Ice%Time = Ice%Time + Ice%Time_step_fast ! advance time

  end subroutine update_ice_model_fast_old

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! update_ice_model_slow - do ice dynamics, transport, and mass changes         !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine update_ice_model_slow (Ice, runoff, calving, &
                                         runoff_hflx, calving_hflx, p_surf)

    type (ice_data_type),             intent(inout)        :: Ice
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: runoff, calving
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: runoff_hflx, calving_hflx
    real, dimension(:,:,:),           intent(in), optional :: p_surf ! obsolete

    real, dimension(isc:iec,jsc:jec)      :: fx_ice, fy_ice
    real, dimension(isc:iec,jsc:jec)      :: fx_cor, fy_cor
    real, dimension(isc:iec,jsc:jec)      :: fx_wat, fy_wat
    real, dimension(isc:iec,jsc:jec)      :: hi_change, h2o_change, bsnk, x
    real, dimension(isc:iec,jsc:jec,2:km) :: snow_to_ice
    real, dimension(isc:iec,jsc:jec,km)   :: part_save, part_save_uv
    real, dimension(isc:iec,jsc:jec)      :: dum1, Obs_h_ice ! for qflux calculation
    real, dimension(isc:iec,jsc:jec,2)    :: Obs_cn_ice      ! partition 2 = ice concentration
    real, dimension(isd:ied,jsd:jed)      :: tmp1, tmp2
    real, dimension(2:km)                 :: e2m
    integer                               :: i, j, k, l, sc, dy, iyr, imon, iday, ihr, imin, isec
    real                                  :: dt_slow, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic, bablt
    real                                  :: heat_limit_ice, heat_res_ice
    real                                  :: tot_heat, heating, tot_frazil
    logical                               :: sent
    real, parameter                       :: LI = hlf

    call get_time(Ice%Time_step_slow,sc,dy); dt_slow = 864e2*dy+sc
    !
    ! Fluxes
    !
    if (present( runoff)) then ! save liquid runoff for ocean
       do j = jsc, jec
          do i = isc, iec
             Ice % runoff(i,j)  = runoff(i,j)
          enddo
       enddo
       if (id_runoff>0) &
            sent = send_data( id_runoff, Ice%runoff, Ice%Time, mask=Ice%mask )
    else
       do j = jsc, jec
          do i = isc, iec
             Ice % runoff(i,j) = 0.0
          enddo
       enddo
    endif

    if (present(calving)) then ! save frozen runoff for ocean
       do j = jsc, jec
          do i = isc, iec
             Ice % calving(i,j) = calving(i,j)
          enddo
       enddo
       if (id_calving>0) &
            sent = send_data(id_calving, Ice%calving, Ice%Time, mask=Ice%mask )
    else
       do j = jsc, jec
          do i = isc, iec
             Ice % calving(i,j) = 0.0
          enddo
       enddo
    endif

    if (present( runoff_hflx)) then ! save liquid runoff hflx for ocean
       do j = jsc, jec
          do i = isc, iec
             Ice % runoff_hflx(i,j)  = runoff_hflx(i,j)
          enddo
       enddo
       if (id_runoff_hflx>0) &
            sent = send_data( id_runoff_hflx, Ice%runoff_hflx, Ice%Time, mask=Ice%mask )
    else
       do j = jsc, jec
          do i = isc, iec
             Ice % runoff_hflx(i,j) = 0.0
          enddo
       enddo
    endif

    if (present(calving_hflx)) then ! save frozen runoff hflx for ocean
       do j = jsc, jec
          do i = isc, iec
             Ice % calving_hflx(i,j) = calving_hflx(i,j)
          enddo
       enddo
       if (id_calving_hflx>0) &
            sent = send_data(id_calving_hflx, Ice%calving_hflx, Ice%Time, mask=Ice%mask )
    else
       do j = jsc, jec
          do i = isc, iec
             Ice % calving_hflx(i,j) = 0.0
          enddo
       enddo
    endif

    !TOM> assume that open water area is not up to date:
    call mpp_clock_end(iceClock)
    call mpp_clock_end(iceClock2)
    tmp1=1.-max(1.-sum(Ice%part_size(:,:,2:km),dim=3),0.0)
    tmp2=ice_avg(Ice%h_ice,Ice%part_size)
    ! Calve off icebergs and integrate forward iceberg trajectories
    if (do_icebergs) call icebergs_run( Ice%icebergs,Ice%Time,                 &
                      Ice%calving, Ice%u_ocn, Ice%v_ocn, Ice%u_ice, Ice%v_ice, &
                      Ice%flux_u, Ice%flux_v, Ice%sea_lev, Ice%t_surf(:,:,1),  &
                      Ice%calving_hflx, tmp1, tmp2)
    call mpp_clock_begin(iceClock2)
    call mpp_clock_begin(iceClock)

    call avg_top_quantities(Ice) ! average fluxes from update_ice_model_fast

    do k = 1, km
       do j = jsc, jec
          do i = isc, iec
             part_save(i,j,k)    = Ice%part_size(i,j,k)
             part_save_uv(i,j,k) = Ice%part_size_uv(i,j,k)
          enddo
       enddo
    enddo
    !
    ! conservation checks: top fluxes
    !
    call mpp_clock_begin(iceClock7)
    if (conservation_check) then
       h2o(2) = h2o(2)+dt_slow*sum(cell_area*(Ice%runoff+Ice%calving                                  &
               +all_avg(Ice%lprec_top+Ice%fprec_top-Ice%flux_q_top,Ice%part_size(isc:iec,jsc:jec,:))),&
                cell_area > 0 )
       heat(2) = heat(2)+dt_slow*sum(cell_area*(                              &
                 all_avg(Ice%flux_sw_vis_dir_top+Ice%flux_sw_vis_dif_top      &
                +Ice%flux_sw_nir_dir_top+Ice%flux_sw_nir_dif_top,             &
                 Ice%part_size(isc:iec,jsc:jec,:))                            &
                +all_avg(Ice%flux_lw_top, Ice%part_size(isc:iec,jsc:jec,:))   &
                -all_avg(Ice%flux_t_top , Ice%part_size(isc:iec,jsc:jec,:))   &
                -all_avg(Ice%flux_lh_top, Ice%part_size(isc:iec,jsc:jec,:))   &
                -LI*(all_avg(Ice%fprec_top,Ice%part_size(isc:iec,jsc:jec,:))  &
                +Ice%calving)), cell_area > 0 )
       tot_frazil = sum(cell_area*Ice%frazil)
    end if
    call mpp_clock_end(iceClock7)

    ! Dynamics
    !

    call mpp_clock_begin(iceClock4)
    tmp1 = ice_avg(Ice%h_snow,Ice%part_size)
    tmp2 = ice_avg(Ice%h_ice,Ice%part_size)

    call mpp_clock_begin(iceClocka)
    if(reproduce_siena_201303) then
       call ice_dynamics(1-Ice%part_size(:,:,1), tmp1, tmp2, Ice%u_ice, Ice%v_ice,                      &
                         Ice%sig11, Ice%sig22, Ice%sig12, Ice%u_ocn, Ice%v_ocn,                         &
                         ice_avg(Ice%flux_u_top_bgrid(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:) ),  &
                         ice_avg(Ice%flux_v_top_bgrid(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:) ),  &
                         Ice%sea_lev, fx_wat, fy_wat, fx_ice, fy_ice, fx_cor, fy_cor)
    else
       call ice_dynamics(1-Ice%part_size(:,:,1), tmp1, tmp2, Ice%u_ice, Ice%v_ice,                      &
                         Ice%sig11, Ice%sig22, Ice%sig12, Ice%u_ocn, Ice%v_ocn,                         &
                         ice_avg(Ice%flux_u_top_bgrid(isc:iec,jsc:jec,:),Ice%part_size_uv(isc:iec,jsc:jec,:) ),  &
                         ice_avg(Ice%flux_v_top_bgrid(isc:iec,jsc:jec,:),Ice%part_size_uv(isc:iec,jsc:jec,:) ),  &
                         Ice%sea_lev, fx_wat, fy_wat, fx_ice, fy_ice, fx_cor, fy_cor)
    endif
    call mpp_clock_end(iceClocka)

    call mpp_clock_begin(iceClockb)
    if(reproduce_siena_201303) then
       call mpp_update_domains(Ice%u_ice, Ice%v_ice, Domain)
       if(tripolar_grid) then
          call cut_check('u_ice', Ice%u_ice) ! these calls fix round off differences
          call cut_check('v_ice', Ice%v_ice) ! in northernmost velocities over the fold
          call mpp_update_domains(Ice%u_ice, Ice%v_ice, Domain)
       endif
    else
       call mpp_update_domains(Ice%u_ice, Ice%v_ice, Domain, gridtype=BGRID_NE)
    endif
    call mpp_clock_end(iceClockb)

    call mpp_clock_begin(iceClockc)
    !
    ! Dynamics diagnostics
    !
    if (id_fax>0) &
         sent = send_data(id_fax, all_avg(Ice%flux_u_top_bgrid(isc:iec,jsc:jec,:),Ice%part_size_uv), Ice%Time)
    if (id_fay>0) &
         sent = send_data(id_fay, all_avg(Ice%flux_v_top_bgrid(isc:iec,jsc:jec,:),Ice%part_size_uv), Ice%Time)
    if (id_fix>0) sent = send_data(id_fix, fx_ice, Ice%Time)
    if (id_fiy>0) sent = send_data(id_fiy, fy_ice, Ice%Time)
    if (id_fcx>0) sent = send_data(id_fcx, fx_cor, Ice%Time)
    if (id_fcy>0) sent = send_data(id_fcy, fy_cor, Ice%Time)
    if (id_fwx>0) sent = send_data(id_fwx, -fx_wat, Ice%Time) ! water force on ice
    if (id_fwy>0) sent = send_data(id_fwy, -fy_wat, Ice%Time) ! ...= -ice on water

    if (id_strna>0) &
         sent = send_data(id_strna, strain_angle(Ice%u_ice,Ice%v_ice), Ice%Time, mask=Ice%mask)
    if (id_sigi>0)  sent = send_data(id_sigi, sigI(ice_avg(Ice%h_ice(isc:iec,jsc:jec,:),          &
                    Ice%part_size(isc:iec,jsc:jec,:)),1-Ice%part_size(isc:iec,jsc:jec,1),         &
                    Ice%sig11(isc:iec,jsc:jec),Ice%sig22(isc:iec,jsc:jec),Ice%sig12(isc:iec,jsc:jec)), &
                    Ice%Time, mask=Ice%mask)
    if (id_sigii>0) sent = send_data(id_sigii, sigII(ice_avg(Ice%h_ice(isc:iec,jsc:jec,:),        &
                    Ice%part_size(isc:iec,jsc:jec,:)),1-Ice%part_size(isc:iec,jsc:jec,1),         &
                    Ice%sig11(isc:iec,jsc:jec),Ice%sig22(isc:iec,jsc:jec),Ice%sig12(isc:iec,jsc:jec)), &
                    Ice%Time, mask=Ice%mask)
    if (id_stren>0) sent = send_data(id_stren, ice_strength(ice_avg(Ice%h_ice(isc:iec,jsc:jec,:), &
                    Ice%part_size(isc:iec,jsc:jec,:)), 1-Ice%part_size(isc:iec,jsc:jec,1)), Ice%Time, mask=Ice%mask)

    if (id_ui>0) sent = send_data(id_ui, Ice%u_ice(isc:iec,jsc:jec), Ice%Time)
    if (id_vi>0) sent = send_data(id_vi, Ice%v_ice(isc:iec,jsc:jec), Ice%Time)

    do k=2,km
       do j = jsc, jec
          do i = isc, iec
             Ice%flux_u_top_bgrid(i,j,k) = fx_wat(i,j)  ! stress of ice on ocean
             Ice%flux_v_top_bgrid(i,j,k) = fy_wat(i,j)  !
          enddo
       enddo
    end do
    call mpp_clock_end(iceClockc)
    call mpp_clock_end(iceClock4)

    !
    ! Thermodynamics
    !
    call mpp_clock_begin(iceClock5)
    if (id_frazil>0) sent = send_data(id_frazil, Ice%frazil/dt_slow,  Ice%Time, mask=Ice%mask)
    snow_to_ice = 0
    bsnk        = 0

    hi_change  = all_avg(Ice%h_ice(isc:iec,jsc:jec,:), Ice%part_size(isc:iec,jsc:jec,:))
    h2o_change = all_avg(DS*Ice%h_snow(isc:iec,jsc:jec,:)+DI*Ice%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))

    ! get observed ice thickness for ice restoring, if calculating qflux
    if (do_ice_restore) &
         call get_sea_surface(Ice%Time, dum1, Obs_cn_ice, Obs_h_ice)
    do k=2,km
       do j=jsc, jec
          do i=isc, iec 
             if (cell_area(i,j) > 0 .and.Ice%h_ice(i,j,k) > 0) then
                ! reshape the ice based on fluxes

                h2o_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
                call ice3lay_resize(Ice%h_snow(i,j,k), Ice%h_ice(i,j,k),               &
                                    Ice%t_ice1(i,j,k), Ice%t_ice2(i,j,k),              &
                                    Ice%fprec_top(i,j,k) *dt_slow, 0.0,                &
                                    Ice%flux_q_top(i,j,k)*dt_slow,                     &
                                    Ice%tmelt (i,j,k), Ice%bmelt(i,j,k),               &
                                    -MU_TS*Ice%s_surf(i,j),                            &
                                    heat_to_ocn, h2o_to_ocn, h2o_from_ocn,             &
                                    snow_to_ice(i,j,k), bablt                          )

                ! modify above-ice to under-ice fluxes for passing to ocean
                Ice%flux_q_top (i,j,k) = h2o_from_ocn/dt_slow ! no ice, evaporation left
                Ice%flux_lh_top(i,j,k) = hlv*Ice%flux_q_top(i,j,k)
                Ice%flux_lw_top(i,j,k) = 0.0
                Ice%flux_t_top (i,j,k) = Ice%bheat(i,j)-heat_to_ocn/dt_slow
                Ice%flux_sw_vis_dif_top(i,j,k) = (Ice%flux_sw_vis_dir_top(i,j,k)+      &
                      Ice%flux_sw_vis_dif_top(i,j,k)+Ice%flux_sw_nir_dir_top(i,j,k)+   &
                      Ice%flux_sw_nir_dif_top(i,j,k))*Ice%pen(i,j,k)*Ice%trn(i,j,k)
                Ice%flux_sw_nir_dir_top(i,j,k) = 0.0
                Ice%flux_sw_nir_dif_top(i,j,k) = 0.0
                Ice%flux_sw_vis_dir_top(i,j,k) = 0.0
                Ice%fprec_top  (i,j,k) = 0.0
                Ice%lprec_top  (i,j,k) = Ice%lprec_top(i,j,k) + h2o_to_ocn/dt_slow

                bsnk(i,j) = bsnk(i,j) - Ice%part_size(i,j,k)*bablt ! bot. melt. ablation

             end if
             !
             ! absorb frazil in thinest ice partition available
             !
             if (Ice%frazil(i,j)>0.and.Ice%part_size(i,j,1)+Ice%part_size(i,j,k)>0) then
                Ice%h_snow   (i,j,k) = Ice%h_snow   (i,j,k) * Ice%part_size(i,j,k)
                Ice%h_ice    (i,j,k) = Ice%h_ice    (i,j,k) * Ice%part_size(i,j,k)
                Ice%t_surf   (i,j,k) = (Ice%t_surf(i,j,k)    *Ice%part_size(i,j,k) &
                                     + (Tfreeze-MU_TS*Ice%s_surf(i,j))*Ice%part_size(i,j,1))
                Ice%part_size(i,j,k) = Ice%part_size(i,j,k) + Ice%part_size(i,j,1)
                Ice%part_size(i,j,1) = 0.0
                Ice%h_snow   (i,j,k) = Ice%h_snow   (i,j,k) / Ice%part_size(i,j,k)
                Ice%h_ice    (i,j,k) = Ice%h_ice    (i,j,k) / Ice%part_size(i,j,k)
                Ice%t_surf   (i,j,k) = Ice%t_surf   (i,j,k) / Ice%part_size(i,j,k)

                call ice3lay_resize(Ice%h_snow(i,j,k), Ice%h_ice (i,j,k),                &
                                    Ice%t_ice1(i,j,k), Ice%t_ice2(i,j,k), 0.0,           &
                                    Ice%frazil(i,j)/Ice%part_size(i,j,k), 0.0, 0.0, 0.0, &
                                    -MU_TS*Ice%s_surf(i,j),                              &
                                    heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic)
                Ice%frazil(i,j) = 0.0;
                !
                ! spread frazil salinification over all partitions
                !
                Ice%lprec_top  (i,j,:) = Ice%lprec_top(i,j,:) + h2o_to_ocn*Ice%part_size(i,j,k)/dt_slow
             end if

          end do  ! i-loop
       end do   ! j-loop
    end do   ! k-loop
    call mpp_clock_end(iceClock5)

    !
    ! Calculate QFLUX's from (1) restoring to obs and (2) limiting total ice.
    !
    call mpp_clock_begin(iceClock6)
    if (do_ice_restore .or. do_ice_limit) then
       do j = jsc, jec
          do i = isc, iec
             Ice % qflx_lim_ice(i,j) = 0.0
             Ice % qflx_res_ice(i,j) = 0.0
          enddo
       enddo

       do j=jsc, jec
          do i=isc, iec
             heat_res_ice   = 0.0
             heat_limit_ice = 0.0
             !
             ! calculate enthalpy
             !
             if (slab_ice) then
                e2m(2) = Ice%h_ice(i,j,2)*DI*LI
             else
                do k=2,km
                   if ((Ice%part_size(i,j,k)>0.0.and.Ice%h_ice(i,j,k)>0.0)) then
                      e2m(k) = e_to_melt(Ice%h_snow(i,j,k), Ice%h_ice(i,j,k)/2, &
                               Ice%t_ice1(i,j,k), Ice%h_ice(i,j,k)/2,           &
                               Ice%t_ice2(i,j,k)) * Ice%part_size(i,j,k)
                   else
                      e2m(k) = 0.0
                   end if
                end do
             end if
             !
             ! calculate heat needed to constrain ice enthalpy
             !
             if (do_ice_restore) then
                ! TK Mod: restore to observed enthalpy (no change for slab, but for
                !         sis ice, results in restoring toward thickness * concentration      

                ! TK Mod for test 10/18/02
                !   Concentration is 1.0 where there is ice for slab case.
                !   The input field Obs_cn_ice may have values less than 1,
                !     so put in if test...

                if (slab_ice) then
                   heat_res_ice = -(LI*DI*Obs_h_ice(i,j)-sum(e2m)) &
                                  *dt_slow/(86400*ice_restore_timescale)
                else                   
                   heat_res_ice = -(LI*DI*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-sum(e2m)) &
                                  *dt_slow/(86400*ice_restore_timescale)
                end if

                !        heat_res_ice = -(LI*DI*Obs_h_ice(i,j)-sum(e2m)) &
                !                        *dt_slow/(86400*ice_restore_timescale)

             end if

             if (do_ice_limit .and. (sum(e2m) > max_ice_limit*DI*LI)) then
                heat_limit_ice = sum(e2m)-LI*DI*max_ice_limit
                ! should we "heat_ice_res = 0.0" ?
             end if

             !
             ! apply constraining heat to ice
             !
             tot_heat = heat_res_ice+heat_limit_ice
             if (slab_ice) Ice%h_ice (i,j,2) = Ice%h_ice (i,j,2) - tot_heat/(DI*LI)

             if (.not. slab_ice .and. (tot_heat>0.0)) then  ! add like ocean-ice heat
                do k=2,km
                   if (e2m(k) > 0.0) then
                      heating = tot_heat/sum(Ice%part_size(i,j,k:km))
                      if (heating*Ice%part_size(i,j,k) > e2m(k)) then ! cat. melts away
                         Ice%h_ice (i,j,k) = 0.0
                         Ice%h_snow(i,j,k) = 0.0
                         tot_heat = tot_heat - e2m(k)
                      else
                         h2o_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
                         call ice3lay_resize(Ice%h_snow(i,j,k), Ice%h_ice(i,j,k),   &
                                             Ice%t_ice1(i,j,k), Ice%t_ice2(i,j,k),  &
                                             0.0, 0.0, 0.0, 0.0, heating,           &
                                             -MU_TS*Ice%s_surf(i,j),                &
                                             heat_to_ocn, h2o_to_ocn, h2o_from_ocn, &
                                             snow_to_ice(i,j,k), bablt              )
                         tot_heat = tot_heat - heating*Ice%part_size(i,j,k)
                      end if
                   end if
                end do
             end if

             tot_heat = heat_res_ice+heat_limit_ice
             if (.not. slab_ice .and. (tot_heat<0.0)) then ! add like frazil
                do k=2,km
                   if (Ice%part_size(i,j,1)+Ice%part_size(i,j,k)>0) exit
                end do
                ! k is thinnest ice partition that can recieve frazil
                Ice%h_snow   (i,j,k) = Ice%h_snow   (i,j,k) * Ice%part_size(i,j,k)
                Ice%h_ice    (i,j,k) = Ice%h_ice    (i,j,k) * Ice%part_size(i,j,k)
                Ice%t_surf   (i,j,k) = (Ice%t_surf(i,j,k)    *Ice%part_size(i,j,k) &
                                     +(Tfreeze-MU_TS*Ice%s_surf(i,j))*Ice%part_size(i,j,1))
                Ice%part_size(i,j,k) = Ice%part_size(i,j,k) + Ice%part_size(i,j,1)
                Ice%part_size(i,j,1) = 0.0
                Ice%h_snow   (i,j,k) = Ice%h_snow   (i,j,k) / Ice%part_size(i,j,k)
                Ice%h_ice    (i,j,k) = Ice%h_ice    (i,j,k) / Ice%part_size(i,j,k)
                Ice%t_surf   (i,j,k) = Ice%t_surf   (i,j,k) / Ice%part_size(i,j,k)

                call ice3lay_resize(Ice%h_snow(i,j,k), Ice%h_ice (i,j,k),          &
                                    Ice%t_ice1(i,j,k), Ice%t_ice2(i,j,k), 0.0,     &
                                    -tot_heat/Ice%part_size(i,j,k), 0.0, 0.0, 0.0, &
                                    -MU_TS*Ice%s_surf(i,j),                        &
                                    heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic)
             end if

             ! Convert constraining heat from energy (J/m^2) to flux (W/m^2)
             Ice%qflx_lim_ice(i,j) = heat_limit_ice / dt_slow
             Ice%qflx_res_ice(i,j) = heat_res_ice / dt_slow
             !
             ! Check for energy conservation
             !
             if (slab_ice) then
                e2m(2) = e2m(2) - Ice%h_ice(i,j,2)*DI*LI
             else
                do k=2,km
                   if (Ice%part_size(i,j,k)>0.0.and.Ice%h_ice(i,j,k)>0.0)                &
                        e2m(k) = e2m(k)-e_to_melt(Ice%h_snow(i,j,k), Ice%h_ice(i,j,k)/2, &
                                 Ice%t_ice1(i,j,k), Ice%h_ice(i,j,k)/2,                  &
                                 Ice%t_ice2(i,j,k)) * Ice%part_size(i,j,k)
                end do
             end if
             !     if (abs(sum(e2m) - heat_res_ice - heat_limit_ice)>DI*LI*1e-3) &
                  !       print *, 'QFLUX conservation error at', i, j, 'heat2ice=',  &
             !                 tot_heat, 'melted=', sum(e2m), 'h*part_size=',    &
                  !                 Ice%h_ice(i,j,:)*Ice%part_size(i,j,:)

          end do
       end do
    end if
    call mpp_clock_end(iceClock6)
    !
    ! Salt fluxes to ocean
    !
    call mpp_clock_begin(iceClock8)
    hi_change  = all_avg(Ice%h_ice(isc:iec,jsc:jec,:), Ice%part_size(isc:iec,jsc:jec,:))-hi_change
    Ice%flux_salt = ice_bulk_salin*DI*hi_change/dt_slow
    if (id_saltf>0)  sent = send_data(id_saltf, Ice%flux_salt, Ice%Time, mask=Ice%mask)

    h2o_change = all_avg(DS*Ice%h_snow(isc:iec,jsc:jec,:)+DI*Ice%h_ice(isc:iec,jsc:jec,:), &
                 Ice%part_size(isc:iec,jsc:jec,:))-h2o_change

    if (id_lsnk>0) then
       x = h2o_change*864e2*365/dt_slow
       do j = jsc, jec
          do i = isc, iec
             if(x(i,j)>0.0) x(i,j) = 0.0
          enddo
       enddo
       sent = send_data(id_lsnk, x, Ice%Time, mask=Ice%mask)
    end if
    if (id_lsrc>0) then
       x = h2o_change*864e2*365/dt_slow
       do j = jsc, jec
          do i = isc, iec
             if(x(i,j)<0.0) x(i,j) = 0.0
          enddo
       enddo
       sent = send_data(id_lsrc,  x, Ice%Time, mask=Ice%mask)
    end if
    if (id_bsnk>0)  sent = send_data(id_bsnk, bsnk*864e2*365/dt_slow, Ice%Time, mask=Ice%mask)
    if (id_tmelt>0) sent = send_data(id_tmelt, ice_avg(Ice%tmelt,Ice%part_size(isc:iec,jsc:jec,:))/dt_slow,   &
                           Ice%Time, mask=Ice%mask)
    if (id_bmelt>0) sent = send_data(id_bmelt, ice_avg(Ice%bmelt,Ice%part_size(isc:iec,jsc:jec,:))/dt_slow,   &
                           Ice%Time, mask=Ice%mask)
    if (id_sn2ic>0) sent = send_data(id_sn2ic, all_avg(snow_to_ice,Ice%part_size(isc:iec,jsc:jec,:))/dt_slow, &
                           Ice%Time, mask=Ice%mask                     )
    if (id_qflim>0) sent = send_data(id_qflim, Ice%qflx_lim_ice, Ice%Time, mask=Ice%mask)
    if (id_qfres>0) sent = send_data(id_qfres, Ice%qflx_res_ice, Ice%Time, mask=Ice%mask)
    !
    ! Ice transport ... all ocean fluxes have been calculated by now
    !
    h2o_change = all_avg(DS*Ice%h_snow(isc:iec,jsc:jec,:)+DI*Ice%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))
    call transport(Ice)

    x = all_avg(DS*Ice%h_snow(isc:iec,jsc:jec,:)+DI*Ice%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))
    if (id_mi>0) sent = send_data(id_mi, x, Ice%Time, mask = Ice%mask)
    
    do j = jsc, jec
       do i = isc, iec
         Ice%mi(i,j) = x(i,j) 
       enddo
    enddo

    if (do_icebergs) call icebergs_incr_mass(Ice%icebergs, x) ! Add icebergs mass in kg/m^2
    if (id_mib>0) sent = send_data(id_mib, x, Ice%Time, mask = Ice%mask) ! Diagnose total mass
    if (id_slp>0) sent = send_data(id_slp, Ice%p_surf, Ice%Time, mask=Ice%mask)
    do j = jsc, jec
       do i = isc, iec
          if (slp2ocean) then
            Ice%p_surf(i,j) = Ice%p_surf(i,j)-1e5
          else
            Ice%p_surf(i,j) = 0.0
          endif
          Ice%p_surf(i,j) = Ice%p_surf(i,j) + grav*x(i,j)
       enddo
    enddo
    h2o_change = x-h2o_change

    if (spec_ice) then                  ! over-write changes with specifications
         call get_sea_surface(Ice%Time, Ice%t_surf(:,:,1), Ice%part_size(isc:iec,jsc:jec,:), Ice%h_ice (isc:iec,jsc:jec,2))
         call mpp_update_domains(Ice%part_size, Domain)
    endif
    do j = jsd, jed
       do i = isd, ied
          Ice%part_size   (i,j,1) = 1.0
       enddo
    enddo

    do j = jsc, jec
       do i = isc, iec
          Ice%part_size_uv(i,j,1) = 1.0
       enddo
    enddo

    do k=2,km
       do j = jsd, jed
          do i = isd, ied
             Ice%part_size(i,j,1) = Ice%part_size   (i,j,1) - Ice%part_size   (i,j,k)
          enddo
       enddo
    enddo
    do k=2,km
       call t_to_uv(Ice%part_size(:,:,k),Ice%part_size_uv(:,:,k))
       Ice%part_size_uv(:,:,1) = Ice%part_size_uv(:,:,1) - Ice%part_size_uv(:,:,k)
    end do
    call mpp_clock_end(iceClock8)
    !
    ! Thermodynamics diagnostics
    !
    call mpp_clock_begin(iceClock9)
    if (id_cn>0) sent = send_data(id_cn, Ice%part_size(isc:iec,jsc:jec,2:km), Ice%Time, mask=spread(Ice%mask,3,km-1)       )

    ! TK Mod: 10/18/02
    !  if (id_obs_cn>0) sent = send_data(id_obs_cn, Obs_cn_ice(:,:,2), Ice%Time, &
         !                                       mask=Ice%mask       )

    if (id_ext>0) sent = send_data(id_ext, ext(Ice%part_size(isc:iec,jsc:jec,1)), Ice%Time, mask=Ice%mask)
    if (id_wnd>0) sent = send_data(id_wnd, Ice%wnd(isc:iec,jsc:jec,1), Ice%Time, mask=Ice%mask)

    if (id_hs>0) sent  = send_data(id_hs, ice_avg(Ice%h_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)), &
                         Ice%Time, mask=Ice%mask)
    if (id_hi>0) sent  = send_data(id_hi, ice_avg(Ice%h_ice(isc:iec,jsc:jec,:) ,Ice%part_size(isc:iec,jsc:jec,:)), &
                         Ice%Time, mask=Ice%mask)

    ! TK Mod: 10/18/02: (commented out...does not compile yet... add later
    !  if (id_obs_hi>0) &
         !    sent = send_data(id_obs_hi, ice_avg(Obs_h_ice,Ice%part_size), Ice%Time, &
    !                     mask=Ice%mask)

    if (id_ts>0) sent = send_data(id_ts, ice_avg(Ice%t_surf-Tfreeze,Ice%part_size(isc:iec,jsc:jec,:)), Ice%Time, mask=Ice%mask)
    if (id_t1>0) sent = send_data(id_t1, ice_avg(Ice%t_ice1(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)), &
                 Ice%Time, mask=Ice%mask)
    if (id_t2>0) sent = send_data(id_t2, ice_avg(Ice%t_ice2(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)), &
                 Ice%Time, mask=Ice%mask)
    if (id_xprt>0) sent = send_data(id_xprt,  h2o_change*864e2*365/dt_slow, Ice%Time, mask=Ice%mask)
    if (id_e2m>0) then
       x = 0.0
       do k = 2, km
          do j = jsc, jec
             do i=isc, iec
                if (Ice%h_ice(i,j,k)>0.0)                                                 &
                     x(i,j) = x(i,j) + Ice%part_size(i,j,k)*e_to_melt(Ice%h_snow(i,j,k),  &
                              Ice%h_ice(i,j,k)/2, Ice%t_ice1(i,j,k), Ice%h_ice(i,j,k)/2, Ice%t_ice2(i,j,k) )
             end do
          end do
       end do
       sent = send_data(id_e2m,  x, Ice%Time, mask=Ice%mask)
    end if

    call ice_top_to_ice_bottom(Ice, part_save, part_save_uv)
    !
    ! conservation checks:  bottom fluxes and final
    !
    if (conservation_check) then
       h2o(3)  = h2o(3) + dt_slow*sum(cell_area *(Ice%lprec+Ice%fprec-Ice%flux_q+Ice%runoff+Ice%calving),&
                 cell_area > 0 )
       heat(3) = heat(3) + tot_frazil + dt_slow*sum(cell_area*(Ice%flux_sw_vis_dir+Ice%flux_sw_vis_dif+  &
                 Ice%flux_sw_nir_dir+Ice%flux_sw_nir_dif+Ice%flux_lw- &
                 Ice%flux_t-Ice%flux_lh -LI*(Ice%fprec+Ice%calving)), cell_area > 0 )
       h2o(4)  = sum(cell_area*all_avg(DI*Ice%h_ice(isc:iec,jsc:jec,:)+ &
                 DS*Ice%h_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)))
       heat(4) = 0.0

       do k=2,km
          do j=jsc, jec
             do i=isc, iec
                if ((Ice%part_size(i,j,k)>0.0.and.Ice%h_ice(i,j,k)>0.0)) then
                   if (slab_ice) then
                      heat(4) = heat(4) - cell_area(i,j) * Ice%part_size(i,j,k)*Ice%h_ice(i,j,2)*DI*LI
                   else
                      heat(4) = heat(4) - cell_area(i,j) * Ice%part_size(i,j,k)*e_to_melt(Ice%h_snow(i,j,k), &
                                Ice%h_ice(i,j,k)/2, Ice%t_ice1(i,j,k), Ice%h_ice(i,j,k)/2, Ice%t_ice2(i,j,k))
                   end if
                end if
             end do
          end do
       end do
    end if

    call get_date(Ice%Time, iyr, imon, iday, ihr, imin, isec)
    call get_time(Ice%Time-set_date(iyr,1,1,0,0,0),isec,iday)
    if(verbose)  call ice_line(iyr, iday+1, isec, Ice%part_size(isc:iec,jsc:jec,:), &
               Ice%t_surf(:,:,1)-Tfreeze) 

    do j=jsc, jec
       do i=isc, iec 
          if (Ice%mask(i,j).and.(abs(sum(Ice%part_size(i,j,:))-1.0)>1e-2)) &
               print *,'ICE%PART_SIZE=',Ice%part_size(i,j,:), 'DOES NOT SUM TO 1 AT', &
               geo_lon(i,j), geo_lat(i,j)
       end do
    end do
    call mpp_clock_end(iceClock9)


  end subroutine update_ice_model_slow

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_redistribute - a simple ice redistribution scheme from Igor Polyakov     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_redistribute(cn, hs, hi, t1, t2)
!!$    real, intent(inout), dimension(2:km  )           :: cn, hs, hi, t1, t2
    real, intent(inout), dimension(2:)           :: cn, hs, hi, t1, t2

    real    :: cw                                 ! open water concentration
    integer :: k

    cw = 1-sum(cn)
    if (cw<0) cn(2) = cw + cn(2) ! open water has been eliminated by convergence

    do k=2,km-1
       if (cn(k)<0) then
          cn(k+1) = cn(k+1)+cn(k); cn(k) = 0 ! pass concentration deficit up to
          hs(k+1) = hs(k+1)+hs(k); hs(k) = 0 ! next thicker category
          hi(k+1) = hi(k+1)+hi(k); hi(k) = 0
          t1(k+1) = t1(k+1)+t1(k); t1(k) = 0 ! NOTE:  here between the thm_pack and
          t2(k+1) = t2(k+1)+t2(k); t2(k) = 0 ! thm_unpack calls, all quantities are
       endif                                ! extensive, so we add instead of
    end do                                 ! averaging

    do k=2,km-1
       if (hi(k)>hlim(k)*cn(k)) then
          cn(k+1) = cn(k+1)+cn(k); cn(k) = 0 ! upper thickness limit exceeded
          hs(k+1) = hs(k+1)+hs(k); hs(k) = 0 ! move ice up to next thicker category
          hi(k+1) = hi(k+1)+hi(k); hi(k) = 0
          t1(k+1) = t1(k+1)+t1(k); t1(k) = 0
          t2(k+1) = t2(k+1)+t2(k); t2(k) = 0
       endif
    end do

    do k=km,3,-1
       if (hi(k)<hlim(k-1)*cn(k)) then
          cn(k-1) = cn(k-1)+cn(k); cn(k) = 0  ! lower thickness limit exceeded;
          hs(k-1) = hs(k-1)+hs(k); hs(k) = 0  ! move ice down to thinner category
          hi(k-1) = hi(k-1)+hi(k); hi(k) = 0
          t1(k-1) = t1(k-1)+t1(k); t1(k) = 0
          t2(k-1) = t2(k-1)+t2(k); t2(k) = 0
       endif
    end do

  end subroutine ice_redistribute

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! transport - do ice transport and thickness class redistribution              !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine transport (Ice)
    type (ice_data_type), intent(inout) :: Ice
    integer :: i, j, k
    logical :: sent
    real, dimension(size(Ice%t_surf,1),size(Ice%t_surf,2)) :: uf0, uf, vf0, vf
    real, dimension(isd:ied,jsd:jed) :: uc, vc ! Local variables, C-grid transporting velocities
    real, dimension(isd:ied,jsd:jed) :: tmp1 ! Local variables, 2D ice concentration
    real, dimension(isd:ied,jsd:jed) :: ustar, vstar ! Local variables, C-grid transporting velocities
    real, dimension(isd:ied,jsd:jed) :: ustaro, vstaro, ustarv, vstarv ! Local variables, C-grid transporting velocities
    real :: u_visc, u_ocn, cnn, grad_eta ! Variables for channel parameterization

    if (slab_ice) then
       call slab_ice_advect(Ice%u_ice, Ice%v_ice, Ice%h_ice(:,:,2), 4.0)
       call mpp_update_domains(Ice%h_ice(:,:,2), Domain)
       do j = jsd, jed
          do i = isd, ied
             if(Ice%h_ice(i,j,2) > 0.0) then
                Ice%part_size(i,j,2) = 1.0
             else 
                Ice%part_size(i,j,2) = 0.0
             endif
          enddo
       enddo
       return;
    end if

    call thm_pack(Ice%part_size(isc:iec,jsc:jec,2:km), Ice%h_snow(isc:iec,jsc:jec,:), Ice%h_ice(isc:iec,jsc:jec,:), &
                  Ice%t_ice1(isc:iec,jsc:jec,:), Ice%t_ice2(isc:iec,jsc:jec,:) )

    call mpp_update_domains(Ice%part_size, Domain) ! cannot be combined with updates below
    call mpp_update_domains(Ice%h_snow, Domain, complete=.false.)
    call mpp_update_domains(Ice%h_ice, Domain, complete=.false.)
    call mpp_update_domains(Ice%t_ice1, Domain, complete=.false.)
    call mpp_update_domains(Ice%t_ice2, Domain, complete=.true.)

    ! Move transporting flow to the C-grid (from the B-grid)
    ! Note: this block of code was moved here from within s/r ice_advect
    ! where the same calculations and mpp_updates were repeated for each
    ! class and variable. -AJA
    uc = 0.0; vc = 0.0
    do j = jsc, jec
      do i = isc, iec
        uc(i,j) = 0.5 * ( Ice%u_ice(i,j-1) + Ice%u_ice(i,j) )
        vc(i,j) = 0.5 * ( Ice%v_ice(i-1,j) + Ice%v_ice(i,j) )
      enddo
    enddo
    if (channel_viscosity>0.) then
    ! This block of code is a parameterization of either (or both)
    ! i) the pressure driven oceanic flow in a narrow channel, or
    ! ii) pressure driven flow of ice itself.
    ! The latter is a speculative but both are missing due to
    ! masking of velocities to zero in a single-cell wide channel.
      tmp1=1.-max(1.-sum(Ice%part_size(:,:,2:km),dim=3),0.0)
      ustar(:,:)=0.; vstar(:,:)=0.
      ustaro(:,:)=0.; vstaro(:,:)=0.
      ustarv(:,:)=0.; vstarv(:,:)=0.
      do j = jsc, jec
        do i = isc, iec
          if (Ice%u_ice(i,j).ne.0..and.Ice%vmask(i,j)==0.) stop 'Ooops' ! Debug new vmask
          if (Ice%v_ice(i,j).ne.0..and.Ice%vmask(i,j)==0.) stop 'Ooops' ! Debug new vmask
          if (uc(i,j)==0. & ! this is a redundant test due to following line
             .and. Ice%vmask(i,j)+Ice%vmask(i,j-1)==0. &  ! =0 => no vels
             .and. wett(i,j)*wett(i+1,j)>0.) then ! >0 => open for transport
               grad_eta=(Ice%sea_lev(i+1,j)-Ice%sea_lev(i,j))    & ! delta_i eta
                        /(0.5*(dxv(i,j)+dxv(i,j-1)))               ! /dx
               u_visc=-ssh_gravity*((dte(i,j)*dte(i,j))/(12.*channel_viscosity)) & ! -g*dy^2/(12*visc)
                        *grad_eta                                                  ! d/dx eta
               u_ocn=sqrt( ssh_gravity*dte(i,j)*abs(grad_eta)/(36.*smag_ocn) ) ! Magnitude of ocean current
               u_ocn=sign(u_ocn, -grad_eta) ! Direct down the ssh gradient
               cnn=max(tmp1(i,j),tmp1(i+1,j))**2. ! Use the larger concentration
               uc(i,j)=cnn*u_visc+(1.-cnn)*u_ocn
               ! Limit flow to be stable for fully divergent flow
               if (uc(i,j)>0.) then
                 uc(i,j)=min( uc(i,j), chan_cfl_limit*dxt(i,j)/dt_adv)
               else
                 uc(i,j)=max( uc(i,j),(-1*chan_cfl_limit)*dxt(i+1,j)/dt_adv)
               endif
               if (id_ustar>0) ustar(i,j)=uc(i,j)
               if (id_uocean>0) ustaro(i,j)=u_ocn
               if (id_uchan>0) ustarv(i,j)=u_visc
          endif
          if (vc(i,j)==0. & ! this is a redundant test due to following line
             .and. Ice%vmask(i,j)+Ice%vmask(i-1,j)==0. &  ! =0 => no vels
             .and. wett(i,j)*wett(i,j+1)>0.) then ! >0 => open for transport
               grad_eta=(Ice%sea_lev(i,j+1)-Ice%sea_lev(i,j))    & ! delta_i eta
                        /(0.5*(dyv(i,j)+dyv(i,j-1)))               ! /dy
               u_visc=-ssh_gravity*((dtn(i,j)*dtn(i,j))/(12.*channel_viscosity)) & ! -g*dy^2/(12*visc)
                        *grad_eta                                                  ! d/dx eta
               u_ocn=sqrt( ssh_gravity*dtn(i,j)*abs(grad_eta)/(36.*smag_ocn) ) ! Magnitude of ocean current
               u_ocn=sign(u_ocn, -grad_eta) ! Direct down the ssh gradient
               cnn=max(tmp1(i,j),tmp1(i,j+1))**2. ! Use the larger concentration
               vc(i,j)=cnn*u_visc+(1.-cnn)*u_ocn
               ! Limit flow to be stable for fully divergent flow
               if (vc(i,j)>0.) then
                 vc(i,j)=min( vc(i,j), chan_cfl_limit*dyt(i,j)/dt_adv)
               else
                 vc(i,j)=max( vc(i,j),(-1*chan_cfl_limit)*dyt(i,j+1)/dt_adv)
               endif
               if (id_vstar>0) vstar(i,j)=vc(i,j)
               if (id_vocean>0) vstaro(i,j)=u_ocn
               if (id_vchan>0) vstarv(i,j)=u_visc
          endif
        enddo
      enddo
    endif

    if(reproduce_siena_201303) then
       call mpp_update_domains(uc, vc, Domain)
    else
       call mpp_update_domains(uc, vc, Domain, gridtype=CGRID_NE)
    endif

    uf = 0.0; vf = 0.0
    uf0 = 0.0; vf0 = 0.0 !Must have a value in case of no subcycling
    do k=2,km
       call ice_advect(uc, vc, Ice%part_size(:,:,k))
       call ice_advect(uc, vc, Ice%h_snow   (:,:,k), uf0, vf0)
       uf = uf + DS*uf0; vf = vf + DS*vf0
       call ice_advect(uc, vc, Ice%h_ice    (:,:,k), uf0, vf0)
       uf = uf + DI*uf0; vf = vf + DI*vf0
       call ice_advect(uc, vc, Ice%t_ice1   (:,:,k))
       call ice_advect(uc, vc, Ice%t_ice2   (:,:,k))
    end do
    sent = send_data(id_ix_trans,  uf, Ice%Time)
    sent = send_data(id_iy_trans,  vf, Ice%Time)

    do j=jsc, jec
       do i=isc, iec
          if (sum(Ice%h_ice(i,j,:))>0)                                  &
               call ice_redistribute(Ice%part_size(i,j,2:km),              &
               Ice%h_snow(i,j,:), Ice%h_ice (i,j,:), &
               Ice%t_ice1(i,j,:), Ice%t_ice2(i,j,:))
          do k=2,km
             if (Ice%part_size(i,j,k)<1e-10) then
                Ice%h_ice (i,j,k) = 0           ! thm_unpack will zero other quantities
                Ice%t_surf(i,j,k) = Tfreeze-MU_TS*Ice%s_surf(i,j)
             end if
          end do
       end do
    end do

    call thm_unpack(Ice%part_size(isc:iec,jsc:jec,2:km), Ice%h_snow(isc:iec,jsc:jec,:), Ice%h_ice(isc:iec,jsc:jec,:), &
                    Ice%t_ice1(isc:iec,jsc:jec,:), Ice%t_ice2(isc:iec,jsc:jec,:) )

    call mpp_update_domains(Ice%part_size, Domain) ! cannot be combined with the two updates below
    call mpp_update_domains(Ice%h_snow, Domain, complete=.false.)
    call mpp_update_domains(Ice%h_ice, Domain, complete=.true.)

    if (id_ustar >0) sent = send_data(id_ustar , ustar(isc:iec,jsc:jec) , Ice%Time)
    if (id_vstar >0) sent = send_data(id_vstar , vstar(isc:iec,jsc:jec) , Ice%Time)
    if (id_vocean>0) sent = send_data(id_vocean, vstaro(isc:iec,jsc:jec) , Ice%Time)
    if (id_uocean>0) sent = send_data(id_uocean, ustaro(isc:iec,jsc:jec) , Ice%Time)
    if (id_vchan>0)  sent = send_data(id_vchan,  vstarv(isc:iec,jsc:jec) , Ice%Time)
    if (id_uchan>0)  sent = send_data(id_uchan,  ustarv(isc:iec,jsc:jec) , Ice%Time)

    return

  end subroutine transport

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ext - modeled ice is one if ice cover >= 0.15 otherwise zero                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function ext(cw)
    real, dimension(isc:,jsc:), intent(in) :: cw
    real,  dimension(isc:iec,jsc:jec)      :: ext
    integer :: i, j 

    do j = jsc, jec
       do i = isc, iec
          if(cw(i,j) < 0.85) then
             ext(i,j) = 1
          else
             ext(i,j) = 0
          endif
       enddo
    enddo

    return
  end function ext

subroutine ocn_ice_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_ice_boundary_type), intent(in) :: bnd_type
 integer ::   n, m, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(ocean_ice_boundary_type):: ', id, timestep
    write(outunit,100) 'ocn_ice_bnd_type%u        ',mpp_chksum(bnd_type%u        )
    write(outunit,100) 'ocn_ice_bnd_type%v        ',mpp_chksum(bnd_type%v        )
    write(outunit,100) 'ocn_ice_bnd_type%t        ',mpp_chksum(bnd_type%t        )
    write(outunit,100) 'ocn_ice_bnd_type%s        ',mpp_chksum(bnd_type%s        )
    write(outunit,100) 'ocn_ice_bnd_type%frazil   ',mpp_chksum(bnd_type%frazil   )
    write(outunit,100) 'ocn_ice_bnd_type%sea_level',mpp_chksum(bnd_type%sea_level)
!    write(outunit,100) 'ocn_ice_bnd_type%data     ',mpp_chksum(bnd_type%data     )
100 FORMAT("CHECKSUM::",A32," = ",Z20)

    do n = 1, bnd_type%fields%num_bcs  !{
       do m = 1, bnd_type%fields%bc(n)%num_fields  !{
          write(outunit,101) 'oibt%',trim(bnd_type%fields%bc(n)%name), &
               trim(bnd_type%fields%bc(n)%field(m)%name), &
               mpp_chksum(bnd_type%fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ocn_ice_bnd_type_chksum

subroutine atm_ice_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_ice_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(atmos_ice_boundary_type):: ', id, timestep
    write(outunit,100) 'atm_ice_bnd_type%u_flux          ',mpp_chksum(bnd_type%u_flux)          
    write(outunit,100) 'atm_ice_bnd_type%v_flux          ',mpp_chksum(bnd_type%v_flux)
    write(outunit,100) 'atm_ice_bnd_type%u_star          ',mpp_chksum(bnd_type%u_star)
    write(outunit,100) 'atm_ice_bnd_type%t_flux          ',mpp_chksum(bnd_type%t_flux)
    write(outunit,100) 'atm_ice_bnd_type%q_flux          ',mpp_chksum(bnd_type%q_flux)
    write(outunit,100) 'atm_ice_bnd_type%lw_flux         ',mpp_chksum(bnd_type%lw_flux)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dir ',mpp_chksum(bnd_type%sw_flux_vis_dir)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dif ',mpp_chksum(bnd_type%sw_flux_vis_dif)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dir ',mpp_chksum(bnd_type%sw_flux_nir_dir)
    write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dif ',mpp_chksum(bnd_type%sw_flux_nir_dif)
    write(outunit,100) 'atm_ice_bnd_type%lprec           ',mpp_chksum(bnd_type%lprec)
    write(outunit,100) 'atm_ice_bnd_type%fprec           ',mpp_chksum(bnd_type%fprec)
    write(outunit,100) 'atm_ice_bnd_type%dhdt            ',mpp_chksum(bnd_type%dhdt)
    write(outunit,100) 'atm_ice_bnd_type%dedt            ',mpp_chksum(bnd_type%dedt)
    write(outunit,100) 'atm_ice_bnd_type%drdt            ',mpp_chksum(bnd_type%drdt)
    write(outunit,100) 'atm_ice_bnd_type%coszen          ',mpp_chksum(bnd_type%coszen)
    write(outunit,100) 'atm_ice_bnd_type%p               ',mpp_chksum(bnd_type%p)
!    write(outunit,100) 'atm_ice_bnd_type%data            ',mpp_chksum(bnd_type%data)
100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine atm_ice_bnd_type_chksum

subroutine lnd_ice_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(land_ice_boundary_type):: ', id, timestep
    write(outunit,100) 'lnd_ice_bnd_type%runoff  ',mpp_chksum(bnd_type%runoff)
    write(outunit,100) 'lnd_ice_bnd_type%calving ',mpp_chksum(bnd_type%calving)
!    write(outunit,100) 'lnd_ice_bnd_type%data    ',mpp_chksum(bnd_type%data)
100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine lnd_ice_bnd_type_chksum

end module ice_model_mod
