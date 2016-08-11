
module atmosphere_mod

!-----------------------------------------------------------------------
!
!         interface for b-grid dynamical core and physics
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use mpp_mod, only: input_nml_file 
use bgrid_core_driver_mod, only: bgrid_dynam_type,       &
                                 bgrid_core_driver_init, &
                                 bgrid_core_driver,      &
                                 bgrid_core_time_diff,   &
                                 bgrid_core_driver_end,  &
                                 get_bottom_data,        &
                                 put_bottom_data

use   bgrid_prog_var_mod, only: prog_var_type, var_init

use      bgrid_horiz_mod, only: get_horiz_grid_size,        &
                                get_horiz_grid_bound, TGRID

use     time_manager_mod, only: time_type, get_time, operator(+)

use              fms_mod, only: file_exist, open_namelist_file, &
                                error_mesg, FATAL, WARNING,     &
                                check_nml_error, stdlog,        &
                                write_version_number,           &
                                mpp_pe, mpp_root_pe,            &
                                close_file, set_domain,         &
                                mpp_clock_id, mpp_clock_begin,  &
                                mpp_clock_end, MPP_CLOCK_SYNC,  &
                                CLOCK_SUBCOMPONENT, NOTE

use        bgrid_vert_mod, only: compute_height_bottom
use bgrid_change_grid_mod, only: change_grid, WIND_GRID, TEMP_GRID

use     bgrid_physics_mod, only: bgrid_physics_down, bgrid_physics_up,  &
                                 bgrid_physics_init, bgrid_physics_end, &
                                 surf_diff_type

use      mpp_domains_mod, only: domain2d

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index, NO_TRACER
use xgrid_mod,          only: grid_box_type

!-----------------------------------------------------------------------

implicit none
private

public  atmosphere_down,       atmosphere_up,       &
        atmosphere_init,       atmosphere_end,      &
        atmosphere_resolution, atmosphere_boundary, &
        atmosphere_cell_area,  atmosphere_restart,  &
        get_atmosphere_axes,   atmosphere_domain

public  get_bottom_mass,     get_bottom_wind

public  surf_diff_type,      get_stock_pe

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 19.0 2012/01/06 19:52:46 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,1]

   integer, dimension(2) :: physics_window = (/0,1/)

   namelist /atmosphere_nml/ physics_window

!-----------------------------------------------------------------------
!---- private data ----

type (bgrid_dynam_type), save :: Dynam
type    (prog_var_type), save :: Var, Var_dt
type        (time_type) :: Time_step_atmos

real                               :: dt_atmos
real,    dimension(:,:,:), pointer :: omega =>NULL()
integer, dimension(4)              :: atmos_axes

integer :: id_dynam, id_phys_down, id_phys_up
logical :: stock_warning_issued = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine atmosphere_down (Time,    frac_land,                         &
                             t_surf,  albedo, albedo_vis_dir,            &
                             albedo_nir_dir,  albedo_vis_dif,            &
                             albedo_nir_dif,  rough_mom,                 &
                             u_star,  b_star, q_star,                    &
                             dtau_du, dtau_dv, tau_x,   tau_y,           &
                             frac_open_sea,                              &
                             gust, coszen, flux_sw, flux_sw_dir,         &
                             flux_sw_dif,                                &
                         flux_sw_down_vis_dir,   flux_sw_down_vis_dif,   &
                         flux_sw_down_total_dir, flux_sw_down_total_dif, &
                             flux_sw_vis,                                &
                             flux_sw_vis_dir,                            &
                             flux_sw_vis_dif, flux_lw,                   &
                             Surf_diff                       )

!
!        Time = time at the current time level
!

  type(time_type),intent(in)    :: Time

  real, intent(in),    dimension(:,:) :: frac_land,                 &
                                         t_surf, albedo, rough_mom, &
                                   albedo_vis_dir, albedo_nir_dir,  &
                                   albedo_vis_dif, albedo_nir_dif,  &
                                   u_star, b_star, q_star,          &
                                   dtau_du, dtau_dv, frac_open_sea
  real, intent(inout), dimension(:,:) :: tau_x,  tau_y
  real, intent(out),   dimension(:,:) :: gust, coszen, flux_sw,    &
                                                      flux_sw_dir,    &
                                                      flux_sw_dif,    &
                                          flux_sw_down_vis_dir, &
                                         flux_sw_down_vis_dif, &
                                         flux_sw_down_total_dir,   &
                                         flux_sw_down_total_dif,   &
                                         flux_sw_vis,         &
                                         flux_sw_vis_dir,    &
                                         flux_sw_vis_dif, flux_lw
  type(surf_diff_type), intent(inout) :: Surf_diff

  type(time_type) :: Time_prev, Time_next
!-----------------------------------------------------------------------

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- dynamics -----

   call mpp_clock_begin (id_dynam)
   call bgrid_core_driver ( Time_next, Var, Var_dt, Dynam, omega )
   call mpp_clock_end (id_dynam)

!---- call physics -----

   call mpp_clock_begin (id_phys_down)
   call bgrid_physics_down (physics_window, dt_atmos,         &
                            Time_prev, Time, Time_next,       &
                            Dynam%Hgrid, Dynam%Vgrid, Dynam,  &
                            Var,     Var_dt,     frac_land,   &
                            albedo,  albedo_vis_dir, albedo_nir_dir,  &
                            albedo_vis_dif, albedo_nir_dif,  &
                            rough_mom,  t_surf,      &
                            u_star,  b_star, q_star,          &
                            dtau_du, dtau_dv, tau_x, tau_y,   &
                            flux_sw,    &
                            flux_sw_dir,    &
                            flux_sw_dif,    &
                            flux_sw_down_vis_dir,     &
                            flux_sw_down_vis_dif,     &
                            flux_sw_down_total_dir,  &
                            flux_sw_down_total_dif,  &
                            flux_sw_vis,  &
                            flux_sw_vis_dir,  &
                            flux_sw_vis_dif,  &
                            flux_lw, coszen, gust,   &
                            Surf_diff, frac_open_sea          )
   call mpp_clock_end (id_phys_down)

!-----------------------------------------------------------------------

 end subroutine atmosphere_down

!#######################################################################

 subroutine atmosphere_up (Time,  frac_land, Surf_diff, lprec, fprec, gust, &
                           u_star, b_star, q_star )

  type(time_type),intent(in)        :: Time
  real, intent(in),  dimension(:,:) :: frac_land
  type(surf_diff_type), intent(inout) :: Surf_diff
  real, intent(out), dimension(:,:) :: lprec,   fprec
  real, intent(inout), dimension(:,:) :: gust
  real, intent(in),  dimension(:,:) :: u_star, b_star, q_star

  type(time_type) :: Time_prev, Time_next
!-----------------------------------------------------------------------
!------ call physics up ------

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

     call mpp_clock_begin (id_phys_up)
     call bgrid_physics_up (physics_window, dt_atmos,        &
                            Time_prev, Time, Time_next,      &
                            Dynam%Hgrid, Dynam%Vgrid, Dynam, &
                            Var,        Var_dt,     omega,   &
                            frac_land,  Surf_diff,           &
                            lprec,      fprec,      gust,    &
                            u_star,     b_star,     q_star   )
     call mpp_clock_end (id_phys_up)

!------ time differencing and diagnostics -------

     call bgrid_core_time_diff ( omega, Time_next, Dynam, Var, Var_dt )

!-----------------------------------------------------------------------

 end subroutine atmosphere_up

!#######################################################################

 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff, Grid_box)

 type (time_type),     intent(in)    :: Time_init, Time, Time_step
 type(surf_diff_type), intent(inout) :: Surf_diff
 type(grid_box_type),  intent(inout) :: Grid_box

  integer :: unit, sec, ierr, io, logunit

!-----------------------------------------------------------------------
!----- read namelist -----

    if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=atmosphere_nml, iostat=io)
        ierr = check_nml_error(io,'atmosphere_nml')
#else
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read (unit, nml=atmosphere_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'atmosphere_nml')
        enddo
 10     call close_file (unit)
#endif
    endif

!----- write version and namelist to log file -----

    call write_version_number(version, tagname)
    logunit = stdlog()
    if ( mpp_pe() == mpp_root_pe() ) write (logunit, nml=atmosphere_nml)

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize dynamical core -----

   call bgrid_core_driver_init ( Time_init, Time, Time_step,    &
                                 Var, Var_dt, Dynam, atmos_axes )

!----- initialize storage needed for vert motion ----

    omega => var_init (Dynam%Hgrid, Dynam%Vgrid%nlev)

!----- initialize physics interface -----
!----- initialize domains for reading global physics data -----

    call set_domain ( Dynam%Hgrid%Tmp%Domain_nohalo )

    call bgrid_physics_init (atmos_axes, Time, Dynam%Hgrid, Dynam%Vgrid, Dynam, &
                             Var, Surf_diff)

!   ----- use entire grid as window ? -----

    if (physics_window(1) <= 0) physics_window(1) = Dynam%Hgrid%Tmp%ie-Dynam%Hgrid%Tmp%is+1
    if (physics_window(2) <= 0) physics_window(2) = Dynam%Hgrid%Tmp%je-Dynam%Hgrid%Tmp%js+1

!  --- initialize clocks for dynamics, physics_down and physics_up

    id_dynam     = mpp_clock_id ('BGRID: dynamical core',   &
                       flags=MPP_CLOCK_SYNC, grain=CLOCK_SUBCOMPONENT )
    id_phys_down = mpp_clock_id ('Physics_down',   &
                       flags=MPP_CLOCK_SYNC, grain=CLOCK_SUBCOMPONENT )
    id_phys_up   = mpp_clock_id ('Physics_up',   &
                       flags=MPP_CLOCK_SYNC, grain=CLOCK_SUBCOMPONENT )

!-----------------------------------------------------------------------

 end subroutine atmosphere_init

!#######################################################################

 subroutine atmosphere_end (Time, Grid_box)

 type (time_type), intent(in) :: Time
 type(grid_box_type),  intent(inout) :: Grid_box
 integer :: unit


    call bgrid_core_driver_end ( Var, Dynam )

!----- initialize domains for writing global physics data -----

    call set_domain ( Dynam%Hgrid%Tmp%Domain_nohalo )

    call bgrid_physics_end   (Time)                           

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

!#######################################################################
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_resolution (nlon, nlat, global)

  integer, intent(out)          :: nlon, nlat
  logical, intent(in), optional :: global

!---- return the size of the grid used for physics computations ----

    call get_horiz_grid_size (Dynam % Hgrid, TGRID, nlon, nlat, global)

 end subroutine atmosphere_resolution

!#######################################################################
 subroutine atmosphere_cell_area(area_out)
   real, dimension(:,:), intent(out) :: area_out
   integer :: j

   do j=Dynam%Hgrid%Tmp%js,Dynam%Hgrid%Tmp%je
     area_out(:,j-Dynam%Hgrid%Tmp%js+1) = Dynam%Hgrid%Tmp%area(j)
   enddo

 end subroutine atmosphere_cell_area
!#######################################################################
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_boundary (blon, blat, global)

    real,    intent(out)          :: blon(:,:), blat(:,:)
    logical, intent(in), optional :: global
    
    real, dimension(size(blon,1)) :: rlonb
    real, dimension(size(blat,2)) :: rlatb
    integer :: i, j
!----- return the longitudinal and latitudinal grid box edges ----------

    call get_horiz_grid_bound (Dynam % Hgrid, TGRID, rlonb, rlatb, global)
    do i=1,size(blon,1)
       blon(i,:) = rlonb(i)
    end do
    do j=1,size(blat,2)
       blat(:,j) = rlatb(j)
    end do

 end subroutine atmosphere_boundary

!#######################################################################
!    returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

 subroutine atmosphere_domain (Domain)
 type(domain2d), intent(inout) :: Domain

   Domain = Dynam % Hgrid % Tmp % Domain_nohalo

 end subroutine atmosphere_domain

!#######################################################################
!    returns the axis indices associated with the coupling grid

 subroutine get_atmosphere_axes ( axes )

   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----

     if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                           'get_atmosphere_axes in atmosphere_mod', &
                           'size of argument is incorrect', FATAL   )

     axes (1:size(axes(:))) = atmos_axes (1:size(axes(:)))

 end subroutine get_atmosphere_axes

!#######################################################################
! returns temp, tracers, pres, height at the lowest model level
!         and surface pressure and sea level pressure

 subroutine get_bottom_mass (t_bot, tr_bot, p_bot, z_bot, p_surf, slp)

  real, intent(out),   &
        dimension(Dynam%Hgrid%Tmp%is:Dynam%Hgrid%Tmp%ie,Dynam%Hgrid%Tmp%js:Dynam%Hgrid%Tmp%je)   &
        :: t_bot, p_bot, z_bot, p_surf, slp
  real, intent(out),   &
        dimension(Dynam%Hgrid%Tmp%is:Dynam%Hgrid%Tmp%ie,Dynam%Hgrid%Tmp%js:Dynam%Hgrid%Tmp%je,Var%ntrace) :: tr_bot

  real, dimension(Dynam%Hgrid%Tmp%is:Dynam%Hgrid%Tmp%ie,Dynam%Hgrid%Tmp%js:Dynam%Hgrid%Tmp%je) :: q_bot
  integer :: i, j, kb, sphum, n
  integer :: is,ie,js,je


  is = Dynam%Hgrid%Tmp%is;  ie = Dynam%Hgrid%Tmp%ie
  js = Dynam%Hgrid%Tmp%js;  je = Dynam%Hgrid%Tmp%je

  sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )

    do j = js, je
    do i = is, ie
      kb = Dynam%Masks%Tmp%kbot(i,j)
      t_bot(i,j) = Var % t(i,j,kb)
    enddo
    enddo

    if(sphum == NO_TRACER) then
      q_bot = 0.0
    else
      do j = js, je
      do i = is, ie
        kb = Dynam%Masks%Tmp%kbot(i,j)
        q_bot(i,j) = Var % r(i,j,kb,sphum)
      enddo
      enddo
    endif

    do n=1,size(tr_bot,3)
    do j = js, je
    do i = is, ie
      kb = Dynam%Masks%Tmp%kbot(i,j)
      tr_bot(i,j,n) = Var % r(i,j,kb,n)
    enddo
    enddo
    enddo

    p_surf = Var % ps (is:ie,js:je)

    slp = Var % ps (is:ie,js:je)

    call compute_height_bottom ( Dynam%Vgrid, Var%pssl(is:ie,js:je), &
                                 t_bot, q_bot, z_bot, p_bot,         &
                                 Dynam%Masks%Tmp%kbot(is:ie,js:je)   )


 end subroutine get_bottom_mass

!#######################################################################
! returns u and v on the mass grid at the lowest model level

 subroutine get_bottom_wind (u_bot, v_bot)

  real, intent(out),   &
        dimension(Dynam%Hgrid%Tmp%is:Dynam%Hgrid%Tmp%ie,Dynam%Hgrid%Tmp%js:Dynam%Hgrid%Tmp%je) &
        :: u_bot, v_bot

  real, dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub,Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: u,v,uh,vh

!---- compute lowest level winds on mass grid -----

  call get_bottom_data ( Var % u, Var % v,    &
                         u, v, Dynam%Masks%Vel%kbot )

  call change_grid ( Dynam%Hgrid, WIND_GRID, TEMP_GRID, u, v, uh, vh )

  u_bot = uh (Dynam%Hgrid%Tmp%is:Dynam%Hgrid%Tmp%ie,Dynam%Hgrid%Tmp%js:Dynam%Hgrid%Tmp%je)
  v_bot = vh (Dynam%Hgrid%Tmp%is:Dynam%Hgrid%Tmp%ie,Dynam%Hgrid%Tmp%js:Dynam%Hgrid%Tmp%je)

 end subroutine get_bottom_wind

!#######################################################################

 subroutine get_stock_pe(index, value)

 ! This is a dummy routine.
 ! It is neccessary to satisfy revision 13.0.4.3.2.1 of atmos_coupled/atmos_model.f90
 ! Since that revision of atmos_coupled/atmos_model.f90 does nothing with the result,
 ! this routine can be a dummy.
 ! If and when the result is needed, it should be the total water content of the
 ! global atmosphere (Kg), including vapor, liquid and ice.

 integer, intent(in) :: index
 real, intent(out)   :: value

 value = 0.0
 if(.not.stock_warning_issued) then
   call error_mesg('get_stock_pe','Stocks not yet implemented. Returning zero.',NOTE)
   stock_warning_issued = .true.
 endif

 end subroutine get_stock_pe

!#######################################################################

end module atmosphere_mod

