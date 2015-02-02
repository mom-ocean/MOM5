
module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for B-grid dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------

use bgrid_core_driver_mod, only: bgrid_dynam_type,       &
                                 bgrid_core_driver_init, &
                                 bgrid_core_driver,      &
                                 bgrid_core_time_diff,   &
                                 bgrid_core_driver_end,  &
                                 get_bottom_data,        &
                                 put_bottom_data,        &
                                 atmosphere_domain

use   bgrid_prog_var_mod, only: prog_var_type, var_init

use      bgrid_horiz_mod, only: get_horiz_grid_size,        &
                                get_horiz_grid_bound, TGRID

use     time_manager_mod, only: time_type, get_time, operator(+)

use              fms_mod, only: file_exist, open_namelist_file, &
                                error_mesg, FATAL,              &
                                check_nml_error, stdlog,        &
                                write_version_number,           &
                                mpp_pe, mpp_root_pe,            &
                                close_file, set_domain

! routines used by subroutine bgrid_physics
use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID
use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_vert_mod       , only: vert_grid_type, &
                                 compute_pres_full, compute_pres_half
use bgrid_halo_mod       , only: update_halo, UWND, VWND, TEMP, &
                                 NORTH, EAST, WEST, SOUTH
use hs_forcing_mod       , only: hs_forcing_init, hs_forcing

!-----------------------------------------------------------------------

implicit none
private

public  atmosphere,      &
        atmosphere_init, &
        atmosphere_end,  &
        atmosphere_resolution,  &
        atmosphere_boundary,    &
        get_atmosphere_axes,    &
        atmosphere_domain

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmosphere.F90,v 19.0 2012/01/06 19:52:51 fms Exp $'
character(len=128) :: tag = '$Name: tikal $'

!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,0]

   integer, dimension(2) :: physics_window = (/0,0/)

   namelist /atmosphere_nml/ physics_window

!-----------------------------------------------------------------------
!---- private data ----

type (bgrid_dynam_type), save :: Dynam
type    (prog_var_type), save :: Var, Var_dt
type        (time_type) :: Time_step_atmos

real                               :: dt_atmos
real,    dimension(:,:,:), pointer :: omega =>NULL()
integer, dimension(4)              :: atmos_axes

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine atmosphere (Time)
 type (time_type), intent(in) :: Time

  type(time_type) :: Time_prev, Time_next
!-----------------------------------------------------------------------

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- dynamics -----

   call bgrid_core_driver ( Time_next, Var, Var_dt, Dynam, omega )

!---- call physics -----

   call bgrid_physics ( physics_window, dt_atmos, Time_next,  &
                        Dynam%Hgrid,   Dynam%Vgrid,    Dynam,             &
                        Var,     Var_dt                       )

!---- time differencing and diagnostics -----

   call bgrid_core_time_diff ( omega, Time_next, Dynam, Var, Var_dt )

!-----------------------------------------------------------------------

 end subroutine atmosphere

!#######################################################################

 subroutine atmosphere_init ( Time_init, Time, Time_step )

 type (time_type),     intent(in)    :: Time_init, Time, Time_step

 integer :: unit, sec, ierr, io, is, ie, js, je, i, j
 real, allocatable :: blon(:), blat(:), blon_2d(:,:), blat_2d(:,:)

!-----------------------------------------------------------------------
!----- read namelist -----

    if (file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read (unit, nml=atmosphere_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'atmosphere_nml')
        enddo
 10     call close_file (unit)
    endif

!----- write version and namelist to log file -----

    call write_version_number ( version, tag )
    if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=atmosphere_nml)

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

    is = Dynam%Hgrid%Tmp%is
    ie = Dynam%Hgrid%Tmp%ie
    js = Dynam%Hgrid%Tmp%js
    je = Dynam%Hgrid%Tmp%je
    allocate ( blon(ie-is+2), blon_2d(ie-is+2,je-js+2) )
    allocate ( blat(je-js+2), blat_2d(ie-is+2,je-js+2) )
    call atmosphere_boundary(blon, blat, global=.false.)
    do i=1,ie-is+2
      blon_2d(i,:) = blon(i)
    enddo
    do j=1,je-js+2
      blat_2d(:,j) = blat(j)
    enddo
    call hs_forcing_init ( atmos_axes, Time, blon_2d, blat_2d )
    deallocate(blon, blat, blon_2d, blat_2d)

!   ----- use entire grid as window ? -----

    if (physics_window(1) <= 0) physics_window(1) = ie-is+1
    if (physics_window(2) <= 0) physics_window(2) = je-js+1

!-----------------------------------------------------------------------

 end subroutine atmosphere_init

!#######################################################################

 subroutine atmosphere_end

 integer :: unit

    call bgrid_core_driver_end ( Var, Dynam )

 end subroutine atmosphere_end

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
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_boundary (blon, blat, global)

    real,    intent(out)          :: blon(:), blat(:)
    logical, intent(in), optional :: global

!----- return the longitudinal and latitudinal grid box edges ----------

    call get_horiz_grid_bound (Dynam % Hgrid, TGRID, blon, blat, global)

 end subroutine atmosphere_boundary

!#######################################################################
!    returns the axis indices associated with the coupling grid

 subroutine get_atmosphere_axes ( axes )

   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----

     if ( size(axes) < 0 .or. size(axes) > 4 ) call error_mesg (    &
                    'get_atmosphere_axes in atmosphere_mod', &
                           'size of argument is incorrect', FATAL   )

     axes (1:size(axes)) = atmos_axes (1:size(axes))

 end subroutine get_atmosphere_axes

!#######################################################################

subroutine bgrid_physics ( window, dt_phys, Time, Hgrid, Vgrid, &
                           Dynam, Var, Var_dt )

!-----------------------------------------------------------------------
!
!   Time      =  current time (time_type, see time manager)
!
!-----------------------------------------------------------------------
  integer, intent(in)                :: window(2)
  real,    intent(in)                :: dt_phys
       type(time_type),intent(in)    :: Time
type (horiz_grid_type),intent(inout) :: Hgrid
type  (vert_grid_type),intent(in)    :: Vgrid
type(bgrid_dynam_type),intent(in)    :: Dynam
type   (prog_var_type),intent(in)    :: Var
type   (prog_var_type),intent(inout) :: Var_dt

!-----------------------------------------------------------------------
  integer :: j, k, n, is, ie, js, je, i1, i2, j1, j2, nt
  integer :: ix, jx, idim, jdim
!-----------------------------------------------------------------------

   real, dimension(window(1),window(2),Vgrid%nlev) :: p_full, u_dt, v_dt

   real, dimension(window(1),window(2),Vgrid%nlev+1) :: p_half

   real, dimension(Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: uh, vh, uh_dt, vh_dt

   real, dimension(window(1),window(2)) :: pssl_new
!-----------------------------------------------------------------------
!---------------------------- do physics -------------------------------

    idim = window(1)
    jdim = window(2)

!   --- momentum and momentum tendency on mass grid ---

    call update_halo (Hgrid, UWND, Var_dt%u, halos=SOUTH+WEST)
    call update_halo (Hgrid, VWND, Var_dt%v, halos=SOUTH+WEST)

    if (Dynam%Masks%sigma) then
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID,    Var%u,    Var%v, uh   , vh   )
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID, Var_dt%u, Var_dt%v, uh_dt, vh_dt)
    else
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID,   &
                            Var%u,    Var%v, uh   , vh   , mask_inp=Dynam%Masks%Vel%mask)
        call change_grid (Hgrid, WIND_GRID, TEMP_GRID,   &
                          Var_dt%u, Var_dt%v, uh_dt, vh_dt, mask_inp=Dynam%Masks%Vel%mask)
    endif

!   --- loop through physics windows ---

    nt = Var%ntrace
    js = Hgrid%Tmp%js

    do while ( js <= Hgrid%Tmp%je )

       je = min ( js+jdim-1, Hgrid%Tmp%je )
       jx = je-js+1
       is = Hgrid%Tmp%is

    do while ( is <= Hgrid%Tmp%ie )

       ie = min ( is+idim-1, Hgrid%Tmp%ie )
       ix = ie-is+1

!      ---- pass updated surface pressure ----
       pssl_new(1:ix,1:jx) = Var%pssl(is:ie,js:je) + &
                             Var_dt%pssl(is:ie,js:je) * dt_phys

       call compute_pres_full (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_full(1:ix,1:jx,:))
       call compute_pres_half (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_half(1:ix,1:jx,:))


       u_dt(1:ix,1:jx,:) = uh_dt(is:ie,js:je,:)
       v_dt(1:ix,1:jx,:) = vh_dt(is:ie,js:je,:)


!      ---- j-axis indices in the global physics grid ----

       j1 = js-Hgrid%Tmp%js+1; j2 = j1+(je-js)
       i1 = is-Hgrid%Tmp%is+1; i2 = i1+(ie-is)

!-----------------------------------------------------------------------
!-------------------------- call physics -------------------------------
!------------ (need to add leap-frog option for uh,vh) -----------------
!-----------------------------------------------------------------------
  if (.not.Dynam%Masks%sigma) then
!------------ eta coordinate -------------------------------------------

      call hs_forcing ( i1, i2, j1, j2, dt_phys, Time ,&
                            Hgrid%Tmp%alm(is:ie,js:je)    ,&
                            Hgrid%Tmp%aph(is:ie,js:je)    ,&
                            p_half   ( 1:ix, 1:jx,:)      ,&
                            p_full   ( 1:ix, 1:jx,:)      ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            u_dt     ( 1:ix, 1:jx,:)      ,&
                            v_dt     ( 1:ix, 1:jx,:)      ,&
                            Var_dt%t (is:ie,js:je,:)      ,&
                            Var_dt%r (is:ie,js:je,:,:)    ,&
                            mask=Dynam%Masks%Tmp%mask(is:ie,js:je,:) ,&
                            kbot=Dynam%Masks%Tmp%kbot(is:ie,js:je)    )

  else
!------------- sigma coordinate ----------------------------------------

      call hs_forcing ( i1, i2, j1, j2, dt_phys, Time ,&
                            Hgrid%Tmp%alm(is:ie,js:je)    ,&
                            Hgrid%Tmp%aph(is:ie,js:je)    ,&
                            p_half   ( 1:ix, 1:jx,:)      ,&
                            p_full   ( 1:ix, 1:jx,:)      ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            u_dt     ( 1:ix, 1:jx,:)      ,&
                            v_dt     ( 1:ix, 1:jx,:)      ,&
                            Var_dt%t (is:ie,js:je,:)      ,&
                            Var_dt%r (is:ie,js:je,:,:)     )

  endif

      ! compute momentum tendencies on mass grid for physics only
      !      udt(phys) = udt(current) - udt(before physics)
        uh(is:ie,js:je,:) = u_dt(1:ix,1:jx,:) - uh_dt(is:ie,js:je,:)
        vh(is:ie,js:je,:) = v_dt(1:ix,1:jx,:) - vh_dt(is:ie,js:je,:)

        is = is + idim

     enddo

        js = js + jdim

     enddo

!-----------------------------------------------------------------------
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

end subroutine bgrid_physics

!#######################################################################

end module atmosphere_mod

