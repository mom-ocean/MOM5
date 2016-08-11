             
module bgrid_core_mod

!-----------------------------------------------------------------------
!
!              gfdl global b-grid dynamical core
!
!           (with eta and hybrid pressure coordinate)
!
!-----------------------------------------------------------------------
!--------------------------- modules -----------------------------------
!-----------------------------------------------------------------------

use bgrid_prog_var_mod    , only: prog_var_type, var_init,   &
                                  prog_var_times_scalar
use bgrid_horiz_mod       , only: horiz_grid_type
use bgrid_vert_mod        , only: vert_grid_type, compute_pres_depth, &
                                  compute_pressures
use bgrid_masks_mod       , only: grid_mask_type, grid_masks_init
use bgrid_advection_mod   , only: advection_init, advection, advection_end
use bgrid_horiz_diff_mod  , only: horiz_diff_init, horiz_diff
use bgrid_horiz_adjust_mod, only: horiz_adjust_vel, horiz_adjust_mass,  &
                                  press_grad, compute_grad_pres, div_damping, &
                                  press_grad_fv
use bgrid_vert_adjust_mod , only: vert_adjust
use bgrid_polar_filter_mod, only: pfilt_control_type, polar_filter_init, &
                                  polar_filter, polar_filter_wind, TGRID
use bgrid_halo_mod        , only: update_halo, TEMP, UWND, VWND
use bgrid_sponge_mod      , only: sponge_driver, sponge_init

use                fms_mod, only: error_mesg, FATAL, write_version_number,        &
                                  mpp_clock_id, mpp_clock_begin, mpp_clock_end,   &
                                  MPP_CLOCK_SYNC, CLOCK_MODULE_DRIVER, uppercase, &
                                  stdlog
use          constants_mod, only:  CP_AIR, RDGAS, RVGAS

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index

!-----------------------------------------------------------------------

implicit none
private

public  update_bgrid_core, bgrid_core_init, bgrid_core_end

public  bgrid_dynam_type

!-----------------------------------------------------------------------
!  ----- defined data types -----
!    Hgrid = horizontal grid constants
!    Vgrid = vertical grid constants
!    Masks = eta coordinate topography masks and indices
!    Pfilt = polar filter constants
!
!  ----- 2-dimensional (nlon,nlat) fields -----
!    fis  = geopotential height of the surface
!    fisl = geopotential height at eta=1. (for eta coord = 0.0,
!    res  = reciprical of eta at the surface
!
!  ----- time step terms ----
!    nt_adv = no. of advection time steps per atmosphere step (integer)
!    nt_adj = no. of adjustment time steps per advection step (integer)
!    dt_adj = adjustment time step in seconds (real)
!
!  ----- miscellaneous ----
!    fopt        = filtering option [integer]
!    pgf_method  = pressure gradient algorithm [integer]
!    sphum       = tracer index for specific humidity [integer]
!    coeff_ddamp = coefficient for divergence damping [real]
!    avg_omega   = omega averaging flag [logical]
!    verbose     = verbose flag [integer]
!-----------------------------------------------------------------------

type bgrid_dynam_type
   type(horiz_grid_type), pointer  :: Hgrid => NULL()
   type (vert_grid_type), pointer  :: Vgrid => NULL()
   type (grid_mask_type)           :: Masks
   type (pfilt_control_type)       :: Pfilt
   real, pointer, dimension(:,:)   :: fis  => NULL(), &
                                      fisl => NULL(), &
                                      res  => NULL()
   real                            :: dt_adj
   real                            :: coeff_ddamp
   integer                         :: nt_adv, nt_adj
   integer                         :: fopt, verbose
   integer                         :: pgf_method
   integer                         :: sphum
   logical                         :: avg_omega
end type bgrid_dynam_type

!-----------------------------------------------------------------------
!------- internal options ---------!  alpha_implicit determines how the
                                   !  coriolis and press grad force
   real  :: alpha_implicit = 0.5   !  terms are solved
                                   !    = 0.5  trapezoidal implicit
                                   !    = 1.0        fully implicit


! parameters for pressure gradient scheme
   integer, parameter :: SIMMONS_BURRIDGE=0, FINITE_VOLUME=1

!---- internal parameters ----

   real, parameter :: d608 = (RVGAS-RDGAS)/RDGAS

!---- version number ----

   character(len=128) :: version='$Id: bgrid_core.F90,v 19.0 2012/01/06 19:53:27 fms Exp $'
   character(len=128) :: tagname='$Name: tikal $'

!---- performance timing info ----

   integer :: id_advect, id_adjust
   logical :: do_clock_init = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine bgrid_core_init (Dynam, Hgrid, Vgrid, fis, res, dt, ntadj, ntadv,   &
                              pgf_scheme, filter_option, filter_weight,   &
                              ref_lat_filter, ddamp_coeff, avg_omega,     &
                              verbose)

   type(bgrid_dynam_type), intent(inout)     :: Dynam
   type(horiz_grid_type), intent(in), target :: Hgrid
   type (vert_grid_type), intent(in), target :: Vgrid
   real, intent(in), dimension(:,:),  target :: fis, res
   real,            intent(in)               :: dt
   integer,         intent(in)               :: ntadj, ntadv

!          ---- optional arguments ----

   character(len=*),intent(in), optional :: pgf_scheme
   integer,         intent(in), optional :: filter_option,   &
                                            filter_weight,   &
                                            verbose
   real,            intent(in), optional :: ref_lat_filter,   &
                                            ddamp_coeff
   logical,         intent(in), optional :: avg_omega

!-----------------------------------------------------------------------
!
!    performs initialization for b-grid dynamics type
!
! input:  Hgrid      horizontal grid constants
!         Vgrid      vertical grid constants
!         fis        geopotential height of the surface
!         res        reciprocal of eta at the surface
!
!         dt         adjustment time step in seconds
!         ntadj      number of adjustment time steps for each
!                    advective time step [integer]
!         ntadv      number of advection time steps for each
!                    update call [integer]
!
!         IMPORTANT:  The input arguments (Hgrid, Vgrid, fis, res)
!                     must have space in memory for the entire
!                     duration of the model integration
!
! input (optional):
!
!  filter_option       Determines how polar filtering is performed.
! 
!  filter_weight       Weight applied to the polar filter that will
!                      increase (or decrease) the strength of the standard
!                      polar filter response function.
!
!  ref_lat_filter      The reference latitude at which polar filtering
!                      (in each hemisphere) will begin to be applied.
!
!  ddamp_coeff         Coefficient for divergence damping.
!
!  verbose             Flag that control additional printed output
!                      Currently, this option is not being used.
!
!  avg_omega           return the omega diagnostic averaged over all
!                      adjustment time steps
!
! NOTE: also see bgrid_core_driver for description of optional arguments
!
!-----------------------------------------------------------------------
   integer :: logunit
!  ---- required time step arguments -----

   Dynam % nt_adj = ntadj
   if (Dynam % nt_adj <= 0) call error_mesg ('bgrid_core_init',  &
                            'input argument ntadj must be >= 1', FATAL)

   Dynam % nt_adv = ntadv
   if (Dynam % nt_adv <= 0) call error_mesg ('bgrid_core_init',  &
                            'input argument ntadv must be >= 1', FATAL)

   Dynam % dt_adj = dt / float(Dynam%nt_adv*Dynam%nt_adj)
   if (Dynam % dt_adj <= 0.0) call error_mesg ('bgrid_core_init',  &
                             'input argument dt must be > 0.', FATAL)

!  ---- optional arguments ----

   Dynam % fopt = 2
   if (present(filter_option)) Dynam % fopt = filter_option

   Dynam % verbose = 0
   if (present(verbose)) Dynam % verbose = max(0, verbose)

   Dynam % avg_omega = .false.
   if (present(avg_omega)) Dynam % avg_omega = avg_omega

!  initialize performance clock 
   if (do_clock_init) then
     ! initialize performance timing
       id_advect = mpp_clock_id ('BGRID: advect loop', &
                           flags=MPP_CLOCK_SYNC, grain=CLOCK_MODULE_DRIVER)
       id_adjust = mpp_clock_id ('BGRID: adjust loop', &
                           flags=MPP_CLOCK_SYNC, grain=CLOCK_MODULE_DRIVER)
       do_clock_init = .false.
   endif
!-----------------------------------------------------------------------
!     ------ pointers ------

   Dynam % Hgrid => Hgrid
   Dynam % Vgrid => Vgrid
   Dynam % fis   => fis
   Dynam % res   => res

!-----------------------------------------------------------------------
!     ------- allocate space for other variables -------

   Dynam%fisl => var_init (Dynam%Hgrid)

!-----------------------------------------------------------------------
!  ---- eta coordinate masks ----

   call grid_masks_init ( Dynam%Hgrid, Dynam%Vgrid, Dynam%res, Dynam%Masks )

!  ---- define sea level geop height ----
   if (Dynam%Masks%sigma) then
        Dynam%fisl = Dynam%fis
   else
        Dynam%fisl = 0.0
   endif

!-----------------------------------------------------------------------
!       ---- initialize polar filtering -----

   call polar_filter_init ( Dynam%Pfilt, Dynam%Hgrid, Dynam%Vgrid%nlev,    &
                            reflat=ref_lat_filter, weight=filter_weight,   &
                            sigma=Dynam%Masks%sigma, verbose=Dynam%verbose )

!-----------------------------------------------------------------------
!------- initialization of other bgrid_core modules ----------

     call advection_init  ( Dynam%Hgrid )
     call horiz_diff_init ( Dynam%Hgrid )
     call sponge_init     ( Dynam%Hgrid )

!-----------------------------------------------------------------------
!------- save other settings -----

   ! divergence damping
     Dynam % coeff_ddamp = ddamp_coeff     
   ! tracer index for specific humidity
     Dynam % sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
   ! pressure gradient scheme option
     if (uppercase(trim(pgf_scheme)) == 'SIMMONS_BURRIDGE') then
            Dynam % pgf_method = SIMMONS_BURRIDGE
     else if (uppercase(trim(pgf_scheme)) == 'FINITE_VOLUME') then
            Dynam % pgf_method = FINITE_VOLUME
            ! cannot use FV pgf with shallow water version (nlev=1)
            if (Dynam%Vgrid%nlev == 1) call error_mesg ('bgrid_core_init', &
                                    'invalid PGF scheme when nlev=1', FATAL)
     else if (uppercase(trim(pgf_scheme)) == 'DEFAULT') then
            Dynam % pgf_method = SIMMONS_BURRIDGE
     else
            call error_mesg ('bgrid_core_init',  &
                             'invalid PGF scheme', FATAL)
     endif

!-----------------------------------------------------------------------

     call write_version_number(version, tagname)

     logunit = stdlog()
     write (logunit,10) dt, dt/real(Dynam%nt_adv), Dynam%dt_adj
  10 format (/,'dynamical core time step (seconds) = ',f7.2, &
             /,'     advective time step (seconds) = ',f7.2, &
             /,'    adjustment time step (seconds) = ',f7.2,/)

!-----------------------------------------------------------------------

 end subroutine bgrid_core_init

!#######################################################################

 subroutine update_bgrid_core (Var, Var_dt, Dynam, omega, div, mfew, mfns )

!-----------------------------------------------------------------------
! Var       = prognostic variables
! Var_dt    = prognostic variable tendencies
! Dynam     = data type for dynamical core constants and data
! omega     = omega (vertical velocity) diagnostic (Pa/s)
!-----------------------------------------------------------------------
   type (prog_var_type), intent(in)    :: Var
   type (prog_var_type), intent(inout) :: Var_dt
   type(bgrid_dynam_type), intent(inout) :: Dynam
   real, intent(out), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                                Dynam%Hgrid%jlb:Dynam%Hgrid%jub, &
                                Dynam%Vgrid%nlev) :: omega, div, mfew, mfns
!-----------------------------------------------------------------------

real, dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: psdt

real, dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: pssl

real, dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub, Dynam%Vgrid%nlev) :: &
                dpde, few, fns, divp, pgfew, pgfns, dpde_old,          &
                pfull, wta, wtb, cew, cns, u, v, tq, uo, vo, to, omgalf

real, dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub, Dynam%Vgrid%nlev+1)  &
                :: phalf, etadot

   real    :: fadv, tdt_adj, scale
   integer :: i, k, n, m, nt
!-----------------------------------------------------------------------
!    ---- definition of local variables ----
!
!    psdt  = surface pressure tendency (Pa/s)
!    pssl  = surface pressure adjusted to eta=1.
!    tq    = (virtual) temperature
!    divp  = mass divergence (Pa/s)
!    pgfew = zonal pressure gradient force (at v pts) (m/s2)
!    pgfns = meridional pressure gradient force (at v pts) (m/s2)
!    dpde  = pressure thickness of model layers (Pa)
!    pfull = pressure at full model levels
!    phalf = pressure at half model levels (between full levels)
!    wta,
!      wtb = weights for obtaining the pressure at full levels (no units)
!    cew,
!      cns = grad(p)/p term (zonal and meridional components) (no units)
!    few,
!      fns = mass fluxes time summed over advection interval (Pa-m2/s)
!            (zonal and meridional components)
!   etadot = vertical mass flux summed over advection interval (Pa/s)
!   omgalf = thermodynamic (omega-alpha) term (Pa/s)
!   uo,vo,to = momentum and temperature at the start of an advective time step
!   dpde_old = pressure thickness of model layers at the start of an advective time step
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------- number of tracers and variable (time) levels ---------

    nt = Var%ntrace

!   ---- set-up time step ----

    tdt_adj = Dynam % dt_adj

!   ------ pressure variables ------

    pssl = Var % pssl + tdt_adj * Var_dt % pssl

    call compute_pressures ( Dynam%Vgrid, pssl,           &
                             phalf, pfull, dpde, wta, wtb )

    call compute_grad_pres ( Dynam%Hgrid, Dynam%Vgrid%nplev, &
                             phalf, dpde, wta, wtb, cew, cns )

!   ------ zero fluxes and diagnostic output -----

    few     = 0.0;  fns     = 0.0;  etadot     = 0.0

    omega = 0.0;  div = 0.0;  mfew = 0.0;  mfns = 0.0

!-----------------------------------------------------------------------
!************** start of bgrid_core time step loop *********************
call mpp_clock_begin (id_advect)

 do m = 1, Dynam % nt_adv

!   ------ save this time level for advection ------

    dpde_old = dpde

    to = Var % t + tdt_adj * Var_dt % t
    uo = Var % u + tdt_adj * Var_dt % u
    vo = Var % v + tdt_adj * Var_dt % v

!   ------ variables at current time level ------

    tq = to
    call compute_virt_temp ( Dynam%sphum, tq, Var%r, Var_dt%r, tdt_adj )
    u  = uo
    v  = vo

!-----------------------------------------------------------------------
!*************** start of adjustment time step loop ********************
call mpp_clock_end   (id_advect)
call mpp_clock_begin (id_adjust)

 do n = 1, Dynam % nt_adj

!-----------------------------------------------------------------------
!------------------------- compute mass divergence ---------------------
!--------------------- horizontal alpha-omega term ---------------------

   call horiz_adjust_mass ( Dynam%Vgrid%nplev, Dynam%Hgrid, &
                            Dynam%Masks,       u,           &
                            v,                 dpde,        &
                            cew,               cns,         &
                            few,               fns,         &
                            divp,              omgalf       )


!    ---- polar filtering -----

     if (Dynam%fopt >= 1) then
       call polar_filter (Dynam%Pfilt, divp, omgalf, TGRID, Dynam%Masks%Tmp%mask)
     endif

     call update_halo (Dynam%Hgrid, TEMP, divp  )
     call update_halo (Dynam%Hgrid, TEMP, omgalf)

!--- save divergence diagnostic (units: 1/sec) ---
     div = div + divp/dpde

!-----------------------------------------------------------------------
!------- compute sfc pres tendency, vert. alpha-omega, & vert vel.------

    call vert_adjust ( Dynam%Vgrid, Dynam%res, divp, wta, wtb,    &
                       Dynam%Masks%Tmp%mask, omgalf, etadot, psdt )     

!-----------------------------------------------------------------------
!-------- update surface pressure tendency -----------------------------

    Var_dt % ps   = Var_dt % ps   + psdt
    Var_dt % pssl = Var_dt % pssl + psdt * Dynam % res

!-------- recompute pressure variables at next time level -----------

    pssl = Var % pssl + tdt_adj * Var_dt % pssl
    call compute_pressures ( Dynam%Vgrid, pssl, phalf, pfull, dpde, &
                             wta, wtb )

!----- do not execute the following code -----
!   with the shallow water version (nlev=1)
    if (Dynam%Vgrid%nlev > 1) then
!-------- update thermodynamic tendency ----------

       omgalf = omgalf/dpde
       Var_dt % t = Var_dt % t + omgalf * tq * RDGAS/CP_AIR

!  ---- compute omega diagnostic (Pa/s) ----

       if (Dynam%avg_omega) then
          omega = omgalf * pfull + omega
       else
          omega = omgalf * pfull
       endif

    endif

!-----------------------------------------------------------------------
!------------ compute geopotential height and rt/p ---------------------
!     ----- (use smoothed value for pssl if leapfrog) ----

    tq = Var % t + tdt_adj * Var_dt % t
    call compute_virt_temp ( Dynam%sphum, tq, Var%r, Var_dt%r, tdt_adj )

    call compute_grad_pres ( Dynam%Hgrid, Dynam%Vgrid%nplev, &
                             phalf, dpde, wta, wtb, cew, cns )


  select case (Dynam%pgf_method)
    case (SIMMONS_BURRIDGE)
      call press_grad ( Dynam%Hgrid, Dynam%Vgrid, Dynam%Masks,    &
                        Dynam%fisl, tq, dpde, wta, wtb, cew, cns, &
                        pgfew, pgfns )
    case (FINITE_VOLUME)
      call press_grad_fv ( Dynam%Hgrid, Dynam%Vgrid, Dynam%Masks, &
                           Dynam%fisl, tq, phalf, wta, wtb,       &
                           pgfew, pgfns )
  end select

!------------------- adjustment of wind components ---------------------

    call horiz_adjust_vel ( Dynam%Hgrid, Dynam%Masks,                &
                                        tdt_adj, pgfew, pgfns,       &
                            u, v, Var_dt%u, Var_dt%v, alpha_implicit )

  ! polar filtering of momentum (new scheme)
    if (Dynam%fopt == 2)  &
    call prog_var_filter_vel ( Dynam%Pfilt, tdt_adj, Dynam%Masks%Vel%mask, &
                               Var%u, Var%v, Var_dt%u, Var_dt%v            )

!-----------------------------------------------------------------------
!----------------------- advection -------------------------------------

    if (n ==  Dynam%nt_adj) then
        call mpp_clock_end   (id_adjust)
        call mpp_clock_begin (id_advect)

        ! accumulate mass fluxes for diagnostic output
              mfew = mfew + few
              mfns = mfns + fns

         call advection ( Dynam%Pfilt,                           &
                          Dynam%Hgrid, Dynam%Vgrid, Dynam%Masks, &
                          Dynam%fopt, tdt_adj, dpde_old, dpde,   &
                          few, fns, etadot,                      &
                          uo, vo, to,   Var, Var_dt              )

        call mpp_clock_end   (id_advect)
        call mpp_clock_begin (id_adjust)
    endif

!-----------------------------------------------------------------------

  ! polar filtering of momentum (old scheme)
    if (Dynam%fopt == 1)  &
    call prog_var_filter_vel ( Dynam%Pfilt, tdt_adj, Dynam%Masks%Vel%mask, &
                               Var%u, Var%v, Var_dt%u, Var_dt%v            )

  ! damping using filtered divergence
    if (abs(Dynam % coeff_ddamp) > 1.e-8) then
      call div_damping ( Dynam%Hgrid, Dynam%Vgrid, Dynam%Masks, &
                         tdt_adj, Dynam % coeff_ddamp,          &
                         dpde, divp, Var_dt%u, Var_dt%v )
    endif

!   ---- recompute momentum at next time level ----

    call update_halo (Dynam%Hgrid, UWND, Var_dt%u )
    call update_halo (Dynam%Hgrid, VWND, Var_dt%v )

    u = Var % u + tdt_adj * Var_dt % u
    v = Var % v + tdt_adj * Var_dt % v

!-----------------------------------------------------------------------

 enddo
 call mpp_clock_end   (id_adjust)
 call mpp_clock_begin (id_advect)

!-----------------------------------------------------------------------
!----------------- horizontal diffusion --------------------------------

    call horiz_diff ( Dynam%Hgrid, Dynam%Masks, Dynam%Vgrid%nplev, &
                      tdt_adj, dpde, pfull, Var, Var_dt            )

!-----------------------------------------------------------------------
!----------------- sponge at top of model -------------------

    call sponge_driver ( Dynam%Hgrid, Dynam%Vgrid%nplev, &
                         tdt_adj, dpde, Var, Var_dt      )

!    ---------------- halo updates --------------------
!    skip on the last pass will update last when needed

    if (m == Dynam % nt_adv) cycle
    call update_halo (Dynam%Hgrid, TEMP, Var_dt%t )
    call update_halo (Dynam%Hgrid, TEMP, Var_dt%r )
    call update_halo (Dynam%Hgrid, UWND, Var_dt%u )
    call update_halo (Dynam%Hgrid, VWND, Var_dt%v )

!-----------------------------------------------------------------------

 enddo
 call mpp_clock_end   (id_advect)

!-----------------------------------------------------------------------
!  ---- scale tendencies before physics for time step difference ----

   scale = 1.0 / float(Dynam%nt_adj*Dynam%nt_adv)
   call prog_var_times_scalar (Var_dt, scale)

  !-- return time averaged diagnostic quantities --
   if (Dynam%avg_omega) omega = omega * scale
   div  = div  * scale
   mfew = mfew * scale
   mfns = mfns * scale

!-----------------------------------------------------------------------

 end subroutine update_bgrid_core
        
!#######################################################################

 subroutine bgrid_core_end (Dynam)

!-----------------------------------------------------------------------
   type(bgrid_dynam_type), intent(inout) :: Dynam
!-----------------------------------------------------------------------
! Dynam = data type for dynamical core constants and data

   call advection_end

!-----------------------------------------------------------------------

 end subroutine bgrid_core_end

!#######################################################################

 subroutine compute_virt_temp ( sphum, tq, r, rdt, dt )

  integer, intent(in) :: sphum
  real, intent(inout) :: tq(:,:,:)
  real, intent(in)    :: r (:,:,:,:), rdt(:,:,:,:), dt

  real, dimension(size(r,1),size(r,2),size(r,3)) :: q

! compute virtual temperature
! using temperature and specific humidity (in the tracer array)

    if ( size(r,4) == 0 ) return
    if ( sphum <= 0 ) return

    q  = r(:,:,:,sphum) + dt * rdt(:,:,:,sphum)
    tq = tq * ( 1. + d608 * q )

 end subroutine compute_virt_temp

!#######################################################################
!-------- polar filter routine for momentum components ---------

 subroutine prog_var_filter_vel ( Pfilt, dt, mask, u, v, udt, vdt )

 type(pfilt_control_type), intent(in)    :: Pfilt
 real,                     intent(in)    :: dt
 real, dimension(:,:,:),   intent(in)    :: mask
 real, dimension(:,:,:),   intent(in)    :: u, v
 real, dimension(:,:,:),   intent(inout) :: udt, vdt

   real, dimension(size(u,1),size(u,2),size(u,3)) :: ut, vt

!  ---- momentum ----
     ut = u + dt*udt
     vt = v + dt*vdt
     call polar_filter_wind ( Pfilt, ut, vt, mask )
     udt = (ut-u)/dt
     vdt = (vt-v)/dt

 end subroutine prog_var_filter_vel

!#######################################################################
!#######################################################################

end module bgrid_core_mod

